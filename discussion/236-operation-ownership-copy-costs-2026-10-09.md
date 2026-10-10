# 236 — Operation ownership and copy costs

Status: Proposed
Date: 2026-10-09
Relates: [181](181-python-boundary-ownership-2026-08-03.md),
[168](168-api-hygiene-2026-07-27.md),
[213](213-editor-overlay-storage-2026-08-27.md)

## Scope

This document owns the ownership contracts of graph-IR and graph operations:
which operations consume, borrow, or mutate their input, how they are named,
and the copy costs those contracts imply. It also owns copy-on-write storage for
the three graph-IR aggregates.

[181](181-python-boundary-ownership-2026-08-03.md) remains the plan of record
for the molecule lifecycle, recovery, and Python ownership. Its plan follows the
decisions recorded here.

The operations are reviewed in eight steps: leaf forms and the normal-form
tower; two-operand operations; constraint containers; entries, frame actions,
and batches; Molecule as a value; Molecule as a published object; graph
engines; and Python. All eight are decided.

## Decided: the normal-form tower

1. normalize, reframe, canonicalize, and their tracked forms consume their
   input. They are logical analogs and share one contract. FrameTransport's
   reframe_by also consumes.
2. A consuming operation moves its input into its result. It does not
   internally clone an input it already owns.
3. The equalities normalized_eq, framed_eq, and canonical_eq borrow their
   operands. A default may clone and run the consuming operation.
   Implementations override it with a borrowing computation where cost matters.
4. Normalize::normalized is removed. A caller holding a reference writes
   `x.clone().normalize()`. The scalar forms override normalized_eq with a
   borrowing comparison of their trivially normal variants. This supersedes
   168's to_normalized rename.
5. canonical_hash borrows. Hashing does not transform its input, so a consuming
   signature does not fit its meaning. It consumes today only because it was
   added as a default on top of the consuming canonicalize.

### Removal of normalized: measurement

A worktree experiment at f3d902498 removed normalized and its 11 overrides and
migrated every call site to `clone().normalize()`. That includes the trait
defaults, the scalar meet/join/matches implementations, the entity-entry
normalize, and the 39 canonicalization key sites. It did not add the scalar
normalized_eq overrides.

With all features, 8,787 graph-IR library tests, 107 lattice property tests,
and 1,993 graph library tests passed. One property assertion was removed: it
checked only that normalized agreed with normalize.

| Benchmark | Cases | Median change | Mean change | Outside noise threshold |
| --- | ---: | ---: | ---: | --- |
| graph-IR canonicalize | 111 | +0.10% | −0.12% | 17 faster, 13 slower, within ±4% except one −9.6% case |

Canonicalize cases that do not call normalized moved by −3.8% to +2.9%, which
is the noise level of these runs. The substructure bench was not compared: for
literal patterns its atom and bond comparison reaches normalized only through
expression patterns. The composite resolver calls meet on scalar forms, the
changed path. The load average was about 3 during the
baseline and 5–6 during the comparison.

### Full-level canonicalization clone

Canonicalize::canonicalize for Molecule passes `&self` to the level dispatch. At
the Full level, canonicalize_full_with_options receives `&Molecule` and clones
it to obtain an owned value for normalization. The key computation, partition,
search, and every leaf candidate read that normalized molecule. Because `self`
is still alive, normalizing the clone copies the atom and bond tables and every
overlay table that normalization writes. The consuming route already owns the
molecule, so that copy is avoidable. The other levels have no such clone.

Sources: [traits](../umol-graph-ir/src/ir/traits.rs),
[canonicalization](../umol-graph-ir/src/ir/canonicalize.rs).

## Decided: two-operand operations

Two-operand operations borrow both operands. The exception is an infallible
operation whose result extends the receiver: it consumes the receiver and reuses
its storage, as `String + &str` and Vec::extend do.

| Operation | Contract |
| --- | --- |
| Molecule::combine | `combine(self, other: &Molecule) -> Molecule`: infallible, the result is the receiver with other appended, so it consumes. Python clones the receiver and consumes the clone. |
| Molecule::combine_from | In-place counterpart, unchanged signature; combine is combine_from followed by returning the receiver. |
| Molecule::combine_all | n-ary without a receiver; borrows every input and builds once. |
| meet, join, frame-action compose | The result does not extend the receiver; borrow both. |
| meet_pushout, tracked_meet_pushout | Keep the receiver's ids but can fail; borrow both, so failure does not lose the receiver. |
| Reaction::compose | Returns several reactions; borrows both. |

combine_from adds other's entities one at a time, and each add_bond rebuilds
the graph. The batch methods add_atoms, add_bonds, and the overlay batch
methods already append with one graph extension.

## Decided: constraint containers

Collection methods on the per-entity `*ConstraintsForm` types and on Constraints
take `&mut self`, following std collections.

update has two meanings, each uniform within its kind. Every form's
`update(&self, &…Update) -> Self` applies an update record and returns a new
form; its production callers in Edits hold a borrowed stored form. Every
`*ConstraintsForm::update(&mut self, &other)` merges a same-type container in
place by key, as Python's dict.update does.

extend takes `impl IntoIterator` on every container.

compact consumes on every public type: the per-entity constraint forms,
Constraint, RelationalConstraint, MoleculeConstraint, the per-entity containers,
and Constraints. tracked_compact follows compact's ownership. restore stays
`&mut self`.

## Decided: entries, frame actions, and batches

No changes. Entity entries are crate-private views whose tower operations
consume. The frame-action algebra borrows: identity, inverse, and compose.
Edits is a builder whose methods take `&mut self` and which application
consumes. Deltas is a collection whose methods take `&mut self` and whose
normalize consumes. inverse consumes on all eight field-change types, the nine
per-entity delta types, and Delta, and borrows on all five frame-action types.

## Decided: transport operations

map, try_map, remap, and try_remap restate references under a correspondence or
remapping. They are logical analogs of reframe_by and compact and consume on
every graph-IR type: the entity collections, Molecule, ReactionSpan,
StereoLigand, Constraint, MoleculeConstraint, and RelationalConstraint.
MoleculeMetadata::remap takes a correspondence and is renamed map, in Rust and
Python.

The five graph-core relation-set types keep their borrowing map, try_map,
remap, and try_remap. They rebuild the set and clone each payload. A consuming
version would save only the payload clones, one per overlay entity, which is
not a strong enough argument for changing graph-core's API.

## Decided: constraint placement

lift_constraints and inline_constraints restate where a molecule's constraints
are stored. They consume: `lift_constraints(self) -> Self` and
`inline_constraints(self) -> Result<Self, Contradiction>`. A contradiction is
intrinsic to the value. inline moves the constraints instead of cloning them
into a validation plan. 166 still owns lift's treatment of an undetermined
stereo kind.

## Decided: derivations

extract, tracked_extract, edits, split, and tracked_split keep borrowing. They
derive results from a source that stays the primary object; HMO perception, the
Kekulizer, and Python keep the host after extracting.

## Decided: Molecule editing

Molecule::edit, apply, and tracked_apply consume, as they do today. The editor
owns its working molecule and writes to its storage directly; a caller that
keeps the molecule clones it first, as reaction application does. In-place
editing exists only with recovery, through edit_into and apply_into, because a
borrowed editor without rollback would leave a published Molecule in an
unchecked state after failure. into_editor is not added; it would duplicate
`edit(self)`.

Python edit, apply, and tracked_apply bind as clone followed by the Rust
method, so the source and its views stay valid; Edits are consumed. Python's
in-place route with recovery is transact and tracked_transact.

## Decided: graph engine naming

The bare verb consumes its molecule. `_into` mutates a borrowed molecule with
recovery. A caller that keeps its molecule clones it explicitly. This covers the
composite resolve, resolve_with_report, and project, the individual resolvers'
resolve and project, Transformer::transform, and their `_into` counterparts.
Python binds the bare verb as clone followed by the Rust method, so the Python
source survives every outcome.

## Decided: planning and the remaining engines

Planned edits apply to the molecule they were planned from. A plan that fails to
apply is a planner defect, not an operation outcome. The five individual
resolvers' application of their own plans, and the composite resolver's
application of phase plans, therefore become plan invariants. The Apply variants
of AromaticityError, AromaticityProjectError, BondsError, IsotopeError,
IsotopeProjectError, MulticenterBondsError, StereoError, and StereoProjectError
are produced only there and are removed.

AromaticityPerceiver::add_systems consumes: `add_systems(&self, molecule:
Molecule, systems) -> Molecule`. Its edits are applied through an editor; it
cannot fail, so it needs no transaction. It has no production callers.

Kekulizer establishes validity during planning. A valid aromatic system admits
a matching of the size its electron count and charge require; planning already
checks this. Planning also checks that the charge, lone pairs, and unpaired
electrons it proposes for each matched and exposed atom are compatible with the
atom's current values. The whole-molecule validation after publication, and its
four PostLocalization error variants, are removed.

Because every transformer failure occurs during planning, transform_into takes
no RollbackMethod and applies the plan in place without a transaction. The same
holds for the individual resolvers' resolve_into and project_into: their
underdetermination, contradiction, and errors are decided during planning.

Consuming operations do not hand their input back on failure. A caller that may
need the input afterward clones it first. Handing the input back would widen
every result type and fix a contract the future may not want. The
composite resolver's `_into` methods keep RollbackMethod, because its later
phases can fail after earlier ones have changed the molecule.

The hydrogen transformers get no separate graph-IR functions. They follow the
transformer pattern: a borrowed, fallible plan and its application. 166's public
fold_explicit_hydrogens and unfold_implicit_hydrogens are withdrawn.

Composite resolution and Molecule apply genuinely fail after mutation: each
later resolution phase plans against the result of the earlier ones, and
Molecule apply executes caller-supplied edits.

## Copy-on-write storage

Molecule's graph, atom and bond tables, and six overlay collections are
Arc-backed. A clone increments reference counts, and the first write to a shared
table copies that table. Three storages do not follow this rule. Cloning them
copies their contents.
[213](213-editor-overlay-storage-2026-08-27.md#remaining-costs) recorded the
Molecule constraint copy as a remaining cost.

| Aggregate | Field | Current type | Shared on clone |
| --- | --- | --- | --- |
| Molecule | constraints | Constraints, wrapping `Vec<Constraint>` | No |
| Reaction | lhs | Molecule | As Molecule |
| Reaction | deltas | Deltas, wrapping `Vec<Delta>` | No |
| ReactionSpan | graph | Graph | Yes |
| ReactionSpan | atoms, bonds | `Vec<EntitySpan<AtomForm>>`, `Vec<EntitySpan<BondForm>>` | No |
| ReactionSpan | six overlay span sets | Newtypes wrapping their relation set directly | No |
| ReactionSpan | constraints | `Vec<ConstraintSpan>` | No |

Known clone sites that copy this storage:

- Every Molecule clone copies its constraint list. Most molecules have none,
  and cloning an empty Vec allocates nothing.
- Reaction::apply and Reaction::tracked_apply snapshot the reaction into the
  returned iterator, copying its Deltas.
- Reaction and ReactionSpan use the default framed_eq, which clones both
  operands.
- Python ReactionSpan canonicalize and tracked_canonicalize clone the wrapped
  span.

The application iterator also clones its prepared Deltas for every match. Each
match reframes its own copy, so that copy remains.

### Design

The layout follows Molecule's existing pattern for each storage role.

- Plain tables and public containers get an Arc at the aggregate field, as
  Molecule's atom and bond tables do: Molecule `constraints: Arc<Constraints>`;
  Reaction `deltas: Arc<Deltas>`; ReactionSpan
  `atoms: Arc<Vec<EntitySpan<AtomForm>>>`,
  `bonds: Arc<Vec<EntitySpan<BondForm>>>`, and
  `constraints: Arc<Vec<ConstraintSpan>>`.
- ReactionSpan's six overlay span newtypes wrap their relation set in an Arc, as
  Molecule's DativeBonds and its peers do.

Reads go through Deref. Mutating sites use Arc::make_mut. Consuming paths that
rebuild a table use Arc::unwrap_or_clone, which moves a uniquely owned table and
copies a shared one. Prefer Arc::make_mut in place to mem::take, because taking
a default Arc allocates an empty one.

Known costs:

- Every construction allocates one Arc per wrapped table, including empty
  tables: one more per Molecule, one per Reaction, and nine per ReactionSpan.
  Molecule already pays this for its overlays.
- The first write to a shared table makes two allocations and copies the table,
  as 213 measured for atom and bond tables.
- Reads add one pointer indirection.

Rejected alternative: an Arc inside Constraints and Deltas. It changes public
container types beyond the three aggregates.

The change is internal. No public signature, constructor, error, or Python API
changes.

### Measurement before production

A worktree experiment at f3d902498 implemented the design for all three
aggregates. With all features, 8,787 graph-IR library tests, 417 graph-IR
property tests, 1,993 graph library tests, and 4,275 IO library tests passed.
Twelve micro-benchmarks measured construction and clone before and after:

| Case | Before | After | Change |
| --- | ---: | ---: | ---: |
| Molecule construction, benzene | 728 ns | 742 ns | +1.4% |
| Molecule construction, 100-atom chain | 2.19 µs | 2.22 µs | +1.6% |
| ReactionSpan construction, 6 atoms | 3.04 µs | 3.07 µs | +2.7% |
| ReactionSpan construction, 100 atoms | 17.7 µs | 17.8 µs | +0.8% |
| Reaction::new, 100 atoms | 29 ns | 35 ns | +19% |
| Molecule clone, benzene | 38 ns | 34 ns | −11% |
| Molecule clone, 100-atom chain | 38 ns | 34 ns | −11% |
| Molecule clone, 100 molecule-level constraints | 777 ns | 34 ns | −96% |
| ReactionSpan clone, 6 atoms | 473 ns | 34 ns | −93% |
| ReactionSpan clone, 100 atoms | 3.84 µs | 34 ns | −99% |
| Reaction clone, 6 and 100 atoms | 65 ns | 37 ns | −42% |

The construction cost is one Arc allocation per wrapped table, about 6 ns for
Reaction::new. The load average was about 10 at the start of the comparison and
about 3 during the baseline, so the small construction changes are within
noise.

## Outside this document

Calling engines with the molecule as receiver, so that engine calls chain left to
right, needs a separate design.

## Implementation plan

This plan runs before 181's. Both plans together are the 0.10 release blocker.

No production change starts before its subitem is written here. Steps 2–8 gain
subitems only after their decisions.

### S0 — Normal-form tower

#### S0a — Remove Normalize::normalized **Done**

**Modules:** graph-IR traits, the 11 scalar forms, lattice implementations,
entity entries, canonicalization, property tests. **Breaking, red→green within
the subitem.** [dep: none]

Remove normalized from Normalize and its 11 overrides. Migrate every caller to
`clone().normalize()`. Override normalized_eq in the scalar forms to compare
trivially normal variants without cloning. The stereo forms delegate to
coset_normalized_eq, a kind-parameterized coset comparison beside coset_meet.
Remove the property assertion that compared normalized with normalize. Update
rustdoc that refers to normalized.

**Tests:** normalized_eq overrides agree with the default for every variant
pair of each scalar form, including unsatisfiable values.

**Verification:** the canonicalize and resolve benches against a baseline taken
before the change; the composite resolver calls meet on scalar forms.

#### S0a1 — Normalized reaction-span constraint keys **Done**

**Modules:** graph-IR canonicalization and its tests. **Red→green
within the subitem.** [dep: none]

ReactionSpan canonical_eq is not invariant under renumbering when the span
carries a constraint with an id list. Molecule-level subsets (ChargeSum,
UnpairedElectronCoupling, BondOrderSum, Connected) and relational atom lists
(DativeBondDonors, DativeBondContainsAllDonors) are sets whose normal form is
sorted and deduplicated. Each span search leaf keys
`span.remap(&remapping).reframe_by(&action)`, which leaves the remapped lists in
permuted order, and the molecule-level constraint block of
reaction_span_comparison_key reads them in stored order. Molecule is
unaffected, and the span's own final canonicalization step normalizes.

reaction_span_comparison_key keys each constraint span after normalizing it, as
it already keys normalized entity values. Spans without constraints take the
unchanged path.

**Tests:** canonical_eq of a span with a molecule-level subset constraint and of
a span with a relational atom-list constraint against a renumbered copy; the
same pairs through Reaction canonical_eq.

**Verification:** a short run of the reaction canonical_eq cases against a
baseline taken before the change.

Molecule-level constraints are a set under normalized identity: normalized-equal
stored constraints, such as `c` and `And([c])`, count as one constraint, and every
constraint delta acts on all of them. Materialization, ReactionSpan::superimpose
(and with it from_sides and composition), Molecule normalization, and canonical
comparison keys follow this definition, as the data-type guide states.
Application and reaction-span integrity do not; S0a2 and S0a3 bring them in line.

#### S0a2 — Conflicting constraint spans in reaction-span integrity **Done**

**Modules:** graph-IR reaction_span and its tests. **Red→green within the
subitem.** [dep: none]

ReactionSpan integrity accepts entries that give normalized-equal constraints two
span states, such as `Unchanged(c)` with `Added(c)`. Each pairing of different states
asserts the constraint both present and absent on one side, and to_reaction turns
such a span into a reaction that does not materialize. check_integrity rejects
entries in which two constraint spans of different states carry
normalized-equal constraints, with the new variant
`ReactionSpanIntegrityError::ConflictingConstraintSpans { constraint:
Constraint }`. Normalized-equal constraints under one state remain valid.

**Tests:** try_from_entries rejects `Unchanged` with `Added`, `Unchanged` with
`Removed`, and `Added` with `Removed` of one constraint, also when the second
constraint is a non-normal equivalent such as `And([c])`, and accepts
normalized-equal constraints under one state. These cases replace
test_reaction_span_to_reaction_constraint_conflict, whose spans become invalid.

#### S0a3 — Constraint removal acts on all normalized-equal constraints **Done**

**Modules:** graph-IR molecule apply and reaction tests. **Red→green within the
subitem.** [dep: none]

RemoveMoleculeConstraint removes the one stored constraint that is structurally
equal to the edit's. A reaction whose lhs stores `And([c])` and whose delta
removes `c` therefore materializes but fails to apply with MissingEntry, and a
host holding two equivalent copies keeps one. RemoveMoleculeConstraint removes
every stored constraint normalized-equal to the edit's and fails with
MissingEntry when there is none; the undo payload records each removed constraint
with its position, as it already supports.

**Tests:** test_reaction_apply_at_non_normal_constraint asserts the exact
product; a molecule-apply case removes two equivalent stored copies through one
normal-form edit and its undo restores both. The `equivalent` case of
test_molecule_apply_edit_molecule_constraint_error moves to the positive cases.

#### S0a4 — Set-like constraints in reaction-span construction **Done**

**Modules:** graph-IR reaction_span, the reaction property tests. **Red→green
within the subitem.** [dep: S0a2]

to_reaction_span and ReactionSpan::superimpose build the constraint spans from
the normalized constraint sets of the two sides, in constraint order, so a
span carries no normalized-equal duplicates. Their rustdoc states the set-like
semantics: the lhs projection has the same structural entities and canonically
equal constraint-set semantics. A redundant or non-normal stored lhs constraint
therefore does not return exactly through a span. test_reaction_span_roundtrip
and the two `span.lhs() == host` assertions of the application properties
compared the lhs exactly; they compare it with its molecule-level constraints in
normal form. to_reaction emits one delta per normal form.

**Tests:** test_reaction_to_reaction_span_constraint covers a non-normal stored
constraint removed by its stored form and `c` beside `And([c])` under one
Remove, each with one Removed span in normal form.

#### S0a5 — Top-level conjunctions stay one constraint **Done**

**Modules:** graph-IR constraint and reaction_span and their tests; the
data-type guide. **Red→green within the subitem.** [dep: S0a4]

Constraints normalization flattened a top-level `And` into its conjuncts, while
it kept a top-level `Or` as one constraint, and every constraint operation
identifies a stored constraint by its own normal form. A canonical reaction then
removed the conjuncts separately and no longer applied to a host that stores the
conjunction. Constraints normalization normalizes each stored constraint on its
own, so a top-level `And` stays one constraint, as an `Or` does; empty `And` and
`Or` are dropped, and the rest sorts and deduplicates. Reaction-span
normalization derives each constraint span's state from the two normalized
sides, ordered by constraint, so the span's lhs and rhs equal the normalized
sides. Materialization emits the constraint spans in the same constraint order, so a
normalized span rematerializes in its own order. The data-type guide
states the rule.

**Tests:** test_constraints_normalize keeps a two-conjunct top-level `And`;
test_reaction_span_normalize_constraints covers added and removed conjunctions
beside an unchanged conjunct and an empty conjunction;
test_reaction_span_reframe_projection and test_reaction_canonicalize_application
hold on the S0a6 generators.

#### S0a6 — Molecule-level constraints in aggregate generators **Done**

**Modules:** graph-IR property-test strategies and the reaction property tests.
**Additive, green.** [dep: S0a1–S0a5]

The Molecule properties draw molecule-level constraints from
molecule_with_constraints_strategy; the ReactionSpan and Reaction generators
produced none. materializable_reaction_strategy gives a reaction one or two lhs
constraint_strategy trees with probability 1/20, each removed by a constraint
delta with probability 1/2, and independently one added tree with probability
1/20; reaction_span_strategy gives an entries-built span one or two constraint
spans with probability 1/10, in constraint order. Trees range over the value's
own entity counts, and a draw that fails molecule integrity, span integrity, or
materialization keeps the unconstrained value. Matching rejects patterns with
molecule-scope constraints until doc 195's evaluation lands, so the application
properties that match (test_reaction_apply, test_react_react,
test_react_react_slice) draw from matchable_reaction_strategy, which adds
constraints but gives the lhs none; the apply_at properties keep them.
molecule_with_removals_strategy no longer removes a bond that connects a stereo
site to one of its atom ligands, since topology removal cascades a stereo entity
only when a participant is removed.

Measured over 2,000 draws each: 3.8% of materializable reactions carry lhs
constraints and 5.2% constraint deltas; 2.4% of matchable reactions carry
constraint deltas and none lhs constraints; 5.4% of reaction spans carry
constraint spans. The graph-IR property suite runs in 37.5 s, against 33.8 s
before.

**Tests:** the existing aggregate properties on the extended generators.

#### S0a7 — Canonical equality contract properties **Done**

**Modules:** graph-IR property tests for molecule, reaction span, and reaction
canonicalization. **Additive, green.** [dep: S0a6]

For Molecule, ReactionSpan, and Reaction, test_<aggregate>_canonical_eq_contract
asserts that canonical_eq agrees with equality of the canonical forms, with two
contradictions equal. molecule_canonical_pair_strategy derives `b` from a
generated molecule by renumbering, normalization, reframing, setting one atom's
charge, or adding an atom. reaction_canonical_pair_strategy derives it from a
generated reaction by renumbering its span, normalization, reframing,
rematerializing its span, a no-op charge modification of one atom, or adding an
atom; reaction_span_canonical_pair_strategy takes the spans of those pairs.

Measured over 1,000 pairs each, as identical / canonically equal but different /
not equal: molecules 17% / 42% / 41%, reaction spans 62% / 16% / 22%, reactions
35% / 35% / 30%. The properties pass at 256 and at 2,048 cases.

#### S0b — Canonicalization helpers **Done**

**Modules:** graph-IR canonicalization and its tests. **Green.**
[dep: none]

Replace the private canonicalization helper layers with two private methods per
aggregate, defined in the canonicalization module:

```rust
impl Molecule {        // likewise ReactionSpan and Reaction
    fn tracked_canonicalize_at(self, level: DescriptionLevel, context: &CanonicalizeContext)
        -> Result<(Self, MoleculeRemapping), Self::Error>;
    fn canonical_key_at(&self, level: DescriptionLevel, context: &CanonicalizeContext)
        -> Result<CanonicalComparisonKey, Self::Error>;
}
```

The trait methods compute the level as today and call these methods.
canonicalize drops the remapping from tracked_canonicalize_at;
tracked_canonicalize returns it; canonical_eq compares canonical_key_at
results. Error conversion uses `?`. Users still cannot choose the level.

- Molecule: tracked_canonicalize_at matches on the level and calls the level
  workers canonicalize_topology, canonicalize_constitution,
  canonicalize_structure, and canonicalize_full. These are today's
  `_with_options` functions under the bare names; they take the owned molecule
  where they normalize it. The Full level normalizes the owned molecule and no
  longer clones it. canonical_key_at absorbs canonical_key_by and
  canonical_key_by_full.
- ReactionSpan: tracked_canonicalize_at normalizes the owned span instead of a
  clone. canonical_key_at replaces canonical_reaction_span_key.
- Reaction: tracked_canonicalize_at converts to a span, calls the span method,
  and converts back. canonical_key_at converts to a span and returns the span's
  key, so Reaction canonical_eq compares span keys instead of canonicalizing
  both reactions.
- Delete the forwarding and conversion functions:
  canonicalize_molecule_by_effective,
  canonicalize_molecule_with_remapping_by_effective,
  canonicalize_reaction_span_by, canonicalize_checked_reaction_span_by,
  canonicalize_checked_reaction_span_with_remapping_by,
  canonicalize_reaction_span_by_effective,
  canonicalize_reaction_span_with_remapping_by_effective,
  canonicalize_reaction_by, canonicalize_reaction_by_effective, and
  canonicalize_reaction_with_remapping_by_effective. Delete the test-only
  per-level wrappers; tests call the level workers or the trait methods.
- Rename canonicalize_level_with_constraints to level_with_constraints.
- The span method's error converts into the reaction error through
  `impl From<ReactionSpanCanonicalizeError> for ReactionCanonicalizeError`.

**Tests:** existing canonicalization, tracked remapping, and canonical_eq
results are unchanged.

**Verification:** the canonicalize bench, against a baseline taken before the
change. Over the complete, scaling canonicalize and tracked_canonicalize, and
reaction canonical_eq cases, median changes were within ±3% except Reaction
canonical_eq, which is 16–18% faster at the median and up to 29% faster on small
reactions; single-case outliers did not reproduce on rerun.

#### S0c — Borrowing canonical_hash **Done**

**Modules:** graph-IR canonicalization; Python bindings that call it.
**Breaking, red→green within the subitem.** [dep: S0b]

Change the signature to `canonical_hash(&self, context) -> Result<u64,
Self::Error>`. Each aggregate hashes its canonical_key_at result, the value
canonical_eq compares. Equal comparison keys then imply equal hashes.

Update the Canonicalize trait listing in the data-type guide.

No Python binding calls canonical_hash.

**Tests:** hash equality for canonically equal inputs under remapping, and a
value and its canonical form hash alike; the canonical equality contract
properties also assert that the hashes agree with canonical_eq on their pairs,
with two errors equal.

**Verification:** the canonicalize bench canonical_hash cases, against a
baseline taken before the change: 4–7% faster, since the hash no longer builds
the canonical form.

#### S0d — Consuming tower implementations move their storage **Done**

**Modules:** graph-IR aromatic, dative, multicenter, noncovalent, stereo,
reaction. **Green.** [dep: none]

Molecule's six overlay collections and the six span types clone each entity's
attributes inside normalize and reframe_by:
`*set.data_mut(id) = set.data(id).clone().normalize()?`. Move the stored value
out with mem::take on data_mut, normalize or transport it, and write it back;
the forms implement Default. Reaction normalize passes its deltas to
normalize_reaction_deltas by reference and rebuilds them. Change that function
to consume the deltas. Results are unchanged.

The span types store `EntitySpan<form>`, so EntitySpan gains
`impl<T: Default> Default for EntitySpan<T>`, returning `Unchanged(T::default())`.
normalize_reaction_deltas also moves the attributes it transports. Its three
borrowing callers (to_reaction_span, apply_at, prepare_deltas) clone the deltas
they pass, where the function used to clone internally.

**Tests:** existing normalize and reframe results for each type and for
Reaction.

**Verification:** the reframe and canonicalize bench cases, against a baseline
taken before the change: reframe cases 0.5–5.6% faster, the canonicalize and
reaction canonical_eq cases unchanged within ±4%.

### S1 — Combination, constraints, and components

#### S1a — Batch combine_from **Done**

**Modules:** graph-core id, graph, and relation sets; graph-IR molecule, editor,
and overlay collections; the graph-core and reaction benches. **Breaking,
red→green within the subitem.** [dep: none]

Rewrite combine_from to append other's entities with add_atoms, add_bonds, and
the overlay batch methods, so that one combine extends the graph once. Ids,
order, attributes, and constraints are unchanged.

The batch-append methods returned their appended ids as an opaque
`impl ExactSizeIterator`, which inherits Iterator's `must_use` although the
mutation has already happened; every caller that does not need the ids had to
discard them. They return `IdRange<Id>` instead: the contiguous positions
`start..end`, iterating as ids through `Id: From<usize>`, with
Iterator, DoubleEndedIterator, ExactSizeIterator, and FusedIterator and no
`must_use`. This covers Graph::extend_nodes, extend_edges, and extend, the five
relation sets' extend, Molecule's and MoleculeEditor's batch add methods, and the
overlay collections' extend. The graph-core ids NodeId, EdgeId, and RelationId
move with IdRange into a new id module; the crate root re-exports the same
names. The discards of appended ids in graph-core, its benches, and graph-IR are
removed.

**Tests:** existing combine and combine_from results, including every overlay
kind and molecule-level constraints; IdRange iteration in both directions, its
length, and its start-after-end panic.

**Verification:** the reaction bench's new combine_from/path_pair cases, against
a baseline taken before the change: 36% faster on 8 atoms and 86% faster on 64
atoms (81.6 to 11.7 µs), within 10% of combine_all.

#### S1b — Consuming combine **Done**

**Modules:** graph-IR molecule and fragment; Python molecule binding; tests and
benches that call combine. **Breaking, red→green within the subitem.**
[dep: S1a]

Change the signature to `combine(self, other: &Molecule) -> Molecule`,
implemented as combine_from followed by returning the receiver. Fragment
assembly passes its owned bodies unchanged. Callers that keep the receiver write
`a.clone().combine(&b)`: the molecule tests, the meet_pushout and structure
property tests, the whitepaper test, and the reaction bench. The Python combine
binding clones the receiver and consumes the clone, so the Python source is
unchanged.

**Tests:** test_molecule_combine_receiver: combine of an owned and of a cloned
receiver gives identical results, and the original of a cloned receiver is
unchanged. The reaction bench's combine case clones its receiver in the batch
setup, outside the measurement.

#### S1c — Consuming Constraints::compact

**Modules:** graph-IR constraint/molecule and molecule. **Breaking, red→green
within the subitem.** [dep: none]

Change the signatures to `compact(self, compaction: &MoleculeCompaction) ->
Self` and `tracked_compact(self, compaction: &MoleculeCompaction) -> (Self,
CascadedConstraints)`. The Molecule callers compact through constraints_mut
with mem::take. restore is unchanged. Change Constraints::extend to take
`impl IntoIterator<Item = Constraint>`, as the per-entity containers do;
callers passing a Vec are unchanged.

**Tests:** existing compaction, tracked-compaction, and restoration results.

#### S1d — Consuming constraint placement

**Modules:** graph-IR molecule and its tests. **Breaking, red→green within the
subitem.** [dep: none]

Change the signatures to `lift_constraints(self) -> Self` and
`inline_constraints(self) -> Result<Self, Contradiction>`. inline takes the
molecule-level list, meets each entry into its entity store, and returns the
contradiction on failure; it no longer clones the constraints into a plan.
Migrate the test callers in molecule/tests and tests/property/edit.

**Tests:** existing lift, inline, and round-trip results, including collisions
and contradictions.

#### S1e — Batch component building in split

**Modules:** graph-IR molecule. **Green.** [dep: none]

tracked_split builds each component with one add_atom, add_bond, or overlay add
per entity, and each add_bond extends the component's graph. Build each
component with the editor batch methods add_atoms, add_bonds, and the overlay
batch methods, so that each component extends its graph once. Components,
their order, and the returned correspondences are unchanged.

**Tests:** existing split and tracked_split results, including every overlay
kind.

**Verification:** the split case in the reaction bench, against a baseline
taken before the change.

### S2 — Transport

#### S2a — Consuming entity-collection transport

**Modules:** graph-IR aromatic, dative, multicenter, noncovalent, stereo.
**Breaking, red→green within the subitem.** [dep: none]

Change map and try_map, and the crate-private remap, on Molecule's six overlay
collections and the six span types to take `self`. Callers that keep the
collection clone it. S4d touches the same span types; either order works.

**Tests:** existing transport results for each type.

#### S2b — Consuming Molecule and ReactionSpan remap

**Modules:** graph-IR molecule/remap, reaction_span, canonicalization, and
callers. **Breaking, red→green within the subitem.** [dep: S2a]

Change remap and try_remap on Molecule and ReactionSpan to take `self`. Permute
the owned tables instead of cloning them, and remap the overlay collections with
their consuming remap. Canonicalization's final remap consumes the owned
normalized value; callers that keep the source write `m.clone().remap(&r)`.

**Tests:** existing remapping identity, inverse, and composition properties.

**Verification:** the canonicalize bench, against a baseline taken before the
change.

#### S2c — MoleculeMetadata::map

**Modules:** graph-IR dsl/metadata; Python metadata binding.
**Breaking, red→green within the subitem.** [dep: none]

Rename MoleculeMetadata::remap to map, keeping its consuming signature, and
rename the Python method to match.

**Tests:** existing metadata remapping tests under the new name.

### S3 — Engines

#### S3a — Consuming add_systems

**Modules:** graph aromaticity and its tests. **Breaking, red→green within the
subitem.** [dep: none]

Change the signature to `add_systems(&self, molecule: Molecule, systems: Vec<(Vec<AtomId>,
AromaticSystemForm)>) -> Molecule`. Apply its edits through an editor instead of a
transaction. Migrate the test callers in the aromaticity tests and
tests/kekulization.

**Tests:** existing perception and system-addition results.

#### S3b — Resolver plan application as an invariant

**Modules:** graph resolve, the five individual resolvers, and their tests.
**Breaking, red→green within the subitem.** [dep: none]

Replace post-plan `apply(edits)?`, `finish()?`, and `transact([edits])?` in the
individual resolvers, and the composite resolver's phase-plan application, with
expect and an invariant message, as the transformers do. Remove the Apply
variants from the eight error types and the composite ResolveError mappings that
wrap them.

**Tests:** existing resolution and projection results; remove tests that
construct the removed variants.

#### S3c — Kekulizer validity during planning

**Modules:** graph transform/kekulizer and the kekulization tests. **Breaking,
red→green within the subitem.** [dep: none]

In planning, check for each matched and exposed atom that the proposed charge,
lone pairs, and unpaired electrons are compatible with its current values, and
report an incompatibility as the planning error
`KekulizeError::IncompatibleLocalizedAtom { system, atom }`. Remove validate_localized_candidate
and the PostLocalizationValenceInvariant, PostLocalizationValenceInvariantError,
PostLocalizationSpinInvariant, and PostLocalizationSpinInvariantError variants.
transform_into applies the plan in place without a transaction.

**Tests:** the cyclopentadienyl-anion case with a conflicting literal charge and
lone-pair count fails during planning with IncompatibleLocalizedAtom; existing
kekulization results are unchanged.

#### S3d — In-place transform_into without a transaction

**Modules:** graph transform/aromatizer and transform/delocalize_charge.
**Green.** [dep: none]

Aromatizer and DelocalizeCharge transform_into apply their plan in place instead
of through transact. Planning failures leave the receiver unchanged.

**Tests:** existing transform_into results, and an unchanged receiver on
Aromatizer planning failure.

#### S3e — Individual resolver _into without a transaction

**Modules:** graph resolve isotope, aromaticity, stereo, bonds, multicenter.
**Green.** [dep: S3b]

The eight resolve_into and project_into methods apply their plan in place
instead of through transact. Outcomes decided during planning leave the receiver
unchanged.

**Tests:** existing _into results, and an unchanged receiver on
underdetermination and contradiction.

#### S3f — Transformer rule in the data-type guide

**Modules:** docs/development/data-types.md. **Documentation.** [dep: S3c, S3d, S3e]

Replace the transformer paragraph, which says transform_into uses a transaction
and restores on error, with the rule: a transformer derives a complete plan from
its input; every rejection is reported during planning, before mutation; a plan
applies to the molecule it was derived from, so a failed application is a defect;
transform consumes and applies the plan to its input, and transform_into applies
it in place without a transaction or rollback method. State the same rule for
the individual resolvers' resolve_into and project_into. Check that the
nomenclature guide's Transformation entry still agrees.

### S4 — Copy-on-write storage

#### S4a — Molecule constraints

**Modules:** graph-IR molecule, molecule/editor, molecule/pushout,
molecule/remap, molecule/tests. **Green.** [dep: none]

Change the field to `Arc<Constraints>`. constraints_mut uses Arc::make_mut, and
every mutating call (take, push, extend, remove_at, compact, tracked_compact,
restore) goes through constraints_mut. normalize uses Arc::make_mut followed by
mem::take on the inner Constraints. reframe_by and reframe use
Arc::unwrap_or_clone. Wrap at construction in Molecule::from_entries, remapping,
and pushout. MoleculeEditor::from_parts takes `Arc<Constraints>`. Readers that
need an owned list call constraints().clone().

**Tests:** migrate the struct literals in molecule/tests that set the field.
Add clone-independence cases: after cloning a molecule with molecule-level
constraints, a constraint edit, compaction, transaction rollback, and normalize
on one copy leave the other unchanged.

#### S4b — Reaction deltas

**Modules:** graph-IR reaction. **Green.** [dep: none]

Change the field to `Arc<Deltas>`. Reaction::try_new wraps the supplied Deltas.
normalize wraps its rebuilt Deltas. reframe_by, reframe, and into_parts use
Arc::unwrap_or_clone. The application iterator keeps its own Deltas field.

**Tests:** migrate the test literals and assignments that set the field.
Add clone-independence cases through normalize and reframe, and into_parts of a
shared and of a uniquely owned reaction.

#### S4c — ReactionSpan tables

**Modules:** graph-IR reaction_span. **Green.** [dep: none]

Change atoms, bonds, and constraints to `Arc<Vec<…>>`. Wrap at construction in
from_entries, try_remap, and reverse_sides. normalize, reframe_by, and reframe
rebuild each table with Arc::unwrap_or_clone. entries and try_remap clone the
inner Vec where they need an owned table.

**Tests:** clone-independence cases through normalize and reframe.

#### S4d — ReactionSpan overlay span newtypes

**Modules:** graph-IR dative, aromatic, multicenter, noncovalent, stereo.
**Green.** [dep: none]

Wrap the relation set of DativeBondSpans, AromaticSystemSpans,
MulticenterBondSpans, NoncovalentBondSpans, StereoAtomSpans, and
StereoBondSpans in an Arc. new, try_map, and remap wrap their result.
into_entries uses Arc::unwrap_or_clone. The Normalize, FrameTransport, and
reframe implementations write through Arc::make_mut.

**Tests:** migrate the unit tests that permute participants directly. Add a
clone-independence case per type through normalize and reframe.

#### S4e — Benchmark cases

**Modules:** graph-IR benches. **Additive, green.** [dep: S4a–S4d]

Add the aggregate_storage bench target with the twelve measured cases.


### S5 — Development guides

#### S5a — Ownership rules in the guides

**Modules:** docs/development/data-types.md and nomenclature.md.
**Documentation.** [dep: S0–S4]

Add to the data-type guide:

- under the operation taxonomy, the ownership rule: the bare verb consumes;
  `_into` mutates in place with recovery; a caller that keeps its input clones
  it; a consuming operation uses its ownership; nothing is handed back on
  failure; logical analogs share one contract, namely the normal-form tower and
  the transport operations map, remap, compact, and reframe_by;
- under raw representation and semantic normal form, that the equalities borrow,
  a default may clone, and implementations override with borrowing computations;
- the two-operand rule: borrow both operands, except an infallible operation
  whose result extends its receiver, which consumes it; in-place counterparts
  are named `_from` or `_with`;
- under Molecule mutation and publication, that every fallible operation on a
  published Molecule is atomic and in-place editing exists only with recovery;
- under representation ownership, the copy-on-write storage of the aggregates:
  a clone increments reference counts, writes go through Arc::make_mut, and
  consuming rebuilds use Arc::unwrap_or_clone.

Add to the nomenclature guide an `_into` suffix entry and an update entry: a
form's update applies an update record and returns a new form; a container's
update merges a same-type container in place by key. Note in the Combine entry
that combine consumes its receiver.

#### S5b — Whitepaper code against the library

**Modules:** whitepaper/main.tex and umol-graph/tests/whitepaper.rs.
**Documentation.** [dep: S0–S4]

Check every code listing and inline code reference in the whitepaper against
the library after this plan: names, signatures, ownership, and results. Update
the listings that differ, and keep umol-graph/tests/whitepaper.rs compiling the
Rust-primer appendix listings exactly as printed.

**Tests:** the whitepaper test.

### S6 — Verification

#### S6a — Final gate

**Modules:** workspace. **Verification.** [dep: S0–S5]

Run workspace tests and lints with all features, and the Python suite in the
repository's Python 3.13 environment. Rerun benchmarks only for cases whose
production code differs from what was measured.

**Critical path:** within S0, S0b before S0c; within S1, S1a before S1b; within S2, S2a before S2b; S4e after S4a–S4d; S5a after S0–S4; S6a last. The stages S0–S4 are otherwise independent of each other.
