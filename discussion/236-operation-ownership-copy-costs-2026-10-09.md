# 236 — Operation ownership and copy costs

Status: Proposed
Date: 2026-10-09
Relates: [181](181-python-boundary-ownership-2026-08-03.md),
[168](168-api-hygiene-2026-07-27.md),
[213](213-editor-overlay-storage-2026-08-27.md)

## Scope

This document owns the ownership contracts of graph-IR and graph operations:
which operations consume, borrow, or mutate their input, how they are named, and
the copy costs those contracts imply. It also owns the planning of resolvers and
transformers, with the shared aromaticity machinery, the constraint readings and
discharge, and copy-on-write storage for the three graph-IR aggregates.

[181](181-python-boundary-ownership-2026-08-03.md) remains the plan of record
for the molecule lifecycle, recovery, and Python ownership. Its plan follows the
decisions recorded here.

The operations are reviewed in eight steps: leaf forms and the normal-form
tower; two-operand operations; constraint containers; entries, frame actions,
and batches; Molecule as a value; Molecule as a published object; resolvers
and transformers; and Python.

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

## Decided: resolver and transformer naming

The bare verb consumes its molecule. `_into` mutates a borrowed molecule with
recovery. A caller that keeps its molecule clones it explicitly. This covers the
composite resolve, resolve_with_report, and project, the individual resolvers'
resolve and project, Transformer::transform, and their `_into` counterparts.
Python binds the bare verb as clone followed by the Rust method, so the Python
source survives every outcome.

## Decided: aromatic system names

Identifiers name an aromatic system in full, as `aromatic_system` and
`aromatic_systems`, also where aromaticity perception supplies the context. A
module name or a very long compound may shorten the name, and then to
`aromatic`, never to `system`. This plan applies the rule to the aromaticity
perception, resolver, and transformer code it reworks; the rest of the workspace
is separate work.

## Decided: planning in resolvers and transformers

Planned edits apply to the molecule they were planned from. A plan that fails to
apply is a planner defect, not an operation outcome. The five individual
resolvers' application of their own plans, and the composite resolver's
application of phase plans, therefore become plan invariants. The Apply variants
of AromaticityError, AromaticityProjectError, BondsError, IsotopeError,
IsotopeProjectError, MulticenterBondsError, StereoError, and StereoProjectError
are produced only there and are removed.

Every failure of an individual resolver or transformer occurs during planning. A
planner detects during planning that its plan would violate integrity, and fails
then. Property tests check, over generated inputs, that every plan these
operations determine applies to its input and that the result passes the
integrity check. The data-type guide states the rule.

The aromatic and stereo mismatch policies are removed: a disagreement between an
assertion and an entity is always a contradiction. AromaticityMismatchPolicy,
AromaticBondConstraintMismatchPolicy, and StereoMismatchPolicy go with their
configuration fields and Python bindings. The failure policies stay.
AromaticityResolveConfig keeps `perception: AromaticityConfig`,
`aromatic_valence_failure: AromaticityFailurePolicy`, and
`aromatic_system_failure: AromaticityFailurePolicy`, and gains
`atom_selection_algorithm: AtomSelectionAlgorithm`. StereoResolveConfig keeps
`tetrahedral_stereo_failure`, `stereo_atom_failure`, `cis_trans_stereo_failure`,
and `stereo_bond_failure`, each a StereoFailurePolicy.

The only integrity violation an aromaticity plan can produce is
AromaticSystemsOverlap. derive's stored-system rule excludes it: a stored system
that no perceived system reproduces with equal members is an
AromaticSystemFailure, and so is every stored system that a perceived system
overlaps without equal members. The selection applies `aromatic_system_failure`
to it, and under Keep the stored system's component gains no system. The
aromatizer plans nothing when systems exist, and the kekulizer and projection
only remove systems.

The aromatic perception machinery is shared and crate-private, and no individual
operation duplicates it. The perceiver dispatches to its three rules, and its
methods give what each rule supplies: its candidate rings, its perception from
the atoms' contributions, its components, its candidate systems, and its
acceptance test. `candidate_components(&self, molecule: &Molecule, rings:
&RingSet, can_contribute: impl Fn(AtomId) -> bool) -> Vec<BTreeSet<AtomId>>`
returns the components that the selection searches one at a time: for Hückel and
Clar the connected unions of rings whose atoms can contribute, for HMO the π
components that contain a ring atom. `candidate_aromatic_systems(&self,
molecule: &Molecule, rings: &RingSet, atom: AtomId, can_contribute: impl
Fn(AtomId) -> bool) -> Vec<BTreeSet<AtomId>>` returns the atom sets the rule
could perceive as one system containing `atom`: Hückel the rings and ring unions
that contain it, Clar the rings that contain it, and HMO the atom's π component
over the atoms that can still contribute, when that component has three or more
atoms and a ring atom. An atom can still contribute when some choice that
remains gives it an aromatic contribution. `accepts_contributions(&self,
members: &[&[AromaticValence]]) -> bool` returns false only when no choice of
the members' possible contributions is accepted by the rule: Hückel tests
whether some choice of an aromatic contribution per member totals 4n+2, and HMO
and Clar, which supply no such bound, return true.

Perception takes each atom's contribution through `impl Fn(AtomId) ->
AromaticValence`: NotAromatic excludes the atom, and Aromatic(n) contributes n.
Completions are always literal. An atom without completions whose
`evidenced_complete` aromatic valence is aromatic with a non-literal count makes
the result Underdetermined before perception; every other atom's contribution is
that reading's literal, through `as_lit`. No contribution is lowered to a bare
integer. derive is the policy-free core of every aromaticity path: for given
contributions it perceives the systems and classifies the inconsistencies,
including the stored-system rule and the failure checks of the next section,
once. The selection searches the completions and applies the failure policies to
the derivation of each combination. The standalone resolver runs the selection
with no completions, and the conformance validator reports the derivation's
first inconsistency, as stereo perception's derive serves the stereo resolver
and validator. `AromaticityPerceiver::plan_aromatic_system(&self, molecule:
&Molecule, atoms: Vec<AtomId>, form: AromaticSystemForm) -> Edits` adds a system
and discharges its members' aromatic valence; the standalone resolver, the
composite resolver, and the aromatizer use it.
AromaticityPerceiver::add_systems, which no production code calls, is removed.

The perceiver, its derivation, the three rule types with their errors and
output, and electrons_from_kekule are crate-private. So are the composite
resolver's phases and the state they pass: the selection, select and
select_with_report; the admit methods of ValenceResolver, AtomTypingValence, and
CountsValence; and ResolveState and AtomCompletions. Outside the crate only
tests use them, and those tests move to the public operations. Outside the
crate, aromaticity is reached through the resolver's resolve and project and
their plans, the aromatizer, the kekulizer, and the conformance validator.

The aromatizer reads `evidenced`, not `evidenced_complete`: its input's aromatic
facts are by purpose unstated, and the closure would read every unmarked atom as
non-aromatic. It runs derive with each atom's Kekulé contribution, and an atom
whose evidence is non-aromatic is excluded from perception. A member of a
perceived system whose evidence gives another contribution is an
AromaticValenceMismatch, identified by the id that system receives when the plan
applies. An atom that asserts aromaticity and that no perceived system claims is
an AromaticValenceFailure, and a bond that asserts aromaticity outside every
perceived system is an AromaticBondConstraintFailure. Transformers have no
policies, so each of these is a contradiction. Bonds that assert non-aromaticity
restrict the aromatizer's candidates as they restrict the resolvers'.

Kekulizer establishes validity during planning. A valid aromatic system admits a
matching of the size its electron count and charge require; planning already
checks this. Planning also builds each system atom's planned form: the atom as
the plan leaves it, with its planned charge and lone pairs and without the
aromatic valence assertion the plan removes, and with its incidence values met
into its constraints, namely the localized valence of its planned bonds, a
non-aromatic valence, and its unchanged multicenter valence and dative pairs. A
meet to ⊥, such as an asserted `#v` that the planned bonds change, is
`KekulizeError::ContradictedAssertion { aromatic_system, atom, constraint:
AtomConstraintForm }`. Otherwise ValenceInvariantsValidator::validate_atom
decides electron conservation from the planned form, and a mismatch is
`KekulizeError::ValenceInvariant { aromatic_system, mismatch: ValenceMismatch
}`. Kekulization never changes unpaired electrons, so the spin half of the
whole-molecule validation has no planning counterpart. The whole-molecule
validation after publication, and its four PostLocalization error variants, are
removed.

The selection's completion search keeps its pruning, which highly connected
systems require. Its enumeration bound counts the perception calls the search
makes, not the product of the flexible atoms' options, so the pruning extends
the systems the search handles; a search that exceeds the bound makes the result
Underdetermined. An aromatic-only atom is an atom whose every completion asserts
an aromatic valence. The pruning cuts a subtree when, for some aromatic-only
atom, every candidate system that contains it fails the rule's acceptance test
for every choice left below it. The cut applies under the Error
`aromatic_valence_failure` policy, under which an unclaimed aromatic-only atom
invalidates the combination; under Keep, the pruned search visits what the full
enumeration visits. `AtomSelectionAlgorithm::{FullDfs, BacktrackingDfs}` chooses
between the full enumeration and the pruned search; it is the field
`atom_selection_algorithm` of AromaticityResolveConfig and defaults to
BacktrackingDfs. Cross-checks compare the two algorithms through
Resolver::resolve_with_report, where both stay within the bound, and the
test-side copy of the selection goes.

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

## Decided: constraint readings and discharge

Every constraints view has six readings of a key: `asserted`, `derived`, and
`evidenced`, and their complete variants. Each reading is the meet of its
sources, starting from top. `asserted` meets the entity's own stored value,
`derived` the value derived from the relations the entity takes part in, and
`evidenced` all the evidence the molecule holds: the own assertion, the derived
value, and, for an atom's aromatic valence, each incident bond's assertion of
aromaticity, read as `#a+`. Each reading returns `Option<F>`: `None` is ⊥, as in
`meet`, and the key's vacuous form is no information, so every value has one
encoding. Top met with a value is that value's normal form, and a vacuous entry
and an absent one both give top, so every reading returns the same value on a
molecule and on its normal form. On a molecule that normalizes, only `evidenced`
reaches ⊥, from conflicting evidence.

A complete variant replaces no information with the key's definite negative
where the key has one. For an assertion, no information is a top value: a
vacuous entry states nothing. For a relation, it is the absence of the overlay:
an overlay entity states that it exists even when its value is open, as a stereo
atom of undetermined kind, whose derived tetrahedral stereo is top, does not
close to `#T!`. `evidenced_complete` returns `evidenced` where it is not top and
`derived_complete` otherwise. Its claim is that the molecule states every
positive overlay fact, as an assertion or as a relation, without either source
alone being complete; `asserted_complete` claims that the assertions are
complete, and `derived_complete` that the relations are. Topology keys always
have a derived value, so the readings agree on them. Every reading that derives
a ring key requires ring context.

Each constraint form type gains `undetermined(key)`, its vacuous form at a key,
which a reading returns for a key without a source. The containers' typed
getters return the stored value, owned, and the vacuous value for an absent key,
as `aromatic()`, `electron_count()`, `intramolecular()`, and `stereogenicity()`
already do. `get`, `contains`, and `iter` remain the storage accessors, and the
views' inherited getters follow the containers'.

`asserted_complete` reads only the entity's own storage and the absence cell. An
incident bond's assertion of aromaticity moves into the evidenced readings,
where it is met with the atom's own assertion instead of yielding to it: an atom
that asserts non-aromaticity next to a bond that asserts aromaticity reads ⊥,
not non-aromatic. The constraints views' `is_compatible` compares against
`evidenced`.

Valence admission, under atom typing and under counts, reads
`evidenced_complete`, and so do aromaticity perception, for each atom's
contribution, and stereo perception. Admission builds the ring set when a
registry row states a ring key. An atom whose evidence meets to ⊥ admits no
candidate, so atom typing and counts report their existing NoMatch. derive
reports an atom whose aromatic evidence meets to ⊥ as an AromaticValenceMismatch
when the conflict is with its system, and as an AromaticValenceBondMismatch when
the conflict is with an incident bond. Three outcomes change. Counts reads an
atom that asserts non-aromaticity inside a stored system as ⊥, not aromatic.
Perception reads an atom that asserts m electrons inside a stored system that
gives it n ≠ m as ⊥, not m. An atom without aromatic evidence outside every
system is excluded from perception instead of contributing zero electrons.
Stereo perception's outcomes are unchanged.

Substructure matching, and reaction matching through it, asks the constraints
views' `satisfies` whether a host entity meets a pattern entity's constraints.
`satisfies` states that logical requirement and answers it algebraically: every
pattern entry must be satisfied by the host's `evidenced_complete` value at its
key, and ⊥ satisfies no entry. On a host whose relations realize its positive
assertions, the result is unchanged. A host atom whose aromatic assertion no
system realizes, or whose only aromatic evidence is an incident bond, now
matches an aromatic pattern entry.

TableIR raising writes `#a!` on an atom only when the table marks the atom
non-aromatic and one of its bonds aromatic. Elsewhere, the absence of aromatic
evidence already reads non-aromatic under `evidenced_complete`. Raising writes
no `#m!`.

The composite resolver has no closing discharge pass. Discharge is part of each
resolver's plan: the plan removes the assertions that its own results realize.
Valence discharges `#v` on each atom it completes. Aromaticity discharges the
aromatic valence of each member of a system it adds, and the aromatic constraint
of each bond outside the stored systems that its result confirms; it writes no
aromatic constraint on bonds. Stereo discharges the tetrahedral or cis-trans
constraint at each site where it adds an entity. `reset_aromatic_valence` and
`reset_stereo_constraints` are removed. Assertions on entities that resolution
neither adds nor completes stay. Resolution no longer evaluates the other
assertions or the molecule-level list; the validators report a violated one.

Each check that the closing pass made for a discharged key moves into the
planning of the resolver that discharges the key. Valence admission already
checks `#v`, and stereo planning already applies the stereo failure policies and
reports each mismatch as a contradiction. Aromaticity planning gains two failure
checks in derive and one restriction in perception. An atom that asserts
aromaticity and that no system claims is an AromaticValenceFailure. A bond that
asserts aromaticity outside every system is an AromaticBondConstraintFailure.
Both follow `aromatic_valence_failure`. A bond that asserts non-aromaticity
restricts perception: no system the plan adds contains it, because a candidate
that would contain it is not considered. Under Error, such a bond leads to a
contradiction only through the failure checks, when no valid candidate remains.
A stored system that contains such a bond remains an
AromaticBondConstraintMismatch. AromaticityInconsistency gains
`AromaticBondConstraintFailure { bond: BondId }` and
`AromaticValenceBondMismatch { atom: AtomId, bond: BondId }`.

Projection is unchanged. The aromaticity and stereo projection plans meet each
stored assertion with the projected value and fail during planning; they use no
policy and no reset, and the composite resolver's projection runs them without a
discharge pass. Resolving a projected molecule consumes the projected
assertions: aromaticity discharges them on each system it adds, and stereo at
each site where it adds an entity, so resolution after projection leaves none
behind. Export output is unchanged. Matching on a pattern's constraint entries
is invariant under projection: under `evidenced_complete`, an atom or bond of
the host has the same aromatic and stereo values before projection, from its
entities, and after, from the projected assertions. A pattern that contains an
aromatic system or a stereo entity as an entity still needs that entity in the
host.

## Decided: invariants

These changes establish the following invariants. Each is a property in the
plan.

1. A plan determined by an individual resolver or transformer applies to its
   input, and the result passes the integrity check.
2. Aromaticity resolution, standalone or composite, never produces overlapping
   aromatic systems; a perceived system that overlaps a stored one makes the
   stored system an AromaticSystemFailure during planning.
3. Admission, aromaticity perception for each atom's contribution, stereo
   perception, and matching on the host side read every derivable constraint
   they use through `evidenced_complete`, so they agree on its value.
4. An incident bond's assertion of aromaticity is the same evidence as the
   atom's own `#a+`; conflicting evidence is ⊥.
5. Resolution discharges the assertions its results realize: the aromatic
   valence of each member of a system it adds, the aromatic constraint of each
   bond outside the stored systems that its result confirms, the tetrahedral or
   cis-trans constraint at each site where it adds an entity, and `#v` on each
   atom it completes.
6. Every other assertion, and the molecule-level list, leaves resolution
   unchanged.
7. A SMILES or MOL input that resolves under the Error policies carries no
   aromatic or stereo assertion afterward.
8. Resolving a projected molecule leaves no projected assertion and, under the
   same model, reproduces the projected systems.
9. A host and its projection satisfy the same pattern constraint entries.
10. `AtomSelectionAlgorithm::FullDfs` and `BacktrackingDfs` give the same
    resolution wherever both stay within the bound.
11. No system that the resolver or the aromatizer adds contains a bond that
    asserts non-aromaticity or an atom whose evidence disagrees with its
    contribution.
12. Every atom the kekulizer changes satisfies electron conservation and its own
    assertions.
13. Every constraint reading returns the same value on a molecule and on its
    normal form.

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

Calling resolvers and transformers with the molecule as receiver, so that their
calls chain left to right, needs a separate design.

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

#### S1c — Consuming Constraints::compact **Done**

**Modules:** graph-IR constraint/molecule and molecule. **Breaking, red→green
within the subitem.** [dep: none]

Change the signatures to `compact(self, compaction: &MoleculeCompaction) ->
Self` and `tracked_compact(self, compaction: &MoleculeCompaction) -> (Self,
CascadedConstraints)`. The Molecule callers compact through constraints_mut
with mem::take. restore is unchanged. Change Constraints::extend to take
`impl IntoIterator<Item = Constraint>`, as the per-entity containers do;
callers passing a Vec are unchanged.

**Tests:** existing compaction, tracked-compaction, and restoration results.

#### S1d — Consuming constraint placement **Done**

**Modules:** graph-IR molecule and its tests. **Breaking, red→green within the
subitem.** [dep: none]

Change the signatures to `lift_constraints(self) -> Self` and
`inline_constraints(self) -> Result<Self, Contradiction>`. inline takes the
molecule-level list, meets each entry into its entity store, and returns the
contradiction on failure; it no longer clones the constraints into a plan.
Migrate the test callers in molecule/tests and tests/property/edit.

**Tests:** existing lift, inline, and round-trip results, including collisions
and contradictions. Two list entries with the same key on one entity combine
by meet: test_molecule_inline_constraints_meets_list_collision and
test_molecule_inline_constraints_contradiction_on_list_collision replace the
last-wins test, which asserted that the later entry overwrote the earlier.

#### S1e — Batch component building in split **Done**

**Modules:** graph-IR molecule. **Green.** [dep: none]

tracked_split builds each component with one add_atom, add_bond, or overlay add
per entity, and each add_bond extends the component's graph. Build each
component with the editor batch methods add_atoms, add_bonds, and the overlay
batch methods, so that each component extends its graph once. Components,
their order, and the returned correspondences are unchanged.

**Tests:** existing split and tracked_split results, including every overlay
kind.

**Verification:** the split case in the reaction bench, against a baseline
taken before the change: split and tracked_split are 60% faster at size 64 and
unchanged at size 8.

#### S1e1 — Mixed topology batch add **Done**

**Modules:** graph-IR molecule and editor. **Green.** [dep: S1e]

Add `add_topology(&mut self, atoms: Vec<AtomForm>, bonds: Vec<([AtomId; 2],
BondForm)>) -> (IdRange<AtomId>, IdRange<BondId>)` to Molecule and the
molecule editor, over Graph::extend. Bond endpoints refer to the resulting
atom space, including the appended atoms. add_atoms and add_bonds delegate to
it. combine_from and tracked_split add their atoms and bonds with one
add_topology call, so the separate node extension no longer copies a shared
adjacency that the edge rebuild then replaces. tracked_split records each
atom's id within its component during the partition, so bond endpoints are
known before the call.

**Tests:** test_molecule_editor_add_topology, for bonds between appended and
existing atoms on owned and shared sources, and
test_molecule_editor_add_topology_identity, for empty input; existing combine
and split results.

**Verification:** the combine and split cases in the reaction bench, against a
baseline taken before the change: combine and combine_from are 7–10% faster at
size 8 and 3% faster at size 64; split and tracked_split are 3–4% faster at
both sizes.

### S2 — Transport

#### S2a — Consuming entity-collection transport **Done**

**Modules:** graph-IR aromatic, dative, multicenter, noncovalent, stereo.
**Breaking, red→green within the subitem.** [dep: none]

Change map and try_map, and the crate-private remap, on Molecule's six overlay
collections and the six span types to take `self`. Callers that keep the
collection clone it. S5d touches the same span types; either order works.

**Tests:** existing transport results for each type.

#### S2b — Consuming Molecule and ReactionSpan remap **Done**

**Modules:** graph-IR molecule/remap, reaction_span, canonicalization, and
callers. **Breaking, red→green within the subitem.** [dep: S2a]

Change remap and try_remap on Molecule and ReactionSpan to take `self`. Permute
the owned tables instead of cloning them, and remap the overlay collections with
their consuming remap. Canonicalization's final remap consumes the owned
normalized value; callers that keep the source write `m.clone().remap(&r)`.
The topology, constitution, and structure searches take the molecule by value.
The full-level search remaps the normalized molecule with the remapping it
returns, instead of rebuilding the complete candidate for the selected order.
The reaction-span search remaps the normalized span with the remapping it
already computed, and reaction_span_from_candidate is removed. The per-leaf
remaps in complete_candidate and the reaction-span leaf candidate clone their
borrowed source.

**Tests:** existing remapping identity, inverse, and composition properties.

**Verification:** the canonicalize bench, against a baseline taken before the
change. The remapping cases are 1–9% faster, with the clone moved into batch
setup. Complete canonicalization is within ±4% for most of the corpus and
16–24% faster for overlay_heavy and frame_relative_stereo_constraint, which
skip the rebuilt complete candidate. The scaling canonicalize and
tracked_canonicalize cases are within ±2%.

#### S2c — MoleculeMetadata::map **Done**

**Modules:** graph-IR dsl/metadata; Python metadata binding.
**Breaking, red→green within the subitem.** [dep: none]

Rename MoleculeMetadata::remap to map, keeping its consuming signature, and
rename the Python method to match.

**Tests:** existing metadata remapping tests under the new name.

### S3 — Constraint readings and discharge

#### S3a — Aromatic skeleton generator **Done**

**Modules:** umol-graph tests/property/strategies.rs. **Additive, green.** [dep:
none]

Add a strategy for aromatic skeletons: single, fused, and coupled five- and
six-membered rings over C, N, O, and S; on each atom no aromatic evidence,
`#a+`, `#a<n>`, or `#a!`; on each bond none, `#a`, or `#a!`; undetermined or
literal hydrogen counts; optionally stored systems and stereo marks. Its
documentation states the domain it emits. Add the SMILES and MOL renderings of a
generated skeleton.

**Tests:** a coverage test that every evidence kind, topology, and option
occurs; properties that the SMILES and MOL renderings raise to the skeleton's
elements, hydrogen counts, and bonds.

#### S3a1 — Constraint forms at a key **Done**

**Modules:** graph-IR constraint/atom, bond, dative, aromatic, multicenter,
noncovalent, and stereo; umol-graph-ir tests/property. **Additive, green.**
[dep: none]

Add `undetermined(key)` to the eight constraint form types: the vacuous form at
`key`, keeping the ring scope, ligand pair, or ligand permutation that the key
carries.

**Tests:** per form type, the vacuous form at each key, compared with a literal;
a property over generated constraint forms of the eight types:
`undetermined(x.key())` has the key of x, is vacuous, equals
`x.as_undetermined()`, and meets x to its normal form, or to ⊥ where x does not
normalize.

#### S3a2 — Owned typed getters on constraint containers

**Modules:** graph-IR constraint/atom, bond, and dative, and view/constraints;
their callers in graph-IR, graph, io, and umol-py. **Breaking, red→green within
the subitem.** [dep: none]

The typed getters that return `Option<&T>` return `T`, the stored value or the
vacuous value for an absent key: the atom getters (`valence`,
`aromatic_valence`, `multicenter_valence`, `tetrahedral_stereo`, `degree`,
`total_degree`, `total_valence`, `ring_degree`, `ring_valence`,
`total_hydrogens`, `donated_pairs`, `accepted_pairs`, `ring_count`, and
`ring_size_count`), `cis_trans_stereo` on bonds, and `ring_count` and
`ring_size_count` on bonds and dative bonds. The views' inherited getters
follow. Callers match the value instead of an option, and matches that cover
both `None` and the vacuous value collapse.

**Tests:** per getter, a stored value and an absent key, which reads the vacuous
value; existing results.

**Verification:** the resolve and substructure bench groups, against a baseline
taken before the change.

#### S3b — Constraint readings as meets

**Modules:** graph-IR view/constraints and the entity view modules, and their
callers in graph-IR, graph, io, and umol-py; umol-graph-ir tests/property;
docs/development/nomenclature.md. **Breaking, red→green within the subitem.**
[dep: S3a1]

Give the eight constraints views the six readings of the constraint-readings
decision, each folding `meet` over its sources from `undetermined(key)` and
returning `Option<F>`, with `None` as ⊥ and the vacuous form as no information.
`asserted`, `asserted_complete`, `derived`, and `derived_complete` change from
`Option<&F>` and `Option<F>` with `None` for absence; `evidenced` and
`evidenced_complete` are added. `asserted_complete` changes only its encoding
and keeps the incident-bond merge, which S3g removes. Callers migrate.

In the nomenclature guide's Derived and asserted entry, define the six readings
as meets, the closure of each complete variant, and the claim behind it, and
list them under In code. Replace the sentence that each consumer reads exactly
one flat accessor: a pattern's own entries read `asserted`, resolution and
matching hosts read `evidenced_complete`, and aromatization reads `evidenced`.
Correct the doc comment of atom_derived_constraint, which says that a bond
asserting aromaticity yields `Aromatic(Undetermined)` on the derived side.

**Tests:** per view and reading, a source present, absent, and vacuous, and, for
`evidenced`, two sources agreeing and conflicting (`None`); each complete
variant closing where its source carries no information and keeping top for an
overlay entity with an open value, such as a stereo atom of undetermined kind; a
topology key; on atoms, aromatic valence evidenced by an incident bond's `#a`
alone, an own assertion met with it, a non-aromatic assertion next to it
(`None`), an incident bond's `#a!`, which is no evidence for the atom, and a
ring key with ring context. A property over generated molecules, for invariants
3, 4, and 13: `evidenced_complete` equals the corresponding flat reading or
their meet, an atom's own `#a+` and an incident bond's `#a` read alike, and
every reading returns the same value on a molecule and on its normal form.

**Verification:** the resolve and substructure bench groups, against a baseline
taken before the change.

#### S3c — Valence admission reads evidenced_complete

**Modules:** graph valence/atom_typing and valence/counts. **Green.** [dep: S3b]

Atom typing matches each registry-row entry against the atom's
`evidenced_complete` value, replacing the inline composition, and keeps a ground
atom under resolution when that reading's aromatic valence is `#a+`, replacing
the system-membership test. It builds the ring set when a row states a ring key.
Counts applies the same ground-atom test and reads aromatic evidence and the
candidate aromatic valences from `evidenced_complete`. An atom whose evidence
meets to ⊥ admits no candidate and reports NoMatch.

**Tests:** existing admission results; NoMatch under both sources for an atom
that asserts non-aromaticity next to a bond that asserts aromaticity, and under
counts for an atom that asserts non-aromaticity inside a stored system; a
registry row with a ring key admits by the atom's ring membership.

**Verification:** the resolve bench groups resolve_ingest_smiles and
smiles_roundtrip/resolve, against a baseline taken before the change.

#### S3d — Perception reads typed contributions from evidenced_complete

**Modules:** graph aromaticity, its rule modules, resolve/aromaticity, and
transform/aromatizer. **Green.** [dep: S3b]

Perception takes each atom's contribution through `impl Fn(AtomId) ->
AromaticValence` instead of `Fn(AtomId) -> Option<u8>`: find_systems, each
rule's find_from_rings, and the selection's contributions; electrons_from_kekule
returns `AromaticValence`, NotAromatic where the Kekulé structure gives no
contribution. An atom without completions whose evidence is aromatic with a
non-literal count makes the result Underdetermined before perception. derive and
the selection read each atom's contribution from `evidenced_complete` through
`as_lit`, replacing derive's contribution closure and stored_contribution.
Report ⊥ as AromaticValenceMismatch when the conflict is with the atom's system,
and as the new AromaticityInconsistency variant `AromaticValenceBondMismatch {
atom: AtomId, bond: BondId }` when it is with an incident bond.

**Tests:** existing derivation and selection results; six atoms asserting one
electron each and a methylene carbon without aromatic evidence, closing a
seven-membered ring, no longer form a system, so the six atoms are an
AromaticValenceFailure in standalone resolution; an atom that asserts
non-aromaticity next to a bond that asserts aromaticity is an
AromaticValenceBondMismatch in standalone resolution.

#### S3e — Stereo perception reads evidenced_complete

**Modules:** graph stereo. **Green.** [dep: S3b]

Read each site's tetrahedral or cis-trans evidence from `evidenced_complete`. A
stereo value derives or confirms the site's entity, ⊥ is the existing mismatch,
and the negative leaves the site without an entity.

**Tests:** existing stereo derivation and resolution results, including failures
and mismatches.

#### S3f — Matching reads evidenced_complete

**Modules:** graph-IR view/constraints and the substructure tests;
docs/development/nomenclature.md. **Green.** [dep: S3b, S3a]

Restate the constraints views' `satisfies` over `evidenced_complete`: every
pattern entry is satisfied by the host's value at its key, and ⊥ satisfies no
entry. Substructure matching keeps calling it. Restate the views'
`is_compatible` over `evidenced`. Restate `satisfies` in the nomenclature
guide's Constraints view entry and list the evidenced readings there.

**Tests:** existing substructure and reaction matching results; a host atom
whose aromatic assertion no system realizes, and one whose only aromatic
evidence is an incident bond, match an aromatic pattern entry; the views'
satisfies and is_compatible cases; a property over generated hosts, for
invariant 3: a single pattern entry matches a host entity exactly when the
entity's `evidenced_complete` value satisfies it; a property in umol-graph
tests/property/project.rs over hosts resolved from the skeleton generator, for
invariant 9: a host and its projection satisfy the same generated aromatic and
stereo pattern entries.

**Verification:** the substructure bench groups, against a baseline taken before
the change.

#### S3g — asserted_complete without the incident-bond merge

**Modules:** graph-IR view/atom and the view tests;
docs/development/nomenclature.md. **Green.** [dep: S3c]

The atom's `asserted_complete` closes `asserted` alone, where it is top, to the
absence cell, and no longer reads incident bonds. Move the guide's sentence on
the two aromatic placements to the evidenced readings.

**Tests:** an atom whose only aromatic evidence is an incident bond reads
non-aromatic under `asserted_complete`; the existing cases with own assertions.

#### S3h — Raising writes #a! only next to an aromatic bond

**Modules:** umol-io table_ir/raise and its tests; graph ingest tests.
**Green.** [dep: S3c, S3d]

Raise an atom that the table marks non-aromatic with `#a!` only when one of its
table bonds is aromatic.

**Tests:** raising `CCO` writes no aromatic valence; raising `c1ccccc1` writes
`#a+` on every atom and the aromatic constraint on every bond; raising an
aliphatic atom joined by an explicit aromatic bond writes `#a!` on it, and
ingesting that input reports NoMatch; resolved ingestion results are unchanged.

#### S3h1 — Raising writes #a! on a single bond between aromatic atoms

**Modules:** umol-io table_ir/raise and its tests; graph ingest tests.
**Green.** [dep: S3h, S3l]

Raise a single table bond with `#a!` when the table marks both of its atoms
aromatic. In SMILES this is an explicit single bond between aromatic atoms, such
as the bond joining the rings of `c1ccccc1-c1ccccc1`, which SMILES defines as
not aromatic.

**Tests:** raising `c1ccccc1-c1ccccc1` writes `#a!` on the bond joining the
rings and `#a` on every ring bond; raising `c1ccccc1C` and `CC` writes no
aromatic constraint on their single bonds; ingesting `c1ccc-cc1` fails with an
AromaticValenceFailure; other resolved ingestion results are unchanged.

#### S3i — Aromaticity discharges what its systems realize

**Modules:** graph resolve/aromaticity and resolve; umol-py resolve bindings and
Python tests; docs/development/nomenclature.md. **Breaking, red→green within the
subitem.** [dep: none]

plan_system removes the aromatic valence assertion of each member of the system
it adds and writes no aromatic constraint. The aromaticity plan removes the
aromatic constraint of each bond outside the stored systems that its result
confirms. Remove `reset_aromatic_valence` from AromaticityResolveConfig and from
the Python binding's keyword, getter, and repr. Remove its name from the
nomenclature guide's Reset entry.

**Tests:** standalone resolution removes those assertions for each added system
and keeps them for a stored system; `#a!` on the bond joining the two rings of
biphenyl is removed; composite results are unchanged.

#### S3j — Stereo discharges what its entities realize

**Modules:** graph resolve/stereo; umol-py resolve bindings and Python tests;
docs/development/nomenclature.md. **Breaking, red→green within the subitem.**
[dep: none]

The stereo plan removes the tetrahedral or cis-trans constraint at each site
where it adds an entity. Remove `reset_stereo_constraints` from
StereoResolveConfig and from the Python binding's keyword, getter, and repr.
Remove the guide's Reset entry with its last code names.

**Tests:** standalone and composite resolution remove the constraint at each
site that gains an entity, including one whose configuration is undetermined,
and keep it at a site with an existing entity.

#### S3k — Valence discharges #v

**Modules:** graph resolve. **Green.** [dep: none]

The composite resolver's constitution plan removes `#v` from each atom whose
fields it completes.

**Tests:** the constitution plan of an atom asserting `#v` that it completes
removes the assertion, and that of an atom it does not complete leaves it.

#### S3l — Aromaticity planning makes the closing pass's checks

**Modules:** graph aromaticity, its rule modules, and resolve/aromaticity.
**Green.** [dep: S3d, S3a]

In derive and in the selection, an atom that asserts aromaticity and that no
system claims is an AromaticValenceFailure, which the selection now also reports
for atoms without completions. A bond that asserts aromaticity outside every
system is the new AromaticityInconsistency variant
`AromaticBondConstraintFailure { bond: BondId }`. Both follow
`aromatic_valence_failure`. No rule's candidates and no selection combination
include a system that contains a bond asserting non-aromaticity. A stored system
that contains one remains an AromaticBondConstraintMismatch.

**Tests:** each check in standalone and composite resolution under the Error and
Keep failure policies; a benzene ring with one bond asserting non-aromaticity
forms no system, so its aromatic atoms are an AromaticValenceFailure under Error
and keep their assertions under Keep; naphthalene with such a bond on the
periphery of one ring forms a system on the other ring; a property over the
skeleton generator, for invariant 11: no system the resolver adds contains a
bond asserting non-aromaticity or an atom whose evidence disagrees with its
contribution.

#### S3m — Closing discharge pass removed

**Modules:** graph resolve and ingest tests; docs/development/nomenclature.md.
**Breaking, red→green within the subitem.** [dep: S3h, S3h1, S3i, S3j, S3k, S3l,
S3a]

Remove plan_discharge and its call sites, DischargeContradiction,
ResolveContradiction::Discharge, and ResolveError::DischargeEvaluation and
Discharge. Rewrite the nomenclature guide's Discharge entry: discharge is part
of each resolver's plan and removes the assertions that its own results realize.
Adjust the Concrete entry's reference to the closing pass.

**Tests:** composite resolution keeps assertions on entities it neither adds nor
completes and keeps the molecule-level list unevaluated; tests of the removed
variants go; properties over the skeleton generator: for invariants 5 and 6, the
resolved assertions in standalone and composite resolution are the input's
without the discharged ones, and the molecule-level list is unchanged; for
invariant 7, the renderings that resolve under the Error policies carry no
aromatic or stereo assertion; for invariant 8, resolving a projected molecule
leaves no projected assertion and reproduces its systems under the same model.

**Verification:** the resolve bench groups, against a baseline taken before the
change.

### S4 — Resolvers and transformers

#### S4a — Resolution and transform benchmark cases

**Modules:** umol-graph benches/resolve.rs, a new benches/transform.rs target,
and umol-graph Cargo.toml. **Additive, green.** [dep: none]

Add to the resolve bench a selection group over fused aza-aromatics whose
flexible atoms give option products from about 10² to beyond 4096, under the
Hückel, HMO, and Clar rules, and a group for AromaticityResolver::resolve over
the existing corpus. Add a transform bench target with the kekulizer over
aromatic inputs and the aromatizer over Kekulé inputs. Input construction and
cloning stay outside timing.

**Tests:** each case's setup asserts that its input resolves or transforms.

#### S4b — Shared system addition

**Modules:** graph aromaticity, resolve/aromaticity, resolve,
transform/aromatizer, and tests/kekulization. **Breaking, red→green within the
subitem.** [dep: S3i]

Move plan_system to `AromaticityPerceiver::plan_aromatic_system(&self, molecule:
&Molecule, atoms: Vec<AtomId>, form: AromaticSystemForm) -> Edits`, keeping the
discharge of S3i. The standalone resolver, the composite resolver's constitution
plan, and the aromatizer's plan use it; the aromatizer's inline copy goes, so it
writes no aromatic constraint on system bonds. Remove
AromaticityPerceiver::add_systems; tests/kekulization builds its fixture from
the DSL.

**Tests:** aromatizer results without aromatic constraints on system bonds;
existing resolution results.

#### S4c — Resolver plan application as an invariant

**Modules:** graph resolve, the five individual resolvers, and their tests.
**Breaking, red→green within the subitem.** [dep: S3a]

Replace post-plan `apply(edits)?`, `finish()?`, and `transact([edits])?` in the
individual resolvers, and the composite resolver's phase-plan application, with
expect and an invariant message, as the transformers do. Remove the Apply
variants from the eight error types and the composite ResolveError mappings that
wrap them.

**Tests:** existing resolution and projection results; remove tests that
construct the removed variants; a property over the skeleton generator, for
invariant 1, replacing the four fixed-SMILES integrity properties: every
determined result of the aromaticity resolver's resolve and project, the
kekulizer, and the aromatizer passes the integrity check, and the consuming and
in-place forms agree.

#### S4d — Kekulizer validity during planning

**Modules:** graph transform/kekulizer and the kekulization tests. **Breaking,
red→green within the subitem.** [dep: S3a]

In planning, build each system atom's planned form: the atom as the plan leaves
it, with its planned charge and lone pairs and without the aromatic valence
assertion the plan removes, and with the localized valence of its planned bonds,
a non-aromatic valence, and its unchanged multicenter valence and dative pairs
met into its constraints. Report a meet to ⊥ as
`KekulizeError::ContradictedAssertion { aromatic_system, atom, constraint:
AtomConstraintForm }`, and a mismatch from
ValenceInvariantsValidator::validate_atom on the planned form as
`KekulizeError::ValenceInvariant { aromatic_system, mismatch: ValenceMismatch
}`. Remove validate_localized_candidate and the
PostLocalizationValenceInvariant, PostLocalizationValenceInvariantError,
PostLocalizationSpinInvariant, and PostLocalizationSpinInvariantError variants.
transform_into applies the plan in place without a transaction.

**Tests:** the cyclopentadienyl-anion case with a conflicting literal charge and
lone-pair count fails during planning with ValenceInvariant, and a benzene
carbon asserting `#v2` with ContradictedAssertion; existing kekulization results
are unchanged; a property over resolved skeletons with `#v` and charge
assertions, for invariant 12: every atom the kekulizer changes satisfies
electron conservation and its own assertions, or kekulization fails during
planning with ContradictedAssertion or ValenceInvariant.

**Verification:** the transform bench's kekulizer cases, against a baseline
taken before the change.

#### S4e — In-place transform_into without a transaction

**Modules:** graph transform/aromatizer and transform/delocalize_charge.
**Green.** [dep: none]

Aromatizer and DelocalizeCharge transform_into apply their plan in place instead
of through transact. Planning failures leave the receiver unchanged.

**Tests:** existing transform_into results, and an unchanged receiver on
Aromatizer planning failure.

#### S4f — Individual resolver _into without a transaction

**Modules:** graph resolve isotope, aromaticity, stereo, bonds, multicenter.
**Green.** [dep: S4c]

The eight resolve_into and project_into methods apply their plan in place
instead of through transact. Outcomes decided during planning leave the receiver
unchanged.

**Tests:** existing _into results, and an unchanged receiver on
underdetermination and contradiction.

#### S4g — Mismatch policies removed

**Modules:** graph resolve/aromaticity, resolve/stereo, and resolve; umol-py
resolve bindings, the Python package exports, and Python tests. **Breaking,
red→green within the subitem.** [dep: none]

Remove AromaticityMismatchPolicy, AromaticBondConstraintMismatchPolicy, and
StereoMismatchPolicy, with the fields aromatic_valence_mismatch,
aromatic_bond_constraint_mismatch, tetrahedral_stereo_mismatch, and
cis_trans_stereo_mismatch. Every mismatch is a contradiction: the aromaticity
resolver's Keep, RemoveConstraint, and ReplaceEntity branches and the stereo
resolver's Keep, RemoveConstraint, ReplaceEntity, and RemoveBoth branches go.
Remove the Python classes, keywords, getters, and repr fields.

**Tests:** each mismatch is a contradiction in standalone and composite
resolution; tests of the removed branches go.

#### S4h — Rule-supplied components and candidate systems

**Modules:** graph aromaticity, its rule modules, and resolve/aromaticity.
**Green.** [dep: S3d]

Add `candidate_aromatic_systems(&self, molecule: &Molecule, rings: &RingSet,
atom: AtomId, can_contribute: impl Fn(AtomId) -> bool) -> Vec<BTreeSet<AtomId>>`
and `accepts_contributions(&self, members: &[&[AromaticValence]]) -> bool` to
the perceiver, dispatching to its rules. candidate_components becomes the
perceiver method `candidate_components(&self, molecule: &Molecule, rings:
&RingSet, can_contribute: impl Fn(AtomId) -> bool) -> Vec<BTreeSet<AtomId>>`,
dispatching to its rules: for Hückel and Clar the connected unions of rings
whose atoms can contribute, for HMO the π components that contain a ring atom.
An atom can still contribute when some choice that remains gives it an aromatic
contribution. The pruning cuts a subtree when, for some aromatic-only atom,
every candidate system that contains it fails accepts_contributions over its
members' remaining possible contributions. Remove claim_candidates and
accepts_range.

**Tests:** per rule, candidate_aromatic_systems and accepts_contributions;
existing selection results.

**Verification:** the resolve bench's selection group, against a baseline taken
before the change.

#### S4i — AtomSelectionAlgorithm and the search bound

**Modules:** graph resolve/aromaticity and resolve; umol-py resolve bindings,
the Python package exports, and Python tests. **Breaking, red→green within the
subitem.** [dep: S4h, S3a]

Add `pub enum AtomSelectionAlgorithm { FullDfs, BacktrackingDfs }`, defaulting
to BacktrackingDfs, and the field `atom_selection_algorithm` to
AromaticityResolveConfig. FullDfs visits every combination of a component;
BacktrackingDfs applies the cut under the Error `aromatic_valence_failure`
policy and otherwise visits what FullDfs visits. The enumeration bound counts
the perception calls of the search instead of the product of the flexible atoms'
options, and a search that exceeds it makes the result Underdetermined. Bind the
enum and the keyword in Python. Add AtomSelectionAlgorithm to the In code line
of the nomenclature guide's Algorithm entry. Add the algorithm as an axis of the
resolve bench's selection group.

**Tests:** the indole skeleton with a flexible ring nitrogen under HMO at
thresholds 0.2 and 0.1 resolves alike under both algorithms; a component whose
option product exceeds the bound but whose pruned search stays within it
resolves under BacktrackingDfs; a property over skeletons with undetermined
hydrogen counts under the Hückel, HMO, and Clar rules, for invariant 10: FullDfs
and BacktrackingDfs give the same resolution through
Resolver::resolve_with_report wherever both stay within the bound.

**Verification:** the selection group under both algorithms, and under
BacktrackingDfs against a baseline taken before the change.

#### S4j — One aromaticity resolution path

**Modules:** graph resolve/aromaticity and aromaticity. **Green.** [dep: S3l,
S4b, S4g, S4h, S3a]

AromaticityResolver::plan_resolve runs the selection on the molecule's own
atoms, with no completions, and plans the accepted systems with
AromaticityPerceiver::plan_aromatic_system. The selection applies the failure
policies to the derivation of each combination, so derive alone classifies the
inconsistencies, including the stored-system rule and the failure checks of S3l,
and the conformance validator keeps reporting its first inconsistency. derive's
stored-system rule reads the systems perceived from all atoms' contributions: a
stored system that no perceived system reproduces with equal members is an
AromaticSystemFailure, which replaces the separate perception of each stored
system from its own contributions. The selection applies
`aromatic_system_failure` to it; under Keep, the stored system's component gains
no system.

**Tests:** existing standalone resolution results; a stored benzene-ring system
inside naphthalene whose ten atoms assert one electron each is an
AromaticSystemFailure under Error and in conformance validation and gains no
system under Keep, and resolution never fails integrity; existing conformance
validation results; a property over skeletons with stored systems and aromatic
evidence on neighbouring rings, for invariant 2: a perceived system that
overlaps a stored one is an AromaticSystemFailure under Error, and under Keep
the stored systems are unchanged.

**Verification:** the resolve bench's AromaticityResolver::resolve group,
against a baseline taken before the change.

#### S4k — Aromatizer reads evidenced

**Modules:** graph transform/aromatizer and aromaticity. **Green.** [dep: S3d,
S3l, S4b, S4j, S3a]

The aromatizer runs derive with each atom's Kekulé contribution and reads each
atom's `evidenced` aromatic valence. An atom whose evidence is non-aromatic is
excluded from perception. A member of a perceived system whose evidence gives
another contribution is an AromaticValenceMismatch, identified by the id that
system receives when the plan applies. An atom that asserts aromaticity and that
no perceived system claims is an AromaticValenceFailure, and a bond that asserts
aromaticity outside every perceived system an AromaticBondConstraintFailure;
each is an AromatizeError::Contradiction. The rules do not consider a candidate
system that contains a bond asserting non-aromaticity.

**Tests:** a Kekulé benzene with one atom asserting non-aromaticity forms no
system; one atom asserting two electrons is an AromaticValenceMismatch; a Kekulé
naphthalene with a peripheral bond asserting non-aromaticity forms a system on
the other ring; a cyclohexane carbon asserting aromaticity is an
AromaticValenceFailure, and a cyclohexane bond asserting aromaticity an
AromaticBondConstraintFailure; the invariant 11 property of S3l, run on the
aromatizer.

**Verification:** the transform bench's aromatizer cases, against a baseline
taken before the change.

#### S4l — Crate-private aromaticity machinery and composite phases

**Modules:** graph aromaticity, transform/aromatizer, resolve,
resolve/aromaticity, resolve/valence, valence/atom_typing, valence/counts, and
valence/completion; tests/property/resolve.rs, tests/property/atom_typing.rs,
and tests/matching. **Breaking, red→green within the subitem.** [dep: S4i, S4j]

Make crate-private: AromaticityPerceiver and AromaticityDerivation; the hueckel,
hmo, and clar modules with their rule types, errors, and output;
electrons_from_kekule; select and select_with_report; the admit methods of
ValenceResolver, AtomTypingValence, and CountsValence; and ResolveState and
AtomCompletions. Remove the test-side exhaustive selection and select_scenario,
which S4i's property replaces; the other tests that called them use the public
operations.

**Tests:** the property suite and tests/matching through the public operations.

#### S4m — Planning rule in the data-type guide

**Modules:** docs/development/data-types.md. **Documentation.** [dep: S4d, S4e,
S4f]

Replace the transformer paragraph, which says transform_into uses a transaction
and restores on error, with the rule: a transformer derives a complete plan from
its input; every rejection is reported during planning, before mutation,
including a plan that would violate integrity; a plan applies to the molecule it
was derived from, so a failed application is a defect; transform consumes and
applies the plan to its input, and transform_into applies it in place without a
transaction or rollback method. State the same rule for the individual
resolvers' resolve_into and project_into. Check that the nomenclature guide's
Transformation entry still agrees.

#### S4n — Aromatic system names in the reworked code

**Modules:** graph aromaticity and its rule modules, resolve,
resolve/aromaticity, transform/aromatizer, transform/kekulizer, and
transform/delocalize_charge, with their tests and the call sites the renames
reach. **Breaking, red→green within the subitem.** [dep: S3m, S4d, S4e, S4f,
S4k, S4l]

In this code, identifiers that say `system` for an aromatic system say
`aromatic_system`, public ones included: find_systems becomes
find_aromatic_systems, the `systems` fields of AromaticityDerivation and
ResolveState become `aromatic_systems`, and the `system` fields of
AromaticityInconsistency and KekulizeError and the locals follow, as in
`stored_aromatic_systems`. AromaticityProjectError's NonConcreteSystem,
ChargedSystem, and SystemSpin become NonConcreteAromaticSystem,
ChargedAromaticSystem, and AromaticSystemSpin, with `aromatic_system` fields.
DelocalizationPlan's `system` field and the parameter of its derive become
`aromatic_system`, and the test cases assertion_without_system in the
aromaticity validator and in_system in counts become
assertion_without_aromatic_system and in_aromatic_system.

**Tests:** existing results under the new names.

### S5 — Copy-on-write storage

#### S5a — Molecule constraints

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

#### S5b — Reaction deltas

**Modules:** graph-IR reaction. **Green.** [dep: none]

Change the field to `Arc<Deltas>`. Reaction::try_new wraps the supplied Deltas.
normalize wraps its rebuilt Deltas. reframe_by, reframe, and into_parts use
Arc::unwrap_or_clone. The application iterator keeps its own Deltas field.

**Tests:** migrate the test literals and assignments that set the field.
Add clone-independence cases through normalize and reframe, and into_parts of a
shared and of a uniquely owned reaction.

#### S5c — ReactionSpan tables

**Modules:** graph-IR reaction_span. **Green.** [dep: none]

Change atoms, bonds, and constraints to `Arc<Vec<…>>`. Wrap at construction in
from_entries, try_remap, and reverse_sides. normalize, reframe_by, and reframe
rebuild each table with Arc::unwrap_or_clone. entries and try_remap clone the
inner Vec where they need an owned table.

**Tests:** clone-independence cases through normalize and reframe.

#### S5d — ReactionSpan overlay span newtypes

**Modules:** graph-IR dative, aromatic, multicenter, noncovalent, stereo.
**Green.** [dep: none]

Wrap the relation set of DativeBondSpans, AromaticSystemSpans,
MulticenterBondSpans, NoncovalentBondSpans, StereoAtomSpans, and
StereoBondSpans in an Arc. new, try_map, and remap wrap their result.
into_entries uses Arc::unwrap_or_clone. The Normalize, FrameTransport, and
reframe implementations write through Arc::make_mut.

**Tests:** migrate the unit tests that permute participants directly. Add a
clone-independence case per type through normalize and reframe.

#### S5e — Benchmark cases

**Modules:** graph-IR benches. **Additive, green.** [dep: S5a–S5d]

Add the aggregate_storage bench target with the twelve measured cases.


### S6 — Development guides

#### S6a — Ownership rules in the guides

**Modules:** docs/development/data-types.md and nomenclature.md.
**Documentation.** [dep: S0–S5]

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

#### S6b — Whitepaper code against the library

**Modules:** whitepaper/main.tex and umol-graph/tests/whitepaper.rs.
**Documentation.** [dep: S0–S5]

Check every code listing and inline code reference in the whitepaper against
the library after this plan: names, signatures, ownership, and results. Update
the listings that differ, and keep umol-graph/tests/whitepaper.rs compiling the
Rust-primer appendix listings exactly as printed.

**Tests:** the whitepaper test.

### S7 — Verification

#### S7a — Final gate

**Modules:** workspace. **Verification.** [dep: S0–S6]

Run workspace tests and lints with all features, and the Python suite in the
repository's Python 3.13 environment. Rerun benchmarks only for cases whose
production code differs from what was measured.

**Critical path:** S3a1 → S3b → S3d → S3l → S4j → S4l → S4n → S6a → S7a.