# 236 — Operation ownership and copy costs

Status: Completed
Date: 2026-10-09
Relates: [181](181-python-boundary-ownership-2026-08-03.md),
[168](168-api-hygiene-2026-07-27.md),
[213](213-editor-overlay-storage-2026-08-27.md),
[166](166-molecule-ops-2026-07-27.md),
[147](147-kekulization-single-structure-2026-07-13.md)

## Scope

This document owns the ownership contracts of graph-IR and graph operations:
which operations consume, borrow, or mutate their input, how they are named, and
the copy costs those contracts imply. It also owns the planning of resolvers and
transformers, with the shared aromaticity machinery, the constraint readings and
discharge, their ring context and Python interface, and copy-on-write storage
for the three graph-IR aggregates.

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

Placement is explicit. Composite resolution does not normalize the molecule-level
list or move its entity leaves into entity stores. Entity assertions remain
available as operands of logical constraints; accepting those expressions does
not require automatic placement. The MOL/SMILES ingestion path already writes
its assertions directly on entities and produces an empty molecule-level list.
Callers that want to move bare leaves can invoke inline_constraints themselves;
that operation leaves logical subtrees in place.

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

The aromatic perception machinery is shared, and no individual operation
duplicates it. The perceiver dispatches to its three rules, and its
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

The perceiver, its derivation, the hueckel, hmo, and clar modules with their
rule types, errors, and output, and ResolveState and AtomCompletions remain
public. AromaticityPerceiver's new, find_systems, and derive and the admit
methods of ValenceResolver, AtomTypingValence, and CountsValence also remain
public: they establish their ring context or candidate sets from the supplied
molecule. AtomCompletions is also part of the public ResolveReport and its
Python interface.

Restrict only the currently public methods that rely on unchecked relationships
between independently supplied inputs. The three rule types' find_from_rings
methods become crate-private because they assume the supplied RingSet belongs
to the molecule. AromaticityResolver's select and select_with_report become
crate-private because they assume candidate sets admitted for that molecule and
compatible pending aromatic systems in ResolveState. Public access to those
consumers would require additional input checks; the data carriers themselves
do not need restricted visibility. Their external test and benchmark callers
move to the public perceiver or resolver operations. Existing internal helpers,
including electrons_from_kekule, retain their visibility.

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

## Constraint readings and discharge

The S3b reading contract and measured implementation are accepted as of
2026-10-12. All six semantic readings on all eight constraints views return
normalized owned Option<F>, with meet's interpretation: None is bottom and
Some(undetermined(key)) is top. Storage absence contributes top, not a separate
reading outcome. The ten numeric projections listed below likewise return
Option<NumForm> with None as bottom. Interpretation is local to contributing
sources; storage getters still borrow, storage is not normalized eagerly, and
attribute mutation remains infallible. This settles invariant 13.

The optimized trial and retained caller optimizations justify accepting the
performance trade-off on the measured workloads. Small standalone owned-read
costs and rare-form copying remain as recorded; acceptance does not claim zero
overhead or universal speedups. The diagnostic contracts and Python reading
interface are also settled below. Shared container readings and diagnostic
formatting are implemented in S3b, owned/borrowed ring data in S3b1, and Python
ring bindings in S3b2. S3b3 implements Python constraint readings and completes
the numeric projections. Later perception and discharge work retains its
separately scheduled scope.

Every constraints view has six readings of a key: `asserted`, `derived`, and
`evidenced`, and their complete variants. Each reading is the meet of its
sources, starting from top. `asserted` meets the entity's own stored value,
`derived` the value derived from the relations the entity takes part in, and
`evidenced` all the evidence the molecule holds: the own assertion, the derived
value, and, for an atom's aromatic valence, each incident bond's assertion of
aromaticity, read as `#a+`. Each reading returns `Option<F>`: `None` is ⊥, as in
`meet`, and the key's vacuous form is no information, so every value has one
encoding. Top met with a value is that value's normal form, and a vacuous entry
and an absent one both give top, so every reading agrees exactly before and
after successful molecule normalization in the same entity and ring context.
On a molecule that normalizes, only the evidenced readings reach ⊥, from
conflicting evidence.

A complete variant replaces no information with the key's definite negative
where the key has one. For an assertion, no information is a top value: a
vacuous entry states nothing. For a relation, it is the absence of the overlay:
an overlay entity states that it exists even when its value is open, as a stereo
atom of undetermined kind, whose derived tetrahedral stereo is top, does not
close to `#T!`. `evidenced_complete` returns `evidenced` where it is not top and
`derived_complete` otherwise. Its claim is that the molecule states every
positive overlay fact, as an assertion or as a relation, without either source
alone being complete; `asserted_complete` claims that the assertions are
complete, and `derived_complete` that the relations are. Where a topology key
has a projection, `derived` and `derived_complete` agree. The own assertion
can be absent or disagree, and meeting it into evidence can yield ⊥. Every
reading that derives a ring key requires ring context.

Each constraint form type gains `undetermined(key)`, its vacuous form at a key,
which a reading returns for a key without a source. The constraint container of
each of the eight entity kinds has a named getter for every key, and its
constraints view has the same getter. Each returns `&T`, for `Copy` and
non-`Copy` values alike: the stored value, or a static vacuous value for an
absent key, so a vacuous entry and an absent key read alike. The getter of a key
with a parameter (a ring scope, a ligand permutation, or a ligand pair) takes the
parameter and returns the inner value; an iterator over the key's entries yields
references to the stored outer forms, so the entries can be listed without
knowing their parameters. Getters and iterators keep doc 122's rule that read
paths borrow and clone nothing. `get`, `contains`, and `iter` are unchanged.

`asserted_complete` reads only the entity's own storage and the absence cell. An
incident bond's assertion of aromaticity moves into the evidenced readings,
where it is met with the atom's own assertion instead of yielding to it: an atom
that asserts non-aromaticity next to a bond that asserts aromaticity reads ⊥,
not non-aromatic. The constraints views' `is_compatible` compares against
`evidenced`.

Valence admission and implicit hydrogen inference, under atom typing and under
counts, read `evidenced_complete`, and so do aromaticity perception, for each
atom's contribution, and stereo perception. Admission builds the ring set when a
registry row states a ring key. An atom whose evidence meets to ⊥ admits no
candidate, so atom typing and counts report their existing NoMatch. derive
reports an atom whose aromatic evidence meets to ⊥ as an AromaticValenceMismatch
when the conflict is with its system, and as an AromaticValenceBondMismatch when
the conflict is with an incident bond. Four outcomes change. Counts reads an
atom that asserts non-aromaticity inside a stored system as ⊥, not aromatic.
Perception reads an atom that asserts m electrons inside a stored system that
gives it n ≠ m as ⊥, not m. An atom without aromatic evidence outside every
system is excluded from perception instead of contributing zero electrons.
Implicit hydrogen inference, which export uses to omit stored H counts, infers
no count for an atom whose stored valence or aromatic valence conflicts with its
derived value, where the stored value took precedence, or whose only aromatic
evidence is an incident bond's `#a`, where the atom read as non-aromatic.
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

The composite resolver has no opening placement or closing discharge pass.
Discharge is part of each resolver's plan: the plan removes the assertions that
its own results realize.
Valence discharges `#v` on each atom it completes. Aromaticity discharges the
aromatic valence of each member of a system it adds, and the aromatic constraint
of each bond outside the stored systems that its result confirms; it writes no
aromatic constraint on bonds. Stereo discharges the tetrahedral or cis-trans
constraint at each site where it adds an entity. `reset_aromatic_valence` and
`reset_stereo_constraints` are removed. Assertions on entities that resolution
neither adds nor completes stay. Resolution no longer evaluates the other
assertions or the molecule-level list; the validators report a violated one.
It does not normalize or inline that list. Existing failure policies that
explicitly request removal still remove the selected assertions or entities;
entity removal retains its required constraint-reference compaction.

Standalone stereo resolution completes the structural step without requiring a
literal configuration. Positive open assertions `#T+` and `#C+` can become
entities with undetermined cosets, and are then discharged. Vacuous `#T*` and
`#C*` do not require entities. The composite resolver retains its final
concreteness check: an open entity leaves the result Underdetermined, and
in-place resolution preserves the input.

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

Explicit non-aromatic bonds delimit candidate aromatic systems under all three
rules, including HMO. Two aromatic regions joined by such a bond remain separate
systems, each evaluated under the selected rule. Giving the localized bond
assertion priority over HMO's grouping by connected pi components is a deliberate
pragmatic modeling choice. Within a fused region, the assertion excludes every
candidate containing both endpoints; admissible neighboring rings remain
candidates. Candidate construction and subsequent merging must both preserve
these boundaries, since a system contains every bond between its member atoms.
If merging or region recovery still produces a candidate containing both
endpoints of an excluded bond, perception rejects that candidate. It does not
search alternative partitions. This conservative loss of possible subregions is
accepted for the edge case; other independent candidates remain available.
HMO can numerically omit a coupling while retaining both endpoints, but the IR
does not store coupling exclusions within a system. Such a region is therefore
rejected before HMO evaluation rather than published with a contradictory bond
reading.

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

### S3b premise review (2026-10-11)

This review describes the pre-trial implementation and the questions that led
to the accepted signature and closure contract below. Its ownership and
normalization questions were settled on 2026-10-12.

The proposal conflates three decisions: which evidence contributes to a
reading; whether its result must be in exact normal form; and whether the API
returns borrowed or owned data. Semantic consistency and avoiding eager
rewrites both govern the decision. Neither preserving the current borrowed
lookup nor requiring a normalized result is justified by that principle alone.

Avoiding eager rewrites is a presumption, not an absolute prohibition. As with
the trade-off behind integrity checks, an earlier guarantee needs a concrete
account of the downstream simplifications it buys and the additional work it
imposes. Establishing normal storage at publication, normalizing a temporary
reading, and evaluating a semantic query without rewriting are different
choices: they affect different callers, preserve different raw observations,
and pay their costs at different frequencies. The review must identify which
consumers need each guarantee, what checks or special cases it removes, and
what copying or repeated computation it adds before settling the boundary.

- **Stored assertions already borrow.** All eight views' asserted accessors
  return Option<&F> directly from the constraint container. The typed getters
  and iterators also borrow and remain borrowed even in the written S3b plan.
  No source combination occurs in asserted. The claim that all assertion
  access must now return owned values is unsupported.
- **Exact normalization offers a downstream simplification.** A stored
  Valence(LitSet({2})) becomes Valence(Lit(2)) under normalization. A reference
  to the former cannot be a reference to the latter; a new value is needed if
  the accessor promises that exact result. The two stored readings already
  agree under normalized_eq, but that alone does not make them interchangeable
  for consumers: NumForm::as_lit and is_ground distinguish them by design.
  Returning normal forms is therefore a substantive API guarantee, not merely
  cosmetic equality. A borrowed raw form would require consumers to normalize
  it or restrict themselves to operations respecting normalized equivalence.
  Raw get/contains/iter can expose storage distinctions. Whether removing that
  burden from reading consumers justifies normalization at every reading, or
  at an earlier boundary, remains part of the trade-off.
- **Normalizing the result alone is insufficient.** A molecule with one bond
  whose order is LitSet({1}) gives derived Valence(Undetermined), whereas its
  normal form gives Valence(Lit(1)). The valence sum uses NumForm arithmetic,
  which returns Undetermined for non-literal operands. Meeting the result with
  top cannot recover the lost information. Also, an isolated atom's explicitly
  vacuous aromatic assertion currently gives asserted_complete top, but after
  normalization removes that entry it gives NotAromatic. These are semantic
  changes, not just different encodings. S3b's invariance requires examining
  derivation inputs and closure decisions, not only wrapping existing results
  in meet. Both discrepancies are reproduced by the public-API probe.
- **Some computed results require construction.** For example,
  {1, 2, 3} meet {2, 3, 4} is {2, 3}, which is neither operand. A derived count,
  an absent parameterized key's outer form, or a normal form different from
  storage cannot generally be returned as a plain reference to an existing
  form. This justifies construction on those paths, not cloning every stored
  assertion. Cow is not carried forward: it reintroduces the borrow-or-normalize
  API complexity removed with Normalize::normalized, without establishing a
  sufficient benefit.
- **Folding from top is not cost-neutral.** NumForm::meet clones and normalizes
  both operands, then clones the surviving operand in its top case. For a
  multi-element LitSet, even top meet x copies the boxed set twice. Expressions
  can copy trees; nested forms delegate to these operations. Existing
  matches/is_compatible implementations have borrowing shortcuts, while
  current view satisfies and valence admission already construct some meets.
  Account for copies in those actual callers; neither zero added cost nor a
  workload-wide slowdown follows from the signatures alone.
- **The consumers distinguish absence from failure.** Incidence validation's
  keyed asserted lookup currently treats None as no assertion to validate.
  A change to None meaning ⊥ must change that behavior deliberately, or it can
  report a contradiction as success. The new evidence reading can give None
  that meaning without requiring a storage lookup to use the same encoding.
  Stereo constraints currently have no derived projection at all, so their
  stored assertions are a particularly clear case without a merge to perform.
- **Two subsidiary claims need correction.** Only the two derived readings
  agree on derivable topology keys: an isolated atom without a degree
  assertion currently has asserted None and derived Degree(Lit(0)). Also, the
  earlier S3b sequence retained asserted_complete's incident-bond fallback
  until S3g while describing every reading as a meet. That fallback gives an
  own assertion precedence; it is not a meet with bond evidence. The revised
  migration below removes this intermediate contract.

The reproduced discrepancies identify work needed for the proposed invariance
law; they do not establish where normalization should occur or justify eager
rewriting by themselves. A plain reference to arbitrary stored F cannot
promise a normalized F, but that observation alone settles neither the result
contract nor ownership. Compare the guarantees, costs, and downstream
simplifications together before implementing S3b. At this review, only the
meet-compatible result encoding was settled; the later trial and measurements
support the accepted scope of owned readings and local normalization.

Evidence: [view/constraints](../umol-graph-ir/src/ir/view/constraints.rs),
[NumForm](../umol-graph-ir/src/ir/num.rs),
[atom typing](../umol-graph/src/ops/valence/atom_typing.rs),
[incidence validation](../umol-graph/src/ops/validate/constraint/incidence.rs),
and [stereo views](../umol-graph-ir/src/ir/view/stereo.rs).
The public-API probe in scratch/236-constraint-readings reproduces the stereo
matching discrepancy and the raw/normalized assertion distinction. The 51
existing constraint-container matching unit cases pass; they omit the failing
vacuous stereo cases.

### S3b producer and consumer analysis (2026-10-11)

This source analysis preceded the trial and acceptance. Its conclusions require
neither running tests nor benchmarking; they identify the benefit, dependency
boundary, and copying costs of normalized readings computed on demand.

**Consumer benefit.** Returning a normal form makes structural literal and
groundness inspection independent of the source's raw encoding. It does not
prove that the whole consuming operation is normalization-invariant: consumers
also inspect other raw attributes, such as element and charge.

| Consumer | Current observation | What the proposed readings supply |
| --- | --- | --- |
| Incidence validation; top-level stereo-constraint evaluation | is_undetermined and is_ground select success or underdetermination before matches. | A singleton set or reducible expression is classified by its normal form. A contradictory reading is distinct from an absent assertion. |
| Atom-typing admission and classification | Repeated per-row branches combine asserted and derived values, clone assertion-only values, and choose closure fallbacks. | One evidence/closure contract replaces the repeated source-selection branches. This simplification comes from evidence composition, independently of normalization. |
| Implicit-H inference in both valence implementations | Typed assertion getters take precedence over derived values; as_lit requires structural literals. | Normalized evidence supplies literal interpretation and detects disagreement instead of selecting one source. |
| Aromaticity derive | Raw aromatic assertions are tested with as_lit and matched against Aromatic(Lit); absence selects the stored-system contribution. | Normalized evidence removes encoding-dependent literal checks and centralizes the assertion/system/bond-source decision. |
| Stereo perception's skipped-bond selector | It checks only the outer Stereo variant, which normalization preserves. | This particular predicate gains nothing from normalizing or owning the coset. |
| Constraint-view satisfies and is_compatible | All eight views already construct combined host forms, with separate branches for missing sources. | Evidence readings centralize composition; these callers already pay some ownership costs. Borrowing matches/is_compatible shortcuts do not themselves require materialized normal forms. |
| Typed getters, get, and iter | Storage inspection and borrowing of assertions. | No normalized owned result is needed merely to inspect storage. |

The current resolver discharge also uses derived.is_ground after a meet and
can distinguish singleton sets from literals. Its removal is already planned;
its work is not a permanent downstream saving attributable to S3b. Python's
constraint views currently expose storage operations and typed getters, not
the keyed semantic readings. Changing Rust readings alone does not establish
Python parity.

Evidence: [incidence evaluation](../umol-graph/src/ops/validate/constraint/incidence.rs),
[stereo-constraint evaluation](../umol-graph/src/ops/validate/constraint.rs),
[atom typing](../umol-graph/src/ops/valence/atom_typing.rs),
[counts](../umol-graph/src/ops/valence/counts.rs),
[aromaticity](../umol-graph/src/ops/aromaticity.rs),
[constraint views](../umol-graph-ir/src/ir/view/constraints.rs), and
[Python constraint views](../umol-py/src/constraint/atom.rs).

**Producer dependencies.** Molecule normalization changes attributes and
constraint values, not topology, relation membership, or entity ids. Constraint
normalization preserves keys and removes vacuous entries. Therefore source
interpretation can be local to the requested key. The table covers the derived
keys of all eight constraint views; asserted reads need only their own entry.

| Derived key or group | Inputs requiring interpretation |
| --- | --- |
| Atom degree, ring degree, ring membership; bond ring membership | Topology and supplied ring context only; results are already literal. |
| Atom valence and ring valence | Orders of the contributing localized bonds; normalize before arithmetic loses non-literal information. |
| Atom donated and accepted pairs | Orders of the contributing dative bonds. Preserve the current limits on multi-donor projection; normalization does not define that missing chemistry. |
| Atom total degree | Topology counts plus the atom's implicit-hydrogen form. |
| Atom total hydrogens | The atom's implicit-hydrogen form and neighboring element forms; normalize before literal hydrogen classification. |
| Atom aromatic and multicenter valence | Overlay membership and the relevant positional electron-count cells. ElectronCountsForm normalization is identity; missing cells retain their current undetermined meaning. |
| Atom total valence | Localized valence, implicit hydrogens, and aromatic/multicenter contributions from the preceding rows. |
| Atom tetrahedral and bond cis-trans stereo | Overlay presence and configuration kind; for the matching kind, interpret the coset under that kind. Normalize the projected stereo form. A present unknown configuration remains open under closure. |
| Bond aromatic; binary dative aromatic | Overlay membership only; results are Boolean literals or top. |
| Dative ring membership; non-binary dative aromatic | No defined projection; reads top, including under complete derivation. |
| Aromatic-system and multicenter-bond electron count | Sum the stored electron vector by reference. No normalization or vector copy is needed. |
| Noncovalent intramolecularity | Localized-bond connectivity only. Its graph traversal is independent of attribute normalization. |
| Stereo-atom and stereo-bond constraint keys | No derived projection; all derived readings are top. |

Evidence: [atom projections](../umol-graph-ir/src/ir/view/atom.rs),
[ring projections](../umol-graph-ir/src/ir/view/ring.rs),
[bond projections](../umol-graph-ir/src/ir/view/bond.rs),
[dative projections](../umol-graph-ir/src/ir/view/dative.rs),
[aromatic projections](../umol-graph-ir/src/ir/view/aromatic.rs),
[multicenter projections](../umol-graph-ir/src/ir/view/multicenter.rs),
[noncovalent projections](../umol-graph-ir/src/ir/view/noncovalent.rs), and
[stereo projections](../umol-graph-ir/src/ir/view/stereo.rs).

For asserted_complete, normalize the selected assertion before deciding whether
it states anything. Close top only at a key with a definite negative. For
derived_complete, retain overlay presence separately from its open value;
normalization does not remove the overlay. The atom's additional aromatic
evidence needs only the incident bonds' Boolean aromatic assertions, not their
orders or whole bond forms. Evidenced readings combine these relevant sources;
their implementation need not construct each public reading separately and
then meet the resulting owned values.

The proposed failure boundary is the relevant source interpretation. An empty
numeric domain contributing to a valence sum yields bottom; an unrelated
contradictory attribute does not. An undetermined summand must not hide a later
contradictory contributing input. Electron-vector length and stereo frame
agreement remain separate contextual questions, not newly imposed normalization
checks. Attribute writes stay infallible because all interpretation is on the
read path and writes establish no new normal-form invariant.

The public projections need the same semantic contract. The
constraint producers share valence, total, dative-pair, and ring-valence helpers
with public entity-view accessors returning NumForm. Their current arithmetic
collapses non-literal inputs to top, including inputs that normalize to a
literal or contradiction. Normalizing that result cannot recover the lost
information. A semantic constraint reader must interpret inputs before that
collapse. The signature contract below follows that requirement;
it does not establish a normal-form invariant on storage.

**Copy classification.** Fix the input ownership, output type, and semantics
before judging a copy. With borrowed source storage, owned F, and the existing
Box/BTreeSet/String payloads, retained output data needs independent storage.
That establishes a lower bound on result construction, not a requirement to
clone the whole input. A copy is proven avoidable only when a concrete path
preserves the same contract without it. Normalization working space is a third
category: it may be justified by the chosen algorithm without being required
by output ownership.

| Case | Required by the result | Avoidable work demonstrated by the code |
| --- | --- | --- |
| Borrowed numeric set with several members, returned unchanged | Independent storage for the returned set while the original remains stored. | Top meet x currently clones the set into a normalized temporary, then clones it again; the temporary can be moved. |
| Borrowed singleton numeric set, returned as a literal | Copy the one integer; no set in the result. | Cloning the entire set before inspecting its cardinality; a borrowed cardinality/value check produces the same result. |
| Two numeric literal sets intersected | Construct the returned intersection if it remains a set. | meet first clones both operands, although their intersection can be computed from borrowed sets and only the result collected. |
| Already owned normal result | Transfer its ownership to the caller. | Cloning it merely to pass through a one-source or top-identity branch. |
| A symbolic expression normalizes to an expression | Independent representation of the retained result. | No general claim that all temporary work is avoidable; inspect the consuming normalizer before classifying each copy. |

Two additional implementation costs are visible without measurement.
stereo::canon_coset consumes its input but clones a plain LitSet while matching
by reference; it can move that set. The finite-domain stereo constraint
normalizer constructs sets even for Lit and Undetermined; those variants can
return directly. Thus even a scalar result can currently incur temporary
allocations, although output ownership does not require them. These are
analysis findings, not authorization to broaden S3b with normalization refactors.

**Optimization and measurement priorities.** Undetermined and Lit are by far
the most frequent variants in the code base, as supplied by the user; sets and
expressions are relevant but rare. Prioritize direct handling of Undetermined
and scalar Lit, with no unnecessary normalization or allocation. A straightforward
normalization fallback is acceptable for rare forms. The avoidable-copy examples
above identify possibilities, not a requirement to specialize every rare variant.

Performance comparisons must preserve the input distribution of the selected
workloads. Report targeted set/expression cases separately; do not combine
equally weighted variant cases into a headline result or inflate rare-form
frequencies to stress normalization in the primary comparison. No numerical
variant frequencies have been established. State the workloads and their
composition, and distinguish their results from general runtime claims.
Correctness obligations still cover every variant regardless of frequency.

This analysis supports local interpretation without a normalized-storage
invariant. Exact normal-form results have a concrete consumer benefit, but the
cost of retained owned payloads remains a trade-off. The following contract
specifies the signatures and closure boundary accepted after measurement.

### S3b signature and closure contract (accepted 2026-10-12)

**Settled result encoding.** Constraint readings share meet's Option semantics:
None denotes bottom, and Some contains a form, including top. This applies to
all six semantic readings and the numeric projections returning Option<NumForm>.
Unknown or absent evidence is Some of the appropriate undetermined form;
closure supplies the definite negative where its rule applies. This preserves
the existing open/closed distinction while removing None's former use for
missing evidence. For example, with no incident dative relation, derived at
DonatedPairs changes from None to Some(DonatedPairs(Undetermined));
derived_complete remains Some(DonatedPairs(Lit(0))). Ordinary storage lookups
still report entry presence; they are not semantic readings.

**Accepted output guarantee.** Each semantic reading returns a normal form:
Some contains the interpreted value, including top, and None means bottom.
Interpret only the contributing sources, before arithmetic or classification
can discard information. Use direct Undetermined/Lit paths and consuming
normalization as the fallback for rarer representations. This promises neither
a more precise symbolic arithmetic nor normalization of unrelated attributes.
Storage getters, constructors, attribute mutation, NumForm arithmetic, and
AsLit retain their contracts.

| Public surface | Return type | Relevant input |
| --- | --- | --- |
| All eight entity constraint containers: asserted, asserted_complete | Option<F>, with the entity's constraint form F and existing key argument | The own stored assertion; no molecule or ring context. |
| All eight constraints views: asserted, asserted_complete, derived, derived_complete, evidenced, evidenced_complete | Option<F>, with the entity's constraint form F and existing key argument | Sources selected by the reading and key. |
| AtomView: valence, donated_pairs, accepted_pairs | Option<NumForm> | Contributing bond or dative orders. |
| AtomView: total_degree | Option<NumForm> | Implicit hydrogens plus topology counts. |
| AtomView: total_hydrogens | Option<NumForm> | Implicit hydrogens and neighboring elements. |
| AtomView: total_valence, covalence | Option<NumForm> | Localized valence and implicit hydrogens, plus their existing overlay contributions. |
| AtomView: heavy_atom_degree, heavy_atom_valence | Option<NumForm> | Neighboring elements and, for valence, orders of contributing bonds. |
| RingAtomView: ring_valence | Option<NumForm> | Orders of incident bonds selected by the supplied ring context. |
| AtomView: degree, multicenter_degree, aromatic_valence, aromatic_covalence, multicenter_valence | NumForm, unchanged | Topology or electron-vector cells; no fallible normalization source. |
| RingAtomView and RingBondView: ring_degree where provided, ring_membership, ring_count, ring_size_count | NumForm, unchanged | Supplied ring topology. |
| Aromatic-system and multicenter-bond views, including editor and mutable views: electron_count | NumForm, unchanged | Borrowed electron-vector cells. |

The heavy-atom properties are included because their element filter also changes
under normalization: a singleton hydrogen set must filter like literal H.
Keep their existing policy that a genuinely unresolved neighbor counts as
heavy. A contradictory element gives bottom; a bond to a known hydrogen does
not contribute to heavy_atom_valence, so its order need not be interpreted.
Preserve total_hydrogens' existing arithmetic precision: with no neighbors its
implicit-hydrogen form survives, rather than being collapsed by adding zero.

Every reading and listed projection must agree exactly before and after a
successful molecule normalization, with the same entity and ring context.
Every Some result is normal. Relevant bottom propagates even if another input is
top; unrelated bottom does not poison a reading. No new electron-vector length,
stereo-frame, range, or mutation checks follow from this contract. The separate
contextual integrity contracts still apply.

**Complete readings.** Let A be normalized own assertion, D the open relation
projection, and E their evidence meet, including incident-bond aromatic evidence
on atoms. Bottom is preserved in every complete variant.

| Reading | Closure rule |
| --- | --- |
| asserted_complete | Return A unless it is top; then use the key's definite negative if one exists. Never read incident bonds. |
| derived_complete | Close absence of an overlay to its definite negative; preserve a present overlay's open value. For topology projections, return D. An unsupported projection stays top. |
| evidenced_complete | Return E unless it is top; then return derived_complete. |

The definite negatives are atom NotAromatic, NotMulticenter, NotStereo, and
zero donated/accepted pairs; bond false aromaticity and NotStereo; and dative
false aromaticity. Other keys have no closure negative. Nonbinary dative
aromaticity has no defined projection and remains top on the derived side.
An incident bond's true aromatic assertion contributes atom #a+ to E; false or
absent contributes no positive evidence. A present stereo entity of unknown
kind stays top under derived_complete; a known matching kind with an open coset
retains its positive stereo form. Missing ring context retains the existing
precondition and must not be encoded as bottom.

**Caller migration.** The Option change is semantic even where the Rust type
is unchanged. Existing source-selection matches use None for absence; retaining
them after the change would misinterpret bottom and omit closure for Some(top).

| Caller | Required migration |
| --- | --- |
| Constraint-view satisfies and is_compatible | Use evidenced_complete and evidenced respectively; bottom satisfies no entry and is incompatible. Include S3f in S3b. |
| Atom-typing and counts admission; implicit-H inference | Use the S3c evidence contract in the same change, including NoMatch/no inferred count for bottom. Do not flatten bottom into existing unknown-value defaults. |
| Incidence, ring, relational, and stereo constraint validation | Distinguish top from bottom before groundness tests. Use existing constraint-contradiction diagnostics or false predicate results for bottom. |
| Valence invariant checking and enumeration | Handle bottom before literal extraction; enumeration yields no candidates. The check uses the diagnostic contract below. |
| Counts classification | A contradictory projection cannot become zero through as_lit().unwrap_or(0); report CountsMismatch as specified below. |
| Current resolver discharge | Migrate None handling while this code exists; None can no longer mean no derived source to check. Its later removal is not a saving credited to S3b. |
| ECFP and Morgan extraction | Adapt property extraction under the existing molecule-is-concrete precondition, which excludes contradictory contributing attributes. Do not add a new public failure solely for the Option wrapper. |

**Settled diagnostics (2026-10-12).** These changes preserve the existing
operation-specific contradiction channels.

- CountsMismatch uses valence: Option<i64>, with Some only when a literal is
  available. Display the integer for Some and "valence unavailable" for None,
  without inventing zero. This field reports diagnostic availability, not the
  lattice reading. CountsError::NoMatch remains the admission outcome.
- ValenceInvariantsValidator already reports an electron-conservation failure
  as ValenceMismatch::OrbitalCount { atom_id, orbital_count, electron_count };
  the composite validator wraps it in ValidateContradiction::ValenceInvariant.
  Retain that existing diagnostic and terminology. It cannot describe a
  contradictory projection before those counts exist, so retain the trial's
  ValenceMismatch::ContradictoryConstraint { atom_id: AtomId,
  key: AtomConstraintKey }. Keep ordinary open inputs underdetermined. The
  standalone public validator cannot assume constraint validation already ran;
  mapping bottom to underdetermination would lose the new distinction.

**Settled standalone and Python interface (2026-10-12).** The two own-assertion
readings belong on each Rust entity constraint container. For example,
AtomConstraintsForm gains asserted(&self, key: AtomConstraintKey) and
asserted_complete(&self, key: AtomConstraintKey), both returning
Option<AtomConstraintForm>. The other seven containers use their corresponding
key and form types. These methods own the normalization and own-assertion closure
logic; molecule-backed constraints views delegate to them. Storage getters,
container lattice comparisons, and mutation retain their existing contracts.

Python keeps one constraints-view class per entity kind. Its existing backing
storage enum distinguishes a molecule entity from a standalone form; no separate
storage-proxy class is needed. Bind all six readings under their Rust names with
the corresponding key argument and F | None result. The dispatch contract is:

| View backing or context | Available operation | Outcome |
| --- | --- | --- |
| Molecule entity or standalone form | asserted, asserted_complete | Delegate to the Rust constraint container; None is bottom. asserted yields top for an absent assertion; asserted_complete closes top to the key's definite negative where defined. |
| Molecule entity | derived, derived_complete, evidenced, evidenced_complete | Delegate to the Rust semantic view with the same evidence and closure rules. |
| Standalone form | derived, derived_complete, evidenced, evidenced_complete, satisfies, is_compatible | Raise ValueError because molecule context is missing. |
| Molecule entity without required ring context | A reading or view comparison that needs ring derivation | Raise ValueError; missing context is not bottom. |

The view's satisfies reads evidenced_complete and is_compatible reads evidenced.
Their molecule-context requirement applies even when a query could finish
without examining a graph-dependent key. Raw container comparisons remain
separate lattice operations. asserted and asserted_complete need no ring context,
including at ring keys. A standalone form is not interpreted as an isolated
molecule atom with known-empty incidence.

Bind the ring-context creation and borrowing interface settled below, and the
atom/bond constraints views' with_rings interface. with_rings returns another
view of the same molecule entity with explicit ring context; it does not rewrite
storage or clone the molecule.
Semantic reading results are independent normal values. Ordinary attribute and
entity-constraint mutation stays visible through live views, while consumption
and structural mutation retain the existing accessor-invalidation contract.
Existing consumed-owner and invalidated-view exceptions remain unchanged.
The ten changed numeric projections bind as NumForm | None; typed storage
getters retain their raw-storage meaning.

S3b publishes the final asserted_complete rule together with its callers:
S3c, S3f, and S3g are included in that migration, retaining their detailed scope
below. No intermediate assertion reading keeps the incident-bond fallback.
The reading signatures, closure rules, diagnostics, Python reading interface,
and measured performance trade-off are accepted. S3b supplies the shared
container implementation, diagnostic formatting, and Rust verification;
S3b1 implements ring ownership, and S3b2 supplies the Python ring interface.
S3b3 binds the constraint readings and completes the numeric projections. The
comparisons below emphasize common variants and keep rare-form stress results
separate.

### Ring data and view ownership (accepted 2026-10-12)

RingSet owns the computed ring topology and retains no molecule reference.
Ring membership, sizes, and relations remain reusable while the covalent
topology and its atom/bond ids stay unchanged. Changing attributes, including
bond orders, does not require enumerating the rings again. RingViews borrows
the Molecule and internally holds Cow<'a, RingSet>; molecule-dependent readings
such as ring valence use the molecule's current attributes. The view exposes
immutable access, and Cow remains an internal ownership detail.

The ring accessor produces a separate computed data object. This distinguishes
it from atoms(), bonds(), and the other accessors that only expose views over
data held elsewhere. Computing rings at molecule construction or rebuild was
considered and rejected. Perception remains explicit, with the model and
algorithms selected at the call, and there is no molecule-level ring cache.

- Molecule::rings(&self, model, config) -> RingViews<'_> retains its current
  signature and computes the set once for that call.
- RingViews::new(molecule: &'a Molecule, rings: RingSet) -> Self remains
  crate-private. It consumes the set into Cow::Owned without cloning.
- RingViews::from_ring_set(rings: &'a RingSet, molecule: &'a Molecule) -> Self is
  public and stores Cow::Borrowed. Reconstruction neither allocates, clones,
  nor repeats perception.
- RingViews::into_ring_set(self) -> RingSet remains public and calls into_owned:
  it moves an owned set or clones a borrowed set into an independent value.
  The borrowed case's copying cost is explicit in its contract.

A caller can extract the owned set, mutate molecule attributes, and reconstruct
a borrowed view over the same set. Contextual use requires that the set
describes the supplied molecule's topology and id space. Mismatched inputs must
not panic but have no specified semantic result. Local guards protect contextual
access; reconstruction does not revalidate the entire topology or repeat
perception.

RingSet keeps its existing topology queries and remains operation-produced;
this change adds no arbitrary-parts constructor. The atom/bond constraints
views' with_rings methods continue to borrow RingSet. Existing owning RingViews
consumers and rings(model, config).into_ring_set() callers retain their routes.
Consumers holding an independent RingSet can reconstruct temporary borrowed
views when needed.
Ring models, enumeration algorithms, query meanings, and the current single-ring
member accessors retain their contracts.

The Python surface covers RingSetKind, RingModel, RingRelation, RingConnection,
RingGraph, RingSet, RingViews, RingView, RingAtomView, and RingBondView, alongside
the existing RingConfig and RingScope. Bind the corresponding Rust operations
and the established collection protocols. A freshly computed Python RingViews
owns the computed set and retains the Python Molecule owner; a reconstructed
view retains the Python RingSet and Molecule owners. Each call reconstructs
temporary Rust borrows without copying the ring data. Bind into_ring_set with
the same consuming semantics and move-versus-clone costs; consumption invalidates
the view's dependent accessors. Nested access and iteration retain their required
owners and remain lazy. Attribute mutations remain visible through
molecule-backed views; molecule consumption and structural operations invalidate
those views under the existing accessor-counter contract. An independently held
RingSet remains available for its own topology queries. RingSet::graph retains
its current Rust contract of returning an owned graph copy, with no additional
binding-side graph clone.

The scope is this ring-ownership change and these bindings. The broader
view-argument audit and ring-member-accessor redesign are outside this work.

### S3b Rust trial (2026-10-12)

The authorized in-tree trial implements the six normalized Option<F> readings
on all eight views and the ten Option<NumForm> projections. Raw getters still
borrow stored forms; mutation remains infallible and stores the supplied forms.
Normalization is local to reads. Undetermined and Lit numeric/element sources
avoid normalization and heap cloning; evidence composition moves a top identity
operand's counterpart instead of cloning it. Finite stereo relation forms also
return their already-normal Undetermined/Lit values directly.

The trial includes the coupled S3c/S3f/S3g caller migration, existing discharge,
constraint validation, and fingerprint extraction. Its CountsMismatch and
ValenceMismatch payload changes are now accepted. The Python workspace still
compiles. At this trial snapshot the shared container readings and accepted
semantic-reader binding surface were not yet implemented; S3b completion below
supplies the Rust containers and diagnostic formatting. S3b1–S3b3 supply ring
ownership and bindings. S3d/S3e perception work and the separately verified
S3a3/S3a4 items are not included.

**Correctness.** Fixed cases distinguish absence, explicit top, literal,
singleton, contradiction, top followed by bottom, and bottom followed by top.
They check local interpretation without storage rewriting, infallible mutation,
heavy/H classification, closure, and incident-bond aromatic evidence. Properties
compare all six readings across successful whole-molecule normalization and
check numeric evidence composition against an independent meet calculation.
The reader property visits all eight entity families when generated, rather
than assuming every generated molecule contains every family.

Three matching expectations change with the settled evidence closure: a marked
aromatic bond supplies positive atom evidence, excludes a negative atom query,
and is not contradicted merely by the absence of an aromatic overlay. Resolver
regression cases preserve stereo-entity assertions now that unsupported
projection returns Some(top), rather than None. Bottom aromatic evidence on a
ground atom produces NoMatch in both admission policies. A bottom localized
valence exercises the new invariant diagnostic.

The two generated H-elision suites retain their independent allowed-H oracles,
registry permutations, duplicate rows, and generated valences. Their fixtures
now realize a nonzero asserted valence with a bond; the supporting isotope-
specified atom retains H explicitly. The former isolated-atom fixture
contradicted its own nonzero valence assertion under the meet contract. Separate
inference cases retain contradictory assertion/topology inputs and expect no
inferred count.

Verification passed:

- Graph-IR tests with proptest: 9,791 unit tests and 430 property-target tests,
  plus all enabled integration targets; four existing ignored tests remain.
- Graph tests with proptest and conformance: 1,994 unit tests, 171 property-target
  tests, 683 resolution conformance cases, and the other integration targets.
- Workspace/all-targets check using Python 3.13.15; Clippy for both changed
  crates, all targets, with the property and conformance features and warnings
  denied; formatting and diff checks.

**Measurement method.** The committed baseline is
24bb200dd7eab8df2c1fdcbc061f338812c2086b, built with rustc 1.96.0
(ac68faa20 2026-05-25). Saved baseline executables were measured before changes
and again after the initial candidate. The final candidate includes the
correctness fixes. Runs are sequential, with no concurrent builds or tests.
Criterion uses 20 samples, 0.2 s warmup, and 0.5 s requested measurement time
(long-running matching cases use Criterion's longer estimate). Reported times
use Criterion's slope estimate when available and mean otherwise. The two
baseline runs show drift separately from candidate differences; these short
runs do not establish a production workload distribution.

The new constraint_readings bench measures atom asserted, derived,
derived_complete, valence, and satisfies calls. Absent, Undetermined, and Lit
cases are reported separately from singleton sets, larger sets, and expressions;
there is no equal-weight aggregate. Allocation counting is a separate run of
1,000 calls per case, outside timed runs. Bytes count allocation/reallocation
requests, not peak live memory. The harness includes common dispatch and
black_box overhead. Baseline asserted borrows; candidate asserted owns and
normalizes, so these are costs of the changed contract, not equivalent raw
storage accessors.

Saved executables, commands, per-case estimates and confidence intervals,
allocation counts, and full comparison CSVs are in scratch/236-read-costs.
The tracked bench is umol-graph-ir/benches/constraint_readings.rs. Reproduce the
saved comparison with measure.py baseline-repeat baseline and measure.py
candidate-final, sequentially, then summarize.py candidate-final.

**Measured results.** Every common microcase (absent, Undetermined, Lit,
across all five measured calls) makes zero allocations in both versions.
Returning an owned scalar form does not introduce a heap clone. It does have a
CPU cost in asserted/evidence paths. The two candidate runs agree on that
micro-level direction and scale. Selected final-candidate times in ns/call:

| Case | Baseline | Baseline repeat | Final trial |
| --- | ---: | ---: | ---: |
| Absent asserted | 4.36 | 4.18 | 6.09 |
| Undetermined asserted | 5.74 | 5.94 | 15.44 |
| Lit asserted | 5.72 | 5.71 | 16.19 |
| Lit derived | 12.16 | 12.22 | 10.67 |
| Lit valence | 10.49 | 10.44 | 9.13 |
| Absent satisfies | 19.52 | 19.23 | 26.82 |
| Undetermined satisfies | 32.59 | 31.48 | 41.18 |
| Lit satisfies | 31.47 | 31.16 | 43.98 |

Common derived/derived_complete/valence cases are approximately 11–15% faster
than the repeated baseline. Common satisfies cases are 31–41% slower, or
approximately 8–13 ns more per call. The much larger asserted percentage
compares a cheap borrowed lookup with an owned semantic read; the unchanged
raw getters remain available for callers that need storage access.

Whole-operation examples, in microseconds per operation:

| Case | Baseline | Baseline repeat | Final trial | Increase against the two baselines |
| --- | ---: | ---: | ---: | ---: |
| Ingest benzene | 19.03 | 18.70 | 19.92 | 4.7–6.5% |
| Ingest octane | 12.76 | 12.62 | 13.42 | 5.2–6.4% |
| Resolve branched | 11.98 | 11.79 | 12.57 | 4.9–6.6% |
| Resolve chain_64 | 112.76 | 110.42 | 119.00 | 5.5–7.8% |
| Resolve aromatic_fused | 43.84 | 42.16 | 46.25 | 5.5–9.7% |
| Resolve shared_chain | 12.99 | 12.86 | 14.38 | 10.7–11.9% |

The seven ingest cases and sixteen resolution-only cases are retained
individually in the CSV. Their changes are not combined into a workload average.
The resolution-only bench checks that each case determines before timing it.
Alkene_four has 9.4% baseline drift and opposite signs against the two baselines;
tetra_ring is effectively unchanged against the repeat. Those cases do not
support a precise regression estimate. Other cases and the initial candidate
also vary, so the table is evidence of this trial's cost, not a universal rate.

The six existing OpenSMILES/VF2 corpus matching runs were also retained. Their
patterns carry no constraint entries, and corpus loading silently skips ingest
failures without recording workload parity. They therefore do not establish
the cost of semantic constraint reads or equivalent end-to-end matching work.
The direct satisfies benchmark does exercise a literal valence constraint.

Rare-form allocation results are separate. Each cell is baseline → trial,
written as allocations/requested bytes per call:

| Source form | asserted | valence | satisfies |
| --- | --- | --- | --- |
| Singleton set | 0/0 → 2/128 | 2/128 → 2/128 | 4/256 → 4/256 |
| Two-member set | 0/0 → 3/152 | 2/128 → 3/152 | 7/408 → 6/304 |
| Constant expression | 0/0 → 3/224 | 2/96 → 3/224 | 5/320 → 6/448 |

Final asserted timings for these rare cases are 52, 68, and 91 ns, respectively;
valence is 51, 65, and 87 ns. The baseline borrowed assertions take about 6 ns.
The rare derived results can differ semantically from the baseline because
normalization now precedes information-losing arithmetic. These are not simply
faster/slower implementations of an unchanged raw read.

**Initial-trial assessment.** The semantic scheme is implementable without eager
storage normalization or fallible attribute mutation. The tested common paths
avoid new allocations, but the added CPU cost is measurable. These results do
not establish an unavoidable cost of ownership or Option<F>: the trial also
changes evidence composition, closure, and callers, and retains some repeated
reads. Rare-form materialization has not been optimized. The optimization
measurements below support acceptance of the revised implementation's cost;
the trial does not justify declaring all assertion accessors owned.

### S3b cost attribution (2026-10-12)

The trial raised two questions: how much additional work the revised semantics
require, and how much of the measured slowdown comes from avoidable work in
the implementation. Allocation-free readings alone answer neither question;
the following audit and caller experiments establish the measured acceptance
trade-off without claiming an exact lower bound on semantic cost.

The expected baseline for centralization is roughly unchanged total work when
the core methods replace equivalent work formerly performed by callers. A
check moved into a reading must not also remain in its caller without a reason.
Centralization can improve cost by sharing results or simplifying branches; it
does not inherently require more meets or more traversals. Function boundaries,
enum materialization, dispatch, and inlining can still affect CPU time even
without allocation. Those are implementation effects, not new semantic
obligations or demonstrated unavoidable costs of Option<F>.

Some obligations are genuinely stronger. A non-normal source must be
interpreted before arithmetic loses its information, relevant bottom must
survive other top inputs, and evidence from independently supplied sources must
agree. For example, H inference previously accepted a literal assertion in
preference to derivation; its new reading must also examine relation evidence
to detect a conflict. However, Undetermined and Lit are already normal: their
contract does not require a general normalization traversal. Nor does a check
newly added to a reader establish net new work for a whole operation if that
operation already performed the check elsewhere. Reuse must respect intervening
mutation and each consumer's required result and diagnostic behavior.

Source comparison against the measured committed baseline identifies concrete
repetition in the trial:

| Path | Finding and baseline distinction |
| --- | --- |
| All views' evidenced_complete | Calls evidenced, which derives the key, then calls derived_complete when the result is top. For topology keys those derivations are identical, so a top valence can traverse the same bonds twice. Overlay closure can also repeat presence queries; eliminating that repetition must preserve the distinction between an absent overlay and a present overlay with an open value. |
| CountsValence::admitted_completions → CountsInput::for_molecule_atom | The new upfront aromatic evidenced_complete result is discarded on the ordinary non-ground path and recomputed by CountsInput. The valence precheck followed by another valence derivation also repeats work, but that duplication already existed in the baseline. It is an optimization opportunity, not solely a new cost attributable to centralization. |
| AtomTypingValence admission | The trial reads aromatic evidence before row selection and reads it again for rows carrying the aromatic key. The old preliminary assertion reading was conditional on the atom being ground. Per-row recomputation of other keys existed before the trial; those readings are unchanged between rows for a fixed molecule and context, though their comparison with each row still differs. |
| asserted → meet_readings → NumForm::meet | asserted clones and normalizes the own assertion. When both evidence operands are non-top, NumForm::meet clones and normalizes both again, including the already-normal derived value. The baseline meet already normalized its operands; the trial adds normalization before that meet. The top-identity path in meet_readings avoids this second meet, so it is not a cost of every read. |

In particular, the baseline satisfies implementation already met asserted and
derived evidence and then compared the host reading with the pattern. The
literal satisfies slowdown cannot be explained simply by saying that meeting
was newly required. The trial adds preparation of an owned normal assertion,
which the subsequent numeric meet prepares again. In the top numeric case,
the complete reading can additionally derive the same value twice. These are
source-level findings; compiler optimization and measured contribution still
need to be distinguished from the apparent work in the source.

Meeting the host's evidence and comparing that host with a query are different
obligations, so they are not automatically duplicate checks. The common numeric
query comparison already uses a direct matches/is_compatible path rather than
another general meet. Removing repeated preparation or repeated derivation
must preserve both obligations, bottom propagation, and the closure policy.

The trial therefore contains avoidable work; the recorded timings do not say
how much of the slowdown it explains. The discriminating comparison is between
implementations with the same new semantics, removing identified repetition
and specializing already-normal common forms, then measuring the same common
paths and whole operations against the saved results. Rare forms remain
separate. The original-versus-trial comparison alone cannot isolate the cost
of stronger semantics from the cost of this implementation. This attribution
review records findings only; it does not change the trial code or claim that
all of the slowdown can be recovered.

### S3b implementation cost audit (2026-10-12)

The optimized trial keeps the same semantic contract. All six semantic readers
return normal Option<F> with None as bottom; storage remains raw, its getters
borrow, and mutation remains infallible. No public symbols, unchecked reading
methods, or caller-supplied trust flags were added in this optimization pass.

**Work removed.** Complete evidence now derives each key once. Own and
incident-bond evidence determine whether that derivation uses closure: a
non-top meet operand cannot produce a top result, so closure cannot subsequently
be needed in that branch. With top evidence, the complete derivation already
preserves present open overlays and closes absent ones. The property suite
checks this against the definition: complete evidence equals open evidence
unless it is top, in which case it equals the complete derivation.

Numeric assertions use the existing numeric reader directly, avoiding the
enclosing constraint's clone/normalize path for Undetermined and Lit. Absent or
top non-aromatic assertions go directly to derivation without materializing an
intermediate top constraint. An already-positive own aromatic reading needs no
incident-bond scan: positive bond marks can add no information to it, and
negative marks do not negate atom aromaticity. The own reading is normalized
before that shortcut; contradictory own evidence still returns bottom.

Composition moves top and identical normal readings directly. NumForm::meet
handles scalar top/literal pairs directly, and aromatic/multicenter meets let
their numeric meet normalize the payload once instead of preparing it again
at the enclosing level. Numeric and atom normalization have inline hints.
Counts admission reuses its upfront aromatic reading; atom typing reuses it
for aromatic row constraints and computes groundness once. All reuse is local
to an immutable molecule/context.

**Comparison.** The committed baseline remains 24bb200dd7eab8df2c1fdcbc061f338812c2086b.
Saved binaries isolate the original trial, an intermediate optimization, and
the final optimization. The first series ran baseline → trial → intermediate
→ intermediate → trial → baseline. The final series ran final → baseline →
final; its repeated reading timings cover the common cases, with rare cases
timed separately in the first final run. Each case uses 30 samples, 0.3 s
warm-up, and 1 s measurement. No builds or correctness checks ran alongside
timing. The tables give arithmetic means of two per-case point estimates,
using the last two baseline runs; they do not average different workloads or
assign a variant frequency.

| Common numeric reading, ns | Baseline | Original trial | Optimized trial |
| --- | ---: | ---: | ---: |
| Absent assertion: asserted | 4.25 | 5.97 | 6.02 |
| Literal assertion: asserted | 5.80 | 15.14 | 10.06 |
| Top assertion: asserted | 5.86 | 15.10 | 9.37 |
| Absent assertion: satisfies | 19.61 | 26.28 | 18.30 |
| Literal assertion: satisfies | 31.51 | 43.12 | 30.51 |
| Top assertion: satisfies | 31.87 | 40.34 | 20.07 |

The final common numeric valence and derived readings also remain faster than
the committed baseline in these point estimates. Standalone owned assertion
reads retain an absolute cost of roughly 2–4 ns over the old borrowed result;
the composed matching path recovers that cost through less work elsewhere.
All 15 tested common microcases still make zero allocation/reallocation
requests. This is not a whole-workflow or all-form allocation claim.

| Whole operation, µs | Baseline | Original trial | Optimized trial |
| --- | ---: | ---: | ---: |
| Ingest octane | 12.74 | 13.58 | 12.04 |
| Ingest benzene | 19.13 | 19.87 | 18.58 |
| Resolve chain_64 | 111.76 | 113.57 | 108.48 |
| Resolve aromatic_fused | 43.07 | 43.46 | 43.39 |
| Resolve shared_chain | 13.09 | 13.22 | 13.11 |
| Resolve shared_branch | 15.08 | 15.31 | 15.46 |

Across all 23 ingest/resolution cases, the optimized per-case changes range
from −5.5% to +2.5% against baseline. Eleven are within ±2%, eleven are more
than 2% faster, and one is more than 2% slower. These are descriptive counts,
not significance tests or workload weights. Baseline repeats drift by up to
4.7%; the two final shared_branch runs differ by 6.8%, exceeding its +2.5%
midpoint difference. Small differences do not establish a remaining semantic
penalty. The earlier trial's precise percentage regressions are likewise not
stable across timing rounds.

Rare forms remain separate. Final asserted readings take about 50 ns for a
singleton set, 64 ns for a two-member set, and 90 ns for the constant expression;
they still request 2/128, 3/152, and 3/224 allocations/bytes respectively.
Their satisfies allocation counts remain 4/256, 6/304, and 6/448. These costs
are not demonstrated lower bounds: for example, a singleton set could be
recognized without first cloning its allocation.

**Attribution and limits.** The original slowdown was substantially avoidable.
The optimized implementation reaches roughly the committed baseline's cost
on the measured whole operations while preserving the new semantics. This
does not prove zero marginal semantic cost: primitive meet improvements also
benefit code that existed before the trial, and this is not a comparison with
an independently optimized old implementation. Nor do these measurements
cover every constraint kind or atom-typing row workload. The corpus matching
controls excluded above still supply no evidence for constraint-reading cost.

Two further repetitions were identified for the follow-up below. Counts
admission's raw-derived
valence precheck is distinct from testing whether combined evidence is literal:
an own literal assertion can make the latter literal while derivation remains
unknown. Removing the precheck changes admission; avoiding its second traversal
requires reusing the derivation while retaining that policy. Atom typing can
still reread non-aromatic keys across rows. Neither repetition is an inherent
semantic cost or a reason to introduce unchecked public APIs. No cache spanning
mutation or broad replacement of semantic reads with raw getters was added.

For the measured common paths, this evidence and the caller follow-up below
support the accepted semantic contract and optimized implementation. They supply
no performance justification for eager storage normalization, weaker bottom
propagation, or a new unchecked reading surface. Performance acceptance and the
diagnostic and binding contracts are settled. The S3b completion below records
the final Rust work; ring ownership and bindings follow in S3b1–S3b3.

**Verification and artifacts.** Graph-IR tests with proptest pass (9,794 unit
tests and 430 property-target tests, plus integrations; four ignored). Graph
tests with proptest and conformance pass (1,994 unit tests, 171 property-target
tests, 683 resolution conformance cases, plus integrations). Relevant-crate
Clippy with those features and all targets, workspace all-target checks with
Python 3.13.15, formatting, and diff checks pass. Added aromatic cases cover
positive literal, singleton, and contradictory own evidence; the closure law
is checked through public readings. The existing mutation and bottom tests
remain intact.

The binaries, hashes, source snapshots, allocation counts, logs, per-run
confidence intervals, and runner are under scratch/236-read-costs:
optimization-experiment.json, optimized-final-source.patch,
optimized-final-added-files, final-comparison.csv, the audit-* result files,
and measure-focused.py. These preserve intermediate measurements as well as
the final comparison.

### S3b related caller optimizations (2026-10-12)

This follow-up tests the remaining caller repetitions as performance work,
independently of whether they arose with the semantic change. The selected
implementation preserves the reading contract, raw storage, infallible mutation,
and existing admission policy. It adds no public API or unchecked path.

Atom typing reuses the registry lookup slice and caches requested complete
readings within one atom's admission scan. Cached bottom rejects only rows
requesting that key; it does not reject later rows that omit it. The cache
includes the ring scope in its key, uses the same immutable ring context, and
stores the five default-registry keys inline. Registries requesting more keys
can spill to the heap. No cache survives the call or spans mutation.

Counts reuse was tested in two forms. Both retained the admission precondition
that derived valence be literal, including the skip when an own assertion could
make unknown derivation literal. Meeting an assertion with an already-literal
derivation can only return that same literal or bottom. The first version
constructed a normal owned assertion and checked compatibility; its admission
cost was essentially unchanged. The second used NumForm::is_compatible directly,
whose contract is equivalence to canonicalizing meet existence, avoiding both
the repeated traversal and assertion materialization. Direct admission improved
about 1–3%, but shared-chain resolution was 7.9% slower than the control midpoint,
with both repeats slower. Partial-triene resolution also spiked in one run.
These observations do not establish a cause or an aggregate workload cost, but
the benefit was insufficient to retain the rewrite. Counts production code was
restored to the previous optimized trial; the added semantic regression cases
remain. The derivation reuse is not an inherent semantic impossibility.

The literal-derived guard skips a non-ground atom with unknown incident bond
orders; it does not report NoMatch for that unknown value. Direct counts
admission can therefore finish with no completion for that atom. A literal
own assertion does not override the guard; the surrounding valence resolver
can also return Underdetermined earlier when it cannot verify that assertion
against incidence. Bond defaults do not infer orders. An otherwise ordinary
acyclic input with unknown orders remains non-concrete and overall resolution
returns Underdetermined. Bottom at the derived-valence check instead produces
NoMatch. This admission restriction predates the reading changes.

Two atom-typing experiments were rejected. Extending this cache to ground-atom
classification increased its measured cost by 4.9–11.5%; attribute filtering
already avoids most repeated reads there. Classification keeps direct readings
and the single registry lookup. Unconditionally testing row constraints before
attributes accelerated open-H admission but slowed methyl-radical admission
and resolution. The selected ordering tests attributes first when the atom's
implicit-H form is constrained, and cached constraints first when it is top.
Each conjunction is unchanged, and attribute compatibility is evaluated at most
once per row. This is an internal performance choice, not a new input contract.

The new valence_readings benchmarks use the default registry and ordinary
raised/resolved SMILES inputs: octane, a 64-carbon chain, benzene, pyridine,
ammonium, and methyl radical. They measure counts admission, typing admission,
typing classification, and complete typing resolution with MostSaturated.
Parsing, raising, registry initialization, and resolver setup occur outside
measurement; mutable resolution also excludes input cloning. Setup checks
successful admission, classification, and resolution. These workloads exercise
literal/undetermined forms; they imply no population frequencies and contain
no synthetic set/expression weighting.

**Final measurements.** Each case uses 30 samples, 0.3 s warm-up, and 1 s
measurement, without concurrent builds or tests. The first experiment used
control → candidate → candidate → control. Subsequent runs separately preserved
the rejected filter ordering, the counts compatibility trial, two fresh runs
of the committed baseline, a third control, and two final runs. The manifest
records the complete order. The table averages each case's last two control
point estimates and two final point estimates; it does not pool workloads.
“Control” is the optimized S3b trial before this follow-up, with the same new
benchmark functions. All times are microseconds.

| Typing workload | Admit control | Admit final | Resolve control | Resolve final |
| --- | ---: | ---: | ---: | ---: |
| Octane | 12.933 | 10.338 | 23.469 | 20.479 |
| 64-carbon chain | 105.954 | 83.228 | 192.994 | 175.592 |
| Benzene | 6.699 | 5.031 | 24.302 | 23.310 |
| Pyridine | 6.952 | 5.389 | 25.208 | 23.962 |
| Ammonium | 0.611 | 0.619 | 1.689 | 1.720 |
| Methyl radical | 0.819 | 0.815 | 1.911 | 1.901 |

Admission improves 20–25% and complete typing resolution 4–13% on the four
chain/aromatic inputs. Ground typing classification improves 2–5% after retaining
only lookup reuse there. Methyl radical stays essentially unchanged. Ammonium
admission is 1.4% slower (about 8 ns), and resolution is 1.8% slower (about 31 ns);
the latter is within its final repeat spread. The result is not a uniform
speedup across every input.

With counts production code restored, direct counts admission differs by
−0.6% to +1.7% from control. The 23 existing ingest/counts-resolution cases
range from −3.4% to +2.5%; the shared-chain case returns to control cost
(12.526 → 12.505 µs). Relative to the two fresh committed-baseline runs, those
same cases range from −7.8% to −0.5%. Baseline repeats themselves differ by up
to 5.5%, and control/final repeats also vary. These small percentages are
measurements, not proof of exact semantic overhead or universal gains. The
large typing-admission reductions survive the repeated controls clearly.

The raw storage and six semantic readers are unchanged by this follow-up.
Default-registry reading reuse adds no heap allocation for the cache itself;
this is not a claim that admission as a whole is allocation-free. The final
tree retains the typing changes and restores the counts implementation.
The subsequent S3b completion is recorded below.

**Verification and artifacts.** The final graph suite passes with proptest and
conformance: 2,010 unit tests, 171 property-target tests, and 683 resolution
conformance cases, plus integrations. Graph all-target Clippy with those
features and warnings denied passes. Sixteen new exact cases cover counts'
unknown-derived skip, assertion conflicts and bottom, repeated registry keys,
top/bottom readings, a later row omitting a bottom key, and distinct ring scopes
with more than five cached keys. Formatting and diff checks pass.

The experiment manifest, immutable binaries, release-source snapshots, full
source patch and added files, per-run confidence intervals, comparisons, and
logs are under scratch/236-read-costs. related-experiment.json maps every
variant; related-accepted-resolve is the final binary. The final comparisons
are related-accepted-summary.json and the CSV files comparing related-control-b,
related-control-c, related-accepted-a, and related-accepted-b, plus the separate
committed-baseline comparison. Earlier binaries remain as rejected experiments.

### S3b Rust completion (2026-10-12)

All eight Rust constraint containers now supply asserted and asserted_complete;
the views delegate their own-assertion readings to them. The shared atom reader
retains the common-variant fast paths and evidence's single storage lookup.
Raw storage and infallible attribute mutation are unchanged. CountsMismatch
prints a literal valence as its integer and an unavailable valence as
“valence unavailable”; public classification preserves None for a bottom
projection. The accepted ValenceMismatch variants and coupled S3c/S3f/S3g
migrations are retained. S3b is complete; S3b1–S3b3 follow separately.

**Correctness and verification.** Added 102 container cases cover absence, top,
normalization, bottom, closure, scoped keys, and unchanged storage. The existing
normalization property now also checks container/view equality for both own
readings across all eight families. Six diagnostic cases cover display and
public counts classification without substituting zero for unavailable valence.
The all-feature workspace test run passes 36,927 cases, including 430 graph-IR
properties, 171 graph properties, and 683 resolution conformance cases. Three
subsequently added classification cases pass in the focused ten-case run.
Workspace all-target, all-feature Clippy with warnings denied, graph-IR rustdoc
with warnings denied, formatting, and diff checks pass. Workspace commands use
Python 3.13.15. The source and API review found no changes outside S3b's agreed
Rust surface; the new Python methods are still scheduled in S3b2–S3b3.

**Performance.** Repeated comparisons use the saved accepted reader/caller
binaries for the completion changes, and the committed pre-S3b binary for the
twelve fixed matching fixtures. Per-case timing changes are:

| Cases | Count | Control | Change range |
| --- | ---: | --- | ---: |
| Common readings | 15 | Accepted reader implementation | −8.2% to +0.5% |
| Rare set/expression readings | 15 | Accepted reader implementation | −1.6% to +1.4% |
| Ingest and resolution | 23 | Accepted caller implementation | −1.0% to +1.9% |
| Direct valence workloads | 24 | Accepted caller implementation | −1.4% to +1.1% |
| Fixed overlay matching | 12 | Committed pre-S3b baseline | −0.8% to +1.7% |

The common cases remain allocation-free; all rare-case allocation counts are
unchanged. One common-reading repeat drifted by up to 195% on the same binary;
it is retained, and a fresh complete common-group comparison agrees with the
earlier stable run. The table uses the median of three run medians for each
common case and the mean of two run medians elsewhere. Whole-operation repeat
drift reaches 2.9%; small differences do not establish universal speedups or
exact overhead. No workload average mixes the common and rare cases. The
OpenSMILES corpus matching runs remain excluded for the workload-parity reason
recorded above; no further investigation was performed.

Commands, source snapshots, immutable binaries and hashes, all timing runs,
allocation counts, and comparison CSVs are in scratch/236-s3b-completion.
verification.json, build-manifest.json, benchmark-manifest.json,
benchmark-summary.json, and common-repeat-summary.json identify the final
checks and results.

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
6. Resolution preserves every other assertion and the molecule-level list,
   without normalization, placement, or evaluation. Explicit removal policies
   retain their requested removals and the associated reference compaction.
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
13. Every constraint reading and the numeric projections listed in S3b agree
    exactly before and after successful molecule normalization, with the same
    entity and ring context. Every Some result is normal; None denotes bottom.
    Only contributing sources are interpreted, without rewriting storage or
    changing attribute-mutation fallibility.

## Copy-on-write storage

Molecule's graph, atom and bond tables, six overlay collections, and molecule-level
constraints are Arc-backed, as are Reaction's deltas and ReactionSpan's atom, bond,
constraint, and six overlay span collections. A clone increments reference counts,
and the first write to shared storage copies that storage.
[213](213-editor-overlay-storage-2026-08-27.md#remaining-costs) recorded the
Molecule constraint copy as a remaining cost; S5a removes it.

| Aggregate | Field | Current type | Shared on clone |
| --- | --- | --- | --- |
| Molecule | constraints | `Arc<Constraints>` | Yes |
| Reaction | lhs | Molecule | As Molecule |
| Reaction | deltas | `Arc<Deltas>` | Yes |
| ReactionSpan | graph | Graph | Yes |
| ReactionSpan | atoms, bonds | `Arc<Vec<EntitySpan<AtomForm>>>`, `Arc<Vec<EntitySpan<BondForm>>>` | Yes |
| ReactionSpan | six overlay span sets | Newtypes wrapping their relation set in Arc | Yes |
| ReactionSpan | constraints | `Arc<Vec<ConstraintSpan>>` | Yes |

ReactionSpan's default framed_eq clones both operands, and Python canonicalize
and tracked_canonicalize clone the wrapped span. These clones share all storage;
the subsequent transformations copy shared storage when they write or rebuild it.

Reaction application iterators share the reaction snapshot's delta storage after
S5b. Each iterator still owns its prepared Deltas and clones them for every match.
Each match reframes its own copy, so that copy remains.

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

### S5 production measurement (2026-10-13)

The aggregate_storage benchmark compares 888d7884c, immediately before S5a,
with 2793335a5, after S5d. Both revisions run the same benchmark source. This is
a fresh comparison of the twelve case categories above; the earlier experiment's
absolute timings are not the baseline for these results.

The fixtures are a six-carbon benzene with an aromatic overlay and a 100-carbon
chain. The constrained molecule carries 100 distinct molecule-level ChargeSum
constraints over individual atoms. Each reaction has one element-change delta;
the spans have unchanged entities. All forms are literal or undetermined.
Molecule and span construction measure from_entries; reaction construction
measures new. Criterion batches input preparation outside construction timing
and drops the results afterward. Clone timing includes dropping each clone while
retaining its source. No parsing is timed.

Four sequential runs used before/after/after/before order, 50 samples per case,
1-second warmup, 2-second measurement, and 20,000 resamples, with no concurrent
builds or tests. The ranges below are the two run estimates, not confidence
intervals; changes compare their means.

| Case | Before, two-run range | After, two-run range | Change |
| --- | ---: | ---: | ---: |
| Molecule construction, benzene | 735.24–753.21 ns | 732.48–741.06 ns | −1.0% |
| Molecule construction, 100-atom chain | 2.215 µs | 2.155–2.278 µs | +0.0% |
| Reaction::new, 100-atom chain | 48.19–48.45 ns | 32.71–33.47 ns | −31.5% |
| ReactionSpan construction, benzene | 4.159–4.206 µs | 4.099–4.153 µs | −1.4% |
| ReactionSpan construction, 100-atom chain | 17.766–17.990 µs | 17.819–17.990 µs | +0.1% |
| Molecule clone, benzene | 38.21–38.85 ns | 33.24–33.31 ns | −13.6% |
| Molecule clone, 100-atom chain | 38.29–38.81 ns | 33.25–33.32 ns | −13.7% |
| Molecule clone, 100 molecule-level constraints | 2.462–2.471 µs | 33.27–33.31 ns | −98.7% |
| Reaction clone, benzene | 69.40–70.28 ns | 36.45–36.95 ns | −47.5% |
| Reaction clone, 100-atom chain | 68.74–69.14 ns | 36.44–36.70 ns | −46.9% |
| ReactionSpan clone, benzene | 572.91–574.07 ns | 33.27–33.29 ns | −94.2% |
| ReactionSpan clone, 100-atom chain | 3.799–3.810 µs | 33.25–33.39 ns | −99.1% |

The clone reductions repeat across both comparisons. No construction slowdown
repeats across both: the molecule-chain result varies from +2.8% to −2.7%, so its
near-zero mean is not evidence that the added Arc allocation is free. These
microbenchmarks measure construction and clone/drop, not mutation after sharing
or end-to-end operation throughput. Raw Criterion samples, confidence intervals,
the baseline snapshot, binaries, commands, and run metadata are in
scratch/236-s5e/.

## Outside this document

Calling resolvers and transformers with the molecule as receiver, so that their
calls chain left to right, needs a separate design.

## Implementation plan

This plan closes after S7a, before 181's plan begins.
Both plans together are the 0.10 release blocker.

No production change starts before its subitem is written here. Subitem labels
remain stable; implementation follows document order and the explicit
dependencies, with the agreed exception that S3b1–S3b3 finish before returning
to S3a3–S3a4. S3c, S3f, and S3g specify parts of S3b, as marked
below. S3b is the Rust reading migration; S3b1–S3b3 sequence ring ownership,
ring bindings, and semantic-reading bindings after it.
S3d includes S4g's aromaticity mismatch-policy removal: reading contradictory
evidence as bottom cannot supply a ReplaceEntity candidate. S3e includes the
stereo mismatch-policy removal for the same reason.

Use focused checks for individual subitems. Run workspace tests after a set of
related changes and at the final gate.

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

#### S3a2 — Entity constraint getters return references **Done**

**Modules:** graph-IR constraint/atom, bond, dative, aromatic, multicenter,
noncovalent, and stereo, and view/constraints; their callers in graph-IR, graph,
io, and umol-py. **Breaking, red→green within the subitem.** [dep: none]

Every named getter on the eight constraint containers returns `&T`, the stored
value or a static vacuous value for an absent key:

- the getters that return `Option<&T>`: on atoms `valence`, `aromatic_valence`,
  `multicenter_valence`, `tetrahedral_stereo`, `degree`, `total_degree`,
  `total_valence`, `ring_degree`, `ring_valence`, `total_hydrogens`,
  `donated_pairs`, `accepted_pairs`, `ring_count`, and `ring_size_count`;
  `cis_trans_stereo` on bonds; and `ring_count` and `ring_size_count` on bonds
  and dative bonds;
- the getters that return owned values: `aromatic` on bonds and dative bonds,
  `electron_count` on aromatic systems and multicenter bonds, `intramolecular` on
  noncovalent bonds, and `topicity` and `stereogenicity` on stereo atoms and
  stereo bonds;
- `ligand_symmetry` and `fluxionality` on stereo atoms and stereo bonds, which
  return the whole form and instead return its `BooleanForm` value, as `topicity`
  and the ring getters return the value at their parameter.

`ring_memberships` on atoms, bonds, and dative bonds becomes public and yields
`&RingMembershipForm`, as the stereo iterators yield references to their forms.
The views' inherited getters follow and return `&'a T`; atom views gain
`ring_size_count`, the atom, bond, and dative bond views gain `ring_memberships`,
and the stereo atom and stereo bond views gain `ligand_symmetry`. Callers match
the value instead of an option: matches that cover both `None` and the vacuous
value collapse, `unwrap_or` with the vacuous value goes, and a caller that needs
an owned value clones it where it does. A stored vacuous value reads as an absent
key, so `infer_implicit_hydrogens` derives the valence of an atom that stores
`#v*` instead of inferring no count. The umol-py getters follow: those that return
`None` for an absent key return `Undetermined`, `ligand_symmetry` and
`fluxionality` return the value, as `topicity` already does, and the atom, bond,
and dative bond constraints gain `ring_memberships`, a list of ring membership
forms, as `ligand_symmetries`, `fluxionalities`, and `topicities` list theirs.

**Tests:** per getter, on containers and views, a stored value and an absent key,
which reads the vacuous value, including `ligand_symmetry` and `fluxionality` at a
permutation with and without an entry; `ring_memberships` on containers, views,
and the umol-py constraints yields each stored ring membership; existing results,
except the two `open_valence` inference cases, which expected no count for `#v*`
and are removed, since `filtered_nonliteral_h` under atom typing and
`unresolved_h_lp` under counts cover an atom without `#v`, and the umol-py tests
that read an absent key as `None`, which read `Undetermined`.

**Verification:** the resolve and substructure bench groups, against a baseline
taken before the change.

#### S3a3 — Vacuous stereo constraint matching **Done**

**Modules:** graph-IR constraint/stereo. **Green.** [dep: S3a2]

In the shared stereo constraint-container matches implementation, compare each
ligand-symmetry entry's invariant against the target's ligand_symmetry getter
at that permutation, and each fluxionality entry's active value against the
target's fluxionality getter. Absence then reads as top, as it already does for
topicity and stereogenicity. Keep matching borrowed and pattern-driven.

**Tests:** on both containers, cover all four keys with a vacuous pattern and
an absent target entry; a target carrying only another parameterized key;
matching and conflicting determined values; and a determined pattern against
an absent or vacuous target. Check exact booleans against the meet-based
definition and retain the existing normalization cases.

**Completion (2026-10-12):** The shared matches implementation reads the
target's ligand-symmetry invariant and fluxionality active value through the
borrowed getters at the pattern's full key. An absent key therefore contributes
top, including when the target contains only another permutation, orientation,
or ligand pair. Matching remains pattern-driven and adds no cloning,
allocation, or normalization. Public signatures and storage are unchanged.

Both containers gain 28 exact cases covering all four key families, with
checks against the meet definition and normalized forms in the matching tables.
The new raw-input cases exposed ten failures before the correction; all 107
targeted constraint-container matching cases pass afterward. Existing matching
and normalization cases are retained.

**Verification:** Workspace tests with all features pass (37,001 tests, ten
ignored), including the stereo lattice properties and normalization cases.
The full Python suite passes (2,204 tests, two skipped) after rebuilding under
Python 3.13.15. Workspace Clippy with all targets/features and warnings denied,
formatting, and diff checks also pass. Logs are in scratch/236-s3a3.

#### S3a4 — Borrow stereo validation entries **Done**

**Modules:** graph ops/validate/stereo. **Green.** [dep: S3a2]

Pass the topicity and ligand-symmetry iterators directly from both stereo
validation callers. Change the private validate_symmetry parameters to accept
iterators of references, removing the two collected Vecs per entity. Clone a
non-Copy asserted payload only when constructing a contradiction that owns it.
Preserve validation order and all existing outcomes.

**Tests:** migrate the private helper's existing cases to borrowed iterators;
run the public stereo-conformance cases for determined, underdetermined, and
contradictory results, including invalid topicity positions and permutation
degrees. No new public API or result type.

**Completion (2026-10-12):** Both stereo validation callers pass their borrowed
topicity and ligand-symmetry iterators directly to the private validate_symmetry
helper. This removes two Vec collections and their entry clones per entity.
The validation loops, order, and outcomes are unchanged; non-Copy asserted
payloads are cloned only when constructing an owned contradiction. Public APIs
and stored constraints are unchanged.

**Verification:** All 44 focused stereo-validation cases pass before and after
the refactor. Existing helper cases now use borrowed iterators. Five added
public cases cover underdetermined configurations and topicity assertions for
both entity kinds, plus preservation of a non-Copy assertion in a contradiction.
The rebuilt Python suite passes (2,204 tests, two skipped) under Python 3.13.15.
Workspace Clippy with all targets/features and warnings denied, formatting,
and diff checks pass. The workspace test run was stopped at the user's request
to reserve it for batches of related changes; it is not a completed gate for
this subitem. Logs are in scratch/236-s3a4.

#### S3b — Constraint readings as meets **Done**

**Completed (2026-10-12).** The six semantic readings return normalized owned
forms with None as bottom, and the
ten public numeric projections return Option<NumForm> under the same convention.
The signature and closure contract above settles invariant 13 and the complete
readings' source boundaries. The optimized trial and retained caller changes
are accepted with their recorded measurement limits. Storage remains raw and
attribute mutation infallible. The shared standalone-container implementation,
diagnostic formatting, and Rust verification are complete. The settled
Python interface follows in S3b2–S3b3 after the ring-ownership change in S3b1.

**Modules:** graph-IR entity constraint containers, view/constraints and the entity
view modules, and their callers in graph-IR, graph, io, and umol-py;
umol-graph-ir tests/property;
docs/development/nomenclature.md. **Breaking, red→green within the subitem.**
[dep: S3a, S3a1]

Give the eight constraints views the six readings of the constraint-readings
decision, semantically meeting their sources from undetermined(key) and
returning Option<F>, with None as ⊥ and the vacuous form as no information.
This specifies the result, not a requirement to invoke meet on every source:
handle common Undetermined and Lit inputs directly and move owned intermediates.
Interpret projection inputs before their information is lost in arithmetic.
`asserted`, `asserted_complete`, `derived`, and `derived_complete` change from
`Option<&F>` and `Option<F>` with `None` for absence; `evidenced` and
`evidenced_complete` are added. asserted_complete closes only the own assertion
from this subitem onward. Include the caller migrations detailed under S3c and
S3f and the closure change under S3g in this breaking unit; those labels retain
their scope below but are not subsequent intermediate contracts. Migrate the
remaining consumers listed in the signature contract, including diagnostics.

Finish the shared asserted and asserted_complete methods on all eight Rust
constraint containers before binding them. Have the molecule-backed views
delegate to these methods, preserving the accepted common-variant fast paths.
S3b3 binds these shared implementations through each Python constraints view's
backing-storage enum after S3b2 supplies the ring-context interface.

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

Verify container/view equality for both assertion readings. Check the optional
numeric projection results and diagnostic display without substituting zero
for an unavailable literal. S3b2–S3b3 carry the Python cases.

**Verification:** the resolve and substructure bench groups, against a baseline
taken before the change. Preserve the selected workloads' variant distribution;
report rare set/expression stress cases separately. Include projection
invariance, bottom propagation, and unchanged infallible mutation in correctness
coverage; performance acceptance does not replace these contracts.

#### S3b1 — Owned and borrowed ring data in RingViews **Done**

**Modules:** graph-IR ir/molecule, ir/ring, and ir/view/ring and their tests.
**Additive, green.** [dep: S3b]

Store Cow<'a, RingSet> inside RingViews alongside its molecule borrow. Preserve
Molecule::rings(model, config) and crate-private RingViews::new, which consumes
the computed set into Cow::Owned. Add public RingViews::from_ring_set(&RingSet,
&Molecule) using Cow::Borrowed. Keep into_ring_set and implement it through
into_owned, documenting the move for owned data and clone for borrowed data.
Preserve the existing query implementations and use the same methods for both
ownership modes. Enumeration occurs only at the explicit perception call.

Guard contextual access so mismatched sets and molecules cannot cause a panic,
without promising meaningful results or adding full topology validation.
Existing owning consumers, including ConstraintEvaluation, keep their current
storage and call signatures. The general view-argument cleanup is not part of
this subitem.

**Tests:** retain exact expected members, ids, relations, and readings for the
existing owning route; exercise reconstruction and extraction for the borrowed
route against independent expected values. Reuse one RingSet across a bond-order
mutation, recreating the borrowed view: membership and sizes stay the same and
ring valence reads the new order. Verify standalone topology queries after the
source molecule is dropped and independence of a set extracted from a borrowed
view after the source set is dropped. Cover mismatched contextual inputs without
requiring a semantic result. Run the existing matching, fingerprint, aromaticity,
and constraint-validation suites that exercise owning ring context.

**Completed 2026-10-12.** RingViews now stores Cow internally and exposes the
borrowed reconstruction method. The two molecule-dependent ring readings guard
absent atom ids before accessing adjacency. Existing owning consumers and query
signatures are unchanged. Tests check borrowed storage identity, ownership-aware
extraction, current attribute readings, and incompatible context; two new
properties cover reconstruction equivalence and panic freedom under independent
context.

**Verification:** all 55 ring-view unit cases and four ring properties pass.
The workspace test suite with all features passes (36,943 tests, 10 ignored),
as do workspace Clippy with all targets/features and warnings denied, graph-IR
rustdoc with warnings denied, formatting, and diff checks. Logs are in
scratch/236-s3b1.

#### S3b2 — Python ring data and views **Done**

**Modules:** umol-py ring and molecule bindings, extension registration, Python
package exports, and ring tests. **Additive, green.** [dep: S3b1]

Bind the ring family inventoried in the ownership decision, including
Molecule.rings(model, config) -> RingViews,
RingViews.from_ring_set(rings, molecule) -> RingViews, and
RingViews.into_ring_set() -> RingSet. Preserve the distinction between owned
and borrowed ring data and the consuming extraction's move-versus-clone costs.
Retain the required Python owners and reconstruct Rust borrows per call. Bind
the existing ring-set, ring-graph, individual-ring, ring-atom, and ring-bond
queries, including RingAtomView.ring_valence -> NumForm | None. Keep borrowed
collection access read-only and iteration lazy. Match Rust construction and
field access, and use the existing enum, integer-id, indexing, and exception
conventions. RingSet and RingGraph acquire no arbitrary-parts constructors;
the crate-private consuming RingViews::new is not a Python constructor.

**Tests:** exact ring membership, sizes, relations, graph connections, and
atom/bond readings on representative ring configurations; indexing and absence;
owner retention after dropping caller references; current attribute reads;
structural invalidation and molecule consumption for nested views and iterators;
and continued standalone RingSet access. Cover both ownership routes through
into_ring_set, including consumption and dependent-accessor invalidation, and
reconstruction with incompatible context without a panic. Distinguish bottom
from top in ring valence. Verify that creating or partially advancing an iterator
does not materialize the remaining entries.

**Completion (2026-10-12):** the ten ring types and Molecule.rings are bound and
exported. Reconstruction retains the RingSet and Molecule owners without
copying or perception. Consuming extraction moves an owned set or clones a
borrowed one; storage-identity tests verify both costs. Rust's writable record
fields remain writable, while borrowed graph edges and ring-member sequences
are read-only. RingView.id is a writable label, independent of its member
slices, matching Rust.

Consumed RingViews raises ConsumedError; its dependent views, sequences, and
iterators raise InvalidatedViewError. Molecule consumption and structural
operations invalidate the entire molecule-backed view tree, including exhausted
iterators. Attribute changes remain visible, and independent RingSet access
remains available. All exceptions use the existing shared types.

Collection adapters retain owners and convert only requested entries. Atom/bond
rings() reconstructs the Rust filtered iterator and advances with nth, so each
step rescans its visited prefix; it neither buffers the remaining rings nor
retains a Rust borrow across calls. No Rust public API was added for this
binding adaptation.

**Verification:** 71 Python ring cases pass, including ownership, lifetime,
invalidation, incompatible context, lazy access, and bottom-versus-top readings.
The full Python suite passes (2,057 tests, two skipped) after rebuilding under
Python 3.13.15. Workspace tests with all features pass (36,945 tests, ten
ignored), as do workspace Clippy with all targets/features and warnings denied,
formatting, and diff checks. Logs are in scratch/236-s3b2.

#### S3b3 — Python constraint readings and numeric projections **Done**

**Modules:** umol-py entity constraint containers and views, atom/bond bindings,
and Python tests. **Additive, green.** [dep: S3b, S3b2]

Bind asserted and asserted_complete on all eight constraint containers and the
six readings on their constraints views, using each family's existing key and
form types. Bind view satisfies and is_compatible against the corresponding
constraint container. Dispatch through the existing backing-storage enum and
the shared Rust implementations, with the settled ValueError boundaries.

Bind AtomConstraintsView.with_rings and BondConstraintsView.with_rings; the new
view retains the RingSet owner and the same molecule entity. Add AtomView's nine
optional numeric projections: valence, donated_pairs, accepted_pairs,
total_degree, total_hydrogens, total_valence, covalence, heavy_atom_degree, and
heavy_atom_valence. Together with S3b2's ring_valence, these are the ten settled
NumForm | None projections. Typed storage getters keep their raw meaning.

**Tests:** all eight entity families; standalone and molecule-backed assertion
readings; ValueError for missing molecule or required ring context, including
the settled empty-comparison case; successful explicit ring-context reads;
bottom versus top; optional numeric results; and unchanged live-view mutation
and invalidation. Verify Rust/Python reading parity and independent returned
forms without introducing a Python-only reading implementation or a separate
storage-proxy class.

**Completion (2026-10-12):** All eight Python constraint containers expose the
two own-assertion readings, and their views expose all six readings plus
satisfies and is_compatible. They delegate to the shared Rust implementations;
None denotes bottom and undetermined denotes top. Standalone assertion reads
work without molecule context; graph-dependent operations reject missing
context with ValueError, including empty comparisons. Atom/bond ring keys
require explicit ring context for derivation, evidence, and comparisons.

with_rings retains the supplied RingSet owner and the same molecule entity
without copying. The resulting view shares the existing invalidation contract;
attribute and constraint mutation remain visible. Returned readings are
independent normal values, and storage getters retain their raw meaning.
AtomView's nine numeric properties complete the ten projections with S3b2's
RingAtomView.ring_valence. No additional Python types are introduced.

**Verification:** 147 new Python cases cover all eight entity families,
normalization without storage rewrites, bottom/top, closure, scoped keys,
missing context, explicit ring ownership, live mutation, invalidation, and the
numeric projections. The full Python suite passes (2,204 tests, two skipped)
after rebuilding under Python 3.13.15. Workspace tests with all features pass
(36,945 tests, ten ignored), as do workspace Clippy with all targets/features
and warnings denied, formatting, and diff checks. Logs are in scratch/236-s3b3.

#### S3c — Valence admission and implicit hydrogen inference read evidenced_complete **Done**

**Completed with S3b.** Modules: graph valence/atom_typing and
valence/counts. Retained here to specify its behavior and tests.

Atom typing matches each registry-row entry against the atom's
`evidenced_complete` value, replacing the inline composition, and keeps a ground
atom under resolution when that reading's aromatic valence is `#a+`, replacing
the system-membership test. It builds the ring set when a row states a ring key.
Counts applies the same ground-atom test and reads aromatic evidence and the
candidate aromatic valences from `evidenced_complete`. An atom whose evidence
meets to ⊥ admits no candidate and reports NoMatch.

`infer_implicit_hydrogens`, under both sources, reads the atom's valence and
aromatic valence from `evidenced_complete`, replacing the rule that a stored
assertion takes precedence over the derived value. An assertion that conflicts
with its derived value meets to ⊥ and infers no count; an incident bond's `#a`
reads as `#a+` and infers no count, as the atom's own `#a+` does. Both doc
comments follow.

**Tests:** existing admission and export results; NoMatch under both sources for
an atom that asserts non-aromaticity next to a bond that asserts aromaticity,
and under counts for an atom that asserts non-aromaticity inside a stored
system; a registry row with a ring key admits by the atom's ring membership.
For `infer_implicit_hydrogens` under both sources: the stored-valence and
stored-aromatic cases, whose assertions conflict with their derived values,
infer no count; an atom outside the stored systems next to a bond that asserts
aromaticity infers no count; under counts, the cases that assert a nonzero `#v`
on an atom without bonds, and the table case, gain bonds that realize it and
keep their results.

**Verification:** the resolve bench groups resolve_ingest_smiles and
smiles_roundtrip/resolve, against a baseline taken before the change.

#### S3d — Perception reads typed contributions from evidenced_complete **Done**

**Modules:** graph aromaticity, its rule modules, resolve/aromaticity, resolve,
and transform/aromatizer; umol-py resolve bindings, exports, and tests.
**Breaking, red→green within the subitem.** [dep: S3b]

Remove AromaticityMismatchPolicy and AromaticBondConstraintMismatchPolicy,
including the aromatic_valence_mismatch and aromatic_bond_constraint_mismatch
config fields, recovery branches, Python classes, keywords, getters, and repr
fields. Every aromaticity mismatch is a contradiction. This portion of S4g
moves here because bottom evidence cannot produce a replacement candidate.

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
Cover mismatch contradictions in standalone and composite resolution; remove
tests of the retired recovery policies and retain the Python config parity tests.

**Completion (2026-10-12):** Perception callbacks, completion selection, and
electrons_from_kekule now use AromaticValence. NotAromatic excludes an atom;
Aromatic(0) remains eligible. Numeric consumers use infallible casts and
ordinary arithmetic under the repository size assumption; integer-representation
bounds do not determine aromatic eligibility. derive reads complete evidence
once per atom; selection does so for atoms without completions. Both retain
scalar contributions for subsequent perception calls.
Positive evidence with a nonliteral count leaves the result Underdetermined;
bottom reports the corresponding system or localized-bond mismatch. Stored
electron-vector length failures retain their existing classification.

Both aromaticity mismatch policies and their recovery branches are removed
from Rust and Python, including config fields, exports, constructor keywords,
getters, and repr fields. The remaining failure policies cannot permit a
mismatch. The stereo policy removal is covered by S3e.

**Verification:** All focused checks pass: 331 aromatic unit cases, 187
composite-resolver cases, 13 aromatizer cases, the exhaustive-selection
reference property, 34 aromatic conformance fixtures, and 27 matching and
kekulization integration cases. The missing-evidence regression fails before
the change and passes afterward; its explicit zero-contribution counterpart
remains aromatic. Tests also cover evidence refinement, both new and existing
mismatch paths, and unchanged inputs after rejected resolution. The reference
property retains its generator and search law.

After rebuilding with Python 3.13.15, all 60 Rust binding cases and 225 focused
Python cases pass. Clippy for umol-graph and umol-py with all targets/features
and warnings denied, formatting, and diff checks pass. Workspace tests were
not rerun for this subitem. Logs are in scratch/236-s3d.

After removing the numeric conversion guards, the aromatic unit cases,
selection reference property, aromatic conformance fixtures, and umol-graph
Clippy with all targets/features pass again.

**Performance:** Recorded before removing the numeric conversion guards.
A before/after comparison of the existing resolve_ingest_smiles
aromatic cases used the s3d-before baseline, 20 samples, one-second warm-up,
and two-second measurement per case. No builds or tests ran concurrently.

| Case | Criterion change | 95% interval |
| --- | ---: | ---: |
| benzene | -0.20% | -0.54% to +0.09% |
| pyridine | -0.65% | -1.56% to +0.04% |
| quinoline | -0.33% | -0.77% to +0.06% |
| purine | +1.92% | +1.08% to +2.87% |

Criterion detects no significant change for the first three cases and flags
a regression for purine. This single short comparison does not isolate its
cause or establish performance beyond these cases. Both benchmark logs are
in scratch/236-s3d.

#### S3e — Stereo perception reads evidenced_complete **Done**

**Modules:** graph-IR stereo, view/atom, view/bond, and view/constraints; graph stereo,
resolve/stereo, and resolve; umol-py resolve bindings, exports, and tests.
**Breaking, red→green within the subitem.**
[dep: S3b]

Read each site's tetrahedral or cis-trans evidence from `evidenced_complete`. A
stereo value derives or confirms the site's entity, ⊥ is the existing mismatch,
and the negative leaves the site without an entity.
Independent constraint and entity realizability failures retain their existing
diagnostics and failure policies. Bottom in an own assertion without an entity
is a constraint failure.
The resolver's partial-input gate reads the same evidence, so a stored entity
can refine an open assertion before that gate classifies it as underdetermined.

Transport stored stereo configurations into the assertion frame in the shared
graph-IR projections before meeting them with assertions. Actual ligands follow
neighbor order, with virtual ligands last; cis-trans uses that order within each
endpoint block. Exchanging whole endpoint blocks preserves cis-trans configuration.
The previous projections copied the stored coset without transport and could report
false bottom when stored ligands are reordered. This prerequisite fix is included
in S3e and also corrects the Python readings.

Remove StereoMismatchPolicy, the tetrahedral_stereo_mismatch and
cis_trans_stereo_mismatch config fields, recovery branches, and corresponding
Python classes, keywords, getters, repr fields, and exports. Every mismatch is
a contradiction. This remaining portion of S4g moves here because bottom
evidence cannot supply a ReplaceEntity candidate.

**Tests:** existing stereo derivation and resolution results, including failures
and mismatches; frame transport in shared readings, including virtual ligands;
refinement of open assertions by stored entities; unconditional
mismatch rejection in standalone and composite resolution; Python config parity.

**Completion (2026-10-12; performance accepted):** perception reads complete evidence, uses site lookup
instead of scanning every stereo entity, and builds each candidate frame once.
The partial-input gate allows a stored entity to refine an open assertion.
Ground assertions or configurations skip that reading because their evidence
is ground or bottom; absent assertions and entities also cannot produce partial
evidence. For consistent tetrahedral evidence, entity validation checks the
distinct ligand set without transporting the configuration again. Cis-trans
retains its block-aware check, with the existing same-frame shortcut.
Independent entity checks still reject an unrealizable virtual ligand without
inventing a replacement entity. The small-ring assertion exclusion is retained.
StereoMismatchPolicy and its Rust/Python config and export surface are removed,
completing S4g.

The shared projections now transport through a four-position stack array and the
existing coset action, without rewriting storage. Six frame regressions and two
open-assertion refinement regressions failed before their fixes and pass after.
Molecule construction and attribute mutation retain their existing contracts.

**Verification:** 280 graph stereo unit cases (including composite resolution
and validation), 144 graph-IR constraint-view cases, three constraint-reading
properties including reframing invariance, 36 stereo resolution conformance
cases, two stereo projection properties, 50 Rust binding cases, and 373 focused
Python cases passed after rebuilding the extension with Python 3.13.
Clippy passed with warnings denied for graph-IR, graph, and Python across all
targets, including the affected property/conformance features. Formatting and
git diff --check passed. No workspace tests were run for this subitem.

The four existing smiles_roundtrip/resolve stereo cases were compared with the
s3e-before baseline using 20 samples, one-second warm-up, and two-second
measurement, with no simultaneous builds or tests:

| Case | Time change | 95% interval |
| --- | ---: | ---: |
| tetra_four | +2.49% | -0.67% to +5.06% |
| tetra_ring | +5.61% | +2.63% to +7.90% |
| alkene_four | +4.68% | -2.12% to +11.01% |
| shared_chain | +4.31% | +3.65% to +4.98% |

Criterion flags regressions for tetra_ring and shared_chain and detects no
significant change for the other two. The baseline had outliers; this short
comparison records the measured cost of the subitem without isolating its
cause or establishing broader performance. Logs are in scratch/236-s3e.

The user accepts approximately 1% overhead but finds 5% too costly. Performance
acceptance was reopened to attribute and reduce these regressions. The user
accepted the optimized implementation and the measurements below on 2026-10-12,
completing S3e.

**Performance attribution (2026-10-12):** rebuilt the committed baseline
cb63b560b773de7be6b2af119564854e7ca516d5 and saved separate executables for
the original S3e implementation and a diagnostic variant with its consumer
changes but the old projections. The diagnostic variant deliberately lacks the
frame correction and is only a cost-isolation experiment. Eight interleaved
rounds of a separate public-API harness gave these median times, in microseconds:

| Case | Baseline | S3e consumers, old projections | Original S3e |
| --- | ---: | ---: | ---: |
| tetra_four | 7.797 | 7.725 | 8.096 |
| tetra_ring | 12.712 | 13.096 | 13.503 |
| alkene_four | 11.253 | 11.452 | 11.671 |
| shared_chain | 13.387 | 13.598 | 13.920 |

This diagnostic timing includes input clone and result destruction, unlike the
existing Criterion resolve group. It isolates the projection change at roughly
0.37–0.41 microseconds for the tetrahedral cases and 0.22–0.32 for the cis-trans
cases. The remaining consumer cost includes evidence reads on sites without
stereo and repeated readings in the partial-input gate. A branched molecule
without stereo gained approximately 0.31 microseconds from those consumer changes.

The tetrahedral ring has noncanonical stored ligand order. Original S3e
transported its configuration in the partial-input gate, again in perception's
evidence reading, and again in the independent entity check. Only the evidence
reading needs that transport: a ground configuration cannot yield partial
evidence, and consistent tetrahedral evidence leaves the entity check with the
ligand-set comparison. Other inputs also paid general coset transport for an
identity permutation. None of this rewrites the stored molecule.

The implementation now skips those redundant reads and transports. Literal
coset normalization returns the existing normal value, and literal/literal
meet compares the values directly; both avoid temporary sets. Identity
transport retains the existing invalid-index and degree behavior. Complete
evidence semantics, partial forms, mismatch diagnostics, and infallible attribute
mutation are unchanged.

Allocation counts from one warmed resolution, excluding input construction:

| Case | Baseline | Original S3e | Optimized S3e |
| --- | ---: | ---: | ---: |
| tetra_four | 45 | 43 | 41 |
| tetra_ring | 67 | 67 | 63 |
| alkene_four | 99 | 99 | 88 |
| shared_chain | 134 | 134 | 112 |

The original regression therefore did not come from increased total allocation
counts. Removing the temporary literal sets additionally improves existing
work in the closure checks. The complete stereo reading sweep allocates zero
times on these four resolved inputs after optimization.

Source snapshots, saved executables, build records, allocation counts, and
interleaved timing logs are in scratch/236-s3e/attribution.

**Final performance comparison:** three repetitions of the four existing
Criterion resolve cases, using separate saved baseline and final executables,
40 samples, one-second warm-up, and two-second measurement. Baseline and final
runs were adjacent for each case; the second repetition reversed both variant
and case order. No builds or tests ran concurrently. These are Criterion's
reported change estimates, with the reversed comparison expressed as final
relative to baseline; negative means faster:

| Case | Run 1 | Run 2 | Run 3 | Median estimate |
| --- | ---: | ---: | ---: | ---: |
| tetra_four | +0.58% | -2.15% | -6.26% | -2.15% |
| tetra_ring | -2.11% | -3.56% | -0.81% | -2.11% |
| alkene_four | -2.26% | -3.18% | -1.81% | -2.26% |
| shared_chain | -1.26% | +0.24% | -4.52% | -1.26% |

The repeated 4–6% regression is no longer observed. All positive point estimates
are below 1%, but this is not a strict 1% statistical bound: the noisy tetra_four
first run has a 95% interval of -2.39% to +5.80%, and shared_chain's second run
has -1.83% to +2.96%. The other comparisons show improvements. The median column
summarizes the three estimates; it is not a pooled confidence interval, and
these four cases do not establish performance on other workloads.

**Follow-up verification:** 3,081 graph-IR stereo unit cases, 282 graph stereo
unit cases, 144 constraint-view cases, 115 lattice property cases, three
constraint-reading properties, 36 stereo conformance cases, and two stereo
resolver properties passed. After rebuilding the extension with Python 3.13,
254 constraint-reading and resolver Python cases passed. Clippy passed with
warnings denied for graph-IR, graph, and Python across all targets and the
affected property/conformance features. Formatting and git diff --check passed.
No workspace tests were run.

#### S3f — Matching reads evidenced_complete **Done**

**Completed with S3b.** Modules: graph-IR view/constraints and the
substructure tests; docs/development/nomenclature.md. Retained here to specify
its behavior and tests; requires the completed S3a matching contract.

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

#### S3g — asserted_complete without the incident-bond merge **Done**

**Completed with S3b.** Modules: graph-IR view/atom and the view tests;
docs/development/nomenclature.md. Retained here to specify its behavior and tests.

The atom's `asserted_complete` closes `asserted` alone, where it is top, to the
absence cell, and no longer reads incident bonds. Move the guide's sentence on
the two aromatic placements to the evidenced readings.

**Tests:** an atom whose only aromatic evidence is an incident bond reads
non-aromatic under `asserted_complete`; the existing cases with own assertions.

#### S3h — Raising writes #a! only next to an aromatic bond **Done**

**Modules:** umol-io table_ir/raise and its tests; graph ingest tests.
**Green.** [dep: S3b, S3d]

Raise an atom that the table marks non-aromatic with `#a!` only when one of its
table bonds is aromatic.

**Tests:** raising `CCO` writes no aromatic valence; raising `c1ccccc1` writes
`#a+` on every atom and the aromatic constraint on every bond; raising an
aliphatic atom joined by an explicit aromatic bond writes `#a!` on it, and
ingesting that input reports NoMatch; resolved ingestion results are unchanged.

**Completion (2026-10-12):** molecule raising adds the negative assertion during
its existing bond traversal, only for a table atom marked non-aromatic and
incident to a definite aromatic table bond. Standalone atom raising has no bond
context and writes only positive aromatic assertions. No additional adjacency
structure is built. Invalid bond references still report MoleculeEntries errors.

The 267 raising tests, 242 ingest tests, and three focused raising properties
pass, including exact unchanged ethanol and benzene results and C:C reporting
CountsError::NoMatch. Scoped Clippy with warnings denied, formatting, and diff
checks pass. Workspace tests remain for the related-change batch gate.

The existing smiles_roundtrip/raise cases were compared against commit
204b66a5b with 40 samples, one-second warm-up, and two-second measurement,
without concurrent builds or tests:

| Case | Time change | 95% interval |
| --- | ---: | ---: |
| chain_64 | -6.49% | -6.77% to -6.21% |
| branched | -2.88% | -3.60% to -1.76% |
| aromatic_fused | -0.73% | -1.21% to -0.26% |

Criterion reports improvements for the aliphatic cases and a change within its
noise threshold for aromatic_fused. These measure raising only. The baseline
executable and logs are saved in scratch/236-s3h.

#### S3i — Aromaticity discharges what its systems realize **Done**

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

**Completion (2026-10-12):** system addition now removes each member's own
aromatic-valence assertion and writes no bond assertion. Standalone aromaticity
planning and composite constitution planning share the bond discharge planner:
outside pre-existing systems, it removes assertions satisfied by membership in
the accepted systems, including negative assertions on linking bonds. It uses
the existing derivation or selection without perceiving again. It allocates the
membership lookup only when a relevant bond assertion exists. Standalone
resolution preserves member and internal-bond assertions of stored systems.
Discharge preserves unrelated constraint keys and inherent fields.

Removed reset_aromatic_valence from the Rust config and its Python constructor,
getter, and repr, and from the nomenclature guide. The shared planner adds no
public API or failure variant. Composite resolution still has its closing
discharge pass until S3m.

Both new discharge fixtures failed before implementation. All 1,999 graph unit
tests, four resolution properties, 683 resolution conformance cases, and 50
Rust binding tests pass. The generated discharge property covers separate or
linked rings with mixed stored/new systems, absent/vacuous/positive internal
bond assertions, preservation of unrelated keys, and idempotence. Composite
conformance snapshots are unchanged. After rebuilding with Python 3.13, the
affected Python suites pass with 407 tests and two existing skips. Scoped
Clippy with warnings denied, formatting, and diff checks pass. Workspace tests
remain for the related-change batch gate.

The existing smiles_roundtrip/resolve cases were compared against commit
6a47de37f with 40 samples, one-second warm-up, and two-second measurement,
without concurrent builds or tests:

| Case | Initial time change | 95% interval |
| --- | ---: | ---: |
| chain_64 | +2.39% | +2.03% to +2.79% |
| aromatic_fused | -5.62% | -8.03% to -3.15% |
| aromatic_lone_pair | -4.04% | -4.99% to -3.01% |

Short repeats of the saved executables exposed drift: a forward chain
comparison detected no change, while a reverse comparison favored the baseline
by 6.40%. Four further pairs, alternating execution order, gave a change of
+0.49% in the average of their mean times (82.21 to 82.61 µs), with individual
pairs ranging from -0.99% to +1.58%. These controls do not reproduce the larger
chain slowdown, but do not establish a strict one-percent bound. The new bond
scan performs no allocation on this assertion-free input. The aromatic cases
avoid the former per-system bond collection and assertion writes. Executables,
logs, and alternating-run estimates are saved in scratch/236-s3i.

#### S3j — Stereo discharges what its entities realize **Done**

**Modules:** graph resolve/stereo, resolve, and stereo; umol-py resolve bindings
and Python tests; docs/development/nomenclature.md. **Breaking, red→green within the subitem.**
[dep: none]

The stereo plan removes the tetrahedral or cis-trans constraint at each site
where it adds an entity. Remove the partial-input gate so structurally
realizable positive evidence can produce an entity with an open configuration;
the composite resolver's final concreteness check remains. Remove
`reset_stereo_constraints` from StereoResolveConfig and from the Python
binding's keyword, getter, and repr.
Remove the guide's Reset entry with its last code names.

**Tests:** standalone resolution removes the constraint at each site that gains
an entity, including one whose configuration is undetermined, and keeps it at
a site with an existing entity. Vacuous stereo assertions add no entity. The
composite resolver uses this plan and preserves its existing concrete results;
open configurations still give Underdetermined with unchanged input. Composite
preservation of existing-site assertions follows when S3m removes the closing
discharge pass.

**Completion (2026-10-12):** stereo planning unconditionally removes the
tetrahedral or cis-trans assertion at each newly added entity's site, using
direct constraint-removal edits. Retained entities do not trigger discharge;
other keys and inherent fields are preserved by discharge. The existing failure
policies and small-ring assertion exclusions remain.

The partial-input gate prevented the planned open-configuration case from
reaching perception. Its removal was accepted during implementation:
standalone resolution may produce Th* or Ct* from positive open evidence and
discharge the assertion without choosing a coset. Vacuous assertions add no
entity, and existing open entities are retained. Unrealizable open evidence now
reaches the same realizability checks and failure policies as literal evidence.
The composite resolver's final concreteness check and preservation of the input
on Underdetermined are unchanged. Its closing discharge pass remains until S3m.

Removed reset_stereo_constraints from the Rust config and Python constructor,
getter, and repr, and removed the nomenclature guide's Reset entry. No public
symbol or error variant was added.

The new discharge and open-configuration fixtures failed before their changes.
All 2,005 graph unit tests, five resolution properties, 683 resolution conformance
cases, and 50 Rust binding tests pass. The discharge property compares exact
independent results for mixtures of stored/new tetrahedral and cis-trans
entities with literal or undetermined configurations, including unrelated-key
preservation and idempotence. Composite conformance snapshots are unchanged.
After rebuilding with Python 3.13, all 287 focused Python tests pass. Scoped
Clippy with warnings denied, formatting, and diff checks pass. Workspace tests
remain for the related-change batch gate.

The existing smiles_roundtrip/resolve cases were compared against commit
83e1f2ae5 with 40 samples, one-second warm-up, and two-second measurement,
without concurrent builds or tests:

| Case | Time change | 95% interval |
| --- | ---: | ---: |
| chain_64 | -0.24% | -0.73% to +0.25% |
| tetra_four | -0.98% | -1.69% to -0.30% |
| alkene_four | -11.41% | -11.99% to -10.83% |
| shared_chain | -11.79% | -13.49% to -10.60% |

Criterion detects no change for chain_64, reports tetra_four within its noise
threshold, and reports improvements for the two cis-trans cases. Removing the
gate eliminates its preliminary atom/bond scans and the duplicate
skipped_stereo_bonds query, including ring enumeration when needed. These
measurements cover the listed literal-configuration inputs. Executables and
logs are saved in scratch/236-s3j.

#### S3k — Valence discharges #v **Done**

**Modules:** graph resolve. **Green.** [dep: none]

The composite resolver's constitution plan removes `#v` from each atom whose
fields it completes.

**Tests:** the constitution plan of an atom asserting `#v` that it completes
removes the assertion, and that of an atom it does not complete leaves it.

**Completion (2026-10-12):** the shared constitution planner removes the valence
assertion from each selected atom completion before computing its attribute
update. The existing edit projection emits the checked removal. Other stored
assertions and atoms without completions are preserved by this plan. Admission
already checks valence assertions; no additional validation pass or public API
change is needed. The closing discharge pass remains until S3m, so final
composite preservation of assertions on untouched atoms still awaits that step.

Seven exact planner cases cover literal, set-valued, explicitly vacuous, and
absent valence assertions; preservation of other keys and an untouched atom in
the same molecule; and empty, complete, and uncompleted inputs without selected
completions. The three discharge cases failed before the implementation change.
All 2,012 graph unit tests, five resolution properties, and 683 resolution
conformance cases pass; snapshots are unchanged. Scoped Clippy with warnings
denied, formatting, and diff checks pass. Workspace tests remain for the
related-change batch gate.

The existing smiles_roundtrip/resolve cases were compared against commit
d60dab38b with 40 samples, one-second warm-up, and two-second measurement,
without concurrent builds or tests:

| Case | Time change | 95% interval |
| --- | ---: | ---: |
| chain_64 | -1.97% | -5.41% to +0.93% |
| aromatic_fused | -0.63% | -3.11% to +1.45% |

Criterion detects no performance change in either case. These inputs carry no
valence assertions, so the comparison checks overhead on ordinary resolution
paths, not discharge savings. Executables and logs are saved in scratch/236-s3k.

#### S3l — Aromaticity planning makes the closing pass's checks **Done**

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

Explicit non-aromatic bonds delimit candidates for HMO as well as Hückel and
Clar. Candidate construction may recover several systems within a connected pi
component; each candidate is evaluated separately under the selected rule.
Any combined candidate that still contains a non-aromatic bond is rejected,
without searching alternative partitions. This conservative handling of the
overlap edge case is accepted for all three rules.

**Tests:** each check in standalone and composite resolution under the Error and
Keep failure policies; a benzene ring with one bond asserting non-aromaticity
forms no system, so its aromatic atoms are an AromaticValenceFailure under Error
and keep their assertions under Keep; naphthalene with such a bond on the
periphery of one ring forms a system on the other ring; under each rule, two
benzene rings joined by a bond asserting non-aromaticity form separate systems;
a coronene hub-bond exclusion rejects the conflicting combined region while a
separate benzene remains perceived; a property over the skeleton generator, for
invariant 11: no system the resolver adds contains a bond asserting
non-aromaticity or an atom whose evidence disagrees with its contribution.

**Completion (2026-10-12):** derive and completion selection report unclaimed
positive bond assertions as AromaticBondConstraintFailure under the existing
aromatic-valence failure policy. Selection also checks positive atom assertions
without completion carriers. Stored systems are assessed independently of the
negative-bond exclusions, preserving AromaticBondConstraintMismatch for a
conflicting stored system.

Hückel filters rings, unions, and merged candidates. HMO removes excluded bonds
from component construction and recovers admissible ring regions when a fused
component still contains both endpoints. Clar filters sextets and separates the
selected regions at excluded bonds. All three reject any combined region that
still contains an excluded bond, preserving independent candidates. HMO builds
the component graph directly over eligible atoms without copying the molecule.
The public rule methods and dispatcher enforce the same restriction.

All 2,035 graph unit tests, six resolution properties, and 683 resolution
conformance cases pass. The coronene regression failed under all three rules
before the final rejection check. The property covers both the original open
skeletons and complete carbon renderings with the same topology and bond
assertions. Scoped Clippy with proptest and conformance enabled and warnings
denied, formatting of changed files, and diff checks pass. Workspace and Python
checks remain at the related-change and final gates.

Seven existing default-rule workloads were measured against the saved executable
from commit 0df8e928b. Two paired runs used 40 samples, one-second warm-up, and
two-second measurement, without concurrent builds or tests. The second run
reversed both executable and case order. Columns report observed changes in
mean time:

| Case | Before then after | After then before |
| --- | ---: | ---: |
| ingest benzene | +1.64% | +2.47% |
| ingest naphthalene | +1.32% | +1.53% |
| ingest quinoline | +0.48% | +0.59% |
| ingest purine | +2.41% | +1.18% |
| resolve chain_64 | -1.06% | -0.05% |
| resolve aromatic_fused | +1.17% | -0.07% |
| resolve aromatic_lone_pair | +0.81% | +1.93% |

These runs show a small added cost on the aromatic workloads. The closing
discharge pass still runs until S3m; its removal is not included in these
measurements. Executables, runner, raw samples, and logs are in scratch/236-s3l.

#### S3h1 — Raising writes #a! on a single bond between aromatic atoms **Done**

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

**Completion (2026-10-12):** molecule raising emits the negative aromatic bond
assertion for a Single table bond whose two atoms are explicitly marked
aromatic. Ordinary and explicitly Shared covalent bonds use the same path.
Unknown endpoint flags do not imply aromaticity, and invalid endpoints still
return the existing InvalidReference error. The standalone bond conversion has
no endpoint context and keeps its existing behavior.

Exact cases cover the biphenyl bridge and all its ring bonds, both orders of
the aromatic/aliphatic neighbor, ethane, unknown flags, Shared bonds, and an
invalid endpoint. The new ingestion rejection and three cases requiring `#a!`
failed before the implementation change. All 4,157 IO unit tests, 241
ingestion tests, 683 resolution conformance cases, and 26 SMILES property tests
pass, including the existing exact resolved biphenyl fixture. Scoped Clippy
with warnings denied, formatting of the two changed Rust files, and diff checks
pass. Workspace tests and full-workspace formatting remain at their planned
batch and final gates.

The existing smiles_roundtrip/raise cases were compared with commit 8cceaccdd
using saved executables, 40 samples, one-second warm-up, and two-second
measurement. Two passes used opposite execution orders, without concurrent
builds or tests. Changes below compare medians of per-iteration sample times:

| Case | Before then after | After then before |
| --- | ---: | ---: |
| chain_64 | +0.34% | -0.33% |
| branched | +1.62% | -0.55% |
| aromatic_fused | -0.14% | +0.77% |
| aromatic_lone_pair | -0.87% | -0.16% |

Criterion's slope estimates ranged from -1.18% to +2.24% across these passes;
outliers made its sample-mean comparisons noisier. The branched slowdown did
not repeat, and these short runs do not establish a sub-percent cost. These
inputs check overhead without adding a negative assertion; they do not measure
the cost of storing the new bridge assertion. Executables, raw samples,
estimates, and logs are saved in scratch/236-s3h1.

#### S3m — Automatic placement and closing discharge removed **Done**

**Modules:** graph resolve, property and resolution conformance tests;
docs/development/nomenclature.md.
**Breaking, red→green within the subitem.** [dep: S3h, S3h1, S3i, S3j, S3k, S3l,
S3a]

Remove plan_discharge and its call sites, DischargeContradiction,
ResolveContradiction::Discharge, and ResolveError::DischargeEvaluation and
Discharge. Also remove plan_placement and its four call sites,
PlacementContradiction, ResolveContradiction::Placement, and
ResolveError::Placement. Molecule::inline_constraints remains an explicit
operation. Rewrite the nomenclature guide's Discharge entry: discharge is part
of each resolver's plan and removes the assertions that its own results realize.
Adjust the Concrete entry's reference to the closing pass.

**Tests:** composite resolution keeps assertions on entities it neither adds nor
completes and keeps the molecule-level list unchanged, including bare entity
leaves, singleton wrappers, logical subtrees, and contradictory assertions;
tests of the removed variants go. List preservation is checked under Error
policies; explicit Remove policies retain their existing compaction behavior.
Properties over the skeleton generator: for invariants 5 and 6, the
resolved assertions in standalone and composite resolution are the input's
without the discharged ones, and the molecule-level list is unchanged; for
invariant 7, the renderings that resolve under the Error policies carry no
aromatic or stereo assertion; for invariant 8, resolving a projected molecule
leaves no projected assertion and reproduces its systems under the same model.

**Verification:** the resolve bench groups, against a baseline taken before the
change.

**Completion (2026-10-12):** all four composite entry points omit automatic
placement and closing discharge. They preserve the molecule-level list exactly
under Error policies, including bare entity leaves, singleton logical wrappers,
and contradictions. The final integrity and concreteness checks remain.
Molecule::inline_constraints and constraint validation remain explicit operations.
A separately invoked general constraint-discharge transformer is tracked with
the hydrogen transformers in
[166](166-molecule-ops-2026-07-27.md#explicit-constraint-discharge).

The conformance inputs now rely on the resolver's closed aromaticity reading:
369 fixtures no longer request the blanket `:aromatic-valence :not-aromatic`
default, and 15 omit redundant explicit atom `#a!` marks; one fixture needed
both changes. All 683 cases pass against their original, unchanged snapshots.
Without positive atom, incident-bond, or stored-system evidence, an atom reads
NotAromatic and is excluded from resolution's aromaticity perception. S3h and
S3h1 already cover the boundary cases that require explicit exclusions.
The DSL default retains its meaning, including its use for registry patterns;
deliberately supplied negative assertions retain their meaning and preservation
tests. No general negative-assertion discharge is added.

The related-change gate passed 37,045 workspace tests with all features enabled,
including 177 graph properties, with 10 ignored; workspace Clippy with all
targets/features and warnings denied; and 2,200 Python tests with 2 skipped
after rebuilding against the repository's Python 3.13 environment. That gate
preceded the fixture correction; the corrected corpus then passed its 683
conformance cases without snapshot updates. Formatting and diff checks pass.

The 33 existing resolve workloads were compared with the saved executable from
commit 63e8152ca, using 40 samples, one-second warm-up, and two-second measurement.
Initial apparent regressions prompted paired repeats of six affected cases and
two controls, with executable order reversed in the second pair. These are
changes in mean time, with negative values meaning faster:

| Case | Before then after | After then before |
| --- | ---: | ---: |
| ingest methane | -2.64% | -6.32% |
| ingest benzene | +0.15% | -1.91% |
| ingest pyridine | -1.02% | -2.23% |
| ingest naphthalene | -1.25% | -0.22% |
| counts chain_64 | -8.61% | +1.00% |
| counts aromatic_fused | -12.06% | -7.48% |
| typing chain_64 | -1.55% | -2.58% |
| dense88 unique ownership | -2.18% | +7.25% |

Three further dense88 pairs, using 80 samples and four-second measurement,
gave +0.21%, +2.40%, and -3.29%. The larger initial slowdowns were not sustained;
the variation does not establish a uniform speedup or a precise small-effect
bound. These benchmarks do not consume the migrated conformance inputs, and
the measured production code is unchanged by that correction. Executables,
samples, and verification logs are retained in scratch/236-s3m.

### S4 — Resolvers and transformers **Done**

#### S4a — Resolution and transform benchmark cases **Done**

**Modules:** umol-graph benches/resolve.rs, a new benches/transform.rs target,
and umol-graph Cargo.toml. **Additive, green.** [dep: none]

Add to the resolve bench a Hückel selection group over fused aza-aromatics
whose flexible atoms give option products from 128 through 4096, and a group
for AromaticityResolver::resolve over aromatic cases from the existing corpus.
Add a transform bench target with the kekulizer over aromatic inputs and the
aromatizer over Kekulé inputs, through both consuming and in-place operations.
Input construction and cloning stay outside timing.

The selection scope was narrowed on 2026-10-12. Clar requires benzenoid inputs,
so its rejection of aza-aromatics is not a useful selection baseline. HMO work
is excluded from this subitem. The current selector returns Underdetermined for
a component whose option product exceeds 4096 before searching; timing that
exit would not measure selection work. Successful cases above that product
belong with S4i's changed search bound. These deliberately flexible inputs are
search stress cases, not a representative distribution of molecular inputs.

**Tests:** each case's setup asserts that its input resolves or transforms.

**Completed (2026-10-12).** The selection group times Resolver::resolve with
atom typing, MostSaturated, and natural isotopes; admission checks outside timing
verify the flexible-atom count and option product. The standalone group uses
benzene, pyridine, pyrrole, furan, thiophene, naphthalene, quinoline, purine,
biphenyl, and anthracene. Its setup projects aromaticity out of the resolved
input and checks that standalone resolution restores the same framed molecule.

The transform target uses the existing paired aromatic/Kekulé fixtures for
benzene, pyridine, pyrrole, furan, thiophene, azulene, and indole with a prescribed
donor. Setup checks exact kekulized outputs, perceived system members and forms,
and agreement between consuming and in-place operations. No rejected input or
no-op transformation is timed.

The 42 new workloads have saved Criterion baselines against production commit
810941e0a8ac2cc7ffc3b131332d5618e65854e4, using 40 samples, one-second warm-up,
and a two-second measurement target. Criterion extended the 4096-option case to
about 3.4 seconds to collect 40 samples. Mean selection times were:

| Fixture | Flexible atoms | Option product | Mean time |
| --- | ---: | ---: | ---: |
| Heptaazanaphthalene | 7 | 128 | 0.411 ms |
| Octaazanaphthalene | 8 | 256 | 0.822 ms |
| Decaazaanthracene | 10 | 1024 | 7.461 ms |
| Dodecaazatetracene | 12 | 4096 | 85.277 ms |

Standalone aromaticity resolution ranges from 7.5 to 28.5 µs. Kekulization
ranges from 5.2 to 9.3 µs consuming and 6.1 to 10.6 µs in place; aromatization
from 7.0 to 15.6 µs consuming and 7.3 to 16.2 µs in place. These are starting
measurements, not before/after comparisons. Executables, source copies, raw
samples, confidence intervals, commands, and verification logs are retained in
scratch/236-s4a.

An untimed input probe also found that SMILES-resolved pyrrole, furan, thiophene,
and purine fail the baseline Kekulizer with PostLocalizationValenceInvariant,
carrying OrbitalCount with orbital_count 6 and electron_count 8. The transform
bench inherited invalid aromatic donor inputs from the existing paired fixtures:
they already stored the donor lone pair locally and asserted a one-electron
contribution despite their two-electron system entries. Their passing transforms
do not qualify those donor timings as correctness baselines. The failures
reproduce on 2026-09-29 commit 6de803627e7621920e0af743b54c8202248d9e59;
Kekulizer is unchanged between that commit and the S4a production baseline.
S4a1 corrects the inputs and implementation before further transformer work.

All 110 resolve benchmark checks and 28 transform benchmark checks pass, as do
Clippy for both targets with warnings denied, touched-file formatting, and the
diff check. No workspace gate was rerun for this benchmark-only subitem.

#### S4a1 — Kekulizer donor localization and fixed exposures **Done**

**Modules:** graph transform/kekulizer, tests/kekulization and its fixtures,
benches/transform.rs. **Breaking, red→green within the subitem.** [dep: S4a]

Restore one local lone pair for every stored two-electron aromatic contributor
when removing its system. Preserve the local lone-pair count for zero- and
one-electron contributors, and preserve atom charge, hydrogens, and spin for
fixed exposures. Retain the separate existing localization of a mobile system
charge; an anionic mobile exposure receives exactly one pair.

Remove the first-release cap on the number of fixed zero-/two-electron sites.
Exclude all such sites from double-bond matching and require a perfect matching
on the remaining one-electron sites. This changes the fixed-site boundary of
147; multiple mobile exposures, mixed fixed/mobile demand, and open-shell system
support are unchanged. A fixed assignment without a matching still fails with
NoMatching; this subitem does not change aromaticity selection or hydrogen
placement to make it localizable.

**Public contract:** Kekulizer and KekulizeConfig construction and the Transformer
method signatures are unchanged. Remove KekulizeError::MultiplePrescribedExposures.
Restoration requires a literal local lone-pair count and reuses
UndeterminedExposedAtomLonePairs on failure. transform consumes its input;
transform_into preserves its input on failure; transform_iter omits rejected
candidates. No new public constructors, carriers, or Python methods are added.
Post-localization validation and transaction removal remain S4d's work.

**Tests:** correct the aromatic donor fixtures and their contribution assertions;
validate inputs as well as outputs. Cover N, O, S, and locally charged C donors,
zero-electron acceptors, multiple fixed sites, and existing mobile charge cases.
Run the aromatic resolution corpus under atom typing and counts, asserting
successful localization or the exact NoMatching result for F420's selected
assignment. Check consuming/in-place agreement and unchanged input on failure.
Qualify the corrected benchmark inputs and capture fresh kekulization timings;
do not compare valid donor inputs against the invalid S4a donor measurements.

**Verification:** focused Kekulizer and integration tests, aromatic resolution
conformance, graph crate Clippy, transform benchmark setup checks and timings,
touched-file formatting, and diff review. No workspace gate for this subitem.

**Completed (2026-10-12).** The plan retains each exposed atom's zero-/two-electron
state instead of losing it after matching. The existing fixed-site vector moves
into the system plan; no second id vector is allocated for those sites. Fixed
donors and mobile anionic exposures use one pair-restoration path. The
MultiplePrescribedExposures variant and its rejection are removed.

All 34 aromatic conformance fixtures were checked under atom typing and counts
with MostSaturated. Each model now localizes 33 successfully, with valid inputs
and outputs; this includes the six previously rejected, localizable multiple-site
cases. F420's selected assignment returns NoMatching, preserving the input for
the in-place operation. The N/O/S donor fixtures and borepin's zero-electron
assertion are corrected, and benchmark setup now validates the paired inputs.

F420's failure concerns its selected hydrogen/contribution assignment. Fixture
atom N20 contributes one electron but has only the zero-electron carbonyl
carbons C18 and C21 as neighbors. Removing fixed sites leaves N20 isolated and
an odd nine-atom component; matching covers eight of ten required atoms.
Specifying N17 H0 and N20 H1 before resolution swaps their donor roles: N17
contributes one, N20 two. That assignment has a perfect residual matching and
localizes successfully under both models. Both assignments have 14 aromatic
electrons; the present Hückel ring/union count check does not establish
localizability. S4a1 leaves this selection behavior unchanged.

Verification passes: 41 Kekulizer unit cases, 98 integration cases (including
the corpus and both explicit-hydrogen F420 cases), the existing public Kekulizer
integrity property, all 34 aromatic conformance cases with unchanged snapshots,
graph crate Clippy over all targets with conformance/proptest and warnings
denied, and all 28 transform benchmark setup checks. Diagnostic sources,
before/after logs, matching evidence, executables, and benchmark samples are
retained in scratch/236-kekulizer-diagnosis.

Fresh valid-donor baseline means are 5.48–8.53 µs consuming and 6.38–10.52 µs
in place. The unchanged benzene, pyridine, and azulene inputs were compared
against the saved S4a executable in both run orders (40 samples, one-second
warm-up, two-second measurement). Changes in mean time were:

| Case | Before then after | After then before |
| --- | ---: | ---: |
| benzene consuming | +40.32% | +4.09% |
| benzene in place | +17.06% | +0.62% |
| pyridine consuming | +7.69% | +1.27% |
| pyridine in place | +0.55% | +3.32% |
| azulene consuming | +3.14% | -0.15% |
| azulene in place | +0.15% | -2.55% |

The large initial benzene difference did not persist. These variable runs do
not establish a precise small-effect bound. Donor measurements are new
baselines, not comparisons with the invalid former inputs. No workspace test
or formatting gate was run; touched files are formatted and the diff checked.

#### S4b — Shared system addition **Done**

**Modules:** graph aromaticity, resolve/aromaticity, resolve,
transform/aromatizer, and tests/kekulization. **Breaking, red→green within the
subitem.** [dep: S3i, S4a1]

Move plan_system to `AromaticityPerceiver::plan_aromatic_system(&self, molecule:
&Molecule, atoms: Vec<AtomId>, form: AromaticSystemForm) -> Edits`, keeping the
discharge of S3i. The standalone resolver, the composite resolver's constitution
plan, and the aromatizer's plan use it; the aromatizer's inline copy goes, so it
writes no aromatic constraint on system bonds. Remove
AromaticityPerceiver::add_systems; tests/kekulization builds its fixture from
the DSL.

**Tests:** aromatizer results without aromatic constraints on system bonds;
existing resolution results.

**Completed (2026-10-12).** The crate-private
AromaticityPerceiver::plan_aromatic_system now owns system insertion and member
assertion discharge. The standalone resolver, composite constitution plan, and
aromatizer all use it. The public add_systems method is removed; no new public
constructor, error, or Python interface is introduced. Resolver bond-assertion
discharge remains unchanged. The C60 kekulization input is a DSL fixture that
preserves the previous MOL-derived input, including its bond assertions.

Verification passed: 50 shared-aromaticity unit cases, 615 resolver unit cases,
14 aromatizer cases, 98 kekulization integration cases, 34 aromatic resolution
conformance snapshots, and graph-crate all-target Clippy with conformance and
proptest enabled. Touched Rust files are formatted and the full diff checked.
No workspace gate or new timing comparison was run.

**Outstanding correctness finding.** Before this change, at commit cb28e47b8,
Aromatizer::transform on the valid benzene Kekulé fixture returned a molecule
that failed ValenceInvariantsValidator with OrbitalCount at atom 0:
orbital_count 10, electron_count 9. It adds the aromatic electron contributions
but retains the Kekulé double bonds. S4b leaves that accounting unchanged;
its output tests verify the insertion and discharge contract, not electron
conservation. The benchmark setup likewise checks the perceived systems but
does not validate the complete aromatized output. S4b1 corrects this before S4d;
S4k retains its separate evidence-reading scope. The baseline probe and input
export are retained in scratch/236-s4b.

#### S4b1 — Aromatizer contribution recovery and electron transfers **Done**

**Modules:** graph transform/aromatizer and its unit tests, aromaticity where
contribution context requires it, benches/transform.rs.
**Breaking, red→green within the subitem.** [dep: S4b]

Correct contribution recovery from the localized representation. Distinguish
double bonds that supply the candidate system from exocyclic double bonds;
recover zero-electron carbonyl sites and two-electron anionic carbon donors.
Use the atom's charge and available lone pairs when identifying donors instead
of assigning two electrons from the element alone.

For each selected system, plan the corresponding electron transfers: change
its contributing localized double bonds to single bonds and withdraw one local
lone pair from each two-electron donor. Zero-electron sites withdraw no pair;
external double bonds remain localized. Preserve charge, hydrogens, spin, and
unaffected atom and bond state. These transfers belong to Aromatizer's plan,
not the shared system-insertion helper: resolver inputs already account for
their aromatic electrons differently.

Establish the planned atoms' electron/valence validity before applying edits.
Keep the Transformer entry points, ownership behavior, and existing-system
no-op contract. Move electrons_from_kekule's agreed crate-private visibility
forward from S4l and reuse the perception's ring context for contribution
recovery. Add `AromatizeError::ContradictedAssertion { aromatic_system:
AromaticSystemId, atom: AtomId, constraint: AtomConstraintForm }` and
`AromatizeError::ValenceInvariant { aromatic_system: AromaticSystemId,
mismatch: ValenceMismatch }`, matching S4d's vocabulary. The system id is the
one the planned insertion receives. Missing information needed for a transfer
uses the existing Underdetermined variant; incomplete unrelated state does not
require whole-molecule groundness. Every failure precedes mutation. S4e still
owns transaction removal, and S4k owns evidence handling.
HMO/Clar changes and resolver selection/backtracking remain outside this fix.

**Investigation (2026-10-12).** At commit 7173a14f1, the aromatic conformance
corpus gives 34 testable Kekulé inputs per valence model, using atom typing and
counts with MostSaturated and Daylight Hückel. This includes F420 with prescribed
hydrogens; the original F420 assignment retains S4a1's NoMatching result.
Aromatizer produces 28 electron-invalid outputs and six unchanged inputs per
model. Diagnostic transfers restore conservation, but contribution recovery
still misses six systems and only recovers F420's benzene ring. Supplying the
reference contributions and membership plus those transfers reproduces all
34 references; that control uses information unavailable to a real reader and
is not an implementation. The manifest, probe, results, and faithful input/output
exports are retained in scratch/236-aromatizer-diagnosis.

**Tests:** use inline parameterized success, identity, and error cases for
one-electron bond transfer, two-electron N/O/S and anionic-C pair transfer,
zero-electron sites retaining external double bonds, unchanged ineligible or
already-aromatic input, required versus unrelated undetermined data, assertion
conflicts, and electron-count mismatches. Check exact outputs, ownership-path
agreement, and unchanged input even when a later candidate fails. Do not add
external fixture files or a permanent suite copied from the conformance corpus.

The conformance-derived inputs remain a diagnostic sweep in scratch/236-s4b1,
under both valence models, including F420 with prescribed hydrogens. Validate
complete inputs and outputs. Distinguish conservation from exact system recovery:
porphine-annulene permits a larger electron-valid perceived system after
localization, so a universal exact-membership roundtrip is not the contract.

**Verification:** focused transform tests and the diagnostic sweep, aromatic resolution
conformance, graph crate Clippy, touched-file formatting, and diff review.
Strengthen transform benchmark setup to validate complete aromatized outputs
and capture fresh consuming/in-place timings. Earlier timings of electron-invalid
outputs are not equivalent-correctness baselines. No workspace gate for this
subitem.

**Completed (2026-10-12).** Aromatizer reuses one candidate ring set and reads
localized contributions once. Its plan transfers internal double bonds and
donor lone pairs, retains external carbonyl bonds, and checks the planned atoms
before mutation. The approved ContradictedAssertion and ValenceInvariant errors
identify the affected planned system and atom; unrelated underdetermined state
is preserved. electrons_from_kekule is crate-private. Existing-system no-op and
Transformer ownership contracts remain unchanged.

Permanent coverage consists of three inline behavior tables: nine successful
transformations, five identities, and five failures, plus the existing
configuration test. Success cases include N/O/S and anionic-C donor variants,
zero-electron sites, assertion discharge with unrelated-assertion preservation,
and unrelated partial state. The later-system failure checks that the input is
unchanged. No external fixture files or corpus integration tests were added.

Verification passed: 70 transform unit cases, 98 existing kekulization cases,
34 aromatic resolution conformance snapshots, the Aromatizer integrity property,
and graph-crate all-target Clippy with conformance and proptest enabled. The
scratch diagnostic passed all 68 runs across the 34 inputs and two valence models,
validating complete inputs and outputs and comparing ownership paths. All 28
transform benchmark setups pass with complete output validation. The final
20 Aromatizer unit cases, touched-file formatting, and diff checks pass. No
workspace gate was run.

Fresh timings cover seven fixtures through consuming and in-place Aromatizer
paths. Final means span 7.45–16.84 µs. A single sequential comparison against the
first corrected implementation measured 4.6–13.4% lower means after omitting
redundant zero-default constraints from temporary validation forms, keeping the
common form within its inline capacity. This is an exploratory optimization
comparison, not a precise regression estimate against the earlier invalid-output
implementation. Commands, manifests, diagnostic results, binaries, and raw timing
samples are retained in scratch/236-s4b1.

#### S4d — Kekulizer validity during planning **Done**

**Modules:** graph transform/kekulizer and the kekulization tests. **Breaking,
red→green within the subitem.** [dep: S3a, S4a1, S4b1]

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

**Implementation, 2026-10-13.** Kekulizer now checks the planned localized
valence and retained multicenter/dative assertions, then electron conservation,
before applying any edits. The two planned error variants replace the four
PostLocalization variants. Both ownership paths use the same plan; transform_into
moves the receiver's storage through the consuming application path after planning
succeeds, without a transaction or recovery copy. Unchanged spin data and atoms
outside the systems are not subjected to whole-molecule validation.

The temporary validation form copies the intrinsic fields and constructs only
the incidence constraints the conservation check uses. Absent aromatic incidence
and zero multicenter/dative incidence use the bare validator's zero defaults.
When both a relation and its own assertion are absent, no temporary reading or
meet is needed. Exposed atoms are consumed once in participant order, avoiding
a repeated search for each system atom. Unknown incidence remains underdetermined,
rather than becoming zero; a bottom reading fails planning. These shortcuts do
not rewrite stored assertions or add a whole-molecule groundness requirement.

Inline unit cases cover incompatible valence and donation assertions, a later
system's mismatch after an earlier partial atom, preserved spin and overlay
data, unknown external bond order, and normalized or bottom external order.
The existing cyclopentadienyl error cases now use inline inputs and check both
ownership paths and iterator omission. No external fixture files or corpus
tests were added. The new property constructs complete single and coupled
C/N/O/S rings with fixed contributions, then varies a charge and an optional
valence assertion. It checks output conservation, retained assertions, exact
ownership-path agreement, and unchanged input on failure, without skipped cases.

Focused verification passed: 76 transform unit cases, 98 existing kekulization
cases, both Kekulizer properties at 512 cases each, and graph-crate all-target
Clippy with conformance and proptest enabled. Touched-file formatting and diff
checks passed. No workspace gate was run.

The final implementation was compared with the clean pre-S4d binary at
7432a0fb9b11dc28dec8cd22279e2f1c39e9c349. Three paired passes used the same seven
fixtures and both ownership paths, with 40 samples, a one-second warmup, and a
two-second measurement per case. The second and third passes reversed the first
pass's before/after order. Each entry below is the median of the three paired
percentage changes in the sample means; negative means faster.

| Fixture | Consuming | In place |
| --- | ---: | ---: |
| Benzene | +1.1% | −5.0% |
| Pyridine | +1.7% | −7.7% |
| Pyrrole | +5.0% | −5.4% |
| Furan | +2.5% | −7.9% |
| Thiophene | +1.4% | −7.2% |
| Azulene | +3.7% | −11.0% |
| Indole with prescribed donor | +3.5% | −11.4% |

There is material timing variation. The second pass's old pyrrole measurement
and new furan measurement had large timing outliers, so pooled means would imply
a misleading pyrrole improvement and exaggerate the furan regression. All runs,
sample distributions, paired changes, source copies, binaries, and commands are
retained in scratch/236-s4d; no run was discarded. The measurements support a
remaining consuming-path cost, not a one-percent bound or a claim that the cost
is unavoidable. The in-place improvement is consistent across all three passes.

**Planned-order table comparison, 2026-10-13.** Two direct-valence trials remove
the whole-molecule Vec<Option<i64>>. Internal bonds contribute one each, with
one added for a matched atom; external bond orders retain the existing reading.
The trials differ only in how they identify internal neighbors: the existing
aromatic-system incidence index, or contains on the plan's existing atom list.
The temporary atom construction, assertion meets, conservation check, matching,
and edit application are unchanged. Both trials pass the same focused tests,
512-case properties, and Clippy checks.

Each trial was compared with both the table implementation and the pre-S4d
binary in three passes, rotating execution order. The seven fixtures, both
ownership paths, 40 samples, one-second warmup, and two-second measurement are
unchanged. These ranges span the fixtures' median paired changes in sample
means relative to the table implementation; positive means slower.

| Direct-valence membership | Consuming | In place |
| --- | ---: | ---: |
| Aromatic-system incidence index | +1.0% to +4.7% | +1.1% to +3.0% |
| Existing plan atom list | +1.2% to +3.3% | −0.6% to +1.6% |

The table implementation is retained: neither replacement recovers the
consuming regression. The atom-list trial's apparent in-place pyridine gain
varies in sign across passes and also changes sign when using sample medians.
Both replacements introduce membership work, so these timings do not isolate
the allocation's cost or establish that the table is optimal. The scope remains
the seven small aromatic fixtures. Source snapshots, exact diffs, binary hashes,
commands, full timing ranges, and raw samples are retained in scratch/236-s4d
under the orders and orders-list comparisons. The production source exactly
matches the tested table version.

**Accepted and completed (2026-10-13).** The measured trade-off is accepted.
The consuming regression is relative to its own pre-S4d baseline; the retained
consuming implementation remains about 3–6% faster than the in-place version.
The latest table measurements in the orders-list comparison give the following
median run means across three passes:

| Fixture | Consuming (µs) | In place (µs) |
| --- | ---: | ---: |
| Benzene | 6.467 | 6.860 |
| Pyridine | 6.552 | 6.918 |
| Pyrrole | 5.458 | 5.706 |
| Furan | 5.493 | 5.692 |
| Thiophene | 5.436 | 5.723 |
| Azulene | 9.386 | 9.652 |
| Indole with prescribed donor | 8.232 | 8.474 |

Both paths share planning and consuming edit application. The in-place path
also constructs and discards the empty replacement used by mem::take; its
isolated cost was not measured, and the Criterion batching loops differ.

#### S4e — In-place transform_into without a transaction **Done**

**Modules:** graph transform/aromatizer and transform/delocalize_charge.
**Green.** [dep: S4b1]

Aromatizer and DelocalizeCharge transform_into apply their plan in place instead
of through transact. Planning failures leave the receiver unchanged.

**Tests:** existing transform_into results, and an unchanged receiver on
Aromatizer planning failure.

**Completed (2026-10-13).** Both transform_into implementations derive their
complete plan before moving the receiver through Molecule::apply. Application
uses the consuming editor without an undo journal or recovery copy, retaining
the final integrity check. Empty plans leave the receiver untouched. Planning
errors leave Aromatizer's receiver unchanged; DelocalizeCharge remains
infallible. The planners, public signatures, and chemistry behavior are unchanged.

Existing inline success, identity, and error cases already cover both ownership
paths, including an Aromatizer failure in a later system. Reused that coverage
without adding tests or fixtures. Verification passed: all 76 transform unit
cases, both affected integrity properties at 512 cases each, graph-crate
all-target Clippy with conformance and proptest enabled, touched-file formatting,
and diff checks. No workspace gate or new timing comparison was run. The contract
and verification logs are retained in scratch/236-s4e.

#### S4g — Mismatch policies removed **Done**

**Completed with S3d and S3e.** The aromaticity policy removal is included in
S3d, and the stereo policy removal is included in S3e.

**Modules:** graph resolve/stereo and resolve; umol-py
resolve bindings, the Python package exports, and Python tests. **Breaking,
red→green within the subitem.** [dep: none]

Remove StereoMismatchPolicy, with the fields tetrahedral_stereo_mismatch and
cis_trans_stereo_mismatch. Every stereo mismatch is a contradiction: the
resolver's Keep, RemoveConstraint, ReplaceEntity, and RemoveBoth branches go.
Remove the Python classes, keywords, getters, and repr fields.

**Tests:** each mismatch is a contradiction in standalone and composite
resolution; tests of the removed branches go.

#### S4h — Rule-supplied components and candidate systems **Done**

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

**Completed (2026-10-13).** The perceiver dispatches component discovery,
containing candidates, and contribution acceptance to the three rules.
Hückel's acceptance uses the actual alternatives, retaining the greatest
reachable total for each residue modulo four in fixed stack storage. Zero is
an aromatic contribution; NotAromatic is not. HMO and Clar retain their
unconditional contribution bound. The old claim_candidates and accepts_range
methods are removed; the new machinery is crate-private.

The selection includes flexible non-ring atoms in HMO components. It reuses
candidate topology while eligibility is unchanged and recomputes it when a
selected NotAromatic alternative changes eligibility. Candidates are
deduplicated and visited by increasing size. Each acceptance check covers all
required member atoms together, preserving the per-atom pruning criterion
without repeating that check for every member. Keep still disables pruning;
the enumeration bound and public algorithm selection remain S4i's scope.

Verification passed: 282 focused aromaticity cases, including all three rules'
candidates and a flexible non-ring HMO atom in both option orders; an exhaustive
comparison of 33,825 small contribution-choice domains with direct enumeration;
the existing selection-equivalence and aromatic-constraint properties at 512
cases each; 34 aromatic resolution snapshots; graph-crate all-target Clippy with
conformance and proptest enabled; touched-file formatting and diff review.
No external test fixtures, workspace gate, or Python changes were added.

The initial implementation repeated candidate acceptance for different required
atoms and measured 8.5–26.7% slower than the baseline. Removing that repetition
and trying smaller candidates first recovered the cost while retaining the
exact contribution predicate. The final consuming-resolution selection timings
against pre-S4h commit 28173f19ceaa027138a809f3276f3878ac0741bb are:

| Fixture | Before (ms) | After (ms) | Change |
| --- | ---: | ---: | ---: |
| Heptaazanaphthalene, 128 combinations | 0.401 | 0.280 | −30.3% |
| Octaazanaphthalene, 256 combinations | 0.803 | 0.556 | −30.9% |
| Decaazaanthracene, 1024 combinations | 7.255 | 5.525 | −23.8% |
| Dodecaazatetracene, 4096 combinations | 82.869 | 77.023 | −6.9% |

Absolute times are median run means from three paired passes; changes are the
medians of the paired percentage changes. Each case used 40 samples, a one-second
warmup, and a two-second requested measurement, extended by Criterion when needed.
The middle pass reversed before/after order. Every pass improved on every fixture;
no run was discarded. These are Hückel selection stress cases, not measurements
of HMO, Clar, or resolution workloads generally. Source snapshots, binary hashes,
commands, full timing ranges, and raw samples are retained in scratch/236-s4h.

#### S4i — AtomSelectionAlgorithm and the search bound **Done**

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

**Completed (2026-10-13).** AtomSelectionAlgorithm is available in Rust and
Python, including the config keyword, getter, and representation;
BacktrackingDfs remains the default. Validity checks apply independently of
algorithm choice. The bound allows 4096 search perception calls per component;
an attempted additional call returns Underdetermined. The option-product
precheck and its temporary allocation are removed.

Public resolver tests cover a full search at the bound and a 14-membered ring
with thirteen flexible nitrogens: its 8192 assignments exceed FullDfs's bound,
but BacktrackingDfs resolves it. Under Keep, both algorithms reach the bound.
The public-resolution property replaces the copied test-side selector and
compares outcomes and reports over small skeletons with undetermined hydrogen
counts, both valence sources, both tie-breaks, both failure policies, and the
Hückel, HMO, and Clar rules, without filtering outcomes.

The flexible-nitrogen indole case gives the same HmoInvalidInput contradiction
under both algorithms at thresholds 0.2 and 0.1. One admitted completion has
nine pi electrons; HMO rejects odd electron counts and selection propagates
that contradiction. This existing behavior is preserved and tested exactly;
these cases do not demonstrate successful HMO resolution.

**Verification:** 622 resolver unit cases, 4096 generated search-equivalence
cases (seed 20261013), 34 aromatic resolution snapshots, 52 Rust binding cases,
and 220 Python cases passed. Graph and Python all-target Clippy passed with
property/conformance features and warnings denied. Touched Rust files were
formatted and the full diff reviewed. Workspace tests and the full formatting
pass remain at their planned gates.

**Performance:** baseline commit 278432e482988543ca3dcc69c5e9cdd13ffda938 was
built and measured before edits. Three passes rotate the baseline,
BacktrackingDfs, and FullDfs order, with 40 samples, one-second warmup, and
two-second requested measurement time. No builds or tests ran concurrently.
Times are medians of the three run means; changes are paired percentages.

| Case | Before (ms) | BacktrackingDfs (ms) | FullDfs (ms) | Median default change | Paired range |
| --- | ---: | ---: | ---: | ---: | ---: |
| heptaazanaphthalene, 128 assignments | 0.283 | 0.280 | 0.368 | −1.04% | −2.39% to −0.68% |
| octaazanaphthalene, 256 assignments | 0.560 | 0.558 | 0.752 | −0.40% | −1.21% to +0.29% |
| decaazaanthracene, 1024 assignments | 5.568 | 5.299 | 7.177 | −5.18% | −5.44% to −4.82% |
| dodecaazatetracene, 4096 assignments | 78.533 | 70.170 | 84.854 | −10.65% | −11.30% to −10.07% |

The additional tridecaazaannulene case, with 8192 assignments, resolves under
BacktrackingDfs in 65.279 ms (64.400–66.691 ms). Its baseline and FullDfs
outcomes are Underdetermined, so those exits are not timed. These Hückel search
stress cases do not establish a workload distribution or HMO/Clar performance.
The measurements do not isolate the cause of the larger cases' improvements.
Commands, binaries, source hashes, and all raw samples are in scratch/236-s4i.

**Deduplication follow-up (2026-10-13).** Direct timing of the completion
membership scan found 36.75 ms, or 52.3% of resolution time, in the largest
fixture. Its 1712 distinct perceived combinations required 1,464,616 failed
equality comparisons. The later restriction scan received one survivor and
performed no equality comparisons. The diagnostic is in scratch/236-s4i-scaling.

Removing early deduplication preserved the checked results and reduced the
largest case by 53.2%, but retaining repeated combinations made the
4096-assignment duplicate fixture 8.0% slower. Insertion-ordered hash sets
retained early compaction, but hashing complete forms increased the smallest
search fixture by 7.9% with the default hasher and 11.4% with XXH3.

The retained implementation uses BTreeMap for full-combination identity,
recording each distinct combination's first discovery position and restoring
that order before validity and selection. With P perceived combinations and K
distinct results, collection uses O(P log(K + 1)) key comparisons and order
restoration uses O(K log K) integer comparisons. It retains K combinations and
introduces no form cloning or dependency. The later restriction scan uses a
HashSet of borrowed references, retaining first-occurrence order; a single
survivor needs no hash set. Search order, the perception-call bound, failure
policies, and tie-break semantics are unchanged.

Three alternating paired passes use the same 40-sample, one-second-warmup,
two-second-measurement protocol. The baseline is the working S4i implementation
before this follow-up. Times below are medians of run means; percentages are
medians of paired changes, not ratios of the displayed medians.

| Case | Before (ms) | Ordered map (ms) | Median paired change | Paired range |
| --- | ---: | ---: | ---: | ---: |
| heptaazanaphthalene, 128 assignments | 0.291 | 0.300 | +1.28% | −1.31% to +2.88% |
| octaazanaphthalene, 256 assignments | 0.572 | 0.563 | −1.55% | −1.82% to −0.52% |
| decaazaanthracene, 1024 assignments | 5.548 | 4.659 | −16.40% | −16.58% to −15.75% |
| dodecaazatetracene, 4096 assignments | 73.422 | 38.804 | −46.62% | −47.45% to −46.53% |

Synthetic selection fixtures make 16, 256, or 4096 assignments produce the same
single ring and restriction; their median changes are +0.48%, −0.62%, and
−0.07%. They exercise duplicate handling, not a chemistry workload distribution.
End-to-end benzene and pyridine ingestion changed by +0.41% and +1.10%
respectively. No builds or tests ran concurrently with measurements, and no
runs or samples were discarded. All four trial implementations, raw samples,
commands, binaries, and source snapshots are retained in scratch/236-s4i-dedup.

**Verification:** 628 resolver cases, including three duplicate-selection
behaviors under both algorithms, passed. Those cases check a shared restriction,
selection between distinct restrictions with its exact report, and retained
ambiguity under Strict. The 4096-case search-equivalence property (seed
20261013), 34 aromatic resolution snapshots, and graph all-target Clippy with
property/conformance features also passed. The incremental diff was reviewed;
workspace tests and the full formatting pass remain at their planned gates.

#### S4j — One aromaticity resolution path **Done**

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

**Completion (2026-10-13).** Standalone resolution now calls the same selection
as composite resolution, with an empty completion carrier. The perceiver owns
contribution reading and inconsistency classification. Stored systems are checked
against the whole perceived partition; their separate perception is removed.
Keep excludes additions throughout a failed stored system's search component,
while independent components can still resolve. Mismatches remain contradictions.
Literal atom assertions can determine contributions when the stored electron
vector is undetermined; genuinely nonliteral complete evidence still gives
Underdetermined. Existing stored systems are preserved.

Selection with no flexible contributions derives once for the whole molecule.
Membership queries binary-search the already sorted member list for one system.
Multiple systems use a sorted SmallVec of atom/system-index pairs, with inline
storage through 32 perceived atoms. This avoids quadratic membership scans;
it does not change perception or search complexity. S4i's ordered combination
and restriction deduplication remain in place.

The overlap regression was run against the pre-S4j production code: Error and
Keep both reached Apply(Integrity(AromaticSystemsOverlap)), while conformance
validation accepted the stored benzene inside naphthalene. The new code gives
AromaticSystemFailure under Error and in validation, and preserves the molecule
under Keep. Exact cases also cover stored systems outside the search components,
aromatic bonds between systems, mismatches before pruning, electron-vector
length failures, and report preservation. The Hückel skeleton property covers
single, fused, and linked five- and six-membered rings, comparing resolution
with and without the stored first ring. It checks exact resulting systems,
Error/Keep, integrity, and agreement of consuming and in-place resolution.

Validation passed: 305 aromaticity unit cases, 643 resolver unit cases,
21 resolve-related property-target checks with PROPTEST_CASES=4096 and seed
20261013, 34 aromatic resolution conformance cases, and graph Clippy with all
targets, proptest/conformance, and warnings denied. Touched Rust files were
formatted with nightly rustfmt; the full workspace format remains at S7a.
No workspace or Python gate ran for this subitem.

**Measurements.** Baseline: committed pre-S4j code at 52d1ef4e. Three before/after
pairs per fixture used alternating order, 40 samples, one second of warmup, and
two seconds of requested measurement, with no concurrent compilation or tests.
The table gives the median paired percentage change and the range of the three
pairs; negative is faster. These are per-fixture results, without workload
weighting. Setup and input cloning are outside standalone resolution timing.

| Standalone aromaticity resolution | Median change | Paired range |
| --- | ---: | ---: |
| benzene | +0.56% | -0.33% to +0.97% |
| pyridine | -0.10% | -1.52% to +0.46% |
| pyrrole | -0.39% | -0.76% to +0.13% |
| furan | -0.49% | -0.52% to +0.33% |
| thiophene | -0.32% | -0.99% to +0.56% |
| naphthalene | -0.03% | -2.01% to +0.68% |
| quinoline | -0.52% | -0.71% to +0.35% |
| purine | +0.36% | +0.22% to +1.16% |
| biphenyl | -0.69% | -1.72% to -0.49% |
| anthracene | +0.51% | +0.10% to +0.53% |

| Additional comparison | Median change | Paired range |
| --- | ---: | ---: |
| SMILES ingest: benzene | -20.40% | -21.09% to -19.76% |
| SMILES ingest: pyridine | +1.81% | +1.72% to +2.38% |
| Hückel selection: 128 choices | +1.41% | +0.25% to +2.18% |
| Hückel selection: 4096 choices | +0.45% | +0.36% to +2.09% |

Standalone median changes are within 0.7% of baseline. The benzene ingest saving
comes with the path that avoids search setup for fixed contributions; pyridine
ingest remains 1.81% slower in this comparison. The 4096-choice case is 0.45%
slower, retaining S4i's large-search improvement. These changes are reported
separately rather than combined into an overall speedup.

Commands, regression logs, source snapshots, frozen binaries, individual
Criterion results, and the comparison script are in scratch/236-s4j/;
final/ contains the accepted timing comparison. The earlier linear-scan and
always-indexed trials are retained as experiment evidence.

#### S4c — Resolver plan application as an invariant **Done**

**Modules:** graph resolve, the five individual resolvers, and their tests.
**Breaking, red→green within the subitem.** [dep: S3a, S4j]

S4j establishes the stored-system overlap check during planning before this
subitem removes application errors. Replace post-plan `apply(edits)?`,
`finish()?`, and `transact([edits])?` in the individual resolvers, and the
composite resolver's phase-plan application, with expect and an invariant
message, as the transformers do. Remove the Apply variants from the eight error
types and the composite ResolveError mappings that wrap them.

**Tests:** existing resolution and projection results, including S4j's overlap
cases under Error and Keep through resolve and resolve_into; remove tests that
construct the removed variants; a property over the skeleton generator, for
invariant 1, replacing the four fixed-SMILES integrity properties: every
determined result of the aromaticity resolver's resolve and project, the
kekulizer, and the aromatizer passes the integrity check, and the consuming and
in-place forms agree.

**Completion (2026-10-13):** the eight individual Apply variants are removed.
Application, probes, finalization, and commits of the resolvers' own plans now
assert the planner's applicability and integrity guarantees. ResolveError's
Isotope, Commit, and Apply variants and ProjectError::Apply are removed too.
BondsError, IsotopeError, MulticenterBondsError, and StereoError remain empty
enums; operation signatures and planning outcomes are unchanged.

The composite borrowed operations retain Transaction::run. Its checked
MoleculeApplyError boundary is asserted separately from the returned planning
Result; only a determined result commits, so planning rejection still rolls
back. Externally applying saved Edits remains checked, and the stale-plan tests
remain. No tests constructed the removed resolver error variants. Rustdoc now
states the trusted-plan guarantee and lists only planning errors. Python has no
matches on the removed variants and needs no binding changes.

The four fixed-SMILES integrity properties are replaced by three property groups
over raw C/N/O/S ring skeletons, stored first-ring skeletons, and independently
concrete saturated-carbon skeletons. They check standalone aromaticity resolution
and projection, Kekulizer, Aromatizer, DelocalizeCharge, and composite resolution:
successful publication passes the public editor integrity gate, consuming and
in-place outcomes agree, and rejection preserves the input. Successful Kekulé
outputs additionally exercise aromatization with localized double bonds. The
separate Kekulizer valence-conservation property is unchanged.

Validation passed: 643 resolver unit cases, the four publication properties,
22 resolver-filtered property-target checks, 159 projection-filtered checks,
and 34 aromatic resolution conformance cases. Property runs used 4,096 cases
and seed 20261013. Graph-crate all-target clippy with proptest and conformance
features, touched-file nightly formatting, and git diff --check passed. No
workspace test run was added. The contract inventory and check logs are in
scratch/236-s4c/.

#### S4f — Individual resolver _into without a transaction **Done**

**Modules:** graph resolve isotope, aromaticity, stereo, bonds, multicenter.
**Green.** [dep: S4c]

The eight resolve_into and project_into methods apply their plan in place
instead of through transact. Outcomes decided during planning leave the receiver
unchanged.

**Tests:** existing _into results, and an unchanged receiver on
underdetermination and contradiction.

**Completion (2026-10-13):** all eight methods finish planning before moving
the receiver through Molecule::apply. They use no undo journal or recovery copy;
empty plans leave the receiver untouched. Nonempty plans retain the publication
integrity check and S4c's asserted applicability guarantee. Planning errors,
underdetermination, and contradiction leave the input exactly unchanged.
Rustdoc states this contract; signatures, planning, and Python bindings are
unchanged. The composite resolver's transactions remain for its later subitem.

Existing aromaticity and multicenter planning cases now also assert receiver
preservation through resolve_into. Validation passed: 643 resolver unit cases,
four publication properties, 140 isotope-filtered property-target checks, and
159 projection-filtered checks. Property runs used 4,096 cases and seed
20261013. Graph-crate all-target clippy with proptest and conformance features,
touched-file nightly formatting, and git diff --check passed. No workspace
test run or timing comparison was added. The contract inventory and check logs
are in scratch/236-s4f/.

#### S4k — Aromatizer reads evidenced **Done**

**Modules:** graph transform/aromatizer and aromaticity. **Green.** [dep: S3d,
S3l, S4b1, S4j, S3a]

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

**Completion (2026-10-13):** the aromatizer reads each atom's evidenced aromatic
valence once, before recovering its localized contribution. Negative evidence
excludes the atom; top evidence skips a redundant meet. Only conflicting atom
ids are retained for comparison with selected systems, and a mismatch reports
the id the system would receive on application. Nonliteral count assertions
remain admissibility constraints: a set containing the recovered contribution
is accepted; a disjoint set is an AromaticValenceMismatch.

The shared derive_from_rings path checks required atoms and bonds without
repeating ring perception. Bottom classification is shared with the complete
reader; the latter's success path and signature are unchanged. Inputs without
rings still undergo the evidence and claim checks, then return without running
the rule or allocating a contribution vector. This preserves acyclic no-ops,
including donors under Clar. All rejections precede mutation. The existing
no-op for inputs already carrying aromatic systems remains unchanged.

The original implementation fails nine new regression cases. Final verification
passes 33 aromatizer cases, 418 aromatic-filtered unit cases (with some overlap),
five publication properties at 4096 cases each and seed 20261013, 34 aromatic
resolution conformance cases, graph all-target Clippy with proptest/conformance
and warnings denied, targeted formatting, and diff checks. The new invariant 11
property uses independently localized carbon skeletons, with and without varied
atom/bond assertions, under Hückel, HMO, and Clar. Exact perception expectations
state their model; the property does not require equal outcomes between rules.
No external fixture files or public API changes were added. The workspace and
Python gates were not rerun for this subitem.

The pre-change transform executable was built from committed
6c317534b1ff2a6ae02544db38dbbad672c27cea. The final comparison uses three
alternating before/after pairs per case, 40 samples, one second of warmup, and
two seconds requested measurement, with no concurrent builds or tests. Inputs
are cloned outside timing. Changes below are medians of paired percentage
changes; parentheses give the full paired range. Positive means slower.

| Fixture | Consuming change | In-place change |
| --- | ---: | ---: |
| benzene | +0.48% (+0.18% to +0.91%) | +2.14% (+1.04% to +4.50%) |
| pyridine | +0.82% (+0.72% to +1.76%) | +0.79% (-0.62% to +1.36%) |
| pyrrole | +1.74% (+0.47% to +2.77%) | +0.87% (+0.60% to +2.20%) |
| furan | +3.04% (+0.80% to +8.89%) | +3.42% (+0.88% to +3.66%) |
| thiophene | +0.90% (-1.21% to +1.45%) | -0.23% (-0.97% to +1.53%) |
| azulene | +0.53% (+0.47% to +0.94%) | +1.23% (-1.83% to +2.73%) |
| indole, prescribed donor | +2.14% (+1.88% to +2.74%) | -0.01% (-1.49% to +9.00%) |

Eight of fourteen medians are within one percent of baseline; the larger ones
remain explicit above. The ranges do not establish a one-percent bound. These
are fixed Hückel fixtures with unasserted aromatic valences, not a workload
distribution or timings for HMO, Clar, or set/expression-heavy evidence.

The final implementation avoids meeting top evidence, retaining full evidence
forms, and repeating perception. The added work still includes evidence reads
and shared claim checks; these timings do not isolate their individual costs.
Source snapshots, frozen executables, exact commands, all raw timings, and
verification logs are in scratch/236-s4k/. Its final/ directory contains the
final comparison; the root-level timing files retain the initial implementation.

#### S4l — Restrict methods with unchecked contextual inputs **Done**

**Modules:** graph aromaticity and its rule modules, resolve/aromaticity;
tests/matching and benches/resolve.rs. **Breaking, red→green within the
subitem.** [dep: S4i, S4j]

Make only these five currently public methods crate-private:

- HueckelAromaticity::find_from_rings, HmoAromaticity::find_from_rings, and
  ClarAromaticity::find_from_rings: their RingSet and molecule must agree.
- AromaticityResolver::select and select_with_report: their ResolveState must
  contain candidates admitted for the supplied molecule and compatible pending
  aromatic systems.

Keep the modules, rule types, errors, HmoOutput, AromaticityPerceiver,
AromaticityDerivation, ResolveState, and AtomCompletions public. Keep the
perceiver's new, find_systems, and derive and all three admit methods public.
ResolveReport and its Python interface retain the same candidate carrier.
Existing crate-private helpers retain their visibility; S4b1 already changed
electrons_from_kekule.

Migrate tests/matching's direct Clar call to AromaticityPerceiver::find_systems
under the Clar model, preserving the fixture and expected result. Adapt the
direct-selection duplicate-combination benchmark to a public Resolver entry
point with candidates produced by admission; verify that admission supplies the
intended alternatives during setup. Preserve its duplicate-combination coverage
and the existing admission and public-operation selection benchmarks. S4i has
already replaced the test-side exhaustive selector with its public-resolution
property; no further property removal is needed.

**Tests:** the existing graph property suite and tests/matching through the
public operations; build the resolve benchmark target.

**Completed (2026-10-13).** Only the five named methods became crate-private;
their Assumes sections state contextual validity prerequisites. The modules,
data types, admission methods, public perception operations, and report
interfaces retain their visibility. The perceiver now calls the rule methods
directly, replacing its redundant dispatch wrapper; each perception still
computes bond exclusions once, with unchanged algorithms and error mapping.

The Clar matching fixture uses the public perceiver and checks its complete
18-member result. The duplicate-combination benchmark now obtains its candidates
through admission and times Resolver::resolve. Setup checks the exact admitted
alternatives, verifies the same perceived system for every one of the 16, 256,
and 4096 assignments, and checks the unresolved report. Its case name is now
resolve_same_partition because composite-resolution timings are not comparable
with the former selection-only timings. The other benchmarks and property laws
are unchanged.

Verification passed: the matching case, 418 aromaticity-filtered unit cases,
179 checks in the graph property target with 256 cases and seed 20261013,
all three migrated Criterion cases in test mode, and graph all-target Clippy
with proptest and conformance enabled. Focused formatting and git diff --check
passed. No timing comparison or workspace/Python gate was run for this subitem.
The contract and verification logs are in scratch/236-s4l/.

#### S4m — Planning rule in the data-type guide **Done**

**Modules:** docs/development/data-types.md, docs/development/nomenclature.md,
and the graph Transformer trait rustdoc. **Documentation.** [dep: S4d, S4e, S4f]

Replace the transformer paragraph, which says transform_into uses a transaction
and restores on error, with the rule: a transformer derives a complete plan from
its input; every rejection is reported during planning, before mutation,
including a plan that would violate integrity; a plan applies to the molecule it
was derived from, so a failed application is a defect; transform consumes and
applies the plan to its input, and transform_into applies it in place without a
transaction or rollback method. State the same rule for the individual
resolvers' resolve_into and project_into. Check that the nomenclature guide's
Transformation entry still agrees.

**Completed (2026-10-13).** The data-type guide, Transformation entry, and
Transformer rustdoc now state complete planning before mutation, rejection
without changing the input, and application/publication failure as a defect.
The guide applies the same rule to individual resolvers and explains why the
composite Resolver retains its transaction across phases.

#### S4n — Aromatic system names in the reworked code **Done**

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

**Completed (2026-10-13).** The public methods, result fields, diagnostic
variants and fields, private methods and locals, test names, and benchmark
callers use the aromatic-system names. The private Kekulizer SystemPlan is now
AromaticSystemPlan. Valence admission also constructs the renamed ResolveState
field. The changes preserve operation behavior, test inputs, property laws, and
diagnostic display text; derived Debug output uses the renamed fields and variants.

**S4 verification (2026-10-13).** All checks passed:

- Graph property target: 179 checks, with 256 generated cases and seed 20261013.
- Graph resolution conformance: all 683 cases, with unchanged snapshots.
- Full workspace tests: 21,962 passed, 9 ignored, no failures.
- Workspace all-target Clippy with graph proptest/conformance enabled and warnings denied.
- Python 3.13 extension rebuilt with maturin; pytest: 2,204 passed, 2 skipped.
- Touched-file nightly formatting, git diff --check, and full diff review.

The contract inventory and logs are in scratch/236-s4mn/. The full workspace
formatting pass remains at S7a. S4 is complete.

### S5 — Copy-on-write storage

#### S5a — Molecule constraints **Done**

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

**Completed (2026-10-13).** Molecule stores its constraints in Arc, and all
mutating access detaches shared storage through constraints_mut. Normalization
takes the inner Constraints without replacing the Arc; reframing and remapping
move uniquely owned constraints with Arc::unwrap_or_clone. Both compaction paths
obtain mutable storage once. MoleculeEntries and the public constraint accessors
retain their owned and borrowed types, respectively; no public signature or
failure behavior changes.

Clone-independence tests cover editor constraint mutation, topology compaction,
transaction rollback with snapshots of both the entry and intermediate states,
and normalization. Contradictory normalization is checked with both unique and
shared storage.

**Verification:** the graph-IR unit, integration, property, and doc tests passed
(10,470 passed, 7 ignored), including 433 property-target checks with 256
generated cases and seed 20261013. The final expansion of the normalization
error cases passed all eight focused checks. Graph-IR all-target Clippy with
proptest and warnings denied, touched-file nightly formatting, git diff --check,
and full diff review passed. The contract inventory and logs are in
scratch/236-s5a/.

#### S5b — Reaction deltas **Done**

**Modules:** graph-IR reaction. **Green.** [dep: none]

Change the field to `Arc<Deltas>`. Reaction::try_new wraps the supplied Deltas.
normalize wraps its rebuilt Deltas. reframe_by, reframe, and into_parts use
Arc::unwrap_or_clone. The application iterator keeps its own Deltas field.

**Tests:** migrate the test literals and assignments that set the field.
Add clone-independence cases through normalize and reframe, and into_parts of a
shared and of a uniquely owned reaction.

**Completed (2026-10-13).** Reaction stores its deltas in Arc. normalize,
reframe_by, reframe, and into_parts use Arc::unwrap_or_clone, preserving the move
path for uniquely owned deltas and copying shared deltas before transformation
or extraction. Borrowing application paths explicitly clone the inner Deltas
when preparing them. ReactionApplicationIter retains its owned prepared Deltas.
Constructors, accessors, errors, and Python component ownership are unchanged.

The six overlay-removal normalization cases and the reframing change-chain case
now check both unique and shared delta storage, including preservation of a
retained clone. into_parts checks exact parts, buffer preservation when unique,
and independence when shared; its empty case remains covered.

**Verification:** 186 reaction unit tests passed, with 3 ignored. The full
graph-IR unit, integration, property, and doc suites passed: 10,483 passed,
7 ignored, no failures. The 433 property-target checks used 256 generated cases
and seed 20261013. Graph-IR all-target Clippy with proptest and warnings denied,
touched-file nightly formatting, git diff --check, and full diff review passed.
The contract inventory and logs are in scratch/236-s5b/.

#### S5c — ReactionSpan tables **Done**

**Modules:** graph-IR reaction_span. **Green.** [dep: S2b]

Change atoms, bonds, and constraints to `Arc<Vec<…>>`. Wrap at construction in
from_entries, try_remap, and reverse_sides. normalize, reframe_by, reframe, and
the consuming try_remap rebuild each table with Arc::unwrap_or_clone, moving
uniquely owned tables and cloning shared ones. The borrowing entries accessor
clones the inner Vec where it needs an owned table.

**Tests:** clone-independence cases through normalize and reframe.

**Completed (2026-10-13).** ReactionSpan's atom, bond, and constraint tables are
Arc-backed. Consuming normalization, reframing, and remapping move uniquely owned
tables with Arc::unwrap_or_clone and copy shared tables before rebuilding them.
reframe keeps constraints owned through normalization and transport, wrapping
only the final result. reframe_by rebuilds constraints while retaining the atom
and bond tables. Construction and side reversal wrap their output tables; the
borrowing entries accessor clones the inner vectors. ReactionSpanEntries,
public slice accessors, errors, and Python ownership are unchanged.

Existing normalization, frame-transport, and reframing cases now check unique
and shared storage against independent originals. They cover actual atom, bond,
and constraint rewrites, shared table buffers, and preservation of the atom
buffer through frame transport. Existing remapping and inverse-remapping tests
cover both shared and unique consumption.

**Verification:** all 135 reaction-span unit tests passed. The full graph-IR
unit, integration, property, and doc suites passed: 10,486 passed, 7 ignored,
no failures. The 433 property-target checks used 256 generated cases and seed
20261013. Graph-IR all-target Clippy with proptest and warnings denied,
touched-file nightly formatting, git diff --check, and full diff review passed.
The contract inventory and logs are in scratch/236-s5c/.

#### S5d — ReactionSpan overlay span newtypes **Done**

**Modules:** graph-IR dative, aromatic, multicenter, noncovalent, stereo,
reaction_span.
**Green.** [dep: none]

Wrap the relation set of DativeBondSpans, AromaticSystemSpans,
MulticenterBondSpans, NoncovalentBondSpans, StereoAtomSpans, and
StereoBondSpans in an Arc. new, try_map, and remap wrap their result.
into_entries uses Arc::unwrap_or_clone. The Normalize, FrameTransport, and
reframe implementations write through Arc::make_mut.

**Tests:** migrate the unit tests that permute participants directly. Add a
clone-independence case per type through normalize and reframe.

**Completed (2026-10-13).** All six overlay span collections wrap their relation
sets in Arc. normalize, reframe_by, and reframe obtain mutable storage once per
operation, copying only when shared. into_entries moves uniquely owned relation
sets through Arc::unwrap_or_clone. new and the mapping methods wrap their results;
mapping retains graph-core's existing borrowing rebuild without a preliminary
relation-set copy. ReactionSpan::reframe moves its overlay fields directly,
avoiding six temporary default Arc allocations. Public interfaces, validation,
failure behavior, and Python ownership are unchanged.

The existing tests now check clone independence for every collection through
normalization, direct and tracked reframing, and frame transport. Normalization
covers both unique and shared storage with full expected collection values;
reframing also checks initial storage sharing and retains independent originals.

**Verification:** all 68 focused span tests passed. The full graph-IR unit,
integration, property, and doc suites passed: 10,492 passed, 7 ignored, no
failures. The 433 property-target checks used 256 generated cases and seed
20261013. Graph-IR all-target Clippy with proptest and warnings denied,
touched-file nightly formatting, git diff --check, and full diff review passed.
The contract inventory and logs are in scratch/236-s5d/.

#### S5e — Benchmark cases **Done**

**Modules:** graph-IR benches. **Additive, green.** [dep: S5a–S5d]

Add the aggregate_storage bench target with the twelve measured cases.

**Completed (2026-10-13).** Added aggregate_storage with all twelve construction
and clone cases, using inline Rust fixtures and public constructors. Input
preparation and result destruction are excluded from constructor timing; clone
timing includes destruction of the clone. The production comparison against the
commit before S5a, including both run estimates, is recorded in the copy-on-write
storage section above. No production implementation or public API changed.

**Verification:** all twelve cases built and executed on both revisions, followed
by four full timing runs. The S5 workspace test gate passed with graph-IR's
proptest feature enabled: 22,421 passed, 10 ignored, no failures. The 433 graph-IR
property-target checks used 256 generated cases and seed 20261013. Workspace
all-target Clippy with the same feature and warnings denied passed. After
maturin develop in the repository's Python 3.13.16 environment, Python passed
2,204 tests with 2 skipped. Touched-file nightly formatting, git diff --check,
and full diff review passed. Artifacts and commands are in scratch/236-s5e/.


### S6 — Development guides

#### S6a — Ownership rules in the guides **Done**

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

**Completed (2026-10-13).** The data-type guide records consuming and borrowing
contracts, two-operand exceptions, borrowed semantic equalities, atomic molecule
operations, and aggregate copy-on-write storage. Recovery covers both
transactions and complete planning before mutation. The nomenclature guide
defines the _into suffix and both update meanings, and records combine's
consuming receiver.

**Verification:** checked the rules and named examples against the current Rust
signatures and implementations. Full diff review and git diff --check passed.
Documentation only; no code or tests changed.

#### S6b — Moved to 181 S12b

Manuscript listing tests and the read-only inline API review are part of
[181 S12b](181-python-boundary-ownership-2026-08-03.md#s12b--final-integration-and-closeout).

### S7 — Verification

#### S7a — Final gate **Done**

**Modules:** workspace. **Verification.** [dep: S0–S5, S6a]

Run `cargo +nightly fmt --all` over the full workspace in this final pass.
Formatting-only changes outside the implementation's scope are acceptable here.
Run workspace tests and lints with all features, and the Python suite in the
repository's Python 3.13 environment. Rerun benchmarks only for cases whose
production code differs from what was measured.

**Completed (2026-10-13).** Full workspace nightly formatting grouped imports in
six graph-core relation modules. No behavior changed, so the recorded benchmark
results were retained without rerunning them.

**Verification:** cargo test --locked --workspace --all-features --no-fail-fast
passed: 37,239 tests, 10 ignored, no failures, including the feature-gated
property and conformance suites and doctests. The run set PROPTEST_CASES=256
and PROPTEST_RNG_SEED=20261013. Workspace all-target Clippy with all features
and warnings denied passed. After maturin develop --locked --all-features,
the Python 3.13.16 suite passed: 2,204 tests, 2 skipped. The final nightly format
check, git diff --check, and full diff review passed. The Rust toolchain was
1.96.0; source revision f41db411e plus the import formatting. Logs, commands,
test totals, and environment details are in scratch/236-s7a/.

**Critical paths:** S3a1 → S3b → S3d → S3l → S4j, then S4c → S4f → S4n
and S4l → S4n. Both join S6a → S7a, after the remaining stage dependencies.
S4a → S4a1 → S4b → S4b1 establishes the corrected transformer baseline
before S4d, S4e, and S4k.
The binding path is S3b → S3b1 → S3b2 → S3b3 → S7a.
The cross-document sequence is S7a → close 236 → 181.
