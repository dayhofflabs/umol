# 209 — Normalization and canonical semantics

Status: In Progress
Date: 2026-08-25
Relates: [168](168-api-hygiene-2026-07-27.md),
[186](186-molecule-canonicalization-2026-08-05.md),
[208](208-canonicalization-scaling-2026-08-24.md),
[210](210-relation-frame-storage-2026-08-25.md),
[211](211-relation-frames-and-api-2026-08-26.md),
[data-type guide](../docs/development/data-types.md),
[nomenclature guide](../docs/development/nomenclature.md)

## Purpose

Canonicalization optimization exposed unresolved relationships among structural hashing,
normalization, equivalence, canonical keys, and description levels. This document separates public
semantic changes from optimizations that must preserve the current semantics. The public boundary
is settled below; allocation and key representation remain non-blocking implementation questions.

The immediate scope is the shared `Canonicalize` implementation for `Molecule`, `Reaction`, and
`ReactionSpan`. This includes private effective-level inspection for their stored forms; an
independent redesign of any of the three representations remains out of scope. Normalization-driven
changes to relation storage are tracked in doc 211, which owns the relation API and the
frame-transport operation used here and supersedes doc 210. This document defines normalization so
that it is valid both under the current eager-ordering storage and after that migration.

`Molecule: Hash` is present on the development trunk but has not appeared in a release. The hash
scope must therefore be settled and implemented before the next release. Neither its current
implementation nor its numeric output carries a released compatibility obligation.

## Settled public boundary

Canonicalization is complete-only in the public API. The public `DescriptionLevel`,
`Molecule::description_level()`, `canonicalize_by`, `canonical_hash_by`, and `canonical_eq_by`
introduced during doc 208 are to be removed from Rust and Python. No public `normalize_by`,
`equiv_by`, or `equiv_under_by` operations are added.

The public `Canonicalize` contract retains:

```text
canonicalize
canonicalize_with_correspondence
canonical_hash
canonical_eq
```

The canonicalizer retains a private `CanonicalizeLevel` with `Topology`, `Constitution`,
`Structure`, and `Full`. Private inspection selects the lowest level required by each stored
aggregate:

- a molecule requires the greatest level of its populated entity families and constraints;
- a reaction requires the greater level of its lhs and every delta;
- a reaction span requires the greatest level of its populated entity-span families and constraint
  spans.

Atoms and localized bonds have base level `Topology`; non-stereo overlays have base level
`Constitution`; stereo entities have base level `Structure`. An explicit constraint delta or span,
a `ModifyConstraint` delta, or any carried entity form with a non-empty inline constraint store
requires `Full`. An entity span inspects every carried side, including both sides of `Modified`.
These inspectors are canonicalization machinery, not properties or operations of the graph-IR
representations. There is no public way to request, inspect, or force a level.

Each aggregate's `canonicalize` and `canonicalize_with_correspondence` use the private level selected
for their operand. `canonical_eq` uses the greater private level required by its two operands. These
reductions are exact only when they preserve the complete operation's value and failure behavior.

This also confines normalization and equivalence to their complete existing meanings. Complete
canonicalization normalizes the complete returned aggregate. `equiv` and `equiv_under` continue to
compare complete graph-IR values; they do not acquire level-dependent variants.
`canonical_eq` retains its existing complete comparison, contradiction, and integrity behavior;
this document reopens its allocation strategy, not its relation.

For `Molecule`, `canonicalize_with_correspondence` must satisfy:

```text
canonicalize_with_correspondence(x).0 == canonicalize(x)
reframe(x).remap(canonicalize_with_correspondence(x).1) == canonicalize(x)
```

The second law relates `canonicalize` to the level below it, not to the raw input. The earlier form
`x.remap(c) == canonicalize(x)` was false whenever `x` was not already reduced and reframed, because
`remap` relabels ids and does neither. S5a already stated the correct relation against a normalized
source; the discrepancy was a product of the partial nesting corrected below.

It need not reproduce the correspondence that an inaccessible forced-`Full` path would select
among symmetry-equivalent alternatives. That comparison is useful only as internal evidence.

## Stability boundary

Canonicalization is deterministic within a fixed umol version and context, but canonical forms,
selected correspondences, and canonical hash values are not stable between versions. They are not
persistent identifiers.

The private level selection, typed comparison schema, and search implementation may therefore
change between versions without violating the public contract. A future level-selecting API, if
ever justified, would be a new operation with newly stated semantics; the complete-only API makes
no advance promise about such projections.

## Canonical hash

`canonical_hash` remains public and hashes the complete canonical aggregate. The accepted baseline
is explicit complete canonicalization followed by the aggregate's structural `Hash`. It may
materialize the canonical aggregate.

Avoiding that materialization later is an optimization, not a new operation. A virtual canonical
view or streaming hash is acceptable only if it preserves the result and failure behavior of
canonicalizing and then hashing under the same release, context, and hasher. The hash remains a
non-persistent, collision-prone value with no cross-release numeric stability guarantee.

There is no public level-constrained structural or canonical hash. Replacing the derived structural
`Hash` implementation is unnecessary unless later measurements justify it for the complete
operation. Allocation, key representation, and complete-operation performance belong to doc 208
and do not block this API correction.

## Normalization design reopened

The complete-operation properties exposed a missing semantic layer between normalization of an
entity's attribute form and canonicalization of an indexed aggregate. The existing private
`normalize_molecule` is not an adequate implementation of that layer: it clones and normalizes
every attribute form, then writes the results back family by family. It neither defines entity
normal form nor coordinates participant-frame changes with external constraints. Its placement in
`canonicalize.rs` also conceals a general graph-IR operation inside one canonicalization path.

Three existing operations must be distinguished:

- `Normalize` folds leaf values, expressions, attribute forms, and constraint collections;
- graph-core relation construction canonicalizes unordered participant factors and transports
  position-indexed relation data through the same permutation; and
- the private `normalize_reaction_span_entries` only collapses an equivalent `Modified` span to
  `Unchanged`. It does not normalize the values and should eventually be named for that reduction.

`Deltas::normalize` retains its existing fold and normal form, but not its fixed-frame premise. That
premise holds today only because `Delta::remap` re-sorts an aromatic or multicenter atom list and
permutes the attributes through the resulting order (`delta.rs:2631`, `:2645`, `:2672`, `:2686`).
Doc [211](211-relation-frames-and-api-2026-08-26.md) S5a makes `Delta::remap` frame-preserving, so
frame selection for a delta becomes this operation's responsibility, discharged by S3a below. The
fold itself is unchanged: it preserves input order
within each entity's operation chain while folding that chain, then canonically orders the
independent results. Entity `Add`/`Remove` lifecycles retain their existing cancellation and
contradiction rules. Repeated molecule-level constraint additions and removals are deduplicated,
and a matching addition and removal cancel. This work does not redesign those semantics; it only
requires aggregate normalization to express participant-bearing delta values in the coordinated
normalized frames before applying the existing fold.

### Nesting correction

Doc [211](211-relation-frames-and-api-2026-08-26.md) settles the operation layer this document sits
in. Only prefixes of the pipeline are meaningful, so the three operations are nested:

```text
normalize     = reduce
reframe       = reduce + frame
canonicalize  = reduce + frame + id
```

`Normalize` is the **reduction only**: it folds value expressions, normalizes set representations,
and deduplicates constraints. It does not select participant frames. Frame selection is `reframe`,
whose surface and laws are in doc 211. `canonicalize` is unchanged as the outermost composite and
keeps its familiar name.

`equiv` and the `Equiv` trait retire. `equiv` mixed the reduction and the frame quotient under one
name, which is why it means something different on a form than on a molecule. Its replacement is
`normalized_eq`, answering the reduction question only, with `framed_eq` and `canonical_eq` for the
outer levels.

`Normalized<T>` is removed in the same pass. It has no construction site and no use anywhere in the
workspace: six mentions total, being its definition, its impl block, a module doc line, the
re-export, and two unrelated comments containing the word. The semantic-deduplication key it was
built for was never taken up. No `Canonical<T>` counterpart is added; its only justification was the
analogy to a type that never earned its keep.

Where the sections below say that entity or aggregate normalization selects a participant frame,
read `reframe`. The semantics recorded there — the residual stabilizer, the admissible frame
actions, the coordinated transport of inline and molecule-level constraints — are unchanged; only
which operation owns them moves. Every subitem naming `Normalize` for frame work is respecified
accordingly.

### Normalization decomposition

The following decomposition is settled for implementation.

**Leaf normalization** retains the current `Normalize` semantics: fold expressions, normalize
set-like values and constraints, and report contradictions. It has no entity ids or participant
frames.

**Entity normalization** operates on a complete entity entry, not only its attribute form. It
preserves the entity id, selects a deterministic local participant presentation wherever order is
not semantic, transports every position-dependent value through that change, and then normalizes
the attributes. The family-specific consequences are:

- atoms normalize `AtomForm`;
- localized bonds retain graph-core's canonical endpoint order and normalize `BondForm`;
- dative bonds canonicalize their donor set and normalize their form;
- aromatic systems and multicenter bonds sort their participants, apply the same permutation to
  their electron-count vectors, and normalize their forms;
- noncovalent bonds canonicalize their endpoint pair and normalize their form; and
- stereo atoms and stereo bonds normalize the complete ligand-frame-relative assertion, including
  configuration and frame-relative constraints.

The current `RelationData` implementations show no additional entity family with a
position-dependent payload. The same difficult entities recur in deltas, reactions, and reaction
spans because those aggregates carry alternative representations of the entries.

**Aggregate normalization** preserves entity ids and graph topology while normalizing complete
entries and aggregate constraints together. A stereo entity normalizer must expose the selected
local frame permutation so that `Molecule` can apply it to molecule-level constraints referring to
that entity. This local frame action is not an entity-id correspondence and does not by itself
justify extending `MoleculeCorrespondence`.

`ReactionSpan` must apply one shared local frame action to every carried side of an `EntitySpan`
and to its constraint spans; independently normalizing lhs and rhs values would be incorrect.
`Reaction` likewise requires coordinated normalization of its lhs and deltas. Existing entities in
the delta sequence use lhs frames, while added entities carry new frames, so independently
normalizing `lhs` and `deltas` is insufficient.

### Owned normalization

`Normalize::normalize(self)` already has the right ownership contract for an efficient
implementation. A consumed aggregate may be changed progressively; if normalization finds a
contradiction, the partially processed value is dropped and no rollback is required. Molecule
normalization should live with the molecule implementation and operate directly on its
copy-on-write stores. `Arc::make_mut` then clones an entity-family store only when that store is
shared, rather than cloning every form into temporary vectors unconditionally.

A public fallible `normalize_in_place(&mut self)` is not currently proposed. Such an operation
would require rollback or would expose partial mutation after an error. Consuming normalization
provides the relevant no-unconditional-clone path without that contract.

### Stereo-frame integrity boundary

Normalization may operate only on a stereo frame whose participants can be interpreted relative to
its site. Under the current compact ligand representation, this is part of `Molecule`
representation integrity rather than chemistry validation:

- for a stereo atom, every actual-atom ligand is a localized neighbor of the site and every
  implicit-hydrogen or lone-pair ligand is borne by the site atom; and
- for a stereo bond, each of the two consecutive ligand blocks belongs to one endpoint of the site
  bond. An actual-atom ligand is a localized neighbor of that endpoint other than the opposite site
  endpoint, and a virtual ligand is borne by that endpoint. The blocks may be exchanged as a whole;
  moving one ligand between them is not an admissible frame action.

`Molecule::check_integrity` must enforce these conditions through the checked and asserted
construction paths. A violation is a `MoleculeIntegrityError`, not a lattice `Contradiction`, and
normalization must not repair it. Constructors do not require a normalized entity frame.

This scope retains the existing atom-or-virtual-ligand representation and does not reopen the
general relation between topology and overlay entities. In particular, aromatic, dative,
multicenter, haptic, and other non-localized attachment semantics remain outside doc 209.

### Stereo normalization boundary

The current stereo-atom and stereo-bond relation storage remains unchanged in this scope. Stereo
bond endpoint blocks stay a semantic partition of the flat ordered ligand frame rather than becoming
a new stored or public type. Wrapping a ligand frame and form in another public aggregate would not
solve the external-constraint action and is not required for the integrity checks above.

The proposed virtual-ligand quotient needs a bond-specific qualification. For a stereo atom,
same-kind virtual ligands share one bearing site. For a stereo bond, a lone permutation of virtual
ligands borne by opposite endpoints is not generally an admissible cis/trans or axial frame action.
The admissible quotient is therefore the intersection of same-kind virtual permutations, the
structural-incidence stabilizer, and the stereo kind's parent group, rather than an unconditional
product of symmetric groups across the complete bond frame.

The admissible frame-action group belongs to the stereo kind's `CosetSpace`, not to `StereoAtom`,
`StereoBond`, or normalization-specific helpers. `Tetrahedral`, `SquarePlanar`,
`TrigonalBipyramidal`, and `Octahedral` admit their full symmetric groups; `CisTrans` and `Axial`
admit the eight actions of `S2 wr S2`: swaps within either two-ligand block and exchange of the two
blocks. The existing `CosetSpace` parent group stores this data. `umol-perm` publicly exposes:

```rust
pub fn allows(&self, permutation: Permutation) -> bool;
pub fn normalizer<T: Ord>(&self, frame: &[T]) -> Option<Permutation>;
```

`allows` returns `false` for the wrong degree or an action outside the kind's parent group and does
not fail or panic. `normalizer` returns `None` exactly when the frame length differs
from the space degree. A matching frame always has an allowed normalizing permutation. Graph IR
uses that successful path after its stereo-frame integrity gate. Neither operation exposes group
enumeration for normalization. The four unrestricted kinds sort the complete frame once and derive
that sorting permutation, without a group-membership test.
`CisTrans` and `Axial` sort within each two-position block and then order the two blocks, again
deriving one action rather than searching their eight-element parent group. An undetermined
configuration does not skip frame normalization: the ligand frame is normalized by the same
operation, the undetermined configuration remains fixed, and all frame-relative constraints are
transported. With no asserted kind, no kind-specific restriction is applied. Kinded stereo
entities delegate to their `CosetSpace` and do not branch on individual kinds.

Every frame-relative constraint must be invariant individually under this residual group. Other
constraints cannot collectively make a non-invariant constraint representable. Complete
stereo-entry normalization reframes the ligand list, transports the configuration and each
constraint through the same action, and rejects any constraint not fixed by every residual action.
Residual invariance is checked using generators such as adjacent swaps of equal ligands, not by
walking every element of the stabilizer. For this scope, that failure is reported as
[`Contradiction`]. Replacing the normalization error type and correcting failure comparison belong
to doc [168](168-api-hygiene-2026-07-27.md).

Normalization evidence must guarantee a nonidentity admissible frame action and a nonuniform
position-sensitive payload. The required laws include inverse frame roundtrips, convergence under
normalization, idempotence, and coordinated transport of inline and molecule-level stereo
constraints. General molecule generators that sort participants or use uniform payloads do not
exercise this domain.

### Deferred relation-storage redesign

Doc [211](211-relation-frames-and-api-2026-08-26.md) owns the migration from eager aromatic and
multicenter participant ordering in graph-core to explicit frame transport in graph-IR, absorbing
what doc [210](210-relation-frame-storage-2026-08-25.md) previously scoped. Its callback removal,
cross-operation audit, and raw structural-equality consequences are prerequisites for S2b rather
than completion conditions here. The normalization contract settled in this document is valid
before and after that migration: under eager ordering an aromatic frame is already sorted, and
after the migration this document's frame selection sorts it, so the normal form is the same.

### Settled normalization laws

Every stereo entity receives one deterministic ligand-frame presentation, including an entity with
an undetermined configuration. The residual stabilizer contains exactly the actions that permute
equal-kind virtual-ligand occurrences, preserve structural incidence, and are allowed by every
asserted `StereoKind` on that frame. Actual-atom ligands remain distinguished by id. Individual
frame-relative constraints must be invariant under a generating set of this stabilizer.

For supported aggregates `x` and `y` whose normalization succeeds:

```text
normalize(normalize(x)) == normalize(x)
normalized_eq(x, y) iff normalize(x) == normalize(y)

reframe(reframe(x)) == reframe(x)
framed_eq(x, y) iff reframe(x) == reframe(y)
```

Both equality laws use structural equality of the respective forms. Inputs differing only by an
admissible local participant-frame restatement are `framed_eq` but not `normalized_eq`, because the
reduction does not select a frame. Normalization
preserves entity ids and graph topology. It reports `Contradiction` when carried value normalization
fails or a frame-relative constraint is not representable under the residual stabilizer.

`Molecule` applies each selected local frame action to the complete entity entry and every
molecule-level constraint referring to it. `Deltas` retains its existing fixed-frame normal form.
`Reaction` coordinates lhs frame actions with every delta using those frames and normalizes an added
entity in its own frame before applying the existing delta fold. `ReactionSpan` selects one frame
action for every entity span and applies it to every carried side and associated constraint span
before reducing a normalized equivalent `Modified` span to `Unchanged`.

This work is not retroactively part of completed S0c. The revised plan begins its normalization
scope at S1.

## Staged implementation plan

Every stage ends with a green workspace. This plan removes the accidental level surface and installs
the private dispatch needed by doc 208, defines complete aggregate normalization, and then removes
the accidental public level surface. It does not select or implement a hash-specific optimization.

### S0 — Establish private aggregate dispatch

#### S0a — Restore the private canonicalization level and leaf inspection **Done**

**Module:** `umol-graph-ir/src/ir/canonicalize.rs` and its unit tests.

Add private `CanonicalizeLevel::{Topology, Constitution, Structure, Full}` and private
inspection for deltas and entity spans. Give each entity family its base level, inspect both carried
forms of a modified span, and raise any inline or explicit constraint change to `Full`. Keep the
public `DescriptionLevel` surface temporarily so this subitem is additive and green.

**Tests and evidence:** Use module-local `rstest` tables covering every delta family and operation,
every entity-span position, both sides of `Modified`, all inline constraint stores, and explicit
constraint deltas. Assert the private containment order without adding a public query.

**Change class:** additive private infrastructure (green).

**Dependencies:** [dep: none]

#### S0b — Derive the effective level of each aggregate **Done**

**Module:** `umol-graph-ir/src/ir/canonicalize.rs` and its unit tests.

Add private aggregate inspectors. A molecule takes the maximum required by its entity collections
and molecule constraints. A reaction takes the maximum of its lhs and ordered deltas. A reaction
span takes the maximum across all entity-span collections and its constraint spans. Binary
selection takes the maximum required by the two operands. Inspection is total and does not
validate, normalize, apply deltas, or materialize a reaction span.

**Tests and evidence:** Cover the four levels for each aggregate, including a reaction whose lhs
and delta require different levels and a modified entity span whose two sides require different
levels. Cover explicit constraint-span presence here. Assert that dense id remapping and delta
inversion preserve the selected level.

**Change class:** additive private infrastructure (green).

**Dependencies:** [dep: S0a]

#### S0c — Route complete aggregate operations through the private level **Done**

**Module:** `umol-graph-ir/src/ir/canonicalize.rs`, its unit tests, and aggregate canonicalization
properties.

Route `canonicalize` and `canonicalize_with_correspondence` for `Molecule`, `Reaction`, and
`ReactionSpan` through their selected private levels. Route each `canonical_eq` through the greater
level required by its operands. `canonical_hash` continues to canonicalize and then structurally
hash, so it inherits unary dispatch without a separate path.

**Tests and evidence:** Compare selected and forced-full internal results for representative
molecules, reactions, and reaction spans at every private level. Assert exact canonical aggregates,
unchanged contradiction and integrity behavior, aggregate-specific correspondence transport laws,
complete canonical-hash invariance, and asymmetric `canonical_eq` cases where only one operand
requires the higher level. Do not freeze the forced-full correspondence itself. Assert
`canonical_hash(x) == hash(canonicalize(x))` directly for all three aggregates.

**Change class:** semantics-faithful private dispatch (green).

**Dependencies:** [dep: S0b]

**Current state:** Effective-level routing, forced-`Full` representative comparisons, aggregate
hash checks, and current correspondence checks are implemented and green. The external property
suite cannot yet assert exact structural transport from the normalized source because complete
aggregate `Normalize` implementations do not exist. S5a retains that stronger property as an
explicit closeout condition after S2 and S3 provide them.

### S1 — Establish stereo frame actions and integrity

#### S1a — Add the public stereo-kind frame operations **Done**

**Module:** `umol-perm/src/coset.rs` and its unit tests.

Add `CosetSpace::allows` and `CosetSpace::normalizer`. `allows` checks degree and parent-group
membership. `normalizer` returns `None` only for a frame-degree mismatch;
otherwise it derives one allowed action by a complete-frame sort for unrestricted kinds, or by
within-block sorting followed by block ordering for `CisTrans` and `Axial`. Neither path enumerates
the parent group.

**Tests and evidence:** Use exact `rstest` tables for every kind, wrong-degree `allows` and
normalization inputs, nonidentity normal ordering, unrestricted-kind sorting, and the two restricted
block actions. Assert that every returned normalizing action is accepted by `allows` and that
selecting again from the resulting frame returns the identity action.

**Change class:** additive public operations (green).

**Dependencies:** [dep: S0c]

#### S1b — Enforce stereo-site incidence integrity **Done**

**Module:** `umol-graph-ir/src/ir/molecule/integrity.rs`, its public error type, and construction
tests.

Extend `Molecule::check_integrity` to require stereo-atom ligands to be borne by or adjacent to the
site and stereo-bond ligand blocks to belong to the corresponding bond endpoints. Preserve whole
block exchange as admissible and reject individual cross-endpoint movement. Route checked and
asserted construction through the same existing integrity implementation.

**Tests and evidence:** Cover actual and virtual ligands for stereo atoms and bonds, both endpoint
block orientations, a misplaced actual ligand, a wrongly anchored virtual ligand, and a ligand
moved across bond endpoint blocks. Assert the exact `MoleculeIntegrityError` variants and retain all
existing integrity cases.

**Change class:** strengthened representation-integrity contract (green).

**Dependencies:** [dep: S0c]


### S2 — Implement complete molecule normalization

#### S2a — Take the relation surface from doc 211

**Current state:** Settled. Doc [211](211-relation-frames-and-api-2026-08-26.md) selects the
relation contract this subitem was waiting for, and doc 210 is superseded there.

The withdrawal of `permute_participants_with` is **lifted**, narrowed, and renamed. Doc
[211](211-relation-frames-and-api-2026-08-26.md) S4a adds `permute_with`, and `permute_1_with` /
`permute_2_with` on the birelation shapes: one entry, a validated permutation, no payload access.
Because the multiset cannot change, incidence stays valid and is not rebuilt. The original proposal
was a general public mutation family put forward before the relation API had been reviewed; this is
the reviewed and much narrower form, and it is required because removing eager sorting removes the
only path that could reorder a stored frame. Complete entry normalization therefore applies the
selected action in place rather than through owned reconstruction, which would cost one heap
allocation per entity out of a flat CSR plus two index rebuilds.

The operations this subitem depends on are:

- the form-level `reframe_to` methods on all six forms, inherent rather than a trait, from doc 211
  S3b;
- `StereoAtomForm::select_frame` and `reframe_by` for the stereo case, from doc 211 S3b; and
- frame-preserving `new` and `into_entries`, from doc 211 S5b.

`Reframe` declines an ambiguous frame change rather than selecting a representative. Resolving that
ambiguity is this document's work, through `CosetSpace::normalizer` and generator-based
residual-invariance checking; doc 211 does not do it.

Doc 211 S5b ends with thirteen enumerated canonicalization and hash failures which this document
closes. They are `test_canonicalize_constitution_family_minimum`,
`test_canonicalize_constitution_participant_order`, `test_kindless_stereo_atom_frame_order`,
`test_kindless_stereo_bond_frame_order`, `test_minimum_kinded_stereo_frames` cases 1 to 4,
`reaction_span::test_reaction_span_canonicalize::case_2_constitution`, and the property tests
`reaction::canonicalize::test_reaction_canonical_eq_by`,
`reaction::canonicalize::test_reaction_canonical_hash`,
`reaction::span::canonicalize::test_reaction_span_canonical_hash`, and
`reaction::span::canonicalize::test_reaction_span_canonicalize`. Restoring them green is the exit
condition for S2 and S3 of this document.

Molecule-level stereo constraint transport under pushout is absorbed here from doc 210 and belongs
to S2b's aggregate constraint handling.

**Change class:** dependency resolution; no implementation in this subitem.

**Dependencies:** [dep: S1a, S1b, doc 211 S3g, doc 211 S5b]

#### S2b — Implement complete molecule normalization

**Module:** `umol-graph-ir/src/ir/molecule.rs`, a focused `molecule/normalize.rs` implementation
module if needed for file size, `stereo.rs`, constraints, and their unit tests.

**Respecified by the nesting correction.** This subitem implements two operations, not one:
`Normalize for Molecule` is the reduction over every entity family and every constraint store, and
`Reframe for Molecule` is the frame prefix that reduces and then selects each entity's frame. The
frame semantics below — the residual stabilizer, the admissible actions, the coordinated transport
of inline and molecule-level constraints — belong to the second, and the selected action reaches
molecule-level constraints through `reframe_with_action`.

Implement both directly as aggregate operations over every entity family.
Keep the aggregate control flow in that trait implementation; a separate source module may keep the
file readable, but must not introduce free-function kernels or a semantic split between stereo and
non-stereo entities. Preserve entity ids and topology, mutate copy-on-write stores through
`Arc::make_mut`, and drop the consumed partially normalized value if a leaf reports
`Contradiction`; no rollback is required.

Normalize ordinary entity forms in their stored frames. Under the current eager-ordering storage,
aromatic and multicenter entries already have a selected participant frame; normalize their
position-sensitive electron-count payload in that frame. Doc 210 later moves frame selection into
this operation and uses the relation contract selected through S2a. For stereo entries now, select
one admissible normalizing action and use the owned reconstruction seam selected in S2a to apply
that action to the ligand list while transporting the configuration and inline constraints. Apply
the same selected action to molecule-level constraints referring to the entry, then normalize the
complete constraint store. Undetermined configurations undergo the same frame normalization.

For stereo, compute the residual stabilizer from equal-kind virtual ligands, structural incidence,
and every asserted kind. Check each frame-relative constraint independently with stabilizer
generators and report `Contradiction` when it is not invariant. Do not enumerate a complete
permutation group during normalization.

**Tests and evidence:** Cover every entity family, nonuniform aromatic and multicenter electron
counts, already-normal inputs, and a leaf contradiction. Exercise nonidentity stereo frame changes
with nonuniform position-sensitive constraints, inverse frame roundtrips, every unrestricted kind,
both restricted kinds, undetermined configuration, same-kind virtual ligands, coordinated inline
and molecule-level constraints, and a non-invariant constraint. Assert that topology and ids are
unchanged, shared input stores are not mutated, and uniquely owned stores are not cloned solely for
normalization. Do not use uniform payloads or pre-sorted stereo frames as the only evidence. Add
focused benchmarks for already-normal and non-normal owned molecules.

**Change class:** additive public trait implementation (green).

**Dependencies:** [dep: S1a, S1b, S2a]

#### S2b.1 — Unify `Molecule::equiv_under` across entity kinds

**Module:** `umol-graph-ir/src/ir/molecule.rs` and its unit tests.

`Molecule::equiv_under` implements frame transport twice, by hand, with different shapes per entity
kind. Doc [211](211-relation-frames-and-api-2026-08-26.md) S4d replaces the machinery underneath
both but does not unify them, because the aggregate has no frame members to unify them onto until
S2b supplies `Reframe for Molecule`.

What the two paths are today:

- **The four distinct-participant families.** Map the participants through the correspondence, derive
  the single order from mapped-to-stored, restate, compare normal forms. After 211 S4d this reads as
  `reframe_to` composed with `equiv`. One target, one action.
- **The two stereo families.** `Permutation::between_all(mapped_ligands, stored_ligands)` enumerates
  every bijection, keeps those under which the configurations agree, and carries the surviving *set*
  forward in `stereo_frames` so `constraints_equiv_under_stereo_frames` can check molecule-level
  constraints against the candidates consistently across entities.

The difference is not arbitrary. The mapped-to-stored bijection is unique exactly when the
participants are distinct, which molecule integrity guarantees for the four; repeated ligands leave a
residual stabilizer, so for stereo it is a set. And only stereo carries frame-relative constraints —
the other four families' inline constraint forms are not position-indexed — so only stereo needs the
candidate set to reach the constraint store.

So the two shapes are both correct and neither generalizes the other by accident: the unique
bijection is the degenerate case of the candidate set. What is wrong is that they are two hand-written
transports rather than one operation with a degenerate case, which is the same defect 211 removed
everywhere else.

**The work.** Once `Reframe for Molecule` exists, express `equiv_under` as the aggregate's
supplied-target member — the `_to` reading of the frame quotient, where `reframe_with_action` is the
selected-target reading — with the candidate set as the general case and the unique bijection as its
degeneracy. Whether the name survives is settled here too and not before: `Molecule` gains `equiv`
from the blanket impl over `Normalize` in S2b, so the inherent `Molecule::equiv` must go, and only
then is it clear what the witness-taking operation should be called beside it.

**A pre-existing narrowing to resolve with it.** `Molecule::equiv_under`'s stereo path short-circuits
when the mapped ligand frame already equals the stored one and the configurations agree: it
`continue`s, which skips pushing the entity onto `stereo_frames`. That list is what
`constraints_equiv_under_stereo_frames` recurses over, trying every candidate permutation against the
molecule-level constraints. An entity omitted from it has its constraints compared under the identity
alone.

Falling through instead would offer identity *plus the frame's residual stabilizer* — the further
permutations under which the configurations still match, which exist exactly when ligands repeat. So
the short-circuit can reject a pair the full path accepts: identical ligand frames, matching
configurations, and molecule-level frame-relative constraints agreeing only under a non-identity
stabilizer element.

Not a regression — the branch predates doc 211 — and the failure is a false negative: it rejects a
pair it should accept, never the reverse.

**Constructed and confirmed**, rather than left as a conjecture.
`test_molecule_equiv_under_stereo_stabilizer_constraint` builds a stereo atom with two implicit
hydrogens at positions 0 and 1, an undetermined configuration so the swap preserves it, and molecule-
level `Topicity` constraints on the two sides differing exactly by that swap. `equiv_under` returns
`false`. Deleting the short-circuit and nothing else makes the same test pass, which isolates the
mechanism to that branch. The test is `#[ignore]`d, naming this subitem; unignore it when the
unification lands.

**Tests and evidence:** Retain every `equiv_under` case. Assert that the unified operation agrees with
both replaced paths on generated inputs, including a stereo entry with repeated virtual ligands where
the candidate set has more than one member, and a non-stereo entry where it has exactly one.

**Change class:** simplification with caller migration (green).

**Dependencies:** [dep: S2b, doc 211 S4d]

#### S2c — Integrate and verify molecule normalization

**Module:** `umol-graph-ir/src/ir/canonicalize.rs`, molecule property tests, and normalization
benchmarks.

Replace the private canonicalization-only molecule normalizer with the operations from S2b. Remove
the inherent `Molecule::equiv`: it compares field by field rather than comparing normal forms, so it
answers a different question under the same name, and `normalized_eq` and `framed_eq` replace it.
Correct its retired doc comment's claim of equality "in the current id and participant frame". Do
not add a second public normalization entry point or a free-function alias.

**Tests and evidence:** State and test idempotence and the equality law at both levels:
`normalize(normalize(x)) == normalize(x)` with `normalized_eq`, and `reframe(reframe(x)) ==
reframe(x)` with `framed_eq`. Assert the containment `normalized_eq => framed_eq => canonical_eq`. Assert structural equality after normalizing an
admissible reframing, contradiction behavior, and canonicalization convergence. Retain the stronger
structural-equality assertion rather than only rechecking `equiv`.

**Change class:** semantics-preserving integration and property evidence (green).

**Dependencies:** [dep: S2b]

### S3 — Extend normalization to reactions

#### S3a — Normalize reactions in coordinated frames

**Module:** `umol-graph-ir/src/ir/reaction.rs`, `delta.rs`, molecule normalization support, and
reaction property tests.

Implement `Normalize` for `Reaction`. Reuse the lhs molecule's selected frame actions for existing
entity deltas, normalize added entities in their own frames, transport frame-relative constraint
deltas, and then invoke the existing ordered `Deltas::normalize` fold. Do not change delta sequence,
cancellation, or contradiction semantics.

**Tests and evidence:** Cover existing and added stereo entities, coordinated constraint deltas,
noncommuting entity-delta chains, repeated add/remove normalization, idempotence, and equality of
normal forms for admissibly reframed reactions.

**Change class:** additive public trait implementation (green).

**Dependencies:** [dep: S2c]

#### S3b — Normalize reaction spans in one shared frame

**Module:** `umol-graph-ir/src/ir/reaction_span.rs`, molecule normalization support, and
reaction-span property tests.

Implement `Normalize` for `ReactionSpan`. Select one frame action per entity span and apply it to
every carried side and associated constraint span. Normalize carried values before reducing an
equivalent `Modified` span to `Unchanged`; rename the existing private reduction helper to describe
that narrower role.

**Tests and evidence:** Cover every span tag, a modified stereo entry whose sides require one shared
nonidentity action, coordinated constraint spans, idempotence, normalized structural equality, and
the existing span/reaction conversion laws.

**Change class:** additive public trait implementation and private-helper correction (green).

**Dependencies:** [dep: S2c]

#### S3c — Route aggregate canonicalization through complete normalization

**Module:** `umol-graph-ir/src/ir/canonicalize.rs`, its unit tests, and aggregate canonicalization
properties.

Remove the remaining canonicalization-local normalization paths. Route molecule, reaction, and
reaction-span canonicalization through their owning `Normalize` implementations while retaining
private effective-level selection and the existing public canonicalization errors.

**Tests and evidence:** Re-run exact canonical aggregates, correspondence transport, renumbering,
contradiction, integrity, canonical equality, and canonical-hash laws for all three aggregates.
Include nonidentity stereo-frame and nonuniform-constraint cases.

**Change class:** semantics-preserving integration (green).

**Dependencies:** [dep: S3a, S3b]

### S4 — Remove the public level surface

#### S4a — Retire the Python consumers

**Module:** `umol-py/src/canonicalize.rs`, molecule bindings, package exports, and Python tests.

Remove Python `DescriptionLevel`, `Molecule.description_level`, and every `canonicalize_by` and
`canonical_eq_by` method from molecule, reaction, and reaction-span bindings. Remove their package
exports and signature inventory entries. The Rust provider still exists during this subitem, so the
Python package can return to green before the Rust surface is removed.

**Tests and evidence:** Rebuild the extension with the repository Python 3.13 environment and run
the focused import, molecule, reaction, and reaction-span tests. Assert the complete operations and
their existing exceptions rather than replacing removed projection tests.

**Change class:** breaking Python API removal (green after its caller and test migration).

**Dependencies:** [dep: S3c]

#### S4b — Retire the Rust level API

**Module:** `umol-graph-ir/src/ir/canonicalize.rs`, `molecule.rs`, and the graph-IR root exports.

Remove `canonicalize_by`, `canonical_hash_by`, and `canonical_eq_by` from `Canonicalize` and all
aggregate implementations. Remove public `DescriptionLevel`, `Molecule::description_level`, and the
root re-export. Convert every retained internal level parameter to private `CanonicalizeLevel` and
remove reaction-projection helpers that existed only for the public reduced operations.

Complete the trait's quotient shape in the same pass, rather than editing it twice. `Canonicalize`
is the entity-id quotient, and a group quotient has four members: act by a supplied witness, select
a representative, select and expose the witness, and equality modulo the group.
`canonicalize` and `canonicalize_with_correspondence` are already the middle two.

- **Absorb `remap` as the act member.** It is inherent today only because the witness types were
  never layered. Note the consequence: `remap` exists on `Molecule` and `ReactionSpan` but **not on
  `Reaction`**, so `Reaction` gains one, transporting its normalized materialized span and converting
  back, which is the shape S5a already assumes.
- **Give `canonical_eq` a default body.** It is a required member today, so its meaning lives in
  three separate implementations, and `Molecule::canonical_eq` (`canonicalize.rs:5211`) does not
  canonicalize and compare at all: it short-circuits on `==` and then compares `canonical_key_by` at
  a selected level. That is a legitimate optimization, but it should override a default body that
  states the meaning, exactly as `canonical_hash` already does.
- **Derive `equiv_under`.** It is act-by-witness followed by equality one level down —
  `x.equiv_under(y, w)` is `x.remap(w)` compared to `y` — so it becomes a provided method rather
  than a hand-written one. `Molecule::equiv_under` has no non-test callers, is absent from the
  Python bindings, and appears only in canonicalization and molecule test assertions, so nothing
  depends on its hand-optimized form. It currently spans two quotient levels at once, fixing the id
  witness while *searching* frame witnesses for stereo through `Permutation::between_all`; under the
  layering that decomposes into act by the id witness, then frame-level equality.

The context stays a concrete `&CanonicalizeContext` parameter on this trait only. `Normalize` and
`Reframe` take no context because they are deterministic computations rather than searches; the
parameter enters exactly where configurable tie-breaking does. No associated context type is
introduced, because all three implementations use the same context. Splitting that type's semantic
and algorithmic halves is recorded in doc [168](168-api-hygiene-2026-07-27.md) and is deliberately
out of scope here.

The act member is `remap` alone. `try_remap` is removed at this layer and at the relation-set layer
together, so the two stay aligned. `remap` has six non-test call sites here and thirty-four at the
relation-set layer; `try_remap` has none at either. Its preconditions stay individually checkable
through `Molecule::check_integrity`, `Correspondence::is_total`, and the entity counts, so what is
removed is a pre-bundled check rather than a capability. Closure under laws is not an argument for
keeping it: that argument covers algebraic members, and checked-versus-asserted is an ergonomics
pair. Doc [211](211-relation-frames-and-api-2026-08-26.md) records the matching relation-set
removal.

**Tests and evidence:** Compile the graph-IR library and confirm that `Canonicalize` exposes only
the four quotient members plus `canonical_hash`. Assert the trait laws directly: idempotence of
`canonicalize`; that transporting through the exposed correspondence reproduces the canonical form;
that `canonical_eq` agrees with comparing canonical forms; and that `canonicalize` is invariant
under `remap`. Because `canonical_eq` gains a default body that `Molecule` overrides, assert on
generated inputs that the override and the default agree. The next
subitem migrates test and benchmark callers, so this breaking subitem may be red within S4.

**Change class:** breaking Rust API removal (red until S4c).

**Dependencies:** [dep: S4a]

#### S4c — Migrate Rust tests, properties, and benchmarks

**Module:** graph-IR canonicalization unit and property tests and
`umol-graph-ir/benches/canonicalize.rs`.

Remove representation-level and public projection tests. Retain complete canonicalization, hash,
equality, renumbering, contradiction, integrity, and correspondence laws. Keep forced-level
comparisons only inside the canonicalization module as private implementation evidence. Replace
benchmark calls to removed methods with the complete public operations needed for the target to
compile; doc 208 owns the expanded performance matrix.

**Tests and evidence:** Run graph-IR unit tests, the feature-gated canonicalization properties, and
the canonicalization benchmark build. Search Rust source outside the canonicalization module for
the removed enum and methods.

**Change class:** caller and verification migration (restores Rust green).

**Dependencies:** [dep: S4b]

#### S4d — Align the living development guides

**Module:** `docs/development/data-types.md`, `docs/development/nomenclature.md`, and
`docs/development/property-tests.md`.

Remove the description-level entry, suffix inventory member, retired-name rule, public trait
methods, and level-specific property claims. Describe aggregate canonicalization, canonical hash,
and canonical equality only through their complete operations. Retain structural domain and
incidence-level terminology where it remains independently meaningful.

Correct the four places where the guides state that normalization does not touch participant frames,
which this document reverses:

- `nomenclature.md:307` — "**Not:** *normalize*, which operates within an existing id and participant
  frame."
- `nomenclature.md:718` — "`equiv` — equality of normalized forms in the current id and participant
  frame."
- `nomenclature.md:1160` — "**Normalize** puts a form into a deterministic normal form without
  changing entity ids or participant frames."
- `nomenclature.md:1172` — "**Not:** aggregate canonicalization, which selects an entity and
  participant frame."

The correction is a scoping one, not a reversal. A form cannot see a participant frame, because the
frame lives with the participants in the relation set, so each statement stays true at the form
level and must say so. A relation set and an aggregate can see frames and select them; only
`Canonicalize` additionally selects entity ids, so the `canonicalize` entry at `:291` should stop
double-counting the frame.

Correct `Molecule::equiv`'s doc comment (`molecule.rs:469`), "Complete semantic equality in the
current id and participant frame". S2c makes it compare normal forms, so the participant clause is
wrong; the id clause remains. This also restores the law that `equiv` agrees with
`equiv_under` under the identity correspondence, which
`test_molecule_equiv_under_identity_reduces_to_equiv` asserts but cannot currently exercise.

**Tests and evidence:** Search the living guides for `DescriptionLevel`, `description_level`, and
the removed `_by` names; every remaining occurrence must describe another current API rather than
the retired surface. Search for "participant frame" and confirm every remaining occurrence names its
carrier. Run `git diff --check`.

**Change class:** documentation migration (green).

**Dependencies:** [dep: S4c]

### S5 — Verify and close the semantic correction

#### S5a — Assert normalized-source correspondence transport

**Module:** molecule, reaction, and reaction-span canonicalization property tests.

For every generated source whose normalization and canonicalization succeed, obtain
`(canonical, correspondence)` from `canonicalize_with_correspondence` and assert exact structural
equality after normalizing the source and transporting it through that correspondence. Use direct
`Molecule::remap` and `ReactionSpan::remap`; for `Reaction`, transport its normalized materialized
span and convert the result back to a reaction. Also retain
`canonicalize_with_correspondence(source).0 == canonicalize(source)` for all three aggregates.

**Tests and evidence:** The asserted relation is structural equality, not `equiv`, `canonical_eq`,
or equality after another canonicalization. Generated inputs include non-normal carried values and
nonidentity entity correspondences.

**Change class:** closeout property evidence (green).

**Dependencies:** [dep: S4d]

#### S5b — Run the cross-language verification gate

**Module:** the graph-IR canonicalization surface, `umol-py`, and all workspace callers.

Run formatting, graph-IR unit and feature-gated property tests, the canonicalization benchmark
build, graph-IR linting, the Python 3.13 build and tests, and workspace test and lint gates. Audit
the public Rust and Python inventories for the complete-only surface. `CanonicalizeLevel` may occur
only in canonicalization-private implementation and module-local tests.

**Tests and evidence:** Every gate passes; complete canonicalization and transport properties remain
green; no public documentation, export, binding, example, or benchmark requires a level selector.

**Change class:** verification only (green).

**Dependencies:** [dep: S5a]

#### S5c — Reconcile the discussion records

**Module:** docs 208 and 209 and `discussion/000-status.md`.

Record the implemented private-dispatch and API-removal outcome. Resume doc 208 from the
complete-only API, with canonical-hash measurement and any allocation work remaining there. Mark
doc 209 `Completed` only after S5b and keep cross-version stability explicitly unsupported.

**Tests and evidence:** Discussion links and statuses agree, the doc-208 next action names current
private and public surfaces, and `git diff --check` passes.

**Change class:** closeout documentation (green).

**Dependencies:** [dep: S5b]

### Dependency summary

S0 is complete, as are S1a and S1b, which completes S1. The two stereo-integrity rules that make
frame selection unambiguous moved to doc [211](211-relation-frames-and-api-2026-08-26.md) as its S3f
and S3g, so that document is implementable end to end without returning here; S2a depends on them
there. S2a also depends on doc 211 S5b, which must land before S2b begins. After S2c, S3a and S3b branch and join at S3c. The critical path from S2b is
`S2b -> S2c -> S3a/S3b -> S3c -> S4a -> S4b -> S4c -> S4d -> S5a -> S5b -> S5c`.
Canonical-hash allocation, key allocation, orbit pruning, and prefix pruning remain in doc 208 and
are not completion conditions here.

Doc 211 S5b leaves the workspace red on thirteen canonicalization and hash tests, enumerated in S2a.
S2 and S3 of this document restore them. That is a deliberate red period across the two documents,
not an unplanned regression.
