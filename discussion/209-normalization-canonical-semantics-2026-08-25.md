# 209 — Normalization and canonical semantics

Status: In Progress
Date: 2026-08-25
Relates: [168](168-api-hygiene-2026-07-27.md),
[186](186-molecule-canonicalization-2026-08-05.md),
[208](208-canonicalization-scaling-2026-08-24.md),
[210](210-relation-frame-storage-2026-08-25.md),
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
changes to relation storage are tracked separately in doc 210. This document defines normalization
against the current storage while avoiding assumptions that would prevent that later migration.

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
x.remap(canonicalize_with_correspondence(x).1) == canonicalize(x)
```

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

`Deltas::normalize` retains its existing fixed-frame fold and normal form. It preserves input order
within each entity's operation chain while folding that chain, then canonically orders the
independent results. Entity `Add`/`Remove` lifecycles retain their existing cancellation and
contradiction rules. Repeated molecule-level constraint additions and removals are deduplicated,
and a matching addition and removal cancel. This work does not redesign those semantics; it only
requires aggregate normalization to express participant-bearing delta values in the coordinated
normalized frames before applying the existing fold.

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

Doc [210](210-relation-frame-storage-2026-08-25.md) owns the separate migration from eager
aromatic and multicenter participant ordering in graph-core to explicit frame transport in
graph-IR. Its callback removal, cross-operation audit, raw structural-equality consequences, and
correspondence-equivalence naming are not completion conditions for doc 209. The normalization
contract settled here must remain valid before and after that migration.

### Settled normalization laws

Every stereo entity receives one deterministic ligand-frame presentation, including an entity with
an undetermined configuration. The residual stabilizer contains exactly the actions that permute
equal-kind virtual-ligand occurrences, preserve structural incidence, and are allowed by every
asserted `StereoKind` on that frame. Actual-atom ligands remain distinguished by id. Individual
frame-relative constraints must be invariant under a generating set of this stabilizer.

For supported aggregates `x` and `y` whose normalization succeeds:

```text
normalize(normalize(x)) == normalize(x)
x.equiv(y) iff normalize(x) == normalize(y)
```

The second law uses structural equality of the normal forms. Inputs differing only by an admissible
local participant-frame restatement therefore normalize to the same stored value. Normalization
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

#### S1b — Enforce stereo-site incidence integrity

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

#### S2a — Normalize complete non-stereo entity entries

**Module:** new `umol-graph-ir/src/ir/molecule/normalize.rs` and its unit tests.

Add the owned molecule-normalization kernel beside `Molecule`. Normalize complete atom, bond,
dative, aromatic, multicenter, and noncovalent entries rather than cloning attribute forms into
temporary vectors. Preserve entity ids and topology, use the current relation frames, and mutate
copy-on-write stores only when required.

**Tests and evidence:** Cover every non-stereo family, nonuniform aromatic and multicenter electron
counts, already-normal inputs, and a leaf contradiction. Assert topology and ids are unchanged,
normalization is idempotent, and shared input stores are not mutated. Add focused normalization
benchmarks for already-normal and non-normal owned molecules.

**Change class:** additive private normalization kernel (green).

**Dependencies:** [dep: S0c]

#### S2b — Normalize complete stereo entries

**Module:** `umol-graph-ir/src/ir/molecule/normalize.rs`, `stereo.rs`, stereo constraints, and their
unit tests.

Extend the kernel to select one stereo frame action, reframe the ligand list, and transport the
configuration and inline constraints together. Compute the residual stabilizer from equal-kind
virtual ligands, structural incidence, and every asserted kind. Check each frame-relative
constraint independently with stabilizer generators and report `Contradiction` when it is not
invariant. Apply the same frame normalization to undetermined configurations.

**Tests and evidence:** Exercise nonidentity frame changes with nonuniform position-sensitive
constraints, inverse frame roundtrips, each unrestricted kind, both restricted kinds, undetermined
configuration, same-kind virtual ligands, and a non-invariant constraint. Do not use uniform
payloads or pre-sorted frames as the only evidence.

**Change class:** additive private stereo-entry normalization (green).

**Dependencies:** [dep: S1a, S1b, S2a]

#### S2c — Publish molecule normalization and coordinated constraints

**Module:** `umol-graph-ir/src/ir/molecule.rs`, `molecule/normalize.rs`, `canonicalize.rs`, and
molecule property tests.

Implement `Normalize` for `Molecule`. Apply every selected stereo frame action to the entity entry
and all molecule-level constraints that refer to it, then normalize the complete constraint store.
Replace the private canonicalization-only molecule normalizer with this operation and make
`Molecule::equiv` compare the resulting complete normal forms in the current id frame.

**Tests and evidence:** State and test `normalize(normalize(x)) == normalize(x)` and
`x.equiv(y) iff normalize(x) == normalize(y)`. Cover coordinated inline and molecule-level stereo
constraints, structural equality after normalizing an admissible reframing, contradiction behavior,
and canonicalization convergence. Retain the stronger structural-equality assertion rather than
only rechecking `equiv`.

**Change class:** public trait implementation and semantics-preserving integration (green).

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

**Tests and evidence:** Compile the graph-IR library and confirm that `Canonicalize` exposes only
`canonicalize`, `canonicalize_with_correspondence`, `canonical_hash`, and `canonical_eq`. The next
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

**Tests and evidence:** Search the living guides for `DescriptionLevel`, `description_level`, and
the removed `_by` names; every remaining occurrence must describe another current API rather than
the retired surface. Run `git diff --check`.

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

After completed S0, S1a, S1b, and S2a are independent additive prerequisites that join at S2b.
S3a and S3b then branch from S2c and join at S3c. The remaining critical path is
`S2b -> S2c -> S3a/S3b -> S3c -> S4a -> S4b -> S4c -> S4d -> S5a -> S5b -> S5c`.
Canonical-hash allocation, key allocation, orbit pruning, and prefix pruning remain in doc 208 and
are not completion conditions here. The relation-storage redesign remains independently deferrable
in doc 210.
