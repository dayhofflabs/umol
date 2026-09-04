# 211 — Relation frames and the relation API

Status: Completed
Date: 2026-08-27
Relates: [168](168-api-hygiene-2026-07-27.md),
[166](166-molecule-ops-2026-07-27.md),
[204](204-reaction-application-redesign-2026-08-19.md),
[208](208-canonicalization-scaling-2026-08-24.md),
[209](209-normalization-canonical-semantics-2026-08-25.md),
[210](210-relation-frame-storage-2026-08-25.md),
[212](212-remapping-layer-2026-08-26.md),
[213](213-editor-overlay-storage-2026-08-27.md),
[214](214-aggregate-frame-semantics-2026-08-28.md),
[215](215-integrity-minimization-2026-08-28.md),
[data-type guide](../docs/development/data-types.md),
[nomenclature guide](../docs/development/nomenclature.md)

## Supersession notice — 2026-08-28

Doc [214](214-aggregate-frame-semantics-2026-08-28.md) supersedes only this completed document's
decision to permit repeated virtual ligands and to interpret a stored coset through their residual
occurrence stabilizer. Current integrity requires pairwise-distinct, `MAX_DEGREE`-bounded stereo
frames; the orbit search, virtual-block swaps, `FrameAction`, and `find_reframed` machinery are
scheduled to be unwound. Frame-preserving relation storage, entry identity, participant transport,
and the remainder of the completed relation API work stay in force.

### Completion addendum — 2026-08-29

Doc 214 completed that unwind. Published stereo frames are now pairwise-distinct and bounded;
occurrence-orbit search, `FrameAction`, `find_reframed`, and their repeat-specific permutation
machinery are removed. Current graph-IR transport uses `FrameTransport` and `Reframe`, and current
fixed-frame comparison uses `normalized_eq`. The historical repeat-valid analysis below remains as
the record of the superseded decision.

## Purpose

The graph-core relation types carry construction-time participant ordering, a payload reindex
callback, participant lookup, participant alignment, graph-id transport, incidence, and relation
algebra on one set of storage types. Three mutually incompatible participant-alignment rules exist,
and the payload protocol that was added to repair construction-time sorting is implemented
non-trivially by two of six payload types while the hardest frame problem bypasses it entirely.

This document replaces the earlier API accounting with one selected design, absorbs the storage
migration previously scoped in doc [210](210-relation-frame-storage-2026-08-25.md), and carries the
staged implementation plan.

## What the current design complects

One concept — a **frame** (a relation's ordered participant presentation) and a **reframing** (a
bijection between two presentations of the same participant multiset) — is spread over three
carriers and four owners.

| carrier | degree | location |
| --- | --- | --- |
| `Vec<ParticipantPosition>` | unbounded | `FactorOrdering::canonicalize_positions` |
| `Permutation` | at most 6, `Copy` | `umol-perm`, stereo configuration |
| implicit, from both sides having been sorted at construction | — | aromatic and multicenter pushout |

The action has four owners: `RelationData::on_permutation`, `RelationEquiv::equiv_under`,
`StereoAtomForm::transform_frame_by`, and `ElectronCountsForm::permute`. Alignment is derived by
three different rules: sorting in `FixedRelationSet::participant_permutation`, multiset search in
the storage participant matcher, and first-position search in the two private `MoleculeEditor`
copies. A
reordered query can therefore be found by one and rejected by another.

### Evidence

Six entity families use three of the five storage shapes.

| family | shape | factor 1 | factor 2 | payload position-sensitive |
| --- | --- | --- | --- | --- |
| aromatic system | `VarRelationSet` | `NodeId` `Unordered` | — | yes, `electrons` |
| multicenter bond | `VarRelationSet` | `NodeId` `Unordered` | — | yes, `electrons` |
| noncovalent bond | `FixedRelationSet<2>` | `NodeId` `Unordered` | — | no, `on_permutation` is empty |
| dative bond | `FixedVarBirelationSet` | `NodeId` `Ordered` 1 | `NodeId` `Unordered` | no, `on_permutation` is empty |
| stereo atom | `FixedVarBirelationSet` | `NodeId` `Ordered` 1 | `StereoLigand` `Ordered` | no; frame transport bypasses the protocol |
| stereo bond | `FixedVarBirelationSet` | `EdgeId` `Ordered` 1 | `StereoLigand` `Ordered` | no; frame transport bypasses the protocol |

`RelationData` and `BiRelationData` therefore exist to serve construction-time sorting, which is the
behaviour this design removes. Stereo ligands are marked `Ordered`, so every stereo
`on_permutation` receives the identity and does nothing; stereo frame transport runs through
`transform_frame_by` instead.

`FixedFixedBirelationSet` and `VarVarBirelationSet` have no consumer outside their own property
tests. They are retained: they are the closure of the arity and factor-count axes, and a foundation
crate keeps a complete uniform surface.

`try_remap` is **not** retained on that argument, and is removed at both the relation-set and
aggregate layers. Closure under laws covers algebraic members; checked-versus-asserted is an
ergonomics pair and does not get that argument. The evidence agrees: `remap` has 34 non-test call
sites at the relation-set layer and 6 at the aggregate layer, while `try_remap` has none at either.
Its preconditions remain individually checkable through existing public predicates —
`Molecule::check_integrity`, `Correspondence::is_total`, and the entity counts — so removing it takes
away a pre-bundled check rather than a capability. As it stands it also returns a bare `Option` that
conflates a failing integrity contract with a correspondence that does not fit, which doc
[168](168-api-hygiene-2026-07-27.md) records as a failure-expression question in its own right.

Both layers keep exactly one member, and they keep the same one.

`has_incident` and `has_incident_edge` are retained. They have eighteen call sites across the six
entity views and atom-constraint evaluation, and they are a `binary_search` where `incident` is two
`partition_point` calls.

`find_by_participants` is **removed as a name**. It conflates two operations whose keys differ:
lookup, which names an entity by its constituents, and coincidence, which decides whether two entries
from two sides denote the same relation. They part on derivability.

**Coincidence stays in graph-core**, as `coincident`. Its key is the full participant multiset, which
every relation set already holds, so nothing family-specific enters. It is also required there:
`pushout` and `pullback` join on it, and an operation those two are defined in terms of cannot be an
implementation detail of a caller. `participants_match` is its mechanism and is retained. The six
entity families expose it as `coincident_id`, in each one's own vocabulary, so a family-level caller
reaches it without going through storage.

**Lookup moves to the entity families.** Its key is the family's uniqueness key — any member atom for
an aromatic system, the site for a stereo entity — which no storage shape can state, for the same
reason frame structure is not derivable from one.

Also renamed: `relation_ids` is `ids`. On a relation set the qualifier repeats the type.

Aromatic and multicenter pushout is correct today only because construction sorted both sides into
the same frame. Removing eager ordering without explicit transport would let `combine` meet
misaligned electron vectors. This is confirmed by measurement below.

## The identity rule

`Ordered` and `Unordered` encode whether storage sorts. They describe the semantics wrongly in both
directions.

An aromatic system is `Unordered`, meaning order is not the datum, but once the electron counts
moved into the payload its order became the coordinate frame of that vector. Sorting a frame whose
payload rides along is what required `on_permutation`. Stereo ligands are `Ordered`, meaning order
is the datum, but the order is not the datum either: it is the coordinate frame the coset is read
against, and any reordering carried by a matching coset transform denotes the same entity. This is
why the storage matcher ignores the marker and compares stereo ligands as a multiset.

One rule replaces the axis for all six families:

> The participant multiset is the relation's identity. The stored frame is the coordinate system its
> payload is expressed in.

Every factor is a frame. Removing the markers changes no comparison: the storage matcher already
ignored them, which was itself the contradiction between what the marker claimed and what matching
did.

The rule also settles why electron counts are payload rather than part of the participant key.
Identity must be invariant under resolution. `ElectronCountsForm::Undetermined` is filled in later
by resolution and by format raising; if counts were part of the key, an aromatic system's identity
would change when its counts were resolved.

**Stereo asymmetry to preserve.** Integrity establishes stereo uniqueness by *site alone*. Lookup
and pushout coincidence key on *site and ligand multiset*. `of_id(site, ligands)` therefore does not
find a stereo entity that exists on that site under a different ligand set, and a same-site
different-ligand collision remains two distinct entries which checked publication rejects. The
lookup key is a strict superset of the uniqueness key, for stereo only, deliberately.

## Comparison semantics

Every relation comparison is a composition of three independent choices. Today each site makes all
three by hand, per family, which is why no two of them agree. The composition is normative:
**identity is the participant multiset, frame transport is `reframe`, and the value relation is the
caller's. No site derives its own alignment.**

| operation | identity | frame transport | value relation |
| --- | --- | --- | --- |
| relation-set `PartialEq` | stored sequence equal | none | `==` |
| family lookup | the family's uniqueness key | none — the key has no frame in it | none |
| pushout coincidence | full participants | into the retained left frame | `meet` |
| `Molecule::equiv` | stored sequence equal | none | normalized equality |
| `Molecule::equiv_under(c)` | mapped through `c` | `reframe`; enumerate for stereo | normalized equality |
| `superimpose`, `difference_to` | participant multiset | rhs into the lhs frame | span or delta construction |
| `canonical_eq` | canonicalize both | selected by normalization | structural `==` |

**Lookup is not a quotient operation.** Its key is whatever integrity establishes as the family's
uniqueness key, and every one of those is order-free, so no frame enters and no reframing is
involved. That is deliberate rather than incidental, and it is why lookup sits outside the three
nested equalities entirely:

| family | uniqueness key | integrity rule |
| --- | --- | --- |
| aromatic system | its atom set, and in fact any single member atom | `AromaticSystemsOverlap` |
| multicenter bond | its atom set | `MulticenterBondsIdentical` |
| noncovalent bond | its unordered pair | `NoncovalentBondsParallel` |
| dative bond | acceptor and donor set | `DativeBondsParallel` |
| stereo atom | the site atom alone | `StereoAtomSitesDuplicate` |
| stereo bond | the site bond alone | `StereoBondSitesDuplicate` |

Ligands are therefore **not** part of stereo lookup. A site bears at most one stereo entity, so
`StereoAtomViews::of_id(site, ligands)` and its stereo-bond counterpart carry an argument that cannot
change the answer; the ligand argument is removed. Ligands could not have served as a key in any
case, because a ligand atom belongs to every stereo entity it participates in and adjacent
stereocentres routinely make each other ligands.

Coincidence is the separate operation, and it does compare full participants: two overlays on one
site with different ligand sets must stay **distinct**, so the glued molecule carries both and
checked publication reports the over-coordination. Merging them instead would ask a frame transport
to relate two different ligand sets.

The key is currently computed in three places — a sort inside the storage matcher, seven view
`of_id` methods translating each family's natural key, and the DSL namespace's `HashMap`, whose
entries it already calls the canonical participant key. The family type becomes the single owner of
both the key and the index that serves it. Index shape differs per family and is a cost question
recorded in doc [208](208-canonicalization-scaling-2026-08-24.md); graph-core's union incidence
index remains, since `has_incident` and incident iteration have their own consumers.

The last row states a relation, not an algorithm. `canonical_eq` is equality modulo entity
relabeling, frame selection, and value normalization; canonicalizing both operands is how it is
decided, not what it means.

## Layer boundary

**graph-core** owns dense relation ids, participant and payload storage, the derived incidence
index, structural equality and hashing, graph-id transport, and relation algebra under explicit
entry evidence. It never inspects a payload.

**graph-IR** owns frame selection, payload frame transport, semantic relation identity,
family-specific uniqueness, and symmetry-aware alignment.

Under this boundary graph-core has no frame operation at all. `remap` and `compact` relabel
participant ids and preserve sequence; `pushout` and `pullback` hand entries to a caller-supplied
closure.

## graph-core surface

### Id transport

`Correspondence` is currently the only one of the three id-transport concepts with a complete
layering. The two missing single-id-space layers are where the duplication sits.

| concept | single id space | graph, node and edge | molecule, eight families |
| --- | --- | --- | --- |
| partial bijection | `Correspondence<Id>` | `GraphCorrespondence` | `MoleculeCorrespondence` |
| removal and dense shift | absent | `Compaction` | `IdCompaction` |
| total relabel | absent | `Remapping` | `IdRemapping` |

`compact_relation`, `uncompact_dense`, and `normalize_removed` in `ir/remap.rs` reimplement
`Compaction::compact_node`, `Compaction::uncompact_node`, and the sort-and-dedup in
`Compaction::new` over a different id type. Filling the compaction layer deletes them.

The selected layering, regular on both axes:

```rust
pub struct Compaction<Id> { .. }        // removal-driven dense renumbering, one id space
pub struct GraphCompaction { .. }       // Compaction<NodeId> + Compaction<EdgeId>
pub struct MoleculeCompaction { .. }    // GraphCompaction + six Compaction<..Id>
```

`GraphCompaction` remains two-space because `RelationParticipant::compact` needs both spaces for a
participant that may be a node or an edge. The generic layer is extracted from inside it rather than
imposed on it. `MoleculeCompaction` holds typed per-family compactions rather than untyped
`Vec<RelationId>`, parallel to `MoleculeCorrespondence` holding typed correspondences. Every entity
id is a `u32` newtype with bidirectional `From`, so the `Compaction<RelationId>` to
`Compaction<DativeBondId>` conversion at the six family boundaries is a zero-cost typed map.

The remapping row is out of scope here and is recorded in doc
[212](212-remapping-layer-2026-08-26.md). It is the same concept at three layers, but unlike
compaction it is not duplicated code: `Remapping` is dense positional
(`self.nodes.get(old.0 as usize)`) while `IdRemapping` is eight `HashMap<Id, Id>`. Extracting
`Remapping<Id>` therefore requires a representation decision that `Compaction<Id>` does not, and
nothing in this document depends on it.

### Relation storage

The implemented storage contract is retained except where the identity rule changes it. Relation ids
are dense and follow entry order; the participant and payload vectors have equal length; incidence
is derived, deduplicated per relation and reference, and excluded from equality and hashing; indexed
access panics for an out-of-range `RelationId`; structural equality and hashing compare the stored
participant sequences and payloads; and construction does not reject two relations with the same
participant multiset. The type stays a dense collection of relation instances rather than a
participant-keyed map: graph-IR integrity establishes family-specific uniqueness, graph-core does
not.

The two clauses that change are the ordering-dependent ones. `Unordered` construction sorting
participants and transporting `D` through `RelationData::on_permutation`, and `Ordered` construction
preserving the supplied sequence, both become one rule: construction preserves the supplied
sequence and never inspects `D`.

Shown for `VarRelationSet`; the other four shapes carry the identical surface with their own
participant and arity parameters.

```rust
#[derive(Clone, Debug, Default)]
pub struct VarRelationSet<P, D> { .. }
// PartialEq, Eq, and Hash are hand-written over participants and payloads; incidence is derived
// and excluded.

impl<P: RelationParticipant, D> VarRelationSet<P, D> {
    pub fn new(entries: Vec<(Vec<P>, D)>) -> Self;
    pub fn into_entries(self) -> Vec<(Vec<P>, D)>;

    pub fn count(&self) -> usize;
    pub fn contains(&self, id: RelationId) -> bool;
    pub fn ids(&self) -> impl ExactSizeIterator<Item = RelationId>;
    pub fn participants(&self, id: RelationId) -> &[P];
    pub fn data(&self, id: RelationId) -> &D;
    pub fn data_mut(&mut self, id: RelationId) -> &mut D;

    pub fn iter(&self) -> impl ExactSizeIterator<Item = (RelationId, &[P], &D)>;
    pub fn iter_mut(&mut self) -> impl ExactSizeIterator<Item = (RelationId, &[P], &mut D)>;

    pub fn incident(&self, node: NodeId) -> &[RelationId];
    pub fn incident_edge(&self, edge: EdgeId) -> &[RelationId];
    pub fn has_incident(&self, node: NodeId) -> bool;
    pub fn has_incident_edge(&self, edge: EdgeId) -> bool;

    /// Permute one entry's participants in place. `order` must be a permutation of
    /// `0..arity`, which is validated; the participant multiset is therefore unchanged and the
    /// incidence index stays valid. Birelation shapes carry `permute_1_with` and `permute_2_with`.
    pub fn permute_with(&mut self, id: RelationId, order: &[ParticipantPosition]);

    pub fn remap(&self, remapping: &GraphRemapping) -> Self;
    pub fn compact(&self, compaction: &GraphCompaction)
        -> (Self, Compaction<RelationId>);

    pub fn pushout(
        &self,
        right: &Self,
        combine: impl FnMut((&[P], &D), (&[P], &D)) -> Option<D>,
    ) -> Option<RelationPushout<Self>>;

    pub fn pullback(
        &self,
        right: &Self,
        combine: impl FnMut((&[P], &D), (&[P], &D)) -> Option<D>,
    ) -> Option<RelationPullback<Self>>;
}
```

`iter_mut` yields immutable participants and a mutable payload. Participants own the incidence
index, so mutating them in place would desynchronise it; the only way to change a frame is to
reconstruct through `into_entries` and `new`. That states the frame-preserving contract at the type
level.

`combine` receives the borrowed form of the `(participants, data)` entry currency that `new` and
`into_entries` already use. Birelation shapes pass the three-element tuple. The object retains the
left frame, so `combine` must return a payload expressed in the left entry's frame; that is a
documented precondition rather than an inferred one. This replaces the current signature, whose
correctness depends on both sides having been sorted at construction.

`new` preserves the supplied participant sequence and treats `D` as opaque, and loses its `D: Clone`
bound, which is already unused: the constructor moves each payload out of the supplied entries and
never clones one. `Default` constructs the empty parallel storage and incidence directly and is
retained on all five shapes.

`compact` and `remap` relabel participants and preserve sequence, and must not be collapsed into one
operation. Compaction is partial and renumbers relation ids; remapping is total over the asserted
source range and preserves them. `RelationParticipant::compact`, `uncompact`, and `remap` likewise
remain three distinct graph-id operations: `uncompact` exists for editor rollback and is not the
inverse surface of relation-set `compact`.

`count`, `contains`, `ids`, `participants`, and `data` answer different direct collection
questions and are not redundant. No public relation-view layer is introduced to combine them, and no
entry-view wrapper is introduced for the entry iterators; either would add indirection without
resolving a semantic problem.

No `permute_participants_with` operation or factor-specific variant is added. Owned entry
reconstruction through `into_entries` and `new` is the transformation seam. A general consuming
entry transformation would have to be assessed against that pair and is not presumed here.

Removed, each because the redesign removes its referent:

| removed | reason |
| --- | --- |
| `FactorOrdering`, `Ordered`, `Unordered` | encode a distinction that is false in both directions |
| `RelationData`, `BiRelationData` | exist only to repair construction-time sorting |
| `RelationEquiv`, `BiRelationEquiv` | subsumed by the form's `reframe_to` composed with `normalized_eq` |
| `find_by_participants` | conflates lookup with coincidence; the lookup key is not derivable from a storage shape, the coincidence key is and stays as `coincident` |
| `FixedRelationSet::participant_permutation` and the two `MoleculeEditor` copies | three incompatible implementations of one concept |
| `ParticipantAnchor`, `RelationParticipant::anchor` | a second routing protocol beside `refs`, which already carries the information |
| `data_iter_mut` | becomes `iter_mut`, which also carries the id and participants |

### In-place participant permutation

Removing eager sorting removes the only path that could reorder a stored frame, so one must replace
it. `permute_with`, and `permute_1_with` / `permute_2_with` on the birelation shapes, permute a
single entry's participants in place.

Taking a permutation rather than a closure over the slice is what makes the operation safe: the
multiset cannot change, so each relation's participant *set* and its relation id are both preserved,
and **the incidence index stays valid without being rebuilt**. The contract and the safety condition
are the same statement.

`ParticipantPosition` is therefore retained, no longer as the callback protocol's `σ` but as this
operation's position argument. `RelationData`, `BiRelationData`, and `FactorOrdering` are still
removed.

The order is a one-line image under the convention `new[i] = old[order[i]]` — the same convention as
today's `canonicalize_positions` and as `Permutation::act` in umol-perm, so no second convention is
introduced. A stereo caller writes its degree-6 `Permutation` into that form, leaving `Permutation`
bounded and `Copy`.

The caller must build the order in a **buffer reused across entries**, one per `reframe` call. A
fresh vector per entity would reintroduce the per-entity allocation that reconstruction was rejected
for.

Graph core moves participants and nothing else; the caller transports the payload itself through the
form's own method. Two in-place steps, no allocation, and the layer boundary is unchanged.

The alternative — reconstruction through `into_entries` and `new` — is what doc
[209](209-normalization-canonical-semantics-2026-08-25.md) S2a settled on when this operation was
withdrawn, and it is expensive for a reason that only became visible once the storage was examined.
`VarRelationSet` is CSR: one `offsets` vector and one flat `participants` vector, with no per-entity
allocation. `into_entries` explodes that into one heap allocation per entity, and `new` then rebuilds
both the CSR and the incidence index — the latter provably unnecessary, since reframing cannot change
what is incident on what. That also contradicts this document's own owned-normalization intent of
mutating copy-on-write stores through `Arc::make_mut` rather than cloning into temporary vectors.

The operation withdrawn in doc 209 was a general public mutation family proposed before the relation
API had been reviewed. What is reinstated here is narrower: one entry, a validated permutation, no
payload access, and an invariant that is checked rather than promised.

Raw relation-set and aggregate `Eq`, `Ord`, and `Hash` remain representation-sensitive. Two values
whose participant sequences differ are structurally distinct before normalization. Faithful boundary
serialization preserves the supplied frame. No uniqueness or canonical-storage guarantee is
introduced.

### Private implementation

The graph-core `Incidence` representation is a coherent derived index and stays private.
`remap_factor`, coverage checking, and participant matching stay implementation helpers rather than
public concepts; their duplicated traversals consolidate under the public contracts above.

`MoleculeEditor` keeps its private `FixedSetStorage`, `VarSetStorage`, and `FixedVarSetStorage`
wrappers. Their copy-on-write materialization, append, removal, rollback extraction, and publication
roles are legitimate editor machinery. Three duplications inside them are not, and this work removes
all three: the separate participant-alignment implementations, the removal-discovery traversal
beside the compaction path, and the difference in participant-frame behaviour before and after
publication, which disappears once construction preserves frames. The wrappers remain thin storage
adapters; no public relation-set editor type is justified by the present evidence.

## graph-IR surface

### Quotient structure

Three things act on graph-IR values, and only two of them are groups.

| | value normalization | `G_frame` | `G_id` |
| --- | --- | --- | --- |
| group action | no; idempotent confluent rewriting | yes | yes |
| witness | none | permutation, or the frame pair | correspondence |
| transport by a witness | meaningless | `reframe` | `remap` |
| select a representative | `normalize` | frame selection | `canonicalize` |
| equality modulo it | normalized equality | this document | `canonical_eq` |

Value normalization cannot join the quotient shape: folding `1 + 1` to `2` is a reduction, not a
group element, so there is nothing to transport by. `Normalize` paired with a blanket equality is
the right structure for equality modulo a confluent rewriting, and its separation from the frame
and id operations is not an accident.

The other two columns are the same shape at two carriers. A group quotient has four members: act by
a supplied witness, select a representative, select a representative and expose the witness that
reaches it, and equality modulo the group. `Canonicalize` is already three of the four at the id
level; the frame level currently has none of them collected in one place. Completing the id level —
`Canonicalize` absorbing `remap`, defaulting `canonical_eq`, and deriving `equiv_under` — belongs to
doc [209](209-normalization-canonical-semantics-2026-08-25.md) S4b, which already edits that trait;
this document owns the frame level.

**Frame selection does not belong under `Normalize`.** A confluent rewriting and the selection of an
orbit representative are different operations: the rewriting has no witness, no inverse, and no
composition law, while the orbit selection has all three. Placing them under one name — whether by
carrier-relative scoping or otherwise — leaves the combined operation with neither algebra cleanly,
so nothing about it can be stated as a law. That is the same conflation as `Ordered`/`Unordered`
meaning two things, and it is what makes properties over the combined operation degenerate into
laws that cannot fail.

The frame quotient therefore needs its own operation set and its own name. Neither is settled in
this document.

The equalities are nested because the quotients are applied in order, each reading the output of
the one before it:

```text
==  subset of  modulo reduction  subset of  modulo reduction and frame  subset of  canonical_eq
```

That containment must be a consequence of the separated operations, not a property secured by
defining the middle term loosely enough to make it hold.

**Why these operations are scattered today.** The trait chain is derivational and terminates before
the aggregates. `Equiv` is a blanket impl over `Normalize`; `Molecule` does not implement
`Normalize`, so it gets no `Equiv` either, and `Molecule::equiv` (`molecule.rs:473`) is an inherent
method shadowing a trait the type does not implement, comparing field by field rather than comparing
normal forms. `equiv_under` and `remap` are inherent for the same reason. Each trait
was therefore scoped to wherever its derivation happened to work rather than to a carrier level, and
the relation set — which falls between form and aggregate — belongs to no trait at all. That is why
its frame operations became storage callbacks and free functions. Doc
[209](209-normalization-canonical-semantics-2026-08-25.md) S2b repairs the chain by implementing
`Normalize` for `Molecule`; the inherent `Molecule::equiv` must then be removed rather than
re-pointed, because it answers a different question under the same name.

### Reframe

The frame quotient. `Reframe` selects a determinate participant frame and restates the
frame-relative payload accordingly.

#### The three nested operations

Only prefixes of the pipeline are meaningful, so the design is nested rather than a free pipeline.
Frame selection reads reduced values to break ties, and id selection reads both, so no stage stands
alone.

```text
normalize     = reduce
reframe       = reduce + frame
canonicalize  = reduce + frame + id
```

| level | operation | equality | quotient |
| --- | --- | --- | --- |
| reduction | `normalize` | `normalized_eq` | value rewriting, confluent, no group |
| frame | `reframe` | `framed_eq` | `G_frame` |
| id | `canonicalize` | `canonical_eq` | `G_id` |

Each operation contains the one above it, so the containment of the equalities is a consequence of
the definitions rather than a property secured by defining a middle term loosely:

```text
==  subset of  normalized_eq  subset of  framed_eq  subset of  canonical_eq
```

`equiv` and the `Equiv` trait retire. `equiv` was a compromise that mixed the reduction and the
frame quotient under one name, which is why it means something different on a form than on a
molecule and why the properties written over it degenerate. Its replacement is `normalized_eq`,
derived from `Normalize` exactly as `Equiv` is blanket-derived today, and answering only the
reduction question.

#### Where the semantics live

Reframe semantics cannot be derived from a storage shape. A bound such as `D: Reframe<L2>` on
`FixedVarBirelationSet` decides that factor 2 bears the frame and factor 1 does not, which is entity
semantics that no storage type states: the shape admits `N1 > 1`, and the retained
`FixedFixedBirelationSet` and `VarVarBirelationSet` could have both factors frame-bearing.

Each entity family therefore gets a type that carries its own frame structure, wrapping its storage
shape on the `Graph(Arc<Csr>)` precedent and replacing the corresponding `Arc<..>` field on
`Molecule`:

| family | wraps | frame-bearing factor | site |
| --- | --- | --- | --- |
| aromatic systems | `VarRelationSet<NodeId, _>` | members | — |
| multicenter bonds | `VarRelationSet<NodeId, _>` | members | — |
| noncovalent bonds | `FixedRelationSet<NodeId, _, 2>` | the pair | — |
| dative bonds | `FixedVarBirelationSet<NodeId, 1, NodeId, _>` | donors | acceptor |
| stereo atoms | `FixedVarBirelationSet<NodeId, 1, StereoLigand, _>` | ligands | atom |
| stereo bonds | `FixedVarBirelationSet<EdgeId, 1, StereoLigand, _>` | ligands | bond |

The family type owns the quotient members — select, select-with-action, and frame-equality — because
it is the value that knows both which factor is a frame and what the payload means. graph-core owns
storage shapes and knows nothing about frames.

#### The per-form methods are inherent, not a trait

There is no trait for the form-level methods. The only site that would need generic dispatch is
`EntitySpan<T>`, and that is served by a closure-taking `try_map` over its four variants rather than
a bound; `EntitySpan` currently has only `lhs`, `rhs`, and `superimpose`, so the combinator is new
but small.

The trait's apparent benefit — forcing an explicit frame decision at compile time — was never real.
A trait forces a decision when a new form *type* appears, not when a new *field* appears on an
existing form. What catches a new field is exhaustive destructuring in the method body, and that
works with or without a trait.

Selection is form-dependent only for stereo. Everywhere else it is "sort the participants", which
the family type does without consulting the payload:

| family | selects the frame | the form's method |
| --- | --- | --- |
| aromatic, multicenter | family type sorts participants | `reframe_to(from, to)` — reindex `electrons`; `Undetermined` unchanged |
| noncovalent, dative | family type sorts participants | `reframe_to(from, to)` — returns `self`; payload is frame-invariant |
| stereo atom, bond | the form, via `CosetSpace::normalizer` under the asserted kind | `select_frame(current)` and `reframe_by(Permutation)` |

The two frame-invariant families keep a method that does nothing, deliberately, and it must
destructure exhaustively so that adding a position-indexed field fails to compile here rather than
being silently left unframed. Today's equivalent — `DativeBondForm`'s empty `on_permutation` — has
no such guard, which is why the ceremony is worth its two bodies:

```rust
impl DativeBondForm {
    /// Frame-invariant: no field is position-indexed, so a frame change carries the form unchanged.
    /// Destructured exhaustively on purpose — a new positional field must fail to compile here.
    pub fn reframe_to(self, _from: &[AtomId], _to: &[AtomId]) -> Option<Self> {
        let Self { order, constraints } = self;
        Some(Self { order, constraints })
    }
}
```

The multiset precondition is checked once by the family type rather than repeated in each form.

`reframe_by` exists exactly where the action is self-contained, which is stereo; a permutation stands
alone, whereas a target frame means nothing to a form that does not carry its source. It is
`transform_frame_by` under its settled name.

The action is per entry and keyed by the family's own id type rather than positionally, and it is
present only on the two stereo families — only they carry frame-relative content outside the payload
and therefore need the action to escape to molecule-level constraints.

#### Carriers that hold more than one form per frame

`Deltas` and `ReactionSpan` both hold frame-bearing values whose shape differs from a molecule's,
and neither bends the nesting.

**`Deltas`** runs reduce and reframe but never canonicalize: its deltas *reference* entity ids
without owning them, as a relation set references participants. Its reduction is frame-insensitive,
so `reduce` then `frame` is well defined here. `fold_group` keys on entity id, `fold_created` seeds
from `Add { atoms, attributes }`, and `EntityOp::Remove { .. }` ignores its payload entirely, so an
`Add`/`Remove` cancellation never compares frames. The sorting in `Delta::remap` was therefore never
serving the fold; it was making structural equality of `Deltas` frame-insensitive, exactly the role
storage sorting played for relation sets, and it gets the same treatment in S5a.

**`ReactionSpan`** stores `EntitySpan<Form>` against a *single* participant list, so a `Modified`
span carries two forms in one frame. The family type carries every side through one frame change
using `EntitySpan::try_map`, declining if any side declines. `reframe_with_action`
returns a single action, which is correct because it is genuinely shared. Doc 209's requirement that
one action reach every carried side stops being a special rule and becomes a consequence of the
representation.

Selection then needs one rule, stated once for every carrier:

> The frame action is selected in the **intersection of the admissible groups over all carried
> sides**.

| carried sides | intersection |
| --- | --- |
| none kinded | the full symmetric group on the frame |
| one kinded, one undetermined | the kinded side's parent group |
| both kinded, same kind | that kind's parent group |
| both kinded, different kinds | rejected by integrity; see below |

With one carried side this degenerates to that side's own group, so `Molecule` and `Deltas` use the
same rule as `ReactionSpan` rather than a simplification of it.

The last row is a real hazard, not a hypothetical. `Axial` is a stereo-**atom** kind and its parent
group is restricted, while `Tetrahedral` at the same degree is unrestricted, so a `Modified`
stereo-atom span could carry sides wanting different normalizers. Doc
[209](209-normalization-canonical-semantics-2026-08-25.md) adds the integrity rules that exclude it.

#### The ambiguity boundary

Where participants are distinct — every family except stereo — integrity makes the frame change
unique and `reframe_to` is total on multiset-equal inputs.

Where participants repeat there is no single restatement, so the two stereo forms do not carry
`reframe_to` at all. `reframe_by` is the whole of their frame action, and every stereo caller
derives its own candidates through `Permutation::between_all`. Resolving the ambiguity by selecting
a normalizing action is `reframe` itself, through `CosetSpace::normalizer` and generator-based
residual-invariance checking. Resolving it by search is the following section.

The four non-stereo forms keep a `reframe_to` method, including the two whose payload is
frame-invariant, and each body destructures exhaustively. That is what forces an explicit decision
when a position-sensitive field is added — not a trait, which would only force one when a new form
*type* appeared. For stereo the same guard is `reframe_by`, which destructures
`Self { configuration, constraints }`. It is the one property of `RelationData` worth carrying
forward, and today's `DativeBondForm::on_permutation` does not have it.

The family-level members are retained even where one has no current consumer, because the members of
an algebraically closed set constrain one another through laws; omitting one leaves a hole in the law
set. That is not the situation of a speculative feature, which has no logical connection to the rest
of the surface.

#### A stored coset denotes its orbit

A frame's repeated ligands generate its residual stabilizer. The participants determine that
stabilizer, and the stabilizer together with one coset index determines the whole orbit, so a stored
`Lit` is a representative and storing only it loses nothing. The storage shape is unchanged.

> A stereo entry's configuration denotes the orbit of the stored coset under the residual stabilizer
> of its stored frame.

Orbits partition the coset space, so for determinate cosets the meet of two entries has two
outcomes: the same orbit, where either representative may be kept, or different orbits, which is
bottom. No set-valued result is required, and which representative survives does not matter. With no
repeated ligand the stabilizer is trivial, every orbit is a singleton, and this is today's
comparison exactly.

**The case it describes.** A tetrahedral site bearing F, Cl and two implicit hydrogens. The two
hydrogens are equal `StereoLigand` values, so the transposition of their frame positions lies in the
stabilizer. The tetrahedral coset space has index two and an odd permutation exchanges its cosets,
so that one stabilizer element carries `Lit(0)` to `Lit(1)`. Both records describe the same achiral
molecule. The site is prochiral and its hydrogens are enantiotopic: R/S cannot be stated without
choosing which hydrogen is which, and the ligand labels do not carry that choice.

A frame repeat is narrower than prochirality. The two methyls of an isopropyl group are also
enantiotopic, but they are distinct atoms and therefore distinct `StereoLigand` values. Repeats
arise only for `ImplicitHydrogen` and `LonePair`, which carry the site's own atom id, so two of one
kind on one site are literally equal.

**Two defects this closes, both reachable on legal data.**

`Permutation::between` declines whenever a frame repeats a ligand, including when the two frames are
identical, so the stereo `reframe_to` declines there too. Two stereo entries that coincide and agree
are then rejected: `glue`'s combine returns `None` and `pushout` fails the whole molecule glue
rather than that one entry, and the reaction application arms return
`ApplyError::StereoFrameMismatch`.

Behind it, and masked by it, the meet is itself wrong. `coset_meet(Lit(0), Lit(1), kind)` intersects
`{0}` with `{1}`, and `canon_coset` returns `Contradiction` on the empty set. Two records related by
a stabilizer element are bottom to each other. This becomes observable the moment transport is made
total by picking one witness, which is why making transport total is not the repair.

**Two approaches not taken.**

Forbidding repeated virtual ligands by integrity — at most one `ImplicitHydrogen` and at most one
`LonePair` per stereo atom, and per endpoint for a stereo bond — would be exactly sufficient, since
atom-kind ligands are already required to be distinct and virtual ligands carry the site's own atom
id. It reverses a settled decision rather than adding a rule. Doc
[103](103-stereochemistry-overlay-and-ports-2026-05-28.md) has distinctness "neither required nor
asserted", and the frame-selection section above states that a repeated ligand frame is not an error
state. It would also leave `select_frame`'s orbit-representative machinery without a purpose, since
a trivial stabilizer reduces selection to a sort, and it would require an explicit hydrogen before
any prochiral stereo assertion.

Canonicalizing each side to its orbit representative and comparing the representatives is not a
quotient. On an `[F, Cl, H, H]` frame `Lit(0)` minimizes to `Lit(0)`, while `LitSet({0, 1})` is
already stabilizer-invariant and minimizes to itself. The two values denote the same fact and the
representative does not identify them. Least-representative and orbit-closure agree on singletons
and part company on sets.

**The uniform shape.** Every entry-level stereo operation searches the candidate actions and keeps
the first that satisfies its own relation:

```rust
Permutation::between_all(source_frame, target_frame)
    .into_iter()
    .find_map(|action| /* restate under `action`, then the site's relation */)
```

Taking the first success is sound because the successful actions are a subset of one coset of the
stabilizer. Any two of them therefore differ by a stabilizer element, and so do the results they
produce, which under the rule above denote the same arrangement. That independence comes from the
orbit reading of a stored value rather than from the relation, so it holds for whichever relation a
site asks.

Searching with `equiv` supplies orbit equality without introducing a second relation. A host holding
`Lit(1)` against a delta asserting `Lit(0)` over a frame with two equal hydrogens fails under the
identity and succeeds under the transposition, which is the answer the rule requires. `equiv` is
already what every delta old-state check in `transact.rs` and every `difference_to` field comparison
uses, so no site changes the relation it asks.

| site | relation searched | what rides on the winning action |
| --- | --- | --- |
| `StereoAtoms::glue`, `StereoBonds::glue` | `meet` returns `Some` | the met value |
| `StereoAtomDelta::Remove` and its bond twin | `equiv` against the host entry | the restated attributes |
| `StereoAtomDelta::ModifyField` and its bond twin | `equiv` of `old` against the host entry | `new.apply(action)` |
| `Molecule::difference_to`, `ReactionSpan::superimpose` | `equiv` of the rhs against the lhs entry | the transported rhs form |

`Molecule::equiv_under` and the editor's `stereo_atom_equiv` and `stereo_bond_equiv` already have
this shape. The `ModifyField` arms need one action shared by `old` and `new`, and searching supplies
it: the action is the one under which `old` agreed.

`ModifyField` is the one site where the action itself escapes the search rather than only its
result, because a second value rides on it. Returning it together with the restated `old` keeps that
value from being recomputed:

```rust
let (action, reframed_old) = Permutation::between_all(&before, &after)
    .into_iter()
    .find_map(|action| {
        let restated = old.apply(action)?;
        restated.equiv(host_configuration).then_some((action, restated))
    })
    .ok_or(/* .. */)?;
*old = reframed_old;
*new = new.apply(action).ok_or(/* .. */)?;
```

Searching also moves what a failure means. `ApplyError::StereoFrameMismatch` currently reports that
two frames could not be aligned; after the change an exhausted candidate set reports that the rule's
old state is not the one the host holds, which is what `TransactionError::OldStateMismatch` already
names one phase later. Whether the two become one error, and in which phase it is raised, is settled
in S4d.2 rather than here.

**Where the closure cannot reach.** The constraint-side meets, `TetrahedralStereo` in
`constraint/atom.rs` and `CisTransStereo` in `constraint/bond.rs`, hold no participants. They cannot
compute the stabilizer and cannot close under it. A frame-blind `meet` on a form remains the finer
relation, consistent with `==` ⊆ `normalized_eq` ⊆ `framed_eq` being a nesting rather than three
spellings of one thing.

A `Permutation::between_one` returning a single witness was written and then removed: one witness is
precisely what the search must not take, since the result would depend on which one the scan
reached. `Permutation::between` is retained, for the sites where a unique relabelling is the question
rather than an incidental requirement.

#### Laws

Within the frame quotient:

```text
reframe(reframe(x)) == reframe(x)
reframe_with_action(x) == (y, a)  =>  reframe_by(normalize(x), a) == y
framed_eq(x, y)  <=>  reframe(x) == reframe(y)
framed_eq(x, reframe_to(x, f, g))            for every admissible f, g
reframe_to(x, f, f) == x
reframe_by(reframe_by(x, a), b) == reframe_by(x, a . b)
reframe_by(reframe_by(x, a), inverse(a)) == x
```

Between adjacent levels, each law relating a level to the one below it rather than to the raw input:

```text
normalized_eq(x, y)  =>  framed_eq(x, y)  =>  canonical_eq(x, y)
reframe(x).remap(c) == canonicalize(x)       where c is canonicalize's own witness
```

The second is the correct form of the transport law. Doc 209 currently states it two ways: its
settled-boundary section has `x.remap(c) == canonicalize(x)`, which is false whenever `x` is not
already reduced and reframed, while S5a states the correct version against a normalized source. The
discrepancy is a direct product of the previous partial nesting.

### Frame selection

Frame selection is `reframe` and belongs here. Doc
[209](209-normalization-canonical-semantics-2026-08-25.md) owns what an *aggregate* does with the
selected action: applying it to molecule-level constraints referring to the entity, coordinating one
action across both sides of an entity span, and the aggregate normal-form laws. The boundary is
entries against aggregates — this document owns selecting, transporting, and comparing a relation
entry's frame; doc 209 owns carrying the exposed action through the rest of an aggregate.

Doc 209's S2b currently specifies `Normalize for Molecule` as selecting participant frames *and*
normalizing attributes. Under the nesting above, `Normalize` is the reduction only and frame
selection is `reframe`, so that subitem is respecified there.

## Operational audit

Absorbed from doc 210 and reconciled against current code. Every site that combines or compares two
relation entries must make the frame relationship explicit.

| operation | current handling | required |
| --- | --- | --- |
| construction, compaction, id remapping | sorts `Unordered` factors and transports the payload through the callback | preserve the supplied sequence; leave the payload untouched |
| `Molecule::equiv_under` | `participant_permutation` plus `RelationEquiv::equiv_under` for four families; `Permutation::between_all` filtered by `transform_frame_by` for stereo | `Reframe` for the four; the existing enumerate-and-filter retained for stereo |
| `Molecule::meet_pushout` | `glue_var_overlays` relies on both sides being sorted; `stereo_glue_entries` pre-aligns the right side in a separate find-and-rebuild pass | one entry-passing `pushout` per family, reframing the right entry into the retained left frame inside `combine` |
| editor equality and alignment | two private `participant_permutation` copies whose first-position search can reuse a query position and need not be a permutation | `Reframe` on the one shared path |
| editor removal | `birelation_removed`, `var_relation_removed`, `fixed_relation_removed` rediscover removed ids in a second traversal | consume `Compaction<RelationId>` returned by `compact` |
| `ReactionSpan::superimpose`, `Molecule::difference_to` | remaps rhs participant ids but places matched rhs forms in the lhs frame without transforming them | explicit `Reframe` of the rhs entry into the lhs frame |
| substructure matching | maps the pattern ligand frame and asks the host for the coset in that frame | retain the explicit overlay alignment without a payload callback |
| reaction application of stereo configuration transport | `Permutation::between` in the two `ModifyField` arms and `reframe_to` in the two `Remove` arms, both declining on a repeated ligand frame | search the candidate actions and keep the first satisfying `equiv` against the host entry |
| reaction application of frame-relative constraint changes | not covered for every constraint kind | out of scope here; flagged for doc [204](204-reaction-application-redesign-2026-08-19.md) |
| molecule-level stereo constraints under pushout | id-remapped only, not frame-transported | doc 209, which owns aggregate constraint transport |

The editor, superposition, difference, and pushout entries in this table are corrections, not only
migrations. They are sites this work touches regardless, and touching them correctly is the fix.

## Migration evidence

The behavioural blast radius of frame-preserving storage was measured, not estimated.
`Unordered::canonicalize_positions` was made an identity — no type changes — and every suite was run.

| crate and target | passed | failed |
| --- | --- | --- |
| `umol-graph-core` lib | 681 | 14 |
| `umol-graph-core` property and integration | 170 | 0 |
| `umol-graph-ir` lib | 6170 | 14 |
| `umol-graph-ir` canonicalization | 14 | 1 |
| `umol-graph-ir` property | 324 | 4 |
| `umol-graph-ir` other integration | 24 | 0 |
| `umol-graph`, including resolution, kekulization, fingerprint | 1725 | 0 |
| `umol-io`, including SMILES, MOL, SDF parsing | 15989 | 0 |

Thirty-three failures in roughly 25,130 tests, in three groups and no others.

**Nineteen assert the removed behaviour by name** and change with it: graph-core
`test_unordered_canonicalize::{case_2_reversed, case_3_shuffled}`, `test_*_participants_sorted`,
`test_*_into_entries::case_1_canonical_entries`, `test_*_remap`; graph-IR
`test_remap_delta::{case_03_dative_resort, case_04_aromatic_resort_permute, case_05_aromatic_remove,
case_06_multicenter_resort_permute, case_07_noncovalent_resort}`.

**Thirteen are canonicalization frame selection**, repaired by doc 209:
`test_canonicalize_constitution_family_minimum`, `test_canonicalize_constitution_participant_order`,
`test_kindless_stereo_atom_frame_order`, `test_kindless_stereo_bond_frame_order`,
`test_minimum_kinded_stereo_frames` cases 1 to 4,
`reaction_span::test_reaction_span_canonicalize::case_2_constitution`, and the property tests
`reaction::canonicalize::test_reaction_canonical_eq_by`,
`reaction::canonicalize::test_reaction_canonical_hash`,
`reaction::span::canonicalize::test_reaction_span_canonical_hash`, and
`reaction::span::canonicalize::test_reaction_span_canonicalize`.

**One is the predicted correctness dependency**, `test_molecule_meet_pushout_overlays`, confirming
that aromatic pushout is correct today only because both sides were sorted. Explicit transport in S5
repairs it.

Resolution, kekulization, SMILES, MOL, and SDF parsing are unaffected, so the change does not reach
the I/O layer.

## Absorbed and superseded scope

Doc [210](210-relation-frame-storage-2026-08-25.md) is superseded by this document. Its storage
semantics, operational audit, structural-equality consequences, and evidence boundary are absorbed
above. Its `equiv_under` name collision closes when `RelationEquiv` is removed, so
`Molecule::equiv_under` is retained unchanged. Its molecule-level stereo constraint transport
belongs to doc 209 and its reaction-application item to doc 204.

Doc 209 keeps its scope and its own plan. Its S3a is revised from a blocked reconstruction
prerequisite to a dependency on S5b of this document.

## Open item

**Naming.** `Compaction<Id>`, `GraphCompaction`, `MoleculeCompaction`, `Reframe`, `reframe`, and
`reframe_by` are proposed and are subject to the nomenclature guide. Nothing else in this document
introduces a public name.

## Staged implementation plan

Every stage ends green except S5, which ends with the thirteen enumerated canonicalization failures
and is closed by doc 209. That exception is deliberate and its exact contents are listed above.

### S0 — Establish the evidence boundary

The property suite constrains the operations this work leaves alone and is vacuous on the ones it
changes. That must be corrected before any surface changes, not delivered alongside them.

Two measurements establish the position. `test_molecule_equiv_under_identity_reduces_to_equiv`
(`tests/property/molecule/comparison.rs:127`) asserts exactly the law most at risk,
`equiv_under(identity) == equiv`, but builds its second operand as a clone with one atom's charge
changed. The two participant frames are therefore identical by construction, and the law cannot
fail on a frame difference. `test_molecule_equiv_agrees_with_equality_for_normalized_molecules`
at `:120` holds trivially while storage sorts every frame.

`stereo_frame_permutation_strategy` exists and is correct, generating permutations inside the
kind's parent group. Its consumers are `stereo/semantics.rs`, `molecule/meet_pushout.rs`, and
`reaction/application.rs`. No comparison or canonicalization property uses it. No aromatic or
multicenter frame-variation strategy exists at all, because eager sorting makes one unwritable
until S5b.

#### S0a — Reach comparison and canonicalization with reframed pairs **Done**

**Module:** `umol-graph-ir/tests/property/strategies.rs`,
`tests/property/molecule/comparison.rs`, and `tests/property/molecule/canonicalize.rs`.

Add a strategy producing a molecule together with an admissible stereo reframing of it, built from
the existing `stereo_frame_permutation_strategy` so the reframing stays inside the kind's parent
group and carries configuration and frame-relative constraints. Extend the comparison and
canonicalization properties to draw from it.

**Done.** `stereo_reframed_molecule_pair_strategy` in `tests/property/strategies.rs` yields a
tetrahedral stereo atom over four element-distinguishable ligands and its reframing under a
nonidentity parent-group action, with the configuration and constraints carried through
`transform_frame_by`. `test_molecule_equiv_under_reframed` and
`test_molecule_canonicalize_reframed` consume it.

The clone-and-perturb operand in `test_molecule_equiv_under_identity_reduces_to_equiv` was
**kept rather than replaced**: it covers the frame-identical domain, which the reframed strategy
does not, and the property-test guide treats distinct operational domains as distinct evidence.
Adding the second domain achieves what this subitem needs — the law becomes capable of failing —
without discarding the first.

The aromatic and multicenter half of this generator is not writable here. It arrives with S5b,
whose evidence requirements already state it.

**Tests and evidence:** The strategy must guarantee a nonidentity frame action and a nonuniform
position-sensitive payload; a generator that can emit only the identity action or uniform payloads
does not establish the boundary. Assert the reframing is admissible by construction.

**Change class:** additive evidence; expected to turn existing laws red (green is not the success
condition for this subitem).

**Dependencies:** [dep: none]

#### S0b — Classify and record the resulting failures **Done**

**Module:** this document.

Run the graph-IR unit and property suites under the strengthened generators and classify every
failure as one of: the law is wrong, the current implementation is wrong, or the difference is an
intended consequence of this work. Record the list here, as the frame-preserving measurement above
is recorded. An unexplained failure stops the plan rather than being carried into S1.

The expected finding is that `equiv` and `equiv_under` diverge on an admissibly reframed stereo
pair, since `Molecule::equiv` compares stored ligand frames directly (`molecule.rs:540`, `:551`)
while `Molecule::equiv_under` documents that "stereo ligand frames may differ by any admissible
permutation". That divergence exists today and is not introduced by this work.

**Result.** Under the strengthened generators the graph-IR suites report **one** failure, and it is
the new property.

| target | passed | failed |
| --- | --- | --- |
| `umol-graph-ir` lib | 6184 | 0 |
| `tests/canonicalization.rs` | 15 | 0 |
| `tests/property.rs` | 329 | **1** |
| `frag_macro`, `mol_macro`, `mol_macro_ui`, `reaction_span` | 24 | 0 |

`molecule::comparison::test_molecule_equiv_under_reframed` — **the current implementation is wrong**,
in the sense that the two operations disagree about a case both document. Reduced to a minimal
example, an admissibly reframed tetrahedral pair gives `equiv_under(identity) == true` and
`equiv == false`. `Molecule::equiv_under` documents that "stereo ligand frames may differ by any
admissible permutation"; `Molecule::equiv` compares stored ligand frames directly (`molecule.rs:540`,
`:551`). The divergence is pre-existing and is repaired by doc 209 S2c, which makes `equiv` compare
normal forms.

`molecule::canonicalize::test_molecule_canonicalize_reframed` **passes**, which is the informative
negative: `canonicalize`, `canonical_eq`, and `canonical_hash` already agree across an admissible
stereo reframing. The defect is confined to `equiv`, not general to the comparison surface.

No unexplained failure arose, so the plan continues to S1.

This evidence covers stereo only. The aromatic and multicenter half of the generator remains
unwritable until S5b removes eager sorting, and S5b's evidence requirements state it.

**Tests and evidence:** The classification above, with each failure named.

**Change class:** verification and record only.

**Dependencies:** [dep: S0a]

### S1 — Extract the single-id-space compaction layer

#### S1a — Add `Compaction<Id>` and re-express the graph compaction **Done**

**Module:** `umol-graph-core/src/graph.rs` and its unit tests.

Add generic `Compaction<Id>` carrying the removal list, dense forward shift, and reverse lookup.
Re-express the existing two-space type as `GraphCompaction` holding `Compaction<NodeId>` and
`Compaction<EdgeId>`, preserving `compact_node`, `compact_edge`, `uncompact_node`, `uncompact_edge`,
`compact_node_vec`, and `compact_edge_vec` behaviour exactly. Migrate graph-core callers, including
`RelationParticipant` and the five relation shapes.

**Tests and evidence:** Table tests for `Compaction<Id>` covering removal at the front, middle, and
end, a removed id, an id beyond the range, empty removals, and the forward-then-reverse roundtrip.
Retain every existing graph compaction assertion against `GraphCompaction`.

`Compaction<Id>` bounds `Id: Copy + Ord + Add<usize, Output = Id> + Sub<usize, Output = Id>` and
stores `Vec<Id>`, so no new trait is introduced: the shift arithmetic uses `std::ops`, and `NodeId`
and `EdgeId` gained those two impls. `GraphCompaction::new` now takes typed id vectors rather than
`Vec<u32>`. That widens the id types' public surface — `NodeId + 5` becomes legal — which is the
arithmetic `compact_node` performed by hand before.

**Change class:** breaking rename with an additive generic (red until S1b).

**Dependencies:** [dep: S0b]

#### S1b — Migrate graph-IR to `GraphCompaction` **Done**

**Module:** `umol-graph-ir/src/ir/remap.rs`, `molecule/editor.rs`, `ligand.rs`, and the six
constraint modules.

Rename every `Compaction` reference to `GraphCompaction`. No semantic change.

**Tests and evidence:** Existing graph-IR unit and property suites pass unchanged.

**Change class:** caller migration (restores green).

**Dependencies:** [dep: S1a]

### S2 — Type the molecule-level compaction

#### S2a — Replace `IdCompaction` with `MoleculeCompaction` **Done**

**Module:** `umol-graph-ir/src/ir/remap.rs`, `molecule/editor.rs`, and their unit tests.

Hold a `GraphCompaction` plus six typed `Compaction<..Id>` in place of six `Vec<RelationId>`. Delete
`compact_relation`, `uncompact_dense`, and `normalize_removed`. Preserve `UndoCompaction` behaviour
over the new representation. Add the zero-cost `Compaction<RelationId>` to `Compaction<..Id>` typed
map used at the six family boundaries.

**Tests and evidence:** Retain every existing compaction and rollback assertion. Add cases where a
removed entity of one family does not shift another family's ids, and where graph removal cascades
into a relation removal.

`MoleculeCompaction::new` and `relations` take typed per-family `Vec<..Id>` rather than
`Vec<RelationId>`, so the twelve `transact.rs` call sites and the six editor removal primitives lost
the `.into()` down-conversion they previously needed. `Compaction::new` performs the sort and dedup
that `normalize_removed` did. The `Compaction<RelationId>` to `Compaction<..Id>` map is not added
here: nothing produces a `Compaction<RelationId>` until S4b returns one from relation-set `compact`,
so the typed conversion currently happens at the editor's removal site instead. `Compaction<Id>`
required `Add<usize>` and `Sub<usize>` on the eight entity ids, added to `define_id!`.

**Change class:** breaking type replacement with caller migration (red then green within the stage).

**Dependencies:** [dep: S1b]

### S3 — Establish the entity families and the reframe operations

#### S3a — Introduce the six entity-family types **Done**

**Module:** `umol-graph-ir/src/ir/molecule.rs`, the six family modules, and their unit tests.

Add one type per entity family, each wrapping its storage shape as `Graph` wraps `Arc<Csr>`, and
replace `Molecule`'s six `Arc<..>` fields with them. Each type states which factor bears the frame
and which, if any, is a site, and speaks only graph-IR ids.

**Tests and evidence:** Assert delegation for count, indexed access, incidence, and lookup on every
family, and that `Molecule`'s existing entity-family behaviour is unchanged.

Named by the plurality rule in the nomenclature guide: a bare trailing `s` marks the container, and
the guide already forbids replacing it with `Store`, `Set`, or another suffix. The six types are
added to that entry as its example.

The surface is uniform across the six and carries no storage vocabulary. Every family has `new`,
`count`, `contains`, `ids`, `attributes`, `attributes_mut`, `incident_ids`, `has_incident`, and
`into_entries` public, and `attributes_iter_mut`, `remap`, `into_arc`, and `glue` crate-visible.
The factor accessors are per-family and match the DSL vocabulary rather than a shared invented one:

| family | factor accessors |
| --- | --- |
| `AromaticSystems` | `atoms` |
| `MulticenterBonds` | `atoms` |
| `NoncovalentBonds` | `atoms`, returning `[AtomId; 2]` by value |
| `DativeBonds` | `acceptor`, `donors` |
| `StereoAtoms` | `site`, `ligands` |
| `StereoBonds` | `site`, `ligands`, plus `incident_bond_ids` and `has_incident_bond` for the edge site |

Three families keep one crate-visible escape hatch — `AromaticSystems::atom_nodes`,
`MulticenterBonds::atom_nodes`, `DativeBonds::acceptor_node` and `donor_nodes` — so a borrowed view
field stays a borrow rather than becoming an owned conversion. `NoncovalentBonds` and the two stereo
families need none: the first returns its pair by value, and `StereoLigand` is already a graph-IR
type.

`RelationId`, `NodeId`, and `EdgeId` do not appear in any signature. `data` is `attributes`,
`relation_ids` is `ids`, and `participants_1`/`participants_2` are gone. Two methods remain
crate-visible and carry a removal note naming the subitem that retires them:
`find_by_participants`, which the lookup relocation replaces with `of_id` keyed on the family's
uniqueness key, and `participant_permutation`, which S4d replaces with `reframe_to`. Neither stereo
family carries `participant_permutation`: `equiv_under` enumerates ligand frames through
`Permutation::between_all` instead, so the mechanical translation left it dead and it was dropped.
`compact` is likewise absent from all six: molecule compaction runs through the editor's storage
enums, not the family types.

Two operations moved onto the family types rather than staying generic over storage, because their
genericity no longer typechecks across six distinct newtypes. `glue_var_overlays` and
`stereo_glue_entries` are deleted from `pushout.rs`, and all six families now carry the same
`glue(&self, right, remapping) -> Option<Self>`: relabel `right` into this molecule's id space, then
meet coinciding entries and carry the rest. The stereo implementations do one thing more — reframe a
coinciding right configuration onto the retained left ligand frame before the meet — which is
family-owned knowledge that the shared signature does not have to express. `Molecule::pushout`
therefore reads as six identical `glue(..)?.into_entries()` calls, replacing about thirty lines of
stereo-specific inline reconstruction. This is the consolidation S4 was expected to perform,
arriving early because the type change forced it.

`Molecule::try_from_entries` now routes all six families through their own constructors. It
previously reimplemented five of the six conversions inline, so the family constructors were
unreachable on the only construction path that mattered.

**Change class:** additive types with a field-type change on `Molecule` (green).

**Dependencies:** [dep: S2a]

#### S3b — Add the form-level reframe methods **Done**

**Module:** `umol-graph-ir/src/ir/aromatic.rs`, `multicenter.rs`, `noncovalent.rs`, `dative.rs`,
`stereo.rs`, `electrons.rs`, and their unit tests.

Add `reframe_to` to all six forms as an inherent method — no trait. `AromaticSystemForm` and
`MulticenterBondForm` reindex `electrons` by participant, leaving `Undetermined` unchanged.
`NoncovalentBondForm` and `DativeBondForm` return `self`, **destructuring exhaustively** so a future
position-indexed field fails to compile there, with a comment saying why the no-op is written the
long way. Stereo adds `select_frame(current)` and renames `transform_frame_by` to `reframe_by`.

**Tests and evidence:** Table tests over nonuniform electron vectors with a nonidentity frame change,
`Undetermined`, mismatched lengths, a frame pair that is not a reordering of one multiset, the
identity frame, and the inverse roundtrip. For stereo cover every kind, both restricted kinds, an
undetermined configuration, a frame change outside the parent group, and repeated virtual ligands
where the change is ambiguous. Do not use uniform payloads as the only evidence.

**Change class:** additive (green).

**Dependencies:** [dep: none]

**Settled surface.** Every form takes `self` by value and returns `Option<Self>`, so the
reconstruction is consuming and the exhaustive destructure is the compile-time guard doc
[210](210-relation-frame-storage-2026-08-25.md) asked for:

```rust
impl AromaticSystemForm  { pub fn reframe_to(self, from: &[AtomId], to: &[AtomId]) -> Option<Self> }
impl MulticenterBondForm { pub fn reframe_to(self, from: &[AtomId], to: &[AtomId]) -> Option<Self> }
impl NoncovalentBondForm { pub fn reframe_to(self, _from: &[AtomId], _to: &[AtomId]) -> Option<Self> }
impl DativeBondForm      { pub fn reframe_to(self, _from: &[AtomId], _to: &[AtomId]) -> Option<Self> }
impl StereoAtomForm      { pub fn reframe_to(self, from: &[StereoLigand], to: &[StereoLigand]) -> Option<Self> }
impl StereoBondForm      { pub fn reframe_to(self, from: &[StereoLigand], to: &[StereoLigand]) -> Option<Self> }
impl StereoAtomForm      { pub fn reframe_by(self, permutation: Permutation) -> Option<Self> }
impl StereoBondForm      { pub fn reframe_by(self, permutation: Permutation) -> Option<Self> }
impl ElectronCountsForm  { pub fn reframe_to(self, from: &[AtomId], to: &[AtomId]) -> Option<Self> }
```

The two electron-bearing forms delegate their one position-indexed field to
`ElectronCountsForm::reframe_to`, which derives the reordering itself rather than through
`umol_perm::Permutation`: an aromatic system may exceed the fixed permutation degree, as doc 210
recorded. Stereo keeps `Permutation::between`, whose fixed degree is right for its kinds. Both
derivations decline on a repeated participant, on a length disagreement, and when `to` is not a
reordering of `from`.

Stereo's `transform_frame` is `reframe_to` under its settled name, and `transform_frame_by` is
`reframe_by`. Both changed from `&self` to `self`; the borrowing call sites clone explicitly, which
is what the old signature did internally anyway.

None of the four non-stereo families' inline constraint forms are position-indexed —
`AromaticSystemConstraintForm::ElectronCount`, `MulticenterBondConstraintForm::ElectronCount`,
`NoncovalentBondConstraintForm::Intramolecular`, `DativeBondConstraintForm::{Aromatic,
RingMembership}` — so `constraints` carries unchanged in all four bodies. Only stereo's
frame-relative constraints move with the frame.

**`select_frame`.** `select_frame(&self, current: &[StereoLigand]) -> Option<Permutation>` returns
one action: the admissible permutation that presents the frame sorted and, among those, presents
this form least.

```rust
let sorted = match self.configuration.kind() {
    Some(kind) => kind.class_key().space().normalizer(current)?.act(current),
    None => /* plain sort, the full symmetric group */,
};
Permutation::between_all(current, &sorted)
    .into_iter()
    .filter_map(|action| Some((self.clone().reframe_by(action)?, action)))
    .min()
    .map(|(_, action)| action)
```

An unkinded configuration has no parent group, so the admissible group is the full symmetric group
and the frame is sorted outright — the intersection table's "none kinded" row.

Ligands that compare equal leave a residual stabilizer, so several admissible actions present the
frame sorted. The selected one is the **orbit representative** under that stabilizer: the action
whose reframed form is least. This is the same shape as `SymmetryCarrier::is_stereogenic`
(`symmetry.rs:394`), which identifies cosets a local symmetry cannot distinguish through
`CosetSpace::orbit_reps`. Minimizing over the whole form rather than the coset alone carries the
frame-relative constraints into the same representative; `orbit_reps` canonicalizes the coset index
only and would leave a `LigandSymmetry` or `Topicity` constraint in whichever presentation the
normalizer happened to pick.

`reframe_by` already rejects actions outside the parent group, so admissibility needs no separate
filter. The result is total: there is always at least one candidate, and no decline path.

A repeated ligand frame is **not** an error state. Doc
[103](103-stereochemistry-overlay-and-ports-2026-05-28.md) settles it: the element is an arrangement
record, ordered ligands plus coset, with distinctness "neither required nor asserted", and a stored
coset is a faithful labeled-arrangement fact, never vacuous. Stereogenicity is a derived predicate
plus an assertable constraint, not a storage policy. Selection must therefore accept such a frame
rather than decline on it, which is what the orbit representative delivers.

**Evidence.** `repeated_ligand_stereo_atom_strategy` generates a form whose frame repeats a virtual
ligand by construction, together with a parent-group action restating it, and asserts the repeat
before asserting anything else. Its law is that the canonical value does not depend on which
presentation of the arrangement selection started from, plus convergence. Unit tables cover a
tetrahedral site with two implicit hydrogens, with three lone pairs, and with all four ligands
equal; the restricted `Axial` parent; and a 1,1-disubstituted alkene under the partitioned cis/trans
parent, which is the common non-stereogenic bond case.

**Found while testing, and fixed: a reachable panic on the coset action.** `StereoKind::act`
`expect`ed `CosetSpace::reindex`, whose `None` covers both an inadmissible permutation and an
out-of-range coset index. `reframe_by`'s guard tested only the former, and through a hardcoded index
`0` rather than the configuration's own. Form constructors are permissive, so
`StereoAtomForm::new(StereoKind::Tetrahedral, 2u32)` builds a form whose coset is outside its kind's
two-element coset space, and any frame action on it panicked.

The validation belongs in the first operation that requires the property to hold, which is `act` —
it is what needs the index to be in range. So `act` returns `Option<u32>` and the chain above it
propagates:

```rust
StereoKind::act(self, index: u32, permutation: Permutation) -> Option<u32>
StereoCoset::{apply, swap, mirror}      -> Option<Self>   // via map_index over literal indices
StereoConfigurationForm::{apply, swap, mirror} -> Option<Self>   // via map_kinded
StereoAtomForm/StereoBondForm::{apply, swap, mirror} -> Option<Self>
```

`reframe_by` then propagates with `?` and its own guard keeps only what *it* requires — that the
permutation match the kind's degree and lie in the parent group, which the constraint transport
needs whether or not a coset is present. The guard now says so directly through `CosetSpace::allows`
instead of `reindex(0, ..)`; `canonicalize::stereo_frame_permutations` used the same fake-index idiom
and was changed with it.

Callers: `canonicalize` maps the decline to `Contradiction` — a coset index outside its kind's space
denotes no arrangement, which is bottom — and reaction application maps it to
`ApplyError::StereoFrameMismatch`, already the local idiom for an inapplicable frame action.

#### S3c — Add `EntitySpan::try_map` **Done**

**Module:** `umol-graph-ir/src/ir/delta.rs` and its unit tests.

Add a closure-taking fallible map over the four `EntitySpan` variants, so a span can carry every side
through one frame change without a trait bound on its payload.

**Tests and evidence:** Cover all four variants, including a `Modified` span where the closure fails
on one side only.

**Change class:** additive (green).

**Dependencies:** [dep: none]

#### S3d — Add in-place participant permutation **Done**

**Module:** `umol-graph-core/src/relation.rs` and its unit tests.

Add `permute_with` to the single-factor shapes and `permute_1_with` / `permute_2_with` to the
birelation shapes. Validate that `order` is a permutation of `0..arity` and panic otherwise; the
participant multiset is then unchanged by construction, so the incidence index is left untouched
rather than rebuilt.

**Tests and evidence:** Assert the permuted participant sequence exactly under the
`new[i] = old[order[i]]` convention, that the payload is untouched, that incidence answers identically before and after, that the identity order is a no-op,
and that a non-permutation order is rejected. Cover both factors independently on a birelation.

**Change class:** additive (green).

**Dependencies:** [dep: S2a]

#### S3e — Add the family-level reframe operations **Done**

**Module:** the six family modules and their unit tests.

Implement `reframe`, `reframe_with_action`, and `framed_eq` on each family type, as the
`Reframe` trait in `traits.rs`: `reframe_with_action` is the required member and the other two are
laws defaulted over it, with an associated `Action` — a position order for the four
distinct-participant families, a `Permutation` for the two stereo families. Reduce first, then
select: the four frame-invariant and electron-bearing families sort their frame-bearing factor
without consulting the payload; the two stereo families ask the form under the intersection rule.
Return one action per entry, keyed by the family's own id type; only the stereo families carry a
meaningful action.

Apply the selected order **in place** — `permute_with` for participants, the form's own method for
the payload — rather than through `into_entries` and `new`. The reconstruction route costs one heap
allocation per entity out of a flat CSR, plus a CSR rebuild and an incidence rebuild that reframing
cannot invalidate.

**Tests and evidence:** Assert idempotence, that `framed_eq` agrees with comparing `reframe` results,
that the returned action carries the reduced value into the selected frame, and the inverse
roundtrip. Include a stereo frame with repeated ligands, where the action is the orbit
representative rather than a unique reordering.

**Change class:** additive (green).

**Dependencies:** [dep: S3a, S3b, S3d]

#### S3f — Require a stereo kind admissible for its site type **Done**

**Module:** `umol-graph-ir/src/ir/molecule/integrity.rs`, its public error type, and construction
tests.

Integrity validates a stereo entry's ligand arity (`check_stereo_frame_arity`) and its coset index,
but never that the asserted `StereoKind` is admissible for the site it sits on. `Tetrahedral`,
`CisTrans`, `Axial`, and `SquarePlanar` all have degree 4, so the arity check passes for a stereo
bond asserting `Tetrahedral`. Add the check with a new `MoleculeIntegrityError` variant naming the
entity, the asserted kind, and the site type.

The admissible table is:

| site | kinds |
| --- | --- |
| stereo atom | `Tetrahedral`, `SquarePlanar`, `TrigonalBipyramidal`, `Octahedral`, `Axial` |
| stereo bond | `CisTrans`, `Axial` |

`Axial` is admissible on both, since axial chirality arises at an allene's central atom and about an
atropisomeric biaryl bond. `CisTrans` is bond-only and the other four are atom-only.

This rule is a property of a single molecule, so it propagates without further code:
`ReactionSpanIntegrityError::{Lhs, Rhs}` wrap `MoleculeIntegrityError` per side, and
`Reaction::check_integrity` calls `self.lhs.check_integrity()` first.

**Tests and evidence:** Cover every kind on both site types, asserting the exact error variant for
each inadmissible pairing and success for each admissible one. Include a degree-4 inadmissible case
so the new check is shown to catch what the arity check cannot. Retain all existing integrity cases.

The check runs before the arity check, so an inadmissible pairing reports the site mismatch rather
than a ligand count, and the site/kind table is matched exhaustively: a new stereo kind must decide
its site in `integrity.rs` or fail to compile.

One canonicalization fixture had to change. `para_stereo_canonicalization_molecule` gave its eight
outer stereo atoms a cycle of four distinct kinds — `Tetrahedral`, `CisTrans`, `Axial`,
`SquarePlanar` — using the kind purely as a label to make the sites distinguishable during partition
refinement. Only three kinds are admissible on an atom at degree 4, so the fourth descriptor is now
a second `Tetrahedral` coset instead. The sites stay pairwise distinguishable and
`test_structure_partition`'s round counts are unchanged.

**Change class:** strengthened representation-integrity contract (green).

**Dependencies:** [dep: none] — its prerequisite, stereo-site incidence integrity, is complete

#### S3g — Prohibit a stereo kind change within one entity **Done**

**Module:** `umol-graph-ir/src/ir/reaction_span.rs`, `umol-graph-ir/src/ir/reaction/integrity.rs`,
their public error types, and construction tests.

Unlike S3f this is a property of a *pair*, so each side can be individually valid while the pairing
is not, and it needs a check at both entry points rather than riding the per-side delegation:

- `ReactionSpan::check_integrity` — reject an `EntitySpan::Modified { lhs, rhs }` whose two sides
  assert different `StereoKind` values. A new `ReactionSpanIntegrityError` variant.
- `Reaction::check_integrity`, inside `ReactionIntegrityCheck` — reject a `ModifyField` delta that
  replaces `StereoConfigurationForm::Kinded(k1, _)` with `Kinded(k2, _)` for `k1 != k2`. A new
  `ReactionIntegrityError` variant.

`StereoConfigurationForm::Undetermined` on either side remains admissible and contributes no
restriction.

**Nothing becomes inexpressible.** Converting a stereocenter's geometry class — an sp3 center
becoming an allene axial center — is written as removal plus addition instead of as a modification.
The rule constrains the encoding within one entity, not the set of reactions the representation can
state. That encoding also matches the chemistry: the stereogenic unit itself changes rather than its
configuration, and the two entities carry different ids with neither side ending on a duplicate
site.

Together with S3f this makes frame selection unambiguous at every carrier under the intersection
rule above: the only ambiguous case is two carried sides asserting different kinds, which these two
rules exclude.

**Tests and evidence:** Cover a modified span whose sides differ in kind, one whose sides share a
kind, one with an undetermined side, and the equivalent delta cases. Assert the exact error
variants, and assert that the removal-plus-addition encoding of a kind change is accepted.

Cover the from-sides route explicitly, since that is where the rule is observable. Two molecules
differing only in the stereo kind of one site are each valid alone, so the outcome turns entirely on
the supplied correspondence: matching the two sites 1:1 asserts one entity changed geometry and
`ReactionSpan::superimpose` and `Molecule::difference_to` both decline; leaving that family
unmatched asserts two entities and superposition yields the removal and the addition. The two cases
share one fixture so the molecules are provably the same in both.

`ReactionSpan::superimpose` built its result with `from_entries`, which panics on an integrity
failure, so the rejection escaped as a panic rather than the `None` its signature promises. It now
uses `try_from_entries(..).ok()`, matching how it already declines an incompatible correspondence.
Its doc comment names the new `None` case.

**Change class:** strengthened representation-integrity contract; breaking for any caller that
encoded a kind change as a modification rather than as removal plus addition (green).

**Dependencies:** [dep: S3f]

#### S3h — Extend the family-level reframe to the span-bearing families **Done**

**Module:** `umol-graph-ir/src/ir/reaction_span.rs`, the six family modules, and their unit tests.

Six span types — `AromaticSystemSpans`, `MulticenterBondSpans`, `NoncovalentBondSpans`,
`DativeBondSpans`, `StereoAtomSpans`, `StereoBondSpans` — each wrapping the storage shape its
`Molecule` peer wraps, with `EntitySpan<Form>` as the payload, and each implementing `Reframe`.
`ReactionSpan`'s six fields become them. One action carries every carried side through S3c's
combinator, and the span declines whole if any side declines.

**The surface is duplicated, not shared.** A payload parameter on the six `Molecule` families would
have let the implementations be reused literally, but it complicates the primary carrier to serve a
downstream one. Each span type lives beside its peer so an entity family's frame structure is stated
in one file.

**The `Modified` representative.** The four distinct-participant families select by sorting without
consulting the payload, so two carried sides cannot disagree. The two stereo families ask the form:
S3f and S3g guarantee both sides assert the same kind, hence the same parent group, hence one
candidate set — but where equal ligands leave a residual stabilizer the set holds more than one
element and the sides can each prefer a different member. The selected action minimises the reframed
span as a whole, so the representative is the pair's, not either side's. This required `EntitySpan`
to derive `PartialOrd, Ord`, which its sibling `ConstraintSpan` already did.

Retyping `ReactionSpan`'s accessors moved `canonicalize.rs`, `incidence.rs`, and
`dsl/reaction_span.rs` off the storage vocabulary as well: they now read `ids`, `attributes`,
`atoms`, `site`, `ligands`, `donors`, and `acceptor` against typed family ids instead of
`RelationId`.

Split from S3e because only this path needs the intersection rule. A `Modified` span carries two
forms against a single participant list, so without S3f and S3g its two sides can assert different
stereo kinds and the rule has no unambiguous answer. S3e carries one side per entry and degenerates
to that side's own group, so it does not need them.

**Tests and evidence:** Cover a `Modified` span whose two sides need the same nonidentity action, and
one where a frame-relative constraint on a single side forces the whole span to decline.

**Change class:** additive (green).

**Dependencies:** [dep: S3c, S3e, S3f, S3g]

### S4 — Move entry comparison and composition onto the new surface

#### S4a — Add `iter` and `iter_mut` to the five relation shapes **Done**

**Module:** `umol-graph-core/src/relation.rs`, its unit tests, and `umol-graph-ir/src/ir/molecule.rs`.

Replace `data_iter_mut` with `iter` and `iter_mut` yielding the id, immutable participants, and the
payload. Migrate the six `Molecule::modify_*` callers.

**Tests and evidence:** Assert exact yielded tuples in relation-id order for all five shapes, the
`ExactSizeIterator` length, and the empty set. Retain the existing `modify_*` assertions.

`iter_mut` yields the participants immutably beside a mutable payload — the two live in different
fields, so the borrows are disjoint. Keeping participants immutable is the point: changing them
would invalidate the incidence index, and `permute_with` is the one operation allowed to leave it
intact.

The six family types keep `attributes_iter_mut` and now build it by dropping the id and participants
from `iter_mut`, so `Molecule::modify_*` is unchanged. The relation property suite's
`assert_data_iter_mut` keeps its assertions and takes the mapped payload iterator.

**Change class:** breaking replacement with caller migration (green within the stage).

**Dependencies:** [dep: S2a]

#### S4b — Return the relation compaction from `compact` **Done**

**Module:** `umol-graph-core/src/relation.rs`, its unit tests, and
`umol-graph-ir/src/ir/molecule/editor.rs`.

Return `(Self, Compaction<RelationId>)` from all five shapes. Delete `fixed_relation_removed`,
`var_relation_removed`, and `birelation_removed` and the traversal in `MoleculeEditor::remove`.

**Tests and evidence:** Assert the surviving set and the returned compaction together, including a
relation dropped because one participant was removed, a relation dropped because a second-factor
participant was removed, and an empty compaction leaving ids unchanged.

`RelationId` gained `Add<usize>` and `Sub<usize>`, which `Compaction<Id>` requires, on the
`NodeId`/`EdgeId` precedent from S1a.

Each `compact` body stopped short-circuiting through `filter_map` and `?`: that form discarded the
drop instead of recording it. They now match on the mapped participants and push the relation id
onto the removal list when any factor is gone.

The editor's three `*Storage::compact` wrappers report the same way. Their `Shared` arm forwards
what the relation set returns; their `Mutable` arm, which walks its own `Vec`, collects the ids
itself. `MoleculeEditor::remove` then reads the six compactions its own calls return, which is what
retires the separate traversal.

**Change class:** breaking return-type change with caller migration (green within the stage).

**Dependencies:** [dep: S2a]

#### S4c — Pass entries to `pushout` and `pullback` **Done**

**Module:** `umol-graph-core/src/relation.rs`, its unit tests, and
`umol-graph-ir/src/ir/molecule/pushout.rs`.

`combine` receives both entries as a tuple per side — `(&[P], &D)` on the single-factor shapes,
`(&[L1], &[L2], &D)` on the birelations. No new types: the alternative was a named borrowed view per
shape, which would have had to carry `participants_1` / `participants_2` as public field names,
putting storage vocabulary into a public type.

Each of the six family `glue` methods is now one `pushout` call that reframes the right entry into
the retained left frame inside `combine`, using the bare act member `reframe_to`; the entries are
already reduced, so the prefix member is not used here. The two stereo methods lose their pre-pass
entirely — they had been finding the coincidence with `coincident_id`, reframing, rebuilding an
entry vector, and then calling `pushout`, which looked the same coincidence up a second time.

The retained-left-frame precondition is documented: the object keeps the left entry's participants,
so `combine` must return a payload in the left frame. graph-core cannot enforce it without knowing
what a frame is, which is the thing it must not know.

The earlier text here named `glue_var_overlays` and `stereo_glue_entries`. Both were already deleted
in S3a; what S4c collapses is the double lookup inside the six `glue` methods that replaced them.

**Tests and evidence:** Cover a coincidence whose two sides carry different frames, a right-only
entry, an inadmissible meet, and the existing pushout and pullback correspondence laws. Add an
aromatic case with nonuniform electron counts supplied in different frames on the two sides.

**Settled by measurement, not argument.** Three facts were recorded as graph-core tests rather than
reasoned about, after several rounds of my asserting them wrongly:

- An `Ordered` factor's two sides reach `combine` in **different frames today**, with no dependence
  on S5b: `new` preserves the supplied order, so a coincidence found by multiset can hold two
  presentations. This is why stereo `glue` already reframes by hand outside `pushout`. An
  `Unordered` factor's sides arrive in one frame because `new` sorts them, so the same defect is
  latent for the other four families until S5b.
- Two right entries coinciding with one left entry is **rejected**, not merged: the right
  coprojection must be injective and `Correspondence` asserts it.
- Consequently `pushout` reading its output buffer and `pullback` reading the source were
  indistinguishable. Both now read the source, so the two operations agree.

**Deferred to measurement.** Three allocation costs, none with semantic content, left for
adversarial review to surface rather than pre-empted here:

- `pushout` clones every payload of `self` to seed its output vector and then rebuilds the CSR and
  incidence index through `new`, for an operation that changes no participant. The same
  reconstruction cost as S3e's.
- `glue`'s combine reframes into an owned value and then calls `meet`, which builds a third. The
  `Lattice` in-place member is `narrow_from`, which cannot serve: its default calls `meet` anyway,
  and its `bool` return collapses ⊥ into "unchanged", which `glue` must distinguish. Saving the
  allocation would need a new member consuming an owned receiver and yielding `Option<Self>`.
- `Permutation::between_all` returns a `Vec<Permutation>`, and every entry-level stereo comparison
  now calls it. Where the frame does not repeat the vector holds exactly one 24-byte element and the
  caller's `find_map` consumes it immediately, so the cost is one heap allocation per comparison over
  a search that does no branching. Its recursion already produces candidates one at a time: a variant
  handing each to a visitor as `visit` reaches it would drop the vector and let a caller stop at the
  first success without materialising the rest. The public name is not settled and does not have to
  mirror the private `visit`.

**Change class:** breaking signature change with caller migration (green within the stage).

**Dependencies:** [dep: S3b, S3c, S4a]

#### S4d — Move `Molecule::equiv_under` onto `Reframe` **Done**

**Module:** `umol-graph-ir/src/ir/molecule.rs`, `molecule/editor.rs`,
`umol-graph-core/src/relation.rs`, and their unit tests.

Replace `participant_permutation` with the bare act member `reframe_to` in `Molecule::equiv_under`
for the four distinct-participant families, retaining the stereo enumerate-and-filter path. Delete
`participant_permutation` from all five graph-core shapes and the four family wrappers.
`ParticipantAnchor` and `RelationParticipant::anchor` were deleted with S4c, where the
`coincident` / `coincident_edge` split replaced them. Lookup moves to the entity-family types keyed
by their uniqueness keys; coincidence stays in graph-core as `coincident` and is exposed by the
families as `coincident_id`, both landed with S4c.

`RelationEquiv` and `BiRelationEquiv` cannot be deleted here: the editor is their only remaining
consumer, and it moves in S4d.1.

**Tests and evidence:** Retain every `equiv_under` case. Add a correspondence under which the mapped
left frame differs from the stored right frame for each of the four families. Assert family lookup
for all six families, including that a stereo site bearing one entity is found without reference to
its ligands, and that a ligand atom shared by two adjacent stereocentres does not confuse either
lookup.

The four families' replacement is behavior-preserving and was checked as such rather than assumed.
`participant_permutation` canonicalized the mapped query and required it to equal the stored
sequence, returning the sort's σ; for these four factors canonicalization is a sort, so the
precondition is multiset equality — the same one `reframe_to` derives its bijection from. The
`is_permutation_invariant` fast path agrees with `ElectronCountsForm::Undetermined` passing through
`reframe_to` unchanged. Repeated participants would separate them, and molecule integrity forbids
those for all four.

**The rename is not decided here, and doc 210's reasoning about it stands.** `equiv_under` names
three operations on two levels: `RelationEquiv::equiv_under` and `BiRelationEquiv::equiv_under`,
blanket over every form, which reindex a payload by a position order and then compare normal forms;
and `Molecule::equiv_under`, inherent, which verifies complete equivalence under a supplied
`MoleculeCorrespondence`. The first two are frame transport composed with `equiv` — `reframe_to` plus
`equiv` replaces them exactly — and deleting them leaves one `equiv_under` in the tree.

`framed_eq` does not compete for the name: it takes no witness, where `Molecule::equiv_under` takes
one and verifies it. So the collision doc [210](210-relation-frame-storage-2026-08-25.md) recorded
does dissolve when the payload protocol goes, as 210 said.

What is left open is a different question, and it is not S4d's. `Molecule` implements neither
`Normalize` nor `Equiv`, so `Molecule::equiv` is an inherent method shadowing a trait the type does
not implement. Doc [209](209-normalization-canonical-semantics-2026-08-25.md) S2b repairs that by
implementing `Normalize` for `Molecule`, at which point the blanket impl supplies `equiv` and the
inherent one must go. Only then is it clear what the witness-taking operation should be called
relative to what `equiv` means on the same receiver. S4d proceeds without renaming.

**Change class:** breaking removal with caller migration (green within the stage).

**Dependencies:** [dep: S3e, S4c]

#### S4d.1 — Move the editor's participant alignment onto the frame surface **Done**

**Module:** `umol-graph-ir/src/ir/molecule/editor.rs`, `umol-graph-ir/src/ir/traits.rs`, and their
unit tests.

Split from S4d because it is not a mechanical replacement. The editor defines its own
`participant_permutation` on each of the three storage wrappers and it is a **different function**:

```rust
stored.iter().map(|s| query.iter().position(|q| q == s).map(ParticipantPosition)).collect()
```

It takes the first `query` position matching each stored participant with no used-tracking, so two
stored positions can map to the same query position and the result need not be a permutation at all.
`Permutation::between`, which `reframe_to` uses, marks positions used and declines. The two agree
exactly on distinct participants and diverge on a repeated frame — which for the editor means stereo
ligands, the one place repeats are legal, reached through two of the six call sites.

So the migration carries a semantic decision: whether accepting a non-permutation there was ever
correct. Keep the editor's implementation as a private function for the duration, assert on generated
inputs that the derived result agrees, and record precisely where it does not, rather than assuming
it does. The divergence is expected, not hypothetical.

Once the six call sites are migrated, `RelationEquiv` and `BiRelationEquiv` lose their last consumer
and are deleted with their re-exports.

**Found by the harness, before any migration.** `participant_permutation` did two jobs: derive the
alignment *and* establish participant identity, since it returned `None` when the participants
differed. `reframe_to` does only the first, and for `NoncovalentBondForm` and `DativeBondForm` it
does neither — their bodies read neither frame and always return `Some`, because the payload is
frame-invariant. Aromatic and multicenter keep the check by accident, through
`ElectronCountsForm::reframe_to` validating the multiset in order to derive its reordering.

S4d had already shipped that hole into `Molecule::equiv_under` for those two families: a
correspondence pairing entities whose participants disagree compared equivalent. No existing test
covered it. Both sites now check identity explicitly through `same_participants`, and
`test_molecule_equiv_under_rejects_mismatched_participants` covers the case that was missing.

The lesson for the remaining migration: **`reframe_to` is not a replacement for
`participant_permutation`** — it replaces the alignment half only, and every call site must be
checked for whether it was also relying on the identity half.

**The two stereo sites do not use `reframe_to`.** They take the search shape of the section above —
`Permutation::between_all` filtered by `reframe_by` and `equiv` — because the editor's two stereo
call sites are exactly where a repeated ligand frame is legal, and `reframe_to` declines there. The
four distinct-participant sites use `reframe_to` as described.

**And the non-uniformity is per value, not per family.** The first reading of this was that aromatic
and multicenter keep the identity check because `ElectronCountsForm::reframe_to` must match the
multiset to derive its reordering. That holds only for `Lit`: the `Undetermined` arm returns
`Some(Undetermined)` without reading either frame, exactly as the frame-invariant forms do. So the
same defect was live for those two families whenever the electron vector was undetermined, and was
found by testing the claim rather than by reasoning about it.

Identity is therefore asked for at **every** site, in all four families, through `is_coincident`.
Transport never establishes it reliably, so no caller may infer it from a `Some`. Pinned at three
levels: `test_*_form_reframe_to_does_not_establish_identity` states it at the form, where a
determinate vector rejects a non-reordering and an undetermined one accepts even disjoint frames;
`test_molecule_editor_*_equiv_undetermined_electrons` and
`test_molecule_equiv_under_rejects_mismatched_overlay_participants` cover the two consumers.

This is the answer to whether `reframe_to`'s contract should be made uniform by validating in every
impl. It should not: validation there would be per-value anyway — the `Undetermined` arm has nothing
to reorder and no reason to read the frames — so the honest split is transport in `reframe_to`,
identity in `is_coincident`, and callers asking for both.

**Recorded divergence.** The stereo differential pins both readings side by side:

| offered frame | coset | current | derived |
| --- | --- | --- | --- |
| stored | 0 | true | true |
| stored | 1 | false | false |
| transposed | 0 | **true** | **false** |
| transposed | 1 | **false** | **true** |
| multiset differs | 0 | false | false |

The two middle rows are inverted, which is what a transposition does to a tetrahedral coset. The
derived column is correct. `test_molecule_editor_stereo_atom_equiv_reordered_frame` asserts the
first of them directly and is `#[ignore]`d until the migration lands.

The differential covers input classes rather than generated inputs: the `*_equiv` methods are
crate-visible and the property target cannot reach them.

**Tests and evidence:** An editor comparison against a reordered frame for each storage wrapper, and
the differential above over generated participant lists including repeated stereo ligands. Retain
every editor and transaction assertion.

**Change class:** breaking removal with caller migration, carrying a semantics decision (green).

**Dependencies:** [dep: S4d]

#### S4d.2 — Put the remaining stereo transport sites on the search shape **Done**

**Module:** `umol-graph-ir/src/ir/stereo.rs`, `reaction.rs`, `umol-perm/src/permutation.rs`, and
their unit tests.

**Module:** also `umol-graph-ir/src/ir/substructure.rs`, `view/stereo.rs`, `traits.rs`,
`ir/reframe.rs`, and `umol-perm/src/permutation.rs`.

A workspace sweep found **fourteen** sites, not six. Beyond the six scoped here —
`StereoAtoms::glue` and `StereoBonds::glue` through `reframe_to`, the two `Remove` arms of
`reframe_stereo` through `reframe_to`, and its two `ModifyField` arms through
`Permutation::between` — the same shape sits in `Molecule::equiv_under` (already searching), the two
stereo arms of the matcher's `verify_overlays`, the editor's two `*_equiv`, and
`StereoAtomView::coset_for` with its bond twin.

They differed only in the relation and in what the caller took. That is now one operation:
`FrameAction` and `find_reframed` in `ir/reframe.rs`, over `Permutation::visit_between`. Each site
is one closure — `meet` for `glue`, `matches` for the matcher, `equiv` for the editor and the
reaction arms, and the coset projection for `coset_for`. `Molecule::equiv_under` keeps its own
enumeration, because it hands the actions to a molecule-level constraint check rather than selecting
one.

`Permutation::between_all` is renamed `enumerate_between` and joined by `visit_between`, the
`ControlFlow` primitive it now collects, on the `Graph::visit_paths` / `enumerate_paths` precedent.
Nothing collects a candidate vector any more, which resolves S4c's third deferred allocation.
`reframe_to` is deleted from the two stereo forms; the four distinct-participant families keep
theirs.

**The matcher was silently dropping stereo constraints.** Its two arms compared `kind` and then a
coset through `coset_for` and `coset_matches`, so a pattern asserting `Topicity`, `LigandSymmetry`,
`Fluxionality` or `Stereogenicity` matched a host that contradicted it — against the module's own
rule that an unevaluated construct fails loudly, which is what a molecule-scope constraint does.
Comparing whole forms closes it, and the explicit kind guard goes with it, since
`StereoConfigurationForm::meet` already rejects a kind mismatch. `coset_matches` keeps one caller,
`matches_value` on the kinded constraint types, where the kind is a macro constant and there is no
configuration to hand to `Lattice::matches`.

`permutation_for` and `permutation_for_ligands` are deleted. `coset_for` is retained: its consumer
is `umol-graph`'s stereo perception, comparing an entity's coset against a `#T` / `#C` constraint,
which carries no frame of its own and so has no form to compare against.

Two sites are left untraced. `umol-io/src/table_ir/raise.rs:290` and `:383` call
`Permutation::between(..).expect("validated .. frames contain the same ligands")`, so a repeated
frame panics there rather than declining. Whether one can reach them depends on what that upstream
validation checks. `symmetry.rs:561` also uses `between`, but is guarded: `all_distinct` at line 190
excludes a repeated frame before the branch that calls it, so the symmetry machinery declines to
reason about a prochiral centre rather than being wrong about it.

**Tests and evidence:** Glue two coincident stereo atoms over a frame with two implicit hydrogens
whose stored cosets differ by the stabilizer transposition, and assert the glue succeeds rather than
failing the whole molecule. The same for a stereo bond over a 1,1-disubstituted alkene. Apply a
reaction whose rule removes, and one whose rule modifies, a prochiral stereo centre against a host
storing the other representative. Assert that a genuinely different orbit still fails, separating
the two outcomes that `None` conflates today. Retain every existing application and pushout
assertion: with no repeated ligand the candidate set holds one action and behaviour is unchanged.

**Change class:** correctness change with caller migration (green).

**Dependencies:** [dep: S4d.1]

#### S4e — Add explicit transport to superposition, difference, and matching **Done**

**Module:** `umol-graph-ir/src/ir/reaction_span.rs`, `molecule.rs`, `substructure.rs`, and their
unit tests.

`Molecule::difference_to` delegates to `ReactionSpan::superimpose`, so there is one site. A
`Modified` span carries two forms against a single participant list — the lhs one — and the rhs form
arrived from the remapped rhs set still in its own frame, ids relabelled and sequence preserved.
Each of the four matched-pair loops now restates it into the lhs frame first: `reframe_to` for the
four distinct-participant families, `find_reframed` for the two stereo families. A frame that admits
no restatement declines the superposition rather than silently comparing across frames.

The matching clause was satisfied by S4d.2, which put both stereo arms of `verify_overlays` on the
same search.

**Tests and evidence:** `test_reaction_span_superimpose_stereo_reframed` states one stereocentre on
the two sides in transposed ligand orders with the correspondingly flipped coset — the same
arrangement — and asserts the span is `Unchanged` and the difference empty. Before the change it
recorded `Modified { Lit(0), Lit(1) }`, a stereo inversion that is not there.
`test_reaction_span_superimpose_aromatic_reframed` is its aromatic counterpart. That one passes
today, because construction still sorts an `Unordered` factor and transports the payload with it, so
both sides reach the span in one frame; it is the guard for S5b, after which only the explicit
restatement keeps them aligned.

**Change class:** correctness change with caller migration (green).

**Dependencies:** [dep: S3h, S4d]

### S5 — Make storage frame-preserving

#### S5a — Take the graph-IR sorting sites off `Unordered` **Done**

**Module:** `umol-graph-ir/src/ir/delta.rs`, `canonicalize.rs`.

Nine sites sort through `Unordered`, which S5b deletes:

| site | what it sorts |
| --- | --- |
| `delta.rs:2605`, `:2620` | `Delta::remap`, dative donors |
| `delta.rs:2647`, `:2661` | `Delta::remap`, aromatic atoms, payload permuted to match |
| `delta.rs:2688`, `:2702` | `Delta::remap`, multicenter atoms, payload permuted to match |
| `delta.rs:2731`, `:2744` | `Delta::remap`, noncovalent atoms |
| `canonicalize.rs:3576` | `sort_ligand_frame` |

Replace each with `canonicalize_positions`'s body: a stable sort of indices by value, yielding the
sorted values and the position order.

Whether `Delta::remap` should sort at all is doc
[209](209-normalization-canonical-semantics-2026-08-25.md) S3a's.

**Tests and evidence:** No behaviour changes; every existing assertion stands.

**Change class:** mechanical (green).

**Dependencies:** none.

#### S5b — Remove the ordering markers and the payload callback **Done**

**Module:** `umol-graph-core/src/relation.rs`, `lib.rs`, their unit and property tests, and every
graph-IR type annotation carrying an ordering marker.

Delete `FactorOrdering` with `canonicalize` and `canonicalize_positions`, its implementors `Ordered`
and `Unordered`, and `RelationData` / `BiRelationData` with `on_permutation` and
`is_permutation_invariant`. In graph-IR delete the six impls and `EntitySpan`'s two forwarding impls
at `delta.rs:1114`–`1137`. Retain `ParticipantPosition` as `permute_with`'s position argument. Remove the `O` parameters from all five shapes, the `D: RelationData` bounds,
and the unused `D: Clone` bound on `new`. Stop sorting in `new`.

`overlay_matches` (`substructure.rs:571`) calls `on_permutation` for real frame transport of the four
non-stereo overlays during matching. Move it to `reframe_to` first.

`ElectronCountsForm::permute` stays: `Delta::remap` still permutes through
`AromaticSystemForm::permute` and `MulticenterBondForm::permute` until doc 209 S3a.
Update the nineteen tests that assert the removed behaviour.

**Tests and evidence:** Assert that `new` preserves a supplied unsorted frame for every shape, that
`remap` and `compact` preserve sequence, and that structural equality distinguishes two frames of
the same multiset. Convert the resort cases in `test_remap_delta` to frame-preservation cases.

**Also rewrite S3e's `reframe` fixtures here.** The four `Unordered` families — aromatic,
multicenter, noncovalent, and dative donors — sort at construction today, so a freshly built family
is already in its selected frame and `reframe` is a no-op on it. Their S3e tests therefore build a
family and then reach into the wrapped storage with `permute_with` to manufacture an unselected
frame, which is fixture surgery standing in for the storage this subitem delivers. Once `new` stops
sorting, each of those four fixtures constructs the unsorted frame directly and the
`Arc::make_mut(..).permute_with(..)` line goes away. The two stereo families need no change: their
ligand factor is already `Ordered` and preserves the supplied frame.

**Change class:** breaking removal (red; the thirteen canonicalization failures remain).

**Dependencies:** [dep: S5a]

**Exit state:** The workspace compiles and every suite passes except the thirteen enumerated
canonicalization and hash tests, which require frame selection from doc 209. Doc 209 S2 resumes
here.

### S6 — Align the guides

#### S6a — Update the nomenclature and development guides **Done**

**Module:** `docs/development/nomenclature.md` and `docs/development/data-types.md`.

Add a *reframe* entry. Rewrite the relation-set entry, whose three axes still name `Unordered` and
`Ordered` as what controls canonicalization. Remove `ParticipantPosition` and `ParticipantAnchor`
from the participant and incidence entries. Rewrite the equivalence entry, which describes the
frame-aware `equiv_under` traits and permutation-invariance skipping. Restate the electron-counts
entry in frame terms. Record the compaction layering and the participant-multiset identity rule.

**Tests and evidence:** Search both guides for every removed name; each remaining occurrence must
describe another current API. Run `git diff --check`.

**Change class:** documentation migration (green).

**Dependencies:** [dep: S5b]

### Dependency summary

The critical path is
`S0a -> S0b -> S1a -> S1b -> S2a -> S3d -> S3e -> S4d -> S4e -> S5a -> S5b -> S6a`.

`S2a -> S4a -> S4c -> S4d` is a parallel branch of comparable length, and `S3f -> S3g -> S3h` gates
`S4e` alongside `S4d`. `S4d.1 -> S4d.2` extends the critical path between `S4d`
and `S4e`: S4e's stereo half uses the search shape that S4d.2 establishes. `S4b` hangs off `S2a` and
joins nothing until S5. The editor storage collapse is doc [213](213-editor-overlay-storage-2026-08-27.md).

S3b and S3c are additive and were built at any point before their consumers. S3a introduces the
family types and gates S3e. S3d supplies the in-place `permute_with` that S3e applies, which is why
it precedes S3e rather than sitting in S4 with the rest of the graph-core work: grouping graph-core
subitems by theme put a consumer before its foundation.

**No subitem depends on another document.** S3f and S3g are the two stereo-integrity rules that make
the intersection rule unambiguous; they were specified in doc
[209](209-normalization-canonical-semantics-2026-08-25.md) and moved here, because S3h cannot be
built without them and this plan has to be executable as one list. Doc 209's S2a depends on them at
their new labels.

S5b is the only subitem that ends red. That is an exit state, not a dependency: it deliberately
leaves the thirteen enumerated canonicalization and hash tests failing, and doc 209 S2 resumes from
there. Every stage is required; nothing here is deferrable.

S0 is a prerequisite for the whole sequence and not only for this document: the same laws guard doc
209's normalization and its completion of the `Canonicalize` quotient. It is deliberately expected
to end with failures rather than green, and S0b's recorded classification is what makes the rest of
the plan safe to execute.
