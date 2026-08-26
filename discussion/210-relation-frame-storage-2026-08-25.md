# 210 — Relation frame storage

Status: Proposed
Date: 2026-08-25
Relates: [209](209-normalization-canonical-semantics-2026-08-25.md),
[data-type guide](../docs/development/data-types.md)

## Purpose

Graph-core currently sorts participants of unordered relations during construction and transports
position-sensitive graph-IR payloads through callbacks. Aromatic systems and multicenter bonds
therefore acquire a canonical participant frame at storage time, while stereo relations preserve
their supplied ligand frame and require explicit semantic transport elsewhere.

This document scopes a separate migration toward frame-preserving relation storage and explicit
graph-IR frame transport. It follows doc 209: normalization must first define the canonical
representative that storage will cease to choose implicitly. The migration is a separate work unit,
not part of the canonicalization-level and hash correction.

No staged implementation plan is appropriate until the relation semantics, alignment carrier, and
operation boundaries below are settled.

## Boundary and dependency

Doc 209 owns complete normalization semantics for leaves, entity entries, `Molecule`, `Deltas`,
`Reaction`, and `ReactionSpan`, together with stereo-frame integrity. It must produce a normal form
that is independent of whether an aromatic or multicenter participant sequence was sorted at
construction.

This document owns:

- removal of eager aromatic and multicenter participant ordering from graph-core storage;
- removal or reduction of graph-core payload callback protocols;
- explicit participant alignment and graph-IR frame transport at comparison and composition
  boundaries;
- the consequences for raw structural equality, hashing, ordering, and faithful serialization;
- the cross-operation audit required to prevent omitted frame actions; and
- the naming distinction between aggregate correspondence verification and local frame transport.

The migration should install the graph-IR destination before removing the graph-core source.
Frame-transport code established for normalization may temporarily implement the existing callback
protocol; the later removal must not introduce a second independent permutation implementation.

## Proposed storage semantics

Every graph-core relation set preserves the supplied participant sequence and treats its payload as
opaque. Construction, compaction, and id remapping relabel participants without selecting a new
frame or interpreting relation data. Constructors do not require a normalized frame.

For aromatic systems and multicenter bonds, graph-IR normalization selects a deterministic
participant sequence and transports the electron-count vector through the same alignment. This is
the same ownership pattern already required by stereo: storage preserves a frame; graph-IR decides
which changes of frame preserve meaning.

Raw relation-set and aggregate `Eq`, `Ord`, and `Hash` remain representation-sensitive. Two values
whose participant sequences differ are structurally distinct before normalization, even when one
is a semantic reframing of the other. Normalized semantic equality and canonical hashing converge.
Faithful boundary serialization preserves the supplied frame rather than silently emitting a
graph-core-selected frame. No uniqueness or canonical-storage guarantee is introduced.

This direction could remove `RelationData`, `BiRelationData`, `on_permutation`,
`is_permutation_invariant`, and the relation-payload `RelationEquiv::equiv_under` and
`BiRelationEquiv::equiv_under` protocol from graph-core. It does not require storing atom/electron
pairs, which would duplicate relation participants and would not represent a whole-vector
undetermined value.

It remains open whether `Ordered`, `Unordered`, and `FactorOrdering` continue as relation-semantic
markers once they no longer control storage-time ordering.

## Participant alignment

Any operation that matches two stored relation entries must make their participant alignment
visible before combining their payloads. A structural match should carry both the relation id and a
right-to-left alignment; a combining operation must not receive two apparently aligned payloads
without that information. Graph-core may derive and carry the alignment but must not interpret or
transform graph-IR payloads.

The alignment has a named direction and arbitrary degree. Its defining law is:

```text
target[i] = source[alignment.source_position(i)]
```

It cannot generally be represented by `umol_perm::Permutation`, whose fixed degree is appropriate
for supported stereo kinds but not for arbitrary aromatic and multicenter systems. A validated
alignment may be converted to `Permutation` when a stereo operation requires its group action.

Graph-IR needs one exhaustive, consuming frame-transport operation per position-sensitive form. An
illustrative internal shape, not a selected public API, is:

```rust
fn reframe_by(self, alignment: &ParticipantAlignment) -> Result<Self, FrameMismatch>
```

For `AromaticSystemForm` and `MulticenterBondForm`, the implementation reconstructs the form by
naming every field, applies the alignment only to `electrons`, and carries other fields unchanged.
Adding a field then requires an explicit frame-sensitivity decision at compile time.
`ElectronCountsForm::Undetermined` is unchanged; a literal vector is reordered. Degree mismatch is
an error rather than a silent no-op.

The form operation does not update relation participants or molecule-level constraints. An
entry-level operation derives the alignment between source and target participant sequences,
transports the form, and returns the target frame and form together. Aggregate operations coordinate
external frame-relative constraints. An `EntitySpan` applies one selected alignment to every
carried side; it must not independently select an alignment for each side.

## Relation lookup and composition

`find_by_participants` already identifies relations by participant multiset rather than stored
order, including relation sets whose factors are marked `Ordered`. Deferred ordering therefore does
not require a new lookup identity. Lookup finds the relation id; the caller derives the alignment
between the two stored sequences.

Stereo pushout is a partial precedent. After remapping the right molecule into the common id space,
`stereo_glue_entries` finds a coincident left entry, restates the right configuration in the retained
left ligand frame, and then combines the forms. Aromatic and multicenter pushout can use the same
shape: find the coincident entry, transport the right electron-count vector into the retained left
frame, and combine aligned forms. Right-only entries retain their supplied frame until aggregate
normalization.

The generic relation-set pushout may remain payload-opaque if coincident payloads are pre-aligned to
the retained left frame. Alternatively, pushout could supply the participant alignment to its
combining closure. That interface decision remains open.

## Operational audit

Removing eager ordering affects every operation that constructs, remaps, compares, or combines a
participant-bearing entry. Existing stereo handling identifies the operation sites but is not a
complete implementation to copy.

| Operation | Current stereo handling | Aromatic and multicenter consequence |
| --- | --- | --- |
| Relation construction, compaction, and id remapping | Preserve ligand sequence and leave the form in that frame. | Preserve participant sequence and leave electron counts in that frame. |
| Aggregate normalization and canonicalization | Explicitly select and transport stereo frames, including molecule-level stereo constraint leaves in canonicalization. | Explicitly sort participants and transport electron counts; do not rely on graph-core reconstruction. |
| `Molecule::equiv_under` | Map participants through the molecule correspondence, enumerate admissible ligand-frame permutations, transform forms, and coordinate molecule constraints. | Derive the atom-participant alignment and transport electron counts before comparison. |
| Substructure matching | Map the pattern ligand frame and ask the host for the coset in that frame. | Retain explicit overlay alignment without graph-core payload callbacks. |
| Editor equality checks | Mutable storage preserves supplied frames, but stereo equality currently uses the invariant-payload path and does not transform a reordered frame. | Retain direct participant alignment and transport electron counts explicitly. |
| Molecule pushout | `stereo_glue_entries` aligns the right inline form to the retained left frame before `meet`; molecule-level stereo constraint leaves are currently only id-remapped. | Apply the same left-frame alignment to electron counts before `meet`. |
| `ReactionSpan::superimpose` and `Molecule::difference_to` | Remap rhs ligand ids but currently place matched rhs forms in the lhs frame without transforming them. | Replace implicit alignment through unordered reconstruction with explicit rhs-to-lhs electron transport. |
| Reaction application | Reframe stereo configuration changes and removal forms into the matched host frame; position-sensitive constraint changes are not all covered. | Reframe electron-count field changes into the matched host overlay frame. |
| Delta remapping and composition | Relabel ligand ids in place and retain their sequence. | Stop sorting aromatic and multicenter participants; preserve the frame until normalization. |

The uncovered stereo omissions are independently relevant correctness work: editor comparison,
reaction-span superposition, reaction application of frame-relative constraints, and pushout of
molecule-level stereo constraints must not become the model for other entity families. The storage
migration must establish one explicit frame-transport rule and apply it at every operation boundary.

## Evidence boundary

Frame-action properties must guarantee a nonidentity participant alignment and a nonuniform
position-sensitive payload. General molecule strategies that sort participants, construction that
erases unordered input order, atom-only correspondence scenarios, and uniform electron vectors do
not make an omitted action observable.

Small arities should enumerate every admissible permutation. Required laws include inverse frame
roundtrips, convergence under normalization, preservation of correspondence-relative equivalence
and substructure results, pushout invariance after normalization, and difference/application
reconstruction across independently framed sides. Minimized failures remain regressions attached to
the violated operation law.

Compile-time exhaustiveness and operation-local contracts are as important as property evidence.
Every position-sensitive form transport must reconstruct the form exhaustively, and every
payload-combining boundary must either receive aligned payloads under an explicit precondition or
receive the alignment itself.

## Correspondence-equivalence naming

The inherent `Molecule::equiv_under` and relation-payload
`RelationEquiv::equiv_under` operations have different meanings. The molecule operation verifies
complete equivalence under a supplied `MoleculeCorrespondence`; the payload operation compares two
forms after a local participant-frame permutation. Removing the payload protocol eliminates the
name collision but does not settle whether the molecule operation should be renamed to emphasize
correspondence verification.

## Open semantic questions

Before implementation planning, settle:

1. whether graph-core retains any semantic distinction between `Ordered` and `Unordered` factors
   after both preserve supplied participant sequences;
2. the concrete participant-alignment carrier, its direction, validation, and visibility;
3. whether relation pushout exposes alignment to the combining closure or requires pre-aligned
   payloads;
4. the exact ownership and visibility of graph-IR frame-transport operations;
5. which existing stereo omissions are prerequisites for the migration and which remain separate
   corrections; and
6. whether `Molecule::equiv_under` is retained or renamed after the payload protocol is removed.
