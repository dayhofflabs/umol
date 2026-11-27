# 255 — Reaction-span bond normalization

Status: Proposed
Date: 2026-11-26
Relates: [233](233-graph-ir-delta-review-2026-09-23.md),
[238](238-non-cascading-removal-2026-10-16.md),
[254](254-edge-colored-automorphisms-2026-11-26.md),
[integrity guide](../docs/development/integrity.md)

## Purpose

Molecule representation integrity guarantees that every localized bond and overlay is unique for
its participant key ([integrity.md](../docs/development/integrity.md)): `{a, b}` for `Bond` and
`NoncovalentBond`, `(acceptor, donor multiset)` for `DativeBond`, participant set for
`AromaticSystem` and `MulticenterBond`, and site for `StereoAtom` and `StereoBond`.

As a prerequisite to making `umol_graph_core::Graph` strictly simple in
[254](254-edge-colored-automorphisms-2026-11-26.md) (`S0b`), `ReactionSpan` adds
`ReactionSpanIntegrityError::ParallelBonds` (rejecting two union-frame bonds on the same `{a, b}`,
including `Removed` + `Added`) and `Reaction` application adds a provisional
`bond_endpoints_are_disjoint` check.

However, the same structural-identity tension exists across **all** overlays:
[238](238-non-cascading-removal-2026-10-16.md) introduced edit-handle "replacement identity"
allowing explicit `Remove` + `Add` over identical participants (for example `AromaticSystem` over
the same atom set, or same-kind stereo over the same frame) to remain unmatched in `witness` and
`ReactionSpan` and drop host-only constraints. That directly contradicts participant-key structural
identity (where only atoms have independent token identity, and bonds and non-kind-changing
overlays are uniquely identified by their participant keys — note `MoleculeCorrespondence::induce`
already distinguishes mismatched determined stereo kinds like `Tetrahedral -> Axial`).

`test_reaction_tracked_apply_at_aromatic_replacement` in `umol-graph-ir/src/ir/reaction.rs` was
removed in [254](254-edge-colored-automorphisms-2026-11-26.md) `S0b` because it asserted same-key
`Remove` + `Add` replacement identity for an internal bond and an aromatic system over identical
participants; the DPO aromatic-bond dangling checks remain covered by
`test_reaction_apply_at_aromatic_bond` and stereo kind replacement remains covered by
`test_reaction_tracked_apply_at_stereo_kind` / `test_reaction_tracked_apply_at_stereo_replacement`.

## Scope to settle

1. **Structural identity vs. same-key replacement across bonds and overlays:** reconcile
   participant-key uniqueness ([integrity.md](../docs/development/integrity.md)) with
   [238](238-non-cascading-removal-2026-10-16.md)'s replacement identity — decide whether union-frame
   same-key `Removed` + `Added` is rejected uniformly across bonds and non-kind-changing overlays,
   and whether host-only-constraint reset via same-key replacement is retired.
2. **Delta normalization and consistency ([233](233-graph-ir-delta-review-2026-09-23.md)):** settle how
   normalization/consistency disposes of same-key `Remove` + `Add` deltas (rejecting rather than
   silently coalescing into `Modify`, which would alter DPO dangling semantics) and revisit the
   provisional `bond_endpoints_are_disjoint` check in `Reaction::prepare_deltas` /
   `Reaction::apply_deltas_at`, along with same-ID `Add` + `Remove` cancellation (which can hide
   transient parallel-bond adds or leave dangling created references).
3. **Sequential composition (`Reaction::compose`):** settle the interaction in `compose_overlap`
   when reaction `A` deletes a bond on `{a, b}` (preserving `a, b`) and reaction `B` creates a bond
   on `{a, b}`: composing witnesses leaves the two bonds unmatched, so `ReactionSpan::superimpose`
   returns `None` (dropping that overlap because a single simple-graph DPO rule cannot express the
   intermediate bond-deletion dangling check).
4. **Stereo and constraint interactions:** preserve legitimate kind-changing or frame-changing stereo
   replacements (`StereoKindModified`) and verify constraint-reference behavior under the settled
   union-frame uniqueness rule.
5. **Integrity check disposition:** align `ReactionSpanIntegrityError` across `Bond` and overlays
   once the structural-identity and normalization contracts are settled.
