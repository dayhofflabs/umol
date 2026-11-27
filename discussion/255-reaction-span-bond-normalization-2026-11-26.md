# 255 — Reaction-span bond normalization

Status: Proposed
Date: 2026-11-26
Relates: [233](233-graph-ir-delta-review-2026-09-23.md),
[254](254-edge-colored-automorphisms-2026-11-26.md),
[integrity guide](../docs/development/integrity.md)

## Purpose

Localized bond identity in `umol` is defined uniquely by unordered endpoint pair `{a, b}`
([integrity.md](../docs/development/integrity.md)). As a prerequisite to making
`umol_graph_core::Graph` strictly simple in [254](254-edge-colored-automorphisms-2026-11-26.md),
`ReactionSpan` adds `ReactionSpanIntegrityError::ParallelBonds` rejecting a `Removed` bond and an
`Added` bond on the same unordered endpoint pair `{a, b}`.

That integrity check is a provisional boundary pending resolution of the broader `Reaction` /
`ReactionSpan` delta-normalization contract in [233](233-graph-ir-delta-review-2026-09-23.md).

## Scope to settle

1. **Delta-to-span conversion and normalization:** define how a `Reaction` delta sequence that
   removes a bond on `{a, b}` and adds a bond on `{a, b}` folds into `EntitySpan::Modified` (or
   `EntitySpan::Unchanged` when both normalized forms agree) when materializing or normalizing a
   `ReactionSpan`.
2. **Overlay interactions:** verify how `StereoBond` (and any bond-referencing constraints) attached
   to a bond on `{a, b}` behave when removed+added bond deltas fold into a single bond span.
3. **Integrity check disposition:** once normalization/materialization handles removed+added pairs
   consistently, decide whether `ReactionSpanIntegrityError::ParallelBonds` remains as the permanent
   union-graph uniqueness check on `ReactionSpanEntries` or is adjusted alongside the settled
   normalization contract.
