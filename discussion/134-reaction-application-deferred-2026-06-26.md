# 134 — Reaction-application: deferred items

Two reaction-application items were left to work out while landing the `Edit`/`Delta`/`Undo` vocabulary refactor
(`Add`/`Modify`/`Remove` ↔ `Added`/`Modified`/`Removed`, key-based `ModifyConstraint`, by-value
molecule `Add`/`Remove`, `CascadedConstraints`) and designing the reaction EDN surface (doc 133).
The vocabulary/transact groundwork is landed and green; these are the remaining application
features, both increment-2. Captured here so they are not lost; to be tackled after the reaction
DSL (133) is implemented.

## 1 — molecule-level constraints don't apply in reactions

**State.** `ReactionAst::apply_at` has `Delta::Constraint(_) => {}` — a molecule-level
`ConstraintDelta` (add/remove on the flat `Vec<Constraint>`) is **silently dropped** on apply.
(The span side is resolved separately: `ReactionSpanAst` gains a `constraints: Vec<ConstraintSpan>`
field and `to_reaction_span` populates it — doc 133, decision 14. This item is the *operational*
`apply_at` path only.)

**Mechanism is ready.** The transact layer now supports general by-value
`Edit::Add/RemoveMoleculeConstraint` (true multiset; `Remove` is last-match-by-value with the
position captured for rollback via `CascadedConstraints`). So this is *closable* — only the
lowering is missing.

**Blocker.** Applying a rule-frame constraint to the host needs its atom/bond refs remapped through
the match `m` (lhs atom → host atom, created atoms → appended). `IdRemapping` / `Constraint::remap`
only express *removal-compaction* (old→new after deleting ids); the match is an arbitrary
injection, so there is no existing remap that fits.

**Proposed solution.**
- Add a total match-based ref-map: `Constraint::map_topology_refs(self, atom: impl Fn(AtomId) ->
  AtomId, bond: impl Fn(BondId) -> BondId) -> Self`, mirroring `Constraint::remap`'s tree traversal
  (entity-leaves, `Relational`, `Molecule`-scope `:atoms`/`:bonds` lists and sub-pattern anchors,
  `And`/`Or`/`Not`). It is total (no drops). Overlay entity-leaf refs (dative/aromatic/multicenter/
  noncovalent/stereo) are out of increment-1 — the embedding maps only atoms/bonds — so those error
  (and are subsumed by item 2).
- Wire `apply_at`: for `Delta::Constraint(ConstraintDelta::Add/Remove(c))`, emit
  `Edit::Add/RemoveMoleculeConstraint { c.map_topology_refs(lhs→host via m, created→appended) }`,
  ordered **before** `RemoveTopology` so transact's own renumbering (`remap_with_update`) finishes
  the job (and drops the constraint if it referenced a removed atom).
- Test: a reaction with a molecule-level constraint delta applies and round-trips through rollback.

Smaller of the two: one remap method + the `apply_at` arm.

## 2 — overlay entity `Delta` infrastructure

**State.** The AST `Delta` is `Atom(AtomDelta) | Bond(BondDelta) | Constraint(ConstraintDelta)`
only. The six overlay entities — dative bond, aromatic system, multicenter bond, noncovalent bond,
stereo atom, stereo bond — have **no `*Delta`**. So a reaction cannot add / modify / remove an
overlay entity at all.

**Mismatch.** The reaction EDN surface (doc 133) and the `apply` / `to_reaction_span` design are
written for all eight entities, but the AST `Delta` covers three (atom / bond / constraint). This is
the increment-1 / increment-2 split in the localized-topology scoping, left to work out.

**Proposed solution.**
- Add the six overlay `*Delta` enums on the `EntityDelta` pattern: `Add` / `ModifyField` /
  `ModifyConstraint` / `Remove`, extend `Delta`, and the `apply_at` / `to_reaction_span` lowering.
- The `EntityDelta` generic path (the `diff` / `deltas_from_states` default methods + the
  `field_ops!` macro) means each new entity is mostly mechanical — per entity only: its `*Delta`
  enum, the `EntityDelta` impl (the `(variant ⇒ field)` map), the `Atoms` type, `into_delta`, the
  `MoleculeAst` field, the `Edit` variant, and its `EntitySpan` span slot.
- **Span generalization (tied to 133 work item 2).** Every `MoleculeAst` entity collection is
  parameterized by its AST payload type `R`: atoms/bonds as `Vec<R>`, the six overlays as `*Relation`
  / `*Birelation` sets that split a topological part (participants + incidence over `NodeId` /
  `EdgeId`) from a `data: Vec<R>` payload column. `ReactionSpanAst` is the identical structure with
  that payload lifted `R → EntitySpan<R>` in every collection — topology shared, only the per-entry
  value becomes a span; `EntitySpan<T>` itself is unchanged. It's done for the two vecs; lifting the
  six overlay relation sets' `data` columns to `EntitySpan<…>` is the remaining work.

Larger of the two: six overlay deltas + the span generalization.

## Sequencing

Both are increment-2, after the reaction DSL (133) lands. Item 1 is independent and small. Item 2 is
the larger one and subsumes item 1's overlay-ref limitation. Neither blocks the 133 surface
design, which is for atom/bond/constraint (localized topology) and complete on its own.
