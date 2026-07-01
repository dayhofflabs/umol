# 134 — Reaction-application: deferred items

Five items were left to work out while landing the `Edit`/`Delta`/`Undo` vocabulary refactor
(`Add`/`Modify`/`Remove` ↔ `Added`/`Modified`/`Removed`, key-based `ModifyConstraint`, by-value
molecule `Add`/`Remove`, `CascadedConstraints`) and designing the reaction EDN surface (doc 133) —
two reaction-application features, plus structural entity refs (item 3) and a constraint-store
cleanup (item 4), which surfaced while implementing 133, plus a serde streaming-parser / naming
review (item 5) that surfaced while wiring the R12 tests. The vocabulary/transact groundwork is
landed and green; all are increment-2, captured here so they are not lost, to be tackled after the
reaction DSL (133) is implemented.

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
- **Span DSL for overlays (render + parse).** `dsl/reaction_span.rs` covers atoms/bonds/constraints
  only. Adding overlay spans needs the overlay entry parsers/renderers shared the way atoms/bonds are:
  atoms/bonds reuse `parse_atom_entry`/`parse_bond_entry` and the factored `render_atom_value` /
  `render_bond_entry` (molecule.rs). The six molecule overlay renderers (`render_dative` … etc.) are
  still inline (TODO in `molecule.rs`) — factor each into a shared `render_<entity>_entry`
  parameterized by the rendered value, then the span wraps them with the `{:add|:modify|:remove}`
  verbs, mirroring `render_bond_span_entry`.
- **Remapping split + reindex (graph-core; settled 2026-06-29, design in doc 131 #7/#8).** Overlay
  deltas and lhs overlays must re-anchor through frame changes. Split `Remapping` into two disjoint
  types: `RemovalRemapping` (today's `Remapping`, monotonic compaction, `apply_removal_remapping ->
  Self`, never reindexes) and a new general `Remapping` (total relabel, `apply_remapping -> (Self,
  Vec<ParticipantPosition>)`). `ParticipantPosition(u32)` is a new graph-core newtype (positions are
  implicit slice indices today). Only `VarRelationSet` (aromatic/multicenter `electrons`) consumes
  the permutation, via `AromaticSystemAst`/`MulticenterBondAst::reindex`; dative/noncovalent are
  scalar and stereo is `Ordered` (coset reframe is the separate #8 op), so they ignore it. `apply_*`
  added uniformly across all five relation-set types. Compose's `remap_delta` HashMaps unify onto the
  general `Remapping`.

Larger of the two: six overlay deltas + the span generalization + the remapping split.

## 3 — structural entity refs

**State.** Every `<entity>-ref` is `int | keyword` (position or id) — reaction surface (doc 133) and
spec (§7.9). A bond/overlay with no id can only be named by position; to reference it otherwise you
must give it an `:id`.

**Want.** Name an entity by its constituents instead — a bond by its endpoints, an aromatic /
multicenter system by its members, a dative bond by donors+acceptor, a stereo element by
site+ligands. (Atoms are the base; no structural form.)

**Form — explicit map, uniform.** `<entity>-ref ::= int | keyword | <structural-map>`, the §4 entry
form minus `:type`/`:id`:

| entity | structural ref |
|---|---|
| bond, noncovalent-bond | `{:atoms [atom-ref atom-ref]}` |
| aromatic-system, multicenter-bond | `{:atoms [atom-ref+]}` |
| dative-bond | `{:donors [atom-ref+] :acceptor atom-ref}` |
| stereo-atom | `{:site atom-ref [:ligands [atom-ref+]]}` |
| stereo-bond | `{:site bond-ref [:ligands [...]]}` |

Map form (not the bare `[atom-ref …]` vector) so a structural ref is self-delimiting and never
collides with the vectors it nests inside (anchor pairs `[[ref ref]+]`, relational `[ref target]`).
Checked: the forms paired with a ref never key on `:atoms`/`:donors`/`:site`, and stereo forms are
vectors, so `[<ref> <form>]` stays unambiguous by shape/key.

**Sites.** One production, so it reaches every non-atom ref at once: reaction `:remove`/`:modify`;
entity constraints; relational leading refs plus the embedded `bond-ref` in `:dative-bond-parallels`
/ `:stereo-bond-site`; `:bond-order-sum :bonds`; anchor pairs; stereo-bond entry `:site`. In an
entity constraint a structural ref is the *query* form (point at an id-less entity by its parts); the
relational leaves still assert membership of a named entity.

Extending a production propagates everywhere automatically:

  - A. Reaction deltas (doc 133, not yet in spec): :remove <ref>, :modify <ref> …. ← the proposed site.
  - B. Entity constraints (§7.9 677–683): {:bond [bond-ref form]} + the 6 other non-atom entities.
  - C. Relational constraints (§7.9 686–717): the leading <entity>-ref of every leaf, plus two embedded cross-refs to a bond-ref: :dative-bond-parallels [dative-bond-ref bond-ref] and
  :stereo-bond-site [stereo-bond-ref bond-ref].
  - D. Molecule-scope (§7.9 722): :bond-order-sum {:bonds [bond-ref+] …}.
  - E. Anchor-spec (§7.9 794–800): pair-lists [[<entity>-ref <entity>-ref]+] for bonds, dative-bonds, aromatic-systems, multicenter-bonds, noncovalent-bonds, stereo-atoms, stereo-bonds.
  - F. Stereo-bond entry :site (§4 line 148): a bond-ref to an existing bond.

**Cost.** Not a `define_ref` tweak. The macro emits a uniform `enum { Index, Id }` resolving over an
`id → index` map. The structural variant has a per-entity payload (`[AtomRef;2]` / `Vec<AtomRef>` /
donors+acceptor / site+ligands) and resolves by matching the host's entity *collections* by
constituents — a different payload and resolution per entity. §4.1 uniqueness (no two bonds on a
pair, aromatic systems disjoint, …) makes the match ≤1 hit, so the semantics are clean; the code
shape is the work.

**Decided (2026-06-29): no parallel overlays of any kind — structural refs are total and
unambiguous.** Every overlay family — including **noncovalent bonds** — forbids two entries with the
same constituents/roles. This drops the one prior exception (spec §4.1 currently lets a noncovalent
pair carry two different-kind interactions); the clearer semantics outweigh that flexibility. So a
structural ref names **any** overlay (incl. noncovalent) by its constituents with ≤1 hit, and the
overlay-composition correspondence is by structure = by id (doc 131 #7 settled). **Spec §4.1 must
change:** rewrite the `:noncovalent-bonds` clause to "at most one interaction per unordered atom
pair" (drop the different-kind-coexist allowance), and add a `:multicenter-bonds` clause (no two with
the same atom set) — §4.1 currently omits multicenter.

## 4 — entity constraints: uniqueness-by-key explicit

**State.** Each entity constraint store (`AtomConstraints`, plus the bond / overlay analogs) is
addressed two ways: by `kind` (`AtomConstraintKind`), which is **multi-valued** — `RingMembership`
has one entry per size — and by `key` (`AtomConstraintKey`), which is **always unique** (a kind's
subkeys distinguish its entries). `add` calls this distinction `is_unique` ("unique kinds replace,
ring appends") — misleading: by key *everything* is unique; "non-unique" only means a kind has
subkeys.

**Want.** Redesign the constraint-store API so uniqueness-by-key is the clear primary invariant —
every entry is addressable by a unique key, a kind is just a group of keys — and retire the
kind-centric `is_unique`/"unique" nomenclature. Applies to all entity constraint stores.

**Now.** `AtomAst::update` (reaction `:modify` resolution) already works per key via `remove_by_key`;
this item is the broader API/nomenclature cleanup, not blocking.

## 5 — FromEdn streaming-parser audit + parse-direction naming review

Two related reviews surfaced while wiring the R12 tests (doc 133).

**5a — streaming `from_edn_str` audit.** Per dsl-serialization, every EDN-shaped `*Dsl` must override
`FromEdn::from_edn_str` to drive its hand-written streaming `read_*` — never the trait default, which
delegates to the tree (`read_string` + `from_edn`), violating "streaming must not delegate to tree."
Today only the top-level / full-entity types do (`MoleculeDsl`, `ReactionDsl`, `AtomDsl`, `BondDsl`,
`DativeBondDsl`, … via `read_subgrammar_all`); the constituent boundary types inherit the
tree-delegating default. Audit every `FromEdn` impl and add the override where missing. The full
type list and per-type wiring notes are in doc 133's R12 deferred note; in summary:

- Most delegate to an existing single-arg `read_<type>_dsl` (mechanical); a few readers
  (`read_topicity`, `read_stereogenicity`, the stereo constraint readers) are private and need
  `pub(super)`; the ref and stereo constraint `*Dsl` macros need the override in the macro body.
- `read_molecule_constraint_dsl(de, key)` and `read_relational_constraint_dsl(de, key)` take a
  pre-consumed dispatch key (the `{`+key is consumed by the umbrella `read_constraint_dsl`), so they
  have **no** standalone reader — a thin standalone reader or a small restructure is needed, a design
  decision rather than a mechanical change.
- `fuzz_constraints` calls `ConstraintDsl`/`ConstraintsDsl::from_edn_str` directly, so its "streaming"
  arm is currently the tree path — fixing the override makes that target exercise the real streaming
  reader.

**5b — parse-direction naming.** Review the naming across the two parse directions and decide a
consistent scheme: string → AST (`FromStr`/`parse_<type>` → AST via defaults) versus the EDN-DSL
boundary (`FromEdn::from_edn` tree / `from_edn_str` streaming, `read_<type>` / `read_edn_<type>`,
`*Dsl` vs `*Ast` trait targets). Confirm the `parse_`/`read_`/`read_edn_` split reads unambiguously
for "parse a string into an AST" vs "parse into the EDN/DSL form," and settle the prospective
`read_edn_` → `pick_` rename (dsl-serialization, not yet adopted).

## Sequencing

All five are increment-2, after the reaction DSL (133) lands. Item 1 is independent and small. Item
2 is the larger one and subsumes item 1's overlay-ref limitation. Item 3 rides with item 2 — six of
its seven structural forms are overlay entities and share the resolution rework; bond is the only
localized entity gaining one, and it is unblocked (name an existing bond by `:id` or index until
then). Item 4 is independent API/naming cleanup. Item 5 is independent serde cleanup (5a mechanical
plus the two dispatch-key readers; 5b a naming pass). None blocks the 133 surface design, which is
atom/bond/constraint with index|id refs and complete on its own.

## 6 — Increment-2 implementation plan

Dependency-ordered plan for the overlay-composition increment: items 2 + 3 above, the
overlay-composition semantics settled in doc 131 (#7/#8), and item 1's overlay-ref limitation that
item 2 subsumes. Each subitem lists its location and `[deps]`. Items 4 and 5, and umol-graph
no-parallel-overlay enforcement, are independent and out of this plan's critical path.

**By location (i).**

| location | subitems |
|---|---|
| graph-core | I0a–I0e |
| umol-ast/ast (entity ASTs) | I1a, I6b |
| umol-ast/ast/delta.rs | I1b, I1c, I1d, I6a |
| umol-ast reaction.rs (apply) | I2a, I2b, I6c |
| umol-ast reaction_span.rs | I2c–I2f, I7a |
| umol-ast compose.rs | I3a, I3b, I3c |
| umol-ast/dsl (tree / stream / format / traits) | I4a–I4e, I5a, I6d |

**Dependency-ordered items (ii–iv).**

**I0 — graph-core: remapping & participant positions** `[—]` **Done**
- I0a — graph-core: `ParticipantPosition(u32)` newtype. `[—]` **Done** (relation.rs, beside `RelationId`; re-exported from lib.rs).
- I0b — graph-core: rename `Remapping`→`RemovalRemapping`; relation-set `apply_remapping`→`apply_removal_remapping` on all five types (no reindex — monotonic); migrate callers (`MoleculeBuilder`, the `IdRemapping` wrapper); move data-column compaction to free `remove_node_vec`/`remove_edge_vec`. `[—]` **Done** (79 type + 22 method renames; `apply_to_*_vec`→free fns over `map_node`/`map_edge`; `IdRemapping`/`UndoRemapping` untouched; graph-core 296 + umol-ast 3914 tests green).
- I0c — graph-core: new general `Remapping` — total relabel (`Vec<NodeId>`/`Vec<EdgeId>`), `map_node`/`map_edge`/`unmap_*`/constructors. `[—]` **Done** (`new`/`map_node`/`map_edge` + `#[rstest]` tests; graph-core 301 green). `unmap_*` **deferred** — partial inverse for a non-bijective injection with no current consumer (forward covers apply/remap_delta/reindex; `reverse()` builds a fresh forward map; rollback uses `Edit`/`Undo`); add when a consumer appears.
- I0d — graph-core: `apply_remapping(&Remapping) -> (Self, …σ)` on all five relation-set types (argsort → per-relation participant permutation). `[I0a, I0c]` **Done.** Per-factor σ (non-conflated): single-factor sets return `(Self, Vec<ParticipantPosition>)`, birelations `(Self, Vec<ParticipantPosition>, Vec<ParticipantPosition>)`. Added `RelationParticipant::remap(&Remapping) -> Self` (renamed the removal pair → `remap_removal`/`unmap_removal`) and `FactorOrdering::canonicalize_positions` (argsort returning σ; Ordered = identity). graph-core 306 + umol-ast 3914 green.
- I0e — ~~graph-core: common-subgraph enumeration over the incidence (Levi) graph~~ **Void — folded into I3b.** No graph-core work: `enumerate_common_subgraphs` is already generic over any `Graph` + predicates, and umol-ast already builds the Levi graph with overlay pseudonodes (`MoleculeAst::incidence_graph`, `IncidenceNodeSelection::{topological,constitution,full}`, and `substructure_matches_incidence`). "Overlay-aware overlap" is a umol-ast/compose change (I3b), not graph-core. (Superseded 2026-06-30: I3b does **not** route through `incidence_graph` after all — it keeps the simple-graph overlap and derives overlay correspondence post-hoc; see I3b point (1). The "no graph-core work" conclusion stands all the more.)

**I1 — overlay reindex + the four uniform overlay deltas** `[I0]` **Done**
- I1a — umol-ast/ast: `AromaticSystemAst::permute` / `MulticenterBondAst::permute` — reorder `electrons` by a `&[ParticipantPosition]`, both delegating to the shared `ElectronCountsAst::permute` (`Undetermined` unchanged; charge/spin/constraints positionless). `[I0a]` **Done** (only the two positional families have it — no blanket needed; umol-ast 3918 green).
- I1b — delta.rs: the four non-stereo overlay `*Delta` (`DativeBondDelta`, `AromaticSystemDelta`, `MulticenterBondDelta`, `NoncovalentBondDelta`) + `EntityPatch`/`EntityFold` + extend the `Delta` sum + `Delta::inverse`. `[—]` **Done** (DSL-shaped fields — dative `donors`/`acceptor`, aromatic/multicenter `atoms: Vec`, noncovalent `atoms: [_;2]`; macros except noncovalent's hand-written patch (uninhabited `*Constraint`); added `PartialOrd,Ord` to the four `*FieldChange`s; `From<usize>` already via `define_id!`). Crate intentionally **red** until I1d/I2 (verified the only errors are the expected non-exhaustive `Delta` matches).
- I1c — delta.rs: `Deltas::canonicalize` over the four families (generic `EntityFold` fold + `fold_group`; mechanical). `[I1b]` **Done** — four id-keyed group-maps + dispatch arms + `fold_group::<…>` loops; removed-atom dangling check extended to overlay `Add` participants (dative: `donors`+`acceptor`). Closes the canonicalize non-exhaustive error.
- I1d — delta.rs: `remap_delta` over overlay deltas — re-anchor participants (atom map; aromatic/multicenter re-sort + `permute` via I1a) **and** the overlay's own relation id (bond-analogous, via per-overlay-id maps the callers build). **Reordered after I2c** (decided 2026-06-29): `reverse`/`compose` build the overlay-id maps from the span's overlay columns (I2c), so I1d depends on I2c. Note: overlay relation ids are not `NodeId`/`EdgeId`, so the general `Remapping` (I3a) won't carry them — a separate per-overlay-id map is needed. `[I1a, I1b, I0c, I2c]` **Done** — new `IdRemapping` bundle (remap.rs) carries the six total-relabel maps (atom/bond + four overlays) with `map_*` accessors; `remap_delta(delta, &IdRemapping)` relabels overlay ids and atoms, re-sorts aromatic/multicenter participants (`Unordered::canonicalize_positions`) and `permute`s electrons to match (dative/noncovalent relabel only, no positional payload). `reverse` builds all six reversed-frame maps via one `reversed_remapping` helper (replaced the inline atom/bond duplication). compose passes empty overlay maps + TODO(I3) — it has no overlay handling yet. Naming reworked (sign-off 2026-06-30): graph-core `RemovalRemapping`→`Compaction` (`remap_removal`/`unmap_removal`→`compact`/`uncompact`); ast `IdRemapping`→`IdCompaction`, `UndoRemapping`→`UndoCompaction`; `Remapping`/`IdRemapping` now name the general case. Tests — `test_remap_delta` (per-family id/atom relabel + aromatic/multicenter re-sort+permute + overlay ModifyField + Constraint passthrough) and `test_reversed_remapping` (identity/removed/created/mixed); **both pass** (umol-ast green) once the lib was bridged to link. The per-kind helper is named `reversed_remapping` (past participle — the remapping for the reversed reaction, not a verb; builds one kind's old→new map feeding `IdRemapping`; "remap" not "reindex"/"relabel" — reindex is the coset op).

**I2 — apply + span for the four overlays** `[I1]` **Done**
- I2a — reaction.rs: `apply_at` overlay lowering arms (overlay `Delta` → existing overlay `Edit`); fold in item 1 (molecule-constraint lowering via `Constraint::map_topology_refs`). `[I1b]` **§2 (overlay arms) Done** — four overlay families lowered: `ModifyField`/`ModifyConstraint`→`sets` (host ref via `m.host_*`); `Add`/`Remove` in a post-`atom_ref` second pass (dative `atoms = [donors…, acceptor]` per transact's `split_last`; noncovalent `[a,b]`), `overlay_adds` after topology adds, `overlay_removes` before `RemoveTopology`; noncovalent `ModifyConstraint` is a no-op (uninhabited `NoncovalentBondConstraint`). `apply_at` compiles. **§1 (constraint lowering) Done** (sign-off 2026-06-30): no `map_topology_refs`; instead a parallel total `Constraint::remap(&IdRemapping)` mirroring `Constraint::compact(&IdCompaction)` top-to-bottom (no inlining/coupling — each in-flow inner constraint got its own `remap`: Constraint/MoleculeConstraint/RelationalConstraint/SubPatternAnchor/the four no-op overlay singulars/stereo + `remap_atom_subset`/`remap_bond_subset`; Atom/Bond pass `c` through, mirroring compact). `IdRemapping` extended with the two stereo maps (full 8-entity counterpart of `IdCompaction`). `apply_at` builds the match `IdRemapping` (lhs→host via `m.host_*`, created atoms/bonds→appended; overlay/stereo lhs-only — a created-overlay constraint ref is unsupported) and lowers `Delta::Constraint` to `Edit::Add/RemoveMoleculeConstraint` before `RemoveTopology`. Tests: `test_constraint_remap` (7 cases) + `test_reaction_ast_apply_at_molecule_constraint` (refs re-anchor `[0,1]`→`[1,2]` through the match). umol-ast 3939 green. **I2a complete.** **Compaction-naming rename done (both crates, sign-off 2026-06-30)**: every compaction-flow symbol carries `compact_`/`uncompact_` (`Compaction::compact_node`/`compact_edge`/`uncompact_*`, `compact_node_vec`/`compact_edge_vec`, `IdCompaction::compact_*`, `UndoCompaction::uncompact_*`, `undo_compaction`, `Constraint::compact`) — disjoint from the remapping flow (`map_node`/`map_edge`, `IdRemapping::map_*`, `Constraint::remap`); graph-core green + 306 tests.
- I2b — reaction.rs: DPO dangling check extended to overlay incidence (a deleted atom's overlay participations, not only `bond_ids()`). `[I1b]` **Done** — four `removed_host_<overlay>` sets collected from the overlay `Remove` deltas (`m.host_*(id)`) in the delta scan; the gluing loop now also rejects a deleted host atom carrying a `dative_bond_ids`/`aromatic_system_id`/`multicenter_bond_ids`/`noncovalent_bond_ids` membership not in the corresponding removed set. **Stereo deferred to I6c** (no stereo deltas yet → removal inexpressible; `AtomView` exposes only the stereo *site*, not ligand incidence). Tests: `test_reaction_ast_apply_at_overlay_dangling` (noncovalent kept → `Dangling`) + `test_reaction_ast_apply_at_overlay_removed_not_dangling` (rule removes the bond → ok). umol-ast 3941 green.
- I2c — reaction_span.rs: lift the six overlay relation-set `data` columns to `EntitySpan<…>` (all six incl. stereo, uniform). `[I1b]` **Done** — six `EntitySpan<D>` overlay relation-set fields on `ReactionSpanAst` (no `Arc`, unlike `MoleculeAst` — the span isn't shared) + `from_parts` gains the six params. Confirmed compiles (the generic relation sets take `D = EntitySpan<…>` directly — no new machinery). `from_parts` now 10-arg → callers (`to_reaction_span` [I2d], `ReactionSpanDsl` `FromAst`/`IntoAst` [I4], tests) are red until updated; overlay accessors added when consumed (I2e/I2f).
- I2d — reaction_span.rs: `to_reaction_span` folds overlay deltas onto `lhs` overlays (an `apply_*_change` per family). `[I2c, I1b]` **Done** — four `apply_*_change` helpers in delta.rs (dative/aromatic/multicenter/noncovalent, mirroring `apply_atom_change`); `to_reaction_span` accumulates removed/added/changed per overlay id, folds each `lhs` overlay (Removed/Modified/Unchanged) and appends created ones, participants mapped to the union frame via `atom_index`. Stereo (both): no deltas yet (I6) → all `Unchanged` from `lhs`, site/ligand ids identity-mapped. All ten `from_parts` columns now built. `to_reaction_span` compiles; residual reds are I4 (DSL/test `from_parts` callers), I1d (`remap_delta`), I2a (`apply_at`).
- I2e — reaction_span.rs: `left()`/`right()` carry unchanged overlays through `from_parts`. `[I2c]` **Done** — `project` refactored to take a `Side` (Left/Right) enum + generic `entity_side` helper (closures couldn't be generic over the six overlay payloads); builds `compacted`/`compacted_bonds` maps then projects all six overlay columns into the `from_parts` vecs (pick side value via `entity_side`, compact participants/ligands; an overlay is dropped if its side value is absent or any participant is dropped). The dropped overlays' union ids per family now feed the constraints' `IdCompaction` (was empty lists — a latent miscompaction of overlay-referencing molecule constraints, fixed). Tests: `..._project_overlay_unchanged` (carried on both sides) + `..._project_overlay_added` (right only). umol-ast 3943 green.
- I2f — reaction_span.rs: `to_reaction` (span→deltas) for overlays (`EntityDelta::diff`/`deltas_from_states` defaults). `[I2c, I1b]` **Done** — `to_reaction` extends the atom/bond `deltas_from_states` pattern to the four non-stereo overlays: per family, the relation set's `EntitySpan<…>` `data` column is collected into a `Vec` (no slice accessor on the relation sets) and fed to `<Delta>::deltas_from_states` with a participant closure reading the union-frame ids (`participants`/`participants_1`+`_2`) — dative `(donors, acceptor)` from `participants_2`/`participants_1[0]`, aromatic/multicenter `Vec` from `participants`, noncovalent `[a,b]`. `Added`/`Removed` → structural `Add`/`Remove`, `Modified` → the field/constraint `diff` default; all reuse the existing infrastructure (no new fold code). Stereo skipped (no `StereoAtomDelta`/`EntityFold` yet — I6a; its span columns are all `Unchanged`). Removed the stale "overlays dropped on conversion" module TODO. Tests: six overlay roundtrip cases folded into `test_reaction_span_ast_to_reaction` (`dative_add`/`aromatic_add`/`multicenter_add`/`noncovalent_add`/`_remove`/`_modify`) covering each family's participant closure + Add/Remove/Modify states; umol-ast 3949 green.

**I3 — composition for overlays** `[I1, I2, I0e]`
- I3a — compose.rs: unify `remap_delta`/frame algebra onto the general `Remapping`; extend the four-class composite frame to overlay relation ids. `[I1d, I0c, I0d]` **Scope settled 2026-06-30 (two forks).** (1) *Vehicle* — keep the `IdRemapping` HashMap bundle as `remap_delta`'s input; the "unify onto general `Remapping`" clause is superseded by I1d (graph-core's dense `Remapping` can't carry overlay relation ids, and its `apply_remapping` over relation sets is only needed for the lhs-overlay carry, which is out of scope here). No `remap_delta` signature change. (2) *Boundary* — **created overlays only**: I3a fills composite-frame classes (3) (A-created) + (4) (B-created) overlay ids — disjoint by construction — and guards by bailing the whole composition when **either reactant lhs carries any overlay** (`a.lhs.has_overlays() || b.lhs.has_overlays()`). The guard is on the lhs, not the deltas: `lhs_c` is built from `raw_graph()` only, so *any* lhs overlay — modified, removed, **or untouched** — would be silently dropped; the lhs check subsumes a delta-scan (you can't modify/remove an overlay absent from the lhs) and created overlays live in the deltas so they still compose. Classes (1) (lhs_A) + (2) (L_B context), `lhs_C` overlay materialization, and the overlap-region overlay correspondence move to I3b. (A first cut guarded only overlay `ModifyField`/`ModifyConstraint`/`Remove` deltas — unsound: an untouched lhs overlay passed and was dropped. Corrected to the lhs guard 2026-06-30, with a `untouched_overlay` test case.)
- I3b — compose.rs: overlay correspondence + lhs-overlay carry, on the **simple-graph** overlap. `[I1, I2, I0d]`
  **Approach settled 2026-06-30 (two points).**
  - **(1) No Levi switch — simple-graph enumeration + post-hoc overlay correspondence.** The existing `enumerate_common_subgraphs` over `raw_graph()` (atoms as nodes, bonds via `edge_match`) already yields the atom/bond overlap; overlays do **not** go through `incidence_graph` (the Levi/pseudonode form). Rationale: `common_subgraph` enumeration is *unconstrained* — it has **no connected variant**, deliberately, precisely so overlay bridges (B₂H₆ multicenter bonds joining the two halves, H-bonded base pairs) never need to be structural edges. With no connectivity gate, covalently-disjoint matches are first-class, so an overlay never carries a match — overlays are pure attributes whose correspondence is a deterministic function of the atom map, derivable post-hoc. The Levi form would only inflate the association graph (`|V|`: `n`→`n+m`) on the Bron–Kerbosch bottleneck for zero gain (the gain it was meant to give — overlay-only connectivity — is exactly the connected-MCS feature this module omits). This supersedes the I0e/old-I3b "call enumeration over `incidence_graph(constitution())`" framing.
  - **(2) Correspondence is well-defined because parallel overlays are excluded.** Per overlap, R_A overlay `X` ↔ L_B overlay `Y` iff `participants(X)` maps onto `participants(Y)` under the atom map and kinds meet. This is a partial *function* (at most one `Y` per `X`) because the no-parallel-overlay invariant (spec §4.1) forbids **complete typed participant overlap**: two same-typed overlays may share *some* participants (partial overlap is fine) but never the *whole* participant set. The Levi form would need the identical uniqueness (else a pseudonode pairs ambiguously), so this is a shared precondition, not a tradeoff.
  
  Sub-pieces (incl. those deferred from I3a):
  - (a) **overlay correspondence (post-hoc)** — per point (2), from each overlap's atom map.
  - (b) **lhs-overlay carry** — materialize `lhs_C`'s overlays (widen the `from_atoms_and_bonds` build to `from_parts`): class (1) lhs_A overlays (identity atom relabel) and class (2) L_B *context* overlays (participants via `db_atom`), relabeled through graph-core `Remapping` + `apply_remapping` (the σ reorders aromatic/multicenter electrons).
  - (c) **overlap-region overlays** — a B lhs-overlay whose participants *all* lie in the R_A∩L_B overlap maps (via (a)) onto its matching A-product overlay, reusing that composite overlay id (no fresh one). A *partial*-overlap B overlay (some participants in the overlap, some in context) is a context overlay → class (2) (b), not a correspondence.
  - (d) **drop the I3a guard** — once (a)–(c) exist, replace the `a.lhs.has_overlays() || b.lhs.has_overlays()` bail so reactants carrying overlays compose. Overlay admissibility (R_A must supply an overlay B requires on the overlap) is I3c.

  **Impl plan settled 2026-06-30.** Verified: `from_parts` → relation `new` **preserves relation order** (canonicalizes participants within a relation, never reorders relations) ⇒ a composite overlay's relation id = its push position, so the four-class order is just push order. But `new` sorts each relation's participants **without permuting its data** (`relation.rs:407–411`) ⇒ a non-monotonic atom relabel (class (2) `db_atom` scatter) breaks aromatic/multicenter electron alignment, so positional class-(2) overlays need the σ-permute (`apply_remapping`/`permute`) — class (1) is identity (exempt), dative/noncovalent are positionless. Mechanism:
  - **Four-class overlay frame, per family**, built **inside** the overlap loop ((2) and offsets are overlap-dependent — I3a's overlap-independent created maps move in). Composite id order (1) lhs_A · (2) L_B context · (3) A-created · (4) B-created; (1)(2) live in `lhs_c`, (3)(4) are delta-created. `da_<F>`: lhs_A id → (1) index, A-created delta id → `|(1)|+|(2)|+rank`. `db_<F>`: L_B context id → `|(1)|+rank`, overlap-region id → its corresponding A composite id, B-created → `|(1)|+|(2)|+|(3)|+rank`.
  - **Correspondence is post-hoc in composite-atom space (span-free):** build the A-side-in-R_A index `{lhs_A F-overlays, id ∉ da-Remove} ∪ {da-Add}` keyed `(kind, sorted composite participants) → composite id`; for each L_B F-overlay map participants via `db_atom` — all-in-overlap ⇒ look up (hit = correspondence, miss = B needs an overlay R_A lacks ⇒ skip overlap), else ⇒ context/class (2).
  - **Guard (decision):** drop the broad `lhs.has_overlays()` bail in I3b, keep a **narrow interim guard** bailing only compositions that delete an overlay-bearing atom (the dangling risk); removed when I3c lands.
  - **RcAnchored (decision):** A's overlay `ModifyField`/`ModifyConstraint`/`Remove` deltas **extend the reaction center** (their participant atoms join `rc_a`), so RcAnchored keeps overlaps chaining an overlay-only A edit.

  **Done 2026-06-30 (DAMN; stereo is I6).** A generic `place_overlays<I, K, E>(a_lhs, a_removed, a_created, lb, b_created, relabel) -> Option<(da, db, lc)>` owns the four-class id assignment + correspondence + carry; it is called once per overlay kind in DAMN order. The key `K` and entry `E` are opaque and `relabel: FnMut(E) -> E` carries a class-2 entry into the composite frame — so dative's acceptor and noncovalent's pair are just entry shapes, aromatic/multicenter fold the electron `permute` into their `relabel`, and the design absorbs stereo at I6 (its bond-anchored key carries a bond + atoms, its `relabel` is the coset reindex). Correspondence is post-hoc on the sorted composite participant key (`§4.1` makes it a function). `lhs_c` widened to `from_parts`; the broad guard replaced by a stereo bail + a narrow `deletes_overlay_atom` (delete-an-overlay-bearing-atom) interim guard. RcAnchored extended: **any** overlay delta adds its participant atoms to `rc_a` (`Add`/`Remove` carry their atoms; the `Modify` variants name only the overlay, so its atoms are looked up in `lhs_A`) — `Add` anchors too (decided 2026-06-30: an overlay-creation is a chain just like a modify/remove). Tests: `carry` / `remove_carried` / `correspondence` / `required_absent` (skip) / `aromatic_carry` (positional) / `rc_modify` / `rc_remove` / `rc_add` / `rc_aromatic_remove`. umol-ast 3959 green, clippy clean. These are point cases; the combinatorial space needs a property test — see **I3-prop** below.
- I3c — compose.rs: admissibility (boundary-bond / combined-frame dangling) extended over overlay incidence. `[I3a, I2b]` **Done 2026-06-30** — the `db_removed_atoms` combined-frame dangling loop now also checks, for each B-deleted shared atom, that every R_A overlay incident to it (all participants in the overlap, B removes the correspondent) is co-deleted, matched by sorted L_B participant set (`b_removed_<kind>`); else the overlap is skipped. The blunt interim `deletes_overlay_atom` guard is gone — compose keeps only the stereo bail (I6) and the (now overlay-aware) combined-frame check, consistent with how it already assumes each input reaction is individually well-formed for bonds. A dedicated compose case that *trips* overlay dangling is deferred to **I3-prop** (the property test exercises it).

**Tier-2 DPO validator — Done 2026-06-30** `[I3c discussion]`. Prompted by the "A's own validity" observation: dangling-freedom is a **tier-2 DPO invariant** (model-free graph invariant), so it belongs in an external, *requested* validator, not enforced at construction — a dangling reaction is constructible, tier-2-invalid until validated (lazy, as for molecule tier-2/3). Added `umol-ast` `ast::validate::dpo::{DpoValidator, DpoContradiction, DpoError}`, same struct+`validate`+`Solution` pattern as the tier-1 validators. Two entry points (repetition accepted to avoid a conversion): `validate_reaction(&ReactionAst)` reads each deleted atom's lhs incidences against the delta Remove sets (no span build); `validate_reaction_span(&ReactionSpanAst)` scans a `Removed` atom for a surviving incidence (added `pub(crate)` overlay accessors on `ReactionSpanAst`). Bonds + the four DAMN overlays; stereo deferred to I6. Compose's during-check stays (it's the *generator* — Kekulization analogy: emit only valid reactions; lazy validation governs checking, not generation). Two distinct "dangling" notions clarified: rule-intrinsic dangling-freedom = tier-2 invariant (this validator, and ≡ compose's combined-frame check on the built composite); the classic DPO *dangling condition* at `apply_at` (I2b) is host-relative/match-time and genuinely separate. **Architectural note (open):** tier-2 validators are model-free and could all move umol-graph→umol-ast (crate boundary = model-independence), and model-free checks (tiers 1–2) could become methods *on* the AST while model-carrying tier-3 stays external — deferred; the molecule tier-1/2 relocation gets its own pass + discussion doc.

**I3-prop — compose property test** `[I3a, I3b, I3c]` **Done**

**Done 2026-07-01** — in `umol-ast/tests/property/reaction.rs` (feature `proptest`), on a new `overlay_reaction_strategy` (DAMN-overlay lhs via `from_parts`; `build_reaction` extended to co-remove a deleted atom's incident overlays, keeping the generated reactions DPO-valid). Properties landed and green: `compose_sound_overlay` (P1 soundness), `compose_dangling_free` (P3 dangling, via `DpoValidator`), `compose_rc_anchored_subset` (P2), `apply_reproduces_right_overlay` (probe: plain overlay reaction `apply == right()`), `span_roundtrip` (reaction↔span fidelity), and `compose_well_formed_overlay` (every composite's `apply == right()`). The last was `#[ignore]`d on the first cut (2026-06-30); it surfaced the three apply/remap bugs listed below, and is enabled and green now they are fixed. Coverage completed 2026-07-01 for the cheap properties: `compose_determinism` (P4), `compose_canonical_deltas` (P3), `compose_distinct_overlays` (P6), and deterministic `disjoint_sum` / `disjoint_rc_anchored` cases in `test_reaction_ast_compose` (P5). See the per-property coverage line under **Properties** below. **P1 completeness (`Full`: seq ⊆ composed) is the one deferred property**, tracked in doc 135 (needs monomorphism overlaps + `meet` interface + seam delta-rebasing); not blocking.

**Overlap enumeration split + compose rewire — done; P1 completeness itself deferred to doc 135 (not blocking).** Root of P1 incompleteness: compose enumerated R_A ∩ L_B overlaps with Bron–Kerbosch *maximal*-clique enumeration, which drops every non-maximal partial overlap and the empty one — yet each distinct overlap is a distinct composite. graph-core now separates the two tasks (mirroring the MCS split): `maximal_common_subgraphs` + `MaximalCommonSubgraphAlgorithm::BronKerbosch` (former behaviour, renamed) and a new complete `enumerate_common_subgraphs` + `CommonSubgraphEnumerationAlgorithm::Backtracking` (every clique of the modular product; shared `modular_product` + `subgraphs_from_cliques`; materialized). Compose rewired to the complete one; `test_reaction_ast_compose` split — constructive cases assert *containment* of their specific composite, `test_reaction_ast_compose_disjoint` keeps the disjoint cases exact. All compose unit + property tests green. **P1 completeness (`seq ⊆ composed`) stays `#[ignore]`d** (`compose_complete_overlay`): closing it needs three coordinated compose changes (monomorphism overlaps, `meet` interface, seam delta-rebasing) that are genuine rule-composition design work, not I3-prop follow-through — **tracked in doc 135**. I3 is otherwise complete and this does not block the rest of 134.

**Bugs the property test surfaced:**
1. **Overlay pushout-complement (fixed).** A context overlay whose overlap participant is A-created was placed into `lhs_c` referencing a non-`lhs_c` (created) atom → out-of-bounds in `apply`. Added the overlay analog of the context-bond `is_ra_created` admissibility in `compose_all` (skip the overlap).
2. **Composite `apply() != right()` → root cause: an `apply_at` id-staleness bug, not compose. Fixed 2026-07-01.** Some composites yielded **no** `apply` product. Diagnosis: **not compose-specific.** The minimal reproducer is a plain hand-written reaction removing ≥2 overlays of the *same kind* (`ReactionAst::new` + two `AromaticSystemDelta::Remove`, no compose) — `apply` errored `Transaction(IdOutOfRange("aromatic system"))` (or `OldStateMismatch` for embeddings where the shifted slot stays in range). `apply` was silently embedding-order-dependent: of N self-matches, only those mapping the removed overlays into descending removal order succeeded (this is why an early probe read `products = 1` — a lucky embedding). Compose merely makes same-kind multi-removes routine. Full analysis + fix in **Edit algebra** below.
3. **Constraint edits emitted after overlay removals. Fixed 2026-07-01.** `Constraint` carries overlay ids; `apply_at` extended `constraint_edits` after `overlay_removes` with original-space ids, so a constraint referencing a surviving overlay whose lower-id sibling was removed got a stale id. Fix: emit `constraint_edits` before all removals (analysis under **Edit algebra**). Untested previously (reaction constraints over overlays are rare).
4. **`remap_delta` didn't canonicalize unordered participants for dative/noncovalent. Fixed 2026-07-01.** After the two fixes above, `compose_well_formed_overlay` shrank to a further case: composing a plain molecule with `remove-atom + its dative + its noncovalent`. `remap_delta` re-canonicalizes participant order for aromatic/multicenter (both `Unordered`, via `canonicalize_positions` because they also permute the electron ast) but the **dative (donors)** and **noncovalent (atoms)** arms — also `Unordered` — only mapped the ids. So a relabel like `[0,1]→[2,0]` left the delta at `[2,0]` while storage canonicalizes the unordered set to `[0,2]`, and `apply`'s order-sensitive `OldState` check mismatched. Fix: call `FactorOrdering::canonicalize(&mut …)` on donors / atoms after remap (no per-participant ast to permute, so no `canonicalize_positions`); a no-op for `Ordered`, so it is safe to apply uniformly. Dative-donors was the same latent bug, masked in the minimal case by a single donor.

All four fixed; `compose_well_formed_overlay` is enabled and green (property suite 101/101, umol-ast lib 3970/3970, workspace green, clippy clean).

**Edit algebra (apply_at lowering) — id-space renumbering analysis, 2026-07-01**

`apply_at` lowers a reaction to a `Vec<Edit>` and `transact`s it. Every edit carries **pre-resolved ids** (resolved once, against the pre-apply state, via `m.host_*`). Correctness reduces to one question: *which edits renumber which id space, and was each edit's ids resolved before or after such a renumbering.* An edit is wrong iff an **earlier** edit renumbered an id space it references.

Renumbering table (verified against `builder.rs` / `transact.rs`). "append" = grows at the end, existing ids unchanged (non-renumbering); "compact" = dense swap-remove + truncate → renumbers every id above the removed one.

| Edit | Atom | Bond | its overlay type X | other overlay types | constraints |
|---|---|---|---|---|---|
| `AddAtoms` / `AddBonds` | append | append | — | — | — |
| `Add<Overlay X>` | — | — | append | — | — |
| `Modify*` (any) | — | — | — | — | — |
| `Remove<Overlay X>` | — | — | **compact** | — | rewrites X-refs |
| `RemoveTopology` | **compact** | **compact** | **compact** (cascade) | **compact** (cascade) | rewrites all refs |
| `Add/RemoveMoleculeConstraint` | — | — | — | — | insert/remove one |

Two facts drive everything: `Remove<X>` compacts **only X's own Vec** (`remove_aromatic_systems` touches `self.aromatic_systems` + constraint refs, nothing else); `RemoveTopology` compacts **every** space, because `self.remove` cascades `apply_compaction` over every overlay set (drops incident, remaps participants → dense).

Verdicts on the four interaction groups:

- **(iv) different overlay types — fully commute.** `Remove aromatic` never touches the multicenter Vec. All cross-type overlay ops (add/remove) are independent; batching is *per type*, never across types. Correct today.
- **(iii) overlay vs topology — one ordering constraint, already met.** `RemoveTopology` renumbers all overlay spaces, so every overlay-by-id edit must precede it; overlay *adds* reference atom ids, so `AddAtoms` must precede them. `apply_at` already orders `AddAtoms → overlay_adds → … → overlay_removes → RemoveTopology`. Correct today.
- **(i) topology among itself — already the batched form.** `RemoveTopology` is a single set op (`graph.remove(&nodes, &edges)`, one compaction, order-free, auto-takes incident bonds); adds are append-only, placed first. No intra-topology ordering problem. Correct today.
- **(ii) same overlay type — was the bug, fixed.** N single-id `Remove<X>` edits, each compacting X, so the 2nd..Nth carried stale ids. The only broken group under this analysis.

**Batching principle (the fix, implemented 2026-07-01).** No permutation rules are needed. For each compacted id space, compact it **exactly once**, by a single set op, with all inputs resolved against the pre-op state — which is what `RemoveTopology` already does for atoms + bonds + cascade. Extended to overlays: the single-id `Edit::Remove<Overlay>` variants (DAMN) were replaced by per-type batched `Edit::Remove<Overlay>s { removes: Vec<(ref, atoms, ast)> }` (stereo stays single-id until I6 adds stereo deltas). `apply_at` collects removed ids per kind and emits one batched edit each; transact resolves all ids against the pre-removal state, `OldState`-checks each, then calls `remove_<type>(&ids)` once. The `Undo::RestoreRemoved<Overlay>` variants became plural (`Vec<Removed…>` + one batch `undo_compaction` from `IdCompaction::relations`), restored via the existing plural `restore_<type>s` primitives; `IdCompaction`'s `compact_relation`/`uncompact_dense` already handle a multi-id removed set (sort+dedup) consistent with `remove_relations`. Same-type adds append and never interact with same-type removes within a reaction (a rule never adds an overlay then removes *that* overlay), so adds + removes of one type need no ordering. ("Descending id order" also fixes (ii) — the array-index-delete idiom on one Vec — but encodes correctness as an implicit ordering invariant; the set-op form is order-free and matches the crate's own topology-removal discipline.)

**Constraint-ordering finding (confirmed real, fixed 2026-07-01).** `Constraint` carries overlay ids (`Constraint::AromaticSystem(AromaticSystemId, …)`, molecule.rs:42-47). `apply_at` extended `constraint_edits` **after** `overlay_removes` with ids in the original space, but `overlay_removes` already compacted those spaces → a constraint referencing a surviving overlay whose lower-id sibling was removed got a stale id (then `RemoveTopology` compacted it a second time). Untested because reaction constraints over overlays are rare. Invariant: **all constraint edits must precede all removals (topology and overlay), so each removal's `constraints.compact` cascades into them** — or else be resolved in the post-removal space. Fix: `constraint_edits` now emitted ahead of `overlay_removes`. Constraint atom/bond refs stay correct (overlay removes don't touch those spaces; `RemoveTopology`, still last, compacts them).

**Participant-order axis (distinct from id-renumbering; fixed 2026-07-01).** Orthogonal to *which id* a delta names is *in what order* it lists an unordered relation's participants. An unordered overlay stores its participants canonically (sorted); `apply`'s `OldState` check compares the stored participant vector against the delta's order-sensitively. So a delta must carry canonical participant order, and `remap_delta` (used by compose/reverse) must re-canonicalize after relabeling. Aromatic/multicenter did (via `canonicalize_positions`, needed anyway to permute their electron ast); dative-donors and noncovalent-atoms — also `Unordered`, but with no per-participant ast — did not. Fix: `FactorOrdering::canonicalize(&mut …)` on both after remap (no-op for `Ordered`, so uniform and safe). This is a property of the *delta representation*, not the edit sequence, so it sits beside — not inside — the renumbering analysis.

**Related redundancy (noted, not folded in).** `RemoveTopology`'s cascade already removes overlays incident to a removed atom, so for those the explicit `Remove<Overlay>` is redundant with the cascade (and is where (ii) failed). Splitting removes into "atom-incident → cascade" vs "independent → explicit batch" would drop the redundancy but trade away the explicit `OldState` verification; separate axis.

**Direct transact property coverage (added 2026-07-01).** The `apply`/compose properties reach `transact` only indirectly, and the pre-existing `transaction_case_strategy` emitted single-edit cases (no overlay removes), so it would not have caught (ii). Added `overlay_transaction_strategy` (strategies.rs): a fixed 6-carbon-path base carrying two overlays of each DAMN kind, plus a valid multi-edit sequence — atom appends, an atom-charge modify, one batched `Remove<Overlay>s` per kind over a chosen subset (exercises ≥2 same-kind), and a topology removal of a chosen atom subset (cascade). Ordered adds → modify → overlay removes → topology, so every id resolves against the pre-removal base and `transact` succeeds. `transaction_edits_strategy` `prop_oneof`s it with the old single-edit cases; both feed the existing `test_molecule_builder_transact_rollback` (success round-trip via `Transaction::rollback`) and `..._unchecked` (checked vs unchecked converge). This catches the (ii) class directly and locks the edit algebra.

**Overlay-remove constraint cascade — fixed 2026-07-01.** The overlay-remove path (`remove_<overlay>s` → `constraints.compact`) *drops* a molecule-level `Constraint::<Overlay>(id, …)` whose id is removed (and remaps surviving-but-shifted refs), but the checked undo (`RestoreRemoved<Overlay>s`) captured no `CascadedConstraints` (unlike `RemoveTopology`), so `Transaction::rollback` of such a removal did not restore them. Fixed by mirroring `RemoveTopology`: the four DAMN checked `apply_edit` arms now clone the pre-removal constraints, remove, then `compact_with_update` the snapshot with the forward `IdCompaction` to capture a `CascadedConstraints`, stored on the (new) `cascade` field of each `Undo::RestoreRemoved<Overlay>s`; `apply_undo` restores the overlays then `cascade.rollback_into` re-inserts dropped / reverts remapped constraints. Verified load-bearing: `test_molecule_builder_transact_remove_aromatic_system_restores_constraint` fails when the cascade restore is removed. The overlay-transaction generator now also adds molecule constraints referencing overlays before the removals, so the round-trip properties exercise drop+restore and remap+restore.

The example cases (14 in `test_reaction_ast_compose` + `…_apply_equivalence`) are points in a large space: 4 DAMN kinds × {4 composite classes + correspondence} × {4 delta types per side} × {overlap geometry} × {scope} × error cases. A generator-based property test (proptest) over small random `(A, B, H)` is required. Below are the properties to verify and the full variant space the generator must sample.

**Properties.** Let `seq(H) = ⋃_{H' ∈ A.apply(H)} B.apply(H')` (sequential application), and products compared as canonical `MoleculeAst` (canonicalize-on-compare).

Coverage as of 2026-07-01 (tests in `umol-ast/tests/property/reaction.rs` unless noted): **P1 soundness** `compose_sound`/`_overlay` (at `host = a.lhs`); **P1 self-consistency** `compose_well_formed`/`_overlay` (`apply(C.lhs) == C.right()`); **P1 completeness deferred** (own item below). **P2** subset half `compose_rc_anchored_subset` (apply-⊆ half transitive). **P3** canonical-deltas `compose_canonical_deltas`, dangling-free `compose_dangling_free`; validity/participant-refs implicit via `apply`; dense ids structural (overlay stores are `Vec`s). **P4** `compose_determinism`. **P5** deterministic `disjoint_sum` / `disjoint_rc_anchored` cases in `test_reaction_ast_compose` (compose.rs). **P6** no-parallel-overlays `compose_distinct_overlays`. **P7** obsolete (guard removed in I3c; subsumed by P1 + `compose_dangling_free`).

- **P1 — apply-equivalence (master, `Full`).** For every host `H`: `⋃_{C ∈ compose(A,B,Full)} C.apply(H)` equals `seq(H)` **as a set** of canonical products. (Soundness: every composite product is a sequential product. Completeness: every sequential product is some composite's product.) Multiset/multiplicity equality is a stronger variant to attempt; if duplicate overlaps or symmetric automorphisms break exact multiplicity, fall back to set equality and assert multiplicity separately. **Covered:** soundness + self-consistency; **completeness deferred to doc 135.**
- **P2 — `RcAnchored` is a sound filter.** `compose(A,B,RcAnchored) ⊆ compose(A,B,Full)`, and `⋃ compose(A,B,RcAnchored).apply(H) ⊆ seq(H)`. RcAnchored is **not** complete w.r.t. `seq` (it drops composites whose overlap misses A's reaction center — the non-chaining ones); that loss is intended, so P1's completeness is asserted only for `Full`. **Covered** (subset half direct; apply-⊆ half transitive via P1 soundness).
- **P3 — well-formed composite.** Every `C ∈ compose(…)`: deltas are in canonical normal form (`C.deltas == C.deltas.canonicalize()`); `lhs_c` is a valid `MoleculeAst`; every overlay/atom/bond participant references an existing `lhs_c` (or created) id; per-kind overlay ids are dense `0..n`; no parallel overlays (§4.1); no dangling reference. **Covered:** canonical-deltas + dangling direct; validity/refs implicit via `apply`; dense ids are structural.
- **P4 — determinism.** `compose(A,B,scope)` returns an identical `Vec` on repeated calls, and is invariant under pre-canonicalizing A/B. **Covered.**
- **P5 — empty overlap = disjoint sum (`Full`).** For **disjoint** `A`/`B` (no matchable atom) the sole overlap is the empty one and `Full` yields `A ⊔ B` (atom/bond/overlay id spaces concatenated, both delta sets relabeled, no interaction). Note the empty overlap is enumerated **only when maximal** — Bron-Kerbosch returns maximal common subgraphs, so with any non-empty overlap the empty one is absent; the disjoint sum is not a member alongside non-empty overlaps. Hence tested deterministically on constructed disjoint reactants, not as a filter over arbitrary composites. **Covered** (`Full` disjoint sum + `RcAnchored` → `[]`).
- **P6 — correspondence reuses, never duplicates.** An overlap-region L_B overlay that matches an A-side overlay maps to that A composite id (no fresh id, no duplicate overlay in `lhs_c`); a created/context overlay gets a fresh id; the union of carried + created ids is exactly `0..n` per kind. **Covered:** the no-duplicate-overlay invariant is asserted per kind; dense id union is structural.
- **P7 — guard soundness (interim, until I3c).** When the narrow guard bails (a deleted atom carries an overlay), `compose` returns `[]` — a sound under-approximation; once I3c lands, P1 must hold for these inputs instead. **Obsolete** — the interim guard was removed in I3c; those inputs now fall under P1 soundness + `compose_dangling_free`.

**Variant space (generator dimensions).** The generator samples each independently; stereo is **excluded until I6**.

1. **Topology of A and of B (per atom / per bond):** unchanged · `Add` (created) · `ModifyField` · `ModifyConstraint` · `Remove`.
2. **Overlap geometry (R_A ∩ L_B):** empty · single-atom · partial · full; over **symmetric** molecules (automorphism multiplicity) and **asymmetric** ones (unique overlap).
3. **Overlay kind:** dative · aromatic · multicenter · noncovalent (each generated on A.lhs, on L_B, and as a create).
4. **Overlay composite class:** (1) lhs_A · (2) L_B context (≥1 participant off-overlap) · (3) A-created · (4) B-created · (5) L_B entirely on overlap (correspondence candidate); plus the **straddling** case (some participants on overlap, some off → class 2).
5. **Overlay delta (per side):** none (carried) · `Add` · `ModifyField` · `ModifyConstraint` (only where inhabited — dative/aromatic/multicenter; noncovalent's is uninhabited) · `Remove`.
6. **Correspondence outcome (overlap-region L_B overlay):** hit (A-side match → reuse) · miss (required-absent → skip overlap) · both-sides-edit (A modifies the overlay *and* B modifies its correspondent — the delta-interaction case).
7. **`RcAnchored` anchoring (per A overlay delta):** any overlay delta (`Add`/`Modify`/`Remove`) with participants on the overlap (anchored) · same off the overlap (not anchored, that overlap dropped).
8. **Positional realignment (aromatic/multicenter):** class-2 carry with a non-monotonic `db_atom` relabel (electron `permute` exercised) · class-1 identity carry (no permute) · `ModifyField{Electrons}` re-sort in `remap_delta`.
9. **Admissibility / error:** admissible · inadmissible (boundary bond on an A-created overlap atom) · composite-frame dangling (B deletes a shared atom whose R_A bond it cannot see) · overlay required-absent (6) · narrow-guard bail (delete an overlay-bearing atom — P7) · stereo present → bail (I6).
10. **Multiplicity / cross-kind:** multiple overlays of one kind (id density) · A creates kind K₁ while B creates K₂ (disjoint ids) · several kinds present at once.

**I3c** unblocks P7→P1 for the dangling cases; **I6** adds stereo (two more kinds, the bond-anchored site, the coset-reindex realignment) to dimensions 3–8.

**I4 — overlay span DSL** `[I1b, I2c]`
- I4a — dsl/molecule.rs: factor the six inline overlay renderers (`render_dative` …) into shared `render_<entity>_entry` (the existing TODO). `[—]` **Done 2026-07-01** — six `render_<entity>_entry(id, participants, type_edn, meta)` (dative / aromatic / multicenter / noncovalent / stereo-atom / stereo-bond) mirroring `render_bond_entry`; each `render_<entity>` is now a thin caller passing the ast's rendered `:type`. The parameter is `type_edn` (the already-rendered `:type` Edn — a `[left right]` / op-wrapped vector for spans, not an ast; `render_bond_entry`'s `value` renamed to match). Behaviour-preserving: DSL lib + molecule roundtrip tests green.
- I4b — dsl/reaction_span.rs (tree): `parse_<entity>_span_entry` over the shared entry parsers, `{:add|:modify|:remove}`-wrapped. `[I4a, I2c]`
- I4c — dsl/reaction_span.rs (format): `render_<entity>_span_entry` over I4a, span-wrapped. `[I4a, I2c]`
- I4d — dsl/reaction_span.rs (stream): streaming `read_*` for overlay span entries. `[I4b]`
- I4e — dsl/reaction_span.rs (traits): extend `ReactionSpanDsl` `FromEdn`/`ToEdn`/`FromStr`/`Display` over overlays. `[I4b, I4c, I4d]`

**I5 — structural entity refs (item 3)** `[I2c, I4]`
- I5a — dsl: the uniform `<entity>-ref` structural-map variant + per-entity resolution by constituents (item 3); six of its seven forms are overlays and ride here. `[item 3]`

**I6 — stereo deltas (the novel part)** `[I1, I2a]`
- I6a — delta.rs: `StereoAtomDelta`/`StereoBondDelta` — membership `Add`/`Remove` plus relative ops `Permute`/`Mirror`/`Swap`; lower at apply to `SetStereoAtomField::Configuration{matched, perm·matched}` via the umol-perm coset algebra. `[I1b, I2a]`
- I6b — umol-ast/ast: `TransformFrameStereoAtom`/`Bond` (#8) — coset action under a reindexing's induced ligand-frame permutation (deliberate frame op; self-inverting). `[I6a]`
- I6c — reaction.rs / reaction_span.rs: stereo apply lowering + stereo span fold (columns already in I2c). `[I6a, I2a, I2c]`
- I6d — dsl: stereo delta / relative-op surface. `[I4, I6a]`

**I7 — standalone fragment diff** `[I2f]`
- I7a — reaction_span.rs: lift `to_reaction`'s AST-diff into a standalone `diff(MoleculeAst, MoleculeAst, correspondence) -> Deltas` — the model-blind substrate for umol-graph de-aromatization (doc 131). `[I2f]`

**Out of critical path:** item 4 (constraint-store uniqueness-by-key cleanup), item 5 (FromEdn streaming audit + parse-direction naming), and umol-graph `validate` enforcement of the no-parallel-overlay invariant (spec §4.1) — all independent.
