# 134 — Reaction-application: overlays

**Metadata terminology and APIs.** The implemented persistent-metadata and
parse-time-context model, including the replacement of the concrete namespace
names used below, is specified by [doc 169](169-dsl-metadata-context-2026-07-27.md).

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
- I2f — reaction_span.rs: `to_reaction` (span→deltas) for overlays (`EntityDelta::diff`/`deltas_from_states` defaults). `[I2c, I1b]` **Done** — `to_reaction` extends the atom/bond `deltas_from_states` pattern to the four non-stereo overlays: per family, the relation set's `EntitySpan<…>` `data` column is collected into a `Vec` (no slice accessor on the relation sets) and fed to `<Delta>::deltas_from_states` with a participant closure reading the union-frame ids (`participants`/`participants_1`+`_2`) — dative `(donors, acceptor)` from `participants_2`/`participants_1[0]`, aromatic/multicenter `Vec` from `participants`, noncovalent `[a,b]`. `Added`/`Removed` → structural `Add`/`Remove`, `Modified` → the field/constraint `diff` default; all reuse the existing infrastructure (no new fold code). Stereo skipped (no `StereoAtomDelta`/`EntityFold` yet — I6a; its span columns are all `Unchanged`). Removed the stale "overlays dropped on conversion" module TODO. Tests: six overlay roundtrip cases folded into `test_reaction_span_ast_to_reaction` (`dative_add`/`aromatic_add`/`multicenter_add`/`noncovalent_add`/`_remove`/`_modify`) covering each family's participant closure + Add/Remove/Modify states; umol-ast 3949 green. **Stereo `to_reaction` closed 2026-07-01** (I6a-c landed the stereo deltas): the two stereo families are recovered without `EntityFold` (they have none — the relative ops have no `EntityOp` image) via a direct column walk — `Added`/`Removed` carry site + ligand frame, `Modified` uses `EntityPatch::diff` (absolute `ModifyField`); mirrors the render logic. Discovered via the `span_roundtrip` property (stereo spans silently lost their deltas before). Tests: `stereo_atom_add`/`_remove`/`_modify` + `stereo_bond_add` `to_reaction` roundtrip cases.

**I3 — composition for overlays** `[I1, I2, I0e]`
- I3a — compose.rs: unify `remap_delta`/frame algebra onto the general `Remapping`; extend the four-class composite frame to overlay relation ids. `[I1d, I0c, I0d]` **Scope settled 2026-06-30 (two forks).** (1) *Vehicle* — keep the `IdRemapping` HashMap bundle as `remap_delta`'s input; the "unify onto general `Remapping`" clause is superseded by I1d (graph-core's dense `Remapping` can't carry overlay relation ids, and its `apply_remapping` over relation sets is only needed for the lhs-overlay carry, which is out of scope here). No `remap_delta` signature change. (2) *Boundary* — **created overlays only**: I3a fills composite-frame classes (3) (A-created) + (4) (B-created) overlay ids — disjoint by construction — and guards by bailing the whole composition when **either reactant lhs carries any overlay** (`a.lhs.has_overlays() || b.lhs.has_overlays()`). The guard is on the lhs, not the deltas: `lhs_c` is built from `raw_graph()` only, so *any* lhs overlay — modified, removed, **or untouched** — would be silently dropped; the lhs check subsumes a delta-scan (you can't modify/remove an overlay absent from the lhs) and created overlays live in the deltas so they still compose. Classes (1) (lhs_A) + (2) (L_B context), `lhs_C` overlay materialization, and the overlap-region overlay correspondence move to I3b. (A first cut guarded only overlay `ModifyField`/`ModifyConstraint`/`Remove` deltas — unsound: an untouched lhs overlay passed and was dropped. Corrected to the lhs guard 2026-06-30, with a `untouched_overlay` test case.) **Done 2026-06-30** — created-overlay boundary (classes 3/4) + the lhs guard; both then subsumed by I3b's `place_overlays` (which materializes classes 1/2 and drops the broad guard).
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

**I4 — overlay span DSL** `[I1b, I2c]` **Done**
- I4a — dsl/molecule.rs: factor the six inline overlay renderers (`render_dative` …) into shared `render_<entity>_entry` (the existing TODO). `[—]` **Done 2026-07-01** — six `render_<entity>_entry(id, participants, type_edn, meta)` (dative / aromatic / multicenter / noncovalent / stereo-atom / stereo-bond) mirroring `render_bond_entry`; each `render_<entity>` is now a thin caller passing the ast's rendered `:type`. The parameter is `type_edn` (the already-rendered `:type` Edn — a `[left right]` / op-wrapped vector for spans, not an ast; `render_bond_entry`'s `value` renamed to match). Behaviour-preserving: DSL lib + molecule roundtrip tests green.
- I4b — dsl/reaction_span.rs (tree): `parse_<entity>_span_entry` over the shared entry parsers, `{:add|:modify|:remove}`-wrapped. `[I4a, I2c]` **Done 2026-07-01** — six `parse_<entity>_span_entry` mirroring `parse_bond_span_entry`: a `full` closure reuses the shared `parse_<entity>_entry` for `None`/`:add`/`:remove`, and `:modify` re-parses participants + a `[left right]` `:type` pair. `SpanInput` gained the six overlay `Vec` fields wired into `parse_span_input`.
- I4c — dsl/reaction_span.rs (format): `render_<entity>_span_entry` over I4a, span-wrapped. `[I4a, I2c]` **Done 2026-07-01** — six `render_<entity>_span_entry` reuse `render_<entity>_entry` + `single_key_map` (Unchanged bare, add/remove/modify wrapped; modify's `:type` a `[left right]` vector), walked in `render_span_edn` via each overlay set's `relation_ids()`/participants/`data`. `render_stereo_ligand` generalized `StereoLigandView`→`StereoLigand` (view carried no data the renderer used) so molecule + span share it; molecule callers switched to `ligand_frame()`. Added `pub(crate)` `stereo_atoms()`/`stereo_bonds()` accessors on `ReactionSpanAst`.
- I4d — dsl/reaction_span.rs (stream): streaming `read_*` for overlay span entries. `[I4b]` **Done 2026-07-01** — the span entry grammar is tree-only (reuses the molecule entry parsers), so `read_span_input` buffers each overlay section element to an `Edn` and dispatches to the tree `parse_<entity>_span_entry`, same as atoms/bonds.
- I4e — dsl/reaction_span.rs (traits): extend `ReactionSpanDsl` `FromEdn`/`ToEdn`/`FromStr`/`Display` over overlays. `[I4b, I4c, I4d]` **Done 2026-07-01** — `SpanInput::into_ast` registers overlay `:id`s, resolves refs (`.into_ast`), builds the six relation-set span columns, enforces per-side participant consistency (the bond check generalized to overlay participants; stereo-bond also checks its site bond), and fills `EntityCounts`. The `FromAst`/`IntoAst` config bridge maps each overlay set's `EntitySpan` payload through the overlay `*Dsl` (three generic `map_*_span` rebuild helpers, one per relation-set shape; participants carry through). Tests: six `test_parse_<entity>_span_entry`, six `test_render_<entity>_span_entry`, overlay cases across both roundtrips, a dative `into_ast` build case, and two per-side-consistency error cases. umol-ast 4048 lib green; workspace builds clean.

**I5 — structural entity refs (item 3)** `[I2c, I4]` — same kernel as I7; design in §7. **Not built; reopened and tracked in doc 135.**
- I5a — dsl: the uniform `<entity>-ref` structural-map variant + per-entity resolution by constituents (item 3); six of its seven forms are overlays and ride here. Resolution = the §7.3 `find_by_participants` kernel. `[item 3]`

**I6 — stereo deltas (the novel part)** `[I1, I2a]` **Done**
- I6a — delta.rs: `StereoAtomDelta`/`StereoBondDelta` — membership `Add`/`Remove` plus relative ops `Permute`/`Mirror`/`Swap`; lower at apply to `SetStereoAtomField::Configuration{matched, perm·matched}` via the umol-perm coset algebra. `[I1b, I2a]` **Core done 2026-07-01; Delta-sum wiring → I6c.** Settled shape (sign-off 2026-07-01): the relative ops are **eager coset-algebra methods on the AST**, not lazy operator chains and not via `canonicalize` — `StereoCosetAst::{apply(kind, σ), swap(kind), mirror(kind)}` map through the variants (`Lit`/`LitSet` move eagerly via `StereoKind::act` = `space.reindex`; `Undetermined` fixed; an open `Term(Var)` keeps one operator layer), lifted through `StereoConfigurationAst` and `StereoAtomAst`/`StereoBondAst` (constraints untouched — the aromatic/multicenter `permute` parallel). Named `apply` (the `^` op), **not** `permute` (that's the electron-reorder word). `Swap` = the kind involution (`StereoKind::involution`, π-independent, no arg — not a `{a,b}` transposition); `Mirror` = enantiomer (`mirror_permutation`, π-independent). `StereoKind::act` made `pub` (`#[allow(unused)]` dropped). Delta enums = `Add`/`Remove`/`ModifyField`/`ModifyConstraint` (the DAMN set-ops) + the relative ops `Apply{kind, permutation}`/`Swap{kind}`/`Mirror{kind}` (no pre-state) with `inverse` (`Apply(σ)→Apply(σ⁻¹)`, `Swap`/`Mirror` self-inverse). The relative ops **carry `kind: StereoKind`** (sign-off 2026-07-01: the coset algebra is parametrized by kind, so a relative op is uninterpretable without it, and the absolute arms carry the kind inside their config/ast) — this makes every stereo delta kind-determinate and lets `canonicalize` fully compose them. Field named `permutation`, not `perm`/`relabeling` (no abbreviations; "relabel" is the id-remap word — sign-off 2026-07-01; `Coset::reindex`'s parameter renamed to match). `StereoAtomFieldChange`/`StereoBondFieldChange` gained `PartialOrd, Ord` (needed for the delta's `Ord`; the DAMN field-changes got it in I1b). 30 unit tests, umol-ast 4091 green. **Not** wired into the `Delta` sum: the DAMN deltas integrate via the `EntityPatch`/`EntityFold` traits (canonicalize-fold / `remap_delta` / diff / apply all dispatch through them), built around exactly the four DAMN arms; the relative ops don't fit, so the sum extension + `From` + the relative-op flow through canonicalize/remap/apply is I6c's substance and lands there.
- I6b — umol-ast/ast: `TransformFrameStereoAtom`/`Bond` (#8) — coset action under a reindexing's induced ligand-frame permutation (deliberate frame op; self-inverting). `[I6a]` **Done 2026-07-01 (AST-level only, per sign-off).** `StereoAtomAst`/`StereoBondAst::transform_frame(before: &[StereoLigand], after: &[StereoLigand]) -> Self` (via the `stereo_element!` macro) = `self.apply(Permutation::between(before, after))` — the induced frame permutation (`between` asserts the two frames are the same ligand multiset reordered; genuine ligand-set changes are membership, not a frame permutation) fed to I6a's `apply`. Self-inverting falls out (`between(after, before) = between(before, after)⁻¹`). Carries frames (not just atom indices), so it covers virtual↔explicit position swaps. Tests: identity / transposition (Lit 0→1) / even-cycle (0→0) / virtual↔explicit swap / self-inverse round-trip. The `Edit::TransformFrameStereoAtom`/`Bond` variant that carries the frames + its `transact`/`Undo` (doc 131's "Edit variant pair") rides with I6c, not here.
- I6c — reaction.rs / reaction_span.rs: stereo apply lowering + stereo span fold (columns already in I2c). `[I6a, I2a, I2c]` **Done 2026-07-01 (all sites green; umol-ast 4101).** Wired the I6a delta enums into `Delta` + the flow through every dispatch site (design below realized in full):
  - **sum** — add `StereoAtom`/`StereoBond` to `Delta`, `From` impls, `Delta::inverse` arms (mechanical).
  - **four DAMN arms** — `EntityPatch`/`EntityFold` via the existing `diff_field_ops!`/`fold_field_ops!` macros with `{ Configuration => configuration }`; so `to_reaction` **diff produces only absolute `ModifyField`** (relatives come from authoring/compose, never from diffing concrete configs).
  - **relative ops in `canonicalize` — option (b), enabled by the delta's `kind`**: normalize `Swap{k}`→`Apply(k, k.involution())`, `Mirror{k}`→`Apply(k, k.mirror_permutation())`; compose consecutive `Apply(k,σ)∘Apply(k,τ)`→`Apply(k, σ.compose(τ))` (drop if identity — a right action, `act(act(c,σ),τ)=act(c,σ∘τ)`); fold a relative after an absolute `ModifyField`/`Add` by evaluating it onto that config (kind from the config). A **bespoke stereo fold path** (the four arms still ride `fold_group`; the relative composition is stereo-specific, since the coset algebra doesn't fit the entity-agnostic `EntityOp`).
  - **`remap_delta` — id relabel only, no coset transform.** The ligand factor is `Ordered` (frame not re-sorted on remap), so a coset stays valid under an atom-id relabel; `Apply`'s permutation is position-space (untouched). Map `map_stereo_atom`/`map_stereo_bond`, site, ligand atom ids.
  - **apply lowering (`apply_at`) — same-frame, no reconciliation** (mirrors aromatic/multicenter `ModifyField`, which pass through with just the host id; "same incidence frame", doc 104). Four arms → the existing `Add`/`Remove`/`ModifyStereoAtomField` edits; relatives → resolve against the matched host config via I6a `apply`/`swap`/`mirror` → `ModifyStereoAtomField::Configuration{matched, op(matched)}`.
  - **span fold (`to_reaction_span`)** — `apply_stereo_atom_change`/`apply_stereo_bond_change` (the I2d analog): fold each stereo delta onto the lhs span column (relatives via I6a `apply` on the lhs config).
  - **DPO stereo-dangling** — extend the deleted-host-atom incidence check (reaction.rs:192 TODO) to stereo site + ligand incidence.
  - **compose / dpo / dsl arms** — compose already bails on `lhs` stereo; extend the bail to cover any stereo *delta* (an Add-only rule has no lhs stereo) so compose's stereo arms are honestly unreachable; the dpo `validate_reaction`/`_span` and dsl/reaction Delta arms get real (small) handling. Compose *aliasing* stereo (the only `TransformFrame` consumer — apply is same-frame, remap is `Ordered`-preserving) stays out; the `Edit::TransformFrameStereoAtom`/`Bond` variant + transact lands with compose-stereo, not I6c.

  **Done 2026-07-01 (all sites green; umol-ast 4101 tests, clippy clean).** Landed: `Delta` sum + `From`(direct) + `inverse` arms; `EntityPatch` for the two stereo deltas (four arms, `{ Configuration => configuration }`; `apply_constraint` via `remove_by_key`); `remap_delta` stereo arms (id-relabel only); **two new `Edit::ModifyStereo{Atom,Bond}Constraint`** variants (sign-off) with unchecked/checked `transact` + `Undo::ApplyEdit` inverse + builder `apply_modify_stereo_{atom,bond}_constraint`. Unit tests: stereo `diff`/`apply_field`/`apply_field_error`.

  **Discovered mid-build — stereo removal must batch** (sign-off): singular `Edit::RemoveStereoAtom`/`RemoveStereoBond` hit the same dense-compaction stale-id bug the DAMN batched-remove fix solved (multiple same-kind removals). Refactor to `RemoveStereoAtoms`/`RemoveStereoBonds { removes: Vec<…> }` + `Undo::RestoreRemovedStereo{Atom,Bond}s { removed: Vec<…>, undo_compaction, cascade }`, mirroring `RemoveDativeBonds` (transact checked ~812: resolve+validate each, `IdCompaction::relations` with stereo ids in the 5th/6th slot, clone `pre_constraints`, `remove_stereo_{atoms,bonds}(&ids)`, `pre_constraints.compact_with_update(&forward)` → `cascade`). **Done 2026-07-01:** `edit.rs` `Edit`+`Undo` variant shapes; transact `apply_edit_unchecked` + checked-with-cascade + `apply_undo` (`cascade.rollback_into`); builder plural `restore_stereo_{atoms,bonds}` made `pub(super)`, singular wrappers dropped; 3 tests migrated to the batched shape. Batching compiles; only the 4 Delta dispatch sites remain red.

  **`apply_at` — done 2026-07-01.** First pass: four set-ops → the `ModifyStereo*Field`/`ModifyStereo*Constraint` edits; relative ops read the matched host coset (`host.stereo_{atom,bond}(host_id).coset()` + the delta's `kind` → I6a `apply`/`swap`/`mirror` on `StereoConfigurationAst`) → absolute `Configuration` edit; `Remove` tracks `removed_host_stereo_{atom,bond}`. Second pass: `Add`→`Edit::AddStereo{Atom,Bond}` (new `bond_ref` helper for the bond-sited stereo bond), `Remove`→batched `RemoveStereo{Atoms,Bonds}`. **DPO stereo-dangling done** — a deleted atom's `host.stereo_{atoms,bonds}().incident_ids(atom)` (site or ligand) must be co-deleted. (Corrected a wrong assumption: stereo incidence *is* exposed — on the stereo views `StereoAtom/BondViews`, not `AtomView`.)

  **`to_reaction_span` span-fold — done.** Stereo accumulation maps + fold section (aromatic-style Removed/Modified/Unchanged + appended Added) via `apply_stereo_{atom,bond}_change` (= `apply_aromatic_change` + the three relative-op arms `*ast = ast.apply/swap/mirror`); added-entity site/ligand ids mapped to the union frame (new `bond_index` for the bond-sited stereo bond).

  **`(b)` canonicalize — done.** Bespoke `fold_stereo_{atom,bond}_group` (stereo isn't `EntityFold` — the relative ops have no `EntityOp` image): *created* path seeds from `Add` and applies each op via `apply_stereo_*_change` (`Add`+`Remove` cancel); *preserved* path splits config ops (folded by `fold_stereo_config`) / constraints (by key) / `Remove`. `fold_stereo_config` runs rules ii–vi: relatives compose by permutation (`σ∘τ`), a set absorbs a leading relative by `apply(old, σ⁻¹)` (vi) and transforms `new` through a trailing relative (v), two sets fuse (`new₁==old₂`); the net is classified by the factored `StereoKind::canonicalize_permutation(g) -> Option<CosetOp>` (priority Mirror > Swap > Apply — factored out of `canon_coset` per rule iv, reused by both; `CosetOp` = `Swap`/`Mirror`/`Apply(Permutation)`, the priority chosen by the fn not intrinsic to the permutation). A `Remove` reverts the net action + constraints onto the removed (original) ast. Tests: `swap∘swap`/`mirror∘mirror`/`apply²`/`swap∘mirror` cancel, `swap`→`mirror` (chiral priority), `Add`+swap → transformed `Add`, swap+`Remove` → reverted `Remove`.

  **dsl/reaction.rs render** — stereo added to the existing overlay-delta `todo!` group (the `render_deltas` `todo!` originally read "I4: overlay delta DSL rendering" — a **mislabel**: I4 was the overlay *span* DSL (`reaction_span.rs`), not the operational reaction DSL's overlay *deltas* (`reaction.rs`), which neither 133 nor I4 implemented for any overlay). Folded into I6d's expanded scope below.

  **Batched-remove refactor (2026-07-01):** singular `RemoveStereo{Atom,Bond}` → batched `RemoveStereo{Atoms,Bonds} { removes: Vec<…> }` (type aliases `Stereo{Atom,Bond}Removal`) + `Undo::RestoreRemovedStereo{Atom,Bond}s { …, cascade }`, mirroring `RemoveDativeBonds` (fixes the same multi-remove stale-id hazard). **DPO stereo-dangling done** — via `host.stereo_{atoms,bonds}().incident_ids(atom)` (stereo incidence is on the stereo views, not `AtomView`).
- I6d — dsl: **operational reaction-DSL overlay-delta surface, all six overlays (DAMNSS)** + the stereo relative-op syntax. `[I4, I6a]` **Scope expanded 2026-07-01** (was "stereo delta / relative-op surface"): `render_deltas` + `parse_delta_input` in `dsl/reaction.rs` handle only `:atom`/`:bond`/`:constraint`; every overlay delta is `todo!`/"unknown reaction delta". Since the render + parse dispatch is one shared surface, cutting it to stereo-only would leave DAMN behind the same `todo!` — the exact deferral-then-bite pattern. So I6d covers all six: `{:<overlay> {:add <entry> | :remove <ref> | :modify [<ref> <partial>]}}` (reusing the I4a `render_<entity>_entry` + the entry parsers), plus the stereo relative ops as verbs `:swap`/`:mirror`/`:apply` (permutation in disjoint-cycle notation, reusing `perm_cycles`/`Permutation` `Display`).

  **Done 2026-07-01 (umol-ast 4146, workspace clippy clean).** Singular delta keys `:dative-bond`/`:aromatic-system`/`:multicenter-bond`/`:noncovalent-bond`/`:stereo-atom`/`:stereo-bond` (matching `:atom`/`:bond`, not the plural collection keys). Landed:
  - `update()` on all six overlay ASTs (partial-merge for `:modify` diff, mirroring `AtomAst::update`; aromatic/multicenter merge spin field-wise); `EntityCounts` overlay allocators; `ReactionMetadata` six overlay id-maps (created-id roundtrip).
  - 24 `DeltaInput` variants + resolution arms (reuse the molecule entry inputs/parsers + the `define_ref!` overlay refs); `Partial<Overlay>Dsl` per overlay for `:modify`.
  - Tree `parse_delta_<overlay>_input` + streaming `read_delta_<overlay>_input` + dispatch; render arms; `render_deltas` extended.
  - **stereo relative ops carry an explicit `kind`** (sign-off 2026-07-01): `:swap [ref kind]` / `:mirror [ref kind]` / `:apply [ref kind "(0,1)"]` — positional vector, kind keyword. Unlike deriving it from the lhs ref, the explicit kind decouples the delta from the lhs (works when the lhs config is `Undetermined`), is symmetric with the resolved `StereoXDelta` which already carries `kind`, and lets `:apply` carry a real `Permutation` (`parse_permutation(input, kind.degree())`, a new cycle-string wrapper mirroring `parse_stereo_coset`) instead of degree-agnostic cycles. No lhs kind lookup at resolution (tier-1 only; no kind-vs-lhs consistency check).
  - Tests: `test_reaction_dsl_from_edn_to_edn_roundtrip` gains a case per overlay × add/remove/modify + stereo swap/mirror/apply.
- I6e — constraint removal via `:modify` (all six overlays). `[I6d]` **Done 2026-07-01 (umol-ast 4146).** A `ModifyConstraint{new:None}` renders by putting an *undetermined* constraint in the partial (`old.as_undetermined()`), which the partial DSL must render and the parser read back.
  - **DAMN** already round-tripped: dative `#a*`/`#R*`, aromatic/multicenter `#e*` (added `as_undetermined()` to `AromaticSystem`/`MulticenterBond` constraints; dative/atom/bond already had it).
  - **Stereo constraint DSL redesign** (prerequisite, sign-off 2026-07-01): `LigandSymmetryAst.member: MemOp` → `present: BooleanAst` and `FluxionalityAst` gained `present: BooleanAst` (giving both an undetermined form; `MemOp::In`→`Lit(true)`, `NotIn`→`Lit(false)`). Lattice `meet`/`merge_same_key`/`is_ground` lifted to three-valued via `BooleanAst`; `as_undetermined()` added to the stereo constraint enums. Parametrized notation (`#tag <parameter> <trailing>`): `!` moved from a leading member-marker to a trailing `<bool>` (`""`/`+`→true, `!`→false, `*`→undetermined) on `#p`/`#f`; `#o` renders the ligand pair before the (unchanged) topicity relation; `#g` unchanged. Topicity's 3-literal lattice confirmed to need only single-literal / not-literal (`!x`) / `*`. `:present` EDN key mirrors the `#a` `BooleanAst` form; fluxionality EDN became a map to carry it. Downstream `umol-graph` validator updated (`present` → `Undetermined` is vacuous).
  - **Stereo `:modify` kind-carry** (the hard case: stereo constraints render/parse against the config's kind, which a constraint-only modify leaves unchanged): a *separate* partial parser (`partial_stereo_{atom,bond}` via `partial_stereo_configuration`) makes the coset **optional** (omitted = unchanged) and the kind mandatory once anything else appears (`"*"` alone, or `"Th#o(0,1)="`; `*#o…` rejected). `StereoConfigurationAst::update` merges config field-wise (an undetermined coset in the partial keeps the base coset), so restating the coset isn't required. The partial renderer omits the unchanged coset and renders *all* constraints (incl. undetermined). Tests: stereo topicity removal/change + ligand-symmetry removal round-trips, `StereoConfigurationAst::update` (5 cases), partial parser/render (14 cases).
  - **`kind` on `ModifyConstraint`** (sign-off 2026-07-01, superseding the render-time lhs borrow): stereo constraints render/parse against the config's kind, and a constraint-only modify's `ModifyConstraint` carried neither the kind nor a config to derive it from — so the first cut borrowed it from the lhs at render (`render_deltas(lhs, …)` with `lhs.stereo_atom(id).kind()`, which `.expect`s a concrete kind → **panics** on an `Undetermined`-geometry center with a kind-free `Topicity`/`Stereogenicity` constraint). Instead `Stereo{Atom,Bond}Delta::ModifyConstraint` now carries `kind: Option<StereoKind>` (like the relative ops carry `kind`, but `Option` — a kind-free constraint on an open geometry genuinely has none). `EntityPatch::modify_constraint` can't supply it (no config in the signature), so the two stereo deltas **override `EntityPatch::diff`** to stamp `kind` from the entity's config; the canonicalize fold threads it through the constraint group, `inverse`/`remap` carry it. Render reads it off the delta — no lhs param, no panic. `apply`/`canonicalize`/`diff` ignore it (serialization-only). The `Undetermined`+kind-free-constraint case degrades (constraint elided, no round-trip) instead of panicking — full support would ungate the kind-free constraint renderer.

**I7 — standalone fragment diff** `[I2f]` — unified with I5; design in §7.
- I7a — reaction_span.rs: `MoleculeAst::diff(&other, &MoleculeCorrespondence) -> Deltas` (= `align` → `to_reaction`) — the model-blind substrate for umol-graph de-aromatization (doc 131). `[I2f]`
- I7b — graph-core: the correspondence carrier + `find_by_participants` kernel (§7.3). `[I2f]`
- I7c — umol-ast: `MoleculeCorrespondence` (base + `Total` wrapper + embedding view), `induce`/`from_maps`, surface syntax (§7.4). `[I7b]`
- I7d — the `align` span constructor + `Sp::correspondence()` + `ReactionAst` import constructor (§7.5). `[I7c]`

## 7 — Correspondence & diff: I5 + I7 unified (design)

I5 (structural entity refs) and I7 (standalone diff) are one kernel seen twice: **resolve an entity by its constituent atoms**. I5 does it in one molecule (a `<entity>-ref` names a bond by its endpoints); I7 does it across two molecules under an atom mapping (a left bond corresponds to the right bond over the mapped endpoints). `substructure::verify_overlays` already runs this kernel for the embedding case: given a subiso `atom_map`, it matches each pattern overlay to a host overlay via `<collection>.connecting(mapped_atoms)`. So I7 generalizes existing code and I5 re-exposes it as a DSL surface.

### 7.1 Shared model — the correspondence, its producers, and the span

graph-core already names the object: `CommonSubgraph` is *"a vertex correspondence between two graphs"* — `mapping: Vec<(NodeId, NodeId)>` + `edge_count`. Every common-subgraph task produces one:

| producer | shape | correspondence character |
|---|---|---|
| `subgraph_isomorphisms` (VF2/Ullmann/RI/…) | `Vec<Vec<usize>>` | total on the query, injective (an **embedding**) |
| `maximum_common_induced/edge_subgraph[s]` (McGregor) | `CommonSubgraph` | maximal-size partial bijection |
| `maximal_common_subgraphs` (Bron–Kerbosch) | `Vec<CommonSubgraph>` | every non-extendable partial bijection |
| `enumerate_common_subgraphs` (backtracking) | `Vec<CommonSubgraph>` | every partial bijection |
| (future) MCS-derived atom map | `CommonSubgraph` | the reaction atom map when none is supplied |

All five yield **a partial bijection between two node id spaces**, differing only in totality / injectivity / maximality, not in kind. They return different shapes today (`CommonSubgraph` from MCS, raw `Vec<usize>` from subiso); the unification is one carrier — `Correspondence` (§7.3) — that all producers emit (MCES wrapping it in `Mcs` to cache its edge-count objective), plus the §7.2 reads.

The umol-ast **span** (`EntitySpan`, `ReactionSpanAst`) is *not* a fourth producer — it is the correspondence **lifted with values and a direction**:

| correspondence (graph-core: valueless, undirected) | span (umol-ast: AST-valued, L→R) |
|---|---|
| mated pair `(l, r)` | `Unchanged` (values equal) / `Modified{left,right}` (differ) |
| left-exposed `l` | `Removed` |
| right-exposed `r` | `Added` |

`align(L, R, corr)` = **lift** (refine each mated pair by comparing AST values, orient each exposed node); `Sp::correspondence()` = **forget** (drop values, merge Unchanged/Modified back to mated). The correspondence is the transformation-free substrate; the span is its reaction-specific enrichment. This is why the span carries `Added`/`Removed` (a direction) and the correspondence must not.

### 7.2 Vocabulary for the two-sided membership — matching-theoretic (settled)

`Presence` / `classify` are too generic and lack two-sidedness. **Frame the correspondence as a matching in the bipartite graph `L ⊔ R`** (left nodes on one side, right nodes on the other, an edge per corresponding pair). Then the standard matching vocabulary names all three classes of `L ⊔ R`:

- **mated** — paired with a partner on the other side (the DPO interface K);
- **left-exposed** / **right-exposed** — unpaired, present on only one side. *Exposed* is the matching-theory term for an unmatched vertex; the bipartite framing makes it read naturally.

**Matching subsumes the partial-function view rather than competing with it.** A correspondence *is* a partial bijection `L ⇀ R`, and `mate(x) -> Option<NodeId>` *is* that partial function — `Some` = in the domain, `None` = exposed. So the natural "apply a partial map" intuition survives (via `mate`), and matching *adds* the two nouns partial-fn lacks (`exposed`, split by side); partial-fn would name only `domain`/`image` — two projections of the mated class — and force `L∖domain` / `R∖image` set-diffs for the two classes a diff most cares about. No tradeoff; matching dominates. Consistent with graph-core's existing `Matching` (`mate`/`is_matched`).

**Excluded: `deleted` / `created`** (DPO). They name the same three classes but presuppose *direction* — they are the span's `Removed` / `Added`. The correspondence is transformation-free (§7.1), so its class names must be direction-free; `exposed` says "unpaired," not "deleted."

**Surface — three accessors, mapped 1:1 to the span build** (the primary form; the per-element enum `{ Mated { left, right }, LeftExposed(NodeId), RightExposed(NodeId) }` / `mate(x)` is the secondary "what is this id" query):

| accessor | yields | span lift (`align`) |
|---|---|---|
| `mates()` | `&[(left, right)]` | compare AST values → `Unchanged` / `Modified` |
| `left_exposed()` | `&[NodeId]` | `Removed` |
| `right_exposed()` | `&[NodeId]` | `Added` |

`classify` (the membership walk) is dropped in favor of these three accessors.

### 7.3 graph-core primitive (point 4)

Two reusable additions, consumed by I5, I7 `induce`, and `verify_overlays`:

1. **structural-match kernel — `find_by_participants`, per relation-set type and per factor.** Relations (one factor, fixed or variable arity): `find_by_participants(&self, &[P]) -> Option<RelationId>`. Birelations (two factors): `find_by_participants(&self, &[P1], &[P2]) -> Option<RelationId>` — one slice per factor, in factor order (dative `(acceptor, donors)`, stereo `(site, ligands)`). Participant type is generic per factor (`P` = `NodeId` / `EdgeId` / `StereoLigand`; stereo-atom is `(&[NodeId], &[StereoLigand])`). Uniqueness (§4.1) ⇒ ≤1 hit. Inherent methods per relation-set type, called from the caller's per-family arms (`verify_overlays` / `induce`) — no unifying trait. Ad hoc today as `<collection>.connecting(atoms)` in umol-ast; lift to graph-core.

   - **Per-factor, never the flat union.** Matching the union of all participants across factors is unsound — it can't tell donors from acceptor, `A→B` from `B→A`. The two-factor structure must be respected.
   - **Each factor matched as an unordered *multiset* — key on identity, not the storage marker.** Every current overlay identifies by its participant multiset: the five `Unordered` factors trivially, and stereo (the one `Ordered`-with-multiplicity factor) too — its ligand order is the coset **frame** (a *value*), the stereo atom's identity being the ligand **bag** (§4.1). So the kernel keys on the participant *multiset*, **deliberately ignoring** the factor's `Ordered`/`Unordered` marker. **Multiset, not set:** virtual ligands (implicit H, lone pairs) legitimately repeat — duplicate virtuals give a valid coset (vacuous in Th, not Ct); only duplicate *explicit* ligands are an error (caught at stereo construction). Multiplicity is real data; a set would collapse it. **Comparison = sort-both-and-compare** (= multiset equality): canonicalize a *copy* of the query and of the stored participants by the factor's total order (`StereoLigand: Ord`) and compare element-wise — exactly `Unordered::canonicalize` (`sort_unstable`), applied to a copy so the stored `Ordered` frame stays intact. Node/edge factors never have duplicates (a repeated node is an error, like a repeated explicit ligand), so there multiset = set — one uniform rule, no dedup. **No `as_ordered` argument:** it is redundant with the marker, *and* following the marker would wrongly order-match stereo (a query in a different frame order would miss the center). Ordered *match* (order-is-identity) is the faithful lookup only for a relation whose order genuinely is its identity (a directed / sequenced hyperedge) — umol has none, so it is not provided.
   - **Distinct from the reverse read** `incident_ids(node) → relations` (unordered union — removal cascades). Three directions total: `participants_1/2(rel) → nodes`, `incident_ids(node) → relations`, `find_by_participants(factors) → Option<RelationId>`.

   *Note (resolved):* `Ordered::canonicalize` is a no-op ("position is the datum", relation.rs), so the stereo relation set keys by exact ligand order rather than structurally enforcing §4.1 set-uniqueness — but that is not a problem: the multiset match is frame-order-independent (finds the center regardless of stored frame order), duplicate *explicit* ligands are a construction error, and duplicate *virtual* ligands (implicit H / lone pairs) are valid (vacuous in Th, not Ct).
2. **`Correspondence` carrier + reads** — its own module (not limited to common subgraphs): the vertex-level partial bijection over two graphs, the §7.2 mated/exposed membership, the induced edge correspondence (an edge is mated iff its endpoints are mated to an edge on the other side — `find_by_participants` under the node map), and derived reads `node_count() = mates().len()`, `shared_edge_count(&a, &b)` (edge count is a projection onto the graphs — line 767's loop — not stored). Producers:
   - **subiso** → `Vec<Correspondence>` — the total-on-left, injective case (`left_exposed()` empty, `right_exposed()` = host outside the match);
   - **MCIS**, **maximal / complete enumeration** → `Correspondence` / `Vec<Correspondence>` — objective is vertex count (`= mates().len()`, already on the correspondence) or none;
   - **MCES** → `Mcs { correspondence, edge_count }` — the *only* producer whose objective isn't free on the correspondence; branch-and-bound already holds `edge_count`, so `Mcs` caches it (delegates `mates()`/`node_count()`, adds `edge_count()`);
   - an MCS-based atom-map derivation → `Correspondence`.

   `CommonSubgraph` **retires** — split into `Correspondence` (the general common subgraph) and `Mcs` (the MCES objective-result). `Matching` is *not* a producer — a single-graph matching (`edges()`, `is_perfect`), vocabulary donor only. Keep `Correspondence` lightweight (pairing + derived-on-demand reads) so subiso emits one per match on the hot enumeration path.

### 7.4 `MoleculeCorrespondence` (umol-ast)

Name: `MoleculeCorrespondence`.

**Shape.** A per-entity partial bijection between two `MoleculeAst` id spaces — atoms + bonds + six overlays, exactly `MoleculeEmbedding`'s eight `host_*` fields but *partial on both sides* (not total-on-sub). The **atom** part is a node-level vertex correspondence — a set of `(AtomId, AtomId)` pairs, the same shape as `CommonSubgraph::mapping()`; the seven entity parts are **induced** from it by structural match. `MoleculeCorrespondence` layers the atom-map derivation and the entity induction on top of the graph-core carrier.

**One type, three cases** (point iv — combine, don't fork; point v — a marker wrapper, not a type per case, per the MoleculeAst-wrapper experience):

- `MoleculeCorrespondence` — general, partial on both sides.
- **Total is a predicate, not a wrapper** (settled 2026-07-02). `Correspondence::is_total()` — `mate_count == left_count == right_count`, no exposed on either side, O(1) — marks the balanced case: same atoms, so a diff needs **no union frame** (R's atoms map 1:1 into L's; bonds/overlays still move within the shared frame — de-aromatization, doc 131). The general path subsumes it (the exposed branches are O(1) dead checks); the union-frame skip is a fast-path branch keyed on `atoms().is_total()`. A `Total(_)` wrapper would add only a *compile-time* no-Add/Remove guarantee for a caller that needs it — none does yet, so a predicate suffices; add the wrapper later if one appears.
- the **embedding** (injective, total-on-one-side) case is what `MoleculeEmbedding` already is; it folds in as a view/wrapper, retiring the standalone struct (subiso + `verify_overlays` build a `MoleculeCorrespondence` directly).

**Constructors.**
- `induce(L, R, atom_map) -> MoleculeCorrespondence` — derive the seven entity maps from the atom pairs by structural match (the I5 / `verify_overlays` kernel). SMIRKS / GML / map-number **import entry point**.
- `from_maps(…)` — the fully-explicit per-entity form `induce` materializes.

**Surface syntax (settled).** EDN-shaped (top-level key), following the span precedent but collapsing to **atom mated-pairs only**:

```clojure
{:correspondence {:atoms [[0 3] [1 1] [:c4 :c9]]}}
```

- A vec of `[left-ref right-ref]` **mated** pairs; each ref is the standard `<entity>-ref` (`int | keyword`, mixable). The `:atoms` wrapper is kept (explicit, room to grow) even though the surface is atoms-only.
- **No exposed entries, no per-family keys, no `nil`.** The whole surface is the atom map because everything else is derived at resolution: exposed atoms are inferred (every non-mated atom *is* exposed), and bond/overlay maps are `induce`d (determined by the atom map — a bond mates iff its mapped endpoints mate).
- **Two-phase, dissolving the parse-time L/R dependence** (the `*Input` pattern, cf. `SpanInput::into_ast`): parse → `MoleculeCorrespondenceInput { atoms: Vec<(AtomRef, AtomRef)> }` (no L/R needed) → `resolve(L, R)` does keyword→id, infers exposed, and induces every entity map. Resolution is where L/R legitimately enter.
- **Structural-ref pairs (once I5) are allowed but tautological.** Atoms are the base (no structural form), so a structural-map ref could only name an *entity-level* pair (bond/overlay by constituents) — which `induce` already derives, and which (endpoints being unordered) can't even fix the atom bijection, only restate it. If written, `resolve` treats it as a **consistency assertion** against the induced correspondence (a contradicting one is an error), never a silent override. Nothing to build beyond the uniform ref grammar.

### 7.5 The primitive set and its operations (point iii)

Six objects:

- **`MoleculeAst`** `M` — a graph state.
- **`MoleculeCorrespondence`** `C` — a partial bijection `M ⇀ M` (+ induced entity maps), valueless. Atom slice = a graph-core `Correspondence` (accessor `atom_correspondence()`).
- **`Deltas`** `D` — the difference: self-contained, invertible, canonical.
- **`ReactionAst`** `Rx` = `(lhs: M, deltas: D)` — the lhs-anchored operational store.
- **`ReactionSpanAst`** `Sp` — the hub / generalized CGR: the superimposed union graph with `EntitySpan` columns.
- **`ReactionDerivation<'a>`** `Der` = `(host: &'a M, product: M, comap: C)` — a production fired once; `apply`'s codomain (§7.7 #5). Host borrowed, product/comap owned.

Names: `superimpose` / `difference_to` are proposed (2026-07-02, §7.5 review); `atom_map` was a placeholder for the atom `Correspondence`.

| op | signature | status |
|---|---|---|
| induce | `(M, M, atom: Correspondence) → C` | new (I5/I7 kernel) |
| `superimpose` (was `align`) | `Sp::superimpose(L: M, R: M, C) → Sp` | new — lift correspondence + values → span |
| `difference_to` (was `diff`) | `L: M::difference_to(R: &M, &C) → D` | new = `superimpose(L,R,C).to_reaction().deltas()`; receiver = source (lhs), fixed by the apply law |
| `from_sides` (was `from_lhs_rhs_map`) | `Rx::from_sides(L: M, R: M, atom: Correspondence) → Rx` | new = `Rx::new(L, L.difference_to(R, induce(L,R,atom)))` — SMIRKS import |
| replay | `Rx::to_reaction_span() → Sp` | done — store → hub |
| read-off | `Sp::to_reaction() → Rx` | done (I2f/I6c) — hub → store |
| project | `Sp::left()/right() → M` | done |
| recover corr | `Sp::correspondence() → C` | new — forget values (K interface) |
| apply | `Rx::apply(host, match) → Der` | change (was `→ M`; the comap is `Der`'s second half) |
| reverse | `Rx::reverse() → Rx` = `Sp` swap sides → read-off | (131) |
| compose (rule∘rule) | `Rx::compose(&Rx) → impl Iterator<Rx>` | (131) — via the hub; the *reserved* `compose` |
| inverse / compose | `D::inverse()`, `D::compose` (`D∘D`) | done |
| derivation accessors | `Der::product() → M`, `Der::comap() → C` | new |
| derivation abstract | `Der::to_reaction() → Rx` (= `difference_to` under `comap`) | new |
| correspondence compose | `C::compose(&C) → C` | new — relational `A⇌B` ∘ `B⇌C`; substrate for `Der::chain`; receiver disambiguates from `Rx::compose` |
| derivation chain | `Der::chain(Der) → Der` (comaps `compose` = pathway atom-map) | new |
| derivation reverse | `Der::reverse() → Der` | new |

`Sp` is the pivot `apply` and `difference_to` route through, and `Sp ↔ (L, R, C)` is a two-way bridge — `superimpose` builds `Sp` from a pair + correspondence; `correspondence()`+`left()`+`right()` read them back. `C` sits on **both sides** of the algebra: consumed by `superimpose`/`difference_to`, produced by `induce`, `Sp::correspondence()`, and `apply` (the `Der` comap). The two import paths share `induce`: `induce`→`difference_to`→`Rx` (deltas store) and `induce`→`superimpose`→`Sp` (CGR / SMIRKS working form).

**Laws** (diff ⊣ apply): `apply(L, L.diff(R, C)) ≡ R` under `C`; `L.diff(L, identity) = ∅`; `L.diff(R, C).inverse() ≡ R.diff(L, C⁻¹)` (= `reverse`); `diff = canonicalize ∘ raw`, so `L.diff(apply(L, d)) = canonicalize(d)`, identity only up to normalization.

**Semantic factoring (where the seams are).** Denotationally — independent of how anything is stored — the five are **four kinds** (DPO / patch-theory):

- **State** — `MoleculeAst`; an object.
- **Correspondence** — `MoleculeCorrespondence`; a valueless partial iso (`K↪L`, `K↪R`), a *relation*.
- **Patch** — `Deltas`; a self-contained, invertible arrow carrying `old→new`, its implicit domain only its *read-set*.
- **Production** — a matchable rule `L←K→R` with `L ⊇` read-set (carries preserved application-context).

`ReactionAst` and `ReactionSpanAst` **denote the same kind — a production** — in two presentations: operational-oriented `(lhs, deltas)` vs the symmetric attributed span (the CGR). Mutually determined (the lossless hub round-trip, 131); every difference — orientation, `K` explicit-vs-derived, which native op is cheap (apply vs compose/CGR) — is *representational*.

The three distinctions among the other kinds are semantic and firm: **Correspondence vs Patch** (valueless relation vs valued arrow); **Patch vs Production** (`Deltas ⊊ Production` — a production carries an arbitrary `L ⊇` read-set as application-context and is *matched* into a host, not applied in a known frame); **State** is the base object.

This **keeps all five** — the `ReactionAst`/`ReactionSpanAst` split is justified operationally (131: apply-native store vs compose/CGR working form) — but the seam is that that split is representational, not semantic; if operational needs ever changed, demoting `ReactionSpanAst` to a non-public transient is where the pressure would legitimately land. If *direction* is taken as semantic (a reaction *means* reactants→products), `ReactionAst` = *oriented* production and `ReactionSpanAst` = *symmetric* production — but even then one is the other decorated with a side, not an independent kind.

**The other direction — what's missing is `apply`'s discarded output.** A correspondence sits on *both* sides of the algebra: supplied to `diff` (`diff(L, R, C) → Deltas`), and *produced* by `apply` (`apply(P, host, match) → product` **plus the host↔product comap** — preserved atoms mated, deleted left-exposed, created right-exposed). The set makes the input side first-class (`MoleculeCorrespondence` feeds `diff`) but **drops the output side**: `apply` returns a bare `MoleculeAst`, discarding the comap it just built. That comap *is* the per-step atom map; 131 recomputes atom maps lazily "from (educts, rule, match)", but they are produced *during* `apply` (like MCES's `edge_count` during search) — re-deriving is the same waste.

So the missing representation is a **`Derivation` = `(product: MoleculeAst, comap: MoleculeCorrespondence)`** — the DPO direct derivation `G⇒H` reduced to its externally-useful data. It is the **instance** of a `Production` (rule : derivation ∷ function : one evaluation): a production is a matchable schema, a derivation is one concrete firing carrying its *ground-truth* map (its value over a post-hoc `diff(host, product)` is that `apply` *knows* the exact map — it created the atoms — where a later diff would have to reconstruct the correspondence). Its data is a *composition* of two existing kinds (State + Correspondence), but it denotes a role the four kinds don't name (the applied-rule 2-cell). It round-trips the rule layer: `ReactionAst::apply → Derivation` instantiates, `Derivation::to_reaction → ReactionAst` abstracts back (host as lhs, `diff(host, product, comap)`).

Compositions that are real but belong **one layer up** (the reaction-network / DG layer — built *from* the five, not peers): the **hyperedge** `(educts, rule, match)` (the compact form a `Derivation` expands from), the **pathway** (ordered productions, comaps composed = multi-step atom-map propagation), **molecule collections** / **rule systems**. Not missing from the algebra; the network on top of it.

Confirmed not-missing: match = injective `Correspondence`; CGR = `ReactionSpanAst`; multi-component = disconnected `MoleculeAst`; pattern = lattice `MoleculeAst`; reverse = `Production` reverse; symmetry = graph-core `Automorphism` (§7.6).

### 7.6 Kernel vs aggregation — where automorphism fits

The §7.1 producers share one primitive — a structure-preserving node map, lifted to entities by `find_by_participants` (§7.3) — but diverge in how node maps are **aggregated** and whether the result carries a transformation:

| aggregation | producer | object | transformation? |
|---|---|---|---|
| one partial bijection | MCS / enumeration / diff | correspondence → span | diff only (`align`) |
| a set of embeddings | subiso | `Vec<embedding>` | no |
| a group | automorphism (`auto.rs`, nauty) | `Automorphism` (generators, orbits, canonical labeling, order) | no |

**Automorphism** is the group aggregation over `L = L`. A single element is the endomorphic, total, structure-preserving special case of a correspondence (`Total` + endo + iso), and its `Vec<NodeId>` image is the `(i, perm[i])` pair form — but `Automorphism` deliberately stores the *group* (generators + orbits), not a bag of elements; folding it into the correspondence carrier would discard that. It stays its own type. It **reuses the kernel** — the `find_by_participants` induction lifts an atom permutation to a consistent bond/overlay permutation, needed the day symmetry must act on overlays (stereo coset frames, a full-molecule canonical key; `auto.rs` is node/edge-only today) — but **not** `align`/diff: a structure-preserving map has ∅ self-diff, and the §7.2 mated/exposed membership degenerates (one id space, all mated).

Shared model: **the node-map + induction kernel is common; the correspondence carrier, the embedding set, and the automorphism group are three aggregations over it, and only the correspondence carries a transformation (the span).**

### 7.7 Open decisions

1. ~~graph-core carrier~~ **settled** (§7.3): `CommonSubgraph` retires → `Correspondence` (its own module; subiso / MCIS / maximal / complete enumeration emit it) + `Mcs { correspondence, edge_count }` (MCES only, caching its objective). `edge_count`/`node_count` are derived `Correspondence` reads; kept lightweight for the enumeration path; `Matching` donates vocabulary only. Membership matching-theoretic mated/exposed (§7.2).
2. ~~`align` / import names~~ **settled** (§7.5): `superimpose` (span constructor), `difference_to` (diff, receiver = source), `from_sides` (import), `Correspondence::compose`, `Der::chain`.
3. ~~`MoleculeCorrespondence` surface syntax~~ **settled** (§7.4): atom mated-pair vec under `:atoms` (`int | keyword` refs); `MoleculeCorrespondenceInput` → `resolve(L, R)` (keyword→id, infer exposed, induce entities); structural-ref pairs allowed-but-tautological (consistency-checked).
4. ~~`MoleculeEmbedding` retirement~~ **settled**: yes — it retires into the injective `MoleculeCorrespondence` view (subiso emits a graph-core `Correspondence`; `verify_overlays` becomes the entity-map induction). No standalone `MoleculeEmbedding` struct.
5. ~~`ReactionDerivation`~~ **settled** — added; `apply`'s codomain is `ReactionDerivation`, not a bare product (the comap is computed during `apply` regardless — 131's lazy-map storage is orthogonal).
   - **name — `ReactionDerivation`**, no `Ast` (the `Ast` discriminator is **lattice structure**, not DSL-presence — `MoleculeAst` is lattice-valued, `Deltas`/`MoleculeCorrespondence`/`ReactionDerivation` are not; which of the six keep `Ast` is a separate audit pass). *Not* `ProductDerivation` (misweights the product half).
   - **shape — borrowed** `ReactionDerivation<'a> { host: &'a MoleculeAst, product: MoleculeAst, comap: MoleculeCorrespondence }` (the `MoleculeEmbedding<'a>` precedent). Only the host borrows; `product`/`comap` are owned, so `reverse`/`to_reaction` are self-contained while the borrow lives, and persisting past the host is a drop to owned `(product, comap)`. Lifetime-bound — fits the transient-right-after-`apply` pattern; the network stores the `(educts, rule, match)` hyperedge and recomputes, so it never holds a `ReactionDerivation` long-term.
   - **API beyond accessors** (`product()`/`comap()`/`atom_map()`): `chain` (`G⇒H` ∘ `H⇒J` → `G⇒J`, comaps `compose` = pathway atom-map propagation), `reverse` (`H⇒G`, comap inverse), `to_reaction`/`to_reaction_span` (abstract back to the rule layer).

### 7.8 Naming principle (settled)

**`Name = [primary-noun]? + concept-word [+ Ast if lattice-valued]`.** The qualifier appears **iff** the bare concept-word is generic or collides (with another layer, or another umol domain); it is the concept's **primary noun** — what the type is fundamentally about — never decorative or for family-grouping.

The six under the rule:

| type | why the qualifier (or not) |
|---|---|
| `MoleculeAst`, `ReactionAst` | concept-word *is* the noun; `Ast` = lattice-valued (§7.7 #5) |
| `MoleculeCorrespondence` | `Correspondence` collides with graph-core `Correspondence` and is generic → `Molecule` (a rule-free relation *between molecules*) |
| `ReactionSpanAst` | bare `Span` collides broadly → `Reaction`; `Ast` per lattice |
| `ReactionDerivation` | `Derivation` is generic and collides with spectroscopic *derivatives* → `Reaction` (the *firing of a rule*); non-lattice, no `Ast` |
| `Deltas` | specific plural, field-faithful (`Constraints` precedent), and *not* reaction-specific (`diff` returns it bare) → **bare** |

Corollary: `Correspondence`→`Molecule` but `Derivation`→`Reaction` because their primary nouns differ (a relation is about its operands; an instance is about its rule) — one rule, different outputs, not an inconsistency. `ReactionDeltas` declined: family-grouping is not a reason to qualify, and it would mislabel a fundamental cross-layer type as reaction-specific.

Two seams deliberately *not* name-encoded: `span ⊃ correspondence` is a derivation, expressed by `Sp.correspondence()` (they are different *kinds* — valueless relation vs valued production — not two sizes of one thing); `Deltas` is off the categorical *Patch* vocabulary but on the codebase's own (`Delta` enum, `diff`→`Deltas`).

### 7.9 Staged implementation plan

Modules top-down: **graph-core** (foundation) → **umol-ast `ast`** → **umol-ast `dsl`** (surface). Additive subitems keep the tree green; only breaking ones (signature / return-type change, type retirement) go red, and every stage ends green.

**S0 — graph-core foundation** (additive) **Done**
- **S0a** graph-core: `find_by_participants` on `FixedRelationSet` / `VarRelationSet` / `FixedVarBirelationSet` — per-factor (`&[P]` / `&[P1], &[P2]`), unordered-multiset via sort-both-compare (§7.3.1). `[dep: —]` **Done**
- **S0b** graph-core: `Correspondence` module — struct + `mates` / `left_exposed` / `right_exposed` / `node_count` / `shared_edge_count` / induced-edge reads (§7.2, §7.3.2). `[dep: —]` **Done**

**S1 — `MoleculeCorrespondence` + constructors** (additive) **Done**
- **S1a** ast: `MoleculeCorrespondence` struct (8 families) + `atom_correspondence()` + reads. `[dep: S0b]` **Done**
- **S1b** ast: `induce(L, R, atom: Correspondence)` — entity-map induction via `find_by_participants`. `[dep: S1a, S0a, S0b]` **Done**
- **S1c** ast: `from_maps`. `[dep: S1a]` **Skipped**
- **S1d** ast: `MoleculeCorrespondence::compose`. `[dep: S1a]` **Done**
- **S1e** graph-core: `Correspondence::is_total()` predicate (`mate_count == left_count == right_count`) instead of a `Total(_)` wrapper — the fast path branches on `atoms().is_total()`; union-frame skip is a later optimization. `[dep: S0b]` **Done**

**S2 — diff / span / import ops** (additive) **Done**
- **S2a** ast: `ReactionSpanAst::superimpose(L, R, C)` — build the union-frame span from a pair + correspondence. `[dep: S1a]` **Done**
- **S2b** ast: `ReactionSpanAst::correspondence()` — recover `C` (forget values). `[dep: S1a]` **Done**
- **S2c** ast: `MoleculeAst::difference_to(&R, &C)` = `superimpose` → `to_reaction` → `deltas`. `[dep: S2a]` **Done** — inherent method in reaction_span.rs (not molecule.rs, avoiding a foundational-module upward import; same split-impl pattern as substructure/symmetry/incidence).
- **S2d** ast: `ReactionAst::from_sides(L, R, atom)` = `induce` + `difference_to` + `new`. `[dep: S1b, S2c]` **Done**
- **S2e** ast (property test, feature `proptest`): cross-validate the two span constructions — `superimpose(L, R, C)` (Strategy A, direct) **==** the delta path (`to_reaction_span` of the deltas). Off a generated reaction span take `L = lhs`, `R = right()`, `C = correspondence()`, and assert `superimpose(L, R, C) == span`. A mismatch flags a `diff`-completeness or frame gap — the whole point of building A independently (two paths, assert equality; cf. doc 135, where testing an unrelated method surfaced a real compose gap). `[dep: S2a, S2b]` **Done** — plain + overlay variants (`test_reaction_span_ast_superimpose_matches_delta_path{,_overlay}`); both paths agree on all generated reactions, no gap surfaced.

**S3 — `ReactionDerivation` + `apply` codomain**
- **S3a** ast: `ReactionDerivation<'a>` struct + `lhs` / `rhs` / `comap` / `atom_map` + `to_reaction` / `reverse` / `chain`. Additive. `[dep: S1a, S1d, S2c]` **Done** — new `reaction_derivation.rs`; needed `Correspondence::reverse` (graph-core) + `MoleculeCorrespondence::reverse` (comap inversion).
- **S3b** ast: change `ReactionAst::apply` / `apply_at` codomain `M → ReactionDerivation`; migrate callers to `.rhs()`. **red→green**. `[dep: S3a]` **Done** — `apply_at` builds the `lhs↔rhs` comap during transform (host→product atom map by survivor-rank, then `induce` for bonds/overlays); all `apply`/`apply_at` callers read `.rhs()`.

**Naming resolution (during S3, 2026-07-02): the two reaction sides are `lhs` / `rhs` throughout.** `ReactionDerivation` fields/accessors are `lhs` (borrowed) / `rhs` (owned) — *not* `host` / `product` (a graph-nomenclature + chemistry chimera; chemical `educt`/`product` was declined too since a side is one *graph* that may hold several species, so species-plural names are awkward). Carried to `ReactionSpanAst::lhs()` / `rhs()` (was `left()`/`right()`), `EntitySpan` `Modified { lhs, rhs }` + `lhs()`/`rhs()` (was `left`/`right`), `ConstraintSpan`, `EntityPatch::diff*` args, and the `superimpose`/`difference_to`/`induce`/`from_sides` molecule args. Graph-core `Correspondence` keeps `left`/`right` (domain-neutral: it serves subiso/MCS/enumeration, not just reactions).

**S4 — retire `MoleculeEmbedding`** **Done**
- **S4a** ast: `MoleculeCorrespondence` is the one injective-view type — a subiso match and an induced subgraph are the *same* object (host↔sub injective map), not two things. So both `verify_overlays`/`substructure_matches` **and** `induced_subgraph` return `MoleculeCorrespondence`; `apply_at(&self, host, &MoleculeCorrespondence)` (host threaded explicitly — the correspondence carries no host ref; `host_atom(id)` → `atoms().right_of(…).expect(total-on-pattern)`); `extract` / `edits` become ops over (host, correspondence) on `MoleculeAst` (they only need the atom subset — `extract` ignores the overlay maps today). Drop the `MoleculeEmbedding` struct; migrate `substructure.rs` + `reaction.rs` + `molecule.rs` (`induced_subgraph`) + tests. **red→green**. `[dep: S1b, S3b]` **Done** — added `Correspondence::from_images(images, right_count)` to graph-core (dense-left constructor a match/embedding induces); `induced_subgraph`/`substructure_matches`/`verify_overlays` all build it; `MoleculeAst::{extract,edits}(&sub)`; `apply_at(host, &corr)`. Scope was larger than "test-only": umol-graph's `hmo.rs`/`kekulizer.rs`/`fingerprint/pattern.rs` also consumed it — migrated their `host_atoms`/`host_bonds`/`sub_atom` reads to the native `Correspondence` API (`.atoms().mates()` / `.left_of()`), no `host_*` accessors added to the general type.
- **S4b** graph-core (sequence after S4a lands): same unification one layer down, but graph-core's subgraph is node **and** edge (`edge_induced_subgraph` selects an explicit edge set — chords excluded — that a bare `Correspondence<NodeId>` can't hold; consumers read `host_edge()` too). So the mirror is a **two-family bundle** `GraphCorrespondence { nodes: Correspondence<NodeId>, edges: Correspondence<EdgeId> }` — the graph-core base that `MoleculeCorrespondence` (atoms + bonds + 6 overlays) extends. `induced_subgraph` (edges = all among the node set) and `edge_induced_subgraph` (edges = the explicit set) both return it; `Graph::extract(&self, sub) -> Graph`; `host_node`/`host_edge` reads become `sub.nodes()/.edges()` `right_of`/`left_of`. Retire `Embedding`; migrate `cycles.rs` (graph-core), `ring.rs` (umol-ast), `fingerprint/substructure.rs` (umol-graph). Graph-core `Correspondence` stays `left`/`right`. **red→green**. `[dep: S4a]` **Done** — `GraphCorrespondence { nodes, edges }` + `Graph::extract(&self, sub)`; both producers return it; consumers read `.nodes()/.edges()` `mates`/`right_of`/`left_of`. This is also the S6 target (see S6a): every graph-core subgraph producer converges on `GraphCorrespondence`.

**S5 — DSL surface** (additive) **Dropped 2026-07-03** — a `MoleculeCorrespondence` never serializes standalone. It is a *computational intermediary* (produced by `induce` / `substructure_matches` / `ReactionSpanAst::correspondence()`; consumed by `superimpose` / `difference_to` / `apply_at` / `from_sides`), same serde-free status as graph-core `Correspondence` / `GraphCorrespondence`. Its round-trip surface **is the span**: `superimpose(L, R, C) → ReactionSpanAst` fuses both molecules + the pairing + values, and `ReactionSpanDsl` (I4) already serializes it losslessly (`lhs()`/`rhs()`/`correspondence()` recover the three parts). Foreign-format import (SMIRKS/GML map numbers) feeds `induce(L, R, atom_map)` programmatically inside the umol-io reader — the map lives in the *foreign* syntax, not a native `:correspondence` EDN. A standalone `:correspondence` form would only add a redundant second way to author a both-sides reaction (which entities are `Unchanged`/`Modified`/`Added`/`Removed` **is** the pairing) — against not proliferating surfaces.
- ~~**S5a** dsl: `MoleculeCorrespondenceInput { atoms }` + `FromEdn` / `ToEdn` + `resolve(L, R)`; `:correspondence {:atoms […]}`.~~ Dropped (see above).

**S6 — graph-core producer unification** **Done**
- **S6a** graph-core: subgraph producers converge on the two carriers. **Done 2026-07-03.** MCIS / MCES / maximal / complete→`GraphCorrespondence` (the search builds it directly — no interim `CommonSubgraph` twin, node/edge counts come from the two graphs it holds); this **subsumes the `Mces { correspondence, edge_count }` wrapper** — each objective is a family `mate_count()` (MCES = `edges().mate_count()`, MCIS = `nodes().mate_count()`), no separate `edge_count`. **subiso→`Vec<Correspondence<NodeId>>`** (node-only; edges are derived on demand via `GraphCorrespondence::induced(left, right, nodes)`, which runs `edge_mates` — kept off the hot enumeration path). Added `GraphCorrespondence::induced`; retired `CommonSubgraph`, `Mces`, the three singular `maximum_common_*_subgraph` methods (a plural + `.into_iter().next()` covers them), and the `Enumerate` enum (the search always enumerates all maxima). ast side: per-family `induced_{bonds,dative_bonds,aromatic_systems,multicenter_bonds,noncovalent_bonds}` in `correspondence.rs` feed both `MoleculeCorrespondence::induce` and `substructure::verify_overlays` (now taking `Correspondence<NodeId>`, layering the predicate + pattern-total gate over the inducer pairing, stereo bespoke coset); the `GraphView` adapter keeps its `Vec<Vec<AtomId>>` public contract, adapted internally. Migrated callers: `common_subgraph.rs`/`subiso.rs`/`lib.rs` (graph-core), `substructure.rs`/`compose.rs`/`view/graph.rs` (umol-ast). `[dep: S0b, S4b; sequence after S4]`

Critical path **S0 → S1 → S2 → S3 → S4**. S5 rides after S1 (parallel to S2–S4). S6 is graph-core-only; place last.
