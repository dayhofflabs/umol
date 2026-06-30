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

**I0 — graph-core: remapping & participant positions** `[—]`
- I0a — graph-core: `ParticipantPosition(u32)` newtype. `[—]` **Done** (relation.rs, beside `RelationId`; re-exported from lib.rs).
- I0b — graph-core: rename `Remapping`→`RemovalRemapping`; relation-set `apply_remapping`→`apply_removal_remapping` on all five types (no reindex — monotonic); migrate callers (`MoleculeBuilder`, the `IdRemapping` wrapper); move data-column compaction to free `remove_node_vec`/`remove_edge_vec`. `[—]` **Done** (79 type + 22 method renames; `apply_to_*_vec`→free fns over `map_node`/`map_edge`; `IdRemapping`/`UndoRemapping` untouched; graph-core 296 + umol-ast 3914 tests green).
- I0c — graph-core: new general `Remapping` — total relabel (`Vec<NodeId>`/`Vec<EdgeId>`), `map_node`/`map_edge`/`unmap_*`/constructors. `[—]` **Done** (`new`/`map_node`/`map_edge` + `#[rstest]` tests; graph-core 301 green). `unmap_*` **deferred** — partial inverse for a non-bijective injection with no current consumer (forward covers apply/remap_delta/reindex; `reverse()` builds a fresh forward map; rollback uses `Edit`/`Undo`); add when a consumer appears.
- I0d — graph-core: `apply_remapping(&Remapping) -> (Self, …σ)` on all five relation-set types (argsort → per-relation participant permutation). `[I0a, I0c]` **Done.** Per-factor σ (non-conflated): single-factor sets return `(Self, Vec<ParticipantPosition>)`, birelations `(Self, Vec<ParticipantPosition>, Vec<ParticipantPosition>)`. Added `RelationParticipant::remap(&Remapping) -> Self` (renamed the removal pair → `remap_removal`/`unmap_removal`) and `FactorOrdering::canonicalize_positions` (argsort returning σ; Ordered = identity). graph-core 306 + umol-ast 3914 green.
- I0e — ~~graph-core: common-subgraph enumeration over the incidence (Levi) graph~~ **Void — folded into I3b.** No graph-core work: `enumerate_common_subgraphs` is already generic over any `Graph` + predicates, and umol-ast already builds the Levi graph with overlay pseudonodes (`MoleculeAst::incidence_graph`, `IncidenceNodeSelection::{topological,constitution,full}`, and `substructure_matches_incidence`). "Overlay-aware overlap" is just compose calling the existing enumeration over `incidence_graph(constitution())` with kind-matching predicates — a umol-ast/compose change (I3b), not graph-core.

**I1 — overlay reindex + the four uniform overlay deltas** `[I0]`
- I1a — umol-ast/ast: `AromaticSystemAst::permute` / `MulticenterBondAst::permute` — reorder `electrons` by a `&[ParticipantPosition]`, both delegating to the shared `ElectronCountsAst::permute` (`Undetermined` unchanged; charge/spin/constraints positionless). `[I0a]` **Done** (only the two positional families have it — no blanket needed; umol-ast 3918 green).
- I1b — delta.rs: the four non-stereo overlay `*Delta` (`DativeBondDelta`, `AromaticSystemDelta`, `MulticenterBondDelta`, `NoncovalentBondDelta`) on the `EntityDelta` pattern (`Add`/`ModifyField`/`ModifyConstraint`/`Remove`, atom-set structural payload, stable id); their `Atoms` types + `into_delta`; extend the `Delta` sum. `[—]`
- I1c — delta.rs: `Deltas::canonicalize` over the four families (generic `EntityDelta` fold + `field_ops!`; mechanical). `[I1b]`
- I1d — delta.rs: `remap_delta` over overlay deltas — re-anchor participants through the general `Remapping`, reindex positional payloads via I1a. `[I1a, I1b, I0c]`

**I2 — apply + span for the four overlays** `[I1]`
- I2a — reaction.rs: `apply_at` overlay lowering arms (overlay `Delta` → existing overlay `Edit`); fold in item 1 (molecule-constraint lowering via `Constraint::map_topology_refs`). `[I1b]`
- I2b — reaction.rs: DPO dangling check extended to overlay incidence (a deleted atom's overlay participations, not only `bond_ids()`). `[I1b]`
- I2c — reaction_span.rs: lift the six overlay relation-set `data` columns to `EntitySpan<…>` (all six incl. stereo, uniform). `[I1b]`
- I2d — reaction_span.rs: `to_reaction_span` folds overlay deltas onto `lhs` overlays (an `apply_*_change` per family). `[I2c, I1b]`
- I2e — reaction_span.rs: `left()`/`right()` carry unchanged overlays through `from_parts`. `[I2c]`
- I2f — reaction_span.rs: `to_reaction` (span→deltas) for overlays (`EntityDelta::diff`/`deltas_from_states` defaults). `[I2c, I1b]`

**I3 — composition for overlays** `[I1, I2, I0e]`
- I3a — compose.rs: unify `remap_delta`/frame algebra onto the general `Remapping`; extend the four-class composite frame to overlay relation ids. `[I1d, I0c, I0d]`
- I3b — compose.rs: overlay-aware overlap — call the existing generic `enumerate_common_subgraphs` over `MoleculeAst::incidence_graph(constitution())` (both sides) with overlay-pseudonode kind-matching predicates, so overlay-only connectivity is visible and the overlay correspondence falls out of the Levi overlap. `[I1, I2 — both incidence builder + enumeration already exist]`
- I3c — compose.rs: admissibility (boundary-bond / combined-frame dangling) extended over overlay incidence. `[I3a, I2b]`

**I4 — overlay span DSL** `[I1b, I2c]`
- I4a — dsl/molecule.rs: factor the six inline overlay renderers (`render_dative` …) into shared `render_<entity>_entry` (the existing TODO). `[—]`
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
