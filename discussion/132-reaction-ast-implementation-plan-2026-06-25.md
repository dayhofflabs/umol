# 132 — Reaction AST: implementation plan (increment 1)

Implements doc 131 for **localized topology**: `ReactionAst = MoleculeAst (lhs) + Deltas`,
`apply`, minimal `compose`, the supporting graph-core primitives, and retirement of the
doc-127 interim. Overlays + stereo + overlay composition (#7–#9) are **increment 2**,
out of scope here. Molecules keep their overlays in `lhs` (unchanged, carried through);
rules in this increment only add/remove/modify atoms, bonds, and constraints.

Terminology: *topology* = atoms + **localized** bonds (single, double, triple); *overlays*
are the **non-localized** entities (aromatic systems, dative, multicenter, noncovalent) plus
stereo. The diff between topology and overlay is localization — that is the increment 1 /
increment 2 boundary.

All design decisions are settled in doc 131 (§"Algorithm decisions and open points"):
DPO + injective matches (#1); the delta reduction system (#2); created-entity numbering via
canonical labeling (#3 ⊂ #6); overlap enumeration (#4); composite = glue+concat+canonicalize
(#5); canonical labeling = nauty on the incidence graph (#6).

## Verified layering (placement is sound)

- `umol-ast` depends on `umol-graph-core` → `apply`/`compose` can call the generic graph-core
  primitives (subiso, common-subgraph, canonical labeling) on a `Graph` view of the AST, with
  an `AtomAst`/`BondAst` compatibility predicate. No `umol-graph` dependency needed.
- `nauty-Traces-sys` is already a `umol-graph-core` dep (used in `algorithms/auto.rs`).
- `transact`/`Undo` live in `umol-ast` (`ast/molecule/transact.rs`), so `apply`'s lowering
  path is in-crate.
- `algorithms/mcs.rs` (maximum only: `McisAlgorithm`/`McesAlgorithm`/`McGregor`) is **merged
  into** a new `algorithms/common_subgraph.rs` that also holds the enumeration; callers use
  the crate-root re-exports, so the move is path-only (`pub mod`/`pub use` + any in-crate
  `algorithms::mcs` references). `algorithms/subiso.rs` (`SubgraphIsomorphismAlgorithm`) is
  reused for apply.

## W1 — graph-core primitives **Done**

New/extended in `umol-graph-core/src/algorithms/`:

1. **Common-subgraph module** (`common_subgraph.rs`, merging the current `mcs.rs`). Both
   operations on the shared modular-product foundation, sharing the `CommonSubgraph` result:
   - **maximum** (moved from `mcs.rs`): `Mcis`/`Mces`, variant `McGregor`, `McsConnectivity`,
     methods `maximum_common_induced_subgraph[s]` / `…_edge_subgraph[s]`.
   - **enumeration** (new): `CommonSubgraphEnumerationAlgorithm` (cf.
     `MatchingEnumerationAlgorithm`), variant `BronKerbosch` (full-name, not `Bk`), method
     `enumerate_common_subgraphs(…, alg) -> Vec<CommonSubgraph>` — modular product +
     Bron–Kerbosch maximal-clique enumeration → every maximal common *induced* subgraph
     (disconnected included), under caller `node_match`/`edge_match` predicates. **No
     connectivity parameter** — disconnected is both required and complete for overlaps
     (a connected `L_B` can still match `R_A` in disjoint regions, e.g. bimolecular
     consumption of two A-produced fragments; connected-only would drop those), and rule
     sides are small so there is no pruning motivation. `McsConnectivity` therefore stays
     MCS-only (the maximum path); the enumeration never needs it.
   Operation names (maximum vs enumeration) disambiguate — no maximum/maximal collision. The
   merge is path-only (callers use the crate-root re-exports).
2. **Canonical labeling + generic incidence** (extend `auto.rs`). nauty's `densenauty`
   already yields the canonical labeling of a vertex-colored graph; expose
   `canonical_form(&Graph, coloring) -> Remapping`. Add a **generic incidence transform**
   (edge-color→vertex-color) covering **topology only** — atoms + localized bonds (single,
   double, triple) → a vertex-colored graph. *By design this is the topology portion only.*
   The **full ABDAMNSS incidence** (overlays, multicenter/dative/noncovalent, stereo as
   colored structure) stays in `umol-graph`, the only layer with the full entity accounting,
   and calls graph-core's `canonical_form`. For the localized-topology reactions of this
   increment the generic transform suffices; the full incidence is increment 2.
   Consumer note: `canonical_form` serves reaction iso-dedup (equality *up to atom
   renumbering*), which is a `umol-graph` generation concern (follow-on). `umol-ast` does
   **not** call nauty in this increment — see W2.3.
3. **Directed graph + deterministic topological sort.** Directed graphs are a graph-core
   foundation type, not a domain-side or per-caller construction. New `DiGraph` (general
   directed graph, CSR out-adjacency, `digraph.rs`, reusing `NodeId`) — one model serves
   the delta-dependency DAG now and cyclic reaction/derivation networks later. New
   `algorithms/toposort.rs`: `impl DiGraph::topological_order<K: Ord>(key, alg)` via Kahn's
   algorithm draining a **key-ordered ready set** (ties by `NodeId`) → a unique order;
   `None` on a cycle (so the "DAG" case needs no separate type). Used to sequence deltas
   when lowering.

Tests: small hand-built graphs with enumerated expected overlaps / canonical labels /
orderings; assert exact result sets, not counts.

## W2 — umol-ast: resolved delta + reaction AST

The resolved-delta vocabulary is **molecule-level, not reaction-scoped** — it is the
`Delta` counterpart of the deferred `Edit` in `ast/edit.rs`, reusable for base+delta
molecule storage and MMP (the reason `ReactionEdit` was rejected in doc 131). So it lives
at **`ast/delta.rs`, a sibling of `ast/edit.rs`**; only the reaction-specific types
(`ReactionAst`, `CondensedReactionAst`) use it. The doc-127 interim `ast/reaction.rs` is
replaced in W4.

1. **Per-family delta enums — done** (`ast/delta.rs`, sibling to `ast/edit.rs`).
   Increment-1 families: `AtomDelta` and `BondDelta`, each with `Add`, `Remove`, `SetField`
   (reusing the existing `AtomFieldChange`/`BondFieldChange` payloads — so atom element and
   **bond order** changes ride here), and `SetConstraint { old, new }` for inline per-entity
   constraints (keyed old→new: `(None,Some)` add, `(Some,None)` remove, `(Some,Some)`
   modify — old/new sharing a key, which is what lets canonicalization chain-fuse and detect
   contradictions, like `SetField`). Identity is the **stable per-family id** — uniform
   across all families including overlays, where an atom set can't identify the entity
   (multicenter spans >2 atoms; multicenter/noncovalent allow parallels); `BondDelta`
   `Add`/`Remove` additionally carry `endpoints` as structural payload, not as identity.
   `ConstraintDelta { Add(Constraint), Remove(Constraint) }` for molecule-level constraints,
   with **multiset** semantics (the store is not force-deduped — dedup is only the lazy
   canonical form), so `Add`/`Remove` are genuine inverses. Sum:
   `enum Delta { Atom, Bond, Constraint }`; per-family `inverse()` (closure under inversion:
   `Add`↔`Remove`, `SetField`/`SetConstraint` swap old/new).
   *(Dative/aromatic/multicenter/noncovalent/stereo families are increment 2.)*
2. **`Deltas`** newtype (`Vec<Delta>`) — **flat container** for increment 1. `impl
   Canonicalize` = the #2 reduction system: group by target `(family,id)` / `(family,id,
   field)`, fold each (fuse chains, cancel `Add`/`Remove`, drop no-ops, `Remove` subsumes
   prior sets, `Add` absorbs sets); entity `SetConstraint` folds as a keyed old→new chain by
   `(entity, constraint-key)`; molecule-level `ConstraintDelta` cancels `Add`↔`Remove`
   one-for-one **preserving multiplicity** (multiset, *not* dedup — corrects doc 131's
   "duplicate `Add C` dedups"); then the cross-entity dangling check on the folded set;
   `Err(Contradiction)` on the conditions in doc 131. Order is **not** stored. (Per-entity split container is the one deferred data-structure decision — a
   non-breaking later refactor since both shapes assemble over the same per-family enums.)
3. **`ReactionAst { lhs: MoleculeAst, deltas: Deltas }`** (`reaction.rs`). `impl Canonicalize`
   composing `lhs`' and `deltas`' canonical forms; lazy `equiv`. Equality here is
   **value-level in a fixed atom frame** (`lhs` + delta canonical forms); equality *up to
   atom renumbering* (reaction iso-dedup) is a separate `umol-graph` operation via
   `canonical_form` on the condensed graph — `umol-ast` cannot reach the full incidence, so
   it does not attempt iso-canonicalization (forced by layering). `reverse()` (map `inverse`
   over deltas, re-anchor on the rhs). Migrate `ast/molecule/rewrite.rs` to this type.
4. **`CondensedReactionAst`** + **`to_condensed()`** (`condensed.rs`). Replay `deltas` on
   `lhs` → the superimposed union graph with per-element `(left,right)` values + atom map
   (= the generalized CGR). Membership derived from each pair. Convenience: `right()`,
   `membership()`, `atom_map()`.
5. **`apply_at(&self, host, match) -> Result<MoleculeAst, ApplyError>`**: remap `deltas` onto
   `match` (K ids → host ids, created → `New`) → a `Vec<Edit>` → existing `transact`. DPO
   **dangling check** before transacting (reject if a deleted host atom carries unmatched
   bonds); matches assumed injective.
6. **`apply(&self, host) -> impl Iterator<MoleculeAst>`**: enumerate injective matches of
   `lhs` into `host` via graph-core `subiso` (compatibility = `AtomAst`/`BondAst` lattice
   meet) over a `Graph` view of each (reuse `ast/embedding.rs`), `× apply_at`. Lazy.

Tests: `Deltas::canonicalize` table cases (each fuse/cancel/contradiction row from doc 131) +
idempotence (`canonicalize∘canonicalize == canonicalize`) via proptest; `to_condensed` /
`reverse` round-trip; `apply` on concrete localized-only reactions (bond make, bond break,
order change, atom add/remove) asserting the exact product `MoleculeAst`; dangling rejection.

## W3 — umol-ast: minimal compose

`ReactionAst::compose(&self, other) -> impl Iterator<ReactionAst>` (`compose.rs`):

1. `R_A = self.to_condensed().right()`; `L_B = other.lhs`.
2. Enumerate overlaps via the W1 all-maximal-common-subgraph primitive on `(R_A, L_B)`,
   compatibility = lattice meet.
3. **RC-anchored filter** (default): keep overlaps that touch `self`'s reaction center (the
   changed elements of `self.deltas` projected onto `R_A`). A `CompositionScope::Full` flag
   drops the filter (the free rule-algebra sum, for algebra/CRN work).
4. Per admissible overlap (DPO dangling check on the combined frame): **glue** (id-identify
   per the overlap map — assemble-along-a-partial-map, id reuse, no node merge) → **concat**
   `self.deltas ++ remapped other.deltas` → **`canonicalize`** (#2, which performs the
   create-then-delete cancellation and keep-then-modify fusion) → `ReactionAst(lhs_C,
   deltas_C)`, where `lhs_C = self.lhs + (L_B \ overlap)`.

Tests: compose two concrete localized bonding rules; assert
`compose(A,B).apply(H) == B.apply(A.apply(H))` on hand-built `H`; RC-anchored vs `Full`
result sets; admissibility rejection; associativity on a small `(A,B,C)` triple (inherited
from #2, verified empirically).

## W4 — retire interim + migrate callers

- Remove the doc-127 correspondence macro/types from `ast/reaction.rs`.
- Migrate `ast/molecule/rewrite.rs` and `umol-graph/fingerprint/reaction.rs` (reaction
  fingerprint) to the new `ReactionAst`.

## Out of scope (increment 2 / follow-on)

- Overlay + stereo delta families; overlay composition (#7); stereo `TransformFrame` (#8);
  saturation for iv (#9); the per-entity-split `Deltas` container; canonical encoding of
  overlays into the incidence graph.
- Exporters (SMIRKS/GML/CGR) as `umol-io` boundary types and `ReactionDsl` (EDN) — they read
  `to_condensed`; separate workstream, not blocking the core.
- Network-generation (iii) APIs over molecule collections → `umol-graph` ops, drive
  `apply`/`compose`; follow-on.

## Risks / verification

- **nauty canonical_form stability** — guard with fixed-input canonical-label tests.
- **BK enumeration blowup** — bounded in practice by RC-anchoring; `Full` mode is opt-in.
- **Incidence-encoding placement — decided.** Generic topology-only incidence
  (atoms + localized bonds) in graph-core; the full ABDAMNSS (non-localized) incidence in
  umol-graph. Not a risk, a layering decision.
- **AST→`Graph` view** — confirm `ast/embedding.rs` exposes what subiso/common-subgraph need
  (atoms→nodes, bonds→edges, compatibility hook); add the bridge if absent.
- **Confluence (#2)** — covered by the idempotence proptest, not assumed.
- Whole workspace green + `cargo clippy --workspace --all-targets`.
