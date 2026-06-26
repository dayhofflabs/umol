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

## W2 — umol-ast: resolved delta + reaction AST **Done**

The resolved-delta vocabulary is **molecule-level, not reaction-scoped** — it is the
`Delta` counterpart of the deferred `Edit` in `ast/edit.rs`, reusable for base+delta
molecule storage and MMP (the reason `ReactionEdit` was rejected in doc 131). So it lives
at **`ast/delta.rs`, a sibling of `ast/edit.rs`**; only the reaction-specific types
(`ReactionAst`, `ReactionSpanAst`) use it. The doc-127 interim `ast/reaction.rs` is
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
2. **`Deltas` newtype + `Canonicalize` — done** (`Vec<Delta>`, flat container). The #2
   reduction is **generic over a `DeltaFamily` trait** (the created/preserved state machine
   written once) with the uniform per-variant field ops (`fuse_field`/`field_is_identity`/
   `apply_field`) generated by a `field_ops!` **macro** from the `(variant ⇒ ast field)` map —
   so atoms, bonds, and the six increment-2 overlay families share one fold. Per entity: fuse
   `SetField` chains per field (continuity-checked), `Add` absorbs sets, `Add`+`Remove`
   cancels, `Remove` subsumes prior sets and carries the *original* value (changes reverted on
   the removed ast); entity `SetConstraint` folds as a keyed old→new chain by `(entity,
   constraint-key)`. Molecule-level `ConstraintDelta` is **multiset**-folded (net multiplicity
   per constraint, `Add`↔`Remove` cancel one-for-one — *not* dedup; corrects doc 131). Then
   the cross-entity dangling check (an added bond's endpoint must not be a net-removed atom),
   then a stable `sort()` (the new `Ord` on the delta enums) — order is **not** stored.
   `Err(Contradiction)` on the doc-131 conditions. Confluence checked by an idempotence
   proptest. (Flat-vs-split container is the one deferred data-structure decision — a
   non-breaking later refactor since both shapes assemble over the same per-family enums.)
3. **`ReactionAst { lhs: MoleculeAst, deltas: Deltas }` — done** (`reaction.rs`). The type +
   `impl Canonicalize`. `canonicalize` reduces `deltas` (the #2 reduction) and **passes `lhs`
   through**: `MoleculeAst` has structural `PartialEq`/`Eq` (cache-excluded) but no whole-molecule
   `Canonicalize`, so there is no canonical form to compose — `ReactionAst` derives `PartialEq`
   and reduces only the deltas. Equality here is **value-level in a fixed atom frame** (`lhs` +
   delta canonical forms); equality *up to atom renumbering* (reaction iso-dedup) is a separate
   `umol-graph` operation via `canonical_form` on the condensed graph — `umol-ast` cannot reach
   the full incidence, so it does not attempt iso-canonicalization (forced by layering).
   `right()`/`left()` are projections of the span (item 4, value-accumulation — no lowering);
   `reverse()` is done (invert each delta + re-anchor to the product frame, item 4); only `apply`
   onto a host (items 5–6) needs the delta→`Edit` lowering + the `transact` engine. The interim
   doc-127 `ReactionAst`/`Assignment`/`Stereo{Atom,Bond}Correspondence`, `RewriteError`, and
   `ast/molecule/rewrite.rs::apply_rule` (no production callers, superseded by item 5 `apply`) are
   **removed**. `umol-graph/fingerprint/reaction.rs` is left red until it migrates to `right()`
   (W4).
4. **`ReactionSpanAst` + `to_reaction_span()` — done** (`reaction_span.rs`). The materialized
   superimposed `L ∪_K R` graph (the DPO rule span; deliberately **not** "condensed" — it is an
   *expansion* of, not a condensation of, the operational `ReactionAst`). A graph-core `Graph`
   over the `lhs` frame (deleted elements kept as nodes/edges, created elements appended) +
   parallel `Vec<Change<AtomAst>>` / `Vec<Change<BondAst>>`, where `Change = Unchanged | Modified
   { left, right } | Added | Removed` carries DPO membership + the per-side value(s). Built by
   **value-accumulation**: canonicalize the deltas, annotate each `lhs` element, and compute a
   `Modified` element's right value by applying its field/constraint changes at the value level
   (reusing the canonicalize fold's apply via `pub(crate)` `apply_atom_change` / `apply_bond_change`
   in `delta.rs`) — **no RHS materialized, no `transact`, no back-map**. `left()` / `right()`
   project a side back to a `MoleculeAst` (survivors renumbered); the DPO span `L ←K─ R` is the
   tag partition. **Scope:** localized topology only — molecule-level constraints and overlays are
   not represented here (carried by the operational form). **No exporters** (CGR / MOD): those are
   `umol-io` boundary types via `FromAst`, out of scope; GML in/out is out of scope.
   **`reverse()` — done** (same module): invert each delta (`Delta::inverse`) and re-anchor its
   ids / bond endpoints to the product's compacted frame (survivors take `right()`'s frame;
   deleted elements become created and take fresh ids). `reverse().to_reaction_span()` swaps the
   span's sides. No lowering — pure id remap over the canonical deltas.
5. **`apply_at(&self, m: &MoleculeEmbedding) -> Result<MoleculeAst, ApplyError>` — done**
   (`reaction.rs`). Reuses `MoleculeEmbedding` as the match (`host_atom` / `host_bond` maps).
   Canonicalize `deltas`; **DPO gluing check** (a deleted host atom keeps no localized bond the
   rule does not also delete → `ApplyError::Dangling`); lower the match-remapped deltas to a
   `Vec<Edit>` in emit order **AddAtoms → AddBonds → Set\* → RemoveTopology (last)**, created
   atoms taking `New(0..k)` (flat counter) and preserved/removed referencing `Id(host …)`; then
   **checked** `transact` on `host.edit()` (precondition holds for concrete reactions — the match
   forces `host == old`; a genuine mismatch surfaces as `ApplyError::Transaction`). Molecule-level
   constraints deferred. `ApplyError` is flat: `Dangling { host_atom } | Inconsistent |
   Transaction(_)`.
6. **`apply(&self, host, subiso) -> impl Iterator<MoleculeAst>` — done** (`reaction.rs`).
   `lhs.substructure_matches(host, GraphAndOverlays, subiso)` × `apply_at`; the `subiso`
   algorithm is an explicit argument (transparency). DPO-invalid (dangling) matches are skipped
   (`filter_map`).

Tests: `Deltas::canonicalize` table cases (each fuse/cancel/contradiction row from doc 131) +
idempotence (`canonicalize∘canonicalize == canonicalize`) via proptest; `to_reaction_span` /
`reverse` round-trip; `apply` on concrete localized-only reactions (bond make, bond break,
order change, atom add/remove) asserting the exact product `MoleculeAst`; dangling rejection.

## W3 — umol-ast: minimal compose **Done**

Implemented in `compose.rs`: `ReactionAst::compose(&self, other, scope) -> Vec<ReactionAst>` +
`CompositionScope { RcAnchored, Full }`; the frame remap reuses `pub(crate) delta::remap_delta`
(also now backing `reverse`). Both gluing conditions are enforced (boundary-bond /
pushout-complement rejection during `lhs_C` construction; combined-frame dangling for B's
deletion of shared atoms). Tests: exact composite, `apply`-equivalence, an A-created overlap
(create-then-modify fusion across the seam), and a boundary-bond rejection (`N→N-C` ∘ `N-C-O`
yields none).

`ReactionAst::compose(&self, other, scope) -> Vec<ReactionAst>` (`compose.rs`).
Sequential composition A;B: applying the composite equals `B.apply(A.apply(H))`. A = `self`,
B = `other`.

**Overlap enumeration.** `R_A = self.to_reaction_span().right()`, `L_B = other.lhs`. Enumerate
maximal common *induced* subgraphs of `(R_A, L_B)` via the W1 `enumerate_common_subgraphs`
(`BronKerbosch`); compatibility is **symmetric** — `Lattice::meet(…).is_some()` on the
`AtomAst` / `BondAst` pairs (not the asymmetric `matches`). Each overlap `E` is the matter A
produces and B consumes, as a node correspondence `R_A ↔ L_B` (`CommonSubgraph::mapping`).

**Scope.** `CompositionScope::RcAnchored` (default) keeps overlaps touching A's reaction center
(the changed elements of `self.deltas`, projected onto `R_A`); `Full` keeps all (the free
rule-algebra sum, for algebra / CRN work).

**Composite frame.** `lhs_C` and `deltas_C` share one frame with four id classes, allocated in
order: (1) `lhs_A` atoms `0..n_A`; (2) `L_B \ E` atoms appended `n_A..n_A+e`; (3) A-created
atoms (from `deltas_A`), shifted to `≥ n_A+e`; (4) B-created atoms after those. Maps:

- **R_A → composite**: read off A's span survivors — R_A id `k` is the k-th survivor of A's span
  (`lhs_A` non-removed in place, A-created appended) → its `lhs_A` id (class 1) or its shifted
  A-created id (class 3).
- **L_B → composite**: an overlap atom → its R_A node → composite (above); a non-overlap atom →
  its class-2 id; a B-created atom → class 4.

`deltas_C = shift_created(deltas_A) ++ remap(deltas_B)`, then `canonicalize` (#2 — fires the
create-then-delete cancellation and keep-then-modify fusion across the A|B seam). `lhs_C =
lhs_A + (L_B \ E)` with the incident L_B context bonds.

**Admissibility = the DPO gluing conditions (genuine, not heuristics).** An overlap is rejected
iff the sequential composite provably does not exist:

- **Pushout-complement (boundary bonds).** An `L_B` context bond incident to an overlap atom
  that A *created* and to a non-overlap atom cannot sit in `lhs_C` (its created endpoint does not
  exist before the composite runs). No `H` realizes such a match — after A that bond is absent
  (A did not create it; it was not in `H`), so B cannot match. Reject. When the overlap atom is
  `lhs_A`-preserved the bond *does* sit in `lhs_C` → admissible.
- **Dangling.** B deleting a shared atom that retains an undeleted incident bond in the combined
  frame → reject (the `apply` gluing check, on the combined frame).

Each rejected overlap has no `B.apply(A.apply(H))` witness; the conditions are silent in the
common (admissible) case.

Tests: compose two concrete localized bonding rules; assert
`compose(A,B).apply(H) == B.apply(A.apply(H))` on hand-built `H` (the oracle — exercises
admissible composites and confirms rejected overlaps have no witness); RC-anchored vs `Full`
result sets; admissibility rejection; associativity on a small `(A,B,C)` triple (inherited from
#2, verified empirically).

## W4 — retire interim + migrate callers **Done**

- doc-127 correspondence macro/types, `Assignment`, `RewriteError`, and
  `ast/molecule/rewrite.rs` (`apply_rule`) removed in W2.3.
- `umol-graph/fingerprint/reaction.rs` migrated: the product side is derived via
  `to_reaction_span().right()` (was the stored `rhs`); the three tests build
  `ReactionAst::new(lhs, deltas)`. `FingerprintError::Inconsistent` added for the
  derive-failure path. Workspace green (`cargo build`/`clippy --all-targets`); no interim
  references remain.

## Out of scope (increment 2 / follow-on)

- Overlay + stereo delta families; overlay composition (#7); stereo `TransformFrame` (#8);
  saturation for iv (#9); the per-entity-split `Deltas` container; canonical encoding of
  overlays into the incidence graph.
- Exporters (SMIRKS/GML/CGR) as `umol-io` boundary types and `ReactionDsl` (EDN) — they read
  `to_reaction_span`; separate workstream, not blocking the core.
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
