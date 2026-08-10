# Maximum common subgraph: API and McGregor implementation plan

Status: McGregor implemented in `umol-graph-core/src/algorithms/mcs.rs` (MCIS + MCES +
seeded MCES, both connectivity modes, both enumerate modes). Future section
(RASCAL / ModularProductClique, the `maximum_independent_sets` primitive, seeded MCIS,
ITS atom-mapping) outstanding.

## Scope

New module `umol-graph-core/src/algorithms/mcs.rs` providing maximum common subgraph
between two `Graph`s, in two flavors sharing one result type:

- **MCIS** — maximum common *induced* subgraph (maximize mapped vertices; an edge is
  shared iff present in both).
- **MCES** — maximum common *edge* subgraph (maximize shared edges; missing edges
  allowed, classical isomorphic edge subgraph).

Both take caller `node_match`/`edge_match` predicates (the chemistry seam, identical in
spirit to `subgraph_isomorphisms`), `Connected`/`Disconnected` connectivity, and
`First`/`AllMaximum` enumeration. MCES additionally has a seeded entry point with a
forced `anchor` correspondence and a warm-start `hint`.

One algorithm is implemented now: **McGregor backtracking** (vertex-mapping variant),
the general solver covering the whole config space — the role VF2 plays for `subiso` and
branch-and-bound plays for `mis`. The algorithm enums therefore ship single-variant
(`McGregor`), matching `SubgraphIsomorphismAlgorithm::Vf2`,
`MaxIndependentSetAlgorithm::BranchAndBound`, `AutomorphismAlgorithm::Nauty`.

### Out of scope (named here so the enums have room to grow, not built now)

- **RASCAL** (MCES via max-weight clique on the line/edge product) and
  **ModularProductClique** (MCIS via max clique on the modular product) — faster routes,
  added as new enum variants when built. Both reuse/extend `mis.rs` (clique = MIS on the
  product complement); note `maximum_independent_set` returns a single maximum, so the
  clique route naturally serves `First` and would need an all-maximal-cliques pass for
  `AllMaximum` — another reason McGregor (which enumerates natively) is the first cut.
- **Graded / order-changing commonality** (a 2→1 bond change counted as a partial match):
  rejected for MCS — it is not classical MCES. The minimal-chemical-distance / ITS model
  belongs to a separate future algorithm built *on top* of MCES (seed + edit cost), not
  crammed in here.
- **Seeded MCIS** — the `anchor`/`hint` plumbing is kind-agnostic, so an
  `..._induced_subgraph_seeded` is a one-line wrapper if wanted; not exposed now.
- Chemistry-level predicates (element / bond-order match) live in the caller
  (`umol-graph`), not here.

## API

Singular vs plural is encoded by the **method name**, not a flag — matching
`subgraph_isomorphism(s)` / `maximum_independent_set(s)`. Singular returns
`CommonSubgraph` (a maximum always exists, empty at worst — so no `Option`), plural
returns `Vec<CommonSubgraph>` (all maxima). There is no `McsEnumerate`; the singular
path uses the tie-pruning bound (`bound <= best`), the plural keeps ties
(`bound < best`). With only `connectivity` left as a knob there is no `McsConfig` —
`McsConnectivity` is passed directly (matching `subiso`/`mis`, which take no config).

```rust
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum McsConnectivity { Connected, Disconnected }

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum McisAlgorithm { McGregor }   // + ModularProductClique later

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum McesAlgorithm { McGregor }   // + Rascal later

/// One common subgraph: vertex correspondence and number of shared edges.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CommonSubgraph {
    mapping: Vec<(NodeId, NodeId)>, // (self node, other node), sorted by self node
    edge_count: usize,
}
impl CommonSubgraph {
    pub fn mapping(&self) -> &[(NodeId, NodeId)];
    pub fn node_count(&self) -> usize; // mapping.len()
    pub fn edge_count(&self) -> usize; // shared edges under the mapping
    pub fn is_empty(&self) -> bool;
}

impl Graph {
    // induced (MCIS)
    pub fn maximum_common_induced_subgraph(..)  -> CommonSubgraph;       // one maximum
    pub fn maximum_common_induced_subgraphs(..) -> Vec<CommonSubgraph>;  // all maxima

    // edge (MCES)
    pub fn maximum_common_edge_subgraph(..)  -> CommonSubgraph;
    pub fn maximum_common_edge_subgraphs(..) -> Vec<CommonSubgraph>;

    // seeded edge (anchor forced + never skipped; hint warm-starts the bound)
    pub fn maximum_common_edge_subgraph_seeded(..)  -> CommonSubgraph;
    pub fn maximum_common_edge_subgraphs_seeded(..) -> Vec<CommonSubgraph>;
}
// every method: (&self, other, node_match, edge_match, [anchor, hint,] connectivity, alg)
```

The plural is the combinatorial path (symmetric scaffolds and degenerate matches
multiply the count); prefer the singular unless ties are needed.

## Algorithm: McGregor backtracking

Reference: J. J. McGregor, "Backtrack search algorithms and the maximal common
subgraph problem", *Software: Practice and Experience* 12 (1982) 23–34. Vertex-mapping
formulation (extend a partial vertex correspondence, allow vertices to be skipped,
branch-and-bound on the objective) — the variant used in the RASCAL paper's baseline.

State (parallel to `Vf2State`):
- `mapping: Vec<Option<u32>>` over `self` nodes → `other` node, `reverse: Vec<bool>` over
  `other`.
- `current` correspondence list; `best: Vec<CommonSubgraph>` incumbent(s).
- For `Connected`: a frontier of `self` nodes adjacent to the mapped set (after the first
  pair, candidates are drawn only from the frontier).

Search step:
1. **Bound / prune.** MCIS objective is vertex count: prune if
   `mapped + remaining_self_unprocessed ≤ best_size`. MCES objective is edge count: prune
   if `current_edges + upper_bound_remaining_edges ≤ best_edges`, where the bound counts
   edges incident to still-mappable vertices (a loose degree bound; RASCAL's tighter
   bound is the later optimization).
2. **Pick** the next unprocessed `self` node (frontier-restricted under `Connected`).
3. **Branch — map** it to each unused `other` node `v` with `node_match(u, v)` true and
   the pair *consistent* (below); recurse.
4. **Branch — skip** it (the MCS degree of freedom VF2 lacks): leave `u` unmapped,
   recurse. Skipping is forbidden for `anchor` nodes.
5. At a leaf (all `self` nodes processed) compare the realized objective to `best`:
   strictly better replaces; equal appends iff `AllMaximum`; under `First` an `≥` prune in
   step 1 keeps a single incumbent.

Consistency of adding `(u → v)` given already-mapped pairs `(u' → v')`:
- **MCIS (induced):** for every mapped `u'`, `edge(u,u')` present in `self` ⟺
  `edge(v,v')` present in `other`; when both present, `edge_match` must hold on the two
  edge ids.
- **MCES (edge):** only the forward direction — when `edge(u,u')` present in `self` *and*
  `edge(v,v')` present in `other`, `edge_match` must hold; absence on either side is
  allowed (subgraph need not be induced). `edge_count` increments per matched edge.

Edges via `Graph::find_edge` (used the same way in `Vf2State::feasible`).

Seeding:
- **anchor:** pre-place each pair (set `mapping`/`reverse`, seed the connected frontier
  from their neighbors, as `Vf2State::seed_anchor` does); reject up front if any pair
  fails `node_match` or the pairs are mutually inconsistent; anchor nodes are never
  skipped or reassigned.
- **hint:** compute the objective of the largest consistent sub-correspondence within the
  hint and install it as the initial incumbent `best`, so branch-and-bound prunes from a
  good lower bound. Purely a warm start — never constrains the result.

## Files

- `umol-graph-core/src/algorithms/mcs.rs` — new module (types, methods, `McGregor` impl,
  tests).
- `umol-graph-core/src/algorithms.rs` — `pub mod mcs;` between `matching` and `mis`.
- `umol-graph-core/src/lib.rs` — re-export `CommonSubgraph, McesAlgorithm, McisAlgorithm,
  McsConnectivity` (between the `matching` and `mis` lines).

## Tests (`mcs.rs`, `#[rstest]` table tests, specific mappings/counts asserted)

- **MCIS, identical small graphs** — triangle vs triangle: `node_count == 3`,
  `edge_count == 3`, one mapping under `First`.
- **MCIS, induced strictness** — path P3 vs triangle K3: max induced common is an edge
  (`node_count == 2`), not the open path embedded in the triangle (which would be a
  non-induced match), confirming the induced consistency check.
- **MCES vs MCIS divergence** — P4 vs C4 (or the P3/K3 pair): MCES recovers more shared
  edges than MCIS recovers as an induced match; assert the differing `edge_count`.
- **Connectivity** — two graphs whose common part is two disjoint edges:
  `Disconnected` gives `node_count == 4`/`edge_count == 2`; `Connected` gives a single
  edge.
- **Enumerate** — symmetric case (e.g. P3 self-MCS) where `First` returns one mapping and
  `AllMaximum` returns the full symmetric set; assert exact mapping lists.
- **node_match / edge_match filters** — a labeling that forbids one correspondence shrinks
  the result to the asserted mapping (mirrors the `subiso` filter cases).
- **Seeded MCES** — `anchor` forces a pair into the result (assert it is present in every
  returned mapping; an anchor that excludes the global optimum yields the constrained
  optimum); `hint` produces the same maxima as the unseeded call (assert mapping-set
  equality — hint changes only search, not answer); empty `anchor`+`hint` equals the
  plain `maximum_common_edge_subgraph`.
- **Empty / oversize** — empty `other` → `vec![]` of an empty subgraph or `vec![]`
  (pick one and pin it); single isolated nodes → vertex-only maxima.

## Verification

- `cargo test -p umol-graph-core --lib mcs`
- `cargo clippy -p umol-graph-core --lib`
- `cargo build -p umol-graph` (downstream still builds with the new re-exports).

## Future (separate go-aheads)

1. `ModularProductClique` / `Rascal` variants over the product graph + `mis.rs`, with
   RASCAL's edge bound.

   The clique route's `AllMaximum` needs an enumeration `mis.rs` does not yet provide.
   Note the distinction: **maximum** ISs/cliques = all sets of the *largest* size (a
   handful — what MCS `AllMaximum` wants, i.e. all common subgraphs of the optimal size);
   **maximal** ISs/cliques = all *non-extendable* sets of any size (Bron–Kerbosch,
   exponential superset, Moon–Moser 3^(n/3)) — needed only for a Koch-style *connected*
   MCS-via-cliques formulation, not for plan 107. So the clique route needs all-maximum,
   not all-maximal; Bron–Kerbosch is a separate heavier primitive, build only if that
   formulation is adopted.

   Decided API for the all-maximum primitive (additive, nothing existing changes):

   ```rust
   pub fn maximum_independent_set(&self, alg) -> Vec<NodeId>;        // keep: single, cheap, common path
   pub fn maximum_independent_sets(&self, alg) -> Vec<Vec<NodeId>>;  // new: all maxima
   ```

   - New sibling method, not a changed return type — the single-result case is dominant
     (a clique-route `First` calls it); plural is opt-in like `subgraph_isomorphisms`.
   - Same `BranchAndBound` variant: the B&B already walks the whole tree; collecting all
     size-maxima is the same `First`/`AllMaximum` leaf change (keep ties when
     `current.len() == best.len()` instead of pruning). No new algorithm ⇒ no new variant.
   - Lives in `mis.rs` (general primitive); MCS calls it on the product complement.
   - An `MisEnumerate` param on one method was rejected: the return type must change with
     it (`Vec<NodeId>` vs `Vec<Vec<NodeId>>`), so it cannot be one signature.
   - Output: each inner set `sort_unstable`'d (as today); outer `Vec` sorted for
     determinism.
2. Seeded MCIS wrapper if a use case appears.
3. Minimal-chemical-distance / ITS atom-mapping algorithm consuming MCES (seed) + a bond
   edit cost — the graded model deliberately excluded from MCS.
