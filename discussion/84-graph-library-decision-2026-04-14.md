# Graph Library Decision

## Context

Step 7 of the unified AST migration (doc 80) requires ring enumeration, aromaticity perception, and neighbor iteration over `MoleculeAst`. Planning this step revealed that `MoleculeAst` has no native graph topology — `bonds: Vec<BondTuple>` is a flat list requiring O(|bonds|) scans for every neighbor query.

The indirection analysis in doc 82 ("Indirection analysis" section) traced the data path for four algorithms:

| Algorithm | Current indirection | With native topology |
|---|---|---|
| Ring enumeration | 4–5 conversions (petgraph → HashMap → BTreeMap → DenseProjection) | 0–1 (induced subgraph if filtering) |
| VF2 subgraph iso | 1 (build petgraph::Graph from flat bonds) | 1 (build petgraph) or 0 (native VF2) |
| Morgan/WL | 1 (build dense adjacency from flat bonds) | 0 (borrow slices) |
| DPO reactions | N/A | 0 for reads, new graph for output |

Every algorithm builds its own adjacency from scratch. `MorganTarget` already implements CSR internally. `AtomAdjacency` + `DenseProjection` exist solely to bridge petgraph's linked-list representation to the dense arrays the algorithms actually need.

## Three options considered

### Option 1: CSR in MoleculeAst

Add CSR (offsets + neighbors arrays) to `MoleculeAst`, built from `Vec<BondTuple>` at parse time. Bond attributes in parallel `IndexVec<BondIdx, BondAst>`.

- Fast to implement
- Good read performance, eliminates most indirection
- Immutable topology — attribute mutation (solver narrowing) works, but topology mutation does not
- Creates an asymmetry: mutable predicates, immutable topology
- Essentially "another table molecule with more pieces"
- No further optimization path

### Option 2: Keep petgraph

Continue using `StableGraph` in `MoleculeBuilder` and building `petgraph::Graph` in the matcher.

- `StableGraph` supports deletion but petgraph's VF2 requires `Graph` (not `StableGraph`) — can't have stable indices AND VF2 in the same representation
- Workaround: clone `StableGraph` → `Graph` for VF2
- Every graph algorithm extracts adjacency from petgraph into HashMap/BTreeMap/Vec<Vec<>> — petgraph stores adjacency internally but its linked-list representation is incompatible with what algorithms need
- Wall looming: any topology mutation conflicts with VF2's requirements
- petgraph contributes exactly one algorithm we use: vanilla VF2 (~2.7x slower than RDKit's VF2+)

### Option 3: New graph library crate

Build `Graph<N, E>` as a separate crate with native adjacency, stable indices, and algorithmic primitives.

- Significant upfront cost (~1 week)
- Solves the fundamental problem: one topology, all algorithms borrow from it
- Mutation + cache efficiency
- Can implement VF2+ directly, removing the petgraph algorithm dependency
- Similar precedent: EDN parser was a 1-week detour that became a solid foundation

## Decision: Option 3

New crate `umol-graph-core` (name TBD). Current `umol-graph` becomes `umol-discrete` (name TBD) — the discrete molecular structure layer (AST, solver, perception, matching, fingerprints).

## Mutation model

The argument for mutable topology was examined:

- **Solver**: never mutates topology. Narrows attributes only.
- **DPO reactions**: by definition produces a new graph. Input is immutable. The pushout constructs output from interface graph + RHS.
- **Reaction network construction**: append-mostly. Species and reactions are added; deletion is rare.

None of these require in-place topology mutation in the strong sense. But petgraph's `StableGraph`/`Graph` split (stable indices OR VF2, not both) is a real wall. And append-only growth is a legitimate mutation operation that CSR handles poorly.

The primary motivation for option 3 is not mutation per se but eliminating the indirection tax: one shared topology that all algorithms can borrow from without conversion.

## Design constraints

- `Vec<Vec<(NodeId, EdgeId)>>` adjacency (Graphs.jl style), not LEDA half-edge lists. Max degree ≤ 6 for molecular graphs; `Vec::retain` on the neighbor list is faster than any linked structure at bounded degree.
- Stable indices via free lists (reuse deleted slots).
- Do not target billion-node graphs. Molecular graphs (< 200 atoms) and reaction networks (< 10^6 species) are both in-memory with `Vec<Vec<>>`. Billion-node processing is a fundamentally different problem.
- Generic `Graph<N, E>` — node type N, edge type E. Not molecule-specific. Molecule-specific algorithms (Morgan, aromaticity, ring perception) stay in `umol-discrete`; generic graph algorithms (BFS, DFS, connected components, biconnected components, VF2) live in `umol-graph-core`.

## Scope

| Component | Estimate | Notes |
|---|---|---|
| `Graph<N, E>` struct + mutation | ~500 lines | add/remove node/edge, neighbor iteration, stable indices |
| BFS, DFS | ~100 lines | adapt from existing `algorithms/` |
| Connected components | ~50 lines | adapt from existing `algorithms/` |
| Biconnected components | ~100 lines | adapt from existing `algorithms/` |
| VF2 | ~500 lines | vanilla first, VF2+ later |
| `MoleculeAst` migration | ~300 lines | replace `Vec<BondTuple>` with `Graph` topology |
| Ring enumeration on `Graph` | ~200 lines | eliminates `DenseProjection`, `AtomAdjacency` |
| Morgan on `Graph` | ~100 lines | borrow adjacency directly |
| **Total** | **~1850 lines** | |

## RelationSet: typed hyperedge collections

`MoleculeAst` has multiple relation types beyond localized bonds: dative bonds, noncovalent bonds, aromatic systems, multicenter bonds. These are not graph edges — they don't need adjacency traversal. But they reference atoms by `NodeId` and must stay consistent when nodes are removed.

`RelationSet<R>` is a typed collection of relations (hyperedges) over a shared `NodeId` space. Each relation has an ordered participant list (`Vec<NodeId>`) and typed data (`R`). Per-node incidence lists enable O(incident) cascade removal when a node is deleted.

### Mapping to molecular relations

| Relation | participants | R (data) |
|---|---|---|
| Dative bond | `[donor, acceptor]` (ordered) | `BondAst` |
| Noncovalent bond | `[source, target]` (ordered) | `NoncovalentBondAst` |
| Aromatic system | `[atom0, ..., atomN]` | `AromaticSystemAst` (charge, spin, electron count) |
| Multicenter bond | `[atom0, ..., atomN]` | `MulticenterBondAst` (charge, spin, electron count) |

Participants carry the topology (who is connected to whom). `R` carries the chemistry (attributes of the relation). This mirrors the `Graph<AtomAst, BondAst>` separation where adjacency is topology and `N`/`E` are attributes.

### MoleculeAst composition

```
MoleculeAst {
    graph: Graph<AtomAst, BondAst>,              // primary topology with adjacency
    dative_bonds: RelationSet<BondAst>,           // directed binary, no adjacency
    noncovalent_bonds: RelationSet<NoncovalentBondAst>,
    aromatic_systems: RelationSet<AromaticSystemAst>,
    multicenter_bonds: RelationSet<MulticenterBondAst>,
    constraints: Vec<MoleculeConstraint>,
}
```

### Why not parallel arrays

The previous design stored each relation type as a `Vec<BondTuple>` or `Vec<AromaticSystem>` with atom indices. Node removal required manual cascade cleanup across every vec — one `retain` call per relation type, in user code. `RelationSet` internalizes this: `remove_participant(node)` drains the incidence list for that node, removes each incident relation, and cleans co-participants' incidence entries. No user-side bookkeeping.

### Incidence queries

- `aromatic_systems.has_incident(atom)` — O(1), replaces the current O(|systems| x |atoms_per_system|) scan in `is_in_aromatic_system`
- `aromatic_systems.incident(atom)` — returns the `RelationId`s for all systems containing that atom
- `dative_bonds.incident(atom)` — all dative bonds involving that atom

### Current AST type changes needed

`AromaticSystem { atoms: Vec<AtomIdx> }` splits: participants absorbed by `RelationSet`, data becomes `AromaticSystemAst` (charge, spin, electron count — to be defined). Same for `MulticenterBond`. `BondTuple { source, target, bond }` splits: source/target become participants, bond data stays as `BondAst`.

## What this unblocks

- Step 7 (aromaticity perception on solver) — ring enumeration and neighbor iteration work natively on `MoleculeAst`
- Future VF2+ implementation — no petgraph API constraint
- DPO graph rewriting — read from `Graph`, construct output as new `Graph`
- WL iteration, other graph algorithms — borrow adjacency directly

## VF2 implementation notes

### Current: vanilla VF2 (Cordella et al., 2004)

Implemented in `umol-graph-core/src/algorithms/vf2.rs`. Operates on undirected `Graph` with caller-supplied `node_match` and `edge_match` closures. Same algorithm as petgraph's `subgraph_isomorphisms_iter`, without the directed-graph-with-doubled-edges workaround (native undirected graph halves the per-feasibility-check iteration cost).

### RDKit's VF2 differences

RDKit uses a modified VF2 ("VFLib") with domain-specific optimizations:

- **Fingerprint pre-screening**: Morgan-based quick-reject before entering VF2. Most non-matches never reach the search.
- **Degree-based candidate ordering**: highest-degree query node first (most constrained), narrowing the search tree. Closer to VF2+ than vanilla.
- **Chirality and ring membership as first-class feasibility predicates**: checked inside VF2, not as post-filters.
- **Lazy iteration**: stateful iterator returning one match at a time. Avoids collecting all matches for containment checks.

The reported ~2.7x gap vs petgraph is primarily pre-screening and ordering, not asymptotic.

### VF2+ (Jüttner & Madarasi, 2018)

Three additions over vanilla:

- **Static candidate ordering**: precomputed query node ordering based on domain sizes and connectivity, not dynamic smallest-index selection.
- **Label-frequency cutting rules**: compare neighbor label histograms at each candidate pair; prune if the target can't satisfy the query's distribution.
- **Backjumping**: skip backtracking levels that can't affect the current failure.

2-10x over vanilla VF2 on sparse molecular graphs. The current `Vf2State` struct is structured for this: `next_query_node` and the candidate loop are separate concerns that VF2+ would replace.

### Closure architecture: algorithm vs semantics boundary

The `subgraph_isomorphisms` function takes two closure parameters: `node_match` and `edge_match`. These are the stable boundary between topology-level search (inside the algorithm) and chemistry-level semantics (provided by the caller).

VFLib's molecule-specific heuristics split cleanly across this boundary:

**Inside the algorithm** (topology heuristics, modify `Vf2State` internals):
- Degree-based or domain-size-based candidate ordering (VF2+ `next_query_node`)
- Label-frequency cutting rules (VF2+ feasibility)
- Backjumping (VF2+ search control)

**Via the closures** (semantic predicates, caller-supplied):
- Element, charge, hydrogen matching
- Ring membership, ring count, ring size (via `DerivedCatalog`)
- Chirality
- Any future `DerivedPred`

**Pre-screening** is a third closure, evaluated once before VF2 enters the search tree. It takes the query and target as a whole and returns whether VF2 should run at all. Fingerprint-based quick-reject (Morgan feature subset check) is the primary use case. The signature:

```rust
pub fn subgraph_isomorphisms(
    query: &Graph,
    target: &Graph,
    node_match: &mut impl FnMut(NodeId, NodeId) -> bool,
    edge_match: &mut impl FnMut(EdgeId, EdgeId) -> bool,
    pre_screen: Option<&mut impl FnMut(&Graph, &Graph) -> bool>,
) -> Vec<Vec<usize>>
```

Or equivalently, the caller checks pre-screening before calling `subgraph_isomorphisms` — the algorithm doesn't need to know about it. Either way, pre-screening is caller logic, not algorithm logic.

The three concerns are orthogonal: vanilla VF2 with rich closures (current state + catalog pushdown), VF2+ with simple closures, VF2+ with rich closures + pre-screening. Any combination works without modifying the others.

### Bounded-treewidth alternatives

Molecular graphs have treewidth ≤ 3 (acyclics) to ≤ 6 (most ring systems). This is exploitable in principle:

- **Tree decomposition + DP**: subgraph iso on treewidth-*k* graphs in O(n^{k+1} · p). For *k* = 3-6, polynomial with manageable exponent. Linear-time decomposition exists (Bodlaender) but the constant is large.
- **Color-coding** (Alon, Yuster, Zwick): randomized O(2^p · m) for pattern size *p*. Competitive for p ≤ 15 (most SMARTS queries).

For single-molecule, single-query matching (< 200 atoms), VF2+ with pre-screening is hard to beat. Tree decomposition becomes relevant in two scenarios:

- **Many queries against one target**: amortize the decomposition cost. This is the reaction pattern matching case — many rewrite templates applied to the same molecule. Precomputing the target's tree decomposition once and running each pattern query as DP over it can outperform repeated VF2 calls.
- **Large graphs**: reaction networks (10^6+ nodes), protein interaction graphs.

## Derived property integration with VF2

### Problem

SMARTS predicates like `R3` (in exactly 3 SSSR rings) and `r6` (in a ring of size 6) are derived from topology, not intrinsic atom properties. They should not be materialized on `AtomAst`, which represents intrinsic properties. But they must be evaluated inside VF2's feasibility check for effective pruning — a `[r6]` query atom tried against every non-ring target atom wastes entire subtrees.

### Architecture: predicate pushdown with demand-driven materialization

Three layers, mirroring SQL index scans / Datalog IDB materialization / Cypher property indexes:

**Layer 1 — predicate specification (MoleculeAst)**

`AtomAst` stays clean. Ring constraints live on `MoleculeConstraint` as `DerivedPred` variants, scoped to specific atoms via `RelationRefs`:

```rust
DerivedPred::RingCount(ValueAst::Lit(3))    // R3
DerivedPred::InRingOfSize(6)                 // r6
```

This is the WHERE clause. The AST carries *what* to check, not *how*.

**Layer 2 — materialized indexes (MatchTarget)**

`MatchTarget` precomputes derived properties from graph topology into a catalog:

```rust
struct MatchTarget<'a> {
    ast: &'a MoleculeAst,
    catalog: DerivedCatalog,
}

struct DerivedCatalog {
    ring_info: Option<RingInfo>,
    // future: distance_matrix, bcc, ...
}

struct RingInfo {
    ring_count: Vec<u16>,                    // per atom: SSSR ring count
    ring_sizes: Vec<SmallVec<[u8; 4]>>,      // per atom: sorted ring sizes
}
```

Computed **on demand**: if no query references ring membership, `ring_info` stays `None`. This is the database deciding which indexes to build for a given query plan.

**Layer 3 — predicate pushdown (MatchQuery + node_match)**

`MatchQuery` compiles per-atom derived predicates from the query's `MoleculeConstraint` vec. During VF2, `node_match` checks both intrinsic properties and compiled derived predicates:

```
node_match(q_node, t_node):
    q_ast.atom(q_node).matches_ground(t_ast.atom(t_node))
    && query.derived_predicates(q_node)
           .all(|pred| target.catalog.satisfies(t_node, pred))
```

### Integration with VF2: no algorithm changes needed

VF2 already takes external `node_match` and `edge_match` closures. The algorithm is topology-only — all semantic matching is delegated to these closures. The current matcher already closes over `q_ast` and `t_ast` to check element/charge/etc. Adding derived predicate checks extends the closure body, not the algorithm:

```rust
// Current
let mut node_match = |q: NodeId, t: NodeId| {
    q_ast.atom(q.into()).matches_ground(t_ast.atom(t.into()))
};

// With derived predicates
let mut node_match = |q: NodeId, t: NodeId| {
    q_ast.atom(q.into()).matches_ground(t_ast.atom(t.into()))
        && query.check_derived(q, &target.catalog, t)
};
```

VF2 calls `node_match` at the top of its feasibility check, before edge consistency and look-ahead. Derived predicates evaluated here prune at the earliest possible point — a rejected candidate at depth 0 eliminates an entire subtree.

The same closure mechanism works for VF2+: the algorithm changes candidate ordering and cutting rules internally, but `node_match` and `edge_match` remain the external predicate interface.

### What pushes down vs what stays as post-filter

| Predicate | Scope | Evaluation |
|---|---|---|
| Element, charge, H-count | per atom, intrinsic | `node_match` (current) |
| Ring count, ring size | per atom, derived | `node_match` via catalog (pushdown) |
| Bond order | per edge, intrinsic | `edge_match` (current) |
| Dative/aromatic subset | multi-atom relational | post-filter (current) |
| Total charge, total spin | whole-match aggregate | post-filter |

Boundary: per-atom predicates (intrinsic or derived) push into `node_match`. Multi-atom relational predicates stay as post-filters on the complete assignment. This matches the database analogy — single-table predicates push into the scan operator, join predicates stay in the join operator.

## Implementation status

### Complete

| Component | Lines | Location |
|---|---|---|
| `Graph` struct (CSR, Arc/CoW) | ~330 | `umol-graph-core/src/graph.rs` |
| Mutations + `Remapping` | ~130 | `umol-graph-core/src/graph.rs` |
| `FixedRelationSet<R, N>` | ~120 | `umol-graph-core/src/relation.rs` |
| `VarRelationSet<R>` | ~120 | `umol-graph-core/src/relation.rs` |
| Connected components | ~35 | `umol-graph-core/src/algorithms/connected.rs` |
| Biconnected components | ~100 | `umol-graph-core/src/algorithms/bcc.rs` |
| Cycle enumeration | ~90 | `umol-graph-core/src/algorithms/cycles.rs` |
| Maximum independent set | ~70 | `umol-graph-core/src/algorithms/mis.rs` |
| VF2 subgraph isomorphism | ~170 | `umol-graph-core/src/algorithms/vf2.rs` |
| `UnionFind` | ~35 | `umol-graph-core/src/union_find.rs` |
| `MoleculeAst` migration | ~400 | `umol-graph/src/ast/molecule.rs` |
| Matcher on native VF2 | ~100 | `umol-graph/src/ast/matcher.rs` |

### petgraph removal from AST layer

petgraph removed from `matcher.rs` (replaced by native VF2), `morgan.rs` (petgraph benchmark variant deleted), `hueckel_rule.rs` (replaced by `UnionFind`). petgraph remains only in `graph_ir/` (`molecule.rs`, `molecule_builder.rs`), which is being superseded by MoleculeAst + solver.

### Remaining

- VF2+ candidate ordering and cutting rules
- Fingerprint pre-screening for matcher
- `graph_ir/` petgraph removal — not planned, GraphIR is going away
- Old standalone algorithms (`umol-graph/src/algorithms/{bcc,cycles,mis}.rs`) — kept for GraphIR consumers, removed when GraphIR is removed

## Relationship to other docs

- Doc 80: step 7 suspended pending this work
- Doc 82: indirection analysis motivating this decision
- Doc 83: step 7 plan summary and suspension note
