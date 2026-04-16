# 82 — Graph Representation for Molecular Workloads

## Context

During Morgan fingerprint optimization (doc 80, step 5), CSR views and petgraph were benchmarked for ECFP computation over `MoleculeAst`. The results prompted a review of whether petgraph is the right graph substrate for umol.

## Benchmark evidence

### ECFP4 (radius 2), 9,120 molecules, pre-built views

| Implementation | Time | us/mol | vs RDKit |
|---|---|---|---|
| direct (MoleculeAst) | 128 ms | 14.0 | 1.1x faster |
| MorganTarget (CSR) | 98 ms | 10.7 | 1.4x faster |
| petgraph (Graph) | 53 ms | 5.8 | 2.5x faster |
| MorganTargetOpt (CSR + u64) | 55 ms | 6.0 | 2.5x faster |
| RDKit (C++) | 132 ms | 14.7 | baseline |

### ECFP6 (radius 3), 9,120 molecules, pre-built views

| Implementation | Time | us/mol | vs RDKit |
|---|---|---|---|
| petgraph (Graph) | 73 ms | 8.0 | 2.3x faster |
| MorganTargetOpt (CSR + u64) | 74 ms | 8.1 | 2.3x faster |
| RDKit (C++) | 166 ms | 18.5 | baseline |

### Key finding

petgraph and CSR views are indistinguishable for molecular-sized graphs. At < 200 atoms, the entire adjacency fits in L1 cache regardless of memory layout. The bottleneck is hashing and duplicate removal, not neighbor iteration.

## petgraph assessment

petgraph is the default Rust graph library. Benchmarks show it is sufficient for molecular workloads.

**Representation:**
- `Graph` uses a linked-list-of-edges adjacency structure (inherited from Boost BGL). Cache-hostile in theory, but irrelevant for molecular-sized graphs where the working set fits in L1.
- `StableGraph` adds indirection for stable indices after removal.
- No first-class CSR representation.
- Edge weights stored inline in the edge array, mixing topology and data.

**What it provides:**
- VF2 subgraph isomorphism (used in the matcher, doc 80 step 4)
- BFS/DFS, connected components, toposort
- Correct, tested implementations

**What it costs:**
- API friction: `NodeIndex`/`EdgeIndex` wrapper types, `Directed`/`Undirected` type parameter
- No measurable performance cost for molecular graphs (< 200 atoms)

## Better designs in the literature

- **Compressed Sparse Row (CSR/CSC)**: optimal for static graphs with iteration-heavy workloads. One contiguous slice per vertex for neighbor access. This is what the Morgan benchmark validated.
- **LEDA / GraphTool (C++)**: separate topology from attributes. GraphTool uses contiguous arrays internally (relies on Boost for VF2). https://leda.uni-trier.de /, https://graph-tool.skewed.de
- **Ligra / GBBS (MIT, Shun & Blelloch)**: parallel graph algorithms over compressed sparse representations. Edge-map/vertex-map primitives. https://github.com/ParAlg/gbbs
- **Graphs.jl (Julia)**: clean interface/storage separation. `SimpleGraph` is `Vector{Vector{Int}}` with no abstraction tax.
- **ECS pattern (game engines)**: topology in one array, atom data in parallel array, bond data in parallel array. Same index space. This is effectively what `MorganTargetOpt` does.
- **BYO**: benchmark suite, not an impl. https://github.com/wheatman/BYO

## What molecular graphs need

- Static after construction (no insert/delete)
- Small (< 200 atoms, < 1000 bonds typically)
- Iteration-heavy (neighbor traversal, BFS, subgraph matching)
- Multiple edge types (localized bonds, dative, noncovalent) — ideally separate topologies, not a tagged union per edge
- Attribute access by index (atom/bond properties in parallel arrays)

## VF2 substructure matching benchmark

Three query patterns over 9,120 molecules, element-only atom matching, any-bond matching:

| Pattern | umol (ms) | RDKit (ms) | Ratio | Hits |
|---|---|---|---|---|
| branched (5 atoms, C(C)C(C)N) | 131 | 51 | 2.6x slower | 3979 / 3990 |
| phenol (7 atoms, 6-ring + O) | 239 | 82 | 2.9x slower | 2071 / 2065 |
| bicyclic (9 atoms, fused 5-6) | 137 | 51 | 2.7x slower | 171 / 175 |

Hit count differences are from RDKit rejecting ~150 molecules during sanitization.

petgraph's VF2 is ~2.7x slower than RDKit's C++ VF2. This is within "proceed and optimize later" territory.

### CSR VF2 re-run (2026-04-15)

After the AST migration: bonds moved into a CSR `Graph` in `MoleculeAst`, VF2 switched from petgraph to the native implementation in `umol-graph-core` operating directly on CSR neighbor slices. Same three patterns, same 9,120-molecule corpus.

| Pattern | CSR VF2 (ms) | petgraph VF2 (ms) | Speedup | RDKit (ms) | Hits (CSR) |
|---|---|---|---|---|---|
| branched | 108 | 131 | 1.21x | 51 | 4000 |
| phenol | 194 | 239 | 1.23x | 82 | 2076 |
| bicyclic | 111 | 137 | 1.23x | 51 | 178 |

~20% uniform speedup from removing the petgraph build step and the `Directed` double-edge overhead. Gap to RDKit narrows from 2.6–2.9x to 2.1–2.4x; the remaining gap is algorithmic (vanilla VF2 vs VF2+).

Hit counts shifted up by 21/5/7. Unrelated to the CSR switch — the `AromaticValenceAst::matches` function previously treated `Undetermined` as exact-match-only instead of a wildcard, inconsistent with the sibling `ElementAst`/`IsotopeAst`/`HydrogenAst`. Fixing that lets element-only patterns correctly match aromatic-ring atoms they always should have matched.

### Optimization path

**Preprocessing around petgraph (hours, estimated ~30-40% improvement):**
- Switch matcher from `Directed` to `Undirected` — currently builds 2 edges per bond, doubling VF2's search space
- Add `has_match` using `.next()` on the iterator instead of collecting all assignments
- Pre-filter: reject target atoms that can't match any query atom (wrong element, insufficient degree)

**Replace the algorithm (1-2 weeks, parity or better):**
- petgraph implements vanilla VF2 (Cordella 2004). VF2+ (Jüttner & Madarasi 2018) uses better candidate ordering, up to 10x improvement on some graph classes. ~500-800 lines.
- VF2+ can run on petgraph's `Graph` type via the trait API — replace the algorithm, not the data structure.
- Chemistry-aware pruning (ring membership, aromatic matching, degree bounds) is what RDKit layers on top. Requires correct molecular semantics first.

The gap is in algorithm quality, not graph representation.

## Reaction networks: a different workload

petgraph is sufficient for molecular graphs but will not scale to reaction networks (100k–1B+ nodes). At that scale, the concerns change entirely:

- **Cache locality dominates**: petgraph's linked-list edge traversal causes cache misses on every neighbor access. CSR gives sequential reads.
- **Memory overhead**: petgraph stores prev/next pointers per edge. At 1B nodes with avg degree ~4, that's billions of wasted pointer words — potentially the difference between fitting in RAM or not.
- **Parallelism**: Ligra/GBBS-style edge-map primitives are designed for parallel BFS/SSSP over billion-node graphs. petgraph is single-threaded.
- **Compression**: at 1B nodes, WebGraph-style gap encoding may be needed just to fit the adjacency in memory.

Molecular graphs and reaction networks are fundamentally different workloads. If umol needs reaction network path-finding, it will require a separate graph substrate (CSR at minimum, likely with parallel traversal and compression).

## Indirection analysis: four algorithms over MoleculeAst

### What the algorithms need

All four algorithms operate on the same primitive: "give me the neighbors of atom i". The differences are in metadata carried alongside traversal.

| Algorithm | Core topology query | Per-edge metadata | Per-node metadata |
|---|---|---|---|
| Ring enumeration | neighbors(i) → \[j\] | bond index (to build Ring) | aromatic hint (for filtering) |
| VF2 subgraph iso | neighbors(i) → \[j\] | bond attributes (match predicate) | atom attributes (match predicate) |
| Morgan/WL | neighbors(i) → \[(j, bond_idx)\] | bond order | atom invariants |
| DPO reactions | neighbors(i) → \[j\] | bond identity (for deletion) | atom identity (for mapping) |

The topology query is identical. Metadata is always available from parallel arrays indexed by `AtomIdx`/`BondIdx`.

### Current data paths

**Ring enumeration** (from `MoleculeBuilder`):

1. `StableGraph<AtomPattern, BondPattern>` (petgraph, linked-list edges)
2. → `adjacency_list()` → `HashMap<AtomIndex, Vec<AtomIndex>>`
3. → `AtomAdjacency::from_map` → `BTreeMap<AtomIndex, Vec<AtomIndex>>`
4. → optionally `induced()` → new `BTreeMap` (filtered)
5. → `to_dense()` → `DenseProjection { atoms: Vec<AtomIndex>, adj: Vec<Vec<usize>> }`
6. → `enumerate_simple_cycles` / `biconnected_components`

4–5 conversions before the algorithm runs. petgraph stores adjacency internally, but its linked-list representation isn't compatible with what cycle enumeration needs (dense contiguous indices), so we extract → hash → sort → reindex every time.

From `MoleculeAst` (not yet implemented), there would be an additional flat `Vec<BondTuple>` → adjacency step at the start.

**VF2** (from `MoleculeAst`):

1. `Vec<BondTuple>` (flat)
2. → `build_graph()` → `petgraph::Graph<usize, usize, Directed>`
3. → `subgraph_isomorphisms_iter` (petgraph's VF2)

1 conversion. petgraph is the destination — built because `subgraph_isomorphisms_iter` requires it.

**Morgan direct** (from `MoleculeAst`):

1. `Vec<BondTuple>` (flat)
2. → build `Vec<Vec<(usize, usize)>>` adjacency (one pass over bonds)
3. → `ecfp_loop` iterates over slices

1 conversion to dense adjacency, then direct iteration.

**MorganTarget** (from `MoleculeAst`):

1. `Vec<BondTuple>` (flat)
2. → build CSR: `adj: Vec<(usize, usize)>` + `offsets: Vec<usize>`
3. → `ecfp_loop` iterates over slices

1 conversion. `MorganTarget` already implements CSR. If the AST had CSR natively, this constructor would be a borrow.

### Where the indirection comes from

The ring enumeration path is the worst case and reveals the structural problem. petgraph stores edges as linked lists indexed by `NodeIndex(u32)` — each node has a "first outgoing edge" pointer, edges form a linked list per node. This gives O(degree) neighbor iteration but with pointer chasing and no contiguous memory layout.

Cycle enumeration (Johnson's) needs dense contiguous indices — atom 0..n-1, each with a contiguous neighbor array. `DenseProjection` exists to renumber atoms and pack adjacency into `Vec<Vec<usize>>`. petgraph can't provide this directly, so we go through `HashMap` → `BTreeMap` → `Vec<Vec<usize>>`.

The `BTreeMap` intermediate (`AtomAdjacency`) exists because the ring enumerator needs `induced()` — filtering to a subset of atoms. This is a legitimate operation (aromatic subgraph extraction), but the data structure choice is incidental. CSR can produce an induced subgraph just as easily, and the output is already dense.

### Per-algorithm evaluation with CSR in AST

**Ring enumeration:**

CSR → (filter to aromatic atoms if needed → build induced CSR) → `enumerate_simple_cycles` / `biconnected_components`

`DenseProjection` becomes unnecessary — CSR is already dense. `AtomAdjacency` intermediate goes away. The only remaining step is optional aromatic filtering, which is inherent to the problem. Eliminates 3 intermediate data structures (HashMap, BTreeMap, DenseProjection).

**VF2:**

Two options:

- (a) Build `petgraph::Graph` from CSR → use `subgraph_isomorphisms_iter`. Same cost as current (iterate neighbors, add edges). Marginal improvement — iterating CSR neighbors is a slice scan vs iterating flat bonds.
- (b) Implement VF2 directly on CSR. Neighbor queries become slice lookups instead of petgraph iterator traversal. Since petgraph's VF2 is vanilla VF2 (~2.7x slower than RDKit's VF2+), a custom implementation is the long-term path regardless.

Option (a) preserves compatibility. Option (b) is separate work.

**Morgan/WL:**

`morgan_direct` currently builds `Vec<Vec<(usize, usize)>>` from bonds. With CSR in the AST, this is a borrow — no construction. `MorganTarget` already IS a CSR; with CSR in the AST, `MorganTarget::new` copies atom invariants but borrows topology.

WL iteration (same access pattern) would iterate `csr.neighbors(i)` directly.

**DPO reactions:**

DPO = find pattern (VF2 on input graph) → delete matched edges/nodes → add new edges/nodes → output graph.

Input graph is read-only — CSR is ideal for pattern-finding. Output is a new graph constructed from the modified bond list, then CSR-built. Graph surgery works on a mutable intermediate (bond list), not on the source CSR. Natural flow: read from CSR, compute changes, write new CSR.

### Indirection summary

| Algorithm | Current from builder | Current from AST | With CSR in AST |
|---|---|---|---|
| Ring enum | 4 conversions | 5 (extra flat→adj) | 0–1 (induced subgraph if filtering) |
| VF2 | N/A (builder graph) | 1 (build petgraph) | 1 (build petgraph) or 0 (native VF2) |
| Morgan | N/A | 1 (build dense adj) | 0 (borrow slices) |
| DPO | N/A | N/A | 0 for reads, new CSR for output |

Every algorithm currently builds its own adjacency representation from scratch. With CSR in the AST, most algorithms borrow directly. The only remaining conversion is VF2 if petgraph's implementation is kept, and optional aromatic-subgraph induction for ring enumeration.

### Structural change to MoleculeAst

`bonds: Vec<BondTuple>` splits into topology (CSR) and bond attributes (`IndexVec<BondIdx, BondAst>`). Source/target are encoded in the CSR; bond attributes live alongside indexed by `BondIdx`.

Flat `Vec<BondTuple>` remains as the parser's intermediate output and the serialization format. CSR is built from it at parse time. Serialization reconstructs the bond list from CSR + attributes.

### Mutation model

`MoleculeAst` is mutable during solving — the solver narrows `Undetermined` fields to `Lit` values. But this is attribute mutation only, never structural mutation. No atoms or bonds are added or removed during solving. CSR is fully compatible with attribute mutation via parallel arrays. Topology is fixed at parse time and never changes.

Construction happens once (parser → build CSR from bond list). The bond list is the parser's natural output; CSR construction is a single O(|bonds|) pass with a sort.

## Current status

CSR is the right topology substrate for `MoleculeAst`. Every algorithm either borrows CSR directly or builds a transient petgraph `Graph` for VF2 compatibility. petgraph remains available as a library dependency for its VF2 implementation until a native VF2+ replaces it. Reaction networks are a separate concern requiring a different solution.
