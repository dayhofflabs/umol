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
- **LEDA / GraphTool (C++)**: separate topology from attributes. GraphTool uses contiguous arrays internally.
- **Ligra / GBBS (MIT, Shun & Blelloch)**: parallel graph algorithms over compressed sparse representations. Edge-map/vertex-map primitives.
- **Graphs.jl (Julia)**: clean interface/storage separation. `SimpleGraph` is `Vector{Vector{Int}}` with no abstraction tax.
- **ECS pattern (game engines)**: topology in one array, atom data in parallel array, bond data in parallel array. Same index space. This is effectively what `MorganTargetOpt` does.

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

## Current status

petgraph is confirmed sufficient for molecular workloads and remains the graph substrate for molecule-level operations (VF2 matcher, Morgan fingerprints). No CSR escape hatch needed for molecules. Reaction networks are a separate concern requiring a different solution.
