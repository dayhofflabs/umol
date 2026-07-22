# Numerical and geometric primitives

Catalog of general-purpose numerical primitives that umol uses or may need.
Not a roadmap — an inventory of techniques, when each is appropriate, and the
state of available implementations. Items here are pulled into specific module
roadmaps when a concrete need arises.

## Linear assignment

Match elements of two sets to minimize total cost. Standard formulation:
given an N×N cost matrix C, find permutation π minimizing Σᵢ C(i, π(i)).

### Methods

| Method | Complexity | Notes |
|---|---|---|
| Greedy | O(N²) | Picks smallest entries iteratively. Not optimal, may pick poorly. |
| Hungarian (Munkres / Kuhn-Munkres) | O(N³) | Optimal. Standard choice. |
| Jonker-Volgenant | O(N³) but faster constants | Modern variant of Hungarian, often 2-10× faster in practice. |
| Auction algorithm | O(N² log N) average | Better for sparse cost matrices. |

### Use cases in umol

- **CSM permutation search** — atom-to-image matching under symmetry operations
- **RMSD with optimal atom assignment** — comparing structures with no a priori mapping
- **Atom mapping in reactions** — reactant-to-product correspondence
- **Trajectory frame matching** — MD analysis, conformer alignment
- **Equivalence set membership** — when atoms are nearly but not exactly equivalent

### Available crates

- `pathfinding` — actively maintained, includes Hungarian (Kuhn-Munkres)
- `lap` / `lap_jv` — Jonker-Volgenant variant
- `linfa-clustering` has it as a dependency for some clustering methods

Decision: use external crate, do not reimplement.

## Sphere sampling (S²)

Generate points (or directions) on the 2-sphere. Different methods optimize
different criteria: uniformity, integration accuracy, hierarchical refinement,
or specific symmetry.

### Methods

| Method | Property | Construction | Use case |
|---|---|---|---|
| Fibonacci sphere (golden angle) | Quasi-uniform, deterministic | `θ = 2πk/φ`, `z = 1 − 2k/(N−1)` | General direction sampling, low N |
| HEALPix | Equal-area pixelization, hierarchical | Recursive subdivision | Cosmology, hierarchical refinement |
| Lebedev quadrature | Optimal polynomial integration up to degree L | Octahedral symmetry constraints | Orientation-averaged observables |
| Spherical t-designs | Exact polynomial integration up to degree t | Variational construction | Higher-order moment integration |
| Repulsion (Thomson) | Local minimum of Coulomb energy | Iterative minimization | Aesthetic uniformity, expensive |
| Random uniform | IID samples | Inverse CDF | Monte Carlo |

### Use cases in umol

- **CSM direction initialization** — Fibonacci, ~10–100 points
- **Orientation-averaged tensor properties** (Raman, ROA, hyperpolarizabilities, ⟨α⟩, ⟨β²⟩, ⟨γ²⟩) — Lebedev is the QC standard, exact for low-degree polynomials
- **Powder pattern simulation** — needs uniform crystal orientations
- **Visualization sampling** — surface point clouds, isosurface generation
- **Initial guesses for SO(3) optimization** — Fibonacci or Hopf as starting set

### Available crates

- `lebedev_laikov` — Lebedev quadrature points (the standard implementation port)
- `healpix` — HEALPix pixelization
- Fibonacci sphere is trivial (~10 lines), implement inline

## Rotation sampling (SO(3))

Sample rigid rotations. Distinct from S² sampling because each direction has a
circle of rotations about it.

### Methods

| Method | Property | Notes |
|---|---|---|
| Hopf fibration (Yershova et al.) | Incremental, multiresolution, near-uniform | S³ → S² × S¹ factorization. Allows refinement on demand. |
| Shoemake (quaternion) | Random uniform | Standard random method, IID quaternions |
| Successive orthogonal images (SOI) | Deterministic uniform | Static grid |
| 24-cell, 600-cell deterministic | Highly uniform low-N grids | Polytope vertices in S³ |
| Euler angle grid | Trivial | Non-uniform near poles, generally avoid |

### Use cases in umol

- **Conformer / pose alignment** — full orientation matters
- **Docking-like problems** — searching rigid-body orientations
- **Full O(3) integration** for tensor invariants when Lebedev is insufficient
- **MD / Monte Carlo initial conditions** — random rigid orientations

### Available crates

Limited. Likely need to implement Hopf fibration sampling from the Yershova
2010 paper if and when needed. Quaternion uniform random is in `nalgebra` /
`rand_distr`.

## Graph algorithms

Algorithms operating on the molecular graph. Live in `umol-graph-core`. All dispatch through
enum arguments (`Graph::method(&self, ..., alg: AlgorithmEnum)`).

### Connected components

| Algorithm | Complexity | Reference |
|---|---|---|
| BFS flood fill | O(V+E) | Standard BFS traversal |

Implemented: `Graph::connected_components(alg: ConnectedComponentsAlgorithm)`.

### Biconnected components

| Algorithm | Complexity | Reference |
|---|---|---|
| Tarjan DFS | O(V+E) | Tarjan 1972 "Depth-first search and linear graph algorithms" |

Implemented: `Graph::biconnected_components(alg: BiconnectedComponentsAlgorithm)`.

### Cycle enumeration

| Algorithm | Output | Complexity | Use case |
|---|---|---|---|
| Vismara relevant cycles | All cycles not decomposable as XOR of shorter cycles | O(V·E·C) where C = relevant cycle count | SMARTS `R` (ring count), ring-system classification, aromaticity |
| BFS shortest cycle through edge | Smallest ring containing a given edge | O(V+E) per query | SMARTS `r` (smallest ring size), ring membership bit |
| BFS shortest cycle through node | Min over incident edges | O(V+E) per query | SMARTS `r` per atom |

Implemented: `Graph::enumerate_cycles(max_cycle_size, alg: CycleEnumerationAlgorithm)`,
`Graph::shortest_cycle_through_edge(edge, alg: ShortestCycleAlgorithm)`,
`Graph::shortest_cycle_through_node(node, alg: ShortestCycleAlgorithm)`.

SSSR/MCB not implemented — relevant cycles subsume them without the non-uniqueness problem
(see discussion/57).

- Essential cycles (subset of relevant that appear in every MCB) — not needed; relevant cycles
  are cheap enough and essential is too sparse for symmetric-equivalent rings.
- All simple cycles — exponential, not needed for any current use case.

### Matching

| Algorithm | Output | Complexity | Use case |
|---|---|---|---|
| Edmonds' blossom (Gabow simplification) | Maximum matching | O(V³) | Kekulization (single), radical detection (unmatched = radical) |
| Branch-and-bound with Edmonds oracle | All perfect/maximum matchings | Exponential worst-case | Kekulization (all representations), tautomer enumeration |

Implemented: `Graph::maximum_matching(node_order, alg: MaximumMatchingAlgorithm) -> Result<Matching, MaximumMatchingError>`,
`Graph::enumerate_perfect_matchings(alg: MatchingEnumerationAlgorithm)`,
`Graph::enumerate_maximum_matchings(alg: MatchingEnumerationAlgorithm)`.

References: Edmonds 1965, Gabow 1976 simplification. Ref impl: cp-algorithms.com.

### Maximum independent set

| Algorithm | Complexity | Use case |
|---|---|---|
| Branch-and-bound with greedy upper bound | Exponential worst-case | Clar's rule (maximum number of disjoint aromatic sextets) |

Implemented: `Graph::maximum_independent_set(alg: MaximumIndependentSetAlgorithm)`.

### Subgraph isomorphism

| Algorithm | Complexity | Use case |
|---|---|---|
| VF2 | O(V!·V) worst-case, fast on sparse molecular graphs | SMARTS matching, substructure search |

Implemented: `Graph::subgraph_isomorphisms(query, node_match, edge_match, alg: SubgraphIsomorphismAlgorithm)`,
`Graph::subgraph_isomorphisms_at(query, anchor, node_match, edge_match, alg: SubgraphIsomorphismAlgorithm)`.
Target graph is the receiver; query is the argument. Match predicates are pairwise closures
(lazy evaluation — most candidate pairs never tested).

Reference: Cordella et al. 2004 "A (sub)graph isomorphism algorithm for matching large graphs".

### Automorphism and canonical labeling

| Algorithm | Output | Use case |
|---|---|---|
| nauty (sparse) | Orbit partition, canonical labeling, group order | Symmetry-equivalent atom detection, canonical SMILES |

Implemented: `Graph::automorphisms(node_color, alg: AutomorphismAlgorithm)`.
Coloring function `Fn(NodeId) -> C where C: Ord + Copy` encodes vertex partition (nauty
requires a total ordering). Returns `Automorphism` with orbit queries, canonical labeling,
and `AutoGroupOrder` (exact u32 or approximate f64).

Reference: McKay & Piperno 2014 "Practical graph isomorphism, II". Impl: `nauty-Traces-sys` FFI.

## Other potentially relevant primitives

Tracked here for completeness, not yet needed:

### Quasi-Monte Carlo sequences
- Sobol, Halton, Niederreiter sequences for low-discrepancy sampling in
  arbitrary dimension. Useful for high-dimensional integration (e.g., conformer
  ensemble averaging, free-energy calculations).

### Optimization on manifolds
- Riemannian gradient descent on SO(3), Stiefel manifold for orbital
  optimization, Grassmannian for subspace problems. Crate: `manopt` analogs.

### Linear assignment generalizations
- **Earth mover's distance** / Wasserstein for soft matching when atom
  identities aren't fixed
- **Optimal transport** for comparing distributions of atomic positions

### Spatial data structures
- KD-tree, ball tree, octree for nearest-neighbor queries. Crate: `kiddo`,
  `kdtree`. Already needed for non-bonded interaction lists, equivalence set
  detection at scale.
