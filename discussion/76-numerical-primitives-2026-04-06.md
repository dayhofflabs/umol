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

Algorithms operating on the molecular graph. Live in `umol-graph-core`.

### Cycle enumeration

| Algorithm | Output | Complexity | Use case |
|---|---|---|---|
| Vismara relevant cycles | All cycles not decomposable as XOR of shorter cycles | O(V·E·C) where C = relevant cycle count | SMARTS `R` (ring count), ring-system classification, aromaticity |
| BFS shortest cycle through edge | Smallest ring containing a given edge | O(V+E) per query | SMARTS `r` (smallest ring size), ring membership bit |

Decision: implement both in `umol-graph-core`. Removes the naive DFS `enumerate_simple_cycles`
currently there. SSSR/MCB not implemented — relevant cycles subsume them without the
non-uniqueness problem (see discussion/57).

- Essential cycles (subset of relevant that appear in every MCB) — not needed; relevant cycles
  are cheap enough and essential is too sparse for symmetric-equivalent rings.
- All simple cycles — exponential, not needed for any current use case.

### Matching

| Algorithm | Output | Complexity | Use case |
|---|---|---|---|
| Edmonds' blossom | Maximum matching in general graphs | O(V³) | Kekulization (single), radical detection (unmatched = radical) |
| Uno's enumeration | All perfect/maximum matchings | O(V·E) per matching, Edmonds as subroutine | Kekulization (all representations), tautomer enumeration |

Decision: implement both in `umol-graph-core`.

- Kekulization: run Edmonds on the π-subgraph. Perfect matching → Kekulizable;
  unmatched vertices → radical centers.
- Tautomer enumeration: Uno enumeration on the mobile-H subgraph where vertices are
  mobile-H sites and edges are allowed H-shifts. Each matching assigns double-bond
  positions; unmatched vertices receive the mobile proton.

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

