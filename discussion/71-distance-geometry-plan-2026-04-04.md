# Distance geometry implementation plan

## Scope

Core DG algorithms only: bounds matrix, triangle smoothing, distance sampling, coordinate embedding. No force-field refinement (ETKDG torsion potentials, chirality enforcement, etc.).

## References

- Crippen & Havel, "Distance Geometry and Molecular Conformation" (1988), pp. 252-253 (triangle smoothing), pp. 312-313 (eigendecomposition embedding)
- Riniker & Landrum, J. Chem. Inf. Model. 2015, 55, 2562-2574 (ETKDG)
- Wang, Witek, Landrum, Riniker, J. Chem. Inf. Model. 2020, 60, 2044-2058 (ETKDGv3, small rings, macrocycles)
- RDKit source: `Code/DistGeom/` (BoundsMatrix, TriangleSmooth, DistGeomUtils)

## Module location

`umol-geometric::algorithms::distance_geometry`

All components are pure numerical algorithms with no chemistry dependencies. Bounds matrix *construction* from molecular topology (bond lengths -> 1-2 bounds, angles -> 1-3 bounds, torsions -> 1-4 bounds) is a separate concern that bridges graph IR and geometric models.

## Components

### 1. BoundsMatrix

Symmetric N x N matrix storing lower bounds in the lower triangle and upper bounds in the upper triangle.

```rust
pub struct BoundsMatrix {
    n: usize,
    data: Vec<f64>,  // row-major N x N
}
```

Methods:
- `new(n) -> Self` (initialized to 0)
- `upper_bound(i, j) -> f64`
- `lower_bound(i, j) -> f64`
- `set_upper_bound(i, j, val)`
- `set_lower_bound(i, j, val)`
- `set_upper_bound_if_tighter(i, j, val)` — only if new value < current upper and > current lower
- `set_lower_bound_if_tighter(i, j, val)` — only if new value > current lower and < current upper

### 2. Triangle smoothing

O(N^3) enforcement of triangle inequality on the bounds matrix. Floyd-Warshall-style triple loop.

```rust
pub fn triangle_smooth(bounds: &mut BoundsMatrix, tol: f64) -> Result<(), SmoothingError>
```

For each triple (i, j, k):
- Upper bound: U(i,j) = min(U(i,j), U(i,k) + U(k,j))
- Lower bound: L(i,j) = max(L(i,j), L(i,k) - U(k,j), L(j,k) - U(i,k))
- Error if L(i,j) > U(i,j) beyond tolerance

### 3. Distance sampling

Sample a concrete distance matrix from smoothed bounds. Each distance drawn uniformly from [L(i,j), U(i,j)].

```rust
pub fn sample_distances(
    bounds: &BoundsMatrix,
    rng: &mut impl Rng,
) -> DMatrix<f64>
```

Returns symmetric N x N distance matrix. Uses nalgebra `DMatrix`.

### 4. Embedding: eigendecomposition (Crippen & Havel)

Classical metric matrix embedding. Given sampled distances D:

1. Compute squared distances: D^2
2. Double-centering: T(i,j) = 0.5 * (d0i^2 + d0j^2 - D^2(i,j)), where d0i^2 = (1/N) * sum_j D^2(i,j) - (1/N^2) * sum_{i,j} D^2(i,j)
3. Eigendecompose T, take top 3 eigenvectors
4. Coordinates: x_i = sqrt(lambda_k) * v_k(i) for k = 1..3

```rust
pub fn embed_eigen(
    distances: &DMatrix<f64>,
    dim: usize,                    // typically 3
    rng: &mut impl Rng,
    config: &EigenEmbedConfig,
) -> Result<DMatrix<f64>, EmbedError>
```

Config:
- `random_negative_eigenvalues: bool` — if true, assign random coords for components with negative eigenvalues; if false, fail
- `max_zero_eigenvalues: usize` — fail if this many or more eigenvalues are zero (default 2)

Returns N x dim coordinate matrix.

### 5. Embedding: random coordinates

Place atoms randomly in a cubic box, to be refined by external minimization.

```rust
pub fn embed_random(
    n: usize,
    dim: usize,
    box_size: f64,
    rng: &mut impl Rng,
) -> DMatrix<f64>
```

Returns N x dim coordinate matrix.

### 6. Error types

```rust
pub enum SmoothingError {
    InfeasibleBounds { i: usize, j: usize, lower: f64, upper: f64 },
}

pub enum EmbedError {
    NegativeEigenvalue { index: usize, value: f64 },
    TooManyZeroEigenvalues { count: usize, max: usize },
    DimensionMismatch,
}
```

## Dependencies

- `nalgebra` (already in umol-geometric): `DMatrix`, eigendecomposition via `SymmetricEigen`
- `rand` (new dependency): `Rng` trait, `Uniform` distribution

## Not in scope

- **Bounds matrix construction from molecular topology**: requires graph IR, bond/angle/torsion knowledge. Lives in `umol-convert` or a dedicated bounds-building module.
- **Force-field refinement**: ETKDG torsion potentials, chirality enforcement, planarity restraints. All ETKDG versions (v1, v2, v3, srETKDGv3) differ only in refinement parameters, not in the core DG algorithm.
- **Conformer ensemble generation**: multiple embeddings + RMSD pruning.

## Test plan

- Triangle smoothing on hand-constructed bounds matrices (feasible and infeasible)
- Eigendecomposition embedding on exact distance matrices (recover known geometries: triangle, square, tetrahedron)
- Round-trip: known coordinates -> distance matrix -> embedding -> compare (up to rotation/translation)
- Random embedding: verify output dimensions and box containment
