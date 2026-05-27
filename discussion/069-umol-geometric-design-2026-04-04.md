# umol-geometric design

## Molecule struct

```rust
pub struct Molecule<G: PointGroup = C1> {
    elements: Vec<Element>,
    coords: Coordinates<G>,
    charge: i32,
    spin: SpinMultiplicity,
}
```

- `Molecule` (no prefix within the module).
- Electron count N_e derived from elements + charge.
- `SpinMultiplicity` from umol-data.
- Default type parameter `G = C1` — writing `Molecule` gives an asymmetric molecule with no overhead.

## Coordinates

```rust
enum Coordinates<G: PointGroup> {
    Cartesian(DMatrix<f64>),
    Symmetric {
        group: G,
        unique_atoms: Vec<usize>,
        full_coords: DMatrix<f64>,
    },
}
```

- Only `Cartesian` implemented. `Symmetric` is a provision for future work.
- Internal coordinates, redundant internal coordinates, and distance matrix representations are future variants.
- nalgebra types throughout.

## Symmetry

Two distinct layers:

1. **Point group symmetry** (geometric): parameterized by `G`. Affects coordinate representation. Detected from geometry or imposed during construction.
   - Detection: full coordinates → identify point group + atom orbits.
   - Construction: asymmetric unit atoms at general/special positions + group operations → full molecule (e.g., C at (x, 0, 0) + D6h → benzene ring).
   - Character table data available in `materials/character_tables/table_data/` (922 `.lis` files, ad hoc ASCII).

2. **Permutational symmetry** (nuclear): partition of atom indices into orbits of identical nuclei. Derived from `elements` alone — not stored, computed on demand. Independent of geometry.

Symmetry-adapted coordinates (linear combinations transforming as irreducible representations) are future work, relevant for geometry optimizers.

## Bond perception

### Formulation

Given N nuclei with elements and 3D positions, assign bond orders b_i ∈ {0, 1, 2, 3} to all atom pairs.

- **Scoring**: a function f(d, Z_a, Z_b) → [p_0, p_1, ..., p_k] maps interatomic distance and element pair to likelihoods for each bond order. No hard connectivity cutoff — bond existence is implicit (p_0 ≈ 1 for distant pairs).
- **Constraints**: Σ b_i = v_a for each atom a (valence sum).
- **Objective**: maximize Π f_i(b_i), equivalently maximize Σ log f_i(b_i).
- **Computational cutoff**: prune pairs with d > 2× max covalent radius sum — optimization, not a model parameter.

### Scoring function (BondDistanceModel)

Encapsulates the distance → bond order likelihood mapping.

Signature: f(d, (Z_a, q_a), (Z_b, q_b)) → [p_0, p_1, ..., p_k], where q is an optional charge (relevant for transition metals with oxidation-state-dependent radii; absent/zero for most organic elements). The (Element, Option<i8>) charge annotation lives on BondDistanceModel, not on the molecule struct. **Charge-dependent radii are not yet implemented.**

**Model — Pyykko & Atsumi radii + Gaussian/sigmoid:**

Parameters per element: r_cov^(1), r_cov^(2), r_cov^(3) (single, double, triple bond covalent radii). Source: Pyykko & Atsumi (2009), verified against primary paper (J. Phys. Chem. A 2015, 119, 2326). Stored in `umol-params/src/covalent_radii.rs`.

Expected bond length for order k: μ_k(a, b) = r_cov^(k)(Z_a) + r_cov^(k)(Z_b).

Likelihoods:
- k ≥ 1: p_k ∝ exp(−(d − μ_k)² / 2σ_k²)
- k = 0: p_0 = 1 / (1 + exp(−α(d − μ_1 − δ)))

Default parameters: σ = 0.10 Å (all orders), α = 15.0, δ = 0.30 Å.

### Configuration (BondPerceptionConfig)

```rust
pub struct BondPerceptionConfig {
    pub model: BondDistanceModel,
    pub target_valences: Option<Vec<u8>>,  // None → octet-rule default
    pub max_iter: usize,                   // default 200
    pub step_scale: f64,                   // default 0.5
}
```

Target valences default to the octet rule: ve ≤ 4 → ve, else 8 − ve.

### Solver

**Implemented: Lagrangian relaxation** with subgradient method on the dual of the valence constraints. Generic solver in `algorithms/optimization.rs`, independent of chemistry. Each subproblem picks the best of ≤4 values per bond. When the Lagrangian solution is infeasible (degenerate cases like benzene where symmetric dual costs prevent convergence), a greedy primal recovery step assigns bond orders by marginal log-likelihood gain subject to remaining valence capacity.

**Future: Belief propagation** — max-product BP on the molecular graph for marginal probability distributions P(b_i = k) per bond.

### Feasibility

Solution is not guaranteed to exist (infeasible valence constraints). Not unique in general (resonance structures). The solver handles both: reports `feasible` flag and per-atom `valence_residuals`.

## Parameters (umol-params)

- **Covalent radii**: Pyykko & Atsumi values. Three radii per element (all 118): r_cov^(1), r_cov^(2), r_cov^(3) for single, double, triple bonds. Model parameters, not fundamental element data.
- **Sigmoid and width parameters**: currently hardcoded defaults in `BondDistanceModel::default()`. Will move to umol-params when tuned.
- Future: trained/fitted parameters from structural databases, charge-dependent radii for transition metals.

## Conversions

- `Molecule` (geometric) → `umol_graph::Molecule`: bond perception gives `BondPerceptionResult`; **conversion to graph `Molecule` not yet implemented**.
- `umol_graph::Molecule` → `Molecule` (geometric): distance geometry embedding (deferred).
- `umol_graph::Molecule` + external coordinates → `Molecule` (geometric): simple assembly (deferred).

## Implementation status

| # | Item | Status |
|---|------|--------|
| 1 | `Molecule<G>` with `Coordinates::Cartesian` | Done |
| 2 | Pyykko & Atsumi covalent radii in umol-params | Done |
| 3 | `BondDistanceModel` scoring function | Done |
| 4 | Lagrangian relaxation solver (`algorithms/optimization`) | Done |
| 5 | `BondPerceptionResult` → `umol_graph::Molecule` conversion | Next |
| 6 | Tests on simple organic molecules | Done (5 molecules) |
| 7 | Belief propagation solver | Future |
| 8 | Charge-dependent radii | Future |
| 9 | Symmetry detection/construction | Future |
