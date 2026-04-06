# Symmetry module roadmap

## Current state (2026-04-06)

Implemented:
- Point group detection via libmsym (finite groups + C∞v/D∞h)
- Character tables with class representatives and irrep characters
- Character table display (`CharacterTableDisplay`)
- Irrep algebra: direct product, reduce, symmetric square, antisymmetric square
- Selection rules: electric dipole, magnetic dipole, Raman, contains_totally_symmetric
- Translation/rotation/quadratic irreps
- SALCs for arbitrary basis functions
- Symmetry coordinates (3N DOF → trans/rot/vib by irrep)
- Equivalence sets and atom permutations
- SymmetryOp: Display (`E`, `C3²`, `σv`, `S4³`, `i`), `transform_point`, `is_proper`
- PointGroup queries: `is_chiral`, `has_inversion`, `is_abelian`, `is_cyclic`, `is_cubic`, `has_complex_irreps`, `principal_axis_order`, `totally_symmetric_irrep`
- `Irrep::is_gerade()` for both finite and linear centrosymmetric groups
- `SymmetryOp::class` as `usize` with documented indexing semantics
- `generate_symmetry_images` (build molecule from asymmetric unit + group)
- Real representation treatment: complex conjugate irreps fused into real 2D reps (standard QC convention)

## Near-term: display tables

### Operation multiplication table
Display the Cayley table for a point group: rows and columns are operations, cells show the product. Requires operation composition (matrix multiply + identify resulting operation). Useful for pedagogy and verification.

### Irrep direct product table
Display all pairwise direct products of irreps. Compact table format: row × column → decomposition. All data is already available via `direct_product()`.

## Near-term: continuous symmetry measures and symmetrize_to

### Motivation
`symmetrize()` detects the group and snaps to it. `symmetrize_to()` should snap to a *user-specified* group. This requires measuring how far a structure is from a given symmetry — Avnir's Continuous Symmetry Measures (CSM).

### Two functions

- `symmetry_measure(centers, group) -> f64` — CSM distance to nearest G-symmetric structure. No modification.
- `symmetrize_to(centers, group, threshold) -> Result<SymmetryResult>` — project onto nearest G-symmetric structure. Fail if CSM > threshold.

### Algorithm

1. Assign atom permutations under each operation of G (approximate nearest-neighbor matching, not exact).
2. For each atom i, compute symmetry-averaged position: p̃_i = (1/h) Σ_R R⁻¹ p_{R(i)}.
3. CSM = Σ |p_i - p̃_i|² / Σ |p_i - center|² (normalized).
4. The p̃_i are the nearest G-symmetric structure.

### Relation to existing code

- `symmetrize()` uses libmsym (exact mapping from detected group) — unchanged.
- `symmetrize_to()` uses CSM projection (approximate mapping to specified group) — new algorithm.
- `generate_symmetry_images()` builds a molecule from an asymmetric unit — unchanged.

## Medium-term: subgroups and correlation

Source: Katzer's data (subgroup chains and correlation tables for all common point groups).

- Subgroup enumeration for a given group
- Correlation tables: irrep mapping between group and subgroup
- Symmetry descent chains

Applications: Jahn-Teller analysis, crystal field splitting, working with molecules whose actual symmetry is lower than idealized.

## Medium-term: tensor symmetry and selection rules

### Motivation
Selection rules for spectroscopic properties beyond electric dipole require decomposing higher-rank tensors under point group operations. Current infrastructure handles rank-1 (dipole) and rank-2 symmetric (Raman/quadrupole). Need general machinery.

### Target tensor properties

| Property | Tensor | Rank | Index symmetry |
|---|---|---|---|
| Electric dipole (μ) | polar vector | 1 | — |
| Magnetic dipole (m) | axial vector | 1 | — |
| Electric quadrupole (Θ) | symmetric | 2 | [α²] |
| Polarizability (α) | symmetric | 2 | [α²] |
| Optical rotation tensor (G') | general | 2 | α ⊗ α (no index symmetry) |
| Hyperpolarizability (β) | rank-3, partial | 3 | symmetric in last two indices |
| ROA tensors | mixed | 2-3 | products of translation, rotation, quadrupole |
| Second hyperpolarizability (γ) | rank-4, partial | 4 | intrinsic permutation symmetry |

### Required infrastructure
- Tensor product decomposition of arbitrary irrep combinations
- Symmetric and antisymmetric parts of higher powers (rank > 2)
- Mixed-symmetry tensor decomposition (e.g., β symmetric in last two indices only)
- Selection rules: does a given tensor irrep decomposition contain the product of initial and final state irreps?

### Building blocks already in place
- `direct_product(a, b)` — rank-2 products
- `symmetric_square(a)` / `antisymmetric_square(a)` — rank-2 symmetric/antisymmetric
- `reduce(characters)` — decomposition into irreps
- `translation_irreps()` / `rotation_irreps()` / `quadratic_irreps()` — rank-1 and rank-2 bases

## Long-term: double groups

For systems with spin-orbit coupling. The point group is extended by a 2π rotation (which maps spinors to −spinors), doubling the number of operations. Same algebraic framework but with half-integer angular momentum representations.

Applications: heavy-element chemistry, spin-forbidden transitions, relativistic electronic structure.

## Not planned
- Space groups (not relevant for molecular symmetry)
