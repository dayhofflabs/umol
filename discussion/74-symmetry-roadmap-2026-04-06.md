# Symmetry module roadmap

## Current state (2026-04-06)

Implemented:
- Point group detection via libmsym (finite groups + C∞v/D∞h)
- Character tables with class representatives and irrep characters
- Irrep algebra: direct product, reduce, symmetric square, antisymmetric square
- Selection rules: electric dipole, magnetic dipole, Raman, contains_totally_symmetric
- Translation/rotation/quadratic irreps
- SALCs for arbitrary basis functions
- Symmetry coordinates (3N DOF → trans/rot/vib by irrep)
- Equivalence sets and atom permutations

## Near-term: structural improvements

### SymmetryOp labels
SymmetryOp carries kind/order/power but has no Display. Need string labels: `E`, `C3²`, `σv`, `σh`, `S4³`, `i`. Required for character table formatting.

### Character table formatter
Pretty-print a character table given the group. Depends on SymmetryOp labels.

### Chirality query
`PointGroup::is_chiral() -> bool` — true iff the group contains no improper operations. Trivial.

### SymmetryOp::orientation cleanup
`Horizontal`/`Vertical`/`Dihedral`/`None` is a libmsym-ism attached to all operations but only meaningful for reflections. Consider restricting to reflection operations or documenting the semantics.

### SymmetryOp::class typing
Currently `i32`. Should be a newtype or at minimum documented as indexing into `class_sizes`/`class_reps`/`Irrep::characters`.

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
