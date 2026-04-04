# libmsym integration: API design

## Status: draft (revised after design discussion)

## Design principles

1. **Semantics first, implementation second.** The public API models the domain (point groups, molecular symmetry, symmetry-adapted bases). The C library is an implementation detail hidden behind the API boundary.
2. **Every molecule has symmetry.** C1 is the trivial case, not the absence of symmetry. There is no separate "symmetric molecule" type.
3. **Molecules are immutable.** Transformations (symmetrization, perception) produce new molecules.
4. **PointGroup is a value type**, not a type parameter. A `Vec<Molecule>` with mixed symmetries works without `Box<dyn>`. A registry can be added later if needed; going the other direction is harder.
5. **Heavy computations (SALCs, projectors) are methods that return owned results.** The molecule stores the cheap ingredients (equivalence sets, atom permutations). Expensive derived data is computed on demand, owned by the caller.

## Domain model

Three distinct concepts:

### 1. PointGroup — abstract algebraic object

One C2v exists. It has no knowledge of any molecule.

```rust
struct PointGroup {
    kind: PointGroupKind,   // enum: Ci, Cs, Cn, Cnv, ...
    n: i32,                 // principal axis order
    name: String,           // "C2v", "Td", etc.
    order: usize,           // number of group elements
    operations: Vec<SymmetryOperation>,
    character_table: CharacterTable,
}
```

Contains: symmetry operations (group elements) with 3×3 matrices, character table, irreps, class structure. Can answer purely algebraic questions: direct products, selection rules, subgroup relationships.

Constructible by name (`PointGroup::from_schoenflies("C2v")`) without any molecular data. For C1: one operation (identity), one irrep (A), trivial character table.

### 2. Molecule — always carries symmetry data

```rust
struct Molecule {
    elements: Vec<Element>,
    coords: Coordinates,
    charge: i32,
    multiplicity: SpinMultiplicity,

    // Symmetry data — always present, trivial for C1
    group: PointGroup,
    equivalence_sets: Vec<Vec<usize>>,  // atom index orbits
    atom_permutations: Vec<Vec<usize>>, // one per operation
    orientation: Orientation,            // transform to standard orientation
}
```

A freshly constructed molecule from Cartesian coordinates gets `PointGroup::c1()` with identity permutations and singleton equivalence sets (one per atom). This is cheap: a few small `Vec<usize>`.

Perception produces a new `Molecule` with the discovered group and associated data:
```rust
impl Molecule {
    fn perceive_symmetry(&self, thresholds: &Thresholds) -> Result<Molecule, Error>;
    fn symmetrize(&self, thresholds: &Thresholds) -> Result<Molecule, Error>;
}
```

Both return new `Molecule` values. `perceive_symmetry` discovers the group and populates equivalence sets / permutations. `symmetrize` additionally snaps coordinates to exact symmetry.

### 3. Symmetry-adapted bases — computed on demand

A symmetry basis is a representation of the point group on a vector space derived from the molecule. The molecule stores the raw ingredients; the analysis is a method that returns an owned result.

```rust
impl Molecule {
    /// SALCs of AO basis functions via libmsym.
    fn salcs(&self, basis: &[BasisFunction]) -> Result<Vec<Salc>, Error>;


    /// Symmetry-adapted coordinates: Γ_3N → Γ_trans + Γ_rot + Γ_vib.
    /// Pure Rust, using atom_permutations + character table.
    /// Returns SALCs of atomic displacements classified by irrep and category.
    fn symmetry_coordinates(&self, weighting: MassWeighting) -> SymmetryCoordinates;
}
```

The caller owns the result. If they need it multiple times, they hold onto it. The molecule is not burdened with cached derived data.

## PointGroup type

```rust
#[derive(Debug, Clone, PartialEq)]
struct PointGroup {
    kind: PointGroupKind,
    n: i32,
    name: String,
    order: usize,
    operations: Vec<SymmetryOp>,
    character_table: CharacterTable,
}

impl PointGroup {
    /// Construct the trivial group.
    fn c1() -> Self;

    // ... add labels for molecular point groups with n <= 8 via macro

    /// Construct by type and principal axis order.
    fn new(kind: PointGroupKind, n: i32) -> Result<Self, Error>;

    /// Construct by Schoenflies symbol.
    fn from_schoenflies(name: &str) -> Result<Self, Error>;

}
```

### SymmetryOp type

```rust
struct SymmetryOp {
    kind: SymmetryOpKind,  // Identity, ProperRotation, ImproperRotation, Reflection, Inversion
    order: i32,
    power: i32,
    orientation: SymmetryOpOrientation,  // None, Horizontal, Vertical, Dihedral
    vector: [f64; 3],             // axis or plane normal
    class: i32,                   // index into character table classes
    matrix: Matrix3<f64>,         // the 3×3 transformation matrix
}
```

Note: the 3×3 matrix is stored on the operation. libmsym provides axis + type + order; the matrix is derived from those. Having it pre-computed avoids reconstructing it every time we build representation matrices.

### CharacterTable

```rust
struct CharacterTable {
    irreps: Vec<Irrep>,
    class_sizes: Vec<i32>,
    class_operations: Vec<SymmetryOp>,  // one representative per class
    characters: Vec<Vec<f64>>,  // [irrep_index][class_index]
    order: usize,
}

#[derive(Clone, PartialEq, Eq, Hash)]
struct Irrep {
    name: String,
    dimension: i32,
    index: usize,
}
```

Methods on CharacterTable (purely algebraic, no molecular data):

```rust
impl CharacterTable {
    /// Decompose Γ_a ⊗ Γ_b into irreps.
    fn direct_product(&self, a: &Irrep, b: &Irrep) -> Vec<(Irrep, u32)>;

    /// Does Γ_a ⊗ Γ_b ⊗ Γ_c contain the totally symmetric irrep?
    fn contains_totally_symmetric(&self, a: &Irrep, b: &Irrep, c: &Irrep) -> bool;

    /// Reduce an arbitrary representation (characters per class) into irreps.
    fn reduce(&self, characters: &[f64]) -> Vec<(Irrep, u32)>;

    /// Irreps spanned by (x, y, z) — the vector representation.
    fn translation_irreps(&self) -> Vec<(Irrep, u32)>;

    /// Irreps spanned by (Rx, Ry, Rz) — the pseudovector representation.
    fn rotation_irreps(&self) -> Vec<(Irrep, u32)>;
}
```

`translation_irreps` and `rotation_irreps` are computed from the 3×3 matrix characters of the operations (trace of Cn = 1 + 2cos(2π/n), etc.) and their symmetric/antisymmetric parts. No external data needed, but Katzer tables serve as validation.

## Molecule construction and perception

```rust
impl Molecule {
    /// From Cartesian coordinates. Symmetry defaults to C1.
    fn from_cartesian_angstrom(
        elements: Vec<Element>,
        coords: &[f64],
        charge: i32,
        spin: SpinMultiplicity,
    ) -> Self;

    /// Detect point group symmetry. Returns a new molecule with the
    /// discovered group, equivalence sets, atom permutations, and
    /// standard orientation.
    fn perceive_symmetry(&self, thresholds: &Thresholds) -> Result<Molecule, Error>;

    /// Symmetrize: perceive + snap coordinates to exact symmetry.
    /// Returns a new molecule with exact symmetry.
    fn symmetrize(&self, thresholds: &Thresholds) -> Result<Molecule, Error>;

    /// Generate full molecule from asymmetric unit + group name.
    fn symmetrize_to(
        point_group_kind: &PointGroupKind,
        elements: Vec<Element>,
        coords: &[f64],
    ) -> Result<Molecule, Error>;
}
```

Internally, `perceive_symmetry` and `symmetrize`:
1. Create a transient libmsym context (not exposed)
2. Set atoms, call `msymFindSymmetry` (and `msymSymmetrizeElements` for symmetrize)
3. Extract: group, operations, character table, equivalence sets
4. Compute atom permutations (see open question below)
5. Construct and return the new `Molecule`
6. Drop the C context

## Displacement analysis (pure Rust)
```rust
struct SymmetryCoordinates {
    gamma_total: Vec<(Irrep, u32)>,
    gamma_trans: Vec<(Irrep, u32)>,
    gamma_rot: Vec<(Irrep, u32)>,
    gamma_vib: Vec<(Irrep, u32)>,
    bases: Vec<IrrepBasis>,  // each with Displacement terms and CoordinateCategory
    weighting: MassWeighting,
}

// IrrepBasis (defined in SALC section) contains:
//   irrep, partners: Vec<Salc>, category, weighting

impl Molecule {
    fn symmetry_coordinates(&self, weighting: MassWeighting) -> SymmetryCoordinates;
}
```

Algorithm:

1. **Build 3N×3N representation matrices** for each operation, using `atom_permutations` and the 3×3 operation matrices stored on `SymmetryOp`.
2. **Compute Γ_3N characters**: χ(R) = Σ_{atoms fixed by R} tr(3×3 matrix of R). Only atoms mapped to themselves contribute.
3. **Reduce Γ_3N** via `CharacterTable::reduce`.
4. **Identify Γ_trans and Γ_rot** via `CharacterTable::translation_irreps()` and `rotation_irreps()`.
5. **Subtract** to get Γ_vib.
6. **Project** via P_μ = (l_μ/h) Σ_R χ_μ(R)* · D_3N(R) applied to trial vectors. Orthogonalize within each irrep subspace (SVD). Classify each resulting coordinate as translation, rotation, or vibration.

Steps 1-5 use only characters (cheap). Step 6 needs the full 3N×3N matrices (allocated transiently). nalgebra handles all linear algebra.

Mass-weighting: the projection in step 6 can use mass-weighted coordinates (for normal mode analysis) or unweighted (for symmetry classification). Both are useful; the method should accept an option.

## Selection rules

Methods on `CharacterTable`, not on `Molecule`, because these are purely algebraic:

```rust
impl CharacterTable {
    /// Electric dipole allowed? Checks Γ_i ⊗ Γ(x,y,z) ⊗ Γ_f ⊃ A1.
    fn electric_dipole_allowed(&self, initial: &Irrep, final_: &Irrep) -> bool;

    // Electric quadrupole allowed? (Raman + s)
    fn electric_quadrupole_allowed(&self, initial: &Irrep, final_: &Irrep) -> bool;

    /// Raman allowed? Checks against irreps of (x², y², z², xy, xz, yz).
    fn raman_allowed(&self, initial: &Irrep, final_: &Irrep) -> bool;

    /// Magnetic dipole allowed? Checks against irreps of (Rx, Ry, Rz).
    fn magnetic_dipole_allowed(&self, initial: &Irrep, final_: &Irrep) -> bool;
}
```

Convenience methods on `SymmetryCoordinates`:

```rust
impl SymmetryCoordinates {
    /// Which vibrational irreps are IR-active?
    fn ir_active(&self, ct: &CharacterTable) -> Vec<&Irrep>;

    /// Which vibrational irreps are Raman-active?
    fn raman_active(&self, ct: &CharacterTable) -> Vec<&Irrep>;
}
```

## Unified SALC type

All symmetry-adapted linear combinations — of atoms, displacements, or orbitals — share the same structure:

```rust
/// What kind of object is being combined.
enum SalcBasis {
    /// s-type: atoms themselves (permutation representation).
    Atom,
    /// p-type: Cartesian displacement component on an atom.
    Displacement(CartesianAxis),
    /// Spherical harmonic centered on an atom.
    /// l and m determine transformation properties; n is irrelevant for symmetry.
    SphericalHarmonic { l: i32, m: i32 },
}

enum CartesianAxis { X, Y, Z }

/// One weighted contribution to a SALC.
struct SalcTerm {
    atom: usize,
    basis: SalcBasis,
    coefficient: f64,
}

/// A symmetry-adapted linear combination.
struct Salc {
    terms: Vec<SalcTerm>,
}

/// A complete irrep basis: one or more partner SALCs spanning an irrep.
struct IrrepBasis {
    irrep: Irrep,
    partners: Vec<Salc>,                   // length = irrep.dimension
    category: Option<CoordinateCategory>,  // Translation/Rotation/Vibration (p-type only)
    weighting: MassWeighting,              // how the SALC was constructed
}

enum CoordinateCategory {
    Translation,
    Rotation,
    Vibration,
}

enum MassWeighting {
    Unweighted,
    MassWeighted,
}
```

This type is used uniformly:
- `Molecule::symmetry_coordinates()` returns `Vec<IrrepBasis>` with `Displacement` terms
- `Molecule::salcs(basis)` returns `Vec<IrrepBasis>` with `SphericalHarmonic` terms (via libmsym)
- s-type reduction returns `Vec<IrrepBasis>` with `Atom` terms

The `Salc` type is naturally sparse (only atoms in the relevant equivalence sets have non-zero coefficients). Internally during projection, DMatrix representations are used for linear algebra; the `Salc` term list is the public output format.

### AO SALCs via libmsym

```rust
impl Molecule {
    /// Compute SALCs for AO basis functions.
    /// The input specifies which (l, m) shells to include on which atoms.
    /// Uses libmsym internally (msymGetSALCs).
    fn ao_salcs(&self, basis: &[OrbitalIndex]) -> Result<Vec<IrrepBasis>, Error>;
}

struct OrbitalIndex {
    atom_index: usize,
    n: i32,  // needed by libmsym for bookkeeping, not for symmetry
    l: i32,
    m: i32,
}
```

## Crate structure and dependencies

```
umol-geometric
  ├── umol-data       (Element, atomic masses, units)
  └── umol-msym       (PointGroup, CharacterTable, perception/symmetrization impl)
        └── umol-msym-sys
              └── libmsym (C, git submodule)
```

- `umol-msym` depends on `umol-data` (needs `Element` for mass lookup in conversions to libmsym atoms).
- `umol-geometric` depends on `umol-msym` (needs `PointGroup`, `CharacterTable`, etc. as fields and return types).
- The `Molecule` methods (`perceive_symmetry`, `displacement_analysis`, etc.) live in `umol-geometric` since `Molecule` is defined there.
- `umol-msym` provides: `PointGroup`, `CharacterTable`, `Irrep`, `SymmetryOperation`, `Thresholds`, `Error`, and the internal `Context` wrapper for the C API.
- `umol-msym` does NOT depend on `umol-geometric`. No circular dependency.

## What libmsym provides vs. what we build

| Need | Source | Notes |
|---|---|---|
| Perception (atoms → group) | libmsym via `umol-msym` | `msymFindSymmetry` |
| Symmetrization | libmsym via `umol-msym` | `msymSymmetrizeElements` |
| Character table | libmsym via `umol-msym` | `msymGetCharacterTable` |
| Symmetry operations | libmsym via `umol-msym` | `msymGetSymmetryOperations` |
| Equivalence sets | libmsym via `umol-msym` | `msymGetEquivalenceSets` |
| Generation from asymmetric unit | libmsym via `umol-msym` | `msymGenerateElements` |
| AO SALCs | libmsym via `umol-msym` | `msymGetSALCs` |
| 3×3 operation matrices | Rust, derived from operation type/axis/order | Stored on `SymmetryOperation` |
| Atom permutations | Rust, from operation matrices + atom positions | Resolved: compute in Rust |
| Direct products, reduction | Rust, from character table | `CharacterTable` methods |
| Selection rules | Rust, from character table | `CharacterTable` methods |
| Translation/rotation irreps | Rust, from 3×3 matrix characters | `CharacterTable` methods |
| Displacement analysis + SACs | Rust, from permutations + 3×3 matrices + character table | `Molecule::displacement_analysis` |

## Current implementation state

- `umol-msym-sys`: complete. All FFI bindings, libmsym builds via `cc`.
- `umol-msym`: scaffolded. `Context` wrapper with 7 unit tests + 18 Katzer comparison tests (25 total, all passing). Covers: perception, character tables, direct products, selection rules, validated against Katzer data for 18 point groups.
- `umol-geometric::Molecule`: currently has `<G: PointGroup>` type parameter — to be replaced with `PointGroup` value field.

## Open questions

1. ~~**Atom permutations.**~~ **Resolved.** Computed in Rust. For each operation, apply the 3×3 matrix to each atom position, find the nearest atom within threshold. libmsym computes this internally (`findPermutation` in `permutation.c`) but stores it per-equivalence-set and does not expose it through the public API. Our Rust implementation produces global permutations (atom i → atom j across the whole molecule), which is the format we need for building 3N×3N representation matrices. The algorithm is O(N²) per operation — negligible for molecular sizes.

2. ~~**Katzer data as test oracle.**~~ **Resolved.** Automated comparison implemented in `umol-msym/tests/katzer_comparison.rs`. 18 tests cover C2, Cs, Ci, C2v, C2h, D2d, D2h, C3, C6, S4, D2, T, Th, Td, O, Oh, I, Ih. Character tables from libmsym match Katzer data (comparison is order-independent via character fingerprints, handles different class orderings and irrep naming conventions). Irrational characters (Ih golden ratio) match to 3 decimal places.

3. ~~**Mass-weighting in displacement analysis.**~~ **Resolved.** Enum parameter `MassWeighting` on `symmetry_coordinates()`, saved on the `SymmetryCoordinates` result type so downstream code knows what it's working with.

4. ~~**PointGroup construction without molecules.**~~ **Resolved.** All point groups treated uniformly via libmsym — `PointGroup::from_schoenflies("C7v")` works the same as `from_schoenflies("C2v")`. No special-casing for crystallographic groups. libmsym handles arbitrary n.

5. ~~**Coordinates variant.**~~ **Resolved.** The `Symmetric` variant is removed from the `Coordinates` enum. The enum itself stays — it distinguishes coordinate types (Cartesian, internal, redundant internal, etc.), which is a separate concern. The `Symmetric` variant was orthogonal to that distinction and is no longer needed because: (a) every molecule now carries its `PointGroup` + equivalence sets + permutations, (b) the asymmetric unit is a derived view (one representative per equivalence set), (c) symmetry-adapted coordinates are produced on demand by `symmetry_coordinates()`, which returns both SALC definitions and evaluated values.
