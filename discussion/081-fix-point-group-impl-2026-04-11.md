# Fix PointGroup: one libmsym-backed source of truth, clean abstract / detection split

**Status:** implemented. 203 umol-msym + 42 umol-geometric tests pass against the rewritten types. Remaining warnings are pre-existing `dead_code` on unused `Context` methods.

## 1. Problem

Two tangled defects in `umol-msym`:

1. **`PointGroup` leaks a frame-dependent 3D matrix representation through its allegedly abstract API.** The singleton cached under a `SchoenfliesLabel` is populated by whoever asks first — either by libmsym reading a real molecule (`from_context`) or by `from_schoenflies::construct` feeding libmsym a hardcoded pair of fake atoms. The character table is the same either way, but the operation matrices and axis vectors are not: they belong to some arbitrary reference frame. Every downstream caller that reaches `group.ops()[i].matrix` is silently using matrices from an unrelated molecule.

2. **`SymmetryOp` is underdefined as a label.** In C₂ᵥ, the two mirror planes σᵥ and σᵥ′ share the same `(kind, order, power, orientation) = (Reflection, 1, 1, Vertical)`. They are distinct group elements (each a singleton class in the order-4 group), but nothing in the current `SymmetryOp` type distinguishes them. A `SymmetryOp` value therefore cannot be used as a unique key into anything — not a multiplication table, not a list of representation matrices, not a character row index.

`Molecule.ops: Vec<SymmetryOp>` was added as a local workaround for #1 and must be removed: frame-dependent matrices do not belong on `Molecule`, and the workaround only papered over #2 by caching the "right" matrices alongside ambiguous labels.

## 2. Mathematical foundation

- **Group identity is positional, not descriptive.** Elements of a finite group G are distinguished by their position in the group, not by any tuple of attributes. Human-readable labels like "σᵥ" are annotations on top of an index, not substitutes for one. In C₂ᵥ, σᵥ and σᵥ′ are distinct elements; the only reason they share a descriptive tuple is that our descriptive vocabulary is too coarse.

- **Orientation-independent vs orientation-dependent data must not be mixed.**

  Orientation-independent (lives on the abstract group):
  - Cardinality, conjugacy-class partition, multiplication table
  - Character table, irrep labels and dimensions
  - Which irreps carry translation / rotation / quadratic basis functions
  - Subgroup lattice, correlation tables given a subgroup embedding map

  Orientation-dependent (lives in a 3D realization):
  - 3×3 matrices implementing each group element
  - Rotation axis / reflection normal unit vectors
  - The permutation of atoms induced by each operation (depends on atom labelling)
  - Cartesian projection operators used in symmetry-adapted coordinate analysis

The job of the type system is to put the first set on `PointGroup` and the second set on a separate per-molecule object.

## 3. libmsym already provides the abstract data

libmsym has an element-free path for constructing a point group and emitting its abstract structure. The relevant calls:

- `msymSetPointGroupByName(ctx, "C2v")` — `libmsym/src/msym.c:83`. Takes a context and a name string. Never calls `ctxGetElements`. Internally calls `generatePointGroupFromName` — `libmsym/src/point_group.c:106` → `generatePointGroupFromStruct` — `libmsym/src/point_group.c:118` → `generateSymmetryOperations(type, n, order, &sops)` — `libmsym/src/point_group.c:1186`. `generateSymmetryOperations` is a pure function of `(type, n, order)`; it produces the canonical operation list in a canonical orientation, with no molecule, no equivalence sets, no geometry.

- `msymGetSymmetryOperations(ctx, &len, &sops)` — `libmsym/src/context.c:398`. Guards only on `ctx->pg->sops`. Works immediately after `msymSetPointGroupByName`.

- `msymGetCharacterTable(ctx, &ct)` — `libmsym/src/context.c:323`. Lazily calls `generateCharacterTable(pg->type, pg->n, pg->order, pg->sops, &pg->ct)`. No element dependency. Returns irrep names, dimensions, class counts, and the full character table.

The current `from_schoenflies::construct` trick that feeds libmsym fake atoms `[1.0, 0.3, 0.7]`, `[0.5, 0.8, 0.2]` to coax a character table out is simply not the supported way. The supported way is the name-based API above. It takes five libmsym calls and returns everything the abstract layer needs.

### Canonical op ordering is preserved through detection

In the detection path, `msymFindSymmetry` → `findPointGroup` — `libmsym/src/point_group.c:410` → `generatePointGroup` — `libmsym/src/point_group.c:330`. `generatePointGroup` at line 349 calls **the same** `generateSymmetryOperations(type, n, pg->order, &pg->sops)`. Lines 360–365 then rotate the axis vectors into the molecule frame via the precomputed `transform`, but the op list itself is regenerated canonically in canonical order.

Consequence: for fixed `(type, n)`, libmsym produces ops in a deterministic order regardless of whether you came in via a name or via molecule detection. Canonical construction and per-molecule detection yield ops at the *same indices*; they differ only in the axis vectors and the 3×3 matrices. Matching a `SymmetryOp` in the abstract layer to its per-molecule matrix is therefore an index lookup, with a debug assertion checking `(type, order, power, class)` agreement at each index.

## 4. Layer split

The split is module-level policy, not a second crate. Both layers link against `umol-msym-sys`.

**Abstract layer** — `point_group.rs`, `subgroup.rs`, `linear.rs`:
- Owns `PointGroup`, `SymmetryOp`, `Irrep`, the `REGISTRY` of cached singletons.
- Owns multiplication, reduction, direct products, symmetric/antisymmetric squares, correlation tables, chirality/gerade queries.
- Calls libmsym in exactly one place: `PointGroup::from_schoenflies(label)`, which opens a minimal context, calls `msymSetPointGroupByName`, extracts ops + character table into owned Rust data, releases the context. After that call returns, no libmsym pointer is retained anywhere on the singleton.
- Never holds a molecule, an equivalence set, or any frame-dependent data.

**Detection layer** — `context.rs`, `detect.rs`, `matrix_rep.rs`, `basis.rs`:
- Owns `Context` (the FFI wrapper), `SymmetryResult`, `SymmetryDescentResult`, `MatrixRep`, SALC construction.
- Calls libmsym for molecule perception, subgroup descent, SALC/basis-function machinery.
- Produces `SymmetryResult` values that reference `&'static PointGroup` from the abstract layer plus a per-molecule `MatrixRep`.
- `Context` owns its `ffi::msym_context` and never exposes it. Callers obtain a `MatrixRep` via `Context::symmetry_representation(group)`, which reads libmsym's op list and hands back owned Rust data. No FFI handle leaks across the module boundary.

Linear groups (`C∞v`, `D∞h`) stay in `linear.rs` with their existing hand-rolled tables; libmsym has no finite character table for them and there is nothing to extract.

### 4.1 Why the registry is correct once its contents are abstract

The registry itself is a reasonable idea: there *is* one C₂ᵥ. Every construction of the abstract C₂ᵥ — whether triggered by detecting water, detecting formaldehyde, or descending from C₂ᵥₕ — produces byte-identical `op_data`, `classes`, `mul_table`, `irrep_data`. Caching one `&'static PointGroup` per `SchoenfliesLabel` is the correct representation of that fact.

What makes the current registry broken is not the pattern, it's the contents. Today, a `PointGroup` stores `ops[i].matrix` and `ops[i].vector` from whichever libmsym context built it first. Those fields are frame-dependent, so first-writer-wins *does* leak into observable behaviour: subsequent callers of a given label get matrices from an unrelated molecule's orientation. The `select_subgroup` path hits this particularly hard, because descent starts from a libmsym subgroup context whose axis alignment depends entirely on the parent molecule.

Once matrices and axis vectors are removed from `PointGroup` (§5.2) and relocated to per-molecule `MatrixRep` (§5.3), the cached content is invariant across all construction sites. Whichever molecule happens to populate the singleton first, every other caller would have populated it identically. The race becomes harmless by construction, the subgroup descent bug disappears, and the lock is contended only on first-registration of each label (≤ 40 labels per process, one lock acquisition each).

Two notions of "same group" remain, and they answer different questions:

- **Same abstract group** — `a.group.label() == b.group.label()`. True for water's and formaldehyde's C₂ᵥ handles, because the registry hands both molecules the same singleton. Used when asking whether a character table, irrep symbol set, or reduction formula from one analysis applies to another.
- **Same singleton instance** — `std::ptr::eq(a.group, b.group)`. Equivalent to label equality in the presence of the registry — this is exactly what "one C₂ᵥ" means — and remains the fast path for `SymmetryOp::eq`. Used as the basis for op-index identity when matching `SymmetryOp` values to `MatrixRep` slots.

## 5. Type design

### 5.1 `SymmetryOp` is a handle into its parent group

```rust
#[derive(Clone, Copy)]
pub struct SymmetryOp {
    group: &'static PointGroup,
    index: usize,
}
```

`SymmetryOp` is to `PointGroup` what `Irrep` already is: a lightweight handle carrying a pointer to the cached singleton plus a positional index. Equality is `ptr::eq` on the group plus `==` on the index. This makes `SymmetryOp` a legitimate unique key:

- σᵥ in C₂ᵥ and σᵥ in D₂ₕ compare unequal because their `group` pointers differ.
- σᵥ and σᵥ′ in C₂ᵥ compare unequal because their indices differ, even though every descriptive attribute agrees.

Descriptive attributes are accessors that read from the parent group's op data:

```rust
impl SymmetryOp {
    pub fn group(&self)       -> &'static PointGroup { self.group }
    pub fn kind(&self)        -> SymmetryOpKind      { ... }
    pub fn order(&self)       -> i32                 { ... }
    pub fn power(&self)       -> i32                 { ... }
    pub fn orientation(&self) -> SymmetryOpOrientation { ... }
    pub fn class(&self)       -> usize               { ... }
    pub fn character(&self, irrep: Irrep) -> f64     { ... }
}

impl PartialEq for SymmetryOp {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.group, other.group) && self.index == other.index
    }
}
```

No `matrix`, no `vector`, no `transform_point` on `SymmetryOp`. Those require a 3D realization and live on the type that owns one (§5.3).

The `usize` inside is an opaque identity key. External code obtains `SymmetryOp` values from `group.ops()`, `group.class_reps()`, semantic accessors on `PointGroup` (§5.2), or multiplication-table lookups; it does not construct `SymmetryOp` by writing an integer literal.

### 5.2 `PointGroup` keeps its ops

`PointGroup` stores the full abstract group. Nothing gets removed.

```rust
pub struct PointGroup {
    label: SchoenfliesLabel,
    order: usize,
    op_data: Vec<OpData>,              // index -> kind/order/power/orientation/class
    classes: Vec<Vec<usize>>,           // class_index -> sorted op indices
    class_rep_indices: Vec<usize>,      // class_index -> chosen op index
    mul_table: Vec<Vec<usize>>,         // [i][j] = index of op_i · op_j
    irrep_data: Vec<IrrepData>,         // symbol, dimension, character row (by class)
}
```

Public accessors:

- `ops(&'static self) -> impl Iterator<Item = SymmetryOp>` — every element, once.
- `op(&'static self, index: usize) -> SymmetryOp` — internal use and testing.
- `class_reps(&'static self) -> impl Iterator<Item = SymmetryOp>` — one op per class. Ordinary `SymmetryOp` values, same type as anything from `ops()`. Used for character-table display and reduction formulas.
- `multiply(a: SymmetryOp, b: SymmetryOp) -> SymmetryOp` — reads `mul_table`. Both arguments must belong to `self` (debug assertion).
- `irreps(&'static self) -> Vec<Irrep>` — unchanged.
- `reduce`, `direct_product`, `symmetric_square`, `antisymmetric_square` — unchanged in meaning; characters reached through the character table indexed by class.
- Semantic accessors: `identity()`, `inversion()` (where present), `principal_rotations()` (ops with kind=Proper, order=n_max), etc., as needed by callers currently grepping through ops for a specific element.

### 5.3 `MatrixRep` is a separate per-molecule object

A (faithful) 3D matrix representation of G is a homomorphism ρ : G → O(3). Concretely it is a vector of 3×3 matrices and axis vectors, one per group element, indexed by `SymmetryOp::index`. "Matrix representation" is the standard group-theory name for this, and it denotes the *whole homomorphism*, not a single matrix. The type name is shortened to `MatrixRep` at call sites.

```rust
pub struct MatrixRep {
    group: &'static PointGroup,
    matrices: Vec<Matrix3<f64>>,   // matrices[op.index()] = ρ(op)
    axes: Vec<Vector3<f64>>,       // axes[op.index()] = rotation axis / reflection normal
}

impl MatrixRep {
    pub fn group(&self) -> &'static PointGroup { self.group }
    pub fn order(&self) -> usize { self.matrices.len() }
    pub fn matrix(&self, op: SymmetryOp) -> &Matrix3<f64> {
        assert!(std::ptr::eq(op.group(), self.group));
        &self.matrices[op.index()]
    }
    pub fn axis(&self, op: SymmetryOp) -> &Vector3<f64> { ... }
    pub fn transform_point(&self, op: SymmetryOp, p: Vector3<f64>) -> Vector3<f64> { ... }
    pub fn matrices(&self) -> &[Matrix3<f64>] { ... }
    pub fn axes(&self) -> &[Vector3<f64>] { ... }
    pub fn iter(&self) -> impl Iterator<Item = (SymmetryOp, &Matrix3<f64>, &Vector3<f64>)> + '_ { ... }
    pub fn identity_only(group: &'static PointGroup) -> Self { ... }
}
```

Construction happens inside the detection layer via `Context::symmetry_representation(group)`, which calls `msymGetSymmetryOperations` and reconstructs each 3×3 matrix from the libmsym sop tuple `(type, order, power, axis)`. Because libmsym generates ops canonically (§3), the `i`-th libmsym op corresponds to the `i`-th `op_data` entry in the abstract singleton. Construction asserts the op count matches `group.order()` and copies matrices and axes into owned storage. `MatrixRep` holds no FFI pointer.

`identity_only` produces a C₁ realization directly from a `&'static PointGroup`, used by the `Molecule::from_parts` default constructor and by the C₁ descent path.

## 6. Data locations on the caller side

### 6.1 `SymmetryResult` bundles the result of perception

```rust
pub struct SymmetryResult {
    pub group: &'static PointGroup,
    pub representation: MatrixRep,
    pub equivalence_sets: Vec<EquivalenceSet>,
    pub centers: Vec<SymmetryCenter>,
}
```

`ops: Vec<SymmetryOp>` goes away from `SymmetryResult`; the same information is reachable as `group.ops()` (abstract handles) + `representation.matrix(op)` / `representation.axis(op)` (frame-dependent data).

`SymmetryDescentResult` gains `child_representation: MatrixRep` in place of `child_ops: Vec<SymmetryOp>`.

### 6.2 `Molecule` carries `MatrixRep` instead of `Vec<SymmetryOp>`

```rust
pub struct Molecule {
    elements: Vec<Element>,
    coordinates: Coordinates,
    charge: i32,
    multiplicity: SpinMultiplicity,
    group: &'static PointGroup,
    representation: MatrixRep,
    equivalence_sets: Vec<Vec<usize>>,
    atom_permutations: Vec<Vec<usize>>,   // indexed by op.index()
}
```

`ops: Vec<SymmetryOp>` is gone. `representation: MatrixRep` replaces it and is the sole carrier of frame-dependent data. `atom_permutations` stays: once computed, the permutation `atom i → atom π(i)` under a given op is invariant under rigid-body rotation and is a property of the labelled atoms, not of the coordinate frame. The permutation vector is indexed by op index — the same index that selects matrices from `representation`.

`symmetry_coordinates` stays on `Molecule` and reads `self.representation`. The earlier sketch considered moving it to `SymmetryResult`; in practice the method needs `{group, representation, atom_permutations, coordinates}` as a bundle, all of which are already owned by `Molecule`, and the analysis is naturally phrased as a method on a symmetrized molecule. Moving it would require `SymmetryResult` to duplicate coordinates and atom permutations or take them as parameters, which is less clean.

## 7. What the abstract layer gets from libmsym, concretely

`PointGroup::from_schoenflies(label)` does:

1. `msymCreateContext()`, `msymSetPointGroupByName(ctx, label_string)`.
2. `msymGetSymmetryOperations(ctx, &len, &sops)` → canonical ops. For each `i in 0..len`, copy `(type, order, power, orientation, class)` into `op_data[i]`.
3. Bucket op indices into `classes` by the `class` field. Pick one index per class for `class_rep_indices`.
4. Compute `mul_table` by matrix-multiplying the canonical 3×3 matrices libmsym handed us in-context, matching the product to an existing canonical op by `(type, order, power, orientation, class)` plus axis-vector agreement. Write back the index. Matrices are scratch; they are discarded after the table is built.
5. `msymGetCharacterTable(ctx, &ct)` → irrep names, dimensions, class counts, character values. Copy into `irrep_data`.
6. `msymReleaseContext(ctx)`.

After step 6 the `PointGroup` owns all its data and never touches libmsym again. There is one source of truth (libmsym's generator tables); the Rust side caches its output.

## 8. What goes away

- `PointGroupKind::ops` and `class_reps` as separately stored `SymmetryOp` values — replaced by `op_data` indexed by position, with `ops()` and `class_reps()` returning handle views.
- `SymmetryOp { vector, matrix }` fields.
- `SymmetryOp::compute_matrix`.
- `from_schoenflies::construct`'s hardcoded seed atoms `[1.0, 0.3, 0.7]`, `[0.5, 0.8, 0.2]`. Replaced by the name-based libmsym API.
- `PointGroup::ops()` as a source of frame-dependent matrices (it becomes purely abstract).
- `Molecule.ops`.
- `SymmetryResult.ops` (becomes `representation`).
- `SymmetryDescentResult.child_ops` (becomes `child_representation`).
- The identity-descent / C₁-descent paths in `lower_symmetry` that clone `group.ops()` into frame-dependent matrices — they now clone the parent `representation` instead, or construct an identity `MatrixRep` from the molecule's frame.

## 9. Implementation strategy

Rewrite the core in place; preserve the periphery that already works. Not a full crate rip, not an incremental migration with parallel types behind a flag.

### 9.1 Rewrite from scratch

- **`types.rs`** — redefine `SymmetryOp` as `{ group: &'static PointGroup, index: usize }`. Collapse `PointGroupKind::{ops, class_reps}` into a single `op_data: Vec<OpData>` store plus `class_rep_indices: Vec<usize>`. Drop `SymmetryOp::{matrix, vector}` and `compute_matrix`.
- **`point_group.rs`** — `PointGroup` internals, `REGISTRY`, `from_schoenflies` rewritten on top of `msymSetPointGroupByName`. Delete the seed-atom `construct`. Port `reduce`, `direct_product`, `symmetric_square`, `antisymmetric_square`, `translation_irreps`, `rotation_irreps`, `quadratic_irreps`, `r_squared_classes`, `is_chiral`, `has_inversion`, `Irrep::is_gerade` to read class-indexed characters — none need matrices. Rewrite `CharacterTableDisplay` on the handle view of `class_reps()`.
- **`MatrixRep`** — new type in its own module `matrix_rep.rs`. FFI-free: construction takes `(group, matrices, axes)` and is driven by `Context::symmetry_representation(group)`.

### 9.2 Mechanical rewrite — same logic, new types

- **`detect.rs`** — `SymmetryResult.ops` → `representation: MatrixRep`; same for `SymmetryDescentResult.child_ops` → `child_representation`. Identity-descent and C₁-descent paths build `MatrixRep::identity_only(group)` instead of cloning matrices off `group.ops()`. The `lower_symmetry` class-mapping loop iterates `child_group.ops()` and matches each `child_representation.matrix(child_op)` against the parent-op matrices carried on `SubgroupInfo.parent_ops`, which is now `Vec<(Matrix3<f64>, usize)>`.
- **`context.rs`** — the old `symmetry_operations()` and `character_table()` methods are gone. Their replacement is `Context::symmetry_representation(&self, group: &'static PointGroup) -> Result<MatrixRep, Error>`, which reads libmsym's op list via `msymGetSymmetryOperations`, reconstructs the 3×3 matrices via `compute_op_matrix`, and returns a `MatrixRep`. `Context` never exposes its `ffi::msym_context`.
- **`umol-geometric/molecule.rs`** — remove `ops` field; add `representation: MatrixRep`. `compute_atom_permutations` takes `&MatrixRep` and iterates `representation.matrices()`. `symmetry_coordinates` reads `self.representation` and accesses matrices via `self.representation.matrix(op)` and classes via `op.class()` (now a method, not a field).

### 9.3 Preserve — do not touch beyond what the type rename forces

- **`linear.rs`** — libmsym has no finite character table for C∞ᵥ / D∞ₕ; the hand-rolled Λ-indexed machinery is correct and cannot be replaced by libmsym calls.
- **`subgroup.rs`** — `CorrelationTable` algebra is already abstract; only the types it imports change.
- **`basis.rs`** — SALC / `BasisFunction` types are orthogonal to the point-group cleanup.
- Existing test fixtures and `#[rstest]` case tables — expected values (class sizes, irrep counts, characters, reduction outputs) are correct and become the regression safety net for the rewrite.

### 9.4 Atomicity

The core type change lands in one commit series, not behind a feature flag. No "old `SymmetryOp` and new `SymmetryOp` coexist" — that doubles the surface area for no benefit. The ~2-week broken-tree window swallows the intermediate states. Tasks #8 (`SubgroupInfo` move) and #9 (`orbital_name` move) run as pre-rewrite cleanup so `context.rs` is already clean when the mechanical rewrites touch it.
