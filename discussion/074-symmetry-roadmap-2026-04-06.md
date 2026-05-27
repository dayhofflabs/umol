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

### Existing implementations (surveyed 2026-04-06)

| Project | Language | License | Notes |
|---|---|---|---|
| [continuous-symmetry-measure/csm](https://github.com/continuous-symmetry-measure/csm) | Python | GPL-2 | Official Avnir-group implementation |
| [abelcarreras/symgroup](https://github.com/abelcarreras/symgroup) | Python | — | Independent CSM implementation |
| [abelcarreras/posym](https://github.com/abelcarreras/posym) | Python | — | Point symmetry analysis using CSM (normal modes, wave functions) |
| [abelcarreras/WFNSYM](https://github.com/abelcarreras/WFNSYM) | C/Python | — | CSM of electronic wave functions |
| [cosymlib](https://cosymlib.readthedocs.io/) | Python | — | CSM + continuous shape measures |

**License constraint**: official csm is GPL-2, incompatible with permissive licensing. Algorithms must be reimplemented from papers, not ported.

**Key papers**:
- [CSM Software, J. Chem. Inf. Model. 2024](https://pubs.acs.org/doi/10.1021/acs.jcim.4c00609) ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11267602/)) — current Avnir-group implementation, revised algorithms for speed/accuracy
- [Approximate algorithms for large structures, 2023 (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10636902/) — relevant if exact permutation search is too expensive

### Algorithm details (from Avnir group papers, 2023–2024)

Core loop is **direction-permutation iteration**:
1. Initialize symmetry direction `ν_sym` (Cartesian axes or Fibonacci sphere sample)
2. Iterate until convergence:
   a. Build distance matrix `A_ij = ||T(Q_i) − Q_j||²` for symmetry operation T
   b. Find permutation π via Hungarian / greedy / structure-preserving search
   c. Analytically optimize `ν_sym` for current π
3. `CSM = 100 · M(G) / D` where `M(G) = (1/2n) · min Σᵢ Σₖ ||Tⁱ(Q_k) − Q_{πᵢ(k)}||²`

| Permutation search | Complexity | Structure preservation |
|---|---|---|
| Greedy | `O(N²)` per iteration | 63–90% |
| Hungarian (Munkres) | `O(N³)` per iteration | 86–90% |
| Structure-preserving recursive | `O(N³ · branches)` with pruning | 100% |
| Exact (all permutations) | exponential | 100% |

**Key observations:**
- For low CSM (< 2 units, typical use case), all methods converge to identical results — simple version is fine.
- Official Avnir software handles only cyclic groups (Cs, Ci, Cn, Sn) — **not Td/Oh/Ih**. Polyhedral group → subgroup is a niche umol could fill.
- For our primary use case (lower a detected high-symmetry structure to a specified subgroup), the orientation is already fixed by the parent group, so the orientation refinement loop may be unnecessary.

### Implementation scope

Phased approach:
1. **Phase 1**: greedy nearest-neighbor + same-element matching, no orientation optimization, parent-group axes assumed correct. Sufficient for the "lower symmetry to a known subgroup" use case. Reuses existing `compute_atom_permutations` machinery.
2. **Phase 2**: Hungarian assignment for guaranteed bijection.
3. **Phase 3**: Fibonacci direction sampling + iterative refinement, for the general case where the molecule is not pre-aligned.

## Cross-cutting: numerical primitives

The CSM algorithm motivates two general-purpose primitives that should not live inside the symmetry module:

### Hungarian algorithm (linear assignment)
Used by CSM for atom-to-image matching, but also relevant to:
- RMSD with optimal atom assignment (general structural comparison)
- Atom mapping between reactant/product in reaction coordinates
- Trajectory frame matching (MD analysis)
- Conformer alignment

Crate options: `pathfinding`, `lap`. Decision: use external crate, do not reimplement.

### Fibonacci sphere sampling
Used by CSM for direction initialization, but also relevant to:
- Solid angle integration (numerical quadrature)
- Initial guesses for any SO(3) optimization
- Powder spectrum simulation (crystal orientation averaging)
- Docking / pose sampling

Trivial to implement (≈10 lines). Decision: implement as a utility, location TBD (`umol-numerics`? new `umol-geometry-utils`?).

### Where to put them
Both primitives are general numerical tools, not symmetry-specific. Options:
- New crate `umol-numerics` for shared numerical utilities
- Add to existing `umol-msym/src/linear.rs` (currently linear-algebra helpers) and rename
- Inline in symmetry module initially, extract when a second consumer appears

Recommend: inline in symmetry module for now (YAGNI), extract when needed.

## Medium-term: subgroups and correlation

### Goals

The subgroup machinery covers four interlocking questions, not just lattice membership:

1. **Lattice structure** — which groups are subgroups of which, edge types (normal, setting, notation, basis change), embedding multiplicity.
2. **Site symmetries** — for a molecule with parent group G, the stabilizer subgroup of each atom (or each Wyckoff-like position) and the orbit decomposition.
3. **Irrep correlation under descent** — for G ⊃ H, how each parent irrep λ ∈ Irr(G) decomposes into Irr(H). The inverse direction (induction Ind_H^G) is also useful.
4. **Distortion-driven descent** — given a vibrational or electronic mode of irrep μ in G, find the maximal subgroup H ⊂ G in which μ becomes totally symmetric. This is the *epikernel* of μ.

Applications: Jahn–Teller analysis, soft-mode phase transitions, ligand-field splitting, NMR equivalence on lowered symmetry, isomer search, vibronic coupling, crystal-field analysis, descent under nuclear displacement.

### Reference: Altmann & Herzig, *Point-Group Theory Tables* (2nd ed.)

Local copy: `materials/symmetry/Altmann S. L., Herzig, Point-Graph Theory Tables - 2nd ed.pdf`

This is the comprehensive reference. Chapter 9 contains 12 subgroup-lattice graphs covering all common point groups (Cn, Sn, Dn, Dnh, Dnd, Cnv, Cnh up to n=10, plus T/Td/Th/O/Oh, I/Ih). Chapter 16 documents the per-group tables T 22 – T 75 with subgroup elements and correlation tables.

**Edge types in the graphs:**

| Style | Meaning |
|---|---|
| Solid | Invariant (normal) subgroup |
| Dash-dot | Change of setting (axis relabelling) |
| Dotted | Change of notation |
| Dashed | Subduction with basis change |
| Thick double | **Subduction may fail** — only in graphs 11 (cubic) and 12 (icosahedral) |

**Graph index:**

| Graph | Coverage |
|---|---|
| 1 | C6 / D6h family |
| 2 | C7 family |
| 3 | C8 / D8h family |
| 4 | C9 family |
| 5 | C10 / D10h family |
| 6 | D6d, S12 |
| 7 | D7d, S14 |
| 8 | D8d, S16 |
| 9 | D9d, S18 |
| 10 | D10d, S20 |
| 11 | Cubic: T, Td, Th, O, Oh |
| 12 | Icosahedral: I, Ih |

### "Subduction may fail" — the embedding ambiguity

For polyhedral groups, the same abstract subgroup can have inequivalent embeddings. Example: D2 ⊂ O has two distinct realizations (along C4 axes vs. along C2′ axes). A single correlation table cannot encode both. Altmann–Herzig flags these cases with the thick-double edge style and provides separate tables for each embedding when it matters.

**Design decision needed:** subgroup descent must accept an embedding selector for cubic and icosahedral parents. API sketch:

```rust
enum Embedding { Default, Named(&'static str) }
fn correlation(parent: PointGroup, child: PointGroup, embed: Embedding) -> CorrelationTable
```

For non-polyhedral parents, `Default` is unique.

### Encoding strategy

Three viable approaches, in increasing order of preference:

**(A) Hand-encode from Altmann–Herzig.** Type ~3000–4000 correlation entries from tables T 22 – T 75. Authoritative but typo-prone, ~7–9 days, high verification burden.

**(B) Compute via character restriction.** Encode each maximal subgroup as the index set of its parent operations. Then correlation tables are derived: restrict each parent irrep's character vector to the subgroup classes and `reduce()` against the subgroup's character table. We already have all character tables via libmsym. Polyhedral embedding ambiguity is handled naturally — different operation index sets for the same abstract subgroup yield different correlation tables. ~5–6 days.

**(C) Port Gernot Katzer's algorithms.** See next section. ~5–7 days but with validated subgroup enumeration logic instead of hand-derived rules.

Common type definitions for all three:

```rust
enum SubductionKind {
    Invariant,        // normal subgroup, identity basis
    Setting,          // axis relabelling needed
    Notation,         // pure notation difference
    BasisChange,      // explicit basis transformation
    MayFail,          // multiple inequivalent embeddings
}

struct SubgroupEdge {
    parent: PointGroup,
    child: PointGroup,
    kind: SubductionKind,
    basis_transform: Option<Mat3>,
    embedding: Option<&'static str>,
}

struct CorrelationTable {
    parent: PointGroup,
    child: PointGroup,
    rows: &'static [(Irrep, &'static [Irrep])],  // parent → reduction in child
}
```

### Reference implementation: Gernot Katzer's character tables

Local mirror: `materials/character_tables/` (HTML pages, `.lis` files, and `ptgroup.js`).

**Surprising structural fact**: Katzer's site does not store subgroup or correlation tables anywhere. The static `.html` and `.lis` files contain only character tables, classes, and l-symmetry / multipole decompositions. Everything subgroup-related is **computed at page load by `ptgroup.js`** running in the browser.

This is excellent news: the algorithms are already written, validated by ~25 years of chemist use, and small enough to port directly.

**Key functions in `ptgroup.js`:**

| Function | Lines | What it does |
|---|---|---|
| `make_subgroup(grpname)` | 3245– | Enumerates subgroups from the parent Schoenflies symbol via regex rules (e.g. `Cnh → Cmh` for m\|n; `Dnh → Dmd` under specific divisor conditions) |
| `make_subgroup_oct` | 3127– | Special-case enumeration for cubic groups (T, Td, Th, O, Oh) |
| `make_subgroup_ico` | 3197– | Icosahedral — explicitly **not implemented** (bails out at line 3131) |
| `calc_symmetry_reduction(subgroups)` | 4108– | Character restriction: matches subgroup ops to parent class indices, then finds subgroup irreps whose characters match the parent's restricted vector |
| `calc_distortion_subgroup` | called 4173 | Distortion-driven descent (epikernel enumeration) |

The character restriction loop in `calc_symmetry_reduction` is exactly approach (B). Katzer's added value over (B) is the **subgroup enumeration regex rules** in `make_subgroup`, which have been debugged against thousands of cases and are unlikely to have systematic gaps for finite Cn/Sn/Dn families.

**Effort, porting Katzer:**

| Task | Effort |
|---|---|
| Port `make_subgroup` regex rules to Rust (finite Cn/Sn/Dn families) | 1 d |
| Operation matching: subgroup ops → parent class indices | 0.5 d |
| Port `calc_symmetry_reduction` over our `CharacterTable` | 0.5 d |
| Port `calc_distortion_subgroup` (epikernel enumeration) | 1 d |
| Port `make_subgroup_oct` (cubic groups) | 0.5 d |
| **Icosahedral fallback** — hand-encode I/Ih edges from Altmann–Herzig | 1–2 d |
| Site symmetries (independent of Katzer) | 0.5–1 d |
| Tests, cross-check against A–H worked examples | 1 d |

**Total: 5–7 days.** The gain over approach (B) is risk reduction — Katzer's `make_subgroup` validated against ~25 years of use removes the "did I miss a subgroup" failure mode.

**Caveat — distortion descent ambiguity in `ptgroup.js`** (line 253–263 of `index.html`): for truly degenerate irreps in non-abelian groups, the two E components may distort to different (or the same) subgroups. Katzer documents three notations — `→ D3 C3v` (different), `→ [ 2 D3 ]` (same type, ambiguous embedding, linked to 4n-fold axes), and the standard single-target case. Our implementation needs the same three-way distinction in the `DescentBranch` API.

### Site symmetries (atomic stabilizers)

For each atom i in a molecule with point group G, the **site symmetry** S_i ⊆ G is the subgroup of operations that fix atom i (the stabilizer in group-theoretic terms). The orbit-stabilizer theorem gives

  |G| = |orbit(i)| · |S_i|

so site symmetries refine the equivalence-set partition we already compute.

**Current state.** `Molecule::atom_permutations()` gives the orbit decomposition (atoms in the same orbit are equivalent). What is missing: the *labelled* site subgroup for each orbit.

**Required addition.** For each orbit, compute the stabilizer of one representative atom by intersecting the operations that map it to itself, then identify which subgroup of G that stabilizer is isomorphic to (and which embedding, if the parent is polyhedral). API sketch:

```rust
struct AtomSite {
    orbit_id: usize,
    representative_atom: usize,
    stabilizer: PointGroup,
    embedding: Option<&'static str>,
}

impl Molecule {
    pub fn atom_sites(&self) -> Vec<AtomSite>
}
```

**Why it matters.**
- **NMR / EPR equivalence under descent**: when G drops to H, an orbit of G splits into orbits of H. The split is determined by how S_i ∩ H sits inside H. Predicts NMR multiplet patterns under symmetry breaking without redoing the full perception.
- **Local-mode analysis**: vibrational modes localized on a site transform as Ind_{S_i}^G(local irreps). For Jahn–Teller and vibronic coupling.
- **Functional group symmetry**: a substituent's local irreps determine its coupling to the rest of the molecule.

### Distortion-driven descent (epikernels)

Given an irrep μ of G with `quadratic_irreps()`, `translation_irreps()`, or vibrational classification, ask: *which subgroup H ⊂ G is preserved by a distortion along a basis vector of μ?*

**Definition.** The **kernel** of μ is the subgroup K(μ) = { R ∈ G : Γ_μ(R) = 1 } — operations that act as identity on the entire μ representation. The **epikernel** of a specific direction d in μ is the larger subgroup E(μ, d) = { R ∈ G : Γ_μ(R)·d = d } — operations that fix that particular direction.

**Operational meaning.** If a molecule distorts along a normal mode of irrep μ, the resulting structure has symmetry equal to the epikernel of the chosen direction within μ. Different directions in a multidimensional irrep give different (possibly inequivalent) epikernels — the basis for Jahn–Teller distortion paths and isomer enumeration.

**Algorithmic content.**
1. For each subgroup H ⊂ G, check whether the subduction Sub_H^G(μ) contains the totally symmetric irrep A1(H). If yes, H is an epikernel candidate for some direction in μ.
2. Maximal such H are the chemically interesting "first descent" subgroups.
3. For multidimensional μ, enumerate the inequivalent maximal epikernels — these correspond to distinct distortion isomers (Jahn–Teller minima, soft-mode wells).

**Why it's straightforward once correlation tables exist.** The check "does Sub_H^G(μ) contain A1(H)?" is one row lookup in the correlation table. Epikernel enumeration is one pass over the subgroup lattice. No new representation theory needed.

**Reference.** Jarić & Birman's epikernel principle (1977 onward) — maximal-epikernel directions are physically realized minima for soft-mode and Jahn–Teller systems. Altmann–Herzig's correlation tables provide all the data needed; the principle just selects rows.

API sketch:

```rust
struct DescentBranch {
    via_irrep: Irrep,
    direction_label: Option<&'static str>,  // e.g. "T2g→tetragonal" vs "T2g→trigonal"
    child: PointGroup,
    embedding: Option<&'static str>,
}

impl PointGroup {
    /// All maximal subgroups reachable by distortion along the given irrep.
    pub fn epikernels(&self, mode: Irrep) -> Vec<DescentBranch>;

    /// Full descent tree under all non-trivial irreps (Jahn–Teller landscape).
    pub fn distortion_descent(&self) -> Vec<DescentBranch>;
}
```

### Subduction from O(3)

Chapter 12 of Altmann–Herzig gives closed-form Wigner-D matrix elements (eqs. 25, 27, 30) and characters

  χʲ(φ) = sin((j+½)φ) / sin(φ/2)

This means subduction *from* O(3) — i.e., reducing arbitrary l harmonics under any point group — does not need per-l hard-coded tables. Compute χʲ at each class angle, then `reduce()` against the point group's character table. Useful for:

- Crystal-field splitting of arbitrary l manifolds (not just d and f)
- Spherical-tensor decomposition of multipole moments
- Selection rules for higher-l states (relativistic, Rydberg)

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

For systems with spin-orbit coupling. The point group is extended by Ẽ (a 2π rotation, Ẽ² = E) which commutes with all g ∈ G, doubling the order. Spinor representations are the additional irreps in G̃ beyond the vector irreps of G.

Applications: heavy-element chemistry, spin-forbidden transitions, relativistic electronic structure, magnetic anisotropy in single-molecule magnets.

### Construction recipe (Altmann–Herzig Chapter 10)

Opechowski's theorem on class structure:

| Class type in G | Effect in G̃ |
|---|---|
| Regular | Splits: C(g) → C̃(g) + C̃(g̃) |
| Irregular (bilateral binary rotations + one of orthogonal-mirror pair) | Stays single: C̃(g) ≡ C̃(g̃) |

**Number of spinor irreps** = number of regular classes of G.
**Sum of squared dimensions** = |G| (same as vector irreps).

This is a pure character-table operation built on the existing `CharacterTable` infrastructure — no libmsym extension needed. Steps:

1. Classify each operation as regular or irregular (geometric test on the operation matrix).
2. Build G̃'s class structure from G's via Opechowski.
3. Construct projective representation of G with factor system from SU(2) (Chapter 11).
4. Diagonalize / reduce to get spinor irreps.

Two equivalent framings (Chapter 10 §2):
- **Double-group**: doubled operation set, vector representations of G̃. Simpler conceptually, multiplication rules not uniquely defined.
- **Projective**: original group, matrices that close up to a phase factor. Cleaner algebra, requires projective-rep machinery.

Altmann–Herzig recommends projective; both yield the same spinor characters and selection rules.

### SU(2) ↔ SO(3) (Chapter 11)

Cayley–Klein parameters give the explicit SU(2) matrix for any rotation R(φ n̂):

  a = cos(φ/2) − i nz sin(φ/2),  b = −(ny + i nx) sin(φ/2)

Pauli gauge: i ∈ O(3) maps to the identity SU(2) matrix. This fixes the sign convention for spinor operations across all centrosymmetric groups.

### Selection rules with spinors

Once double-group character tables are built, spin-forbidden transitions, ZFS tensors, and SOC matrix elements reduce to the same `direct_product` / `reduce` machinery already in place. The `is_gerade` infrastructure for centrosymmetric groups carries over (eq. 33 of Chapter 13: complex-conjugate spinor bases are u when ordinary spinors are g, and vice versa — already aligned with our real-irrep treatment).

## Not planned
- Space groups (not relevant for molecular symmetry)
