# Aromaticity Perception Design

Date: 2026-03-11

## 1. Survey of existing approaches

### 1.1 RDKit

RDKit's `setAromaticity` (in `Aromaticity.cpp`) dispatches on an `AromaticityModel` enum with six
variants: `AROMATICITY_RDKIT` (alias `AROMATICITY_DEFAULT`), `AROMATICITY_SIMPLE`,
`AROMATICITY_MDL`, `AROMATICITY_MMFF94`, `AROMATICITY_CUSTOM`, and an implicit default.

All four concrete models follow the same algorithmic structure:

1. Compute SSSR via `symmetrizeSSSR`.
2. Classify each ring atom's electron donor type (`ElectronDonorType`: Vacant, One, Two, OneOrTwo,
   Any, No) using `getAtomDonorTypeArom` and filter candidates via `isAtomCandForArom`.
3. For each ring (and, when enabled, fused ring combinations built via `makeRingNeighborMap` /
   `pickFusedRings` / `applyHuckelToFused`), sum the min/max electron contributions and check
   whether any value in the range satisfies 4n+2.
4. Mark atoms and bonds as aromatic.

The models differ only in their parameterization of this pipeline:

- **DEFAULT/RDKIT**: Elements through row 3 plus Se and Te. Exocyclic double bonds to
  more-electronegative atoms "steal" one electron from the ring atom's contribution. Fused ring
  combinations checked up to `maxNumFusedRings`. No ring size limits. At most one
  `AnyElectronDonorType` atom per ring.
- **SIMPLE**: Same algorithm restricted to ring sizes 5-6 with fused combinations disabled.
- **MDL**: C and N only. One-electron donors only. No third-row elements, no triple bonds, no
  exocyclic multiple bonds. Minimum ring size 6. Maximum 6 fused rings with maximum bond overlap 1.
- **MMFF94**: Force-field-specific. Uses sp2 hybridization checks, iterative convergence for fused
  systems, N/O/S special handling. Same 4n+2 gate at the end.

The `CUSTOM` variant accepts a function pointer, delegating entirely to user code.

### 1.2 CDK

CDK's modern aromaticity API (post-May & Steinbeck 2014 refactor) factors the problem into two
orthogonal strategy objects:

```java
new Aromaticity(ElectronDonation.daylight(), Cycles.all())
```

- `ElectronDonation` — how many pi-electrons each atom contributes. Implementations: `cdk()`,
  `daylight()`, `piBonds()` (most conservative: only count electrons from explicit pi bonds, ignore
  lone pairs).
- `CycleFinder` — which rings to check. Implementations: `Cycles.all()`, `Cycles.relevant()`,
  `Cycles.mcb()`, with composable fallback via `Cycles.or(a, b)`.

The 4n+2 check is fixed in the `Aromaticity` class. The factoring is clean, but the
`ElectronDonation` implementations are ad-hoc functions that recompute electron counts from atomic
properties, duplicating logic that belongs in a valence model.

Reference: May & Steinbeck, "Efficient ring perception for the Chemistry Development Kit," J.
Cheminform. 6:3 (2014). doi:10.1186/1758-2946-6-3.

### 1.3 OEChem / PubChem

OEChem provides five aromaticity models: `OEAroModelMDL`, `OEAroModelTripos`, `OEAroModelMMFF`,
`OEAroModelDaylight`, `OEAroModelOpenEye`. All are ring-based with 4n+2 electron counting; the
variation is in candidate filtering and electron counting rules.

OEChem is explicitly "SSSR-free" — it uses a pi-subgraph approach: extract the subgraph of atoms
that could participate in aromaticity, find connected components, detect cycles locally within each
component. This decouples aromaticity from global ring perception.

PubChem uses OEChem's `OEAroModelOpenEye` with max aromatic cycle path length 40 for structure
standardization. Aromaticity is perceived after tautomer canonicalization, followed by Kekulization.
PubChem's comparison of five aromaticity models (Hahnke, Kim, Bolton 2018, Fig. 5) shows significant
disagreement across models for many structures — confirming that model choice is not cosmetic.

Reference: Hahnke, Kim, Bolton, "PubChem chemical structure standardization," J. Cheminform. 10:36
(2018). doi:10.1186/s13321-018-0293-8.

### 1.4 Observation: one algorithm, three genuinely different approaches

Every mainstream cheminformatics aromaticity implementation is a parameterization of the same
algorithm: filter candidate atoms, obtain a ring set, check 4n+2. The 4n+2 check itself is
invariant across all models. The variation is in three dimensions: (a) which atoms are candidates,
(b) how many electrons each contributes, (c) which rings are tested. This algorithm will be called
**Hueckel rule aromaticity** to distinguish it from Hueckel MO theory, whose name it borrows without
performing any MO calculation.

Two genuinely different approaches exist but are not widely implemented in cheminformatics toolkits:

- **Hueckel Molecular Orbital (HMO) aromaticity**: an actual eigensolve on the pi-adjacency matrix
  (Hueckel Hamiltonian). Produces orbital energies, fractional bond orders, and delocalization energy
  as a quantitative aromaticity descriptor. Not ring-based — operates on the entire conjugated
  subgraph at once, handling fused systems naturally without ring-combination enumeration. The
  computational cost (O(n^3) eigensolve, n = number of pi-atoms) is negligible for real molecules
  (sub-millisecond for n < 30). Heteroatom handling through parametric Coulomb and resonance
  integrals (Van-Catledge 1980). The implementation complexity is comparable to or less than the
  thicket of empirical rules in RDKit's `isAtomCandForArom` + `getAtomDonorTypeArom` +
  `countAtomElec` — a 20-line eigensolve replaces hundreds of lines of special-case logic.
- **Clar's pi-sextet rule** (Liu & Green, Proc. Combust. Inst. 2019): a global ILP optimization
  maximizing disjoint pi-sextets rather than per-ring 4n+2 checking. Applicable to benzenoid PAHs
  only. Produces differential aromaticity across rings in fused systems.

All three approaches are planned as implementations of a common `AromaticityModel` trait:

| Model                      | What it computes                | Ring detection needed? | Scope              |
| -------------------------- | ------------------------------- | ---------------------- | ------------------- |
| `HueckelRuleAromaticity`   | 4n+2 electron count per ring    | Yes (cycles)           | All ring systems    |
| `HmoAromaticity`           | Eigensolve, bond orders, dE     | No (pi-subgraph only)  | All conjugated      |
| `ClarAromaticity`          | Maximum disjoint pi-sextets     | Yes (ring topology)    | Benzenoid PAHs only |

## 2. Ring scope

### 2.1 SSSR is the wrong abstraction

The Smallest Set of Smallest Rings (SSSR), also known as the Minimum Cycle Basis (MCB), is
non-unique. In cubane, 6 equivalent square faces exist but an MCB must pick 5 of 6; no canonical
choice exists without breaking symmetry. The choice of MCB can affect downstream aromaticity results
when the selected basis happens to omit a ring that should be checked.

No downstream chemical task actually requires a minimum cycle basis. The questions chemistry asks are
always more specific: "is this bond in a ring?", "what is the smallest ring through this atom?",
"which atoms share a delocalized pi system?" All are answerable with simpler primitives.

RDKit nominally computes `symmetrizeSSSR` but then builds fused ring combinations on top of it,
effectively working around the SSSR's limitations. OEChem is explicitly SSSR-free. The CDK paper
demonstrates that uniquely defined cycle sets (essential and relevant cycles) are computable in time
comparable to MCB.

Reference: Berger, Flamm, Gleiss, Leydold, Stadler, "Counterexamples in chemical ring perception,"
J. Chem. Inf. Comput. Sci. 44:323 (2004). doi:10.1021/ci030405d.

### 2.2 Two viable strategies

Two ring detection strategies have distinct trade-offs. Both will be implemented behind a
`RingStrategy` enum and compared empirically.

**PiSubgraph**: Extract the induced subgraph of atoms whose `AtomTypeSpec` candidates include at
least one with `aromatic_valence > 0`. Find connected components of this subgraph. Enumerate cycles
locally within each (typically small) component. This skips cycle enumeration in non-aromatic ring
systems entirely — a real saving for molecules like paclitaxel with complex fused saturated rings.

**GlobalAllCycles**: Compute biconnected components (Hopcroft-Tarjan). Enumerate all elementary
cycles per component up to `max_ring_size`. One computation serves aromaticity, stereochemistry, and
future descriptor consumers.

Comparison molecules: paclitaxel, steroids (fused non-aromatic); coronene, porphyrins (large fused
aromatic); camptothecin (mixed); cubane (pathological).

### 2.3 Shared ring primitives

Both strategies share a set of ring primitives in `rings.rs`:

- `biconnected_components(graph)` — Hopcroft-Tarjan. Partitions into ring systems. Answers
  `is_cyclic` in O(V+E).
- `smallest_ring_size(bond)` — local BFS. Needed for stereo-in-ring checks.
- `enumerate_cycles(subgraph, max_size)` — used locally by PiSubgraph on pi-components, or globally
  by GlobalAllCycles on biconnected components.

## 3. The three aromaticity models

### 3.1 Hueckel rule (4n+2 electron counting)

The Hueckel 4n+2 rule states that a planar, fully conjugated monocyclic system is aromatic if the
number of pi-electrons is 4n+2 for some non-negative integer n (i.e., 2, 6, 10, 14, ...).

In cheminformatics practice, the rule is applied as follows:

1. For each candidate ring (and, when configured, fused combinations of rings), determine the
   pi-electron contribution of each atom. Contributions depend on element, charge, bond environment,
   and lone pair availability.
2. Sum the contributions. If multiple atoms have ambiguous contributions (e.g., a nitrogen that could
   donate 1 or 2 electrons), enumerate the feasible combinations.
3. If any feasible sum satisfies 4n+2, the ring (or fused system) is aromatic.

The check is: `(electron_count >= 2) && ((electron_count - 2) % 4 == 0)`.

All cheminformatics aromaticity models implement this identical check. The variation across models is
entirely in step 1 (which atoms, how many electrons) and in which rings reach step 2.

The name `HueckelRuleAromaticity` distinguishes the simple counting rule from Hueckel Molecular
Orbital theory. The rule is a consequence of HMO theory applied to monocyclic systems, but the
cheminformatics implementation performs no MO calculation — it counts electrons and checks
divisibility.

### 3.2 Hueckel Molecular Orbital (HMO) aromaticity

An actual semiempirical quantum-chemical calculation on the pi-subgraph:

1. **Build the Hueckel Hamiltonian** H on the pi-adjacency graph. Diagonal elements: alpha (or
   alpha + h_X * beta for heteroatom X). Off-diagonal: beta for bonded C-C pairs (or k_CX * beta
   for C-X bonds). Parameters from Van-Catledge (1980) for C, N, O, S.
2. **Diagonalize** H. Eigenvalues = orbital energies. O(n^3) where n = number of pi-atoms.
3. **Fill electrons** by aufbau into lowest orbitals. Electron count from candidate `a` values.
4. **Compute delocalization energy**: dE = E_pi - E_ref, where E_ref = n_double * 2 * beta
   (isolated double bonds). Significantly negative dE indicates aromatic stabilization.
5. **Compute pi-bond orders**: B_ij = 2 * P_ij from the density matrix P. Aromatic rings have
   bond orders near 1.5; non-aromatic conjugated systems have alternating values.
6. **Classify**: dE per pi-electron exceeds threshold tau -> aromatic.

The HMO approach handles fused systems globally without ring-combination enumeration — the
eigensolve operates on the entire conjugated subgraph. It produces continuous descriptors (bond
orders, delocalization energy) rather than binary yes/no, and handles heteroatoms through
parameterized integrals rather than ad-hoc element lists.

The computational cost (eigensolve of an n x n matrix) is negligible for real molecules. Benzene:
n=6. Coronene: n=24. Porphyrin: n~20. All sub-millisecond. The implementation is no more complex
than the empirical rule thicket it replaces.

HMO bond orders also provide a natural guide for Kekulization: assign double bonds to edges with
highest pi-bond order.

### 3.3 Clar aromaticity (pi-sextet optimization)

Maximizes the number of disjoint pi-sextets across a fused benzenoid ring system. ILP formulation
(Hansen & Zheng 1994): maximize sum of y_r (ring r has a sextet) subject to: for each atom a, sum
of x_b (bond b is a double bond) over bonds incident to a + sum of y_r over rings containing a = 1.
Variables are binary.

Produces differential aromaticity: in phenanthrene, the two outer rings get sextets while the
center ring has a localized double bond. This matches experimental observations (outer rings are
more aromatic). The Hueckel rule cannot express this — it reports all three rings as aromatic.

Limited to benzenoid hydrocarbons (all-carbon, all 6-membered rings). Extension to heterocycles is
an open research problem.

Reference: Liu & Green, Proc. Combust. Inst. 37:575 (2019); Hansen & Zheng, J. Math. Chem. 15:93
(1994).

## 4. Folding structural element perception into the valence registry

### 4.1 The `v` and `a` fields

The `AtomTypeSpec` in the valence registry uses two fields relevant to aromaticity:

- `v` (valence): the bond order sum of **localized** bonds. This excludes aromatic (`a`), dative
  (`>`/`<`), and multicenter (`m`) contributions, and implicit hydrogens (`H`).  `v` counts bond
  orders as they are — a C=C double bond contributes 2, a C-C single bond contributes 1.
  It is not restricted to sigma bonds.
- `a` (aromatic_valence): the number of electrons donated by this atom to an aromatic pi-system.

The total electron budget of an atom is `v + 2*> + a + m + H`, which must
equal the atom's outer-shell electron count adjusted for charge and lone pairs.

Example: benzene carbon in aromatic SMILES `c1ccccc1`. Each carbon has two aromatic ring bonds
(each contributing 1 to the sigma skeleton) and one implicit hydrogen. `v = 2` (2 ring sigma
bonds at order 1 each), `H = 1` (one implicit H), `a = 1` (one pi-electron donated to the aromatic
system). Total 4 = carbon's normal valence.

### 4.2 Electron contributions from candidate sets

The per-atom electron contribution for aromaticity perception is derived directly from the `a` field
on `AtomTypeSpec` candidates, already populated by the valence phase. No separate `ElectronDonorType`
enum is needed.

For each ring atom, the contribution range is `[min(a), max(a)]` across surviving candidates. The
4n+2 check tries assignments within these ranges. Ambiguity is resolved by the check itself:
whichever assignment satisfies 4n+2 survives; candidates with incompatible `a` values are narrowed
out.

This maps cleanly to every RDKit `ElectronDonorType`:

| RDKit ElectronDonorType  | Equivalent in candidate set                      |
| ------------------------ | ------------------------------------------------ |
| `NoElectronDonorType`    | No candidate has `a > 0`                         |
| `VacantElectronDonorType`| Candidate with `a = 0` (empty p-orbital)         |
| `OneElectronDonorType`   | All aromatic candidates have `a = 1`             |
| `TwoElectronDonorType`   | All aromatic candidates have `a = 2`             |
| `OneOrTwoElectronDonorType` | Candidates with both `a = 1` and `a = 2`     |
| `AnyElectronDonorType`   | Dummy atom; handled by element scope exclusion   |

### 4.3 Structural exclusions replaced by registry design

RDKit's `isAtomCandForArom` function applies a series of ad-hoc structural exclusion checks. Each
of these is replaceable by appropriate registry design — if the registry does not define aromatic
candidates (`a > 0`) for a given (element, charge, v) combination, that atom will not be an aromatic
candidate. No separate exclusion flags are needed.

**Multiple unsaturations** (RDKit: atom with more than one double or triple bond is excluded). In our
design: an atom with two double bonds has a high `v` value. If the registry does not include an
`a > 0` candidate at that `v`, the atom is not an aromatic candidate. Example: the central carbon in
allene (C=C=C) has `v = 4` from two double bonds. If `{C+0v4a1}` is not in the registry, this atom
has no aromatic candidate.

**Triple bonds in rings** (RDKit MDL: excluded). Same mechanism: atoms with in-ring triple bonds
have high `v`. If the registry does not include aromatic candidates at that valence, they are
excluded.

**Exocyclic multiple bonds** (RDKit MDL: excluded; RDKit DEFAULT: allowed but electrons are
"stolen"). See section 4.4.

**Radicals on heteroatoms** (RDKit: heteroatoms or charged carbons with radical electrons are
excluded). In our design: candidates with `unpaired_electrons > 0` and `a > 0` would need to be
explicitly present in the registry to allow aromatic radicals. If they are absent, radical atoms have
no aromatic candidate.

**High coordination** (RDKit: degree > 3 excluded). Atoms with degree > 3 have high `v`. If no
aromatic candidate exists at that `v`, the atom is excluded.

### 4.4 Exocyclic bonds and the no-mutation principle

RDKit's DEFAULT model includes an "electronegativity stealing" rule: when a ring atom has an
exocyclic double bond to a more-electronegative atom, one electron is subtracted from the ring
atom's pi-contribution. This is an implicit charge separation performed during aromaticity
perception — a mutation of the electronic bookkeeping.

Our design eliminates this rule. Validation does not mutate the structure, nor does it anticipate
mutation. The input representation determines the result:

**Tropone** (2,4,6-cycloheptatrienone):
- Kekule form `O=C1C=CC=CC=C1`: The C=O carbon has `v = 4` (2 ring sigma + 1 C=O double bond).
  Registry candidate `{C+0v4a0}` gives `a = 0`. The 7-membered ring has 6 remaining pi-electrons
  from the other carbons. 4n+2 with n=1. Aromatic.
- Charge-separated form `[O-][C+]1C=CC=CC=C1`: The C+ has `v = 3`, charge +1. Registry candidate
  `{C+1v3a0}` gives `a = 0`. Same result: 6 pi-electrons, aromatic.

Both representations work without a stealing rule. No mutation needed.

**Fulvene** (methylenecyclopentadiene):
- Neutral form `C=C1C=CC=C1`: The exocyclic-adjacent ring carbon has candidates with both `a = 0`
  and `a = 1`. With `a = 0`: ring has 4 pi-electrons (not 4n+2). With `a = 1`: ring has 5
  pi-electrons (not 4n+2). Not aromatic in either case.
- Charge-separated form `[CH2+][C-]1C=CC=C1`: The ring is cyclopentadienyl anion (Cp-) with 6
  pi-electrons. Aromatic. Whether this interpretation applies is a property of the input
  representation.

This is the correct behavior: aromaticity of a structure as drawn, not aromaticity after implicit
normalization. If a standardization step wants to produce the charge-separated form, it operates
explicitly on an already-validated structure and triggers re-validation.

## 5. Presets and separation of spec from implementation

### 5.1 Presets defined by specification

Presets correspond to chemical specifications, not to toolkit-internal model names:

**Daylight** — the normative aromaticity model for SMILES. Lowercase atom symbols (`c`, `n`, `o`,
`s`) in SMILES mean "aromatic under Daylight rules." Needed for SMILES round-trip correctness.
- `ElementScope`: C, N, O, S, Se, As.
- `RingScope`: no ring size limits; fused combinations checked.
- `RingStrategy`: PiSubgraph (default) or GlobalAllCycles.

**MDL** — the normative aromaticity model for MOL/SDF files. Needed for MOL round-trip correctness.
- `ElementScope`: C, N only.
- `RingScope`: minimum ring size 6; fused combinations up to 6 rings.
- `RingStrategy`: PiSubgraph (default) or GlobalAllCycles.

**Permissive** — broad perception for general use. Defined by "be as inclusive as the chemistry
allows," not by mimicking any specific toolkit.
- `ElementScope`: any element with aromatic valence states in the registry.
- `RingScope`: no ring size limits; fused combinations checked.
- `RingStrategy`: PiSubgraph (default) or GlobalAllCycles.

Force-field-specific models (MMFF94) and toolkit-internal convenience presets (RDKit SIMPLE) are
excluded from the preset list — they do not correspond to specifications and are trivially
constructible from the config struct if ever needed.

### 5.2 Separation of concerns

- **Specification** defines what the aromaticity model *means*: which elements participate, which
  rings are considered, what the pi-electron requirement is (for HueckelRule); or what stabilization
  threshold qualifies as aromatic (for HMO); or what constitutes a benzenoid system (for Clar).
- **Algorithm** implements the check: candidate filtering, ring enumeration and 4n+2 summation
  (HueckelRule); eigensolve, electron filling, and energy/bond-order computation (HMO); ILP
  formulation and solution (Clar).
- **Code** provides the data structures, configuration, and integration into the resolution pipeline.

The specification is captured in the preset constructors and in the atom type registry (which defines
the valid `a` values for each element/charge/valence combination). All three algorithms implement the
`AromaticityModel` trait. The code is organized in `rings.rs`, `aromaticity.rs` (trait +
model-specific submodules), `kekulize.rs`, and modifications to `resolver.rs`, `config.rs`, and
`error.rs`.

## 6. Architecture

### 6.1 AromaticityModel trait

```rust
pub trait AromaticityModel {
    fn detect_aromatic_systems(
        &self,
        builder: &MoleculeBuilder,
        ring_info: &RingInfo,
    ) -> Result<Vec<AromaticSystemCandidate>, ResolutionError>;
}
```

Three implementations:

- `HueckelRuleAromaticity` — 4n+2 electron counting on rings and fused combinations.
- `HmoAromaticity` — Hueckel MO eigensolve on pi-subgraph. Produces bond orders and delocalization
  energy. Classifies aromaticity by stabilization energy threshold.
- `ClarAromaticity` — ILP pi-sextet optimization on benzenoid ring systems.

### 6.2 Pipeline

The aromaticity phase is the third phase of the resolution pipeline, following topology and valence.
It absorbs the existing `kekulize` stub.

1. **Ring detection**: compute biconnected components on the full graph. Then, per `RingStrategy`:
   extract pi-subgraph and enumerate cycles locally (PiSubgraph), or enumerate all cycles per
   biconnected component (GlobalAllCycles).
2. **Element filtering**: restrict candidate rings to those where all atoms pass `ElementScope`.
3. **Aromatic system detection**: for each ring and (when `include_fused` is set) fused combination
   of rings: derive `[min(a), max(a)]` per atom from candidates, enumerate feasible assignments,
   check 4n+2. The combinatorial cost is bounded: typical organic ring atoms have unambiguous `a`,
   so the number of ambiguous atoms per ring is usually 0-2.
4. **Candidate narrowing**: for atoms in detected aromatic systems, eliminate `AtomTypeSpec`
   candidates whose `a` is inconsistent with the detected contribution.
5. **Kekulization**: assign definite single/double bond orders to aromatic bonds such that each
   atom's valence constraints are satisfied. Initial implementation: backtracking DFS. Future:
   augmenting-path matching.
6. **AromaticSystem construction**: build `AromaticSystem` objects from the detected systems and
   attach to `MoleculeBuilder`.

### 6.3 Data flow

```
MoleculeBuilder (from valence phase)
  |
  v
BiconnectedComponents (rings.rs)
  |
  +--[PiSubgraph]--> extract pi-candidates --> connected components --> local enumerate_cycles
  |
  +--[GlobalAllCycles]--> enumerate_cycles per biconnected component
  |
  v
Candidate rings (filtered by ElementScope, RingScope)
  |
  v
AromaticityModel dispatch (HueckelRule / HMO / Clar)
  |
  v
Narrow candidates (remove AtomTypeSpec entries with incompatible 'a')
  |
  v
Kekulize (assign bond orders)
  |
  v
Build AromaticSystem objects, attach to MoleculeBuilder
```

## 7. Configuration

### 7.1 Per-model configuration

**HueckelRuleConfig** — three axes, deliberately minimal:

```rust
pub struct HueckelRuleConfig {
    pub element_scope: ElementScope,
    pub ring_scope: RingScope,
    pub ring_strategy: RingStrategy,
}
```

No `StructuralExclusions` — handled by registry design (section 4.3). No `ExocyclicElectronRule` —
the representation determines the result (section 4.4). No `ElectronDonorType` — derived from
candidate `a` values (section 4.2).

**HmoConfig:**

```rust
pub struct HmoConfig {
    pub element_scope: ElementScope,
    pub stabilization_threshold: f64,  // dE per pi-electron, in units of |beta|
    pub ring_strategy: RingStrategy,
}
```

**ClarConfig:**

```rust
pub struct ClarConfig {
    pub ring_strategy: RingStrategy,
    // Benzenoid-only; no element scope needed (implicitly C only).
}
```

**Shared types:**

```rust
pub enum ElementScope {
    Any,
    AllowList(Vec<Element>),
}

pub struct RingScope {
    pub min_ring_size: usize,
    pub max_ring_size: usize,
    pub include_fused: bool,
    pub max_fused_combination: usize,
}

pub enum RingStrategy {
    PiSubgraph,
    GlobalAllCycles,
}
```

### 7.2 AromaticityResolveConfig

```rust
pub enum AromaticityModelChoice {
    HueckelRule(HueckelRuleConfig),
    Hmo(HmoConfig),
    Clar(ClarConfig),
}

pub struct AromaticityResolveConfig {
    pub enabled: bool,
    pub model: AromaticityModelChoice,
    pub kekulize: bool,
}
```

### 7.3 Error variants

```rust
AromaticityNoKekulization(String),   // no valid Kekule structure exists
AromaticityInconsistent(String),     // aromatic hint but no aromatic system detected
```

## References

- Van-Catledge, "A Pariser-Parr-Pople-based set of Hueckel molecular orbital parameters," J. Org.
  Chem. 45:4801 (1980). doi:10.1021/jo01311a060.
- RDKit aromaticity: `Aromaticity.cpp`, commit a0725e3.
- May & Steinbeck, "Efficient ring perception for the CDK," J. Cheminform. 6:3 (2014).
- Hahnke, Kim, Bolton, "PubChem chemical structure standardization," J. Cheminform. 10:36 (2018).
- Liu & Green, "Capturing aromaticity in automatic mechanism generation software," Proc. Combust.
  Inst. 37:575 (2019).
- Grinberg Dana, Liu, Green, "Automated chemical resonance generation and structure filtration for
  kinetic modeling," Int. J. Chem. Kinet. 51:1 (2019).
- Artemova, Jaillet, Redon, "Automatic molecular structure perception for the Universal Force
  Field," J. Comput. Chem. 37:1191 (2016).
- Berger et al., "Counterexamples in chemical ring perception," J. Chem. Inf. Comput. Sci. 44:323
  (2004).
- Discussion documents: 14-aromaticity-2025-03-08.md, 41-opensmiles-aromaticity-2025-10-01.md,
  53-molecule-validation-scheme-2026-02-17.md, 57-sssr-needed-2026-03-11.md.
