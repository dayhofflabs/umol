# Prompt

Please review the plan in @/Users/dr/.cursor/plans/aromaticity_perception_framework_5071f98f.plan.md and the design document in @discussion/58-aromaticity-perception-2026-03-11.md for implementing aromaticity perception in umol resolution pipeline. Questions to address:
1. Completeness of design
2. Configurability (integration with @umol-models-graph/src/graph_ir/config.rs configuration) and presets (incl. naming)
3. Lower-level algorithmic prerequisites (graph for computing cycle basis, linear algebra for HMO, ILP for Clar algorithm). Which crates can be usefully employed here besides petgraph for graph algorithms and nalgebra for linear algebra?
4. Configuration of the lower-level algorithms.
5. Data structures and API design for the lower-level algorithms
6. Interaction with the registry design
7. Overall critique of the underlying science
8. Critique design and implementation
9. Testing and comparison with existing implementations (RDKit, CDK, RMG)

# Premium

## Findings

- Critical — current default registry does not encode aromatic electron donation (a > 0) at all, which blocks
  the proposed “registry-driven aromaticity” design. The plan/design assume aromatic candidates are already
  represented in AtomTypeSpec.a, but the default registry currently has only valence/spin/lone-pair style entries.
- Critical — ValenceStrategy::Counts currently forces aromatic_valence = 0, so aromaticity cannot work under
  that strategy without additional inference logic. This directly conflicts with the design’s candidate
  narrowing based on a.
- High — pipeline phase order in current code is incompatible with the proposed architecture.
  kekulize() is called before aromaticity now; the design says aromaticity phase should perform detection,
  narrowing, then kekulization.
- High — planned configurability does not yet fit ResolveConfig shape and defaulting behavior. Current aromatic
  config is only enabled: bool, so model choice/presets/thresholds cannot be expressed until config API is expanded.
- Medium — HMO classification by a single global threshold (tau) is scientifically fragile across heterocycles/charged
  systems. Good as a descriptor, risky as sole yes/no gate without calibration set and fallback logic.
- Medium — Clar scope is correctly limited to benzenoids, but implementation plan should explicitly define
  failure/fallback semantics. If non-benzenoid input is encountered, behavior should be deterministic (error vs
  skip vs fallback model), not implicit.

## Additional crates

- Graph/cycle work:
  - fixedbitset (fast cycle/ring membership bitsets)
  - smallvec (already used)
  - indexmap (stable deterministic iteration)
- Numerics:
  - faer (robust eigensolvers, good performance)
  - optionally ndarray/ndarray-linalg if you want LAPACK-backed decompositions
- ILP (Clar):
  - good_lp as modeling API with backend choice (coin_cbc or highs)
  - or direct backend crate if you want fewer abstractions.
  - optionally sprs if sparse HMO becomes relevant for large fused systems

## Configuration of lower-level algorithms

- Add explicit config structs for algorithm internals:
  - Cycle/rings: max_cycle_len, max_cycles_per_component, time_budget_ms, dedup_mode.
  - HMO: parameter table version (VanCatledge1980 enum), occupation scheme, degeneracy tolerance,
    energy normalization, aromaticity decision mode (Threshold, Ranked, DescriptorOnly).
  - Clar: solver backend, mip gap, time limit, deterministic seed, non-benzenoid policy.

## Data structure design

- Recommend splitting “compute” from “classify”: PiSystem, RingInfo, CycleSet, HmoResult, ClarResult, AromaticAssignment.
- Trait should return rich evidence (energies, bond orders, sextet placement), not just final aromatic systems.
- Keep inputs immutable (&MoleculeBuilder is okay short-term), but introduce graph/candidate views to decouple from resolver internals.

## Interaction with registry design

- Define registry validation checks:
  - disallow contradictory aromatic states,
  - assert expected aromatic-capable element coverage per preset,
  - maintain versioned registry hash for reproducibility.

## Design critique

- Biggest practical risk is introducing high algorithmic complexity before nailing base contracts.
  Suggested implementation order tweak:
  (a) finalize registry + config + ring primitives + HueckelRule end-to-end first,
  (b) then HMO descriptor mode,
  (c) then Clar ILP and optional HMO classification.

# GPT-5.3 Codex

## Findings

- High: cycle enumeration strategy as written is at risk of exponential blow-up.
  “All elementary cycles per biconnected component” is not bounded well enough for polycyclic graphs
  unless you define hard computational budgets and pruning semantics.
- Medium: HMO classification is underspecified for degeneracy/open-shell cases.
  The “20-line eigensolve” claim is directionally attractive but scientifically risky without explicit
  rules for orbital degeneracy, spin state handling, and numerical tolerances.
- Medium: current trait shape forces ring infrastructure even for HMO.
  ring_info is mandatory in the trait API, but HMO is explicitly described as not needing ring enumeration.
  This may force unnecessary coupling/cost.

## Design gaps

- Missing hard acceptance criteria for: performance budgets, numerical stability, fallback behavior,
  and unsupported chemistry behavior (especially Clar).

## Additional crates

- Linear algebra: **nalgebra-lapack** (or **ndarray-linalg**) if you need robust/fast eigensolvers beyond pure-Rust paths.
- ILP/Clar: good_lp as modeling layer; backend via HiGHS or CBC (coin_cbc) depending deployability.

## Data structures and API design

- Keep AromaticityModel but pass a richer context, not always RingInfo.
- Suggested split: AromaticityContext { pi*subgraph, ring_data: Option<*>, atom_candidates }.
- Output should carry both binary and quantitative data:
  - aromatic membership
  - chosen per-atom a
  - optional bond order scores (HMO)
  - provenance (model, parameters, confidence flags)

## Interaction with registry

- Clarify precedence when registry and ElementScope disagree (registry should remain source of truth;
  scope should be a model filter).

## Design critique

- Biggest implementation risk is complexity concentration in ring enumeration/fused-combo handling.
- De-risk by delivering in slices:
  - robust HueckelRule + PiSubgraph first
  - deterministic Kekulization
  - HMO quantitative outputs
  - Clar as opt-in advanced module
  - Add observability from day 1 (debug trace of candidate elimination decisions).

# Gemini 3 Pro

## Algorithmic comment

- ILP (Clar):
  - Avoid heavy ILP solvers (like coin_cbc or minilp) if possible.
  - Mapping to Max Independent Set: The Clar problem (max disjoint sextets) on a benzenoid system is equivalent
    to finding the Maximum Independent Set on the dual graph of the benzenoid rings.
  - Since benzenoid dual graphs are planar and typically small, a simple recursive branch-and-bound algorithm
    implemented in pure Rust will be faster and lighter than linking a generic LP solver.

## Configuration

- Timeout/Limits: GlobalAllCycles needs a hard cutoff (e.g., max_ring_size: 14 or max_cycles: 2000)
  to prevent stalling on highly fused structures like fullerene or infinite lattices if the user provides weird input.

## Data structures and API design

The rings.rs module should expose a clean intermediate representation:

```rust
pub struct RingSystem {
    pub atom_indices: Vec<AtomIndex>,
    pub bonds: Vec<BondIndex>,
    pub rings: Vec<Vec<AtomIndex>>, // Elementary cycles
}

pub struct RingInfo {
    pub systems: Vec<RingSystem>, // Disconnected or biconnected components
}
```

This decouples the detection (topology) from the model (electron counting).

## Science critique

- Exocyclic Bonds: The "no mutation" stance (handling tropone/fulvene via input representation) is theoretically pure
  but puts burden on the input source.
  - Risk: Users providing SMILES O=C1C=CC=CC=C1 might expect it to be aromatic (via "stealing" rules). If umol treats
    it as non-aromatic (because O=C consumes the valence), users might be confused.
  - Mitigation: Detailed error messages or a "Diagnostic" that suggests "Did you mean [O-]...?" if a ring is close
    to aromaticity.

## Algorithm critique

- Kekulization: The plan suggests "Backtracking DFS". For purely aromatic systems, this is fine. For larger conjugated
  systems, it can be slow (exponential).
- Suggestion: Use Edmonds' Blossom algorithm (maximum matching) if performance becomes an issue, but DFS is likely
  fine for the size of typical aromatic systems.

## Test suite additions

- Add:
  - Quinone (non-aromatic)
  - Borazine (inorganic benzene equivalent - tests element scope)
  - Porbsite/Porphyrin (large macrocycles - tests cycle finding limits)

## Implementation adjustments

1. rings.rs: Implement biconnected_components (custom DFS) and enumerate_cycles (Tiernan/Johnson).
2. config.rs: Update AromaticityResolveConfig.
3. aromaticity.rs: Define the trait.
4. HueckelRule: Implement first as it blocks most standard use cases.
5. Integration: Wire into resolver.rs.
6. HmoAromaticity: Implement using nalgebra

# Composer 1.5

## Design gaps

- **System charge**: AromaticSystem (aromatic.rs) has no charge, but discussion-58
  expects total_π = Σ aromatic_valence − charge for systems like tropylium. That
  needs a charge field and electron_count() updated.
- **Resolver order**: In resolver.rs the sequence is kekulize →
  resolve_aromaticity. The plan says aromaticity should absorb kekulize. The order
  should be: ring detection → model dispatch → narrow → kekulize → AromaticSystem
  construction.
- **Default registry**: default-registry.toml has no specs with aromatic_valence
  > 0. Aromaticity models assume such candidates exist. Add at least C, N, O, S
  >    aromatic specs (e.g. {C+0v2a1}, {N+0/1v2a1}, etc.).
- **Aromatic hint precedence**: For aromatic_hint()-marked atoms/bonds, how and when
  are these used vs. model detection? The interaction with candidate narrowing is
  unclear.

## Additional crates

- Check: **biconnected-components** implements Hopcroft–Tarjan; it works with
  UnGraph-like types via traits. Worth checking petgraph compatibility (e.g. Graph, StableGraph).
- **branch-and-bound** (docs.rs): generic B&B with Subproblem; could host the Hansen–Zheng formulation.
- **lpsolve**: external dependency; stronger MILP support but adds a C dependency.
- **Custom LP-only**: Abeledo & Atkinson showed the Clar problem’s LP relaxation
  gives integral solutions (unimodularity). So a pure LP solver (e.g. lp-modeler +
  good_lp, or simplex) could suffice. That simplifies both implementation and
  dependencies.

## Configuration of lower-level algorithms

### Ring enumeration

- max_ring_size in enumerate_cycles — likely 20–24 for practical chemistry.
- RingScope::min_ring_size, max_ring_size — should feed enumerate_cycles and max_fused_combination.
- Consider a max_cycles_per_component cap to avoid blow-up in pathological graphs.

### HMO

- stabilization_threshold (dE per π-electron in |β| units): ~0.1–0.2 β for typical aromatics; make it configurable.
- Van–Catledge parameters: either fixed table or HmoConfig::heteroatom_params: Option<HeteroatomParams> for overrides.
- Electron filling: closed-shell by default; radicals need open-shell handling.

## Data structures and API design

- Ring info

```rust
// rings.rs
pub struct RingInfo { ... }  // or similar aggregate

pub fn biconnected_components<G>(graph: &G) -> Vec<Vec<NodeId>>;
pub fn smallest_ring_size<G>(graph: &G, bond: EdgeIndex) -> Option<usize>;
pub fn enumerate_cycles<G>(subgraph: &G, max_size: usize) -> Vec<Vec<NodeId>>;
pub fn is_cyclic<G>(graph: &G) -> bool;
```

- Aromaticity model

```rust
pub trait AromaticityModel {
    fn detect_aromatic_systems(
        &self,
        builder: &MoleculeBuilder,
        ring_info: &RingInfo,
    ) -> Result<Vec<AromaticSystemCandidate>, ResolutionError>;
}
```

- AromaticSystemCandidate should carry: atoms (and per-atom contributions), charge, and, for HMO,
  bond orders and delocalization energy.

### Method considerations

- **Hückel rule**:
  - Planarity and conjugation: design avoids explicit planarity checks. For macrocycles or twisted rings
    this can mislabel; acceptable for many cheminformatics uses, but worth documenting.
- **HMO**:
  - Threshold τ: inherently empirical; should be tunable.
  - Heteroatoms: parameters exist for common cases (C, N, O, S); others need defaults or a clear extension path.
  - Open-shell: design mentions radicals; HMO should specify how unpaired electrons are placed.
- **Clar**:
  - Multiple Clar structures: ILP gives one optimum; other structures may exist; document that the choice is
    arbitrary unless tie-breaking is defined.
- **Kekulization**:
  - Kekulization: DFS backtracking can be slow for large systems. Plan mentions augmenting-path matching;
    it would be good to scope that for later.

## Error handling

- Add AromaticityModelError(String) for model-specific failures (HMO convergence, ILP infeasibility, etc.).

## Testing

- Add
  Porphyrin — large fused aromatic system.
  Coronene — compact fused benzenoid.
  Pentalene — 4n antiaromatic; should be rejected.
  Biphenylene — 4- and 6-rings; Clar does not apply; HueckelRule edge case.

# Opus 4.6

## Design gaps

- Fused-ring combination enumeration is underspecified. The plan says "fused combinations" are checked for HueckelRule
  but doesn't specify the algorithm. RDKit uses makeRingNeighborMap + pickFusedRings to build fused systems by iteratively
  merging rings that share edges. The design needs to specify: (a) how fused systems are constructed (union of rings
  sharing at least one bond vs. one atom), (b) the traversal/combination strategy, (c) the bound on combination size
  (max_fused_combination in RingScope). The combinatorial explosion is real: phenanthrene has 3 individual rings, 3
  pairwise fusions, and 1 triple fusion = 7 checks. Coronene has 7 rings and 127 subsets. The bound matters.
- `AromaticSystemCandidate` is undefined. The trait returns Vec<AromaticSystemCandidate> but this type is never
  specified. It needs to bridge between model-specific output (HueckelRule: which atoms, what a values; HMO: atom set +
  bond orders + dE; Clar: sextet assignment) and the generic post-processing steps (narrowing, kekulization). This type's
  design is central to the architecture and needs explicit attention.
- Post-detection divergence across models is unaddressed. The pipeline assumes a uniform detect -> narrow -> kekulize ->
  build flow, but the three models produce fundamentally different information:
  - HueckelRule: boolean aromaticity per ring, discrete a assignments.
  - HMO: continuous bond orders, delocalization energy per pi-system. The "narrowing" step is different —
    it's threshold-based on dE, and bond orders guide kekulization rather than being derived from it.
  - Clar: sextet vs. non-sextet rings, localized double bonds already determined by the ILP.

## Required lower-level algorithm implementations

- Biconnected components (Hopcroft-Tarjan): ~80 lines of Rust. No crate needed; the algorithm is simple and
  well-documented. Implementing it directly on petgraph's StableGraph avoids abstraction mismatch.
- Cycle enumeration: The plan says enumerate_cycles(subgraph, max_size) but doesn't specify the algorithm.
  Options:
  - Horton's algorithm for minimum cycle basis — O(V·E²), overkill and gives MCB not all cycles.
  - Johnson's algorithm for all elementary cycles — the standard choice. O((V+E)(C+1)) where C = number of cycles.
    Available via petgraph's ecosystem? No. Must be implemented.
  - Tiernan's algorithm — simpler but less efficient.
  - Relevant cycles (Vismara 1997) — union of all MCBs, unique. More complex to implement but chemically better
    motivated per discussion-57.
- Recommendation: implement Johnson's algorithm bounded by max_ring_size for cycle enumeration. For the PiSubgraph
  strategy, the subgraphs are small enough that the exponential worst case is irrelevant.
- Perfect matching (Kekulization): The pi-subgraph is not necessarily bipartite (azulene, [10]annulene), so
  Hopcroft-Karp doesn't apply. Options:
  - Backtracking DFS (plan's initial choice): correct, simple, adequate for molecules up to ~50 pi-atoms.
  - Edmonds' blossom algorithm for general perfect matching: O(V³). Correct for all cases. The blossom-rs
    crate or similar could be useful.
  - HMO-guided greedy + backtracking: assign double bonds in decreasing order of pi-bond order, backtrack
    on conflict. Nearly linear in practice for molecules with a clear HMO solution.

## Crates to consider

| Crate       | Purpose                                     | Recommendation                                                                    |
| ----------- | ------------------------------------------- | --------------------------------------------------------------------------------- |
| good_lp     | Unified LP/ILP interface (HiGHS, CBC, etc.) | Consider for Clar. Avoids implementing branch-and-bound. But adds FFI dependency. |
| highs       | HiGHS solver bindings                       | Alternative to good_lp. HiGHS is Apache-licensed and handles MIP.                 |
| minilp      | Pure Rust LP solver                         | Insufficient — no integer variables, can't solve Clar ILP directly.               |
| coin_cbc    | CBC solver bindings                         | Alternative. CBC is well-tested for MIP.                                          |
| sprs        | Sparse matrices                             | Not needed. HMO matrices are dense and small (n < 50). nalgebra dense is fine.    |
| pathfinding | Graph search algorithms                     | Not directly useful. Has A\*, Dijkstra, but no cycle enumeration or matching.     |

Recommendation for Clar ILP: For real benzenoid PAHs, the problem is tiny (coronene: 7 ring variables, 24 atom constraints; circumcoronene: ~19 rings).

A hand-written branch-and-bound in pure Rust (~100 lines) is preferable to an FFI dependency for this problem
size. If generality matters, good_lp with highs feature is the cleanest option. Note that for benzenoid systems,
the LP relaxation of the Hansen-Zheng formulation often gives integer solutions (the constraint matrix has
favorable structure), so solving the LP with minilp and checking integrality may suffice in practice —
falling back to branch-and-bound only on fractional solutions.

## Configuration of lowe-level algorithms

- HMO eigensolve has no tunable parameters — it's a fixed computation. The only configuration is the Van-Catledge
  parameter table (h_X, k_CX values). This should be a compile-time constant table, not user-configurable Van-Catledge
  doesn't give a threshold; typical values from Aihara's circuit resonance energy work suggest ~0.01 |beta| per
  pi-electron as a reasonable cutoff, but this is somewhat arbitrary for the binary classification use case. The design
  should document the default and its provenance.
- Kekulization backtracking should have a max_backtrack_steps bound to fail gracefully on pathological inputs rather
  than spinning. This is missing from the design.

## Data structures and API design

- The aromatic system candidates are model-dependent. Suggestions: enum or trait object

```rust
pub struct AromaticSystemCandidate {
    pub atoms: Vec<AtomIndex>,
    pub electron_assignments: Vec<(AtomIndex, u8)>,
    pub rings: Vec<Vec<AtomIndex>>,
    pub metadata: AromaticityMetadata,
}

pub enum AromaticityMetadata {
    HueckelRule,
    Hmo { bond_orders: BTreeMap<(AtomIndex, AtomIndex), f64>, dE: f64 },
    Clar { sextet_rings: Vec<Vec<AtomIndex>> },
}
```

- RingInfo object

```rust
pub struct RingInfo {
    pub biconnected_components: Vec<Vec<AtomIndex>>,
    pub cycles: Vec<Vec<AtomIndex>>,  // all cycles up to max_size
    pub atom_ring_membership: BTreeMap<AtomIndex, Vec<usize>>,  // atom -> cycle indices
    pub bond_ring_membership: BTreeMap<(AtomIndex, AtomIndex), Vec<usize>>,
}
```

- Aromaticity model trait

```rust
pub trait AromaticityModel {
    fn detect_aromatic_systems(
        &self,
        builder: &MoleculeBuilder,
        ring_info: &RingInfo,
    ) -> Result<Vec<AromaticSystemCandidate>, ResolutionError>;
}
```

One concern: HMO doesn't need ring_info for its eigensolve (it operates on the pi-subgraph adjacency). It uses ring_info
only to report which rings are aromatic post-hoc. The trait could take ring_info: Option<&RingInfo> or the HMO
implementation could simply ignore it. The current signature is acceptable — passing unused data is preferable to
complicating the trait.

- AromaticSystem struct extension
  The existing AromaticSystem stores only Vec<AromaticContribution> (atom + a value). The design notes that ring
  membership (rings: Vec<Vec<AtomIndex>>) should be discussed. My recommendation: add it. Clar needs it (to distinguish
  sextet vs. non-sextet rings). HMO benefits from it (bond orders per ring). Even HueckelRule produces it naturally. The
  cost is minimal. Also consider adding electron_count as a cached field rather than recomputing from contributions each
  time. For HMO, per-bond pi-bond orders are a high-value output. Consider a separate PiBondOrderMap attached to
  MoleculeBuilder alongside the AromaticSystem objects, rather than cramming it into AromaticSystem.

## Interaction with Registry Design

The design's strongest architectural decision is deriving electron contributions from registry candidate a values rather
than a separate ElectronDonorType enum. This eliminates an entire class of synchronization bugs between the valence
model and the aromaticity model.

Critical dependency: registry population. The current default registry has zero aromatic entries. The atom type
specifications needed for the core aromatic elements are at minimum:

| SMILES atom           | Registry entry            | Notes                      |
| --------------------- | ------------------------- | -------------------------- |
| c (benzene)           | {Cv3Ha1}, {Cv3a1}         | 1 pi-electron              |
| [nH] (pyrrole N)      | {N/1v2Ha2}                | 2 pi-electrons (lone pair) |
| n (pyridine N)        | {N/1v2a1}                 | 1 pi-electron              |
| o (furan O)           | {O/1v2a2}                 | 2 pi-electrons (lone pair) |
| s (thiophene S)       | {S/1v2a2}                 | 2 pi-electrons (lone pair) |
| [se] (selenophene Se) | {Se/1v2a2}                | 2 pi-electrons             |
| c+ (tropylium C)      | {C+v3a0}                  | 0 pi-electrons (empty p)   |
| c- (Cp- C)            | {C-/1v3a1} or {C-/1v2Ha2} | depends on representation  |

These need careful curation per preset. The ElementScope in HueckelRuleConfig must be consistent with which elements
have a > 0 entries in the registry. A mismatch (element in scope but no aromatic candidates in registry) should be
detected at config validation time, not silently produce empty results.

v semantics for aromatic bonds. The design says v = localized bond order sum. For benzene carbon from aromatic SMILES
(c1ccccc1), each C has two aromatic bonds (each order 1 in BondBuilder) and one implicit H. So v = 2 + 1 = 3 (two ring
bonds + H). But wait — the from_table_bond code sets aromatic bonds to order 1. So v is computed from the order field of
BondBuilder, which is 1 for aromatic bonds. This means v for benzene carbon is 3 (1+1+1 = two ring bonds + one H),
giving candidate {Cv3Ha1}. This is consistent but requires that the valence phase runs before aromaticity and correctly
accounts for aromatic bonds at order 1. The ordering in resolver.rs already does this (valence before aromaticity), so
this is fine.

Narrowing direction. The aromaticity phase narrows candidates by removing those with incompatible a values. This is the
same compositional pattern as the valence phase. But there's a subtlety: if a candidate has a = 1 and the aromaticity
check determines the atom contributes 2 electrons, that candidate is eliminated. The surviving candidates must have a = 2.
But what if no candidate with a = 2 exists? This is the AromaticityInconsistent error — the molecule claims
aromaticity but the registry doesn't support the required electron count. The design handles this correctly in principle
but should specify the error message content clearly.

## Science Critique

- HMO: scientifically sound for the intended purpose. The Hückel approximation (no electron-electron repulsion, only
  nearest-neighbor interactions) is crude by quantum-chemical standards but exactly right for aromaticity classification.
  The key scientific consideration is the Van-Catledge parameterization:
  The parameters are from 1980 and calibrated against experimental ionization energies. They're adequate for C, N, O, S.
  For elements beyond these four (Se, As, Te, P in aromatic rings), Van-Catledge parameters may not be available. The
  design should address this: either extend the parameter table from more recent literature, or fall back to HueckelRule
  for elements without HMO parameters. The delocalization energy threshold (tau) for binary classification introduces a
  discontinuity in what is otherwise a continuous measure. The design correctly identifies this as a cost. For practical
  use, tau should be calibrated against known aromatic/non-aromatic pairs (e.g., benzene vs. cyclobutadiene vs.
  cyclooctatetraene).
- Clar: correct but narrow. The Clar model is scientifically sound for benzenoid PAHs and produces genuinely useful
  differential aromaticity information that the other models cannot. The limitation to all-carbon, all-6-membered ring
  systems is inherent to Clar's theory. The Hansen-Zheng ILP formulation is the standard approach. The scientific risk is
  low.
- Open question: antiaromaticity. The design doesn't address 4n (antiaromatic) systems. Cyclobutadiene and
  cyclooctatetraene are listed as negative tests, which is correct — they should not be classified as aromatic. But should
  the framework explicitly detect and report antiaromaticity? HMO can do this naturally (positive dE = destabilization).
  HueckelRule can check for 4n. Consider adding an Antiaromatic classification to the output, at least optionally.
- Open question: Möbius aromaticity. Twisted annulenes follow a 4n rule (aromatic when 4n pi-electrons). This is
  esoteric but the HMO framework could handle it by adjusting the topology of the Hamiltonian. Not needed for v1 but
  worth noting as a future extension.
- Electron filling for open-shell systems. The HMO design assumes closed-shell aufbau filling. For odd-electron systems
  (e.g., radicals), the filling rule needs clarification: do you fill alpha electrons only? Use fractional occupation? The
  design says electron count comes from candidate a values, which should always give integer total electron counts, but an
  odd number of electrons in a pi-system creates a half-filled orbital. The density matrix calculation changes. This edge
  case needs attention.

## Design critique

- The AromaticityModel trait may need builder mutation access. The current signature takes &MoleculeBuilder (immutable).
  But HMO may want to annotate bond orders on BondBuilders, and candidate narrowing modifies AtomBuilder candidate sets.
  Two options:
  - The trait returns candidates only; the caller (resolver) does the mutation. Cleaner separation.
  - The trait takes &mut MoleculeBuilder. More direct but couples the model to mutation.
- BondBuilder::aromatic_hint lifecycle. The aromatic_hint: Option<bool> is set during topology resolution from input and
  consumed during kekulization. After kekulization, aromatic bonds should have definite orders (1 or 2) and the hint
  becomes irrelevant. The build() method on BondBuilder drops the hint, which is correct. But there's no assertion that
  all Some(true) hints have been resolved before build() is called. Adding a check in MoleculeBuilder::build() that no
  aromatic_hint == Some(true) bonds remain unresolved would catch pipeline ordering bugs.

## Testing

- Quantitative HMO validation against textbook values. Streitwieser's "Molecular Orbital Theory for Organic Chemists"
  (1961) and Yates' HMO tables provide reference eigenvalues, total pi-energies, and bond orders for benzene, naphthalene,
  azulene, etc. These are exact (analytic) for homoatomic systems. The test suite should verify: Benzene: eigenvalues {2,
  1, 1, -1, -1, -2} in units of beta, E_pi = 8 beta, dE = 2 beta, all bond orders = 1.667. Naphthalene: E_pi = 13.683
  beta (10 electrons), dE = 3.683 beta. Azulene: E_pi = 13.364 beta, dE = 3.364 beta, non-uniform bond orders.
- RDKit comparison methodology. For HueckelRule with Daylight preset, the results should match RDKit's AROMATICITY_RDKIT
  model on the common subset of molecules where the electron counting rules agree. Differences will arise from: (a) the
  no-stealing rule for exocyclic bonds, (b) ring set differences (SSSR vs. pi-subgraph). Document these expected
  divergences explicitly in tests with comments explaining why they differ.
- CDK comparison. CDK's ElectronDonation.daylight() + Cycles.all() is the closest analog to Daylight preset +
  GlobalAllCycles. CDK is open-source Java; reference values can be extracted from CDK unit tests (AromaticityTest.java).
- RMG (Reaction Mechanism Generator) comparison. RMG uses Liu-Green's Clar/HMO approach. Their rmgpy.molecule.aromatic
  module provides reference implementations. Particularly useful for Clar and HMO validation on PAHs. RMG's test suite
  (test_aromatic.py) has cases for phenanthrene, pyrene, coronene with expected sextet counts.

# Comments

- A seed set of aromatic atom types for C, N, O, S, Se have been added to the default registry, covering aromaticity in
  benzene, tropylium, cyclopentadienyl, pyrrole, pyridine, pyridinium, furan, pyrylium, thiophene, and thiopyrylium,
  selenophene. This should suffice for algorithm development. Additional elements will be added in a later pass.
- The perception algorithm should not mutate, instead, explicit standardization based on a configurable set of rules
  should be performed prior to validation (separation of concerns).
- Kekulization is **not** part of the resolution process, in contrast to RDKit. GraphIR stores `AromaticSystem` as a
  native property, Kekulization should be implemented as an independent process
  based on a validated aromatic system.
- Input handling: Both aromatic and Kekule structures should be processed by the resolver. Intermediate conversion of
  aromatic systems to Kekule structures with subsequent aromatization (RDKit's sanitization approach) is explicitly
  **not** the model for the implementation here. This means that we should reason through the following cases: Kekule-like
  input, aromaticity perception on or off, aromatic input, aromaticity perception on or off. In either of these four
  cases, narrowing of the options by converting the input to a different representation (e.g., aromatic -> Kekule) is not
  acceptable.
- Hückel parameter sets:
  Review newer papers on the Pariser--Parr--Pople (PPP) method
  https://pubs.aip.org/aip/jcp/article/134/2/024114/965341/Excitation-energy-calculation-of-conjugated
  https://pubs.rsc.org/en/content/articlehtml/2026/dd/d5dd00445d
  for new parameter sets
- Diagnostics is not an error type; it exists for a different purpose (out of scope for this work). Don't touch it.
- Question in need of additional exploration: How should the aromaticity perception interact with the `aromatic_hint`
  fields in `AtomBuilder` and `BondBuilder`? Should the treatment of aromatic hints be configurable?
