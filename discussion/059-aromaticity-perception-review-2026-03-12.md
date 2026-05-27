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

## Premium

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

## GPT-5.3 Codex

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

## Gemini 3 Pro

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

## Composer 1.5

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

## Opus 4.6

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

## Comments

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

## Prompt

Ok, let's review status of the aromatic flags in the inputs to aromaticity perception carefully. We talked about that a bit in @discussion/58-aromaticity-perception-2026-03-11.md but it is not fully worked out, seems like.

1. The TableIR input provides 2 types of aromaticity flags that come from parsing SMILES / MOL inputs:
   a. atom aromatic_hint (obtained from lowercase element names in SMILES (no equivalent in MOL)
   b. bond aromatic_hint (either explicit ":" or implicitly between lowercase atoms in SMILES, bond type "aromatic" in MOL)
2. The atom type valence validator assigns valence types with aromatic valences (`a`) based on the input. At the moment, the matching does not take the aromatic_hints into account. This needs to be built.
3. The counts validator (RDKit-like) only applies a simple valence counting rule based on max valence. There is no handling of aromatic hints in this approach at all. Needs to be built. The RDKit approach is that everything is kekulized, then an aromaticity perception is run. I'd like to get a clearer workflow of their approach. However, I do not want to emulate this (s. below).
   Let me now state the basic semantic principles used throughout umol:
   A. Immutability should be employed wherever possible. That applies both to the key objects (e.g., Atom, Bond, Molecule) should be immutable. But this also applies more broadly to the semantic concepts, meaning that the structures should not be modified by default (no "sanitization" in the RDKit parlance) during parsing (-> TableIR) or resolution (-> GraphIR). This design makes it necessary that the resolution be quite permissive (important to avoid RDKit's infuriating failures to parse fe-porphyrin etc.), which is quite a heavy design constraint. Concretely, this also means that Kekule structures remain in the Kekule formulations -- no aromaticity should be layered onto them in the standard resolution process, not even with aromaticity resolution = on. That is important. The Kekule structure is a valid GraphIR representation. So is an aromatic representation.
   B. Any modifications should be expliciti and operate on a well-defined input after resolution. This means that need to be able to perceive aromatic systems from Kekule structures (fully resolved = GraphiIR Molecule structure). Conversely, we need to be able to kekulize aromatic structures (again from GraphIR representations). This means that the aromaticity perception primitives need to work both on the partially resolved structures (MoleculeBuilder objects produced by the valence validator) using aromatic hints and on fully resolved structures (Molecule objects), generating MoleculeBuilder objects for subsequent resolution.
   We didn't put enough thought into this aspect of the design but it is clearly blocking further progress since this distinction is not well worked out.
   One more note on the transformations, which is not something that needs to be addressed immediately but should be included in the design considerations. Aromaticity perception / Kekulization are two types of transformations that can be applied to GraphIR. Other examples include
4. stereochemistry perception (same principle apply: no sanization, no mandatory perception steps even if stereo resolution is on) / stereo removal (racemization if you will)
5. perception of dative bonds and multicenter bonds (not supported by SMILES and MOL but really important in order to support chemistry properly).
6. "sanitization" steps, for example, fixing of incorrect structures (hypervalent N is an example, CN(=O)=O)
7. ionization / deionization
8. tautomerization
9. reactions
   The preliminary idea is that these transformations do not modify the GraphIR structure but generate a new MoleculeBuilder, which can then be subject to the same resolution phases.
   The two key principles are (they are slightly in tension but that is the nature of the domain):
10. Invalid molecules should not be representable as GraphIR Molecule
11. The definition of what is valid should be quite permissive. Warning is fine, erroring out on hypervalent N or silently modifying it is not.
    Please start by reviewing this outline, the current state of implementation (see @discussion/53-molecule-validation-scheme-2026-02-17.md and @discussion/54-molecule-validation-phases-2026-02-17.md . The aromaticity perception discussion is in @discussion/58-aromaticity-perception-2026-03-11.md and @discussion/59-aromaticity-perception-review-2026-03-12.md ).
    Then make a suggestion for the implementation plan, identify missing pieces, current blockers, restructuring needs in the current umol-models-graph codebase. Ask clariying questions as necessary.

## Composer 1.5

## 2. Current implementation status

### 2.2 Bond aromatic hints — **not used**

The design (discussion-58 §4.5) says: "If a bond has `aromatic_hint = Some(true)` but neither endpoint has aromatic candidates, this is a mismatch." The valence phase does not currently check bond hints. This check belongs in the **aromaticity phase** as a validation step when hints are present, not in valence. Valence only needs to narrow atom candidates; bond-hint consistency can be verified after aromaticity detection (if hinted bonds are not covered by any AromaticSystem → error/warning).

-> This logic is flawed. Since there is no backtracking, if the atom has no aromaticity hint, then aromatic specs will not be picked up by the valence phase.
**TODO**: When setting the aromatic hint on atoms, check incident bonds for aromatic hints. If any bond has an aromatic hint, then the atom is aromatic.

### 2.3 Counts validator — **blocks aromaticity**

```rust:88:99:umol-models-graph/src/graph_ir/resolver/valence.rs
        match AtomTypeSpec::new(
            ...
            AromaticValence::None,
            0,
        ) {
```

The counts validator always uses `AromaticValence::None`. **Both strategies (AtomTyping and Counts) must support aromatic atoms**—this is standard in both MOL and SMILES and cannot be out of scope. Counts needs inference rules: given aromatic hint (atom or bond), element, charge, explicit valence → infer `a` (e.g. C with v=3 and hint → a=1; N with v=2 and hint → a=1 or 2 depending on H). Design and implement as part of Phase 1/5.

## 4. Implementation plan

### Phase 1: Aromatic hint handling in valence (atom typing)

1. **Extend `AtomTypeQuery`** with an explicit aromatic hint constraint:
   - `aromatic_valence_query: Option<AromaticValenceQuery>` where `AromaticValenceQuery` is `RequireAromatic | RequireNonAromatic | Exact(AromaticValence)`.
   - **Inference rule (no mutation):** Use either atom hint or bond hint. `RequireAromatic` when: `atom.aromatic_hint() == Some(true)` OR any incident bond has `aromatic_hint == Some(true)`. `RequireNonAromatic` when: `atom.aromatic_hint() == Some(false)` AND all incident bonds have `aromatic_hint == Some(false)`. Otherwise unconstrained. No changes to atoms; both sources feed the query.

-> **TODO**: Instead mirror the AromaticValence enum but add a new variant: `Any`. Use `Option<AromaticConstraint>` in line with other `Option<T>` fields.

- `Some(AromaticConstraint::Valence(n))` → constrains to candidates with `a = Valence(n)`
- `Some(AromaticConstraint::Any)` → constrains to candidates with `a = Valence(_)`
- `Some(AromaticConstraint::None)` → constrains to candidates with `a = None`
- `None` → no constraint

2. **Update `AtomTypeRegistry::candidates_for`** (or equivalent) to apply the aromatic filter:
   - `RequireAromatic` → retain only specs with `aromatic_valence.is_aromatic()`.
   - `RequireNonAromatic` → retain only specs with `aromatic_valence == None`.
   - `Exact(x)` → retain only specs matching `x`.

-> **TODO**: Update `AtomTypeRegistry::candidates_for` to apply the aromatic filter.

3. **Tests:** Benzene from aromatic SMILES (`c1ccccc1`) should resolve with `AtomTyping` and default registry. Each carbon should get candidates narrowed to aromatic (e.g. `{Cv2a1H}`, `{Cv3a1}` for bridgehead).

### Phase 2: Bond hint validation in aromaticity phase

1. **After** `AromaticityModel::aromatic_systems` returns, for each bond with `aromatic_hint = Some(true)`:
   - Check that both endpoints are in some `AromaticSystem`.
   - If not, apply policy: `AromaticityMatchPolicy::Error` → return `AromaticityInconsistent`; `Ignore` → continue.

2. **Config:** Add `AromaticityResolveConfig.aromatic_hint_policy: AromaticityHintPolicy` with variants `Strict` (error on mismatch), `Lenient` (warn), `Ignore`.

### Phase 3: Transformations (Molecule → MoleculeBuilder)

Resolution uses `resolve_aromaticity_with`, which **mutates** `MoleculeBuilder` in place. Transformations are different: they operate on **immutable** `Molecule` and produce `MoleculeBuilder`.

1. **`kekulize(mol: &Molecule) -> MoleculeBuilder`**
   - Input: `Molecule` with `AromaticSystem` objects.
   - Output: `MoleculeBuilder` with explicit single/double bonds (aromatic bonds replaced by 1 or 2).
   - Module: `umol_models_graph::transform::kekulize`.
   - Algorithm: Backtracking DFS or matching; optionally HMO bond order guidance.
   - Contract: Does not modify input; produces new builder for re-resolution.

2. **`aromatize(mol: &Molecule) -> MoleculeBuilder`** (or similar name)
   - Input: `Molecule` (Kekule form, fully resolved).
   - Output: `MoleculeBuilder` with topology unchanged but `AromaticSystem` objects attached (and bond orders adjusted to aromatic where applicable, per representation).
   - Runs aromaticity perception on the resolved molecule; attaches detected systems to new builder.

-> **TODO**: Need to add a no-op path for resolution if aromatic systems are already present. Similar for valence phase.

### Phase 4: Counts + aromatic (required)

- Design inference rules: given aromatic hint (atom or bond), element, charge, explicit valence → infer `a`.
- Implement for `ValenceStrategy::Counts`. Both strategies must support aromatic SMILES and MOL.

## 5. Restructuring needs in umol-models-graph

### 5.1 Suggested layout

```
graph_ir/
  atom.rs, atom_type.rs, bond.rs, ...
  resolver/
    valence.rs
    aromaticity/
      mod.rs           # AromaticityModel dispatch
      hueckel_rule.rs
      hmo.rs
      clar.rs
  transform/           # NEW: post-resolution transformations
    mod.rs
    aromaticity.rs     # perception on Molecule
    kekulize.rs
    ...
```

-> **TODO**: The aromaticity rules are broader than the resolver and should be in a separate module. One needs to be able to run aromaticity perception on MoleculeBuilder or Molecule
and generate the results according to the semantics of the step (resolution or transformation).

### 5.2 Data flow clarity

| Input type                    | Valence phase                                          | Aromaticity phase                                                                  |
| ----------------------------- | ------------------------------------------------------ | ---------------------------------------------------------------------------------- |
| **Aromatic** (hints)          | Constrain to aromatic candidates via `RequireAromatic` | Validate hinted systems; narrow to single `a` per atom                             |
| **Kekule** (no hints)         | No aromatic constraint                                 | Detect 4n+2; narrow candidates; build AromaticSystem                               |
| **Molecule** (transformation) | N/A                                                    | Perception from resolved `a`; produce AromaticSystem for attachment to new builder |

### 5.3 Config additions

```rust
pub struct AromaticityResolveConfig {
    // enabled removed — see §6.2
    pub model: AromaticityModelChoice,
    pub aromatic_hint_policy: AromaticityHintPolicy,  // Strict | Lenient | Ignore
    // ...
}
```

## Premium

## Proposed Architecture Restructure

- Keep phase pipeline unchanged for resolution, but split aromatic functionality into two reusable layers:
  - Core aromatic primitives that work on a graph view + atom candidate/hint metadata.
  - Transformation adapters:
    - builder adapter for in-pipeline aromaticity phase,
    - molecule adapter for post-resolution transforms that output MoleculeBuilder.
- Introduce an explicit Transformation module family (aromaticity perception, kekulization first) under
  graph_ir, with a stable contract: Molecule -> MoleculeBuilder (or BuilderView -> BuilderDelta internally).

## Implementation Phases

### Phase 1: Enforce strict aromatic hints in valence

- Update atom-typing query construction in [/Users/dr/.cursor/worktrees/umol/uff/umol-models-graph/src/graph_ir/atom_type.rs](/Users/dr/.cursor/worktrees/umol/uff/umol-models-graph/src/graph_ir/atom_type.rs):
  - `aromatic_hint=true` constrains to aromatic-capable candidates (a > 0).
  - `aromatic_hint=false` remains non-aromatic constraint (a=None).
- Add hint mismatch policy wiring for explicit diagnostics (atom-level and bond-level), reusing policy patterns from valence ambiguity/no-match.
- Define bond-hint consistency check: if bond is aromatic-hinted, at least one aromatic-feasible assignment must exist on both endpoints after valence narrowing.

-> **TODO**: Generally agree, a0 does not mean non-aromatic (tropylium). See `AromaticConstraint` discussion under Composer 1.5.

### Phase 2: Make Counts aromatic-capable

- Extend counts_candidates in valence.rs to emit aromatic-capable candidates when constraints permit.
- Add a lightweight aromatic donation mapping for Counts (initially C/N/O/S/Se with charge/valence-conditioned possibilities) so Counts can generate a=1/a=2 alternatives where chemically plausible.

-> **TODO**: yes, please spare me `lightweight` and just do it.

- Ensure implicit-H inference under Counts respects aromatic hints (e.g., fused aromatic carbons with H0 vs edge carbons H1).
- Add ambiguity policy behavior for multi-candidate Counts output (fail/ignore as configured).

### Phase 3: Separate aromatic primitives from resolver internals

- Extract aromatic perception/kekulization core interfaces from aromaticity.rs into reusable modules under a new transformation namespace.

-> **TODO**: Aromaticity primitives should be in a separate module from the resolver.

- Keep resolver aromatic phase as a thin adapter over this core.

-> **TODO**: Yes, but please spare me LLM-babble `thin adapter`.

### Phase 4: Workflow and config consolidation

- Extend resolve/config docs and types so behavior is explicit for:
  - aromatic-hint strictness,
  - mismatch policy,
  - Counts aromatic candidate generation rules.
- Clarify and enforce invariant: no automatic sanitization in parse/resolve; transformations are explicit and produce new builders.
- Align discussion/design docs with actual API contracts in:
  - 53-molecule-validation-scheme-2026-02-17.md
  - 54-molecule-validation-phases-2026-02-17.md
  - 58-aromaticity-perception-2026-03-11.md
  - 59-aromaticity-perception-review-2026-03-12.md

## GPT-5.3 Codex

### Implementation Plan

1. Rework aromaticity config from boolean gate to behavior contract

- Update config.rs
  - Remove/deprecate AromaticityResolveConfig.enabled.
  - Add input_mode (how to treat aromatic hints in resolved output intent).
  - Add hint_policy (strict/warn/ignore).
  - Keep model selection (HueckelRule/Hmo/Clar) + ring strategy; model always available.
- Update resolver flow in resolver.rs so aromaticity phase always executes and chooses behavior by config mode.

-> **TODO**: Fine, let's discuss naming first. `input_mode` is too vague.
What is the difference between `input_mode` and `hint_policy`? There are also chirality hints.

2. Make aromatic hints first-class constraints in valence matching
   -> **TODO**: See `AromaticConstraint` discussion under Composer 1.5.

3. Add bond-hint consistency validation

- Add a small aromatic-hint consistency pass (new helper module under resolver, e.g. resolver/hints.rs):
  For bonds with aromatic_hint == Some(true), require chemically coherent aromatic endpoints/candidates or aromatic-system membership depending on phase.
- Enforce strict/warn/ignore policy.
- Keep hints non-persistent in final Atom/Bond (current behavior remains aligned with your immutability semantics).

-> **TODO**: Fine

4. Define counts-strategy aromatic behavior explicitly

Extend valence.rs with one of two explicit modes:

- CountsAromaticMode::Disabled (counts can parse aromatic input but treats aromatic participation as unresolved; requires later transform to aromatize), or
- CountsAromaticMode::InferFromHints (minimal inference from aromatic hints/ring context for candidate a values).

-> **TODO**: Makes no sense to ignore, same argument as turning off aromaticity perception.

5. Separate “perception” and “transformation” APIs

- Add explicit transformation API surface (new module graph_ir/transform/):
  - perceive_aromaticity_from_kekule(&Molecule, cfg) -> Result<MoleculeBuilder, ...>
  - kekulize_aromatic(&Molecule, cfg) -> Result<MoleculeBuilder, ...>
- Reuse existing ring/model code rings.rs, resolver/aromaticity/\*.rs) but avoid in-place mutation of existing Molecule.

-> **TODO**: Fine, shorter names, `write_a_tolstoy_novel_just_to_check_the_compiler_limits_and_make_myself_look_very_sophisticated`
is not my style.

6. Upgrade kekulization output to builder-producing transform

- Keep existing solver core in kekule.rs.
- Add wrapper that applies assignment onto a cloned builder, clears aromatic hints/systems as configured, and returns `MoleculeBuilder` for re-resolution.

-> **TODO**: Needs deeper analysis that this shooting from the hip.

7. Add lifecycle invariants and tests for the four input/representation cases

Add tests under tests/resolution/ and aromaticity unit suites:

- Kekule input -> aromatic system perceived (no bond-order mutation in resolver path).
- Aromatic input -> valence constrained by hints; mismatch policy behavior verified.
- Aromatic representation retained as aromatic GraphIR when configured.
- Explicit transform round-trips: aromatic -> kekule transform; kekule -> aromatic transform.
- Add regression cases: tropylium, Cp-, fulvene/tropone representations, porphyrin-like permissive handling.

-> **TODO**: Fine

## Claude Opus 4.6

## Analysis: Aromatic Flags in Resolution

### Current State Summary

The pipeline currently flows as:

SMILES/MOL → TableIR (aromatic_hint on atoms/bonds) → Topology → Valence → Aromaticity → Stereo → Build → Molecule

Here's what works and what doesn't at each stage.

1. TableIR → AtomBuilder/BondBuilder (works correctly)

Aromatic hints propagate cleanly:

SMILES c → TableAtom.aromatic = Some(true) → AtomBuilder.aromatic_hint = Some(true)
Between two aromatic atoms: BondOrder::Aromatic → BondBuilder(order=1, aromatic_hint=Some(true))
MOL aromatic bond type → same path 2. Valence matching (partially broken)

2. Valence matching (partially broken)

```rust
    let aromatic_valence = match atom.aromatic_hint() {
    Some(false) => Some(AromaticValence::None),
    - => None,
      };
```

This is asymmetric: aromatic_hint = Some(false) correctly constrains to non-aromatic specs, but aromatic_hint = Some(true) leaves aromatic_valence = None in the query, meaning it matches both aromatic and non-aromatic specs. For benzene C from c1ccccc1, the query ?{Cv2} matches both {Cv2a1H} (aromatic, 1 pi-electron, 1 implicit H) and potentially {Cv2H2} (non-aromatic, 2 implicit H) if that spec exists. This causes ambiguity rather than correct aromatic-constrained matching.

3. Counts validator (broken for aromatic input)

```rust
AromaticValence::None,
0,
```

Hardcodes AromaticValence::None. Aromatic atoms cannot be represented at all under Counts strategy.

4. Aromaticity phase (incomplete)

The phase detects aromatic systems and adds them to the builder, but does not narrow atom candidates. After valence + aromaticity phases, atoms with aromatic_hint = Some(true) can still have multiple candidates (both aromatic and non-aromatic), causing build() to fail on ambiguity.

5. Build phase (missing validation)

BondBuilder::build() silently drops aromatic_hint without asserting that all aromatic_hint = Some(true) bonds were resolved into AromaticSystem membership.

6. Registry (populated correctly)

The default registry contains the necessary aromatic entries: {Cv2a1H}, {Cv3a1}, {N/1v2a1}, {Nv2a2H}, {O/1v2a2}, {S/1v2a2}, {Se/1v2a2}, etc.

## Key Architectural Clarification: Resolution vs. Transformation

Your message establishes a critical principle that partially supersedes the input handling matrix in discussion/58:

> Kekule structures remain in the Kekule formulations — no aromaticity should be layered onto them in the standard resolution process,
> not even with aromaticity resolution = on.

| Input                             | Resolution behavior                                                                                                                                                      |
| --------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Kekule (no aromatic hints)        | Valence matches non-aromatic specs. No AromaticSystem created. Molecule is Kekule. Always.                                                                               |
| Aromatic (aromatic hints present) | Valence constrained to aromatic specs. Aromaticity phase validates ring membership + Hückel (or configured model). AromaticSystem objects created. Molecule is aromatic. |
| Mixed (e.g. c1ccccc1CC)           | Aromatic atoms → aromatic specs; aliphatic atoms → non-aromatic specs. AromaticSystem for the aromatic part only.                                                        |

And separately, transformations operate on fully resolved Molecule objects:

```rust
Molecule (Kekule) → perceive_aromaticity() → MoleculeBuilder → resolve → Molecule (aromatic)
Molecule (aromatic) → kekulize() → MoleculeBuilder → resolve → Molecule (Kekule)
```

This means the aromaticity resolution phase during the standard pipeline only needs to handle aromatic-hinted input. It does not detect aromaticity in Kekule structures — that is a separate, explicit transformation.

## The Four Concrete Gaps to Close

### Gap 1: Aromatic hint constraining in valence matching (AtomType strategy)

`aromatic_hint = Some(true)` must constrain the query to `aromatic_valence = Some(AromaticValence::Valence(*))` — meaning only aromatic specs match. This is the single most important fix. Currently lines 421–424 of atom_type.rs need to become:

```rust
let aromatic_valence = match atom.aromatic_hint() {
Some(true) => None, // match any Valence(*), exclude None
Some(false) => Some(AromaticValence::None),
None => None,
};
```

But this alone isn't sufficient — None in the query means "match anything", which includes non-aromatic specs too. What we actually need is a richer query constraint. The query's `aromatic_valence: Option<AromaticValence>` field can't express "must be `Valence(*)` for any value of the inner u8." This requires either:

- A new query variant (e.g. `AromaticValenceConstraint::MustBeAromatic` | `MustBeNonAromatic` | `Any`)
- Or changing the matching logic to filter differently

-> **TODO**: See discussion regarding `AromaticConstraint` in Composer 1.5 addressing this issue.

### Gap 2: Aromatic hint handling in Counts strategy

The Counts validator needs to be extended. When aromatic_hint = Some(true):

- The atom is in an aromatic ring, contributing some pi-electrons
- Bond order sum from graph edges gives the σ-skeleton (aromatic bonds at order 1)
- Need to infer aromatic valence from element + charge + σ-skeleton
- Implicit H inference must account for the pi-electron budget

This is non-trivial. The RDKit approach is to kekulize first (assign alternating single/double bonds), then run aromaticity perception on the kekulized form. Since we explicitly reject that approach, the Counts strategy needs its own aromatic inference logic, or it should be clearly documented as not supporting aromatic input (requiring the AtomType strategy for aromatic SMILES).

-> **TODO**: Cannot reject this strategy because this is too hard. Boo-hoo. The kekulization is a lazy shortcut, not a scientific necessity. The fact that there is usually only one spec in the registry per atom / charge state (which are known) tells us that this is solvable.
What we need to add are calculations for:

- Number of unpaired electrons (need the equivalent to "radical assignment" in RDKit), a bit rough but acceptable. This is the essence of the counts method -- a bit imprecise but does not need a registry (see the "lightweight aromatic donation mapping" in Premium).
- Number of implicit hydrogens (naphthalene bridging = H0 vs alpha/beta carbons = H1)? This is analogous to the aliphatic case once the aromatic valence is established. Needs to be done but there is plenty of precedent to fall back on.

### Gap 3: Aromaticity phase candidate narrowing

Currently resolve_aromaticity_with adds AromaticSystem objects but doesn't narrow atom candidates. It needs to:

1. For each atom with aromatic_hint = Some(true): verify that the atom is in a detected AromaticSystem
2. Narrow candidates to those whose aromatic_valence is consistent with the AromaticSystem's pi-electron assignment
3. Handle the error case: aromatic-hinted atom not in any detected system → AromaticityInconsistent

### Gap 4: Bond aromatic hint validation at build time

MoleculeBuilder::build() should verify that every bond with aromatic_hint = Some(true) has both endpoints in some AromaticSystem. If not, this is a resolution pipeline bug (aromatic bond not resolved).
