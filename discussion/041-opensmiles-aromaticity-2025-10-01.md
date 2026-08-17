### What RDKit does today (three aromaticity models)

- RDKit default (Daylight-like, in Aromaticity.cpp) [source]
  - Core idea: identify conjugated ring atoms, attempt Kekulization, then mark rings as “aromatic” when the ring can be assigned 4n+2 π-electrons with all ring atoms in an sp2-like state. Bond symbols may remain single in Kekulé form, but atoms are marked aromatic.
  - Mechanics (high level):
    - Find ring membership and a small ring basis.
    - For each candidate ring system, count π-electrons based on atom type/charge rules, check planarity-ish constraints (via atom types), require conjugation continuity.
    - If feasible, mark atoms/bonds aromatic; otherwise revert to non-aromatic Kekulé.
  - Strengths:
    - Fast; works for typical organic rings.
    - Good compatibility with legacy SMILES semantics.
  - Limitations (well known):
    - Fused/polycyclic systems: 4n+2 on a ring-by-ring basis is not globally consistent; results depend on chosen cycle basis.
    - Heteroaromatics and “borderline” cases can flip across models; behavior is model-dependent.
    - Heuristic and not formulated as a global optimization, so ambiguous systems get arbitrary choices.
  - Reference: RDKit’s aromaticity perception implementation in `Aromaticity.cpp` [RDKit source link].

- MDL aromaticity model (CTFile compatibility) [source]
  - Conservative, bond-type oriented view: only certain 5/6-member ring patterns with allowed heteroatom valences are marked aromatic. It largely mirrors historical SDfile semantics rather than electronic structure.
  - Pros/cons:
    - Pros: stable interchange for MDL/SDfile roundtrips.
    - Cons: intentionally limited; not a physics-based criterion; under-calls aromaticity in many fused systems.

- MMFF aromaticity (force-field oriented) [source]
  - Encodes the MMFF94 notion of aromaticity to get the right parameters (charges, torsions). Focused on rings where MMFF has parameters (5/6-member rings with defined patterns).
  - Pros/cons:
    - Pros: matches MMFF parameterization goals.
    - Cons: narrower scope; not appropriate as a canonical cheminformatics definition.

[RDKit source link]: https://github.com/rdkit/rdkit/blob/23ffd85f60d5cbedc86c698933f0fbaeabc81437/Code/GraphMol/Aromaticity.cpp

### Liu–Green (Clar-structure-based) approach

- Core idea: formulate aromaticity in fused polycyclic systems by selecting a maximum set of disjoint π-sextets (Clar sextets), then placing remaining double bonds. This is a global optimization, not ring-by-ring heuristics.
- Algorithmic shape:
  - Build an optimization problem on the cycle graph (typically benzenoid hexagons): binary variables for selecting sextets; constraints ensure disjointness and valence feasibility; objective maximizes number (or weight) of sextets.
  - Yields a Clar cover: which rings are sextetic and where localized double bonds sit.
  - From the Clar cover one can infer fractional bond orders (≈1.5 on sextet edges; ≈1/2 weighting on edges shared by multiple resonance contributors, etc.) and a resonance-energy proxy.
- Strengths:
  - Global, graph-optimization view; stable for large PAHs; captures resonance distribution across fused systems.
  - Produces unambiguous sextet allocations and a principled bond-order picture for benzenoids.
- Limitations:
  - Natively targets benzenoid (hexagonal) systems; generalizing to 5-member hetero rings or general conjugated networks requires extra rules beyond Clar’s original scope.
- Citation: Liu & Green, a Clar-structures-based aromaticity algorithm (ScienceDirect link).

[ScienceDirect link]: https://www.sciencedirect.com/science/article/pii/S1540748918301883

### A graph-based, formal path forward (answering your three questions)

- Q1: Assign fractional bond orders and, where relevant, fractional charge/spin on atoms
  - Option A (preferred general method): Hückel Molecular Orbital (HMO) on the π-graph
    - Build the π-adjacency matrix A on the conjugated subgraph (edges between p-orbital centers), choose simple parameters (α=0, β=−1 scaling), fill electrons by valence rules/charge, compute density matrix P from occupied MOs, and derive bond orders B_ij = 2·P_ij for conjugated edges.
    - Output: fractional bond orders per bond and π-charge/spin per atom. This directly tells you which bonds “demand” non-integer order (e.g., ≈1.5), when the effect is negligible (e.g., butadiene), and can inform implicit-H choices without heuristics.
    - Pros: fully graph-based, global, handles heterocycles, substituent effects, charges, radicals. Clear numeric outputs.
    - Cons: requires small linear algebra (eig/SVD). Still fast for typical SMILES (graphs <10^4 atoms are easy; typical molecules are <<1000).
  - Option B (benzenoids): Clar-ILP (Liu–Green)
    - Use when the conjugated core is essentially hexagonal. Produces sextet assignment and localized double bonds. From this, assign bond orders (sextet edges ~1.5, localized ~2.0 or ~1.0).
    - Pros: crisp for PAHs; matches chemical intuition and modern theory for benzenoids.
    - Cons: needs extension for 5-member rings, heteroatoms, and non-benzenoid graphs.
  - Option C (Kekulé ensemble averaging)
    - Enumerate (or sample) Kekulé structures (perfect matchings on the π-graph) and average bond orders across them; optionally weight by an energy model.
    - Pros: directly tied to resonance structures; naturally delivers fractional bond order.
    - Cons: counting/enumerating perfect matchings is #P-hard in general. Feasible for small/planar cases; can be done with Pfaffians on planar graphs, or via MCMC sampling otherwise.

- Q2: Symmetry-equivalent structures under aromaticity
  - Compute the automorphism group of the conjugated subgraph with constraints that preserve the chosen delocalization descriptor:
    - For Clar: automorphisms that preserve the Clar cover (or any maximum Clar cover if multiple exist) define position orbits. Substitutions within the same orbit are symmetry-equivalent; e.g., 2- and 3- on naphthalene collapse to one orbit under the sextet/double-bond pattern.
    - For HMO: automorphisms that preserve the vector of fractional bond orders (or density matrix invariants) define equivalence classes.
  - Implementation outline: canonical labeling (e.g., partition refinement) using edge/vertex invariants from the chosen aromatic model; compute vertex orbits to report substituent-equivalence classes.

- Q3: A clean, topology-first algorithm (minimize heuristics)
  - Use a two-tier model:
    - Tier 1 (general): HMO on the π-subgraph to produce fractional bond orders and π-charge/spin; accept as the canonical “aromaticity descriptor”. This avoids cycle-basis dependence, handles heteroatoms/charges/radicals, and provides exactly what you need for implicit-H logic.
    - Tier 2 (specialized): Clar-ILP for benzenoids as an optional override when its preconditions hold. This yields crisp sextet assignments and often produces the most interpretable ring-local picture.
  - Verification rather than invention:
    - If SMILES uses lowercase/aromatic tokens, verify they are consistent with the chosen model (e.g., ring system supports delocalization with non-integer bond orders above a threshold).
    - If SMILES uses Kekulé (uppercase, explicit single/double), compute delocalization and optionally suggest normalization (style lint) without rewriting semantics.

### Practical outputs for our pipeline

- From HMO:
  - Per-bond π-bond order b_ij ∈ [0, 2].
  - Per-atom π-charge ρ_i and spin if radicals present.
  - Flags for “significantly non-integer” bonds, e.g., |b_ij − round(b_ij)| ≥ τ (τ ≈ 0.2), to decide when aromatic treatment is semantically essential vs. cosmetic.

- From Clar-ILP (when applicable):
  - Set of sextetic rings, set of localized double bonds.
  - Derived bond-order schematic (sextet edges ~1.5).
  - A consistent Kekulé for export when needed, and symmetry orbits for substitution positions.

- Symmetry/equivalence classes:
  - Orbits of vertices (and edges) under automorphisms that preserve the chosen aromatic descriptor (HMO B_ij or Clar cover). Directly answers when constitutional isomers collapse under aromatic symmetry assumptions.

### Where RDKit’s models fit in

- RDKit default/MDL/MMFF can still be mapped as “profiles” for compatibility:
  - We can run our HMO/Clar verification first and then verify RDKit/MDL/MMFF constraints for reporting compatibility diagnostics or emitting style hints.
  - This preserves interoperability while keeping our core logic principled and global.

### Integration plan (post-parse checker)

- Detection and model selection
  - Detect benzenoid cores (all cycles are 6 and embedded on hexagonal lattice) → enable Clar-ILP fast path.
  - Otherwise, run HMO on the π-subgraph.

- Outputs to STIR/GIR side-channels
  - Store fractional bond orders, π-charges, optional sextet flags.
  - Store symmetry orbits for atoms/positions.

- Diagnostics/lints (examples)
  - AROM_INCONSISTENT_LOWERCASE: lowercase atoms used but no significant delocalization detected.
  - AROM_SUGGEST_LOWERCASE: significant delocalization detected but Kekulé-only SMILES provided (style).
  - AROM_CLAR_CONFLICT: multiple incompatible sextet allocations or sextet fails electron constraints (when Clar profile chosen).
  - AROM_DOUBLE_BOND_IN_AROMATIC_RING: explicit “-” required between two aromatic atoms when model predicts non-aromatic edge (verification).

### Computational notes

- HMO complexity: O(n^3) for eigensolve on the π-subgraph; n is number of conjugated atoms. Typical molecules are small; this is easily sub-millisecond to a few milliseconds.
- Clar-ILP: small 0–1 ILPs over the ring-interference graph; tractable for realistic PAHs. We can also use maximum independent set on the ring graph when constraints match Clar’s classical assumptions.
- Kekulé ensemble: exact counting is expensive; fall back to sampling for larger graphs or planar-Pfaffian tricks for special cases.

### Recommendation

- Make HMO the default verification/annotation engine for all conjugated systems.
- When the molecule is benzenoid, also run Clar-ILP and prefer it for human-facing “aromatic ring” narratives and symmetry classes.
- Keep RDKit/MDL/MMFF as compatibility profiles that translate our numeric picture into their legacy yes/no decisions for interchange and spec-comparison.

- References:
  - RDKit aromaticity implementation: `Aromaticity.cpp` [RDKit source link].
  - Liu & Green Clar-structure algorithm [ScienceDirect link].

- RDKit source link: https://github.com/rdkit/rdkit/blob/23ffd85f60d5cbedc86c698933f0fbaeabc81437/Code/GraphMol/Aromaticity.cpp
- ScienceDirect link: https://www.sciencedirect.com/science/article/pii/S1540748918301883

- In short:
  - Use HMO for general fractional bond orders/π-charges; Clar-ILP for benzenoids; compute symmetry orbits against those descriptors; treat RDKit/MDL/MMFF as compatibility layers rather than the source of truth.
