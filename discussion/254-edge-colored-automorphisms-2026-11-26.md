# 254 — Edge-colored graph automorphisms

Status: In Progress
Date: 2026-11-26
Relates: [120](120-molecular-graph-generation-2026-06-21.md),
[143](143-vendored-nauty-integration-2026-07-11.md),
[176](176-ast-naming-2026-07-31.md),
[216](216-canonicalization-performance-2026-08-30.md),
[233](233-graph-ir-delta-review-2026-09-23.md),
[238](238-non-cascading-removal-2026-10-16.md),
[245](245-canonical-eq-cost-2026-11-12.md),
[246](246-stereo-derived-readings-2026-11-12.md),
[255](255-reaction-span-bond-normalization-2026-11-26.md),
[data types guide](../docs/development/data-types.md),
[integrity guide](../docs/development/integrity.md),
[nomenclature guide](../docs/development/nomenclature.md),
[operations guide](../docs/development/operations.md),
[python api guide](../docs/development/python-api.md)

## Purpose

`umol-graph-core` currently delegates graph automorphism and canonical-labeling computation to
vendored sparse nauty (`umol-nauty-sys`). Because nauty accepts only vertex colors, every caller
that distinguishes edge types (`Graph::canonical_key`, `umol-graph-ir` topology/constitution/
structure/full canonicalization, and `GraphSymmetry`) must subdivide colored edges into extra
vertices.

Exploratory prototyping and benchmarking in `scratch/automorphism-edges` established that:
1. Native edge colors plus a pure-Rust trail-backed individualization--refinement engine
   (`AutomorphismAlgorithm::TrailIr`) beat subdivided nauty across the full 39-case benchmark suite
   (33 Tier 1--4 cases + 6 Tier-5 multi-ligand chemical complexes) and beat weighted C nauty
   (`dispatch_swg`) on 38/39 cases, delivering 1.1×--2.6× end-to-end speedups on symmetric
   `Molecule::canonicalize`, up to 1.19× on `ReactionSpan::canonicalize`, and 1.48×--29.3× on
   heavy stereo-symmetry derived readings.
2. Both `Nauty` and `TrailIr` require simple graphs. Eliminating the two structural non-simple
   `Graph` representations in `umol-graph-ir` (`IncidenceGraph` virtual stereo-ligand backlinks
   and `ReactionSpan` removed+added bonds on one atom pair), enforcing localized-bond simplicity
   during editing, and deduplicating two higher-crate edge-list builders (`umol-graph` HMO
   ring-bond unions and `umol-io` SMILES stereo-marker components) allow `Graph` and `DiGraph`
   themselves to become strictly simple.

This document records the settled design, data-type contracts, exploratory evidence, and the
staged production implementation plan on `feature/automorphism-edge-colors`.

### Priorities and objectives

- **Primary objective — pure performance.** Keeping `umol-nauty-sys` (extended with weighted
  `dispatch_swg`) as an explicit secondary backend (`AutomorphismAlgorithm::Nauty`) preserves an
  independent C differential oracle and fallback for extreme permutation-group combinatorial graphs,
  while `AutomorphismAlgorithm::TrailIr` becomes the primary pure-Rust backend and higher-layer
  default.
- **Primary target workloads:**
  1. Complete canonicalization of `Molecule`, `Reaction`, and `ReactionSpan` across `Topology`,
     `Constitution`, `Structure`, and `Full` prefixes (including reaction-network construction).
  2. Automorphism and symmetry queries supporting stereochemistry (`GraphSymmetry` and derived
     readings in doc [246](246-stereo-derived-readings-2026-11-12.md)).
- **Secondary workloads:** `Graph::canonical_key` (`SubstructureFeaturizer`) and `GraphView`
  atom/bond automorphism queries.

## Settled design and data-type contracts

### 1. Strictly simple `Graph` and `DiGraph` (prerequisite cleanup)

`umol_graph_core::Graph` and `umol_graph_core::DiGraph` become strictly simple: no out-of-range
endpoints, no self-loops, and no parallel edges.

```text
Type and role: GraphError (error type for Graph and DiGraph construction and edge mutation)
Open carrier or operation-issued value: error enum
Intrinsic representation invariants: N/A
Variants:
  - NodeOutOfBounds { node: u32, node_count: usize }
  - SelfLoop { node: u32 }
  - ParallelEdge { endpoints: [u32; 2] }
Public constructors: Graph::try_new / Graph::new, Graph::try_add_edge / Graph::add_edge,
  Graph::try_extend_edges / Graph::extend_edges, Graph::try_extend / Graph::extend,
  DiGraph::try_new / DiGraph::new
Failure, absence, and panic behavior: checked try_* methods return Result<_, GraphError> without
  mutating the receiver on error; asserted counterparts check the same contract and panic on
  violation. Graph::restore retains its existing matching-history contract.
Rust/Python boundary: Rust-only (Graph and DiGraph are not bound in umol-py).
```

Retired remnants and multigraph fallbacks removed in the same cleanup:
- `Graph::{is_simple, is_dense, node_bound, edge_bound, subdivide_edges}` (all `node_bound` /
  `edge_bound` internal call sites switch to `node_count` / `edge_count`).
- `SubdividedGraph`, `SubdivisionNodeSource`, `NonSimpleGraphError`.
- The 8 multigraph cycle fallback / restricted-dispatch methods in `algorithms/cycles.rs`
  (`try_visit_simple_cycles`, `visit_simple_cycles_fallback`, `try_enumerate_simple_cycles`,
  `enumerate_simple_cycles_fallback`, `try_visit_relevant_cycles`, `visit_relevant_cycles_fallback`,
  `try_enumerate_relevant_cycles`, `enumerate_relevant_cycles_fallback`), keeping only the 4 direct
  methods (`visit_simple_cycles`, `enumerate_simple_cycles`, `visit_relevant_cycles`,
  `enumerate_relevant_cycles`) and dropping loop/subdivision paths in Horton MCB and Vismara.
- Parallel-edge ambiguity checks in `Correspondence::edge_matched_pairs` / `shared_edge_count`
  (`umol-graph-core/src/correspondence.rs`).
- Because `Molecule` stores `graph: Graph` directly and `Transaction` (`RollbackMethod::Journal`)
  mutates `&mut Molecule` in place, transient self-loop and parallel bonds are no longer
  representable during editing: `Edit::AddBonds` rejects self-loops and parallel bonds before
  mutating `Graph` with `MoleculeApplyError::Integrity(DuplicateAtom / ParallelBonds)` (leaving
  `Reaction::apply_deltas_at`'s `Ok(None)` path and Python's `InvalidStructureError` intact without
  adding duplicate `TransactionError` variants), and `MoleculeEditor::{add_bond, add_bonds,
  add_topology}` remain asserted (panicking on self-loops and parallel bonds, matching out-of-range
  atoms). `Molecule::try_from_entries` and `ReactionSpan::try_from_entries` retain their exact error
  contracts.

### 2. Simple `IncidenceGraph` (`StereoAnchor`) and `ReactionSpan` bond uniqueness

- **`IncidenceGraph` (`umol-graph-ir/src/ir/incidence.rs`):** eliminates parallel edges uniformly
  for both `StereoAtom` and `StereoBond` with zero extra vertices (`node_count == sum(entity_counts)`
  and total `entity(NodeId) -> Entity` preserved):
  - `StereoAtom`: one `Incidence::StereoAnchor(StereoVirtualLigands)` edge to `site_atom` + one
    `Incidence::StereoLigand` edge per actual-atom ligand.
  - `StereoBond`: one `Incidence::StereoSite` edge to `site_bond` (connecting to the `Bond` entity
    node) + one `Incidence::StereoAnchor(StereoVirtualLigands)` edge to each endpoint atom `a`, `b`
    + one `Incidence::StereoLigand` edge per actual-atom ligand.
  - `StereoVirtualLigands`: `bitflags!` `u8` carrier (`IMPLICIT_HYDROGEN`, `LONE_PAIR`, `VIRTUAL`)
    deriving `Copy + Clone + Debug + Default + PartialEq + Eq + PartialOrd + Ord + Hash`. Because
    representation integrity forbids anchor atoms from appearing as actual ligands of the same
    stereo entity, every `(stereo_entity, anchor_atom)` edge is unique and `IncidenceGraph::graph()`
    is always simple.
- **`ReactionSpan` (`umol-graph-ir/src/ir/reaction_span.rs`):** because localized bond identity is
  its unordered endpoint pair `{a, b}` ([integrity.md](../docs/development/integrity.md)),
  `ReactionSpanIntegrityError` gains `ParallelBonds { atoms: [AtomId; 2] }`, rejected in
  `check_reaction_span_entries` and `ReactionSpan::check_integrity` (and inherited by
  `Reaction::to_reaction_span`). Settling `Reaction` / `ReactionSpan` delta normalization of
  removed+added bonds on the same endpoint pair (and revisiting this check) is tracked separately in
  [255](255-reaction-span-bond-normalization-2026-11-26.md) alongside
  [233](233-graph-ir-delta-review-2026-09-23.md).

### 3. Edge-colored automorphism and canonical-labeling API (`umol-graph-core`)

#### Algorithm selector

```rust
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AutomorphismAlgorithm {
    Nauty,
    TrailIr,
}
```

- **`TrailIr`:** pure-Rust trail-backed edge-colored IR engine (zero-copy `Graph` CSR + `&[u32]`
  edge colors, `uniform_edges` fast path, undo-trail `Partition`, hybrid $O(h)$ sparse / 2-color
  zero-sort Hopcroft refinement over the full `u32` edge-color range, exact false-twin and
  monochromatic true-twin collapse, O(1)-gated recursive NUCR at `nnt >= 8`, `FirstSmallest`
  target-cell heuristic, flat DFS stacks, and `get_canon: bool`). The prefixed `Ir` abbreviation is
  accepted as a concise, unambiguous discriminator in this context (following the discriminator
  rationale in [176](176-ast-naming-2026-07-31.md)). No separate `Bliss` variant is added because
  `TrailIr` with gated recursive NUCR dominates both NUCR-free `TrailIr` on large multi-ligand
  complexes and twin-free Bliss on H-/methyl-rich chemistry.
- **`Nauty`:** vendored sparse nauty extended with native weighted `dispatch_swg` (`umol-nauty-sys`,
  supporting both `get_canon = true` and `get_canon = false`). Retained as the independent C
  differential oracle and fallback for extreme permutation-group combinatorial graphs.

#### Operations on `Graph`

All operations are infallible on strictly simple `Graph`:

```rust
impl Graph {
    pub fn automorphisms<NC: Ord, EC: Ord>(
        &self,
        node_color: impl Fn(NodeId) -> NC,
        edge_color: impl Fn(EdgeId) -> EC,
        alg: AutomorphismAlgorithm,
    ) -> AutomorphismGroup;

    pub fn automorphisms_with_scratch<NC: Ord, EC: Ord>(
        &self,
        node_color: impl Fn(NodeId) -> NC,
        edge_color: impl Fn(EdgeId) -> EC,
        alg: AutomorphismAlgorithm,
        scratch: &mut AutomorphismScratch,
    ) -> AutomorphismGroup;

    pub fn canonical_labeling<NC: Ord, EC: Ord>(
        &self,
        node_color: impl Fn(NodeId) -> NC,
        edge_color: impl Fn(EdgeId) -> EC,
        alg: AutomorphismAlgorithm,
    ) -> CanonicalLabeling;

    pub fn canonical_labeling_with_scratch<NC: Ord, EC: Ord>(
        &self,
        node_color: impl Fn(NodeId) -> NC,
        edge_color: impl Fn(EdgeId) -> EC,
        alg: AutomorphismAlgorithm,
        scratch: &mut AutomorphismScratch,
    ) -> CanonicalLabeling;

    pub fn canonical_key(
        &self,
        node_color: impl Fn(NodeId) -> Vec<u8>,
        edge_color: impl Fn(EdgeId) -> Vec<u8>,
        alg: AutomorphismAlgorithm,
    ) -> Vec<u8>;
}
```

- Dropping `+ Copy` (requiring only `Ord`) lets both `Copy` ranks and rich `Ord` color types
  (`Vec<u8>`, composite tuples) be passed directly without manual pre-ranking.
- `Graph::canonical_key` calls `self.canonical_labeling` directly (zero edge subdivision) and
  encodes canonical node colors followed by the sorted canonical `(pos[u], pos[v], edge_color)` list.

#### Data-type contract sheets

```text
Type and role: AutomorphismGroup (symmetry result of Graph::automorphisms / automorphisms_with_scratch)
Open carrier or operation-issued value: closed, operation-issued value (replaces AutomorphismOutput)
Intrinsic representation invariants:
  - orbits.len() == node_count; each orbit representative is the minimum NodeId in its orbit.
  - orbit_count equals the number of distinct values in orbits.
  - each generator has length node_count and is a bijection over 0..node_count preserving graph
    adjacency, node colors, and edge colors.
  - the subgroup generated by generators has order group_order and orbit partition orbits.
Contextual properties and supplied context: NodeId queries (orbit_of, same_orbit) assume node ids
  belong to the graph that produced this value.
Public constructors: none (issued only by Graph::automorphisms* and CanonicalLabeling).
Public accessors: node_count(&self) -> usize, orbit_count(&self) -> usize,
  orbit_of(&self, node: NodeId) -> NodeId, same_orbit(&self, first: NodeId, second: NodeId) -> bool,
  group_order(&self) -> AutomorphismGroupOrder, generators(&self) -> &[Vec<NodeId>].
Failure, absence, and panic behavior: orbit_of / same_orbit panic on out-of-bounds NodeId (matching
  Graph::degree and Graph::neighbors).
Algebraic / semantic properties:
  - Orbits, orbit_count, group_order, and the generated subgroup are invariant across
    AutomorphismAlgorithm variants (Nauty and TrailIr agree).
  - The specific generating basis in generators() is algorithm-dependent.
  - automorphisms and automorphisms_with_scratch produce identical results for every input.
Rust/Python boundary: Rust-only (umol-py exposes only AutomorphismAlgorithm and higher-layer configs).
```

```text
Type and role: CanonicalLabeling (canonical node order plus AutomorphismGroup)
Open carrier or operation-issued value: closed, operation-issued value
Intrinsic representation invariants:
  - labels.len() == automorphisms.node_count() and labels is a bijection over 0..node_count.
  - automorphisms satisfies all AutomorphismGroup invariants.
Public constructors: none (issued only by Graph::canonical_labeling*).
Public accessors: labels(&self) -> &[NodeId], automorphisms(&self) -> &AutomorphismGroup,
  into_automorphisms(self) -> AutomorphismGroup.
Algebraic / semantic properties:
  - For a fixed AutomorphismAlgorithm and color ordering, labels() induces an identical canonical
    colored adjacency matrix across all vertex renumberings of isomorphic graphs.
  - Tie-break ordering of labels() is algorithm-dependent (Nauty vs TrailIr need not emit identical
    permutations, though both induce valid canonical keys).
  - canonical_labeling and canonical_labeling_with_scratch produce identical results.
Rust/Python boundary: Rust-only.
```

```text
Type and role: AutomorphismScratch (reusable opaque scratch storage for automorphism/canonical calls)
Open carrier or operation-issued value: opaque reusable scratch carrier
Intrinsic representation invariants: no caller-observable state; safe to reuse across arbitrary
  graphs, sizes, colorings, and AutomorphismAlgorithm variants.
Public constructors: AutomorphismScratch::new() and Default (plus Clone, Debug).
Failure, absence, and panic behavior: infallible.
Algebraic / semantic properties: prior calls on a scratch instance never affect the output of
  subsequent calls; *_with_scratch produces results identical to the corresponding no-scratch call.
Rust/Python boundary: Rust-only (not bound in umol-py or GraphView).
```

#### Repo-wide `*Scratch` / `*_with_scratch` admission policy (to record in `nomenclature.md`)

- Add a `*Scratch` type and `*_with_scratch` method counterpart only when benchmarks establish that
  per-call buffer allocation is a measurable share of runtime in repeated-call loops (search trees,
  stabilizers, batch calls).
- Always keep the no-scratch operation as the primary entry point with the same receiver (`Graph`).
- Keep `*Scratch` opaque, problem-scoped (algorithm-independent), reusable across arbitrary inputs
  with output identical to the no-scratch method, and Rust-only unless binding benchmarks justify
  otherwise.

### 4. Downstream alignment (`umol-graph-ir`, `umol-graph`, `umol-py`, docs)

- **`GraphView` (`umol-graph-ir/src/ir/view/graph.rs`):**
  - `GraphView::automorphisms<AC: Ord, BC: Ord>(&self, atom_color: impl Fn(AtomId) -> AC, bond_color: impl Fn(BondId) -> BC, alg: AutomorphismAlgorithm) -> AtomAutomorphismGroup` (replaces `AtomAutomorphism`, exposing `atom_count`, `orbit_count`, `orbit_of`, `same_orbit`, `group_order`, `generators`).
  - `GraphView::canonical_labeling<AC: Ord, BC: Ord>(&self, atom_color: impl Fn(AtomId) -> AC, bond_color: impl Fn(BondId) -> BC, alg: AutomorphismAlgorithm) -> AtomCanonicalLabeling` (exposing `labels(&self) -> Vec<AtomId>` and `automorphisms(&self) -> &AtomAutomorphismGroup`).
- **Internal `umol-graph-ir` canonicalization and `GraphSymmetry`:**
  - `AutomorphismAdapter::new` wraps `incidence_graph.graph().clone()` with `colors` and
    `edge_colors = initial_colors.incidences.clone()` (zero subdivision, zero `Graph::new`).
  - `topology_search_input` (both `Molecule` and `ReactionSpan`) uses `molecule.raw_graph().clone()`
    / `span.graph().clone()` directly with `edge_colors` (zero modal bond subdivision).
  - `OrderedPartition::refine` uses sparse run-length-encoded `(cell_col, edge_color, count)`
    signatures in exact descending dense-count-vector order ($O(\sum \deg)$ per non-singleton cell).
  - `canonical_search_impl` reuses one stack-local `AutomorphismScratch` across `search_partition`
    (`canonical_labeling_with_scratch`); `visit_remapping_candidates` reuses one stack-local
    `AutomorphismScratch` (`automorphisms_with_scratch`); `Molecule::graph_symmetry` calls
    `Graph::automorphisms` (`get_canon = false`) and stores `Option<Box<AutomorphismGroup>>`.
  - Unit tests assert algorithm-independent invariants (canonicalization parity between `Nauty` and
    `TrailIr`, relabeling invariance, sound pruning bounds) rather than hard-coding backend-specific
    generating bases or tie-broken labelings.
- **Higher-layer config defaults (Rust + Python + `operations.md`):**
  - Switch default `automorphism_algorithm` from `AutomorphismAlgorithm::Nauty` to
    **`AutomorphismAlgorithm::TrailIr`** in `GraphSymmetryConfig` (`umol-graph-ir`),
    `CanonicalizeConfig`, `AddStereoConfig`, `RemoveStereoConfig`, `EnumerateStereoConfig`,
    `StereoValidateConfig`, `StructuralFingerprintConfig` (`umol-graph`), and all corresponding
    `umol-py` configs + `umol_py::algorithm::AutomorphismAlgorithm::{Nauty(), TrailIr()}` +
    `umol-py/python/umol/_umol.pyi` + [operations.md](../docs/development/operations.md).

## Exploratory evidence (`scratch/automorphism-edges`, 2026-11-26)

All experiments below ran in the isolated `scratch/automorphism-edges` worktree on `rustc 1.99.1
(58cf8f9c3 2026-10-13)`, Apple M1 Pro `arm64`, Criterion `--warm-up-time 0.2 --measurement-time 0.5
--sample-size 20` (medians in µs unless noted).

### 10--20 atom FFI and setup breakdown (unsubdivided topology floor)

`rust_prep`: color sort + CSR build + `NautyInput::try_new`; `c_unique`: `umol_nauty_sys::run` on
all-unique colors (C validation + 9 `malloc`/`free` calls + 1-pass `sparsenauty` floor);
`c_actual`: `umol_nauty_sys::run` with actual node colors; `total`: `Graph::automorphisms`:

| Case | $n$ / $m$ | `rust_prep` (µs) | `c_unique` (µs) | `c_actual` (µs) | `total` (µs) | `(prep+unique) / total` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `propane_uniform` | 11 / 10 | 0.174 | 0.509 | 5.213 | 5.574 | 12.3% |
| `naphthalene_kekule` | 10 / 11 | 0.167 | 0.500 | 1.545 | 1.838 | 36.3% |
| `dimethylbutanol_stereo_levi` | 20 / 23 | 0.261 | 0.594 | 2.848 | 3.404 | 25.1% |
| `dimethylbutanol_stereo_direct_bonds` | 11 / 14 | 0.185 | 0.530 | 2.004 | 2.369 | 30.2% |
| `mol_gen_12_r1` | 12 / 12 | 0.200 | 0.524 | 0.673 | 1.047 | 69.1% |
| `mol_gen_16_r2` | 16 / 17 | 0.249 | 0.557 | 0.768 | 1.216 | 66.3% |
| `mol_gen_20_r3` | 20 / 22 | 0.283 | 0.581 | 0.820 | 1.364 | 63.3% |
| `cubane_2color` | 8 / 12 | 0.157 | 0.494 | 2.468 | 2.728 | 23.9% |

### Head-to-head core benchmark matrix (33 Tier 1--4 cases)

Differential verification passed across all 13,598 simple graphs on $n \le 8$ (`simple-through-8.g6`)
with 2-, 3-, and 4-color edges and all 33 corpus cases (8× relabelings + site stabilizers).

| Case | `full_tot` | `mod_tot` | `S1 swg` | `TrailIr (canon)` | `TrailIr (auto)` | `full / IR` | `mod / IR` | `S1 / IR` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `tier1_ir/propane_uniform` | 8.12 | 6.22 | 6.03 | 0.903 | 0.798 | **8.99×** | 6.89× | **6.68×** |
| `tier1_ir/propane_disconnected` | 7.43 | 6.34 | 6.25 | 0.802 | 0.704 | **9.26×** | 7.91× | **7.79×** |
| `tier1_ir/ethane_radicals` | 5.94 | 5.73 | 5.74 | 0.280 | 0.259 | **21.18×** | 20.46× | **20.47×** |
| `tier1_ir/naphthalene_kekule` | 3.38 | 2.87 | 2.17 | 1.448 | 1.329 | **2.33×** | 1.98× | **1.50×** |
| `tier1_ir/disconnected_rings_kekule` | 5.84 | 5.02 | 4.15 | 2.692 | 2.580 | **2.17×** | 1.86× | **1.54×** |
| `tier1_ir/dimethylbutanol_stereo_levi` | 6.68 | 4.74 | 3.68 | 2.126 | 1.952 | **3.14×** | 2.23× | **1.73×** |
| `tier1_ir/dimethylbutanol_stereo_direct_bonds` | 4.50 | 3.70 | 2.67 | 0.823 | 0.692 | **5.47×** | 4.50× | **3.24×** |
| `tier1_ir/overlay_heavy_incidence` | 4.52 | 3.65 | 1.87 | 0.439 | 0.425 | **10.30×** | 8.31× | **4.26×** |
| `tier1_ir/porphyrin_core_3color` | 5.78 | 4.84 | 3.08 | 3.079 | 2.835 | **1.88×** | 1.57× | **1.00×** |
| `tier2_gen/mol_gen_12_r1` | 2.26 | 1.77 | 1.19 | 0.331 | 0.321 | **6.83×** | 5.35× | **3.60×** |
| `tier2_gen/mol_gen_16_r2` | 3.09 | 2.01 | 1.38 | 0.460 | 0.442 | **6.72×** | 4.37× | **3.00×** |
| `tier2_gen/mol_gen_20_r3` | 3.60 | 2.48 | 1.57 | 0.541 | 0.526 | **6.66×** | 4.58× | **2.90×** |
| `tier2_gen/mol_gen_40_r4` | 6.97 | 4.08 | 3.24 | 1.123 | 1.099 | **6.21×** | 3.63× | **2.89×** |
| `tier2_gen/mol_gen_80_r6` | 14.79 | 9.56 | 7.23 | 2.644 | 2.618 | **5.59×** | 3.62× | **2.73×** |
| `tier2_gen/cubic_16_2color` | 3.42 | 2.43 | 1.63 | 0.732 | 0.688 | **4.67×** | 3.32× | **2.23×** |
| `tier2_gen/cubic_24_3color` | 4.78 | 4.02 | 2.47 | 1.369 | 1.358 | **3.49×** | 2.94× | **1.80×** |
| `tier2_gen/cubic_60_2color` | 11.39 | 8.56 | 5.70 | 3.196 | 3.174 | **3.56×** | 2.68× | **1.78×** |
| `tier3_cage/cubane_2color` | 4.56 | 3.54 | 3.00 | 1.713 | 1.606 | **2.66×** | 2.07× | **1.75×** |
| `tier3_cage/adamantane_2color` | 4.21 | 3.29 | 2.40 | 1.270 | 1.161 | **3.32×** | 2.59× | **1.89×** |
| `tier3_cage/dodecahedron_2color` | 12.24 | 8.43 | 6.97 | 5.895 | 5.704 | **2.08×** | 1.43× | **1.18×** |
| `tier3_cage/c60_kekule` | 29.96 | 19.29 | 16.87 | 15.684 | 15.049 | **1.91×** | 1.23× | **1.08×** |
| `tier3_cage/c60_3color` | 11.38 | 9.50 | 5.24 | 3.220 | 3.196 | **3.53×** | 2.95× | **1.63×** |
| `tier3_cage/c70_kekule` | 32.52 | 21.96 | 15.66 | 11.700 | 10.884 | **2.78×** | 1.88× | **1.34×** |
| `tier4_torture/grid_8x8_axes` | 20.68 | 15.88 | 11.32 | 10.705 | 9.887 | **1.93×** | 1.48× | **1.06×** |
| `tier4_torture/grid_8x8_3color` | 13.67 | 11.06 | 5.80 | 3.534 | 3.519 | **3.87×** | 3.13× | **1.64×** |
| `tier4_torture/hypercube_6_dims` | 63.43 | 54.99 | 23.10 | 21.019 | 19.712 | **3.02×** | 2.62× | **1.10×** |
| `tier4_torture/hypercube_6_2color` | 74.89 | 49.58 | 31.07 | 33.273 | 31.889 | **2.25×** | 1.49× | 0.93× |
| `tier4_torture/petersen_2color` | 5.58 | 4.30 | 3.42 | 2.215 | 2.089 | **2.52×** | 1.94× | **1.54×** |
| `tier4_torture/four_hexagons_kekule` | 18.44 | 15.56 | 12.62 | 9.699 | 9.423 | **1.90×** | 1.60× | **1.30×** |
| `tier4_torture/k10_2factor` | 8.79 | 4.93 | 4.13 | 3.030 | 2.723 | **2.90×** | 1.63× | **1.36×** |
| `tier4_torture/rook_4x4_axes` | 17.73 | 13.27 | 8.48 | 7.786 | 7.444 | **2.28×** | 1.70× | **1.09×** |
| `tier4_torture/shrikhande_2color` | 11.85 | 7.27 | 5.30 | 4.488 | 4.155 | **2.64×** | 1.62× | **1.18×** |
| `tier4_torture/miyazaki_ladder_20` | 9.09 | 6.22 | 4.49 | 3.868 | 3.632 | **2.35×** | 1.61× | **1.16×** |

Three-site stabilizer sequences (µs, 3 calls):

| Case | `full_sub` | `mod_sub` | `S1 swg` | `TrailIr (canon)` | `TrailIr (auto)` | `full / auto` | `S1 / auto` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `dimethylbutanol_stereo_levi` | 13.72 | 9.00 | 8.55 | 4.67 | 4.09 | **3.35×** | **2.09×** |
| `dimethylbutanol_stereo_direct_bonds` | 9.90 | 7.61 | 6.90 | 2.02 | 1.75 | **5.65×** | **3.94×** |
| `cubane_2color` | 6.07 | 4.58 | 4.41 | 1.87 | 1.50 | **4.05×** | **2.94×** |
| `dodecahedron_2color` | 14.16 | 9.18 | 7.91 | 5.05 | 4.40 | **3.22×** | **1.80×** |
| `c60_kekule` | 41.94 | 27.20 | 22.59 | 17.55 | 15.06 | **2.78×** | **1.50×** |
| `four_hexagons_kekule` | 34.12 | 28.29 | 24.98 | 18.45 | 17.65 | **1.93×** | **1.42×** |

### Tier-5 complex chemistry stress tests (multi-ligand NUCR vs Bliss vs nauty)

`bliss_nucr`: recursive NUCR without twin collapse; `trail_no_nu`: twin collapse without NUCR;
`trail_nucr`: twin collapse + gated recursive NUCR; `trail_auto`: `trail_nucr` with `get_canon = false`.
Verified across 3× relabelings (exact orbit and normalized group-order agreement up to
$|G| \approx 1.52 \times 10^{51}$ on `cdse_t3_stbu16`):

| Case | $n$ / $m$ | `mod_tot` | `S1 swg` | `bliss_nucr` | `trail_no_nu` | `trail_nucr` | `trail_auto` | `S1 / trail_nucr` | `no_nu / nucr` | `bliss / nucr` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `decapeptide_mixed_explicit_h` | 177 / 180 | 386.9 | 335.1 | 22.9 | 44.5 | **19.1** | 17.9 | **17.53×** | **2.33×** | **1.20×** |
| `decapeptide_cyclic_sym` | 180 / 185 | 385.6 | 333.6 | 45.8 | 76.8 | **36.6** | 34.9 | **9.12×** | **2.10×** | **1.25×** |
| `macrocycle_octaphenyl_porphyrin` | 116 / 128 | 54.4 | 38.0 | 31.3 | 40.9 | **31.9** | 30.9 | **1.19×** | **1.28×** | 0.98× |
| `cdse_t3_sph16` | 206 / 232 | 286.5 | 178.0 | 81.9 | 166.9 | **84.2** | 82.5 | **2.11×** | **1.98×** | 0.97× |
| `cdse_doped_sph16` | 206 / 232 | 220.9 | 136.8 | 38.6 | 129.7 | **38.8** | 36.9 | **3.53×** | **3.34×** | 0.99× |
| `cdse_t3_stbu16` | 238 / 248 | 5739.2 | 5157.4 | 284.3 | 541.6 | **141.8** | 139.0 | **36.38×** | **3.82×** | **2.01×** |

### Sources of `TrailIr`'s speedup relative to nauty

1. **Algorithmic differences:** native edge-colored refinement (`1.14×--2.38×` over modal
   subdivision); exact false-/true-twin collapse (`2.6×--8.9×` on H-/radical-rich molecules,
   `1.20×--2.01×` on top of NUCR on twin-rich complexes); gated recursive NUCR (`1.28×--3.82×` over
   NUCR-free `TrailIr`, up to `36.4×` over weighted nauty on multi-ligand complexes); `FirstSmallest`
   target-cell (`1.06×--1.33×` over `bestcell`/`fsm`); hybrid $O(h)$ sparse `inv_lab` + 1-/2-color
   zero-sort refinement (`1.26×--2.67×` on medium/large and cage graphs); automorphism-only mode
   (`1.03×--1.25×`).
2. **Data layout and memory management:** zero-copy `Graph` CSR reads, LIFO undo trail, flat DFS
   arenas (zero per-node heap allocations), and scratch reuse (`~0.35--1.0 µs` saved on 12--20 atom
   molecules, `1.3×--2.1×`).
3. **FFI / C-boundary elimination:** removes the `0.65--0.86 µs` per-call Rust CSR rebuild + C
   validation + 9 `malloc`/`free` floor (`63%--69%` of total time on 12--20 atom asymmetric graphs).

### End-to-end `umol-graph-ir` trial (`canonicalize` & `constraint_readings`)

Wiring `TrailIr` + native edge colors + sparse RLE `OrderedPartition::refine` + scratch reuse into
`umol-graph-ir` passed `577/577` `proptest` property tests, `1/1` doctest, Clippy `-D warnings`,
and `13,183/13,188` unit tests (the 5 deltas were all solver-basis / tie-break pins: ethane
canonical-label tie-break, `case_3_symmetric_stereo` dropping `24 -> 4` leaves via 3-cycle
generators, and three $S_3$ `ReactionSpan` cases shifting `3 -> 5` leaves).

`umol-graph-ir/benches/canonicalize.rs` (medians in µs):

| Benchmark | Baseline (nauty + sub) | Trial (`TrailIr` + edges) | Change | Speedup |
| --- | ---: | ---: | ---: | ---: |
| `complete/canonicalize/symmetry_heavy_radicals` | 25.51 | 9.91 | −61.5% | **2.57×** |
| `complete/canonicalize/feature_free_disconnected` | 29.53 | 15.08 | −48.5% | **1.96×** |
| `complete/canonicalize/feature_free_connected` | 25.57 | 14.06 | −44.9% | **1.82×** |
| `complete/canonicalize/para_stereo_cascade` | 250.93 | 151.01 | −39.9% | **1.66×** |
| `complete/canonicalize/symmetric_stereo_butanol` | 56.22 | 45.74 | −18.5% | **1.23×** |
| `complete/canonicalize/disconnected_rings` | 23.12 | 19.26 | −16.6% | **1.20×** |
| `complete/canonicalize/symmetric_stereo_dimethylbutanol` | 149.28 | 125.59 | −16.3% | **1.19×** |
| `complete/canonicalize/overlay_heavy` | 18.24 | 16.41 | −10.2% | **1.11×** |
| `complete/canonicalize/mixed_atom_and_bond_stereo` | 14.17 | 12.88 | −8.5% | **1.10×** |
| `complete/canonicalize/overlay_unconstrained` | 13.02 | 12.01 | −8.0% | **1.08×** |
| `complete/canonicalize/frame_relative_stereo_constraint` | 8.55 | 8.02 | −6.0% | **1.07×** |
| `complete/canonicalize/tetrahedral_stereo` | 8.07 | 7.61 | −6.0% | **1.06×** |
| `complete/canonicalize/cis_trans_stereo_bond` | 8.51 | 8.01 | −6.0% | **1.06×** |
| `complete/canonicalize/ordinary_naphthalene` | 12.12 | 11.56 | −5.3% | **1.05×** |
| `complete/canonicalize/aromatic_electron_counts` | 11.91 | 11.32 | −4.9% | **1.05×** |
| `complete/canonicalize/para_stereo_trichloropentane` | 26.68 | 25.90 | −3.6% | **1.03×** |
| `complete/canonicalize/symmetric_charge_sum` | 21.66 | 21.02 | −2.8% | **1.03×** |
| `complete/canonicalize/asymmetric_charge_sum` | 134.80 | 132.50 | −1.9% | **1.02×** |
| `complete/canonicalize/symmetric_stereo_triisopropanolamine` | 35.40 | 35.24 | −0.8% | 1.00× |
| `complete/canonicalize/meso_dichlorobutane` | 19.30 | 19.77 | +2.7% | 0.98× |
| `scaling/canonical_eq_equal/feature_free_connected` | 12.61 | 12.49 | −1.2% | 1.01× |
| `scaling/canonical_eq_equal/feature_free_disconnected` | 12.04 | 11.79 | −1.9% | 1.02× |
| `scaling/canonical_eq_equal/symmetry_heavy_radicals` | 7.56 | 7.59 | +0.5% | 1.00× |
| `reaction_span/canonicalize/multicenter_nonuniform` | 7.64 | 6.40 | −15.5% | **1.19×** |
| `reaction_span/canonicalize/raw_forms` | 22.51 | 20.47 | −8.5% | **1.10×** |
| `reaction_span/canonicalize/entity_constraints` | 21.78 | 20.11 | −8.4% | **1.08×** |
| `reaction_span/canonicalize/aromatic_nonuniform` | 12.68 | 12.41 | −3.3% | **1.02×** |
| `reaction_span/canonicalize/stereo_frame_constraints` | 9.54 | 9.33 | −2.3% | **1.02×** |
| `reaction_span/canonicalize/dative_roles` | 4.77 | 4.69 | −2.2% | **1.02×** |

`umol-graph-ir/benches/constraint_readings.rs` (`stereo_derived_readings/*`):

| Benchmark | Baseline (nauty + sub) | Trial (`TrailIr` + edges) | Change | Speedup |
| --- | ---: | ---: | ---: | ---: |
| `square_planar/topicity` | 7.61 µs | 259.6 ns | −96.6% | **29.3×** |
| `square_planar/stereogenicity` | 15.95 µs | 9.75 µs | −39.5% | **1.64×** |
| `trigonal_bipyramidal/stereogenicity` | 273.57 µs | 184.52 µs | −32.9% | **1.48×** |
| `trigonal_bipyramidal/topicity` | 18.38 µs | 15.98 µs | −12.8% | **1.15×** |
| `meso/stereogenicity` | 12.5 ns | 11.9 ns | −5.1% | **1.05×** |
| `cis_trans/topicity` | 2.91 µs | 2.83 µs | −3.2% | **1.03×** |
| `meso/topicity` | 4.6 ns | 4.5 ns | −2.1% | 1.02× |
| `meso/ligand_symmetry` | 385 ns | 379 ns | −2.2% | 1.02× |
| `cis_trans/ligand_symmetry` | 4.09 µs | 4.07 µs | −0.4% | 1.00× |
| `chlorofluoropropane/ligand_symmetry` | 392 ns | 393 ns | −0.1% | 1.00× |
| `square_planar/ligand_symmetry` | 285 ns | 287 ns | +0.4% | 0.99× |
| `chlorofluoropropane/topicity` | 695 ns | 704 ns | +1.1% | 0.99× |
| `chlorofluoropropane/stereogenicity` | 1.835 µs | 1.867 µs | +1.6% | 0.98× |
| `cis_trans/stereogenicity` | 6.47 µs | 6.65 µs | +2.8% | 0.97× |
| `trigonal_bipyramidal/ligand_symmetry` | 359 ns | 374 ns | +3.9% | 0.96× |

## Staged implementation plan

Every stage ends green (`cargo test`, `cargo clippy --all-targets -- -D warnings`, and Python tests
from `S6b` onward). Breaking public API changes are paired with their caller migrations within the
same stage.

### S0 — `umol-graph-ir` simple-graph prerequisites

Eliminate the two structural non-simple `Graph` representations in `umol-graph-ir` before enforcing
simplicity on `Graph`.

#### S0a — Uniform `StereoAnchor` encoding in `IncidenceGraph` [Done]

**Module:** `umol-graph-ir/src/ir/{incidence.rs, symmetry.rs, molecule/symmetry.rs,
reaction_span/symmetry.rs}`, affected `umol-graph-ir` tests; **breaking (red→green)** (internal
`Incidence` migration within `umol-graph-ir`). `[dep: none]`

- Add `StereoVirtualLigands` (3-bit flag value for `ImplicitHydrogen`, `LonePair`, `Virtual`).
- Replace `Incidence::StereoLigand(StereoLigandKind)` with `Incidence::StereoAnchor(StereoVirtualLigands)`
  (on `(StereoAtom, site_atom)` and `(StereoBond, endpoint_atom)`), unit `Incidence::StereoSite` (on
  `(StereoBond, site_bond)`), and unit `Incidence::StereoLigand` (on actual-atom ligands only).
- Update `Molecule::incidence_graph`, `ReactionSpan::incidence_graph`, incidence color readers, and
  `incidence.rs` / `canonicalize` / `symmetry` unit tests; assert every `IncidenceGraph` is simple.
- **Outcome:** Landed as planned with no extra edge cases; all `IncidenceGraph` instances are simple
  with zero added vertices.

#### S0b — Unique bond-endpoint integrity check on `ReactionSpan` [Done]

**Module:** `umol-graph-ir/src/ir/{reaction_span.rs, reaction.rs, reaction_span/canonicalize.rs,
reaction_span/canonicalize/tests.rs}`, `docs/development/integrity.md`; **breaking (red→green)**
(integrity error variant + rule-level bond-endpoint uniqueness guard + affected tests updated
together). `[dep: none]`

- Add `ReactionSpanIntegrityError::ParallelBonds { atoms: [AtomId; 2] }`, enforced in
  `check_reaction_span_entries` (inherited by `ReactionSpan::from_entries`, `superimpose`, and
  `Reaction::to_reaction_span`).
- Add provisional rule-level `bond_endpoints_are_disjoint` check in `Reaction::prepare_deltas`
  (`ApplyPreconditionError::InconsistentReaction`) and `Reaction::apply_deltas_at`
  (`ApplyError::Inconsistent`) with a TODO to revisit once delta normalization disposes of
  same-endpoint remove+add combinations (leaving `Deltas::normalize` untouched).
- Update [integrity.md](../docs/development/integrity.md), append dated addendum to
  [238](238-non-cascading-removal-2026-10-16.md), and update affected tests; follow-up delta
  normalization is tracked in [255](255-reaction-span-bond-normalization-2026-11-26.md).
- **Findings:** `Reaction::test_reaction_tracked_apply_at_aromatic_replacement` asserted same-endpoint
  `RemoveBond`+`AddBond` replacement identity, contradicting structural bond uniqueness in
  [integrity.md](../docs/development/integrity.md); it was removed and the broader reconciliation
  between structural participant-key identity and doc [238](238-non-cascading-removal-2026-10-16.md)
  same-key `Remove`+`Add` semantics across bonds and overlays was recorded in
  [255](255-reaction-span-bond-normalization-2026-11-26.md).

### S1 — Strictly simple `Graph` and `DiGraph` (`umol-graph-core`)

Retire obsolete multigraph/remnant APIs and close `Graph` and `DiGraph` to simple graphs.

#### S1a — Remove multigraph fallbacks, `SubdividedGraph`, `NonSimpleGraphError`, and bound remnants [Done]

**Module:** `umol-graph-core/src/{graph.rs, lib.rs, correspondence.rs, algorithms/}`,
`umol-graph-core/{tests/property/, benches/algorithms.rs}`, `umol-graph-ir/src/ir/{symmetry.rs,
molecule/canonicalize.rs, reaction_span/canonicalize.rs}`, `docs/development/{nomenclature.md,
property-tests.md}`; **breaking (red→green)** (independently green commit). `[dep: S0a, S0b]`

- Replace all `node_bound()` / `edge_bound()` calls inside `umol-graph-core` with `node_count()` /
  `edge_count()`, and remove `Graph::{is_dense, node_bound, edge_bound, subdivide_edges}`,
  `SubdividedGraph`, `SubdivisionNodeSource`, and `NonSimpleGraphError`.
- Remove the 8 `try_*` / `*_fallback` cycle methods in `cycles.rs` (keeping the 4 direct
  `visit_*` / `enumerate_*` methods), and drop loop/subdivision paths in Read--Tarjan, Horton MCB,
  Vismara, and Kolodzik URF.
- Simplify `Correspondence::edge_matched_pairs` and `shared_edge_count` by using `Graph::find_edge`
  directly and dropping the parallel-edge ambiguity checks and multigraph test cases.
- Replace `SubdivisionNodeSource` in `umol-graph-ir/src/ir/symmetry.rs` with a temporary private
  `AdapterNodeSource` enum until `S5` removes subdivision completely.
- Migrate `tests/property/{cycles.rs, cycles/captured.rs, cycles/exhaustive.rs, correspondence.rs}`
  and `benches/algorithms.rs` to simple graphs, remove `multigraph-through-4-edges-5.tsv`, and
  update [nomenclature.md](../docs/development/nomenclature.md) and
  [property-tests.md](../docs/development/property-tests.md).
- **Findings:** Reordered ahead of `S1b` so both `S1a` and `S1b` remain independently green commits
  without interleaving obsolete multigraph/subdivision APIs with strict simplicity enforcement;
  removing the multigraph cycle corpus also retired `generate_captured_results.py`'s multigraph
  fixture generation.

#### S1b — `GraphError`, checked/asserted simplicity on `Graph` and `DiGraph`, and editing alignment [Done]

**Module:** `umol-graph-core/src/{graph.rs, digraph.rs, lib.rs,
algorithms/{topological_sort.rs, connectivity/components.rs, traversal.rs, matching/maximum.rs}}`,
`umol-graph-core/{tests/property/, benches/restore.rs}`, `umol-graph-ir/src/ir/{molecule.rs,
reaction_span.rs, incidence.rs, symmetry.rs, molecule/, reaction_span/}`,
`umol-graph-ir/tests/property/edit.rs`, `umol-graph/src/ops/aromaticity/hmo.rs`,
`umol-io/src/smiles/render/stereo.rs`, `umol-py/{src/, tests/}`, `docs/development/integrity.md`;
**breaking (red→green)** (independently green commit). `[dep: S1a]`

- Add `GraphError::{NodeOutOfBounds, SelfLoop, ParallelEdge}` (re-exported at crate root).
- Add `Graph::{try_new, try_add_edge, try_extend_edges, try_extend}` and `DiGraph::try_new`; make
  asserted `Graph::{new, add_edge, extend_edges, extend}` and `DiGraph::new` enforce the same
  contract (rejecting out-of-range endpoints, self-loops, and parallel edges within batches or
  against existing edges), keep `Graph::restore` panic-free on malformed inputs, and remove
  `Graph::is_simple`.
- Migrate `topological_sort.rs` self-loop tests to 2-cycles and migrate remaining `umol-graph-core`
  unit tests, `benches/restore.rs`, and property tests (`graph.rs`, `restore.rs`, `connectivity.rs`,
  `traversal.rs`) to simple graphs; add `rstest` unit and `proptest` property tests covering all
  `GraphError` variants and asserted panic boundaries on `Graph` and `DiGraph`.
- Align `umol-graph-ir` and `umol-py`: preserve `Molecule::try_from_entries` and
  `ReactionSpan::try_from_entries` error contracts; make `Edit::AddBonds` reject self-loops and
  parallel bonds before mutating `Graph` with `MoleculeApplyError::Integrity(DuplicateAtom /
  ParallelBonds)`; keep `MoleculeEditor::{add_bond, add_bonds, add_topology}` asserted; update
  [integrity.md](../docs/development/integrity.md) and transient-parallel/self-loop tests in
  `umol-graph-ir` and `umol-py`.
- **Findings:**
  - Enforcing simplicity inside `Graph::try_build_csr`'s existing degree pass and neighbor sort adds
    zero extra heap allocation; `DiGraph::try_new` preserves input successor order while permitting
    opposite-direction arcs `[a, b]` and `[b, a]`.
  - Removing transient parallel/self-loop localized bonds during `Molecule` editing eliminated the
    now-redundant bond simplicity scan in `Molecule::check_integrity`, required `Molecule::apply_edit`
    and `apply_edit_with_undo` to return `Result<_, MoleculeApplyError>` (with direct `.into()`
    conversions at `TransactionError` return sites in `apply.rs` rather than broadening scope to
    redesign internal edit errors), and migrated deferred-integrity tests in `umol-graph-ir` and
    `umol-py` to transient noncovalent-bond violations.
  - Strict simplicity exposed two non-deduplicated edge-list builders in higher crates, both fixed
    by canonicalizing and deduplicating edges before `Graph::new`:
    `HmoAromaticity::recognize_with_exclusions` (`umol-graph/src/ops/aromaticity/hmo.rs`, shared
    bonds across fused rings) and `derive_marker_components`
    (`umol-io/src/smiles/render/stereo.rs`, shared candidate bonds across double-bond stereo sites).

### S2 — Native weighted sparse nauty (`umol-nauty-sys`)

#### S2a — Weighted sparse `dispatch_swg` in C [Done]

**Module:** `umol-nauty-sys/src/umol_nauty_swg.{c,h}`, `umol-nauty-sys/build.rs`; **additive
(green)**. `[dep: none]`

- Implement `isautom_swg`, `testcanlab_swg`, `updatecan_swg`, `init_swg`, conservative
  `cheapautom_swg`, `refine_swg`, `targetcell_swg`, and `dispatch_swg` over `sparsegraph.w`
  without modifying vendored upstream nauty C files. Include a brief file-header notice recording
  its derivation from upstream `nausparse.c` (Apache-2.0 §4(b)). Compile cleanly under
  `-Wall -Wextra`.
- **Findings:**
  - Hardened four edge cases while porting `umol_nauty_swg.c` from `scratch/automorphism-edges`:
    guarded `umol_sparsenauty_weighted`'s self-loop scan on `sg->nde > 0` using indexed access
    (avoiding `NULL + 0` pointer arithmetic UB on edgeless graphs), switched `MASH` arithmetic to
    `unsigned long` (avoiding signed 32-bit `long` overflow on LLP64 platforms when
    `wt == INT_MAX - 1`), guarded `init_swg` on `hin != NULL`, and unified the single-/multi-weight
    non-singleton refinement loop without sorting single-weight splitters.
  - Confirmed zero warnings under C11 `-Wall -Wextra`, untouched `umol-nauty-sys/nauty/`, and exact
    symbol isolation (`libumol_nauty_shim.a` exports only `umol_dispatch_swg` and
    `umol_sparsenauty_weighted` from `umol_nauty_swg.o`, with all upstream nauty references
    `umol_nauty_`-prefixed via `build.rs`).

#### S2b — Safe Rust weighted nauty boundary with `get_canon` control [Done]

**Module:** `umol-nauty-sys/{include/umol_nauty.h, src/umol_nauty.c, src/lib.rs}`; **additive
(green)**. `[dep: S2a]`

- Expose `umol_nauty_run_weighted` (taking `bool get_canon`) and safe Rust `WeightedNautyInput` /
  `run_weighted` (with full CSR, weight, and partition validation).
- Add `rstest` unit tables covering empty/singleton graphs, error variants, edge-color symmetry
  breaking, `get_canon = false`, stabilizers, group order, and multi-threaded concurrency.
- **Findings:**
  - Settled `run_weighted(&WeightedNautyInput, bool) -> Result<NautyOutput, NautyError>` keeping
    `NautyOutput.canonical_labels: Vec<u32>` unchanged (empty, zero-allocation when
    `get_canon = false`, passing `NULL` across FFI for `canonical_labels` and `canong`), preserving
    strictly additive `S2b` scope while `S4a` separates `AutomorphismGroup` and `CanonicalLabeling`
    in `umol-graph-core`.
  - Factored the shared validation, allocation, TLS-guard, and cleanup flow in `umol_nauty.c` into
    `static run_sparse_nauty` shared by `umol_nauty_run` and `umol_nauty_run_weighted`, and added
    `_Static_assert(sizeof(sg_weight) == sizeof(int))`.
  - Strengthened unit assertions over the scratch prototype to verify exact directed edge-weight
    preservation on every reported generator, canonical-relabeling invariance across vertex
    renumberings, the `i32::MAX as u32 - 1` weight boundary, `get_canon = true/false` parity across
    all test tables, and mixed-mode multi-threaded concurrency (`69` passing `umol-nauty-sys` unit
    tests).

#### S2c — `umol-nauty-sys` provenance README [Done]

**Module:** `umol-nauty-sys/{README.md, Cargo.toml}`; **additive (green)**. `[dep: S2b]`

- Add a crate-level `umol-nauty-sys/README.md` (and set `readme = "README.md"` in
  `umol-nauty-sys/Cargo.toml`) documenting the vendored upstream nauty `2.9.3` minimal file subset,
  build-time symbol prefixing in `build.rs`, and the umol-authored `umol_nauty.c` /
  `umol_nauty_swg.c` boundary (including `umol_nauty_swg.c`'s adaptation of `nausparse.c`'s sparse
  dispatch callbacks while keeping `umol-nauty-sys/nauty/` unmodified).
- **Findings:**
  - Switched `umol-nauty-sys/Cargo.toml` from `readme.workspace = true` (which previously pointed at
    the workspace root `../README.md`) to crate-local `readme = "README.md"` and verified package
    inclusion with `cargo package -p umol-nauty-sys --list --allow-dirty`.
  - Kept all README links crate-local so the published crate remains self-contained without
    broken workspace-relative paths.

### S3 — Pure-Rust `TrailIr` engine (`umol-graph-core` private modules)

#### S3a — `trail_ir` internal engine and exhaustive differential suite

**Module:** `umol-graph-core/src/algorithms/automorphism/trail_ir.rs` and child modules
`trail_ir/{partition, preprocess, refine, components, search}.rs` (using `module.rs`, not
`mod.rs`); **additive (green)**. `[dep: S1b, S2b]`

- Port the verified pure-Rust trail-backed edge-colored IR engine (zero-copy `Graph` CSR view,
  `uniform_edges` fast path, undo-trail `Partition`, hybrid $O(h)$ sparse / 2-color zero-sort
  Hopcroft refinement handling the full `u32` edge-color range, exact false-/true-twin collapse,
  O(1)-gated recursive NUCR at `nnt >= 8`, `FirstSmallest` target-cell heuristic, flat DFS stacks,
  and `get_canon: bool`).
- Carry unit and differential tests verifying exact orbit, group-order, generator-subgroup,
  relabeling-invariance, and `get_canon = false` parity against weighted `Nauty` across all 13,598
  simple graphs on $n \le 8$ (`simple-through-8.g6`) with multi-color edges and representative
  molecular/cage/complex fixtures.

### S4 — Public edge-colored automorphism API (`umol-graph-core` & `GraphView`)

#### S4a — `AutomorphismAlgorithm::TrailIr`, `AutomorphismGroup`, `CanonicalLabeling`, `AutomorphismScratch`, and `Graph` methods

**Module:** `umol-graph-core/src/{algorithms/automorphism.rs, lib.rs}`,
`docs/development/nomenclature.md`; **breaking (red→green)**. `[dep: S3a]`

- Add `AutomorphismAlgorithm::TrailIr`.
- Replace `AutomorphismOutput` with closed operation-issued `AutomorphismGroup` and
  `CanonicalLabeling`, plus opaque algorithm-independent `AutomorphismScratch`.
- Redesign `Graph::{automorphisms, automorphisms_with_scratch, canonical_labeling,
  canonical_labeling_with_scratch}` over `<NC: Ord, EC: Ord>`, dispatching both `Nauty` (via
  `run_weighted`) and `TrailIr` (reusing `AutomorphismScratch` buffers for both backends).
- Rewrite `Graph::canonical_key` directly over `Graph::canonical_labeling` (eliminating edge
  subdivision).
- Document all public rustdoc contracts (per `data-type-contracts`) and record the uniform
  `*Scratch` / `*_with_scratch` admission policy in
  [nomenclature.md](../docs/development/nomenclature.md).
- Add `rstest` tables covering both `Nauty` and `TrailIr`, scratch reuse across differing graph
  sizes/colorings/algorithms, non-`Copy` color types, `canonical_key` isomorphism/discrimination,
  and concurrency.

#### S4b — `GraphView` alignment and caller migration

**Module:** `umol-graph-ir/src/ir/{view/graph.rs, symmetry.rs, molecule/symmetry.rs,
molecule/tests.rs}`; **breaking (red→green)** (restores workspace green after `S4a`). `[dep: S4a]`

- Replace `AtomAutomorphism` with `AtomAutomorphismGroup` (adding `bond_color` and `generators`)
  and `AtomCanonicalLabeling` (`labels`, `automorphisms`) on `GraphView::{automorphisms,
  canonical_labeling}`.
- Migrate all existing `umol-graph-ir` internal call sites and tests to the new signatures.

### S5 — Native edge-colored `umol-graph-ir` search and stereo symmetry

#### S5a — Zero-subdivision `AutomorphismAdapter`, sparse edge-colored `OrderedPartition::refine`, and `AutomorphismScratch` reuse

**Module:** `umol-graph-ir/src/ir/{symmetry.rs, symmetry/tests.rs, molecule/symmetry.rs}`;
**additive/internal (green)**. `[dep: S4b]`

- Replace `AutomorphismAdapter::new`'s subdivision builder with direct `source.clone()` +
  `edge_colors = initial_colors.incidences.clone()`, removing all subdivision state and
  `node_sources`.
- Replace dense `O(n·c)` `OrderedPartition::refine` with sparse RLE `(cell_col, edge_color, count)`
  signatures in exact descending dense-count-vector order (including `rstest` tables for edge-color
  splitting and active-vs-all equivalence).
- Thread one stack-local `AutomorphismScratch` across `search_partition`
  (`canonical_labeling_with_scratch`) and `visit_remapping_candidates`
  (`automorphisms_with_scratch`), honoring the caller's selected `AutomorphismAlgorithm`.
- Switch `Molecule::graph_symmetry` to `Graph::automorphisms` (`AutomorphismGroup`, `get_canon =
  false`).

#### S5b — Zero-subdivision `topology_search_input` and algorithm-independent test contracts

**Module:** `umol-graph-ir/src/ir/{molecule/canonicalize.rs, molecule/canonicalize/tests.rs,
reaction_span/canonicalize.rs, reaction_span/canonicalize/tests.rs}`; **additive/internal
(green)**. `[dep: S5a]`

- Remove modal bond subdivision from both `topology_search_input` functions, using
  `molecule.raw_graph().clone()` / `span.graph().clone()` directly with `edge_colors`.
- Replace the 5 backend-basis / tie-break assertions (`test_automorphism_adapter_automorphisms`,
  `test_canonicalize_structure_visited_leaves::case_3_symmetric_stereo`, and
  `test_canonical_candidate::{case_12_lhs_configuration, case_13_rhs_configuration,
  case_16_asymmetric_constraint}`) with algorithm-independent semantic invariants and explicit
  `Nauty` vs `TrailIr` canonicalization parity tests.

### S6 — Switch higher-layer defaults to `TrailIr` (Rust, `umol-py`, and docs)

#### S6a — Rust config defaults and `operations.md`

**Module:** `umol-graph-ir/src/ir/molecule/symmetry.rs`, `umol-graph/src/{ops/canonicalize.rs,
ops/transform/stereo.rs, ops/validate/stereo.rs, fingerprint/substructure.rs}`, affected tests,
`docs/development/operations.md`; **breaking (red→green)** (default value change + default-checking
tests updated together). `[dep: S5b]`

- Change default `automorphism_algorithm` from `AutomorphismAlgorithm::Nauty` to
  `AutomorphismAlgorithm::TrailIr` across all higher-layer configs and update
  [operations.md](../docs/development/operations.md).

#### S6b — Python bindings and stub parity (`umol-py`)

**Module:** `umol-py/src/{algorithm.rs, canonicalize.rs, symmetry.rs, transform.rs,
fingerprint/config.rs}`, `umol-py/python/umol/_umol.pyi`, `umol-py/tests/`; **breaking
(red→green)**. `[dep: S6a]`

- Expose `AutomorphismAlgorithm.TrailIr()` alongside `AutomorphismAlgorithm.Nauty()`, update all
  Python config defaults to `AutomorphismAlgorithm.TrailIr()` in lockstep with Rust (per
  [python-api.md](../docs/development/python-api.md)), update `.pyi` stubs and `pytest` / Rust
  binding tests, and verify in `umol-py/.venv` (Python 3.13) with `maturin develop` and
  `pytest -q umol-py/tests`.

### S7 — Final verification gate and deferrable benchmark retention

#### S7a — Full workspace verification and Criterion confirmation

**Module:** workspace-wide verification; **additive (green)**. `[dep: S6b]`

- Run `cargo +nightly fmt --all`, `cargo test --workspace`,
  `cargo test -p umol-graph-core -p umol-graph-ir --features proptest --test property`,
  `cargo clippy --workspace --all-targets -- -D warnings`, and `umol-py` `maturin develop` +
  `pytest -q umol-py/tests`.
- Re-run `umol-graph-ir` `canonicalize` and `constraint_readings` Criterion benchmarks to confirm
  production numbers match exploratory evidence; close out doc `254` in `discussion/000-status.md`.

#### S7b — Optional permanent edge-colored benchmark suite (deferrable)

**Module:** `umol-graph-core/benches/`; **additive (green)**. `[dep: S7a]`

- Deferrable: promote a compact subset of the Tier 1--5 edge-colored benchmark corpus from
  `scratch/automorphism-edges` into `umol-graph-core/benches/` if ongoing tracking of `TrailIr` vs
  `Nauty` is desired. Core functionality and end-to-end benchmarks work without it.

### Critical path and deferrable stages

- **Critical path:** `S0a/S0b -> S1a -> S1b -> S3a -> S4a -> S4b -> S5a -> S5b -> S6a -> S6b -> S7a`
  (`S2a -> S2b` runs independently of `S0`/`S1` and joins at `S3a`).
- **Deferrable:** `S7b` (permanent core edge-colored benchmark suite).
- **Tracked follow-up (separate doc):** [255](255-reaction-span-bond-normalization-2026-11-26.md)
  covers `Reaction` / `ReactionSpan` normalization of removed+added bonds on the same atom pair and
  revisiting `ReactionSpanIntegrityError::ParallelBonds` once [233](233-graph-ir-delta-review-2026-09-23.md)
  settles delta normalization semantics.
