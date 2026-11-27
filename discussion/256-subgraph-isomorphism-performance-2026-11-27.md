# 256 — Subgraph isomorphism performance

Status: In Progress
Date: 2026-11-27
Relates: [151](151-python-molecule-workflows-2026-07-13.md),
[162](162-common-subgraph-algs-2026-07-25.md),
[207](207-reaction-network-spike-2026-08-24.md),
[254](254-edge-colored-automorphisms-2026-11-26.md),
[255](255-reaction-span-bond-normalization-2026-11-26.md),
[data types guide](../docs/development/data-types.md),
[nomenclature guide](../docs/development/nomenclature.md),
[operations guide](../docs/development/operations.md),
[python api guide](../docs/development/python-api.md)

## Purpose

`umol-graph-core` exposes six subgraph-isomorphism (monomorphism) selectors on `Graph`:
`Vf2`, `Ullmann`, `Ri`, `ArcMatch { path_length }`, `Vf2Rdkit`, and `RayKirsch`. Higher-level
workflows (`PatternFingerprinter` in `umol-graph`, `SubstructureSearchConfig`,
`ReactionApplicationConfig`, and `PatternFingerprintConfig` in `umol-py`) currently default to
`SubgraphIsomorphismAlgorithm::Vf2Rdkit`.

An empirical investigation in `scratch/subgraph-isomorphism` (branched from `6ad23e5d3`) measured
the six existing solvers, 17 experimental solver variants, `umol-graph-core` primitives,
`umol-graph-ir` boundary overhead, and prepared-pattern / prepared-host reuse across:

1. **Core graph regimes (`umol-graph-core`):** 33 multi-regime cases covering molecular cages
   (`naphthalene`, `cubane`, `c60`, `c60_kekule`, `c70`), combinatorial / symmetric / unsatisfiable
   graphs (`c60/C7`, `c60/fused55`, `grid8x8`, `grid10x10`, `hypercube_6`, `rook_4x4`,
   `4hex/2hex`), large labeled chemical complexes (`decapeptide` 177n, `octaphenyl_porphyrin` 116n,
   `cdse_t3_sph16` 206n), Levi incidence graphs (`levi52/pyr10`), and dense Erdős–Rényi graphs
   (`er50_0.4`, `er100_0.25`), across `count`, `first`, and anchored `_at` modes.
2. **Real OpenSMILES corpus (`umol-graph`):** 8,540 parsed molecules (median 17 atoms, p95 39, max
   721; **99.2% have ≤ 64 atoms**; 30.2% carry aromatic overlays) across 8 patterns (`branched`,
   `phenol`, `bicyclic`, `ester`, `sulfonam`, `steroid` [0 hits], `disc_r6c` [disconnected],
   `wild_c6` [all-wildcard 6-ring]), for both `GraphAndOverlays` and `Incidence`, all-embeddings
   and first-match screening.
3. **Many-patterns-against-one-host workloads:** `PatternFingerprinter` (9 wildcard templates × 300
   molecules) and a 24-pattern functional-group / reaction-LHS library (including 6 patterns with
   `RingMembership` constraints) × 500 molecules.

All key algorithmic references were accessible (PDFs in `materials/subgraphs/` plus open-access
papers for Glasgow, VF3L, DAF, and VEQ).

## Empirical findings (2026-11-27, Apple M1 Pro `arm64`)

All 17 experimental solvers passed differential embedding-set parity against `Vf2Rdkit` across 400
random labeled/unlabeled graph pairs (`simple-through-8.g6`, $n \le 6$), empty/singleton/
disconnected/anchored/`first` cases, `grid81` ($65 \le n \le 128$), and `grid144` ($n > 128$).

### 1. Core `umol-graph-core` multi-regime benchmarks (closure boundary)

Selected representative rows from the 33-case matrix (medians, baseline `Graph::find_edge`):

| Case | Hits | `Vf2Rd` | `RK` | `Ri` | `Vf2` | `Ullm` | `Arc6` | `RKO` | `ArcF0` | `Hyb-0` | `Hyb-LA` | `RiL-LA` | `GlasF` | Best vs `Vf2Rd` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| `naphthalene/path3/cnt` | 28 | 1.20 µs | 960 ns | 1.53 µs | 2.57 µs | 6.92 µs | 10.7 µs | **610 ns** | 2.05 µs | 654 ns | 666 ns | 720 ns | 888 ns | `RKO` 1.97× |
| `naphthalene/path3/1st` | 1 | **129 ns** | 179 ns | 304 ns | 194 ns | 671 ns | 6.78 µs | 139 ns | 1.61 µs | 218 ns | 225 ns | 202 ns | 172 ns | `Vf2Rd` 1.00× |
| `cubane/square/cnt` | 48 | 3.33 µs | 2.90 µs | 4.01 µs | 7.69 µs | 22.4 µs | 74.2 µs | 1.88 µs | 4.31 µs | **1.65 µs** | 1.90 µs | 2.15 µs | 2.21 µs | `Hyb-0` 2.02× |
| `c60/hex/cnt` | 240 | 107.4 µs | 93.1 µs | 144 µs | 384 µs | 1.10 ms | 1.79 ms | 54.2 µs | 60.0 µs | 30.5 µs | **27.5 µs** | 32.0 µs | 54.2 µs | `Hyb-LA` 3.91× |
| `c60/hex/1st` | 1 | **432 ns** | 458 ns | 909 ns | 815 ns | 8.54 µs | 536 µs | 490 ns | 18.2 µs | 1.08 µs | 1.09 µs | 788 ns | 1.06 µs | `Vf2Rd` 1.00× |
| `c60_kek/hex_kek/cnt` | 120 | 29.6 µs | 20.9 µs | 32.9 µs | 114 µs | 531 µs | 618 µs | **10.1 µs** | 23.2 µs | 11.4 µs | 12.2 µs | 13.5 µs | 23.1 µs | `RKO` 2.92× |
| `c70/hex/cnt` ($n=70$) | 300 | 117.5 µs | 94.0 µs | 162 µs | 418 µs | 1.39 ms | 2.29 ms | 60.2 µs | 67.2 µs | **38.2 µs** | 38.9 µs | 44.2 µs | 71.8 µs | `Hyb-0` 3.08× |
| `c60/C7/UNSAT` | 0 | 196.5 µs | 162 µs | 255 µs | 741 µs | 1.08 ms | 3.07 ms | 106 µs | 104 µs | 55.2 µs | **45.4 µs** | 53.9 µs | 98.9 µs | `Hyb-LA` 4.33× |
| `c60/fused55/UNSAT` | 0 | 63.0 µs | 51.3 µs | 77.7 µs | 173 µs | 605 µs | 6.76 ms | 31.0 µs | 51.3 µs | 18.8 µs | **15.3 µs** | 18.0 µs | 45.6 µs | `Hyb-LA` 4.12× |
| `grid8x8/C5/UNSAT` | 0 | 102.4 µs | 85.6 µs | 127 µs | 321 µs | 146 µs | 1.43 ms | 56.6 µs | 55.2 µs | 22.4 µs | **21.1 µs** | 24.6 µs | 35.8 µs | `Hyb-LA` 4.85× |
| `grid10x10/C5/UNSAT` | 0 | 168.3 µs | 149 µs | 201 µs | 606 µs | 344 µs | 3.80 ms | 95.9 µs | 89.8 µs | **41.1 µs** | 42.4 µs | 48.0 µs | 67.0 µs | `Hyb-0` 4.10× |
| `q6/K4/UNSAT` | 0 | 32.0 µs | 28.4 µs | 33.5 µs | 15.9 µs | 123 µs | 844 µs | 22.3 µs | 43.4 µs | 4.73 µs | **2.66 µs** | 2.88 µs | 7.26 µs | `Hyb-LA` 12.0× |
| `rook4x4/K5/UNSAT` | 0 | 35.1 µs | 30.0 µs | 41.0 µs | 3.15 µs | 208 µs | 5.15 ms | 22.4 µs | 28.6 µs | 8.57 µs | **862 ns** | 914 ns | 11.8 µs | `Hyb-LA` 40.8× |
| `4hex/2hex/cnt` | 1728 | 254.3 µs | 150 µs | 309 µs | 591 µs | 5.03 ms | 1.62 ms | **92.0 µs** | 106 µs | 109 µs | 120 µs | 145 µs | 220 µs | `RKO` 2.76× |
| `deca177/tyr_ph_o/cnt` | 2 | 5.21 µs | 3.69 µs | 12.5 µs | 43.1 µs | 27.4 µs | 39.1 µs | 4.15 µs | 19.9 µs | 2.71 µs | 2.71 µs | 6.31 µs | **2.65 µs** | `GlasF` 1.97× |
| `deca177/indole/UNSAT` | 0 | 8.78 µs | 5.90 µs | 13.7 µs | 75.7 µs | 7.72 µs | 24.0 µs | 4.11 µs | 2.95 µs | 2.56 µs | 2.65 µs | 6.22 µs | **2.55 µs** | `GlasF` 3.44× |
| `porph116/paraOH/UNSAT` | 0 | 5.32 µs | 3.73 µs | 22.5 µs | 33.8 µs | 6.41 µs | 17.1 µs | 6.65 µs | **1.04 µs** | 1.75 µs | 1.72 µs | 11.9 µs | 1.79 µs | `ArcF0` 5.13× |
| `cdse206/SPh/cnt` | 32 | 27.1 µs | 19.1 µs | 41.9 µs | 294 µs | 337 µs | 422 µs | 11.9 µs | 23.9 µs | 5.12 µs | 5.17 µs | 18.9 µs | **4.97 µs** | `GlasF` 5.45× |
| `cdse206/biph/UNSAT` | 0 | 26.2 µs | 19.2 µs | 4.74 µs | 289 µs | 10.7 µs | 34.5 µs | **2.31 µs** | 4.03 µs | 5.28 µs | 5.24 µs | 3.10 µs | 5.22 µs | `RKO` 11.3× |
| `levi52/pyr10/cnt` | 8 | 4.36 µs | 2.91 µs | 3.94 µs | 13.7 µs | 67.8 µs | 105 µs | **1.65 µs** | 12.5 µs | 2.42 µs | 2.50 µs | 3.43 µs | 2.88 µs | `RKO` 2.64× |
| `er50_0.4/K5/cnt` | 16080 | 12.76 ms | 13.0 ms | 13.4 ms | 26.6 ms | 32.3 ms | 421 ms | 12.8 ms | 6.99 ms | **1.06 ms** | 1.21 ms | 1.34 ms | 1.23 ms | `Hyb-0` 12.1× |
| `er100_0.25/K6` | 1440 | 43.25 ms | 41.7 ms | 42.1 ms | 110 ms | 87.8 ms | 4817 ms | 41.3 ms | 19.2 ms | **2.20 ms** | 2.43 ms | 2.71 ms | 2.66 ms | `Hyb-0` 19.6× |

Key algorithmic conclusions:
- `Ullmann`, `Vf2`, `Ri`, and `ArcMatch` win on zero regimes (`Ullmann` loses 33/33 by 5×–60×; `Vf2` loses 33/33 to `Hyb-LA`; `ArcMatch(6)` is 10×–70,000× slower, and even CSR-rewritten `ArcF6`/`ArcF0` loses 32/33 to `Hyb-LA` with path reduction never paying for itself).
- Pre-search arc consistency (`B-Tp`, `B-Ed`) and DAF failing-set backjumping (`B-TpFS`, `B-EdFS`) lost to `B-No` on 28/33 core cases and 8/8 OpenSMILES patterns.
- `Graph::find_edge`: switching to min-degree endpoint + linear scan for `deg <= 8` regressed uniform-degree graphs by 20%–56% (`c60`, `q6`, `rook4x4`) and `PatternFingerprinter` by 12%; keep `find_edge`'s branchless `binary_search_by_key` unchanged and avoid unnecessary calls via `nb.edge` / bitset adjacency.
- Allocation-only scratch saves only ~3%–5% on $n_2 \le 128$ (where bitset engines stay stack-allocated) and ~15%–20% on $n_2 > 128$, whereas **prepared pattern / prepared host** reuse yields **4×–22×**.

### 2. End-to-end OpenSMILES, `Incidence`, `PatternFingerprinter`, and 24-pattern library

- **OpenSMILES 8,540 molecules (all embeddings, ms, with zero-`BTreeSet` `Correspondence`):**

| Pattern | Hits | `ProdFull` | `FastVer` (`Vf2Rd`) | `PrepFull` (`PrepPat` + fast verify) | `HybLA-Raw` | `PrepP-Raw` | `Both-Raw` | End-to-End Speedup |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `branched` | 15,845 | 53.9 | 41.1 | **11.0** | 15.0 | 8.3 | **2.9** | **4.9×** |
| `phenol` | 7,012 | 77.2 | 77.1 | **11.4** | 20.3 | 9.4 | **4.1** | **6.8×** |
| `bicyclic` | 624 | 39.1 | 41.4 | **7.6** | 24.3 | 7.3 | **3.1** | **5.1×** |
| `ester` | 1,147 | 17.2 | 16.4 | **6.8** | 8.8 | 6.2 | **2.2** | **2.5×** |
| `sulfonam` | 632 | 8.3 | 7.7 | **1.8** | 1.6 | 1.6 | **0.4** | **4.7×** |
| `steroid` | 0 | 39.2 | 43.1 | **5.9** | 27.0 | 5.7 | **3.2** | **6.6×** |
| `disc_r6c` | 32,220 | 244.5 | 191.5 | **25.5** | 31.5 | 18.9 | **15.2** | **9.6×** |
| `wild_c6` | 147,768 | 300.8 | 185.4 | **42.6** | 67.8 | 19.7 | **17.5** | **7.1×** |

- **OpenSMILES first-match screening (ms, raw search):** `PrepP` achieves **2.2×–12.8×** and `PrepPat + PrepHost` achieves **5.9×–21.6×** (`sulfonam` 7.6 → 0.35 ms, `disc_r6c` 141.4 → 7.3 ms, `phenol` 54.6 → 3.5 ms, `steroid` 45.5 → 3.3 ms); one-time `PreparedHost` build for all 8,540 molecules costs **5.47 ms (0.64 µs/mol)**.
- **`PatternFingerprinter` (9 templates × 300 mols, bit-exact):** 103.4 ms → **24.4 ms (4.2×)** with fast `verify_overlays` + zero-`BTreeSet` `Correspondence` alone → **5.38–6.05 ms (17×–19×)** with `PreparedHost` + `PreparedPattern` + borrowed-embedding visitor.
- **24-pattern library (6 with `RingMembership` constraints) × 500 mols:** 67.4 ms → **12.4 ms (5.4×)** full `Vec<MoleculeCorrespondence>` (first-match **10.1 ms**, **6.7×**).

## Settled decisions (SIN)

1. **Semantics and algorithmic transparency:**
   - Monomorphism (non-induced) semantics, deterministic traversal, and explicit algorithm selection at `umol-graph-core` and `umol-graph-ir` are preserved.
   - Subgraph-isomorphism solvers never invoke `AutomorphismAlgorithm` implicitly; symmetry constraints are passed explicitly as data values in `umol-graph-core` and configured via `automorphism_algorithm` on `SubstructureMatchConfig` in `umol-graph-ir`.
2. **`SubgraphIsomorphismAlgorithm` portfolio (4 variants):**
   - **Retire:** `Vf2`, `Ullmann`, `Ri`, `ArcMatch { path_length }`, and `ARCMATCH_DEFAULT_PATH_LENGTH`.
   - **Retain / upgrade:** `Vf2Rdkit` (reference; cleaned static BFS order + flat parent CSR + direct parent edge), `RayKirsch` (upgraded RI-guided `Seed`/`Extend`/`Close` plan + `deg(q) <= deg(t)`).
   - **Add:** `BitsetRi` (RI-DS rarity+connectivity ordering + `u64`/`u128`/wide bitsets + lazy per-row adjacency + direct sparse single-parent scan + 2-instruction bit-parallel monomorphism 1-lookahead) and `Glasgow` (compact-trail bit-parallel forward checking + dynamic MRV + unit propagation + Hall/AllDifferent pruning).
   - **Higher-level defaults:** switch `PatternFingerprinter`, `SubstructureSearchConfig`, `ReactionApplicationConfig`, and `PatternFingerprintConfig` from `Vf2Rdkit` to `BitsetRi`.
3. **Core primitives & predicate boundary (`umol-graph-core`):**
   - Add `Graph::max_degree(&self) -> usize` (cached `u32` on `Csr`). Keep `Graph::find_edge` unchanged.
   - Optimize `Correspondence::new` / `Correspondence::from_images` (zero-heap inline validation for `len <= 32`, empty fast-return, skip sort when sorted).
   - Introduce `SubgraphIsomorphismPredicate` trait (using `undetermined`, not `any`: `all_edges_undetermined`, `edge_is_undetermined`, query node equivalence classes, bulk `u64`/`u128` class candidate fill, optional edge-neighbor bitmask hook) plus an internal `pub(crate)` closure adapter forwarding `Graph`'s closure methods into the same engine.
   - Add `SubgraphIsomorphismScratch`, `PreparedSubgraphQuery` (alg-tagged plan + `QuerySymmetryConstraints`), and `PreparedSubgraphTarget` (alg-independent degree/adjacency bitsets).
   - Add mode splits: `has_*`, `first_*`, `count_*` (with leaf `POPCNT` shortcut for undetermined leaf edges) alongside `visit_*` and `enumerate_*` (including `_at` variants on `Graph` and scratch-taking methods on `PreparedSubgraphQuery`, without duplicating `_with_scratch` closure overloads on `Graph`).
4. **Symmetry reduction (push-down + exact):**
   - Two orthogonal booleans on `SubstructureMatchConfig` (and Python configs): `reduce_pattern_symmetry: bool` and `reduce_host_symmetry: bool`, alongside `automorphism_algorithm: AutomorphismAlgorithm`.
   - Query/reaction symmetry ($G_Q$): Grochow–Kellis stabilizer-chain strict inequalities ($\phi(b_i) < \phi(u)$) pushed down into all 4 solvers; deterministic base priority `(Reverse(deg), q)` makes the constraint set and emitted canonical orbit representatives strictly invariant across `TrailIr`/`Nauty` and all 4 `SubgraphIsomorphismAlgorithm`s.
   - Host symmetry ($G_H$): push-down root orbit-rep restriction + shallow point-stabilizer orbit restriction while stabilizer is non-trivial (zero search work on asymmetric hosts $|G_H| = 1$).
   - Combined $G_Q \times G_H$: simultaneous push-down + leaf-only exact canonicality guard active only when both $|G_Q| > 1$ and $|G_H| > 1$.
   - Proper stereo-, overlay-, and constraint-preserving automorphisms via existing `umol-graph-ir` `molecule/symmetry.rs` and `reaction_span/symmetry.rs` (never improper mirror symmetries).
   - On `Reaction::apply*` / `React::react`, `reduce_pattern_symmetry = true` uses $\text{Aut}(\text{Reaction})$ via `Reaction::to_reaction_span()` projected to LHS atoms (conservative `{1}` fallback if span conversion fails), never bare $\text{Aut}(\text{LHS})$.
5. **Prepared carriers (`umol-graph-ir` and `umol-py`):**
   - `PreparedSubstructurePattern` (immutable snapshot, `try_new` validates `MoleculeScopeConstraints` once, binds `SubstructureMatchConfig`, deduplicates `AtomForm` classes, precomputes `u128`/`u8` element/bond/overlay masks, `QuerySymmetryConstraints`, and optional Levi prep).
   - `PreparedSubstructureHost` (immutable `&self`, `Send + Sync`, frozen in Python; eagerly builds ~0.6 µs `PreparedSubgraphTarget` + element/bond/overlay bitmasks; lazily populates `OnceLock` slots for `IncidenceGraph(Constitution)`, `RingSet(Relevant, 22)` per `RelevantCycleEnumerationAlgorithm`, and stereo-preserving `TargetSymmetry` per `AutomorphismAlgorithm`).
   - Add `Reaction::apply_to_host`, `tracked_apply_to_host`, `apply_to_host_to_reaction`, `apply_to_host_to_reaction_span`, and `impl React for PreparedSubstructureHost`.
   - `PreparedReaction` (caching pre-normalized `Deltas`) is explicitly out of scope and noted in [255](255-reaction-span-bond-normalization-2026-11-26.md).
6. **Caller optimization & Python defaults:**
   - Fast `induced_bonds` (`find_edge`) + empty-overlay short-circuit in `verify_overlays` + borrowed-embedding visitor `visit_substructure_embeddings`.
   - `PatternFingerprinter` precompiles static templates and uses `PreparedSubstructureHost` + `visit_substructure_embeddings` (always `reduce_pattern_symmetry = false, reduce_host_symmetry = false`, RDKit `uniquify = false` bit-exact).
   - Python defaults (`SubstructureSearchConfig`, `ReactionApplicationConfig`): `subgraph_isomorphism_algorithm = BitsetRi()`, `automorphism_algorithm = TrailIr()`, `reduce_pattern_symmetry = False`, `reduce_host_symmetry = False` (conservative defaults for now; revisiting `(True, False)` / `(True, True)` deferred until usage patterns are clearer).
   - Python bindings expose updated configs, `Molecule::{has,first,count}_substructure_match*`, frozen `PreparedSubstructurePattern` and `PreparedSubstructureHost`, and `PreparedSubstructureHost.react` / `Reaction.apply_to_host*`.

## Data-type contract sheets

### `QuerySymmetryConstraints` (`umol-graph-core`)
- **Type and role:** Operation-issued / validated carrier of Grochow–Kellis strict-order symmetry-breaking inequalities over a query node space `0..node_count`.
- **Open carrier or operation-issued value:** Validated value (constructed by `none`, `from_graph*`, or checked `try_from_generators` for external proper subgroups such as stereo-aware `Molecule`/`ReactionSpan` groups).
- **Intrinsic representation invariants:** Every constraint `(u, v)` satisfies `u.index() < node_count`, `v.index() < node_count`, `u != v`, and forms a strict DAG ordered along the base priority `(Reverse(deg), q)`. Stored generator list (used only for combined $G_Q \times G_H$ leaf canonicality) contains valid bijections on `0..node_count`.
- **Contextual properties and supplied context:** Pairing with a `query: &Graph` requires `node_count == query.node_count()` and that represented permutations preserve `query` edges and predicates.
- **Public constructors:** `none(node_count)`, `from_graph` / `from_graph_with_scratch`, `try_from_generators(query, generators)` (rejects out-of-range, non-bijective, or non-edge-preserving generators).
- **First public consumer requiring context:** `PreparedSubgraphQuery::with_symmetry` / `try_with_symmetry` (validates `node_count` match; mismatch returns error or panics on asserted constructor without indexing panic).
- **Semantic properties:** For fixed `(query, node_color, edge_color)`, `from_graph` yields identical constraints under `TrailIr` and `Nauty`. Enforcing these inequalities emits exactly one canonical embedding per $G_Q$-orbit.
- **Rust/Python boundary:** Rust-only (`umol-graph-core` internal/public primitive; Python configures symmetry via booleans + `AutomorphismAlgorithm` on configs).

### `TargetSymmetry` (`umol-graph-core`)
- **Type and role:** Operation-issued host symmetry carrier wrapping a target's `AutomorphismGroup` (plus bitmask of minimum orbit representatives and stabilizer helper data).
- **Open carrier or operation-issued value:** Operation-issued value (`from_graph` / `from_graph_with_scratch` / `from_automorphism_group`).
- **Intrinsic representation invariants:** `orbits` has length `node_count`; `orbit_reps_mask` marks every node `v` with `orbits[v] == v`; every generator is a bijection on `0..node_count`.
- **Contextual properties and supplied context:** Valid for the target `Graph` (and coloring/stereo context) that produced it; checked for `node_count == target.node_count()` on entry to match operations so mismatched context cannot panic.
- **Semantic properties:** Orbit representatives and emitted canonical representatives under $G_H$-reduction are invariant across `TrailIr` and `Nauty`.
- **Rust/Python boundary:** Rust-only (hidden behind `PreparedSubstructureHost` in `umol-graph-ir` / `umol-py`).

### `PreparedSubgraphQuery` & `PreparedSubgraphTarget` & `SubgraphIsomorphismScratch` (`umol-graph-core`)
- **Type and role:**
  - `PreparedSubgraphQuery`: immutable alg-tagged query plan (`Graph` snapshot + degree histogram + alg-specific plan + `QuerySymmetryConstraints`).
  - `PreparedSubgraphTarget`: immutable alg-independent target index (`Graph` snapshot + degree histogram + `deg_ge` / `adj` bitsets).
  - `SubgraphIsomorphismScratch`: opaque stateless reusable buffer carrier (`Clone + Debug + Default`), obeying the `nomenclature.md` `*Scratch` policy.
- **Intrinsic invariants:** Snapshotted `Graph` (`Arc`-shared) guarantees stored degrees, plans, and bitsets always match their owned `Graph` — impossible to desync by mutating the source graph.
- **Semantic properties:** Every `PreparedSubgraphQuery` × `PreparedSubgraphTarget` operation produces the exact same embedding set / boolean / count as the corresponding `Graph` operation with identical algorithm, predicate, and symmetry constraints; `*_with_scratch` is strictly independent of prior scratch calls.
- **Rust/Python boundary:** Rust-only.

### `PreparedSubstructurePattern` (`umol-graph-ir`, bound in `umol-py`)
- **Type and role:** Immutable prepared pattern carrier for 1-pattern × many-hosts and many-patterns × 1-host substructure matching.
- **Open carrier or operation-issued value:** Closed validated carrier (`Clone + Debug + Send + Sync`).
- **Intrinsic representation invariants:** Holds an `Arc`-shared `Molecule` snapshot with zero molecule-scope constraints, its `SubstructureMatchConfig`, consistent deduplicated `AtomForm` classes, `u128`/`u8` pre-screen masks, `PreparedSubgraphQuery` (with proper stereo+overlay-preserving `QuerySymmetryConstraints` if `reduce_pattern_symmetry`), and prebuilt Levi `IncidenceGraph` + query if `match_algorithm == Incidence`.
- **Public constructors:** `try_new(pattern: &Molecule, config: SubstructureMatchConfig) -> Result<Self, SubstructureMatchError>` (and `#[new]` in Python accepting `(pattern, *, config=None)`, raising `ValueError` on `MoleculeScopeConstraints`).
- **Failure behavior:** Construction rejects molecule-scope constraints once; all subsequent matching calls (`matches`, `has_match`, `first_match`, `count_matches`, `visit_*`, `*_in_host`) are infallible.
- **Semantic properties:** `PreparedSubstructurePattern::try_new(p, cfg)?.matches(h)` equals `p.substructure_matches(h, cfg)?` for every `(p, h, cfg)`. Cloning is $O(1)$ / independent; mutating the original `Molecule` after construction does not affect the prepared pattern.
- **Rust/Python boundary:** Immutable (`#[pyclass(frozen)]`) owned object in `umol-py`; `matches`/`has_match`/`first_match`/`count_matches` accept either `Molecule` or `PreparedSubstructureHost` (justified Python union adaptation avoiding duplicated `*_in_host` names in Python).

### `PreparedSubstructureHost` (`umol-graph-ir`, bound in `umol-py`)
- **Type and role:** Immutable prepared host carrier for many-patterns × 1-host and $R \times M$ matching and reaction application.
- **Open carrier or operation-issued value:** Closed immutable snapshot (`Clone + Debug + Send + Sync`) with transparent `OnceLock` lazy caches.
- **Intrinsic representation invariants:** Holds an `Arc`-shared `Molecule` snapshot, matching `PreparedSubgraphTarget`, exact `present_elements: u128`, `present_bond_orders: u8`, `present_overlays: u8`, per-element & per-bond-order bitmasks, and `OnceLock` slots for `IncidenceGraph(Constitution)` (+ Levi `PreparedSubgraphTarget`), `RingSet(Relevant, 22)` per `RelevantCycleEnumerationAlgorithm`, and stereo-preserving `TargetSymmetry` per `AutomorphismAlgorithm`.
- **Public constructors:** `new(host: &Molecule) -> Self` (infallible; `#[new]` in Python taking `molecule: Molecule`).
- **Failure behavior:** Infallible construction; read-only `molecule(&self) -> &Molecule` (independent `Molecule` clone in Python); `React::react` / `Reaction::apply_to_host*` preserve exact `ApplyPreconditionError` / `ApplyError` channels.
- **Semantic properties:** Matching any pattern or reaction against `PreparedSubstructureHost::new(h)` produces the exact same match set, count, existence, and reaction applications as matching directly against `h`. Lazy caches are semantically transparent (observably pure) and thread-safe. Mutating the original `Molecule` after construction does not affect the prepared host.
- **Rust/Python boundary:** Immutable (`#[pyclass(frozen)]`) owned object in `umol-py`; exposes `molecule` property and `react(reaction, *, config=None)`. `Reaction.apply*` methods accept `Molecule | PreparedSubstructureHost` (or `apply_to_host` — matching Rust parity cleanly).

## Staged implementation plan

Every stage ends **green** (`cargo test --workspace` + clippy `-D warnings`, and `maturin develop` + `pytest` once `umol-py` is touched).

### Stage S0 — Zero-API-churn foundations & primitives (additive/internal, green)
- **S0a (`umol-graph-core`) — Done:** Cache `max_degree: u32` in `Csr` and expose `Graph::max_degree(&self) -> usize`; replace `BTreeSet` validation in `Correspondence::new` and `Correspondence::from_images` with inline zero-heap validation for `len <= 32`, empty fast-return, and skip-sort-if-sorted. Unit + property tests. **Additive (green).** `[dep: none]`
  - *Implementation notes:* `Csr::max_degree` is populated in `try_build_csr` and needs no update in `extend_nodes` (appending isolated nodes leaves `max_degree` unchanged). In `Correspondence::new`'s inline `len <= 32` scan, recording `duplicate_right` across `0..i` rather than returning immediately on the first right duplicate preserves exact `DuplicateLeftId`-before-`DuplicateRightId` error precedence (fixing a subtle precedence divergence present in the scratch prototype). Verified with `umol-graph-core` unit tests and `proptest` (`graph::test_graph_max_degree`, `correspondence::test_correspondence_new`, `correspondence::test_correspondence_new_bijection`, and widened `0..48` `test_correspondence_from_images`).
- **S0b (`umol-graph-ir`) — Done:** Rewrite `induced_bonds` to use `right.raw_graph().find_edge` directly (no `BTreeMap`) and short-circuit empty pattern overlay kinds in `verify_overlays`. Unit tests covering overlay-free and overlay-bearing parity. **Additive (green).** `[dep: S0a]`
  - *Implementation notes:* `induced_bonds` validates `atoms` count compatibility (preventing out-of-bounds `find_edge` panics if called directly with mismatched counts) and indexes `atoms.matched_pairs()` directly in $O(1)$ when `atoms.is_total_on_left()` instead of binary-searching `right_of` per bond endpoint. `induce_by_key` also gained an `left_count == 0 || right_count == 0` fast-return, and `verify_overlays` short-circuits both empty pattern overlay kinds (preserving exact `right_count`) and insufficient host overlay counts (`host.<kind>().count() < pattern.<kind>().count()`). Verified green across the full workspace (`umol-graph-ir` `correspondence`/`substructure` proptests + `cargo test --workspace` + `clippy --workspace --all-targets -- -D warnings`).

### Stage S1 — Core predicate trait, symmetry constraints, scratch, and prepared query/target (additive, green)
- **S1a (`umol-graph-core`) — Done:** Define `SubgraphIsomorphismPredicate` trait (using `undetermined` edge hooks, query node equivalence classes, bulk `u64`/`u128` class candidate fill, optional edge-neighbor bitmask hook) and internal `pub(crate)` `ClosureSubgraphIsomorphismPredicate` adapter. Unit tests. **Additive (green).** `[dep: S0a]`
  - *Implementation notes & public-surface tightening:* Placed in `algorithms/subgraph_isomorphism/predicate.rs` (establishing the `subgraph_isomorphism/` submodule layout for `S1`–`S2`). Bulk candidate fills (`node_candidates_u64`, `node_candidates_u128`, `fill_node_candidates_words`) take an `eligible` bitmask so callers can intersect precomputed masks (e.g., `deg_ge[deg(q)]` & element masks) directly and default fallbacks only probe set bits; `edge_neighbors_u64` / `edge_neighbors_u128` return `Option` exact neighbor bitmasks so bond-order bitsets bypass `find_edge` + `match_edge` and preserve leaf `POPCNT`. After reviewing all workspace callers (`umol-graph-ir` `substructure.rs` implements the trait directly; `GraphView` only proxies `Graph`'s closure methods as an escape hatch and can wire additional methods additively later), tightened the public surface by making `ClosureSubgraphIsomorphismPredicate` `pub(crate)` (stripped down to the minimal `(node_match, edge_match)` bridge) and dropping `_with_scratch` closure overloads on `Graph` (keeping `SubgraphIsomorphismScratch` strictly on `PreparedSubgraphQuery`). Verified with 12 unit tests and `clippy --all-targets --features proptest -- -D warnings`.
- **S1b (`umol-graph-core`) — Done:** Implement `QuerySymmetryConstraints` (Grochow–Kellis stabilizer chain with deterministic `(Reverse(deg), q)` base; `none`, `from_graph`, `from_graph_with_scratch`, `try_from_generators`) and `TargetSymmetry` (`none`, `from_graph`, `from_graph_with_scratch`, `from_automorphism_group`, `try_from_generators`, shallow point-stabilizer orbit helper, combined $G_Q \times G_H$ exact leaf canonicality guard). Exhaustive differential tests on `simple-through-8.g6` verifying exact orbit counts and `TrailIr` == `Nauty` invariance. **Additive (green).** `[dep: S0a]`
  - *Implementation notes & findings:* Added `algorithms/subgraph_isomorphism/symmetry.rs` exporting `QuerySymmetryConstraints`, `TargetSymmetry`, and shared `SymmetryGeneratorError` (`LengthMismatch`, `NodeOutOfBounds`, `DuplicateImage`, `EdgeNotPreserved`). Added `TargetSymmetry::none` and `TargetSymmetry::try_from_generators` so `S4b` `PreparedSubstructureHost` can pass proper stereo-preserving host subgroups (which cannot be built as `AutomorphismGroup`). Both carriers compute exact point stabilizers via one-pass Schreier's Lemma + sparse Sims filter + exact disjoint-support commutation pruning (`Support`), requiring only a single `AutomorphismAlgorithm` call per graph and $O(n)$ Schreier-tree storage per active level. Aligning `TargetStabilizerChain` push-down and `is_canonical_double_orbit_embedding` onto the exact same `(Reverse(deg), q)` `query_priority_order` as `QuerySymmetryConstraints` prevents cross-base orbit loss (verified against the 4-node `H = {0-2, 1-3}, Q = {0-1}` counterexample) and guarantees identical canonical representative sets across all four `SubgraphIsomorphismAlgorithm`s and both `AutomorphismAlgorithm`s. Reusing the stabilizer Schreier trees avoids BFS during combined leaf recursion, and `suffix_is_global_min[r]` prevents factorial branching on symmetric cliques (`K_10` in `K_12`). Deliberately omitted `PartialEq` on both carriers so equality does not conflate generator bases or distinct subgroups with matching constraints. Verified in 0.76 s (debug): 13,598/13,598 `simple-through-8.g6` graphs (`TrailIr` == `Nauty`), 52×52 pairs through order 5 (exact orbit counts and unique lex-min representatives across $G_Q$, $G_H$, and $G_Q \times G_H$), $K_{10} \hookrightarrow K_{12}$, 12 disjoint triangles ($|G| = 6^{12}$), wide $n = 144$ bitmasks, and clean `clippy --all-targets --features proptest -- -D warnings`.
- **S1c (`umol-graph-core`) — Done:** Implement `SubgraphIsomorphismScratch`, `PreparedSubgraphTarget` (`t_hist`, `deg_ge`, `adj` `u64`/`u128`/wide), and `PreparedSubgraphQuery` (alg-tagged plan enum + `QuerySymmetryConstraints`). Unit tests. **Additive (green).** `[dep: S1a, S1b]`
  - *Implementation notes & findings:* Added `algorithms/subgraph_isomorphism/prepared.rs` exporting `SubgraphIsomorphismScratch`, `PreparedSubgraphTarget`, and `PreparedSubgraphQuery`. `PreparedSubgraphTarget` stores suffix-sum cumulative `degree_hist` (`hist[d]` = vertices with `deg >= d`) and width-specialized `TargetBitsets` (`U64` for `<= 64`, `U128` for `65..=128`, `Wide` for `> 128`), eagerly materializing wide adjacency matrices up to `WIDE_EAGER_ADJ_LIMIT = 1024` and falling back to lazy per-row scratch adjacency on larger sparse targets to avoid quadratic memory allocation. `PreparedSubgraphQuery` precomputes flat CSR symmetry-bound tables (`sym_lower_bounds` / `sym_upper_bounds` for $O(1)$ bitmask/scalar enforcement), `StaticOrderPlan` for `Vf2Rdkit`, and `RkStep` (`Seed`/`Extend`/`Close`) for `RayKirsch` starting at `priority_order[0]` (`(Reverse(deg), q)`) so `RayKirsch` gets immediate depth-0 host orbit-rep pruning. `try_with_symmetry -> Option<Self>` (paired with asserted `with_symmetry`) validates both `node_count` and `priority_order` compatibility for non-trivial symmetry without introducing a single-variant error type. Verified with 9 unit tests (including snapshot mutation isolation, `n = 0/4/70/144/1050` bitset widths, symmetry CSR bounds, plan invariants, and mismatch rejection) plus full Stage S1 workspace gate (`cargo test --workspace`, `cargo clippy --workspace --all-targets -- -D warnings`, and `clippy -p umol-graph-core --all-targets --features proptest -- -D warnings`).

### Stage S2 — Core engines (`BitsetRi`, `Glasgow`, upgraded `RayKirsch`, cleaned `Vf2Rdkit`) + mode splits (additive, green)
- **S2a (`umol-graph-core`):** Add `BitsetRi` variant and implementation (`u64`/`u128`/wide bitsets, degree-sequence filter, descending-degree class candidate fill, RI-DS rarity+connectivity order, lazy per-row `adj`, direct sparse single-parent scan, 2-instruction bit-parallel 1-lookahead, push-down $G_Q$/$G_H$ symmetry, anchored `_at`, leaf `POPCNT` for `count`). **Additive (green).** `[dep: S1c]`
- **S2b (`umol-graph-core`):** Add `Glasgow` variant and implementation (compact $n_1$-word flat trail `u64`/`u128`/wide, bit-parallel FC, dynamic MRV, unit propagation, Hall/AllDifferent wipeout pruning, push-down $G_Q$/$G_H$ symmetry, anchored `_at`, leaf `POPCNT` for `count`). **Additive (green).** `[dep: S1c]`
- **S2c (`umol-graph-core`):** Upgrade `RayKirsch` (`RkStep` RI-guided plan + `deg(q) <= deg(t)` + symmetry constraints) and clean up `Vf2Rdkit` (static BFS order + flat parent CSR + direct parent edge + symmetry constraints); expose `has_subgraph_isomorphism`, `first_subgraph_isomorphism`, `count_subgraph_isomorphisms` (+ `_at`) on `Graph` and scratch-taking `visit_matches`/`has_match`/`first_match`/`count_matches` on prepared query/target. **Additive (green).** `[dep: S1c]`
- **S2d (`umol-graph-core`):** Exhaustive differential & property test suite verifying exact embedding-set parity across `[Vf2Rdkit, RayKirsch, BitsetRi, Glasgow]`, both `[TrailIr, Nauty]`, all 4 `(pattern_sym, host_sym)` combinations, `visit`/`has`/`first`/`count`, anchored `_at`, and $n \le 64$, $65..=128$, $n > 128$. **Additive (green).** `[dep: S2a, S2b, S2c]`

### Stage S3 — Retire dominated `SubgraphIsomorphismAlgorithm` variants & migrate workspace (breaking red→green)
- **S3a (`umol-graph-core`):** Remove `Vf2`, `Ullmann`, `Ri`, `ArcMatch { path_length }`, and `ARCMATCH_DEFAULT_PATH_LENGTH` from `SubgraphIsomorphismAlgorithm` and `lib.rs`. **Breaking (red until S3b).** `[dep: S2d]`
- **S3b (`umol-graph-ir`, `umol-graph`, `umol-py`, benches, tests, fuzz targets, `README.md`):** Migrate all callers, tests, fuzz targets (`fuzz_reaction.rs`), `umol-py` `SubgraphIsomorphismAlgorithm` bindings, and `umol-graph-core/benches/subgraph_isomorphism.rs` (expanding benchmark coverage to include UNSAT/hard/large/symmetry-reduced cases across the 4 retained algorithms). Workspace green. **Migration (green).** `[dep: S3a]`

### Stage S4 — `umol-graph-ir` config symmetry fields, prepared carriers, and mode splits (breaking config + migration)
- **S4a (`umol-graph-ir`, `umol-graph`, `umol-py`, benches, tests, fuzz):** Add `automorphism_algorithm: AutomorphismAlgorithm`, `reduce_pattern_symmetry: bool`, and `reduce_host_symmetry: bool` to `SubstructureMatchConfig`; migrate all struct literals across the workspace. **Breaking → green.** `[dep: S3b]`
- **S4b (`umol-graph-ir`):** Implement `PreparedSubstructureHost` (immutable snapshot, `OnceLock` lazy `IncidenceGraph`/`RingSet`/stereo-preserving `TargetSymmetry`, `u128`/`u8` masks, per-element & per-bond-order bitsets) and `PreparedSubstructurePattern` (deduplicated `AtomForm` classes, `undetermined` bond masks, `u128`/`u8` pre-screen, proper stereo+overlay-preserving `QuerySymmetryConstraints`, Levi prep). **Additive (green).** `[dep: S4a]`
- **S4c (`umol-graph-ir`):** Add `Molecule::{visit_substructure_embeddings, has_substructure_match, first_substructure_match, count_substructure_matches}` and `GraphView` mode splits; wire `Molecule::visit_substructure_matches` over the fast prepared/predicate engine; add comprehensive unit & property tests (`GraphAndOverlays` == `Incidence`, stereo-preserving symmetry, all 4 symmetry modes, prepared vs direct parity, all 4 algorithms). **Additive (green).** `[dep: S4b]`

### Stage S5 — `Reaction` symmetry soundness & `PreparedSubstructureHost` integration (additive, green)
- **S5a (`umol-graph-ir`):** Wire `Reaction::apply*` / `React::react` so `reduce_pattern_symmetry = true` builds reaction-preserving `QuerySymmetryConstraints` via `Reaction::to_reaction_span()` + `reaction_span/symmetry.rs` (preserving deltas & stereo; conservative `{1}` fallback if span conversion fails), never bare `lhs` symmetry. **Additive (green).** `[dep: S4c]`
- **S5b (`umol-graph-ir`):** Add `Reaction::{apply_to_host, tracked_apply_to_host, apply_to_host_to_reaction, apply_to_host_to_reaction_span}` and `impl React for PreparedSubstructureHost`; add tests for regioisomer preservation on symmetric LHS (e.g. benzene monosubstitution on toluene), chiral stereocenter preservation, symmetric bimolecular `A + A` deduplication, symmetric host deduplication (benzene, biphenyl, adipic acid), and parity with `Reaction::apply`. **Additive (green).** `[dep: S5a]`

### Stage S6 — `umol-graph` `PatternFingerprinter`, defaults, and benchmarks (additive/default switch, green)
- **S6a (`umol-graph`):** Switch `PatternFingerprinter::new()` default to `BitsetRi`; precompile templates and use `PreparedSubstructureHost` + `visit_substructure_embeddings` across both `GraphAndOverlays` and `Incidence` and all 4 `SubgraphIsomorphismAlgorithm` variants; verify 100% bit-exact RDKit 2026.03.3 test parity. **Additive/default switch (green).** `[dep: S4c]`
- **S6b (`umol-graph`, `umol-graph-ir`):** Update `umol-graph/benches/substructure.rs` (4 algorithms + prepared pattern/host + screening + symmetry-reduced) and `umol-graph-ir/benches/reaction.rs` (`BitsetRi` + `apply_to_host` library benchmark); record before/after benchmark verification. **Additive (green).** `[dep: S5b, S6a]`

### Stage S7 — `umol-py` bindings, defaults, documentation, and verification (additive/default switch, green)
- **S7a (`umol-py`):** Update `SubstructureSearchConfig`, `ReactionApplicationConfig`, and `PatternFingerprintConfig` defaults (`BitsetRi()`, `TrailIr()`, `reduce_pattern_symmetry=False`, `reduce_host_symmetry=False`); expose `Molecule::{has_substructure_match, first_substructure_match, count_substructure_matches}`. **Additive/default switch (green).** `[dep: S5b, S6a]`
- **S7b (`umol-py`, docs):** Bind frozen `PreparedSubstructureHost` and `PreparedSubstructurePattern` (supporting `Molecule | PreparedSubstructureHost` on `PreparedSubstructurePattern` matching methods and `Reaction.apply*` / `PreparedSubstructureHost.react`); add Rust PyO3 tests + Python `pytest` suite (`umol-py/.venv`, Python 3.13); update `docs/development/{nomenclature,operations,python-api}.md` and `README.md`. **Additive (green).** `[dep: S7a]`

### Critical path & verification gate
- **Critical path:** `S0a → S1a/S1b → S1c → S2a/S2b/S2c → S2d → S3a → S3b → S4a → S4b → S4c → S5a → S5b → S7a → S7b` (`S6a/S6b` run in parallel with `S5`/`S7a` after `S4c`).
- **Final verification gate:** `cargo +nightly fmt --all -- --check`, `cargo clippy --workspace --all-targets -- -D warnings`, `cargo test --workspace`, feature-gated `proptest` suites for `umol-graph-core` and `umol-graph-ir`, `maturin develop` + `pytest -q umol-py/tests` under `umol-py/.venv` (Python 3.13), and targeted benchmark regression checks.

## Out-of-scope follow-up ideas

Orthogonal extensions and internal refinements that build cleanly on the settled `256` interfaces without changing them (not part of `S0`–`S7`):

1. **`PreparedReaction` (reaction-side preparation for $1 \times M$ and $R \times M$):** Cache both pre-normalized `Deltas` (deferred to [255](255-reaction-span-bond-normalization-2026-11-26.md)) and the reaction's `PreparedSubstructurePattern` + `ReactionSpan`-derived `QuerySymmetryConstraints` together so applying one reaction across many hosts avoids rebuilding the LHS pattern and span symmetry per host.
2. **Indexed / multiplexed batch search (`SubstructurePatternSet` and host screening indices):** An inverted feature/element index over large rule libraries ($R \gg 10^2$–$10^3$) to skip even the ~1–2 ns bitmask pre-screen per rule, and a fingerprint / structural-key pre-screen index over large molecule databases ($M \gg 10^4$).
3. **Component-wise multi-reactant filtering (`A + B` bimolecular enumeration):** Match disconnected reaction LHS components independently against per-molecule `PreparedSubstructureHost` instances before forming reactant pairs, avoiding $O(M^2)$ `Molecule::combine_all` constructions.
4. **Parallel batch wrappers and Python GIL release:** Leverage `PreparedSubstructurePattern` / `PreparedSubstructureHost`'s `Send + Sync` / frozen immutability for `rayon` batch methods in Rust and GIL-releasing single/batch screening calls in `umol-py`.
5. **Extended pre-screen bitmasks:** Add formal-charge, isotope, radical, stereo, or ring-presence summary masks inside `PreparedSubstructurePattern` / `PreparedSubstructureHost` if corpus profiling identifies queries that pass element/bond/overlay screening and reject late.
6. **Optional symmetry-breaking for `has` / `first` and revisiting Python defaults:** Re-evaluate conservative `(False, False)` Python defaults toward `(True, False)` for `SubstructureSearchConfig` and `(True, True)` for `ReactionApplicationConfig` once downstream usage patterns stabilize, and assess whether hard symmetric UNSAT `has`/`first` workloads warrant opt-in symmetry acceleration.
