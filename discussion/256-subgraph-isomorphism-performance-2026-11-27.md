# 256 — Subgraph isomorphism performance

Status: Completed
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

`umol-graph-core` exposes six subgraph-isomorphism (monomorphism) selectors on `Graph`: `Vf2`,
`Ullmann`, `Ri`, `ArcMatch { path_length }`, `Vf2Rdkit`, and `RayKirsch`. Higher-level workflows
(`PatternFingerprinter` in `umol-graph`, `SubstructureSearchConfig`, `ReactionApplicationConfig`,
and `PatternFingerprintConfig` in `umol-py`) currently default to
`SubgraphIsomorphismAlgorithm::Vf2Rdkit`.

An empirical investigation in `scratch/subgraph-isomorphism` (branched from `6ad23e5d3`) measured
the six existing solvers, 17 experimental solver variants, `umol-graph-core` primitives,
`umol-graph-ir` boundary overhead, and prepared-pattern / prepared-host reuse across:

1. **Core graph regimes (`umol-graph-core`):** 33 multi-regime cases covering molecular cages
   (`naphthalene`, `cubane`, `c60`, `c60_kekule`, `c70`), combinatorial / symmetric / unsatisfiable
   graphs (`c60/C7`, `c60/fused55`, `grid8x8`, `grid10x10`, `hypercube_6`, `rook_4x4`, `4hex/2hex`),
   large labeled chemical complexes (`decapeptide` 177n, `octaphenyl_porphyrin` 116n,
   `cdse_t3_sph16` 206n), Levi incidence graphs (`levi52/pyr10`), and dense Erdős–Rényi graphs
   (`er50_0.4`, `er100_0.25`), across `count`, `first`, and anchored `_at` modes.
2. **Real OpenSMILES corpus (`umol-graph`):** 8,540 parsed molecules (median 17 atoms, p95 39, max
   721; **99.2% have ≤ 64 atoms**; 30.2% carry aromatic overlays) across 8 patterns (`branched`,
   `phenol`, `bicyclic`, `ester`, `sulfonam`, `steroid` [0 hits], `disc_r6c` [disconnected],
   `wild_c6` [all-wildcard 6-ring]), for both `GraphAndOverlays` and `Incidence`, all-embeddings and
   first-match screening.
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
- `Ullmann`, `Vf2`, `Ri`, and `ArcMatch` win on zero regimes (`Ullmann` loses 33/33 by 5×–60×; `Vf2`
  loses 33/33 to `Hyb-LA`; `ArcMatch(6)` is 10×–70,000× slower, and even CSR-rewritten
  `ArcF6`/`ArcF0` loses 32/33 to `Hyb-LA` with path reduction never paying for itself).
- Pre-search arc consistency (`B-Tp`, `B-Ed`) and DAF failing-set backjumping (`B-TpFS`, `B-EdFS`)
  lost to `B-No` on 28/33 core cases and 8/8 OpenSMILES patterns.
- `Graph::find_edge`: switching to min-degree endpoint + linear scan for `deg <= 8` regressed
  uniform-degree graphs by 20%–56% (`c60`, `q6`, `rook4x4`) and `PatternFingerprinter` by 12%; keep
  `find_edge`'s branchless `binary_search_by_key` unchanged and avoid unnecessary calls via
  `nb.edge` / bitset adjacency.
- Allocation-only scratch saves only ~3%–5% on $n_2 \le 128$ (where bitset engines stay
  stack-allocated) and ~15%–20% on $n_2 > 128$, whereas **prepared pattern / prepared host** reuse
  yields **4×–22×**.

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

- **OpenSMILES first-match screening (ms, raw search):** `PrepP` achieves **2.2×–12.8×** and
  `PrepPat + PrepHost` achieves **5.9×–21.6×** (`sulfonam` 7.6 → 0.35 ms, `disc_r6c` 141.4 → 7.3 ms,
  `phenol` 54.6 → 3.5 ms, `steroid` 45.5 → 3.3 ms); one-time `PreparedHost` build for all 8,540
  molecules costs **5.47 ms (0.64 µs/mol)**.
- **`PatternFingerprinter` (9 templates × 300 mols, bit-exact):** 103.4 ms → **24.4 ms (4.2×)** with
  fast `verify_overlays` + zero-`BTreeSet` `Correspondence` alone → **5.38–6.05 ms (17×–19×)** with
  `PreparedHost` + `PreparedPattern` + borrowed-embedding visitor.
- **24-pattern library (6 with `RingMembership` constraints) × 500 mols:** 67.4 ms → **12.4 ms
  (5.4×)** full `Vec<MoleculeCorrespondence>` (first-match **10.1 ms**, **6.7×**).

## Settled decisions (SIN)

1. **Semantics and algorithmic transparency:**
   - Monomorphism (non-induced) semantics, deterministic traversal, and explicit algorithm selection
     at `umol-graph-core` and `umol-graph-ir` are preserved.
   - Subgraph-isomorphism solvers never invoke `AutomorphismAlgorithm` implicitly; symmetry
     constraints are passed explicitly as data values in `umol-graph-core` and configured via
     `automorphism_algorithm` on `SubstructureMatchConfig` in `umol-graph-ir`.
2. **`SubgraphIsomorphismAlgorithm` portfolio (4 variants):**
   - **Retire:** `Vf2`, `Ullmann`, `Ri`, `ArcMatch { path_length }`, and
     `ARCMATCH_DEFAULT_PATH_LENGTH`.
   - **Retain / upgrade:** `Vf2Rdkit` (reference; cleaned static BFS order + flat parent CSR +
     direct parent edge), `RayKirsch` (upgraded RI-guided `Seed`/`Extend`/`Close` plan + `deg(q) <=
     deg(t)`).
   - **Add:** `BitsetRi` (RI-DS rarity+connectivity ordering + `u64`/`u128`/wide bitsets + lazy
     per-row adjacency + direct sparse single-parent scan + 2-instruction bit-parallel monomorphism
     1-lookahead) and `Glasgow` (compact-trail bit-parallel forward checking + dynamic MRV + unit
     propagation + Hall/AllDifferent pruning).
   - **Higher-level defaults:** switch `PatternFingerprinter`, `SubstructureSearchConfig`,
     `ReactionApplicationConfig`, and `PatternFingerprintConfig` from `Vf2Rdkit` to `BitsetRi`.
3. **Core primitives & predicate boundary (`umol-graph-core`):**
   - Add `Graph::max_degree(&self) -> usize` (cached `u32` on `Csr`). Keep `Graph::find_edge`
     unchanged.
   - Optimize `Correspondence::new` / `Correspondence::from_images` (zero-heap inline validation for
     `len <= 32`, empty fast-return, skip sort when sorted).
   - Introduce `SubgraphIsomorphismPredicate` trait (using `undetermined`, not `any`:
     `all_edges_undetermined`, `edge_is_undetermined`, query node equivalence classes, bulk
     `u64`/`u128` class candidate fill, optional edge-neighbor bitmask hook) plus an internal
     `pub(crate)` closure adapter forwarding `Graph`'s closure methods into the same engine.
   - Add `SubgraphIsomorphismScratch`, `PreparedSubgraphQuery` (alg-tagged plan +
     `QuerySymmetryConstraints`), and `PreparedSubgraphTarget` (alg-independent degree/adjacency
     bitsets).
   - Add mode splits: `has_*`, `first_*`, `count_*` (with leaf `POPCNT` shortcut for undetermined
     leaf edges) alongside `visit_*` and `enumerate_*` (including `_at` variants on `Graph` and
     scratch-taking methods on `PreparedSubgraphQuery`, without duplicating `_with_scratch` closure
     overloads on `Graph`).
4. **Symmetry reduction (push-down + exact):**
   - Two orthogonal booleans on `SubstructureMatchConfig` (and Python configs):
     `reduce_pattern_symmetry: bool` and `reduce_host_symmetry: bool`, alongside
     `automorphism_algorithm: AutomorphismAlgorithm`.
   - Query/reaction symmetry ($G_Q$): Grochow–Kellis stabilizer-chain strict inequalities
     ($\phi(b_i) < \phi(u)$) pushed down into all 4 solvers; deterministic base priority
     `(Reverse(deg), q)` makes the constraint set and emitted canonical orbit representatives
     strictly invariant across `TrailIr`/`Nauty` and all 4 `SubgraphIsomorphismAlgorithm`s.
   - Host symmetry ($G_H$): push-down root orbit-rep restriction + shallow point-stabilizer orbit
     restriction while stabilizer is non-trivial (zero search work on asymmetric hosts $|G_H| = 1$).
   - Combined $G_Q \times G_H$: simultaneous push-down + leaf-only exact canonicality guard active
     only when both $|G_Q| > 1$ and $|G_H| > 1$.
   - Proper stereo-, overlay-, and constraint-preserving automorphisms via existing `umol-graph-ir`
     `molecule/symmetry.rs` and `reaction_span/symmetry.rs` (never improper mirror symmetries).
   - On `Reaction::apply*` / `React::react`, `reduce_pattern_symmetry = true` uses
     $\text{Aut}(\text{Reaction})$ via `Reaction::to_reaction_span()` projected to LHS atoms
     (conservative `{1}` fallback if span conversion fails), never bare $\text{Aut}(\text{LHS})$.
5. **Prepared carriers (`umol-graph-ir` and `umol-py`):**
   - `PreparedSubstructurePattern` (immutable snapshot, `new` panics on molecule constraints, binds
     `SubstructureMatchConfig`, deduplicates `AtomForm` classes, precomputes `u128`/`u8`
     element/bond/overlay masks, `QuerySymmetryConstraints`, and optional Levi prep).
   - `PreparedSubstructureHost` (immutable `&self`, `Send + Sync`, frozen in Python; eagerly builds
     ~0.6 µs `PreparedSubgraphTarget` + element/bond/overlay bitmasks; lazily populates `OnceLock`
     slots for `IncidenceGraph(Constitution)`, `RingSet(Relevant, 22)` per
     `RelevantCycleEnumerationAlgorithm`, and stereo-preserving `TargetSymmetry` per
     `AutomorphismAlgorithm`).
   - Add `Reaction::apply_to_host`, `tracked_apply_to_host`, `apply_to_host_to_reaction`,
     `apply_to_host_to_reaction_span`, `impl React for PreparedSubstructureHost`, and
     `PreparedReaction` (pre-validated/normalized `Deltas` + `PreparedSubstructurePattern` with
     reaction-preserving `QuerySymmetryConstraints`, infallible `apply*` / `apply_to_host*`
     iterators, and `React::react_prepared`) so repeated reaction application obeys the same
     symmetric "prepare whatever you reuse" model as substructure matching (when
     [255](255-reaction-span-bond-normalization-2026-11-26.md) lands later, it updates delta
     normalization behind `PreparedReaction` without changing its interface).
6. **Caller optimization & Python defaults:**
   - Fast `induced_bonds` (`find_edge`) + empty-overlay short-circuit in `verify_overlays` +
     borrowed-embedding visitor `visit_substructure_embeddings`.
   - `PatternFingerprinter` precompiles static templates and uses `PreparedSubstructureHost` +
     `visit_substructure_embeddings` (always `reduce_pattern_symmetry = false, reduce_host_symmetry
     = false`, RDKit `uniquify = false` bit-exact).
   - Python defaults (`SubstructureSearchConfig`, `ReactionApplicationConfig`):
     `subgraph_isomorphism_algorithm = BitsetRi()`, `automorphism_algorithm = TrailIr()`,
     `reduce_pattern_symmetry = False`, `reduce_host_symmetry = False` (conservative defaults for
     now; revisiting `(True, False)` / `(True, True)` deferred until usage patterns are clearer).
   - Python bindings expose updated configs, `Molecule::{has,first,count}_substructure_match*`,
     frozen `PreparedSubstructurePattern`, `PreparedSubstructureHost`, and `PreparedReaction`, and
     `PreparedSubstructureHost.react` / `Reaction.apply*` / `PreparedReaction.apply*`.

## Data-type contract sheets

### `QuerySymmetryConstraints` (`umol-graph-core`)
- **Type and role:** Operation-issued / validated carrier of Grochow–Kellis strict-order
  symmetry-breaking inequalities over a query node space `0..node_count`.
- **Open carrier or operation-issued value:** Validated value (constructed by `none`, `from_graph*`,
  or checked `try_from_generators` for external proper subgroups such as stereo-aware
  `Molecule`/`ReactionSpan` groups).
- **Intrinsic representation invariants:** Every constraint `(u, v)` satisfies `u.index() <
  node_count`, `v.index() < node_count`, `u != v`, and forms a strict DAG ordered along the base
  priority `(Reverse(deg), q)`. Stored generator list (used only for combined $G_Q \times G_H$ leaf
  canonicality) contains valid bijections on `0..node_count`.
- **Contextual properties and supplied context:** Pairing with a `query: &Graph` requires
  `node_count == query.node_count()` and that represented permutations preserve `query` edges and
  predicates.
- **Public constructors:** `none(node_count)`, `from_graph` / `from_graph_with_scratch`,
  `try_from_generators(query, generators)` (rejects out-of-range, non-bijective, or
  non-edge-preserving generators).
- **First public consumer requiring context:** `PreparedSubgraphQuery::with_symmetry` /
  `try_with_symmetry` (validates `node_count` match; mismatch returns error or panics on asserted
  constructor without indexing panic).
- **Semantic properties:** For fixed `(query, node_color, edge_color)`, `from_graph` yields
  identical constraints under `TrailIr` and `Nauty`. Enforcing these inequalities emits exactly one
  canonical embedding per $G_Q$-orbit.
- **Rust/Python boundary:** Rust-only (`umol-graph-core` internal/public primitive; Python
  configures symmetry via booleans + `AutomorphismAlgorithm` on configs).

### `TargetSymmetry` (`umol-graph-core`)
- **Type and role:** Operation-issued host symmetry carrier wrapping a target's `AutomorphismGroup`
  (plus bitmask of minimum orbit representatives and stabilizer helper data).
- **Open carrier or operation-issued value:** Operation-issued value (`from_graph` /
  `from_graph_with_scratch` / `from_automorphism_group`).
- **Intrinsic representation invariants:** `orbits` has length `node_count`; `orbit_reps_mask` marks
  every node `v` with `orbits[v] == v`; every generator is a bijection on `0..node_count`.
- **Contextual properties and supplied context:** Valid for the target `Graph` (and coloring/stereo
  context) that produced it; checked for `node_count == target.node_count()` on entry to match
  operations so mismatched context cannot panic.
- **Semantic properties:** Orbit representatives and emitted canonical representatives under
  $G_H$-reduction are invariant across `TrailIr` and `Nauty`.
- **Rust/Python boundary:** Rust-only (hidden behind `PreparedSubstructureHost` in `umol-graph-ir` /
  `umol-py`).

### `PreparedSubgraphQuery` & `PreparedSubgraphTarget` & `SubgraphIsomorphismScratch` (`umol-graph-core`)
- **Type and role:**
  - `PreparedSubgraphQuery`: immutable alg-tagged query plan (`Graph` snapshot + degree histogram +
    alg-specific plan + `QuerySymmetryConstraints`).
  - `PreparedSubgraphTarget`: immutable alg-independent target index (`Graph` snapshot + degree
    histogram + `deg_ge` / `adj` bitsets).
  - `SubgraphIsomorphismScratch`: opaque stateless reusable buffer carrier (`Clone + Debug +
    Default`), obeying the `nomenclature.md` `*Scratch` policy.
- **Intrinsic invariants:** Snapshotted `Graph` (`Arc`-shared) guarantees stored degrees, plans, and
  bitsets always match their owned `Graph` — impossible to desync by mutating the source graph.
- **Semantic properties:** Every `PreparedSubgraphQuery` × `PreparedSubgraphTarget` operation
  produces the exact same embedding set / boolean / count as the corresponding `Graph` operation
  with identical algorithm, predicate, and symmetry constraints; `*_with_scratch` is strictly
  independent of prior scratch calls.
- **Rust/Python boundary:** Rust-only.

### `PreparedSubstructurePattern` (`umol-graph-ir`, bound in `umol-py`)
- **Type and role:** Immutable prepared pattern carrier for 1-pattern × many-hosts and many-patterns
  × 1-host substructure matching.
- **Open carrier or operation-issued value:** Closed validated carrier (`Clone + Debug + Send +
  Sync`).
- **Intrinsic representation invariants:** Holds an `Arc`-shared `Molecule` snapshot with zero
  molecule-scope constraints, its `SubstructureMatchConfig`, consistent deduplicated `AtomForm`
  classes, `u128`/`u8` pre-screen masks, `PreparedSubgraphQuery` (with proper
  stereo+overlay-preserving `QuerySymmetryConstraints` if `reduce_pattern_symmetry`), and prebuilt
  Levi `IncidenceGraph` + query if `match_algorithm == Incidence`.
- **Public constructors:** `new(pattern: &Molecule, config: SubstructureMatchConfig) -> Self` (and
  `#[new]` in Python accepting `(pattern, *, config=None)`).
- **Failure behavior:** Panics on construction if `pattern` carries molecule constraints (matching
  does not evaluate molecule constraints; no error type is encoded into the public API for this
  unimplemented capability). All matching calls (`matches`, `has_match`, `first_match`,
  `count_matches`, `visit_*`, `*_in_host`) are infallible.
- **Semantic properties:** `PreparedSubstructurePattern::new(p, cfg).matches(h)` equals
  `p.substructure_matches(h, cfg)` for every `(p, h, cfg)`. Cloning is $O(1)$ / independent;
  mutating the original `Molecule` after construction does not affect the prepared pattern.
- **Rust/Python boundary:** Immutable (`#[pyclass(frozen)]`) owned object in `umol-py`; exposes
  `pattern()` method returning a detached `Molecule` copy (preserving the `python-api.md`
  property-mutation contract), immutable `config` property, and
  `matches`/`has_match`/`first_match`/`count_matches` accepting either `Molecule` or
  `PreparedSubstructureHost` via `SubstructureHostLike` (justified Python union adaptation avoiding
  duplicated `*_in_host` names in Python).

### `PreparedSubstructureHost` (`umol-graph-ir`, bound in `umol-py`)
- **Type and role:** Immutable prepared host carrier for many-patterns × 1-host and $R \times M$
  matching and reaction application.
- **Open carrier or operation-issued value:** Closed immutable snapshot (`Clone + Debug + Send +
  Sync`) with transparent `OnceLock` lazy caches.
- **Intrinsic representation invariants:** Holds an `Arc`-shared `Molecule` snapshot, matching
  `PreparedSubgraphTarget`, exact `present_elements: u128`, `present_bond_orders: u8`,
  `present_overlays: u8`, per-element & per-bond-order bitmasks, and `OnceLock` slots for
  `IncidenceGraph(Constitution)` (+ Levi `PreparedSubgraphTarget`), `RingSet(Relevant, 22)` per
  `RelevantCycleEnumerationAlgorithm`, and stereo-preserving `TargetSymmetry` per
  `AutomorphismAlgorithm`.
- **Public constructors:** `new(host: &Molecule) -> Self` (infallible; `#[new]` in Python taking
  `molecule: Molecule`).
- **Failure behavior:** Infallible construction; read-only `molecule(&self) -> &Molecule`
  (independent `Molecule` clone in Python); `React::react` / `Reaction::apply_to_host*` preserve
  exact `ApplyPreconditionError` / `ApplyError` channels.
- **Semantic properties:** Matching any pattern or reaction against
  `PreparedSubstructureHost::new(h)` produces the exact same match set, count, existence, and
  reaction applications as matching directly against `h`. Lazy caches are semantically transparent
  (observably pure) and thread-safe. Mutating the original `Molecule` after construction does not
  affect the prepared host.
- **Rust/Python boundary:** Immutable (`#[pyclass(frozen)]`) owned object in `umol-py`; exposes
  `molecule()` method returning a detached `Molecule` copy (preserving the `python-api.md`
  property-mutation contract), `react(reaction, *, config=None)`, and `react_prepared(reaction)`.
  `Reaction.apply*` and `PreparedReaction.apply*` accept `Molecule | PreparedSubstructureHost` via
  `SubstructureHostLike`.

### `PreparedReaction` (`umol-graph-ir`, bound in `umol-py`)
- **Type and role:** Immutable prepared reaction carrier for repeated reaction application ($1
  \times M$ and $R \times M$).
- **Open carrier or operation-issued value:** Closed validated value (`Clone + Debug + Send +
  Sync`).
- **Intrinsic representation invariants:** Holds an `Arc`-shared `Reaction` snapshot, pre-validated
  and pre-normalized `Arc<Deltas>` satisfying every reaction-wide application precondition
  (`ApplyPreconditionError`: valid stereo delta domains, normalized deltas, disjoint bond endpoints,
  rule-local `check_reaction_dpo`), and a precompiled `PreparedSubstructurePattern` for `lhs` under
  `config` (zero molecule-scope constraints on `lhs`; reaction-preserving `QuerySymmetryConstraints`
  via `Reaction::to_reaction_span()` when `config.reduce_pattern_symmetry` is `true`).
- **Public constructors:** `new(reaction: &Reaction, config: SubstructureMatchConfig) ->
  Result<Self, ApplyPreconditionError>` (no `try_` prefix per `nomenclature.md`; `#[new]` in Python
  taking `(reaction, *, config=None)`).
- **Accessors:** `reaction(&self) -> &Reaction`, `config(&self) -> SubstructureMatchConfig`.
- **Failure behavior:** `new` validates reaction-wide preconditions first (returning
  `ApplyPreconditionError`, or `InvalidStructureError` in Python; `ValueError` in Python if live
  `Reaction` components fail integrity) and then compiles `lhs` (panicking if `reaction.lhs()`
  carries molecule constraints, matching `Reaction::apply`). Application methods (`apply`,
  `tracked_apply`, `apply_to_reaction`, `apply_to_reaction_span`, their `*_to_host` counterparts,
  and `React::react_prepared`) are infallible at entry (returning `ReactionApplicationIter<T>` /
  `ReactionProductsIter` directly); per-match execution failures remain `ApplyError` iterator items.
- **Semantic properties:** Applying `PreparedReaction::new(r, cfg)?` against `h` (or
  `PreparedSubstructureHost::new(h)`) produces the exact same match order, products,
  correspondences, reactions, reaction spans, split product sets, and `ApplyError` outcomes as
  `r.apply*(h, cfg)?` / `h.react(&r, cfg)?`. Cloning is $O(1)$; mutating the source `Reaction` or
  host after preparation has no effect.
- **Rust/Python boundary:** Immutable (`#[pyclass(frozen)]`) owned object in `umol-py`; exposes
  `reaction()` method returning a detached `Reaction` copy (preserving the `python-api.md`
  property-mutation contract), immutable `config` property,
  `apply`/`tracked_apply`/`apply_to_reaction`/`apply_to_reaction_span` accepting `Molecule |
  PreparedSubstructureHost` via `SubstructureHostLike`, and `react_prepared` / `react_all_prepared`
  parity on `Molecule` and `PreparedSubstructureHost`.

## Staged implementation plan

Every stage ends **green** (`cargo test --workspace` + clippy `-D warnings`, and `maturin develop` +
`pytest` once `umol-py` is touched).

### Stage S0 — Zero-API-churn foundations & primitives (additive/internal, green)
- **S0a (`umol-graph-core`) — Done:** Cache `max_degree: u32` in `Csr` and expose
  `Graph::max_degree(&self) -> usize`; replace `BTreeSet` validation in `Correspondence::new` and
  `Correspondence::from_images` with inline zero-heap validation for `len <= 32`, empty fast-return,
  and skip-sort-if-sorted. Unit + property tests. **Additive (green).** `[dep: none]`
  - *Implementation notes:* `Csr::max_degree` is populated in `try_build_csr` and needs no update in
    `extend_nodes` (appending isolated nodes leaves `max_degree` unchanged). In
    `Correspondence::new`'s inline `len <= 32` scan, recording `duplicate_right` across `0..i`
    rather than returning immediately on the first right duplicate preserves exact
    `DuplicateLeftId`-before-`DuplicateRightId` error precedence (fixing a subtle precedence
    divergence present in the scratch prototype). Verified with `umol-graph-core` unit tests and
    `proptest` (`graph::test_graph_max_degree`, `correspondence::test_correspondence_new`,
    `correspondence::test_correspondence_new_bijection`, and widened `0..48`
    `test_correspondence_from_images`).
- **S0b (`umol-graph-ir`) — Done:** Rewrite `induced_bonds` to use `right.raw_graph().find_edge`
  directly (no `BTreeMap`) and short-circuit empty pattern overlay kinds in `verify_overlays`. Unit
  tests covering overlay-free and overlay-bearing parity. **Additive (green).** `[dep: S0a]`
  - *Implementation notes:* `induced_bonds` validates `atoms` count compatibility (preventing
    out-of-bounds `find_edge` panics if called directly with mismatched counts) and indexes
    `atoms.matched_pairs()` directly in $O(1)$ when `atoms.is_total_on_left()` instead of
    binary-searching `right_of` per bond endpoint. `induce_by_key` also gained an `left_count == 0
    || right_count == 0` fast-return, and `verify_overlays` short-circuits both empty pattern
    overlay kinds (preserving exact `right_count`) and insufficient host overlay counts
    (`host.<kind>().count() < pattern.<kind>().count()`). Verified green across the full workspace
    (`umol-graph-ir` `correspondence`/`substructure` proptests + `cargo test --workspace` + `clippy
    --workspace --all-targets -- -D warnings`).

### Stage S1 — Core predicate trait, symmetry constraints, scratch, and prepared query/target (additive, green)
- **S1a (`umol-graph-core`) — Done:** Define `SubgraphIsomorphismPredicate` trait (using
  `undetermined` edge hooks, query node equivalence classes, bulk `u64`/`u128` class candidate fill,
  optional edge-neighbor bitmask hook) and internal `pub(crate)`
  `ClosureSubgraphIsomorphismPredicate` adapter. Unit tests. **Additive (green).** `[dep: S0a]`
  - *Implementation notes & public-surface tightening:* Placed in
    `algorithms/subgraph_isomorphism/predicate.rs` (establishing the `subgraph_isomorphism/`
    submodule layout for `S1`–`S2`). Bulk candidate fills (`node_candidates_u64`,
    `node_candidates_u128`, `fill_node_candidates_words`) take an `eligible` bitmask so callers can
    intersect precomputed masks (e.g., `deg_ge[deg(q)]` & element masks) directly and default
    fallbacks only probe set bits; `edge_neighbors_u64` / `edge_neighbors_u128` return `Option`
    exact neighbor bitmasks so bond-order bitsets bypass `find_edge` + `match_edge` and preserve
    leaf `POPCNT`. After reviewing all workspace callers (`umol-graph-ir` `substructure.rs`
    implements the trait directly; `GraphView` only proxies `Graph`'s closure methods as an escape
    hatch and can wire additional methods additively later), tightened the public surface by making
    `ClosureSubgraphIsomorphismPredicate` `pub(crate)` (stripped down to the minimal `(node_match,
    edge_match)` bridge) and dropping `_with_scratch` closure overloads on `Graph` (keeping
    `SubgraphIsomorphismScratch` strictly on `PreparedSubgraphQuery`). Verified with 12 unit tests
    and `clippy --all-targets --features proptest -- -D warnings`.
- **S1b (`umol-graph-core`) — Done:** Implement `QuerySymmetryConstraints` (Grochow–Kellis
  stabilizer chain with deterministic `(Reverse(deg), q)` base; `none`, `from_graph`,
  `from_graph_with_scratch`, `try_from_generators`) and `TargetSymmetry` (`none`, `from_graph`,
  `from_graph_with_scratch`, `from_automorphism_group`, `try_from_generators`, shallow
  point-stabilizer orbit helper, combined $G_Q \times G_H$ exact leaf canonicality guard).
  Exhaustive differential tests on `simple-through-8.g6` verifying exact orbit counts and `TrailIr`
  == `Nauty` invariance. **Additive (green).** `[dep: S0a]`
  - *Implementation notes & findings:* Added `algorithms/subgraph_isomorphism/symmetry.rs` exporting
    `QuerySymmetryConstraints`, `TargetSymmetry`, and shared `SymmetryGeneratorError`
    (`LengthMismatch`, `NodeOutOfBounds`, `DuplicateImage`, `EdgeNotPreserved`). Added
    `TargetSymmetry::none` and `TargetSymmetry::try_from_generators` so `S4b`
    `PreparedSubstructureHost` can pass proper stereo-preserving host subgroups (which cannot be
    built as `AutomorphismGroup`). Both carriers compute exact point stabilizers via one-pass
    Schreier's Lemma + sparse Sims filter + exact disjoint-support commutation pruning (`Support`),
    requiring only a single `AutomorphismAlgorithm` call per graph and $O(n)$ Schreier-tree storage
    per active level. Aligning `TargetStabilizerChain` push-down and
    `is_canonical_double_orbit_embedding` onto the exact same `(Reverse(deg), q)`
    `query_priority_order` as `QuerySymmetryConstraints` prevents cross-base orbit loss (verified
    against the 4-node `H = {0-2, 1-3}, Q = {0-1}` counterexample) and guarantees identical
    canonical representative sets across all four `SubgraphIsomorphismAlgorithm`s and both
    `AutomorphismAlgorithm`s. Reusing the stabilizer Schreier trees avoids BFS during combined leaf
    recursion, and `suffix_is_global_min[r]` prevents factorial branching on symmetric cliques
    (`K_10` in `K_12`). Deliberately omitted `PartialEq` on both carriers so equality does not
    conflate generator bases or distinct subgroups with matching constraints. Verified in 0.76 s
    (debug): 13,598/13,598 `simple-through-8.g6` graphs (`TrailIr` == `Nauty`), 52×52 pairs through
    order 5 (exact orbit counts and unique lex-min representatives across $G_Q$, $G_H$, and $G_Q
    \times G_H$), $K_{10} \hookrightarrow K_{12}$, 12 disjoint triangles ($|G| = 6^{12}$), wide $n =
    144$ bitmasks, and clean `clippy --all-targets --features proptest -- -D warnings`.
- **S1c (`umol-graph-core`) — Done:** Implement `SubgraphIsomorphismScratch`,
  `PreparedSubgraphTarget` (`t_hist`, `deg_ge`, `adj` `u64`/`u128`/wide), and
  `PreparedSubgraphQuery` (alg-tagged plan enum + `QuerySymmetryConstraints`). Unit tests.
  **Additive (green).** `[dep: S1a, S1b]`
  - *Implementation notes & findings:* Added `algorithms/subgraph_isomorphism/prepared.rs` exporting
    `SubgraphIsomorphismScratch`, `PreparedSubgraphTarget`, and `PreparedSubgraphQuery`.
    `PreparedSubgraphTarget` stores suffix-sum cumulative `degree_hist` (`hist[d]` = vertices with
    `deg >= d`) and width-specialized `TargetBitsets` (`U64` for `<= 64`, `U128` for `65..=128`,
    `Wide` for `> 128`), eagerly materializing wide adjacency matrices up to `WIDE_EAGER_ADJ_LIMIT =
    1024` and falling back to lazy per-row scratch adjacency on larger sparse targets to avoid
    quadratic memory allocation. `PreparedSubgraphQuery` precomputes flat CSR symmetry-bound tables
    (`sym_lower_bounds` / `sym_upper_bounds` for $O(1)$ bitmask/scalar enforcement),
    `StaticOrderPlan` for `Vf2Rdkit`, and `RkStep` (`Seed`/`Extend`/`Close`) for `RayKirsch`
    starting at `priority_order[0]` (`(Reverse(deg), q)`) so `RayKirsch` gets immediate depth-0 host
    orbit-rep pruning. `try_with_symmetry -> Option<Self>` (paired with asserted `with_symmetry`)
    validates both `node_count` and `priority_order` compatibility for non-trivial symmetry without
    introducing a single-variant error type. Verified with 9 unit tests (including snapshot mutation
    isolation, `n = 0/4/70/144/1050` bitset widths, symmetry CSR bounds, plan invariants, and
    mismatch rejection) plus full Stage S1 workspace gate (`cargo test --workspace`, `cargo clippy
    --workspace --all-targets -- -D warnings`, and `clippy -p umol-graph-core --all-targets
    --features proptest -- -D warnings`).

### Stage S2 — Core engines (`BitsetRi`, `Glasgow`, upgraded `RayKirsch`, cleaned `Vf2Rdkit`) + mode splits (additive, green)
- **S2a (`umol-graph-core`) — Done:** Add `BitsetRi` variant and implementation (`u64`/`u128`/wide
  bitsets, degree-sequence filter, descending-degree class candidate fill, RI-DS rarity+connectivity
  order, lazy per-row `adj`, direct sparse single-parent scan, 2-instruction bit-parallel
  1-lookahead, push-down $G_Q$/$G_H$ symmetry, anchored `_at`, leaf `POPCNT` for `count`).
  **Additive (green).** `[dep: S1c]`
  - *Implementation notes & findings:* Added `algorithms/subgraph_isomorphism/bitset_ri.rs` and
    `SubgraphIsomorphismAlgorithm::BitsetRi` (plus exhaustive match arms in
    `umol-py/src/algorithm.rs` and the two benchmark `subiso_name` helpers to keep the workspace
    green). Monomorphized the `u64`/`u128` core over `WordAdj` (`PrecomputedWordAdj` vs.
    `LazyWordAdj`), `LOOKAHEAD`, and `COUNT_ONLY`, keeping direct `Graph` calls zero-heap for `n2 <=
    128` and `m1 <= 128` (`STACK_PARENTS_LIMIT = 128`) and deduplicating `query_node_class_count()`
    on stack buffers whenever `num_classes <= N` so prepared callers incur zero per-match heap
    allocations. In `Wide` mode (`n2 > 128`), capped lazy dense scratch adjacency at
    `WIDE_LAZY_ADJ_LIMIT = 4096` (falling back to `Graph::find_edge` on larger sparse targets so
    memory never grows quadratically). Added `is_canonical_double_orbit_with_built_chain` in
    `symmetry.rs` so combined $G_Q \times G_H$ leaves reuse the active push-down
    `TargetStabilizerChain` without rebuilding stabilizers; static `prefix_len_after[d]` handles
    arbitrary RI-DS vertex orders even when rare vertex labels force `order[0] !=
    priority_order[0]`, skips building the final leaf Schreier tree in host-only mode (`!has_q_sym
    && r + 1 == n1`), and enables leaf `POPCNT` counting whenever all leaf parent edges are
    bitmask-resolved (`edge_is_undetermined` or `edge_neighbors`) and no active $G_Q \times G_H$
    check or non-trivial host stabilizer remains. Also skipped lookahead `adj.row(v)` at `depth ==
    n1 - 1` (where `q_term == q_rem == 0`), resolved a `prepared.rs` test `&[]` type-inference
    ambiguity with `serde_json` (`&[0u32; 0]`), and removed now-used
    `#[allow(dead_code)]`/`#[allow(unused_imports)]` attributes in `predicate.rs`, `symmetry.rs`,
    and `subgraph_isomorphism.rs`. Verified across `u64`/`u128`/wide (`144`, `1050`, `4200`),
    high-degree `qd >= 15` exact filtering, anchored `_at`, predicate hooks (class deduplication +
    `edge_neighbors` + zero `match_edge` on undetermined/masked edges), scratch reuse after early
    `Break`, non-priority start-order symmetry, and all 52×52 small graph pairs × 4 `(G_Q, G_H)`
    modes (`TrailIr` == `Nauty`).
- **S2b (`umol-graph-core`) — Done:** Add `Glasgow` variant and implementation (compact $n_1$-word
  flat trail `u64`/`u128`/wide, bit-parallel FC, dynamic MRV, unit propagation, Hall/AllDifferent
  wipeout pruning, push-down $G_Q$/$G_H$ symmetry, anchored `_at`, leaf `POPCNT` for `count`).
  **Additive (green).** `[dep: S1c]`
  - *Implementation notes & findings:* Extracted shared private bitset primitives (`WordMask`,
    `WordAdj`, `WideAdj`, `degree_sequence_ok`, `init_word_candidates`, `init_wide_candidates`) into
    `algorithms/subgraph_isomorphism/bitset.rs` to eliminate duplication between `BitsetRi` and
    `Glasgow`. Added `algorithms/subgraph_isomorphism/glasgow.rs` and
    `SubgraphIsomorphismAlgorithm::Glasgow` (plus exhaustive match arms in
    `umol-py/src/algorithm.rs` and the two benchmark `subiso_name` helpers). In `u64`/`u128` mode,
    `Glasgow` uses a compact `(n1 + 1) * n1`-word flat trail (stack-allocated for `n1 <= 16` in
    direct calls; reused `scratch.trail_u64`/`trail_u128` in prepared calls), tracks verified
    forward-checked edges (`edge_is_undetermined` or `edge_neighbors == Some`) in a value-passed
    `verified_edges: u128` bitmask so verified edges skip `find_edge`/`match_edge` and unlock leaf
    `POPCNT`, propagates forced singletons (`unit_done` bitmask passed by value so each singleton
    propagates at most once along any root-to-leaf path), and runs the Glasgow bit-parallel
    Hall/AllDifferent sweep gated by an exact `small_count >= 2` theorem (`|dom[q]| < k_nu`) so
    unconstrained branches skip sorting completely while Hall-discovered singletons loop back into
    unit propagation. In `Wide` mode (`n2 > 128`), used an $O(n_1 \cdot \text{words})$ domain bitset
    with live `dom_cnt` (`O(n1)` dynamic MRV and `O(1)` leaf `POPCNT`) and a LIFO sparse word-level
    undo trail instead of copying `(n1 + 1) * n1 * words`, keeping memory linear in $n_1 \cdot n_2 /
    64$ and falling back to per-row `target.neighbors` materialization above `WIDE_LAZY_ADJ_LIMIT =
    4096`. Push-down `G_Q` Grochow–Kellis inequalities narrow peer domains directly during FC and
    unit propagation; dynamic `G_H` tracks `prefix_len` along `priority_order`, filters
    `dom[p_{next_prefix_len}]` to current stabilizer orbit minima immediately after pushing a
    stabilizer level, and reuses the built `TargetStabilizerChain` at leaves for combined $G_Q
    \times G_H$ checks. Verified across `u64`/`u128`/`Wide` (`144`, `1050`, `4200`), Hall-to-unit
    cascades, predicate hooks, anchored `_at`, dynamic non-priority MRV symmetry, scratch reuse, and
    all 52×52 small graph pairs × 4 symmetry modes (`TrailIr` == `Nauty`).
- **S2c (`umol-graph-core`) — Done:** Upgrade `RayKirsch` (`RkStep` RI-guided plan + `deg(q) <=
  deg(t)` + symmetry constraints) and clean up `Vf2Rdkit` (static BFS order + flat parent CSR +
  direct parent edge + symmetry constraints); expose `has_subgraph_isomorphism`,
  `first_subgraph_isomorphism`, `count_subgraph_isomorphisms` (+ `_at`) on `Graph` and
  scratch-taking `visit_matches`/`has_match`/`first_match`/`count_matches` on prepared query/target.
  **Additive (green).** `[dep: S1c]`
  - *Implementation notes & findings:* Extracted cleaned `Vf2Rdkit`
    (`algorithms/subgraph_isomorphism/vf2_rdkit.rs`) and upgraded `RayKirsch`
    (`algorithms/subgraph_isomorphism/ray_kirsch.rs`) into dedicated submodules matching `BitsetRi`
    and `Glasgow`'s `visit_direct`/`count_direct`/`visit_prepared`/`count_prepared` shape and
    removing the legacy inline implementations from `subgraph_isomorphism.rs`. `Vf2Rdkit` consumes
    `PreparedQueryPlan::Vf2Rdkit(StaticOrderPlan)` directly, checks the primary parent edge from CSR
    `nb.edge` without `find_edge`, enforces `deg(q) <= deg(t)`, and pushes down static `G_Q`/`G_H`
    symmetry; `RayKirsch` executes `PreparedQueryPlan::RayKirsch(Vec<RkStep>)` directly with `deg(q)
    <= deg(t)`, `edge_is_undetermined` skipping, and step-indexed `G_Q`/`G_H` push-down (handling
    `Close` steps cleanly without redundant stabilizer work). Exposed
    `Graph::{has,first,count}_subgraph_isomorphism{,_at}` (returning `bool`,
    `Option<Correspondence<NodeId>>`, and `u64`) and
    `PreparedSubgraphQuery::{visit_matches,has_match,first_match,count_matches}` (validating
    `target_symmetry.node_count() == target.node_count()` up front and temporarily supporting legacy
    selectors via `RefCell` + leaf canonicality until `S3a` retires them). Trimmed 5 unused fields
    (`cand_u64`, `cand_u128`, `class_u64`, `class_u128`, `trail_wide`) from opaque
    `SubgraphIsomorphismScratch` and unused `priority_rank` from `PreparedSubgraphQuery`, leaving
    zero `#[allow(dead_code)]` attributes across `subgraph_isomorphism/`. Verified with 200
    `subgraph_isomorphism` unit tests (including cross-engine parity across `[Vf2Rdkit, RayKirsch,
    BitsetRi, Glasgow]` × all 4 symmetry modes and all `Graph`/`PreparedSubgraphQuery` mode splits)
    plus full workspace clippy and test pass.
- **S2d (`umol-graph-core`) — Done:** Exhaustive differential & property test suite verifying exact
  embedding-set parity across `[Vf2Rdkit, RayKirsch, BitsetRi, Glasgow]`, both `[TrailIr, Nauty]`,
  all 4 `(pattern_sym, host_sym)` combinations, `visit`/`has`/`first`/`count`, anchored `_at`, and
  $n \le 64$, $65..=128$, $n > 128$. **Additive (green).** `[dep: S2a, S2b, S2c]`
  - *Implementation notes & findings:* Added independent definition-level monomorphism, colored
    automorphism-group, subgroup-closure, and canonical lex-min orbit-representative calculations in
    `tests/property/subgraph_isomorphism/exhaustive.rs`, and expanded
    `tests/property/subgraph_isomorphism.rs` into 27 differential, exhaustive, boundary, and
    generative tests: (1) direct `Graph` proptests over random, non-induced planted, and anchored
    `_at` labeled graphs comparing all 8 currently-present selectors across
    `enumerate`/`visit`/`has`/`first`/`count`/early `Break` directly against
    `exhaustive_monomorphisms` (removing VF2/Ullmann as ground truth ahead of `S3a`); (2) random and
    non-induced planted `U64` proptests exercising overlapping lattice-style bitmask predicates
    (`PointwiseMaskPredicate` vs. `BitsetMaskPredicate` with `all_edges_undetermined`,
    `edge_is_undetermined`, `query_node_class_count`, bulk
    `node_candidates_{u64,u128}`/`fill_node_candidates_words`, and `edge_neighbors_{u64,u128}` in
    `AllSome`/`Mixed`/`AllNone` modes), `TrailIr` == `Nauty`, all 4 `(pattern_sym, host_sym)` modes,
    and all 4 core engines against independent `canonical_orbit_representatives`; (3) exhaustive
    `simple-through-8.g6` corpus checks over all 52 × 52 pairs through order 5 (2,704 pairs) and all
    18 × 156 order-1..=4 queries in order-6 targets (2,808 pairs) across all 4 core engines × 4
    symmetry modes × `[TrailIr, Nauty]`; (4) `U128` (`65..=128`) and `Wide` (`129..=196`) generative
    proptests over sparse labeled graphs with repeated symmetric motifs, shuffled vertex IDs,
    non-priority rare-leaf branching motifs (forcing `BitsetRi` RI-DS and `Glasgow` MRV to branch
    outside `priority_order[0]` under active `G_Q`/`G_H` symmetry), UNSAT cases, direct/prepared
    parity, anchored `_at`, and independent `is_valid_monomorphism` validation; (5) explicit
    `#[rstest]` coverage for empty graphs (`0/0`, `0/8`), exact word boundaries (`64`, `65`, `128`,
    `129`), `n1 = 18` (> `Glasgow` stack-trail cap) across all three widths, high-degree wheels
    (`W_17` with `W_16` decoy at `64`/`128`/`140`), sparse wide targets above both
    `WIDE_EAGER_ADJ_LIMIT` (`1050`) and `WIDE_LAZY_ADJ_LIMIT` (`4120`), proper subgroups via
    `try_from_generators` (`C_4`, `A_4`, disconnected pair), and `SubgraphIsomorphismScratch` reuse
    across regimes after forced early `Break`. All 231 tests in `umol-graph-core --test property
    --features proptest` passed (0 failures, 0 clippy warnings) with zero defects uncovered in
    `Vf2Rdkit`, `RayKirsch`, `BitsetRi`, `Glasgow`, or symmetry reduction.

### Stage S3 — Retire dominated `SubgraphIsomorphismAlgorithm` variants & migrate workspace (breaking red→green)
- **S3a (`umol-graph-core`) — Done:** Remove `Vf2`, `Ullmann`, `Ri`, `ArcMatch { path_length }`, and
  `ARCMATCH_DEFAULT_PATH_LENGTH` from `SubgraphIsomorphismAlgorithm` and `lib.rs`. **Breaking (red
  until S3b).** `[dep: S2d]`
  - *Implementation notes & findings:* Removed `Vf2`, `Ullmann`, `Ri`, `ArcMatch { path_length }`,
    and `ARCMATCH_DEFAULT_PATH_LENGTH` from `SubgraphIsomorphismAlgorithm` and
    `umol-graph-core/src/lib.rs`, deleting ~1,800 lines of retired VF2/Ullmann/RI/ArcMatch
    implementation and ArcMatch-specific tests from `subgraph_isomorphism.rs`. Removed
    `PreparedSubgraphQuery::visit_matches_legacy`, `RefCell`, and wildcard dispatch arms in
    `prepared.rs`, and narrowed `is_canonical_double_orbit_embedding` in `symmetry.rs` to
    `#[cfg(test)]`. Updated `umol-graph-core/tests/property/subgraph_isomorphism.rs` and the
    selector list in `umol-graph-core/benches/subgraph_isomorphism.rs` so `umol-graph-core` compiles
    and passes `--all-targets --features proptest` cleanly (0 clippy warnings, 2,172 unit tests +
    231 property tests + doc tests passing; downstream crates intentionally red until `S3b`).
- **S3b (`umol-graph-ir`, `umol-graph`, `umol-py`, benches, tests, fuzz targets, `README.md`) —
  Done:** Migrate all callers, tests, fuzz targets (`fuzz_reaction.rs`), `umol-py`
  `SubgraphIsomorphismAlgorithm` bindings, and `umol-graph-core/benches/subgraph_isomorphism.rs`
  (expanding benchmark coverage to include UNSAT/hard/large/symmetry-reduced cases across the 4
  retained algorithms). Workspace green. **Migration (green).** `[dep: S3a]`
  - *Implementation notes & findings:* Migrated all remaining references to `Vf2`, `Ullmann`, `Ri`,
    `ArcMatch`, and `ARCMATCH_DEFAULT_PATH_LENGTH` across `umol-graph-ir` (`compose.rs`,
    `molecule/tests.rs`, `reaction.rs`, `substructure.rs`, `tests/integrity.rs`,
    `tests/property/{substructure,reaction/*}.rs`, `benches/reaction.rs`,
    `fuzz/fuzz_targets/fuzz_reaction.rs`), `umol-graph` (`src/fingerprint/pattern.rs`,
    `benches/substructure.rs`), `umol-py`
    (`src/{algorithm,fingerprint/config,molecule,reaction,substructure}.rs` and
    `tests/test_{algorithm,fingerprint,reaction,substructure,workflow}.py`), and
    `docs/development/property-tests.md` (`README.md` already used
    `SubgraphIsomorphismAlgorithm.RayKirsch()` and needed no edits). Expanded
    `umol-graph-core/benches/subgraph_isomorphism.rs` into four Criterion groups across `[Vf2Rdkit,
    RayKirsch, BitsetRi, Glasgow]`: core labeled matrix (using `has_subgraph_isomorphism` and
    `count_subgraph_isomorphisms` for `first`/`count`), `unsat_hard` (`petersen/square`,
    `dodecahedron/square`, `c60/square`, `c60/k4`, `c70/square`, `grid_8x8/triangle`,
    `k10_minus_m5/k6`), `large` (`c70/hexagon`, `grid_12x12/{square,grid_3x3}`,
    `grid_16x16/grid_3x3` via prepared carriers + scratch), and `symmetry_reduced`
    (`none`/`pattern`/`host`/`both` on `cubane/square`, `dodecahedron/pentagon`, `c60/hexagon`,
    `c70/hexagon`). Verified zero remaining retired selector occurrences outside
    `discussion/`/`scratch/`/`materials/`, clean `cargo clippy --workspace --all-targets -- -D
    warnings`, clean `cargo clippy -p umol-graph-core -p umol-graph-ir --all-targets --features
    proptest -- -D warnings`, full `cargo test --workspace`, `umol-graph-ir --test property
    --features proptest substructure` (34 passed), `maturin develop` + `pytest -q umol-py/tests`
    (3,654 passed, 2 skipped under Python 3.13.18), and `cargo bench -p umol-graph-core --bench
    subgraph_isomorphism -- --test`.

### Stage S4 — `umol-graph-ir` config symmetry fields, prepared carriers, and mode splits (breaking config + migration)
- **S4a (`umol-graph-ir`, `umol-graph`, `umol-py`, benches, tests, fuzz) — Done:** Add
  `automorphism_algorithm: AutomorphismAlgorithm`, `reduce_pattern_symmetry: bool`, and
  `reduce_host_symmetry: bool` to `SubstructureMatchConfig`; migrate all struct literals across the
  workspace. **Breaking → green.** `[dep: S3b]`
  - *Implementation notes & findings:* Added `automorphism_algorithm: AutomorphismAlgorithm`,
    `reduce_pattern_symmetry: bool`, and `reduce_host_symmetry: bool` to `SubstructureMatchConfig`
    in `umol-graph-ir/src/ir/substructure.rs` (keeping no `Default` impl at the graph-IR layer so
    all algorithm and symmetry choices remain explicit). Migrated all workspace struct literals
    across `umol-graph-ir` (`substructure.rs`, `reaction.rs` including the `React` doctest,
    `compose.rs`, `benches/reaction.rs`, `fuzz/fuzz_targets/fuzz_reaction.rs`, `tests/integrity.rs`,
    `tests/property/{substructure,reaction/*}.rs`), `umol-graph` (`src/fingerprint/pattern.rs`
    keeping both reduction flags `false` for RDKit `uniquify = false` parity,
    `benches/substructure.rs`), and `umol-py` (`src/{substructure,reaction,molecule}.rs`, keeping
    the Python config surface unchanged until `S7a` exposes the new options once `S4c` and `S5a`
    wire them). Verified green with `cargo +nightly fmt --all`, `git diff --check`, `cargo clippy
    --workspace --all-targets -- -D warnings`, `cargo clippy -p umol-graph-ir --all-targets
    --features proptest -- -D warnings`, `cargo test --workspace`, and `maturin develop` + `pytest
    -q umol-py/tests` (3,654 passed, 2 skipped under Python 3.13.18).
- **S4b (`umol-graph-ir`) — Done:** Implement `PreparedSubstructureHost` (immutable snapshot,
  `OnceLock` lazy `IncidenceGraph`/`RingSet`/stereo-preserving `TargetSymmetry`, `u128`/`u8` masks,
  per-element & per-bond-order bitsets) and `PreparedSubstructurePattern` (deduplicated `AtomForm`
  classes, `undetermined` bond masks, `u128`/`u8` pre-screen, proper stereo+overlay-preserving
  `QuerySymmetryConstraints`, Levi prep). **Additive (green).** `[dep: S4a]`
  - *Implementation notes & findings:* Added generic canonical bottom-up Schreier–Sims subgroup
    extraction (`extract_subgroup_generators`, `project_generators`) in
    `umol-graph-ir/src/ir/symmetry.rs` (reusable by `ReactionSpan` in `S5a`), reconciled
    molecule-level proper (`Orientation::Proper`) generator extraction with `GraphSymmetry` in
    `umol-graph-ir/src/ir/molecule/symmetry.rs` (`pub(crate)`
    `Molecule::proper_incidence_generators(&self, &GraphSymmetry)` over
    `IncidenceLevel::Constitution` via `remapping_matches_stereo`, plus
    `Molecule::proper_atom_generators(&self, algorithm, include_constraints)` taking a direct
    `raw_graph().automorphisms(...)` fast path for overlay-free, stereo-free molecules and otherwise
    projecting `proper_incidence_generators`), and implemented `PreparedSubstructureHost` and
    `PreparedSubstructurePattern` in `umol-graph-ir/src/ir/substructure/prepared.rs` (re-exported
    from `ir::substructure` and `ir`). Canonical coset-image sorting (`g_cur(x)`) and `Stab(base)`
    budget fallback (`SUBGROUP_SEARCH_NODE_BUDGET = 8192`) preserve strict `TrailIr` == `Nauty`
    invariance even on non-trivial proper stereo subgroups (e.g., tetrahedral `CH4` `Th0` yielding
    exact order-12 `A4 < S4`). `PreparedSubstructureHost` stores compact bit-ranked per-element atom
    bitmasks (`u64`/`u128`/`Wide`), bit-ranked per-bond-order neighbor bitsets (`u64`/`u128`),
    `u128`/`u8`/overlay-count screening summaries with explicit literal/satisfiability guards for
    symbolic or unsatisfiable raw host forms, and `Arc`-shared `OnceLock` lazy caches for
    `IncidenceGraph(Constitution)` (+ precomputed per-edge electron contributions),
    `RingSet(Relevant, 22)` per `RelevantCycleEnumerationAlgorithm`, and proper atom/incidence
    `TargetSymmetry` per `AutomorphismAlgorithm`. `PreparedSubstructurePattern` deduplicates
    atom/bond predicate classes (`Any`/`LitElementOnly`/`ElementSetOnly`/`Filtered`,
    `Undetermined`/`OrderMask`/`General`), reuses thread-local `SubgraphIsomorphismScratch` +
    atom-slice buffers with re-entrancy fallback, and calls
    `PreparedSubgraphQuery::{has_match,count_matches}` directly (unlocking leaf `POPCNT`) whenever
    `present_overlays == 0`. Verified with `cargo +nightly fmt --all`, `git diff --check`, `cargo
    clippy -p umol-graph-ir --lib --tests -- -D warnings`, `cargo test -p umol-graph-ir --lib
    substructure` (94 passed), and `cargo test -p umol-graph-ir --lib
    test_molecule_proper_generators` (13 passed, including direct agreement with
    `Molecule::proper_orbit_of`).
- **S4c (`umol-graph-ir`) — Done:** Add `Molecule::{visit_substructure_embeddings,
  has_substructure_match, first_substructure_match, count_substructure_matches}` and `GraphView`
  mode splits; wire `Molecule::visit_substructure_matches` over the fast prepared/predicate engine;
  add comprehensive unit & property tests (`GraphAndOverlays` == `Incidence`, stereo-preserving
  symmetry, all 4 symmetry modes, prepared vs direct parity, all 4 algorithms). **Additive
  (green).** `[dep: S4b]`
  - *Implementation notes & findings:* Rewired `Molecule::visit_substructure_matches` over
    `PreparedSubstructurePattern` (retiring the legacy per-run `host_ring_context`,
    `constraint_tables`, `visit_substructure_matches_graph_and_overlays`, and
    `visit_substructure_matches_incidence` helpers), removed `SubstructureMatchError` and
    `ApplyPreconditionError::Match` completely (switching `PreparedSubstructurePattern::new` and all
    `Molecule::*substructure_*` methods to infallible signatures that panic upfront when a pattern
    carries molecule constraints rather than encoding an unimplemented capability into public
    `Result` types, including `umol-py` `Molecule.substructure_matches`), added
    `Molecule::{visit_substructure_embeddings, has_substructure_match, first_substructure_match,
    count_substructure_matches}` and `GraphView::{has,first,count}_subgraph_isomorphism{,_at}`
    (`first_*` returning `Option<Vec<AtomId>>` matching
    `GraphView::enumerate_subgraph_isomorphisms`'s `Vec<Vec<AtomId>>` element type), and expanded
    `umol-graph-ir/tests/property/substructure.rs` with stereo-symmetry fixtures
    (`test_substructure_stereo_symmetry`), full mode-split and prepared-carrier parity in
    `test_substructure_cross_validation{,_planted}`, symmetry-mode cross-validation
    (`test_substructure_symmetry_modes_planted`), and independent exhaustive `n!`
    proper-automorphism orbit validation (`test_substructure_symmetry_exact_orbits`).
- **S4d (`umol-graph-ir`) — Done:** Make `IncidenceGraph`, `Incidence`, `IncidenceLevel`,
  `StereoVirtualLigands`, and `{Molecule, ReactionSpan}::incidence_graph` `pub(crate)`; keep
  localized `Bond`s as direct edges (`0..bond_count`) rather than pseudonodes and drop
  `Incidence::BondEndpoint` and `Incidence::StereoSite` (stereo bonds attach via `StereoAnchor` on
  both endpoint atoms); address `IncidenceGraph` through its own `#[repr(transparent)]` domain IDs
  (`IncidenceNodeId`, `IncidenceEdgeId`, `IncidenceEdge::{Bond, Incidence}`) with `raw_graph() ->
  &Graph` and zero-cost `From` conversions at the `umol-graph-core` boundary; update
  `molecule/symmetry.rs`, `molecule/canonicalize.rs`, `reaction_span/symmetry.rs`,
  `reaction_span/canonicalize.rs` (retiring `TopologySearchInput`/`topology_search_input`),
  `substructure/prepared.rs` (enabling `edge_is_undetermined`, `all_edges_undetermined`, and
  `edge_neighbors_{u64,u128}` on localized bonds in `SubstructureMatchAlgorithm::Incidence`), tests,
  benches, and `docs/development/nomenclature.md`. **Internal cleanup + perf (green).** `[dep: S4c]`
  - *Implementation notes & findings:* Narrowed `IncidenceGraph`, `Incidence`, `IncidenceLevel`,
    `StereoVirtualLigands`, and `{Molecule, ReactionSpan}::incidence_graph` to `pub(crate)`
    (dropping the re-exports from `ir.rs` and introducing crate-private `#[repr(transparent)]`
    `IncidenceNodeId`/`IncidenceEdgeId` and `IncidenceEdge::{Bond(BondId), Incidence(Incidence)}`
    with `raw_graph() -> &Graph`, mirroring `Molecule`'s `AtomId`/`BondId`/`raw_graph()` boundary).
    Localized bonds remain direct atom–atom edges (`0..bond_count`) with an `O(1)`
    `raw_graph().clone()` fast path when no overlays are included, and `incidences: Vec<Incidence>`
    stores only overlay–atom incidences (zero bytes on overlay-free molecules). Dropped
    `Incidence::BondEndpoint`, `Incidence::StereoSite`, `TopologySearchInput`,
    `topology_search_input`, `OrderedPartition::from_color_ranks`, and
    `canonical_search_from_partition`; added `OrderedPartition::individualize_pair` for fixing a
    `StereoBond`'s unordered endpoint pair during `stereogenic_site_kinds`; upgraded
    `IncidenceSubstructurePredicate` in `substructure/prepared.rs` to match localized bonds directly
    with `all_edges_undetermined`, `edge_is_undetermined`, and `edge_neighbors_{u64,u128}` (widening
    `HostBondOrderBitsets::U64` in `u128` for `atom_count <= 64 < incidence_node_count <= 128`,
    verified by `test_prepared_substructure_pattern_incidence_u64_atoms_u128_nodes`); and fixed a
    pre-existing bug in `visit_remapping_candidates` (`ir/symmetry.rs`) that checked `same_cells`
    only after refinement (allowing incompatible initial individualizations to collapse to identical
    discrete partitions). Updated tests, `benches/canonicalize.rs`, `tests/integrity.rs`,
    `substructure.rs` rustdoc, and `docs/development/nomenclature.md`. Verified green with `cargo
    +nightly fmt --all`, `git diff --check`, `cargo clippy -p umol-graph-ir --all-targets --features
    proptest -- -D warnings`, `cargo check -p umol-graph -p umol-io -p umol-geometric-graph`,
    `umol-graph-ir` unit tests (`incidence`: 136, `symmetry`: 697, `canonicalize`: 501,
    `substructure`: 111), `--test integrity` (19), and `--test property --features proptest`
    (`substructure`: 42, `canonicalize`: 25, `symmetry`: 10, `stereo`: 158).
  - *Before (`ddfa0c61f`) vs after (`S4d`) benchmarks:* Pure `S4d` improved 43/43 measured
    `canonicalize` cases with 0 regressions (topology-only 1.01–1.08×; `aromatic_electron_counts`
    1.21×; `para_stereo_cascade` 1.31× mol / 1.36× span; `symmetric_stereo_triisopropanolamine`
    1.36× / 1.32×; `symmetric_stereo_butanol` 1.84× / 1.78×; `symmetric_stereo_dimethylbutanol`
    2.40× / 2.35×) and 6/6 `constraint_readings` cases (`chfclbr/stereogenic_site_kinds` 1.33×,
    `chfclbr/prepare_and_stereogenic_site_kinds` 1.83×, `stereogenic_site_kinds/*` 1.60–1.92×).
    Combined `S4c`+`S4d` one-shot `umol-graph/benches/substructure.rs` (`BitsetRi`, rebuilding
    prepared pattern/host per call) improved OpenSMILES ~9k `Incidence` by 1.87–3.48× (`branched`
    44.4→23.7 ms, `phenol` 59.2→24.4 ms, `bicyclic` 89.9→25.9 ms — within 15–23% of
    `GraphAndOverlays`) and `aromatic_ring*/incidence` by 1.33–1.46×, while 0-bond 3-atom one-shot
    `multicenter` (+0.4–0.5 µs) and one-shot `branched/graph_overlays` (18.2→20.3 ms, +~230 ns/mol)
    reflect `S4c` per-call `PreparedSubstructure{Pattern,Host}::new` overhead when not reusing
    prepared carriers.

### Stage S5 — `Reaction` symmetry soundness & `PreparedSubstructureHost` integration (additive, green)
- **S5a (`umol-graph-ir`) — Done:** Wire `Reaction::apply*` / `React::react` so
  `reduce_pattern_symmetry = true` builds reaction-preserving `QuerySymmetryConstraints` via
  `Reaction::to_reaction_span()` + `reaction_span/symmetry.rs` (preserving deltas & stereo;
  conservative `{1}` fallback if span conversion fails), never bare `lhs` symmetry. **Additive
  (green).** `[dep: S4c]`
  - *Implementation notes & findings:* Added `ReactionSpan::proper_incidence_generators` in
    `umol-graph-ir/src/ir/reaction_span/symmetry.rs` over `IncidenceLevel::Full` using
    `initial_colors(..., DescriptionLevel::Structure, true)` and
    `ir::symmetry::extract_subgroup_generators` with a deterministic base covering all
    `StereoAtom`/`StereoBond` nodes, sites, and ligands plus any entities referenced by
    `ConstraintSpan`s (reusing `pub(crate) edit::collect_constraint_entities` and `pub(super)
    reaction_span::canonicalize::remapping_from_images`), falling back to `None` (`{1}`) on
    malformed coloring, non-normalizable stereo records, or non-reframable constraint spans. Added
    `Reaction::proper_lhs_generators` in `umol-graph-ir/src/ir/reaction.rs` projecting
    `ReactionSpan` generators onto `lhs.raw_graph()` (`GraphAndOverlays`) or
    `lhs.incidence_graph(Constitution).raw_graph()` (`Incidence`, mapping each LHS entity node
    through `span_incidence.entity -> lhs_incidence.node_of` so added atoms/overlays never misalign
    indices), short-circuiting empty deltas to `Molecule` proper symmetry and falling back
    conservatively to `{1}` when `to_reaction_span()` fails. Refactored
    `PreparedSubstructurePattern::new` onto a shared internal `build` helper with
    `PreparedSubstructurePattern::for_reaction` and wired `ReactionApplicationIter::new` (covering
    `apply`, `tracked_apply`, `apply_to_reaction`, `apply_to_reaction_span`, and `React::react`)
    through it. Verified green with `cargo +nightly fmt --all`, `git diff --check`, `cargo clippy -p
    umol-graph-ir --all-targets --features proptest -- -D warnings`, `umol-graph-ir --lib reaction`
    (1,437 passed, including 8 new `test_reaction_span_proper_incidence_generators` cases and 7 new
    `test_reaction_apply_symmetry_*` cases across both strategies, both automorphism algorithms, and
    all 4 subgraph-isomorphism algorithms), `--lib substructure` (111 passed), and `--test property
    --features proptest reaction` (105 passed).
- **S5b (`umol-graph-ir`) — Done:** Add `Reaction::{apply_to_host, tracked_apply_to_host,
  apply_to_host_to_reaction, apply_to_host_to_reaction_span}` and `impl React for
  PreparedSubstructureHost`; add tests for regioisomer preservation on symmetric LHS (e.g. benzene
  monosubstitution on toluene), chiral stereocenter preservation, symmetric bimolecular `A + A`
  deduplication, symmetric host deduplication (benzene, biphenyl, adipic acid), and parity with
  `Reaction::apply`. **Additive (green).** `[dep: S5a]`
  - *Implementation notes & findings:* Replaced `PreparedSubstructurePattern::for_reaction` with a
    decoupled `pub(crate) PreparedSubstructurePattern::with_pattern_symmetry(pattern, config, impl
    FnOnce(&Molecule, Option<&IncidenceGraph>) -> Vec<Vec<NodeId>>)` hook (removing
    `substructure/prepared.rs`'s dependency on `Reaction` with zero dispatch overhead) and private
    `Reaction::prepare_lhs` / `Reaction::proper_lhs_generators`. Added
    `ReactionApplicationIter::new_in_host`, `Reaction::{apply_to_host, tracked_apply_to_host,
    apply_to_host_to_reaction, apply_to_host_to_reaction_span}`, and `impl React for
    PreparedSubstructureHost` in `umol-graph-ir/src/ir/reaction.rs`, reusing
    `PreparedSubstructurePattern::matches_in_host` and `host.molecule().clone()`. Added
    `test_reaction_apply_to_host_symmetry_and_parity` (6 cases: benzene monosubstitution, biphenyl
    C–H monosubstitution, adipic acid monodeprotonation, symmetric bimolecular 3-copy `C–O + C–O`,
    `meso` host preserving both enantiomers, and chiral `C2` host reducing equivalent sites, swept
    across both match strategies, all 4 symmetry modes, both automorphism algorithms, and all 4
    subgraph-isomorphism algorithms with exact `Reaction::apply` parity),
    `test_reaction_apply_to_host_outputs_isolation_and_preconditions` (exact parity across all 5
    entry points, post-prepare host mutation isolation, and `ApplyPreconditionError` parity), and
    `test_reaction_apply_to_host_panics_on_lhs_molecule_constraints`. Verified green with `cargo
    +nightly fmt --all`, `git diff --check`, `cargo clippy -p umol-graph-ir --all-targets --features
    proptest -- -D warnings`, `umol-graph-ir --lib reaction` (1,445 passed, 3 ignored), `--lib
    substructure` (111 passed), `--test property --features proptest reaction` (105 passed, 1
    ignored), and `--doc reaction` (1 passed).

### Stage S6 — `umol-graph` `PatternFingerprinter`, defaults, and benchmarks (additive/default switch, green)
- **S6a (`umol-graph`) — Done:** Switch `PatternFingerprinter::new()` default to `BitsetRi`;
  precompile templates and use `PreparedSubstructureHost` + `visit_substructure_embeddings` across
  both `GraphAndOverlays` and `Incidence` and all 4 `SubgraphIsomorphismAlgorithm` variants; verify
  100% bit-exact RDKit 2026.03.3 test parity. **Additive/default switch (green).** `[dep: S4c]`
  - *Implementation notes & findings:* Switched `PatternFingerprinter::new()` default from
    `Vf2Rdkit` to `BitsetRi` in `umol-graph/src/fingerprint/pattern.rs`. Added a lock-free
    `[OnceLock<Vec<PreparedPatternTemplate>>; 8]` static cache indexed by
    `(SubstructureMatchAlgorithm, SubgraphIsomorphismAlgorithm, RelevantCycleEnumerationAlgorithm)`
    (keeping `reduce_pattern_symmetry = false, reduce_host_symmetry = false` for RDKit `uniquify =
    false` parity), prepared one `PreparedSubstructureHost` per molecule, pre-extracted ground
    `atom_numbers` and `bond_types` once per molecule after `is_concrete()`, visited matches via
    `PreparedSubstructurePattern::visit_embeddings_in_host`, and set bits directly in `BitFp`
    (narrowly widening `BitFp::{zeros,set}` to `pub(super)` in `bit_fp.rs`, removing per-call
    `Vec<u64>` sorting/dedup while preserving `NotConcrete`-before-`ZeroWidth` error precedence).
    Expanded `test_pattern_fingerprinter_fingerprint_algorithm` to the full `2 × 4 = 8`
    `(match_algorithm, subgraph_isomorphism_algorithm)` matrix across all 6 pinned RDKit 2026.03.3
    molecules plus empty-molecule and error-precedence cases, and aligned the two
    `GraphPatternFingerprinter::new()` roundtrip test rows in `umol-py/src/fingerprint/config.rs`
    (leaving Python `PatternFingerprintConfig()`'s default switch for `S7a`). Verified green with
    `cargo +nightly fmt --all`, `git diff --check`, `cargo clippy -p umol-graph -p umol-py
    --all-targets -- -D warnings` (under `umol-py/.venv` Python 3.13.18), `umol-graph --lib
    fingerprint` (108 passed), `umol-graph --test fingerprint` (26 passed), `umol-py --lib
    fingerprint` (156 passed), and `cargo bench -p umol-graph --bench fingerprint -- "pattern"`
    (`fingerprint_fixture/pattern` 978 ns; `fingerprint_structural/pattern/300` **4.42 ms**, −95.68%
    / **23.1×** faster vs prior ~102 ms baseline).
- **S6b (`umol-graph`, `umol-graph-ir`) — Done:** Update `umol-graph/benches/substructure.rs` (4
  algorithms + prepared pattern/host + screening + symmetry-reduced) and
  `umol-graph-ir/benches/reaction.rs` (`BitsetRi` + `apply_to_host` library benchmark); record
  before/after benchmark verification. **Additive (green).** `[dep: S5b, S6a]`
  - *Implementation notes & findings:* Expanded `umol-graph/benches/substructure.rs` with
    `substructure/prepared` (`build_hosts/corpus`, `pattern/*`, `both/*` across all 4 algorithms × 2
    strategies), `substructure/screening` (`direct_has`, `prepared_pattern_has`,
    `prepared_both_has`, `prepared_both_count` across `branched`, `phenol`, `bicyclic`, and 0-hit
    `steroid`), and `substructure/symmetry` (`ring6_corpus` and `ring6_coronene` across
    `none`/`pattern`/`host`/`both`). Expanded `umol-graph-ir/benches/reaction.rs` with
    `reaction/apply_to_host/all_matches` and a `reaction_library` group (16 realistic reactions
    covering acyclic, `RingScope::All`-constrained, symmetric-LHS, and rare-element pre-screened
    rules × 6 hosts, comparing `apply` vs `apply_to_host_per_run` vs `apply_to_host_prebuilt` under
    `Vf2Rdkit`/`BitsetRi` and unsym/sym modes). Verified `cargo +nightly fmt --all`, `git diff
    --check`, `cargo clippy -p umol-graph -p umol-graph-ir --benches -- -D warnings`, and `--
    --test` on both benches.
  - *Benchmark verification (Apple M1 Pro, `--warm-up-time 1 --measurement-time 2 --sample-size
    20`):*
    - **OpenSMILES ~9k all-matches (`GraphAndOverlays`):** `branched` direct `Vf2Rdkit` 33.8 ms →
      direct `BitsetRi` 19.7 ms (1.72×) → `prepared/pattern` 14.2 ms (2.39×) → `prepared/both`
      **5.07 ms** (**6.67×**; `Incidence` `prepared/both` 5.32 ms); `phenol` 52.1 → 20.6 → 13.6 →
      **4.97 ms** (**10.5×**; `Incidence` 5.31 ms); `bicyclic` 39.4 → 21.2 → 13.2 → **4.81 ms**
      (**8.20×**; `Incidence` 5.74 ms); `build_hosts/corpus` **11.9 ms** (~1.3 µs/mol).
    - **OpenSMILES ~9k screening (`GraphAndOverlays`, `BitsetRi`):** `prepared_both_has` **1.90 ms**
      (`branched`, vs `Vf2Rdkit` 7.84 ms, **4.12×**), **2.87 ms** (`phenol`, vs 23.3 ms, **8.13×**),
      **4.38 ms** (`bicyclic`, vs 21.7 ms, **4.95×**), **4.57 ms** (`steroid` 0-hit, vs 23.0 ms,
      **5.04×**); `prepared_both_count` 2.27 / 3.41 / 4.69 / 4.64 ms.
    - **Symmetry (`BitsetRi`, `GraphAndOverlays`, `PreparedSubstructureHost`):** `ring6_coronene`
      16.5 µs (`none`) → **3.16 µs** (`pattern`, **5.22×**), **6.15 µs** (`host`, **2.68×**), **5.88
      µs** (`both`, **2.81×**). On `ring6_corpus` (~9k organic mols, `count_matches_in_host`),
      `none` 10.3 ms → `pattern` **7.39 ms** (**1.39×**); `host` initially measured 138.5 ms (`both`
      81.5 ms) due to fixed-point `SchreierTree`/generator-cloning allocations in
      `TargetStabilizerChain::try_push_orbit_min`, unfiltered unrelated-substituent host generators,
      and per-call `SearchSymmetry` `Vec` allocations. Optimizing `TargetStabilizerChain`
      (zero-allocation fixed-point step sharing, single-involution fast path, lazy `SimsFilter`,
      `keep_tree` flag, and exact candidate orbit-closure generator filtering via precomputed
      `TargetSymmetry` masks) plus reusing `SearchSymmetry` buffers in `SubgraphIsomorphismScratch`
      cut `ring6_corpus` `host` from 138.5 ms → **21.4 ms** (**−84.6% / 6.48×**), `both` from 81.5
      ms → **18.2 ms** (**−77.7% / 4.47×**), and `ring6_coronene` `host`/`both` by 23–29%.
    - **Reaction library (16 rules × 6 hosts):** `apply/vf2rdkit_unsym` 695 µs →
      `apply/bitset_ri_unsym` 663 µs → `apply_to_host_per_run` 534 µs (1.30×) →
      `apply_to_host_prebuilt` **482 µs** (**1.44×**); with `(true, true)` symmetry: `apply` 1.26 ms
      → `apply_to_host_per_run` 973 µs (1.29×) → `apply_to_host_prebuilt` **837 µs** (**1.50×**). On
      the 1-shot 6-atom `reaction/apply/all_matches` fixture, `apply` is 4.56 µs vs `apply_to_host`
      4.24 µs (reflecting per-call reaction LHS preparation before `S6d` `PreparedReaction`).
- **S6c (`umol-graph-core`) — Done:** Eliminated remaining `TargetStabilizerChain` search-time
  allocations: zero-copy `Level0 { active_mask, dropped_supp }` borrowing from `&TargetSymmetry`
  for $n \le 128$ and $|S| \le 64$, exact disjoint-support isolated-involution fast path for
  multi-generator stabilizers ($\text{Stab}_{\langle g_m \rangle \times H}(v) = H$ via precomputed
  `gen_supports`, `gen_is_involution`, and `gen_disjoint_from` bitmasks, zero
  BFS/`SimsFilter`/`compute_orbits`/derived groups), compact inline
  `TargetSchreierTree::Involution` (zero heap allocation even with `keep_tree = true`), and
  reusable `SimsFilter` buffers across calls (`dirty_buckets`, `reset`). Threaded
  `target_sym: &TargetSymmetry` through `is_orbit_min_at_current` and
  `try_push_orbit_min_with_tree` across all four engines (`BitsetRi`, `Glasgow`, `RayKirsch`,
  `Vf2Rdkit`) and expanded unit coverage across `U64`/`U128`/`Wide` domains. **Internal
  optimization (green).** `[dep: S6b]`
  - *Verification:* `cargo +nightly fmt --all`, `git diff --check`, `cargo clippy -p
    umol-graph-core -p umol-graph-ir --all-targets --features proptest -- -D warnings` clean;
    `umol-graph-core` `--lib subgraph_isomorphism` (186 passed) and `--test property --features
    proptest subgraph_isomorphism` (27 passed); `umol-graph-ir` `--lib substructure` (111 passed),
    `--lib test_reaction_apply` (103 passed, 1 ignored), and `--test property --features proptest
    substructure` (42 passed).
  - *Benchmark verification (`BitsetRi`, `GraphAndOverlays`, Apple M1 Pro, `--warm-up-time 1
    --measurement-time 2 --sample-size 20`):*
    - `ring6_corpus`: `none` 10.5 ms, `pattern` **7.30 ms** (−0.9%), `host` 21.4 ms → **18.2 ms**
      (**−15.1%** vs post-S6b; **−86.9% / 7.63×** vs initial 138.5 ms), `both` 18.2 ms →
      **14.8 ms** (**−18.3%** vs post-S6b; **−81.9% / 5.51×** vs initial 81.5 ms).
    - `ring6_coronene`: `none` 16.6 µs, `pattern` **3.15 µs**, `host` 6.15 µs → **5.74 µs**
      (**−6.3%**, **2.89×** faster than `none`), `both` 5.88 µs → **5.43 µs** (**−8.8%**,
      **3.06×** faster than `none`).
- **S6d (`umol-graph-ir`, benches) — Done:** Added `PreparedReaction` (`Clone + Debug + Send +
  Sync`, `Arc`-backed $O(1)$ clone; `new(reaction, config) -> Result<Self,
  ApplyPreconditionError>`, `reaction`, `config`, and infallible `apply*` / `apply_to_host*`
  iterators), re-exported from `umol_graph_ir::ir`; added `React::react_prepared(&self, reaction:
  &PreparedReaction) -> ReactionProductsIter` for `Molecule`, `PreparedSubstructureHost`, and
  `[Molecule]`. Refactored `ReactionApplicationIter<T>` to store `Option<ReactionApplicationState>`
  with `Arc<Deltas>` (zero input clones on 0-match pre-screens; immediate state drop on error or
  exhaustion), moved `bond_endpoints_are_disjoint` to `tracked_apply_at` (already checked in
  `prepare_deltas`), and made `apply_deltas_at` borrow `&Deltas` via `Cow::Borrowed` when all
  deltas are `Atom`/`Bond`. Added `#[rstest]` unit, parity (across all 2 × 4 × 2 × 4 strategy /
  symmetry / automorphism / subiso combinations), isolation, error-precedence, span-fallback, and
  `React::react_prepared` tests plus `PreparedReaction` benchmarks in
  `umol-graph-ir/benches/reaction.rs`. **Additive (green).** `[dep: S5b, S6c]`
  - *Verification:* `cargo +nightly fmt --all`, `git diff --check`, `cargo clippy -p umol-graph-ir
    --all-targets --features proptest -- -D warnings` clean; `umol-graph-ir` `--lib
    ir::reaction::tests` (278 passed, 3 ignored) and `--doc` (2 passed, 3 ignored).
  - *Benchmark verification (`BitsetRi`, `GraphAndOverlays`, Apple M1 Pro, `--warm-up-time 1
    --measurement-time 2 --sample-size 20`):*
    - Single 6-atom aromatic rule (`reaction/*/all_matches`): `apply` 4.68 µs → `apply_to_host`
      4.20 µs → `prepared_apply` **2.73 µs** (1.71×) → `prepared_apply_to_host` **2.28 µs**
      (**2.05×**).
    - `reaction_library` (16 rules × 6 hosts), unsym: `apply/vf2rdkit_unsym` 662 µs →
      `apply/bitset_ri_unsym` 632 µs → `prepared_rxn` 516 µs → `apply_to_host_per_run` 492 µs →
      `apply_to_host_prebuilt` 467 µs → `prepared_both_per_run` **380 µs** (**1.74×** vs
      `Vf2Rdkit`) → `prepared_both_prebuilt` **351 µs** (**1.88×** vs `Vf2Rdkit`).
    - `reaction_library` (16 rules × 6 hosts), `(true, true)` sym: `apply/bitset_ri_sym` 1.20 ms →
      `apply_to_host_per_run` 897 µs → `apply_to_host_prebuilt` 809 µs → `prepared_rxn` 516 µs
      (2.33×) → `prepared_both_per_run` **250 µs** (**4.81×** vs `apply_sym`; **1.52×** faster than
      unsym `prepared_both_per_run`) → `prepared_both_prebuilt` **175 µs** (**6.86×** vs
      `apply_sym`; **2.00×** faster than unsym `prepared_both_prebuilt`; **3.78×** faster than
      `apply/vf2rdkit_unsym`).

### Stage S7 — `umol-py` bindings, defaults, documentation, and verification (additive/default switch, green)
- **S7a (`umol-py`) — Done:** Switched `SubstructureSearchConfig`, `ReactionApplicationConfig`, and
  `PatternFingerprintConfig` defaults to `SubgraphIsomorphismAlgorithm.BitsetRi()`, exposed
  `automorphism_algorithm` (default `AutomorphismAlgorithm.TrailIr()`), `reduce_pattern_symmetry`
  (default `False`), and `reduce_host_symmetry` (default `False`) on `SubstructureSearchConfig` and
  `ReactionApplicationConfig`, and exposed `Molecule::{has_substructure_match,
  first_substructure_match, count_substructure_matches}` (`umol-py/src/{substructure,reaction,
  molecule,fingerprint/config}.rs`). **Additive/default switch (green).** `[dep: S5b, S6a]`
- **S7b (`umol-py`, docs) — Done:** Bound frozen `PreparedSubstructureHost`,
  `PreparedSubstructurePattern`, and `PreparedReaction` (`umol-py/src/{substructure,reaction,
  molecule,lib}.rs`, `umol-py/python/umol/__init__.py`) with `SubstructureHostLike` (`Molecule |
  PreparedSubstructureHost`) on `PreparedSubstructurePattern` matching methods, `Reaction.apply*`,
  and `PreparedReaction.apply*`; exposed detached-copy `molecule()` / `pattern()` / `reaction()`
  methods (preserving `docs/development/python-api.md`'s property-mutation contract) and read-only
  `config` properties; exposed `Molecule::{react_prepared, react_all_prepared}` and
  `PreparedSubstructureHost::{react, react_prepared}`. Added Rust PyO3 tests and Python `pytest`
  coverage (`test_{substructure,reaction,fingerprint,import,class_constructor_signatures}.py`);
  updated `docs/development/{nomenclature,operations,python-api}.md` and `README.md`. **Additive
  (green).** `[dep: S6d, S7a]`
  - *Final verification gate (`umol-py/.venv`, Python 3.13.18):* `cargo +nightly fmt --all --
    --check`, `git diff --check`, `cargo clippy --workspace --all-targets -- -D warnings` (0
    warnings), `cargo test --workspace` (all workspace crates + doc tests green, including 1,792
    `umol-py` Rust unit tests), `umol-graph-core --test property --features proptest
    subgraph_isomorphism` (27 passed), `umol-graph-ir --test property --features proptest`
    (`substructure`: 42 passed; `reaction`: 105 passed, 1 ignored), `maturin develop` + `pytest -q
    umol-py/tests` (**3,726 passed, 2 skipped**), and README Python snippet verification.

### Critical path & verification gate
- **Critical path:** `S0a → S1a/S1b → S1c → S2a/S2b/S2c → S2d → S3a → S3b → S4a → S4b → S4c → S5a →
  S5b → S6c → S6d → S7b` (`S6a/S6b/S7a` run in parallel once their upstream deps land).
- **Final verification gate:** `cargo +nightly fmt --all -- --check`, `cargo clippy --workspace
  --all-targets -- -D warnings`, `cargo test --workspace`, feature-gated `proptest` suites for
  `umol-graph-core` and `umol-graph-ir`, `maturin develop` + `pytest -q umol-py/tests` under
  `umol-py/.venv` (Python 3.13), and targeted benchmark regression checks.

## Out-of-scope follow-up ideas

Orthogonal extensions and internal refinements that build cleanly on the settled `256` interfaces
without changing them (not part of `S0`–`S7`):

1. **Indexed / multiplexed batch search (`SubstructurePatternSet` and host screening indices):** An
   inverted feature/element index over large rule libraries ($R \gg 10^2$–$10^3$) to skip even the
   ~1–2 ns bitmask pre-screen per rule, and a fingerprint / structural-key pre-screen index over
   large molecule databases ($M \gg 10^4$).
2. **Component-wise multi-reactant filtering (`A + B` bimolecular enumeration):** Match disconnected
   reaction LHS components independently against per-molecule `PreparedSubstructureHost` instances
   before forming reactant pairs, avoiding $O(M^2)$ `Molecule::combine_all` constructions.
3. **Parallel batch wrappers and Python GIL release:** Leverage `PreparedSubstructurePattern` /
   `PreparedSubstructureHost` / `PreparedReaction`'s `Send + Sync` / frozen immutability for `rayon`
   batch methods in Rust and GIL-releasing single/batch screening calls in `umol-py`.
4. **Extended pre-screen bitmasks:** Add formal-charge, isotope, radical, stereo, or ring-presence
   summary masks inside `PreparedSubstructurePattern` / `PreparedSubstructureHost` if corpus
   profiling identifies queries that pass element/bond/overlay screening and reject late.
5. **Optional symmetry-breaking for `has` / `first` and revisiting Python defaults:** Re-evaluate
   conservative `(False, False)` Python defaults toward `(True, False)` for
   `SubstructureSearchConfig` and `(True, True)` for `ReactionApplicationConfig` once downstream
   usage patterns stabilize, and assess whether hard symmetric UNSAT `has`/`first` workloads warrant
   opt-in symmetry acceleration.
