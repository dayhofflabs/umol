# Experiment A: matching-enumeration correctness implementation plan

## Objective and boundary

Implement the correctness baseline from discussion 145 before changing kekulization chemistry or
attempting Uno:

- replace the unsound greedy branch-pruning test with an exact Edmonds extension oracle;
- validate perfect/maximum enumeration against an independent exhaustive oracle on all small graphs;
- expose bounded-memory, early-stoppable enumeration while preserving the current eager APIs;
- establish relabeling invariance and prescribed-hole correctness at graph-core level;
- add planar FKT/Pfaffian counting as a genuinely independent count oracle;
- add matching transport through `GraphCorrespondence`, the localized correspondence fix required
  when matching work proceeds;
- benchmark and record the baseline that later Uno work must beat.

This plan does not implement chemical electron-demand derivation, charge localization, Uno, general
maximal-matching enumeration, symmetry quotienting, or an automatic planarity algorithm. FKT accepts
an explicit validated planar embedding; supplying/deriving chemical embeddings beyond the shared
fixtures is separate work.

The workspace is green after every subitem. All test edits follow the test-writing conventions:
`rstest`, inline case construction, definition-parallel order, specific structural assertions, and no
helper constructors. Existing nonconforming tests in a module touched by a subitem are normalized in
that subitem rather than copied.

## S0 — Trusted inputs and reference vocabulary

### S0a — Normalize matching tests and add exhaustive reference enumeration

**Module:** `umol-graph-core/src/algorithms/matching.rs` test module

**Kind:** additive test foundation (green)

**Dependencies:** [dep: none]

Normalize the existing matching tests in definition order: convert remaining bare `#[test]` cases to
`#[rstest]`, replace the `ids` constructor helper with inline `NodeId` vectors, keep assertion-only
helpers only where they express a reusable invariant, and consolidate cases by the public method
under test.

Add a test-only exhaustive oracle that scans edge subsets and returns canonical sorted edge-id sets
for:

- every matching;
- all matchings of an exact cardinality;
- perfect matchings;
- maximum matchings.

The oracle must not call `maximum_matching`, `enumerate_*`, `can_reach`, or any production matching
helper. Validate the oracle itself against closed-form/hand-listed cases: empty graph, isolated
vertices, paths/cycles, triangle, `K4`, disconnected edges, and a graph with maximal matchings of
different cardinalities. Add explicit tests demonstrating that greedy matching cardinality is not an
upper bound on maximum cardinality; this documents the production bug without asserting buggy
production output.

**Verification:** `cargo test -p umol-graph-core algorithms::matching::tests`; Clippy on the touched
test target with warnings denied.

### S0b — Consolidate reusable matching graph fixtures

**Module:** `umol-graph-core/test-data/matching/`, graph-core test/benchmark fixture loaders,
`umol-graph-core/benches/algorithms.rs`, and the existing coronene test consumer

**Kind:** additive fixture consolidation (green)

**Dependencies:** [dep: S0a]

Create shared declarative graph-only fixture assets needed by Experiment A: benzene, naphthalene,
coronene with its face embedding, azulene, C60 with its face embedding, disconnected cycles,
ladders, and small grids. Use a minimal stable edge/face-list format that both graph-core and the
higher-crate tests can consume with local test-only loaders; do not add a production dependency or a
test-support feature to graph-core. Move the existing C60 benchmark topology into the shared asset
rather than copying it and adapt the automorphism benchmark loader without changing its measured
graph. Make the existing coronene aromaticity test consume the same topology asset while retaining
its molecule-specific atom construction locally.

Fixtures return topology/embedding data only and contain no expected algorithm outputs. Validate
node/edge counts, face-edge incidence, Euler characteristic, and selected degrees so a malformed
fixture cannot become a shared false oracle.

**Verification:** existing graph-core benchmarks compile with `cargo bench -p umol-graph-core
--bench algorithms --no-run`; focused fixture tests and existing aromaticity fixture consumers remain
green.

**Stage exit:** the exhaustive oracle and shared fixtures are trusted independently of production
enumeration; the workspace remains green.

## S1 — Correspondence-aware matching transport

### S1a — Add fallible matching transport through `GraphCorrespondence`

**Module:** `umol-graph-core/src/algorithms/matching.rs`,
`umol-graph-core/src/correspondence.rs`, and `umol-graph-core/src/lib.rs`

**Kind:** additive API (green)

**Dependencies:** [dep: S0a]

Add a small graph-core transport API with an explicit error type. A `Matching` on the
correspondence's left graph can be transported to the right graph only when:

- every matched left edge has a right edge mate;
- both endpoints have node mates consistent with that edge mate;
- transported edges remain pairwise vertex-disjoint;
- graph/count/id bounds agree.

The result is a right-side `Matching` whose mate vector spans the right graph; right-only vertices
remain exposed. Return specific errors for exposed matched edges/nodes, inconsistent endpoint/edge
images, and a non-matching image. Do not silently drop entities and do not reinterpret `Matching` as
a `Correspondence<NodeId>`.

Prefer an inherent operation such as
`Matching::transport(&GraphCorrespondence, left: &Graph, right: &Graph)` so matching semantics and
construction invariants remain owned by `Matching`. Re-export only the error type needed by callers.

Tests cover total isomorphism, induced-subgraph→host transport, right-only exposed vertices,
relabeling, each error variant, and round-trip transport through a total correspondence and its
reverse. Normalize every edited correspondence/matching test to the test-writing conventions and
keep tests parallel to definition order.

**Verification:** focused matching/correspondence tests, `cargo test -p umol-graph-core`, and Clippy
with warnings denied.

**Stage exit:** matching results can cross extraction correspondences without manual id vectors;
there is no behavior change to enumeration.

## S2 — Exact extension oracle and correct eager enumeration

### S2a — Add residual constrained-graph construction

**Module:** `umol-graph-core/src/algorithms/matching.rs`

**Kind:** additive internal primitive (green)

**Dependencies:** [dep: S0a, S1a]

Introduce one internal state representation for binary partition search: included edges, excluded
edges, covered vertices, and included cardinality. Centralize include/undo and exclude/undo operations
so conflicting included edges and all edges incident to newly covered vertices are handled once.

Add a residual-graph builder that compacts uncovered vertices, retains exactly the nonexcluded edges
between them, and returns its residual→original `GraphCorrespondence`. Keep compaction deterministic
in original node/edge order. Validate the builder against hand-written states and exhaustive small
states: residual topology, correspondence images, covered/exposed sets, and restoration after undo.

Do not wire this primitive into public enumeration yet.

**Verification:** focused internal tests and full graph-core tests.

### S2b — Add an exact Edmonds extension decision

**Module:** `umol-graph-core/src/algorithms/matching.rs`

**Kind:** additive internal primitive (green)

**Dependencies:** [dep: S2a]

Implement `can_extend_to(target_size)` by running Edmonds maximum matching on the compact residual
graph and testing
`included_size + residual_maximum.size() >= target_size`. Preserve cheap sound prechecks (included
size, uncovered vertex capacity, and impossible conflicts) ahead of Edmonds. Remove the greedy count
from the decision path; a greedy result may be retained only as a lower bound or initial matching,
never as a pruning upper bound.

Cross-check the decision against the exhaustive oracle for every valid include/exclude state of all
simple graphs through the practical exhaustive bound (recommended `V <= 5`), including disconnected
and non-bipartite graphs. Add a fixed counterexample where the old greedy decision rejects a state
that the exact oracle accepts.

**Verification:** exhaustive state test completes within an explicit practical test-time budget;
focused tests, full graph-core tests, and Clippy are green.

### S2c — Rewire eager perfect/maximum enumeration to the exact oracle

**Module:** `umol-graph-core/src/algorithms/matching.rs`

**Kind:** correctness fix (green)

**Dependencies:** [dep: S2b]

Replace `can_reach` and the duplicated branch bookkeeping in `enumerate_rec` with the search state and
exact extension decision. Preserve the public `Vec<Matching>` methods and their deterministic edge
branch order. Emit only at exact target cardinality and construct each result from the included-edge
state.

Cross-check complete canonical edge-set equality—not just counts—between production and exhaustive
enumeration for every simple graph through `V <= 5`. Test both perfect and maximum outputs,
including empty/isolated graphs, and assert validity, cardinality, uniqueness, and completeness.
Run property tests on sampled larger graphs and edge insertion orders to extend coverage without
making the exhaustive suite unbounded.

**Verification:** focused exhaustive/property tests, `cargo test -p umol-graph-core --features
proptest`, and Clippy with warnings denied.

**Stage exit:** current eager enumeration is complete and correct against an independent oracle for
the exhaustive bound; no greedy pruning remains.

## S3 — Streaming, early stop, and relabeling invariance

### S3a — Add a bounded-memory visitor enumeration API

**Module:** `umol-graph-core/src/algorithms/matching.rs` and `umol-graph-core/src/lib.rs`

**Kind:** additive API (green)

**Dependencies:** [dep: S2c]

Generalize the recursive exact search to accept a visitor returning `ControlFlow<B>`. Add public
`visit_perfect_matchings` and `visit_maximum_matchings` methods (names may follow the crate's existing
visitor convention if one exists) that pass each owned `Matching` to the visitor and immediately
propagate `Break`. Return the visitor's `ControlFlow` so callers can attach an early-stop value.

The traversal keeps only search state/undo logs plus the current output; it must not allocate a
hidden output `Vec`. Document deterministic traversal as an implementation property, not a canonical
ordering contract.

Tests cover zero outputs, one output, full traversal, stop before/after the first output, stop after
`k`, and equivalence with the eager set. A visitor that retains nothing demonstrates that result
storage does not grow with output count; use structural/state assertions rather than a fragile
allocator test.

**Verification:** focused visitor tests, full graph-core tests, and Clippy.

### S3b — Implement eager APIs as visitor collectors

**Module:** `umol-graph-core/src/algorithms/matching.rs`

**Kind:** internal rewire (green)

**Dependencies:** [dep: S3a]

Make `enumerate_perfect_matchings` and `enumerate_maximum_matchings` thin compatibility collectors
over the visitor APIs. Delete the separate eager recursion so there is one enumeration engine and
one extension decision. Assert byte-for-byte equality with the pre-rewire deterministic result order
on stable fixtures, while completeness tests continue to compare sets.

**Verification:** focused eager/visitor equivalence tests and full graph-core tests.

### S3c — Add relabeling and hole-deletion conformance

**Module:** `umol-graph-core/tests/matching_enumeration.rs` and shared fixtures from S0b

**Kind:** additive integration conformance (green)

**Dependencies:** [dep: S0b, S3b]

For deterministic and property-generated relabelings, enumerate the original and relabeled graphs,
map relabeled matching edges back through `GraphCorrespondence`, and compare canonical solution sets.
Cover perfect and maximum enumeration on bipartite/non-bipartite and disconnected graphs.

Model prescribed holes at graph-core level by extracting/deleting chosen vertices, enumerating
perfect matchings of the remainder, transporting results back, and asserting exactly those vertices
are exposed while every retained vertex is covered. Include one mobile-hole reference check by
enumerating each single-vertex deletion and comparing its union with the graph's near-perfect
maximum matchings. This validates the formulation for Experiment B without adding chemical demand
types.

**Verification:** integration tests with `proptest`, full graph-core tests, and workspace check.

**Stage exit:** enumeration is complete, bounded-memory when visited, early-stoppable, invariant
under relabeling, and validated for fixed/mobile one-hole graph formulations.

## S4 — Independent planar Pfaffian count oracle

### S4a — Add and validate an explicit planar embedding carrier

**Module:** new `umol-graph-core/src/algorithms/planar_matching_count.rs`, algorithms module exports,
and shared planar fixtures

**Kind:** additive API foundation (green)

**Dependencies:** [dep: S0b]

Add a `PlanarEmbedding` carrier whose face boundary walks reference graph node/edge IDs and identify
the outer face. Construction validates bounds, consecutive adjacency, edge incidence (twice for a
connected sphere embedding, with explicit bridge handling if supported), consistent orientation,
connectedness assumptions, and Euler characteristic. Reject malformed/noncellular embeddings with a
specific error enum.

The API deliberately accepts an embedding and does not claim to test planarity or discover one.
Document this boundary. Tests cover cycle, `K4`, cube, coronene, and C60 embeddings plus each invalid
shape. Keep embedding tests independent of matching enumeration.

**Verification:** focused embedding tests and full graph-core tests.

### S4b — Implement exact Kasteleyn signing and Pfaffian arithmetic

**Module:** `umol-graph-core/src/algorithms/planar_matching_count.rs` and
`umol-graph-core/Cargo.toml`

**Kind:** additive algorithm (green)

**Dependencies:** [dep: S4a]

Add the minimal pure-Rust exact-integer dependency required for counts beyond machine width
(`num-bigint`/workspace-consistent equivalent). Construct Kasteleyn face-parity equations over
`GF(2)`, solve for edge signs deterministically while omitting the one redundant outer-face equation,
and build the signed skew-symmetric matrix. Compute its Pfaffian with fraction-free/exact arithmetic;
return a nonnegative arbitrary-precision count.

Keep the signing solver and Pfaffian kernel separately testable. Kernel tests use explicit
skew-symmetric matrices with hand-computed Pfaffians, row/column permutation sign changes, zero
matrices, and odd dimensions. Signing tests assert every bounded face's required parity and remain
independent of the final count.

Expose a narrow method such as `Graph::count_perfect_matchings_planar(&PlanarEmbedding)` returning
`Result<BigUint, PlanarMatchingCountError>`. Do not add it to `MatchingEnumerationAlgorithm`: counting
is not enumeration.

**Verification:** focused arithmetic/signing tests, full graph-core tests, Clippy, and dependency-tree
review confirming only pure-Rust additions.

### S4c — Cross-validate FKT counts and hole-count identities

**Module:** `umol-graph-core/tests/matching_count.rs`, shared fixtures, and doc 145

**Kind:** additive independent validation (green)

**Dependencies:** [dep: S3c, S4b]

Compare FKT counts with:

- closed-form/hand counts for even cycles, disconnected cycle products, `K4`, and cube;
- exhaustive subset counts on all planar graphs in the practical small bound for which an embedding
  fixture/generator is available;
- distinct exact-oracle enumeration counts for benzene, naphthalene, coronene, and C60;
- prescribed-hole counts `PM(G-H)` on heterocycle-shaped fixtures;
- the one-mobile-hole identity `sum_v PM(G-v)` versus the number of near-perfect matchings, ensuring
  each output is counted once at its exposed vertex.

Known coronene/C60 counts are cross-checks, not the sole oracle. Test disconnected multiplication and
arbitrary-precision overflow beyond `u128`. State explicitly that counts are labeled; symmetry
quotienting is out of scope.

**Verification:** focused count integration tests, full graph-core tests with `proptest`, workspace
check, and Clippy.

**Stage exit:** planar counts are independent of enumeration and agree across exhaustive, known, and
chemical stress fixtures.

## S5 — Baseline performance and acceptance record

### S5a — Add output-sensitive enumeration benchmarks

**Module:** `umol-graph-core/benches/algorithms.rs`, optional focused matching benchmark module, and
`umol-graph-core/benches/README.md`

**Kind:** additive benchmark (green)

**Dependencies:** [dep: S3c, S4c]

Benchmark the exact-oracle visitor on the Experiment A corpus. Separate:

- preprocessing/first output (visitor immediately breaks);
- bounded prefixes (`k = 1, 10, 100` where available);
- full streaming fold with no retained outputs;
- eager collection compatibility cost;
- FKT count time;
- total outputs and per-output throughput.

Criterion cannot directly summarize inter-output delay inside one iteration; add a diagnostic runner
or measurement mode that records median/p95/max delay outside Criterion's timed loop without
committing raw machine data. Include coronene and C60, nonalternant graphs, prescribed-hole deletions,
disconnected products, grids/ladders, and dense bipartite graphs. Construction and embedding
validation stay outside timed closures.

**Verification:** benchmark targets compile and a reduced sample run completes.

### S5b — Run and record the Experiment A baseline

**Module:** doc 145 or an adjacent benchmark note

**Kind:** additive verification (green)

**Dependencies:** [dep: S5a]

Run the corpus on a recorded host/compiler/profile. Record preprocessing, first-output latency,
prefix/full throughput, delay distribution, peak search-state/output-retention behavior, FKT time,
and independently verified output counts. Treat C60's 12,500 labeled matchings and coronene's known
count as acceptance checks only after FKT/exhaustive cross-validation succeeds.

Acceptance requires:

- exact equality with exhaustive solution sets through the small-graph bound;
- no relabeling or hole-conformance failures;
- visitor early stop returns without exploring the remaining enumeration tree;
- FKT and enumeration counts agree on every planar fixture;
- no benchmark case exhibits an unexplained super-polynomial gap between outputs at the tested
  sizes (total exponential output volume is expected).

Run `cargo fmt --all`, graph-core tests with `proptest`, workspace tests (subject to the documented
local Python-link limitation), Clippy with warnings denied, and `git diff --check`.

**Stage exit:** Experiment A is the trusted general-graph reference and benchmark baseline for Uno;
the workspace is green.

## Critical path

`S0a → S1a → S2a → S2b → S2c → S3a → S3b → S3c → S4c → S5a → S5b`, with
`S0a → S0b → S4a → S4b → S4c` as the independent counting branch joining at S4c.

S0b may proceed in parallel with S1/S2 after S0a. S4a/S4b may proceed in parallel with S2/S3 once
the shared fixtures exist. S4c is the join point because it compares the independent counter with
the corrected enumerator.

## Deferrable work

- Automatic planar embedding/planarity testing is deferrable; explicit embeddings suffice for the
  Experiment A oracle and chemical fixtures.
- FKT acceleration through modular arithmetic/CRT, rank-one monomer updates, or reused factorizations
  is deferrable until counting benchmarks justify it.
- Full layered Hopcroft–Karp is deferrable; Edmonds supplies the correctness oracle.
- Symmetry-inequivalent counting via Burnside's lemma is deferrable and must not alter labeled counts.
- Uno 1997/2001, general maximal matching, chemical `MatchingDemand`, charge localization, and the
  kekulizer migration belong to later experiments/stages.
