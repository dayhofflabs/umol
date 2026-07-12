# 143 — Vendored nauty integration — 2026-07-11

## Decision

Replace the `nauty-Traces-sys` dependency of `umol-graph-core` with a workspace crate,
`umol-nauty-sys`, that vendors the required nauty C sources, compiles them with `cc`, and exposes a
small handwritten interface through an opaque C shim. This removes the build-time `bindgen` and
LLVM/libclang dependency that complicates Python source distributions and wheel builds.

A pure-Rust individualization/refinement implementation is not part of this work. The goal is to
simplify distribution without replacing a mature graph-canonization implementation. It can be
reconsidered only if the vendored nauty integration proves inadequate.

## Nauty, sparse nauty, and Traces

Nauty and Traces are alternative canonical-labeling/automorphism search engines shipped in the same
source package. They solve substantially the same problem but use different search strategies.
"Sparse nauty" is not Traces: it is nauty operating through its sparse-graph representation and
dispatch, invoked by `sparsenauty`. Traces has its own entry point and search implementation. Traces
accepts sparse undirected graphs, but unlike nauty it does not presently accept digraphs.

The current `umol-graph-core` backend calls `sparsenauty`, not `Traces`. Consequently the first
vendored crate should contain only the nauty sources needed by `sparsenauty` and its dependencies;
the Traces and gtools sources should not be compiled unless a later benchmark establishes a reason
to expose Traces as a second backend. This reduces native build time and the maintained C surface.

This distinction should also be reflected in naming: the crate is `umol-nauty-sys`, and the existing
algorithm variant remains `AutomorphismAlgorithm::Nauty`.

## Current contract

`Graph::automorphisms` currently supplies all of the following:

- a canonical labeling;
- the vertex orbit partition;
- automorphism-group order;
- a generating set of vertex permutations;
- caller-provided vertex colors;
- repeated stabilizer calculations by assigning a site a unique color.

All are required. In particular, generators are consumed by stereochemical orientation grading and
site stabilizers, while canonical labeling is consumed by exact graph keys. The replacement must not
quietly reduce the result to canonicalization or orbit discovery alone.

## Rust API cleanup

Introduce a backend result named `AutomorphismOutput` and write out `canonical_labels` in full:

```rust
pub struct AutomorphismOutput {
    orbits: Vec<NodeId>,
    canonical_labels: Vec<NodeId>,
    node_count: usize,
    orbit_count: usize,
    group_order: AutomorphismGroupOrder,
    generators: Vec<Vec<NodeId>>,
}
```

`Graph::automorphisms` should return `AutomorphismOutput`. If retaining `Automorphism` as a public
name is desirable for source compatibility, it can temporarily be a deprecated type alias; otherwise
the rename can be made directly while the API is still internal to the workspace. Rename
`canonical_lab` to `canonical_labels` internally and expose
`AutomorphismOutput::canonical_labels()`.

### Group order

The current `AutoGroupOrder::{Exact(u32), Approx(f64)}` has two problems: `u32` is unnecessarily
small, and collapsing nauty's `(grpsize1, grpsize2)` scientific representation into one `f64` loses
the representation supplied by the solver.

Retype it as a nauty-independent public value that preserves that representation, for example:

```rust
pub enum AutomorphismGroupOrder {
    Exact(u128),
    Scientific { mantissa: f64, exponent: i32 },
}
```

Use `Exact` only when the solver result can be converted without loss; otherwise retain its mantissa
and base-10 exponent. Provide conversion/display helpers rather than making callers reconstruct the
numeric convention. `u128` is ample for ordinary molecular use but is not presented as an unbounded
integer. If exact orders beyond `u128` become a requirement, add an optional big-integer conversion
or compute the order from the returned generators; neither is needed for the distribution fix.

## `umol-nauty-sys` design

Follow the workspace pattern established by `umol-msym-sys`, with one important constraint: do not
hand-copy nauty's internal structs into Rust. Keep `sparsegraph`, `optionblk`, `statsblk`, allocation,
macros, and version checks entirely on the C side.

The crate should contain:

```text
umol-nauty-sys/
├── Cargo.toml
├── build.rs
├── include/umol_nauty.h
├── src/lib.rs
└── nauty/                 # vendored upstream source and license
```

`build.rs` compiles the shim and the minimal nauty source closure using `cc::Build`. It defines
`USE_TLS` so separate Rust threads may execute nauty independently. Vendored version and license
must be recorded explicitly; nauty 2.6 and later is Apache-2.0.

The C API should be narrow and owned by umol rather than mirroring upstream ABI. Conceptually:

```c
typedef void (*umol_nauty_generator_fn)(
    void *context, const unsigned int *permutation, unsigned int n);

int umol_nauty_run(
    unsigned int n,
    const size_t *offsets,
    const unsigned int *neighbors,
    const unsigned int *colors,
    unsigned int *canonical_labels,
    unsigned int *orbits,
    umol_nauty_generator_fn report_generator,
    void *context,
    double *group_mantissa,
    int *group_exponent);
```

The exact integer widths should be finalized against nauty's supported vertex range and the Rust
graph limits. The wrapper must validate conversions before entering C.

The shim is responsible for:

1. constructing nauty's `lab`/`ptn` ordered partition from the supplied color ranks;
2. constructing or borrowing a nauty `sparsegraph` from CSR input;
3. selecting sparse undirected options, canonical-form generation, and the generator callback;
4. calling `nauty_check` and `sparsenauty`;
5. copying canonical labels, orbits, and group-size fields to caller-owned outputs;
6. freeing every nauty allocation before returning;
7. translating nauty failure states into a small stable umol error enum.

Passing an opaque callback context removes the current Rust `thread_local!` generator accumulator.
Rust can pass a pointer to a scoped `Vec<Vec<NodeId>>`; the callback copies each permutation while it
is valid. The safe wrapper owns all buffers and contains the raw FFI call and callback plumbing.

## Source selection

Start from the upstream source version already exercised by the workspace, but audit for the newest
Apache-licensed release before landing. Determine the minimal source closure empirically from
`sparsenauty` and link validation rather than copying the full source list from
`nauty-Traces-sys`. In particular, exclude `traces.c`, gtools programs, file-format utilities, and
unrelated command-line support unless required by the linker.

Keep the vendored upstream files unmodified where practical. Put portability fixes and the stable
API in `umol_nauty.c`/`umol_nauty.h`, so an upstream refresh is a replace-and-test operation.

## Integration steps

1. Establish the benchmark corpus and record a baseline through the existing `nauty-Traces-sys`
   integration. This lands first, so the wrapper replacement cannot erase its own comparison point.
2. Add `umol-nauty-sys` to the workspace and vendor the selected upstream release plus its license.
3. Implement and unit-test the C shim independently of `umol-graph-core`.
4. Add the handwritten Rust declarations and a safe crate-local wrapper around `umol_nauty_run`.
5. Replace `nauty-Traces-sys` in `umol-graph-core/Cargo.toml`.
6. Refactor `auto.rs` to produce `AutomorphismOutput`, `canonical_labels`, and the new group-order
   representation while preserving `AutomorphismAlgorithm::Nauty`.
7. Run the same benchmark corpus through the vendored integration and investigate material
   regressions in graph conversion, callback collection, or nauty configuration before removing the
   old dependency.
8. Remove bindgen from this dependency path and confirm with `cargo tree -i bindgen`.
9. Validate Rust tests and Python builds on Linux, macOS, and Windows without LLVM/libclang.

## Performance workstream

Performance evaluation is a prerequisite and migration gate, not post-hoc polish. Sparse nauty is
the only production algorithm initially, so evaluate it against umol's workload and compare the new
shim with the existing wrapper around the same solver. This isolates integration overhead and
configuration mistakes even though it is not an algorithm-to-algorithm comparison.

Extend the existing `umol-graph-core` Criterion group before implementing the replacement. Cover:

- representative molecules and full incidence graphs across realistic sizes;
- rigid graphs and highly symmetric graphs;
- disconnected repeated components;
- uniform, chemically representative, and unique vertex colorings;
- ordinary automorphism calculation;
- uniquely colored site-stabilizer runs, including a sequence of stabilizers for one graph;
- `canonical_key`, whose edge subdivision increases the nauty graph size;
- scaling/stress families: paths, cycles, complete graphs, grids, hypercubes, and selected strongly
  regular or otherwise refinement-resistant graphs.

Keep construction outside the timed region unless graph construction itself is intentionally under
test. Report latency/throughput by graph size and workload class; retain Criterion distributions so
outliers on difficult graphs are visible. Record the compiler profile, target, CPU, nauty version,
and relevant compile definitions with each durable baseline.

Define acceptance in two layers:

1. **Migration parity:** the shim should introduce no material regression relative to
   `nauty-Traces-sys` on the same `sparsenauty` workload. Any threshold should allow normal benchmark
   noise but be fixed before measuring the new implementation.
2. **Workload fitness:** representative molecule, incidence, stabilizer, and canonical-key calls
   should meet explicit latency targets chosen for the Python-facing workflows. This establishes
   whether sparse nauty is adequate even without a competing backend.

Keep the benchmark input representation backend-neutral. If difficult families or real workloads
later fail the fitness target, add an opaque `umol_traces_run` beside `umol_nauty_run`, feed both the
same CSR/color inputs, and add `AutomorphismAlgorithm::Traces`. The common
`AutomorphismOutput` remains unchanged. Traces-specific options, statistics, callback plumbing, and
source files stay in the C shim. Canonical permutations and generator sets may differ between the
engines, so compare canonical graph/key, orbit and group semantics, correctness, and runtime rather
than raw output arrays.

### S5 benchmark result (2026-07-11)

The post-migration corpus was run on the same baseline host and compiler recorded in
`umol-graph-core/benches/README.md`, with Criterion's default 3-second warm-up, 100 samples, and
5-second measurement period. The pre-migration measurements were retained in Criterion's `base`
dataset and the new measurements were saved as `vendored-nauty`. Raw Criterion data remains local.

Across 69 comparable cases, the median change was a 1.2% slowdown. The distribution was mixed:
14 cases were more than 5% slower, four were more than 10% slower, and 15 were more than 5% faster.

| Representative case | Old binding | Vendored | Change |
|---|---:|---:|---:|
| ordinary path 64, degree colors | 7.07 µs | 8.67 µs | +22.7% |
| ordinary grid 8×8, degree colors | 7.22 µs | 8.44 µs | +16.9% |
| C60, three site stabilizers | 19.51 µs | 21.49 µs | +10.1% |
| dodecahedron, three site stabilizers | 8.33 µs | 9.26 µs | +11.2% |
| C60 incidence graph | 27.47 µs | 28.09 µs | +2.3% |
| C60 canonical key | 31.84 µs | 30.56 µs | −4.0% |
| C70 canonical key | 36.61 µs | 36.04 µs | −1.6% |
| C60, unique colors | 3.70 µs | 2.78 µs | −25.0% |

The nauty version, `WORDSIZE=64`, `USE_TLS`, and sparse entry point are unchanged. The principal
adapter-level difference is extra boundary work: graph-core ranks colors and builds CSR, while the C
shim validates and copies CSR, allocates nauty arrays, and sorts the already-ranked colors again.
This is consistent with fixed overhead being visible in some small and medium cases, but the mixed
results do not isolate one cause.

This run does not formally close the performance gate: the design required a regression threshold
fixed before measurement, but no numeric threshold was recorded before S5. The four regressions
above 10% should therefore be reviewed or rerun after eliminating redundant color sorting/CSR work
before sparse nauty is declared accepted for the Python-facing workload.

### S5c no-sort result (2026-07-11)

S5c passes graph-core's Rust-computed vertex partition order through `NautyInput` and the private C
ABI. The safe Rust boundary validates its length, vertex bounds, uniqueness, and nondecreasing color
order; the C boundary repeats the memory-safety-critical checks. The shim now initializes `lab` and
`ptn` directly and contains no `qsort` or colored-vertex comparison machinery.

The identical 69 cases were saved as `vendored-nauty-no-sort`. Relative to the old-binding `base`
dataset, the median change improved from +1.2% to +0.8%. Cases slower by more than 5% fell from 14
to five, no case remained slower by more than 10%, and 19 cases were faster by more than 5%.
Relative to the first vendored run, the median case improved by 1.0%.

| Representative case | Old binding | First vendored | No-sort | No-sort vs old | No-sort vs first |
|---|---:|---:|---:|---:|---:|
| ordinary path 64, degree colors | 7.07 µs | 8.67 µs | 6.72 µs | −4.9% | −22.5% |
| ordinary grid 8×8, degree colors | 7.22 µs | 8.44 µs | 7.42 µs | +2.7% | −12.1% |
| C60, three site stabilizers | 19.51 µs | 21.49 µs | 18.41 µs | −5.6% | −14.3% |
| dodecahedron, three site stabilizers | 8.33 µs | 9.26 µs | 8.55 µs | +2.7% | −7.7% |
| C60 incidence graph | 27.47 µs | 28.09 µs | 27.77 µs | +1.1% | −1.1% |
| C60 canonical key | 31.84 µs | 30.56 µs | 28.88 µs | −9.3% | −5.5% |
| C70 canonical key | 36.61 µs | 36.04 µs | 33.55 µs | −8.4% | −6.9% |
| grid 8×8 canonical key | 38.20 µs | 37.63 µs | 34.99 µs | −8.4% | −7.0% |

All four previous regressions above 10% disappeared. The largest remaining old-binding regression
was C60 with degree colors at +6.6%; the other four cases above 5% ranged from +5.1% to +6.1%.

**Acceptance decision (2026-07-12):** the S5 performance gate is accepted. Although no numeric
threshold was fixed before measurement, the accepted basis is the no-sort result: +0.8% median
change across the full corpus, no regression above 10%, only five regressions above 5%, elimination
of all four initially material regressions, and representative incidence, stabilizer, and
canonical-key latencies at or near old-binding parity. Sparse nauty is therefore considered adequate
for the Python-facing workloads in scope. S5b's LLVM-free platform/package matrix remains a separate
distribution acceptance gate.

## Verification

Retain the existing graph cases and add boundary tests at the shim seam:

- empty and single-vertex graphs (handled in Rust if nauty requires nonempty input);
- uniform and fully distinct color partitions;
- path, cycle, complete graph, cubane, Petersen, and fullerene cases;
- every generator is a bijection preserving colors and adjacency;
- returned generators induce the returned orbit partition;
- canonical keys are invariant under many input relabelings;
- a uniquely colored site is fixed by every returned stabilizer generator;
- parallel calls do not share generator state;
- group orders exercise exact and scientific representations;
- invalid sizes and failed integer conversions return errors rather than truncating.

For migration confidence, freeze representative outputs from the current binding before replacing
it. Generator sets themselves are not unique, so compare the generated group semantics, group order,
orbits, and canonical graph/key rather than requiring byte-identical generator lists.

## Non-goals

- implementing Bliss or a new individualization/refinement engine in Rust;
- exposing the complete nauty/Traces API;
- using bindgen during normal builds;
- adding Traces merely because its sources ship beside nauty;
- guaranteeing that canonical label permutations remain numerically identical across future nauty
  versions. Stored canonical keys that require long-term stability must pin the vendored version.
