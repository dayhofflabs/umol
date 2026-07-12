# 144 — Vendored nauty implementation plan — 2026-07-11

Implements the settled design in doc 143: replace `nauty-Traces-sys` with a vendored
`umol-nauty-sys` crate and an opaque C shim around `sparsenauty`; rename the result to
`AutomorphismOutput`, write out `canonical_labels`, and preserve nauty's group-order representation.
A pure-Rust canonization engine and Traces are not required for the core deliverable.

The stages are cut so the workspace is green after every stage. New benchmarks, the native crate,
the safe wrapper, and the new Rust vocabulary land additively before the one breaking backend/API
migration.

## S0 — Baseline and migration gates

The existing `nauty-Traces-sys` backend remains untouched. This stage makes its behavior and cost
observable before there is a replacement to bias the comparison.

### S0a — Backend-neutral benchmark corpus

**Module:** `umol-graph-core/benches/algorithms.rs`  
**Kind:** additive (green)  
**Dependencies:** none

Extend the Criterion automorphism group with reusable, backend-neutral benchmark cases:

- ordinary automorphism calls on representative rigid and symmetric molecular graphs;
- uniform, chemically representative, and unique color partitions;
- disconnected repeated components;
- sequences of uniquely colored stabilizer calls on one graph;
- `canonical_key` calls, including the cost of the edge-subdivision representation;
- scaling families: paths, cycles, complete graphs, grids, hypercubes, and a selected
  refinement-resistant/strongly-regular family;
- representative full incidence-graph sizes, constructed directly in graph-core benchmark terms so
  the foundation crate does not gain a dependency on `umol-ast`.

Construct graphs and color vectors outside the timed closure. Group benchmark IDs by workload,
graph, size, and coloring so the same corpus can later compare Nauty and Traces. Document the
baseline command, compiler profile, target, CPU, nauty version, and compile definitions in the
benchmark module or an adjacent checked-in README; do not commit Criterion's machine-specific raw
measurements.

**Verification:** run `cargo bench -p umol-graph-core --bench algorithms automorphism` and the new
canonical-key/stabilizer filters; retain one named Criterion baseline from the existing binding for
the post-migration comparison.

### S0b — Backend-semantic conformance cases

**Module:** `umol-graph-core/src/algorithms/auto.rs` tests and, where appropriate,
`umol-graph-core/tests/property.rs`  
**Kind:** additive (green)  
**Dependencies:** none

Strengthen the current backend tests before changing it:

- assert every returned generator is a bijection preserving adjacency and input colors;
- assert generator-induced components reproduce the returned orbit relation;
- assert canonical keys/forms are invariant under input relabeling;
- assert every generator from a uniquely colored stabilizer fixes that site;
- cover disconnected components and a graph whose group order does not fit `u32`;
- exercise concurrent calls so generator capture cannot leak between threads.

Do not freeze a particular generator set: generating sets are non-unique. Likewise, compare the
canonical graph/key rather than requiring a backend-independent numeric canonical permutation. New
or edited tests use `rstest`, table cases, specific assertions, and definition-parallel ordering per
the test-writing conventions.

**Verification:** `cargo test -p umol-graph-core` and
`cargo test -p umol-graph-core --features proptest --test property`.

**Stage exit:** graph-core tests are green; the old backend has a durable semantic oracle and a
recorded performance baseline.

## S1 — Vendored native foundation

This stage adds a workspace crate that builds but has no consumer. It cannot affect graph-core's
existing backend.

### S1a — Crate, upstream sources, and licensing

**Module:** workspace `Cargo.toml`; new `umol-nauty-sys/Cargo.toml`, `build.rs`, and
`umol-nauty-sys/nauty/`  
**Kind:** additive (green)  
**Dependencies:** none

Add `umol-nauty-sys` as a workspace member. Vendor nauty 2.9.3, its Apache-2.0 license, version
record, and only the source/header closure needed by `sparsenauty`. Keep upstream files unchanged.
Compile with `cc::Build`, `USE_TLS`, warnings suppressed for upstream code only, and the same
`WORDSIZE` behavior currently selected by `nauty-Traces-sys`. Emit precise `rerun-if-changed`
directives.

Establish the minimal closure by link validation. Exclude `traces.c`, Traces-only support, gtools,
command-line programs, and unrelated utilities unless a symbol dependency proves necessary.

**Verification:** `cargo check -p umol-nauty-sys` on the local target; inspect the build log/source
list to confirm Traces and gtools are absent.

### S1b — Opaque C shim

**Module:** new `umol-nauty-sys/include/umol_nauty.h` and
`umol-nauty-sys/src/umol_nauty.c`  
**Kind:** additive (green)  
**Dependencies:** [dep: S1a]

Define the stable umol-owned C ABI and its error enum. `umol_nauty_run` accepts validated CSR
offsets/neighbors, ranked vertex colors, caller-owned canonical-label/orbit buffers, a generator
callback plus opaque context, and group-order mantissa/exponent outputs.

The implementation owns all nauty details: construct `sparsegraph`, construct the ordered
`lab`/`ptn` partition, configure sparse undirected canonicalization and the automorphism callback,
call `nauty_check` and `sparsenauty`, copy outputs, and release allocations on every path. Convert
nauty status into stable shim errors. Add compile-time assertions for the integer-width assumptions
shared with the header.

Keep the function name backend-specific. A future Traces backend gets `umol_traces_run`, not a
mode flag inside `umol_nauty_run`.

**Verification:** C compilation with warnings enabled for the shim; a small Rust-side invocation is
deferred to S2, where outputs can be asserted safely.

**Stage exit:** the new native archive and shim build on their own; the workspace remains green and
graph-core still uses `nauty-Traces-sys`.

## S2 — Safe `umol-nauty-sys` Rust surface

### S2a — Raw declarations and owned input/output types

**Module:** `umol-nauty-sys/src/lib.rs`  
**Kind:** additive (green)  
**Dependencies:** [dep: S1b]

Handwrite only the declarations for the umol shim, never nauty's structs. Add crate-private raw FFI
types and public/safe Rust-owned input and output types for CSR topology, ranked colors, canonical
labels, orbits, generators, and `(mantissa, exponent)`. Define a typed error corresponding to every
shim status. Validate:

- CSR offset length, monotonicity, terminal offset, and neighbor bounds;
- color length;
- vertex/directed-edge counts and all `usize`→C integer conversions;
- output lengths before entering C.

**Tests:** table-test valid and invalid input shapes and each conversion boundary reachable without
impractical allocations. Use specific error-variant assertions.

### S2b — Safe run wrapper and callback collector

**Module:** `umol-nauty-sys/src/lib.rs`  
**Kind:** additive (green)  
**Dependencies:** [dep: S2a]

Add the single safe `run` operation. It allocates output buffers, passes a scoped pointer to an owned
generator collector as the callback context, copies each callback permutation while valid, invokes
the shim, and returns owned results. Contain all `unsafe` in this crate and document the callback
lifetime and thread-safety invariants.

**Tests:** empty/singleton policy, same- and different-color edge, path, cycle, complete graph,
disconnected graph, canonical-label validity, orbit values, known group-order fields, generator
permutations, stabilizer coloring, error propagation, and parallel independent calls. Tests use
`rstest`; group semantics are asserted rather than exact generator lists except where the group has
one possible non-identity generator.

**Stage exit:** `cargo test -p umol-nauty-sys` is green; graph-core still uses the old dependency.

## S3 — Additive graph-core result vocabulary

This stage creates the final public types without yet changing `Graph::automorphisms`.

### S3a — `AutomorphismGroupOrder`

**Module:** `umol-graph-core/src/algorithms/auto.rs`, re-export in
`umol-graph-core/src/lib.rs`  
**Kind:** additive (green)  
**Dependencies:** none

Add:

```rust
pub enum AutomorphismGroupOrder {
    Exact(u128),
    Scientific { mantissa: f64, exponent: i32 },
}
```

Add checked construction from solver mantissa/exponent, display, and checked exact conversion. Use
`Exact` only when conversion is lossless; retain the normalized scientific pair otherwise. Do not
add a big-integer dependency.

**Tests:** table-test exact zero-exponent values, exactly representable positive exponents,
scientific fallback, overflow, invalid/non-finite mantissas if the constructor is fallible, display,
and exact-conversion behavior.

### S3b — `AutomorphismOutput`

**Module:** `umol-graph-core/src/algorithms/auto.rs`, re-export in
`umol-graph-core/src/lib.rs`  
**Kind:** additive (green)  
**Dependencies:** [dep: S3a]

Add the final owned result with fields `orbits`, `canonical_labels`, `node_count`, `orbit_count`,
`group_order`, and `generators`. Add the final query surface:

- `node_count`;
- `orbit_count`;
- `orbit_of` and `same_orbit`;
- `canonical_labels`;
- `group_order`;
- `generators`.

Keep the old `Automorphism` and `AutoGroupOrder` operational until S4. A crate-private constructor
allows S4 to assemble a validated output without exposing mutable fields.

**Tests:** table-test orbit queries and accessors through structurally meaningful outputs; avoid
tautological field-only assertions.

**Stage exit:** graph-core exports both old and new vocabularies; all existing callers and tests are
unchanged and green.

## S4 — Backend and API migration

This is the only intentionally breaking stage. It changes one apply surface,
`Graph::automorphisms`, then migrates every workspace caller before the stage ends.

### S4a — Vendored backend adapter beside the old backend

**Module:** `umol-graph-core/Cargo.toml`; `umol-graph-core/src/algorithms/auto.rs`  
**Kind:** additive (green)  
**Dependencies:** [dep: S0b, S2b, S3b]

Add the path dependency on `umol-nauty-sys`. Implement a private graph-core adapter that:

- ranks the caller's ordered colors exactly as today;
- converts `Graph` adjacency into the sys crate's CSR input without changing graph semantics;
- handles the empty graph in Rust;
- maps sys node indices, canonical labels, orbits, and generators into `NodeId`;
- maps mantissa/exponent into `AutomorphismGroupOrder`;
- returns `AutomorphismOutput`.

Leave the public dispatch on the old backend temporarily. Add differential conformance cases that
run old and new adapters on the S0 corpus and compare canonical graph/key semantics, orbits, group
order, and generated-group semantics—not raw generator vectors.

**Verification:** `cargo test -p umol-graph-core`; run the S0 benchmark corpus once with an internal
temporary selection mechanism if necessary, without exposing a second public Nauty variant.

### S4b — Rewire `Graph::automorphisms` and `canonical_key`

**Module:** `umol-graph-core/src/algorithms/auto.rs` and `src/lib.rs`  
**Kind:** breaking (red→green within S4)  
**Dependencies:** [dep: S4a]

Change `Graph::automorphisms` to return `AutomorphismOutput`, dispatch
`AutomorphismAlgorithm::Nauty` to the vendored adapter, and update `canonical_key` to call
`canonical_labels()`. Retire `Automorphism`, `AutoGroupOrder`, the old FFI imports, callback, and
thread-local accumulator. Update graph-core tests to the new type and method names, preserving the
S0 semantic cases.

The workspace is temporarily red until S4c migrates downstream callers.

### S4c — Migrate workspace consumers

**Module:** `umol-ast/src/ast/symmetry.rs`, `umol-ast/src/ast/view/graph.rs`, AST tests,
`umol-graph/src/fingerprint/substructure.rs`, and every remaining `Automorphism`/
`AutoGroupOrder`/`canonical_labeling` use found by `rg`  
**Kind:** breaking caller migration (red→green)  
**Dependencies:** [dep: S4b]

Replace type imports with `AutomorphismOutput`, method calls with `canonical_labels` and
`group_order`, and old group-order pattern matches with the new representation. Preserve stereo
generator grading, orbit fixpoints, stabilizer behavior, and canonical substructure keys. Do not
alter algorithm selection: callers continue to request `AutomorphismAlgorithm::Nauty`.

Update affected tests in definition order and under the test-writing conventions. Assert the same
domain behavior rather than wrapper implementation details.

### S4d — Remove the old dependency and differential scaffolding

**Module:** `umol-graph-core/Cargo.toml`, `Cargo.lock`, temporary differential-test code  
**Kind:** breaking cleanup (red→green)  
**Dependencies:** [dep: S4c]

Remove `nauty-Traces-sys`, its old adapter, and any temporary dual-backend selector. Keep the
backend-neutral S0 conformance suite and benchmarks. Regenerate the lockfile and verify that native
`bindgen` is no longer pulled through graph-core.

**Stage exit:** `cargo test --workspace` is green; `rg` finds no old public type/method names or
`nauty_Traces_sys`; `cargo tree -p umol-graph-core -i bindgen` reports no matching dependency.

## S5 — Performance and distribution acceptance

### S5a — Post-migration benchmark comparison

**Module:** benchmark results/process from S0; no production API change  
**Kind:** additive verification (green)  
**Dependencies:** [dep: S4d]

Run the identical S0 corpus against the vendored backend using the same target, profile, CPU, and
Criterion settings. Compare against the named old-binding baseline. Investigate material regressions
in CSR conversion, color ranking, callback collection, allocation, or nauty compile options. Record
the accepted comparison and the predeclared noise/regression threshold in doc 143 or an adjacent
benchmark note.

Also record workload-fitness latency for representative Python-facing molecule, incidence,
stabilizer, and canonical-key paths. This is the gate for declaring sparse nauty adequate; it is not
an algorithm comparison.

### S5b — LLVM-free package/build matrix

**Module:** workspace CI/package configuration and `umol-py` build documentation as needed  
**Kind:** additive verification (green)  
**Dependencies:** [dep: S4d]

Build/test on Linux, macOS, and Windows with a C compiler but without LLVM/libclang installed.
Exercise at least `cargo build -p umol-graph-core`, workspace tests, the maturin source-distribution
path, and wheel builds for supported Python targets. Confirm the vendored license is included in
source distributions as required and that no build step searches for libclang.

### S5c — Remove duplicate color sorting and remeasure

**Module:** `umol-graph-core/src/algorithms/auto.rs`, `umol-nauty-sys/src/lib.rs`,
`umol-nauty-sys/include/umol_nauty.h`, `umol-nauty-sys/src/umol_nauty.c`, and the S5 benchmark note

**Kind:** internal boundary optimization (red→green)

**Dependencies:** [dep: S5a]

Extend the private `NautyInput`/C-shim boundary with the vertex partition order already computed
while graph-core ranks generic Rust colors. Validate that this order is a permutation of all
vertices and is nondecreasing by ranked color. Initialize nauty's `lab` and `ptn` arrays directly
from that order, then remove `umol_colored_vertex`, its comparison function, and the C-side `qsort`.
Keep CSR validation, canonical-label, orbit, group-order, and generator behavior unchanged.

Update the sys-crate boundary tests for valid and invalid partition orders and retain graph-core's
backend-semantic conformance suite. Run formatting, Clippy, sys/core tests, and the workspace check
before benchmarking so the subitem ends green.

Rerun the identical 69-case S0/S5a Criterion corpus on the same host, compiler, profile, warm-up,
sample count, and measurement period. Save it under a new named baseline, compare it both with the
old-binding `base` dataset and the first `vendored-nauty` run, and update doc 143 with the complete
summary and representative molecule, incidence, stabilizer, and canonical-key latencies. Explicitly
report whether the four regressions above 10% shrink, persist, or move; do not declare the
performance gate accepted without a recorded acceptance threshold.

**Stage exit:** the performance gate is accepted; the supported native and Python build matrix is
green without LLVM/libclang; and the duplicate color sort is absent with the remeasurement recorded.
The core deliverable is complete.

## S6 — Traces backend (deferrable)

This entire stage is optional and is not needed when S5 shows sparse nauty meets workload targets.

### S6a — Traces native entry point

**Module:** `umol-nauty-sys` vendored sources, build script, C header/shim  
**Kind:** additive (green)  
**Dependencies:** [dep: S5c]

Only after a benchmark case justifies it, add the minimal Traces source closure and a separate
`umol_traces_run` entry point over the same umol-owned CSR/color/result convention. Keep
Traces-specific options, statistics, and callback translation internal. Add sys-crate conformance
tests parallel to the nauty cases.

### S6b — Public Traces selection and comparison

**Module:** `umol-graph-core/src/algorithms/auto.rs`, benchmarks, conformance tests  
**Kind:** additive (green)  
**Dependencies:** [dep: S6a]

Add `AutomorphismAlgorithm::Traces`, adapt its sys output to the existing
`AutomorphismOutput`, and run the same S0 benchmark/conformance corpus for both algorithms. Compare
canonical graph/key, orbit and group semantics, and performance; do not compare raw canonical-label
or generator arrays. Document that Traces supports only the undirected input currently in scope.

**Stage exit:** both backends are green and benchmarked; default selection changes only through a
separate design decision.

## Critical path

`S0a/S0b → S1a → S1b → S2a → S2b → S3a → S3b → S4a → S4b → S4c → S4d → S5a → S5c`.

S5b can proceed in parallel after S4d and joins S5c at the S5 stage exit.

S1 and S3 may proceed in parallel after S0 because the native foundation and graph-core vocabulary
are independent. S4a is the join point. S4b–S4d must land as one stage/green commit series because
the return-type rename and caller migration are inseparable.

## Deferrable work

- **S6 (Traces):** entirely deferrable; the core distribution objective is complete after S5.
- Exact arbitrary-precision group order is not staged; `u128` plus scientific fallback fulfills the
  settled contract.
- A pure-Rust individualization/refinement engine is outside this plan.
- Further allocation reduction or caching is deferred unless the S5 benchmark comparison identifies
  it as necessary for migration parity.
