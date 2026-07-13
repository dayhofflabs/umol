# Graph-core benchmark baselines

The automorphism migration corpus lives in `algorithms.rs` under three Criterion groups:

- `automorphism` — ordinary uniform/degree/unique colorings plus incidence graphs;
- `automorphism_stabilizer` — three uniquely colored site stabilizers per iteration;
- `canonical_key` — edge-subdivided, node/edge-colored canonical keys.

Save the pre-migration baseline with:

```sh
cargo bench -p umol-graph-core --bench algorithms -- automorphism --save-baseline nauty-traces-sys
cargo bench -p umol-graph-core --bench algorithms -- canonical_key --save-baseline nauty-traces-sys
```

After replacing the binding, compare on the same machine and toolchain with:

```sh
cargo bench -p umol-graph-core --bench algorithms -- automorphism --baseline nauty-traces-sys
cargo bench -p umol-graph-core --bench algorithms -- canonical_key --baseline nauty-traces-sys
```

Record `rustc -vV`, the operating-system/CPU description, and any non-default build flags alongside
durable comparison results. Criterion's machine-specific raw results remain local and are not
committed.

## Initial binding

- profile: workspace `bench` profile (`release` optimization plus debug information);
- binding: `nauty-Traces-sys` 0.11.0;
- bundled solver: nauty 2.9.3, sparse undirected entry point (`sparsenauty`);
- definitions: `USE_TLS`; `WORDSIZE=64` on this 64-bit target; no `native`, `popcnt`, or `lzc`
  feature;
- baseline host: `aarch64-apple-darwin`, Darwin 24.6.0;
- baseline compiler: `rustc 1.98.0-nightly (b354133fb 2026-06-03)`.

## Matching enumeration

The output-sensitive matching corpus is part of the single-file `algorithms` benchmark suite.
It covers molecular
fixtures (benzene, naphthalene, coronene, azulene, and C60), disconnected components, a ladder,
an odd grid using maximum matching, a prescribed-hole residual path, and dense bipartite K6,6.
Fixture parsing, graph construction, embedding validation, and reference output counts happen
before Criterion's timed closures.

The groups distinguish:

- `matching_first_output`: time to the first visitor result;
- `matching_visit_prefix`: visitor prefixes of 1, 10, and 100 available outputs;
- `matching_visit_full`: full streaming traversal retaining only a scalar count;
- `matching_eager_collection`: full traversal including collection allocation;
- `matching_fkt_count`: independent planar perfect-matching counts.

The full visitor and eager groups label throughput in enumerated outputs. A short smoke
measurement can be run with:

```sh
cargo bench -p umol-graph-core --bench algorithms -- matching_first_output/benzene --sample-size 10 --measurement-time 1
```

Inter-output delay diagnostics are deliberately outside Criterion's timed loops. Enable a single
diagnostic traversal per corpus case with:

```sh
UMOL_MATCHING_DELAY_DIAGNOSTICS=1 cargo bench -p umol-graph-core --bench algorithms -- --test
```

The diagnostic TSV reports the observed first, median, 95th-percentile, and maximum delay.
These scheduler-sensitive observations and Criterion's machine-specific raw results remain local
and are not committed.
