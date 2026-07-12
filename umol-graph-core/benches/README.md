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

