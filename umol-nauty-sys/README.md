# umol-nauty-sys

Vendored native foundation for umol's sparse nauty integration.

This crate compiles a minimal vendored subset of upstream nauty and exposes a
small, safe Rust boundary (`NautyInput`, `WeightedNautyInput`, `NautyOutput`,
`NautyGroupOrder`, `NautyError`, `run`, `run_weighted`, and `NAUTY_VERSION`).
Upstream nauty data structures, options, statistics, and dynamic-allocation
macros remain private to the C layer.

## Vendored upstream nauty (`nauty/`)

- **Upstream release:** nauty and Traces `2.9.3` (Brendan McKay and Adolfo
  Piperno; <https://users.cecs.anu.edu.au/~bdm/nauty/>,
  <https://pallini.di.uniroma1.it/>).
- **License:** Apache-2.0; retained in [`nauty/COPYRIGHT`](nauty/COPYRIGHT) and
  [`nauty/LICENSE-APACHE`](nauty/LICENSE-APACHE).
- **Minimal source closure:** only the files required by sparse nauty are
  vendored, not the full nauty and Traces distribution:
  - C sources: `nauty.c`, `nautil.c`, `nausparse.c`, `schreier.c`, `naurng.c`
  - Headers and templates: `nauty.h`, `nausparse.h`, `schreier.h`, `naurng.h`,
    `naututil.h`, `sorttemplates.c`
  - Upstream metadata and licenses: `COPYRIGHT`, `LICENSE-APACHE`, `README`,
    `This_is_nauty2_9_3.txt`

Every file under `nauty/` is kept unmodified from upstream. All umol-authored C
and Rust code lives in `include/`, `src/`, and `build.rs`.

## Build configuration and symbol namespacing (`build.rs`)

[`build.rs`](build.rs) compiles two static C archives via `cc`:

1. `umol_nauty` — the five vendored upstream C files (C99).
2. `umol_nauty_shim` — [`src/umol_nauty.c`](src/umol_nauty.c) and
   [`src/umol_nauty_swg.c`](src/umol_nauty_swg.c) (C11, `-Wall -Wextra`).

Both builds define:

- `USE_TLS` for thread-local workspace isolation across concurrent threads.
- `WORDSIZE=32` on Windows or 32-bit targets and `WORDSIZE=64` otherwise.
- Compile-time symbol prefixing (`-D<symbol>=umol_nauty_<symbol>`) for all 71
  global symbols exported by the vendored upstream source closure
  (`PREFIXED_SYMBOLS` in [`build.rs`](build.rs)), preventing link-time
  collisions if another copy of nauty is linked into the same process.

The shim archive exports only `umol_`-prefixed symbols (`umol_nauty_run`,
`umol_nauty_run_weighted`, `umol_dispatch_swg`, and
`umol_sparsenauty_weighted`).

## C boundary and weighted sparse dispatch (`include/`, `src/`)

- [`include/umol_nauty.h`](include/umol_nauty.h) and
  [`src/umol_nauty.c`](src/umol_nauty.c) implement the C boundary used by
  [`src/lib.rs`](src/lib.rs). They validate CSR offsets, vertex/degree bounds,
  vertex partitions, and edge weights (`0..INT_MAX`), reject same-thread
  reentrant calls, manage thread-local generator callbacks, and dispatch both
  unweighted (`umol_nauty_run`) and edge-weighted (`umol_nauty_run_weighted`,
  with `bool get_canon` to toggle canonical-labeling work) runs.
- [`src/umol_nauty_swg.h`](src/umol_nauty_swg.h) and
  [`src/umol_nauty_swg.c`](src/umol_nauty_swg.c) implement a weighted sparse
  dispatch vector (`umol_dispatch_swg`) and entry point
  (`umol_sparsenauty_weighted`) adapted from upstream `nausparse.{c,h}` to
  operate directly on non-negative edge weights in `sparsegraph.w` without
  modifying `nauty/`. They provide weight-aware `isautom_swg`,
  `testcanlab_swg`, `updatecan_swg`, `init_swg`, `refine_swg`,
  `targetcell_swg`, and conservative `cheapautom_swg` (`k <= nnt + 1`) for
  simple undirected edge-colored graphs.

## License

Licensed under the Apache License, Version 2.0 (see
[`LICENSE-APACHE`](LICENSE-APACHE), [`nauty/COPYRIGHT`](nauty/COPYRIGHT), and
[`nauty/LICENSE-APACHE`](nauty/LICENSE-APACHE)).
