# Vendored third-party code

The native sources below are tracked as ordinary repository files so that
repository archives, Cargo packages, and Python source distributions are
self-contained. Upstream code and licenses remain in their respective system
crates; umol-specific interfaces are kept outside the upstream source trees.

## libmsym

- **Location:** `umol-msym-sys/libmsym`
- **Upstream:** <https://github.com/mcodev31/libmsym>
- **Revision:** `85e47232376a8e735c2a7b5283f40b59b3953db1`
- **Revision date:** 2024-06-15
- **Description:** the complete upstream tree at the recorded revision
- **License:** MIT; retained in `umol-msym-sys/libmsym/LICENSE`

`umol-msym-sys/build.rs` selects and compiles the C source files used by the
Rust wrapper. To refresh libmsym, replace the directory with a clean checkout
of the chosen upstream revision, update the revision above, retain the
upstream license, and verify the packaged `umol-msym-sys` crate. A successful
workspace checkout alone does not establish package completeness.

## nauty

- **Location:** `umol-nauty-sys/nauty`
- **Upstream release:** nauty and Traces 2.9.3
- **Upstream sites:** <https://users.cecs.anu.edu.au/~bdm/nauty/> and
  <https://pallini.di.uniroma1.it/>
- **Imported into umol:** commit
  `8524073db96fe6587fba06e68639ec35aad9a54d` on 2026-07-11
- **Description:** the upstream source closure required by `sparsenauty`, not
  the complete nauty and Traces distribution
- **License:** Apache-2.0; retained in `umol-nauty-sys/nauty/COPYRIGHT` and
  `umol-nauty-sys/nauty/LICENSE-APACHE`

`umol-nauty-sys/build.rs` records the selected upstream source and header
files. The stable umol interface and portability configuration live in
`umol-nauty-sys/src`, `umol-nauty-sys/include`, and the build script. The
vendored directory contains the selected upstream source closure. An upstream
refresh should preserve that separation and verify the packaged
`umol-nauty-sys` crate.

## CoordGen

- **Location:** `umol-coordgen-sys/coordgen`
- **Upstream release:** CoordGen 3.0.2
- **Upstream:** <https://github.com/schrodinger/coordgenlibs>
- **Revision:** `c4dd5b0e1f1971c06c7ab85725c185e47211814e`
- **Revision date:** 2023-02-01
- **Description:** the library source and header closure, including generated
  built-in templates; tests, examples, and optional MAE-parser support are not
  included
- **License:** BSD-3-Clause; retained in
  `umol-coordgen-sys/coordgen/LICENSE`

`umol-coordgen-sys` compiles this snapshot only when its `native` feature is
enabled. That path requires a C++11 compiler; the default feature set does not
invoke a C++ compiler. The focused native boundary is verified with
`cargo test -p umol-coordgen-sys --features native`. To refresh CoordGen,
replace the vendored source closure from a clean upstream revision, retain the
upstream license, update the revision above, and verify both the feature-gated
crate and its packaged contents.
