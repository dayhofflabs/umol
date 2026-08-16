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
upstream license, and verify the packaged `umol-msym-sys` crate rather than
only the workspace checkout.

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
`umol-nauty-sys/src`, `umol-nauty-sys/include`, and the build script rather
than in the vendored directory. An upstream refresh should preserve that
separation and verify the packaged `umol-nauty-sys` crate.
