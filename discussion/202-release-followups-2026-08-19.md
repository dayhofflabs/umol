# 202 — Release follow-ups

Status: In Progress
Date: 2026-08-31
Relates: [163](163-release-preparation-2026-07-26.md),
[168](168-api-hygiene-2026-07-27.md)

## Scope

This document tracks operational work intentionally deferred after the successful 0.6.0 release.
It does not reopen the released API, packaging metadata, or publication graph recorded in doc 163.

## Trusted publishing **Done**

The first release after 0.6.0 should exercise the crates.io OIDC path now configured in
`release.yml`. All eighteen existing crates have a trusted publisher for repository
`dayhofflabs/umol`, workflow `release.yml`, and environment `crates-io`; the bootstrap token has
been revoked and the `CARGO_REGISTRY_TOKEN` was removed from the `crates-io` environment.

## 0.8.0 release preparation

The workspace's nineteen packages and internal dependency requirements are set to 0.8.0.
RELEASE_NOTES.md replaces the previous notes with changes since 0.7.0, including breaking
API migrations. The README adds a verified primer and file-based SVG output.

The new umol-coordgen-sys package includes its MIT/Apache license texts alongside the vendored
CoordGen BSD license. Native tests pass (17 tests), and `cargo package --locked --offline
-p umol-coordgen-sys --features native --allow-dirty` successfully builds the extracted package.
The archive contains 50 files, approximately 132 KiB compressed. The local Python extension
build and its suite pass (1,479 tests, two skipped); workspace checking also passed.

### First CoordGen publication remains manual

crates.io trusted publishing cannot perform a crate's first publication; see the
[crates.io development update](https://blog.rust-lang.org/2025/07/11/crates-io-development-update-2025-07/).
After reviewing and committing the prepared tree, run locally with normal Cargo credentials:

```sh
cargo publish --locked -p umol-coordgen-sys --features native --dry-run
cargo publish --locked -p umol-coordgen-sys --features native
```

Then configure that crate's trusted publisher for repository `dayhofflabs/umol`, workflow
`release.yml`, environment `crates-io`. The release workflow already lists umol-coordgen-sys
before its consumers and skips versions already published. No workflow change is needed for
this bootstrap. The normal release dispatch requires the reviewed `v0.8.0` tag to exist.

No crate was uploaded, no publisher settings were changed, and no commit or tag was created
during this preparation. The dirty-tree allowance above was for local package verification only.

## Whitepaper link

The Python primer and Rust appendix have been updated for 0.8.0. The README contains the
adapted tour; the paper source remains author-maintained.

The README currently links the bundled `docs/umol-whitepaper.pdf`. Replace that target with the
permanent arXiv page once the paper is published; keep the bundled PDF available until then.

## Published Rust documentation

Review the generated docs.rs pages for all published workspace crates from a new user's point of
view. Add concise crate- and module-level background where the generated API listings do not explain
the crate's responsibility, its relationship to the other umol models, or the intended entry points.
Clean up stale descriptions, weak navigation, broken or unclear cross-links, and public items whose
rendered documentation does not state enough context to use them correctly.

Treat the published pages as the review surface, while editing rustdoc in the source. Documentation
corrections stay in this follow-up; findings that require changing visibility, re-exports, module
boundaries, or public names belong in doc 168.

## Python documentation site

Set up and publish an online documentation site for `umol-py`. Select Sphinx or another maintained
generator after checking how well it documents the compiled PyO3 surface; the choice should support
runtime API introspection, cross-linked reference pages, search, and ordinary prose guides without
requiring a parallel hand-maintained signature inventory.

The initial site should include installation and first-use guidance, model background, the main
molecule and reaction workflows, and a navigable API reference. Add an automated documentation
build and deployment path, then update the package's `Documentation` URL and repository links to the
published site.

## Completion criteria

- A release publishes the Rust crates through OIDC without the bootstrap credential.
- The obsolete GitHub secret is deleted and the trusted-publishing-only decision is recorded.
- Public-runner CI timing is measured and any justified timeout or cache change is applied.
- The README points to the permanent arXiv page when it exists.
- The docs.rs pages have been reviewed and the source rustdoc supplies the missing orientation,
  navigation, and public-API context.
- The `umol-py` documentation site is generated from maintained prose and the live extension API,
  deployed automatically, and linked from the package metadata.
