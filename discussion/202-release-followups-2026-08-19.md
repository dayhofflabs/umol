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

## Whitepaper link

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
