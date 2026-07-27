# Praparations for Rust and Python package release

Status: **Proposed**
Date: 2026-07-26
Relates: [151](151-python-molecule-workflows-2026-07-13.md),

## Scope

This document covers verification and release steps for the following Rust crates:

- `umol-ast`
- `umol-ast-macros`
- `umol-chem`
- `umol-edn`
- `umol-edn-macros`
- `umol-geometric`
- `umol-geometric-core`
- `umol-geometric-graph`
- `umol-graph`
- `umol-graph-core`
- `umol-io`
- `umol-msym`
- `umol-msym-sys`
- `umol-nauty-sys`
- `umol-params`
- `umol-perm`
- `umol-py`
- `umol-utils`

The crates should be tagged as revision 0.6.0, this version reflects prior internal iteration. The CI/CD pipeline should be set up. The python package should be prepared and published to the pypi server.

## Rust Additional Steps

1. The repo needs a README.md document with the Getting Started section (corresponds to the Primer section of the whitepaper).
2. Need to carefully consider if umol-ast and umol-io dependencies should be re-exported from umol-graph.
3. The workspace needs to set version.workspace = true in individual crates and [workspace.package] version = "0.6.0" in the top-level Cargo.toml.
4. Need to check which other fields should be set in the Cargo.toml files.
5. Check if all crates need to be published now, umol-geometric*, umol-msym* are not required by the graph infrastructure.

## Python Additional Steps

1. Check which additional fields need to be added to the pyproject.toml file.
2. CI/CD pipeline setup for building Python wheels (linux, macos-arm).
