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

## Registry name availability (checked 2026-08-02)

Availability is time-sensitive; re-check immediately before publishing.

**crates.io — `umol` is claimed; the other eighteen are free.**

`umol` was published as a `0.0.0` placeholder on 2026-08-02 (https://crates.io/crates/umol): three
files, no functionality, no dependencies, a README stating that it is a name reservation for a
library under active development. It commits no API, so it does not constrain doc
[176](176-ast-naming-2026-07-31.md); the naming decision remains open until the real 0.6.0 ships.
The placeholder source is at `~/Source/rust/umol-placeholder`, outside the workspace.

Note for the real release: a version number can never be reused, so 0.6.0 must be published fresh;
`0.0.0 -> 0.6.0` is a valid increase.

The remaining eighteen were free as of the same date: `umol-ast`, `umol-ast-macros`, `umol-chem`,
`umol-edn`, `umol-edn-macros`, `umol-geometric`, `umol-geometric-core`, `umol-geometric-graph`,
`umol-graph`, `umol-graph-core`, `umol-io`, `umol-msym`, `umol-msym-sys`, `umol-nauty-sys`,
`umol-params`, `umol-perm`, `umol-py`, `umol-utils`. Each returned HTTP 404 from `https://crates.io/api/v1/crates/<name>`. Nothing blocks the Rust side.
These carry far less squatting risk than the bare four-letter name did, so they can wait for the real
release rather than being reserved individually — a placeholder for a crate that never ships is the
practice this project should avoid.

**PyPI — `umol` is taken, and this blocks the plan above as written.**

- One release, `0.1.0`, carrying **zero files**. `pip install umol` cannot succeed; there is no
  artifact to install.
- Owner listed as Steven (Yuhang) Wang, summary "umol: molecular data analytics suite".
- Linked repository `github.com/UMOL/umol-py` was created 2017-06-14 and last pushed 2017-06-22 —
  eight days of activity, nine years ago, zero stars, not archived.

**Decision (author, 2026-08-02): publish as `umol-py`.** Do not pursue the name.

**The distribution name and the import name are independent**, which is why this costs almost
nothing: `scikit-learn` imports as `sklearn`, `pillow` as `PIL`, `attrs` as `attr`. Every listing in
the whitepaper and the documentation continues to read `from umol import ...` whatever the wheel is
called. The whitepaper's Availability section names the package by how it is imported rather than
where it is published, so it is already correct and needs no revision.

**Why not PEP 541**, recorded so the option is not reopened without new information. The policy
declares a project abandoned only when *all three* hold: owner unreachable, no releases in twelve
months, no activity on the home page. PyPI attempts contact at least three times and stops after six
weeks. Requests are filed at `github.com/pypi/support` after the requester has tried the owner
directly by email and repository issue.

Two tracks exist, and ours is the harder one. *Continuation* — taking over the same project — is
lenient. *Reuse* — a different project claiming the name — additionally requires meeting notability
requirements and showing download statistics prove the incumbent is unused. \umol\ is
unambiguously reuse.

The case would probably succeed: no release since 2017, home page untouched since 2017-06-22, and a
package with zero files has no downloads to speak of. But six weeks of mandated waiting precede any
consideration, on a volunteer-staffed tracker, for a discretionary outcome — realistically months —
in exchange for a purely cosmetic gain. Reconsider only if `umol` on PyPI becomes a genuine source of
user confusion.

**Do not plan a later migration to `umol`, and if one ever happens, obey one rule.** The classic
painful renames (`BeautifulSoup` to `bs4`, `pycrypto` to `pycryptodome`) hurt because the *import*
name changed and every dependent's source had to change. That does not apply here: `import umol` is
correct from the first release under any distribution name, so the distribution name appears in a
`pip install` line and nowhere else. Distribution and import names differing is ordinary — `PyYAML`
imports as `yaml`, `python-dateutil` as `dateutil`, `Pillow` as `PIL`, `attrs` as `attr`.

If a switch is ever made, the mechanism is a shim: the real code moves to the new distribution and
the old one remains as a metadata-only package whose entire content is a dependency on the new one.
Existing installs and pins keep working; releases are never deleted, and yanking is the advisory tool
if one must be withdrawn.

**The rule: the old distribution must ship zero modules after a switch.** If both distributions
contain code, both install a `umol/` directory into the same `site-packages`, and which one wins
depends on install order. That is the one failure mode in this area that is genuinely hard to
diagnose.

## Rust Additional Steps

1. The repo needs a README.md document with the Getting Started section (corresponds to the Primer section of the whitepaper).
2. Need to carefully consider if umol-ast and umol-io dependencies should be re-exported from umol-graph.
3. The workspace needs to set version.workspace = true in individual crates and [workspace.package] version = "0.6.0" in the top-level Cargo.toml.
4. Need to check which other fields should be set in the Cargo.toml files.
5. Check if all crates need to be published now, umol-geometric*, umol-msym* are not required by the graph infrastructure.

## Python Additional Steps

1. Check which additional fields need to be added to the pyproject.toml file.
2. CI/CD pipeline setup for building Python wheels (linux, macos-arm).
