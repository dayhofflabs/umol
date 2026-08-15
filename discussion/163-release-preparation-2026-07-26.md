# Praparations for Rust and Python package release

Status: **In Progress** — guiding document for the 0.6.0 release
Date: 2026-07-26
Relates: [151](151-python-molecule-workflows-2026-07-13.md),

## Scope

This document covers verification and release steps for the workspace crates
(names as of the current tree; the former `umol-ast`/`umol-ast-macros` are
`umol-graph-ir`/`umol-graph-ir-macros`):

- `umol-graph-ir`
- `umol-graph-ir-macros`
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

## Semantic versioning policy (0.6 line)

Cargo's compatibility convention for pre-1.0 versions applies: the leftmost
nonzero component is the breaking boundary, so `0.6.x` releases must be
compatible with `0.6.0`, and any breaking change ships as `0.7.0`. Consequences:

- Conformance-suite growth, atom-typing registry rows, and living
  default-valence-table additions are compatible (`0.6.x`); resolution outcomes
  may change through added candidate rows — that is the documented update
  policy (candidate sets are the monotone observable), state it in the release
  notes. Frozen preset tables never change in place; a preset revision is a new
  named file and at least a minor bump.
- Internal cleanup that does not move public paths is `0.6.x`. Module
  restructuring that changes public paths — including the docs 164–168 API
  worklists — is `0.7.0` material, not a patch.
- `umol-py` follows the same version as the workspace; the wheel and the crates
  are tagged together.

## Functional-gap assessment for 0.6 (2026-08-15)

The whitepaper is the functionality specification; every concrete claim in its
listings has been verified against the live code (doc 194 S6e6). Two code-level
items must land before the tag:

1. `AromaticityRule::Clar { ring_limits }` advertises a payload that
   `ClarAromaticity` never reads (the ring request is hardcoded at six).
   Advertised-but-ignored configuration is exactly the "behaves unexpectedly"
   class; either the payload is read or it is removed — removal is a type
   change, so it must precede the tag.
   **Done 2026-08-15:** removed — no `RingLimits` field applies to Clar (the
   sextet size is inherent; the union limits govern the Hückel mechanism);
   `Clar` is a unit variant in Rust and a no-argument constructor in Python.
2. `umol-py/pyproject.toml` still declares `name = "umol"`; the 2026-08-02
   decision below publishes as `umol-py` (import name stays `umol`).
   **Done 2026-08-15:** renamed to `umol-py`; the wheel installs as
   `umol-py`, the module imports as `umol`. The dormant PyPI `umol` takeover
   stays unpursued per the PEP 541 analysis below; revisit in a year if it
   ever becomes a source of confusion.

Assessed and **not** blocking, with the reason recorded:

- **Doc 195 (molecule-scope constraint matching)**: a pattern carrying a
  molecule-scope constraint is refused loudly with an error naming the
  construct (the doc 194 S1a gate) instead of matching as if unconstrained.
  The whitepaper neither documents molecule-scope constraints in patterns nor
  promises their evaluation, so this is a documented limitation for the
  release notes, not a gap against the specification. Same for the gated
  constraint-remove reaction path (six `#[ignore]`d tests).
- **Doc 193 (recursive subpattern constraints)**: unimplemented proposal; not
  in the whitepaper.
- **Doc 149 Part B (hashing tiers)**: the whitepaper presents hashability as a
  consequence of canonical labeling and demonstrates `canonical_eq` only;
  `canonical_hash` stays deferred until a consumer (network dedup) exists.
- **Docs 164–168 (cleanup worklists)**: no blockers. The two once-functional
  items in doc 166 are resolved or out of scope — the validator completion
  landed with the doc 194 S6a family (spin invariants included), and the
  implicit/explicit hydrogen transformer pair is unimplemented but not
  whitepaper-promised (a good `0.6.x` addition). The API reshapes in these
  docs are `0.7.0` material under the versioning policy.
- **Doc 56 (registry rows)** and conformance growth: additive data, `0.6.x`.
- **Doc 175 (`Ground<T>` API)**: unimplemented; likely superseded by the
  concreteness API (`is_concrete`/`into_concrete`, doc 194 S6e5) — review and
  close, no release impact.
- **Doc 180 (facade crate)**: the whitepaper's Rust primer instructs depending
  on `umol-graph`/`umol-graph-ir`/`umol-graph-core` directly, so 0.6 ships
  without a facade; the crates.io `umol` placeholder stays as the name
  reservation. A later facade that only re-exports is additive (`0.6.x` or
  `0.7`); closing off the individual crates behind it would be breaking and
  runs against the multi-model composition philosophy — the individual crates
  stay public regardless.
- **Doc 182 (Python resolution exposure)**: reviewed 2026-08-15 — the gap is
  real and current. Python has no `resolve` method anywhere: a molecule built
  by editing or `parse` cannot be resolved except by round-tripping through
  ingest, so the refine loop (assert, re-resolve) is not executable from
  Python. Not a tag blocker — the whitepaper shows no `resolve()` listing,
  and adding a method is additive (`0.6.x`) — but it is the most user-visible
  API hole in this assessment, and it carries doc 182's open design question
  (verdict value versus exception) to settle before implementing. The doc's
  Rust reference section predates the S4 rework (`Resolver::with_config`,
  `ResolverError` spellings) and needs refreshing against `ResolveConfig`
  stored on `Resolver` and the `Resolve*` names.
- Behaviors to state in the release notes rather than change: higher stereo
  kinds (allene, square-planar, octahedral) are staged off by the default
  `StereoModel`; the whitepaper demonstrates tetrahedral and cis/trans only.

## Release collateral (missing items)

1. **README.md** at the repository root — Getting Started from the whitepaper
   primer (step already listed below); include the `pip install umol-py` /
   `import umol` distinction.
2. **Release notes for 0.6.0** — the whitepaper feature set, the known
   limitations above, and the data-update policy statement.
3. **License files** — `pyproject.toml` declares `MIT OR Apache-2.0`, but the
   repository has no `LICENSE-MIT`/`LICENSE-APACHE` files and the crate
   manifests carry no `license` field; both are required for crates.io.
4. **CI** — no `.github/` exists. Needed: a test workflow (workspace build,
   `--all-features --tests` so the conformance and proptest suites run,
   clippy, fmt) and a release workflow (crates.io publish in dependency
   order; maturin wheel builds for linux x86_64 and macos-arm — the binding
   uses abi3, so one wheel per platform — plus sdist, publish to PyPI).
5. **Workspace version** — `[workspace.package] version = "0.6.0"` plus
   `version.workspace = true` in member crates (step 3 below), and the
   remaining Cargo metadata fields (description, license, repository,
   keywords, readme) per crate.
6. **Name availability re-check** — the 2026-08-02 check covered the old
   `umol-ast`/`umol-ast-macros` names; `umol-graph-ir` and
   `umol-graph-ir-macros` have not been checked. Re-check all names
   immediately before publishing.

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
2. Need to carefully consider if umol-graph-ir and umol-io dependencies should be re-exported from umol-graph.
   Standing state (2026-08-15): `umol-graph`'s root is pub modules only with zero
   re-exports, and the whitepaper's listings use module paths (`ops::resolve::Resolver`,
   `ingest::ingest_smiles`) consistently with it — that is the de-facto export policy
   and single types get no one-off root exceptions.
   **Decided (author, 2026-08-15): no cross-crate re-exports**, proc macros excepted.
   Consumers depend on the defining crates directly; the whitepaper's three-crate
   dependency block is the documented experience.
   **Still open: the per-crate root-export policy** — what each crate's `lib.rs`
   re-exports of its own modules. Today `umol-graph` and `umol-graph-ir` are
   pub-modules-only (paths like `ops::resolve::Resolver`, `ir::Molecule`;
   `#[macro_export]` macros land at the root by the mechanism). Whether that stays the
   uniform rule or roots curate a selection is a doc 168 api-hygiene decision; settle
   it before the README teaches the import style.
3. The workspace needs to set version.workspace = true in individual crates and [workspace.package] version = "0.6.0" in the top-level Cargo.toml.
4. Need to check which other fields should be set in the Cargo.toml files.
5. Check if all crates need to be published now, umol-geometric*, umol-msym* are not required by the graph infrastructure.

## Python Additional Steps

1. Check which additional fields need to be added to the pyproject.toml file.
2. CI/CD pipeline setup for building Python wheels (linux, macos-arm).
