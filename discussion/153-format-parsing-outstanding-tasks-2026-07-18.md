# SMILES / MOL / SDF parsing outstanding tasks

Status: **Active task inventory**
Date: 2026-07-18
Relates: [047](047-smiles-conformance-suite-2026-01-21.md), [048](048-smiles-parser-configuration-2026-01-23.md), [100](100-table-ir-raise-ast-2026-05-27.md), [112](112-ctfile-winnow-migration-2026-06-13.md), [151](151-python-molecule-workflows-2026-07-13.md), [152](152-basic-molecule-wildcards-2026-07-18.md)

## Purpose

This document collects the known outstanding tasks around SMILES, MOL, and SDF parsing after the OpenSMILES wildcard work in doc 152. It is an inventory, not a staged implementation plan. Numbering is for stable reference only and does not imply implementation order.

The shared direction is:

- external formats should have explicit boundary objects rather than leaking parser internals into graph-level APIs;
- ordinary SMILES and CXSMILES should be separated at the format-boundary level;
- MOL and SDF need the same boundary cleanup, but their parser implementation is also due for a coordinated migration;
- parser result types should move toward semantically named TableIR structures instead of basic/extended variants chosen for parser convenience.

## Closed baseline

The following items are already handled by doc 152 and should not be reopened in this workstream:

- `*` is part of ordinary OpenSMILES parsing.
- `Molecule` can represent wildcard atoms.
- SMILES wildcard atoms raise to `ElementAst::Undetermined`.
- The old `WILDCARDS` and `BASIC_OPENSMILES` parse-flag split has been removed.
- `SmilesIoConfig::basic_opensmiles()` has been removed.
- Basic/OpenSMILES conformance classification has been collapsed.
- `EXTENDED_AROMATICS` and `EXTENDED_BONDS` remain ordinary SMILES acceptance-policy flags. They are not CXSMILES-specific merely because they previously lived near extended parser paths.

Doc 151 also records a benchmark gate that rejected direct replacement of compact `Molecule` by current `ExtendedMolecule`: direct replacement has unacceptable size and parse-time cost. The remaining design target is therefore a compact semantic superset, likely with cold extension records or side tables.

## Outstanding tasks

### T0 — Separate ordinary SMILES from CXSMILES

The current `SmilesIoConfig` still contains CX/ChemAxon-specific flags and presets. These should move behind a separate future CXSMILES boundary.

Required cleanup:

- Introduce or plan an explicit `CxSmiles` boundary type separate from ordinary `Smiles`.
- Remove CX-specific flags from ordinary `SmilesIoConfig`:
  - `CHEMAXON_EXTENSIONS`;
  - `SKIP_UNKNOWN_CHEMAXON_TAGS`;
  - `CHEMAXON`;
  - `SmilesIoConfig::chemaxon()`.
- Preserve `EXTENDED_AROMATICS` and `EXTENDED_BONDS` as ordinary SMILES acceptance-policy flags.
- Decide the public fate of `parse_extended_smiles*` and `parse_extended_reaction_smiles*` once the CXSMILES boundary exists.
- Update diagnostic and conformance tooling that still reports categories such as `basic_chemaxon` and `chemaxon`.
- Define conversion semantics:
  - `Smiles -> CxSmiles` should be lossless;
  - `CxSmiles -> Smiles` should be fallible when CX-only payload is present.

### T1 — Define CXSMILES semantic payload ownership

The CX parser already has payloads that exceed ordinary `Molecule`, but the semantic ownership is still unresolved.

Required decisions and work:

- Inventory which CXSMILES entries can be represented by compact `Molecule` plus existing AST raise semantics and which require extended payload storage.
- Resolve the current CX TODOs around fragment groups, stereo groups, relative stereo, ligand order, and related molecule-level payloads.
- Define how CX coordination and multicenter constructs map into TableIR and AST semantics.
- Define coordinate/property ownership for CX payloads and whether these are shared with MOL/SDF records.
- Add CXSMILES tests, conformance fixtures, fuzz coverage, and rendering/round-trip cases after the boundary is explicit.

### T2 — Add MOL and SDF boundary wrappers

SMILES now has an explicit parsed-format wrapper direction. MOL and SDF need the same treatment instead of graph-level APIs directly exposing parser helper names.

Required work:

- Add explicit `Mol` and `Sdf` boundary objects analogous to `Smiles`.
- Keep syntax parsing separate from model ingestion.
- Define ingestion APIs from `Mol`/`Sdf` into graph-level molecule types.
- Decide whether SDF exposes streaming/iterator-style records, eager record vectors, or both.
- Define render APIs and per-operation configs for MOL and SDF.
- Retire, rename, or redirect old helper APIs such as direct `parse_mol*_to_table_ir*` entry points after wrappers exist.

### T3 — Replace basic/extended TableIR parser result split

The current `Molecule` / `ExtendedMolecule` split was introduced for parser efficiency and now leaks an implementation concern into semantic design.

Required work:

- Design a compact semantic superset for molecule records that does not carry the current `ExtendedMolecule` size and latency penalties on the common path.
- Split `ExtendedAtom`, `ExtendedBond`, SGroup, RGroup, and related structures into semantically named records instead of parser-convenience containers.
- Decide whether cold extension records, side tables, or another compact layout should hold rarely used payloads.
- Remove parser result type choices where possible so each external format has one parsed boundary representation.
- Review the public localized-bond endpoint accessor on TableIR. Keep the
  endpoint pair directly accessible if that is the stable semantic shape;
  avoid wrapping it only for parser convenience.
- Add whole-record MOL/SDF benchmarks before replacing representation internals.
- Keep parse latency, allocation count, and retained size as explicit gates.

### T4 — Relocate or redesign `ChiralityFrame`

`ChiralityFrame` is still a molecule-level TableIR field even though it reflects source-format stereo conventions.

Required work:

- Redesign `ChiralityFrame` so raw stereo descriptors cannot be interpreted without the required source frame.
- Avoid polluting ordinary test construction with source-format chirality details when tests do not care about raw descriptor frames.
- Preserve the distinction between SMILES `FirstNeighborToward` and CTFile `LastNeighborAway` semantics.
- Add targeted tests across SMILES, MOL, SDF, and eventually CX enhanced stereo.
- Coordinate the redesign with TableIR-to-AST raise semantics.

### T5 — Close TableIR-to-AST raise gaps for parsed formats

Several parsed-format features still need semantic audits or stronger tests before they can be considered stable model-ingestion behavior.

Required work:

- Audit CTAB `vvv` mapping into AST valence constraints.
- Represent query-MOL `hhh >= 2` as minimum-hydrogen pattern semantics rather than a literal hydrogen count.
- Add stronger corpus/test pressure for `M HYD`.
- Define CX coordination and multicenter raise semantics.
- Decide coordinate ownership on the molecule wrapper and raise behavior for coordinate-bearing records.
- Re-export or expose `RaiseError` variants if stricter raise checks become public API.
- Close the existing `LitSet` TODO in raise code where appropriate.
- Add or update resolution fixtures for bond-only aromatic constraints.

### T6 — Migrate MOL/SDF parser implementation from nom to winnow

Doc 112 remains the active migration record for CTFile parsing internals.

Required work:

- Establish baseline MOL/SDF parser tests and benchmarks before changing parser internals.
- Inventory remaining nom-specific parser code.
- Build shared winnow parser infrastructure for fixed-width CTFile records.
- Replace the nom-specific `ParseError` layer with winnow-compatible errors.
- Preserve current behavior except where error classification, location reporting, or parser commitment is intentionally improved.
- Rework SDF data-block dispatch with deterministic line-oriented parsing where that is cleaner than directly translating the nom parser.
- Keep intermediate states controlled; the final migration boundary should leave MOL/SDF parsing green.

### T7 — Expand conformance, fuzzing, and benchmarks

Parsing coverage needs to track the boundary split and TableIR representation changes.

Required work:

- Add representative whole-record MOL and SDF benchmarks, not only primitive atom/bond line microbenchmarks.
- Add CXSMILES conformance categories and snapshots once CX has an explicit boundary.
- Add no-panic fuzz targets for CXSMILES parsing and for MOL/SDF parsing.
- Add parse-plus-raise fuzz or property coverage where the AST semantics are well-defined.
- Move corpus-validity checks out of parser unit tests and into dedicated
  integration tests so unit tests do not repeatedly validate fixture setup.
- Audit conformance-test feature gates. Expensive or externally sourced
  conformance suites must be gated consistently, while ordinary semantic
  regressions remain in the default test surface.
- Refresh classification tools after the ordinary-SMILES/CXSMILES split.
- Preserve regression fixtures for known parser and raise issue classes.

### T8 — Add Python-facing format APIs after Rust boundaries stabilize

The current Python workflow round focuses on resolved SMILES. MOL, SDF, and CXSMILES should not be forced through the same API before the Rust boundary types are settled.

Required work:

- Bind ordinary SMILES config and errors for the current Python workflow.
- Interpret reaction SMILES as `ReactionAst` through the Rust format and graph
  boundaries as tracked in
  [doc 170](170-reaction-smiles-python-2026-07-28.md).
- Defer Python `CxSmiles`, `Mol`, and `Sdf` APIs until their Rust boundary objects and configs exist.
- Keep parsing methods operation-specific, with separate config types for SMILES, MOL, SDF, and CXSMILES.
- Avoid generic format-polymorphic parsing APIs unless a later design shows concrete value.

## Immediate decision points

These are the main ordering and design choices that remain open:

- Whether `CxSmiles` initially wraps current `ExtendedMolecule` as an interim measure or waits for the compact semantic superset.
- Whether `Mol`/`Sdf` wrappers should land before or after the CTFile winnow migration.
- Whether `ChiralityFrame` relocation should be bundled with the TableIR semantic-superset work or handled earlier as a targeted cleanup.
- What benchmark and conformance evidence is required before old direct parser helpers are removed or redirected.

## Non-goals

- Do not reintroduce a `basic_opensmiles` mode.
- Do not reintroduce a `WILDCARDS` parse flag.
- Do not move `EXTENDED_AROMATICS` or `EXTENDED_BONDS` to CXSMILES solely because they previously lived beside extended-parser code.
- Do not propagate CXSMILES, MOL, or SDF Python APIs before the Rust boundary types are explicit.
