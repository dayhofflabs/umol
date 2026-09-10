# SMILES / MOL / SDF parsing outstanding tasks

Status: Proposed
Date: 2026-07-18
Relates: [047](047-smiles-conformance-suite-2026-01-21.md), [048](048-smiles-parser-configuration-2026-01-23.md), [100](100-table-ir-raise-ast-2026-05-27.md), [112](112-winnow-unification-2026-06-13.md), [151](151-python-molecule-workflows-2026-07-13.md), [152](152-basic-molecule-wildcards-2026-07-18.md),
[206](206-umol-perm-review-2026-08-21.md),
[217](217-rhea-participant-failures-2026-08-30.md),
[224](224-smiles-ring-closure-frame-2026-09-08.md),
[225](225-depiction-problems-2026-09-09.md)

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

### T4 — Extend stereo format coverage

[224](224-smiles-ring-closure-frame-2026-09-08.md) implements the targeted frame replacement:
SMILES finalization publishes complete stereo-atom frames, raising maps them directly to graph IR,
and ChiralityFrame is removed. MOL parity remains unread; wedge and coordinate interpretation is
unchanged. Limited atom-frame and directional-bond permutations are covered through parsing,
raising, resolution, and supported Python ingestion.

Remaining work:

- Define operative frames for additional supported stereo kinds and CX enhanced stereo when their
  format semantics are settled.
- Extend MOL/SDF stereo coverage as those boundaries develop; raw parity does not currently define
  an operative frame.
- Coordinate future format coverage with TableIR-to-graph-IR raising.

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

### T6 — Unify parser implementations on Winnow

Doc 112 owns direct nom removal from the workspace at the `0.8.0` compatibility
boundary.

Required work:

- Update and verify the Winnow 1.x dependency used by `umol-edn`.
- Remove the unused nom dependency from `umol-geometric`.
- Port the private CXSMILES parser without folding in the future CXSMILES boundary
  redesign.
- Establish baseline MOL/SDF parser tests and benchmarks before changing parser internals.
- Build shared winnow parser infrastructure for fixed-width CTFile records.
- Replace the nom-specific CTfile `ParseError` layer with parser-library-neutral errors.
- Preserve current behavior except where error classification, location reporting, or parser commitment is intentionally improved.
- Rework SDF data-block dispatch with deterministic line-oriented parsing where that is cleaner than directly translating the nom parser.
- Remove all direct nom dependencies and source references from the workspace.

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

The doc 224 catalog audit also found that the existing cis/trans-decalin resolution inputs each
resolve to the opposite catalog SMILES entry. Visual review confirmed that their names were reversed.
On 2026-09-09 the fixture names were swapped together with their matching snapshots, preserving
the inline constraint inputs and resolution settings. General randomized permutation testing has
no selected scope in doc 224.

#### Stereo fixture verification — 2026-09-09

Depiction defects and geometry ownership are tracked separately in
[doc 225](225-depiction-problems-2026-09-09.md); fixture identity review remains here.

The table tracks identity review of current resolution fixtures, separately from correction of an
identified defect. Reviewed outcomes record the user's visual review; pending rows await review.
For all rows, agreement with a catalog SMILES alone is not independent verification of its chemical name.

Depictions use the current EDN inputs and configuration overrides, resolved with the conformance
suite's counts model and MostSaturated tie-break, then rendered by Depict::depict and
Depiction::render_svg. Sheet composition changes only placement, scale, labels, and SVG identifiers.
The fixture filenames label the drawings; they do not assert that the names are correct.

| Fixture | Verification | Remaining check or action |
| --- | --- | --- |
| cis-decalin | Reviewed — name corrected | Formerly named trans-decalin; fixture and matching snapshot renamed together. |
| trans-decalin | Reviewed — name corrected | Formerly named cis-decalin; fixture and matching snapshot renamed together. |
| alpha-d-glucopyranose | Reviewed — appears correct | None identified in visual review. |
| cis-1-2-dichlorocyclohexane | Reviewed — appears correct | None identified in visual review. |
| trans-1-2-dichlorocyclohexane | Reviewed — appears correct | None identified in visual review. |
| r-methyloxirane | Reviewed — correct | Confirmed in user visual review. |
| l-ascorbic-acid | Reviewed — correct | Confirmed in user visual review. |
| 2r3r-dichlorobutane | Reviewed — correct | Confirmed in user visual review. |
| 2s3s-dichlorobutane | Reviewed — correct | Confirmed in user visual review. |
| meso-dichlorobutane | Reviewed — correct | Confirmed in user visual review. |
| 2-3-4-trichloropentane | Reviewed — correct structure | User identifies (2S,4S)-2,3,4-trichloropentane; C3 is not stereogenic. The current resolved output nevertheless retains a stereo-atom entity at C3 (site 3); investigate that discrepancy. |
| 2r3e-pent-3-en-2-ol | Reviewed — correct | Confirmed in user visual review. |
| e-cyclooctene | Pending — depiction unavailable | Resolution succeeds; CoordGen ignores ring E/Z below nine atoms and returns geometry inconsistent with the requested OppositeSide relation at bond 2. |
| z-cyclooctene | Pending | Named E/Z assignment. |
| e-2-3-difluorobut-2-ene | Pending | Named E/Z assignment. |
| z-2-3-difluorobut-2-ene | Pending | Named E/Z assignment. |
| z-2-fluorobut-2-ene | Pending | Named E/Z assignment. |
| e-azomethane | Pending | Named E/Z assignment. |
| z-azomethane | Pending | Named E/Z assignment. |
| z-butan-2-one-oxime | Pending | Named E/Z assignment. |
| 2e4e-hexa-2-4-diene | Pending | Named E/Z assignments. |
| 2e-hexa-2-4-diene-partial | Pending | First double bond E; second unspecified. |
| cyclohexene-asserted | Pending | Verify why the supplied #C1 assertion leaves no stereo entity in the successful resolved output; depiction alone cannot check this policy. |

The E-cyclooctene depiction failure is localized to the vendored CoordGen backend. Its MACROCYCLE
threshold is nine atoms; sketcherMinimizerBond::isStereo excludes bonds in smaller rings, and
setAbsoluteStereoFromStereoInfo uses that predicate before applying the supplied cis/trans relation.
The fixture resolves with an E stereo bond and the IO adapter supplies OppositeSide. Returned
coordinates fail the Rust boundary's independent relative-side check, which raises
CisTransGeometryMismatch instead of publishing the wrong drawing. This is a layout limitation,
not evidence that the fixture resolved as Z. A bounded trial changing MACROCYCLE from nine to eight
makes E-cyclooctene depict successfully but makes Z-cyclooctene fail the SameSide geometry check.
The other 21 previously rendered fixtures produce identical SVGs. The macro also selects the ring
layout algorithm, so this is not an isolated stereo-recognition switch. The trial was rejected and
the nine-atom definition restored; no fixture or geometry validator was changed permanently.
A second trial bypassed the output geometry check at the eight-atom threshold to inspect both
returned drawings. E has well-separated opposite-side substituents. Z has one ring bond nearly
collinear with the double bond (about 0.005 degrees from its line, normalized cross product magnitude
0.0000866). Its small deviation is on the opposite side and exceeds the validator's 0.000001
relative tolerance. Thus the reported Z mismatch comes from nearly degenerate generated geometry,
not a cleanly drawn E arrangement. Both temporary changes were restored after capture.

These files are in the stereo_tetrahedral and stereo_cis_trans resolution fixture directories.
The 35-entry catalog comparison covered all of them except cyclohexene-asserted; after the three
corrections in doc 224, only the two decalins differed from their catalog references. This selected
review list does not claim independent identity verification of every stereo fixture in the repo.

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

### T9 — Verify OpenSMILES arrangement numbering against umol-perm cosets

Moved from the umol-perm review ([206](206-umol-perm-review-2026-08-21.md), open item 5).

- Pin the TH and AL index-to-arrangement-number correspondence with fixtures; check first
  whether the existing SMILES conformance suite already exercises it end to end.
- Verify the TB/OH enantiomer pairing against the OpenSMILES `@`/`@@` numbering (the unverified
  note in `ClassKey::build`); the specification copy is in `materials/formats/opensmiles`, and
  RDKit can serve as cross-validation.
- Pin the fixed-point-freeness of the TB/OH axial swaps once the pairing is verified; the
  review's exact scan (moved cosets TH 2/2, AX 2/2, TB 20/20, OH 30/30; CT 0/2, SP 0/3) is the
  expected evidence.

## Immediate decision points

These are the main ordering and design choices that remain open:

- Whether `CxSmiles` initially wraps current `ExtendedMolecule` as an interim measure or waits for the compact semantic superset.
- Whether `Mol`/`Sdf` wrappers should land before or after the CTFile winnow migration.
- What benchmark and conformance evidence is required before old direct parser helpers are removed or redirected.

## Non-goals

- Do not reintroduce a `basic_opensmiles` mode.
- Do not reintroduce a `WILDCARDS` parse flag.
- Do not move `EXTENDED_AROMATICS` or `EXTENDED_BONDS` to CXSMILES solely because they previously lived beside extended-parser code.
- Do not propagate CXSMILES, MOL, or SDF Python APIs before the Rust boundary types are explicit.
