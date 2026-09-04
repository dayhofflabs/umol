# 112 — Winnow parser unification

Status: Completed
Date: 2026-06-13
Relates: [153](153-format-parsing-outstanding-tasks-2026-07-18.md)

## Purpose

Replace direct nom use throughout the workspace with Winnow and remove nom at the
`0.8.0` compatibility boundary. The migration should leave each parser idiomatic for
its input format rather than introduce a common parser framework.

This supersedes the original CTfile-only plan. That plan incorrectly described CTfile
as the last production nom user and combined a useful technical design with an obsolete
implementation sequence.

## Current evidence

As of 2026-09-03:

- `umol-edn` and `umol-graph-ir` depend directly on Winnow 1.x. The lockfile resolves
  Winnow 1.0.4.
- Winnow 1.0.4 is current. Its changes since 1.0.1 concern `seq!` field capacity,
  `ascii::float`, `binary::bits`, and documentation; none applies to the API used by
  `umol-edn`.
- The CTfile parser and the private CXSMILES parser in `umol-io` use nom 8;
  `umol-io` also depends on Winnow 1.x for the staged migration.
- The unused nom 8 dependency has been removed from `umol-geometric`; `umol-io` is
  now the only direct nom user in the workspace.
- `ctfile::ParseError` publicly exposes `nom::error::ErrorKind`, and its public
  conversion functions accept nom error types. Removing nom therefore changes public
  API even if accepted CTfile input remains unchanged.
- `ctfile::parser` also exposes CTAB block parsers and line parsers whose return types
  name nom's parser API. The line parsers are exported only for Criterion benchmarks,
  and no parser-returning function has a workspace consumer outside `umol-io`.
- The current MOL benchmark measures individual counts, atom, bond, atom-list, and
  property records. There is no full-file SDF benchmark, and the extended-SMILES
  benchmark contains no actual CX annotation.

## Compatibility boundary

Complete removal of nom is a `0.8.0` change. The migration will not retain deprecated
nom-bearing variants, conversion functions, optional compatibility features, or other
temporary shims for a `0.7.x` release.

The breaking surface comprises the CTfile error API and removal of the public
parser-returning functions. Parser-library replacement does not otherwise authorize
changes to public format boundaries, accepted syntax, or constructed TableIR values.

The existing whole-MOL and whole-SDF entry points remain public. Make `ctab_block`,
`extended_ctab_block`, `counts_input`, the basic and extended atom and bond inputs,
`legacy_atom_list_input`, and the basic and extended property inputs private. Do not
replace them with a new public benchmark-support API.

## Scope

- Update the locked Winnow 1.x release and verify `umol-edn` without redesigning its
  parser.
- Remove the unused nom dependency from `umol-geometric`.
- Port the private CXSMILES parser in `umol-io` from nom to Winnow while preserving its
  current public entry points and errors.
- Port MOL/SDF CTfile parsing from nom to Winnow.
- Replace the nom-bearing CTfile error surface with parser-library-neutral errors.
- Retire the record-level CTfile microbenchmarks and benchmark a small set of public
  full-file parsing operations instead.
- Make the parser-returning CTfile implementation functions private.
- Remove nom from workspace manifests and the lockfile.
- Preserve conformance behavior and check parsing performance against a baseline
  recorded before implementation.

## Boundaries

This work does not:

- redesign the ordinary SMILES/CXSMILES boundary;
- introduce `Mol` or `Sdf` boundary objects;
- replace the current basic/extended TableIR split;
- change TableIR-to-graph-IR ingestion semantics;
- centralize the independent EDN, graph-IR, SMILES, and CTfile error types; or
- require those parsers to share helpers merely because they use the same library.

Those format and representation questions remain in doc 153.

## Winnow version update

The direct dependency declarations should remain `winnow = "1.0"`; the manifest need
not pin a patch release. Update the lockfile to the current compatible release and run
the complete `umol-edn` and `umol-graph-ir` feature surfaces. Source changes are
required only if compilation or tests identify an incompatibility.

The older Winnow 0.7 lockfile entry is transitive through test infrastructure. It is not
part of the direct parser migration and need not be eliminated by adding dependency
overrides.

## CXSMILES migration

The CXSMILES parser is private parser machinery behind the existing SMILES APIs. Port
its combinators and numeric parsing directly to Winnow, preserving:

- the existing CX fields and returned payloads;
- accepted and rejected inputs;
- current public `smiles::ParseError` values; and
- current configuration behavior.

Do not combine this port with the future separation of ordinary SMILES and CXSMILES.

## CTfile parser design

### Input and parser shape

Use Winnow's mutable-input model and `.parse_next()` internally. Stateful block parsers
should remain direct loops: atom, bond, property, and SDF blocks carry counts, physical
line numbers, feature flags, and cross-line state that are clearer in imperative code
than in nested combinators.

Counts, atoms, bonds, and most property records are fixed-column physical lines. Parse a
line boundary once, then use a line-local `LocatingSlice<&[u8]>`. The local offset is the
column, while the block parser supplies the physical line number. Fixed-width helpers
must preserve the field-start column when parsing a detached field.

Direct indexing remains appropriate in optimized fixed-column parsers. Migration to
Winnow does not require expressing every field as a combinator expression.

### Dispatch and commitment

Use deterministic prefix dispatch for record kinds with unique tags, especially CTfile
property records. Use `alt` only where more than one grammar branch can own the same
input.

MOL is a fixed-column, block-structured format. Once block context selects a counts,
atom, bond, or other fixed record, a malformed field is a committed failure and normally
does not require backtracking. The main exceptions are block-boundary recognition and
the multiline property records whose ownership can depend on a following physical line.

`Backtrack` therefore means that the parser has not yet determined whether the input
belongs to the current block or multiline construct. At clear failure sites the port may
replace a recoverable error with `Cut`, including semantic construction and consistency
failures after a construct has been identified. This is permitted cleanup rather than a
requirement to convert every existing recoverable error.

For SDF data blocks, prefer record-oriented line dispatch between a data field beginning
with `>` and the required `$$$$` delimiter. A malformed recognized field or delimiter
must not be reinterpreted as the other record kind.

### CTfile error contract

The replacement error surface must not expose Winnow or nom types. It should retain
domain errors that identify format concepts, including unexpected end of a counted
block, a missing `M  END`, a missing SDF delimiter, invalid property accumulation, and
inconsistent SGroups.

Remove:

- `Incomplete`, because CTfile entry points parse complete inputs;
- `NomError` and `NomErrorKind`;
- the nom `ParseError` implementation and `From<nom::Err<_>>`; and
- `from_nom` and the block-specific `*_from_nom` functions.

Syntax failures should identify the actual failed format condition and its line and
column where available. Derive the flat public variants from the current parser,
conversion, and accumulator failure sites. Preserve reachable semantic errors, remove
unused variants, and do not introduce nested record-kind or syntax-kind taxonomies merely
to translate parser-library errors.

Do not add a generic external-error conversion that erases domain conversion failures.
Use explicit conversion errors and apply `Backtrack` or `Cut` at the grammar site where
commitment becomes known.

## Public contract

`ctfile::ParseError` remains the public failure type of the existing MOL and SDF entry
points. It reports the format or construction condition that prevented the requested
TableIR value and does not expose nom or Winnow values. Complete-input parsing has no
`Incomplete` outcome, and parsing does not panic on malformed input.

The existing whole-file parsing functions retain their names, inputs, configuration,
output types, and accepted syntax. No CTfile API is added to Python by this work. The
parser-returning CTAB and record functions are removed from the public surface rather
than translated to Winnow return types.

## Behavioral and performance evidence

Before the corresponding parser migration, record:

- the complete `umol-edn` test result against its current locked Winnow version;
- `umol-io` unit, integration, and feature-gated conformance results;
- representative CXSMILES benchmark results before the CXSMILES port;
- representative full-file MOL and SDF benchmark results against the nom implementation
  before the CTfile port; and
- the complete inventory of direct nom dependencies and source references.

During migration, preserve accepted input and constructed values unless a separately
recorded correction is required. Add focused evidence for commitment and location
behavior, including malformed fixed-width fields, recognized property bodies, truncated
counted blocks, malformed SDF headers and delimiters, zero-field SDF records, and LF/CRLF
locations.

The completed work must leave no direct nom dependency or nom source reference in the
workspace and no unexplained material full-file parsing regression.

## Implementation plan

### S0 — Dependency and evidence baseline

- **S0a — Current parser verification.** **Done.** (`umol-edn`, `umol-graph-ir`,
  `umol-io`): record the existing test, integration, and applicable feature-gated
  conformance results before changing dependencies. Additive, green. [dep: none]

  Baseline on 2026-09-03:

  - `cargo test -p umol-edn --all-features`: 1,218 passed;
  - `cargo test -p umol-graph-ir --all-features`: 6,933 passed, 7 ignored;
  - `cargo test -p umol-io`: 3,393 passed;
  - `cargo test -p umol-io --features conformance`: 16,082 passed; and
  - `cargo test -p umol-io --features proptest`: 3,401 passed.

  All commands completed with no failures.
- **S0b — Winnow patch update.** **Done.** (`Cargo.lock`, `umol-edn`,
  `umol-graph-ir`): retain `winnow = "1.0"`, update the locked compatible patch release,
  and run both crates' complete feature surfaces. Additive, green. [dep: S0a]

  Winnow 1.x was updated from 1.0.1 to 1.0.4 without changing either manifest.
  `cargo test -p umol-edn --all-features` passed 1,218 tests, and
  `cargo test -p umol-graph-ir --all-features` passed 6,933 tests with 7 ignored.
- **S0c — Unused dependency removal.** **Done.** (`umol-geometric/Cargo.toml`): remove
  nom and verify `umol-geometric`. Additive, green. [dep: S0a]

  `cargo test -p umol-geometric --all-features` passed 218 tests. `cargo tree -i nom`
  now reports only `umol-io` as a direct workspace dependency.
- **S0d — `umol-io` migration dependency and CX benchmark.** **Done.**
  (`umol-io/Cargo.toml`, `umol-io/benches/smiles_parsing.rs`): add Winnow alongside nom
  and add a small set of genuine CX-annotation cases to the existing top-level
  extended-SMILES benchmark. Record the nom baseline. Additive, green. [dep: S0b]

  The nom baseline on 2026-09-03 was:

  - coordinates: 263.77–267.76 ns;
  - ferrocene multicenter bonds: 750.11–753.19 ns; and
  - SGroup hierarchy: 302.30–305.83 ns.

  `cargo test -p umol-io` passed 3,393 tests after the dependency and benchmark
  additions.

The workspace is green at the end of S0.

### S1 — Private CXSMILES parser

- **S1a — CXSMILES parser port.** **Done.**
  (`umol-io/src/smiles/parser/cx.rs`): port the private CX parsing functions and their
  tests to Winnow, preserving accepted input, returned payloads, configuration behavior,
  unknown-tag handling, and public `smiles::ParseError` results. Additive with respect to
  public API, green. [dep: S0d]

  The private parser now uses Winnow's mutable byte-slice input and its backtrack/cut
  distinction; no CX parser details or new error types were added to the public API.
  `cargo test -p umol-io cx` passed 196 focused tests, and `cargo test -p umol-io`
  passed 3,393 tests in total (3,387 unit tests and 6 integration tests).
  `cargo clippy -p umol-io --all-targets -- -D warnings` also passed.
- **S1b — CXSMILES evidence.** **Done.** (`umol-io` tests and `smiles_parsing`
  benchmark): run the existing CX tests and the S0d top-level CX cases, and account for
  any material regression. Additive, green. [dep: S1a]

  The Winnow measurements on 2026-09-03 were:

  - coordinates: 256.34–262.33 ns;
  - ferrocene multicenter bonds: 739.49–746.18 ns; and
  - SGroup hierarchy: 296.27–298.10 ns.

  All three intervals are slightly lower than the recorded nom baseline, so there is no
  material regression. `cargo test -p umol-io cx` passed 196 focused tests.

The workspace is green at the end of S1; CTfile still uses nom.

### S2 — CTfile parser and error migration

- **S2a — Full-file CTfile benchmark baseline and privacy cleanup.** **Done.**
  (`umol-io/benches/mol_parsing.rs`, `umol-io/src/ctfile/parser.rs` and record modules):
  replace record-level microbenchmarks with a few representative top-level basic MOL,
  extended MOL, and SDF parsing cases; make the CTAB block and record parser functions
  private; and record the full-file baseline against the unchanged nom implementation.
  Breaking public-surface cleanup, red to green within the subitem. [dep: S0a]

  The benchmark now measures the public whole-file parsers using the existing caffeine
  MOL, RDKit copolymer/SGroup MOL, and ten-component wwPDB SDF conformance inputs. The
  unchanged nom implementation measured 2.7175–2.7352 us, 4.6568–4.6641 us, and
  78.328–78.446 us, respectively, on 2026-09-03. The benchmark-only public exposure of
  the CTAB block and record parsers was removed; the public whole-file parsers are
  unchanged. `cargo test -p umol-io` passed 3,393 tests (3,387 unit tests and 6
  integration tests), and `cargo clippy -p umol-io --all-targets -- -D warnings`
  passed.
- **S2b — Parser-neutral errors and input foundation.** **Done.**
  (`umol-io/src/ctfile/error.rs`, `umol-io/src/ctfile/parser/utils.rs`): replace the
  nom-bearing public error machinery with the cause-based contract above and port line,
  location, fixed-width, integer, float, and field helpers to Winnow. Breaking
  coordinated migration. [dep: S2a]

  `ctfile::ParseError` no longer exposes nom types, streaming-incomplete states, or
  generic parser errors. Its remaining variants describe reachable CTfile syntax and
  construction conditions; counted-record and SDF-header errors carry physical
  locations, and invalid references distinguish atom and bond indices. The private
  parser foundation now uses line-local Winnow `LocatingSlice` inputs and preserves
  field-start columns while retaining CTfile's zero-padded decimal and LF/CRLF
  behavior. An isolated build of the foundation passed 182 focused tests.

  As specified for the coordinated S2b–S2f rewire, `cargo check -p umol-io --lib`
  remains red at the unmigrated record and composition modules: they still import the
  removed line iterator, call the former helper signatures and nom conversions, and
  require nom's error traits. No compatibility adapters were added; S2c owns the first
  consumer migration.
- **S2c — Fixed CTfile records** (`header.rs`, `counts.rs`, `atom.rs`, `bond.rs`,
  `legacy_atom_list.rs`, `rgroup.rs`): port each record family and its exact success and
  committed failure cases, preserving physical line and column reporting. Breaking
  coordinated migration. [dep: S2b]

  The fixed header, counts, atom, bond, legacy atom-list, and R-group parsers now use
  Winnow directly. Counted blocks consume physical lines through the shared line reader;
  missing records are committed `UnexpectedEof` failures, while malformed owned records
  are committed record-specific failures at the field-start or first trailing byte
  column. Enumerated atom and bond fields are checked before conversion, so malformed
  codes report syntax errors rather than reaching infallible conversion branches. The
  existing success partitions were retained and the failure tables now assert exact
  Winnow columns; explicit counted-block EOF and formerly panic-prone code cases were
  added. An isolated fixed-record harness passed 802 tests on Rust 1.87, and the same
  harness passed Clippy with warnings denied on the current toolchain.

  As expected for the coordinated S2b-S2f migration, the crate remains red at the
  unmigrated property, accumulation, and top-level composition modules. No Nom adapter
  was added; S2d owns the next consumers.
- **S2d — SGroups and properties** (`sgroup.rs`, `properties.rs`). **Done.** Port property
  dispatch and bodies, use `Cut` after recognized tags where ownership is clear, and
  retain the necessary backtracking at multiline-property boundaries. Breaking
  coordinated migration. [dep: S2b, S2c]

  SGroup helpers and basic and extended property parsers now use Winnow directly.
  Property dispatch remains backtrackable until a supported tag owns the line; failures
  in recognized property bodies are committed and the block reports the exact physical
  line and byte column. Atom aliases and legacy group abbreviations retain their
  two-line structure, including committed unexpected-EOF errors when the continuation
  line is absent. The original property test matrix was retained and exact Winnow error
  modes and locations replaced nom error kinds. An isolated harness passed 656 tests on
  Rust 1.87, including 297 property-parser and 56 SGroup tests, and passed Clippy with
  warnings denied on the current toolchain.

  As expected for the coordinated S2b-S2f migration, the crate remains red only at the
  unmigrated accumulation and top-level composition boundaries. No Nom adapter was
  added; S2e owns the next consumers.
- **S2e — Property accumulation and conversion** (`accumulator.rs`, `context.rs`,
  `convert.rs`). **Done.** Route duplicate properties, invalid codes and isotopes, invalid
  references, and SGroup consistency failures through the parser-neutral error contract.
  Breaking coordinated migration. [dep: S2d]

  Property accumulation now uses the settled atom- and bond-specific out-of-bounds
  errors and the consistently named unterminated SGroup-data error. Conversion already
  returned the parser-neutral invalid-code and invalid-isotope variants and required no
  implementation change; its tests now assert their complete fields. Accumulator tests
  likewise assert complete duplicate-property, reference, isotope, and SGroup errors,
  including type constraints, missing data content and termination, and mismatched data
  indices. `Context` remains the private state carrier and required no semantic change.
  An isolated harness passed 910 tests on Rust 1.87, including 133 accumulator and 121
  conversion tests, and passed Clippy with warnings denied on the current toolchain.

  As expected for the coordinated S2b-S2f migration, the crate remains red only in the
  unmigrated top-level composition. No Nom adapter was added; S2f owns the remaining
  CTAB, MOL, and SDF consumers.
- **S2f — CTAB, MOL, and SDF composition** (`parser.rs`, `sdf_data.rs`). **Done.** Port the
  top-level composition, retain backtracking where block ownership is not yet known, use
  deterministic record-oriented SDF dispatch, preserve all public whole-file entry
  points, and migrate their tests. This completes the breaking rewire and restores a
  green tree. [dep: S2c, S2d, S2e]

  CTAB composition now drives the migrated counted blocks through Winnow's mutable
  input directly, and the existing whole-MOL and whole-SDF functions retain their
  public signatures. SDF data parsing dispatches deterministically between `>` fields
  and the required `$$$$` delimiter. Recognized malformed headers report
  `InvalidSdfDataHeader`; absent or malformed delimiters report `MissingDelimiter`.
  Zero-field records, LF and CRLF input, inter-record whitespace, and the legacy
  repeated-terminator case have focused coverage. No Nom adapter or new public parser
  surface was added.

  The 313 MOL and 43 SDF snapshots affected by the parser-neutral error migration were
  audited and updated. Constructed values and success categories are unchanged; the
  differences are field columns, corrected cross-record line accounting, and the
  settled delimiter classification. `cargo test -p umol-io` passed 3,416 tests,
  `cargo test -p umol-io --features conformance` passed 16,105 tests, and
  `cargo test -p umol-io --features proptest` passed 3,424 tests. Clippy passed for all
  `umol-io` targets with warnings denied.
- **S2g — Full-file regression check.** **Done.** (`mol_parsing` benchmark): rerun the unchanged
  S2a benchmark set against Winnow and account for any material regression. Additive,
  green. [dep: S2f]

  The initial integer-field port converted matched bytes through UTF-8 and `FromStr` and
  caused a material regression. Integer fields now use `atoi` to parse their fixed-width
  ASCII bytes directly. The final measurements on 2026-09-03 were:

  - caffeine MOL: 2.9984–3.4649 us;
  - copolymer/SGroup MOL: 4.0682–4.1991 us; and
  - ten-component SDF: 74.959–76.963 us.

  The caffeine case remains a fraction of a microsecond slower than its nom baseline and
  was accepted as immaterial. The extended MOL and SDF cases are faster than their nom
  baselines. `cargo test -p umol-io` passed 3,416 tests, and Clippy passed for all
  `umol-io` targets with warnings denied.

S2b through S2f are one coordinated breaking migration: the shared parser result and
error types cross the CTfile modules, so the tree is required to be green at the stage
boundary rather than through artificial dual-parser adapters. The workspace is green at
the end of S2.

### S3 — Nom removal and final verification

- **S3a — Direct nom removal.** **Done.** (`umol-io/Cargo.toml`, workspace sources,
  `Cargo.lock`): remove the remaining dependency and imports, then confirm that no direct
  nom dependency or source reference remains. The transitive Winnow 0.7 test-infrastructure
  dependency is not part of this cleanup. Breaking cleanup, green. [dep: S1b, S2g]

  The final direct nom dependency was removed from `umol-io` and the lockfile. No nom
  package remains in the dependency graph, and no workspace source reference remains.
- **S3b — Verification gate.** **Done.** (workspace): run formatting; the default, integration,
  conformance, and property-test surfaces affected by the parser changes; workspace
  tests and lint; and the final full-file benchmarks. Additive, green. [dep: S3a]

  Nightly formatting passed. The final focused results were 1,218 `umol-edn` tests;
  6,933 `umol-graph-ir` tests with 7 ignored; and 3,416 default, 16,105 conformance,
  and 3,424 property-enabled `umol-io` tests. `cargo test --workspace` and
  `cargo clippy --workspace --all-targets -- -D warnings` passed under the repository
  Python 3.13 environment. The final full-file benchmark results are recorded in S2g.

The critical path is S0b → S0d → S1 → S3 and S0a → S2 → S3. No stage is deferrable:
the CXSMILES port, CTfile port, dependency removal, and regression evidence are all
required to complete workspace-wide nom removal.
