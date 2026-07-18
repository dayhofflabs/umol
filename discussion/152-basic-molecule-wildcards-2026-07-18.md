# OpenSMILES wildcards in basic TableIR

Status: **Active implementation plan**

Date: 2026-07-18

Relates: 047 (SMILES conformance suite), 048 (SMILES parser configuration),
100 (TableIR raise), 151 (Python molecule workflows and the TableIR
representation decision)

## Goal

The ordinary `Smiles` boundary must accept the complete OpenSMILES grammar,
including the `*` atom, without routing ordinary input through
`ExtendedMolecule`. The basic `table_ir::Molecule` representation remains the
compact result type. A wildcard atom raises to `ElementAst::Undetermined`; a
workflow that requires a determined molecule may subsequently report
underdetermination, but parsing itself succeeds.

The implementation includes unit tests, semantic property tests, fuzz coverage,
and paired basic/extended parsing benchmarks. It also removes the artificial
configuration and conformance-classification split between "basic OpenSMILES"
and OpenSMILES once the basic parser supports `*`.

CXSMILES, MOL/SDF boundary unification, and the eventual compact semantic
superset for all TableIR formats remain outside this task.

## Representation

The basic atom stores:

```rust
pub element: Option<Element>
```

`Some(element)` denotes a concrete element and `None` denotes the OpenSMILES
`*`. This avoids adding the much larger extended `AtomSymbol` representation to
every basic atom and maps directly to `ElementAst::Undetermined`.

A wildcard atom has `aromatic: None`: `*` does not constrain aromaticity. Other
atom fields remain independent, so bracketed wildcard forms may still carry an
isotope, charge, hydrogen count, chirality, class, span, or other field already
representable by basic `Atom`.

## Staged implementation plan

### S0 — Basic TableIR representation

- **S0a — compact wildcard atoms** **Done**
  (`umol-io/src/table_ir/atom.rs` and direct `Atom.element` consumers): change
  `Atom.element` from `Element` to `Option<Element>`. Existing element
  constructors continue to accept `Element` and store `Some(element)`; add
  wildcard constructors, including a span-preserving form. Convert
  `Atom { element: None, .. }` to
  `AtomSymbol::WildcardAtom(WildcardAtom::Any)` and accept that exact extended
  symbol in `TryFrom<ExtendedAtom> for Atom`, while continuing to reject other
  wildcard kinds and extended symbols. Migrate concrete-element consumers in
  the CTFile and SMILES implementations without unchecked unwraps. Record the
  resulting `Atom` layout and require that the optional element does not
  materially increase the current 104-byte representation.

  Rust table tests cover concrete and wildcard constructors, concrete and
  wildcard basic/extended roundtrips, and continued rejection of heavy,
  heteroatom, halogen, metal, atom-list, and other extended symbols.
  **Breaking representation migration (red to green).** `[dep: none]`

  **Implemented verification:** `Atom` remains 104 bytes. Constructor and
  conversion tests cover concrete elements, the `Any` wildcard, and rejection
  of every other wildcard category and extended atom-list symbol.

- **S0b — molecule and AST semantics** **Done**
  (`umol-io/src/table_ir/molecule.rs`, `table_ir/raise.rs`): include basic
  wildcards in `Molecule::sum_formula`, using the existing extended convention
  (`*`, `*2`, and so on). Raise `Some(element)` to `ElementAst::Lit(element)`
  and `None` to `ElementAst::Undetermined`. A wildcard does not receive the
  aromatic-heteroatom implicit-hydrogen rule. Isotope, charge, hydrogen, spin,
  class, stereo, and other independent fields remain preserved.

  Rust table tests cover basic formulas such as `C3O2*` and `C*2`, exact
  `AtomAst` values for bare and bracketed wildcards, and basic-to-extended-to-
  basic molecule roundtrips containing wildcards. **Additive behavioral
  completion (green).** `[dep: S0a]`

  **Implemented verification:** basic formulas count one or multiple
  wildcards; raising preserves independent isotope, charge, hydrogen,
  lone-pair, and spin fields, including independent multiplicity; aromatic
  wildcards retain undetermined implicit hydrogens. A molecule roundtrip also
  preserves wildcard class, stereo, source span, label, and value fields.

S0 ends with the `umol-io` unit suite green.

### S1 — Basic SMILES parsing

- **S1a — bare and bracketed wildcard grammar** **Done**
  (`umol-io/src/smiles/parser.rs`, `parser/utils.rs`, `parser/builder.rs`): add
  a basic-builder wildcard path that creates an atom with `element: None` and
  `aromatic: None`. Parse bare `*` through the ordinary attachment, branch,
  component, and ring machinery. Change the basic bracket parser to return an
  optional element so bracketed wildcard fields remain representable. Support
  the OpenSMILES wildcard shapes already accepted by the extended parser. `*`
  is a core grammar token rather than behavior controlled by an optional flag.

  Rust table cases assert exact TableIR atoms, bonds, fields, and spans for at
  least `*`, `C*C`, `C(*)C`, `C-*`, `*.*`, `*1CC1`, `[*]`, `[*:1]`, and
  representative bracket fields accepted by the extended grammar. **Additive
  parser support (green).** `[dep: S0a]`

  **Implemented verification:** nine exact TableIR rows cover bare, chain,
  branch, explicit-bond, disconnected-component, ring, bracket, class, and
  combined bracket-field forms. Wildcards carry `element: None` and
  `aromatic: None`; all atom and bond spans are asserted structurally.

- **S1b — format, reaction, and raise integration** **Done**
  (`umol-io/src/smiles.rs`, the SMILES reaction parser tests, and TableIR raise
  tests): verify that the ordinary `Smiles` wrapper accepts wildcard input,
  that reaction SMILES accepts wildcards on reactant, agent, and product sides,
  and that `Smiles -> Molecule -> MoleculeAst` produces undetermined elements.
  At the graph workflow boundary, pin the semantic distinction:
  `Smiles::parse("*")` succeeds, while ingestion as a determined molecule
  reports `Underdetermined` rather than a syntax or model-conversion error.

  Rust table tests cover each boundary and exact error/result variants.
  **Additive integration coverage (green).** `[dep: S0b, S1a]`

  **Implemented verification:** the ordinary `Smiles` wrapper returns exact
  wildcard TableIR; a single exact `*>*>*` reaction assertion covers reactant,
  agent, and product molecules, including side-local spans and format
  metadata; and the wrapper-to-TableIR-to-AST path yields an undetermined
  element with the expected IO ground defaults. Both valence planners classify
  a non-literal element as `Underdetermined` with an empty edit plan; direct
  valence resolution and the composite resolver return without mutation or
  running later stages. Consequently, `Smiles::ingest` and the text convenience
  API report their exact `Underdetermined` variants without an ingestion-layer
  precheck.

- **S1c — wildcard property tests**
  (`umol-io/tests/smiles_property.rs`): add `*` to the existing SMILES-biased
  arbitrary-input alphabet. Add generated valid `C`/`*` chains containing at
  least one wildcard and assert atom count, wildcard positions, bond endpoints,
  and source spans. Parse generated valid wildcard structures through both the
  basic and extended parsers and require the basic result, converted to
  `ExtendedMolecule`, to equal the extended result. Raise generated chains and
  require `C -> ElementAst::Lit(Element::C)` and
  `* -> ElementAst::Undetermined` at every position.

  The differential property is the primary migration oracle: it compares the
  new implementation with the existing independent extended implementation
  rather than only restating aggregate expectations. **Additive property
  coverage (green).** `[dep: S0b, S1a]`

S1 ends with the unit, integration, and property suites green.

### S2 — Remove the wildcard configuration and classification split

- **S2a — consolidate SMILES configuration**
  (`umol-io/src/smiles/config.rs`, parser defaults, and `Smiles` defaults):
  remove `SmilesParseFlags::WILDCARDS`, because `*` is part of OpenSMILES.
  Remove `BASIC_OPENSMILES` and `SmilesIoConfig::basic_opensmiles()`. Remove or
  consolidate presets whose only distinction was wildcard acceptance instead
  of retaining identical aliases. The ordinary default becomes
  `SmilesIoConfig::opensmiles()`. Extended-aromatic and extended-bond flags
  remain independent acceptance-policy options.

  Rust table tests cover the remaining named configurations, flag composition,
  and display behavior. **Breaking public configuration change (red until
  callers migrate).** `[dep: S1b]`

- **S2b — classification and parsing-conformance migration**
  (`umol-graph/src/bin/classify_smiles_strings.rs`, related diagnostic tools,
  `umol-io/tests/smiles_parsing`, and their fixtures/snapshots): migrate both
  consumers of the old split. In the `umol-graph` classification tools, remove
  the separate "basic OpenSMILES without wildcards" result, category, counts,
  reports, and ordering assumptions. In the `umol-io` parsing-conformance
  suite, remove the corresponding parse result and category, update shortcut
  logic that previously inferred OpenSMILES success from basic-OpenSMILES
  success, and update expected category metadata and snapshots to the single
  OpenSMILES classification. Migrate `test_smiles` and every other workspace
  caller of the retired configuration surface in the same subitem.

  Conformance tests continue to exercise the same source corpora, now with
  wildcard acceptance as part of the ordinary OpenSMILES result. Exact
  classification/report tests pin the revised category structure and prevent
  reintroduction of the result-type-driven parser distinction. **Breaking
  caller and conformance migration, restoring green.** `[dep: S2a]`

S2 ends with the workspace green and no wildcard capability switch or
basic-versus-OpenSMILES conformance category.

### S3 — Fuzz and performance gates

- **S3a — wildcard fuzz coverage**
  (`umol-io/fuzz`): add corpus seeds for a bare wildcard, a mixed wildcard
  chain, a branch with a bracketed/classed wildcard, a wildcard ring, and
  representative bracket fields. On a successful parse, extend the existing
  fuzz target to attempt TableIR-to-`MoleculeAst` raising inside the same panic
  boundary. Retain the general arbitrary-byte target rather than creating a
  wildcard-only target. Compile the fuzz crate, replay the seed corpus, and run
  a bounded fuzz session. **Additive fuzz coverage (green).**
  `[dep: S1b, S2b]`

- **S3b — paired basic/extended benchmarks**
  (`umol-io/benches/smiles_parsing.rs`): activate a representative wildcard
  corpus containing a single bare wildcard, a long bare-wildcard chain, a
  bracketed wildcard chain, a mixed element/wildcard chain, and wildcard
  branches or rings. Benchmark every input through both ordinary
  `Smiles::parse_bytes` and `parse_extended_smiles_bytes`. Retain the existing
  non-wildcard basic benchmark groups as regression controls.

  Record the basic `Atom` layout, paired wildcard latency, and representative
  ordinary chain, ring, and bracket latency. The basic atom remains compact;
  ordinary non-wildcard parsing must not regress beyond measurement noise; and
  wildcard parsing should follow the basic-path cost profile rather than the
  current extended-path profile. **Additive benchmark gate (green).**
  `[dep: S1a, S2b]`

S3 ends with the fuzz target compiling and replaying its seeds, the bounded fuzz
run clean, and the benchmark comparison recorded.

## Dependencies and completion

The critical path is:

```text
S0a -> S0b -> S1a -> S1b -> S2a -> S2b -> S3
```

S1c may proceed after S1a and S0b. S3a and S3b may proceed independently once
their listed dependencies are complete, but both are part of the completion
gate. No stage in this plan is deferrable. CXSMILES/MOL/SDF unification and the
larger TableIR representation spike remain explicitly outside the task.
