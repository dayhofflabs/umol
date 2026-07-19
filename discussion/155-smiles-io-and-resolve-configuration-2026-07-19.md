# SMILES IO and molecular resolution configuration

Status: **Active implementation plan**
Date: 2026-07-19
Relates: [094](094-dsl-ast-io-ergonomics-2026-05-07.md),
[100](100-table-ir-raise-ast-2026-05-27.md),
[151](151-python-molecule-workflows-2026-07-13.md),
[153](153-format-parsing-outstanding-tasks-2026-07-18.md)

## Purpose

The resolved SMILES workflow crosses three distinct policy boundaries:

1. the external SMILES representation accepted and emitted by `umol-io`;
2. the chemistry model used to interpret and resolve the raised `MoleculeAst`;
3. operational choices governing how resolution modifies its input.

These policies should remain distinct even when a higher-level Rust or Python
operation composes them. This document defines their responsibilities and the
round-trip contract required of the SMILES boundary before the Python workflow
surface is completed.

The design is not constrained by the currently implemented parser or the
absence of a renderer. The purpose is to establish the complete shape that the
eventual two-way format API must satisfy.

## 1. `SmilesIoConfig` is bidirectional

`SmilesIoConfig` configures the paired parse and render operations for ordinary
SMILES. It is not merely a parser-options container. Using one configuration in
both directions follows the precedent of the contextual `FromAst` / `IntoAst`
conversions: the configuration defines a compatible pair of transformations
rather than two unrelated collections of defaults.

A useful conceptual decomposition is:

```rust
pub struct SmilesIoConfig {
    pub syntax: SmilesSyntaxFlags,
    pub parse: SmilesParsePolicy,
    pub render: SmilesRenderPolicy,
}
```

The exact policy fields remain part of the renderer design, but their ownership
is fixed:

- `SmilesSyntaxFlags` selects grammar capabilities that are meaningful in both
  directions. Extended aromatic element symbols and extended bond symbols are
  examples: the parser accepts them and the renderer may emit them.
- `SmilesParsePolicy` holds genuinely input-directional behavior such as
  acceptance strictness, recovery, checking, and diagnostics.
- `SmilesRenderPolicy` holds genuinely output-directional behavior such as
  surface normalization, explicitness, and traversal choices.

The aggregate config is bidirectional; every individual setting need not be.
`SmilesIoFlags` would obscure this distinction. The shared bitset is therefore
named `SmilesSyntaxFlags`, while directional policy retains directional names.

### 1.1 Surface-normal round trips

Parsing and rendering are not required to reproduce every accepted spelling.
Let

```text
N_c(source) = render(parse(source, c), c)
```

for an IO configuration `c`. The required laws are:

```text
N_c(source) = source             when source is syntax-normal under c
N_c(N_c(source)) = N_c(source)   for every accepted source
```

For a non-normal input, the rendered spelling may differ. It must still encode
the same ordered molecular representation and must be accepted by the same IO
configuration.

This is deliberately weaker than byte-for-byte replay of arbitrary input. The
boundary therefore does not need to retain the original source or a full
concrete syntax tree merely to satisfy the round-trip law. It must retain enough
ordering information for rendering to select one stable surface-normal spelling.

### 1.2 Surface normalization is not graph canonicalization

The term *canonical SMILES* normally implies canonical graph labeling and atom
ordering. That is not the operation defined here. Documentation and APIs should
use *syntax-normal* or *surface-canonical* when referring to the fixed points of
`N_c`.

The parsed representation preserves input order relevant to rendering:

- atom encounter order;
- component order;
- bond or neighbor encounter order where branch and ring-closure choices need
  it.

Rendering may normalize redundant syntax within that order, including optional
bond notation, bracket use, numeric spelling, ring labels, and branch or closure
punctuation. It does not relabel atoms, reorder components, or search for a
graph-canonical traversal. Canonical atom ordering, if provided, is a separate
and explicitly requested graph operation.

The renderer design must verify that the ordered format value retains enough
information to make syntax-normal rendering deterministic. If it does not, the
boundary should gain the minimal missing ordering metadata rather than silently
perform graph canonicalization.

### 1.3 Boundary operations

The intended format-level shape is:

```rust
Smiles::parse_with(source, &io_config)
smiles.render_with(&io_config)
```

Conversion between `Smiles` and graph models remains a separate boundary.
Constructing a `Smiles` value from a graph model is fallible when the model is
not representable under the chosen IO configuration. The resulting format value
then renders to syntax accepted by that same configuration.

## 2. CXSMILES remains outside the Python round

Ordinary SMILES and CXSMILES remain separate future format boundaries. The
current Python workflow must not expose CX-specific capabilities merely because
the existing Rust parser configuration still carries them.

The ordinary Python surface hides:

- `CHEMAXON_EXTENSIONS`;
- `SKIP_UNKNOWN_CHEMAXON_TAGS`;
- the `CHEMAXON` preset;
- `chemaxon()`.

Its ordinary lenient preset must not enable CX syntax implicitly. Extended
aromatic elements and extended bond symbols remain ordinary SMILES syntax
capabilities and stay visible.

This round does not pull the full `CxSmiles` / `CxSmilesIoConfig` split from doc
153 into the Python workflow. That work still owns the CX boundary type, payload
representation, conversion semantics, renderer, and round-trip tests. The
Python binding simply does not stabilize the existing CX-specific Rust members.

## 3. Fixed raise semantics

TableIR-to-`MoleculeAst` raising remains fixed. It performs the format-to-model
interpretation recorded in doc 100 and has no public configuration object.
Failures are model-conversion errors, not evidence that a user-selectable raise
policy is missing.

## 4. `ChemistryModel` is semantic configuration

`ChemistryModel` selects the chemical interpretation used by resolution and
validation. It includes valence, aromaticity, and stereochemistry models and
their model data. It does not include operational decisions about consuming or
retaining source constraints.

The Python binding must expose the complete usable model vocabulary rather than
only `ChemistryModel::default()`:

- `ChemistryModel`;
- both `ValenceModel` variants, including owned `AtomTypeRegistry` and
  `ValenceTable` values;
- all `AromaticityModel` variants;
- `ElementScope` and `RingLimits`;
- `StereoModel`, `StereoKindModel`, and `InconsistencyPolicy`.

Python wrappers own their model data. Conversion to Rust may use owned `Cow`
variants internally; Rust borrowing and storage choices do not restrict the
Python model surface.

## 5. `ResolveConfig` is operational configuration

Resolution also has operational policy that is not part of the chemistry model.
The existing aromaticity and stereo configs govern such behavior as charge
delocalization and removal of source constraints. Removing those controls would
conflate a convenient default with the complete operation.

The composite resolve operation therefore receives a top-level config:

```rust
pub struct ResolveConfig {
    pub aromaticity: AromaticityResolveConfig,
    pub stereo: StereoResolveConfig,
}
```

Only resolver stages with operational policy need fields. New stage configs join
`ResolveConfig` when their operations acquire genuine policy; empty placeholder
configs are unnecessary.

The naming rule is operation-centered:

- `ResolveConfig`, not `ResolverConfig`;
- `AromaticityResolveConfig`, not `AromaticityResolverConfig`;
- `StereoResolveConfig`, not `StereoResolverConfig`.

The config describes the `resolve` operation, not the `Resolver` object that
executes it. `Resolver` may retain constructors such as `new(model)` for default
operational policy and `with_config(model, config)` for explicit policy without
changing that naming rule.

The graph ingestion boundary must accept both axes explicitly: a
`ChemistryModel` for semantic selection and a `ResolveConfig` for operational
behavior. Its unconfigured convenience path may supply higher-level defaults for
both.

## 6. Python resolved-SMILES operation

The complete Python operation has three independent keyword-only configuration
arguments:

```python
MoleculeAst.from_smiles(
    source,
    *,
    io_config=None,
    chemistry_model=None,
    resolve_config=None,
)
```

Omission selects the higher-level ordinary defaults:

- ordinary OpenSMILES IO policy;
- `ChemistryModel.default()`;
- `ResolveConfig.default()`.

The names expose the policy boundaries rather than flattening their fields into
method keyword arguments. Raise remains an internal fixed step between parsing
and resolution.

## 7. Consequences for the Python workflow plan

The remaining S4 work in doc 151 must be revised before `MoleculeAst.from_smiles`
is implemented:

1. update the current SMILES flag/config bindings to the paired IO design and
   hide CX-only members;
2. bind the complete chemistry-model vocabulary;
3. add and bind the top-level `ResolveConfig` vocabulary, including the
   operation-centered rename of stage configs;
4. extend Rust ingestion to accept explicit model and resolve configs;
5. only then add the configured Python `from_smiles` operation.

The staged implementation is maintained below; doc 151 forwards its remaining
resolved-SMILES work here.

## 8. Staged implementation plan

Every subitem ends with its affected Rust and Python suites green. Breaking
subitems include all caller migrations required to restore the build; none leave
a stage boundary red.

### S0 — Rust model and resolve-operation foundations

- **S0a — flatten `ValenceModel`** **Done**
  (`umol-graph/src/ops/model.rs`, `ops/valence/{atom_typing,counts}.rs`,
  `ops/{resolve,validate}/valence.rs`, and callers): remove `AtomTypingModel` and
  `CountsModel`. Replace them with the named-field variants
  `ValenceModel::AtomTyping { registry: Cow<'static, AtomTypeRegistry> }` and
  `ValenceModel::Counts { table: Cow<'static, ValenceTable> }`.
  `AtomTypingValence` borrows `AtomTypeRegistry` directly and `CountsValence`
  borrows `ValenceTable` directly. Migrate defaults, fixtures, integration tests,
  and every constructor without compatibility aliases. Exact dispatch and
  default-model tests cover borrowed defaults and owned custom data. **Breaking
  Rust API migration (red → green).** `[dep: —]`
- **S0b — operation-name the aromaticity config** **Done**
  (`umol-graph/src/ops/resolve/aromaticity.rs` and callers): rename
  `AromaticityResolverConfig` to `AromaticityResolveConfig`, preserving
  `delocalize_charge`, `reset_aromatic_valence`, defaults, and behavior. Migrate
  every constructor, import, and test. **Breaking Rust API migration (red →
  green).** `[dep: —]`
- **S0c — operation-name the stereo config** **Done**
  (`umol-graph/src/ops/resolve/stereo.rs` and callers): rename
  `StereoResolverConfig` to `StereoResolveConfig`, preserving
  `reset_stereo_constraints`, its default, and behavior. Migrate every
  constructor, import, and test. **Breaking Rust API migration (red → green).**
  `[dep: —]`
- **S0d — composite `ResolveConfig`** **Done**
  (`umol-graph/src/ops/resolve.rs`): add `ResolveConfig { aromaticity:
  AromaticityResolveConfig, stereo: StereoResolveConfig }`, `Default`, and
  structural equality. `Resolver::new(model)` delegates to default operational
  policy; `Resolver::with_config(model, config)` passes each branch to the
  corresponding stage. Tests prove exact default equivalence and independent
  propagation of every operational field. **Additive (green).** `[dep: S0b,
  S0c]`

S0 ends with a smaller semantic-model vocabulary and one complete operational
config for the composite resolve operation.

### S1 — Rust SMILES IO and ingestion migration

- **S1a — shared SMILES syntax configuration** **Done**
  (`umol-io/src/smiles/config.rs`, parser entry points, conformance tools,
  fuzz/benchmark callers, and workspace consumers): rename `SmilesParseFlags`
  to `SmilesSyntaxFlags`, `parse_flags` to `syntax_flags`, and
  `with_parse_flags` to `with_syntax_flags`. `OPENSMILES` remains zero;
  `LENIENT` becomes the ordinary-SMILES union of extended aromatic and extended
  bond syntax and no longer enables CX. Existing explicit CX bits and presets
  remain Rust-only pending doc 153 rather than being moved into a new boundary
  here. Tests pin bit positions, ordinary presets, OR composition, display, and
  explicit CX opt-in. **Breaking Rust API and preset migration (red → green).**
  `[dep: —]`
- **S1b — explicit resolve policy at interpretation and ingestion boundaries** **Done**
  (`umol-graph/src/ingest.rs`, current resolved MOL helpers, and callers): make
  `Interpret::interpret` accept both `&ChemistryModel` and `&ResolveConfig`;
  rename `MoleculeIngestError` to `MoleculeInterpretationError`; extend
  `ingest_smiles_with` / `ingest_smiles_bytes_with` and every public resolved
  format path that already accepts a chemistry model with explicit resolve
  policy. Unconfigured convenience operations use `ChemistryModel::default()`
  and `ResolveConfig::default()`. Migrate all workspace callers in the same
  subitem. End-to-end tests distinguish semantic-model selection from each
  operational resolve option and prove default parity. **Breaking Rust API
  migration (red → green).** `[dep: S0a, S0d, S1a]`

S1 ends with Rust exposing every policy required by the eventual Python input
operation while retaining the paired `SmilesIoConfig` direction.

### S2 — Rust model-data value semantics

- **S2a — `AtomTypeRegistry` structural equality**
  (`umol-graph/src/ops/valence/registry.rs`): add exact structural equality over
  registry entries; the content hash remains metadata and is not used as an
  equality substitute. Tests distinguish equal reconstruction, differing
  patterns, and differing element/charge buckets. **Additive (green). Done.**
  `[dep: —]`
- **S2b — `ValenceEntry` / `ValenceTable` structural equality**
  (`umol-graph/src/ops/valence/table.rs`): add exact value equality for entries
  and tables, including target covalences and aromatic valences. Tests
  distinguish equal reconstruction, missing elements, and changed entry
  values. **Additive (green). Done.** `[dep: —]`
- **S2c — chemistry-model structural equality**
  (`umol-graph/src/ops/model.rs`): add structural `PartialEq` throughout
  `ValenceModel`, `AromaticityModel`, `StereoModel`, and `ChemistryModel`, using
  the exact registry/table equality from S2a/S2b and ordinary floating-point
  equality for HMO thresholds. Tests vary every branch and aggregate field.
  **Additive (green). Done.** `[dep: S0a, S2a, S2b]`

### S3 — corrected Python SMILES IO values

- **S3a — `SmilesSyntaxFlags`** (`umol-py/src/smiles.rs`, registration,
  exports, and installed tests): replace the current `SmilesParseFlags` binding
  with immutable `SmilesSyntaxFlags`. Expose `EXTENDED_AROMATICS`,
  `EXTENDED_BONDS`, `OPENSMILES`, and the ordinary `LENIENT` preset with bitwise
  OR, validated bit construction, equality, and repr. CX bits and presets are
  rejected by the Python constructor and are neither registered nor exported.
  Conversion tests cover every exposed bit and preset without implying that the
  Rust CX members were removed.
  **Breaking Python API correction (red → green). Done.**
  `[dep: S1a]`
- **S3b — corrected `SmilesIoConfig` binding**
  (`umol-py/src/smiles.rs`, registration, exports, and installed tests): expose
  immutable `SmilesIoConfig` through `opensmiles()`, ordinary `lenient()`, and
  keyword-only `with_syntax_flags(syntax_flags=...)`, plus a detached read-only
  `syntax_flags` value, structural equality, and repr. Remove Python
  `chemaxon()`, `parse_flags`, and `with_parse_flags`; do not expose lint or
  future render-policy internals before their contracts are settled. Rust/Python
  conversions cover all ordinary configs and an arbitrary OR composition.
  **Breaking Python API correction (red → green). Done.** `[dep: S3a]`

### S4 — Python valence-model data

Each public value introduced in this stage is registered in `_native` and
exported from `umol` in its own subitem.

- **S4a — `AtomTypeRegistry`**
  (`umol-py/src/model/valence.rs`, registration, and exports): bind an immutable
  owned registry with `default()`, `from_atoms(...)`, `from_toml(...)`, exact
  equality, content-hash inspection, read-only pattern lookup, and
  stable repr. Invalid TOML raises `ValueError`; returned atoms are
  detached Python values. Conversion tests cover borrowed-default and owned
  custom registries. **Additive (green). Done.** `[dep: S2a]`
- **S4b — `ValenceEntry`** (`umol-py/src/model/valence.rs`): bind immutable
  target-covalence and aromatic-valence sequences with a keyword-only
  constructor, detached getters, normalization matching Rust table insertion,
  equality, and repr. Tests cover empty and multi-state entries. **Additive
  (green). Done.** `[dep: S2b]`
- **S4c — `ValenceTable`**
  (`umol-py/src/model/valence.rs`, registration, and exports): bind an immutable
  owned table with `default()`, keyword-only construction from an element-to-
  entry mapping, `from_toml(...)`, exact equality, content-hash inspection,
  read-only lookup, and repr. Invalid TOML raises `ValueError`; returned entries
  are detached. Conversion tests cover borrowed-default and owned custom tables.
  **Additive (green). Done.** `[dep: S2b, S4b]`
- **S4d — `ValenceModel`** (`umol-py/src/model/valence.rs`): bind only the
  direct `AtomTyping { registry }` and `Counts { table }` variants, with
  keyword-only payloads, structural equality, variant repr, and separate
  `from_rust` / `to_rust` implementations. No Python or Rust intermediate model
  wrappers are reintroduced. Tests cover owned round trips for both variants.
  **Additive (green). Done.** `[dep: S0a, S2c, S4a, S4c]`

### S5 — Python aromaticity and stereo models

Each public value introduced in this stage is registered in `_native` and
exported from `umol` in its own subitem.

- **S5a — `ElementScope`** (`umol-py/src/model.rs`): bind `Any()` and
  `AllowList(elements)` as immutable values with detached element sequences,
  equality, repr, and Rust round trips. **Additive (green). Done.** `[dep: —]`
- **S5b — `RingLimits`** (`umol-py/src/model/aromaticity.rs`): bind every field
  with a keyword-only constructor carrying the Rust defaults, immutable getters,
  equality, repr, and Rust round trips. Boundary tests cover zero/invalid Python
  integer conversion and nondefault fused-ring limits. **Additive (green).
  Done.** `[dep: —]`
- **S5c — `AromaticityModel`** (`umol-py/src/model/aromaticity.rs`): bind
  `HueckelRule`, `Hmo`, and `Clar` with keyword-only variant fields, plus
  `daylight()`, `mdl()`, and `permissive()` presets. Implement structural
  equality, variant repr, and separate conversions; tests cover every variant,
  preset, scope, ring-limit, and HMO threshold path. **Additive (green). Done.**
  `[dep: S2c, S5a, S5b]`
- **S5d — `InconsistencyPolicy`** (`umol-py/src/model/stereo.rs`): bind `Keep`,
  `Strip`, and `Error` as a fieldless immutable enum with equality, hashing,
  repr, and Rust round trips. **Additive (green). Done.** `[dep: —]`
- **S5e — `StereoKindModel`** (`umol-py/src/model/stereo.rs`): bind `scope` and
  `fluxionality` with a keyword-only constructor, detached ownership, equality,
  repr, and Rust round trips. **Additive (green). Done.** `[dep: S5a]`
- **S5f — `StereoModel`** (`umol-py/src/model/stereo.rs`): bind the per-kind
  model mapping, `para_stereo`, `max_iterations`, and inconsistency policy with
  keyword-only construction and `default()`. Missing stereo kinds map to Rust
  `None`; getters return detached mappings. Tests cover the exact default kind
  set, every supported `StereoKind`, disabled kinds, all scalar fields,
  equality, repr, and Rust round trips. **Additive (green). Done.** `[dep: S2c,
  S5d, S5e]`

### S6 — Python aggregate semantic and operational configs

- **S6a — `ChemistryModel`** (`umol-py/src/model.rs`, registration, and
  exports): bind `valence`, `aromaticity`, and `stereo` as an immutable config
  with keyword-only construction, `default()`, detached getters, structural
  equality, repr, and owned Rust round trips. Tests pin the full default and
  independently replace every branch. **Additive (green). Done.** `[dep: S2c,
  S4d, S5c, S5f]`
- **S6b — `AromaticityResolveConfig`** (`umol-py/src/resolve.rs`): bind
  `delocalize_charge` and `reset_aromatic_valence` with a keyword-only
  constructor carrying Rust defaults, immutable getters, equality, repr, and
  Rust round trips. **Additive (green). Done.** `[dep: S0b]`
- **S6c — `StereoResolveConfig`** (`umol-py/src/resolve.rs`): bind
  `reset_stereo_constraints` with a keyword-only constructor carrying the Rust
  default, immutable getter, equality, repr, and Rust round trips. **Additive
  (green). Done.** `[dep: S0c]`
- **S6d — `ResolveConfig`** (`umol-py/src/resolve.rs`, registration, and
  exports): bind the aromaticity and stereo branches with keyword-only
  construction, `default()`, detached getters, structural equality, repr, and
  Rust round trips. Tests pin default composition and independent replacement of
  both branches. **Additive (green).** `[dep: S0d, S6b, S6c]`

### S7 — fully configured Python SMILES ingestion

- **S7a — `MoleculeAst.from_smiles`**
  (`umol-py/src/molecule.rs`, `src/error.rs`, installed tests, and exports): add
  `from_smiles(source, *, io_config=None, chemistry_model=None,
  resolve_config=None)`. Omitted values select the three higher-level defaults;
  explicit values lower through their owned wrappers and call the fully
  configured Rust ingestion path. Preserve the settled syntax,
  model-conversion, contradiction, underdetermination, and execution-error
  mapping. Rust/PyO3 and installed tests assert exact determined molecules,
  ordinary syntax configurations, both valence models, representative
  aromaticity/stereo model changes, every operational resolve field, every
  reachable error category, keyword-only rejection, and detached ownership.
  This consumes the already-complete semantic exception mapping from doc 151
  S4b. **Additive (green).** `[dep: S1b, S3b, S6a, S6d]`

S7 is the completed resolved-SMILES Python deliverable.

### Critical path and deferred follow-ons

The three prerequisite branches join at S7a:

```text
S1a → S3a → S3b ───────────────────────────────┐
S0a + S2a/S2b → S2c → S4/S5 → S6a ───────────┼→ S7a
S0b/S0c → S0d → S1b and S6b/S6c → S6d ───────┘
```

No stage in this plan is deferrable relative to the configured Python ingestion
contract. The SMILES renderer and concrete `SmilesParsePolicy` /
`SmilesRenderPolicy` fields remain a separate follow-on design; no empty policy
types are introduced here. The `CxSmiles` boundary and `CxSmilesIoConfig` remain
in doc 153, while the present Python surface deliberately hides CX-only Rust
members.
