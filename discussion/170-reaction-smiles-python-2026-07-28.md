# Reaction SMILES import in Python

Status: **Completed**
Date: 2026-07-28
Relates: [151](151-python-molecule-workflows-2026-07-13.md),
[153](153-format-parsing-outstanding-tasks-2026-07-18.md),
[155](155-smiles-io-and-resolve-configuration-2026-07-19.md),
[162](162-common-subgraph-algs-2026-07-25.md)

## Scope

This document tracks two remaining Python reaction-interface changes:

1. ingest reaction SMILES into `ReactionAst`;
2. replace the raw algorithm argument of Python reaction composition with a
   Python-only operation config.

The reaction SMILES parser already has an explicit external-format boundary,
but the parsed value cannot yet be interpreted as the graph reaction model.
The composition change is smaller and independent. It is recorded here because
both changes complete high-level Python reaction workflows without changing
the algorithmically transparent Rust APIs.

## Existing reaction SMILES boundary

`ReactionSmiles` is the parsed semantic value of reaction SMILES and privately
wraps `table_ir::Reaction`. It supports configured text and byte parsing and
read-only or consuming TableIR access:

```text
reaction SMILES text
    -> ReactionSmiles
    -> table_ir::Reaction
```

The TableIR value contains reactants, products, agents, and atom-map classes.
The pipeline stops there. There is no `Interpret for ReactionSmiles`, no
`ingest_reaction_smiles*` operation in `umol-graph`, and therefore no
model-ingestion operation for Python to bind.

`ReactionDsl` is not an intermediate representation for this conversion. It is
the faithful internal DSL boundary for an existing `ReactionAst`. Its
restructuring simplifies exact construction, rendering, and tests, while
external reaction SMILES still enters through TableIR.

## Interpretation into `ReactionAst`

Add a graph-layer interpretation path parallel to molecule SMILES:

```text
ReactionSmiles
    -> interpret reactants as MoleculeAst
    -> interpret products as MoleculeAst
    -> derive atom Correspondence from atom-map classes
    -> ReactionAst::from_sides(reactants, products, correspondence)
```

Both molecule sides use the same `ChemistryModel` and `ResolveConfig`. Raising
and resolution do not renumber atoms, so the atom ids collected by the parser
remain valid after side interpretation. The conversion should nevertheless
test this contract directly rather than relying on it implicitly.

Implement `Interpret for ReactionSmiles` with `ReactionAst` as its output.
Add the same four graph-layer convenience shapes as molecule SMILES:

```rust
ingest_reaction_smiles
ingest_reaction_smiles_bytes
ingest_reaction_smiles_with
ingest_reaction_smiles_bytes_with
```

The unconfigured forms use `SmilesIoConfig::opensmiles()`,
`ChemistryModel::default()`, and `ResolveConfig::default()`. The configured
forms require all three Rust configuration references. No parsing,
model-conversion, or resolution logic is duplicated in `umol-py`.

The layers remain independently usable and visible in Rust:

1. `ReactionSmiles::parse*` performs syntax-only parsing;
2. `Interpret for ReactionSmiles` converts the parsed value under an explicit
   chemistry model and resolve config;
3. `ingest_reaction_smiles*` composes parsing and interpretation for the common
   text-to-model workflow;
4. Python `ReactionAst.from_reaction_smiles` lowers its keyword-only configs
   and calls `ingest_reaction_smiles_with`.

The top-level Python method must not duplicate or hide a second conversion
pipeline. Errors retain their phase and source across all four layers.

## Atom-map semantics

`table_ir::Reaction::atom_mapping` groups the reactant and product atom ids for
each numerical atom-map class. Convert it to the partial atom correspondence
used by `ReactionAst::from_sides` as follows:

- exactly one reactant and one product atom with a class produce one matched
  pair;
- a class present on only one side produces no pair, so that atom remains an
  explicit deletion or addition;
- more than one atom with the same class on either side is valid reaction
  SMILES but is an error for the strict single-`ReactionAst` projection because
  it does not select one partial bijection;
- absence of atom-map classes produces the empty correspondence and therefore
  a reaction containing deletions and additions rather than an inferred map.

The operation must not infer atom mappings or select a mapping algorithm
silently. Mapping inference is a separate future operation with its own
algorithm configuration. Atom-map class numbers are boundary annotations used
to construct the correspondence; they are not retained as atom fields in
`ReactionAst`.

Reaction-SMILES map numbers denote equivalence classes rather than necessarily
unique atom pairs. Repeated class members can express chemically equivalent
atoms, alternative mechanisms, or incomplete knowledge about atom provenance.
The TableIR vectors preserve this set-valued format semantics.

Repeated atom-map classes are rejected in the strict reaction
raise/model-conversion phase, before either side is resolved and before deltas
are derived. They are not parser syntax errors: the reaction SMILES and its
TableIR representation remain valid, but one `ReactionAst` requires one
partial-bijection atom map.

A separate future plural projection may enumerate every partial bijection
compatible with the reaction-SMILES equivalence classes and return the
corresponding vector of reactions. This is a lossless expansion of ambiguity,
not mapping inference: it retains all compatible choices instead of selecting
or optimizing one. The exact treatment of unequal class cardinalities and
partially specified classes must be designed with that operation. It is not
part of `from_reaction_smiles`, whose singular return type retains the strict
interpretation.

## Agents

`ReactionAst` represents a left-hand molecule plus deltas and has no agent
channel. Agents cannot be folded into the left-hand side without changing
their meaning, and silently discarding them would make import lossy without an
observable diagnostic.

Reject a nonempty reaction-SMILES agent section in the same
raise/model-conversion phase, before resolution. A future reaction-record
boundary may preserve agents and other record-level data, but that type is not
part of this work. Empty agent sections, including the common
`reactants>>products` spelling, remain accepted.

## Errors

Reaction ingestion must preserve the existing Python error taxonomy while
adding side context:

| Failure | Python exception |
| --- | --- |
| reaction SMILES syntax | `ParseError` |
| reactant/product TableIR-to-AST conversion | `ModelConversionError` |
| reactant/product resolver contradiction | `ContradictionError` |
| reactant/product remains underdetermined | `UnderdeterminedError` |
| reactant/product resolver execution failure | `RuntimeError` |
| ambiguous multi-atom map class | `ModelConversionError` |
| nonempty agents | `ModelConversionError` |

The diagnostic for a side interpretation failure identifies whether the
reactants or products failed while preserving the underlying error as its
source. Rust should use one reaction-input error rather than duplicating every
molecule interpretation category into reactant and product variants.

The public Rust errors are `ReactionInterpretationError` for parsed-value
interpretation and `ReactionSmilesInputError` for the combined text-to-model
operation. The strict mapping error describes an ambiguous class that cannot
be projected into one `ReactionAst`, rather than describing valid
reaction-SMILES syntax as a duplicate or malformed map.

## Python API

Expose the configured convenience operation directly on `ReactionAst`:

```python
ReactionAst.from_reaction_smiles(
    source,
    *,
    io_config=None,
    chemistry_model=None,
    resolve_config=None,
)
```

The three options have the same Python wrapper types, keyword-only behavior,
ownership semantics, and defaults as `MoleculeAst.from_smiles`. The more
explicit method name retains “reaction SMILES” as the external format name and
does not conflict with the molecule operation.

No Python `ReactionSmiles` boundary wrapper is required for this workflow. It
can be added later only if Python needs syntax-only parsing or repeated
interpretation of the same parsed input under different chemistry models.
The absence of that Python wrapper does not narrow the Rust surface:
`ReactionSmiles::parse*`, `Interpret`, and `ingest_reaction_smiles*` remain
public operations.

## Reaction composition configuration

Rust keeps the existing algorithmically transparent signature:

```rust
ReactionAst::compose(
    &self,
    other: &ReactionAst,
    algorithm: CommonSubgraphEnumerationAlgorithm,
)
```

Python replaces the method-level `algorithm` keyword with the Python-only
`ReactionCompositionConfig`:

```python
config = ReactionCompositionConfig(
    common_subgraph_enumeration_algorithm=(
        CommonSubgraphEnumerationAlgorithm.DirectBacktracking()
    )
)

first.compose(second, config=config)
```

`ReactionAst.compose(other, *, config=None)` uses
`ReactionCompositionConfig.default()`. The default is
`DirectBacktracking()`. The paired benchmark cases in doc 162 measured the
direct implementation approximately 1.7--2.3 times faster than the modular
product implementation while returning the same complete ordered results.
These preliminary measurements justify the current high-level default but are
not a permanent benchmark policy. A later pass must review all Python and
other high-level defaults against representative benchmarks.

This supersedes only the direct Python selector decision in doc 151 S9q.
`CommonSubgraphEnumerationAlgorithm` remains exported so callers can construct
an explicit config, and the Rust operation continues to require a selector.
The config name, Python-only placement, keyword-only method argument, and
`DirectBacktracking()` default are settled.

This is also the binding rule for additional algorithm-selecting operations
from `umol-ast`: as those operations gain Python bindings, each receives a
dedicated Python operation-config wrapper instead of exposing a bare algorithm
keyword. The Rust APIs remain algorithmically transparent and continue to
accept their selectors directly; no config is added to `umol-ast` solely for
Python uniformity.

## Verification requirements

The reaction ingestion work must cover:

- mapped preserved atoms and attribute changes;
- unmapped deletions and additions;
- one-sided atom-map classes;
- duplicate classes on each side;
- completely unmapped reactions;
- rejection of nonempty agents;
- syntax and each reactant/product interpretation error category;
- propagation of non-default IO, chemistry, and resolve configuration;
- the exact keyword-only Python signature and detached returned objects.

The composition migration must cover:

- exact config construction, equality, getters, and repr;
- `DirectBacktracking()` as the omitted-config default;
- explicit propagation of both enumeration algorithms;
- rejection of the removed method-level `algorithm` keyword;
- unchanged complete, deterministic composition results.

## Staged implementation plan

### S0 — Rust interpretation foundation **Done**

S0 leaves the existing molecule-SMILES API and behavior unchanged while
factoring the complete TableIR-molecule interpretation sequence for reuse by
the two reaction sides.

#### S0a — Share the TableIR molecule interpretation kernel **Done**

- **Module:** `umol-graph/src/ingest.rs`
- **Work:** Move the existing TableIR-molecule raise-and-resolve sequence into
  one module-local operation used by `Interpret for Smiles` and, subsequently,
  by both sides of `Interpret for ReactionSmiles`. The shared operation covers
  the whole conversion sequence and returns `MoleculeInterpretationError`; do
  not split isolated conversion steps into additional helpers or expose new
  public API.
- **Compatibility:** Internal refactor; no public or behavioral change.
- **Dependencies:** None.
- **Tests:** Keep the existing exact molecule interpretation, configured-model,
  configured-resolution, and error-propagation tests green. Add no test of the
  private extraction itself; exercise it through the existing public
  `Interpret for Smiles` surface.

#### S0b — Add `ReactionInterpretationError` **Done**

- **Module:** `umol-graph/src/ingest.rs`
- **Work:** Add the approved public error type for interpreting a parsed
  `ReactionSmiles`. It must retain reactant/product context around an
  underlying `MoleculeInterpretationError` and represent the two strict
  reaction-model conversion failures: a nonempty agent section and an atom-map
  class that cannot be projected into one partial bijection. Keep the latter
  failures distinct from reaction-SMILES syntax errors. The variants are
  `Reactants`, `Products`, `AmbiguousAtomMapClass`, and `AgentsUnsupported`.
- **Compatibility:** Additive.
- **Dependencies:** None.
- **Tests:** Use table cases to verify exact display text, error sources, side
  context, and the classification of agent and ambiguous-map failures as
  interpretation rather than syntax failures.

#### S0c — Add `ReactionSmilesInputError` **Done**

- **Module:** `umol-graph/src/ingest.rs`
- **Work:** Add the approved public error type for the composed parse-and-
  interpret operation. It must preserve a reaction parser error as the syntax
  branch and `ReactionInterpretationError` as the interpretation branch,
  including their source chains.
- **Compatibility:** Additive.
- **Dependencies:** S0b.
- **Tests:** Verify exact display text and sources for representative syntax,
  side-model, side-resolution, ambiguous-map, and agent failures. Test the
  public error categories rather than duplicating parser or resolver internals.

### S1 — Rust reaction ingestion **Done**

#### S1a — Implement strict `Interpret for ReactionSmiles` **Done**

- **Modules:** `umol-graph/src/ingest.rs`,
  `umol-graph/tests/` if the existing ingest test module becomes too large.
- **Work:** Implement `Interpret<Output = ReactionAst>` for `ReactionSmiles`.
  Reject nonempty agents and non-singular atom-map classes in the strict
  raise/model-conversion phase before either side is resolved. Interpret the
  reactants and products through the shared S0a operation with the same
  `ChemistryModel` and `ResolveConfig`, preserving side context on failure.
  Convert every one-to-one map class to a pair, leave one-sided and absent
  classes unmatched, construct the partial atom `Correspondence`, and call
  `ReactionAst::from_sides`. Do not infer, optimize, or select atom mappings.
- **Compatibility:** Additive.
- **Dependencies:** S0a, S0b.
- **Tests:** Compare against independently written exact reaction-AST/DSL
  expectations for mapped preserved atoms, mapped attribute and bond changes,
  one-sided classes, unmapped deletions and additions, and completely unmapped
  reactions. Add separate exact failures for multiplicity on the reactant
  side, product side, and both sides of a map class; nonempty agents; and every
  reactant/product interpretation category. Verify explicitly that parser atom
  ids still identify the same atoms after side interpretation. Expected values
  must not be constructed by the implementation path under test.

#### S1b — Add public `ingest_reaction_smiles*` operations **Done**

- **Module:** `umol-graph/src/ingest.rs`
- **Work:** Add `ingest_reaction_smiles`,
  `ingest_reaction_smiles_bytes`, `ingest_reaction_smiles_with`, and
  `ingest_reaction_smiles_bytes_with`. Match the molecule-SMILES division
  exactly: unconfigured operations select the documented high-level defaults,
  while configured operations require explicit references to
  `SmilesIoConfig`, `ChemistryModel`, and `ResolveConfig`. Compose the existing
  parser with S1a and return `ReactionSmilesInputError`; do not duplicate
  interpretation logic.
- **Compatibility:** Additive.
- **Dependencies:** S0c, S1a.
- **Tests:** Verify text/byte parity, default/convenience parity, explicit
  propagation of non-default IO/model/resolve settings, and representative
  syntax and interpretation failures through each applicable entry point.

### S2 — Python reaction-SMILES ingestion **Done**

#### S2a — Map reaction input errors to Python exceptions **Done**

- **Module:** `umol-py/src/error.rs`
- **Work:** Add the conversion from `ReactionSmilesInputError` to the existing
  Python exception hierarchy. Preserve reactant/product context in messages
  while mapping the underlying category to `ModelConversionError`,
  `ContradictionError`, `UnderdeterminedError`, or `RuntimeError`. Map syntax
  to `ParseError`, and map ambiguous-class and agent failures to
  `ModelConversionError`.
- **Compatibility:** Additive.
- **Dependencies:** S0c.
- **Tests:** Use table cases covering every row of the error table above,
  including both sides where side context matters. Assert the exact Python
  exception type and message.

#### S2b — Add `ReactionAst.from_reaction_smiles` **Done**

- **Modules:** `umol-py/src/reaction.rs`,
  `umol-py/tests/test_reaction.py`, `umol-py/tests/test_import.py`
- **Work:** Add the static method with the settled keyword-only
  `io_config`, `chemistry_model`, and `resolve_config` parameters. Lower the
  optional Python wrappers in the same way as `MoleculeAst.from_smiles`, call
  `ingest_reaction_smiles_with`, and return a detached Python-owned
  `ReactionAst`. Do not add a Python `ReactionSmiles` wrapper or another
  conversion pipeline.
- **Compatibility:** Additive.
- **Dependencies:** S1b, S2a.
- **Tests:** Exercise mapped, one-sided, and unmapped reactions through Python;
  non-default values for each of the three configs; every public exception
  category; the exact keyword-only signature; rejection of positional config
  arguments; and ownership after the source/config objects leave scope. Update
  the import/signature inventory with a fixed expected signature.

### S3 — Python reaction composition configuration **Done**

This stage is independent of S0--S2 and may be implemented in parallel with
the reaction-ingestion path.

#### S3a — Add `ReactionCompositionConfig` **Done**

- **Modules:** `umol-py/src/reaction.rs`, `umol-py/src/lib.rs`,
  `umol-py/tests/test_reaction.py`, `umol-py/tests/test_import.py`
- **Work:** Add and register a frozen, equality-comparable Python config with
  the single `common_subgraph_enumeration_algorithm` field. Provide the usual
  keyword-only constructor, getter, `default()`, and repr. Implement the
  complete `to_rust`/`from_rust` boundary conversion separately for this type.
  Its default is
  `CommonSubgraphEnumerationAlgorithm.DirectBacktracking()`.
- **Compatibility:** Additive.
- **Dependencies:** None.
- **Tests:** Verify construction, default construction, equality/inequality,
  getter, exact repr, both Rust conversion directions, module export, and the
  direct-backtracking default.

#### S3b — Migrate `ReactionAst.compose` to the config **Done**

- **Modules:** `umol-py/src/reaction.rs`,
  `umol-py/tests/test_reaction.py`, `umol-py/tests/test_import.py`, and all
  Python call sites found by repository-wide search.
- **Work:** Replace the Python-only
  `compose(other, *, algorithm=...)` signature with
  `compose(other, *, config=None)`. Lower the supplied or default config and
  pass its selector to the unchanged, algorithmically transparent Rust
  `ReactionAst::compose`. Remove the old method-level keyword rather than
  supporting two configuration paths.
- **Compatibility:** Breaking Python signature migration; all in-repository
  callers change in this subitem so the stage ends green.
- **Dependencies:** S3a.
- **Tests:** Verify omitted-config behavior uses direct backtracking, explicit
  configs propagate both algorithms, the removed `algorithm` keyword is
  rejected, the public signature is exact, and complete deterministic results
  and existing snapshots are unchanged apart from the selected default where
  applicable.

### Build and verification order

The reaction-ingestion critical path is
`S0a + S0b -> S0c -> S1a -> S1b -> S2a -> S2b`. The independent composition
path is `S3a -> S3b`. Each subitem includes its focused Rust or Python tests;
each stage ends with formatting and the affected crate suites green. After
both paths complete, activate `umol-py/.venv`, rebuild the extension, run the
full Python suite, then run workspace tests and clippy.

The plural ambiguity-preserving reaction projection and the benchmark-driven
review of all high-level algorithm defaults remain explicitly deferred. No
subitem in the plan above is deferrable for the singular
`ReactionAst.from_reaction_smiles` or configured Python composition surfaces.

## References

- Daylight, [*SMILES Theory: Reaction Atom
  Maps*](https://daylight.com/dayhtml/doc/theory/theory.smiles.html#RTFToC40):
  reaction maps are equivalence classes with no completeness or uniqueness
  requirement.
- OpenSMILES, [*OpenSMILES
  specification*](https://opensmiles.org/opensmiles.pdf): atom classes are
  application-defined integers and may label multiple atoms.
