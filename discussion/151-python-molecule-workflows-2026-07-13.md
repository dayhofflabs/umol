# Python molecule workflows: SMILES, fingerprints, and reactions

Status: **Active design / general implementation plan**

Date: 2026-07-13

Updated: 2026-07-15

Relates: 126 (fingerprints), 131–135 (reaction AST and application), 137 and 140
(Python binding conventions), 148 (operation boundaries and validation), 150
(completed reaction bindings)

## Goal

The first workflow-oriented Python deliverable makes the existing `MoleculeAst`
binding useful without reimplementing chemistry in Python. It has three required
capabilities:

1. parse a SMILES string into a resolved `MoleculeAst`;
2. compute, inspect, compare, and export molecule and reaction fingerprints as
   separate result families;
3. parse and apply reactions while retaining derivations and correspondences.

The reaction value and application surface is implemented. The remaining round
adds resolved SMILES parsing and molecule/reaction fingerprints, refines workflow
configuration, and closes the operation-boundary error gaps shared by all three
capabilities.

`ReactionSpanAst`, the span DSL, unresolved low-level SMILES parsing, general
chemistry-model configuration, DRFP, and BRIDGIT are outside this round.

## Layering: algorithm transparency and workflow configuration

Low-level algorithmic APIs follow the algorithmic-transparency policy: when an
operation has multiple implementations, the caller passes the algorithm enum
explicitly. `umol-graph-core` is the canonical layer for this policy. Algorithmic
entry points in other low-level modules follow the same rule where algorithm
selection is itself part of the operation. These layers do not acquire a default
algorithm merely to shorten higher-level calls.

Workflow APIs in `umol-graph` and format APIs in `umol-io` instead accept
configuration objects. A configuration object owns the selected algorithms and
their parameters together with other policy that belongs to the operation. This
is the appropriate higher-level boundary for two reasons:

1. new settings can be added without lengthening every method signature;
2. the configuration can define a reasonable workflow default without imposing
   that default on `umol-graph-core` or another algorithmic primitive layer.

The crate boundary is not the only distinction: a low-level algorithmic primitive
inside a higher-level crate may still take an explicit algorithm enum, while the
crate's public workflow entry point accepts a config and performs the explicit
dispatch internally. Within `umol-io`, this likewise distinguishes any
algorithm-selecting parser primitive from the configured format workflow that
contains it.

The Python API is a workflow boundary and follows the configuration-object side
of this policy. In particular:

- SMILES parsing accepts `SmilesIoConfig` and uses its default when omitted;
- molecular fingerprint generation accepts a fingerprint configuration that
  selects the featurizer and its parameters;
- reaction application accepts an application configuration whose default selects
  the standard matching behavior.

The underlying algorithm enums remain available where needed to construct or
inspect configs, and are still passed explicitly from the binding implementation
to the Rust algorithmic methods. Python does not create a second, unrelated set
of algorithm choices.

### Computation API organizing principle

The public unit is an operation: one semantic computation over one kind of input
with one coherent configuration space. Python computation APIs are organized by
the following rules:

1. Distinct operations receive distinct method invocations and distinct config
   types. SMILES, MOL, and SDF parsing are separate operations even though all
   produce molecules; hashed, pattern, and structural molecular fingerprints are
   likewise separate operations.
2. Alternative algorithms remain inside one config when they implement the same
   operation with the same input and result contract. WL, ECFP, and Morgan are
   variants selected by `HashedFingerprintConfig`, not three methods and three
   unrelated configs.
3. A config contains the algorithm selection, algorithm parameters, and policy
   that may expand together. The method does not grow one keyword argument per
   new setting, and low-level algorithm enums do not acquire defaults for the
   convenience of the workflow layer.
4. A method has a fixed result type where practical. When two result encodings
   arise from the same configured computation, separate methods may share one
   config: `hashed_fingerprint` and `counted_hashed_fingerprint` both consume
   `HashedFingerprintConfig` but return different feature-set types.
5. When one operation intrinsically has variant result semantics, it returns one
   dedicated operation-specific sum type rather than a broad cross-operation
   union. `combined_fingerprint` returns `ReactionCombinedFingerprint` with
   Difference and DisjointUnion variants.
6. Molecule and reaction computations have separate configs and result types even
   when they reuse the same lower-level featurizer or storage substrate. DRFP and
   BRIDGIT therefore become separate future reaction operations, not variants of
   a molecule hashed-fingerprint config.
7. A config is optional only when its default represents the operation's minimal
   or ordinary baseline and other settings layer behavior onto that baseline.
   `SmilesIoConfig` satisfies this rule. Sibling fingerprint definitions such as
   WL, ECFP, and Morgan do not form such a hierarchy, so
   `HashedFingerprintConfig` is required even though each selected variant may
   default its own parameters.

New functionality joins an existing method/config only when it preserves that
operation's input, semantics, and result contract. Otherwise it adds a new
operation pair. This keeps expansion systematic without creating one method per
minor algorithm option or one universal config/result union.

## Existing Python foundation

`umol-py/src/molecule.rs` wraps `umol_ast::ast::MoleculeAst` by value and exposes:

- construction from every entity family;
- live atom, bond, overlay, and stereo views;
- value equality and a compact representation;
- crate-private immutable and mutable access to the Rust AST;
- a crate-private `from_inner` constructor for sibling bindings that receive an
  owned Rust `MoleculeAst`.

Python binding types own their corresponding Rust values. Entity views retain
their molecule owner and route mutations back through it. New public classes are
registered in `umol-py/src/lib.rs`, re-exported from `python/umol/__init__.py`, and
listed in `__all__`.

## Current reaction deliverable

The reaction half of the workflow is implemented through doc 150:

- the complete resolved delta hierarchy and `Deltas` container;
- live `ReactionAst` lhs and delta components;
- reaction DSL parse/render, canonicalization, reverse, `from_sides`, and compose;
- return-only `Correspondence`, `MoleculeCorrespondence`, and
  `ReactionDerivation` values;
- all eight molecule-entity correspondence families plus the atom map;
- reverse, chain, and conversion of a derivation back to `ReactionAst`;
- an owned one-shot application iterator.

`ReactionAst.apply` snapshots the reaction and host, eagerly enumerates
graph-and-overlay correspondences, and lazily constructs one owned
`ReactionDerivation` per successful `next`. Successive derivations and every
observed molecule, correspondence, and recovered reaction are detached values.
The iterator class itself is not a package export.

The current Python signature passes `SubgraphIsomorphismAlgorithm` directly. The
workflow-configuration policy above refines it to:

```python
reaction.apply(host, config=None)

ReactionApplicationConfig(
    strategy=SubstructureMatchAlgorithm.GraphAndOverlays(),
    algorithm=SubgraphIsomorphismAlgorithm.Vf2Rdkit(),
)
```

Omitting the config selects these high-level workflow defaults. Both fields
remain configurable: `GraphAndOverlays` is the faster ordinary strategy, while
`Incidence` is required when connectivity is carried only by multicenter or
other overlays. `Vf2Rdkit` is the standard subgraph-isomorphism backend for the
Python workflow. This does not introduce a low-level default: the config still
passes its selected strategy and algorithm explicitly to the Rust matcher.

One reaction correctness issue remains. Application currently skips every
`ApplyError`, conflating expected rejection of one embedding with invalid reaction
input and internal lowering/transaction failures. This round must preserve lazy
emission while distinguishing those outcomes; it does not add `apply_at` or a
public diagnostic report.

## Remaining SMILES surface

### Resolved boundary

There are two relevant Rust paths:

- `umol_io::smiles::parse_smiles_to_ast` parses and raises to AST but does not
  resolve the result;
- `umol_graph::parse::parse_smiles` parses, raises, and runs `Resolver` with
  `SmilesIoConfig::basic_opensmiles()` and `ChemistryModel::default()`.

Python exposes the resolved operation:

```python
mol = MoleculeAst.from_smiles("c1ccccc1")
mol = MoleculeAst.from_smiles(source, config=SmilesIoConfig.opensmiles())
```

Omitting `config` uses the higher-level default. The binding calls the explicit
Rust configured path internally; it does not reproduce parsing or resolution in
Python.

### Configuration scope

This round binds `SmilesIoConfig` with its named presets, including
`basic_opensmiles`, `opensmiles`, `basic_max`, `strict`, `extended`, `lenient`,
and the ChemAxon presets. It also binds `SmilesParseFlags` with bitwise OR so a
caller can construct an arbitrary supported parser-capability combination and
pass it through `SmilesIoConfig.with_parse_flags(...)`. Both values are owned and
immutable on the Python side.

`SmilesLintFlags` and `SmilesLintConfig` are not part of this round. The Rust
parsing workflow does not currently consume them, and lint behavior is not yet
sufficiently specified to expose as an effective Python control. They should be
bound only after the Rust workflow applies a settled lint configuration.

`ChemistryModel` is not part of this round. Binding it would also require public
configuration for valence registries/tables, aromaticity models and scopes, ring
limits, and stereo perception policy. Parser/format configuration and chemistry
resolution policy remain separate axes; adding a model argument later does not
change the ordinary `from_smiles(source, config=None)` path.

Arbitrary lint-name configuration is also deferred unless its ownership is first
made suitable for external callers: `SmilesLintConfig` currently stores
`Vec<&'static str>`, which cannot directly own Python strings. Named IO presets do
not have this problem.

The unresolved `umol-io` AST parser remains an advanced, later surface. This
round's constructor promises a determined molecule or a typed operation error.

### Resolved-SMILES error

The current Rust resolved parser returns `Box<dyn UmolError>` across a fixed
parse → raise → resolve pipeline. Before binding it, the operation receives one
compact concrete error that preserves the categories callers can act on:

- SMILES syntax failure;
- TableIR-to-model conversion failure;
- resolver contradiction;
- underdetermined result;
- resolver execution failure.

The Python mapping uses semantic categories rather than mirroring Rust crate
boundaries. A SMILES syntax failure uses the same public `ParseError` category as
other syntax parsers; contradiction continues to use `ContradictionError`;
`RaiseError` maps to `ModelConversionError`; and underdetermination remains a
separate `UnderdeterminedError`. A fixed operation does not expose a dynamically
boxed error merely because it crosses crates.

### Resolution atomicity

Resolved parsing builds a fresh local `MoleculeAst`, so a failed parse operation
cannot expose the resolver's partially mutated value. The general Rust resolver
must nevertheless have the stronger operation contract: `Resolver::resolve`
retains all accepted narrowing on `Determined` or `Underdetermined`, and leaves
its input unchanged on contradiction or execution failure.

The resolver implements this contract through the existing `Edit`/`Undo`
transaction path rather than whole-molecule clone-and-publish. Each resolver
stage first computes its complete result read-only and emits an edit batch. The
batch is then applied immediately, because the following stage must observe the
fully materialized intermediate `MoleculeAst`:

```text
state0 --valence edits--> state1
state1 --aromaticity edits--> state2
state2 --stereo edits--> state3
state3 --remaining default edits--> final state
```

These are separate transaction batches under one logical resolver transaction.
Their undo journals are appended in chronological order. A later contradiction
or execution failure rolls the combined journal back in reverse order, so later
stages are restored before earlier stages. If applying the current batch fails,
`transact` first restores that batch and the resolver then restores the earlier
combined journal. A rollback failure preserves both the original resolver cause
and the rollback cause. Success retains the materialized state. Public individual
stage resolvers use the same planners with a single atomic batch.

The required shared transaction addition is therefore only chronological
composition of owned `Transaction` journals. It is not an explicit borrowed
transaction scope, arbitrary historical rollback, or savepoint facility.

### Entity update values

Read-only resolver planning needs a neutral way to describe attribute updates
before they are wrapped as host-relative `Edit`s. The existing reaction-update
DSL already has that role, but currently represents an update by placing
`Undetermined` values in the corresponding complete entity AST. That overloads
`Undetermined`: in a complete AST it is a lattice value, while in a partial DSL
field it means "not specified." It also makes an explicit update from a
determined field to `Undetermined` impossible to represent.

Replace the overloaded representation with dedicated update values for all
eight entity families: `AtomUpdate`, `BondUpdate`, `DativeBondUpdate`,
`AromaticSystemUpdate`, `MulticenterBondUpdate`, `NoncovalentBondUpdate`,
`StereoAtomUpdate`, and `StereoBondUpdate`. This is the atom, localized-bond,
and DAMNSS relation family in full; the migration must not leave some entity
update DSL types on the old representation.

An update uses one `Option` per independently parsed ordinary leaf. `None` means
"leave unchanged"; `Some(value)`, including `Some(Undetermined)`, means "set
this leaf to exactly this value." Flat entity fields carry that `Option`
directly. Composite fields preserve their independently addressable AST leaves
rather than replacing the complete composite value. The fixed-field
representation is preferable to a `Vec` of field variants because it excludes
duplicate updates and gives deterministic field order without another keyed
container.

Spin is the shared composite case. `SpinStateAst` contains independently parsed
unpaired-electron and multiplicity `ValueAst` leaves, so atom, localized-bond,
aromatic-system, and multicenter-bond updates reuse:

```rust
pub struct SpinStateUpdate {
    pub unpaired: Option<ValueAst>,
    pub multiplicity: Option<ValueAst>,
}
```

An omitted `#u` or `#s` preserves that component; an explicit `#u*` or `#s*`
sets only that component to `Undetermined`. Thus applying the update `#s1` to
`#u2#s3` produces `#u2#s1`, not `#u*#s1`. Complete-entity parsing still starts
from a fresh `SpinStateAst` and may subsequently apply configured defaults, but
update parsing does not default a missing component before merging it with the
current entity. Edit and delta projection materialize the resulting complete
`SpinStateAst` and continue to emit one whole-spin field change with exact
`old` and `new` values.

Each update reuses its entity's existing constraint container. An empty
container means that no constraint key changes. A determined constraint sets or
replaces its key; an undetermined constraint removes that key. Removal and
setting a constraint to undetermined are the same canonical operation, so there
is no separate constraint-removal variant.

Each complete entity AST exposes both directions of the public update algebra:

```rust
let update: AtomUpdate = atom.difference_to(&other);
assert!(atom.update(&update).canonical_eq(&other));
```

`difference_to` is total because the outer field `Option` distinguishes an
omitted field from an explicit undetermined target. The name is intentionally
shared with `MoleculeAst::difference_to`: both methods return the natural
directed difference at their structural level. An entity has only attribute
differences and returns its `*Update`; a molecule also has structural changes
and returns `Deltas` under a correspondence.

The update is also the canonical input for constructing attribute operations in
both mutation vocabularies. The constructors require the current entity so they
can populate exact `old` values and omit canonical no-ops:

```rust
let edits: Vec<Edit> =
    Edit::for_atom_update(AtomHandle::Id(id), current, &update);
let deltas: Vec<AtomDelta> = AtomDelta::for_update(id, current, &update);
```

Repeat `Edit::for_*_update` and entity-`Delta::for_update` for all eight entity
families. Existing full-state entity delta diffing becomes composition through
the update value:

```rust
let update = lhs.difference_to(rhs);
let deltas = AtomDelta::for_update(id, lhs, &update);
```

Structural delta additions and removals do not pass through an update. Explicit
stereo `Apply`, `Swap`, and `Mirror` deltas also remain direct operations because
an absolute target state does not retain which relative operation produced it.
The two projections may share crate-internal traversal code, but that code does
not introduce another stored change representation and neither projection is
implemented by converting through the other vocabulary.

Rename the existing
`PartialAtomDsl`, `PartialBondDsl`, `PartialDativeBondDsl`,
`PartialAromaticSystemDsl`, `PartialMulticenterBondDsl`,
`PartialNoncovalentBondDsl`, `PartialStereoAtomDsl`, and `PartialStereoBondDsl`
wrappers to the corresponding `AtomUpdateDsl`, `BondUpdateDsl`,
`DativeBondUpdateDsl`, `AromaticSystemUpdateDsl`, `MulticenterBondUpdateDsl`,
`NoncovalentBondUpdateDsl`, `StereoAtomUpdateDsl`, and `StereoBondUpdateDsl`.
These types are reaction modification payloads; no present operation gives
"partial entity" a meaning distinct from "entity update." The renamed wrappers
parse and render the dedicated update values directly. `Edit` and the
corresponding entity `Delta` independently project the same update into their
own operation vocabularies, filling `old` values from the current entity. There
is no generic `AttributeChange` layer and no partial `Edit`-to-`Delta`
conversion.

The current resolver edits are direct projections of their computed candidates:

- The counts and atom-typing valence planners inspect all atoms before mutation.
  For every selected candidate they emit `ModifyAtomField` changes for each
  resolver-owned field that changes, including implicit hydrogens, lone pairs,
  and spin, and keyed `ModifyAtomConstraint` changes for valence,
  aromatic-valence, and other derived constraints. Each edit records the exact
  current value as `old`; unchanged fields and constraint keys emit nothing.
- The aromaticity planner identifies systems from materialized valence results
  and emits `AddAromaticSystem` with the final system AST. Charge-delocalization
  policy may additionally emit atom charge and aromatic-valence changes. It also
  emits the localized-bond aromatic constraint changes currently performed by
  aromatic-system insertion. Resetting the source aromatic-valence constraints
  is opt-in and emits keyed constraint removals.
- The stereo planner identifies sites from materialized aromaticity results,
  computes ligand frames, and emits `AddStereoAtom` and `AddStereoBond`. Resetting
  the source tetrahedral or cis-trans constraint is opt-in and emits the
  corresponding keyed atom- or bond-constraint removal.
- The localized-bond and multicenter default stages follow the same field-diff
  pattern. Their charge and spin field-change variants and transactional
  application/undo support already exist.

A consumed keyed constraint is removed with `old: Some(constraint), new: None`.
An explicit undetermined constraint is vacuous and canonicalizes to absence, so
the planners do not emit it as an intermediate spelling of removal.

Candidate selection may later retain several states rather than requiring one
selected state. That changes the candidate value or `Underdetermined` outcome,
not the edit-diff and transaction protocol.

This round does not add `transact_validated`, generic validator composition,
pre-transaction hooks, resolver savepoints, or fallback policy. A post-apply
acceptance gate becomes necessary only when a resolver postcondition cannot be
decided during read-only planning. Savepoints become necessary only for a
concrete fallback that must retry from an earlier materialized stage.

## Remaining fingerprint surface

### Rust algorithms and output shapes

The molecular fingerprint facility lives in `umol-graph::fingerprint`:

| Rust entry point | Configuration | Output |
|---|---|---|
| WL | rounds + frozen refinement scheme | `FeatureSet<u64>` / `CountedFeatureSet<u64>` |
| ECFP | radius + seed | `FeatureSet<u64>` / `CountedFeatureSet<u64>` |
| Morgan | radius | `FeatureSet<u64>` / `CountedFeatureSet<u64>` |
| `PatternFingerprinter` | width, default 2048 | `BitFp` |
| `SubstructureFeaturizer` | maximum bond count | `FeatureSet<Vec<u8>>` |

The Python workflow accepts operation-specific fingerprint configs rather than
adding each algorithm parameter directly to `MoleculeAst`. Dispatch from each
config to the underlying algorithm remains explicit.

The hashed operation has two methods with fixed return types and one shared
configuration:

```python
mol.hashed_fingerprint(config) -> HashedFeatureSet
mol.counted_hashed_fingerprint(config) -> CountedHashedFeatureSet
```

`HashedFingerprintConfig` selects WL, ECFP, or Morgan and carries the selected
algorithm's parameters. Binary versus counted output does not duplicate that
configuration: it is an aggregation/result choice over the same featurization.
The config argument is required because these algorithms are sibling fingerprint
definitions, not optional capabilities layered over one minimal default. Each
variant may still provide defaults for its own radius, rounds, or named scheme.

The initial config exposes frozen named hash schemes, not raw hash seeds or
aggregation controls. A named scheme defines the reproducible fingerprint
identity and carries its identifier width; the resulting hashed feature set
reports that width as metadata. WL and ECFP use algorithm-specific named scheme
types because their recipes are not interchangeable. Morgan's scheme is fixed by
the pinned RDKit-compatible definition, so its ordinary config carries no scheme
field. Raw custom schemes may be added later if a concrete interoperability need
outweighs the loss of a named reproducibility contract.

The initial public scheme identities are:

```python
WlHashScheme.Xxh3Sorted64V1()
EcfpHashScheme.Xxh3_64V1()
```

These names identify the frozen recipe and version used by the initial binding.
Choosing descriptive names here does not establish a general rule against stable
code names for future schemes; it only avoids promoting the current placeholder
bird names into public compatibility commitments.

The hashed config variants have the following parameter policy:

- `Morgan(radius=2)` uses the conventional radius-2 default;
- `Ecfp(radius=2, scheme=<frozen default>)` uses the conventional radius-2
  default and the frozen named ECFP scheme;
- `Wl(rounds, scheme=<frozen default>)` requires an explicit
  `RefinementRounds` value because `Fixed(n)` and `ToFixpoint` define materially
  different fingerprints; its named scheme may have a workflow default.

`RefinementRounds` is wrapped directly, preserving its `Fixed(n)` and
`ToFixpoint` variants rather than encoding fixpoint as `None` or a sentinel
integer. Public scheme names are stable and versioned rather than inheriting the
current placeholder family names.

Pattern and structural fingerprints remain separate operations with separate
configs and fixed result types:

```python
mol.pattern_fingerprint(
    config=PatternFingerprintConfig(width=2048),
) -> BitFp
mol.structural_fingerprint(
    StructuralFingerprintConfig(max_bonds),
) -> StructuralFeatureSet
```

`PatternFingerprintConfig` is optional at the method boundary because the
2048-bit pattern fingerprint is an established ordinary baseline; specifying a
config layers width and future pattern-library/matching policy onto the same
operation. Width must be positive.

`StructuralFingerprintConfig` is required because `max_bonds` directly defines
the feature universe and computational bound and has no neutral conventional
default. Zero is valid and produces atom features only. The config can later
carry enumeration-algorithm, label-policy, and related bounds without changing
the method signature.

All implemented molecular families belong to this round, staged as the common
hashed featurizers, pattern fingerprints, then exact substructure features.
Reaction fingerprints form a separate operation, config, and result family;
they are not another output mode of a molecular fingerprint call. The current
Difference/DisjointUnion pipeline is exposed as:

```python
reaction.combined_fingerprint(config)
```

`ReactionCombinedFingerprintConfig` is an enum-like config whose variant selects
the reaction combination and whose payload is the required molecule-side hashed
config:

```python
ReactionCombinedFingerprintConfig.Difference(
    molecule=HashedFingerprintConfig.Morgan(),
)
ReactionCombinedFingerprintConfig.DisjointUnion(
    molecule=HashedFingerprintConfig.Morgan(),
)
```

Neither the variant nor the molecule fingerprint definition has a default. The
workflow config lowers internally to the explicit Rust reaction combinator;
`ReactionCombinator` is not a separate ordinary Python workflow argument.
Difference computes counted molecule features before subtraction, while
DisjointUnion computes binary molecule features before tagging them by reaction
side. The operation returns a dedicated `ReactionCombinedFingerprint`:

```python
ReactionCombinedFingerprint.Difference(
    features=SignedHashedFeatureSet(...),
)
ReactionCombinedFingerprint.DisjointUnion(
    features=RoleTaggedHashedFeatureSet(...),
)
```

`SignedHashedFeatureSet` stores `(identifier, signed_count)` entries;
`RoleTaggedHashedFeatureSet` stores `(reaction_side, identifier)` entries. Both
retain hash-width metadata and remain distinct from molecule fingerprint result
types. DRFP and BRIDGIT are deferred, genuinely separate reaction operations with
future `DrfpConfig`/`BridgitConfig` values; the current surface must leave room
for them without treating them as variants of `HashedFingerprintConfig`.

### Python fingerprint values

Fingerprint results remain dedicated immutable Python values rather than being
reduced immediately to tuples, mappings, or bytes:

- `HashedFeatureSet` for sorted sparse integer identifiers;
- `CountedHashedFeatureSet` for sorted `(identifier, count)` entries;
- `BitFp` for fixed-width folded bits;
- `StructuralFeatureSet` for exact `Vec<u8>` canonical keys.

The hashed wrappers preserve the Rust identifier width rather than hardcoding
the current `u64` producers. Internally each is a tagged specialization over the
supported `u32`, `u64`, and `u128` `FeatureSet`/`CountedFeatureSet` forms; current
featurizers exercise the `u64` form. Python identifiers remain ordinary `int`,
while an `id_width` (or equivalently named dtype metadata) property retains the
hash-width contract. Equality, similarity, and subset operations require
compatible widths. This width metadata is part of the wrapper design even though
native NumPy export is deferred.

`StructuralFeatureSet` remains a distinct class rather than another hashed-width
variant. Its keys export as detached Python `bytes`, the native immutable
representation for an opaque byte vector.

The wrappers retain the Rust operations:

- `HashedFeatureSet`: ids, length/iteration, Tanimoto, Dice, subset, and fold;
- `CountedHashedFeatureSet`: entries, length/iteration, and count lookup;
- `BitFp`: width, indexed access, population count, Tanimoto, Dice, and subset;
- structural feature sets: keys, length/iteration, and subset.

Native exports are detached snapshots. The binding does not require NumPy and
does not invent a second implementation of similarity or folding.

Native NumPy export is deferred. The ordinary exports already permit
`numpy.asarray(features.ids, dtype=numpy.uint64)` when NumPy is present, while
the Rust-backed Tanimoto, Dice, subset, count, and fold operations avoid an array
round-trip. A later optional integration can add `__array__` or `to_numpy` for
hashed identifiers using the retained width metadata and a deliberately specified
packed/unpacked representation for `BitFp`; `u128` requires a nonstandard NumPy
representation, and variable-length structural keys remain better represented
as `list[bytes]` than an object-dtype array.

### Fingerprint operation safety

Fingerprint generation has one shared molecular precondition: the molecule must
be ground. Every public concrete and enum-dispatched Rust entry point must check
that precondition consistently before Python can call it. Concrete featurizers
must not reach `expect("ground atom")` through a public path.

Dynamic fingerprint arguments are ordinary API errors:

- fold width must be nonzero;
- bit lookup must bounds-check;
- operations on two runtime-width bit fingerprints must reject unequal widths or
  define explicit unequal-width semantics.

These conditions must return errors suitable for `ValueError`/`IndexError`, not
panic or assert. A caller-visible ground-molecule witness type is not required.

## Shared Python operation errors

Python exposes a small semantic taxonomy shared across these workflows:

- `ParseError` for textual syntax failures;
- `ModelConversionError` when an accepted source representation cannot be
  represented in the requested molecular model;
- `InvalidStructureError` when an already constructed model value fails an
  operation's structural preconditions;
- `ContradictionError` for a chemistry contradiction;
- `UnderdeterminedError` when an operation requires a determined molecule;
- built-in `RuntimeError` for a failure that should be impossible after accepted
  input;
- built-in `ValueError` and `IndexError` for ordinary dynamic arguments.

These are boundary categories, not a Python copy of every Rust error enum. They
do not require prepared-operation types, validation witnesses, checked/unchecked
method pairs, or application reports.

For fingerprints, a non-ground molecule maps to `UnderdeterminedError`; a proven
inconsistency maps to `ContradictionError`. For resolved SMILES, the concrete
operation error maps syntax to `ParseError`, TableIR-to-`MoleculeAst` raising to
`ModelConversionError`, contradiction to `ContradictionError`,
underdetermination to `UnderdeterminedError`, and an unexpected resolver failure
to `RuntimeError`.

For reaction application, invalid reaction or host structure is reported as
`InvalidStructureError` when the iterator is created, where it can be checked
once. Expected match-local rejection remains filtered. An internal application
failure is raised as `RuntimeError` from `next` and permanently terminates the
one-shot iterator.

`ModelConversionError` describes conversion between molecular representations
or models. Rust-to-Python and Python-to-Rust wrapper conversion is binding
machinery, not model conversion; failures there continue to use the appropriate
ordinary `TypeError`, `ValueError`, or `RuntimeError`.

## Scope taken from doc 148

This round adopts the operation-boundary conclusions of doc 148, not its complete
transaction and validator architecture.

| Doc 148 area | This round | Reason |
|---|---:|---|
| direct public operation shape | yes | preserves the ordinary `from_smiles`, fingerprint, and `apply` workflows |
| concrete resolved-SMILES error | yes | required for a stable Python boundary |
| uniform fingerprint ground check | yes | prevents public paths to unchecked literal access |
| safe fingerprint width/index arguments | yes | caller-controlled values must not panic |
| reaction application error classification | yes | expected rejection must not hide internal failure |
| dedicated entity update values | yes | removes the overloaded partial-AST meaning of `Undetermined` before resolver edit planning |
| resolver edit planning and whole-pipeline rollback | yes | makes the public resolver atomic while preserving materialized stage dependencies |
| owned transaction-journal composition | yes | combines sequential resolver batches under one rollback boundary |
| generic validated transactions | no | broader infrastructure project |
| open validator target/composition design | no | not required to expose these workflows |
| transformer atomicity | no | no transformer is added here |
| resolver fallback or savepoints | no | no concrete alternative-strategy retry requires them yet |
| general product-postcondition framework | no | not required for this binding round |

The narrow reaction classification uses the present `ApplyError` meanings:

- `Dangling` is match-local rejection;
- an embedding-dependent `StructuralConflict` is match-local rejection;
- `Inconsistent` is invalid reaction input and should be detected once;
- `Transaction` is an internal failure and must be surfaced.

This classification can improve the existing lazy iterator without first
implementing `transact_validated`, a molecule read facade, validator combinators,
or resolver savepoints. Broader validation may refine the categories later
without changing the Python happy path.

## Overview implementation plan

The subitem-level plan should preserve green boundaries and order the remaining
work as follows:

1. **Entity updates and atomic resolution.** Replace the overloaded complete-AST
   payloads behind all eight entity update DSL families with dedicated update values;
   add chronological composition of owned transaction journals; refactor resolver
   stages into read-only planning plus edit emission; apply each stage immediately
   so the next observes its materialized result; and roll the complete pipeline
   back on contradiction or execution failure.
2. **Other operation contracts in Rust.** Introduce the concrete resolved-SMILES
   error; make fingerprint precondition and dynamic-argument paths non-panicking;
   classify reaction application outcomes without changing matching algorithms.
3. **Workflow configuration values.** Bind `SmilesIoConfig`; define fingerprint
   and reaction-application configs with high-level defaults and explicit
   algorithm selection internally.
4. **Resolved SMILES.** Add `MoleculeAst.from_smiles`, semantic error mapping,
   preset coverage, and resolved-result conformance tests.
5. **Fingerprint values.** Bind `HashedFeatureSet`,
   `CountedHashedFeatureSet`, `StructuralFeatureSet`, and `BitFp` with their
   operations and detached native exports.
6. **Molecular fingerprint algorithms.** Bind the common hashed algorithms,
   pattern fingerprinting, and finally exact substructure features.
7. **Reaction fingerprints.** Bind the existing Difference/DisjointUnion
   operation through `ReactionCombinedFingerprintConfig` and reaction-specific
   results. Keep DRFP and BRIDGIT as deferred, separate operation/config families.
8. **Reaction refinement.** Change `ReactionAst.apply` to accept the application
   config/default, retain eager matching and lazy derivation emission, filter only
   match rejection, and surface fatal iterator errors.
9. **Public verification.** Register/export every public value, run installed
   Python workflows across SMILES, fingerprints, and reactions, and close with
   the complete Rust/Python suites, workspace clippy, rustfmt, and diff checks.

Benchmarks and parity corpora should be selected before algorithm wrappers are
implemented, so binding overhead and bit/feature identity are measured against
the Rust entry points rather than inferred after the API is complete.

## Staged implementation plan

Every subitem includes focused verification and leaves its crate green unless it
is explicitly marked as a breaking migration. Rust tests use `rstest`, table
cases for variant families, exact result or error assertions, and source-order
placement. Breaking Rust or Python signatures are migrated with all workspace
callers before their stage ends.

### S0 — Measurement, update, and transaction foundations

- **S0a — fingerprint identity and cost baselines**
  (`umol-graph/benches/fingerprint.rs`, existing fingerprint fixtures, and a
  small binding-overhead harness): select the conformance molecules and reaction
  fixtures used throughout the implementation; record exact WL, ECFP, Morgan,
  pattern, structural, Difference, and DisjointUnion results from the current
  Rust entry points; and record Rust computation costs separately from an empty
  Python-call baseline. Reuse the existing SMILES corpus and pinned RDKit
  fixtures rather than introducing a second corpus. Tests assert the exact
  fixture identities; the benchmark records measurements without pass/fail
  thresholds. **Additive (green).** `[dep: —]`
- **S0b — chronological transaction composition**
  (`umol-ast/src/ast/molecule/transact.rs`): add one owned `Transaction`
  composition operation that appends a later undo journal after an earlier one.
  Combined rollback must reverse the complete chronological journal, including
  field, keyed-constraint, and added-overlay edits. Tests use two and three
  materialized batches and assert exact restoration, including failure after a
  later batch. Do not add borrowing scopes, implicit drop rollback, arbitrary
  snapshots, or savepoints. **Additive (green).** `[dep: —]`
- **S0c — `AtomUpdate` migration**
  (`umol-ast/src/ast/spin.rs`, `ast/atom.rs`, `dsl/atom.rs`, `dsl/reaction.rs`,
  delta/edit projection): add shared leaf-wise `SpinStateUpdate`, optional
  ordinary leaves, and `AtomConstraintsAst`; implement
  `SpinStateAst::{update, difference_to}`, `AtomAst::{update, difference_to}`,
  `Edit::for_atom_update`, and
  `AtomDelta::for_update`; route full-state `AtomDelta::diff` through the update;
  rename `PartialAtomDsl` to `AtomUpdateDsl`; and migrate atom reaction modify
  parsing/rendering. Tests cover omitted, determined, and explicit undetermined
  leaves, including independent `#u`/`#s` preservation; empty/set/replace/remove
  constraints; deterministic edit and delta projection; the update-difference
  law; and DSL round trips.
  **Breaking Rust/DSL representation migration (red→green).** `[dep: —]`
- **S0d — `BondUpdate` migration**
  (`umol-ast/src/ast/bond.rs`, `dsl/bond.rs`, `dsl/reaction.rs`, delta/edit
  projection): repeat the S0c contract for order, charge, shared leaf-wise spin,
  and
  `BondConstraintsAst`; add `BondAst::{update, difference_to}`,
  `Edit::for_bond_update`, and `BondDelta::for_update`; route full-state bond
  diffing through the update; rename `PartialBondDsl` to `BondUpdateDsl`; and
  migrate localized-bond reaction modification without changing structural
  add/remove deltas. Table tests cover every field, keyed constraint set/removal,
  difference, projections, and DSL round trips. **Breaking Rust/DSL
  representation migration (red→green).**
  `[dep: S0c]`
- **S0e — `DativeBondUpdate` migration**
  (`umol-ast/src/ast/dative.rs`, `dsl/dative.rs`, `dsl/reaction.rs`, delta/edit
  projection): repeat the update contract for dative-bond fields and constraints,
  add `DativeBondAst::{update, difference_to}`,
  `Edit::for_dative_bond_update`, and `DativeBondDelta::for_update`; route
  full-state diffing through the update; then rename `PartialDativeBondDsl` to
  `DativeBondUpdateDsl` and migrate reaction modification. Tests cover the
  complete field/constraint surface, difference, projections, and DSL round
  trips. **Breaking Rust/DSL representation migration (red→green).** `[dep: S0c]`
- **S0f — `AromaticSystemUpdate` migration**
  (`umol-ast/src/ast/aromatic.rs`, `dsl/aromatic.rs`, `dsl/reaction.rs`,
  delta/edit projection): repeat the update contract for aromatic-system fields,
  shared leaf-wise spin, and constraints; replace the current overloaded
  partial-spin representation with the common `SpinStateUpdate` semantics; add
  `AromaticSystemAst::{update, difference_to}`,
  `Edit::for_aromatic_system_update`, and `AromaticSystemDelta::for_update`;
  route full-state diffing through the update; then rename
  `PartialAromaticSystemDsl` to `AromaticSystemUpdateDsl` and migrate reaction
  modification. Tests cover the complete field/constraint surface, difference,
  projections, and DSL round trips. **Breaking Rust/DSL representation migration
  (red→green).** `[dep: S0c]`
- **S0g — `MulticenterBondUpdate` migration**
  (`umol-ast/src/ast/multicenter.rs`, `dsl/multicenter.rs`, `dsl/reaction.rs`,
  delta/edit projection): repeat the update contract for multicenter-bond fields,
  shared leaf-wise spin, and constraints; replace the current overloaded
  partial-spin representation with the common `SpinStateUpdate` semantics; add
  `MulticenterBondAst::{update, difference_to}`,
  `Edit::for_multicenter_bond_update`, and `MulticenterBondDelta::for_update`;
  route full-state diffing through the update; then rename
  `PartialMulticenterBondDsl` to `MulticenterBondUpdateDsl` and migrate reaction
  modification. Tests cover the complete field/constraint surface, difference,
  projections, and DSL round trips. **Breaking Rust/DSL representation migration
  (red→green).** `[dep: S0c]`
- **S0h — `NoncovalentBondUpdate` migration**
  (`umol-ast/src/ast/noncovalent.rs`, `dsl/noncovalent.rs`, `dsl/reaction.rs`,
  delta/edit projection): repeat the update contract for noncovalent-bond fields
  and constraints; add `NoncovalentBondAst::{update, difference_to}`,
  `Edit::for_noncovalent_bond_update`, and `NoncovalentBondDelta::for_update`;
  route full-state diffing through the update; then rename
  `PartialNoncovalentBondDsl` to `NoncovalentBondUpdateDsl` and migrate reaction
  modification. Tests cover the complete field/constraint surface, difference,
  projections, and DSL round trips. **Breaking Rust/DSL representation migration
  (red→green).** `[dep: S0c]`
- **S0i — `StereoAtomUpdate` migration**
  (`umol-ast/src/ast/stereo.rs`, `dsl/stereo.rs`, `dsl/reaction.rs`, delta/edit
  projection): repeat the update contract for stereo-atom fields and constraints,
  add `StereoAtomAst::{update, difference_to}`, `Edit::for_stereo_atom_update`,
  and `StereoAtomDelta::for_update`; route absolute full-state diffing through
  the update while retaining direct relative delta operations; then rename
  `PartialStereoAtomDsl` to `StereoAtomUpdateDsl` and migrate reaction
  modification. Tests cover absolute and relative configurations, keyed
  constraint set/removal, difference, projections, and DSL round trips.
  **Breaking Rust/DSL representation migration (red→green).** `[dep: S0c]`
- **S0j — `StereoBondUpdate` migration**
  (`umol-ast/src/ast/stereo.rs`, `dsl/stereo.rs`, `dsl/reaction.rs`, delta/edit
  projection): add `StereoBondAst::{update, difference_to}`,
  `Edit::for_stereo_bond_update`, and `StereoBondDelta::for_update`; route
  absolute full-state diffing through the update while retaining direct relative
  delta operations; rename `PartialStereoBondDsl` to `StereoBondUpdateDsl`; and
  migrate its reaction modification paths. Tests cover absolute and relative
  configurations, keyed constraint set/removal, difference, projections, and
  DSL round trips. **Breaking Rust/DSL representation migration (red→green).**
  `[dep: S0c]`
- **S0k — entity-update property suite**
  (`umol-ast/tests/property/update.rs`, `tests/property/strategies.rs`): add
  strategies for all eight `*Update` and `*UpdateDsl` families, including
  omitted leaves, explicit undetermined leaves, independent spin components,
  determined constraints, and undetermined constraint removals. For each family,
  prove the string and EDN
  parse/render round trips over the representable update grammar, and prove
  `lhs.update(&lhs.difference_to(&rhs)).canonical_eq(&rhs)` over generated entity
  AST pairs. Preserve minimized regression cases in the normal proptest
  regression file. **Additive verification gate (green).**
  `[dep: S0c, S0d, S0e, S0f, S0g, S0h, S0i, S0j]`

S0 ends with unchanged public workflows, frozen comparison fixtures, and the
complete eight-family update vocabulary plus the minimal transaction composition
needed by atomic resolution, with the public update algebra and update DSLs
covered by cross-family properties.

### S1 — Edit-planned, atomic resolution

- **S1a — counts-valence edit planner**
  (`umol-graph/src/ops/valence/counts.rs`): separate read-only candidate
  calculation from application. Add the public read-only
  `CountsValence::plan(&MoleculeAst) -> Result<Vec<Edit>, CountsError>`; the edit
  vector is the plan, with no plan wrapper or planner trait. Planning inspects
  every atom before mutation. For each applicable atom it first computes the
  complete selected `AtomAst`, derives an `AtomUpdate` with
  `current.difference_to(&selected)`, and projects that update through
  `Edit::for_atom_update`, thereby recording exact current `old` payloads and
  omitting canonical no-ops. The existing single-atom resolver remains a local
  calculate-then-assign operation.

  The public molecule resolver is the plan-and-apply convenience API. It maps a
  planning contradiction to `Solution::Contradictory`, otherwise creates a
  `MoleculeEditor` from the unchanged source, applies the complete edit vector
  with the existing checked `MoleculeEditor::transact`, and publishes
  `editor.build()` only after that batch succeeds. Its application error remains
  separate from `CountsError`:

  ```rust
  pub fn resolve(
      &self,
      ast: &mut MoleculeAst,
  ) -> Result<Solution<(), CountsError>, TransactionError>;
  ```

  There is no additional `apply`, `rollback`, or "successful plan" helper:
  transaction application and rollback use the existing `MoleculeEditor` and
  `Transaction` APIs directly. Tests assert the complete public plan, successful
  materialization, empty-plan no-op omission, and unchanged input on a
  contradiction discovered after an earlier atom was successfully planned.
  **Public planning API plus atomic resolver migration (green).** `[dep: S0k]`
- **S1b — atom-typing edit planner**
  (`umol-graph/src/ops/valence/atom_typing.rs`): add the public read-only
  `AtomTypingValence::plan(&MoleculeAst) -> Result<Vec<Edit>, AtomTypingError>`
  with the same plan representation and application boundary as counts valence.
  Plan construction first requires a concrete element for every atom, reporting
  the offending `AtomId`. For each eligible atom, clone its current `AtomAst`
  once, extend only that clone with topology-derived
  `AtomView::derive_constraints(false)`, select the preferred compatible
  registry entry according to `compare_valence_preference`, and narrow the clone
  against the borrowed selected entry. Derive an `AtomUpdate` with
  `current.difference_to(&selected)` and project it with
  `Edit::for_atom_update`; the source molecule is never enriched in place.

  The public resolver maps planning failure to `Solution::Contradictory`, applies
  the complete edit vector through one checked transaction, and publishes only
  after successful application:

  ```rust
  pub fn resolve(
      &self,
      ast: &mut MoleculeAst,
  ) -> Result<Solution<(), AtomTypingError>, TransactionError>;
  ```

  Registry candidates remain borrowed and selection uses an iterator rather
  than materializing a candidate vector. The temporary cost is therefore one
  cloned `AtomAst` per eligible atom, not one clone per registry candidate.
  Avoiding that clone would require a merged view over stored and derived
  constraints; defer such an optimization to the constraint-container
  restructuring instead of introducing a temporary matching abstraction here.
  Classification remains distinct: `classify_molecule_atom` uses
  `derive_constraints(true)`, where absent overlays are definite negatives,
  while planning uses the pre-perception `false` form. Tests cover the exact
  plan, derived constraints, selected field values, empty-plan identity,
  concrete-element and no-match errors, successful transactional
  materialization, unchanged input after a late contradiction, and exact
  preservation of unrelated constraints. **Public planning API plus atomic
  resolver migration (green).** `[dep: S0k]`
- **S1c — valence resolver dispatch**
  (`umol-graph/src/ops/resolve/valence.rs`): add the public dispatcher
  `ValenceResolver::plan(&MoleculeAst) -> Result<Vec<Edit>,
  ValenceContradiction>`, mapping the selected counts or atom-typing planner's
  chemistry error into the existing shared contradiction type. The public
  `resolve` method dispatches through `plan` rather than through the engines'
  plan-and-apply convenience methods, applies that edit vector with one checked
  `MoleculeEditor::transact`, and publishes only after success. Planning
  contradiction remains `Solution::Contradictory`; transaction failure remains
  `ValenceError::Transaction`. This makes the valence stage itself the single
  atomic application boundary while retaining the engine-specific public
  convenience APIs. Shared exact-plan, successful-resolution, contradiction,
  and unchanged-input tests cover both model variants and ensure that a failure
  after an earlier plannable atom never exposes a partially resolved prefix.
  **Public dispatch planner plus atomic stage application (green).** `[dep: S1a,
  S1b]`
- **S1d — aromaticity edit planner**
  (`umol-graph/src/ops/aromaticity.rs`, `ops/resolve/aromaticity.rs`, resolver
  policy configuration): add public `AromaticityResolverConfig {
  delocalize_charge, reset_aromatic_valence }`, defaulting to the existing
  behavior (`true`, `false`), and construct either with `new(model)` or
  `with_config(model, config)`. Expose the read-only planner as

  ```rust
  pub fn plan(
      &self,
      ast: &MoleculeAst,
  ) -> Result<Solution<Vec<Edit>, AromaticityContradiction>, AromaticityError>;
  ```

  Perception identifies systems from the materialized valence state. For each
  determined system, derive the final `AromaticSystemAst` before emitting
  `AddAromaticSystem`; homogeneous charge delocalization is a pure calculation
  returning atom updates and the final system charge/electron counts, while
  heterogeneous systems retain localized charges. Merge charge changes and an
  opt-in aromatic-valence reset into one `AtomUpdate` per atom, then project
  them with `Edit::for_atom_update`. A reset uses
  `AromaticValenceAst::Undetermined` in the update and therefore projects to
  `old: Some(constraint), new: None`, never to a stored vacuous constraint.
  Mark every localized bond induced by the system atom set through `BondUpdate`
  and `Edit::for_bond_update`.

  `resolve` dispatches through `plan`, applies a determined edit vector with a
  single checked `MoleculeEditor::transact`, and publishes only after success;
  underdetermination and contradiction do not apply a partial plan. Tests cover
  the exact plan, empty-plan identity, homogeneous delocalization on/off,
  heterogeneous localization, exact bond marking, canonical source-constraint
  removal, and unchanged input on perception contradiction/planning error.
  **Public planner plus atomic resolver migration (green).** `[dep: S0k, S1c]`
- **S1e — stereo edit planner**
  (`umol-graph/src/ops/resolve/stereo.rs`): compute atom/bond sites, canonical
  ligand frames, and final stereo ASTs from the materialized aromaticity state;
  emit `AddStereoAtom`/`AddStereoBond` directly and express opt-in removal of
  consumed tetrahedral/cis-trans constraints through `AtomUpdate`/`BondUpdate`.
  Tests cover both supported stereo kinds, ligand-frame identity, reset on/off,
  no-op sites, and unchanged input on the configured inconsistency path.
  **Additive internal refactor (green).**
  `[dep: S0k, S1d]`
- **S1f — localized-bond default edit planner**
  (`umol-graph/src/ops/resolve/bonds.rs`): replace direct charge/spin defaulting
  with read-only `BondUpdate` construction, `Edit::for_bond_update`, and one
  transaction batch. Tests assert exact old/new fields, preservation of
  determined values, no-op batches, and rollback. **Additive internal refactor
  (green).** `[dep: S0k]`
- **S1g — multicenter default edit planner**
  (`umol-graph/src/ops/resolve/multicenter.rs`): replace direct charge/spin
  defaulting with read-only `MulticenterBondUpdate` construction,
  `Edit::for_multicenter_bond_update`, and one transaction batch. Tests mirror
  S1f over multicenter values and assert complete restoration. **Additive
  internal refactor (green).** `[dep: S0k]`
- **S1h — composite resolver transaction**
  (`umol-graph/src/ops/resolve.rs`): plan and apply valence, aromaticity, stereo,
  localized-bond, and multicenter batches in order. Each batch is materially
  present before the next planner runs; each returned transaction is appended to
  the composite journal. `Determined` and `Underdetermined` retain the final
  state; contradiction or execution failure restores the original input. Extend
  `ResolverError` only as needed to preserve transaction and dual-cause rollback
  failures, migrating exhaustive workspace matches in the same subitem. Tests
  prove stage-to-stage observation, complete reverse rollback from every stage,
  accepted underdetermined narrowing, and both-cause diagnostics. **Breaking
  error-enum migration (red→green).** `[dep: S0b, S1c, S1d, S1e, S1f, S1g]`

S1 is the atomic-resolution milestone. The general resolver is correct without
clone-and-publish, `transact_validated`, validator combinators, or savepoints.

### S2 — Concrete resolved-SMILES contract

- **S2a — resolved-SMILES error value**
  (`umol-graph/src/parse.rs`): add one concrete operation error preserving the
  fixed parse → TableIR conversion → resolve categories: SMILES parse failure,
  `RaiseError`, resolver contradiction, underdetermination, and resolver execution
  failure. Implement source conversions and exact display/source behavior.
  Table tests construct every category directly without depending on the parser
  to manufacture unrelated failures. **Additive (green).** `[dep: S1h]`
- **S2b — resolved parser migration**
  (`umol-graph/src/parse.rs` and workspace callers): change the resolved SMILES
  functions from `Box<dyn UmolError>` to the S2a error, preserving the default
  `SmilesIoConfig::basic_opensmiles()` and explicit configured path. Migrate all
  callers and tests in the same subitem. End-to-end cases distinguish syntax,
  TableIR-to-model conversion, contradiction, underdetermination, and successful
  determined output. **Breaking return-type migration (red→green).** `[dep: S2a]`

### S3 — Safe and reproducible Rust fingerprint contracts

- **S3a — stable named hash schemes**
  (`umol-graph/src/hash.rs`, `fingerprint/wl.rs`, `fingerprint/ecfp.rs`): add the
  frozen WL `Xxh3Sorted64V1` and ECFP `Xxh3_64V1` recipe identities and route the
  ordinary featurizers through them. Exact-ID fixtures pin seed, aggregation,
  width, and version. Placeholder bird constructors may remain internal during
  migration but are not exposed as Python compatibility identities. **Additive
  (green).** `[dep: S0a]`
- **S3b — uniform ground-input contract**
  (`umol-graph/src/fingerprint/{featurizer,wl,ecfp,morgan,pattern,substructure}.rs`):
  make every public concrete and enum-dispatched fingerprint entry point reject
  non-ground molecules before literal extraction. Convert concrete WL, ECFP, and
  Morgan methods to checked results and migrate reaction fingerprinting,
  benchmarks, and workspace callers within the stage. Tests run the same ground,
  non-ground, and contradictory inputs through every concrete and enum path and
  assert `FingerprintError` variants rather than panics. **Breaking signature
  migration (red→green).** `[dep: S3a]`
- **S3c — checked `BitFp` and folding arguments**
  (`umol-graph/src/fingerprint/{bit_fp,feature_set}.rs`): replace panicking bit
  lookup, zero-width folding, and unequal-width similarity/subset paths with
  checked Rust contracts (`Option` for lookup and typed/nonzero or `Result`
  contracts for width-bearing operations). Migrate pattern fingerprinting,
  benches, and tests in the same subitem. Cases cover zero width, last valid and
  first invalid index, equal/unequal widths, empty values, and unchanged valid
  similarity results. **Breaking signature migration (red→green).** `[dep: S0a]`
- **S3d — fingerprint parity gate**
  (`umol-graph` fingerprint tests and benchmark): rerun S0a identities through
  the checked APIs; pin counted/binary agreement, structural byte keys, pattern
  bits, and both reaction combinators; and record any cost change caused by
  checks. No wrapper work begins until this gate is green. **Additive (green).**
  `[dep: S3b, S3c]`

### S4 — Python errors and resolved SMILES

- **S4a — binding dependencies and module skeletons**
  (`umol-py/Cargo.toml`, `src/lib.rs`, `src/smiles.rs`,
  `src/fingerprint.rs`): add optional `umol-graph` and `umol-io` dependencies to
  the existing `graph` feature and create the binding modules without public
  classes. `cargo check` covers default and no-default feature builds. **Additive
  (green).** `[dep: S2b, S3d]`
- **S4b — semantic exception classes**
  (`umol-py/src/error.rs`, `src/lib.rs`, `python/umol/__init__.py`): retain
  `ParseError` and `ContradictionError`; add and export `ModelConversionError`,
  `InvalidStructureError`, and `UnderdeterminedError`; and add mapping helpers
  while reserving built-in `RuntimeError`, `ValueError`, and `IndexError` for the
  settled categories. Rust/PyO3 tests assert exact class mapping and messages;
  installed tests assert imports. **Additive (green).** `[dep: S4a]`
- **S4c — composable `SmilesParseFlags`**
  (`umol-py/src/smiles.rs`): bind every effective parser capability and named
  parse preset, with value equality, repr, and bitwise OR producing another
  immutable flag value. Round-trip tests cover individual bits, representative
  combinations, zero/basic OpenSMILES, and rejection of unknown bits. Do not
  bind lint flags. **Additive (green).** `[dep: S4a]`
- **S4d — `SmilesIoConfig`** (`umol-py/src/smiles.rs`): bind the owned immutable
  config, named presets, `with_parse_flags`, structural equality, and repr.
  Conversion tests cover every public preset and an arbitrary OR-composed flag
  set; no lint or chemistry-model fields appear. **Additive (green).**
  `[dep: S4c]`
- **S4e — resolved `MoleculeAst.from_smiles`**
  (`umol-py/src/molecule.rs`, `src/error.rs`): add
  `from_smiles(source, config=None)` over the explicit configured Rust path;
  omission selects the higher-level default. Map syntax, model conversion,
  contradiction, underdetermination, and unexpected resolver execution to the
  settled Python classes. Installed tests assert exact determined structures,
  preset/custom-flag behavior, each reachable error category, and detached
  ownership. **Additive (green).** `[dep: S2b, S4b, S4d]`

S4 is the first complete Python deliverable: resolved SMILES with effective
parser configuration and typed failures.

### S5 — Python fingerprint configuration values

- **S5a — `RefinementRounds`**
  (`umol-py/src/fingerprint/config.rs`): bind `Fixed(rounds)` and `ToFixpoint()`
  directly, with inherent `from_rust`/`to_rust`, equality, payload access, and
  constructor-shaped repr. Cases cover zero/fixed rounds and fixpoint without a
  sentinel representation. **Additive (green).** `[dep: S4a]`
- **S5b — `WlHashScheme`**
  (`umol-py/src/fingerprint/config.rs`): bind the WL-specific scheme type and
  only its initial `Xxh3Sorted64V1()` identity. Conversion tests pin recipe and
  64-bit metadata; no seed or aggregation constructor is public. **Additive
  (green).** `[dep: S3a, S4a]`
- **S5c — `EcfpHashScheme`**
  (`umol-py/src/fingerprint/config.rs`): bind the ECFP-specific scheme type and
  only its initial `Xxh3_64V1()` identity. Conversion tests pin recipe and
  64-bit metadata; no seed constructor is public. **Additive (green).**
  `[dep: S3a, S4a]`
- **S5d — `HashedFingerprintConfig`**
  (`umol-py/src/fingerprint/config.rs`): bind `Morgan(radius=2)`,
  `Ecfp(radius=2, scheme=<default>)`, and
  `Wl(rounds, scheme=<default>)`. The config is required by computation methods;
  conversion lowers each variant to an explicit Rust featurizer. Tests cover
  defaults, explicit parameters, every lowering path, equality, and repr.
  **Additive (green).** `[dep: S5a, S5b, S5c]`
- **S5e — `PatternFingerprintConfig`**
  (`umol-py/src/fingerprint/config.rs`): bind the optional-method config with
  `width=2048`, positive-width validation, conversion, equality, and repr. Tests
  cover default, custom positive width, and zero/negative rejection. **Additive
  (green).** `[dep: S3c, S4a]`
- **S5f — `StructuralFingerprintConfig`**
  (`umol-py/src/fingerprint/config.rs`): bind the required-method config with
  `max_bonds`, accepting zero. Tests cover zero and positive bounds, conversion,
  equality, and repr. **Additive (green).** `[dep: S4a]`
- **S5g — `ReactionCombinedFingerprintConfig`**
  (`umol-py/src/fingerprint/config.rs`): bind required `Difference(molecule=...)`
  and `DisjointUnion(molecule=...)` variants. Lower both the molecular
  featurizer and explicit Rust `ReactionCombinator`; do not export the combinator
  as a separate workflow argument. Tests cover both variants with every hashed
  config family, equality, and repr. **Additive (green).** `[dep: S5d]`
- **S5h — configuration registration and exports**
  (`umol-py/src/lib.rs`, `python/umol/__init__.py`): register/export every S5
  value and list it in `__all__`. Installed tests construct all variants through
  the public package and assert required versus optional arguments. **Additive
  (green).** `[dep: S5d, S5e, S5f, S5g]`

### S6 — Python fingerprint result values

- **S6a — `HashedFeatureSet`**
  (`umol-py/src/fingerprint/value.rs`): add an immutable return-only value with
  internal `u32`/`u64`/`u128` specializations, ordinary Python-int identifier
  snapshots, `id_width`, length/iteration, Tanimoto, Dice, subset, and fold.
  Operations reject incompatible widths. Rust/PyO3 tests cover all three widths,
  detached exports, exact operations, incompatible widths, equality, and repr.
  **Additive (green).** `[dep: S3c, S4b]`
- **S6b — `CountedHashedFeatureSet`**
  (`umol-py/src/fingerprint/value.rs`): add the parallel three-width return-only
  value with detached `(identifier, count)` entries, `id_width`, length/iteration,
  and count lookup. Tests cover sorting assumptions at the Rust boundary,
  absent/present counts, all widths, detachment, equality, and repr. **Additive
  (green).** `[dep: S6a]`
- **S6c — `BitFp`** (`umol-py/src/fingerprint/value.rs`): wrap the checked Rust
  value with width, indexed access, population count, Tanimoto, Dice, and subset.
  Map invalid indices to `IndexError` and unequal widths to `ValueError`. Tests
  cover empty/nonempty values, boundary indices, width mismatches, exact
  similarities, equality, and repr. **Additive (green).** `[dep: S3c, S4b]`
- **S6d — `StructuralFeatureSet`**
  (`umol-py/src/fingerprint/value.rs`): wrap `FeatureSet<Vec<u8>>` as an immutable
  return-only value exposing detached Python `bytes` keys, length/iteration, and
  subset. Tests cover embedded zero bytes, variable lengths, ordering,
  detachment, subset, equality, and repr. **Additive (green).** `[dep: S4a]`
- **S6e — `ReactionSide`**
  (`umol-py/src/fingerprint/reaction.rs`): bind the Reactant/Product value used
  in role-tagged identifiers, with Rust conversion, equality, and repr. Table
  tests cover both sides. **Additive (green).** `[dep: S4a]`
- **S6f — `SignedHashedFeatureSet`**
  (`umol-py/src/fingerprint/reaction.rs`): bind the three-width return-only value
  over `(identifier, signed_count)`, with `id_width`, entries, length/iteration,
  count lookup, equality, and repr. Tests cover cancellation-free nonzero
  entries, positive/negative counts, ordering, widths, and detached exports.
  **Additive (green).** `[dep: S6b]`
- **S6g — `RoleTaggedHashedFeatureSet`**
  (`umol-py/src/fingerprint/reaction.rs`): bind the three-width return-only value
  over `(reaction_side, identifier)`, with `id_width`, ids, length/iteration,
  equality, and repr. Tests cover both sides, ordering, widths, and detached
  exports. **Additive (green).** `[dep: S6a, S6e]`
- **S6h — `ReactionCombinedFingerprint`**
  (`umol-py/src/fingerprint/reaction.rs`): bind the return-only
  `Difference(features=...)` and `DisjointUnion(features=...)` sum type with
  inherent Rust conversion, equality, payload access, and repr. Tests cover both
  Rust variants and ensure their payload classes cannot be interchanged.
  **Additive (green).** `[dep: S6f, S6g]`
- **S6i — result registration and exports**
  (`umol-py/src/lib.rs`, `python/umol/__init__.py`): register/export every public
  S6 value while keeping result constructors return-only where specified.
  Installed tests assert imports, non-constructibility, iteration, and native
  snapshot types (`int`, `bytes`, tuples). **Additive (green).** `[dep: S6a,
  S6b, S6c, S6d, S6e, S6h]`

Native NumPy integration is deliberately absent from S6; retained width metadata
and ordinary snapshots are the compatibility seam for a later optional layer.

### S7 — Molecular fingerprint operations

- **S7a — hashed and counted methods**
  (`umol-py/src/molecule.rs`, `src/fingerprint`): add
  `hashed_fingerprint(config) -> HashedFeatureSet` and
  `counted_hashed_fingerprint(config) -> CountedHashedFeatureSet` over the same
  required `HashedFingerprintConfig`. Dispatch explicitly to WL, ECFP, or Morgan
  and map non-ground input to `UnderdeterminedError`. Tests compare every result
  to S0a exact Rust identities, cover defaults/explicit parameters, and prove
  returned values are detached. **Additive (green).** `[dep: S3d, S4b, S5d,
  S6a, S6b]`
- **S7b — pattern fingerprint method**
  (`umol-py/src/molecule.rs`): add
  `pattern_fingerprint(config=None) -> BitFp`, using the 2048-bit baseline when
  omitted and explicit configured width otherwise. Tests compare exact S0a bits,
  cover optional/default equivalence and custom width, and map non-ground input
  without panics. **Additive (green).** `[dep: S5e, S6c]`
- **S7c — structural fingerprint method**
  (`umol-py/src/molecule.rs`): add required-config
  `structural_fingerprint(config) -> StructuralFeatureSet`. Tests compare exact
  byte keys, cover `max_bonds=0` atom-only output and a bounded connected case,
  map non-ground input, and prove byte snapshots are detached. **Additive
  (green).** `[dep: S5f, S6d]`
- **S7d — molecular fingerprint workflow gate**
  (`umol-py/tests/test_fingerprint.py`, benchmark harness): exercise all five
  public molecule methods from installed Python, exact configs/results,
  similarities, subsets, folding, invalid arguments, and typed non-ground
  failures. Compare binding overhead with S0a without imposing a timing
  threshold. **Additive (green).** `[dep: S7a, S7b, S7c]`

### S8 — Reaction fingerprints

- **S8a — combined reaction method**
  (`umol-py/src/reaction.rs`, `src/fingerprint/reaction.rs`): add required-config
  `ReactionAst.combined_fingerprint(config) -> ReactionCombinedFingerprint`.
  Difference requests counted molecule features and returns signed product-minus-
  reactant counts; DisjointUnion requests binary features and returns role tags.
  Lower the config to the explicit Rust combinator and map inconsistent or
  non-ground reaction sides through the settled semantic exceptions. Tests cover
  identity and changing reactions, every molecular featurizer family, both
  variants, exact width metadata, and detached results. **Additive (green).**
  `[dep: S5g, S6h, S7a]`
- **S8b — reaction fingerprint workflow gate**
  (`umol-py/tests/test_fingerprint.py`): exercise both combined variants from an
  installed package, assert exact S0a results and payload types, and verify that
  molecular and reaction feature-set classes are not interchangeable. **Additive
  (green).** `[dep: S8a]`

DRFP and BRIDGIT remain separate future operation/config families; they do not
extend `ReactionCombinedFingerprintConfig` in S8.

### S9 — Configured and correctly classified reaction application

- **S9a — Rust application preflight and outcome classification**
  (`umol-ast/src/ast/reaction.rs`, validation modules): validate the reaction and
  host once before matching; separate invalid-structure preconditions,
  match-local `Dangling`/embedding-dependent structural rejection, and internal
  transaction/lowering failures. Replace the broad `filter_map(... .ok())`
  behavior with an internal result-bearing path while preserving eager match
  enumeration. Tests assert preflight-before-enumeration, filtered local
  rejection, and visible internal failure with the original diagnostic.
  Migrate the Rust application return contract and all workspace callers in the
  same subitem. **Breaking Rust application migration (red→green).**
  `[dep: —]`
- **S9b — `SubstructureMatchAlgorithm` binding**
  (`umol-py/src/correspondence.rs`): bind `GraphAndOverlays()` and `Incidence()`
  with inherent `from_rust`/`to_rust`, equality, and repr. Table tests cover both
  variants; keep the already-bound six-way `SubgraphIsomorphismAlgorithm` as the
  backend selector. **Additive (green).** `[dep: S4a]`
- **S9c — `ReactionApplicationConfig`**
  (`umol-py/src/reaction.rs`): add the immutable config with default
  `GraphAndOverlays()` strategy and `Vf2Rdkit()` backend. Convert both fields
  explicitly for the Rust call; tests cover defaults, `Incidence`, all backend
  variants, equality, and repr. **Additive (green).** `[dep: S9b]`
- **S9d — fatal-error-aware lazy iterator**
  (`umol-py/src/reaction.rs`): change the private one-shot iterator so
  `__next__` skips only classified match rejection, raises internal failures as
  `RuntimeError`, permanently terminates after a fatal failure, and continues to
  emit owned derivations lazily from the eager correspondence vector. Tests
  cover zero/multiple successes, rejection between successes, fatal termination,
  repeated exhaustion, stable order, and detached results. **Breaking private
  iterator migration (red→green).** `[dep: S9a]`
- **S9e — configured `ReactionAst.apply` migration**
  (`umol-py/src/reaction.rs`, package callers): replace the required algorithm
  argument with `apply(host, config=None)`. Snapshot reaction and host, run S9a
  preflight at iterator creation, eagerly enumerate using both configured
  strategy and backend, and return S9d. Map preflight failure to
  `InvalidStructureError`; preserve lazy derivation construction. Migrate all
  installed examples/tests in the same subitem. **Breaking Python signature
  migration (red→green).** `[dep: S4b, S9c, S9d]`
- **S9f — application workflow gate**
  (`umol-py/tests/test_reaction.py`): update the complete installed reaction
  workflow for default and explicit configs, both match strategies, representative
  backends, eager correspondence snapshotting, lazy owned results, local
  rejection, invalid-structure creation failure, and fatal iterator error.
  **Additive (green).** `[dep: S9e]`

S9 preserves the existing direct `apply` happy path while correcting error
classification. It adds no `apply_at`, report object, or prepared-reaction type.

### S10 — Public integration and release gate

- **S10a — public surface audit**
  (`umol-py/src/lib.rs`, `python/umol/__init__.py`): audit native registration,
  package imports, `__all__`, constructor visibility, signatures, reprs, and
  exception mapping for every S4–S9 value and method. Installed tests compare the
  exported-name set and verify that no lint, raw-seed, combinator, iterator,
  NumPy, DRFP, or BRIDGIT implementation leaked into this round. **Additive
  (green).** `[dep: S7d, S8b, S9f]`
- **S10b — complete installed workflows**
  (`umol-py/tests`): run three coherent public workflows—configured resolved
  SMILES, molecule/reaction fingerprints, and configured reaction application—
  with exact values, typed errors, source non-mutation, detached results, and
  cross-operation composition. **Additive (green).** `[dep: S10a]`
- **S10c — workspace gate** (workspace): run focused resolver/fingerprint/
  reaction suites, complete crate and installed-Python suites, workspace clippy
  over all targets with warnings denied, rustfmt, benchmark compile/run, and
  `git diff --check`; record final test counts and benchmark environment without
  imposing unstable performance thresholds. **Additive (green).** `[dep: S10b]`

At S10 the required round is complete.

## Critical path and deferrals

The resolution/SMILES path is:

The update foundation begins with `S0c`, then fans out through
`{S0d, S0e, S0f, S0g, S0h, S0i, S0j}` and joins at the S0k property gate.
For resolution, `S0k → {S1a, S1b} → S1c → S1d → S1e`; `S1f`, `S1g`, and `S0b`
proceed alongside those planners and join at
`S1h → S2a → S2b → S4a → S4b → {S4c → S4d} → S4e`.

The fingerprint path is:

`S0a → S3a → S3b`, with `S3c` joining at `S3d → S4a`; the parallel configuration
and result branches `{S5, S6}` join at `S7 → S8`.

The application path is:

`S9a → S9d`, while `S4a → S9b → S9c`; both join with S4b at
`S9e → S9f`. The three paths join at `S10a → S10b → S10c`.

The following are explicitly deferrable and are not on the required critical
path: NumPy integration; `SmilesLintFlags`/`SmilesLintConfig`; Python
`ChemistryModel`; the unresolved `umol-io` AST parser; raw/custom hash schemes;
DRFP and BRIDGIT; generic `transact_validated`; validator combinators;
transaction scopes, savepoints, and resolver fallback; general transformer
atomicity; `apply_at` and diagnostic application reports. The complete S10
workflows do not depend on any of them.
