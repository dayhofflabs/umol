# Python molecule workflows: SMILES, fingerprints, and reactions

Status: **Active design / general implementation plan**

Date: 2026-07-13

Updated: 2026-07-18

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

Configured DSL construction is available through `MoleculeDefaults` and
`ReactionDefaults`. Their ordinary constructors select `Required` for every
configurable field and constraint: DSL-to-AST conversion leaves omitted values
undetermined, while AST-to-DSL conversion preserves every value explicitly.
`ground()` fills ordinary entity fields while leaving omitted constraints
required. `MoleculeAst.parse(text, *, defaults=None)` and
`ReactionAst.parse(text, *, defaults=None)` apply the selected policy during
DSL-to-AST conversion; `None` selects the ordinary no-substitution defaults.
Reaction defaults cover the LHS and
add/remove snapshots for all eight entity families. Partial update payloads are
not defaulted. `zeroed()` is deliberately not exposed on the Python side and is
not used by these construction paths.

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

The Rust API now exposes the two boundaries separately and a compact combined
operation:

- `Smiles::parse*` performs syntax parsing and returns the checked SMILES format
  value backed by TableIR;
- `Interpret::interpret` interprets a borrowed format value under an explicit
  `ChemistryModel` and `ResolveConfig` and returns a determined `MoleculeAst`;
- `umol_graph::ingest::ingest_smiles*` composes both boundaries for ordinary callers,
  using `SmilesIoConfig::opensmiles()`, `ChemistryModel::default()`, and
  `ResolveConfig::default()` in the unconfigured form.

Python exposes the resolved operation:

```python
mol = MoleculeAst.from_smiles("c1ccccc1")
mol = MoleculeAst.from_smiles(
    source,
    io_config=SmilesIoConfig.opensmiles(),
    chemistry_model=ChemistryModel.default(),
    resolve_config=ResolveConfig.default(),
)
```

The three configs are keyword-only; omitting any of them uses its higher-level
default. The binding calls the explicit Rust configured path internally; it does
not reproduce parsing or resolution in Python.

### Configuration scope

The accepted configuration design and its staged implementation are maintained
in [doc 155](155-smiles-io-and-resolve-configuration-2026-07-19.md).
`SmilesIoConfig` is a paired parse/render configuration whose shared
`SmilesSyntaxFlags` are composable. The Python ordinary-SMILES surface hides CX
members without pulling the full CX boundary split into this round. Wildcards
are part of OpenSMILES and do not have a capability flag.

`SmilesLintFlags` and `SmilesLintConfig` are not part of this round because the
lint acceptance, diagnostic, and ownership contracts are not yet sufficiently
specified for an effective Python control. They should be bound after those
contracts are settled.

The same plan binds the complete `ChemistryModel` vocabulary and the separate
operational `ResolveConfig`; neither is reduced to a default-only Python value.

Arbitrary lint-name configuration requires owned strings rather than
`Vec<&'static str>` before it can accept Python-provided names. Named IO presets
do not have this problem.

This round removes the unresolved `umol-io` AST shortcut from the public format
API. Syntax parsing instead produces a `Smiles` format value backed by TableIR;
graph-model construction is the separate, explicit ingestion step. The Python
constructor composes both steps and promises a determined molecule or a typed
operation error.

### SMILES input error

The current Rust resolved parser returns `Box<dyn UmolError>` across a fixed
parse → raise → resolve pipeline. Before binding it, parsing is made strictly
syntax-to-`Smiles`, while the graph-owned `Interpret` trait converts that parsed
format value into a determined `MoleculeAst`. The combined convenience operation
receives one compact `SmilesInputError` preserving the categories callers can
act on:

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
  If any element is non-literal, planning returns
  `Solution::Underdetermined(Vec::new())`: the baseline is molecule-wide and
  emits no edits. No later resolver stage runs. Resolving independently
  plannable atoms may be added later as an edit-carrying underdetermined outcome
  without changing the planner return type.
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
WlHashScheme.Xxh3SortedWidth64V1()
EcfpHashScheme.Xxh3Width64V1()
```

These names identify the frozen recipe and version used by the initial binding.
Choosing descriptive names here does not establish a general rule against stable
code names for future schemes; it only avoids promoting the current placeholder
bird names into public compatibility commitments.

The hashed config variants have the following parameter policy:

- `Morgan(radius=2)` uses the conventional radius-2 default;
- `Ecfp(radius=2, hashing_scheme=<frozen default>)` uses the conventional
  radius-2 default and the frozen named ECFP scheme;
- `Wl(rounds, hashing_scheme=<frozen default>)` requires an explicit
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
    config=StructuralFingerprintConfig(max_bonds=max_bonds),
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
2. **Other operation contracts in Rust.** Introduce the parsed `Smiles` format
   value, graph-owned ingestion trait, and concrete `SmilesInputError`; make
   fingerprint precondition and dynamic-argument paths non-panicking; classify
   reaction application outcomes without changing matching algorithms.
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
  `CountsValence::plan(&MoleculeAst) -> Solution<Vec<Edit>, CountsError>`; the
  edit vector is the plan, with no plan wrapper or planner trait. Planning
  inspects every atom before candidate selection. A non-literal element returns
  `Underdetermined(Vec::new())`; a chemistry contradiction returns
  `Contradictory`, and neither outcome mutates the molecule. For each applicable
  atom in a determined plan it first computes the complete selected `AtomAst`,
  derives an `AtomUpdate` with
  `current.difference_to(&selected)`, and projects that update through
  `Edit::for_atom_update`, thereby recording exact current `old` payloads and
  omitting canonical no-ops. The existing single-atom resolver remains a local
  calculate-then-assign operation.

  The public molecule resolver is the plan-and-apply convenience API. It maps a
  planning outcome directly: underdetermination returns without constructing an
  editor, contradiction remains `Solution::Contradictory`, and only a determined
  plan creates a `MoleculeEditor`, applies the complete edit vector with the
  existing checked `MoleculeEditor::transact`, and publishes `editor.build()`
  after that batch succeeds. Its application error remains separate from
  `CountsError`:

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
  `AtomTypingValence::plan(&MoleculeAst) -> Solution<Vec<Edit>,
  AtomTypingError>` with the same plan representation and application boundary
  as counts valence. Plan construction first requires a concrete element for
  every atom; a non-literal element returns `Underdetermined(Vec::new())`
  without selecting any candidates. For each eligible atom in a determined
  plan, clone its current `AtomAst` once, extend only that clone with topology-derived
  `AtomView::derive_constraints(false)`, select the preferred compatible
  registry entry according to `compare_valence_preference`, and narrow the clone
  against the borrowed selected entry. Derive an `AtomUpdate` with
  `current.difference_to(&selected)` and project it with
  `Edit::for_atom_update`; the source molecule is never enriched in place.

  The public resolver returns an underdetermined or contradictory planning
  outcome without mutation. It applies only a determined edit vector through
  one checked transaction and publishes after successful application:

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
  non-literal-element underdetermination, no-match errors, successful transactional
  materialization, unchanged input after a late contradiction, and exact
  preservation of unrelated constraints. **Public planning API plus atomic
  resolver migration (green).** `[dep: S0k]`
- **S1c — valence resolver dispatch**
  (`umol-graph/src/ops/resolve/valence.rs`): add the public dispatcher
  `ValenceResolver::plan(&MoleculeAst) -> Solution<Vec<Edit>,
  ValenceContradiction>`, mapping the selected counts or atom-typing planner's
  contradiction into the existing shared type while preserving determined and
  underdetermined outcomes. The public `resolve` method dispatches through
  `plan` rather than through the engines' plan-and-apply convenience methods.
  It returns immediately on underdetermination, applies only a determined edit
  vector with one checked `MoleculeEditor::transact`, and publishes only after
  success. Planning contradiction remains `Solution::Contradictory`;
  transaction failure remains `ValenceError::Transaction`. The composite
  resolver likewise stops the chain at valence underdetermination, before
  aromaticity and stereo. Shared exact-plan, partial, successful-resolution,
  contradiction, and unchanged-input tests cover both model variants and ensure
  that neither a later open element nor a later contradiction exposes a
  partially resolved prefix.
  **Public dispatch planner plus atomic stage application (green).** `[dep: S1a,
  S1b]`
- **S1d — aromaticity edit planner**
  (`umol-graph/src/ops/aromaticity.rs`, `ops/resolve/aromaticity.rs`, resolver
  policy configuration): add public `AromaticityResolveConfig {
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
  (`umol-graph/src/ops/resolve/stereo.rs`): add public
  `StereoResolveConfig { reset_stereo_constraints }`, defaulting to `false`,
  and expose

  ```rust
  pub fn plan(
      &self,
      ast: &MoleculeAst,
  ) -> Result<Solution<Vec<Edit>, StereoContradiction>, StereoError>;
  ```

  Planning computes atom/bond sites, canonical ligand frames, and final stereo
  ASTs from the materialized aromaticity state. Realizable assertions emit
  `AddStereoAtom`/`AddStereoBond` directly; opt-in removal of consumed
  tetrahedral/cis-trans constraints is projected through `AtomUpdate` and
  `BondUpdate`. The existing model-level `InconsistencyPolicy` is now enforced:
  `Keep` retains an unrealizable assertion, `Strip` emits its canonical keyed
  removal, and `Error` returns `UnrealizableAtom`/`UnrealizableBond` without
  mutation. `resolve` applies a determined plan in one checked transaction.
  Exact-plan and full-result tests cover both supported kinds and ligand frames,
  reset on/off, no-op sites, all three inconsistency policies, and unchanged
  input on contradiction. **Public planner plus atomic resolver migration
  (green).** `[dep: S0k, S1d]`
- **S1f — localized-bond default edit planner**
  (`umol-graph/src/ops/resolve/bonds.rs`): add public
  `plan(&MoleculeAst) -> Vec<Edit>`. For each localized bond, derive the target
  charge and spin without mutation, express independent spin-component changes
  through `SpinStateAst::difference_to`, project the resulting `BondUpdate`
  through `Edit::for_bond_update`, and apply the complete vector in one checked
  transaction from `resolve`. Completely undetermined spin still defaults to a
  closed shell; partially determined spin retains its known component and uses
  the existing high-spin completion rule. Tests assert exact old/new fields,
  determined-value identity, full resolved output, and complete batch rollback
  after a later stale precondition. **Public planner plus atomic resolver
  migration (green).** `[dep: S0k]`
- **S1g — multicenter default edit planner**
  (`umol-graph/src/ops/resolve/multicenter.rs`): add public
  `plan(&MoleculeAst) -> Vec<Edit>` and mirror S1f over
  `MulticenterBondUpdate`/`Edit::for_multicenter_bond_update`, including the
  independent spin-component semantics and one checked transaction in
  `resolve`. Tests mirror the exact-plan, identity, full-result, and stale-plan
  restoration coverage over multicenter values. **Public planner plus atomic
  resolver migration (green).** `[dep: S0k]`
- **S1h — composite resolver transaction**
  (`umol-graph/src/ops/resolve.rs`): the composite resolver now plans and applies
  valence, aromaticity, stereo, localized-bond, and multicenter batches in
  order. After each successful batch it consumes the editor with `build()` and
  immediately opens `state.edit()` for the next batch. This materially exposes
  the preceding stage to the next read-only planner without cloning a
  materialized editor; the accumulated `Transaction` journal remains valid for
  the equivalent reopened editor state. Each returned transaction is appended
  to that journal.

  A final `Determined` or `Underdetermined` result publishes the final state. A
  contradiction or execution failure reverse-replays the composite journal and
  publishes the restored editor state before returning. `ResolverError` adds
  `RollbackFailed { cause: ResolverRollbackCause, rollback }`; the cause
  distinguishes the original resolver contradiction from the original resolver
  execution error, so rollback failure never erases either diagnostic. Tests
  prove valence-to-aromaticity and aromaticity-to-stereo observation, reverse
  rollback after later aromaticity/stereo contradictions, accepted
  underdetermined narrowing, identity, and both dual-cause diagnostic shapes.
  The `umol-ast` property suite additionally compares a generated multi-edit
  overlay transaction with the same edits applied as individually materialized
  batches, appends every returned journal, and proves that the combined rollback
  reconstructs the original AST. **Composite transactional resolver and
  error-enum migration (green).** `[dep: S0b, S1c, S1d, S1e, S1f, S1g]`

S1 is the atomic-resolution milestone. The general resolver is correct without
clone-and-publish, `transact_validated`, validator combinators, or savepoints.

### S2 — Explicit SMILES representation and concrete ingestion

- **S2a — parsed SMILES value and ingestion contract**
  (`umol-io/src/smiles.rs`, `umol-io/src/smiles/parser.rs`, and
  `umol-graph/src/ingest.rs`): make the boundary between parsing an external
  format and constructing the graph model explicit.

  `umol_io::smiles::Smiles`, with a private `table_ir::Molecule` field, is the
  parsed semantic value of a molecular SMILES representation. It is not the
  source string and does not preserve source spelling. The private wrapper
  establishes that the contained TableIR originated from, or was checked for
  representation by, SMILES; arbitrary TableIR cannot be wrapped unchecked.
  `Smiles::parse`, `parse_bytes`, `parse_with`, and `parse_bytes_with` perform
  syntax parsing only and return `Smiles`. `as_table_ir` and `into_table_ir`
  expose the neutral boundary value read-only or by ownership; there is no
  mutable escape that could invalidate the format invariant. This associated
  API replaces public names such as `parse_smiles_bytes_to_table_ir_with`
  without weakening the distinction between syntax parsing and model
  construction.

  The graph-owned public interpretation trait is:

  ```rust
  pub trait Interpret {
      type Output;
      type Error;

      fn interpret(
          &self,
          model: &ChemistryModel,
          resolve_config: &ResolveConfig,
      ) -> Result<Self::Output, Self::Error>;
  }
  ```

  `Interpret for Smiles` returns `MoleculeAst` and borrows the parsed value so the
  same external representation can be interpreted under more than one
  chemistry model. Its implementation performs TableIR-to-`MoleculeAst`
  conversion and then resolution; those are interpretation details rather than
  parsing semantics. `MoleculeInterpretationError` preserves the four post-parse
  categories: model conversion (`RaiseError`), resolver contradiction,
  underdetermination, and resolver execution failure.

  One flat `SmilesInputError` covers the combined text-to-model convenience
  operation. Its public variants are `Syntax(ParseError)`,
  `ModelConversion(RaiseError)`, `Contradiction(ResolverContradiction)`,
  `Underdetermined(ResolveUnderdetermined)`, and `Execution(ResolverError)`.
  Unprefixed display preserves the underlying diagnostic, and each variant
  exposes its wrapped error directly through `Error::source`; conversions from
  `ParseError` and `MoleculeInterpretationError` remove dynamic error erasure
  without exposing pipeline crate boundaries to Python. Both concrete errors
  implement `UmolError`.

  The reverse path follows the same structural boundary when SMILES output is
  added: `MoleculeAst` is converted, with an explicit output configuration, to
  a checked `Smiles` value, and `Smiles::render` serializes that value to text.
  Rendering does not bypass the format value by returning a string directly
  from `MoleculeAst`. Its configuration and errors remain separate from parse
  configuration and errors. `ReactionSmiles` applies the same private-wrapper
  pattern over reaction TableIR. `Mol` and `Sdf` may later do the same over
  their appropriate TableIR values; they need no common `*Format` suffix or
  marker trait until an operation actually requires format-polymorphic code.

  Rust table tests construct every error category directly and assert exact
  display/source behavior. Wrapper tests prove the syntax-only parse boundary,
  private-format invariant, TableIR access, successful ingestion, and reuse
  under distinct chemistry models. **Additive (green).** `[dep: S1h]`
- **S2b — SMILES ingestion migration**
  (`umol-io/src/smiles.rs`, `umol-graph/src/ingest.rs`, the former
  `umol-graph/src/parse.rs`, and workspace callers): resolved workspace callers
  use either `Smiles::parse*` followed by `Interpret::interpret` or the compact
  `umol_graph::ingest::ingest_smiles*` convenience surface. The latter returns
  `SmilesInputError`, uses `SmilesIoConfig::opensmiles()` plus
  `ChemistryModel::default()` as its unconfigured default, and provides explicit
  text/bytes configured paths. The public unresolved SMILES-to-AST helpers and
  public `parse_smiles*_to_table_ir*` functions are retired; the raw byte parser
  remains crate-private for parser implementation and tests. `parse` therefore
  means syntax-to-format-value rather than syntax-to-graph-model. End-to-end
  cases distinguish syntax, TableIR-to-model conversion, contradiction,
  underdetermination, and successful determined output. Python's later
  `Molecule.from_smiles` composes the same two boundaries without requiring a
  public Python `Smiles` wrapper; a future
  `Molecule.to_smiles` likewise converts through the Rust `Smiles` value before
  rendering. **Breaking API migration (green).** `[dep: S2a]`

#### S2b follow-up — TableIR and format-boundary unification

The wildcard configuration failure exposes a lower-level design problem rather
than a missing validation check. `table_ir::Molecule` predates `MoleculeAst` and
was introduced as an optimized representation of fixed-composition molecules.
`ExtendedMolecule` is a strict representational superset: it replaces `Atom`
with `ExtendedAtom` and `Bond` with `ExtendedBond`, and adds optional CTFile and
CX annotation payloads. The basic/extended distinction then propagates into two
SMILES parsers, two MOL parsers, two SDF parsers, capability presets that select
between them, conversion helpers, tests, and benchmarks.

That optimization boundary no longer matches the semantic boundaries of the
system. OpenSMILES includes the `*` wildcard, and `MoleculeAst` represents it
natively as `ElementAst::Undetermined`; parsing it need not fail merely because
the result is not ground. A workflow that promises a determined molecule may
still report underdetermination after conversion and resolution, while a plain
format-to-model conversion can return the patterned `MoleculeAst`. Groundness
belongs at the model or workflow boundary, not in selection of the TableIR
storage type.

The proposed direction is therefore to collapse `Molecule` and
`ExtendedMolecule` into one TableIR `Molecule`, using the current extended
representation as the semantic superset. `ExtendedAtom` and `ExtendedBond`
likewise become the sole `Atom` and `Bond` representations. Each external
format then has one parser and one TableIR result type. Parser configuration may
control acceptance policy, strictness, vendor deviations, ignored data, and
diagnostics, but it must not select a different result structure. The current
`basic_*`/`extended_*` parser families and `BASIC_MAX`/`EXTENDED_MAX`
representation split can disappear.

The extended atom and bond structures do have a larger inline footprint:
`ExtendedAtom` carries query fields, optional vectors, and a property map, while
`ExtendedBond` carries query fields and a property map. Empty collections avoid
their backing allocations but still increase the containing value's size. That
is a storage-layout question, not a reason to retain two semantic types and two
parser implementations. If measurement shows that ordinary concrete molecules
are materially affected, the single TableIR type can later hide cold fields
behind an optional boxed extension record or side tables. Such an optimization
would preserve the unified API and parser.

The existing `umol-io/benches` suites provide the migration gate. The SMILES
benchmark runs the same chain, branch, bond, ring, component, bracket, and
whitespace corpora through the basic and extended parser paths, with a separate
wildcard corpus for the extended path. It can therefore measure the direct cost
of changing the ordinary result representation while verifying that wildcard
support remains effective. The MOL benchmark compares the basic and extended
atom, bond, and property parsers on overlapping inputs, isolating where the
larger representation enters parsing cost. It is mostly a token/parser-component
benchmark rather than a whole-MOL/SDF workload, so unification should retain
those comparisons and add representative end-to-end records before removing
the old path. The decision rule is not that the unified representation must be
free: it is that measured cost must justify a storage optimization inside the
single type, rather than preservation of duplicate public semantics.

TableIR unification does not require every format to share a boundary wrapper.
The wrappers express grammar provenance and rendering guarantees, while the
unified TableIR expresses what was parsed:

- `Smiles` accepts the complete OpenSMILES grammar, including `*`, and wraps the
  unified `table_ir::Molecule`. `WILDCARDS` and `basic_opensmiles` cease to be
  optional capabilities; they are part of the format.
- A future `CxSmiles` may separately accept OpenSMILES plus CXSMILES annotations
  and wrap the same TableIR type. This is useful because CXSMILES has its own
  parsing and rendering invariant and carries useful additions such as radical
  annotations. `Smiles -> CxSmiles` is lossless; `CxSmiles -> Smiles` is
  fallible when CX-only information would be discarded.
- Other graph-string systems such as HELM, DeepSMILES, BigSMILES, or SELFIES are
  independent formats if they are supported at all; they should not be folded
  into an open-ended SMILES dialect flag set.
- `Mol` and `Sdf` likewise wrap the unified TableIR representation. Query
  features, SGroups, RGroups, and vendor extensions affect accepted syntax and
  retained payload, not whether the parser returns `Molecule` or
  `ExtendedMolecule`.

`EXTENDED_AROMATICS` and `EXTENDED_BONDS` are not CXSMILES-specific. Implementations
such as RDKit support varying subsets of this wider SMILES syntax, including
bond forms such as `$` and `~`. They may remain options of the ordinary SMILES
boundary; CX tag parsing and rendering are the part that belongs specifically
to a future `CxSmiles` boundary and configuration.

The next decision is benchmark-gated rather than an unconditional unification:

```text
Measure the ordinary-input cost of ExtendedMolecule
├── acceptable
│   ├── settle the complete umol-io migration scope
│   │   ├── OpenSMILES plus extended aromatic/bond syntax
│   │   ├── CXSMILES
│   │   └── MOL/SDF
│   └── unify TableIR molecule/atom/bond types and parser result types
└── too high
    ├── add OpenSMILES `*` support to basic Molecule
    ├── audit and remove only Smiles options that truly require ExtendedMolecule
    ├── retain the Molecule/ExtendedMolecule split inside umol-io
    └── schedule a focused representation-design spike
```

In either branch, only ordinary SMILES ingestion propagates into `umol-graph`
and `umol-py` in this implementation round. CXSMILES, MOL, and SDF support stop
at the `umol-io` boundary and are scheduled for the next Python workflow pass.
If the cost is acceptable, they should nevertheless be included far enough in
the Rust-side migration to avoid retaining duplicate TableIR and parser APIs
that would have to be removed immediately afterward.

The cost gate measures both parsing latency and retained representation size.
The SMILES suite already supplies a direct basic/extended comparison over the
same ordinary inputs. The MOL suite supplies component-level atom, bond, and
property comparisons; representative whole-record MOL/SDF cases must be added
or measured before using it to justify removal of the basic path. No arbitrary
percentage defines “acceptable” before seeing the data: report the relative
cost across small, ordinary drug-like, and larger structures, together with
atom/bond size and allocation differences. If a material cost appears, first
consider optimizing cold storage inside one semantic type; preserve duplicate
types and parsers only if the evidence justifies that larger design cost.

The fallback audit is feature-specific. It must not remove
`EXTENDED_AROMATICS` or `EXTENDED_BONDS` merely because they currently travel
through an extended parser. Features that the basic representation already can
hold remain available; only syntax whose result genuinely requires
`ExtendedMolecule` is excluded from the `Smiles` boundary until the design
spike resolves it.

This gate precedes S4c-S4e. Python must not freeze the current capability and
configuration split before the Rust representation decision is made.

##### Benchmark gate result

The current `ExtendedMolecule` representation is not an acceptable direct
replacement for `Molecule`. Measurements on 2026-07-18 used the existing
Criterion suites on the same machine and release profile. Representative quick
runs gave the following parser costs; the 10-carbon chain was also repeated as
a standard 100-sample Criterion run and confirmed the directional result.

| Input | Basic path | Extended path | Extended cost |
|---|---:|---:|---:|
| SMILES linear `C10` | 396 ns | 576 ns | +46% |
| SMILES cyclohexane | 372 ns | 474 ns | +27% |
| SMILES mixed-bond `C10` | 460 ns | 609 ns | +32% |
| SMILES adamantane | 485 ns | 771 ns | +59% |
| SMILES 50 bracket atoms | 1.29 us | 2.15 us | +66% |
| MOL 69-column atom record | 86.5 ns | 97.6 ns | +13% |
| MOL 21-column bond record | 26.3 ns | 34.5 ns | +31% |
| MOL single charge property | 49.4 ns | 52.0 ns | +5% |

The inline representation difference is also substantial:

| Type | Basic | Extended | Ratio |
|---|---:|---:|---:|
| atom | 104 bytes | 288 bytes | 2.77x |
| bond | 40 bytes | 96 bytes | 2.40x |
| molecule header | 200 bytes | 440 bytes | 2.20x |

These sizes do not include backing allocations, but they directly increase the
atom and bond vector allocations and retained cache footprint for ordinary
molecules. The consistent SMILES latency increase shows that this is not only a
dormant-data size concern. The component-level MOL results vary with how much
work is representation-independent, reinforcing the need for whole-record MOL
and SDF benchmarks before a later representation redesign.

The immediate path therefore follows the high-cost branch:

1. Add OpenSMILES `*` to basic `Molecule` and raise it to
   `ElementAst::Undetermined`; do not route ordinary SMILES through
   `ExtendedMolecule` merely to represent one core token.
2. Audit `EXTENDED_AROMATICS` and `EXTENDED_BONDS` against the actual basic atom
   and bond structures. Retain every supported form whose value already fits;
   parser provenance alone is not evidence that a feature requires
   `ExtendedMolecule`.
3. Remove CXSMILES-specific flags and presets from the ordinary `Smiles`
   configuration. Keep the extended parser and `ExtendedMolecule` inside
   `umol-io` for CXSMILES, MOL, and SDF while their boundary design is deferred
   to the next pass.
4. Use one complete OpenSMILES default, including `*`. The distinction between
   ordinary and extended aromatic/bond
   syntax remains an acceptance-policy option within SMILES rather than a
   result-type choice.
5. Schedule a representation-design spike for eventual TableIR unification.
   Its target is a compact semantic superset, potentially with cold extension
   records or side tables, measured against these baselines. It must not begin
   by renaming the current extended structures and accepting their cost.
6. Revisit `ChiralityFrame` placement as part of the CXSMILES/MOL/SDF boundary
   and representation spike. The immediate parser invariant keeps the
   molecule-level field `None` unless an `Atom.chirality` descriptor is present;
   the spike should determine whether the frame belongs directly with that raw
   descriptor so the representation cannot express a descriptor without its
   required source frame. Expand this section after the focused work in 152 is
   complete.

The focused implementation plan for items 1 and 4 is
[152-basic-molecule-wildcards-2026-07-18.md](152-basic-molecule-wildcards-2026-07-18.md).
It includes the compact basic-atom representation migration, basic parser and
raise changes, removal of the wildcard configuration split, migration of both
the `umol-graph` classification tools and the `umol-io` parsing-conformance
suite, and the required unit, property, fuzz, and benchmark gates.

##### Follow-up status after 152

The focused wildcard/OpenSMILES work is complete. `Molecule` accepts `*`,
raises it to `ElementAst::Undetermined`, and uses one default OpenSMILES
configuration without `SmilesSyntaxFlags::WILDCARDS`, `BASIC_OPENSMILES`, or
`SmilesIoConfig::basic_opensmiles()`. The former basic/OpenSMILES
classification split has also been removed; the only remaining
`basic_opensmiles` Rust reference is a conformance-table test that asserts the
retired category name is unknown. `EXTENDED_AROMATICS` and `EXTENDED_BONDS`
remain ordinary SMILES acceptance-policy flags because their retained values
fit the basic `Atom`/`Bond` TableIR structures.

The broader S2b follow-up is not closed. `CHEMAXON_EXTENSIONS`,
`SKIP_UNKNOWN_CHEMAXON_TAGS`, `SmilesSyntaxFlags::CHEMAXON`, and
`SmilesIoConfig::chemaxon()` still live on the ordinary `SmilesIoConfig`
surface, and the diagnostic/conformance tools still distinguish
`basic_chemaxon` from `chemaxon` inside the deferred CXSMILES boundary. The
compact TableIR semantic-superset spike has not been executed, and
`ChiralityFrame` remains a molecule-level TableIR field assigned by the SMILES
and CTFile parsers. These items belong to the next CXSMILES/MOL/SDF boundary
and representation pass, not to the completed wildcard task.

CXSMILES, MOL, and SDF do not propagate into `umol-graph` or `umol-py` in this
round. Their current Rust support remains available inside `umol-io`; unified
boundary types and downstream ingestion belong to the next workflow pass.

### S3 — Safe and reproducible Rust fingerprint contracts

- **S3a — stable named hash schemes**
  (`umol-graph/src/hash.rs`, `fingerprint/wl.rs`, `fingerprint/ecfp.rs`): add the
  frozen WL `Xxh3SortedWidth64V1` and ECFP `Xxh3Width64V1` recipe identities and
  route the ordinary featurizers through them. Exact-ID fixtures pin seed,
  aggregation, width, and version. `WlFeaturizer` and `EcfpFeaturizer` store
  those named schemes; raw `RefinementXxh3Scheme::new` remains the low-level
  refinement recipe constructor for Rust-side experiments and hash tests.
  **Implemented (green).** `[dep: S0a]`
- **S3b — uniform ground-input contract**
  (`umol-graph/src/fingerprint/{featurizer,wl,ecfp,morgan,pattern,substructure}.rs`):
  every public concrete and enum-dispatched molecule fingerprint entry point
  rejects non-ground molecules before literal extraction. Concrete WL, ECFP, and
  Morgan methods now return checked results; enum dispatch delegates to those
  concrete checks; pattern and substructure keep their existing checked shape.
  Reaction fingerprinting propagates molecule-side `NotGround` and maps
  inconsistent deltas to `Inconsistent`. Tests pin exact ground output,
  non-ground errors for concrete and enum paths, and reaction inconsistency.
  **Implemented (green).** `[dep: S3a]`

  Wildcard-aware circular fingerprints are a separate research direction rather
  than an exception to this contract. Doc
  [154](154-lattice-aware-probabilistic-fingerprints-2026-07-18.md) centers that
  work on wildcard-bearing reaction representations for enzyme–reaction
  contrastive learning, with exact lattice encodings, frequency-weighted
  containment, and probabilistic refinement as supporting constructions.
- **S3c — checked `BitFp` and folding arguments**
  (`umol-graph/src/fingerprint/{bit_fp,featurizer,pattern}.rs`): `BitFp::get`
  returns `Option<bool>`; folding over `FeatureSet<u32/u64/u128>` rejects zero
  width with `FingerprintError::ZeroWidth`; and
  `BitFp::{tanimoto,dice,is_subset}` reject unequal widths with
  `FingerprintError::WidthMismatch { left, right }`.
  Pattern fingerprinting propagates the checked fold result. Tests cover zero
  width, the last valid and first invalid indices, equal and unequal widths,
  empty values, collisions, and unchanged valid similarity and subset results.
  **Implemented (green).** `[dep: S0a]`
- **S3d — fingerprint parity gate**
  (`umol-graph` fingerprint tests and benchmark): the S0a WL, ECFP, Morgan,
  pattern, structural, Difference, and DisjointUnion identities remain exact
  through the checked APIs. The counted and binary WL, ECFP, and Morgan paths
  produce the same sorted identifier sets. Relative to the saved S0a Criterion
  fixture baseline, ECFP is the only reported regression: +2.78% with a
  95%-confidence interval of +2.25% to +3.28%. Criterion reports WL, structural,
  Difference, and DisjointUnion within its noise threshold and no detectable
  change for Morgan or pattern. The gate has no timing pass/fail threshold.
  **Implemented (green).** `[dep: S3b, S3c]`

### S4 — Python errors and resolved SMILES

- **S4a — binding dependencies and module skeletons**
  (`umol-py/Cargo.toml`, `src/lib.rs`, `src/smiles.rs`,
  `src/fingerprint.rs`): optional direct `umol-graph` and `umol-io` dependencies
  are enabled by the existing `graph` feature. Private `smiles` and `fingerprint`
  binding modules provide the workflow implementation locations without adding
  Python classes. Default and no-default feature builds are green. **Implemented
  (green).** `[dep: S2b, S3d]`
- **S4b — semantic exception classes**
  (`umol-py/src/error.rs`, `src/lib.rs`, `python/umol/__init__.py`): `ParseError`
  and `ContradictionError` remain public; `ModelConversionError`,
  `InvalidStructureError`, and `UnderdeterminedError` are exported alongside
  them. The concrete `SmilesInputError` returned by
  `umol_graph::ingest::ingest_smiles` maps syntax, model conversion,
  contradiction, and underdetermination to those semantic classes, and maps
  resolver execution failure to built-in `RuntimeError`. Rust/PyO3 tests assert
  exact classes and messages; installed tests assert import identity.
  **Implemented (green).** `[dep: S4a]`
- **S4c — composable `SmilesSyntaxFlags`**
  (`umol-py/src/smiles.rs`): the immutable value binds ordinary-SMILES
  `EXTENDED_AROMATICS` and `EXTENDED_BONDS`, plus the `OPENSMILES` and `LENIENT`
  presets. Validated bit construction, the read-only `bits` property, value
  equality, repr, and bitwise OR preserve the `umol-io` flag semantics; unknown,
  retired, and CX-only bits raise `ValueError`. Rust/Python conversion tests
  cover every public capability, representative combinations, OpenSMILES, and
  unknown bits. CX and lint flags are not exposed. **Implemented (green).**
  `[dep: S4a]`
- **S4d — `SmilesIoConfig`** (`umol-py/src/smiles.rs`): the owned immutable
  config exposes `opensmiles()`, `lenient()`, `with_syntax_flags(...)`, a
  detached read-only `syntax_flags` value, structural equality, and repr.
  Lowering reconstructs the `umol-io` config with internal lint defaults; CX,
  lint, and chemistry-model fields are absent from Python. Conversion tests
  cover every public preset and arbitrary OR-composed flags. **Implemented
  (green).** `[dep: S4c]`

  The complete chemistry-model and resolve-operation configuration bindings and
  their supporting Rust migrations are implemented in
  [doc 155](155-smiles-io-and-resolve-configuration-2026-07-19.md).
- **S4e — resolved `MoleculeAst.from_smiles`**
  (`umol-py/src/molecule.rs`, `src/error.rs`):
  `from_smiles(source, *, io_config=None, chemistry_model=None,
  resolve_config=None)` lowers owned configuration values into the fully
  configured Rust ingestion path, with omission selecting the three high-level
  defaults. It returns an exact determined `MoleculeAst` and preserves the
  syntax, model-conversion, contradiction, underdetermination, and execution
  error taxonomy. Tests cover ordinary syntax configurations, both valence
  models, representative aromaticity and stereo models, every resolve-policy
  field, keyword-only rejection, and detached ownership. **Implemented
  (green).** `[dep: S4b, doc 155 S7a]`

S4 completes the configured resolved-SMILES Python deliverable.

### S5 — Python fingerprint configuration values

- **S5a — `RefinementRounds`**
  (`umol-py/src/fingerprint/config.rs`): bind `Fixed(rounds)` and `ToFixpoint()`
  directly, with inherent `from_rust`/`to_rust`, equality, payload access, and
  constructor-shaped repr. Cases cover zero/fixed rounds and fixpoint without a
  sentinel representation. **Implemented (green).** `[dep: S4a]`
- **S5b — `WlHashScheme`**
  (`umol-py/src/fingerprint/config.rs`): bind the WL-specific scheme type and
  only its initial `Xxh3SortedWidth64V1()` identity. Conversion tests pin recipe
  and 64-bit metadata; no seed or aggregation constructor is public.
  **Implemented (green).** `[dep: S3a, S4a]`
- **S5c — `EcfpHashScheme`**
  (`umol-py/src/fingerprint/config.rs`): bind the ECFP-specific scheme type and
  only its initial `Xxh3Width64V1()` identity. Conversion tests pin recipe and
  64-bit metadata; no seed constructor is public. **Implemented (green).**
  `[dep: S3a, S4a]`
- **S5d — `HashedFingerprintConfig`**
  (`umol-py/src/fingerprint/config.rs`): bind `Morgan(radius=2)`,
  `Ecfp(radius=2, hashing_scheme=<default>)`, and
  `Wl(rounds, hashing_scheme=<default>)`. The config is required by computation
  methods; conversion lowers each variant to an explicit Rust featurizer. Tests
  cover defaults, explicit parameters, every lowering path, equality, and repr.
  **Implemented (green).** `[dep: S5a, S5b, S5c]`
- **S5e — `PatternFingerprintConfig`**
  (`umol-py/src/fingerprint/config.rs`): bind the optional-method config with
  `width=2048`, positive-width validation, conversion, equality, and repr. Tests
  cover default, custom positive width, and zero/negative rejection.
  **Implemented (green).** `[dep: S3c, S4a]`
- **S5f — `StructuralFingerprintConfig`**
  (`umol-py/src/fingerprint/config.rs`): bind the required-method config with
  `max_bonds`, accepting zero. Tests cover zero and positive bounds, conversion,
  equality, and repr. **Implemented (green).** `[dep: S4a]`
- **S5g — `ReactionCombinedFingerprintConfig`**
  (`umol-py/src/fingerprint/config.rs`): bind required `Difference(molecule=...)`
  and `DisjointUnion(molecule=...)` variants. Lower both the molecular
  featurizer and explicit Rust `ReactionCombinator`; do not export the combinator
  as a separate workflow argument. Tests cover both variants with every hashed
  config family, equality, and repr. **Implemented (green).** `[dep: S5d]`
- **S5h — configuration registration and exports**
  (`umol-py/src/lib.rs`, `python/umol/__init__.py`): register/export every S5
  value and list it in `__all__`. Installed tests construct all variants through
  the public package and assert required versus optional arguments. **Implemented
  (green).** `[dep: S5d, S5e, S5f, S5g]`

### S6 — Python fingerprint result values

- **S6a — `HashedFeatureSet`**
  (`umol-py/src/fingerprint/value.rs`): add an immutable return-only value with
  internal `u32`/`u64`/`u128` specializations, ordinary Python-int identifier
  snapshots, `id_width`, length/iteration, Tanimoto, Dice, and subset. Fold is
  wired in S6c alongside its `BitFp` result type.
  Operations reject incompatible widths. Rust/PyO3 tests cover all three widths,
  detached exports, exact operations, incompatible widths, equality, and repr.
  **Implemented (green).** `[dep: S3c, S4b]`
- **S6b — `CountedHashedFeatureSet`**
  (`umol-py/src/fingerprint/value.rs`): add the parallel three-width return-only
  value with detached `(identifier, count)` entries, `id_width`, length/iteration,
  and count lookup. Tests cover sorting assumptions at the Rust boundary,
  absent/present counts, all widths, detachment, equality, and repr.
  **Implemented (green).** `[dep: S6a]`
- **S6c — `BitFp`** (`umol-py/src/fingerprint/value.rs`): wrap the checked Rust
  value with width, indexed access, population count, Tanimoto, Dice, and subset.
  Generalize checked Rust folding across the supported identifier widths and
  expose `HashedFeatureSet.fold(width) -> BitFp` here.
  Map invalid indices to `IndexError` and unequal widths to `ValueError`. Tests
  cover empty/nonempty values, boundary indices, width mismatches, exact
  similarities, equality, and repr. **Implemented (green).** `[dep: S3c, S4b]`
- **S6d — `StructuralFeatureSet`**
  (`umol-py/src/fingerprint/value.rs`): wrap `FeatureSet<Vec<u8>>` as an immutable
  return-only value exposing detached Python `bytes` keys, length/iteration, and
  subset. Tests cover embedded zero bytes, variable lengths, ordering,
  detachment, subset, equality, and repr. **Implemented (green).** `[dep: S4a]`
- **S6e — `ReactionSide`**
  (`umol-py/src/fingerprint/reaction.rs`): bind the Reactant/Product value used
  in role-tagged identifiers, with Rust conversion, equality, and repr. Table
  tests cover both sides. **Implemented (green).** `[dep: S4a]`
- **S6f — `SignedHashedFeatureSet`**
  (`umol-py/src/fingerprint/reaction.rs`): bind the three-width return-only value
  over `(identifier, signed_count)`, with `id_width`, entries, length/iteration,
  count lookup, equality, and repr. Tests cover cancellation-free nonzero
  entries, positive/negative counts, ordering, widths, and detached exports.
  **Implemented (green).** `[dep: S6b]`
- **S6g — `RoleTaggedHashedFeatureSet`**
  (`umol-py/src/fingerprint/reaction.rs`): bind the three-width return-only value
  over `(reaction_side, identifier)`, with `id_width`, ids, length/iteration,
  equality, and repr. Tests cover both sides, ordering, widths, and detached
  exports. **Implemented (green).** `[dep: S6a, S6e]`
- **S6h — `ReactionCombinedFingerprint`**
  (`umol-py/src/fingerprint/reaction.rs`): bind the return-only
  `Difference(features=...)` and `DisjointUnion(features=...)` sum type with
  inherent Rust conversion, equality, payload access, and repr. Tests cover both
  Rust variants and ensure their payload classes cannot be interchanged.
  **Implemented (green).** `[dep: S6f, S6g]`
- **S6i — result registration and exports**
  (`umol-py/src/lib.rs`, `python/umol/__init__.py`): register/export every public
  S6 value while keeping result constructors return-only where specified.
  Installed tests assert imports and non-constructibility. PyO3 tests assert
  iteration and native snapshot types (`int`, `bytes`, tuples); the installed
  workflow checks repeat those assertions once S7 and S8 provide public result
  producers. **Implemented (green).** `[dep: S6a, S6b, S6c, S6d, S6e, S6h]`

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
  returned values are detached. **Implemented (green).** `[dep: S3d, S4b,
  S5d, S6a, S6b]`
- **S7b — pattern fingerprint method**
  (`umol-py/src/molecule.rs`): add
  `pattern_fingerprint(config=None) -> BitFp`, using the 2048-bit baseline when
  omitted and explicit configured width otherwise. Tests compare exact S0a bits,
  cover optional/default equivalence and custom width, and map non-ground input
  without panics. **Implemented (green).** `[dep: S5e, S6c]`
- **S7c — structural fingerprint method**
  (`umol-py/src/molecule.rs`): add required-config
  `structural_fingerprint(config) -> StructuralFeatureSet`. Tests compare exact
  byte keys, cover `max_bonds=0` atom-only output and a bounded connected case,
  map non-ground input, and prove byte snapshots are detached. **Implemented
  (green).** `[dep: S5f, S6d]`
- **S7d — molecular fingerprint workflow gate**
  (`umol-py/tests/test_fingerprint.py`, benchmark harness): exercise all four
  public molecule methods across the five WL, ECFP, Morgan, pattern, and
  structural definitions from installed Python, including the counted path and
  exact configs/results, similarities, subsets, folding, invalid arguments, and
  typed non-ground failures. Compare binding overhead with S0a without imposing
  a timing threshold. **Implemented (green).** `[dep: S7a, S7b, S7c]`

  The release binding harness on macOS 15.7.3 arm64 with CPython 3.13.14 records
  the best of five repeats: 24 ns for an empty Python call, 28 ns for the trivial
  PyO3 call, 2.462 µs for WL, 1.939 µs for ECFP, 1.753 µs for Morgan, 1.788 µs
  for counted Morgan, 4.549 µs for pattern, and 9.343 µs for structural. The
  same-machine 100-sample Criterion fixture intervals are 2.235–2.241 µs,
  1.760–1.767 µs, 1.540–1.546 µs, 4.363–4.383 µs, and 8.986–9.015 µs
  respectively for the five binary definitions. The directional Python
  increment is therefore about 0.17–0.36 µs per fixture call, including result
  wrapping; these measurements are recorded without a pass/fail threshold.

### S8 — Reaction fingerprints

- **S8a — combined reaction method**
  (`umol-py/src/reaction.rs`, `src/fingerprint/reaction.rs`): add required-config
  `ReactionAst.combined_fingerprint(config) -> ReactionCombinedFingerprint`.
  Difference requests counted molecule features and returns signed product-minus-
  reactant counts; DisjointUnion requests binary features and returns role tags.
  Lower the config to the explicit Rust combinator and map inconsistent or
  non-ground reaction sides through the settled semantic exceptions. Tests cover
  identity and changing reactions, every molecular featurizer family, both
  variants, exact width metadata, and detached results. **Implemented (green).**
  `[dep: S5g, S6h, S7a]`
- **S8b — reaction fingerprint workflow gate**
  (`umol-py/tests/test_fingerprint.py`): exercise both combined variants from an
  installed package, assert exact S0a results and payload types, and verify that
  molecular and reaction feature-set classes are not interchangeable.
  **Implemented (green).** `[dep: S8a]`

DRFP and BRIDGIT remain separate future operation/config families; they do not
extend `ReactionCombinedFingerprintConfig` in S8.

### S9 — Configured reaction application and algorithm transparency

#### Algorithm-transparency audit

The repository-wide audit found several production calls that select a
`umol-graph-core` algorithm inside an operation whose public arguments or
configuration do not own that choice. These are silent defaults and must not be
preserved merely because an algorithm enum currently has one implementation.
Low-level algorithmic operations pass selectors explicitly; higher-level
workflows may define defaults, but the selected algorithms remain fields of the
workflow configuration and are passed explicitly at the lower boundary.

| Location | Hidden selection | Required correction |
| --- | --- | --- |
| `umol-ast/src/ast/ring.rs` | `CycleEnumerationAlgorithm::Vismara` behind both `rings` and `rings_with` | Separate ring-set semantics into `RingModel` and the family-specific simple/relevant selectors into `RingConfig`; higher-level ring-consuming workflows own inspectable defaults. This also affects ECFP, Morgan, and aromaticity transitively. |
| `umol-ast/src/ast/symmetry.rs` | `AutomorphismAlgorithm::Nauty` for initial symmetry, fixpoint reruns, and site stabilizers | Add the selection to `GraphSymmetryConfig` and retain it in `GraphSymmetry`; `StereoValidateConfig`, not `StereoModel`, owns the higher-level default. |
| `umol-ast/src/ast/compose.rs` | `CommonSubgraphEnumerationAlgorithm::Backtracking` behind `ReactionAst::compose` | Remove `CompositionScope`; complete composition enumerates every admissible overlap. Require the enumeration algorithm at the Rust boundary and expose it as a visible Python keyword with a high-level default. |
| `umol-ast/src/ast/reaction.rs`, `umol-py/src/reaction.rs` | `SubstructureMatchAlgorithm::GraphAndOverlays` while only the graph-core subisomorphism backend is supplied | Rust application accepts both the structure strategy and subisomorphism backend explicitly. Python owns both in `ReactionApplicationConfig`. |
| `umol-graph/src/fingerprint/pattern.rs` | `GraphAndOverlays` plus `SubgraphIsomorphismAlgorithm::Vf2` | `PatternFingerprinter` / `PatternFingerprintConfig` own both selectors. |
| `umol-graph/src/fingerprint/substructure.rs` | `SubgraphEnumerationAlgorithm::Esu` plus `AutomorphismAlgorithm::Nauty` | `SubstructureFeaturizer` / `StructuralFingerprintConfig` own the enumeration and canonicalization selectors. |
| `umol-graph/src/ops/aromaticity/clar.rs` | `MaximumIndependentSetAlgorithm::BranchAndBound` | `AromaticityConfig`, separate from `AromaticityModel`, owns the independent-set selector. |
| `umol-graph/src/ops/aromaticity/hmo.rs` | `ConnectedComponentsAlgorithm::Bfs` | `AromaticityConfig`, separate from `AromaticityModel`, owns the connected-components selector. |
| `umol-graph/src/ops/transform/kekulizer.rs` | `MaximumMatchingAlgorithm::Edmonds` replaced the configured selection for mobile exposure | The behavioral defect is fixed by doc 158 S5d: both modes honor the operation-level selector and surface non-bipartite Hopcroft--Karp input without fallback. Rename `KekulizationModel` to `KekulizationConfig` in S9t. |

`umol-io/src` contains no graph-core algorithm selection. WL, ECFP, and Morgan
constructing their corresponding refinement variants is not a silent default:
the named higher-level operation already specifies that algorithm. Dispatch
matches, tests, benchmarks, and fuzz targets likewise do not constitute hidden
defaults.

Inside `umol-graph-core`, a named algorithm may specify a subsidiary algorithm
as part of its implementation without adding another selector to the public
API. The current cases are Hopcroft-Karp using BFS bipartition, Vismara using
Tarjan biconnected components, and EC circular refinement using BFS
neighborhood traversal. Their implementation functions remain private: do not
widen them to `pub(crate)` or add visibility solely to bypass the ordinary API.
The implementation instead calls the existing public selector-bearing method
with the chosen enum variant and documents why that subsidiary choice belongs
to the named algorithm.

Algorithm selectors and execution bounds are operational configuration, not
chemistry-model parameters. `AromaticityConfig` owns `RingConfig` plus the
connected-components and maximum-independent-set selectors passed to the shared
perception object; sharing is justified here because perception requires ring
enumeration, connected components, and maximum independent set.
`StereoModel` retains semantic perception choices, while
`StereoValidateConfig` owns graph-automorphism selection and the fixpoint
iteration bound. `InconsistencyPolicy` moves to `StereoResolveConfig`, because
it controls resolver behavior. The parallel naming is deliberate: an `XModel`
describes chemical semantics, while an `XConfig` describes execution of the
operation using that model.

The same pass removes three graph-core naming and organization inconsistencies:
the original `MaxMatchingAlgorithm` rename is superseded by doc 158 S5d's
separate `BipartiteMaximumMatchingAlgorithm` and
`GeneralMaximumMatchingAlgorithm`,
`MaxIndependentSetAlgorithm` becomes `MaximumIndependentSetAlgorithm`, and
`planar_matching_count.rs` becomes `matching_count.rs`. Counting remains
separate from matching construction and enumeration: it has its own planar
embedding, Pfaffian, and error machinery, while the public planar type and
method names remain specific. `NonBipartiteGraphError` makes Hopcroft--Karp's
bipartite precondition a release-mode result rather than a debug assertion.

`CompositionScope` does not describe two algorithms for the same operation.
`Full` is complete reaction composition, whereas `RcAnchored` discards valid
composites using a derived reaction-center filter. The implementation performs
complete overlap enumeration before filtering, so `RcAnchored` provides no
search saving; its derived center also omits stereo and important deletion-only
cases. Remove the enum and define `ReactionAst::compose` as complete composition
only, including the empty overlap. With scope removed, composition has one
algorithm selector and no coherent multi-field configuration space. Python
therefore exposes that selector directly as a keyword with a visible workflow
default rather than introducing a one-field composition config.

- **S9a — Rust application preflight and outcome classification**
  (`umol-ast/src/ast/reaction.rs`, validation modules): validate the reaction and
  host once before matching; separate invalid-structure preconditions,
  match-local `Dangling`/embedding-dependent structural rejection, and internal
  transaction/lowering failures. Replace the broad `filter_map(... .ok())`
  behavior with an internal result-bearing path while preserving eager match
  enumeration. Tests assert preflight-before-enumeration, filtered local
  rejection, and visible internal failure with the original diagnostic.
  Migrate the Rust application return contract and all workspace callers in the
  same subitem. `ReactionAst::apply` now returns preflight failure separately
  from its result-bearing lazy application iterator; the iterator filters only
  classified match-local rejection and terminates after yielding one internal
  failure. Attribute deltas are re-projected through the corresponding
  `*Update` against each matched host entity, so pattern-side undetermined
  values do not become incorrect edit preconditions and independently omitted
  spin leaves remain unchanged. **Implemented (green).**
  The application property suite now covers all eight entity-update families
  with effective pattern differences, host refinement, partial spin, keyed
  constraints, and canonical-equivalence assertions; the reaction properties
  pass a 4,096-case soak.
  `[dep: —]`
- **S9b — `SubstructureMatchAlgorithm` binding**
  (`umol-py/src/correspondence.rs`): bind `GraphAndOverlays()` and `Incidence()`
  with inherent `from_rust`/`to_rust`, equality, and repr. Table tests cover both
  variants; keep the already-bound six-way `SubgraphIsomorphismAlgorithm` as the
  backend selector. **Implemented (green).** `[dep: S4a]`
- **S9c — `ReactionApplicationConfig`**
  (`umol-ast/src/ast/reaction.rs`, `umol-py/src/reaction.rs`): migrate the
  low-level Rust `ReactionAst::apply` operation to accept both
  `SubstructureMatchAlgorithm` and `SubgraphIsomorphismAlgorithm` explicitly,
  and migrate all Rust callers in the same subitem. Add the immutable Python
  config with default `GraphAndOverlays()` strategy and `Vf2Rdkit()` backend.
  Convert both fields explicitly for the Rust call; tests cover defaults,
  `Incidence`, all backend variants, equality, and repr. **Implemented
  (green).** The Rust application signature and all callers now supply both
  selectors; direct application tests cover both structure strategies. Python
  exposes the immutable, keyword-only `ReactionApplicationConfig`, including
  explicit Rust conversions and package-level export, while the Python
  `ReactionAst.apply` signature remains scheduled for S9e. `[dep: S9b]`
- **S9d — fatal-error-aware lazy iterator**
  (`umol-py/src/reaction.rs`): change the private one-shot iterator so
  `__next__` skips only classified match rejection, raises internal failures as
  `RuntimeError`, permanently terminates after a fatal failure, and continues to
  emit owned derivations lazily from the eager correspondence vector. Tests
  cover zero/multiple successes, rejection between successes, fatal termination,
  repeated exhaustion, stable order, and detached results. **Implemented
  (green).** The iterator now filters only `ApplyError::is_match_rejection()`;
  the first other application error becomes a `RuntimeError` and marks the
  iterator terminal before it is returned. `[dep: S9a]`
- **S9e — configured `ReactionAst.apply` migration**
  (`umol-py/src/reaction.rs`, package callers): replace the required algorithm
  argument with `apply(host, config=None)`. Snapshot reaction and host, run S9a
  preflight at iterator creation, eagerly enumerate using both configured
  strategy and backend, and return S9d. Map preflight failure to
  `InvalidStructureError`; preserve lazy derivation construction. Migrate all
  installed examples/tests in the same subitem. **Implemented (green).**
  `ReactionAst.apply(host, *, config=None)` snapshots reaction and host,
  validates both before matching, eagerly captures correspondences using the
  configured selectors, and returns the lazy S9d iterator. Precondition errors
  map to `InvalidStructureError`; installed callers use the default or explicit
  keyword-only config. `[dep: S4b, S9c, S9d]`
- **S9f — application workflow gate**
  (`umol-py/tests/test_reaction.py`): update the complete installed reaction
  workflow for default and explicit configs, both match strategies, representative
  backends, eager correspondence snapshotting, lazy owned results, local
  rejection, invalid-structure creation failure, and fatal iterator error.
  **Implemented (green).** The installed workflow covers the default
  `GraphAndOverlays`/`Vf2Rdkit` pair, explicit `GraphAndOverlays`/`Ullmann` and  `Incidence`/`Vf2` pairs, snapshot and detachment behavior, a rejected match
  between successful derivations, precondition failure at iterator creation,
  and a fatal error raised lazily with permanent exhaustion. `[dep: S9e]`
- **S9g — graph-core terminology and matching-count module**
  (`umol-graph-core/src/algorithms/{matching,mis,matching_count}.rs`,
  `algorithms.rs`, `lib.rs`, all workspace callers): rename
  `MaxMatchingAlgorithm` to `MaximumMatchingAlgorithm` and
  `MaxIndependentSetAlgorithm` to `MaximumIndependentSetAlgorithm`. Rename
  `planar_matching_count.rs` to `matching_count.rs` without moving its contents
  into the already-large matching-construction module; retain
  `PlanarEmbedding`, `PlanarEmbeddingError`, `PlanarMatchingCountError`, and
  `count_perfect_matchings_planar`. Migrate imports, type annotations, tests,
  benchmarks, and documentation in the same subitem. **Breaking public renames
  and module move (red→green).** **Implemented (green).** Graph-core and every
  workspace caller adopted the then-current `MaximumMatchingAlgorithm` and
  `MaximumIndependentSetAlgorithm`; matching-count implementation and exports
  live under `algorithms::matching_count`, while its planar-specific public
  types and operation retain their existing names. Doc 158 S5d subsequently
  replaced the graph-core matching selector with separate bipartite and general
  selectors. `[dep: —]`
- **S9h — fallible maximum matching and genuine Hopcroft-Karp**
  (`umol-graph-core/src/algorithms/matching.rs`,
  `umol-ast/src/ast/{matching,view/graph}.rs`, all callers): add
  `MaximumMatchingError::NonBipartite` and make `Graph::maximum_matching` and
  the typed `GraphView` adapter return `Result`. Replace the current repeated
  augmenting-BFS implementation behind `HopcroftKarp` with layered BFS plus
  batched DFS, preserving caller-supplied node-order determinism. Validate the
  bipartite precondition in all build modes; never substitute Edmonds. Record
  the pre-change matching benchmark baseline, retain the benchmark after the
  replacement, and migrate every caller in the same subitem. Tests cover the
  typed odd-cycle error, Edmonds/Hopcroft-Karp cardinality parity on bipartite
  graphs, deterministic order, and propagation through `GraphView`.
  **Breaking return-contract migration plus algorithm correction (red→green).**
  **Implemented (green).** `Graph::maximum_matching` and `GraphView` now return
  the exported typed error; all callers handle it explicitly, and the
  kekulizer preserves aromatic-system context. Hopcroft-Karp now uses shortest
  augmenting-path layers with batched DFS, with randomized bipartite parity
  against Edmonds. A retained Criterion benchmark records the local
  pre-change → post-change ranges as 343.54–345.70 ns → 243.91–244.58 ns for
  hexagon, 450.81–453.62 ns → 270.20–273.83 ns for cubane, and
  12.693–12.820 µs → 4.2859–4.3043 µs for a 16×16 grid. `[dep: S9g]`
  Doc 158 S5d subsequently replaced this interim API with the direct
  `bipartite_maximum_matching` and `general_maximum_matching` operations plus
  the infallible `bipartite_maximum_matching_or_general` fallback operation,
  and replaced `MaximumMatchingError` with `NonBipartiteGraphError` on the
  direct Hopcroft--Karp path.
- **S9i — graph-core subsidiary algorithms**
  (`umol-graph-core/src/algorithms/{matching,cycles,refine}.rs`): use the
  existing public selector-bearing calls for BFS bipartition inside
  Hopcroft-Karp, Tarjan biconnected components inside Vismara, and BFS
  neighborhood traversal inside EC circular refinement. Document at each call
  why the subsidiary choice is fixed by the named parent algorithm. Do not
  expose private implementation methods or widen them to `pub(crate)`. Focused
  algorithm tests and strict graph-core Clippy remain green. **Additive
  documentation (green).** **Implemented (green).** All three implementations
  already called the public selector-bearing APIs; their call sites now state
  why the subsidiary algorithm is fixed by the parent operation. No private
  visibility was widened and no behavior changed. `[dep: S9h]`
- **S9j — Python graph-core selector values**
  (`umol-py/src/algorithm.rs`, `lib.rs`, `python/umol/__init__.py`): establish
  the singular binding module for algorithm selectors, moving the existing
  `SubstructureMatchAlgorithm` and `SubgraphIsomorphismAlgorithm` wrappers there
  without changing their public names. Bind
  `SimpleCycleEnumerationAlgorithm`,
  `RelevantCycleEnumerationAlgorithm`, `AutomorphismAlgorithm`,
  `CommonSubgraphEnumerationAlgorithm`, `SubgraphEnumerationAlgorithm`,
  `MaximumIndependentSetAlgorithm`, and `ConnectedComponentsAlgorithm` with
  inherent `from_rust`/`to_rust`, equality, and repr. Table tests cover every
  variant and installed tests cover all exports. **Additive public values plus
  internal module move, followed by the family-specific cycle-selector
  migration (green).** **Implemented (green).** The singular binding module
  owns all nine selector wrappers; reaction application imports the moved
  selectors from it without a Python API change. The family-blind cycle
  selector was replaced by the separate Read--Tarjan simple-cycle and Vismara
  relevant-cycle values specified by doc 158; all selectors are registered and
  exported with variant-complete Rust and installed-Python coverage.
  `[dep: S9b, S9g]`
- **S9k — family-specific ring model and configuration**
  (`umol-ast/src/ast/{molecule,ring,view/graph,view/ring}.rs`, all Rust
  callers): replace the intermediate family-blind selector with
  `RingModel { kind, max_ring_size }` and
  `RingConfig { simple_cycle_algorithm, relevant_cycle_algorithm }`.
  Rename `RingFamily` to `RingSetKind`; make `MoleculeAst::rings(model, config)`
  the sole general entry point; and make `RingSet::enumerate` dispatch to the
  family-specific graph-core collector. Remove `rings_with`, `atom_filter`,
  induced-cycle filtering, and endpoint-based bond reconstruction. Migrate all
  Rust callers together. Exact and property tests cover both ring-set kinds,
  bounds, selector routing, edge identity, reindexing, and the fixed Relevant
  projection used by existing ring constraints. **Breaking AST ring API and
  complete caller migration (red→green).** **Implemented (green).**
  `RingModel` and `RingConfig` are public AST values with defaults of
  Relevant/22 and Read--Tarjan/Vismara. `RingSet::enumerate` takes the graph,
  model, and config in that order, consumes edge-aware graph-core `Cycle`
  values, and excludes one- and two-atom cycles only at the chemical-ring
  boundary. The general molecule surface now consists only of
  `MoleculeAst::rings`; every workspace caller supplies `RingModel` and
  `RingConfig` explicitly. Existing ring-membership constraints retain a
  documented fixed Relevant projection. The complete design and validation
  record is in [doc 158](158-ring-model-and-enumeration-2026-07-22.md).
  `[dep: —]`

- **S9l — ring selection in hashed-fingerprint configs**
  (`umol-graph/src/fingerprint/{ecfp,morgan,featurizer}.rs`,
  `umol-py/src/{fingerprint/config,ring}.rs`): add `ring_config` to
  `EcfpFeaturizer` and `MorganFeaturizer`; their high-level constructors retain
  the AST default as an inspectable workflow choice. Bind the keyword-only
  `RingConfig` and extend the ECFP and Morgan variants of
  `HashedFingerprintConfig`, conversions, equality, and repr. Migrate reaction
  fingerprint configs through their nested molecule config. Rust and installed
  Python tests compare default and explicit selectors and preserve exact
  fingerprint identities. **Breaking config-shape migration (red→green).**
  **Implemented (green).** ECFP and Morgan now store and pass the complete
  `RingConfig` into the fixed Relevant/22 fingerprint `RingModel`; WL remains
  unchanged. Python exports a frozen keyword-only `RingConfig` with independent
  simple- and relevant-cycle selectors. Both hashed configurations lower the
  nested value to their concrete Rust featurizers, and reaction difference and
  disjoint-union fingerprints preserve it through their nested molecular
  configuration. Exact molecular and reaction payload tests pass under default
  and explicit configuration. `[dep: S9j, S9k]`
- **S9m — aromaticity model/config separation**
  (`umol-graph/src/ops/aromaticity.rs`, `aromaticity/{clar,hmo}.rs`,
  `resolve/aromaticity.rs`, `transform/aromatizer.rs`,
  `validate/aromaticity.rs`, `umol-py/src/resolve.rs`): add
  `AromaticityConfig { ring_config,
  connected_components_algorithm, maximum_independent_set_algorithm }`, with
  Read--Tarjan/Vismara, BFS, and branch-and-bound as inspectable high-level
  defaults.
  Pass the shared config explicitly to `AromaticityPerception::find_systems`;
  each top-level operation owns its copy rather than adding selectors to the
  model or changing the perception dispatch type. Add the shared config as
  `AromaticityResolveConfig.perception`; add default and configured constructors
  to the aromatizer and aromaticity validator. Bind `AromaticityConfig` in
  Python and extend `AromaticityResolveConfig`, leaving `AromaticityModel`
  unchanged. Migrate all Rust and Python callers together. Tests prove
  default-result parity and selector propagation through resolution,
  aromatization, validation, HMO, and Clar.
  **Additive config followed by breaking resolve-config migration (red→green).**
  **Implemented (green).** Rust `AromaticityConfig` owns `RingConfig`, the
  connected-components selector, and the maximum-independent-set selector,
  with Read--Tarjan/Vismara, BFS, and branch-and-bound defaults.
  `AromaticityPerception` constructs the fixed Relevant ring model at the
  chemistry model's size bound and threads each operational selector to the
  corresponding low-level operation. Resolver configuration nests it as
  `perception`; the aromatizer and validator have default and configured
  constructors. Python exports the frozen keyword-only configuration and nests
  it in `AromaticityResolveConfig`. Exact Hückel-rule, HMO, Clar, resolver,
  aromatizer, validator, conformance, and installed ingestion tests preserve
  established structures while covering explicit selector propagation.
  `[dep: S9g, S9j, S9k]`
- **S9n — configurable graph-symmetry backend**
  (`umol-ast/src/ast/symmetry.rs`, all callers): add
  `automorphism_algorithm` to `GraphSymmetryConfig`, retain it in
  `GraphSymmetry`, and use it for initial symmetry, fixpoint reruns, and site
  stabilizers. Migrate every config literal in the same subitem. Tests pin the
  selected backend at all three call sites and preserve the existing Nauty
  result under an explicit selector. **Breaking config-shape migration
  (red→green).**
  **Implemented (green).** `GraphSymmetryConfig` now carries the explicit
  automorphism selector, and `GraphSymmetry` retains it for site-stabilizer
  reruns. Initial symmetry, fixpoint refinement, and site stabilization all
  pass the selected backend. Every Rust config literal now selects Nauty
  explicitly. Focused symmetry, stereo-view, model-config, and stereo property
  tests preserve the established orbit, chirality, and stereogenicity results
  while exercising the construction, fixpoint, and stabilizer paths under that
  selector. `[dep: —]`
- **S9o — validation operation configs**
  (`umol-graph/src/ops/{validate,validate/stereo}.rs`): add
  `StereoValidateConfig { automorphism_algorithm, max_iterations }`, with Nauty
  and 16 as inspectable defaults, and make `StereoConformanceValidator` consume
  it together with the semantic `StereoModel`. Add
  `ValidateConfig { aromaticity, stereo }`, where the aromaticity branch is the
  shared `AromaticityConfig`; `Validator::new(model)` delegates to defaults and
  `with_config(model, config)` passes both branches explicitly. Tests cover
  equality/defaults, both composite branches, and propagation into
  `GraphSymmetryConfig`. **Additive operation configs and constructor path
  (green).** `[dep: S9m, S9n]`
- **S9p — stereo model/resolve-config migration**
  (`umol-graph/src/ops/{model,resolve/stereo}.rs`,
  `umol-py/src/{model/stereo,resolve}.rs`, all callers): remove
  `max_iterations` and `inconsistency` from `StereoModel`. Move
  `InconsistencyPolicy` to the resolution module and add it to
  `StereoResolveConfig`, preserving `Error` as the high-level default. Move the
  Python wrapper to `resolve.rs`, update `StereoResolveConfig`, and remove both
  operational fields from Python `StereoModel`; the public
  `InconsistencyPolicy` name remains unchanged. Python does not bind the
  validation-only configs until it exposes a validation operation. Tests cover
  model/config equality and repr, Rust/Python conversion, all inconsistency
  policies, default ingestion parity, and removal of the obsolete model fields.
  **Breaking Rust and Python model/config migration (red→green).**
  `[dep: S9o]`
- **S9q — explicit reaction-composition enumeration**
  (`umol-ast/src/ast/compose.rs`, `umol-py/src/reaction.rs`, package callers):
  remove `CompositionScope` from Rust, its Python binding, registration, exports,
  and callers. Make `ReactionAst::compose` enumerate the complete admissible
  overlap set, including the empty overlap, and require
  `CommonSubgraphEnumerationAlgorithm` explicitly at the Rust boundary. Expose
  the selector directly in Python as the keyword-only `algorithm`, defaulting
  visibly to `Backtracking()`, and pass it through without another config type.
  Migrate all Rust and Python callers in the same subitem. Table and property
  tests cover explicit selector propagation, complete composite sets, and
  stereo-only and deletion-only reactions that the former reaction-center
  filter could empty. **Breaking Rust and Python API migration (red→green).**
  `[dep: S9j]`
- **S9r — pattern-fingerprint matching configuration**
  (`umol-graph/src/fingerprint/pattern.rs`,
  `umol-py/src/fingerprint/config.rs`): add `match_algorithm` and
  `subgraph_isomorphism_algorithm` to `PatternFingerprinter` and
  `PatternFingerprintConfig`, with `GraphAndOverlays` and VF2 as inspectable
  high-level defaults. Thread both values to template matching; update
  constructors, conversions, equality, repr, and all callers. Cross-algorithm
  tests require identical fingerprints for representative ground molecules.
  **Breaking config-shape migration (red→green).** `[dep: S9j]`
- **S9s — structural-fingerprint algorithm configuration**
  (`umol-graph/src/fingerprint/substructure.rs`,
  `umol-py/src/fingerprint/config.rs`): add `subgraph_enumeration_algorithm` and
  `automorphism_algorithm` to `SubstructureFeaturizer` and
  `StructuralFingerprintConfig`, with ESU and Nauty as inspectable high-level
  defaults. Pass the canonicalization selector into the canonical-key operation
  instead of closing over a literal. Update constructors, conversions, equality,
  repr, and all callers; tests pin selector propagation and unchanged structural
  feature sets. **Breaking config-shape migration (red→green).** `[dep: S9j]`
- **S9t — configured and fallible kekulization matching**
  (`umol-graph/src/ops/transform/kekulizer.rs`, `ops/transform.rs`, all callers):
  rename `KekulizationModel` to `KekulizationConfig` while retaining its local
  operation-level `MaximumMatchingAlgorithm`. Doc 158 S5d already removed the
  silent Edmonds substitution: both modes dispatch Edmonds to
  `general_maximum_matching` and Hopcroft--Karp to
  `bipartite_maximum_matching`, mapping `NonBipartiteGraphError` to
  `KekulizerError::NonBipartiteMatching(system)` without fallback. Preserve
  that behavior through the config rename and migrate all callers. Tests pin
  the renamed config construction, exact selector propagation, and the
  existing unchanged-input error contract. **Breaking config rename with
  caller migration (red→green).** `[dep: S9g, S9h]`
- **S9u — algorithm-transparency gate** (workspace): search every non-test Rust
  call to a graph-core selector-bearing operation. Require each selection to
  originate in a low-level method argument, a stored operational-config field,
  one of S9i's documented subsidiary choices, or an operation whose public name
  fixes that exact algorithm family and therefore offers no selection. Verify
  that `umol-io` remains free of hidden graph-core choices; run focused matching,
  matching-count, ring, symmetry, composition, fingerprint, aromaticity,
  validation, kekulization, reaction, and Python config suites, strict Clippy,
  formatting, benchmarks, and `git diff --check`. Record any future exception
  beside its call site rather than broadening a visibility boundary. **Additive
  gate (green).**
  `[dep: S9f, S9i, S9l, S9m, S9p, S9q, S9r, S9s, S9t]`

S9 preserves the existing direct `apply` happy path while correcting error
classification, separates chemistry models from operational configuration, and
removes every audited hidden graph-core algorithm selection. Edmonds remains the
general-graph matching default; another non-bipartite algorithm is not required
by S9, and Gabow remains contingent on corpus benchmarks. S9 adds no `apply_at`,
report object, or prepared-reaction type. No S9 subitem is deferrable.

The follow-up comparison vocabulary, panic-removal plan, and property-suite
reorganization are specified in [doc 156](156-ast-comparison-and-property-suite-2026-07-20.md).

### S10 — Public integration and release gate

- **S10a — public surface audit**
  (`umol-py/src/lib.rs`, `python/umol/__init__.py`): audit native registration,
  package imports, `__all__`, constructor visibility, signatures, reprs, and
  exception mapping for every S4–S9 value and method. Installed tests compare the
  exported-name set and verify that no lint, raw-seed, combinator, iterator,
  NumPy, DRFP, or BRIDGIT implementation leaked into this round. **Additive
  (green).** `[dep: S7d, S8b, S9u]`
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
`S1h → S2a → S2b → S4a → S4b → {S4c → S4d}`, then continues through doc 155
S0–S7 to the configured `from_smiles` operation.

The fingerprint path is:

`S0a → S3a → S3b`, with `S3c` joining at `S3d → S4a`; the parallel configuration
and result branches `{S5, S6}` join at `S7 → S8`.

The application path is:

`S9a → S9d`, while `S4a → S9b → S9c`; both join with S4b at
`S9e → S9f`. The transparency foundation is `S9g → S9h → S9i`, with
`S9g → S9j` and `S9k → {S9l, S9m}`. The consumer branches are
`{S9j, S9k} → S9m`, `{S9m, S9n} → S9o → S9p`,
`S9j → {S9q, S9r, S9s}`, and `{S9g, S9h} → S9t`; all join the application
path at `S9u`. The three
deliverable paths join at
`S10a → S10b → S10c`.

The following are explicitly deferrable and are not on the required critical
path: NumPy integration; `SmilesLintFlags`/`SmilesLintConfig`; Python
`ChemistryModel`; the unresolved `umol-io` AST parser; raw/custom hash schemes;
DRFP and BRIDGIT; generic `transact_validated`; validator combinators;
transaction scopes, savepoints, and resolver fallback; general transformer
atomicity; `apply_at` and diagnostic application reports. The complete S10
workflows do not depend on any of them.
