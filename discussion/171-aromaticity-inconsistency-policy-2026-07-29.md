# 171 — Aromaticity resolution and conformance

Status: In Progress
Date: 2026-07-29
Relates: [166](166-molecule-ops-2026-07-27.md),
[177](177-nomenclature-guide-2026-07-31.md)

Aromatic input that the chemistry model declines currently leaves projections behind with no relation
to project from. This document states the problem, the resolver policy that follows from
kekulization being an explicit operation, and the related charge-delocalization and conformance work
that must land as one coherent unit.

## Symptom

`MoleculeAst.from_smiles("c1ccoc1", chemistry_model=<mdl>)` returns a molecule with **no** aromatic
system whose atoms still carry `#a` (`#a2` on the oxygen) and whose bonds still carry `#a`.
`reset_aromatic_valence = true` does not change the output. Verified 2026-07-29, Python bindings:

| SMILES | `daylight()` | `mdl()` |
| --- | --- | --- |
| furan `c1ccoc1` | 1 aromatic system | 0 |
| thiophene `c1ccsc1` | 1 | 0 |
| pyrrole `c1cc[nH]c1` | 1 | 0 |
| pyridine `c1ccncc1` | 1 | 1 |
| benzene `c1ccccc1` | 1 | 1 |

`mdl()` is `HueckelRule` over {C, N} with `min_ring_size: 6`, so furan fails on both the element scope
and the ring size.

Rendered output for furan under `mdl()`:

```clojure
{:atoms ["C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" … "O#i=#c0#h0#n#u0#s#v2#d0#t0#a2#m!" …]
 :bonds [[0 4 "1#c0#u0#s#a"] [0 1 "1#c0#u0#s#a"] …]}
```

## Why it is wrong

Per-atom `#a` (aromatic valence) and per-bond `#a` (aromatic incidence) are projections of an
aromatic-system relation onto its participants. With no system present they project from nothing. This
contradicts the entity model: a derivable constraint whose relation is absent has no derived side.

## Why the obvious repair does not work

Aromaticity is the only overlay carrying information that the localized bonds would otherwise have to
hold. In the output above every ring bond is order 1 and every ring carbon is `#v2 #h1 #a1`. Strip the
`#a` and each ring carbon has three bonds' worth of valence and needs four, so **reset alone yields an
invalid molecule**. Discarding aromaticity requires *supplying* the alternating bond orders it stood
in for.

Stereo does not have the same bond-order dependency. Removing an unrealizable stereo assertion leaves
the constitution and localized bond orders intact, so removing stereo information is a sound explicit
stereo-resolution action. The policy must still say whether it removes the projection, the entity, or
both. A shared policy enum would expose actions that are unsound for aromaticity.

## Constraint: perception is flag-driven

Aromaticity perception in the ingest path is driven by the input's per-atom aromatic flags, not by ring
analysis of a Kekulé structure. Verified: `C1=CC=CC=C1` yields **zero** aromatic systems under every
model, benzene included, while `c1ccccc1` yields one under both `daylight()` and `mdl()`.

Consequence: once the flags are discarded there is no flag-driven path that rediscovers the system.
This rules out `Strip`. It does not rule out `Keep`: retaining the assertions preserves the information
needed by a later operation or a different chemistry model.

## Representation-consistency policy

The base entity and its stored projections are independent inputs. Neither is authoritative merely
because it is materialized or because it appears in a constraint container. Derivation evaluates the
two independently and classifies their relationship; only resolution applies a provenance-sensitive
policy.

Three failure classes are distinct:

| Failure | Meaning |
| --- | --- |
| Constraint failure | A non-vacuous constraint cannot produce a valid entity under the selected topology and model |
| Entity failure | A structurally readable entity is not realizable under the selected topology and model |
| Entity/constraint mismatch | The entity and constraint are independently realizable but disagree |

Malformed representation shape is outside this policy. Dangling participants, participant/data
length mismatches, duplicate sites, and similar conditions are unconditional
`EntityStructureValidator` contradictions. A resolver must not interpret malformed storage in order
to choose which information to retain.

Resolver policies are separate per failure class and per operation. A bare `Strip` does not say which
representation is discarded and is therefore not an adequate policy name.

| Failure | Aromaticity | Stereo |
| --- | --- | --- |
| Constraint failure | error or keep | error, keep, or remove the constraint |
| Entity failure | error or keep | error, keep, or remove the entity |
| Entity/constraint mismatch | error, keep, remove the constraint, or replace the entity from the constraint | error, keep, remove the constraint, replace the entity from the constraint, or remove both |

Replacing an entity is one atomic plan: remove the conflicting relation and materialize the valid
constraint-derived replacement. Removing the entity while retaining the constraint is not a stable
resolution result because a subsequent call would recreate it. Removing a constraint retains the
independently valid entity and clears only the redundant assertion.

Aromaticity admits fewer recovery actions because the aromatic-system relation carries bond
information not otherwise present in the localized representation:

- an unmatched constraint cannot simply be removed;
- an internally invalid entity cannot be removed without a valid replacement;
- an entity/constraint mismatch may replace the entity only when the constraint independently yields
  a valid complete system;
- removing both sides is never sound.

Stereo relations do not carry constitutional bond information, so removing an invalid stereo entity
or removing both sides of a mismatch is sound. This does not make either representation authoritative;
it only makes more explicit recovery choices available.

The same symmetry applies to model-independent projections. Localized topology versus `#v`, dative
relations versus `#d`/`#t`, and multicenter relations versus `#m` must never be silently reconciled.
Their scalar projections cannot reconstruct the contributing bonds, so the generic recovery choices
are error, keep, or remove the projection; removing the underlying relation is a separate structural
operation.

For now, the dative rule applies only to binary donor→acceptor bonds. Multi-donor entries conflate
coordination/haptic bonding with binary dative bonding; their `#d`/`#t` projection remains an explicit
stub pending the entity split and extensibility work in
[doc 117](117-entity-model-extensibility-2026-06-20.md). This deferral limits coverage but does not
change the representation-consistency policy.

The policy covers all entity kinds without requiring a generic repair operation:

- atom and bond constraints are checked against their `AtomView` and `BondView` projections;
- dative, aromatic, multicenter, noncovalent, stereo-atom, and stereo-bond constraints are checked
  against their corresponding relation views;
- entity-local aggregate constraints, such as aromatic-system and multicenter-bond electron counts,
  are checked against the entity fields from which they are derived;
- relational and molecule constraints are evaluated against the complete molecule;
- ring constraints are evaluated against the explicitly selected fixed ring projection.

Where no resolver can reconstruct a sound replacement, this policy adds validation but not an
automatic mutation. A caller may still remove a constraint or structural entity through the explicit
editing API. Resolution policy is added only to an operation that can define a stable result for every
offered variant.

Validators have no inconsistency policy. Every determined constraint failure, entity failure, and
entity/constraint mismatch is `Contradictory`. Absent and `Undetermined` constraints are vacuous. A
non-vacuous assertion that cannot yet be decided is `Underdetermined`.

**Transformations are explicit operations, not policy variants or reaction primitives.** Resolution
and validation never invoke them implicitly. Kekulization and aromatization are the transformations
at issue here; in particular, kekulization has no simple inverse that would make it a sound primitive
reaction operation.

Raw parse/raise has no aromatic-system relation for `Kekulizer` to consume. A caller who wants furan
localized before applying `mdl()` must first resolve it under a model that accepts the source aromatic
form, then run the explicit transformation:

```text
parse/raise → resolve(source model) → kekulize → resolve/validate(target model)
```

The current Python surface does not expose this full chain.

## Derivation and verification

The consistency classification belongs with the derivation that has enough information to evaluate
both representations. Policy is applied only by resolvers:

- the constraint is evaluated without treating an existing entity as evidence that the constraint is
  realizable;
- the entity is evaluated without treating its constraint as evidence that the entity conforms;
- only after both evaluations does derivation classify a constraint failure, entity failure, or
  entity/constraint mismatch;
- The `*Perception` types carry the model and perform the operation. The `*Derivation` values are
  policy-free results consumed by resolvers and validators.
- each inconsistency identifies both the failure class and the affected constraint site or entity ID;
- resolvers apply the separately configured policy for that failure class before constructing one
  atomic plan;
- `AromaticityConformanceValidator` and `StereoConformanceValidator` consume the same policy-free
  derivation results but always report mismatches as contradictions.

The aromatic and stereo derivations distinguish an independently unrealizable constraint or entity
from two independently realizable representations that disagree. Stereo relations are evaluated by
deriving the canonical ligand frame from the molecule and mapping the relation's stored frame into
it; relation presence is not evidence of conformance. An explicit `NotStereo` constraint and a valid
stereo relation therefore form a mismatch.

`Undetermined` constraints are vacuous and do not assert that an entity is missing. Validators ignore
an absent or `Undetermined` constraint. A non-vacuous assertion without a matching relation, or an
explicit constraint that contradicts its relation, is `Contradictory`. A non-ground but non-vacuous
assertion that cannot yet be decided is `Underdetermined`.

`reset_aromatic_valence` clears the aromatic-valence constraints on atoms of newly materialized
systems by setting them to `Undetermined`. It is not the lever for a system the model declines.

### Current perception API

`AromaticityPerception` retains `find_systems` as the low-level algorithmic operation whose caller
supplies the π-electron source. `derive` is the standard AST-facing operation: it reads aromatic
assertions, calls `find_systems`, and compares the accepted systems with existing relations and
non-vacuous constraints.

```rust
pub struct AromaticityDerivation {
    pub systems: Vec<(Vec<AtomId>, AromaticSystemAst)>,
    pub inconsistencies: Vec<AromaticityInconsistency>,
}

pub enum AromaticityInconsistency {
    AromaticValenceFailure { atom: AtomId },
    AromaticSystemFailure { system: AromaticSystemId },
    AromaticValenceMismatch {
        atom: AtomId,
        system: AromaticSystemId,
    },
    AromaticBondConstraintMismatch {
        bond: BondId,
        system: AromaticSystemId,
    },
}

impl AromaticityPerception {
    pub fn new(model: &AromaticityModel) -> Self;

    pub fn find_systems<F>(
        &self,
        ast: &MoleculeAst,
        config: AromaticityConfig,
        electrons_at: F,
    ) -> Result<
        Solution<Vec<(Vec<AtomId>, AromaticSystemAst)>, AromaticityContradiction>,
        AromaticityError,
    >
    where
        F: Fn(&AtomView<'_>) -> Option<u8>;

    pub fn derive(
        &self,
        ast: &MoleculeAst,
        config: AromaticityConfig,
    ) -> Result<
        Solution<AromaticityDerivation, AromaticityContradiction>,
        AromaticityError,
    >;
}
```

`AromaticityResolver` and `AromaticityConformanceValidator` consume `derive`.
`Aromatizer` calls `find_systems` directly with its Kekulé electron source. The current mutating
`AromaticityPerception::add_systems` is removed; resolvers and transformers own materialization.

`StereoPerception::derive_stereo_atom` and `derive_stereo_bond` are the public per-entity operations.
The molecule-wide `derive` calls them for non-vacuous `#T` and `#C` assertions and compares their
results with existing stereo relations.

```rust
pub struct StereoPerception {
    model: StereoModel,
}

pub struct StereoDerivation {
    pub atoms: Vec<(AtomId, Vec<StereoLigand>, StereoAtomAst)>,
    pub bonds: Vec<(BondId, Vec<StereoLigand>, StereoBondAst)>,
    pub inconsistencies: Vec<StereoInconsistency>,
}

pub enum StereoInconsistency {
    TetrahedralStereoFailure { atom: AtomId },
    StereoAtomFailure { stereo_atom: StereoAtomId },
    TetrahedralStereoMismatch {
        atom: AtomId,
        stereo_atom: StereoAtomId,
    },
    CisTransStereoFailure { bond: BondId },
    StereoBondFailure { stereo_bond: StereoBondId },
    CisTransStereoMismatch {
        bond: BondId,
        stereo_bond: StereoBondId,
    },
}

impl StereoPerception {
    pub fn new(model: &StereoModel) -> Self;

    pub fn derive(&self, ast: &MoleculeAst) -> StereoDerivation;

    pub fn derive_stereo_atom(
        &self,
        ast: &MoleculeAst,
        atom: AtomId,
        coset: &StereoCoset,
    ) -> Option<(Vec<StereoLigand>, StereoAtomAst)>;

    pub fn derive_stereo_bond(
        &self,
        ast: &MoleculeAst,
        bond: BondId,
        coset: &StereoCoset,
    ) -> Option<(Vec<StereoLigand>, StereoBondAst)>;
}
```

`StereoResolver` applies the six operation-specific policy fields to `StereoDerivation` before
constructing one edit plan. `Error` returns the exact inconsistency before any edits are applied;
`Remove` removes only the independently invalid side; and mismatch policies can retain both sides,
remove the constraint, replace the relation from the constraint, or remove both. For `NotStereo`,
replacement removes the relation and adds no new relation. `StereoConformanceValidator` integration
is completed separately in S5b. The helper that derives the two ligands at one end of a cis-trans
bond remains private.

### Approved inconsistency and policy names

The settled meanings of constraint, failure, mismatch, inconsistency, policy, contradiction, error,
and the associated operation names are recorded in the repository-wide
[nomenclature guide](177-nomenclature-guide-2026-07-31.md#glossary).

Diagnostic variants and config fields use the concrete constraint and entity names. They do not
introduce a generic `Projection` or `Representation` vocabulary, and they do not force parallel
wording where the underlying constraint and entity names differ. `AromaticBondConstraint` refers to
`BondConstraintAst::Aromatic`; it does not introduce another AST type.

The following blocks fix the public names. Their diagnostic fields show the required identity
information but may acquire further detail when the classifications are implemented.

```rust
pub enum AromaticityInconsistency {
    AromaticValenceFailure { atom: AtomId },
    AromaticSystemFailure { system: AromaticSystemId },
    AromaticValenceMismatch {
        atom: AtomId,
        system: AromaticSystemId,
    },
    AromaticBondConstraintMismatch {
        bond: BondId,
        system: AromaticSystemId,
    },
}

pub enum AromaticityFailurePolicy {
    Error,
    Keep,
}

pub enum AromaticityMismatchPolicy {
    Error,
    Keep,
    RemoveConstraint,
    ReplaceEntity,
}

pub enum AromaticBondConstraintMismatchPolicy {
    Error,
    Keep,
    RemoveConstraint,
}
```

An aromatic bond constraint cannot derive a replacement aromatic system, so its mismatch policy does
not expose `ReplaceEntity`. The dedicated enum keeps every offered recovery action total rather than
making one variant context-dependent.

Electron-contribution failures within an aromatic system are
`AromaticityInconsistency::AromaticSystemFailure`.
`AromaticSystemConstraintAst::ElectronCount` belongs to constraint validation, not aromaticity
resolution.

```rust
pub enum StereoInconsistency {
    TetrahedralStereoFailure { atom: AtomId },
    StereoAtomFailure { stereo_atom: StereoAtomId },
    TetrahedralStereoMismatch {
        atom: AtomId,
        stereo_atom: StereoAtomId,
    },
    CisTransStereoFailure { bond: BondId },
    StereoBondFailure { stereo_bond: StereoBondId },
    CisTransStereoMismatch {
        bond: BondId,
        stereo_bond: StereoBondId,
    },
}

pub enum StereoFailurePolicy {
    Error,
    Keep,
    Remove,
}

pub enum StereoMismatchPolicy {
    Error,
    Keep,
    RemoveConstraint,
    ReplaceEntity,
    RemoveBoth,
}
```

Policy enums are shared where their action sets are identical; the config field identifies the
specific constraint or entity affected by `Remove`, `RemoveConstraint`, `ReplaceEntity`, or
`RemoveBoth`. All policies default to `Error`.

```rust
pub struct AromaticityResolveConfig {
    pub perception: AromaticityConfig,
    pub aromatic_valence_failure: AromaticityFailurePolicy,
    pub aromatic_system_failure: AromaticityFailurePolicy,
    pub aromatic_valence_mismatch: AromaticityMismatchPolicy,
    pub aromatic_bond_constraint_mismatch: AromaticBondConstraintMismatchPolicy,
    pub reset_aromatic_valence: bool,
}

pub struct StereoResolveConfig {
    pub tetrahedral_stereo_failure: StereoFailurePolicy,
    pub stereo_atom_failure: StereoFailurePolicy,
    pub tetrahedral_stereo_mismatch: StereoMismatchPolicy,
    pub cis_trans_stereo_failure: StereoFailurePolicy,
    pub stereo_bond_failure: StereoFailurePolicy,
    pub cis_trans_stereo_mismatch: StereoMismatchPolicy,
    pub reset_stereo_constraints: bool,
}
```

## Resolver boundary

Charge delocalization moves out of `AromaticityResolver` into the explicit transformer specified
below, and `AromaticityResolveConfig::delocalize_charge` is retired in the same public-config
migration that adds the aromatic inconsistency policy.

The lasting resolver invariant is that resolution does not change a localized bond order, add or
remove an atom or localized bond, or otherwise alter the constitution. Kekulization and aromatization
change bond orders, which is why they remain caller-invoked transformations.

## Charge delocalization

Move charge delocalization out of `AromaticityResolver` and into a
model-independent `DelocalizeCharge` transformer under `umol-graph/src/ops/transform`. Retire
`AromaticityResolveConfig::delocalize_charge`.

Delocalization rewrites one resolved representation into another; it does not fill undetermined
state. For an aromatic system whose participants carry literal formal charges, the transformer moves
the summed charge onto the system, sets the contributing atoms to literal zero, and adjusts their
π-electron contributions so the total is preserved. It leaves undetermined or non-literal charges
unchanged and is idempotent.

No `LocalizeCharge` inverse is added: the choice of which atom receives the charge is not canonical,
so localization is not a function of the delocalized structure.

The localized and delocalized representations are observably distinct. A `"C#c-1"` pattern matches
the localized form and not the delocalized form; `"C#a2"` and `"C#a1"` likewise distinguish their
stored contributions. The choice cannot be made silently.

Consequently, `ingest::smiles` returns the localized form for inputs such as `[cH-]1cccc1` unless the
caller explicitly applies `DelocalizeCharge`. Resolution may still add aromatic and stereo
relations. With the explicit projection-reset options disabled, it changes no already determined
field or constraint.

## Aromaticity conformance

`AromaticityConformanceValidator` validates every stored representation of aromaticity against
independent perception:

- stored aromatic-system participant sets match the perceived sets;
- stored per-atom electron contributions match the perceived contributions;
- every non-vacuous localized-bond aromatic assertion agrees with the bonds induced by the stored
  systems;
- every non-vacuous aromatic-valence assertion agrees with the corresponding system contribution;
- an asserted aromatic atom or bond without a matching system, an extra stored system, or an explicit
  negative assertion on a system participant is a contradiction;
- absent and `Undetermined` projections are ignored;
- a non-ground, non-vacuous contribution that cannot be decided yields `Underdetermined`, not a false
  match.

This work was previously listed in doc 166 and belongs here because the resolver and validator must
use the same policy-free constraint comparison.

## Validator ownership

Constraint evaluation and model conformance are separate responsibilities:

- `EntityStructureValidator` rejects malformed entity storage unconditionally.
- `IncidenceConstraintValidator` evaluates entity-local constraints derived from fields and directly
  incident localized bonds or overlay relations. Its focused AST-layer operation receives a
  `ConnectedComponentsAlgorithm` for noncovalent `#I`; every other incidence constraint is
  algorithm-free.
- `RingConstraintValidator` evaluates ring membership, ring degree, and ring valence under the fixed
  Relevant-through-22 semantics. Its focused AST-layer operation receives a
  `RelevantCycleEnumerationAlgorithm` directly.
- `ConstraintValidator` coordinates model-independent, closed-world constraint evaluation over a
  resolved molecule and folds `And`/`Or`/`Not`. A logical tree may mix incidence, ring, relational,
  and molecule-scope leaves, so the complete fold cannot be split into independent validator
  outcomes.
- `AromaticityConformanceValidator` and `StereoConformanceValidator` own model-dependent entity
  realizability and entity/constraint classification.
- the composite `Validator` preserves the integrity → invariants → conformance order.

Thus an atom's `#v`, `#d`, `#t`, `#a`, `#m`, and `#T` are incidence constraints evaluated through
`AtomView`; `#a` and `#C` on localized bonds are evaluated through `BondView`. These values may be
Boolean, counts, or weighted sums: incidence describes their structural source, not their value
type. Ring-derived variants remain stored in the entity constraint containers but dispatch to the
ring evaluator.

Stereo-atom and stereo-bond constraint leaves have two distinct uses. The AST constraint evaluator
can compare their references and stored values while folding a logical constraint tree. Chemical
claims such as ligand symmetry, topicity, and stereogenicity are evaluated only by
`StereoConformanceValidator`, using `StereoValidateConfig` and the selected `StereoModel`; the AST
validator does not duplicate graph-symmetry perception. The same division applies generally:
model-independent constraint evaluation establishes what the AST asserts, while conformance decides
whether a model accepts the assertion.

The contradiction hierarchy follows evaluation mechanism first, then the concrete entity and
constraint name. Ring contradictions therefore do not sit among incidence contradictions merely
because both constraints are stored on an atom or bond. `RelationalConstraint` retains its narrower
meaning of a reference-bearing molecule-scope constraint; it is not another name for an incidence
constraint.

Ring enumeration is performed once and only when the logical tree or inline containers contain a
ring-derived constraint. Dative-bond ring membership remains the explicit unimplemented case pending
the dative versus coordination/haptic topology decision in doc 117.

The direct specialized operations receive their algorithm selectors explicitly. Composite
constraint validation uses structured operation config rather than flattening every subsidiary
selector into one method signature. The AST-layer config has no default because its selectors are
algorithmically transparent:

```rust
pub struct ConstraintValidateConfig {
    pub relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
    pub connected_components_algorithm: ConnectedComponentsAlgorithm,
    pub substructure_match_algorithm: SubstructureMatchAlgorithm,
    pub subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm,
}
```

`ConstraintValidator` stores this config. The higher-level `umol_graph::ValidateConfig` nests it as
`constraint` and supplies documented defaults. A focused incidence validator receives
`ConnectedComponentsAlgorithm` directly for `#I`; a focused ring validator receives
`RelevantCycleEnumerationAlgorithm` directly. The complete constraint validator runs ring
enumeration, connected-components analysis, or substructure matching only when a present non-vacuous
leaf requires it.

The dedicated aromaticity and stereo validators then answer the separate model-conformance question.
Doc 166 now points to this work unit for the complete constraint-integrity implementation.

## Open

- Whether the contradiction should name the offending atoms and the reason (element out of scope
  against ring size), beyond the exact unmatched atom and bond IDs required initially.
- Which high-level Python transformation surface should expose the source-resolve → kekulize → target
  resolve/validate chain.

## Staged implementation plan

Every Rust subitem carries focused `#[rstest]` coverage and leaves its affected crate green; Python
subitems carry focused pytest coverage. Tests use exact `Solution`, contradiction, edit-plan, error,
or molecule equality rather than summary-only assertions.

S0 and S1 record completed foundation work. They remain part of the dependency graph, but their
initial mismatch enums and single-policy resolver configs are not the final public contract. S2
performs those breaking migrations explicitly; the completed subitems are not relabeled as though
they had already implemented the refined classifications.

### S0 — Policy-free foundations

#### S0a — Charge-delocalization transformer **Done**

**Module:** `umol-graph/src/ops/transform/delocalize_charge.rs`,
`umol-graph/src/ops/transform.rs`

**Kind:** additive (green)

**Dependencies:** `[dep: none]`

Add public `DelocalizeCharge` implementing `Transformer`. Move the model-independent charge and
π-contribution calculation behind this transformer without changing the existing resolver path yet.
Use `Infallible` as its error type: systems with non-literal data are left unchanged. Derive every
selected system update from the immutable input before applying the complete transformation.

Focused tables cover cyclopentadienyl anion, tropylium cation, a heterogeneous aromatic system,
non-literal input, and multiple systems. They assert charge conservation, π-electron conservation,
field preservation, and idempotence.

#### S0b — Aromatic perception and derivation **Done**

**Module:** `umol-graph/src/ops/aromaticity.rs`

**Kind:** additive (green)

**Dependencies:** `[dep: none]`

Add `AromaticityPerception::derive`, returning public `AromaticityDerivation`. The result contains
the systems independently accepted by the selected model together with sorted unmatched atom, bond,
and existing-system IDs. Compare non-vacuous atom aromatic-valence assertions, non-vacuous
localized-bond aromatic assertions, stored aromatic-system participant sets, and stored electron
contributions. Existing systems must be re-perceived; their presence does not count as acceptance.
Absent and `Undetermined` projections are ignored.

Focused tables cover MDL furan rejection, Daylight furan acceptance, missing and extra atom/bond
assertions, vacuous projections, contribution mismatch, an already conformant system, and an
existing system rejected by the selected model. They assert the complete `AromaticityDerivation`.

#### S0c — Stereo perception and derivation **Done**

**Module:** `umol-graph/src/ops/stereo.rs`, `umol-graph/src/ops.rs`

**Kind:** additive (green)

**Dependencies:** `[dep: none]`

Add public model-carrying `StereoPerception` and policy-free result `StereoDerivation`. Move the
current atom ligand-frame, bond ligand-frame, and side-ligand calculations into
public `StereoPerception::derive_stereo_atom` and `derive_stereo_bond`; only the side-ligand helper
remains private. The molecule-wide `derive` calls the per-entity methods, leaving the resolver wired
to its existing methods until the new implementation has parity. The result contains realizable
materializations and exact unrealizable or mismatched non-vacuous `#T` and `#C` sites without
choosing `Keep`, `Strip`, or `Error`. It also compares an asserted site already covered by a stereo
relation instead of skipping it.

Ordered table tests cover successful tetrahedral and cis-trans derivation, disabled kinds, element
scope, ligand arity, aromatic exclusion, existing elements, and exact atom/bond contradictions.

**Stage exit:** charge delocalization is available as an explicit transformer, and aromaticity and
stereo can derive and compare their relations without an operational policy.
`cargo test -p umol-graph` passes.

### S1 — Resolver contract migration

#### S1a — Aromatic resolver and Python configuration **Done**

**Module:** `umol-graph/src/ops/aromaticity.rs`,
`umol-graph/src/ops/resolve/aromaticity.rs`, `umol-graph/src/ops/resolve.rs`,
`umol-graph/src/parse.rs`, `umol-graph/src/ingest.rs`, `umol-py/src/resolve.rs`,
`umol-py/src/lib.rs`, and the affected Rust and Python configuration tests

**Kind:** breaking public config migration (red → green)

**Dependencies:** `[dep: S0a, S0b]`

Add `AromaticityInconsistencyPolicy::{Keep, Error}` beside `AromaticityResolveConfig`, defaulting to
`Error`, and add it as the config's `inconsistency` field. In the same public-config migration, remove
`delocalize_charge`; `AromaticityPerception::add_systems` and `AromaticityResolver` must no longer
equalize charge implicitly.

`AromaticityResolver::plan` consumes the S0b derivation before constructing edits. `Error` returns the
exact contradiction without mutation. `Keep` retains unmatched projections and adds only independently
accepted systems that are not already present; it also retains an existing system rejected by the
selected model. Repeated resolution of a conformant molecule is an identity. Add a
composite-resolver case in `ops/resolve.rs` proving that an aromatic contradiction rolls back edits
already applied by valence resolution.

Migrate all Rust config literals and the Python wrapper in the same subitem. Python exposes
`AromaticityInconsistencyPolicy`, makes `inconsistency` keyword-only, removes
`delocalize_charge`, and updates constructor signatures, getters, conversions, equality, repr, import
coverage, and all affected workflow/molecule/reaction tests. The localized post-ingest charge
representation becomes the new exact expected value.

#### S1b — Stereo policy rename **Done**

**Module:** `umol-graph/src/ops/resolve/stereo.rs`, `umol-graph/src/ops/resolve.rs`,
`umol-graph/src/ingest.rs`, `umol-py/src/resolve.rs`, `umol-py/src/lib.rs`, and the affected Rust and
Python configuration tests

**Kind:** breaking public type rename (red → green)

**Dependencies:** `[dep: S0c]`

Rename the existing public `InconsistencyPolicy` to `StereoInconsistencyPolicy` in Rust and Python.
Retain `Keep`, `Strip`, and `Error`, retain `Error` as the default, and change no stereo-resolution
behavior. Migrate every Rust import/config literal and every Python import, annotation, conversion,
repr, and exact enum/config test in the same subitem. Rewire `StereoResolver` to the S0c derivation
through `StereoPerception` without changing its planned edits.

**Stage exit:** the initial operation-specific config migration is green. S2 performs the final
failure-class and recovery-action split.
`cargo test -p umol-graph`, `cargo test -p umol-py`, `maturin develop`, and
`pytest -q umol-py/tests` pass with Python 3.13 active.

### S2 — Final aromaticity and stereo contracts

#### S2a — Aromatic inconsistency classification and policy migration **Done**

**Module:** `umol-graph/src/ops/aromaticity.rs`,
`umol-graph/src/ops/resolve/aromaticity.rs`,
`umol-graph/src/ops/validate/aromaticity.rs`, `umol-py/src/resolve.rs`,
`umol-py/src/lib.rs`, and affected Rust and Python callers

**Kind:** breaking public type and config migration (red → green)

**Dependencies:** `[dep: S0a, S0b, S1a]`

Replace `AromaticityMismatch` with `AromaticityInconsistency` and independently derive constraint
failure, entity failure, and entity/constraint mismatch. An existing system is re-evaluated from its
own participants and contributions; it is not accepted merely because perception finds a system on
the same atoms. Preserve deterministic ordering and include the exact atom, bond, and system IDs
fixed above.

Replace the single resolver policy with `AromaticityFailurePolicy`,
`AromaticityMismatchPolicy`, `AromaticBondConstraintMismatchPolicy`, and the four approved config
fields. The bond-specific policy omits `ReplaceEntity` because an aromatic bond constraint cannot
derive a complete replacement system. Implement `Keep`, `RemoveConstraint`, and atomic
`ReplaceEntity` only for the failure classes whose action tables admit them. `Error` returns the exact
inconsistency as a contradiction before mutation. Remove the old policy type and migrate the Python
enum/config surface, conversions, repr, equality, exports, and callers in the same subitem.

Focused tables assert the complete derivation and complete planned edits for each inconsistency and
policy. They include independently invalid constraints and systems, valid-but-different systems,
atom and bond mismatches, reset interaction, idempotent conformant input, and transaction identity on
error. Python tables assert exact enum/config values and exact resolved or contradictory outcomes.

#### S2b — Stereo inconsistency classification and policy migration **Done**

**Module:** `umol-graph/src/ops/stereo.rs`, `umol-graph/src/ops/resolve/stereo.rs`,
`umol-graph/src/ops/validate/stereo.rs`, `umol-py/src/resolve.rs`,
`umol-py/src/lib.rs`, and affected Rust and Python callers

**Kind:** breaking public type and config migration (red → green)

**Dependencies:** `[dep: S0c, S1b]`

Replace `StereoMismatch` with `StereoInconsistency` and distinguish tetrahedral/cis-trans constraint
failure, stereo-atom/stereo-bond entity failure, and independently realizable mismatches. Existing
relations are assessed through the same public per-entity derivation operations as uncovered
constraint sites; relation presence is not evidence of conformance.

Replace `StereoInconsistencyPolicy` with `StereoFailurePolicy`, `StereoMismatchPolicy`, and the six
approved config fields. Implement `Remove` for the exact failed side identified by its config field,
and implement `RemoveConstraint`, atomic `ReplaceEntity`, and `RemoveBoth` for mismatches. Remove the
old behavior that silently ignored relation mismatches. Migrate the Python enum/config surface and
all callers in the same subitem.

Focused tables assert complete derivations and planned edits for atom and bond cases under every
applicable action, including reset interaction, conformant identity, and unchanged input on error.
Python tables assert exact config conversion and resolution outcomes.

#### S2c — Composite resolver transaction behavior **Done**

**Module:** `umol-graph/src/ops/resolve.rs`

**Kind:** behavioral correction (green)

**Dependencies:** `[dep: S2a, S2b]`

Update the composite resolver to propagate the refined aromaticity and stereo contradictions without
losing their payloads. Verify that a late inconsistency rolls back edits from every earlier resolver
stage, while successful stages are materialized before the next stage derives its plan. Cover both
aromaticity and stereo failures with exact source-molecule equality after rollback and exact final
molecule equality after success.

**Stage exit:** derivation, resolution, validation consumers, Rust callers, and Python bindings use
only the final inconsistency and policy types. The workspace is green; the old mismatch and policy
types no longer exist.

### S3 — Model-independent constraint evaluation

#### S3a — Configuration and contradiction foundations **Done**

**Module:** `umol-ast/src/ast/validate/constraint.rs`,
`umol-ast/src/ast/validate/constraint/`, `umol-ast/src/ast.rs`

**Kind:** additive (green)

**Dependencies:** `[dep: none]`

Add public `ConstraintValidateConfig`, `IncidenceConstraintValidator`,
`IncidenceConstraintContradiction`, `RingConstraintValidator`, and
`RingConstraintContradiction`. Establish the `ConstraintContradiction` hierarchy by evaluation
mechanism and concrete entity/constraint identity. Keep operational failures in `ConstraintError`;
semantic failures remain on the `Solution::Contradictory` side. Do not add defaults at the AST layer
and do not add one-field config wrappers to focused validators.

Focused tables establish construction and exact contradiction/error wrapping. No test merely checks
that an error or contradiction is present.

#### S3b — Entity incidence constraint evaluation **Done**

**Module:** `umol-ast/src/ast/validate/constraint/incidence.rs` and focused evaluator modules under
`umol-ast/src/ast/validate/constraint/`

**Kind:** additive (green)

**Dependencies:** `[dep: S3a]`

Implement the focused incidence evaluator over `AtomView`, `BondView`, and overlay views. Cover atom
valence, donated and accepted pairs, aromatic and multicenter valence, tetrahedral stereo incidence,
degree, total degree, total valence, and total hydrogens; localized-bond aromatic and cis-trans
incidence; binary dative aromatic incidence; noncovalent intramolecular status using an explicit
`ConnectedComponentsAlgorithm`; and the algorithm-free entity aggregates for aromatic-system and
multicenter-bond electron counts. The multi-donor `#d`/`#t` case remains an explicit unimplemented
branch pointing to doc 117 rather than acquiring provisional semantics.

Evaluate literal, finite-set, range, and undetermined values through their lattice semantics. Vacuous
constraints are determined identities; non-vacuous values whose derived side is unavailable are
underdetermined; decided disagreement returns the exact contradiction. Focused tables cover every
constraint variant and all three `Solution` outcomes. Complement-set forms belong to the finite
stereo relation lattices evaluated in the later relational/conformance stages; the incidence leaf
types do not define one.

#### S3c — Ring constraint evaluation **Done**

**Module:** `umol-ast/src/ast/validate/constraint/ring.rs`

**Kind:** additive (green)

**Dependencies:** `[dep: S3a]`

Implement `RingConstraintValidator` for atom ring degree, ring valence, and membership and localized
bond membership under the fixed Relevant-through-22 semantics. Its public focused operation takes
`RelevantCycleEnumerationAlgorithm` directly. Build one ring projection per validation call and
reuse it across every requested size and entity. Do not attempt dative-bond ring membership; return
the explicit unsupported operational error tied to doc 117.

Focused tables cover acyclic, monocyclic, fused, bridged, size-filtered, vacuous, non-ground, and
contradictory cases. The current Vismara case verifies explicit selector forwarding; future
relevant-cycle implementations join the same exact-outcome table.

#### S3d — Relational constraint evaluation **Done**

**Module:** `umol-ast/src/ast/validate/constraint/relational.rs`

**Kind:** additive (green)

**Dependencies:** `[dep: S3b, S3c]`

Evaluate every `RelationalConstraint` variant against the referenced dative bond, aromatic system,
multicenter bond, noncovalent bond, stereo atom, or stereo bond. Exact-set, contains, all, any,
acceptor/site, endpoint, and parallel-bond forms retain their declared ordered or unordered
semantics. Nested atom predicates dispatch through the same atom-constraint evaluation path,
including the shared ring projection when needed.

Focused tables cover every relational variant, bad references as operational errors, determined
truth and contradiction, vacuous nested predicates, and underdetermined nested predicates.

#### S3e — Molecule aggregate and connectivity constraints **Done**

**Module:** `umol-ast/src/ast/validate/constraint/molecule.rs`

**Kind:** additive (green)

**Dependencies:** `[dep: S3a, S3b]`

Evaluate `ChargeSum` and `BondOrderSum` over explicit subsets or the whole-molecule scope,
preserving underdetermination whenever a required value is non-literal. Validate the atom scope of
`UnpairedElectronCoupling`; its vacuous target is determined and every non-vacuous target remains
underdetermined pending the angular-momentum operation specified in doc 173. Evaluate `Connected`
with the supplied `ConnectedComponentsAlgorithm`: selected atoms must belong to one localized-bond
component, paths may pass through unselected atoms, and empty or singleton selections are
connected. References outside the molecule are operational errors rather than false predicates.

Focused tables cover subset and whole-molecule semantics, exact totals, partial values, empty
selections, disconnected selections, and selector parity where more than one implementation exists.

#### S3f — Subpattern constraints **Done**

**Module:** `umol-ast/src/ast/validate/constraint/molecule.rs`,
`umol-ast/src/ast/substructure.rs`, `umol-ast/src/ast/reaction.rs`,
`umol-graph/src/fingerprint/pattern.rs`, `umol-py/src/substructure.rs`,
`umol-py/src/reaction.rs`, `umol-py/src/fingerprint/config.rs`, and affected callers

**Kind:** breaking public matching-config migration (red → green)

**Dependencies:** `[dep: S3c, S3e]`

Add mandatory, default-free `SubstructureMatchConfig` at the AST layer, containing
`SubstructureMatchAlgorithm`, `SubgraphIsomorphismAlgorithm`, and
`RelevantCycleEnumerationAlgorithm`. Use it in `MoleculeAst::substructure_matches` and
`ReactionAst::apply`; higher-level graph and Python configs retain their documented defaults.
Pass the existing `ConstraintValidateConfig` to `MoleculeConstraintValidator::validate` instead of
four loose selectors.

Evaluate `SubPattern` using the supplied matching config, preserving each `SubPatternAnchor` form.
Thread the selected relevant-cycle algorithm into derived host ring constraints so subpattern
evaluation introduces no hidden graph-core algorithm choice. Use the existing substructure-match
result to decide the Boolean constraint; an early-exit search API is not part of this work unit.

Focused tables cover unanchored and anchored success, absence, overlays, ring-constrained patterns,
invalid anchors, and exact equivalence across every supported matching-selector combination.

#### S3g — Recursive constraint coordinator **Done**

**Module:** `umol-ast/src/ast/validate/constraint.rs`, `umol-ast/src/ast.rs`, and all AST-layer callers

**Kind:** breaking validator construction and call migration (red → green)

**Dependencies:** `[dep: S3b, S3c, S3d, S3e, S3f]`

Replace the unit stub with a configured `ConstraintValidator`. Evaluate inline entity containers and
the top-level `Constraints` store as one conjunction, and recursively fold `Constraint::And`, `Or`,
and `Not` without discarding contradiction or underdetermination. Lazily initialize ring,
connectivity, and substructure data only when a present non-vacuous leaf needs it. Stereo-entity
leaves compare references and stored values here; graph-symmetry truth remains solely in
`StereoConformanceValidator`.

Migrate every AST-layer caller to explicit `ConstraintValidateConfig`. Focused tables assert exact
mixed-category logical outcomes. Property tests cover permutation invariance of `And`/`Or`, double
negation, vacuous conjunction, and agreement between inline constraints and their equivalent
top-level entity leaves without using the coordinator to construct expected values.

**Stage exit:** every supported model-independent constraint has explicit determined,
underdetermined, contradictory, and operational-error semantics; the dative-ring stub has an explicit
outcome and doc-117 reference. `ConstraintValidator` is no longer a stub, and
`cargo test -p umol-ast` plus the feature-gated property suite pass.

### S4 — Higher-level validation and resolver preconditions

#### S4a — Composite validation configuration **Done**

**Module:** `umol-graph/src/ops/validate.rs` and affected Rust callers

**Kind:** breaking public config migration (red → green)

**Dependencies:** `[dep: S3g]`

Add `constraint: ConstraintValidateConfig` to `umol_graph::ValidateConfig`, construct the AST
validator from it, and preserve integrity → invariants → conformance ordering. Its documented
higher-level defaults are Relevant/Vismara, connected-components/BFS,
substructure/GraphAndOverlays, and subgraph-isomorphism/VF2-RDKit. Migrate every config literal and
caller in the same subitem. Composite validation must not run an algorithmic evaluator that the
input does not require.

Focused tables assert the complete default config, explicit selector propagation, exact wrapping of
each constraint contradiction/error category, and unchanged outcomes for molecules without
constraints.

#### S4b — Resolver constraint preconditions **Done**

**Module:** `umol-graph/src/ops/resolve/valence.rs`,
`umol-graph/src/ops/resolve/aromaticity.rs`,
`umol-graph/src/ops/resolve/multicenter.rs`, `umol-graph/src/ops/resolve/stereo.rs`,
`umol-graph/src/ops/resolve.rs`

**Kind:** behavioral correction (green)

**Dependencies:** `[dep: S2c, S3b]`

Before candidate selection, each direct resolver checks the incidence constraints that form its
preconditions: localized valence and dative pair counts, aromatic valence and aromatic-bond
incidence, multicenter valence, and tetrahedral/cis-trans stereo incidence. A determined mismatch is
an exact contradiction; a vacuous constraint does nothing; an underdetermined prerequisite stops the
resolution chain without emitting a partial plan. Keep the multidonor dative case as the documented
stub. Constraint derivation preserves non-literal aromatic and multicenter contributions as AST
values rather than assuming that relation fields are literal.

The composite resolver continues to materialize each successful stage before deriving the next and
applies the complete sequence transactionally. Focused tables assert no emitted edits on an
underdetermined stop, exact contradiction payloads, dependence on the materialized intermediate
state, and rollback of every earlier stage after a later precondition failure.

**Stage exit:** higher-level validation exposes documented defaults without hiding AST-layer
selectors, and every resolver rejects the model-independent inconsistencies relevant to its own
candidate derivation. `cargo test -p umol-graph` passes.

### S5 — Conformance completion

#### S5a — Aromatic constraint conformance **Done**

**Module:** `umol-graph/src/ops/validate/aromaticity.rs`, `umol-graph/src/ops/validate.rs`

**Kind:** behavioral correction (green)

**Dependencies:** `[dep: S2a]`

Keep `AromaticityConfig` as the validator's algorithm configuration; do not add a validation
inconsistency policy. Replace the current count-and-atom-set-only comparison with the S2a
classification. An asserted constraint without a matching relation, a stored relation rejected by
perception, or an explicit constraint/entity mismatch is `Solution::Contradictory`. Absent and
`Undetermined` constraints are ignored; a non-ground, non-vacuous contribution that cannot be
decided is `Solution::Underdetermined`. Extend `AromaticityValidatorContradiction` with exact
deterministic payloads and preserve the existing setup-error boundary.

Focused tables cover participant-set mismatch, per-atom contribution mismatch, localized-bond flag
mismatch, aromatic-valence mismatch, vacuous projections, model rejection of a stored system,
non-ground non-vacuous contributions, and a fully conformant existing system.

#### S5b — Stereo constraint conformance **Done**

**Module:** `umol-graph/src/ops/validate/stereo.rs`, `umol-graph/src/ops/validate.rs`

**Kind:** behavioral correction (green)

**Dependencies:** `[dep: S2b]`

Keep `StereoValidateConfig` limited to graph-symmetry algorithms and iteration limits; do not add an
inconsistency field. Before the existing relation and symmetry checks, use `StereoDerivation` to
require every non-vacuous asserted `#T` and `#C` site from `StereoPerception::derive` to have a
realizable, matching stereo relation. Absent and `Undetermined` projections are ignored. Add exact
atom/bond contradiction variants for an unrealizable, absent, or mismatched relation.

Focused tables cover absent and unrealizable tetrahedral/cis-trans relations, relation/assertion
mismatch, conformant existing relations, and preservation of every existing graph-symmetry
validation result.

#### S5c — Composite conformance outcomes **Done**

**Module:** `umol-graph/src/ops/validate.rs`

**Kind:** additive tests (green)

**Dependencies:** `[dep: S4a, S5a, S5b]`

Extend the composite validator tables with exact aromatic and stereo projection contradictions and
an underdetermined aromatic contribution. Assert that the contradiction is wrapped in the correct
`ValidatorContradiction` variant and that validation does not mutate the input.

**Stage exit:** standalone and composite validators reject constraint/entity inconsistencies
without resolver policy or mutation. `cargo test -p umol-graph` passes.

### S6 — Public ingestion acceptance

#### S6a — Rust ingestion propagation **Done**

**Module:** `umol-graph/src/ingest.rs`, `umol-graph/src/parse.rs`,
`umol-graph/tests/resolution/`

**Kind:** additive (green)

**Dependencies:** `[dep: S2a, S4b, S5a]`

Add molecule and reaction SMILES table cases proving that MDL furan, thiophene, and pyrrole return
the exact aromatic contradiction by default, while explicit `Keep` preserves the unmatched
constraints without adding a system. Add the bare-n `c1cccn1` case from doc 174 as an interim
regression: default ingestion returns the exact nitrogen aromatic-valence contradiction instead of
the old closed-shell C4H4N result. Retain positive MDL pyridine/benzene and Daylight
furan/thiophene/pyrrole
references. Update charge-sensitive resolution fixtures and snapshots to the localized
representation; do not normalize them through `DelocalizeCharge` in expected-value construction.

#### S6b — Python ingestion propagation

**Module:** `umol-py/tests/test_molecule.py`, `umol-py/tests/test_reaction.py`,
`umol-py/tests/test_workflow.py`

**Kind:** additive (green)

**Dependencies:** `[dep: S2a, S6a]`

Verify that configured molecule and reaction SMILES ingestion maps default aromatic rejection to the
existing `ContradictionError`, that explicit `Keep` returns the exact preserved representation, and
that the default output retains localized charge because Python has no implicit delocalization
option. No change to `umol-py/src/error.rs` is planned because the existing resolver-contradiction
mapping already has the correct public category. This stage does not add transformers or validators
to the Python surface.

**Stage exit:** Rust and Python ingestion demonstrate default rejection, explicit retention, and the
new localized charge representation. The feature-gated resolution conformance suite and installed
Python suite pass.

## Critical path and deferral

The final aromatic resolver path is `S0a + S0b → S1a → S2a → S2c → S4b → S6a → S6b`;
aromatic conformance is `S2a → S5a → S5c`. The stereo path is
`S0c → S1b → S2b → S2c → S4b → S5b → S5c`. Constraint validation is
`S3a → S3b + S3c → S3d + S3e → S3f → S3g → S4a`.

No stage in this work unit is deferrable. Rich model-rejection reasons, `LocalizeCharge`,
and Python transformation/validation methods remain separate proposed work. Transformations remain
explicit and outside the reaction primitive set.

Final verification is:

1. `cargo fmt --all`
2. `cargo test -p umol-ast`
3. `cargo test -p umol-ast --features proptest --test property --no-fail-fast`
4. `cargo test -p umol-graph`
5. `cargo test -p umol-graph --features conformance --test resolution`
6. `cargo test --workspace`
7. `cargo clippy --workspace --all-targets -- -D warnings`
8. With `umol-py/.venv` active and `python` confirmed as Python 3.13,
   `maturin develop` and `pytest -q umol-py/tests`
9. `git diff --check`
