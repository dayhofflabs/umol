# 177 — Nomenclature guide

Status: **Informational**
Date: 2026-07-31
Relates: [125](125-constraints-as-projections-2026-06-22.md),
[166](166-molecule-ops-2026-07-27.md),
[171](171-aromaticity-inconsistency-policy-2026-07-29.md),
[176](176-ast-naming-2026-07-31.md)

## Purpose

This is the living guide to terms coined or assigned a repository-specific meaning in umol. It is
normative for new public names and explanatory documentation once a term is recorded as settled. It
does not rename existing APIs by itself; a disagreement between the guide and existing code is a
separate migration task.

Prefer the established domain noun over a newly generalized synonym. Public names need not be
artificially parallel when the underlying constraints, entities, or operations have different names.
A shared name is useful only when the semantics and available operations are genuinely shared.

Doc 176 separately considers the unsettled `Ast` suffix and crate naming. This guide does not decide
that proposal.

## How to use this guide

- **Before naming something new**, check *Retired and discouraged* first. It is indexed by the wrong
  word, not the right one.
- **Every glossary entry is self-contained.** Nothing depends on the surrounding section, so landing
  mid-file by search is safe.
- **`Not:` lines name the confusable neighbours verbatim**, so searching for a near-miss term reaches
  the entry that corrects it.
- **`In code:` gives the identifiers**, which is the bridge from this prose to the API surface being
  edited.
- Headings name the **concept**. If a type is renamed, the `In code` line changes and the heading
  does not, so links and citations to this guide stay valid.

## Boundaries

Orientation, read once. The contrasts below are the ones most often collapsed; each participant has
its own glossary entry.

**Operations.** Three disjoint kinds, established in [doc 166](166-molecule-ops-2026-07-27.md):
*resolution* fills undetermined state using a chemistry model; *transformation* rewrites one resolved
representation into another; *validation* checks integrity, model-independent invariants, and
model-dependent conformance without mutation. Kekulization, aromatization, and charge delocalization
are transformations, not resolver behaviour, because they alter determined representation.

**Validation tiers.** *Integrity* (well-formed storage and shape, no model), *invariants*
(model-independent physics), *conformance* (accepted by a selected chemistry model). Distinct tiers,
run in that order.

**Derivation and policy.** *Perception* produces a policy-free *derivation*, which may carry
*inconsistencies*. A *policy* maps a classified inconsistency to a *recovery action*. Policies belong
to resolvers only: perception, chemistry models, and validators have none.

**Choices.** A *model* decides which result is chemically accepted. A *config* decides how an
operation is performed. A *policy* decides what happens after an inconsistency is established. An
*algorithm* selects one implementation of one algorithmic problem. These are four different things
and none is a synonym for another.

**Determinacy.** *Undetermined* is stored lattice state on a value. *Underdetermined* is an operation
outcome. The words differ by two letters and by kind; the glossary keeps them apart deliberately.

## Suffixes

The generative core of this guide. Each family gains members as work lands, so a new type's suffix
should be chosen here rather than by analogy with whichever neighbour was read last. Counts are as of
2026-07-31.

| Suffix | Means | Count | Crates |
| --- | --- | --- | --- |
| `*Algorithm` | selects one implementation of one algorithmic problem | 36 | graph-core (primitives), ast, graph |
| `*Config` | composite operational configuration | 30 | ast (ops), graph, io, py — never graph-core |
| `*Model` | semantic choices deciding chemical acceptance | 13 | graph |
| `*Policy` | maps a classified inconsistency to a recovery action | 11 | edn, graph, py |
| `*Kind` | unit-variant enum discriminating a family | 11 | ast, geometric, graph-core, msym, py |
| `*Constraint` | one assertable predicate over an entity | 6 | ast, py |
| `*Constraints` | the container holding an entity's constraints | — | currently spelled `*ConstraintsAst` |
| `*Key` | identifies a constraint slot within a container | 13 | ast, perm, py |
| `*View` / `*Views` / `*ViewMut` | borrowed accessor into a molecule; plural is the collection; `Mut` is the editing form | 37 / 15 / 16 | ast, io, py |
| `*Delta` | an encoded change belonging to a reaction side | 18 | ast, py |
| `*Update` | a field-level change applied to one entity | 10 | ast |
| `*Defaults` | values used where an input states nothing | 13 | ast, py |
| `*Overrides` | values that replace what an input stated | 10 | ast |
| `*Entry` | one row of a table-shaped registry or format | 41 | several |
| `*Contradiction` | a semantic rejection, the `Solution::Contradictory` payload | 18 | ast, graph |
| `*Mismatch` | two independently meaningful things that disagree | 3 | graph |
| `*Error` | operational or setup failure, the `Err` side | 58 | all |
| `*Validator` | performs one validation tier or one check within it | 12 | ast, graph |
| `*Resolver` | fills undetermined state under a chemistry model | 5 | graph |
| `*Ast` / `*Dsl` | internal representation / boundary surface | 77 / 38 | ast — **under review, doc 176** |

Families with no member yet, named in settled vocabulary and expected to gain one:
`*Failure` (doc 171), transformers (doc 166 — see *Transformer naming* below).

### Constraint singular and plural

`*Constraint` is one assertable predicate. `*Constraints` is the container holding an entity's
constraints. The plural is never a synonym for "several constraints" in a signature; it names the
store.

**Not:** using the plural for a `Vec<Constraint>` parameter.
**In code:** `AtomConstraintAst` against `AtomConstraintsAst`; likewise for every entity kind.
**Settled by:** 138.

### View, Views, ViewMut

A `*View` is a borrowed accessor into a molecule, not an owned value. `*Views` is the collection over
one entity kind. `*ViewMut` is the editing form. Sub-families exist: `*ConstraintsView`,
`*EditorView`, `*EditorViewMut`.

**Not:** the owned representation, which is the `*Ast` type; a view does not survive its molecule.
**In code:** `AtomView`, `AtomViews`, `AtomViewMut`, `AtomEditorViewMut`.
**Settled by:** 086.

### Delta and Update

A `*Delta` is an encoded change belonging to a reaction side. A `*Update` is a field-level change
applied to one entity. They are not interchangeable and exist per entity kind in both families.

**Not:** each other. `AtomDelta` is reaction-side change encoding; `AtomUpdate` is an entity edit.
**In code:** `AtomDelta`, `Deltas`, `ConstraintDelta` against `AtomUpdate`, `AtomFieldChange`.
**Settled by:** 134.

### Defaults and Overrides

`*Defaults` supplies values where an input states nothing. `*Overrides` replaces values the input did
state.

**Not:** each other. The difference is whether the input spoke.
**In code:** `MoleculeDefaults`, `AtomDefaults` against `MoleculeOverrides`, `AtomOverrides`.
**Settled by:** to fill during the sweep.

### Transformer naming — unsettled

Three patterns are in use for the same family, and doc 166 adds two more members:

- agent noun: `Kekulizer`, `Aromatizer`, each having grown `*Error` and `*Config`;
- verb phrase: `DelocalizeCharge`;
- target phrase: `ToExplicitHydrogens`, `ToImplicitHydrogens` (doc 166, not yet written).

The house pattern for a thing that performs an operation is the agent noun — `Resolver`, `Validator`,
`Kekulizer`. The counter-argument is that `Kekulizer` and `Aromatizer` are configured engines with
their own error and config types, while `DelocalizeCharge` and the hydrogen transforms are
parameterless operations, which may justify two conventions rather than one.

This needs a decision before doc 166 lands its three transformers, not after.

## Retired and discouraged

Indexed by the word not to use. Add a row whenever an entry's `Not:` line forbids a specific
spelling.

| Do not write | Write instead | Because |
| --- | --- | --- |
| `projection` for a stored constraint or the incidence category | `constraint`, `incidence constraint` | `projection` names an actual mapping between representations |
| `predicate` or `representation` for the stored object | `constraint` | the repository term for a possibly non-ground assertion |
| `EntityConstraint` for the non-ring subset | `incidence constraint` | entity constraints include ring constraints |
| `relational constraint` for an incidence constraint | `incidence constraint` | `RelationalConstraint` already means a molecule-scope, reference-bearing constraint |
| `undetermined` for an operation outcome | `underdetermined` | stored state versus outcome |
| `underdetermined` for stored lattice state | `undetermined` | as above |
| `reset` for removing or replacing an entity | `remove`, `replace` | reset clears a constraint to its undetermined form |
| `config` for semantic acceptance choices | `model` | config is operational |
| `policy` for chemical acceptance | `model` | policy acts after acceptance is established |

## Glossary

Alphabetical by concept. Entry shape: definition, optional detail, then `Not:`, `In code:`,
`Settled by:`.

### Algorithm

An **algorithm** enum selects a concrete implementation of one algorithmic problem. `*Algorithm` is
the suffix for algorithm selectors: `umol-graph-core` defines the graph-algorithm primitives, and
higher layers follow the same suffix.

**Not:** config (which may *contain* an algorithm selection), model, policy.
**In code:** `*Algorithm` enums in `umol-graph-core`, e.g. `CommonSubgraphEnumerationAlgorithm`,
`RelevantCycleEnumerationAlgorithm`; also higher-layer selectors such as `SubstructureMatchAlgorithm`.
**Settled by:** 171.

### Application

**Application** executes a complete edit plan transactionally and publishes the result only when
every edit succeeds.

**Not:** plan (which is derived without mutating), transformation.
**In code:** `apply`, `apply_at`.
**Settled by:** 166.

### Combine

**Combine** forms the disjoint union of two structures, yielding one structure with two components.

Named `join` originally; renamed to avoid collision with the lattice join. Implemented and unlikely
to gain analogues, so it is recorded rather than generative.

**Not:** *join*, which is the lattice least upper bound.
**In code:** `combine`, `combine_all`, `combine_from`.
**Settled by:** 151.

### Config

A **config** contains operational choices controlling how an operation is performed, including
algorithm selection and iteration limits. `*Config` is the suffix for composite configuration: in
`umol-ast` it belongs to ops rather than to individual methods, and it is used throughout
`umol-graph`, `umol-io`, and `umol-py`. `umol-graph-core` defines no configs; its operations take
algorithm selectors directly.

**Not:** model, policy.
**In code:** `*Config`, e.g. `ResolveConfig`, `ValidateConfig`, `SubstructureSearchConfig`,
`ConstraintValidateConfig`, `SmilesIoConfig`.
**Settled by:** 171.

### Conformance

**Conformance** asks whether a structure is accepted by a selected chemistry model. Conformance
validators may reuse policy-free derivations, but they convert every determined inconsistency into a
contradiction.

**Not:** integrity, invariant. All three are validation tiers and are not interchangeable.
**In code:** `validate_conformance`, `*ConformanceValidator`.
**Settled by:** 171.

### Constraint

A **constraint** is a possibly non-ground assertion represented by the constraint AST types. This is
the public repository term.

**Not:** `projection`, `predicate`, or `representation`, when naming the stored object or an
operation over it.
**In code:** `AtomConstraintAst`, `AtomConstraintsAst`, and the per-entity equivalents.
**Settled by:** 125, 171.

### Contradiction

A **contradiction** is a semantic rejection represented by `Solution::Contradictory`. Validators have
no recovery policy: every determined failure or mismatch in their scope becomes a contradiction.

**Not:** error (operational, outside `Solution`), failure or inconsistency (policy-free
classifications that a resolver may still act on).
**In code:** `Solution::Contradictory`, `*Contradiction`.
**Settled by:** 171.

### Derivation

A **derivation** is the policy-free result of perception, including candidates and exact
inconsistencies.

**Not:** resolution (which applies policy), validation.
**In code:** —
**Settled by:** 171.

### Determined

**Determined** is an operation outcome: the pass produced a fully resolved result.

**Not:** *ground*, which is a property of stored state rather than of an outcome. An operation may
return `Determined` with a payload that is not ground if the operation's own contract is satisfied.
**In code:** `Solution::Determined`, `is_determined`, `into_determined`.
**Settled by:** to fill during the sweep.

### Entity

An **entity** is one of the eight kinds represented by `Entity`: atom, localized bond, dative bond,
aromatic system, multicenter bond, noncovalent bond, stereo atom, or stereo bond. Use the concrete
entity name in diagnostics and action fields when it matters which kind is affected.

**Not:** *kind*, which names what distinguishes one entity family from another rather than an
instance. Not *overlay* either: entity covers all eight kinds, overlay covers the six that are not
topology.
**In code:** `Entity`, `EntityKind`.
**Settled by:** 134, 086. Origin of the approach: 079.

### Entity constraint

An **entity constraint** is the broad category of constraints stored on or addressed to an entity. It
includes ring constraints.

**Not:** a name for the non-ring subset — use *incidence constraint* for that.
**In code:** `EntityConstraint` must not be narrowed to the non-ring subset.
**Settled by:** 171.

### Error

An **error** is an operational or setup failure outside the semantic `Solution`, such as a failed
transaction or unavailable model parameters.

**Not:** contradiction, failure, inconsistency — all of which are semantic.
**In code:** the `Err` side of `Result<Solution<_, _>, _>`. Every module error type implements
`umol_utils::UmolError`, which supplies `as_any` for downcasting; `Box<dyn UmolError>` is the
cross-module boundary form and `?` promotes into it.
**Settled by:** 171.

### Failure

A **failure** is a policy-free derivation classification. A constraint failure means that a
non-vacuous constraint cannot produce a valid entity under the selected topology and model. An entity
failure means that a structurally readable entity is not realizable under them.

**Not:** error (operational, not semantic); contradiction (a failure is not yet one — a resolver may
be configured to retain or remove the affected input).
**In code:** —
**Settled by:** 171.

### Ground term

A **ground term** is a structure in which every inherent field holds a definite value. Groundness is
structural: it says the lattice is resolved to a bottom element, not that the structure satisfies
chemistry invariants or that its entities are mutually consistent.

**Not:** valid, or chemically admissible. Doc 173 separates structural groundness from both.
**In code:** `Lattice::is_ground`, `Ground<T>` (planned, doc 175).
**Settled by:** 173. Whitepaper glossary carries the chemist-facing form.

### Incidence constraint

An **incidence constraint** is an entity-local constraint whose value is derived from the entity's
fields and directly incident localized bonds or overlay relations, without running a separate graph
algorithm. The category includes:

- degree and valence aggregates over incident localized bonds;
- donated- and accepted-pair aggregates over incident dative bonds;
- aromatic-valence and aromatic-bond values from incident aromatic systems;
- multicenter-valence values from incident multicenter bonds;
- the corresponding directly incident stereo relations.

`Incidence` names the structural relationship, not the value type. An incidence-derived value may be
Boolean, a count, or a weighted sum.

**Not:** *relational constraint* (which has a narrower existing meaning); *projection*; a synonym for
the whole entity-constraint category.
**In code:** `IncidenceConstraintValidator`, `IncidenceConstraintContradiction`; the established
`incident*` methods expose the same relationship.
**Settled by:** 171.

### Inconsistency

An **inconsistency** is the policy-free umbrella classification for a constraint failure, entity
failure, or entity/constraint mismatch found during derivation. An inconsistency does not by itself
select an authority or recovery action.

**Not:** contradiction. A resolver policy may retain the inconsistency, remove one or both inputs
where sound, replace the entity, or report a contradiction.
**In code:** `AromaticityInconsistency`, `StereoInconsistency`, which identify the exact constraint
site and entity involved.
**Settled by:** 171.

### Inherent field

An **inherent field** identifies an entity and is always present on it, though its value may be
undetermined.

**Not:** constraint, which is assertable without contributing to identity and may be absent.
**In code:** the non-constraint fields of each `*Ast` entity type.
**Settled by:** 104, 173.

### Integrity

**Integrity** concerns well-formed storage, entity shape, references, and model-independent
constraint consistency.

**Not:** invariant, conformance. All three are validation tiers and are not interchangeable.
**In code:** `validate_integrity`, `EntityStructureValidator`, `ConstraintValidator`.
**Settled by:** 171.

### Invariant

An **invariant** is a model-independent physical or mathematical condition over an otherwise
well-formed structure.

**Not:** integrity, conformance. All three are validation tiers and are not interchangeable.
**In code:** `validate_invariants`, `ValenceInvariantsValidator`, `SpinInvariantsValidator`.
**Settled by:** 171.

### Lattice

**Lattice** is the internal term for the partial order on attribute values under which matching,
meet, and join are defined. The whitepaper calls it the *attribute lattice*; "attribute" is a
clarification added for chemist readers and is not repository terminology.

**Not:** *attribute lattice* in code or internal documentation.
**In code:** `Lattice`, `meet`, `join`, `matches`, `is_compatible`.
**Settled by:** 113.

### Mismatch

A **mismatch** means that a constraint and entity, or a constraint and its derived value, are each
independently meaningful but disagree.

**Not:** failure (where one side is not realizable at all).
**In code:** concrete names in diagnostic variants and config fields, such as
`AromaticValenceMismatch`, `CisTransStereoMismatch`.
**Settled by:** 171.

### Model

A **model** contains semantic choices defining which result is chemically accepted.

**Not:** config, policy, algorithm.
**In code:** `ChemistryModel`, `ValenceModel`, `AromaticityModel`, `StereoModel`, `RingModel`.
**Settled by:** 171.

### Molecular structure

A **molecular structure** is a molecular topology together with its overlay entities.

**Not:** molecular topology, which excludes them.
**In code:** —
**Settled by:** whitepaper glossary; no settled repository term.

### Molecular topology

The **molecular topology** is the attributed undirected simple graph of atoms and localized bonds,
carrying no aromatic, stereo, or coordination information.

**Not:** molecular structure, which includes the overlay entities.
**In code:** —
**Settled by:** whitepaper glossary; no settled repository term.

### Overlay

An **overlay** is one of the six entity kinds that are not molecular topology: dative bond, aromatic
system, multicenter bond, noncovalent bond, stereo atom, stereo bond. Atoms and localized bonds are
the topology and are not overlays.

**Scope conflict, unresolved.** `IncidenceNodeSelection::OVERLAYS` covers only four of the six —
dative, aromatic, multicenter, noncovalent — with `STEREO` as a separate flag, so that
`constitution()` is atoms, bonds and those four while `full()` adds stereo. The narrower sense is
chemically motivated (constitution against configuration) but it makes the same word mean four kinds
in one place and six everywhere else. Either the flag or the term should change.

**Not:** *relation* or *hyperedge*, which are whitepaper framings for the same thing and appear in
source comments descriptively; overlay is the repository term. Not *entity*, which is the umbrella
over all eight kinds.
**In code:** `GraphAndOverlays`, `verify_overlays`, `RemovedOverlays`,
`IncidenceNodeSelection::OVERLAYS` (narrower, see above).
**Settled by:** 134. Whitepaper §5 carries the chemist-facing definition.

### Participant

A **participant** is an entity referenced by an overlay entity — the atoms of an aromatic system, the
donors and acceptor of a dative bond, the site and ligands of a stereo element.

**Not:** member, constituent, or argument.
**In code:** `ParticipantPosition`.
**Settled by:** 134.

### Perception

**Perception** applies a selected chemistry model and algorithms to identify or assess structural
entities. It produces a policy-free derivation.

**Not:** resolution (which applies policy and edits), validation.
**In code:** `AromaticityPerception`.
**Settled by:** 171.

### Plan

A **plan** is the complete edit sequence derived without mutating the source object.

**Not:** application (which executes a plan).
**In code:** `plan`.
**Settled by:** 166.

### Policy

A **policy** is operational configuration mapping a classified inconsistency to a recovery action.
Policies belong to resolvers, not perception, derivation, chemistry models, or validators. A policy
enum is shared when the available action set is identical; the config field supplies the concrete
constraint or entity context.

**Not:** model. A model determines chemical acceptance; a policy determines what an operation does
after acceptance or inconsistency has been established.
**In code:** `AromaticityInconsistencyPolicy`, `StereoInconsistencyPolicy`.
**Settled by:** 171.

### Projection

Use **projection** for an actual mapping from one representation or indexed relation to another, as
in pullback projections, reaction-side projections, or cycle projection from a subdivision graph. It
may be used descriptively when explaining that an entity relation induces participant values.

**Not:** the public name of stored constraints, or of the incidence-constraint category.
**In code:** —
**Settled by:** 125.

### Recovery action

A **recovery action** is a policy variant naming what an operation does about a classified
inconsistency:

- `Keep` retains the affected inputs;
- `Remove` removes the target identified by the config field;
- `RemoveConstraint` retains the independently valid entity;
- `ReplaceEntity` atomically replaces the entity with the valid constraint-derived result;
- `RemoveBoth` removes both members of the pair identified by the mismatch policy;
- `Error` converts the classified inconsistency into a contradiction.

**Not:** separate enums created solely to repeat identical variants with longer target-specific
spellings. Policy enums may be shared when they admit exactly the same action set; the config field
supplies the context.
**In code:** the variants above.
**Settled by:** 171.

### Refinement

**Refinement** names three unrelated operations. Always qualify it.

- **Lattice refinement** — the order relation on attribute values; `b` refines `a` when `a ∧ b = b`.
  Doc 113.
- **Colour refinement** — the automorphism and canonical-labelling procedure. `RefinementAlgorithm`,
  `umol-graph-core/src/algorithms/refinement.rs`.
- **Circular refinement** — the fingerprint iteration over atom environments.
  `CircularRefinementAlgorithm`, `CircularRefinementHash`. Docs 126, 154.

**Not:** any of the three used unqualified where another could be meant.
**In code:** `RefinementAlgorithm`, `CircularRefinementAlgorithm`, `Lattice::matches`.
**Settled by:** 113, 126, 154.

### Relational constraint

A **relational constraint** is a molecule-scope, reference-bearing constraint relating an overlay
entity to atoms, bonds, roles, or predicates.

**Not:** a synonym for *incidence constraint*, merely because the derived value depends on a
relation.
**In code:** `RelationalConstraint`.
**Settled by:** 171.

### Reset

**Reset** clears a constraint by setting it to its undetermined form.

**Not:** a general synonym for removing an entity or replacing determined structural information.
**In code:** `reset_aromatic_valence`, `reset_stereo_constraints`.
**Settled by:** 171.

### Resolution

**Resolution** applies configured recovery policy, constructs one atomic edit plan, and may apply
that plan transactionally.

**Not:** transformation (which rewrites determined representation), validation (which does not
mutate), perception (which is policy-free).
**In code:** `Resolver`, `resolve`, `ResolveConfig`.
**Settled by:** 166, 171.

### Ring constraint

A **ring constraint** is a constraint whose value requires ring enumeration. Ring membership, ring
degree, and ring valence belong here even though they are stored in atom, bond, or dative-bond
constraint containers.

Ring constraints are separate from incidence constraints because their evaluation has an algorithmic
dependency. Their fixed molecular semantics are the Relevant ring set through size 22; the selected
relevant-cycle enumeration algorithm is operational configuration.

**Not:** incidence constraint.
**In code:** `RingConstraintValidator`, `RingConstraintContradiction`.
**Settled by:** 171.

### Solution

**`Solution<T, C>`** is the three-valued outcome of an engine pass: `Determined`, `Underdetermined`,
or `Contradictory(C)` with a typed diagnostic payload. Setup and parameter failures travel separately
in `Result` and never collapse into it, so every umol engine returns `Result<Solution<_, _>, _>`.

An operation must decide which outcomes it treats as success, and the two reasonable answers differ
precisely on `Underdetermined`: a validator accepts it (only `Contradictory` fails), a transformer
does not (it cannot rewrite a representation on partial information). State the choice when adding an
operation; it is the outcome most easily overlooked. `into_observation` and `into_decisive` encode
the two, but with 4 and 1 call sites they are conveniences rather than established vocabulary.

**Not:** `Result`. `Solution` carries the semantic verdict; `Result` carries operational success.
Both appear in one signature and mean different things.
**In code:** `Solution`, `umol_utils::solution`.
**Settled by:** to fill during the sweep.

### Split

**Split** decomposes a structure into its connected components.

Implemented and unlikely to gain analogues; recorded rather than generative.

**Not:** a partition by any other criterion.
**In code:** `split`.
**Settled by:** 086.

### Transformation

A **transformation** explicitly rewrites one valid representation into another. Kekulization,
aromatization, and charge delocalization are transformations because they alter determined
representation rather than fill undetermined state.

**Not:** a resolver policy; not resolution.
**In code:** `umol-graph/src/ops/transform`.
**Settled by:** 166.

### Underdetermined

**Underdetermined** is an operation outcome: the available input does not justify a determined answer
or a contradiction.

**Not:** *undetermined*, which is stored state. Do not use the two interchangeably.
**In code:** `Solution::Underdetermined`.
**Settled by:** 171.

### Undetermined

**Undetermined** is stored lattice state: an AST value asserts no concrete value at that position. An
absent or undetermined constraint is vacuous and does not assert that an entity is missing.

**Not:** *underdetermined*, which is an operation outcome. Do not use the two interchangeably.
**In code:** `ValueAst::Undetermined` and the per-type equivalents.
**Settled by:** 171.

### Validation

**Validation** reports semantic contradictions without repairing or selecting an authoritative
representation.

**Not:** resolution (which repairs), transformation. Integrity, invariants, and conformance are its
three tiers.
**In code:** `Validator`, `validate`.
**Settled by:** 166, 171.

## Maintaining this guide

Add a term when a discussion assigns it a durable meaning that future code or documentation must
preserve. Every entry carries four slots:

- the definition, in one sentence where possible;
- `Not:` — the nearby terms it must not be used to mean, spelled out so a search for the wrong word
  reaches this entry;
- `In code:` — the public type or method names implied by the decision, or `—` if none yet;
- `Settled by:` — the discussion document owning the decision and any implementation work.

Whenever a `Not:` line forbids a specific spelling, add the corresponding row to *Retired and
discouraged*.

Do not turn the guide into an inventory of every exported type. Ordinary Rust or chemistry terms need
entries only when umol narrows, distinguishes, or coordinates their meanings.

Keep the file loadable in one piece — under roughly a thousand lines. If it outgrows that, split by
domain deliberately rather than letting it sprawl past the point where a reader or an agent sees all
of it.

### Provenance still to fill

`Entity` is credited to the sweep rather than to a document. `Derivation`, `Failure`, and
`Projection` have no `In code` names yet. Both gaps are expected to close when the terminology sweep
over the discussion corpus runs.
