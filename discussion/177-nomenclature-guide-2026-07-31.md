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

This initial version records the constraint and resolver vocabulary settled in docs 166 and 171. A
separate later sweep should extract other durable terminology from recent discussion documents; that
sweep is not part of the present work.

## Glossary

Terms are ordered alphabetically by their primary public name.

### Algorithm

An **algorithm** enum selects a concrete implementation of one algorithmic problem.

Algorithmic transparency does not require flattening every subsidiary selector into every high-level
method signature:

- graph-core and focused lower-level operations receive the relevant algorithm enum directly;
- composite operations collect subsidiary choices in operation-specific, nestable config objects;
- higher-level APIs may provide documented defaults;
- an algorithm-dependent evaluator remains separate from evaluators that do not need it;
- unused algorithmic subsystems are not executed merely because their config is present.

### Application

**Application** executes a complete edit plan transactionally and publishes the result only when
every edit succeeds.

### Config

A **config** contains operational choices controlling how an operation is performed, including
algorithm selection and iteration limits. Config is not a synonym for model.

### Conformance

**Conformance** asks whether a structure is accepted by a selected chemistry model. It is a
validation tier, distinct from integrity and invariants. Conformance validators may reuse policy-free
derivations, but they convert every determined inconsistency into a contradiction.

### Constraint

A **constraint** is a possibly non-ground assertion represented by the constraint AST types. This is
the public repository term. Do not replace it with `projection`, `predicate`, or `representation`
when naming the stored object or an operation over it.

### Contradiction

A **contradiction** is a semantic rejection represented by `Solution::Contradictory`. Validators have
no recovery policy: every determined failure or mismatch in their scope becomes a contradiction.

### Derivation

A **derivation** is the policy-free result of perception, including candidates and exact
inconsistencies.

### Entity

An **entity** is one of the eight kinds represented by `Entity`: atom, localized bond, dative bond,
aromatic system, multicenter bond, noncovalent bond, stereo atom, or stereo bond. Use the concrete
entity name in diagnostics and action fields when it matters which kind is affected.

### Entity constraint

An **entity constraint** is the broad category of constraints stored on or addressed to an entity.
It includes ring constraints; therefore `EntityConstraint` must not be used as the name of the
non-ring subset.

### Error

An **error** is an operational or setup failure outside the semantic `Solution`, such as a failed
transaction or unavailable model parameters. Error types occupy the `Err` side of
`Result<Solution<_, _>, _>`.

### Failure

A **failure** is a policy-free derivation classification. A constraint failure means that a
non-vacuous constraint cannot produce a valid entity under the selected topology and model. An entity
failure means that a structurally readable entity is not realizable under them.

Failure is not `Error`: it is a semantic outcome, not an operational failure. It is also not itself a
`Contradiction`; a resolver may be configured to retain or remove the affected input.

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
Boolean, a count, or a weighted sum. The established `incident*` methods expose the same underlying
relationship.

The intended validator names are `IncidenceConstraintValidator` and
`IncidenceConstraintContradiction`.

### Inconsistency

An **inconsistency** is the policy-free umbrella classification for a constraint failure, entity
failure, or entity/constraint mismatch found during derivation. `AromaticityInconsistency` and
`StereoInconsistency` identify the exact constraint site and entity involved.

An inconsistency does not by itself select an authority or recovery action. A resolver policy may
retain it, remove one or both inputs where sound, replace the entity, or report a contradiction.

### Integrity

**Integrity** concerns well-formed storage, entity shape, references, and model-independent
constraint consistency. It is a validation tier, distinct from invariants and conformance.

### Invariant

An **invariant** is a model-independent physical or mathematical condition over an otherwise
well-formed structure. It is a validation tier, distinct from integrity and conformance.

### Mismatch

A **mismatch** means that a constraint and entity, or a constraint and its derived value, are each
independently meaningful but disagree. Use the concrete constraint and entity names in diagnostic
variants and config fields, such as `AromaticValenceMismatch` or `CisTransStereoMismatch`.

### Model

A **model** contains semantic choices defining which result is chemically accepted. Model is not a
synonym for config or policy.

### Perception

**Perception** applies a selected chemistry model and algorithms to identify or assess structural
entities. It produces a policy-free derivation.

### Plan

A **plan** is the complete edit sequence derived without mutating the source object.

### Policy

A **policy** is operational configuration mapping a classified inconsistency to a recovery action.
Policies belong to resolvers, not perception, derivation, chemistry models, or validators. A policy
enum is shared when the available action set is identical; the config field supplies the concrete
constraint or entity context.

Policy is not a synonym for model. A model determines chemical acceptance; a policy determines what
an operation does after acceptance or inconsistency has been established.

### Projection

Use **projection** for an actual mapping from one representation or indexed relation to another, as
in pullback projections, reaction-side projections, or cycle projection from a subdivision graph.
It may be used descriptively when explaining that an entity relation induces participant values, but
it is not the public name of stored constraints or of the incidence-constraint category.

### Recovery action

Use action names for policy variants:

- `Keep` retains the affected inputs;
- `Remove` removes the target identified by the config field;
- `RemoveConstraint` retains the independently valid entity;
- `ReplaceEntity` atomically replaces the entity with the valid constraint-derived result;
- `RemoveBoth` removes both members of the pair identified by the mismatch policy;
- `Error` converts the classified inconsistency into a contradiction.

Policy enums may be shared when they admit exactly the same action set. The config field supplies the
concrete constraint or entity context; separate enums must not be created solely to repeat identical
variants with longer target-specific spellings.

### Relational constraint

`RelationalConstraint` has a narrow existing meaning: a molecule-scope, reference-bearing constraint
relating an overlay entity to atoms, bonds, roles, or predicates. Do not use `relational constraint`
as a synonym for incidence constraint merely because the derived value depends on a relation.

### Reset

**Reset** clears a constraint by setting it to its undetermined form. It is not a general synonym for
removing an entity or replacing determined structural information.

### Resolution

**Resolution** applies configured recovery policy, constructs one atomic edit plan, and may apply
that plan transactionally.

### Ring constraint

A **ring constraint** is a constraint whose value requires ring enumeration. Ring membership, ring
degree, and ring valence belong here even though they are stored in atom, bond, or dative-bond
constraint containers.

Ring constraints are separate from incidence constraints because their evaluation has an algorithmic
dependency. Their fixed molecular semantics are the Relevant ring set through size 22; the selected
relevant-cycle enumeration algorithm is operational configuration.

The intended validator names are `RingConstraintValidator` and `RingConstraintContradiction`.
Dative-bond ring membership remains undefined until the ring topology for dative versus
coordination/haptic entities is settled.

### Transformation

A **transformation** explicitly rewrites one valid representation into another. It is not a resolver
policy. Kekulization, aromatization, and charge delocalization are transformations because they alter
determined representation rather than fill undetermined state.

### Underdetermined

**Underdetermined** is an operation outcome: the available input does not justify a determined answer
or a contradiction. It is represented by `Solution::Underdetermined`.

### Undetermined

**Undetermined** is stored lattice state: an AST value asserts no concrete value at that position. An
absent or undetermined constraint is vacuous and does not assert that an entity is missing.

Do not use `undetermined` and `underdetermined` interchangeably.

### Validation

**Validation** reports semantic contradictions without repairing or selecting an authoritative
representation. Integrity, invariants, and conformance are distinct validation tiers.

## Coordination of constraint validation

Incidence-constraint validation is algorithm-free. Ring-constraint validation receives a
relevant-cycle algorithm. A complete constraint validator coordinates them through structured
configuration rather than a growing flat argument list. The same organization extends to other
algorithm-dependent constraint categories as they are implemented.

The policy and validation vocabulary above is settled by
[doc 171](171-aromaticity-inconsistency-policy-2026-07-29.md).

## Maintaining this guide

Add a term when a discussion assigns it a durable meaning that future code or documentation must
preserve. Record:

- the precise scope of the term;
- nearby terms it must not be used to mean;
- the public type or method names implied by the decision;
- the discussion document owning any implementation work.

Do not turn the guide into an inventory of every exported type. Ordinary Rust or chemistry terms
need entries only when umol narrows, distinguishes, or coordinates their meanings.
