# Data type contracts

## Purpose

This is a normative developer guide for deciding which properties belong in constructors,
converters, validators, and transformations. It applies primarily to aggregate model types such as
`MoleculeAst`, `ReactionAst`, and `ReactionSpanAst`. Small value types may enforce stronger
invariants when those invariants define the value represented by the type.

The central rule is that an aggregate constructor establishes that its input can be represented; it
does not establish every useful property of the resulting value. Semantic properties are lazy and
must be requested through named operations.

An invariant is therefore enforced by the first operation that requires it, not by the earliest
operation capable of checking it. A caller may request an explicit validator earlier. If it does
not, a later conversion or operation reports the failure when that invariant becomes a precondition
of producing its result.

This distinction is already applied to DPO validation: a dangling reaction is constructible but
DPO-invalid until checked. This guide generalizes that decision and makes it available independently
of the reaction implementation history.

## Operation taxonomy

| Operation | Establishes | Does not implicitly do |
| --- | --- | --- |
| Construction | Representation and referential integrity | Semantic validation, canonicalization, resolution, repair |
| Conversion | Faithful representation in the target type | Semantic validation, normalization, silent loss or repair |
| Validation | A named semantic property | Mutation or repair |
| Transformation | An explicitly named change of representation or state | Pretend that the input was already valid or canonical |

### Construction

Construction enforces the minimum invariants needed for the value to be a coherent instance of its
type. For aggregate ASTs, these are representation invariants:

- stored references name entities in the owning namespace;
- participant, site, ligand, anchor, and constraint references are resolvable;
- correspondence pairs are in range and form a partial bijection;
- parallel collections have the shape required by the representation;
- a representation variant contains the data required to interpret that variant.

Construction does not establish model-independent semantics merely because they can be checked
without a chemistry model. In particular, construction does not imply:

- DPO dangling-freedom;
- chemical validity or conformance to a model;
- agreement between independently meaningful entities and constraints;
- groundness;
- canonical form;
- satisfaction of an operation's preconditions.

An asserted constructor such as `from_entries` and a checked constructor such as
`try_from_entries` establish the same invariant. They differ only in how a violated construction
contract is reported. The asserted form is for producers that establish the invariant by
construction; the checked form is for untrusted or independently assembled input.

### Conversion

A conversion preserves every source state representable by its target, including semantically
invalid states. It fails only when the target representation cannot structurally encode the source
or when the source contains contradictory information that cannot be represented as one target
value.

A conversion must not silently make its input acceptable. Dropping dangling entities, removing
constraints, selecting a semantic interpretation, canonicalizing, or resolving undetermined state
are transformations rather than conversion mechanics. If such behavior is useful, expose it as a
separately named operation.

The error belongs to the failed boundary:

- failure to assemble the target's entries uses the target entry error;
- contradictory source information uses the relevant contradiction or a conversion error that
  preserves that cause;
- failure of an optional semantic precondition belongs to the operation requesting that property,
  not to the target constructor.

### Validation

Validation checks one named property without mutation. It may traverse the entire value and may be
more expensive than construction. That cost is appropriate because the caller explicitly requested
the property.

Validation remains separate even when a property is model-independent. DPO dangling-freedom is a
semantic invariant of a rule rather than an invariant required to store the rule. Canonicality is
similarly useful but not required for representation.

### Transformation

Canonicalization, resolution, repair, stripping, cascading removal, and closure under a rewriting
semantics are transformations. Their names and return types must expose the change. A constructor or
ordinary conversion must not perform one as an incidental implementation step.

## Provenance and contextual validity

Whether a consumer must re-establish a contextual property depends on what kind of value it
receives.

An **open data carrier** can be constructed independently of the objects with which it will later be
used. Its intrinsic constructor checks only its own representation, but a contextual consumer checks
the properties required to combine it with those objects. `Correspondence` is such a carrier: atom
pairs may come from SMIRKS or another mapped external format, and Rust and Python both permit direct
construction. A valid partial bijection is not necessarily a correspondence over the particular two
molecules supplied to a later operation. `MoleculeCorrespondence` is likewise not bound to molecule
instances by its ids.

An **operation-issued value** may instead be provenance-bound. A `Transaction` is issued by applying
edits and records how to undo that particular successful application. It has no public constructor
for independently asserting that provenance. Replacing it with another transaction or mutating the
object independently violates the operation contract; rollback must not panic, but it does not owe
correct restoration for the compromised pairing.

The resulting rules are:

- do not add validation merely to defend against swapping an opaque, provenance-bound result into a
  different operation history;
- do validate an open carrier when a public operation legitimately accepts it alongside independent
  objects and needs contextual agreement to produce a correct result;
- an internal producer-consumer path may assert the contextual property established by its producer;
- simple ids and open correspondences remain deliberately unbound to object identity; this avoids
  heavier identity infrastructure but makes contextual checks the responsibility of operations that
  promise a relationship between independently supplied objects.

For molecule correspondences, an atom correspondence read from an external format is a normal input
to induction over a supplied molecule pair. Count agreement and any structural uniqueness required
to derive the remaining entity families are therefore operation preconditions, not defenses against
deliberate tampering. A full correspondence produced for the same molecule pair by a conforming
operation satisfies those checks by provenance. Reusing that result with the same unchanged pair
cannot newly produce a contextual mismatch. Supplying an atom correspondence and molecule pair
independently is different: this is the ordinary bridge from mapped formats such as SMIRKS, and
induction may report incompatible carrier sizes or non-unique entity incidence. Supplying a full
correspondence for another molecule pair is likewise an ordinary public-input case whenever the
consuming API accepts the correspondence and pair independently; it is not treated as tampering
with an opaque operation result.

### Containing fallibility

Fallibility is not propagated merely because an implementation calls another fallible operation.
It belongs on the first public operation whose promised result cannot be produced without the
property:

- use `Option` when there is one ordinary absence condition and the cause carries no useful detail;
- use `Result` when the caller can distinguish or act on failure causes;
- preserve an existing operation-specific error surface when an internal implementation route
  changes;
- assert only on internal paths whose producer establishes the required property.

For example, converting a partial correspondence to a dense remapping may return `None` because no
total-left mapping exists. This can occur for a correspondence correctly produced for its molecule
pair: partiality is part of the correspondence model, whereas a remapping is total on its source.
It does not make correspondence construction, composition, reversal, or unrelated consumers
fallible. Exact error taxonomy remains subject to the repository-wide error review; the
construction/validation boundary does not require introducing a new error type for each method.

## Reaction-span application

`ReactionSpanAst` stores a union-frame structure. Referential integrity is therefore evaluated in
that union namespace. A bond, overlay, stereo entity, or constraint may reference a union entity
that is absent from one side without making the span structurally malformed.

For example, an `Unchanged` bond incident to a `Removed` atom is:

- representable in the union-frame span;
- DPO-invalid;
- not projectable to a referentially intact right-hand `MoleculeAst` without changing its meaning.

For a span materialized from `ReactionAst`, the lhs is already a referentially intact
`MoleculeAst`, while the deltas may be DPO-invalid. `ReactionAst::to_reaction_span` therefore
succeeds: the union-frame span can represent the rule exactly. `ReactionSpanAst::rhs` is the first
operation in this path that requires the rhs to be a referentially intact molecule, so it reports
the unavailable reference. Calling `DpoValidator` before projection reports the same underlying
problem earlier by explicit request.

Projection cannot return a valid molecule without changing the rule. Dropping the surviving bond,
overlay, or constraint would impose cascading or SqPO-style closure rather than faithfully project
the stored side.

Consequently:

- `ReactionSpanAst::try_from_entries` checks union-frame references, not side-level DPO semantics;
- `DpoValidator` owns the dangling-freedom check;
- parsing a structurally intact span does not imply DPO validation;
- `ReactionAst::to_reaction_span` preserves a DPO-invalid but representable reaction;
- projecting either side is fallible;
- projection must not drop surviving incidences or constraints merely because their references are
  absent from that side.

### Projection errors

`ReactionSpanAst::lhs` and `ReactionSpanAst::rhs` should return
`Result<MoleculeAst, MoleculeEntriesError>`. A failed projection has assembled a molecule whose
entries contain an unavailable reference; this is exactly the target constructor's error, not a DPO
contradiction and not a new chemical error.

The shared implementation should retain every entity and constraint present on the selected side,
remap all references that survive compaction, and report `MoleculeEntriesError::InvalidReference`
when a selected entry refers to an entity absent from that side. It must not obtain success by
dropping the referring entry.

This choice keeps the public distinction precise:

- `DpoValidator` answers whether the span satisfies DPO dangling-freedom;
- `lhs` or `rhs` answers whether that side can be represented as a `MoleculeAst` and returns it;
- a future cascading or SqPO projection, if needed, must have an explicit transformation name.

### Downstream API consequences

Making side projection fallible changes the direct span conversion surface:

- `ReactionSpanAst::to_reaction` becomes fallible because it needs a valid lhs molecule. A
  DPO-invalid rhs does not by itself prevent this conversion; `ReactionAst` can represent the same
  dangling rule.
- `ReactionSpanAst::correspondence` remains infallible because it reads span states without
  materializing either molecule.
- `ReactionSpanAst::superimpose` remains an asserted producer: two referentially intact molecules and
  a valid correspondence produce projectable sides by construction.
- property tests over generated valid reactions unwrap the projections as established invariants;
  tests over arbitrary structurally intact spans must admit projection failure.
- Python `lhs`, `rhs`, and `to_reaction` raise `InvalidStructureError` when projection fails. This is
  a failed operation precondition on an already constructed object, rather than invalid Python
  constructor arguments.

This does not by itself change the public contracts of `ReactionAst::reverse`, reaction composition,
or reaction fingerprints. Those are operations on `ReactionAst`, and their current use of
`to_reaction_span().rhs()` is an implementation choice rather than part of their semantics. They
must handle the fallible internal projection within their existing operation-specific result
surface or materialize the product without routing through a public span projection.

`ReactionAst::to_reaction_span` retains its existing `Contradiction` surface in this work. A broader
review may distinguish contradictory deltas from invalid delta references, but the span-projection
change is not a reason to introduce a shared reaction-conversion error hierarchy.

## Review questions

When adding or revising an API, ask in order:

1. Can the type store this value without unresolved references or malformed representation state?
2. If yes, is the rejected property instead a semantic predicate that should have a named
   validator?
3. Does the conversion preserve all information the target can represent?
4. Is any dropping, normalization, repair, or resolution being hidden inside construction or
   conversion?
5. Does each failure report the boundary that actually failed rather than a broader semantic label?

If the answer to question 1 is yes, a semantic failure alone is not a reason to reject construction.
