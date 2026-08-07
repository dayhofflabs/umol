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

## Reaction-span construction

`ReactionSpanAst` stores the union-frame encoding of an actual span of two molecules. The union
namespace is necessary to share entity identity across the sides, but union-frame integrity alone is
not sufficient: the entries selected on the left and right must each form a referentially intact
`MoleculeAst`.

For example, an `Unchanged` bond incident to a `Removed` atom is representable as an arbitrary
annotated union graph, but its right projection contains an edge without both endpoints. It is
therefore not a reaction span. Accepting it and dropping the bond during projection would impose a
cascading or SqPO-style transformation rather than faithfully represent the supplied entries.

This is a representation invariant of `ReactionSpanAst`, not chemical validation. A type named as a
span must contain both objects of the span. Consequently:

- `ReactionSpanAst::try_from_entries` checks the union namespace and the referential integrity of
  each projected side;
- `ReactionSpanAst::from_entries` asserts the same invariant for trusted producers;
- parsing accepts only entries that construct an actual two-sided span;
- `ReactionSpanAst::lhs`, `rhs`, `to_reaction`, and `correspondence` are infallible;
- projection retains every entity and constraint present on the selected side and remaps its
  references without repair or silent loss.

The two reaction representations deliberately have different construction boundaries.
`ReactionAst` is a permissive lhs-plus-deltas carrier and may contain a DPO-invalid rule.
`ReactionAst::to_reaction_span` is therefore fallible: it rejects deltas whose projected product
cannot form the second molecule of an actual span. The conversion retains its existing
`Contradiction` surface pending the repository-wide error review. Conversely, every constructed
`ReactionSpanAst` has a valid lhs, so `ReactionSpanAst::to_reaction` is infallible.

### DPO validation

The current `DpoValidator::validate_reaction_span` check becomes redundant. A removed atom with a
surviving bond or overlay would already violate the right-side construction invariant. The
constructor's symmetric check is stronger: it also rejects an added atom required by a surviving
left-side entity and covers stereo and constraint references.

Remove the reaction-span validator entry point rather than retaining a validator that can only
confirm a type invariant. `DpoValidator::validate_reaction` remains useful for permissive
`ReactionAst` values, and match-dependent dangling checks remain part of reaction application: they
concern the supplied host and match, not the reaction span alone.

This keeps fallibility at the representation boundary. `ReactionAst::reverse`, reaction
composition, and reaction fingerprints retain their operation-specific result surfaces; an
internal route through `ReactionAst::to_reaction_span` does not require side projection itself to
become fallible. Python checked construction reports invalid span entries as `ValueError`, while
`lhs`, `rhs`, and `to_reaction` return their values directly.

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
