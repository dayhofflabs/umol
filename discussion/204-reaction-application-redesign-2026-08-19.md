# 204 — Reaction application result redesign

Status: Superseded
Date: 2026-08-19
Relates: [131](131-reaction-application-design-2026-06-24.md),
[187](187-assembly-disassembly-2026-08-05.md),
[201](201-molecular-data-first-steps-2026-08-19.md),
[203](203-atom-mapping-2026-08-19.md),
[218](218-mutation-witness-2026-08-31.md),
[data-type contracts](../docs/development/data-types.md),
[nomenclature](../docs/development/nomenclature.md)

## Superseded

Doc [218](218-mutation-witness-2026-08-31.md) subsumes this work. It retains
the diagnosis below, removes `ReactionDerivation`, and defines reaction
application in the same covariant-witness vocabulary as other molecule
operations. The historical text remains below as the input that led to that
decision.

## Scope

Reconsider the result contract of `Reaction::apply` and decide whether
`ReactionDerivation` should remain an operation-issued result, be renamed or
reshaped around genuinely application-specific information, or be removed.
The design must establish which information a successful application needs to
retain, including the concrete reaction, host-to-product correspondence,
rule-to-host match, and any useful operation provenance. It must cover the Rust
and Python application iterators, the shape of one item in a multi-result
application, and the product-oriented `React` adapter.

`Reaction` remains the central lhs-plus-deltas representation. `ReactionSpan`
remains its materialized superimposed form, and atom mapping continues to
return `Correspondence<AtomId>` values that may be lifted through
`Reaction::from_sides`. This document does not yet select an application result
type, public API, migration, or staged implementation plan.

## Justification

For every materializable reaction, `Reaction` and `ReactionSpan` encode the
same reaction morphism in operational and superimposed forms. A span is
equivalently the two molecule sides plus their correspondence, modulo
lhs-anchored reindexing and reaction normal form. The current
`ReactionDerivation` stores that same side pair and correspondence, so its only
additional representational distinction is the exact rhs frame.

That distinction does not by itself justify a separate semantic reaction type.
`ReactionDerivation` also does not retain the originating reaction or the
rule-to-host match, so it is not a complete application record. Its public
`chain` operation accepts independently supplied derivations without checking
agreement of the intermediate molecule sides. Before the type is expanded or
used as a persistence boundary, reaction application should be redesigned from
the information its consumers actually need and the redundant type should be
removed if no independent role remains.

The original motivation nevertheless remains substantial. Applying a reaction
or rule to a host is a plural operation: each successful match emits one result
whose essential payload is the product side plus the host-to-product
correspondence. `ReactionDerivation` can be understood as the named composite
return type for that item rather than as another authoritative reaction model.
The redesign must evaluate that narrow role directly. In particular, it should
separate the usefulness of an owned application item from the current choice to
also store the known host and expose general `reverse`, `chain`, and
`to_reaction` operations. Keeping a focused composite result remains a live
option even if the broader derivation abstraction does not.
