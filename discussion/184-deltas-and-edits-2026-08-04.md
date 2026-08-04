# 184 — Deltas and edits: why there are two

Status: Informational
Date: 2026-08-04
Relates: [043](043-mutative-undoable-mutation-2025-12-23.md),
[086](086-molecule-ast-api-2026-04-16.md),
[131](131-reaction-application-design-2026-06-24.md),
[148](148-validated-transactions-operations-2026-07-13.md),
[179](179-python-editing-and-transactions-2026-08-02.md)

`Delta` and `Edit` describe nearly the same kinds of molecular change. This makes the two
hierarchies look redundant, especially because reaction application eventually lowers deltas to
edits. They should nevertheless remain distinct: replacing either one with the other would obscure
which structure owns the operation and would weaken the useful guarantees of the replaced type.

## The answer in one line

**A delta is a rule-relative algebraic change; an edit is an ordered host-relative mutation plan.**

"Declarative versus imperative" is a useful shorthand, but the relative coordinate system is the
more exact distinction. A delta is interpreted in the id space owned by a reaction and can be
reasoned about before any host is selected. An edit names entities in one initial host, plus
same-transaction `New(n)` handles, and acquires its full effect only when applied to that host.

`Undo` is a third object with different semantics again: it records the concrete effect of one
successful edit application, including host-dependent cascades and id compaction.

| | `Delta` / `Deltas` | `Edit` / `Edits` | `Undo` / `Transaction` |
| --- | --- | --- | --- |
| relative to | reaction LHS and rule-owned id space | initial host and earlier operations | one realized application |
| describes | complete rule-side change | requested mutation steps | actual inverse operations |
| order | normalized by `canonicalize` | semantically significant | rollback order |
| inversion | `Delta::inverse` is total | not available from the plan alone | applied by rollback |
| structural removal | every deletion explicit | dependent structure may cascade | captures everything actually removed |
| primary use | reversal and rule composition | molecule construction and mutation | transactional recovery |

## Inversion terminology

A delta is not generally an involution and should not be called "self-inverting." Addition and
removal, for example, are different operations.

The precise statements are:

- the delta vocabulary is **closed under inversion**;
- `Delta::inverse` is total; and
- inversion is an **involution as a function**:
  `delta.clone().inverse().inverse() == delta`.

Every removal carries the removed entity, and every field modification carries both `old` and
`new`, so the inverse delta can be constructed without applying the rule. Reversing an entire
reaction still requires more than mapping `Delta::inverse`: `ReactionAst::reverse` also constructs
the product-side LHS and re-anchors ids into its compacted id space.

An edit often carries similar information for checked application, but that does not make the edit
plan statically invertible. Additions receive concrete host ids only during application; removals
compact host ids; and a topology removal may discard incident bonds, overlays, and constraints that
the edit did not name. Checked application observes these effects and records the corresponding
`Undo` values. That realized journal, rather than the original `Edit`, is what can restore the host.

## The two independent axes

Two distinctions are easy to conflate here.

### Rule-relative versus host-relative

A delta refers to entities in a reaction-owned frame. That frame survives matching, remapping,
canonicalization, reversal, and composition. The same rule can therefore be applied to different
hosts through different matches.

An edit refers to `Id(n)` in the transaction's initial host or to `New(n)`, the nth same-kind entity
created earlier in the same edit sequence. Its order is part of its meaning. Sorting, folding, or
deduplicating arbitrary edits could change which creations later handles name or move an operation
before its prerequisite.

The two addressing schemes are similarly expressive for construction: a reaction over an empty
LHS and an edit sequence can both build a molecule. That does not make them interchangeable. The
important difference is who owns the frame and which operations must remain meaningful before a
host has been chosen.

### DPO versus cascade deletion

DPO versus SqPO is an application policy, not the fundamental difference between `Delta` and
`Edit`.

The current reaction path is strict DPO. `DpoValidator` checks that deleting an LHS atom explicitly
deletes its incident LHS structure, and reaction application performs match-local dangling checks
against the host before lowering the deltas to edits. The resulting plan is then executed by the
same edit engine whose `RemoveTopology` operation is capable of cascading dependent removals. The
DPO preconditions make that cascade unnecessary for a valid application. This combination already
demonstrates that a DPO rule representation and a cascade-capable mutation primitive are compatible.

Conversely, deltas could be interpreted with SqPO semantics without replacing their vocabulary.
Application could admit a match with external incidence and remove the host-only dependent
structure automatically. The cost is semantic: the original deltas would no longer describe the
complete realized transformation, and applying the syntactic inverse rule could not recreate the
host-only structure that was discarded. A realized derivation or undo journal would have to retain
that additional information.

`RemoveTopology` is therefore best described as **cascade deletion** or **SqPO-style dangling
cleanup**, not as a complete implementation of SqPO rewriting. An edit plan has no rule match on
which to define non-injective matching, cloning, or a final pullback complement. Those concerns
belong to the rewrite application layer, not to an imperative mutation vocabulary.

## Why neither hierarchy should replace the other

### Replacing deltas with edits

A host-independent "edit" would need stable rule-owned ids, complete removed values, canonical
folding, inversion, remapping, and composition. It would also have to exclude or separately record
host-dependent cascade effects. Adding those properties recreates `Delta` under a less precise name.

Using current `Edits` directly would instead lose the algebra that reactions need. `New(n)` and
sequence order are appropriate for construction, but they are the wrong basis for canonical rule
identity and composition.

### Replacing edits with deltas

Using deltas for ordinary host mutation would require the caller or a planning pass to materialize
every old value and every entity that a removal will affect, assign a stable frame, and construct a
complete before/after difference before mutation begins. Construction through newly issued handles
would become less direct, and transactional failure would still require a realized undo journal.

That planning layer would recreate `Edits` and `Undo`, while making the ordinary mutation API carry
information that is available more reliably from the host during application.

## Shared vocabulary without a shared hierarchy

The overlap between the two representations is real. Field-change and constraint-change values
with `old` and `new` sides should be shared wherever their semantics are identical, as they already
are in several cases. Common projection or validation code may also be factored when the abstraction
is complete.

This does not justify parameterizing one public top-level mutation hierarchy by its index space.
The important differences are not only the index type: canonicalization, ordering, completeness,
cascade discovery, and inversion all differ. A generic hierarchy would either expose operations
that are invalid for some instantiations or encode the distinction through a second layer of traits
and policies.

## Possible extensions

### Deltas under SqPO application

This is possible, but no new delta variants are required. It would be an alternative reaction
application policy with a richer realized derivation. It would weaken the present and useful
property that the rule explicitly accounts for every structural deletion, while adding no required
functionality at present. Strict DPO should remain the only delta application semantics until an
excision-style consumer justifies the tradeoff.

### Exact removal for edits

Edits already expose the useful cascade-removal operations: `remove_atom`, `remove_bond`, and
`remove_topology` all produce `RemoveTopology`, and removal of affected structure and constraints is
captured during the transaction. There is no missing "full set of SqPO verbs" to add. SqPO is not a
verb family at this layer.

The meaningful complementary extension would be an **exact** or **strict** removal operation that
rejects the edit when unnamed dependent structure remains. It could serve two possible consumers:

- callers that want host mutation to fail rather than cascade; and
- reaction lowering, where it could check at execution time that prior DPO validation and lowering
  really produced a complete removal plan.

Neither consumer currently requires a new public operation. The existing reaction path establishes
DPO preconditions before mutation, while ordinary editing deliberately uses cascade deletion.
Adding a removal-policy enum or parallel variants now would encode an unused distinction and expand
an already large edit vocabulary. The strict operation should be designed when a concrete consumer
can determine whether the distinction belongs in an edit variant, an `Edits` method, or transaction
application configuration.

## Conclusion

The durable model has three parts:

1. **Deltas** are complete, rule-relative changes suitable for algebraic operations.
2. **Edits** are ordered, host-relative plans suitable for construction and mutation.
3. **Undo journals** are realized, host-relative inverses suitable for rollback.

DPO versus SqPO remains an orthogonal choice at application. The current asymmetry—strict DPO for
reaction rules and cascade deletion for direct edits—is deliberate and useful, but it is not the
reason the two mutation hierarchies exist.

## Lineage

- [043](043-mutative-undoable-mutation-2025-12-23.md) — the origin of a modification as a
  serializable value with an inverse, from the JSON-patch model.
- [086](086-molecule-ast-api-2026-04-16.md) — the `MoleculeAst` API within which editing sits.
- [131](131-reaction-application-design-2026-06-24.md) — reaction application, where deltas meet a
  host through a match, and the DPO/SqPO policy distinction.
- [148](148-validated-transactions-operations-2026-07-13.md) — validated transactions and the
  operation lifecycle; why a batch is the unit.
- [179](179-python-editing-and-transactions-2026-08-02.md) — the edit vocabulary and transactions on
  the Python surface, and the standalone edit document.
