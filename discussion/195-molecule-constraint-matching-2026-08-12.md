# 195 — Matching molecule-level constraints

Status: Proposed
Date: 2026-08-12
Relates: [128](128-substructure-derived-predicates-2026-06-23.md),
[165](165-ast-api-worklist-2026-07-27.md),
[193](193-subpattern-constraints-2026-08-09.md),
[194](194-constraint-assertion-semantics-2026-08-10.md)

## Purpose

Design and implement the evaluation of molecule-level constraints during substructure matching.
The matcher today evaluates only per-entity inline constraints; a pattern's `:constraints` list —
`Relational` and `Molecule` leaves and the `And`/`Or`/`Not` combinators — is silently ignored
(`substructure.rs` reads only the atom and bond containers). A pattern carrying
`(charge-sum 0)` matches as if unconstrained.

## Interim policy

Until this design lands, matching rejects a pattern whose molecule-scope `Constraints` list is
non-empty (the doc 194 S1a gate), with an error naming the construct. A query using an
unevaluated construct fails loudly rather than returning silently weakened results; fail →
evaluate is the compatible evolution, while the previous silent ignore admitted false positives.

The gate also closes the matcher path for reactions whose LHS carries a molecule-scope
constraint (every `ConstraintDelta::Remove`-style rule that asserts the removed constraint on
its LHS); `apply_at` with a supplied correspondence is ungated. Six tests exercising
`TransactionError::MissingEntry` through that path are `#[ignore]`d until evaluation lands
(2026-08-13); dropping the gate must un-ignore them:

- `umol-graph-ir/src/ir/reaction.rs`: `test_reaction_application_iter_error`,
  `test_reaction_products_iter_error`, `test_reaction_apply_error` (`case_1_transaction`);
- `umol-graph-ir/tests/property/reaction/malformed.rs`: `test_reaction_apply_error`;
- `umol-py/src/reaction.rs`: `test_reaction_application_iter_error`,
  `test_reaction_products_iter_error`.

## Scope

Doc 194 normalizes constraint placement only through resolution, and patterns are never resolved:
a pattern's list may carry bare entity leaves alongside `Relational` and `Molecule` leaves and
combinators. Matching must therefore evaluate:

- bare entity leaves, bound to host entities through the match correspondence and evaluated like
  inline predicates through the doc 194 constraints views;
- entity leaves carrying ring keys (`#R`, `#x`, `#y`), bare or inside combinators, which need the
  ring context: the `RingSet` built once per match run (fixed Relevant projection, size 22)
  serves the inline predicate path and this evaluation alike;
- `Molecule` leaves (`ChargeSum`, `BondOrderSum`, `UnpairedElectronCoupling`, `Connected`);
- `Relational` leaves, whose entity references are pattern ids and must be bound to host entities
  through the match correspondence;
- combinator trees whose entity leaves are evaluated the same way.

## Required design

- Scoping per predicate kind: whether a `Molecule` predicate constrains the matched image or the
  whole host — `Connected` over the image is a different assertion than over the host, and the
  same question applies to every sum predicate.
- Evaluation point: combinators over entity leaves break per-entity independence, so evaluation
  completes at embedding completion; decide what is soundly prunable earlier and what remains a
  post-filter.
- Reading: host evaluation uses the closure reading (`derived_complete`), mirroring `satisfies`.
- Cost model at reaction-network scale: per-embedding evaluation of relational predicates, gated
  on presence as the ring machinery gates today; one `RingSet` per match run shared between the
  inline predicate path and molecule-scope evaluation, built only when ring keys occur anywhere
  in the pattern.
- Interaction with doc 193: a recursive subpattern constraint is a molecule-scope constraint of
  the same family; the evaluation surface designed here is the one doc 193 plugs into.
- Conformance and property coverage for the complete contract.

An implementation plan belongs here after these questions are settled.
