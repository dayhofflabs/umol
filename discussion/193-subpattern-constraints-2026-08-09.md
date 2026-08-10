# 193 — Recursive subpattern constraints

Status: Proposed
Date: 2026-08-09
Relates: [164](164-dsl-edn-worklist-2026-07-27.md),
[165](165-ast-api-worklist-2026-07-27.md), [176](176-ast-naming-2026-07-31.md),
[192](192-python-api-type-roles-2026-08-09.md)

## Purpose

Design recursive subpattern constraints as one coherent feature before reintroducing them to the
graph IR, DSL, and Python API. Their temporary removal in doc 176 is a postponement, not a rejection
of recursive pattern matching as a useful operation.

The removed implementation allowed a molecule constraint to contain another complete molecule and
an anchor relating outer and inner entities. That recursive carrier reached across the IR, metadata,
parsing, matching, validation, remapping, mutation machinery, and Python ownership model before the
end-to-end contract was settled. Retaining isolated pieces would make those pieces accidental
constraints on the eventual design.

## Required design

The next design pass must settle the feature as a whole:

- define how recursive predicates relate to ordinary substructure search and to the non-recursive
  entity and relational constraint model;
- define the anchor or correspondence representation across all eight entity kinds, including
  structural integrity, cardinality, and the semantics of an unanchored pattern;
- decide whether nested patterns have namespaces and metadata, and how outer and inner references
  are resolved without conflating their id spaces;
- specify the DSL and EDN representation, including named versus anonymous inner patterns;
- separate construction integrity from matching-time semantic validation and state the required
  matching algorithms explicitly;
- define compaction, remapping, ordering, normalization, and eventual canonicalization behavior;
- define how recursive constraints participate in deltas, constraint edits, transactions, and
  removal of referenced entities;
- define Python ownership and mutability without returning disconnected mutable copies from an
  immutable constraint value; and
- provide unit, property, conformance, benchmark, and fuzz coverage for the complete contract.

## Constraints on reintroduction

The feature should not return as a dormant enum variant or as a Rust-only carrier. Reintroduction
must provide a usable path from syntax through validation and matching, with Rust and Python
construction semantics that agree. The public happy path should remain a recursive constraint, not
a collection of exposed lifecycle adapters.

An implementation plan belongs here after these questions are settled.
