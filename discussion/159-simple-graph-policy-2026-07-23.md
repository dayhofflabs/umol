# Simple-graph policy

Status: **Active**
Date: 2026-07-23
Relates: [118](118-validator-architecture-2026-06-20.md),
[124](124-tier1-structural-wellformedness-2026-06-21.md),
[131](131-reaction-application-design-2026-06-24.md),
[136](136-dpo-primitives-2026-07-04.md),
[148](148-validated-transactions-operations-2026-07-13.md),
[158](158-ring-model-and-enumeration-2026-07-22.md)

## Scope

This document records the treatment of self-loops and parallel edges across
graph storage, molecular input, validation, generated structures, and graph
algorithms. It also records the narrow decision needed by doc 158.

The general policy is tentative because its exact reaction semantics remain
open. That does not block ring work: graph-core cycle semantics admit
non-simple graphs, while the molecule validator separately reports self-loops
and parallel localized bonds as structural contradictions.

## Current behavior

### Graph storage is permissive

`umol-graph-core::Graph` can represent self-loops and parallel edges.
`Graph::new` constructs its adjacency representation directly from an edge
list, and `Graph::add_edge` appends an edge without imposing simple-graph
semantics.

This is a property of the storage representation, not evidence that every
graph algorithm has meaningful multigraph semantics. The molecular domain
currently has no identified need for self-loops or parallel localized bonds.
The graph-core consumers that deliberately exercise those cases are tests of
storage behavior or individual algorithms.

Changing `Graph::new` into a fallible simple-graph constructor would be a broad
API migration. It would also prevent graph storage from representing invalid
input long enough to diagnose it. Neither consequence is required for the
ring work.

### Molecular construction preserves invalid structures

The low-level molecule construction paths inherit the storage behavior:

- `MoleculeAst::from_parts` passes localized bond endpoints to `Graph::new`;
- `MoleculeEditor::add_bond` delegates to `Graph::add_edge`;
- `MoleculeBuilder` bond construction delegates to `MoleculeEditor`;
- checked edit transactions validate references and ranges but do not reject a
  topologically conflicting added bond;
- molecule DSL and TableIR lowering ultimately construct a `MoleculeAst` from
  parts.

These paths can therefore represent a localized bond self-loop or parallel
localized bonds. That is useful for input diagnostics, but construction alone
does not certify structural validity.

SMILES parsing already follows this separation. The syntax parser can preserve
a syntactically valid ring closure that produces a self-loop in TableIR.
Topology validation is a later concern. Parsing an external representation and
accepting it as a valid molecule are distinct operations.

### Validation already defines the domain restriction

`EntityStructureValidator` reports localized bond self-loops and parallel
localized bonds as `BondSelfLoop` and `BondsParallel`. `BondViews::has_conflict`
provides the corresponding structural-conflict predicate used by generated
operations.

Reaction application validates its LHS and host, applies the DPO conditions,
and rejects a generated product with structural conflicts. Pushout operations
likewise use the conflict predicate to prevent publication of a structurally
conflicting result. These are existing examples of permissive representation
combined with stricter operational output guarantees.

Not every input path currently performs entity-structure validation before
returning its higher-level result. In particular, SMILES interpretation raises
TableIR and runs resolution without an explicit entity-structure validation
stage. The policy audit below must classify such boundaries before deciding
whether this is a defect in the operation or an intentionally unvalidated
result.

## Tentative policy

Representability, validity, and operational guarantees are separate:

| Boundary | Non-simple structure | Responsibility |
| --- | --- | --- |
| Parsed or otherwise user-provided input | May be represented | A validator reports structural contradictions |
| Low-level AST construction and mutation | May be represented | The caller has no validity guarantee until validation |
| An operation documented to return a valid domain object | Must not publish a non-simple result | The operation validates before publication |
| Algorithmically generated structure | Must not introduce a self-loop or parallel edge | The generating operation rejects the candidate or returns an error |
| A graph algorithm defined only for simple graphs | May receive permissive graph storage | The algorithm checks its precondition and returns an error |

In particular:

- permissive graph storage does not impose multigraph semantics on algorithms;
- user input remains inspectable and diagnosable instead of failing in a
  low-level constructor;
- algorithms do not silently reinterpret or simplify non-simple input;
- generated structures must not introduce a structural contradiction, even if
  their input representation was permitted to contain one;
- direct construction or editing remains lower-level than validation and does
  not acquire an implicit valid-by-construction promise.

“Generated” includes transformation, resolution, graph rewriting, reaction
application, and construction from an algorithm's output. Mere parsing or
faithful lowering of user input is not generation for this purpose.

## Rejected broad migration

Simple-graph enforcement should not be propagated through every constructor.
In particular, the present design does not:

- introduce a special bundle type solely to pass checked graph parts;
- make `Graph::new`, `MoleculeAst::from_parts`, or every editor and builder
  operation return `Result`;
- require every caller to prove simplicity before it can construct an
  inspectable value;
- use graph storage permissiveness as a reason to support multigraph semantics
  in every algorithm.

Such a migration would substantially complicate the ordinary construction API
without solving a current domain requirement. Validation and operation-specific
preconditions provide the narrower separation needed here.

## Ring-algorithm decision

The ring work in doc 158 does not add a public simple-graph precondition:

1. Graph-core cycles preserve edge identity and include self-loops and
   parallel-edge cycles.
2. The edge-aware Read--Tarjan implementation handles non-simple graphs
   directly.
3. Vismara, minimum-cycle-basis, and Unique Ring Family analysis use their
   direct simple-graph path when possible and an internal subdivision fallback
   otherwise.
4. Detection is fused into the mandatory initial graph analysis. An
   unmodified simple-graph algorithm has no reliable failure signal on which to
   base fallback.
5. Cycle visitation, cycle collection, basis computation, URF decomposition,
   `MoleculeAst::rings`, and `RingViews` remain infallible.
6. The molecule layer maps graph cycles to chemical rings and excludes cycles
   with fewer than three distinct atoms. Structural validation remains a
   separate operation.

The subdivision construction belongs in graph-core. The existing
`MoleculeAst::incidence_graph` is the overlay-capable molecule-level Levi graph;
no additional umol-ast subdivision type is required for ring perception.

No graph, molecule, editor, builder, parser, or transaction constructor changes
are prerequisites for doc 158.

## Reaction semantics remain open

The general generated-output rule is clear, but its application to reactions
needs a dedicated decision:

- whether a non-simple reaction LHS is an invalid rule, an input that remains
  representable until validation, or both at different API boundaries;
- whether a non-simple host is rejected before matching or produces an
  operation-level validation error;
- whether a rule capable of introducing a self-loop or parallel bond is invalid
  independently of a match;
- whether a conflict that arises only for a particular correspondence rejects
  that correspondence or fails the entire application;
- whether these failures belong to rule planning, application preconditions,
  or the existing per-application structural-conflict error;
- how the answers differ, if at all, between DPO and SqPO rewriting.

The existing reaction and pushout conflict checks are a sound baseline: they do
not publish generated non-simple products. This document does not strengthen
that baseline into a complete reaction policy.

## Follow-up audit

A later pass should classify every relevant API as input-preserving,
validating, or result-generating. The audit should cover:

- molecule and reaction DSL parsing;
- SMILES, MOL, and SDF parsing, raising, and ingestion;
- Rust and Python constructors;
- `MoleculeAst::from_parts`;
- `MoleculeEditor`, `MoleculeBuilder`, and transactions;
- reaction application and pushouts;
- resolvers and other transformations;
- induced-subgraph, extraction, compaction, and related graph operations.

For each result-generating operation, the audit should verify that structural
conflicts are rejected before the result is published. For input-preserving
operations, it should verify that a documented validation route exists. Test
coverage should be checked explicitly for every input avenue, especially
`MoleculeEditor` and `MoleculeBuilder`, rather than inferred from parser tests.

This audit is broader follow-up work. It does not block the total edge-aware
cycle semantics selected for ring algorithms.
