# Graph algorithm execution APIs

Status: **Proposed**
Date: 2026-07-27
Relates: [105](105-dsl-fixes-2026-06-06.md),
[158](158-ring-model-and-enumeration-2026-07-22.md),
[162](162-common-subgraph-algs-2026-07-25.md)

## Scope

This document tracks the progressive conversion of graph algorithms from eager
result collection to explicit result-delivery APIs. Algorithm selection and
result delivery are independent:

- the algorithm enum selects the algorithm;
- `visit_*`, `enumerate_*`, and `iter_*` describe how results are delivered.

The current graph-core API already has visitor and eager forms for matching and
cycle enumeration. Subgraph isomorphism, common-subgraph enumeration, paths,
and connected-subgraph enumeration remain eager.

## Naming contract

- `visit_*`: callback delivery with `ControlFlow` for early termination;
- `enumerate_*`: eager collection of every result;
- `iter_*`: a resumable Rust iterator with an explicit search cursor.

Do not use “iterative” to describe visitor delivery. An implementation may
remain recursive while visiting results.

Not every operation needs all three forms. `visit_*` plus `enumerate_*` is the
normal intermediate state. Add `iter_*` only when suspension across calls has a
consumer.

## Implementation layering

Before a resumable cursor exists:

```text
visit_* owns the search
enumerate_* collects visit_*
```

After a cursor exists:

```text
iter_* owns the search state
visit_* consumes iter_*
enumerate_* collects iter_*
```

All algorithm variants behind one public selector must obey the same delivery
contract before the visitor is exposed. A visitor that streams for one variant
but materializes for another is not a streaming API.

## Priority

1. Subgraph isomorphisms: highest-value visitor conversion because it enables
   existence checks, limits, reaction application during matching, and reduced
   Python-side materialization.
2. Common and maximal common subgraphs.
3. Paths and connected subgraphs.
4. Resumable iterators where a Python generator or long-lived Rust cursor
   justifies the explicit backtracking state.

Matching and cycle enumeration already establish the visitor/collector
pattern. Their APIs should be used as the convention, not reimplemented through
another abstraction.

## Subgraph-isomorphism cursor

A visitor conversion is local: replace solution-vector pushes at search leaves
with an emitter and thread `ControlFlow` through recursion.

A true iterator is a separate design task. It requires an owned cursor for each
algorithm containing:

- current mapping and inverse mapping;
- candidate position at every depth;
- restoration state for domains and terminal sets;
- static preprocessing results;
- algorithm-specific ordering or reduction data.

Before implementing cursors, decide whether all selector variants expose
separate cursor types behind an enum or share a sufficiently concrete common
state. Do not erase algorithm-specific state behind boxed callbacks merely to
force one type.

## Input-domain dispatch

Simple-graph versus multigraph execution is independent of result delivery.
Cycle enumeration keeps one algorithm selector and dispatches between the
direct simple-graph implementation and its subdivision-based fallback.
`try_*`, explicit fallback, and combined operations retain the naming settled
in doc 158.

## Verification

- Visitor and eager forms return the same normalized result set.
- Early termination stops the search rather than only truncating a materialized
  vector.
- Property tests cover all selector variants.
- Benchmarks distinguish preprocessing cost, time to first result, total
  enumeration time, and peak result-storage cost.
- External implementations may generate one-time fixtures but are never test
  dependencies.

## Completion criteria

- Every multi-result algorithm has an explicit decision about visitor support.
- Eager APIs collect the visitor where a visitor exists.
- No operation named `iter_*` materializes its complete result set.
- Python laziness is backed by resumable state rather than an eagerly filled
  Rust vector.
