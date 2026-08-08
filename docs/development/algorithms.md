# Algorithm execution contracts

This is a normative developer guide for the execution and result-delivery API
of graph algorithms in `umol-graph-core` and of the higher-layer operations
that wrap them.

## Selection and delivery are independent

An `*Algorithm` enum selects one implementation of one algorithmic problem.
How results reach the caller is a separate axis with three forms:

- `visit_*` — callback delivery; the visitor returns `ControlFlow`, and
  `Break` terminates the search;
- `enumerate_*` — eager collection of every result;
- `iter_*` — a resumable iterator with an explicit search cursor.

The naming for these forms is defined in the nomenclature guide under *Result
delivery*.

## Rules

- An operation returning a single value — one output struct, one set, one
  coloring, one count, an `Option` — has no visitor form.
- Where a visitor exists, the eager form collects through it; one operation
  never has two independent search implementations. Before a cursor exists,
  `visit_*` owns the search and `enumerate_*` collects it. Once an `iter_*`
  cursor exists, the cursor owns the search state and both other forms are
  derived from it.
- Every algorithm variant behind one public selector obeys the same delivery
  contract before a visitor is exposed. A visitor that streams for one variant
  but materializes for another is not a streaming API.
- A wrapper layer that mirrors an operation's `enumerate_*` form also mirrors
  its `visit_*` form. Current consumer counts do not gate the mirror; they are
  a poor indicator of usage, and the mirror is a wrapper, not algorithmic
  work.
- Add `iter_*` only when suspension across calls has a consumer. An operation
  named `iter_*` must not materialize its complete result set. The cursor
  borrows the graph it searches; an owning adapter for an FFI boundary is
  layered on top of the cursor, not built into it.

## Visitor payload

Low-level visitors emit borrowed data. The payload borrows the search's
working state, reusing a per-search scratch buffer where the leaf
representation is not already contiguous. Owned, typed result assembly
(`Correspondence`, `GraphCorrespondence`) belongs to the collector and to any
caller that retains results.

Visitors are themselves an efficiency measure, so efficiency takes precedence
for these low-level methods; signature asymmetry across visitor families is
acceptable and is not a reason to forfeit zero-copy delivery. A deviation to
an owned payload is justified only when the borrowed form has bad ergonomics.
Higher-layer operations (molecule-level wrappers) may emit owned typed
results.

## Emission order

Traversal is deterministic for a fixed graph representation, but emission
order is not a contract. Consumers must not depend on the order or on the
identity of an early-terminated prefix. Visitor and eager forms agree as
normalized result sets, and property tests compare accordingly across every
selector variant.

## Streamability

A `visit_*` form requires that result membership is decidable at emission:
every value passed to the visitor is a final member of the result set, never
revised or discarded by later search.

- A leaf-emitting backtracking search satisfies this: a completed leaf is a
  result.
- An exact bound computed before the search satisfies it: every emission
  meeting the bound is final. Matching enumeration streams this way, with the
  exact target cardinality established up front.
- Incumbent-based branch and bound does not satisfy it: the bound rises
  during search and earlier best results are superseded. Such an operation
  stays eager, or gains a separate fixed-target enumeration phase run after
  the optimum is known.

Early termination must stop the search itself, not truncate a materialized
result vector.
