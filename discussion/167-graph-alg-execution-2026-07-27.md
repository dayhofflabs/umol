# Graph algorithm execution APIs

Status: **Active**
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

Settled 2026-08-07:

- An operation returning a single value — one output struct, one set, one
  coloring, a count, an `Option` — takes no visitor form and keeps its
  descriptive name.
- Every eager operation returning a collection of results is named
  `enumerate_*`, whether or not a visitor form exists or is planned. The
  plural-noun eager forms are renamed; see the ledger below.

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

## Streamability requirement

Settled 2026-08-07. A `visit_*` form requires that result membership is
decidable at emission: every value passed to the visitor is a final member of
the result set and is never revised or discarded by later search.

- Leaf-emitting backtracking searches satisfy this: a completed leaf is a
  result. This covers subgraph isomorphism (all six selectors), the
  common-subgraph clique walks, paths, connected subgraphs, and the cycle
  enumerations.
- An exact bound computed before the search satisfies it. Matching enumeration
  streams because Edmonds supplies the exact target cardinality up front, so
  every emitted matching is final.
- Incumbent-based branch and bound does not satisfy it: the bound rises during
  search and a "best so far" may be superseded. `maximum_common_*` (McGregor)
  therefore has no visitor form.

### Streaming and maximum common subgraph

Reviewed 2026-08-07. None of the established maximum-common-subgraph
algorithms — McGregor (1982) backtracking with bounds, the McSplit family
(McCreesh, Prosser, Trimble 2017), clique reductions over the modular product
with similarity-threshold pruning (RASCAL; Raymond, Gardiner, Willett 2002) —
emits final results before the search completes. Anytime branch and bound
emits a stream of improving incumbents, but each emission supersedes the
previous one; that is a different contract from `visit_*`, where every
emission is final, and would need its own name if ever exposed.

The streamable construction is two-phase: run branch and bound to establish
the optimum size, then stream every correspondence of exactly that size with a
fixed-target search — the same shape `visit_maximum_matchings` uses with its
Edmonds bound, at the cost of roughly one additional search pass. Deferred: no
consumer requests it.

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

Settled 2026-08-07: cursors borrow the graph they search, like the existing
search-state types. The Python boundary cannot hold a Rust lifetime, so the
S5 design includes an owning adapter layered on top of the cursor; the cursor
itself does not own or share the graph to serve FFI.

## Input-domain dispatch

Simple-graph versus multigraph execution is independent of result delivery.
Cycle enumeration keeps one algorithm selector and dispatches between the
direct simple-graph implementation and its subdivision-based fallback.
`try_*`, explicit fallback, and combined operations retain the naming settled
in doc 158.

## Delivery decisions

Settled 2026-08-07, satisfying the completion criterion that every
multi-result algorithm has an explicit decision about visitor support.

| operation | decision | basis |
| --- | --- | --- |
| `subgraph_isomorphisms` (+ `_at`), all six selectors | visitor | leaf-push backtracking in every selector variant |
| `enumerate_common_subgraphs`, both selectors | visitor | clique/backtracking walks; leaves are final |
| `maximal_common_subgraphs` (Bron–Kerbosch) | visitor | maximal cliques are final at emission |
| `maximum_common_{induced,edge}_subgraphs` (+ `_seeded`) | no visitor | incumbent bound revises results mid-search |
| `enumerate_connected_subgraphs` (ESU) | visitor | leaf-push; structural fingerprint benefits from peak-memory reduction |
| `enumerate_paths` | visitor, deferrable | leaf-push; no consumer exists |
| `automorphisms`, `canonical_key` | no visitor | nauty owns the search; output is one compact struct |
| `connected_components`, `biconnected_components` | no visitor | one linear pass producing one partition |
| simple/relevant cycles, perfect/maximum matching enumeration | visitor, done | reference pattern |
| single-value operations (`minimum_cycle_basis`, `unique_ring_families`, `maximum_independent_set`, `bipartition`, `neighborhood`, `refine`, `circular_refine`, single matchings, counts, shortest cycles, `topological_order`) | no visitor | single-value rule |

### Renames

Eager collection-returning operations take the `enumerate_*` prefix.
Plural-sounding operations returning one value keep their names:
`automorphisms` (one `AutomorphismOutput`), `unique_ring_families` (one
`UniqueRingFamilies`), `maximum_independent_set` (one set).

| current | renamed |
| --- | --- |
| `Graph::subgraph_isomorphisms` / `_at` | `Graph::enumerate_subgraph_isomorphisms` / `_at` |
| `Graph::maximal_common_subgraphs` | `Graph::enumerate_maximal_common_subgraphs` |
| `Graph::maximum_common_induced_subgraphs` | `Graph::enumerate_maximum_common_induced_subgraphs` |
| `Graph::maximum_common_edge_subgraphs` / `_seeded` | `Graph::enumerate_maximum_common_edge_subgraphs` / `_seeded` |
| `Graph::connected_components` | `Graph::enumerate_connected_components` |
| `Graph::biconnected_components` | `Graph::enumerate_biconnected_components` |
| `GraphView` mirrors (subgraph isomorphism ×2, components ×2) | the same `enumerate_*` names |

### Payload, ordering, wrappers, cursors

Settled 2026-08-07, recorded normatively in the algorithm execution guide:

- Low-level visitors emit borrowed payloads — zero-copy, with a per-search
  scratch buffer where the leaf representation is not already contiguous.
  Owned typed assembly (`Correspondence`, `GraphCorrespondence`) belongs to
  the collectors and to higher layers. Asymmetry across visitor families is
  acceptable; a deviation to an owned payload needs an ergonomic reason,
  efficiency otherwise takes precedence for low-level methods.
- Emission order is deterministic for a fixed representation but is not a
  contract; visitor/eager agreement is normalized-set equality.
- A wrapper layer that mirrors an eager form mirrors the visitor form;
  current consumer counts do not gate parity. `GraphView` gains visitor
  mirrors for every eager form it wraps.
- `iter_*` cursors borrow the graph. The owning adapter needed at the Python
  boundary is part of the S5 design, layered on top of the cursor.
- The molecule-level visitor name is `visit_substructure_matches`.

## Staged migration plan

Recorded 2026-08-07. Stages end green; a breaking subitem and the caller
migration restoring green share a stage. Cross-crate stages verify with
`--all-features --tests` plus clippy.

Uptake motivating the order: no crate outside `umol-graph-core` calls any
`visit_*` form. The full substructure match set is materialized for an
existence check in constraint validation, by `ReactionAst::apply` before its
lazy derivation loop, by every Python `substructure_matches` call, and per
template in the pattern fingerprint.

### S0 — renames (breaking) **Done**

- S0a: subgraph isomorphism ×2 plus `GraphView` mirrors; callers:
  substructure search (2 sites), `GraphView` (2), tests, benches.
- S0b: common-subgraph maximal/maximum ×4; callers are graph-core tests and
  benches only.
- S0c: components ×2 plus `GraphView` mirrors; callers: constraint validation
  (2 sites), HMO aromaticity, tests, benches.

### S1 — subgraph-isomorphism streaming core (green) **Done**

- S1a: convert the six per-selector searches to `ControlFlow` emitters over
  `&[usize]` embeddings, replacing the leaf pushes; internal only. Ullmann
  and RI stop cloning their mapping buffer per result.
- S1b: public `visit_subgraph_isomorphisms` (+ `_at`) emitting the borrowed
  embedding `&[NodeId]` (query position → host node) from a per-search
  scratch buffer; the eager forms collect the visitor and lift each embedding
  to `Correspondence<NodeId>`. Property tests: visitor and eager agree for
  all six selectors, ordinary and anchored; early termination makes strictly
  fewer `node_match` calls than a full search on a multi-match fixture.
  Retires the first-match TODO in the source. [dep: S1a]
- S1c: benches gain time-to-first-result (`Break` after one) and count-only
  peak-storage rows beside the existing totals. [dep: S1b]

### S2 — substructure streaming (green) **Done**

- S2a: both substructure strategies drive the visitor internally with overlay
  verification inside the callback; public behavior unchanged. [dep: S1b]
- S2b: `MoleculeAst::visit_substructure_matches` emitting
  `MoleculeCorrespondence`; `substructure_matches` collects it. [dep: S2a]
- S2c: the constraint-validation existence check breaks on the first
  admissible match. [dep: S2b]
- S2d: the pattern fingerprint consumes the visitor per template. [dep: S2b]
- S2e: `GraphView` visitor mirrors for every eager form it wraps — matchings
  and cycles against the existing visitors, subgraph isomorphism against
  S1b. Wrappers only, no algorithmic work. [dep: S1b]

### S3 — common-subgraph visitors (green)

- S3a: emitters for the modular-product clique walk, Bron–Kerbosch, and
  direct backtracking, emitting the borrowed pair slice `&[(NodeId, NodeId)]`
  from a per-search scratch buffer; the batch clique conversion is removed.
- S3b: public `visit_common_subgraphs` and `visit_maximal_common_subgraphs`
  emitting the borrowed pair slice; eager forms collect and lift each result
  to `GraphCorrespondence`; tests and benches as in S1. [dep: S3a]
- S3c: `ReactionAst::compose` composes per overlap inside the visitor; public
  signature unchanged. [dep: S3b]

### S4 — connected subgraphs and paths (green)

- S4a: ESU emitter and `visit_connected_subgraphs` emitting the borrowed
  working buffer `&[EdgeId]`; the eager form collects by copying.
- S4b: the structural fingerprint computes canonical keys per emitted
  subgraph without materializing the family. [dep: S4a]
- S4c (deferrable): paths visitor; no consumer exists.

### S5 — deferred: cursors and Python laziness

A separate design task (see the cursor section above): the owned per-selector
cursor state and the enum-of-cursors against shared-state decision precede any
code. The concrete consumer is a Python generator for substructure matches and
reaction application; every current Python iterator wraps an already-collected
vector, and visitor delivery cannot suspend across `__next__`. Gate on
committing to that Python surface.

Critical path: S0 → S1 → S2. S3 and S4 depend on S0 plus their own subitems
only. Every stage ends green and delivers value on its own, so the migration
can pause at any stage boundary.

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
