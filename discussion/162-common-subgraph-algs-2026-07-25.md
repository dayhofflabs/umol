# Common-subgraph enumeration alternatives

Status: **Completed**
Date: 2026-07-25
Relates: [108](108-mcs-algorithms-2026-06-09.md),
[135](135-reaction-composition-completeness-2026-07-01.md),
[151](151-python-molecule-workflows-2026-07-13.md),
[161](161-property-tests-as-specs-2026-07-25.md)

## Scope

`Graph::enumerate_common_subgraphs` currently enumerates every admissible
partial node correspondence by building the modular product and enumerating
all of its cliques. `ReactionAst::compose` is its principal consumer and
requires the complete result, including the empty correspondence.

This document tracks two changes:

1. rename the underspecified
   `CommonSubgraphEnumerationAlgorithm::Backtracking` variant to
   `ModularProductBacktracking`;
2. add `DirectBacktracking`, which enumerates the same correspondences directly
   without materializing the modular product.

Both algorithms retain the existing eager `Vec<GraphCorrespondence>` API.
Neither changes the definitions of `EmbeddingKind::Induced` and
`EmbeddingKind::Monomorphism`, the treatment of node and edge predicates, or
the deterministic ordering of the returned set.

The Python names follow the Rust names:

```python
CommonSubgraphEnumerationAlgorithm.ModularProductBacktracking()
CommonSubgraphEnumerationAlgorithm.DirectBacktracking()
```

The visible default of `ReactionAst.compose` remains the existing algorithm,
renamed to `ModularProductBacktracking()`. Adding the direct implementation
does not itself justify a default change.

## Existing modular-product algorithm

For every node pair `(u, v)` accepted by `node_match`, the modular product has
one vertex. Two product vertices are adjacent when they use distinct nodes and
their source-graph relationship is compatible under the selected embedding:

- if both source edges exist, `edge_match` must accept them;
- if neither source edge exists, the pairs are compatible;
- if exactly one source edge exists, the pairs are compatible for
  `Monomorphism` and incompatible for `Induced`.

Every clique is therefore one admissible partial injective node
correspondence. The implementation enumerates all cliques, constructs their
edge correspondences, and sorts and deduplicates the result.

`ModularProductBacktracking` accurately names both parts of this algorithm:
the intermediate graph and the clique walk. The current `Backtracking` name
distinguishes neither from the direct alternative.

If `P` node pairs pass `node_match`, the product adjacency matrix retains
`O(P²)` bits before output materialization. This can dominate memory even when
structural incompatibility later excludes most mappings.

## Direct algorithm

`DirectBacktracking` searches partial injective node mappings without an
intermediate product graph:

1. Visit left-graph nodes in ascending `NodeId` order.
2. For the current left node, branch over every compatible unused right node
   and over leaving the left node unmapped.
3. Before extending with `(u, v)`, compare it with every already selected
   `(u', v')` using exactly the same edge-presence and `edge_match` rules as
   the modular product.
4. After every left node has been considered, emit the node correspondence
   and the correspondence of source edges present on both sides.

Node candidates are computed once from `node_match`. Search state consists of
the current sorted mating, the used-right-node bitset, and the recursion
position. The implementation applies `edge_match` while testing extensions;
predicate-result caching is an optimization to consider only if benchmarks
show repeated expensive edge predicates to be material.

Every partial injective mapping has exactly one search path because each left
node is visited once and each right node can be selected once. The empty
mapping is the leaf reached by taking every unmapped branch. The implementation
still sorts the final vector by node mating so its public ordering is identical
to `ModularProductBacktracking`.

The direct representation removes the product's quadratic retained adjacency
matrix. It does not change the exponential number of valid correspondences or
the eager cost of returning all of them.

## Staged implementation plan

### S0 — Baseline and algorithm name

- **S0a — common-subgraph benchmark baseline**
  (`umol-graph-core/benches/algorithms.rs`): add a
  `common_subgraph_enumeration` Criterion group for complete enumeration under
  both embedding kinds. Include small dense-compatible, structurally
  selective, label-selective, and molecular-graph pairs; pin each case's output
  count outside the timed section so comparisons cannot benchmark incomplete
  work. Record the existing modular-product implementation before adding the
  alternative. **Additive (green).** `[dep: —]`
  **Implemented (green).** The benchmark group covers dense-compatible,
  structurally selective induced and monomorphic, label-selective, and
  benzene--naphthalene cases. Their pinned complete-output counts are 34, 69,
  209, 109, and 1,957. A short baseline run measured approximately 6.18,
  13.9, 44.0, 23.0, and 613 µs respectively; these measurements are comparison
  observations, not acceptance thresholds.
- **S0b — coordinated selector rename**
  (`umol-graph-core/src/algorithms/common_subgraph/enumeration.rs`,
  `umol-ast/src/ast/compose.rs`, `umol-py/src/algorithm.rs`,
  `umol-py/src/reaction.rs`, Rust and Python callers): rename
  `Backtracking` to `ModularProductBacktracking` throughout Rust and Python.
  Update the Python repr, conversion, equality tests, exports, and the visible
  `ReactionAst.compose` default. Exact common-subgraph and composition results
  must remain unchanged. **Breaking Rust and Python naming migration
  (red→green).** `[dep: S0a]`
  **Implemented (green).** The Rust and Python variants are now
  `ModularProductBacktracking` and `ModularProductBacktracking()`. The Python
  repr and visible `ReactionAst.compose` default use the same name; all Rust,
  Python, property-test, benchmark, and discussion-doc callers were migrated
  without changing enumeration or composition results.

### S1 — Direct enumeration

- **S1a — direct search kernel**
  (`umol-graph-core/src/algorithms/common_subgraph/enumeration.rs`): implement
  and test the module-local direct partial-mapping search before extending the
  public enum. Exact tests call the kernel directly and cover empty graphs,
  isolated nodes, incompatible node and edge labels, injectivity, disconnected
  mappings, the empty correspondence, and the induced/monomorphism
  edge-presence distinction. Its sorted output must equal the existing public
  modular-product result for the same fixtures. **Additive (green).**
  `[dep: S0b]`
  **Implemented (green).** The `direct_backtracking` kernel precomputes node
  candidates, walks left nodes in id order, tracks used right nodes with a
  bitset, applies the embedding and edge-predicate checks while extending, and
  constructs sorted complete `GraphCorrespondence` output at the leaves. Both
  kernels use the same node-to-edge correspondence materializer, preventing
  their edge results from drifting. Exact and cross-algorithm tests cover all
  listed cases.
- **S1b — coordinated public selector**
  (`umol-graph-core/src/algorithms/common_subgraph/enumeration.rs`,
  `umol-py/src/algorithm.rs`, `umol-py/tests/test_algorithm.py`,
  `umol-py/tests/test_reaction.py`): add
  `CommonSubgraphEnumerationAlgorithm::DirectBacktracking`, dispatch to the S1a
  kernel, and migrate every exhaustive Rust and Python match in the same
  subitem. Bind `DirectBacktracking()` with exact Rust conversion, equality,
  and repr behavior. Exercise reaction composition explicitly with both
  selectors and assert the same complete, deterministically ordered
  composites. Keep `ModularProductBacktracking()` as the Python default.
  **Additive Rust and Python API with exhaustive-enum caller migration
  (red→green).** `[dep: S1a]`
  **Implemented (green).** `DirectBacktracking` is a public graph-core selector
  and dispatches to the direct kernel without constructing the modular
  product. Python exposes `DirectBacktracking()` with exact conversion,
  equality, and repr behavior. Both selectors produce the same exact ordered
  composites across the Python composition cases, while
  `ModularProductBacktracking()` remains the visible Python default.
- **S1c — generated cross-validation and comparative benchmarks**
  (`umol-graph-core/tests/property/common_subgraph.rs`,
  `umol-graph-core/tests/property.rs`,
  `umol-graph-core/benches/algorithms.rs`): for generated labeled graph pairs,
  compare the complete sorted `Vec<GraphCorrespondence>` from the two
  algorithms under both embedding kinds. Include asymmetric graph sizes,
  empty graphs, and predicates that reject node or edge pairs. Extend the S0
  benchmark group with `DirectBacktracking` and report both runtime and output
  count on the unchanged cases. **Additive verification (green).**
  `[dep: S1b]`
  **Implemented (green).** Generated labeled simple-graph pairs independently
  vary from empty through four left and five right nodes. Equality-label
  predicates exercise rejected node and edge pairs; for both embedding kinds,
  256 generated cases produced identical complete, sorted
  `Vec<GraphCorrespondence>` values from the two algorithms. The unchanged
  benchmark cases retain output counts 34, 69, 209, 109, and 1,957. A short
  paired run measured modular-product/direct times of approximately
  6.20/2.82, 14.1/6.09, 44.6/19.3, 23.4/10.0, and 605/364 µs respectively.
  These measurements are comparison observations, not acceptance thresholds.

At the end of S0 and S1, the workspace test and lint suites must be green. The
critical path is `S0a → S0b → S1a → S1b → S1c`; no stage in this plan is
deferrable.

## Deferred ZDD representation

A zero-suppressed decision diagram is a possible representation for a future
*restricted* composition operation. Compatible node pairs can serve as the
decision variables, while the represented family contains the injective,
embedding-compatible pair sets. Such a representation could support:

- restricting mapped-node cardinality without expanding every smaller
  correspondence;
- enumerating selected cardinality layers;
- applying composition-specific mandatory or forbidden-pair constraints;
- retaining a large overlap family compactly when the caller does not require
  every composite immediately.

This is not useful to the current eager complete-composition contract if every
represented correspondence is immediately expanded into a `ReactionAst`.
Ordering by mapped-node count would also require cardinality analysis of the
ZDD rather than relying on its variable order. If overlap size includes shared
edges rather than only selected node pairs, it is not ordinary set cardinality
and requires additional state.

No ZDD dependency, public type, or algorithm variant is part of S0 or S1. It
should be reconsidered only together with a restricted or lazy composition API
whose filters can be applied before materialization.

## Maximum and maximal common subgraphs

Raised 2026-07-31. Most callers who reach for common-subgraph enumeration want the **maximum** common
subgraph (largest node count) or a **maximal** one (cannot be extended), rather than the complete
family filtered by a criterion they have to supply. Both are reasonable operations to offer. Neither
may be plumbed into composition.

### Why composition cannot use them

`ReactionAst::compose` requires the complete result **including the empty correspondence**, and the
empty correspondence is never maximal. It is the case where the two reactions act on disjoint parts
of the host — two independent transformations — which is the dominant case in reaction-network
generation, not an edge case. Restricting to maximal or maximum overlaps therefore discards the most
common situation first, because a zero-size overlap is minimal by construction.

The composition law degrades in a specific way that is worth stating precisely rather than as "not
necessarily equivalent". With a single maximum overlap the composite explores one identification
where sequential application explores all of them, so

```text
composed ⊆ sequential      (containment, not equality)
```

Equality holds only for the complete enumeration. Any restricted variant must document the
containment direction rather than claim the law.

### Maximal is a different enumerator, not a filter

Maximal common subgraphs are exactly the **maximal cliques** of the modular product, which
Bron–Kerbosch produces natively. So this is a change of clique walk, not a post-filter over the
existing output: extension stops instead of every smaller correspondence being materialised and then
discarded. That is the efficiency the request is really after, and it needs no compressed
representation to obtain.

### Maximum is branch-and-bound

Enumerate-all-then-take-largest is the naive route and is precisely what should be avoided. The
established approach is branch and bound with bounds on the achievable remainder — the McSplit family
and its predecessors. This belongs as its own algorithm variant, not as a criterion applied to
enumeration output.

### Consequence for the deferred ZDD

Neither operation supplies the "restricted or lazy composition API" that the ZDD section above is
waiting for. Maximum wants pruning, not a compact family; maximal wants a different walk. The case
that would justify a ZDD is a **cardinality-threshold** query — overlaps of at least `k` mapped
nodes — or repeated queries against one retained family. That is a third consumer, distinct from both
composition and maximum common subgraph, and should not be conflated with this request.

### Connectedness is the axis callers ask for first

Cheminformatics maximum common subgraph almost always means a single **connected** fragment, while
the current enumeration yields disconnected correspondences as well. Raymond and Willett's MCIS and
MCES distinction maps approximately onto the existing `EmbeddingKind::Induced` and
`EmbeddingKind::Monomorphism`, so that axis is already covered; connected versus disconnected is not,
and it is the parameter users reach for before they reach for maximum versus maximal.

### Proposed shape

Three named operations with visibly different contracts, rather than one call with a mode parameter:

| operation | contract |
| --- | --- |
| `enumerate_common_subgraphs` | complete, including the empty correspondence; the only one composition may use |
| `maximal_common_subgraphs` | maximal cliques; no completeness guarantee |
| `maximum_common_subgraph` | one largest; branch and bound |

The completeness guarantee is what composition depends on, so it should be legible at the call site
rather than hidden in an argument.

## References

- Minato, [*Zero-Suppressed BDDs for Set Manipulation in Combinatorial
  Problems*](https://doi.org/10.1109/DAC.1993.203958), 1993.
- Kawahara et al., [*Frontier-Based Search for Enumerating All Constrained
  Subgraphs with Compressed Representation*](https://doi.org/10.4230/LIPIcs.SEA.2020.9),
  2020.
- Raymond and Willett, [*Maximum common subgraph isomorphism algorithms for the
  matching of chemical structures*](https://doi.org/10.1023/A:1021271615909), 2002.
- Bron and Kerbosch, [*Algorithm 457: finding all cliques of an undirected
  graph*](https://doi.org/10.1145/362342.362367), 1973.
- McCreesh, Prosser and Trimble, [*A Partitioning Algorithm for Maximum Common Subgraph
  Problems*](https://doi.org/10.24963/ijcai.2017/99), 2017.
