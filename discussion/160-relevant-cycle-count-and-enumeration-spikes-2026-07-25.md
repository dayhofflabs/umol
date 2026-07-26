# Relevant-cycle counting and enumeration alternatives

Status: **Active**
Date: 2026-07-25
Relates: [105](105-dsl-fixes-2026-06-06.md),
[158](158-ring-model-and-enumeration-2026-07-22.md),
[159](159-simple-graph-policy-2026-07-23.md)

## Scope

Doc 158 established the cycle-set semantics and implemented bounded
Read--Tarjan simple-cycle enumeration, compact Vismara relevant-cycle analysis,
Horton minimum cycle bases, and Kolodzik Unique Ring Families. This document
tracks three follow-ups:

1. expose the exact number of bounded relevant cycles without materializing
   them;
2. evaluate the Hanser path-graph reduction with the May--Steinbeck
   implementation refinements as a second simple-cycle enumerator;
3. evaluate the Birmelé et al. output-optimal enumerator as a more ambitious
   alternative.

The count operation is required work. The two enumerators are benchmark and
design spikes. Neither becomes a public algorithm value merely because it has
been prototyped.

This work does not add another `RingSetKind`, change the definition of simple
or relevant cycles, or reopen essential cycles, SSSR compatibility, and
triplet-short cycles.

## Existing implementation

The public graph-core cycle surface currently selects:

- `SimpleCycleEnumerationAlgorithm::ReadTarjan`;
- `RelevantCycleEnumerationAlgorithm::Vismara`;
- `MinimumCycleBasisAlgorithm::Horton`;
- `UniqueRingFamilyAlgorithm::Kolodzik`.

`RelevantCycleAnalysis` already stores compact Vismara cycle families. Every
family has one cycle weight and an exact count obtained from the product of the
two compatible shortest-path counts. `UniqueRingFamilies` sums those values
into the existing public `RelevantCycleCount(BigUint)` type.

The missing operation is therefore not a second counting algorithm. It is a
public projection of information already computed by Vismara before the
families are expanded into individual cycles.

## Relevant-cycle count contract

For a graph `graph`, source-edge bound `max_cycle_size`, and relevant-cycle
algorithm `algorithm`, the count operation returns the exact number of cycles
that a complete call to `visit_relevant_cycles` would emit with the same
arguments.

Consequently:

- the result uses `RelevantCycleCount`, not `usize`;
- a family contributes either its complete count or zero because every cycle
  in a Vismara family has the same weight;
- the bound is measured in source-graph edges;
- self-loops contribute one cycle each when the bound is at least one;
- distinct parallel-edge cycles remain distinct;
- subdivision is an implementation detail and does not change the source-edge
  count or cycle identity;
- the count is zero for an acyclic graph or when the bound excludes every
  cycle;
- increasing the bound cannot decrease the count;
- counting must not construct any individual `Cycle`.

The direct, fallback, and combined operations must follow the same split as
relevant-cycle visitation:

- the direct operation checks graph simplicity and returns
  `NonSimpleGraphError`;
- the fallback operation counts loops and uses subdivision for the loopless
  remainder;
- the combined operation chooses the direct or fallback path after the same
  simplicity check used by visitation.

The likely names are:

```rust
impl Graph {
    pub fn try_count_relevant_cycles(
        &self,
        max_cycle_size: usize,
        algorithm: RelevantCycleEnumerationAlgorithm,
    ) -> Result<RelevantCycleCount, NonSimpleGraphError>;

    pub fn count_relevant_cycles_fallback(
        &self,
        max_cycle_size: usize,
        algorithm: RelevantCycleEnumerationAlgorithm,
    ) -> RelevantCycleCount;

    pub fn count_relevant_cycles(
        &self,
        max_cycle_size: usize,
        algorithm: RelevantCycleEnumerationAlgorithm,
    ) -> RelevantCycleCount;
}
```

These names follow the existing `try_visit_*`, `visit_*_fallback`, and
`visit_*` operation families. They remain proposed names until explicitly
approved.

The `RelevantCycleCount` documentation must be widened: it is the exact number
of relevant cycles represented by compact cycle-family state, not a value
specific to one Unique Ring Family.

### Verification

The count implementation must establish:

- equality with the number of visitor invocations for every manageable graph
  in the checked-in graph corpus;
- equality with collected relevant-cycle length for the literature fixtures;
- exact loop, parallel-edge, disconnected, acyclic, and bound-zero behavior;
- equality among direct and combined results on simple graphs;
- equality among fallback and combined results on non-simple graphs;
- monotonicity under increasing `max_cycle_size`;
- counts exceeding `usize::MAX` without individual-cycle emission.

The last property should use compact path-family construction rather than
attempting to enumerate the represented cycles.

## Hanser--May spike

Hanser, Jauffret, and Kaufmann enumerate all elementary cycles by progressively
removing vertices from a path graph. Each retained path-graph edge represents a
simple source-graph path. Combining compatible paths through a removed vertex
creates new path edges, while a path-graph loop represents one completed source
cycle.

May and Steinbeck retain this algorithm but improve its practical
implementation through incidence-list storage, binary path sets, removal
ordering, and a maximum intermediate-degree feasibility measure. The resulting
algorithm is a credible chemical-graph implementation rather than only a
historical alternative.

The spike must answer:

1. Can the May--Steinbeck representation be expressed as a visitor without
   first collecting all cycles?
2. Can `max_cycle_size` prune path construction directly and exactly?
3. How do removal ordering and intermediate path-graph degree affect runtime
   and peak memory on fused, bridged, cage, and macrocyclic graphs?
4. Does it provide a repeatable performance region in which it materially
   improves on Read--Tarjan?
5. Can a promoted implementation preserve graph-core's loop and parallel-edge
   semantics through an edge-aware representation or subdivision fallback?

The feasibility threshold is instrumentation during the spike, not part of the
cycle-set semantics. The spike must never silently return a partial set. If a
later public implementation needs an abort threshold, that requires a separate
explicit incomplete-result or error contract.

The spike starts with simple graphs. Promotion requires a defined non-simple
path consistent with doc 159 and source-edge bounds.

## Birmelé et al. spike

Birmelé et al. reduce undirected cycle listing to optimal listing of simple
`s`--`t` paths through a sequence of biconnected components. Their unbounded
algorithm runs in

```text
O(E + sum(length(cycle) for cycle in cycles))
```

which is the input size plus the size of the emitted output.

This is a stronger asymptotic target than Read--Tarjan and the Hanser path
graph, but the paper's result is for unbounded enumeration of simple undirected
graphs. The graph-core operation has two additional requirements:

- apply `max_cycle_size` during traversal rather than enumerate and discard
  longer cycles;
- preserve the visitor's early-break behavior.

The spike must therefore separate two questions:

1. Can the paper's unbounded algorithm be implemented faithfully with exact
   output parity and the stated amortized organization?
2. Can a depth-bounded adaptation prune the search without invalidating
   correctness, and what complexity and practical performance does that
   adaptation retain?

The unbounded optimality result must not be claimed for the bounded adaptation
without a corresponding argument. Promotion also requires a multigraph path,
most likely loop extraction followed by subdivision and projection, with bounds
translated back into source-edge units.

## Comparative evaluation

The two spikes must share one harness and compare against the completed
Read--Tarjan implementation. Correctness and performance evidence are separate.

Correctness evaluation uses:

- exhaustive small-graph fixtures and generated graphs;
- the checked-in cycle-family corpus;
- the literature graphs in `tests/property/cycles/literature.rs`;
- exact normalized edge-cycle sets at several length bounds;
- visitor/collector equality and early-break behavior.

Performance evaluation must include data independent of the correctness corpus:

- representative molecular structures sampled from a real chemical dataset;
- parameterized fused-ring chains, bridged systems, cages, macrocycles, and
  high-cycle-count structures;
- sparse nonchemical graphs that expose output-sensitive behavior;
- several source-edge bounds, including common chemical bounds and effectively
  unbounded enumeration.

Measure:

- total enumeration time;
- time to first emitted cycle;
- time for early termination after fixed output counts;
- peak retained memory and allocation volume;
- preprocessing time;
- emitted cycle count and total emitted edge count;
- Hanser path-graph maximum intermediate degree;
- scaling with graph size, cycle count, and total emitted edge count.

No external chemistry or graph package becomes a benchmark or test dependency.
One-off external comparisons may be captured as checked-in evidence under the
same policy as doc 158.

## Promotion criteria

An experimental enumerator is eligible for the public algorithm selector only
when:

- it returns exactly the existing simple-cycle set at every tested bound;
- it applies the bound during traversal;
- it supports visitor early termination;
- its non-simple behavior is defined without silently selecting a different
  algorithm;
- its implementation and algorithm citation are maintainable;
- benchmarks demonstrate a material and repeatable advantage over Read--Tarjan
  for an identified workload.

Promotion does not imply replacement. Multiple public implementations are
justified only when their performance regions are meaningfully different. The
public enum variant names require approval after the spikes identify which
algorithms, if any, deserve promotion.

## References

- Hanser, Jauffret, and Kaufmann, [*A New Algorithm for Exhaustive Ring
  Perception in a Molecular Graph*](https://doi.org/10.1021/ci960322f), 1996.
- May and Steinbeck, [*Efficient Ring Perception for the Chemistry Development
  Kit*](https://doi.org/10.1186/1758-2946-6-3), 2014.
- Vismara, [*Union of All the Minimum Cycle Bases of a
  Graph*](https://doi.org/10.37236/1294), 1997.
- Birmelé et al., [*Optimal Listing of Cycles and st-Paths in Undirected
  Graphs*](https://doi.org/10.1137/1.9781611973105.134), 2013.
