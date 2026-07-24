# Ring model and cycle enumeration

Status: **Active**
Date: 2026-07-22
Relates: [057](057-sssr-needed-2026-03-11.md),
[058](058-aromaticity-perception-2026-03-11.md),
[059](059-aromaticity-perception-review-2026-03-12.md),
[086](086-molecule-ast-api-2026-04-16.md),
[149](149-molecule-ring-cache-and-hashing-2026-07-13.md),
[151](151-python-molecule-workflows-2026-07-13.md),
[159](159-simple-graph-policy-2026-07-23.md)

## Scope

The algorithm-transparency work in doc 151 made the cycle-enumeration backend
explicit, but exposed a deeper problem: the current API does not cleanly
separate the requested ring family from the algorithm used to enumerate it.
This document defines that separation before ring configuration is threaded
through fingerprints, aromaticity, and Python.

The design covers:

- the mathematical contract of each ring family;
- the algorithms capable of enumerating each family;
- edge-aware cycle semantics for loops and parallel edges;
- subdivision as the internal adapter for simple-graph algorithms;
- the division between `RingModel` and `RingConfig`;
- bounded enumeration and incomplete-result behavior;
- the current atom-filter facility.

The staged implementation plan follows the settled design.

## Prior decisions

This design does not reopen the general SSSR question.

- Doc 057 rejects a non-unique minimum cycle basis as the default chemical ring
  description. It identifies relevant cycles, essential cycles, and local
  shortest-ring queries as the meaningful alternatives.
- Doc 058 separates exhaustive all-cycle enumeration from relevant-cycle
  perception and from pi-subgraph aromaticity. It also records the CDK split
  between all, relevant, and minimum-cycle-basis finders.
- Doc 059 recommends Johnson enumeration for the then-proposed local pi-subgraph
  strategy. The present design also covers bounded global enumeration on an
  undirected molecular graph, for which an undirected algorithm should be
  considered directly.
- Doc 086 selects bounded Vismara relevant cycles as the ordinary global ring
  set and explicitly accepts differences from SSSR-based SMARTS engines in
  symmetric and fused graphs.

The design therefore determines which additional cycle-set models have an
actual consumer, which lower-level algorithms implement them, and whether
RDKit compatibility justifies exposing RDKit-specific ring constructions.

## Topological cycle model

Ring enumeration operates on the unweighted molecular topology. An edge is
either present or absent; its bond order, aromatic status, direction-like
annotation, and other edge properties do not affect cycle membership, cycle
length, ring-family membership, or cycle identity.

Consequently:

- cycle length is the number of topological edges;
- cycle identity includes edge identity, not only the traversed vertices;
- changing an existing bond's properties does not change the `RingSet`;
- algorithms must not accept edge weights or edge-label predicates as part of
  ring enumeration;
- algorithms selected for the same family must return the same topological
  cycles regardless of bond attributes.

Chemical operations may inspect bond properties after enumeration. Such
classification belongs to the consuming model, not to `RingModel` or the cycle
enumerator.

Graph-core cycle operations are defined on the full storage semantics of
`Graph`, including self-loops and parallel edges. A self-loop is a one-edge
cycle. Two parallel edges form a two-edge cycle, and cycles traversing different
parallel edges remain distinct.

The molecule layer has a narrower chemical contract. `Ring` retains its
existing minimum of three distinct atoms, so one-edge and two-edge graph cycles
are not chemical rings. A molecule containing a localized bond self-loop or
parallel localized bonds remains structurally invalid and is reported by
`EntityStructureValidator`; deriving a ring view does not become another
validation operation. For valid molecules the graph-core and chemical domains
coincide.

## Current API does not express its behavior

`CycleEnumerationAlgorithm` currently contains only `Vismara`, and
`Graph::enumerate_cycles` therefore enumerates relevant cycles. `RingSet` then
interprets the same output in two ways:

- `RingSetKind::Simple` accepts the Vismara result unchanged;
- `RingSetKind::Relevant` filters that result for induced cycles.

Neither branch implements its documented distinction. `Simple` does not return
all simple cycles or one minimum cycle basis. `Relevant` receives relevant
cycles before applying its extra filter.

The current return type, `Vec<Vec<NodeId>>`, also cannot express graph-core's
stored topology. It cannot distinguish two cycles that traverse different
parallel edges, and the current minimum length of three silently omits
self-loop and parallel-edge cycles. These are representation and algorithm
coverage defects, not reasons to make cycle queries fallible.

Vismara defines the relevant cycles as the union of all minimum cycle bases,
equivalently as the cycles that cannot be generated from strictly shorter
cycles. The family is a graph property; Vismara is one algorithm for computing
it. See [Vismara, *Union of all the minimum cycle bases of a graph*
(1997)](https://www.lirmm.fr/~vismara/papers/vismara97.pdf).

Under this unweighted topological model, every relevant cycle is chordless. A
chord splits a cycle into two strictly shorter cycles, whose symmetric
difference is the original cycle. The induced-cycle filter is therefore
redundant. This implication need not hold for arbitrary positive edge weights,
but weighted relevance is outside the ring model.

The existing claim in doc 086 that relevant cycles need not be induced applies
too broadly and must not define the new API.

## Cycle-set families

The family is a semantic choice. The algorithm is an operational choice that
must satisfy the selected family's contract.

| Family | Contract | Recommendation |
| --- | --- | --- |
| Relevant | Union of all minimum cycle bases, bounded by the requested maximum cycle size. | Required; ordinary global ring model. |
| Simple | Every elementary cycle, bounded by the requested maximum cycle size; rotations and reversal identify the same undirected cycle. | Required as an explicit exhaustive model, not as the default. |
| Essential | Intersection of all minimum cycle bases. Unique, polynomial in number, but possibly empty for a cyclic graph. | Mathematically valid but deferred until a consumer needs this exact set. |

At graph-core, these definitions include one-edge and two-edge cycles. The
relevance of longer cycles is computed in the full multigraph cycle space
before the molecule layer removes cycles with fewer than three distinct atoms.
The molecule layer must not first delete or collapse invalid edges and then
recompute a different cycle family.

A single minimum cycle basis is not another cycle-set model for chemical ring
perception. It is one selected basis, can depend on tie-breaking, and has a
different result contract. If it is needed for linear-algebraic work, it should
be a separate basis operation with its own algorithm selector. It must not be
called SSSR unless a compatibility operation deliberately reproduces a
particular SSSR implementation.

Chordless cycles become a separate family only if a consumer requires them
independently of relevant cycles. Under the present unweighted model, adding a
chordless family would not distinguish the current relevant cycles.

The maximum ring size is part of `RingModel` under the present contract because
it changes which rings belong to the returned set. A computational result cap,
time limit, or memory limit belongs in `RingConfig`, but exhausting such a limit
must return an explicit incomplete result or error rather than silently alter
the family.

## Ring concepts that are not `RingSetKind` variants

Several useful ring operations do not return another interchangeable set of
cycles and should not be forced into the family selector.

| Concept | Appropriate contract |
| --- | --- |
| Minimum cycle basis | Separate basis operation returning one linearly independent basis. |
| Smallest cycle through an atom or bond | Local shortest-cycle query. Graph-core already exposes these operations. |
| Ring membership without enumeration | Bridge or biconnected-component query; equivalent in purpose to RDKit `FastFindRings`, not a cycle family. |
| Shortest cycle through each vertex or edge | Batch form of the local queries; add only if repeated consumers require it. |
| Shortest cycle through each vertex triple | Specialized set used by CDK for CACTVS/PubChem keys; add with that fingerprint rather than to obtain generic toolkit parity. |
| Unique ring families | Decomposition of relevant cycles into equivalence classes, with family-level atom/bond membership and relevant-cycle counts. Its result is a collection of families, not another flat cycle set. |

`RingSetKind` keeps this cycle-set selection distinct from the Unique Ring
Family concept. Unique Ring Families have their own result type and operation
rather than becoming `RingSetKind::Unique`.

## RDKit parity

RDKit exposes four materially different ring facilities:

- `GetSSSR` computes one potentially non-unique SSSR using a modified Figueras
  algorithm;
- `GetSymmSSSR` adds rings to produce a symmetric small set of smallest rings;
- `FastFindRings` marks cyclic atoms and bonds without constructing an SSSR;
- `FindRingFamilies` computes Unique Ring Families through
  RingDecomposerLib. `RingInfo` then exposes ring-family membership and the
  number of relevant cycles.

References: [RDKit ring-finding documentation](https://rdkit.org/docs/RDKit_Book.html#ring-finding-and-sssr),
[RDKit `findSSSR` API](https://www.rdkit.org/docs/cppapi/namespaceRDKit_1_1MolOps.html),
[RDKit `FindRingFamilies` API](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html),
and [RDKit `RingInfo`](https://rdkit.org/docs/cppapi/RingInfo_8h_source.html).

Parity should be divided into two policies:

1. **Semantic parity is desirable.** umol should answer invariant questions
   that RDKit users rely on: whether an atom or bond is cyclic, its smallest
   containing cycle size, relevant-cycle enumeration, and—if a consumer needs
   it—Unique Ring Family decomposition. These operations need not reproduce
   RDKit's internal SSSR.
2. **Exact SSSR-output parity is not a default goal.** `GetSSSR` is
   non-canonical, while `GetSymmSSSR` is an RDKit-specific construction rather
   than a general mathematical family. Reproducing either output requires the
   same algorithm, tie-breaking, and ordering conventions, not merely another
   minimum-cycle-basis algorithm.

Exact RDKit parity becomes justified only for a named compatibility consumer,
such as reproducing an RDKit descriptor, SMARTS ring-count predicate, or
aromaticity model whose result demonstrably depends on the symmetrized SSSR. In
that case it belongs in an explicit RDKit compatibility model and is tested
against versioned RDKit reference data. It must not redefine the ordinary
`RingModel`.

The useful RDKit parity additions are therefore:

| Priority | Addition | Placement |
| --- | --- | --- |
| Required | Cyclic atom/bond and smallest-cycle queries | Existing graph-core local operations and ring views; no new family. |
| Required | Relevant and bounded all-simple cycle sets | `RingModel` plus the family-specific enumeration algorithm. |
| Candidate | Unique Ring Family decomposition | Separate operation/result, justified by ring-topology, descriptor, or conformer consumers. |
| Conditional | Exact symmetrized-SSSR output | Explicit RDKit compatibility model only. |
| Not planned | Arbitrary SSSR/MCB as the ordinary ring set | Separate basis operation only if linear-algebraic work needs it. |

The corresponding algorithm audit is:

| Algorithm or implementation | Reason to add |
| --- | --- |
| Vismara relevant-cycle enumeration | Already required by the ordinary umol ring model; not added merely for RDKit parity. |
| Read--Tarjan or another bounded all-simple enumerator | Required by the explicit exhaustive model; RDKit does not require this parity. |
| Unique Ring Family decomposition | The one new algorithmic facility supported by both a canonical graph contract and a current RDKit operation. Add when the separate decomposition API has a consumer. |
| Modified Figueras plus RDKit symmetrization | Add only inside an exact RDKit compatibility implementation. Another MCB algorithm does not provide parity. |
| DFS ring marking | No cycle-enumeration addition. Use the existing topology operations for cyclic membership. |

RDKit's `includeDativeBonds` and `includeHydrogenBonds` switches on
`FindRingFamilies` are not part of the parity target. In umol, the topology
presented to ring perception determines which edges exist; bond attributes do
not opt edges into or out of the ring model.

Unique Ring Families are the strongest parity candidate beyond the two cycle
sets. They extend relevant cycles, are invariant under atom ordering, and are
polynomial in number even when the relevant-cycle set is exponential. See
[Kolodzik, Urbaczek, and Rarey, *Unique Ring Families: A Chemically Meaningful
Description of Molecular Ring Topologies*
(2012)](https://doi.org/10.1021/ci200629w) and the
[RingDecomposerLib API](https://ringdecomposerlib.readthedocs.io/en/latest/RingDecomposerLib_8h.html).

## All-simple-cycle enumeration

The required operation is complete enumeration of every undirected simple
cycle whose length does not exceed the model's maximum. It must emit each cycle
once, treating rotation and reversal as the same cycle, and must apply the
length bound during enumeration rather than enumerate longer cycles and discard
them afterward.

The leading initial algorithm is the undirected Read--Tarjan backtracking
algorithm. It was designed for listing simple cycles in undirected graphs, uses
linear space, and has `O(V + E + E N)` total time for `N` reported objects. See
[Read and Tarjan, *Bounds on Backtrack Algorithms for Listing Cycles, Paths,
and Spanning Trees* (1975)](https://doi.org/10.1002/net.1975.5.3.237). Its path
search admits a direct maximum-depth bound, but the bounded adaptation must be
specified and tested as part of the implementation rather than assumed from
the unbounded proof.

The paper states its algorithm for graphs without loops or multiple edges, but
also states that it is readily modified to support them. The implementation
should make paths edge-aware rather than routing non-simple input through the
subdivision fallback:

- distinguish every traversal by `EdgeId`;
- emit each self-loop directly as a one-edge cycle;
- permit two distinct parallel edges to close a two-edge cycle;
- continue to forbid repeated internal vertices for longer elementary cycles;
- normalize by the ordered node-and-edge traversal, so parallel-edge cycles are
  not deduplicated by their node sequence.

This direct extension keeps bounded all-cycle enumeration in the original graph
and avoids doubling its search depth.

Two alternatives do not change the public model:

- Johnson's algorithm is the established general baseline, with
  `O((V + E)(N + 1))` time and linear space, but is formulated for directed
  graphs. An undirected implementation needs biconnected-component reduction,
  duplicate-orientation handling, and a correct bounded-search adaptation. See
  [Johnson, *Finding All the Elementary Circuits of a Directed Graph*
  (1975)](https://doi.org/10.1137/0204007).
- Birmelé et al. give an asymptotically optimal undirected enumerator whose
  running time is the input size plus the total size of the emitted cycles. It
  is a stronger but more involved implementation candidate if benchmarks show
  Read--Tarjan's per-cycle edge factor to matter. See [Birmelé et al.,
  *Optimal Listing of Cycles and st-Paths in Undirected Graphs*
  (2013)](https://doi.org/10.1137/1.9781611973105.134).

The Gupta--Suzumura bounded directed-cycle algorithm is not a suitable basis
for the first implementation. A later analysis gives counterexamples in which
it omits valid bounded cycles; adopting it would require independently applying
and validating the correction. See [Bauernöppel and Sack, *Finding All
Bounded-Length Simple Cycles in Directed Graphs -- Revisited*
(2025)](https://arxiv.org/abs/2512.08392).

The family-specific graph-core selectors are:

```rust
pub enum SimpleCycleEnumerationAlgorithm {
    ReadTarjan,
}

pub enum RelevantCycleEnumerationAlgorithm {
    Vismara,
}
```

The public operations are likewise separated:

```rust
impl Graph {
    pub fn visit_simple_cycles<B>(
        &self,
        max_cycle_size: usize,
        algorithm: SimpleCycleEnumerationAlgorithm,
        visitor: impl FnMut(Cycle) -> ControlFlow<B>,
    ) -> ControlFlow<B>;

    pub fn enumerate_simple_cycles(
        &self,
        max_cycle_size: usize,
        algorithm: SimpleCycleEnumerationAlgorithm,
    ) -> Vec<Cycle>;

    pub fn visit_relevant_cycles<B>(
        &self,
        max_cycle_size: usize,
        algorithm: RelevantCycleEnumerationAlgorithm,
        visitor: impl FnMut(Cycle) -> ControlFlow<B>,
    ) -> ControlFlow<B>;

    pub fn enumerate_relevant_cycles(
        &self,
        max_cycle_size: usize,
        algorithm: RelevantCycleEnumerationAlgorithm,
    ) -> Vec<Cycle>;
}
```

The operation name determines the cycle-set semantics, and its enum selects an
implementation of those semantics. Invalid family/algorithm combinations are
therefore unrepresentable, and no algorithm-selection error is required.

The related simple-path and spanning-tree algorithms from the same paper are
not part of this work. The simple-path algorithm is a direct adaptation of the
cycle backtracker, while the spanning-tree algorithm uses separate partial-tree
state and bridge/cycle restrictions. Neither operation is scheduled here.

## Subdivision fallback for simple-graph algorithms

Neither Vismara nor Read--Tarjan defines a runtime failure that identifies
non-simple input. Both are stated for graphs without loops or multiple edges,
and an unmodified implementation may silently omit or conflate cycles rather
than return a recognizable failure. The public cycle APIs therefore cannot use
“try the simple algorithm, fall back if it fails” dispatch.

Vismara analysis already requires an initial graph traversal and biconnected
decomposition. That preprocessing should detect a self-loop or repeated
endpoint while it scans adjacency:

```text
initial Vismara / RCF analysis
    |
    +-- no non-simple edge encountered -> continue on the original graph
    |
    `-- non-simple edge encountered
            -> discard the partial linear-time analysis
            -> use subdivision fallback
```

This is optimistic dispatch without a separate ahead-of-time validation pass.
The detection is internal, produces no public error, and at worst discards part
of one `O(V + E)` preprocessing traversal.

The fallback uses the subdivision graph, also called the barycentric
subdivision. Replacing every non-loop edge by a path of length two creates a
bipartite graph whose vertices are the original vertices and one inserted
vertex per original edge. For a loopless multigraph this graph is simple:
parallel original edges have distinct inserted vertices.

Self-loops are handled separately:

- each loop is an independent one-edge cycle and an essential member of every
  cycle basis;
- remove loops from the graph passed to subdivision;
- run the simple-graph analysis on the once-subdivided loopless remainder;
- add each loop back as its own relevant cycle, basis member, and one-cycle
  Unique Ring Family.

Every non-loop cycle has a unique cycle in the subdivision graph with twice its
length. Uniform doubling preserves basis weight comparisons, minimum cycle
bases, relevant-cycle membership, and the edge-set relations used to form ring
families. Bounds supplied in original-edge units must be doubled for the
internal run and divided back when results are reported.

The subdivision construction belongs in graph-core:

```rust
impl Graph {
    pub fn subdivide_edges(&self) -> SubdividedGraph;
}
```

`subdivide_edges()` performs exactly one subdivision of every current edge. It
has no repetition count and no once/twice variants. A caller requiring another
subdivision applies the operation to the returned graph. Ring analysis removes
self-loops before calling it, so its fallback requires only one subdivision.

`SubdividedGraph` owns the resulting graph and retains enough correspondence to
map:

- every original node to its subdivision node;
- every original edge to its inserted node;
- every subdivision incidence edge to its original edge.

The result is a graph to which ordinary graph operations can be applied, while
its mappings translate suitable results back to the source graph. It is not a
general substitute for `Graph`: only algorithms with an explicit subdivision
correspondence should consume it.

No new umol-ast subdivision type is needed. The existing
`MoleculeAst::incidence_graph` is already the molecule-level Levi construction:
`topological()` represents atoms and localized bonds, `constitution()` adds
the constitutional overlays, and `full()` also adds stereo entities. That
representation remains useful for symmetry and incidence-based matching.
Ring perception must continue to use localized topology only; overlay
pseudonodes can create Levi-graph cycles that are not chemical rings.

## Cycle, minimum-cycle-basis, and URF results

Cycle enumeration, minimum cycle bases, and Unique Ring Family decomposition
need a graph-core cycle value that preserves both traversal order and edge
identity:

```rust
pub struct Cycle {
    nodes: Vec<NodeId>,
    edges: Vec<EdgeId>,
}
```

Its fields remain private. `nodes()`, `edges()`, `len()`, and `is_empty()`
provide access while construction preserves the invariant that the two ordered
sequences describe the same cycle. `Graph::enumerate_cycles` should return
`Vec<Cycle>` instead of node sequences alone. Edge identity matters because
graph-core permits parallel edges.

The representation admits:

- a self-loop as one node and one edge;
- a parallel-edge cycle as two nodes and two distinct edges;
- an ordinary cycle as equal-length node and edge sequences of length three or
  greater.

Normalization compares the alternating node-and-edge traversal under every
rotation and reversal. It must not use the node sequence alone. This preserves
the identity of parallel edges while making the starting point and traversal
direction irrelevant.

`Cycle` is intentionally close to the private `Ring` in `ast/ring.rs`, but the
types belong to different API layers. `Cycle` uses `NodeId` and `EdgeId` and is
a pure graph result; `Ring` uses `AtomId` and `BondId` and supports molecule
ring views and indices. `Ring` should be constructed from `Cycle` by mapping
both stored sequences directly. It must not reconstruct bonds from adjacent
nodes with `find_edge`, which discards the enumerator's edge identity. A public
generic cycle type or representation-level slice conversion is not warranted
merely to eliminate the similar fields.

`Ring::new` continues to reject mapped cycles with fewer than three distinct
atoms. This is the deliberate graph-cycle to chemical-ring boundary, not an
enumeration failure.

### Minimum cycle basis

Minimum cycle basis is a separate, unbounded graph operation:

```rust
pub enum MinimumCycleBasisAlgorithm {
    Horton,
}

pub struct MinimumCycleBasis {
    cycles: Vec<Cycle>,
    total_length: usize,
}

impl Graph {
    pub fn minimum_cycle_basis(
        &self,
        algorithm: MinimumCycleBasisAlgorithm,
    ) -> MinimumCycleBasis;
}
```

`MinimumCycleBasis` exposes `dimension()`, `total_length()`, and `iter()`.
The dimension is the number of returned cycles and equals `E - V + C`. A
forest returns an empty basis. The operation has no maximum-cycle-size
argument: imposing one can make a cycle-space basis unavailable. It returns
one minimum basis and does not claim that tied bases or their ordering are
canonical.

The formula and operation include multigraph edges. Every self-loop contributes
one independent basis vector. The subdivision fallback computes the remaining
loopless basis and reports `total_length()` in original-edge units rather than
subdivision-edge units. The operation remains infallible.

The result is not a `RingSetKind` or a `RingSet`. It is a linearly independent
basis selected by the requested algorithm. The first implementation is
Horton's minimum-cycle-basis algorithm; shared shortest-path and candidate-cycle
state with relevant-cycle analysis does not change that selector.

### Unique Ring Families

Unique Ring Families likewise have a separate decomposition operation and
result:

```rust
pub struct UniqueRingFamilyId(pub u32);

pub enum UniqueRingFamilyAlgorithm {
    Kolodzik,
}

pub struct UniqueRingFamily {
    nodes: Vec<NodeId>,
    edges: Vec<EdgeId>,
    weight: usize,
    relevant_cycle_count: RelevantCycleCount,
}

pub struct UniqueRingFamilies {
    families: Vec<UniqueRingFamily>,
    node_to_families: Vec<Vec<UniqueRingFamilyId>>,
    edge_to_families: Vec<Vec<UniqueRingFamilyId>>,
    // Private polynomial decomposition state used for relevant-cycle emission.
}

impl Graph {
    pub fn unique_ring_families(
        &self,
        algorithm: UniqueRingFamilyAlgorithm,
    ) -> UniqueRingFamilies;
}
```

`UniqueRingFamilies` provides `count()`, `ids()`, `iter()`, `get()`, and
reverse node/edge membership queries. Each family exposes its node and edge
unions, common cycle weight, and relevant-cycle count. `RelevantCycleCount` is
a dedicated count type rather than `usize`, since the represented relevant
cycle count can be exponential.

The polynomial decomposition must not eagerly materialize every relevant
cycle. The retained decomposition state instead supports explicit lazy
emission through a visitor:

```rust
pub fn visit_relevant_cycles<B>(
    &self,
    family: UniqueRingFamilyId,
    visitor: impl FnMut(Cycle) -> ControlFlow<B>,
) -> ControlFlow<B>;
```

The operation has no maximum-ring-size argument: a bounded result would not be
the graph's URF decomposition. Relevant-cycle prototypes and the finer RCF
decomposition remain internal until a consumer requires them.

The molecule layer maps these results to `AtomId` and `BondId`, following the
existing `GraphView` adapter pattern. Neither result is represented as a
`RingSet`.

For non-simple input, URF construction uses the same loop extraction and
subdivision mapping as relevant-cycle analysis. Each self-loop forms a
one-cycle family. Node and edge unions, common cycle weight, reverse membership,
and relevant-cycle counts are mapped back to the original graph before the
result is exposed. The operation remains infallible.

### Replace the current eager Vismara backend

The current `enumerate_cycles_vismara` implementation is not a suitable base
for minimum-cycle-basis or URF operations. It eagerly expands every shortest
path with `ShortestPathTree::all_paths_to`, forms path-pair products, filters
the resulting cycles, and deduplicates the materialized output. The compact
Vismara cycle-family information needed by the other operations is discarded.

Replace this backend with the two-phase construction used by
RingDecomposerLib:

```text
eager polynomial analysis on original graph
    BCC traversal + non-simple detection
        |
        +-- simple -> shortest-path DAGs
        |
        `-- non-simple -> loop extraction -> subdivision -> BCCs
                             -> shortest-path DAGs

    shortest-path DAGs -> Vismara cycle families -> relevant cycle families

projections from retained relevant-cycle-family state
    -> visit relevant cycles lazily
    -> count relevant cycles
    -> minimum cycle basis
    -> Unique Ring Families
```

This is an algorithm replacement, not a direct port of RingDecomposerLib's C
architecture. graph-core retains its own `Graph`, identifiers, BCC operation,
subdivision correspondence, and reusable traversal facilities. The internal
relevant-cycle-family state is not a new public ring model.

The retained family state is always expressed in original graph identifiers.
Subdivision nodes and edges must not escape through cycle, basis, or family
results.

The eager `Graph::enumerate_cycles` operation may remain as a collecting
wrapper over the streaming traversal. This follows the existing matching API:
the `visit_*` operation owns the enumeration, while `enumerate_*` collects the
same emissions into a `Vec`. The generic visitor and `ControlFlow` add no
required result cloning and allow early termination.

Cycle visitation and collection remain infallible on both direct and
subdivision paths. Structural validation of a molecule is a separate operation;
constructing `RingViews` does not introduce a validation error into ordinary
ring queries.

This work contributes the ring-specific case to the graph-algorithm streaming
audit in doc 105, Other item 30. A true resumable `Iterator` remains distinct
from a visitor and is not required merely to eliminate eager result storage.

## `RingModel` and `RingConfig`

Both types belong in umol-ast next to the ring representation and views. The
model records the requested meaning:

```rust
pub struct RingModel {
    pub kind: RingSetKind,
    pub max_ring_size: usize,
}
```

The config records how each supported ring set is computed. Its public fields
use the family-specific algorithm enums:

```rust
pub struct RingConfig {
    pub simple_cycle_algorithm: SimpleCycleEnumerationAlgorithm,
    pub relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
}
```

The representation is complete without a map, missing entries, private fields,
or compatibility validation. Selecting `RingSetKind::Simple` reads
`simple_cycle_algorithm`; selecting `RingSetKind::Relevant` reads
`relevant_cycle_algorithm`. Each field's type admits only algorithms
implementing that operation.

This structure permits higher-level defaults without assigning a silent
default in graph-core. A fingerprint or aromaticity config may carry a complete
`RingConfig`; the umol-ast operation selects the corresponding public field and
a direct graph-core call continues to receive the family-specific algorithm
explicitly.

## Atom-filter audit

No production caller currently supplies a selective atom filter to
`MoleculeAst::rings_with`.

- `AromaticityPerception` passes `|_| true` for Hückel-rule, HMO, and Clar
  perception.
- The model-specific implementations apply element scope and aromatic
  eligibility after ring enumeration.
- ECFP and Morgan use the unfiltered canonical ring operation.
- The only selective uses are unit and property tests of the ring API itself.

The filter is therefore not a Clar facility. It is currently an unused public
facility for enumerating rings in an induced atom subgraph.

That distinction matters for relevant cycles: relevance in an induced subgraph
is not necessarily the same as enumerating relevant cycles in the full graph
and discarding cycles containing excluded atoms. The closure changes the graph
against which the family is defined; it is not merely an output predicate.

An opaque closure also cannot be inspected, compared, retained in a model, or
represented directly in Python. The ring model/config API should therefore not
incorporate `atom_filter`. Remove it from the public ring operation. The current
call sites do not justify a new scope type, a field on `RingModel`, or a separate
induced-subgraph ring operation. Such an operation can be designed if a real
consumer appears.

## Consequences for doc 151

S9k made the previously hidden Vismara selection explicit and is complete. S9l
and S9m must not freeze the intermediate, family-blind
`cycle_enumeration_algorithm` field into fingerprint and aromaticity configs.
They should consume the settled `RingModel`/`RingConfig` design from this
document. Existing behavior maps every current production ring request to the
Relevant family until a consumer is deliberately migrated to all-simple
semantics.

The existing ring-membership constraints continue to use relevant cycles as a
projection. This is useful for direct SMARTS mapping but is not a satisfactory
general semantics for complex ring conditions. A future constraint redesign
should prefer an unambiguous boolean ring marker and express more complex ring
requirements as substructure constraints rather than flattened ring counts.
That redesign is explicitly outside this round; the existing constraints remain
available.

No design decisions remain before implementation. The staged plan below places
`RingModel` and `RingConfig` in umol-ast; keeps the kind-specific algorithm
selectors and subdivision operation in graph-core; removes `atom_filter`; and
preserves existing constraint behavior.

## Staged implementation plan

New graph-core APIs coexist with the current cycle API until the Rust and Python
callers have migrated. Every subitem includes its unit or property tests and
leaves the workspace green unless it is explicitly marked as a breaking
red-to-green migration. Test additions follow the workspace test-writing
conventions: `rstest` tables use exact structural expectations, property tests
state independent invariants, and benchmark inputs are not used as correctness
references.

### S0 — Baselines and graph values

- **S0a — cycle corpus and performance baseline** **Done**
  (`umol-graph-core/benches/algorithms.rs`,
  `src/algorithms/cycles.rs`): retain the existing molecular-topology corpus
  and record the current Vismara baseline before replacing it. Establish
  reusable benchmark inputs and naming so S1b, S2b, S2c, and S2e can add
  separate all-simple-cycle, basis, relevant-cycle, and URF groups without
  turning this stage into a placeholder for later operations. Add exact tests
  for the small graphs whose cycle sets are independently known; the
  larger molecular graphs remain benchmark inputs rather than expected-output
  cases. **Additive benchmark/test preparation (green).** `[dep: —]`

  The reusable corpus is exposed to the benchmark module by `cycle_corpus()`.
  The current implementation is recorded under the stable Criterion path
  `relevant_cycles/vismara/<case>` and saved locally as baseline
  `doc158-s0a`. The exact table already covers the independently known acyclic,
  single-ring, fused-ring, complete, prism, cube, and naphthalene cases; no
  expected outputs were inferred from the larger benchmark inputs.

  Baseline command:

  ```text
  cargo bench -p umol-graph-core --bench algorithms -- relevant_cycles \
      --save-baseline doc158-s0a
  ```

  Baseline recorded 2026-07-23:

  | Case | Criterion time interval |
  | --- | ---: |
  | `path_6` | 364.44–381.17 ns |
  | `hexagon` | 11.588–12.187 µs |
  | `naphthalene` | 43.781–46.264 µs |
  | `prismane` | 22.273–23.045 µs |
  | `cubane` | 71.967–73.580 µs |
  | `adamantane` | 70.648–73.472 µs |
  | `dodecahedron` | 330.88–346.97 µs |
  | `icosahedron` | 321.57–331.68 µs |
  | `c60` | 3.6622–3.7993 ms |
  | `c70` | 27.838–29.832 ms |
- **S0b — edge-aware `Cycle`** **Done**
  (`umol-graph-core/src/algorithms/cycles.rs`, `src/lib.rs`): add the public
  `Cycle` value with private node/edge sequences, invariant-preserving internal
  construction, accessors, length, equality, hashing, and normalization under
  rotation and reversal. Admit one-edge loops and two-edge parallel cycles and
  distinguish cycles that traverse different `EdgeId`s. Tests use exact whole
  cycles for loop, digon, triangle, rotated, reversed, and parallel-edge cases.
  **Additive API (green).** `[dep: —]`

  `Cycle` is exported from graph-core with private, corresponding node and edge
  sequences and public `nodes()`, `edges()`, and `length()` accessors. Internal
  construction verifies non-empty equal-length sequences, distinct nodes and
  edges, graph membership, and endpoint incidence before normalizing the joint
  representation under rotation and reversal. The existing node-only Vismara
  result remains unchanged but uses the same normalization path.
- **S0c — `SubdividedGraph`** **Done**
  (`umol-graph-core/src/graph.rs`, `src/lib.rs`,
  `benches/algorithms.rs`): add `Graph::subdivide_edges() ->
  SubdividedGraph`. Perform exactly one subdivision, own the resulting graph,
  and expose the source-node, inserted-edge-node, and incidence-edge mappings
  needed to translate algorithm results. Replace the benchmark-local
  `subdivided` helper with the public operation. Unit tests cover empty,
  isolated, path, cycle, parallel-edge, and self-loop inputs; property tests
  prove `|V'| = |V| + |E|`, `|E'| = 2|E|`, endpoint incidence, and mapping
  totality. **Additive API (green).** `[dep: —]`

  `SubdividedGraph` is exported from graph-core together with
  `SubdivisionNodeSource`. Source nodes retain their ids, inserted nodes are
  indexed by source edge id, and the two incidence edges of each source edge
  are consecutive. `node_source` and `edge_source` map subdivision results
  back; `node_of` and `incidence_edges_of` provide the forward mappings. The
  automorphism benchmark uses the public construction, and unit and property
  tests cover empty, isolated, simple, parallel-edge, and self-loop inputs.
- **S0d — exhaustive cycle reference** **Done**
  (`umol-graph-core/tests/property.rs`,
  `tests/property/cycles.rs`, `tests/property/cycles/exhaustive.rs`): split the
  cycle properties into a
  hierarchical module and add deliberately simple, test-only exhaustive
  enumeration that
  enumerates cycles as connected 2-regular edge subsets. Treat self-loops and
  parallel-edge digons as edge-distinct cycles, and provide independent
  cycle-space rank and linear-independence operations. The exhaustive
  implementation must not call any production cycle-enumeration or cycle-space
  implementation. Exact table tests pin it on small named simple and
  non-simple graphs before property tests use it. **Additive test
  infrastructure (green).** `[dep: S0b]`

  The property target is split into `cycles` and `subiso` modules. The cycle
  reference represents a cycle as its sorted edge-id set and exhaustively
  checks every nonempty subset for connected 2-regular incidence, counting a
  loop twice at its endpoint. It computes cycle-space rank with its own
  component traversal and tests linear independence by dynamic GF(2)
  elimination.
  Exact tables cover empty and acyclic graphs, ordinary cycles, loops,
  edge-distinct parallel digons, disconnected cycles, rank, and dependent and
  independent cycle sets.
- **S0e — exhaustive cycle-family reference** **Done**
  (`umol-graph-core/tests/property/cycles/exhaustive.rs`): extend the test-only
  exhaustive implementation to enumerate all cycle bases, retain every minimum
  basis, derive
  relevant cycles as their union, and derive Unique Ring Families directly
  from the defining relation. Keep the implementation intentionally
  exponential and bound its property-test strategies accordingly. Exact table
  tests pin tied minimum bases, relevant-cycle unions, and URF partitions
  before the production implementations are compared with them. **Additive
  test infrastructure (green).** `[dep: S0d]`

  `enumerate_cycle_bases` exhaustively selects every independent
  cycle-space basis, `minimum_cycle_bases` retains every basis of minimum
  total edge count, and `relevant_cycles` returns their sorted union.
  `unique_ring_families` applies the defining pair relation directly:
  relevant cycles must have equal weight, share an edge, and have a symmetric
  difference in the span of strictly shorter cycles; connected components of
  that relation are the families. Exact tables cover unique and tied minimum
  bases, positive and negative URF relations, loops, parallel edges, and
  disconnected graphs. The property strategy exercising the exponential
  family reference is bounded to five edges.
- **S0f — generated exhaustive graph corpus** **Done**
  (`umol-graph-core/tests/data/simple-through-8.g6`,
  `tests/property/corpus.rs`, `tests/property/strategy.rs`): invoke
  `/opt/homebrew/bin/geng` once
  during implementation to generate a checked-in graph6 corpus of
  non-isomorphic simple graphs through the selected order, and record the
  exact generation command and nauty version. Add a test-only graph6 reader and
  an internal bounded generator for graphs with loops and parallel edges.
  Ordinary property and validation tests read the checked-in corpus or use the
  internal generator: they must never spawn `geng`, search for it, skip tests
  when it is absent, or make a `geng` installation a build, test, or development
  prerequisite. Different algorithms may use explicitly documented corpus
  bounds appropriate to their output size. **Additive, self-contained test
  corpus (green).** `[dep: S0d, S0e]`

  `simple-through-8.g6` contains all 13,598 non-isomorphic simple graphs of
  orders one through eight. It was generated once with nauty 2.9.3
  (Homebrew) using:

  ```text
  for order in {1..8}; do
      /opt/homebrew/bin/geng -q "$order"
  done > umol-graph-core/tests/data/simple-through-8.g6
  ```

  The file is 93,655 bytes with SHA-256
  `bbe34489bc2875f5d29b3f4342f6ab6e04d339105d0d89e139f9078f8f0e10f8`
  and XXH3-64 `c3c03691841d9b70`. The test-only `parse_graph6` reader supports
  the corpus's compact graph6 order encoding, `simple_graphs` iterates the
  checked-in data, and `multigraph` provides bounded random graphs with loops
  and repeated edges. Exact tests pin graph6 decoding, counts by order, and the
  corpus digest. No Rust test invokes or searches for `geng`.
- **S0g — algorithm-neutral graph corpus support** **Done**
  (`umol-graph-core/tests/data/simple-through-8.g6`,
  `tests/property/corpus.rs`, `tests/property/strategy.rs`, `tests/property.rs`,
  `tests/property/cycles.rs`): keep the exhaustive simple-graph corpus and its
  graph6 reader outside any algorithm-specific test hierarchy. Keep checked-in
  graph collections in the property target's `corpus` module and generated
  graph inputs in its `strategy` module; cycle properties consume both without
  owning them. Their items use `pub(super)` to remain visible throughout this
  integration-test crate without becoming workspace test support. This makes
  the corpus directly reusable by other property modules while preserving its
  S0f contents, provenance, integrity checks, and no-`geng` runtime policy.
  **Property-suite reorganization without behavioral changes (green).**
  `[dep: S0f]`

### S1 — Complete bounded cycle enumeration

- **S1a — local shortest-cycle multigraph semantics** **Done**
  (`umol-graph-core/src/algorithms/cycles.rs`): make the existing shortest-cycle
  queries return length one for a self-loop and length two for parallel edges,
  while preserving existing simple-graph results. Table tests assert exact
  edge- and node-local answers for loops, digons, ordinary cycles, bridges, and
  acyclic nodes. **Behavioral correction without a signature change (green).**
  `[dep: S0b]`

  Edge-local BFS handles a self-loop directly before excluding the queried
  edge; for distinct endpoints, the existing edge-exclusion search finds a
  parallel alternative as a one-edge path and therefore returns a two-edge
  cycle. Node-local queries retain their minimum-over-incident-edges semantics.
  Exact tables cover both public queries on loops, digons, ordinary cycles,
  bridges, and acyclic nodes.
- **S1b — edge-aware Read--Tarjan visitor and collector** **Done**
  (`umol-graph-core/src/algorithms/cycles.rs`,
  `src/algorithms/cycles/simple.rs`, `src/lib.rs`,
  `benches/algorithms.rs`): add `SimpleCycleEnumerationAlgorithm::ReadTarjan`,
  `visit_simple_cycles`, and the collecting `enumerate_simple_cycles` wrapper.
  Apply the maximum length during search, use edge-aware paths, emit loops and
  digons directly, normalize once at emission, and propagate visitor
  `ControlFlow` without materializing unused cycles. Exact tests cover bounds,
  disconnected components, fused and bridged graphs, parallel alternatives,
  deterministic collection order, and early termination. Property tests
  compare small random multigraphs against the S0d exhaustive reference and
  verify visitor/collector equality, uniqueness, and relabeling invariance.
  Extend the
  baseline corpus with bounded and unbounded
  all-simple-cycle benchmarks. **Additive API (green).**
  `[dep: S0a, S0b, S0d]`

  `SimpleCycleEnumerationAlgorithm::ReadTarjan` is exported with
  `visit_simple_cycles` and its collecting `enumerate_simple_cycles` wrapper.
  The edge-aware search emits loops directly, retains only extensions with a
  bounded return path that avoids the current path, and advances through a
  single fruitful extension without recursion. The minimum node fixes cycle
  rotation; ordering the first and closing edges fixes reversal while retaining
  every edge-distinct parallel cycle. `Cycle` is constructed and normalized
  only when emitted, and `ControlFlow::Break` propagates immediately.

  Exact tables cover bounds, disconnected and fused components, bridges,
  parallel alternatives, deterministic traversal, and early termination.
  Properties compare bounded results with independent exhaustive edge-subset
  enumeration, prove visitor/collector equality and uniqueness, and preserve
  edge-cycle sets under node relabeling. Criterion groups
  `simple_cycles_bounded_8/read_tarjan/*` and
  `simple_cycles_unbounded/read_tarjan/*` cover the full bounded molecular
  corpus and the repeatable unbounded small/medium corpus, respectively.
- **S1c — simple-cycle differential validation** **Done**
  (`umol-graph-core/tests/data/simple-cycles/`): run the
  completed Read--Tarjan implementation over the S0f corpus and compare it
  first with the S0d exhaustive reference, then independently with NetworkX
  and igraph from the `work2` environment. Restrict NetworkX comparison to
  simple-graph node-cycle semantics; use igraph edge-cycle results for
  non-simple cases only
  after pinning its loop and parallel-edge behavior with exact probes. Store
  normalized failing cases and a validation report, but do not make either
  Python library or the micromamba environment a prerequisite of the Rust
  property suite. **One-off differential validation with checked-in evidence
  (green).** `[dep: S0f, S1b]`

  The exhaustive Rust reference agrees with Read--Tarjan on all 208
  non-isomorphic simple graphs through order six. NetworkX 3.5, igraph 1.0.1,
  and Read--Tarjan agree on all 13,598 simple graphs through order eight,
  containing 1,526,236 cycles. After exact probes fixed igraph's edge-cycle
  treatment of loops and parallel edges, igraph and Read--Tarjan also agree on
  all 3,453 non-simple canonical edge multisets through four nodes and five
  edges. No failing cases were produced.

  The normalized external answers are checked in as
  `simple-through-8.tsv` (13,598 rows, 17,818,141 bytes, SHA-256
  `7ae435ea65b1021a183a30508ec692353b66040a643d1338b6bcc5a3787e0070`)
  and `multigraph-through-4-edges-5.tsv` (3,453 rows, 98,828 bytes, SHA-256
  `287b0297108954bd647840fe87a1c2e56ec60280fd2b950d35a720c5d389b2ef`).
  Rust-only regression tests compare graph-core with these captured results;
  they never import, invoke, search for, or conditionally skip based on
  NetworkX, igraph, Python, or `work2`.

  The one-off generator
  `tests/data/simple-cycles/generate_captured_results.py` uses NetworkX
  3.5 and python-igraph 1.0.0 backed by igraph 1.0.1. It normalizes
  simple-graph node cycles over rotation and reversal, and non-simple edge
  cycles as sorted edge-id sets; `-` denotes an empty cycle set. Before writing
  the captured results, it verifies igraph's treatment of one and multiple
  loops, two and three parallel edges, a parallel edge plus a triangle, and
  combined loops, digons, and triangles. It then requires NetworkX and igraph
  to agree on every simple-graph node-cycle set before capturing the external
  results. Regeneration is deliberately separate from the Rust test suite:

  ```text
  micromamba run -n work2 python \
      umol-graph-core/tests/data/simple-cycles/generate_captured_results.py \
      umol-graph-core/tests/data/simple-through-8.g6 \
      umol-graph-core/tests/data/simple-cycles/simple-through-8.tsv \
      umol-graph-core/tests/data/simple-cycles/multigraph-through-4-edges-5.tsv
  ```

### S2 — Cycle-space, relevant-cycle, and URF analysis

- **S2a — shared cycle-space kernel** **Done**
  (`umol-graph-core/src/algorithms/cycles.rs`,
  `src/algorithms/cycles/basis.rs`,
  `src/algorithms/cycles/relevant.rs`): add the internal edge-vector,
  independence, shortest-path-DAG, and candidate-cycle machinery shared by
  Horton and the Vismara/RCF analysis. Keep subsidiary BCC selection explicit
  inside the named algorithms. Unit tests pin rank, independence, symmetric
  difference, shortest-path alternatives, and candidate reconstruction;
  property tests compare the computed cycle-space rank with `E - V + C`.
  **Additive internals (green).** `[dep: S0b]`

  `EdgeVector` stores edge incidence over GF(2), and `CycleVectorBasis`
  maintains a high-pivot row-echelon basis with insertion reporting linear
  independence. `cycle_space_rank` computes `E - V + C` using the explicit
  `ConnectedComponentsAlgorithm::Bfs` selection. `ShortestPathDag` retains
  every node/edge predecessor alternative, materializes exact `ShortestPath`
  values, and `ShortestPath::cycle_with` reconstructs a `Cycle` from two
  compatible paths and a closing edge. The existing Vismara implementation
  now consumes this path machinery and checks its full, unbounded result
  against the cycle-space rank in debug builds; its public result type and
  algorithm selection are unchanged. Exact tests cover symmetric difference,
  independent and dependent insertion, rank, path alternatives, parallel
  edges, candidate reconstruction, and intersecting-path rejection. A bounded
  multigraph property verifies that the span of all enumerated simple cycles
  has rank `E - V + C`.
- **S2b — Horton minimum cycle basis** **Done**
  (`umol-graph-core/src/algorithms/cycles.rs`,
  `src/algorithms/cycles/basis.rs`, `src/lib.rs`,
  `benches/algorithms.rs`): add `MinimumCycleBasisAlgorithm::Horton`,
  `MinimumCycleBasis`, and `Graph::minimum_cycle_basis`. Extract self-loops as
  independent basis members, use `SubdividedGraph` for the remaining
  non-simple topology, and report cycles and total length in source identifiers
  and source-edge units. Exact tests cover forests, disconnected graphs,
  tied bases, loops, and parallel edges. Property tests verify dimension,
  independence, spanning of the cycle space, and minimal total length against
  the S0e exhaustive bases on small graphs. Add basis benchmarks over the S0
  corpus. **Additive API (green).** `[dep: S0a, S0b, S0c, S0e, S2a]`

  `MinimumCycleBasisAlgorithm::Horton`, `MinimumCycleBasis`, and
  `Graph::minimum_cycle_basis` are exported from graph-core. The result exposes
  `dimension`, `total_length`, and `iter`; it uses source `NodeId` and `EdgeId`
  values and makes no canonical-ordering claim. Horton constructs one
  deterministic shortest-path tree per root, orders candidate cycles by source
  length and identifiers, and greedily selects independent edge vectors until
  the cycle-space rank is reached.

  Self-loops are extracted as mandatory one-edge basis members. If the
  remaining loopless graph has parallel edges, it is subdivided once, solved
  as a simple graph, and mapped back to source identifiers. Reported lengths
  are divided structurally by this mapping rather than retaining subdivision
  edge counts. Exact tests cover forests, disconnected cycles, tied K4 bases,
  loops, digons, and three parallel edges. The bounded multigraph property
  compares dimension, independence, spanning, and total length with exhaustive
  minimum bases; an extended 4,096-case run also passes.

  Criterion group `minimum_cycle_basis/horton/*` covers the S0 corpus. A
  10-sample verification run on 2026-07-23 measured:

  | Case | Criterion time interval |
  | --- | ---: |
  | `path_6` | 574.19–575.96 ns |
  | `hexagon` | 11.251–11.315 µs |
  | `naphthalene` | 36.543–36.608 µs |
  | `prismane` | 21.672–21.701 µs |
  | `cubane` | 38.160–38.225 µs |
  | `adamantane` | 62.441–62.552 µs |
  | `dodecahedron` | 286.62–289.17 µs |
  | `icosahedron` | 188.31–188.60 µs |
  | `c60` | 3.2030–3.2084 ms |
  | `c70` | 4.5588–4.5669 ms |
- **S2c — compact Vismara/RCF analysis on simple graphs** **Done**
  (`umol-graph-core/src/algorithms/cycles/relevant.rs`): replace eager
  shortest-path expansion with retained shortest-path DAGs, Vismara cycle
  families, and relevant-cycle-family state on the direct simple-graph path.
  Keep the new analysis internal until the total public operation is available.
  Exact tests cover odd and even prototypes, multiple shortest paths, fused and
  bridged systems, family counts, and lazy early termination. Property tests
  compare small simple graphs with the S0e definition-level reference.
  Benchmark the new simple path beside the S0 baseline before the legacy
  backend is removed. **Additive replacement internals (green).**
  `[dep: S0a, S0b, S0e, S2a, S2b]`

  `RelevantCycleAnalysis` retains one admissible shortest-path DAG per
  Vismara root and only those cycle families whose prototypes are independent
  of every strictly shorter prototype. Odd and even prototypes retain source
  node and edge identifiers. Alternative shortest paths are traversed directly
  from the DAGs through `ControlFlow`; neither the path alternatives nor their
  Cartesian products are materialized. The existing node-only
  `Graph::enumerate_cycles` operation now collects this traversal while its
  public compatibility surface remains unchanged.

  Exact tests cover odd and even prototypes, an unequal theta graph with
  multiple shortest paths, fused and bridge-connected ring systems, exact
  family prototypes and emitted cycles, and early termination at both the path
  and cycle-family levels. A bounded simple-graph property compares the
  emitted edge sets with the S0e exhaustive union of all minimum cycle bases;
  the extended 4,096-case run passes.

  The existing `relevant_cycles/vismara/*` Criterion group was rerun with ten
  samples against the intervals recorded by S0a:

  | Case | S0a eager interval | S2c compact interval |
  | --- | ---: | ---: |
  | `path_6` | 364.44–381.17 ns | 405.10–406.13 ns |
  | `hexagon` | 11.588–12.187 µs | 5.6489–5.6815 µs |
  | `naphthalene` | 43.781–46.264 µs | 10.656–10.681 µs |
  | `prismane` | 22.273–23.045 µs | 11.408–11.422 µs |
  | `cubane` | 71.967–73.580 µs | 19.908–20.302 µs |
  | `adamantane` | 70.648–73.472 µs | 25.791–25.869 µs |
  | `dodecahedron` | 330.88–346.97 µs | 89.248–92.303 µs |
  | `icosahedron` | 321.57–331.68 µs | 125.13–125.33 µs |
  | `c60` | 3.6622–3.7993 ms | 581.13–582.10 µs |
  | `c70` | 27.838–29.832 ms | 1.0181–1.0672 ms |

  The compact construction is faster on every cyclic corpus case, from about
  twofold on the small prism and hexagon cases to about 27-fold on C70. The
  acyclic path remains sub-microsecond but is approximately 7–11% slower.
- **S2d — total relevant-cycle public API** **Done**
  (`umol-graph-core/src/algorithms/cycles.rs`,
  `src/algorithms/cycles/relevant.rs`, `src/lib.rs`): add
  `RelevantCycleEnumerationAlgorithm::Vismara`,
  `visit_relevant_cycles`, and `enumerate_relevant_cycles`. Fuse non-simple
  detection into initial analysis; continue directly for simple graphs, or
  extract self-loops, subdivide the loopless remainder once, and map the result
  back. Preserve source-edge bounds, cycle identities, collection order, and
  infallibility. Tests cover direct/fallback parity, standalone loops, parallel
  digons, longer parallel alternatives, mixed disconnected graphs, and visitor
  termination. Property tests compare small multigraph results with the S0e
  exhaustive minimum-basis union. **Additive API (green).**
  `[dep: S0b, S0c, S0e, S2c]`

  `RelevantCycleEnumerationAlgorithm::Vismara` now selects the public
  `Graph::visit_relevant_cycles` and `Graph::enumerate_relevant_cycles`
  operations. The size bound counts source edges. Both operations return
  normalized `Cycle` values carrying source `NodeId` and `EdgeId` identities.
  The visitor can terminate traversal through `ControlFlow`; the collector is
  implemented by collecting that traversal.

  The total traversal extracts source self-loops during its first edge pass.
  It runs the compact Vismara analysis directly when the remaining graph is
  simple. Parallel edges instead trigger one subdivision of the loopless
  remainder, followed by source-identifier projection. Graphs containing loops
  but no parallel edges use a compact loopless graph and its source-edge map.
  The legacy node-only Vismara operation delegates to the new collector while
  retaining its former exclusion of loops and digons.

  Exact tests cover the direct and subdivision paths, source-distinct loops,
  parallel digons and longer alternatives, mixed disconnected graphs, source
  size bounds, deterministic collection, and visitor termination. Bounded
  multigraph properties compare the public result with the S0e exhaustive
  definition and prove visitor/collector agreement. The exhaustive comparison
  passes 4,096 generated cases; the full `umol-graph-core` test suite and
  all-targets Clippy gate pass.
- **S2e — Unique Ring Family decomposition** **Done**
  (`umol-graph-core/src/algorithms/cycles.rs`,
  `src/algorithms/cycles/urf.rs`, `src/lib.rs`,
  `benches/algorithms.rs`): add `UniqueRingFamilyAlgorithm::Kolodzik`,
  `RelevantCycleCount`, `UniqueRingFamilyId`, `UniqueRingFamily`, and
  `UniqueRingFamilies`, including reverse node/edge membership and lazy
  per-family relevant-cycle visitation. Retain only polynomial decomposition
  state, use arbitrary-precision counts, add one family per extracted
  self-loop, and expose only source identifiers. Exact tests cover independent,
  fused, bridged, symmetric, loop, and parallel-edge families. Property tests
  compare partitions with S0e and prove that stored counts equal lazy emission
  on tractable graphs, reverse indices agree with family unions, and relabeling
  preserves the decomposition. Add decomposition and lazy emission benchmarks.
  **Additive API (green).**
  `[dep: S0a, S0b, S0c, S0e, S2b, S2c, S2d]`

  `Graph::unique_ring_families` now returns the polynomial
  `UniqueRingFamilies` decomposition selected by
  `UniqueRingFamilyAlgorithm::Kolodzik`. Each `UniqueRingFamily` exposes its
  source-node and source-edge unions, common source-edge weight, and exact
  `RelevantCycleCount`. The collection provides family identifiers,
  descriptor iteration and lookup, reverse node/edge membership, and
  `visit_relevant_cycles` for explicit lazy expansion with `ControlFlow`
  termination.

  The implementation groups retained relevant-cycle families by the defining
  equal-weight relation: their prototypes have equal remainders modulo the
  strictly smaller cycle space and their edge unions intersect. Shortest-path
  DAGs provide polynomial node/edge unions and arbitrary-precision path-product
  counts without materializing the cycles. Self-loops become one-cycle
  families. Parallel-edge inputs use one subdivision and all public
  descriptors and emitted cycles are projected back to source identifiers.

  Exact tests cover independent, fused, bridge-connected, symmetric,
  self-loop, and parallel-edge systems, including a family representing
  multiple lazily emitted cycles and visitor termination. Bounded multigraph
  properties compare full partitions and counts with the S0e exhaustive
  definition, verify reverse membership, and preserve the full decomposition
  under node relabeling. The extended property set passes 4,096 generated
  cases; the full `umol-graph-core` suite and all-targets Clippy gate pass.

  Criterion groups `unique_ring_families/decomposition/*` and
  `unique_ring_families/lazy_emission/*` cover the S0 corpus. A short
  ten-sample verification run on 2026-07-23 measured:

  | Case | Decomposition | Full lazy emission |
  | --- | ---: | ---: |
  | `path_6` | 825.61–832.11 ns | empty |
  | `hexagon` | 6.7556–6.8082 µs | 0.77046–1.6456 µs |
  | `naphthalene` | 13.380–13.460 µs | 1.5600–1.5734 µs |
  | `prismane` | 14.327–14.546 µs | 2.9883–3.0111 µs |
  | `cubane` | 24.162–24.393 µs | 3.6549–3.6848 µs |
  | `adamantane` | 31.974–32.237 µs | 4.6203–4.6645 µs |
  | `dodecahedron` | 106.76–107.90 µs | 8.6328–8.7013 µs |
  | `icosahedron` | 144.29–144.93 µs | 11.853–14.791 µs |
  | `c60` | 690.98–695.02 µs | 24.284–24.492 µs |
  | `c70` | 1.1666–1.1812 ms | 28.098–28.412 µs |
- **S2f — cycle-family differential validation**
  (`umol-graph-core/tests/data/cycles/`): run the
  completed MCB, relevant-cycle, and URF operations over the bounded S0f corpus.
  Compare MCB dimension and total weight with NetworkX, igraph, CDK, and
  RingDecomposerLib without requiring the same member cycles when minimum bases
  are tied; compare normalized relevant-cycle edge sets with CDK and
  RingDecomposerLib; and compare URF partitions, counts, and lazy emission with
  RingDecomposerLib. Run RingDecomposerLib's independent exponential validator
  on the shared tractable cases. Record source revisions, commands, semantic
  exclusions, and normalized failures. These external tools corroborate the
  S0e exhaustive reference and remain one-off validation dependencies, not
  Rust property-test dependencies. **One-off differential validation with
  checked-in evidence
  (green).** `[dep: S0f, S2b, S2d, S2e]`

### S3 — umol-ast ring model and views

- **S3a — `RingSetKind` rename**
  (`umol-ast/src/ast/{ring,molecule,view/ring}.rs`, `src/ast.rs`, all workspace
  Rust callers): rename `RingFamily` to `RingSetKind` and migrate imports,
  `family` fields, accessors, and arguments to `kind`, and fixtures and
  assertions without a compatibility alias.
  Existing unit and property tests retain their exact semantics under the new
  name. **Done.** The enum, API vocabulary, and every workspace Rust caller now
  use `RingSetKind`/`kind`; no compatibility alias remains. **Breaking rename
  and complete caller migration (red→green).**
  `[dep: —]`
- **S3b — `RingModel` and `RingConfig`**
  (`umol-ast/src/ast/ring.rs`, `src/ast.rs`): add the public AST-layer model
  and config. `RingModel` holds `kind` and `max_ring_size`; `RingConfig` exposes
  separate public simple- and relevant-cycle algorithm fields. Add ordinary
  high-level defaults for relevant rings up to size 22 and the settled
  Read--Tarjan/Vismara implementations. Table tests assert exact construction,
  defaults, and equality; no compatibility-validation error is introduced.
  **Done.** `RingModel` and
  `RingConfig` are public AST-layer value types with explicit public fields;
  their defaults select Relevant/22 and Read--Tarjan/Vismara respectively.
  Dispatch remains in S3c because the two algorithm fields have different
  selector types. **Additive API (green).**
  `[dep: S1b, S2d, S3a]`
- **S3c — infallible ring construction and public entry point**
  (`umol-ast/src/ast/{ring,molecule,view/graph,view/ring}.rs`, all workspace
  Rust callers): make `RingSet` consume graph-core `Cycle` values and map both
  nodes and edges directly. Remove `find_edge` reconstruction, the redundant
  induced-cycle filter, `atom_filter`, `rings_with`, and the family-blind cycle
  selector. Make the general molecule ring entry point accept `RingModel` and
  `RingConfig`, dispatch to the kind-specific graph-core operation, and keep
  `RingViews` infallible. Split `GraphView` cycle enumeration into typed simple
  and relevant methods. Migrate every Rust caller in the same subitem using the
  former effective Relevant/Vismara behavior. Exact tests cover both kinds,
  selection of the corresponding `RingConfig` field, bounds, bond identity,
  absence of one- and two-atom chemical rings, and
  structurally invalid but non-panicking views; property tests preserve ring
  reindexing and view/index consistency. **Done.** `RingSet::enumerate` now
  takes the graph first, followed by `RingModel` and `RingConfig`, dispatches
  to the typed graph-core collectors, maps both node and edge identifiers
  directly, and rejects graph cycles shorter than three atoms only at the
  chemical-ring boundary. `MoleculeAst::rings` is the sole general entry
  point; `GraphView` exposes separate simple- and relevant-cycle collectors;
  `rings_with`, `atom_filter`, induced-cycle filtering, and AST-layer use of
  the legacy selector are removed. Every Rust caller uses the former effective
  Relevant/Vismara behavior pending its operation-specific S4 configuration.
  **Breaking ring API and complete Rust caller migration (red→green).**
  `[dep: S1b, S2d, S3b]`
- **S3d — retained ring-constraint behavior**
  (`umol-ast/src/ast/view/ring.rs`,
  `src/ast/constraint/{atom,bond,dative,ring}.rs`, and ring property tests):
  make the ring views used to ground or match existing ring-membership
  constraints use the fixed Relevant projection, without adding a ring-set
  parameter to `RingMembershipAst` or `RingScope`. Tests retain exact
  boolean/count/size behavior and make the fixed Relevant choice explicit.
  The boolean-marker/substructure redesign remains outside this round.
  **Done.** Substructure host-target derivation now materializes only the atom
  and localized-bond ring constraints requested anywhere in the pattern,
  using an explicit Relevant projection through size 22 and the default
  relevant-cycle implementation. K4 matching cases distinguish these counts
  from all-simple-cycle counts for both atoms and bonds. Ring constraint and
  view documentation records the fixed projection, and generated ring
  properties cover exact membership booleans plus total and per-size counts.
  Dative-bond ring membership remains an asserted constraint: deriving it
  requires the separately deferred definition of ring topology containing
  dative overlays rather than the localized graph used here.
  **Internal migration with unchanged public API (green).** `[dep: S3c]`

### S4 — Rust workflow configuration

- **S4a — fingerprint ring configuration**
  (`umol-graph/src/fingerprint/{ecfp,morgan,featurizer,reaction}.rs`): add the
  AST `RingConfig` to ECFP and Morgan featurizers and thread it through enum and
  reaction-featurizer composition. Their ordinary constructors retain the
  inspectable Relevant/Vismara workflow default; struct literals and
  conversions are migrated together. Exact fingerprint fixtures prove default
  identity stability, and configured tests prove the relevant selector reaches
  ring construction. **Breaking config-shape migration (red→green).**
  `[dep: S3c]`
- **S4b — aromaticity ring configuration**
  (`umol-graph/src/ops/aromaticity.rs`,
  `ops/aromaticity/{hueckel_rule,hmo,clar}.rs`,
  `ops/{resolve,transform,validate}` aromaticity callers): introduce the
  operation-level `AromaticityConfig` already specified by doc 151, embedding
  AST `RingConfig` beside connected-components and
  maximum-independent-set selectors. Perception constructs the fixed Relevant
  `RingModel` with each aromaticity model's size bound and passes the selected
  algorithm explicitly. Migrate resolver, aromatizer, validator, fixtures, and
  conformance tests together. Exact tests preserve aromatic systems and
  transactional results under defaults and exercise explicit selector
  propagation through Hückel-rule, HMO, and Clar paths. **Additive config
  followed by breaking caller migration (red→green).** `[dep: S3c]`

### S5 — Python configuration surface

- **S5a — family-specific algorithm values and `RingConfig`**
  (`umol-py/src/{algorithm,ring}.rs`, `src/lib.rs`,
  `python/umol/__init__.py`): replace the Python
  `CycleEnumerationAlgorithm` with separate simple- and relevant-cycle
  algorithm classes and bind keyword-only `RingConfig` fields with the settled
  defaults. Implement inherent `from_rust`/`to_rust`, equality, and repr.
  Variant-complete Rust table tests and installed-package tests cover exports,
  construction, defaults, and repr. **Breaking Python selector replacement and
  export migration (red→green).** `[dep: S1b, S2d, S3b]`
- **S5b — fingerprint config integration**
  (`umol-py/src/fingerprint/config.rs`, fingerprint and reaction bindings,
  `umol-py/tests/test_fingerprint.py`): add keyword-only nested `ring_config`
  to ECFP and Morgan configuration variants, lower it into the Rust
  featurizers, and carry it through reaction fingerprint configs. Installed
  tests compare default and explicit configuration and retain exact fingerprint
  payloads. **Breaking Python config-shape migration (red→green).**
  `[dep: S4a, S5a]`
- **S5c — aromaticity config integration**
  (`umol-py/src/model/aromaticity.rs`, resolver/aromatizer/validator bindings,
  installed aromaticity and resolution tests): bind the Rust
  `AromaticityConfig` with keyword-only ring, connected-components, and
  independent-set selectors; nest it in the existing operation configs and
  migrate every Python caller. Tests cover defaults, explicit nested configs,
  equality, repr, lowering, and unchanged resolved/aromatized structures.
  **Breaking Python config-shape migration (red→green).** `[dep: S4b, S5a]`

### S6 — Legacy removal and release gates

- **S6a — retire the intermediate cycle API**
  (`umol-graph-core/src/algorithms/cycles.rs`, `src/lib.rs`,
  `benches/algorithms.rs`, remaining workspace callers): remove
  `CycleEnumerationAlgorithm`, `Graph::enumerate_cycles`, the eager
  `ShortestPathTree::all_paths_to` Vismara backend, and obsolete benchmark/test
  adapters. The family-specific visitor/collector APIs become the only public
  enumeration surface. Run the exact fixture suites, graph-core and AST
  property suites, Rust workflow tests, installed Python tests, workspace
  clippy, the criterion comparison against S0, and the one-off S1c/S2f
  differential validations before accepting the removal. The ordinary gate
  remains reproducible without `geng` because it consumes the checked-in S0f
  corpus. **Breaking legacy removal with all callers already migrated
  (red→green).**
  `[dep: S1c, S2f, S3c, S4a, S4b, S5a, S5b, S5c]`
- **S6b — documentation and status closure**
  (`discussion/151-python-molecule-workflows-2026-07-13.md`, this document,
  `discussion/000-status.md`): update the superseded S9k/S9l/S9m descriptions
  to the family-specific AST design; record the measured benchmark comparison,
  exhaustive corpus bounds, generation provenance, external implementation
  revisions, validation commands, and deliberate semantic exclusions; and mark
  doc 158 complete only after every S6a gate passes. **Documentation closure
  after a green workspace (green).** `[dep: S6a]`

### Critical path and deferral

The critical path is:

```text
S0b -> S0d -> S0e -> S0f
          \-> S1b ----------------> S1c <- S0f
S0b -> S2a --\
S0e ----------> S2b -> S2c -> S2d -> S2e -> S2f <- S0f
                                 |
{S1b, S2d, S3a} -> S3b -> S3c -> {S4a, S4b}
                                          |
                                          -> {S5b, S5c}
                                                 |
{S1c, S2f, S5b, S5c} -------------------------> S6a -> S6b
```

`S0c` joins the cycle-space path at S2b and S2d; `S3a` can proceed independently
without blocking graph-core development. S1c and S2f are one-off implementation
gates rather than recurring external prerequisites; their checked-in evidence
and the S0f corpus keep subsequent Rust verification self-contained. No
implementation stage is deferrable in this round: Read--Tarjan, Horton,
Vismara/RCF, URF, AST integration, workflow configuration, Python migration,
and independent validation are all part of the settled deliverable.
Essential cycles, exact RDKit SSSR compatibility, general path/tree
enumeration, the ring-constraint redesign, and CX/MOL/SDF work remain outside
this plan rather than becoming deferred stages.
