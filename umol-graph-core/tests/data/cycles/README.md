# Captured cycle-family results

`cycle-families-through-8.tsv` contains independently calculated minimum-cycle-
basis, relevant-cycle, and Unique Ring Family results for every non-isomorphic
simple graph through order eight in `../simple-through-8.g6`.

The Rust property suite reads the captured TSV directly. Python, Java, CDK,
RingDecomposerLib, NetworkX, and igraph are regeneration tools, not test
dependencies.

## TSV schema

Each tab-separated row contains:

1. the source graph in graph6 form;
2. the minimum cycle basis dimension;
3. the minimum cycle basis total length;
4. normalized relevant cycles;
5. normalized Unique Ring Families;
6. RingDecomposerLib exponential-validation status.

A cycle is encoded as its sorted source edge identifiers. Cycles are separated
by `;` and ordered first by length and then lexicographically. `-` represents an
empty cycle set.

Each Unique Ring Family is encoded as:

```text
weight:relevant-cycle-count:nodes:edges:cycles
```

Families are separated by `/`. Within one family, nodes and edges are sorted
identifiers, cycles are separated by `.`, and the family ordering is by weight,
cycles, nodes, and edges.

The validation-status field is:

- `1` when RingDecomposerLib's exponential family validator ran and accepted
  the result;
- `0` when relevant cycles exist but exponential validation was skipped because
  the graph exceeds the configured order limit;
- `-` when the graph has no relevant cycles.

Minimum-cycle-basis member cycles are intentionally absent. A minimum cycle
basis need not select the same members when minimum bases are tied; only its
dimension and total length are compared.

## Independent calculations

`generate_captured_results.py` performs the following comparisons:

- NetworkX and igraph independently calculate minimum-cycle-basis dimension
  and total length.
- CDK calculates the same two basis values and the normalized relevant
  edge-cycle set.
- RingDecomposerLib calculates the basis values, relevant edge-cycle set, and
  complete Unique Ring Family decomposition.
- RingDecomposerLib's eager and iterator relevant-cycle operations must agree.
- Family counts, weights, node unions, edge unions, and emitted cycles must
  agree internally.
- CDK and RingDecomposerLib must return identical normalized relevant-cycle
  sets.
- RingDecomposerLib's definition-level exponential validator is applied to
  cyclic graphs through the configured order limit.

The resulting TSV is written only after these comparisons succeed for every
input graph.

## Last captured run

The checked-in results were generated with:

```text
NetworkX                 3.5
python-igraph            1.0.0
igraph                   1.0.1
CDK revision             7e130959efbc0d3561e3437175ca5da83147c298
RingDecomposerLib        3a7ff93de0d9c4f6a5661508549c6063573f39c7
validation through order 5
```

The run covered:

```text
graphs                                  13,598
normalized relevant cycles             131,589
normalized Unique Ring Families        116,205
exponentially validated cyclic graphs       30
normalized comparison failures               0
```

The captured scope is finite, simple, unweighted, undirected graphs. Multigraph
semantics are verified separately by the Rust property suite.

The input and output hashes are:

```text
simple-through-8.g6
bbe34489bc2875f5d29b3f4342f6ab6e04d339105d0d89e139f9078f8f0e10f8

cycle-families-through-8.tsv
12d68a42fd387ac94a37002ccb7aa82ebbf7b40befd9ce66015e0434169c2769
```

## Regeneration

The generator requires prebuilt CDK and RingDecomposerLib source trees plus a
Python environment containing NetworkX and python-igraph.

```sh
cmake -S "$RDL_SOURCE" -B "$RDL_BUILD" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build "$RDL_BUILD"

mvn -f "$CDK_SOURCE/pom.xml" -pl base/core -am package -DskipTests

micromamba run -n work2 python \
    umol-graph-core/tests/data/cycles/generate_captured_results.py \
    umol-graph-core/tests/data/simple-through-8.g6 \
    umol-graph-core/tests/data/cycles/cycle-families-through-8.tsv \
    --cdk-source "$CDK_SOURCE" \
    --rdl-source "$RDL_SOURCE" \
    --rdl-library "$RDL_BUILD/src/RingDecomposerLib/libRingDecomposerLib.dylib"
```

The generator prints tool versions, revisions, counts, and hashes after a
successful run. Update this README when replacing the captured results.
