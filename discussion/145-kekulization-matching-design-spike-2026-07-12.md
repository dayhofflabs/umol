# Kekulization and matching design spike

## Scope

This spike separates three related projects:

1. make single-structure kekulization work for heterocycles and charged non-benzenoids;
2. evaluate output-sensitive enumeration of all Kekulé structures, especially Uno's algorithms;
3. remove the localized impedance mismatch between matching results and graph correspondences when
   algorithmic work proceeds.

The first is the immediate product improvement. The second is a research/benchmark track. The third
should land with the first matching change rather than remain another manual id-map.

## Existing implementation

`umol-graph/src/ops/transform/kekulizer.rs` extracts each aromatic system, builds a manual host→sub
node map and sub→host bond vector, and asks `Graph::perfect_matching` for a perfect matching. The
only single-matching implementation is deterministic recursive DFS. It therefore assumes every
aromatic-system atom participates in exactly one localized double bond.

`umol-graph-core/src/algorithms/matching.rs` also contains:

- a general-graph Edmonds maximum matching implementation;
- a bipartite augmenting-path implementation named `HopcroftKarp` (currently repeated BFS,
  `O(V(V+E))`, not the full layered `O(E sqrt(V))` algorithm);
- perfect- and maximum-matching enumeration by binary include/exclude branch-and-bound;
- eager `Vec<Matching>` enumeration results.

The enumerator is currently not a safe reference implementation. Its `can_reach` function is
described as an Edmonds oracle but computes a greedy residual matching and rejects a branch when the
greedy matching is too small. A greedy matching is a lower bound on the residual maximum, not an
upper bound, so this can prune a branch that contains valid target-size matchings. This must be
replaced by an exact extension oracle or removed before enumeration results are trusted.

## Chemical problem statement

Let the aromatic-system graph be `G=(V,E)`. A localized double-bond assignment is a matching `M`.
An atom incident to `M` has one localized π bond; an exposed atom has no localized π bond and must
carry the appropriate local 0e/2e state.

The aromatic system already stores positional per-atom electron contributions. This provides more
information than the current kekulizer uses:

- contribution `1`: atom must be covered by a matching edge;
- contribution `2`: prescribed exposed 2e donor (pyrrole/furan/thiophene class);
- contribution `0`: prescribed exposed 0e acceptor (borepin/carbocation class).

For a heterogeneous system these locations are chemically meaningful and charge is normally already
localized on the heteroatom. If the exposed set `H` is prescribed, the feasibility problem is simply
a perfect matching of the induced graph `G-H`. Calling this a maximal-matching problem is misleading:
maximality would additionally require the exposed vertices to form an independent set, which is not
part of the electron-demand constraint.

Homoelement charged systems are different. Charge equalization deliberately removes an arbitrary
atom-local charge and stores it on the aromatic system. Cp− and tropylium therefore contain only
per-atom `1` contributions but require one unspecified exposed atom. Their immediate problem is a
maximum matching whose deficiency `|V|-2|M|` equals the required hole count. Kekulization must also
move the aromatic-system charge back onto a selected exposed atom before deleting the system entry;
the current transformer deletes the system without performing that localization.

The general form is a degree-constrained matching with per-vertex demand in `{0,1}`:

- fixed `0`: prescribed hole;
- fixed `1`: must be covered;
- flexible `{0,1}`: a hole location chosen by the algorithm;
- an exact total matching cardinality/electron balance.

For the first implementation, the chemically important cases can use two simpler paths:

1. delete prescribed holes, then require a perfect matching on the remainder;
2. when holes are delocalized and their count equals the graph's matching deficiency, use a general
   maximum matching and take its exposed vertices as the localization sites.

Cases requesting more holes than the minimum deficiency, mixed prescribed/flexible holes, magnitude
greater than one, or open-shell systems should be rejected explicitly until a degree-constrained or
weighted matching formulation is designed. They should not silently accept an arbitrary maximal
matching.

## Recommended single-structure algorithm

Introduce a backend-neutral request rather than another kekulizer-specific algorithm enum:

```text
MatchingInput {
    required_covered: vertex set,
    required_exposed: vertex set,
    exposed_count: exact total,
}
```

For the initial closed-shell cases:

1. derive prescribed 0e/2e holes from `AromaticSystemAst::electrons`;
2. remove those vertices from the system graph;
3. if every remaining vertex is required, run general Edmonds and require a perfect result;
4. if the system has delocalized charge, run general Edmonds, verify the deficiency equals the
   requested hole count, and choose the deterministic result using canonical node/neighbor order;
5. map matched edges and exposed atoms back through the extraction correspondence;
6. write single/double bond orders and localize system charge/electron-pair state at exposed atoms;
7. validate the resulting localized molecule before removing the aromatic-system record.

Edmonds is the correct default because azulene, odd-ring fused systems, and C60 are non-bipartite.
The bipartite implementation is an optional fast path after its algorithm/name contract is resolved.
It is not necessary for the first correctness improvement.

## Enumeration: what Uno does and does not solve

The supplied papers support three distinct conclusions.

### Bipartite perfect and maximum matchings

Uno 1997 enumerates perfect and maximum matchings of a bipartite graph from an initial matching using
alternating cycles/paths and binary subproblems. The 1997 perfect-matching result is `O(V)` time per
output after initial matching; the improved 2001 algorithm gives `O(E sqrt(V))` preprocessing and
`O(log V)` amortized time per perfect matching, with `O(E+V)` space. The improvements depend on the
bipartite orientation `D(G,M)`, SCC trimming, careful branch-edge choice, and incremental state.
Implementing only the visible alternating-cycle recursion does not obtain the advertised bound.

This is directly applicable to benzenoids such as coronene after prescribed holes are deleted, but
not to general Kekulé graphs. In particular, C60 contains pentagons and is non-bipartite.

### Non-bipartite maximal matchings

Uno 2001 gives reverse-search enumeration of all maximal matchings in a general graph in
`O(E+V+Delta*N)` total time and `O(E+V)` space. Those outputs are inclusion-maximal, not necessarily
maximum, perfect, or compliant with a hole demand. Filtering its output for Kekulé structures can
have unacceptable delay because exponentially many smaller maximal matchings may occur between
valid outputs. It is scientifically interesting but is not the appropriate C60 Kekulé enumerator.

### General-graph Kekulé enumeration

For general graphs, retain binary partition/flashlight search but use an exact Edmonds extension
oracle. At a state with included/excluded edges:

1. reject conflicting included edges;
2. remove vertices covered by included edges and all excluded edges;
3. compute a maximum matching on the residual graph;
4. recurse only if `included_size + residual_maximum_size >= target_size` and all required-covered
   vertices can still be covered.

With a polynomial-height include/exclude tree and an exact polynomial extension oracle, this is a
polynomial-delay baseline. It will be slower asymptotically than bipartite Uno, but it is general,
much easier to validate, and reuses graph-core's Edmonds implementation. The current greedy oracle
should be replaced by this baseline before attempting Uno.

Enumeration should expose a lending iterator, callback, or visitor rather than returning only a
`Vec`. C60 has 12,500 Kekulé structures; eager collection obscures delay and imposes output-sized
memory even when a caller wants a prefix, count limit, or streaming fold.

## Independent counting without enumeration

For a planar aromatic-system graph, the number of perfect matchings can be computed independently
by the Fisher–Kasteleyn–Temperley (FKT) method. Given a planar embedding, construct a Pfaffian
(Kasteleyn) orientation and its signed skew-symmetric adjacency matrix. The absolute Pfaffian—or an
equivalent determinant computation—gives the number of perfect matchings, typically in `O(V^3)`
arithmetic time. It does not generate the matchings and does not depend on either branch-and-bound
or Uno.

FKT is not restricted to bipartite graphs. It therefore supplies an independent counting oracle for
both coronene and C60 even though Uno's bipartite perfect-matching algorithm does not apply to C60.
For planar bipartite graphs the same construction can be expressed as a determinant of a signed
Kasteleyn bipartite adjacency matrix; the unsigned permanent is not sufficient computationally.

Hole constraints reduce to related counts:

- prescribed holes `H`: compute the perfect-matching count of `G-H` once;
- one mobile hole: sum `PM(G-v)` over allowed vertices `v`; each near-perfect matching is counted
  once at its unique exposed vertex;
- a fixed small number `k` of mobile holes: sum `PM(G-H)` over allowed `k`-vertex subsets `H`;
- mixed prescribed/mobile holes: delete the prescribed set first, then sum over the allowed mobile
  subsets.

This remains independent of matching enumeration, although mobile-hole counting enumerates hole
sets. It is polynomial for fixed small `k` (`O(V^k)` Pfaffian evaluations naively), but does not make
the general monomer–dimer problem easy. Counting matchings of arbitrary cardinalities, or allowing
an unbounded number of mobile holes, is #P-hard even in important planar settings.

For a nonplanar graph, exact independent formulations still exist but lose FKT's polynomial bound:

- the permanent counts perfect matchings in a bipartite graph;
- the hafnian of the symmetric adjacency matrix counts perfect matchings in a general graph;
- deletion–contraction and treewidth dynamic programming compute matching-polynomial coefficients.

These are exponential in the general case (or exponential in width), but remain useful validation
oracles for small nonplanar fixtures. A Tutte-matrix determinant is useful for testing existence; it
must not be mistaken for an exact count because symbolic terms/cancellations carry more information
than an ordinary numeric determinant.

FKT counts labeled matchings. Counting only symmetry-inequivalent Kekulé structures is a different
problem: use the automorphism group with Burnside's lemma and count matchings fixed by each group
element (or treat symmetry reduction as a post-enumeration layer). The benchmark oracle must state
which count it expects.

## Experimental sequence

### Experiment A — correctness baseline

The detailed staged implementation plan is in discussion 146.

- Replace the unsound greedy extension test with exact residual Edmonds.
- Cross-check all graphs through a practical small-graph bound against exhaustive edge-subset
  enumeration.
- Assert uniqueness, matching validity, target cardinality, prescribed holes, and required coverage.
- Add relabeling invariance of the solution set (compare canonical edge sets, not output order).
- Add early-stop streaming tests.
- For planar fixtures, compare the number of distinct enumerated outputs with an independently
  implemented FKT count; retain exhaustive edge-subset counts only for small graphs.

This baseline is the reference against which Uno is compared, not the existing branch-and-bound.

### Experiment B — single kekulization with holes

Use at least:

- benzene and pyridine: no holes, perfect matching;
- pyrrole, furan, thiophene: one prescribed 2e heteroatom hole;
- borepin/boratabenzene: one prescribed 0e boron hole;
- Cp− and tropylium: one delocalized charged hole, with charge localized in the output;
- azulene: non-bipartite perfect matching;
- fused heterocycles with both prescribed and ordinary vertices;
- invalid electron vectors, impossible prescribed holes, unsupported multiple/flexible holes, and
  open-shell systems: specific errors.

Round-trip checks must compare electron accounting and charge, not demand restoration of the same
aromatic-system identity.

### Experiment C — enumeration performance

Benchmark the exact-oracle baseline and bipartite Uno separately. Report preprocessing, time to first
result, median/p95/max inter-output delay, total time, peak memory, and output count. A total-runtime
benchmark alone is insufficient for an output-sensitive algorithm.

Corpus:

- benzene (2), naphthalene (3), linear acenes, phenanthrene, pyrene, and coronene (20 expected);
- representative prescribed-hole heterocycles and charged rings, enumerating both allowed hole
  locations and matchings where applicable;
- azulene and other nonalternant systems;
- C60 (12,500 expected) as the principal non-bipartite/output-volume stress case;
- disconnected repeated aromatic components to expose Cartesian-product behavior;
- synthetic grids/ladders and dense bipartite graphs to test delay independently of chemistry.

Use FKT as the primary independent count oracle for planar cases, including prescribed-hole vertex
deletions and C60. Independently sourced known counts and exhaustive subsets on small graphs remain
cross-checks for the FKT implementation itself. Coronene and C60 topology fixtures already exist in
the repository and should be shared rather than copied into each benchmark module.

### Experiment D — Uno prototype gate

Implement the 1997 bipartite alternating-cycle algorithm first, preserving an initial matching and
incremental include/exclude state. Compare its solution set against the exact-oracle baseline before
attempting the 2001 SCC trimming/amortized refinements. Proceed to the 2001 algorithm only if the
1997 prototype is correct and benchmark evidence shows its delay matters on realistic benzenoids.

Do not present the non-bipartite maximal-matching paper as a C60 fallback. A general-graph Uno-like
Kekulé enumerator would require a separate derivation around blossoms or another exact extension
oracle.

## Correspondence integration

The kekulizer currently reconstructs host/sub maps even though `extract` already supplies a
`GraphCorrespondence`. When Experiment A or B lands, add one small transport operation at the
graph-core layer: map a `Matching` on the left graph through a `GraphCorrespondence` to its right-side
matched edges and exposed-node images (or construct the corresponding right-side `Matching` when a
right graph is supplied). The API must fail if a matched edge or required node is exposed by the
correspondence; silently dropping it would corrupt a matching.

Use that operation in kekulization and remove `host_to_sub`/`host_bonds`. Do not make matching itself
a `Correspondence<NodeId>`: an undirected self-matching and a left↔right partial bijection have
different semantics, and conflating them would make exposure and cardinality confusing.

## Decisions and open questions

Recommended decisions:

- Treat prescribed-hole kekulization as perfect matching after vertex deletion.
- Treat the common delocalized one-hole ions as maximum matching plus explicit charge localization.
- Use Edmonds for the first general implementation; keep bipartite dispatch as an optimization.
- Repair exact enumeration before evaluating Uno.
- Evaluate Uno only for bipartite perfect matching; do not filter general maximal matchings.
- Make enumeration streaming and bounded-memory.
- Add planar Pfaffian counting as an independent validation/counting path, not as an enumeration
  backend.
- Land correspondence transport with the first algorithmic change.

Questions to settle before an implementation plan:

1. For system charge magnitude greater than one, is every unit of charge intended to create one
   exposed 0e/2e site, or can the aromatic electron vector encode a different localization demand?
2. Should a canonical single Kekulé form choose the exposed site jointly with matching edges under
   the existing canonical node order, or should chemistry/HMO bond-order scores rank candidates?
3. Is the first enumeration API required to enumerate hole placements and bond matchings together,
   or only bond matchings for a fixed hole set?
4. Should output canonicalization quotient symmetry-equivalent Kekulé forms, or enumerate every
   labeled matching and leave symmetry reduction to a separate layer?

## Sources reviewed

- Takeaki Uno, *Algorithms for Enumerating All Perfect, Maximum and Maximal Matchings in Bipartite
  Graphs* (1997).
- Takeaki Uno, *A Fast Algorithm for Enumerating Bipartite Perfect Matchings* (2001).
- Takeaki Uno, *A Fast Algorithm for Enumerating Non-Bipartite Maximal Matchings* (2001).
- Takeaki Uno, *A New Approach to Efficient Enumeration by Push-out Amortization* (2014).
- Conte et al., *Designing Output Sensitive Algorithms for Subgraph Enumeration* (2025).
- Fisher–Kasteleyn–Temperley planar dimer/perfect-matching counting via Pfaffian orientation.
- Existing umol discussions 53, 58, 85, 87, and 93 and the current matching/kekulizer code.
