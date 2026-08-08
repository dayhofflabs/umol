# Enumeration algorithm candidates

Status: **Informational**
Date: 2026-08-08
Relates: [167](167-graph-alg-execution-2026-07-27.md),
[158](158-ring-model-and-enumeration-2026-07-22.md),
[162](162-common-subgraph-algs-2026-07-25.md)

## Scope

Follow-up to [167](167-graph-alg-execution-2026-07-27.md). With the delivery
migration complete, this document inventories algorithm variants and related
algorithms not yet present in `umol-graph-core`, within the already-implemented
families and in adjacent families with defined chemical applications. Each
candidate carries its chemical application, its delivery classification under
the algorithm execution guide's streamability rule, and, where relevant, its
suitability for a resumable cursor (the deferred 167 S5 work).

Sources: the Yoshida–Fukuda enumeration table
(`materials/enumeration/Yoshida and Fukuda 1995 - Table of enumeration
algorithms.pdf`) and Wasa's maintained catalogue at
<https://kunihirowasa.github.io/enum/problem_list>.

Nothing here is scheduled. Selection happens per candidate when a consumer
appears; the one standing priority is noted below.

## Delivery classes

Assessed against the streamability requirement (result membership decidable at
emission; doc 167):

- **streams** — a leaf- or output-once enumerator; a `visit_*` form is direct.
- **ordered emission** — streams, and results arrive in cost order, each final
  when emitted; a stronger contract than the plain visitor.
- **two-phase** — streams after an exact bound is computed first (the
  `visit_maximum_matchings` construction: bound, then fixed-target
  enumeration).
- **no visitor** — incumbent branch and bound or single-result; eager only.

Reverse-search algorithms (Avis–Fukuda 1996) are additionally the most
cursor-amenable class: their entire state is the current object plus a local
successor rule, so a resumable `iter_*` needs no backtracking stack. If S5 is
ever picked up, a reverse-search operation is the cheapest first cursor.

## Candidates within implemented families

### Matchings

Implemented: Hopcroft–Karp, Edmonds, perfect-matching existence,
visit/enumerate for perfect and maximum matchings via branch and bound with
the exact Edmonds bound.

- **Uno enumeration** (Uno 1997): perfect/maximum bipartite matchings by
  alternating-cycle exchanges over one maximum matching; amortized per-output
  cost far below the current branch and bound, which re-runs Edmonds per
  branch node. A second `MatchingEnumerationAlgorithm` variant behind the
  existing visitors; the uniform delivery contract holds. **Streams.**
- **Maximal matchings** (Uno 1997): an operation the crate lacks entirely
  (maximal ≠ maximum; a different enumerator, not a filter — the same
  distinction doc 162 records for common subgraphs). **Streams.**
- **K-best weighted matchings** (Murty 1968; Chegireddy–Hamacher 1987):
  ranked Kekulé structures under bond-weight models. Requires edge weights,
  which `Graph` does not carry; a larger addition than a selector variant.
  **Ordered emission.**

### Cycles

Implemented: Read–Tarjan simple cycles, Vismara relevant cycles, Horton
minimum cycle basis, Kolodzik unique ring families, subdivision fallbacks,
shortest cycle through node/edge.

- **Directed elementary cycles** (Johnson 1975): nothing on `DiGraph` beyond
  topological order today. Cycle detection in reaction and metabolic network
  digraphs (futile and autocatalytic cycles). **Streams.**
- **Chordless cycles** (Uno–Satoh 2014): induced cycles — rings without
  transannular shortcuts; a ring-perception variant. **Streams.**
- **Hanser all-rings** (Hanser–Jauffret–Kaufmann 1996): the cheminformatics
  classic as an alternative `SimpleCycleEnumerationAlgorithm` selector.
  **Streams.**
- **Faster minimum cycle basis** (de Pina 1995; Berger et al. 2004):
  alternative `MinimumCycleBasisAlgorithm` selectors. Single result —
  **no visitor**, by the single-value rule.

### Paths

Implemented: bounded-length simple paths (all endpoint pairs, canonical
orientation), visitor and eager.

- **K-shortest simple paths** (Yen 1971; Eppstein 1998 for walks): pathway
  ranking. **Priority: this is the standing candidate for the planned
  reaction-network work** — ranked route enumeration over large reaction
  networks. Requires arc weights on the network graph, not the molecular
  graph. **Ordered emission**; Eppstein's structure is itself organized as a
  lazy iterator, so it is also cursor-amenable.
- **s–t paths** (Read–Tarjan 1975): endpoint-constrained variant of the
  existing enumeration (linker and chain analysis). **Streams.**

### Connected subgraphs

Implemented: connected edge subgraphs (ESU adaptation, bounded size), visitor
and eager.

- **Node-induced connected subgraphs** (Wernicke 2006 — the original ESU is
  node-induced; Avis–Fukuda 1996 by reverse search): atom-centered fragments,
  motif counting, Rücker-style subgraph-count descriptors (Rücker–Rücker
  2000). **Streams**; the reverse-search form has polynomial delay, bounded
  memory, no size cap, and is cursor-amenable.

### Maximal cliques and common subgraphs

Implemented: all cliques and Bron–Kerbosch with pivoting over the modular
product, direct backtracking, McGregor maximum common subgraph.

- **Worst-case-optimal pivoting** (Tomita–Tanaka–Takahashi 2006) and
  **degeneracy-ordered Bron–Kerbosch** (Eppstein–Löffler–Strash 2010):
  selector variants for the maximal-clique walk. **Stream.**
- **Koch c-cliques** (Koch 2001): Bron–Kerbosch modified to enumerate
  *connected* maximal common subgraphs directly — realizes the
  `McsConnectivity::Connected` axis that doc 162 identifies as the parameter
  chemists reach for first, currently absent from the maximal enumeration.
  **Streams.**
- **McSplit** (McCreesh–Prosser–Trimble 2017): maximum-common-subgraph
  variant beside McGregor. Incumbent branch and bound — **no visitor**
  (doc 167, streamability section).

### Independent sets

Implemented: one maximum independent set by branch and bound.

- **All maximal independent sets** (Tsukiyama–Ide–Ariyoshi–Shirakawa 1977):
  polynomial delay. Chemically: all Clar covers of the sextet-conflict graph,
  where the aromaticity code currently takes a single maximum set.
  **Streams.**
- **All maximum independent sets**: **two-phase** (independence number first,
  then fixed-target enumeration). For bipartite conflict graphs the bound
  phase is polynomial via König, so streaming is nearly free (Kashiwabara et
  al. 1992).

## Adjacent families with defined chemical applications

Not variants of implemented families; listed for completeness from the same
sources.

- **Spanning trees and arborescences** (Shioura–Tamura–Uno 1997;
  Gabow–Myers 1978): King–Altman enumeration of directed spanning trees of
  kinetic diagrams for steady-state enzyme rate laws (King–Altman 1956).
  **Streams**; reverse-search variants are cursor-amenable.
- **Minimal cut-sets** (Tsukiyama et al. 1980): minimal cut sets of
  biochemical reaction networks (Klamt–Gilles 2004). **Streams** (linear time
  per cutset).
- **Matching and independent-set counting** (Hosoya index; Merrifield–Simmons
  index): single-value counting operations, delivery-exempt like
  `count_perfect_matchings_planar`.

## References

- D. Avis, K. Fukuda. Reverse search for enumeration. Discrete Appl. Math.
  65 (1996) 21–46.
- F. Berger, P. Gritzmann, S. de Vries. Minimum cycle bases for network
  graphs. Algorithmica 40 (2004) 51–62.
- C. R. Chegireddy, H. W. Hamacher. Algorithms for finding K-best perfect
  matchings. Discrete Appl. Math. 18 (1987) 155–165.
- J. C. de Pina. Applications of shortest path methods. PhD thesis,
  University of Amsterdam, 1995.
- D. Eppstein. Finding the k shortest paths. SIAM J. Comput. 28 (1998)
  652–673.
- D. Eppstein, M. Löffler, D. Strash. Listing all maximal cliques in sparse
  graphs in near-optimal time. ISAAC 2010, LNCS 6506, 403–414.
- H. N. Gabow, E. W. Myers. Finding all spanning trees of directed and
  undirected graphs. SIAM J. Comput. 7 (1978) 280–287.
- T. Hanser, P. Jauffret, G. Kaufmann. A new algorithm for exhaustive ring
  perception in a molecular graph. J. Chem. Inf. Comput. Sci. 36 (1996)
  1146–1152.
- D. B. Johnson. Finding all the elementary circuits of a directed graph.
  SIAM J. Comput. 4 (1975) 77–84.
- T. Kashiwabara, S. Masuda, K. Nakajima, T. Fujisawa. Generation of maximum
  independent sets of a bipartite graph and maximum cliques of a circular-arc
  graph. J. Algorithms 13 (1992) 161–174.
- E. J. King, C. Altman. A schematic method of deriving the rate laws for
  enzyme-catalyzed reactions. J. Phys. Chem. 60 (1956) 1375–1378.
- S. Klamt, E. D. Gilles. Minimal cut sets in biochemical reaction networks.
  Bioinformatics 20 (2004) 226–234.
- I. Koch. Enumerating all connected maximal common subgraphs in two graphs.
  Theor. Comput. Sci. 250 (2001) 1–30.
- C. McCreesh, P. Prosser, J. Trimble. A partitioning algorithm for maximum
  common subgraph problems. IJCAI 2017, 712–719.
- K. G. Murty. An algorithm for ranking all the assignments in order of
  increasing cost. Oper. Res. 16 (1968) 682–687.
- G. Rücker, C. Rücker. Automatic enumeration of all connected subgraphs.
  MATCH Commun. Math. Comput. Chem. 41 (2000) 145–149.
- A. Shioura, A. Tamura, T. Uno. An optimal algorithm for scanning all
  spanning trees of undirected graphs. SIAM J. Comput. 26 (1997) 678–692.
- E. Tomita, A. Tanaka, H. Takahashi. The worst-case time complexity for
  generating all maximal cliques and computational experiments. Theor.
  Comput. Sci. 363 (2006) 28–42.
- S. Tsukiyama, M. Ide, H. Ariyoshi, I. Shirakawa. A new algorithm for
  generating all the maximal independent sets. SIAM J. Comput. 6 (1977)
  505–517.
- S. Tsukiyama, I. Shirakawa, H. Ozaki, H. Ariyoshi. An algorithm to
  enumerate all cutsets of a graph in linear time per cutset. J. ACM 27
  (1980) 619–632.
- T. Uno. Algorithms for enumerating all perfect, maximum and maximal
  matchings in bipartite graphs. ISAAC 1997, LNCS 1350, 92–101.
- T. Uno, H. Satoh. An efficient algorithm for enumerating chordless cycles
  and chordless paths. Discovery Science 2014, LNCS 8777, 313–324.
- S. Wernicke. Efficient detection of network motifs. IEEE/ACM Trans.
  Comput. Biol. Bioinform. 3 (2006) 347–359.
- J. Y. Yen. Finding the K shortest loopless paths in a network. Manage.
  Sci. 17 (1971) 712–716.
