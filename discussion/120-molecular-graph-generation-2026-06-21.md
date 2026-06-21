# 120 — Molecular graph generation: a palette for sampling and fuzzing (2026-06-21)

## Goal

A general-purpose facility for constructing molecular graphs from well-defined
random models. It is not specific to property tests: the same generators feed
algorithm fuzzing, benchmark-corpus construction, statistical studies, and fixture
generation. The immediate motivation is **molecular fuzzing** — driving the graph,
perception, and matching algorithms with large numbers of inputs drawn from a
*stated* distribution so that coverage and bias are known quantities rather than
artifacts of an ad-hoc generator.

The facility follows the same transparency rule as the rest of the codebase
(cf. the subgraph-isomorphism selector): not one opaque "use this" generator, but a
**palette of named, cited algorithms**. The caller chooses the molecule class and
the distribution shape; the facility exposes which algorithms realize that choice,
with each algorithm's guarantee (exact / asymptotic / Markov-chain) and cost stated.

## Scope of this document

**Topology only** — the bare graph skeleton (vertices and edges), no atom/bond
labels, no overlays, no stereo. Decoration (element, bond order, charge, hydrogens)
and the overlay/stereo layers are deferred; each is its own distribution problem and
its own palette, built on top of the skeleton layer. This document fixes the
skeleton layer and its architecture.

## Why MoleculeAst-as-data suits this

Generation is just construction of the same value every algorithm already consumes.
A generator has the shape `(params, &mut Rng) -> Graph` (lifted to `MoleculeAst` by
the decoration layer). Because the molecule is plain data, a generator needs no
coupling to any test harness, and generators compose: a skeleton sampler followed by
a decoration sampler followed by an overlay sampler. The proptest `Strategy` becomes
one *adapter* over a generator, not the generator itself.

## The three-axis decomposition

The design separates three concerns that ad-hoc generators entangle. Keeping them
orthogonal is what makes the distribution well-defined.

1. **Class** — the support set, expressed as structural constraints:
   - always: simple (no parallel edges, no self-loops), undirected;
   - connectivity (connected, or allow components);
   - maximum degree Δ (the molecular hallmark — the valence cap);
   - vertex count *n* (fixed or ranged);
   - edge count *m*, or equivalently the cyclomatic number *m − n + 1* (number of
     independent rings); acyclic (tree) is the *m = n − 1* corner;
   - optional further restrictions: *d*-regular, bipartite, girth bound.

2. **Distribution** — the probability measure over that class:
   - uniform over **labeled** graphs (vertices distinguishable);
   - uniform over **unlabeled** graphs (isomorphism classes);
   - a prescribed **degree sequence** (uniform among graphs with that sequence);
   - **Boltzmann** (size becomes a tunable random variable, exact-uniform given
     the realized size);
   - an explicitly parameterized model (e.g. a fixed-*p* edge model).

3. **Algorithm** — the concrete sampler realizing a (class, distribution) pair, with
   its guarantee and cost. Some pairs admit several algorithms with different
   tradeoffs; some admit none — that is surfaced honestly rather than approximated
   silently.

The caller states class + distribution; the facility offers the algorithms that
realize it. This is the transparency requirement made concrete: the algorithm is a
named choice, not a hidden default.

## Subtleties the facility must expose, not bury

- **Labeled vs. unlabeled uniformity.** Sampling labeled-uniform and forgetting the
  labels gives each isomorphism class probability proportional to *n! / |Aut(G)|*:
  asymmetric graphs are overrepresented, highly symmetric ones suppressed. For large
  random graphs almost everything is asymmetric and the distinction vanishes; for the
  *small, symmetric* graphs typical of molecules it is large. The choice must be
  explicit, because it determines which algorithms are even applicable (the cheap
  labeled methods vs. the canonical-form / enumeration methods).

- **Connectivity conditioning.** Most labeled-uniform methods can emit disconnected
  graphs. Above the connectivity threshold a constant fraction are connected, so
  rejection is cheap; below it, a connectivity-aware (block-decomposition) sampler is
  needed. The chosen handling is part of the algorithm's contract.

- **Degree cap is central, not optional.** The bounded maximum degree is what
  distinguishes molecular skeletons from generic graphs. The Erdős–Rényi family has
  an unbounded (Poisson) degree distribution and is a poor structural fit; rejecting
  down to bounded degree distorts it further. It is retained only as a labeled
  baseline, documented as such.

- **Planarity is not assumed.** Cage and cluster topologies are non-planar; the
  facility must not silently restrict to planar graphs.

- **Shrinking is a sampler property.** The proptest adapter needs a generator that
  shrinks toward minimal counterexamples. Not every sampler shrinks naturally; this
  is recorded per algorithm, and is a reason the tree-plus-augmentation construction
  is attractive for the fuzzing use.

## Algorithm palette (topology layer)

| algorithm | class served | distribution | guarantee | cost | source |
|---|---|---|---|---|---|
| Edge model G(n,p) / G(n,m) | all graphs / fixed-edge | labeled-uniform | exact | O(n²) / O(m) | Gilbert 1959; Erdős–Rényi 1960 |
| Configuration model + non-simple rejection | fixed degree sequence | labeled-uniform (asymptotic for bounded Δ) | asymptotic | O(m), bounded rejection | Bollobás 1980 |
| Double-edge-swap (switch chain) MCMC | fixed degree sequence, connectivity optionally preserved | labeled-uniform at stationarity | Markov chain | poly mixing (regular / bounded) | Fosdick et al. 2018; Greenhill 2015 |
| Random *d*-regular | *d*-regular | labeled-uniform | asymptotic | moderate *d* | McKay–Wormald 1990; Wormald 1999 |
| Prüfer-sequence tree | trees | labeled-uniform | exact | O(n) | Prüfer 1918 |
| Tree + *k* non-tree edges (under Δ-cap) | trees + fixed cyclomatic number | tree backbone exact; ring placement parameterized | exact backbone | O(n + k) | (composition of the above) |
| Recursive method | any decomposable class via a grammar | uniform, fixed size | exact | big-integer count tables, then ~O(n log n)/draw | Nijenhuis–Wilf 1978; Flajolet–Zimmermann–Van Cutsem 1994 |
| Boltzmann sampler | decomposable class | uniform given realized size, size tunable | exact (conditioned on size) | O(n) | Duchon–Flajolet–Louchard–Schaeffer 2004 |
| Canonical-augmentation / orderly generation | connected, Δ-bounded chemical graphs | unlabeled (isomorph-free) | exact, no duplicates | grows fast; small *n* | McKay 1998; Kerber–Laue (MOLGEN) |
| Enumerate-then-sample (nauty `geng`) | connected, Δ-bounded | unlabeled-uniform | exact for small *n* | exhaustive | McKay–Piperno 2014 |

Historical anchor: Pólya's enumeration theorem (Pólya 1937) was developed to count
chemical isomers — unlabeled structures under a symmetry group. Combined with the
recursive method it yields uniform generation over isomorphism classes, the
principled route to the unlabeled-uniform target.

## Architecture sketch (proposal — not yet a plan)

- **Layering.** The skeleton sampler produces a `Graph` and lives with the graph
  type (umol-graph-core owns `Graph`/`Csr`). A separate decoration palette lifts
  `Graph → MoleculeAst` by assigning elements/bond orders under a stated
  distribution — out of scope here but the reason the skeleton output is `Graph`,
  not a half-built `MoleculeAst`.

- **Dispatch.** Enum-of-algorithms mirroring `SubgraphIsomorphismAlgorithm`: each
  variant carries its parameters (e.g. degree sequence, *n* range, cyclomatic
  number, swap count). A class is a small constraint struct; the sampler validates
  that its parameters are consistent with the class and returns the algorithm's
  stated distribution.

- **Randomness.** A `rand::Rng` is injected so runs are reproducible from a seed; no
  global RNG.

- **Proptest adapter.** A thin wrapper turns a sampler into a `Strategy`, providing
  shrinking where the algorithm supports it (drop augmentation edges, then leaves,
  for the tree-plus-*k* construction). This is the only point the test framework
  enters; everything else is library code usable from any binary.

## Validation of the samplers

A sampler is only as trustworthy as the evidence that it hits its stated
distribution.

- For small *n*, enumerate the entire class with nauty `geng -c -D<Δ>` and run a χ²
  test of empirical isomorphism-class frequencies against the target (uniform → equal
  frequencies; labeled-uniform → proportional to *n! / |Aut|*, with |Aut| from
  nauty). nauty is a development-only, non-permanent oracle here — the same role the
  RDKit reference plays for the matchers, not a shipped dependency.
- Assert class membership (connected, Δ-bound, simple) on every draw.
- Report realized degree-sequence and ring-size histograms so the distribution is
  documented rather than assumed.

## Open decisions (require sign-off before implementation)

1. Output type and layering: `Graph` in graph-core with a separate decoration layer,
   vs. a single MoleculeAst-producing facility. (Leaning to the former — separation
   of skeleton from decoration.)
2. Whether labeled-vs-unlabeled is a first-class axis the caller selects, or a
   property fixed per algorithm.
3. Which algorithms seed the initial palette. Candidate starter set: the edge-model
   baseline; configuration-model + switch-chain MCMC (degree-sequence control);
   Prüfer tree + *k*-edge augmentation (ring control, shrinkable); plus the nauty
   validation harness. The recursive / Boltzmann and canonical-augmentation methods
   are larger and can follow.
4. Crate/module placement and the dev-only nauty dependency.
5. Shrinking contract for the proptest adapter.

## Out of scope (separate future palettes)

- Decoration: element, bond order, charge, implicit-hydrogen assignment over a given
  skeleton — its own distribution question (e.g. element frequencies, valence
  feasibility).
- Overlays: aromatic systems, dative/multicenter/noncovalent relations.
- Stereo configuration.

Each builds on the skeleton layer and is recorded separately when taken up.

## References

- Gilbert, "Random Graphs," *Ann. Math. Statist.* 1959.
- Erdős & Rényi, "On the Evolution of Random Graphs," 1960.
- Bollobás, "A probabilistic proof of an asymptotic formula for the number of
  labelled regular graphs," *Eur. J. Combin.* 1980 (configuration model).
- Fosdick, Larremore, Nishimura & Ugander, "Configuring Random Graph Models with
  Fixed Degree Sequences," *SIAM Review* 2018.
- Greenhill, "The switch Markov chain for sampling irregular graphs," *SODA* 2015.
- McKay & Wormald, "Uniform generation of random regular graphs of moderate degree,"
  *J. Algorithms* 1990; Wormald, "Models of random regular graphs," 1999.
- Prüfer, "Neuer Beweis eines Satzes über Permutationen," 1918.
- Nijenhuis & Wilf, *Combinatorial Algorithms*, 2nd ed., 1978 (recursive method).
- Flajolet, Zimmermann & Van Cutsem, "A calculus for the random generation of
  labelled combinatorial structures," *Theoret. Comput. Sci.* 1994.
- Duchon, Flajolet, Louchard & Schaeffer, "Boltzmann samplers for the random
  generation of combinatorial structures," *Combin. Probab. Comput.* 2004.
- Pólya, "Kombinatorische Anzahlbestimmungen für Gruppen, Graphen und chemische
  Verbindungen," *Acta Math.* 1937.
- McKay, "Isomorph-free exhaustive generation," *J. Algorithms* 1998.
- McKay & Piperno, "Practical graph isomorphism, II," *J. Symbolic Comput.* 2014
  (nauty/`geng`).
- Kerber, Laue et al., MOLGEN structure generator (orderly generation of chemical
  graphs).
