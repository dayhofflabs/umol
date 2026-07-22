# Kekulization algorithms

## Problem statement

Given a molecular graph where some bonds are marked aromatic (bond order 1.5 / unspecified), assign integer bond orders (single=1, double=2) such that each atom's valence is satisfied. This is equivalent to finding a perfect matching on the subgraph induced by aromatic bonds: matched edges become double bonds, unmatched edges become single bonds.

Three levels of capability:

1. **Single Kekulé structure** — find one valid assignment. Required for any toolkit that reads SMILES.
2. **Deterministic Kekulé structure** — the same molecule always produces the same assignment regardless of input atom/bond ordering. Required for canonical SMILES and InChI round-tripping.
3. **All Kekulé structures** — enumerate every valid assignment. Required for resonance energy calculations, Clar's rule, reaction mechanism enumeration, and e-graph representations.

## The graph theory problem

Kekulization is a perfect matching problem on the aromatic subgraph. The aromatic subgraph contains only the atoms and bonds in aromatic rings. A perfect matching M ⊆ E covers every vertex exactly once; matched edges become double bonds.

| Graph class | Single matching | Count all | Enumerate all |
|---|---|---|---|
| Bipartite (most aromatics: benzene, naphthalene, ...) | O(E√V) Hopcroft-Karp | O(V³) FKT / Pfaffian | O(V) per matching, Uno 1997 |
| General (azulene, odd-ring fused systems) | O(EV²) Edmonds blossom | #P-complete | O(V) per matching via blossom tree |
| Planar (all 2D molecular graphs) | O(E√V) | O(V³) FKT | Efficient via planar decomposition |

Most molecular aromatic systems are bipartite (alternant hydrocarbons). Non-bipartite cases arise from odd-membered rings (azulene, tropylium, corannulene) and certain heterocyclic fusions.

## Existing implementations

### RDKit (C++, BSD)

**Algorithm**: greedy DFS with backtracking.

- Start from each unmatched atom, try to assign a double bond to one of its aromatic neighbors
- If stuck, backtrack: undo assignments made since the last branch point, try the next candidate
- Runs in the aromatic subgraph only (pre-extracted)

**Determinism fix** (PR #9125, 2025): before DFS, sort candidate atoms by canonical rank. Process atoms in canonical order, sort neighbors by rank at each step. This ensures the same Kekulé structure is produced regardless of input ordering. The canonical ordering comes from Morgan/SMILES canonicalization, which is already available when Kekulization runs.

**Limitations**: produces only one structure. The backtracking search has exponential worst case but is fast in practice for molecular-sized graphs.

### CDK (Java, LGPL)

**Algorithm**: Edmonds blossom algorithm via John May's implementation.

- Handles general (non-bipartite) graphs correctly
- Produces a single maximum matching
- Used for both Kekulization and electron-pair assignment

### Open Babel (C++, GPL)

**Algorithm**: alternating path augmentation.

- Similar to RDKit but with a `FindPath()` function that searches for augmenting paths
- Greedy initial matching, then augmentation via path search

### InChI (C, MIT)

InChI does NOT kekulize in the traditional sense. Its normalization layer:

1. Identifies "mobile" and "fixed" hydrogen positions via tautomer perception
2. Represents delocalized bonds abstractly in the connection layer — the `/b` layer records double bond geometry, and the `/h` layer separates mobile vs fixed H
3. The canonical form deliberately avoids choosing a specific Kekulé structure for tautomeric systems

InChI's approach to aromaticity is fundamentally different from SMILES-based toolkits: instead of picking one resonance form, it encodes the entire equivalence class. This is closer to level 3 (all structures) but compressed into a canonical descriptor rather than explicit enumeration.

For non-tautomeric aromatic systems (e.g. benzene with no mobile H), InChI does assign double bonds internally during normalization, using an augmenting-path approach on the delocalized bond network. The specific implementation is in the normalization pass of the InChI generation algorithm (source: `INCHI_BASE/src/` in the IUPAC-InChI repository).

## Algorithms for single matching

### Greedy DFS with backtracking (RDKit-style)

```
kekulize(G_aromatic):
    M = {}
    for each atom a in canonical order:
        if a is unmatched:
            if not try_match(a, M):
                return FAILURE
    return M

try_match(a, M):
    for each aromatic neighbor b of a (in canonical order):
        if b is unmatched:
            add (a,b) to M
            return true
        if (b, c) in M and try_match(c, M \ {(b,c)}):
            remove (b,c) from M
            add (a,b) to M
            return true
    return false
```

Simple, fast for molecular graphs, but exponential worst case. Sufficient for level 1 and 2.

### Hopcroft-Karp (bipartite only)

O(E√V). Finds maximum matching by finding maximal sets of shortest augmenting paths simultaneously. Optimal for alternant hydrocarbons. Would need a bipartite check + fallback to Edmonds for non-bipartite cases.

### Edmonds blossom (general)

O(EV²). Handles odd cycles by contracting blossoms. Correct for all molecular graphs. Available as the `blossom5` library (Kolmogorov, 2009) or simpler implementations. Micali-Vazirani achieves O(E√V) for general graphs.

### Recommendation for single matching

Greedy DFS with backtracking, processing atoms in canonical order (Morgan rank or DFS discovery order on the canonical graph). This matches RDKit's proven approach, is simple to implement, handles molecular-sized graphs in microseconds, and gives determinism for free.

Edmonds blossom is the correct general algorithm but overkill for molecular graphs. If we encounter non-bipartite aromatic systems frequently, a simple blossom implementation (not Blossom V) would be worth adding.

## Algorithms for all Kekulé structures

### Uno's algorithm (bipartite, 1997–2001)

Enumerates all perfect matchings of a bipartite graph with O(V) delay per matching (1997) or O(log V) delay (2001). Total time O(E√V + N·V) where N = number of matchings.

**How it works**: Start from one perfect matching M₀. Build a directed graph D(G, M) where matched edges go one direction and unmatched edges go the other. Every other perfect matching M' differs from M₀ by a set of alternating cycles in D. The algorithm systematically finds these cycles by DFS, generating one new matching per cycle found.

This is the natural fit for molecular Kekulé enumeration:
- Most aromatic systems are bipartite
- The number of Kekulé structures is typically small (3 for naphthalene, 5 for anthracene, 20 for coronene)
- O(V) per structure means enumeration is effectively free for molecular graphs

### Adaptation for general graphs

For non-bipartite graphs, the alternating-cycle approach generalizes but requires blossom contraction. The delay per matching increases. In practice, non-bipartite aromatic systems in real molecules are rare enough that a slower fallback is acceptable.

### FKT counting (planar graphs)

If only the COUNT of Kekulé structures is needed (e.g. for Hosoya index, topological resonance energy), the FKT algorithm computes this in O(V³) via Pfaffian of the adjacency matrix, without enumerating the structures. All molecular graphs are planar, so FKT always applies.

### Recommendation for all structures

Uno's bipartite enumeration with blossom fallback for non-bipartite cases. The bipartite check is O(V+E) via 2-coloring. For the e-graph representation, each Kekulé structure becomes a node in the e-graph, with the aromatic system as the e-class.

## Implementation plan for umol-graph-core

### Phase 1: single deterministic Kekulé structure

New module `algorithms/matching.rs`:

- `fn perfect_matching(graph: &Graph, node_order: &[NodeId]) -> Option<Vec<EdgeId>>`
- Greedy DFS with backtracking, processing nodes in the given order
- Returns `None` if no perfect matching exists (non-kekulizable)
- The canonical node order comes from the caller (umol-graph's SMILES canonicalization or Morgan algorithm)

**Status (2026-04-29):** landed. `Graph::perfect_matching(node_order, PerfectMatchingAlgorithm::BacktrackingDfs) -> Option<Matching>` in `umol-graph-core/src/algorithms/matching.rs`. Returns the existing `Matching` struct rather than `Vec<EdgeId>` directly. Companion: `BipartitionAlgorithm::Bfs` in `algorithms/coloring.rs` and `MaximumMatchingAlgorithm::HopcroftKarp`. As of 2026-07-22, Hopcroft-Karp uses layered BFS plus batched DFS and reports non-bipartite input as `MaximumMatchingError::NonBipartite`.

### Phase 2: all Kekulé structures

- `fn enumerate_perfect_matchings(graph: &Graph) -> Vec<Vec<EdgeId>>`
- Bipartite check → Uno's algorithm or blossom-based enumeration
- Returns all perfect matchings
- For the e-graph: each matching is one term in the aromatic e-class

**Status (2026-04-29):** Uno's algorithm has been hard to extract from the original papers and from public reference implementations. Until that lands, `Kekulize::generate_all` returns a single-element iterator (the canonical structure produced by `BacktrackingDfs`). The trait shape stays uniform; enumeration capability lights up later. Alternatives if Uno remains unavailable: (a) backtracking enumeration (DFS that emits each completed matching, then backtracks — exponential worst case but fine at molecular sizes), (b) alternating-cycle enumeration from one matching (Uno's idea without polynomial-delay machinery; constructions in Lovász & Plummer, *Matching Theory*).

### Phase 3 (optional): matching count

- `fn count_perfect_matchings(graph: &Graph) -> u64`
- FKT algorithm for planar graphs (all molecular graphs)
- Useful for Hosoya index and resonance energy without full enumeration

## References

- Edmonds, J. (1965). Paths, trees, and flowers. *Canadian Journal of Mathematics*, 17, 449–467.
- Uno, T. (1997). Algorithms for enumerating all perfect, maximum and maximal matchings in bipartite graphs. *ISAAC*, LNCS 1350, 92–101.
- Uno, T. (2001). A fast algorithm for enumerating bipartite perfect matchings. *ISAAC*, LNCS 2223, 367–379.
- Kolmogorov, V. (2009). Blossom V: a new implementation of a minimum cost perfect matching algorithm. *Mathematical Programming Computation*, 1(1), 43–67.
- May, J. (2014). *Efficient cheminformatics algorithms and tools*. PhD thesis, University of Cambridge.
- RDKit PR #9125: Deterministic kekulize, independent of atom and bond order (2025).
- Kasteleyn, P. (1961). The statistics of dimers on a lattice. *Physica*, 27(12), 1209–1225.
