# Graph Library Decision

## Context

Step 7 of the unified AST migration (doc 80) requires ring enumeration, aromaticity perception, and neighbor iteration over `MoleculeAst`. Planning this step revealed that `MoleculeAst` has no native graph topology — `bonds: Vec<BondTuple>` is a flat list requiring O(|bonds|) scans for every neighbor query.

The indirection analysis in doc 82 ("Indirection analysis" section) traced the data path for four algorithms:

| Algorithm | Current indirection | With native topology |
|---|---|---|
| Ring enumeration | 4–5 conversions (petgraph → HashMap → BTreeMap → DenseProjection) | 0–1 (induced subgraph if filtering) |
| VF2 subgraph iso | 1 (build petgraph::Graph from flat bonds) | 1 (build petgraph) or 0 (native VF2) |
| Morgan/WL | 1 (build dense adjacency from flat bonds) | 0 (borrow slices) |
| DPO reactions | N/A | 0 for reads, new graph for output |

Every algorithm builds its own adjacency from scratch. `MorganTarget` already implements CSR internally. `AtomAdjacency` + `DenseProjection` exist solely to bridge petgraph's linked-list representation to the dense arrays the algorithms actually need.

## Three options considered

### Option 1: CSR in MoleculeAst

Add CSR (offsets + neighbors arrays) to `MoleculeAst`, built from `Vec<BondTuple>` at parse time. Bond attributes in parallel `IndexVec<BondIdx, BondAst>`.

- Fast to implement
- Good read performance, eliminates most indirection
- Immutable topology — attribute mutation (solver narrowing) works, but topology mutation does not
- Creates an asymmetry: mutable predicates, immutable topology
- Essentially "another table molecule with more pieces"
- No further optimization path

### Option 2: Keep petgraph

Continue using `StableGraph` in `MoleculeBuilder` and building `petgraph::Graph` in the matcher.

- `StableGraph` supports deletion but petgraph's VF2 requires `Graph` (not `StableGraph`) — can't have stable indices AND VF2 in the same representation
- Workaround: clone `StableGraph` → `Graph` for VF2
- Every graph algorithm extracts adjacency from petgraph into HashMap/BTreeMap/Vec<Vec<>> — petgraph stores adjacency internally but its linked-list representation is incompatible with what algorithms need
- Wall looming: any topology mutation conflicts with VF2's requirements
- petgraph contributes exactly one algorithm we use: vanilla VF2 (~2.7x slower than RDKit's VF2+)

### Option 3: New graph library crate

Build `Graph<N, E>` as a separate crate with native adjacency, stable indices, and algorithmic primitives.

- Significant upfront cost (~1 week)
- Solves the fundamental problem: one topology, all algorithms borrow from it
- Mutation + cache efficiency
- Can implement VF2+ directly, removing the petgraph algorithm dependency
- Similar precedent: EDN parser was a 1-week detour that became a solid foundation

## Decision: Option 3

New crate `umol-graph-core` (name TBD). Current `umol-graph` becomes `umol-discrete` (name TBD) — the discrete molecular structure layer (AST, solver, perception, matching, fingerprints).

## Mutation model

The argument for mutable topology was examined:

- **Solver**: never mutates topology. Narrows attributes only.
- **DPO reactions**: by definition produces a new graph. Input is immutable. The pushout constructs output from interface graph + RHS.
- **Reaction network construction**: append-mostly. Species and reactions are added; deletion is rare.

None of these require in-place topology mutation in the strong sense. But petgraph's `StableGraph`/`Graph` split (stable indices OR VF2, not both) is a real wall. And append-only growth is a legitimate mutation operation that CSR handles poorly.

The primary motivation for option 3 is not mutation per se but eliminating the indirection tax: one shared topology that all algorithms can borrow from without conversion.

## Design constraints

- `Vec<Vec<(NodeId, EdgeId)>>` adjacency (Graphs.jl style), not LEDA half-edge lists. Max degree ≤ 6 for molecular graphs; `Vec::retain` on the neighbor list is faster than any linked structure at bounded degree.
- Stable indices via free lists (reuse deleted slots).
- Do not target billion-node graphs. Molecular graphs (< 200 atoms) and reaction networks (< 10^6 species) are both in-memory with `Vec<Vec<>>`. Billion-node processing is a fundamentally different problem.
- Generic `Graph<N, E>` — node type N, edge type E. Not molecule-specific. Molecule-specific algorithms (Morgan, aromaticity, ring perception) stay in `umol-discrete`; generic graph algorithms (BFS, DFS, connected components, biconnected components, VF2) live in `umol-graph-core`.

## Scope

| Component | Estimate | Notes |
|---|---|---|
| `Graph<N, E>` struct + mutation | ~500 lines | add/remove node/edge, neighbor iteration, stable indices |
| BFS, DFS | ~100 lines | adapt from existing `algorithms/` |
| Connected components | ~50 lines | adapt from existing `algorithms/` |
| Biconnected components | ~100 lines | adapt from existing `algorithms/` |
| VF2 | ~500 lines | vanilla first, VF2+ later |
| `MoleculeAst` migration | ~300 lines | replace `Vec<BondTuple>` with `Graph` topology |
| Ring enumeration on `Graph` | ~200 lines | eliminates `DenseProjection`, `AtomAdjacency` |
| Morgan on `Graph` | ~100 lines | borrow adjacency directly |
| **Total** | **~1850 lines** | |

## RelationSet: typed hyperedge collections

`MoleculeAst` has multiple relation types beyond localized bonds: dative bonds, noncovalent bonds, aromatic systems, multicenter bonds. These are not graph edges — they don't need adjacency traversal. But they reference atoms by `NodeId` and must stay consistent when nodes are removed.

`RelationSet<R>` is a typed collection of relations (hyperedges) over a shared `NodeId` space. Each relation has an ordered participant list (`Vec<NodeId>`) and typed data (`R`). Per-node incidence lists enable O(incident) cascade removal when a node is deleted.

### Mapping to molecular relations

| Relation | participants | R (data) |
|---|---|---|
| Dative bond | `[donor, acceptor]` (ordered) | `BondAst` |
| Noncovalent bond | `[source, target]` (ordered) | `NoncovalentBondAst` |
| Aromatic system | `[atom0, ..., atomN]` | `AromaticSystemAst` (charge, spin, electron count) |
| Multicenter bond | `[atom0, ..., atomN]` | `MulticenterBondAst` (charge, spin, electron count) |

Participants carry the topology (who is connected to whom). `R` carries the chemistry (attributes of the relation). This mirrors the `Graph<AtomAst, BondAst>` separation where adjacency is topology and `N`/`E` are attributes.

### MoleculeAst composition

```
MoleculeAst {
    graph: Graph<AtomAst, BondAst>,              // primary topology with adjacency
    dative_bonds: RelationSet<BondAst>,           // directed binary, no adjacency
    noncovalent_bonds: RelationSet<NoncovalentBondAst>,
    aromatic_systems: RelationSet<AromaticSystemAst>,
    multicenter_bonds: RelationSet<MulticenterBondAst>,
    constraints: Vec<MoleculeConstraint>,
}
```

### Why not parallel arrays

The previous design stored each relation type as a `Vec<BondTuple>` or `Vec<AromaticSystem>` with atom indices. Node removal required manual cascade cleanup across every vec — one `retain` call per relation type, in user code. `RelationSet` internalizes this: `remove_participant(node)` drains the incidence list for that node, removes each incident relation, and cleans co-participants' incidence entries. No user-side bookkeeping.

### Incidence queries

- `aromatic_systems.has_incident(atom)` — O(1), replaces the current O(|systems| x |atoms_per_system|) scan in `is_in_aromatic_system`
- `aromatic_systems.incident(atom)` — returns the `RelationId`s for all systems containing that atom
- `dative_bonds.incident(atom)` — all dative bonds involving that atom

### Current AST type changes needed

`AromaticSystem { atoms: Vec<AtomIdx> }` splits: participants absorbed by `RelationSet`, data becomes `AromaticSystemAst` (charge, spin, electron count — to be defined). Same for `MulticenterBond`. `BondTuple { source, target, bond }` splits: source/target become participants, bond data stays as `BondAst`.

## What this unblocks

- Step 7 (aromaticity perception on solver) — ring enumeration and neighbor iteration work natively on `MoleculeAst`
- Future VF2+ implementation — no petgraph API constraint
- DPO graph rewriting — read from `Graph`, construct output as new `Graph`
- WL iteration, other graph algorithms — borrow adjacency directly

## Relationship to other docs

- Doc 80: step 7 suspended pending this work
- Doc 82: indirection analysis motivating this decision
- Doc 83: step 7 plan summary and suspension note
