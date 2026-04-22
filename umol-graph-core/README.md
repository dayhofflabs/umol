# umol-graph-core

CSR (compressed sparse row) graph and relation-set primitives underlying the chemistry layers.

Three data structures, all `Arc`-shared with copy-on-write mutation:

| Type | Purpose |
|---|---|
| `Graph` | undirected graph, topology only (no edge data) |
| `FixedRelationSet<R, N>` | N-ary relation with compile-time arity (e.g. dative bonds: N=2) |
| `VarRelationSet<R>` | N-ary relation with runtime arity (e.g. aromatic systems) |

The chemistry layer stacks edge/relation data parallel to these structures; umol-graph-core itself holds only `NodeId` / `EdgeId` / `RelationId` identifiers.

## Identifiers

- `NodeId(u32)` — dense positional index; stable across mutations only when accompanied by a `Remapping`.
- `EdgeId(u32)` — dense positional index into each graph; semantics as above.
- `RelationId(u32)` — dense positional index per relation set.

`Remapping` is returned from destructive mutations (`Graph::remove`, relation-set removal). Consumers that hold ids across mutation MUST translate them through the returned `Remapping`; there is no stability guarantee otherwise.

## Graph invariants

**Storage.** `Graph` stores a `Topology` carrying
- `offsets: Vec<u32>` (length `node_count + 1`) — CSR prefix sum.
- `neighbors: Vec<Neighbor>` — flat CSR of `(NodeId, EdgeId)` entries.
- `endpoints: Vec<[NodeId; 2]>` — one pair per edge.

**Guaranteed:**

- `node_count` and `edge_count` are the lengths of the respective id domains.
- `offsets` is monotonically non-decreasing; `offsets[0] == 0`, `offsets[node_count] == neighbors.len()`.
- For every node `n`, `neighbors[offsets[n]..offsets[n+1]]` is sorted ascending by `Neighbor::node`. Consumers MAY use `binary_search_by_key(&n, |nb| nb.node)` for O(log deg) lookup.
- `edge_endpoints(eid)` returns `[a, b]` with `a ≤ b`. Self-loops return `[n, n]`. Canonicalization happens once at construction; consumers MAY rely on min-first for dedup patterns and for canonical equality on structures parallel to `EdgeId`.
- Every edge appears in the neighbor slices of both of its endpoints (undirected representation). Self-loops appear once.
- Parallel edges are permitted; they appear as repeated entries in the neighbor slice and as distinct `EdgeId`s.
- `find_edge(a, b)` returns `Some(e)` iff `a` and `b` are adjacent; on multiple parallel edges, any one of them MAY be returned.

**Not guaranteed:**

- The global order of edges beyond construction order.
- That `Remapping` preserves dense packing in any specific way beyond what its methods document.

## Relation-set invariants

Both `FixedRelationSet<R, N>` and `VarRelationSet<R>` carry
- `data: Vec<R>` parallel to the `RelationId` domain.
- an incidence CSR built once at construction.

**Guaranteed on both:**

- `RelationId(i).index()` is the position in `data`; `contains(id)` is a bounds check.
- `incident(node)` returns a slice of `RelationId`s sorted ascending (they are drawn from a single CSR segment, and segments are built in relation-id order).
- `has_incident(node)` is O(log N) on the flat incidence-node array.
- `data_mut(id)` and participant accessors are O(1) index lookups.

**`FixedRelationSet` specifics:**

- Participants stored as `[NodeId; N]` per relation, **sorted ascending by `NodeId` at construction**. `participants(id)` returns the sorted array; accessors MAY assume `participants[0] ≤ … ≤ participants[N-1]`.
- Duplicates in the caller's input are preserved (no dedup).
- Directional consumers — dative bonds are the sole example — carry direction in their `R` payload, not in participant order. `R`'s constructor is responsible for flipping the direction marker when the caller-supplied pair required a swap.
- `data` alone determines `PartialEq`; canonical participant order makes structural equality between sets with permuted input order hold automatically.

**`VarRelationSet` specifics:**

- Participants of one relation are stored in a sliced `Vec<NodeId>` framed by `offsets`, **sorted ascending by `NodeId` at construction** (`new` sorts the caller's vec in place before copying). Accessors return a sorted slice; duplicates in the caller's input are preserved.
- A node appears multiple times in one relation iff the caller's input did.

## Algorithms

Each algorithm in `algorithms/` is an enum dispatch over named implementations. Consumers pick an algorithm variant explicitly; there is no silent default. See the per-module docs. Algorithms operate on `Graph` alone; relation sets are not traversed by the algorithms module.

## Mutation semantics

- `Graph` mutations rebuild the `Topology` via `Arc::make_mut` semantics: other clones continue to see the old topology. There is no in-place edit path.
- Relation sets currently have no public mutation API beyond `data_mut`; structural edits go through the chemistry layer's builders and produce new relation sets.

## What this crate does not do

- No chemistry-specific types. `R` is opaque.
- No edge data on `Graph`; edge properties live alongside `EdgeId` in consumer structures.
- No string parsing, no serialization. `Debug` is available; no `Display`, no EDN.
- No thread-level concurrency primitives beyond `Arc`. All types are `Send + Sync` when `R` is.
