# 216 — Canonicalization performance

Status: Proposed
Date: 2026-08-30
Relates: [109](109-permutation-infrastructure-2026-06-09.md),
[110](110-molecular-symmetry-structure-2026-06-11.md),
[186](186-molecule-canonicalization-2026-08-05.md),
[208](208-canonicalization-scaling-2026-08-24.md)

## Purpose

Completed doc [208](208-canonicalization-scaling-2026-08-24.md) reduced complete canonicalization
of its retained feature-free reaction-network products from milliseconds to tens of microseconds by
selecting the lowest sufficient private description level and restoring sound automorphism-orbit
pruning. It did not re-run the ethane and propane network workloads because those sources and
reporting paths remain on the atom-mapping branch.

This document owns that deferred end-to-end measurement and the subsequent decision about whether
canonicalization needs another optimization stage. It also records candidate improvements that can
now be evaluated against the post-208 implementation. It does not select a backend, authorize a new
canonicalization carrier, or add an implementation plan before the measurements identify the
remaining cost.

## Inherited evidence

The pre-208 network profile reported:

| Case | Products canonicalized | Total | Product canonicalization |
| --- | ---: | ---: | ---: |
| Ethane | 855 | 0.368 s | 0.349 s |
| Propane | 17,929 | 143.173 s | 142.791 s |

Canonicalization therefore accounted for 99.73% of the propane run, at about 8 ms per produced
derivation. The final uninstrumented doc-208 benchmark measures its three retained feature-free
cases at about 47-96 us per complete canonicalization. This is a two-order-of-magnitude per-call
improvement on the retained workload, but only the network rerun can establish the resulting
end-to-end bottleneck.

Doc 208 also evaluated typed prefix pruning. The exact differential fixture reduced visited leaves
from 24 to 6 under the favorable branch order, but the production-shaped benchmark regressed by
12-15%. The guaranteed prefix became informative only after all atom images were fixed, at which
point constructing it duplicated most of the leaf-key work. Prefix pruning and its purpose-built
machinery were therefore removed. A future search design may reuse incrementally constructed key
state, but it must not reinstate the removed approach without new evidence.

## Current cost structure

The graph-IR canonicalizer performs a library-ordered typed individualization/refinement search.
The nauty result supplies automorphism orbits and generators for sound branch pruning and canonical
labels for branch order; the backend labeling convention does not define the accepted graph-IR
representative.

The topology-level `IncidenceGraph` has one node for every atom and localized bond and two incidence
edges per bond. `AutomorphismAdapter` already keeps unique ordinary bond-endpoint incidences as
direct edges, so it does not subdivide that topology a second time. Role- or value-bearing
incidences still require colored occurrence nodes at the vertex-colored nauty boundary.

The current graph-core nauty path constructs an owned CSR input for each call. The C shim allocates
and copies its partition, orbit, and sparse-graph buffers, requests a canonical graph unconditionally
with `options.getcanon = TRUE`, and frees those buffers after the call. Within one graph-IR
canonicalization, the adapter topology is fixed while successive search nodes vary the partition
colors. The final doc-208 feature-free cases visit one typed leaf but make five or six backend
calls, so backend setup and repeated stabilizer discovery remain concrete measurement targets.

## Required measurement

After the completed canonicalization branch is merged into the atom-mapping branch:

1. Re-run the ethane and propane extended-rule closures with the same seeds, rules, and limits as
   the inherited profile. Confirm identical flask, adjacency, transformation, and product counts
   before interpreting timings.
2. Record total network time, reaction application time, product-canonicalization time, number of
   canonicalization calls, and time per produced derivation.
3. Retain representative products from any newly dominant slow class and decompose their
   canonicalization into incidence construction, typed refinement and leaf-key work, backend calls,
   and final remapping. Record carrier sizes, residual cells, leaves, orbit-pruned branches, and
   backend-call counts.
4. Measure the backend itself separately from Rust-side adapter construction and result projection.
   Where a prototype changes the backend request or carrier, compare exact canonical aggregates
   and correspondence transport under dense renumbering before comparing time.

The first decision is whether the post-208 network performance is already sufficient. A larger
network case is justified only if propane no longer identifies the remaining bottleneck.

## Candidate optimizations

### Request only the backend result that the search consumes

Topology and constitution search require automorphism orbits; structure search additionally
filters projected generators. They do not require nauty's canonical graph to define the graph-IR
minimum. Full search currently disables orbit pruning, so without prefix pruning its backend calls
provide branch order only and cannot reduce the exhaustive leaf set.

Measure an automorphism-only nauty request that returns the required orbits and generators without
constructing a canonical graph. For a full search with neither orbit nor prefix pruning, also
measure deterministic local branch order with no backend call. These are private operational
changes and must preserve the exact typed minimum.

### Reuse fixed adapter and native storage

Measure a reusable backend session that retains the fixed CSR topology and capacity across the
partitions visited by one canonicalization. Partition colors, orbit output, and generator output
still change per call. This isolates allocation and copying from nauty's search cost without
changing graph semantics.

### Derive child stabilizers from one generated group

The search currently asks the backend for automorphisms of successive individualized partitions.
A generated permutation group and exact stabilizer chain could instead derive point stabilizers and
their orbits from a root generating set. This is a concrete consumer for the generated-group and
BSGS work anticipated by docs [109](109-permutation-infrastructure-2026-06-09.md) and
[110](110-molecular-symmetry-structure-2026-06-11.md). It should proceed only if repeated backend
calls remain material after the network and backend-mode measurements.

### Reduce the vertex-colored carrier when bond markers dominate

Molecule integrity makes a localized bond unique by its unordered atom endpoints. A compact exact
topology carrier can therefore represent one invariantly selected bond-value class as direct atom
edges and introduce a colored marker vertex only for bonds outside that class. For `n` atoms, `m`
bonds, and `m_nondefault` marked bonds, the candidate carrier has `n + m_nondefault` vertices and
`m + m_nondefault` edges instead of the topology incidence carrier's `n + m` vertices and `2m`
edges. A deterministic typed tie-break is required when several bond-value classes have the same
size.

This is a prototype question, not a selected replacement. It must preserve bond identity and every
selected bond distinction, integrate with the shared incidence facility rather than create a
second public molecular model, and be compared with the existing selectively subdivided adapter.
An edge-colored individualization/refinement backend remains the more general way to consume typed
incidences without vertex gadgets, but a less mature solver is useful only if the smaller carrier
wins end to end.

### Reduce canonicalization calls above the single-call algorithm

If one canonicalization is no longer the dominant cost but the network still spends materially on
the aggregate call count, measure reaction-application orbit reduction, pre-canonical duplicate
detection, canonical-result caching, and reuse of source refinement after local edits. These are
distinct from optimizing one canonicalization and must not make molecular identity depend on
derivation history.

## Decision boundary

Use the measured dominant term to select at most the next narrow optimization:

- backend canonical-graph or allocation cost: backend request modes or session reuse;
- repeated child backend calls: generated groups and stabilizers;
- adapter vertex count and backend search: compact or edge-colored carrier work;
- Rust typed search or leaf construction: revise that search directly;
- repeated network calls with acceptable single-call time: optimize the reaction-network path.

Any selected implementation must preserve complete canonical aggregates, canonical equality and
hash behavior, and operation-issued correspondence transport. Benchmarks and external solvers are
evidence, not runtime dependencies. Add a staged implementation plan only after the measurements
settle the target and its exact representation.
