# 216 — Canonicalization performance

Status: In Progress
Date: 2026-08-30
Relates: [109](109-permutation-infrastructure-2026-06-09.md),
[110](110-molecular-symmetry-structure-2026-06-11.md),
[186](186-molecule-canonicalization-2026-08-05.md),
[207](207-reaction-network-spike-2026-08-24.md),
[208](208-canonicalization-scaling-2026-08-24.md)

## Purpose

Completed doc [208](208-canonicalization-scaling-2026-08-24.md) reduced complete canonicalization
of its retained feature-free reaction-network products from milliseconds to tens of microseconds by
selecting the lowest sufficient private description level and restoring sound automorphism-orbit
pruning. It did not re-run the ethane and propane network workloads because those sources and
reporting paths remain on the atom-mapping branch.

This document owns that deferred end-to-end measurement and the evidence-gated implementation of
the next narrow general canonicalizer improvements. It does not authorize every candidate below:
each experimental path has an explicit selection or disposition gate, and broader backend,
generated-group, or reaction-network changes move to separate work only if their measurements
justify them.

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

## Measurement program

The investigation proceeds from representative use to the general operation, without treating the
reaction-network implementation as the optimization boundary.

1. Re-run several reaction networks under both the normal-polarity and extended rule catalogs.
   Include a range of molecule sizes and symmetry profiles, and record whether each run reaches
   closure or an explicit bound. For each run, report generation wall time, matching, application,
   product canonicalization, interning, other time, canonicalization calls, time per derivation, and
   canonicalization as a percentage of generation time. Do not interpret differences between the
   two catalogs as a controlled algorithm comparison unless the seeds and admitted chemistry are
   also comparable.
2. Retain representative products from the network runs and decompose the general molecule
   canonicalization operation into:
   - effective description-level selection;
   - incidence-graph construction;
   - normalized entity and incidence key construction;
   - color ranking and automorphism-adapter construction;
   - partition-descriptor construction and initial refinement;
   - recursive refinement, backend calls and projection, orbit selection, branch ordering, and
     leaf-key construction and comparison;
   - correspondence construction; and
   - final remapping and reframing, or complete-candidate materialization at the full level.
3. Record carrier sizes, residual cells, refinement and backend-call counts, visited leaves,
   orbit-pruned branches, and allocation-sensitive timings. Measure the backend separately from
   Rust-side adapter construction and result projection. The existing Criterion measurements of
   incidence construction, remapping, reframing, and complete canonicalization provide partial
   controls but do not replace measurements of the unmeasured middle stages.
4. Use that decomposition to evaluate optimizations of the general canonicalizer first. A
   reaction-network-specific reduction in canonicalization calls remains useful evidence and a
   possible later optimization, but the network is a training workload rather than the sole
   consumer or objective.

Where a prototype preserves the current typed order, compare exact canonical aggregates and
correspondence transport under dense renumbering before comparing time. A carrier may deliberately
define a different canonical numbering: canonical representatives and keys are deterministic for a
fixed release and context but are not stable identifiers under the current 0.x semantic-versioning
contract. Such a candidate instead must be exactly idempotent and invariant under valid dense
renumbering, preserve the represented aggregate through its returned correspondence, and induce the
same canonical-equivalence classes. Record whether it changes the current representative, but do
not reject it for that reason alone.

The network rerun first establishes whether canonicalization performance still blocks doc 207 and
then supplies evidence for the continuing optimization discussion. It does not by itself establish
that canonicalization performance is satisfactory. The catalog comparisons and corpus work in doc
207 provide additional representative workloads when broader or larger calibration cases are
needed.

## Post-208 network measurement

The merged atom-mapping branch was measured with the same extended carbon--hydrogen catalog, case
manifest, and bounds as the inherited profile:

```text
cargo run --release -p umol-reaction-network -- \
  --rule-catalog experimental/reaction-network/data/extended-carbon-hydrogen-rules.edn \
  --case-manifest experimental/reaction-network/data/extended-carbon-hydrogen-cases.edn \
  --case CASE --max-flasks 100000 --max-generations 64
```

Five release executions per case produced the following invariant counts:

| Case | Flasks | Directed adjacencies | Undirected adjacencies | Transformations |
| --- | ---: | ---: | ---: | ---: |
| Ethane | 101 | 402 | 201 | 855 |
| Propane | 1,230 | 8,760 | 4,380 | 17,929 |

Every run reached complete closure. The full reversibility checks also completed: ethane produced
432 identity and 423 source-automorphism reversals, while propane produced 9,420 identity and 8,509
source-automorphism reversals. Neither case had a missing reversal.

The median generation time and the millisecond-resolution internal phase counters were:

| Case | Generation | Matching | Application | Canonicalization | Interning | Other | Canonicalization per derivation |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Ethane | 0.052 s | 0.010 s | 0.004 s | 0.035 s | 0.002 s | 0.002 s | 40.9 us |
| Propane | 1.475 s | 0.129 s | 0.082 s | 1.181 s | 0.049 s | 0.031 s | 65.9 us |

The rounded internal counters need not sum exactly to the separately measured generation wall
time. Compared with the inherited profile, ethane generation is 7.1 times faster and product
canonicalization is 10.0 times faster. Propane generation is 97.1 times faster and product
canonicalization is 120.9 times faster. Canonicalization remains the largest propane phase, but its
share falls from 99.73% to about 80%, and the complete closure now takes about 1.5 seconds rather
than 143 seconds. Full runner time including the reversibility diagnostic is about 4.1 seconds.

No newly dominant slow product class is visible in the aggregate counters: the network-wide 65.9 us
per propane derivation lies within the tens-of-microseconds range represented by the retained
doc-208 cases. This does not complete measurement item 3. Canonicalization still takes about 80% of
propane generation, and the backend sample below leaves most of that time attributed only to the
combined Rust-side canonicalization path. The next investigation must decide whether retained
products from doc 207, more detailed phase instrumentation, or measurements above the single-call
canonicalizer best explain that remaining cost.

### Expanded network matrix

Exploratory release executions then applied the normal-polarity oxygen catalog and the extended
carbon--hydrogen catalog to the checked-in `network-cases.edn` seeds. These are workload
measurements, not a controlled comparison of the catalogs: their admitted electronic states and
element coverage differ. Every execution reached complete closure.

| Catalog | Seed | Flasks | Transformations | Generation | Canonicalization | Share |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Normal | `C2H4O` | 49 | 222 | 0.014 s | 0.008 s | coarse |
| Normal | `C3H6O` | 435 | 3,121 | 0.237 s | 0.189 s | 79.9% |
| Normal | `C4H6O` | 4,537 | 40,374 | 3.298 s | 2.701 s | 81.9% |
| Extended | `CH4O` | 6 | 22 | 0.004 s | 0.001 s | coarse |
| Extended | `C2H4O` | 51 | 312 | 0.020 s | 0.010 s | coarse |
| Extended | `C3H6O` | 1,057 | 12,135 | 0.837 s | 0.616 s | 73.6% |
| Extended | `C4H6O` | 10,472 | 159,043 | 12.098 s | 9.424 s | 77.9% |

The internal phase counters have one-millisecond resolution, so the smallest closures do not
support a stable percentage. All four closures longer than 0.2 seconds spend 73.6-81.9% of
generation time in canonicalization. The earlier extended propane result is therefore not an
isolated outcome, and both rule catalogs provide useful training workloads for the general
canonicalizer investigation.

## Backend sample

A one-second, one-millisecond-interval sample used the repository's release-equivalent `profiling`
profile so that the unchanged Rust and native call stacks remained symbolized. Sampling began
after input loading and covered 829 main-thread samples during propane generation. The measured
generation took 1.564 seconds under the profiler, compared with the 1.475-second unprofiled median.

The relevant inclusive sample counts are nested rather than additive:

| Scope | Samples | Share of canonicalization | Share of generation sample |
| --- | ---: | ---: | ---: |
| Product canonicalization | 681 | 100% | 82.1% |
| `AutomorphismAdapter::automorphisms_for_partition` | 124 | 18.2% | 15.0% |
| `umol_nauty_run` C shim and native call | 93 | 13.7% | 11.2% |
| `umol_nauty_sparsenauty` | 80 | 11.7% | 9.7% |

The 31 samples between the complete adapter path and the C shim cover Rust-side graph input
construction and result projection. The other 557 canonicalization samples remain in incidence and
key construction, typed partition refinement and leaf work, correspondence construction, remapping,
reframing, and their allocations. The backend is therefore a minority of both one canonicalization
and the complete network workload. Even eliminating the full native C call would bound the
generation improvement at about 1.13 times; eliminating the complete adapter path would bound it at
about 1.18 times.

## General canonicalization phase sample

A second profiling build temporarily prevented inlining at the private phase boundaries. A
five-second, one-millisecond-interval sample of the extended `C4H6O` generation attributed 3,154
main-thread samples to topology-level product canonicalization. The temporary boundaries increased
generation time by about 2%, so the counts are diagnostic proportions rather than replacement
Criterion results. The annotations were removed after sampling.

| Phase | Samples | Share of canonicalization |
| --- | ---: | ---: |
| Effective description-level selection | 2 | 0.1% |
| Incidence-graph construction | 36 | 1.1% |
| Normalized entity and incidence keys | 159 | 5.0% |
| Initial color ranking | 451 | 14.3% |
| Automorphism-adapter construction | 44 | 1.4% |
| Partition descriptors | 97 | 3.1% |
| Canonical search | 1,665 | 52.8% |
| Correspondence construction | 33 | 1.0% |
| Final remapping, reframing, and associated destruction | 667 | 21.1% |

The search samples break down further:

| Search work | Samples | Share of search | Share of canonicalization |
| --- | ---: | ---: | ---: |
| Initial partition construction from descriptors | 308 | 18.5% | 9.8% |
| Initial and recursive partition refinement | 724 | 43.5% | 23.0% |
| Automorphism calls and result projection, including nauty | 299 | 18.0% | 9.5% |
| of which the native nauty call | 194 | 11.7% | 6.2% |
| Typed leaf-candidate construction | 219 | 13.2% | 6.9% |
| Branch, orbit, and remaining search mechanics | 115 | 6.9% | 3.6% |

The nauty row is contained in the complete automorphism row; those two rows are not additive. The
other search rows partition the sampled search time.

The profile corresponds to visible allocation-heavy code paths. Initial colors clone the complete
recursive keys into a `BTreeSet`, move them into a `BTreeMap`, and then look up every original key.
Initial partition construction clones descriptors into another `BTreeMap`. Each refinement round
reconstructs dense cell indices, allocates a dense signature for every node, and groups each cell
through a `BTreeMap`. Final construction first remaps the complete molecule and then reframes the
owned intermediate. These paths, rather than incidence construction or the native backend alone,
are the first general-operation targets that the next controlled benchmarks need to discriminate.

A focused Criterion rerun keeps the retained complete canonicalizations at 95.2 us for the connected
case, 90.9 us for the disconnected case, and 46.8 us for the symmetry-heavy radical case. Remapping
the same cases takes only 3.28 us, 3.23 us, and 2.03 us respectively. Separate reframing controls
take 4.11 us for the overlay-heavy case and 10.6 us for a 128-participant aromatic system. The
reframing controls do not share the retained topology corpus and therefore cannot be subtracted
from complete canonicalization, but they confirm that the sampled final-construction term must not
be described as remapping cost alone. A fused-construction experiment needs a same-case breakdown
of remapping, normalization, frame selection, frame transport, and aggregate destruction.

## Settled execution strategy

The post-208 result removes the immediate 143-second propane blocker and demonstrates that the
completed work was effective. It does not establish that no further work is warranted.
Canonicalization remains the dominant measured generation phase at about 80%, and doc 207 already
provides concrete network workloads on which to continue the investigation. No hypothetical future
workload is required before this document can proceed.

The backend sample rules out only one broad interpretation: the native nauty call is not itself the
dominant remaining cost of the measured propane closure. The next work therefore starts with a
common controlled baseline and three independent investigation tracks: dense initial ranking and
partition construction, final-result construction, and a compact localized-bond carrier.

Those tracks are independent enough to prototype in parallel but not to interpret independently.
Carrier size changes the work performed by color ranking, partition construction, refinement, and
the backend. The compact-carrier and refinement results must therefore be compared separately and
together before either gain is treated as additive. Final remapping and reframing remain independent
of the search carrier and may proceed in parallel throughout.

The first compact-carrier prototype is private and topology-only. It compares the current carrier,
one exact ordinary-single-bond form represented by direct atom edges, and an invariant modal
bond-form class with a typed tie-break. Every other bond retains a marker. The public
`IncidenceGraph` remains unchanged during the experiment; a production carrier is selected only
after exact differential comparison and end-to-end measurement. Structure-level extension must
retain a marker for any stereo-bond site unless a later design proves an exact edge-site encoding.

All experimental worktrees start from one recorded commit and toolchain. Agents may edit separate
worktrees concurrently, but performance and allocation runs are serialized on the same machine.
No worktree may check out, reset, merge into, or otherwise mutate the protected `main` branch.
Unselected prototypes and profiling-only instrumentation are removed after their evidence is
recorded.

## Candidate optimizations

The phase sample gives the first three general candidates below the largest measured ceilings, but
does not yet select among them. The earlier backend and carrier candidates remain open where they
also reduce refinement or search work, but a wholesale backend replacement cannot be justified by
treating the native call itself as the dominant cost.

### Rank colors and construct the initial partition without cloning keys

Initial color ranking and initial partition construction account for about 24% of the sampled
canonicalization. Both currently clone ordered recursive keys into tree collections. Measure an
index-based ordering that sorts references to the existing keys, assigns dense ordered colors, and
forms partition cells from those colors. At topology level, partition descriptors can consume the
already ordered dense colors instead of cloning `InitialColorKey` values. Any generalization to
higher description levels must preserve their independently constructed descriptor order.

### Reduce partition-refinement allocation

Partition refinement accounts for about 23% of the sampled canonicalization, excluding initial
partition construction. Measure reusable cell-index and signature storage, sort-and-group
refinement in place of per-cell tree maps, and a worklist refinement that does not rebuild every
cell after every split. The exact ordered equitable partition is semantic search state: an
optimization must preserve its greater-signature-first order as well as its cell membership.

### Fuse final remapping and reframing

Final result construction accounts for about 21% of the sampled canonicalization. Measure a fused
internal construction path that writes normalized forms, remapped ids, representative participant
frames, and transported constraints into the final aggregate without first building a complete
remapped intermediate. This must implement the existing normalization--reframing--canonicalization
pipeline, not create a second representative or weaken the operation-issued correspondence.

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
derivation history. A cheap exact structural comparison or invariant screen before canonicalization
may suppress an important subset of redundant products, but it benefits this calling path rather
than the general canonicalization operation and is therefore not the first optimization target.

## Decision boundary

Use the measured dominant term and the carrier/refinement interaction matrix to select each narrow
optimization independently. The candidate list is not an implementation package:

- initial color and partition setup: reference-based dense ranking and partition construction;
- partition refinement: reusable storage or a worklist refinement with the same ordered result;
- final result construction: fused remapping and reframing;
- backend canonical-graph or allocation cost: backend request modes or session reuse;
- repeated child backend calls: generated groups and stabilizers;
- carrier size across initial setup, refinement, and backend search: compact or edge-colored
  carrier work;
- Rust typed search or leaf construction: revise that search directly;
- repeated network calls with acceptable single-call time: optimize the reaction-network path.

Any selected implementation must preserve complete canonical aggregates, canonical equality and
hash behavior, and operation-issued correspondence transport. Benchmarks and external solvers are
evidence, not runtime dependencies.

## Staged work plan

Every subitem below adds no public API unless it explicitly says otherwise. A prototype that would
require a new public graph-core or graph-IR seam stops at its disposition gate until that contract
is settled in this document or a linked follow-up. A stage ends with the active implementation
green; experimental worktrees may contain unselected alternatives until their disposition subitem.

### S0 — Establish the shared baseline and isolated workspaces

#### S0a — Freeze the semantic and performance corpus

**Module:** Rhea participant-corpus tooling, neutral fixture support, semantic tests, and a dedicated
Criterion target in `experimental/atom-mapping`; the existing CTfile boundary in `umol-io` and
faithful `MoleculeRecord` storage in `umol-store`; and self-contained constructed coverage in
`umol-graph-ir/benches/canonicalize.rs`; additive (green); no public API. [dep: none]

Freeze two explicitly different cohorts. The workload cohort retains the reaction-network products
that motivated this work and measures their topology path. The higher-description cohort supplies
deliberately selected constitution, structure, para-stereo, and full cases. Do not describe the
latter as reaction-network-derived or infer its frequency from the available reaction rules.

The starting complete-canonicalization benchmark is not sufficient by itself:

| Effective path | Existing complete cases | Coverage limit |
| --- | --- | --- |
| Topology | ordinary naphthalene, disconnected rings, and three retained network products | Good coverage of the current workload and topology scaling |
| Constitution | none | `overlay_heavy` selects Full; `large_aromatic` currently measures reframing only |
| Structure | tetrahedral stereo, meso dichlorobutane, and two para-stereo atom cases | No structure-level stereo-bond case and no imported chemical-size distribution |
| Full | `overlay_heavy` | One synthetic aggregate; its constraints are not frame-relative |

The incidence-construction groups force several `IncidenceLevel` values over this corpus, but that
does not exercise the corresponding complete effective canonicalization paths.

Use the ChEBI participant structures distributed with Rhea, not molecules reconstructed from the
reaction-SMILES source described by doc [203](203-atom-mapping-2026-08-19.md) and not the differently
standardized RXNMapper output. Rhea publishes one MDL Molfile per participant in
`rhea-mol.tar.gz`, the same structures in `rhea.sdf.gz`, and the small-molecule participant ids in
`chebiId_name.tsv`. This is the release-synchronized ChEBI subset that supplies Rhea's molecular
participants. It avoids reaction direction, stoichiometry, repeated reaction occurrences, and
lhs/rhs pairing while retaining the chemical population relevant to the reaction-network workload.
The complete ChEBI distribution is a possible later expansion, not part of the frozen baseline.

Freeze one Rhea release before the census and record the release identifier, source URLs, byte
sizes, and content fingerprints of `rhea-mol.tar.gz` and `chebiId_name.tsv`. Use the individual
Molfile records in stable ChEBI-id order so that one malformed participant does not prevent the
remaining records from being assessed. Join through `chebiId_name.tsv` and admit every
structure-bearing ChEBI participant to the census; do not pre-filter by size, connectivity,
perceived chemical validity, or expected parser support. One ChEBI participant is one source record
even when its represented graph is disconnected, as for an ion pair. Do not split such a record
into connected components.

For each participant, run the explicit CTfile parse--TableIR raise--resolution path under one
recorded chemistry model and `ResolveConfig`. Resolution is a named operation, not hidden inside
faithful ingestion. Record parse, raise, underdetermined, contradiction, execution, and unexpected
non-ground outcomes separately. Store every successfully resolved molecule exactly through
`MoleculeRecord`, preserving its entity ids and participant frames without normalization,
aromatization, or canonicalization. The source provenance is the Rhea release and ChEBI id; no
synthetic reaction or lhs/rhs pair is introduced.

Build a separate derived cohort by applying `Aromatizer` with
`AromaticityModel::daylight()` and `AromaticityConfig::default()` to each successfully resolved
participant. `Aromatizer` is a ground-to-ground transformation: it replaces the localized double
bonds of each perceived conjugated system with an aromatic-system overlay on a localized
single-bond scaffold and does not create `#a` assertions. Its output therefore selects Constitution
without a second resolution pass. Verify groundness, store the derived molecule beside its ChEBI
source provenance, and never mix its counts with the untransformed resolved cohort. Rhea does not
supply graph-IR constraints or the other overlay kinds systematically, so retain constructed cases
for dative, multicenter, noncovalent, para-stereo, and frame-relative full-constraint coverage.

The local `/Users/dr/Dayhoff/azoreductases/rhea-chebi-smiles.tsv` file may support an exploratory
dry run, but it is an older RDKit-derived beta representation and is not the frozen source unless
its release provenance and content fingerprint are selected explicitly. The Molfile corpus remains
the authoritative S0a input.

Select fixed ChEBI-participant cases only after the census. The combined benchmark surface must
include small, middle, and large resolved molecules for each populated effective path; at structure
level it must include atom stereo, bond stereo, multiple stereo sites, and a molecule carrying both
site kinds. At constitution level it must include one aromatic system, multiple or fused systems,
and a larger aromatic participant frame from the derived aromatized cohort. Preserve distinct ChEBI
records in census statistics when they resolve to equal canonical structures, but select named
imported cases by ChEBI id. Keep the existing para-stereo cases because imported stereo does not
establish that the para-stereo fixpoint changes a partition.

Record the source release and content fingerprints, repository commit, Rust toolchain, build
profile, and machine state. For each fixed case, record the exact canonical aggregate, verify that
its returned correspondence transports the source to that aggregate after reframing, and retain the
current dense-renumbering property tests as the semantic comparison surface. The S0a baseline table
reports complete canonicalization time and source graph and entity counts. S0c adds adapter sizes,
initial key and cell counts, refinement and search counters, phase timings, and allocation
measurements after the necessary private instrumentation exists. The existing canonicalization
unit, integration, and property tests must pass before timing begins.

**Evidence.** The frozen source is Rhea release 141 dated 2026-06-10:

| Source | URL | Bytes | XXH3-128 content fingerprint |
| --- | --- | ---: | --- |
| `rhea-mol.tar.gz` | `https://ftp.expasy.org/databases/rhea/ctfiles/rhea-mol.tar.gz` | 5,172,329 | `ef5cb65c0becec63f39da049e3639566` |
| `chebiId_name.tsv` | `https://ftp.expasy.org/databases/rhea/tsv/chebiId_name.tsv` | 589,145 | `2e530735f3d6e0b185a3db1e71929bb7` |

The builder uses `CtfileIoConfig::basic()`, the MDL valence model, the default stereo model and
`ResolveConfig`, followed independently by the Daylight aromaticity model and default aromatization
configuration. Its stable ChEBI-id census is:

| Population | Count |
| --- | ---: |
| Listed participants | 12,950 |
| Molfile records | 14,048 |
| Joined participants | 12,942 |
| Listed participants without a Molfile | 8 |
| Molfile records without a listed participant | 1,106 |
| Stored resolved records | 10,107 |
| Rejected listed participants, including missing Molfiles | 2,843 |
| Stored aromatized records | 10,107 |
| Aromatization failures | 0 |

The rejection ledger contains 8 source-read failures, 1,464 parse failures, 56 raise failures, 378
underdetermined resolutions, and 937 resolution contradictions. It contains no resolution-execution
or unexpected non-ground outcome. Doc [217](217-rhea-participant-failures-2026-08-30.md) owns the
independent cause analysis; these rejections are census results, not a requirement to turn S0a into
a complete Rhea ingestion project.

The populated effective paths before and after the separate aromatization transformation are:

| Cohort | Topology | Constitution | Structure | Full |
| --- | ---: | ---: | ---: | ---: |
| Resolved | 3,439 | 0 | 6,668 | 0 |
| Aromatized | 1,850 | 1,589 | 6,668 | 0 |

The fixed imported cases live under `experimental/atom-mapping/fixtures/canonicalization`, with
shared non-production metadata under `experimental/atom-mapping/support`. Their semantic tests and
Criterion target are also owned by `experimental/atom-mapping`. `umol-graph-ir` does not include or
load these files; its general benchmark retains self-contained constructed molecules only.

The imported complete-canonicalization baseline is in microseconds. Each entry is Criterion's point
estimate from 20 samples after a 0.5 s warm-up and 1 s measurement period; batched input cloning is
outside the measured routine.

| Effective path | Fixed case | Atoms | Bonds | Selected overlays | Without para stereo | With para stereo |
| --- | --- | ---: | ---: | --- | ---: | ---: |
| Topology | CHEBI:15379, O2 | 2 | 1 | none | 10.046 | 9.953 |
| Topology | CHEBI:2453, acyclovir | 16 | 17 | none | 73.083 | 73.230 |
| Topology | CHEBI:46245, ubiquinone-10 | 63 | 63 | none | 777.71 | 783.05 |
| Constitution | CHEBI:16150, benzoate | 9 | 9 | 1 aromatic system | 65.726 | 65.598 |
| Constitution | CHEBI:17097, biphenyl | 12 | 13 | 2 aromatic systems | 135.80 | 136.87 |
| Constitution | CHEBI:57306, protoporphyrin IX | 42 | 46 | 1 aromatic system, 24 participants | 325.38 | 326.17 |
| Structure | CHEBI:10983, (R)-3-hydroxybutanoate | 7 | 6 | 1 stereo atom | 81.413 | 81.853 |
| Structure | CHEBI:15903, beta-D-glucose | 12 | 12 | 5 stereo atoms | 256.95 | 257.57 |
| Structure | CHEBI:57287, CoA | 48 | 50 | 6 stereo atoms | 752.93 | 755.14 |

The release-141 participant set produces atom stereo but no stereo-bond entity. The self-contained
general benchmark therefore supplies the missing stereo-bond, mixed stereo-site, para-stereo, other
overlay, and full frame-relative-constraint paths:

| Constructed case | Atoms | Bonds | Relevant entities or constraints | Without para stereo | With para stereo |
| --- | ---: | ---: | --- | ---: | ---: |
| `cis_trans_stereo_bond` | 6 | 5 | 1 stereo bond | 52.017 | 52.239 |
| `mixed_atom_and_bond_stereo` | 11 | 9 | 1 stereo atom, 1 stereo bond | 110.48 | 110.66 |
| `frame_relative_stereo_constraint` | 5 | 4 | stereo-atom topicity in the form and molecule constraint store | 83.976 | 83.919 |
| `overlay_heavy` | 8 | 8 | all six overlay kinds and constraints | 164.08 | 164.10 |
| `para_stereo_trichloropentane` | 8 | 7 | 3 stereo atoms | 207.83 | 264.99 |
| `para_stereo_cascade` | 14 | 40 | 10 stereo atoms | 778.39 | 689.91 |

Measurements use canonicalizer base commit `6f234b92c9742213a75c3f7269f83688cf53953f`,
`rustc 1.96.0 (ac68faa20 2026-05-25)`, the optimized benchmark profile with debug information, and
an Apple M1 Pro MacBookPro18,1 with 10 cores and 32 GiB memory on AC power. Performance runs were
serial.

Focused verification passed: the corpus builder's 6 tests; 18 imported-corpus exact aggregate,
correspondence, and identity tests; 257 graph-IR canonicalization unit tests; 15 graph-IR
canonicalization integration tests; the 7 dense-remapping canonicalization properties; both
Criterion targets; and Clippy for both touched packages, all targets, and the property feature with
warnings denied.

**Done.**

#### S0b — Create the experimental worktrees

**Module:** repository operations; additive (green); no source or public API change. [dep: S0a]

Create detached or temporary-branch worktrees for dense setup, final construction, and compact
carrier work from the recorded baseline. Record their paths and commits in the working notes.
Agents do not edit the shared checkout concurrently. Do not create, check out, reset, merge, or
commit on `main`; selected work returns only to the active feature branch. Run all Criterion,
allocation, and reaction-network measurements serially.

**Working notes.** The isolated worktrees are clean detached checkouts of S0a commit
`ad6e5332cf53f4fc301957388df8629ec9e16db2` (`Add canonicalization corpus`):

| Investigation | Worktree |
| --- | --- |
| Dense setup | `/private/tmp/umol-216-dense-setup` |
| Final construction | `/private/tmp/umol-216-final-construction` |
| Compact carrier | `/private/tmp/umol-216-compact-carrier` |

The active feature checkout remains `/Users/dr/Dropbox/Source/rust/umol` on
`feature/atom-mapping`. No worktree is attached to `main`. Each investigation owns only its named
worktree, and performance and allocation measurements remain serial.

**Done.**

#### S0c — Establish phase and allocation measurements

**Module:** profiling-only graph-IR canonicalization instrumentation and doc 216; additive (green);
no public API. [dep: S0a, S0b]

Use symbolized sampling and the macOS Allocations instrument to record allocation count, allocated
bytes, and peak live bytes for initial key construction, color ranking, partition construction,
refinement, backend preparation and projection, leaf construction, correspondence construction,
remapping, normalization, frame selection and transport, and aggregate destruction. If exact phase
snapshots require a counting allocator or non-inlining boundaries, keep them in an isolated
profiling worktree and remove them after recording the evidence. Timing and allocation measurements
are separate runs.

**Evidence.** The symbolized timing half is the 3,154-sample general phase profile above. Allocation
measurement used an optimized test executable at S0a commit
`ad6e5332cf53f4fc301957388df8629ec9e16db2` in a fourth detached disposable worktree. The macOS
Allocations instrument was attempted first, but the local Instruments Devices plug-in reported a
missing weak symbol before launching the executable and produced an empty trace. The permitted
fallback therefore wrapped `System` with a profiling-only counting allocator and placed nested
scopes in the existing canonicalization path rather than copying the algorithm. Every case ran in
a separate process, outside Criterion, and the resulting canonical aggregate was compared with the
frozen corpus aggregate; the small Full case was compared with an uninstrumented result produced
before the counters were enabled.

Each table cell is `allocation or reallocation calls / gross allocated bytes / peak live-byte
increase from phase entry`. Counts and gross bytes exclude nested child scopes; the live-byte
increase includes values returned by a child that remain live in its parent. A reallocation counts
the new allocation and its complete new size. These are Rust-global-allocator measurements: native
nauty allocations are excluded, while Rust-side backend preparation and projection are included.

| Phase | Topology large | Constitution large | Structure large | Full constraint |
| --- | ---: | ---: | ---: | ---: |
| Normalization | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 40 / 4,390 / 3,060 |
| Incidence construction | 18 / 16,740 / 10,684 | 18 / 16,056 / 10,000 | 20 / 28,860 / 16,660 | 12 / 1,884 / 1,228 |
| Initial keys | 894 / 93,744 / 81,840 | 677 / 74,928 / 63,024 | 1,429 / 170,112 / 74,016 | 71 / 11,520 / 5,136 |
| Color ranking | 889 / 96,328 / 93,744 | 676 / 75,596 / 70,416 | 1,426 / 161,752 / 78,048 | 72 / 12,200 / 5,248 |
| Adapter construction | 35 / 13,812 / 12,716 | 31 / 16,680 / 13,136 | 67 / 30,824 / 15,420 | 28 / 2,880 / 1,488 |
| Partition construction | 1,915 / 150,844 / 75,600 | 1,263 / 106,280 / 53,856 | 4,018 / 367,404 / 143,328 | 441 / 43,348 / 22,516 |
| Refinement | 11,601 / 3,535,240 / 14,280 | 3,610 / 1,081,324 / 12,848 | 5,053 / 1,683,760 / 15,096 | 195 / 37,560 / 1,420 |
| Backend preparation and projection | 21 / 11,622 / 9,456 | 0 / 0 / 0 | 44 / 25,324 / 10,712 | 0 / 0 / 0 |
| Leaf construction | 1,023 / 103,004 / 85,328 | 744 / 83,376 / 66,584 | 909 / 96,408 / 79,024 | 78 / 8,980 / 10,036 |
| Correspondence | 25 / 3,776 / 2,924 | 22 / 2,700 / 2,032 | 24 / 3,112 / 2,304 | 36 / 1,104 / 224 |
| Remapping | 48 / 47,988 / 27,016 | 81 / 36,372 / 18,320 | 132 / 45,348 / 20,800 | 154 / 15,368 / 3,584 |
| Frame selection | 0 / 0 / 0 | 0 / 0 / 0 | 576 / 21,888 / 120 | 192 / 7,296 / 120 |
| Frame transport | 11 / 1,172 / 212 | 20 / 2,088 / 768 | 48 / 2,372 / 292 | 68 / 4,384 / 560 |
| Aggregate destruction | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 | 0 / 0 / 0 |
| **Measured total** | **16,480 / 4,074,270 / —** | **7,142 / 1,495,400 / —** | **13,746 / 2,637,164 / —** | **1,387 / 150,914 / —** |

Aggregate destruction freed 2,731 allocations and 258,720 bytes for Topology large, 1,970 and
197,420 bytes for Constitution large, 1,525 and 171,120 bytes for Structure large, and 127 and
17,144 bytes for Full constraint. It allocated nothing, so a live-byte increase is not applicable.

Refinement alone accounts for 86.8%, 72.3%, and 63.8% of measured gross Rust allocation bytes in
the three large cases. Partition construction plus refinement accounts for 90.5%, 79.4%, and 77.8%
respectively. The largest transient increases differ: color ranking leads the large Topology and
Constitution cases, while partition construction leads Structure and the small Full case. Frame
selection performs many short-lived allocations at Structure and Full. Remapping and frame
transport have modest gross byte volume even though the symbolized profile attributes 21.1% of
Topology time to final construction and associated destruction; allocation bytes therefore do not
stand in for timing, and S1b must still decompose that path. Constitution large and Full constraint
were already discrete after initial refinement and made no backend call, which explains their zero
backend rows rather than indicating that those description levels never use the backend.

The exact focused release probe passed for all four cases. The profiling-only source changes and
worktree were removed after these results were recorded.

**Done.**

### S1 — Run the independent prototypes

#### S1a — Prototype dense color ranking and initial partition construction

**Module:** `umol-graph-ir/src/ir/canonicalize.rs` and its module-local tests; additive (green); no
public API. [dep: S0b, S0c]

Replace recursive-key cloning through `BTreeSet` and `BTreeMap` with reference-index ordering that
assigns the same dense ordered colors. Form the topology initial partition from those colors without
cloning `InitialColorKey` values. Preserve the current implementation as a test-only comparator for
the duration of the selection stage.

Add focused `rstest` tables for empty, repeated, and ordered-distinct key sets and a property that
the prototype returns the same colors, ordered partition, and canonical aggregate under dense
renumbering and that its correspondence transports the source to that aggregate. Benchmark the
phase and complete operation on every S0 case.

**Prototype.** The detached dense-setup worktree remains based on S0a commit
`ad6e5332cf53f4fc301957388df8629ec9e16db2`; its uncommitted changes are isolated at
`/private/tmp/umol-216-dense-setup`. `rank_initial_colors` now sorts indices whose keys are borrowed
from the existing entity and incidence slices, then writes each ordered color rank directly into
the result vectors. The topology path constructs `OrderedPartition` directly from the adapter's
color ranks and removes ranks absent from that adapter. It does not construct or sort recursive
`InitialColorKey` partition descriptors. Constitution, Structure, and Full retain their typed
partition descriptors. The cloned-key ranking and topology descriptor path remain test-only
comparators.

The focused phase probe ran 2,000 iterations per case in separate release processes. Times are
microseconds per call; `—` means that the case does not select the topology partition path. Every
new result was compared with the retained reference before timing.

| Case | Color ranking, reference -> prototype | Speedup | Topology partition, reference -> prototype | Speedup |
| --- | ---: | ---: | ---: | ---: |
| `ordinary_naphthalene` | 4.982 -> 1.190 | 4.2x | 6.253 -> 0.217 | 28.8x |
| `disconnected_rings` | 5.546 -> 1.366 | 4.1x | 7.091 -> 0.218 | 32.5x |
| `overlay_heavy` | 7.132 -> 2.627 | 2.7x | — | — |
| `tetrahedral_stereo` | 2.911 -> 0.388 | 7.5x | — | — |
| `cis_trans_stereo_bond` | 3.766 -> 1.124 | 3.4x | — | — |
| `mixed_atom_and_bond_stereo` | 6.982 -> 2.437 | 2.9x | — | — |
| `frame_relative_stereo_constraint` | 2.841 -> 0.376 | 7.6x | — | — |
| `meso_dichlorobutane` | 3.769 -> 1.776 | 2.1x | — | — |
| `para_stereo_trichloropentane` | 5.123 -> 2.902 | 1.8x | — | — |
| `para_stereo_cascade` | 20.348 -> 8.712 | 2.3x | — | — |
| `feature_free_connected` | 7.406 -> 1.548 | 4.8x | 9.986 -> 0.193 | 51.7x |
| `feature_free_disconnected` | 7.369 -> 1.406 | 5.2x | 9.805 -> 0.218 | 45.0x |
| `symmetry_heavy_radicals` | 3.571 -> 0.589 | 6.1x | 5.078 -> 0.126 | 40.3x |
| `topology_small_chebi_15379` | 1.245 -> 0.158 | 7.9x | 1.554 -> 0.062 | 25.1x |
| `topology_middle_chebi_2453` | 14.420 -> 4.971 | 2.9x | 15.620 -> 0.422 | 37.0x |
| `topology_large_chebi_46245` | 57.491 -> 22.971 | 2.5x | 61.554 -> 0.883 | 69.7x |
| `constitution_one_system_chebi_16150` | 8.079 -> 2.872 | 2.8x | — | — |
| `constitution_multiple_systems_chebi_17097` | 11.612 -> 5.313 | 2.2x | — | — |
| `constitution_large_frame_chebi_57306` | 43.639 -> 15.700 | 2.8x | — | — |
| `structure_one_atom_site_chebi_10983` | 6.312 -> 1.471 | 4.3x | — | — |
| `structure_multiple_atom_sites_chebi_15903` | 12.394 -> 5.501 | 2.3x | — | — |
| `structure_large_chebi_57287` | 49.935 -> 16.057 | 3.1x | — | — |

The complete-operation Criterion comparison uses the frozen S0a baseline and reports relative mean
time; negative values are improvements. All confidence intervals exclude zero. The
`para_stereo_trichloropentane` result with para-stereo enabled is below Criterion's 1% practical
significance threshold; every other case is reported as an improvement.

| Constructed case | Para stereo off | Para stereo on |
| --- | ---: | ---: |
| `ordinary_naphthalene` | -16.49% | -15.50% |
| `disconnected_rings` | -12.14% | -12.15% |
| `overlay_heavy` | -5.07% | -5.07% |
| `tetrahedral_stereo` | -8.27% | -8.30% |
| `cis_trans_stereo_bond` | -9.30% | -9.71% |
| `mixed_atom_and_bond_stereo` | -7.24% | -7.54% |
| `frame_relative_stereo_constraint` | -7.06% | -5.33% |
| `meso_dichlorobutane` | -3.92% | -3.94% |
| `para_stereo_trichloropentane` | -3.66% | -0.76% |
| `para_stereo_cascade` | -2.82% | -2.62% |
| `feature_free_connected` | -18.24% | -20.05% |
| `feature_free_disconnected` | -18.78% | -19.55% |
| `symmetry_heavy_radicals` | -17.64% | -17.70% |

| Imported case | Para stereo off | Para stereo on |
| --- | ---: | ---: |
| `topology_small_chebi_15379` | -27.72% | -27.96% |
| `topology_middle_chebi_2453` | -35.48% | -35.64% |
| `topology_large_chebi_46245` | -10.36% | -10.74% |
| `constitution_one_system_chebi_16150` | -7.80% | -8.44% |
| `constitution_multiple_systems_chebi_17097` | -5.24% | -8.63% |
| `constitution_large_frame_chebi_57306` | -10.84% | -9.77% |
| `structure_one_atom_site_chebi_10983` | -13.12% | -14.31% |
| `structure_multiple_atom_sites_chebi_15903` | -7.59% | -6.04% |
| `structure_large_chebi_57287` | -11.80% | -10.62% |

Focused verification passed: 264 canonicalization unit tests, including the new three-case ranking
and partition tables and the dense-renumbering differential property; 15 canonicalization
integration tests; the 7 feature-gated molecule-canonicalization properties; 18 imported-corpus
tests; and graph-IR Clippy for all targets with the property feature and warnings denied. The
prototype exposes no public API and no correctness, phase, or complete-operation regression was
observed.

**Done.**

#### S1b — Decompose final-result construction

**Module:** graph-IR molecule remapping, aggregate reframing, and canonicalization profiling;
additive (green); no public API. [dep: S0b, S0c]

Measure remapping validation and correspondence conversion, atom and bond reordering, graph
reconstruction, overlay remapping, normalization, representative-frame selection, constraint frame
transport, final allocation, and destruction of the intermediate. Do not implement fusion in this
subitem. Record which work is duplicated by the consecutive `remap` and `reframe` calls and whether
the duplication remains material on topology as well as overlay-, stereo-, and constraint-bearing
cases.

**Evidence.** Profiling-only scopes and a counting allocator were added to the clean detached S1b
worktree, exercised through ignored white-box tests, and then removed. Each source was first
canonicalized normally; its returned correspondence was then passed through
`try_remap(...).reframe()`, and the result was compared exactly with the canonical aggregate before
measurement. The three large cases ran 10,000 final constructions and the three smaller cases ran
20,000, serially, in the release profile with para-stereo enabled.

The scopes add a diagnostic cost. Compared with the unchanged Criterion remapping controls, the
instrumented remap total was about 13% higher for ordinary naphthalene, 14% higher for the retained
network product, and 20% higher for the small constraint case. The following microseconds per call
and their shares of the frozen S0a complete-operation means are therefore conservative upper bounds,
not attainable speedups.

| Case | Effective path | Remap | Reframe | Result destruction | Total | Share of complete |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `feature_free_connected` | Topology | 3.752 | 1.766 | 0.490 | 6.008 | 6.1% |
| `ordinary_naphthalene` | Topology | 3.809 | 1.722 | 0.503 | 6.034 | 9.3% |
| `topology_large_chebi_46245` | Topology | 12.655 | 5.855 | 1.018 | 19.529 | 2.5% |
| `constitution_large_frame_chebi_57306` | Constitution | 11.876 | 5.146 | 0.927 | 17.948 | 5.5% |
| `structure_large_chebi_57287` | Structure | 15.913 | 10.254 | 1.033 | 27.201 | 3.6% |
| `frame_relative_stereo_constraint` | Full | 4.542 | 3.155 | 0.583 | 8.280 | 9.9% |

The remap decomposition is:

| Case | Validation | Correspondence conversion | Atom and bond reordering | Graph reconstruction | Overlay remapping | Constraint id transport | Result integrity |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `feature_free_connected` | 0.020 | 0.370 | 0.928 | 0.330 | 1.196 | 0.022 | 0.575 |
| `ordinary_naphthalene` | 0.020 | 0.366 | 0.915 | 0.382 | 1.197 | 0.022 | 0.599 |
| `topology_large_chebi_46245` | 0.020 | 1.428 | 4.309 | 1.104 | 1.211 | 0.023 | 4.226 |
| `constitution_large_frame_chebi_57306` | 0.020 | 1.127 | 3.115 | 0.869 | 2.421 | 0.023 | 3.927 |
| `structure_large_chebi_57287` | 0.020 | 1.354 | 3.463 | 0.931 | 4.021 | 0.022 | 5.727 |
| `frame_relative_stereo_constraint` | 0.020 | 0.257 | 0.560 | 0.298 | 2.160 | 0.069 | 0.851 |

The reframe decomposition is:

| Case | Action domain | Normalization | Representative-frame selection | Overlay frame transport | Constraint frame transport |
| --- | ---: | ---: | ---: | ---: | ---: |
| `feature_free_connected` | 0.025 | 0.999 | 0.113 | 0 | 0.024 |
| `ordinary_naphthalene` | 0.025 | 0.959 | 0.112 | 0 | 0.024 |
| `topology_large_chebi_46245` | 0.026 | 5.081 | 0.113 | 0 | 0.025 |
| `constitution_large_frame_chebi_57306` | 0.026 | 3.765 | 0.334 | 0.207 | 0.024 |
| `structure_large_chebi_57287` | 0.026 | 5.681 | 0.486 | 2.508 | 0.025 |
| `frame_relative_stereo_constraint` | 0.056 | 1.447 | 0.170 | 0.442 | 0.157 |

The difference between each total and its listed scopes is probe and control-flow overhead. Gross
Rust allocation counts and bytes close exactly against the S0c allocation measurements:

| Case | Remap allocations / bytes | Reframe allocations / bytes | Result deallocations / bytes |
| --- | ---: | ---: | ---: |
| `feature_free_connected` | 44 / 9,872 | 11 / 1,172 | 19 / 4,636 |
| `ordinary_naphthalene` | 45 / 9,836 | 11 / 1,172 | 19 / 4,552 |
| `topology_large_chebi_46245` | 48 / 47,988 | 11 / 1,172 | 19 / 20,764 |
| `constitution_large_frame_chebi_57306` | 81 / 36,372 | 20 / 2,088 | 24 / 15,324 |
| `structure_large_chebi_57287` | 132 / 45,348 | 48 / 2,372 | 24 / 17,680 |
| `frame_relative_stereo_constraint` | 77 / 7,684 | 34 / 2,192 | 25 / 3,228 |

The code and measurements correct the premise that reframing destroys a separately allocated
complete intermediate molecule. `Molecule::reframe` consumes the remapped molecule, mutates its
unique atom and bond arcs in place, and moves each overlay store through `mem::take`. There is one
result aggregate, whose ordinary destruction is reported above. The avoidable destruction is of
temporary vectors and relation stores inside the two operations, not of a second complete molecule.

The actual overlap and boundaries are:

- Atom and bond remapping clones each source vector and then constructs both an `Option`-wrapped
  source vector and an equally sized target vector. The target allocation is necessary; the
  temporary source representation and its traversal are not intrinsic to a known total bijection.
- Each overlay is first rebuilt by its relation-set `remap`, consumed into entries, reordered by
  entity id, and rebuilt again. Reframing then traverses the rebuilt store to normalize its payload,
  select the representative participant frame, and permute the store. The topology cases still
  spend about 1.2 us and 22 allocations / 2,344 bytes rebuilding six empty overlay stores during
  remap.
- Fused aggregate reframing performs another 11 allocations / 1,172 bytes even on those topology
  cases. These are empty replacement stores introduced by moving all six aggregates through
  `mem::take`, not chemically meaningful frame work.
- The operation-issued correspondence is materialized first as eight `Correspondence` values, then
  copied into graph-core `Remapping` vectors and an `IdRemapping`. This is representation conversion,
  not a second semantic transport. It grows with molecule size and is visible on the large cases.
- Public `try_remap` must validate an independently supplied correspondence and return an
  integrity-valid closed molecule. Canonicalization nevertheless sends it an operation-issued
  total correspondence and then pays the complete result integrity audit. That audit costs about
  0.6 us on the small topology cases and 3.9--5.7 us on the large higher-description cases.
- Id transport and participant-frame transport are not duplicate work. Constraints require both,
  in that order. Normalization and representative-frame selection are also semantically required;
  only their scheduling with cloning, remapping, and store construction can be fused. At Full, the
  source was already normalized before search, so the final normalization is additionally a known
  repeated reduction.

The avoidable work is real, but a canonicalization-only fused constructor has an unfavorable
boundary. The complete measured final path, including required result allocation, normalization,
frame work, and inflated probe overhead, is only 2.5--9.9% of the complete operation on these cases.
A fusion could remove only part of that ceiling, while duplicating the construction and integrity
logic of `Molecule::try_remap` across atoms, bonds, six overlay kinds, and recursive constraints.
S1b therefore does not support the broad fused path proposed by S2c. The narrower opportunities are
general remapping/reframing improvements: direct target construction for a proven total bijection,
empty-overlay fast paths, cheaper move-out of empty aggregate stores, and a carefully delimited
operation-issued integrity path. They should not be smuggled into a private canonicalizer copy.

The profiling source and ignored measurement tests were removed. The S1b worktree is again exactly
the clean S0a commit and exposes no public API.

**Done.**

#### S1c — Prototype a compact topology carrier

**Module:** private graph-IR canonicalization carrier/adapter code, module-local tests, and the
canonicalization benchmark; additive (green); no public API. [dep: S0b, S0c]

Keep the public `IncidenceGraph` as the semantic source. Prototype a private topology search carrier
that represents the selected exact normalized bond-form class as direct atom edges and every other
bond as a colored marker vertex. Compare both the complete ordinary neutral single-bond form and the
invariant modal exact bond-form class; resolve equal modal counts by the frozen typed key order.
Derive bond images and bond correspondence from mapped unordered atom endpoints. Keep leaf-key
construction and final molecule construction unchanged.

Add differential tests against the current carrier for all bounded simple graphs already used by
canonicalization, mixed single and non-single bonds, electronically distinct order-one bonds,
disconnected graphs, and dense entity renumberings. Each case compares the complete canonical
aggregate across the original input, dense renumberings, the current representative, and the
candidate representative; it also checks exact idempotence and the correspondence transport law.
Equality with the current representative is recorded but is not required. Record carrier sizes,
partition statistics, allocations, and phase and complete-operation timings.

**Evidence.** The isolated prototype started from commit
`ad6e5332cf53f4fc301957388df8629ec9e16db2` in
`/private/tmp/umol-216-compact-carrier`. The fully compact version uses the compact carrier for its
ordered initial partition, every recursive refinement, and backend automorphism work. Ordered color
refinement depends on carrier adjacency, so this version selects a different exact representative
for imported CHEBI:2453. Its returned correspondence remains valid and equivalent inputs converge
to the same compact representative. The difference is permitted: the existing `Canonicalize`
contract and development guides already state that canonical keys and representatives are
deterministic within a fixed release and context but are not stable across 0.x releases or suitable
as unversioned persistent identifiers. The initial rejection of this candidate imposed a stronger
compatibility requirement than the repository's actual contract.

A representative-preserving hybrid was measured while isolating the cause: it retained the current
incidence graph for ordered refinement and projected its cell ids onto a compact backend adapter.
That path produced only modest setup and backend gains because refinement still operated on the
complete carrier. It is not retained as the primary candidate now that representative stability is
known not to be required.

The retained test-only prototype keeps the public `IncidenceGraph` as the complete semantic source
but performs topology search on the compact carrier. Search individualizes atom nodes and marked
bonds remain carrier nodes. The unchanged typed leaf key orders every bond from its mapped unordered
atom endpoints and exact normalized bond fields, and the unchanged correspondence builder derives
every bond image from those endpoints. Both policies and all prototype entry points remain
module-private and test-only; production canonicalization still uses the current carrier.

The final release-mode carrier and search measurements were:

| Case | Policy | Same current representative | Carrier nodes / edges | Direct bonds | Initial cells | Refinements | Backend calls | Leaves | Orbit-pruned |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Ordinary naphthalene | Current | Yes | 21 / 22 | — | 7 | 3 | 2 | 1 | 2 |
| Ordinary naphthalene | Ordinary single | Yes | 21 / 22 | 0 | 7 | 3 | 2 | 1 | 2 |
| Ordinary naphthalene | Modal | Yes | 10 / 11 | 11 | 3 | 3 | 2 | 1 | 2 |
| CHEBI:2453 | Current | Yes | 33 / 34 | — | 33 | 1 | 0 | 1 | 0 |
| CHEBI:2453 | Ordinary single | No | 20 / 21 | 13 | 20 | 1 | 0 | 1 | 0 |
| CHEBI:2453 | Modal | No | 20 / 21 | 13 | 20 | 1 | 0 | 1 | 0 |
| CHEBI:46245 | Current | Yes | 126 / 126 | — | 124 | 2 | 1 | 1 | 1 |
| CHEBI:46245 | Ordinary single | Yes | 77 / 77 | 49 | 76 | 2 | 1 | 1 | 1 |
| CHEBI:46245 | Modal | Yes | 77 / 77 | 49 | 76 | 2 | 1 | 1 | 1 |

The constructed naphthalene uses an incomplete order-one form, so it deliberately does not match
the complete ordinary-single key. Ground imported molecules do. Unlike the representative-preserving
hybrid, the fully compact carrier also reduces the initial partition domain and every refinement
domain.

Gross Rust allocation counts and direct phase probes show where that reduction is realized:

| Case | Policy | Setup | Search | Complete | Allocations / bytes |
| --- | --- | ---: | ---: | ---: | ---: |
| Ordinary naphthalene | Current | 13.840 us | 44.161 us | 64.693 us | 1,421 / 231,558 |
| Ordinary naphthalene | Ordinary single | 13.145 us | 40.593 us | 60.674 us | 1,356 / 223,506 |
| Ordinary naphthalene | Modal | 11.057 us | 17.415 us | 34.630 us | 670 / 98,728 |
| CHEBI:2453 | Current | 32.222 us | 30.277 us | 71.204 us | 1,695 / 223,556 |
| CHEBI:2453 | Ordinary single | 25.095 us | 15.696 us | 50.040 us | 1,119 / 150,516 |
| CHEBI:2453 | Modal | 27.576 us | 15.582 us | 52.417 us | 1,216 / 158,060 |
| CHEBI:46245 | Current | 131.419 us | 612.108 us | 768.296 us | 16,492 / 4,076,938 |
| CHEBI:46245 | Ordinary single | 97.582 us | 228.004 us | 353.049 us | 7,782 / 1,506,645 |
| CHEBI:46245 | Modal | 109.385 us | 214.812 us | 348.234 us | 8,155 / 1,534,061 |

The fully compact carrier materially reduces recursive refinement as well as adapter setup. The
release probe reduces complete time by 6% for the no-compaction ordinary-single naphthalene control,
30% for CHEBI:2453, and 54% for CHEBI:46245; modal reduces those cases by 46%, 26%, and 55%.
Allocations and gross bytes fall correspondingly. Modal is much stronger for incomplete forms that
cannot match the complete ordinary-single key, but its counting and selection cost remains visible
when it chooses the same class as ordinary-single. S2a owns the serialized Criterion matrix and the
combined comparison with dense setup; these direct probes establish that both fully compact policies
remain eligible.

Differential coverage includes every simple graph through four atoms with incomplete and complete
order-one forms, original and reverse dense ids, the current representative, the compact
representative, and both policies. Each is exactly invariant and idempotent under the compact
policy, and every returned correspondence transports its source to that representative. Focused
cases add mixed bond orders, electronically distinct order-one bonds, and a disconnected graph. A
temporary imported-corpus probe established the same convergence and transport properties from the
source, current representative, and compact representative for all three topology cases under both
policies; CHEBI:2453 is the sole case among those three whose representative changes. The retained
reference route passes 265 canonicalization unit tests, 15 compatibility tests, seven
molecule-canonicalization properties, all 18 imported-corpus cases, and warning-denying Clippy with
the `proptest` feature.

**Done.**

### S2 — Select and integrate the first independent gains

#### S2a — Compare the prototype matrix

**Module:** canonicalization benchmark results and doc 216; additive (green); no public API.
[dep: S1a, S1c]

Compare the unchanged baseline, dense setup alone, each compact-carrier choice alone, and dense setup
combined with each carrier. A candidate is eligible only after exact semantic comparison passes.
That comparison requires invariance, idempotence, correspondence transport, and unchanged
canonical-equivalence classes; it does not require equality with the current representative.
Timing selection requires separated Criterion confidence intervals on at least two representative
topology cases, including one larger mostly-single-bond case, without a systematic regression on
the higher-description corpus. Record allocation and search-statistic changes even when wall time is
neutral. Select ordinary-single, modal, or no compact carrier explicitly; do not retain two
production carrier policies.

**Evidence.** The matrix used the frozen S0a Criterion baseline, serialized runs, a 0.5 s warm-up,
a 1 s measurement, and 20 samples. The table reports mean change from that baseline; negative values
are improvements. The compact-only columns use the S1c implementation. The combined columns share
the dense entity-color ranks between carrier selection and adapter construction instead of ranking
the same exact bond forms again.

| Case | Baseline | Dense setup | Compact ordinary-single | Compact modal | Dense + ordinary-single | Dense + modal |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Constructed naphthalene | 0.00% | -16.49% | -6.93% | -47.27% | -17.82% | -54.60% |
| Constructed disconnected rings | 0.00% | -12.14% | -5.21% | -43.80% | -14.09% | -50.30% |
| Feature-free connected | 0.00% | -18.24% | -41.85% | -39.77% | -48.95% | -49.29% |
| Feature-free disconnected | 0.00% | -18.78% | -37.04% | -35.28% | -44.56% | -44.92% |
| Symmetry-heavy radicals | 0.00% | -17.64% | -5.53% | -8.17% | -18.15% | -22.37% |
| CHEBI:15379 | 0.00% | -27.72% | -5.56% | -15.11% | -26.70% | -31.38% |
| CHEBI:2453 | 0.00% | -35.48% | -32.73% | -28.06% | -46.70% | -47.37% |
| CHEBI:46245 | 0.00% | -10.36% | -54.86% | -53.19% | -60.49% | -60.57% |

Dense setup and compact search are complementary. Dense setup removes key-cloning and partition
construction work whether or not the carrier shrinks; the compact carrier reduces the domain of
every topology refinement and backend call. Sharing dense ranks also removes the setup disadvantage
seen in the isolated modal probe. The final combined means were:

| Case | Dense + ordinary-single | Dense + modal | Modal relative to ordinary-single |
| --- | ---: | ---: | ---: |
| Constructed naphthalene | 53.889 us | 29.899 us | -44.52% |
| Constructed disconnected rings | 77.995 us | 45.139 us | -42.13% |
| Feature-free connected | 48.480 us | 48.282 us | -0.41% |
| Feature-free disconnected | 50.439 us | 50.081 us | -0.71% |
| Symmetry-heavy radicals | 38.305 us | 36.426 us | -4.90% |
| CHEBI:15379 | 7.328 us | 6.820 us | -6.93% |
| CHEBI:2453 | 39.057 us | 38.451 us | -1.55% |
| CHEBI:46245 | 306.870 us | 305.550 us | -0.43% |

The ordinary-single and modal timing intervals are disjoint for constructed naphthalene,
disconnected rings, CHEBI:15379, and CHEBI:2453. They overlap narrowly for CHEBI:46245, where both
combined candidates are about 60.5% faster than baseline. The compact carrier is not selected at
higher description levels. Their final combined-modal measurements therefore exercise dense setup
alone: every constructed Constitution or Structure control improved by 3.81--10.61%, and every
imported Constitution or Structure control improved by 6.73--12.95%. There is no systematic
higher-description regression.

The shared ranks also remove modal's allocation penalty. Counts include the complete operation and
are gross allocations / bytes:

| Case | Current | Compact ordinary-single | Compact modal | Dense + ordinary-single | Dense + modal |
| --- | ---: | ---: | ---: | ---: | ---: |
| CHEBI:2453 | 1,695 / 223,556 | 1,119 / 150,516 | 1,216 / 158,060 | 813 / 121,584 | 808 / 121,416 |
| CHEBI:46245 | 16,492 / 4,076,938 | 7,782 / 1,506,645 | 8,155 / 1,534,061 | 6,603 / 1,398,681 | 6,598 / 1,399,257 |

Dense ranking preserves the exact initial partition and therefore does not change the S1c search
statistics. Modal reduces the constructed naphthalene carrier from 21 nodes / 22 edges to 10 / 11;
both compact policies reduce CHEBI:2453 from 33 / 34 to 20 / 21 and CHEBI:46245 from 126 / 126 to
77 / 77. The combined candidate still passes the bounded simple-graph differential through four
atoms, the mixed-order, electronic-bond-form, and disconnected cases, exact idempotence, invariance,
and correspondence transport.

**Selection.** Dense setup proceeds to S2b. Modal exact bond form is the sole selected compact
topology carrier policy: it directly represents a largest normalized bond-form class, avoids the
ordinary-single no-compaction case for incomplete inputs, is allocation-equivalent after rank
sharing, and is no slower in the measured topology matrix. Ordinary-single remains useful evidence
but is rejected as a production policy. S2d still owns the selected carrier's private name and
representation before S3a integrates it.

**Done.**

**Delivery route.** The selected general canonicalization changes are integrated on
`feature/canonicalization-performance`, an independent branch from `upstream/main`. The branch
contains no atom-mapping or reaction-network implementation. A fresh benchmark-only main baseline
is frozen before applying the selected changes. After the remaining refinement and backend gates,
the branch proceeds through a protected-main pull request and the `v0.7.0` release. Released `main`
is then merged into `feature/atom-mapping` for S5; this keeps the production patch independent while
still validating it on the CHEBI and network workloads before doc 216 closes.

#### S2b — Integrate or reject dense setup

**Module:** graph-IR canonicalization search and tests; additive (green); no public API.
[dep: S2a]

If selected, apply the dense-ranking and initial-partition implementation to the independent
canonicalization branch, preserve the exact partition order, and run its focused unit, integration,
property, and benchmark gates. Remove the test-only old implementation after the differential
evidence is retained. If it is not selected, record the result and remove the prototype worktree
without changing production code.

**Implementation.** `rank_initial_colors` now orders references to the existing initial-color keys,
assigns dense ranks, and constructs the initial `OrderedPartition` directly from those ranks. The
canonical search accepts that preconstructed partition. The cloned-key implementation remained only
through the differential gate; its unit tables and the dense-renumbering property established
identical ranks, ordered cells, and topology canonicalization before the comparator was removed.

The integration passed the focused rank and partition tables, the generated differential property,
all 15 canonicalization integration tests, and clippy with all graph-IR targets and the
`proptest` feature. Against the fresh `upstream/main` baseline, the complete-operation benchmark
improved by approximately 13--19% for feature-free and topology-heavy controls, 5--10% for ordinary
overlay and stereo controls, and was within the benchmark noise threshold for the two para-stereo
controls. Both para-stereo modes showed the same pattern and no case regressed.

**Done.**

#### S2c — Integrate or reject fused final construction

**Module:** graph-IR molecule remapping, reframing, constraints, and canonicalization; additive
(green); no public API. [dep: S1b]

Only if S1b identifies material duplicated work, implement one private fused path that constructs
the normalized, id-remapped, representative-framed aggregate and transports constraints without a
complete intermediate molecule. It must return the same operation-issued correspondence as the
existing path and preserve the normalization--reframing--canonicalization pipeline. Add focused
examples for topology and frame-sensitive aggregates and properties comparing the fused and
composed paths under dense renumbering. Otherwise record the decomposition and reject fusion.

**Decision.** Reject the canonicalization-only fused constructor. S1b found no separately allocated
complete molecule between remapping and reframing: reframing consumes the remapped aggregate and
mutates or moves its stores in place. Avoidable work remains in general remapping and reframing, but
the complete measured final path accounts for only 2.5--9.9% of the operation and a private fused
path would duplicate molecule construction and integrity logic. The narrower opportunities remain
follow-up improvements to the shared operations rather than a canonicalization-specific path.

**Done.**

#### S2d — Record the carrier decision

**Module:** doc 216 and status note; additive (green); no public API. [dep: S2a]

Record the selected carrier policy or the evidence for retaining the current carrier. If a compact
carrier is selected, settle its private representation and name before S3a. If the measured winner
would require changing the public `IncidenceGraph`, stop and prepare the required data-type contract
and plan correction before implementation.

**Decision.** The selected private representation is `CompactTopologyCarrier`. It contains the
`AutomorphismAdapter` used by search and its initial `OrderedPartition`. Atom vertices occupy the
adapter's source prefix. Bonds in the largest normalized exact bond-form class are represented as
direct edges; ties select the lowest dense semantic color. Every other bond is represented by a
colored subdivision vertex. The existing `SubdivisionNodeSource` and general canonical search are
reused, so production has no carrier-policy enum, new node-source vocabulary, or separate compact
search implementation.

The public `IncidenceGraph` remains the complete semantic representation and supplies the typed leaf
key and correspondence. The compact carrier affects only private topology search: the leaf key
still orders all bonds by mapped endpoints and normalized fields, and the returned correspondence
still maps every atom and bond.

**Done.**

### S3 — Integrate the selected carrier and optimize refinement

#### S3a — Integrate or reject the compact topology carrier

**Module:** private graph-IR canonicalization carrier/adapter code, correspondence construction,
tests, and benchmarks; additive (green); no public API. [dep: S2b, S2d]

If selected, integrate exactly the policy recorded by S2d and retain the public `IncidenceGraph` as
the complete semantic incidence representation. Search only the atom and marked-bond vertices;
induce direct-class bond images from atom images. Preserve complete bond ordering in the typed leaf
key and complete bond correspondence in the returned witness. Re-run the S1c differential suite
after integration. If no compact policy was selected, close this subitem by removing the prototype.

**Implementation.** Molecule topology canonicalization and topology canonical-key comparison now
construct `CompactTopologyCarrier` from the shared dense entity colors. The carrier uses modal-class
direct edges and subdivisions for every other bond class, then passes its initial partition to the
existing canonical search. Automorphism projection and backend branch ordering retain only the atom
source prefix. Complete bond images continue to be induced by the unchanged typed leaf candidate
and correspondence builder.

The integrated carrier passed its three-case construction table, focused mixed-order, electronic,
and disconnected relabeling cases, and bounded exhaustive validation for every incomplete and
complete order-one simple graph through four atoms. Those checks establish exact idempotence,
dense-renumbering invariance, correspondence transport, integrity, and agreement between pruned and
unpruned search. The old incidence-carrier comparator found the permitted representative change on
a five-atom disconnected case; both correspondences transported their sources correctly, and all
equivalent compact-carrier inputs converged. The comparator was then removed. All 269
canonicalization unit tests, 15 canonicalization integration tests, and seven feature-gated molecule
canonicalization properties passed.

Against the fresh `upstream/main` baseline, complete topology canonicalization improved by about
55% for constructed naphthalene, 50--51% for disconnected rings, 45--49% for the two retained
feature-free controls, and 22% for the symmetry-heavy radical control in both para-stereo modes.
Higher-description controls continued to receive the dense-setup gain without a systematic
regression; the para-stereo trichloropentane control was statistically unchanged with para stereo
enabled and every other measured control improved.

**Done.**

#### S3b — Reduce refinement allocation on the selected carrier

**Module:** `OrderedPartition` construction and refinement with module-local tests and benchmarks;
additive (green); no public API. [dep: S3a]

First measure reusable cell-index and signature storage plus sort-and-group splitting against the
selected carrier. Preserve a test-only reference refinement until the prototype proves identical
ordered cells for empty, discrete, unsplit, multiply split, disconnected, and individualized
partitions. Add a property comparing exact ordered refinement and complete canonicalization under
dense renumbering. Integrate only the narrow changes whose allocation reduction produces no
systematic benchmark regression.

**Implementation.** `OrderedPartition::refine` now reuses one node-to-cell vector, one flat
node-by-cell signature vector, and the outer next-cell vector across fixed-point rounds. Each round
computes all exact neighbor-count signatures once, orders a cell by descending signature with the
existing node-id tie break, and groups adjacent equal signatures. Unsplit cells move directly into
the next partition; split groups receive their required cell vectors. This retains the previous
asymptotic signature storage while removing the per-node signature vectors and per-cell tree maps.

The prototype retained the previous map-based implementation as a test-only reference through its
disposition gate. It produced identical ordered cells for the six prescribed partition shapes and
through recursive individualization for every simple graph on four nodes. The permanent bounded
exhaustive property compares exact recursive refinement under reverse dense renumbering and checks
complete topology canonicalization, correspondence transport, and convergence for both incomplete
and complete order-one forms. The temporary reference and allocation probe were then removed.

The allocation probe measured refinement alone with an initially unsplit partition. Constructed
naphthalene fell from 53 allocations / 4,120 gross bytes to 11 / 608, disconnected six-membered
rings from 18 / 776 to 3 / 120, and a 77-node path from 4,845 / 759,868 to 95 / 51,480. Against the
fresh integrated-S3a Criterion baseline, every complete-operation case improved in both para-stereo
modes: 8.8--25.1% without para stereo and 6.9--19.8% with para stereo. The lower bounds were the
small frame-relative constraint case; constructed naphthalene improved by about 20%, disconnected
rings by about 19--20%, the retained feature-free controls by about 12--15%, and the 77-atom
topology corpus case by about 25% without para stereo. No confidence interval admitted a
regression.

All 275 canonicalization unit tests, 15 canonicalization integration tests, and seven
feature-gated molecule canonicalization properties passed, as did warning-denying Clippy with the
`proptest` feature.

**Done.**

#### S3c — Evaluate worklist refinement

**Module:** private graph-IR partition refinement, tests, and benchmarks; additive (green); no public
API. [dep: S3b]

This subitem is conditional and deferrable. Proceed only if refinement remains a dominant measured
term after S3b. Prototype a worklist refinement that reaches exactly the same greater-signature-first
ordered equitable partition. Differentially compare the final stable ordered partition and canonical
result with the retained reference. Integrate it only if the larger change materially improves at
least the carrier-size or symmetry-scaling cases; otherwise remove it and record the result.

#### S3d — Re-measure carrier/refinement interactions

**Module:** canonicalization benchmarks and doc 216; additive (green); no public API.
[dep: S3b; S3c if executed]

Report the factorial comparison of selected setup, carrier, and refinement changes, including
carrier sizes, allocation counts and bytes, refinement rounds, backend calls, leaves, and complete
time. Do not report the sum of isolated gains as the combined improvement.

### S4 — Re-evaluate backend work against the selected search

#### S4a — Measure backend request modes and local branch order

**Module:** profiling-only graph-core/nauty and graph-IR search experiments; additive (green); no
production public API. [dep: S3d]

Measure automorphism-only backend requests for topology and constitution search and deterministic
local branch order where the backend result cannot prune. Compare exact canonical aggregates,
correspondence transport, backend calls, search nodes, and complete time. Keep any cross-crate
request seam inside the experimental worktree until its public or private contract is explicitly
settled.

#### S4b — Dispose backend session and stabilizer candidates

**Module:** doc 216 and, when warranted, a new proposed discussion document; additive (green); no
production API. [dep: S4a]

If backend allocation and repeated calls remain material, record separate estimates for reusable
native storage and generated-group stabilizers. Session ownership or a graph-core request API, and
the generated-group/BSGS work related to docs 109 and 110, require their own settled contract before
implementation. Otherwise reject or defer them with the measured ceiling. This subitem does not
introduce either mechanism under doc 216.

#### S4c — Release the selected general improvements

**Module:** repository release metadata and operations; additive (green); no public API beyond the
new package version. [dep: S3d, S4b]

Run the independent candidate's complete semantic, property, benchmark, workspace-test, lint, and
applicable Python gates. Reconcile its public-symbol inventory, prepare the 0.7.0 workspace version
and internal dependency requirements as a separate commit, push the branch to `upstream`, and merge
it through the protected-main pull-request path after CI passes. Tag the resulting `main` commit as
`v0.7.0`, dispatch the release workflow, and verify the Rust and Python publication results before
S5 consumes the release.

### S5 — Return to the reaction-network workload

#### S5a — Re-run the network matrix

**Module:** `experimental/reaction-network` measurements and doc 216; additive (green); no public
API. [dep: S2c, S4c]

Merge released `main` into `feature/atom-mapping`, then run the same normal-polarity and extended
closures, bounds, and phase accounting used above. Re-run the retained CHEBI canonicalization
corpus on that merged tree. Record exact flask and transformation counts, generation and
canonicalization time, canonicalization per derivation, and canonicalization share. Run performance
cases serially and compare against the recorded pre-plan baseline rather than an intermediate
worktree.

#### S5b — Measure the opportunity for pre-canonical duplicate screening

**Module:** reaction-network diagnostic instrumentation; additive (green); no production public API.
[dep: S5a]

Add diagnostic-only accounting for products rejected by a cheap exact structural comparison or
invariant screen before complete canonicalization. The diagnostic must still canonicalize every
candidate through the existing path and compare the answer, so it measures a possible saved-call
set without changing molecular identity or network output. Record hit rate, false-positive work,
screen cost, and attainable generation-time ceiling. Add a focused `rstest` table covering an exact
raw duplicate, a differently numbered but canonically equal product, and an unequal product, with
specific expected diagnostic counts.

#### S5c — Dispose reaction-network-specific optimization

**Module:** doc 216 and, when warranted, a new proposed discussion document; additive (green); no
production API. [dep: S5b]

If the screen has a material positive ceiling, move its exact semantics, cache lifetime, and
reaction-network integration into a separate proposed work unit. Otherwise remove the diagnostic
and record the negative result. Doc 216 does not implement a derivation-history-dependent identity
path.

### S6 — Verification, cleanup, and closeout

#### S6a — Run the final semantic and lint gates

**Module:** `umol-graph-ir` and affected workspace callers; additive (green); no public API.
[dep: S2c, S3d, S4b, S5c]

Run formatting, focused canonicalization unit and integration tests, the feature-gated
canonicalization property suite, the canonicalization benchmark target, and graph-IR clippy while
iterating. At the final gate run the applicable workspace all-target/all-feature tests and lint once;
run the Python 3.13 build and canonicalization tests if Rust changes affect the bound operation.
Reconcile the public-symbol inventory and confirm that no experimental helper, option, carrier, or
request mode escaped its planned visibility.

#### S6b — Reconcile permanent documentation and final evidence

**Module:** `docs/development/data-types.md`, `docs/development/nomenclature.md`, doc 216, and
`discussion/000-status.md`; additive (green); no public API. [dep: S6a]

Update the permanent guides only for durable representation or operational facts that changed.
Record before/after phase, allocation, Criterion, and reaction-network tables; distinguish isolated
from combined gains. Move any selected but deferred backend or reaction-network work into linked
proposed documents.

#### S6c — Remove experimental worktrees and close the work unit

**Module:** repository operations, doc 216, and the status index; additive (green); no source or
public API change. [dep: S6b]

Resolve and verify the exact temporary worktree paths and confirm that every selected general change
is present in released `main` and the merged atom-mapping branch. Remove the profiling and unselected
experimental worktrees and their temporary branches. Confirm the affected worktrees are clean
except for intended doc-216 closeout changes, set doc 216 and its index row to `Completed`, clear the
index note, and record the closeout date.

### Critical path and deferrable work

The critical path is
S0a--S0b--S0c--S1a/S1c--S2a--S2b/S2d--S3a--S3b--S3d--S4a--S4b--S4c--S5a--S6.
S1b and S2c run in parallel with the search/carrier path. Performance executions are serialized
even when code work is parallel.

S3c is deferrable when narrow refinement changes remove the measured allocation term. S4b's backend
session and stabilizer implementations and S5c's reaction-network optimization are explicitly
outside the core deliverable; this document completes once each has either been rejected or moved
to a linked proposed work unit.
