# 207 — Lightweight reaction-network generator

Status: Completed
Date: 2026-08-24
Relates: [200](200-molecular-data-substrate-2026-08-19.md),
[201](201-molecular-data-first-steps-2026-08-19.md),
[205](205-mapping-test-corpus-2026-08-20.md),
[216](216-canonicalization-performance-2026-08-30.md),
[187](187-assembly-disassembly-2026-08-05.md),
[218](218-path-constraints-2026-09-01.md)

## Purpose

This document scopes a small native reaction-network generator whose immediate purpose is to
produce larger operation-witnessed atom-mapping cases for doc 205. It should establish whether the
existing reaction DSL, reaction application, canonicalization, correspondence composition, and
corpus tooling are already sufficient for useful network generation before any of those operations
are optimized for this workload.

This is intentionally a lightweight experiment and must not grow into a free-standing polished
product. It is also the second downstream expansion from the molecular-data experiment--doc 201 led
to the atom-mapping corpus in doc 205, which now motivates reaction-network generation. That
provenance is useful, but it is also a reason to keep the new scope narrow.

The required result is an in-memory network generator and quasireaction-subgraph (QRS) extractor
that can produce larger controlled atom-mapping evidence, together with persistent QRS GraphML and
mapping artifacts that can replace the old corpus. It is not a general persistent network service,
a database workload, or the billion-node architecture considered in doc 200.

## Reference inputs

The existing Python projects are behavioral references only:

- `/Users/dr/Source/python/colibri2` provides the earlier breadth-first reaction-network generation
  behavior;
- `/Users/dr/Source/python/pheasant2` provides the quasireaction-subgraph selection behavior; and
- `/Users/dr/Projects/chemical-network/smarts-rules/rules-o.toml` and `rules-n.toml` provide example
  primitive reaction rules.

The generator must not depend on these projects, RDKit, their serialized SMILES, GraphML, or their
data infrastructure. The SMIRKS rules are small mapped two-site graph edits: atom-charge changes
combined with localized bond addition, removal, or order change. They fit native umol `Reaction`
values directly and do not justify adding general SMIRKS support.

The reference files are not normative rule data. In particular, `rules-o.toml` assigns the name
`CO_12_10` to two different rules. Native transcription must give every rule an unambiguous name
and verify the intended forward/reverse pairs.

The native rule catalog is one EDN document. EDN provides the outer catalog envelope--including
rule names and inverse links--and each rule value uses the existing `Reaction` EDN DSL directly.
This is not SMIRKS-in-TOML, a second rule notation, or a reason to add reaction serialization
outside graph IR. The experimental crate only reads the catalog envelope and delegates reaction
values to the existing graph-IR EDN decoder.

The network cases are data in a second EDN document. Each entry carries its name, formula, native
`Molecule` EDN seed, optional reference counts, and QRS parameters. Neither the rules nor the cases
are Rust constants or case-specific branches. The readers accept arbitrary documents of the same
shape, and the command-line runner takes their paths explicitly; the five checked-in entries are
reference inputs to the generic operations.

## Initial reproduction experiment

The first target is deliberately concrete: reproduce one existing oxygen-rule network at each size
from two through six non-hydrogen atoms. Each run uses a native reaction-DSL transcription of the
normal-polarity, fully reversible rule set represented by `rules-o.toml`, together with the same
starting composition and seed as the corresponding reference network.

The five reference networks are:

| Size | Formula | Seed | Nodes | Undirected edges |
| --- | --- | --- | ---: | ---: |
| 2 atoms | CH4O | `CO` | 6 | 6 |
| 3 atoms | C2H4O | `CC=O` | 47 | 69 |
| 4 atoms | C3H6O | `CCC=O` | 430 | 959 |
| 5 atoms | C4H6O | `CC=CC=O` | 4,513 | 13,584 |
| 6 atoms | C5H6O | `CCC#CC=O` | 49,637 | 194,750 |

The node and edge counts are the corresponding fields of the historical
`classification/network-stats.csv`; the seeds come from `classification/networks.csv`. In that
statistics script, edges are counted after conversion to a simple undirected graph. Its ring- and
strain-filtered neutral counts are not reproduction objectives for this experiment.

For each of the five networks, the experiment records:

- canonical flask count;
- simple directed source/target adjacency count;
- simple undirected adjacency count, directly comparable to `network-stats.csv`;
- complete native transformation count, which may exceed the simple adjacency count when several
  rules or occurrences connect the same flasks;
- neutral endpoint count under the QRS charge predicate;
- extracted QRS count and path multiplicities; and
- single-core generation and extraction time.

Count comparison is diagnostic rather than a conformance assertion. The reference implementation
canonicalizes through SMILES, discards parallel transformations in a simple directed graph, and may
contain defects in the distributed code. A difference must be assigned to native-versus-SMILES
representation, canonical identity, simple-edge collapse, QRS selection, a reference defect, or a
native-generator defect; exact equality is not assumed in advance.

Reproducing these five networks is a calibration step, not the endpoint of the experiment. Once
the count differences and reversibility checks are understood, the same generator must be usable
with additional formulas and seeds, larger bounded network closures, and QRS extraction settings
beyond the historical shortest-path-only run. The five networks, their seeds, and the historical
path-length and slack parameters therefore remain explicit inputs rather than being hardcoded into
the generator or extractor.

The useful new output is operation-witnessed evidence that the historical GraphML does not retain:
rule applications, edge correspondences, composed endpoint correspondences, and the paths
supporting them. Persisted QRS values and their non-symmetry-equivalent endpoint mapping classes are
durable corpus artifacts, not values generated only to be counted and discarded. A QRS GraphML file
preserves the graph record; faithful endpoint mappings and their evidence may be stored separately
through the doc-201 columnar substrate or another database and joined to that QRS. Subsequent atom-
mapping work should read those artifacts instead of reconstructing evidence from the historical
network archive.

The reversible rule set supplies a stronger invariant than aggregate counts. Every native rule must
have its intended inverse, and every generated directed edge must be recoverable by applying that
inverse at the target and reaching the canonical source flask. The check also composes the forward
and reverse atom correspondences. Literal identity demonstrates labeled reversal; a nontrivial
source automorphism is recorded separately from both identity and failure to recover the source.

## Settled first-cut approach

### Native rules and application

Rules are written in the umol reaction DSL and lowered to ordinary graph-IR `Reaction` values before
generation. The first cut uses `Reaction::apply` with an explicit substructure-matching
configuration. It does not introduce a second primitive-rule representation, a specialized
two-atom executor, or a lower-level graph-core rewrite path.

This deliberately exercises the operations that have just been built. Each successful application
produces a `ReactionDerivation`, including the concrete source, product, and operation-issued
correspondence. If application performance is inadequate, the first alternative to investigate is
a borrowed or streaming form of the same reaction operation. A compiled two-site executor is a
later optimization candidate, not part of the initial semantics.

A network state is one `Molecule`, which may be disconnected. Spectator components remain part of
the host and product. A disconnected rule lhs can therefore describe an intercomponent association
or an intracomponent closure through the same graph-rewrite semantics; no separate flask or
cyclization representation is introduced.

This complete molecular state is called a flask in the network vocabulary. `FlaskId` addresses its
canonical `Molecule`; the experiment does not add a second `Flask` wrapper around that value. A
`Transformation` is one concrete rule application from a source flask to a target flask, retaining
the rule identity and source-to-target correspondence.

### Network construction

The first implementation is single-threaded and bounded. It starts from a canonical seed state,
applies every supplied rule to each frontier state, canonicalizes every generated product, interns
the resulting canonical molecule, and continues breadth-first with newly discovered flasks. Bounds
on generations and network size are explicit experiment inputs so an unexpectedly large closure
does not consume the machine.

Canonicalization is the closure boundary. A raw reaction product is never looked up in the flask
table or admitted to the frontier. It is first canonicalized under the one context fixed for the
run; only that canonical value is compared with existing flasks, interned if new, and scheduled for
rule application. Consequently, isomorphic products reached with different atom ids or insertion
orders are one flask rather than separate states.

Canonical equality identifies flasks under one explicit canonicalization context. The
canonical value and any derived key are scoped to the producing umol version; no persistent stable
identifier claim is made. Discovering an existing flask records another transformation rather than
discarding the application.

For each application, the source-to-raw-product correspondence from `ReactionDerivation` is
composed with the raw-product-to-canonical-product correspondence returned by canonicalization.
The resulting source-flask-to-target-flask correspondence is retained on the directed
transformation together with the rule identity and generation provenance. Thus canonical
deduplication does not erase atom identity, and two paths to the same flask may carry different
endpoint maps.

Complete network construction remains in memory and does not require a persistent network store.
Emitted QRS values are different: the corpus producer writes their graph records as
GraphML and materializes their endpoint mappings through the doc-201 storage substrate and the
experimental atom-mapping corpus path, or through another explicit database boundary. This does not
add resumability, database coordination, or a storage schema for complete reaction networks.

### QRS extraction

QRS extraction follows the useful semantics of `pheasant2` while retaining information that its
GraphML representation discarded:

- eligible endpoint pairs are neutral states;
- a retained path has no eligible neutral internal state, so it cannot be decomposed at an earlier
  neutral endpoint;
- for endpoints at shortest-path distance `d`, retain every simple path of length at most
  `d + slack`, subject to an explicit global path-length or case-count bound; and
- the displayed QRS is the induced subgraph over the union of retained path flasks, while the
  ordered paths themselves remain available as evidence.

The historical extraction used a maximum path length of eight, zero slack, at most three rings in
either endpoint, and rejection of endpoints containing a recognized strain-inducing group. The
native experiment retains the first two parameters but deliberately drops the ring and strain
filters: they do not strengthen atom-mapping evidence, and retaining the excluded cases makes the
resulting corpus more challenging. Zero slack selects all shortest paths and no longer paths. The
batch size of 1,000 only partitioned the written GraphML files; it was not a limit on QRS cases or
paths and has no semantic counterpart in the in-memory first cut. Maximum length eight and zero
slack configure the calibration run only. Later runs may increase the maximum length, admit paths
above the shortest length through positive slack, or select new networks without changing the QRS
representation.

As in `pheasant2`, path discovery uses an undirected projection of the network. The native network
nevertheless retains the direction of every transformation. Traversing a transformation against
its generation direction uses the reversed correspondence when composing the initial-to-final
mapping.

For reference parity, the initial neutral predicate follows `pheasant2`: the sum of absolute
atom-local formal charges is zero, with its explicit charged-triple-bond exceptions. This is a
quasireaction selection predicate, not a statement that all other net-neutral or charge-separated
molecules are chemically invalid. A later experiment may replace it with another explicitly named
endpoint predicate.

The extractor composes edge correspondences along each retained path, deduplicates identical full
endpoint correspondences, and retains the provenance of every path supporting each result. It can
therefore distinguish path-stable and path-dependent endpoint maps without reconstructing rule
applications after serialization.

### Atom-mapping evidence

An extracted endpoint correspondence is exact evidence for a transformation sequence that the
native generator performed. It establishes recovery of a known network path, not chemical
mechanism and not optimality under doc 205's A0 or A1 objective. The exact mapper or an independent
formulation must still certify the objective value where that claim is needed.

The producer does not select a curated endpoint sample. For a network population and extraction
domain declared by the caller, it can emit every eligible QRS and every distinct endpoint mapping
class. Molecular size, path length, rule domain, symmetry, and mapping stability are downstream
stratification dimensions for doc 205 rather than filters imposed by this operation.

## Performance and concurrency

The initial measurement is deliberately direct: generation time, number of canonical flasks,
number of transformations, and the size and path multiplicity of extracted cases for the five
reference networks. A repeatable end-to-end benchmark records those values from the first complete
generator. Fine-grained phase or kernel benchmarks are unnecessary unless the aggregate counts and
timings fail to explain observed runtime.

A single-core runtime of ten to twenty minutes per network is already useful; up to sixty minutes
is acceptable for the first cut. Classification networks are independent jobs, so they can be
scheduled as separate single-threaded processes without an internal generation work queue. The
population, resource preflight, and resulting external concurrency belong to doc 205. If a later
workload cannot be decomposed by network, concurrency remains a separate execution-design question
rather than a new network model.

Likewise, reaction-application specialization is justified only by evidence that application or
substructure matching, rather than canonicalization or closure size, is the limiting operation.
The first cut must not optimize hypothetical bottlenecks.

## Boundaries

This work does not include:

- a general reaction-network crate or stable public network API;
- SMIRKS parsing or runtime dependence on the Python reference projects;
- persistent, distributed, resumable, or billion-node generation;
- energy calculations, kinetic models, chemical ranking, or mechanistic claims;
- a general reaction-rule library or literature-template extraction;
- polished network visualization or an interactive network product;
- replacement of `Reaction::apply` with corpus-specific graph-IR code; or
- changes to graph-core or graph-IR unless the direct experiment demonstrates a generally useful
  missing operation.

## Repository boundary

Non-published research consumers belong under the top-level `experimental/` directory. The intended
layout is:

```text
experimental/
  atom-mapping/
  reaction-network/
  reaction-templates/   # if doc 201 S8 proceeds
```

The reaction-network generator is a separate `publish = false` workspace crate at
`experimental/reaction-network`. The existing `corpora/atom-mapping` crate should move mechanically
to `experimental/atom-mapping` before or alongside creation of the new crate. Its corpus artifacts
can remain purpose-named inside that experiment; the top-level directory names the lifecycle and
authority of the crates rather than presuming that each crate is itself only a corpus.

The move also renames the package from `umol-atom-mapping-corpus` to `umol-atom-mapping`. The crate
now covers source ingestion, persistent evidence, annotation, inspection, and algorithm development;
`corpus` is too narrow for that responsibility. The `experimental/` location, rather than the
package name, communicates its current lifecycle.

The move must update workspace membership, repository instructions, checked-in commands, notebook
paths, and query-script paths together. It changes organization, not atom-mapping semantics. No
other experimental crate is renamed or generalized merely to share infrastructure.

Promotion of any experimental code into the permanent umol surface requires a separate disposition
decision based on demonstrated reuse.

The first cut exports the generated directed multigraph as GraphML for inspection and use by
ordinary graph tools. GraphML already supplies nodes, directed edges, parallel edge ids, and typed
application data, so there is no reason to define another graph container. The export contains the
native molecule EDN for each flask node and the rule name for each transformation edge. It is an
inspectable topology projection, not a faithful persistence or reload format. The classification
campaign additionally writes one durable GraphML artifact per retained QRS. Exact endpoint mappings
are separate faithful records, which may be stored through the doc-201 columnar substrate and the
atom-mapping relations or in another database and joined to the QRS artifact. The runner may write
plain `.graphml` or compressed `.graphml.bz2`; compression does not change the graph format.

## First-cut construction contracts

The experimental crate needs a small library surface because its command-line runner, integration
tests, and benchmarks are separate crate consumers. That surface is explicitly experimental; it is
not a commitment to a permanent reaction-network API.

The input carriers are open values:

- `ReactionCatalogEntry` contains a unique name, the name of its inverse, and a `Reaction`;
- `NetworkCase` contains a name, formula, concrete seed, optional reference metadata, and QRS
  parameters;
- `GenerateNetworkConfig` contains the explicit matching and canonicalization choices plus flask
  and generation limits; and
- `ExtractQuasireactionConfig` contains the maximum path length, slack, and an explicit path bound.

These values have public fields and no validating constructor. Neither configuration implements
`Default`: algorithm selection and resource limits remain visible at each call site. Catalog-wide
conditions such as name uniqueness and symmetric inverse links are checked when the catalog enters
generation. `read_reaction_catalog` reads the experiment-local outer EDN envelope and returns a
typed envelope/reaction-decoding error; it delegates every reaction value to graph IR.
`read_network_cases` does the same for the network-case envelope and delegates every seed value
to graph IR. Neither reader performs canonicalization or supplies omitted experiment parameters.

Flasks and transformations use transparent zero-based physical ids. Generated networks,
transformations, reversibility results, QRS values, paths, and endpoint correspondences are
operation-issued values with private fields, no public from-parts constructors, and read-only
accessors. Their integrity comes from `generate_network`, `check_reversibility`, and
`visit_quasireaction_subgraphs`, not from callers assembling parallel vectors.

`write_network_graphml` is the only first-cut graph-export operation. It writes to a caller-supplied
byte sink and has a typed I/O/encoding error; it introduces no public GraphML data model or reader.

Resource limits yield an explicit incomplete status together with the useful partial network; they
must never masquerade as completed closure. Invalid catalog metadata, canonicalization failures,
reaction-application failures, and QRS values whose literal charge cannot be determined are typed
operation errors. In particular, QRS extraction does not use unchecked literal extraction or
silently classify a non-literal charge as neutral or non-neutral.

## Staged implementation plan

The critical path is S0 -> S1 -> S2 -> S3 -> S4 -> S5 -> S6 -> S7. Every stage ends with a green
workspace; the large calibration runs are evidence gates, not additions to the ordinary test suite.

### S0 — Establish the experimental workspace boundary

#### S0a — Move and rename atom-mapping experiment [breaking, red -> green] [dep: none] **Done**

Implementation status: complete. The experiment now lives at `experimental/atom-mapping`, and its
Cargo package is `umol-atom-mapping`. The persisted
`umol-atom-mapping-corpus.experimental.1` schema identifier remains unchanged because it identifies
the existing corpus schema rather than the Rust package. No Rust public API changed.

Mechanically move `corpora/atom-mapping` to `experimental/atom-mapping`. Update the workspace
member, repository instructions and crate map, checked-in commands, notebook paths, query-script
paths, and path-sensitive tests in the same subitem. Rename the Cargo package from
`umol-atom-mapping-corpus` to `umol-atom-mapping`, including package-targeting commands. Preserve
the source, fixtures, notebook contents, and all atom-mapping semantics. Workspace members remain
alphabetically ordered.

Verification:

- the existing atom-mapping Rust, query, and notebook checks pass from the new path;
- a repository search finds no live `corpora/atom-mapping` path or
  `umol-atom-mapping-corpus` package reference; and
- the move introduces no content changes beyond path references.

#### S0b — Create the reaction-network experiment crate [additive, green] [dep: S0a] **Done**

Implementation status: complete. The private module skeleton is registered as the non-published
`umol-reaction-network` package at `experimental/reaction-network`. It adds no public API.

Create the `publish = false` crate at `experimental/reaction-network`, add it to the workspace in
alphabetical order, and establish the module skeleton `catalog`, `network`, `reversibility`,
`quasireaction`, `graphml`, and `reference`. Add only the existing umol crate dependencies needed by
the settled design; do not add Python, RDKit, SMIRKS, a database backend, or concurrency
dependencies.

Verification:

- workspace metadata lists the new package exactly once at the intended path; and
- the new crate compiles and its empty test target passes without placeholder public APIs.

S0 gate: complete. The experimental directory is the only home of both experiment crates, their
existing checks pass, and the workspace is green.

### S1 — Transcribe the native reference inputs

#### S1a — Add the EDN rule-catalog reader and input [additive, green] [dep: S0b] **Done**

In `catalog.rs`, add open `ReactionCatalogEntry` and generic `read_reaction_catalog`. Add the
normal-polarity oxygen rules as a checked-in EDN input under the experiment's data directory. The
outer envelope carries unique names and inverse names; each rule value is decoded through the
existing `Reaction` EDN implementation. Do not add a new DSL, SMIRKS parser, TOML schema,
experiment-specific reaction serializer, or Rust match over the reference rules.

Catalog loading reports malformed envelope data and reaction decoding separately. Generation later
checks catalog-wide uniqueness, inverse resolution, and inverse-link symmetry. The loader itself
does not assert that the catalog is chemically complete.

Tests:

- exact examples cover the envelope, every rule family, and a malformed reaction value;
- table-driven tests establish unique names, resolved inverse names, symmetric inverse links, and
  the expected catalog size; and
- forward/reverse entries are checked as inverse graph edits rather than by their names alone.

Implementation status: complete. The generic reader preserves catalog order, distinguishes EDN,
envelope, and embedded-reaction failures, and delegates each reaction value to graph IR. The
checked-in catalog contains all 22 normal-polarity oxygen rules. The source TOML used the name
`CO_12_10` for two distinct rules; the carbon(+1)-oxygen(-1) inverse is named `CO_12_00` in the
EDN catalog so that both inverse pairs have unique, symmetric names. Exact graph-side tests cover
every entry and establish that each linked pair is structurally inverse.

#### S1b — Add the EDN network-case reader and input [additive, green] [dep: S1a] **Done**

In `reference.rs`, add open `NetworkCase` and generic `read_network_cases`. Add the five cases as a
checked-in EDN input under the experiment's data directory. Each entry contains a native concrete
`Molecule` seed with explicit hydrogens, formula, optional reference node and edge
counts, maximum path length, and slack. The SMILES spellings remain optional descriptive metadata;
they are not parsed at runtime. The cases are values passed to the generic generator, not branches
inside it.

Tests:

- each seed passes graph-IR integrity checking and has only literal entity forms required by the
  experiment;
- formula, non-hydrogen atom count, charge state, and seed connectivity match the table above; and
- the checked-in manifest contains exactly one case for each size from two through six, while a
  separately constructed manifest with another case is accepted by the same reader.

Implementation status: complete. The ordered reader distinguishes EDN, envelope, and seed-decoding
failures and delegates every native seed to graph IR. The five checked-in cases use concrete
explicit-hydrogen molecules and retain the supplied formulas, reference SMILES, reference counts,
maximum path length, and slack. Reference fields are optional; the same reader accepts an unrelated
case without them. Exact tests cover every seed's composition, heavy-atom skeleton, hydrogen
placement, neutrality, connectivity, and integrity.

S1 gate: complete. The EDN catalog and all five concrete seeds load independently of the Python
reference projects, and `cargo test` for the experimental crate is green.

### S2 — Generate the canonical reaction network

#### S2a — Define the network result and configuration [additive, green] [dep: S1a] **Done**

In `network.rs`, add the public input and result surface established above:

- open `GenerateNetworkConfig` with explicit `SubstructureMatchConfig`,
  `CanonicalizeContext`, maximum flasks, and maximum generations;
- transparent zero-based `FlaskId` and `TransformationId`;
- closed, operation-issued `ReactionNetwork`, `Transformation`, and
  `GenerateNetworkReport`;
- `GenerateNetworkStatus::{Complete, FlaskLimit, GenerationLimit}`; and
- `GenerateNetworkError` for catalog, canonicalization, and application failures.

The network uses a graph-core `Graph` as its undirected topology. Its dense node and edge positions
align exactly with `FlaskId` and `TransformationId`; canonical flasks and complete directed
transformations are the corresponding external payload tables. The network exposes the graph,
simple directed and undirected adjacency counts, and id-based payload lookup. `Transformation`
exposes its source and target `FlaskId`, rule name, and full source-to-target
`MoleculeCorrespondence`. It does not expose mutable parallel storage or a public constructor.

Tests:

- compile-time and runtime surface tests exercise every public accessor and enum case;
- private white-box tests reject unresolved and asymmetric inverse catalogs at the operation
  boundary; and
- limit statuses are distinguishable from both completion and errors.

Implementation status: complete. The configuration is an open field-only carrier with explicit
matching, canonicalization, flask-limit, and generation-limit inputs. Transparent physical ids
address closed operation-issued network and transformation values. A graph-core `Graph` is the
authoritative topology: `NodeId(i)` addresses `FlaskId(i)` and `EdgeId(i)` addresses
`TransformationId(i)`, including self-loops and parallel transformations. The read-only surface
exposes that graph, canonical flasks, directed applications, rule identity, generation provenance,
full entity correspondences, id lookup, and simple directed and undirected adjacency counts. The
report keeps a complete or bounded partial network distinct from its completion status. Catalog
validation rejects duplicate names, unresolved inverses, and asymmetric inverse links;
canonicalization, application-precondition, and application-item failures retain separate error
variants.

#### S2b — Implement single-threaded breadth-first closure [additive, green] [dep: S1b, S2a] **Done**

Implement `generate_network` as breadth-first closure over one possibly disconnected `Molecule`.
Canonicalize the seed before creating `FlaskId(0)`. Canonicalize every generated product with a
correspondence before any membership lookup or frontier insertion, intern only the canonical
molecule, and compose the application correspondence with the raw-product-to-canonical
correspondence before recording the transformation. Use the canonical `Molecule` itself as the
internal hash key; canonicalization defines flask identity, while ordinary equality and hashing
only perform lookup over values already admitted through that boundary.

Apply every catalog reaction with the configuration's explicit substructure-matching algorithm.
Each successful `ReactionDerivation` becomes one `Transformation` after its product is canonicalized
and its correspondence is moved into the canonical target frame. Record distinct rule applications
or occurrences even when their source and target flasks are already known. Discovering an existing
flask does not end transformation recording.
Stop only at closure or an explicit flask/generation limit and return the corresponding status.
Accumulate transformation endpoints during generation and construct the graph-core CSR once when
the complete or bounded `ReactionNetwork` is issued; do not rebuild the CSR for each discovery.

Tests:

- exact examples cover seed canonicalization, new-flask discovery, canonical deduplication,
  self-loops, parallel transformations, and disconnected products;
- isomorphic products with different atom ids and insertion orders intern as one flask while both
  transformations remain recorded;
- correspondence composition is checked against independently induced endpoint entities;
- zero generations returns only the canonical seed with `GenerationLimit` status;
- flask and generation limits retain a valid partial network without claiming completion; and
- catalog, application-precondition, application-item, and canonicalization failures retain their
  distinct causes.

Implementation status: complete. `generate_network` validates the supplied catalog, canonicalizes
the seed, and expands newly interned flasks in breadth-first order. Every raw product is
canonicalized with its correspondence before lookup; the canonical `Molecule` is the hash key. The
application and product-canonicalization correspondences are composed before each transformation is
recorded. Existing targets retain parallel applications and self-loops, while new targets receive
dense ids and enter the frontier once. Zero and reached resource bounds return explicit partial
reports. The graph topology is constructed once from the transformation table when any report is
issued. Exact tests cover canonical seed handling, new flasks, products reached through different
atom orders, self-loops, parallel transformations, disconnected products, full endpoint
correspondences, graph/payload alignment, both limits, and the operation's error boundaries.

#### S2c — Establish generator laws and the first benchmark [additive, green] [dep: S2b] **Done**

Add property tests over small concrete reversible rule systems. They establish that flask and
transformation ids resolve in both the payload tables and aligned graph topology, canonical flask
values are unique, every transformation correspondence belongs to its endpoint frames, and dense
renumbering of the seed does not change the canonical flask or adjacency sets.
Property generators may use ordinary `#[test]` entry points with behavior-naming carve-outs, as
permitted by the testing guide.

Add a repeatable benchmark for complete generation of the two- and three-atom reference networks.
Record flask, simple-adjacency, transformation, and elapsed-time measurements together so a
faster incomplete result cannot look like an improvement.

Implementation status: complete. The feature-gated property target generates complete closures of
one- to three-atom paths under reversible charge and bond-order rules. Separate properties cover
dense graph/payload id alignment, unique canonical flasks, correspondence endpoint frames, and
invariance of canonical flask and adjacency sets under dense seed remapping. No reference network
implementation is embedded in the property suite.

`cargo bench -p umol-reaction-network --bench network` performs a complete preflight and includes
its counts in each Criterion id; every timed iteration rejects an incomplete result. The optimized
run produced:

| Case | Flasks | Directed adjacencies | Undirected adjacencies | Transformations | Criterion time |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2 atoms | 6 | 12 | 6 | 19 | 2.1186–2.1366 ms |
| 3 atoms | 49 | 144 | 72 | 222 | 23.792–23.894 ms |

These native closure counts are evidence, not assertions against the historical
SMILES-canonicalized counts. Their attribution belongs to the calibration stage.

S2 gate: complete. The generator reaches complete closure for the two- and three-atom reference
cases, all four generated network laws pass, and the initial benchmark records result size with
elapsed time.

### S3 — Verify rule-set reversibility

#### S3a — Add transformation-level reversal checking [additive, green] [dep: S2b] **Done**

In `reversibility.rs`, implement `check_reversibility`. For each recorded transformation, apply its
declared inverse rule to the canonical target, canonicalize each result, and look for recovery of
the canonical source. Compose the forward and recovered reverse correspondences and classify the
result as identity, a source automorphism, or missing. Do not infer reversibility merely from the
presence of an opposite adjacency.

The per-transformation check and aggregate report are operation-issued values with accessors and no
public constructors. Failures to apply or canonicalize the inverse remain typed operational errors;
an inverse that applies but cannot recover the source is a `Missing` outcome rather than an error.

Tests:

- exact cases distinguish identity, nontrivial source automorphism, missing inverse recovery, and
  several inverse candidates;
- a parallel opposite edge cannot substitute for the declared inverse application; and
- correspondence composition induces the classified permutation on every entity family present.

Implementation status: complete. `check_reversibility` applies only the declared inverse to each
canonical target, canonicalizes every candidate, and records identity, nonidentity source
automorphism, or missing recovery. The composed forward and reverse atom map induces the complete
source correspondence across all entity families; this avoids treating a removed and recreated
bond's intermediate id as persistent identity. Identity takes precedence when several inverse
candidates recover the source. The per-transformation values and aggregate report are closed,
read-only operation results, and catalog, application, and canonicalization failures remain typed
errors. Exact tests cover all three outcomes, multiple candidates, complete induced entity-family
maps, and the rule that an unrelated opposite network edge is not reversal evidence.

#### S3b — Check the small reference closures [additive, green] [dep: S2c, S3a] **Done**

Add integration tests that generate the complete two- and three-atom networks and require every
transformation to recover its source through the declared inverse. Report identity and automorphism
counts separately so symmetry remains visible.

Implementation status: complete. The current checked-in two-atom closure has 19 transformations:
12 reverse through identity and 7 through a nonidentity source automorphism. The three-atom closure
has 222 transformations: 144 reverse through identity and 78 through a source automorphism. Neither
closure has a missing reversal. The integration test reports the transformation id, source and
target flasks, rule, and declared inverse for every missing result.

S3 gate: complete. Every transformation in both small reference closures is reversibly witnessed,
with identity and source-automorphism counts asserted separately.

### S4 — Extract quasireaction subgraphs and endpoint maps

#### S4a — Select endpoints and bounded simple paths [additive, green] [dep: S2a] **Done**

In `quasireaction.rs`, add open `ExtractQuasireactionConfig` and the closed QRS result types. Traverse
the `ReactionNetwork`'s graph-core `Graph` directly as undirected topology; each parallel edge
retains its aligned `TransformationId` as a distinct path support. Visit each unordered eligible
neutral endpoint pair once; zero-based `FlaskId` order chooses only its presentation direction and
has no chemical meaning. Compute the shortest distance and visit all simple paths no longer than
`shortest + slack`, subject to maximum length and an explicit per-case bound on supporting
transformation paths. Paths with a neutral internal flask are ineligible. No ring or strain filter
is introduced. Implement this exact BFS/bounded-DFS selection inside the experimental crate; do not
add a generic graph-core search API before reuse is demonstrated.

The neutral predicate is exact for the first cut: all atom-local charges must be literal; otherwise
extraction returns an error identifying the flask and atom. It implements the historical sum
of absolute charges and the two explicit charged-triple-bond exceptions without changing molecule
construction or graph-IR literal APIs.

`visit_quasireaction_subgraphs` uses visitor delivery rather than eagerly collecting the whole
experiment. Early visitor termination is distinct from an exhausted search. Reaching the path bound
reports an incomplete case and does not emit it as complete evidence.

Tests:

- table-driven cases cover neutral and charged endpoints, both exceptions, a non-literal charge,
  neutral internal nodes, maximum length, positive slack, and the path bound;
- graph fixtures distinguish shortest paths from longer paths admitted by slack; and
- visitor termination stops cleanly without altering already emitted cases.

Implementation status: complete. `visit_quasireaction_subgraphs` preflights every atom-local charge,
applies the two historical charged-triple-bond exceptions, and visits each retained endpoint pair
in increasing `FlaskId` order. It computes shortest distances and enumerates paths directly over the
network's graph-core topology, excludes neutral internal flasks, and treats every aligned parallel
edge as a distinct transformation support. The open config makes maximum length, slack, and the
per-case transformation-path bound explicit. Closed subgraph, path, and oriented-traversal values
retain the induced flask set and supporting transformation sequence without exposing constructors.
A bounded case is emitted with `PathLimit` status only after another eligible support is detected;
it cannot appear complete. Non-literal atom charge is a typed preflight error, and visitor
termination remains distinct from exhausted enumeration through `ControlFlow`.

#### S4b — Compose and deduplicate endpoint correspondences [additive, green] [dep: S3a, S4a] **Done**

For every retained path, orient each traversed transformation with the path. Use its stored
correspondence in the generation direction and its reverse in the opposite direction. Compose the
ordered sequence into a full endpoint correspondence. Deduplicate equal endpoint correspondences
within a QRS while retaining every supporting ordered path and its chosen transformation sequence.

`QuasireactionSubgraph`, `QuasireactionPath`, traversal records, and
`EndpointCorrespondence` remain operation-issued. Their accessors expose endpoint ids, shortest
length, induced flask ids, ordered traversals, the full composed correspondence, and supporting path
ids; they do not permit callers to manufacture mismatched frames.

Tests:

- exact paths cover all-forward, all-reverse, and mixed-direction traversal;
- parallel transformations yielding the same endpoint map are deduplicated with both supports
  retained;
- genuinely path-dependent endpoint maps remain distinct; and
- reversing a complete path yields the reverse endpoint correspondence.

Implementation status: complete. Every retained path has a zero-based local
`QuasireactionPathId`. Its transformation correspondences are composed in traversal order, using
the stored direction for forward traversal and the reversed correspondence otherwise.
`EndpointCorrespondence` retains one distinct full `MoleculeCorrespondence` and the ids of every
supporting path; exact equality performs the first-cut deduplication. The subgraph resolves path ids
and exposes the distinct endpoint correspondences without adding public constructors for paths,
traversals, endpoint correspondences, or subgraphs. Tests cover forward, reverse, and mixed paths;
equal maps from parallel transformations; distinct path-dependent maps; and reversal of a complete
path.

#### S4c — Add QRS laws and extraction benchmarks [additive, green] [dep: S3b, S4b] **Done**

Add property tests that every emitted path is simple, begins and ends at the advertised flasks,
contains no neutral internal flask, satisfies the length window, and induces its advertised
endpoint correspondence. Add extraction benchmarks for the complete two- and three-atom networks
with `max_length = 8` and `slack = 0`, recording endpoint-pair, QRS, path, distinct-map, and support
multiplicities with elapsed time.

Implementation status: complete. The feature-gated property target checks generated complete
networks through the public extraction API. It establishes endpoint order and neutrality, simple
paths, endpoint and traversal continuity, exclusion of neutral internal flasks, the configured
length window, the induced flask union, graph-edge direction, complete support partitioning, and
fresh composition and induction of every advertised endpoint correspondence.

`cargo bench -p umol-reaction-network --bench quasireaction` generates each network once and times
only extraction with `max_length = 8`, `slack = 0`, and no effective path bound. Every emitted QRS
must have `Complete` status. The optimized run produced:

| Case | Represented endpoints | QRS pairs | Paths | Distinct maps | Supports | Criterion time |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 atoms | 0 | 0 | 0 | 0 | 0 | 395.27–399.37 ns |
| 3 atoms | 4 | 6 | 47,486 | 67 | 47,486 | 214.88–218.61 ms |

One-shot release probes established the immediate larger-network boundary. These runs used complete
network closures but explicitly bounded paths per QRS, so their extraction counts are incomplete
where reported as limited and are not Criterion measurements:

| Case | Flasks | Transformations | Generation | Path bound | QRS pairs | Complete / limited | Retained paths | Extraction |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 4 atoms | 435 | 3,121 | 0.905 s | 1,000 | 119 | 53 / 66 | 79,392 | 0.426 s |
| 5 atoms | 4,537 | 40,374 | 12.689 s | 1,000 | 2,717 | 661 / 2,056 | 2,202,664 | 14.510 s |
| 6 atoms | 49,948 | 537,524 | 166.619 s | 100 | 57,548 | 1,723 / 55,825 | 5,727,636 | 284.430 s |

The six-atom result shows that bounded extraction itself is already material at the intended scale;
it does not justify an optimization in S4, but it provides the evidence needed to interpret the S6
calibration runs.

S4 gate: small-network QRS extraction is exact, reversible path composition is covered, and the
benchmark reports both work size and time.

### S5 — Add GraphML output and the inspection runner

#### S5a — Add the GraphML topology projection [additive, green] [dep: S2b] **Done**

In `graphml.rs`, add `write_network_graphml`. Write one directed GraphML node per canonical flask
and one edge with a unique id per recorded transformation, preserving parallel transformations.
Declare typed keys for the native molecule EDN and rule name, escape their text through the XML
writer, and write to a caller-supplied `Write`. Do not add a new GraphML reader or claim that the
projection can reconstruct the stored correspondences.

Tests:

- exact output cases cover an empty transformation set, a self-loop, parallel transformations, and
  XML text requiring escaping;
- parsing the output with the existing XML library recovers the expected node, edge, source,
  target, molecule, and rule values; and
- an injected sink failure is returned as the typed graph-output error.

Implementation status: complete. `write_network_graphml` writes standard directed GraphML through
`quick-xml` to a caller-supplied byte sink. Dense graph nodes become `n...` elements with native
molecule EDN in string-valued `molecule` data, and every aligned transformation becomes a uniquely
identified `e...` edge with string-valued `rule` data. Self-loops and parallel transformations are
preserved. The public surface adds only the writer and its I/O error; it adds no GraphML value type,
reader, correspondence encoding, or reload claim. Exact and parsed-output tests cover empty,
self-loop, parallel, and escaped-text cases, and an injected sink failure retains its `io::Error`.

#### S5b — Add a direct command-line runner [additive, green] [dep: S4c, S5a] **Done**

Add a binary that accepts an explicit native EDN rule-catalog path and either a network-case
manifest path plus case name or a native molecule DSL seed. It runs generation, reversibility
checking, and optional QRS extraction. The runner supplies the ordinary application defaults
GraphAndOverlays, VF2, Vismara, Nauty, and no para-stereo refinement; each remains an explicit
command-line override. Flask/generation, path-length, slack, and path-count bounds remain explicit
inputs; manifest values are visible rather than silently overridden.
Output a stable, human-readable tabular summary to stdout. An optional output path writes
`.graphml` or `.graphml.bz2` through S5a; it does not create a persistent dataset schema.

Tests:

- table-driven argument tests cover named cases, custom seeds, explicit algorithm choices, limits,
  and invalid combinations;
- output tests require completion status, all network counts, reversibility classes, QRS
  multiplicities, and elapsed phases to have unambiguous labels; and
- file tests distinguish plain and bzip2-compressed GraphML and reject unsupported output suffixes.

Implementation status: complete. The `umol-reaction-network` binary accepts either a named case
from an explicit EDN manifest or a molecule DSL seed, together with an explicit rule catalog and
generation bounds. Its application-level defaults are GraphAndOverlays, VF2, Vismara, Nauty, and
no para-stereo refinement. Every selector remains overridable; ArcMatch additionally requires its
query-path length. QRS extraction is optional; a direct seed requires explicit path length and
slack, while a named case may inherit those two values from the manifest. The tab-separated summary
shows both manifest and effective values and labels command-line versus manifest sources, so an
override is never silent.

The same summary reports generation status, flask, directed-adjacency, undirected-adjacency, and
transformation counts; all three reversibility classes; complete and path-limited QRS, path, map,
support, and represented-endpoint counts; optional GraphML output; and input, generation,
reversibility, QRS, GraphML, and total elapsed time. A `.graphml` suffix selects plain output and
`.graphml.bz2` selects bzip2 compression through the S5a writer; other suffixes are rejected before
a file is created. The command adds no library API or persistent dataset type. Table-driven tests
cover both input modes, every conditional argument relationship, the stable summary labels, and
plain, compressed, and rejected file outputs. A direct run of the two-atom manifest case reproduced
6 flasks, 6 simple edges, 19 transformations, and no missing reversals.

Long-running commands also emit progress to stderr without changing the stable stdout summary.
The stream identifies input, generation, reversibility, QRS, and GraphML phase boundaries.
Generation reports the current breadth-first generation, expanded and queued flasks, retained
flasks and transformations, application count, source flask, and active rule at approximately
one-second work intervals. Reversibility reports checked and total transformations, the three
accumulated result classes, the active transformation, and inverse candidates considered. These
events are internal logging from the existing operations; they add no callback variant or progress
type to the public library API.

#### S5c — Exercise the complete small workflow [additive, green] [dep: S5b] **Done**

Add an end-to-end integration test for the two-atom reference case through the same library calls as
the runner. It checks complete generation, the historical node and simple-undirected-edge counts,
full reversibility, and internally consistent QRS/path/map totals. Keep the three-atom case in the
benchmark and optional local evidence run unless its test runtime is comfortably small.

Implementation status: complete. The ordinary `workflow` integration target reads the checked-in
rule catalog and case manifest through their public readers, selects the two-atom case, and runs the
same generation, reversibility, and QRS operations as the command-line runner. It requires complete
closure with 6 flasks, 12 directed adjacencies, 6 simple undirected adjacencies, and 19 retained
transformations; the flask and simple-edge counts also agree with the manifest reference values.
All 19 transformations have declared-inverse witnesses, split into 12 identities and 7 source
automorphisms. Complete extraction under the case's length and slack values emits no two-atom QRS
pair, matching the established S4 result; its represented-endpoint, path, map, and support totals
are correspondingly zero. The nonzero three-atom extraction remains in the benchmark and property
suite rather than entering the ordinary end-to-end test.

S5 gate: complete. A user can run the experiment without an environment variable, Python program,
database, or hand-written identifier and can see whether closure was complete.

### S6 — Run the calibration and extend the evidence

#### S6a — Run all five reference calibrations [additive evidence, green] [dep: S5c] **Done**

Run the five named networks in increasing size and record in this document the exact command,
completion status, canonical flask count, all three edge counts, neutral endpoint count,
reversibility classes, QRS/path/map multiplicities, phase timings, and comparison with the
historical node and undirected-edge counts. A mismatch is investigated and explained; it is not
patched over by hardcoding the historical value.

The two-atom end-to-end test remains the ordinary regression gate. The larger runs are reproducible
evidence commands rather than long-running default tests. Stop a run only at an explicit limit and
record its incomplete status.

Implementation status: complete. Each calibration used the application defaults established in
S5b and the manifest's `max-path-length = 8` and `slack = 0`. The exact command was:

```text
cargo run --release -q -p umol-reaction-network -- \
  --rule-catalog experimental/reaction-network/data/normal-polarity-oxygen-rules.edn \
  --case-manifest experimental/reaction-network/data/network-cases.edn \
  --case CASE --max-flasks 100000 --max-generations 64 \
  --extract-qrs --max-paths-per-case PATH_BOUND
```

`CASE` and `PATH_BOUND` were respectively `2-atoms`/`1000000`, `3-atoms`/`1000000`,
`4-atoms`/`1000`, `5-atoms`/`1000`, and `6-atoms`/`100`. Every network reached complete closure;
only QRS path enumeration was bounded. The two- and three-atom QRS runs were complete. The larger
QRS totals are retained evidence up to their explicit per-pair path bounds, not complete corpus
sizes. The runner now reports `qrs.neutral_endpoint_count` using the same private endpoint predicate
as QRS selection, so this input population is reproducible without a second analysis program.

| Case | Reference flasks | Native flasks | Delta | Reference edges | Native edges | Delta |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 atoms | 6 | 6 | 0 | 6 | 6 | 0 |
| 3 atoms | 47 | 49 | +2 | 69 | 72 | +3 |
| 4 atoms | 430 | 435 | +5 | 959 | 968 | +9 |
| 5 atoms | 4,513 | 4,537 | +24 | 13,584 | 13,662 | +78 |
| 6 atoms | 49,637 | 49,948 | +311 | 194,750 | 196,293 | +1,543 |

| Case | Directed adjacencies | Transformations |
| --- | ---: | ---: |
| 2 atoms | 12 | 19 |
| 3 atoms | 144 | 222 |
| 4 atoms | 1,936 | 3,121 |
| 5 atoms | 27,324 | 40,374 |
| 6 atoms | 392,586 | 537,524 |

The historical values are calibration landmarks, not conformance values. Exact agreement holds
for the two-atom case, after which the native closures have small and increasing surpluses. The old
pipeline cannot support unique attribution of those differences: it admitted products through
configured charge limits and RDKit reaction/SMILES behavior, canonicalized identity through
SMILES, and stored the closure without transactional insertion. Consequently, a missing historical
flask can represent a deliberate admission decision, RDKit behavior, or a dropped concurrent
insertion. A diagnostic application of the current `colibri2` component charge test did not
reproduce the historical totals consistently, so no such filter was added to the native generator
and no historical count was hardcoded.

| Case | Identity reversals | Source-automorphism reversals | Missing reversals | Neutral endpoints |
| --- | ---: | ---: | ---: | ---: |
| 2 atoms | 12 | 7 | 0 | 1 |
| 3 atoms | 144 | 78 | 0 | 4 |
| 4 atoms | 1,936 | 1,185 | 0 | 17 |
| 5 atoms | 27,300 | 13,074 | 0 | 89 |
| 6 atoms | 392,358 | 145,166 | 0 | 483 |

All 581,260 native transformations recover their canonical source through the declared inverse;
159,510 do so through a nonidentity source automorphism. This is a direct check of the generated
rules and correspondences and does not depend on agreement with the historical closure size.

| Case | Represented endpoints | QRS pairs | Complete / limited | Paths | Distinct maps | Supports |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 atoms | 0 | 0 | 0 / 0 | 0 | 0 | 0 |
| 3 atoms | 4 | 6 | 6 / 0 | 47,486 | 67 | 47,486 |
| 4 atoms | 17 | 119 | 53 / 66 | 79,392 | 1,085 | 79,392 |
| 5 atoms | 89 | 2,717 | 661 / 2,056 | 2,202,664 | 22,804 | 2,202,664 |
| 6 atoms | 483 | 57,548 | 1,723 / 55,825 | 5,727,636 | 204,135 | 5,727,636 |

| Case | Input | Generation | Reversibility | QRS | GraphML | Total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 atoms | 0.006075 s | 0.004268 s | 0.002768 s | 0.000001 s | 0 s | 0.013113 s |
| 3 atoms | 0.006216 s | 0.026335 s | 0.037176 s | 0.212962 s | 0 s | 0.282698 s |
| 4 atoms | 0.006440 s | 0.886970 s | 1.943624 s | 0.423196 s | 0 s | 3.260326 s |
| 5 atoms | 0.006121 s | 12.502553 s | 31.251388 s | 14.379625 s | 0 s | 58.141659 s |
| 6 atoms | 0.006356 s | 165.664390 s | 411.322894 s | 282.055558 s | 0 s | 859.134601 s |

The largest calibration completed in 14 minutes 19 seconds on one core. This establishes that
independent single-process runs are usable; it does not trigger an internal work queue. The broader
classification campaign is now gated by the separate canonicalization investigation in doc 208.

**Post-0.7.0 rerun (2026-08-31).** After completing docs 208 and 216, the same five commands were
run serially at atom-mapping commit `301bf740f31cc07c018780228c072b23ce2e3795` with
`rustc 1.96.0 (ac68faa20 2026-05-25)`. Every network again reached complete closure. Flask,
directed- and undirected-adjacency, transformation, neutral-endpoint, represented-endpoint, QRS,
complete/limited, path, support, and missing-reversal counts exactly reproduce the original tables.

The original timing table above remains the pre-optimization record. The updated phase timings are:

| Case | Input | Generation | Reversibility | QRS | GraphML | Total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 atoms | 0.006781 s | 0.002822 s | 0.000777 s | 0.000004 s | 0.000002 s | 0.010599 s |
| 3 atoms | 0.000427 s | 0.008233 s | 0.009621 s | 0.227410 s | 0.000003 s | 0.245782 s |
| 4 atoms | 0.000449 s | 0.131536 s | 0.239172 s | 0.451319 s | 0.000002 s | 0.822751 s |
| 5 atoms | 0.000417 s | 1.807144 s | 3.469268 s | 15.088502 s | 0.000002 s | 20.368901 s |
| 6 atoms | 0.000432 s | 25.726392 s | 53.702460 s | 283.153811 s | 0.000002 s | 362.630484 s |

The material generation phases are 6.44-6.92 times faster for the four- through six-atom cases,
and reversibility is 7.66-9.01 times faster. QRS enumeration is essentially unchanged: the
six-atom phase moved from 282.056 to 283.154 seconds. The complete six-atom calibration consequently
falls from 859.135 seconds (14 minutes 19 seconds) to 362.630 seconds (6 minutes 3 seconds), a
2.37-fold end-to-end speedup. The four- and five-atom totals improve by 3.96- and 2.85-fold.
QRS extraction is now the dominant remaining cost rather than network generation or reversibility.

Two witness-sensitive presentation counts changed. Identity and source-automorphism reversals are
respectively 27,324 / 13,050 for five atoms and 392,586 / 144,938 for six atoms, moving 252
reversals from the latter bucket to the former while preserving their sum and zero missing
reversals. Exact endpoint-correspondence counts are 22,667 and 203,541, decreases of 137 and 594;
all supporting path counts are unchanged. These counters compare physical entity ids or raw
`MoleculeCorrespondence` values. They can therefore change when the normalization, reframing, and
canonicalization pipeline selects different valid correspondence witnesses. They do not establish
a change in the non-symmetry-equivalent endpoint-mapping classes.

Doc 216 has now resolved the canonicalization gate. The six-atom run remains well below the
60-minute work-queue trigger considered by the original plan, so the corpus producer needs no
internal generation work queue.

#### S6b — Exercise beyond historical reproduction [additive evidence, green] [dep: S6a] **Done**

Do not begin by increasing slack over the calibration networks. Positive slack primarily multiplies
path-enumeration work, and the bipartite network topology means that slack one admits no additional
paths. Coverage and throughput for slack of at least two are a separate experiment after the value
of shortest-path evidence is established.

Instead, add selected stoichiometries whose endpoints can distinguish chemically meaningful atom
provenance. For each selected endpoint pair, retain the shortest rule paths and full endpoint
correspondences, quotient the maps by independent endpoint symmetry, and compare every class with
A0. Record separately whether a class is an A0 optimum and whether it is supported by a named
chemical interpretation. A network path establishes its generated derivation; it does not by
itself establish that the corresponding mechanism occurs experimentally.

The first case is neutral formic acid plus methanol, `C2H6O3`, with methyl formate plus water as the
selected product. `A_AC_2` is addition--elimination and `A_AC_1` is elimination--addition; despite
their mechanistic difference, both predict that the methanol oxygen becomes the ester oxygen and
the formic-acid hydroxyl oxygen becomes water. `A_AL_1` predicts the opposite oxygen provenance.
The native neutral graph-edit network is used to witness those two endpoint-map classes, not to
model acid catalysis or claim one of the three mechanisms.

The seed lives in `data/atom-mapping-cases.edn`, separate from the five historical calibration
cases. Select a small, varied continuation by size, path length, rule family, symmetry, mapping
stability, and A0 agreement. Retain rule paths and full endpoint correspondences; do not introduce
a doc-205 corpus schema in this crate.

Validation compares every selected endpoint map with fresh composition of its supporting path and
spot-checks the resulting cases through the existing atom-mapping inspection workflow.

Implementation status: complete. The esterification closure completed with 987 flasks, 5,382
directed adjacencies, 2,691 undirected adjacencies, and 8,856 transformations. Methyl formate plus
water is flask 167 in this version-scoped run and lies at shortest distance four from the seed.
Slack remained zero. Complete extraction retained 288 shortest transformation paths and four full
labeled endpoint maps. Two maps with 144 path supports have the acyl-cleavage oxygen provenance;
two maps with 144 supports have the alkyl-cleavage provenance. There is no third provenance class.

The example acyl-provenance path is `0->4->22->74->167`; the example alkyl-provenance path is
`0->3->22->74->167`. Both traverse rule families `CO_10_00`, `CO_01_10`, `OH_10_10`, and
`OH_01_00`, but the concrete matches carry the two hydroxyl oxygens through different histories.
The acyl path removes the formic-acid C-O bond before adding the methanol O-C bond, so its graph-
edit order is elimination--addition, consistent with `A_AC_1`. The alkyl path likewise removes the
methanol C-O bond before adding the formic-acid O-C bond, consistent with `A_AL_1` provenance and
order.

A separate targeted length-six search establishes an addition--elimination path rather than
inferring it from the endpoint map:

```text
0->5->29->73->22->74->167
CO_21_00 · CO_01_10 · CO_10_00 · CO_12_00 · OH_10_10 · OH_01_00
```

The first edit changes the carbonyl `C=O` to `C+(-O-)`; the second adds methanol oxygen to that
carbon while the formic-acid C-O bond remains present, yielding the tetrahedral carbon bonded to
carbonyl `O-`, methanol `O+`, and the original acid oxygen. Only the third edit removes the acid
C-O bond. This is direct operation evidence for the `A_AC_2` addition--elimination ordering within
the neutral graph-edit model. It is two steps longer than the shortest path and therefore would
first appear at slack two. The targeted existence check is not a measurement or endorsement of
exhaustive positive-slack extraction. Slack is therefore not merely a throughput or completeness
knob in this use case: zero-slack extraction omits a chemically meaningful, and for ordinary
non-activated acids and alcohols usually preferred, reaction ordering. This particular omission
does not change the endpoint atom-mapping class because the `A_AC_2` and `A_AC_1` paths compose to
the same correspondence.

Both A0 implementations agree on 24 labeled optima. Canonical mapped-reaction equality places
them in exactly the same two non-symmetry-equivalent classes, split twelve and twelve by oxygen
provenance. Every optimum matches all eleven atoms and induces two localized-bond additions, two
deletions, and no order modifications. This case therefore demonstrates chemically meaningful A0
degeneracy, not a chemically supported mapping outside the A0 minimum. The search for the latter
is not needed to complete this subitem. The length-six `A_AC_2` path composes to the same acyl-
cleavage class as the length-four `A_AC_1` path and is also an A0 optimum; the endpoint map cannot
encode their different intermediate order.

This establishes a set-valued relation at both boundaries: a reactant/product pair can admit
multiple non-symmetry-equivalent endpoint mappings, and one endpoint mapping can admit multiple
path or mechanism classes. A fully atom-tracked concrete path induces one endpoint mapping, but
neither the endpoint pair nor the mapping determines the preceding member uniquely.

The implementation added no Rust public API. The reusable neutral-endpoint predicate remains a
private implementation detail shared by extraction and the command-line summary; the only new
surface is the explicitly labeled `qrs.neutral_endpoint_count` output row and the checked-in
esterification case.

#### S6c — Exercise an extended carbon--hydrogen rule set [additive evidence, green] [dep: S6b] **Done**

Add a separate checked-in carbon--hydrogen catalog that includes normal and inverted heterolysis,
homolysis, and neutral singlet carbon-sextet states. Do not modify the S6a calibration catalog or
use the Python rule-set generator as the semantic authority. Its output may identify candidate
cases, but the native catalog is derived and checked independently.

The first cut admits only the following exact local atom states:

| Element | State | Atom DSL |
| --- | --- | --- |
| C | ordinary closed shell | `C#c0#n0#u0#s` |
| C | cation | `C#c+#n0#u0#s` |
| C | anion | `C#c-#n#u0#s` |
| C | radical | `C#c0#n0#u#s2` |
| C | neutral singlet sextet | `C#c0#n#u0#s` |
| H | covalently bound | `H#c0#n0#u0#s` |
| H | proton | `H#c+#n0#u0#s` |
| H | hydride | `H#c-#n#u0#s` |
| H | radical | `H#c0#n0#u#s2` |

All hydrogens are explicit, so the rules need not spell out the invariant `#h0`. No atom has
formal-charge magnitude above one, more than one unpaired electron, a simultaneous charge and
radical, or a charged sextet state. Triplet carbenes are outside this experiment.

Derive every forward rule from a decrement of exactly one localized C-C or C-H bond order. Allocate
the released electron pair as `2+0`, `0+2`, or `1+1`; retain a transition only when both resulting
atom states occur in the table. The three allocations yield the two heterolytic orientations and
homolysis without treating them as unrelated hand-written families. In particular, allocation to
an already charged or radical carbon can produce the allowed singlet-sextet state, which is why
the three groups cannot be cleanly separated. Pair every retained rule with the structural inverse
that increments the same bond order.

Every LHS atom and bond states its charge, lone-pair count, unpaired-electron count, and
multiplicity explicitly. Every rule is checked independently for total charge conservation,
localized electron accounting, membership of all four endpoint atom states in the table, and an
exact inverse graph edit. These checks, rather than agreement with the Python generator, define the
catalog.

The expanded catalog increases alternative transformations and path multiplicity, but it does not
by itself break network bipartiteness: every rule changes the total localized bond order by one, so
bond-order parity still colors the graph and slack one still admits no paths. A future
bond-parity-preserving electronic-state rule would be a separate semantic extension, not something
introduced merely to obtain odd cycles.

Run the existing generator on the two- and three-carbon alkane cases, ethane and propane. Ethane is
the smallest case with both C-C and C-H choices; propane adds two symmetry-related C-C bonds without
introducing oxygen or nitrogen. Compare the normal-polarity and extended catalogs on the same seeds.
Record closure status, flasks, transformations, simple adjacencies, parallel transformation
multiplicity, electronic-state composition, reversibility classes, and runtime. Select concrete
endpoint pairs only after inspecting that network; do not perform an indiscriminate positive-slack
QRS run or change the existing endpoint predicate in this subitem. The result should show how much
branching and mapping multiplicity the enlarged electronic-state domain creates before oxygen or
nitrogen is added.

The first native catalog contains 54 named rules in
`experimental/reaction-network/data/extended-carbon-hydrogen-rules.edn`; the ethane and propane
seeds are checked in separately in `extended-carbon-hydrogen-cases.edn`. Parsing and reciprocal
inverse-link validation pass for the full catalog. The current rules remain bipartite because every
application changes total localized bond order by one. A singlet-to-triplet carbene state change
would preserve bond-order parity and could break that property, but triplet carbenes and such a
state-only rule are deliberately absent from this first catalog.

Initial release runs produced:

- Ethane, `C2H6`: complete closure with 101 flasks and 855 transformations in 0.372 s. The
  complete 0.900 s reversibility check found 432 identity reversals, 423 source-automorphism
  reversals, and none missing. The closure contains 19 neutral flasks and therefore 171 unordered
  neutral-flask pairs.
- Propane, `C3H8`: complete closure with 1,230 flasks and 17,929 transformations in 141.3 s. The
  intentionally partial reversibility check covered 17,556 transformations: 9,239 identity,
  8,317 source automorphisms, and none missing. The closure contains 120 neutral flasks and
  therefore 7,140 unordered neutral-flask pairs.

The decreasing rate of new-flask discovery and continuing addition of transformations is ordinary
late closure densification, not by itself evidence of symmetry overhead. A second profiling run
separated the eager match search, lazy reaction construction, product canonicalization, and
interning/recording seams already present in the implementation:

| Case | Total | Matching | Application | Canonicalization | Interning | Other |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Ethane | 0.368 s | 0.011 s | 0.004 s | 0.349 s | 0.002 s | 0.002 s |
| Propane | 143.173 s | 0.164 s | 0.107 s | 142.791 s | 0.074 s | 0.037 s |

Product canonicalization accounts for 99.73% of the profiled propane closure. Matching accounts for
0.11%, so making the existing VF2 enumeration faster cannot materially improve this case. Symmetry
reduction can still be valuable if it removes application-orbit members before their products are
canonicalized. Query-only symmetry reduction must use automorphisms of the complete rule, including
its deltas, rather than automorphisms of its LHS alone. Equivalent host sites require the separate
host-automorphism action; a complete classification is an orbit of embeddings under both rule and
host automorphisms. This has not yet been measured, so no fraction of the 17,929 applications is
currently claimed to be symmetry-redundant.

The profile does not trigger an internal reaction-network work queue. Independent classification
networks can already run as separate single-threaded processes, while symmetry-adapted rule-
application enumeration is a separate algorithmic question. The canonicalization result instead
opens doc 208: full product canonicalization itself must be understood before deciding the size and
scheduling of a population campaign in doc 205. Repeating the full reversibility diagnostic is not
required to establish the network closure. At that point, S6c remained open pending the requested
network-composition comparison.

Post-doc-208 reruns and backend sampling are recorded in
[doc 216](216-canonicalization-performance-2026-08-30.md). Five-run medians preserve the exact
closures while reducing ethane generation to 0.052 s and propane generation to 1.475 s. The full
propane reversibility diagnostic now also covers all 17,929 transformations with none missing.
These results remove canonicalization performance as the immediate blocker to the network work but
did not by themselves complete S6c's separate network-composition comparison. Further
canonicalization-performance investigation remains open in doc 216.

The remaining comparison used the 22-rule normal-polarity catalog and the 54-rule extended catalog
with the same ethane and propane manifests and generation bounds. Five-run median generation times
are paired with the exact closure and full reversibility counts:

| Catalog | Case | Flasks | Directed adjacencies | Undirected adjacencies | Transformations | Identity reversals | Source-automorphism reversals | Generation |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Normal polarity | Ethane | 4 | 6 | 3 | 17 | 6 | 11 | 0.0017 s |
| Extended | Ethane | 101 | 402 | 201 | 855 | 432 | 423 | 0.052 s |
| Normal polarity | Propane | 17 | 46 | 23 | 117 | 46 | 71 | 0.012 s |
| Extended | Propane | 1,230 | 8,760 | 4,380 | 17,929 | 9,420 | 8,509 | 1.475 s |

Every closure and reversal check completed, with no missing reversal. The extended catalog increases
the ethane flask count by 25.3 times and transformation count by 50.3 times. For propane, the
corresponding increases are 72.4 and 153.2 times. Mean directed adjacencies per flask rise from 1.50
to 3.98 for ethane and from 2.71 to 7.12 for propane; mean transformations per flask rise from 4.25
to 8.47 and from 6.88 to 14.58, respectively.

Temporary GraphML projections preserve one edge per transformation, so grouping those edges by
directed source/target pair gives the parallel-transformation comparison:

| Catalog | Case | Parallel pairs | Mean transformations per directed adjacency | Maximum multiplicity |
| --- | --- | ---: | ---: | ---: |
| Normal polarity | Ethane | 4 / 6 | 2.83 | 6 |
| Extended | Ethane | 238 / 402 | 2.13 | 12 |
| Normal polarity | Propane | 33 / 46 | 2.54 | 6 |
| Extended | Propane | 4,981 / 8,760 | 2.05 | 24 |

The extended domain therefore creates many more parallel transformations and raises the observed
maximum, but the faster growth of distinct directed adjacencies lowers mean multiplicity and the
fraction of adjacency pairs that are parallel. The enlargement is primarily a state-space and
branching effect, not a uniform increase in duplicate applications per adjacency.

Counting each local atom-state occurrence across the canonical flasks gives the carbon composition
below. Parentheses give the percentage of all carbon occurrences in that network:

| Catalog | Case | Ordinary | Cation | Anion | Radical | Neutral singlet sextet |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Normal polarity | Ethane | 3 (37.5%) | 1 (12.5%) | 4 (50.0%) | 0 | 0 |
| Extended | Ethane | 36 (17.8%) | 40 (19.8%) | 40 (19.8%) | 44 (21.8%) | 42 (20.8%) |
| Normal polarity | Propane | 21 (41.2%) | 5 (9.8%) | 25 (49.0%) | 0 | 0 |
| Extended | Propane | 850 (23.0%) | 716 (19.4%) | 716 (19.4%) | 752 (20.4%) | 656 (17.8%) |

The corresponding hydrogen composition is:

| Catalog | Case | Ordinary | Proton | Hydride | Radical |
| --- | --- | ---: | ---: | ---: | ---: |
| Normal polarity | Ethane | 21 (87.5%) | 3 (12.5%) | 0 | 0 |
| Extended | Ethane | 210 (34.7%) | 108 (17.8%) | 108 (17.8%) | 180 (29.7%) |
| Normal polarity | Propane | 116 (85.3%) | 20 (14.7%) | 0 | 0 |
| Extended | Propane | 3,888 (39.5%) | 1,618 (16.4%) | 1,618 (16.4%) | 2,716 (27.6%) |

This is a historical-catalog comparison, not a strict subset experiment. In particular, the
normal-polarity catalog's carbon anion has zero lone pairs (`#n0`), while the extended catalog's
admitted carbon anion has one (`#n`). The observed expansion therefore measures both the additional
heterolytic orientation, homolysis, and sextet transitions and the extended catalog's explicit
electron-state model. Together with the independent catalog audit, these results complete S6c.

The extended catalog now has an executable audit in the reaction-network catalog integration test.
It checks all 54 rules against the nine admitted C/H atom states, independently enumerates the 27
allowed bond-decrement transitions, verifies the `2+0`, `0+2`, or `1+1` allocation of each released
electron pair together with the corresponding endpoint formal-charge changes, and requires every
named inverse to be an exact reciprocal graph edit. This completes the catalog-state,
charge-conservation, electron-allocation, and inverse checks without making the catalog reader or
network generator responsible for this reference-input-specific derivation.

S6 gate: at least one network at each feasible size has a reproducible result, the reason for every
incomplete size is known, and doc 205 has larger exact operational witnesses rather than duplicate
primitive steps. The continuation includes at least one mechanism-discriminating endpoint case
and reports A0 agreement per non-symmetry-equivalent mapping class. The extended C/H catalog has
independent state, conservation, inverse, and small-closure evidence. Positive-slack scaling is not
a prerequisite for this gate.

### S7 — Persist QRS corpus artifacts and dispose

#### S7a — Add QRS corpus output [additive, green] [dep: S6a] **Done**

The internal work-queue branch is not triggered. The post-0.7.0 six-atom calibration completes in
about six minutes on one core, and independent network jobs can be scheduled externally. Do not add
an internal generation work queue, symmetry-adapted match enumeration, or a canonicalization
workaround in this subitem. Corpus population and campaign scheduling belong to doc 205.

Every QRS emitted by a corpus job is a long-lived artifact, not an aggregate count or temporary
visitor value. Its induced graph is written once as QRS GraphML. Each included source-network
transformation's full operation-issued `MoleculeCorrespondence` is persisted once and joined to
GraphML by the source-network transformation id. For a complete QRS, the graph and its edge
correspondences allow the eligible supporting paths to be derived without path-node or path-step
rows. For a path-limited QRS, the graph is the induced envelope of the retained prefix and the
aggregate support facts remain durable, but the individual prefix paths are intentionally not
persisted and need not be recoverable across enumeration-order changes.

Endpoint atom-mapping evidence is set-valued. Extraction composes path correspondences, projects
each result to its heavy-atom pairs, collapses full correspondences that differ only in hydrogen
assignments, and quotients the remaining heavy-atom mappings by independent endpoint
automorphisms. One deterministic heavy-atom representative is persisted per observed orbit.
Distinct heavy-atom orbits remain distinct evidence. Each representative retains its minimum
supporting path length, support counts by path length, observed heavy-mapping count, and collapsed
full-correspondence count. Completion or path-bound status remains explicit; these observations are
evidence rather than an intrinsic correctness score.

The GraphML location is the source case's explicit record location in the atom-mapping dataset.
QRS metadata, heavy-atom representatives, support aggregates, and transformation-to-correspondence
joins are optional corpus-specific relations beside the existing faithful `umol-store` records.
They do not enlarge `MoleculeCorrespondenceRecord`, make GraphML a faithful molecular format, or add
reaction-network schema to `umol-store`.

`write_quasireaction_graphml` writes exactly the induced flasks and transformations with source-
network physical ids, endpoint labels, shortest length, completion status, and the transformation-
id join. The existing full-network GraphML bytes remain unchanged.
`CorpusIngest::ingest_graph_ir_quasireaction` checks the closed QRS against its supplied source
network before atomically recording the endpoint mapping input, one faithful full correspondence
per observed source-network transformation, heavy-atom orbit representatives, and path-length
support aggregates. Tests cover exact hydrogen-placement collapse and the endpoint-automorphism
quotient separately.

The four QRS relations are optional additions to the experimental atom-mapping dataset:
`quasireaction_subgraphs`, `quasireaction_transformations`,
`quasireaction_mapping_classes`, and `quasireaction_mapping_class_supports`. Existing datasets
without them remain readable. Both nested and coordinated faithful-record layouts round-trip the
new relations and their transformation correspondences.

The `build_quasireaction_corpus` command runs one named network job. It accepts explicit generation
and extraction bounds plus optional unordered endpoint restrictions, writes compressed QRS GraphML
as values are emitted, writes the atom-mapping dataset after extraction, and writes
`quasireaction-corpus.json` last as the completion marker. The manifest contains source paths and
xxHash values, algorithms, bounds, endpoint restrictions, network census, per-QRS artifact
inventory, dataset census, and elapsed time. An incomplete network closure, missing requested
endpoint, ingestion failure, or artifact error leaves no completion manifest. With no endpoint
restriction, the command persists every eligible QRS in the declared network and extraction
domain.

Release-mode preflights used the normal-polarity catalog, the native nested layout, full network
closure, and the previously established QRS bounds:

| Case and path bound | QRS | Complete / limited | Paths | Full maps | Heavy classes | Time | Artifact size |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 atoms / 100,000 | 6 | 6 / 0 | 47,486 | 67 | 6 | 0.59 s | 276 KiB |
| 4 atoms / 1,000 | 119 | 53 / 66 | 79,392 | 1,085 | 127 | 2.40 s | 928 KiB |
| 5 atoms / 1,000 | 2,717 | 661 / 2,056 | 2,202,664 | 22,667 | 2,979 | 63.05 s | 16 MiB |

The five-atom run retained 20,674 distinct source-network transformation correspondences and
peaked at 411,271,168 bytes RSS with no swapping. These runs validate the artifact shape, checked
ingestion, output transaction boundary, and individual-job I/O behavior. They are technical
preflights rather than the QRS population defined by doc 205.

#### S7b — Reconcile and close [breaking, red -> green] [dep: S7a] **Done**

Inventory the implemented public surface against the construction contracts and remove accidental
constructors, duplicated free-function surfaces, hidden visibility hierarchies, and unused helpers.
Confirm that no corpus-specific operation entered graph IR, graph core, umol-store, umol-io, or the
Python package. Record separately evidenced candidates for reusable reaction-application,
canonicalization, path-search, or correspondence kernels; promotion requires its own design
decision.

The experimental reaction-network crate remains the generic, non-published producer used by doc
205. Its public operation-issued network, transformation, reversibility, QRS, and path values are
consumed by its runner, benchmarks, property suite, and atom-mapping persistence boundary. Reducing
the crate to commands and fixtures would either duplicate those operations in binaries or force the
atom-mapping experiment to reconstruct closed state. Promotion remains unwarranted: network
generation, bounded QRS search, and heavy-mapping evidence each still have only this experimental
consumer.

The public-surface audit found no independently constructible route into closed `ReactionNetwork`,
`Transformation`, `QuasireactionSubgraph`, path, traversal, or endpoint-correspondence values. The
new reaction-network surface is limited to `write_quasireaction_graphml` and its typed error. The
atom-mapping side adds the open `GraphIrQuasireactionRecord`, four open query-row types, their
dataset fields and aggregate accessors, and the checked
`CorpusIngest::ingest_graph_ir_quasireaction` consumer. The consumer, rather than the open carrier,
checks QRS/network contextual agreement and commits atomically. No constructor, visibility change,
or test-only accessor was added to the closed producer types.

No corpus-specific operation or schema entered graph IR, graph core, `umol-store`, `umol-io`, or
the Python package. QRS search remains in the reaction-network experiment; heavy-map projection and
orbit aggregation remain in the atom-mapping experiment. No development comparator or reusable
algorithm is stranded in test support. The repository crate map already records both experimental
crates and their responsibilities. The logical-schema module documentation now points to its
experimental physical realization without claiming stable keys or schema identity.

Final verification:

- `cargo +nightly fmt --all`;
- crate tests, feature-gated property tests, and both benchmark entry points;
- activate `umol-py/.venv` and confirm that `python` is Python 3.13;
- with that environment active, `cargo test --workspace`;
- with that environment active, `cargo clippy --workspace --all-targets -- -D warnings`; and
- `git diff --check`.

All final verification gates passed. This closes the lightweight reaction-network producer and its
QRS persistence boundary. Population-scale generation, rule-domain expansion, and aromatic and
stereo corpus construction remain active proposed work in doc 205 rather than extending this
implementation plan.

## Addendum — 2026-09-01 reversibility census

The identity/source-automorphism split recorded in the initial reproduction experiment predates
the later canonicalization and frame-transport work. A release-mode rerun with the current
implementation reached the same complete closures and reproduced every recorded flask, directed
adjacency, undirected adjacency, and transformation count. It also reported no missing reversals.
The current reversibility classification is:

| Case | Identity reversals | Source-automorphism reversals | Missing reversals |
| --- | ---: | ---: | ---: |
| 2 atoms | 12 | 7 | 0 |
| 3 atoms | 144 | 78 | 0 |
| 4 atoms | 1,936 | 1,185 | 0 |
| 5 atoms | 27,324 | 13,050 | 0 |
| 6 atoms | 392,586 | 144,938 | 0 |

Thus all 581,260 transformations still recover their canonical source through the declared
inverse; 159,258 do so through a nonidentity source automorphism. This table supersedes only the
earlier reversibility classification, not the network census or QRS evidence. The five calibration
runs used `normal-polarity-carbon-hydrogen-oxygen-rules.edn`, `max-flasks = 100000`,
`max-generations = 64`, and no QRS extraction. As a control, the historical 22-entry
`normal-polarity-oxygen-rules.edn` produced the same updated 5-atom split under the current binary,
so the shift is implementation-version evidence rather than a rule-catalog difference.
