# 200 — A molecular data substrate

Status: Informational
Date: 2026-08-19
Relates: [201](201-molecular-data-first-steps-2026-08-19.md),
[205](205-mapping-test-corpus-2026-08-20.md),
[207](207-reaction-network-spike-2026-08-24.md),
[015](015-four-domains-semantic-model-2025-03-10.md),
[131](131-reaction-application-design-2026-06-24.md),
[139](139-mutability-hashability-equality-2026-07-09.md),
[180](180-umol-facade-crate-2026-08-02.md),
[190](190-enumeration-algorithm-candidates-2026-08-08.md),
[191](191-graph-internal-coordinates-2026-08-08.md)

## Purpose

This document proposes a set of extensions for molecular I/O, storage, querying, reaction
networks, geometry, ontologies, and interactive analysis. It is a statement of architectural
possibility, not a compatibility promise and not an implementation plan. umol is young enough that
the right question is not how to serialize its present Rust structs with minimal disruption. The
right question is what data substrate would make the library unusually capable at molecular and
reaction scale, then how the in-memory types should participate in it.

A durable stable-identifier design remains desirable but should not block exploration. The first
iterations can use structural equality for exact representation identity and canonical equality for
semantic deduplication. A corpus may derive a content key from the canonical form as long as that key
is explicitly scoped to the producing umol version and canonicalization context and can be rebuilt.
This is a working identity proxy, not a disguised cross-version identifier.

[201](201-molecular-data-first-steps-2026-08-19.md) separately scopes the semantic questions and
experiments needed before an implementation plan can be written.

## Thesis

umol should make the same molecular semantics available in six physical situations:

- an ordinary owned `Molecule`, `Reaction`, `ReactionSpan`, or geometric value;
- a compact binary record;
- a row projected from a memory-mapped or remotely ranged corpus;
- a columnar batch processed without constructing one Rust or Python object per record;
- a database relation participating in chemical query operators;
- an interactive notebook object exposing structure, provenance, and witnesses.

These are not six unrelated export features. They are physical realizations of a shared logical
model. The owned graph object remains the natural unit for editing and sophisticated graph
algorithms. The batch, corpus, and relational forms become peers optimized for movement, screening,
aggregation, and network-scale computation. Neither is the one true molecular model, and neither
collapses the graph, geometric, and external-format model classes distinguished by docs 015 and 180.

The opportunity is larger than fast persistence. A native substrate can make chemical structure and
transformation first-class data-processing concepts rather than opaque blobs passed to scalar UDFs.

## Two complementary binary forms

### Portable records

`Molecule`, `Reaction`, and `ReactionSpan` need a versioned logical wire schema independent of Rust
layout and allocator details. It should support:

- bounded mechanical decoding before fields are accessed;
- selective access to headers and common fields without reconstructing the full value;
- deterministic encoding where useful for caching and reproducible artifacts;
- explicit schema and semantic versions;
- preservation of underdetermined values, constraints, relation kinds, deltas, and provenance;
- evolution without treating the current enum layout as eternal.

A row-oriented encoding is useful for messages, caches, individual database values, and random
record retrieval. FlatBuffers, Cap'n Proto, postcard, or a purpose-built encoding are implementation
candidates, not the architecture. A benchmark may show that a compact serde codec is the right
small-record form while a schema-addressable codec is the right random-access form. There is no
reason to force one encoding to win every workload.

Fast deserialization need not mean transmuting bytes into the current owned type. A borrowed record
view can answer simple questions and can materialize an owned value when an algorithm needs it. The
stored, borrowed, and owned forms must expose the same molecular semantics. Decoding rejects bytes
that cannot be safely interpreted, but chemistry, coherence, normalization, and other semantic
properties remain lazy and operation-specific just as they are for the owned container. Storage
must not introduce a fortified parallel molecule API.

### Native columnar batches

At one million to one billion records, per-object allocation and per-row decoding become the wrong
default. umol should investigate first-class `MoleculeBatch`, `ReactionBatch`, and
`ReactionSpanBatch` representations built from offsets, dense entity columns, validity bitmaps,
dictionaries, and relation tables. Arrow compatibility is highly desirable, but Arrow is the memory
interchange mechanism rather than the chemical design.

A molecule batch could hold, schematically:

- molecule offsets into contiguous atom, bond, and higher-relation arrays;
- entity-local identifiers and endpoint columns;
- compact discriminants plus value columns for literals, sets, ranges, variables, and constraints;
- optional canonical identifiers and cached feature summaries;
- separate child arrays for variable-length relation members and provenance.

Reaction batches should avoid repeating whole graph payloads when the domain model already provides
better factorization. Depending on the workload, a reaction record can reference educt structures,
a transformation or rule, products, participants, conditions, and evidence. `ReactionSpan` has an
especially natural analytical layout: one aligned entity frame with left/right membership and value
columns. It can expose bond changes, atom changes, preserved context, and reaction centers without
diffing two separately decoded molecules for every query.

Batch forms should support batch-native operations and cheap row projections. Constructing an owned
object is an explicit materialization boundary, not the hidden inner loop of every scan.

## A chemical content store

The on-disk target should be more than “a Parquet file containing serialized molecules.” A useful
corpus separates four related stores:

1. **Structure store** — molecules, patterns, fragments, canonical identities, and model class.
2. **Transformation store** — deltas, reaction rules, reaction spans, and composed rules.
3. **Record store** — observed reactions, participant roles, stoichiometry, conditions, outcomes,
   literature references, source locations, and curation state.
4. **Knowledge store** — names, external identifiers, ontology assertions, computed classifications,
   evidence, and derivation lineage.

This separation permits deduplication without erasing context. Ten observations may refer to one
chemical transformation; one rule may generate millions of derivations; one molecular identity may
have many source records, patterns, charge states, tautomers, or geometric realizations. Identity,
observation, transformation, and geometry must not be collapsed into one oversized row.

A native corpus format could use immutable segments, manifests, footer statistics, bitmap and
fingerprint indexes, and background compaction. It should permit:

- append and parallel ingestion without rewriting the corpus;
- exact lookup by stable identifier;
- projection of only the columns needed by a workload;
- remote range reads and memory mapping;
- deterministic derived artifacts with source and algorithm lineage;
- sharding and partition pruning without changing logical query semantics;
- recovery of a corpus snapshot from its manifest.

Parquet is a strong physical component for scalar and nested columns. It need not be the entire
format. Dense graph topology, large relation arrays, fingerprints, and specialized indexes may be
better stored as coordinated segments referenced by the same manifest. The important design choice
is a stable logical corpus and replaceable physical components, not allegiance to one container.

## A structured reaction template library

A reaction template library is a nearer-term use case that exercises almost the whole substrate
without first requiring billion-node network execution. It should connect five kinds of fact without
collapsing them:

- the concrete reaction reported or curated from a source;
- the generalized reaction rule or template inferred from or illustrated by it;
- stable references to substrate, product, reagent, catalyst, solvent, and other molecular records;
- structured conditions and outcomes;
- literature evidence and the provenance of every normalization, mapping, and generalization.

The central relation is many-to-many. A rule can be supported by many concrete reactions; one
concrete reaction can illustrate multiple rule abstractions at different radii or levels of
generality; a paper can report many reactions; the same molecule can occupy different roles in
different records. The library should store these relations explicitly rather than embedding an
example reaction as an annotation on a template.

A concrete reaction record also should not pretend that the source, curated interpretation, and
computed transformation are identical. It may need to preserve:

- the source reaction as reported, including partial or ambiguous structures;
- normalized molecular references and participant roles;
- asserted versus computed atom correspondence;
- the derived `Reaction` or `ReactionSpan` and the algorithm/model version that produced it;
- one or more extracted rules, their generalization policy, and their supporting witnesses;
- yield, selectivity, conversion, and other outcome measurements with units and uncertainty;
- typed conditions such as temperature, pressure, time, solvent, catalyst, atmosphere, and
  concentration, plus source text for information that has not been normalized;
- bibliographic identifiers and a precise locator for the table, scheme, example, or passage that
  supports the record.

This makes the library useful for precedent retrieval, rule inspection, condition comparison,
template extraction, reaction classification, and eventually reaction prediction. A query should be
able to move in both directions: from a concrete reaction to every rule and source assertion derived
from it, or from a proposed rule to its examples, counterexamples, conditions, and literature.

The template library naturally feeds reaction-network generation. Curated rules provide the rewrite
vocabulary; concrete examples test applicability and selectivity; condition and evidence records can
annotate generated derivations without being confused with topological possibility. Network results
in turn expose where a rule is too broad, too narrow, or missing contextual constraints.

## Reaction networks as a scale-driving workload

Docs 131 and 139 already identify reaction networks as the 100k–1B-node consumer for canonical,
deduplicated molecular values. This is not merely an eventual application of the storage work. It is
an excellent workload for shaping it.

A persistent reaction-network representation should treat the result as a derivation hypergraph:

- nodes refer to stable molecular identities;
- hyperedges record educt and product sets, rule or transformation identity, and application
  provenance;
- frontier membership and generation round are explicit operational data;
- repeated discovery of an existing product records a new derivation rather than cloning the node;
- atom correspondence and full products may be derived lazily when their witnesses make that safe;
- network snapshots can be resumed, queried, compared, and extended.

This realizes the distinction in doc 131 between a canonical-keyed node and the many deltas or paths
that reach it. It also creates a natural home for ranked routes, cycles, cut sets, motifs, and the
network algorithms surveyed in doc 190.

The bulk execution model should be able to apply a rule against candidate batches, screen impossible
educt combinations, canonicalize products, deduplicate them against an on-disk index, and append new
nodes and derivations without loading the network into a single process. A billion-node ambition
changes the abstraction before it changes the optimizer: resumable frontiers, partitionable work,
idempotent writes, and explicit provenance become basic semantics.

## Chemical relational algebra

Arrow, Parquet, DuckDB, and Polars provide excellent transport, execution, and ecosystem access. A
distinctive umol integration should go beyond registering scalar functions over binary blobs. It
should define chemical logical operators whose physical implementations can exploit molecular
indexes and batch layout:

- substructure semi-join and substructure join;
- exact or strength-parameterized molecular identity join;
- reaction-rule application over educt relations;
- reaction-center extraction and grouping;
- fragmentation and fragment containment;
- canonical grouping and deduplication;
- ontology classification with an executable witness;
- neighborhood, path, and motif queries over derivation networks.

A substructure join is a pipeline, not one predicate call per row:

1. prune segments using statistics and coarse chemical features;
2. apply bitmap and fingerprint screens;
3. run exact graph matching only on survivors;
4. optionally return the embedding or another witness, not just a Boolean.

The query optimizer should understand which stages are necessary and which indexes satisfy them. A
DuckDB extension, Polars expression namespace, DataFusion physical operator, or another engine can
host this work; the durable contribution is the logical and physical chemical operator contract.

### Querying underdetermined chemistry

umol's graph IR can represent variables, sets, ranges, and constraints. Flattening these to nulls or
rejecting them from the analytical substrate would discard one of the library's strongest ideas.
Queries over patterns and underdetermined structures need explicit semantics.

For many predicates the useful answer is three-valued:

- **definitely** — every admissible resolution satisfies the predicate;
- **possibly** — at least one admissible resolution satisfies it, but not all do;
- **cannot** — no admissible resolution satisfies it.

Indexes can expose lower and upper feature bounds to prune these cases safely. An exact operator can
return a constraint witness or matching embedding explaining the result. This gives patterns a
native analytical role rather than treating them as an exceptional string format.

## `ReactionSpan` as an analytical primitive

`ReactionSpan` is not only a visualization convenience. Its aligned before/after frame makes it a
candidate columnar fact model for transformations. Across a reaction corpus it can support:

- counts and distributions of created, removed, and modified entities;
- reaction-center fingerprints;
- retrieval by changed bond or atom environment;
- clustering and similarity of transformations;
- rule generalization from observed examples;
- comparison of authored rules with observed reactions;
- direct left/right projection when complete molecules are required.

The stored fact may still be `lhs + deltas`, an observed lhs/rhs record, or a derivation reference;
the span need not become the only authoritative reaction representation. The point is that a batch
of spans can be a primary analytical projection with its own indexes and operators.

## Geometry and 3D are peer data, not optional decoration

The data substrate should preserve doc 015's model-class distinction. Coordinates are not optional
fields that complete a graph molecule. A graph-model identity may be associated with zero, one, or
many geometric-model instances through an explicit bridge and atom correspondence.

A geometry store should therefore represent:

- a stable geometric or conformer identity and the structure identity it is associated with;
- coordinate arrays, units, coordinate frame, and atom correspondence;
- conformer ensembles rather than only a preferred conformer;
- generation method, optimization state, energy and model, uncertainty, and source provenance;
- periodic cell, trajectory, vibration, or transition-structure data when the model requires it.

Coordinates themselves are ideal dense columnar data: conformer offsets, point arrays, and mapping
arrays can be scanned or transferred to numerical tools without constructing Python point objects.
The bridge remains explicit because a geometric model can admit multiple graph interpretations and a
graph model can admit multiple geometries.

This becomes particularly interesting in reaction networks. A node may have a conformer and energy
ensemble; a derivation may carry a transition structure, trajectory, or computed barrier; a route
query may combine chemical steps with geometric or energetic evidence. The network remains a
chemical derivation graph, while geometry supplies additional peer records and edge evidence rather
than silently redefining node identity.

Exports should cover both archival and visualization uses. Conventional 3D formats remain explicit
boundary representations with documented loss. Columnar numerical export should preserve the richer
internal relation between structures, geometries, ensembles, and provenance.

## An executable ontology service

Ontology integration should not stop at attaching ChEBI, RXNO, GO, or Rhea identifiers to rows. umol
can connect symbolic knowledge to executable chemical semantics. An internal concept can combine:

- external ontology identifiers and versioned source assertions;
- an intensional definition, such as a molecular pattern, reaction-span pattern, rule, or network
  predicate;
- extensional members asserted by a source or computed by umol;
- typed relations to molecular classes, substructures, transformations, pathways, and processes;
- evidence, confidence, curation status, and provenance;
- an executable witness explaining why a record was classified.

ChEBI concepts can be related to molecular and substructure patterns. RXNO concepts can be related to
reaction-span patterns and rule families. GO biological-process terms and Rhea reactions can connect
chemical transformations to pathway and biological context. These links form a coherent graph only
if asserted knowledge, computed classification, and molecular identity remain distinguishable.

The service should support a productive discrepancy loop: evaluate an ontology definition against a
corpus, inspect false positives and negatives with structural witnesses, refine either the executable
definition or the mapping, and record the evidence. RDF/OWL may be an interchange and publication
surface. Internally, typed relational tables plus executable patterns may be more effective for bulk
classification and provenance.

This is also a route to explainable chemical data. A result can carry not only “is an amide-forming
reaction” but the matched reaction center, rule, ontology path, source version, and confidence that
justify the claim.

## Jupyter as a semantic debugger

SVG export is the minimum notebook integration, not the endpoint. umol objects should expose stable
custom MIME representations and linked tabular views so a notebook can inspect their semantics:

- 2D structure and reaction drawings with selectable entities;
- interactive 3D conformers and ensembles;
- toggles for constraints, relation layers, atom correspondence, and reaction changes;
- highlighted substructure matches and query witnesses;
- ontology paths and classification evidence;
- canonicalization labels, automorphism orbits, and alternative resolutions;
- linked selection between a dataframe row, a network neighborhood, and a molecular view;
- query-plan inspection showing which screens and exact operators admitted a candidate.

The visualization payload should reference stable identities and structured scene data rather than
requiring notebooks to scrape SVG element names. Static SVG and image fallback remain important for
export, publication, and environments without the interactive frontend.

The notebook then becomes more than a renderer. It is a semantic debugger for graph transformations,
queries, ontological definitions, canonicalization, and network generation.

## Cross-cutting semantic requirements

The substrate succeeds only if the representations agree on several points:

- **Identity is explicit.** A stable identifier, a structural equality relation, a canonical key,
  an external accession, and a corpus row location are different things.
- **Validation is preserved.** Borrowed and columnar views cannot create states that the owned public
  type would reject without labeling them as unchecked boundary data.
- **Loss is named.** Conversion between graph, geometric, ontology, and external-format records is
  fallible or lossy where the models differ.
- **Provenance is queryable.** Derived values record source records, algorithm and model versions,
  parameters, and witnesses.
- **Materialization is visible.** APIs distinguish cheap projection from reconstruction of an owned
  graph or geometry.
- **Schema evolution is semantic.** Adding an enum variant or changing a validation rule is not
  disguised as a byte-layout concern.
- **Algorithm selection remains explicit.** Storage and query statistics may choose physical
  execution strategies, but chemistry and graph-algorithm choices that affect semantics remain
  operation configuration.

## What this document does not decide

This document does not select FlatBuffers over postcard, Arrow over another columnar ABI, Parquet over
another segment format, or DuckDB over Polars. It does not assign the work to a new crate, define a
public batch API, settle which reaction projection is authoritative on disk, or specify the stable
identifier assumed as a prerequisite.

It also does not schedule the work. The next task is to settle scope, logical records, query and
validation semantics, and the evidence expected from representative experiments. That is the purpose
of doc 201. A staged implementation plan would be premature until those questions have answers.
