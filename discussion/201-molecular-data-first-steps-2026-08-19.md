# 201 — Molecular data: scope and first experiments

Status: Completed
Date: 2026-08-19
Relates: [200](200-molecular-data-substrate-2026-08-19.md),
[220](220-readable-depiction-2026-09-02.md),
[205](205-mapping-test-corpus-2026-08-20.md),
[207](207-reaction-network-spike-2026-08-24.md),
[131](131-reaction-application-design-2026-06-24.md),
[139](139-mutability-hashability-equality-2026-07-09.md),
[181](181-python-boundary-ownership-2026-08-03.md),
[199](199-open-container-integrity-2026-08-18.md)

## Purpose

This document scopes an exploratory first cut into faithful molecular boundary values, columnar
storage, and visualization. The atom-mapping corpus supplied the first stored values; molecule
generation with database-backed deduplication is the storage experiment that decides whether the
substrate is useful. The work is intended to reveal the actual representation problems. It is not a
product schema, an implementation plan, or a prerequisite-clearing exercise. Row-oriented binary
codecs remain part of the broader direction in doc 200, but they are not the first storage path
explored here.

The first iteration is intentionally provisional. Its artifacts and surface API may be replaced
after the experiment, but its implementation is not exempt from the repository's semantic and
architectural rules. It may record version-specific assumptions, use a rebuildable identity proxy,
and keep metadata opaque. Its useful output is evidence: which values round-trip, which layouts are
awkward, where allocation and decoding costs occur, and which questions are real after data has
passed through the system once.

## Discovery order

The broader questions in doc 200 remain meaningful. This document neither answers nor dismisses
them; it moves them from prerequisites of the first experiment to subjects of later, evidence-backed
design work. The intended order is:

1. build the thinnest end-to-end data path using current semantics;
2. collect concrete examples of friction, loss, duplication, cost, and useful access patterns;
3. restate the larger questions in terms of those observations;
4. compare alternatives against representative cases;
5. settle semantics and public contracts;
6. only then write an implementation plan for a durable feature.

The first spike must not zero-shot identity, lifecycle, metadata, query, geometry, ontology, or
schema-evolution contracts. Failure to answer one of those questions is not a deficiency of the
spike. Discovering which version of the question matters is part of its result.

## Experimental implementation policy

The experimental status relaxes how completely the first public surface must be designed; it does
not relax semantic or architectural discipline. The first cut may be a narrow vertical slice with
provisional names and operations, followed by an explicit cleanup pass once the experiment has
shown which concepts deserve a durable API.

The following constraints remain firm:

- preserve the distinctions between graph IR, owned graph objects, external-format boundaries,
  columnar projections, and geometric models;
- use the public operations of existing data types rather than reaching into their representation
  or adding storage-specific access paths to them;
- place new functionality in modules owned by the appropriate crate instead of routing around
  established crate boundaries;
- preserve existing test organization and semantic properties; experimental comparisons and data
  fixtures must have explicit homes and must not hide independent algorithms inside test support;
- keep mechanical decoding failure separate from lazy molecular and operation-specific semantics;
- avoid introducing a second molecular model merely because a physical encoding needs a convenient
  layout; a peer semantic model requires an explicit capability that graph IR does not provide.

The first cut need not implement every eventual method, codec, batch operation, query adapter, or
visualization feature. It should nevertheless be extension-friendly along the axes already known to
vary. In particular, the generic boundary and storage design must not bake the first test corpus, one
codec, one mapper score, or the initial metadata shortcuts into closed assumptions that require a
redesign to add another boundary family, derived projection, provenance field, or physical backend.

Extension-friendly does not mean building a generic framework in advance. The experiment should
introduce only the abstractions required by the first complete data path, keep experiment-specific
facts separate from reusable molecular boundaries, and leave unsupported semantics explicit. The
cleanup pass may rename, consolidate, or remove provisional surface area, but it must not be relied
upon to repair boundary violations or semantic shortcuts in the initial implementation.

## Working decisions

### Identity is proxied, not settled

Do not design the durable stable identifier in this work. Use the equality ladder already present:

- structural equality identifies an exact stored representation in its current frame;
- canonical equality identifies the same complete graph-IR assertion across admissible remappings;
- an experiment-local content key may hash a canonical encoding for lookup and deduplication.

Any such content key is scoped to the umol version and canonicalization context that produced it. It
may change when canonicalization changes and must be cheap to rebuild. Corpus-local record ids may be
used for joins inside one experimental dataset, but they make no cross-corpus or cross-version
promise.

This lets storage and query work expose what a stable identifier would actually need to do before a
durable scheme is designed.

### Faithful boundary values are distinct types

The experiment requires a distinct owned external boundary family using the `*Record` suffix:
`MoleculeRecord`, `ReactionRecord`, `ReactionSpanRecord`, and, if needed as an independent value,
`MoleculeCorrespondenceRecord`. Here `record` names the materialized logical boundary value, not a
byte buffer, Arrow row, or Parquet layout. Each value is self-contained, belongs to one dataset
schema, and remains distinct from graph IR, `*Dsl`, TableIR, geometric values, columnar batches, and
corpus entries.

The boundary has an external representation, not independent molecular semantics. Its logical
schema may differ from the Rust structs, but the value represented is the graph-IR value. Postcard,
FlatBuffers, Arrow, or another physical encoding does not define that meaning.

Self-contained applies to each boundary family independently. Under the working name, a
`ReactionRecord` contains the complete `Reaction` value, including the molecular lhs and its deltas,
and can be decoded without resolving references to separately stored molecule records. A
store-relative representation that instead links molecule records, correspondences, and reaction
data would be a different boundary type with a separate design; it is not an alternative physical
layout of `ReactionRecord`.

### Boundary values faithfully represent graph IR

The first cut uses faithful graph-IR boundary values. Each value is an external representation of
one graph-IR value. It owns a logical schema and a representation-integrity contract independent of
the Rust structs, but it preserves every graph-IR state and entity frame.
The defining roundtrip is exact structural equality:

```text
ir -> boundary -> ir == ir
```

This is the direct analogue of `MoleculeDsl` without the DSL's default-elision purpose. It gives the
external path a simple contract, complete fidelity, and direct property laws. The schema
consequently tracks graph-IR capabilities; that is an accepted cost of this representation.

Corpus observations, adjudication, provenance, conditions, literature data, and other semantics
outside graph IR join from separate storage layers. They do not enlarge `MoleculeRecord` or
`ReactionRecord` into a second molecular model. In particular, the current graph-IR model is the
semantic authority for atom mapping. If the atom-mapping corpus exposes molecular or reaction
semantics that it cannot express, that is evidence to adjust graph IR, not a reason to introduce a
lossily convertible atom-mapping molecule type.

This decision is scoped to the faithful storage boundary. It does not collapse genuinely distinct
models into graph IR. TableIR remains a source-format model with a fallible raise into graph IR, and
geometric and graph molecules remain peer model instances connected through
`umol-geometric-graph`.

For the first experiment, a record is not an openly assembled carrier. Records and columnar batches
are produced only by inserting representation-integrity-valid graph-IR values, and persisted data
is assumed to have been written through that same path under the matching experimental schema
generation. The provisional API therefore needs no public entry constructor, mutable record fields,
or checked/asserted constructor pair. Conversion from graph IR and reconstruction back into graph
IR preserve the producer-established integrity contract and are infallible under that contract.

Mechanical Arrow, Parquet, and schema-mismatch failures remain ordinary decoding errors. Checking
arbitrary externally assembled columns, corrupted payloads, or independently constructed record
entries before publishing a record is a later integrity boundary. The first experiment must not
silently present its closed-producer assumption as a general-purpose trusted-data decoder.

Schema identity, if it needs an explicit marker at all, belongs once in dataset metadata. It is not
a record field or a per-row dispatch tag. One experimental dataset has one schema throughout; the
first cut neither reads mixed-schema datasets nor defines compatibility or migration semantics
between schema generations.

The owned records are inspectable typed values and define the meaning shared by every physical
adapter. Under the closed producer contract, they implement the repository's infallible `FromIr`
and `IntoIr` conversion vocabulary with unit context. They need not be allocated on every hot write:
an Arrow builder may encode a graph-IR value directly, and selected Arrow rows may reconstruct graph
IR directly. Those paths are optimizations of the record conversion, not alternative schemas, and
must produce the same columns and satisfy the same exact structural-roundtrip properties.

Arrow batches and Parquet datasets are separate physical values rather than aliases or backing
storage for `*Record`. A later row-oriented codec may encode the same records without changing their
meaning. The first cut does not add borrowed record views or make a record a handle into a batch.

### `Reaction` is the central owned reaction representation

Within graph IR, represent concrete reactions and reaction rules as `Reaction`: an lhs molecule
plus `Deltas`.
`lhs + deltas` expresses disparate molecular transformations directly and readably, and it is the
representation around which reaction application, composition, and the homoiconic molecule/rule
relationship were designed.

`ReactionSpan` remains a useful derived representation for aligned before/after analysis,
visualization, interchange, or a specialized columnar projection. It is not the authoritative stored
reaction in this exploration.

Making `ReactionSpan` central later would require evidence that it materially improves important
operations or storage without obscuring the transformation vocabulary or weakening `Reaction` as
the ordinary API. Easier canonicalization and compatibility with existing ecosystems are not, by
themselves, sufficient reasons.

The claim that `Reaction` cannot be independently canonicalized is not assumed here. The experiment
should record what current reaction canonicalization actually does and any failure or ambiguity it
encounters; it should not settle representation priority by assertion.

### Conversion and validation boundaries

A boundary value does not acquire different chemical meaning because it is serialized, memory
mapped, or projected from a columnar batch. Boundary conversion is representation conversion, not
model conversion. Partial or lossy conversion remains appropriate for genuinely distinct models
such as TableIR or geometric molecules, not for this storage path.

The storage boundary may reject mechanically undecodable bytes: truncated buffers, impossible tags,
or offsets that cannot be accessed safely. That is not a new validation model. Aggregate coherence,
chemistry, normalization, and operation preconditions remain lazy and are checked only by the
operation that requires them, following the open-container direction of doc 199.

The exploration must not introduce opaque molecule handles, `Validated<T>` wrappers, checked and
unchecked API families, fallible accessors for ordinary readable fields, or a required finalization
ritual before routine use. The external representation may reveal improvements to owned
construction and editing, but both should remain open and directly usable models.

### Queries start with concrete structures

Only concrete molecule and reaction-rule matching is in scope. Use the existing graph-IR matching
semantics, embedding kinds, and explicit algorithm selection. A boundary value may be converted and
materialized before matching in the first experiment. A native operation over a faithful boundary
or columnar projection must produce the same matches.

Pattern truth over underdetermined values, three-valued query results, optimizer-visible chemical
operators, and witness-carrying database plans are later work.

## Corpus ownership

Create a top-level `corpora/` area, beginning with `corpora/atom-mapping/`. This is a research and
test corpus, not a conformance suite: there is no settled external-boundary specification or mapping
objective against which conformance could yet be claimed.

The checked-in corpus definition should own source manifests and checksums, provenance, curated
annotations, small diagnostic cases, and reproducible population or sampling definitions. Large
imported sources and generated Arrow or Parquet artifacts need not be committed and should be
rebuildable into a caller-selected data location. Minimized cases that later become durable
operation regressions move into the owning crate's tests; the full corpus does not become test
support owned by one crate.

`materials/` remains the home of papers, external code, and research inputs. It is not the runtime
or generated corpus location.

The reusable store owns only faithful molecular record families, their Arrow batch realizations,
and Parquet persistence. A corpus may eventually own a separate schema for its contextual facts,
but the atom-mapping experiment currently has only provisional in-memory relations. Their physical
schema is deferred until the study concepts are settled. Annotation relations remain consumer-owned.
A reaction-template corpus may later define its own condition, bibliography, and rule relations
over the same store substrate.

Do not introduce generic mapping-result, evidence, observation, or experiment records in
`umol-store`. This first cut does not justify a reusable corpus abstraction; it keeps mapping-
specific tables in `corpora/atom-mapping/` and shares only the molecular boundary and storage
machinery.

Within the atom-mapping corpus, the chemical substrate and the mapping experiment are separate.
`MappingInput` identifies one ordered lhs/rhs pair of `MoleculeRecord` keys; each molecule fixes one
local entity-id frame. `AtomMapping` identifies that input and carries only an atom correspondence
over its frames. It is neutral with respect to how the mapping was obtained, whether it is optimal,
and whether it is chemically correct. The referenced molecules supply its carrier counts.

Different raw or standardized readings of one source case are different mapping inputs. An explicit
alignment between two such inputs contains separate lhs-to-lhs and rhs-to-rhs correspondences. A
`Reaction` may be derived from a materializable input and atom mapping for chemical operations,
comparison, or depiction. Whether to materialize that faithful value is a physical-layout choice;
it does not become a second semantic authority. `ReactionSpan` is likewise a derived view in this
experiment.

The experiment records concrete sources of atom mappings and internal algorithm executions as
facts joined to `AtomMapping`. An imported mapper result, an exact network reconstruction, and an
internal mapping run are different sources and retain different information instead of being
forced through one generic evidence record. Chemical-input provenance, such as the Rhea reaction
row, is also distinct from mapping-result provenance, such as an RXNMapper result row or a
quasireaction reconstruction witness.

`MappingRun` records an internal execution: its input, selected algorithm and complete applicable
parameters, execution outcome, timing, and emitted atom mappings. An objective is an algorithm
parameter when the algorithm has one; the first cut does not define an independent objective
catalog. Optimality, enumeration completeness, optimum multiplicity, confidence, and similar facts
belong to a concrete algorithm result only when that algorithm can report them. They are not
universal run fields and are not evidence of chemical correctness.

`EditAnalysis` is a pure derived value computed from a well-formed `AtomMapping` and its input. It
reports compatibility and localized-bond additions, deletions, and modifications for inspection and
visualization. It does not search for another mapping or make a minimum claim, and it is not a
stored corpus artifact. Equivalence grouping is also a derived comparison result; it may later be
materialized as a cache if corpus measurements justify doing so.

Chemical correctness is reported separately from the chemical substrate, algorithm results, and
derived edit analysis. Reference support, manual adjudication, ambiguity, synthesized conclusions,
and explanatory notes may assess one atom mapping or a set of mappings. A proved optimum under an
algorithm's objective does not fill or imply that assessment. The first S6d persistence cut does
not define empty assessment or synthesis relations in advance; doc 205 owns their design together
with the concrete annotation workflow that will produce and consume them.

## Experimental crate boundary

Use a provisional `umol-store` crate for the faithful storage boundary, columnar batches, Arrow
realizations, and Parquet persistence. The crate depends directly on `umol-graph-ir`; the corpus
harness or a later facade composes it with parsing, chemistry-aware operations, and visualization.
The `store` name describes the broader role without claiming that the first experiment implements a
database engine or that Arrow is the permanent physical representation.

This seam does not settle the meaning or eventual dependency direction of `umol-io`. In particular,
it is compatible with the working hypothesis that `umol-io` may later depend on `umol-graph`, but it
does not require that redesign or place storage behind `umol-io`. The first experiment must not use
crate placement to decide that separate architectural question implicitly.

## Columnar first cut

The first storage experiment is deliberately columnar-first. At small scale, a simpler row codec
would likely be adequate in absolute terms even if it won a microbenchmark. That result would not
exercise the representation expected to matter for larger scans, joins, and batch operations.
Conversely, adding a compact row codec later is straightforward once the faithful boundary exists,
whereas allowing blob storage to shape ingestion and query APIs would make the columnar design
harder to recover. The experiment therefore biases its limited evidence toward the expected
scale-driving workload.

The semantic decomposition is already known and should not be rediscovered through arbitrary
payload columns. The logical columnar model has molecular values, eight typed entity families and
their local ids, participants and references, constraints, reactions and deltas, correspondences,
and separate corpus observations, runs, provenance, and annotations.

The first physical realization should use molecule or reaction offsets into dense child arrays and
typed relation data. It must not place serialized Postcard, FlatBuffers, or other record payloads in
a primary binary column. Row-oriented codecs, point stores, and blob-plus-index baselines are
follow-on experiments rather than conveniences in the initial schema.

Compare two Arrow realizations while keeping the molecular leaf encoding identical:

- a nested aggregate schema in which one top-level molecule or reaction row owns `List` and
  `Struct` children; and
- coordinated molecule, atom, bond, entity, member, delta, and correspondence tables connected by
  record-local keys or stored ranges.

A dataset artifact contains a caller-selected nonempty subset of the four record families; it does
not imply that molecules, reactions, reaction spans, and molecule correspondences are all present.
Calling a family write operation with an empty iterator records a present but empty family, while
not calling it leaves the family absent. A stateful artifact writer owns the new root, partition and
row-group policy, per-family key allocation, and actual physical-table inventory until a consuming
`finish` writes the manifest last. Each selected family's keys begin independently at zero, and each
family may be supplied once as an iterator so persistence does not require complete in-memory
collections.

Readers remain per record type because reconstruction of one family requires no coordinated
all-family read. They first consult the manifest inventory and distinguish a family that was not
included from a declared table whose Parquet parts are missing. `ColumnarLayout::table_names()`
describes the physical tables supported by a schema realization; the artifact manifest records only
the physical tables actually present. In the nested layout, one selected record family contributes
one physical table. In the coordinated layout, it contributes its complete typed table set. This is
artifact assembly, not a new aggregate semantic value over the four record families.

Both realizations are columnar. Arrow nested values are themselves offsets plus contiguous child
arrays; the comparison is principally between an aggregate-oriented logical schema and a normalized
relational schema. Holding the leaf representation constant isolates ownership, projection,
reconstruction, and query ergonomics from separate questions about how individual graph-IR forms
are encoded.

Entity identity remains molecule-anchored exactly as in graph IR. `AtomId`, `BondId`, and the other
typed entity ids name entries in their owning molecular frame; a query-visible cross-record
reference pairs the parent record key or row with that local id. Absolute positions in Arrow child
arrays are derived physical addresses only. They may accelerate traversal and joins, but filtering,
rechunking, concatenation, or compaction must be free to change them without changing molecular
identity or correspondence semantics.

The first realization stores explicit typed local-id columns. Nested entity structs retain the id
even where list position could derive it; coordinated entity tables retain both the parent row and
the local id. This makes filtered query results, correspondence endpoints, and witnesses
self-describing without recovering their original array positions, and it keeps the leaf schema
comparable between the two layouts. Parent ranges and absolute child offsets remain available as
physical accelerators. The storage and scan cost of the redundant monotone id columns is measured;
omitting them later is an optimization rather than a change to identity semantics.

### Concrete logical schema sketch

The first experiment should use the following concrete widths and container shapes. These are
schema choices for generated artifacts, not new graph-IR semantics. `record_key` is an artifact-
local `UInt64` assigned by the writer and is not part of `MoleculeRecord`, `ReactionRecord`, or
`ReactionSpanRecord` equality. Molecular local ids, list positions, delta ordinals, constraint-node
ids, and sparse-payload ids are `UInt32`, matching graph IR's typed ids. A reference to one of those
values in a coordinated table is therefore the composite `(record_key, local_id)` or
`(record_key, payload_id)`, never the local integer alone.

Semantic identity is distinct from physical addressing. Typed graph-IR entity ids and the stored
sequence or tree relations belong to the faithfully represented value. `record_key`, `payload_id`,
constraint-node ids, corpus relation keys, Arrow offsets, and absolute child positions only connect
parts of one generated artifact. They are unique only within the scope required to resolve those
references; they make no global, content-identity, or canonical-allocation guarantee. A writer may
reassign them during rechunking, compaction, or re-encoding without changing the represented record.
Equality, hashing, and deduplication of records must therefore not depend on those addresses.

Correspondence carrier counts are `UInt64`, not local ids. This preserves the current
`Correspondence<T>` count value independently of the `UInt32` endpoint representation; converting
an artifact whose count exceeds the platform's `usize` remains a mechanical decoding failure.

Use ordinary Arrow `List`, with 32-bit offsets, for record-owned sequences and partition batches
before any child array approaches the signed 32-bit offset limit. Billion-record scale is a
dataset-partitioning concern, not a reason for a single `LargeList` array. `LargeList` remains a
physical alternative if measurement finds a real need. Ordered values additionally carry an
explicit `UInt32` ordinal when they are placed in a relation table; the nested realization gets the
same order from the list and need not duplicate the ordinal inside each child struct.

Stable schema tags use non-null `UInt8` code columns. Dictionary encoding may be applied by Arrow or
Parquet as a physical optimization, but dictionary indices are not the semantic enum codes. The
first schema also avoids Arrow `Union`: a stable tag plus typed nullable payload columns is easier to
project through Parquet and SQL engines and makes non-applicable payload nulls explicit. Null never
means graph-IR `Undetermined`.

| Semantic value | Arrow leaf shape | Faithful cases |
| --- | --- | --- |
| typed local id or position | `UInt32` | Exact wrapper value; its field name supplies the Rust id type. |
| dataset record key | `UInt64` | Writer-assigned storage identity only. |
| enum discriminator or operator | `UInt8` | Schema-defined codebook; unknown codes are decoding errors. |
| `NumForm` | `Struct<tag: UInt8, lit: Int64?, payload_id: UInt32?>` | `Undetermined` is tag-only, `Lit` uses `lit`, and sets, ranges, arithmetic trees, and predicate trees use a typed sparse payload. |
| `ElementForm` | `Struct<tag: UInt8, atomic_number: UInt8?, payload_id: UInt32?>` | Atomic numbers `1..=118` are literals; sets, complement sets, and variables use typed payloads. `*` is the `Undetermined` tag, not element zero or null. |
| `IsotopeMassForm` | `Struct<tag: UInt8, mass_number: UInt32?, payload_id: UInt32?>` | `Natural` and `Undetermined` are distinct tag-only cases; literal mass, sets, and variables remain distinct. |
| `ElectronCountsForm` | `Struct<tag: UInt8, counts: List<Int64>?>` | The literal vector is positional and may be empty; `Undetermined` is tag-only. |
| `BooleanForm` | `Struct<tag: UInt8, lit: Boolean?>` | The form tag distinguishes `Undetermined` from a literal boolean. |
| closed literal-enum forms | `Struct<tag: UInt8, lit: UInt8?>` | The form tag distinguishes `Undetermined` from a literal; the literal enum has its own schema codebook. |
| stereo configuration | `Struct<tag: UInt8, kind: UInt8?, coset: Struct<tag: UInt8, lit: UInt32?, payload_id: UInt32?>?>` | `StereoConfigurationForm::Undetermined`, kinded undetermined cosets, literal cosets, sets, and operator terms remain distinct. |
| text carried by a form or provenance value | `Utf8` | Exact UTF-8 contents; no interning is semantic. |

`payload_id` is scoped to the owning top-level record. In the nested realization that record owns
typed payload lists. In the coordinated realization every payload row carries `record_key` and
`payload_id`. This is stable under batch splitting and Parquet rewriting while costing nothing on
the dominant `Lit` and `Undetermined` paths beyond the nullable reference.

Sparse payloads are typed rather than a generic value table. The required first-cut families are:

- ordered set/domain member rows `(payload_id, ordinal, value)` for numeric, element, isotope, and
  stereo-coset sets;
- one numeric payload header for `LitSet`, `RangeFrom`, `RangeTo`, `ArithExpr`, or `PredExpr`, with
  the range bound inline and a root node id for a recursive expression;
- arithmetic and predicate node rows with a node tag, applicable literal, variable, or operator
  fields, plus ordered child-edge rows `(payload_id, parent_node_id, ordinal, child_node_id)` and
  typed membership-set rows;
- element and isotope variable rows containing the exact name and optional membership operator and
  domain; and
- stereo-term node rows, ordered child edges, and permutation-image rows, preserving operator-tree
  shape and permutation order.

No payload family accepts values belonging to another form type. A later implementation may
generate these repetitive typed schemas from the graph-IR enums, but the generated columns remain
ordinary public Arrow fields rather than an opaque serializer embedded in Arrow.

The molecule-like entity decomposition is:

| Entity family | Identity and structural columns | Attribute columns |
| --- | --- | --- |
| atom | `atom_id: UInt32` | element, isotope mass, charge, implicit hydrogens, lone pairs, unpaired-electron count and multiplicity, plus sparse atom constraints |
| localized bond | `bond_id`, `atom_0_id`, `atom_1_id: UInt32` | order, charge, unpaired-electron count and multiplicity, plus sparse bond constraints |
| dative bond | `dative_bond_id`, `acceptor_atom_id: UInt32`; donor members | order plus sparse dative constraints |
| aromatic system | `aromatic_system_id: UInt32`; atom members | positional electron counts, charge, unpaired-electron count and multiplicity, plus sparse aromatic constraints |
| multicenter bond | `multicenter_bond_id: UInt32`; atom members | positional electron counts, charge, unpaired-electron count and multiplicity, plus sparse multicenter constraints |
| noncovalent bond | `noncovalent_bond_id`, `atom_0_id`, `atom_1_id: UInt32` | interaction-kind form plus sparse noncovalent constraints |
| stereo atom | `stereo_atom_id`, `site_atom_id: UInt32`; ordered stereo ligands | stereo configuration plus sparse stereo-atom constraints |
| stereo bond | `stereo_bond_id`, `site_bond_id: UInt32`; ordered stereo ligands | stereo configuration plus sparse stereo-bond constraints |

For localized and noncovalent bonds, `atom_0_id` and `atom_1_id` are the two scalar columns of one
fixed-arity, unordered endpoint pair. The numeric qualifiers name stored positions only; they do not
assign chemical roles or direction to the endpoints.

Each stereo-ligand member is `Struct<ordinal: UInt32, atom_id: UInt32, kind: UInt8>` in a relation
table, or the same fields without a redundant ordinal in an ordered nested list. Dative donors,
aromatic members, and multicenter members similarly use typed member rows; an ordinal is retained
even where the relation's semantic participant order is canonical or unordered, because exact
structural reconstruction must not depend on a query engine's row order. Positional electron counts
remain one attribute vector and must align with the corresponding stored member order.

Entity constraints are rows or nested entries of
`Struct<constraint_tag: UInt8, key_parameter: UInt8/UInt32?, value: typed form>`. The exact value
shape is selected by the entity family and tag; it is not a generic binary or string value. Molecule-
scope constraints instead use record-local typed nodes, an explicit ordered root list, and ordered
parent-child edges. Entity and relational leaves contain their complete typed ids and payloads.

### Two physical realizations of the same frame

The nested molecule batch has one row per stored value:

| Top-level column | Arrow type |
| --- | --- |
| `record_key` | `UInt64` |
| `atoms`, `bonds`, `dative_bonds`, `aromatic_systems`, `multicenter_bonds`, `noncovalent_bonds`, `stereo_atoms`, `stereo_bonds` | one `List<Struct<...>>` per typed entity family, using the shapes above |
| `constraint_roots` | `List<UInt32>` |
| `constraint_nodes` | `List<Struct<node_id: UInt32, tag: UInt8, typed leaf payload...>>` |
| `constraint_edges` | `List<Struct<parent_node_id: UInt32, ordinal: UInt32, child_node_id: UInt32>>` |
| typed non-literal payload families | one sparse `List<Struct<...>>` per payload family |

The coordinated realization uses the same leaf structs but places them in tables:

| Table family | Primary columns | Cardinality |
| --- | --- | --- |
| `molecule_records` | `record_key: UInt64` | one row per molecule |
| each typed entity table | `record_key`, typed `local_id`, structural columns, attribute forms | zero or more rows per molecule |
| dative/aromatic/multicenter member tables | `record_key`, typed relation id, `ordinal`, `atom_id` | zero or more rows per relation |
| stereo-ligand tables | `record_key`, typed stereo id, `ordinal`, `atom_id`, `kind` | zero or more rows per stereo entity |
| each typed entity-constraint table | `record_key`, typed entity id, constraint key and typed value | sparse; zero or more rows per entity |
| molecule constraint roots/nodes/edges | `record_key` plus the node or ordinal columns above | sparse; exact ordered forest per molecule |
| each typed form-payload table | `record_key`, `payload_id`, typed payload columns | sparse; zero or more rows per record |

This is coordinated rather than entity-attribute-value storage: there are eight typed entity table
families and typed child tables, not one `entities(kind, payload)` table. Implementations may expose
common builders internally, but table names and schemas preserve the entity type at the boundary.
The molecule, reaction lhs, and reaction-span union frame instantiate this same family under
separate schema namespaces or table prefixes; they are not mixed behind an `owner_kind` tag.

### Literal-first form columns

The expected corpus distribution is dominated by `Lit`, with occasional `Undetermined`; sets,
ranges, variables, arithmetic expressions, and predicates are rare, and recursive numeric
expressions are not expected from imported external data. The first realization therefore stores
each form-valued field as an explicit discriminator, a nullable dense literal column, and an
optional record-local reference to sparse typed non-literal payload data.

`Undetermined` has its own discriminator and is never represented by Arrow null alone. Null in the
literal or payload column means only that the payload does not apply to that tagged variant.
Non-literal payloads retain sets, complement sets, variable domains, range bounds, and recursive
expression nodes without normalization or opaque serialization. Recursive node tables preserve the
exact variant, tree shape, and child order needed for structural roundtrip.

The same leaf encoding is used in nested and coordinated layouts. Imported corpus measurements are
expected to exercise the literal fast path; small generated fixtures separately exercise every
non-literal form so absence from the first corpus does not narrow the faithful-boundary contract.
A later non-literal-heavy consumer may still reveal a different distribution.

### Sparse constraint columns

Constraints use typed sparse storage rather than wide nullable structs or opaque payloads. Each
entity family has sparse constraint rows keyed by the parent record, the entity's typed local id,
and the complete constraint key, including any key parameter such as `RingScope`. No row means that
the constraint is absent. A present constraint whose value is `Undetermined` remains an explicit
tagged value and is not collapsed into absence. The encoding preserves the sorted, unique-by-key
entry contract of the graph-IR entity constraint containers.

Molecule-scope `Constraints` use typed recursive nodes. Root position and child position preserve
the exact top-level sequence and the child order of `And` and `Or`; repeated structurally equal
nodes remain repeated. `Not` retains its single child, and each leaf retains its complete typed
payload. Entity and relational references use the same parent record plus typed local-id scheme as
participants. Conversion performs no normalization, flattening, sorting, or deduplication.

The nested and coordinated Arrow realizations may place these rows differently, but they encode the
same constraint sequence and node relations. An unconstrained record therefore pays little payload
cost while every graph-IR constraint remains representable and directly inspectable.

### Reaction and reaction-span columns

`ReactionRecord` contains a complete `MoleculeRecord` lhs followed by one logical sequence of
resolved deltas. `Deltas` is flat, not recursive, and is not ordering-invariant. Every stored delta
therefore has its explicit sequence position as well as its entity family, delta variant, referenced
entity id, and complete variant payload. Physical variant- or family-specific payload tables may use
that position to recover the one sequence; they must not replace it with independently ordered
per-family collections.

The only recursive value a delta may contain is a `Constraint` carried by `ConstraintDelta`; that
payload reuses the recursive constraint representation above. This does not make `Delta` or
`Deltas` recursive. Old and new values, complete removed values, addition participants, and every
other variant payload remain explicit so conversion neither normalizes nor reconstructs missing
delta information.

`ReactionSpanRecord` has the same structural decomposition as `MoleculeRecord`, matching the
in-memory `ReactionSpan` and `ReactionSpanEntries`: one union-frame topology, the same eight typed
entity sequences and participants, and the molecule constraint sequence. Each entity's attributes
are lifted through `EntitySpan`, with an explicit `Unchanged`, `Modified`, `Added`, or `Removed` tag;
`Modified` carries complete lhs and rhs attributes. Constraint entries use their corresponding span
tags. The record is not lowered into independent lhs and rhs molecule records.

Concretely, the nested `ReactionRecord` row contains
`Struct<record_key: UInt64, lhs: Struct<molecule frame>, deltas: List<Struct<ordinal: UInt32,
family_tag: UInt8, variant_tag: UInt8, entity_id: UInt32?, typed payload...>>, typed sparse
payloads...>`. The coordinated realization has a `reaction_records(record_key)` table, a distinct
reaction-lhs instantiation of the molecular-frame tables, one ordered `reaction_deltas` table with
the discriminator and referenced id, and family/variant-specific payload tables keyed by
`(record_key, ordinal)`. Add and remove payloads contain complete attributes and participants;
modify-field payloads contain both old and new values; modify-constraint payloads preserve both
optional sides. The central ordered delta table is never replaced by nine independently ordered
family tables.

The nested `ReactionSpanRecord` row reuses the molecule-frame shape but replaces every entity's
attribute struct by
`Struct<span_tag: UInt8, unchanged_or_one_sided: attributes?, lhs: attributes?, rhs: attributes?>`.
The coordinated realization similarly reuses a distinct span-frame table namespace with these
columns. `Unchanged`, `Added`, and `Removed` use the one-sided attribute column; `Modified` uses both
`lhs` and `rhs`. Participants and sites occur once in the union frame. Constraint-span roots use
their actual three-way `Unchanged`, `Added`, or `Removed` tag and one complete constraint tree;
`ConstraintSpan` has no `Modified` variant.

### Correspondence columns

`MoleculeCorrespondenceRecord` contains the eight typed component correspondences already present
in `MoleculeCorrespondence`; the atom correspondence is the first component and is not stored a
second time. Each component faithfully stores its left and right carrier counts plus its matched
local-id pairs in left-id order. The counts are part of `Correspondence<T>` rather than a
denormalized summary: paired ids alone do not identify the unmatched members or even the sizes of
the two carriers.

An external source that supplies only atom pairs is an atom-mapping corpus fact. It may use the same
atom-component column shape without introducing a second atom correspondence beside a full
`MoleculeCorrespondenceRecord` or requiring a general-purpose `AtomCorrespondenceRecord` in the
first cut. Links identifying the two corpus molecular frames related by either value belong to the
corpus layer; they are not fields of the self-contained correspondence value.

The nested correspondence row has eight named component structs, one for each typed entity family:

`Struct<left_count: UInt64, right_count: UInt64,
pairs: List<Struct<left_id: UInt32, right_id: UInt32>>>`.

The coordinated realization uses eight typed component tables
`<family>_correspondences(record_key, left_count, right_count)` and eight typed pair tables
`<family>_correspondence_pairs(record_key, ordinal, left_id, right_id)`. One component row exists
even when its pair table is empty, because the carrier counts are semantic data. Pair ordinals
preserve the required left-id order and make malformed out-of-order artifacts detectable rather
than silently sorted during decoding.

### Atom-mapping experiment model

The corpus model is organized around a reusable chemical substrate and mapping-experiment facts
joined to it. The chemical substrate consists of source-associated lhs/rhs `MappingInput` values,
explicit alignments between alternative input frames, and neutral `AtomMapping` values with their
ordered atom pairs. Full molecule correspondences used by input alignments remain ordinary
`MoleculeCorrespondenceRecord` values. The input supplies each atom mapping's carrier frames and
counts; pair ordinals preserve correspondence order.

`AtomMapping` replaces the provisional candidate terminology. Candidate, optimum, selected result,
and chemically correct are roles or assessments, not properties of the stored correspondence.
Materialized `Reaction` and `ReactionSpan` values remain derived chemical views unless a later
physical-layout decision demonstrates a reason to store them.

Experiment facts retain their concrete scope. Imported mapper output records source provenance and
producer-reported values such as confidence. Network reconstruction records its exact witness.
Internal `MappingRun` values record execution facts and link their ordered outputs to atom mappings.
Algorithm-specific results retain only the values and guarantees that the selected algorithm can
actually report; proof, enumeration completeness, and multiplicity therefore do not become
mandatory fields of every run. Chemical-correctness assessment remains a separate reportable fact.

The edit decomposition becomes an on-demand `EditAnalysis` used by inspection and visualization.
It is neither an objective catalog nor a stored candidate evaluation. After S6d reconciliation, the
first cut does not retain `Objective`, generic `CandidateEvaluation`, or generic evidence relations.
The working Rust values remain provisional and do not yet declare Arrow or Parquet table names,
field codes, codebooks, or compatibility guarantees.

The first realization is a faithful record representation without denormalized molecular
projections. It does not duplicate a reaction's derived rhs, `ReactionSpan`, element composition,
entity counts, canonical form, fingerprint, or similar values merely to accelerate scans. Such
values are computed from the faithful columns or after reconstructing the selected record. Separate
corpus observations, mapping results, and provenance remain separate stored facts. Adjudications
will join in the same way when doc 205 settles their workflow and representation; none are fields or
cached projections of the molecular record. If measurements later justify
materialized molecular projections, their authority and invalidation semantics are a subsequent
design question.

The experiment should measure physical questions: nested arrays against coordinated tables where
both are plausible, scan and projection cost, materialization and selected-row reconstruction,
dictionary effectiveness, Parquet compression and row-group behavior, and the mapping and reaction
queries required by doc 205. It should inspect whether DuckDB or Polars can express those queries
directly and whether useful leaf projections and predicates reach Parquet. It should not use those
measurements to decide whether graph entities deserve explicit columns.

## Depiction and notebook output

Visualization may not depend on supplied coordinates. The first cut includes an explicit 2D layout
algorithm; the existing RDKit- or OpenBabel-class layout is sufficient for this experimental stage.
Layout remains separate from SVG encoding so that improving the algorithm does not change the
depiction or notebook contracts.

TableIR's optional positions are source-format evidence. They participate in interpreting formats
such as MOL, including resolving depiction-dependent stereochemistry during the fallible raise into
graph IR. TableIR is therefore neither the output of graph layout nor an intermediate to which a
graph-IR molecule should be lowered merely to render it.

An automatically generated 2D depiction is instead a derived presentation coupled to graph-IR
entities. Coordinate generation does not itself produce that presentation. It produces a
`MoleculeLayout`, which depiction construction combines with the graph-IR value to produce drawing
items, typed entity references, and bounds. Multiple layouts and depictions may correspond to
the same graph-IR value, and neither reconstructs or redefines that value. SVG is a rendering of the
depiction, and the Jupyter MIME adapter is a consumer of the SVG path rather than a second renderer.

The explicit coupling operation resembles the graph/geometric bridge, but the represented thing is
different. A geometric molecule is a standalone molecular-model instance with its own elements,
coordinates, charge, multiplicity, and operations; bond perception performs a genuine model
conversion back to graph IR. A molecule layout is an editable coordinate assignment in one graph-IR
atom frame, and a generated depiction is a disposable presentation of a graph-IR value. Their typed
references depend on that value, and neither has a reverse conversion to a graph model.

The first cut should therefore treat `Depiction` as a dedicated visualization boundary type in
`umol-io`, alongside but separate from TableIR, rather than reproduce the
`umol-geometric-core` / `umol-geometric` / `umol-geometric-graph` crate structure. It must not reuse
`table_ir::Molecule`: TableIR carries source-format state and participates in raising, while
`Depiction` carries graph-derived presentation state and participates in rendering. The dependency
seam is:

```text
graph IR --layout projection + algorithm--> molecule layout
graph IR + molecule layout --depiction construction--> depiction
depiction --SVG renderer--> SVG
```

Separate `umol-io` modules should own molecule layouts, the format-neutral depiction, explicitly
selected layout algorithms, depiction construction, and SVG rendering. Keeping `Depiction` as a
named public boundary rather than an SVG-private helper leaves a clean extraction path. If depiction
later becomes an independently constructed, persisted, or edited model with its own operations and
several model bridges, it can move to a representation crate plus coupling crates without changing
the layout or renderer contracts.

The first SVG renderer may build one deterministic document in a `String` and use a fixed style in
nominal-bond-length units. These are provisional renderer choices rather than `Depiction` semantics.
If measurements justify it, the same item traversal may write incrementally to a caller-supplied
text or XML sink. A later coherent SVG configuration may own font choices, stroke widths,
multiple-bond spacing, marker and arrow geometry, and viewport padding or sizing policy; these
remain renderer concerns rather than fields on `Depiction`. Molecule size determines content bounds
but should not by itself cause labels and strokes to shrink.

### Molecule layout contract

`MoleculeLayout` is a small open carrier rather than an opaque backend result or another molecular
model. It contains one finite two-dimensional position for every atom in its own dense `AtomId`
frame. Coordinates are dimensionless with a nominal bond length of one and use a mathematical
y-up convention; an SVG renderer performs the output-axis transform.

Public construction establishes only the layout's intrinsic integrity, including finite
coordinates. It cannot establish agreement with a separately supplied molecule. The first operation
that combines independently supplied values checks that the molecule's atom count equals the
layout's frame size. A layout returned for that molecule by a conforming algorithm satisfies the
same condition by construction. The carrier remains directly editable so callers may adjust atom
positions before depiction.

`MoleculeLayout` contains no atom labels, bond styles, colors, chemical interpretation, drawing
bounds, or reverse conversion to graph IR. Its atom positions are the reusable result needed by
reaction alignment and composition. Reaction layout may rotate, reflect, translate, or constrain
component layouts before depiction without manipulating SVG constructs.

CoordGen consumes an operational layout projection rather than a replacement molecule. The
projection contains topology, literal elements and bond orders where useful to the backend, generic
atoms or connectivity for underdetermined forms, and optional fixed-coordinate constraints. The
original graph-IR forms remain authoritative and are read again during depiction construction where
the first depiction projection supports them. An underdetermined or set-valued element may influence
layout as a generic atom without requiring its complete form to appear in the depiction. This
projection is a heuristic input to coordinate generation, not a lossy model conversion.

The first depiction projection includes literal element, isotope, charge, implicit hydrogen, and
localized bond order. Aromatic systems project to aromatic markers on their member atoms and induced
localized bonds. Stereo entities project to markers on their site atoms or bonds; deriving wedges,
hatches, or configuration-dependent geometry is later work. Dative, multicenter, and noncovalent
entities and arbitrary constraints may be omitted. Two graph-IR values that differ only outside this
projection may therefore have the same depiction.

The CoordGen adapter consumes the same atom/bond projection. Later layout algorithms may consume
additional overlays as hints without changing `MoleculeLayout`.

Layout preserves the supplied entity frame and does not canonicalize implicitly. The initial
determinism contract is scoped to the same graph-IR representation, algorithm, backend version, and
layout inputs. A caller that wants a canonical-frame depiction canonicalizes explicitly before
layout.

### Reaction layout

The first reaction depiction uses indexed sides. Materialization supplies the two molecular sides
and their correspondence; the molecule-layout algorithm positions each side independently; the
reaction compositor places them around a reaction arrow; and matching atom pairs receive the same
depiction-local map index. These displayed indices identify correspondence pairs. They are not
`AtomId` values or persistent molecular identifiers.

`Depict` is the end-to-end operation shared by `Molecule` and `Reaction`. Its `depict_with` method
takes the explicitly selected molecule-layout algorithm. A molecule is laid out and lowered
directly. A reaction first materializes its `ReactionSpan`, then derives the two molecular sides and
their atom correspondence from that span. Failure to materialize an internally inconsistent
reaction is reported at this lazy depiction boundary.

Independently supplied molecular sides and an atom correspondence are not themselves a reaction.
They remain useful for corpus candidates and other mapped-pair evidence, but their public operations
are named `depict_from_sides` and `depict_from_sides_with`. The latter selects the molecule-layout
algorithm explicitly. Neither operation is presented as depiction of a `Reaction` value.

Map indices are assigned consecutively from zero in the correspondence's left-id order. Unmatched
atoms receive no map index. A boundary format whose specification fixes another base converts the
indices at that boundary; it does not change the internal zero-based convention.

This mode does not require the correspondence to influence molecular coordinates. It is robust
under rearrangements and component merging or splitting, and it keeps the first reaction compositor
independent of a mapped-substructure alignment algorithm.

Correspondence-aware alignment is a useful later reaction-layout algorithm. It keeps lhs and rhs on
separate sides of the arrow but uses their atom correspondence to place conserved substructures in
similar orientations. A minimal form lays out both sides independently and applies a rigid
translation and rotation to one side. This is technically straightforward but cannot reconcile
different internal layouts of the same substructure. A stronger form constrains regeneration of one
side from coordinates on the other. It must select conserved anchor subgraphs, resolve symmetry and
competing anchors, and handle additions, removals, rearrangements, and component merging or
splitting. That selection makes it materially more involved than indexed sides.

Correspondence-aware alignment remains an explicitly selected reaction-layout algorithm rather than
silent polishing of indexed layout. It may reuse any molecule-layout backend. Overlaying both sides
in one union-frame layout and drawing correspondence connectors between independently laid-out sides
remain other possible analytical presentations; neither is part of the first cut.

### Layout and viewing backends

The ecosystem divides by responsibility rather than into interchangeable whole-system backends.
Algorithm selection should remain explicit at each applicable level:

- molecule coordinate generation;
- reaction composition and correspondence-aware alignment;
- depiction rendering; and
- interactive two- or three-dimensional viewing.

A molecule layout algorithm may therefore be used inside more than one reaction-layout algorithm.
For example, a correspondence-aware reaction compositor may arrange molecule coordinates produced
by CoordGen or Indigo. A single flat selector that treats coordinate generation, reaction alignment,
SVG rendering, and notebook interaction as variants of the same operation would obscure these
boundaries.

[CoordGen](https://github.com/schrodinger/coordgenlibs) is the leading focused molecule-layout
candidate. It is BSD-3-Clause, accepts coordinate constraints, and does not impose a second
chemical toolkit on the depiction model. It does not provide reaction composition or SVG rendering.
Its C++ API requires a small C wrapper for Rust, and adopting it would add a C++ compiler to a build
that currently requires only a C compiler.

[Indigo](https://github.com/epam/Indigo) is the strongest second native candidate. It is Apache-2.0,
has a C API, and directly lays out and renders both molecules and reactions. It is also a much larger
toolkit. Its reaction layout supplies component and arrow composition, but the inspected
implementation does not use atom correspondences to align conserved substructures. It should be
compared with a compositor that does. Indigo layout and Indigo rendering are separate experimental
roles; using one does not require making the other authoritative.

[CDK](https://github.com/cdk/cdk) is the principal reference for mapped-reaction layout. Its
`StructureDiagramGenerator` automatically aligns mapped reactant and product substructures and
allows that behavior to be selected explicitly. The JVM boundary and LGPL license make CDK
unattractive as a shipped umol backend, but its algorithm and output are important external
evidence for designing and testing correspondence-aware composition.

[CDPKit](https://github.com/molinfo-vienna/CDPKit) is another reference implementation rather than
a proposed shipped backend. Its C++ molecule and reaction layout and rendering facilities broaden
the comparison beyond CDK, while its LGPL license and full-toolkit build make it a poor initial
distribution dependency.

[OpenChemLib](https://github.com/Actelion/openchemlib) and its JavaScript build are worth retaining as
browser-side layout and editing references. [Ketcher](https://github.com/epam/ketcher) is a promising
reaction editor and possible later notebook interface, but its automatic layout is provided by
Indigo through a service or WebAssembly and is not an independent layout algorithm. Mol* and
Jmol/JSmol are candidates for later interactive 3D consumption of geometric records, not 2D
structure-diagram generation. Open Babel may remain an external comparison, but its GPL-2.0 license
makes it unsuitable as a linked or vendored umol backend. MolView is neither a reusable layout
library nor redistributable under ordinary open-source terms and is excluded.

If a native external backend survives exploration, its source should be vendored into the
repository rather than attached as a git submodule. The nauty and msym integrations establish the
desired distribution property: a checkout contains the source needed for a self-contained build
without a second repository operation. Vendoring does not remove the C++ toolchain cost, so build
portability, compiler availability, binary size, and maintenance burden must be considered alongside
layout quality before a backend becomes part of an ordinary build.

The feature-branch exploration may use C++, JVM, WebAssembly, subprocesses, or unbundled local
checkouts freely. Those experiments provide evidence and do not commit the durable crate graph or
default build to their dependencies.

The SVG payload may retain structured references for the atoms, bonds, correspondence pairs, and
changes present in the first depiction projection. Other overlay and constraint references are not
required until their depiction is attempted. Notebook code must not recover semantics by scraping
presentation-generated element names.

## Deferred reaction-template counterexample

A small reaction-template library remains a plausible second consumer, but it is no longer part of
this first experiment or S8. The sketch below preserves the boundary of that possible follow-up so
it does not become mixed into the atom-mapping and annotation work.

That consumer would connect:

- source-specific `ReactionFact` values containing concrete reactions;
- directed reaction rules represented by the same graph-IR `Reaction` type;
- the many-to-many `RuleApplication` relation between supplied rules and facts;
- substrates and products derived from each concrete reaction;
- an optional opaque condition value;
- a stable literature document identifier and a location within that document.

It would be deliberately independent of the atom-mapping experiment. Its purpose would not be
to validate or generalize the mapping study's source-case, mapping-input, run, objective,
evaluation, evidence, or assessment relations. It should reuse only the faithful molecular records,
storage machinery, graph-IR operations, and depiction facilities that prove independently useful.
No shared corpus abstraction is authorized by the present atom-mapping experiment.

`ReactionFact` is source-specific rather than a unique chemical-reaction identity. It refers to one
self-contained `ReactionRecord`, carries the optional condition value, and identifies one document
and source location. Multiple facts may contain structurally or canonically equal reactions under
different conditions or in different sources. A source-frame reaction remains authoritative for
the fact; its canonical form or experiment-local canonical key is derived for lookup and comparison.

A rule is another self-contained `ReactionRecord` referenced in the rule role. No graph-IR
`ReactionRule` type or corpus rule wrapper is required merely to assign that role. Rules are directed;
a retained reverse is a separate rule. Rule extraction and generalization are outside this
experiment: rules arrive as supplied inputs.

`RuleApplication` is the two-key join between a rule record and a reaction fact. It asserts that the
fact is an application of the rule but initially stores no match correspondence, extraction record,
or embedded application result. Insertion checks only that both references exist. A separate lazy
inspection operation applies the rule to the fact's lhs and reports whether any resulting concrete
reaction is equal to the fact reaction under the selected reaction-equivalence semantics. If later
examples demonstrate a need to distinguish embeddings, a correspondence may be added then.

Substrate and product molecules are the connected components derived from the concrete reaction's
lhs and rhs. They are not duplicated as authoritative participant rows in the first cut. Catalysts,
solvents, and other non-transformed source information remain in the provisional condition value.
Canonical forms, component decompositions, application checks, and depictions are derived results,
not stored chemical authorities.

Conditions and selectivity are not designed here. Preserve one optional free-format string or opaque
structured payload and let the experiment reveal representative condition data without claiming to
normalize it. Do not introduce a curation workflow, user-management model, assertion history,
participant-role ontology, rule-extraction lineage, or multi-layer evidence system.

Bibliography is similarly shallow: retain a stable document identifier plus a simple location
reference such as a page, scheme, table, example, or record label. Nothing more is required to test
the data path.

A minimal corpus decomposition therefore needs only:

- faithful self-contained reaction records containing `lhs + deltas`;
- `ReactionFact(fact_key, reaction_record_key, condition, document_id, document_location)`;
- `RuleApplication(rule_record_key, reaction_fact_key)`.

This is a test fixture for the substrate, not a proposed final schema. Its value is precisely that it
places a very small, structurally different consumer beside atom mapping and exposes which proposed
abstractions were actually use-case-specific.

## First experimental slice

The first slice should answer one practical question: can the atom-mapping corpus move through a
coherent external data path from molecular values and correspondences, through a faithful boundary
and Arrow/Parquet storage, to direct queries and inspectable SVG output?

The slice uses the faithful-boundary contract and exact graph-IR structural equality. The columnar
projection implements that settled contract rather than defining it.

The slice should include representative `Molecule`, `Reaction`, `ReactionSpan`, and atom-
correspondence values; the complete logical columnar decomposition needed to represent them;
source-specific corpus provenance; one concrete mapping query; deterministic 2D layout; SVG output;
and a Jupyter display path. Adjudication is deferred to the annotation workflow in doc 205 rather
than represented by empty tables in this slice. The template-library records remain a second nearby
consumer that checks whether the substrate has been specialized too narrowly to atom mapping.

Measurements cover batch-construction time, Arrow and Parquet size, write and read time, scan and
projection time, allocations, selected-row reconstruction, materialization cost, layout and
rendering cost, and the hard-tail populations already identified in doc 205. The useful unit is
corpus construction, persistence, retrieval, query, and visual inspection through one shared
representation path.

## Parked questions

The following are important parts of the direction in doc 200, but the first exploration does not
settle or prototype them:

- a durable stable identifier or cross-version canonicalization profile;
- a complete condition, selectivity, outcome, or bibliography schema;
- user management, curation workflow, or assertion history;
- geometric records or 3D export;
- ontologies or executable classification;
- underdetermined query semantics;
- native DuckDB, Polars, or DataFusion chemical operators;
- billion-record corpus organization, compaction, remote access, or schema evolution;
- reaction-network persistence;
- a store-relative reaction map linking separately stored molecules and correspondences;
- row-oriented Postcard, FlatBuffers, or purpose-built codecs and point stores;
- a public borrowed-view or columnar-batch API;
- replacement of `Reaction` by `ReactionSpan` as the central representation.

The Arrow or Parquet schema may identify generated experimental artifacts adequately. If an
explicit format version proves necessary, one dataset-level marker is sufficient. Compatibility
policy is not required for generated experimental artifacts.

Parked means awaiting concrete evidence, not rejected or considered secondary. Each question should
return when the experimental path supplies representative states, failure modes, or competing
operations against which a design can be judged.

### Promotion guard for schema coverage

Schema-change detection is not part of the experiment, but it must be settled before promoting
`umol-store` from provisional feature work to the development trunk. The required property is that
adding a graph-IR entity family, faithfully stored field, entity-constraint variant, or molecule-
constraint variant cannot leave the record or columnar schema silently incomplete.

The handwritten S0 and S1a implementation provides only partial compile-time detection. Exhaustive
record-conversion matches reject many new enum variants, and explicit struct construction may expose
some new fields or aggregate families. The record types, storage codebooks, `LeafKind`,
`EntityConstraintKind`, `SparsePayloadKind`, and Arrow schema constructors are nevertheless separate
catalogues. After the immediate record-layer errors are handled, omitting a corresponding codebook
or Arrow field can still compile, while exact schema tests only freeze the catalogue already named
by the test.

Before promotion, choose a coverage mechanism that makes every relevant graph-IR extension fail at
compile time or in one dedicated completeness test until its record conversion, codebook, leaf or
sparse schema, aggregate layouts, and roundtrip coverage have been updated. A single declarative
inventory or code generation and exhaustive adapters plus cross-inventory tests are possible
approaches. This document does not select or schedule either mechanism.

## Findings needed for the next iteration

At the end of the spike, record what was learned without requiring every question to have an answer:

1. Which graph-IR states or entity frames, if any, make faithful boundary or columnar encoding
   awkward or costly?
2. Does exact structural roundtrip hold across every representative boundary value, Arrow
   projection, Parquet roundtrip, and corpus case?
3. Does a canonical content key provide useful deduplication within one version and context?
4. Does graph-IR `Reaction` encode both concrete examples and rules without storage-specific
   distortion?
5. Is keeping `ReactionSpan` derived inconvenient for any measured operation?
6. Which physical columnar layout best realizes the already-settled logical entity decomposition?
7. Does materializing owned molecules dominate query cost at the explored scale?
8. Does the depiction scene cleanly separate layout, SVG, and notebook behavior, and does any
   observed use justify extracting it from `umol-io` as an independent model?
9. Did an existing metadata schema help, or did free-format conditions and simple references suffice?
10. Which question, if any, is now concrete enough to deserve a second experiment or a design
    decision?

These are observation prompts, not acceptance criteria. The findings determine which semantic or
implementation question becomes concrete enough to study next. The boundary semantics, physical
comparison, depiction projection, and first consumers are now settled sufficiently to sequence the
experiment without treating its outcome as predetermined.

[205](205-mapping-test-corpus-2026-08-20.md) applies this exploratory direction to atom mapping.
It treats the mapping test corpus as the first demanding case study for the same columnar, query,
and notebook path rather than constructing a separate research-only data stack.

## Staged implementation plan

The plan has two initially independent foundation branches. S0-S3 establish faithful records and
the two columnar realizations. S4-S5 establish layout, depiction, SVG, and notebook display. S6 joins
them in the atom-mapping corpus and S7 records its measurements; S8 tests molecule generation,
database-backed deduplication, and retrieval at representative scale. Mandatory S9 then assigns
every provisional surface an explicit disposition. Every physical-layout or layout-algorithm
choice remains explicit at the operation boundary.

Unless a subitem says otherwise, it is additive and the tree remains green after the subitem. Each
subitem includes its tests rather than leaving verification to a later hardening stage. Property
tests assert exact structural laws; example tables cover named variants and failures. Benchmarks use
fixed constructed values or caller-supplied generated artifacts, never a hidden runtime dependency
on `materials/`.

### S0 — Faithful owned records

S0 establishes the semantic boundary before any Arrow schema can accidentally define it.

- **S0a — `umol-store::record::form`: form records and conversion vocabulary.** Add the new
  `umol-store` workspace crate, its `record` module, the owned tagged forms needed by graph IR, and
  the unit-context `FromIr`/`IntoIr` implementations. The records represent `Undetermined`, literal,
  set, range, variable, arithmetic, predicate, and stereo-term forms without normalization. Exact
  variant tables and property tests cover `ir -> record -> ir == ir`, including every rare recursive
  form. Add the first conversion benchmark cases here. **Additive (green). Done.** [dep: none]
- **S0b — `umol-store::record::molecule`: entity and molecule records.** Add the eight typed entity
  record families, scalar `atom_0_id`/`atom_1_id` endpoint fields, typed members, entity constraints,
  molecule constraint trees, and `MoleculeRecord`. Preserve local ids, entity order, participant
  order, repeated constraints, and exact graph-IR forms. Example tests exercise all eight entity
  families and generated property tests assert exact molecule roundtrip. **Additive (green).**
  **Done.** [dep: S0a]
- **S0c — `umol-store::record::reaction`: reaction records.** Add a self-contained, closed
  `ReactionRecord` with public read-only accessors for its lhs and ordered deltas. Preserve one flat
  ordered delta sequence, complete add/remove/modify payloads, and recursive constraints only where
  carried by `ConstraintDelta`. The initially added public aggregate fields were removed after the
  design-fidelity review; no public assembly constructor is exposed. Tests distinguish reordered
  deltas, exercise the accessors, and cover every delta family and variant; property tests assert
  exact reaction roundtrip. **Additive (green). Done.** [dep: S0b]
- **S0d — `umol-store::record::reaction_span`: span records.** Add the molecule-like union-frame
  `ReactionSpanRecord`, retaining each entity span tag and the exact one-sided or lhs/rhs attributes.
  Tests cover all entity families, all valid span tags, and constraint-span ordering; property tests
  assert exact span roundtrip. **Additive (green). Done.** [dep: S0b]
- **S0e — `umol-store::record::correspondence`: correspondence records.** Add
  `MoleculeCorrespondenceRecord` with eight named typed components, `UInt64`-domain carrier counts,
  and pairs in left-id order. Tests cover empty carriers, unmatched entities, non-total mappings, and
  all eight components; property tests assert exact correspondence roundtrip. **Additive (green).
  Done.** [dep: S0b]
- **S0f — `umol-store/benches/record.rs`: record baselines.** Benchmark conversion for small literal,
  constraint-heavy, recursive-form, reaction, span, and correspondence fixtures. Record allocations
  and throughput separately so later direct-to-Arrow paths can be judged as optimizations of the
  same semantics. A benchmark smoke test verifies every fixture roundtrips before timing it.
  **Additive (green). Done.** [dep: S0b, S0c, S0d, S0e]

**S0 gate:** `cargo test -p umol-store`, its feature-gated property suite, and the record benchmark
smoke path pass. No Arrow or Parquet dependency is needed to use the owned records.

### S1 — Shared Arrow leaves and nested batches

S1 implements one complete Arrow realization while separating the shared leaf encoding from its
aggregate ownership shape.

Nested aggregate batches are producer-issued physical values in this experiment. Their initially
added public `try_new(RecordBatch)` adoption constructors were removed after the design-fidelity
review because schema equality alone cannot establish record integrity. `from_records` is the public
production path, and wrapper-level zero-copy `slice` operations preserve the producer-established
contract. Admission of arbitrary Arrow columns remains a later, separately designed reader boundary.

- **S1a — `umol-store::arrow::schema`: codebooks and shared leaf schemas.** Add explicit `UInt8`
  codebooks, field names, widths, nullability, typed sparse-payload schemas, and one schema-construction
  API shared by both layouts. Do not derive semantic codes from Rust discriminants or dictionary
  indices. Exact schema tests assert every field's name, type, nullability, and code assignment.
  **Additive (green). Done.** [dep: S0a]
- **S1b — `umol-store::arrow::leaf`: form and constraint columns.** Add builders and readers for the
  dense literal/undetermined path and typed sparse payload tables, including ordered expression and
  constraint trees. Decoder errors are limited to mechanical schema, offset, tag, and range failures;
  they do not normalize or chemically validate. Variant tables, malformed-column examples, and
  property tests assert `record -> Arrow -> record == record`. **Additive (green). Done.** [dep: S1a]
- **S1c — `umol-store::arrow::nested::molecule`: nested molecule batches.** Add one-row-per-record
  nested molecule encoding with explicit typed ids and record-owned child lists. Tests cover empty
  and mixed-sized batches, multiple records with repeated local ids, nonliteral payloads, and exact
  selected-row reconstruction. **Additive (green). Done.** [dep: S0b, S1b]
- **S1d — `umol-store::arrow::nested::reaction`: nested reaction batches.** Add the self-contained
  lhs struct and one ordered delta list with typed variant payloads. Tests prove that delta order and
  complete payloads survive batch slicing and reconstruction. **Additive (green). Done.** [dep: S0c, S1c]
- **S1e — `umol-store::arrow::nested::{reaction_span, correspondence}`: remaining nested values.**
  Add union-frame span columns and eight named correspondence component structs. Tests cover all span
  tags, empty pair lists with nonzero counts, and exact reconstruction after batch slicing.
  **Additive (green). Done.** [dep: S0d, S0e, S1c]
- **S1f — `umol-store/benches/arrow.rs`: nested construction and materialization.** Benchmark batch
  construction, leaf projection, selected-row reconstruction, whole-batch materialization, and
  allocation counts over the S0 fixture families and scalable literal-only values. **Additive
  (green). Done.** [dep: S1c, S1d, S1e]

**S1 gate:** nested batches preserve exact records and exact graph-IR roundtrip across the example
and property suites; the Arrow benchmark target builds and its smoke fixtures pass.

### S2 — Coordinated Arrow tables

S2 realizes the same leaves as typed coordinated tables. It does not introduce an entity-value
table, owner-kind tag, or a second logical schema.

- **S2a — `umol-store::arrow::coordinated::molecule`: coordinated molecule tables.** Add molecule,
  eight entity, member, constraint, and typed payload table builders/readers keyed by
  `(record_key, local_id)` or `(record_key, payload_id)`. Explicit ordinals reconstruct every stored
  sequence without relying on engine row order. Tests permute physical rows where permitted and
  assert exact semantic reconstruction. **Additive (green). Done.** [dep: S0b, S1b]
- **S2b — `umol-store::arrow::coordinated::reaction`: coordinated reaction tables.** Add a distinct
  reaction-lhs frame, central ordered delta table, and typed variant payload tables keyed by
  `(record_key, ordinal)`. Tests reject missing or duplicate ordinal references mechanically and
  preserve exact delta order. **Additive (green). Done.** [dep: S0c, S2a]
- **S2c — `umol-store::arrow::coordinated::{reaction_span, correspondence}`: remaining table
  families.** Add distinct span-frame tables plus the eight component and eight pair table families.
  Tests cover all span tags, carrier counts independent of pair counts, and out-of-order pair
  detection. Factored entity-constraint rows identify their `one_sided`, `lhs`, or `rhs` span
  branch with a stable code; pair rows retain explicit ordinals and decoding also checks strictly
  increasing left-id order. **Additive (green). Done.** [dep: S0d, S0e, S2a]
- **S2d — `umol-store::arrow::equivalence`: cross-layout laws.** Add a property suite that assigns
  fresh physical keys independently, encodes the same record through both layouts, and proves that
  both reconstruct the same record and graph-IR value. Assertions also prove that equality, hashing,
  and deduplication ignore reassigned physical keys. **Additive (green). Done.** [dep: S1e, S2c]
- **S2e — `umol-store/benches/arrow.rs`: coordinated comparison.** Add matching construction,
  projection, scan, selected-row reconstruction, whole-batch materialization, and allocation cases
  using the exact S1 fixtures and leaf encoding. **Additive (green). Done.** [dep: S1f, S2d]

**S2 gate:** both Arrow layouts satisfy one cross-layout roundtrip suite, and the benchmark harness
can compare them without changing input generation or leaf semantics.

### S3 — Parquet datasets and direct scans

S3 turns each Arrow realization into an explicitly selected on-disk artifact. It does not add a
row codec, transparent layout default, or mixed-version reader.

- **S3a — `umol-store::dataset`: dataset selection and metadata.** Add an explicit
  `ColumnarLayout::{Nested, Coordinated}` selector, one experimental schema identity if the reader
  needs it, dataset table naming, and mechanical `StoreError` cases. Physical keys are allocated per
  artifact and remain outside record equality. Tests cover deterministic table discovery, explicit
  layout selection, and rejection of absent or mismatched experimental metadata. **Additive
  (green). Done.** [dep: S1e, S2c]
- **S3b — `umol-store::parquet::nested`: nested artifact writer and per-type readers.** Add a
  stateful writer that creates one artifact, accepts any caller-selected nonempty subset of record
  families through independent iterator operations, allocates keys per family, and writes the
  manifest with the actual table inventory only on consuming `finish`. An explicitly selected empty
  family remains distinct from an absent family. Keep readers per record type and distinguish an
  unselected family from a declared table with missing Parquet parts. Tests cover one-family and
  multi-family artifacts, empty and absent families, duplicate selection, temporary multi-file
  datasets, row-group boundaries, selected-leaf projection, and exact reconstruction for all four
  record types. **Additive (green). Done.** [dep: S3a]
- **S3c — `umol-store::parquet::coordinated`: coordinated dataset writer and reader.** Persist every
  typed table belonging to each caller-selected record family under the same artifact-writer and
  manifest-inventory contract as S3b, and recover records through per-type readers using explicit
  keys and ordinals. Tests roundtrip temporary multi-file datasets, reorder table fragments, and
  verify that compaction-style physical key reassignment leaves reconstructed records unchanged.
  **Additive (green). Done.** [dep: S3a]
- **S3d — `umol-store::scan`: concrete projection and retrieval operations.** Add only the operations
  needed by the first experiment: enumerate artifact-local record keys, project selected literal
  fields, and reconstruct a selected molecule, reaction, span, or correspondence. Matching still
  materializes graph IR and uses existing graph-IR semantics. The same operation table tests both
  layouts. **Additive (green). Done.** [dep: S3b, S3c]
- **S3e — `umol-store/benches/parquet.rs`: persistence and scan baselines.** Measure Arrow and Parquet
  size, write/read time, leaf projection, selected-row retrieval, full materialization, allocations,
  compression, and row-group sensitivity for both layouts. Benchmark setup validates exact
  roundtrip before collecting timings. **Additive (green). Done.** [dep: S3d]

**S3 gate:** temporary Parquet artifacts in both layouts reconstruct identical records; the complete
storage benchmark matrix runs from one input generator and reports physical size as well as time.

### S4 — Format-neutral depiction and SVG

S4 is independent of S0-S3. It establishes the graph-IR-coupled presentation path before adopting a
native layout backend.

- **S4a — `umol-io::layout`: `MoleculeLayout`.** Add the open dense `AtomId`-frame coordinate carrier,
  y-up and nominal-bond-length conventions, intrinsic finite-coordinate construction, direct editing,
  and the contextual frame-size check used when combining it with a molecule. Tests cover zero atoms,
  nonfinite input, post-construction edits, and frame mismatch without adding chemistry validation.
  **Additive (green). Done.** [dep: none]
- **S4b — `umol-io::depiction`: depiction scene vocabulary.** Add format-neutral atom, bond, text,
  marker, arrow, bounds, and typed graph-reference items plus `Depiction`. Construction remains
  independent of SVG element names. Exact example tests cover bounds and typed references.
  **Additive (green). Done.** [dep: S4a]
- **S4c — `umol-io::depiction::molecule`: graph-IR lowering.** Combine a molecule and layout into the
  settled first projection: literal element, isotope, charge, implicit hydrogen, localized bond
  order, aromatic markers, and stereo site markers. Omit unsupported overlays and constraints
  explicitly. Tests use paired molecules that differ inside and outside the projection and assert
  the documented visible/omitted behavior. **Additive (green). Done.** [dep: S4b]
- **S4d — `umol-io::svg`: SVG renderer.** Render `Depiction` with deterministic item ordering,
  y-axis conversion, escaping, view bounds, and structured references for represented entities.
  Tests parse the output as XML and make exact assertions about items, references, and escaped
  text rather than relying only on whole-file snapshots. Add rendering benchmarks. **Additive
  (green). Done.** [dep: S4b, S4c]
- **S4e — `umol-io::depiction::reaction`: indexed-side composition.** Add the explicitly selected
  `depict_from_sides` reaction compositor over already supplied molecule layouts and an atom
  correspondence. Place sides around an arrow and assign zero-based consecutive map indices in
  left-id order. Tests cover unmatched atoms, component changes, reordered ids, and deterministic
  index assignment. **Additive (green). Done.** [dep: S4c, S4d]

**S4 gate:** callers can construct/edit a layout, lower the supported molecule or indexed reaction
projection, and render deterministic inspectable SVG without CoordGen or Python.

### S5 — CoordGen layout and Jupyter display

S5 supplies the first automatic layout algorithm while keeping the backend projection, presentation,
and notebook surface separate.

- **S5a — `umol-coordgen-sys`: vendored CoordGen interface.** Add a vendored source snapshot and its
  BSD-3-Clause notice, a minimal C ABI wrapper, and the C++ build integration needed only by the
  experimental layout feature. Do not use a git submodule. ABI smoke tests cover an empty graph, one
  atom, and a bonded graph; CI/build documentation records the C++ compiler requirement.
  **Additive (green). Done.** [dep: none]
- **S5b — `umol-io::layout::coordgen`: projection and adapter.** Add the operational topology/literal
  projection and `MoleculeLayoutAlgorithm::CoordGen`; callers pass the selector explicitly. Generic
  atoms represent unsupported element forms only for layout input, while graph IR remains
  authoritative for depiction. Tests assert exact atom-frame preservation, finite normalized
  coordinates, deterministic results for a fixed backend, and no implicit canonicalization.
  **Additive (green). Done.** [dep: S4a, S5a]
- **S5c — `umol-io/benches/layout.rs`: layout evidence.** Benchmark representative acyclic, cyclic,
  aromatic, disconnected, underdetermined, and mapping hard-tail molecules. Store backend/version
  information with results and validate frame integrity before timing. **Additive (green). Done.**
  [dep: S5b]
- **S5d — `umol-io::depiction`: automatic molecule and reaction paths.** Add `Depict::depict_with`
  for `Molecule` and `Reaction`, taking the explicit molecule-layout selector. Reaction depiction
  materializes its span and derives its own sides and correspondence. Retain
  `depict_from_sides_with` for independently supplied mapped sides; correspondence affects its
  labels, not the generated coordinates. Tests compare automatic paths with separately laid-out
  composition, assert the same map-index semantics, and cover reaction-materialization failure.
  **Additive (green). Done.** [dep: S4e, S5b]
- **S5e — `umol-py::depiction`: SVG display boundary.** Bind explicit molecule/reaction layout and
  depiction entry points and return a small SVG display value whose `_repr_svg_` exposes already
  rendered SVG. Do not put a hidden layout choice inside `Molecule.__repr__` or `Reaction.__repr__`.
  The Python surface is gated by the non-default `depiction` feature, leaving the ordinary graph
  binding build independent of CoordGen and its C++ toolchain. Rust/Python parity tests cover
  explicit selection, errors, and exact MIME output; the Python 3.13 build and pytest gate is
  mandatory. **Additive (green). Done.** [dep: S4d, S5d]

**S5 gate:** a graph-IR molecule or mapped reaction reaches Jupyter SVG through an explicitly chosen
CoordGen layout and indexed-side reaction composition. The ordinary non-C++ workspace path remains
available when the experimental layout feature is disabled.

### S6 — Atom-mapping corpus vertical slice

S6 joins storage and depiction in `corpora/atom-mapping`. The corpus is a non-published research
consumer, not a new runtime data directory or a second molecular model.

- **S6a — `corpora/atom-mapping::schema`: provisional corpus relation values.** Introduced the
  source-case, input, alignment, mapping-assertion, pair, objective, evaluation, class, run, and
  statistic values initially required by S6b and S6c. The initial Arrow tables, stable table names
  and status codes, unresolved provenance/preprocessing catalogs, generic annotation/evidence rows,
  candidate-owned producer and confidence, duplicated carrier counts, and run-owned generic best
  score were removed after S6d exposed that the study model was not settled. The remaining test
  surface was provisional logical Rust data at that stage. S6d subsequently replaced candidate
  terminology with neutral atom mappings, added the concrete physical corpus schema, and removed the objective
  catalog, stored evaluations and equivalence classes, sparse run statistics, and universal proof
  and enumeration status before persistence. **Breaking (red→green). Done.** [dep: S0e]
- **S6b — `corpora/atom-mapping::ingest`: source and atom-mapping ingestion.** Ingest the existing
  graph-IR/Rhea fixtures and external mapping results into explicit source frames, inputs,
  alignments, and atom mappings. Imported confidence remains producer output rather than an
  objective or constraint; S6d retains it in the concrete RXNMapper observation relation. Tests
  pin frame transformations and the known Rh_20309, Rh_63116,
  Rh_12116, hydrogen, and parent Diels-Alder diagnostic cases without treating one arbitrary
  symmetry representative as uniquely correct. **Additive (green). Done.** The ingestion aggregate
  keeps source and mapper frames distinct and accepts only explicit frame alignments. The saved
  Rhea/RXNMapper CSV population currently
  yields 15,937 raised source frames and 14,418 atom mappings; the other 12 mapper rows refer to
  source rows that do not raise through the current boundary. The saved mapper CSV has no confidence
  column. When present in another imported file, confidence remains available to the inspection
  report as producer-reported output metadata but is not attached intrinsically to a mapping or
  treated as chemical-correctness evidence.

  The classification population is a retained reconstruction operation, not another downloaded
  mapper result. Its index reader selects the 194,778 nonempty `valid` or `reducible` endpoint
  subgraphs, returns to each original directed network, follows every shortest path present in the
  archived endpoint GraphML, and replays the rule in its recorded direction or reversed direction
  as required. The old rule strings are interpreted as mapped queries: outer component grouping and
  SMARTS `#1` hydrogen spelling are lowered at this boundary, and an omitted bracket-atom hydrogen
  term does not become a zero-hydrogen constraint.

  Node SMILES are resolved and every implicit hydrogen is expanded into an ordinary graph-IR atom
  before the mapping input is stored. Replay uses an operation-local element/charge/bond projection
  of those same atom frames so current derived valence fields do not redefine the historical rule
  state. Every exact intermediate frame alignment is explored during composition but is not stored
  as a corpus atom mapping. All shortest paths contribute to one endpoint set, which is deduplicated
  by exact initial-to-final atom-pair sequence before insertion. Differently labeled endpoint maps
  remain distinct even when molecular automorphisms make them symmetry-equivalent; no symmetry
  representative is selected. The endpoint subgraph remains the source case and contains the path
  graph that produced the set. S6d retains each supporting shortest path and ordered archived-rule
  sequence as a self-contained reconstruction witness joined to the neutral endpoint mapping.

  The checked-in command accepts explicit `--start <INDEX>` and `--count <COUNT>` options so
  reconstruction is restartable before S6d supplies persistent dataset output. Its optional
  `--output <REPORT.html>` projection analyzes every selected endpoint mapping and renders its
  indexed-side SVG, exact atom-id pairs, edit decomposition, and source provenance for human
  spot-checking; the report is a presentation artifact, not a corpus storage representation. The
  first ten selected subgraphs yield 226 distinct endpoint mappings and 1,854 endpoint atom-pair
  rows. In the first five-atom case, the 20 mappings witnessed over five shortest paths collapse to
  four distinct initial-to-final correspondences.

  The companion `inspect_rhea_rxnmapper` command accepts caller-selected Rhea source and RXNMapper
  result CSVs with the same explicit start, count, and HTML-output controls. It renders the stored
  mapper frame through the shared report path and includes the original and mapped reaction SMILES,
  external confidence when present, exact atom pairs, and edit decomposition. Rows that cannot
  pass the existing source or mapper ingestion boundary are reported and counted rather than
  silently entering the review set. [dep: S6a]
- **S6c — `corpora/atom-mapping::evaluate`: edit analysis and run capture.** Reused the current
  mapping and algorithm APIs to check mapping compatibility, compute edit components, and capture
  exact runs and their complete labeled outputs. Compatibility means that the atom correspondence uses
  the mapping input's frames and every matched atom pair has compatible element forms. Ingestion
  already establishes the carrier counts, id ranges, and partial-bijection shape. A compatible atom
  mapping is not thereby maximum-cardinality, globally optimal, chemically correct, or a selected
  symmetry representative. Global optimality can be established only by an applicable completed
  exact run; chemical correctness remains a separate assessment.

  An initial implementation also stored induced-reaction canonical-equivalence classes. S6d removed
  that derived relation; a later comparison may compute it without discarding the underlying labeled
  mappings. No second hidden mapping algorithm is added under the corpus tests. Example and property
  tests compare complete labeled optimum sets where exact enumeration is available. **Additive
  (green). Done.** The corpus-local operation decomposes localized-bond edits without making an
  optimality claim; this research operation does not add a graph-IR API.

  The implementation initially persisted these derived values as candidate evaluations and attached
  proof and enumeration status to every run. S6d supersedes that storage shape: edit decomposition
  becomes an on-demand `EditAnalysis`, while optimality, enumeration, and multiplicity move to the
  specific internal algorithm result that can establish them. The completed example and
  feature-gated property tests remain evidence for compatibility, edit decomposition, symmetric
  labeled enumeration, and agreement of the complete labeled optimum sets produced by
  branch-and-bound and exhaustive search over generated small graph-IR molecule pairs. [dep: S6b]
- **S6d — `corpora/atom-mapping::{schema,dataset}`: atom-mapping study artifact.** Persist the
  atom-mapping evidence and algorithm-development use case rather than a generic collection of
  all four `umol-store` record families. A study is organized around a specific source reaction and
  its explicit lhs/rhs molecular frames. The chemical substrate consists of `MappingInput`, explicit
  frame alignments, neutral `AtomMapping` values, and faithful molecule or correspondence records.
  `Reaction` and `ReactionSpan` remain derived chemical views unless later measurements justify
  materializing them.

  Reconcile the existing S6a-S6c Rust values before persistence. Rename candidate keys, rows, pairs,
  and accessors to atom-mapping terminology. Remove `ObjectiveRow`, `ObjectiveKey`, stored candidate
  evaluations, and their keys. Compute the concrete `EditAnalysis` on demand from a well-formed atom
  mapping and its input for inspection and visualization. It reports compatibility and localized-
  bond edits; it does not search, claim a minimum, or become a stored artifact.

  Retain `MappingRunRow` for internal mapping executions without assuming one objective model or one
  kind of result. A run records its input, algorithm and complete applicable parameters, execution
  outcome, timing, and emitted atom mappings. Objective, optimality, enumeration completeness,
  multiplicity, confidence, and similar values belong to a concrete algorithm-specific result only
  when applicable. In particular, proof and enumeration status are not mandatory universal run
  fields.

  Imported RXNMapper, Indigo, or other mapper output is a concrete source observation rather than an
  internal run. Producer-reported confidence belongs to that output observation and is not intrinsic
  to the atom mapping. Exact quasireaction reconstruction retains its own witness and remains
  distinguishable from provisional external output. Do not introduce a generic source, evidence, or
  annotation framework merely to place these unlike facts in one table.

  Chemical correctness is a separate reportable assessment, not algorithm evidence. Reference
  support, manual adjudication, ambiguity, synthesized conclusions, and development notes may target
  an atom mapping or a set of mappings. Neither an edit score nor a proved algorithmic optimum
  implies the assessment. A synthesis may retain more than one labeled mapping or a result stated
  under a named equivalence relation. Their concrete relations and physical storage are deferred to
  the annotation workflow in doc 205; S6d does not create speculative empty tables for them.

  **Implemented reconciliation and persistence:** the speculative physical corpus schema and
  metadata attached to the wrong facts are gone. The Rust surface now uses `AtomMapping` and `MappedAtomPair`, derives
  `EditAnalysis` on demand, retains generic execution facts in `MappingRunRow`, and records optimum
  score and complete labeled multiplicity only in `MinimumEditRunResult`. The objective catalog,
  stored evaluation and equivalence relations, universal proof/enumeration status, and generic run
  statistics have been removed. `RheaReactionSource` retains the source row; `RxnMapperOutput`
  retains the producer version, raw mapped reaction SMILES, and optional producer confidence.
  Quasireaction reconstruction stores every distinct supporting shortest path as ordered node ids
  plus the ordered archived rule strings and application directions. Multiple witnesses may join
  to one neutral endpoint mapping.

  The unmapped subject remains an explicit lhs/rhs frame because graph-IR `Reaction` already
  includes a transformation. An atom mapping may induce a `Reaction` and `ReactionSpan` for
  analysis and depiction, but neither derived value becomes an independent stored authority
  without a separate decision. Write the corpus relations and only the faithful record families
  actually referenced by them through both explicitly selected `umol-store` layouts. Generated
  artifacts go only to caller-selected paths. Temporary-dataset tests verify exact reference
  resolution, preserve the distinctions among source observations, reconstruction witnesses, and
  algorithm-specific results, cover multiple equivalent labeled mappings, and reconstruct every
  selected faithful record exactly. This subitem does not add chemical-assessment or synthesis
  schema; those concerns remain in doc 205. Reaction-template, condition, and literature-extraction
  schemas remain deferred.

  The physical artifact places queryable, partitioned Parquet relation tables beside a delegated
  `umol-store` record artifact. Nested and coordinated record layouts are both explicit options;
  molecule correspondences are omitted when no relation references them. A small corpus-specific
  JSON manifest records only schema generation, selected record layout, table inventory, and
  optional-family presence. It is written last as the completion marker; it is not a row store,
  chemistry authority, transaction log, or replacement for the Parquet columns. Temporary-artifact
  tests force multiple files, cover both layouts and all three evidence kinds, preserve labeled
  alternatives and reconstruction witnesses, resolve every key, reconstruct every faithful record
  exactly, and verify omission of an unreferenced correspondence family. **Breaking (red→green).
  Done.** [dep: S3d, S6c]
- **S6e — `corpora/atom-mapping::query`: concrete comparison queries.** Add reproducible DuckDB and
  Polars query scripts for mapping disagreements, edit-analysis distributions, derived symmetry
  classes, algorithm-specific optimality or multiplicity where applicable, and hard-tail
  populations. The scripts operate directly on Parquet columns; tests compare their
  result tables with the Rust analysis path on a small committed fixture, without
  making either engine a runtime dependency of `umol-store`.

  The completed query project keeps DuckDB and Polars as isolated, locked Python dependencies of
  the corpus experiment. Both engines normalize RXNMapper atom ids through the stored lhs and rhs
  input alignments before comparing the resulting source-frame pair set with every labeled output
  of an applicable exact run. Unaligned outputs and aligned inputs without an exact run remain
  visible with explicit non-comparison statuses. They also query the algorithm-specific optimum
  and multiplicity relation and rank runs by observed runtime. Edit decomposition remains derived:
  each engine joins
  the coordinated atom, localized-bond, input, and mapped-pair tables and reproduces the Rust
  `EditAnalysis` distribution without introducing a candidate-evaluation relation. This first
  columnar expression supports the expected `Undetermined` and `Lit` element and bond-order leaves
  and rejects other variants rather than assigning them approximate query semantics.

  “Symmetry class” is not left unnamed. The checked projection uses full induced-reaction canonical
  equality, one of the distinct relations identified in doc 205. Rust derives its membership rows
  from graph-IR reactions into a disposable Parquet query input; those rows are not part of the
  atom-mapping dataset or a replacement for the retained labeled mappings. The committed fixture
  contains both physical record layouts, one symmetric RXNMapper result that belongs to the exact
  optimum set, one three-edit result that disagrees with the exact two-edit family, and an empty
  correspondence at the zero-compatible-atom boundary. Unaligned and no-exact-run rows exercise
  both explicit non-comparison statuses. A single fixture command regenerates both layouts, the
  named class projection, and `expected-query-results.json`.
  DuckDB and Polars tests compare all five result tables with those expected results. **Additive
  (green). Done.** [dep: S6d]
- **S6f — `umol-io::depiction::molecule`: imported-aromatic projection repair.** The Rhea reports
  expose a depiction defect at the S6 boundary: reaction-SMILES/TableIR raise represents imported
  aromatic atoms and bonds through definite aromatic constraints, while S4c emits aromatic markers
  only for `AromaticSystem` overlays. Constraint-only aromatic systems consequently appear as
  ordinary localized single-bond structures. Extend the first projection to mark definite aromatic
  atom and bond constraints as well as aromatic-system members, without resolving aromaticity,
  changing graph IR, or duplicating markers when both representations apply. Exact examples cover
  constraint-only, overlay-only, combined, nonaromatic, and underdetermined cases. This changes no
  public API.

  The depiction marks an atom as aromatic when its aromatic-valence constraint says that it is
  aromatic, even if the exact valence is unknown. It marks a bond as aromatic only when the bond's
  aromatic constraint is `true`. Constraints that say nonaromatic or leave aromaticity
  undetermined do not produce markers. Markers are ordered by atom or bond id. If an atom or bond is
  aromatic both through a constraint and through an `AromaticSystem` overlay, the depiction draws
  one marker linked to both. It reads the stored constraints and overlays directly; it does not run
  aromaticity resolution or mutate the molecule. **Additive (green). Done.** [dep: S4c, S6b]
- **S6g — `umol-py::store` and `corpora/atom-mapping/notebooks`: atom-mapping case inspection.** Add
  the narrow read-only Python adapter needed to open an explicitly selected experimental layout and
  reconstruct a molecule or reaction by record key. The notebook queries Parquet, retrieves the
  case's molecular sides, and displays the induced reaction through S5. Python tests use temporary
  datasets; a notebook execution test runs the small fixture without network access. This remains
  read-only inspection; the annotation workflow and the assessment/synthesis relations it requires
  belong to doc 205.

  The completed Python surface is intentionally smaller than the storage crate. `ColumnarLayout`
  names the required nested or coordinated selection, and `StoreDataset.open(root, layout)` returns
  a read-only handle after checking the artifact's completion manifest, schema generation, and
  selected layout. Its only value operations are `read_molecule(record_key)` and
  `read_reaction(record_key)`. They reconstruct the existing Python graph-IR values; they do not
  expose record constructors, writers, scans, or query-engine wrappers. A key absent from a present
  family raises `KeyError`; artifact, layout, schema, and decoding failures raise `RuntimeError`.
  Record keys remain artifact-local physical addresses.

  `inspect_atom_mapping_case.ipynb` queries the existing atom-mapping relation tables with DuckDB
  for a named source case and one concrete mapping-run output. It follows the returned molecule
  record keys, reconstructs the lhs and rhs through `StoreDataset`, builds the atom correspondence
  from the queried pair rows, induces the graph-IR `Reaction`, and displays the SVG produced by the
  explicit CoordGen selector. The committed execution check runs those ordinary code cells directly
  against the small query fixture and asserts the queried keys, mapping pairs, reconstructed carrier
  sizes, reaction arrow, and correspondence labels. DuckDB stays in the corpus query project's
  locked environment rather than becoming an `umol-py` dependency. The Python package tests copy
  the record fixture to temporary directories and cover both physical layouts; Rust-side binding
  tests additionally write temporary molecule-and-reaction datasets under both layouts and read
  both value families through the Python adapter. **Additive (green). Done.** [dep: S5e, S6d, S6e,
  S6f]

  The Python storage adapter is temporary experimental scaffolding, not a surface to expand into a
  corpus or database API. Generic faithful storage remains in scope for `umol-store`, while corpus
  generation, source-specific relations, query tooling, annotations, and notebooks remain outside
  umol proper. Before this document closes, `umol-py::store` must be removed or moved behind a
  separately scoped distribution boundary; doc 205 may use it during the experiment but must not
  add corpus-specific bindings to `umol-py`.

**S6 gate:** one command can generate both corpus layouts at caller-selected locations; DuckDB and
Polars can query the settled facts directly; selected cases reconstruct exactly and render in the
executed notebook, including imported aromatic systems represented by definite constraints.

### S7 — Deferred comparative experiment and evidence

S7 is not the next implementation stage. Doc 205 now uses the S6 path to generate the substantive
atom-mapping corpus, improve automatic evidence, conduct a bounded annotation campaign, and develop
the mapping algorithm. Those uses may restructure the questions and measurements below. S7 resumes
only after doc 205 has supplied real corpus access patterns and must not substitute storage tuning
for that work.

- **S7a — `experimental/atom-mapping/benches/columnar_layouts.rs`: layout comparison.** The
  matrix was narrowed to the question the disposition needs: store 10,000 molecules and 10,000
  atom mappings under each physical layout and report time and space. The population is the
  first 10,000 Rhea reactions with saved RXNMapper mappings under `materials/atom_mapping`
  (6 skipped at the ingestion boundary), ingested through the corpus ingest in 1.5 s; molecule
  records are the mapped-frame lhs and rhs sides, correspondence records are the RXNMapper atom
  mappings induced over each pair. Every cell is one wall-clock pass, the fastest of three when a
  pass takes under 20 s; exact roundtrip is checked on every materialized result; writer
  parameters are 1,024 records per file and 512 per row group. Run with
  `cargo bench --bench columnar_layouts` in the crate. **Additive (green). Done.**
  [dep: S3e, S6b]

  | 10,000 molecules | Nested | Coordinated |
  | --- | ---: | ---: |
  | Parquet on disk | 10.8 MB | 15.3 MB |
  | Arrow in memory | 96.0 MB | 105.6 MB |
  | Arrow construction | 2.59 s | 2.10 s |
  | Parquet write | 2.64 s | 2.35 s |
  | Materialize all records | 2.29 s | 1,082 s |
  | One record by key | 254 ms | 368 ms |
  | Key listing | 9 ms | under 1 ms |

  | 10,000 atom mappings | Nested | Coordinated |
  | --- | ---: | ---: |
  | Parquet on disk | 2.1 MB | 6.4 MB |
  | Arrow in memory | 8.7 MB | 19.5 MB |
  | Arrow construction | 0.19 s | 0.12 s |
  | Parquet write | 0.21 s | 0.21 s |
  | Materialize all records | 0.10 s | 52.7 s |
  | One record by key | 15 ms | 56 ms |

  The current state is unambiguous in favor of the nested layout. Space favors it throughout: 30%
  smaller on disk for molecules, three times smaller for correspondences, and smaller in memory.
  Writing favors the coordinated layout by 10-20% for molecules and is even for correspondences.
  Reconstruction differs by a factor of about 500 in both families because the coordinated reader
  decodes every table into a generic value vector and then scans the complete child tables once
  per record key; that is quadratic in the population and a property of the reader, not of the
  layout, but a linear reader would still pass through the dynamic value layer that the nested
  reader does not need. Reading one record by key is a full decode under both layouts, so neither
  offers a cheap keyed lookup. Key listing is the only operation where the coordinated layout is
  ahead, and the difference is milliseconds. What the measurement does not show is the
  circumstance under which the coordinated layout is useful: analytical queries over its flat
  tables from DuckDB or Polars were not timed, and its case rests entirely on that axis.

  Decision: the nested layout is the only physical layout. The coordinated layout, its selector,
  and the comparison benchmark are removed; the query scripts read the nested record files
  directly; the artifact manifest keeps its layout marker with the constant value `nested`.
  The crate keeps its name and is not split: its public surface is the flat `umol_store::graph`
  module over sealed submodules, the way `umol_graph_ir::ir` is built, with `StoreError` alone at
  the crate root. It exposes the faithful records, the dataset writer and readers, keyed reads and
  key listings, manifest validation, and the decode error; the Arrow batches, leaf codec, schema
  kinds and codebooks, dataset metadata pieces, and the element-literal projection are not exposed
  now, and the projection, the Arrow benchmark, and the batch, leaf, and schema test suites were
  removed with that decision.
- **S7b — workspace verification.** Run formatting, workspace tests, feature-gated property suites,
  the experimental CoordGen feature, clippy, the Python 3.13 build/pytest gate, query-script tests,
  and notebook smoke execution. Record any intentionally non-default native or Python dependency in
  the reproducibility instructions. **Additive (green). Done.** [dep: S7a] Run after the coordinated
  layout was removed: nightly formatting; `cargo test --workspace --all-features --tests` under the
  Python 3.13 venv (77 suites, 29,955 tests); `cargo clippy --workspace --all-targets
  --all-features -- -D warnings`; tests and clippy with all features in both experimental crates;
  `maturin develop` with `store` and `depiction` followed by pytest (1,324 passed, the two doc 195
  skips); the DuckDB and Polars query tests against the regenerated fixture; the notebook
  execution check; `git diff --check`. No non-default native or Python dependency was needed
  beyond the CoordGen feature the notebook already requires.

**S7 gate:** the first experimental slice is complete and reproducible, with evidence rather than
intuition deciding which physical layout and crate surfaces survive.

### S8 — Database-backed molecule interning experiment

S8 asks the storage question that the completed Parquet roundtrip did not answer: can canonical
molecules produced by a generation workflow be interned, deduplicated, reopened, and retrieved at
representative scale, and does the nested Arrow representation provide a useful advantage over a
simple canonical-DSL payload? It does not attempt another atom-mapping review workflow.

The experiment stores only molecules. Reaction correspondences, transformations, provenance,
annotation, generic corpus metadata, Polars object integration, and chemical query operators are out
of scope. The current `Canonicalize::canonical_hash` under one explicit canonicalization context is
the experiment-local unique key. It is sufficient for this rebuildable measurement and makes no
cross-version identity promise.

- **S8a — `umol-store::graph`: public aggregate Arrow batches.** **Done 2026-09-02.** Reverse the
  premature S7 sealing decision for the four aggregate batches: `MoleculeBatch`, `ReactionBatch`,
  `ReactionSpanBatch`, and `MoleculeCorrespondenceBatch`. Each exposes infallible construction from
  its closed faithful records with caller-supplied physical-key offset, checked adoption of an
  arbitrary `RecordBatch`, length and emptiness queries, borrowed batch access, physical-key
  access, and exact record reconstruction. Checked adoption
  accepts only the exact Arrow schema; physical value errors remain lazy `DecodeError` results from
  key or record access. Neither path canonicalizes, repairs, or validates chemistry. Keep leaf
  codecs, builders, schema codebooks, and component-table details sealed. Exact example tests cover
  records -> batch -> records, nonzero key offsets, absent rows, and rejected schemas for all four
  families. This is the only public-API change in S8, and it is reconciled against the faithful
  record contract before completion. **Additive (green).**
  [dep: S1e, S7b]
- **S8b — `experimental/reaction-network::molecule_store`: DuckDB molecule table.** **Done
  2026-09-02.** Add an experiment-local DuckDB adapter whose native table carries a database-local
  id, the canonical hash, and the existing nested molecule fields. Ingest and retrieve
  `MoleculeBatch` values without adding DuckDB as an `umol-store` dependency or an unconditional
  reaction-network dependency. A staging relation handles duplicates within and across batches;
  interning returns the existing or inserted id for every input row. Exact tests cover all-new,
  all-existing, mixed, and within-batch duplicates, close/reopen, point retrieval, and
  reconstruction of the original canonical molecules. Add no public API to a published crate.
  **Additive (green).** [dep: S8a]
- **S8c — `experimental/reaction-network`: generated workload and baselines.** **Done
  2026-09-02.** Add a feature-gated `compare-molecule-storage` runner over the built-in
  normal-polarity oxygen catalog. The natural replay is the canonical seed followed by each
  recorded transformation target in application order. It uses the same graph-and-overlays, VF2,
  Vismara, and Nauty configuration as the network benchmark, with para-stereo disabled, and sends
  the same candidate order through the existing `HashMap<Molecule, _>` interner, the nested
  DuckDB table, and a DuckDB table containing the canonical hash plus canonical molecule DSL.
  Every replayed first-seen id, point-read molecule, and completely materialized molecule is
  checked for exact agreement; inserted counts are checked on the mixed and all-existing passes.
  New ids are explicitly assigned in first-seen order;
  `MoleculeStore::enumerate_molecules` returns exact records in database-local id order. A
  committed two-atom case covers the comparison, while the caller selects the built-in case,
  candidate count, transaction batch size, point-read count, and output directory. The runner
  reports the catalog, configuration, replay definition, generated network shape, representation
  preparation, mixed ingest, an all-existing duplicate replay, reopening, sampled point reads,
  complete materialization, and database file size. **Additive (green).** [dep: S8b]
- **S8d — `discussion/201`: storage disposition evidence.** **Done 2026-09-02.** The measured
  natural workloads were 3,122 candidates and 435 unique molecules for the four-atom case and
  40,375 candidates and 4,537 unique molecules for the five-atom case. A 100,000-candidate run
  cycled the five-atom application-order replay and therefore still contained 4,537 unique
  molecules; it is duplicate-heavy replay evidence, not a claim about 100,000 unique molecules.
  All three stores returned exactly the same ids and molecule values.

  On that 100,000-candidate release run, the in-memory baseline processed mixed ingest at about
  614,000 candidates/s and the all-existing pass at about 698,000 candidates/s. The nested DuckDB
  path processed them at about 12,100 and 13,000 candidates/s, respectively; the DSL-payload
  DuckDB path processed them at about 316,000 and 341,000 candidates/s. One hundred exact point
  reads took 0.580 s from the nested table and 0.023 s from the DSL table. Complete materialization
  of 4,537 molecules took 0.361 s and 0.058 s, respectively. After the two replay passes, the
  nested database occupied 2,371,584 bytes and the DSL database 1,847,296 bytes.

  The nested record codec is not the bottleneck: preparing all 4,537 `MoleculeRecord` values took
  0.0067 s, while rendering their molecule DSL took 0.0279 s. The cost appears after that boundary,
  when DuckDB ingests the wide nested rows and when exported rows are normalized back to the exact
  Arrow schema before reconstruction. The nested representation improved no measured operation.
  The hash-plus-DSL table demonstrates practical persistent deduplication and retrieval, but is
  deliberately equivalent to storing an opaque faithful payload; it demonstrates no
  `umol-store`-specific database capability.

  These network counts reproduce the native normal-polarity census already recorded in doc 207:
  the four-atom closure has 435 molecules and 968 undirected adjacencies rather than the historical
  Python pipeline's 430 and 959, and the five-atom closure has 4,537 and 13,662 rather than 4,513
  and 13,584. Doc 207 records the attribution boundary and a 2026-09-01 control in which the
  historical 22-entry oxygen catalog produced the same updated five-atom split under the current
  binary. The difference is therefore implementation-version evidence, not a rule-catalog change.
  The storage runner reports the actual authoritative network shape instead of silently treating
  the manifest's descriptive reference fields as current generator counts.

  Peak resident memory was not attributed per backend because the combined runner retains the
  generated network and workload across all three measurements; its process peak would not be a
  backend comparison. Runs at one million and ten million repeated candidates were also stopped
  after the natural and 100,000-candidate results established the same stable order-of-magnitude
  difference. The largest unique set measured is 4,537 molecules, so the experiment says nothing
  about ten million unique molecules, database chemical queries, concurrent writers, schema
  evolution, or long-term identity.

  The evidence supports retaining the faithful records and aggregate Arrow batches as a working
  low-level codec through S9, while keeping `umol-store` provisional and making no claim that its
  nested representation is a useful DuckDB schema. Keep the existing Parquet/scan layer unchanged
  through S9: this experiment supports neither removing it nor extending it. The evidence does not
  justify promoting the nested DuckDB adapter or retaining the temporary `umol-py::store` bridge.
  The adapter and comparison runner are experiment apparatus, not a general storage API. A future
  nested database design requires a named structural query that can demonstrate an advantage over
  the opaque-payload baseline; candidate count alone is not that use case. Verification passed the
  default and feature-enabled reaction-network test suites, Clippy with warnings denied, rustdoc
  with warnings denied, formatting, and `git diff --check`. **Additive (green).** [dep: S8c]

**S8 gate:** satisfied. The same canonical molecule workload has exact, reproducible in-memory,
nested-DuckDB, and DSL-payload-DuckDB results through the largest useful duplicate-heavy bound, with
the unique-population and memory limitations stated explicitly. The result is sufficient for the
S9 retention decision without pretending to establish large-unique-corpus behavior.

### S9 — Storage teardown and retained-surface cleanup

S9 applies the disposition established by S7 and S8 by extracting two independent branches from
the same upstream-main commit rather than promoting or tearing down the provisional substrate in
place. One branch contains `umol-coordgen-sys` and the `umol-io`/`umol-py` depiction surface. The
other contains the mapping implementation and the standalone atom-mapping and reaction-network
experiments. Neither branch contains `umol-store`, its Python bridge, or the storage-specific corpus
workflow built to exercise it. The experimental branch retains only the scientific work after its
dependencies on the discarded storage layer and storage-only apparatus are removed.
The faithful-roundtrip design and the comparative measurements remain in this document; they
describe the approach well enough to reconstruct it against a future graph IR without preserving
the crate or its exact experimental Arrow schema. Git history remains the archaeological source for
exact codebooks and field layouts. No storage stub or speculative schema-evolution follow-up is
retained.

The depiction branch is independent of that result. `umol-coordgen-sys`, the `umol-io` layout,
depiction, and SVG path, and the explicit Python SVG display boundary remain. S9 cleans and verifies
that retained surface after removing the storage consumers.

- **S9a — `discussion/201`: final storage disposition.** **Done 2026-09-02.** Record the decision to
  remove `umol-store` rather than retain it as a stub. Preserve the structural-roundtrip contract,
  codec approach, S7 layout comparison, S8 DuckDB comparison, and their measurement limits here.
  State the reusable conclusion: faithful graph-IR records and nested Arrow encoding are feasible
  and the codec is fast, but neither the Parquet scans nor nested DuckDB representation demonstrated
  a useful storage or query capability beyond an opaque canonical payload. Future storage work must
  begin with a concrete database-side operation that benefits from structural projection. This
  subitem changes documentation only and adds no public API. **Additive (green).** [dep: S7b, S8d]
- **S9b — `experimental/reaction-network`: remove the S8 apparatus.** **Done 2026-09-02 on the
  experimental branch.** Delete the DuckDB molecule
  adapter, comparison runner, `molecule-store` feature, and their optional DuckDB and `umol-store`
  dependencies. Preserve the reaction-network generator, rule catalogs, native census, QRS
  producer, and their existing tests and benchmarks. The reaction-network default and property
  suites and Clippy restore the standalone crate to green. This removes experimental code and adds
  no public API. **Breaking (red -> green).** [dep: S9a]
- **S9c — `experimental/atom-mapping`: detach the storage layer while preserving mapping work.**
  **Done 2026-09-02 on the experimental branch.**
  Keep the experimental crate and its source ingestion, mapping inputs, alignments, mapping and QRS
  evidence, evaluation, spot checks, rule/network integration, and canonicalization fixtures,
  tests, and benchmarks. Remove the store-backed dataset reader/writer, generated query fixture,
  DuckDB and Polars query project, inspection notebook and its launch/environment helpers, and
  binaries or tests whose purpose is only to produce or consume that persisted corpus. Migrate the
  surviving in-memory aggregate from `MoleculeRecord` and `MoleculeCorrespondenceRecord` to the
  existing graph-IR `Molecule` and `MoleculeCorrespondence` values. `Molecule` remains a closed
  integrity-valid aggregate; `MoleculeCorrespondence` remains an open carrier whose agreement with
  the referenced frames is established by the existing ingest operation that first combines them.
  `CorpusIngest` remains operation-populated with no public parts constructor, performs no implicit
  canonicalization or chemistry validation, and retains its existing `IngestError`, ordinary
  absence, and no-panic contracts. No Python boundary replaces it. Rename storage-specific keys or
  accessors only where they continue to describe in-memory scientific values. Retain the
  selected-EDN and provenance work of the canonicalization corpus builder while removing its
  unconsumed full-population Parquet output and store-relative keys. The notebook's usability or any
  replacement annotation interface is a separate decision, not part of this teardown. This removes
  storage-only experimental surfaces, preserves the scientific public surface, and adds no public
  API. The atom-mapping tests, property suite, and retained benchmarks restore the standalone crate
  to green. **Breaking (red -> green).** [dep: S9a]
- **S9d — `umol-py`: remove the temporary storage bridge.** **Done 2026-09-02 by branch
  extraction.** Remove the `store` feature,
  `umol-store` dependency, `StoreDataset` module and export, Python wrapper plumbing, and storage
  tests. Retain the independent `depiction` feature, `MoleculeLayoutAlgorithm`, `Svg`, and the
  explicit molecule and reaction `depict_with` methods. Reconcile the remaining Python export list
  and feature matrix, then rebuild under the repository's Python 3.13 environment and run the
  focused import and depiction tests. This removes the provisional `StoreDataset` public class and
  adds no public API. **Breaking (red -> green).** [dep: S9a, S9c]
- **S9e — workspace and repository surface: remove `umol-store`.** **Done 2026-09-02 by branch
  extraction.** After all consumers are green,
  delete the crate and remove it from the workspace, dependency metadata, and lockfiles. Remove
  generated store artifacts and update `AGENTS.md`, current development guides, doc 205's open
  storage assumptions, and the status index to describe the remaining repository accurately.
  Preserve completed docs 207 and 216 as historical records; add a dated clarification only where
  the removal would otherwise make a current claim ambiguous. This removes the crate's complete
  public API and adds no replacement. **Breaking (red -> green).** [dep: S9b, S9c, S9d]
- **S9f — `umol-io`, `umol-coordgen-sys`, and `umol-py`: retain and reconcile depiction.** **Done
  2026-09-02.** The extracted branch preserves explicit CoordGen selection and the separation
  between layout, format-neutral depiction, SVG rendering, and Python rich display. The public
  surface matches the S4-S5 contract: `MoleculeLayout` has checked construction and frame
  agreement; `Depiction` remains operation-issued with no public aggregate constructor; rendering
  consumes that closed scene; and Python exposes only the explicit layout selector, frozen `Svg`
  value, and molecule and reaction `depict_with` operations. The vendored native boundary is
  feature-gated and reports input, allocation, backend, and non-finite-output failures without
  exposing raw pointers. No implicit algorithm default, configuration layer, or additional public
  seam was introduced by the split. A later visual-quality pass is a separate product process, not
  unfinished scope of this storage experiment. **Additive (green).** [dep: S9d]
- **S9g — repository verification and closeout.** **Done 2026-09-02.** Both extracted branches are
  rooted directly at upstream main commit `4053d820054302c13c7598f619c6ab810c063a3f` and have no
  storage dependency. The depiction branch contains no experimental or `umol-store` path; the
  experimental branch contains no CoordGen, depiction, DuckDB, Parquet, notebook, or `umol-store`
  dependency. Its two crates are standalone `publish = false` workspaces excluded from the root
  workspace, and Cargo lockfiles are no longer generally ignored.

  Post-split verification passed the default and `proptest` suites plus strict all-target Clippy in
  both experimental crates; the CoordGen ABI tests; `umol-io` with CoordGen, including the layout
  suite; Rust tests and strict Clippy for the depiction-enabled Python binding; the Python 3.13
  pytest suite with depiction enabled; local packaging of `umol-coordgen-sys`; formatting; and
  `git diff --check`. The original combined branch remains archival working state rather than a
  merge dependency. The first crates.io publication of `umol-coordgen-sys` is a separate release
  operation required before trusted publishing can manage later releases; it is not a completion
  condition for this experiment. **Additive (green).** [dep: S9e, S9f]

**S9 gate:** satisfied. `umol-store` and every storage-only consumer are absent, both experimental
crates and their retained scientific work are green, the retained depiction and Python SVG
surfaces are reconciled and verified, and this document preserves the useful approach and negative
evidence without leaving an unowned storage stub.

### Dependency summary

The storage branch is `S0 -> S1 -> S2 -> S3`, after which the atom-mapping path reaches S6 and
continues in doc 205. S7 records the first layout evidence; S8 reuses the molecule record and Arrow
path in a database-backed generation and deduplication experiment; both terminate in the S9 branch
split. The depiction branch receives its retained-surface reconciliation in S9f. The independent
experimental branch preserves the mapping and network work without Store or depiction. The
reaction-template counterexample is deferred rather than a dependency of this experiment.

Direct graph-IR-to-Arrow builders, direct Arrow-row-to-IR readers, correspondence-aware reaction
alignment, external layout comparisons, row codecs, borrowed views, materialized chemical
projections, and native database chemical operators are not scheduled stages. They return only if
future concrete evidence justifies them; any reconstruction must preserve faithful graph-IR
roundtrip but need not reproduce this experiment's exact records or Arrow schema.
