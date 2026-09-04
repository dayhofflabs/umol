# 212 — The remapping id-transport layer

Status: In Progress
Date: 2026-08-26
Relates: [211](211-relation-frames-and-api-2026-08-26.md),
[214](214-aggregate-frame-semantics-2026-08-28.md),
[nomenclature guide](../docs/development/nomenclature.md)

## Purpose

Doc [211](211-relation-frames-and-api-2026-08-26.md) completes the compaction row of the
id-transport types because relation-set `compact` returns a relation compaction. That work exposed
the same missing layer in the remapping row. Remapping is not on 211's critical path and its
extraction turns on a representation question that compaction does not raise, so it is recorded
separately.

This document records the settled representation, public surface, and audit of existing
construction sites. Implementation is in progress.

## Finding

There are three id-transport concepts. Each has three layers: one id space, the graph's node and
edge spaces, and the molecule's eight entity families.

| concept | single id space | graph, node and edge | molecule, eight families |
| --- | --- | --- | --- |
| partial bijection | `Correspondence<Id>` | `GraphCorrespondence` | `MoleculeCorrespondence` |
| removal and dense shift | `Compaction<Id>` | `GraphCompaction` | `MoleculeCompaction` |
| total relabel | `Remapping<Id>` | `GraphRemapping` | `MoleculeRemapping` |

The aggregate names state the layer they cover. The legacy sparse `IdRemapping` remains only until
its consumers are classified and migrated later in this plan.

## How remapping differs from compaction

Doc 211 extracts the compaction layer because the arithmetic is duplicated: `compact_relation`,
`uncompact_dense`, and `normalize_removed` in `ir/compact.rs` are line-for-line reimplementations of
`Compaction::compact_node`, `Compaction::uncompact_node`, and the sort-and-dedup in
`Compaction::new`. That extraction is mechanical and deletes code.

Remapping has no such duplication. The two layers use different representations for the same
concept.

```rust
pub struct Remapping {                  // graph.rs:643
    nodes: Vec<NodeId>,                 // dense, indexed by source id
    edges: Vec<EdgeId>,
}
// try_map_node(old) = self.nodes.get(old.0 as usize).copied()

pub struct IdRemapping {                // ir/remap.rs:235
    atom: HashMap<AtomId, AtomId>,      // sparse
    bond: HashMap<BondId, BondId>,
    dative: HashMap<DativeBondId, DativeBondId>,
    // five more
}
```

Extracting a shared `Remapping<Id>` therefore requires choosing a representation first. That is a
design decision, not an extraction.

## Representation decision

All three remapping layers use a dense positional representation. The image vector defines the
source domain: source id `i` maps to the value at position `i`, so construction establishes that
every source id has an image. The image need not be dense, injective, or surjective; a remapping may
embed its dense source into a sparse or larger ambient id space.

The vector is stored as `index_vec::IndexVec<Id, Id>`. The existing id newtypes implement the
crate's `Idx` indexing contract, which confines conversion to the machine index required by vector
storage. They do not acquire general `Into<usize>` conversions. This preserves the distinction
between an id as a typed table index and an arbitrary integer conversion.

The extracted single-id-space layer is therefore `Remapping<Id>`, graph core holds one such value
for nodes and one for edges, and the molecule layer holds one for each of its eight entity kinds.
The current `IdRemapping` hash maps do not establish this contract. They have no declared source
domain, so their effective guarantee is only that every lookup made by a particular consumer happens
to have an entry.

## Public contract

The three public layers are:

- `Remapping<Id>`, containing one typed dense image vector;
- `GraphRemapping`, containing node and edge remappings; and
- `MoleculeRemapping`, containing a `GraphRemapping` and remappings for dative bonds, aromatic
  systems, multicenter bonds, noncovalent interactions, stereo atoms, and stereo bonds.

Construction supplies the complete image vector for every represented source id space. The
aggregate types expose their component remappings. `Remapping<Id>::map` returns the image of a
source id and asserts that the id is in the declared source domain; `try_map` returns `None` for an
out-of-domain source id. The aggregate types provide the corresponding typed pairs, such as
`map_node` and `try_map_node`, for every component id space.

`GraphCorrespondence::to_remapping` and `MoleculeCorrespondence::to_remapping` continue to return
`Option`: a correspondence has a remapping exactly when it is total on its dense left-hand id
spaces. Their successful results construct the corresponding dense aggregate directly.

## Operation naming and failure boundary

`map_*` names lookup of one id. `remap` names transport of a complete reference-bearing value. A
public whole-value operation that accepts an independently supplied remapping provides both forms:

- `try_remap` returns `None` when any referenced id lies outside the remapping's source domain; and
- `remap` asserts complete coverage and documents that requirement.

The paired whole-value surface applies to `Constraint`, `RelationalConstraint`,
`MoleculeConstraint`, and `Delta`. The current free `remap_delta` operation becomes the inherent
`Delta::remap` and `Delta::try_remap` pair. Existing molecule and reaction-span transport remains
paired in the same way even where the transport carrier is a correspondence rather than a
remapping.

## Adjacent naming corrections

Two existing operations do not have remapping semantics and should be renamed as part of this
cleanup:

- `MoleculeMetadata::remap` retains only metadata whose referenced entities are present in a
  correspondence. It is projection, so the Rust and Python operation becomes `project`.
- `Sgroup::remap_indices` translates external row indices through an independently supplied lookup
  and returns `None` on missing coverage. It becomes `try_reindex`.

## Construction-site findings

### Correspondence conversion

The two `to_remapping` operations are the same operation at two layers, and both are built from a
correspondence that is total on the left over a dense left id space.

`GraphCorrespondence::to_remapping` (`correspondence.rs:365`) checks `is_total_on_left`, then
collects only the right id of each matched pair into a positional vector. That is valid precisely
because the left space is dense and `matched_pairs` is in left order.

`MoleculeCorrespondence::to_remapping` (`ir/correspondence.rs:246`) performs the identical check and
then collects the whole pair into a hash map, discarding the density it just established.

Both conversions therefore become direct construction of dense image vectors.

### Constraint edits

`ConstraintEdit` currently uses `IdRemapping` for two different operations.

`ConstraintEdit::new` collects only the entity ids referenced by the supplied constraint, maps them
to normalized indices in private per-kind handle vectors, and rewrites the constraint to those
indices. Its source is the finite set of references encountered in that constraint, not a declared
dense id space. This operation is not a remapping and must no longer be represented by the remapping
type.

`ConstraintEdit::resolve` maps every dense private handle index to the entity id obtained by
resolving that handle. Its source domain is exactly each private vector's `0..len`, with one image
for every source id. This is a remapping and admits the dense representation directly.

### Molecule splitting

`Molecule::split` returns each compact component together with a `MoleculeCorrespondence` from the
component to the original molecule. The current implementation inverts the matched pairs into a
sparse original-to-component `IdRemapping` solely to rewrite constraints routed to that component.
That mapping is partial on the original molecule and is therefore not a remapping. Constraint
references should instead be transported through the existing correspondence in the reverse
direction. The operation remains a split producing a correspondence; no compaction result is added.

### Reaction reversal

`Reaction::reverse` currently materializes `self.to_reaction_span()?`, but uses the span only to
obtain its right projection as the new left-hand molecule. It separately collects removed and
created ids, builds sparse maps over every id the inverted deltas will reference, and uses those
maps to re-anchor the inverted deltas. This works because the hash-map representation checks
coverage only when each reference is looked up. In particular, ids introduced by `Add` deltas need
not form a dense source space, so these maps do not satisfy the remapping contract.

The selected replacement uses the reaction span as the complete intermediate:

```text
Reaction -> ReactionSpan -> reverse sides -> Reaction
```

Reversing the span exchanges `Added` with `Removed`, exchanges the left and right values of
`Modified`, and leaves `Unchanged` unchanged; constraint spans reverse in the same way. The span
remains in its dense union id space. `ReactionSpan::to_reaction` then reanchors that union to the
new left-hand side by assigning its present entities the dense prefix and its additions the
following ids. Its union-to-reaction map is therefore a genuine dense remapping.

This is an implementation decision for the existing `Reaction::reverse` operation and does not
change its public API. It replaces the current sparse delta-remapping path rather than adding that
path to the remapping abstraction.

## Scope boundary

This document owns the remapping row only. The compaction row, the relation-set surface, frame
transport, and the `Reframe` operation belong to doc
[211](211-relation-frames-and-api-2026-08-26.md) and are not reopened here. Nothing in doc 211
depends on this document, and doc 211 may be implemented and closed while this remains proposed.

Doc 211 has landed the compaction row and supplies the naming and layering precedent. The existing
construction sites have now been classified: only genuine total relabelings move to the dense
layer; operation-local bookkeeping and correspondence transport do not migrate merely because they
currently use `IdRemapping`.

## Public-symbol inventory

The implementation adds:

- `Compaction<Id>::compact_vec` for applying one id-space compaction to its aligned data column;
- `Remapping<Id>` with `new`, `map`, and `try_map`;
- `index_vec::Idx` implementations for the existing graph-core and graph-IR id newtypes;
- `GraphRemapping` with `new`, `nodes`, `edges`, and the typed node/edge `map_*` and `try_map_*`
  methods;
- `MoleculeRemapping` with `new`, `graph`, six overlay-component accessors, and typed `map_*` and
  `try_map_*` methods for all eight molecule entity kinds; and
- `try_remap` alongside `remap` on `Constraint`, `RelationalConstraint`, `MoleculeConstraint`, and
  `Delta`.

The typed molecule methods use the complete entity names: `map_dative_bond`,
`map_aromatic_system`, `map_multicenter_bond`, `map_noncovalent_bond`, `map_stereo_atom`, and
`map_stereo_bond`, with matching `try_map_*` names. Atom and bond lookup remains `map_atom` and
`map_bond`.

The implementation changes or retires:

- the graph-specific free functions `compact_node_vec` and `compact_edge_vec`, replaced by the
  generic `Compaction<Id>` method;
- the current graph aggregate `Remapping`, renamed to `GraphRemapping` so the generic type can own
  the unqualified name;
- `IdRemapping`, replaced by `MoleculeRemapping`;
- the return types of `GraphCorrespondence::to_remapping` and
  `MoleculeCorrespondence::to_remapping`;
- free `remap_delta`, replaced by `Delta::remap` and `Delta::try_remap`;
- Rust and Python `MoleculeMetadata::remap`, renamed to `project`; and
- `Sgroup::remap_indices`, renamed to `try_reindex`.

`Reaction::reverse`, `Molecule::{remap, try_remap}`, and
`ReactionSpan::{remap, try_remap}` retain their current public signatures. No remapping carrier is
added to the Python API.

## Implementation plan

### S0 — Establish the performance baseline

#### S0a — Benchmark reaction reversal

**Module:** `umol-graph-ir/benches/reaction.rs`.

Add one top-level `Reaction::reverse` benchmark using a reaction with additions, removals,
modifications, overlays, and constraints, so it exercises the current sparse re-anchoring path.
Record the baseline result before changing the implementation and retain the same case for the final
comparison.

**Tests and evidence:** Run the existing reaction reversal tests and the filtered reaction benchmark;
record the benchmark result in this document.

**Change class:** additive evidence (green).

**Dependencies:** none.

### S1 — Extract the graph-core remapping row

#### S1a — Add `Remapping<Id>` and re-express the graph aggregate **Done**

**Module:** `umol-graph-core/src/remap.rs`, `compact.rs`, `graph.rs`, `relation.rs`,
`correspondence.rs`, `lib.rs`, and their unit tests.

Add the dense `Remapping<Id>` image vector and its `new`, `map`, and `try_map` methods, backed by
`index_vec::IndexVec<Id, Id>`. Implement `index_vec::Idx` for `NodeId`, `EdgeId`, and `RelationId`;
do not add general id-to-`usize` conversions. Rename the existing node/edge aggregate to
`GraphRemapping`, store one generic remapping for each id space, and preserve its typed lookup
methods while adding `nodes` and `edges` accessors. Update graph-core relation transport and
`GraphCorrespondence::to_remapping` to the new aggregate. Name the private operation modules
`remap` and `compact` consistently. Move aligned-column application from the graph-specific free
functions to `Compaction<Id>::compact_vec`.

**Tests and evidence:** Table tests cover an empty image vector, first and last images, a sparse
target image, a repeated target image, checked out-of-domain lookup, and the asserted lookup panic.
Retain the existing node/edge aggregate and correspondence-conversion cases under
`GraphRemapping`. Exact `compact_vec` cases cover identity, removed positions, scattered removals,
and removals outside the supplied column. Run `cargo test -p umol-graph-core`.

**Change class:** breaking rename with an additive generic (red until S1b).

**Dependencies:** [dep: S0a]

The generic remapping is backed by `index_vec::IndexVec`, and the graph-core id newtypes implement
its `Idx` contract without gaining general integer conversions. `cargo test -p umol-graph-core`
passes.

#### S1b — Migrate graph-IR graph remapping consumers **Done**

**Module:** graph-remapping imports and call sites in `umol-graph-ir/src/ir`, including the relation
families, ligand transport, molecule remapping, pushout, and reaction spans.

Replace uses of the former graph aggregate name with `GraphRemapping`. This is a type-name migration
only; participant and ligand transport semantics remain unchanged.

**Tests and evidence:** Existing graph-IR molecule, relation, pushout, and reaction-span tests pass
unchanged. Run `cargo test -p umol-graph-ir`.

**Change class:** caller migration (restores green).

**Dependencies:** [dep: S1a]

All graph-IR consumers use `GraphRemapping`. `cargo test -p umol-graph-ir` passes.

### S2 — Add the molecule remapping aggregate

#### S2a — Add `MoleculeRemapping` **Done**

**Module:** `umol-graph-ir/src/ir/remap.rs`, `compact.rs`, `id.rs`, `ir.rs`, and unit tests.

Add `MoleculeRemapping` alongside the existing `IdRemapping`. It contains a `GraphRemapping` plus
six `Remapping<..Id>` overlay components. Its constructor takes the graph aggregate and six dense
image vectors; its component accessors and typed `map_*`/`try_map_*` pairs expose the contract
listed above. Implement `index_vec::Idx` for the eight molecule id types. Preserve the old
aggregate's empty `Default` behavior for the replacement type. Keep molecule compaction and
remapping in the separate private `compact` and `remap` modules.

**Tests and evidence:** Exact table cases exercise successful and out-of-domain lookup for every
entity kind, component access, sparse and repeated target images, and the asserted failure boundary.
Run `cargo test -p umol-graph-ir ir::compact` and `cargo test -p umol-graph-ir ir::remap`.

**Change class:** additive (green).

**Dependencies:** [dep: S1b]

`MoleculeRemapping` now carries the graph remapping and six typed dense overlay remappings, with
checked and asserted lookup for every entity kind. All eight graph-IR id newtypes implement
`index_vec::Idx`; no general id-to-`usize` conversion was added. `cargo test -p umol-graph-ir
ir::compact` passes all 31 selected cases and `cargo test -p umol-graph-ir ir::remap` passes all 25
selected cases.

### S3 — Remove uses that are not remappings

#### S3a — Keep constraint-edit handle substitution local

**Module:** `umol-graph-ir/src/ir/edit.rs` and its `ConstraintEdit` tests.

Replace the `IdRemapping` built by `ConstraintEdit::new` with direct substitution from each
referenced entity id to the normalized slot for its resolved handle. Keep those sparse per-kind
lookups local to `ConstraintEdit`; do not introduce a transport carrier or a new public helper.
`ConstraintEdit::resolve` remains a true dense-remapping consumer and is not changed in this
subitem.

**Tests and evidence:** Retain exact construction and resolution cases. Add cases for sparse source
ids, repeated use of one referenced entity, a missing handle, and a handle-kind mismatch. Run the
focused `ConstraintEdit` tests.

**Change class:** internal representation change (green).

**Dependencies:** [dep: S2a]

#### S3b — Transport split constraints through the correspondence

**Module:** `umol-graph-ir/src/ir/molecule.rs` and split tests.

Delete the partial original-to-component `IdRemapping`. For each constraint already routed to a
component, translate its references through the reverse direction of the component-to-original
`MoleculeCorrespondence`. Use the existing `ConstraintEdit` handle boundary for the translation;
do not add a second public constraint-transport surface. Preserve the returned component molecule
and component-to-original correspondence exactly.

**Tests and evidence:** Exact split cases cover constraints referring to atoms, bonds, every overlay
kind, relational constraints, and molecule constraints; assert both the complete component and the
returned correspondence. Run the focused molecule split tests.

**Change class:** internal transport correction (green).

**Dependencies:** [dep: S3a]

#### S3c — Reverse reactions through the reaction span

**Module:** `umol-graph-ir/src/ir/reaction_span.rs` and its reaction reversal tests.

Replace the sparse reverse bookkeeping with the settled path
`Reaction -> ReactionSpan -> reverse sides -> Reaction`. Side reversal is private implementation:
swap `Added`/`Removed`, swap the values of `Modified`, preserve `Unchanged`, and apply the same side
swap to constraint spans. Delete `reversed_remapping`; do not add a public reaction-span reversal
method.

**Tests and evidence:** Exact cases cover every entity family and every span state, including
constraints and stereo frames. Retain the reverse-twice and span-side equivalence assertions. Rerun
the S0a benchmark and record the comparison.

**Change class:** internal algorithm replacement with unchanged public API (green).

**Dependencies:** [dep: S3b]

### S4 — Move true remapping consumers onto the dense layer

The public signature changes in this stage may leave the workspace red between subitems. The stage
ends with every producer migrated and `IdRemapping` removed.

#### S4a — Add checked constraint transport

**Module:** `umol-graph-ir/src/ir/constraint/molecule.rs`, `relational.rs`, the leaf constraint
modules, and their unit tests.

Change constraint remapping to accept `MoleculeRemapping`. Add `try_remap` to `Constraint`,
`RelationalConstraint`, and `MoleculeConstraint`; it returns `None` when any referenced id is outside
the corresponding source domain. Make `remap` the asserted counterpart with a documented coverage
requirement. Leaf forms that contain no entity references remain infallible and do not acquire a
meaningless checked variant.

**Tests and evidence:** Exhaustive table cases cover the id-bearing variants and recursive
`And`/`Or`/`Not`, with exact successful values and missing coverage in each id family. Keep each
no-reference variant as an exact identity case. Run the focused constraint tests once the stage is
green.

**Change class:** breaking public signature and additive checked operations (red until S4d).

**Dependencies:** [dep: S2a, S3c]

#### S4b — Move delta transport onto `Delta`

**Module:** `umol-graph-ir/src/ir/delta.rs` and its unit tests.

Replace free `remap_delta` with inherent `Delta::try_remap` and `Delta::remap` over
`MoleculeRemapping`. The checked method declines on the first uncovered entity or participant; the
asserted method documents and asserts complete coverage. Preserve participant sequence and stereo
frame order.

**Tests and evidence:** Exhaustive table cases cover all delta families, entity ids, participants,
stereo sites and ligands, and nested constraint deltas. Assert exact `None` results for missing
coverage and retain the inverse-remapping roundtrip as a targeted unit test. Do not add a property
test for this finite variant surface.

**Change class:** breaking free-to-inherent API migration (red until S4d).

**Dependencies:** [dep: S4a]

#### S4c — Convert correspondence and molecule producers

**Module:** `umol-graph-ir/src/ir/correspondence.rs`, `edit.rs`, `molecule.rs`,
`molecule/remapping.rs`, `molecule/pushout.rs`, and their unit tests.

Change `MoleculeCorrespondence::to_remapping` to construct dense image vectors and return
`MoleculeRemapping`. Convert the true remapping producers in `ConstraintEdit::resolve`, molecule
combination, molecule remapping, and pushout from hash maps to dense vectors. Update their
constraint transport to the checked or asserted operation justified by each producer's coverage.
Delete the now-unused offset hash-map helper.

**Tests and evidence:** Retain correspondence totality, combination, molecule remapping, edit
resolution, and pushout cases. Add exact coverage assertions where a producer establishes the dense
source domains. Run the focused tests for those modules once the stage is green.

**Change class:** breaking producer and caller migration (red until S4d).

**Dependencies:** [dep: S4b]

#### S4d — Convert reaction-span producers and retire `IdRemapping`

**Module:** `umol-graph-ir/src/ir/reaction_span.rs`, `remap.rs`, `ir.rs`, and their unit tests.

Convert the genuine remappings used by superimposition, side projection, and `to_reaction` to dense
image vectors and `MoleculeRemapping`. Replace direct hash-map indexing with the typed asserted
lookups whose domains those producers establish. Remove `IdRemapping`, its hash-map imports, and all
remaining old-name references.

**Tests and evidence:** Retain exact superimposition, lhs/rhs projection, span-to-reaction,
roundtrip, and invalid-reference cases. Search the workspace for `IdRemapping` and `remap_delta`;
only historical discussion text may remain. Run `cargo test -p umol-graph-ir` and rerun both the S0a
reaction benchmark and the existing molecule-remapping benchmark.

**Change class:** final caller migration and type retirement (restores green).

**Dependencies:** [dep: S4c]

### S5 — Correct adjacent operation names

This stage is deferrable relative to the dense-remapping deliverable, but required before this
document can close.

#### S5a — Rename metadata projection in Rust and Python

**Module:** `umol-graph-ir/src/dsl/metadata.rs`, `umol-py/src/metadata.rs`, Rust tests, and
`umol-py/tests/test_metadata.py`.

Rename `MoleculeMetadata::remap` to `project` on both language surfaces and migrate all callers. The
operation continues to retain matched keywords, omit unmatched keywords, and preserve atom aliases.

**Tests and evidence:** Exact Rust and Python cases cover a total correspondence, unmatched keyword
omission, alias preservation, and composition. For Python verification, activate `umol-py/.venv`,
confirm Python 3.13, rebuild with `maturin develop`, and run the focused metadata tests.

**Change class:** breaking rename with caller migration (red then green within the subitem).

**Dependencies:** [dep: S4d]

#### S5b — Rename S-group index translation

**Module:** `umol-io/src/table_ir/sgroup.rs`, `ctfile_data.rs`, `cx_data.rs`, and their unit tests.

Rename `Sgroup::remap_indices` to `try_reindex` and migrate both TableIR callers. Preserve its
checked lookup behavior and the callers' existing omission of an S-group whose required indices
cannot be translated.

**Tests and evidence:** Exact cases cover every index-bearing S-group field and missing atom and
bond coverage, including connecting bonds. Run `cargo test -p umol-io`.

**Change class:** breaking rename with caller migration (red then green within the subitem).

**Dependencies:** [dep: S5a]

### S6 — Align documentation and close

#### S6a — Update the guides and complete verification

**Module:** `docs/development/nomenclature.md`, `docs/development/data-types.md`, crate rustdoc,
this document, and `discussion/000-status.md`.

Record the three-layer remapping vocabulary, dense source-domain contract, `map`/`remap` distinction,
and checked/asserted naming rule. Replace stale names and examples. Reconcile every changed public
symbol against the inventory above, record the benchmark comparison and verification results, then
mark this document Completed only when the implementation and all required checks have landed.

**Tests and evidence:** Run `cargo +nightly fmt --all`, `cargo test --workspace`, and
`cargo clippy --workspace --all-targets -- -D warnings`. Run the Python gate from S5a. Search code
and current guides for every retired symbol, and run `git diff --check`.

**Change class:** documentation migration and verification (green).

**Dependencies:** [dep: S5b]

### Dependency summary

The critical path is
`S0a -> S1a -> S1b -> S2a -> S3a -> S3b -> S3c -> S4a -> S4b -> S4c -> S4d -> S5a -> S5b -> S6a`.

The dense-remapping deliverable is green at S4d. S5 is independently deferrable if only that core
deliverable is required, but the adjacent public names remain part of this document's accepted
scope and the document remains Proposed until S5 and S6 are complete.
