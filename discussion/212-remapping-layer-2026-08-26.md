# 212 — The remapping id-transport layer

Status: Completed
Date: 2026-08-26
Relates: [211](211-relation-frames-and-api-2026-08-26.md),
[214](214-aggregate-frame-semantics-2026-08-28.md),
[218](218-mutation-witness-2026-08-31.md),
[nomenclature guide](../docs/development/nomenclature.md)

## Purpose

Doc [211](211-relation-frames-and-api-2026-08-26.md) completes the compaction row of the
id-transport types because relation-set `compact` returns a relation compaction. That work exposed
the same missing layer in the remapping row. Remapping is not on 211's critical path and its
extraction turns on a representation question that compaction does not raise, so it is recorded
separately.

This document records the dense typed carrier extracted from that audit, the removal of two
non-remapping uses, the reaction-reversal cleanup, and the remapping-module name alignment. The
wider operation-witness work continues in doc 218.

## Finding

There are three id-transport concepts. Each has three layers: one id space, the graph's node and
edge spaces, and the molecule's eight entity families.

| concept | single id space | graph, node and edge | molecule, eight families |
| --- | --- | --- | --- |
| partial bijection | `Correspondence<Id>` | `GraphCorrespondence` | `MoleculeCorrespondence` |
| removal and dense shift | `Compaction<Id>` | `GraphCompaction` | `MoleculeCompaction` |
| dense total map | `Remapping<Id>` | `GraphRemapping` | `MoleculeRemapping` |

The aggregate names state the layer they cover. The legacy sparse `IdRemapping` remains because
migrating its consumers requires the wider operation-witness design now assigned to doc 218 on the
`feature/mutation-witness` branch.

## Semantic classification

For operation witnesses, the intended concepts are narrower than the current carriers alone
express:

- a correspondence is a partially bijective mapping with explicit left and right counts;
- a remapping is the efficient special case of a total bijection with a dense left-hand side; and
- a compaction is the order-preserving special case with no unmatched right-hand elements, so the
  right-hand side is no larger than the left-hand side.

The special carriers should be used only when those properties are part of an operation's
semantics. Otherwise, correspondence is the default witness. This distinction matters for witness
direction and composition, not merely for storage.

The `Remapping<Id>` implemented in this work establishes a dense, total source domain, but its
constructor currently permits sparse or repeated target images. It therefore implements the
representation needed by a remapping without yet enforcing the full bijection contract above.
Tightening that contract cannot be done by adding a repeated-image check alone: existing producers
first need to be classified by operation semantics, and multi-output operations need a target id
space capable of distinguishing their outputs. That semantic work is transferred to doc 218.

## How remapping differs from compaction

Doc 211 extracts the compaction layer because the arithmetic is duplicated: `compact_relation`,
`uncompact_dense`, and `normalize_removed` in `ir/compact.rs` are line-for-line reimplementations of
`Compaction::compact_node`, `Compaction::uncompact_node`, and the sort-and-dedup in
`Compaction::new`. That extraction is mechanical and deletes code.

Before S1a, remapping had no such duplication. The two layers used different representations for
the same concept.

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

The extracted carriers use a dense positional representation. The image vector defines the source
domain: source id `i` maps to the value at position `i`, so construction establishes that every
source id has an image.

The vector is stored as `index_vec::IndexVec<Id, Id>`. The existing id newtypes implement the
crate's `Idx` indexing contract, which confines conversion to the machine index required by vector
storage. They do not acquire general `Into<usize>` conversions. This preserves the distinction
between an id as a typed table index and an arbitrary integer conversion.

The extracted single-id-space layer is therefore `Remapping<Id>`, graph core holds one such value
for nodes and one for edges, and the molecule layer holds one for each of its eight entity kinds.
The current `IdRemapping` hash maps do not establish a declared source domain, so their effective
guarantee is only that every lookup made by a particular consumer happens to have an entry.

## Public contract

The public additions completed here are:

- `Remapping<Id>`, containing one typed dense image vector;
- `GraphRemapping`, containing node and edge remappings; and
- `MoleculeRemapping`, containing a `GraphRemapping` and remappings for dative bonds, aromatic
  systems, multicenter bonds, noncovalent interactions, stereo atoms, and stereo bonds.

Construction supplies the complete image vector for every represented source id space. The
aggregate types expose their component remappings. `Remapping<Id>::map` returns the image of a
source id and asserts that the id is in the declared source domain; `try_map` returns `None` for an
out-of-domain source id. The aggregate types provide the corresponding typed pairs, such as
`map_node` and `try_map_node`, for every component id space.

`GraphCorrespondence::to_remapping` now constructs `GraphRemapping` directly when the
correspondence is total on its dense left-hand node and edge spaces. The molecule-level conversion
still returns `IdRemapping`; its migration is part of the deferred semantic work.

## Construction-site findings

### Correspondence conversion

`GraphCorrespondence::to_remapping` (`correspondence.rs:365`) checks `is_total_on_left`, then
collects only the right id of each matched pair into a positional vector. That is valid precisely
because the left space is dense and `matched_pairs` is in left order.

`MoleculeCorrespondence::to_remapping` (`ir/correspondence.rs:431`) performs the identical check and
then collects the whole pair into a hash map, discarding the density it just established.

Only the graph-core conversion moved in this work. Whether and how the molecule conversion should
produce a remapping is coupled to the witness-direction and composition questions assigned to doc
218.

### Constraint edits

Before S3a, `ConstraintEdit` used `IdRemapping` for two different operations.

`ConstraintEdit::new` collected only the entity ids referenced by the supplied constraint, mapped
them to normalized indices in private per-kind handle vectors, and rewrote the constraint to those
indices. Its source was the finite set of references encountered in that constraint, not a declared
dense id space. S3a replaced the remapping carrier with a local sparse entity map.

`ConstraintEdit::resolve` maps every dense private handle index to the entity id obtained by
resolving that handle. Its source domain is exactly each private vector's `0..len`, with one image
for every source id. It remains an `IdRemapping` consumer until the operation-witness migration.

### Molecule splitting

`Molecule::split` returns each compact component together with a `MoleculeCorrespondence` from the
component to the original molecule. Before S3b, the implementation inverted the matched pairs into
a sparse original-to-component `IdRemapping` solely to rewrite constraints routed to that
component. That mapping was partial on the original molecule and therefore was not a remapping. S3b
now transports those constraint references through the existing correspondence in the reverse
direction. The operation remains a split producing a correspondence; no compaction result was
added.

### Reaction reversal

Before S3c, `Reaction::reverse` materialized `self.to_reaction_span()?`, but used the span only to
obtain its right projection as the new left-hand molecule. It separately collected removed and
created ids, built sparse maps over every id the inverted deltas would reference, and used those
maps to re-anchor the inverted deltas. This worked because the hash-map representation checked
coverage only when each reference was looked up. In particular, ids introduced by `Add` deltas did
not need to form a dense source space, so those maps did not satisfy the remapping contract.

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
depends on this document.

Doc 211 has landed the compaction row and supplies the naming and layering precedent. The existing
construction sites have now been classified: only genuine total relabelings move to the dense
layer; operation-local bookkeeping and correspondence transport do not migrate merely because they
currently use `IdRemapping`.

The following work is transferred to doc 218 on the `feature/mutation-witness` branch:

- consistently generating covariant source-to-result witnesses;
- representing multi-input and multi-output witnesses, including tagged output id spaces for
  `Molecule::split`;
- composing correspondences, remappings, and compactions without collapsing every operation onto
  one carrier;
- classifying and migrating the remaining `IdRemapping` producers and consumers, including
  constraints, deltas, molecule combination, pushout, and reaction-span operations; and
- enforcing the bijection semantics of remapping and retiring `IdRemapping` when its consumers have
  valid replacements.

For example, each component correspondence returned by `Molecule::split` currently points from the
component back to the source molecule. A covariant witness for the whole split instead maps each
source atom to a tagged target `(component_id, local_atom_id)`. Repeated bare local ids are therefore
a codomain-representation defect, not evidence that remapping permits repeated images.

## Public-symbol inventory

The implementation adds:

- `Compaction<Id>::compact_vec` for applying one id-space compaction to its aligned data column;
- `Remapping<Id>` with `new`, `map`, and `try_map`;
- `index_vec::Idx` implementations for the existing graph-core and graph-IR id newtypes;
- `GraphRemapping` with `new`, `nodes`, `edges`, and the typed node/edge `map_*` and `try_map_*`
  methods;
- `MoleculeRemapping` with `new`, `graph`, six overlay-component accessors, and typed `map_*` and
  `try_map_*` methods for all eight molecule entity kinds.

The typed molecule methods use the complete entity names: `map_dative_bond`,
`map_aromatic_system`, `map_multicenter_bond`, `map_noncovalent_bond`, `map_stereo_atom`, and
`map_stereo_bond`, with matching `try_map_*` names. Atom and bond lookup remains `map_atom` and
`map_bond`.

The implementation changes or retires:

- the graph-specific free functions `compact_node_vec` and `compact_edge_vec`, replaced by the
  generic `Compaction<Id>` method;
- the current graph aggregate `Remapping`, renamed to `GraphRemapping` so the generic type can own
  the unqualified name;
- the return type of `GraphCorrespondence::to_remapping`; and
- the private graph-core operation modules, named `remap` and `compact`.

`Reaction::reverse`, `Molecule::{remap, try_remap}`, and
`ReactionSpan::{remap, try_remap}` retain their current public signatures. No remapping carrier is
added to the Python API. `IdRemapping`, molecule-level correspondence conversion, constraint and
delta transport, metadata projection, and S-group reindexing retain their current APIs in this
document's scope.

## Implementation plan

### S0 — Establish the performance baseline

#### S0a — Benchmark reaction reversal **Done**

**Module:** `umol-graph-ir/benches/reaction.rs`.

Add one top-level `Reaction::reverse` benchmark using a reaction with additions, removals,
modifications, overlays, and constraints, so it exercises the current sparse re-anchoring path.
Record the baseline result before changing the implementation and retain the same case for the final
comparison.

**Tests and evidence:** Run the existing reaction reversal tests and the filtered reaction benchmark;
record the benchmark result in this document.

**Change class:** additive evidence (green).

**Dependencies:** none.

The retained mixed reaction case includes added, removed, modified, and unchanged entities, an
aromatic-system modification, and added, removed, and unchanged constraints. Before S3c, the
sparse re-anchoring implementation measured 14.121–14.166 µs.

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

#### S3a — Keep constraint-edit handle substitution local **Done**

**Module:** `umol-graph-ir/src/ir/edit.rs` and its `ConstraintEdit` tests.

Replace the `IdRemapping` built by `ConstraintEdit::new` with one local sparse mapping from each
referenced `Entity` to the normalized `Entity` slot for its resolved handle. Apply that mapping to
the entity references in the constraint; do not introduce a transport carrier or a new public
helper.
`ConstraintEdit::resolve` remains a true dense-remapping consumer and is not changed in this
subitem.

**Tests and evidence:** Retain exact construction and resolution cases. Add cases for sparse source
ids, repeated use of one referenced entity, a missing handle, and a handle-kind mismatch. Run the
focused `ConstraintEdit` tests.

**Change class:** internal representation change (green).

**Dependencies:** [dep: S2a]

`ConstraintEdit::new` now maps referenced entities through one local sparse entity map while its
typed handle vectors retain the normalized slots; `ConstraintEdit::resolve` remains unchanged as
the dense-remapping consumer. The focused `ConstraintEdit` suite passes all 23 selected tests,
including sparse ids, repeated/shared handles, all entity kinds, missing handles, and handle-kind
mismatches. Package clippy passes with warnings denied.

#### S3b — Transport split constraints through the correspondence **Done**

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

`Molecule::split` now obtains component handles by looking up each routed constraint reference in
the reverse direction of the returned component-to-original correspondence, then resolves the
result through the ordinary `ConstraintEdit` application path. The partial original-to-component
`IdRemapping` and its construction helper are gone. All six focused split tests pass, including an
exact rich-component case covering every entity kind, relational and molecule constraints, and the
complete returned correspondences. Package clippy passes with warnings denied.

#### S3c — Reverse reactions through the reaction span **Done**

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

`Reaction::reverse` now materializes the reaction span, reverses its sides privately, and converts
the reversed span back to a reaction. The sparse reverse bookkeeping and `reversed_remapping` are
gone. The former `remap_delta` implementation remains test-only; its disposition is part of the
delta witness work transferred to doc 218.

Exact tests cover all four entity-span states, all three constraint-span states, every entity
family, stereo frames, side equivalence, and reversal as an involution after span normalization.
The existing feature-gated span-side property also passes. The full `umol-graph-ir` suite passes
6,555 unit tests plus its integration and documentation tests, and package clippy passes with
warnings denied. The retained benchmark now measures 13.462–13.534 µs, a statistically significant
4.0–4.9% improvement over the S0a baseline.

### S4 — Align the molecule remapping module name

#### S4a — Rename the private molecule remapping modules **Done**

**Module:** `umol-graph-ir/src/ir/molecule/remapping.rs`,
`umol-graph-ir/tests/property/molecule/remapping.rs`,
`umol-graph-ir/tests/property/reaction/span/remapping.rs`, and their parent module declarations.

Rename the molecule implementation module and the molecule and reaction-span property-test modules
from `remapping` to `remap`, matching the operation noun used by graph core and the rest of graph
IR. This is a private module rename only; do not change remapping carriers, producers, consumers,
operation direction, or public APIs.

**Tests and evidence:** Run `cargo test -p umol-graph-ir` and
`cargo test -p umol-graph-ir --features proptest --test property`.

**Change class:** private naming correction (green).

**Dependencies:** [dep: S3c]

The molecule implementation module and the molecule and reaction-span property-test modules are now
named `remap`; no public symbol or behavior changed. The graph-IR test suite and its feature-gated
property target pass.

### S5 — Align documentation and close

#### S5a — Complete verification and close the record **Done**

**Module:** `docs/development/nomenclature.md`, `docs/development/data-types.md`, crate rustdoc,
this document, and `discussion/000-status.md`.

Reconcile the public additions against the inventory above, preserve the explicit boundary between
the implemented total-map carrier and the deferred remapping semantics, record final verification,
and mark this document Completed only after S4a and the checks have landed.

**Tests and evidence:** Run `cargo +nightly fmt --all`, `cargo test -p umol-graph-core`,
`cargo test -p umol-graph-ir`, `cargo test -p umol-graph-ir --features proptest --test property`,
package clippy for both crates with warnings denied, and `git diff --check`.

**Change class:** documentation alignment and verification (green).

**Dependencies:** [dep: S4a]

The nomenclature and data-type guides now distinguish the intended total-bijection semantics from
the weaker contract currently enforced by the dense carriers, and they record the sparse
`IdRemapping` as legacy transport. `cargo +nightly fmt --all -- --check`, both package test suites,
the feature-gated graph-IR property target, and package clippy with warnings denied pass.

### Dependency summary

S0 through S5 are complete. Witness semantics and migration were not prerequisites for closing
this document; they are now one coherent work item in doc 218.

## Closeout

This work extracted the generic dense carrier, renamed the graph aggregate, added the corresponding
molecule aggregate, removed two improper sparse-remapping uses, and simplified reaction reversal
through `ReactionSpan`. It does not claim that the current dense constructors enforce remapping's
full bijection semantics, nor does it retire `IdRemapping`. Those coupled changes remain with the
witness design in doc 218.
