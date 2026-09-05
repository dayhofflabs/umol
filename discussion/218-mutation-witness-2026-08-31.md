# 218 — Mutation witnesses for the editing surface

Status: In Progress
Date: 2026-08-31
Relates: [179](179-python-editing-and-transactions-2026-08-02.md),
[184](184-deltas-and-edits-2026-08-04.md),
[199](199-open-container-integrity-2026-08-18.md),
[204](204-reaction-application-redesign-2026-08-19.md),
[212](212-remapping-layer-2026-08-26.md),
[data type contracts guide](../docs/development/data-types.md)

## Purpose

External consumers key their own data by entity ids — per-atom provenance
registries, correspondence caches, any structure held outside the `Molecule`.
Every mutation that moves the entity id space silently invalidates those keys:
after a removal compacts the survivors, "atom 7" is atom 6, and an external
record keyed by the old id lies. This is the accepted cost of ordinary
numerical indices over bound references (docs 114 and 139 record that
background); the payment is that **every operation that moves the id space
must be able to surrender the move as a value** — a mutation witness.

The triggering request is a downstream repair/assemble workflow: repair
(delete a misread atom, fix a bond order, correct a stereocentre — the editor
path) and assemble (build from several sources — `combine`, which already
returns a correspondence). The transactional design is not in question; what
the consumer needs becomes knowable exactly at commit and is currently not
returned.

## Current state

- The editor returns piecewise information in Rust: additions return their
  new ids, and `MoleculeEditor::remove` returns a pre-to-post
  `MoleculeCompaction`. Finalization returns only the molecule, so an editing
  session has no aggregate witness.
- `Molecule::apply(edits)` returns only the new molecule. Contrary to the
  first version of this document, it does not construct and discard a
  `MoleculeCorrespondence`; the edit machinery has operation-local handles
  and compactions but no session-level correspondence to reuse.
- Transactions retain undo payloads but expose no forward or rollback
  witness.
- `canonicalize_with_correspondence` returns a source-to-canonical
  correspondence even though canonicalization has the stronger remapping
  semantics.
- Reaction application returns `ReactionDerivation`, which duplicates the two
  sides and correspondence already representable by `Reaction` and
  `ReactionSpan` without being a complete application record.
- `combine_all` returns one input-to-combined correspondence per input.
  `combine` and `combine_from` return only the correspondence from `other`
  into the result. These mappings restate a deterministic append layout and
  are not nontrivial operation witnesses.
- `split` returns one component-to-source correspondence per component. This
  is the reverse of a covariant source-to-result operation witness.
- `MoleculeRemapping`, added by doc 212, currently has no production
  consumers. Constraint transport now consumes `MoleculeCorrespondence`;
  S1e removed `IdRemapping` and migrated its live consumers.

## Settled semantic foundation

All operation witnesses are **covariant**: they map ids in an operation's
source object or objects to ids in its result object or objects. Provenance in
the reverse direction may be useful, but it is not the operation witness.

The carrier follows the operation's semantics:

- `MoleculeCorrespondence` is the general single-source, single-result case:
  a partial bijection with explicit left and right counts for all eight entity
  families.
- `MoleculeRemapping` is the efficient special case of a total bijection with
  a dense source. It is appropriate when renumbering is the operation, as in
  canonicalization.
- `MoleculeCompaction` is the order-preserving removal case with no unmatched
  result entities. Each component declares its source count and derives its
  result count, so the witness is complete without either molecule. It is
  appropriate when selection or removal is the operation.

The special carriers are used only when their stronger property is part of
the operation's semantics; correspondence is the default otherwise. The
current dense `Remapping` and `MoleculeRemapping` constructors enforce only a
total dense source map and still permit repeated or sparse target images.
They do not yet enforce the semantic remapping contract.

### Dense-carrier audit

The current `GraphRemapping` producers confirm that the weaker contract is in
active use rather than merely permitted by its constructor:

- `GraphCorrespondence::to_remapping` requires only totality on the left, so
  the result may be an injection into a larger right-hand space;
- molecule combination maps the appended operand into an offset range of the
  larger combined molecule;
- reaction-span superimposition maps the right molecule into the larger union
  frame;
- pushout transports one operand through its coprojection into the pushout;
- split maps source atoms to component-local ids, which repeat across
  components because the component tag is absent; and
- graph-core relation transport accepts the aggregate as an arbitrary total
  participant-id map.

Only the complete molecule `remap` path first establishes totality on both
sides and therefore uses the carrier as a semantic remapping. Tightening the
current constructor in place would break legitimate injection and transport
uses. Split's aggregate map also lacks a component tag and is replaced by
separate source-to-component correspondences rather than forced into a
remapping. These transport uses migrate to correspondences with explicit
source and target counts. `MoleculeRemapping` becomes canonicalization's
witness only after its bijection contract is enforced.

The constructor-check experiment on 2026-09-04 required every image to be
in range and unique. Library tests for graph-core, graph-IR, and graph yielded
8,193 passes and 90 failures: 52 from non-bijective test fixtures and 38 from
production paths through four constructors. Reaction composition,
application, and unmapped reaction-SMILES failures traced back to those same
sites. The temporary assertions were removed after the experiment.

### Migration of non-bijective transport

- **Split:** construct one source-to-component correspondence per component.
  Use it to transport that component's participants and constraints, and
  return it from the witness-bearing form. The component collection supplies
  component identity; local ids need no component tag.
- **Pushout:** consume the existing input-to-pushout correspondences directly.
  Remove the conversion of the right correspondence into `GraphRemapping`.
- **Combination:** construct the appended operand's correspondence into the
  combined id space from its offsets and the full result counts. Use it
  internally. The public operation remains witness-free because the append
  layout determines the mapping.
- **Superimposition:** construct an RHS-to-union correspondence from the
  existing union assignment and explicit union counts. Use it to transport
  RHS entities and references.

Relation participants, stereo ligands, and constraints gain
correspondence-based transport. Every source id referenced by the transported
value must have an image. An unmatched reference is a transport failure,
never an instruction to delete the referencing entity or constraint. Each
of these four producers establishes coverage for the values it transports;
split transports only the entities selected for the current component.

Retain remapping-based transport for bijective renumbering. Both routes share
the underlying id-lookup mechanics without duplicating participant and
constraint traversal. Correspondence-based transport uses `map` and
`try_map`; bijective transport retains `remap` and `try_remap`. For graph-core
relation sets the correspondence route has this shape:

```rust
pub fn map(&self, correspondence: &GraphCorrespondence) -> Self;
pub fn try_map(&self, correspondence: &GraphCorrespondence) -> Option<Self>;
```

`try_map` returns `None` when any referenced participant lacks an image;
`map` asserts the same coverage condition. Mapping retains relation rows,
participant order, and payloads while replacing participant ids. Graph-IR
stereo ligands and entity sets also accept `&GraphCorrespondence`. The
inherent ligand methods delegate to `RelationParticipant`; they do not add
a separate atom-correspondence API. Constraints follow the same naming and
coverage contract with `MoleculeCorrespondence` for their entity references.
These methods consume a correspondence; the `with_correspondence` suffix is
reserved for operations that return a witness.

### Remapping construction and consumption

`Remapping<Id>` is an open carrier for a permutation of one dense id space.
If its image vector has length `n`, every image must lie in `0..n` and occur
exactly once. Construction therefore rejects an out-of-range or repeated
image. No separate missing-image case exists: `n` distinct in-range images
necessarily cover the target space. `Remapping::new` returns
`Result<Remapping<Id>, RemappingError<Id>>`; `len` and `is_empty` expose the
shared source and target size.

The image storage is a private `Vec<Id>`. Construction, lookup, and vector
reordering require `Id: Copy + Into<usize>`; operations that enumerate
source ids additionally require `From<usize>`. The public lookup still
accepts and returns `Id`, so different entity-id types cannot be mixed.
ID newtypes implement `From<Id> for usize`; their existing `From<usize>`
and inherent `index` methods remain. `index_vec` and its `Idx` implementations
are removed: typed lookup does not require a typed storage dependency, and
the trait never prevented integer extraction. This does not change the
existing narrowing behavior of `usize`-to-ID conversion.

The aggregate constructors receive already-valid components and are
infallible. `GraphRemapping::new` takes the node and edge `Remapping` values,
and `MoleculeRemapping::new` takes the graph aggregate and the six overlay
remappings. Validation remains at the single-id-space boundary rather than
being repeated with aggregate-specific error types.

`Remapping::remap_vec<T>(&self, values: Vec<T>) -> Vec<T>` consumes a vector
and places each source value at its image position, without requiring
`T: Clone`. `try_remap_vec` returns `Option<Vec<T>>`, with `None` exactly when
the vector length differs from the remapping length. The asserted form
panics on the same mismatch; both consume the input even on failure.
This single-space operation owns vector reordering for molecule and span
renumbering rather than separate private free functions.

Agreement with an independently supplied molecule or reaction span is
contextual, not intrinsic to the carrier. `Molecule::remap` and
`ReactionSpan::remap` accept `&MoleculeRemapping`; their checked companions
verify that all eight component lengths equal the receiver's entity counts,
while their asserted companions document and panic on the same mismatch.
The carrier itself already establishes bijectivity, so consumers do not
recheck it.

There is no correspondence-to-remapping conversion in the intended public
surface. The current `GraphCorrespondence::to_remapping` has no production
consumer and admits total-left injections, while
`MoleculeCorrespondence::to_remapping` exists only to construct the legacy
sparse `IdRemapping` inside the current correspondence-based molecule and
reaction-span remapping implementations. Remove both methods. Canonicalization
constructs `MoleculeRemapping` directly from its selected entity orders and
uses the infallible remapping-to-correspondence conversion when it needs the
general carrier. `Molecule::framed_eq_under` likewise accepts a
`&MoleculeRemapping`, because complete framed equality under a supplied id
relation requires exactly a dense total bijection. A checked narrowing can be
added later if a concrete consumer requires it; theoretical convertibility is
not enough to retain a public seam.

For a single-input, single-output edit that returns a correspondence, removed
entities are left-unmatched and added entities are right-unmatched. Additions
can therefore be recovered from the correspondence itself:

> **Addition law.** The post-state entities without a pre-image — the
> correspondence's right-unmatched entities, taken in id order — are the
> script's surviving additions, in creation order.

This holds for the current `Edits` semantics because additions append,
compaction preserves order, and an entity added and then removed in the same
script never materializes. A caller that removed some of its own `New(n)`
handles knows which additions survived, so `New(n)` resolution is
caller-computable from the correspondence alone. This law still requires
tests, including an add-then-remove case.

### Transformation classification

The same classification applies to transformations:

- **Identity class** — id-preserving transformations (`Normalize`, `Reframe`,
  attribute-only edits, bond and multicenter-bond resolution, and charge
  delocalization): the witness is the identity, stated as a law rather than
  returned.
- **Append class** — additions preserve every existing id and allocate new ids
  contiguously after the old range. That layout is part of the operation's
  contract, so pre- and post-state counts determine the mapping and the new
  ids. Primitive editor methods may still return the entity they create when
  construction needs it; that return is not a mutation witness.
- **Compaction class** — pure removals and selections, where the move is
  order-preserving and dense: the witness is `MoleculeCompaction`, the
  stronger type whose shape carries that guarantee. `MoleculeEditor::remove`
  already conforms.
- **Remapping class** — dense bijective renumbering, where renumbering is the
  operation: the witness is `MoleculeRemapping`.
- **Correspondence class** — general additions, removals, and replacements:
  the witness is `MoleculeCorrespondence`.

Class membership is assigned per operation by its effect on the id space,
regardless of the owning layer, and the effect is judged per family. The
chemistry-layer transformers currently move at most one family: the
`Aromatizer` appends aromatic-system entities and the `Kekulizer` removes all
of them; their bond-order and constraint changes are attribute changes.
Neither needs a returned witness. The former has the append law. In the
latter, the affected family has no survivors to transport and every other
family is unchanged. Doc 166's hydrogen transformations have the same split:
unfolding has the append law, while folding selected hydrogens has a
nontrivial `MoleculeCompaction`.

`BondsResolver` and `MulticenterBondsResolver` only modify attributes and
preserve every id. `AromaticityResolver` and `StereoResolver` may replace
overlay entities, and the aggregate `Resolver` runs both structural phases.
That replacement is real, but resolver results are deliberately excluded from
the witness-bearing API in this work rather than expanding every successful
solution with transport data.

### Public operation audit

The audit covers public operations whose input contains one or more existing
molecule id spaces and whose result changes or materializes a molecule id
space. It excludes construction from no molecule source (`MoleculeBuilder`,
parsers, and external-format conversion), read-only selection and planning
(`induced_subgraph`, `edits`, resolver plans, and perception), and operations
that transform a `Reaction` or `ReactionSpan` without issuing a molecule
result. These exclusions have no source molecule ids to transport or belong
to a different aggregate's witness contract.

| Operation family | Shape and id effect | Current return | Required witness conclusion |
| --- | --- | --- | --- |
| Direct molecule attribute and constraint mutation; `Normalize`; `Reframe`; `DelocalizeCharge`; bond and multicenter-bond resolution | one to one, every entity id preserved | value or `()` | Identity law; no returned carrier. |
| Editor `add_*`; `AromaticityPerceiver::add_systems` | one to one, append | new id for individual editor additions; `add_systems` returns `()` | Existing ids and allocation order are contractual, so counts determine all new ids. Keep primitive created-id returns needed for continued construction; do not add witness or created-id returns to higher-level append operations. |
| `Aromatizer`; `Fragment::finish_open`; hydrogen unfolding | one to one, high-level append transformation | bare result | Append law; no witness-bearing form. |
| Editor relation-family `remove_*` | one to one, order-preserving removal in a selected relation family | `()` | Keep plain methods witness-free; add `with_compaction` forms where surviving ids can move. Each method already constructs the required family compaction internally. |
| `Kekulizer` | one to one, removes every aromatic system and preserves every other family | bare result | No returned witness: the affected family has no survivors and every other id is unchanged. |
| `MoleculeEditor::remove` | one to one, order-preserving removal with cascades | `MoleculeCompaction` | Keep `remove` witness-free and provide `remove_with_compaction`. |
| `Molecule::extract` | one to one, selection followed by order-preserving removal | molecule only | Keep `extract` returning the molecule and provide `extract_with_compaction`. The supplied subgraph correspondence is a sub-to-host selection descriptor, not the result witness. |
| `Molecule::apply`; editor `apply`, `transact`, `snapshot`, and finalization; transaction rollback | one to one, arbitrary additions, removals, and modifications | molecule, editor, transaction, or `()` without an aggregate witness | Witness-bearing application returns a source-to-result correspondence; rollback returns its inverse. Witness-bearing snapshot/finalization exposes the session correspondence when the editor originated from an existing molecule. Plain variants return no witness. |
| Aromaticity, stereo, and aggregate resolution | one to one, replacement-capable in overlay families | successful `Solution` without a witness | No witness-bearing form in this work. |
| `Transformer::transform` and `generate_all` | one source to one result, or to alternative results; effect depends on transformer | bare molecule values | No blanket witness method. Add an operation-specific form only when a transformer has nontrivial transport that callers need. |
| `canonicalize_with_correspondence`; `Molecule::remap`; `ReactionSpan::remap`; `Molecule::framed_eq_under` | one to one, dense bijective renumbering | canonical value plus correspondence; remapped value, remapped span, or Boolean from a caller-supplied correspondence | Canonicalization's witness-bearing form becomes `canonicalize_with_remapping`; the other operations accept `MoleculeRemapping`. Remapping operations return no witness because the caller supplied it. |
| `Reaction::apply_at` and `Reaction::apply` | one host to one product per application | `ReactionDerivation` with host-to-product `comap` | Remove `ReactionDerivation`. Offer product-only, product-with-correspondence, realized-`Reaction`, and realized-`ReactionSpan` application forms. |
| `combine`, `combine_all`, and `combine_from` | one or many inputs to one result by disjoint append | operand-to-result or all input-to-result correspondences | Append order determines every input mapping. Return only the combined molecule, with preservation and ordering made contractual. |
| `meet_pushout` | two inputs to one result with possible identification | `MoleculePushout` containing both input-to-result correspondences | The correspondences are nontrivial. Keep a witness-bearing result; align the plain and `with_correspondence` forms with the naming rule. |
| `split` | one source to several simultaneous component results | one component-to-source correspondence per result | Plain `split` returns components. A witness-bearing form returns each component with its source-to-component correspondence. The collection supplies component identity. |
| `Fragment::attach` and fragment `Add` | two fragment bodies to one result body by append | fragment only; the internal `combine` correspondence is discarded | Append law; no witness-bearing form and no extra created-id return. `finish` alone is identity. |
| `React` for `Molecule` and `[Molecule]` | one or many reactants to alternative sets of product components | component molecules only | Keep the convenience operation product-only. Callers needing transport compose the explicit combine, reaction-application, and split operations. |

An operation witness explains the source-to-result entity pairing so callers
can follow atom provenance through a sequence of operations. It is not an
audit trace and does not promise operational persistence through every
internal removal and addition. The operation constructs its atom
correspondence; internal remove/add mechanics alone do not determine the
pairing. Incidence-based induction of bonds and overlays may pair equivalent
entities across an internal remove/add sequence. That is consistent with this
contract and is not a defect requiring an action journal.

### Reaction application

This document subsumes doc
[204](204-reaction-application-redesign-2026-08-19.md). `ReactionDerivation`
is removed. It stores the same side pair and correspondence already
representable as a `Reaction` or `ReactionSpan`, while retaining neither the
originating reaction nor the rule-to-host match needed for a complete
application record. Its exact rhs frame does not justify another semantic
reaction type, and its unconstrained `chain` operation does not check that the
independently supplied intermediate sides agree.

Reaction application instead exposes the useful results directly. A
successful application may be requested as:

- the product molecule alone;
- the product and its covariant host-to-product
  `MoleculeCorrespondence`;
- the realized `Reaction`; or
- the realized `ReactionSpan`.

These are methods on `Reaction`, with the following names:

| Primary result | Supplied match | Iterate over matches |
| --- | --- | --- |
| Product molecule | `apply_at` | `apply` |
| Product and correspondence | `apply_at_with_correspondence` | `apply_with_correspondence` |
| Realized reaction | `apply_at_to_reaction` | `apply_to_reaction` |
| Realized reaction span | `apply_at_to_reaction_span` | `apply_to_reaction_span` |

The correspondence form returns `(Molecule, MoleculeCorrespondence)`. The
reaction and span forms return those primary objects directly. Iterators
follow the same item contracts, including operation errors, without a
replacement result wrapper for `ReactionDerivation`. Rust and Python expose
the same names and result shapes. The `React` trait methods remain separate,
providing the product-component convenience workflow.

The graph-core surface supports the same distinctions rather than supplying
one universal answer. Graph additions return their new ids; removals return a
`GraphCompaction`; pushout returns both input-to-object coprojections; and
subdivision has an operation-specific result because source edges become
different result entities. Conversely, the `PushoutComplement::context` arrow
is context-to-host because that is its categorical meaning. It is not a
covariant host-to-context mutation witness and should not determine the public
molecule-operation direction.

### Multiple inputs and outputs

The single-pair aggregate types are insufficient by themselves when an
operation changes object arity.

The public representation should retain the operation's individual
correspondences rather than flatten every object into one tagged carrier when
the mappings are nontrivial. `combine`, `combine_all`, and `combine_from` are
deterministic append layouts, so input counts and ordering recover every
input-to-result mapping without returning a carrier. `meet_pushout` is
different: overlapping input entities may be identified, so its
witness-bearing result retains the distinct left-to-result and right-to-result
correspondences.

For `split`, each output component can carry its own source-to-component
correspondence. The result collection identifies the component, so the entity
ids do not need to contain a component tag. Considering all correspondences together
distinguishes an entity routed to another component from an entity absent from
every result.

This leaves no present need for a universal multi-object witness type. A
future operation requiring a flattened view can add one from demonstrated
composition requirements rather than making every current operation expose
tagged ids.

### Composition across carrier types

Composition must preserve covariant direction while allowing the strongest
operational carrier at each step. For example, an editor removal should not
stop returning a compact `MoleculeCompaction` merely so it can compose with a
later general edit.

The general composition vocabulary is correspondence. A compaction carries
its pre-state counts and can therefore be lowered losslessly: enumerate its
survivors as matched pairs, use removed ids as the unmatched left side, and
set the right count to the survivor count. A semantic remapping can be lowered
directly to a total correspondence with equal carrier counts.
Primitive append information can be lowered from the preservation law and the
pre- and post-state counts: existing ids map identically and created ids are
right-unmatched. Composition may then return the general correspondence even
when its inputs were specialized carriers or did not allocate a witness at
their direct API boundary.

Multiple-object composition reconstructs append mappings from the input
layout, then composes the individual correspondences through the shared
object. A caller that needs end-to-end transport across combination, reaction
application, and splitting performs those explicit operations and composes
their mappings. `React` does not expose that richer result. A pushout retains
its distinct input-to-object correspondences through later composition; it is
not collapsed into a relation that permits arbitrary many-to-many
"identity".

Subtype-to-correspondence conversion follows the information carried by each
type. Once `MoleculeRemapping` enforces its bijection contract, lowering it is
infallible and context-free, so `From<&MoleculeRemapping> for
MoleculeCorrespondence` is the idiomatic baseline. An owned `From` is added
only if an owned call site needs it.

The inverse narrowing is not exposed. Neither current `to_remapping` method
produces a needed semantic remapping after the remapping consumers and
canonicalization use `MoleculeRemapping` directly.

`Compaction<Id>` is a complete witness between two finite dense id spaces. It
stores the source count and the sorted, deduplicated removed ids; the result
count is the source count minus the number removed. Its construction surface
is:

```rust
pub enum CompactionError<Id> {
    RemovedIdOutOfRange { id: Id, source_count: usize },
}

pub fn new(
    source_count: usize,
    removed: Vec<Id>,
) -> Result<Compaction<Id>, CompactionError<Id>>;

pub fn identity(source_count: usize) -> Compaction<Id>;
pub fn source_count(&self) -> usize;
pub fn result_count(&self) -> usize;
```

`new` rejects a removed id outside `0..source_count`; input order and repeated
removed ids do not matter. There is no separate asserted public constructor:
operation producers already know their source counts and handle an impossible
construction error internally. A context-free `Default` or `empty` cannot
mean identity over an undeclared id space and is removed; callers use
`identity(source_count)`.

Compaction consumers respect these declared bounds. `compact` returns `None`
for a removed or out-of-source-range id. `uncompact` asserts that its input
is below `result_count`; `try_uncompact` returns `None` otherwise.
`compact_vec` requires exactly `source_count` input elements and asserts
that condition; `try_compact_vec` returns `None` on a length mismatch.

The aggregate constructors accept already-valid components and are
infallible. `GraphCompaction::new` takes the node and edge `Compaction`
values, while `MoleculeCompaction::new` takes the graph aggregate and the six
typed overlay compactions. Validation is not repeated in aggregate-specific
error types. Graph, relation-set, and molecule-editor removal producers all
know the pre-removal table sizes required to construct these components.

Lowering is consequently context-free and infallible:
`From<&Compaction<Id>> for Correspondence<Id>`,
`From<&GraphCompaction> for GraphCorrespondence`, and
`From<&MoleculeCompaction> for MoleculeCorrespondence`. Each component
correspondence uses the stored source count as its left count, the derived
result count as its right count, and maps every survivor to its compacted id.
Owned `From` implementations are added only if an owned call site needs them.

No conversion from a list of added entities is needed. `MoleculeCorrespondence`
is an open carrier, and callers can construct the append correspondence
directly from the known counts and allocation order.

Composition remains an operation on correspondences rather than accepting a
generic conversion argument. Specialized witnesses are lowered explicitly,
then `MoleculeCorrespondence::compose(&self, next: &MoleculeCorrespondence)`
returns `Result<MoleculeCorrespondence, _>`. It verifies, for every entity
kind, that `self`'s right count equals `next`'s left count. Equal counts cannot
prove that the two open carriers refer to the same intermediate molecule;
that provenance is intentionally not represented.

The generic `Correspondence<Id>::compose` enforces the same intermediate-count
condition. `compose_all` returns `Result<Option<Self>, _>`: `None` denotes an
empty input, while `Err` denotes incompatible consecutive carriers. An
`Into<&MoleculeCorrespondence>` argument cannot represent conversions that
materialize a new value, while `Into<MoleculeCorrespondence>` would hide the
allocation performed by lowering. The replacement of non-semantic
`IdRemapping` consumers follows the migration above; exact composition error
variants are implementation details within the settled count-agreement
contract.

### API direction

An operation must not leave id behavior implicit. Identity and append
operations state their layouts as contracts. General mutations expose a
witness only when the mapping cannot be recovered from the operation and the
pre- and post-state shapes.

The method without a `with_*` suffix never returns a witness. A
witness-bearing companion names its carrier precisely:

- `with_correspondence` for `MoleculeCorrespondence`;
- `with_compaction` for `MoleculeCompaction`; and
- `with_remapping` for `MoleculeRemapping`.

This requires breaking alignment of existing methods whose plain form
currently returns a witness. The intended pairs include
`MoleculeEditor::remove` / `remove_with_compaction`, `Molecule::extract` /
`extract_with_compaction`, and the applicable relation-family removal methods.
`Molecule::apply`, editor finalization, and the transaction path gain
correspondence-bearing companions while retaining their existing primary
results in the plain forms. The witness-bearing canonicalization method is
`canonicalize_with_remapping`, not `canonicalize_with_correspondence`.

`split` returns component molecules; `split_with_correspondence` returns each
component paired with its covariant source-to-component correspondence.
Transaction `rollback_with_correspondence` returns the inverse witness while
plain `rollback` retains its ordinary result.

The molecule pushout pair is:

```rust
meet_pushout(...) -> Option<Molecule>
meet_pushout_with_correspondence(...)
    -> Option<(Molecule, MoleculePushoutCorrespondence)>
```

`MoleculePushoutCorrespondence` contains named `left` and `right`
`MoleculeCorrespondence` components, each directed from its input to the
result. Their right counts agree for all eight entity kinds. It contains no
molecule and replaces the existing `MoleculePushout` result container. The
two component correspondences can each be composed with a subsequent
result-to-result correspondence.

The same naming and output/witness separation applies across graph-core and
graph-IR. Graph pushout uses `GraphPushoutCorrespondence`; relation-set
pushout uses `RelationPushoutCorrespondence`. Their bare methods return the
graph or relation set, and `pushout_with_correspondence` returns the output
paired with the two input-to-result correspondences. Graph pullback uses
`PullbackCorrespondence`, and graph pushout complement uses
`PushoutComplementCorrespondence`, through the corresponding
`*_with_correspondence` methods. Their categorical directions are preserved:
pullback `left` and `right` map the result to the inputs; complement `context`
maps context to host and `interface` maps interface to context. These
categorical correspondences are not reinterpreted as covariant mutation
witnesses. Graph removals use `*_with_compaction` companions.

Chemistry-layer operation API policy in `umol-graph` is settled separately.
This work adapts its callers of changed graph-core and graph-IR APIs without
adding chemistry-operation witness variants. The transformer audit above is
context for that later work, not an instruction to implement those variants.

Identity operations, append-only operations, caller-supplied remapping,
resolver operations, and the `React` convenience API have no witness-bearing
variant. In particular, `combine`, `combine_all`, `combine_from`, fragment
assembly, aromaticity perception, aromatization, hydrogen unfolding, and
Kekulization return no transport data beyond primitive ids needed to continue
construction. Their id behavior is contractual and caller-computable.

A blanket witness method on `Transformer` is not justified. A structural
transformer may gain an operation-specific pair only when it has nontrivial
transport that callers need.

**Laws.**

- Output equivalence: for identical inputs, bare and witness-bearing variants
  have the same success/failure semantics and produce the same ordinary
  output. Discarding the witness recovers the bare result; mutating variants
  also leave the same resulting state, and iterators yield corresponding
  outputs in the same order. The bare implementation need not allocate and
  discard a witness to satisfy this law.
- Composition: the witness of a sequence of mutations is the composition of
  the step witnesses; the editor's session correspondence equals the
  composition of its piecewise witnesses.
- Inversion: rollback's witness is the inverse partial bijection; entities
  created by the rolled-back mutation are absent from it.
- Identity: an attribute-only mutation yields the identity correspondence,
  so consumers need no moved-or-not case split.

**Python parity and shape.** Python should catch up to the resulting Rust
surface: `MoleculeCompaction` (recently redesigned; the bindings are behind)
and the applicable variant pairs, mirrored by name. The
scikit-learn keyword shape
(`with_correspondence=True` selecting a richer return) was considered and
rejected: a flag-dependent return type needs `Literal` overloads to type
precisely and degrades to a union under a runtime flag, while method pairs
type exactly — decisive for a surface that maintains a signature inventory.
The existing `canonicalize_with_correspondence` binding migrates with Rust to
the remapping-bearing name and return type. Python exposes frozen `Remapping`
values constructed from image lists, and frozen `MoleculeRemapping` values
constructed from eight remapping components named `atoms`, `bonds`, and the
six overlay kinds. It does not expose `GraphRemapping`. Python mirrors the accepted
method pairs and `(output, witness)` shapes, including
`MoleculePushoutCorrespondence`. Existing APIs may break to follow these
contracts; no return-type-changing witness flag is introduced.

## Related work

- Doc [204](204-reaction-application-redesign-2026-08-19.md) is superseded by
  this document. Its diagnosis of `ReactionDerivation` and reaction
  application result requirements is incorporated above.
- Doc [212](212-remapping-layer-2026-08-26.md): the id-transport
  foundation is complete. It introduced the dense carriers but deliberately
  left their full bijection contract, the remaining `IdRemapping` migration,
  and cross-witness composition to this document. There is no existing
  `Molecule::apply` correspondence to reuse.

## Boundaries

- No stable ids, bound references, or identity infrastructure.
- Participant-frame positions are outside the witness: correspondences move
  entity ids, and frame-action transport remains the reframe and
  canonicalization machinery's concern.
- No unconditional witness returns on the plain operations.
- No witness-bearing variants for operations whose identity behavior is fully
  expressed by an identity or append-preservation law.
- No migration from `IdRemapping` merely because a consumer needs an id map;
  the producer must first establish the selected witness semantics.
- No universal multi-object carrier without a demonstrated operation that
  cannot be represented clearly by an operation-specific collection of
  correspondences.
- No editor API redesign beyond the witness-bearing variants and bindings.

## Remaining implementation work

- Add source counts and checked construction to compaction, migrate its
  producers, and implement the settled infallible remapping-to-correspondence
  and compaction-to-correspondence conversions. Implement checked
  correspondence composition and define its exact error types.
- Enforce permutation construction in `Remapping`, compose
  `GraphRemapping` and `MoleculeRemapping` from valid components, migrate the
  remapping consumers, and remove correspondence-to-remapping narrowing.
- Migrate the four identified non-bijective transport producers and their
  consumers to correspondence-based `map`/`try_map`; replace remaining
  `IdRemapping` uses consistently with their operation semantics.
- Record and test existing-id preservation and allocation order for primitive
  additions, `combine`, `combine_all`, `combine_from`, and fragment assembly;
  add the covariant witness-bearing `split` form.
- Implement source-to-result correspondence accumulation for `Molecule::apply`,
  editor finalization, and transactions, and the inverse witness for
  `rollback_with_correspondence`.
- Implement witness-bearing variants for direct relation removals and extraction
  in graph-core and graph-IR; chemistry-layer transformer APIs remain separate.
- Replace `MoleculePushout` with the accepted plain molecule result and
  `(Molecule, MoleculePushoutCorrespondence)` witnessed result.
- Replace `ReactionDerivation` in the Rust and Python application iterators
  with the direct product, product-with-correspondence, `Reaction`, and
  `ReactionSpan` forms.
- Verify output equivalence for every bare/witnessed pair and mirror the
  settled Rust surface in Python.

## Implementation plan

All stages begin unchecked. Each subitem includes its focused tests, rustdoc,
and affected caller updates. An additive item keeps the tree green; a
breaking item may be red during editing but includes the migrations needed
to restore it. Every stage ends green. Public names and semantics above are
the contract; stop if implementation requires a new public abstraction or a
different failure or construction boundary.

Python compile adaptations accompany the Rust change that requires them;
the later Python stage completes exposure and behavioral parity. Changes in
umol-graph and umol-io are limited to adapting consumers. No mutating git
operations are part of this plan.

The remaining order follows the transport dependencies. First separate edit
handle resolution from molecular correspondence, then migrate constraint
transport and its callers together. Do not introduce an adapter between
`IdRemapping` and correspondence to preserve the old constraint API.
Strict remapping construction follows the four non-bijective producer fixes;
it does not depend on compaction. True constraint renumbering uses the
approved remapping-to-correspondence conversion once that conversion exists.
Completed S0/S1a/S1b items retain their identifiers. Dependencies below name
prerequisites, not merely the preceding item in the chosen execution order.

### S0 — Regression cases and measurement baseline

- [x] **S0a — Transport regressions.** Graph-core relation/rewriting tests and
  graph-IR molecule/reaction-span tests. Additive (green). [dep: none]
  Retain focused cases for combination offsets, split into multiple
  components, partial-overlap pushout, and RHS-to-union transport. Check
  complete output and expected pairings, including stereo participants and
  constraint references. Distinguish these valid operations from constructor
  fixtures that currently permit non-bijective remappings. Do not alter the
  constructor yet.
  Completed 2026-09-04: full-output and correspondence assertions cover
  combination offsets, a split with stereo and constraints in the later
  component, partial-overlap graph/molecule pushout, and reaction-span
  superimposition. Relation transport uses the injection produced by a
  partial-overlap pushout, checking participant order and payload retention.
  Production code and non-bijective constructor fixtures remain unchanged.
  Verification: `cargo test -p umol-graph-core -p umol-graph-ir --lib --no-fail-fast`
  passed (756 core, 6,555 IR; 3 ignored); `cargo +nightly fmt --all --check`
  and `git diff --check` passed.
- [x] **S0b — Operation benchmarks.** Existing graph-core `algorithms` and
  graph-IR `reaction`/`canonicalize` benchmark targets. Additive (green).
  [dep: S0a] Add focused removal, edit-batch, combination, split, and pushout
  measurements beside the existing reaction application and renumbering
  measurements. Keep setup outside timing, record input sizes and selected
  algorithms, and save a Criterion baseline before implementation. Later
  add witnessed variants on the same cases. These fixtures supply regression
  evidence, not claims about production scale.
  Completed 2026-09-04: added 14 operation cases and saved `mutation-s0`
  across 30 cases, including existing reaction and renumbering measurements.
  S0 gates passed: core/IR library, integration, and doc tests (7,404 passed,
  6 ignored); property targets with `PROPTEST_CASES=256` (483 passed,
  1 ignored); `cargo check --workspace` with Python 3.13.15 activated;
  `cargo +nightly fmt --all --check`; and `git diff --check`.

#### S0 benchmark fixtures

The new `graph_mutation` and `molecule_mutation` groups use paths of 8 and
64 nodes/atoms with 7 and 63 bonds. Removal deletes the middle node/atom
and its two incident bonds. The three-edit batch removes the last carbon,
adds oxygen, and bonds it to the preceding carbon. Combination concatenates
two equal-sized paths; pushout identifies one endpoint of each, producing
15/127 atoms and 14/126 bonds. Split starts with 8/64 atoms in two equal
paths (6/62 bonds). These inputs have no overlays or constraints.

Fixture construction, edit-batch cloning, and mutable input preparation are
outside timing. Removal measures the graph/editor mutation and returned
compaction, without publishing a molecule. `Molecule::apply` includes its
editor construction, checked application, and publication. The other cases
measure the complete existing operation, including any witness it currently
returns. Future bare/witnessed variants must use these same inputs.

The retained `reaction` group uses its existing six-atom/six-bond aromatic
application case and five-atom reaction-span reversal fixture. Matching
selects GraphAndOverlays, VF2, and Vismara. New mutation cases expose no
algorithm choice; molecule splitting uses union-find. The retained
`canonicalize/remapping` group measures supplied reverse renumberings over
its thirteen existing cases, including stereo, overlays, and constraints; it does
not run canonical labeling. Its benchmark ids record all eight entity counts.

Baseline recorded on 2026-09-04, macOS arm64, Rust 1.96.0, optimized bench
profile with debug information. Criterion used 30 samples, one second of
warmup, and two seconds of measurement per case. The three targets ran
sequentially after builds and verification finished. Saved estimates and
samples are under `/Users/dr/.cargo-target/criterion`, in each case's
`mutation-s0` directory. This local artifact is not tracked; the commands
below recreate it. Point estimates in microseconds:

| Operation | 8 atoms/nodes | 64 atoms/nodes |
| --- | ---: | ---: |
| Graph removal | 0.220 | 1.083 |
| Graph pushout | 2.146 | 21.152 |
| Molecule editor removal | 1.516 | 6.017 |
| Molecule three-edit application | 2.725 | 12.282 |
| Molecule combination | 2.711 | 19.330 |
| Molecule split | 4.646 | 51.909 |
| Molecule pushout | 7.810 | 67.701 |

The existing reaction case measured 1.695 µs for matching, 4.321 µs for
matched application, and 13.756 µs for reversal.
Retained renumbering estimates range from 2.032 to 15.612 µs; the eight-atom
overlay-heavy case measured 7.849 µs. Individual estimates (µs):

| Renumbering case | Time |
| --- | ---: |
| Ordinary naphthalene | 3.369 |
| Disconnected rings | 3.519 |
| Overlay-heavy | 7.849 |
| Tetrahedral stereo | 3.457 |
| Cis/trans stereo bond | 3.651 |
| Mixed atom and bond stereo | 5.489 |
| Frame-relative stereo constraint | 3.625 |
| Meso dichlorobutane | 4.528 |
| Para-stereo trichloropentane | 5.538 |
| Para-stereo cascade | 15.612 |
| Feature-free connected | 3.247 |
| Feature-free disconnected | 3.198 |
| Symmetry-heavy radicals | 2.032 |

```sh
cargo bench -p umol-graph-core --bench algorithms -- graph_mutation --save-baseline mutation-s0 --warm-up-time 1 --measurement-time 2 --sample-size 30
cargo bench -p umol-graph-ir --bench reaction -- --save-baseline mutation-s0 --warm-up-time 1 --measurement-time 2 --sample-size 30
cargo bench -p umol-graph-ir --bench canonicalize -- canonicalize/remapping --save-baseline mutation-s0 --warm-up-time 1 --measurement-time 2 --sample-size 30
```

For later comparisons, replace `--save-baseline mutation-s0` with
`--baseline mutation-s0`; do not overwrite the baseline after changing the
implementation. These are bounded regression measurements, not corpus-scale
estimates.

### S1 — Correspondence composition and transport

- [x] **S1a — Composition.** Graph-core and graph-IR correspondence modules.
  Breaking (red→green). [dep: S0a] Make single-space, graph, and molecule
  composition check intermediate counts. Migrate every caller, including
  rewriting and reaction composition; retain errors only at public boundaries
  where independently supplied counts can disagree. Implement the settled
  `compose_all` empty-input contract. Test exact count mismatches, identity,
  associativity on compatible carriers, and unmatched entities through a
  sequence. Preserve the prohibition on object-identity checks.
  Completed 2026-09-04. All three correspondence layers now return `Result`
  from `compose`, and `Result<Option<Self>, _>` from `compose_all`.
  Constructors, conversions, and open-carrier status are unchanged; count
  agreement is checked at composition, not construction. New public errors:
  `CorrespondenceComposeError` carries both intermediate counts;
  `GraphCorrespondenceComposeError` identifies nodes or edges;
  `MoleculeCorrespondenceComposeError` carries the entity kind and count error.
  `Graph::pullback` and the temporarily retained `ReactionDerivation::chain`
  expose these errors. `pushout_complement` retains `Option`, rejecting
  incompatible intermediate/host counts; reaction composition retains its
  existing surface because its applications share a producer-established host.
  Python composition/chaining raises `ValueError` for count mismatches.
  Identity and associativity properties use compatible carriers; exact tests
  cover both count directions, empty carriers, every aggregate component,
  and mismatches after an initially compatible composition.
  Verification: core/IR unit tests (7,329 passed, 3 ignored), properties with
  `PROPTEST_CASES=256` (483 passed, 1 ignored), Python-binding Rust tests
  (1,634 passed, 2 ignored), rebuilt-extension molecule/reaction pytest
  (178 passed, 2 skipped), workspace all-targets check, core/IR all-targets
  clippy with `-D warnings`, nightly formatting, and `git diff --check` passed.
- [x] **S1b — Graph participant transport.** Graph-core relation participants
  and all five relation-set families. Additive (green). [dep: S1a]
  Add `map`/`try_map` using `GraphCorrespondence`; share traversal with
  remapping transport. Test covered subsets, missing images, atom/edge
  factors, nonidentity mappings, preserved row order, and positional payloads.
  Completed 2026-09-04. `RelationParticipant` and all five relation-set families
  expose `map`/`try_map` over `GraphCorrespondence`. Mapping requires images
  only for referenced ids, preserves rows/frames/payloads, and rebuilds
  incidence indexes. Mapping and remapping share private participant traversal;
  constructors and public remapping contracts are unchanged. The existing
  `StereoLigand` trait implementation was adapted, including virtual anchors;
  inherent molecular reference transport remains S1d.
  Exact tests cover partial coverage, missing/out-of-range references in both
  factors, asserted failures, identity, inverse, composition, and incidence.
  Generated permutation cases check composition, inverse, and remapping parity
  for all five families. Verification: core/IR unit tests (7,377 passed,
  3 ignored), relation properties with `PROPTEST_CASES=256` (10 passed),
  workspace all-targets check with Python 3.13.15 activated, core/IR all-targets
  clippy with `-D warnings`, nightly formatting, and `git diff --check` passed.
- [x] **S1c — Edit handle resolution.** Graph-IR `ConstraintEdit` construction
  and resolution. Internal rewiring (green). [dep: S0a]
  Remove this consumer of `IdRemapping` before changing constraint transport.
  Resolve the existing per-kind handle indices as an edit operation, using the
  existing edit-local entity traversal rather than manufacturing a molecular
  witness. Preserve handle interning, kind checks, resolution failures, and
  nested constraints. Test all entity kinds and repeated handle references.
  No new public lookup abstraction or witness API is part of this item.
  Completed: resolution uses the existing edit-local entity traversal and a
  `HashMap<Entity, Entity>` containing only the distinct references in the
  constraint. Public signatures and resolution errors are unchanged.
  Exact tests cover all eight entity kinds, repeated references, nested
  constraints, subsets, and each resolver's failure. Verification: 32 focused
  constraint-edit tests and 6,580 graph-IR library tests passed (3 ignored);
  graph test-target checks with `proptest`, nightly formatting, and
  `git diff --check` passed.
- [x] **S1d — Ligand and entity-set transport.** Graph-IR ligand and
  entity-set methods. Additive (green). [dep: S1b]
  Add correspondence-based `map`/`try_map` with the accepted coverage
  contract, using the graph-core participant transport already implemented.
  Test atom/bond references, virtual ligand anchors, preserved entity rows,
  payloads, and frames. An unmapped reference must fail rather than disappear.
  Keep frame transport distinct from id transport.
  Completed: `StereoLigand` and all twelve ordinary/span entity-set wrappers
  expose `map`/`try_map` over `GraphCorrespondence`, delegating to existing
  graph-core transport. Constructors and remapping methods are unchanged.
  Exact tests cover nonidentity mappings, unmatched unused ids, missing
  participant images, virtual anchors, atom/bond sites, positional payloads,
  span attributes, row/frame preservation, and asserted failures.
  Verification: graph-IR library tests (6,659 passed, 3 ignored), graph-IR
  all-targets clippy with `-D warnings`, graph/IO all-targets checks, nightly
  formatting, and `git diff --check` passed.
- [x] **S1e — Constraint transport and legacy removal.** Graph-IR constraint,
  molecule, reaction-span, correspondence, and delta modules. Breaking
  (red→green). [dep: S1c, S1d]
  Replace sparse `IdRemapping` constraint transport with correspondence-based
  `map`/`try_map`. Migrate every live caller in this same item: pushout uses
  its existing correspondence; span union/side assignments use explicit
  source/result counts; molecule/span renumbering uses its already-supplied
  correspondence. Check coverage only for references actually transported,
  including every entity kind, nested constraints, and optional subsets;
  preserve predicates and frame positions. Do not retain the old constraint
  signature through a shim or a second traversal.
  Remove the unused `remap_delta` and its orphaned tests, both
  correspondence-to-remapping conversions, and `IdRemapping` after its last
  caller is migrated. Remove value-only no-op transport methods where they
  exist solely for that legacy traversal. Keep the current public molecule/
  span renumbering signatures until S2b. Verify live edit, projection,
  superimposition, pushout, and renumbering tests and absence of retired names.
  Completed: `Constraint`, `MoleculeConstraint`, and `RelationalConstraint`
  expose correspondence-based `map`/`try_map`; value-only constraint payloads
  are preserved directly. Combination and pushout reuse their correspondences;
  span transport declares the relevant id-space counts. Side projection still
  assigns absent entities beyond the valid side prefix, preserving integrity
  diagnostics; that internal assignment permutes the full union-sized domain,
  not the selected molecule's domain.
  Removed `IdRemapping`, both narrowing conversions, the unused delta helper
  and its tests, and orphaned no-op/offset helpers. The composition identity
  property now checks correspondence images directly without changing its law.
  Live guides are synchronized. Verification: core/IR library tests (7,471
  passed, 3 ignored), graph/IO library tests (4,297 passed), core/IR properties
  with `PROPTEST_CASES=256` (487 passed, 1 ignored), workspace all-targets
  check with Python 3.13.15, core/IR all-targets clippy with `-D warnings`,
  nightly formatting, and `git diff --check` passed.
- [x] **S1f — Four graph transport producer migrations.** Graph-IR molecule combination,
  split, pushout, and reaction-span superimposition. Breaking internal
  rewiring (red→green). [dep: S1d, S1e] Replace their non-bijective
  `GraphRemapping` construction with the correspondences specified above.
  Reuse the constraint correspondences established in S1e where applicable;
  do not reconstruct a legacy sparse carrier.
  Keep existing public return shapes until their dedicated stages. Use
  actual source/result counts and the existing pushout correspondences;
  split constructs a separate correspondence for each selected component.
  Run S0a regressions and downstream composition/application/ingest tests.
  Completed: combination and RHS-to-union participant transport use counted
  graph correspondences; pushout passes its existing right correspondence to
  overlay gluing. Split constructs separate source-to-component atom, bond,
  graph, and molecule correspondences, using them for participants and
  constraints. Its public result still reverses the molecule correspondence
  to the existing component-to-source direction until S5a. Removed the
  aggregate split remapping and duplicate local bond lookup. The four paths
  no longer construct `GraphRemapping`; genuine renumbering is unchanged.
  Added an exact interleaved-id split regression covering three components,
  repeated component-local ids, constraint routing, and complete outputs.
  S1 gates passed: core/IR library, integration, and doc tests (7,565 passed,
  6 ignored); core/IR properties with `PROPTEST_CASES=256` (487 passed,
  1 ignored), with both split properties rechecked after the final change;
  graph/IO library tests (4,297 passed); explicit MOL/SDF/SMILES conformance
  (12,689 passed); workspace all-targets check with Python 3.13.15;
  core/IR all-targets clippy with `-D warnings`; nightly formatting; and
  `git diff --check`.

### S2 — Bijective remapping and renumbering

- [x] **S2a — Strict remapping construction.** Graph-core remap and graph-IR
  remap modules. Breaking (red→green). [dep: S1f]
  Enforce dense bijective images, assemble graph/molecule aggregates from
  validated components, and add infallible widening to correspondences.
  Replace non-bijective lookup fixtures with valid permutations; retain
  invalid images as exact constructor rejection cases. Migrate all remaining
  constructors in the same item. The original four production regressions
  must now pass with the permanent construction checks.
  Completed 2026-09-04: construction checks images with a temporary typed
  boolean vector and returns `RemappingError::ImageOutOfRange` or
  `DuplicateImage`. Storage remains private; `len` and `is_empty` expose its
  size. Aggregate constructors accept validated components. Borrowed `From`
  conversions widen single-space, graph, and molecule remappings to
  correspondences; no owned conversion or inverse narrowing was added.
  Constructor callers and fixtures are migrated, including the pushout
  transport test's direct use of its correspondence. Exact rejection and
  widening cases pass, as does a sorted-reference check of all 701 image
  vectors of lengths zero through four with images in `0..=length`.
  Core/IR library tests passed (7,499 passed, 3 ignored), including the four
  production regressions; core/IR properties with `PROPTEST_CASES=256`
  (487 passed, 1 ignored); graph/IO library tests (4,297 passed); workspace
  all-targets check with Python 3.13.15; core/IR all-targets clippy with
  `-D warnings`; nightly formatting; and `git diff --check`.
  Molecule/span operation signatures and canonicalization remain for S2b.
- [x] **S2b — Molecule/span renumbering and canonicalization.** Graph-IR
  molecule remap, constraints, reaction span, canonicalization, traits, and
  dependent bindings. Breaking (red→green). [dep: S2a, S1e]
  Accept `MoleculeRemapping` in renumbering and `framed_eq_under`, and return
  it from `canonicalize_with_remapping`. Canonicalization constructs it
  directly. Route genuine constraint renumbering through the approved
  remapping-to-correspondence conversion and S1e's traversal, not a new lookup
  adapter. Migrate consumers, benchmarks, and tests together. Check exact
  identity/inverse/composition, all-family reference transport, framed
  equivalence, and equality of bare/witnessed canonical outputs.
  Completed 2026-09-04: molecule/span renumbering and molecule framed
  comparison accept `MoleculeRemapping`, checking all eight component counts.
  Molecule, reaction, and span canonicalization construct and return dense
  remappings directly; constraint transport widens to correspondence.
  `Remapping::remap_vec` and `try_remap_vec` own vector transport without a
  `Clone` bound. Python exposes the approved frozen `Remapping` and
  `MoleculeRemapping` types and `canonicalize_with_remapping` methods.
  Consumers, benchmarks, and normative documentation are migrated.
  Inverse test references sort source/target pairs independently of
  `remap_vec`; composition references apply the two image lookups directly.
  Existing generated domains and exact assertions are preserved.
  Core/IR library tests passed (7,527 passed, 3 ignored), as did integration
  and doc tests. Core/IR properties passed with `PROPTEST_CASES=256`
  (488 passed, 1 ignored), including the independent vector-remapping reference.
  The corrected IR property suite and affected canonicalization unit test
  were rerun successfully. Graph/IO library tests passed (4,297);
  rebuilt Python tests passed (1,350 passed, 2 skipped, Python 3.13.15).
  Workspace all-targets check and clippy, core/IR property-enabled clippy
  with `-D warnings`, nightly formatting, and `git diff --check` passed.

### S3 — Complete compaction carriers

- [x] **S3a — Single-space and graph compaction.** Graph-core compact, graph,
  and relation modules. Breaking (red→green). [dep: S0a, S1a]
  Add checked source-count construction, declared-count identity, and count
  accessors. Build graph aggregates from valid components; remove unbounded
  empty/default identities. Migrate graph and relation producers, including
  affected graph-IR construction sites, in this subitem. Add infallible
  conversions to single-space and graph correspondences. Test out-of-range
  removals, duplicate/order handling, empty and fully removed domains,
  survivor ordering, and compact/uncompact roundtrips.
  Include the dependent molecule component constructor and remove `relations`,
  `empty`, and raw-vector construction. Migrate editor and transaction
  producers with all eight pre-removal counts, preserving them through
  `UndoCompaction`. These construction changes are required for S3a to compile.
  Completed 2026-09-04: single-space compactions carry source counts, reject
  out-of-range removals, and expose bounded compact/uncompact and vector
  operations with the approved checked companions. Graph and molecule
  aggregates accept validated components; unbounded defaults and molecule
  `relations`/`empty` constructors are removed. Graph, relation-set, editor,
  and transaction producers preserve pre-removal counts in every component.
  Borrowed single-space and graph correspondence conversions preserve the
  complete survivor pairings and both counts. Rollback checks all eight
  result counts before inverse transport and reports `RollbackStateMismatch`
  on disagreement. Tests state the error behavior without historical notes.
  Core/IR library tests passed (7,572 passed, 3 ignored), plus integration
  and doc tests. Core/IR properties passed at `PROPTEST_CASES=256`
  (489 passed, 1 ignored), including independent survivor enumeration.
  Graph/IO library tests passed (4,297); rebuilt Python tests passed
  (1,350 passed, 2 skipped, Python 3.13.15). Workspace all-targets check,
  core/IR property-enabled all-targets clippy with `-D warnings`, nightly
  formatting, and `git diff --check` passed.
- [ ] **S3b — Molecule compaction.** Graph-IR compact, editor, transaction,
  and constraint consumers. Breaking (red→green). [dep: S3a]
  Add infallible molecule-correspondence conversion over the count-bearing
  components introduced in S3a. Test cascaded removals and exact all-family pairings against
  the source/result tables, plus existing rollback laws.

### S4 — Graph-core output/witness separation

- [ ] **S4a — Removal and relation compaction.** Graph-core graph and relation
  modules. Breaking (red→green). [dep: S3a]
  Give graph removal methods and relation-set compaction plain output-only
  forms and `*_with_compaction` companions. Preserve cascading versus
  dangling-condition behavior. Migrate editor and rewriting callers that
  need the compaction. Test identical resulting graphs/sets and failures for
  each pair, plus the witness's expected survivor images.
- [ ] **S4b — Pushout families.** Graph-core rewriting and relation modules.
  Breaking (red→green). [dep: S1a, S1b]
  Replace `Pushout` and `RelationPushout<S>` output containers with the
  accepted graph and relation pushout correspondence types. Bare methods
  return the object; witnessed methods return `(object, correspondence)`.
  Preserve named left/right mappings and their common result counts.
  Migrate graph-IR gluing/composition callers. Test both input mappings and
  bare/witnessed output equivalence, including coincidences.
- [ ] **S4c — Pullback and pushout complement.** Graph-core rewriting.
  Breaking (red→green). [dep: S1a]
  Introduce the accepted `PullbackCorrespondence` and
  `PushoutComplementCorrespondence` result separation and method pairs.
  Preserve categorical directions and existing admissibility conditions.
  Migrate reaction composition and rewriting tests; check the defining
  commutative mappings and output equivalence. Keep subdivision's distinct
  graph representation and its cross-entity source accessors intact.

### S5 — Graph-IR molecule operation returns

- [ ] **S5a — Combination and split.** Graph-IR molecule and fragment callers.
  Breaking (red→green). [dep: S1f]
  Make combination methods witness-free and state append order. Make split
  return components, with `split_with_correspondence` returning component/
  source-to-component pairs. Migrate Rust/Python callers, including `React`
  implementation plumbing. Test empty inputs, multiple components, all
  entity families, component ordering, append layout, and exact output
  equivalence. Keep chemistry-layer public signatures unchanged.
- [ ] **S5b — Molecule pushout.** Graph-IR molecule pushout and composition
  callers. Breaking (red→green). [dep: S1f, S4b]
  Replace `MoleculePushout` with `MoleculePushoutCorrespondence` and the
  accepted method pair. Keep left/right correspondence access independent
  of the molecule and preserve common result counts. Test overlay meets,
  stereo frames, constraints, inadmissible input, and output equivalence.
- [ ] **S5c — Removal and extraction.** Graph-IR editor removal families and
  molecule extraction. Breaking (red→green). [dep: S3b, S4a]
  Add the `with_compaction` variants and make plain methods witness-free.
  Preserve selection order, cascading behavior, constraints, and transaction
  use of removal. Test exact all-family compactions and equality of resulting
  state across each method pair. Adapt callers rather than adding hydrogen
  or other chemistry transformer variants.

### S6 — Editor and transaction correspondences

- [ ] **S6a — Editor session witness.** Graph-IR molecule editor.
  Additive (green). [dep: S1a, S3b, S5c]
  Track source-to-current pairings across direct additions, removals,
  modifications, and restoration. Expose witnessed snapshot and checked/
  asserted build companions using existing publication boundaries. The
  correspondence may describe transient editor id spaces; it does not
  require publication of an intermediate molecule. Test multi-step sessions,
  all entity families, repeated snapshots, and publication failures. Compare
  session pairings with composed operation pairings.
- [ ] **S6b — Batch application and transactions.** Graph-IR molecule apply,
  editor transact/apply, and transaction rollback. Additive (green).
  [dep: S6a] Add correspondence companions retaining the ordinary molecule,
  editor, transaction, or unit output. Forward witnesses cover the particular
  operation; session witnesses cover the editor's source. Rollback returns
  the inverse operation correspondence. Preserve transaction append and
  rollback behavior. Test failures, appended batches, add/remove sequences,
  atom provenance, inverse direction, and bare/witnessed output equivalence.
  Do not impose operational persistence on induced non-atom pairings.

### S7 — Reaction application results

- [ ] **S7a — Application at a supplied match.** Graph-IR reaction application.
  Breaking (red→green). [dep: S6b, S5b]
  Implement the four accepted `apply_at` result forms. Keep host-to-product
  correspondence direction, existing preconditions, and stereo/frame
  handling. Migrate direct consumers of the old derivation return. Test
  identical bare/witnessed products and failures, realized reaction/span
  consistency, and expected atom pairings.
- [ ] **S7b — Iteration and derivation retirement.** Graph-IR reaction
  iterators and their Rust/Python consumers. Breaking (red→green).
  [dep: S7a, S5a] Implement the four accepted iterative method forms, preserving
  captured-input ownership, match order, lazy result production, and terminal
  error behavior. Remove `ReactionDerivation`, its exports and bindings, and
  update all consumers. Keep `React` methods separate and verify their
  products still match explicit combine/application/split. No new result
  container duplicates the returned molecule/reaction/span.

### S8 — Python parity and integrated verification

- [ ] **S8a — Python exposure.** umol-py correspondence, molecule, edit,
  reaction, exports, and API inventory. Additive/breaking (red→green).
  [dep: S2b, S3b, S4c, S5a, S5b, S5c, S6b, S7b]
  Complete bindings for the accepted witness types, conversions,
  and operation pairs beyond earlier compile adaptations. Mirror Rust names,
  construction constraints, errors, and tuple/iterator shapes. Rebuild the
  native extension before tests. Verify direct calls outside notebooks,
  composition, bare/witnessed equivalence, and absence of flag-dependent
  witness returns and retired types.
- [ ] **S8b — End-to-end laws and measurement.** Existing operation property
  suites and S0 benchmarks. Additive (green). [dep: S8a]
  Exercise mixed compaction/remapping/correspondence sequences and
  input-to-pushout-to-split composition. Preserve the existing semantic laws
  while updating their APIs. Compare final timings with S0 and report both
  bare-operation regression and witness cost for the measured inputs.
- [ ] **S8c — Documentation and closeout.** Normative development guides,
  public documentation, doc 218, and status index. Documentation (green).
  [dep: S8b] Reconcile the final exported API against this design, update
  obsolete normative descriptions, and remove retired names from live
  examples and inventories. Preserve historical closed discussion records.
  Run the final gates below; record results before marking the work complete.

### Verification gates and critical path

Each subitem runs focused tests for its changed behavior. Each stage runs
library and integration tests for graph-core and graph-IR, their explicitly
enabled property targets, affected downstream tests, and a workspace compile
check including bindings. Use `PROPTEST_CASES=256` for the stage property
runs. Python-related commands run with `umol-py/.venv` activated and Python
3.13 confirmed; rebuild the extension before any affected pytest suite.

Final gates:

```sh
cargo +nightly fmt --all
source umol-py/.venv/bin/activate
python --version
cargo test --workspace
cargo test -p umol-graph-core -p umol-graph-ir --features proptest --test property
cargo test -p umol-graph --features proptest --test property
cargo clippy --workspace --all-targets -- -D warnings
maturin develop --manifest-path umol-py/Cargo.toml
pytest -q umol-py/tests
git diff --check
```

Run these in the activated shell; separate tool invocations must reactivate
it. Run affected conformance targets explicitly where consumer changes reach
their code paths. Re-run S0 benchmark commands with the saved baseline and
include fixture sizes and absolute timings in the closeout.

Transport dependency path: baseline → graph participant transport → edit
handle separation and molecular reference transport → legacy removal and
four producer fixes → strict remapping/conversions → renumbering consumers.
The independent compaction path is baseline/composition → counted compaction
→ removal/extraction pairs → editor witnesses → reaction results. Output/
witness separation for pushout, pullback, and complement builds on the
correspondence foundation. These paths join at Python parity and integrated
verification; compaction is not a prerequisite for the remapping fix.
No implementation stage above is optional for this work's completion.
Chemistry-layer witness API design remains separate and is not a deferred
stage of this plan.
