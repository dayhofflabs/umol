# 218 — Mutation witnesses for the editing surface

Status: Proposed
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
  consumers. `IdRemapping` remains the sparse transport used by constraints,
  deltas, molecule construction, and reaction-span operations.

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
  result entities. It is appropriate when selection or removal is the
  operation.

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
remapping. The existing generic and graph aggregate therefore need either
non-semantic map names or replacement by the appropriate correspondence or
operation-specific carrier. `MoleculeRemapping` becomes canonicalization's
witness only after its bijection contract is enforced.

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
| `Molecule::apply`; editor `apply`, `transact`, `snapshot`, and finalization; transaction rollback | one to one, arbitrary additions, removals, and modifications | molecule, editor, transaction, or `()` without an aggregate witness | Forward application needs a source-to-result correspondence accumulated from the executed edit identities. Rollback needs its inverse. Plain snapshot/finalization only needs the session witness when the editor originated from an existing molecule. |
| Aromaticity, stereo, and aggregate resolution | one to one, replacement-capable in overlay families | successful `Solution` without a witness | No witness-bearing form in this work. |
| `Transformer::transform` and `generate_all` | one source to one result, or to alternative results; effect depends on transformer | bare molecule values | No blanket witness method. Add an operation-specific form only when a transformer has nontrivial transport that callers need. |
| `canonicalize_with_correspondence`; `Molecule::remap` | one to one, dense bijective renumbering | source-to-canonical correspondence; caller-supplied correspondence respectively | Canonicalization's witness-bearing form returns `MoleculeRemapping`; `remap` returns no witness because the caller supplied it. |
| `Reaction::apply_at` and `Reaction::apply` | one host to one product per application | `ReactionDerivation` with host-to-product `comap` | Remove `ReactionDerivation`. Offer product-only, product-with-correspondence, realized-`Reaction`, and realized-`ReactionSpan` application forms. |
| `combine`, `combine_all`, and `combine_from` | one or many inputs to one result by disjoint append | operand-to-result or all input-to-result correspondences | Append order determines every input mapping. Return only the combined molecule, with preservation and ordering made contractual. |
| `meet_pushout` | two inputs to one result with possible identification | `MoleculePushout` containing both input-to-result correspondences | The correspondences are nontrivial. Keep a witness-bearing result; align the plain and `with_correspondence` forms with the naming rule. |
| `split` | one source to several simultaneous component results | one component-to-source correspondence per result | Plain `split` returns components. A witness-bearing form returns each component with its source-to-component correspondence. The collection supplies component identity. |
| `Fragment::attach` and fragment `Add` | two fragment bodies to one result body by append | fragment only; the internal `combine` correspondence is discarded | Append law; no witness-bearing form and no extra created-id return. `finish` alone is identity. |
| `React` for `Molecule` and `[Molecule]` | one or many reactants to alternative sets of product components | component molecules only | Keep the convenience operation product-only. Callers needing transport compose the explicit combine, reaction-application, and split operations. |

The audit also corrects an assumption about post-hoc reconstruction. An
operation witness records operational persistence, not merely structural
similarity between the finished molecules. An edit may remove an overlay and
add a new overlay with the same incidence. `MoleculeCorrespondence::induce`
matches bonds and overlays from mapped atom incidence, so it would describe
that replacement as preservation. Editor, resolver, and transformer witnesses
must therefore be accumulated from the executed operations rather than
inferred from the result.

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

The plain application name returns the product, while the correspondence form
uses the `with_correspondence` suffix. The explicitly named reaction and span
forms are different primary results, not witness variants of the product
method. The plural application iterators follow the same item contracts rather
than wrapping each item in a replacement for `ReactionDerivation`. Rust and
Python expose the same distinctions. `React` remains a product-only
convenience API.

The correspondence form must also correct a defect in the current
implementation. Reaction application constructs atom persistence
operationally and then calls `MoleculeCorrespondence::induce` for the other
entity families. Incidence-equivalent removal and addition can therefore be
misreported as persistence. The new host-to-product correspondence must be
accumulated from the executed operation identities rather than reconstructed
from structural similarity.

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

The general composition vocabulary is correspondence. A compaction can be
lowered losslessly when supplied the pre-state counts: enumerate its
survivors as matched pairs, use removed ids as the unmatched left side, and
set the right count to the survivor count. A semantic remapping can be
lowered directly to a total correspondence with equal carrier counts.
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

`MoleculeCompaction` deliberately stores only the compact arithmetic and not
the source counts. Its lowering therefore needs the source molecule (or
equivalent counts) and can fail if the objects disagree. Ordinary `From`
cannot express that context. A tuple-shaped `TryFrom` could, but obscures the
operation; prefer a named fallible conversion on either
`MoleculeCompaction` or `MoleculeCorrespondence`. Which type owns that method
remains an API-shape decision. Do not add counts to `MoleculeCompaction`
solely to make the conversion context-free.

The composition methods themselves and the replacement of non-semantic
`IdRemapping` consumers remain to be designed.

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

`split` likewise has a plain component result and a correspondence-bearing
form with covariant source-to-component mappings. `meet_pushout` retains its
nontrivial input-to-result correspondences in a witness-bearing form; the
exact plain return and result-type names remain to be settled.

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
the remapping-bearing name and return type. Operation-specific multi-object
results must be settled in Rust before they are bound.

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

## Open items

- Define lowering and composition among `MoleculeRemapping`,
  `MoleculeCompaction`, append results, and `MoleculeCorrespondence`, including
  the named contextual conversion from compaction to correspondence.
- Decide how semantic bijectivity is established by `Remapping` and
  `MoleculeRemapping`, and then classify every remaining `IdRemapping`
  producer before retiring the sparse carrier.
- Record and test existing-id preservation and allocation order for primitive
  additions, `combine`, `combine_all`, `combine_from`, and fragment assembly;
  add the covariant witness-bearing `split` form.
- Determine how `Molecule::apply`, editor finalization, and transactions
  accumulate a session-level source-to-result correspondence; none exists in
  the current implementation.
- Decide whether transaction rollback returns the inverse witness directly or
  callers invert the forward witness.
- Define witness-bearing variants for direct relation removals and extraction,
  and add operation-specific transformer variants only where the audit shows
  nontrivial useful information.
- Settle the plain and witness-bearing `meet_pushout` result shapes.
- Replace `ReactionDerivation` in the Rust and Python application iterators
  with the direct product, product-with-correspondence, `Reaction`, and
  `ReactionSpan` forms.
- Correct reaction application's non-atom witness construction so that
  remove/add replacement is not inferred as preservation from incidence.
- Bind the settled Rust surface in Python only after these contracts are
  fixed.
