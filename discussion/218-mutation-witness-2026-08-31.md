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
the consumer needs exists exactly at commit and is currently dropped.

## Current state

- The **editor** returns witnesses piecewise in Rust: `add_atom` and its
  siblings return real ids immediately, and `MoleculeEditor::remove` returns
  `MoleculeCompaction`. Nothing reaches Python, and a session of edits yields
  no single pre-to-post correspondence.
- The **`Edits` script path** returns none: `Molecule::apply(edits)` returns
  only the new molecule, although its implementation constructs a
  `MoleculeCorrespondence` internally to transport constraints and drops it.
  `New(n)` handles are necessarily unresolved before commit and are never
  resolved for the caller after it. The transaction path behaves the same.
- `combine` and `canonicalize_with_correspondence` already surrender their
  moves; reaction application carries a correspondence in its derivation.

## Settled design

**The witness is `MoleculeCorrespondence`, unchanged.** A pre-to-post partial
bijection over all eight entity families: removed entities are absent,
survivors map through the compaction. No composite witness type is
introduced, because additions are recoverable from the correspondence
itself:

> **Addition law.** The post-state entities without a pre-image — the
> correspondence's right-unmatched entities, taken in id order — are the
> script's surviving additions, in creation order.

This holds because additions append, compaction preserves order (so
additions remain above survivors throughout), and an entity added and then
removed in the same script never materializes. `Edits` already forbids the
operations that would break handle order (sorting, deduplication,
concatenation). A caller that removed some of its own `New(n)` handles knows
which additions survived, so `New(n)` resolution is caller-computable from
the correspondence alone. The law is stated as a `# Semantic properties`
entry and pinned by tests, including an add-then-remove case.

### Transformation classification

The same principle covers every transformation, in three classes:

- **Identity class** — id-preserving transformations (`Normalize`, `Reframe`,
  resolution, attribute-only edits): the witness is the identity, stated as a
  law rather than returned.
- **Compaction class** — pure removals and selections, where the move is
  order-preserving and dense: the witness is `MoleculeCompaction`, the
  stronger type whose shape carries that guarantee. `MoleculeEditor::remove`
  already conforms.
- **Correspondence class** — everything that may add or renumber: the witness
  is `MoleculeCorrespondence`.

Class membership is assigned per operation by its effect on the id space,
regardless of the owning layer, and the effect is judged per family. The
chemistry-layer transformers move exactly one family: the `Aromatizer` adds
aromatic-system entities and the `Kekulizer` removes them; every other
family is identity, and their bond-order and constraint changes are
attribute changes, invisible to the witness. The `Kekulizer`'s witness is
therefore exactly a `MoleculeCompaction` (empty in seven components), while
the `Aromatizer` needs the correspondence form, since compactions cannot
express additions. Doc 166's hydrogen transformations will split the same
way (explicitation adds atoms and bonds; collapse removes them). The
`Transformer` trait carries the variant as a required
`transform_with_correspondence` member: required rather than defaulted, so no
transformer can silently claim the identity witness, and uniform across
classes, so a transformer pipeline's total witness is the composition of its
steps — identity-class transformers return the identity witness explicitly,
and compaction-class transformers speak correspondence through the
conversion (their stronger compaction may additionally be exposed
inherently).

Two independent questions govern the correspondence class, and both are
answered uniformly. First, witness versus documented layout: no operation's
result id layout is part of its contract, so a witness is the only sanctioned
mapping between input and output ids — documenting the layout merely to make
the mapping caller-computable would freeze an implementation detail into the
contract. Second, unconditional versus variant: variant, everywhere. A caller
holding no external ids is indifferent to where entities landed, so no
signature bakes in the assumption that every caller wants the mapping.
`combine`, `combine_from`, and `split`
([142](142-join-split-2026-07-10.md)) currently return correspondences
unconditionally and migrate to the same variant pair as the rest of the
surface.

**Compaction-to-correspondence conversion.** `Compaction<Id>` stores only the
sorted removed-id set; it is deliberately size-free, with survivor images
computed arithmetically. The conversion to `MoleculeCorrespondence` therefore
takes the pre-state entity counts (or the pre-state molecule) as context and
is lossless given them: it enumerates the survivor pairs, and the original
compaction is recovered from the result's unmatched left ids. Composition
across witness classes happens in the correspondence vocabulary through this
conversion.

**Variant pairs, not unconditional returns**, for single-value
transformations, per the `canonicalize` / `canonicalize_with_correspondence`
precedent: most call sites have no external id holder, and the plain
operation's shape should not bake that assumption in.

- `Molecule::apply` / `Molecule::apply_with_correspondence`.
- The transaction path gains the equivalent variant.
- The editor gains **both** forms: a finalization variant returning the
  session correspondence (the composition of the session's steps), and the
  piecewise route — `MoleculeCompaction` becomes Python-visible with a
  conversion to `MoleculeCorrespondence`, so callers may compose step
  witnesses themselves.

**Laws.**

- Composition: the witness of a sequence of mutations is the composition of
  the step witnesses; the editor's session correspondence equals the
  composition of its piecewise witnesses.
- Inversion: rollback's witness is the inverse partial bijection; entities
  created by the rolled-back mutation are absent from it.
- Identity: an attribute-only mutation yields the identity correspondence,
  so consumers need no moved-or-not case split.

**Python parity and shape.** Python catches up to the current Rust surface:
`MoleculeCompaction` (recently redesigned; the bindings are behind) and the
variant pairs, mirrored by name. The scikit-learn keyword shape
(`with_correspondence=True` selecting a richer return) was considered and
rejected: a flag-dependent return type needs `Literal` overloads to type
precisely and degrades to a union under a runtime flag, while method pairs
type exactly — decisive for a surface that maintains a signature inventory.
`canonicalize_with_correspondence` is already bound as a pair, so the idiom
is uniform. All eight families are exposed, as `MoleculeCorrespondence`
already provides.

## Alignment, not scope

- Doc [204](204-reaction-application-redesign-2026-08-19.md): the redesigned
  application result should speak this same vocabulary — a host-to-product
  `MoleculeCorrespondence` obeying the same addition law — so that consumers
  handle edit-apply and reaction-apply witnesses identically. That is an
  input to 204's design, decided there.
- Doc [212](212-remapping-layer-2026-08-26.md): the id-transport
  representation decision concerns the same data; the implementation should
  reuse the correspondence `apply` already constructs internally rather than
  adding a second computation, and 212's representation choice should be
  settled or explicitly deferred when this lands.

## Boundaries

- No stable ids, bound references, or identity infrastructure.
- Participant-frame positions are outside the witness: correspondences move
  entity ids, and frame-action transport remains the reframe and
  canonicalization machinery's concern.
- No unconditional witness returns on the plain operations.
- No editor API redesign beyond the added variants and bindings.

## Open items

- Whether the transaction's rollback returns the inverse witness directly or
  callers invert the forward witness.
- The exact reuse point of `apply`'s internal correspondence (implementation
  detail, settled during the plan).
