# 179 — Expose editing and transactions on the Python surface

Status: Completed
Date: 2026-08-02
Relates: [178](178-python-lattice-ops-2026-08-01.md),
[043](043-mutative-undoable-mutation-2025-12-23.md),
[199](199-open-container-integrity-2026-08-18.md)

Doc 178 exposed the lattice operations because a central part of the model was unreachable from the
interface most users have. The same problem applies to molecule mutation: Python exposes reaction
`Delta` values, but it cannot construct and apply the host-specific `Edit` values used to modify a
particular molecule.

The distinction is semantic, not merely one of API layers. An `Edit` names entities in a host
molecule, may refer to entities created earlier in the same edit sequence, and carries the old-state
preconditions required for checked application. It is a reified mutation — an in-memory data value
rather than an operation — and is replayable against that host state, but it is not a portable
transformation that can be applied to an arbitrary molecule. A reaction instead combines a pattern
with `Deltas`; applying it to a host requires a match that anchors the pattern's entities in that host.
Construction can also be expressed as a pattern plus deltas, but it need not be understood as a
reaction: an ordered edit sequence can directly describe an independent generation step such as
attaching a methyl group to ammonia.

Earlier discussions used *serializable* loosely for this reified, data-form property. Actual external
serialization is a separate capability and is also required here: `Edits` parses and renders through
the existing EDN DSL. This does not introduce serde, pickle support, or another serialization
framework.

## Justification

Two reasons, and the second is the stronger.

**The whitepaper.** Section 9 (Mutation and Delta Encoding) describes both host-specific editing and
pattern-anchored transformation. Its listings are Python, like every other section. Without this work,
the section cannot show construction as a sequence of serializable edits or distinguish that sequence
from the pattern-plus-deltas representation of a reaction. The editing vocabulary would instead have
to be described in prose or shown in Rust — either of which makes the section the odd one out.

**Transactionality is a claim the model makes, not a convenience.** A batch of edits either applies
wholly or not at all. `Edit` is the caller-facing requested mutation vocabulary; checked application
records the realized `Undo` journal needed for exact rollback and restores that journal immediately if
an edit fails. An operational property that cannot be exercised from the primary interface is a claim
the reader has to take on trust.

The Python surface is additive. The Rust work makes `Edits` the sole semantic representation of an
ordered edit batch and moves every producer and transaction boundary away from interchangeable
`Vec<Edit>` values.

## Scope

**In:**

- `Edits`, the ordered container corresponding to `Deltas`. Order is semantic because later edits
  may refer to entities created by earlier edits; the container does not canonicalize or deduplicate.
  It is also the complete construction surface: additions, checked updates, batched removals,
  topology removal, and molecule-constraint changes append directly to the container.
- A full Rust migration from semantic `Vec<Edit>` accumulators to `Edits`, including transaction
  entry points, resolver plans, molecule operations, reaction application, DSL resolution, tests,
  and property generators. `Edit` remains the public raw entry enum, but production code constructs
  executable sequences through `Edits`.
- EDN `parse` and `render` operations for the edit sequence, using the existing DSL infrastructure.
- A normative standalone-edit grammar in `umol-ast/spec/umol-dsl-spec.md`, adjacent to the existing
  reaction-map and delta grammar in §8. The specification must define the ordered vector shape,
  handle syntax, checked `:expect` / `:update` modification form, removal preconditions, and batch
  semantics rather than leaving them only in this implementation discussion.
- Entity handles (`AtomHandle`, `BondHandle`, and the rest), resolved either as an `Id` naming an
  entity in the transaction's initial host or as a same-sequence `New` creation ordinal. Python uses
  bare non-negative integers for initial-host ids and one generic immutable `New` wrapper for
  creation ordinals; it does not expose eight typed handle classes.
- Stable handles for initial-host and same-batch entities, with an invalid `Id(n)` or an unissued or
  forward `New(n)` reported through `TransactionError::HandleOutOfRange`. Without creation handles an
  edit batch cannot bond an atom it just added, which is most of what a batch is for.
- `MoleculeAst::apply(edits)` as the ordinary immutable-style operation. It creates an editor,
  applies the edits atomically, consumes the editor to build the result, and leaves the source
  molecule unchanged.
- `MoleculeEditor`, `transact`, consuming `build`, and non-consuming `snapshot`, together with the
  `Transaction` returned by `transact` and its `rollback` operation.
- The `*Update` types the `Edits::update_*` methods take.
- `TransactionError` as the Python exception for application and rollback failures. DSL failures
  remain `ParseError`; use of a consumed editor or transaction raises `RuntimeError`.

**Out, unless the implementer finds a reason:**

- `transact_unchecked` — the unchecked path is a Rust-internal optimisation and exposing it invites
  callers to skip the validation that makes the facility worth having.
- `Transaction::append`.
- Identifier compaction and remapping (`IdCompaction`, `UndoCompaction`, `remap_delta`). Internal
  bookkeeping; no reader-facing question depends on it.
- A parallel serde, JSON, pickle, or other serialization framework. The external representation is
  the umol EDN DSL.
- `Vec<Edit>` overloads, `IntoIterator<Item = Edit>` transaction entry points, or other parallel
  batch representations. Raw `Edit` values may be inspected and appended, but executable batches
  are `Edits`.

## Settled semantics

- A failed `transact` leaves the structure untouched. This is the property the paper asserts, so it
  is the property the binding must not weaken. Application remains optimistic: it mutates the live
  editor, records each realized `Undo`, and replays that journal if a later edit fails. An undo
  generated by the same application must always succeed; `RollbackFailed` remains the defensive
  error for a broken internal invariant, and generated valid-prefix/failing-tail tests must establish
  that it is unreachable in ordinary checked application.
- `MoleculeAst::apply` is the ordinary path when no intermediate editor state or explicit rollback
  is needed. It returns the edited molecule and retains the original.
- `MoleculeEditor::build` is finalization and consumes the editor. `snapshot` materializes the
  editor's current state without consuming it; it is reserved for inspection followed by further
  editing or rollback, not required in the ordinary construction flow. This use of *snapshot* follows
  the repository nomenclature in [doc 177](177-nomenclature-guide-2026-07-31.md).
- The returned `Transaction` is the detached journal defined by the Rust API, not a borrowed guard.
  Python mirrors `transaction.rollback(editor)`. Rollback consumes the transaction; a second call
  raises `RuntimeError`. Applying the transaction returned by successful application to the
  corresponding post-application editor, or applying an appended journal to the end state of that
  consecutive transaction chain, restores the original molecule exactly. Python does not expose
  `append` in this round.
- Rust `Transaction` remains `Clone`: cloning a detached journal is valid, and applying a clone to
  an unrelated state is governed by the rollback guarantee boundary below. Python deliberately does
  not expose that cloning surface and retains consuming wrapper semantics.
- A detached journal is not bound to an editor. Applying an unrelated transaction to an editor has
  no semantic result guarantee and may return a transaction error after partial mutation, but it
  must not panic. Undo application therefore checks the bounds and structural assumptions required
  for panic-free reconstruction and returns `TransactionError::RollbackStateMismatch` instead of
  indexing or unwrapping invalid state; it does not clone the molecule, retain an application
  snapshot, compensate a failed rollback, or attempt to prove that the journal and editor share
  provenance. This follows the same deliberate tradeoff as accepting a valid `AtomId` from one
  molecule when indexing another.
  `Transaction` has no public constructor from arbitrary `Undo` values. Successful application is
  its ordinary source; a valid empty transaction remains available for accumulation.
- Handles accept the Rust `Id` and `New` forms. `Id(n)` stably names entity `n` of the handle's kind
  in the transaction's initial host, irrespective of later dense-id compaction. `New(n)` stably names
  the nth entity of that kind created earlier in the edit sequence: `AtomHandle::New(0)` and
  `BondHandle::New(0)` may both occur in one sequence. The handle type supplies the entity kind, so
  no cross-kind type-mismatch state or `RefTypeMismatch` error is needed.
- Python handle arguments accept `int | New`. `0` denotes entity zero in the initial host;
  `New(0)` denotes the first entity of the relevant kind created by this edit sequence. The argument
  position supplies the kind, as it does for a bare integer.
- `TransactionError::HandleOutOfRange { kind, index, count }` means that an `Id(index)` did not exist
  in the initial host or a typed `New(index)` ordinal has not yet been issued for that entity kind;
  `count` is the size of the relevant namespace. `HandleRemoved { kind, index }` means that the
  handle was valid but its entity was subsequently removed directly or by cascade. The error does
  not repeat whether the caller used `Id` or `New`; the caller already owns that handle. The current
  `IdOutOfRange` variant is subsumed by these handle-resolution errors.
- A topology or overlay-removal batch rejects a repeated entity before mutating the editor with
  `TransactionError::DuplicateRemoval { kind: EntityKind }`. This is an operation precondition, not
  container deduplication: `Edits` continues to preserve repeated entries and their order. Add
  `Display` for `EntityKind` using the established entity names so the error remains typed without
  exposing Rust variant spelling in its message.
- The existing `Edit::for_*_update` functions move to mutating `Edits::update_*` methods without
  replacements on `Edit`. Their signature is `(&mut self, id, current_entity_ast, update)` and they
  append zero or more checked field/constraint edits to the existing sequence. The current entity
  AST supplies the checked old values; no molecule context is required. Overlay-removal methods take
  complete batches of recorded entity handles, participant handles, and entity ASTs. A stereo-atom
  site is an atom handle, a stereo-bond site is a bond handle, and each ligand carries an atom handle
  plus its `StereoLigandKind`; Python represents every such handle as `int | New`.

### Rollback guarantee boundary

The transaction API makes three deliberately different guarantees:

1. If checked application rejects an edit, its internally generated journal restores the exact
   pre-application molecule.
2. Applying a successful transaction to its corresponding post-application state restores the exact
   pre-application molecule. The same holds for an appended journal and the end state of its
   consecutive transaction chain.
3. Applying a transaction to any other editor must not panic. A structurally incompatible state
   returns `RollbackStateMismatch`; a structurally plausible but unrelated state may accept the
   journal and mutate unrelated entities. Its result is deliberately unspecified.

The third guarantee does not attempt provenance checking. Entity handles and realized undo ids use
the repository's ordinary dense integer newtypes, so a numerically valid id can name an unrelated
entity just as an `AtomId` from one molecule can name an atom in another. Preventing this would
require binding the journal to an editor, persistent lineage/version state, or retaining a molecule
snapshot; none is justified for this facility. Rustdoc on `Transaction`, `Transaction::rollback`,
and `RollbackStateMismatch` must state this boundary directly.

## Standalone edit DSL

An `Edits` document is a bare ordered vector. Individual edits use the singular entity keys already
used by reaction deltas. Existing host entities are written as positional integers and creations as
`{:new n}` in every position that accepts a handle. For example:

```clojure
[{:atom {:add "C#h3"}}
 {:bond {:add [0 {:new 0} "1"]}}]
```

A checked modification keeps the entity handle positional and records its precondition separately
from its update:

```clojure
{:atom
 {:modify [0 {:expect "#h3"
              :update "#h2"}]}}
```

Both values use the corresponding partial entity DSL. They must address the same field or constraint
keys. An undetermined constraint in either position denotes absence, so `:expect "#v*"` followed by
`:update "#v4"` adds a valence constraint.

Singleton additions and removals use singular forms. Simultaneous atom/bond removal is the distinct
`:topology` operation because atom removal can cascade incident bonds and overlays before a separately
recorded bond removal validates. The combined operation validates and captures one pre-removal state
and compacts once.
Batched overlay removals use the corresponding plural entity form and retain each removed entity's
participants and AST. These forms preserve the batching that affects edit semantics; implementation
batching that does not affect semantics need not appear in the surface representation.

`Edits.parse(text, *, defaults=None)` and `edits.render(*, defaults=None)` use `MoleculeDefaults`,
following molecule parsing. Defaults apply to full entity definitions and recorded removal state,
never to partial `:expect` / `:update` values.

## Handle allocation and resolution

Symbolic handles belong to construction of the edit sequence; concrete ids belong to its application.
The two responsibilities remain separate:

- `Edits` is the handle issuer and complete sequence builder. Its private storage is an ordered
  `Vec<Edit>`, accompanied by one private creation counter per entity kind. The singular and batched
  `add_*` operations append entries, advance the corresponding counts, and return the typed `New(n)`
  handles. The eight `update_*` operations append every derived field/constraint edit. Topology,
  batched overlay-removal, and molecule-constraint methods append their corresponding entries.
  Construction from an iterator and EDN parsing derive all eight counters in one pass.
- The creation counters assign stable `New(n)` ordinals only. They never contain host ids, do not
  change when an entity is removed, and do not determine whether a created entity remains alive.
  The `Edits` rustdoc must state both handle meanings and this ownership boundary explicitly:
  `Id(n)` names an initial-host entity, `New(n)` names a same-kind creation, and neither is rewritten
  inside the edit sequence when application compacts concrete ids.
- The mutation surface of `Edits` is append-only. It exposes iteration and appending operations that
  preserve already-issued handles, but not mutable iteration, insertion, removal, or reordering that
  could change the ordinal of an earlier creation. Two independently constructed `Edits` values are
  not naively concatenated: the second sequence's `New` handles would need per-kind rebasing, and its
  host ids may have been invalidated by structural edits in the first sequence.
- `MoleculeEditor` allocates the concrete per-family ids while applying the sequence.
- For each entity kind, the transaction applicator maintains an initial-host table indexed by
  `Id(n)` and a creation table indexed by `New(n)`. Each entry is either the entity's current concrete
  id or a tombstone. Its module and application-state rustdoc must state that this host-specific
  realization, including liveness and compaction, belongs exclusively to transaction application.
- After every removal, the applicator applies its `IdCompaction` to both tables for all eight entity
  kinds. Surviving concrete ids are updated and entities removed directly or by cascade become
  tombstones without changing their handle identities. Resolving a tombstone is a transaction error
  rather than an alias to whatever entity later occupies the same dense id. Topology removal can
  affect all eight kinds because atom or bond removal can cascade through every overlay family;
  direct overlay removal affects its corresponding kind. The current implementation passes `Id`
  handles through unchanged and does not update the created-id list, so either handle form can retain
  a stale id or silently resolve to a different entity after compaction.
- The handle-realization tables are application state, not part of `Transaction` and not a public
  combined creation namespace. If cross-family creation chronology or a post-application result map
  becomes a real consumer need, it can be recorded independently without changing typed-handle
  semantics.

`Deltas` provides the repository precedent for this boundary. Public reaction, difference,
validation, and composition APIs use `Deltas`, not `Vec<Delta>`. Two remaining semantic internal
accumulators currently use `Vec<Delta>`: `EntityFold::deltas_from_states`/reaction-span lowering and
the reaction property generator. They should be migrated to `Deltas` in the same work so the stated
container rule holds internally as well as publicly. Vectors used only as table-test inputs or as
private scratch storage inside a container algorithm are not competing semantic representations.

The EDN form is correspondingly local to the surrounding entity family: an integer denotes initial-
host `Id(n)` and `{:new n}` denotes `New(n)`. No molecule metadata is required to interpret either
form.

Keyword ids and structural refs are not part of the first standalone-edit grammar. Unlike a reaction
map, an edit vector contains neither its host molecule nor the host's persistent DSL metadata.
Supporting those forms therefore requires a separate host-aware resolution design, including the
meaning of a structural ref after preceding edits have changed or compacted the molecule. The
initial grammar is deliberately extensible to those reference forms without accepting them yet.

## Verification

Follow 178: keep the algebraic property tests in Rust, and have the Python tests verify method
availability and representative cross-boundary results. Add at minimum a round trip in which a batch
is applied, the resulting structure is compared against the expected value, the transaction is rolled
back, and the structure is compared against the original.

Do not attempt to generate arbitrary `(MoleculeAst, Edits)` pairs. A generator that reproduces
checked old-state payloads, cascades, handle liveness, and dense compaction would duplicate too much
of the transaction implementation. Use narrow generated traces over small distinguishable fixtures:
choose an entity kind, initial and created counts, removal subsets, a surviving or removed target,
and the position of an invalid batched entry. Compute expected identity by filtering independently
labeled entities rather than by calling `IdCompaction` or transaction helpers.

The Rust property suite must state and verify the resulting laws: initial-host `Id` handles remain
stable across compaction; surviving `New` handles retain their identities; removed handles become
tombstones; creation ordinals are never reused; creation namespaces are independent by entity kind;
a failed batched edit is atomic; apply followed by rollback is identity even when later edits use
handles across compaction; checked and unchecked application agree for the same valid traces; and
collecting or parsing an edit sequence recovers the next per-kind creation ordinals exactly.

The transaction matrix is broader than those new handle laws. Prepared valid-edit strategies must
span every `Edit` and `Undo` family, including both stereo overlays, all eight field and constraint
receivers, molecule constraints, topology cascades, and the six direct overlay removals. For each
family, checked apply followed by rollback must restore the exact ordered AST, and checked and
unchecked application must agree on the valid input. A second strategy inserts a rejected edit after
a generated valid prefix and verifies the exact original editor, the primary application error, and
the absence of `RollbackFailed`; vary the rejected member across the first, middle, and last position
of every batched edit.

Constraint-compaction properties must cover removed and remapped references for every entity kind,
including mixed updates, duplicate values, and exact order restoration. Molecule-constraint tests
must likewise exercise duplicate values: addition appends, removal selects the last matching value,
and rollback restores the original multiset order. Correct the existing comment that says removal
selects the first match.

`Transaction::append` composes journals from transactions applied consecutively to successive
states; it does not imply that an arbitrary `Edits` batch can be split into independently resolved
single-edit batches. Replace the existing over-broad materialization property with that journal law,
including empty-journal identity and exact reverse rollback order. At the resolver boundary, retain
successful multi-stage materialization tests and cover underdetermined, contradictory, and error
exits after an earlier stage has committed, always restoring the original molecule exactly.

Use direct unit tables for the cases that are cumbersome or less informative to generate: exact
`HandleOutOfRange` and `HandleRemoved` payloads for both handle forms; a three-entity initial-host
alias regression; add-two/remove-first/use-second and add-remove-add concrete-id reuse; topology
cascade tombstones for every overlay family; interleaved `New(0)` handles for all eight kinds; and an
invalid first, middle, or last member of `AddBonds` and each batched overlay removal with exact
failure atomicity; repeated atom, bond, and overlay removals with exact `DuplicateRemoval` values;
empty transaction and empty edit-batch identities; every currently uncovered field/constraint
receiver; exact molecule-constraint duplicate ordering; and cross-application of journals to
unrelated editor shapes to establish that every undo path returns normally or with an error rather
than panicking. The state produced by those deliberately mismatched applications is not asserted.

Separate the mismatched-journal tests into two cases. A structurally impossible journal/editor pair
asserts the exact `RollbackStateMismatch` value. A structurally plausible but semantically unrelated
pair calls rollback directly and asserts only by returning from the test without panic; it must not
compare the result or resulting editor. Pinning a particular wrong molecule would turn accidental
behavior into an API contract and obstruct later validation improvements.

Do not add artificial tests merely to mention error variants. Remove the unused `KindShapeMismatch`
and `DuplicateEntry`; add the typed `DuplicateRemoval` for the actual repeated-removal failure.
`RollbackFailed` remains a defensive internal-invariant error and is tested by establishing its
absence from valid-prefix/failing-tail application, not by constructing an invalid public journal.

The DSL implementation must also test the normative grammar directly: parse/render round trips for
all entity families, existing and `New` handles, checked modifications, removal preconditions,
constraint edits, and the simultaneous-removal forms whose semantics cannot be reproduced by a
sequence of independently applied removals.

An inventory check must establish that no public Rust signature accepts or returns `Vec<Edit>` and
that no production or property-generation accumulator representing an executable batch uses it.
The corresponding check for `Vec<Delta>` must leave only test-case data and private scratch storage
inside `Deltas` algorithms.
The reader-facing check is Section 9's listings, which must execute against the built module before
the section ships.

## Staged implementation plan

### S0 — Establish the semantic containers

- **S0a — Complete the `Deltas` precedent internally.** Change the crate-private entity-fold
  recovery helper to append directly to a `Deltas` accumulator, make
  `ReactionSpanAst::to_reaction` and the reaction property generator construct `Deltas` directly,
  and leave vectors used only for table-test parameters or scratch storage inside
  `Deltas::canonicalize` alone. Test reaction-span recovery and the generated reaction invariants
  against exact `Deltas` values, then verify that no production or property-generation accumulator
  representing a delta collection uses `Vec<Delta>`. **Breaking (red→green).** [dep: none]
  **Done.**
- **S0b — Add the `Edits` container and structural construction surface.** In
  `umol-ast/src/ast/edit.rs`, add `Edits` with private ordered `Vec<Edit>` storage and eight private
  creation counters. Provide `new`, immutable slice/length/iteration access, consuming iteration,
  `FromIterator<Edit>`, and raw single-entry `push`; update counters on `push` and derive them in one
  pass for collected entries. Add singular and batched addition methods returning every allocated
  typed `New(n)` handle, simultaneous topology removal, the six batched overlay-removal methods, and
  molecule-constraint add/remove. Do not add mutable iteration, insertion, removal, reordering,
  sequence concatenation, or an exposed `Vec<Edit>`. Re-export `Edits` from `ast.rs`. Test exact
  emitted entries, interleaved creations, batched atom/bond additions, independent ordinals for all
  eight kinds, complete removal preconditions, and counter recovery through `FromIterator`. Its
  rustdoc must define `Id(n)` as a stable initial-host identity and `New(n)` as a stable per-kind
  creation ordinal, while concrete ids, liveness, and compaction belong to transaction application.
  Add a property over arbitrary interleavings of the eight addition families: collecting or pushing
  the same entries must recover the same next per-kind `New(n)` values as uninterrupted construction.
  **Additive (green).** [dep: none]
  **Done.**

### S1 — Migrate every Rust edit producer and transaction boundary

- **S1a — Move checked update expansion onto `Edits`.** Move the bodies of
  all eight `Edit::for_*_update` functions to mutating `Edits::update_*` methods that append zero or
  more entries and return nothing; remove the old functions without aliases. Remove the trivial
  `Edit::add_atom`, `add_bond`, `remove_atom`, and `remove_bond` helpers so `Edit` remains the raw
  entry representation rather than a parallel construction API. Convert the focused `edit.rs` tests
  to the container API and assert for every entity family that `Edits::update_*` emits nothing for an
  empty update and that applying its output reproduces `current.update(update)`. Callers elsewhere in
  the workspace remain red until the rest of S1 completes; no permanent compatibility API is added.
  **Breaking (red within S1).** [dep: S0b]
  **Done.**
- **S1b — Move the transaction boundary to `Edits`.** Change `MoleculeEditor::transact` and
  `transact_unchecked` to accept concrete `Edits`, then migrate the transaction unit tests, edit
  property strategies, fixtures, examples, and benchmarks from vector batches and removed `Edit`
  helpers to `Edits`. Preserve order and duplicates exactly. A minor temporary internal adapter is
  acceptable if it makes this subitem independently verifiable, but it must be removed before S1i;
  the stage is otherwise allowed to remain red between subitems. **Breaking (red within S1).**
  [dep: S1a]
  **Done.**
- **S1c — Realize stable handles in per-kind namespaces.** In
  `umol-ast/src/ast/molecule/transact.rs`, replace the combined created-entity list with an
  initial-host table and a created-entity table for each kind. Initialize the former from the host's
  starting id spaces and append concrete ids to the latter as additions succeed; both tables retain
  tombstones after removal. Replace the older handle errors with
  `HandleOutOfRange { kind: EntityKind, index: usize, count: usize }` and
  `HandleRemoved { kind: EntityKind, index: usize }`; remove `RefTypeMismatch` and `IdOutOfRange`.
  Resolve every typed `Id(n)` against its initial-host table and every `New(n)` against its same-kind
  creation table. Replace the older
  *ref* terminology in handle rustdoc at the same boundary. Apply every removal's `IdCompaction` to
  both tables for all eight kinds, updating surviving ids and retaining tombstones for direct and
  cascaded removals. Add module and application-state rustdoc stating the stable meaning of both
  handle variants and that transaction application alone owns their host-specific realization,
  liveness, and compaction. Resolve and validate every fallible input to a batched edit
  before mutating the editor; in particular, fix `AddBonds` so failure in a later entry cannot leave
  earlier entries from the same edit applied without an undo. Add `Display` for `EntityKind`, remove
  the unused `DuplicateEntry`, add `DuplicateRemoval { kind: EntityKind }`, and reject repeated
  resolved entities in topology and overlay-removal batches before mutation. Add exact transaction
  unit tables for both errors and handle forms; the three-entity initial-host alias regression;
  add-two/remove-first/use-second; add-remove-add concrete-id reuse; topology-cascade tombstones for
  each overlay family; interleaved per-kind `New(0)` handles; and invalid first, middle, and last
  entries in `AddBonds` and every batched overlay removal. Every failure case must compare the whole
  editor state with its pre-transaction value.
  **Breaking (red within S1).** [dep: S1b]
  **Done.**
- **S1d — Add stable-handle transaction properties.** In
  `umol-ast/tests/property/strategies.rs`, add constrained trace strategies over small fixtures whose
  entities carry independently distinguishable values. Generate entity kind, initial and created
  counts, removal subsets, surviving or removed targets, interleaved creation order, and invalid
  batch position; compute expected identity by filtering labeled entities rather than using
  `IdCompaction` or transaction helpers. In `tests/property/edit.rs`, verify stable initial-host `Id`
  handles, stable surviving `New` handles, tombstones, non-reused creation ordinals, independent
  per-kind namespaces, and failed-batch atomicity. Extend the existing valid transaction strategy so
  apply-rollback identity and checked-unchecked agreement include compaction followed by later handle
  use. **Additive tests (red within S1).** [dep: S1c]
  **Done.**
- **S1e — Make rollback structurally fallible.** Remove `Transaction::new(Vec<Undo>)` entirely and
  construct the transaction fields directly inside checked application; retain `Clone` and the valid
  empty transaction used for accumulation. Remove the
  unused `KindShapeMismatch`. Add `TransactionError::RollbackStateMismatch` and make every undo path
  structurally fallible: check the entity counts, indices, compaction dimensions, and reconstruction
  slots it requires before indexing or unwrapping, and return
  `TransactionError::RollbackStateMismatch` for incompatible state. Do not clone the editor, retain a
  post-application snapshot, generate compensating undos, or validate unaffected state. Correct
  molecule-constraint removal documentation to its last-match multiset semantics. Use focused unit
  tables to exercise successful application of every `Undo` family, structurally incompatible
  counts/indices/compactions/reconstruction slots, empty journals, duplicate removals, uncovered
  field and constraint receivers, and exact molecule-constraint duplicate ordering. Add rustdoc to
  `Transaction`, `rollback`, and `RollbackStateMismatch` stating the three-part guarantee boundary
  above. **Breaking implementation and unit tests (red within S1).** [dep: S1d]
  **Done.**
- **S1f — Complete the transaction-law properties.** Add prepared valid-edit strategies spanning
  every `Edit` and `Undo` family, including all field-change variants, all eight constraint receivers,
  both stereo overlay families, molecule constraints, direct removals, and topology cascades. Verify
  exact apply-rollback identity and checked-unchecked agreement for each family. Generate a valid
  prefix followed by a rejected edit and assert the primary error, absence of `RollbackFailed`, and
  exact original editor state. Add constraint-compaction properties for removed and remapped
  references of every kind with duplicates and order. Replace the current split-batch append property
  with the actual law: journals from transactions applied consecutively to successive states append
  and undo in reverse order, with the empty journal as identity. Cross-apply independently generated
  valid journals and editor states to verify that mismatched rollback never panics, without inspecting
  the result or resulting state. **Additive tests (red within S1).** [dep: S1e]
  **Done.**
- **S1g — Migrate molecule and reaction producers.** Change `MoleculeAst::edits` and every
  `ReactionAst` executable-batch accumulator to `Edits`. Reaction application must emit its globally
  scheduled additions, updates, removals, and constraints through the complete container surface
  without losing batching. Migrate the corresponding AST unit tests, property generators, fuzz
  targets, fixtures, examples, and benchmarks. **Breaking (red within S1).** [dep: S1f]
  **Done.**
- **S1h — Migrate valence and structurally narrow graph plans.** Change the valence, atom-typing,
  bond, and multicenter resolver/planner results from `Vec<Edit>`/`Solution<Vec<Edit>, _>` to
  `Edits`/`Solution<Edits, _>`, construct their batches through `Edits`, and migrate their tests and
  fixtures. **Breaking (red within S1).** [dep: S1g]
  **Done.**
- **S1i — Migrate aromaticity and stereo plans and close the inventory.** Change the aromaticity and
  stereo planner results and the top-level resolution coordinators to `Edits`, migrate their tests
  and fixtures, and remove any temporary migration adapters. At the top-level resolver boundary,
  verify exact restoration on underdetermined, contradictory, and error exits after one or more
  earlier stages have committed. Verify by source inventory that no
  public signature and no production, test fixture, or property generator representing an executable
  batch uses `Vec<Edit>`, then run the complete workspace tests to restore green. Direct `Edit`
  construction may remain only inside `Edits` and representation-oriented matching/rendering.
  **Breaking (red→green).** [dep: S1h]
  **Done.**
- **S1j — Add snapshot and immutable apply.** In `ast/molecule/editor.rs`, add
  `MoleculeEditor::snapshot(&self) -> MoleculeAst` without changing the consuming `build`. In
  `ast/molecule.rs`, add `MoleculeAst::apply(&self, edits: Edits) -> Result<MoleculeAst,
  TransactionError>` using the checked transaction path. Test that snapshots do not consume or
  detach subsequent editor changes, successful apply leaves the source unchanged, and failed apply
  returns the transaction error while leaving the source unchanged. **Additive (green).**
  [dep: S1i]
  **Done.**

### S2 — Constraint edit representation, standalone edit DSL, and normative grammar

- **S2a — Add the standalone handle codec.** Add `umol-ast/src/dsl/edit.rs` with the shared
  surface representation and typed conversions for existing non-negative integer handles and
  `{:new n}` handles. Reject negative indices, keywords, structural references, and malformed
  `:new` maps at this boundary. Add direct parse/render cases for both forms in every typed handle
  position and exact parse-error cases for the rejected forms. **Additive (green).** [dep: S1i]
  **Done.**
- **S2b — Encode atom, bond, topology, and molecule-constraint edits.** In `dsl/edit.rs`, implement
  the singular atom/bond add and remove forms, simultaneous `:topology` removal, checked
  `:expect`/`:update` modifications, entity-constraint changes, and molecule-constraint add/remove.
  Full additions and recorded removal state participate in `MoleculeDefaults`; partial expect and
  update values do not. Enforce that the expect and update sides address the same field or
  constraint key. Test defaults in both directions, constraint addition/removal via undetermined,
  and a topology removal whose meaning cannot be reproduced by sequential removal edits.
  The lowering path must append through `Edits` rather than constructing an intermediate vector.
  **Additive (green).** [dep: S1i, S2a]
  **Done.**
- **S2c — Encode dative, aromatic, multicenter, and noncovalent edits.** Add the corresponding
  forms, including singular additions, batched plural removals with participant and AST
  preconditions, checked field updates, and constraint updates. Test every operation for each family
  with both existing and `New` participant/entity handles and exact parse/render output. **Additive
  (green).** [dep: S1i, S2a]
  **Done.**
- **S2d — Encode stereo overlay edits.** Add stereo-atom and stereo-bond additions, batched
  removals, checked field updates, and constraint updates. Preserve the distinction between atom and
  bond sites and render each ligand as an atom handle plus `StereoLigandKind`. Test both site kinds,
  mixed existing/`New` ligand frames, removal preconditions, and every stereo update form.
  **Additive (green).** [dep: S1i, S2a]
  **Done.**
- **S2e — Add the normalized `ConstraintEdit` representation.** In
  `umol-ast/src/ast/edit.rs`, add the public `ConstraintEdit` value used only by molecule-constraint
  edit entries. Keep `Constraint`, `ConstraintDelta`, and the inline entity-constraint ASTs
  unchanged. Store one private `Constraint` tree whose target-molecule ids are normalized slots,
  together with eight private typed handle vectors (`AtomHandle` through `StereoBondHandle`); an id
  of kind `K` and index `i` in the private tree indexes handle vector `K[i]`. Intern repeated
  references once per kind, reuse the slot throughout logical combinators, and map only the target
  member of each `SubPatternAnchor` pair—the nested pattern and its ids retain their independent
  coordinate space. Provide an infallible conversion from a concrete `Constraint` that maps every
  referenced entity to its same-id `Id` handle, plus a checked public construction path from a
  complete, kind-correct set of per-entity handle mappings for reaction lowering, DSL parsing, and
  callers that need `New` references. Do not expose the normalized constraint or its slot ids as a
  second public constraint model. Test a single entity leaf, repeated references in nested
  `And`/`Or`/`Not`, explicit relational references, quantified relational predicates that carry no
  atom identity, molecule-wide constraints with `None` subsets, mixed-kind subsets, and anchored and
  unanchored subpatterns. Verify that equivalent repeated handles share one slot and that converting
  an ordinary constraint preserves every exact predicate and maps every outer reference to `Id`.
  **Additive (green).** [dep: S1i]
  **Done.**
- **S2f — Make constraint edits handle-aware end to end.** Change
  `Edit::{AddMoleculeConstraint, RemoveMoleculeConstraint}` to carry `ConstraintEdit`, and change
  `Edits::add_molecule_constraint` / `remove_molecule_constraint` to accept and append that concrete
  value. Keep the ordinary concrete-`Constraint` path available through an explicit infallible
  conversion at the call site. In
  `ast/molecule/transact.rs`, resolve all eight handle vectors through `ApplicationState` before
  mutation, construct the total `IdRemapping` required by the normalized tree, and materialize an
  ordinary concrete `Constraint` through `Constraint::remap`; a forward, out-of-range, or removed
  handle must return the existing transaction error before the constraint list changes. Store only
  the realized concrete constraint in the undo journal so existing constraint compaction and
  reverse-order rollback remain authoritative. In reaction application, build the complete
  reaction-entity-to-edit-handle mapping while scheduled additions issue their actual per-kind
  `New(n)` handles: preserved LHS entities map to matched host `Id` handles and created entities map
  to those issued `New` handles. Lower every `ConstraintDelta` through that mapping, retaining the
  established additions → updates → constraints → removals schedule. Correct `remap_delta` so both
  constraint-delta variants remap their contained `Constraint` rather than passing it through.
  Migrate all raw edit constructors and tests in the same subitem. Add exact application and
  rollback cases for initial and created references of every kind; mixed relational references;
  compaction before later constraint addition/removal; duplicate-value last-match removal; removed
  and forward handles; target-only subpattern-anchor remapping; constraint deltas over created atoms
  and overlays; and direct `remap_delta` coverage. **Breaking (red→green).** [dep: S2e]
  **Done.**
- **S2g — Extend molecule-constraint edit parsing to handles.** Add a private parallel surface in
  `dsl/edit.rs`. Reuse the ordinary constraint value DSLs,
  while parsing every target-molecule reference as a typed edit handle and building
  `ConstraintEdit` directly: integers become `Id`, `{:new n}` becomes `New`, repeated handles share
  one normalized per-kind slot, and subpattern pattern-side refs remain local to the nested pattern.
  Render typed handles directly from the normalized constraint without exposing its slots.
  Cover single entity leaves; logical trees over existing and created entities; explicit relational
  references of every participating kind; quantified predicates that require no extra atom handle;
  molecule subsets; and anchored subpatterns. Parsing and `IntoAst` enforce only the structural
  integrity of the normalized handle representation; handle liveness and chemical semantics remain
  transaction and validator concerns. **Additive (green).** [dep: S2a, S2b, S2f]
  **Done.**
- **S2h — Expose `EditsDsl` and specify the grammar.** Add the public `EditsDsl` root with
  `FromStr`, `Display`, `FromEdn`, and `ToEdn`, plus `IntoAst<Edits>`/`FromAst<Edits>` under
  `MoleculeDefaults`; re-export it from `dsl.rs`. Parsing must rebuild the eight `Edits` counters in
  one pass through the `Edits` construction surface, while rendering must preserve edit order and
  duplicates. Rendering a `ConstraintEdit` must substitute its typed handles for normalized slots
  without exposing those slots. Add generated round-trip tests spanning every `Edit` variant,
  constraint trees referencing every handle family, and direct conformance cases for the examples in
  this document. Add the normative standalone-edit grammar adjacent to the reaction/delta grammar
  in `umol-ast/spec/umol-dsl-spec.md`, including vector ordering, handles inside constraints,
  defaults, checked updates, removal preconditions, and batching. **Additive (green).**
  [dep: S2b, S2c, S2d, S2g]
  **Done.**

### S3 — Python update values

- **S3a — Bind `UnpairedElectronsUpdate`.** In `umol-py/src/spin.rs`, expose the independent
  optional count and multiplicity updates with complete `from_rust`/`to_rust` conversion. Test empty,
  single-component, two-component, and explicit-undetermined values. **Additive (green).**
  [dep: S1i]
  **Done.**
- **S3b — Bind `StereoConfigurationUpdate`.** In `umol-py/src/stereo.rs`, expose `Unchanged`,
  `Undetermined`, and `Kinded { kind, coset }`, preserving the distinction between an omitted coset
  and an undetermined coset. Test every variant and its Rust round trip. **Additive (green).**
  [dep: S1i]
  **Done.**
- **S3c — Bind `AtomUpdate`.** In `umol-py/src/atom.rs`, expose every optional atom field, the
  nested unpaired-electron update, and `AtomConstraintsAst`. Test the empty value, a mixed field and
  constraint update, and explicit constraint removal. **Additive (green).** [dep: S3a]
  **Done.**
- **S3d — Bind `BondUpdate`.** In `umol-py/src/bond.rs`, expose order, charge, unpaired-electron,
  and constraint updates. Test independent spin-component preservation and constraint removal.
  **Additive (green).** [dep: S3a]
  **Done.**
- **S3e — Bind `DativeBondUpdate`.** In `umol-py/src/dative.rs`, expose order and constraint
  updates and test empty, field, and constraint-removal values. **Additive (green).** [dep: S1i]
  **Done.**
- **S3f — Bind `AromaticSystemUpdate`.** In `umol-py/src/aromatic.rs`, expose electron-count,
  charge, unpaired-electron, and constraint updates. Test mixed updates and independent spin
  components. **Additive (green).** [dep: S3a]
  **Done.**
- **S3g — Bind `MulticenterBondUpdate`.** In `umol-py/src/multicenter.rs`, expose electron-count,
  charge, unpaired-electron, and constraint updates with the same exact-value tests as the Rust
  surface. **Additive (green).** [dep: S3a]
  **Done.**
- **S3h — Bind `NoncovalentBondUpdate`.** In `umol-py/src/noncovalent.rs`, expose kind and
  constraint updates and test explicit undetermined values separately from omitted fields.
  **Additive (green).** [dep: S1i]
  **Done.**
- **S3i — Bind `StereoAtomUpdate`.** In `umol-py/src/stereo.rs`, expose configuration and
  stereo-atom constraint updates, with cases for unchanged, kind-only, absolute, cleared, and
  constraint-removal behavior. **Additive (green).** [dep: S3b]
  **Done.**
- **S3j — Bind `StereoBondUpdate`.** In `umol-py/src/stereo.rs`, expose configuration and
  stereo-bond constraint updates with the same variant inventory and exact conversion checks.
  **Additive (green).** [dep: S3b]
  **Done.**
- **S3k — Rename Python coercion unions from `*Arg` to `*Like`.** Apply the conventional Python
  suffix uniformly to the 21 internal PyO3 input unions: `ValueLike`, `BooleanLike`,
  `ElectronCountsLike`, `ElementLike`, `IsotopeMassLike`, `NoncovalentBondKindLike`,
  `AromaticValenceLike`, `MulticenterValenceLike`, `TetrahedralStereoLike`,
  `CisTransStereoLike`, `TopicityRelationLike`, `StereoConfigurationLike`, `ConstraintsLike`, and
  the eight entity-specific `*ConstraintsLike` types. Update every constructor, setter, helper,
  macro parameter, conversion test, and explanatory comment in `umol-py`; retain `*Input` for DSL
  boundary representations. Verify that the complete Python constructor signatures and coercion
  behavior are unchanged. **Breaking (red→green).** [dep: S3j]
  **Done.**

### S4 — Python `Edit` and `Edits` values

- **S4a — Add the generic `New` handle.** Add `umol-py/src/edit.rs` with the immutable,
  value-equal `New(index)` Python class and internal conversion adapters that accept `int | New` and
  produce the typed Rust handle required by the surrounding argument position. Reject negative and
  overflowing integers. Test equality, repr, immutability, and conversion into all eight typed
  handles without exposing typed Python handle classes. **Additive (green).** [dep: S0b]
  **Done.**
- **S4b — Bind the raw `Edit` enum.** In `umol-py/src/edit.rs`, expose every Rust `Edit` variant
  with exact `from_rust`/`to_rust`, value equality, and repr, reusing the existing field-change,
  entity AST, constraint, and stereo-ligand wrappers. Expose `ConstraintEdit` as the handle-aware
  value carried by molecule-constraint edit variants. Do not add a second static add/remove/update
  construction API on `Edit`. Test every variant in both directions. **Additive (green).**
  [dep: S1i, S4a]
  **Done.**
- **S4c — Bind `Edits` storage and additions.** Expose value equality, construction from entries,
  raw `append(Edit)`, length, indexing, snapshot iteration, and every singular/batched `add_*`
  method with its returned generic Python `New` handle. Test order and duplicate preservation,
  negative indexing, counter recovery from constructor entries, batched results, and independent
  per-kind `New(0)` values. Do not expose `extend`, mutable iteration, insertion, removal,
  reordering, or a list-returning conversion. **Additive (green).** [dep: S4a, S4b]
  **Done.**
- **S4d — Bind `Edits` updates, removals, constraints, and DSL.** Expose all eight mutating
  `update_*` methods, simultaneous topology removal, the six plural overlay-removal methods,
  molecule-constraint add/remove, and `parse(text, *, defaults=None)`/`render(*, defaults=None)`.
  Test zero- and multi-entry updates, exact old-state capture, complete removal preconditions,
  all-family DSL round trips, defaults, and parse-error classification. **Additive (green).**
  [dep: S2h, S3c, S3d, S3e, S3f, S3g, S3h, S3i, S3j, S4c]
  **Done.**
- **S4e — Register and inventory the editing data surface.** Register `New`, `ConstraintEdit`,
  `Edit`, `Edits`, both leaf updates, and all eight entity updates in `umol-py/src/lib.rs`; add
  import-surface assertions for the complete inventory. Typed handles follow the ordinary Python
  convention of using integer entity ids and are covered by the exact package export inventory.
  **Additive (green).** [dep: S4d]
  **Done.**

### S5 — Python editor and transaction lifecycle

- **S5a — Map transaction failures.** In `umol-py/src/error.rs`, add and register
  `TransactionError`, mapping every Rust application and rollback failure to it while leaving DSL
  syntax failures as `ParseError`. Test `HandleOutOfRange`, `HandleRemoved`, representative old-state,
  `DuplicateRemoval`, malformed-edit, and `RollbackStateMismatch` failures and their exact exception
  class. **Additive (green).**
  [dep: S1i]
  **Done.**
- **S5b — Bind `MoleculeEditor` inspection and finalization.** In
  `umol-py/src/transaction.rs`, add and register the Python-owned editor wrapper with non-consuming
  `snapshot()` and consuming `build()`. Store the Rust editor as consumable state so any operation
  after `build` raises `RuntimeError`. Test that repeated snapshots track current state, snapshot
  does not consume the editor, build does, and the built molecule is detached. **Additive (green).**
  [dep: S1j, S4d]
  **Done.**
- **S5c — Bind checked transactions and rollback.** In `umol-py/src/transaction.rs`, add and
  register the consumable `Transaction` wrapper, `MoleculeEditor.transact(edits)`, and
  `transaction.rollback(editor)`. Preserve Rust's detached
  journal semantics, consume the transaction on rollback, and map use of a consumed editor or
  transaction to `RuntimeError`. Test success, failure atomicity, exact rollback, rollback against
  an incompatible editor as an ordinary Python exception rather than a panic, and the rejected
  second rollback. Do not assert the resulting editor state for the deliberately mismatched case.
  Do not expose `transact_unchecked` or `Transaction::append`. **Additive (green).** [dep: S5a, S5b]
  **Done.**
- **S5d — Expose the immutable molecule path.** Add `MoleculeAst.edit()` and
  `MoleculeAst.apply(edits)` in `umol-py/src/molecule.rs`. Test that `apply` returns the expected
  molecule without changing its receiver, surfaces `TransactionError` on failure, and covers the
  methyl-to-ammonia construction that motivates the work. **Additive (green).** [dep: S5c]
  **Done.**

### S6 — Cross-boundary and reader-facing verification

- **S6a — Verify the complete workflow and whitepaper example.** In the Python integration suite,
  exercise one edit sequence through construction, EDN render/parse, immutable apply, editor
  transact, snapshot, rollback, and build; assert exact molecules at every boundary and retain the
  original after both successful and failed operations. Reproduce Section 9's host-specific edit
  listing solely through the public Python API. Run formatting, clippy, the focused and full
  `umol-ast` unit/property suites, the `umol-py` Rust tests, rebuild the extension under the Python
  3.13 virtual environment, and run the complete Python suite. **Additive (green).**
  [dep: S4e, S5d]
  **Done.**

**Critical path:** S0b → S1a → S1b → S1c → S1d → S1e → S1f → S1g → S1h → S1i →
S2a → S2b/S2c/S2d → S2e → S2f → S2g → S2h → S4a → S4b → S4c → S4d → S4e → S5a → S5b →
S5c → S5d → S6a. S1j can proceed after S1i and converges at S5b. S0a is an independent
precedent cleanup. The S3 update wrappers can proceed after S1i and converge at S4d.

No implementation stage is deferrable. S0a makes the container rule honest rather than merely
aspirational; the Rust `Edits` migration and per-kind resolution establish one sound batch boundary;
the DSL makes edits externally serializable; the update, edit, editor, and transaction wrappers form
the promised Python surface; and S6 closes the reader-facing requirement that motivates the work.

## Notes

- Naming is unaffected by doc [176](176-ast-naming-2026-07-31.md): these are new names on the Python
  surface, and if `*Ast` later becomes `*Def` the classes move but the methods do not.

## 2026-08-18 addendum — immutable application failure

The completed S5d surface originally mapped only edit-execution failure from
`Molecule.apply`. Publicly constructed edits can execute successfully and then
fail the checked molecule-integrity publication gate. Rust now reports these
causes through `MoleculeApplyError`; Python retains `TransactionError` for edit
execution and reports failed publication as `InvalidStructureError`. In both
cases the source molecule remains unchanged.
