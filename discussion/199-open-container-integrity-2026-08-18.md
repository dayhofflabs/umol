# 199 — Open-container integrity and operation boundaries

Status: Proposed
Date: 2026-08-18
Relates: [148](148-validated-transactions-operations-2026-07-13.md),
[179](179-python-editing-and-transactions-2026-08-02.md),
[198](198-resolver-performance-2026-08-17.md),
[201](201-molecular-data-first-steps-2026-08-19.md),
[data-type contracts](../docs/development/data-types.md),
[nomenclature](../docs/development/nomenclature.md)

## Motivation

`Molecule`, `Reaction`, and `ReactionSpan` are intended to be relatively open
scientific data containers. They support direct inspection and meaningful
editing, retain non-normal and non-canonical representations, and defer
normalization, canonicalization, resolution, and semantic validation to named
operations. Ordinary use should not require a heavily guarded mutation
protocol or a fallible setter for every field.

The containers must nevertheless have well-defined behavior. Some stored
states cannot safely support every operation: references may not resolve,
participant-indexed data may not match its frame, two relations may conflict
under the relation family's intended semantics, or a reaction may not
materialize a two-sided span. The design question is where each property is
required and checked, not whether integrity is valuable in the abstract.

The current development guides take a stronger position: representation
integrity is established by construction and preserved by every public
mutation that produces an aggregate graph-IR value. The implemented APIs also
provide direct live mutation, a deliberately permissive `Reaction`, and
operation-specific checks. This document reopens that policy. The guides remain
the current record until the discussion settles, then must be revised to match
the selected model.

Resolver profiling prompted the review because repeated editor finalization and
whole-molecule integrity checks are visible costs. Resolver performance is one
consumer of the result; it must not determine the general container semantics
by itself.

## Scope

This discussion covers:

- the minimum properties required for safe storage and ordinary access;
- the broader coherence properties reported by `check_integrity`;
- which properties construction establishes and direct mutation preserves;
- the preconditions and failure behavior of operations that need stronger
  coherence;
- the meanings of `build`, `try_build`, `snapshot`, `apply`, and `transact`;
- the distinction between edit-application failure and a post-application
  integrity failure;
- the corresponding Rust and Python behavior;
- the implications for molecule and reaction operations, including resolver
  phase publication.

This discussion does not select chemistry validation policy, require eager
normalization or canonicalization, remove direct mutation as a design goal, or
introduce `Validated<T>`, typestate, witness values, or another public wrapper
layer. It does not contain an implementation plan.

## Starting position

The following principles are established requirements for this discussion:

1. Aggregate graph-IR types remain usable as open, editable containers.
2. Direct mutation is an important public capability in both Rust and Python.
3. Normalization, canonicalization, resolution, and semantic validation remain
   explicit and lazy.
4. Fallibility belongs at an operation boundary where the caller can act on it;
   it is not propagated through unrelated methods solely because an internal
   implementation can fail.
5. A checked convenience route may coexist with open editing, but it must not
   make the ordinary path ceremonial.
6. Performance may influence where equivalent checks occur, but it does not
   define which states the containers represent.

The unresolved issue is how much of the current `check_integrity` predicate is
an invariant of every live container and how much is a condition checked only
by operations that require it.

## Current surfaces

### Molecule construction and mutation

`Molecule::try_from_entries` constructs the graph and relation stores and runs
`Molecule::check_integrity`. `Molecule::from_entries` uses the same path and
asserts success. The current integrity check covers storage/table agreement,
references, participant uniqueness, cross-entity relation uniqueness,
participant/electron-count alignment, stereo frame and value domains, and
constraint references.

After construction, `Molecule` exposes mutable entity views, whole-family
attribute modification, and mutable molecule constraints. These are useful
editing surfaces. Some changes preserve every current integrity condition;
others can create a state that a later `check_integrity` rejects, such as a
participant-indexed electron vector of another length, a stereo value outside
the stored ligand frame, or a constraint naming an unavailable entity.

The Python bindings intentionally expose corresponding live views and setters.
Any selected policy must therefore be coherent across both languages; changing
only Rust construction or only Python setters would leave two different data
models.

### Editor finalization

`MoleculeEditor` is explicitly transient and permits multi-step construction.
Its low-level additions and mutable views need not make every intermediate
state suitable for all molecule operations.

Currently:

- `try_build` consumes the editor, materializes a `Molecule`, runs the complete
  integrity check, and returns `MoleculeIntegrityError` on failure;
- `build` delegates to `try_build` and asserts success;
- `snapshot` clones the editor and uses the checked `try_build` path.

This gives `build` the meaning "asserted checked finalization." An open-container
model could instead give `build` the meaning "ordinary finalization" and reserve
`try_build` for an explicitly checked result. Neither meaning is selected yet.

### Edit application and transactions

`MoleculeEditor::apply` consumes an editor and an `Edits` batch. It checks
handles, old-state preconditions, edit shapes, and the evolving batch state. On
failure the partially modified editor is inaccessible and is dropped. On
success it returns another transient editor and does not itself promise that
the complete `Molecule::check_integrity` predicate holds.

`MoleculeEditor::transact` performs the same edit application while retaining a
realized undo journal. A failed batch restores the borrowed editor. A successful
batch returns a `Transaction` that may later be rolled back. Transaction
atomicity concerns application and restoration; it does not by itself settle
the integrity status of the resulting editor.

`Molecule::apply` is the immutable-style convenience operation. Its current
sequence is:

```text
source molecule
    -> edit
    -> consuming editor apply
    -> checked editor try_build
    -> result molecule
```

Its `MoleculeApplyError` distinguishes edit-execution `TransactionError` from
post-application `MoleculeIntegrityError`. The operation leaves its source
unchanged under either failure and does not publish the transient editor.
This checked boundary is selected for the 0.6 API without settling the broader
meaning of ordinary editor finalization or the complete open-container model.

### Reaction and reaction span

`Reaction::new` deliberately accepts an lhs and deltas without checking whether
the deltas can be interpreted in every later context. Operations such as span
materialization and application establish the properties they require.

`ReactionSpan::try_from_entries` currently takes a stronger construction
position: both side projections must form integral molecules. Its projections
and conversion back to `Reaction` are consequently infallible. Whether that
stronger boundary remains appropriate should be assessed from the meaning of a
span and the usefulness of representing incomplete spans, not forced to match
either `Molecule` or `Reaction` mechanically.

### Operation checks

Canonicalization, remapping, reaction operations, boundary conversion, and
resolution currently invoke integrity checks at different boundaries. Some
operations rebuild complete aggregates through checked constructors; some
prove transport properties and use internal construction directly; some accept
permissive inputs and check only when a stricter representation is requested.

The inventory is useful evidence, but existing placement does not settle the
general policy. Each operation must state the subset of properties required to
produce its promised result and what happens when those properties do not hold.

## Property spectrum to classify

The current single integrity category contains properties with different
operational consequences. The discussion must classify them before deciding
where checks belong.

| Property group | Examples | Question to settle |
| --- | --- | --- |
| Storage and indexing shape | graph node/atom and edge/bond table agreement; dense table bounds | Must this always hold so ordinary access remains defined? |
| Stored reference resolution | overlay participants; stereo sites and ligands; constraint entity ids | May an open container temporarily carry unresolved references, and which accessors remain available? |
| Entity-local interpretability | participant/electron-count length; stereo frame arity and coset/permutation domain | Is this required to inspect the entity, or only for operations using the affected field? |
| Relation-family coherence | duplicate participants, parallel pairs, overlapping aromatic systems, duplicate stereo sites | Is this part of what the container stores, or a condition imposed by consumers using relation-set semantics? |
| Contextual compatibility | a correspondence against supplied molecules; a reaction against a host | Check at the first combining operation. |
| Semantic validity | physical invariants and chemistry-model conformance | Continue to use explicit validators. |
| Normal form | normalized values and constraints; canonical entity frame | Continue to use explicit transformations. |

The names of these groups are descriptive placeholders. The discussion may
retain one integrity predicate, divide it into smaller predicates, or keep one
diagnostic operation while allowing only part of it to be an always-held
container invariant.

## Contract sheets to settle

### Molecule

```text
Type and role: editable aggregate molecular graph IR
Open carrier or operation-issued value: open carrier
Intrinsic representation invariants: unresolved; at minimum, safe storage and indexing
Contextual properties and supplied context: operation-specific models, algorithms, correspondences, and hosts
Semantic predicates and validators: invariants and conformance in umol-graph
Public constructors: new, builder, from_entries, try_from_entries
Conversions and preserved information: DSL and external-format boundaries; no implicit normalization or repair
Explicit transformations: normalize carried forms, canonicalize, resolve, and chemistry-aware transforms
First public consumer requiring each contextual property: to be inventoried per operation
Failure, absence, and panic behavior: open; especially build, apply, remap, and canonicalize
Algebraic, preservation, or roundtrip properties: existing equality, remapping, editing, and rollback laws remain relevant
Rust/Python boundary: the same openness and failure categories in both languages
```

### Reaction

```text
Type and role: lhs molecule plus resolved delta sequence
Open carrier or operation-issued value: open carrier
Intrinsic representation invariants: unresolved; construction is currently permissive
Contextual properties and supplied context: span materializability, host and match compatibility, DPO conditions
Semantic predicates and validators: operation-specific reaction and chemistry checks
Public constructors: new and from_sides
Conversions and preserved information: to/from ReactionSpan where materializable
Explicit transformations: reverse, compose, canonicalize, and application-derived operations
First public consumer requiring each contextual property: span conversion, composition, canonicalization, or application
Failure, absence, and panic behavior: preserve operation-specific distinctions
Algebraic, preservation, or roundtrip properties: side/span and reverse/composition laws
Rust/Python boundary: permissive construction and contextual failures should agree
```

### ReactionSpan

```text
Type and role: union-frame representation of two reaction sides
Open carrier or operation-issued value: currently a checked aggregate; desired openness is unresolved
Intrinsic representation invariants: union storage plus the disputed requirement that both projections form molecules
Contextual properties and supplied context: dense remappings and source correspondences
Semantic predicates and validators: no chemistry validation during construction
Public constructors: from_entries, try_from_entries, and superimpose
Conversions and preserved information: lhs, rhs, to_reaction, and remapping
Explicit transformations: normalization and canonicalization
First public consumer requiring each contextual property: projection, conversion, remapping, or canonicalization
Failure, absence, and panic behavior: depends on whether projections remain guaranteed
Algebraic, preservation, or roundtrip properties: side projection and span/reaction normal-form laws
Rust/Python boundary: checked construction and projection behavior must remain aligned
```

## `apply` and finalization alternatives

The apply question has two independent axes:

1. What does successful edit application establish?
2. What does converting an editor into a `Molecule` establish?

The following alternatives are open.

### A. Open-result application

`MoleculeEditor::apply` continues to mean successful execution of the edit
batch. Ordinary `build` finalizes the current state without a complete
integrity check, while `try_build` remains the explicit checked finalization.
`Molecule::apply` uses ordinary finalization and reports only transaction
failures.

This aligns closely with open-container editing and keeps the ordinary path
small. Operations receiving the result must check the properties they require.
The design must specify which minimum storage properties edit application and
ordinary finalization still guarantee.

### B. Checked molecule application

`MoleculeEditor::apply` still returns a transient editor, but
`Molecule::apply` performs checked finalization. Its result type exposes both
edit-application and integrity failures through an operation-specific error.
`build` may remain asserted checked finalization.

This gives the immutable-style convenience a stronger postcondition, at the
cost of a broader fallible surface and a complete check even when the caller
does not require it. Python must expose the additional failure category.

**Selected for the 0.6 `Molecule::apply` boundary (2026-08-18).** Rust returns
`MoleculeApplyError::Transaction` or `MoleculeApplyError::Integrity`. Python
maps those causes onto the existing `TransactionError` and
`InvalidStructureError` exception classes. This decision does not settle
`build`, `try_build`, direct mutation, or which integrity properties an open
container must always preserve.

### C. Integrity-preserving edit application

The edit engine checks enough local and cross-entity conditions that applying a
batch to a suitable source establishes the selected integrity postcondition.
`Molecule::apply` can then retain its current compact result type, and
finalization may rely on the edit operation's proof.

This may support incremental checks over touched entities, but it increases the
semantic responsibility of edit application and may reject edit sequences that
are useful while constructing an open intermediate state. The postcondition
would need to differ between `MoleculeEditor::apply` and `Molecule::apply`, or
the source editor would need a documented precondition.

### D. Asserted post-application finalization

Keep the current shape: `Molecule::apply` reports transaction failures and
asserts that a successfully applied batch also passes final integrity checking.
This treats a post-application integrity failure as a broken producer
precondition.

This is compact and may be acceptable for trusted plans, but public caller-built
`Edits` can reach the same method. The documentation would need to state the
precondition precisely and decide whether a panic is an appropriate response
for ordinary caller input.

### E. Parallel checked and open operations

Expose separate checked and open variants of application or finalization.

This is mechanically flexible but adds names and makes callers choose a
lifecycle policy at routine call sites. It is not the default direction unless
distinct recurring consumers justify both routes.

These alternatives may be combined selectively. For example, `build` could be
ordinary finalization while `Molecule::apply` remains checked, or transaction
application could preserve a minimal safe-storage subset while leaving broader
coherence to `check_integrity`.

## Operation-boundary questions

The complete audit should record, for each operation:

| Operation | Current boundary to inspect | Decision required |
| --- | --- | --- |
| Direct entity and constraint mutation | Live mutable form or container | Which properties must the mutation preserve, and may later consumers reject the state? |
| `from_entries` / `try_from_entries` | Complete integrity check | Do both continue to accept exactly the same states, and what initial guarantee do they establish? |
| `build` / `try_build` / `snapshot` | Complete integrity check | Ordinary versus checked versus asserted finalization semantics |
| `MoleculeEditor::apply` | Edit handles, shapes, and old-state checks | Whether success also establishes any aggregate coherence |
| `MoleculeEditor::transact` | Atomic application with undo journal | Keep rollback semantics independent from aggregate integrity, or add an acceptance gate? |
| `Molecule::apply` | Transaction or integrity failure through checked `try_build` | Retain the selected bounded behavior while the wider container model is reviewed |
| `combine` / `combine_from` | Rebuild and checked publication | Whether disjoint concatenation can rely on its input properties |
| `remap` / `try_remap` | Source integrity and correspondence checks | Which source properties the transport actually requires |
| Canonicalization | Full integrity precondition | Whether to require all current integrity conditions or an operation-specific subset |
| Resolution | Repeated checked intermediate publication | Entry, phase, and final-output requirements |
| Reaction span materialization | Reaction and projected-side checks | Whether a span continues to guarantee two usable projections |
| Boundary conversion | Checked target construction | Which failures belong to conversion and which remain representable target states |

This matrix must include current Rust and Python entry points, public rustdoc,
tests, property suites, examples, benchmarks, and the specification where it
promises behavior.

## Performance considerations

The current ChEBI profile shows repeated checked publication as a meaningful
part of resolver time. A single complete integrity pass over the corpus is much
smaller than the repeated resolver cost, so lifecycle repetition and
copy-on-write materialization need to be separated from the cost of the
predicate itself.

Possible optimizations depend on the selected semantics:

- operations may check only the properties they require;
- a trusted transformation may reuse properties established by its source and
  its own edits;
- an incremental check may examine touched entities and affected incidence;
- an ordinary open finalization may skip broader coherence entirely;
- the complete explicit `check_integrity` operation may still be optimized as a
  diagnostic and boundary tool.

No resolver-only unchecked path should be selected before the general meanings
of open state, finalization, and operation preconditions are settled.

## Documentation consequences

Once the design is settled, update together:

- `docs/development/data-types.md`, especially the categorical statements that
  every public mutation preserves integrity and every aggregate publication
  runs the complete check;
- the `Integrity`, `Integrity check`, `Application`, and `Transaction` entries
  in `docs/development/nomenclature.md`;
- the `data-type-contracts` skill, whose current classification assumes the
  stronger construction guarantee;
- Rust and Python public documentation for construction, editing,
  finalization, application, canonicalization, remapping, and reaction spans;
- specifications and tests that currently imply a stronger or weaker boundary.

Doc 148 remains the broader discussion of transaction acceptance, validator
composition, transformer atomicity, and resolver lifecycle. This document must
settle the underlying aggregate-container semantics before those transaction
and lifecycle proposals are revised. Doc 198 remains the performance evidence
and should consume the result without independently defining it.

## Open questions

1. Which properties are indispensable for safe ordinary access to a live
   `Molecule`?
2. Which current `MoleculeIntegrityError` conditions may be ordinary temporary
   states of an open container?
3. Does `check_integrity` remain one broad diagnostic predicate, or should the
   always-required and operation-specific subsets receive distinct names?
4. What initial guarantee should `from_entries` and `try_from_entries`
   establish, and must they accept the same input domain?
5. What should `build`, `try_build`, and `snapshot` each mean?
6. **Answered for the 0.6 boundary:** `Molecule::apply` uses checked
   finalization and reports transaction and integrity failures separately.
   The wider container decision may revisit this in a later breaking release.
7. Does transaction success say anything about aggregate coherence beyond
   successful edit execution?
8. Which operations require the complete current integrity predicate, and
   which need only a subset?
9. Should `ReactionSpan` retain a stronger construction boundary than
   `Molecule` and `Reaction` because its name promises two projectable sides?
10. How should live Python setters and collection views reflect the selected
    policy without making routine editing cumbersome?
11. After semantics are settled, which repeated checks can be removed, shared,
    or made incremental without changing observable behavior?
