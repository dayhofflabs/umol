# Validated transactions and operation lifecycle

## Status and scope

This discussion was prompted by kekulization, but the issue is broader than one transformer. A
mutation may be locally well-formed as a sequence of edits while its completed result violates a
structural, physical, or model-dependent invariant. Conversely, a resolver deliberately starts from
an incomplete value, so underdetermination cannot generally be treated as transaction failure.

The design therefore has four coupled parts:

1. transaction acceptance hooks and rollback;
2. an extensible validator abstraction over more than one target type;
3. transformer atomicity and selection of postconditions;
4. resolver atomicity, partial information, and fallback.

This is a design discussion, not an implementation plan. The working direction below is intended to
make the semantic boundaries explicit before APIs are selected.

## Current implementation

### Transactions

`MoleculeEditor::transact(Vec<Edit>)` applies edits immediately. It records a realized `Undo` journal
and reverse-replays that journal if an edit cannot be applied. On success it returns a `Transaction`
that may later be rolled back explicitly. There is no separate begin/commit phase and no validation
of the completed molecular state.

This distinction matters for terminology. A returned `Transaction` is currently an applied rollback
token, not an open transaction awaiting commit. A "pre-commit hook" cannot simply be added to
`Transaction::commit`, because no such operation exists.

Doc 86 specified tier-1 validation at commit and a pre-commit validation hook, but that portion of
the design did not land. Its cost model remains useful: one validation per natural edit batch and an
undo journal proportional to changed entities, not a whole-molecule rollback snapshot. The proposed
API in that document needs to be reconciled with the immediate-application transaction model rather
than copied literally.

### Validators

Validation is closed on two independent axes.

First, validation targets are hard-coded. Individual validators accept `MoleculeAst`; two invariant
validators additionally expose a separate `validate_atom(&AtomAst)` method. `impl
AsRef<MoleculeAst>` permits wrappers around an already-built AST but cannot expose an uncommitted
`MoleculeEditor` state.

Second, composition is hard-coded. `umol_graph::ops::validate::Validator` owns one field for every
known validator, calls them in a fixed tier order, and has manually enumerated
`ValidatorContradiction` and `ValidatorError` unions. Adding a validator requires editing the
composite type, its constructor, its dispatch, and both diagnostic enums. Adding a target requires
new ad hoc methods on the individual and composite validators.

The current three-way result remains valuable:

- `Determined(())`: the target satisfies the check and contains enough information to decide it;
- `Underdetermined(())`: no contradiction is known, but the target is incomplete;
- `Contradictory(c)`: the available information disproves the condition;
- outer `Err(e)`: the validator itself could not execute or was misconfigured.

Transaction acceptance must not flatten these cases prematurely.

### Transformers

`Transformer::transform_into(&mut MoleculeAst)` has no trait-level atomicity contract. Its provided
`transform(&MoleculeAst)` clones the input and delegates to `transform_into`, which protects the
caller's source only for the by-value-result form.

The two current transformers have different failure shapes:

- `Aromatizer` computes perceived systems before writing them. Its fallible planning occurs before
  mutation, but the applied systems are not passed through a general postcondition gate.
- `Kekulizer` can prove the matching requirements during planning, but localized bond orders,
  charge, lone pairs, and the surrounding structure must be checked together afterward. It
  currently creates a candidate AST, validates that candidate, and overwrites the input only on
  success.

The latter is logically transactional but uses a private draft/finalize implementation rather than
the shared `Edit`/`Undo` machinery.

### Resolvers

Resolution currently mutates a `MoleculeAst` directly and sequentially:

1. valence;
2. aromaticity;
3. stereo;
4. localized-bond defaults;
5. multicenter-bond defaults.

This is not atomic today. The counts and atom-typing engines narrow atoms one at a time and can fail
after earlier atoms have been modified. The composite resolver can likewise return a contradiction
after one or more earlier resolver stages have succeeded and mutated the input.

That behavior has not been stated as a deliberate partial-progress contract. It is especially
problematic for fallback: retrying another strategy against a partially narrowed value is not
equivalent to retrying it against the original input.

## Terminology and semantic separation

The following concepts should remain distinct:

- An **edit precondition** establishes that an `Edit` applies to the state it names, for example an
  `old` field value matching the current value. The transaction engine already checks these.
- An **operation input contract** establishes that a transformer or resolver is applicable, for
  example kekulization receiving supported aromatic-system electron demands.
- A **transaction acceptance gate** inspects the completed tentative batch and decides whether it
  may remain applied.
- A **validator** is a read-only computation returning determined, underdetermined, contradictory,
  or execution-error status. It does not decide transaction policy by itself.
- An **acceptance policy** maps validation results to accept/reject. Transformers will usually reject
  contradictions and unexpected underdetermination. Resolvers normally accept underdetermination
  and reject only contradiction or execution failure.
- A **fallback policy** decides what state to restore and which alternative operation to try after
  rejection. It belongs above an individual validator.

Calling all of these "validation hooks" obscures important differences.

## Transaction lifecycle design

### Required post-application gate

The essential primitive is a read-only gate after the complete edit batch has been applied but
before `transact` reports success:

```rust
pub fn transact_validated<E>(
    &mut self,
    edits: Vec<Edit>,
    accept: impl FnOnce(&MoleculeEditor, &TransactionSummary) -> Result<(), E>,
) -> Result<Transaction, ValidatedTransactionError<E>>;
```

The algorithm is:

1. apply edits while recording the existing `Undo` journal;
2. expose the tentative state read-only to `accept`;
3. return the rollback-capable `Transaction` if accepted;
4. reverse-replay the journal if rejected;
5. preserve both causes if rollback after rejection also fails.

The callback must not mutate the editor: otherwise its writes are absent from the journal. It runs
once per batch, never once per edit. `TransactionSummary` is only a placeholder in this sketch.
Local validators may need changed entities and affected neighborhoods, while genuinely global
validators can scan the complete target. The repository already has three relevant representations:
symbolic caller-facing `Edit`, resolved logical `Delta`, and physical rollback `Undo`. The design
must first determine whether a resolved `Delta` sequence or a minimal touched-entity projection is
sufficient. It must not introduce a fourth parallel change vocabulary merely to feed hooks, nor
require validators to reverse-engineer physical `Undo` variants.

The transaction layer should be independent of chemistry and `Solution`. Its generic gate accepts or
rejects; a caller-specific policy performs the mapping from validation outcomes.

### Do we need a pre-application hook?

A symmetric pair of hooks is not automatically desirable. Most current pre-application work is one
of the following:

- edit preconditions, already owned by `transact`;
- operation planning, which should complete before constructing an edit batch;
- an input validator that can be called explicitly before `transact_validated`.

A generic pre-hook earns its place only if it provides functionality that cannot be expressed cleanly
by those mechanisms. Plausible cases are:

- capture a compact baseline used by the post-gate, such as total charge or a resolution progress
  measure;
- validate a relation between before and after states, such as monotone lattice narrowing;
- observe the exact realized transaction boundary uniformly for instrumentation or policy.

These needs may be served without a transaction-owned pre-hook. The caller can calculate context
before invoking `transact_validated`, and old/new edit payloads plus the realized transaction summary
may be enough for relational checks. A closure can capture the context.

Working direction: make the post-application acceptance gate the required primitive. Do not add a
pre-hook merely for lifecycle symmetry. Reconsider a paired before/after protocol only if resolver
monotonicity or another concrete consumer cannot be expressed without copying the baseline state.

### One-shot gate versus explicit transaction scope

There are two materially different APIs:

1. **One-shot validated transaction.** Apply one batch, run one gate, then return success or restore.
   This fits current transformers and can be implemented around the existing journal.
2. **Explicit transaction scope.** `begin`, apply several batches, create savepoints, validate, then
   commit or roll back. This better supports resolver pipelines and alternative search but changes
   borrowing, nesting, and panic/drop behavior substantially.

An explicit scope should not be introduced solely to make the word "commit" fit doc 86. It becomes
justified if resolver fallback needs nested or staged rollback. Until then, one-shot acceptance is
the smaller load-bearing primitive.

## Opening the validator design

### Validation as a target-parameterized capability

The base abstraction should make the target a type parameter rather than encode targets as method
names:

```rust
pub trait Validate<T: ?Sized> {
    type Contradiction;
    type Error;

    fn validate(
        &self,
        target: &T,
    ) -> Result<Solution<(), Self::Contradiction>, Self::Error>;
}
```

This opens both extension axes:

- a new validator implements `Validate<ExistingTarget>` without changing a central composite;
- an existing validator can implement `Validate<NewTarget>` without adding another named method to
  a top-level validator API.

Not every validator must support every target. Unsupported combinations simply have no trait
implementation.

### What is the molecule target?

Parameterizing `Validate` does not by itself let the same algorithm inspect both `MoleculeAst` and
`MoleculeEditor`. Three approaches need comparison:

1. **Build a temporary `MoleculeAst`.** Minimal API work, but may copy materialized editor storage and
   defeats the main efficiency objective.
2. **A single borrowed molecule facade.** Both containers produce a `MoleculeRef<'_>` with the
   read-only atom, bond, overlay, constraint, and topology queries validators need. This avoids
   duplicating validator implementations but requires the facade to abstract the editor's shared and
   mutable storage variants.
3. **Read-capability traits.** Validators are generic over narrowly scoped traits such as atom,
   topology, and overlay access. This is the most open design but can become a large trait hierarchy,
   especially because current view and iterator types are tied to `MoleculeAst`.

The target abstraction should be driven by a method-usage audit of the validators, not by copying the
entire `MoleculeAst` API into one omnibus trait. A borrowed facade is the likely pragmatic starting
point; narrow capabilities remain preferable if the audit reveals clean boundaries.

### Composition

`Validator` should cease being the only way to compose validation. There are three useful levels:

- individual typed validators;
- typed combinators or tuples for a fixed operation-specific pipeline;
- an application convenience bundle representing the standard full chemistry validation sequence.

The current top-level bundle may remain as convenience, but adding an independent validator must not
require editing it. Operation-specific validation should not need to instantiate the full bundle and
then manually select fields, as kekulization currently does.

Typed composition preserves concrete diagnostic types but produces nested sum types unless callers
map diagnostics into an operation error. Dynamic composition permits a registry but requires a
common erased diagnostic and allocation. The likely balance is typed composition at the core,
explicit diagnostic mapping at subsystem boundaries, and no global runtime registry until a real
plugin/configuration consumer appears.

Tier labels—integrity, invariant, and conformance—remain useful metadata and standard groupings. They
should not be the closed dispatch mechanism.

## Transformers and transaction acceptance

A transformer should have a clear atomicity contract: on `Err`, the caller's input is unchanged; on
success, the result satisfies the transformer's declared postconditions. This is stronger and more
useful than allowing each implementation to choose partial mutation semantics.

The natural implementation separates planning from application:

```text
immutable input -> fallible plan -> Edit batch -> tentative apply -> selected postconditions
                -> accept or rollback
```

Input contracts belong in planning or an explicit pre-validation call. Postconditions belong in the
transaction acceptance gate. The transaction engine should not know which transformer is running.

Initial consumer matrix:

| Transformer | Fallible work before edits | Required acceptance gate | Pre-hook need |
|---|---|---|---|
| Aromatizer | perception and system construction | entity structure; aromatic-system consistency; selected invariants/conformance to be decided | none identified |
| Kekulizer | matching, hole selection, localization plan | entity structure, valence invariants, spin invariants; possibly constraints once implemented | none identified |

The matrix must be completed as transformers are added. "Run the global validator" is not an
adequate default: aromatization and kekulization change representation and may intentionally make a
representation-specific conformance check inapplicable during the transition.

Ideally a transform plan lowers entirely to `Vec<Edit>`. Direct field writes bypass the undo journal
and make operation-specific snapshotting necessary. The existing `Edit` vocabulary covers the
kekulizer's atom-field and bond-field changes, atom and bond constraint removal, and
aromatic-system removal; the remaining question is how cleanly the plan emits their required old/new
payloads.

## Resolvers, incomplete inputs, and fallback

Resolvers cannot use the same acceptance policy as transformers. Their input is intentionally
underdetermined, and a successful pass may remain underdetermined. For ordinary resolution:

- `Determined` accepts;
- `Underdetermined` also accepts, possibly with a progress/no-progress distinction;
- `Contradictory` rejects and restores the chosen resolver boundary;
- execution `Err` rejects and restores.

Three additional questions must be settled.

### What is the atomic boundary?

Candidates are:

1. one atom or entity;
2. one resolver stage, such as valence or aromaticity;
3. the complete composite resolver invocation;
4. one fixpoint iteration if resolution becomes iterative.

Per-entity commits preserve useful progress but make strategy fallback and cross-entity checks
difficult. Whole-pipeline atomicity gives the strongest caller guarantee but can discard valid,
expensive narrowing when a later optional stage fails. Stage-level transactions are a plausible
middle ground, provided the public composite contract states whether earlier stages remain visible on
later failure.

The current accidental mixture—per-atom mutation inside stages and persistent earlier stages in the
composite—should not become the contract by default.

### What must a resolver preserve?

Resolution should normally be monotone in the AST lattice: it may replace undetermined information
with narrower information but must not widen or contradict the input. That is a relational
before/after property, not an ordinary validator over the final target alone.

Possible enforcement sources are:

- the `Edit` batch itself, whose old/new values can be checked for narrowing;
- a resolved transaction summary, preferably reusing `Delta` or a projection of it;
- a paired before/after validation protocol if neither is sufficient.

This is the strongest concrete candidate for a pre-transaction context, but it does not yet prove
that a general pre-hook is necessary.

### How does fallback work?

Fallback must retry against a defined baseline. If strategy B is intended as an alternative to
strategy A, B must not observe A's rejected partial narrowing unless that dependency is explicit.

A one-shot stage can be retried by rolling back its transaction and applying the next strategy. A
pipeline that needs to retain earlier common work while trying several later branches needs
savepoints or a stack of transactions rolled back strictly in reverse order. Existing rollback
tokens are not safe as arbitrary historical snapshots after unrelated later transactions compact
dense id spaces.

Therefore resolver fallback may be the consumer that eventually justifies an explicit transaction
scope with nested savepoints. It should be designed from concrete fallback cases rather than assumed
as part of the first validated-transaction change.

## Provisional architecture

The working decomposition is:

```text
operation planner
    -> Vec<Edit>
    -> MoleculeEditor transaction + Undo journal
    -> read-only tentative molecule target
    -> operation-selected validator composition
    -> acceptance policy
    -> accept Transaction | rollback with preserved diagnostic
```

For resolvers, add a separately stated atomic boundary and fallback controller around that core.

This keeps responsibilities narrow:

- transactions provide tentative application and exact restoration;
- molecule read abstractions expose tentative state without materialization;
- validators report facts and uncertainty;
- transformers and resolvers select validators and acceptance policy;
- fallback orchestration chooses retry baselines and alternatives.

## Questions to resolve before implementation planning

1. Can a borrowed `MoleculeRef` cover all current validator queries without allocation, including
   editor storage that has been materialized into mutable vectors?
2. Does an acceptance gate need the resolved `Delta` sequence, only a touched-entity projection, or
   no change summary at all? What stable entity-id semantics can it promise across removals and
   compaction?
3. Is post-application acceptance sufficient for the first implementation, with captured caller
   context for relational checks, or does a concrete resolver require a paired hook protocol?
4. Should `Transformer::transform_into` guarantee input preservation on every error, or should the
   mutating method be replaced by a plan/apply API that makes this automatic?
5. Which exact validators are postconditions of aromatization and kekulization? In particular, when
   should representation-specific conformance run?
6. What is the public atomicity contract of `Resolver::resolve`: per stage or whole composite call?
7. Which real resolver alternatives require fallback, and at what shared baseline? This determines
   whether one-shot transactions suffice or savepoints are required.
8. Should underdetermined-but-no-progress resolution be a successful `Underdetermined`, a separate
   status, or an orchestration-level fixpoint condition?
9. How are typed validator diagnostics composed without recreating a manually closed global union?
10. Which current direct mutations cannot yet be represented faithfully as `Edit`s?

## Suggested design-spike sequence

Before producing a staged implementation plan:

1. audit read methods used by every validator and prototype an allocation-free view over both
   `MoleculeAst` and `MoleculeEditor`;
2. prototype one-shot `transact_validated` using the existing journal and a generic closure error;
3. lower kekulization to an edit batch and use it as the post-acceptance reference consumer;
4. construct a resolver failure case that currently leaves partial mutation, then compare whole-call,
   stage, and per-entity rollback semantics;
5. use that resolver case to decide whether captured context is enough or explicit savepoints are
   required;
6. only then settle the validator composition API and operation trait changes.
