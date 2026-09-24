# 213 — Molecule and reaction mutation

Status: In Progress
Date: 2026-08-27
Relates: [117](117-entity-model-extensibility-2026-06-20.md),
[166](166-molecule-ops-2026-07-27.md),
[211](211-relation-frames-and-api-2026-08-26.md),
[214](214-aggregate-frame-semantics-2026-08-28.md),
[228](228-python-api-parity-2026-09-21.md),
[229](229-aggregate-integrity-review-2026-09-22.md),
[data-type guide](../docs/development/data-types.md),
[nomenclature guide](../docs/development/nomenclature.md)

## Design status — 2026-09-24

This document owns the molecule/reaction mutation redesign. S0 and S1a are
implemented. Graph-core mutation and restoration are complete in
[166](166-molecule-ops-2026-07-27.md); editor integration remains here. After
that integration, return to 166 for the operation changes and hydrogen folding.
Doc 228 is unchanged by this review and its withdrawn ownership migration is not
an implementation dependency.

| Area | Status | Concrete position |
| --- | --- | --- |
| Storage delegation, participant methods, Edit/Delta/Undo variants, local getters | Settled design; S1a complete | Use the existing typed entity sets and graph-core mutation/restoration; contracts below. |
| Editing and recovery | Settled design | Owning, destructive editor; separate borrowed, scoped transaction. Editor and Transaction probe check integrity and return an immutable Molecule borrow; no probe callback. |
| resolve/project/transform consumers | Settled design; integration work remains | resolve/project consume destructively; resolve_into/project_into mutate borrowed inputs with recovery. Consuming resolution uses Solution<Molecule, C, ()>; reporting is explicit. Ingest uses report-free resolution. Transformer signatures follow the same ownership naming. |
| Molecule attribute methods | Settled: retain current placement | Keep the eight integrity-preserving mutable methods and the nine checked callbacks on Molecule, with their existing guarantees. The editor retains the complete mutation vocabulary. |
| Transaction correspondence | Settled design | tracked_commit returns the whole transaction's correspondence. Omit Transaction::tracked_apply unless a concrete need for intermediate tracking arises. |
| Python bindings | Prepared-batch transaction interface settled; migration remains | Molecule.transact and tracked_transact submit prepared Edits; Rust applies and commits within one borrowed transaction. No interactive Python Transaction or scoped TLS dependency. Molecule and Edits input-transfer changes remain; the editor already supports consumption. |
| Edits and multiple batches | Settled | Edits accumulates one sequence. Multiple batches execute through separate Transaction::apply calls under one commit/rollback boundary. No independent-batch composition API on Edits. |
| Mutation errors | Settled design | Retain application/integrity categories and chemistry outcomes; add Aborted and remove obsolete rollback failures. ResolveError::Apply and ProjectError::Apply carry MoleculeApplyError. |

The staged implementation plan below sequences these contracts and integration
obligations. S0 records the baseline and additive Solution type; S1a adds aromatic
and multicenter set mutation. Editor and molecule mutation API changes have not
started. The lift_constraints defect and its undetermined-stereo policy are a
separate focused correction, recorded under
[other operations](#other-moleculereaction-operations).

## Editor and transaction API

```rust
pub struct MoleculeEditor { molecule: Molecule }
pub struct Transaction<'scope> { /* private borrowed state and journal */ }

impl Molecule {
    pub fn edit(self) -> MoleculeEditor;
    pub fn apply(self, edits: Edits) -> Result<Self, MoleculeApplyError>;
    pub fn tracked_apply(
        self, edits: Edits,
    ) -> Result<(Self, MoleculeCorrespondence), MoleculeApplyError>;
    pub fn transact(
        &mut self,
        batches: impl IntoIterator<Item = Edits>,
    ) -> Result<(), MoleculeApplyError>;
    pub fn tracked_transact(
        &mut self,
        batches: impl IntoIterator<Item = Edits>,
    ) -> Result<MoleculeCorrespondence, MoleculeApplyError>;
}

impl MoleculeEditor {
    // All direct entity/constraint mutation remains available here.
    pub fn apply(self, edits: Edits) -> Result<Self, MoleculeApplyError>;
    pub fn tracked_apply(
        self, edits: Edits,
    ) -> Result<(Self, MoleculeCorrespondence), MoleculeApplyError>;
    pub fn probe(&self) -> Result<&Molecule, MoleculeIntegrityError>;
    pub fn finish(self) -> Result<Molecule, MoleculeIntegrityError>;
}

impl Transaction<'_> {
    pub fn run<T, E>(
        molecule: &mut Molecule,
        f: impl for<'scope> FnOnce(Transaction<'scope>) -> Result<T, E>,
    ) -> Result<T, E>
    where E: From<MoleculeApplyError>;

    pub fn apply(&mut self, edits: Edits) -> Result<(), MoleculeApplyError>;
    pub fn probe(&self) -> Result<&Molecule, MoleculeIntegrityError>;
    pub fn commit(self) -> Result<(), MoleculeApplyError>;
    pub fn tracked_commit(self)
        -> Result<MoleculeCorrespondence, MoleculeApplyError>;
    pub fn rollback(self);
}
```

The block expresses the accepted editor and transaction interface. Currently
Molecule::edit and Molecule::apply borrow &self and create an independent result;
editor transact returns a detached journal; snapshot/try_build/build publish the
working state. Those lifecycle methods are replaced by the surface above.
MoleculeBuilder retains its asserted build for fresh construction.

Transaction::run is the scoped execution entry point, not a constructor returning
a transaction. It passes the borrowed handle to the callback and retains the
recovery guard until that callback returns. Molecule::transact and
tracked_transact are prepared-batch conveniences over run: apply each batch,
then commit or tracked_commit. A single batch is passed as [edits]. The callback
entry point is Rust-only; both prepared-batch methods are exposed in Python under
the same names. No Molecule::transaction or with_transaction method is added.

| Method/path | Success | Error, cancellation, or drop | Copies and integrity |
| --- | --- | --- | --- |
| edit and direct methods | One owned editor; direct changes and batches may mix. | No recovery promise. | No implicit clone, journal, or correspondence. Intermediate state may be incomplete. |
| Editor apply | Returns the same continuing owner. | Drops the editor, including earlier changes. | Checks edit preconditions; no full integrity check per batch. |
| probe | Returns an integrity-checked immutable borrow, tied to &self. | Integrity error; no borrow is returned; editing can continue. | No implicit copy. A deliberate clone of the checked molecule is allowed. |
| finish | Returns a constructor-equivalent Molecule. | Drops the editor. | One integrity gate. |
| Molecule apply | edit, editor apply, finish. | Drops the input. | Same checks and no recovery storage. |
| Molecule transact / tracked_transact | Applies all batches and commits; tracked form returns the overall correspondence. | Restores the receiver on failure. | Journal-based recovery, no molecule clone; one integrity gate at commit. |
| Transaction apply | Applies immediately; later batches can use changed state. | Restores the whole transaction and prevents further mutation or commit. | Same preconditions, plus private undo recording. |
| commit | Checks integrity and requests acceptance. | Failed integrity restores and aborts. | Acceptance requires both successful commit and callback Ok. |
| rollback or callback without commit | Restores transaction-entry state. | Callback Err is returned after restoration; unwinding also restores. | No recovery clone or empty replacement molecule. |

The probe borrow prevents mutation or consumption while the reference
remains in use. Once its last use ends, editing continues normally. A transaction
reference cannot escape its scope. Integrity is checked before returning the
reference, so an explicit clone has the same constructor guarantee as finish.
The integrity check establishes the observation barrier; a probe callback adds
no protection against observing an invalid molecule. Transaction::run retains
its scoped callback to keep responsibility for restoration inside the library.

### One Edits or separate apply calls

- Appending edits to one Edits creates one continuous sequence. Id refers to an
  entity in the molecule before that sequence starts; New numbering continues
  throughout the sequence, separately for each entity kind.
- Two separate apply calls execute two sequences. The second call's Id refers to
  the molecule after the first call, and its New numbering starts over. These
  rules also apply to the editor's tracked_apply.

For example, two separately applied Edits can each create an atom called New(0):
those handles refer to different atoms. Transaction processing applies the batches
separately and gives them one commit/rollback boundary. Edits::push appends an
entry written for the receiving sequence; it does not rebase another batch.

### Transaction application and probe errors

If apply fails, the entire transaction is rolled back, including
earlier successful calls. Further application and commit fail. Ignoring the error
does not allow the transaction to succeed: the transaction call returns Aborted
unless its callback returns its own error.

If probe fails its integrity check, it returns no reference and does not undo any
changes. The transaction stays active, so later edits can repair the state. After
an application failure has already rolled back the transaction, probe can only
observe the restored original.

### Restoration and ownership

Restoration is modulo normalized_eq; exact non-normal form encodings need not
survive. Participant sequences retain their separately settled exact comparison.
Manipulated undo data must not panic and has no correctness guarantee. Undo stays
public; the journal and scope guard are private. Forgetting the supplied
Transaction cannot forget the library's restoration guard. An apply panic caught
inside the callback also requires restoration before another call; the
implementation obligations below cover this case.

There is one owning editor, without a mode or borrowed alternative. Transaction
borrows Molecule storage directly and shares private mutation kernels with the
editor. It exposes batches, not unrestricted mutable views. Transaction::run
borrows the receiver, so errors do not need to carry the original molecule back
to the caller. An editor does not provide transaction entry.

### Transaction correspondences

Use ordinary apply for every batch. tracked_commit has commit's integrity,
acceptance, and error behavior, and returns the entity-id mapping from transaction
entry to the final state. It constructs that correspondence from the created ids
and removal compactions already recorded in the private undo journal. The
transaction does not maintain an additional correspondence throughout execution;
ordinary commit constructs none.

The journal remains available for rollback until the callback succeeds. A caller
may return the mapping through that callback; a callback error returns no
successful result. Transaction::tracked_apply is omitted: no current consumer
needs an intermediate batch correspondence. Add it only if that need arises.

Editor tracked_apply returns only its batch's correspondence; preceding/following
direct changes are outside that mapping. No path adds participant-position
tracking: arbitrary replacement is not a frame permutation. Correspondence
materialization can require receiver-sized storage; ordinary direct mutation does
not pay for it.

## Resolver and transformer interfaces

The resolve/project ownership names are settled:

| Method inputs | Ownership and failure contract | Execution |
| --- | --- | --- |
| resolve(Molecule), project(Molecule, ProjectFlags) | Consume the input; return the changed molecule on accepted completion. Rejection drops the input without recovery. | Owning editor; no implicit copy or undo journal. |
| resolve_into(&mut Molecule), project_into(&mut Molecule, ProjectFlags) | Mutate the borrowed input; restore it on rejection. | Scoped transaction; no recovery copy. |

Both routes establish constructor-equivalent integrity before returning a molecule
or accepting a change. Extend Solution with a separately typed Underdetermined
payload, defaulting to the Determined payload type:

```rust
pub enum Solution<T, C, U = T> {
    Determined(T),
    Underdetermined(U),
    Contradictory(C),
}
```

Existing two-parameter uses retain their meaning. Methods that operate on a shared
payload (data, into_data, map) retain their existing behavior for Solution<T, C, T>;
predicates and methods that do not require equal payload types generalize to U.
No new mapping-method family is required by this decision.

Consuming resolution returns a molecule only in Determined; Underdetermined
discards it. Reporting remains explicit through the with_report suffix; ordinary
resolution does not construct a ResolveReport. The corresponding signatures are:

```rust
impl Resolver<'_> {
    pub fn resolve(&self, molecule: Molecule)
        -> Result<Solution<Molecule, ResolveContradiction, ()>, ResolveError>;
    pub fn resolve_with_report(&self, molecule: Molecule)
        -> Result<Solution<(Molecule, ResolveReport), ResolveContradiction, ResolveReport>, ResolveError>;
    pub fn project(&self, molecule: Molecule, flags: ProjectFlags)
        -> Result<Solution<Molecule, ProjectContradiction, ()>, ProjectError>;

    pub fn resolve_into(&self, molecule: &mut Molecule)
        -> Result<Solution<(), ResolveContradiction>, ResolveError>;
    pub fn resolve_into_with_report(&self, molecule: &mut Molecule)
        -> Result<Solution<ResolveReport, ResolveContradiction>, ResolveError>;
    pub fn project_into(&self, molecule: &mut Molecule, flags: ProjectFlags)
        -> Result<Solution<(), ProjectContradiction>, ProjectError>;
}
```

resolve_with_report returns the report alongside the molecule when determined,
or the report alone when underdetermined. Contradiction and execution-error
categories remain unchanged. No copying convenience is introduced.

The Transformer interface follows the same ownership naming. Rename generate_all to
transform_iter; this method borrows the source and yields independent outputs:

```rust
pub trait Transformer {
    type Error;
    fn transform(&self, molecule: Molecule) -> Result<Molecule, Self::Error>;
    fn transform_into(&self, molecule: &mut Molecule) -> Result<(), Self::Error>;
    fn transform_iter<'a>(&'a self, molecule: &'a Molecule)
        -> impl Iterator<Item = Molecule> + 'a;
}
```

Current transform borrows and clones its input before calling transform_into.
The consuming interface instead uses the owning editor. transform_iter retains the
borrowed input and creates candidates on demand, preserving the current empty
iterator on failure. All three current generate_all implementations eagerly
transform before returning a boxed zero-or-one-element iterator. Defer that work
until iteration; removing Box alone does not make it lazy.

Use impl Iterator: each implementation can return its own concrete iterator
without a required heap allocation or virtual next call. This makes Transformer
unavailable as a dyn trait with this method; the workspace has no such callers
or bindings. No current consumer justifies retaining the boxed return.

The caller of resolve_into/project_into keeps its molecule on every outcome.
Determined installs the change; Underdetermined, Contradictory, and Err restore
the entry state modulo normalized_eq. Default resolve_into returns Determined(()) or
Underdetermined(()); contradictions and execution errors retain their causes.
resolve_into_with_report additionally returns unresolved candidates and tie-break
information in ResolveReport. Neither returns a partially resolved molecule.
Discarding the optional report gives the same outcome, errors, and final molecule
as ordinary resolve_into. Borrowed projection retains its unit success payload
and does not require concreteness.

Ordinary resolution must skip report construction and bookkeeping used only by the
report. Keep the candidate sets needed for resolution itself. Replace the current
to_report().unresolved.is_empty() decision with a direct test of those candidates;
do not clone unresolved alternatives or collect/sort tie-break records merely to
discard them. All resolution entry points share the chemistry routines and retain
the same phase order and acceptance conditions.

Current [ingest](../umol-graph/src/ingest.rs) and
[MOL parse](../umol-graph/src/parse.rs) place ResolveReport in their
ResolveUnderdetermined errors. In the redesign, these boundaries use ordinary
report-free consuming resolution. They return a molecule only when determined;
underdetermination drops the candidate and returns a payload-free
ResolveUnderdetermined. Contradiction and execution failures retain their causes.
Python ingestion likewise raises UnderdeterminedError without a report attribute.
Explicit resolver reporting remains available to callers that request it.
The current ingestion tests asserting report contents and the Python exception
adapter must change with this contract.

resolve_into, project_into, and transform_into update the caller's molecule on
success and restore it on rejection. No recovery copy is required. Ingest, parse,
and export already own their candidates and can call consuming resolve/project
without a recovery journal.

### Ownership, cost, and result at each call site

| Caller or entry point | Selected execution | Success and failure | Additional storage/work |
| --- | --- | --- | --- |
| Public composite resolve_into(&mut M) | One transaction across all phases. | Determined commits; all other outcomes restore. Unit payload by default; resolve_into_with_report returns the optional report. Contradiction/error causes survive. | No recovery clone or default report construction; journal retained across seven batches; four probes plus final commit. |
| Public standalone resolve(M) | Plan first; apply an accepted plan through the owning editor, then finish. | Keep each resolver's existing Solution policy; rejection drops the owned input. | One final integrity gate; no implicit candidate copy or recovery journal. |
| Public standalone resolve_into(&mut M) | Plan first; apply an accepted plan in one transaction, then commit. | Keep each resolver's existing Solution policy; rejection preserves the original molecule. | One final integrity gate; undo retained until acceptance, then discarded; no recovery copy. |
| Public project_into(&mut M, flags) | One transaction; plan each stage against the preceding stage's result. | Only complete Determined commits. | No outer candidate copy or nested stage journals. Probes before later planning and one commit gate. |
| Ingest and MOL parse | Pass the freshly raised M to report-free consuming resolution. | Finish, then accept only a concrete M. Underdetermination becomes a payload-free boundary error; contradiction and execution errors retain their causes. Rejected owned state is dropped. | No recovery clone, journal, or report construction; three probes plus finish. Raise's construction check remains. |
| Export/convey | Clone the retained source once; pass the candidate to consuming project; lower the result. | Return boundary representation; any failure drops the candidate. | One source-preserving candidate, with COW on changed tables. No further recovery copy or journal. |
| Transformer::transform_into | One transaction around the operation's batch/checks. | Existing success/error contract; receiver restored on rejection. | Undo instead of recovery copies. DelocalizeCharge keeps Infallible as described below. |
| Built-in transform(M) | Pass the input to the owning editor; finish after the operation. | Return the transformed molecule; failure drops the owned state. | No implicit candidate copy or recovery journal. |
| transform_iter(&M) | Create source-preserving candidates on demand; owning editor and finish for each. | Independent outputs; retain the current empty iterator on failure. | Copies for independent candidates, no recovery journal or required iterator box. |
| Reaction::apply_at and product iterators | Keep host; host.clone().edit().apply(edits), then finish. | Keep Result<Option<Molecule>, ApplyError>, including ordinary non-applicability. | One product candidate; remove discarded journal. Preserve reaction-specific correspondence. |

The source evidence is [resolve/project](../umol-graph/src/ops/resolve.rs),
[Transformer](../umol-graph/src/ops/transform.rs),
[ingest](../umol-graph/src/ingest.rs), [parse](../umol-graph/src/parse.rs),
[export](../umol-graph/src/export.rs), and
[reaction lowering/application](../umol-graph-ir/src/ir/reaction.rs).

### Resolution and projection phase mapping

Use the editor and transaction directly. Share the chemistry routines; repeat
their sequencing and outcome handling in the two entry points. There is no
BatchSession, ownership adapter, or generic phase runner.

The current composite resolution has seven batches. Preserve their exact inputs
and rejection order:

| Step | Read and plan | Apply or reject |
| --- | --- | --- |
| Opening | plan_placement and isotope.plan_resolve read the original intact molecule. | Apply placement, then isotope edits. Composite resolution accepts isotope Underdetermined edits here. |
| Constitution | probe the placed state; run valence.admit, aromaticity.select, the final atom tie-break, and constitution edit construction. | Stop at the first unaccepted outcome; otherwise apply the constitution batch. |
| Remaining entities | probe the post-constitution state. Run stereo.plan_resolve; only if Determined, compute bonds.plan_resolve and multicenter_bonds.plan_resolve from that same state. | Apply stereo, then bonds; handle the saved multicenter outcome, then apply its edits if Determined. This preserves the current failure priority. |
| Discharge | probe the resulting state; run plan_discharge. | Reject its contradiction or apply its edits. |
| Completion | Consuming resolve calls finish, then checks is_concrete on the returned molecule. Borrowed resolve_into checks is_concrete through probe, then commits only if concrete. | Return the accepted molecule or commit the change. Otherwise drop the owned state or roll back the entire transaction. |

The consuming path uses three probes and finish: four integrity checks, matching
the current four build boundaries, with no intermediate copy. The borrowed path
uses four probes and commit: five checks. The last probe and commit check the same
state; this is a concrete cost of the selected transaction interface. No validity
cache or unchecked commit is proposed to avoid it. Reporting remains opt-in.

IsotopeResolver, AromaticityResolver, StereoResolver, BondsResolver, and
MulticenterBondsResolver currently expose plan. Rename these methods to
plan_resolve, retaining their existing return types and semantics, so resolution
and projection planning use plan_resolve and plan_project. The composite Resolver
has no public whole-operation plan method. Its plan_discharge and the other
specifically named planning routines retain their names.

Valence admission and aromatic selection already return owned ResolveState values.
The final atom tie-break and constitution edit construction are currently inline
in Resolver::resolve; extract those as ordinary private routines in resolve.rs so
the two entry points share their chemistry rules. Plans returned from probe's
borrow own their data, so that borrow ends before apply changes the molecule.

Projection first checks localized/multicenter charge and spin, even for empty
flags. It then runs the selected stereo, aromaticity, valence, and isotope phases
in that order. Valence projection is currently a no-op. Each other phase already
constructs an Edits value and only then opens an editor, transacts, and builds.
Extract the planning part; the composite operation must not call those mutating
phase wrappers and create separately committed transactions.

The projection counterparts to plan_resolve are:

```rust
impl StereoResolver {
    pub fn plan_project(&self, molecule: &Molecule)
        -> Result<Edits, StereoProjectError>;
}
impl AromaticityResolver {
    pub fn plan_project(&self, molecule: &Molecule)
        -> Result<Edits, AromaticityProjectError>;
}
impl IsotopeResolver {
    pub fn plan_project(&self, molecule: &Molecule)
        -> Result<Edits, IsotopeProjectError>;
}
```

These additions let both composite and standalone projection use the same plans.
They borrow without mutation and preserve the existing planning errors. None of
these three current planning bodies produces Underdetermined or Contradictory;
their extracted result therefore needs Edits and the existing phase error, not
another Solution layer. These public methods complement plan_resolve.

Plan the first selected changing phase against the original input. For each later
phase, read the changed state through probe and apply that phase's plan. Keep all
applications in one editor or one transaction; finish/commit after the last phase.
Thus a late isotope error drops the owned candidate or undoes all preceding stereo
and aromaticity edits. Projection does not require a concrete final molecule.

Both routes retain the same phase order, plan inputs, and rejection conditions.
Repeating that orchestration adds source maintenance, not a second execution at
runtime. Consuming execution has no recovery journal or molecule copy. Borrowed
execution records undo across all phases. Ingest/export use the public consuming
entry points; export's one source-preserving copy remains at that boundary.

Standalone isotope, aromaticity, stereo, and multicenter resolution execute only
Determined plans. Bonds resolution has an infallible plan; valence admission is
immutable and needs no session. Preserve these differences from composite
resolution, which also applies partial isotope defaults.

Standalone projection follows the same ownership split where available:
project consumes the input, applies plan_project through the owning editor, and
calls finish; project_into applies the plan through a transaction and commits.
Both preserve the phase's existing error semantics. Rejection drops the consumed
input or preserves the borrowed receiver, respectively.

### Transformation-specific mapping

| Operation | Steps retained or changed |
| --- | --- |
| Aromatizer | Perceive/select before writes. Derive covered bond ids from input topology and planned systems; put additions and bond assertions in one batch. One commit/finish gate replaces editor/build followed by direct writes. |
| DelocalizeCharge | Plan all changes first; batch atom charges/assertions and same-length system electron vectors/charge. Its producer establishes shape preservation; retain Infallible and assert that producer contract. One gate replaces per-system candidate copies/checks. Journaling is extra work over raw assignments. |
| Kekulizer | Plan matchings and unchanged-field admission first; batch field changes and system removal; probe for current valence/spin checks; commit/finish. Keep present diagnostics. The known late spin failure also fails on the input; it does not prove conformant kekulization inherently fallible. |
| HydrogenFolder/HydrogenUnfolder | Use the ordinary Transformer surface; the topology, reference, and stereo rewrite remains owned by graph IR as specified in 166. Reuse participant replacement Edits. Fold replaces ligands before deleting H; unfold adds atom/bond handles before replacing virtual ligands/counts. Eligibility is planning; incidence/frame integrity is completion. |

Single-batch transformations share their plan producer between the two ownership
wrappers. Kekulizer also shares its post-execution checks. Existing chemistry
errors need no blanket From<MoleculeApplyError> requirement: an internal scope can
return a domain outcome without commit and assert only the producer-established
execution contract. Unexpected execution failure must not become a new chemistry
contradiction. This preserves Infallible and avoids forcing new error variants on
every Transformer implementation.

### Other molecule/reaction operations

extract/tracked_extract use an independent candidate and retain their compaction
result. split builds fresh components; combine/combine_all build fresh aggregates.
None needs recovery of a rejected local candidate. combine_from should replace
its unguarded move-out with recoverable append/count checkpoints; this adds
recovery work but no initial clone. inline_constraints already plans meets before
writes and needs no journal merely for centralization. lift_constraints has a
reproduced drain-before-panic defect; move its admission before writes and settle
its undetermined-stereo policy as a focused operation correction.

Normalize, reframe, remap, and canonicalize retain their established contracts.
Reaction/ReactionSpan definition editing continues through existing parts/entries
and checked construction; no current consumer requires a new ReactionEditor.
Reaction product integrity conflicts retain their existing Ok(None) versus error
classification, and product correspondence keeps its own pairing semantics.

## Additional API contracts

### Molecule mutable methods — retain current placement

Keep the current method placement and guarantees. The editor retains the complete
mutable vocabulary; these existing methods also remain on Molecule:

| Family | Signature shape and guarantee |
| --- | --- |
| atom_mut, bond_mut, dative_bond_mut, noncovalent_bond_mut; modify_atoms, modify_bonds, modify_dative_bonds, modify_noncovalent_bonds | &mut self returns a form view, or takes FnMut(Form) -> Form. Permitted writes preserve current integrity; no finish required. |
| try_modify_aromatic_system(s), try_modify_multicenter_bond(s), try_modify_stereo_atom(s), try_modify_stereo_bond(s), try_modify_constraints | &mut self plus callback -> Result<(), MoleculeIntegrityError>. Preserve the receiver on rejection; accept only integrity-valid changes. The current implementation checks a private candidate. |

Participant mutation remains on the editor's mutable entity views, as already
settled. This work does not relocate or redesign the existing Molecule methods.

### Python bindings: consuming inputs and prepared-batch transactions

| Rust operation in the design | Required Python behavior | Required migration |
| --- | --- | --- |
| Molecule::edit(self), apply(self, Edits) | Move the molecule; later owner access raises ConsumedError. Owner-backed child access raises InvalidatedViewError. Failed application does not restore the owner. | Current Python Molecule stores a live value directly; edit/apply preserve it. Root access must become fallible before this can change. |
| Editor::apply(self, Edits), finish(self) | Move the editor, including on failure; return the continuing editor or completed molecule on success. | Current editor already has consumed-state storage. Application clones Edits; change the wrapper to transfer the batch. |
| Resolver::resolve_into(&mut M) | Mutate that Python molecule on Determined; return the outcome, with no second molecule or default report. Reporting must be requested explicitly. | Current Python resolve copies and returns a molecule and report inside Determined. Its name and result shape would change. |
| Resolver::resolve(M), project(M, flags) | Consume the Python molecule under the same owner/view invalidation contract as edit/apply; no implicit recovery copy. | Use the same Molecule input-transfer change as edit/apply and the result signatures above. |
| Molecule::transact, tracked_transact | Submit prepared Edits, mutate the receiver on success, restore it on failure; tracked_transact returns the whole transaction's correspondence. | Bind the Rust conveniences over Transaction::run. No interactive Python Transaction is exposed. |

These are specific migration obligations, not authorization for a blanket
Option-based conversion of forms or entries. Doc 228 is unchanged.

Python does not expose editor probe. Remove snapshot and tracked_snapshot;
finish returns the integrity-checked molecule. Rust retains probe for multi-phase
operations; exposing intermediate whole-molecule inspection in Python has no
demonstrated consumer requirement.

Python transaction execution accepts one prepared batch or several prepared
batches. Rust opens a borrowed transaction, applies each batch separately, checks
integrity at commit, and returns success or an error. Any application or commit
failure restores the transaction-entry molecule. Separate batches retain their
own New namespaces and use the molecule at each batch's entry for Id references;
they are not concatenated. Tracked execution returns the correspondence for the
whole transaction on success.

The molecule stays mutably borrowed throughout execution. There is no recovery
clone, empty replacement, or borrowed transaction handle crossing into Python.
Python cannot inspect intermediate states, plan further batches inside the
transaction, or explicitly commit/rollback a retained handle. High-level Rust
operations can still perform their own multi-phase transactions internally.

This deliberately exposes a subset of Rust's transaction capabilities rather
than binding the interactive callback. It avoids the scoped dispatcher, handle
invalidation machinery, and scoped-tls-hkt dependency demonstrated below. Python
calls molecule.transact([edits]) or molecule.tracked_transact([first, second]).
The batch list is accepted before mutation begins; no Python iterator runs
between applications. Rust callers use the same prepared-batch methods or
Transaction::run(&mut molecule, callback) for interactive execution.

### Edits — one accumulated sequence

```rust
impl Edits {
    pub fn new() -> Self;
    pub fn push(&mut self, edit: Edit);
    pub fn add_atom(&mut self, attributes: AtomForm) -> AtomHandle;
    // Existing typed creation/update/removal methods remain.
}
```

Retain one accumulated sequence, one initial-host namespace, and one creation
namespace per kind. There is no independent-batch composition, append(Edits),
extend(Edits), or caller-offset API. Multiple batches are handled by transaction
processing under the separate-apply rules above. This is settled, not deferred.

Reaction lowering currently puts fragments already written in the final namespace
into temporary Edits and then pushes them into the final sequence. That constructs
one sequence and does not require an independent-batch composition operation.

### Error boundary

Keep the current public categories:

```rust
pub enum MoleculeApplyError {
    Transaction(TransactionError),
    Integrity(MoleculeIntegrityError),
}
```

Retain handle, old-state, missing-entry, duplicate-removal, and malformed-edit
errors. Add TransactionError::Aborted for a failed transaction that a callback
tries to continue. Transaction application and commit failures restore the
molecule before returning their original error. If the callback catches an
application or commit failure and returns success, the scope returns Aborted;
if it returns its own error, that error is preserved after restoration.

Remove RollbackFailed and RollbackStateMismatch: the private journal stays paired
with its receiver, and rollback uses the storage restoration primitives. This
depends on that ownership contract, not merely on panic freedom. No aggregate
undo validator or new mutation-mode error is added.

For the newly checked aggregate probe/completion paths, add
ResolveError::Apply(MoleculeApplyError) and
ProjectError::Apply(MoleculeApplyError), with From<MoleculeApplyError> adapters
for the scoped API. Preserve the existing phase-specific chemistry errors,
contradiction payloads, and explicitly requested ResolveReport. Integrity failures
never become Underdetermined or Contradictory. Reaction application retains its separate
product-integrity classification. This deliberately leaves the existing nested
error vocabulary intact; a general error redesign is not needed here.

Python retains TransactionError for application failures and InvalidStructureError
for integrity failures. Prepared-batch calls return the original failure; they
expose no interactive handle through which to continue an aborted transaction.

## Settled mutation vocabulary

The following contracts are settled design and are not implemented in graph IR.

### Participant mutation

Settled 2026-09-19: participant-list mutation belongs on the entity mutable view.
Its methods delegate through the owning typed set to graph-core participant
mutation so that incidence remains synchronized. A writable participant slice
alone cannot provide that contract.

The mutable overlay view holds a mutable borrow of its owning typed entity set
and the selected entity id, both private. It does not retain participant slices,
copied mutable site fields, or a separate mutable attribute reference. Attribute
and participant accessors borrow through the view for the duration of each access;
attributes and attributes_mut expose the stored payload. Mutation methods use the
selected id to delegate through the set to graph-core. Writes take effect
immediately; dropping the view performs no deferred work. Rust prevents retaining
a participant borrow across a mutation that may relocate its storage.

StereoAtomEditorViewMut and StereoBondEditorViewMut expose
ligand(position: StereoLigandPosition) -> StereoLigand. The accessor returns the
stored value by copy and panics for an out-of-range position. Callers need not
index a participant slice themselves. This follows the existing ligand(position)
name on StereoAtomView and StereoBondView; those molecule-backed views return
StereoLigandView, whereas the editor accessor requires no molecule-wide context.

Participant replacement and attribute assignment are independent: for the same
supplied replacement and attribute values, either order produces the same final
state. Replacement preserves attributes and constraints; attribute assignment
preserves participants. Either order may temporarily violate molecule integrity,
which remains the publication boundary's responsibility. This commutation does
not apply to frame-preserving permutations, which also transport attributes and
frame-relative constraints.

Separately review whether immutable entity views can use a container borrow plus
id as well, accounting for methods that require molecule-wide context. That
review is not a prerequisite for settling the mutable editor view.

All methods below take &mut self and return (). The view supplies the entity id.

| Entity mutable view | Mutation signatures |
| --- | --- |
| Aromatic system and multicenter bond | replace_atoms(atoms: &[AtomId]); replace_atom(position: AtomPosition, atom: AtomId); insert_atom(position: AtomPosition, atom: AtomId); remove_atom(position: AtomPosition) |
| Noncovalent bond | replace_atoms(atoms: [AtomId; 2]); replace_atom(position: AtomPosition, atom: AtomId) |
| Dative bond | replace_donors(donors: &[AtomId]); replace_acceptor(acceptor: AtomId); replace_donor(position: AtomPosition, donor: AtomId); insert_donor(position: AtomPosition, donor: AtomId); remove_donor(position: AtomPosition) |
| Stereo atom | replace_ligands(ligands: &[StereoLigand]); replace_site(site: AtomId); replace_ligand(position: StereoLigandPosition, ligand: StereoLigand); insert_ligand(position: StereoLigandPosition, ligand: StereoLigand); remove_ligand(position: StereoLigandPosition) |
| Stereo bond | The same ligand methods as stereo atoms; replace_site(site: BondId) uses a bond site. |

Introduce AtomPosition for non-ligand atom positions, including dative donors and
noncovalent endpoints. Retain StereoLigandPosition for ligand positions. Both are
positions within the selected entity's sequence, not entity ids. ParticipantPosition
remains graph-core vocabulary; translate to it inside the delegation rather than
expose it in these graph-IR methods. Update StereoLigandPosition documentation to
describe the current stored frame in an editor, whose length may temporarily
differ from the configuration kind's degree.

Whole replacement preserves supplied order. Single replacement preserves length
and other positions. Insertion places the new value before the supplied position;
the current length is the append position. Removal preserves survivor order.
Replacement/removal panics unless position < length; insertion panics unless
position <= length. Fixed arity is enforced by the array, with no insertion/removal
for fixed components. Each factor-specific method preserves the other factor.
Changing both uses two calls; intermediate draft inconsistency is permitted.
Do not add combined direct methods. A possible reduction in incidence updates
does not justify an additional public surface without demonstrated need.

These operations preserve entity ids, all attributes and constraints, and other
rows, while maintaining incidence immediately. Empty variable sequences,
duplicates, and references outside the current graph are admitted as intermediate
draft state; publication establishes aggregate integrity. Stereo insertion/removal
does not infer geometry or adjust configuration. No undo journal or identity
correspondence is produced by these participant operations.

Frame-preserving permutation is outside this work. It transports attributes and
frame-relative constraints as well as participants, requiring coordination beyond
the borrowed entity container. Moving higher-level permutation coordination down
can be considered later; it must not expand this scope. Supplying reordered
participants to replacement retains ordinary replacement semantics and does not
claim to preserve the represented configuration.

Localized-bond endpoint replacement is another gap. Graph has no edge-endpoint
replacement operation, so retaining BondId while rewiring a bond would require
additional graph-storage design beyond the completed relation surface.

### Edit variants

Whole replacement of each named graph-IR component is sufficient for batch
mutation. Each Edit variant has exactly id, old, and new fields:

| Variant | id type | old and new type |
| --- | --- | --- |
| ReplaceAromaticSystemAtoms | AromaticSystemHandle | Vec<AtomHandle> |
| ReplaceMulticenterBondAtoms | MulticenterBondHandle | Vec<AtomHandle> |
| ReplaceNoncovalentBondAtoms | NoncovalentBondHandle | [AtomHandle; 2] |
| ReplaceDativeBondDonors | DativeBondHandle | Vec<AtomHandle> |
| ReplaceDativeBondAcceptor | DativeBondHandle | AtomHandle |
| ReplaceStereoAtomSite | StereoAtomHandle | AtomHandle |
| ReplaceStereoAtomLigands | StereoAtomHandle | Vec<(AtomHandle, StereoLigandKind)> |
| ReplaceStereoBondSite | StereoBondHandle | BondHandle |
| ReplaceStereoBondLigands | StereoBondHandle | Vec<(AtomHandle, StereoLigandKind)> |

Graph-IR names distinguish atoms, donors, acceptors, sites, and ligands rather
than exposing graph-core's generic participant terminology. Single-position
replacement, insertion, and removal need no separate Edit variants: each can be
expressed by replacing the whole list. Changing both distinguished components uses
two edits. Execution delegates to the same storage-owned mutation used by the
entity views and preserves the other component, entity id, attributes, and
constraints.

Swapping old and new participant states defines the inverse replacement. Realized
undo records the actual previous stored participants with resolved ids, so exact
restoration does not depend on retaining batch-local handles. No separate
positional inverse algebra is needed.

Forward replacement compares the named old component exactly after resolving all
handles. Lists compare in stored order, including each stereo ligand's atom and
kind; sites and acceptors compare by id. A mismatch returns the existing
OldStateMismatch before the replacement mutates anything. Invalid handles retain
their existing errors. Delta comparison uses the reaction's frame; application
retains rule-to-host frame alignment rather than requiring equal raw storage order
between independently framed rule and host. Undo restores its saved value without
an expected-post-value comparison.

### Delta and Undo variants

Use nine factor-specific Delta variants and nine corresponding Undo variants.
Every variant carries the affected entity's typed id. Delta variants additionally
carry old and new values; Undo variants carry only the named saved value, with
the same type shown in the table.

| Delta variant | old and new type | Undo variant | Saved field |
| --- | --- | --- | --- |
| AromaticSystemDelta::ReplaceAtoms | Vec<AtomId> | RestoreAromaticSystemAtoms | atoms |
| MulticenterBondDelta::ReplaceAtoms | Vec<AtomId> | RestoreMulticenterBondAtoms | atoms |
| NoncovalentBondDelta::ReplaceAtoms | [AtomId; 2] | RestoreNoncovalentBondAtoms | atoms |
| DativeBondDelta::ReplaceDonors | Vec<AtomId> | RestoreDativeBondDonors | donors |
| DativeBondDelta::ReplaceAcceptor | AtomId | RestoreDativeBondAcceptor | acceptor |
| StereoAtomDelta::ReplaceSite | AtomId | RestoreStereoAtomSite | site |
| StereoAtomDelta::ReplaceLigands | Vec<StereoLigand> | RestoreStereoAtomLigands | ligands |
| StereoBondDelta::ReplaceSite | BondId | RestoreStereoBondSite | site |
| StereoBondDelta::ReplaceLigands | Vec<StereoLigand> | RestoreStereoBondLigands | ligands |

Delta inversion swaps old and new. Undo captures the actual previous stored value
with resolved ids and exact sequence order, then restores it through ordinary
storage replacement. Local target-access guards provide the settled no-panic undo
contract; undo does not validate an expected post-value. These operations preserve
the other component, attributes, and constraints. No positional or combined-component
variants are needed.

The replacement undos do not use graph-core restore or restore_participants:
those undo row removal and reference compaction, respectively. These new undos
restore an arbitrary saved component on an existing row without changing its id.

Extend the existing Delta normalization and composition rules to the named
components, rather than introducing another change algebra:

| Sequence | Folded result |
| --- | --- |
| A to B, then B to C | A to C |
| A to B, then C to D, with B unequal to C | Existing contradiction result |
| A to A | Identity change eliminated |
| Add followed by replacement | Final component absorbed into Add |
| Replacement followed by removal | Remove carries the original entry |
| Add, changes, then removal | Cancel as in the existing created-entity fold |

Component continuity uses exact sequence equality. Do not sort component lists
or implicitly transport attributes during this folding. Context-free Deltas
normalization retains replacement variants; lowering to complete removal/addition
entries uses the reaction lhs to supply the complete entity state.

### Reaction-span representation

A replacement Delta may target an existing entity id without requiring that entity
to be preserved in the resulting reaction span. Follow the existing stereo-kind
precedent: changes incompatible with a Modified entry's shared incidence are
represented as removal plus addition, with each entry carrying its own frame and
attributes. The span correspondence leaves these removed and added entities
unpaired. Editor identity preservation does not impose reaction-span identity
preservation.

Keep ReactionSpan's shared-incidence representation. Do not introduce separate
left/right participant lists within one preserved span entity merely to support
replacement deltas. Lower the new Delta variants into the existing span
representation, including the resulting correspondence. This does not broaden
the deferred permutation work.

Current implementation distinction: correspondence induction leaves stereo
entities of different determined kinds unmatched, and superimposition represents
them as Removed and Added. An incompatible Modified span is rejected. Reaction
application also rejects a configuration ModifyField that changes determined
kind; it does not automatically rewrite that modification into removal/addition.
The existing precedent establishes the representation and its consumers, not a
general automatic conversion of arbitrary replacements.

Settled conversion approach: materialize the before/after states for replacement
deltas and reuse the existing correspondence induction, superimposition, and
reference-mapping operations. Shared-incidence and stereo-kind compatibility use
the existing rules; incompatible entities become removal/addition entries. Span
conversion back to a reaction may therefore express the replacement as removal
and addition rather than recover the original Delta spelling.

Incidence-changing replacements are expected to be rare. Favor this straightforward
path over specialized change-history tracing or optimization without demonstrated
need. Connecting replacement deltas to these existing operations is implementation
work, not a requirement for another reaction representation or public API.

### Local getters — existing surface and proposed additions

Provide the same local read conveniences on immutable and mutable editor entity
views, using the molecule-backed immutable counterparts' names and meanings. Participant
inspection should not require callers to index slices or unpack attributes merely
because the entity is being edited.

MoleculeEditor already provides atom, bond, dative_bond, aromatic_system,
multicenter_bond, noncovalent_bond, stereo_atom, and stereo_bond, their mutable
counterparts, entity counts, and constraints access. Those entry points are not
additions. The additions below are methods on the returned editor entity views,
both immutable and mutable; they are not yet implemented.

All current editor entity views expose id and attributes as public fields.
Mutable views hold a mutable attributes reference. Their remaining local read
surface and the proposed additions are:

| Entity kind | Already available on editor views | Proposed getter additions |
| --- | --- | --- |
| Atom | id and attributes fields only | element, isotope_mass, charge, implicit_hydrogens, lone_pairs, unpaired_electrons |
| Localized bond | atoms field | atom_ids, order, charge, unpaired_electrons |
| Dative bond | atom_ids() | donor_ids, acceptor_id, donor_count, atom_count, order |
| Aromatic system and multicenter bond | atom_ids() | atom_count, electrons, electron_count, charge, unpaired_electrons |
| Noncovalent bond | atoms field | atom_ids, kind |
| Stereo atom and stereo bond | site and ligands fields | site_id, ligands, ligand_count, ligand, ligand_position, ligand_frame, ligand-kind filters, ids and counts; configuration accessors |

The proposed views replace public field access with id() and attributes(); mutable
views additionally provide attributes_mut(). These methods are also additions.
References returned by mutable-view accessors are tied
to the accessor borrow. atom_ids() returns [AtomId; 2] for localized and noncovalent
bonds; ordinary variable atom/donor collections use exact-size iterators of AtomId.
Stereo site_id() returns AtomId or BondId according to the entity kind. ligands()
returns an exact-size iterator of copied StereoLigand values. ligand_frame()
returns the owned ordered frame. ligand_position(atom: AtomId) returns
Option<StereoLigandPosition> for the first matching actual-atom ligand, following
the existing getter's lookup semantics. Ligand-kind filters return stored ligand
values; their id and count companions keep the existing names and meanings.

Stereo configuration() returns &StereoConfigurationForm, kind() returns
Option<StereoKind>, and coset() returns Option<&StereoCoset>. The optional results
follow the underlying configuration accessors and admit an unfinished undetermined
configuration without a panic. Other attribute getters borrow the corresponding
stored forms. Include local payload predicates such as is_ground and
is_undetermined where their existing meanings apply. constraints() returns a
reference to the stored entity constraint container, not a molecule-backed
constraints view.

Methods resolving other entity views or deriving topology-dependent quantities
(neighbors, aromatic system bonds, valence, or stereo-bond site endpoint atoms)
need context beyond the owning overlay set. Identify these separately rather than
expanding the new mutable view's ownership merely for mechanical parity. These
context-dependent methods remain outside the settled local getter surface.

## Storage integration and implementation obligations

Replace FixedSetStorage, VarSetStorage, and FixedVarSetStorage with the existing
typed entity sets. MoleculeEditor owns Molecule directly; no second mutable-vector
representation, OverlayEditor trait, or new public construction seam is needed.
Graph core owns incidence maintenance, compaction, and restoration; graph IR owns
attributes, constraints, molecular frames, and cross-entity coordination.

| Typed set | Existing relation storage | Factors |
| --- | --- | --- |
| AromaticSystems, MulticenterBonds | VarRelationSet | Atoms |
| NoncovalentBonds | FixedRelationSet | Two atoms |
| DativeBonds | FixedVarBirelationSet | Acceptor; donors |
| StereoAtoms | FixedVarBirelationSet | Atom site; ligands |
| StereoBonds | FixedVarBirelationSet | Bond site; ligands |

### Storage delegation

Retain graph-IR coordination while replacing the editor's storage implementations.
Overlay operations delegate through the typed entity sets that own copy-on-write
relation storage. This inventory records the implementation obligations for the
API defined above; it is not a staged implementation plan.

| Editor surface | Current implementation and required change |
| --- | --- |
| add_atom, add_bond | Already delegate to Graph::add_node/add_edge. Retain parallel attribute-array updates; direct additions record no correspondence. |
| add_dative_bond, add_aromatic_system, add_multicenter_bond, add_noncovalent_bond, add_stereo_atom, add_stereo_bond | Replace mutable-entry-vector insertion with relation add through the typed set; direct additions record no correspondence. |
| atom_mut, bond_mut | Already mutate copy-on-write attribute arrays. Retain graph-IR ownership of these arrays. |
| Mutable attribute views for all six overlay kinds | Replace entry-vector materialization with typed-set access to relation payload mutation. Keep attributes and participants accessible in the entity view. |
| push_constraint, constraints_mut, inline constraints through attribute views | Retain graph-IR mutation; graph-core does not interpret molecular constraints. |
| Each overlay remove_* / tracked_remove_* pair | Use relation tracked_remove and its returned row compaction instead of editor-owned row removal and separately constructed compaction. Retain constraint compaction; remove public tracked direct mutation in favor of tracked batch execution. |
| Topology remove / tracked_remove | Retain Graph::tracked_remove_cascading; replace wrapper-specific overlay compaction with relation tracked_compact through typed sets. Retain atom/bond attribute compaction, cascade coordination, and molecule-wide constraint updates. Public tracking belongs to batch execution. |
| Internal undo-addition removal and restore_* methods | Route removal through the same primitives. Replace graph/relation reconstruction with storage restore_participants and restore in the appropriate order. Retain attribute and constraint restoration with the settled local guards below. |
| apply, transact, and tracked counterparts | Replace the current lifecycle with the editor/transaction API above. Both execution paths share graph-IR mutation operations; handle resolution and forward preconditions remain batch concerns, and undo capture/replay remains transactional. |
| snapshot, try_build, build, and tracked counterparts | Remove relation-row rebuilding at publication. Replace editor publication with finish, Molecule::apply, or scoped commit, and replace snapshot with checked probe access. Fresh MoleculeBuilder retains asserted build. |
| Overlay reads and views; six internal *_equiv methods | Remove storage-wrapper dispatch. Use typed-set accessors and graph-core participant comparison where applicable; retain graph-IR frame transport and payload comparison. |

Plain editor removal still needs compaction internally even when it returns no
witness. Atom/bond addition already delegates to Graph, and entity-set frame
transformations already use storage participant permutation. Undo of additions
can reuse removals; undo of removals needs the restoration capability below.

Existing offered-old comparisons align participants and transport payloads into
the stored frame. Preserve the direction to[i] = from[action[i]]; ordinary
unordered factors use DynPermutation and stereo uses bounded Permutation.
Stereo-bond alignment may reorder within endpoint blocks or exchange whole
blocks, but cannot move one ligand between blocks. Their internal placement is
not part of the settled public API; no public comparator or generic trait is
proposed.

### Removal rollback: current implementation

The current editor saves removed entries, compactions, and cascaded constraints,
then reconstructs graph/relation rows during undo. validate_undo precedes those
reconstruction methods. The changes below replace this reconstruction and checker
with storage restoration and local access guards; they are not yet implemented.

### Restoration and local guards

Rollback follows the same history-bound contract as graph-core restoration:

- Rolling back a transaction from its matching post-state restores the pre-state
  modulo normalized_eq, consistent with forward old-value matching. Structural
  equality of non-normal form encodings is not required.
- Mismatched or manipulated history must not cause a panic. The resulting state
  is unspecified; rollback need not detect misuse, report it, or leave the editor
  unchanged.

An accepted equivalent old value may be retained in Undo; capturing its original
encoding is not required. Participant sequence comparison and graph-core storage
restoration retain their exact contracts.

This changes the current structural-rejection contract. Remove
RollbackStateMismatch and the Result used solely to report rollback mismatches.
RollbackFailed is likewise unnecessary once internal undo follows this contract.
Forward edit errors and old-state checks remain: forward edits are independently
supplied requests, whereas undo restores changes already recorded by the operation.

Remove validate_undo without introducing another aggregate checker. Put necessary
guards directly in the methods that consume the saved data, and delegate guards
already owned by graph-core to its restoration operations:

| Consumer | Responsibility |
| --- | --- |
| Graph and relation restoration | Delegate to storage-owned restore methods and their panic guards. |
| Atom/bond attribute restoration | Guard sizes, indices, and slot access in the attribute restoration methods. |
| Dative restoration | Guard acceptor extraction from a saved participant sequence. |
| Undo of additions and field changes | Guard target access in the consuming removal or field-restoration method. Restore recorded field values without requiring equality with the recorded post-value. |
| Constraint restoration | Guard sizes, positions, and iteration locally; restore recorded values without requiring equality with the recorded post-value. |
| Correspondence on failure | Return no correspondence on failure. Editor tracked_apply computes its batch mapping; transaction tracked_commit derives the overall mapping from the journal. Neither requires an unconditional session accumulator. |

Do not retain blanket checks of all entity counts, equality between forward and
inverse compaction records, or external graph membership of saved overlay
participants. Storage delegation removes the editor's relation-row reconstruction
checks. Attribute and constraint restoration must still protect their own accesses;
they do not need to prove that the supplied history matches the editor.

The consuming implementations must become panic-free before validate_undo is
removed: deleting the checker alone would expose existing indexing, unwrap, and
expect paths. These restoration semantics are implementation prerequisites for
the public transaction lifecycle above, not current implementation behavior.

### Cost and memory decisions

The comparisons below separate planning/output costs from recovery costs.
They do not treat every allocation or existing Arc detachment as an editor copy.

| Route | Additional recovery storage | Success work | Failure/cancellation work |
| --- | --- | --- | --- |
| Owning direct editor | None. | Direct kernels and final gate; checks at requested probes. | Drop owned storage; no recovery. |
| Owning batch apply | None. | Batch realization/preconditions and mutation; probe/finish gates. Molecule apply includes finish. | Drop editor and successful prefix. |
| Borrowed transaction | Undo entries plus removal/cascade metadata across all phases. | Batch realization, mutation, checked probes, final gate; discard journal on acceptance. | Reverse the whole operation. |
| Independent product | One output candidate, required because source survives success. | Owning direct/batch route plus COW from source sharing. | Drop candidate; source was never changed. |

Batch realization is not free. Currently each ApplicationState allocates initial
handle vectors for all entity kinds, and every editor allocates full session
correspondence. The design removes unconditional session correspondence. Keep
initial batch handles as identity-by-count until a removal requires translation;
keep created-handle vectors only for actual creations. At the first affected
compaction, materialize the required mapping for that entity kind and update it.
Thus field-only batches need no receiver-sized handle vector; removal can still
require receiver-sized work. This preserves handle semantics without pretending
that compaction and incidence maintenance are proportional only to removed rows.

The existing measurements below give a useful sparse-change comparison: at 88
atoms the copied atom/bond tables requested about 16.5/7.6 KiB, while an isolated
Undo entry occupied 752 bytes. They do not establish that the journal is always
smaller. At the current enum size, 88 separate field records already require
66,176 bytes for entries alone, before nested payloads and compactions. That is an
arithmetic illustration, not a measured full resolver workload. Dense resolution
and kekulization must not be described as necessarily cheaper than copying.

Use the existing Undo representation initially so the public mutation/restoration
contract is not coupled to a new compressed journal design. Its large enum layout
is a recorded performance liability. Avoid avoidable payload clones by consuming
Edits and moving old values where their existing comparison contract permits;
do not reconstruct a full molecule to save an undo. An implementation benchmark
must report actual journal entry count and peak retained bytes for dense resolver
fixtures before making a performance claim. This is a validation obligation, not
another speculative benchmark campaign or an automatic recovery-strategy switch.

Likewise, relation participant replacement/add/remove currently rebuild incidence,
and Graph additions rebuild CSR. Undo of removal can require whole-set compaction
and restoration work. Those costs remain under both ownership routes. The design
selects journaled recovery to avoid mandatory whole-receiver recovery copies for
borrowed recovering operations, and selects no journal for destructive editing or
independent products. No speedup over raw direct mutation is claimed for a
transaction; its additional work buys recovery.

### Failure, integrity, and unwind implementation obligations

Retain existing edit precondition errors and the two MoleculeApplyError categories
for application and integrity. Add the latched TransactionError::Aborted needed
by the scoped lifecycle. There is no mutation-path error. Remove the obsolete
rollback-mismatch/error protocol as already settled above. Undo remains public;
the transaction's journal, state, guard, and any internal partial-edit records are
private and cannot be substituted by the caller.

Field execution accepts normalized equality with the offered old value. Reversing
that accepted change is consistent with rollback modulo normalized_eq: restoring
a literal in place of an equivalent singleton literal set is permitted. Retain
this comparison and undo behavior; no extra capture of the original encoding is
required. This is a guarantee for the operation's own matching history, not a
correctness guarantee for manipulated journals.

The scope alone is not enough for panic recovery. Current transact records Undo
after apply_edit_with_undo returns. AddAtoms/AddBonds can mutate several rows before
that happens; topology removal also mutates before assembling all undo metadata.
The implementation must change this ordering:

1. Resolve handles, check offered old state and edit shape, and obtain all fallible
   preparation data before the corresponding writes.
2. Reserve journal capacity and prepare restoration data before mutation. For a
   multi-row operation, register progress at each completed storage step, so an
   interruption can restore its completed prefix as well as earlier edits.
3. Keep the assignment and inverse registration free of intervening fallible calls;
   otherwise register the inverse before assignment.
4. For compaction across graph, attributes, relations, and constraints, capture
   removed data and the required maps before their respective steps; update the
   private progress record before another potentially panicking step. Delegate
   restoration to the settled graph-core operations.
5. Audit the concrete graph-IR participant/form implementations used by these
   kernels. Generic graph-core callbacks are not evidence that arbitrary user code
   is present in this fixed graph-IR path. A primitive must either leave storage
   intact on unwind or have a recovery record covering its partial write state.
6. Finish correspondence extraction and other return bookkeeping before marking
   the private guard accepted. A commit request must not disarm restoration while
   an operation-owned step can still unwind.

A per-apply guard restores and latches Aborted if the batch unwinds, including
when that panic is caught inside the callback. The outer scope guard covers
callback abandonment, error, and unwinding after commit was requested. Ordinary
unwinding out of the scope restores before the caller can read the receiver.
Do not catch and convert panics into chemistry errors. Process abort, including
allocation failure that aborts the process, has no return-state
guarantee. Capacity/index failures reachable from edit inputs are checked before
writes; library-issued rollback must not itself panic. This is necessary work for
the borrowed design, not a property established by today's implementation.

probe failures are representation failures, not chemical contradictions. Final
acceptance additionally uses each consumer's existing model-dependent condition.
Tracked and untracked execution have the same acceptance and mutation result.
Owned and borrowed execution have the same accepted molecule when given the same
initial state and plans; their intentional difference is what survives rejection.

## Implementation verification

Preserve the existing laws for handle namespaces, apply/transact result agreement,
rollback, failure recovery, constraint compaction, publication, and correspondence
composition. Relevant suites include
[edit properties](../umol-graph-ir/tests/property/edit.rs),
[publication](../umol-graph-ir/tests/property/molecule/publication.rs),
[compaction](../umol-graph-ir/tests/property/molecule/compaction.rs), and the
[reaction properties](../umol-graph-ir/tests/property/reaction).

For rollback, state the matching-history restoration law using normalized_eq,
including accepted old values with different equivalent encodings. Do not require
structural equality of those encodings. Update misuse coverage to require freedom
from panics, without asserting a particular mismatch
error, restored result, or unchanged editor for manipulated inputs.

Extend coverage from the settled contract: all six overlay kinds plus atoms and
localized bonds, participant and payload coordination, legal and illegal frame
alignment, cascade restoration, reaction reference updates, and Python failure
ownership. Include position-sensitive payloads so that frame tests prove transport.
Benchmark representative direct,
apply, and transact paths from the beginning, including copy-on-write sharing and
repeated relation mutation; do not assume the current fixtures establish scale.

## Supporting evidence

The results below preserve the evidence needed after deleting scratch. Baseline
API names and behavior describe the measured code, not the selected interface.
Unselected experiments are identified explicitly; none adds implementation scope.

### Python scoped transaction binding — 2026-09-22

**Feasible alternative, not selected.** The selected Python interface submits
prepared batches and keeps transaction handling entirely inside Rust. This
experiment's scoped access machinery and dependency are not implementation
requirements; the results are retained as evidence for the choice.

A standalone PyO3 0.29.0 module ran under Python 3.13.15: **35 pytest cases
passed**, and strict Clippy passed. It binds a Rust transaction borrowing a model
molecule (a value vector plus declared count), with real Python callbacks and
exceptions. It does not yet bind graph-IR transactions or the full Molecule read
surface. A separate negative compilation confirmed that a borrowed transaction
cannot be stored directly in a lifetime-parameterized pyclass.

The binding uses [scoped-tls-hkt 0.1.5](https://docs.rs/scoped-tls-hkt/0.1.5/scoped_tls_hkt/)
to register a stack-local dispatcher during the callback. Python handles contain
only a unique scope identifier and thread identifier. The Rust transaction stays
in a stack-local RefCell<Option<Transaction>>; the library owns the recovery
guard and journal. The experiment forbids unsafe code; the dependency encapsulates
the scoped-reference machinery. No molecule or batch copying occurs.

The tests cover multiple batches, commit, rollback, uncommitted return, partial
application failure, failed integrity checks, repair after a rejected probe,
exceptions and injected Rust panics (including after commit request), caught
application failures that still abort, escaped handles/views, source aliases,
and nested transactions on different molecules. Probe holds a shared borrow:
mutation and finalization are rejected during its callback. Source aliases are
blocked by PyO3's exclusive borrow. Scope exit invalidates retained handles;
later scopes cannot reactivate them. Ordinary Python exception identity survives.

For model inputs containing 2 and 100,000 values, the same three writes recorded
three undo entries on success, rollback, and exception. The original vector's
data address stayed unchanged. On this machine each Python handle's Rust payload
is 16 bytes, the stack transaction slot is 32 bytes, and each model undo is 24
bytes. These exclude Python object overhead, journal capacity, and dispatcher
stack storage; they are not total-memory or runtime benchmarks.

The costs are a new dependency, scoped dispatch and runtime borrow checks per
call, and Python handle allocation. Handles work only on the callback thread;
moving an active handle to another thread is rejected. Python checked reads use
probe(callback), not a freely retained Molecule reference. This demonstrates a
feasible binding without recovery copies; it does not establish that this is the
only or simplest implementation.

Reproduction used scratch/213-python-transaction after activating umol-py/.venv:
cargo run --offline --manifest-path scratch/213-python-transaction/Cargo.toml
--target-dir /Users/dr/.cargo-target/213-python-transaction
--bin python-transaction-probe (the executable registers the module and runs
pytest). Clippy used the same manifest, target directory, and binary with
-- -D warnings. Checking --bin borrowed_class intentionally fails with
“#[pyclass] cannot have lifetime parameters.” These results are retained here so
the scratch directory can be deleted.

### Checked-borrow lifetime check — 2026-09-22

A standalone Rust model checked the proposed probe(&self) -> Result<&Molecule, E>
signature with rustc 1.96.0, edition 2024. Planning through the reference and then
applying after its last use compiled for both the owning editor and a scoped
borrowed transaction. Three negative cases were rejected: returning the reference
from the transaction callback (E0515 and an incompatible lifetime); applying while
the reference is still in use (E0502); and consuming the editor while its reference
is still in use (E0505). No adapter or callback-based probe was needed.

The commands used rustc --edition=2024 --emit=metadata on
scratch/213-probe-borrow/main.rs, with separate escape, transaction_write, and
editor_move cfgs for the negative cases. This checks accessor lifetimes only;
integrity checking and rollback were not implemented in that model. The result
supports replacing the probe callback without changing the integrity barrier.
These recorded results are sufficient to delete the scratch files.

### Ownership and execution prototype — 2026-09-22

A standalone std-only prototype compiled without warnings using rustc 1.96.0.
Its model molecule owns a value vector and declared count; complete state was
compared in each case. The selected shape has one owning MoleculeEditor and a
transaction borrowing the existing molecule. Both use the same private edit
kernel. There is no copied draft, move-out, empty placeholder, or ownership mode.

Eleven transaction cases passed: accepted commit; ordinary domain rejection;
explicit rollback; later batch error; failed final integrity; rejected probe
followed by repair; forgotten capability; callback error after commit request;
callback panic after commit request; rollback after an application failure without
clearing the abort latch; and a two-write edit panic caught inside the callback.
The last case restores before another probe, rejects commit, and returns Aborted.
It established the need for per-apply cleanup as well as the outer scope guard.

An earlier code-sharing experiment ran three paired cases on Transaction and an
Option<MoleculeEditor> adapter: two dependent successful batches, domain rejection
after the first, and a second-batch precondition error. Both accepted paths produced
values [2, 7, 9] and declared count 3. Rejection restored the borrowed input and
discarded the owned candidate. Destructive execution recorded zero undos. The
adapter was subsequently rejected; these results establish only the tested
ownership behavior, not the selected implementation structure.

Five owned-editor cases passed: direct/batch/direct interleaving; repair by a batch
after a direct change failed probe; dropping a directly changed editor; batch
failure discarding earlier direct changes; and final integrity failure discarding
the input. Whole-model Clone instrumentation counted zero calls across all cases.
Fourteen borrowed transaction entries incurred zero measured scope-setup
allocations, excluding batch execution, journal growth, observation, and mutation.
That limited measurement establishes no production transaction timing.

Four compile-fail attempts were rejected: returning Transaction from its scope;
returning a borrowed slice through probe; accessing the source during its mutable
borrow; and calling transaction on the owning editor. The first two failed for
lifetimes, source access with E0502, and the absent editor method with E0599.

Commands used rustc --edition=2024 on scratch/213-owned-editor/main.rs, execution
of the binary, and separate compilations with escape_transaction, escape_view,
source_access, and editor_transaction configuration flags. The scratch source and
logs are disposable; the cases, results, and limits above are the durable record.
The model establishes ownership, cancellation, and shared phase execution. It
does not establish production graph/relation interruption safety, chemistry
correctness, performance, tracking behavior, or Python feasibility. No production
source or workspace test was changed for the experiment.

### Consumer audit probes — 2026-09-21

Three public-API consumer probes passed: input spin rejection equals kekulizer's
late contradiction; projection mutates earlier stages before isotope rejection
while preserving the composite receiver; resolution performs nonempty isotope
planning before later underdetermination and preserves its input. The separate
lift_constraints probe constructed its input successfully and observed
`lift_panicked=true, molecule_unchanged=false`. A standalone safe-Rust guard probe
confirmed the mem::forget counterexample.

The consumer probe reused these committed fixtures at revision 1fd0c3aa7fbda4e47ec122d0222dd5829ab772a9:

| Probe | Reproduction input and configuration | Additional observation |
| --- | --- | --- |
| Kekulization | kekulizer.rs, test_kekulizer_transform_into_error::spin_invariant; default KekulizeConfig, atom order 0..6. | Run SpinInvariantsValidator on the input; its contradiction exactly equals the payload of PostLocalizationSpinInvariant. |
| Projection | resolve.rs, test_resolver_project_phase_error::isotope; default ChemistryModel, Natural isotope policy, ProjectFlags::all(). | Run stereo and aromatic projection separately first: both return Determined and remove their respective overlays. Isotope projection then reports NonGroundIsotope at AtomId(0). Composite projection leaves the source equal to its original. |
| Resolution | resolve.rs, test_resolver_resolve_later_underdetermined; default ChemistryModel and Natural isotope policy. | IsotopeResolver::plan returns a nonempty Determined Edits; complete resolution returns Underdetermined and preserves exact input. |

The lift probe used Molecule::try_from_entries with one default AtomForm plus
AtomConstraintForm::degree(0), one stereo entry `(AtomId(0), vec![],
StereoAtomForm::default())`, and all other entries default. catch_unwind around
lift_constraints observed the panic and unequal receiver. The guard probe used a
borrowed state with a degree matching its vector length, changed degree through
the guard, then forgot the guard whose Drop would have restored it; reading the
state afterward observed the mismatch.

The probes used only standalone scratch crates/programs; no production tests,
Python build, workspace gate, or new benchmark campaign was run. The earlier cost
measurements below remain the performance evidence. Their conclusions and the
probe inputs/results are recorded here so deleting scratch loses no design result.
The replacement editor itself was not implemented or benchmarked by these probes.

### Feasibility measurements — 2026-09-20

#### Decision supported

**Interpretation revised 2026-09-21:** these successful-path measurements support
the potential cost savings of removing ownership copies and journal work, but
do not override the revised recovery contract. The consuming variant discards
input on failure and is not an equivalent replacement for an operation promising
unchanged input on error.

Across these fixtures, projection median time fell 28–37% with
uniquely owned input and 22–28% with a retained shared source. Requested allocation
bytes fell 58–69% and 50–55%, respectively. Removing journals alone reduced
projection time 16–25%; consuming ownership produced an additional reduction.
The benefit therefore survives the intentional source copy used by export.

Resolution gains are smaller: 3–7% median-time reduction for unique input and
1–4% for shared input, with 4–9% and 2–4% fewer requested bytes. The composite
resolver already uses journal-free application; its changes here remove
intermediate ownership copies. These small timing differences are directional
observations from a local feasibility run, not guaranteed speedups. They support
passing ownership through resolution without justifying additional special cases.

No change to chemistry selection, projection meaning, or final acceptance
criteria was needed for the successful workloads. That does not justify giving
up recovery. The selected design uses destructive consuming methods and journaled
borrowed methods. The clone measurements below document the alternative's costs;
they do not leave the recovery mechanism undecided.

#### Method and scope

Source baseline: 1c2d0482dec5f459fc5d75e51ecd04f3e84e6b68. An isolated source copy
under scratch held the experiment; production Rust sources were not changed.
Environment: macOS arm64, rustc 1.96.0 (ac68faa20, 2026-05-25), Cargo release
profile, offline dependencies. Only the graph dependency closure was built;
Python and workspace-wide tests were not involved.

The three variants were:

1. **Current:** the existing resolver and projection implementations.
2. **No journal:** projection retains current cloning/ownership and replaces each
   standalone phase's transact(edits) with consuming apply(edits), followed by
   the same checked publication. This is an isolation measurement, not another
   proposed public API. There is no corresponding resolution variant because
   the aggregate resolver already uses apply.
3. **Consuming:** move all Molecule fields into the existing editor constructor;
   compute each phase's plans before moving the input, then move accepted results
   between phases. Projection additionally consumes its initial candidate and
   uses journal-free phase application. No recovery clone or empty-placeholder
   transfer is included in the final measured variant.

All variants retain the same chemistry algorithms, edit preconditions, aggregate
integrity checks, and phase boundaries. Current editor storage wrappers and eager
session-correspondence construction remain, so these are measurements of the
ownership/journal changes, not of the complete storage-delegation redesign.
Publication retains the current checked/asserted calls; the draft's Result
adapters and active Transaction type are not implemented by this experiment.
Hydrogen transformations have no implementation to measure yet.

**Unique** means each input was independently parsed/raised for resolution or
independently ingested for projection. **Shared** means the input was cloned from
a retained source. Input creation, input cloning, parsing, and resolver setup are
outside timing. Output destruction is included identically across variants.
Thus the shared rows measure the operation after an intentional source copy,
as retained for export; they do not include the cost of making that outer copy.
They do not prescribe implicit copying in the Python bindings.

Each row is the median of 11 batch means. Batches contain 8–256 operations,
selected once per workload/ownership pair toward approximately 2 ms from a warm
baseline operation. Variant order rotates between samples. There are no
CPU-affinity or machine-load controls. A counting system allocator is enabled
only for a separate untimed operation after warm-up. Allocation calls include
reallocations; requested bytes sum requested allocation/reallocation sizes.
These are allocation-traffic measurements, not peak live memory or process RSS.
No expensive precision reruns were used to refine small resolution differences.

Model: ValenceModel::smiles(), default chemistry model otherwise,
IsotopePolicy::Natural, and ProjectFlags::all(). Fixtures:

| Fixture | Input SMILES |
| --- | --- |
| octane | CCCCCCCC |
| chain64 | 64 consecutive C characters |
| naphthalene | c1ccc2ccccc2c1 |
| stereo | `N[C@H](F)/C=C/C` |
| combined | `[13CH3][C@H](F)/C=C/c1ccccc1` |
| combined8 | Eight copies of combined, separated by dots |

combined8 is a synthetic disconnected scaling case, not a claim about production
molecule sizes. The fixture selection covers ordinary chains, aromaticity,
stereo, and their combination; it is not a corpus-weighted throughput estimate.

#### Correctness checks

Before timing, the consuming resolver matched the current resolver's complete
Molecule output and report for all six fixtures. Both projection variants matched
the current complete Molecule output. Additional checks matched wildcard
underdetermination, a late discharge contradiction for C#c0#h4#n0#u0#s#D5, and a
late isotope-projection error after earlier stereo/aromatic phases. The baseline
projection preserved its source on that error; the consuming operation returned
the same diagnostic and discarded its owned state, as intended. These checks
establish equivalence for the measured cases, not general algorithm verification.

#### Resolution measurements

| Fixture | Input ownership | Current µs | No journal µs | Consuming µs | Requested KiB, current → consuming | Allocation calls, current → consuming |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| octane (8 atoms) | unique | 13.32 | — | 12.52 | 39.4 → 35.6 | 91 → 85 |
| octane (8 atoms) | shared | 12.81 | — | 12.29 | 39.4 → 37.8 | 91 → 89 |
| chain64 (64 atoms) | unique | 95.92 | — | 91.14 | 346.2 → 316.7 | 411 → 405 |
| chain64 (64 atoms) | shared | 94.89 | — | 93.51 | 346.2 → 334.2 | 411 → 409 |
| naphthalene (10 atoms) | unique | 36.54 | — | 34.95 | 86.5 → 80.0 | 644 → 614 |
| naphthalene (10 atoms) | shared | 35.54 | — | 34.74 | 86.5 → 83.3 | 644 → 629 |
| stereo (6 atoms) | unique | 15.23 | — | 14.15 | 44.8 → 41.5 | 194 → 184 |
| stereo (6 atoms) | shared | 14.89 | — | 14.26 | 44.8 → 43.2 | 194 → 189 |
| combined (11 atoms) | unique | 45.15 | — | 43.83 | 104.4 → 97.8 | 889 → 867 |
| combined (11 atoms) | shared | 44.87 | — | 43.79 | 104.4 → 101.1 | 889 → 878 |
| combined8 (88 atoms) | unique | 402.01 | — | 387.82 | 1297.6 → 1245.8 | 6378 → 6258 |
| combined8 (88 atoms) | shared | 391.66 | — | 384.16 | 1297.6 → 1271.7 | 6378 → 6318 |

#### Projection measurements

| Fixture | Input ownership | Current µs | No journal µs | Consuming µs | Requested KiB, current → consuming | Allocation calls, current → consuming |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| octane (8 atoms) | unique | 1.94 | 1.54 | 1.32 | 10.9 → 3.4 | 13 → 8 |
| octane (8 atoms) | shared | 1.59 | 1.19 | 1.14 | 10.9 → 4.9 | 13 → 10 |
| chain64 (64 atoms) | unique | 12.77 | 9.84 | 8.48 | 95.3 → 35.2 | 20 → 15 |
| chain64 (64 atoms) | shared | 11.42 | 8.68 | 8.48 | 95.3 → 47.3 | 20 → 17 |
| naphthalene (10 atoms) | unique | 8.87 | 6.84 | 5.61 | 61.7 → 26.1 | 111 → 61 |
| naphthalene (10 atoms) | shared | 8.20 | 6.33 | 5.94 | 61.7 → 30.5 | 111 → 77 |
| stereo (6 atoms) | unique | 5.97 | 5.00 | 4.33 | 17.9 → 6.6 | 76 → 53 |
| stereo (6 atoms) | shared | 5.42 | 4.53 | 4.25 | 17.9 → 8.4 | 76 → 58 |
| combined (11 atoms) | unique | 12.80 | 10.60 | 8.84 | 54.6 → 20.9 | 162 → 106 |
| combined (11 atoms) | shared | 12.16 | 10.11 | 8.99 | 54.6 → 24.9 | 162 → 118 |
| combined8 (88 atoms) | unique | 88.15 | 72.75 | 61.26 | 436.1 → 183.3 | 659 → 414 |
| combined8 (88 atoms) | shared | 82.96 | 68.50 | 61.73 | 436.1 → 214.9 | 659 → 482 |

#### Remaining costs

Molecule clones increment Arc counts for graph and entity storage, but copy the
owned global Constraints vector. Subsequent mutations can copy shared tables;
consuming ownership avoids that when no other owner remains. It does not remove
algorithmic work, storage-index rebuilding, or copies required by an intentionally
retained external source.

Editor construction currently creates identity correspondence vectors for all
entity kinds, even for untracked use. That cost is retained in every measured
variant. The selected design removes unconditional correspondence accumulation;
these measurements do not quantify that change. Nor do they measure rollback or
predict the performance of the proposed recovery implementation.

The experiment's input definitions, method, checks, results, and interpretation
are retained here so scratch sources and logs can be deleted without losing the
design evidence.

### Clone and first-write measurements — 2026-09-21

#### Decision supported

The initial Molecule clone costs 18–22 ns and allocates
nothing on these fixtures. Copying shared tables on first mutation is the
material cost: at 88 atoms, the atom and bond probes individually cost 2.63 µs /
16.5 KiB and 1.99 µs / 7.6 KiB. A subsequent write to an already detached table
allocates nothing. Retain one candidate through the operation rather than
repeatedly sharing and detaching modified tables between phases.

These results alone do not establish acceptable clone-and-apply overhead. They
establish that the existing representation already makes the initial clone
inexpensive, while recovery by copying is not free. The earlier successful-path
projection timings cannot justify discarding the original where the public
operation promises recovery. The selected design uses journaling for
that contract and an owning editor for explicitly destructive execution.
No production API was changed by this experiment.

#### Method and limits

Source: 1fd0c3aa7fbda4e47ec122d0222dd5829ab772a9, with only discussion documents
modified. macOS arm64, rustc 1.96.0 (ac68faa20, 2026-05-25), release profile,
offline dependencies. The standalone Rust probe used the current graph,
graph-IR, and graph-core crates directly. It was run with
`cargo run --release --offline --manifest-path scratch/213-clone-cost/Cargo.toml`.
No production source changes, Python builds, or workspace-wide tests were needed.

Reuse the six SMILES fixtures listed in the preceding experiment, constructed
with ingest_smiles before timing. Each row is the median of 11 batch means,
2,048 operations per batch. Setup and destruction are outside timing:

- Clone: preallocate 2,048 empty Option<Molecule> slots, then fill each with a
  clone of the retained source. The source is passed through black_box.
- Atom/bond first write: prepare 2,048 source clones, then assign NumForm::Lit(1)
  to the first atom/bond's charge through atom_mut/bond_mut. Keeping the source
  alive forces copy-on-write. This is a storage probe, not a chemical operation.
- Detached write: prepare the same clones and obtain each target mutable view
  once before timing, forcing detachment. Time the same charge assignment.
- Aromatic storage: reconstruct an Arc<VarRelationSet<NodeId, AromaticSystemForm>>
  from the molecule's public aromatic-system views, preserving participants and
  attributes. Time Arc::make_mut followed by data_mut on the first relation and
  the same charge assignment. The detached control calls make_mut during setup.
  This uses the storage type held by AromaticSystems; it excludes the checked
  public mutation wrapper's candidate creation and integrity validation.

An untimed operation counts allocation/reallocation calls and requested bytes
with a system-allocator wrapper. Counting is disabled during timing. Requested
bytes describe allocation traffic, not peak live memory or RSS; the retained
source remains alive, and the candidate adds the copied table. Cloning shares
graph and six relation sets as well as the atom and bond tables. Only the
explicitly targeted table detaches in each probe. Independent assertions check
the assigned fields and compare the retained molecule to a separately ingested
original; the relation probe also checks its retained source is unchanged.

All fixtures have **zero global constraints**. Constraints owns its vector and
clones its entries eagerly, so these results do not establish constant-cost
cloning for constraint-bearing molecules. Nested atom/bond forms can likewise
increase table-copy costs beyond these fixtures. The aromatic probe covers one
variable relation-set shape, not every overlay. Graph topology mutation is not
measured. combined8 is the same synthetic disconnected scaling fixture as
before; the observed range stops at 88 atoms. Batch working sets and allocator
reuse affect timings; these are bounded cost probes, not application throughput
predictions or a new journal comparison.

#### Results

First-write times include table detachment and the field assignment. All clone
and detached-write rows requested zero allocations.

| Fixture | Atoms / bonds | Aromatic / stereo atom / stereo bond rows | Clone ns | Atom first µs / bytes | Bond first µs / bytes |
| --- | ---: | ---: | ---: | ---: | ---: |
| octane | 8 / 7 | 0 / 0 / 0 | 22.4 | 0.166 / 1,576 | 0.194 / 656 |
| chain64 | 64 / 63 | 0 / 0 / 0 | 18.7 | 1.948 / 12,328 | 1.375 / 5,584 |
| naphthalene | 10 / 11 | 1 / 0 / 0 | 17.6 | 0.204 / 1,960 | 0.261 / 1,008 |
| stereo | 6 / 5 | 0 / 1 / 1 | 18.1 | 0.141 / 1,192 | 0.139 / 480 |
| combined | 11 / 11 | 1 / 1 / 1 | 18.1 | 0.228 / 2,152 | 0.273 / 1,008 |
| combined8 | 88 / 88 | 8 / 8 / 8 | 19.0 | 2.628 / 16,936 | 1.988 / 7,784 |

Each atom/bond first write made two allocations: the new Arc payload and its
table buffer. Detached atom writes took 3.9–20.6 ns; detached bond writes took
5.4–12.9 ns across these fixtures.

| Fixture | Aromatic storage first µs | Allocation calls | Requested bytes | Detached write ns |
| --- | ---: | ---: | ---: | ---: |
| naphthalene | 0.097 | 7 | 488 | 3.6 |
| combined | 0.099 | 7 | 408 | 3.7 |
| combined8 | 0.278 | 14 | 1,948 | 3.4 |

The fixture definitions, reconstruction recipe, measurement boundaries, results,
checks, and interpretation above remain valid records after deleting scratch.

### Journal cost for field mutation — 2026-09-21

#### Decision and scope

Journaling remains the preferred recovery approach: its entries grow with the
operations and captured payloads rather than requiring copies of whole affected
tables. The motivation for journal-free mutation is its advantage over recording
anything, not a preference for cloning over journaling. Undo coverage is not an
open obstacle. Preparation followed by infallible commit is not an independent
design option without a specific proven transformation demonstrating it.

The measurements below isolate that recording cost and distinguish it from the
current batch implementation's additional bookkeeping. They do not implement
a future transaction redesign or measure a complete transformation.

#### Method

Same source, compiler, allocator probe, fixtures, 11 samples, and 2,048 operations
per sample as the clone experiment. Prepare editors before timing and detach the
target atom/bond table or materialize the aromatic table through its mutable
accessor. Thus no measured path copies a shared table. Each operation changes the
first entity's charge to NumForm::Lit(1):

- **Direct:** assign through the editor's mutable view.
- **Record:** the same assignment plus an isolated recording prototype using the
  actual Undo variant. Clone the old charge, form its inverse field change,
  allocate Vec<Undo>::with_capacity(1), and retain the entry. This isolates
  recording; it is not a new transaction implementation and performs no checks.
- **Apply:** execute the corresponding preconstructed single-edit Edits through
  the current editor apply method, including its preconditions and handle setup.
- **Transact:** execute the same Edits through the current transact method and
  retain its returned journal.

Editor/Edits construction, initial detachment, output destruction, journal
destruction, and final integrity validation are outside timing. Apply and
transact consume the prepared Edits during timing. Untimed checks verify the
changed field, the exact issued Undo entry, full original-molecule restoration
through Transaction::rollback, and the published apply result for every case.
The recording prototype constructs the same inverse entry checked against the
production journal. The aromatic editor uses its current materialized entry
vector, whereas the preceding aromatic copy probe used relation storage.

Scratch binaries are journal and journal_relation in the preceding probe crate,
run with `cargo run --release --offline --manifest-path
scratch/213-clone-cost/Cargo.toml --bin journal` (and --bin journal_relation).
No production code changed. Timing is sensitive to the prepared working set:
the larger editor objects and correspondence buffers make direct-write timings
different from the earlier Molecule probes. Use the paired columns here for
journal overhead; do not divide these by earlier direct-write timings.

#### Results

All times are ns per successful operation, excluding teardown as above.

| Fixture | Field | Direct | Record | Apply | Transact |
| --- | --- | ---: | ---: | ---: | ---: |
| octane | atom | 4 | 51 | 124 | 221 |
| octane | bond | 6 | 91 | 137 | 260 |
| chain64 | atom | 68 | 199 | 276 | 517 |
| chain64 | bond | 15 | 107 | 226 | 443 |
| naphthalene | atom | 4 | 43 | 135 | 233 |
| naphthalene | bond | 6 | 78 | 147 | 277 |
| naphthalene | aromatic | 4 | 41 | 144 | 246 |
| stereo | atom | 5 | 39 | 148 | 272 |
| stereo | bond | 6 | 81 | 148 | 304 |
| combined | atom | 4 | 56 | 167 | 310 |
| combined | bond | 6 | 97 | 170 | 341 |
| combined | aromatic | 4 | 55 | 161 | 331 |
| combined8 | atom | 110 | 248 | 403 | 757 |
| combined8 | bond | 25 | 147 | 328 | 658 |
| combined8 | aromatic | 6 | 103 | 261 | 601 |

Direct writes allocate nothing. Every isolated record makes one 752-byte
allocation: size_of::<Undo>() is 752 bytes on this target because the enum must
accommodate its largest variant. This is the current representation's cost for
one small field entry, not an inherent minimum for journaling. Exact-capacity
reservation matches current transact. These probes do not measure vector growth
over multiple edits or heap-owning field values.

Allocation totals are identical for each field within a fixture:

| Fixture | Apply calls / requested bytes | Transact calls / requested bytes |
| --- | ---: | ---: |
| octane | 2 / 120 | 5 / 992 |
| chain64 | 2 / 1,016 | 5 / 2,784 |
| naphthalene | 3 / 176 | 7 / 1,104 |
| stereo | 4 / 104 | 9 / 960 |
| combined | 5 / 200 | 11 / 1,152 |
| combined8 | 5 / 1,600 | 11 / 3,952 |

The current batch paths are not operation-sized overall: ApplicationState::new
allocates initial handle tables for all entity kinds in both paths, and transact
also clones the session correspondence before execution. At combined8, the
1,600 bytes of handle tables plus another 1,600 bytes for the saved
correspondence plus the 752-byte entry explain the transact total. That
receiver-sized bookkeeping is distinct from the journal and must not be
presented as an inherent journaling cost.

For this single-field workload, isolated recording adds roughly 35–140 ns over
direct writes and a fixed 752 bytes; the existing batch transaction adds more
work. The earlier atom/bond table copies were 1.99–2.63 µs and 7.6–16.5 KiB at
88 atoms. This supports the operation-sized journal approach for sparse changes,
while exposing current bookkeeping and entry-size costs. It does not establish
whole-transformation speedups or require a broader benchmark campaign.

## Staged implementation plan

Graph-core mutation and restoration from 166 and the aggregate integrity gate
optimized in 229 are prerequisites already implemented. Each subitem includes
focused tests of its stated behavior, using public operations for property tests.
Run affected-crate tests and checks as each stage closes; every stage ends green.
Only breaking signature changes and rewires may leave the tree temporarily red
within a stage. Python checks use `umol-py/.venv` with Python 3.13. S0 records
the benchmark baseline before those changes.

### S0 — Baseline and additive types

- **S0a — completed 2026-09-24** (`umol-graph-ir/benches`,
  `umol-graph/benches`, existing external
  tests; additive) Record current direct mutation, batch apply, transaction,
  resolution, and projection costs on sparse and dense fixtures, with unique
  and shared inputs. Count allocations and journal entries/peak retained bytes
  where a journal exists. Add exact public-behavior cases for failure recovery,
  publication, and phase rejection before changing the implementation.
  [dep: none]

  Baseline at `3c54cff68`, Rust 1.96.0 on aarch64-apple-darwin. The retained
  [editor benchmark](../umol-graph-ir/benches/editor.rs) uses an 8-atom chain
  with one atom-field edit and an 80-atom chain with eight edits, eight disjoint
  aromatic systems, and eight dative bonds. Direct means `edit` + field writes
  + `try_build`; apply means `Molecule::apply`; transact means `edit` +
  `transact` + `try_build`. The retained
  [resolve benchmark](../umol-graph/benches/resolve.rs) uses octane and eight
  disconnected copies of `[13CH3][C@H](F)/C=C/c1ccccc1`. Unique inputs are
  independently constructed; shared inputs are cloned from a retained source.
  Setup and edit-list construction are outside timing. All paths include their
  current integrity publication checks. Each time is Criterion's point estimate
  from 10 samples (0.1 s warmup, 0.2 s measurement), in microseconds.

  | Fixture | Input | Direct | Apply | Transact | Resolve | Project |
  | --- | --- | ---: | ---: | ---: | ---: | ---: |
  | chain8 / octane | unique | 0.355 | 0.448 | 0.575 | 12.172 | 1.317 |
  | chain8 / octane | shared | 0.339 | 0.411 | 0.553 | 11.906 | 1.249 |
  | overlays80 / combined8 | unique | 2.595 | 2.938 | 3.430 | 356.410 | 80.134 |
  | overlays80 / combined8 | shared | 2.598 | 2.823 | 3.378 | 368.580 | 74.634 |

  A separate untimed counting-System-allocator probe measured allocation calls
  and requested bytes for one operation after setup. The probe code and text
  output were removed after recording the results; it did not affect the
  retained benchmark. Unique/shared allocation counts agreed for each fixture.

  | Fixture | Direct calls / bytes | Apply calls / bytes | Transact calls / bytes | Resolve calls / bytes | Project calls / bytes | Journal entries / inline bytes |
  | --- | ---: | ---: | ---: | ---: | ---: | ---: |
  | chain8 / octane | 4 / 1,696 | 6 / 1,816 | 9 / 2,688 | 83 / 39,808 | 11 / 11,024 | 1 / 752 |
  | overlays80 / combined8 | 7 / 17,208 | 11 / 18,608 | 16 / 26,024 | 6,053 / 1,304,656 | 620 / 438,396 | 8 / 6,016 |

  The returned journal length is the peak entry count for these successful
  one-batch operations. Inline bytes are `len * size_of::<Undo>()` (752 bytes
  per entry here); these field-only undos have no nested heap payload. The
  figure excludes the 24-byte Vec header and any spare capacity. Thus the
  retained journal footprint requested by the current exact-capacity reservation
  is 776 / 6,040 bytes for these two cases, excluding allocator metadata.
  Allocation bytes are request traffic, not peak live memory or RSS. The current
  `Molecule::edit` shares storage with its input even when that input is unique;
  its first write detaches the atom table. This explains why unique and shared
  allocation counts match, and does not predict the proposed owning editor.
  Exact external cases now pin recovery after a valid edit followed by an
  invalid handle, publication rejection of parallel localized bonds, and
  source preservation after late isotope projection on an ingested
  stereo/aromatic molecule. The existing composite phase-error table also
  covers stereo and aromatic rejection order.

- **S0b — completed 2026-09-24** (`umol-utils::solution`; additive) Add
  `Solution<T, C, U = T>`, keeping the existing two-parameter behavior and
  equal-payload methods. Test both equal and distinct payload types and the
  existing conversion laws. [dep: none]

  `Underdetermined` now carries `U`; the two-parameter form still uses `T`.
  Predicates, contradiction mapping, and conversion methods accept either
  payload shape. `data`, `into_data`, and `map` retain their equal-payload
  signatures. Both shapes and their conversion laws pass `umol-utils` tests;
  the dependent Rust and Python binding crates compile with the default form.

### S1 — Typed-set storage delegation

All six typed sets use `restore_topology_ids` to translate the atom and bond ids
stored in surviving entries back to their original ids. `restore` then
reinstates saved entity entries at their original ids. These methods delegate
to graph-core's `restore_participants` and `restore`, respectively.

- **S1a — completed 2026-09-24** (`ir::aromatic`, `ir::multicenter`; additive) Give the owning typed
  sets the needed add, remove, restore, and atom-list mutation methods, delegating
  incidence and row compaction to graph-core. Test ordered atoms,
  incidence, compaction, and matching-history restoration. [dep: S0a]

  The following shorthand applies once for each `(Set, Id, Form)` pair:
  `(AromaticSystems, AromaticSystemId, AromaticSystemForm)` and
  `(MulticenterBonds, MulticenterBondId, MulticenterBondForm)`. `position` is a
  zero-based offset in the selected atom sequence; S2c gives public editor views
  `AtomPosition` and converts it at this boundary.

  ```text
  Set::add(&mut self, atoms: &[AtomId], attributes: Form) -> Id
  Set::remove(&mut self, ids: &[Id])
  Set::tracked_remove(&mut self, ids: &[Id]) -> Compaction<Id>
  Set::restore(&mut self, rows: &Compaction<Id>, removed: Vec<(Id, Vec<AtomId>, Form)>)
  Set::restore_topology_ids(&mut self, graph: &GraphCompaction)
  Set::replace_atoms(&mut self, id: Id, atoms: &[AtomId])
  Set::replace_atom(&mut self, id: Id, position: usize, atom: AtomId)
  Set::insert_atom(&mut self, id: Id, position: usize, atom: AtomId)
  Set::remove_atom(&mut self, id: Id, position: usize)
  Set::compact(&self, graph: &GraphCompaction) -> Self
  Set::tracked_compact(&self, graph: &GraphCompaction) -> (Self, Compaction<Id>)
  ```

  These methods are `pub(crate)`; S2 exposes public mutation through the editor
  views. They delegate to `VarRelationSet`, retaining the owning set's
  copy-on-write storage. Focused tests cover order, draft duplicates, incidence,
  attribute preservation, compaction, and matching-history restoration.
  The implementation matches the interfaces above. Focused tests, strict
  graph-IR Clippy, nightly formatting, and `git diff --check` pass.
- **S1b** (`ir::dative`, `ir::noncovalent`; additive) Add the same ownership
  operations for distinguished acceptor/donors and fixed endpoints. Test each
  factor independently, including duplicate and temporarily invalid draft
  atom references. [dep: S0a]

  `DativeBonds` takes `(donors, acceptor, attributes)` on addition and restores
  the same fields with their original id. The acceptor is the fixed factor;
  only donors admit insertion and removal.

  ```text
  DativeBonds::add(&mut self, donors: &[AtomId], acceptor: AtomId, attributes: DativeBondForm) -> DativeBondId
  DativeBonds::remove(&mut self, ids: &[DativeBondId])
  DativeBonds::tracked_remove(&mut self, ids: &[DativeBondId]) -> Compaction<DativeBondId>
  DativeBonds::restore(&mut self, rows: &Compaction<DativeBondId>, removed: Vec<(DativeBondId, Vec<AtomId>, AtomId, DativeBondForm)>)
  DativeBonds::restore_topology_ids(&mut self, graph: &GraphCompaction)
  DativeBonds::replace_acceptor(&mut self, id: DativeBondId, acceptor: AtomId)
  DativeBonds::replace_donors(&mut self, id: DativeBondId, donors: &[AtomId])
  DativeBonds::replace_donor(&mut self, id: DativeBondId, position: usize, donor: AtomId)
  DativeBonds::insert_donor(&mut self, id: DativeBondId, position: usize, donor: AtomId)
  DativeBonds::remove_donor(&mut self, id: DativeBondId, position: usize)
  DativeBonds::compact(&self, graph: &GraphCompaction) -> Self
  DativeBonds::tracked_compact(&self, graph: &GraphCompaction) -> (Self, Compaction<DativeBondId>)
  ```

  `NoncovalentBonds` uses a fixed ordered pair; there is no endpoint insertion
  or removal. Its row methods follow the same remove/restore/compact contracts.

  ```text
  NoncovalentBonds::add(&mut self, atoms: [AtomId; 2], attributes: NoncovalentBondForm) -> NoncovalentBondId
  NoncovalentBonds::remove(&mut self, ids: &[NoncovalentBondId])
  NoncovalentBonds::tracked_remove(&mut self, ids: &[NoncovalentBondId]) -> Compaction<NoncovalentBondId>
  NoncovalentBonds::restore(&mut self, rows: &Compaction<NoncovalentBondId>, removed: Vec<(NoncovalentBondId, [AtomId; 2], NoncovalentBondForm)>)
  NoncovalentBonds::restore_topology_ids(&mut self, graph: &GraphCompaction)
  NoncovalentBonds::replace_atoms(&mut self, id: NoncovalentBondId, atoms: [AtomId; 2])
  NoncovalentBonds::replace_atom(&mut self, id: NoncovalentBondId, position: usize, atom: AtomId)
  NoncovalentBonds::compact(&self, graph: &GraphCompaction) -> Self
  NoncovalentBonds::tracked_compact(&self, graph: &GraphCompaction) -> (Self, Compaction<NoncovalentBondId>)
  ```

  These are proposed crate-private set methods. S2c exposes the domain-level
  editor-view methods with `AtomPosition` after S1d replaces the storage wrappers.
- **S1c** (`ir::stereo`; additive) Add site and ligand operations to the stereo
  typed sets, preserving stored frames and payloads during ordinary replacement.
  Test atom and bond sites, virtual ligands, incidence, and restoration.
  [dep: S0a]

  The shorthand applies to `(Set, Id, Site, Form)` equal to
  `(StereoAtoms, StereoAtomId, AtomId, StereoAtomForm)` or
  `(StereoBonds, StereoBondId, BondId, StereoBondForm)`.

  ```text
  Set::add(&mut self, site: Site, ligands: &[StereoLigand], attributes: Form) -> Id
  Set::remove(&mut self, ids: &[Id])
  Set::tracked_remove(&mut self, ids: &[Id]) -> Compaction<Id>
  Set::restore(&mut self, rows: &Compaction<Id>, removed: Vec<(Id, Site, Vec<StereoLigand>, Form)>)
  Set::restore_topology_ids(&mut self, graph: &GraphCompaction)
  Set::replace_site(&mut self, id: Id, site: Site)
  Set::replace_ligands(&mut self, id: Id, ligands: &[StereoLigand])
  Set::replace_ligand(&mut self, id: Id, position: usize, ligand: StereoLigand)
  Set::insert_ligand(&mut self, id: Id, position: usize, ligand: StereoLigand)
  Set::remove_ligand(&mut self, id: Id, position: usize)
  Set::compact(&self, graph: &GraphCompaction) -> Self
  Set::tracked_compact(&self, graph: &GraphCompaction) -> (Self, Compaction<Id>)
  ```

  StereoBonds restoration translates its bond sites and the atom references in
  its ligands. These are proposed crate-private set methods. S2c exposes
  `StereoLigandPosition` on the public editor views. Site and ligand replacement
  leave the other factor and payload unchanged; frame-preserving permutation
  remains a separate operation.
- **S1d** (`ir::molecule::editor`; internal rewire, red→green) Replace the
  `*SetStorage` wrappers with the typed sets for reads, additions, removals,
  compaction, and restoration. Keep the current public editor lifecycle and
  undo checker for this stage while routing row restoration through typed sets.
  Remove the temporary dead-code expectations on the new typed-set methods.
  Test all six overlay kinds against the S0 behavior and rerun the
  storage-sensitive benchmarks.
  [dep: S1a, S1b, S1c]

### S2 — Editor entity views

- **S2a** (`ir::view`, `ir::molecule::editor`; breaking, red→green) Make each
  mutable overlay editor view borrow its owning typed set and id. Replace public
  fields on immutable and mutable editor views with short-lived accessors,
  including `attributes_mut`; migrate their Rust and Python callers in this
  subitem. Test immediate writes and the specified commutation of participant
  replacement with attribute assignment. [dep: S1d]
- **S2b** (`ir::view`; additive) Add the local immutable and mutable editor-view
  getters listed above, including `ligand(position)` and exact-size participant
  iterators. Test values, ordering, and out-of-range behavior against immutable
  Molecule views. [dep: S2a]
- **S2c** (`ir::view`, `ir::id`, typed sets; additive) Add `AtomPosition` and the
  factor-specific replacement, insertion, and removal methods on mutable views;
  clarify the current-frame meaning of `StereoLigandPosition`. Test position
  boundaries, unaffected factors and attributes, immediate incidence
  maintenance, and constructor-equivalent integrity at the current editor
  publication gate, not at each draft write.
  [dep: S2a, S2b]

### S3 — Batch mutation and reaction vocabulary

- **S3a** (`ir::edit`; breaking, red→green) Add the nine whole-component Edit
  variants and matching saved-value Undo variants. Extend Edits construction
  without independent-batch composition. Test construction and handle
  namespaces; execution comparisons belong to S3b. [dep: S1d]
- **S3b** (`ir::molecule::transact`; breaking, red→green) Realize those edits
  through typed-set mutation and their undos through saved components. Keep
  unrelated factors, attributes, constraints, and ids unchanged. Test each
  forward/undo pair, old-state errors before mutation, and frame alignment.
  [dep: S2c, S3a]
- **S3c** (`ir::delta`; breaking, red→green) Add the nine Delta variants and
  extend frame transport, inversion, normalization, composition, and Add/Remove
  folding under the settled exact component comparison. Test continuity,
  contradiction, identity, and created-entity cancellation laws.
  [dep: S0a]
- **S3d** (`ir::reaction`, `ir::reaction_span`; breaking, red→green) Lower
  replacement deltas through before/after materialization and existing
  correspondence induction and superimposition. Test compatible preserved
  entities, incompatible removal/addition, application, and roundtrips without
  asserting that Delta spelling survives span conversion. [dep: S3b, S3c]
- **S3e** (`umol-py::edit`, delta/reaction bindings; breaking, red→green)
  Extend the existing Python variant families for the new changes, with
  Rust-equivalent construction and failure behavior; migrate exhaustive matches
  and add parity cases. Do not add unrelated Rust API coverage merely for parity.
  This closes the stage after the enum changes. [dep: S3a, S3c, S3d]

### S4 — Recovery machinery before the public lifecycle switch

- **S4a** (`ir::molecule::editor`, `ir::molecule::transact`; additive internal
  work) Route topology and relation restoration through the existing graph-core
  operations and add local attribute, constraint, and target-access guards.
  Keep the existing detached journal surface until S5. Test matching-history
  recovery under `normalized_eq` and panic freedom for manipulated undo data
  without asserting its result.
  [dep: S1d, S3b]
- **S4b** (`ir::molecule::transact`; internal rewire, red→green) Separate shared
  graph-IR mutation kernels from handle realization; prepare fallible data before
  writes, record progress for multi-row edits and cascades, and avoid initial
  receiver-sized handle tables until compaction needs them. Test failures and
  injected unwinds at each completed storage step. [dep: S4a]
- **S4c** (`ir::molecule::transact`; additive internal work) Build the private
  scoped recovery guard and journal on those kernels. Verify callback error,
  cancellation, `mem::forget` of the supplied handle, caught apply panic,
  commit failure, and restoration before the receiver is observable. Record
  journal count and peak retained bytes for the dense resolver fixture.
  [dep: S4b]

### S5 — Borrowed transaction API

- **S5a** (`ir::molecule::transact`, `ir::error`; breaking, red→green) Replace
  the detached Transaction with `Transaction::run`, immediate `apply`, checked
  `probe`, `commit`, `tracked_commit`, and `rollback`; add Molecule's prepared-batch
  `transact` and `tracked_transact`. Retire detached-journal `undos`, `append`,
  `rollback`, and `tracked_rollback`, plus editor `transact`/`tracked_transact`,
  without adding Transaction::tracked_apply. Remove `validate_undo` and obsolete
  rollback errors after S4's local guards. Retain the drafted MoleculeApplyError
  return type. Test separate batch handle namespaces, whole-transaction
  correspondence, abort latching, callback value return after no-commit
  rollback, and checked acceptance. [dep: S4c]
- **S5b** (`umol-graph-ir` reaction and molecule callers; breaking, red→green)
  Migrate uses of detached journals. Preserve reaction `Ok(None)` versus error
  classification and host-to-product correspondence.
  Test failed applications and product integrity. [dep: S5a]
- **S5c** (`umol-graph::ops`; breaking, red→green) Migrate existing borrowed
  resolve/project execution to scoped transactions, sharing the planned
  batches and using checked probes between phases. Keep their present public
  signatures in this stage. Test late rejection restores the entry molecule and
  preserves chemistry diagnostics. [dep: S5a]
- **S5d** (`umol-py::transaction`, `umol-py::molecule`; breaking, red→green)
  Replace Python's detached journal with prepared-batch `transact` and
  `tracked_transact`; give Edits the consumed-input behavior needed to transfer
  all prepared batches before mutation and expose no Python transaction handle.
  Test independent New namespaces, rollback on late failure, and tracked result
  parity. [dep: S5a, S5b]

### S6 — Owning editor and publication

- **S6a** (`ir::molecule`, `ir::molecule::editor`; breaking, red→green) Make
  `edit(self)` own its Molecule, make Molecule `apply`/`tracked_apply` consume,
  and replace editor snapshot/build methods with checked `probe`/`finish`.
  Retain MoleculeBuilder's asserted build and the existing integrity-preserving
  Molecule mutable methods. Test direct/batch interleaving, invalid probe then
  repair, destructive failure, and constructor-equivalent publication.
  [dep: S1d, S5a]
- **S6b** (`umol-graph-ir` callers, including reaction application; breaking,
  red→green) Migrate editor construction/publication and remove session
  correspondence accumulation and public tracked direct removal. Use one
  source-preserving product candidate where the host survives; keep batch-only
  tracking distinct from whole-transaction tracking. Test product and
  correspondence laws. [dep: S6a]
- **S6c** (`umol-graph`, `umol-io` callers; breaking, red→green) Migrate remaining
  editor publication sites to `finish` or borrowed transactions without changing
  their chemistry or boundary outcomes. Test published outputs and rejection
  behavior. [dep: S6a, S5c]
- **S6d** (`umol-py::molecule`, `umol-py::transaction`, `umol-py::edit`; breaking,
  red→green) Transfer Molecule inputs on consuming calls and reuse S5d's Edits
  transfer; raise `ConsumedError` for the owner and `InvalidatedViewError` for
  owner-backed accessors. Replace Python snapshot/build with finish, without a blanket
  consumed-state migration of unrelated forms. Test aliases, nested views,
  successful and failed consumption, and explicit copies. [dep: S6a, S6b]

### S7 — Resolution, projection, and boundaries

- **S7a** (`umol-graph::ops::resolve` and phase modules; breaking, red→green)
  Extract shared planning, rename phase `plan` to `plan_resolve`, add the three
  `plan_project` methods, and implement consuming `resolve`/`project` alongside
  recovering `resolve_into`/`project_into`. Preserve phase order, rejection
  priority, and standalone-versus-composite isotope policy. Test accepted
  outputs, each later-phase rejection, both ownership contracts, and application
  versus integrity errors through MoleculeApplyError in phase-specific errors.
  [dep: S0b, S5c, S6a]
- **S7b** (`umol-graph::ops::resolve`; breaking, red→green) Make reports opt-in
  with `resolve_with_report` and `resolve_into_with_report`; avoid default
  report-only collection while retaining resolution candidates. Test equal
  outcomes and final molecules with and without reports, including
  underdetermination. [dep: S7a]
- **S7c** (`umol-graph::ingest`, `parse`, `export`, `umol-io` boundaries; breaking,
  red→green) Pass owned ingest/parse candidates to report-free resolution; keep
  one intentional source copy for export projection. Remove report payloads from
  default underdetermination errors. Test boundary diagnostics and output
  integrity. [dep: S7a, S7b]
- **S7d** (`umol-py::resolve` and boundary adapters; breaking, red→green)
  Mirror the consuming/borrowed names, explicit reporting, and payload-free
  ingestion underdetermination. Test Python ownership, result shapes, and error
  parity. [dep: S6d, S7a, S7b, S7c]

### S8 — Transformations and remaining operations

- **S8a** (`umol-graph::ops::transform`; breaking, red→green) Change Transformer
  to consuming `transform`, recovering `transform_into`, and lazy
  `transform_iter` returning `impl Iterator`; migrate all trait callers. Test
  failure ownership, independent iterator outputs, and deferred execution.
  [dep: S5c, S6a]
- **S8b** (`umol-graph::ops::transform::{aromatizer,delocalize_charge,kekulizer}`;
  breaking, red→green) Share each operation's planning between the two routes;
  use one publication gate and retain existing chemistry errors, including
  DelocalizeCharge's `Infallible` contract. Test each operation's success,
  rejection, and checked output. [dep: S8a]
- **S8c** (`ir::molecule`; additive) Make `combine_from` recover its append/count
  checkpoints on unwind without an initial whole-molecule copy. Test unchanged
  receiver after a caught unwind. The separate lift_constraints correction
  is not part of this subitem. [dep: S6a]

### S9 — Contract and performance closeout

- **S9a** (`docs/development`, public rustdoc, examples; additive) Reconcile
  the living data-type, nomenclature, integrity, and Python API guides with the
  implemented lifecycle; remove stale snapshot/detached-journal descriptions.
  Document precise failure, ownership, and integrity boundaries without citing
  discussion records from source. Check links and examples. [dep: S5a, S6d, S7d, S8b]
- **S9b** (workspace gates and benchmarks; additive) Compare the final direct,
  apply, and transaction paths with S0, including dense resolver journal size
  and the source-preserving export path. Review the complete diff against 213,
  run formatting, workspace tests, strict Clippy and rustdoc, explicit feature
  suites, Python 3.13 build/tests, and the pinned Rust 1.87 gate once at final
  closeout. Record results and update the discussion status only after the full
  scope passes. [dep: S8c, S9a]

**Critical path:** S0 → S1 → S2/S3 → S4 → S5 → S6 → S7/S8 → S9.
S3c can proceed alongside S2; S3b, S3d, and S3e wait for S2c as recorded in
their dependencies. No speculative optimization stage is required.
Immutable-view simplification, graph-core bond
endpoint rewiring, intermediate transaction tracking, an interactive Python
transaction, Undo compression, and the hydrogen operations in 166 remain outside
this plan; the 213 mutation surface enables the latter work.
