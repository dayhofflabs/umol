# 213 — Molecule and reaction mutation

Status: In Progress
Date: 2026-08-27
Relates: [117](117-entity-model-extensibility-2026-06-20.md),
[166](166-molecule-ops-2026-07-27.md),
[211](211-relation-frames-and-api-2026-08-26.md),
[214](214-aggregate-frame-semantics-2026-08-28.md),
[228](228-python-api-parity-2026-09-21.md),
[229](229-aggregate-integrity-review-2026-09-22.md),
[231](231-view-arguments-2026-09-25.md),
[232](232-view-accessors-2026-09-26.md),
[data-type guide](../docs/development/data-types.md),
[nomenclature guide](../docs/development/nomenclature.md)

## Design status — 2026-09-29

This document owns the molecule/reaction mutation redesign. S0a–S0b, S1a–S1c,
S2a, the revised S2b, S2c, S2d, S2g, S2h, and S2i1 are implemented. The previous S2b mutable-view
attempt was reverted. S2i2–S2i4 are complete; the migration compiles and its
verification passes. S2i5 is complete: mutable molecule/editor views are separate
types. All eight entity-view families expose private ids and attribute borrows
through matching accessors. Stereo views use the owning sets for site and ligand
access; their ligand frames are borrowed. S2j is complete for all eight entity
families: matching local getters and editor-only structural mutation are implemented.
S2k1's in-place DSL conversion, S2k2's Rust callback caller migration, and S2l's
Python assignment and read-only molecule constraint access are implemented.
S2m removes the remaining mutation callbacks and closes S2; S2f is cancelled.
S3a1–S3b are implemented. Replacement Deltas are withdrawn; the nine replacement
Edits and their Undo variants remain. S3c/S3d record the selective removal and
retained reaction integration; both are verified. S3e completes the Python Edit
migration. S3f's graph-core bulk additions, S3g's relation-set bulk additions,
and S3h's typed-overlay extend methods are implemented. S3i's Molecule/editor bulk
additions and S3j's mutable correspondence methods are implemented. S3k's
index-arithmetic cleanup is complete across graph-core, graph-ir, and graph.
S4a is complete: undo restoration calls Molecule and constraint storage methods.
Editor batch loops remain under the editor module; single-edit execution and
handle state are in molecule::apply. Fields remain private and internal Molecule
mutation methods use pub(crate). S4b is complete: S4b9 consolidates the entity-set
implementations and names their mapping-returning compaction methods compact.
S4d is complete: comparable single-entity entries live in their owning entity
modules, and both Edit execution paths use specialized framed_eq implementations.
S5a is implemented: scoped transactions provide immediate batch application,
checked probes and commit, rollback, and optional correspondence. Molecule's
prepared-batch conveniences use the same lifecycle. S5b's graph-ir caller and
test migration is complete; graph-ir passes its checks. S5c's borrowed graph
operations use scoped transactions and pass the graph checks. S5d1 implements
Python Edits consumption. S5d2 implements accessor counters and Storage names.
S5d3 implements prepared-batch Python transactions; S5 is complete.
S6a implements consuming edit/application and checked probe/finish publication.
S6b1 migrates graph-ir callers and removes editor session correspondence. S6b2
rewrites combine_from to append through Molecule methods. S6c1 implements the
transformation plans and borrowed execution. S6c2 completes the remaining Rust
caller migration. S6d1–S6d3 complete Python Molecule consumption, fallible owner
access, and editor finish; the S6 workspace gate passes. S7a implements consuming
resolve/project and recovering resolve_into/project_into, with shared phase plans.
S7b is next: separate opt-in reporting from the composite resolver's borrowed path.
Graph-core mutation and restoration are complete in
[166](166-molecule-ops-2026-07-27.md); editor integration remains here. After
that integration, return to 166 for the operation changes and hydrogen folding.
Doc 228 is unchanged by this review and its withdrawn ownership migration is not
an implementation dependency.

[232](232-view-accessors-2026-09-26.md) separately reviews id names, collection
return shapes, accessor allocations, and graph algorithms in views. It does not
reopen S2i or block S2j.

| Area | Status | Concrete position |
| --- | --- | --- |
| Storage delegation, participant methods, Edit/Undo variants, local getters | Settled design; S1a–S1c complete | Use the existing typed entity sets and graph-core mutation/restoration; contracts below. |
| Editing and recovery | Settled design | Owning, destructive editor; separate borrowed, scoped transaction. Editor and Transaction probe check integrity and return an immutable Molecule borrow; no probe callback. |
| resolve/project/transform consumers | Settled design; integration work remains | resolve/project consume destructively; resolve_into/project_into mutate borrowed inputs with recovery. Consuming resolution uses Solution<Molecule, C, ()>; reporting is explicit. Ingest uses report-free resolution. Transformer signatures follow the same ownership naming. |
| Molecule attribute methods | Implemented; S2 complete | Mutable borrows expose every entity attribute and entity-level constraint in Molecule and MoleculeEditor. Rust and Python retain simple assignment, including aromatic/multicenter/stereo. Mutation callbacks are removed. |
| Entity-view structures and API | S2i5 and S2j complete | Molecule uses *View / *ViewMut; editor uses *EditorView / *EditorViewMut. Corresponding molecule/editor methods have identical signatures and semantics. All attributes remain freely mutable; structural mutation is editor-only. |
| Molecule-level constraint mutation | Implemented; S2 complete | Molecule::constraints provides reads; the editor exposes &mut Constraints. Python molecule constraint entries and iteration are lazy and read-only. try_modify_constraints is removed. |
| Transaction correspondence | Settled design | tracked_commit returns the whole transaction's correspondence. Omit Transaction::tracked_apply unless a concrete need for intermediate tracking arises. |
| Python bindings | S5d and S6d complete | Molecule.transact and tracked_transact submit prepared Edits under a borrowed transaction. Molecule edit/apply/tracked_apply consume; editor finish publishes. Owners and views enforce consumption/invalidation. Resolver ownership and reporting changes remain in S7d. |
| Edits and multiple batches | Settled | Edits accumulates one sequence. Multiple batches execute through separate Transaction::apply calls under one commit/rollback boundary. No independent-batch composition API on Edits. |
| Mutation errors | Settled design | Retain application/integrity categories and chemistry outcomes; add Aborted and remove obsolete rollback failures. ResolveError::Apply and ProjectError::Apply carry MoleculeApplyError. |

The staged implementation plan below sequences these contracts and integration
obligations. S0a–S0b record the baseline and additive Solution type; S1a–S1c add
typed-set mutation for all six overlay kinds. The interrupted S1d changes were
reverted. The S0c attempt is also reverted: removing Molecule mutation removed
useful editing APIs without an adequate replacement, particularly in Python.
The revised access model retains Molecule mutable views. S2 sequences the
replacement APIs, caller migration, and callback removal. The previous S2b
implementation, field interfaces, added errors, and checked entity-constraint
API were reverted.
The [bounded study](#fieldframe-agreement-first-use-study--2026-09-24) records
the consumer changes enabling unchecked attribute assignment. The revised design
keeps molecule/editor views separate for both immutable and mutable access.
S2i5 removed the mutable-view const parameter without changing typed editor storage
or attribute assignment. S2j implements the local-getter inventory and editor-only
structural mutations below for all eight entity families.
The lift_constraints defect and its undetermined-stereo policy are a separate
focused correction, recorded under
[other operations](#other-moleculereaction-operations).

Replacement verbs and payloads in the Edit DSL are approved in S3.
Python consumption, counter-based accessor invalidation, and storage names are
approved below. S2b is complete: Rust's unit error is NoJoinError and Python
join raises NoJoinError. S2c's bounded coset-operation fixes and S2d's role-only
incidence/count-aware consumers are complete. S2f is cancelled. S2g's frame-consumer
checks, S2h's aggregate-integrity changes, and S2i1–S2i5 are complete. S2j is
complete. S2k1, S2k2, S2l, and S2m are implemented; S2 is complete. S3a1–S3a3
and S3b are implemented. S3c/S3d remove replacement Deltas while retaining the
approved reaction names, semantics, and dative-factor migration. S3e–S3k, S4a,
and S4b are complete. S4d's caller migration and comparison optimization are
complete; S4b9 closes the strict lint gate. S5a–S5d3 are complete, including
Python runtime verification. S6a–S6b2 implement the owning editor, migrate
graph-ir callers, and rewire combine_from. S6c1–S6c2 complete the remaining Rust
caller migration. S6d1–S6d3 complete Python ownership and the S6 gate; S7a is next.

## Editor and transaction API

```rust
pub struct MoleculeEditor { molecule: Molecule }
pub struct Transaction<'scope> {
    molecule: &'scope mut Molecule,
    journal: &'scope mut Vec<Undo>,
    status: &'scope mut TransactionStatus,
}

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

The transaction surface and owning editor interfaces above are implemented.
Editor publication uses probe/finish, and Molecule edit/application consume their
input. S6b1 removes editor session correspondence and migrates graph-ir callers;
S6c–S6d migrate chemistry, format, and Python callers.
MoleculeBuilder retains its asserted build for fresh construction.

Transaction::run is the scoped execution entry point, not a constructor returning
a transaction. The approved ownership arrangement is a private guard local to
run, owning the mutable molecule borrow, journal, and transaction status. The
Transaction passed to the callback borrows access to those fields; it does not
own the guard. The guard remains alive until the callback returns. Its Drop
restores unaccepted changes, so forgetting the handle cannot disable cleanup.
Successful commit requests acceptance; run accepts only on callback Ok as well.
Molecule::transact and
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
Transaction cannot forget the library's restoration guard. Callback unwinding
between completed edits restores their journal. Recovery from allocation
failure or internal panics during an edit is not part of the contract.

There is one owning editor, without a mode or borrowed alternative. Transaction
borrows Molecule directly and calls the same crate-private mutation methods as the
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

The accepted Transformer interface follows the same ownership naming. Rename generate_all to
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
None needs recovery of a rejected local candidate. combine_from must migrate
from its move-out/editor route without an initial recovery clone; append/count
checkpoints for internal panics are not required. inline_constraints already plans meets before
writes and needs no journal merely for centralization. lift_constraints has a
reproduced drain-before-panic defect; move its admission before writes and settle
its undetermined-stereo policy as a focused operation correction.

Normalize, reframe, remap, and canonicalize retain their established contracts.
Reaction/ReactionSpan definition editing continues through existing parts/entries
and checked construction; no current consumer requires a new ReactionEditor.
Reaction product integrity conflicts retain their existing Ok(None) versus error
classification, and product correspondence keeps its own pairing semantics.

## Additional API contracts

### Molecule mutation boundary

Retain ordinary entity-field and entity-level constraint mutation uniformly
across all eight entity kinds, in both Rust and Python. Molecule and MoleculeEditor
both expose all entity attributes through mutable borrows of their forms.
Simple assignment applies to aromatic systems, multicenter bonds, stereo atoms,
and stereo bonds just as it does to atoms and localized, dative, and noncovalent
bonds. It includes electron counts, complete stereo configurations, and direct
entity-level constraint mutation. No assignment-time integrity check, special
setter, callback, per-field mutable-accessor family, or attribute permission
is needed. Python property and nested-constraint writes change the stored value
through the corresponding Rust mutable borrow, without a copy/edit/publication
cycle or detached mutable copy.

Remove all Molecule modify_* and try_modify_* callbacks after migrating their
callers; these methods are still present in the restored source. This includes
private bulk callbacks, try_modify_checked, and try_modify_constraints.
Molecule-level constraints are read-only on Molecule and mutable through the
editor. Remove public Molecule::constraints_mut and ConstraintsViewMut implemented in
S2a; retain the editor’s &mut Constraints access. Entity-level constraints remain
directly mutable through entity attributes.

The previous S2b implementation was reverted in full. Its per-field mutable accessors and checked
entity-constraint methods are withdrawn. The study below establishes the required
consumer changes for removing electron-count length and coset-index agreement
from aggregate integrity checks. Implement those changes with the constructor
changes across molecule, reaction, and reaction-span paths before exposing the
uniform mutation surface. No integrity checks have yet been removed in code.

Charge delocalization and MoleculeDsl conversion mutate entity attributes through
Molecule's mutable views. Neither requires editor publication or an integrity
check. Python preserves field assignment and direct entity-constraint mutation;
molecule-level constraint writes use the editor.
Remove the molecule-backed top-level constraint mutation paths and their binding
callback plumbing.

### Field/frame agreement: first-use study — 2026-09-24

**Finding:** moving these checks to the operations that interpret the fields is
feasible. The checks are small and concentrated. It requires explicit failure
handling in count-dependent canonicalization/matching, stereo normalization and
transport;
deleting the constructor checks alone would leave panics and incorrect results.
No integrity check or replacement mutation API has been implemented by this study.

Scope is electron-count/atom-list agreement, stereo configuration/frame agreement,
and the frame-relative entity constraints affected by direct constraint mutation.
Topology references, incidence, uniqueness, and the storage bound on stereo-frame
size are unchanged. Raw field reads, assignments, cloning, and faithful DSL
serialization need no new checks.

#### Electron counts

For a literal vector the check is counts.len() == atom_count. Undetermined has
no vector to check. The direct consumers are:

| Location | Present behavior on disagreement | First-use consequence |
| --- | --- | --- |
| Molecule::incidence_graph and ReactionSpan::incidence_graph / electron_span in [incidence.rs](../umol-graph-ir/src/ir/incidence.rs) | Short vectors panic while copying counts[position] into edge labels; long vectors lose their extra entries in those labels. Connectivity does not require electron counts. | Remove electron-count payloads from incidence labels and their construction. Keep both incidence builders infallible; check agreement in consumers that interpret counts. |
| aromatic_valence and multicenter_valence in [view/atom.rs](../umol-graph-ir/src/ir/view/atom.rs) | Checked access makes a missing cell Undetermined; existing cells still produce literal contributions and excess cells are ignored. Length disagreement causes no indexing panic. | Accepted as-is. Keep NumForm return types and existing behavior; add no row-length check to these getters. |
| AromaticSystemView::electron_count and MulticenterBondView::electron_count in [view/aromatic.rs](../umol-graph-ir/src/ir/view/aromatic.rs) / [view/multicenter.rs](../umol-graph-ir/src/ir/view/multicenter.rs) | Sum every stored entry, including entries with no member atom. A short vector contributes only its stored values. | Accepted as-is. Keep NumForm return types and existing behavior; add no atom-count agreement check to these sums. |
| initial_color_keys, reaction_span_entity_keys, and constitution_candidate::electron_occurrence_fields in [canonicalize.rs](../umol-graph-ir/src/ir/canonicalize.rs) | Coloring reads count-bearing incidence labels; zip truncates a mismatched atom/count pair when building a canonical key. | Read contributions from entity attributes when building count-dependent colors/keys. Check agreement before that interpretation; no repeated length check inside each candidate is needed. |
| ElectronCountsForm::reframe_by and ::permute in [electrons.rs](../umol-graph-ir/src/ir/electrons.rs), with aromatic/multicenter form and set delegation | reframe_by already returns None on a length mismatch. The older permute method silently leaves the counts unchanged. | The active set reframe paths already propagate failure. Retain that behavior; the public permute method needs an explicit unsuccessful outcome if it is to promise frame transport. Its form wrappers have the same issue. |
| MatchingInput::from_system, DelocalizationPlan::derive, AromaticityPerceiver::derive in [kekulizer.rs](../umol-graph/src/ops/transform/kekulizer.rs), [delocalize_charge.rs](../umol-graph/src/ops/transform/delocalize_charge.rs), [aromaticity.rs](../umol-graph/src/ops/aromaticity.rs) | Already check lengths: respectively ElectronCountMismatch, no applicable plan, and AromaticSystemFailure. | Keep these existing local checks. No new validation layer is needed for these operations. |

Required incidence changes:

- Remove electron contributions from molecule and reaction-span incidence edge
  labels, including the count-specific span payloads and electron_span extraction.
  Use atom-role names AromaticAtom and MulticenterAtom instead of
  AromaticParticipant and MulticenterParticipant. Participant is graph-core
  terminology, not the role name for these graph-IR atoms.
- Preserve count-dependent canonicalization by reading the entity attributes at
  coloring/key construction, with agreement checked there.
- Incidence-based substructure matching currently uses edge counts for early
  filtering, then verify_overlays compares the transported entity attributes.
  Remove its dependence on count-bearing incidence labels. Any retained
  count-based filtering must read entity attributes; agreement is required when
  interpreting those counts. Preserve valid-input matching results for both
  matching algorithms. Removing the early filter may affect search cost.
- Graph symmetry does not read these electron-count edge labels; removing them
  creates no electron-length failure boundary for graph_symmetry.
- Update the incidence, canonicalization, and matching tests to reflect role-only
  incidence labels while preserving the operations' valid-input semantics.

The four derived-value getters remain unchanged. For atoms [a, b, c], counts
[1, 2] yield per-atom contributions 1, 2, Undetermined and total 3; counts
[1, 2, 0, 5] yield contributions 1, 2, 0 and total 8. These results on malformed
input are acceptable. Check length only in operations whose work requires
complete atom/count correspondence, not merely because a derived result could
be chemically meaningless on invalid input.

Reframing, gluing, and reaction removal transport already use the checked
FrameTransport path. Stored-vector equality, update, and serialization do not
interpret a cell as an atom contribution and need no length check.

#### Stereo cosets and frames

CosetSpace already checks indices when retrieving or transforming arrangements.
StereoKind::count gives the number of cosets; tetrahedral indices are 0 and 1.
Do not add blanket range validation to normalization, lattice comparisons, raw
access, assignment, or serialization. Meaningful algebraic results are not
required for malformed indices. Keep Rust lattice behavior unchanged; S2b only
renames the unit error from NoJoin to NoJoinError, including the derive macro's
return type. The proposed JoinError enum migration is withdrawn.
Python join raises NoJoinError for Rust's NoJoinError instead of returning None (S2b).
Python meet retains None for bottom; normalization behavior is unchanged.

The coset-operation changes are limited to action compatibility, failure
propagation through existing Option returns, and correct domain simplification.
Reject an incompatible supplied permutation in StereoCoset::apply, including
symbolic and undetermined branches. compose_term returns None rather than
composing incompatible actions; canon_coset then returns the original term
unchanged. Do not turn that case into a new normalization failure.

Plain literals and sets retain their existing normalization behavior. A variable
domain may become unrestricted only when it actually covers the valid range;
otherwise retain the supplied domain, without rejecting or trimming it. Empty
domains remain Contradiction. Existing failure when evaluating a literal action
through checked reindex remains unchanged. coset_apply_permutation propagates
that evaluation failure as None instead of Some(Undetermined).

The concrete gaps and required changes are:

| Location | Present behavior / dependency | First-use consequence |
| --- | --- | --- |
| CosetSpace::unindex, ::reindex, ::enantiomer, ::observable_coset in [coset.rs](../umol-perm/src/coset.rs) | Already use checked lookup and return None for an invalid index; reindex also checks the allowed action. | Reuse these operations. No change to the permutation library is required for index bounds. |
| canon_coset in [stereo.rs](../umol-graph-ir/src/ir/stereo.rs) | Variable domains are erased solely because their size equals kind.count(); tetrahedral {0, 9} becomes unrestricted. | Erase only a domain covering the valid range. Retain other domains. If compose_term cannot compose an action, return the original term unchanged. No new range rejection. |
| compose_term in [stereo.rs](../umol-graph-ir/src/ir/stereo.rs) | A wrong-degree Apply panics in Permutation::compose. | Return None for an incompatible explicit action before composition; canon_coset preserves the input term. |
| StereoCoset::apply, StereoConfigurationForm transport, and StereoAtomView::coset_for / StereoBondView::coset_for | Literal index failures already return None; symbolic and undetermined apply branches bypass action checks. | Check the supplied action against the known kind in StereoCoset::apply. Preserve existing checked literal transport and failure propagation; symbolic construction remains lazy. |
| coset_apply_permutation and its caller symmetry::reexpress | The term branch normalizes, then replaces an error with Undetermined; reexpress subsequently treats it as no literal result. | Return None on normalization failure through the existing Option path. This does not require a new failure return from the public symmetry methods. |
| stereo_refinement_descriptor and canonical_kinded_stereo_frame in [canonicalize.rs](../umol-graph-ir/src/ir/canonicalize.rs) | Permutation::act asserts that ligand length equals the declared kind's degree. | Check that agreement before applying kind-derived permutations. These functions already return Result. |
| Molecule::graph_symmetry, observable_descriptor, has_oriented_center, and project_stereo in [symmetry.rs](../umol-graph-ir/src/ir/symmetry.rs) | An invalid coset can affect the resulting chirality classification. Separately, a kind/frame mismatch can panic while generating ligand swaps, and local projection assumes a determined kind. | Invalid coset indices do not justify making graph_symmetry or stereo_atom_symmetry/stereo_bond_symmetry fallible. No correct chirality classification is required for an invalid index. Keep the separate panic hazards visible without treating them as approval for a fallible symmetry API. |
| StereoConformanceValidator::validate_symmetry and StereoSymmetry::topicity | An out-of-range first topicity position panics; an out-of-range second position can give a false classification. A wrong-degree ligand-symmetry permutation is treated as non-membership, which can satisfy a negative assertion. | Check pair positions and permutation degree before interpreting the entity constraint. Constraint FrameTransport already checks them and returns None. Fluxionality is transported, but this validator does not currently evaluate it. |
| StereoResolver::project; StereoPerception::derive; export::lower_stereo_atom / lower_stereo_bond; Reaction::prepare_deltas | Projection reports failed transport; perception records failed entity candidates; export rejects unsupported literal indices. Reaction application explicitly checks stereo domains. | Retain existing consumer behavior; do not introduce a shared normalization range check or additional validation pass. |

**Failure interfaces:** all public signatures remain unchanged. compose_term is
the only private signature change: it returns Option<(&StereoTerm, Permutation)>.
canon_coset retains its existing Result and preserves the original term on a
composition failure. apply, swap, mirror, reframe_by, and coset_for retain their
Option boundaries. No new lattice or symmetry errors are introduced.
Depiction/CoordGen error-propagation changes (S2f) are cancelled.

Opening the complete configuration field also permits changing kind. Accordingly,
kind/site admissibility and kind/ligand-count agreement require separate review
where consumers rely on them, especially before applying permutations to ligands.
Coset-range handling alone does not resolve these frame-related panic hazards.
Likewise, direct entity-constraint mutation requires the local checks in the
topicity/symmetry consumers above; it does not require checked collection setters.

Reaction constructors repeat count-length checks for Add/Remove and both sides
of ModifyField, and use the shared stereo integrity functions for their payloads
and constraints ([reaction/integrity.rs](../umol-graph-ir/src/ir/reaction/integrity.rs)).
ReactionSpan construction checks its two molecule projections. A subsequent
implementation must adjust these corresponding construction paths too; otherwise
moving a newly admissible molecule payload into a reaction can reject it solely
because it crossed that boundary. The existing reaction application preconditions
and checked frame transport remain at use. This does not reopen unrelated
reaction integrity checks.

#### Verification and remaining interface work

A disposable public-API probe ran against the restored code. It confirmed:

- Tetrahedral Lit(9) normalizes unchanged, but application of an identity frame
  permutation returns None.
- Tetrahedral variable domain {0, 9} normalizes to an unrestricted variable.
- A tetrahedral term containing a degree-three Apply panics during normalization.
- Two electron counts reframed by a degree-three action return None; permute
  with a three-position order leaves them unchanged.

The molecule-level indexing hazards above are source traces, not
claims that current public constructors admit those states. No visibility or
integrity bypass was added to manufacture a molecule for the probe. Its output
is recorded here; scratch/213-field-checks can be deleted.

The feasible direction is direct open field mutation with checks in the named
consumers that require agreement. The four derived electron getters and public
symmetry methods retain their current signatures. Existing Option/Result
boundaries cover normalization, canonicalization, and transport. S2c specifies
the bounded coset-operation fixes; S2g specifies frame-consumer changes.
Rust lattice operations and depiction/CoordGen handling remain unchanged.
S2b aligns the Rust unit-error name and Python join's exception behavior.
The legacy permute signature is a separate unresolved API question; current
active consumers use the already fallible FrameTransport path.
No new field setter family, view permission machinery
for these fields, or whole-molecule recheck is required by this approach. The
revised signatures must be shown before code changes; this study does not approve
them implicitly. The data-type and integrity guides must change with the eventual
contract, not in advance of it.

### Molecule-level constraint mutation

Top-level constraint mutation belongs in the editor. Molecule exposes read-only
access; entity-level constraints remain freely mutable through entity views.
Remove public Molecule::constraints_mut, ConstraintsViewMut, and
Molecule::try_modify_constraints. The checked view implemented in S2a is no
longer part of the intended API.

```rust
impl Molecule {
    pub fn constraints(&self) -> &Constraints;
    pub(crate) fn constraints_mut(&mut self) -> &mut Constraints;
}

impl MoleculeEditor {
    pub fn constraints(&self) -> &Constraints;
    pub fn constraints_mut(&mut self) -> &mut Constraints;
}
```

The crate-private Molecule accessor supplies the editor's mutable borrow by delegation;
it is not part of Molecule's public mutation surface.
The editor borrow supports the existing Constraints operations and whole-collection
assignment. It does not check writes individually or clone the collection.
Remove MoleculeEditor::push_constraint; use constraints_mut().push instead.
The crate-private Molecule::push_constraint remains for Edit execution.
Intermediate references may be invalid; publication checks the resulting molecule.
S6a changes edit to consume the molecule and move its constraint collection into
the editor; this decision needs no additional ownership mechanism.

Constraints owns the following collection primitives:

```rust
pub fn push(&mut self, constraint: Constraint);
pub fn extend(&mut self, constraints: Vec<Constraint>);
pub fn remove_at(&mut self, position: usize) -> Constraint;
pub fn retain(&mut self, predicate: impl FnMut(&Constraint) -> bool);
pub fn clear(&mut self);
pub fn take(&mut self) -> Vec<Constraint>;
```

Only extend is new in this list. It moves the supplied constraints into the
collection in order. These operations preserve duplicates and do not normalize
or validate references. Constraints have list positions, not entity ids;
addition returns no id iterator. Whole-collection replacement is assignment
through constraints_mut, not another method. Batch execution resolves handles
before insertion; explicit removal finds the last equal entry or reports
MissingEntry, then calls remove_at. Journaled removal saves its returned value
and original position. No tracked_remove is needed for that operation.
Compaction, tracked_compact, and restore are specified under Storage delegation;
restore also reinstates entries saved by explicit removal.

Python molecule-level constraint access is read-only, including nested entries.
Remove the molecule constraints property setter and mutation methods on its
ConstraintsView. Standalone Constraints remains mutable. Top-level changes use
the editor/batch editing surface; do not retain a hidden copy/edit/replace callback.
This restriction does not apply to entity-level constraint views or setters.

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

### Python accessor invalidation

**Approved.** S2l implements the counter and invalidation for molecule constraint
accessors at combine_from. The remaining accessor migration belongs to S5d;
consumption belongs to S6d. The Option
wrapper for consuming inputs is retained. Use counter for the invalidation
field and methods. Invalidate at the public operation boundary, without
classifying edits or maintaining per-entity counters.

| Operation | Existing molecule-backed accessors |
| --- | --- |
| transact / tracked_transact | Invalidate immediately before invoking Rust execution, after Python argument preparation succeeds. This includes empty and field-only batches. |
| resolve_into, resolve_into_with_report, project_into, transform_into | The same rule, including no-op, underdetermined, contradictory, and error outcomes. |
| combine_from | The same rule, even though successful append preserves existing ids. |
| Ordinary entity field/form setters and constraint collection mutation | Remain valid and observe the current value at the same entity or collection. A rejected setter leaves the value and accessor validity unchanged. |
| Consuming edit/apply/resolve/project/transform | Existing consumption contract: the owner becomes consumed and every owner-backed accessor becomes invalid. |
| Read-only operations and creation of independent results | Remain valid. |

Once execution starts, rollback restores the molecule but does not revive old
accessors. Argument-conversion, unavailable-input, or borrow failures before
execution do not invalidate them. Fresh accessors obtained after the call see
the accepted result or restored molecule. No comparison with the original
molecule is needed to decide whether a no-op or rollback should revive views.

All eight entity views, entity collections, their owner-backed iterators, and
molecule/entity constraint accessors obey this rule. Descendants inherit their
parent's captured counter; they never adopt the owner's current counter
to revive a stale parent. Check validity before reads, writes, indexing,
iteration, conversion, id access, or creating descendants. Invalid access raises
the existing InvalidatedViewError, rather than IndexError or access to a shifted
entity. The check precedes Rust indexing. Explicit independent copies remain
usable. S2l makes molecule constraint entries live read-only accessors; other
existing owned return values retain their copying policy.

Use one private u64 counter in the Python Molecule wrapper and a captured
counter in each owner-backed accessor. Advance it with ordinary addition under
the data-type guide's lifecycle-counter size assumption. Check owner availability
and counter under the same short
PyO3 borrow used for access. Borrow the receiver exclusively and advance its
counter after preparing Python inputs, then retain that borrow throughout
Rust execution. No Python callback runs between invalidation and completion.
An unwind after execution starts also leaves old accessors invalid.

Concrete wrapper layout (final shape after S6d; S2l adds counter while the
value is still non-optional):

```rust
#[pyclass]
pub struct Molecule {
    value: Option<GraphIrMolecule>,
    counter: u64,
}

#[pyclass]
pub struct AtomView {
    owner: Py<Molecule>,
    counter: u64,
    id: GraphIrAtomId,
}

pub(crate) enum AtomConstraintsStorage {
    Molecule {
        owner: Py<Molecule>,
        counter: u64,
        id: GraphIrAtomId,
    },
    Atom(Py<AtomForm>),
}
```

**Approved replacement names for existing Python enums:**

| Existing types | Approved names |
| --- | --- |
| AtomConstraintsBacking, BondConstraintsBacking, DativeBondConstraintsBacking | AtomConstraintsStorage, BondConstraintsStorage, DativeBondConstraintsStorage |
| AromaticSystemConstraintsBacking, MulticenterBondConstraintsBacking, NoncovalentBondConstraintsBacking | AromaticSystemConstraintsStorage, MulticenterBondConstraintsStorage, NoncovalentBondConstraintsStorage |
| StereoAtomConstraintsBacking, StereoBondConstraintsBacking | StereoAtomConstraintsStorage, StereoBondConstraintsStorage |
| AtomRingSizeBacking, BondRingSizeBacking, DativeBondRingSizeBacking | AtomRingSizeStorage, BondRingSizeStorage, DativeBondRingSizeStorage |

These enums select where a view reads and writes its data: an entity in a
molecule or a standalone form. Rename their containing field from backing to
storage and the macro parameter from $backing to $storage, including both
stereo-generated types. The Py<Molecule> field remains owner. These names apply
to the existing types; their source migration belongs to S5d.

Apply the same extra field to the other seven entity views, the eight entity
collections and their molecule-backed iterators, molecule ConstraintsView,
and the Molecule variants of entity-constraint and ring-size storage enums.
Standalone storage variants are unchanged. Molecule ConstraintIter retains its
collection and checks its counter on access. Entity-constraint key/item iterators
currently own copied entries or keys; they need no counter because they do not
consult the molecule. Their separate copying policy is unchanged.

The complete additional crate-visible methods on the Python Molecule wrapper
are below; fields remain private. Existing to_rust/to_rust_mut become fallible
in S6d and report ConsumedError. No public Python counter API is introduced.

```text
view_counter(&self) -> PyResult<u64>
check_access(&self, expected: u64, accessor: &'static str) -> PyResult<()>
advance_counter(&mut self) -> PyResult<()>
```

view_counter checks that the owner is available, then returns the counter.
check_access returns InvalidatedViewError if the owner is consumed or
the counter differs; the message names the accessor and Molecule. Otherwise it
returns (). advance_counter checks owner availability, then increments
counter. S5d's checks have
no consumed branch until the separately scheduled S6d ownership migration.

Every view getter/setter takes a try_borrow/try_borrow_mut of its owner, calls
check_access, then accesses Rust storage under that same borrow. For
example, after converting the Python assignment value, AtomView::set_element
writes through the entity view:

```rust
let mut owner = self.owner.bind(py).try_borrow_mut()?;
owner.check_access(self.counter, "AtomView")?;
owner.to_rust_mut()?.atom_mut(self.id).attributes_mut().element = value;
```

Root collection creation captures view_counter. Collection indexing and
iteration check the collection's saved counter, then copy it into the
returned view or iterator. Entity constraint/ring-size getters likewise copy
the validated parent's saved counter into the storage enum. __next__ checks
before advancing, including when already exhausted. Id and repr access become
fallible too. These paths must not refresh a stale counter from the owner.

In Molecule.transact, first convert and take the prepared Edits batches and
acquire the receiver's exclusive borrow. With those Rust batches ready, the
execution boundary is:

```rust
owner.advance_counter()?;
owner.to_rust_mut()?.transact(batches).map_err(molecule_apply_error)
```

The same two-step boundary applies to tracked_transact and each listed borrowed
aggregate operation. It is not placed in to_rust_mut: ordinary view setters use
that method and must not invalidate themselves. The Rust journal sees only the
contained GraphIrMolecule, so rollback cannot reset the wrapper's counter.
Explicit equality/copy implementations use molecular values, not derived equality
of the wrapper fields; input extraction must not introduce automatic copies.

This adds one integer per owner/accessor, a comparison per access, and one
increment per aggregate mutation call. It adds no molecular copy, per-entity
map, or Rust graph-IR version field. Its deliberate cost to callers is that even
a failed or field-only transaction requires obtaining fresh views; ordinary
property editing does not. Neither counter nor consumed-state metadata
participates in molecular equality or serialization. Independent copies have
independent owner lifecycles.

Acceptance cases: deleting an earlier entity never retargets an old view;
field-only and empty transactions invalidate; rollback and late integrity
failure leave old views invalid and fresh ones usable; pre-execution rejection
preserves views; ordinary setters keep views live; stale nested access and
iterator advancement fail; independent copies survive.

### Edits — one accumulated sequence

```rust
impl Edits {
    pub fn new() -> Self;
    pub fn push(&mut self, edit: Edit);
    pub fn add_atom(&mut self, attributes: AtomForm) -> AtomHandle;
    pub fn add_dative_bond(
        &mut self, donors: Vec<AtomHandle>, acceptor: AtomHandle,
        attributes: DativeBondForm,
    ) -> DativeBondHandle;
    // Other existing typed creation/update/removal methods remain.
}
```

Retain one accumulated sequence, one initial-host namespace, and one creation
namespace per kind. There is no independent-batch composition, append(Edits),
extend(Edits), or caller-offset API. Multiple batches are handled by transaction
processing under the separate-apply rules above. This is settled, not deferred.

Edits constructs a sequence without executing it or validating it against a
molecule. Id names an entity at this batch's application entry, not at the
start of an enclosing multi-batch transaction. New names a same-kind creation
within this sequence. update_* derives old/new entries from the supplied
current form; it does not simulate earlier entries. Raw push and FromIterator
preserve supplied handles and do not rebase independently constructed batches.

Dative addition and removal payloads distinguish donors from acceptor, as the
editor and DativeBondDelta already do. Replace the combined atoms vector whose
last element denotes the acceptor. Keep donor vectors owned in Edits; direct
editor operations borrow slices. S3a–S3e migrate construction, execution, saved
undo records, lowering, and bindings together.

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

## Mutation vocabulary and mutable views

The participant, Edit, Delta, and Undo contracts below remain in place.
S2i1 removed the public checked molecule-level constraint view. S2i2–S2i4
implemented uniform attribute access and typed editor storage. S2i5 replaced the
shared mutable-view types with separate molecule/editor types; it preserved that
storage and the completed caller migrations.

### Mutable-view structures and access

Molecule and MoleculeEditor expose mutable entity views for all eight kinds:
atom, localized bond, dative bond, aromatic system, multicenter bond,
noncovalent bond, stereo atom, and stereo bond. Every view provides a mutable
borrow of the complete entity form; all attribute fields and entity-level
constraints remain directly assignable. The same assignment semantics apply in
Python. There are no exceptions for electron counts or stereo configuration.

Use four separate view families:

| Access | Molecule | MoleculeEditor |
| --- | --- | --- |
| Immutable | *View<'a> | *EditorView<'a> |
| Mutable | *ViewMut<'a> | *EditorViewMut<'a> |

Within each row, editor read methods are a subset of the molecule interface,
with identical signatures, lifetimes, and semantics for shared methods. Mutable
editor views additionally provide structural mutation. Both mutable families
expose the complete form through attributes_mut. There is no const parameter,
runtime permission flag, checked field setter, or modify/try_modify callback.
Do not introduce individual charge_mut, element_mut, electrons_mut, or coset_mut
accessor families to replace ordinary field assignment.

Attribute access covers these exact payloads:

| Entity | Mutable attribute borrow |
| --- | --- |
| Atom | &mut AtomForm |
| Localized bond | &mut BondForm |
| Dative bond | &mut DativeBondForm |
| Aromatic system | &mut AromaticSystemForm |
| Multicenter bond | &mut MulticenterBondForm |
| Noncovalent bond | &mut NoncovalentBondForm |
| Stereo atom | &mut StereoAtomForm |
| Stereo bond | &mut StereoBondForm |

A view borrows its owner, is neither Clone nor Copy, owns no draft or journal,
and does no work on drop. Views are receivers, not method or callback arguments.
Invalid entity ids continue to panic at the accessor entry points. Existing
copy-on-write storage behavior remains; the view adds no recovery copy.
Entity-specific read access follows [local getters](#local-getters).

Each mutable overlay view stores its entity id and a mutable borrow of its
owning typed set. The molecule and editor types have the same fields but distinct
method sets. Atom views borrow their form; localized-bond views also store the
endpoint ids. S2i5 gives the exact structures and accessors.

Structural changes can break reference, uniqueness, or incidence integrity, so
the editor publication boundary checks them. Attribute assignment does not change
those structural links. No conversion from a molecule mutable view to an editor
mutable view is exposed. Internal batch/undo execution obtains editor mutable
views through crate-private Molecule accessors.

### Participant mutation

Settled 2026-09-19: participant-list mutation belongs on the entity mutable view.
The methods below remain editor-only in the public API; internal batch and undo
execution use the same methods. Attribute access is unrestricted on both
molecule and editor views; structural mutation is a separate boundary.
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

StereoAtomViewMut and StereoBondViewMut expose
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
their existing errors. Undo restores its saved value without an expected-post-value
comparison.

### Undo variants

Each replacement Undo carries the affected entity's typed id and its actual
previous stored value, with resolved ids and exact sequence order.

| Undo variant | Saved field | Saved type |
| --- | --- | --- |
| RestoreAromaticSystemAtoms | atoms | Vec<AtomId> |
| RestoreMulticenterBondAtoms | atoms | Vec<AtomId> |
| RestoreNoncovalentBondAtoms | atoms | [AtomId; 2] |
| RestoreDativeBondDonors | donors | Vec<AtomId> |
| RestoreDativeBondAcceptor | acceptor | AtomId |
| RestoreStereoAtomSite | site | AtomId |
| RestoreStereoAtomLigands | ligands | Vec<StereoLigand> |
| RestoreStereoBondSite | site | BondId |
| RestoreStereoBondLigands | ligands | Vec<StereoLigand> |

Replay uses ordinary storage replacement. Local target-access guards provide the
settled no-panic undo contract; undo does not validate an expected post-value.
The other factor, attributes, and constraints are preserved. These undos do not
use graph-core row restoration or reference uncompaction: they restore a saved
value on an existing row without changing its id.

### Delta and reaction scope

Structural replacement belongs to Edits and direct editor mutation. Deltas retain
addition, removal, field modification, and constraint modification. Reactions
express changes to atom lists, donors, acceptors, sites, and ligands through
explicit Add/Remove entries. The two entries are not inferred to be a replacement.

Replacement Deltas were withdrawn on 2026-09-26 because no consumer currently needs
them. Their normalization rules are tractable, and a representation need not
round-trip every operation; neither issue alone warrants removing the vocabulary.
The demonstrated need is replacement Edits for transformers, including hydrogen
folding/unfolding. That work does not require equivalent Delta variants.

Keep ReactionSpan's shared-incidence representation and existing conversion
semantics. Correspondence induction leaves stereo entities of different determined
kinds unmatched; superimposition represents them as Removed and Added. An
incompatible Modified span is rejected. Reaction application rejects a configuration
ModifyField that changes determined kind; it does not infer removal/addition.

The approved reaction method names and other semantic changes remain. S3c/S3d
remove only replacement Delta support; replacement Edits, Undo, and Edit DSL verbs
are unaffected.

### Local getters

The separate families and matching-interface rule are settled. The getters below
are implemented for all eight entity families in S2j.
Existing immutable Molecule view methods retain their signatures and semantics,
except for the already agreed borrowed ligand_ids return.

AtomView and AtomEditorView have private id and attributes fields, exposed by
id() -> AtomId and attributes() -> &'a AtomForm. AtomView stores molecule first,
then id, then the attribute borrow supplied at construction. AtomEditorView stores
id and the attribute borrow; its crate-private new takes those two arguments.
Both mutable atom views retain id(), attributes(), and attributes_mut(). Matching
atom-view accessors use #[inline] consistently. Python keeps its read-only id
property and existing attribute properties; ownership and copying are unchanged.

The four localized-bond views use the same id/attribute accessor pattern.
BondView stores molecule, id, atoms, and the supplied attribute borrow, in that
order. BondEditorView stores id, atoms, and the attribute borrow; its crate-private
new takes those three arguments. All four expose atom_ids() -> [AtomId; 2],
copying only the two ids, without allocation. All fields are private, and matching
accessors use #[inline]. Python properties and setters are unchanged.

DativeBondViews stores only molecule; its constructor takes that borrow alone.
DativeBondViews and DativeBondView access the set through raw_dative_bonds().
DativeBondView stores molecule and id. DativeBondEditorView stores dative_bonds:
&DativeBonds and id; both mutable views store dative_bonds: &mut DativeBonds and
id. Fields and constructor arguments put the owning borrow first. All fields are
private. All four expose id(), attributes(), donor_ids(), acceptor_id(), and
atom_ids(), with #[inline] on matching accessors. Donor and combined-atom getters
return lazy exact-size iterators in stored donor order; atom_ids appends the
acceptor. Readonly references/iterators retain 'a; mutable-view accessor borrows
last for the method borrow. Both mutable views expose attributes_mut(). Public
entity lookup preserves invalid-id panics; the namespace's get remains optional.

AromaticSystemViews stores only molecule; its constructor takes that borrow alone.
AromaticSystemViews and AromaticSystemView access the set through raw_aromatic_systems().
AromaticSystemView stores molecule and id. AromaticSystemEditorView stores
aromatic_systems: &AromaticSystems and id; both mutable views store
aromatic_systems: &mut AromaticSystems and id. All fields are private, with the
owning borrow first in fields and constructors. All four expose id(), attributes(),
and atom_ids(), with matching #[inline] annotations. atom_ids delegates to the
set's lazy exact-size iterator. Readonly borrows retain 'a; mutable-view accessor
borrows last for the method borrow. Both mutable views retain unrestricted
attributes_mut(). Molecule-dependent queries retain their existing behavior.

MulticenterBondViews stores only molecule; its constructor takes that borrow alone.
MulticenterBondView stores molecule and id. Both access the set through
raw_multicenter_bonds(). MulticenterBondEditorView stores multicenter_bonds:
&MulticenterBonds and id; both mutable views store multicenter_bonds:
&mut MulticenterBonds and id. All fields are private, with the owning borrow first
in fields and constructors. All four expose id(), attributes(), and lazy exact-size
atom_ids(), with matching #[inline] annotations. Readonly borrows retain 'a;
mutable-view accessor borrows last for the method borrow. Both mutable views
retain unrestricted attributes_mut(). Public lookup preserves invalid-id panics.

NoncovalentBondViews stores only molecule; its constructor takes that borrow alone.
NoncovalentBondView stores molecule and id. Both access the set through
raw_noncovalent_bonds(). NoncovalentBondEditorView stores noncovalent_bonds:
&NoncovalentBonds and id; both mutable views store noncovalent_bonds:
&mut NoncovalentBonds and id. All fields are private, with the owning borrow first
in fields and constructors. All four expose id(), attributes(), and
atom_ids() -> [AtomId; 2], with matching #[inline] annotations. Readonly attribute
borrows retain 'a; mutable-view accessor borrows last for the method borrow.
Both mutable views retain unrestricted attributes_mut(). Public lookup preserves
invalid-id panics; atom_ids preserves stored endpoint order without allocation.

StereoAtomViews and StereoBondViews store only molecule. Their singular molecule
views store molecule and id, and read through raw_stereo_atoms() and
raw_stereo_bonds(). StereoAtomEditorView stores stereo_atoms: &'a StereoAtoms
followed by id: StereoAtomId; StereoBondEditorView stores stereo_bonds:
&'a StereoBonds followed by id: StereoBondId. Both mutable families hold the
corresponding mutable set borrow followed by id. Fields are private and
constructors take the owning borrow first. All four views expose id(),
attributes(), site_id(), and ligand_ids(); both mutable families retain
attributes_mut() and constraints(). Readonly attributes/frames retain 'a;
mutable-view accessor borrows last for the method borrow. Invalid-id lookup
panics are preserved. Macros share repeated namespace, readonly, and mutable
implementations; entity-specific navigation and incidence remain explicit.
The additional stereo getters for editor/mutable views are implemented in S2j.

The readonly storage design retains a separately supplied attribute borrow for
atoms and localized bonds. Relation views retrieve their attributes from the
owning entity set through molecule and id. All eight entity families now follow
these structures.

Mutable molecule/editor views both provide id(), attributes(), attributes_mut(),
and constraints() with the exact S2i5 signatures. Their accessor borrows last for
the method borrow. Immutable getters return references tied to the stored borrow
'a, matching the existing immutable molecule views; do not shorten them to &self.
Iterator lifetimes follow the same rule: +'a for immutable views and +'_ for
mutable views.

Both mutable families and immutable editor views expose these local methods.
Every method takes &self. In the table, &Form denotes &'a Form for
immutable views and &Form tied to &self for mutable views.

| Entity | Getter signatures |
| --- | --- |
| Atom | element() -> &ElementForm; isotope_mass() -> &IsotopeMassForm; charge() -> &NumForm; implicit_hydrogens() -> &NumForm; lone_pairs() -> &NumForm; unpaired_electrons() -> &UnpairedElectronsForm |
| Localized bond | atom_ids() -> [AtomId; 2]; order() -> &NumForm; charge() -> &NumForm; unpaired_electrons() -> &UnpairedElectronsForm |
| Dative bond | donor_ids() -> impl ExactSizeIterator<Item = AtomId>; acceptor_id() -> AtomId; atom_ids() -> impl ExactSizeIterator<Item = AtomId>; donor_count() -> usize; atom_count() -> usize; order() -> &NumForm |
| Aromatic system and multicenter bond | atom_ids() -> impl ExactSizeIterator<Item = AtomId>; atom_count() -> usize; electrons() -> &ElectronCountsForm; electron_count() -> NumForm; charge() -> &NumForm; unpaired_electrons() -> &UnpairedElectronsForm |
| Noncovalent bond | atom_ids() -> [AtomId; 2]; kind() -> &NoncovalentBondKindForm |
| Stereo atom and stereo bond | site_id() -> AtomId / BondId respectively; kind() -> StereoKind; coset() -> &StereoCoset; ligand_count() -> usize; ligand_position(id: AtomId) -> Option<StereoLigandPosition>; ligand_ids() -> &[StereoLigand]; atom_ligand_ids() -> impl Iterator<Item = AtomId>; implicit_hydrogen_atom_ids() -> impl Iterator<Item = AtomId>; lone_pair_atom_ids() -> impl Iterator<Item = AtomId>; atom_ligand_count() -> usize; implicit_hydrogen_count() -> usize; lone_pair_count() -> usize |

Dative atom_ids yields donors in stored order followed by the acceptor.
Aromatic/multicenter electron_count sums the stored literal contributions,
otherwise returns NumForm::Undetermined. It does not check contribution-vector
length. Stereo kind/coset match the existing immutable molecule views: they panic
for an undetermined configuration. Optional access remains available through
attributes.configuration.kind()/coset(); no new configuration() getter is needed.
ligand_position takes id: AtomId and finds the first actual-atom ligand with that id.
Filters preserve stored order; virtual-ligand id accessors yield their bearing ids.

The matching-interface decisions below are implemented, including stereo
site/frame access and the additional getters for all entity families:

| Earlier editor proposal | Existing immutable Molecule getter | Implemented decision |
| --- | --- | --- |
| kind() -> Option<StereoKind>; coset() -> Option<&StereoCoset> | kind() -> StereoKind; coset() -> &'a StereoCoset | Match the molecule signatures and undetermined-configuration panic. |
| constraints() -> &EntityConstraintsForm on an immutable editor view | constraints() -> EntityConstraintsView<'a> | Do not add this immutable-editor method. Read stored constraints through attributes.constraints. Both mutable families retain constraints() -> &EntityConstraintsForm. |
| ligand()/ligands() and ligand-kind filters returning StereoLigand | Corresponding methods return StereoLigandView<'a> | Do not add raw-value methods under these names. Raw ordered ligands remain available through borrowed ligand_ids(); retain the id/count methods listed above. |

The existing molecule ligand/constraint view APIs are retained.
Resolved entity views, neighbors, valence, induced bonds, stereo-bond endpoint
lookup, and frame-transport queries stay outside the editor getter additions.
No new checked access, error type, or publication requirement is introduced.

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

Molecule owns seven topology operations because Graph cannot maintain the
parallel atom/bond attributes: add_atom, add_atoms, add_bond, add_bonds,
remove_topology, tracked_remove_topology, and restore_topology. The other
internal Molecule mutations delegate one operation to one owning component.
Editor and transaction execution use these private methods whether they own or
borrow the molecule. Delegation does not combine overlay, topology, and
constraint changes, resolve Edit handles, check offered old values, or record
undo. Their composition remains explicit in edit execution.
All editor and transaction writes, including undo, go through Molecule mutation
methods or mutable views obtained from Molecule. Molecule owns view construction;
these consumers do not borrow its fields to construct views or directly mutate
Graph, attribute vectors, typed sets, or copy-on-write storage. Structural
replacement uses *EditorViewMut, which delegates to the typed sets.
This boundary adds no second mutation vocabulary or public mutable access.

Graph retains add_node/add_edge for individual additions and uses extend for
bulk addition, matching relation sets: extend_nodes adds isolated nodes,
extend_edges adds edges between existing nodes, and extend adds nodes and edges
together. The bare verb covers both topology components, consistently with the
existing graph nomenclature. Edge endpoints supplied to extend use the
resulting graph's node ids, including the appended nodes; there is no separate
handle namespace. Existing ids stay unchanged and each added kind occupies a
contiguous block. Return owned, allocation-free, exact-size iterators over those
blocks, following the existing typed node_ids/edge_ids iterator interfaces:

```rust
impl Graph {
    pub fn extend_nodes(
        &mut self, count: usize,
    ) -> impl ExactSizeIterator<Item = NodeId> + use<>;
    pub fn extend_edges(
        &mut self, edges: &[[NodeId; 2]],
    ) -> impl ExactSizeIterator<Item = EdgeId> + use<>;
    pub fn extend(
        &mut self, node_count: usize, edges: &[[NodeId; 2]],
    ) -> (
        impl ExactSizeIterator<Item = NodeId> + use<>,
        impl ExactSizeIterator<Item = EdgeId> + use<>,
    );
}
```

The returned iterators own their bounds and borrow neither the graph nor the
input edge slice; use<> makes that absence of captures explicit. The mutation
is complete before return, irrespective of whether the caller consumes the
iterators. Empty additions return empty iterators for the corresponding kind.

No tracked addition variants are needed: existing ids map identically and new
ids have no predecessor. Batch execution records added ids for undo and derives
correspondence when requested. AddAtoms/AddBonds can use the respective bulk
operations; combine_from can use extend. Node-only addition extends
adjacency offsets without rebuilding existing edges. Combined addition builds
the final adjacency once, avoiding an intermediate copy-on-write detachment for
node addition followed by an edge rebuild. Single additions delegate to their
bulk counterparts. Measurements are recorded under S3f.

Relation sets retain individual add and gain public extend. Each extend takes
&mut self plus the batch below and returns an owned, allocation-free
`impl ExactSizeIterator<Item = RelationId>`. Returned iterators borrow neither
the receiver nor the supplied participant slices. Precise captures must exclude
those lifetimes while accounting for the enclosing type/const parameters.

| Relation set | extend batch argument |
| --- | --- |
| FixedRelationSet<P, D, N> | Vec<([P; N], D)> |
| VarRelationSet<P, D> | Vec<(&[P], D)> |
| FixedFixedBirelationSet<L1, N1, L2, N2, D> | Vec<([L1; N1], [L2; N2], D)> |
| FixedVarBirelationSet<L1, N1, L2, D> | Vec<([L1; N1], &[L2], D)> |
| VarVarBirelationSet<L1, L2, D> | Vec<(&[L1], &[L2], D)> |

Move payloads into storage; copy variable participant slices into the packed
buffers. Preserve existing rows, supplied row/factor order, and multiplicity;
new ids form one contiguous block in input order. Append the batch and rebuild
incidence once, completing mutation before returning. An empty batch leaves
storage unchanged and returns an empty iterator. Match individual add's
construction contract, including no external graph-membership validation.
No new entry types, IntoIterator inputs, tracked additions, or add_many methods.

All six typed overlay sets gain crate-private extend with the same domain
arguments as their corresponding Molecule bulk additions below, returning an
owned exact-size iterator of their typed ids. They convert graph-IR ids and
delegate to relation extend; do not loop over single additions and rebuild
incidence for each row. Individual set add remains available.

Constraints owns compaction and restoration, aligned with the other storage APIs:

```rust
impl Constraints {
    pub fn compact(&mut self, compaction: &MoleculeCompaction);
    pub fn tracked_compact(
        &mut self,
        compaction: &MoleculeCompaction,
    ) -> CascadedConstraints;
    pub fn restore(&mut self, changes: &CascadedConstraints);
}
```

Rename the existing compact_with_update to tracked_compact without changing its
semantics: translate surviving entity ids, drop entries referencing removed
entities (including a whole compound entry when a subtree references one), and
record removed or rewritten constraints at their original positions. Keep the
existing CascadedConstraints representation and surviving list order.

Move transaction-level restored_constraints into Constraints::restore. Restore
removed entries and old values at their original positions, retaining local
size/access guards but removing comparison against the recorded post-value.
Matching history restores the pre-state; manipulated history must not panic but
has no correctness guarantee. No tracked_restore is needed: restoration uses
the saved changes, and any entity-id translation is already supplied by the
original compaction.

Plain removal calls compact. Transactional removal calls tracked_compact once
on the actual collection and saves the returned changes in the bundled undo.
Remove the current whole-collection recovery clone and second compaction used
to obtain those changes. Undo restores entity storage before restoring constraints.

Retain graph-IR coordination while replacing the editor's storage implementations.
Overlay operations delegate through the typed entity sets that own copy-on-write
relation storage. This inventory records the implementation obligations for the
API defined above; it is not a staged implementation plan.

| Editor surface | Current implementation and required change |
| --- | --- |
| add_atom, add_bond | Move implementations to crate-private Molecule methods; editor methods delegate and batch execution calls the same implementations. Retain Graph addition and parallel attribute-array updates; direct additions record no correspondence. |
| add_dative_bond, add_aromatic_system, add_multicenter_bond, add_noncovalent_bond, add_stereo_atom, add_stereo_bond | Move implementations to crate-private Molecule methods; editor methods delegate and batch execution calls the same implementations. Replace mutable-entry-vector insertion with relation add through the typed set; direct additions record no correspondence. |
| atom_mut, bond_mut | Already mutate copy-on-write attribute arrays. Retain graph-IR ownership of these arrays. |
| Mutable attribute views for all six overlay kinds | Replace entry-vector materialization with typed-set access to relation payload mutation. Keep attributes and participants accessible in the entity view. |
| constraints_mut, inline constraints through attribute views | Remove editor push_constraint; use constraints_mut().push. Keep crate-private Molecule::push_constraint for Edit execution. Retain editor-only top-level mutable access and shared entity-attribute access; graph-core does not interpret molecular constraints. |
| Each overlay remove_* / tracked_remove_* pair | Crate-private Molecule methods remove entries from that owning set only; tracked forms return its typed compaction. Complete editor/batch removal separately assembles MoleculeCompaction and compacts constraints. Use the set's removal implementation; remove public tracked direct mutation. |
| Topology remove / tracked_remove | Rename the editor operation to remove_topology. Crate-private Molecule remove_topology / tracked_remove_topology mutate Graph and atom/bond attributes only; tracked removal returns GraphCompaction. Editor/batch execution separately compacts each overlay set, assembles MoleculeCompaction, and compacts constraints. Graph nomenclature is unchanged. |
| Internal undo-addition removal and restore_* methods | Undo additions with crate-private Molecule untracked removal only, without overlay or constraint compaction. Undo removals with restore_topology for Graph and atom/bond attributes, separate overlay topology-id/row restoration, then constraint restoration. S4a lists the restoration interfaces. |
| Overlay and constraint compaction | Crate-private Molecule compact_<overlay> methods install the returned set and return its row mapping. Constraint delegates compact_constraints / tracked_compact_constraints mutate the collection in place. S4b lists their interfaces. |
| apply, transact, and tracked counterparts | Replace the current lifecycle with the editor/transaction API above. Both execution paths share graph-IR mutation operations; handle resolution and forward preconditions remain batch concerns, and undo capture/replay remains transactional. |
| snapshot, try_build, build, and tracked counterparts | Remove relation-row rebuilding at publication. Replace editor publication with finish, Molecule::apply, or scoped commit, and replace snapshot with checked probe access. Fresh MoleculeBuilder retains asserted build. |
| Overlay reads and views; six internal *_equiv methods | Remove storage-wrapper dispatch. Use typed-set accessors and graph-core participant comparison where applicable; retain graph-IR frame transport and payload comparison. |

The crate-private Molecule addition interfaces are:

```rust
pub(crate) fn add_atom(&mut self, attributes: AtomForm) -> AtomId;
pub(crate) fn add_bond(&mut self, first: AtomId, second: AtomId, attributes: BondForm) -> BondId;
pub(crate) fn add_dative_bond(&mut self, donors: &[AtomId], acceptor: AtomId, attributes: DativeBondForm) -> DativeBondId;
pub(crate) fn add_aromatic_system(&mut self, atoms: &[AtomId], attributes: AromaticSystemForm) -> AromaticSystemId;
pub(crate) fn add_multicenter_bond(&mut self, atoms: &[AtomId], attributes: MulticenterBondForm) -> MulticenterBondId;
pub(crate) fn add_noncovalent_bond(&mut self, atoms: [AtomId; 2], attributes: NoncovalentBondForm) -> NoncovalentBondId;
pub(crate) fn add_stereo_atom(&mut self, site: AtomId, ligands: &[StereoLigand], attributes: StereoAtomForm) -> StereoAtomId;
pub(crate) fn add_stereo_bond(&mut self, site: BondId, ligands: &[StereoLigand], attributes: StereoBondForm) -> StereoBondId;
pub(crate) fn push_constraint(&mut self, constraint: Constraint);
```

Bulk additions have identical interfaces on crate-private Molecule methods and public
MoleculeEditor methods. Each takes &mut self and the listed owned batch vector;
the return is `impl ExactSizeIterator<Item = Id> + use<>`, with Id from the last
column. The editor delegates to Molecule. Variable participant lists remain
borrowed slices, fixed-size lists remain arrays, and attributes move into storage.
The returned iterator owns its bounds and retains no receiver/input borrow.

| Method | Batch argument | Returned Id |
| --- | --- | --- |
| add_atoms | Vec<AtomForm> | AtomId |
| add_bonds | Vec<([AtomId; 2], BondForm)> | BondId |
| add_dative_bonds | Vec<(&[AtomId], AtomId, DativeBondForm)> | DativeBondId |
| add_aromatic_systems | Vec<(&[AtomId], AromaticSystemForm)> | AromaticSystemId |
| add_multicenter_bonds | Vec<(&[AtomId], MulticenterBondForm)> | MulticenterBondId |
| add_noncovalent_bonds | Vec<([AtomId; 2], NoncovalentBondForm)> | NoncovalentBondId |
| add_stereo_atoms | Vec<(AtomId, &[StereoLigand], StereoAtomForm)> | StereoAtomId |
| add_stereo_bonds | Vec<(BondId, &[StereoLigand], StereoBondForm)> | StereoBondId |

Molecule topology bulk additions update Graph and the corresponding attribute
vector together. Overlay bulk additions delegate to the owning set's extend.
Existing ids remain unchanged; returned ids follow input order. These are eager,
sharp mutations with the same checks as the corresponding individual addition;
they do not acquire Edit preconditions, a journal, or publication checks.
The plural Edits constructors continue to construct their existing variants;
these storage additions do not authorize new bulk Edit variants or coalescing
separate edits into a different execution unit.

Each row below specifies a crate-private Molecule removal pair. Both methods take
&mut self and the listed arguments. The bare method returns (); the tracked_
form returns the listed mapping. Neither form mutates top-level constraints.

| Bare method | Arguments after &mut self | Mutated storage | tracked_ return |
| --- | --- | --- | --- |
| remove_topology | atoms: &[AtomId], bonds: &[BondId] | Graph and atom/bond attribute vectors, including incident-bond removal | GraphCompaction |
| remove_dative_bonds | ids: &[DativeBondId] | DativeBonds | Compaction<DativeBondId> |
| remove_aromatic_systems | ids: &[AromaticSystemId] | AromaticSystems | Compaction<AromaticSystemId> |
| remove_multicenter_bonds | ids: &[MulticenterBondId] | MulticenterBonds | Compaction<MulticenterBondId> |
| remove_noncovalent_bonds | ids: &[NoncovalentBondId] | NoncovalentBonds | Compaction<NoncovalentBondId> |
| remove_stereo_atoms | ids: &[StereoAtomId] | StereoAtoms | Compaction<StereoAtomId> |
| remove_stereo_bonds | ids: &[StereoBondId] | StereoBonds | Compaction<StereoBondId> |

Private topology removal leaves overlays to their separate compaction operations.
The editor's public remove_topology remains a complete cascading removal: call
tracked_remove_topology, call the six Molecule compact_* overlay
delegates with its GraphCompaction, assemble MoleculeCompaction::new from the
graph and six row mappings, then call compact_constraints. Batch execution uses
that same sequence and retains the mappings it needs. Explicit overlay removal instead
assembles MoleculeCompaction from the one changed row mapping and identity
mappings for the unchanged kinds. No other overlay set changes in that path.
The typed-set compact interfaces are recorded in S1a–S1c. Their returned row
mappings are required to update dependent constraints and handles.

These compositions may remain explicit in their consumers. Do not introduce
combined helpers, output parameters, or recording modes solely to share the
sequencing. Untracked-to-tracked delegation is the current implementation, not
a settled efficiency requirement. Remove the seven remove_added_* adapters;
each Undo match arm extracts saved ids and calls the crate-private Molecule untracked
removal for the added kind directly. Under matching
history, reverse replay has undone later dependencies and returned added entries
to their appended positions; earlier ids do not change. No separate overlay or
constraint compaction, combined MoleculeCompaction, or correspondence update is
needed for addition undo. Retain local panic guards for manipulated history,
without validating that it is matching history or requiring a trailing block.

Restoration follows the same component boundaries:

```rust
// Crate-private Molecule method.
pub(crate) fn restore_topology(
    &mut self,
    compaction: &GraphCompaction,
    atoms: Vec<RemovedAtom>,
    bonds: Vec<RemovedBond>,
);
```

This method calls Graph::restore and restores the atom/bond attribute vectors.
It does not restore overlays or constraints. Undo execution then translates
surviving overlay topology ids and restores removed rows through the separate
Molecule delegates listed in S4a, then calls restore_constraints. For
overlay-only removal, topology ids did not change, so only the affected set's
restore and constraint restoration are needed. The existing typed-set restore
signatures are recorded in S1a–S1c; add no combined overlay-restoration helper.
Retain the matching-history and panic-free misuse contracts below.

Internal Molecule mutation methods and mutable accessors use pub(crate); storage
fields stay private. Ordinary internal visibility is private or pub(crate), not
pub(super) or path-restricted visibility. Module-test support may use pub(super);
other exceptions require a specific, narrow justification.

Single-entry execution belongs to crate-private Molecule methods: apply_edit,
apply_edit_with_undo, and apply_undo, retaining the existing names. Edit and Undo
remain data enums; Edits retains its construction API. Editor and Transaction
own their batch loops and call these methods on the owned or borrowed Molecule.
Transaction owns journal storage and reverse iteration.

Edit-specific checks belong in batch execution for every variant: resolve
handles and check all preconditions before writing, then use the shared mutation
surface. Remove the sixteen apply_modify_*_field / apply_modify_*_constraint
helpers from the edit interpreter; do not move or recreate them on Molecule.
Field changes use ordinary assignment through attributes_mut after checking the
offered old value. Entity-level constraints use the stored constraint collection.
Journaled execution additionally retains the inverse change. This follows the
existing separation between edit checks and mutation for additions/removals;
no new per-field mutation methods are required.

The following covers every current non-constraint Edit variant and all nine
planned structural replacement variants:

| Edit family | Mutation implementation |
| --- | --- |
| AddAtoms, AddBonds, six overlay additions | Shared crate-private Molecule add_* methods, with bulk topology addition for the existing AddAtoms/AddBonds variants. |
| RemoveTopology, six overlay removals | Crate-private Molecule removal primitives followed by the separate overlay/constraint compaction steps described above; batch execution retains the combined mapping for handle updates. |
| Modify*Field for all eight entity kinds | Assignment through the existing/shared mutable attribute access. |
| Nine Replace* variants for atoms, donors, acceptors, sites, and ligands | Structural methods on *EditorViewMut obtained through crate-private Molecule access; those views delegate to their owning sets. |

The shared attribute access and planned structural replacements complete this
coverage; no additional Molecule modification family is needed. Neither optional
mutable undo-output parameters nor an attributes-and-overlays grouping is approved.
Constraint undo capture is the ordinary returned value of tracked_compact, called
after assembling MoleculeCompaction; it is not an argument to a removal primitive.

Plain editor removal still needs compaction internally even when it returns no
witness. Atom/bond addition already delegates to Graph, and entity-set frame
transformations already use storage participant permutation. Undo of additions
uses only the private untracked removal primitives; undo of removals needs the
restoration capability below.

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
| Dative restoration | Pass the separately stored donors and acceptor to the owning set; no extraction from a combined sequence remains. |
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

Retain the existing Undo representation except for the planned structural
replacement variants and replacement of ApplyEdit with nine explicit constraint
undo variants in S4b. A compressed journal redesign remains separate. The large
enum layout is a recorded performance liability. Avoid avoidable payload clones by consuming
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
the transaction's journal, status, and scope guard are private and cannot be
substituted by the caller. There is no partial-edit recovery record.

Field execution accepts normalized equality with the offered old value. Reversing
that accepted change is consistent with rollback modulo normalized_eq: restoring
a literal in place of an equivalent singleton literal set is permitted. Retain
this comparison and undo behavior; no extra capture of the original encoding is
required. This is a guarantee for the operation's own matching history, not a
correctness guarantee for manipulated journals.

The edit is the execution and undo unit. Resolve all handles, check the whole
edit's preconditions, and prepare its required input data before its first write.
Execute the bundled mutation, then record its bundled Undo. An entity-constraint
None-to-None edit makes no change and records no undo entry. If a later edit or
commit fails, reverse the completed edits. Do not decompose a bundled edit into
separately checked mutations, partial undos, or per-storage-step progress records.

The inspected AddBonds and removal paths already check the complete edit before
mutation; field and constraint changes compare before assignment. No reachable
partial-edit failure was established by that inspection. Allocation failure and
internal programming panics during mutation do not carry a restoration guarantee.
Do not add per-apply unwind guards, reserve every allocation for panic recovery,
or inject panics between storage writes to justify additional machinery. A newly
found input-dependent failure belongs before the write, or in a specific corrected
operation; it does not authorize a general partial-write recovery design.

The scope guard covers explicit rollback, checked application/commit failure,
callback abandonment, callback error, and callback unwinding between completed
edits, including after commit was requested. Accept only after successful commit
and callback Ok; retain the journal until that decision. Do not convert panics
into chemistry errors. Process abort has no return-state guarantee. The separate
settled requirement that restoration with manipulated undo data must not panic
remains unchanged.

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
transaction borrowing the existing molecule. Both use the same crate-private edit
kernel. There is no copied draft, move-out, empty placeholder, or ownership mode.

Eleven transaction cases passed: accepted commit; ordinary domain rejection;
explicit rollback; later batch error; failed final integrity; rejected probe
followed by repair; forgotten capability; callback error after commit request;
callback panic after commit request; rollback after an application failure without
clearing the abort latch; and a two-write edit panic caught inside the callback.
The last case restores before another probe, rejects commit, and returns Aborted.
That injected partial-edit panic tested a stronger contract than the selected
design requires. It does not establish a production need for per-apply cleanup;
that mechanism is excluded. Callback recovery through the outer scope guard
remains part of the selected contract.

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

Numbered children (for example S4b1) divide an unfinished subitem into executable
units. Its unsuffixed label denotes the complete group, not an additional task;
a dependency on that label requires every child. Completed and reverted records
keep their original labels. Each child inherits its group's explicit contracts
and carries focused verification of its own change. Where an enum or signature
migration must span children, the green boundary is stated explicitly; do not
insert compatibility APIs or incomplete match arms between them.

### S0 — Baseline and solution vocabulary

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

- **S0c — reverted.** The attempted removal of Molecule mutation is undone.
  Its replacement is S2 below, which supplies the new views and migrates callers
  before removing callbacks. The former proposed S0d/S0e work is incorporated
  into S2i/S2j; those labels are not executable subitems.

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
  zero-based offset in the selected atom sequence; S2j gives public editor views
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
  Set::compact(&self, graph: &GraphCompaction) -> (Self, Compaction<Id>)
  ```

  These methods are `pub(crate)`; S2j exposes public mutation through the editor
  views. They delegate to `VarRelationSet`, retaining the owning set's
  copy-on-write storage. Focused tests cover order, draft duplicates, incidence,
  attribute preservation, compaction, and matching-history restoration.
  The implementation matches the interfaces above. Focused tests, strict
  graph-IR Clippy, nightly formatting, and `git diff --check` pass.
- **S1b — completed 2026-09-24** (`ir::dative`, `ir::noncovalent`; additive) Add the same ownership
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
  DativeBonds::compact(&self, graph: &GraphCompaction) -> (Self, Compaction<DativeBondId>)
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
  NoncovalentBonds::compact(&self, graph: &GraphCompaction) -> (Self, Compaction<NoncovalentBondId>)
  ```

  These are crate-private set methods. S2j exposes the domain-level
  editor-view methods with `AtomPosition` after S2i replaces the storage wrappers.

  Implemented the recorded interfaces through `FixedVarBirelationSet` and
  `FixedRelationSet`, preserving copy-on-write ownership. Tests cover independent
  donor/acceptor changes, ordered and duplicate atoms, empty donors, fixed
  endpoints, incidence, attribute preservation, and removal/compaction followed
  by restoration. All 44 focused tests, strict graph-IR Clippy, nightly formatting,
  and `git diff --check` pass.
- **S1c — completed 2026-09-24** (`ir::stereo`; additive) Add site and ligand operations to the stereo
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
  Set::compact(&self, graph: &GraphCompaction) -> (Self, Compaction<Id>)
  ```

  StereoBonds restoration translates its bond sites and the atom references in
  its ligands. These are crate-private set methods. S2j exposes
  `StereoLigandPosition` on the public editor views. Site and ligand replacement
  leave the other factor and payload unchanged; frame-preserving permutation
  remains a separate operation.

  Implemented the recorded interfaces through `FixedVarBirelationSet`, preserving
  copy-on-write ownership. Tests cover atom and bond sites, all ligand kinds,
  duplicates and empty frames, shared anchors, independent site/ligand changes,
  unchanged configuration and constraints, incidence, and restoration after
  entity removal and topology compaction. All 56 focused tests, strict graph-IR
  Clippy, nightly formatting, and `git diff --check` pass.
- **S1d — reverted; replaced by S2i.** Storage integration and mutable-view
  restructuring proceed together, without temporary view factories or adapters.
  Completed S1a–S1c are retained.

### S2 — Uniform mutable views and removal of callback mutation

Implement replacements first, migrate their callers, then remove all callback
mutation. Existing callbacks may remain temporarily while their callers are
migrated; add no new callback API or compatibility layer. S2a is complete, but its
checked molecule-level constraint mutation surface is removed in S2i1.
The earlier S2b attempt remains reverted. S2a's completed record is retained;
S2b–S2m below replace the former unimplemented S2b–S2f worklist. Each subitem
states semantics, interface/nomenclature, and verification. S2b covers the Rust
unit-error rename and Python join error mapping; the JoinError enum is withdrawn. S2f is
cancelled. S2c's bounded operation fixes and S2g's frame-consumer policy are approved.

- **S2a — completed 2026-09-24; API superseded, removal in S2i1** (`ir::molecule::constraints`, `ir` re-export; additive, green)
  Add the checked molecule-level ConstraintsViewMut and expose it through
  Molecule::constraints_mut. Keep the editor's existing &mut Constraints.
  [dep: S0a]

  ```text
  Molecule::constraints(&self) -> &Constraints
  Molecule::constraints_mut(&mut self) -> ConstraintsViewMut<'_>
  MoleculeEditor::constraints_mut(&mut self) -> &mut Constraints

  ConstraintsViewMut::push(&mut self, Constraint) -> Result<(), MoleculeIntegrityError>
  ConstraintsViewMut::extend(&mut self, Constraints) -> Result<(), MoleculeIntegrityError>
  ConstraintsViewMut::replace(&mut self, Constraints) -> Result<(), MoleculeIntegrityError>
  ConstraintsViewMut::remove_at(&mut self, usize) -> Constraint
  ConstraintsViewMut::clear(&mut self)
  ConstraintsViewMut::{len, is_empty, as_slice, iter}(&self)
  ```

  The view has one private &mut Molecule field, no public constructor, and no
  mutable collection escape. Define the view and Molecule::constraints_mut in a
  private molecule child module and re-export ConstraintsViewMut through ir.
  This permits access to private molecule storage without adding an unchecked
  mutation gateway or widening its fields. Reuse the existing incoming-reference
  and stereo integrity checks; do not validate unrelated storage or clone the receiver.
  Check every incoming entry before extend or replace writes anything. Preserve
  order and duplicates. Remove/clear are unchecked; invalid removal positions
  panic as on Constraints. Replace the test-only raw Molecule::constraints_mut
  access with this public contract, moving deliberately invalid setup into an
  editor or constructor input as appropriate. Keep try_modify_constraints until
  its bindings and remaining callers have migrated.

  Test valid/invalid references, nested constraints, stereo-frame failures,
  unchanged collections after a late extend failure, duplicate/order preservation,
  replacement, removal, and clearing through public entry points.

  Implemented the listed public surface with a private owner borrow and no public
  constructor or mutable escape. Reused reference checks and made the existing
  stereo-wrapper check crate-visible for the sibling view module. Checks inspect
  incoming constraints only; writes move entries without copying the molecule or
  allocating a second input buffer. Migrated the test-only raw accessor callers;
  editor access and callback APIs remain unchanged pending their planned stages.
  Updated the nomenclature guide for the new molecule-level view.

  Verification: 800 focused molecule, canonicalization, and DSL unit cases pass,
  including 33 new accessor/mutation cases. Three publication-agreement properties
  pass at 128 cases each, covering success and unchanged state on rejection.
  Graph-IR strict Clippy (all targets, proptest enabled), rustdoc with warnings
  denied, nightly formatting, and diff review pass. No workspace or MSRV gate
  was run for this subitem.

- **S2b — completed 2026-09-26: NoJoinError in Rust and Python** (`umol-graph-ir::ir`,
  `umol-graph-ir-macros`, `umol-py::{lattice,error}`, extension registration and
  Python package exports; breaking name and Python behavior, red→green). [dep: S0a]

  **Semantics.** Python join returns a form on success and raises NoJoinError
  when Rust join returns NoJoinError. Absence of a join is an operation failure, not
  the ordinary bottom result represented by None in meet. Leave Python meet and
  normalize unchanged. Rename Rust's NoJoin unit error to NoJoinError in its
  definition, re-export, trait/implementation signatures, derive macro, callers,
  tests, and documentation. Keep its display message and behavior unchanged;
  no compatibility alias, JoinError enum, or InvalidTerm variant.

  **Interfaces and nomenclature.** Rust exposes the unit struct NoJoinError;
  join returns Result<Self, NoJoinError> and widen_with returns
  Result<bool, NoJoinError>. Export umol.NoJoinError, derived from
  Exception, alongside the existing binding exceptions. Preserve the Rust
  NoJoinError display message. In the shared impl_py_lattice macro, change:

  ```diff
  -fn join(&self, py: Python<'_>, other: &Self) -> PyResult<Option<Self>>;
  +fn join(&self, py: Python<'_>, other: &Self) -> PyResult<Self>;
  ```

  The macro uses its existing concrete Python type for other. Map the existing
  Rust error directly to NoJoinError; retain the existing successful conversion.
  Register the exception in the extension and Python package exports. Update
  join documentation and any affected annotations; do not add a second join
  method or reuse ContradictionError for the distinct NoJoinError failure.

  **Verification.** Python tests cover a successful join, a pair of distinct
  constraint kinds whose Rust join returns NoJoinError, the exception's public import
  and message, and an incompatible meet still returning None. Update existing
  tests that expect None from failed join. Rebuild with the repository Python
  3.13 environment and run focused binding tests. Run focused Rust join tests
  and compile the renamed public error through its consumers; no workspace or
  MSRV gate.

  Implemented the Rust unit-error rename throughout graph-IR, its Lattice derive,
  tests, and the nomenclature guide. Python join now returns PyResult<Self> and
  raises the exported NoJoinError with Rust's display message. No lattice
  algorithms, normalization behavior, or meet behavior changed.

  Verification: 162 focused Rust tests and 175 focused Python tests pass after
  rebuilding the extension with Python 3.13. Coverage includes successful joins,
  different-key/scope failures, public exception exports, meet bottom, and
  writable join results from read-only entity attributes. Nightly formatting,
  diff checks, and full diff review pass. No workspace or MSRV gate was run.

- **S2c — completed 2026-09-26: Coset action handling and domain simplification** (`ir::stereo`;
  behavior changes, green). [dep: S0a]

  **Semantics.** Preserve unrestricted coset assignment and existing lattice
  behavior. Add no blanket coset-range checks to normalization, meet, join,
  matches, or is_compatible. Invalid indices need not produce meaningful
  algebraic results. Retain checked lookup when an operation interprets an index.
  Apply these four bounded changes:

  1. StereoCoset::apply checks that the supplied permutation is allowed by the
     known kind before dispatching on the coset variant. Return None for an
     incompatible action, including symbolic and undetermined cosets. Do not
     scan stored indices. swap and mirror generate compatible actions themselves.
  2. compose_term checks each explicit Apply action before composing it. Return
     None for wrong-degree or disallowed actions, including nested ones. This
     covers directly constructed terms as well as terms built through methods.
  3. canon_coset returns the original term unchanged when compose_term returns
     None. A variable domain becomes unrestricted only if it covers exactly
     0..kind.count(); retain other domains without error or trimming. Plain
     literals/sets and existing empty-domain contradictions are unchanged.
  4. coset_apply_permutation maps failed symbolic evaluation to None through its
     existing return type, rather than Some(Undetermined). Its symmetry::reexpress
     caller already handles None. Existing checked literal/set branches remain.

  **Interfaces and nomenclature.** Only the private composition signature changes:

  ```rust
  fn compose_term(term: &StereoTerm, kind: StereoKind)
      -> Option<(&StereoTerm, Permutation)>;
  ```

  The returned term borrows the input. Retain canon_coset's existing
  Result<StereoCoset, Contradiction> and coset_apply_permutation's
  Option<StereoCoset>. No new normalization error variant or blanket range
  failure is added. Existing checked reindex failures during literal action
  evaluation keep their current behavior. All public normalize, apply, swap,
  mirror, reframe_by, and coset_for signatures remain unchanged. No new public
  validator, lookup type, checked setter, or binding interface is introduced.

  **Verification.** Test incompatible supplied actions on literal, symbolic, and
  undetermined cosets; nested incompatible actions preserve the original term
  during normalization without panic. Test that complete variable domains become
  unrestricted and same-size domains containing an invalid index remain intact.
  Check existing literal/set transport failures and symbolic evaluation returning
  None rather than Undetermined. Preserve valid normalization idempotence and
  frame-transport laws. Do not require new lattice rejection or correct symmetry
  classification for malformed indices. Run focused stereo tests and relevant
  properties; no workspace or MSRV gate at this subitem.

  Implemented the four changes in ir::stereo. Incompatible explicit actions
  return None from compose_term; canon_coset preserves the original term.
  Complete-domain folding now checks both domain size and its largest index,
  without allocating or scanning a second set. Symbolic evaluation failures in
  coset_apply_permutation return None. Public signatures and lattice algorithms
  remain unchanged; no new normalization error or range-validation pass was added.

  Verification: 1,541 stereo-related unit tests and 62 stereo-related property-suite
  tests pass (PROPTEST_CASES=128), including valid normalization and
  frame-transport laws. New exact cases cover incompatible actions, nested
  preservation, incomplete domains, and existing evaluation failures. Nightly
  formatting, diff checks, and full scope review pass. No workspace, Python,
  or MSRV gate was run for this Rust-only subitem.

- **S2d — completed 2026-09-26: Role-only incidence and count-aware consumers** (`ir::incidence`,
  `ir::canonicalize`, `ir::substructure`; breaking, red→green).
  [dep: S0a, S2c]

  **Semantics.** Incidence construction records connectivity and chemical roles,
  not electron contributions. It must not inspect electron-vector length.
  Canonicalization still distinguishes contributions attached to different
  atoms; read them from the owning forms at color/key construction. Check
  counts.len() == atom_count before associating literal counts with atoms.
  Return the existing Contradiction on failure and check once before candidate
  enumeration, not inside each candidate. Undetermined counts remain allowed.
  Apply this to molecule and reaction-span sides.

  **Interfaces and nomenclature.** The complete replacement edge-label enum is:

  ```rust
  pub enum Incidence {
      BondEndpoint,
      DativeDonor,
      DativeAcceptor,
      AromaticAtom,
      MulticenterAtom,
      NoncovalentEndpoint,
      StereoSite,
      StereoLigand(StereoLigandKind),
  }
  Molecule::incidence_graph(&self, IncidenceLevel) -> IncidenceGraph
  ReactionSpan::incidence_graph(&self, IncidenceLevel) -> IncidenceGraph
  ```

  Remove AromaticParticipant, MulticenterParticipant, their count-bearing Span
  variants, and electron_span extraction. Retain the existing initial_color_keys
  and reaction_span_entity_keys Result boundaries and return shapes. Their
  aromatic/multicenter occurrence colors read entity attributes; role-only labels
  must not remove chemical information from canonical comparison. Do not replace
  counts with positional edge colors, which would constrain allowed atom permutations.

  **Concrete before/after mapping.** Today an aromatic/multicenter incidence edge
  stores a copied electron contribution, or an EntitySpan of contributions for a
  reaction span. After this change the edge stores only AromaticAtom or
  MulticenterAtom. ReactionSpan continues to own its existing EntitySpan<Form>
  attributes; no entity spans, lhs/rhs values, or change categories are removed
  from that storage.

  In initial_color_keys, use the incidence edge's endpoints and
  IncidenceGraph::entity to identify the owning aromatic system or multicenter
  bond and its incident atom. Find that atom's position in the owning set's
  stored atom sequence. Read the contribution at that position from the form's
  electrons field, or use Undetermined when the field is undetermined. Construct
  the existing contribution-bearing InitialColorKey::Incidence in canonicalize,
  rather than putting the contribution back into IncidenceGraph. Position is
  only a lookup index; it is not included in the color.

  reaction_span_entity_keys does the same lookup against the span's owning set.
  Match its existing EntitySpan<Form>: Unchanged/Added/Removed supplies one
  contribution; Modified supplies lhs and rhs contributions at the same stored
  atom position. Encode the existing span category and corresponding value(s)
  directly in the canonical key. Retain the existing entity colors that also
  encode entity-span categories. No new *Span incidence variant, persistent
  contribution map, public lookup type, or change to ReactionSpan's shape is
  needed; the existing color vectors hold the derived keys.

  Check each literal vector against its owner's atom count once before these
  lookups and candidate enumeration, including both Modified values and all
  other present span values. Every count-dependent canonicalization route must
  reach that check, including constitution_candidate's electron_occurrence_fields
  path; do not rely on zip truncation or repeat the check for each candidate.
  Preserve that path's association of atom images with contributions when sorting
  occurrences. Final molecule/span comparison continues reading the owning forms.

  Incidence matching uses role labels and reads aromatic/multicenter contributions
  from owning forms during edge comparison. Check each literal vector's length
  before looking up the incident atom's contribution; mismatch on either side
  rejects the tentative pairing. Undetermined pattern counts impose no contribution
  restriction. The existing verify_overlays comparison of transported attributes
  remains authoritative for completed candidates, including length disagreement
  through its Option path. Both matching algorithms must return identical
  valid-input match sets. Public signatures and SubstructureMatchError are unchanged.

  **Verification at S2d.** Exact incidence labels and unchanged topology for
  constructible molecules/spans, canonical keys distinguishing atom/count
  association, reframing invariance, and matching-algorithm agreement on valid
  count differences. Implement the count checks here; public-API tests for
  short/long counts and malformed matching candidates belong to S2h, after
  construction admits those inputs. Do not bypass construction or widen
  visibility for tests. Run existing
  focused canonicalization/matching benchmarks before and after: removing early
  pruning can affect search work. Record a decision-relevant regression, not
  precision measurements of invalid-input behavior. Symmetry receives no new
  count checks or error return.

  Implemented the role-only enum and both incidence builders. The entity-key
  pass checks each selected literal count vector before incidence lookup and
  candidate enumeration; modified spans check both sides. Contribution-bearing
  keys retain the stored atom/count association and existing span normalization.
  Matching now prunes count-incompatible pairings during search by reading owning
  forms, and still checks both matched forms before transport through verify_overlays.
  Operation signatures, matching errors, symmetry checks, and normalization are unchanged.
  Malformed aggregate tests remain in S2h, as specified above.

  **Measured tradeoff.** The initial role-only implementation removed count-based
  edge pruning and slowed incidence matching when counts distinguished otherwise
  interchangeable atoms. The follow-up below restores pruning in the matcher
  while retaining role-only incidence storage. Both strategies remain explicit
  caller choices; no public API was added.

  Criterion, VF2, 20 samples, 1 s warm-up, 2 s requested measurement (corpus
  measurements automatically ran longer); times rounded:

  | Matching case | Strategy | Before | After |
  | --- | --- | ---: | ---: |
  | Six-carbon ring, counts [2,1,0,1,0,1], self-match | Incidence | 21.4 µs | 32.3 µs |
  | Same ring | GraphAndOverlays | 20.5 µs | 20.2 µs |
  | Three-carbon multicenter, counts [2,0,1], self-match | Incidence | 1.59 µs | 4.78 µs |
  | Same multicenter | GraphAndOverlays | 3.81 µs | 3.78 µs |
  | Existing phenol pattern over OpenSMILES corpus | Incidence | 727 ms | 707 ms |
  | Same corpus/pattern | GraphAndOverlays | 236 ms | 239 ms |

  The new bounded cases live in the existing umol-graph substructure benchmark;
  the corpus pattern has no count-bearing overlay and cannot expose this loss
  of pruning. The corpus changes were not statistically distinguishable.
  Incidence construction also carries smaller labels and no copied contributions.
  The existing canonicalize benchmark was run before/after for ordinary_naphthalene
  and overlay_heavy: complete canonicalization measured 48.5→22.6 µs and
  194→137 µs; constitution incidence construction measured 1.59→0.525 µs and
  1.64→0.562 µs. The initial canonicalization run overlapped another crate's
  compilation, so those differences do **not** establish a speedup. Matching
  measurements above ran after compilation without concurrent benchmark work.

  Verification: 39 incidence-related unit tests, 286 canonicalization tests,
  and 55 matching tests pass. The new exact cases cover span categories,
  nonuniform contributions, reordered atom lists, and undetermined counts.
  Canonicalization and matching property filters pass 24 and 2 tests respectively
  (PROPTEST_CASES=128), including reframing/remapping laws and algorithm agreement.
  Rustdoc with warnings denied, nightly formatting, diff checks, and full scope
  review pass. No workspace, Python, or MSRV gate was run.

  **Early filtering — tested and adopted 2026-09-26.** The edge comparator reads
  the contribution from each edge's owning form: identify the overlay and atom,
  check the literal vector length against the atom list, then locate that atom's
  position. A length mismatch rejects the tentative edge match. Undetermined
  pattern counts impose no contribution restriction. No incidence payload,
  contribution cache, heap allocation, or public API was added; verify_overlays
  still checks completed candidates. This restores early rejection of atom
  permutations with incompatible contributions without changing incidence storage.

  Repeated the bounded VF2 cases with uniform and undetermined controls, using
  the same Criterion settings and sequential, isolated measurements:

  | Incidence matching case | Late check only | Adopted early filter |
  | --- | ---: | ---: |
  | Aromatic ring, nonuniform counts | 31.61 µs | 20.97 µs |
  | Aromatic ring, uniform counts | 32.59 µs | 34.34 µs |
  | Aromatic ring, undetermined counts | 31.80 µs | 32.49 µs |
  | Multicenter, nonuniform counts | 4.71 µs | 1.61 µs |
  | Multicenter, uniform counts | 5.17 µs | 5.37 µs |
  | Multicenter, undetermined counts | 4.75 µs | 5.02 µs |

  The nonuniform cases improve by 34% and 66%; controls pay 2–6% for lookups
  that cannot prune. GraphAndOverlays control timings remain within about 1.3%.
  These examples support matcher-side early filtering, with a small measured
  cost when contributions do not distinguish atoms. They do not establish a
  workload-wide gain. The tested comparator is now implemented in
  visit_substructure_matches_incidence, with no new helper or cache. The expanded
  benchmark and exact permutation tests remain. All 63 matching unit tests pass
  with either implementation, including all six graph algorithms; early filtering
  also passes both matching properties with PROPTEST_CASES=128.
  Short/long-vector aggregate tests remain scheduled for S2h,
  when public construction admits those values. Results here survive deletion of
  the temporary prototype and logs in scratch/s2d-early-filter.

- **S2e — moved to S2h.** Existing electron-use boundaries and getter behavior
  require verification after aggregate construction admits the relevant forms.
  They do not require a separate implementation subitem.

- **S2f — cancelled.** No depiction/CoordGen error-handling, error-enum, or
  private-signature changes are included in this work.

- **S2g — completed. Frame-dependent attribute consumers** (`ir::canonicalize`,
  `ir::stereo` transport, `umol-graph::ops::validate::stereo`;
  behavior changes and enum extension, red→green). [dep: S2c]

  **Semantics.** Uniform attribute assignment includes changing the stereo kind,
  symbolic actions, and entity constraints. Check ligand-count/kind agreement
  before canonicalization applies a kind-sized permutation to ligand storage;
  return its existing Contradiction. Retain FrameTransport's existing degree,
  allowed-action, and constraint-position checks and Option failure. Validate
  both positions of a topicity pair before querying the symmetry group. A
  wrong-degree ligand-symmetry assertion must not become a satisfied negative
  assertion merely because group membership returns false.

  **Interfaces and nomenclature.** Canonicalization and transport signatures
  stay unchanged. Public graph_symmetry, stereo_atom_symmetry,
  stereo_bond_symmetry, and their result queries stay infallible; no promise of
  correct chirality classification on invalid coset indices is added.

  **Approved failure behavior.** For a kind/ligand-count
  mismatch in observable_descriptor, return its existing None before creating
  kind-sized swaps. This keeps the symmetry API infallible and prevents that
  specific indexing/expect hazard; it makes no claim that the resulting symmetry
  is meaningful for the malformed frame. Existing preconditions unrelated to
  these attribute changes are not redesigned here.

  For validation of a non-undetermined topicity assertion, add:

  ```rust
  StereoConformanceContradiction::TopicityPositionOutOfRange {
      pair: StereoLigandPair,
      degree: usize,
  }
  ```

  Return Solution::Contradictory with that payload before querying either
  position. For a wrong-degree ligand-symmetry assertion, use the existing
  LigandSymmetryViolation { asserted } rather than allowing a negative assertion
  to pass. Keep validate_symmetry's existing Solution return and the public
  validate Result<Solution<...>, StereoConformanceError> signature. Do not invent
  a derived topicity for an invalid pair, add errors to symmetry itself, or add
  a separate validator type.

  **Verification.** Kind-sized permutations versus shorter/longer ligand frames,
  both topicity positions, wrong-degree positive and negative symmetry
  assertions, and existing valid frame-transport/canonicalization properties.
  Tests needing aggregate states currently rejected by construction run in S2h;
  standalone form/action cases and valid aggregate cases run here. Tests for
  malformed coset indices must not demand correct symmetry output.

  **Implemented and verified — 2026-09-26.** Both kind-sized canonicalization
  permutation paths reject shorter/longer frames with Contradiction.
  observable_descriptor returns None on kind/ligand-count disagreement.
  Conformance validation checks both topicity positions and rejects wrong-degree
  ligand-symmetry assertions, including negative assertions. Existing transport
  checks and all operation signatures remain unchanged; the only new public
  symbol is TopicityPositionOutOfRange.

  Passed 297 canonicalization, 19 symmetry, 350 stereo, and 31 stereo-conformance
  unit tests; 22 canonicalization and 31 frame property tests at 128 cases each;
  rustdoc for graph-ir/graph with warnings denied, nightly formatting, and
  git diff --check. Malformed aggregate tests remain in S2h.

- **S2h — Attribute agreement leaves aggregate integrity**
  (`ir::{molecule::integrity,stereo::integrity,reaction::integrity,reaction_span}`,
  constructor/publication callers, electron-use/getter documentation, consumer
  regression tests, and guides; breaking accepted-input change, red→green).
  [dep: S2a, S2c, S2d, S2g]

  **Semantics.** Public construction, editor publication/probe, reaction
  construction, and both reaction-span projections use the same remaining
  integrity contract. Do not allow a value through one path and reject it
  elsewhere solely for the removed attribute agreement. Construction stores
  supplied forms faithfully: no normalization, count padding/truncation,
  coset repair, or stereo erasure.

  Remove literal electron-count/atom-count agreement and coset-index bounds.
  To permit the complete mutable configuration and entity-constraint borrows,
  also move per-entity kind/site, kind/ligand-count, term-permutation-degree,
  and constraint-position/degree requirements to the consumers in S2c/S2g.
  These are attribute agreements, not changes to the stored topology. Keep
  valid entity references, unique/parallel relation rules, ligand incidence,
  duplicate-ligand rejection, and the storage maximum on ligand-frame size.
  Top-level constraint integrity checks remain at construction/publication.
  Top-level writes use the editor; Molecule exposes read-only constraints.

  **Interfaces and nomenclature.** Constructor and publication signatures stay
  unchanged. Remove ElectronCountLengthMismatch and StereoCosetOutOfRange from
  MoleculeIntegrityError/ReactionIntegrityError and their conversion/mapping
  arms. Remove corresponding shared error/check code when no retained consumer
  uses it. Other stereo error variants still used by the separate top-level
  constraint checks remain; do not remove them mechanically by name.
  Reaction Add/Remove and ModifyField old/new payloads follow the same policy.
  Retain reaction reference/removal-incidence checks and the reaction-span rule
  requiring matching determined kinds for a preserved stereo entity.

  **Existing electron-use boundaries and unchanged getters (moved from S2e).**
  Retain ElectronCountsForm::reframe_by returning None on a
  literal-vector/action length mismatch. Retain existing checks in
  MatchingInput::from_system, DelocalizationPlan::derive, and
  AromaticityPerceiver::derive, with ElectronCountMismatch, no applicable plan,
  and AromaticSystemFailure respectively. Do not add a whole-molecule precheck.

  No signature changes to those operations or these four getters:

  ```rust
  AtomView::aromatic_valence(&self) -> NumForm
  AtomView::multicenter_valence(&self) -> NumForm
  AromaticSystemView::electron_count(&self) -> NumForm
  MulticenterBondView::electron_count(&self) -> NumForm
  ```

  Preserve checked per-position access: missing cells give Undetermined; extra
  cells are ignored by per-atom getters and included by totals. No length check
  is added to the getters. The older permute methods are not used by the active
  transport paths; their separate replacement/signature question is not a
  prerequisite for these count-integrity changes. Do not silently redesign them
  here. Remove documentation promising later constructor rejection.

  **Verification.** Constructor/editor publication/Reaction/ReactionSpan agreement on
  newly admissible count and coset payloads; faithful storage and valid-input serialization;
  continued rejection of structural defects. Exercise the first-use rejections
  through public APIs after construction, and the unchanged getter behavior.
  Malformed inputs require panic freedom, not correct results or faithful text
  round-tripping.
  This subitem owns the deferred aggregate cases from S2d/S2g and all
  former S2e coverage. Construct inputs through the public constructors now
  admitting them, using the current editor publication API until S6 introduces
  probe/finish. In the same subitem, verify:

  - Incidence construction remains infallible for short/long counts in molecules
    and spans; count-dependent canonicalization rejects them and matching rejects
    malformed candidates through its existing failure path (S2d).
  - Canonicalization, transport, and stereo validation handle the newly admitted
    frame/constraint disagreements as specified in S2g, without demanding correct
    symmetry classification for malformed configurations.
    Cover the construction-dependent expect calls around complete_candidate in
    canonicalize_full_with_options and reframe_by in reaction_span_canonical_candidate:
    newly admitted constraint/frame disagreements must reach the existing
    Contradiction path before those infallible search callbacks.
  - Standalone electron transport still returns None on a degree mismatch;
    the existing chemistry boundaries retain ElectronCountMismatch, no applicable
    plan, and AromaticSystemFailure. Getter cases use [1, 2] and [1, 2, 0, 5]
    for a three-atom entity and assert the accepted per-atom and total behavior.

  S2h is not complete when checks are merely deleted: these consumer regressions
  and the retained structural-rejection cases must pass in the same subitem.
  Update data-types.md and integrity.md with the actual retained checks now;
  do not leave normative guides contradicting the implementation until S9.

  **Implemented and verified — 2026-09-26.** Removed attribute-agreement checks
  and the two obsolete error variants. Structural and top-level constraint
  checks remain. Constructors and current editor publication preserve supplied
  forms; reaction payloads and both span projections follow the same contract.
  Canonicalization checks stereo frame transport before its infallible search
  callbacks. Existing electron-use failures and getter behavior are unchanged.
  The normative guides and affected Rust/Python tests now reflect this contract.

  Passed 7,087 graph-IR unit tests (3 ignored), 19 public construction/consumer
  regressions, 417 property tests at 128 cases (1 ignored), 1,959 graph tests,
  259 TableIR raise tests, and 259 Python tests
  against a rebuilt Python 3.13 extension. Strict Clippy and rustdoc passed for
  graph-ir/graph/io; nightly formatting and git diff --check passed.

- **S2i — Uniform mutable attribute access and typed editor storage**
  (`ir::{view,molecule,molecule::editor}`, ir exports; group; S2i1–S2i5 complete).
  [dep: S1a, S1b, S1c, S2h]

  S2i1–S2i5 are complete. S2i5 separates the mutable view types; S2j adds getters
  and structural methods.

- **S2i1 — Retire the public checked constraint view**
  (`ir::{molecule,view}`, graph-IR callers/tests; breaking, red→green).
  [dep: S2a, S2h]

  Remove public Molecule::constraints_mut, ConstraintsViewMut, its re-export,
  and the private molecule::constraints module. Migrate the remaining checked
  accessor callers in DSL/canonicalization tests to the existing editor's
  constraints_mut and publication path. Preserve constructor/publication tests
  for invalid top-level references and stereo constraints; retire only tests
  of the removed checked-view API. Python uses try_modify_constraints instead
  and its caller migration remains S2l.

  Add crate-private Molecule::constraints_mut(&mut self) -> &mut Constraints in the
  owning molecule module. S2i3 uses it for the editor delegate; no public
  top-level constraint borrow is added. Remove the checked-view description
  from the nomenclature guide. Keep try_modify_constraints until S2k/S2l migrate
  its consumers, then remove it in S2m. The editor delegates to the private
  accessor rather than accessing Molecule fields directly.

  Verify editor top-level mutation and checked publication, retained read-only
  Molecule access, and absence of the checked view/type export.

  **Implemented and verified — 2026-09-26.** Removed ConstraintsViewMut,
  its exports/module, and its API-specific tests. Test setup uses the editor
  and checked publication. Molecule retains a crate-private constraints_mut accessor,
  used by the existing checked callback; editor delegation follows in S2i3.
  Constructor rejection cases remain, and editor publication tests cover invalid
  references, stereo kinds, frame sizes, permutation degrees, and positions.
  Removed the checked-view description from the nomenclature guide.

  Passed 7,059 graph-IR unit tests (3 ignored) and all 5 publication properties
  at 128 cases. Strict Clippy, rustdoc with warnings denied, nightly formatting,
  and git diff --check passed.

- **S2i2 — Shared mutable-view types and Molecule access** (`ir::{view,molecule}`, exports; breaking, green at S2i4). [dep: S1a, S1b, S1c, S2i1]

  **Implemented record; type sharing is superseded by S2i5.** The signatures
  below describe the completed implementation, not the target view design.

  **Semantics.** Molecule and editor expose the complete form of each entity
  through a mutable borrow, including its constraints. Attribute assignment is
  immediate, unchecked, and identical in capability for all eight entity kinds.
  No journal, cloning for recovery, setter-time checks, validation on drop,
  per-field mutable-accessor family, or runtime permission flag. The const
  parameter restricts structural methods only; it never restricts attributes.
  Existing copy-on-write detachment remains an ordinary storage implementation
  detail. Atom/ligand/site structural replacement remains editor-only.

  **Interfaces and nomenclature.** Use one *ViewMut family for all eight kinds,
  with const EDITOR: bool = false. Molecule accessors return <false>; editor
  accessors return <true>. No distinct *EditorViewMut structures remain.

  ```rust
  pub struct AromaticSystemViewMut<'a, const EDITOR: bool = false> {
      id: AromaticSystemId,
      set: &'a mut AromaticSystems,
  }

  impl<const EDITOR: bool> AromaticSystemViewMut<'_, EDITOR> {
      pub fn id(&self) -> AromaticSystemId;
      pub fn attributes(&self) -> &AromaticSystemForm;
      pub fn attributes_mut(&mut self) -> &mut AromaticSystemForm;
      pub fn constraints(&self) -> &AromaticSystemConstraintsForm;
      // Other read methods from the local-getter inventory.
  }

  impl AromaticSystemViewMut<'_, true> {
      pub fn replace_atoms(&mut self, atoms: &[AtomId]);
      // Remaining structural methods are listed in S2j.
  }

  // Molecule:
  pub fn aromatic_system_mut(
      &mut self, id: AromaticSystemId,
  ) -> AromaticSystemViewMut<'_, false>;
  // MoleculeEditor:
  pub fn aromatic_system_mut(
      &mut self, id: AromaticSystemId,
  ) -> AromaticSystemViewMut<'_, true>;

  // Assignment on either view:
  view.attributes_mut().charge = value; // an atom view
  view.attributes_mut().constraints.set(constraint);
  ```

  Apply the same const parameter, common accessors, and entry-point return types
  to Atom, Bond, DativeBond, MulticenterBond, NoncovalentBond, StereoAtom, and
  StereoBond, using each entity's concrete form/constraint/id types. Atom and
  localized bond views have no structural methods in this scope. Add the missing
  molecule StereoAtomViewMut and StereoBondViewMut access through
  Molecule::stereo_atom_mut/stereo_bond_mut. Existing *_mut entry points keep
  their entity-id arguments and panic for invalid ids.

  **Structures.** All fields are private; the listed borrows carry lifetime 'a.
  Each row defines one structure shared by both const specializations:

  | View | Fields |
  | --- | --- |
  | AtomViewMut | id: AtomId; attributes: &mut AtomForm |
  | BondViewMut | id: BondId; atoms: [AtomId; 2]; attributes: &mut BondForm |
  | DativeBondViewMut | id: DativeBondId; set: &mut DativeBonds |
  | AromaticSystemViewMut | id: AromaticSystemId; set: &mut AromaticSystems |
  | MulticenterBondViewMut | id: MulticenterBondId; set: &mut MulticenterBonds |
  | NoncovalentBondViewMut | id: NoncovalentBondId; set: &mut NoncovalentBonds |
  | StereoAtomViewMut | id: StereoAtomId; set: &mut StereoAtoms |
  | StereoBondViewMut | id: StereoBondId; set: &mut StereoBonds |

  No public constructors. Each view has a crate-private new taking exactly its
  listed fields. Molecule owns the calls that construct views, establishes that
  the entity id exists, and selects the const specialization. Public Molecule
  accessors return <false>; crate-private access supplies the existing <true> views to
  editor, batch, and undo execution. Editor accessors delegate to that private
  access rather than borrowing Molecule fields to construct views themselves.
  Spell that crate-private access as the following Molecule methods, each taking
  &mut self and the listed id and returning the listed view. Public *_mut
  methods call the <false> specialization; editor access calls <true>.

  | Private method, each with const EDITOR: bool | Id | Return |
  | --- | --- | --- |
  | atom_view_mut | AtomId | AtomViewMut<'_, EDITOR> |
  | bond_view_mut | BondId | BondViewMut<'_, EDITOR> |
  | dative_bond_view_mut | DativeBondId | DativeBondViewMut<'_, EDITOR> |
  | aromatic_system_view_mut | AromaticSystemId | AromaticSystemViewMut<'_, EDITOR> |
  | multicenter_bond_view_mut | MulticenterBondId | MulticenterBondViewMut<'_, EDITOR> |
  | noncovalent_bond_view_mut | NoncovalentBondId | NoncovalentBondViewMut<'_, EDITOR> |
  | stereo_atom_view_mut | StereoAtomId | StereoAtomViewMut<'_, EDITOR> |
  | stereo_bond_view_mut | StereoBondId | StereoBondViewMut<'_, EDITOR> |

  These methods select and construct the existing view, not a second mutation
  API. They have the same invalid-id panic as the public accessors. Their bodies
  live in the molecule module and use pub(crate) visibility. Top-level mutable
  constraint access uses the crate-private Molecule::constraints_mut() -> &mut Constraints introduced in S2i1.
  Internal attribute-only writes can use the ordinary <false> views. No new view
  family, public const-selection parameter, or false-to-true conversion is added.
  Export the shared family through ir. id copies the id; attributes,
  attributes_mut, constraints, and structural getters borrow for the method
  call, not 'a. No public conversion changes false to true. EDITOR affects
  available methods only, with no stored flag or runtime branch.

  Add focused cases for uniform attribute access and const-specialized capability.

  **Completed — 2026-09-26; group verification recorded in S2i4.** All eight
  shared mutable-view types have the specified private fields, crate-private
  constructors, and common accessors. Molecule's public accessors return the
  false specialization; its private view constructors support either value.
  Removed the separate mutable editor-view types. Added 16 attribute/constraint
  cases across both specializations and eight invalid-id cases.

- **S2i3 — Editor storage and view delegation** (`ir::molecule::editor`; breaking, green at S2i4). [dep: S2i2]

  Store the editor's draft in its Molecule field, using the six existing typed
  overlay sets rather than *Store wrappers. This private storage restructuring
  precedes S3/S4's shared mutation access; it does not change edit's public
  ownership contract, which changes in S6. An editor overlay view needs only its
  typed entity id and a mutable borrow of that set, with short accessor borrows.
  Molecule views must expose
  attributes without exposing set mutation. Localized atom/bond forms borrow
  their existing attribute tables; localized-bond rewiring is not added.
  Forward existing addition/removal/reframe operations to typed sets. Preserve
  intermediate-state access and current publication behavior until S6.

  Migrate editor construction/read/write paths and test typed-set incidence and
  preservation of current publication behavior. No new public names.

  **Completed — 2026-09-26; group verification recorded in S2i4.** MoleculeEditor
  holds a private Molecule draft and its existing correspondence. The three
  storage-wrapper types and their Arc conversions are removed. All eight mutable
  accessors delegate to Molecule's crate-private accessors with EDITOR=true;
  constraints_mut delegates to Molecule's crate-private accessor. Overlay additions,
  removals, compaction, and restoration use the typed sets. Equivalence checks
  retain frame transport and normalized comparison. Snapshot and build check
  integrity before returning a Molecule; edit's ownership contract is unchanged.

  Six new tests cover live incidence after addition, dense removal, and restoration
  for every overlay kind, plus snapshot independence and whole-molecule publication.
  The wrapper-materialization tests are removed with their types; existing
  publication/error and frame-equivalence tests remain.

- **S2i4 — Slice additions and caller migration** (editor callers across Rust/Python; breaking, red→green). [dep: S2i3]

  Direct editor additions use slices for variable-length atom/donor/ligand
  lists, consistently with the typed sets and the replacement methods in S2j.
  Change these existing arguments; all other arguments and return types stay
  unchanged:

  | MoleculeEditor method | Current argument | Planned argument |
  | --- | --- | --- |
  | add_dative_bond | donors: Vec<AtomId> | donors: &[AtomId] |
  | add_aromatic_system | atoms: Vec<AtomId> | atoms: &[AtomId] |
  | add_multicenter_bond | atoms: Vec<AtomId> | atoms: &[AtomId] |
  | add_stereo_atom | ligands: Vec<StereoLigand> | ligands: &[StereoLigand] |
  | add_stereo_bond | ligands: Vec<StereoLigand> | ligands: &[StereoLigand] |

  Fixed-size lists retain arrays; individual ids and owned forms remain values.
  The packed relation storage copies list elements into its own buffers; it
  cannot adopt a separate per-row vector. AtomId-to-NodeId conversion may still
  require a temporary vector. Migrate every caller of these five signatures in
  this subitem, borrowing existing vectors or locally collecting generated
  sequences as needed. MoleculeBuilder keeps its current public iterator inputs;
  only its calls into the editor change. Edits/Delta payloads retain owned vectors.

  **Compilation boundary.** S2i4 closes this group green. In this subitem migrate every
  caller affected by private fields, accessor borrows, const-specialized return
  types, removed editor-view types, typed storage, and slice addition arguments.
  This includes editor application/undo and existing callback implementations,
  Rust/DSL/graph/IO consumers, Python binding internals, tests, examples, and
  benchmarks. Supply the existing read-access equivalents needed by those
  callers here; S2j adds the remaining local-getter inventory and structural
  methods. The structural signatures shown above describe the resulting API;
  their implementation and tests belong to S2j. Do not add compatibility fields,
  view aliases, or adapters to postpone these caller migrations. Existing
  callback entry points remain until S2k/S2l migrate their uses and S2m removes
  them; update their internals as needed to compile against the new views.

  **Verification.** Ordinary and whole-form assignment for every entity kind,
  constraints mutation, short/long electron vectors, invalid cosets, and access
  after successive attribute changes through the editor. Verify no deferred work
  runs on drop, typed-set incidence remains synchronized after existing add/remove
  operations, and attribute assignment is available on both const specializations.
  S2j tests structural mutation on editor views and its absence from molecule views.
  Review editor mutable accessors to confirm they obtain views from Molecule;
  storage borrows and copy-on-write mutation remain inside Molecule and views.

  **Completed — 2026-09-26.** The five additions take slices. Rust/Python callers,
  benchmarks, and tests use the shared mutable-view accessors. MoleculeBuilder's
  iterator inputs and Edits' owned vectors are unchanged. Bond, dative, aromatic,
  and multicenter mutable views supply atom_ids for existing callers. All eight
  editor mutable accessors delegate to Molecule's view constructors.

  Tests exercise public Molecule/editor attribute and constraint assignment,
  short and long electron vectors, invalid cosets, successive access, editor drop,
  and typed-set incidence after addition/removal/restoration. Removal panic tests
  expect the delegated graph-core diagnostic. The Python reaction-SMILES case
  C[S@]C>> expects stereo-resolution contradiction under the S2g/S2h contract.

  Verification:

  - Workspace all-target compilation and strict Clippy pass.
  - Unit suites: graph 1,959 and IO 4,137 pass. Graph-IR's full run passed 7,064
    cases, with 24 stale panic-message assertions; after correcting those
    assertions, all 61 removal cases pass (three unrelated tests remain ignored).
  - With PROPTEST_CASES=128: edit 20, molecule 79, and project 157 properties pass.
  - Python 3.13 extension rebuild succeeds. Pytest passed 1,573 cases, with one
    stale error assertion and two skips; after correcting that assertion, all
    ten reaction-SMILES error cases pass.
  - Graph-IR rustdoc with warnings denied and doc-tests pass. Nightly formatting,
    full diff review, and git diff --check pass.

- **S2i5 — Separate molecule and editor mutable views**
  (`ir::{view,molecule,molecule::editor}`, exports and explicit view-type callers;
  breaking, red→green). [dep: S2i4]

  **Semantics.** Preserve uniform unchecked attribute/constraint assignment and
  typed editor storage. Replace the const-generic mutable family with distinct
  molecule/editor types. This subitem changes type selection only; structural
  methods and new local getters remain in S2j. Immutable views remain separate.
  No recovery copy, journal, runtime flag, or conversion between view families.

  **Structures and names.** Each row defines two public structures with identical
  private fields. Borrows carry 'a. Each has a crate-private new taking those
  fields, no public constructor, and no Clone/Copy implementation.

  | Molecule type | Editor type | Fields |
  | --- | --- | --- |
  | AtomViewMut<'a> | AtomEditorViewMut<'a> | id: AtomId; attributes: &'a mut AtomForm |
  | BondViewMut<'a> | BondEditorViewMut<'a> | id: BondId; atoms: [AtomId; 2]; attributes: &'a mut BondForm |
  | DativeBondViewMut<'a> | DativeBondEditorViewMut<'a> | dative_bonds: &'a mut DativeBonds; id: DativeBondId |
  | AromaticSystemViewMut<'a> | AromaticSystemEditorViewMut<'a> | aromatic_systems: &'a mut AromaticSystems; id: AromaticSystemId |
  | MulticenterBondViewMut<'a> | MulticenterBondEditorViewMut<'a> | multicenter_bonds: &'a mut MulticenterBonds; id: MulticenterBondId |
  | NoncovalentBondViewMut<'a> | NoncovalentBondEditorViewMut<'a> | noncovalent_bonds: &'a mut NoncovalentBonds; id: NoncovalentBondId |
  | StereoAtomViewMut<'a> | StereoAtomEditorViewMut<'a> | stereo_atoms: &'a mut StereoAtoms; id: StereoAtomId |
  | StereoBondViewMut<'a> | StereoBondEditorViewMut<'a> | stereo_bonds: &'a mut StereoBonds; id: StereoBondId |

  **Methods.** Both types in each row retain identical id(&self) -> EntityId,
  attributes(&self) -> &EntityForm, attributes_mut(&mut self) -> &mut EntityForm,
  and constraints(&self) -> &EntityConstraintsForm. Returned form borrows last
  for the accessor borrow. Preserve the existing atom_ids() methods on both
  types: [AtomId; 2] for Bond; impl ExactSizeIterator<Item = AtomId> + '_ for
  DativeBond, AromaticSystem, and MulticenterBond.

  Molecule's eight public *_mut(&mut self, id) methods return *ViewMut<'_> without
  a const argument. The synonymous editor methods return *EditorViewMut<'_>.
  Molecule's existing crate-private atom_view_mut, bond_view_mut, dative_bond_view_mut,
  aromatic_system_view_mut, multicenter_bond_view_mut, noncovalent_bond_view_mut,
  stereo_atom_view_mut, and stereo_bond_view_mut keep their names and id arguments,
  lose the const parameter, and return the corresponding *EditorViewMut<'_>.
  Public Molecule accessors construct molecule views directly; editor accessors
  delegate to those private methods. Batch/undo execution uses the same private
  access. No public structural route is added to Molecule.

  Preserve invalid-id panics. Export all eight separate mutable editor types;
  migrate explicit generic view annotations and turbofish calls in this subitem.
  No aliases, forwarding view wrappers, new traits, or implementation-generating
  macros. Ordinary implementation duplication is acceptable.

  Within each entity kind, order declarations and their implementations as
  *Views, *View, *EditorView, *ViewMut, *EditorViewMut; supporting functions follow,
  then tests. Stereo-atom and stereo-bond declarations form separate groups;
  the existing shared stereo-query macro follows those groups. Use the entity
  collection name for each mutable overlay's stored borrow, and acceptor for
  the acceptor field in both readonly dative views. Public acceptor_id() is unchanged.

  **Verification.** Preserve the completed attribute/constraint and incidence
  tests while changing their view types. Verify both families expose the same
  attribute behavior, typed editor storage remains in use, and no const-generic
  mutable-view references remain. Run the affected graph-IR tests and compile
  downstream callers, including Python with its prescribed environment. Finish
  green without depending on S2j's additions.

  **Completed — 2026-09-26.** All eight mutable molecule/editor pairs are separate
  public types with matching existing methods and unchanged stored borrows.
  Removed the const parameter, rewired public/crate-private accessors, exported the
  editor mutable types, and updated the nomenclature guide. Existing assignment
  tests name both public return types explicitly; their assertions are unchanged.
  No immutable-view API changes or S2j getter/structural additions are included.

  Graph-IR unit tests: 7,088 passed, three ignored. Workspace all-target strict
  Clippy (including Python bindings under Python 3.13), graph-IR rustdoc with
  warnings denied, nightly formatting, and git diff --check pass. Reviewed the
  full code diff and confirmed matching method bodies for every mutable pair;
  no const-generic mutable-view references remain in code. Field names and
  declaration ordering are aligned across the seven entity-view modules; the
  ordering correction passes all 418 view tests. Function-body comparison
  confirms only field renames, with test contents unchanged.

  **Atom accessor alignment — completed 2026-09-26.** AtomView and
  AtomEditorView expose private id/attributes through id() and attributes();
  all four atom views use #[inline] on matching accessors. AtomView retains its
  molecule, id, and separately supplied attribute borrow in that order.
  Rust callers use the accessors; Python properties and mutation are unchanged.
  Graph-IR unit tests pass (7,088 passed, three ignored), as do all 71 Python
  atom tests after rebuilding the extension. Workspace all-target compilation
  includes the graph-IR/graph/IO property targets and Python depiction feature;
  strict Clippy passes for the four affected crates with those targets/features.
  Graph-IR rustdoc with warnings denied, nightly formatting, and git diff --check
  pass. Caller changes preserve existing operations and test assertions.

  **Bond accessor alignment — completed 2026-09-26.** The four bond views expose
  private ids and attribute borrows through matching accessors with #[inline].
  BondEditorView also exposes its private endpoints through atom_ids(). Rust
  callers are migrated; Python properties and setters are unchanged. Graph-IR
  unit tests pass (7,088 passed, three ignored), as do all 78 Python bond tests
  after rebuilding. Feature-enabled workspace all-target compilation and
  graph-IR rustdoc with warnings denied pass. Strict Clippy for the four affected
  crates, nightly formatting, and git diff --check pass. Reviewed the full diff;
  callers preserve their operations and assertions.

  **Dative accessor alignment — completed 2026-09-26.** The four dative views
  expose matching id, attribute, donor-id, acceptor-id, and combined-atom accessors.
  Callers are migrated; invalid-id panics and Python properties are preserved.
  The 34 added cases cover row selection, readonly borrow lifetimes, donor order,
  empty donor lists, and write-through mutation. Graph-IR unit tests pass
  (7,122 passed, three ignored), along with 37 Python dative tests after rebuilding.
  Feature-enabled workspace all-target compilation, affected-crate strict Clippy,
  graph-IR rustdoc, nightly formatting, and diff checks pass. Full diff reviewed.

  **Aromatic accessor alignment — completed 2026-09-26.** The four aromatic
  views use the structures and matching accessors listed above. Rust callers
  are migrated; Python properties, unrestricted attribute mutation, and invalid-id
  panics are preserved. Twenty added cases cover system selection, stored atom
  order, readonly borrow lifetimes, and write-through mutation. Graph-IR unit
  tests pass (7,142 passed, three ignored), along with 43 Python aromatic tests
  after rebuilding. Feature-enabled workspace all-target compilation, affected-crate
  strict Clippy, graph-IR rustdoc, nightly formatting, and diff checks pass.
  Full diff reviewed.

  **Multicenter accessor alignment — completed 2026-09-26.** The plural view
  holds only molecule; molecule-backed views use raw_multicenter_bonds(). The four
  entity views expose matching id, attribute, and atom-id accessors. Callers are
  migrated; Python properties, unrestricted attribute mutation, and invalid-id
  panics are preserved. Twenty added cases cover bond selection, stored atom
  order, readonly borrow lifetimes, and write-through mutation. Graph-IR unit
  tests pass (7,162 passed, three ignored), along with 43 Python multicenter tests
  after rebuilding. Feature-enabled workspace all-target compilation, affected-crate
  strict Clippy, graph-IR rustdoc, nightly formatting, and diff checks pass.
  Full diff reviewed.

  **Noncovalent accessor alignment — completed 2026-09-26.** The plural view
  holds only molecule; molecule-backed views use raw_noncovalent_bonds(). All four
  entity views expose matching id/attribute accessors and atom_ids() -> [AtomId; 2].
  Callers are migrated; Python properties, unrestricted attribute mutation, and
  invalid-id panics are preserved. Twenty added cases cover bond selection,
  readonly attribute-borrow lifetimes, endpoint order, and write-through mutation.
  Graph-IR unit tests pass (7,182 passed, three ignored), along with 40 Python
  noncovalent tests after rebuilding. Feature-enabled workspace all-target
  compilation, affected-crate strict Clippy, graph-IR rustdoc, nightly formatting,
  and diff checks pass. Full diff reviewed.

  **Stereo view alignment completed — 2026-09-26.** Both stereo families use the
  structures above. Basic id/attribute/site/frame access is aligned across all
  four views, with borrowed ligand frames and migrated callers. Shared macros
  cover repeated implementations; no additional S2j getters or structural
  mutation methods were added. Ten added cases cover borrowed access,
  write-through mutation, stored frame order, and invalid editor ids. Graph-IR
  unit tests pass (7,192 passed, three ignored), as do 62 stereo properties and
  117 Python stereo/constraint tests after rebuilding. Feature-enabled all-target
  strict Clippy, graph-IR rustdoc, nightly formatting, and diff checks pass.
  The full scoped diff was reviewed.

- **S2j — completed 2026-09-26 — Local getters and editor structural mutation**
  (`ir::{id,view}`, typed sets, ligand_ids callers and bindings; additive
  structural methods plus breaking getter return change, red→green). [dep: S2i5]

  **Status.** Complete for all eight entity families. Matching local getters,
  borrowed stereo ligand frames, and editor-only structural mutations are implemented.

  **Atom and localized-bond scope completed — 2026-09-26.** AtomEditorView and
  both mutable atom views expose element(), isotope_mass(), charge(),
  implicit_hydrogens(), lone_pairs(), and unpaired_electrons(). BondEditorView and
  both mutable bond views expose order(), charge(), and unpaired_electrons().
  All match the readonly molecule getters, with 'a on readonly editor borrows
  and method borrows on mutable views. Inline annotations match. No structural
  mutation methods are added for these entities. Six new cases cover readonly
  borrow lifetimes and reads after successive attribute writes. The focused atom
  and bond view suites pass 77 and 23 cases, respectively; graph-IR library/test
  Clippy with warnings denied, nightly formatting, and diff checks pass. The full
  scoped diff is reviewed.

  **Dative scope completed — 2026-09-26.** DativeBondEditorView and both mutable
  dative views expose order(), donor_count(), and atom_count(), matching the
  immutable molecule view. Readonly order() retains 'a; mutable access borrows
  through &self. DativeBondEditorViewMut delegates replace_donors,
  replace_acceptor, replace_donor, insert_donor, and remove_donor to DativeBonds.
  AtomPosition is defined in ir::id and re-exported through ir, with the traits,
  index(), and From<usize> convention of StereoLigandPosition. No molecule-view
  structural mutation is exposed. The focused dative view suite passes 79 cases,
  including ordering, empty/end positions, preserved forms and constraints,
  incidence updates, and publication rejection of invalid atom lists. Graph-IR
  all-target Clippy with proptest enabled and warnings denied, strict rustdoc,
  nightly formatting, and the scoped diff check pass. The scoped diff is reviewed.

  **Aromatic scope completed — 2026-09-26.** AromaticSystemEditorView and both
  mutable aromatic views expose electrons(), charge(), unpaired_electrons(),
  electron_count(), and atom_count(), matching AromaticSystemView. Readonly
  references retain 'a; mutable-view references borrow through &self.
  AromaticSystemEditorViewMut delegates replace_atoms, replace_atom, insert_atom,
  and remove_atom to AromaticSystems. Attributes and constraints are preserved;
  atom-list and attribute writes can be performed in either order. The focused
  aromatic-view suite passes 72 cases, including getter lifetimes, counts, stored
  order, bounds, incidence, preservation of other systems, and publication checks.
  Strict graph-IR library/test Clippy, rustdoc, nightly formatting, and diff checks
  pass. The full scoped diff is reviewed.

  **Multicenter scope completed — 2026-09-26.** MulticenterBondEditorView and both
  mutable multicenter views expose electrons(), charge(), unpaired_electrons(),
  electron_count(), and atom_count(), matching MulticenterBondView. Readonly
  references retain 'a; mutable-view references borrow through &self.
  MulticenterBondEditorViewMut delegates replace_atoms, replace_atom, insert_atom,
  and remove_atom to MulticenterBonds, preserving attributes and constraints.
  The focused multicenter-view suite passes 66 cases covering getter lifetimes,
  counts, ordering, bounds, incidence, preserved attributes and other bonds,
  either order of atom-list/attribute writes, and publication checks. Overlapping
  bonds remain accepted; identical atom sets are rejected. Strict graph-IR
  library/test Clippy, rustdoc, nightly formatting, and diff checks pass. The full
  scoped diff is reviewed.

  **Noncovalent scope completed — 2026-09-26.** NoncovalentBondEditorView and both
  mutable noncovalent views expose kind(), matching NoncovalentBondView. Readonly
  references retain 'a; mutable-view references borrow through &self.
  NoncovalentBondEditorViewMut delegates replace_atoms and replace_atom to
  NoncovalentBonds, preserving attributes and constraints. The focused view suite
  passes 44 cases covering getter lifetimes and reads after writes, endpoint
  order, both positions, bounds, incidence, preserved attributes and other bonds,
  and publication rejection of missing atoms, self-loops, and parallel bonds.
  Strict graph-IR library/test Clippy, rustdoc, nightly formatting, and diff
  checks pass. The full scoped diff is reviewed.

  **Stereo scope completed — 2026-09-26.** Both stereo families expose the
  approved local getters on readonly editor views and both mutable views.
  ligand_position(&self, id: AtomId) uses id consistently, including the readonly
  molecule views. Shared macros implement the getters and editor-only replace_site,
  replace_ligands, replace_ligand, insert_ligand, and remove_ligand delegates.
  Getters borrow stored attributes and ligand slices; structural writes preserve
  configuration and constraints and update the owning sets' incidence.
  The stereo-view suite passes 151 cases, including actual/virtual ligand filters,
  borrow lifetimes, undetermined-configuration panics, position bounds, empty
  lists, stored order, site/ligand preservation, either attribute-write order,
  incidence, and structural publication rejection. The S2j closeout passes all
  693 view tests, 17 stereo property tests, graph-IR all-target Clippy with proptest
  enabled and warnings denied, strict rustdoc, nightly formatting, and diff checks.
  The rebuilt Python 3.13 extension passes all 173 stereo/molecule tests.
  The full scoped diff is reviewed.

  **Readonly stereo editor storage — implemented.** StereoAtomEditorView stores
  stereo_atoms: &'a StereoAtoms followed by id: StereoAtomId. StereoBondEditorView
  stores stereo_bonds: &'a StereoBonds followed by id: StereoBondId. All fields are
  private; each crate-private new takes the set borrow and id in that order.
  attributes(), site_id(), and ligand_ids() read the owning set. Readonly
  attribute/frame borrows retain 'a. Public site/ligands field callers are migrated.
  Editor lookup checks the id before constructing the view and preserves the
  immediate invalid-id panic. Both mutable families also hold their owning sets.

  **Semantics.** Read-only getters borrow existing data; filters are lazy and
  preserve stored order. ligand_ids borrows the stored ordered slice, with no
  allocation or reconstruction through ligand views. Structural
  replacements preserve attributes/constraints and the other relation factor.
  Single replacement preserves length; insertion/removal preserve unaffected
  order. Writes delegate through the typed set to graph-core and update incidence
  immediately. Do not add mutable participant slices or combined replacements.

  **Interfaces and nomenclature.** Add AtomPosition(pub u32), index(self) -> usize,
  and From<usize>, matching StereoLigandPosition's traits/conversion convention.
  AtomPosition names positions in non-ligand atom lists, including donors;
  StereoLigandPosition names ligand positions. ParticipantPosition stays inside
  graph-core delegation. The exact public editor methods, all &mut self -> (), are:

  | Entity | Methods |
  | --- | --- |
  | AromaticSystem / MulticenterBond | replace_atoms(&[AtomId]); replace_atom(AtomPosition, AtomId); insert_atom(AtomPosition, AtomId); remove_atom(AtomPosition) |
  | NoncovalentBond | replace_atoms([AtomId; 2]); replace_atom(AtomPosition, AtomId) |
  | DativeBond | replace_donors(&[AtomId]); replace_acceptor(AtomId); replace_donor(AtomPosition, AtomId); insert_donor(AtomPosition, AtomId); remove_donor(AtomPosition) |
  | StereoAtom | replace_ligands(&[StereoLigand]); replace_site(AtomId); replace_ligand(StereoLigandPosition, StereoLigand); insert_ligand(StereoLigandPosition, StereoLigand); remove_ligand(StereoLigandPosition) |
  | StereoBond | The same ligand methods; replace_site(BondId) |

  replacement/removal require position < length; insertion allows position ==
  length. Violations panic. Fixed-array length is enforced by the type. Atom
  ids are not validated at these sharp editor operations; publication checks
  remaining structural integrity. No localized-bond endpoint rewiring.
  Implement the exact local-getter inventory above on the corresponding views,
  without substituting raw values or Option results for existing molecule-view
  return types. Count getters stay
  infallible with the accepted malformed-input behavior from S2h.

  Implemented on both stereo families, including immutable editor views and both
  separate mutable families:

  ```rust
  pub fn ligand_ids(&self) -> &[StereoLigand];
  ```

  Immutable views return &'a [StereoLigand] from the owning set; mutable views
  borrow through the owning set and return &[StereoLigand] tied to &self. Both
  prevent structural mutation while the slice is in use. No mutable frame slice is
  exposed. Migrated callers include: canonicalization, correspondence,
  reactions, reaction spans, reaction integrity, substructure matching, DSL,
  Python bindings, and tests/properties. Read-only consumers use the slice or
  iter().copied(); use to_vec only where a consumer needs owned storage or must
  sort/change the frame. Do not blanket-copy at every call site to restore the
  old return type. Python's existing ligand-property collection is built directly
  from the borrowed slice without an intermediate Rust Vec; its public return
  shape is unchanged by this getter correction. Complete these signature-driven
  caller changes here rather than leaving the build broken until S2k/S2l.
  S2j finishes green, including Python binding compilation and affected tests;
  it does not depend on the later callback or property-behavior migrations.

  **Verification.** Borrowed frames preserve exact stored order and include actual
  and virtual ligands. Retain caller-level canonicalization, matching, reaction,
  DSL, and Python outcomes after migration; review every introduced to_vec for
  an actual ownership requirement. Ordering, boundaries, unaffected fields/factors, immediate
  incidence, arbitrary attribute/sequence assignment order, and retained
  structural publication rejection. No test should expect publication failure
  solely for an electron-length or coset-index disagreement.

- **S2k — Rust callback caller migration** (`dsl::molecule`,
  `umol-graph::ops::transform::delocalize_charge`, remaining Rust callers;
  breaking rewire, red→green). [dep: S2i, S2j]

  Signature and field-access migrations required to compile are already complete
  in S2i/S2j. S2k1 covers DSL conversion; S2k2 covers charge delocalization and
  remaining Rust callers. Each finishes green.

- **S2k1 — completed 2026-09-26 — In-place DSL defaults conversion** (`dsl::{molecule,atom,bond,
  aromatic,multicenter,stereo}`; internal rewire, green).
  [dep: S2i, S2j]

  **Semantics.** MoleculeDsl conversion mutates attributes through Molecule's
  mutable views. FromIr clones the source once; IntoIr moves self.molecule and
  mutates it directly. Preserve defaults, metadata, ordering, and faithful
  conversion. Neither path uses an editor or performs an integrity check.
  Replace the eight entity conversion calls in each direction as follows:

  | Entity | FromIr | IntoIr |
  | --- | --- | --- |
  | Atom | atom::lower_atom | atom::raise_atom |
  | Localized bond | bond::lower_bond | bond::raise_bond |
  | Aromatic system | aromatic::lower_aromatic_system | aromatic::raise_aromatic_system |
  | Multicenter bond | multicenter::lower_multicenter_bond | multicenter::raise_multicenter_bond |
  | Dative bond, noncovalent bond, stereo atom, stereo bond | No mutation | No mutation |

  The first four kinds use the existing in-place functions through
  attributes_mut(). The other four conversions leave their forms unchanged, so
  remove those molecule-level passes entirely. No temporary entity DSL wrapper,
  per-form clone, or mem::take is needed. Keep the standalone entity DSL
  conversions; reduce StereoAtomDsl::from_ir's two consecutive clones to one.

  **Interfaces and nomenclature.** Keep MoleculeDsl's FromIr/IntoIr signatures.
  Make lower_aromatic_system, raise_aromatic_system, lower_multicenter_bond, and
  raise_multicenter_bond pub(crate), matching the existing atom/bond functions.
  Their signatures remain (&mut Form, &Defaults) -> (). Add no functions or
  public API. These conversions do not depend on the editor lifecycle or its
  S6 migration; IntoIr needs neither an editor storage copy nor drop(self).

  **Verification.** DSL defaults and roundtrips, source preservation, metadata
  behavior, exact entity and frame order, unchanged forms for the four no-op
  kinds, and standalone stereo conversion. Preserve property laws and malformed
  attribute cases; do not introduce correctness requirements for manipulated
  inputs.

  **Verification results — 2026-09-26.** The four lower/raise pairs now mutate
  molecule attributes directly; the eight no-op passes are removed.
  StereoAtomDsl::from_ir clones once. DSL unit tests pass (2,418 cases), as do
  all 11 molecule serialization/defaults properties. Added independent expected
  values cover both conversion directions across all entity kinds and source
  preservation. Nightly formatting and diff checks pass; full diff reviewed.
  Strict Clippy is blocked only by the now-unused crate-private
  modify_aromatic_systems, modify_multicenter_bonds, modify_stereo_atoms, and
  modify_stereo_bonds methods. Their removal remains in S2m; no lint allowances
  were added.

- **S2k2 — completed 2026-09-26 — Charge delocalization and remaining Rust callers**
  (`umol-graph::ops::transform::delocalize_charge`, Rust tests/fixtures;
  internal rewire, green). [dep: S2k1]

  **Semantics.** Charge delocalization applies its plans directly to the borrowed
  molecule through mutable entity views, preserving its planning checks and
  Infallible result. Attribute writes require no editor, publication, or
  integrity check.

  **Interfaces and nomenclature.** Keep DelocalizeCharge's current transformation
  signature. Its rename to
  ChargeDelocalizer is tracked separately in
  [166](166-molecule-ops-2026-07-27.md#charge-delocalization-transformer-name).
  DelocalizationPlan::apply takes &mut Molecule and uses the approved mutable
  view names/accessors from S2i. No mem::take of a borrowed caller's molecule, temporary
  empty receiver, new modify helper, or transitional compatibility callback.
  Migrate ordinary calls and callback method names supplied through macros.
  Move remaining try_modify_constraints callers to editor mutation and
  publication. S2i1 already migrated the public checked-view callers; retain
  reference/frame integrity cases at constructor/editor publication
  boundaries, including deliberately invalid inputs and order preservation.

  **Verification.** Charge delocalization outcomes and tests/fixtures that used
  callbacks. Preserve
  property laws and invalid-input cases; do not make fixtures valid merely to
  avoid the newly explicit first-use failure behavior.

  **Verification results — 2026-09-26.** Delocalization applies its plans through
  molecule mutable views. Remaining Rust consumers use mutable entity
  attributes or editor constraints; the callback methods' own tests remain until
  S2m removes those methods. Inputs, assertions, and property laws are unchanged.
  All 808 selected transformation/resolver/validator tests, 14 resolver/isotope
  properties, seven focused canonicalization/editor-publication tests, and 19
  integrity integration tests pass. All-target Clippy completes with only the
  four unused callback warnings recorded in S2k1; warnings-denied Clippy remains
  blocked by those methods pending S2m. Nightly formatting and diff checks pass;
  full diff reviewed. Python callers and the transformer rename are unchanged.

- **S2l — Python assignment for every entity kind**
  (`umol-py::{atom,bond,dative,aromatic,multicenter,noncovalent,stereo,molecule,constraint}`;
  binding rewire, red→green). [dep: S2a, S2i, S2j]

  Binding compilation against S2i/S2j is already restored. This subitem changes
  assignment behavior and removes binding callbacks, then finishes green.

  **Semantics.** Property assignment and nested entity-constraint mutation write
  through to the stored Rust form in molecule-backed and standalone access.
  Aromatic systems, multicenter bonds, stereo atoms, and stereo bonds have the
  same assignment capability as the other entities. Assigning an incompatible
  electron vector or coset succeeds; the later operation requiring agreement
  owns any failure. No copy/edit/publication cycle or detached mutable form.

  **Interfaces and nomenclature.** Preserve existing Python property names and
  setter syntax, including whole attributes assignment. Delegate to the whole
  mutable-form borrow from S2i; do not create Python set_electrons/set_coset or
  try_modify methods. Entity-level constraints mutate their stored collection.
  Remove the Molecule.constraints setter and mutating methods on the
  molecule-backed top-level ConstraintsView; nested entries must not provide a
  write path into that collection. Reads remain live and read-only. Migrate
  top-level mutation callers to editor/batch operations and publication. Keep
  standalone Constraints and entity-level constraint setters mutable. Invalid
  entity/collection indices remain IndexError.
  Remove each molecule-backed with_mut callback as its callers migrate.
  Prepared transactions and consumption remain in S5/S6. No interactive Python
  transaction or new editor-view binding family is introduced.

  **Read-only constraint entries.** Preserve Constraint.Atom, Constraint.And,
  and the other existing variants, constructors, tuple fields, equality, and
  pattern matching. Constraint privately stores either an owned Rust value or
  a molecule owner, captured counter, collection position, and child-index path.
  ConstraintsView stores either the molecule/counter or a parent Constraint for
  an And/Or child sequence. ConstraintIter retains that sequence, a cursor, and
  its end position. Each step yields one constraint; compositions remain whole
  entries. Child access is explicit and read-only. No constraint tree is copied
  to create a molecule accessor or iterator. Leaf payloads use the existing
  immutable Python value forms, converted when requested. Constraint.copy()
  copies the selected subtree into an independent value; standalone Constraints
  remains mutable and returns independent values.

  Add the already specified Molecule counter and check_access/advance_counter/
  view_counter methods for these accessors. combine_from advances the counter
  immediately before Rust execution; ordinary entity assignment does not.
  Collections, entries, child sequences, and exhausted iterators reject stale
  access with InvalidatedViewError. Conversion into edits/deltas propagates this
  error. S5d2 extends this mechanism to the remaining views and operation
  boundaries; S6d adds owner consumption.

  **Verification.** Python 3.13 assignment tests for all eight kinds, whole-form
  replacement, nested constraints, short/long electron vectors, out-of-range
  cosets, and continued access to backing storage. Check that molecule-level
  constraints cannot be written through the property, collection, or nested
  accessors, while standalone/entity-level constraints remain mutable. Preserve
  invalid top-level reference/frame rejection at editor publication.
  Build the extension before tests and exercise public Python access rather
  than test-only Rust mutation paths.

  **Implemented and verified — 2026-09-26.** All entity property and nested
  constraint writes use direct Rust mutable views; Python with_mut/try_modify
  callers are removed. Molecule constraints use the read-only entry/sequence
  design above. The rebuilt Python 3.13 extension passes 748 focused Python
  tests; all 1,639 binding unit tests pass (two ignored), including editor
  publication rejection of invalid constraint references and stereo frames.
  Variant matching, nested sequences, explicit copying, invalidation, and
  ordinary setter usability are covered. All-target Clippy adds no warnings;
  the four unused Rust callback methods still block warnings-denied lint until
  S2m removes them. Nightly formatting and diff checks pass; full diff reviewed.

- **S2m — Remove callback mutation and close the stage**
  (`ir::molecule`, remaining Rust/Python callers, public documentation;
  breaking removal, red→green). [dep: S2a, S2k, S2l]

  **Semantics.** Attribute assignment has one direct mutable-borrow path for
  every entity kind. Remove the old copy/check/replace callback mechanisms;
  do not retain private equivalents or compatibility shims.

  **Interfaces and nomenclature.** Remove Molecule::modify_atoms,
  modify_bonds, modify_dative_bonds, modify_aromatic_systems,
  modify_multicenter_bonds, modify_noncovalent_bonds, modify_stereo_atoms,
  modify_stereo_bonds, try_modify_aromatic_system, try_modify_aromatic_systems,
  try_modify_multicenter_bond, try_modify_multicenter_bonds,
  try_modify_stereo_atom, try_modify_stereo_atoms, try_modify_stereo_bond,
  try_modify_stereo_bonds, try_modify_constraints, and try_modify_checked.
  This includes private callback methods. S2i1 already removed the public
  checked constraints_mut accessor, ConstraintsViewMut, its module and export,
  and its dedicated API tests. Retain the private mutable accessor used by the
  editor, Molecule::constraints, and MoleculeEditor::constraints_mut. Direct entity
  attribute access remains; only the separate editor mutable views expose
  structural methods.
  Update guides, rustdoc, examples, and bindings to describe current behavior,
  without discussion-doc citations in source.

  **Verification.** Search both call syntax and macro-supplied identifiers to
  establish complete removal. Run the affected graph-IR/graph suites and
  explicit property suites, Python 3.13 build/tests, strict affected lint and
  rustdoc, nightly formatting, and full diff review. Compare the existing S0
  mutation measurements and the S2d matching/canonicalization measurements only
  where the change can affect cost. Full-workspace and Rust 1.87 gates remain
  at S9b, not after every subitem. This stage closes only after all listed
  migrations are green.

  **Implemented and verified — 2026-09-26.** Removed all 18 listed callback
  methods, including private implementations, and their obsolete tests. No
  production callers remained. Existing mutable-view tests cover every entity
  kind, electron-count mismatches, coset values, entity constraints, and invalid
  ids; stereo kind-assignment cases are retained through those views. Editor
  publication tests retain invalid molecule-constraint rejection. Molecule
  rustdoc and the integrity guide describe the direct entity-mutation and
  editor-publication boundaries.

  Graph-IR and graph tests pass with both proptest features and graph conformance
  enabled: 10,712 passed, eight ignored, including doctests. All 1,639 Python
  binding tests pass (two ignored); the rebuilt Python 3.13 extension passes
  1,589 Python tests (two skipped). All-target Clippy and rustdoc pass with
  warnings denied for graph-IR, graph, and Python bindings. Nightly formatting,
  removal searches, diff checks, and full diff review pass. No remaining
  production execution path changes in S2m, so the S0 and S2d measurements were
  not rerun. Full-workspace and Rust 1.87 gates remain S9b.

### S3 — Batch mutation and reaction vocabulary

S3a–S3e form one breaking enum migration; the green boundary is S3e, including
Rust execution, DSL conversion, and Python exhaustive matches. Do not make an
intermediate subitem compile by adding wildcard rejection, dropped variants,
or unimplemented execution branches.

Approved DSL verbs match the Rust method names, using hyphens in EDN. The
Edit DSL uses the existing entity key to select the entity kind:

| Entity keys | Rust methods | DSL verbs |
| --- | --- | --- |
| :aromatic-system, :multicenter-bond, :noncovalent-bond | replace_atoms | :replace-atoms |
| :dative-bond | replace_donors, replace_acceptor | :replace-donors, :replace-acceptor |
| :stereo-atom, :stereo-bond | replace_site, replace_ligands | :replace-site, :replace-ligands |

The payload follows the existing Edit :modify convention:
[handle {:expect old :update new}]. Both values are supplied explicitly. Execution
resolves their handles and compares old with the stored value.

Atom and donor lists are EDN vectors, preserving order. Noncovalent atom vectors
contain exactly two entries. A site or acceptor is one atom/bond handle. Ligand
lists use the existing ligand encoding. Both :expect and :update carry these
values directly; they do not use attribute-form strings.

```edn
{:aromatic-system
 {:replace-atoms [0 {:expect [0 1 2] :update [0 1 3]}]}}
```

Retain the Edit DSL's existing handle resolution. Single-position replacement,
insertion, and removal need no additional Edit variants or DSL operations.
The reaction DSL retains its addition, removal, and modification vocabulary.

- **S3a — implemented** (`ir::edit`, `dsl::edit`; group; breaking, green at S3e)
  [dep: S2i]

- **S3a1 — Structural Edit/Undo variants and dative factors** (`ir::edit`; breaking, green at S3e). [dep: S2i]

  **Semantics.** Add batch representations of structural replacement, which
  Edit currently cannot express. Each edit replaces one complete named component;
  it preserves the entity id, other components, attributes, and constraints.
  Single-position changes use whole-list replacement; no additional single-position
  variants or independent-batch composition API are added. Execution, old-state
  comparison, and undo application belong to S3b.

  **Interfaces and nomenclature.** Add these nine Edit/Undo pairs. Every Edit
  has exactly id, old, and new fields; every Undo has id and the saved field
  listed below. For entity kind X, Edit id is XHandle and Undo id is XId.

  | Edit variant | old/new type | Undo variant | Saved field and type |
  | --- | --- | --- | --- |
  | ReplaceAromaticSystemAtoms | Vec<AtomHandle> | RestoreAromaticSystemAtoms | atoms: Vec<AtomId> |
  | ReplaceMulticenterBondAtoms | Vec<AtomHandle> | RestoreMulticenterBondAtoms | atoms: Vec<AtomId> |
  | ReplaceNoncovalentBondAtoms | [AtomHandle; 2] | RestoreNoncovalentBondAtoms | atoms: [AtomId; 2] |
  | ReplaceDativeBondDonors | Vec<AtomHandle> | RestoreDativeBondDonors | donors: Vec<AtomId> |
  | ReplaceDativeBondAcceptor | AtomHandle | RestoreDativeBondAcceptor | acceptor: AtomId |
  | ReplaceStereoAtomSite | AtomHandle | RestoreStereoAtomSite | site: AtomId |
  | ReplaceStereoAtomLigands | Vec<(AtomHandle, StereoLigandKind)> | RestoreStereoAtomLigands | ligands: Vec<StereoLigand> |
  | ReplaceStereoBondSite | BondHandle | RestoreStereoBondSite | site: BondId |
  | ReplaceStereoBondLigands | Vec<(AtomHandle, StereoLigandKind)> | RestoreStereoBondLigands | ligands: Vec<StereoLigand> |

  For example, the additions to the existing enums are:

  ```rust
  // Edit
  ReplaceAromaticSystemAtoms {
      id: AromaticSystemHandle,
      old: Vec<AtomHandle>,
      new: Vec<AtomHandle>,
  }
  // Undo
  RestoreAromaticSystemAtoms {
      id: AromaticSystemId,
      atoms: Vec<AtomId>,
  }
  ```

  Undo saves the actual previous component with resolved ids, not batch-local
  handles. Existing Edits::push accepts the new variants; no new convenience
  method family is introduced.

  Correct the existing dative construction/removal shapes in the same stage:

  ```rust
  // Edits: replace the combined atoms argument with donors and acceptor.
  pub fn add_dative_bond(
      &mut self, donors: Vec<AtomHandle>, acceptor: AtomHandle,
      attributes: DativeBondForm,
  ) -> DativeBondHandle;
  pub fn add_dative_bonds(
      &mut self,
      bonds: impl IntoIterator<Item = (Vec<AtomHandle>, AtomHandle, DativeBondForm)>,
  ) -> Vec<DativeBondHandle>;
  pub fn remove_dative_bonds(
      &mut self,
      removes: Vec<(DativeBondHandle, Vec<AtomHandle>, AtomHandle, DativeBondForm)>,
  );

  // Edit: replace atoms with the two named factors.
  AddDativeBond {
      donors: Vec<AtomHandle>, acceptor: AtomHandle, attributes: DativeBondForm,
  }
  RemoveDativeBonds {
      removes: Vec<(DativeBondHandle, Vec<AtomHandle>, AtomHandle, DativeBondForm)>,
  }
  ```

  AddedDativeBond and RemovedDativeBond likewise replace atoms: Vec<AtomId>
  with donors: Vec<AtomId> and acceptor: AtomId; id and attributes stay unchanged.
  Preserve donor order and existing comparison semantics. An acceptor is now
  structurally present; this adds no donor-count or chemistry validation.
  S3a3 adapts DSL conversion/rendering to these fields without changing the
  existing dative syntax. Update Edit/Edits rustdoc to the construction and batch-namespace
  contract in “Edits — one accumulated sequence”, including update_* and raw
  push/collection semantics; remove the transaction-only and detached-journal
  descriptions.

  Verify all nine payload shapes, exact list order, Id/New handle namespaces,
  fixed arity, and dative donor/acceptor construction. Execution cases belong to S3b.

  **Implemented — 2026-09-26.** Edit and Undo carry the nine named `Replace*`
  variants; dative construction, removal, and saved entries carry donors and
  acceptor separately. Edits keeps its append-only handle accounting. Local
  construction tests and rustdoc reflect the new payloads and batch namespace;
  the tests cannot run until the S3 consumer migrations compile. Formatting and
  diff checks pass. Compilation is red at the planned S3 boundary:
  editor/transaction and reaction consumers still use the former dative payload
  or lack the new Undo arms. S3b, S3d, and S3e migrate those consumers before
  the S3 green gate.

- **S3a2 — Stereo update construction** (`ir::{edit,stereo}`; additive cleanup, within the S3 migration). [dep: S3a1]

  **Stereo update construction cleanup** (`ir::{edit,stereo}`). In both
  Edits::update_stereo_atom and update_stereo_bond, replace
  current.update(update) with computation of the configuration alone. Reuse
  StereoConfigurationUpdate::apply_to(&self, &StereoConfigurationForm) ->
  StereoConfigurationForm; change this existing method from private to
  pub(crate) so ir::edit can call it. No new method or public signature is added.
  Preserve configuration-update semantics, kind selection for constraint edits,
  normalized equality comparisons, emitted entry order, and old/new payloads.
  Keep the existing per-constraint construction; do not clone and update the
  complete constraints collection only to discard it. Ordinary form update
  methods retain their complete-form behavior.

  Verify exact edits for configuration-only, constraint-only, combined, clearing,
  and unchanged updates. Cover kind-only updates with matching and differing
  current kinds. Review the construction path for discarded whole-constraint
  updates; no benchmark campaign or timing threshold is required.

  **Implemented — 2026-09-26.** Both Edits constructors compute only the updated
  configuration; `StereoConfigurationUpdate::apply_to` is crate-visible. Exact
  edit cases cover the listed update shapes and kind selection. Form updates
  retain their complete-form behavior. Formatting and diff checks pass. The
  S3a1 consumer errors still prevent compilation until the S3 migrations.

- **S3a3 — Edit DSL `replace-*` operations** (`dsl::edit`; breaking, green at S3e). [dep: S3a1, S3a2]

  Extend the existing EditInput grammar and EditsDsl parse/render/conversion
  paths with :replace-atoms, :replace-donors, :replace-acceptor, :replace-site,
  and :replace-ligands under the applicable entity keys. Each uses
  [handle {:expect old :update new}], for example:

  ```edn
  {:aromatic-system
   {:replace-atoms [0 {:expect [0 1 2] :update [0 1 3]}]}}
  ```

  Lists are vectors in stored order; noncovalent atoms have fixed arity two.
  Sites/acceptors use existing handle encoding, and ligand lists use existing
  ligand encoding. EditsDsl::from_ir is exhaustive: support every new variant
  without omission or fallback. Delta and reaction DSL additions belong to later
  S3 subitems, not this Edit DSL change.

  Verify tree/streaming paths where present, ordered old/new vectors, existing
  handle and ligand encodings, dative syntax preservation, and roundtrips for
  every variant. Update Edit/Edits construction and namespace rustdoc.

  **Implemented — 2026-09-26.** EditInput parsing, rendering, and IR conversion
  handle all nine operations beside the existing edits for each entity. Separate
  EditInput and EditsDsl roundtrip cases cover stored vector order, handles, and
  ligand encoding. Negative cases cover missing update, noncovalent arity, and
  unknown ligand kind. Dative add/remove syntax is unchanged. Formatting and diff
  checks pass. Tests cannot run yet because the
  planned S3 consumer migration still leaves the crate uncompilable.

- **S3b** (`ir::molecule::transact`; breaking, green at S3e) Realize those edits
  through the structural methods on *EditorViewMut, obtained
  from Molecule's crate-private access introduced in S2i. Undo uses those same methods
  with the saved components. Do not reach into typed sets from edit execution
  or construct views there from Molecule fields. Keep unrelated factors,
  attributes, constraints, and ids unchanged. Test each
  forward/undo pair, old-state errors before mutation, and frame alignment.
  Migrate dative addition/removal and undo capture/replay to separate donors
  and acceptor; remove combined-vector joining/splitting and pass donor slices
  to the editor. Preserve existing failure and rollback semantics.
  [dep: S2j, S3a]

  **Implemented — 2026-09-26.** All nine edits check resolved old values before
  mutating through editor views and record their exact previous components for
  undo. Dative add/remove and restoration use separate donor and acceptor
  fields. Focused forward/undo, mismatch, invalid-target undo, and dative tests
  pass. S3d retains the dative reaction call-site migration and updates the
  property-test inputs to the separate donor/acceptor signature. The Python
  consumer migration remains in S3e.

  **Structural Edit execution and undo.** Each row below uses the existing
  *EditorViewMut for the resolved entity id, obtained through private
  Molecule access (S2i). Resolve the target and every handle in old/new before
  writing. Compare the old component exactly through view getters: sequences
  in stored order, stereo ligands including atom and kind, and scalar sites or
  acceptors by id. Reject mismatch with OldStateMismatch. Capture the actual
  previous component with resolved ids, then make the listed call. These edits
  neither renumber entities nor change attributes, constraints, or other factors.

  | Edit | View | Mutation call | Undo variant and replay on the same view |
  | --- | --- | --- | --- |
  | ReplaceAromaticSystemAtoms | AromaticSystemViewMut | replace_atoms(&new) | RestoreAromaticSystemAtoms: replace_atoms(&atoms) |
  | ReplaceMulticenterBondAtoms | MulticenterBondViewMut | replace_atoms(&new) | RestoreMulticenterBondAtoms: replace_atoms(&atoms) |
  | ReplaceNoncovalentBondAtoms | NoncovalentBondViewMut | replace_atoms(new) | RestoreNoncovalentBondAtoms: replace_atoms(atoms) |
  | ReplaceDativeBondDonors | DativeBondViewMut | replace_donors(&new) | RestoreDativeBondDonors: replace_donors(&donors) |
  | ReplaceDativeBondAcceptor | DativeBondViewMut | replace_acceptor(new) | RestoreDativeBondAcceptor: replace_acceptor(acceptor) |
  | ReplaceStereoAtomSite | StereoAtomViewMut | replace_site(new) | RestoreStereoAtomSite: replace_site(site) |
  | ReplaceStereoAtomLigands | StereoAtomViewMut | replace_ligands(&new) | RestoreStereoAtomLigands: replace_ligands(&ligands) |
  | ReplaceStereoBondSite | StereoBondViewMut | replace_site(new) | RestoreStereoBondSite: replace_site(site) |
  | ReplaceStereoBondLigands | StereoBondViewMut | replace_ligands(&new) | RestoreStereoBondLigands: replace_ligands(&ligands) |

  Replay guards target access, then writes the saved component without an
  expected-post-value comparison. It does not call row restoration or directly
  access a set. This table supplies the nine replacement rows of S4b's complete
  Edit inventory; S4b supplies the remaining variants and transaction integration.

- **S3c — replacement Delta vocabulary withdrawn; removal complete**
  (`ir::delta`, `ir::canonicalize`; breaking, graph-ir verification in S3d).
  [dep: S3b]

  **Semantics and interface.** Remove the nine ReplaceAtoms, ReplaceDonors,
  ReplaceAcceptor, ReplaceSite, and ReplaceLigands variants from the six overlay
  Delta enums. Remove their inverse, normalization/folding, frame transport, and
  canonicalization branches. Retain the existing Add/Remove and field/constraint
  semantics, including normalized equality and stereo-kind handling. Add no new
  Delta variants or conversions. Replacement Edits and Undo remain implemented.

  Remove cases specific to the deleted vocabulary; preserve the existing laws
  and the attribute/constraint cases from S3d1. The latter now test attributes
  directly through the functions specified in S3d.

- **S3d — selective removal and retained reaction integration; complete**
  (`ir::delta`, `ir::reaction`, `ir::reaction::integrity`, `ir::reaction_span`,
  `dsl::reaction`; breaking, graph-ir green checkpoint). [dep: S3b, S3c]

  This replaces the former S3d1–S3d5 replacement-integration sequence.

  **Semantics.** Remove replacement Delta lowering, reference visitation, and
  reaction DSL parsing/rendering. Keep direct span construction and ordinary
  reaction application. Remove the early created-id scans introduced solely for
  replacement lowering; use the existing collected additions to assign handles.
  Do not infer replacement from Add/Remove or route application through a span.
  Keep dative Edit addition/removal inputs split into donors and acceptor.

  **Modification interfaces.** Each function takes only the form and its Delta,
  returning Result<(), Contradiction>. They apply ModifyField/ModifyConstraint;
  Add/Remove leave the form unchanged. Existing old-value and constraint-key
  checks remain. The functions retain pub(crate) visibility and full entity names:

  ```rust
  pub(crate) fn apply_atom_modification(
      attributes: &mut AtomForm, delta: &AtomDelta,
  ) -> Result<(), Contradiction>;
  pub(crate) fn apply_bond_modification(
      attributes: &mut BondForm, delta: &BondDelta,
  ) -> Result<(), Contradiction>;
  pub(crate) fn apply_dative_bond_modification(
      attributes: &mut DativeBondForm, delta: &DativeBondDelta,
  ) -> Result<(), Contradiction>;
  pub(crate) fn apply_aromatic_system_modification(
      attributes: &mut AromaticSystemForm, delta: &AromaticSystemDelta,
  ) -> Result<(), Contradiction>;
  pub(crate) fn apply_multicenter_bond_modification(
      attributes: &mut MulticenterBondForm, delta: &MulticenterBondDelta,
  ) -> Result<(), Contradiction>;
  pub(crate) fn apply_noncovalent_bond_modification(
      attributes: &mut NoncovalentBondForm, delta: &NoncovalentBondDelta,
  ) -> Result<(), Contradiction>;
  pub(crate) fn apply_stereo_atom_modification(
      attributes: &mut StereoAtomForm, delta: &StereoAtomDelta,
  ) -> Result<(), Contradiction>;
  pub(crate) fn apply_stereo_bond_modification(
      attributes: &mut StereoBondForm, delta: &StereoBondDelta,
  ) -> Result<(), Contradiction>;
  ```

  Span construction and stereo normalization call these functions. Retain the
  approved reaction interfaces and removal of the public check_preconditions
  wrapper:

  ```rust
  impl Reaction {
      fn prepare_deltas(&self) -> Result<Deltas, ApplyPreconditionError>;
      fn apply_deltas_at(
          &self, host: &Molecule,
          correspondence: &MoleculeCorrespondence, deltas: Deltas,
      ) -> Result<Option<(Molecule, MoleculeCorrespondence)>, ApplyError>;
  }
  fn reframe_deltas(
      deltas: Deltas, lhs: &Molecule, host: &Molecule,
      correspondence: &MoleculeCorrespondence,
  ) -> Result<Deltas, ApplyError>;
  ```

  **Verification.** Check graph-ir tests compile; run delta, reaction, reaction-span,
  canonicalization, reaction-DSL, and replacement Edit/Undo cases. Run feature-gated
  delta, reaction, and Edit properties without changing their laws. Migrate the
  dative property inputs to separate donors/acceptor while preserving their bonds.
  Check graph-ir Clippy and rustdoc with warnings denied, nightly formatting,
  and the full diff. No replacement-Delta benchmark is needed. Workspace/Python
  enum migration closes in S3e; Rust 1.87 remains a final gate.

  **Checked — 2026-09-26.** All 1,413 selected tests pass: Delta, reaction,
  reaction-span, canonicalization, reaction DSL, Edit/Undo execution, and the
  feature-gated Delta/reaction/Edit properties, plus the selected integration
  targets. The dative property inputs in edit.rs and strategies.rs use separate
  donors and acceptor; their assertions and generated bonds are unchanged.
  Graph-ir all-target Clippy with proptest and rustdoc pass with warnings denied.
  Nightly formatting, removal searches, and full diff review pass. Replacement
  Edit/Undo and Edit DSL source files are unchanged by this removal. Workspace
  and Python verification is recorded under S3e.

- **S3e — implemented** (`umol-py::edit`; breaking, red→green)
  Expose the nine replacement Edit variants with the Rust payloads and failure
  behavior; migrate exhaustive matches and add parity cases.
  Migrate existing dative Edit variants and Edits addition/removal methods to
  separate donor/acceptor inputs, including plural entry tuples. Delta bindings
  retain their existing vocabulary. Do not add unrelated API coverage.
  This closes the enum migration. [dep: S3a, S3b, S3d]

  **Implemented — 2026-09-26.** Python Edit exposes all nine Replace variants,
  ordered with their entity kinds and carrying id/old/new as in Rust. Dative
  additions take donors, acceptor, attributes; removals take
  (id, donors, acceptor, attributes). The singular/plural Edits methods use those
  same factors. Binding conversion preserves Id/New handles, list order,
  noncovalent pairs, and ligand kinds. Rust execution supplies the existing
  precondition errors.

  **Checked — 2026-09-26.** The rebuilt Python 3.13 extension passes 1,638 tests
  (two skipped); the Rust binding suite passes 1,637 tests (two ignored).
  All nine variants have conversion, Python construction/rendering, application,
  old-value mismatch, and rollback coverage; a same-batch New-handle case checks
  dative additions and both replacements together. Workspace all-target checking,
  binding all-target Clippy, and binding rustdoc pass with warnings denied for
  lint and documentation. Nightly formatting and full diff review pass.
  Full-workspace tests and Rust 1.87 remain at S9b.

- **S3f — implemented** (`umol-graph-core::graph`; additive, green) Implement Graph::extend_nodes,
  extend_edges, and extend with the exact interfaces under Storage delegation. Return
  owned, allocation-free exact-size id iterators; mutation is eager and the
  iterators borrow neither receiver nor inputs. Preserve existing ids and append
  each kind contiguously. Endpoints supplied to extend use the resulting node space;
  invalid endpoints retain the graph addition panic contract. Single additions
  delegate to bulk operations. Node-only addition extends adjacency offsets;
  combined addition builds final adjacency once. Verify empty additions,
  appended ids and endpoints, loops/parallel edges, independence of shared
  clones, and continued graph mutation while returned iterators are held.
  Establish correctness cases and a focused native comparison of repeated versus
  bulk addition before replacing the implementation; use independently built
  expected graphs rather than single additions that delegate to the subject.
  [dep: none]

  **Implemented — 2026-09-26.** All three bulk methods return owned exact-size
  iterators, and add_node/add_edge delegate to them. Node-only addition grows
  offsets, copying shared storage once when necessary; batches containing edges
  rebuild adjacency once. Empty batches preserve shared storage. Invalid endpoints
  panic.

  **Measurements — 2026-09-26.** The graph benchmark has a setup function and
  separate extend_nodes, extend_edges, and extend benchmark functions. Removal
  and pushout benchmarks also reside there, in separate groups, retaining their
  fixtures and timing boundaries from algorithms.rs. Addition fixtures
  are paths of 64 or 1,024 nodes, with batches of one or 16 additions and unique
  or shared storage. Construction and final destruction are excluded; detachment
  is timed. Combined batches add equally many nodes and edges extending the path.
  Times below are microseconds per batch of 16. The original column measures
  repeated single additions before S3f; the other columns use the implementation
  above on the same fixtures.

  | Nodes / storage | Addition | Original singles | Current singles | Bulk |
  | --- | --- | ---: | ---: | ---: |
  | 64 / unique | nodes | 13.784 | 0.189 | 0.101 |
  | 64 / unique | edges | 13.764 | 12.883 | 0.858 |
  | 64 / unique | nodes + edges | 33.644 | 14.193 | 0.923 |
  | 64 / shared | nodes | 13.407 | 0.345 | 0.265 |
  | 64 / shared | edges | 13.728 | 12.558 | 0.778 |
  | 64 / shared | nodes + edges | 28.558 | 14.786 | 0.836 |
  | 1,024 / unique | nodes | 130.300 | 0.344 | 0.335 |
  | 1,024 / unique | edges | 133.880 | 130.450 | 9.361 |
  | 1,024 / unique | nodes + edges | 263.460 | 132.070 | 9.141 |
  | 1,024 / shared | nodes | 129.380 | 2.259 | 2.354 |
  | 1,024 / shared | edges | 133.390 | 130.010 | 9.175 |
  | 1,024 / shared | nodes + edges | 262.880 | 132.090 | 9.197 |

  These results support bulk calls for known edge batches: 16-edge additions
  take about 14–16 times less time than current repeated singles. Node-only gains
  mainly come from avoiding adjacency reconstruction and also benefit add_node.
  Shared node-only addition still copies the existing storage; at 1,024 nodes
  batching does not materially reduce that cost. Single-edge batches remain
  roughly equal to add_edge. Edge addition still rebuilds the full adjacency;
  these synthetic paths do not establish application-wide speedups. No allocation
  counts were measured.

  Reproduce with `cargo bench -p umol-graph-core --bench graph -- graph_extend
  --warm-up-time 0.1 --measurement-time 0.2 --sample-size 20 --noplot`.

  **Checked — 2026-09-26.** All 2,085 graph-core tests pass with proptest enabled,
  including independent topology/adjacency expectations, loops and parallel edges,
  clone independence, invalid endpoints, and iterator lifetimes.
  All-target Clippy with proptest and rustdoc pass with warnings denied. The graph
  benchmark's 48 addition cases pass; the eight moved removal/pushout cases pass
  in Criterion test mode. Both benchmark targets pass Clippy with warnings denied.
  Nightly formatting and full diff review pass.
  Full-workspace tests and Rust 1.87 remain at S9b.

- **S3g — implemented** (`umol-graph-core::relation::{fixed,var,fixed_fixed,fixed_var,var_var}`;
  additive, green) Add public extend to all five sets, retaining individual add.
  The precise batch argument for each set is in the relation-extend table under
  Storage delegation: an owned Vec of its factor/payload tuples, with variable
  factors borrowed as slices and fixed factors passed as arrays. Every method
  takes &mut self and returns an owned, allocation-free
  `impl ExactSizeIterator<Item = RelationId>` without receiver/input-lifetime
  captures. Move payloads, preserve row/factor order and multiplicity, and append
  all rows before rebuilding incidence once. Empty batches are no-ops. Match
  add's external-reference and size contracts; add no tracked_extend, new entry
  types, IntoIterator inputs, or add_many alias. Establish independent row-model
  cases and a focused native repeated-add/bulk-extend comparison before
  implementation. Check returned ids, both incidence factors, coinciding rows,
  empty variable factors, payload moves, and mutation while the returned iterator
  remains live. Do not use a delegating individual add as the sole oracle.
  [dep: none]

  **Implemented — 2026-09-26.** All five sets expose extend with the batch shapes
  above. Fixed factors use arrays; variable factors take slices copied into packed
  storage; payloads move without a Clone bound. Each nonempty call reserves room
  for the batch in each column/buffer, appends its rows, and rebuilds union incidence
  once. Empty batches leave storage unchanged. Precise captures exclude receiver
  and input-slice lifetimes; the exact-size id iterator owns its bounds. Individual
  add implementations, constructors, and their external-reference/size contracts
  are unchanged. The five methods are the only added public symbols.

  **Measurements — 2026-09-26.** Separate per-set benchmarks in benches/relation.rs
  compare repeated add with extend. Each row contains four overlapping node ids;
  birelations additionally contain one edge id. Inputs contain 64 or 1,024 rows;
  batches contain one or 16 rows. Setup clones and batch construction, and final
  set destruction, are excluded. Batch consumption, storage growth, and incidence
  rebuilding are timed. Microseconds per batch of 16:

  | Set | Existing rows | add before S3g | add after S3g | extend |
  | --- | ---: | ---: | ---: | ---: |
  | FixedRelationSet | 64 | 53.989 | 57.252 | 3.834 |
  | FixedRelationSet | 1,024 | 880.080 | 858.140 | 56.338 |
  | VarRelationSet | 64 | 55.083 | 56.088 | 4.160 |
  | VarRelationSet | 1,024 | 921.940 | 843.200 | 58.169 |
  | FixedFixedBirelationSet | 64 | 60.606 | 64.545 | 4.380 |
  | FixedFixedBirelationSet | 1,024 | 856.090 | 935.070 | 58.431 |
  | FixedVarBirelationSet | 64 | 63.421 | 62.934 | 4.607 |
  | FixedVarBirelationSet | 1,024 | 874.850 | 884.590 | 62.012 |
  | VarVarBirelationSet | 64 | 73.547 | 73.669 | 5.500 |
  | VarVarBirelationSet | 1,024 | 999.470 | 1,056.900 | 71.754 |

  Known batches should use extend: the measured 16-row calls take about 13–16
  times less time than repeated add. One-row calls show no consistent gain;
  extend ranged from about 4% faster to 12% slower in this short comparison.
  Individual add remains appropriate for isolated additions. Both operations
  rebuild the full incidence index; batching reduces repeated work without
  changing storage or adding persistent indexes. No allocation counts were
  measured, and these synthetic fixtures do not establish application-wide gains.

  Reproduce with `cargo bench -p umol-graph-core --bench relation -- _extend/
  --warm-up-time 0.1 --measurement-time 0.2 --sample-size 20 --noplot`.

  **Checked — 2026-09-26.** All 2,164 graph-core tests pass with proptest enabled,
  including 74 new unit cases and five new properties. Unit cases cover zero/empty
  factors, coinciding rows, sparse references, non-Clone payloads, exact-size
  iteration, and receiver/input-borrow independence. Sequential-batch properties
  use independent rows and direct incidence scans, including references shared
  across factors and participants reporting both node and edge ids. Removing all
  appended ids restores the original rows and index. The 40 benchmark cases pass;
  all-target Clippy with proptest, rustdoc with warnings denied, nightly formatting,
  and full diff review pass. Workspace and Rust 1.87 gates remain at S9b.

- **S3h — completed 2026-09-28** (`ir::{aromatic,multicenter,dative,noncovalent,stereo}`;
  additive, green) Add crate-private extend to AromaticSystems, MulticenterBonds,
  DativeBonds, NoncovalentBonds, StereoAtoms, and StereoBonds. Each takes &mut self
  and the exact Vec batch specified for the corresponding plural Molecule
  addition under Storage delegation; it returns an owned exact-size iterator
  of that set's domain ids, without borrowing the set or participant inputs.
  Preserve individual add. Convert ids and delegate the batch to relation
  extend once, retaining copy-on-write storage and moving attributes. Verify
  order, ids, incidence, unchanged shared originals, empty batches, and domain
  distinctions: dative donors/acceptor and stereo site/ligands.
  [dep: S1a, S1b, S1c, S3g]

  All six sets implement the specified crate-private extend. Each nonempty
  batch calls its relation set's extend once. Aromatic, multicenter, and dative
  atom-id conversion uses one contiguous temporary node vector; stereo ligand
  slices pass directly to graph-core. Attributes move into storage. Empty batches
  retain the shared allocation, and returned id iterators retain no borrows.

  **Checked — 2026-09-28.** The 36 new unit cases cover empty/populated sets,
  unique/shared ownership, empty batches, row/factor order, repeated atoms,
  independent incidence scans, dative roles, stereo sites and all ligand kinds,
  and iterator use after inputs and receivers are dropped. All 50 tests selected
  by `cargo test -p umol-graph-ir --lib _extend` pass, along with graph-ir all-target
  Clippy, rustdoc with warnings denied, nightly formatting, and diff checks/review.
  Graph-core's S3g properties cover the delegated batch algorithm; Molecule/editor
  integration is S3i. Workspace and Rust 1.87 gates remain at S9b.

- **S3i — completed 2026-09-28** (`ir::molecule`, `ir::molecule::editor`; additive, green) Add crate-private
  Molecule and public editor add_atoms, add_bonds, add_dative_bonds,
  add_aromatic_systems, add_multicenter_bonds, add_noncovalent_bonds,
  add_stereo_atoms, and add_stereo_bonds. Use the full bulk-addition interface
  table under Storage delegation: &mut self, owned Vec batches containing
  borrowed variable participant slices and moved forms, and
  `impl ExactSizeIterator<Item = Id> + use<>` for each concrete domain Id.
  Editor methods delegate. Molecule add_atoms/add_bonds call Graph bulk addition
  and extend the matching attribute vector; overlay methods call the owning
  set's extend. Preserve existing ids and input order, with eager mutation and
  no retained input/receiver borrow. Add no journal, correspondence, publication
  check, bulk Edit variant, or edit-coalescing behavior. Verify ids and complete
  stored entries against independent expected values, graph/attribute alignment,
  empty batches, and use through editor publication. [dep: S2i, S3f, S3h]

  All eight crate-private Molecule additions and public editor delegates are implemented.
  Graph and its atom/bond attribute vectors grow together; each overlay addition
  calls its typed set's extend. Forms move into storage, and the returned iterators
  own their bounds. Empty batches preserve shared storage.

  The existing editor still owns session correspondence until S6b1. Nonempty
  additions extend its target counts once per batch so current tracked publication
  and subsequent removals remain correct. This increments one count through a
  mutable borrow without cloning pair vectors. The crate-private Molecule
  methods do not maintain correspondence; S6b1 removes this editor bookkeeping
  together with the existing single-addition bookkeeping.

  **Checked — 2026-09-28.** The 24 new unit cases cover every addition with unique
  and shared storage, complete published entries, id order, empty batches, and
  iterator use after input slices and the editor are gone. A public property
  reconstructs generated MoleculeEntries through consecutive bulk additions of
  varying size, checking ids, publication, and the current session correspondence.
  `cargo test -p umol-graph-ir --features proptest` passes: 7,951 tests, including
  415 properties; seven tests are ignored. Strict graph-ir Clippy with all targets
  and proptest, rustdoc with warnings denied, nightly workspace formatting, and
  `git diff --check` pass. The full diff was reviewed against S3i and local
  conventions. Workspace and Rust 1.87 gates remain at S9b.

- **S3j — completed 2026-09-28** (`umol-graph-core::correspondence`,
  `ir::correspondence`, `ir::molecule::{editor,transact}`; breaking, green after
  caller migration). [dep: S3i]

  Change the three mutation methods on Correspondence, GraphCorrespondence, and
  MoleculeCorrespondence from consuming receivers to mutable borrows. A
  correspondence stored in another structure must be directly mutable without
  moving it out and installing an empty placeholder. Retain public visibility,
  names, arguments, existing generic bounds, and error types. The nine signatures
  are:

  ```rust
  impl<Id> Correspondence<Id> {
      pub fn extend_right(&mut self, count: usize);
      pub fn compact_right(
          &mut self, compaction: &Compaction<Id>,
      ) -> Result<(), CorrespondenceComposeError>;
      pub fn uncompact_right(
          &mut self, compaction: &Compaction<Id>,
      ) -> Result<(), CorrespondenceComposeError>;
  }

  impl GraphCorrespondence {
      pub fn extend_right(&mut self, nodes: usize, edges: usize);
      pub fn compact_right(
          &mut self, compaction: &GraphCompaction,
      ) -> Result<(), GraphCorrespondenceComposeError>;
      pub fn uncompact_right(
          &mut self, compaction: &GraphCompaction,
      ) -> Result<(), GraphCorrespondenceComposeError>;
  }

  impl MoleculeCorrespondence {
      pub fn extend_right(&mut self, kind: EntityKind, count: usize);
      pub fn compact_right(
          &mut self, compaction: &MoleculeCompaction,
      ) -> Result<(), MoleculeCorrespondenceComposeError>;
      pub fn uncompact_right(
          &mut self, compaction: &MoleculeCompaction,
      ) -> Result<(), MoleculeCorrespondenceComposeError>;
  }
  ```

  **Semantics and failure.** extend_right increases the selected right counts;
  all pairs remain unchanged. compact_right discards pairs whose right ids were
  removed and translates survivors. uncompact_right expands surviving right ids;
  restored positions remain unmatched, and discarded pairs are not recreated.
  All three preserve left counts and reuse pair-vector allocations. Zero extension
  and identity compaction/uncompaction leave the value unchanged.

  compact_right requires each right count to equal its compaction's source count;
  uncompact_right requires the result count instead. On a mismatch, return the
  existing error and leave the entire receiver unchanged. GraphCorrespondence and
  MoleculeCorrespondence check every applicable count before mutating any pairs.
  Preserve existing error selection: nodes before edges, and molecule entity-kind
  order. These checks require no clone or journal. Preserve the existing algorithms
  and the atom/bond removed-id adaptation; this is a receiver-contract change.

  **Migration.** Update graph-core first, then MoleculeCorrespondence, then editor
  additions/removals and undo restoration. Replace the correspondence-specific
  mem::replace/assignment sequences with direct method calls. Migrate chained
  calls in tests, properties, and the existing correspondence benchmark. Keep
  current session tracking until S6b1; S3j changes how its value is updated.
  Construction, conversion, induction, queries, compose, compose_all, and reverse
  retain their contracts. The Python bindings do not expose these three mutation
  methods and need no new surface. Update rustdoc and the data-type guide's
  consuming-correspondence section to describe mutable updates and error behavior.

  **Verification.** Retain exact-result, identity, pair-allocation reuse, and
  composition-equivalence coverage. Add unchanged-receiver assertions on count
  errors at all three levels, including a later mismatching component when earlier
  components would change. Exercise direct mutation through a containing field.
  Run graph-core and graph-ir tests with proptest, strict all-target Clippy, and
  rustdoc with warnings denied; adapt and compile the existing benchmark without
  starting a new measurement campaign. Run nightly formatting and diff checks,
  and review the full change against these nine signatures. Workspace and
  Rust 1.87 gates remain at S9b.

  All nine methods use mutable receivers with the signatures above. Graph and
  molecule compaction/uncompaction check all counts before updating pairs.
  Editor additions/removals and undo restoration call the methods directly;
  their 30 move-out/assignment sequences are removed. Existing session tracking
  remains until S6b1. Rustdoc and the data-type guide describe the new contracts.

  **Checked — 2026-09-28.** Both crates pass with proptest: 10,115 tests passed,
  including 624 properties; seven tests are ignored. The 26 count-error cases
  assert unchanged receivers, with nonidentity compactions exercising late
  mismatches. Existing editor tests cover correspondence mutation through its
  containing field. Strict all-target Clippy, including the migrated benchmark,
  rustdoc with warnings denied, nightly formatting, diff checks, and full diff
  review pass. No new performance measurements were needed.

- **S3k — completed 2026-09-28** (index and count arithmetic; group, green).
  [dep: S3j]

  Indices, positions, offsets, and collection counts are assumed to fit their
  representations. The u32 range is much larger than the molecules that can
  reasonably be represented and manipulated. Integer conversions for these
  values are infallible; use ordinary arithmetic. Do not add checked or saturating
  arithmetic, maximum-value comparisons, or conversion errors to guard index
  overflow.

  Retain bounds checks against actual collections, endpoint/reference checks,
  and compatibility checks between supplied correspondence/compaction domains.
  Their purpose is different from testing integer representability. The Python
  accessor counter retains its separately settled overflow behavior. Chemical
  values and external-format value ranges are not indices and retain their
  existing conversion contracts. No new APIs, numeric wrappers, or generic
  conversion helpers are needed.

- **S3k1 — completed 2026-09-28**
  (`umol-graph-core::graph`, `docs/development/data-types.md`; cleanup, green).
  [dep: S3j]

  Record the rule above in the data-type guide and remove its instruction to
  guard index-size arithmetic and capacity limits during restoration. Keep the
  matching-history and manipulated-history contracts, with index arithmetic
  governed by the common representable-size assumption.

  In Graph::extend_nodes and Graph::extend, replace the four checked additions
  with ordinary additions and remove the three u32-limit assertions. In
  Graph::restore, remove the saturated subtraction and u32/isize capacity
  comparisons. Preserve the empty-node/nonempty-edge guard and the existing
  local id and endpoint checks. All method names, visibility, arguments, and
  return types remain unchanged; add_node/add_edge/extend_edges keep delegating.
  Update their Panics sections and the MoleculeEditor bulk-addition rustdoc to
  describe actual endpoint/index bounds, without an integer-overflow contract.

  Remove the graph extension capacity-panic test and the three restoration
  capacity cases. Retain ordinary addition/restoration fixtures, clone
  independence, invalid endpoints, and mismatched-history coverage. Run the
  focused graph tests and graph/restoration properties; no new benchmark campaign.

  Graph uses ordinary arithmetic for additions and has no index-overflow guards
  in restoration. The unused edge-count calculation is removed with its guard.
  Endpoint/domain checks remain. Graph and MoleculeEditor rustdoc and the data-type
  guide follow the index/count rule. The four extension capacity cases and three
  restoration capacity cases are removed; all other tests and property laws remain.

  **Checked — 2026-09-28.** All 175 Graph unit tests and 13 properties covering
  graph operations, restoration, and common-subgraph consumers pass. Nightly
  formatting, diff checks, and review against the S3k1 scope pass.

- **S3k2 — completed 2026-09-28** (`umol-graph-core::relation`;
  cleanup, green). [dep: S3k1]

  Apply the same rule to restore and restore_participants on FixedRelationSet,
  VarRelationSet, FixedFixedBirelationSet, FixedVarBirelationSet, and
  VarVarBirelationSet. Replace checked row/participant-count additions and
  try_fold accumulation with ordinary addition and fold. Remove u32-limit
  comparisons and the saturated source-count checks in participant restoration.
  Use ordinary allocation/reservation instead of the try_reserve early-return
  branches. Keep the existing restoration algorithms and allocation sizes.

  Preserve all ten signatures, payload movement, reference translation, current
  surviving values, and incidence rebuilding. Keep guarded row/slot accesses
  and reference-domain checks. Remove only oversized-count/capacity cases from
  the five relation test modules; retain sparse ids, missing/duplicate/out-of-range
  rows, matching-history results, and the independent restoration properties.
  Run the focused relation restoration tests/properties and graph-core Clippy.

  All five sets use ordinary count arithmetic and allocation in restoration.
  Row/slot guards, reference-domain checks, signatures, and restoration algorithms
  are unchanged. The 15 oversized-count test cases are removed; other cases and
  property laws are unchanged.

  **Checked — 2026-09-28.** All 221 focused Graph/relation restoration unit tests
  and 38 restoration properties pass. Graph-core all-target Clippy with proptest
  and warnings denied, nightly formatting, diff checks, and scope review pass.

- **S3k3 — completed 2026-09-28**
  (`ir::{aromatic,dative,multicenter,noncovalent,stereo,molecule::transact}`;
  cleanup, green). [dep: S3k2]

  Replace all 16 u32::try_from(position).expect(...) conversions with the direct
  ParticipantPosition(position as u32) construction. The affected methods are
  replace_atom/insert_atom/remove_atom on AromaticSystems and MulticenterBonds;
  replace_donor/insert_donor/remove_donor on DativeBonds; replace_atom on
  NoncovalentBonds; and replace_ligand/insert_ligand/remove_ligand on StereoAtoms
  and StereoBonds. Retain their current usize arguments, visibility, storage
  delegation, and bounds behavior for positions in the supported index range.

  Replace the four checked count additions in molecule::transact::reconstruction_fits,
  molecule::transact::restored_constraints, and the RestoreRemovedTopology arm of
  MoleculeEditor::validate_undo with ordinary addition. Keep this edit limited to
  arithmetic; S4 still owns removal of the old undo-validation machinery and
  restoration rewiring. Do not introduce overflow checks into that replacement.

  Run the existing affected entity-mutation and rollback tests, the S4a1
  restoration tests, and graph-ir Clippy. Preserve invalid-id and invalid-position
  tests that check actual collection bounds. Public view signatures are unchanged.

  All 16 position conversions use ParticipantPosition(position as u32); the four
  undo count calculations use ordinary addition. Signatures, bounds checks,
  storage delegation, and tests are unchanged.

  **Checked — 2026-09-28.** All 429 selected entity-mutation, transaction/rollback,
  and Molecule restoration tests pass. Graph-IR all-target Clippy with proptest
  and warnings denied, nightly formatting, diff checks, and scope review pass.

- **S3k4 — completed 2026-09-28** (`umol-graph::export`; breaking enum
  cleanup, green). [dep: S3k3]

  In the reaction Convey implementation, use pairs.len() as u32 to generate the
  same one-based labels. Remove ReactionConveyError::AtomMapLabel { index: usize },
  whose sole producer is the count conversion. The enum retains Materialization,
  Reactants, and Products; the convey signature and all other errors stay intact.
  Update the method's Errors section. Retain checks on chemical attribute values
  and other external-format representation limits.

  Run the existing reaction convey tests and graph Clippy; verify exhaustive
  consumers compile. Finish this group with warnings-denied rustdoc for the
  affected Rust crates, nightly formatting, and a diff review against the scope
  above. Workspace/Python rebuilds and Rust 1.87 remain at the final gate.

  Reaction convey uses pairs.len() as u32 for the label count. AtomMapLabel is
  removed from ReactionConveyError and the Errors section describes the three
  remaining variants. No other code referenced the removed variant; labels,
  mapping order, chemical-value checks, and tests are unchanged.

  **Checked — 2026-09-28.** All 36 reaction convey/export tests pass. Graph
  all-target Clippy with conformance/proptest and warnings denied passes.
  Warnings-denied rustdoc passes for graph-core, graph-ir, and graph. Nightly
  formatting, diff checks, and S3k scope review pass; S3k is complete.

### S4 — Recovery machinery before the public lifecycle switch

- **S4a — completed 2026-09-28** (`ir::constraint::molecule`, `ir::molecule`, `ir::molecule::editor`,
  `ir::molecule::editor::transact`; group; public rename and restoration rewire, green at S4a2). [dep: S2i, S3b]

- **S4a1 — completed 2026-09-28** (`ir::molecule`; additive, green). [dep: S2i, S3b]

  Route topology and relation restoration through the existing graph-core
  operations and add local attribute and target-access guards. Constraint
  restoration and its guards belong to S4a2.
  Guard missing or out-of-range entries; do not add index-overflow checks.
  Implement crate-private Molecule::restore_topology with
  `(&mut self, &GraphCompaction, Vec<RemovedAtom>, Vec<RemovedBond>) -> ()`.
  Its scope is Graph plus the atom/bond attribute vectors. S4a2 removes the editor's
  combined topology-and-overlay restoration body; undo execution composes this
  primitive with separate crate-private Molecule delegates to each set's existing
  restore_topology_ids and restore, followed by restore_constraints.
  Overlay-only undo calls the affected row-restoration delegate without
  topology-id restoration. Keep these as independent operations.

  Each row defines two crate-private Molecule methods. The row-restoration method
  takes `(&mut self, rows: &Compaction<Id>, removed: Vec<Entry>) -> ()`.
  The topology-id method takes `(&mut self, topology: &GraphCompaction) -> ()`.
  Id is the first member of Entry. Each method forwards to just its owning set;
  removed entries use the pre-removal topology ids, as required by set restore.

  | Row-restoration method | Entry | Topology-id method |
  | --- | --- | --- |
  | restore_dative_bonds | (DativeBondId, Vec<AtomId>, AtomId, DativeBondForm) | restore_dative_bond_topology_ids |
  | restore_aromatic_systems | (AromaticSystemId, Vec<AtomId>, AromaticSystemForm) | restore_aromatic_system_topology_ids |
  | restore_multicenter_bonds | (MulticenterBondId, Vec<AtomId>, MulticenterBondForm) | restore_multicenter_bond_topology_ids |
  | restore_noncovalent_bonds | (NoncovalentBondId, [AtomId; 2], NoncovalentBondForm) | restore_noncovalent_bond_topology_ids |
  | restore_stereo_atoms | (StereoAtomId, AtomId, Vec<StereoLigand>, StereoAtomForm) | restore_stereo_atom_topology_ids |
  | restore_stereo_bonds | (StereoBondId, BondId, Vec<StereoLigand>, StereoBondForm) | restore_stereo_bond_topology_ids |

  Verify each method against matching removals, surviving topology-id translation,
  unchanged-topology restoration, and manipulated-input panic freedom.

  All 13 crate-private methods are implemented. restore_topology delegates topology
  restoration to Graph::restore and restores only the atom/bond attribute vectors.
  Saved attributes move into their original positions; surviving attributes keep
  their current values and move when their storage is uniquely owned. Each overlay
  delegate calls only its owning set's restore or restore_topology_ids method.

  **Checked — 2026-09-28.** All 79 focused restoration tests pass, covering matching
  removals, changed surviving attributes, shared-storage independence, identity,
  separate overlay operations, and manipulated history. Strict graph-ir all-target
  Clippy with proptest, rustdoc with warnings denied, nightly formatting, and diff
  checks pass. The full diff was reviewed against the 13 approved signatures.

- **S4a2 — completed 2026-09-28** (`ir::constraint::molecule`, `ir::molecule::{editor,editor::transact}`; breaking, green). [dep: S4a1, S3k]

  Move the current molecule/transact.rs to molecule/editor/transact.rs and
  declare it as a private child of editor. Keep MoleculeEditor.molecule and
  MoleculeEditor.correspondence private. The current editor batch methods and
  detached-journal rollback then access editor state within its owning module.
  Re-export Transaction and TransactionError through editor and molecule,
  preserving their public paths; migrate internal imports to those re-exports.
  Move the correspondence-allocation test into the editor tests without importing
  another test module. Make the seven editor remove_added_* methods private;
  their callers are now descendants of editor. Use pub(crate) for from_parts and
  the existing internal Molecule mutable accessors, additions, and restorations.
  This changes placement and visibility, not lifecycle or mutation semantics.

  Add public `Constraints::extend(&mut self, constraints: Vec<Constraint>) -> ()`,
  moving entries in order without normalization or reference checks. Retain the
  existing push, remove_at, retain, clear, and take interfaces listed under
  Molecule-level constraint mutation. Explicit constraint-removal edits save
  the value returned by remove_at with its original position; restoration uses
  Constraints::restore, with no new tracked-removal method. Test duplicate/order
  preservation, empty extension, and explicit removal followed by restoration.
  Rename Constraints::compact_with_update to
  `tracked_compact(&mut self, &MoleculeCompaction) -> CascadedConstraints` and
  migrate its callers. Add public
  `Constraints::restore(&mut self, &CascadedConstraints)` in place of
  transaction-level restored_constraints; preserve recorded positions and old
  values, guard local accesses, and remove post-value equality requirements.
  Retain `compact(&mut self, &MoleculeCompaction)`; add no tracked_restore.
  Add crate-private Molecule
  `restore_constraints(&mut self, changes: &CascadedConstraints) -> ()`,
  delegating only to Constraints::restore. These delegates do not add old-value
  checks, journals, or publication checks. Migrate undo callers to them; preserve
  the storage restoration contracts, including panic freedom for manipulated
  history. Remove S4a1's temporary per-method dead_code expectations when these
  methods gain production callers. Test topology-id restoration and row restoration
  together and independently where topology is unchanged, covering all six overlay kinds.
  Keep the existing detached journal surface until S5. Test matching-history
  recovery under `normalized_eq` and panic freedom for manipulated undo data
  without asserting its result. Cover constraint compaction/restoration with
  removed and rewritten entries, duplicate entries, and preserved list order.

  Constraint storage exposes extend, tracked_compact, and restore. Undo removal
  branches call the Molecule restoration methods; constraint reconstruction and
  its local guards are owned by Constraints. Explicit removal saves remove_at's
  returned value and position. Editor restoration bodies and all 13 restoration
  dead-code expectations are removed. Both editor fields are private;
  editor::transact owns the current batch and detached-journal implementation.
  Shared execution moves to molecule::apply in S4b3; borrowed transaction
  ownership moves to molecule::transact in S5. Aggregate undo validation remains
  until its scheduled removal in S4b8.

  **Checked — 2026-09-28.** All 603 focused constraint, editor, transaction,
  tracked-publication, and Molecule restoration unit cases pass. All 26 affected
  constraint/edit property tests pass, including independently generated
  restoration histories. Graph-ir all-target Clippy with proptest and warnings
  denied, warnings-denied rustdoc, nightly formatting, and diff checks/review
  pass. Workspace/Python rebuilds and Rust 1.87 remain at S9b.

- **S4b — Molecule delegation and Edit/Undo execution** (group; breaking,
  green at S4b8; cleanup at S4b9). [dep: S4a, S3i]

  Execute S4b1–S4b9 in order. The complete editor inventory is in S4b2;
  the execution signatures are in S4b3 and every Edit/Undo call is mapped below.
  Changes to journal return types and Undo variants may temporarily break callers
  within this group; migrate all of them by S4b8 without fallback match arms.

- **S4b1 — completed 2026-09-28** (`ir::molecule`; breaking rewire, green at S4b8). [dep: S4a, S3i]

  Move the eight add_* methods,
  push_constraint, and all seven untracked/tracked removal pairs to crate-private
  Molecule methods, with the interfaces recorded under Storage delegation.
  Editor additions delegate; complete editor removals and batch execution compose
  the same primitives.
  Remove MoleculeEditor::push_constraint and migrate its callers to
  editor.constraints_mut().push(constraint). Retain crate-private
  Molecule::push_constraint for Edit execution. Do not replace the removed editor
  method with another top-level constraint convenience delegate.
  Reuse S3i's bulk add_atoms/add_bonds for the corresponding Edit variants;
  preserve individual overlay Edit execution and its existing undo boundaries.
  Rename editor remove to remove_topology and crate-private tracked_remove to
  tracked_remove_topology, migrating callers. The crate-private topology pair mutates
  only Graph and atom/bond attributes and returns GraphCompaction when tracked;
  each overlay pair mutates only its owning set and returns Compaction<Id> when
  tracked. Public editor removal still completes all cascading updates.

  Add the following crate-private Molecule compaction methods. Each takes
  `(&mut self, topology: &GraphCompaction)` and returns the listed mapping.
  Each installs the set returned by its owning set's compact operation. They
  mutate no other component and perform no extra validation or journal recording.
  The row mapping is required for dependent updates; there is no separate
  mapping-discarding method.

  | Method | Return |
  | --- | --- |
  | compact_dative_bonds | Compaction<DativeBondId> |
  | compact_aromatic_systems | Compaction<AromaticSystemId> |
  | compact_multicenter_bonds | Compaction<MulticenterBondId> |
  | compact_noncovalent_bonds | Compaction<NoncovalentBondId> |
  | compact_stereo_atoms | Compaction<StereoAtomId> |
  | compact_stereo_bonds | Compaction<StereoBondId> |

  Complete the crate-private Molecule constraint delegation surface:

  ```rust
  pub(crate) fn push_constraint(&mut self, constraint: Constraint);
  pub(crate) fn extend_constraints(&mut self, constraints: Vec<Constraint>);
  pub(crate) fn remove_constraint_at(&mut self, position: usize) -> Constraint;
  pub(crate) fn compact_constraints(&mut self, compaction: &MoleculeCompaction);
  pub(crate) fn tracked_compact_constraints(
      &mut self, compaction: &MoleculeCompaction,
  ) -> CascadedConstraints;
  // restore_constraints is supplied by S4a.
  ```

  Delegate to push, extend, remove_at, compact, and tracked_compact respectively.
  Preserve their order, duplicate, and failure semantics; remove_constraint_at
  retains remove_at's out-of-range panic. Batch execution finds the last equal
  constraint and handles absence before calling it. No tracked constraint
  removal is needed. The existing editor constraints_mut access continues to
  support retain, clear, take, and whole-collection assignment.

  Verify component boundaries, untracked/tracked storage equivalence, bulk
  delegation, and unchanged unrelated components.

  Molecule owns the recorded individual additions, component-only removal and
  compaction pairs, and constraint delegates, all pub(crate). Editor additions
  delegate; its removals compose these primitives and retain cascading updates.
  Rust/Python topology removal is named remove_topology/tracked_remove_topology;
  editor push_constraint is removed. Existing tracked editor removals remain
  public until S4b2 removes that direct surface.

  **Checked — 2026-09-28.** All 1,067 focused molecule/editor/correspondence unit
  cases and 24 affected property tests pass. The 57 new unit cases cover isolated
  storage changes, tracked/untracked equivalence, unchanged shared originals,
  and constraint order/duplicates. The rebuilt Python extension passes all 164
  transaction/import tests. Warnings-denied rustdoc, nightly formatting, and
  diff checks/review pass. All-target graph-ir Clippy with proptest completes
  with only dead_code warnings for new primitives awaiting their planned
  production callers; no suppression attributes were added.

- **S4b2 — completed 2026-09-28** (`ir::molecule::editor`; breaking rewire, green at S4b8). [dep: S4b1]

  **Complete editor direct-mutation inventory.** The following is the resulting
  public surface, including mutation reached through returned views/borrows.
  Addition signatures are in Storage delegation; all take &mut self and
  delegate to synonymous crate-private Molecule methods.

  | Individual addition | Bulk addition |
  | --- | --- |
  | add_atom | add_atoms |
  | add_bond | add_bonds |
  | add_dative_bond | add_dative_bonds |
  | add_aromatic_system | add_aromatic_systems |
  | add_multicenter_bond | add_multicenter_bonds |
  | add_noncovalent_bond | add_noncovalent_bonds |
  | add_stereo_atom | add_stereo_atoms |
  | add_stereo_bond | add_stereo_bonds |

  Removal methods take &mut self and return (). They compose Molecule removal
  and compaction methods, rather than merely delegating to the synonymous
  component-only removal. No public tracked direct-removal methods remain.

  | Editor removal | Arguments after &mut self |
  | --- | --- |
  | remove_topology | atoms: &[AtomId], bonds: &[BondId] |
  | remove_dative_bonds | ids: &[DativeBondId] |
  | remove_aromatic_systems | ids: &[AromaticSystemId] |
  | remove_multicenter_bonds | ids: &[MulticenterBondId] |
  | remove_noncovalent_bonds | ids: &[NoncovalentBondId] |
  | remove_stereo_atoms | ids: &[StereoAtomId] |
  | remove_stereo_bonds | ids: &[StereoBondId] |

  Mutable accessors are atom_mut, bond_mut, dative_bond_mut,
  aromatic_system_mut, multicenter_bond_mut, noncovalent_bond_mut,
  stereo_atom_mut, and stereo_bond_mut. Each takes its typed entity id and
  returns the corresponding *EditorViewMut<'_> obtained from Molecule (S2i5).
  Every view exposes attributes_mut() -> &mut Form, including direct mutation
  of entity-level constraints. Structural mutation uses these view methods,
  all taking &mut self and returning () (S2j):

  | View | Structural methods |
  | --- | --- |
  | Atom / Bond | None; localized-bond endpoint rewiring is outside this scope |
  | AromaticSystem / MulticenterBond | replace_atoms(&[AtomId]); replace_atom(AtomPosition, AtomId); insert_atom(AtomPosition, AtomId); remove_atom(AtomPosition) |
  | NoncovalentBond | replace_atoms([AtomId; 2]); replace_atom(AtomPosition, AtomId) |
  | DativeBond | replace_donors(&[AtomId]); replace_donor(AtomPosition, AtomId); insert_donor(AtomPosition, AtomId); remove_donor(AtomPosition); replace_acceptor(AtomId) |
  | StereoAtom | replace_ligands(&[StereoLigand]); replace_ligand(StereoLigandPosition, StereoLigand); insert_ligand(StereoLigandPosition, StereoLigand); remove_ligand(StereoLigandPosition); replace_site(AtomId) |
  | StereoBond | The same ligand methods; replace_site(BondId) |

  Top-level constraints use constraints_mut(&mut self) -> &mut Constraints:
  push(Constraint), extend(Vec<Constraint>), remove_at(usize) -> Constraint,
  retain(impl FnMut(&Constraint) -> bool), clear(), take() -> Vec<Constraint>,
  and whole-collection assignment. The other listed collection methods return
  (). There is no editor push_constraint, extend_constraints, or
  remove_constraint_at delegate. Molecule supplies the mutable access; editor
  code does not reach into its fields.

  Undo-only restoration and internal compaction methods are not editor entry
  points. Batch apply/tracked_apply and checked probe/finish remain the separate
  lifecycle surface specified above. Verify this inventory against the final
  editor API and migrate all removed push_constraint callers, including bindings
  and tests where present.

  The seven tracked direct-removal methods are crate-private, and their Python
  bindings are removed. Public removals return () and retain cascading updates.
  The inventory matches all 16 additions, seven removals, eight mutable entity
  accessors, constraints_mut, and the structural view methods above. Additions
  and mutable access already delegate to Molecule; removal composes its storage
  primitives. No view types or mutation semantics changed. External compaction
  properties use public transaction results; their input domains and inverse/
  composition assertions are preserved. The retired tracked-removal benchmark
  is removed; direct-removal and batch benchmarks remain.

  **Checked — 2026-09-28.** All 972 molecule/editor unit cases, three affected
  compaction properties, and 157 Python transaction/import cases pass. The
  Python extension was rebuilt with Python 3.13. Warnings-denied rustdoc,
  nightly formatting, and diff checks/review pass. All-target graph-ir Clippy
  with proptest reports only S4b1's unused primitives awaiting planned callers.

- **S4b3 — completed 2026-09-28** (`ir::molecule::{apply,editor::transact}`; breaking, green at S4b8). [dep: S4b2]

  Move apply_edit, apply_edit_with_undo, and apply_undo from impl MoleculeEditor
  to impl Molecule in a private ir::molecule::apply module. Move ApplicationState,
  HandleTable, and execution-only functions with them. Keep editor batch methods
  and detached Transaction lifecycle in editor::transact until S5; they call the
  shared execution methods rather than access one another's fields. S5 moves the
  borrowed Transaction implementation to molecule::transact.

  The three Molecule methods, ApplicationState, and its new constructor are
  pub(crate). ApplicationState fields, HandleTable, and handle-processing methods
  remain private to apply. Both forward methods take one Edit and the batch's
  mutable ApplicationState; apply_edit_with_undo returns an optional Undo.
  apply_undo takes one Undo. Initialize ApplicationState from Molecule accessors.
  Editor and Transaction retain their batch loops; Transaction retains journal
  storage and reverse replay. Add no execution methods to Edit, Edits, or Undo
  and no separate execution type. Execution calls Molecule mutation methods and
  accessors; its placement does not authorize direct storage access.

  ```rust
  // In ir::molecule::apply, on Molecule:
  pub(crate) fn apply_edit(&mut self, edit: Edit, state: &mut ApplicationState)
      -> Result<(), TransactionError>;
  pub(crate) fn apply_edit_with_undo(&mut self, edit: Edit, state: &mut ApplicationState)
      -> Result<Option<Undo>, TransactionError>;
  pub(crate) fn apply_undo(&mut self, undo: Undo);
  // Existing batch state; fields stay private:
  pub(crate) struct ApplicationState { /* eight HandleTable fields */ }
  impl ApplicationState {
      pub(crate) fn new(molecule: &Molecule) -> Self;
  }
  ```

  Execution returns the existing per-Edit error category. Public lifecycle
  methods wrap it in MoleculeApplyError::Transaction; integrity belongs to
  probe/finish/commit. Undo has no Result or expected-post-value check. The shared
  crate-private method interface keeps editor and transaction fields private.

  Retain ApplicationState's eight existing HandleTable fields, one per entity
  kind, and its typed lookup/push/compact methods. Replace eager initial-id
  vectors with this private layout:

  ```rust
  struct HandleTable<I> {
      initial_count: usize,
      initial: Option<Vec<Option<I>>>,
      created: Vec<Option<I>>,
  }
  // I: Copy + From<usize>
  HandleTable::new(initial_count: usize) -> Self;
  ```

  Before compaction, initial lookup checks initial_count and returns I::from
  the original index; None here means identity, not a removed entity. On the
  first compaction, materialize the original-id mapping and retain None slots
  for removals. Later compactions update that mapping and created. Existing
  initial/created lookup error distinctions and New numbering stay unchanged.
  Test lookups before/after removal, additions preceding removal, and removed
  handles; do not add another state type or public constructor.

  Keep execution on Molecule single-entry; add no internal batch wrappers.
  Editor::apply initializes batch handle state and calls apply_edit for each
  entry. Transaction::apply initializes batch handle state and calls
  apply_edit_with_undo, appending each returned Undo to its journal. Rollback
  traverses that journal in reverse and calls apply_undo. This keeps completed
  undo entries available to Transaction when a later Edit fails. A journaled
  batch wrapper would need to accept the journal or return completed entries on
  failure; its benefit would be encapsulating iteration, not reducing storage
  work. Plural Edit variants already use the corresponding bulk mutation
  methods where provided by this plan.

  After forward topology removal, call all six Molecule compact_* overlay
  delegates and assemble MoleculeCompaction. After explicit overlay removal,
  assemble that mapping with identity components for unchanged kinds. Then
  call compact_constraints or tracked_compact_constraints. Addition undo instead
  calls only the crate-private Molecule untracked removal, as listed below. Remove the
  seven remove_added_* adapters and make those calls in the Undo match arms,
  without overlay/constraint compaction or
  correspondence updates. Keep forward sequencing explicit; add no combined
  coordination helper or recording-mode parameter solely to avoid repetition.
  Remove all sixteen apply_modify_*_field / apply_modify_*_constraint helpers.
  Keep handle resolution and old-value/precondition checks in batch execution,
  and perform field and entity-constraint writes through mutable views obtained
  from Molecule. Top-level constraint changes use the crate-private Molecule delegates.
  Structural Edit replacements and their undos use *EditorViewMut methods wired
  in S3b. Add no Molecule modification family. Prepare fallible data before
  writes, preserve each bundled edit as one
  execution/undo unit, and avoid initial receiver-sized handle tables until
  compaction needs them. Plain removal uses
  compact_constraints; transactional removal uses tracked_compact_constraints
  once on the actual collection and records its result in the bundled undo.
  Remove the whole-collection recovery clone and duplicate compaction; replay
  constraints through restore_constraints after restoring entity storage.
  Constraint capture is the returned value of tracked_compact; no optional
  mutable output parameter is added to removal. Check all current
  non-constraint Edit families and the nine planned replacements against the
  mutation-coverage table. Test that a rejected
  member of a bundled edit causes no writes from that edit, and that a later
  edit failure restores earlier completed edits. Verify that bare and tracked
  compaction install equivalent sets, returned mappings describe the installed
  rows, and each delegate leaves other components unchanged. Review all editor,
  batch, and undo mutation paths: calls must use Molecule methods or views,
  with no direct writes to Molecule storage, direct owning-set mutation, view
  construction from its fields, or copy-on-write manipulation in those callers.
  Keep multi-component sequencing explicit; this access boundary does not
  justify a combined helper. No partial-write records or injected internal-panic
  recovery tests.

  **Edit-to-mutation inventory.** The tables here and S3b enumerate all 33
  current Edit variants and nine planned replacements. Calls without a view
  receiver are crate-private Molecule methods; m is the owned or borrowed Molecule.
  Accessors and views supply all reads, including preconditions and undo capture.
  No execution or replay branch borrows Molecule fields, constructs a view from
  those fields, or directly mutates Graph, attribute vectors, typed sets, or
  copy-on-write storage. The sequences below are instructions for the match
  arms, not additional bundled helper methods.

  Resolve all handles and check all preconditions for the whole Edit before its
  first write, including every member of a plural edit. Invalid/dead handles
  retain their existing errors; repeated removal targets return DuplicateRemoval.
  Store resolved ids in undo, not batch-local handles. Publish one bundled undo
  for the completed Edit, except the entity-constraint None-to-None no-op, which
  produces no journal entry. Replay completed edits in reverse order. Journal-free
  execution makes the same forward calls without constructing undo payloads.
  Neither path adds per-Edit aggregate integrity checks. During replay, retain
  local access guards for manipulated undo data and omit offered-old checks;
  this is not permission to re-enter the checked forward interpreter.

  Verify batch-local handles, completed-Edit journal entries, and no entry for
  a None-to-None constraint edit. Retain the existing public lifecycle until S5.

  **Implemented — single-entry execution and batch ownership.** The private
  molecule::apply module owns ApplicationState, HandleTable, the three
  crate-private execution methods, and existing execution-only functions.
  Initial handle maps allocate only on compaction; New handles retain their
  per-kind ordinals. Editor loops call Molecule and append only Some(Undo);
  all eight None-to-None entity-constraint edits produce no entry. The six
  frame-comparison methods and their tests moved with execution and are private.
  Structural writes use Molecule's editor-view accessors. Fields remain private;
  no execution type, batch wrapper, or visibility exception was added.

  **Verification — 2026-09-28.** Nightly formatting and diff checks pass.
  cargo check -p umol-graph-ir --all-targets reaches library and unit-test
  targets but fails in the scheduled removal/replay migrations: component
  compactions still need S4b5's complete mappings; addition undos still call
  the adapters replaced in S4b4/S4b8; the moved replay body still needs S4b8's
  infallible implementation. Nineteen new handle/no-op cases and the migrated
  comparison cases are not yet executable. No passing test result is claimed.
  S4b4–S4b7 retain the per-family rewiring and S4b8 closes the migration,
  including the retained editor tracked batch/rollback behavior: Molecule
  execution no longer updates the editor's private session correspondence.

- **S4b4 — completed 2026-09-28** (`ir::molecule::apply`; rewire, green at S4b8). [dep: S4b3]

  **Additions.** Resolve the listed topology handles before calling the method;
  attribute forms are moved unchanged. Register returned ids in the batch's New
  namespace in input order. Capture the corresponding existing Added* records,
  with resolved ids, stored frames, and attributes, for the listed Undo variant.

  | Edit | Preconditions beyond the common rules | Molecule mutation | Undo and removal call |
  | --- | --- | --- | --- |
  | AddAtoms | No topology handles | add_atoms(atoms) | RemoveAddedTopology with AddedAtom entries: remove_topology(&ids, &[]) |
  | AddBonds | Resolve both endpoints of every bond | add_bonds(bonds) | RemoveAddedTopology with AddedBond entries: remove_topology(&[], &ids) |
  | AddDativeBond | Resolve donors and acceptor separately | add_dative_bond(&donors, acceptor, attributes) | RemoveAddedDativeBond: remove_dative_bonds(&[id]) |
  | AddAromaticSystem | Resolve all atoms | add_aromatic_system(&atoms, attributes) | RemoveAddedAromaticSystem: remove_aromatic_systems(&[id]) |
  | AddMulticenterBond | Resolve all atoms | add_multicenter_bond(&atoms, attributes) | RemoveAddedMulticenterBond: remove_multicenter_bonds(&[id]) |
  | AddNoncovalentBond | Resolve the two atoms | add_noncovalent_bond(atoms, attributes) | RemoveAddedNoncovalentBond: remove_noncovalent_bonds(&[id]) |
  | AddStereoAtom | Resolve atom site and all ligand atoms | add_stereo_atom(site, &ligands, attributes) | RemoveAddedStereoAtom: remove_stereo_atoms(&[id]) |
  | AddStereoBond | Resolve bond site and all ligand atoms | add_stereo_bond(site, &ligands, attributes) | RemoveAddedStereoBond: remove_stereo_bonds(&[id]) |

  The last column is the entire mutation needed for each addition undo. These
  are crate-private Molecule primitives, not the editor's complete cascading removal
  operations. Reverse replay has undone later references and compactions, so
  matching added entries again occupy trailing positions and earlier ids stay
  unchanged. Do not compact overlays or constraints, assemble MoleculeCompaction,
  update correspondence, or create new undo. The Undo match arms
  only extract ids and call these methods, with local guards for panic freedom
  on manipulated input; they neither validate history nor inspect saved payloads.

  Verify ids, handle registration, saved entries, and reverse removal after
  later dependent edits have been undone.

  **Implemented — addition execution.** AddAtoms and AddBonds use one bulk
  Molecule call per Edit in both execution paths. Bond handles are all resolved
  before insertion. Journaled additions save the stored attributes and resolved
  frames, with one bundled Undo per Edit; returned ids populate each kind's New
  namespace in order. Stereo forms move into storage and their resolved ligand
  vectors are retained for undo without another frame allocation.

  The seven addition-undo branches call only untracked Molecule removal, with
  local bounds guards. All seven editor remove_added_* adapters are removed;
  their topology regression now exercises Molecule::apply_undo directly.
  Tests cover all eight additions, full stored results and undo payloads, New
  handles, reverse replay of dependent additions, rejected bond handles at each
  position, and panic freedom for out-of-range added-entry undos.

  **Verification — 2026-09-28.** Nightly formatting and diff checks/review pass.
  The all-target graph-ir check still fails on S4b5's removal-compaction mapping
  and S4b8's fallible replay body; addition-adapter errors are gone. The new
  execution cases cannot run until that migration compiles. Aggregate undo
  validation remains for removal in S4b8; no passing test result is claimed.

- **S4b5 — completed 2026-09-28** (`ir::molecule::apply`; rewire, green at S4b8). [dep: S4b4]

  **Topology removal.** RemoveTopology resolves atom/bond targets and rejects
  duplicate ids within each list. No offered-old payload is present. Journaled
  execution captures, before writing, removed atoms, explicit and incident bonds, and all overlays
  that will lose a required topology reference, using molecule/entity accessors.
  Capture stored order and attributes. Stereo capture includes sites and virtual
  ligand-bearing atoms; a stereo-bond site drops if its bond is explicitly or
  incidentally removed. Do not clone all top-level constraints for capture.

  Forward calls, in order:

  1. tracked_remove_topology(&atoms, &bonds) returns GraphCompaction and updates
     Graph plus the atom/bond attributes together.
  2. Call each overlay delegate in the following table with that GraphCompaction.
     Each installs its replacement set and returns the corresponding row mapping.
  3. Assemble MoleculeCompaction::new from the graph and six row mappings.
  4. Call tracked_compact_constraints(&compaction) for journaled execution, or
     compact_constraints(&compaction) for journal-free execution.
  5. Update the batch handle mappings with the combined compaction. Record
     RestoreRemovedTopology with the captured entries, compaction,
     undo_compaction, and returned CascadedConstraints.

  | Forward compaction call | Undo: surviving topology-id call | Undo: removed-row call |
  | --- | --- | --- |
  | compact_dative_bonds | restore_dative_bond_topology_ids | restore_dative_bonds |
  | compact_aromatic_systems | restore_aromatic_system_topology_ids | restore_aromatic_systems |
  | compact_multicenter_bonds | restore_multicenter_bond_topology_ids | restore_multicenter_bonds |
  | compact_noncovalent_bonds | restore_noncovalent_bond_topology_ids | restore_noncovalent_bonds |
  | compact_stereo_atoms | restore_stereo_atom_topology_ids | restore_stereo_atoms |
  | compact_stereo_bonds | restore_stereo_bond_topology_ids | restore_stereo_bonds |

  Undo first calls restore_topology(compaction.graph(), atoms, bonds). Then,
  for each row above, restore surviving topology ids with compaction.graph()
  and restore removed rows with that kind's row compaction and saved entries
  using S4a's signatures. Finally call restore_constraints(&cascade). Do not
  combine these calls into another Molecule restoration operation.

  **Explicit overlay removal.** Resolve every target and offered topology
  handle; compare every offered entry before any removal. Preserve the existing
  comparison: align the offered frame to the stored frame and compare transported
  attributes by normalized_eq. Ordinary unordered factors use DynPermutation;
  stereo uses Permutation, with stereo-bond endpoint blocks kept intact. Failed
  alignment or unequal payload yields OldStateMismatch. Journaled execution
  captures actual stored entries in their original order, not the offered frame.

  | Edit | Entry comparison | Molecule removal | Undo variant and row restoration |
  | --- | --- | --- | --- |
  | RemoveDativeBonds | Same acceptor; donor multiset and aligned attributes | tracked_remove_dative_bonds(&ids) | RestoreRemovedDativeBonds: restore_dative_bonds(rows, entries) |
  | RemoveAromaticSystems | Atom multiset and aligned attributes | tracked_remove_aromatic_systems(&ids) | RestoreRemovedAromaticSystems: restore_aromatic_systems(rows, entries) |
  | RemoveMulticenterBonds | Atom multiset and aligned attributes | tracked_remove_multicenter_bonds(&ids) | RestoreRemovedMulticenterBonds: restore_multicenter_bonds(rows, entries) |
  | RemoveNoncovalentBonds | Unordered atom pair and aligned attributes | tracked_remove_noncovalent_bonds(&ids) | RestoreRemovedNoncovalentBonds: restore_noncovalent_bonds(rows, entries) |
  | RemoveStereoAtoms | Same atom site; aligned full ligand frame and attributes | tracked_remove_stereo_atoms(&ids) | RestoreRemovedStereoAtoms: restore_stereo_atoms(rows, entries) |
  | RemoveStereoBonds | Same bond site; aligned full ligand frame and attributes | tracked_remove_stereo_bonds(&ids) | RestoreRemovedStereoBonds: restore_stereo_bonds(rows, entries) |

  After each removal call, assemble MoleculeCompaction with its returned row
  mapping and identity mappings for unchanged kinds. Call
  tracked_compact_constraints once, record its cascade with the saved entries
  and undo_compaction, and update batch handle mappings. Journal-free execution
  calls compact_constraints instead. Undo obtains the saved row mapping from
  undo_compaction.forward(), calls the listed row restoration, then
  restore_constraints(&cascade). No topology-id restoration is required here.

  Verify whole-Edit precondition rejection, all overlay cascades, constraint
  positions, and batch handle updates from the composed mappings.

  **Implemented — removal execution.** Both execution paths compose the graph
  and six overlay compactions before compacting stored constraints and batch
  handles. Explicit overlay removals use identity mappings for unchanged kinds.
  Journaled removal records the cascade returned by tracked_compact_constraints;
  it no longer clones the complete constraint list or compacts a second copy.
  Target resolution, duplicate rejection, frame alignment, and stored-entry
  capture all precede removal. Replay uses the existing per-kind restoration
  calls followed by constraint restoration.

  Thirty added cases cover both execution paths: topology cascades, explicit
  bond-site removal, empty topology removal, all six overlay removals, reordered
  offered frames, constraint positions and duplicates, initial/New handle
  remapping, whole-Edit rejection before writes, and matching-history replay.

  **Verification — 2026-09-28.** Nightly formatting and diff checks/review pass.
  The all-target graph-ir check reports only the ten existing errors from `?`
  in the infallible apply_undo body, scheduled for S4b8. Removal-compaction
  errors are resolved. Tests remain unexecuted until that migration compiles.

- **S4b6 — completed 2026-09-28** (`ir::molecule::apply`; rewire, green at S4b8). [dep: S4b5]

  **Attribute changes.** Resolve the target, read the selected field through
  its view, and require normalized_eq with the offered old value before writing.
  A mismatch is OldStateMismatch. Assign new through the access below. Capture
  the resolved id and inverse FieldChange in the corresponding Undo::Modify*Field
  variant. Replay assigns the inverse change's new value through the same view,
  without comparing its old value. Equivalent offered old forms may be retained,
  consistently with restoration modulo normalized_eq.

  | Edit | Mutable access | FieldChange variant → assigned field | Undo variant |
  | --- | --- | --- | --- |
  | ModifyAtomField | m.atom_mut(id).attributes_mut() | Element → element; IsotopeMass → isotope_mass; Charge → charge; ImplicitHydrogens → implicit_hydrogens; LonePairs → lone_pairs; UnpairedElectrons → unpaired_electrons | ModifyAtomField |
  | ModifyBondField | m.bond_mut(id).attributes_mut() | Order → order; Charge → charge; UnpairedElectrons → unpaired_electrons | ModifyBondField |
  | ModifyDativeBondField | m.dative_bond_mut(id).attributes_mut() | Order → order | ModifyDativeBondField |
  | ModifyAromaticSystemField | m.aromatic_system_mut(id).attributes_mut() | Electrons → electrons; Charge → charge; UnpairedElectrons → unpaired_electrons | ModifyAromaticSystemField |
  | ModifyMulticenterBondField | m.multicenter_bond_mut(id).attributes_mut() | Electrons → electrons; Charge → charge; UnpairedElectrons → unpaired_electrons | ModifyMulticenterBondField |
  | ModifyNoncovalentBondField | m.noncovalent_bond_mut(id).attributes_mut() | Kind → kind | ModifyNoncovalentBondField |
  | ModifyStereoAtomField | m.stereo_atom_mut(id).attributes_mut() | Configuration → configuration | ModifyStereoAtomField |
  | ModifyStereoBondField | m.stereo_bond_mut(id).attributes_mut() | Configuration → configuration | ModifyStereoBondField |

  Bind each view locally while borrowing attributes; the table abbreviates the
  access path, not Rust temporary-lifetime handling. No per-field setters or
  apply_modify_* wrappers are introduced.

  Cover every FieldChange member and equivalent offered-old forms.

  **Implemented — attribute execution.** The eight field wrappers are removed.
  Both Edit executors select and assign fields directly through the listed
  mutable views. Journaled execution clones only the new field value and moves
  the offered old value into the inverse change. Field replay assigns the saved
  value without comparing old state, with local target bounds guards.

  Fifty-four cases cover every FieldChange member in both execution paths,
  equivalent offered-old forms, whole-molecule results, repeated-edit rejection,
  and restoration modulo normalized_eq. They include assignments of electron
  counts with differing lengths and out-of-range cosets. Eight further cases
  exercise panic freedom for field undo with missing targets.

  **Verification — 2026-09-28.** Nightly formatting and diff checks/review pass.
  The all-target graph-ir check reports only two existing errors: aggregate undo
  validation and ApplyEdit recursion still use `?` in infallible apply_undo.
  Those are removed in S4b7/S4b8. Tests remain unexecuted until the migration
  compiles.

- **S4b7 — completed 2026-09-28** (`ir::{edit,molecule::apply}`; breaking, green at S4b8). [dep: S4b6]

  **Entity-level constraints.** Resolve the target. If old/new are both Some,
  require equal keys. At that key, require current and old to be both absent or
  normalized_eq; otherwise return OldStateMismatch before writing. Both None
  is a no-op with no undo entry. The internal journaled executor returns an
  optional Undo, and its caller appends only a present entry; do not introduce
  a dummy undo or retain ApplyEdit for this case. These are the existing
  compare_and_set semantics, expressed at
  the Edit boundary followed by the collection's ordinary set/remove operations.

  | Edit | Collection reached through mutable view | Mutation | Undo variant |
  | --- | --- | --- | --- |
  | ModifyAtomConstraint | m.atom_mut(id).attributes_mut().constraints | set(new) when Some; remove(key) when None | RestoreAtomConstraint |
  | ModifyBondConstraint | m.bond_mut(id).attributes_mut().constraints | set(new) when Some; remove(key) when None | RestoreBondConstraint |
  | ModifyDativeBondConstraint | m.dative_bond_mut(id).attributes_mut().constraints | set(new) when Some; remove(key) when None | RestoreDativeBondConstraint |
  | ModifyAromaticSystemConstraint | m.aromatic_system_mut(id).attributes_mut().constraints | set(new) when Some; remove(key) when None | RestoreAromaticSystemConstraint |
  | ModifyMulticenterBondConstraint | m.multicenter_bond_mut(id).attributes_mut().constraints | set(new) when Some; remove(key) when None | RestoreMulticenterBondConstraint |
  | ModifyNoncovalentBondConstraint | m.noncovalent_bond_mut(id).attributes_mut().constraints | set(new) when Some; remove(key) when None | RestoreNoncovalentBondConstraint |
  | ModifyStereoAtomConstraint | m.stereo_atom_mut(id).attributes_mut().constraints | set(new) when Some; remove(key) when None | RestoreStereoAtomConstraint |
  | ModifyStereoBondConstraint | m.stereo_bond_mut(id).attributes_mut().constraints | set(new) when Some; remove(key) when None | RestoreStereoBondConstraint |

  Remove Undo::ApplyEdit(Box<Edit>). Add the eight variants above, each with
  fields id, key, and constraint. For entity prefix X, their types are XId,
  XConstraintKey, and Option<XConstraintForm>, respectively; all are existing
  domain types. For example:

  ```rust
  RestoreAtomConstraint {
      id: AtomId,
      key: AtomConstraintKey,
      constraint: Option<AtomConstraintForm>,
  }
  ```

  Record the resolved id, the key derived from old or new, and the accepted old
  optional constraint. Retaining an equivalent offered old value follows the
  settled normalized_eq restoration contract. Replay uses the listed mutable
  view: Some(saved) calls constraints.set(saved); None calls
  constraints.remove(key). Guard target access without checking an expected
  post-value or rerunning compare_and_set. Stereo kind remains forward Edit DSL
  context, not an application precondition, and is absent from these undos.

  **Top-level constraints.** Resolve all handles inside ConstraintEdit,
  including nested entries, before mutation. Constraints is an ordered list
  preserving duplicates; these operations do not perform an integrity check.

  | Edit | Additional precondition | Molecule mutation | Undo capture and replay |
  | --- | --- | --- | --- |
  | AddMoleculeConstraint | None | push_constraint(constraint) | RemoveAddedMoleculeConstraint { position }, recording m.constraints().len() before insertion. Replay guards the position, then calls remove_constraint_at(position). |
  | RemoveMoleculeConstraint | Find the last exact-equal entry through m.constraints(); MissingEntry if absent | remove_constraint_at(position) | RestoreMoleculeConstraints containing RemovedConstraint with that position and the returned stored value. Replay calls restore_constraints(&changes). |

  The ninth replacement variant is:

  ```rust
  RemoveAddedMoleculeConstraint {
      position: usize,
  }
  ```

  Reverse replay restores the append position before this undo runs, including
  when later edits removed or compacted constraints. No saved constraint payload,
  equality search, or forward RemoveMoleculeConstraint processing is needed.
  Guard an out-of-range position for panic freedom on manipulated history.
  Whole-entry restoration continues to use restore_constraints; inline
  constraint restoration uses the owning entity view. Migrate all nine former
  ApplyEdit producers and all Undo matches/tests in this subitem. No replacement
  boxed Edit, generic inverse interpreter, or ApplyEdit compatibility variant.

  Test saved optional constraints, keys, duplicate top-level entries, position
  restoration, and the no-op journal case.

  **Implemented — constraint execution.** Both Edit executors check keys and
  offered old values inline, then call set/remove through entity mutable views.
  The eight constraint wrappers and Undo::ApplyEdit are removed. Eight explicit
  constraint undos retain the resolved id, key, and accepted old optional value;
  replay restores through the same views with local target guards. Molecule
  constraints use push_constraint/remove_constraint_at and restore_constraints.
  Addition undo records and guards the insertion position.

  Added cases cover all eight entity kinds in both execution paths, optional
  values, equivalent offered-old forms, key/value mismatches, and no-op journals.
  Further cases cover target resolution for no-ops, duplicate molecule constraints,
  exact-equality removal, nested handle failures before writes, and restoration
  of an insertion position after later removal/compaction is undone. The former
  boxed-Edit receiver test is replaced by nine explicit-undo guard cases.

  **Verification — 2026-09-28.** Nightly formatting and diff checks/review pass.
  The all-target graph-ir check reports only the remaining aggregate undo
  validation `?` in infallible apply_undo, removed in S4b8. Tests remain unexecuted
  until that migration compiles.

- **S4b8 — completed 2026-09-28** (`ir::molecule::apply`, Rust consumers; breaking, red→green). [dep: S4b7]

  Implement Molecule::apply_undo through the calls mapped in S3b and S4b4–S4b7.
  Remove the seven remove_added_* adapters, forward Edit recursion, old-value
  comparisons, and aggregate undo validation from replay. Migrate detached
  rollback callers to this infallible private replay while retaining their
  public lifecycle until S5. No Python Undo wrapper or exhaustive match exists
  to migrate; the Python Transaction lifecycle changes in S5d.

  **Undo coverage.** These tables and S3b cover all 41 resulting Undo variants:
  seven undo additions, seven restore removals, eight restore fields, nine
  restore structural components, eight restore entity constraints, and two
  restore top-level constraint changes. No additional Molecule mutation method
  is needed for replay. Match arms select the operation or field, extract saved
  data and compactions, convert saved entries to the method argument tuples,
  and retain local guards for panic freedom on manipulated inputs. Storage
  restoration and constraint reconstruction belong to the called methods;
  replay performs no compaction algorithm, correspondence update, old-value
  comparison, aggregate undo validation, or forward Edit dispatch.

  **Inventory verification.** Cover every table row with forward mutation and
  matching-history undo, checking unaffected entities/components too. Exercise
  each listed FieldChange member, each constraint option shape, duplicate
  top-level constraints, and topology removal cascading to all six overlay kinds.
  Cover all nine new constraint Undo variants, no entry for None-to-None,
  equivalent old forms, target/position guards, and restoration of a saved
  constraint-addition position after later removal/compaction has been undone.
  Verify complete removal of Undo::ApplyEdit and forward Edit dispatch from replay.
  For addition undo, cover later dependencies and intervening removals undone
  first, then removal of the appended entries with earlier storage preserved.
  Verify removal of all seven remove_added_* adapters. Review the addition-undo
  branches for untracked removal only, without
  overlay/constraint compaction or correspondence bookkeeping.
  Rejected later members of plural edits leave that entire Edit unwritten;
  transaction failure after earlier successful edits restores transaction entry.
  Check source matches against this inventory without wildcard fallbacks. Audit
  every mutation and capture path for accessor-only Molecule access; passing
  behavior tests alone does not establish this boundary.

  **Implementation and verification.** All 41 Undo variants call Molecule
  mutation methods or mutable views. Replay has no forward Edit dispatch,
  remove_added_* adapters, aggregate validation, or old-value comparison.
  The detached Transaction lifecycle remains until S5.

  restore_topology returns before mutation when the attribute counts differ
  from the compaction's survivor counts. This prevents unrelated-history
  restoration from creating a length mismatch that makes a later removal or
  field undo panic. Saved entries are sorted by original id and merged with
  surviving attributes directly into the final vectors. Each affected kind
  allocates one final attribute vector; there is no Option vector or default
  filling. Forms move, with the existing Arc::make_mut copying only when the
  attribute storage is shared.

  All 1,198 molecule unit cases and 52 affected public Edit/frame/reframe
  properties pass. Six regression cases replay an addition/removal journal on
  empty, smaller, and larger unrelated molecules, through both rollback paths.
  Matching restoration covers reversed saved entries, changed survivor values,
  and shared storage. Malformed-entry tests exercise the merge as well as count
  mismatches. The unrelated-history property checks rollback panic freedom; its
  subsequent try_tracked_build comparison was removed because manipulated
  history has no usable-result guarantee. Matching-history correctness checks
  remain. Nightly formatting and diff review pass. S4b9 retains the unused-method
  cleanup and strict lint gate.

  **Overlay restoration.** VarRelationSet and FixedVarBirelationSet restore by
  sorting saved rows and merging them with surviving row ranges into final
  columns. This covers aromatic systems, multicenter bonds, dative bonds, stereo
  atoms, and stereo bonds. It removes the Option row buffers and the intermediate
  copy of saved participants through the old packed buffer. Payloads move without
  a Clone bound. Input buffers are released before incidence rebuilding. The
  trade-off is allocating final payload/offset/fixed columns rather than reusing
  their retained capacity. Noncovalent restoration retains its in-place reordering.

  Native before/after comparison: existing restore benchmark, sparse width-four
  rows, usize payloads, reversed saved entries; 0.1 s warm-up, 0.2 s measurement,
  10 samples. Times in microseconds; each cell is before → after. Retained cases
  preserve column capacity from removal; single removes the middle row.

  | Storage / rows | Single | Single, retained | Interleaved | Interleaved, retained |
  | --- | ---: | ---: | ---: | ---: |
  | Var / 64 | 2.02 → 1.92 | 1.82 → 1.98 | 2.34 → 2.10 | 2.12 → 2.15 |
  | Var / 1024 | 21.82 → 20.59 | 20.54 → 20.28 | 25.78 → 27.93 | 25.10 → 24.18 |
  | FixedVar / 64 | 3.58 → 3.42 | 3.32 → 3.43 | 3.88 → 3.70 | 3.64 → 3.72 |
  | FixedVar / 1024 | 38.18 → 37.57 | 37.14 → 37.12 | 42.99 → 41.76 | 42.54 → 42.03 |

  Retain the direct merge for its simpler storage construction and removal of
  intermediate copies. It is not uniformly faster: the 64-row retained Var
  single-removal case costs about 0.15 µs more (8%). The 1024-row ordinary Var
  interleaved estimate has a wide interval and no detected timing difference.
  No additional precision study is needed for this change. Allocation effects
  above come from code inspection; peak memory was not measured.

  Verification: 221 graph-core restoration unit cases, 38 restoration properties,
  1,198 molecule unit cases, and 21 public Edit properties pass. Graph-core strict
  all-target Clippy and private-item rustdoc pass. Existing tests cover unsorted
  entries, empty factors, non-Clone payloads, incidence, and manipulated history;
  their assertions and generators are unchanged.

- **S4b9 — completed 2026-09-28**
  (`ir::{aromatic,dative,multicenter,noncovalent,stereo,molecule}`; cleanup and internal rename, green).
  [dep: S4b8]

  Remove the six impl-level `#[cfg_attr(not(test), expect(dead_code, ...))]`
  attributes. The six owning overlay sets expose
  `compact(&self, &GraphCompaction) -> (Self, Compaction<Id>)`; their Molecule
  delegates expose `compact_<entity>(&mut self, &GraphCompaction) -> Compaction<Id>`.
  Remove the mapping-discarding methods and rename the mapping-returning methods
  to these bare names. The mappings are required for dependent constraint and
  handle updates. Keep graph-core's public compact/tracked_compact pairs.
  Migrate both Edit execution paths, direct editor removal, and tests. Preserve
  assertions on resulting storage, returned mappings, incidence, and restoration;
  remove only the duplicate cases and comparisons for the deleted counterparts.
  Use extend_constraints for lift_constraints' accumulated additions.

  Merge each owning set's inherent impl blocks in this order: constructor,
  immutable getters, mutable getters, mutators, destructuring methods
  (into_entries), helpers. Place its existing reframe_*_with free function
  directly after the inherent impl. Preserve method bodies and visibility.
  Verify the non-test build, strict all-target Clippy, the affected removal and
  compaction tests, and rustdoc. Add no dead-code allowances.

  **Implementation and verification.** The six owning sets each have one ordered
  inherent impl. Their compact methods and Molecule delegates return the required
  row mappings; the twelve mapping-discarding methods are removed. No new method,
  visibility change, or dead-code allowance is introduced. lift_constraints
  passes its accumulated additions to extend_constraints with order and duplicates
  preserved. All 2,785 focused entity/molecule unit cases pass. The non-test build,
  strict all-target Clippy, and private-item rustdoc with warnings denied pass.
  Nightly formatting and full diff review pass. This also closes S4d7's lint gate.

- **S4c — moved to S5a1.** Introduce the scope guard with Transaction::run,
  which owns it. No unused recovery machinery or temporary public API is added
  at the S4 green boundary.

- **S4d — Single-entity entries and comparison** (group; additive types,
  then caller migration; green). [dep: S1a, S1b, S1c]

  S4d1–S4d6 are additive to the owning sets and can precede S4b closeout.
  S4d7 retains S4b9 as a prerequisite.

  Each owning entity set gains entry(id), returning a complete entity value
  with its atom/site/ligand ids and attributes. The set's own entity id selects
  the entry and is absent from the returned value. The six entry types live in
  their respective entity modules, with Normalize, FrameTransport, and Reframe
  implementations. Comparison uses the actual normalized_eq and framed_eq trait
  methods. The six Molecule *_equiv methods disappear; no molecule comparison
  submodule or renamed comparison-wrapper family is introduced.

  **Common interface and contract.** Types, fields, and entry getters use
  pub(crate), retaining the agreed crate-internal comparison surface. Construct
  old entries with struct literals; add no constructor family, re-export,
  Python type, storage mutation, or view change. Derive Clone, Debug, PartialEq,
  and Eq; equality compares values independently of Cow's borrowed/owned state.
  These are open entry values, not checked molecule constructors. Entry access
  panics on an unavailable entity id, consistently with the existing set getters.

  Normalize changes only attributes and retains the supplied frame. It delegates
  attribute contradictions to the existing Normalize implementation; it adds no
  electron-count length, coset-range, or molecule-reference validation.
  FrameTransport checks action degree against the entry's actual frame, transports
  both that frame and its attributes, and returns None on incompatibility.
  Distinguished acceptors/sites remain fixed. Reframe uses the existing per-kind
  representative-action functions and normalize/transport/normalize semantics.
  Specialize framed_eq to transport one entry's attributes into the other's frame
  for comparison. Preserve the trait's equality contract, including two failed
  reframings comparing equal; use the definition when direct alignment does not
  establish equality. Keep the structural-equality shortcut.

  Cow permits an entry to borrow input and own a transformed result. The getter
  borrows attributes without cloning them. Dative, aromatic, and multicenter
  atom lists currently require NodeId-to-AtomId collection; record that allocation
  rather than changing graph-core storage or adding unsafe slice conversions.
  Noncovalent atom pairs are copied; stereo ligand slices remain borrowed.
  Transformation leaves the source set untouched.

  Place each entry definition and its trait implementations immediately before
  its owning set's definition. Place entry(id) with that set's immutable getters.
  Keep each entry's tests together, before its owning set's tests, with no imports
  between test modules. Each subitem covers exact getter contents, unavailable-id
  panic, borrowed/owned equality, source preservation, normalization idempotence,
  frame-action identity/inverse/composition, and reframing idempotence using valid
  frames. Existing public property suites remain outside src; do not widen these
  interfaces solely to expose them to property tests.

- **S4d1 — completed 2026-09-28** (`ir::dative`; additive, green). [dep: S1b]

  ```rust
  pub(crate) struct DativeBondEntry<'a> {
      pub(crate) donors: Cow<'a, [AtomId]>,
      pub(crate) acceptor: AtomId,
      pub(crate) attributes: Cow<'a, DativeBondForm>,
  }
  // On DativeBonds:
  pub(crate) fn entry(&self, id: DativeBondId) -> DativeBondEntry<'_>;
  ```

  FrameTransport::Action is DynPermutation on donors; the acceptor is fixed.
  Reuse dative_bond_representative_action. Test reordered donors, changed donor
  membership, different acceptors, and normalized attribute equality. Transfer
  the corresponding Molecule comparison cases to entry-level tests while keeping
  the old callers until S4d7.

  Before implementing entry comparisons, extend benches/editor.rs with focused
  overlay-removal cases through the public Edits API, covering all six kinds in
  stored and reordered frames. Use normal fixture/setup functions and capture
  the current comparison cost for S4d7. Keep this bounded to the replacement of
  these comparisons; do not start another general editor benchmark study.

  **Verification — 2026-09-28.** The 18 entry/getter cases pass, including
  ownership, frame transport, reframing, and the transferred comparison cases.

  **Removal baseline — 2026-09-28.** benches/editor.rs now covers all six
  kinds in stored and reordered frames. Each fixture has five atoms, four
  localized bonds, and one of each overlay. Molecule construction and Edit
  preparation are outside timing; the timed operation is editor.apply of one
  removal, without publication. Criterion used 20 samples, 0.2 s warm-up and
  0.5 s measurement per case:

  | Removed entity | Stored frame, ns | Reordered frame, ns |
  | --- | ---: | ---: |
  | Dative bond | 739 | 717 |
  | Aromatic system | 794 | 783 |
  | Multicenter bond | 772 | 771 |
  | Noncovalent bond | 595 | 605 |
  | Stereo atom | 961 | 967 |
  | Stereo bond | 952 | 995 |

  Command: `cargo bench -p umol-graph-ir --bench editor --
  molecule_editor/removal --warm-up-time 0.2 --measurement-time 0.5
  --sample-size 20 --save-baseline s4d-before`.
  These are whole-removal estimates before caller migration, not isolated
  comparison timings. They let S4d7 assess the cost of using entries in the
  actual caller. S4d1–S4d6 do not change that caller, so no second timing run
  is needed at this boundary.

- **S4d2 — completed 2026-09-28** (`ir::aromatic`; additive, green). [dep: S1a]

  ```rust
  pub(crate) struct AromaticSystemEntry<'a> {
      pub(crate) atoms: Cow<'a, [AtomId]>,
      pub(crate) attributes: Cow<'a, AromaticSystemForm>,
  }
  // On AromaticSystems:
  pub(crate) fn entry(&self, id: AromaticSystemId) -> AromaticSystemEntry<'_>;
  ```

  FrameTransport::Action is DynPermutation on atoms. Reuse
  aromatic_system_representative_action; electron contributions move with their
  atoms through the existing form transport. Test equivalent reordered
  atom/count pairs, reordered atoms with unchanged counts, different atom sets,
  and undetermined contributions. Transfer the existing aromatic comparison
  cases. Preserve the established first-use behavior for mismatched count lengths.

  **Verification — 2026-09-28.** The 24 entry/getter cases pass, including
  ownership, frame transport, reframing, and the transferred comparison cases.

- **S4d3 — completed 2026-09-28** (`ir::multicenter`; additive, green). [dep: S1a]

  ```rust
  pub(crate) struct MulticenterBondEntry<'a> {
      pub(crate) atoms: Cow<'a, [AtomId]>,
      pub(crate) attributes: Cow<'a, MulticenterBondForm>,
  }
  // On MulticenterBonds:
  pub(crate) fn entry(&self, id: MulticenterBondId) -> MulticenterBondEntry<'_>;
  ```

  FrameTransport::Action is DynPermutation on atoms. Reuse
  multicenter_bond_representative_action and the existing form transport.
  Transfer the multicenter comparison cases and cover atom/count alignment,
  undetermined contributions, changed membership, and incompatible action degree.
  This has the same count-length contract as S4d2, without sharing entity storage
  through a new generic entry family.

  **Verification — 2026-09-28.** The 24 entry/getter cases pass, including
  ownership, frame transport, reframing, and the transferred comparison cases.

- **S4d4 — completed 2026-09-28** (`ir::noncovalent`; additive, green). [dep: S1b]

  ```rust
  pub(crate) struct NoncovalentBondEntry<'a> {
      pub(crate) atoms: [AtomId; 2],
      pub(crate) attributes: Cow<'a, NoncovalentBondForm>,
  }
  // On NoncovalentBonds:
  pub(crate) fn entry(&self, id: NoncovalentBondId) -> NoncovalentBondEntry<'_>;
  ```

  FrameTransport::Action is a degree-two DynPermutation. Reuse
  noncovalent_bond_representative_action and transport the full form. Transfer
  the stored/reversed/different-pair cases; include changed kind and inline
  constraints. Getter construction allocates neither an atom vector nor a form.

  **Verification — 2026-09-28.** The 13 entry/getter cases pass, including
  ownership, frame transport, reframing, and the transferred comparison cases.

- **S4d5 — completed 2026-09-28** (`ir::stereo`; additive, green). [dep: S1c]

  ```rust
  pub(crate) struct StereoAtomEntry<'a> {
      pub(crate) site: AtomId,
      pub(crate) ligands: Cow<'a, [StereoLigand]>,
      pub(crate) attributes: Cow<'a, StereoAtomForm>,
  }
  // On StereoAtoms:
  pub(crate) fn entry(&self, id: StereoAtomId) -> StereoAtomEntry<'_>;
  ```

  FrameTransport::Action is Permutation on complete StereoLigand values; the site
  is fixed. Reuse stereo_atom_representative_action and the full form transport,
  including configuration and inline constraints. Transfer the stereo-atom frame
  comparison cases and cover changed site, virtual ligands, and frame-relative
  constraints. Entry access borrows both ligands and attributes.

  representative_action retains Reframe's admissible-frame precondition. The
  fallible reframe path uses the existing Option-returning action derivation;
  comparison must not turn an oversized supplied frame into an assertion panic.
  This adds no coset checks to normalization and no new error type.

  **Verification — 2026-09-28.** The 21 entry/getter cases pass, including
  ownership, frame transport, reframing, and the transferred comparison cases.

- **S4d6 — completed 2026-09-28** (`ir::stereo`; additive, green). [dep: S1c]

  ```rust
  pub(crate) struct StereoBondEntry<'a> {
      pub(crate) site: BondId,
      pub(crate) ligands: Cow<'a, [StereoLigand]>,
      pub(crate) attributes: Cow<'a, StereoBondForm>,
  }
  // On StereoBonds:
  pub(crate) fn entry(&self, id: StereoBondId) -> StereoBondEntry<'_>;
  ```

  FrameTransport::Action is Permutation restricted to the existing stereo-bond
  endpoint-block group; the site bond is fixed. Reuse
  stereo_bond_representative_action and the full form transport. Transfer the
  within-block, complete-block-swap, across-block, changed-ligand, and changed-site
  comparison cases. Include determined and undetermined configurations and
  frame-relative constraints. Getter ownership and malformed-frame handling
  follow S4d5; no copy is needed to obtain the ligand slice.

  **Verification — 2026-09-28.** The 33 entry/getter cases pass, including
  ownership, frame transport, reframing, and the transferred comparison cases.

  **Combined S4d1–S4d6 checks.** All 133 new entry/getter cases pass;
  the combined `cargo test -p umol-graph-ir --lib entry` run passes 229 cases.
  The five affected entity-module suites passed during implementation.
  Private-item rustdoc with warnings denied, nightly formatting, and diff checks
  pass. Clippy over library, tests, and benches reports only dead_code: the new
  entries/getters await S4d7, and seven existing Molecule methods await their
  planned callers. No lint suppression was added; the strict gate remains S4d7.
  All six types, fields, and entry getters are pub(crate), as specified.
  Dative/aromatic/multicenter getters allocate their atom-id vectors;
  noncovalent/stereo getters borrow attributes without allocation, and stereo
  getters also borrow ligand slices.

- **S4d7 — completed 2026-09-28** (`ir::molecule::apply`,
  entity tests, `benches/editor.rs`; rewire, green).
  [dep: S4b9, S4d1, S4d2, S4d3, S4d4, S4d5, S4d6]

  In apply_edit and apply_edit_with_undo, obtain the stored entry from the owning
  set, construct the old entry from the already-resolved ids and attributes,
  and use framed_eq. Remove all six Molecule *_equiv methods and their imports;
  S4d1–S4d6 have transferred their tests to the entity modules. Retain no
  delegating comparison wrappers. This does not change Edit/Undo payloads or
  mutation sequencing.

  Preserve removal's separate structured-incidence precondition. Ordinary
  overlays use the existing set is_coincident methods. Stereo entries use site
  equality and Permutation::between on complete ligand values; stereo bonds also
  require that permutation to belong to the existing endpoint-block group.
  This matters when both forms are contradictory: the trait equates two
  contradictions, but removal must still refer to the supplied atoms/site/ligands.
  Keep these preconditions in the Edit branches before writes; do not modify
  framed_eq's contract or introduce another comparison method to conceal them.

  Run the existing removal/old-state and matching-history rollback cases in both
  execution paths, plus cases with contradictory attributes and mismatched
  structured incidence. Run the affected public Edit/reframe properties. Compare
  the same removal benchmarks captured in S4d1 and record the time/allocation
  consequences of entry construction and reframing; the wrapper alone does not
  establish that comparison became cheaper. Update the frame-carrier description
  in the data-type guide, then run affected-crate checks, strict Clippy, rustdoc,
  nightly formatting, and full diff review. Workspace/Python/MSRV gates remain S9b.

  **Implementation.** Both Edit execution paths
  construct old entries and call framed_eq after the incidence check. The six
  Molecule comparison methods are removed. Each entry specializes framed_eq:
  structural equality returns immediately; otherwise it transports one form
  into the other's frame and compares normalized attributes. If this does not
  establish equality, it compares the two reframing results, preserving the
  trait's treatment of contradictions. No comparison wrapper or public API is
  added. Matching-history undo still restores the actual stored frame.

  The default framed_eq implementation doubled reordered stereo removal time by
  reframing both complete entries. The specialization removes that work. Stereo
  incidence uses the derived permutation directly, avoiding the four temporary
  vectors used by the general relation-set coincidence query.

  Same S4d1 fixtures and Criterion settings (0.2 s warm-up, 0.5 s measurement,
  20 samples; saved baseline s4d-before). Central estimates, rounded to ns:

  | Removed entity | Baseline stored | Current stored | Baseline reordered | Current reordered |
  | --- | ---: | ---: | ---: | ---: |
  | Dative bond | 739 | 718 | 717 | 746 |
  | Aromatic system | 794 | 711 | 783 | 818 |
  | Multicenter bond | 772 | 708 | 771 | 834 |
  | Noncovalent bond | 595 | 560 | 605 | 611 |
  | Stereo atom | 961 | 651 | 967 | 1,014 |
  | Stereo bond | 952 | 699 | 995 | 1,053 |

  Retain the specialization: it removes the large regression while keeping the
  comparison in the owning entity type. Stored-frame removals improve; reordered
  removals retain a small cost (about 1–8% in these central estimates). This is
  not a claim that the migration improves every path. No further timing study
  is needed for this decision.

  Allocation accounting from the implementations, not allocator instrumentation:
  entry retrieval allocates one atom-id vector for dative/aromatic/multicenter
  entries and none for noncovalent/stereo entries. Successful alignment creates
  no transported frame vectors and transports only one attribute form. Dynamic
  alignment allocates its image and used-position vectors; bounded stereo
  alignment allocates one image vector. Stereo removal derives alignment once
  for incidence and, when structural equality does not apply, once for framed_eq.
  This second derivation is one additional small allocation over the baseline.
  Attribute allocations depend on the form; the general definition remains the
  fallback when the direct comparison does not establish equality.

  **Verification.** All 2,177 entry and molecule unit cases pass, including 764
  new definition-comparison cases and 24 removal-incidence cases. The affected
  public Edit/frame/reframe property run passes all 52 properties after S4b8's
  replay correction. S4b9 closes the strict all-target Clippy gate and verifies
  the non-test build. Private-item rustdoc with warnings denied, nightly
  formatting, and full diff review pass. No lint suppression was added.
  S4b8 records the unrelated-history property's corrected scope.

### S5 — Borrowed transaction API

S5a removes APIs used by S5b–S5d; those migrations are required before the stage
returns green. S5d's Python invalidation and sequential input-consumption contracts are recorded above.

- **S5a — Scoped Rust transaction lifecycle** (group; breaking, green at
  S5d3). [dep: S4b, S4d7]

- **S5a1 — completed 2026-09-28** — Guard, borrowed handle, and scoped run
  (`ir::molecule::transact`; breaking, green at S5d3). [dep: S4b, S4d7]

  ```rust
  struct TransactionGuard<'a> {
      molecule: &'a mut Molecule,
      journal: Vec<Undo>,
      status: TransactionStatus,
  }

  #[derive(Clone, Copy, PartialEq, Eq)]
  enum TransactionStatus {
      Active,
      CommitRequested,
      Accepted,
      RolledBack,
      Aborted,
  }
  ```

  Create molecule/transact.rs for the borrowed Transaction lifecycle, independent
  of editor. Retire the detached Transaction implementation in editor::transact;
  editor apply/tracked_apply stay under editor. Update the Transaction re-export
  to the new module. Transaction methods call Molecule's shared execution methods
  in apply and never access editor fields.

  Both types and all fields are private to transact. Construct the guard directly
  in Transaction::run with an empty journal and Active status; no public
  constructor, separate journal object, or recovery clone. The Transaction
  handle borrows its three fields separately, using the single-lifetime shape
  at the start of this document. It does not own or borrow the whole guard.

  The only guard method is Drop::drop(&mut self): for Active or CommitRequested,
  pop the journal in reverse and call molecule.apply_undo for each entry. The
  other statuses require no restoration. Do not call forward Edit execution,
  validate undo, update correspondence, or install an empty molecule. These
  operations rely on the completed-Edit journal contract, not partial-write
  recovery. No per-apply guard or catch_unwind inside Edit execution.

  Transaction::run creates the guard on its stack and lends its molecule,
  journal, and status fields to the callback handle. Public methods use only
  these borrows. No extra lifetime, RefCell, TLS, unsafe code, or dependency.

  | Event | State and result |
  | --- | --- |
  | Successful apply | Append completed undos; remain Active |
  | Failed apply | Reverse and drain the journal; set Aborted; return the original error |
  | probe | Check integrity and return the borrow or error; do not change status |
  | Successful commit/tracked_commit | Check integrity; set CommitRequested; keep the journal until run decides acceptance |
  | Failed commit | Reverse and drain the journal; set Aborted; return Integrity(error) |
  | Explicit rollback | Reverse and drain the journal; set RolledBack unless already Aborted, which remains latched |
  | Callback Ok with CommitRequested | Set Accepted, then return the callback value |
  | Callback Ok with Aborted | Return TransactionError::Aborted through MoleculeApplyError and E |
  | Callback Ok without commit | Restore through guard destruction; return the callback value |
  | Callback Err or unwind | Preserve the callback error/unwind; guard restores any unaccepted completed changes |

  apply/commit on an Aborted handle return Aborted without writing. probe may
  read the restored molecule. Consuming commit/rollback prevent later use of
  those handles. Forgetting a handle leaves the guard owned by run. Test that
  case, callback unwinding between completed edits and after a commit request,
  ignored apply errors (including explicit rollback afterward), explicit
  rollback, and integrity rejection. These are
  the integration cases for this scoped API.

  **Implementation and verification.** The guard, borrowed handle, and run live
  in molecule::transact. TransactionError moves there and gains Aborted. The
  detached journal and editor transaction entry points are removed; editor
  apply/tracked_apply retain their implementations. S5a2 exercises the guard
  through public lifecycle methods. The stage gate remains S5d3.

- **S5a2 — completed 2026-09-28** — Batch application, completion, and Molecule conveniences
  (`ir::molecule::transact`, `ir::error`; breaking, green at S5d3). [dep: S5a1]

  Complete the detached Transaction replacement with immediate `apply`, checked
  `probe`, `commit`, `tracked_commit`, and `rollback`; add Molecule's prepared-batch
  `transact` and `tracked_transact`. Retire detached-journal `undos`, `append`,
  `rollback`, and `tracked_rollback`, plus editor `transact`/`tracked_transact`,
  without adding Transaction::tracked_apply. Remove obsolete rollback errors;
  S4b8 already removes validate_undo. Retain the drafted MoleculeApplyError
  return type. Test separate batch handle namespaces, whole-transaction
  correspondence, abort latching, callback value return after no-commit
  rollback, and checked acceptance.

  Implement the exact public signatures at the start of the document. Each
  apply call creates fresh ApplicationState and records only completed edits;
  commit checks the same integrity gate as construction. tracked_commit derives
  whole-transaction correspondence from the retained journal only when requested.
  No unconditionally maintained correspondence field is added to the guard.
  Prepared Molecule calls apply each batch, then commit once. Verify matching
  ordinary/tracked results, separate namespaces, and the S5a1 lifecycle table.

  **Implementation and verification.** All seven methods are implemented in
  molecule::transact. Application and commit errors drain the full journal and
  latch Aborted. probe and commit call Molecule::check_integrity. Ordinary commit
  constructs no correspondence; tracked_commit derives starting counts by reading
  the journal backward, then applies its recorded additions and compactions to
  the correspondence in execution order. It does not replay edits or mutate the
  molecule. RollbackFailed and RollbackStateMismatch are removed; Undo rustdoc and
  the data-type/nomenclature guides describe the scoped lifecycle.

  The 32 public-API unit cases cover callback return/error/unwind, forgotten
  handles, abort latching, separate batch namespaces, probe repair, constructor-
  equivalent integrity rejection, explicit rollback, and correspondence across
  topology and all six overlay kinds. All pass with the S5b caller migration.
  The workspace stage gate remains S5d3.

- **S5b — completed 2026-09-28** (`umol-graph-ir` reaction and molecule callers; breaking, green at S5d)
  Migrate uses of detached journals. Preserve reaction `Ok(None)` versus error
  classification and host-to-product correspondence.
  Test failed applications and product integrity. [dep: S5a]

  **Implementation and verification.** Reaction application uses editor.apply on
  its product candidate, retaining the existing integrity-error classification
  and host-to-product correspondence. Transaction callers and benchmarks use the
  scoped API. Batch/rollback unit cases live under molecule::transact; manipulated
  Undo cases call Molecule::apply_undo internally. All 135 cases from the former
  editor transaction module are retained.

  Public properties cover batch namespaces, rollback, and correspondence
  composition. The detached-journal concatenation and unrelated-receiver
  properties are retired with that public interface; multi-batch rollback and
  internal no-panic Undo cases cover the retained contracts. Constraint-compaction
  properties read transient constraints through the editor and separately verify
  transaction rollback, without requiring transient stereo incidence to pass probe.

  Graph-ir verification passes: 8,775 library tests (3 ignored); 44 selected
  edit, molecule-compaction, and reaction-application property-suite tests with
  PROPTEST_CASES=256; doctests (1 passed, 3 ignored); benchmark compilation;
  all-target strict Clippy with proptest; warnings-denied rustdoc; nightly
  formatting and diff review. The workspace gate remains S5d3.

- **S5c — completed 2026-09-28** (`umol-graph::ops`; breaking, green at S5d) Migrate existing borrowed
  resolve/project execution to scoped transactions, sharing the planned
  batches and using checked probes between phases. Keep their present public
  resolve/project signatures in this stage. Extract these three approved
  methods before composite projection uses them:

  ```rust
  StereoResolver::plan_project(&self, molecule: &Molecule)
      -> Result<Edits, StereoProjectError>;
  AromaticityResolver::plan_project(&self, molecule: &Molecule)
      -> Result<Edits, AromaticityProjectError>;
  IsotopeResolver::plan_project(&self, molecule: &Molecule)
      -> Result<Edits, IsotopeProjectError>;
  ```

  These use the existing error types; no shared phase-error type is introduced.
  Do not nest separately committed phase transactions. Extract shared composite
  chemistry into these accepted private Resolver methods in resolve.rs:

  ```rust
  fn select_atom_completions(&self, state: &mut ResolveState);
  fn plan_constitution(&self, molecule: &Molecule, state: &ResolveState) -> Edits;
  ```

  select_atom_completions applies the current final atom tie-break: select a
  unique best completion, retain unresolved ties, and record the resolved atom
  ids in state.tie_breaks. Before plan_constitution, retain the existing
  rejection of plural completion entries. plan_constitution contains the current
  atom-field and aromatic-system batch construction, preserving stored
  constraints and duplicate-system handling. Neither routine executes edits,
  owns an editor/transaction, or constructs a ResolveReport. Keep plan_placement
  and plan_discharge unchanged. No generic phase runner is introduced.
  Wire these concrete error variants at this boundary; each Apply variant
  carries MoleculeApplyError and provides From<MoleculeApplyError> for scoped
  execution and commit:

  | Error enum | Change |
  | --- | --- |
  | ResolveError | Add Apply(MoleculeApplyError) for probe/commit; change existing Isotope, Commit, Placement, and Discharge payloads from TransactionError to MoleculeApplyError, preserving their phase attribution |
  | ProjectError | Add Apply(MoleculeApplyError) |
  | IsotopeError | Replace Transaction(TransactionError) with Apply(MoleculeApplyError) |
  | IsotopeProjectError | Replace Transaction(TransactionError) with Apply(MoleculeApplyError) |
  | AromaticityError (ops::aromaticity) | Replace Transaction(TransactionError) with Apply(MoleculeApplyError) |
  | AromaticityProjectError | Replace Transaction(TransactionError) with Apply(MoleculeApplyError) |
  | StereoError | Replace Transaction(TransactionError) with Apply(MoleculeApplyError) |
  | StereoProjectError | Replace Transaction(TransactionError) with Apply(MoleculeApplyError) |
  | BondsError | Replace Transaction(TransactionError) with Apply(MoleculeApplyError) |
  | MulticenterBondsError | Replace Transaction(TransactionError) with Apply(MoleculeApplyError) |

  Preserve the existing chemistry variants, contradictions, and aggregate
  chemistry-error wrappers. Migrate affected error conversions, matches, Python
  bindings, and tests in the S5 migration. ValenceError and ValenceProjectError
  gain no Apply variant: valence admission produces no edits, and its current
  projection does not mutate. No Phase error enum or placeholder is introduced.
  S7 changes ownership names and reporting, not this prerequisite.
  Test late rejection restores the entry molecule and preserves chemistry
  diagnostics. [dep: S5a]

  **Implemented.** Composite resolve and project each use one Transaction::run,
  checked intermediate probes, and one commit on a determined result. Other
  outcomes restore the entry molecule. Standalone phases use Molecule::transact;
  projection planning uses the three methods above. Resolver and projector
  execution no longer copies whole molecules. Valence projection remains a no-op.
  The private constitution methods and all specified error payloads are wired.
  Python has no matches on the changed chemistry-error variants to migrate;
  its transaction lifecycle migration remains S5d.

  Added projection-planner cases and extended late projection rejection coverage
  to include restoration of removed overlays and their molecule constraints.
  Verification passes: 1,968 graph library tests, 168 property tests with
  PROPTEST_CASES=256, 683 resolution-conformance cases, all-target strict Clippy
  with both features, warnings-denied rustdoc, nightly formatting, and diff review.
  The workspace gate remains S5d3.

- **S5d — Python prepared transactions and ownership** (group; breaking,
  closes S5 green at S5d3). [dep: S5a, S5b, S5c]

  Execute S5d1–S5d3 in order. No interactive Python transaction handle is added.

- **S5d1 — completed 2026-09-28: Edits consumption and iterator access** (`umol-py::edit` and its consumers; breaking, green at S5d3). [dep: S5a, S5b, S5c]

  **Edits ownership interface.** Replace the current tuple wrapper with:

  ```rust
  #[pyclass]
  pub struct Edits {
      value: Option<GraphIrEdits>,
  }

  impl Edits {
      pub(crate) fn from_rust(value: GraphIrEdits) -> Self;
      pub(crate) fn to_rust(&self) -> PyResult<&GraphIrEdits>;
      pub(crate) fn to_rust_mut(&mut self) -> PyResult<&mut GraphIrEdits>;
      pub(crate) fn take(&mut self) -> PyResult<GraphIrEdits>;
  }
  ```

  from_rust stores Some. Accessors and take report ConsumedError on None; take
  moves the value and leaves None. There is no refill operation or automatic
  clone. Keep the existing append/add/update construction methods, making their
  access fallible; unrelated forms retain their current ownership. Replace
  derived Python equality with value comparison through the fallible accessors;
  two consumed wrappers do not compare equal by virtue of both storing None.
  Indexing, length, repr/render, conversions, and iteration also check access.
  Current yielded Edit values are independent copies and remain usable.

  Keep EditIter's existing owner: Py<Edits>, position: usize, and end: usize.
  __next__(&mut self, py: Python<'_>) -> PyResult<Option<Py<Edit>>> borrows the
  owner and checks for a present value before testing exhaustion. A consumed
  owner raises InvalidatedViewError, including on an already-exhausted iterator.
  No counter is needed: consumed Edits is never refilled. Do not widen this work
  into a new Edit/view representation or change existing entry-copy semantics.

  Migrate all to_rust users to fallible access and use take at consuming Rust
  Edits boundaries, including editor apply. Do not preserve an implicit batch
  clone until S6. Test direct aliases, copying by explicit user action where
  already supported, exhausted/live iterators, equality, and consumed access.

  **Implemented.** Edits uses the Option storage and four access/transfer methods
  above. Existing Molecule/editor batch consumers take their Edits input.
  All 34 construction mutators check access, as do reads, rendering, equality,
  and iterator creation. Repr uses the same entry-list format as Deltas.
  EditIter checks owner availability before exhaustion; yielded entries retain
  their independent-copy semantics. Tests explicitly copy a batch when reusing it.

  Added Rust transfer/iterator tests and Python cases for aliases, successful
  and failed consumption, explicit copies, all mutators, equality, and live,
  partial, empty, and exhausted iterators. Nightly formatting, Python syntax
  checks, and diff review pass. Python 3.13 all-target Cargo checking reports
  the same six errors as before this subitem: the detached Transaction lifetime,
  removed editor transact methods, and retired rollback-error variants.
  Runtime verification is recorded under S5d3.

- **S5d2 — completed 2026-09-28: Molecule accessor counters and Storage names** (`umol-py::molecule`, entity/constraint/collection views; breaking, green at S5d3). [dep: S5d1]

  Implement the [approved accessor invalidation contract](#python-accessor-invalidation).
  Current AtomView stores an owner and a dense id; deleting an earlier atom can make it
  silently address another atom. Owner consumption is not involved in transact,
  so S6d's consumed-state checks alone cannot fix this. Apply the specified
  triggers, failure behavior, nested access, and *Storage names, including the
  later *_into consumers. Whole-molecule mutation invalidates all owner-backed
  accessors at execution entry, including on no-op and rollback; ordinary view
  setters retain access.

  Reuse the Molecule counter methods and molecule constraint accessors from S2l.
  Apply the remaining field and method table under Python accessor invalidation.
  Test all eight entity kinds, nested constraints and ring sizes, exhausted
  iterators, and ordinary setter usability. Molecule storage is not Option yet;
  S6d adds consumption without changing these counter rules.

  All eight entity views, collections, and owner-backed iterators now check their
  captured counter under the storage borrow. Constraint and ring-size accessors
  inherit it; stale parents cannot create fresh descendants. Id, repr, lookup,
  iteration (including exhaustion), assignment, and conversion paths check access.
  The eleven Backing enums and their fields/macro parameters use the approved
  Storage names. Ordinary setters do not advance the counter. Existing owned
  constraint-entry and ring-size iterators retain their copying semantics.

  Python tests cover all eight entity kinds, empty/nonempty combine_from,
  nested ring-size access, live/exhausted iterators, rejected arguments,
  independent values, and ordinary setters. Rust tests use the fallible binding
  accessors. Nightly formatting, Python syntax checking, and diff review pass.
  Python 3.13 all-target Cargo checking reports only the same six existing
  transaction-binding/test errors recorded in S5d1. Runtime verification,
  including transaction invalidation and rollback, is recorded in S5d3. Later *_into
  bindings must advance the same counter at their execution boundary as specified
  above; those bindings are not present yet.

- **S5d3 — completed 2026-09-28: Prepared-batch transaction bindings** (`umol-py::{molecule,transaction}`; breaking, red→green). [dep: S5d2]

  Transfer prepared Edits batches sequentially in input order, with no
  availability precheck, duplicate-object scan, or input recovery. Each transfer
  takes the stored Rust value or raises ConsumedError if it is already absent.
  If a later transfer fails, earlier batches remain consumed and their transferred
  values are dropped; batches not yet visited remain untouched. Thus
  [first, already_consumed] consumes first before failing, and [edits, edits]
  consumes the object on its first occurrence and fails on its second.
  The molecule and its accessor counter remain unchanged: execution starts only
  after all transfers succeed. No implicit clone substitutes for a second
  transfer. Test both failure cases and their input-consumption effects.
  Test Edits iterators retained across consumption as well as direct aliases.

  The Python entry signatures are:

  ```rust
  fn transact(slf: Py<Self>, py: Python<'_>, batches: &Bound<'_, PyAny>)
      -> PyResult<()>;
  fn tracked_transact(slf: Py<Self>, py: Python<'_>, batches: &Bound<'_, PyAny>)
      -> PyResult<MoleculeCorrespondence>;
  ```

  These are methods of the Python Molecule wrapper. Accept any Python iterable
  and collect its Edits objects into Vec<Py<Edits>> before taking their values;
  iteration or element-type errors therefore consume no batches. PyO3's direct
  Vec extraction only accepts sequences, so the binding performs this iteration
  explicitly. Use a short try_borrow_mut for each Edits transfer; after all
  transfers, exclusively borrow the molecule, advance its
  counter, and call the Rust operation. A later transfer, receiver-borrow, or
  counter-overflow error does not restore already consumed batches. No Python
  iteration/callback occurs during Rust application. The returned correspondence
  is the existing Python wrapper, constructed from the Rust result.

  Remove the detached Python Transaction class and editor transact methods.
  Test independent New namespaces, no-op invalidation, rollback after a later
  failure, tracked correspondence, and no invalidation before execution starts.

  **Implemented.** Molecule.transact and tracked_transact accept arbitrary
  iterables, transfer batches sequentially, then advance the receiver counter
  and call the Rust operation under one exclusive borrow. The detached Python
  Transaction class and editor transaction methods are removed, including their
  exports and tests. Transaction tests now exercise the Molecule entry points.

  Tests cover independent New namespaces, commit-only integrity checking,
  later-batch and integrity-error rollback, participant replacement and new
  handles, composed correspondence, consumed/aliased batches, iterator failures,
  and receiver/batch borrow failures. All eight entity-view families and nested
  accessors invalidate on empty transactions; field changes and rollback also
  invalidate views. Extraction and transfer failures preserve receiver access.
  S5d2's stereo constraint tests use iteration; those views have no asdict method.

  **Verification.** Workspace tests pass: 20,498 passed, 9 ignored, including
  doctests. After the iterable-input adjustment, all 1,639 Python-binding Rust
  tests pass again (2 ignored). The rebuilt Python 3.13 extension passes 1,808
  pytest cases (2 skipped), including S5d1/S5d2 runtime coverage. Workspace
  all-target strict Clippy, warnings-denied graph-ir/graph/Python rustdoc,
  nightly formatting, and full diff review pass. S5b/S5c record the feature-gated
  property and conformance checks. S6a is next.

### S6 — Owning editor and publication

S6a–S6d are one ownership migration, returning green after S6d. Migrate every
consumer of the changed signatures, including tests and benchmarks; retaining
temporary cloning adapters is not a way to close an earlier subitem.

- **S6a — completed 2026-09-28** (`ir::molecule`, `ir::molecule::editor`; breaking, green at S6d) Make
  `edit(self)` own its Molecule, make Molecule `apply`/`tracked_apply` consume,
  and replace editor snapshot/build methods with checked `probe`/`finish`.
  Retain MoleculeBuilder's asserted build and the S2 uniform mutable attribute access.
  Move the owned constraint collection with the molecule; no clone is needed.
  Test direct/batch interleaving, invalid probe then repair, destructive failure, and
  constructor-equivalent publication.
  [dep: S2m, S5d]

  **Implemented.** Molecule.edit moves every storage field, including Constraints,
  into the editor. Molecule apply/tracked_apply consume and finish; the tracked
  form obtains its correspondence from the editor's tracked_apply. Editor batch
  methods return MoleculeApplyError. They retain transient state until probe or
  finish invokes the same check_integrity used by construction. The six editor
  snapshot/build methods are removed. MoleculeBuilder.build calls finish and
  asserts its producer contract. Session correspondence remains until S6b1.

  Publication tests now cover borrowed probes, repair after rejection, constructor
  agreement, and moving the atom/bond/constraint allocations through editing.
  Added direct/batch interleaving and consuming-failure cases; updated molecule
  application tests to the consuming contract.

  **Verification boundary.** Nightly formatting and diff review pass.
  cargo check -p umol-graph-ir --lib reports seven caller-migration errors:
  fragment/extract/combine/split still call editor build, and reaction application
  still uses try_build and the previous editor error type. These belong to
  S6b1/S6b2. The added tests cannot execute until the graph-ir callers compile;
  no compatibility adapters were added. The workspace returns green at S6d.

- **S6b — Graph-IR ownership migration** (group; breaking, green at S6d).
  [dep: S6a]

- **S6b1 — completed 2026-09-28 — Graph-IR callers and correspondence** (`umol-graph-ir` reaction/molecule callers; breaking, green at S6d). [dep: S6a]

  Migrate editor construction/publication and remove session
  correspondence accumulation and public tracked direct removal. Use one
  source-preserving product candidate where the host survives; keep batch-only
  tracking distinct from whole-transaction tracking. Test product and
  correspondence laws.
  Keep reaction product failure classification and one intentional host copy.
  A rejected product candidate is dropped; it needs no recovery transaction.

  **Implemented.** MoleculeEditor holds only its Molecule. Direct additions,
  removals, and plain apply no longer accumulate correspondence. tracked_apply
  constructs one mapping for its own batch. Tracked direct removal remains
  crate-private for extraction and internal compaction; no public tracked direct
  removal is exposed. Reaction application uses one explicit host clone and
  consuming Molecule.apply, retaining product-conflict rejection and execution
  errors. Fragment, extraction, split, tests, and benchmarks use probe/finish.
  combine_from's publication call is migrated; its storage rewrite remains S6b2.

  Tests retain direct mutation, cascading removal, restoration, and publication
  assertions. Retired session-mapping assertions are replaced by batch identity,
  batch composition versus transaction correspondence, and mixed remapping /
  compaction / batch correspondence checks. The retired asserted editor-build
  panic test is removed; finish/probe rejection and repair remain covered.
  Benchmark consuming inputs are prepared outside timed apply/tracked_apply;
  editor_session retains its source-preserving copy inside the timed operation.

  **Verification.** Graph-ir library/tests/benches compile with proptest enabled.
  Focused nextest selection: 3,119 passed, covering molecule/editor/transaction,
  reaction, correspondence, affected views and construction tests, integrity,
  and edit/publication/reaction properties (64 generated cases per property).
  Strict graph-ir all-target Clippy with proptest, rustdoc, nightly formatting,
  and diff checks pass. No workspace or MSRV gate was run for this subitem.

- **S6b2 — completed 2026-09-28 — Disjoint append on a borrowed molecule** (`ir::molecule`; breaking rewire, green at S6d). [dep: S6b1, S4b1]

  Migrate combine_from's current mem::take/from_parts
  route here because the owning editor change affects it now, not in S8.
  Preserve its append semantics and borrowed receiver without an initial recovery
  clone. No append/count checkpoints or recovery tests for internal panics between
  storage writes are required. This moves the former S8c migration here. The
  separate lift_constraints correction remains excluded.

  Keep combine_from(&mut self, other: &Molecule) -> () and its existing
  disjoint-append semantics. Save the eight original entity counts solely as
  offsets for other, then use crate-private Molecule additions directly: atoms,
  bonds with shifted endpoints, and the six overlay kinds with shifted sites
  and atoms. Retain the existing correspondence-based ligand and constraint
  translation; append mapped constraints through extend_constraints. Only
  other's borrowed payloads are copied. Do not move self into an editor, call
  mem::take(self), clone self for recovery, or retain offsets as rollback
  checkpoints. Finish with the existing asserted producer integrity guarantee;
  no new public method or error type is needed.

  Verify original-id prefixes, every overlay factor, inline/top-level constraints,
  empty operands, shared source storage, and the unchanged other molecule.

  **Implemented.** combine_from calls Molecule's add methods for all eight
  entity kinds and extend_constraints for mapped constraints. It retains the
  existing offset correspondences and checks integrity directly after appending.
  The receiver stays in place throughout; no recovery copy or journal is added.

  **Verification.** All 13 focused unit cases pass, covering exact appended
  entries, every entity kind's constraint references, inline constraints, empty
  operands, shared storage with an unchanged source, and unique attribute-storage
  reuse. Four combination properties pass with 256 cases each. Strict graph-ir
  library/test Clippy with proptest, nightly formatting, and diff checks pass.
  The S6 workspace gate remains at S6d.

- **S6c — Chemistry and format caller migration** (group; breaking, green
  at S6d). [dep: S6a, S5c]

- **S6c1 — completed 2026-09-28 — Transformation plans and borrowed execution** (`umol-graph::ops`; breaking rewire, green at S6d). [dep: S6a, S5c]

  Migrate graph operation publication sites to `finish` or borrowed transactions
  without changing
  their chemistry or boundary outcomes. This includes AromaticityPerceiver's
  add_systems and the three existing transform_into implementations. Construct
  their batches and shared chemistry checks according to the transformation
  mapping now: borrowed callers cannot simply consume their receiver.
  Preserve their current public signatures and error contracts. S8 reuses these
  plans for the consuming transformation and lazy iterator; it does not supply a
  missing prerequisite for S6. Test published outputs, late rejection recovery,
  and DelocalizeCharge's asserted producer contract.

  The accepted shared single-operation plans are private inherent methods:

  ```rust
  Aromatizer::plan_transform(&self, molecule: &Molecule) -> Result<Edits, AromatizeError>;
  DelocalizeCharge::plan_transform(&self, molecule: &Molecule) -> Edits;
  Kekulizer::plan_transform(&self, molecule: &Molecule) -> Result<Edits, KekulizeError>;
  ```

  These contain only the existing chemistry planning plus Edit construction
  specified in the transformation mapping. Reuse DelocalizationPlan::derive
  and Kekulizer::plan_systems; keep validate_localized_candidate as the shared
  kekulization postcheck. Plans do not execute, clone the whole molecule, or
  introduce a new plan type. Borrowed execution applies the batch in one
  transaction; consuming execution in S8 applies it in one editor. The
  existing AromaticityPerceiver::add_systems signature remains unchanged and
  performs its own Edits construction and transaction, without a new public
  planning method or nesting inside an Aromatizer transaction.

  Verify each existing transformation and AromaticityPerceiver::add_systems,
  including unchanged-input and late-rejection cases.

  **Implemented.** The three private plan_transform methods build Edits from
  the existing chemistry plans. transform_into executes one transaction;
  Kekulizer probes for its valence/spin checks and commits only on acceptance.
  Its recovery candidate copy is removed. DelocalizeCharge retains Infallible
  and asserts its execution contract. AromaticityPerceiver::add_systems batches
  additions and bond assertions independently. Existing error types and public
  signatures remain unchanged. Graph-operation test publication uses finish.

  Planning relies on aromatic systems' disjoint membership guaranteed by
  molecule integrity. It emits edits directly per system. The additional work
  is Edit construction and undo recording for affected fields/entities;
  DelocalizeCharge previously used raw attribute assignment. Final performance
  comparisons remain in S9b.

  **Verification.** 860 focused tests pass across transformations, aromatic
  insertion, resolver/constraint consumers, and kekulization fixtures. These
  include exact output, disjoint systems, unchanged input, and late valence/spin
  rejection with receiver restoration. Graph library checking and focused strict
  Clippy pass; nightly formatting and diff checks pass. The S6 workspace gate
  remains at S6d; remaining boundary/property-fixture callers belong to S6c2.

- **S6c2 — completed 2026-09-28 — Remaining ingest, parse, and export callers** (`umol-graph`, `umol-io`; breaking, green at S6d). [dep: S6c1]

  Migrate remaining edit/build/snapshot call sites to the owning editor and
  finish while retaining their current public signatures. Preserve each source
  copy already required for an independent output; do not add recovery copies.
  Resolver names and report removal change in S7, not here. Verify boundary
  outputs, error categories, and unchanged retained sources.

  **Implemented.** The remaining editor publication calls were in graph
  projection/publication properties and the C60 depiction fixture; all use
  finish. Projection fixtures read atom counts from the owning editor.
  Publication properties retain their whole-molecule equality assertions with
  explicit test copies. Production ingest and parse already use the migrated
  resolver; export retains its independent projection copy. No production
  signatures, error categories, or copying behavior change in this subitem.

  **Verification.** 610 selected tests pass: ingest, parse, export, projection
  and publication properties, and the C60 depiction fixture, with proptest and
  depiction features enabled. These retain exact boundary/error assertions and
  source-preservation checks. Focused strict Clippy, nightly formatting, and
  diff checks pass. No retired editor build/try_build/snapshot calls remain in
  umol-graph or umol-io. The workspace gate remains at S6d3.

- **S6d — completed 2026-09-28 — Python Molecule ownership migration** (group; breaking, closes
  S6 green at S6d3). [dep: S5d, S6a, S6b, S6c]

- **S6d1 — completed 2026-09-28 — Molecule Option storage and accessors** (`umol-py::molecule`; breaking, green at S6d3). [dep: S5d, S6a, S6b, S6c]

  The final Python Molecule fields and counter methods are specified under
  Python accessor invalidation. Complete their ownership interface as follows:

  ```rust
  pub(crate) fn from_rust(value: GraphIrMolecule) -> Self;
  pub(crate) fn to_rust(&self) -> PyResult<&GraphIrMolecule>;
  pub(crate) fn to_rust_mut(&mut self) -> PyResult<&mut GraphIrMolecule>;
  pub(crate) fn take(&mut self) -> PyResult<GraphIrMolecule>;
  ```

  from_rust initializes Some(value) and counter zero. take leaves None and
  raises ConsumedError if already empty; it needs no counter increment because
  check_access rejects consumed owners. No refill path exists. Explicit copies
  copy molecular storage deliberately and start a new wrapper at counter zero.
  Existing editor Option storage follows the same destructive failure rule;
  retain its current private storage rather than introduce a second wrapper.

  Test move-out, repeated take, explicit copies, and equality based on molecular
  values. Fields and accessor signatures are fixed by the tables above.

  **Implemented.** Molecule stores Option<GraphIrMolecule>; all four accessors
  have the specified signatures. Molecule.copy creates an independent owner
  with counter zero. Fallible Python equality compares molecular values.
  Clone and automatic by-value Python extraction are removed from the wrapper.
  Counter arithmetic follows the data-type guide's ordinary-arithmetic rule.

- **S6d2 — completed 2026-09-28 — Fallible owner access throughout bindings** (Molecule/view consumers in `umol-py`; breaking, green at S6d3). [dep: S6d1]

  Migrate every read/write of the wrapper through its fallible accessors,
  including getters, repr/equality, conversions, collections, entity attributes,
  and nested constraints. Reuse S5d2 check_access and counter propagation.
  Do not convert unrelated forms to Option or add automatic clones. Verify
  consumed roots and invalidated children across every entity collection.

  Solution.Determined stores its molecule as Py<Molecule>. Construction retains
  the supplied Python object and its getter returns that same object; neither
  copies molecular storage. Consuming it is visible through every alias,
  including the Solution. Repr and molecular equality propagate ConsumedError.
  This replaces the generated by-value constructor/getter's implicit copies.
  Reaction.lhs already retains a Python owner; its access propagates the same
  consumed-state error. Existing copying constructors and independently produced
  operation results keep their current contracts.

  **Implemented.** Root operations, all eight entity families, collections,
  nested constraints, and Molecule conversions propagate fallible owner access.
  Tests cover consumed roots, invalidated descendants and exhausted iterators,
  independent copies, and shared molecules held by Solution and Reaction.

- **S6d3 — completed 2026-09-28 — Consuming entry points and editor publication** (`umol-py::{molecule,transaction,edit}`; breaking, red→green). [dep: S6d2]

  Transfer Molecule inputs on consuming calls and reuse S5d's Edits
  transfer; raise `ConsumedError` for the owner and `InvalidatedViewError` for
  owner-backed accessors. Replace Python snapshot/build with finish, without a blanket
  consumed-state migration of unrelated forms. Test aliases, nested views,
  successful and failed consumption, and explicit copies. Apply fallible owner
  access throughout the binding consumers, not only mutation entry points;
  getters, repr/equality, conversions, collections, and nested setters must not
  bypass the consumed/invalidation checks.

  **Implemented.** Molecule.edit/apply/tracked_apply transfer storage without
  cloning. Editor apply/tracked_apply return another owning editor; finish
  consumes it and checks integrity. The four Python snapshot/build methods
  are removed; no Python probe is exposed. Transfers are sequential: an
  unavailable Edits input leaves an already-transferred receiver consumed.
  Python argument type errors occur before transfer. Existing error categories
  are preserved. Tests explicitly copy inputs only when retaining them is part
  of the tested workflow.

  **Verification.** 20,510 workspace tests and six doctests pass (nine skipped
  or ignored in total). The all-feature Python-binding suite passes 1,654 Rust
  tests; the rebuilt Python 3.13 extension passes 1,864 pytest cases (two skipped).
  Workspace all-target strict Clippy, all-feature graph-ir/graph/Python
  warnings-denied rustdoc, nightly formatting, and full diff review pass.
  The final test cleanup passes 74 focused Rust cases and all-feature Python
  all-target strict Clippy. S6c records the feature-gated property/conformance
  checks for its Rust changes. S7a is next; MSRV remains a final S9 gate.

### S7 — Resolution, projection, and boundaries

S7a–S7d form one public signature/result migration, returning green at S7d.

- **S7a — completed 2026-09-29** (`umol-graph::ops::resolve` and phase modules; breaking, green at S7d)
  Reuse S5c's shared planning and plan_project methods, rename phase `plan` to
  `plan_resolve`, and implement consuming `resolve`/`project` alongside
  recovering `resolve_into`/`project_into`. Preserve phase order, rejection
  priority, and standalone-versus-composite isotope policy. Test accepted
  outputs, each later-phase rejection, both ownership contracts, and application
  versus integrity errors through MoleculeApplyError in phase-specific errors.
  S0b already implemented Solution's separate underdetermination payload; do not
  repeat that type change. [dep: S0b, S5c, S6d]
  The standalone ownership signatures follow the same substitution for each
  resolver below, with its existing contradiction C and error E:

  ```rust
  fn resolve(&self, molecule: Molecule) -> Result<Solution<Molecule, C, ()>, E>;
  fn resolve_into(&self, molecule: &mut Molecule) -> Result<Solution<(), C>, E>;
  ```

  | Resolver | C | E |
  | --- | --- | --- |
  | IsotopeResolver | IsotopeContradiction | IsotopeError |
  | AromaticityResolver | AromaticityContradiction | AromaticityError |
  | StereoResolver | StereoContradiction | StereoError |
  | BondsResolver | BondsContradiction | BondsError |
  | MulticenterBondsResolver | MulticenterBondsContradiction | MulticenterBondsError |

  Isotope/aromaticity/stereo also expose project/project_into with the same
  ownership and result shapes, substituting their corresponding ProjectError
  types for E. Their existing contradiction types remain C. No standalone
  report API is added. Valence admission/projection retain their existing
  signatures and immutable/no-op behavior; there is no invented valence edit
  path. An underdetermined consuming operation returns (), dropping its input;
  borrowed underdetermination preserves the receiver as already specified.

  **Implementation and verification:** All five standalone resolvers and the
  composite resolver expose the ownership split above; isotope/aromaticity/stereo
  and composite projection expose both routes. Consuming execution uses one owning
  editor without a recovery copy or journal; borrowed execution uses a transaction.
  Every successful route finishes/commits, including empty plans. The first selected
  projection phase plans against the original input; later phases use checked probes.
  Existing callers and benchmarks use the renamed borrowed methods. S7b owns the
  report split: composite resolve_into still returns its existing ResolveReport;
  the new consuming resolve returns Molecule or payload-free underdetermination.
  S7c/S7d own boundary and Python lifecycle changes.

  Both ownership routes are covered by expected-result and rejection cases, including
  late-phase failures, stale plans, all projection flag combinations, and output
  integrity. Error conversions preserve MoleculeApplyError through each phase's error
  type; no integrity predicate or error type changed. Validation: 1,979 graph unit
  tests, 683 resolution conformance cases, and 14 resolver/projection/isotope
  properties pass. Graph all-feature/all-target strict Clippy and warnings-denied
  rustdoc pass; the Python caller compiles under Python 3.13. Nightly formatting and
  diff review pass. Timing comparisons remain in S9b; scratch is empty.

- **S7b** (`umol-graph::ops::resolve`; breaking, green at S7d) Make reports opt-in
  with `resolve_with_report` and `resolve_into_with_report`; avoid default
  report-only collection while retaining resolution candidates. Test equal
  outcomes and final molecules with and without reports, including
  underdetermination. [dep: S7a]

  Default paths inspect ResolveState.completions directly for entries with
  more than one candidate; they do not construct a report to inspect
  report.unresolved. Call to_report only in the explicit reporting methods.
  Preserve the current report payloads and rejection order. Sharing chemistry
  methods does not require a reporting-mode type or generic execution adapter.
- **S7c** (`umol-graph::ingest`, `parse`, `export`, `umol-io` boundaries; breaking,
  green at S7d) Pass owned ingest/parse candidates to report-free resolution; keep
  one intentional source copy for export projection. Remove report payloads from
  default underdetermination errors. Test boundary diagnostics and output
  integrity. [dep: S7a, S7b]
- **S7d** (`umol-py::resolve` and boundary adapters; breaking, red→green)
  Mirror the consuming/borrowed names, explicit reporting, and payload-free
  ingestion underdetermination. Test Python ownership, result shapes, and error
  parity. [dep: S6d, S7a, S7b, S7c]

### S8 — Transformations and remaining operations

S8a and S8b are one trait/implementation change, returning green at S8b. No
temporary clone-based default transform implementation is introduced between them.

- **S8a** (`umol-graph::ops::transform`; breaking, green at S8b) Change Transformer
  to consuming `transform`, recovering `transform_into`, and lazy
  `transform_iter` returning `impl Iterator`; migrate all trait callers. Test
  failure ownership, independent iterator outputs, and deferred execution.
  [dep: S5c, S6d]
- **S8b** (`umol-graph::ops::transform::{aromatizer,delocalize_charge,kekulizer}`;
  breaking, red→green) Reuse S6c's planning and borrowed execution for the new
  consuming route;
  use one publication gate and retain existing chemistry errors, including
  DelocalizeCharge's `Infallible` contract. Test each operation's success,
  rejection, and checked output. [dep: S8a]
- **S8c — moved to S6b.** combine_from caller migration is required during the
  owning editor change, not a later independent optimization.

### S9 — Contract and performance closeout

- **S9a** (`docs/development`, public rustdoc, examples; additive) Reconcile
  the living data-type, nomenclature, integrity, and Python API guides with the
  implemented lifecycle; remove stale snapshot/detached-journal descriptions.
  Document precise failure, ownership, and integrity boundaries without citing
  discussion records from source. Check links and examples. [dep: S5a, S6d, S7d, S8b]
- **S9b** (workspace gates and benchmarks; additive) Compare the final direct,
  apply, and transaction paths with S0, including dense resolver journal size
  and the source-preserving export path. Compare the existing reaction-application
  benchmark cases against the S0 revision using the same fixtures, toolchain,
  and configuration; record the timing and allocation changes. S5b compiled
  these benchmarks but did not measure them. Review the complete diff against 213,
  run formatting, workspace tests, strict Clippy and rustdoc, explicit feature
  suites, Python 3.13 build/tests, and the pinned Rust 1.87 gate once at final
  closeout. Record results and update the discussion status only after the full
  scope passes. [dep: S6b, S9a]

**Execution order:** completed S0a–S0b → completed S1a–S1c → S2 → S3 →
S4 → S5 → S6 → S7/S8 → S9. S2a is implemented; its API is removed in S2i1.
Within the revised S2:

- S2b, S2c, and S2d are complete; S2f is cancelled.
  S2c → S2d establishes coset-operation and count-use behavior. S2e is folded
  into S2h; it is not an executable prerequisite.
- S2g depends on S2c; its interfaces and failure behavior are approved.
- S2h removes aggregate attribute checks only after S2c/S2d/S2g
  consumer changes are complete, then tests the newly admissible inputs through
  public construction and those consumers before closing the subitem.
- S2i1–S2i5 are complete, including separate mutable molecule/editor view types
  with typed storage and unrestricted attributes. S2j is complete: matching local
  getters and editor structural operations compile and pass their focused gates
  without depending on S2k/S2l.
- S2k and S2l migrate Rust and Python callers; both precede removal in S2m.
- S3a depends on S2i; S3b requires S2j. S3c/S3d remove replacement Deltas and
  retain reaction integration; S3e closes the Python Edit migration.
- S3f and S3g supply graph-core bulk additions; S3g → S3h supplies typed-set
  extend, then S3f/S3h → S3i supplies Molecule/editor bulk additions. S3j changes
  correspondence mutation to mutable borrowing and migrates its callers.
  S3k1–S3k4's index-overflow cleanup, S4a, S4b, and S4d are complete.
  S4b uses the additions and the component
  removal/restoration interfaces.
- S4a, S4b, and S4d are complete. S4c is incorporated in S5a's guard and
  scoped run. S5a–S5d3 are complete; S5's build and test gate passes.
  S6a–S6d3 are complete; S7a is next.
- S5d1–S5d3 complete Python ownership, counters, and prepared transactions.
- S6b1/S6b2 separate caller migration from combine_from; S6c1/S6c2 separate
  chemistry and format callers. S6d1–S6d3 close the Python owning migration.

S0c/S1d remain reverted, and the former S0d/S0e work is represented by S2i/S2j.
Former S8c is included in S6b. The DSL payloads and Python consumption,
invalidation, and storage interfaces remain approved. No code or source API
change is authorized merely by writing a proposed interface in this plan.
No speculative optimization stage is required.
Immutable-view simplification, graph-core bond
endpoint rewiring, intermediate transaction tracking, an interactive Python
transaction, Undo compression, and the hydrogen operations in 166 remain outside
this plan; the 213 mutation surface enables the latter work.
