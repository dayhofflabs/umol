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
[data-type guide](../docs/development/data-types.md),
[nomenclature guide](../docs/development/nomenclature.md)

## Design status — 2026-09-26

This document owns the molecule/reaction mutation redesign. S0a–S0b, S1a–S1c,
S2a, the revised S2b, and S2c are implemented. The previous S2b mutable-view
attempt was reverted. S2d is next; S2f is cancelled and the remaining S2 work is
unimplemented. Graph-core mutation and restoration are complete in
[166](166-molecule-ops-2026-07-27.md); editor integration remains here. After
that integration, return to 166 for the operation changes and hydrogen folding.
Doc 228 is unchanged by this review and its withdrawn ownership migration is not
an implementation dependency.

| Area | Status | Concrete position |
| --- | --- | --- |
| Storage delegation, participant methods, Edit/Delta/Undo variants, local getters | Settled design; S1a–S1c complete | Use the existing typed entity sets and graph-core mutation/restoration; contracts below. |
| Editing and recovery | Settled design | Owning, destructive editor; separate borrowed, scoped transaction. Editor and Transaction probe check integrity and return an immutable Molecule borrow; no probe callback. |
| resolve/project/transform consumers | Settled design; integration work remains | resolve/project consume destructively; resolve_into/project_into mutate borrowed inputs with recovery. Consuming resolution uses Solution<Molecule, C, ()>; reporting is explicit. Ingest uses report-free resolution. Transformer signatures follow the same ownership naming. |
| Molecule attribute methods | Uniform unchecked attribute mutation settled; implementation remains | Mutable borrows expose every entity attribute and entity-level constraint in Molecule and MoleculeEditor. Rust and Python retain simple assignment, including aromatic/multicenter/stereo. Remove modify/try_modify callbacks. |
| Mutable-view structures and API | Settled design; S2i remains to implement | One const-generic mutable-view family. EDITOR controls structural mutation only; all attributes are freely mutable for both values. |
| Molecule-level constraint mutation | Editor-only design settled; S2a API scheduled for removal | Retain Molecule::constraints for reads and editor &mut Constraints for writes. S2i1 removes the public checked constraint view; S2k/S2l migrate callback callers and S2m removes try_modify_constraints. |
| Transaction correspondence | Settled design | tracked_commit returns the whole transaction's correspondence. Omit Transaction::tracked_apply unless a concrete need for intermediate tracking arises. |
| Python bindings | Prepared-batch transactions, consumption, and accessor invalidation settled; implementation remains | Molecule.transact and tracked_transact submit prepared Edits; Rust applies and commits within one borrowed transaction. No interactive Python Transaction or scoped TLS dependency. Molecule and Edits input-transfer changes remain; the editor already supports consumption. |
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
the consumer changes enabling unchecked attribute assignment. The const-generic
view design gates only editor structural mutation; attribute mutation is identical
for both values of EDITOR.
The lift_constraints defect and its undetermined-stereo policy are a separate
focused correction, recorded under
[other operations](#other-moleculereaction-operations).

Replacement verbs and payloads in the Edit/reaction DSLs are approved in S3.
Python consumption, counter-based accessor invalidation, and storage names are
approved below. S2b is complete: Rust's unit error is NoJoinError and Python
join raises NoJoinError. S2c's bounded coset-operation fixes are complete;
S2d is next and S2f is cancelled. S2g's frame-consumer decisions remain approved
and unimplemented.

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

The block expresses the accepted editor and transaction interface. Currently
Molecule::edit and Molecule::apply borrow &self and create an independent result;
editor transact returns a detached journal; snapshot/try_build/build publish the
working state. Those lifecycle methods are replaced by the surface above.
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
borrows Molecule directly and calls the same private mutation methods as the
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

Required caller migrations remain: charge delocalization and MoleculeDsl
conversion use an editor; Python preserves field assignment and direct
entity-constraint mutation; molecule-level constraint writes use the editor.
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
| StereoResolver::project; StereoPerception::derive; export::lower_stereo_atom / lower_stereo_bond; Reaction::application_deltas | Projection reports failed transport; perception records failed entity candidates; export rejects unsupported literal indices. Reaction application explicitly checks stereo domains. | Retain existing consumer behavior; do not introduce a shared normalization range check or additional validation pass. |

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
    fn constraints_mut(&mut self) -> &mut Constraints;
}

impl MoleculeEditor {
    pub fn constraints(&self) -> &Constraints;
    pub fn constraints_mut(&mut self) -> &mut Constraints;
}
```

The private Molecule accessor supplies the editor's mutable borrow by delegation;
it is not part of Molecule's public mutation surface.
The editor borrow supports the existing Constraints operations and whole-collection
assignment. It does not check writes individually or clone the collection.
Remove MoleculeEditor::push_constraint; use constraints_mut().push instead.
The private Molecule::push_constraint remains for Edit execution.
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

**Approved; implementation remains in S5d/S6d.** The Option
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
usable. Existing owned return values are not converted to live views by this work.

Use one private u64 counter in the Python Molecule wrapper and a captured
counter in each owner-backed accessor. Advancing it uses checked addition;
overflow raises Python OverflowError before execution, never wraps to revive an
old accessor. Check owner availability and counter under the same short
PyO3 borrow used for access. Borrow the receiver exclusively and advance its
counter after preparing Python inputs, then retain that borrow throughout
Rust execution. No Python callback runs between invalidation and completion.
An unwind after execution starts also leaves old accessors invalid.

Concrete wrapper layout (final shape after S6d; S5d adds counter while the
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
Standalone storage variants are unchanged. Current constraint/key/item iterators
own copied entries or keys; they need no counter because they do not consult
the molecule. Their separate copying policy is unchanged.

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
returns (). advance_counter checks owner availability, then assigns
counter.checked_add(1), mapping overflow to OverflowError. S5d's checks have
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

The participant, Edit, Delta, and Undo contracts below remain in place. S2a's
molecule-level ConstraintsViewMut is implemented but scheduled for removal in
S2m. The previous mutable-view
implementation remains reverted; S2i must implement the uniform attribute access below after the required
first-use checks and aggregate-integrity changes.

### Mutable-view structures and access

Molecule and MoleculeEditor expose mutable entity views for all eight kinds:
atom, localized bond, dative bond, aromatic system, multicenter bond,
noncovalent bond, stereo atom, and stereo bond. Every view provides a mutable
borrow of the complete entity form; all attribute fields and entity-level
constraints remain directly assignable. The same assignment semantics apply in
Python. There are no exceptions for electron counts or stereo configuration.

Use const EDITOR: bool = false solely to control structural mutation: public
Molecule access returns views with EDITOR = false; MoleculeEditor returns
EDITOR = true. Internal batch and undo execution can obtain the same <true>
views through private Molecule accessors. Both
expose the same unrestricted mutable attribute borrow. Do not add a runtime
permission flag, checked field setters, or
modify/try_modify callbacks. Do not introduce individual charge_mut,
element_mut, electrons_mut, or coset_mut accessor families to replace ordinary
field assignment. Attribute reads and writes must have the same shape across
entity kinds and across molecule/editor access.

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

The same view structure serves Molecule and MoleculeEditor. Common read and
attribute methods are defined for either value of EDITOR; atom/donor/ligand/site
replacement is defined only for EDITOR = true. Structural changes can break
reference, uniqueness, or incidence integrity, so the editor publication boundary
must check them. Attribute assignment does not change those structural links.
EDITOR enables the structural capability used by the editor and internal batch
execution; both kinds of view support attribute editing. The const parameter
adds no stored flag or runtime branch. Do not introduce a
second *EditorViewMut family or a conversion enabling structural mutation on a
molecule view.

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

### Local getters

The following local getters are approved for both the unified mutable views
and the existing immutable editor views.
Immutable molecule views retain their current API. Editor getter methods already
exist; the local methods on their return values are mostly additions.

Currently all editor views expose id and attributes fields. Bond and noncovalent
views also expose atoms; stereo views expose site and ligands. Dative, aromatic,
and multicenter editor views already have atom_ids(). Replace editor-view public
fields with accessors, including the common methods above (attributes_mut only on
mutable views). Use the same read method signatures on both editor view families.

Immutable editor views keep their existing stored data; their fields become
private. Construct them through the following crate-private new methods on
impl<'a>, each returning Self. These constructors take storage data, not views.
The dative/aromatic/multicenter signatures already exist; the others replace
struct literals. No new factory methods are added to the typed sets.

| Immutable editor view | new arguments |
| --- | --- |
| AtomEditorView | id: AtomId, attributes: &'a AtomForm |
| BondEditorView | id: BondId, atoms: [AtomId; 2], attributes: &'a BondForm |
| DativeBondEditorView | id: DativeBondId, donors: &'a [NodeId], acceptor: AtomId, attributes: &'a DativeBondForm |
| AromaticSystemEditorView | id: AromaticSystemId, atoms: &'a [NodeId], attributes: &'a AromaticSystemForm |
| MulticenterBondEditorView | id: MulticenterBondId, atoms: &'a [NodeId], attributes: &'a MulticenterBondForm |
| NoncovalentBondEditorView | id: NoncovalentBondId, atoms: [AtomId; 2], attributes: &'a NoncovalentBondForm |
| StereoAtomEditorView | id: StereoAtomId, site: AtomId, ligands: &'a [StereoLigand], attributes: &'a StereoAtomForm |
| StereoBondEditorView | id: StereoBondId, site: BondId, ligands: &'a [StereoLigand], attributes: &'a StereoBondForm |

Every signature below takes &self. Ordinary getters return borrowed stored forms;
ids and StereoLigand values are copied. Iterator results are lazy and borrow the
view. ligand_frame borrows the stored slice without allocating; callers that
need ownership explicitly copy it.

| Entity view | Getter signatures |
| --- | --- |
| Atom | element() -> &ElementForm; isotope_mass() -> &IsotopeMassForm; charge() -> &NumForm; implicit_hydrogens() -> &NumForm; lone_pairs() -> &NumForm; unpaired_electrons() -> &UnpairedElectronsForm |
| Localized bond | atom_ids() -> [AtomId; 2]; order() -> &NumForm; charge() -> &NumForm; unpaired_electrons() -> &UnpairedElectronsForm |
| Dative bond | donor_ids() -> impl ExactSizeIterator<Item = AtomId> + '_; acceptor_id() -> AtomId; atom_ids() -> impl ExactSizeIterator<Item = AtomId> + '_; donor_count() -> usize; atom_count() -> usize; order() -> &NumForm |
| Aromatic system and multicenter bond | atom_ids() -> impl ExactSizeIterator<Item = AtomId> + '_; atom_count() -> usize; electrons() -> &ElectronCountsForm; electron_count() -> NumForm; charge() -> &NumForm; unpaired_electrons() -> &UnpairedElectronsForm |
| Noncovalent bond | atom_ids() -> [AtomId; 2]; kind() -> &NoncovalentBondKindForm |
| Stereo atom and stereo bond | site_id() -> AtomId / BondId respectively; configuration() -> &StereoConfigurationForm; kind() -> Option<StereoKind>; coset() -> Option<&StereoCoset>; plus the ligand methods below |

Dative atom_ids yields the donors in stored order followed by the acceptor.
Aromatic/multicenter electron_count follows the existing getter: sum the stored
literal contributions, otherwise return NumForm::Undetermined. It does not check
whether the contribution vector fits the atom sequence. Stereo kind/coset follow
the configuration form's optional accessors, so an undetermined draft is readable.

Both stereo editor view families provide exactly these ligand methods:

```rust
pub fn ligands(&self) -> impl ExactSizeIterator<Item = StereoLigand> + '_;
pub fn ligand_count(&self) -> usize;
pub fn ligand(&self, position: StereoLigandPosition) -> StereoLigand;
pub fn ligand_position(&self, atom: AtomId) -> Option<StereoLigandPosition>;
pub fn ligand_frame(&self) -> &[StereoLigand];
pub fn atom_ligands(&self) -> impl Iterator<Item = StereoLigand> + '_;
pub fn implicit_hydrogen_ligands(&self) -> impl Iterator<Item = StereoLigand> + '_;
pub fn lone_pair_ligands(&self) -> impl Iterator<Item = StereoLigand> + '_;
pub fn atom_ligand_ids(&self) -> impl Iterator<Item = AtomId> + '_;
pub fn implicit_hydrogen_atom_ids(&self) -> impl Iterator<Item = AtomId> + '_;
pub fn lone_pair_atom_ids(&self) -> impl Iterator<Item = AtomId> + '_;
pub fn atom_ligand_count(&self) -> usize;
pub fn implicit_hydrogen_count(&self) -> usize;
pub fn lone_pair_count(&self) -> usize;
```

ligand panics out of range. ligand_position returns the first matching actual-atom
ligand, not a virtual ligand bearing that atom id. Filters preserve stored order;
virtual-ligand id accessors return their bearing atom ids. constraints returns the
stored ConstraintsForm, not a molecule-backed constraint view.

Neighbors, valence, induced bonds, resolved entity views, stereo-bond endpoint
lookup, and frame-transport queries require more than this local borrow. They
remain on immutable Molecule views and are not added to either editor view family.
No additional graph or Molecule borrow is stored merely for getter parity.

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
replacement uses the existing <true> views, which delegate to the typed sets.
This boundary adds no second mutation vocabulary or public mutable access.

Graph bulk addition is approved as three mutable operations: add_nodes for
isolated nodes, add_edges for edges between existing nodes, and add for nodes
and edges together. The bare verb covers both topology components, consistently
with the existing graph nomenclature. Combined-add edge endpoints use the
resulting graph's node ids, including the appended nodes; there is no separate
handle namespace. Existing ids stay unchanged and each added kind occupies a
contiguous block. Return owned, allocation-free, exact-size iterators over those
blocks, following the existing typed node_ids/edge_ids iterator interfaces:

```rust
impl Graph {
    pub fn add_nodes(
        &mut self, count: usize,
    ) -> impl ExactSizeIterator<Item = NodeId> + use<>;
    pub fn add_edges(
        &mut self, edges: &[[NodeId; 2]],
    ) -> impl ExactSizeIterator<Item = EdgeId> + use<>;
    pub fn add(
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
operations; combine_from can use combined add. Node-only addition can extend
adjacency offsets without rebuilding existing edges. Combined addition builds
the final adjacency once, avoiding an intermediate copy-on-write detachment for
node addition followed by an edge rebuild. Single additions delegate to their
bulk counterparts; no timings are claimed by this design decision.

Relation sets retain individual add and gain public extend. Each extend takes
&mut self plus the batch below and returns an owned, allocation-free
`impl ExactSizeIterator<Item = RelationId>`. Returned iterators borrow neither
the receiver nor the supplied participant slices. Precise captures must exclude
those lifetimes while accounting for the enclosing type/const parameters.

| Relation set | extend batch argument |
| --- | --- |
| FixedRelationSet<P, N, D> | Vec<([P; N], D)> |
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
| add_atom, add_bond | Move implementations to private Molecule methods; editor methods delegate and batch execution calls the same implementations. Retain Graph addition and parallel attribute-array updates; direct additions record no correspondence. |
| add_dative_bond, add_aromatic_system, add_multicenter_bond, add_noncovalent_bond, add_stereo_atom, add_stereo_bond | Move implementations to private Molecule methods; editor methods delegate and batch execution calls the same implementations. Replace mutable-entry-vector insertion with relation add through the typed set; direct additions record no correspondence. |
| atom_mut, bond_mut | Already mutate copy-on-write attribute arrays. Retain graph-IR ownership of these arrays. |
| Mutable attribute views for all six overlay kinds | Replace entry-vector materialization with typed-set access to relation payload mutation. Keep attributes and participants accessible in the entity view. |
| constraints_mut, inline constraints through attribute views | Remove editor push_constraint; use constraints_mut().push. Keep private Molecule::push_constraint for Edit execution. Retain editor-only top-level mutable access and shared entity-attribute access; graph-core does not interpret molecular constraints. |
| Each overlay remove_* / tracked_remove_* pair | Private Molecule methods remove entries from that owning set only; tracked forms return its typed compaction. Complete editor/batch removal separately assembles MoleculeCompaction and compacts constraints. Use the set's removal implementation; remove public tracked direct mutation. |
| Topology remove / tracked_remove | Rename the editor operation to remove_topology. Private Molecule remove_topology / tracked_remove_topology mutate Graph and atom/bond attributes only; tracked removal returns GraphCompaction. Editor/batch execution separately compacts each overlay set, assembles MoleculeCompaction, and compacts constraints. Graph nomenclature is unchanged. |
| Internal undo-addition removal and restore_* methods | Undo additions with private Molecule untracked removal only, without overlay or constraint compaction. Undo removals with restore_topology for Graph and atom/bond attributes, separate overlay topology-id/row restoration, then constraint restoration. S4a lists the restoration interfaces. |
| Overlay and constraint compaction | Private Molecule compact_* / tracked_compact_* methods each delegate to one component. Overlay delegates install the returned replacement set; tracked forms also return its row mapping. Constraint delegates mutate the collection in place. S4b lists their interfaces. |
| apply, transact, and tracked counterparts | Replace the current lifecycle with the editor/transaction API above. Both execution paths share graph-IR mutation operations; handle resolution and forward preconditions remain batch concerns, and undo capture/replay remains transactional. |
| snapshot, try_build, build, and tracked counterparts | Remove relation-row rebuilding at publication. Replace editor publication with finish, Molecule::apply, or scoped commit, and replace snapshot with checked probe access. Fresh MoleculeBuilder retains asserted build. |
| Overlay reads and views; six internal *_equiv methods | Remove storage-wrapper dispatch. Use typed-set accessors and graph-core participant comparison where applicable; retain graph-IR frame transport and payload comparison. |

The private Molecule addition interfaces are:

```rust
fn add_atom(&mut self, attributes: AtomForm) -> AtomId;
fn add_bond(&mut self, first: AtomId, second: AtomId, attributes: BondForm) -> BondId;
fn add_dative_bond(&mut self, donors: &[AtomId], acceptor: AtomId, attributes: DativeBondForm) -> DativeBondId;
fn add_aromatic_system(&mut self, atoms: &[AtomId], attributes: AromaticSystemForm) -> AromaticSystemId;
fn add_multicenter_bond(&mut self, atoms: &[AtomId], attributes: MulticenterBondForm) -> MulticenterBondId;
fn add_noncovalent_bond(&mut self, atoms: [AtomId; 2], attributes: NoncovalentBondForm) -> NoncovalentBondId;
fn add_stereo_atom(&mut self, site: AtomId, ligands: &[StereoLigand], attributes: StereoAtomForm) -> StereoAtomId;
fn add_stereo_bond(&mut self, site: BondId, ligands: &[StereoLigand], attributes: StereoBondForm) -> StereoBondId;
fn push_constraint(&mut self, constraint: Constraint);
```

Bulk additions have identical interfaces on private Molecule methods and public
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

Each row below specifies a private Molecule removal pair. Both methods take
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
tracked_remove_topology, call the six Molecule tracked_compact_* overlay
delegates with its GraphCompaction, assemble MoleculeCompaction::new from the
graph and six row mappings, then call compact_constraints. Batch execution uses
that same sequence and retains the mappings it needs. Explicit overlay removal instead
assembles MoleculeCompaction from the one changed row mapping and identity
mappings for the unchanged kinds. No other overlay set changes in that path.
The typed-set compact/tracked_compact interfaces are recorded in S1a–S1c.

These compositions may remain explicit in their consumers. Do not introduce
combined helpers, output parameters, or recording modes solely to share the
sequencing. Untracked-to-tracked delegation is the current implementation, not
a settled efficiency requirement. Remove the seven remove_added_* adapters;
each Undo match arm extracts saved ids and calls the private Molecule untracked
removal for the added kind directly. Under matching
history, reverse replay has undone later dependencies and returned added entries
to their appended positions; earlier ids do not change. No separate overlay or
constraint compaction, combined MoleculeCompaction, or correspondence update is
needed for addition undo. Retain local panic guards for manipulated history,
without validating that it is matching history or requiring a trailing block.

Restoration follows the same component boundaries:

```rust
// Private Molecule method.
fn restore_topology(
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

Single-entry execution belongs to private Molecule methods: apply_edit,
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
| AddAtoms, AddBonds, six overlay additions | Shared private Molecule add_* methods, with bulk topology addition for the existing AddAtoms/AddBonds variants. |
| RemoveTopology, six overlay removals | Private Molecule removal primitives followed by the separate overlay/constraint compaction steps described above; batch execution retains the combined mapping for handle updates. |
| Modify*Field for all eight entity kinds | Assignment through the existing/shared mutable attribute access. |
| Nine Replace* variants for atoms, donors, acceptors, sites, and ligands | Structural methods on the existing <true> mutable views obtained through private Molecule access; those views delegate to their owning sets. |

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
transaction borrowing the existing molecule. Both use the same private edit
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
  Set::compact(&self, graph: &GraphCompaction) -> Self
  Set::tracked_compact(&self, graph: &GraphCompaction) -> (Self, Compaction<Id>)
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
  Set::compact(&self, graph: &GraphCompaction) -> Self
  Set::tracked_compact(&self, graph: &GraphCompaction) -> (Self, Compaction<Id>)
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

- **S2d — Role-only incidence and count-aware consumers** (`ir::incidence`,
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

  Incidence matching uses role labels and the existing verify_overlays check of
  transported attributes. Remove the edge-payload count filter; do not introduce
  a public matching error merely for a rejected overlay candidate. At that
  candidate check, literal count/frame disagreement on either matched entity
  rejects the candidate through the existing Option path. Both matching
  algorithms must return identical valid-input match sets. Public matching
  signatures and SubstructureMatchError remain unchanged.

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

- **S2e — moved to S2h.** Existing electron-use boundaries and getter behavior
  require verification after aggregate construction admits the relevant forms.
  They do not require a separate implementation subitem.

- **S2f — cancelled.** No depiction/CoordGen error-handling, error-enum, or
  private-signature changes are included in this work.

- **S2g — Frame-dependent attribute consumers** (`ir::canonicalize`,
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
  newly admissible count and coset payloads; faithful storage and serialization;
  continued rejection of structural defects. Exercise the first-use rejections
  through public APIs after construction, and the unchanged getter behavior.
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
  - Standalone electron transport still returns None on a degree mismatch;
    the existing chemistry boundaries retain ElectronCountMismatch, no applicable
    plan, and AromaticSystemFailure. Getter cases use [1, 2] and [1, 2, 0, 5]
    for a three-atom entity and assert the accepted per-atom and total behavior.

  S2h is not complete when checks are merely deleted: these consumer regressions
  and the retained structural-rejection cases must pass in the same subitem.
  Update data-types.md and integrity.md with the actual retained checks now;
  do not leave normative guides contradicting the implementation until S9.

- **S2i — Uniform mutable attribute access and typed editor storage**
  (`ir::{view,molecule,molecule::editor}`, ir exports; group; breaking, green at S2i4).
  [dep: S1a, S1b, S1c, S2h]

  Execute S2i1–S2i4 in order; all S2i contracts apply to the group.

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

  Add private Molecule::constraints_mut(&mut self) -> &mut Constraints in the
  owning molecule module. S2i3 uses it for the editor delegate; no public
  top-level constraint borrow is added. Remove the checked-view description
  from the nomenclature guide. Keep try_modify_constraints until S2k/S2l migrate
  its consumers, then remove it in S2m. The editor delegates to the private
  accessor rather than accessing Molecule fields directly.

  Verify editor top-level mutation and checked publication, retained read-only
  Molecule access, and absence of the checked view/type export.

- **S2i2 — Shared mutable-view types and Molecule access** (`ir::{view,molecule}`, exports; breaking, green at S2i4). [dep: S1a, S1b, S1c, S2i1]

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
  accessors return <false>; private access supplies the existing <true> views to
  editor, batch, and undo execution. Editor accessors delegate to that private
  access rather than borrowing Molecule fields to construct views themselves.
  Spell that private access as the following Molecule methods, each taking
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
  live in the molecule module so its editor/transact children can call them
  without pub(crate) visibility. Top-level mutable constraint access uses the
  private Molecule::constraints_mut() -> &mut Constraints introduced in S2i1.
  Internal attribute-only writes can use the ordinary <false> views. No new view
  family, public const-selection parameter, or false-to-true conversion is added.
  Export the shared family through ir. id copies the id; attributes,
  attributes_mut, constraints, and structural getters borrow for the method
  call, not 'a. No public conversion changes false to true. EDITOR affects
  available methods only, with no stored flag or runtime branch.

  Add focused cases for uniform attribute access and const-specialized capability.

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
  Structural-view mutation and its exclusion from <false> are tested in S2j.
  Review editor mutable accessors to confirm they obtain views from Molecule;
  storage borrows and copy-on-write mutation remain inside Molecule and views.

- **S2j — Local getters and editor structural mutation**
  (`ir::{id,view}`, typed sets, ligand_frame callers and bindings; additive
  structural methods plus breaking getter return change, red→green). [dep: S2i]

  **Semantics.** Read-only getters borrow existing data; filters are lazy and
  preserve stored order. ligand_frame borrows the stored ordered slice, with no
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
  including ligand(StereoLigandPosition) -> StereoLigand and optional stereo
  kind/coset access for an undetermined editor configuration. Count getters stay
  infallible with the accepted malformed-input behavior from S2h.

  Change the existing immutable stereo views and give both const specializations
  of the mutable stereo views the same borrowed-frame contract:

  ```diff
  -pub fn ligand_frame(&self) -> Vec<StereoLigand>;
  +pub fn ligand_frame(&self) -> &[StereoLigand];
  ```

  Immutable views borrow their stored ligand slice; mutable views borrow through
  the owning set. The returned slice is tied to the accessor borrow, preventing
  structural mutation while that slice is in use. No mutable frame slice is
  exposed. Migrate all callers in this subitem: canonicalization, correspondence,
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
  in S2i/S2j. This subitem changes the remaining callback-based mutation flows
  and finishes green.

  **Semantics.** MoleculeDsl conversion uses one editor for the full pass and
  publishes once. Preserve FromIr's intentional source copy, defaults, metadata,
  ordering, and faithful conversion; IntoIr does not gain an extra recovery copy.
  Charge delocalization uses one editor for its complete mutation pass, preserving
  its planning checks and Infallible result. Direct attribute writes replace
  whole-form callback roundtrips.

  **Interfaces and nomenclature.** Keep existing FromIr/IntoIr signatures and
  DelocalizeCharge's current transformation signature. Its rename to
  ChargeDelocalizer is tracked separately in
  [166](166-molecule-ops-2026-07-27.md#charge-delocalization-transformer-name).
  Use the approved mutable
  view names/accessors from S2i and current edit/build lifecycle here; S6 owns
  consuming edit/finish. No mem::take of a borrowed caller's molecule, temporary
  empty receiver, new modify helper, or transitional compatibility callback.
  Migrate ordinary calls and callback method names supplied through macros.
  Move remaining try_modify_constraints callers to editor mutation and
  publication. S2i1 already migrated the public checked-view callers; retain
  reference/frame integrity cases at constructor/editor publication
  boundaries, including deliberately invalid inputs and order preservation.

  **Verification.** DSL roundtrips, preserved source semantics, charge
  delocalization outcomes, and tests/fixtures that used callbacks. Preserve
  property laws and invalid-input cases; do not make fixtures valid merely to
  avoid the newly explicit first-use failure behavior.

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
  Prepared transactions, consumption/counters, and their scheduled S5 work are
  not redesigned here; no interactive Python transaction or new editor-view
  binding family is introduced.

  **Verification.** Python 3.13 assignment tests for all eight kinds, whole-form
  replacement, nested constraints, short/long electron vectors, out-of-range
  cosets, and continued access to backing storage. Check that molecule-level
  constraints cannot be written through the property, collection, or nested
  accessors, while standalone/entity-level constraints remain mutable. Preserve
  invalid top-level reference/frame rejection at editor publication.
  Build the extension before tests and exercise public Python access rather
  than test-only Rust mutation paths.

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
  attribute access remains; the const generic gates structural methods only.
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

### S3 — Batch mutation and reaction vocabulary

S3a–S3e form one breaking enum migration; the green boundary is S3e, including
Rust execution, DSL conversion, and Python exhaustive matches. Do not make an
intermediate subitem compile by adding wildcard rejection, dropped variants,
or unimplemented execution branches.

Approved DSL verbs match the Rust method names, using hyphens in EDN. Both the
Edit and reaction DSLs use the existing entity key to select the entity kind:

| Entity keys | Rust methods | DSL verbs |
| --- | --- | --- |
| :aromatic-system, :multicenter-bond, :noncovalent-bond | replace_atoms | :replace-atoms |
| :dative-bond | replace_donors, replace_acceptor | :replace-donors, :replace-acceptor |
| :stereo-atom, :stereo-bond | replace_site, replace_ligands | :replace-site, :replace-ligands |

Approved payloads follow the respective existing :modify conventions:

| DSL | Payload | Old value |
| --- | --- | --- |
| Edits | [handle {:expect old :update new}] | Supplied explicitly; realized and compared by batch execution. |
| Reaction | [target new] | Derived from the lhs component, as for existing :modify. |

Atom and donor lists are EDN vectors, preserving order. Noncovalent atom vectors
contain exactly two entries. A site or acceptor is one atom/bond handle or target
in the corresponding DSL's existing encoding. Ligand lists are vectors using
the existing ligand encoding. Both :expect and :update carry these component
values directly; they do not use attribute-form strings.

```edn
;; Edits
{:aromatic-system
 {:replace-atoms [0 {:expect [0 1 2] :update [0 1 3]}]}}

;; Reaction
{:aromatic-system
 {:replace-atoms [0 [0 1 3]]}}
```

Retain each DSL's existing handle/target resolution conventions. Single-position
replacement/insertion/removal has no separate Edit or Delta variant in this
design and gains no DSL operation here.

- **S3a — planned; interfaces approved** (`ir::edit`, `dsl::edit`; group; breaking, green at S3e)
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

- **S3a3 — Edit DSL replacement syntax** (`dsl::edit`; breaking, green at S3e). [dep: S3a1, S3a2]

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

- **S3b** (`ir::molecule::transact`; breaking, green at S3e) Realize those edits
  through the structural methods on the existing <true> mutable views, obtained
  from Molecule's private access introduced in S2i. Undo uses those same methods
  with the saved components. Do not reach into typed sets from edit execution
  or construct views there from Molecule fields. Keep unrelated factors,
  attributes, constraints, and ids unchanged. Test each
  forward/undo pair, old-state errors before mutation, and frame alignment.
  Migrate dative addition/removal and undo capture/replay to separate donors
  and acceptor; remove combined-vector joining/splitting and pass donor slices
  to the editor. Preserve existing failure and rollback semantics.
  [dep: S2j, S3a]

  **Structural Edit execution and undo.** Each row below uses the existing
  <true> mutable view for the resolved entity id, obtained through private
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

- **S3c — planned; interfaces approved** (`ir::delta`; breaking, green at S3e)
  [dep: S0a]

  **Semantics.** Add whole-component replacement to the six existing entity
  Delta enums. Unlike Edit, these values use resolved typed ids, not Id/New
  handles. A replacement changes only its named component; attributes,
  constraints, and other components require their own deltas. List order is
  significant. No single-position variants or combined-component variants.

  **Interfaces and nomenclature.** Every new variant has exactly id, old, and
  new fields:

  | New Delta variant | id type | old/new type |
  | --- | --- | --- |
  | AromaticSystemDelta::ReplaceAtoms | AromaticSystemId | Vec<AtomId> |
  | MulticenterBondDelta::ReplaceAtoms | MulticenterBondId | Vec<AtomId> |
  | NoncovalentBondDelta::ReplaceAtoms | NoncovalentBondId | [AtomId; 2] |
  | DativeBondDelta::ReplaceDonors | DativeBondId | Vec<AtomId> |
  | DativeBondDelta::ReplaceAcceptor | DativeBondId | AtomId |
  | StereoAtomDelta::ReplaceSite | StereoAtomId | AtomId |
  | StereoAtomDelta::ReplaceLigands | StereoAtomId | Vec<StereoLigand> |
  | StereoBondDelta::ReplaceSite | StereoBondId | BondId |
  | StereoBondDelta::ReplaceLigands | StereoBondId | Vec<StereoLigand> |

  For example, add to AromaticSystemDelta:

  ```rust
  ReplaceAtoms {
      id: AromaticSystemId,
      old: Vec<AtomId>,
      new: Vec<AtomId>,
  }
  ```

  Extend existing id access, inversion, frame transport, normalization, and
  composition matches; do not add a parallel change algebra or new public
  method family. inverse(self) -> Self swaps old/new without changing id.
  Component continuity compares exact ids and stored sequence order, including
  both atom id and kind for each StereoLigand. Do not sort lists or implicitly
  transport attributes to make successive changes agree.

  Apply the existing folding rules to each named component:

  | Sequence | Normalized result |
  | --- | --- |
  | A → B, then B → C | A → C |
  | A → B, then C → D, with B != C | Existing contradiction result |
  | A → A | Identity eliminated |
  | Add followed by replacement | New component absorbed into Add |
  | Replacement followed by Remove | Remove contains the original component |
  | Add, changes, then Remove | Created entity cancels under the existing rules |

  Context-free Deltas normalization retains replacement variants when they do
  not fold away; it cannot invent the remaining entity attributes required for
  complete Remove/Add entries. S3d owns that lowering using the reaction lhs.
  Retain the existing FrameTransport failure boundary: one frame action cannot
  act on differently sized old/new lists. Incompatible replacements use S3d's
  before/after lowering before host alignment, rather than forcing a shared
  permutation or introducing separate left/right frames into ReactionSpan.

  The reaction DSL uses the replacement verbs listed at the start of S3, with
  [target new] payloads and old values derived from lhs, as for :modify. Its
  parsing/rendering and reaction integration are implemented in S3d; Python
  variant bindings follow in S3e. Undo variants belong to S3a, not this item.

  **Verification.** All nine variants: construction, id preservation, inverse
  involution, exact list/kind ordering, continuity and discontinuity, identity
  removal, Add/Remove folding, and created-entity cancellation. Verify supported
  frame transport and rejection of incompatible action degrees without silently
  changing either list. Reaction application and span-conversion tests belong
  to S3d. The complete enum migration reaches its green boundary at S3e.

- **S3d** (`ir::reaction`, `ir::reaction::integrity`, `ir::reaction_span`,
  `dsl::reaction`; breaking, green at S3e) Lower
  replacement deltas through before/after materialization and existing
  correspondence induction and superimposition. Test compatible preserved
  entities, incompatible removal/addition, application, and roundtrips without
  asserting that Delta spelling survives span conversion. [dep: S3b, S3c]

  Pass dative donors and acceptor separately when lowering existing Add/Remove
  deltas to Edits; stop joining them into one atom vector.

  S2h has already removed count/frame and per-entity stereo-attribute agreement
  from aggregate integrity. S3d does not reintroduce those checks for replacement
  deltas, either at reaction construction or final product publication.

  Extend reference visitation, remapping, and application-domain collection to
  every new old/new component. Extend the retained reaction reference checks to
  their entity ids, atom ids, bond sites, and ligand atoms. Preserve the existing
  removal/source incidence contract. Published molecules and reaction-span sides
  still satisfy the structural contract retained by S2h: valid references,
  relation uniqueness, ligand incidence, duplicate-ligand rejection, and the
  storage frame bound. Keep the retained reaction-span agreement of determined
  kinds for a preserved stereo entity; incompatible changes use Remove/Add.

  Before/after materialization remains necessary to obtain complete entities for
  span conversion and incompatible-change lowering. A single frame action cannot
  transport differently sized old/new lists. Use the agreed Remove/Add lowering
  before host alignment where required; first-use frame transport retains its
  own Option/Result failure behavior. Do not impose attribute agreement on
  intermediate edits or on construction merely because transport may need it.

  Test size-changing atoms/electron counts and ligands/configuration together,
  including nonidentity host alignment. Distinguish retained structural errors
  from count/coset disagreement: the latter alone must not fail construction or
  publication, but an operation requiring normalization or frame transport still
  reports its existing failure when the supplied attributes cannot be interpreted.

  Extend ReactionDsl's tree and streaming readers, rendering, and IR conversion
  with the approved replacement syntax, deriving old component values from lhs
  as for :modify. Test both parsing paths, aliases, ordering,
  roundtrips, and rejection of unavailable entities. Do not encode an
  identity-preserving replacement as remove/add merely to avoid designing its
  DSL representation; span conversion's separate contract remains unchanged.

- **S3e** (`umol-py::edit`, delta/reaction bindings; breaking, red→green)
  Extend the existing Python variant families for the new changes, with
  Rust-equivalent construction and failure behavior; migrate exhaustive matches
  and add parity cases. Do not add unrelated Rust API coverage merely for parity.
  Migrate existing dative Edit variants and Edits addition/removal methods to
  the separate donor/acceptor inputs from S3a, including plural entry tuples.
  This closes the enum migration. [dep: S3a, S3c, S3d]
- **S3f** (`umol-graph-core::graph`; additive, green) Implement Graph::add_nodes,
  add_edges, and add with the exact interfaces under Storage delegation. Return
  owned, allocation-free exact-size id iterators; mutation is eager and the
  iterators borrow neither receiver nor inputs. Preserve existing ids and append
  each kind contiguously. Combined-add endpoints use the resulting node space;
  invalid endpoints retain the graph addition panic contract. Single additions
  delegate to bulk operations. Node-only addition extends adjacency offsets;
  combined addition builds final adjacency once. Verify empty additions,
  appended ids and endpoints, loops/parallel edges, independence of shared
  clones, and continued graph mutation while returned iterators are held.
  Establish correctness cases and a focused native comparison of repeated versus
  bulk addition before replacing the implementation; use independently built
  expected graphs rather than single additions that delegate to the subject.
  [dep: none]
- **S3g** (`umol-graph-core::relation::{fixed,var,fixed_fixed,fixed_var,var_var}`;
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
- **S3h** (`ir::{aromatic,multicenter,dative,noncovalent,stereo}`;
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
- **S3i** (`ir::molecule`, `ir::molecule::editor`; additive, green) Add private
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

### S4 — Recovery machinery before the public lifecycle switch

- **S4a** (`ir::constraint::molecule`, `ir::molecule`, `ir::molecule::editor`,
  `ir::molecule::transact`; group; public rename and restoration rewire, green at S4a2). [dep: S2i, S3b]

- **S4a1 — Molecule topology and overlay restoration** (`ir::molecule`; additive, green). [dep: S2i, S3b]

  Route topology and relation restoration through the existing graph-core
  operations and add local attribute and target-access guards. Constraint
  restoration and its guards belong to S4a2.
  Implement private Molecule::restore_topology with
  `(&mut self, &GraphCompaction, Vec<RemovedAtom>, Vec<RemovedBond>) -> ()`.
  Its scope is Graph plus the atom/bond attribute vectors. S4a2 removes the editor's
  combined topology-and-overlay restoration body; undo execution composes this
  primitive with separate private Molecule delegates to each set's existing
  restore_topology_ids and restore, followed by restore_constraints.
  Overlay-only undo calls the affected row-restoration delegate without
  topology-id restoration. Keep these as independent operations.

  Each row defines two private Molecule methods. The row-restoration method
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

- **S4a2 — Constraint storage and undo restoration wiring** (`ir::constraint::molecule`, molecule/transact; breaking, red→green). [dep: S4a1]

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
  Add private Molecule
  `restore_constraints(&mut self, changes: &CascadedConstraints) -> ()`,
  delegating only to Constraints::restore. These delegates do not add old-value
  checks, journals, or publication checks. Migrate undo callers to them; preserve
  the storage restoration contracts, including panic freedom for manipulated
  history. Test topology-id restoration and row restoration together and
  independently where topology is unchanged, covering all six overlay kinds.
  Keep the existing detached journal surface until S5. Test matching-history
  recovery under `normalized_eq` and panic freedom for manipulated undo data
  without asserting its result. Cover constraint compaction/restoration with
  removed and rewritten entries, duplicate entries, and preserved list order.

- **S4b — Molecule delegation and Edit/Undo execution** (group; breaking,
  green at S4b8). [dep: S4a, S3i]

  Execute S4b1–S4b8 in order. The complete editor inventory is in S4b2;
  the execution signatures are in S4b3 and every Edit/Undo call is mapped below.
  Changes to journal return types and Undo variants may temporarily break callers
  within this group; migrate all of them by S4b8 without fallback match arms.

- **S4b1 — Molecule addition, removal, and compaction methods** (`ir::molecule`; breaking rewire, green at S4b8). [dep: S4a, S3i]

  Move the eight add_* methods,
  push_constraint, and all seven untracked/tracked removal pairs to private
  Molecule methods, with the interfaces recorded under Storage delegation.
  Editor additions delegate; complete editor removals and batch execution compose
  the same primitives.
  Remove MoleculeEditor::push_constraint and migrate its callers to
  editor.constraints_mut().push(constraint). Retain private
  Molecule::push_constraint for Edit execution. Do not replace the removed editor
  method with another top-level constraint convenience delegate.
  Reuse S3i's bulk add_atoms/add_bonds for the corresponding Edit variants;
  preserve individual overlay Edit execution and its existing undo boundaries.
  Rename editor remove to remove_topology and private tracked_remove to
  tracked_remove_topology, migrating callers. The private topology pair mutates
  only Graph and atom/bond attributes and returns GraphCompaction when tracked;
  each overlay pair mutates only its owning set and returns Compaction<Id> when
  tracked. Public editor removal still completes all cascading updates.

  Add the following private Molecule compaction pairs. Each method takes
  `(&mut self, topology: &GraphCompaction)`. The bare method returns (); the
  tracked method returns the listed mapping. Both install the replacement set
  returned by that owning set's compact/tracked_compact. They mutate no other
  component and perform no extra validation or journal recording.

  | Bare method | Tracked method | Tracked return |
  | --- | --- | --- |
  | compact_dative_bonds | tracked_compact_dative_bonds | Compaction<DativeBondId> |
  | compact_aromatic_systems | tracked_compact_aromatic_systems | Compaction<AromaticSystemId> |
  | compact_multicenter_bonds | tracked_compact_multicenter_bonds | Compaction<MulticenterBondId> |
  | compact_noncovalent_bonds | tracked_compact_noncovalent_bonds | Compaction<NoncovalentBondId> |
  | compact_stereo_atoms | tracked_compact_stereo_atoms | Compaction<StereoAtomId> |
  | compact_stereo_bonds | tracked_compact_stereo_bonds | Compaction<StereoBondId> |

  Complete the private Molecule constraint delegation surface:

  ```rust
  fn push_constraint(&mut self, constraint: Constraint);
  fn extend_constraints(&mut self, constraints: Vec<Constraint>);
  fn remove_constraint_at(&mut self, position: usize) -> Constraint;
  fn compact_constraints(&mut self, compaction: &MoleculeCompaction);
  fn tracked_compact_constraints(
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

- **S4b2 — Editor direct-mutation surface** (`ir::molecule::editor`; breaking rewire, green at S4b8). [dep: S4b1]

  **Complete editor direct-mutation inventory.** The following is the resulting
  public surface, including mutation reached through returned views/borrows.
  Addition signatures are in Storage delegation; all take &mut self and
  delegate to synonymous private Molecule methods.

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
  returns the corresponding *ViewMut<'_, true> obtained from Molecule (S2i).
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

- **S4b3 — Single-entry execution and batch ownership** (`ir::molecule::transact`; breaking, green at S4b8). [dep: S4b2]

  Move apply_edit, apply_edit_with_undo, and apply_undo from impl MoleculeEditor
  to impl Molecule in ir::molecule::transact, retaining their private visibility
  and names. Both forward methods take one Edit and the batch's mutable
  ApplicationState; apply_edit_with_undo returns an optional Undo as specified
  below. apply_undo takes one Undo. Initialize ApplicationState from Molecule
  accessors rather than an editor. Editor and Transaction retain their batch
  loops; Transaction retains journal storage and reverse replay. Add no execution
  methods to Edit, Edits, or Undo and no separate execution type. These Molecule
  methods obey the same accessor/mutation-method boundary as the inventory below;
  moving their impl does not authorize direct storage access.

  ```rust
  // Private methods in ir::molecule::transact; editor batch bodies stay here too.
  fn apply_edit(&mut self, edit: Edit, state: &mut ApplicationState)
      -> Result<(), TransactionError>;
  fn apply_edit_with_undo(&mut self, edit: Edit, state: &mut ApplicationState)
      -> Result<Option<Undo>, TransactionError>;
  fn apply_undo(&mut self, undo: Undo);
  // Existing private batch state, with a changed receiver source:
  ApplicationState::new(molecule: &Molecule) -> ApplicationState;
  ```

  Execution returns the existing per-Edit error category. Public lifecycle
  methods wrap it in MoleculeApplyError::Transaction; integrity belongs to
  probe/finish/commit. Undo has no Result or expected-post-value check. Keeping
  execution and editor batch bodies in the same module makes these private
  methods accessible without widening visibility.

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

  After forward topology removal, call all six Molecule tracked_compact_* overlay
  delegates and assemble MoleculeCompaction. After explicit overlay removal,
  assemble that mapping with identity components for unchanged kinds. Then
  call compact_constraints or tracked_compact_constraints. Addition undo instead
  calls only the private Molecule untracked removal, as listed below. Remove the
  seven remove_added_* adapters and make those calls in the Undo match arms,
  without overlay/constraint compaction or
  correspondence updates. Keep forward sequencing explicit; add no combined
  coordination helper or recording-mode parameter solely to avoid repetition.
  Remove all sixteen apply_modify_*_field / apply_modify_*_constraint helpers.
  Keep handle resolution and old-value/precondition checks in batch execution,
  and perform field and entity-constraint writes through mutable views obtained
  from Molecule. Top-level constraint changes use the private Molecule delegates.
  Structural Edit replacements and their undos use the <true> view methods wired
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
  receiver are private Molecule methods; m is the owned or borrowed Molecule.
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

- **S4b4 — Addition execution** (`ir::molecule::transact`; rewire, green at S4b8). [dep: S4b3]

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
  are private Molecule primitives, not the editor's complete cascading removal
  operations. Reverse replay has undone later references and compactions, so
  matching added entries again occupy trailing positions and earlier ids stay
  unchanged. Do not compact overlays or constraints, assemble MoleculeCompaction,
  update correspondence, or create new undo. The Undo match arms
  only extract ids and call these methods, with local guards for panic freedom
  on manipulated input; they neither validate history nor inspect saved payloads.

  Verify ids, handle registration, saved entries, and reverse removal after
  later dependent edits have been undone.

- **S4b5 — Removal and cascading compaction execution** (`ir::molecule::transact`; rewire, green at S4b8). [dep: S4b4]

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
  | tracked_compact_dative_bonds | restore_dative_bond_topology_ids | restore_dative_bonds |
  | tracked_compact_aromatic_systems | restore_aromatic_system_topology_ids | restore_aromatic_systems |
  | tracked_compact_multicenter_bonds | restore_multicenter_bond_topology_ids | restore_multicenter_bonds |
  | tracked_compact_noncovalent_bonds | restore_noncovalent_bond_topology_ids | restore_noncovalent_bonds |
  | tracked_compact_stereo_atoms | restore_stereo_atom_topology_ids | restore_stereo_atoms |
  | tracked_compact_stereo_bonds | restore_stereo_bond_topology_ids | restore_stereo_bonds |

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

- **S4b6 — Attribute edit execution** (`ir::molecule::transact`; rewire, green at S4b8). [dep: S4b5]

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

- **S4b7 — Constraint edits and explicit constraint Undo variants** (`ir::{edit,molecule::transact}`; breaking, green at S4b8). [dep: S4b6]

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
  | RemoveMoleculeConstraint | Find the last exact-equal entry through m.constraints(); MissingEntry if absent | remove_constraint_at(position) | ApplyCascadedConstraints containing RemovedConstraint with that position and the returned stored value. Replay calls restore_constraints(&changes). |

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

- **S4b8 — Undo replay and migration closure** (`ir::molecule::transact`, Rust consumers; breaking, red→green). [dep: S4b7]

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
  replay performs no compaction algorithm, correspondence update, offered-old
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

- **S4c — moved to S5a1.** Introduce the scope guard with Transaction::run,
  which owns it. No unused recovery machinery or temporary public API is added
  at the S4 green boundary.

### S5 — Borrowed transaction API

S5a removes APIs used by S5b–S5d; those migrations are required before the stage
returns green. S5d's Python invalidation and sequential input-consumption contracts are recorded above.

- **S5a — Scoped Rust transaction lifecycle** (group; breaking, green at
  S5d3). [dep: S4b]

- **S5a1 — Guard, borrowed handle, and scoped run**
  (`ir::molecule::transact`; breaking, green at S5d3). [dep: S4b]

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

- **S5a2 — Batch application, completion, and Molecule conveniences**
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

- **S5b** (`umol-graph-ir` reaction and molecule callers; breaking, green at S5d)
  Migrate uses of detached journals. Preserve reaction `Ok(None)` versus error
  classification and host-to-product correspondence.
  Test failed applications and product integrity. [dep: S5a]
- **S5c** (`umol-graph::ops`; breaking, green at S5d) Migrate existing borrowed
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
- **S5d — Python prepared transactions and ownership** (group; breaking,
  closes S5 green at S5d3). [dep: S5a, S5b, S5c]

  Execute S5d1–S5d3 in order. No interactive Python transaction handle is added.

- **S5d1 — Edits consumption and iterator access** (`umol-py::edit` and its consumers; breaking, green at S5d3). [dep: S5a, S5b, S5c]

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

- **S5d2 — Molecule accessor counters and Storage names** (`umol-py::molecule`, entity/constraint/collection views; breaking, green at S5d3). [dep: S5d1]

  Implement the [approved accessor invalidation contract](#python-accessor-invalidation).
  Current AtomView stores an owner and a dense id; deleting an earlier atom can make it
  silently address another atom. Owner consumption is not involved in transact,
  so S6d's consumed-state checks alone cannot fix this. Apply the specified
  triggers, failure behavior, nested access, and *Storage names, including the
  later *_into consumers. Whole-molecule mutation invalidates all owner-backed
  accessors at execution entry, including on no-op and rollback; ordinary view
  setters retain access.

  Apply the complete field and method table under Python accessor invalidation.
  Test all eight entity kinds, nested constraints and ring sizes, exhausted
  iterators, and ordinary setter usability. Molecule storage is not Option yet;
  S6d adds consumption without changing these counter rules.

- **S5d3 — Prepared-batch transaction bindings** (`umol-py::{molecule,transaction}`; breaking, red→green). [dep: S5d2]

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
  fn transact(slf: Py<Self>, py: Python<'_>, batches: Vec<Py<Edits>>)
      -> PyResult<()>;
  fn tracked_transact(slf: Py<Self>, py: Python<'_>, batches: Vec<Py<Edits>>)
      -> PyResult<MoleculeCorrespondence>;
  ```

  These are methods of the Python Molecule wrapper. Finish Python iterable
  extraction before taking batches. Use a short try_borrow_mut for each Edits
  transfer; after all transfers, exclusively borrow the molecule, advance its
  counter, and call the Rust operation. A later transfer, receiver-borrow, or
  counter-overflow error does not restore already consumed batches. No Python
  iteration/callback occurs during Rust application. The returned correspondence
  is the existing Python wrapper, constructed from the Rust result.

  Remove the detached Python Transaction class and editor transact methods.
  Test independent New namespaces, no-op invalidation, rollback after a later
  failure, tracked correspondence, and no invalidation before execution starts.

### S6 — Owning editor and publication

S6a–S6d are one ownership migration, returning green after S6d. Migrate every
consumer of the changed signatures, including tests and benchmarks; retaining
temporary cloning adapters is not a way to close an earlier subitem.

- **S6a** (`ir::molecule`, `ir::molecule::editor`; breaking, green at S6d) Make
  `edit(self)` own its Molecule, make Molecule `apply`/`tracked_apply` consume,
  and replace editor snapshot/build methods with checked `probe`/`finish`.
  Retain MoleculeBuilder's asserted build and the S2 uniform mutable attribute access.
  Move the owned constraint collection with the molecule; no clone is needed.
  Test direct/batch interleaving, invalid probe then repair, destructive failure, and
  constructor-equivalent publication.
  [dep: S2m, S5d]
- **S6b — Graph-IR ownership migration** (group; breaking, green at S6d).
  [dep: S6a]

- **S6b1 — Graph-IR callers and correspondence** (`umol-graph-ir` reaction/molecule callers; breaking, green at S6d). [dep: S6a]

  Migrate editor construction/publication and remove session
  correspondence accumulation and public tracked direct removal. Use one
  source-preserving product candidate where the host survives; keep batch-only
  tracking distinct from whole-transaction tracking. Test product and
  correspondence laws.
  Keep reaction product failure classification and one intentional host copy.
  A rejected product candidate is dropped; it needs no recovery transaction.

- **S6b2 — Disjoint append on a borrowed molecule** (`ir::molecule`; breaking rewire, green at S6d). [dep: S6b1, S4b1]

  Migrate combine_from's current mem::take/from_parts
  route here because the owning editor change affects it now, not in S8.
  Preserve its append semantics and borrowed receiver without an initial recovery
  clone. No append/count checkpoints or recovery tests for internal panics between
  storage writes are required. This moves the former S8c migration here. The
  separate lift_constraints correction remains excluded.

  Keep combine_from(&mut self, other: &Molecule) -> () and its existing
  disjoint-append semantics. Save the eight original entity counts solely as
  offsets for other, then use private Molecule additions directly: atoms,
  bonds with shifted endpoints, and the six overlay kinds with shifted sites
  and atoms. Retain the existing correspondence-based ligand and constraint
  translation; append mapped constraints through extend_constraints. Only
  other's borrowed payloads are copied. Do not move self into an editor, call
  mem::take(self), clone self for recovery, or retain offsets as rollback
  checkpoints. Finish with the existing asserted producer integrity guarantee;
  no new public method or error type is needed.

  Verify original-id prefixes, every overlay factor, inline/top-level constraints,
  empty operands, shared source storage, and the unchanged other molecule.

- **S6c — Chemistry and format caller migration** (group; breaking, green
  at S6d). [dep: S6a, S5c]

- **S6c1 — Transformation plans and borrowed execution** (`umol-graph::ops`; breaking rewire, green at S6d). [dep: S6a, S5c]

  Migrate graph operation publication sites to `finish` or borrowed transactions
  without changing
  their chemistry or boundary outcomes. This includes AromaticityPerceiver's
  add_systems and the three existing transform_into implementations. Construct
  their batches and shared chemistry checks according to the transformation
  mapping now: borrowed callers cannot simply consume their receiver, and the
  DelocalizeCharge editor path introduced in S2k cannot wait until S8 to migrate.
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

- **S6c2 — Remaining ingest, parse, and export callers** (`umol-graph`, `umol-io`; breaking, green at S6d). [dep: S6c1]

  Migrate remaining edit/build/snapshot call sites to the owning editor and
  finish while retaining their current public signatures. Preserve each source
  copy already required for an independent output; do not add recovery copies.
  Resolver names and report removal change in S7, not here. Verify boundary
  outputs, error categories, and unchanged retained sources.

- **S6d — Python Molecule ownership migration** (group; breaking, closes
  S6 green at S6d3). [dep: S5d, S6a, S6b, S6c]

- **S6d1 — Molecule Option storage and accessors** (`umol-py::molecule`; breaking, green at S6d3). [dep: S5d, S6a, S6b, S6c]

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

- **S6d2 — Fallible owner access throughout bindings** (Molecule/view consumers in `umol-py`; breaking, green at S6d3). [dep: S6d1]

  Migrate every read/write of the wrapper through its fallible accessors,
  including getters, repr/equality, conversions, collections, entity attributes,
  and nested constraints. Reuse S5d2 check_access and counter propagation.
  Do not convert unrelated forms to Option or add automatic clones. Verify
  consumed roots and invalidated children across every entity collection.

- **S6d3 — Consuming entry points and editor publication** (`umol-py::{molecule,transaction,edit}`; breaking, red→green). [dep: S6d2]

  Transfer Molecule inputs on consuming calls and reuse S5d's Edits
  transfer; raise `ConsumedError` for the owner and `InvalidatedViewError` for
  owner-backed accessors. Replace Python snapshot/build with finish, without a blanket
  consumed-state migration of unrelated forms. Test aliases, nested views,
  successful and failed consumption, and explicit copies. Apply fallible owner
  access throughout the binding consumers, not only mutation entry points;
  getters, repr/equality, conversions, collections, and nested setters must not
  bypass the consumed/invalidation checks.

### S7 — Resolution, projection, and boundaries

S7a–S7d form one public signature/result migration, returning green at S7d.

- **S7a** (`umol-graph::ops::resolve` and phase modules; breaking, green at S7d)
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
  and the source-preserving export path. Review the complete diff against 213,
  run formatting, workspace tests, strict Clippy and rustdoc, explicit feature
  suites, Python 3.13 build/tests, and the pinned Rust 1.87 gate once at final
  closeout. Record results and update the discussion status only after the full
  scope passes. [dep: S6b, S9a]

**Execution order:** completed S0a–S0b → completed S1a–S1c → S2 → S3 →
S4 → S5 → S6 → S7/S8 → S9. S2a is implemented; its API is removed in S2i1.
Within the revised S2:

- S2b and S2c are complete. S2d is next; S2f is cancelled.
  S2c → S2d establishes coset-operation and count-use behavior. S2e is folded
  into S2h; it is not an executable prerequisite.
- S2g depends on S2c; its interfaces and failure behavior are approved.
- S2h removes aggregate attribute checks only after S2c/S2d/S2g
  consumer changes are complete, then tests the newly admissible inputs through
  public construction and those consumers before closing the subitem.
- S2i → S2j supplies the shared const-generic views and editor structural
  operations after S2h; attribute mutation is unrestricted on both specializations.
  S2i1 closes green before the shared-view migration; S2i2–S2i4 close at S2i4.
  S2j closes its own getter/structural migration green;
  neither waits for S2k/S2l to restore compilation.
- S2k and S2l migrate Rust and Python callers; both precede removal in S2m.
- S3a depends on S2i; S3b requires S2j. S3c is independent of view changes.
- S3f and S3g supply graph-core bulk additions; S3g → S3h supplies typed-set
  extend, then S3f/S3h → S3i supplies Molecule/editor bulk additions. S4b uses
  those additions and the component removal/restoration interfaces.
- S4a closes at S4a2; S4b closes at S4b8. S4c is incorporated in S5a1.
- S5a1–S5a2 introduce the guard and public lifecycle together; S5d1–S5d3
  complete Python ownership, counters, and prepared transactions.
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
