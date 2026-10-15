# 181 — Python boundary ownership and copy costs

Status: In Progress
Date: 2026-08-03
Relates: [179](179-python-editing-and-transactions-2026-08-02.md),
[151](151-python-molecule-workflows-2026-07-13.md),
[199](199-open-container-integrity-2026-08-18.md),
[213](213-editor-overlay-storage-2026-08-27.md),
[231](231-python-api-parity-2026-09-21.md),
[236](236-operation-ownership-copy-costs-2026-10-09.md),
[Python API guide](../docs/development/python-api.md),
[data-type guide](../docs/development/data-types.md)

## Molecule lifecycle

Discussion scope agreed 2026-10-08. The eight sections below identify the decisions
and measurements to work through; they are not an implementation plan.

The consuming molecule API can leave the original Python object unusable after
failure. Avoiding that requires the caller to copy in advance. This lifecycle is
the release blocker; changing only the recovery mechanism inside resolve_into
does not resolve it.

[213's recovery measurements](213-editor-overlay-storage-2026-08-27.md#addendum--recovery-strategy-handoff-2026-09-29)
found clone-based recovery faster and smaller than journaling on the tested
resolution fixtures, with little extra cost over consuming execution. These
measurements include copy-on-write detachment but do not cover failed resolution
or establish the best mechanism for every operation.

## 1. Consuming Rust APIs, including traits

### Current inventory — 2026-10-08

Checked against code at ec82fcdcc. These are implemented interfaces and caller
relationships, not proposed replacements. Existing simple attribute mutation is
outside this ownership review.

### Molecule and editor

Signatures below omit method bodies only. Correspondence in tracked editor
application covers that batch, not all preceding editor changes.

```rust
// Molecule
fn edit(self) -> MoleculeEditor;
fn apply(self, edits: Edits) -> Result<Self, MoleculeApplyError>;
fn tracked_apply(self, edits: Edits)
    -> Result<(Self, MoleculeCorrespondence), MoleculeApplyError>;

// MoleculeEditor
fn apply(self, edits: Edits) -> Result<Self, MoleculeApplyError>;
fn tracked_apply(self, edits: Edits)
    -> Result<(Self, MoleculeCorrespondence), MoleculeApplyError>;
fn probe(&self) -> Result<&Molecule, MoleculeIntegrityError>;
fn finish(self) -> Result<Molecule, MoleculeIntegrityError>;
```

| Operation | Success | Failure / abandonment |
| --- | --- | --- |
| Molecule edit | Moves storage into an owning editor; no recovery copy or journal | Dropping the editor drops its molecule |
| Editor direct mutation | Changes the owned molecule through mutable editor access; may leave intermediate state | No rollback |
| Editor apply / tracked_apply | Checks edit preconditions and returns a transient editor; no integrity check or journal | Drops the editor, including earlier changes; returns no molecule |
| Editor probe | Checks integrity and borrows the current molecule without copying | Keeps the editor available |
| Editor finish | Checks integrity and moves out the molecule | Drops the editor and molecule |
| Molecule apply / tracked_apply | Calls edit, editor application, then finish | Drops input on application or integrity failure |

Python edit consumes the Molecule owner. Python editor apply consumes the old
editor and returns a new editor object even on success. Both molecule and editor
application consume their receiver before transferring Edits; a batch-transfer
failure therefore also loses the receiver. Python finish consumes the editor on
success and failure. These are separate lifecycle decisions, not one necessary
consequence of using an editor.

Sources: [Molecule](../umol-graph-ir/src/ir/molecule.rs),
[editor](../umol-graph-ir/src/ir/molecule/editor.rs),
[editor application](../umol-graph-ir/src/ir/molecule/editor/transact.rs),
[Python molecule](../umol-py/src/molecule.rs),
[Python editor](../umol-py/src/transaction.rs).

### Existing borrowed recovery interfaces

Molecule transact and tracked_transact take &mut self and an IntoIterator of
owned Edits. They return Result<(), MoleculeApplyError> and
Result<MoleculeCorrespondence, MoleculeApplyError>, respectively. They apply all
batches and commit once.

Transaction::run borrows &mut Molecule and lends Transaction to a callback.
Transaction apply borrows the handle mutably and consumes one Edits batch.
Commit, tracked_commit, and rollback consume the handle, not the molecule.
Commit checks integrity; keeping changes also requires the callback to return
Ok. Application failure, commit failure, cancellation, or callback unwinding
between completed edits restores transaction-entry state under normalized_eq,
preserving participant order. These routes currently record undo.

Python exposes prepared-batch transact and tracked_transact, not the scoped
Transaction handle. The molecule remains usable after failure; batches already
transferred remain consumed. Existing views are invalidated before execution,
including when execution subsequently restores the molecule.

Source: [transaction](../umol-graph-ir/src/ir/molecule/transact.rs).

### Resolution and projection

All listed resolver methods borrow the resolver itself. The ownership distinction
concerns their molecule argument.

| Resolver | Consuming methods | Borrowed counterparts |
| --- | --- | --- |
| Resolver | resolve, resolve_with_report, project | resolve_into, resolve_into_with_report, project_into |
| IsotopeResolver | resolve, project | resolve_into, project_into |
| AromaticityResolver | resolve, project | resolve_into, project_into |
| StereoResolver | resolve, project | resolve_into, project_into |
| MulticenterBondsResolver | resolve | resolve_into |
| BondsResolver | resolve | resolve_into |

Consuming resolve/project take Molecule and return
Result<Solution<Molecule, C, ()>, E>, with each resolver's contradiction and error
types. Composite project additionally takes ProjectFlags. Borrowed counterparts
take &mut Molecule and return Result<Solution<(), C>, E>.

- Consuming execution uses an editor without recovery. Only Determined returns a
  molecule; Underdetermined, Contradictory, and Err do not return the input.
- Borrowed execution uses transactions and retains changes only on Determined.
  Composite resolution/projection share one transaction across their phases;
  individual resolvers plan against the input and transact their batch.
- resolve_with_report returns
  Result<Solution<(Molecule, ResolveReport), ResolveContradiction, ResolveReport>, ResolveError>.
  Its borrowed counterpart returns
  Result<Solution<ResolveReport, ResolveContradiction, ResolveReport>, ResolveError>.
  Reports do not provide recovery of the consumed molecule.
- ValenceResolver is an exception to this method family: project already takes
  &mut Molecule plus ValenceTieBreak, returns
  Result<Solution<(), ValenceContradiction>, ValenceProjectError>, and is a no-op.
  It has no public consuming resolve counterpart; composite resolution uses its
  admission and planning operations.

Python Molecule exposes the four composite resolve methods, transferring the
molecule for resolve and resolve_with_report and borrowing it for the two
_into methods. There are no direct Python projector or individual-resolver
bindings in this surface.

Sources: [composite resolver](../umol-graph/src/ops/resolve.rs),
[individual resolvers](../umol-graph/src/ops/resolve/),
[Python molecule](../umol-py/src/molecule.rs).

### Transformer trait

```rust
fn transform(&self, molecule: Molecule) -> Result<Molecule, Self::Error>;
fn transform_into(&self, molecule: &mut Molecule) -> Result<(), Self::Error>;
fn transform_iter<'a>(&'a self, molecule: &'a Molecule)
    -> impl Iterator<Item = Molecule> + 'a;
```

The implementations are Kekulizer, Aromatizer, and DelocalizeCharge. transform
uses journal-free editor application; transform_into uses transactions.
transform_iter lazily clones its source and invokes transform, yielding one
independent result on success and none on rejection. DelocalizeCharge has
Infallible as its error type. No direct Python Transformer binding is present.

Source: [Transformer and its implementations](../umol-graph/src/ops/transform.rs).

### Value-operation traits

| Trait | Consuming interface | Existing borrowed interface |
| --- | --- | --- |
| Normalize | normalize(self) -> Result<Self, Contradiction> | normalized(&self) -> Result<Cow<Self>, Contradiction>; normalized_eq |
| FrameTransport | reframe_by(self, &Self::Action) -> Option<Self> | No borrowed transport counterpart |
| Reframe | reframe(self) -> Result<Self, Contradiction>; tracked_reframe(self) -> Result<(Self, Self::Action), Contradiction> | representative_action; framed_eq |
| Canonicalize | canonicalize(self, &CanonicalizeContext) -> Result<Self, Self::Error>; tracked_canonicalize returning (Self, MoleculeRemapping); canonical_hash returning u64 | canonical_eq |

Failure consumes the value in every consuming entry above. canonical_hash also
consumes it on success although its result contains no aggregate.

Normalize and FrameTransport reach well beyond Molecule: forms, constraints,
entity entries and storage, spans, deltas, and reaction representations use them.
Reframe applies to entity entries/sets/spans and the three aggregates Molecule,
Reaction, and ReactionSpan. Canonicalize has exactly those three aggregate
implementations. Changing these traits is therefore broader than changing the
molecule editor.

There are concrete implementation differences:

- Normalize's default normalized clones and calls normalize; Molecule does not
  override it. Molecule normalize moves fields and uses copy-on-write mutable
  access for its atom/bond tables.
- Frame transport and reframing move and transform owned fields. Default
  framed_eq clones its inputs when structural equality does not suffice.
- All three aggregate canonicalize and tracked_canonicalize implementations
  accept self but immediately pass &self to internal canonicalization routines.
  Their consuming signatures do not currently enable an ownership transfer into
  those routines.
- Python canonicalize and tracked_canonicalize already preserve their source.
  Molecule and ReactionSpan bindings clone; Reaction converts to an owned Rust
  value. Python form normalization and Deltas normalization also preserve their
  receiver. The aggregate wrappers do not expose normalize, reframe, reframe_by,
  or canonical_hash methods.

Sources: [traits](../umol-graph-ir/src/ir/traits.rs),
[aggregate canonicalization](../umol-graph-ir/src/ir/canonicalize.rs),
[Python canonicalization](../umol-py/src/canonicalize.rs),
[Python normalization](../umol-py/src/lattice.rs).

### Production caller relationships

| Caller | Current ownership path |
| --- | --- |
| SMILES ingestion and MOL parsing | Raise a new molecule, then consume it through composite resolve; unresolved outcomes become input errors |
| Molecule export | Clone the borrowed source, then consume that copy through project |
| Reaction export | Obtain side molecules, then consume them through project |
| Supplied-match reaction application | Clone the retained host, then use Molecule apply to produce an independent product |
| Composite resolution/projection | Apply successive planned batches to one editor, probing between phases and finishing once; borrowed counterparts use one transaction |
| Individual resolvers and transformers | Plan edits, then use editor application or transact according to the public ownership route |
| Molecule extraction and construction | Extraction clones before editing; builder and split use editors for newly assembled values |
| Python molecule/editor application | Transfer the owned Rust value out of its Option and wrap the successful result in a new Python owner |
| Python canonicalization | Preserve Python ownership and supply an owned Rust copy/conversion to the consuming trait |

Sources: [ingest](../umol-graph/src/ingest.rs),
[MOL parsing](../umol-graph/src/parse.rs),
[export](../umol-graph/src/export.rs),
[reaction application](../umol-graph-ir/src/ir/reaction.rs),
[molecule construction/extraction](../umol-graph-ir/src/ir/molecule.rs),
[builder](../umol-graph-ir/src/ir/molecule/build.rs).

### Proposed Molecule / MoleculeEditor API — 2026-10-08

The editor owns one working molecule in every route. It has no copying/journaling
mode and does not retain a destination borrow. Rust and Python lifecycles can be
decided independently; this proposal concerns Rust.

#### Entry and completion

```rust
// Molecule: move the molecule into an editor, or apply a batch to it.
fn edit(self) -> MoleculeEditor;
fn apply(self, edits: Edits) -> Result<Molecule, MoleculeApplyError>;
fn tracked_apply(self, edits: Edits)
    -> Result<(Molecule, MoleculeCorrespondence), MoleculeApplyError>;

// MoleculeEditor
fn apply(&mut self, edits: Edits) -> Result<(), MoleculeApplyError>;
fn tracked_apply(&mut self, edits: Edits)
    -> Result<MoleculeCorrespondence, MoleculeApplyError>;
fn probe(&self) -> Result<&Molecule, MoleculeIntegrityError>;
fn finish(self) -> Result<Molecule, MoleculeIntegrityError>;
```

edit, apply, and tracked_apply consume the molecule, as they do today;
[236](236-operation-ownership-copy-costs-2026-10-09.md) owns their ownership contract. edit moves the storage into the editor
without copying; a caller that keeps the molecule clones it first. No route
records undo. Molecule apply and tracked_apply are conveniences over edit,
batch application, and finish; on failure the consumed input is dropped.

Editor apply checks edit preconditions but provides no rollback. It mutably
borrows the editor rather than consuming it. Failure leaves completed changes
in the editor, which remains available for inspection, further mutation, or
discard. Direct mutation likewise permits intermediate states. probe checks
integrity before lending the molecule; failure retains the editor. finish checks
integrity before moving out the molecule; failure drops the working molecule.

#### In-place mutation

Add source-preserving execution that installs the successful working molecule
into a mutably borrowed receiver:

```rust
molecule.edit_into(|editor| {
    editor.atom_mut(id).attributes_mut().charge = charge;
    editor.apply(edits)?;
    Ok(())
})?;

molecule.apply_into(edits)?;
```

edit_into and apply_into are settled. edit_into clones the receiver once,
lends mutable editor access to the callback, and checks integrity after a
successful callback return. Only then does it replace the receiver. Callback
error, integrity failure, or unwinding discards the working molecule and leaves
the receiver unchanged. apply_into is the prepared-batch convenience with the
same receiver guarantee. The caller does not call finish inside the callback.
The settled generic callback signature is given in the migration interface review.

This is in-place mutation at the public API: the existing receiver changes.
Cloning and installation are its implementation, not a different public
ownership contract. Preserving the receiver by discarding a failed working copy
does not give editor apply rollback semantics.

#### Tracking

Editor tracked_apply returns the correspondence from immediately before its
batch to immediately afterward. The editor does not accumulate it. Callers can
compose correspondences from consecutive successful batches. Failure returns no
correspondence and leaves partial editor changes, just as apply does. A mapping
returned before that failed batch does not describe the resulting partial state.
Molecule tracked_apply returns its correspondence only after successful finish.

#### No-recovery route

```rust
let mut editor = molecule.edit();
editor.atom_mut(id).attributes_mut().charge = charge;
editor.apply(edits)?;
let correspondence = editor.tracked_apply(more_edits)?;
let molecule = editor.finish()?;
```

This route makes no recovery copy and records no journal. Dropping the editor
loses the original molecule. Batch failure retains the partial editor only if
the caller handles the error rather than exiting; publication still requires
integrity. No-recovery does not mean publishing an integrity-invalid Molecule.

The transaction API is retained with the proposed implementation changes in
point 2. Value-trait ownership remains a separate decision.

### Resolution / projection / transformation API

[236](236-operation-ownership-copy-costs-2026-10-09.md) settles the naming and ownership. The bare verb consumes its molecule;
`_into` mutates a borrowed molecule with recovery; a caller that keeps its
molecule clones it explicitly. Python binds the bare verb as clone followed by
the Rust method.

| Operation | Molecule ownership and result |
| --- | --- |
| resolve / resolve_with_report / project | Consume the molecule; return it on Determined |
| transform | Consume the molecule; return the transformed molecule on success |
| composite resolve_into / resolve_into_with_report / project_into | Mutable borrow with RollbackMethod last; retain changes only on Determined |
| individual resolvers' resolve_into / project_into, transform_into | Mutable borrow; failures are decided during planning, so the plan is applied in place without a selector |
| transform_iter | Borrow the source and lazily produce independent results |

Existing result and diagnostic types are retained. Ordinary resolve does not
construct a report. Ingestion, parsing, and export pass their owned molecules to
the consuming methods without another copy or recovery journal.
ValenceResolver's existing no-op project remains unchanged.

### Other inventory boundaries

Molecule remap/try_remap, extract, combine, and split already borrow their inputs
and produce results. Reaction application also preserves its host. They belong
in the copy-producing review in point 3, not a list of consuming molecule methods.
Builder completion and consuming conversions from construction data are distinct
from losing an existing molecule and are not proposed for redesign here.

The normative guides must follow the settled decisions: data-types.md and
python-api.md describe consuming mutation, while nomenclature.md's Resolution
entry still describes source-preserving resolution without distinguishing the
consuming route. No normative contract is changed by this inventory.

## 2. Mutating APIs: journals versus cloning

### Transaction interface and rollback selection — revised 2026-10-14

Select tracking at transaction entry through separate scoped handle types.
Transaction exposes run, apply, probe, commit and rollback; TrackedTransaction
exposes the same method names. Transaction::commit returns unit, while
TrackedTransaction::commit returns MoleculeCorrespondence. Remove
Transaction::tracked_commit. Molecule transact/tracked_transact retain their
plain/tracked results and select the corresponding handle internally.

Tracking and recovery are independent choices. Both types provide both recovery
implementations through an explicit selector:

```rust
pub enum RollbackMethod {
    Clone,
    Journal,
}
```

Transaction::run and TrackedTransaction::run take RollbackMethod after the
molecule argument and before the callback. Selection applies to the entire
transaction, not to each batch or phase.
Clone uses one working copy and installs it on success; Journal records undo
while modifying the receiver. Both preserve the original on failure and retain
the same callback, commit/rollback, integrity-checking, error, and correspondence
contracts.

Higher-level operations select a default and allow callers to override it. This
is an execution option, separate from chemistry parameters, following the
explicit algorithm-selection pattern in graph-core. The outer operation owns
the selection; phase planners do not open independent recovery scopes. No
automatic selection heuristic or universal winner is proposed. The composite
_into methods take RollbackMethod last; Python resolve_into and
resolve_into_with_report default to Clone.

```rust
Transaction::run(&mut molecule, RollbackMethod::Clone, |mut tx| {
    tx.apply(first_edits)?;
    let next_edits = plan(tx.probe()?)?;
    tx.apply(next_edits)?;
    tx.commit()
})?;
```

The Clone implementation executes as follows:

1. run clones the receiver once into working storage, leaving the receiver
   unchanged throughout the callback.
2. apply executes each batch immediately on that working molecule. Each batch
   retains its own handle namespace.
3. probe checks integrity before lending a shared molecule borrow.
4. commit consumes the handle, checks integrity, and requests acceptance.
5. Only after a successful callback return does run install an accepted result
   into the receiver.

An application error aborts the transaction; subsequent application or commit
cannot succeed. Failed commit, callback error, explicit rollback, cancellation
without commit, or callback unwinding discards the working copy. Forgetting the
handle cannot install an unaccepted result: run retains responsibility for
publication. No undo replay is needed to preserve the receiver.

### Editor relationship and convenience methods

Both clone-based transaction types hold their working molecule in
TransactionStorage::Clone(Option<Molecule>) and apply batches through the
molecule's edit application methods. Neither contains an editor. Both expose
batches and checked observation, not direct mutable entity views.
Editor apply failure retains partial working state for further editing;
transaction apply failure aborts the entire transaction.

edit_into supports direct mutation followed by checked publication. Transaction
supports edit batches with explicit commit/rollback and abort-on-application-error
semantics. Both preserve the original receiver on failure. Molecule transact
continues to apply prepared batches and commit once; tracked_transact returns
the corresponding whole-operation mapping.

### Correspondence and implementation requirements

Clone-based execution replaces receiver mutation and undo replay with one
working copy, discarded on failure and installed only after commit and callback
success. Transaction states, error behavior, checked probe access, and batch
handle namespaces stay the same.

TrackedTransaction maintains one correspondence from transaction entry to the
state at commit, independently of Undo, with either recovery method. Initialize
it to identity at tracked run entry and update it through ApplicationState after
each successful batch. TrackedTransaction::commit returns it; acceptance still
requires a successful callback return. Transaction maintains no aggregate
correspondence. Its per-batch ApplicationState remains necessary for handle
resolution, but it performs no correspondence construction or updates.

S0e established independent tracking and measured its unconditional cost. The
2026-10-14 decision replaces that unconditional requirement with these separate
types. S0e1/S0e2 introduce the split and repeat the Journal comparison; S0f adds
and measures both recovery methods on both types.

Both implementations are retained. The selected Undo/Edit boxing and batching
work below improves their costs without choosing between them. Destructive
journal-free editing remains a separate contract: it provides no recovery
guarantee and is not a RollbackMethod variant.

## 3. Copy-producing APIs

Re-review operations that preserve their input and return independent results,
including reaction application and transform_iter. Determine whether the required
copy occurs once, and identify repeated copies across public methods, editors,
resolver phases, and Python conversions. Reaction application currently clones
the host and applies edits without a recovery journal.

Keep independent-result semantics distinct from receiver mutation. A cheaper
implementation must still let callers retain the source and the result.

## 4. Benchmarks before and after changes

Record a baseline before implementation and compare the changed APIs for:

- resolution and projection;
- each transformer;
- reaction application with independently retained products;
- sparse and extensive changes, with both uniquely owned and shared input storage;
- bounded failure cases that exercise recovery.

Measure time, allocation traffic, and peak extra memory. Compare equivalent
input/output and failure guarantees, including copy-on-write detachment,
publication checks, and recovery. Python measurements must include the actual
binding path. Use these results to choose the mechanism per workload rather than
assuming either journals or cloning wins universally.

### Recovery comparison: molecule size versus journal size

Vary molecule size and amount of mutation independently. A size series with the
same fraction of atoms changed would confound those two costs. Include ordinary
small molecules below 100 atoms and larger cases, without treating the larger
cases as typical.

| Dimension | Required comparisons |
| --- | --- |
| Molecule size | Hold the edit sequence's size and kind fixed while increasing unaffected molecular storage |
| Journal size | Hold the molecule fixed while increasing effective edits from one/few to extensive changes |
| Edit kind | Attribute changes; additions; removals with saved payloads and cascades; atom/ligand-list replacements |
| Storage touched | Concentrated changes versus changes spread across entity sets, to expose copy-on-write detachment costs |
| Input ownership | Unique storage versus storage shared with a retained molecule |
| Outcome | Successful commit; early and late application failure; final integrity failure; explicit rollback |
| Tracking | Plain and tracked execution, including correspondence construction |

Use a bounded set of these comparisons, not the full Cartesian product. Record
actual molecule/entity counts, effective edit count, Undo entry count, journal
capacity, and owned payload memory: equal entry counts need not mean equal
journal sizes. Separate allocated journal capacity from memory occupied by
entries and their payloads.

Compare complete clone/apply/install and journal/apply/commit paths on the same
inputs and edit sequences. Include final checks and dropping recovery storage;
on failure include discard or rollback. Check successful result agreement and
the promised receiver preservation before interpreting timings. Journal-free
consuming execution can show recovery overhead, but cannot win a comparison
requiring receiver preservation because it lacks that contract.

The result should identify where each mechanism is preferable and whether the
difference matters for resolution, reaction application, and transformations.
It must inform the recovery choice; it does not authorize automatic size-based
dispatch beyond the agreed RollbackMethod selector. Reaction measurements must retain
their source/product ownership requirements rather than timing a different
mutation contract.

### S0a baseline — 2026-10-14

**Completed before production changes.** Source 7fa2bd203b12c406b738a49c3df99db72d6365eb,
after completed 236; the only tracked Rust change is the stereo extension to
the maintained constraint_readings benchmark. The earlier dated studies below
remain historical evidence, not measurements of this revision.

Rust/Cargo 1.96.0, aarch64-apple-darwin, macOS 15.7.3; Python 3.13.16 from the
repository environment. Criterion uses the workspace bench profile (optimized,
debug information enabled); the recovery and Python harnesses use release defaults.
No compilation overlapped timing; measurements ran sequentially. Named Criterion
baseline: `181-s0a-20261014`, 30 samples, 250 ms warmup, 1 s measurement,
10,000 resamples, no plots. Tables report medians, not Criterion slope estimates.
Some short cases have wide confidence intervals: sparse8 unique apply is
0.104 µs with a 95% median interval recorded in the raw data; its interval width
is 26% of the estimate. This baseline alone cannot establish percent-level changes.

Coverage: 348 Criterion cases (editor 24, reaction 32, resolve 118,
constraint_readings 174); 41 recovery fixtures producing 298 timing rows and
298 memory rows, plus four reaction memory rows; 48 Python cases. All benchmark
fixture checks, recovery/result checks, and the three Python harness tests passed.
The Python harness calls actual production methods through PyO3 dispatch in a
fresh release test binary; the root extension was also freshly built and imported.
Only cfg(test) measurement modules differ in the isolated Python source copy.

The current measurements use these boundaries:

| Measurement | Preparation | Included execution and destruction |
| --- | --- | --- |
| Criterion editor/ownership | Input and prepared Edits outside the batched timer | Operation and publication; final returned/input storage destruction deferred by the existing batched harness |
| Criterion reaction/readings | Fixed operands prepared outside | Operation and returned-value destruction inside iter |
| Recovery | Fresh unique input or a clone sharing retained storage; Edits outside | Source copy when required, application, publication, restoration/discard and replaced-receiver destruction inside; final returned values outside |
| Python | Input owner, batch and model outside | Dispatch, conversion, operation, input-owner and output destruction inside |

Recovery timing retains 20 sample means of 32 executions with rotating route
order. Python uses 15 sample means of 32 prepared calls or 1,024 inverse calls,
after five warmups. Allocation measurements are separate runs: calls include
reallocations, bytes sum requested sizes, and peak is extra live requested Rust
heap above operation entry, excluding allocator overhead, CPython heap and RSS.
Stereo allocation rows average 1,000 calls. Reaction memory keeps the returned
products alive until counting stops; its Criterion time includes their destruction.

#### Editor, reaction and composite operations

Editor sparse8 is an eight-atom chain with one charge edit; dense80 is an
80-atom chain with eight charge edits, eight dative bonds and eight aromatic
systems. Time is µs, including editor creation/publication.

| Fixture / ownership | Direct | Apply | Transact |
| --- | --- | --- | --- |
| sparse8 / unique | 0.057 | 0.104 | 0.185 |
| sparse8 / shared | 0.268 | 0.276 | 0.345 |
| dense80 / unique | 0.723 | 0.896 | 1.269 |
| dense80 / shared | 2.382 | 2.479 | 2.900 |

The five-atom overlay-removal fixture also covers all six relation kinds in
stored/reordered frames. Median times range from 0.366–0.523 µs stored and
0.402–0.846 µs reordered; each kind remains a separate saved case.

Reaction application retains the six-atom/six-bond source and produces an
independent product, updating aromatic electron counts across a rotated frame.
Enumeration has one match. Memory entries are calls / bytes / peak bytes.

| Operation | Time, µs | Memory |
| --- | --- | --- |
| Supplied match | 3.279 | 61 / 7,898 / 3,232 |
| Tracked supplied match | 3.274 | 61 / 7,898 / 3,232 |
| Enumerated | 4.337 | 85 / 10,410 / 4,712 |
| Tracked enumerated | 4.350 | 85 / 11,690 / 4,712 |

The other maintained reaction cases include editor/application/tracking,
combine/split/meet_pushout and reversal. The resolver ownership cases below use
the actual in-place composite methods: sparse8 is an eight-carbon chain;
dense88 contains eight copies of [13CH3][C@H](F)/C=C/c1ccccc1. These use the
SMILES valence model and Natural isotope policy. Time is µs.

| Fixture / ownership | Resolve into | Project into |
| --- | --- | --- |
| sparse8 / unique | 10.585 | 0.998 |
| sparse8 / shared | 10.956 | 1.182 |
| dense88 / unique | 270.551 | 55.104 |
| dense88 / shared | 275.294 | 63.226 |

#### Recovery and transformer paths

The bounded fixtures retain the definitions in the historical ownership study
below: active-component / unaffected-chain-atom pairs (1,0), (1,64), (1,256),
(4,64), (16,0), both ownership modes, and the same controlled failure cases.
Resolve/project components have 11 atoms, benzene components six, and the
charge-delocalization component five. Counts and bonds are saved per row.

Consume has no receiver-preservation guarantee; Copy retains the source and
returns a result. Install is the existing harness's clone/apply/install reference,
not the future RollbackMethod::Clone implementation. InPlace calls current public
operations: composite Resolver and prepared transactions use journals; the three
transformers apply complete plans without journals. Transformers therefore have
no Install-versus-journal comparison. The following times are µs, unique input.

| Operation | Groups / extra atoms | Consume | Copy | Install | InPlace |
| --- | --- | --- | --- | --- | --- |
| Resolve | 1 / 64 | 45.714 | 48.386 | 49.545 | 50.160 |
| Resolve | 16 / 0 | 462.370 | 468.877 | 479.759 | 610.900 |
| Project | 1 / 64 | 9.779 | 12.344 | 13.361 | 15.443 |
| Project | 16 / 0 | 81.303 | 90.675 | 96.310 | 122.501 |
| Kekulize | 1 / 64 | 11.385 | 12.774 | — | 11.947 |
| Kekulize | 16 / 0 | 205.377 | 207.362 | — | 206.016 |
| Aromatize | 1 / 64 | 17.703 | 19.029 | — | 18.210 |
| Aromatize | 16 / 0 | 157.520 | 158.410 | — | 162.161 |
| Delocalize | 1 / 64 | 1.003 | 2.169 | — | 1.318 |
| Delocalize | 16 / 0 | 9.310 | 11.238 | — | 9.605 |

Memory for unique input with 16 active components and no extra chain:

| Operation | Consume: calls / bytes / peak | InPlace: calls / bytes / peak |
| --- | --- | --- |
| Resolve | 7,553 / 2,662,334 / 796,794 | 7,618 / 5,742,530 / 1,751,366 |
| Project | 612 / 351,130 / 83,374 | 690 / 1,127,962 / 472,006 |
| Kekulize | 1,697 / 262,606 / 24,358 | 1,713 / 264,006 / 24,358 |
| Aromatize | 3,050 / 571,220 / 404,076 | 3,066 / 572,620 / 404,076 |
| Delocalize | 86 / 35,270 / 17,730 | 102 / 36,670 / 19,070 |

Prepared charge batches vary molecule and mutation sizes independently.
Each effective edit produces one scalar Undo; no owned payload allocation is
needed for that Undo. Edit is 256 bytes and Undo 752 bytes. A separate Vec push
probe matching the journal's growth records 1/8/64 entries at capacities 4/8/64:
3,008/6,016/48,128 buffer bytes, of which 752/6,016/48,128 hold entries.
The following operation peaks also include application scratch storage.

| Atoms / edits / ownership | Install, µs | Journal, µs | Install peak, B | Journal peak, B |
| --- | --- | --- | --- | --- |
| 8 / 1 / unique | 0.438 | 0.271 | 1582 | 3014 |
| 96 / 1 / unique | 2.462 | 0.552 | 18478 | 3014 |
| 384 / 1 / unique | 12.175 | 1.452 | 73774 | 3014 |
| 96 / 8 / unique | 2.557 | 1.127 | 18478 | 6022 |
| 96 / 64 / unique | 3.305 | 5.202 | 18478 | 48134 |
| 384 / 64 / unique | 13.132 | 6.113 | 73774 | 48134 |
| 96 / 1 / shared | 1.815 | 1.863 | 18478 | 21486 |
| 96 / 64 / shared | 2.634 | 6.263 | 18478 | 66606 |

At 96 unique atoms, journal allocation traffic for 1/8/64 changes is
2/3/6 calls and 3,014/9,030/93,254 bytes; Install is three calls and 18,478 bytes
in each case. Tracked 1/64-edit times are 2.634/3.501 µs Install and
0.727/5.028 µs journal. Both tracked correspondence checks passed.

Failure medians, unique input, include recovery. The four edit cases use
64 prepared changes on 96 atoms; the integrity case adds a duplicate bond.

| Outcome | Install/discard, µs | Journal/recovery, µs |
| --- | --- | --- |
| Isotope underdetermination | 1.821 | 2.092 |
| Stereo contradiction | 1.426 | 1.518 |
| Unsupported bond charge | 0.311 | 0.296 |
| Early invalid handle | 0.587 | 0.590 |
| Late invalid handle | 3.294 | 6.708 |
| Final integrity failure | 6.569 | 15.260 |
| Explicit rollback | 3.286 | 5.427 |

These are bounded baseline cases. Chemistry-operation journal counts/capacities
and a complete edit-kind/payload matrix are not measured here; do not infer them
from allocation totals or the controlled scalar-edit cases.

#### Stereo readings and lattice operations

The maintained benchmark now includes tetrahedral stereo atoms and cis/trans
stereo bonds on six-atom fixtures with actual ligand incidence. Stereogenicity
is the bounded assertion key. All six reading methods and normalized_eq, meet,
join, matches and is_compatible are separate cases, with distinct equal form
operands prepared outside measurement. Internally contradictory forms measure
readings only. Common, rare singleton-set and semantic-change rows stay separate;
there is no assumed frequency-weighted aggregate.

Each cell below is atom / bond median in ns. L = known kind/literal coset,
U = known kind/undetermined coset, K = kindless configuration.

| Case | asserted | evidenced | normalized_eq | meet | join | matches | is_compatible |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Common: L, empty | 5.96 / 6.00 | 10.42 / 10.45 | 3.73 / 3.73 | 99.58 / 102.26 | 271.97 / 271.53 | 264.89 / 265.31 | 4.11 / 4.10 |
| Common: U, empty | 5.94 / 5.92 | 10.36 / 10.40 | 3.73 / 3.73 | 98.23 / 98.34 | 101.85 / 101.76 | 263.40 / 263.47 | 3.94 / 3.94 |
| Common: L, literal assertion | 13.70 / 13.60 | 20.24 / 20.13 | 5.59 / 5.61 | 295.33 / 298.90 | 473.29 / 475.95 | 187.14 / 186.61 | 8.41 / 8.39 |
| Common: K, empty | 5.96 / 5.93 | 10.36 / 10.47 | 2.80 / 2.80 | 16.21 / 16.21 | 15.87 / 15.86 | 147.92 / 149.53 | 2.80 / 2.80 |
| Changed semantics: K, literal assertion | 13.63 / 13.75 | 20.01 / 20.13 | 4.97 / 4.98 | 206.32 / 206.33 | 208.14 / 209.74 | 78.02 / 78.32 | 6.85 / 6.85 |
| Rare: L, singleton assertion | 106.24 / 106.67 | 111.89 / 111.36 | 13.38 / 13.37 | 292.26 / 291.15 | 472.38 / 472.26 | 184.49 / 184.36 | 11.88 / 11.82 |
| Changed semantics: K, bottom assertion | 66.69 / 66.16 | 67.65 / 66.69 | — | — | — | — | — |

Derived and derived_complete return Undetermined on these fixtures and take
2.44–2.48 ns without allocations. The complete asserted/evidenced variants have
the same expected values as their non-closed counterparts here. Allocation
counts/bytes below are identical for both entity families; normalized_eq and
is_compatible allocate nothing on the measured non-bottom forms.

| Case | Asserted/evidenced calls / B | Meet calls / B | Join calls / B | Matches calls / B |
| --- | --- | --- | --- | --- |
| Common: L, empty | 0 / 0 | 0 / 0 | 5 / 240 | 4 / 54 |
| Common: U, empty | 0 / 0 | 0 / 0 | 0 / 0 | 4 / 54 |
| Common: L, literal assertion | 0 / 0 | 7 / 131 | 12 / 371 | 2 / 48 |
| Common: K, empty | 0 / 0 | 0 / 0 | 0 / 0 | 4 / 54 |
| Changed semantics: K, literal assertion | 0 / 0 | 7 / 131 | 7 / 131 | 2 / 48 |
| Rare: L, singleton assertion | 4 / 75 | 7 / 131 | 12 / 371 | 2 / 48 |
| Changed semantics: K, bottom assertion | 2 / 27 | — | — | — |

The baseline still reads a stored kindless literal assertion as literal and a
stored bottom assertion as None. S2g1 changes those semantics to vacuity; their
timings must not be presented as same-result optimizations. No production
normalization, reading or lattice behavior changed in S0a.

#### Python boundary and inverse baselines

Resolved carbon chains have 8/80 atoms and 1/8 prepared charge edits. Edit
includes an unchanged finish; apply/tracked_apply consume, while the transaction
methods mutate. Unique/shared pairs below include final owner destruction, so a
shared case can be faster even when it detaches storage: the retained fixture
keeps other tables alive. These are current consuming-contract measurements;
they do not hide a required retained-source copy in preparation.

Times are µs; every paired cell is unique / shared. Allocation figures count
Rust traffic only.

| Atoms / method | Time | Calls | Bytes | Peak B |
| --- | --- | --- | --- | --- |
| 8 / edit | 0.749 / 0.345 | 0 / 0 | 0 / 0 | 0 / 0 |
| 8 / apply | 0.781 / 0.556 | 0 / 2 | 0 / 1576 | 0 / 1576 |
| 8 / tracked_apply | 0.993 / 0.736 | 2 / 4 | 120 / 1696 | 120 / 1696 |
| 8 / transact | 1.237 / 0.866 | 4 / 6 | 3400 / 4976 | 3360 / 4936 |
| 8 / tracked_transact | 1.404 / 1.006 | 6 / 8 | 3520 / 5096 | 3360 / 4936 |
| 80 / edit | 1.692 / 0.582 | 0 / 0 | 0 / 0 | 0 / 0 |
| 80 / apply | 1.833 / 2.273 | 0 / 2 | 0 / 15400 | 0 / 15400 |
| 80 / tracked_apply | 2.120 / 2.521 | 2 / 4 | 1272 / 16672 | 1272 / 16672 |
| 80 / transact | 2.693 / 2.897 | 5 / 7 | 9416 / 24816 | 6368 / 21768 |
| 80 / tracked_transact | 2.941 / 3.189 | 7 / 9 | 10688 / 26088 | 6368 / 21768 |

Resolution uses resolved chains with only lone-pair counts reset to Undetermined
during preparation; isotope composition and implicit hydrogens remain specified.
Both default-model routes must recover the exact source chain, matching successful
Rust resolution. The model argument is either omitted or a prepared Python conversion
of ChemistryModel::default. Input storage is shared; resolve consumes and
resolve_into mutates. Memory entries are calls / bytes / peak bytes.

| Atoms / model | Resolve, µs | Resolve into, µs | Resolve memory | Resolve into memory |
| --- | --- | --- | --- | --- |
| 8 / omitted_default | 9.112 | 9.358 | 71 / 12,956 / 8,028 | 73 / 21,980 / 14,044 |
| 8 / supplied_default | 247.087 | 250.769 | 2,883 / 919,642 / 906,692 | 2,885 / 928,666 / 906,692 |
| 80 / omitted_default | 82.223 | 86.695 | 662 / 167,668 / 96,884 | 668 / 357,172 / 193,140 |
| 80 / supplied_default | 325.928 | 327.970 | 3,474 / 1,074,354 / 906,692 | 3,480 / 1,263,858 / 906,692 |

Inverse covers both Add and Remove, direct entity and aggregate Delta methods.
Atom payloads contain element-variable strings; aromatic payloads contain atom
ids and electron counts. Sizes are cost probes, not a workload distribution.
The table reports Add medians in µs; Remove has the same memory counts and
separate retained timing rows. Exact inverse values are checked against Rust.

| Payload | Entity inverse | Aggregate inverse | Calls / bytes / peak B |
| --- | --- | --- | --- |
| atom / 16 | 0.278 | 0.348 | 4 / 144 / 144 |
| atom / 256 | 0.284 | 0.340 | 4 / 624 / 624 |
| atom / 4096 | 0.417 | 0.473 | 4 / 8,304 / 8,304 |
| aromatic / 6 | 0.263 | 0.337 | 4 / 144 / 144 |
| aromatic / 64 | 0.331 | 0.406 | 4 / 1,536 / 1,536 |

#### Reproduction and retained evidence

The maintained targets are built with cargo bench --offline --no-run:
graph-IR editor/reaction/constraint_readings and graph resolve. The recovery
harness is built twice with cargo build --offline --release, preserving the
non-instrumented executable before enabling its alloc feature. The isolated
Python harness uses cargo test --offline --release -p umol-py --lib --no-run
with the repository Python environment activated. The root extension is built
with maturin develop --offline --release. Build commands, logs, source hashes,
fixture definitions and exact executable invocations are retained in
scratch/181-s0a/README.md, manifest.json, build-*.log and runs.jsonl.

With builds complete and the Python environment active, run.py executes the
targets sequentially; summarize.py produces criterion.csv, recovery-time.csv,
recovery-alloc.csv and python.csv. readings-alloc.log holds all reading allocation
rows. criterion-baseline.tar.gz preserves every named Criterion baseline directory,
including raw samples and confidence intervals; the named baseline also remains
in the Cargo target directory. The maintained benchmark diff is saved as
benchmark.patch. Keep these inputs and raw data for S2g1/S4a/S12a comparisons.

S0a establishes the pre-change reference only. It introduces no recovery selector,
ownership API, boxing, batching or stereo semantic changes. The documented
figures preserve a useful baseline even if temporary artifacts are later removed;
the raw paired cases are needed for the planned implementation comparisons.

### Bounded ownership experiment — 2026-10-08

**Outcome:** molecule size, number of effective edits, and existing storage
sharing change which recovery mechanism wins. The results do not support
replacing all journaled transactions with cloning. They also show why the earlier
resolution result cannot be generalized to sparse mutations.

This experiment uses the existing consuming operations as the journal-free core,
with experimental copy-producing and clone-and-install adapters. Production APIs
and implementations are unchanged. The retained transaction interface can support
either recovery mechanism.

#### Method and coverage

- Source ec82fcdcc; workspace-locked dependency versions; Rust 1.96.0 release
  defaults on arm64 macOS. The isolated executable does not build Python.
- Four routes: consume without recovery; clone and return a separate result while
  retaining the source; clone and install on success; existing journaled _into
  execution. Controlled Edits use apply versus transact, including tracked pairs.
- Timing: median of 20 sample means, 32 executions per sample, rotating route
  order each sample. Input preparation and final result destruction are outside
  timing. Cloning, copy-on-write detachment, old-receiver destruction on install,
  integrity checks, journal disposal, and failure recovery are inside timing.
- Unique inputs are independently constructed; shared inputs retain another
  owner of the initial storage. No measurements run concurrently.
- A separate build instruments System allocation. Allocation traffic sums
  requested sizes, including the full new size of each reallocation. Peak extra
  heap is maximum live requested bytes above the execution-entry baseline,
  excluding allocator overhead and RSS. Input and prepared-Edits storage are
  already present at that baseline.
- All 41 fixture configurations passed result-agreement and receiver-preservation
  checks before measurement; both ownership modes and all four routes were
  measured (328 rows). Tracked edit cases also checked correspondence agreement.
  Recovery comparisons use normalized_eq, matching the journal contract.

Operation fixtures combine active components with an already resolved linear
carbon chain. Active components are raw/resolved
`[13CH3][C@H](F)/C=C/c1ccccc1` for resolution/projection, benzene with an aromatic
overlay for kekulization, concrete alternating-bond benzene for aromatization,
and cyclopentadienyl with one atom-localized negative charge for charge
delocalization. The latter has electron contributions [2,1,1,1,1] before
transformation and [1,1,1,1,1] with overlay charge -1 afterward. Active-component /
extra-chain-atom pairs are (1,0), (1,64), (1,256), (4,64), and (16,0). These are
controlled compositions, not a claim about typical molecular populations.

Controlled edit fixtures are resolved linear carbon chains. Each edit changes a
different atom's charge from 0 to 1, producing exactly one Undo with scalar
payload. Molecule-size / edit-count pairs are (8,1), (96,1), (384,1), (96,8),
(96,64), and (384,64); tracked cases use (96,1), (96,64), and (384,1).

Temporary harness and raw CSVs are in scratch/181-recovery. The definitions,
measurement conditions, and decision-relevant results are recorded here so that
removing scratch does not remove the findings.

#### Operation timings

All times below are microseconds with unique input. Consume lacks recovery;
Copy retains the source and returns an independent result. Only Install and
Journal are competing implementations of the same receiver-mutation contract.

| Operation | Active groups / extra atoms | Total atoms | Consume | Copy | Install | Journal |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Resolve | 1 / 64 | 75 | 62.004 | 64.197 | 66.106 | 66.562 |
| Resolve | 16 / 0 | 176 | 872.032 | 885.658 | 893.715 | 1004.086 |
| Project | 1 / 64 | 75 | 10.025 | 12.707 | 13.538 | 16.025 |
| Project | 16 / 0 | 176 | 82.044 | 91.979 | 99.800 | 114.119 |
| Kekulize | 1 / 64 | 70 | 13.969 | 15.382 | 15.710 | 14.658 |
| Kekulize | 16 / 0 | 96 | 206.859 | 209.193 | 209.755 | 210.958 |
| Aromatize | 1 / 64 | 70 | 15.227 | 16.969 | 17.183 | 15.921 |
| Aromatize | 16 / 0 | 96 | 157.549 | 159.360 | 159.024 | 172.754 |
| Delocalize | 1 / 64 | 69 | 1.016 | 2.173 | 2.658 | 1.303 |
| Delocalize | 16 / 0 | 80 | 9.340 | 11.246 | 12.229 | 13.188 |

Across the complete operation matrix, the following ratio is **Journal time /
Install time**, shown as unique / shared. Above 1 favors cloning; below 1 favors
journaling. Differences of only a few percent are inconclusive at this sampling
scale and are not a basis for selecting an implementation.

| Operation | 1 / 0 | 1 / 64 | 1 / 256 | 4 / 64 | 16 / 0 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Resolve | 1.08 / 1.11 | 1.01 / 1.06 | 0.96 / 1.04 | 1.05 / 1.09 | 1.12 / 1.11 |
| Project | 1.07 / 1.32 | 1.18 / 1.95 | 1.14 / 1.63 | 1.23 / 1.63 | 1.14 / 1.45 |
| Kekulize | 1.02 / 1.07 | 0.93 / 1.04 | 0.90 / 1.03 | 0.99 / 1.02 | 1.01 / 1.02 |
| Aromatize | 1.02 / 1.05 | 0.93 / 1.03 | 0.88 / 1.02 | 1.00 / 1.05 | 1.09 / 1.10 |
| Delocalize | 0.86 / 1.21 | 0.49 / 1.10 | 0.30 / 1.04 | 0.74 / 1.20 | 1.08 / 1.37 |

Resolution benefits most clearly from cloning when many components need changes;
with one active component and much unaffected storage the timing advantage
shrinks or reverses. Projection favors cloning throughout this matrix.
Kekulization is often close; journals avoid copying untouched storage in the
larger sparse unique cases. Charge delocalization shows the strongest sparse
case: at 261 atoms, Install takes 7.236 µs versus Journal's 2.149 µs. Copy takes
5.740 µs and Consume 1.775 µs, so internal owned execution remains materially
useful for this workload.

#### Molecule size and journal size

Time is in microseconds; peak extra heap is in bytes.

| Atoms | Effective edits / Undo entries | Input | Install time | Journal time | Install peak | Journal peak |
| ---: | ---: | --- | ---: | ---: | ---: | ---: |
| 8 | 1 | unique | 0.468 | 0.292 | 1,582 | 3,014 |
| 96 | 1 | unique | 2.494 | 0.580 | 18,478 | 3,014 |
| 384 | 1 | unique | 9.428 | 1.526 | 73,774 | 3,014 |
| 96 | 8 | unique | 2.606 | 1.148 | 18,478 | 6,022 |
| 96 | 64 | unique | 3.370 | 5.320 | 18,478 | 48,134 |
| 384 | 64 | unique | 10.354 | 6.183 | 73,774 | 48,134 |
| 96 | 1 | shared | 1.827 | 1.869 | 18,478 | 21,486 |
| 96 | 64 | shared | 2.707 | 6.413 | 18,478 | 66,606 |
| 384 | 64 | shared | 7.843 | 11.916 | 73,774 | 121,902 |

Undo is 752 bytes in this build even for a scalar charge change. For the unique
96-atom cases, clone/install makes 3 allocation/reallocation calls and requests
18,478 bytes regardless of whether 1, 8, or 64 atoms change. Journaling makes
2 / 3 / 6 calls and requests 3,014 / 9,030 / 93,254 bytes respectively. Its peak
is lower than total traffic because growing the journal frees earlier buffers.
These figures include operation scratch allocations, not only the journal.
The extra heap for one unique-input edit stays 3,014 bytes as the molecule grows
from 8 to 384 atoms; cloning must detach the touched atom table.

Tracked application preserves this crossover: for unique 96-atom input, one edit
costs 2.692 µs with clone/install versus 0.782 µs with journaling; 64 edits cost
3.601 versus 4.906 µs. Tracking does not make a universal recovery choice valid.
Shared input largely removes journaling's sparse advantage because both routes
must detach the shared table, while journaling additionally stores Undo.

#### Operation allocation costs

Unique input with 16 active components and no extra chain. Each entry is
allocation/reallocation calls / requested bytes / peak extra bytes.

| Operation | Clone/install | Journal |
| --- | ---: | ---: |
| Resolve | 13,163 / 3,774,938 / 830,814 | 13,108 / 6,801,790 / 1,751,366 |
| Project | 654 / 411,406 / 140,538 | 690 / 1,127,962 / 472,006 |
| Kekulize | 1,704 / 272,626 / 29,526 | 1,719 / 356,942 / 71,214 |
| Aromatize | 3,292 / 618,168 / 404,076 | 3,326 / 800,148 / 404,076 |
| Delocalize | 110 / 54,058 / 36,458 | 107 / 129,158 / 65,838 |

Resolution's dense peak falls from about 1.75 MB with journals to 0.83 MB with
cloning. Aromatization has the same peak in this case despite different allocation
traffic, showing that recovery storage is not always what determines the peak.

#### Unsuccessful execution

All recovering routes preserved the input under normalized_eq. Times are
microseconds for unique input. The controlled batch failures use a 96-atom chain
and 64 prepared charge changes.

| Outcome | Clone/discard | Journal/recovery |
| --- | ---: | ---: |
| Isotope underdetermination | 2.015 | 2.272 |
| Stereo contradiction | 1.516 | 1.619 |
| Unsupported bond charge | 0.320 | 0.298 |
| Invalid handle before changes | 0.633 | 0.613 |
| Invalid handle after 64 changes | 3.307 | 6.704 |
| Final integrity failure after changes and duplicate bond addition | 6.704 | 15.064 |
| Explicit cancellation after 64 changes | 3.369 | 6.038 |

The resolver rejection fixtures are respectively two methane components with an
unresolved isotope set on one, and methane with an unsatisfiable tetrahedral
stereo constraint. Projection rejects nonzero localized-bond charge. These are
bounded recovery checks, not an exhaustive sample of chemistry failures.

#### Decisions supported and remaining coverage

- Do not eliminate journals based on the resolution result: sparse changes to
  unique storage have demonstrated time and memory benefits from journaling.
- Source-preserving copy-producing methods and internal consuming execution
  remain useful independently of which mechanism implements _into.
- Clone-based resolution/projection merits a production candidate, but one
  universal mechanism for all mutating operations is not supported. No automatic
  size-based dispatch is designed or introduced by this experiment.
- Smaller Undo representations could move the crossover; the measurements
  describe today's journal, not every possible journal implementation.
- This completes the bounded Rust operation experiment. Reaction application,
  Python binding costs, a full edit-kind/removal-payload matrix, and actual journal
  counts/capacities for chemistry operations remain outside this run. Only the
  controlled charge-edit cases establish exact Undo counts. Those remaining
  measurements must not be inferred from these tables.

### Bounded journal study — boxing and batching

Two independent experiments will test whether a smaller journal changes the
recovery tradeoff. Start with boxing; neither experiment requires journals to
win for every operation or changes the public recovery guarantees.

1. **Boxing:** measure the largest Undo payloads, then compare selective boxing
   of those payloads. Keep ordinary field changes inline. Measure journal size,
   allocation count, peak memory, success time, and rollback time. Include removal
   workloads that pay for the new boxes, alongside the existing sparse and dense
   modification fixtures. Show the proposed production types before implementing
   them; the layout models below are experiments, not approved interfaces.
2. **Batching:** inspect the repeated edits emitted by resolution and projection.
   Compare typed lists for one field across entities and several fields across
   entities, with corresponding batched undo storage. Preserve old-value checks,
   application order, and recovery semantics. A Vec of the existing Edit or Undo
   enum does not remove its per-entry size cost. Define only batches justified by
   these consumers; do not expand the entire edit vocabulary mechanically.

Use the provisionally selected Undo 88 / Edit 104 layout for all three batching
configurations: individual edits, grouping by field/key, and grouping by entity.
Compare journal, clone/install, and no-recovery execution. Stop after these bounded candidates:
report which operations benefit and any regressions, rather than pursue a
universally winning mechanism or introduce adaptive dispatch.

#### Boxing: initial layout measurements — 2026-10-08

Measured on arm64 with Rust 1.96.0 using size_of, with compiler layout output
from Rust 1.98.0-nightly as a cross-check. Undo occupies 752 bytes. These sizes
describe inline storage, excluding allocations owned by vectors and forms.

| Variant payload | Bytes |
| --- | ---: |
| RestoreRemovedTopology | 752 |
| Each of the six RestoreRemoved overlay variants | 328 |
| RemoveAddedStereoAtom / RemoveAddedStereoBond | 168 |
| RemoveAddedAromaticSystem / RemoveAddedMulticenterBond | 128 |
| ModifyStereoAtomField / ModifyStereoBondField | 88 |
| RestoreStereoAtomConstraint / RestoreStereoBondConstraint | 80 |
| ModifyAtomField / ModifyBondField | 72 |

RestoreRemovedTopology contains two 256-byte values (MoleculeCompaction and
UndoCompaction), 144 bytes for RemovedOverlays (six vector headers), 48 bytes for
the removed atom and bond vector headers, and 48 bytes for CascadedConstraints
(two vector headers). UndoCompaction itself contains a MoleculeCompaction. Each
overlay-removal payload contains that 256-byte UndoCompaction, a 24-byte removed
entry vector, and the 48-byte CascadedConstraints.

Scratch enum models preserve every variant and field, changing only selected
payloads to one Box per selected variant:

| Layout | Bytes per Undo entry |
| --- | ---: |
| Current enum; unchanged model matches | 752 |
| Box only RestoreRemovedTopology | 336 |
| Box topology removal and all six overlay-removal payloads | 176 |

The last candidate reduces inline journal storage by about 77%, while adding
one allocation for each recorded removal. It adds no box allocation to ordinary
field changes. The next largest payloads are stereo additions. Removal frequency
in the measured chemistry operations still needs counting; these variants must
not be assumed rare merely because scalar modifications are common. No execution
speed improvement has yet been measured for these layouts.

The reproducible layout probe is scratch/181-recovery/examples/undo_layout.rs;
its size output is undo-payload-sizes.txt in the same experiment directory. The
measurements and conclusions above remain available after scratch is deleted.

#### Boxing: agreed experimental stages

Compare two stages before choosing the production representation:

1. **176-byte candidate:** box the complete RestoreRemovedTopology payload and
   the six RestoreRemoved overlay payloads, one allocation per undo entry.
2. **Below-100-byte candidate:** additionally box the payloads of
   RemoveAddedStereoAtom, RemoveAddedStereoBond, RemoveAddedAromaticSystem, and
   RemoveAddedMulticenterBond. Measure the resulting enum size; the target is not
   yet a verified layout result.

All payloads above 100 bytes belong to structural undos: either restoring removed
entities or removing added entities. Field-modification and constraint undos stay
inline in both candidates. The second stage adds allocations for those four
addition kinds, so its smaller journal entries do not by themselves establish
that it is preferable.

Measure both candidates against the unchanged journal and clone/install on the
same resolution, projection, and sparse-edit fixtures, including rollback. Cover
structural additions and removals to expose the allocation cost of the boxes.
Record variant counts, journal size, allocation count, peak memory, and execution
time. Select the stage from these results; neither stage is a commitment to make
journaling the universal recovery mechanism. Batching remains a separate part of
the study.

#### Boxing results — 2026-10-08

The full 3×3 boxing experiment is complete; the separate batching study follows
below.
Production code is unchanged. An isolated source copy at ec82fcdcc crosses all
three Undo sizes with all three Edit sizes:

- **Undo 752:** unchanged. **Undo 176:** box seven removal-undo payloads.
  **Undo 88:** additionally box four addition-undo payloads.
- **Edit 256:** unchanged. **Edit 184:** box AddMoleculeConstraint and
  RemoveMoleculeConstraint. **Edit 104:** additionally box AddStereoAtom,
  AddStereoBond, AddAromaticSystem, and AddMulticenterBond.

The two molecule-constraint variants have 248-byte payloads: ConstraintEdit
includes eight handle-vector headers occupying 192 bytes. Stereo-addition
payloads occupy 176 bytes, and aromatic-system and multicenter-bond additions
occupy 120 bytes. Boxing additions without the two molecule-constraint variants
leaves Edit at 256 bytes.

The experiment retains the earlier fixtures and adds 1/8 sequential topology
removals from a 96-atom chain, with commit and rollback, and 1/64 molecule-constraint
additions, with commit and rollback. All 49 fixture checks pass in every
configuration: successful results agree across routes, failures and cancellation
restore the receiver, and the tracked-edit cases agree on correspondence.
Undo variant counts and maximum journal capacities agree across all nine
configurations. No edit semantics, preconditions, or recovery policy changed.

Most uses have no molecule constraints. The four explicit molecule-constraint
fixtures are separate cost checks for an uncommon workload; primary conclusions
use resolution, projection, transformations, and ordinary edits. There is no
aggregate score giving all fixtures equal weight. Edit 184 benefits ordinary
plans by shrinking their entries without adding allocations to those plans.

Timing uses Rust 1.96.0 on arm64, release builds, 20 sample means of 32 executions,
rotating route order. Builds and timing runs do not overlap. Whole-operation
timings and allocation measurements include edit construction and disposal;
application-only allocations are recorded separately. Timing builds contain no
allocation or journal-count instrumentation. Input preparation and destruction
of returned outputs remain outside the measured region. Both unique and shared
inputs are measured; the matrices below show unique input. These are controlled
fixtures, not a molecular population survey.

Selected journaled times in microseconds. Dense denotes 16 active groups without
ballast; the ordinary edit cases use a 96-atom chain.

| Undo bytes | Edit bytes | Dense resolution | Dense projection | Dense aromatization | 64 charge edits | Eight topology removals |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 752 | 256 | 993.56 | 120.01 | 177.06 | 7.33 | 76.98 |
| 752 | 184 | 965.11 | 117.71 | 171.70 | 7.19 | 77.04 |
| 752 | 104 | 938.01 | 107.95 | 171.89 | 7.00 | 75.05 |
| 176 | 256 | 873.54 | 91.04 | 161.43 | 4.62 | 76.27 |
| 176 | 184 | 874.74 | 90.01 | 161.59 | 4.50 | 76.89 |
| 176 | 104 | 864.12 | 88.51 | 159.23 | 4.24 | 75.71 |
| 88 | 256 | 866.74 | 92.14 | 161.98 | 4.48 | 74.87 |
| 88 | 184 | 864.02 | 89.03 | 162.86 | 4.22 | 77.80 |
| 88 | 104 | 848.98 | 88.09 | 159.14 | 4.02 | 75.08 |

Peak extra live requested memory in KiB, including edit construction:

| Undo bytes | Edit bytes | Dense resolution | Dense projection | 64 charge edits | Eight topology removals |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 752 | 256 | 1710.32 | 460.94 | 63.01 | 37.20 |
| 752 | 184 | 1674.32 | 442.94 | 58.51 | 36.63 |
| 752 | 104 | 1634.32 | 422.94 | 53.51 | 36.01 |
| 176 | 256 | 867.24 | 173.90 | 27.01 | 37.84 |
| 176 | 184 | 867.24 | 155.90 | 22.51 | 37.28 |
| 176 | 104 | 867.24 | 135.90 | 17.51 | 36.65 |
| 88 | 256 | 825.24 | 129.90 | 21.51 | 37.15 |
| 88 | 184 | 825.24 | 111.90 | 17.01 | 36.59 |
| 88 | 104 | 825.24 | 91.90 | 12.01 | 35.96 |

Dense resolution records 1,168 undos: 848 field modifications, 288 entity-constraint
restorations, and 32 additions (16 aromatic systems and 16 stereo bonds). Its peak
journal capacity is 2,048 entries. Undo 88 adds 32 allocations; Edit 104 adds
another 32. Edit 184 adds none to this fixture. Dense projection records 387 undos
with capacity 512, including three removal undos. Undo 176/88 add three
allocations there; neither boxed Edit layout adds any. The eight-topology-removal
fixture gains eight allocations and sees little timing benefit from Undo boxing.

The rare molecule-constraint case pays for its boxes: 64 additions gain 64
allocations with either Edit 184 or 104. With Undo 88, times for Edit 256/184/104
are respectively 15.18/15.55/15.70 microseconds. This is an explicit cost of the
representation, not the basis for selecting the common-workload layout.

Provisional choice: **Undo 88 bytes / Edit 104 bytes**, carried forward into the
batching study. The drop to Undo 176 delivers most of the speed gain; further
memory reduction remains valuable where timing changes little. Compared with
176/256, 88/104 reduces measured peak extra memory from 173.90 to 91.90 KiB for
dense projection and from 27.01 to 12.01 KiB for 64 charge edits. The additional
boxing allocations are an accepted tradeoff for this provisional choice. This
selects a representation candidate, not a universal recovery mechanism.

Findings supporting the choice:

- Undo boxing is worthwhile. The 176-byte stage captures most of the speed gain;
  88 bytes reduces memory further. Additional boxing does not make every operation
  faster: topology-removal time is largely unchanged and small timing differences
  between compact layouts do not establish a universal winner.
- Edit 184 offers a reduction without new allocations in the ordinary fixtures.
  Edit 104 reduces plan memory further and is useful for the measured ordinary
  workloads. With Undo 88, dense projection falls from 111.90 to 91.90 KiB and the
  64-charge-edit case from 17.01 to 12.01 KiB. Its extra allocations for entity
  additions are visible in the resolution counts above.
- Neither Edit stage is an unconditional speed improvement. For example, with
  Undo 752, projection with 64 unaffected atoms measures 16.42/16.22/21.01
  microseconds for Edit 256/184/104. With Undo 176 those values are
  12.72/12.16/11.91, and with Undo 88 they are 13.64/11.63/11.45. The interaction
  between the two representations matters; choosing them solely from separate
  measurements would miss this.
- The earlier cloning advantage is substantially reduced. With Undo 88/Edit 104,
  unique-input dense resolution is 848.98 microseconds journaled versus 839.95
  with clone/install, and dense projection is 88.09 versus 89.54. Differences
  around 1–2% do not establish a decisive winner. Unique-storage sparse changes
  still favor journals, while repeated topology removal favors clone/install on
  time and journals on peak memory.
- Shared storage remains relevant: with Undo 88/Edit 104, 64 charge edits on 96
  shared atoms take 5.04 microseconds journaled versus 4.25 with clone/install;
  dense projection takes 90.44 versus 83.15. Boxing does not establish a universal
  recovery mechanism.
- Payload-size reduction is demonstrated; reducing journal entry count remains
  the independent batching question. No batching or public lifecycle redesign was
  folded into this experiment.

Reproduction scripts, isolated source, binaries, and complete CSVs are under
scratch/181-boxing. The README records the commands and measurement boundaries.
The results above are the durable record; scratch can be deleted without losing
the conclusions. Production payload type names and interfaces remain to be settled
before implementation.

#### Batching: resolution experiment — 2026-10-08

**Result:** grouping by field/key reduces application cost and journal storage,
but does not establish a material speedup for complete resolution. Grouping by
entity with separate per-entity vectors is worse on time and allocation count.
Neither result establishes a universal recovery mechanism.

All three configurations use the provisional Undo 88 / Edit 104 layout:

- **Individual:** existing single-field and single-constraint edits and undos.
- **By field:** typed vectors of `(handle, old, new)` for each atom/bond field or
  entity-constraint key; undos use concrete ids, inverse field values, or the
  previous optional constraint value.
- **By entity:** rows containing one handle, a vector of field changes, and a
  vector of constraint changes; corresponding undo rows use concrete ids.

Only atom and bond updates emitted by resolver plans are batched. They cover all
1,136 repeated field/constraint modifications in the dense fixture. The 32
aromatic-system/stereo-bond additions remain ordinary edits. Other operations are
ordering barriers. Changes to the same field or constraint key retain their
order; grouping can reorder independent changes. Arbitrary Edits are not
regrouped, and no independent-batch composition API is introduced.

Producers append directly to typed vectors; no completed Edit list is repacked.
Entity rows are constructed once per update call. Edit vectors grow normally;
undo vectors reserve their required lengths. An initial collector that dispatched
through a single-field Edit for every value made bond planning take 18.37 µs.
Direct typed emission reduced it to 10.71 µs; the results below use that version.
This construction cost matters when assessing the representation.

**Verification and measurement.** The scratch source copy uses ec82fcdcc,
Rust 1.96.0, arm64, and release builds. Each timing is the median of 20 sample
means of 32 calls. Timings and allocation/counter measurements use separate
executables; builds and timed runs do not overlap. Each batched executable also
executes the unbatched plans and compares complete outputs exactly. All five
successful resolution fixtures and the underdetermined/contradictory fixtures
pass. Eight-entry field and constraint batches verify success, first/middle/last
old-value mismatch, and explicit rollback. Failed batched application restores
its applied prefix before the transaction restores earlier journal entries.
No allocation-failure or internal-panic recovery machinery is added.

Full resolution includes planning, edit construction/application, integrity checks,
and journal disposal. Input preparation and destruction of the returned result
are excluded, as in the boxing study. Allocation figures are requested bytes,
not RSS. The fixtures use the existing 11-atom active group and resolved carbon
ballast; they contain no molecule-level constraints. Their entity-level constraints
are the actual resolver workload, not extra synthetic assertions.

**Complete journaled resolution, unique input.** Each cell is time in µs /
peak extra KiB. Compare columns within this run, not timings across sessions.

| Active groups / ballast atoms | Individual | By field | By entity |
| --- | ---: | ---: | ---: |
| 1 / 0 | 38.52 / 18.61 | 39.06 / 10.59 | 41.05 / 16.26 |
| 1 / 64 | 65.00 / 28.82 | 64.50 / 27.97 | 65.69 / 28.95 |
| 1 / 256 | 143.93 / 81.16 | 141.59 / 80.30 | 145.05 / 81.29 |
| 4 / 64 | 190.44 / 142.40 | 188.44 / 138.31 | 202.92 / 142.25 |
| 16 / 0 | 857.76 / 825.24 | 858.26 / 808.87 | 910.00 / 824.62 |

**Dense resolution, all recovery mechanisms.** Times in µs; each cell is
no recovery / clone-and-install / journal. Both unique and shared inputs are
measured. Copy-producing execution is also retained in the complete CSVs.

| Representation | Unique input | Shared input |
| --- | ---: | ---: |
| Individual | 844.17 / 855.00 / 857.76 | 849.08 / 848.40 / 860.74 |
| By field | 854.43 / 863.45 / 858.26 | 863.26 / 866.22 / 872.51 |
| By entity | 894.32 / 891.78 / 910.00 | 903.43 / 897.47 / 934.06 |

Field batching leaves these complete-operation timings close to the individual
representation. It does not make journal recovery decisively faster than cloning.
The entity representation adds roughly 1–10% to journaled resolution time across
the measured fixtures, with no whole-operation memory advantage over field grouping.

**Dense journal storage and allocation.** All configurations perform the same
848 field changes, 288 entity-constraint changes, and 32 additions.

| Quantity | Individual | By field | By entity |
| --- | ---: | ---: | ---: |
| Outer Undo records | 1,168 | 40 | 36 |
| Peak journal capacity, records | 2,048 | 64 | 64 |
| Outer journal allocation, KiB | 176 | 5.5 | 5.5 |
| Batched constituent changes | — | 1,136 | 1,136 |
| Entity rows in batched storage | — | — | 640 |
| Edit constituent vector length / summed capacity | — | 1,136 / 1,680 | 1,136 / 2,560 |
| Undo constituent vector length / summed capacity | — | 1,136 / 1,136 | 1,136 / 1,136 |
| Total allocation calls | 13,172 | 13,216 | 14,464 |
| Total requested allocation, KiB | 3,530.34 | 3,153.75 | 3,255.97 |
| Peak extra memory, KiB | 825.24 | 808.87 | 824.62 |

The 5.5 KiB figure excludes typed vectors and their payloads. Reducing outer
records does not eliminate the individual changes. Entity row vectors also have
summed edit/undo capacities of 896/640 rows, in addition to their constituent
vectors. Field batching cuts total requested bytes by about 11%, but whole-operation
peak memory by only 2% in the dense fixture. The 11-atom fixture shows a larger
peak reduction, 18.61 to 10.59 KiB. These are different memory questions.

**Separate phase diagnostics.** The existing isotope and bond resolvers supply
the plans. The isotope diagnostic explicitly makes the active atoms' isotope
mass Undetermined: the ordinary raised resolution fixtures already have determined
isotopes and would otherwise produce an empty isotope plan. There are 176 isotope
changes or 352 bond-field changes. Plan construction and prepared journaled
application are timed separately; application inputs here have unique storage.
Each cell is µs / allocation calls / peak extra KiB.

| Phase | Measurement | Individual | By field | By entity |
| --- | --- | ---: | ---: | ---: |
| Isotope | Plan | 4.75 / 8 / 26.01 | 4.51 / 10 / 12.41 | 7.14 / 186 / 60.41 |
| Isotope | Apply | 5.58 / 9 / 22.01 | 3.06 / 4 / 6.88 | 9.43 / 180 / 9.69 |
| Bonds | Plan | 10.37 / 9 / 52.01 | 10.71 / 18 / 32.41 | 18.52 / 186 / 60.41 |
| Bonds | Apply | 9.97 / 10 / 44.01 | 5.68 / 5 / 7.60 | 13.11 / 180 / 9.76 |

Prepared journaled application improves by 43–45% with field batching; direct
application also improves (isotope 3.28 to 2.52 µs; bonds 6.07 to 4.34 µs).
Planning must be counted alongside that saving. Field-batch planning is close to
the ordinary representation; the entity vectors make both planning and journaled
application substantially more expensive. The phase timings are diagnostics,
not an additive breakdown of complete resolution.

**Resolution disposition:** grouping by field/key proceeds to the projection
comparison below. Do not pursue the tested per-entity vector representation.
Resolution alone supports lower application cost and storage, not an end-to-end
speed claim or a decision to use journals everywhere.
No production API or implementation is changed. The harness, generators, complete
CSVs, and initial-collector comparison are in scratch/181-batching; these findings
are the durable record.

#### Batching: projection experiment — 2026-10-08

**Result:** field/key batching improves complete journaled projection on these
fixtures. Unique-input time falls by about 2–17%, and peak extra memory by 26–61%.
Shared-input time falls by 4–13%, with peak reductions of 16–35%. Cloning remains
faster for shared inputs; batching does not establish a universal recovery choice.
Only individual versus field-grouped edits are tested, both at Undo 88 / Edit 104.
The per-entity vector candidate is excluded.

The source copy, direct typed emission, measurement boundaries, five fixture shapes,
and timing/allocation instrumentation are those of the resolution experiment.
All projection phases are enabled. Unlike resolution, projection changes the
additional resolved carbon atoms: their natural-isotope defaults are elided.
For example, one active group plus 256 carbon atoms produces 266 isotope edits.
These are not untouched ballast for this operation. No molecule-level constraints
are added to the fixtures.

**Complete journaled projection.** Each cell is µs / peak extra KiB.

| Groups / additional carbon atoms | Unique, individual | Unique, by field | Shared, individual | Shared, by field |
| --- | ---: | ---: | ---: | ---: |
| 1 / 0 | 6.50 / 6.83 | 6.35 / 5.06 | 6.55 / 11.14 | 6.29 / 9.37 |
| 1 / 64 | 11.94 / 26.28 | 11.05 / 12.81 | 13.95 / 48.00 | 12.67 / 34.53 |
| 1 / 256 | 28.23 / 98.28 | 23.53 / 38.31 | 35.34 / 172.50 | 30.64 / 112.53 |
| 4 / 64 | 26.50 / 41.20 | 24.88 / 19.54 | 29.13 / 73.84 | 27.29 / 52.18 |
| 16 / 0 | 87.32 / 91.90 | 82.54 / 50.00 | 91.79 / 150.77 | 86.94 / 108.86 |

**Dense projection by recovery mechanism.** Times in µs, for 176 atoms.

| Recovery | Unique, individual | Unique, by field | Shared, individual | Shared, by field |
| --- | ---: | ---: | ---: | ---: |
| None | 79.45 | 79.16 | 84.82 | 84.29 |
| Clone and install | 89.52 | 89.33 | 85.06 | 84.28 |
| Journal | 87.32 | 82.54 | 91.79 | 86.94 |

With field batching, journals are faster than clone/install for unique inputs in
all five fixtures, while clone/install is faster for shared inputs in all five.
The largest unique fixture with additional carbon atoms illustrates the difference:
no recovery / clone-install / journal take 21.57 / 33.66 / 23.53 µs. For shared
input, the corresponding times are 29.67 / 29.27 / 30.64 µs. Copy-producing
measurements are retained in the CSVs; no copying API change is inferred here.

**Dense journal storage and allocation, unique input.** The operation performs
160 field changes, 224 entity-constraint changes, and three removal operations.

| Quantity | Individual | By field |
| --- | ---: | ---: |
| Outer Undo records | 387 | 8 |
| Peak journal capacity, records | 512 | 8 |
| Outer journal allocation, KiB | 44 | 0.69 |
| Batched constituent changes | — | 384 |
| Edit constituent vector length / summed capacity | — | 384 / 544 |
| Undo constituent vector length / summed capacity | — | 384 / 384 |
| Allocation calls | 693 | 706 |
| Total requested allocation, KiB | 271.86 | 153.60 |
| Peak extra memory, KiB | 91.90 | 50.00 |

The eight records are five field/key batches plus the three unchanged removal
undos. Typed vectors and removal payloads remain allocated; 0.69 KiB is only the
outer journal. Field batching makes slightly more allocation calls while reducing
their total requested bytes by 44% and peak extra memory by 46% in this fixture.

**Failure checks.** Both representations pass the five successful fixtures,
whole-result equivalence, the existing early bond-charge rejection, and a late
isotope rejection after stereo and aromatic projection have run. A separate check
asserts IsotopeProjectError::NonGroundIsotope for that late failure and verifies
source restoration with batching disabled and enabled. The existing batched
field/constraint prefix-failure checks also pass. On the 11-atom late-failure
fixture, unique-input journal time is 6.88 / 6.87 µs and peak extra memory is
5.96 / 3.56 KiB for individual / field batching; rollback correctness is preserved
without claiming a failure-path speedup.

**Disposition:** field/key batching is supported by the combined resolution and
projection evidence: it reduces application/storage costs, keeps full resolution
roughly unchanged, and improves complete journaled projection. This supports
carrying that representation into production design, subject to settling its
interfaces; it does not select journals over cloning for every operation or
ownership condition. Reaction workloads remain unmeasured in this batching study.
Production sources are unchanged. Reproduction uses projection.py and
projection-run.py in scratch/181-batching; projection-prefixed CSVs and logs retain
the full measurements. These tables preserve the findings independently of scratch.

## 5. Python consumed-object errors

ConsumedError currently applies to Molecule, MoleculeEditor, Edits, and RingViews.

Python Molecule will hold GraphIrMolecule directly, without Option. Migrate its
current consuming methods—edit, apply, tracked_apply, resolve, and
resolve_with_report—to the proposed source-preserving operations. Remove
Molecule::take, molecule-consumed errors, and the fallibility of internal accessors
whose only failure is a consumed molecule. No empty placeholder is installed.
The destructive Rust editor entry point is not required on the Python surface.

Edits retains consumption on transfer and gains an explicit copy method.
MoleculeEditor retains consumption at finish; application mutably borrows it.
RingViews retains consumption through into_ring_set, including ConsumedError for
the consumed owner and InvalidatedViewError for its dependent accessors.

### Remaining consumed wrappers — review, 2026-10-08

Source inspection at ec82fcdcc found only MoleculeEditor and Edits in addition
to Molecule using consumed-owner state. Independent review confirmed the inventory
and the ownership distinctions below. No runtime measurements were made in this pass.

- **MoleculeEditor:** the current Python apply/tracked_apply consume the editor
  and return a new one. Follow the proposed Rust mutable-borrow methods instead:
  Python apply returns None, tracked_apply returns correspondence, and both retain
  the editor, including partial state after application failure. finish still
  consumes it on success or integrity failure. Keeping Option for this terminal
  operation avoids copying the completed working molecule. Agreed disposition:
  retain consumption only at finish.
- **Edits:** application and transaction methods transfer the batch; receiver
  rollback does not restore it. Rust Edits owns Vec<Edit> and creation counters.
  Preserving it requires copying the sequence and owned payloads; shared form
  storage can remain shared. Molecule cloning measurements do not establish this
  cost. Agreed disposition: retain transfer, with an explicit copy method for
  reuse rather than unconditional copying. Python currently lacks that method;
  reconstructing Edits from its iterator converts entries to Python and back.
  copy delegates to Rust Edits::clone without per-entry Python conversion.
  Batch-copy measurements are recorded in Study 3 below.

Evidence: [editor wrapper](../umol-py/src/transaction.rs),
[batch wrapper and iterator](../umol-py/src/edit.rs),
[Rust Edits](../umol-graph-ir/src/ir/edit.rs), and
[edit application](../umol-graph-ir/src/ir/molecule/apply.rs).
This inventory covers consumed-owner wrappers, not every binding that copies a
Rust value before invoking a consuming operation.

## 6. View invalidation

The existing counter mechanism remains: receiver mutation can invalidate views
even when Molecule itself cannot be consumed. Remove consumed-molecule checks
from view access; retain counter checks. Keep the current invalidation semantics:

- Producing an independent result leaves source views usable.
- Whole-molecule mutating operations advance the counter immediately before
  execution, invalidating all existing views and iterators, including nested
  accessors. Invalidation remains after rollback or a no-op, with the same rule
  for RollbackMethod::Clone and RollbackMethod::Journal.
- Argument or batch extraction failure before execution leaves views valid.
- Ordinary attribute and entity-constraint assignment preserves views; reads
  observe the updated values.

Installing a cloned result still changes the receiver's entity-id spaces, so
clone-based execution does not eliminate the need for invalidation. Copy-producing
operations leave the original owner unchanged. No stable-identity or view-repair
mechanism is introduced. Edits iterators remain invalidated by batch consumption.

## 7. Python interfaces and guarantees

For each affected operation, specify receiver mutation, returned value, whether
inputs remain reusable, failure behavior, and effects on aliases and views.
Justify each difference from Rust, including ownership and copying. Update the
bindings, exceptions, examples, tests, and Python API guide together after these
contracts are settled. The broader Python parity and type-role audit remains in
[231](231-python-api-parity-2026-09-21.md).

## 8. Reconsider the premise about Python cloning

The resolution findings do not support a general presumption against cloning.
A copy that provides useful ownership and failure semantics can be preferable to
transferring Rust's consuming lifecycle into Python. Reconsider the Python API
guide's preference for preserving Rust consumption on that basis.

This does not justify redundant copies or detached mutable property results.
Cheap Molecule cloning also does not establish that arbitrary forms or model
tables are cheap to copy. In particular, Python ValenceCandidateSource::to_rust
clones its registry or table into Cow::Owned; repeated conversion of a supplied
ChemistryModel remains a separate concrete copy-cost investigation.

### Study 1 — ChemistryModel and valence data

**Question:** eliminate repeated registry/table copying while preserving the
immutable Python model and its operation semantics. Scope: ChemistryModel,
ValenceModel, ValenceCandidateSource, AtomTypeRegistry, and ValenceTable.

Current Python argument extraction clones the field-based ChemistryModel wrapper;
to_rust then builds another Rust model, cloning the candidate registry/table into
Cow::Owned. Borrowing a custom table into that newly assembled model is obstructed
by its static Cow lifetime. The alternative is to retain a Rust ChemistryModel
in the Python wrapper and borrow the complete model for the operation. This does
not require changing Rust model lifetimes. A Python-reference extraction of the
existing wrapper, followed by its current conversion, provides a smaller
comparison that removes only the extraction copy.

Measure wrapper construction separately from repeated extraction/conversion and
complete resolution. Include the omitted-model default, explicitly supplied
built-in model, and custom registries/tables at stated sizes. Assert model and
resolution-result equivalence. Record time, Rust allocation calls/requested bytes,
and peak extra Rust allocation; the Rust allocator does not measure the Python
heap. Model getter behavior and construction copying remain part of the review.

#### Results — 2026-10-08

**Decision supported:** storing the complete Rust ChemistryModel and borrowing it
for calls removes substantial recurring work for explicitly supplied models.
Python model objects are frozen and their getters return values, so this is
compatible with current semantics. The experiment covers storage and borrowing, extended below to Python
construction and getters; production migration and repr adaptation remain pending. No Arc or Rust
model-lifetime change is required. Merely borrowing the current Python wrapper
removes the extraction clone but still pays for reconstruction.

Boundary costs per call, including destruction of temporary converted models:

| Supplied model | Current extract + convert, µs | Borrow wrapper + convert, µs | Current Rust allocations / requested KiB | Borrow + convert allocations / KiB |
| --- | ---: | ---: | ---: | ---: |
| Built-in default | 260.03 | 121.98 | 2813 / 885.44 | 1406 / 442.54 |
| Registry, 8 entries | 2.93 | 1.41 | 43 / 14.99 | 21 / 7.31 |
| Registry, 64 entries | 22.99 | 10.23 | 279 / 112.99 | 139 / 56.31 |
| Registry, 256 entries | 109.41 | 50.30 | 1109 / 452.22 | 554 / 225.93 |
| Table, 8 elements | 1.04 | 0.49 | 53 / 2.05 | 26 / 0.84 |
| Table, 64 elements | 7.76 | 3.98 | 407 / 17.13 | 203 / 8.38 |
| Table, 118 elements | 14.55 | 7.26 | 751 / 33.97 | 375 / 16.80 |

Borrowing the stored Rust model allocated nothing in all seven cases; the measured
bare borrow was about 0.004 µs and excludes outer Python method dispatch. The
explicit default costs more than the synthetic 256-entry case. Its current
conversion also loses the Rust default's borrowed static registry, making it
owned. The custom cases use owned registries/tables; their savings do not rely on
static storage. Synthetic registry entries repeat built-in patterns; table entries
use actual element keys and synthetic valences. These test allocation scaling,
not representative chemistry.

Complete one-atom resolution with the default model (same input clone in each
route, exact equal determined results):

| Route | µs | Rust allocations | Requested KiB | Peak extra KiB |
| --- | ---: | ---: | ---: | ---: |
| Explicit model, current conversion | 249.78 | 2827 | 890.69 | 885.44 |
| Explicit model, borrow wrapper + convert | 122.68 | 1420 | 447.79 | 446.31 |
| Explicit stored Rust model, borrow | 2.01 | 14 | 5.25 | 3.77 |
| Model argument omitted | 2.00 | 15 | 5.25 | 3.78 |

This is a supplied-model boundary improvement, not a general resolution speedup:
the omitted-model path already avoids the copying. These measurements use one
atom and the default model; the construction-inclusive custom-model comparison
follows below. Large-molecule workloads were not benchmarked. Converting a borrowed default Rust model into the
current wrapper cost 125.22 µs / 1,408 Rust allocations; cloning it into the stored
Rust prototype cost 0.032 µs / one 6-byte allocation. These are Rust-to-wrapper
construction measurements, not Python constructor timings; an already-owned Rust
model can move into that wrapper.

#### Construction, getters, and single-use check — 2026-10-08

The stored-model prototype now has a Python constructor, default constructor, and
all four component getters. Custom construction borrows the supplied frozen
component wrappers and converts them once; it does not retain Python borrows.
The current constructor instead extracts owned copies of those wrappers. Both
resolution routes are invoked through Python functions with their actual argument
extraction types, and return Python molecule objects. They use the same input
molecule clone and the same Rust resolver. All models, getter results, and resolved
molecules were exactly equal; each resolution was determined.

The measured workflow constructs the registry/table, ValenceModel, and enclosing
ChemistryModel, then resolves once or ten times. Inputs are existing Python atom
forms or element/ValenceEntry mappings; preparing/parsing those inputs and the
unchanged connectivity/aromaticity/stereo components is outside the measurements
for both routes. The custom registry contains the 35 built-in carbon patterns.
Tables cover 8/64/118 element keys, with carbon target four and synthetic unused
entries. These are construction-cost probes, not representative workload claims.

Times in µs, including Python dispatch and temporary/result destruction:

| Model | Construct current / stored | Construct + 1 resolve current / stored | Construct + 10 resolves current / stored |
| --- | ---: | ---: | ---: |
| Explicit default | 120.95 / 0.18 | 343.38 / 2.68 | 2468.10 / 23.86 |
| Custom registry, 35 patterns | 1207.41 / 1213.86 | 1223.96 / 1205.25 | 1345.55 / 1225.83 |
| Custom table, 8 elements | 7.48 / 7.18 | 10.55 / 9.35 | 36.34 / 25.09 |
| Custom table, 64 elements | 219.48 / 219.48 | 230.18 / 220.86 | 320.37 / 239.72 |
| Custom table, 118 elements | 682.00 / 704.27 | 707.94 / 710.91 | 880.54 / 728.80 |

**Decision supported:** there is no large additional model-construction cost
concealed by the borrowing comparison. Creating the explicit default becomes
much cheaper because it keeps the existing Rust static registry. For custom
models, registry/table construction dominates one use and remains largely
unchanged; the 118-element one-use candidate was slightly slower in this run.
Do not claim a universal one-use speedup. Ten uses benefit from avoiding repeated
conversion. Omitting the model remained unchanged, 2.31/2.29 µs per resolution.

There are small costs to retain explicitly:

- Component getters still produce independent Python values. Fetching all four
  default components took 115.83/111.21 µs (current/stored); for the custom registry,
  6.73/7.39 µs. The prototype has one extra Rust allocation in the combined getter
  path, about 160–166 extra requested bytes. Model storage does not make extracting
  a large component free, and repeatedly rebuilding a model from its getters is
  not the reuse path measured above.
- tracemalloc reported peak Python-heap use of 632/808 bytes for construction plus
  one resolution and 760/936 bytes for ten, current/stored, across all five cases:
  176 extra bytes for the stored-model workflow. These figures exclude Rust heap
  allocations and are not RSS.
- Rust allocations for explicit-default construction plus one resolution fell
  from 4,236 calls / 1,333.76 KiB requested to 15 calls / 5.25 KiB. For the custom
  118-element table the corresponding figures were 3,185 / 392.62 KiB and
  2,433 / 358.75 KiB. Custom construction still allocates substantial storage;
  it is paid once rather than reassembled at each resolver call.

One additional focused rstest case passed in both timing and allocation modes.
Timing used the same release build, five warmups and fifteen sample means, with
32 repetitions for construction/getters/reuse and 16 for construct-and-use.
Python peak tracing ran separately with garbage collection before measurement;
tracing was disabled during timing. Source and reproduction:
scratch/181-python-studies/model/single_use.rs, run-single-use.py,
build-single-use.log, and model-single-use-{time,alloc}.log. The combined summary
now contains 276 measurement rows. Production sources remain unchanged.

### Study 2 — Reaction reconstruction and application

**Question:** identify avoidable copying/rechecking in repeated use of one Python
Reaction. Current to_rust clones its independently mutable Python children and
calls Reaction::try_new on the assembled value. This study measured paths retaining
that interface; it did not establish that the interface should be preserved.
The live-child contract differs from closed Rust Reaction storage. Its correction
is scheduled in [231](231-python-api-parity-2026-09-21.md#next-work-after-181--reaction-ownership)
immediately after this document's current work.
Measure reconstruction and integrity checking separately from iterator creation
and full product realization. Inspect property access and matching/application
call paths to count reconstruction boundaries. Compare the existing application
path with borrowing the supplied host until Rust creates the iterator's required
snapshot, removing the redundant binding-side host copy. Both candidates retain
fresh reaction construction and validation. Use stated reaction/delta and host
sizes, assert identical products and source preservation, and retain snapshot
semantics when inputs change after iterator creation. If avoiding reconstruction
itself requires an API change, record that boundary rather than prototype an
unapproved child-view or cache design.

#### Results — 2026-10-08

**Decision supported:** borrowing the host removes an unnecessary clone without
changing behavior, but this is a small cleanup. The storage correction in 231 is
justified by the Rust contract, independently of these timings. Allocation counts and requested bytes were identical between
the two application routes: the extra host clone shares storage. Whole-operation
timing differences were small, and the smallest fixture was slightly slower in
the candidate. Do not infer a reliable percentage speedup from this run.

| Pattern / host atoms; delta count | Clone + validate reaction, µs | Current issue + products, µs | Borrow-host issue + products, µs | Allocations / requested KiB, either route |
| --- | ---: | ---: | ---: | ---: |
| 2 / 8; 1 | 0.08 | 2.99 | 3.03 | 37 / 8.25 |
| 2 / 128; 1 | 0.08 | 25.07 | 24.17 | 165 / 36.93 |
| 16 / 128; 15 | 0.28 | 33.79 | 33.42 | 240 / 83.57 |

Fixtures are oxygen-anchored carbon chains with one exact product, changing one
or fifteen bond orders. Clone+validation takes 0.08–0.28 µs here and copies only
the delta vector (one allocation, 0.20–2.93 KiB). The check on precloned inputs
allocates nothing in these cases. This does not establish costs for large
constraint-rich reactions.

Each apply-family call constructs/checks one reaction snapshot. Rust then owns
reaction/host snapshots in its iterator. Property getters lhs/deltas return the
live children without reconstruction; equality/composition reconstruct both
operands. Fresh checking cannot simply be removed while those children remain
independently mutable. The correction in 231 removes that mutation contract rather
than retaining repeated checks. Direct Rust Reaction storage was not measured in
this experiment.

Both routes produced the independently constructed expected molecule and preserved
inputs. After iterator creation, changing the host and live reaction children did
not change the issued product; fresh conversion rejected the deliberately invalid
new child combination. Timing includes boundary issuance and Rust product
realization/destruction, not Python product wrapping or outer Python dispatch.
The source-level host-borrow change also applies to tracked_apply,
apply_to_reaction, and apply_to_reaction_span; those were not separately timed.

### Study 3 — Edit and Delta conversion

**Question:** separate necessary independent-entry copying from repeated boundary
conversion. Scope: Edit, Edits, Delta, Deltas, and their entity payload conversions.
Edits consumption remains agreed; this study does not reopen it.

Compare reconstructing a batch through Python entries with directly cloning its
Rust container, implementing the approved explicit Edits copy operation only in
the experiment. The analogous Deltas copy method is approved as well. Owned whole-container from_rust conversions already move their input;
indexing a retained container cannot instead move its stored entry.

Measure construction/insertion, indexing, iteration, and application separately,
using small field changes, additions with substantial attributes, and constraints
at stated sequence lengths. Check whether owned intermediate results pass through
borrowed entry conversions and copy unnecessarily. Assert exact converted values
and equivalent application results before measuring. Preserve independent-entry
and iterator semantics; do not introduce entry-view types to reduce a benchmark.

#### Results — 2026-10-08

**Decision supported:** implement the approved Edits copy method by directly
cloning the Rust container. Reconstructing through Python entries is substantially
more expensive. Deltas shows the same effect; its explicit copy method is also approved. Neither current Edits application nor
Reaction conversion reconstructs individual Python entries: application transfers
Rust Edits directly, and Reaction clones Rust Deltas directly.

Times in µs, with both copying routes invoked through Python dispatch:

| Payload / entries | Edits reconstruction | Edits direct copy | Deltas reconstruction | Deltas direct copy |
| --- | ---: | ---: | ---: | ---: |
| field, 32 | 8.20 | 0.53 | 9.83 | 0.54 |
| field, 256 | 62.36 | 3.35 | 74.88 | 3.18 |
| addition, 32 | 9.97 | 2.77 | 10.64 | 2.18 |
| addition, 256 | 85.05 | 26.49 | 84.64 | 19.02 |
| constraint, 32 | 6.56 | 0.54 | 8.21 | 0.59 |
| constraint, 256 | 49.50 | 3.47 | 62.77 | 3.80 |

The field fixture changes scalar charge; the constraint fixture changes scalar
atom valence; additions carry an owned element-variable string. Sizes 1, 32, and
256 were exercised. These are payload cost probes, not a workload distribution.
For 256 additions, Edits reconstruction/direct copy used 1,544/769 Rust allocation
calls and 242.50/140.75 KiB requested bytes. Copy followed by actual Python
Molecule.apply took 93.21/34.41 µs, with exact equal complete molecules. The method
keeps its transfer semantics; copying is explicit.

Indexing, iteration, construction from existing Python entries, and append were
also measured. For 256 additions, Edits construction from Python entries took
34.49 µs, repeated append 63.69 µs, and full traversal 52.99 µs. One indexed entry
cost about 0.19 µs. No move-only alternative preserves independent returned
entries and reusable supplied arguments. These costs do not justify introducing
new entry views in this study. Figures use the current production enum layouts,
before the planned boxing/batching changes.

A distinct removable copy occurs in AtomDelta inverse: Rust creates an owned
inverse, then the current borrowed from_rust conversion copies its attributes
again. The prototype still calls Rust inverse, then moves the resulting Add/Remove
attributes into the same read-only Python form representation.

| Owned string bytes | Add inverse current / moved, µs | Remove inverse current / moved, µs | Rust allocation calls current / moved | Requested KiB current / moved |
| --- | ---: | ---: | ---: | ---: |
| 16 | 0.29 / 0.22 | 0.29 / 0.23 | 4 / 2 | 0.14 / 0.07 |
| 256 | 0.29 / 0.23 | 0.30 / 0.23 | 4 / 2 | 0.61 / 0.30 |
| 4096 | 0.45 / 0.29 | 0.44 / 0.30 | 4 / 2 | 8.11 / 4.05 |

**Proposed correction:** move owned inverse results during binding conversion;
retain copies where the source is genuinely borrowed. Exact inverse values and
unchanged inputs were checked. This experiment covers AtomDelta Add/Remove only;
other entity inverse methods have the same borrowed-result conversion pattern,
but their payload-specific savings were not measured. No Python-side inverse
algorithm or form-consumption mechanism is introduced.

**Execution and verification:** isolated source copy under
scratch/181-python-studies at ec82fcdcc; Rust 1.96.0 release, arm64, Python 3.13.15
from umol-py/.venv. A native Rust test binary embeds Python; the installed extension
was not replaced. Production sources are unchanged. Nine focused rstest cases
passed in each of the timing and allocation runs; independent source checks
reviewed candidate semantics and measurement boundaries.

Timing ran serially after compilation: five warmups, fifteen sample means, median
reported. Per-sample repetitions were 100 for model conversion/resolution (1,000
for bare borrowing), 64 for reaction copying, 16 for isolated validation and matching/products,
128 for entry operations, and 512 for inverse operations. Fixture creation is outside the
measurements; result destruction is included. Separate warm-call allocation runs
count Rust allocation/reallocation calls, requested bytes, and peak extra live
bytes. They exclude CPython's heap and RSS. Precloned/preissued isolated reaction
rows free storage prepared outside measurement, so their peak figures are not
comparable to complete-operation rows.

Reproduction and complete 224-row data: scratch/181-python-studies/README.md,
run.py, summary.csv, the per-study sources/notes, and time/alloc logs. These tables
retain the decision-relevant findings independently of scratch. No broad Python
API audit or production migration is claimed complete by these studies.

### Study 4 — Form operations and preliminary operand copies

**Question:** can Python form operations borrow Rust operands rather than copy or
reconstruct them before invoking borrowed Rust methods? Keep source values usable,
result types, errors, and exact operation semantics unchanged. This is a bounded
conversion study, not a redesign of form construction or a new binding framework.

**Two cases:** all eight entity forms and their eight entity constraint containers
already retain Rust values. The molecule-level Constraints container does too.
Their macro conversions clone those values before borrowed operations. These
17 wrappers can avoid those preliminary clones without changing storage or the
Python interface. Leaf forms and individual constraint forms commonly represent
Python variants instead; avoiding reconstruction there requires a separate
representation decision. A conversion cannot borrow a Rust value that is not stored.

#### Results — 2026-10-08

**Approved:** remove preliminary operand clones from is_ground,
is_undetermined, matches, is_compatible, meet, join, and normalized_eq for the
16 entity form/container wrappers. For Constraints, remove them from normalized_eq.
Keep the copying needed to produce independent results; normalize remains
source-preserving. This recommendation changes binding conversion only, not Rust
lattice operations, normalization, integrity checking, or mutation semantics.

Python-call times in µs and Rust allocation calls, current → borrowed:

| Fixture / operation | Time | Allocations |
| --- | ---: | ---: |
| AtomForm, charge set of 16 / is_ground | 0.182 → 0.059 | 4 → 0 |
| AtomForm, charge set of 16 / normalized_eq | 0.398 → 0.167 | 8 → 0 |
| AtomForm, charge set of 16 / meet | 1.040 → 0.792 | 25 → 17 |
| AromaticSystemForm, six electron counts / normalized_eq | 0.229 → 0.134 | 2 → 0 |
| AromaticSystemForm, six electron counts / meet | 0.271 → 0.219 | 3 → 1 |
| MulticenterBondForm, four electron counts / normalized_eq | 0.198 → 0.134 | 2 → 0 |
| AtomConstraintsForm, valence set of 16 / normalized_eq | 0.372 → 0.155 | 8 → 0 |
| Constraints, charge sum over 16 atoms / normalized_eq | 0.232 → 0.134 | 4 → 0 |

The AtomForm normalized_eq case eliminates 864 requested bytes per call; the
aromatic case eliminates 96 bytes. Small allocation-free/default forms offer much
less benefit, and some timings were slightly worse. These results justify removing
unnecessary operand copies, not a claim that every form operation becomes faster.
Normalizing the populated AtomForm still makes five allocations (456 requested
bytes) in either route; the prototype preserves its independent result.

Experimental stored-Rust wrappers for NumForm, ElementForm, AtomConstraintForm,
and MoleculeConstraint measured the potential benefit of avoiding reconstruction:

| Fixture / operation | Current / stored time, µs | Current / stored allocations |
| --- | ---: | ---: |
| NumForm literal / normalized_eq | 0.130 / 0.130 | 0 / 0 |
| NumForm set of 16 / normalized_eq | 0.343 / 0.157 | 8 / 0 |
| ElementForm variable / normalized_eq | 0.205 / 0.133 | 4 / 0 |
| AtomConstraintForm valence set of 16 / normalized_eq | 0.342 / 0.162 | 8 / 0 |
| MoleculeConstraint charge sum over 16 atoms / normalized_eq | 0.178 / 0.138 | 2 / 0 |

These are not complete replacement bindings: constructors, variant access,
nested mutation, and their costs were not migrated or measured. They establish
where reconstruction costs occur, not that changing every leaf type is worthwhile.
No leaf representation migration is selected. Scalar literals already avoid
allocations; source-preserving normalization kept the same allocation counts.

**Coverage and limits:** all eight entity forms and all eight entity constraint
containers were exercised with defaults; populated cases covered AtomForm,
AromaticSystemForm, MulticenterBondForm, AtomConstraintsForm, and Constraints.
Additional leaf/individual-constraint cases included numeric sets of 16 and 128,
a numeric literal, an element variable, a valence constraint, and a molecule
charge-sum constraint. Binary calls used distinct structurally equal operands.
Checks compared complete Rust results, Boolean results, and ordinary absence.
This does not establish behavior for every unequal operand or failure branch,
nor representative application-level speedups.

**Execution:** isolated scratch workspace at ec82fcdcc; Python 3.13.15, Rust 1.96.0,
release arm64. The focused rstest passed in both timing and allocation runs:
204 operation comparisons, 408 rows per run. Five warmups, fifteen sample means,
1,024 calls per sample; median reported. Measurements include Python method
dispatch and result conversion/destruction, exclude fixture construction, and
count Rust allocations separately from timing. Python heap/RSS were not measured.
Files are under scratch/181-python-studies/forms: generate.py, run.py, README.md,
summary.tsv, build.log, and forms-{time,alloc}.log. Production code and the
installed extension are unchanged.

**Separate finding:** StereoAtomForm and StereoBondForm lack Python value equality;
the other six entity forms implement it. The harness therefore compares complete
Rust results. The binding defect is recorded in
[231](231-python-api-parity-2026-09-21.md#stereo-entity-form-equality--2026-10-08),
not repaired or treated as a conversion-performance issue here.

## Migration interface review — 2026-10-08

The ownership and recovery contracts above remain the basis of the migration.
The following makes the remaining interface questions concrete; proposals in this
section are not approval to implement new names or public signatures.

### Scoped editing

Settled scoped editing signatures:

```rust
// Molecule
pub fn edit_into<T, E>(
    &mut self,
    f: impl FnOnce(&mut MoleculeEditor) -> Result<T, E>,
) -> Result<T, E>
where
    E: From<MoleculeApplyError>;

pub fn apply_into(&mut self, edits: Edits) -> Result<(), MoleculeApplyError>;
```

The callback receives a mutable borrow of the single owning editor. It cannot
consume that editor with finish. After callback success, the method calls finish,
converts integrity failure through MoleculeApplyError, and installs only success.
Return the callback's T after installation. Error or callback unwinding leaves the
receiver untouched. This uses the same error-conversion bound as Transaction::run;
it needs no new error enum or editor mode. The names edit_into and apply_into and
the callback signature are settled.

### Transaction entry and higher-level selection

The already selected RollbackMethod belongs to graph-IR transaction handling.
The settled signatures keep the callback last on run and put batches before
execution policy on the prepared-batch conveniences:

```rust
// Transaction
pub fn run<T, E>(
    molecule: &mut Molecule,
    rollback_method: RollbackMethod,
    f: impl for<'scope> FnOnce(Transaction<'scope>) -> Result<T, E>,
) -> Result<T, E>
where E: From<MoleculeApplyError>;

pub fn commit(self) -> Result<(), MoleculeApplyError>;

// TrackedTransaction
pub fn run<T, E>(
    molecule: &mut Molecule,
    rollback_method: RollbackMethod,
    f: impl for<'scope> FnOnce(TrackedTransaction<'scope>) -> Result<T, E>,
) -> Result<T, E>
where E: From<MoleculeApplyError>;

pub fn commit(self) -> Result<MoleculeCorrespondence, MoleculeApplyError>;

// Both scoped handles
pub fn apply(&mut self, edits: Edits) -> Result<(), MoleculeApplyError>;
pub fn probe(&self) -> Result<&Molecule, MoleculeIntegrityError>;
pub fn rollback(self);

// Molecule
pub fn transact(
    &mut self,
    batches: impl IntoIterator<Item = Edits>,
    rollback_method: RollbackMethod,
) -> Result<(), MoleculeApplyError>;
pub fn tracked_transact(
    &mut self,
    batches: impl IntoIterator<Item = Edits>,
    rollback_method: RollbackMethod,
) -> Result<MoleculeCorrespondence, MoleculeApplyError>;
```

Both handles are operation-issued, with private fields and no public constructor
other than their scoped run entry point. They borrow run-owned state and cannot
replace the receiver or recovery storage. They are not copyable or cloneable and
have no conversion or mutable access to the other handle type. In particular,
tracked application cannot bypass correspondence updates through an exposed
ordinary handle. Export TrackedTransaction beside Transaction through molecule
and ir; add no public tracking trait, mode flag or tracking-disabled error.

Each type preserves the existing callback lifetime, separate batch namespaces,
abort latch and recovery contract. apply checks edit preconditions and can leave
transient state; probe and commit check representation integrity. Probe returns
MoleculeIntegrityError; application and commit failures use MoleculeApplyError,
including the existing TransactionError variants. No new failure condition is
introduced by selecting tracking. Failed probe leaves the transaction active;
failed application or commit restores entry state and aborts. Explicit rollback,
uncommitted return, callback error, forgotten handle and callback unwinding
between completed edits retain the existing restoration guarantees.

run returns the callback's T for both types. A tracked commit returns its mapping
to the callback; it does not independently accept changes. Explicit rollback and
return without commit do not request a correspondence, so run needs no additional
optional mapping result. For equivalent edit sequences, discarding the tracked commit result
preserves the ordinary final state and errors. Correspondence covers all eight
entity-id spaces, composes across batches and does not map participant positions.
Regrouping edits preserves the mapping when handles name the same entities in
each batch's namespace. These scoped handles have no DSL or serialization form.

Python transact/tracked_transact expose RollbackMethod through a keyword-only
rollback_method argument defaulting to Journal. Clone remains selectable. Both
methods provide the same recovery guarantee with either choice. Journal is the
default because its recovery storage scales with changes rather than molecule
size; it is not a claim that journaling is always cheaper.

Each mutating _into method except transform_into and the individual resolvers'
_into methods takes RollbackMethod last, after its existing inputs (including
ProjectFlags). Rust callers select explicitly. Python resolve_into and
resolve_into_with_report expose a keyword-only rollback_method defaulting to
Clone; Journal remains selectable. This is not a universal performance claim.
Python transact retains its separately agreed Journal default.
Python continues to expose prepared-batch operations; neither scoped Rust handle
is introduced as a Python callback interface.

### Correspondence without Undo payload retention

ApplicationState already exists for every apply call: it resolves initial and
new handles within that batch. Its update_correspondence method can incorporate
the batch's removals and surviving additions. Reuse that behavior rather than
reconstructing mappings from old attribute and constraint payloads.

TrackedTransaction maintains one MoleculeCorrespondence directly in both
recovery modes. Initialize it to identity at tracked transaction entry and update
it after each successful batch through ApplicationState::update_correspondence.
Drop each batch's ApplicationState afterward. TrackedTransaction::commit moves
out the accumulated value. Transaction has no aggregate mapping to construct,
update or discard. Neither type retains batch mapping history or uses Undo to
reconstruct correspondence.

S0e measured a real cost for unconditional tracking, including a 32.8% increase
for a plain single charge edit on a unique 96-atom chain and 1,528 bytes of
retained mapping storage. S0e2 compares the split implementation with the saved
pre-S0e Journal baseline and S0e's eager candidate, including identity construction
and mapping updates only on the tracked path. S0f measures both types with Clone
and Journal recovery. The S0e measurements remain historical evidence; they do
not establish performance of the new split.

### Transaction guard

The approved recovery storage is a private enum in molecule/transact:

```rust
enum TransactionStorage {
    Journal(Vec<Undo>),
    Clone(Option<Molecule>),
}
```

The retained ownership split is as follows; the S0e2/S0f costs and decision to
proceed are recorded below.
Both run entry points own a private TransactionGuard retaining the receiver's
mutable borrow, TransactionStorage and TransactionStatus. Tracked run also owns
the accumulated correspondence in an Option, borrowed by its small callback
handle and consumed with take at commit. The run-owned guard restores unaccepted
changes even if a handle is forgotten. Cleanup of deliberately forgotten handle
resources is not part of the contract. Share recovery and acceptance logic
privately; neither callback handle owns recovery. The direct-ownership and lazy
identity experiments below failed the timing gate and were withdrawn.
Journal applies edits to the receiver with undo recording. Clone applies edits
without undo recording to its working copy, leaving the receiver untouched.

On application or commit failure, Journal immediately replays undo; Clone drops
the working copy. Both mark the transaction Aborted. probe then observes the
restored or untouched receiver. Successful commit requests acceptance; only a
successful callback return accepts the Journal changes or moves the Clone working
copy into the receiver. Explicit rollback, an uncommitted return, callback error,
forgotten handle, or unwinding discards unaccepted changes through the guard.
The existing TransactionStatus variants suffice. The Option permits taking the
working molecule on acceptance or discarding it on failure; it never places an
empty molecule in the receiver.

### Consuming execution and value traits

Current ingestion, parsing, and export call sites are in umol-graph's ingest.rs,
parse.rs, and export.rs. They are not cross-crate umol-io callers. Consequently,
retaining no-recovery execution for those owned inputs does not require new public
consuming resolver/projector methods.

The public consuming resolve, resolve_with_report, and project serve those callers
directly; no crate-private consuming functions are needed.

For individual single-plan resolvers and transformers, keep planning against the
source and apply to one editor; do not introduce a second execution abstraction.
transform_iter clones its source lazily, once per result, and calls transform.

Proposed scope disposition for Normalize, Reframe, Canonicalize, builder completion,
and consuming construction-data conversions: retain existing Rust ownership
contracts. The molecule/editor/resolver migration does not require changing these
traits. Python copy-preserving operations remain explicit binding decisions;
no blanket consuming-to-borrowed trait rewrite is included. [236](236-operation-ownership-copy-costs-2026-10-09.md) owns the
ownership contracts of these traits; it removes Normalize::normalized and makes
canonical_hash borrow.

### Kindless stereo constraints — decision, 2026-10-14

The stereo kind is the subtype of the entire stereo entity. Its constraints
derive their meaning from that owning kind; they do not carry an independent
intrinsic kind. For both StereoAtomForm and StereoBondForm, constraints are
semantically vacuous when the configuration is Undetermined (`"*"`). A known
kind with an undetermined coset (`"Th*"`, for example) retains meaningful
constraints and does not fall under this rule.

**Raw storage and mutation.** Attribute mutation and form updates remain
infallible. They may retain constraints while the kind is undetermined; assigning
or clearing a kind does not eagerly normalize or erase the stored payload. Raw
equality and hashing continue to observe it. A later raw kind assignment can
make retained data meaningful; normalization or serialization may already have
discarded that data. Raw mutation is not required to commute with normalization.

**Normalization and lattice operations.** Normalization of a kindless stereo
entity clears its constraints before inspecting or normalizing their contents.
An internally contradictory but vacuous constraint therefore cannot make the
entity contradictory. Known-kind entities retain the existing normalization of
their configuration and constraints. The rule belongs to the owning entity;
a standalone constraint container has no kind context of its own.

normalized_eq compares the resulting normal forms. is_undetermined, is_ground,
meet, join, matches, and is_compatible must use the same vacuity rule. In
particular, meeting a kindless entity with a kinded entity must not transfer the
kindless operand's stored constraints into the result. A join whose resulting
kind is undetermined likewise has vacuous constraints. Results must agree
semantically whether the operands were normalized first or not. Entity-bound
semantic constraint readings likewise ignore kindless stored constraints and
return top at the requested key: Some containing that key's undetermined form.
None continues to mean bottom. Raw storage access remains available for
observing and mutating the retained payload.

**Updates and lowering.** Keep the uniform infallible update interfaces,
including Edits::update_stereo_atom and Edits::update_stereo_bond. No additional
kind argument or stereo-specific rejection path is introduced. An Unchanged
configuration update does not assert an undetermined kind: when a supplied
current form has a kind, that context remains available for constraint changes.
A Kinded update supplies context explicitly; an Undetermined update clears the
configuration, making constraints in the resulting entity vacuous. Lowering
uses the resulting configuration's kind; an explicit Undetermined update does
not fall back to the source kind.

Semantic differences and lowering must respect the applicable source/update
kind context. Two kindless forms differing only in stored constraints have no
semantic difference. Ignore constraint contributions for which there is no
determined kind context, while retaining actual configuration changes and
constraint changes with known context. Ignoring a constraint must not become a
constraint-removal operation. Retained edits still use the stored source values
for their old-state preconditions; semantic vacuity is not permission to replace
a stored old value with absence. This decision adds no kind-presence or kind-match
rejection to edit application.

**Serialization.** Entity forms, attribute updates, and single-entry and batched
edits use the same kind-scoped constraint syntax and vacuity convention. Omit
constraints with no determined kind context instead of returning a formatting
error or reaching an expect/panic during conversion. A kindless entity renders
as `"*"`. An update retains its configuration intent: Unchanged and an explicit
Undetermined update remain distinct after their vacuous constraints are omitted.
Use an available supplied current-form context during lowering; standalone
serialization does not invent or infer a kind from a molecule.

For a kindless constraint-only edit or batch row, omit the entire constraint
operation, including both its expectation and update. Do not encode just its new
value as absent, which would request removal. Preserve unrelated configuration
and other edits and the order of retained operations. Raw Rust/Python value
conversion can retain the payload; DSL export applies the semantic omission.
The entity roundtrip is under normalized_eq for this vacuity, not raw equality;
it does not promise to preserve discarded payload for later raw kind assignment.

Harmonize the existing behaviors: entity serialization already omits kindless
constraints, whereas update formatting and edit conversion currently fail on
them. Apply the same omission before those failure paths. Keep the existing
kind-scoped constraint grammar; do not add a separate raw encoding or require
fallible mutation. This does not weaken the kind requirement of molecule-level
stereo constraint leaves.
The shared LigandSymmetry and Fluxionality batch-key syntax remains a separate
unresolved design point in S2; full permutation-image encoding remains rejected.

## Python getter contract review — 2026-10-14

This review records existing Python getter behavior. The naming and immutable
observation changes are assigned to the separate S11c implementation item. They
are not prerequisites for S2a or the other batch subitems.

**Decision:** the concrete getter changes listed in S11c are approved. Retain the existing variants
and constructors. Added Python payload names are limited to names already
established by the corresponding public Rust fields or accessors; unnamed tuple
fields retain positional access. Every payload must be accessible on Python 3.9
without structural pattern matching. Pattern matching remains an optional
convenience on newer Python versions. S11c implements this alignment; designing
names for unnamed Rust payloads is excluded. S2 reuses the existing Python forms
and their getter behavior without this migration.

Reviewed source: fa093b068cc25668661a232cdfaa353a5b0af166. Python sources were
unchanged; the unfinished S2a changes in graph-IR edit, application, and DSL modules
were excluded. Source inspection covered explicit properties, generated enum
properties, nested return values, and the mutation paths described below. Python
3.13.16 probes used the installed extension and agreed with the source traces; the
binary has no embedded revision stamp, so its exact source revision is not asserted.
This is a getter/mutation review, not a full operation, lifecycle, or parity audit.

### Required behavior and working mutation paths

The guide permits both immutable values and writable accessors. It does not require
every property to be immutable. A property result must either reject mutation or
make supported writes affect the backing object. A detached mutable result that
accepts ineffective writes violates that rule. The absence of a setter on the
outer property does not establish recursive immutability.

| Access | Current intended role and supported writes |
| --- | --- |
| An entity field such as AtomForm.charge or AtomView.charge | Returns a leaf value; the field setter replaces it with another complete form. AtomView writes to its molecule. The returned leaf is not a live mutable field accessor. |
| The constraints property of each entity form or entity view, for all eight kinds | Returns a live constraint view. Its setters and set/pop/update methods write to the backing entity. |
| Atom/Bond/DativeBond constraint ring_size_count | Returns a live indexed accessor; assignment and deletion update the backing constraints. |
| Molecule's eight entity collections | Returns live accessors; item assignment replaces entity attributes, and returned entity views support attribute assignment. |
| Entity forms retained by Edit or Delta | Field and constraint writes reject with TypeError; copy explicitly returns an independent writable form. |
| Molecule.constraints and logical Constraint.And/Or children | Returns read-only ConstraintsView accessors. Nested leaf payloads inherit the leaf issues below. |
| Reaction.lhs/deltas and Solution.Determined.molecule | Returns retained mutable objects; changes are visible on subsequent access. Frozen Python field slots do not make these payloads recursively immutable. |

The same field-replacement mechanism handles literals, sets, ranges, and expressions.
For example, these operations were executed successfully on AtomForm, and the
corresponding assignments through Molecule.atoms[0] changed the molecule:

```python
atom.charge = NumForm.LitSet({-1, 0})
atom.charge = NumForm.LitSet(set(atom.charge._0) | {1})
atom.charge = NumForm.ArithExpr(
    ArithExpr.Sum([ArithExpr.Var("q"), ArithExpr.Lit(1)])
)
atom.constraints.valence = NumForm.LitSet({3, 5})
```

Numeric setters accept int or NumForm, not a bare Python set. Leaf variants have
no supported in-place mutation operation: their payload slots reject assignment.
Where supplied, normalize/meet/join return new values. AtomUpdate is also a value;
Edits.update_atom applies such an update using the supplied current AtomForm.

### Confirmed getter deviations

**Detached mutable collections.** Rust collections exposed through generated or
explicit Python getters become lists, sets, or dictionaries. Mutating the returned
collection succeeds while the originating value stays unchanged. The strongest
defense is storage isolation: these paths do protect the original value. That does
not satisfy the guide's requirement that immutable observations reject mutation.
Independent review confirmed this distinction.

| Family | Affected property inventory |
| --- | --- |
| Leaf forms: 19 getter sites | ElementForm.LitSet/NotSet and Var restrictions; IsotopeMassForm.LitSet and Var restrictions; NumForm.LitSet; ArithExpr.Sum/Product; PredExpr.Mem/And/Or; ElectronCountsForm.Lit; StereoTerm.Var/LitSet; StereoCoset.LitSet; TopicityRelationForm.LitSet/NotSet; StereogenicityForm.LitSet/NotSet. |
| RelationalConstraint/MoleculeConstraint: 14 getter sites | Ten relational participant/predicate arrays and vectors; ChargeSum, BondOrderSum, Connected, and UnpairedElectronCoupling participant lists. |
| Edit: 25 properties across 19 variants | Atom/bond addition rows, topology removal lists, relation participants, removal rows, and old/new participant lists. Nested donor/atom/ligand lists also remain mutable copies. Retained entity forms themselves correctly reject writes. |
| Entity Delta: 10 participant properties | Add/Remove donors for dative bonds, atoms for aromatic systems and multicenter bonds, and ligands for stereo atoms and bonds. |
| Mapping and compaction values | Correspondence.matched_pairs/left_unmatched/right_unmatched, Remapping.images, Compaction.removed. Their molecule-level wrappers inherit these nested observations. |
| Models and reports | ValenceEntry.target_covalences/aromatic_valences/fallback_aromatic_valences; StereoModel.kind_models; ElementScope.AllowList.elements; ResolveReport.tie_breaks. |
| Fingerprint values | HashedFeatureSet.ids, CountedHashedFeatureSet.entries, StructuralFeatureSet.keys, SignedHashedFeatureSet.entries, RoleTaggedHashedFeatureSet.ids. |
| Entity topology views | DativeBondView.donors; StereoAtomView.ligands and StereoBondView.ligands. Other inspected entity atom_ids properties return tuples. |

These counts identify direct getter sites, not every nested path to them. Field
changes, individual constraints, entity attributes, and read-only Edit/Delta
payloads also expose the affected leaves. S11c specifies the collection return
types and the exact boundary of this correction. Reaching a leaf through another
getter does not authorize changing that getter's ownership or mutation contract.

**Detached mutable constraint containers.** All eight entity Update.constraints
getters clone into writable ConstraintsForm objects. For example,
`child = update.constraints; child.valence = 4` changes child but leaves update
unchanged. This is confirmed by source and execution, and independently challenged.
It needs an immutable observation consistent with the Update's value role; treating
the copy as writable is not justified by the current guide. Ordinary entity
constraints getters are different: those already write through correctly.

**Detached mutable aggregate children.** ReactionLayout.lhs/rhs clone into writable
MoleculeLayout objects, and ReactionMetadata.lhs clones into writable
MoleculeMetadata. Their mutators succeed without updating the parent. Comments and
tests describing snapshots explain the implementation but do not justify these
property semantics under the guide. These mutable aggregate ownership questions
are excluded from S11c; this finding is not a prerequisite for that item.

The inspected scalar/value-only properties, read-only ring member/connection
accessors, and MoleculeLayout.positions tuple do not show these collection defects.
Mutable ring models and independently owned ring connections have explicit setters.
Explicit copy/asdict operations and value-producing projection methods are not
classified as defective merely because their results are mutable.

### Coverage and review limits

The source inventory covered all Python binding modules, including generated
properties. Runtime probes exercised all 19 leaf and 14 constraint collection
sites and representative aggregate and live mutation paths. Separate reviewers
covered leaf/constraint values and other aggregates; an independent challenger
checked the findings. Complete inventories, source locations, probe bodies and
results are in [leaf getter review](../scratch/181-s2a/leaf-getter-review.md),
[aggregate getter review](../scratch/181-s2a/aggregate-getter-review.md), and
[entity mutation probes](../scratch/181-s2a/entity-getter-probes.txt). The
[independent challenge](../scratch/181-s2a/getter-challenge.md) distinguishes
source-confirmed findings from individually executed cases.

Two adjacent read-method families were identified by source inspection only:
AtomTypeRegistry.patterns_for_element/patterns_for_element_and_charge and
AtomCompletions.get/items copy borrowed forms into writable objects. Their methods'
contracts need separate reconciliation; they are not included in the property
counts. This review does not approve the distinct RingView.id mutation contract,
nor does it establish performance, full invalidation behavior, or general parity.

The existing test_edit.py test at lines 233–237 explicitly allows clearing a
returned atom list and checks only that the stored Edit is unchanged. This is a
concrete test/review blind spot: storage protection was tested instead of rejection
of mutation through immutable observations. It does not establish the reasoning
used in each earlier review. The review failed to distinguish that broader getter
contract from the batch implementation; its correction is assigned to S11c.

## Implementation plan

### S0 — Molecule, editor, and transaction migration

[236](236-operation-ownership-copy-costs-2026-10-09.md) is completed; this plan
remains the 0.10 release blocker. Individual resolvers and transformers complete
planning before mutation and apply without transactions. AromaticityPerceiver
does not use transactions. Composite Resolver's borrowed operations retain one
transaction across phases because later planning can reject earlier changes.
S0f migrates those transaction calls; S10a exposes their recovery selector.
Individual operations still migrate with editor/application interface changes
in S0b and plan emission changes in S3.

The baseline precedes all production changes. The migration subitems follow.
The migration interface review distinguishes settled interfaces from remaining
proposals; an unresolved interface must be specified before its item starts.
Each implementation subitem includes its callers and focused tests; none is a
separate task merely to write another plan.

#### S0a — Benchmark baseline (done)

**Modules:** existing graph-IR editor/reaction/constraint_readings benches, graph
resolve benches, bounded operation/Python measurements, and this record.
**Measurement, green.**
[dep: 236 complete; precedes every production-change subitem here]

Run the editor, reaction, constraint_readings, and resolve benchmark targets before
changing production code. Record the source revision, toolchain, build settings,
fixture sizes, ownership mode, and measurement boundaries with the results in section 4.
Retain a named Criterion baseline for subsequent comparisons; raw outputs may live
in scratch but are not the only record of the results.

| Existing target | Baseline cases |
| --- | --- |
| graph-IR `editor` | sparse8/dense80, unique/shared, direct/apply/transact, and six overlay-removal kinds in stored/reordered frames |
| graph-IR `reaction` | supplied-match and enumerated application, tracked counterparts, editor/application cases, combine/split/meet_pushout, and reversal |
| graph `resolve` | ingestion, resolution, individual aromaticity/stereo projection, composite projection, and sparse8/dense88 ownership cases |
| graph-IR `constraint_readings` | stereo-atom and stereo-bond readings, normalized_eq, meet, join, matches, and is_compatible on bounded literal/undetermined cases |

Extend the existing constraint_readings timing/allocation harness for the stereo
cases before production changes; reuse its Criterion and allocation-count paths.
Cover known kinds with literal or undetermined cosets, empty constraints, and
literal assertions, plus kindless forms with and without stored payloads. Keep
unusual payloads in separately reported cases; do not combine them into an
assumed workload distribution or let sets/expressions dominate the comparison.
Record time, allocation count, and allocated bytes with fixed preparation and
destruction boundaries. Check results outside timing and distinguish kindless
cases whose expected results change from cases with unchanged semantics. S2g1
owns the comparison against this baseline.

Reuse the bounded recovery study's fixtures for all three transformers, bounded
failure, and separate allocation/peak-memory measurements absent from Criterion.
Adapt its harness to the implemented consuming and in-place complete-plan paths
of individual resolvers and transformers. Scope clone/install-versus-journal
comparisons to composite operations and prepared transactions. Keep historical
measurements tied to their recorded revisions; fresh baselines use the completed
236 paths. Verify the harness's source dependencies against the current production
revision before running; do not rerun the boxing-layout matrix.
Use the existing Python-study harness for actual Molecule edit/apply/tracked_apply,
prepared transactions, and resolve calls, with default and supplied models; add
missing call cases to that harness rather than another benchmark framework.
Include direct entity and aggregate Delta inverse allocation baselines for S6c.
Most fixtures have no molecule-level constraints.

Check result and recovery agreement outside timing. For each row specify whether
input/Edits preparation and final result destruction are timed; keep those
boundaries for the final comparison. Include any clone required by the measured
source-preservation contract inside the operation boundary. Record current
consuming calls separately from source-preserving calls: a changed method name
must not silently change what a benchmark measures. Python uses the repository's
3.13 environment and a freshly built extension. Run measurements sequentially.

**Completion:** actual timing and memory results, commands, and fixture definitions
are recorded durably here, with reusable inputs for S2g1/S4a/S12a. This item produces
the pre-change baseline; S12a performs the final comparison.

Completed 2026-10-14. The [S0a measurement record](#s0a-baseline--2026-10-14)
contains the current timing/allocation results, fixture and ownership boundaries,
verification, and retained baseline artifacts. Production code is unchanged.

#### S0b — Mutable editor batch application (done)

**Modules:** graph-IR molecule/editor/transact, all Rust consumers, Python editor
binding, and affected tests. **Breaking, red→green within the subitem.** [dep: S0a]

Change MoleculeEditor apply to take &mut self and return Result<(),
MoleculeApplyError>; tracked_apply takes &mut self and returns only
MoleculeCorrespondence. Keep direct mutation, checked probe, and consuming finish.
Failure retains the completed prefix without rollback. Each application still
starts a fresh handle namespace. Remove editor reassignments/chaining at callers;
do not replace them with clones or extra publication checks.

Change Python editor apply to return None and tracked_apply to return only the
correspondence. Borrow inner mutably; only finish takes the Option value.
Validate editor availability before transferring Edits. Keep the editor usable
after application failure, including its partial state; consumed Edits stay
consumed. Finish consumes on both success and integrity failure.

**Caller migration in this subitem:**

| Current caller | Required change |
| --- | --- |
| graph-IR Molecule apply/tracked_apply in `ir/molecule.rs` | Hold a mutable editor; apply returns unit, tracked_apply returns only the correspondence; finish that editor. Molecule methods keep consuming. |
| graph composite Resolver resolve/resolve_with_report/project in `ops/resolve.rs` | Remove every `editor = editor.apply(...)` reassignment; retain phase order, probe locations, planning-error mappings, application/publication assertions, and final finish. |
| graph individual resolvers: aromaticity/isotope/stereo resolve and project; bonds/multicenter resolve | Split edit/apply/finish chains into a mutable editor, application, then finish. Preserve current public receivers and result types. |
| graph Kekulizer/Aromatizer/DelocalizeCharge transform | Split editor chains likewise; retain Kekulizer's validation during planning and existing plan-application/publication assertions. |
| Python `src/transaction.rs` | Use Option::as_mut for apply/tracked_apply; return None/correspondence. Only finish uses take. Reject an unavailable editor before consuming Edits. |

Migrate Rust editor application tests in `ir/edit.rs`, `ir/molecule/tests.rs`,
`editor/transact.rs`, transaction tests, `tests/property/edit.rs`, and
`tests/property/molecule/compaction.rs`, plus affected resolver tests. Replace returned-editor destructuring/chaining with inspection of
the same editor. Test retained partial state on error; do not weaken publication
or old-state precondition assertions. Python `test_transaction.py` and
`test_edit.py`, and binding unit tests, must distinguish terminal finish from
nonterminal apply, including Edits-transfer failure. These changes cannot wait
for S11. In `test_molecule.py`, migrate the editor=True branch of
test_molecule_apply_consumed_edits here: a failed Edits transfer leaves the editor
available. The Molecule branch migrates in S0c2.

In `benches/editor.rs`, removal still measures application to an existing editor:
make the input editor mutable and return it after applying, rather than measuring
unit destruction. Update the editor paragraph of the Python API guide and the
changed methods' rustdoc in this same subitem.

**Checks:** direct-plus-batch editing, repeated batches, partial failure followed
by inspection/editing, batch correspondence, publication integrity, and Python
editor/Edits lifecycle. Migrate callers and tests with the signature change.

**Interface status:** the apply/tracked_apply receiver and return changes,
partial-state behavior, checked probe, and terminal finish are specified above.
No additional interface decision is identified for this subitem.

Completed 2026-10-14. Rust and Python editor application now mutates the existing
editor; failure retains completed edits and permits further mutation and batches.
All callers and affected lifecycle tests are migrated, with existing planning,
probe, and publication boundaries preserved. The removal benchmark still returns
the editor from the timed closure, preserving its result-destruction boundary.
The data-type, integrity, and Python API guides describe the implemented contract.

Verification passed: graph-IR molecule tests (1,245), edit tests (86), graph
operation tests (1,603), feature-gated edit properties (21) and compaction
properties (3), Python binding editor unit tests (6), and Python transaction,
edit, and molecule tests (331) against a freshly built extension in Python 3.13.
All-target Clippy for graph-IR, graph, and Python passed with graph-IR proptest
enabled and warnings denied, including compilation of the editor benchmark.
Changed Rust files pass nightly formatting checks; git diff --check passes.
Logs are retained in `scratch/181-s0b/`. Full workspace tests remain at the stage
gate.

#### S0c1 — Dropped

into_editor is not added; owned callers keep `edit(self)`. See [236](236-operation-ownership-copy-costs-2026-10-09.md).

#### S0c2 — Python Molecule editing without consumption (done)

**Modules:** Python Molecule edit/apply, binding tests, and the Python guide.
**Python behavior change; Rust signatures unchanged.** [dep: S0b]

Python edit, apply, and tracked_apply borrow the source through to_rust, clone
it, and call the consuming Rust method. They never take the molecule or advance
its counter. Transfer Edits and return the new editor, molecule, or
correspondence.

Python Molecule still has Option storage because resolve/resolve_with_report
consume until S10a. Keep to_rust fallible and retain take for those consumers;
S11a removes that storage only after the last consuming binding changes. Update
binding unit tests and edit/apply-triggered cases in `test_transaction.py`,
`test_edit.py`, `test_molecule.py`, `test_molecule_view.py`, `test_ring.py`,
`test_constraint_readings.py`, `test_molecule_constraint.py`, `test_depiction.py`,
`test_reaction.py`, and `test_resolve.py`. Cover retained nested views and molecules
obtained from Reaction and Solution. Preserve invalidation coverage under actual
whole-molecule mutations. Remove redundant copy-before-edit/apply in Python
`test_resolve.py` and `test_workflow.py`; retain the source-preservation assertions.
Keep resolution-triggered consumption tests until S10a; edit-triggered cases in
the resolver tests migrate here. Update the Python guide's Molecule editing
paragraph.

**Checks:** source preservation on success/failure, independent results, retained
source views, and correct tracked results. Rust edit/apply/tracked_apply,
resolver/transformer signatures, and ingestion/export calls remain unchanged.

Completed 2026-10-14. Python edit/apply/tracked_apply borrow the source and consume
a Rust copy, preserving source aliases, views, and the accessor counter. Edits
transfer and failure behavior remain explicit: an unavailable source does not
consume Edits; a transferred batch stays consumed on application or publication
failure. Tests cover independent results, retained entity and nested constraint
views, rings, and molecules held by Reaction and Solution. Consumption and
invalidation tests now use operations that still consume or mutate. The Python
guide and binding docstrings describe the implemented behavior.

Verification: 77 Rust molecule-binding tests passed. Across the eleven affected
Python test files, 1,064 cases passed after correcting the new test assertions;
two existing molecule-scope constraint-matching cases remain skipped for doc 195.
The corrected molecule/ring files passed their 208-case rerun. Tests used a fresh
extension in Python 3.13. All-target Clippy for umol-py with warnings denied,
nightly formatting of the changed Rust file, and git diff --check passed. Logs
are in `scratch/181-s0c2/`; full workspace tests remain at the stage gate.

#### S0d — Scoped in-place editing conveniences (done)

**Modules:** graph-IR molecule and tests. **Additive, green.**
[dep: S0b]

Implement edit_into and apply_into using the concrete signatures in the interface
review. One working copy, `self.clone().edit()`, one successful installation, no
journal. Return callback
errors unchanged; convert finish errors through MoleculeApplyError. An editor
application error does not latch an abort: the callback may catch it, continue
editing, and return Ok. Install that result if finish passes integrity checks.
A propagated callback error discards the working copy. Do not add Python callback
editing as part of this item.

**Checks:** return-value propagation, callback/application/integrity failure,
callback unwinding, source preservation, successful installation, and repeated
use. Include a caught application error followed by repair and successful
publication. Keep the no-recovery editor contract distinct from this receiver guarantee.

**Interface status:** edit_into/apply_into and the callback signature in the
migration interface review are settled. The signature follows the selected
clone-once, checked-installation and receiver-preservation semantics.

Completed 2026-10-14. Molecule::edit_into and apply_into use one working copy and
install only after successful publication. Callback errors propagate unchanged;
integrity errors convert through MoleculeApplyError. The tests cover returned
values, custom errors, propagated application errors, publication failure,
unwinding, repeated use, and repair after catching a partial application failure.
The public rustdoc and data-type guide state the structural-equality recovery
contract. Python bindings are unchanged.

Verification passed: 555 focused molecule tests, 24 feature-gated edit properties
(including three new scoped-editing properties), all-target graph-IR Clippy with
proptest enabled and warnings denied, nightly formatting of changed Rust files,
and git diff --check. Logs are in `scratch/181-s0d/`; full workspace tests remain
at the stage gate.

#### S0e — Transaction correspondence independent of Undo (done)

Completed candidate and measurement record. The unconditional tracking described
here is replaced by the 2026-10-14 decision scheduled in S0e1/S0e2; its independent
correspondence mechanism and measurements remain the starting point.

**Modules:** graph-IR molecule apply/transact and correspondence tests.
**Internal implementation change, green; tracked_commit's public result is unchanged.**
[dep: S0a]

Implement transaction-entry-to-result correspondence without retaining attribute
or constraint Undo payloads solely for tracking. Reuse ApplicationState's existing
batch mapping behavior. Preserve per-apply handle namespaces, cascading removals,
surviving additions, and the absence of participant-position mapping. Journaled
rollback continues to use Undo independently of this tracking mechanism.

**Storage and calls:** add correspondence: MoleculeCorrespondence to the private
guard and a mutable borrow of it to Transaction. Construct identity over all eight
entity kinds at run entry. After a successful apply loop, call the existing
state.update_correspondence on it and drop state. tracked_commit checks integrity
and moves out the accumulated correspondence; plain commit discards it. Remove
both Undo-inspection passes from tracked_commit. Failed applications abort the
transaction, so their partial mapping is never returned.

**Cost check:** use S0a's plain/tracked transaction fixtures to compare time,
allocation traffic and peak extra memory with the previous Journal implementation.
Cover empty and attribute-only batches as well as structural changes, small/larger
molecules, and unique/shared storage. Include identity construction and updates
in the measurement. Record the incremental cost on plain commit explicitly. Do
not add a mapping history or deferred construction without measured justification.
The same check for Clone belongs to S0f, where that mode becomes available.

**Checks:** entry-to-result mapping across additions, removals/cascades, deletion
of newly added entities, separate batches, attribute-only batches, empty
transactions, and rollback/abort. Measure retained mapping memory separately from
undo storage; changing batch partition must not change the resulting mapping.

Completed 2026-10-14. The private guard now owns the entry-to-result
MoleculeCorrespondence and Transaction borrows it. run constructs all eight
identity components; successful apply calls update them through ApplicationState.
tracked_commit checks integrity and moves out the mapping. Both journal-inspection
passes are removed. Undo still supplies Journal rollback; public signatures,
errors, acceptance and restoration behavior are unchanged. Rustdoc and the
development guides now describe independent correspondence accumulation.

Verification passed: 1,256 molecule unit tests, 25 feature-gated edit tests,
all-target graph-IR Clippy with proptest and warnings denied, nightly formatting
and git diff --check. Exact cases cover temporary additions removed within a
batch and empty batches, alongside existing all-kind identity, overlays, cascades,
cross-batch creation/removal and recovery cases. The new partition property varies
batch size over descending atom removals and subsequent additions, including
cascaded overlay removals. Full workspace checks remain at the stage gate.

The cost comparison rebuilds the preceding Journal implementation at
7d98315387d1b55c9e484bbd833cdcf766c9cbf4 and this implementation with the same
harness. It uses S0a's resolved chains and charge edits, extended with empty
transactions, empty batches and structural changes. Three sequential before/after
pairs alternate binary order; each row below is the median of three run medians,
each from 20 means of 32 calls. Input preparation and final molecule destruction
are excluded; transaction entry, mapping updates, integrity checking and recovery/
correspondence destruction are included. Allocation measurements use separate
instrumented binaries. All 15 fixtures pass output/correspondence checks; 60
case/ownership/tracking combinations are measured. Sources, binaries, commands,
per-run spread, allocation counts and retained memory probes are in
`scratch/181-s0e/`; these are bounded cases, not a workload mix.

Times in µs, before → after:

| Atoms / operation / ownership | Plain | Tracked |
| --- | --- | --- |
| 96 / no batches / unique | 0.333 → 0.488 (+46.5%) | 0.505 → 0.510 (+0.9%) |
| 8 / 1 charge edit / unique | 0.240 → 0.318 (+32.9%) | 0.305 → 0.342 (+12.1%) |
| 96 / 1 charge edit / unique | 0.520 → 0.690 (+32.8%) | 0.687 → 0.704 (+2.5%) |
| 384 / 1 charge edit / unique | 1.400 → 1.640 (+17.1%) | 1.630 → 1.670 (+2.4%) |
| 96 / 1 charge edit / shared | 1.859 → 1.978 (+6.4%) | 2.005 → 2.003 (−0.1%) |
| 96 / 8 charge edits / unique | 0.950 → 1.169 (+23.1%) | 1.120 → 1.191 (+6.3%) |
| 96 / 64 charge edits / unique | 4.805 → 5.003 (+4.1%) | 5.019 → 5.013 (−0.1%) |
| 384 / 64 charge edits / unique | 5.494 → 5.727 (+4.2%) | 5.783 → 5.783 (0.0%) |
| 384 / structural batch / unique | 72.639 → 73.765 (+1.6%) | 73.832 → 73.841 (0.0%) |

The structural batch removes an initial atom and its incident bond, adds an atom
and bond, then adds and removes a temporary atom. Across both ownership cases and
8/96/384 atoms, plain structural medians increase 1.6–3.0%; tracked medians change
−0.2–+1.9%. Small tracked attribute cases also have measurable overhead: the largest
relative increase is the 8-atom single edit, +37 ns (+12.1%). No frequency-weighted
or blanket performance-neutral conclusion follows from these measurements.

The identity map adds two allocations and 120/1,528/6,136 requested heap bytes
for chains of 8/96/384 atoms. Callback-entry live-memory probes isolate that storage
before any Undo exists. After 1/8/64 scalar changes, retained journal buffers remain
3,008/6,016/48,128 bytes, independently of molecule size. Plain commits now pay the
mapping cost even with no edits; an empty molecule still allocates nothing.
Tracked attribute-only allocation traffic is unchanged, but peak extra memory
increases because the mapping is live during application. For a unique 96-atom,
one-edit operation, plain calls/bytes/peak change from 1/3,008/3,008 to
3/4,536/4,536; tracked changes from 3/4,536/3,424 to 3/4,536/4,536. Structural
tracking adds one allocation and 28 requested bytes beyond the old tracked path,
while reusing the existing batch compaction logic.

This completed the eager tracking candidate and established its plain-commit
cost. The selected follow-up is separate Transaction and TrackedTransaction
handles in S0e1/S0e2. Clone-mode measurements remain part of S0f.

#### S0e1 — TrackedTransaction scoped handle (done)

**Modules:** graph-IR molecule/transact, molecule and ir exports, transaction
unit tests and edit properties. **Additive, green.** [dep: S0e]

Add TrackedTransaction with the settled run, apply, probe, commit and rollback
surface. At this boundary run takes the molecule and callback, matching the
currently implemented Journal-only Transaction::run. S0f adds RollbackMethod to
both entry points together. Export the new type beside Transaction through
molecule and ir. Use the existing errors and keep fields private; add no arbitrary
construction, conversion, mutable dereference to Transaction or Python handle.

TrackedTransaction::commit consumes the handle, checks integrity and returns the
entry-to-commit-state MoleculeCorrespondence. Acceptance still requires callback
success. run owns the restoration guard and returns the callback's T; the handle
only borrows scoped state. apply consumes one batch and returns unit; probe lends
the checked molecule; rollback consumes the handle. Retain every existing abort,
cancellation and restoration behavior. Reuse S0e's correspondence mechanism and
shared private recovery logic without constructing a second mapping.

This additive boundary leaves Transaction::tracked_commit and its existing
tracking behavior temporarily in place; S0e2 removes them and migrates their
callers. Document the new type's contract here without prematurely claiming that
ordinary transactions have no mapping cost.

**Checks:** exact correspondence for all eight entity kinds, additions,
removals/cascades, deletion of newly added entities, separate handle namespaces,
participant replacement and empty/attribute-only batches. Cover successful commit,
first/later application failure, failed commit, failed probe, explicit rollback,
uncommitted return, callback error after commit, unwinding between completed edits,
forgotten handles and probe after abort. Check the same errors and restoration as
Transaction. Extend public-operation properties for batch composition/partition,
rollback and equality of plain/tracked outcomes after discarding correspondence.

Completed 2026-10-14. TrackedTransaction exposes run, apply, probe, commit and
rollback, with private fields and exports through molecule and ir. Its commit
returns MoleculeCorrespondence and still requires callback success for acceptance.
The implementation privately reuses Transaction and its restoration guard; it
constructs no second mapping and exposes no conversion or mutable access to the
ordinary handle. Transaction::tracked_commit and ordinary eager tracking remain
until S0e2.

Exact correspondence cases now exercise the new handle, with lifecycle cases
covering cancellation, first/later application failure, probe repair, aborted and
failed commit, callback errors/unwinding after commit, and forgotten handles.
Properties check plain/tracked outcomes and errors, rollback, and composition and
partition independence, including structural edits and cascades. Verification
passed: 1,272 molecule unit tests, 28 feature-gated edit tests, all-target graph-IR
Clippy with proptest and warnings denied, nightly formatting and git diff --check.
Logs are in `scratch/181-s0e1/`. Cost comparisons remain in S0e2 and workspace
checks remain at the stage gate.

#### S0e2 — Ordinary transactions without correspondence (done)

**Modules:** graph-IR transaction implementation and callers, molecule
conveniences, transaction/property tests, development guides, and bounded recovery
measurements. **Breaking method removal, red→green within the subitem.**
[dep: S0e1]

Remove Transaction::tracked_commit. Remove aggregate correspondence storage,
identity construction and update_correspondence calls from ordinary Transaction
execution. Keep ApplicationState for per-batch handle resolution. Shared recovery
storage and acceptance logic remain independent of tracking; only the tracked
run owns a MoleculeCorrespondence and only tracked apply updates it. Neither
public type exposes mutation through the other. Do not replace the removed
method with a runtime tracking-disabled error or optional commit result.

Migrate Molecule::tracked_transact to TrackedTransaction::run and commit;
Molecule::transact stays on Transaction::run and commit. Their public signatures
and Python results remain unchanged here. Migrate direct tracked_commit callers
in transaction tests and examples at the same boundary. Preserve their callback
error/unwind-after-commit and aborted-commit scenarios with the new tracked type;
do not reduce them to ordinary commit cases. Composite Resolver continues to use
ordinary Transaction. Audit workspace consumers, feature-gated properties and
benchmarks for the retired method.

Update rustdoc and data-types.md, integrity.md and nomenclature.md to describe
the two operation-issued handles and their commit results. Remove tracked_commit
from current API examples and publication tables. Python's prepared-batch
operations continue to delegate to the corresponding Rust convenience; scoped
callback handles remain Rust-only. This is the documentation boundary for the
split, rather than waiting for S11b.

**Checks and cost gate:** retain S0e1's lifecycle and correspondence coverage on
both types, including plain/tracked result/error parity and batch partition
independence. Run focused transaction/molecule and feature-gated edit tests.
Repeat the S0e empty, attribute and structural cases across small/larger molecules
and unique/shared storage. Compare the split with both saved pre-S0e Journal and
S0e eager-tracking binaries; record timing, allocation traffic, peak extra memory
and retained correspondence separately from Undo. Verify by implementation
inspection and allocation probes that ordinary transactions neither construct
nor retain aggregate mapping storage. Include entry/update costs on the tracked
path and check its output before timing. Report residual differences from the
original baseline rather than assuming the split removes all overhead. Record
the result before S0f; retain the saved measurements and harness boundaries.

Completed 2026-10-14. Transaction and its guard contain no correspondence storage.
Only TrackedTransaction::run constructs the identity mapping, only its apply
updates it, and its commit returns it after the shared integrity check. Both
handles initially shared a private batch executor, forced inline into both apply
methods. The later direct-loop revision below removes that helper and all related
inlining attributes while retaining shared recovery logic. Transaction::tracked_commit
is removed; Molecule::tracked_transact now
uses the tracked handle. Composite Resolver remains on Transaction. The Rust and
Python prepared-batch signatures and results are unchanged.

Current source and guide references to tracked_commit are removed. The development
guides describe both operation-issued handles and their distinct commit results.
The four former tracked_commit lifecycle cases already have equivalents on
TrackedTransaction from S0e1; those error, unwind and abort cases remain covered.
Verification passed: 1,268 molecule unit tests, 28 feature-gated edit tests,
all-target graph-IR Clippy with proptest and warnings denied, nightly formatting
and git diff --check. Release benchmark builds also compile the graph consumers;
all benchmark molecule/correspondence checks pass. Workspace checks remain at the
stage gate.

The cost comparison uses the saved original and eager S0e binaries, this split,
and an initial split without forced inlining. It retains S0e's 15 fixtures,
60 case/ownership/tracking combinations and measurement boundaries. All four
versions run sequentially in four rotating orders, after builds and tests finish;
each reported value is the median of four run medians, each from 20 means of 32
calls. The candidate uses the same Rust toolchain and dependency lockfile as S0e.
Sources, binaries, commands, raw spreads and allocation probes are preserved in
`scratch/181-s0e2/`, with the original saved binaries in `scratch/181-s0e/`.

Times in µs; original means the pre-S0e journal-derived implementation:

| Atoms / operation / ownership / mode | Original | Eager | Split | Split vs original |
| --- | ---: | ---: | ---: | ---: |
| 96 / no batches / unique / plain | 0.329 | 0.498 | 0.326 | −1.2% |
| 8 / 1 charge edit / unique / plain | 0.244 | 0.309 | 0.227 | −7.1% |
| 96 / 1 charge edit / unique / plain | 0.533 | 0.694 | 0.523 | −1.8% |
| 384 / 1 charge edit / unique / plain | 1.397 | 1.654 | 1.394 | −0.2% |
| 96 / 1 charge edit / shared / plain | 1.835 | 1.935 | 1.842 | +0.4% |
| 96 / 1 charge edit / unique / tracked | 0.712 | 0.717 | 0.718 | +0.9% |
| 8 / structural batch / unique / tracked | 5.422 | 5.468 | 5.562 | +2.6% |
| 384 / structural batch / unique / plain | 73.506 | 73.695 | 74.062 | +0.8% |
| 384 / structural batch / unique / tracked | 75.358 | 74.050 | 76.134 | +1.0% |

Across edit-bearing cases, the largest median increases over the original are
2.0% plain and 2.6% tracked. Some empty tracked cases retain 6–8 ns of overhead:
up to 11.5% for a zero-atom transaction and 5.5% for an empty batch on a shared
8-atom chain. Larger samples can have broad spreads: the unique 384-atom,
64-edit plain case spans 5.536–7.928 µs originally and 5.558–6.087 µs after the
split. These bounded cases do not establish a workload-weighted or uniformly
cost-neutral result.

Performance acceptance remains open. Continuing implementation does not accept
the measured 2.0% plain or 2.6% tracked increases as negligible. Their causes have
not been isolated, and subsequent comparisons must retain the original baseline.

Inlining removes an avoidable call in the initial split. The unique 8-atom
tracked empty-batch median changes from 199.5 ns (195.2–203.2 across run medians)
to 178.1 ns (175.2–182.9), versus 177.7 ns originally. The retained implementation
shares the source loop without paying that additional call; larger timing
differences remain mixed.

All 30 plain allocation rows exactly match the original's allocation count,
requested bytes and peak extra live heap. All 30 tracked rows match the eager
candidate on those metrics. Callback-entry probes confirm zero ordinary mapping
storage, while tracked chains of 8/96/384 atoms retain 120/1,528/6,136 bytes.
After 1/8/64 scalar edits, both handles retain journal buffers of
3,008/6,016/48,128 bytes, independently of molecule size.

Tracked peak-memory costs relative to the original remain: the unique 96-atom,
one-edit case still allocates 3 times and requests 4,536 bytes, but peaks at
4,536 extra bytes instead of 3,424 because mapping and Undo coexist. The same
plain case returns to 1 allocation, 3,008 requested bytes and 3,008 peak extra
bytes. Structural tracking retains S0e's additional allocation and 28 requested
bytes relative to the original. No mapping history or Undo-based reconstruction
is introduced to remove those tracked costs.

#### S0f — Clone and Journal transaction recovery (done)

**Modules:** graph-IR Transaction and TrackedTransaction recovery, both run
entry points' Rust callers, and tests.
**Breaking run signatures, red→green within the subitem.** [dep: S0e2]

Add RollbackMethod::{Clone, Journal} and the agreed selector position on both
run methods: receiver, rollback method, callback. Clone applies to one working
molecule and installs only after successful commit plus callback return. Journal
retains undo-based restoration. Keep each type's apply/probe/commit/rollback signatures,
separate batch namespaces and abort behavior. TrackedTransaction accumulates
correspondence with either recovery method; Transaction does not.

Migrate run callers together. Existing transact/tracked_transact conveniences
retain their signatures at this boundary and explicitly pass Journal internally;
S0g exposes the selector on those conveniences. Composite Resolver callers also
pass Journal until S10a exposes the operation's recovery selector.

**run caller migration in this subitem:** Molecule transact/tracked_transact in
`ir/molecule/transact.rs` pass Journal to Transaction::run and
TrackedTransaction::run respectively. Composite Resolver resolve_into/
resolve_into_with_report/project_into pass Journal to ordinary Transaction::run.
Their public signatures stay unchanged here. Update graph-IR molecule/transaction
unit tests, edit properties and both run methods' rustdoc examples at the same
boundary. Exercise both recovery methods on both handle types; retain Journal
explicitly where a test specifically exercises undo restoration. No run caller
waits for S10.

**Guard storage:** retain the receiver borrow and status on the shared
TransactionGuard; replace its journal field with the approved TransactionStorage
enum above. Each handle borrows recovery state; the tracked handle also borrows
the run's optional correspondence storage. Both run methods initialize Journal
with an empty Vec or Clone with Some(receiver.clone()). Keep restoration owned
by run; no empty placeholder is installed in the receiver. Application dispatch
selects apply_edit_with_undo for Journal and apply_edit for Clone. Both feed the same
per-batch ApplicationState; only TrackedTransaction uses it to update the
accumulated correspondence.

Probe selects the working copy while Clone is active and the receiver after
abort; Journal always probes the receiver. Commit checks that same active
molecule. Accept Clone only after commit plus callback success, moving its value
from the Option into the receiver. Rollback and guard Drop replay Journal undo
or discard the Clone value. Retain the existing status transitions and scoped
callback lifetime. Repeat S0e2's cost cases for both Clone and Journal, keeping
ordinary and tracked execution separate and including recovery setup/disposal.
Export RollbackMethod beside Transaction and TrackedTransaction through molecule
and ir, following the existing transaction-type exports.

**Checks:** both recovery modes and both handle types on success, first/later
application failure, failed commit, explicit rollback, callback error after commit,
unwinding between completed edits, forgotten handles, and tracked structural
changes. Probe after abort observes the
restored entry value. Extend composition, partition and plain/tracked parity
properties across both recovery modes. Do not add recovery from allocation failure
or internal partial-Edit panics.

**Outcome — 2026-10-14.** RollbackMethod::{Clone, Journal} is exported beside
both scoped handles. Their run methods take receiver, recovery method and callback;
construction and tracking remain separate from recovery selection. Clone edits one
working molecule through apply_edit without Undo, and installs it only after
commit plus callback success. Journal keeps undo restoration. Recovery selection
occurs once per batch, and only TrackedTransaction updates correspondence.

Both methods retain the run-owned guard, abort latch, checked probe and commit,
separate batch namespaces, callback-error/unwind restoration and forgotten-handle
recovery. Every Rust run caller is migrated. Molecule prepared-batch conveniences
and composite Resolver callers select Journal internally; their signatures remain
unchanged until S0g and S10a respectively.

Verification passed: 1,314 molecule unit cases, 28 feature-gated edit tests,
643 resolver unit cases, all-target graph-IR/graph Clippy with graph-IR proptest
and warnings denied, nightly formatting and git diff --check. Lifecycle cases
exercise both handles and both methods; generated tests cover restoration,
plain/tracked parity, correspondence composition and batch partitioning.
No workspace or Python gate was run for this subitem.

The cost comparison retains the original and S0e2 binaries, the initial S0f
Journal implementation, and final Journal/Clone execution. It repeats the same
15 fixtures and 60 cases per method, with five sequential rotating orders after
builds/tests finish. Values are medians of five run medians, each from 20 means
of 32 calls. The dependency lockfile matches S0e2. The candidate uses scoped
callbacks with a runtime selector; the saved baselines use prepared-batch methods
with the same apply/commit sequence. Both include recovery setup/disposal.
Sources, binaries, disassembly, raw spreads, allocation probes and commands are in
`scratch/181-s0f/`; the initial candidate and measurements are retained in its
`initial/` directory.

The initial iterator-based Journal loop introduced an out-of-line IntoIter::try_fold
call. Direct loops remove that call from both apply methods while keeping recovery
dispatch outside the loop. In the final comparison, the unique 96-atom, eight-edit
case changes from 1.087 to 1.036 µs plain and 1.289 to 1.220 µs tracked; S0e2 is
1.038 and 1.219 µs respectively. Other differences remain and are not accepted
merely because this avoidable cost was removed.

Times in µs; original means the pre-S0e journal-derived implementation:

| Atoms / operation / ownership / mode | Original | S0e2 | Journal | Clone | Journal vs original |
| --- | ---: | ---: | ---: | ---: | ---: |
| 96 / 1 charge edit / unique / plain | 0.524 | 0.523 | 0.524 | 2.444 | +0.0% |
| 96 / 1 charge edit / unique / tracked | 0.700 | 0.741 | 0.773 | 2.659 | +10.5% |
| 96 / 1 charge edit / shared / plain | 1.817 | 1.822 | 1.830 | 1.778 | +0.7% |
| 96 / 8 charge edits / unique / tracked | 1.210 | 1.219 | 1.220 | 2.824 | +0.9% |
| 96 / 64 charge edits / unique / plain | 4.697 | 4.774 | 4.431 | 3.387 | -5.7% |
| 8 / structural batch / unique / tracked | 5.440 | 5.473 | 5.528 | 4.239 | +1.6% |
| 96 / structural batch / unique / plain | 21.812 | 21.994 | 21.552 | 17.962 | -1.2% |
| 384 / 1 charge edit / unique / plain | 1.428 | 1.409 | 1.484 | 9.339 | +3.9% |
| 384 / 64 charge edits / unique / plain | 5.524 | 6.564 | 6.108 | 12.596 | +10.6% |
| 384 / 64 charge edits / unique / tracked | 5.869 | 6.442 | 7.127 | 12.925 | +21.4% |

The large unique 384-atom, 64-edit case is unstable across runs. Its tracked
Journal median is 21.4% above the original in this comparison, with trial medians
ranging from 5.468 to 7.658 µs, versus 5.680–7.175 µs originally and
5.487–8.209 µs for S0e2. The plain medians rise 10.6%, with similarly broad
ranges. These observations are retained; they do not establish a stable slowdown
bound or cost-neutrality.

A separate entry-point check in `scratch/181-s0f/boundary/` compares original,
S0e2, final scoped Journal and final prepared-batch Journal in four rotating
orders. The latter uses the same Molecule convenience methods as the frozen
baselines. Among edit-bearing cases, its largest increases over the original
are 3.1% plain (231.1 → 238.25 ns) and 10.1% tracked (294.9 → 324.55 ns), both
for a unique eight-atom, one-edit operation. S0e2 is 231.8 and 320.3 ns in those
cases. The largest prepared-batch increases over S0e2 occur in different cases:
2.9% plain for 64 edits on a shared 384-atom molecule, and 4.9% tracked for one
edit on a shared eight-atom molecule. The unique 384-atom, 64-edit tracked medians
are now 5.871 µs original, 5.584 µs S0e2, 5.483 µs scoped and 5.720 µs prepared.
This check does not isolate a stable entry-point penalty or explain the changing
large-case results. The earlier S0e2 2.0%/2.6% figures are observations from that
comparison, not an established upper bound. Empty transactions also retain overhead:
in the same-entry-point check, a shared zero-atom input changes from 39.1 to
43.3 ns plain (+10.7%) and 67.8 to 85.3 ns tracked (+25.8%); S0e2 is 40.35
and 82.0 ns respectively. Performance acceptance remains open.

Journal matches S0e2's allocation count, requested bytes and peak extra heap in
all 60 cases. Ordinary entry mapping storage remains zero; tracked chains of
8/96/384 atoms retain 120/1,528/6,136 bytes with either recovery method.
Retained scalar-edit Journal buffers remain 3,008/6,016/48,128 bytes for 1/8/64
edits. Clone records no Undo; the first scalar edit detaches 1,576/18,472/73,768
bytes for 8/96/384 atoms, with no further retained copy storage for subsequent
scalar edits. Clone's working copy itself adds no entry heap allocation in these
fixtures because molecule storage is shared until mutation.

For a unique 96-atom molecule, a single plain edit therefore peaks at 18,472
extra bytes with Clone versus 3,008 with Journal; at 64 edits Journal peaks at
48,128 while Clone stays at 18,472. Shared input already requires detachment with
Journal. The unique 96-atom structural batch uses fewer allocations with Clone
(102 versus 114), but higher peak extra heap (65,600 versus 30,584 bytes).
The timings and heap costs support keeping recovery selection explicit; they do
not establish one method as uniformly preferable or accept remaining regressions.

**Review correction — 2026-10-14.** S0e2 and S0f are not complete while their
performance gate remains unresolved. S0g must not advance on the basis of passing
correctness tests. The original pre-S0e Journal implementation remains the
acceptance baseline for the combined change; intermediate versions are diagnostic
comparisons, not progressively relaxed acceptance baselines.

The prepared-batch methods already called Transaction::run before S0e. Changing
the benchmark to call scoped run directly introduced no new transaction semantics
and was unnecessary for the Journal comparison. Any resulting execution overhead
is an implementation or measurement issue, not a semantic justification. S0f
adds selectable Clone recovery, but Journal performs the same edits, integrity
checks and recovery obligations as before. Independent tracked correspondence
supports execution without Undo; it computes the same mapping, and its current
identity construction and per-batch updates are implementation choices whose
cost still needs justification or removal.

The 4.9% tracked increase over S0e2 in the same-entry-point comparison is
399.1 → 418.65 ns, or 19.55 ns. The 2.9% plain increase is
11,166.65 → 11,492.85 ns, or 326.2 ns, in a different case. The cumulative
unique eight-atom tracked increase over the original is 294.9 → 324.55 ns,
or 29.65 ns. These remain unresolved observations; neither call-site terminology,
unchanged allocation counts nor broad timing spreads establishes acceptability.

The corrective work is to hold the entry point and measurement boundaries fixed,
check repeatability with unchanged-binary controls, then isolate recovery dispatch,
acceptance/destruction, and correspondence maintenance. Empty and attribute-only
batches deserve explicit attention: their correspondence is unchanged, while
update_correspondence currently checks all eight removal tables and invokes all
eight right-domain extensions even when every added count is zero. Whether the
compiler removes those operations and how much they cost must be established,
not assumed. Changes must retain all lifecycle, error, mapping and rollback laws.
The gate covers cumulative cost against the original, including empty operations,
shared/unique storage, allocation traffic and peak memory. Inconclusive timing
keeps the gate open; unexplained reproducible regressions do not authorize S0g.

A fallback proposal, if the regressions cannot be removed within this transaction
work, is to restore the pre-S0e Journal execution and correspondence construction,
and reconcile the dependent scoped-handle/caller/test/document changes together.
This would roll back the costly transaction changes rather than only remove S0f
and retain unexplained S0e/S0e2 costs. Earlier lifecycle/editor work would remain.
This paragraph records a fallback plan; no production rollback is performed here.

**Cost investigation — 2026-10-14 (inlining subsequently removed).** The measured candidate removes concrete
implementation overhead without changing transaction semantics or public symbols:

- Tracked commit previously copied out its 320-byte correspondence, then wrote
  another complete empty correspondence into the run-owned slot. This candidate
  used an Option, consumed with take. Generated code replaced the empty-value
  rewrite with a discriminant store; run-owned cleanup covered forgotten handles.
- Attribute-only and empty batches previously entered the full correspondence
  updater, including its register saves, stack frame and zero-count domain
  updates. An inline check of the existing handle tables now skips that routine
  when no ids changed. No extra change flag or mapping allocation is introduced.
- Both apply and commit methods permit inlining. The former out-of-line
  boundaries prevented specialization of the known Journal selection in the
  prepared-batch callers. In the final timing binary, those method calls are
  inlined. Clone and Journal still share the same public handles and lifecycle.

The study retains the original prepared-batch entry points, inputs, timing
boundaries, compiler profile and pre-S0e binary. Attribution candidates were
measured separately; the final candidate uses eight rotating trials of original,
unchanged control, candidate and the same unchanged control again. No builds or
tests overlap timing. All 60 cases remain included. Selected medians, in ns:

| Operation / atoms / ownership / tracking | Pre-S0e | Before these fixes | Final candidate |
| --- | ---: | ---: | ---: |
| Empty / 0 / unique / plain | 39.05 | 41.70 | 39.45 |
| One attribute edit / 8 / unique / plain | 234.00 | 237.60 | 230.45 |
| One attribute edit / 8 / unique / tracked | 295.30 | 322.90 | 286.45 |
| One attribute edit / 8 / shared / tracked | 402.65 | 418.90 | 394.20 |
| 64 attribute edits / 384 / shared / plain | 11,904.00 | 11,424.45 | 11,541.00 |
| Structural / 8 / unique / plain | 5,313.50 | 5,255.60 | 5,415.00 |
| Structural / 8 / unique / tracked | 5,480.20 | 5,484.75 | 5,550.15 |

The before-fix column pools both unchanged-control labels. The small tracked
regression is removed in this comparison. This does not establish a general
performance bound: the largest remaining positive median is the small unique
plain structural case, +101.5 ns (+1.91%); its tracked counterpart is +69.95 ns
(+1.28%). These remain unresolved, not accepted as negligible. The intermediate
wrapper-only candidate also showed a small plain attribute regression and was
not used to claim success. Its data remains alongside the final candidate.

Unchanged-binary controls establish a measurement limitation independently of
implementation differences. In the wrapper comparison, the same control's two
eight-run medians for the large shared plain case were 10,690.75 and 11,179.05 ns
(+488.3 ns, +4.57%). This demonstrates that an apparent difference of that size
can arise without a code change; it does not establish that every observed
regression is noise. The final comparison includes trial ranges and paired
bootstrap intervals, not a claim of equivalence from overlapping ranges.

The fixes preserve the before-fix allocation traffic and peak memory in all
60 cases. Relative to pre-S0e, tracked structural cases still allocate one extra
28-byte temporary, and accumulated mapping storage still overlaps Undo storage.
Those costs are not removed by the timing optimizations. The cost gate remains
open for the residual structural difference and the cumulative resource trade-off;
S0g has not advanced. Evidence, source snapshots, binary hashes, all intermediate
results and the reproducible harness are in `scratch/181-transaction-cost-audit/`.

Verification of that candidate: all 1,314 molecule tests and 28
feature-gated edit properties pass, as does strict graph-IR all-targets Clippy.
Formatting and diff checks pass. Each timing binary checks final molecule and
correspondence results before measurement; the separate final allocation probe
matches the earlier candidate in all 60 cases.

**Direct correspondence ownership — 2026-10-14 (experiment withdrawn).** The
approved ownership proposal puts MoleculeCorrespondence directly in the tracked
handle and keeps recovery in the run-owned guard. Deliberate forgetting may leak
the mapping; preventing that leak is not required. This permits consuming commit
to move the correspondence without a borrowed Option or replacement value.

The straightforward implementation adds a 344-byte handle copy before commit
validation and a redundant 280-byte mapping copy during result propagation in
the generated tracked run. Its stack frame grows by 304 bytes. The small unique
tracked attribute case rises from 289.00 to 319.70 ns (+30.70 ns, +10.62%) against
the preceding implementation. These costs add no semantic work. Explicit field
destructuring reduces the first copy to 320 bytes but leaves the other copy;
that experiment was discarded after assembly inspection.

A second candidate keeps direct ownership and forwards the callback Result
through run rather than extracting its value and wrapping it again. This removes
the redundant mapping round trip in the generated code, but leaves the handle
copy. It does not recover the tracked attribute regression and also regresses
plain attribute cases. Selected unique-storage medians from its separate eight
rotating trials, in ns:

| Operation / atoms / tracking | Pre-S0e | Borrowed slot | Owned, result forwarded |
| --- | ---: | ---: | ---: |
| Empty / 0 / tracked | 68.90 | 67.70 | 73.45 |
| Empty batch / 8 / tracked | 168.60 | 176.40 | 189.10 |
| One attribute edit / 8 / plain | 235.40 | 227.20 | 244.75 |
| One attribute edit / 8 / tracked | 301.10 | 286.75 | 318.30 |
| Structural / 8 / tracked | 5,486.40 | 5,512.10 | 5,541.05 |

Both timed candidates pass the unchanged 1,314 molecule unit tests, 28
feature-gated edit properties and strict graph-IR all-targets Clippy. Allocation
traffic and peak memory match the preceding implementation in all 60 cases.
Both comparisons use the original pre-S0e binary and repeat the same preceding
binary under two control labels; the small attribute regressions exceed that
control variation. No public signatures, inlining attributes, mapping algorithms,
workloads or compiler options change. Individual copy costs have not been timed
separately; their presence does not establish a complete timing attribution.

Both production experiments were withdrawn, restoring the preceding borrowed
slot implementation exactly. Direct ownership is not a required implementation
target. S0e2/S0f remain open for correspondence
transfer costs and the existing cumulative cost questions before S0g. Sources,
binaries, assembly, samples and allocation comparisons are preserved in
`scratch/181-transaction-owned-map/`; the pre-S0e baseline is retained.

**Lazy correspondence — 2026-10-14 (experiment withdrawn).** This candidate
keeps the small borrowed handle and leaves identity implicit until the first
successful structural batch or successful commit. The first structural batch
uses its entry counts; prior attribute or participant changes preserve entity ids.
An identity-only commit validates before constructing the result. Public
signatures, recovery semantics, ordinary transaction code, inlining attributes,
mapping results, compiler options and benchmark workloads remain unchanged.

The resource improvement is concrete. No mapping heap is retained during
attribute-only work, compared with 120, 1,528 and 6,136 bytes for the current
8-, 96- and 384-atom cases. Peak extra requested heap matches pre-S0e in all
60 successful-operation cases. Allocation traffic matches the current control;
tracked structural cases still allocate one extra 28-byte temporary relative
to pre-S0e. The returned identity still allocates at commit and can then overlap
the journal; lazy construction does not eliminate the result's storage.

Two separate eight-trial comparisons retain the original prepared-batch entry
points and rotate original, current control, candidate and the same control
again. Builds and tests finish before timing. Selected unique-storage medians,
in ns; the current column pools both control labels:

| Comparison / operation / atoms / tracking | Pre-S0e | Current | Lazy | Lazy vs current |
| --- | ---: | ---: | ---: | ---: |
| First / empty batch / 8 / tracked | 173.15 | 169.25 | 180.05 | +6.38% |
| Repeat / empty batch / 8 / tracked | 169.25 | 163.65 | 177.75 | +8.62% |
| First / one attribute edit / 8 / tracked | 305.65 | 289.75 | 310.30 | +7.09% |
| Repeat / one attribute edit / 8 / tracked | 295.85 | 289.05 | 299.10 | +3.48% |

The repeated empty-batch cost is also 8.50 ns (+5.02%) above pre-S0e. Larger
cases remain variable: the repeated large shared plain 64-edit case is +4.80%
against pre-S0e, while the identical current control's two medians differ by
5.15%. These observations do not establish a general regression bound or justify
discarding the small-case increases. The heap improvement does not settle the
timing trade-off. No timing cost is justified by a new semantic requirement.

The candidate passes 1,318 molecule unit cases, 28 feature-gated edit properties
and strict all-target graph-IR Clippy. Additional cases cover an attribute-only
prefix before structural changes and failed commit while identity is implicit,
with both recovery modes. The prototype and its test/guide changes were withdrawn,
restoring their pre-experiment contents exactly. Source snapshots, patch, binary
hashes, allocation evidence, both timing series and assembly are retained in
`scratch/181-transaction-lazy-map/`. S0e2/S0f remain open; S0g has not advanced.

**Direct batch loops — 2026-10-14 (retained).** Remove all six inlining
attributes: both apply methods, both commit methods, the shared apply_batch
helper, and ApplicationState::update_correspondence. Remove apply_batch and put
its recovery dispatch, edit loops and abort handling directly in Transaction::apply
and TrackedTransaction::apply. Recovery storage, rollback, acceptance and commit
validation remain shared. No public interface, semantic contract, mapping algorithm
or test changes. The explicit loop duplication replaces reliance on forced inlining;
the compiler remains free to inline by its own heuristics.

The comparison keeps all 60 prepared-batch cases, original timing boundaries and
pre-S0e baseline. Eight rotating trials include the preceding inlining version
twice as an unchanged-binary control. Builds and tests finish before timing.
Selected medians, in ns; the preceding column pools both control labels:

| Operation / atoms / ownership / tracking | Pre-S0e | Preceding inlining version | Direct loops | Direct vs pre-S0e |
| --- | ---: | ---: | ---: | ---: |
| Empty / 0 / shared / tracked | 69.00 | 70.30 | 76.50 | +10.87% |
| Empty / 8 / unique / tracked | 145.50 | 137.35 | 156.85 | +7.80% |
| Empty batch / 8 / unique / plain | 98.00 | 93.70 | 106.80 | +8.98% |
| One attribute edit / 8 / unique / plain | 237.60 | 234.70 | 237.95 | +0.15% |
| One attribute edit / 8 / unique / tracked | 303.45 | 291.00 | 302.70 | -0.25% |
| One attribute edit / 96 / shared / plain | 1,816.10 | 1,824.50 | 1,841.80 | +1.42% |
| Structural / 8 / unique / tracked | 5,428.00 | 5,523.75 | 5,548.50 | +2.22% |

The small tracked structural increase is 120.50 ns over pre-S0e. Empty-operation
increases also remain part of the comparison; they are not excluded to present
only edit-bearing costs. The measurements cover the combined attribute removal
and loop duplication, not the individual contribution of either change.
Allocation traffic and peak extra heap match the preceding version in all 60
cases; the previously recorded extra structural allocation and retained mapping
costs remain. This implements the requested removal of inlining hints without
claiming the cumulative performance gate is satisfied. No further tuning
variants were tried; S0e2/S0f remain open before S0g.

Verification passed: 1,314 unchanged molecule unit cases, 28 feature-gated edit
properties, strict all-target graph-IR Clippy, nightly graph-IR formatting and
git diff --check. Sources, patch, binary hashes, allocation comparison, raw timing
series and uncertainty estimates are in `scratch/181-transaction-direct-loops/`.

**Decision to proceed — 2026-10-14.** Retain the direct-loop implementation and
continue with S0g despite the recorded regressions. S0e2 and S0f are complete for
this plan; further transaction tuning no longer blocks the remaining work.
This is acceptance of the known cost for proceeding, not a finding that the
performance differences were eliminated or explained. The original baseline
and all measurements remain available.

#### S0g — Prepared-batch transaction entry points (done)

**Modules:** graph-IR Molecule transact/tracked_transact; Python molecule,
rollback-selector binding/exports, callers, and tests.
**Breaking selector interfaces, red→green within the subitem.** [dep: S0f]

Expose RollbackMethod through Molecule transact/tracked_transact and Python
prepared-batch calls. Delegate to Transaction::run, apply each batch immediately,
and commit once for transact; tracked_transact delegates to TrackedTransaction::run
and its commit. Preserve the existing whole-operation correspondence result.
Python continues to accept arbitrary iterables, extract the batch objects before
receiver mutation, and transfer Edits sequentially without prechecking the entire
list for consumed or repeated objects. Failure transferring a later batch may
leave earlier batches consumed, but leaves the molecule and its view counter
unchanged. Before execution, advance the molecule counter under the existing rule.

**Interface status:** Rust takes batches followed by rollback_method. Python
exposes RollbackMethod and a keyword-only rollback_method defaulting to Journal;
Clone remains selectable. These signatures and the default are settled. Scoped
Python Transaction or TrackedTransaction handles/callbacks are not introduced.

**Prepared-batch caller migration in this subitem:**

| Caller | Required change |
| --- | --- |
| Python `src/molecule.rs` transact/tracked_transact | Accept the keyword-only RollbackMethod selector and forward it, default Journal; retain current extraction, sequential Edits transfer and counter timing |
| graph-IR editor/reaction benchmarks and molecule/transaction/property tests; graph stale-plan rollback tests | Supply the explicit selector. Keep baseline journal rows on Journal; add Clone comparison cases without relabeling old results. |

Update Python RollbackMethod exports with its binding and both method signatures.
Extend prepared-batch Python tests for explicit Clone/Journal; calls omitting the
keyword continue to exercise the selected Journal default. Update transaction
rustdoc and the Python guide here; S11b is the final consistency pass, not the
first documentation update.

**Checks:** Clone/Journal success and failure parity, plain/tracked calls,
arbitrary iterables, repeated/consumed batch inputs, independent batch namespaces,
source restoration, and unchanged view invalidation timing. Update the operation
docs with the selected Python default and exact receiver/input guarantees.

**Completed — 2026-10-14.** Both Rust conveniences now require RollbackMethod.
Python exports the selector and accepts keyword-only rollback_method with Journal
as its default. Batch extraction, sequential consumption, counter timing and
plain/tracked results are unchanged. Callers and guides are updated; benchmark
Journal rows retain their names and Clone comparison rows are added. Existing
prepared-batch cases exercise both recovery methods, including correspondence,
application and integrity failures, and Python argument/consumption failures.
Six successful-plan applications in graph tests use consuming application;
transactions remain in the three stale-plan tests that explicitly check rollback.

**S0 verification passed:** workspace all-target compilation with graph-IR and
graph property features; 1,320 molecule unit cases, 28 edit properties, three
compaction properties, 1,603 graph operation cases, 11 isotope properties,
1,640 Rust binding cases, and 1,096 Python lifecycle/API cases. Two existing Rust
binding cases remain ignored. Strict Clippy passed for all targets of graph-IR,
graph and Python with the property features. The extension was rebuilt in the
Python 3.13 environment. Nightly formatting checks and git diff --check passed.
Logs are in `scratch/181-s0g/`; no new timing comparison was run for this selector
plumbing.

**Compilation and test boundaries:** S0b, S0c2, S0e1, S0e2, S0f and S0g include their
binding/caller changes in the same subitem. Each ends with affected libraries,
tests and benches compiling and focused tests passing; no compatibility wrapper
or disabled test bridges the interval to S10. Check feature-gated property callers
explicitly. At the S0 stage gate, check workspace targets and run the affected
publication, edit, transaction, resolver/transformer and Python lifecycle suites.
Use the Python 3.13 environment for every PyO3 build. Full workspace test/lint and
MSRV runs remain at S12b; baseline timings are not rerun after every subitem.

**S0 gate:** all listed interface decisions are resolved where marked, affected
Rust/Python callers compile and pass their focused tests, and lifecycle/recovery
contracts are verified before the boxing and batching stages. Transaction exposes
no tracked_commit and performs no aggregate correspondence work; tracked results
come from TrackedTransaction::commit. The S0e2 and S0f cost comparisons are recorded.
Full workspace and MSRV verification remain at the final combined gate.

### S1 — Box the selected large payloads

#### S1a — Undo: 88-byte layout (completed)

**Modules:** graph-IR edit vocabulary, application/undo replay, transaction
correspondence handling and affected consumers.
**Breaking Rust payload shape, red→green within the subitem.** [dep: S0d, S0g]

Box the complete payload once for each of:

- RestoreRemovedTopology;
- RestoreRemovedDativeBonds, RestoreRemovedAromaticSystems,
  RestoreRemovedMulticenterBonds, RestoreRemovedNoncovalentBonds,
  RestoreRemovedStereoAtoms, RestoreRemovedStereoBonds;
- RemoveAddedAromaticSystem, RemoveAddedMulticenterBond,
  RemoveAddedStereoAtom, RemoveAddedStereoBond.

For the last four, box their existing Added* value. The seven removal payloads
retain their removed entries, compactions, and cascaded constraints together in
one box per Undo. Field, participant-list, and constraint undos remain inline.
Change construction and destructuring together; do not add clones or reconstruct
payloads merely to access them. Undo replay, id_compaction, and tracked commit
retain their existing meaning.

**Production declarations:** retain the study's complete
payloads, replacing its Experiment-prefixed types with the following open carriers
in ir/edit. All fields are public; derive Clone, Debug, PartialEq, Eq and export
through ir alongside the existing Added*/Removed* records. They introduce no
constructors, validation, or independent mutation operations.

```rust
pub struct RemovedTopology {
    pub atoms: Vec<RemovedAtom>,
    pub bonds: Vec<RemovedBond>,
    pub overlays: RemovedOverlays,
    pub compaction: MoleculeCompaction,
    pub undo_compaction: UndoCompaction,
    pub cascade: CascadedConstraints,
}
// Undo variant:
RestoreRemovedTopology(Box<RemovedTopology>),
```

Each row below defines a public struct with exactly three fields:
`removed: Vec<Entry>`, `undo_compaction: UndoCompaction`, and
`cascade: CascadedConstraints`. Its Undo variant takes one Box of that struct.

| Carrier | Entry | Undo variant |
| --- | --- | --- |
| RemovedDativeBonds | RemovedDativeBond | RestoreRemovedDativeBonds |
| RemovedAromaticSystems | RemovedAromaticSystem | RestoreRemovedAromaticSystems |
| RemovedMulticenterBonds | RemovedMulticenterBond | RestoreRemovedMulticenterBonds |
| RemovedNoncovalentBonds | RemovedNoncovalentBond | RestoreRemovedNoncovalentBonds |
| RemovedStereoAtoms | RemovedStereoAtom | RestoreRemovedStereoAtoms |
| RemovedStereoBonds | RemovedStereoBond | RestoreRemovedStereoBonds |

The remaining boxed declarations use existing types:

```rust
RemoveAddedAromaticSystem(Box<AddedAromaticSystem>),
RemoveAddedMulticenterBond(Box<AddedMulticenterBond>),
RemoveAddedStereoAtom(Box<AddedStereoAtom>),
RemoveAddedStereoBond(Box<AddedStereoBond>),
```

Port the study's construction/destructuring changes, including borrowed access in
id_compaction. Do not port its old journal-based correspondence inspection:
S0e has already moved correspondence accumulation out of Undo. The study skipped
test-only constructors and Python, so migrate those callers explicitly. These
carriers hold the existing restoration data; manipulated values retain the
existing no-panic/no-correctness-guarantee contract.

**Checks:** addition/removal rollback, compactions and correspondence, cloned Undo
values, and no panics from manipulated history. Correct restoration is required
for genuine histories under the existing normalized equality/frame contract.
Confirm 88 bytes on the measured 64-bit target; do not promise a portable ABI size.

**Implementation — 2026-10-14.** Added and exported the seven specified public
removal carriers, boxed the eleven selected Undo variants, and migrated all
production/test constructors and matches. Replay moves the stored payload out
of its box; id_compaction borrows the forward compaction. No new constructors,
validation, cloning, correspondence work, or recovery behavior was introduced.
Undo has no Python or DSL representation to migrate. Edit remains unchanged.

The measured aarch64 layout is **Undo 752 → 88 bytes**, with Edit still 256 bytes.
Existing addition/removal replay cases now also restore from cloned histories;
the compaction accessor is exercised through a real topology-removal journal.
Graph-IR all-target compilation with proptest, 1,320 molecule unit cases,
28 edit properties, three compaction properties, strict all-target Clippy,
nightly formatting, and git diff --check passed. Manipulated-history cases retain
their no-panic contract and genuine-history assertions remain unchanged.

**Cost comparison.** Saved before/after release binaries use identical workloads,
entry points, dependency versions and measurement boundaries. Seven synthetic
workloads cover sparse/dense field edits, 1/8 topology removals and 1/8 aromatic
system additions, each with unique/shared storage, plain/tracked execution and
commit/cancellation. These 56 cases are separate probes, not frequency weights.
Six trials rotate before/after/before, with both before labels invoking the same
binary. Each trial uses 20 sample means of 32 calls; the table pools both before
labels and compares medians across trials. Input/Edits preparation and final
molecule destruction are excluded; application, publication, mapping, journal
disposal and rollback are included. All builds and tests finished before timing.

| Operation / atoms / ownership / tracking / outcome | Before, µs | Boxed, µs | Change |
| --- | ---: | ---: | ---: |
| 1 charge edit / 8 / unique / plain / commit | 0.228 | 0.211 | −7.31% |
| 1 charge edit / 96 / unique / plain / commit | 0.499 | 0.503 | +0.64% |
| 64 charge edits / 96 / unique / plain / commit | 4.574 | 2.438 | −46.71% |
| 8 topology removals / 96 / unique / plain / rollback | 108.261 | 108.420 | +0.15% |
| 8 topology removals / 96 / unique / tracked / rollback | 109.413 | 123.172 | +12.57% |
| 1 topology removal / 96 / shared / tracked / rollback | 13.046 | 14.915 | +14.33% |
| 1 aromatic-system addition / 96 / shared / tracked / rollback | 0.662 | 0.720 | +8.76% |

Allocation probes are separate instrumented builds. For 64 charge edits on unique
96-atom input, plain commit's peak extra requested heap falls from 48,128 to
5,632 bytes, with five allocations in both versions. Eight topology removals add
eight allocations (396 → 404 for plain commit); their peak barely changes
(35,908 → 35,860 bytes). Eight aromatic-system additions likewise add eight
allocations (101 → 109), while peak extra heap falls from 7,628 to 3,212 bytes.
These figures exclude allocator overhead, stack and RSS.

**Measured trade-off.** Tracked topology cancellation rises 12.29–14.33%
across unique/shared and 1/8-removal cases, or 1.79–13.76 µs. Before and after
trial ranges do not overlap in those cases; unchanged-control median differences
there are −0.72% to +0.82%. The increase is not dismissed as measurement noise.
The source audit found no added checks, clones, mapping work or repeated replay;
allocation changes match one box per affected undo. This establishes the new
allocation work, but does not explain the disproportionate tracked-only timing
increase. The additional 57.95 ns for the single shared aromatic addition's
tracked cancellation also remains a recorded regression. No compiler attributes
or speculative tuning changes were added. The memory benefit is established;
acceptance of these new costs is not inferred from the earlier S0 decision.
Evidence, source snapshots, binary hashes and the complete per-case timing and
allocation results are in `scratch/181-s1a/`.

**Decision — 2026-10-14.** Retain the compact Undo representation, accepting the
measured rollback costs for its storage and dense-edit benefits. Retain journaling
as an explicitly selectable recovery method alongside cloning. Complete S1b and
compare the combined layouts against the fixed pre-S1 baseline, including whole
operations and both recovery methods; preserve the S1a results separately. This
does not authorize further compiler tuning or leave a deferred rollback
investigation.

#### S1b — Edit: 104-byte layout (completed)

**Modules:** graph-IR edit vocabulary, application, DSL, Python edits and producers.
**Breaking Rust payload shape, red→green within the subitem.** [dep: S0d, S0g]

Box the complete payload once for AddMoleculeConstraint, RemoveMoleculeConstraint,
AddAromaticSystem, AddMulticenterBond, AddStereoAtom, and AddStereoBond. Keep all
other variants inline. Migrate every constructor, match, conversion, and test in
the same subitem, including Python and DSL callers. Boxing alone changes neither
operation semantics nor DSL/Python lifecycle, syntax, or copying policy.

**Production declarations:** add four public open carriers
in ir/edit, with public fields and Clone, Debug, PartialEq, Eq, exported through
ir. AddBond supplies the naming precedent. No separate constructor is needed.

```rust
pub struct AddAromaticSystem {
    pub atoms: Vec<AtomHandle>,
    pub attributes: AromaticSystemForm,
}
pub struct AddMulticenterBond {
    pub atoms: Vec<AtomHandle>,
    pub attributes: MulticenterBondForm,
}
pub struct AddStereoAtom {
    pub site: AtomHandle,
    pub ligands: Vec<(AtomHandle, StereoLigandKind)>,
    pub attributes: StereoAtomForm,
}
pub struct AddStereoBond {
    pub site: BondHandle,
    pub ligands: Vec<(AtomHandle, StereoLigandKind)>,
    pub attributes: StereoBondForm,
}
// Edit variants:
AddAromaticSystem(Box<AddAromaticSystem>),
AddMulticenterBond(Box<AddMulticenterBond>),
AddStereoAtom(Box<AddStereoAtom>),
AddStereoBond(Box<AddStereoBond>),
AddMoleculeConstraint { constraint: Box<ConstraintEdit> },
RemoveMoleculeConstraint { constraint: Box<ConstraintEdit> },
```

The constraint variants box their existing type directly; the study's extra
one-field constraint carriers are unnecessary. Existing Edits addition methods
retain their signatures and perform the boxing. Python constructors retain their
existing named arguments; Box and the Rust allocation carriers are not separate
Python objects. DSL syntax is unchanged. Conversion owns/moves each payload where
it already owns the input, and borrows it where the existing conversion borrows.
Handle resolution and final integrity checks remain application responsibilities.

**Checks:** existing edit application, conversion/serialization and Python parity
cases for the six variants; cloning and consuming conversion preserve payloads.
Confirm 104 bytes on the measured target. **S1 gate:** both layouts coexist,
affected callers pass, and the selected representation has no extra semantic work.

**Implementation and verification — 2026-10-14.** Added and exported the four
specified addition carriers, boxed the six selected Edit payloads, and migrated
construction, application, DSL, Python conversion, and test/property consumers.
Owned application moves each payload out of its box; borrowed conversions retain
their existing copying policy. No constructor checks, public methods, Python
variant fields, DSL syntax, or transaction behavior changed. The measured
aarch64 layouts are Edit **256 → 104 bytes** and Undo **752 → 88 bytes**.

Workspace all-target compilation with graph-IR/graph property features passed.
Verification passed 10,075 graph-IR unit cases, 28 edit properties, three
compaction properties, 1,603 graph operation cases, 1,640 Rust binding cases,
and 342 Python edit/transaction cases after rebuilding the extension with
Python 3.13. The five existing ignored Rust cases remain ignored. Strict
all-target Clippy for the three affected crates, nightly formatting checks, and
git diff --check passed. Existing property laws, generators' domains, expected
errors, and roundtrip assertions are preserved. No workspace-wide test run was
repeated at this subitem. Evidence is in `scratch/181-s1/`.

**Combined S1 comparison.** The baseline is the committed pre-S1 tree at
a1c4b62bc906555dced2245d0aaac092c234fc89; the candidate includes both boxed layouts.
Both use identical harness sources and dependency versions. The whole-operation
study reuses the earlier recovery fixtures at one active group, one active group
plus 256 extra atoms, and 16 active groups, together with its edit/failure cases.
Consuming, clone/install, and native mutating routes are measured separately.
Native mutation journals resolution/projection and transactions; transformer
mutation applies complete plans without journals. No fixture-frequency weighting
or aggregate score is inferred.

Four trials rotate before/after/before, with the two before labels using the same
binary. Each process reports 20 sample means of 32 calls. Summaries compare eight
baseline trial medians with four candidate medians and retain unchanged-control
differences and trial ranges. All builds and checks finished before serial timing.
Whole operations include internal edit construction and disposal. Separate edit
probes preserve S1a's application-only boundary and now select Clone or Journal
explicitly; that harness change is identical on both sides of this comparison.
Allocation probes use separate instrumented builds.

Selected unique-input medians in microseconds; dense resolution/projection use
176 atoms and the transaction probes use 96 unless noted:

| Operation / route | Pre-S1, µs | S1, µs | Change |
| --- | ---: | ---: | ---: |
| Dense resolution / consuming | 475.900 | 460.298 | −3.28% |
| Dense resolution / clone-install | 487.505 | 475.309 | −2.50% |
| Dense resolution / journal | 599.950 | 480.972 | −19.83% |
| Dense projection / consuming | 80.764 | 77.375 | −4.20% |
| Dense projection / clone-install | 100.026 | 87.740 | −12.28% |
| Dense projection / journal | 117.336 | 85.256 | −27.34% |
| 64 charge edits / plain Journal commit | 4.351 | 2.091 | −51.94% |
| 64 charge edits / plain Clone commit | 3.318 | 3.167 | −4.54% |
| Eight topology removals / tracked Journal cancellation | 108.775 | 107.436 | −1.23% |
| Eight topology removals / tracked Clone cancellation | 63.193 | 63.419 | +0.36% |

Dense journaled resolution's peak extra requested heap falls from 1,710.32 to
869.24 KiB; projection falls from 460.94 to 91.90 KiB. Consuming projection falls
from 81.42 to 43.42 KiB and clone/install projection from 137.28 to 99.28 KiB.
Dense resolution's consuming and clone/install peaks are unchanged, although
their allocation traffic decreases. Journaled resolution makes 64 additional
allocations, matching the 32 selected Edit and 32 selected Undo additions;
journaled projection makes three additional allocations. These figures exclude
allocator bookkeeping, stack and RSS.

**Remaining costs and limits.** Sparse shared resolution with 256 extra atoms
through clone/install rises from 113.306 to 116.035 µs (+2.41%, +2.729 µs), with
overlapping trial ranges. Early projection rejection through native mutation rises
from 285.1 to 311.5 ns (+9.26%, +26.4 ns), with disjoint trial ranges. One-charge
Clone transactions on eight atoms show increases of 16–33 ns across several
ownership/tracking/outcome combinations, up to +9.31%; several have disjoint trial
ranges. These increases remain part of the result. No further tuning was attempted.

Tracked topology cancellation in the combined Journal matrix is within −1.23%
to −0.29% of baseline across the four 1/8-removal and unique/shared cases. The
earlier S1a-only 12–14% increase remains recorded above. The current comparison
does not explain that earlier behavior or isolate S1b's contribution across the
changed harnesses; it establishes the complete S1 build's cost in this paired
study. No compiler attributes or speculative code rearrangements were added.

**Assessment:** retain the selected compact layouts. Their storage benefits and
complete-operation gains are demonstrated, with the remaining increases explicitly
reported above. Retain both recovery methods: in the current unique-input probes,
Journal is faster for 64 charge edits, while Clone is faster for eight tracked
topology removals followed by cancellation. Neither is selected universally.
The comparison is complete; no further performance investigation is scheduled.
All 186 whole-operation rows and 56 rows per recovery method, allocation results,
source/binary hashes, correctness logs and reproduction instructions are retained
in `scratch/181-s1/`.

### S2 — Field/key batches and corresponding undos

The stereo subitems follow the [kindless stereo constraint decision](#kindless-stereo-constraints--decision-2026-10-14).
S2g1 implements form semantics and readings for both stereo entity kinds; S2g2
implements their update lowering and existing serialization paths. S2g3 and S2h
then add the atom and bond batches. The remaining shared stereo-key syntax must
be settled before S2g3; it does not block S2g1 or S2g2.

**Interfaces:** use the study's typed vectors and direct
application/replay loops. Put each new variant beside the corresponding entity's
single-entry variants in Edit and Undo. Do not retain ExperimentEditBatch,
ExperimentUndoBatch, the experimental grouping switch, or their wrapper variants.
The tables in S2a–S2h instantiate the declarations below; E, H, I, F, K, and C
are notation in this plan, not new Rust generic types or aliases.

For every field table row, add the named variant to both Edit and Undo:

```rust
// Edit: H is the entity handle; F is the table's form type.
Variant { changes: Vec<(H, F, F)> }, // handle, old, new
// Undo: I is the concrete entity id.
Variant { changes: Vec<(I, F, F)> }, // id, inverse old, inverse new
// Edits method named in the table:
pub fn method(&mut self, changes: Vec<(H, F, F)>);
```

The Edits method moves the vector into one Edit via push. It neither combines
with preceding edits nor examines/reorders rows. An empty vector remains an
empty batch. Construction does not access a molecule or check old values.

For each entity E, add these constraint declarations (the concrete names/types
are listed in each item):

```rust
// Edit:
ModifyEConstraints { key: K, changes: Vec<(H, Option<C>, Option<C>)> },
// Undo:
RestoreEConstraints { key: K, constraints: Vec<(I, Option<C>)> },
// Edits:
pub fn modify_e_constraints(
    &mut self, key: K, changes: Vec<(H, Option<C>, Option<C>)>,
);
```

Stereo constraint Edit/method rows instead have
`(H, Option<StereoKind>, Option<C>, Option<C>)`, retaining geometry context per
row. Their Undo rows have the same `(I, Option<C>)` shape as the other entities.
Every present old/new constraint must have the batch key; application rejects a
mismatch with TransactionError::OldStateMismatch before writing that row.
`(None, None)` retains the single-entry no-op behavior after handle resolution. Keys include
their parameters: two RingMembership scopes or two ligand permutations are not
the same key. No new error enum or batch-entry error wrapper is introduced.

**Python signatures:** expose the same named Edit variants and Edits methods.
A field constructor has `Edit.Variant(changes: list[tuple[int | New, F, F]])`;
its method is `edits.method(changes) -> None`. Constraint constructors have
`Edit.ModifyEConstraints(key: K, changes: list[tuple[int | New, C | None,
C | None]])`; the method takes the same key and changes arguments. Stereo rows
insert `StereoKind | None` after the handle. K and C are the existing Python
constraint-key/form classes; F is the existing form class, without new scalar
coercions. Rust binding arguments use Vec of tuples, HandleLike, and Py<F>/Py<C>
for form objects, following existing Edit field constructors. Reuse each key's
existing conversion. Conversion to owned Rust batches copies supplied Python
form values once where ownership is needed; inputs remain usable. Construction
of Python Edit values retains the existing Python Edit conversion lifecycle.
Edits methods reject consumed receivers through the existing exception. No new
consumed state or view-invalidation mechanism is introduced for row values.

The Edit variant's changes property returns a tuple of row tuples. Reuse the
existing key, form, and constraint classes and their getter behavior. Later
mutation of constructor input lists or mutable input values does not alter the
Edit. Shared getter immutability and Rust/Python payload-name alignment belong
to S11c and do not gate any S2 batch subitem.

**Native DSL — stereo constraint key syntax unresolved:** use singular entity
keys for all 19 field-batch operations and :modify-constraints on each of the
eight entity kinds. Each table gives the exact field verb. Existing :modify,
addition, removal, and replacement syntax stays unchanged. For example:

```clojure
[{:atom {:modify-charge [[0 {:expect 0 :update 1}]
                         [{:new 0} {:expect 0 :update -1}]]}}]
```

A field row is `[handle {:expect <old-form> :update <new-form>}]`; each value uses
its existing form DSL. Stereo configurations use a string containing the existing
complete configuration head (`"*"`, `"Th*"`, `"Th0"`, etc.), with no constraint
tail or omitted-coset update shorthand. A constraint batch is `{:atom {:modify-constraints {:key <key> :changes
[[handle {:expect <old-constraint-or-nil> :update <new-constraint-or-nil>}]]}}}`,
with the appropriate singular entity key. Values use the existing entity constraint
DSL. Retained stereo rows carry `:kind <kind>` in the expect/update map for both
constraint operands. Omit kindless constraint-only rows and both their operands
under the vacuity decision above. No kind is inferred from the molecule.

The key encoding is explicit even for empty batches. Unit keys use their existing
constraint keyword (for example :valence, :aromatic, :electron-count,
:intramolecular, :stereogenicity). Parameterized keys use:

- RingMembership(All): `{:ring-membership {}}`;
  RingMembership(Size(n)): `{:ring-membership {:size n}}`.
- Topicity: `{:topicity {:pair [i j]}}`.

**Rejected — 2026-10-14:** full permutation-image encoding for stereo constraint
keys was not approved and is rejected. LigandSymmetry and Fluxionality key syntax,
including context for empty batches, remains to be settled. Kindless-row omission
is settled by the vacuity decision above. No replacement key encoding is selected.
The existing uniform infallible update-helper interfaces are retained.

Reject malformed key encodings at deserialization. Key/value agreement remains
an application precondition, just as independently supplied old/new values do.

Extend EditsDsl's existing trait boundary and per-entity EditInput organization.
Provide tree readers and tree rendering for the new forms, with
read_edn_/render_ naming where free functions are needed. Do not add a
public batch DSL type or a generic parser that obscures the per-entity grammar.
Apart from omission of vacuous kindless stereo-constraint rows, valid batches
preserve their variant, key, retained row order, empty rows, and stereo kind
context through DSL conversion. Rust/Python value conversion preserves the raw
batch payload. The study's expansion into individual edits is not a production
conversion. No Delta variants are added.

**Reuse from the study:** port the final direct loops, not the initial recursive
prototype. Journal-free loops resolve a handle, compare the current field with
old using normalized_eq, and move new into the mutable view. Journaled loops
reserve one typed inverse vector, record only completed rows, and return one
Undo. On failure they replay that vector in reverse before returning the error.
Constraint loops retain optional old values under the shared key. Replay checks
actual id bounds and writes through mutable views. The atom/bond loops are tested
study code; the other six entities extend that same structure with their types.
Add the explicit row-key/batch-key comparison: the study's trusted planner always
provided matching keys, whereas the production Edit is independently constructible.

All S2 batch subitems use the following execution contract. Each batch subitem
includes its own vocabulary, Edits construction, application, undo replay,
DSL/Python conversion, and tests. Entity application loops are independent;
S2h reuses the stereo configuration/key grammar introduced in S2g3.

- Field values use the same form types and old-value comparisons as the existing
  single-field variants. Field undos preserve the current inverse old/new shape;
  this work does not introduce a separate old-only optimization.
- Constraint old/new absence, key agreement, and stored previous values preserve
  the existing single-entry contract. Molecule-level constraint batching is not
  part of this stage.
- Preserve handle namespaces and entry order within a batch. Repeated changes to
  the same field/key execute sequentially; do not sort or deduplicate them.
- Check an entry before writing it. Journaled failure restores any completed
  prefix of the current batch before earlier journal entries are replayed.
  Journal-free execution retains its destructive failure contract. No full-
  molecule recovery copy, per-entry large Undo vector, recursive Edit replay, or
  internal-panic recovery mechanism is introduced.
- Both execution paths operate through the existing molecule mutation methods
  and mutable views. Attribute batches do not change entity-id correspondence.
- Construction and execution use typed vectors directly. Existing single-entry
  operations do not become one-element heap batches as an implementation shortcut.

**Checks in each subitem:** equivalence to the corresponding ordered single-entry
sequence on valid inputs; empty and single-entry batches; mixed initial/new
handles; repeated sequential changes to one entity; missing/stale handles and
old-value mismatch at first/middle/last entries; constraint add/remove/modify and
key mismatch; rollback after preceding successful edits; tracked correspondence;
no panics from manipulated undo payloads. Include DSL/Python construction and
conversion checks with each entity, preserving typed batches rather than expanding
them into individual edits. Python checks also cover independence from constructor
inputs, the tuple-of-tuples shape of changes, and exact key and row values.
Shared getter mutation-rejection checks belong to S11c.
Property tests exercise sequence equivalence and restoration independently of
the batching implementation. Invalid histories need not produce correct values.
Application unit tests live in molecule/apply/tests, grouped by operation with
method/variant ordering within each module. Extend these logical groups rather
than creating test modules for implementation subitems.

#### S2a — Atom field and constraint batches (done)

**Modules:** graph-IR ir/edit, molecule/apply, dsl/edit; Python edit; their tests.
**Additive vocabulary with exhaustive-match/caller migration; green at completion.**
[dep: S1a, S1b]

**Concrete interfaces:** H = AtomHandle; I = AtomId.
DSL entity key: :atom.

| Edit and Undo variant | F | Edits method | DSL verb |
| --- | --- | --- | --- |
| ModifyAtomElements | ElementForm | modify_atom_elements | :modify-element |
| ModifyAtomIsotopeMasses | IsotopeMassForm | modify_atom_isotope_masses | :modify-isotope-mass |
| ModifyAtomCharges | NumForm | modify_atom_charges | :modify-charge |
| ModifyAtomImplicitHydrogens | NumForm | modify_atom_implicit_hydrogens | :modify-implicit-hydrogens |
| ModifyAtomLonePairs | NumForm | modify_atom_lone_pairs | :modify-lone-pairs |
| ModifyAtomUnpairedElectrons | UnpairedElectronsForm | modify_atom_unpaired_electrons | :modify-unpaired-electrons |

Constraint interfaces: Edit::ModifyAtomConstraints,
Undo::RestoreAtomConstraints, Edits::modify_atom_constraints;
K = AtomConstraintKey, C = AtomConstraintForm.

Add typed atom field batches for Element (ElementForm), IsotopeMass
(IsotopeMassForm), Charge, ImplicitHydrogens, LonePairs (NumForm), and
UnpairedElectrons (UnpairedElectronsForm). Each field has a vector of
(AtomHandle, old, new); the corresponding Undo vector carries AtomId and inverse
old/new values. Add keyed constraint batches using AtomConstraintKey and optional
AtomConstraintForm old/new values; undo retains each accepted previous value.

Wire direct Edits construction, both application paths, and reverse-order replay
through atom mutable views. Include the DSL/Python conversions in this subitem.
Keep existing single-atom variants and methods unchanged.

**Checks:** the shared S2 cases for all six fields and atom constraint keys;
include several independent fields on the same atom, repeated writes to one field,
and atoms created earlier in the same Edits sequence. Check exact values, restored
state, unchanged correspondence, and boundary conversions.

**Implementation and verification — 2026-10-14.** All six field batches and the
keyed atom constraint batch are implemented in Edit, Undo, Edits, both application
paths, DSL conversion, and Python. Direct loops preserve row order, repeated
handles, prefix recovery, accepted-old-value replay, and correspondence. Python
changes exposes tuple rows using the existing form and constraint classes.

Passed 10,268 graph-IR unit cases (three existing ignored cases), 36 edit
properties, and 505 Python edit/transaction/delta cases after rebuilding the
extension. Application tests were reorganized into logical files under
molecule/apply/tests; all 431 cases were retained and passed after the move.
Affected-crate strict Clippy, workspace all-target compilation, nightly formatting,
and git diff --check passed. Workspace tests remain at the plan gate; production
performance comparisons remain in S4a. Logs are in scratch/181-s2a.

#### S2b — Bond field and constraint batches (done)

**Modules:** graph-IR ir/edit, molecule/apply, dsl/edit; Python edit; their tests.
**Additive vocabulary with exhaustive-match/caller migration; green at completion.**
[dep: S1a, S1b]

**Concrete interfaces:** H = BondHandle; I = BondId.
DSL entity key: :bond.

| Edit and Undo variant | F | Edits method | DSL verb |
| --- | --- | --- | --- |
| ModifyBondOrders | NumForm | modify_bond_orders | :modify-order |
| ModifyBondCharges | NumForm | modify_bond_charges | :modify-charge |
| ModifyBondUnpairedElectrons | UnpairedElectronsForm | modify_bond_unpaired_electrons | :modify-unpaired-electrons |

Constraint interfaces: Edit::ModifyBondConstraints,
Undo::RestoreBondConstraints, Edits::modify_bond_constraints;
K = BondConstraintKey, C = BondConstraintForm.

Add typed bond field batches for Order and Charge (NumForm), and UnpairedElectrons
(UnpairedElectronsForm). Entries carry BondHandle and typed old/new values; Undo
entries carry BondId and inverse values. Add BondConstraintKey batches with optional
BondConstraintForm values and accepted-previous-value undo.

Wire Edits construction, both application paths, replay through bond mutable
views, and DSL/Python conversion together. Preserve the existing single-bond API;
batching does not change bond endpoints or introduce topology edits.

**Checks:** shared S2 cases for all three fields and bond constraint keys; include
newly added bonds, a bond removed before use of its handle, sequential order
changes, rollback after earlier atom edits, and unchanged entity correspondence.

**Implementation and verification — 2026-10-14.** Added the three typed field
batches and keyed constraint batch throughout Edit, Undo, Edits, direct application,
journaled application, DSL, and Python. Constraint keys retain their ring scope,
including in empty batches. Ordered application and reverse replay use bond mutable
views and preserve endpoints and correspondence. Python uses the existing form
classes and tuple-of-tuples changes properties.

Verification passed 10,381 graph-IR unit cases (three existing ignored cases),
41 edit properties, and 526 Python edit/transaction/delta cases after rebuilding
the extension. Coverage includes repeated initial/new handles, removed handles,
first/middle/last failures, key/scope mismatches, prior atom edits, normalized-equal
old values, and manipulated undo payloads. Tests extend the existing logical groups.
Affected-crate strict Clippy, workspace all-target compilation, nightly formatting,
and git diff --check passed. Logs are in scratch/181-s2b.

#### S2c — Dative bond field and constraint batches

**Modules:** graph-IR ir/edit, molecule/apply, dsl/edit; Python edit; their tests.
**Additive vocabulary with exhaustive-match/caller migration; green at completion.**
[dep: S1a, S1b]

**Concrete interfaces:** H = DativeBondHandle; I = DativeBondId.
DSL entity key: :dative-bond.

| Edit and Undo variant | F | Edits method | DSL verb |
| --- | --- | --- | --- |
| ModifyDativeBondOrders | NumForm | modify_dative_bond_orders | :modify-order |

Constraint interfaces: Edit::ModifyDativeBondConstraints,
Undo::RestoreDativeBondConstraints, Edits::modify_dative_bond_constraints;
K = DativeBondConstraintKey, C = DativeBondConstraintForm.

Add Order batches with Vec<(DativeBondHandle, NumForm, NumForm)> and corresponding
DativeBondId inverse vectors. Add DativeBondConstraintKey batches with optional
DativeBondConstraintForm values and accepted-previous-value undo.

Wire Edits construction, both application paths, dative-bond mutable-view writes,
undo replay, and DSL/Python conversion. Keep donor/acceptor replacement as the
existing separate operations; no relation reconstruction is part of a field batch.

**Checks:** shared S2 cases for order and dative constraint keys; include newly
added dative bonds, failure after a completed order change, retained donor order
and acceptor, and rollback through the existing removal/compaction path.

#### S2d — Aromatic system field and constraint batches

**Modules:** graph-IR ir/edit, molecule/apply, dsl/edit; Python edit; their tests.
**Additive vocabulary with exhaustive-match/caller migration; green at completion.**
[dep: S1a, S1b]

**Concrete interfaces:** H = AromaticSystemHandle; I = AromaticSystemId.
DSL entity key: :aromatic-system.

| Edit and Undo variant | F | Edits method | DSL verb |
| --- | --- | --- | --- |
| ModifyAromaticSystemElectrons | ElectronCountsForm | modify_aromatic_system_electrons | :modify-electrons |
| ModifyAromaticSystemCharges | NumForm | modify_aromatic_system_charges | :modify-charge |
| ModifyAromaticSystemUnpairedElectrons | UnpairedElectronsForm | modify_aromatic_system_unpaired_electrons | :modify-unpaired-electrons |

Constraint interfaces: Edit::ModifyAromaticSystemConstraints,
Undo::RestoreAromaticSystemConstraints, Edits::modify_aromatic_system_constraints;
K = AromaticSystemConstraintKey, C = AromaticSystemConstraintForm.

Add Electrons (ElectronCountsForm), Charge (NumForm), and UnpairedElectrons
(UnpairedElectronsForm) batches with AromaticSystemHandle entries and
AromaticSystemId inverse entries. Add AromaticSystemConstraintKey batches with
optional AromaticSystemConstraintForm values and accepted-previous-value undo.

Wire Edits construction, both application paths, mutable-view writes, replay, and
DSL/Python conversion. Preserve atom order and the existing open attribute
contract; batching adds no electron-count length check or system reconstruction.

**Checks:** shared S2 cases for all three fields and aromatic constraint keys;
include electron-count vectors, independent charge changes, repeated electron
changes, late batch rejection, and unchanged atom lists after application/rollback.

#### S2e — Multicenter bond field and constraint batches

**Modules:** graph-IR ir/edit, molecule/apply, dsl/edit; Python edit; their tests.
**Additive vocabulary with exhaustive-match/caller migration; green at completion.**
[dep: S1a, S1b]

**Concrete interfaces:** H = MulticenterBondHandle; I = MulticenterBondId.
DSL entity key: :multicenter-bond.

| Edit and Undo variant | F | Edits method | DSL verb |
| --- | --- | --- | --- |
| ModifyMulticenterBondElectrons | ElectronCountsForm | modify_multicenter_bond_electrons | :modify-electrons |
| ModifyMulticenterBondCharges | NumForm | modify_multicenter_bond_charges | :modify-charge |
| ModifyMulticenterBondUnpairedElectrons | UnpairedElectronsForm | modify_multicenter_bond_unpaired_electrons | :modify-unpaired-electrons |

Constraint interfaces: Edit::ModifyMulticenterBondConstraints,
Undo::RestoreMulticenterBondConstraints, Edits::modify_multicenter_bond_constraints;
K = MulticenterBondConstraintKey, C = MulticenterBondConstraintForm.

Add Electrons (ElectronCountsForm), Charge (NumForm), and UnpairedElectrons
(UnpairedElectronsForm) batches with MulticenterBondHandle entries and
MulticenterBondId inverse entries. Add MulticenterBondConstraintKey batches with
optional MulticenterBondConstraintForm values and accepted-previous-value undo.

Wire Edits construction, both application paths, mutable-view writes, replay, and
DSL/Python conversion. Preserve the atom sequence and existing electron-count
semantics; neither validation nor structural replacement is added by batching.

**Checks:** shared S2 cases for all three fields and multicenter constraint keys;
include newly added bonds, electron-count payloads, sequential changes, failure
following earlier successful field edits, and exact retained atom order.

#### S2f — Noncovalent bond field and constraint batches

**Modules:** graph-IR ir/edit, molecule/apply, dsl/edit; Python edit; their tests.
**Additive vocabulary with exhaustive-match/caller migration; green at completion.**
[dep: S1a, S1b]

**Concrete interfaces:** H = NoncovalentBondHandle; I = NoncovalentBondId.
DSL entity key: :noncovalent-bond.

| Edit and Undo variant | F | Edits method | DSL verb |
| --- | --- | --- | --- |
| ModifyNoncovalentBondKinds | NoncovalentBondKindForm | modify_noncovalent_bond_kinds | :modify-kind |

Constraint interfaces: Edit::ModifyNoncovalentBondConstraints,
Undo::RestoreNoncovalentBondConstraints, Edits::modify_noncovalent_bond_constraints;
K = NoncovalentBondConstraintKey, C = NoncovalentBondConstraintForm.

Add Kind batches with Vec<(NoncovalentBondHandle, NoncovalentBondKindForm,
NoncovalentBondKindForm)> and corresponding NoncovalentBondId inverse vectors.
Add NoncovalentBondConstraintKey batches with optional NoncovalentBondConstraintForm
values and accepted-previous-value undo.

Wire Edits construction, both application paths, mutable-view writes, replay, and
DSL/Python conversion. Leave the two atom IDs and their replacement API unchanged.

**Checks:** shared S2 cases for kind and noncovalent constraint keys; include
literal/undetermined forms, new and removed handles, repeated kind changes, and
unchanged atom IDs and correspondence through successful execution and rollback.

#### S2g1 — Stereo form semantics and constraint readings

**Modules:** graph-IR ir/stereo, ir/view/stereo, ir/view/constraints; Python
stereo and constraint/stereo; their unit and property tests.
**Semantic implementation changes; signatures unchanged, green at completion.**
[dep: S1a, S1b]

Implement kind-aware Normalize and Lattice behavior for StereoAtomForm and
StereoBondForm at the entity level. Clear constraints before normalizing their
contents when configuration is Undetermined. Retain the existing known-kind
behavior, including a known kind with an undetermined coset. Make normalized_eq,
is_undetermined, is_ground, meet, join, matches, and is_compatible agree with
that normalization. Neither a meet nor a join may make vacuous input constraints
meaningful; a kindless result has vacuous constraints.

Apply the same rule to asserted/asserted_complete and their evidenced/comparison
consumers in StereoAtomConstraintsView and StereoBondConstraintsView. Preserve
the existing derived readings. Cover both molecule-backed and value-backed
Python views through the owning form's kind. Keep raw storage access, equality,
hashing, mutation, and the standalone constraint containers' semantics intact;
do not normalize eagerly or change the general-purpose derives. Update the
affected API rustdoc to state the raw/semantic distinction.

**Checks:** exact kindless normal forms with ordinary and internally contradictory
payloads; meaningful constraints on known-kind/undetermined-coset forms; meet in
both operand orders; joins losing the kind; and all four stereo constraint keys.
Extend the existing stereo/lattice/constraint-reading properties to cover raw
kindless payloads and invariance under prior normalization, retaining the lattice
laws. Check raw payload retention separately from semantic equality and verify
Python parity for owned forms and molecule views. Standalone container readings
must still interpret their contents without an owning kind.

In tests/property/constraint_readings, retain view/container assertion equality
for known-kind stereo and the other entity kinds. For kindless stereo, assert
the key's undetermined form on the entity view independently of the stored
container's reading. Generate raw kindless forms with nonempty and contradictory
payloads; do not normalize them before exercising the operation under test.
Retain the existing normalized lattice laws and add raw-versus-normalized
operation/readings agreement on normalizable inputs, including kindless inputs
whose payload would contradict if interpreted independently.

**Performance:** repeat S0a's bounded stereo cases after these changes and record
time and allocation differences. Inspect added cloning, repeated normalization,
or lost allocation-free paths before accepting a regression on unchanged
semantics. Keep common and unusual cases separate in the results.

#### S2g2 — Stereo updates, lowering, and single-entry serialization

**Modules:** graph-IR ir/stereo, ir/edit, dsl/stereo, dsl/edit and their molecule
serialization consumers; Python stereo/edit; affected unit and property tests.
**Semantic implementation changes; signatures unchanged, green at completion.**
[dep: S2g1]

Keep StereoAtomForm::update and StereoBondForm::update as infallible raw updates.
Make both difference_to methods ignore differences that remain kindless.
In Edits::update_stereo_atom and Edits::update_stereo_bond, derive constraint
context from the resulting configuration: Unchanged retains the supplied source
kind, Kinded supplies its kind, and Undetermined supplies none. Omit constraint
contributions without that context while retaining configuration changes. Use
actual stored old constraints for retained edit preconditions. Keep the existing
signatures and application checks; add no kind-presence or kind-match rejection.

Harmonize StereoAtomDsl, StereoBondDsl, both update DSL types, and the single-entry
EditsDsl conversion. Omit kindless constraints before kind-dependent formatting
or conversion can fail. Omit an entire kindless constraint edit, including its
expectation, rather than emitting a removal. Preserve configuration intent and
the order of retained operations. Use the existing kind-scoped grammar and trait
boundaries; Rust/Python value conversion continues to retain raw payloads.
Document entity roundtrips under normalized_eq and the explicit omission rule
for standalone updates/edits; do not promise context-free update serialization
preserves a constraint change that only acquires a kind when applied to a source.

**Checks:** raw update retention; no semantic difference between kindless forms;
Unchanged with and without a source kind; explicit clearing without fallback to
the old kind; and Kinded updates with retained raw old-state preconditions.
Check lowered edits through both application paths. Test Display,
ToEdn, molecule export, and the public EditsDsl conversion with constructed
kindless payloads, including contradictory ones. Assert exact omission, distinct
Unchanged/Undetermined update heads, retained surrounding edits, and absence of
panics through Python's corresponding export paths. Extend the existing stereo
serialization properties to cover these raw states under the stated equivalence;
retain exact roundtrips on the existing lossless domain.

Add generated difference/update and lowering/application properties for both
stereo families. For normalizable source and target forms, applying
source.difference_to(target) to that raw source must be normalized_eq to the
target. For a single update whose source and direct result normalize, applying
the edits lowered against that same raw source must be normalized_eq to the
direct raw form update, through both application paths. Include raw kindless
payloads, known-kind/undetermined-coset forms, kind assignment, and kind clearing.
Compare the immediate semantic results; these laws do not require raw equality or preservation of latent
payload for a subsequent kind assignment. Preserve the exact stored old-value
preconditions in the generated lowering cases.

#### S2g3 — Stereo atom field and constraint batches

**Modules:** graph-IR ir/edit, molecule/apply, dsl/stereo, dsl/edit; Python edit;
their tests.
**Additive vocabulary with exhaustive-match/caller migration; green at completion.**
[dep: S1a, S1b, S2g2]

**Design prerequisite:** settle the shared stereo-key syntax recorded in S2.

**Concrete interfaces:** H = StereoAtomHandle; I = StereoAtomId.
DSL entity key: :stereo-atom.

| Edit and Undo variant | F | Edits method | DSL verb |
| --- | --- | --- | --- |
| ModifyStereoAtomConfigurations | StereoConfigurationForm | modify_stereo_atom_configurations | :modify-configuration |

Constraint interfaces: Edit::ModifyStereoAtomConstraints,
Undo::RestoreStereoAtomConstraints, Edits::modify_stereo_atom_constraints;
K = StereoAtomConstraintKey, C = StereoAtomConstraintForm.

Add Configuration batches with Vec<(StereoAtomHandle, StereoConfigurationForm,
StereoConfigurationForm)> and corresponding StereoAtomId inverse vectors.
Add StereoAtomConstraintKey batches whose entries retain StereoAtomHandle,
Option<StereoKind>, and optional StereoAtomConstraintForm old/new values. Undo
retains the resolved ID and accepted previous constraint value.

Wire Edits construction, both application paths, stereo-atom mutable-view writes,
replay, and DSL/Python conversion. Keep geometry context per entry: one shared
constraint key does not imply a shared kind. Preserve existing configuration
comparison and restoration semantics, without adding coset normalization/checks.
Site and ligand mutation remain separate existing operations.

Extend S2g2's serialization omission to constraint batch rows: a row with no
kind contributes neither an expectation nor an update. Preserve retained row
order and configuration batches; do not infer context from a molecule or turn
omission into removal. Raw batch construction, Python value conversion,
application, and undo keep the corresponding single-entry contracts.

In dsl/stereo, reuse stereo_configuration and fmt_stereo_configuration to expose
crate-private read_edn_stereo_configuration, read_stereo_configuration, and
render_stereo_configuration for the two edit grammars. The readers require the
complete configuration head and reject trailing constraint text. These are
boundary functions for EditsDsl, not a new public configuration type. The shared
stereo-key grammar remains unresolved as recorded in S2. Settle it before
implementing the constraint readers; the atom/bond batch readers use that grammar
with their corresponding key type.

**Checks:** shared S2 cases for configuration and stereo-atom constraint keys;
include different kinds within a keyed batch, per-entry DSL kind context,
configuration changes followed by rejection, and retained site/ligand order.
Add mixed kinded/kindless and entirely kindless constraint batches, including
contradictory omitted payloads. Assert exact serialized omission and retained
order under the settled key grammar, raw Rust/Python payload preservation, and
that export neither panics nor introduces a constraint removal.

#### S2h — Stereo bond field and constraint batches

**Modules:** graph-IR ir/edit, molecule/apply, dsl/edit; Python edit; their tests.
**Additive vocabulary with exhaustive-match/caller migration; green at completion.**
[dep: S1a, S1b, S2g2, S2g3 (shared stereo DSL grammar)]

**Concrete interfaces:** H = StereoBondHandle; I = StereoBondId.
DSL entity key: :stereo-bond.

| Edit and Undo variant | F | Edits method | DSL verb |
| --- | --- | --- | --- |
| ModifyStereoBondConfigurations | StereoConfigurationForm | modify_stereo_bond_configurations | :modify-configuration |

Constraint interfaces: Edit::ModifyStereoBondConstraints,
Undo::RestoreStereoBondConstraints, Edits::modify_stereo_bond_constraints;
K = StereoBondConstraintKey, C = StereoBondConstraintForm.

Add Configuration batches with Vec<(StereoBondHandle, StereoConfigurationForm,
StereoConfigurationForm)> and corresponding StereoBondId inverse vectors.
Add StereoBondConstraintKey batches whose entries retain StereoBondHandle,
Option<StereoKind>, and optional StereoBondConstraintForm old/new values. Undo
retains the resolved ID and accepted previous constraint value.

Wire Edits construction, both application paths, stereo-bond mutable-view writes,
replay, and DSL/Python conversion. Retain kind context per constraint entry and
the existing configuration semantics. Do not narrow the representation to
coset-only changes or change the separate site/ligand replacement operations.

Reuse the form, reading, lowering, and single-entry serialization semantics from
S2g1/S2g2. Apply S2g3's whole-row omission rule to stereo-bond constraint batch
serialization, preserving raw construction/conversion and application/undo
contracts. Do not add a separate bond-specific kind policy.

**Checks:** shared S2 cases for configuration and stereo-bond constraint keys;
include new and removed handles, per-entry geometry context, sequential
configuration changes, failure after earlier batches, and retained site/ligand order.
Cover mixed kinded/kindless and entirely kindless batches, exact omission without
removal, retained order, raw Rust/Python payload preservation, and panic-free export.

**S2 gate:** all eight entity kinds have the same applicable semantics, all
workspace consumers compile and pass their relevant tests, and the new variants
do not defeat the selected compact enum layouts. If a payload would widen an
enum, settle its representation before proceeding; do not silently undo S1.
Run the affected stereo/lattice/constraint-reading property suites explicitly,
alongside the focused Rust and rebuilt Python suites. Both stereo entity kinds
must obey the vacuity decision across normalization, readings, updates, lowering,
and existing and batched DSL export; raw mutation remains infallible.
The S2g1 time/allocation comparison against S0a must be recorded and reviewed,
with any regression's cause and disposition stated before this gate closes.

#### S2i — Edits DSL fuzz target

**Modules:** umol-graph-ir fuzz and scripts/run-fuzz.sh. **Additive, green.**
[dep: S2a–S2h]

Add `fuzz_edits` beside the existing targets: parsing arbitrary input as
EditsDsl must not panic, and a parsed value converts to Edits, renders back
through EditsDsl, and reparses without panicking. Constructor-based kindless
cases remain covered in S2g2/S2g3/S2h; parser fuzzing alone cannot supply those
raw payloads. EditsDsl has one parser, the tree path behind from_edn_str, so there
is no streaming parity to check. Register it in the fuzz Cargo.toml and the
GRAPH_IR_TARGETS list in scripts/run-fuzz.sh before S2j removes fuzz_constraints.

#### S2j — Remove fuzz_constraints

**Modules:** umol-graph-ir fuzz and scripts/run-fuzz.sh. **Removal, green.**
[dep: S2i]

fuzz_constraints parses ConstraintDsl and ConstraintsDsl through the tree path
only; fuzz_molecule already reaches both constraint parsers and checks their
streaming and tree paths against each other. Remove the target, its Cargo.toml
entry, its corpus and artifacts, and its entry in scripts/run-fuzz.sh. Carry its
seeds over to fuzz_molecule, each wrapped in a molecule that places the
constraint where it occurs.

### S3 — Resolver and projector emission

#### S3a — Resolution plans

**Modules:** graph resolver and its phase planners; settled Edits emission methods.
**Implementation change, green.** [dep: S2a–S2h]

Emit field/key batches while producing updates. Do not construct an individual
Edit list and then pack it. Group only independent changes established by the
planner; preserve order for repeated changes and keep structural operations and
phase boundaries in place. Leave chemistry decisions, reports, public lifecycle,
and recovery choice unchanged. The same plans serve journaled and journal-free
execution. Do not force batching into consumers with isolated modifications.

**Checks:** whole-resolution equivalence, underdetermination/contradiction and
source recovery, phase-plan contents, both ownership conditions, and tracked
results where used. The ordinary raised fixtures have determined isotopes; retain
an explicitly undetermined isotope fixture when measuring isotope edit emission.

#### S3b — Projection plans

**Modules:** stereo, aromaticity and isotope projection; composite projector.
**Implementation change, green.** [dep: S2a–S2h]

Emit typed field/key batches directly, retaining phase order and ordinary removal
operations. Additional resolved carbon atoms undergo isotope-default elision;
do not describe them as untouched input in these benchmarks. Do not change the
projection semantics or its existing error and publication boundaries.

**Checks:** complete-output equivalence; early rejection; late isotope rejection
after stereo/aromatic changes with restoration of the original; both execution
paths and ownership conditions. **S3 gate:** production planners exercise the
batch variants without post-construction packing or temporary per-entry journals.

### S4 — Verification and decision record

#### S4a — Measure the production representations and wiring

**Modules:** existing editor/resolver benchmarks; this discussion record.
**Verification/documentation, green.** [dep: S3a, S3b]

Compare the production implementation with the recorded pre-change measurements
and experimental candidates. Report boxing and batching effects separately where the
reference permits: whole-operation time, allocation count/bytes, peak extra
memory, and Edit/Undo counts/capacities. Measure the implemented complete-plan
paths for individual resolvers and transformers. Compare no-recovery,
clone/install, and journal paths for composite operations and prepared
transactions, including failures and structural edits that incur boxes.
Keep interpretation tied to decisions: a smaller outer journal is not the whole
memory footprint, and unchanged resolution time must not be presented as a speedup.

Run focused checks per subitem and the affected stage gates; reserve the full
workspace/Python checks and MSRV check for the final combined gate, rather than
repeating them for every entity kind. Record any regression and its practical
tradeoff before marking this work complete. Preserve durable findings here before
removing scratch artifacts.

### S5 — Store and borrow the Rust chemistry model in Python

#### S5a — ChemistryModel storage, interfaces, and operation callers

**Modules:** umol-py model, molecule, reaction, and their existing tests.
**Internal Rust signature changes, red→green within the subitem; Python call
shapes preserved.** [dep: S0a; independent of S1–S4]

Replace the four stored Python component wrappers with the complete Rust model:

```rust
#[pyclass(eq, frozen, skip_from_py_object)]
#[derive(Clone, Debug, PartialEq)]
pub struct ChemistryModel(GraphChemistryModel);

// Internal conversion boundary.
pub(crate) fn from_rust(model: GraphChemistryModel) -> Self;
pub(crate) fn to_rust(&self) -> &GraphChemistryModel;

// Existing keyword-only Python constructor; component arguments are borrowed.
fn new(
    connectivity: PyRef<'_, ConnectivityModel>,
    valence: PyRef<'_, ValenceModel>,
    aromaticity: PyRef<'_, AromaticityModel>,
    stereo: PyRef<'_, StereoModel>,
) -> Self;
fn default() -> Self;

// Existing Python properties, returning independent values of the same types.
fn connectivity(&self) -> ConnectivityModel;
fn valence(&self) -> ValenceModel;
fn aromaticity(&self) -> AromaticityModel;
fn stereo(&self) -> StereoModel;
```

Construction converts the supplied components once and owns the resulting Rust
model. Inputs remain usable; no Python borrow is stored. default moves the Rust
default into the wrapper, preserving its static borrowed registry/table data.
from_rust moves an owned model; a genuinely borrowed caller explicitly clones
only when it requires an independent owner. Do not add a second conversion API.
The wrapper remains frozen and value-comparable. Remove automatic value-cloning
argument extraction; retain Rust Clone for explicitly requested internal copies.

Preserve the four property types and their independent-value semantics. Adapt
repr to the stored Rust value, retaining its existing Python spelling. Compare
against the Rust default directly; do not construct a complete Python model just
to recognize it. No live component views, caching, Arc storage, consuming model
operations, or Rust model/Resolver lifetime changes are part of this migration.
The constituent Python model types keep their existing public interfaces/storage.

Change each operation's binding argument from Option<ChemistryModel> to
Option<PyRef<'_, ChemistryModel>>. Keep the existing Python name, optionality,
keyword behavior, and defaults. The complete current production caller list is:

- Molecule: from_smiles, to_smiles, resolve, resolve_with_report, resolve_into,
  resolve_into_with_report;
- Reaction: from_reaction_smiles, to_reaction_smiles.

Each caller borrows the stored GraphChemistryModel through to_rust for the whole
synchronous Rust operation. When omitted, construct the operation's existing
local Rust default and borrow it. Preserve the SMILES-specific valence defaults
and all resolve/IO configuration behavior. Keep the PyRef alive until the Rust
operation returns; neither resolver nor result retains it. Use direct local
selection, without a new shared model-selection helper. Migrate internal Rust
callers and tests with the signature change so this subitem ends green. Recheck
the caller inventory at implementation time for additions since the study.

This changes model passing only: it does not change molecule ownership, rollback,
result publication, reports, errors, or view invalidation. Subsequent lifecycle
work must retain borrowed model passing rather than restore owned extraction.

**Checks:** default/custom model value equality, constructor inputs remaining
usable, frozen properties, exact getter values and types, equality and repr,
repeated supplied-model use, and all eight operation boundaries. Retain existing
chemistry successes/failures and source-lifecycle assertions. Verify omitted
versus explicitly supplied defaults according to each operation's actual default,
including the different SMILES valence model. Follow the Python build skill before
compiling and rebuild the extension before Python integration tests.

#### S5b — Verify production construction/use costs and document the contract

**Modules:** Python model/operation integration checks, Python API guide, this
record; isolated measurement harnesses. **Verification/documentation, green.**
[dep: S5a]

Repeat the Study 1 conversion and construction-inclusive comparisons against the
production wrapper, using the recorded baseline. Include omitted/explicit default,
custom registry/table, construction plus one/ten uses, reuse, and component getters.
Measure Python dispatch and return wrapping, Rust allocations, and separate Python
peak memory with the same boundaries. Do not require every one-use fixture to be
faster: custom construction dominates some cases. Confirm that repeated operation
calls borrow the model without rebuilding it, and report any new cost rather than
inferring success from a smaller wrapper implementation.

Document that each Python ChemistryModel owns its Rust model and operations borrow
it for the call; custom construction and component getters may copy. Update the
Python API guide without imposing this representation on unrelated wrappers.
Record final results here. **S5 gate:** unchanged Python model behavior, no implicit
per-call model copy, green affected Rust/Python tests, and production costs checked
against the bounded study. Broader workspace/MSRV verification remains a final
combined-plan gate, not a per-subitem repetition.

### S6 — Python batch copying and owned inverse conversion

#### S6a — Explicit Edits copy

**Modules:** umol-py edit and its Rust/Python tests. **Additive, green.**
[dep: S0a; independent of S1–S5]

Add the Python method with this binding signature:

```rust
// Edits
fn copy(&self) -> PyResult<Self>;
```

Clone the retained GraphIrEdits directly through the existing to_rust/from_rust
boundary. Preserve entry order, duplicates, handles, and payloads exactly; do not
iterate through Python entries, normalize, rebase, or validate against a molecule.
The copy is an independent, unconsumed Edits. The source remains usable and its
iterators remain valid. Copying a consumed source raises the existing ConsumedError.
Consumption on application remains unchanged; users explicitly copy when needed.

**Checks:** empty and populated batches, all entity kinds, field/addition/constraint
payloads, exact copied values, independent subsequent mutation, consuming either
batch without consuming the other, retained source iteration, and consumed-source
failure. Compare application of independently copied batches on equal starting
molecules. Re-run Study 3 direct-copy and copy-plus-apply measurements using the
production method, recording the enum layout used.

#### S6b — Explicit Deltas copy

**Modules:** umol-py delta and its Rust/Python tests. **Additive, green.**
[dep: S0a; independent of S6a]

```rust
// Deltas
fn copy(&self) -> Self;
```

Clone GraphIrDeltas directly and wrap the owned result. This exposes Rust Clone
as an explicit Python operation. Preserve exact order, duplicates, and payloads;
copying does not normalize or check applicability. Both containers remain usable,
and appending to either does not affect the other. Existing indexing and lazy
iteration retain their independent-entry semantics. No entry-view types are added.

**Checks:** exact copies for empty and populated batches, all entity kinds and
constraint deltas, independent appends, unchanged source iterators, and payload
independence. Re-run the Study 3 Deltas direct-copy comparison through Python
dispatch. Document copy alongside construction and append.

#### S6c — Move owned inverse results into Python forms

**Modules:** umol-py entity, delta, and affected conversion callers/tests.
**Internal conversion signature changes, red→green within the subitem; Python
inverse signatures and semantics unchanged.** [dep: S0a; independent of S6a/S6b]

Use ownership at the existing conversion boundary rather than adding a second
inverse implementation or a parallel conversion family:

```rust
// Readonly<T>: accept the owned Rust form instead of borrowing and cloning it.
pub(crate) fn from_rust(py: Python<'_>, value: T::RustForm) -> PyResult<Self>;

// Each of the eight entity Delta wrappers: accept its owned Rust delta.
// AtomDelta shown; use the corresponding types for the other seven kinds.
pub(crate) fn from_rust(py: Python<'_>, delta: GraphIrAtomDelta) -> PyResult<Self>;

// Aggregate Delta: move entity payloads into the corresponding converters.
pub(crate) fn from_rust(py: Python<'_>, delta: GraphIrDelta) -> PyResult<Self>;

// Existing Python method remains unchanged.
fn inverse(&self, py: Python<'_>) -> PyResult<Py<Self>>;
```

Cover AtomDelta, BondDelta, DativeBondDelta, AromaticSystemDelta,
MulticenterBondDelta, NoncovalentBondDelta, StereoAtomDelta, StereoBondDelta, and
aggregate Delta.
Continue calling Rust inverse. Match its owned result by value and move Add/Remove
forms into the existing read-only Python form representation; move owned atom,
ligand, and site lists through their existing conversions as applicable.
Do not borrow the owned inverse merely to clone its complete form again.

Aggregate Delta.inverse passes its owned Rust result through the outer converter,
which dispatches entity payloads by value without an intervening clone.
Migrate genuinely borrowed producers, including indexing and iteration, to make
their required clone explicit before this boundary.
Migrate all Readonly from_rust callers, including Edit conversion, in the same
subitem. Their required independent copies stay in place; no supplied Python form
is consumed. Keep existing field-change and constraint conversion semantics; this
bounded item does not claim to eliminate their nested conversion copies or change
Rust inverse. No normalization or integrity checks are introduced.

**Checks:** exact inverse values and unchanged inputs for every entity kind and
variant; retained read-only form behavior; exact borrowed-entry conversions and
independence after source/container mutation. Exercise both entity inverse and
aggregate Delta inverse. Repeat the measured AtomDelta Add/Remove cases and add
allocation comparisons for overlay forms with owned lists, covering both direct
entity and aggregate inverse paths against S0a. Confirm removal of the second
whole-form copy in both paths without claiming every nested allocation disappears.

**S6 gate:** focused Rust/Python suites pass after rebuilding with the repository
Python environment; the Python API guide records explicit copy semantics and
unchanged inverse behavior. Record production timing/allocation results beside
Study 3, distinguishing measured improvements from unmeasured payloads. Run the
full workspace and MSRV checks only at the combined final gate. S6 is scheduled
work, not deferred; no production implementation is performed by adding this plan.

### S7 — Borrow stored form operands in Python operations

#### S7a — Shared method generation and the eight entity forms

**Modules:** umol-py lattice, atom, bond, dative, aromatic, multicenter,
noncovalent, stereo, and their tests. **Internal additive macro support and caller
migration, green; no public interface change.** [dep: S0a; independent of S1–S6]

Extend the existing impl_py_normalize and impl_py_lattice macros with an explicit
borrowed-operand arm. Keep the existing owned-conversion arm for representations
that reconstruct Rust values. The borrowed arm obtains references through the
wrappers' existing to_rust methods; no new trait, wrapper, cache, or public method
is needed. Preserve shared generation rather than duplicating the method bodies
for each entity kind.

Migrate AtomForm, BondForm, DativeBondForm, AromaticSystemForm,
MulticenterBondForm, NoncovalentBondForm, StereoAtomForm, and StereoBondForm.
Preserve these existing binding signatures and their visibility:

```rust
fn is_undetermined(&self, py: Python<'_>) -> PyResult<bool>;
fn is_ground(&self, py: Python<'_>) -> PyResult<bool>;
pub(crate) fn matches(&self, py: Python<'_>, target: &Self) -> PyResult<bool>;
fn is_compatible(&self, py: Python<'_>, other: &Self) -> PyResult<bool>;
fn meet(&self, py: Python<'_>, other: &Self) -> PyResult<Option<Self>>;
fn join(&self, py: Python<'_>, other: &Self) -> PyResult<Self>;
fn normalized_eq(&self, py: Python<'_>, other: &Self) -> PyResult<bool>;
fn normalize(&self, py: Python<'_>) -> PyResult<Self>;
```

The first seven methods call the corresponding borrowed Rust operation without
preliminary operand clones. meet and join still convert their owned Rust results
into independent Python forms. normalize explicitly clones the borrowed source
before invoking consuming Rust normalization and returns an independent form.
Do not change storage, constructors, setters, read-only flags, or result mutability.

Preserve meet's ordinary None, join's NoJoinError, normalization exceptions, and
normalized_eq semantics. Do not add integrity checks or change Rust lattice
algorithms. Python names, arguments, and source-preserving behavior remain the
same. The stereo equality defect belongs to 231, not this item.

**Checks:** exercise all eight forms, including populated payloads, equal and
unequal operands, incompatible meets, absent joins where the Rust operation
admits them, normalization successes/failures, and unchanged sources. Compare
complete Rust result values where Python value equality is absent. Check that
returned forms are independent and that the owned-conversion macro arm retains
existing leaf-form behavior. Follow the test-writing and Python-build skills;
run focused Rust tests and rebuilt Python integration tests.

#### S7b — Entity constraint containers and molecule Constraints

**Modules:** umol-py constraint atom, bond, dative, aromatic, multicenter,
noncovalent, stereo, molecule, and their tests. **Internal migration, green;
no public interface change.** [dep: S7a]

Use the borrowed macro arm for AtomConstraintsForm, BondConstraintsForm,
DativeBondConstraintsForm, AromaticSystemConstraintsForm,
MulticenterBondConstraintsForm, NoncovalentBondConstraintsForm,
StereoAtomConstraintsForm, and StereoBondConstraintsForm. Apply the same eight
method signatures and semantics as S7a: seven operations borrow their operands;
normalize retains its independent-result copy.

For the molecule-level Constraints container, migrate only its existing
normalization-method generation: normalized_eq borrows both operands, while
normalize clones its source. Do not add lattice operations to Constraints.
Individual constraint forms and molecule/relational constraint enums keep their
existing conversion paths. This item changes no constraint mutation or view API.

**Checks:** all eight entity containers and Constraints; empty/populated values,
distinct equal and unequal operands, exact results, existing absence/errors,
source preservation, and independent returned containers. Include the Study 4
valence-set and charge-sum fixtures. Keep checks independent of the macro bodies.

#### S7c — Verify production copy elimination

**Modules:** affected Python operation tests, Study 4 measurements, and this
record. **Verification/documentation, green.** [dep: S7a, S7b]

Repeat Study 4's existing-storage comparisons through the production Python
methods. Include complete dispatch and result conversion/destruction, with timing
and allocation counting separate. Verify removal of preliminary operand clones
for the 17 wrappers; do not require normalization or result creation to become
allocation-free. Retain small/default cases alongside populated forms, and report
regressions rather than asserting a universal speedup.

Record the final measurements and scope here. No enum-backed leaf or individual
constraint representation migration, speculative binding framework, equality fix,
or additional production benchmark framework is included. **S7 gate:** focused
Rust/Python tests pass, public contracts are unchanged, and measured allocation
changes match the intended removal. Workspace/MSRV checks remain at the combined
final gate. S7 is scheduled work, not deferred.

### S10 — Resolution, projection, transformations, and owned consumers

At S10 entry, every call to graph-IR editing/transaction methods is migrated. S10
changes the graph operations' own public interfaces and migrates their callers
within each subitem; it does not repair a build left broken by S0.

| Boundary | State before S10 | Change owned by S10 |
| --- | --- | --- |
| Composite resolve/resolve_with_report/project | Consume Molecule; use edit internally | Unchanged Rust signatures |
| Composite _into methods | Existing signatures; one transaction across phases, using Journal after S0f | Expose the approved selector and forward it through that transaction |
| Individual resolvers | Consuming ordinary methods; _into applies a complete plan in place without a transaction | No further ownership/recovery change in S10; implemented by 236 |
| Transformer | Consuming transform; transform_iter clones once per result; transform_into applies a complete plan without a transaction | No further ownership/recovery change in S10; implemented by 236 |
| Python resolve/resolve_with_report | Still take the molecule; Option and consumed checks remain | Clone the source and consume the clone; preserve source/views, enabling S11a to remove Option |

Each S10 item includes its changed Rust/Python callers, tests, benchmark call
sites, and the corresponding guide paragraphs. Owned ingestion/parsing and
projected export intermediates keep calling the consuming methods without an
extra clone.

#### S10a — Composite resolver and ingestion/export

**Modules:** graph resolve, Python resolve bindings, consumers and tests.
**Breaking _into selector, red→green within the subitem.** [dep: S0c2, S0g]

Keep the consuming resolve/resolve_with_report/project signatures and all
Solution, report, and error types. Ordinary resolve builds no report. Ingestion,
parsing, and export keep passing their owned molecules. All phases share the
selected transaction for _into execution; no per-phase recovery scopes or
additional source copy. Only Determined commits. Other outcomes preserve the
receiver and diagnostics.

Preserve the planner applicability and integrity assertions established by
[236](236-operation-ownership-copy-costs-2026-10-09.md) S4c. Returned resolver
errors remain planning errors; application, probe, and publication failures of
internally produced plans are implementation defects. Externally supplied edits
remain checked by the editor and transaction APIs, with failure coverage in S0.

Update Python resolve/resolve_with_report to clone the source, consume the clone,
and preserve the source's views. The two _into methods retain their receiver and invalidation contracts and
expose the approved rollback selector. Keep ChemistryModel borrowing from S5 if
already implemented; do not restore extraction copies.

**Checks:** complete determined results, underdetermination, contradictions,
planning errors, optional report equivalence, late projection rejection,
ingest/export fixtures, and unique/shared inputs. Count copies at call boundaries
rather than inferring their number from API names.

#### S10b — Dropped

The individual resolvers keep their consuming ordinary methods. Their in-place
_into methods already apply complete plans without transactions, implemented in
[236](236-operation-ownership-copy-costs-2026-10-09.md) (S4f).

#### S10c — Dropped

transform keeps consuming and transform_iter keeps its per-result clone.
transform_into's in-place application without a transaction is implemented in
[236](236-operation-ownership-copy-costs-2026-10-09.md) (S4d, S4e).

### S11 — Python storage and API guarantees

#### S11a — Remove consumed Molecule storage

**Modules:** umol-py molecule, all molecule consumers/views, error tests and exports.
**Internal breaking conversion signatures, red→green within the subitem.**
[dep: S0c2, S10a]

Store value: GraphIrMolecule and counter: u64 directly. Remove Molecule.take.
Use these crate-private signatures:

```rust
fn to_rust(&self) -> &GraphIrMolecule;
fn to_rust_mut(&mut self) -> &mut GraphIrMolecule;
fn view_counter(&self) -> u64;
fn advance_counter(&mut self);
fn check_access(&self, expected: u64, accessor: &'static str) -> PyResult<()>;
```

Keep from_rust taking an owned molecule. check_access checks counter equality,
not consumed state. Migrate every binding caller in the same subitem, removing
only fallibility caused by consumed Molecule. PyO3 borrow failures and actual
operation errors remain. Keep ConsumedError for Edits, terminal editor finish,
and consumed RingViews. RingViews retains its optional storage and the
InvalidatedViewError contract for dependent accessors after consumption.
Keep all *Storage types and counter terminology; add no replacement lifecycle
wrapper, weak-reference list, or cache.

**Checks:** source access after every formerly consuming operation; retained
views after copying, invalidation before whole-molecule mutation including rollback/
no-op, extraction failure before execution, ordinary setter visibility, and
existing invalid-access exceptions. Exercise nested views, not just entity roots.
Retain RingViews consumption and dependent-accessor invalidation checks.

#### S11b — Python guarantees and user-facing migration

**Modules:** Python API guide, binding docs/examples, Python integration tests,
and affected Rust data-type/nomenclature guides. **Documentation/tests, green.**
[dep: S0d, S0g, S10a, S11a, S6a]

Document source/result/receiver ownership for each changed operation, terminal
editor finish, explicit Edits copy, rollback selection, and view invalidation.
Remove consumed-Molecule examples and exception claims without deleting the
remaining Edits, editor, and RingViews lifecycle guarantees. Justify Python
defaults and copying against the selected Rust operations. Preserve current-state
documentation rather than narrating the migration. Value-trait ownership and the general operation
ownership rules are documented through [236](236-operation-ownership-copy-costs-2026-10-09.md)'s S6a.

#### S11c — Align Rust-named Python payloads and immutable observations

**Modules:** umol-py value types and getters, their conversion callers and tests,
and binding documentation. **Breaking getter/name changes, red→green with caller
migration in this item.** [dep: S2a–S2h, S5a, S6c, S7b]

**Settled scope:** this is a bounded binding correction, not an instruction to
redesign or certify the entire Python object graph. The changes below are already
authorized. Their implementation must not restart the nomenclature, immutability,
or ownership discussion, or become a prerequisite for earlier subitems.

**Names:** use the named fields in the corresponding Rust enum declarations:

| Python payload | Required names |
| --- | --- |
| All 31 RelationalConstraint variants | The exact fields of the matching Rust variant in ir/constraint/relational.rs, including bond/atoms, system/atoms, predicate/predicates, dative/parallel, and stereo_atom/stereo_bond as declared there. |
| MoleculeConstraint.ChargeSum | atoms, sum |
| MoleculeConstraint.BondOrderSum | bonds, sum |
| MoleculeConstraint.Connected | atoms |
| StereoConfigurationUpdate.Kinded | kind, coset |

Already named fields, including MoleculeConstraint.UnpairedElectronCoupling,
keep their Rust names. Align getters, constructor keywords, repr, and matching
fields for the listed variants, retaining positional construction. Genuine Rust
tuple payloads keep _0/_1/_2. Do not add value/members/terms/factors or other names
without an existing Rust counterpart; do not rename Rust to manufacture one.
Every payload remains accessible on Python 3.9 through properties, without match.

**Getter results:** the collection properties in the review's "Confirmed getter
deviations" table have these concrete replacements. The table's named families
are the scope boundary, not an invitation to trace every reachable object.

| Current observation | Required result |
| --- | --- |
| Set payloads in ElementForm, IsotopeMassForm, NumForm, PredExpr, StereoTerm, StereoCoset, TopicityRelationForm, and StereogenicityForm | frozenset, including sets inside optional restrictions. |
| Sequence payloads in ArithExpr, PredExpr, and ElectronCountsForm | tuple; child forms retain their existing classes and value semantics. |
| RelationalConstraint and MoleculeConstraint participant/predicate sequences | tuple; preserve None where present and the existing child constraint classes. |
| Edit and Delta collection properties, including S2 batch rows | tuple for outer sequences, rows, and nested donor/atom/ligand sequences. Retain the existing read-only entity payloads. |
| Correspondence.matched_pairs/left_unmatched/right_unmatched, Remapping.images, Compaction.removed | tuple, with tuple pair entries. |
| ValenceEntry.target_covalences/aromatic_valences/fallback_aromatic_valences, ElementScope.AllowList.elements, ResolveReport.tie_breaks | tuple. |
| StereoModel.kind_models | types.MappingProxyType over the getter's independently owned dictionary; preserve the existing model value classes. |
| HashedFeatureSet.ids, CountedHashedFeatureSet.entries, StructuralFeatureSet.keys, SignedHashedFeatureSet.entries, RoleTaggedHashedFeatureSet.ids | tuple; retain tuple records and immutable bytes keys. |
| DativeBondView.donors, StereoAtomView.ligands, StereoBondView.ligands | tuple; leave the views' setters and owner relationships unchanged. |

Preserve contents, sequence order, duplicates, optionality, and existing copying
semantics. Constructors continue accepting their current list/set/dict inputs
and accept the corresponding immutable getter results. These are output-container
changes, not new domain types, new mutable leaf operations, or new public wrapper
classes. Fixing a shared leaf getter once also fixes access through fields,
constraints, and edits; it does not require redesigning those enclosing types.

For the eight entity Update.constraints getters, return the existing corresponding
ConstraintsView retaining the Update owner. Extend its internal storage choice;
reads borrow that Update's constraints, and all writes, including indexed ring
count writes, raise TypeError. No new public constraint type is introduced.
Ordinary ConstraintsForm construction remains writable. Existing entity forms
and their live constraint views retain their current mutation behavior.

**Explicit exclusions:** ReactionLayout.lhs/rhs, ReactionMetadata.lhs, Reaction
ownership, Solution.Determined.molecule, RingView.id, registry/completion read
methods, and explicit copy/conversion/operation-result methods. Their audit
findings do not extend this item. Do not freeze Molecule, ordinary entity forms,
live writable views, Edits/Deltas containers, or independently mutable owners.
Do not replace their setters with immutable-value reconstruction. Do not add
public frozen/mutable type pairs, new lifetime protocols, or a general wrapper
framework. Further API-wide review is not a completion requirement.

**Checks:** exact values and Rust roundtrips; the named properties above; access
without match; exact tuple/frozenset/mapping-proxy result types; rejected container
writes and Update.constraints mutations; constructor-input independence and
construction from getter results; unchanged live mutation and explicit copies.
Use each listed getter family and the existing nested paths to it. Successful
completion means these changes and checks, not a claim that every Python object
is transitively immutable. Migrate affected callers and tests in this item and
run focused Rust/Python suites with a rebuilt extension. S12b remains the full
workspace gate. S2a–S2h do not wait for this work.

### S12 — Whole-operation comparison and final verification

#### S12a — Equivalent-contract operation measurements

**Modules:** existing graph-IR/graph benchmarks, isolated Python measurements,
this record. **Verification, green.** [dep: S0–S7, S10–S11]

Complete section 4 coverage after production changes: resolution, projection,
all three transformers, and reaction application with retained independent
products. Compare against S0a baselines and fixtures; no repeated boxing-layout search.
Include sparse/extensive edits, small/larger molecules, unique/shared storage,
bounded failures, and plain/tracked execution. Include structural additions,
removals/cascades and participant-list changes where the earlier charge-only
comparison was insufficient. Most cases have no molecule-level constraints.

Measure whole-operation time, allocation traffic, peak extra memory, Undo count/
capacity/payload sizes, and separate correspondence storage. Count clone/installation,
publication checks, rollback/discard, and result destruction consistently. Include
actual Python calls. For composite operations and prepared transactions,
no-recovery execution is an overhead reference rather than a recovery-equivalent
competitor. Compare Clone/Journal only at those boundaries; individual resolvers
and transformers use their implemented complete-plan paths. Reaction application
keeps its source/product contract; do not insert a recovery journal into product
construction simply to create a comparison.

**Checks:** identical successful results and agreed recovery guarantees before
interpreting timings. Record which workloads favor which selector without an
automatic policy. Do not infer leaf-form or Reaction binding redesigns from these
measurements; those remain assigned to 231.

#### S12b — Final integration and closeout

**Modules:** workspace, Python package, guides, discussion/status records,
whitepaper/main.tex (read-only review), whitepaper/listings/, and whitepaper/tests/.
**Verification/documentation, green.** [dep: S0–S7, S10–S12a]

Audit the final public symbols, production callers, examples, and lifecycle docs
against this plan. Run the full workspace tests/lints, required feature-gated
suites, rebuilt Python suite, formatting, and the final Rust 1.87 gate once.

Run the manuscript's extracted-listing harness against the final library:
from the manuscript repository, run `python tests/run.py /path/to/umol` using
an environment with umol and pytest installed. Check inline API references for
names, signatures, ownership, and results. Report discrepancies to the author;
manuscript prose is author-owned and must not be edited by agents. Listing
corrections require authorization. Keep test scaffolding outside the printed
code; do not expand schematic examples or assume a shared namespace between
independent listings.

Preserve measurements and remaining review items in durable documents before
removing scratch artifacts. Close 181 only when its agreed work is complete;
then start the Reaction ownership work already scheduled in 231.

**Dependencies:** S0a precedes every production change;
S0b → S0c2 and S0b → S0d; S0e → S0e1 → S0e2 → S0f → S0g;
then S1a/S1b → S2a–S2h → S3a/S3b → S4a;
independently, S5a → S5b, S6a/S6b/S6c, and S7a → S7b → S7c.
S11c follows S2a–S2h, S5a, S6c, and S7b, before final S12 verification;
its getter and naming changes do not block the batching work.
S1a/S1b and S2a–S2f are independent within their stages. The stereo sequence is
S2g1 → S2g2 → S2g3 → S2h; S2g3 also requires the shared stereo-key syntax to be
settled. S2g1/S2g2 can proceed independently of S2a–S2f. S2i follows S2a–S2h,
and S2j follows S2i. S5 needs no boxing, batching, or molecule
lifecycle change and is not deferred. S0 assigns lifecycle/recovery; S10–S12
assign operation, Python, and final-measurement work. S8/S9 material has moved
to S0; later labels are retained. Their dependencies are listed per item. The
migration interface review retains the unresolved shared stereo-key syntax.
