# 181 — Python boundary ownership and copy costs

Status: Proposed
Date: 2026-08-03
Relates: [179](179-python-editing-and-transactions-2026-08-02.md),
[151](151-python-molecule-workflows-2026-07-13.md),
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

### Transaction interface and rollback selection — 2026-10-08

Keep Transaction::run, apply, probe, commit, tracked_commit, and rollback, and
Molecule transact/tracked_transact. The current implementation records and
replays undo. The proposed interface provides both implementations through an
explicit selector:

```rust
pub enum RollbackMethod {
    Clone,
    Journal,
}
```

Transaction::run takes RollbackMethod after the molecule argument and before the
callback. Selection applies to the entire transaction, not to each batch or phase.
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

The clone-based transaction holds its working molecule in
TransactionStorage::Clone(Option<Molecule>) and applies batches through the
molecule's edit application methods. It does not contain an editor. Transaction
exposes batches and checked observation, not direct mutable entity views.
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

Both recovery methods maintain one correspondence from transaction entry to the
accepted result, independently of Undo. Initialize it to identity at run entry
and update it through ApplicationState after each successful batch. tracked_commit
returns it; plain commit discards it. S0e specifies the implementation and cost checks.

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

ConsumedError currently applies to Molecule, MoleculeEditor, and Edits.

Python Molecule will hold GraphIrMolecule directly, without Option. Migrate its
current consuming methods—edit, apply, tracked_apply, resolve, and
resolve_with_report—to the proposed source-preserving operations. Remove
Molecule::take, molecule-consumed errors, and the fallibility of internal accessors
whose only failure is a consumed molecule. No empty placeholder is installed.
The destructive Rust editor entry point is not required on the Python surface.

Edits retains consumption on transfer and gains an explicit copy method.
MoleculeEditor retains consumption at finish; application mutably borrows it.

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
pub fn run<T, E>(
    molecule: &mut Molecule,
    rollback_method: RollbackMethod,
    f: impl for<'scope> FnOnce(Transaction<'scope>) -> Result<T, E>,
) -> Result<T, E>
where E: From<MoleculeApplyError>;

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

### Correspondence without Undo payload retention

ApplicationState already exists for every apply call: it resolves initial and
new handles within that batch. Its update_correspondence method can incorporate
the batch's removals and surviving additions. Reuse that behavior rather than
reconstructing mappings from old attribute and constraint payloads.

The selected candidate maintains one MoleculeCorrespondence directly in both
recovery modes. Initialize it to identity at transaction entry and update it after
each successful batch through ApplicationState::update_correspondence. Drop each
batch's ApplicationState afterward. tracked_commit returns the accumulated value;
plain commit discards it. Neither mode retains batch mapping history or uses Undo
to reconstruct correspondence.

This makes ordinary transactions pay for correspondence construction and updates.
Whether that cost is significant has not been established. Measure time and memory
for plain and tracked transactions before introducing delayed construction,
structural records, or another tracking-selection interface. S0e measures the
Journal path; S0f measures Clone when that implementation is available.

### Transaction guard

The approved recovery storage is a private enum in molecule/transact:

```rust
enum TransactionStorage {
    Journal(Vec<Undo>),
    Clone(Option<Molecule>),
}
```

Transaction::run owns TransactionGuard, which retains the receiver's mutable
borrow, TransactionStorage, the accumulated correspondence and TransactionStatus.
The callback's Transaction borrows those fields; it does not own recovery.
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

### Outstanding review before implementation

- Higher-level lifecycle/naming is settled in [236](236-operation-ownership-copy-costs-2026-10-09.md). The Python
  resolution rollback default is Clone; the prepared-transaction default is
  Journal.

S1's carrier names, enum declarations, and 88/104 layouts are settled. S2's
Rust and Python interfaces, field/key grouping, prefix recovery, and DSL syntax
are settled, including singular entity keys, row shapes, and parameterized keys.

The crate-private owned resolver/projector functions and leaving unrelated value
traits unchanged follow the scoped design; they do not require a new mechanism
or an independent lifecycle decision. No further semantic gap is identified in
tracking or the two recovery implementations. This is a design-review inventory,
not a claim that the unimplemented migration has passed its verification gates.

## Implementation plan

### S0 — Molecule, editor, and transaction migration

[236](236-operation-ownership-copy-costs-2026-10-09.md) closes before this plan begins.
Both plans together are the 0.10 release blocker.
After 236, the resolvers, transformers, and AromaticityPerceiver no
longer use transactions, so they need no migration here.

The baseline precedes all production changes. The six migration groups follow.
The migration interface review distinguishes settled interfaces from remaining
proposals; an unresolved interface must be specified before its item starts.
Each implementation subitem includes its callers and focused tests; none is a
separate task merely to write another plan.

#### S0a — Benchmark baseline

**Modules:** existing graph-IR editor/reaction benches, graph resolve benches,
bounded operation/Python measurements, and this record. **Measurement, green.**
[dep: 236 complete; precedes every production-change subitem here]

Run the existing editor, reaction, and resolve benchmark targets before changing
production code. Record the source revision, toolchain, build settings, fixture
sizes, ownership mode, and measurement boundaries with the results in section 4.
Retain a named Criterion baseline for subsequent comparisons; raw outputs may live
in scratch but are not the only record of the results.

| Existing target | Baseline cases |
| --- | --- |
| graph-IR `editor` | sparse8/dense80, unique/shared, direct/apply/transact, and six overlay-removal kinds in stored/reordered frames |
| graph-IR `reaction` | supplied-match and enumerated application, tracked counterparts, editor/application cases, combine/split/meet_pushout, and reversal |
| graph `resolve` | ingestion, resolution, individual aromaticity/stereo projection, composite projection, and sparse8/dense88 ownership cases |

Reuse the bounded recovery study's fixtures and harness for all three transformers,
clone/install versus journal, bounded failure, and separate allocation/peak-memory
measurements absent from Criterion. Verify its source dependencies against the
current production revision before running; do not rerun the boxing-layout matrix.
Use the existing Python-study harness for actual Molecule edit/apply/tracked_apply,
prepared transactions, and resolve calls, with default and supplied models; add
missing call cases to that harness rather than another benchmark framework.
Most fixtures have no molecule-level constraints.

Check result and recovery agreement outside timing. For each row specify whether
input/Edits preparation and final result destruction are timed; keep those
boundaries for the final comparison. Include any clone required by the measured
source-preservation contract inside the operation boundary. Record current
consuming calls separately from source-preserving calls: a changed method name
must not silently change what a benchmark measures. Python uses the repository's
3.13 environment and a freshly built extension. Run measurements sequentially.

**Completion:** actual timing and memory results, commands, and fixture definitions
are recorded durably here, with reusable inputs for S4a/S12a. This item produces
the pre-change baseline; S12a performs the final comparison.

#### S0b — Mutable editor batch application

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
| graph composite Resolver resolve/resolve_with_report/project in `ops/resolve.rs` | Remove every `editor = editor.apply(...)` reassignment, retain phase order, probe locations, error mappings, and final finish. |
| graph individual resolvers: aromaticity/isotope/stereo resolve and project; bonds/multicenter resolve | Split edit/apply/finish chains into a mutable editor, application, then finish. Preserve current public receivers and result types. |
| graph Kekulizer/Aromatizer/DelocalizeCharge transform | Split editor chains likewise; retain existing plan assertions and Kekulizer's post-publication validation. |
| Python `src/transaction.rs` | Use Option::as_mut for apply/tracked_apply; return None/correspondence. Only finish uses take. Reject an unavailable editor before consuming Edits. |

Migrate Rust editor application tests in `ir/edit.rs`, `ir/molecule/tests.rs`,
`editor/transact.rs`, transaction tests, `tests/property/edit.rs`, and
`tests/property/molecule/compaction.rs`, plus affected resolver tests. Replace returned-editor destructuring/chaining with inspection of
the same editor. Test retained partial state on error; do not weaken publication
or old-state precondition assertions. Python `test_transaction.py` and
`test_edit.py`, and binding unit tests, must distinguish terminal finish from
nonterminal apply, including Edits-transfer failure. These changes cannot wait
for S11.

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

#### S0c1 — Dropped

into_editor is not added; owned callers keep `edit(self)`. See [236](236-operation-ownership-copy-costs-2026-10-09.md).

#### S0c2 — Python Molecule editing without consumption

**Modules:** Python Molecule edit/apply, binding tests, and the Python guide.
**Python behavior change; Rust signatures unchanged.** [dep: S0b]

Python edit, apply, and tracked_apply borrow the source through to_rust, clone
it, and call the consuming Rust method. They never take the molecule or advance
its counter. Transfer Edits and return the new editor, molecule, or
correspondence.

Python Molecule still has Option storage because resolve/resolve_with_report
consume until S10a. Keep to_rust fallible and retain take for those consumers;
S11a removes that storage only after the last consuming binding changes. Update
binding unit tests and `test_transaction.py`, `test_edit.py`, `test_molecule.py`
and `test_molecule_view.py` where consumption/invalidation was caused by these
three methods. Remove redundant copy-before-edit/apply in Python
`test_resolve.py` and `test_workflow.py`; retain the source-preservation assertions.
Keep consumption tests for resolution until its migration. Update the Python
guide's Molecule editing paragraph.

**Checks:** source preservation on success/failure, independent results, retained
source views, and correct tracked results. Rust edit/apply/tracked_apply,
resolver/transformer signatures, and ingestion/export calls remain unchanged.

#### S0d — Scoped in-place editing conveniences

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

#### S0e — Transaction correspondence independent of Undo

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

#### S0f — Clone and Journal transaction recovery

**Modules:** graph-IR transaction, Rust Transaction::run callers, and tests.
**Breaking run signature, red→green within the subitem.** [dep: S0e]

Add RollbackMethod::{Clone, Journal} and the agreed selector position on run:
receiver, rollback method, callback. Clone applies to one working molecule and
installs only after successful commit plus callback return. Journal retains
undo-based restoration. Keep apply/probe/commit/tracked_commit/rollback signatures,
separate batch namespaces, and abort behavior. Both modes use S0e's tracking.

Migrate run callers together. Existing transact/tracked_transact conveniences
retain their signatures at this boundary and explicitly pass Journal internally;
S0g exposes the selector on those conveniences. Higher-level operation callers
also pass Journal until their own migration selects otherwise.

**run caller migration in this subitem:** Molecule transact/tracked_transact in
`ir/molecule/transact.rs`, composite Resolver resolve_into/resolve_into_with_report/
project_into pass Journal in the new second argument.
Their public signatures stay unchanged here. Update graph-IR molecule/transaction
unit tests, edit properties, and run's rustdoc examples at the same boundary.
Exercise both methods in transaction recovery tests; retain Journal explicitly
where a test specifically exercises undo restoration. No run caller waits for S10.

**Guard storage:** retain the receiver borrow, correspondence and status on
TransactionGuard; replace its journal field with the approved TransactionStorage
enum above. Transaction borrows those fields. run initializes Journal with an
empty Vec or Clone with Some(receiver.clone()). Keep restoration owned by run;
no empty placeholder is installed in the receiver. Application dispatch selects
apply_edit_with_undo for Journal and apply_edit for Clone. Both feed the same
per-batch ApplicationState and accumulated correspondence.

Probe selects the working copy while Clone is active and the receiver after
abort; Journal always probes the receiver. Commit checks that same active
molecule. Accept Clone only after commit plus callback success, moving its value
from the Option into the receiver. Rollback and guard Drop replay Journal undo
or discard the Clone value. Retain the existing status transitions and scoped
callback lifetime. Measure S0e's direct-correspondence cost for Clone plain/tracked
execution as well. Export RollbackMethod beside Transaction through
molecule and ir, following the existing transaction-type exports.

**Checks:** both modes on success, first/later application failure, failed commit,
explicit rollback, callback error after commit, unwinding between completed edits,
forgotten handles, and tracked structural changes. Probe after abort observes the
restored entry value. Do not add recovery from allocation failure or internal
partial-Edit panics.

#### S0g — Prepared-batch transaction entry points

**Modules:** graph-IR Molecule transact/tracked_transact; Python molecule,
rollback-selector binding/exports, callers, and tests.
**Breaking selector interfaces, red→green within the subitem.** [dep: S0f]

Expose RollbackMethod through Molecule transact/tracked_transact and Python
prepared-batch calls. Delegate to Transaction::run, apply each batch immediately,
and commit once. Preserve the existing whole-operation correspondence result.
Python continues to accept arbitrary iterables, extract the batch objects before
receiver mutation, and transfer Edits sequentially without prechecking the entire
list for consumed or repeated objects. Failure transferring a later batch may
leave earlier batches consumed, but leaves the molecule and its view counter
unchanged. Before execution, advance the molecule counter under the existing rule.

**Interface status:** Rust takes batches followed by rollback_method. Python
exposes RollbackMethod and a keyword-only rollback_method defaulting to Journal;
Clone remains selectable. These signatures and the default are settled. Scoped
Python Transaction handles/callbacks are not introduced.

**Prepared-batch caller migration in this subitem:**

| Caller | Required change |
| --- | --- |
| Python `src/molecule.rs` transact/tracked_transact | Accept the keyword-only RollbackMethod selector and forward it, default Journal; retain current extraction, sequential Edits transfer and counter timing |
| graph-IR editor/reaction benchmarks and molecule/transaction/property tests; graph resolver tests | Supply the explicit selector. Keep baseline journal rows on Journal; add Clone comparison cases without relabeling old results. |

Update Python RollbackMethod exports with its binding and both method signatures.
Extend prepared-batch Python tests for explicit Clone/Journal; calls omitting the
keyword continue to exercise the selected Journal default. Update transaction
rustdoc and the Python guide here; S11b is the final consistency pass, not the
first documentation update.

**Checks:** Clone/Journal success and failure parity, plain/tracked calls,
arbitrary iterables, repeated/consumed batch inputs, independent batch namespaces,
source restoration, and unchanged view invalidation timing. Update the operation
docs with the selected Python default and exact receiver/input guarantees.

**Compilation and test boundaries:** S0b, S0c2, S0f and S0g include their
binding/caller changes in the same subitem. Each ends with affected libraries,
tests and benches compiling and focused tests passing; no compatibility wrapper
or disabled test bridges the interval to S10. Check feature-gated property callers
explicitly. At the S0 stage gate, check workspace targets and run the affected
publication, edit, transaction, resolver/transformer and Python lifecycle suites.
Use the Python 3.13 environment for every PyO3 build. Full workspace test/lint and
MSRV runs remain at S12b; baseline timings are not rerun after every subitem.

**S0 gate:** all listed interface decisions are resolved where marked, affected
Rust/Python callers compile and pass their focused tests, and lifecycle/recovery
contracts are verified before the boxing and batching stages. Full workspace and
MSRV verification remains at the final combined gate.

### S1 — Box the selected large payloads

#### S1a — Undo: 88-byte layout

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

#### S1b — Edit: 104-byte layout

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

### S2 — Field/key batches and corresponding undos

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

**Native DSL — operation list and row shapes approved:** use singular entity
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
DSL. Stereo rows additionally carry `:kind <kind-or-nil>` in the expect/update
map, used by both constraint operands. No kind is inferred from the molecule.

The key encoding is explicit even for empty batches. Unit keys use their existing
constraint keyword (for example :valence, :aromatic, :electron-count,
:intramolecular, :stereogenicity). Parameterized keys use:

- RingMembership(All): `{:ring-membership {}}`;
  RingMembership(Size(n)): `{:ring-membership {:size n}}`.
- LigandSymmetry: `{:ligand-symmetry {:image [p0 p1 ...]
  :orientation :proper}}` (or :improper).
- Fluxionality: `{:fluxionality {:image [p0 p1 ...]}}`.
- Topicity: `{:topicity {:pair [i j]}}`.

The complete permutation image retains its degree without choosing a stereo kind
for the shared key; constraint values retain their existing kind-based syntax.
Reject malformed key encodings at deserialization. Key/value agreement remains
an application precondition, just as independently supplied old/new values do.

Extend EditsDsl's existing trait boundary and per-entity EditInput organization.
Provide tree readers and tree rendering for the new forms, with
read_edn_/render_ naming where free functions are needed. Do not add a
public batch DSL type or a generic parser that obscures the per-entity grammar.
Valid batches preserve their variant, key, row order, empty rows, and stereo kind
context through DSL and Python conversion. The study's expansion into individual
edits is not a production conversion. No Delta variants are added.

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

All S2 subitems use the following execution contract. Each subitem includes its
own vocabulary, Edits construction, application, undo replay, DSL/Python conversion,
and tests. Entity application loops are independent; S2h reuses the stereo
configuration/key grammar introduced in S2g.

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
no panics from manipulated undo payloads. Include DSL/Python construction and conversion checks with each entity,
preserving typed batches rather than expanding them into individual edits.
Property tests exercise sequence equivalence and restoration independently of
the batching implementation. Invalid histories need not produce correct values.

#### S2a — Atom field and constraint batches

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

#### S2b — Bond field and constraint batches

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

#### S2g — Stereo atom field and constraint batches

**Modules:** graph-IR ir/edit, molecule/apply, dsl/edit; Python edit; their tests.
**Additive vocabulary with exhaustive-match/caller migration; green at completion.**
[dep: S1a, S1b]

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

In dsl/stereo, reuse stereo_configuration and fmt_stereo_configuration to expose
crate-private read_edn_stereo_configuration, read_stereo_configuration, and
render_stereo_configuration for the two edit grammars. The readers require the
complete configuration head and reject trailing constraint text. These are
boundary functions for EditsDsl, not a new public configuration type. Implement
the shared stereo-key grammar here as well; the atom/bond batch readers supply
their corresponding key type.

**Checks:** shared S2 cases for configuration and stereo-atom constraint keys;
include different kinds within a keyed batch, per-entry DSL kind context,
configuration changes followed by rejection, and retained site/ligand order.

#### S2h — Stereo bond field and constraint batches

**Modules:** graph-IR ir/edit, molecule/apply, dsl/edit; Python edit; their tests.
**Additive vocabulary with exhaustive-match/caller migration; green at completion.**
[dep: S1a, S1b, S2g (shared stereo DSL grammar)]

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

**Checks:** shared S2 cases for configuration and stereo-bond constraint keys;
include new and removed handles, per-entry geometry context, sequential
configuration changes, failure after earlier batches, and retained site/ligand order.

**S2 gate:** all eight entity kinds have the same applicable semantics, all
workspace consumers compile and pass their relevant tests, and the new variants
do not defeat the selected compact enum layouts. If a payload would widen an
enum, settle its representation before proceeding; do not silently undo S1.

#### S2i — Edits DSL fuzz target

**Modules:** umol-graph-ir fuzz. **Additive, green.** [dep: S2a–S2h]

Add `fuzz_edits` beside the existing targets: parsing arbitrary input as
EditsDsl must not panic, and a parsed value converts to Edits without
panicking. EditsDsl has one parser, the tree path behind from_edn_str, so there
is no streaming parity to check. Register it in the fuzz Cargo.toml.

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
memory, and Edit/Undo counts/capacities. Include no-recovery and clone/install
paths alongside journals, plus failures and structural edits that incur boxes.
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

#### S6c — Move owned entity inverse results into Python forms

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

// Existing Python method remains unchanged.
fn inverse(&self, py: Python<'_>) -> PyResult<Py<Self>>;
```

Cover AtomDelta, BondDelta, DativeBondDelta, AromaticSystemDelta,
MulticenterBondDelta, NoncovalentBondDelta, StereoAtomDelta, and StereoBondDelta.
Continue calling Rust inverse. Match its owned result by value and move Add/Remove
forms into the existing read-only Python form representation; move owned atom,
ligand, and site lists through their existing conversions as applicable.
Do not borrow the owned inverse merely to clone its complete form again.

Migrate genuinely borrowed producers, including indexing, iteration, and aggregate
Delta conversion, to make their required clone explicit before this boundary.
Migrate all Readonly from_rust callers, including Edit conversion, in the same
subitem. Their required independent copies stay in place; no supplied Python form
is consumed. Keep existing field-change and constraint conversion semantics; this
bounded item does not claim to eliminate their nested conversion copies or change
Rust inverse. No normalization or integrity checks are introduced.

**Checks:** exact inverse values and unchanged inputs for every entity kind and
variant; retained read-only form behavior; exact borrowed-entry conversions and
independence after source/container mutation. Exercise both entity inverse and
aggregate Delta inverse. Repeat the measured AtomDelta Add/Remove cases and add
allocation comparisons for overlay forms with owned lists. Confirm removal of
the second whole-form copy without claiming every nested allocation disappears.

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
| Composite _into methods | Existing signatures, one Journal transaction | Expose the approved selector and forward it through that transaction |
| Individual resolvers | Consuming ordinary methods, journaled _into methods | No change in 181; see 236 |
| Transformer | Consuming transform; transform_iter clones once; _into uses Journal | No change in 181; see 236 |
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

Update Python resolve/resolve_with_report to clone the source, consume the clone,
and preserve the source's views. The two _into methods retain their receiver and invalidation contracts and
expose the approved rollback selector. Keep ChemistryModel borrowing from S5 if
already implemented; do not restore extraction copies.

**Checks:** complete determined results, underdetermination, contradictions,
execution failures, optional report equivalence, late projection failure,
ingest/export fixtures, and unique/shared inputs. Count copies at call boundaries
rather than inferring their number from API names.

#### S10b — Dropped

The individual resolvers keep their consuming ordinary methods. Their in-place
_into methods are planned in [236](236-operation-ownership-copy-costs-2026-10-09.md) (S3b, S3e).

#### S10c — Dropped

transform keeps consuming and transform_iter keeps its per-result clone.
transform_into's in-place application is planned in [236](236-operation-ownership-copy-costs-2026-10-09.md) (S3c, S3d).

### S11 — Python molecule storage and guarantees

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
operation errors remain. Keep ConsumedError for Edits and terminal editor finish.
Keep all *Storage types and counter terminology; add no replacement lifecycle
wrapper, weak-reference list, or cache.

**Checks:** source access after every formerly consuming operation; retained
views after copying, invalidation before whole-molecule mutation including rollback/
no-op, extraction failure before execution, ordinary setter visibility, and
existing invalid-access exceptions. Exercise nested views, not just entity roots.

#### S11b — Python guarantees and user-facing migration

**Modules:** Python API guide, binding docs/examples, Python integration tests,
and affected Rust data-type/nomenclature guides. **Documentation/tests, green.**
[dep: S0d, S0g, S10a, S11a, S6a]

Document source/result/receiver ownership for each changed operation, terminal
editor finish, explicit Edits copy, rollback selection, and view invalidation.
Remove consumed-Molecule examples and exception claims without deleting the
remaining Edits/editor lifecycle guarantees. Justify Python defaults and copying
against the selected Rust operations. Preserve current-state documentation rather
than narrating the migration. Value-trait ownership and the general operation
ownership rules are documented through [236](236-operation-ownership-copy-costs-2026-10-09.md)'s S5a.

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
actual Python calls. No-recovery is an overhead reference, not a recovery-equivalent
competitor. Reaction application keeps its source/product contract; do not insert
a recovery journal into product construction simply to create a comparison.

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
S0b → S0c2 and S0b → S0d; S0e → S0f → S0g;
then S1a/S1b → S2a–S2h → S3a/S3b → S4a;
independently, S5a → S5b, S6a/S6b/S6c, and S7a → S7b → S7c.
S1a/S1b and S2a–S2g are independent within their stages; S2h follows S2g
once their prerequisites are met, S2i follows S2a–S2h, and S2j follows S2i. S5 needs no boxing, batching, or molecule
lifecycle change and is not deferred. S0 assigns lifecycle/recovery; S10–S12
assign operation, Python, and final-measurement work. S8/S9 material has moved
to S0; later labels are retained. Their dependencies are listed per item. The
migration interface review has no open interface points.
