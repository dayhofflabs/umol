# 213 — Molecule and reaction mutation

Status: Proposed
Date: 2026-08-27
Relates: [117](117-entity-model-extensibility-2026-06-20.md),
[166](166-molecule-ops-2026-07-27.md),
[211](211-relation-frames-and-api-2026-08-26.md),
[214](214-aggregate-frame-semantics-2026-08-28.md),
[data-type guide](../docs/development/data-types.md),
[nomenclature guide](../docs/development/nomenclature.md)

## Scope — revised 2026-09-18

Review and design the full molecule/reaction mutation API: direct changes, checked
batches, transactions and rollback, participant replacement, publication, and
identity tracking. Include the six overlay entity sets owned by Molecule and the
Rust/Python boundary. The original editor-wrapper consolidation is one part of
this work, not its organizing assumption.

The relation-storage work in [166](166-molecule-ops-2026-07-27.md), S0–S5, is
complete. Its S6–S8 graph-core restoration extension is also complete as of
2026-09-19; the storage prerequisite for restoration delegation is available.
Design and implement the common molecule/reaction mutation infrastructure here,
then return to 166 for
non-transactional editing at operation call sites and
HydrogenFolder/HydrogenUnfolder. Projection requirements inform this design;
their remaining operation work stays in 166. This is a current-state review and
open design, not an implementation plan.

## Current ownership

Molecule owns topology, atom/bond attributes, six overlay entity sets, and global
constraints. Each overlay set already owns an Arc-backed relation set:

| Entity set | Storage | Participants |
| --- | --- | --- |
| AromaticSystems | VarRelationSet | Atoms |
| MulticenterBonds | VarRelationSet | Atoms |
| NoncovalentBonds | FixedRelationSet, arity 2 | Two atoms |
| DativeBonds | FixedVarBirelationSet | Fixed acceptor, variable donors |
| StereoAtoms | FixedVarBirelationSet | Fixed atom site, variable ligands |
| StereoBonds | FixedVarBirelationSet | Fixed bond site, variable ligands |

Their crate-private attribute mutation already uses Arc::make_mut. Their public
surface supplies typed reads and transformations; raw construction, entry
extraction, and mutable attribute access are crate-private.
See [aromatic](../umol-graph-ir/src/ir/aromatic.rs),
[multicenter](../umol-graph-ir/src/ir/multicenter.rs),
[noncovalent](../umol-graph-ir/src/ir/noncovalent.rs),
[dative](../umol-graph-ir/src/ir/dative.rs), and
[stereo](../umol-graph-ir/src/ir/stereo.rs).

MoleculeEditor instead holds three private storage-shape wrappers:
FixedSetStorage, VarSetStorage, and FixedVarSetStorage. Each switches from a shared
relation set to a mutable vector of entries. Attribute mutation and insertion
materialize that vector; publication rebuilds relation storage. These wrappers
also duplicate reads, removal, and compaction behavior. See
[editor storage](../umol-graph-ir/src/ir/molecule/editor.rs).

All five graph-core relation shapes now support participant replacement, add,
remove, and tracked_remove. They own incidence-index maintenance and stable row
compaction. They do not interpret molecular payloads, validate external node/edge
references, or coordinate other entity sets and constraints. The graph itself
already supports add_node, add_edge, remove_cascading, and
tracked_remove_cascading.

The molecule-level wiring to the new relation mutation surface remains to be
designed. The existence of five storage shapes does not require five kinds of
molecular wrapper: the current overlays use only the three shapes above.

## Existing mutation paths

| Surface | Current behavior and boundary |
| --- | --- |
| Molecule direct mutation | Atom, bond, dative, and noncovalent attributes have direct mutable views/bulk methods. Aromatic, multicenter, and stereo attributes and global constraints have public checked modification callbacks. Those checked methods modify a private candidate and publish only after integrity succeeds; they do not create an undo journal. |
| MoleculeBuilder | Fresh construction delegates to MoleculeEditor; build uses the asserted integrity gate. |
| MoleculeEditor primitives | Add atoms, bonds, and all six overlays; mutate attributes and constraints; remove topology with cascading consequences; remove each overlay kind. These are direct operations without an undo journal. Participant slices in mutable views remain read-only. |
| MoleculeEditor::apply | Consumes the editor and applies Edits without constructing Undo. Failure returns no partially modified editor. Success still requires publication. |
| MoleculeEditor::transact | Borrows the editor, applies Edits, and returns a Transaction journal. Application failure triggers rollback; rollback failure has its own error. Successful application does not itself publish a valid Molecule. |
| Molecule::apply | Applies Edits to a private editor and checks publication, preserving the source molecule on failure. Reports batch failures separately from integrity failures. |
| Editor snapshot/build | snapshot and try_build check aggregate integrity; snapshot preserves the editor, try_build consumes it. build is the asserted counterpart. |
| Tracked operations | Return identity correspondences/compactions. They are independent of rollback journals. The editor also maintains an initial-to-current session correspondence during ordinary operations. |

See [Molecule](../umol-graph-ir/src/ir/molecule.rs),
[builder](../umol-graph-ir/src/ir/molecule/build.rs),
[editor](../umol-graph-ir/src/ir/molecule/editor.rs), and
[batch application](../umol-graph-ir/src/ir/molecule/transact.rs).

Edits is an ordered batch with handles for initial entities and entities created
within the batch. Overlay removals carry offered old participants and attributes;
field changes carry old/new values. Both apply and transact resolve handles and
check preconditions. Their dispatch paths are separate, with shared operations
underneath. Neither Edits nor the direct editor currently exposes participant
replacement. Undo records realized changes, including cascades and compactions;
it is not simply an inverse Edit. See [edit carriers](../umol-graph-ir/src/ir/edit.rs).

A Transaction is issued for a particular post-state; its rollback guarantee is
bound to that history. Restoration also has to recover deleted rows at their old
dense positions. Append-only relation insertion does not alone replace this
restoration path. Review that path alongside forward mutation.

Other existing molecule operations include extraction, combination, splitting,
constraint lifting/inlining, remapping, reframing, and canonicalization. Their
named semantics remain relevant consumers/producers; a common mutation foundation
does not by itself justify replacing their public APIs with edit batches.

### Reaction definitions, spans, and application

Reaction owns a Molecule lhs and Deltas. It exposes checked/asserted construction,
read access, into_parts, and named transformations, but no ReactionEditor or
mutable lhs/deltas accessor. Deltas itself is an open mutable carrier; rebuilding
a Reaction establishes aggregate integrity. ReactionSpan likewise has an open
entries carrier and checked/asserted publication. See
[Reaction](../umol-graph-ir/src/ir/reaction.rs),
[Deltas](../umol-graph-ir/src/ir/delta.rs), and
[ReactionSpan](../umol-graph-ir/src/ir/reaction_span.rs).

Editing a reaction definition and applying that reaction to a host are distinct
operations. The former must coordinate lhs ids, delta targets, participant
frames, and constraints; span editing must coordinate both sides and their
correspondence. The latter already has contextual matching and applicability
contracts, including ordinary non-applicability versus errors. Its implementation
lowers to molecule Edits, calls transact, discards the journal, and checks product
integrity. It derives its result correspondence from the host/product relation;
that witness must not be silently replaced by an editor witness with different
pairing semantics.

The review includes both paths. Whether reaction-definition editing needs a
public editor, checked modification methods, or reconstruction through existing
constructors is still open. Delta, Edit, and Undo have different roles; reviewing
them together does not establish that they should share a representation.

### Python ownership

The Python MoleculeEditor wraps an optional Rust editor. apply consumes it,
including on failure; transact borrows it. Build consumes it; snapshot does not.
Python exposes removal methods but not the complete Rust direct-add/mutable-view
surface. Transaction rollback consumes the journal. Review parity deliberately
rather than exposing detached mutable copies or assuming every Rust borrow can
be bound directly. See [bindings](../umol-py/src/transaction.rs).

## What projections expose

AromaticityResolver::project and StereoResolver::project construct field updates
and collect full old overlay entries for removal before calling transact and
discarding the journal. IsotopeResolver::project also uses transact. The aggregate
Resolver::project already works on a private candidate before publishing its
successful result. See [aromaticity](../umol-graph/src/ops/resolve/aromaticity.rs),
[stereo](../umol-graph/src/ops/resolve/stereo.rs),
[isotope](../umol-graph/src/ops/resolve/isotope.rs), and
[resolver](../umol-graph/src/ops/resolve.rs).

Switching transact to apply would remove journal construction, but would retain
the work of describing old state, constructing Edits, and resolving handles.
Direct mutation already exists for some of these operations. The design question
is what the caller must express, and which common operations own the remaining
bookkeeping. Chemical interpretation stays in umol-graph; storage maintenance,
reference compaction, and graph-IR integrity belong below it. These are source-level
observations, not measurements of runtime cost.

## Settled storage and editor ownership — 2026-09-18

Remove FixedSetStorage, VarSetStorage, and FixedVarSetStorage from the design.
MoleculeEditor holds the same storage types as Molecule: Graph, Arc-backed
atom/bond attributes, the six typed overlay entity sets directly, and Constraints.
It also retains its editing-session correspondence. Entity-set mutation delegates
to relation storage, using the sets' existing copy-on-write ownership. There is no
separate mutable-vector representation. Transaction journals remain separate from
the editor's intrinsic state.

Retain MoleculeEditor as a separate owned type. Molecule promises aggregate
representation integrity; the editor may hold unfinished changes whose combined
result must pass the publication gate. For example, replacing stereo participants
and adjusting their attributes may require coordinated changes before the result
is a valid Molecule. The editor can be abandoned without publishing that state.

A mutable Molecule borrow alone does not provide this boundary. When the borrow
ends, including through an early return or panic, the caller again has a Molecule
whose integrity must hold. Such mutation must either preserve integrity at each
operation boundary or use controlled editing that validates the result and
preserves/restores the original on failure.

Molecule and MoleculeEditor therefore use the same storage types with different
aggregate guarantees. Integrity-preserving operations may remain directly on
Molecule; retaining the editor does not require routing every mutation through it.
The exact entity-set mutation methods and their visibility remain to be designed.

## Disposition of the earlier wrapper proposal

The earlier sketch proposed OverlayEditor<O>, an Overlays trait, and a
Shared/Mutable-vector representation. None is implemented. Storage-owned
replacement now makes the mutable-vector representation unsuitable as the
participant-mutation mechanism. The typed entity sets already provide
copy-on-write ownership; a second representation is not needed merely to preserve
sharing.

The sketch also described raw constructors and entry extraction as public; they
are currently crate-private. Publishing them through a trait would broaden the
construction boundary. A generic trait, its name and visibility, constructor
reshaping, and allocating editor views are not approved requirements. The editor
will own the existing typed sets directly; the mutation access they need while
the editor is transient remains to be designed.

The useful part of the earlier proposal is participant-frame ownership. Generic
relation storage knows factors and indices; graph IR knows which participants
carry a frame and how attributes transform with that frame. The editor's existing
six *_equiv methods already align offered participants and transport attributes
for old-state comparison. Moving this knowledge to an entity-set operation is an
ownership proposal, not new comparison semantics.

For pairwise alignment, retain the direction `to[i] = from[action[i]]` and compare
attributes after transport into the stored frame. Ordinary unordered factors use
DynPermutation; stereo factors use bounded Permutation. Stereo-bond alignment
must respect endpoint blocks: reorder within blocks or exchange complete blocks,
not move an individual ligand across the boundary. The location and visibility of
this operation remain open; it does not require a public generic trait.

## Mutation surface under discussion — 2026-09-18

The following records the proposed surface for discussion. It does not settle
method names, signatures, or the remaining integrity and failure contracts.

### Existing methods to reroute

Delegate through the typed entity sets that own copy-on-write relation storage:

| Editor surface | Storage operation and retained coordination |
| --- | --- |
| add_dative_bond, add_aromatic_system, add_multicenter_bond, add_noncovalent_bond, add_stereo_atom, add_stereo_bond | Relation add; editor extends session correspondence. |
| Each overlay remove_* / tracked_remove_* pair | Relation tracked_remove; editor compacts constraints and correspondence using its returned row compaction. |
| Mutable overlay attribute views | Entity-set attribute mutation through relation data_mut. |
| Topology remove / tracked_remove | Existing Graph::tracked_remove_cascading, followed by relation tracked_compact; editor coordinates all entity spaces and constraints. |
| Overlay reads and views | Typed entity-set accessors. |

Plain editor removal still needs compaction internally even when it returns no
witness. Atom/bond addition already delegates to Graph, and entity-set frame
transformations already use storage participant permutation. Undo of additions
can reuse removals; undo of removals needs the restoration capability below.

### Missing editor participant operations

| Entity kind | Mutable participant components |
| --- | --- |
| Aromatic system | Atom sequence |
| Multicenter bond | Atom sequence |
| Noncovalent bond | Fixed atom pair |
| Dative bond | Acceptor and donor sequence, together and separately |
| Stereo atom | Site and ligand sequence, together and separately |
| Stereo bond | Site and ligand sequence, together and separately |

Propose whole replacement, individual replacement, and positioned insertion/removal
for variable sequences. Fixed components do not gain insertion/removal. Molecular
names should identify atoms, donors, ligands, and sites rather than expose factor
numbers. These changes preserve the entity id. Ordinary replacement preserves the
attributes; the draft permits the caller to adjust them separately before publication.
This permits ligand changes without removing and recreating the stereo entity.

Reframing has a separate meaning: transport attributes and affected constraints
along with participant order to preserve the represented value. Reuse the existing
frame machinery; storage permutation alone does not perform that transformation.

Localized-bond endpoint replacement is another gap. Graph has no edge-endpoint
replacement operation, so retaining BondId while rewiring a bond would require
additional graph-storage design beyond the completed relation surface.

### Direct Molecule mutation

Retain the existing direct attribute mutation and checked attribute/constraint
callbacks. Consider direct structural methods for complete additions, removals,
and replacements, with success preserving aggregate integrity and failure leaving
the molecule unchanged. Share underlying mutation operations with the editor;
these methods need neither Edits construction nor an undo journal. A private
candidate followed by checked publication is a possible initial implementation;
local checks may suffice where their preservation guarantee is established.

Removal is not automatically integral. Deleting a bond from a stereo site to an
explicit ligand can preserve all stored atom references while breaking required
ligand incidence. Relation compaction does not detect that missing connecting
bond. A direct Molecule operation must reject the result or have an explicitly
defined broader cascade; this choice remains open.

Some complete replacements must accept coordinated attributes, such as changing
an aromatic participant count and its electron-count vector. The editor can
perform these changes in separate steps; direct Molecule mutation must complete
them within its integrity boundary. How much of the editor's primitive surface
should have Molecule conveniences remains open.

### Removal rollback: current implementation

Transactional removal captures actual stored entries before mutation, including
old entity ids, participants in their stored order, and attributes. Topology
removal also captures cascading bonds and overlays. Undo stores the compaction
and the constraint changes needed for restoration. UndoCompaction is an inverse
view of MoleculeCompaction, not a separate history of payloads.

Rollback replays Undo in reverse order. Before each operation, validate_undo checks
counts, reconstruction positions, and relevant reference coverage. For each
relation family, restore_* then:

1. Allocates slots for the original row count.
2. Places surviving rows at their original ids using inverse row compaction.
3. Expands surviving node/edge references through inverse graph compaction.
4. Places saved removed rows at their original ids. Their participants already
   use the original graph ids and must not be expanded again.
5. Stores the reconstructed entries in the current mutable-vector wrapper.
   Relation storage and incidence indices are rebuilt when materialized for publication.

Topology restoration similarly reconstructs atom/bond arrays and rebuilds Graph
from restored endpoints. The transaction layer restores constraints and expands
the session correspondence after restoring the affected tables. See
[restore methods](../umol-graph-ir/src/ir/molecule/editor.rs),
[undo execution](../umol-graph-ir/src/ir/molecule/transact.rs), and
[UndoCompaction](../umol-graph-ir/src/ir/compact.rs).

### Graph-core restoration dependency

The settled graph-core restoration design and public contracts are
recorded with the earlier storage work in
[166 — Graph and relation restoration](166-molecule-ops-2026-07-27.md#graph-and-relation-restoration).
Completed 2026-09-19: Graph::restore and restore/restore_participants on all five
relation-set shapes are implemented and verified in 166 S6–S8;
uncompact remains non-mutating reference translation. Restoration returns (): matching
removal data recovers the original storage; manipulated inputs must not panic but
have no specified restoration result. This document owns editor
integration, including attribute arrays, constraints, correspondence, and
transaction coordination, and delegates graph/relation reconstruction to those
storage operations. The editor currently retains the reconstruction paths described
above; rewiring them is not part of the completed graph-core work.

## Design questions to settle together

1. **Owning operations.** Which mutations belong on typed entity sets, which need
   molecule-wide coordination, and which published Molecule conveniences should
   remain? Delegate relation changes to storage while keeping constraint updates,
   cascades, and identity tracking with the aggregate that owns those references.
2. **Participants and attributes.** Specify whole replacement, individual
   replacement, addition/removal where meaningful, and frame permutation together.
   Storage replacement preserves the payload; graph IR must define how a changed
   frame or participant count interacts with that payload. Distinguish reframing
   the same entity from changing its participants. Decide when both must be
   supplied together and when temporary inconsistency is permitted in a draft.
3. **Execution and failure.** Define a direct path for callers working on a private
   draft and a transactional path for rollback. Decide which operations they
   share, where preconditions are checked, and what survives failure. Keep
   journal-free execution, atomic publication, and reversible history distinct.
4. **References and witnesses.** Separate current dense ids, batch handles,
   session correspondences, per-operation compactions, and rollback restoration.
   Preserve constraints through compaction. Record when an identity survives
   participant replacement and what tracked results mean for changed frames.
5. **Reaction mutation.** Specify editing lhs/deltas and editing both sides of a
   span, including dependent references and publication. Review reaction
   application as a consumer of the same molecule primitives without changing its
   matching, non-applicability, or correspondence semantics by accident.
6. **Public contract and bindings.** Enumerate the intended constructors,
   conversions, mutators, publication gates, and contextual consumers before
   implementation. State invalid-id panics versus recoverable input errors,
   failure atomicity, and Rust/Python ownership for each surface. Existing low-level
   storage failure rules do not automatically determine the batch or Python rules.

The starting proposal for discussion is storage-owned mutation, typed graph-IR
coordination, and transaction handling layered on the same underlying operations.
Public names, signatures, and any additional type or trait remain to be settled.

## Evidence for the eventual implementation

Preserve the existing laws for handle namespaces, apply/transact result agreement,
rollback, failure recovery, constraint compaction, publication, and correspondence
composition. Relevant suites include
[edit properties](../umol-graph-ir/tests/property/edit.rs),
[publication](../umol-graph-ir/tests/property/molecule/publication.rs),
[compaction](../umol-graph-ir/tests/property/molecule/compaction.rs), and the
[reaction properties](../umol-graph-ir/tests/property/reaction).

Extend coverage from the settled contract: all six entity kinds, participant and
payload coordination, legal and illegal frame alignment, cascade restoration,
reaction reference updates, and Python failure ownership. Include position-sensitive
payloads so that frame tests prove transport. Benchmark representative direct,
apply, and transact paths from the beginning, including copy-on-write sharing and
repeated relation mutation; do not assume the current fixtures establish scale.
