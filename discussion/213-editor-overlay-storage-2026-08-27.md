# 213 — Molecule and reaction mutation

Status: Proposed
Date: 2026-08-27
Relates: [117](117-entity-model-extensibility-2026-06-20.md),
[166](166-molecule-ops-2026-07-27.md),
[211](211-relation-frames-and-api-2026-08-26.md),
[214](214-aggregate-frame-semantics-2026-08-28.md),
[228](228-python-api-parity-2026-09-21.md),
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

Priority update 2026-09-21: pause this work for the urgent Python/Rust semantic
parity review in [228](228-python-api-parity-2026-09-21.md). Preserve the settled
editor/transaction scheme below; independent-batch composition remains unresolved
and is not a prerequisite for that review.

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
As revised on 2026-09-20, it holds neither a session correspondence nor an undo
journal. Entity-set mutation delegates to relation storage, using the sets' existing
copy-on-write ownership. There is no separate mutable-vector representation.
Optional correspondence accumulation belongs to checked batch execution; undo
journals belong to transaction execution.

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

The following records the mutation surface and its settled decisions. Remaining
proposals and open integrity and failure contracts are identified separately.

### Existing editor surface and required delegation changes — 2026-09-19

Retain graph-IR coordination while replacing the editor's storage implementations.
Overlay operations delegate through the typed entity sets that own copy-on-write
relation storage. This inventory records the required changes, not an
implementation plan or a settlement of the execution-method names.

| Editor surface | Current implementation and required change |
| --- | --- |
| add_atom, add_bond | Already delegate to Graph::add_node/add_edge. Retain parallel attribute-array updates; remove session-correspondence extension. |
| add_dative_bond, add_aromatic_system, add_multicenter_bond, add_noncovalent_bond, add_stereo_atom, add_stereo_bond | Replace mutable-entry-vector insertion with relation add through the typed set; remove session-correspondence extension. |
| atom_mut, bond_mut | Already mutate copy-on-write attribute arrays. Retain graph-IR ownership of these arrays. |
| Mutable attribute views for all six overlay kinds | Replace entry-vector materialization with typed-set access to relation payload mutation. Keep attributes and participants accessible in the entity view. |
| push_constraint, constraints_mut, inline constraints through attribute views | Retain graph-IR mutation; graph-core does not interpret molecular constraints. |
| Each overlay remove_* / tracked_remove_* pair | Use relation tracked_remove and its returned row compaction instead of editor-owned row removal and separately constructed compaction. Retain constraint compaction; remove session-correspondence maintenance and the public editor tracked_remove_* variants. |
| Topology remove / tracked_remove | Retain Graph::tracked_remove_cascading; replace wrapper-specific overlay compaction with relation tracked_compact through typed sets. Retain atom/bond attribute compaction, cascade coordination, and molecule-wide constraint updates. Remove session-correspondence maintenance and public editor tracked_remove. |
| Internal undo-addition removal and restore_* methods | Route removal through the same primitives. Replace graph/relation reconstruction with storage restore_participants and restore in the appropriate order. Retain attribute and constraint restoration with the settled local guards below. Failed tracked execution discards its correspondence accumulator. |
| apply, transact, and tracked counterparts | Route both execution paths through the same graph-IR mutation operations. Handle resolution and forward edit preconditions remain batch concerns; undo capture/replay remains transactional. Review names and execution contracts separately. |
| snapshot, try_build, build, and tracked counterparts | Publish from the editor's typed sets without rebuilding relation storage from entry vectors. Remove snapshot/tracked_snapshot and tracked direct publication, and consolidate draft publication; transaction commit includes the integrity gate, as settled below. |
| Overlay reads and views; six internal *_equiv methods | Remove storage-wrapper dispatch. Use typed-set accessors and graph-core participant comparison where applicable; retain graph-IR frame transport and payload comparison. |

Plain editor removal still needs compaction internally even when it returns no
witness. Atom/bond addition already delegates to Graph, and entity-set frame
transformations already use storage participant permutation. Undo of additions
can reuse removals; undo of removals needs the restoration capability below.

### Settled editor participant mutation — 2026-09-19

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

### Whole replacement in edit batches — 2026-09-19

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

### Settled Delta and Undo replacement variants — 2026-09-19

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

### Delta replacement and span representation — 2026-09-19

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

### Settled local editor getters — 2026-09-19

Provide the same local read conveniences on immutable and mutable editor entity
views, using the molecule-backed immutable counterparts' names and meanings. Participant
inspection should not require callers to index slices or unpack attributes merely
because the entity is being edited.

The current editor views expose attributes as fields and only a small subset of
the immutable views' participant getters. Add the following local read surface:

| Entity kind | Local getters |
| --- | --- |
| Atom | element, isotope_mass, charge, implicit_hydrogens, lone_pairs, unpaired_electrons |
| Localized bond | atom_ids, order, charge, unpaired_electrons |
| Dative bond | donor_ids, acceptor_id, donor_count, atom_count, order; retain atom_ids |
| Aromatic system and multicenter bond | atom_count, electrons, electron_count, charge, unpaired_electrons; retain atom_ids |
| Noncovalent bond | atom_ids, kind |
| Stereo atom and stereo bond | site_id, ligand_count, ligand, ligand_position, ligand_frame, ligand-kind filters, ids and counts; configuration accessors |

Each editor view provides id() and attributes(); mutable views additionally
provide attributes_mut(). References returned by mutable-view accessors are tied
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

## Settled rollback contract and guard ownership — 2026-09-19

Rollback follows the same history-bound contract as graph-core restoration:

- Rolling back a transaction from its matching post-state restores the recorded pre-state.
- Mismatched or manipulated history must not cause a panic. The resulting state
  is unspecified; rollback need not detect misuse, report it, or leave the editor
  unchanged.

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
| Correspondence on failure | Discard the failed tracked execution's accumulator; the direct editor no longer has a session correspondence to restore. |

Do not retain blanket checks of all entity counts, equality between forward and
inverse compaction records, or external graph membership of saved overlay
participants. Storage delegation removes the editor's relation-row reconstruction
checks. Attribute and constraint restoration must still protect their own accesses;
they do not need to prove that the supplied history matches the editor.

The consuming implementations must become panic-free before validate_undo is
removed: deleting the checker alone would expose existing indexing, unwrap, and
expect paths. This is a settled design, not implemented behavior. The active
transaction design below replaces the public detached journal and removes its
tracked_rollback interface. Tracked results belong only to batch
application and transaction commit, as settled below.

## Transaction and publication restructuring — decisions 2026-09-19, revised 2026-09-20

Retain direct mutation, journal-free batch application, storage-owned restoration,
and the distinction between an editable draft and an integrity-valid Molecule.
The following lifecycle replaces the current detached-journal design. It is
settled design direction, not implemented behavior or a complete signature list.

Currently, transact rolls back an applied prefix on edit failure and returns a
detached undo journal on success. It does not establish aggregate integrity.
Subsequent build or try_build consumes the editor; publication failure does not
trigger rollback. Callers therefore have to compose execution, acceptance, and
recovery themselves. Transaction describes an undo journal rather than an active
transaction lifetime. This separation is not justified merely by the need for
an integrity boundary.

### Execution and recovery

Keep three execution paths:

- Direct mutation maintains storage invariants but does not check batch
  preconditions or establish aggregate molecule integrity.
- Journal-free batch application checks each edit's preconditions and records no
  undo. It remains useful for already-owned, disposable working state.
- Transactional batch application checks the same preconditions and records undo
  privately so failure can restore the retained receiver.

Use journal-based recovery for transactions. Cloning the receiver and mutating a
candidate was rejected as the transaction strategy: extra allocation then grows
with receiver size, whereas journal storage follows the recorded changes,
including saved removed data. The journal-free path does not authorize replacing
transactional recovery with cloning. Existing copy-on-write storage remains;
this is not a claim that all transaction allocation is independent of receiver
size.

Edit preconditions and final aggregate integrity belong to one transactional
acceptance boundary. Check preconditions before each edit and integrity after the
complete operation. Both kinds of failure restore the state from before the
transaction. Do not check aggregate integrity after each edit or batch:
coordinated changes may require inconsistent intermediate draft states.

### Active transaction and Edits composition

Transaction represents an active operation borrowing the receiver mutably.
Submission queues Edits; execution and private undo recording are deferred until
commit. It is not an open history container or a detached result.

The agreed lifecycle is:

1. Begin a transaction by borrowing the Molecule mutably.
2. Submit one or more Edits batches without execution, validation, or journaling.
3. Commit executes the batches, checks each edit's preconditions, records undo,
   and checks aggregate integrity once at the end.
4. An edit or integrity failure during commit restores the original molecule;
   success accepts the complete result.
5. Explicit rollback or dropping the transaction before commit discards queued
   work and releases the borrow. The molecule has not been changed.

Edits remains incrementally constructible as one sequence before execution.
Independent-batch composition is unresolved, as recorded below.
Submitting batches does not expose an interface for inspecting, rearranging,
appending, or otherwise manipulating the transaction's internal history. Multiple
batches share one commit/rollback boundary. Remove detached-journal return and
composition from the design; retaining undo history after commit is not a
requirement. Every mutation exposed within the transaction must participate in
its recovery guarantee; unrestricted writes that bypass recording cannot silently
form part of that surface.

### Publication and names

A directly edited draft needs one fallible publication operation that produces a
Molecule. A transaction over an existing Molecule commits an integrity-preserving
change and needs no subsequent build call. Keep these responsibilities distinct.
Reserve apply for executing edits, submit for queuing them, and transaction for
the active operation with automatic recovery. The current TransactionError name
is misleading for failures also produced by journal-free application; settle its
replacement with the error surface. Exact publication names and signatures remain
open.

build adds panic-on-error behavior to try_build, not a faster unchecked path.
The pair does not justify retaining two publication operations in the redesign.
Remove snapshot and tracked_snapshot. Their current behavior is to produce an
independent, integrity-checked Molecule while retaining the editor. They neither
compose edits nor provide rollback. They arose in
the [Python editing work](179-python-editing-and-transactions-2026-08-02.md), where
inspection followed by continued editing or rollback was the stated use. No
necessary consumer requirement has been identified for keeping this separate
lifecycle operation; Python exposure alone does not justify it.

The participant mutation methods, position types, local getters, and exclusion of
frame-preserving permutation are settled above. Keep transaction/publication
restructuring separate from that contract.

## Draft editor, transaction, and transformation APIs — 2026-09-20

This is a revised design sketch for review, not an implementation plan. Consuming
journal-free editing and borrowing transactions have different failure outcomes;
the distinction is not an assumed performance advantage of moving over swapping.
The participant mutation and Edit/Undo/Delta contracts above remain settled.
The consuming-edit versus borrowing-transaction distinction and deferred
transaction execution were settled on 2026-09-20; remaining draft names and
transformation result shapes are identified below.

### Editor and journal-free application

Proposed signatures; finish and the error names are draft nomenclature:

```rust
impl Molecule {
    pub fn edit(self) -> MoleculeEditor;
    pub fn apply(self, edits: Edits) -> Result<Molecule, MoleculeApplyError>;
    pub fn transaction(&mut self) -> Transaction<'_>;
}

impl MoleculeEditor {
    pub fn apply(self, edits: Edits) -> Result<Molecule, MoleculeApplyError>;
    pub fn finish(self) -> Result<Molecule, MoleculeIntegrityError>;
}
```

edit moves the molecule's existing storage into an owned editor. Individual
direct mutation methods borrow that editor mutably, including the entity views
settled above. Consumption at editor creation and mutable borrowing for each
mutation are compatible. Empty construction can start from Molecule::new().edit().

Current direct editor mutations and the approved participant methods do not
return Result errors. Invalid ids or positions may panic; temporary aggregate
inconsistency is permitted. finish checks the resulting aggregate after direct
mutation. Direct mutation can also precede apply, whose final integrity check
covers the complete resulting state.

apply checks each edit against the evolving state, executes through the direct
mutation API, and checks aggregate integrity before returning a Molecule. It
records no journal. Any edit or integrity error drops the working state; the
input is consumed and no editor is returned for repair. Molecule::apply delegates
to the same operation without requiring an explicit editor in batch-only code.
finish is the single fallible publication operation for manually edited drafts.
Remove the build/try_build pair and snapshot/tracked_snapshot from the proposed
surface. Applying a batch needs no subsequent finish call.

Retain MoleculeApplyError, with Edit(EditError) and Integrity(MoleculeIntegrityError)
variants. The unqualified ApplyError already names reaction application errors;
do not reuse that name for this surface. EditError retains the applicable handle,
shape, and old-state diagnostics currently named TransactionError; journal
mismatch and rollback-failure variants are removed as already settled. Neither
failure path publishes a partial Molecule.

Returning an editor on error was considered and rejected for parsimony. It adds
a repair lifecycle; when the editor borrows a destination, it also leaves the
problem of abandonment unresolved. A borrowed Molecule must not become empty as
an error sentinel: empty is an ordinary valid molecule. A caller deliberately
preserving a source for consuming execution may clone it explicitly.

This explains the ownership distinction: journal-free failure consumes and drops
the working state; a transaction restores its borrowed receiver. Returning the
restored molecule together with an error from a consuming transaction is possible,
but mutable borrowing restores the caller's existing variable directly, including
on abandonment. Neither ownership mechanism inherently requires copying storage;
their difference is recovery and ergonomics, not memory safety or successful-result
integrity. Cloning or journaling would be needed to promise restoration for the
journal-free path, and neither belongs in that path.

### Borrowing transaction with deferred execution — settled 2026-09-20

```rust
impl<'a> Transaction<'a> {
    pub fn submit(&mut self, edits: Edits);
    pub fn commit(self) -> Result<(), MoleculeApplyError>;
    pub fn rollback(self);
}
```

Molecule::transaction(&mut self) creates the transaction with an exclusive borrow
and queued work. submit takes ownership of an Edits batch and retains submission
order; it does not execute edits or check their preconditions, so it has no
Result error. The original molecule stays unchanged until commit starts.

commit consumes the transaction and is its single fallible execution boundary.
It moves the receiver's storage into an internal editor, executes the queued
batches through the direct mutation primitives, and records realized undo. An
empty valid value may occupy the destination during that exclusive borrow; it
is an internal transfer detail, not an error result. No recovery clone is made.
Each precondition is checked immediately before its edit against the evolving
state. Integrity is checked only after all batches, while recovery still owns
the editor and journal.

On success, commit transfers the accepted storage back and discards the journal.
On any edit or integrity error, it replays the journal and restores the original
before releasing the borrow. MoleculeApplyError distinguishes those failures.
The receiver contains either the accepted result or the original, never an empty
error sentinel. Do not call the public consuming editor.apply or finish in a way
that destroys the working state before failed-commit recovery can run.

Explicit rollback and Drop before commit only discard the queued work; no undo
has been recorded or is needed. There is no failed transaction handle left
callable after an execution error, because commit consumes it.

There is no public arbitrary-parts constructor, Clone, journal accessor, journal
composition, or unrestricted mutable editor access on Transaction.
Edits remains incrementally constructible as one sequence before submission.
Independent-batch composition is unresolved.
Multiple submitted batches share the transaction boundary; their New handles
retain batch-local meaning and each batch starts in the then-current id spaces
when executed during commit. Queuing multiple batches does not flatten their
handle namespaces implicitly.

```rust
let mut transaction = molecule.transaction();
transaction.submit(edits);
transaction.commit()?;
```

The common case is one edit list, although multiple submissions remain supported.
Deferred execution was chosen as a modest API simplification over immediate
fallible apply: it gives one error boundary and needs no aborted-but-still-held
transaction state. Immediate execution could report precondition errors earlier
or support dependent interactive work, but no intermediate-state inspection is
required here. Deferred submission retains Edits until commit and delays those
diagnostics. Commit can consume queued edits as the journal grows; no extra copy
of the batch is required.

There is no demonstrated substantial performance difference between immediate
and deferred execution for a single batch followed by commit. The measurements
below compare journaling and ownership choices, not execution timing. Transaction
submission and journal-free editor.apply therefore deliberately have different
names: submit queues work; apply executes and publishes it.

### Edits accumulation and independent batches — 2026-09-21

Retain Edits as one incrementally accumulated sequence with one initial-host id
frame and one creation namespace per entity kind. Its creation methods issue
handles for subsequent edits in that sequence. Resolver plans already accumulate
edits in this way; another batch-builder type is not needed. Do not replace this
with an immutable container whose caller must manage offsets in a Vec.

Distinguish three cases:

- Adding an edit written for the current sequence uses its existing handles;
  no rebasing is involved.
- Combining independently constructed batches must preserve each batch's handle
  meaning. For batches using the same initial-host ids, incoming New handles need
  per-kind rebasing; generic list extension does not supply that operation.
- Temporary fragments written for the final sequence already use its namespace.
  They are not independent batches and must not be automatically rebased. Use
  ordinary internal collections for such fragments rather than treating each as
  an independently composable Edits.

Independent-batch composition remains unresolved and unavailable in the proposed
surface for now. Do not add a generic extend/concatenation API, a caller-supplied
offset, or a Python-only implementation. Retained incoming handles, references
between batches, and initial-host versus successive-state ids are part of that
unresolved contract. Freely constructible raw handles still permit incorrectly
assembled sequences; container immutability alone cannot establish their intended
namespace.

If this work is resumed, move the sequence-assembly and handle-namespace semantics
currently implemented in Deltas lowering into methods owned by Edits, and
restructure the lowering calls to use them. Reaction lowering currently assigns
creation handles for the eventual sequence, collects several Edits fragments,
and pushes them into the final batch. See
[reaction lowering](../umol-graph-ir/src/ir/reaction.rs). This is the existing
producer to reconcile, not a precedent for independent-batch concatenation.
Edits must own any future rebasing, including references nested in edits and
ConstraintEdit; the Python binding must delegate to the Rust operation.

Keep the settled deferred transaction design. Multiple submitted batches retain
their boundaries and resolve their own handles when executed at commit. Merely
consuming handles during application does not require immediate execution.
Constructing a later batch by inspecting an earlier batch's actual result would
be a separate interactive requirement, not established by ordinary accumulation.

This records the design boundary, not an implementation change. The urgent
Python/Rust review proceeds without settling independent-batch composition.

### Tracking belongs to batch execution — settled 2026-09-20

Direct editor mutation is the fast primitive. It provides no public identity
tracking, atomic batch boundary, or rollback. finish checks aggregate integrity
before publication; it does not return provenance. Remove the editor's session
correspondence, public tracked direct-mutation methods, and tracked publication.
There is no tracked_finish or tracking mode for direct editing.

| Path | Checks and publication | Identity tracking | Recovery on failure |
| --- | --- | --- | --- |
| Direct mutation followed by finish | Direct operations may leave a transient draft; finish checks aggregate integrity. | None. | No rollback; failed finish drops the draft. |
| apply / tracked_apply | Per-edit preconditions and one final integrity check; no intermediate molecule is published. | tracked_apply only. | No journal; failure drops the working state. |
| commit / tracked_commit | Per-edit preconditions and one final integrity check across all submitted batches; no intermediate molecule is published. | tracked_commit only. | Journal-based rollback restores the borrowed receiver. |

Tracked batch signatures use the same draft error names as their ordinary forms:

```rust
impl Molecule {
    pub fn tracked_apply(
        self,
        edits: Edits,
    ) -> Result<(Molecule, MoleculeCorrespondence), MoleculeApplyError>;
}

impl MoleculeEditor {
    pub fn tracked_apply(
        self,
        edits: Edits,
    ) -> Result<(Molecule, MoleculeCorrespondence), MoleculeApplyError>;
}

impl<'a> Transaction<'a> {
    pub fn tracked_commit(self) -> Result<MoleculeCorrespondence, MoleculeApplyError>;
}
```

Molecule::tracked_apply maps its input molecule to the published result.
MoleculeEditor::tracked_apply maps the editor state immediately before the batch
to the published result; preceding direct mutations are outside that witness.
Transaction::tracked_commit maps the original receiver to the committed receiver
across all submitted batches. Deferred execution allows tracking to start at
commit without a transaction configuration flag. No tracked_submit or
tracked_rollback is needed: neither operation changes the molecule.

Correspondence records surviving entity identities and their resulting ids,
not attribute or participant changes. Participant replacement is supported by
tracked batches but produces no participant-mapping witness. General replacement
can change membership and length; it is not a frame permutation. The entity
correspondence does not assert frame equivalence. Edit/Delta old/new values and
saved undo components already describe the change. No consumer requiring a
separate positional mapping has been identified.

Correspondence is distinct from saved state for undo;
tracked_apply accumulates correspondence without constructing an undo journal.
Discarding a successful tracked result's correspondence preserves the ordinary
operation's output and final state. Tracked and ordinary forms have the same
checks, errors, and recovery behavior; failure returns no correspondence.

Untracked execution does not accumulate correspondence. Direct mutation still
computes compactions needed to maintain constraints and other internal references;
the graph-core primitives supplying them remain. This necessary bookkeeping does
not expose a tracked direct-mutation convenience. Batch execution uses the same
direct mutation implementation and accumulates provenance only when requested.

Review and migrate existing Rust and Python callers of tracked direct mutations
and tracked publication. Callers needing correspondence must express their changes
as Edits and use tracked_apply or tracked_commit. This is a settled API change,
not implemented behavior.

### Consuming transformation APIs

Resolve and project are common operations, particularly at SMILES boundaries.
Optimize their ordinary path for owned, journal-free execution. Eliminating
every clone is not a requirement; retaining a source intentionally can justify
one. Repeated recovery copies and journals inside the pipeline are avoidable.
Current transformation ownership signatures and unchanged-on-failure guarantees
are revisable, not compatibility constraints.

Proposed basic Transformer and graph-IR hydrogen signatures:

```rust
trait Transformer {
    type Error;
    fn transform(&self, molecule: Molecule) -> Result<Molecule, Self::Error>;
    // Existing enumeration remains a separate capability.
}

fn fold_explicit_hydrogens(molecule: Molecule)
    -> Result<Molecule, FoldHydrogensError>;
fn unfold_implicit_hydrogens(molecule: Molecule)
    -> Result<Molecule, UnfoldHydrogensError>;
```

The consuming transform is fundamental; the trait does not require a borrowing
transform_into that preserves its receiver on failure. Enumerating independent
alternatives can still require copies. The hydrogen chemistry, eligibility,
stereo substitution, error conditions, and idempotence in doc 166 remain; its
borrowed-input/candidate-replacement signatures are superseded by this draft
proposal if approved. No hydrogen implementation exists to benchmark yet.

Resolver::resolve and Resolver::project likewise take Molecule by value and
return the accepted molecule with their existing diagnostic information. They
can borrow it for planning, then transfer ownership through successive
applications. Later resolution phases depend on earlier results; do not force
the whole algorithm into one immutable-input planning pass. Publishing one
completed phase and moving its result into the next editor needs neither a
clone nor snapshotting, and retains the existing Molecule-based read APIs.

For a concrete result-shape draft that preserves the existing distinction between
chemistry outcomes and execution errors:

```rust
enum ResolveOutcome {
    Determined { molecule: Molecule, report: ResolveReport },
    Underdetermined { report: ResolveReport },
    Contradictory(ResolveContradiction),
}

enum ProjectOutcome {
    Determined(Molecule),
    Underdetermined,
    Contradictory(ProjectContradiction),
}

impl Resolver<'_> {
    fn resolve(&self, molecule: Molecule) -> Result<ResolveOutcome, ResolveError>;
    fn project(&self, molecule: Molecule, flags: ProjectFlags)
        -> Result<ProjectOutcome, ProjectError>;
}
```

These result names and shapes are proposed, not settled. The existing Solution
type uses the same payload type for Determined and Underdetermined; simply
putting a Molecule into that payload would also promise a molecule on
underdetermination. The draft makes no such change: partial-refinement publication
remains the separate resolver discussion in doc 217. Standalone phase methods
need the corresponding owned-result treatment; their planners can retain Edits
and existing Solution results.

### Consumer fit and remaining decisions

| Consumer | Ownership path |
| --- | --- |
| MOL parsing and SMILES ingestion | Move the freshly raised molecule into resolve; no recovery copy is needed. |
| SMILES export | Keep the existing intentional copy of the caller's source, then consume it through projection stages. Remove additional internal candidates and journals. |
| Python resolve | Keep its non-mutating, result-copy behavior: clone once at the wrapper and consume that copy in Rust. |
| Standalone Rust transformations | Transfer ownership through the consuming operation; callers choose explicitly whether to retain a source copy. |

Python editor and transaction wrappers still need exact lifecycle signatures.
A consuming Python wrapper must record consumption rather than leave an empty
Molecule pretending to be the original value. Tracking is confined to batch
application and transaction commit; direct editor correspondence accumulation
is removed by the settled design above.

After the Python/Rust review in 228, resume with Undo visibility, draft names
and errors, transformation outcome types, and Python lifecycle. Independent-batch
composition remains unresolved; if revisited, Edits must own the assembly
semantics currently held in Deltas lowering. The authorized experiment below
measures successful resolve/project execution with current storage and checks;
it does not implement the complete transaction redesign or authorize production
API changes.

## Feasibility measurements — 2026-09-20

### Decision supported

The consuming, journal-free transformation path is worth adopting for routine
projection. Across these fixtures, projection median time fell 28–37% with
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

Do not retain borrowing transformation APIs merely to preserve their current
recovery machinery. The measured consumers can use owned results, while export
and Python deliberately preserve a source at their outer boundary. No change to
chemistry selection, projection meaning, or final acceptance criteria was needed
for the successful workloads. Transactional mutation remains a separate recovery
contract; these timings do not suggest replacing it with destructive application.

### Method and scope

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
Thus the shared rows measure the operation after export/Python's intentional
source copy; they do not include the cost of making that outer copy.

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
| stereo | N[C@H](F)/C=C/C |
| combined | [13CH3][C@H](F)/C=C/c1ccccc1 |
| combined8 | Eight copies of combined, separated by dots |

combined8 is a synthetic disconnected scaling case, not a claim about production
molecule sizes. The fixture selection covers ordinary chains, aromaticity,
stereo, and their combination; it is not a corpus-weighted throughput estimate.

### Correctness checks

Before timing, the consuming resolver matched the current resolver's complete
Molecule output and report for all six fixtures. Both projection variants matched
the current complete Molecule output. Additional checks matched wildcard
underdetermination, a late discharge contradiction for C#c0#h4#n0#u0#s#D5, and a
late isotope-projection error after earlier stereo/aromatic phases. The baseline
projection preserved its source on that error; the consuming operation returned
the same diagnostic and discarded its owned state, as intended. These checks
establish equivalence for the measured cases, not general algorithm verification.

### Resolution measurements

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

### Projection measurements

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

### Remaining costs

Molecule clones increment Arc counts for graph and entity storage, but copy the
owned global Constraints vector. Subsequent mutations can copy shared tables;
consuming ownership avoids that when no other owner remains. It does not remove
algorithmic work, storage-index rebuilding, or copies required by an intentionally
retained external source.

Editor construction currently creates identity correspondence vectors for all
entity kinds, even for untracked use. That cost is retained in every measured
variant. The settled design removes that accumulation from direct editing and
untracked execution; these results do not quantify the benefit of that change.
Nor do they measure rollback or predict performance of the redesigned active
transaction.

The experiment's input definitions, method, checks, results, and interpretation
are retained here so scratch sources and logs can be deleted without losing the
design evidence.

## Design questions to settle together

1. **Owning operations.** Which mutations belong on typed entity sets, which need
   molecule-wide coordination, and which published Molecule conveniences should
   remain? Delegate relation changes to storage while keeping constraint updates,
   cascades, and identity tracking with the aggregate that owns those references.
2. **Participants and attributes.** Editor participant mutation and local getters
   are settled above, including independent attribute mutation and temporary draft
   inconsistency. Direct Molecule conveniences still need to define coordinated
   participant/attribute changes within their integrity boundary. Frame-preserving
   permutation is outside the current editor work.
3. **Execution and failure.** Specify signatures and ownership for the settled
   direct, journal-free, and active-transaction paths. Transactional precondition
   and integrity failures share whole-transaction rollback; direct draft
   publication remains separate. Specify the common underlying operations
   without exposing transaction history manipulation.
4. **References and witnesses.** Separate current dense ids, batch handles,
   batch/transaction correspondences, internal compactions, and rollback restoration.
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
   Carry the settled active-transaction lifecycle through its public surface and
   bindings; use the settled tracked batch/commit surface without retaining the
   detached-journal API.

Storage-owned mutation, typed graph-IR coordination, and the active transaction
lifecycle are settled above. Remaining execution/publication names, signatures,
and ownership must be settled before implementation planning.

## Evidence for the eventual implementation

Preserve the existing laws for handle namespaces, apply/transact result agreement,
rollback, failure recovery, constraint compaction, publication, and correspondence
composition. Relevant suites include
[edit properties](../umol-graph-ir/tests/property/edit.rs),
[publication](../umol-graph-ir/tests/property/molecule/publication.rs),
[compaction](../umol-graph-ir/tests/property/molecule/compaction.rs), and the
[reaction properties](../umol-graph-ir/tests/property/reaction).

For rollback, preserve exact restoration laws for matching history. Update misuse
coverage to require freedom from panics, without asserting a particular mismatch
error, restored result, or unchanged editor for manipulated inputs.

Extend coverage from the settled contract: all six entity kinds, participant and
payload coordination, legal and illegal frame alignment, cascade restoration,
reaction reference updates, and Python failure ownership. Include position-sensitive
payloads so that frame tests prove transport. Benchmark representative direct,
apply, and transact paths from the beginning, including copy-on-write sharing and
repeated relation mutation; do not assume the current fixtures establish scale.
