# 238 — Non-cascading removal and DPO application

Status: In Progress
Date: 2026-10-16
Relates: [184](184-deltas-and-edits-2026-08-04.md),
[213](213-editor-overlay-storage-2026-08-27.md),
[181](181-python-boundary-ownership-2026-08-03.md),
[239](239-kindless-stereo-removal-2026-10-17.md),
[data-type guide](../docs/development/data-types.md),
[integrity guide](../docs/development/integrity.md),
[nomenclature guide](../docs/development/nomenclature.md),
[Python API guide](../docs/development/python-api.md)

## Scope

Provide non-cascading topology removal in MoleculeEditor and Edit execution,
retain the existing cascading behavior under an explicit name, and route all
downstream DPO application paths through the non-cascading implementation.
This separation also applies to relation and overlay compaction: transporting
participant references must not implicitly select additional entity removals.
Complete the host-dependent dangling checks before mutation, including the bond
dependencies exposed by [GitHub issue #31](https://github.com/dayhofflabs/umol/issues/31)
and the corresponding stereo cases.

The application-tracking review also suggested two independent stereo changes:
removing kindless stereo entities, and requiring reactions to preserve the kind
of every preserved stereo entity. The combined experiment is complete; its
findings and potential kindless-removal scope are recorded in
[239](239-kindless-stereo-removal-2026-10-17.md). The decision on 2026-10-17 is to
retain the uniform entity interfaces and use the smallest tracking correction
below, together with
directly related consistency fixes under the existing contracts. Shared Updates,
FieldChange types, and EntityPatch implementations remain; kindless forms and
Default remain. The proposed separate stereo delta field-change types and removal
of the stereo EntityPatch implementations are rejected. Uniform treatment of
entities is a design principle, not incidental machinery to discard to simplify
one restriction. The combined audit below records the earlier experimental scope,
not a commitment to complete that migration.

Within the removal work, applicability checks and mutation have distinct
purposes. Pre-mutation checks establish DPO applicability. The removal primitive
determines what the mutation does. A complete pre-mutation dangling check makes
a cascading primitive usable for DPO: there is then no additional dependent
entity for it to delete. The non-cascading primitive makes the intended deletion
explicit and avoids rediscovering dependent removals after the caller has
already established that none remain.

Implementation began on 2026-10-17 with S1a. The subitems below record completion;
graph-IR removal retains its existing cascading behavior until S2e.

### Exclusions

- Molecule-level constraint compaction and removal retain their existing
  behavior. Constraint handling is not the cascading behavior being changed
  here; extending this work to prohibit automatic constraint removal was not
  approved.
- No change to aromatic-system storage or serialization: systems continue to
  store atom participants, without a new stored bond list.
- No general entity-dependency DAG or participant-algebra redesign.
- No SqPO reaction application, automatic overlay repair, or inferred replacement
  of dependent overlays after topology deletion. Lowering an explicitly requested
  stereo kind change to removal and addition is not such a cascade. Cascading
  editor conveniences do not authorize cascading at the reaction layer.
- No additional restrictions on bond-order changes or bond addition, and no
  chemistry validation in removal primitives.
- No manuscript changes.

## Implementation before this work

The implementation was inspected on 2026-10-16.

| Surface | Current behavior |
| --- | --- |
| MoleculeEditor::remove_topology | Removes requested atoms and bonds, incident bonds, and overlays whose stored participants or sites are removed; compacts constraints. There is no non-cascading counterpart. |
| Edit::RemoveTopology | Performs the same removal through plain or journaled execution. |
| Edits::remove_atom / remove_bond | Construct RemoveTopology; they do not provide another removal behavior. |
| Explicit editor and Edit overlay removals | Remove only the requested overlay entries among the entity families. They also compact molecule constraints. Removing a stereo bond does not remove its site bond or ligand atoms. |
| Internal Molecule topology removal | Calls Graph::tracked_remove_cascading and compacts atom/bond attributes; overlay and constraint handling occurs separately. |
| Graph::try_remove / try_tracked_remove | Check the dangling condition first, then execute tracked_remove_cascading. The check prevents extra edge deletion, but the execution still discovers incident edges again. |

The relevant implementations are in
[editor.rs](../umol-graph-ir/src/ir/molecule/editor.rs),
[apply.rs](../umol-graph-ir/src/ir/molecule/apply.rs),
[edit.rs](../umol-graph-ir/src/ir/edit.rs), and
[graph.rs](../umol-graph-core/src/graph.rs).

During the current cascading topology removal, the GraphCompaction is passed
to each overlay store. Its compact operation transports participant references
and drops rows when any required participant has no surviving image, producing
a separate compaction for that overlay's IDs. For StereoBonds, deleting the
localized site bond therefore drops the stereo-bond row before final molecule
integrity checking can observe its missing site. The graph and overlay
compactions are then assembled into a MoleculeCompaction.

The compaction components remain disjoint by entity kind. Explicit stereo-bond
removal produces a stereo-bond compaction and leaves the graph compaction as
identity. The implicit cascade is in the relation-store compact operation that
derives additional removals from missing participant images, not in the
Compaction or MoleculeCompaction constructors themselves.

This split is not a missed migration from GraphCompaction to MoleculeCompaction.
Commit a2aa3e7cb (2026-08-26) replaced IdCompaction, which already combined a
GraphCompaction with six overlay-removal lists. The same commit updated the
editor's topology-removal path while retaining graph compaction as the input
to cascading overlay removal. Doc 213 subsequently specified that composition
explicitly. This establishes the history, not the desirability of the implicit
cascading contract.

By comparison, an actual stereo ligand's adjacency bond is not a stored BondId
participant. Deleting that bond can leave the stereo row present, and final
integrity checking rejects its ligand incidence. Aromatic systems store atom
participants, so deleting a bond between their members does not remove the
system or trigger an equivalent representation-integrity failure.

Reaction::apply_deltas_at currently checks incidence around deleted host atoms
before lowering removals. It lacks the corresponding checks for deleted host
bonds. That produces different outcomes depending on which dependency is
explicitly stored and which is implicit in incidence.

## Settled decisions

### Separate applicability from mutation

1. Reaction application is strictly DPO. Host-dependent dangling checks run
   against the host, match, and explicit deletion plan before any mutation can
   remove the evidence needed for those checks.
2. MoleculeEditor and Edit are mutation interfaces, not implementations of
   reaction DPO applicability. Where cascading previously prevented missing
   references, exact removal checks the required reference survival before
   mutation and rejects invalid input. The core mutation primitives rely on
   those established prerequisites without repeating the checks.
3. MoleculeEditor::probe and finish retain their representation-integrity gate.
   That gate is not a substitute for the pre-mutation DPO check: a cascade can
   leave an integral product after deleting an entity the rule did not remove.
4. All downstream DPO implementations use the non-cascading removal path after
   their applicability checks. Cascading convenience operations remain available
   for callers that explicitly request their behavior.
5. Constraint compaction remains unchanged in both topology-removal paths.

### Graph-core removal

Add Graph::remove and Graph::tracked_remove as the direct counterparts to the
checked removal operations. These remove exactly the specified nodes and
edges, preserve survivor order, translate surviving endpoints, and rebuild
adjacency. They neither check the dangling condition nor discover additional
incident-edge removals. Every retained edge must have surviving endpoints.
Violating this prerequisite is a programming error and may panic during
translation; normal callers must establish it before invoking the primitive.
The bulk operations consume their removal vectors, as specified in the
interface inventory below.

Existing try_remove and try_tracked_remove check the dangling condition and
then use this backend. Existing cascading operations determine the additional
edge removals and then use the same backend. The try_ prefix follows the
nomenclature guide's checked/asserted counterpart convention, rather than
merely announcing an Option result; try_tracked_remove retains that prefix
order.

### Public topology-removal interface

| Surface | Non-cascading operation | Existing behavior retained as |
| --- | --- | --- |
| MoleculeEditor | remove_topology(atoms, bonds) | remove_topology_cascading(atoms, bonds) |
| Edit | RemoveTopology { atoms, bonds } | RemoveTopologyCascading { atoms, bonds } |
| Edits | remove_topology(atoms, bonds) | remove_topology_cascading(atoms, bonds) |

Bare removal deletes the explicitly supplied atoms and bonds without implicitly
deleting additional bonds or overlays. Cascading removal retains the old
topology and overlay behavior. Both retain the current molecule-constraint
handling.

The combined topology operation is retained so an atom/bond deletion batch can
compact both ID spaces together. Separate RemoveAtoms and RemoveBonds variants
are not needed for this design. Edits::remove_atom and remove_bond continue to
construct the bare RemoveTopology variant and therefore become non-cascading.

Editor methods retain their mutating, journal-free ownership model. Edit batches
retain their existing handle namespaces, old-state checks where applicable,
plain/journaled execution, and recovery contracts. The new distinction does not
make attribute mutation fallible or introduce per-removal DPO validation.

The non-cascading topology path obtains GraphCompaction from graph-core,
compacts atom/bond attributes, and translates the surviving overlay participants.
Graph IR then assembles MoleculeCompaction with that graph mapping and six
count-bearing identity overlay mappings. The identities describe overlay row
IDs, even when the stored participant IDs have shifted. Constraints and edit
handles use the resulting MoleculeCompaction.

An explicit overlay removal instead produces a MoleculeCompaction with its own
row mapping and identities for the graph and other overlay families. The
cascading topology path assembles the graph mapping with the actual mappings
returned by exact removal of the selected dependent overlays. These
compositions use existing Compaction, GraphCompaction, and MoleculeCompaction
constructors; no new mapping type or general-purpose constructor is needed.

### Explicit cascading in relation compaction

The distinction must extend through graph-core relation-set compact and
tracked_compact operations, graph-IR overlay stores, and their molecule
delegates. Adding a non-cascading topology entry point while leaving implicit
row deletion under bare compaction names does not satisfy this design.

Overlay row removal and participant translation are separate responsibilities.
Row removal applies the selected removals in that entity's ID space;
compact_participants translates graph references in surviving rows without
deleting rows. Neither operation discovers additional dependent removals.
Where participant-triggered deletion is retained, expose it as an explicitly
named cascading operation and route cascading callers through it. A returned
relation compaction records the removals performed by the selected operation;
returning that mapping does not itself authorize additional removals.

The removal of implicit cascading, the core programming-error boundary, and
pre-mutation rejection of invalid Edit input are settled. The agreed signatures
and shared error are specified below. Molecule-constraint compaction remains
unchanged.

Applying an already specified MoleculeCompaction and discovering a cascade are
different operations. The former takes overlay row removals from the relevant
component and translates participants through the graph component. The latter
uses graph removal and overlay storage to discover additional row removals,
thereby producing components of MoleculeCompaction. The current overlay compact
methods perform the latter operation. GraphCompaction alone is sufficient for
participant-only translation, but cannot specify independent overlay removals.

A completed MoleculeCompaction already contains the overlay removals selected
by cascading. Applying it does not need the original topology to rediscover
those removals. Consulting the molecule is necessary when deriving that
MoleculeCompaction from graph removal, not when interpreting its completed
overlay components. Constructing the mapping and applying it must not be
conflated.

### Bond dependencies at reaction application

The narrow fix uses the existing encoding. Before deleting a matched host bond,
application must account for every retained overlay that requires that bond:

- **Aromatic systems:** a localized bond whose two endpoints belong to the same
  aromatic system is support for that system. Deleting it while retaining the
  system makes the match inapplicable. This conservative rule does not require
  another ring-perception pass or an attempt to repair the system.
- **Stereo atoms:** the bond connecting the site to an actual-atom ligand is
  required incidence. Virtual hydrogens and lone pairs do not introduce such a
  bond dependency.
- **Stereo bonds:** the site bond and the bonds connecting each endpoint to its
  actual-atom ligands are required incidence. Endpoint-block and virtual-ligand
  semantics remain unchanged.

Explicitly removing the dependent overlay in the rule can satisfy the dangling
condition. Any replacement must also be explicit in the rule. Removing and
adding a bond with the same endpoints does not turn deletion of the original
bond into preservation of that bond.

A bond-order change preserves the bond; adding a bond deletes no entity.
Neither violates the dangling condition, and both remain permitted without
requiring an aromatic-system removal. Whether the resulting bond orders and
topology conform to the asserted aromaticity is a Tier 3, model-dependent
validation question. The deletion restriction does not imply a symmetric
restriction on modification or addition. Reaction application does not perform
that chemistry validation implicitly.

Match-local dangling failure retains the existing Ok(None) application result;
application iterators skip that match. Plain, tracked, realized-reaction, and
reaction-span results must have the same applicability. This work does not turn
ordinary non-applicability into an iteration error.

## Caller audit and migration scope

| Caller or surface | Required treatment |
| --- | --- |
| Graph::try_remove / try_tracked_remove | Retain their checked contract; after the check, execute non-cascading removal through the agreed graph-core interface. |
| Graph::pushout_complement / tracked_pushout_complement | They delegate to checked graph removal and must reach its non-cascading backend. Preserve their current absence and correspondence contracts. |
| Reaction::apply_deltas_at | Complete atom/bond dependency checks before mutation and lower to bare RemoveTopology. Do not rely on a later cascade or integrity gate for DPO applicability. |
| All Reaction application result forms and iterators | They share apply_deltas_at; verify consistent applicability and product identity across the entry points. |
| Reaction::compose | compose_overlap invokes tracked_apply_at. Cover both component applications and the resulting composite, including explicit overlay removal and replacement. |
| Python reaction and molecule reaction methods | Continue delegating to the Rust operations; verify the same DPO behavior without a separate Python implementation. |
| MoleculeEditor, plain Edit execution, and journaled Edit execution | Check exact-removal prerequisites before any mutation by that operation and reject invalid inputs. Implement both removal behaviors consistently; preserve compaction, handle transport, and recovery laws. |
| Molecule::extract / tracked_extract / edits | These currently depend on cascading extraction behavior. Route that behavior explicitly through the cascading operation; do not silently change extraction into DPO application. |
| Undo execution | Audit topology removal used to undo additions separately from forward removal. Preserve the restoration contract, including its behavior for manipulated history. |
| Other topology-removal callers, examples, and generators | Classify whether they request a complete deletion or rely on a cascade. Migrate based on that contract, rather than mechanically renaming every call. |
| Explicit overlay removals used by resolvers, projectors, and Kekulizer | Their behavior is outside this topology-removal split; molecule-constraint cleanup is not redesigned here. |

ReactionSpan side projection constructs the selected side directly rather than
using topology-removal edits. Its existing faithful-projection and integrity
contracts remain applicable; it must not acquire a cascading repair path.

### Application correspondence and replay

The current application result builds an atom correspondence and induces bond
and overlay pairings from incidence. Consequently, explicitly removing and
recreating an aromatic system over the same atom set can be reported as a
preserved system. A realized reaction can then omit the system removal even
though the original rule explicitly included it.

This matters when the new dangling check requires that explicit removal:
replaying the realized reaction must not fail because application erased the
deletion identity from its own witness. Audit the application-issued
correspondence, realized reactions/spans, and composition together. Composition
currently reconstructs its reaction from sides using an atom correspondence,
so delegation to the corrected application path alone does not prove that it
preserves explicit replacement identity.

Faithful deletion/addition identity and replay are settled requirements. Reusing
Molecule::tracked_apply for the complete host-to-product correspondence is
accepted. For composition, the proposed trial reverses and composes complete
MoleculeCorrespondence values, then uses Molecule::difference_to and Reaction::new
to construct the composite. Acceptance of that composition implementation requires
measurements; it is not settled by acceptance of the application change. Compare
it with atom-only composition followed by incidence induction while holding the
application changes fixed, and report the complete composition cost. The combined
stereo experiment does not isolate this change. Any selected implementation must
preserve explicit replacement identity. S0b records the verified producer
prerequisites and failure handling for this trial. General edit identity
and the correspondence carrier's contract remain unchanged.

Stereo kind is part of entity identity in reactions. Changing between different
determined kinds represents removal of the original stereo entity and addition
of its replacement. This rule applies to application, its correspondence and
realized results, direct span construction, reversal, and composition. It is a
settled requirement, not a choice left to implementation. Shared Updates,
FieldChange types, and EntityPatch implementations remain; generic edit tracking
and kindless refinement retain their existing contracts.

The current supplied-match application executes a kind change as a field edit,
then leaves the pair unmatched when inducing its correspondence. Direct
Reaction::to_reaction_span instead emits Modified and rejects the different
determined kinds. Correct both producers: reaction lowering must implement the
same replacement semantics that span construction expresses as Removed and
Added. Keep the LHS entity prefix, assign a fresh ID to the replacement, and
transport RHS references and ligand frames accordingly. ReactionSpan's existing
rejection of a different-kind Modified entry remains correct. Reversal continues
through the resulting span; it needs no new interface.

Replacement follows the existing removal/addition rules for constraints. The
replacement receives the constraints specified by the resulting rule side;
additional constraints found only on the matched host entity are not implicitly
inherited. Molecule-level constraints referring to the removed entity follow the
existing constraint-removal behavior. Constraints explicitly supplied for the
replacement refer to its new ID. Preserving the old field edit's retention of
host-only constraints is not a requirement. Application and span conversion must
agree on replacement behavior, rather than merely report different identity
while continuing the old field mutation.

Use the edit-produced correspondence without re-inducing other entity pairings
or pairing an explicitly removed entity with its replacement. At the reaction
boundary, exclude stereo pairs whose two determined kinds differ. Apply that
rule to the final sides of composition as well: Th -> kindless -> Sp can preserve
a pair at each step, but its final Th -> Sp pair must be unmatched. Perform this
adjustment before difference_to; otherwise a valid sequential transformation can
be lost when span construction rejects the pair. General correspondence
composition remains unchanged. Include this work in the composition measurement.

This route needs focused checks of realized reactions, spans, and composition
as part of the existing tracking trial. It is not an implemented or measured
solution. Related consistency fixes must restore the existing contracts across
their producers and consumers; they do not authorize the stricter no-kind-change
delta semantics or a separate stereo patch interface. Removing the enumeration-only
kind-change prohibition is accepted: prepare_deltas must not require the new
configuration's kind to equal the old or LHS kind. Retain the existing old/LHS
compatibility and coset-range/expression-action checks. Supplied-match application
already permits the kind change; both routes report different determined kinds
as unmatched in the result correspondence. Sharing the other validity checks
between the two entry points remains outside this decision. The broader stereo
redesign is not part of this implementation.

The stereo delta lowering correction is accepted: StereoAtomDelta::for_update and
StereoBondDelta::for_update derive constraint context from the resulting
configuration. An explicit kind clear must not fall back to the source kind or
emit vacuous constraint changes. Retain the configuration change and the existing
infallible signatures, matching the established Edit-lowering contract.

The expected costs, based on the current implementation, are:

- Application already constructs a complete correspondence, including for
  result forms that discard it. The trial keeps the same host clone and final
  publication check; tracked_apply does not introduce a journal or another
  molecule clone. Kind-changing reactions lower to replacement operations, so
  their edit execution is not assumed to have the old field-mutation cost.
- The current reconstruction scans source atoms to construct checked survivor pairs,
  and induces seven other components using incidence keys, BTreeMaps, and
  correspondence validation. Variable participant keys also require allocation
  and sorting. The trial removes that host-to-product reconstruction; validating
  the supplied rule-to-host match remains necessary.
- tracked_apply instead allocates and initializes identity pairs for every
  source entity, up to eight nonempty pair buffers. Removal tracking scans
  existing handle state, constructs removal compactions, and filters/translates
  the pairs in place. It retains the pair-buffer allocations but still incurs
  temporary removal vectors and atom/bond ID conversions. Its identity
  correspondences are dense, unlike the count-only identity compactions. Large
  deletions therefore initialize pairs that will be discarded. Skipping
  unchanged handle families does not by itself eliminate correspondence scans.
- Composition reverses and composes all eight components instead of only the
  atom component. This adds pair allocation, validation/sorting during reversal,
  and lookup work for the other components, while eliminating their subsequent
  incidence induction. Span construction, difference calculation, and required
  validation remain.

For application, avoiding incidence-key reconstruction is a plausible saving;
the dense initialization and composition work prevent a blanket speedup claim.
Use a bounded comparison of complete application and composition to accept this
implementation route. A material regression calls for identifying repeated
work or revising the route, not dropping faithful replacement provenance.

## Stereo redesign — bounded experiment

The combined experiment described below concluded on 2026-10-16. Its results,
failures, and remaining migration scope are recorded in
[239](239-kindless-stereo-removal-2026-10-17.md). This section retains the original
hypotheses and scope audit for that combined trial.

Exploratory design and current-code audit, 2026-10-16. The redesign is a plausible
direction, not an approved migration. Its possible simplifications must be
weighed against the implementation and the semantic/interface problems exposed
by a prototype. This investigation remains separate from the removal and
tracking work above.

The experiment's purpose is to decide whether to adopt or discard the design.
Before code work starts, define its limit and stopping point and preserve a
working baseline. Keep the prototype isolated and disposable. Report which
existing checks pass unchanged, which expectations intentionally change, all
failures, and remaining unknowns. Passing tests alone do not establish that the
design is preferable. New architectural dependencies do not automatically expand
the experiment. Its conclusion may be that the simplification is insufficient
or the migration too costly. No redesign implementation subitems follow until
there is an explicit adoption decision.

### Working hypothesis and constraints

The mandatory-kind and reaction-restriction changes below are tentative. The
constraints on infallible direct mutation, ModifyField naming, and faithful
serialization remain requirements for any experiment.

1. Trial requiring every stereo entity to have a definite StereoKind. Remove the kindless
   StereoConfigurationForm::Undetermined state. A known kind with an undetermined
   StereoCoset remains meaningful: Th* and Ct* remain representable.
2. Trial requiring reactions to preserve the kind of every preserved stereo entity. Changing
   kind requires explicit removal and addition. This rule includes application,
   spans, realized reactions, inversion, and composition. It is independent of
   hypothesis 1; taken together they would remove unknown-kind refinement and identity
   exceptions throughout the pipeline.
3. General edits and direct attribute mutation remain capable of changing one
   known kind to another. Direct attribute mutation remains infallible. This
   does not strengthen molecule integrity to require kind/frame agreement,
   valid coset indices, or consistent entity-level constraints. Existing raw-form
   and operation-specific failure contracts remain in force.
4. Keep ModifyField. Do not replace it with a special ModifyCoset delta variant
   merely because stereo currently has one inherent field. StereoAtomUpdate,
   StereoBondUpdate, and their shared supporting types may be redesigned; their
   current shape is not a constraint on the solution.
5. Remove Default from stereo forms that cannot supply a meaningful kind.
   Empty owning collections retain Default. A no-op update is distinct from a
   kindless entity; its exact representation belongs to the update design.
6. A temporary replacement with the original kind, undetermined coset, and empty
   constraints is acceptable for moving attributes out of a mutable reference.
   StereoAtomForm::new(kind, StereoCoset::Undetermined) and the corresponding bond
   constructor already express it. A new convenience constructor is optional,
   not a prerequisite. No new graph-core storage transformation is required just
   to remove the form's Default implementation.
7. Remove bare kindless stereo syntax and kindless constraint-vacuity rules. Do
   not silently translate old kindless input to a guessed kind or discard an
   entity. Ordinary syntactic elision of an undetermined constraint remains a
   separate matter. Preserve the fixed-kind TetrahedralStereoForm and
   CisTransStereoForm variants, including Undetermined and NotStereo, and their
   atom #T / bond #C spellings.
8. The different kind domains have no common top. Different satisfiable kinds
   meet to bottom and join through the existing NoJoinError channel. Within a
   kind, retain coset and constraint expressivity. The raw contradictory-operand
   case and the is_undetermined predicate need the explicit decisions below.
9. Reject an added degree-zero Any/Null kind whose purpose is to provide Default.
   Also reject automatic changes to mutation fallibility, silent omission of an
   unsupported kind change, or catching panics as a recovery strategy.

The experiment would remove an entity-level pattern saying that a stereo entity
exists with arbitrary geometry, and the corresponding cross-kind generalization.
Omitting a stereo entity is not an equivalent assertion. The historical rationale
in [113](113-ast-canonical-equality-and-lattice-2026-06-14.md)
describes that wildcard. The current Lattice trait already supports types without
a common top; parsing also does not require a kindless entity placeholder.

The kindless treatment recorded and implemented in
[181](181-python-boundary-ownership-2026-08-03.md#kindless-stereo-constraints--decision-2026-10-14)
remains the current behavior unless this redesign is adopted. Preserve that
completed record rather than rewriting its historical decisions.

### Normalization, reframing, and ownership

Normalize and FrameTransport do not require Default. The dependency is in the
collection implementations' mem::take calls, including StereoAtomSpans and
StereoBondSpans through EntitySpan<T>::default.

The proposed replacement changes ownership plumbing only. Derive the action
from the stored ligand frame, move out the complete original attributes, and
normalize/transport those original attributes. The temporary wildcard does not
participate in the operation. Install the successful result before transporting
the stored participants or reporting an action. Both sides of a Modified span
must still receive the same action.

These paths consume their owning collection. Arc::make_mut detaches shared
storage before replacement. Failure drops the working collection, and surviving
clones retain their original attributes and frames. The replacement itself
requires no allocation or payload clone; it copies the existing kind and creates
an empty constraint container. Existing normalization/transport costs and
fallibility are separate. Do not introduce another normalization pass, rebuild
participant storage, or add runtime rollback for the temporary replacement.

### Public surface and consumer audit

This inventory covers direct symbol changes and consumers whose semantics must
be checked even when their signatures need not change. It does not authorize
additional API families or unrelated cleanup.

| Surface and current source | Consequence if adopted or design dependency to investigate |
| --- | --- |
| StereoConfigurationForm; StereoAtomForm and StereoBondForm; StereoConfigurationUpdate, StereoAtomUpdate and StereoBondUpdate in [ir/stereo.rs](../umol-graph-ir/src/ir/stereo.rs) | Remove kindless construction and clearing; settle the replacement shape and constructors. Update kind/coset access, update/difference_to, From conversions, AsLit, structural equality/hash/order, Normalize, Lattice, FrameTransport, and apply/swap/mirror. Keep StereoConfiguration as the concrete kind/literal-coset value. |
| StereoAtoms, StereoBonds, their entry implementations, StereoAtomSpans and StereoBondSpans in ir/stereo.rs; EntitySpan in [ir/delta.rs](../umol-graph-ir/src/ir/delta.rs) | Replace form-dependent mem::take uses without changing collection identity or frame transport. Generic EntitySpan need not acquire a new empty variant. Retain its existing Default only where T itself has Default. |
| StereoAtomConstraintsForm and StereoBondConstraintsForm in [ir/constraint/stereo.rs](../umol-graph-ir/src/ir/constraint/stereo.rs); [ir/view/constraints.rs](../umol-graph-ir/src/ir/view/constraints.rs) | Keep constraints governed by the owning entity's kind. Remove kindless asserted-reading and normalization bypasses. Do not add an intrinsic kind to every individual constraint. Preserve Option-for-bottom reading semantics and raw mutation. |
| Constraint::StereoAtom / StereoBond in [ir/constraint/molecule.rs](../umol-graph-ir/src/ir/constraint/molecule.rs), including logical constraints | These already carry a definite kind for interpretation outside an entity. Keep that context and the existing molecule-level constraint integrity checks; removing kindless entity forms does not remove logical stereo assertions. |
| Stereo views in [ir/view/stereo.rs](../umol-graph-ir/src/ir/view/stereo.rs), atom and bond projections in ir/view/atom.rs and ir/view/bond.rs | Form kind/coset access can become total. View kind/coset already return direct values but currently panic on kindless configuration; eliminate that case. Retain absence for a missing overlay and existing rejection of incompatible projection kinds. |
| StereoAtomFieldChange and StereoBondFieldChange; Edit/Undo configuration and constraint variants; Edits::update_stereo_atom / update_stereo_bond in [ir/edit.rs](../umol-graph-ir/src/ir/edit.rs); [ir/molecule/apply.rs](../umol-graph-ir/src/ir/molecule/apply.rs) | Preserve general configuration replacement, checked old-value application, batching, journal restoration, and ID tracking. Decide how the shared field-change/update types separate this capability from reaction modification. Remove kindless constraint-export cases. |
| StereoAtomDelta and StereoBondDelta, for_update, EntityPatch implementations, delta folding, inversion and frame transport in ir/delta.rs; EntityPatch in [ir/traits.rs](../umol-graph-ir/src/ir/traits.rs) | Make reaction modification preserve kind. Reconcile total attribute-diff laws and required kind context with the chosen update/change representation. The current trait's kind-free modify_constraint fallback cannot remain a producer of kindless stereo deltas. |
| Reaction construction/integrity, prepare_deltas, supplied-match application, reframe_deltas and delta-to-Edit lowering in [ir/reaction.rs](../umol-graph-ir/src/ir/reaction.rs) and [ir/reaction/integrity.rs](../umol-graph-ir/src/ir/reaction/integrity.rs) | Use one fixed-kind reaction contract across all entry points. Remove unknown-kind compatibility branches. Decide the earliest existing failure boundary for independently supplied incompatible changes; preserve checks for malformed coset expressions and actions. |
| ReactionSpan publication, superimpose, projection and to_reaction in [ir/reaction_span.rs](../umol-graph-ir/src/ir/reaction_span.rs); Molecule::difference_to; [ir/compose.rs](../umol-graph-ir/src/ir/compose.rs) | Replace the Some/Some-only kind check with unconditional agreement for a preserved entity. Keep added/removed sides and infallible conversion from a valid span. Compose complete provenance so explicit removal/addition does not become modification. |
| MoleculeCorrespondence::induce in [ir/correspondence.rs](../umol-graph-ir/src/ir/correspondence.rs); [ir/substructure.rs](../umol-graph-ir/src/ir/substructure.rs); [ir/molecule/pushout.rs](../umol-graph-ir/src/ir/molecule/pushout.rs) | Equal kinds can pair; different kinds remain unmatched when inducing side identity. Pattern matching and pushout retain their existing attribute-meet and participant-frame contracts. Do not impose reaction kind preservation on arbitrary correspondence construction or generic edit tracking. |
| [ir/canonicalize.rs](../umol-graph-ir/src/ir/canonicalize.rs), [ir/coloring.rs](../umol-graph-ir/src/ir/coloring.rs), and [ir/symmetry.rs](../umol-graph-ir/src/ir/symmetry.rs) | Remove kindless structural keys, descriptor construction, and missing-kind symmetry branches; coloring hashes a definite kind. Preserve frame-aware distinction within each kind, including unknown cosets. Retain existing handling of raw kind/frame or coset incompatibility. Canonical numeric hashes/representatives need not remain byte-stable across this representation change. |
| [dsl/stereo.rs](../umol-graph-ir/src/dsl/stereo.rs), [dsl/edit.rs](../umol-graph-ir/src/dsl/edit.rs), [dsl/reaction.rs](../umol-graph-ir/src/dsl/reaction.rs), dsl/molecule.rs, dsl/reaction_span.rs and dsl/constraint.rs | Migrate compact syntax, tree and streaming EDN paths, rendering, lowering, and all public Dsl conversions together. Standalone constraint Dsl wrappers already carry a definite kind. Update wrappers and self-contained edit/delta serialization need the kind-context decision below. No panic-based formatter or alternate serialization family. |
| [umol-py/src/stereo.rs](../umol-py/src/stereo.rs), src/constraint/stereo.rs, src/edit.rs, src/delta.rs, src/reaction.rs and src/lib.rs | Mirror settled Rust shapes and construction. Change configuration variants/getters, update constructors and read-only views, field changes, constraint-modification context and exhaustive from_rust/to_rust conversions. Preserve live attribute mutation, readonly forms retained in Edit/Delta, ownership and iterator failure contracts. Structural matching remains optional for Python 3.10+; ordinary construction/access must work on 3.9. |
| [umol-graph/src/ops/stereo.rs](../umol-graph/src/ops/stereo.rs), ops/resolve/stereo.rs, ops/validate/stereo.rs; [umol-io/src/table_ir/raise.rs](../umol-io/src/table_ir/raise.rs), table_ir/stereo.rs, layout/stereo.rs and depiction consumers | Perception/raising already produce named kinds. Migrate optional-kind tests and direct construction, retain unresolved cosets and the external-format boundary semantics. Do not add stereo perception, invent a geometry for unspecified input, or change export policy. Several additional direct references here occur only in tests. |
| Public exports in ir.rs and dsl.rs; development guides and [DSL specification](../umol-graph-ir/spec/umol-dsl-spec.md) | Reconcile every changed symbol and signature, describe subkind identity, remove obsolete kindless and Default claims, and align the formal grammar with both parsers. No manuscript edits. |

No graph-core or umol-perm representation change is required by the selected
wildcard replacement approach. Geometry/model support, the six existing
StereoKind values, and the low-level permutation/coset algebra are not being
expanded. The existing one-stereo-entity-per-site integrity rules remain in
force. Unrelated optional values, such as a missing overlay, a missing
constraint, an absent span side, or a bottom reading, retain their meanings.

### Design risks to investigate within the experiment

These SIN questions expose where the hypothesis is incomplete. They are risks
for the experiment, not prerequisites for the tracking fix or decisions that
must all be settled in advance. A coherent candidate should resolve routine
choices together and identify any remaining semantic trade-offs for review.
The experiment can be rejected without resolving every question.

#### SIN1 — Form shape, constructors, and the per-kind top

The narrow representation is a StereoConfigurationForm struct containing
StereoKind and StereoCoset, retained in the existing configuration field of the
two entity forms. It avoids a one-variant Kinded enum and retains the distinction
from concrete StereoConfiguration. This is a proposal, not an approved public
layout. Settle the public fields/accessors and whether kinded becomes new; existing
entity new(kind, coset) already supplies the wildcard replacement.

Also settle is_undetermined. Th* must not become a cross-kind wildcard, but the
predicate can mean top within its kind, as keyed constraint forms already mean
top within their key, or it can remain a global reading that is always false
for the kind-bearing form. The per-kind reading fits subkind identity; all
consumers must still preserve the kind and the existence of its entity. With
constraints, the entity's predicate must account for their information too.

#### SIN2 — Updates, field changes, and the patch algebra

Both general Edits and reaction deltas currently use StereoAtomFieldChange /
StereoBondFieldChange with complete old/new configurations. They also consume
the same StereoAtomUpdate / StereoBondUpdate. Form difference_to currently
produces a total update between normalizable forms even when their kinds differ.
EntityPatch::diff promises a modification-only reconstruction between its
Attributes values, and stereo for_update has no site or ligand arguments with
which to construct removal/addition.

Those contracts cannot all remain unchanged when reaction modification loses
kind-changing capability. Settle together whether Update remains a general edit
description or becomes a modification within a fixed kind; where its kind is
carried; and how the field-change/patch interfaces express the reaction domain.
General configuration edits must remain available either way. Do not silently
narrow difference_to's law, emit no changes on mismatch, assume a hidden
same-kind precondition on an unrestricted public helper, or automatically add
Result to all producers. Equally, do not preserve an ill-fitting trait by
retaining a kindless fallback. Changing the shared types or trait context is
available; removing public helpers is not a settled decision.

The distinction between omitting a coset update and explicitly assigning the
undetermined coset must survive. The corresponding empty/no-op update syntax,
constraint-only updates, kind-only spelling, and full configuration edits must
all follow from this one design, in Rust and Python.

#### SIN3 — Kind context in constraint changes and serialization

Edit's stereo constraint singletons/batches and the stereo delta constraint
variants currently carry Option<StereoKind>. It supplies DSL context; application
does not compare it with the host's kind, and None can cause an entire edit to
be omitted on export. Update formatting likewise omits constraints when its
configuration supplies no kind. Removing kindless entities does not alone fix
these independently constructible values.

The revised representation must provide definite context for every serialized
stereo constraint change and preserve its effect. Settle where that context lives
with SIN2, including one shared kind for a homogeneous constraint batch. It is
context for the operation, not an independently mutable intrinsic kind on each
constraint. Clarify whether edit-side context remains rendering context or also
becomes an application precondition; an additional host-kind check is not already
approved. Reaction-side context must agree with its preserved source kind.

The generic EntityPatch::modify_constraint signature currently has no context
argument. Its stereo implementations explicitly manufacture kind: None and rely
on an overridden diff to avoid that path. This is a concrete interface gap to
resolve, not a reason to restore kindless forms or silently omit serialization.

#### SIN4 — Where invalid reaction modifications are rejected

Reaction::try_new currently permits old/new configuration kinds to differ.
Enumeration's prepare_deltas rejects such changes through
ApplyPreconditionError::InconsistentReaction; supplied-match application does
not perform that same preparation. ReactionSpan publication already reports
StereoKindModified for different determined sides.

The proposed rule would establish uniform kind preservation. Once the experiment
establishes what raw deltas could express, determine
whether the remaining independently supplied kind mismatch is rejected by
Reaction construction/integrity or by a shared application precondition. This
must cover both field and constraint modifications and their agreement with the
lhs, including raw DSL/Python construction. Prefer the existing failure
boundaries and vocabulary; no new fallible direct attribute setter follows from
this decision. Valid span-to-reaction production must retain its infallible route.

#### SIN5 — Join on contradictory raw operands

Different satisfiable kinds have no join. There is a separate existing edge
case: StereoConfigurationForm::join normalizes each operand and substitutes
the global Undetermined state on contradiction. That fallback disappears with
the state. The per-kind coset join also currently falls back to an undetermined
coset; these are distinct behaviors, not proof that one can mechanically replace
the other.

Settle the stereo form's same-kind and cross-kind behavior when either or both
raw operands are contradictory, including entity constraint contradictions.
Relate it explicitly to the existing common-bottom normalized_eq contract and
join/normalization laws. Do not silently pick a kind, drop constraints, or change
the whole workspace's lattice algebra as a side effect of this migration.
This question concerns join, not the safety of the temporary replacement during
normalization or reframing.

### Coverage and measurement consequences

The existing suites provide the homes for coverage; add behavior cases there,
not new test modules named after implementation subitems or external fixture
collections. The affected surfaces include stereo form/DSL unit modules, edit
and delta units, molecule/reaction/span integrity and application tests,
canonicalization and symmetry tests, and Python stereo/edit/delta/reading tests.

- Replace kindless cases with explicit assertions about the removed syntax/API,
  not arbitrary tetrahedral defaults. Preserve all six kind domains, ordinary
  literal and unknown cosets, and deliberately incompatible raw attributes where
  the open-form contract permits them.
- Preserve normalization idempotence, raw-versus-normalized equality/matching,
  within-kind lattice laws, and explicit cross-kind outcomes. Add the decided
  contradictory-join cases and per-kind top law. Update generators and laws
  together rather than restricting draws until old tests pass.
- Cover update/difference reconstruction under its newly settled domain, direct
  known-kind edits with plain and journaled execution, inversion, and faithful
  serialization of constraint-only and no-op updates.
- Cover known-kind normalization and frame transport in forms, molecules, and
  spans: kind preservation, action/inverse and composition laws, idempotent
  reframing, and agreement of fused and tracked executions. Check surviving
  clones after both success and failure and both sides of Modified spans.
- Exercise same-kind reaction modification and different-kind explicit
  replacement, with compatible-degree kinds as well as different-degree kinds.
  Verify no preserved stereo pair crosses kinds, explicit replacement is not
  re-induced, and realized reactions/spans and compositions replay correctly.
  Cover supplied-match and enumerated paths under both matching algorithms.
- Migrate compact-string, tree, and streaming roundtrips for forms, updates,
  edits, reactions, spans, and standalone constraints. Existing entity-string,
  molecule, edit, reaction, and span fuzz targets must exercise the revised
  grammar. Preserve the roundtrip laws; remove only their obsolete kindless
  payload-loss exceptions.
- Check Python constructor/getter/setter parity, immutable retained Edit/Delta
  payloads, update constraint views, ordinary Python 3.9 access, and optional
  structural matching on supported versions. No Python-only stereo semantics.
- Update the existing constraint-readings benchmark to remove its obsolete
  kindless workload while retaining comparable known-kind cases. Use the bounded
  application/composition experiment already required here for a changed tracking
  route. Measure additional allocations or traversals if the chosen interface
  introduces them; no separate benchmark campaign is required just to replace
  mem::take. Structural simplification is not a measured speedup claim.

Permanent documentation changes belong to implementation: data-types.md for
mandatory kind and open-attribute boundaries, integrity.md for any approved
reaction construction rule, nomenclature.md for kind/configuration/coset and
Update/Delta distinctions, python-api.md for the actual binding changes, and
umol-dsl-spec.md for the complete syntax and lowering contracts. The historical
kindless decisions remain in their completed records. No manuscript edits or
whitepaper listing changes are part of this scope.

## Interfaces and representation contracts

### Dense IDs and transient editor state

Removal compacts dense IDs. A reference to a deleted bond cannot simply retain
its old integer: that integer can now name a different surviving bond. The same
problem applies to atom references. This is an ID-transport question, separate
from whether a higher-level reaction satisfies DPO.

The core prerequisite is settled: exact graph removal requires surviving
endpoints for every retained edge; participant compaction requires a surviving
image for every retained reference. Violating either prerequisite is a
programming error and may panic during required translation, as out-of-range
IDs may already panic. No separate dangling scan is required in these
primitives. They do not drop dependent rows, retarget references, or introduce
an invalid-reference representation.

This asserted contract is acceptable only if normal callers establish its
prerequisites. Checked graph removal does so through its check; cascading
removal through its complete selection; DPO application through pre-mutation
applicability checks and correctly ordered edits. Matching undo relies on its
history and reverse execution order. Each migrated call path must establish
the corresponding property; merely replacing the cascading backend does not
provide that guarantee. Existing non-panicking behavior for manipulated undo
data remains a separate obligation, handled locally.

Catching panics is not part of transaction semantics and must not be introduced
to implement this contract. Ordinary application rejection must not reach a
core assertion; a panic is not a substitute for a documented edit error.

The public boundary is settled. Edits are independently constructible inputs,
not provenance-bound recovery data. Applying an invalid removal must report
rejection; avoiding a panic while producing an unspecified result is
insufficient. Wherever cascading previously prevented missing references, the
updated exact operation checks reference coverage and survival before any mutation
by that operation. Editor drafts may already contain references outside the
current atom/bond domains; exact removal rejects those as InvalidReference,
distinct from DanglingReference for a valid reference broken by this removal.
This applies to direct editor removal and both plain and journaled
Edit execution, with the same semantics exposed in Python. Only after the
checks succeed may execution call the asserted core primitives.

Remove or redirect dependent rows before removing their stored references.
The checks cover all affected stores before any of them is changed; discovering
a dangling reference after deleting topology is too late. The current handle
and duplicate checks are retained but are not sufficient on their own.
Rejection leaves the failing removal unapplied. Earlier completed edits follow
the existing boundary: editor application retains its completed prefix, while
transaction application restores its entry state. This does not require
prevalidating an entire sequence against its initial state.

The shared removal error and method signatures below implement this boundary.
The relaxed misuse contract for derived undo data does not apply to Edit inputs.

This prerequisite concerns stored references, not every molecule integrity
condition: deleting a stereo ligand's adjacency bond can leave all stored atom
references intact, allowing an intermediate incidence inconsistency to be
repaired before publication.

### Agreed primitive inventory and naming

The following records the agreed interface.
It replaces the proposal to add relation compact_cascading and
tracked_compact_cascading. Cascade selection belongs to graph IR for overlays;
the selected rows use the existing exact-removal operations. The ordinary graph
cascade remains a graph-core operation over nodes and edges.

| Name | Meaning |
| --- | --- |
| remove | Delete the selected entities and close gaps in their own ID space. Relation row removal does not translate graph references. |
| tracked_remove | The same mutation, returning its compaction. Tracking does not mean journaling. |
| try_remove / try_tracked_remove | The existing checked graph counterparts: reject a dangling edge before mutation. No checked relation-row counterpart is proposed. |
| remove_cascading / tracked_remove_cascading | Graph removal that selects incident edges as well as the requested entities. |
| remove_topology / remove_topology_cascading | Graph-IR entry points distinguishing exact topology removal from selection of additional bonds and overlays. |
| compact_participants | Translate surviving graph references without removing or renumbering relation rows, reordering participant positions, or changing payloads. |

compact_participants complements restore_participants. It returns no mapping:
row IDs and participant positions are preserved. Existing Compaction::compact
continues to translate a single ID; it does not mutate storage. There is no new
apply_compaction operation or removal overload accepting a precomputed
Compaction or MoleculeCompaction.

#### Graph-core Graph

Bulk signatures:

```rust
pub fn remove(&mut self, nodes: Vec<NodeId>, edges: Vec<EdgeId>);
pub fn tracked_remove(
    &mut self,
    nodes: Vec<NodeId>,
    edges: Vec<EdgeId>,
) -> GraphCompaction;

pub fn try_remove(&mut self, nodes: Vec<NodeId>, edges: Vec<EdgeId>) -> Option<()>;
pub fn try_tracked_remove(
    &mut self,
    nodes: Vec<NodeId>,
    edges: Vec<EdgeId>,
) -> Option<GraphCompaction>;

pub fn remove_cascading(&mut self, nodes: Vec<NodeId>, edges: Vec<EdgeId>);
pub fn tracked_remove_cascading(
    &mut self,
    nodes: Vec<NodeId>,
    edges: Vec<EdgeId>,
) -> GraphCompaction;
```

Owned inputs amend the earlier slice proposal. The removal lists can become
the compaction's backing vectors instead of being copied. Input order and
repetitions retain their current treatment. This changes input ownership, not
the graph's mutable receiver or the meaning of rejection. A failed checked
removal leaves the graph unchanged; the supplied vectors have been consumed.
Do not add parallel borrowed/owned variants. Existing single-entity cascading
conveniences retain their names and delegate to the bulk operation.

The asserted exact operation requires live retained endpoints. It performs no
separate dangling scan. The checked operation establishes that prerequisite;
the cascading operation establishes it by selecting additional edges. Both
then use exact removal without repeating cascade discovery. Existing
out-of-range-ID behavior remains explicit. A missing retained-endpoint image
is an asserted programming error under the core contract above.

#### Graph-core relation storage

Apply the same interface to FixedRelationSet, VarRelationSet,
FixedFixedBirelationSet, FixedVarBirelationSet, and VarVarBirelationSet:

```rust
pub fn remove(&mut self, ids: Vec<RelationId>);
pub fn tracked_remove(&mut self, ids: Vec<RelationId>) -> Compaction<RelationId>;
pub fn compact_participants(&mut self, topology: &GraphCompaction);
```

The first two adapt existing methods to consume their removal lists. They
preserve survivor order, participant values, and payloads, and update relation
incidence. Their bounds and duplicate-ID contracts remain unchanged.
compact_participants mutates references and incidence only. It neither chooses
rows to remove nor adds a checked/tracked counterpart. Covering identity is a
no-op; a missing image is an asserted programming error under the core contract
above. Callers must establish coverage and survival before transport.

Retire the current borrowing compact/tracked_compact pair once its consumers
use explicit row removal and participant transport. Do not retain those
rebuilding cascades under new names by default. A future combined path would
need evidence that composition has a material cost, together with a settled
contract; it is not part of this initial surface.

#### Graph-IR owning storage and Molecule delegates

Keep typed slice inputs at the owning overlay boundary. For each store and ID
pair below, retain crate-private remove and tracked_remove, and replace its
row-dropping compact with crate-private compact_participants. AromaticSystems
illustrates the concrete signatures:

```rust
pub(crate) fn remove(&mut self, ids: &[AromaticSystemId]);
pub(crate) fn tracked_remove(
    &mut self,
    ids: &[AromaticSystemId],
) -> Compaction<AromaticSystemId>;
pub(crate) fn compact_participants(&mut self, topology: &GraphCompaction);
```

| Owning store | ID | Replacement Molecule delegate |
| --- | --- | --- |
| DativeBonds | DativeBondId | compact_dative_bond_participants |
| AromaticSystems | AromaticSystemId | compact_aromatic_system_participants |
| MulticenterBonds | MulticenterBondId | compact_multicenter_bond_participants |
| NoncovalentBonds | NoncovalentBondId | compact_noncovalent_bond_participants |
| StereoAtoms | StereoAtomId | compact_stereo_atom_participants |
| StereoBonds | StereoBondId | compact_stereo_bond_participants |

Each replacement delegate is crate-private, takes &mut self and
&GraphCompaction, and returns (). The explicit participant noun distinguishes
it from deleting overlay rows. Sites and ligand anchors are participants in
this storage sense. Existing Molecule bulk removal delegates retain their typed
slice inputs and return () or a typed Compaction; they affect their owning store
only.

Molecule's crate-private remove_topology and tracked_remove_topology likewise
retain typed atom/bond slices and () / GraphCompaction results. They become
exact graph-and-attribute operations. The topology coordinator performs
overlay reference transport and molecule-constraint compaction separately.
This preserves the distinction between the graph-only mapping produced here
and the complete molecule mapping assembled by the coordinator.

These wrappers already construct owned core-ID vectors; the agreed core
input change lets them transfer those vectors without a second input copy.
It does not eliminate conversion allocations or the separate conversion of a
core row compaction back to a typed graph-IR compaction. Those costs are an
investigation target, not a justification for unreviewed conversion APIs.

#### Graph-IR editor and edit batches

MoleculeEditor keeps typed atom/bond slices and a mutable receiver. All four
topology-removal methods are public:

```rust
pub fn remove_topology(
    &mut self,
    atoms: &[AtomId],
    bonds: &[BondId],
) -> Result<(), TransactionError>;
pub fn tracked_remove_topology(
    &mut self,
    atoms: &[AtomId],
    bonds: &[BondId],
) -> Result<MoleculeCompaction, TransactionError>;

pub fn remove_topology_cascading(&mut self, atoms: &[AtomId], bonds: &[BondId]);
pub fn tracked_remove_topology_cascading(
    &mut self,
    atoms: &[AtomId],
    bonds: &[BondId],
) -> MoleculeCompaction;
```

The exact pair reports failed removal prerequisites before mutation. The
cascading pair retains the existing explicit-cascade contract. Tracking changes
the returned result, not the mutation or its prerequisites; keeping the tracked
methods crate-private would not protect an additional invariant.

Make the six existing tracked overlay-removal methods public as well:

| Method | Input |
| --- | --- |
| tracked_remove_dative_bonds | &[DativeBondId] |
| tracked_remove_aromatic_systems | &[AromaticSystemId] |
| tracked_remove_multicenter_bonds | &[MulticenterBondId] |
| tracked_remove_noncovalent_bonds | &[NoncovalentBondId] |
| tracked_remove_stereo_atoms | &[StereoAtomId] |
| tracked_remove_stereo_bonds | &[StereoBondId] |

Each takes &mut self and returns MoleculeCompaction, with identities for
unchanged families. Preserve the existing mutation and bounds contracts and
the six public untracked counterparts. The owning-store and raw Molecule
helpers described above remain crate-private; public mutation goes through the
editor, which also maintains molecule constraints.

Use the existing TransactionError carrier for exact-removal rejection. Add:

```rust
InvalidReference {
    entity: Entity,
}
DanglingReference {
    removed: Entity,
    dependent: Entity,
}
```

InvalidReference names an unavailable referenced atom or bond, using the same
terminology as MoleculeIntegrityError::InvalidReference. DanglingReference names
a valid removed entity and the dependent entity that would retain a stored
reference to it. IDs refer to the state immediately before the rejected removal.
Direct editor removal returns these errors; Edit execution propagates them through
the existing MoleculeApplyError::Transaction wrapper. This introduces no new
error type or transaction behavior.

Edit::RemoveTopology and Edit::RemoveTopologyCascading share these fields:

```rust
atoms: Vec<AtomHandle>,
bonds: Vec<BondHandle>,
```

Edits::remove_topology and Edits::remove_topology_cascading accept those vectors
and return (). Existing remove_atom/remove_bond conveniences select
RemoveTopology. These methods construct edit data; contextual removal checks
belong to application. No cascading singleton variants are proposed. Plain and
journaled execution share the selected deletion and pre-mutation rejection
semantics; capture of recovery data is an execution concern.

### Graph-IR usage and mapping construction

| Workflow | Information available before removal | Required result or consumer |
| --- | --- | --- |
| Explicit overlay removal | IDs, or handles resolved and old values checked by Edit execution | Row compaction for constraints and subsequent edit handles |
| Exact topology removal | Complete atom/bond removal lists and successful pre-mutation reference checks; DPO callers also establish applicability | GraphCompaction for attributes and surviving participants; identity overlay row mappings |
| Cascading editor removal / extraction | Requested topology IDs; dependent removals require selection | Actual graph and overlay compactions; tracked extraction returns their aggregate |
| Plain Edit execution | The same selections, in the current handle namespace | Mappings are needed for subsequent handles even without returned provenance |
| Journaled execution | Explicit selections, or the selected cascade | Removed payloads before mutation, plus mappings and constraint recovery data |
| Undo of additions | IDs issued by the matching history | Exact removal under the restoration contract; no cascade discovery |

No current forward workflow starts with an independently supplied complete
row compaction. Requiring maps as results does not require precomputing a
complete MoleculeCompaction. For explicit removal, selected IDs already
identify the payloads to capture. For a cascade, select each dependent row
once and reuse that selection for capture and mutation.

The proposed cascade coordinator determines the complete bond removal list
against the original topology, then selects overlay rows through existing
incidence accessors. It must also retain the current selection of rows whose
references are outside the source domains, as recorded in S0b. It can process
overlay families in turn while that topology remains available: select unique
row IDs, capture if journaled, and perform exact row removal. After exact
topology removal, translate surviving participants, assemble the component
mappings, and update constraints and
handles. The operation uses existing typed ID lists and compaction carriers;
no additional public plan type is proposed.

Retain MoleculeCompaction as the common tracked editor result in Rust and
Python. The fixed bookkeeping for unchanged entity kinds is an accepted cost
of the simpler API: callers can translate any entity ID and pass the result to
consumers without assembling missing identity components. Each identity stores
its source count and an empty removal vector, with no heap allocation or work
proportional to the number of entities. The eight descriptors occupy 256 bytes
on a 64-bit target, excluding removal buffers. Lower storage layers retain
their narrower typed compactions.

Nonempty removal vectors still allocate storage. Keep the aggregate for
constraints, handles, recovery, and returned tracking. Do not add separate
mapping-free execution paths merely because a particular direct editor call
has no constraints. Acceptance of the fixed representation cost does not
accept unnecessary size-dependent work in consumers: S2f must still skip
unchanged handle families. Returning MoleculeCompaction does not require
expanding it into a MoleculeCorrespondence.

### Known duplicate work and remaining cost questions

Remove structurally identified duplication before comparing alternative
ownership/conversion implementations:

- Where exact removal replaces a protective cascade, move dependency discovery
  before mutation and reject an unselected dependent entity. Reuse existing
  incidence accessors and established checks instead of checking and then
  rediscovering the same dependencies in the mutation backend. Rejection can
  stop at the first failure. Successful execution still transports surviving
  references; moving the check does not make that separate work disappear.
- Establish source-domain coverage from the largest atom and bond IDs in each
  store's existing sorted incidence indices. This requires a fixed number of
  comparisons per store, without another scan of every participant, a new cache,
  or allocation. Wire access to that information through storage. Check deletion
  dependencies through incidence separately; coverage does not establish survival.
- Undo of added topology needs only the graph-removal panic guard described in
  S0b. Do not reuse the editor's broader reference check or scan overlays that
  undo neither reads nor transports.
- The checked Graph removal currently checks incidence and then rediscovers
  incident-edge removals in the cascading backend. Route it to exact removal.
- Journaled topology removal discovers dependent bonds and overlays for
  capture, then discovers them again during mutation. Reuse selected IDs.
- tracked_extract currently selects incident bonds before invoking a cascade
  that selects them again. Preserve extraction's selection semantics while
  assigning that work to one place.
- Relation row removal currently rebuilds incidence. Following it with another
  rebuilding participant operation would repeat that work. Stable compaction
  permits filtering and translating the existing sorted incidence arrays;
  preserve ordering and uniqueness without reconstructing and sorting them.
  Separate row and participant operations still perform separate traversals.
- Avoid copying an owned removal list only to construct the compaction, and
  skip unchanged handle families instead of materializing or traversing their
  tables under identity mappings. Skip participant writes for unchanged
  reference domains before forcing copy-on-write storage.

Two costs need separate decisions after that work:

1. **Cloning before removal.** Exact mutation through Arc::make_mut can clone
   rows that will immediately be deleted. The current borrowing cascade clones
   surviving rows only. Unique and shared storage therefore require separate
   accounting; an in-place design is not automatically faster.
2. **Allocations for ID conversion.** Trace typed graph-IR IDs to core IDs and
   core compactions back to graph-IR compactions. Count the vectors and
   constructor work actually performed. Consuming the core input removes one
   copy but does not establish that the remaining conversions are free.

The investigations below decide whether to keep the simple implementation or
change it. They do not presume a combined removal primitive, a new conversion
trait, a public Cow, or another ID representation.

### Scope of explicit cascading

The agreed baseline is to retain the existing cascade under an explicit name.
That cascade follows stored participant and site references. Removing a
stereo-bond site deletes its stereo bond; removing a stereo ligand's adjacency
bond can retain the stereo row and fail final incidence integrity. Removing an
internal aromatic bond retains the aromatic system. Thus the existing cascade
does not close over all bond dependencies specified for DPO application above.

Current overlay compaction also drops a row when a stored participant or site
is outside the source atom/bond domain, including in transient editor drafts.
Retain that missing-image behavior in the explicit cascade. Check bounds before
following a stored site; a selection based only on removed IDs would miss these
rows. Exact removal instead reports InvalidReference before transport.

Expanding the explicitly cascading convenience to those additional incidence
dependencies is outside this work. Compaction application must not introduce
that expansion implicitly. Reaction application remains strictly DPO.

### DSL and Python

Both Edit variants require faithful standalone-edit DSL serialization and
deserialization. Existing :topology/:remove and singular atom/bond removal
forms continue to select the exact operation. The cascading form uses :topology
with :remove-cascading. Cover both parsing entry paths and serialization.

Python already exposes MoleculeEditor.remove_topology and the six untracked
overlay-removal methods. Change remove_topology to the exact contract and add
these nine editor methods:

| Method | Input | Result |
| --- | --- | --- |
| remove_topology_cascading | atoms, bonds | None |
| tracked_remove_topology | atoms, bonds | MoleculeCompaction |
| tracked_remove_topology_cascading | atoms, bonds | MoleculeCompaction |
| tracked_remove_dative_bonds | ids | MoleculeCompaction |
| tracked_remove_aromatic_systems | ids | MoleculeCompaction |
| tracked_remove_multicenter_bonds | ids | MoleculeCompaction |
| tracked_remove_noncovalent_bonds | ids | MoleculeCompaction |
| tracked_remove_stereo_atoms | ids | MoleculeCompaction |
| tracked_remove_stereo_bonds | ids | MoleculeCompaction |

Inputs retain the existing integer-ID sequence conversion. These methods mutate
the editor without consuming it. Return the existing Python MoleculeCompaction
wrapper by transferring the Rust result; no new result type is needed. Exact
removal rejects invalid or dangling references through the existing
TransactionError exception. Preserve the existing consumed-editor and ID-range
checks and error translation. The bindings delegate dependency checks and
mutation to Rust.

Also add Edits.remove_topology_cascading and Edit.RemoveTopologyCascading, with
the same atom/bond handle fields as the existing exact variant. Existing
Edits.remove_topology and Edit.RemoveTopology adopt the exact contract. Builder
methods return None; checks occur when applying the batch. Tracking remains an
application operation, not another Edit variant or builder method.

Core graph/relation primitives and crate-private owning-store/transport helpers
do not gain separate Python bindings. Retain current handle coercion, immutable
Edit observations, ownership, and error translation. Add no Python-only removal
policy. Examples and tests must remain usable on Python 3.9; structural pattern
matching cannot be the sole way to use either variant.

## Development-guide updates

Update the normative guides with the implementing subitems, so they describe
available behavior. S2h checks the complete set for stale descriptions. These
are changes to existing sections, not a new guide or a duplicate API reference:

- **Data-type guide:** replace the Compaction section's implicit row-dropping
  compact/tracked_compact contract and private tracked-editor descriptions with
  the agreed removal and participant-transport families. State retained-reference
  and source-domain prerequisites, checked rejection before mutation, the core
  programming-error boundary, and public editor MoleculeCompaction results with
  count-bearing identities. In the mutation/publication and provenance sections,
  distinguish independently supplied edits from manipulated undo, and operation
  tracking from incidence induction. Replace the supplied-match application's
  claim that output pairings are induced: explicit deletion/recreation remains
  unmatched on the two sides, through realized results and composition. State the
  same replacement rule for different determined stereo kinds in application,
  span conversion, and reversal, including ordinary constraint-removal behavior
  and the final-side check in composition. Keep general induction and
  correspondence construction unchanged. Owned core lists,
  borrowed typed editor inputs, and unchanged-on-rejection guarantees belong to
  these method contracts. Owners: S1a/S1b, S2d/S2e, and S2h.
- **Nomenclature guide:** distinguish exact row/topology removal, explicit
  remove_cascading/remove_topology_cascading, and compact_participants as
  reference transport without row selection. A compaction records removals;
  it does not select a cascade. Clarify that a mutating tracked removal returns
  the witness alone because its bare result is unit. Apply the existing try_
  convention to checked/asserted core pairs; the ordinarily fallible editor
  operation keeps remove_topology. Tracking preserves surviving entity identity
  without promising an action journal. Owners: S1a/S1b and S2e/S2h.
- **Python API guide:** extend the editor contract with exact/cascading removal,
  all eight tracked removals returning the existing MoleculeCompaction, and
  untracked removals returning None. Distinguish direct integer entity IDs from
  batch handles, and explain that constructing a removal Edit does not check it
  against an editor. Record InvalidReference and DanglingReference rejection
  through TransactionError before that removal mutates state, retaining the
  existing completed-prefix, transaction recovery, consumed-object, and
  batch-transfer rules. Explain the new Edit variant and builder through Python
  3.9-compatible access; preserve Rust names and semantics without binding
  private core helpers. Owners: S2e/S2g.
- **Integrity guide:** distinguish the stored-reference precondition needed
  before dense-ID compaction from final molecule publication checks and
  host-dependent DPO applicability. Add the aromatic and stereo bond-dependency
  examples to the operation boundary, not the aggregate integrity inventory.
  Bond-order modification and bond addition remain permitted; their chemical
  acceptability is Tier 3 validation. Owners: S2d/S2e.

The existing test-writing and property-test policies already cover the required
verification. Add behavior cases and laws to the owning suites; no new testing
policy is needed. Keep public rustdoc self-contained and update examples with
the method or variant they demonstrate. DSL spelling and roundtrips belong to
the existing Edit serialization surface.

## Verification requirements

Verification is part of the implementation scope:

- Exact removal of explicitly listed atoms/bonds, survivor order, shifted IDs,
  and preservation of unrelated overlay rows, frames, and attributes.
- Relation and overlay compaction distinguish reference transport from cascading
  row deletion. Verify each family's ordinary and tracked operations against
  that contract; their returned mappings must describe the same removals.
- Cascading cases that retain the existing incident-bond, overlay, and
  molecule-constraint results. Existing cascade tests must continue testing
  cascades under the explicit name.
- The agreed boundary for references to deleted entities, including a deleted
  ID that would otherwise alias a different survivor after compaction.
- Invalid independently constructed removal edits return the specified error
  before the removal changes topology, overlays, attributes, constraints, or
  handles. Cover plain, journaled, direct editor, and Python entry points;
  verify the existing completed-prefix and rollback behavior separately.
- Drafts with unavailable overlay participants or stereo-bond sites report
  InvalidReference on exact removal; valid references broken by the requested
  deletion report DanglingReference. Check every incident row, including
  overlapping overlays in a draft. Explicit cascading retains the current
  missing-image row selection. Unrelated transient integrity failures remain
  publication concerns.
- Explicit dependency removal followed by topology removal, with plain/editor,
  tracked, and both transaction recovery paths producing the same result.
- Undo and correspondence laws, including shared storage, newly created
  entities, batch handle namespaces, and restoration after a later failure.
  Manipulated undo includes an incomplete incident-bond removal list and retained
  overlay references, not only out-of-range IDs. These cases require only absence
  of panics, not valid references, successful publication, repair, or a particular
  molecular result. Matching history retains its full restoration laws.
- Pre-mutation DPO cases for aromatic internal bonds, stereo-atom ligand bonds,
  stereo-bond site bonds, and stereo-bond ligand bonds. Include host-only
  dependencies absent from the pattern, explicit co-removal, and unrelated bond
  deletion. Do not replace these with final-product integrity tests.
- Identical applicability across all application result forms and iteration;
  faithful explicit replacement identity, realized-reaction/span replay, and
  reaction-composition behavior.
- Stereo atoms and stereo bonds with different determined kinds become Removed
  and Added in direct and realized spans, with no preserved correspondence pair.
  Cover span/reaction replay, reversal, constraint references to the replacement,
  and hosts with additional constraints that replacement must not inherit.
  Same-kind modifications and kindless refinement retain their contracts.
  Composition covers Th -> kindless -> Sp and explicit same-kind replacement;
  neither may acquire a preserved pair through composition or incidence induction.
- Rust/DSL/Python interface parity and roundtrips for both variants.

The aromatic cases from issue #31 require these explicit acceptance cases:

| Case | Expected application behavior |
| --- | --- |
| Delete an internal benzene bond; the pattern omits the host aromatic system | Reject the match; no products |
| Match the aromatic system but retain it while deleting an internal bond | Reject the match; mentioning the system is not explicit co-removal |
| Explicitly remove the matched aromatic system and its internal bond | Permit the deletion; the product contains neither removed entity |
| Delete and re-add an internal bond with the same endpoints while retaining the aromatic system | Reject the match; the final endpoint pair does not erase the original deletion |
| Delete an exocyclic C-H bond on an aromatic carbon | Permit the deletion; preserve the aromatic system |
| Delete the bridge between the two aromatic systems of biphenyl | Permit the deletion; preserve both aromatic systems |
| Change the order of an internal aromatic-system bond | Permit the modification; preserve the bond's identity and aromatic-system entry |
| Add a bond between two members of one aromatic system | Permit the addition; preserve the aromatic-system entry |

Exercise these cases with both GraphAndOverlays and Incidence matching. For
rejected cases, establish that the pattern has matches, then require Ok(None)
from supplied-match application and an empty product sequence from iteration,
without an error. Successful cases compare against explicit expected products,
including topology, changed attributes, and retained or removed aromatic-system
entries. Include the original unannotated bond pattern: rejection must not
depend on adding an aromaticity assertion or preventing aromatic atoms from
matching. The modification/addition cases assert application behavior, not
Tier 3 conformance of their products.

Extend the existing behavior-oriented unit and public-API property suites. Keep
their laws intact when migrating generators from implicit to explicit cascades.
Do not add a separate issue-specific test module or external fixture collection.
Add generated preservation checks for the new boundaries: participant transport
preserves relation rows and payloads while translating every reported reference;
exact removal agrees with cascading removal when all dependencies are explicitly
selected; a rejected exact removal leaves its entry state unchanged; and an
explicitly replaced entity has no surviving host-to-product pairing, including
through realized-result replay and composition. Exercise these through public
operations with independent expected mappings and the stated prerequisites.

Use focused crate checks while implementing. At the completed related-work gate,
run formatting, workspace tests and Clippy, the affected feature-gated property
suites, and Python tests after rebuilding with the repository's Python 3.13
environment. Full workspace checks do not run after each small subitem.

Performance measurements must answer an acceptance question. Retain the original
implementation baseline for the total cost comparison. Before switching graph-IR
removal execution, establish the complete DPO checks and faithful tracking on the
existing cascading backend, with the new behavior tests passing. Preserve that
green implementation as the reference for comparing removal integration with the
same checks and tracking on both sides. Record any earlier core changes in that
reference; its difference from the original is not solely a semantic cost.
Use previously successful cases that remain applicable for timing comparisons;
newly rejected matches are correctness cases, not faster equivalent workloads.
Include complete application and composition to assess the tracking trial,
and complete removal operations to avoid mistaking a local saving for a total
one. Check for repeated traversal or allocation introduced by reference transport
and witness construction. Broaden measurement only if a material regression or
an unresolved cost justifies it.

## Planning readiness

Exact versus cascading removal, pre-mutation DPO checks, separate row removal
and participant transport, and the migration obligation are settled. The
primitive inventory records agreed names, public/editor versus private/storage
visibility, signatures, ownership, shared error, and DSL/Python surfaces. The
core permits programming-error panics for violated prerequisites
on retained references, conditional on normal callers establishing them;
transaction panic recovery is excluded. Independently supplied edits and direct
editor removals reject invalid removal before mutation; Python follows that
contract. Semantics, interfaces, and nomenclature (SIN) are settled for the
removal scope, including the common MoleculeCompaction return type and its fixed
bookkeeping cost. Faithful replacement tracking remains an implementation
obligation. Existing cascade scope is the baseline; no broader cascade is
implied by the strategy.

The stereo identity rule and its public interfaces and nomenclature are also
settled. S0b records the implementation route for replacement lowering, span
construction, application tracking, and composition; S2d implements and measures
it. The separate stereo experiment is outside the implementation strategy below
and does not block it. No removal of kindless forms, new update types, or
resolution of the historical SIN1-SIN5 is required. Only an explicit decision to adopt the broader redesign
would justify integrating its implementation and dependencies into a plan.

S0's code analysis is complete: caller prerequisites, comparison boundaries,
and replay implementation routes are recorded below. Implementation proceeds
through the subitems below.

## Bottom-up implementation strategy

Introduce the core capabilities first. Then compose them in graph IR and remove
known duplicate work. Investigate cloning and ID conversion against that
implementation, separately, so that each measurement answers a remaining
design question. Do not compare a half-migrated implementation or attribute
the cost of new DPO checks to a storage change.

Each implementation subitem carries focused unit cases and applicable public-API
properties in the existing suites. Additive transport methods coexist with the
old delegates until their consumers can migrate. A breaking signature or variant
change includes its caller migration and exhaustive DSL/Python conversions in
the same subitem; it ends green. Every stage ends green. Full workspace gates
run at S5, not after each subitem. No stage is deferred.

### S0 — Caller prerequisites and comparison boundaries

#### S0a — Core and storage contracts (done) [dep: none]

Completed 2026-10-17 by code analysis. The agreed primitive inventory above fits
the current storage and callers; no additional public primitive or failure type
is needed. This subitem changes documentation only. Measurements and the
graph-IR publication/replay audit remain in their assigned subitems.

**Core prerequisites.** Compaction is an open, count-bearing carrier, not proof
that a separately supplied graph or relation set is compatible with it.

| Operation | Prerequisite and producer guarantee |
| --- | --- |
| Graph exact remove/tracked_remove | Supplied IDs belong to the current graph. Every retained edge has surviving endpoints. Node and edge source counts are captured from that graph before removal; they are different domains. Unsorted and repeated removal IDs retain their existing deduplication behavior. |
| Graph try_remove/try_tracked_remove | The same ID-domain prerequisite applies. The existing incident-edge check establishes retained-endpoint survival before mutation; failure remains None with the graph unchanged. Successful execution can delegate to exact removal. |
| Graph remove_cascading/tracked_remove_cascading | The same ID-domain prerequisite applies. Collecting all edges incident to removed nodes establishes retained-endpoint survival. Pass the completed deletion lists to exact removal without discovering those edges again. |
| Relation remove/tracked_remove, all five families | Selected row IDs belong to the current relation set. The row compaction's source count is that set's pre-removal row count. Removal preserves surviving rows' order, participant values and positions, and payloads; it does not translate their references. No graph is required for row removal. |
| Relation compact_participants, all five families | Every node or edge reference reported by a retained participant is inside the corresponding GraphCompaction source domain and has a surviving image. Transport changes references and incidence only. Row IDs, participant positions and multiplicity, and payloads remain unchanged. |

Relation storage owns no graph or graph-domain counts. Its transport prerequisite
is coverage and survival of the references it actually contains, not equality
between graph counts and relation counts, nor equality with the largest used ID
plus one. Unused graph IDs may be removed. A covering identity is a no-op;
zero-domain identity is not a universal identity. RelationParticipant::refs must
report every transported reference, including both components for a participant
that refers to a node and an edge. StereoLigand transports its anchor atom and
retains its ligand kind.

Compaction::new already validates removal bounds before storage mutation.
The owned Vec inputs can supply its removal storage; ownership does not add a
new validity condition. Graph-core assertions express programming errors at
these boundaries. They do not replace checks required by normal public callers,
and transaction recovery does not catch panics.

**Caller audit.** The following records the production callers of the affected
graph and relation removal/compaction primitives before S1.

| Current caller | Established prerequisite and migration |
| --- | --- |
| [Graph removal delegates](../umol-graph-core/src/graph.rs) | Checked removal tests retained endpoints, then currently calls the cascade and repeats discovery. Cascading removal completes the edge list itself. S1a routes both through exact removal after their respective prerequisite is established. Singleton cascading conveniences continue through the cascade. |
| [Graph::tracked_pushout_complement](../umol-graph-core/src/rewriting.rs) | Correspondence composition checks the shared interface domain; host-side counts are checked against the supplied graph. Deleted IDs come from those validated correspondences. The checked graph removal still establishes the dangling condition. Keep that checked route and move its already-owned deletion vectors into it. pushout_complement delegates here. |
| [Molecule::tracked_remove_topology](../umol-graph-ir/src/ir/molecule.rs) | The sole production graph-removal bridge in graph IR converts typed IDs to node/edge vectors, obtains a compaction from the same graph, and compacts equally sized atom/bond attribute arrays. S1a transfers those vectors into the existing cascade. The later exact switch requires the graph-IR coordinator to supply the complete deletion and establish retained-reference survival first. |
| [MoleculeEditor](../umol-graph-ir/src/ir/molecule/editor.rs), plain and journaled [Edit execution](../umol-graph-ir/src/ir/molecule/apply.rs) | These call the raw Molecule bridge, then all six row-dropping overlay compactions. Edit execution first resolves handles and rejects duplicates; journaled execution additionally captures removed rows. Today the cascades establish survivor validity. S2 must establish it before mutation by checking exact removal or selecting the complete existing cascade, remove overlay rows explicitly, remove topology, then transport surviving participants. |
| Undo::RemoveAddedTopology in the same apply module | Matching reverse history removes later dependent additions before the original added topology. Current replay filters out-of-range IDs and uses the graph cascade. A mechanical switch to asserted exact removal would not cover manipulated history. S0b records the local guard before that caller migrates; ordinary history must not gain panic recovery. |
| The six graph-IR owning overlay stores | Explicit remove/tracked_remove builds RelationId vectors for its own rows; pass those vectors directly to the changed core signatures. Existing compact calls borrowing tracked_compact, which selects surviving rows by successful participant translation. That selection must move to the coordinator before replacing it with compact_participants. The topology compaction retains the original graph domains while the surviving rows still hold original participant references. |

The owning-store mapping is FixedRelationSet for NoncovalentBonds;
VarRelationSet for AromaticSystems and MulticenterBonds; and
FixedVarBirelationSet for DativeBonds, StereoAtoms, and StereoBonds. The latter
stores its BondId site in the edge factor; stereo ligand anchors use the node
factor. FixedFixedBirelationSet and VarVarBirelationSet currently have no owning
molecule store, but retain the same public core contract. Their existing core
tests and benchmarks remain consumers of the signature change.

The graph, relation, rewriting and restoration unit/property suites, graph-IR
storage tests, and existing graph/relation benchmarks also call these methods.
S1 migrates their argument ownership without changing what their tests exercise:
cascade cases still request cascading behavior, checked dangling cases still
check rejection, and participant-only cases must explicitly remove dependent
rows before transport. New exact-operation cases are separate from those
existing laws.

**Reproducible comparison boundary.** The original production baseline is
commit 1f782a5f1147e7c70165bd4e84ba1fbae998d33c. At the S0 audit, working-tree changes
were discussion records, so that commit contains the audited implementation. Apply
the same measurement harness to both revisions, use separate Cargo target
directories, and record the toolchain, build flags and commands with the
results. Do not reuse timings from the broader stereo redesign experiment as
evidence for this narrower change.

Use these bounded fixtures; their sizes are controlled comparison cases from
existing benchmarks, not claims about typical molecules.

| Case | Construction and measured operation |
| --- | --- |
| Exact and cascading topology removal | Reuse the 8- and 64-node paths in [graph.rs](../umol-graph-core/benches/graph.rs), also used as carbon/single-bond paths by benchmark_mutation in [reaction.rs](../umol-graph-ir/benches/reaction.rs). Remove middle node m. Exact removal supplies edges m-1 and m; cascade supplies no edges. Both yield the same survivors. Keep checked-to-checked comparisons distinct from the asserted backend; use core timings only to attribute a whole-operation cost when needed. |
| Overlay removal and participant transport | Reuse removal_entries in [editor.rs](../umol-graph-ir/benches/editor.rs): five atoms, bonds (0,1), (0,2), (0,3), (1,4), and one row of each overlay kind. Remove atom 0. The existing cascade removes bonds 0–2 and the dative, aromatic, multicenter, stereo-atom and stereo-bond rows. The noncovalent row over atoms 2 and 4 survives and must be transported. Exact execution removes those five rows explicitly, then atom 0 and bonds 0–2. Compare complete operations with equal final contents. |
| Journaled cascading removal | Use the same mixed-overlay fixture and deletion through a Journal transaction that commits. Include removed-row capture, removal, reference/constraint/handle transport, journal work and publication. Rollback is a correctness requirement, not an additional timing matrix here. |
| Supplied-match application | Reuse application_case in reaction.rs, including its six-atom aromatic participant-frame change. Reuse the stereo experiment's five-atom carbon star and 128-atom extension: site 0, actual ligands 1–4, LHS tetrahedral undetermined coset, RHS coset 1, host coset 0, and a chain from atom 4 through 127 in the larger host. Supply the fixed identity match on the LHS atoms. Time tracked_apply_at and apply_at_to_reaction_span, including complete result construction. |
| Deletion-heavy application | Use the 64-atom carbon path as both LHS and host with a supplied identity match. Delete atoms 0–31 and bonds 0–31, including the boundary bond (31,32). This previously applicable case checks the cost of initializing tracking pairs that removal then discards. Use the same two application result forms. |
| Composition | Reuse the stereo experiment's C/F/Cl/Br/I star, site 0, ligands 1–4: tetrahedral coset 0 changes to 1, followed by its reverse. Measure complete compose with ModularProductBacktracking, including overlap enumeration and all results. Hold the application implementation fixed when comparing atom-only reconstruction with complete-correspondence composition. |

For removal, distinguish unique storage built afresh in untimed setup from shared
storage obtained by cloning a retained seed. A clone-based setup is not a unique
case even if the cloned Molecule variable has one owner. Use the same ownership
condition on both sides. Time prerequisite checks, removal, required transport
and publication together at graph IR; fixture construction stays outside. The
existing editor benchmark excludes finish, so its old timings are not complete
operation measurements. Keep input-vector construction at the same boundary in
core comparisons; end-to-end graph-IR timings include their ID conversions.

Check expected products and correspondence pairs outside the timed loop before
comparing costs. Different-kind replacement and newly rejected dangling matches
are correctness cases; manipulated undo is a no-panic case. Changed acceptance
or missing results must not count as faster execution. S2d preserves an
additional green reference with corrected DPO checks and tracking on the
cascading backend before the graph-IR removal switch. Subsequent storage comparisons use that reference
as well as the original baseline; any earlier core changes in the reference are
recorded. Run the relevant cases at each decision, not a cross-product of all
fixtures, result forms, relation families and rollback modes. No measurements
were run for S0a.

#### S0b — Graph-IR publication and replay boundaries (done) [dep: S0a]

Completed 2026-10-17 by code analysis against the S0a baseline. The existing
publication boundaries and correspondence interfaces support the agreed design.
The additional TransactionError::InvalidReference variant above was approved
for pre-existing invalid references in editor drafts. No production code,
tests, or benchmarks changed or ran for this subitem.

**Removal and publication.** The relevant owners are ir/molecule/editor.rs,
editor/transact.rs, molecule/apply.rs, molecule/transact.rs, and the Python
transaction/error bindings.

| Entry point | Required behavior |
| --- | --- |
| Direct editor exact removal | Check reference coverage and survival before changing any store. Return TransactionError and retain the editor on rejection. Preserve existing removal-ID bounds and duplicate handling. |
| Plain or tracked editor apply | Resolve handles and reject duplicate removals, then use the same reference check. A rejected removal changes nothing; preceding successful edits remain. Tracking changes only the returned witness. |
| Journal or clone transaction | Use the same check before mutation and before journal capture. An edit failure aborts and restores the transaction's entry state through the selected recovery method. |
| probe / finish / commit | Retain the final integrity gate. probe leaves the editor available; finish consumes it on either outcome; failed commit rolls back. Removal checks do not replace publication validation. |
| Python editor and batch methods | Retain consumed-object, ID-range, and batch-transfer checks. Both new error variants use the existing TransactionError exception; publication integrity failures retain their existing translation. |

The editor/Edit preflight has two parts. Check source-domain coverage from the
largest atom/bond IDs in the existing sorted incidence indices; valid removal
IDs alone do not exclude references to not-yet-added entities in a draft. This
is a fixed number of comparisons per store, not a scan of every stored reference.
Then check dependencies of the requested deletions through incidence. Inspect
stored IDs before following a stereo-bond site. Use all rows returned by
owning-store incidence, not the single aromatic-system lookup valid only after
publication: a draft can temporarily contain overlapping systems. Do not call
full molecule integrity validation here; unrelated incidence or uniqueness
failures can still be repaired before publication. Successful exact execution
uses asserted core removal and participant transport, without cascade discovery
afterward.

Explicit cascading must include out-of-domain participant/site references in
its selected rows, preserving the old compaction behavior. Reuse the selected
IDs for journal capture and removal. Source-domain coverage uses the same index
bounds; only a failing bound requires finding the affected rows. Dependency
selection and surviving-reference transport still cost work. S2e must account
for those traversals rather than claim that moving checks eliminates them.

The new reaction DPO checks inspect dependencies of deleted bonds; they do not
add molecule-wide chemistry validation. The core asserted primitives gain no
separate validation pass, publication retains its existing checks, and Python
translates the Rust failures. Checking before mutation can add traversal to
normal removal even though the old cascade already examined dependencies; no
net cost claim is established until the planned operation measurements.

**Undo.** Undo::RemoveAddedTopology currently filters out-of-range IDs and calls
raw topology removal. Matching reverse history removes later dependent additions
first and restores any later compaction before reaching the original addition.
The added topology is again at the end of its stores; removing it does not shift
retained IDs, so this path needs neither overlay transport nor constraint
compaction. For manipulated history, retain bounds filtering and duplicate
normalization, then check only whether a retained graph bond would reference a
removed atom. Use the deleted atoms' graph incidence; skip the undo entry if
that check fails. This replaces graph cascade discovery already performed by
the current path. Do not run the editor's overlay-reference preflight: raw
topology removal neither reads nor transports those references, so checking
them serves no panic-prevention purpose. Manipulated history need not produce
a correct molecule. Do not validate or repair the history, return a new error,
or catch panics.
RestoreRemovedTopology keeps its existing graph, participant, overlay-row, and
constraint restoration order.

**Reaction producers and replay.** The owners are ir/reaction.rs, ir/delta.rs,
ir/reaction_span.rs, ir/compose.rs, and the molecule correspondence/difference
implementations. S2d follows this order:

1. Fold the rule's resulting stereo form before reframe_deltas substitutes
   concrete host values for old values. Decide replacement from the actual host
   and resulting kinds, including a kindless LHS matched to a determined host.
   Different determined kinds use ordinary remove/add edits with the rule's
   resulting attributes, not a field edit retaining host-only constraints.
   Include those removals in DPO applicability checks. Same-kind modifications,
   kindless refinement, and generic Edit tracking retain their contracts.
2. Give each replacement a fresh created handle. Constraint removals refer to
   the old entity; explicit RHS additions refer to the replacement. A single
   overwritten handle map cannot serve both sides. Keep dependent overlay
   removals before topology removal and preserve the existing constraint
   compaction behavior, including dropping a whole molecule constraint whose
   referenced entity is removed. Do not introduce subtree repair or automatic
   retargeting of old constraints.
3. Use host.clone().tracked_apply for publication and the complete witness.
   ApplicationState already distinguishes surviving original entities from new
   ones in all eight families, including same-incidence explicit replacements.
   Remove output incidence induction; retain supplied-match validation and
   stereo reframing. A host dangling violation returns Ok(None); iteration
   skips it. Preserve existing error propagation for other failures and retain
   the final integrity gate. No additional clone or journal is needed.
4. In direct span construction, retain the LHS prefix as Removed entries and
   append Added replacements, allocating fresh IDs without collisions with
   explicit additions. Construct the RHS reference map before constraint spans;
   unchanged old references follow removal semantics, while explicit new
   constraints use replacement IDs. Realized spans and reactions consume the
   full witness through superimpose/difference_to. Published sides, preserved
   incidence, and the reaction kind rule establish those consumers' prerequisites;
   preserve their existing interfaces. Reversal continues through the span.
5. In composition, reverse the left application's full witness and compose it
   with the right witness. Their intermediate counts agree because both use the
   same glue molecule. Remove final stereo pairs with different determined
   kinds before difference_to, including Th -> kindless -> Sp. Construct the
   result with Reaction::new and retain delta normalization and existing
   inadmissible-overlap handling. General correspondence composition does not
   gain a stereo policy or new fallibility. Measure this route separately with
   the application changes held fixed, as already required by S2d.

The for_update correction takes constraint context from the resulting
configuration. prepare_deltas drops only old/new and LHS/new kind-equality
requirements; old/LHS compatibility, coset domains, expression actions, and
site/arity integrity checks remain. Removing the kind prohibition does not make
an invalid stereo form or site admissible.

**Coverage accounting.** Extend the owning molecule/edit/transaction, reaction,
span, composition, property, and Python suites. The current supplied-match
same-kind explicit replacement test expects an induced preserved pair; change
that expectation to unmatched. Keep generic kind-changing Edit tracking paired.
Split valid reaction kind changes from existing invalid-site/arity cases instead
of making all old negative cases succeed. Add replacement-ID collision cases,
old-reference removal versus RHS addition, host-only constraint removal, and
kindless-LHS application to a determined host. Preserve the existing
normalization, result-form agreement, span/reaction replay, reversal, and
composition laws. The reference-error and manipulated-undo cases are specified
in the verification requirements above. No new test module or public testing
seam is needed.

### S1 — Graph-core primitives

#### S1a — Graph removal family (done) [dep: S0a]

Completed 2026-10-17. Graph now exposes exact remove/tracked_remove; all six bulk
removal methods consume node and edge vectors. Checked removal establishes
endpoint survival before exact execution, and cascading removal supplies the
complete incident-edge list. Exact execution asserts surviving endpoints during
translation without a separate dangling scan. The supplied vectors become the
compaction's removal storage. Callers, examples, benchmarks, and tests use the
owned signatures; graph IR still calls the cascade. The data-type and
nomenclature guides describe these implemented contracts.

Verification passed: 1,906 graph-core unit tests, 11 graph properties, 38
restoration properties, and the graph-core doctest. New exact-removal cases
cover survivor IDs/endpoints, loops and parallel edges, shared storage,
identity, duplicate/unsorted inputs, and programming-error prerequisites.
Existing checked/cascading and restoration laws retain their meaning.
Clippy passed for graph-core with all targets and proptest enabled;
cargo check passed for graph IR with all targets. Targeted nightly formatting
and git diff --check passed. Exact-removal benchmark cases are prepared, with
owned inputs constructed outside timing; no performance measurements were run.

The shared execution route and input-vector reuse also complete that part of
S2b. Its remaining duplicate-work review is still open.

#### S1b — Relation removal and participant transport [dep: S0a]

Modules: the five graph-core relation-set families. Adapt remove/tracked_remove
to the agreed owned input and add compact_participants. Breaking red-to-green:
migrate callers of the changed removal signatures together; retain the old
compact/tracked_compact temporarily until S2 moves their consumers. Each family
keeps its own tests in the existing organization. Exercise row/payload
preservation, fixed and variable factor shapes, participant order and
multiplicity, node/edge references, complete participant values, incidence,
identity, and transport/restoration roundtrips. No new combined public removal
operation is introduced. Cover the asserted missing-image boundary separately
from normal producer/consumer sequences that must satisfy it.

S1 ends with focused graph-core tests and its feature-gated property target
passing, and all mechanically affected callers building. It establishes the
capabilities; S2 completes their integration and removes the remaining repeated
work before the cost investigations.

### S2 — Integration and removal of known duplicate work

#### S2a — Incidence maintenance [dep: S1b]

Module: umol-graph-core relation incidence and its storage callers. Internal
change; green. Filter/renumber incidence entries for row removal and translate
their graph keys for participant compaction, using monotonicity instead of
rebuilding and sorting the index. Verify with independent scans through public
incidence accessors across all five families. Account for the two traversals
when the operations compose; do not describe them as a fused pass.

#### S2b — Graph selection and execution [dep: S1a]

Modules: umol-graph-core graph and rewriting. Internal change; green. Route
checked removal to exact execution after its check; route cascading removal
after its selection. Eliminate repeated cascade discovery and avoidable copies
of the removal lists in that execution. Preserve rejection without mutation,
returned mappings, and pushout-complement behavior. Use the existing graph and
rewriting properties to verify the shared execution path.

#### S2c — Owning overlay stores and Molecule transport [dep: S2a, S2b]

Modules: graph-IR overlay stores and ir/molecule. Add the agreed typed
participant methods and six participant-only Molecule delegates; additive,
green. Retain the old cascading compaction delegates and the raw topology
helper's existing behavior until S2e migrates their consumers together. Keep
exact row removal in the owning stores. Cover typed mappings, source counts,
frame/payload preservation, and unchanged shared source molecules. Wire the
constant-cost source-domain coverage check from the existing sorted incidence
indices for S2e; do not add a participant scan, cache, or allocation for it. At
this point use the simple existing copy-on-write model; do not preselect the
S3 answer.

#### S2d — DPO checks and faithful application tracking [dep: S0b, S2c]

Modules: graph-IR reaction application, delta update lowering, reaction-span
construction and reversal, realized results, and composition.
Internal change; green. Complete the pre-mutation bond dangling checks while
the current cascading graph-IR removal backend is still in use. After these
checks, its cascade must select no additional entity.

Align reaction lowering and direct span construction first. Different determined
stereo kinds use ordinary replacement semantics in application and Removed/Added
entries in spans, with fresh replacement IDs, correct RHS references and frames,
and the existing removal/addition behavior for constraints. Do not inherit extra
constraints from the host's removed stereo entity. Include the accepted delta
for_update correction and removal of the enumeration-only kind-change prohibition.
Retain existing checks for invalid forms and applicability; a valid kind change
must no longer fail direct span construction or reversal solely because the
kinds differ. Keep shared public update/delta interfaces and generic edits.

Use tracked_apply for application correspondence, preserving explicit replacement
identity across all entity kinds and applying the reaction stereo identity rule.
Then trial complete correspondence composition, apply the rule to its final
stereo endpoints, and construct the composite through difference_to and
Reaction::new. Acceptance of the composition implementation requires the separate
comparison already specified. No new public constructor or tracking type is needed.

Cover host-only dependencies, every affected overlay kind, explicit removal
and replacement, whole-match rejection, all application result forms, replay,
direct span conversion, reversal, and composition. Include both stereo entity
kinds, host-only constraints, replacement constraint references, same-kind
modification, kindless refinement, and composition through a kindless intermediate.
Include the issue #31 acceptance cases above under both match
algorithms, including permitted bond-order changes and additions. Their chemical
acceptability remains a Tier 3 validation concern. Add the public-API preservation
properties for replacement identity and replay. These tests must pass before
the removal backend changes; keep them unchanged through that migration.

Record this green implementation as the reference with corrected semantics.
Run the bounded application/composition comparison selected in S0a to assess
the tracking route's complete cost, accounting for the added DPO checks and any
earlier core changes. Identify allocation or traversal regressions before
accepting the route. Preserve both this reference and the original baseline for
the subsequent comparisons; no separate permanent backend is introduced.

#### S2e — Editor and edit execution [dep: S2c, S2d]

Modules: graph-IR molecule/editor, molecule/apply, ir/edit, exhaustive DSL/Python
Edit conversions, and existing callers. Breaking red-to-green within this
subitem. Change the raw Molecule topology helper to exact execution together
with editor, plain/journaled application, extraction, and undo callers; switch
forward removal consumers to participant-only transport. Undo of added topology
retains its narrower path described in S0b. Implement the exact/cascading
topology split and the two Edit variants. Expose both tracked topology methods
and all six tracked overlay-removal methods publicly with the agreed
MoleculeCompaction results.
Check all exact-removal prerequisites before mutation or journal capture.
Report TransactionError::InvalidReference for unavailable stored participants
or sites and DanglingReference for references broken by the requested removal.
Use the incidence-index bounds for source-domain coverage; do not scan all
overlay participants to establish it. Reuse the dependency discovery previously
performed by cascading; do not add a redundant cascade scan after checking.
Select cascaded bonds/overlays once, including rows with out-of-domain references
in editor drafts; reuse selected IDs for mutation and journal capture, and
assemble MoleculeCompaction from the actual results.
Compact constraints once for the completed topology operation;
explicit overlay edits retain their existing separate constraint handling.
Migrate extraction and its generated edits to the cascading operation.
Remove duplicate incident-bond selection between extraction and that operation.
Check plain/editor/journal equivalence, subsequent handles, explicit co-removal,
rejection before mutation, rollback, undo of additions, and manipulated-history
behavior. Invalid Edit input has a reported failure, not the unspecified-result
contract of manipulated undo. Before Undo::RemoveAddedTopology reaches asserted
exact removal, retain ID bounds handling and check only the retained graph
endpoints as described in S0b. A failed check leaves that undo entry unapplied.
Do not add overlay checks, history validation, repair, or fallibility to undo.

The new Edit variant must include its DSL tree/streaming parsing and serialization,
Python Edit construction/conversion/observations, and all exhaustive matches in
this subitem. Adjust the existing Python exact-removal entry points for their
fallible Rust contract and test their rejection behavior. Carry focused variant
roundtrips and existing-surface tests with these changes; S2g adds the remaining
Python methods. Keep examples, generators, and other removal callers on their
intended exact or cascading behavior. Reaction lowering now reaches exact
removal with S2d's applicability checks already established.

#### S2f — Handle transport [dep: S2e]

Module: graph-IR molecule/apply ApplicationState and HandleTable. Internal
change; green. Skip unchanged compaction components instead of materializing
or traversing their handle tables. Preserve ID/New handle semantics, removal
errors, and returned correspondence. Test through edit application and its
public tracking results; do not expose internal tables for testing.

#### S2g — Python removal methods [dep: S2e, S2f]

Modules: Python edit/editor bindings and the Python API guide. Additive, green.
Add the nine editor methods listed above and the cascading Edits builder method;
the new Edit variant and existing exact entry points were migrated in S2e.
Return the existing MoleculeCompaction wrapper from all eight tracked removals
and preserve Rust ownership and failure semantics. Cover direct/batch use,
tracked/untracked equivalence, mapping source counts, and rejection before
mutation for both reference errors. Keep Python 3.9-compatible access paths.
Rebuild the extension before Python tests.

#### S2h — Retire implicit relation cascades [dep: S2e, S2g]

Modules: graph-core relation storage, remaining callers, and normative guides.
Breaking red-to-green. Remove the old compact/tracked_compact pair after all
production callers have migrated. Preserve the laws exercised by cascade and
restoration tests by expressing their selection explicitly; do not weaken them
to ordinary row-removal checks. Complete the development-guide inventory above,
including the data-type, nomenclature, Python API, and integrity guides; reconcile
rustdoc and examples with the final interfaces. Audit for hidden cascade
selection and for new duplicate incidence/capture work.

S2 ends with the full topology-removal distinction implemented across its
consumers and the affected unit, property, and Python suites passing. The two
remaining cost investigations start from this coherent implementation, with
the same DPO checks and cascade scope in both sides of each comparison.

### S3 — Decide the shared-storage cloning tradeoff

#### S3a — Isolate cloning before removal [dep: S2h]

Modules: graph-IR owning storage and relevant existing benchmarks. Investigation;
green. Account for payload clones and allocated storage separately for unique
and shared inputs, at deletion fractions chosen to expose the stated tradeoff.
Distinguish overlay copying from the existing atom/bond attribute copying.
Measure a bounded end-to-end removal/extraction workload as well as the storage
step. A narrowly scoped survivor-copying alternative may be prototyped to test
the hypothesis; it does not gain a public API merely to make the experiment
convenient. Keep selection, mappings, constraints, and output semantics equal.

#### S3b — Keep, change, or reject [dep: S3a]

Record the measured cost and a concrete decision. Keep ordinary copy-on-write
if its cost is acceptable; otherwise select and implement the narrow ownership
change supported by the comparison, or reject the candidate design. Any changed
primitive or failure boundary returns to interface review before implementation.
For an accepted storage change, rerun the affected preservation/shared-source
tests and the paired comparison. Retain the evidence needed for the decision,
not an open-ended benchmark worklist.

### S4 — Decide the ID-conversion allocation tradeoff

#### S4a — Account for conversions [dep: S3b]

Modules: graph-IR owning stores, molecule topology removal, and core Compaction.
Investigation; green. Trace both directions: typed atom/bond/overlay IDs into
core IDs, then core row compactions back into typed mappings. Include their
actual consumers, so an allocation eliminated at one boundary is not merely
recreated at the next. Separate vector allocation/copying from bounds checks,
sorting, and deduplication. Measure the conversion contribution and a bounded
whole operation after S3's decision. Do not assume iterator collection reuses
an allocation on every compiler or that type conversion is allocation-free.

#### S4b — Keep, change, or reject [dep: S4a]

Record which remaining allocations are justified and whether an ownership or
conversion adjustment is worthwhile. Specify any changed constructor,
conversion, visibility, and failure contract before implementing a candidate;
preserve typed ID domains and compaction invariants. Do not assume a new generic
conversion trait, representation cast, or removal variant is necessary. Verify
source/result counts, removed IDs, every survivor image, and roundtrips through
the public consumers, then repeat the bounded whole-operation comparison.
An explicit decision to retain the measured cost completes this investigation.

### S5 — Combined acceptance and closeout [dep: S3b, S4b]

Review the final diff against the agreed symbol inventory, semantic boundaries,
and caller audit. Run one combined acceptance comparison where needed to check
that the selected cloning/conversion changes compose without a material
regression. Report absolute and relative costs and measurement uncertainty;
broaden the experiment only to resolve an identified problem. No further
optimization is scheduled merely to confirm an expected reduction.

Run cargo fmt --all, cargo test --workspace, and cargo clippy --workspace
--all-targets -- -D warnings at this completed-item gate. Run the graph-core
and graph-IR property targets explicitly with --features proptest --test
property, and the affected graph property/conformance suites, including
umol-graph's --features conformance --test resolution when its consumers are
covered by this migration. Activate umol-py/.venv in every shell that builds/tests
Python or includes it in a workspace build; verify Python 3.13 before the first
build. Rebuild with maturin develop before pytest -q umol-py/tests. Reconcile
documentation and completion markers after the gates pass.

Critical path: S0a -> S1 -> S2 -> S3 -> S4 -> S5, with S0b required before
S2d and its dependent surface/consumer changes. The stage order preserves the
requested separation: primitives first, known duplicate work next, then two
explicit cost decisions. No implementation choice is justified solely by a
microbenchmark improvement or by uniformity of the interface.
