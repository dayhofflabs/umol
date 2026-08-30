# 214 — Aggregate frame semantics

Status: In Progress
Date: 2026-08-28
Relates: [103](103-stereochemistry-overlay-and-ports-2026-05-28.md),
[104](104-stereochemistry-implementation-plan-2026-05-31.md),
[111](111-stereo-phase-c-impl-2026-06-11.md),
[134](134-reaction-application-overlays-2026-06-26.md),
[137](137-python-bindings-2026-07-05.md),
[140](140-python-entity-bindings-plan-2026-07-09.md),
[168](168-api-hygiene-2026-07-27.md),
[186](186-molecule-canonicalization-2026-08-05.md),
[208](208-canonicalization-scaling-2026-08-24.md),
[209](209-normalization-canonical-semantics-2026-08-25.md),
[211](211-relation-frames-and-api-2026-08-26.md),
[212](212-remapping-layer-2026-08-26.md),
[213](213-editor-overlay-storage-2026-08-27.md),
[215](215-integrity-minimization-2026-08-28.md),
[data-type guide](../docs/development/data-types.md),
[nomenclature guide](../docs/development/nomenclature.md)

## Scope

This document replaces the unfinished scope of doc
[209](209-normalization-canonical-semantics-2026-08-25.md). It settles the stereo-frame domain that
blocked that plan, records every aggregate publication boundary that must enforce it, removes the
orbit semantics and machinery that repeated virtual ligands required, and carries the remaining
molecule, reaction, reaction-span, canonicalization, and complete-only API work into one staged
implementation plan.

It also supersedes, in part, the linked historical stereo documents wherever they permit the same
`StereoLigand` value to occur more than once in one stored frame or derive an occurrence orbit from
such repetition. Their prochirality, topicity, relation-storage, parser, and canonicalization
discussions remain historical design records and are not replaced wholesale.

Canonical-search scaling remains in doc
[208](208-canonicalization-scaling-2026-08-24.md), the general remapping representation remains in
doc [212](212-remapping-layer-2026-08-26.md), and reaction-application result provenance remains in
doc 204 on the atom-mapping branch. This document settles the frame semantics those consumers use;
it does not absorb their independent scope.

## Inherited completed state

The successor begins after work that remains implemented:

- doc 209 S0 added private effective-level inspection and routed complete aggregate
  canonicalization, equality, and hash operations through it;
- doc 209 S1 added stereo-kind frame operations and stereo-site incidence integrity;
- doc 211 made relation storage preserve supplied participant sequences, placed entry-level frame
  selection and payload transport in graph IR, and supplied in-place participant permutation; and
- doc 209 S2a recorded that relation-surface dependency as settled and introduced no implementation
  of its own.

The public description-level retirement proposed by doc 209 has not happened. It remains in the
plan below rather than being reported as inherited work.

## Current evidence

Review at clean commit `387650962` and the subsequent entry-point audit established:

- molecule integrity deliberately checks uniqueness only for actual-atom stereo ligands, so equal
  virtual ligands currently pass;
- kindless stereo frames may currently exceed `umol-perm`'s private `MAX_DEGREE`, although every
  actual stereo frame operation eventually needs a bounded `Permutation`;
- `Reaction::new` is unchecked, its fields are public, reaction `Add` deltas do not receive the
  molecule stereo-entry checks, and the live Python reaction facade can publish a fresh unchecked
  snapshot;
- molecule-level `Constraint::StereoAtom` and `Constraint::StereoBond` wrappers check their kind's
  degree and payload but not whether that kind is admissible for the referenced atom or bond site;
- `Permutation::visit_between`, `Permutation::between_all`, `FrameAction`, `find_reframed`,
  `virtual_block_swaps`, and the product-candidate portions of molecule equality and canonicalization
  exist specifically because equal frame values make the frame action non-unique;
- `reaction::application::test_reaction_apply_stereo_atom_update` and its bond twin expose a
  separate, older frame transport error: a rule `old` is compared to the concrete host with `equiv`
  rather than `matches`, and the actual host value is not installed as the realized transaction
  pre-state;
- reaction application transports configuration modifications and removals but skips stereo
  `ModifyConstraint` deltas and molecule-level `ConstraintDelta`s that carry frame-relative values;
- reaction application has no corresponding rule-to-host transport for aromatic or multicenter
  electron-count deltas, and its stereo-only helper therefore cannot serve as the complete overlay
  frame boundary;
- reaction integrity compares ordinary overlay-removal participants as unordered incidence, does
  not compare a created entity's `Remove` frame with its `Add` frame, and delta folding or span
  materialization may consequently combine a removal payload with the lhs or added frame without
  first establishing that they are the same frame;
- `RelationalConstraint::NoncovalentBondEndsSatisfy` assigns its two predicates to the stored
  endpoint positions, so sorting or aligning a noncovalent frame without moving that molecule-level
  value changes the assertion;
- stereo-bond integrity gives every frame two consecutive endpoint blocks independently of whether
  its form asserts `CisTrans`, `Axial`, or no kind, so form-selected full symmetric sorting can
  produce an incidence-invalid bond frame; and
- constitution-level canonicalization still fails to select one presentation for aromatic and
  multicenter participant frames.

The baseline was re-run at `d1f63468d`, whose only change after `387650962` was documentation. It
reproduced doc 209's complete fifteen-test semantic failure ledger:

- the graph-IR library passed 6361 tests, failed four, and ignored four;
- the canonicalization integration target passed 14 and failed one; and
- the property target passed 321, failed all ten persisted failure cases, and ignored one in
  708.61 seconds.

The four library failures are `test_canonicalize_constitution_entity_kind_minimum`,
`test_canonicalize_constitution_participant_order`,
`test_aromatic_systems_glue_differing_frames`, and `test_molecule_meet_pushout_overlays`. The
integration failure is `reaction_span::test_reaction_span_canonicalize::case_2_constitution`. The
property failures are `molecule::comparison::test_molecule_equiv_under_reframed`, both
`reaction::application::test_reaction_apply_stereo_*_update` properties, the four reaction
canonicalization equality/hash/by/reversal properties, and the three corresponding reaction-span
hash/canonicalize/reversal properties.

The same audit found one unrelated mechanical failure: the `React` rustdoc example had a malformed
`use` statement introduced by `387650962`. It was repaired with this planning correction and is not
part of the semantic ledger. Formatting and the `umol-graph-core`, `umol-perm`, `umol-graph`, and
`umol-io` test suites were green; Python was not rerun for this baseline.

## Settled semantic contract

### Stored stereo frames are bounded and pairwise-distinct

Every stereo atom and stereo bond frame published as part of a `Molecule`, `Reaction`, or
`ReactionSpan` satisfies both rules:

1. `frame.len() <= umol_perm::MAX_DEGREE`.
2. Every complete `StereoLigand { atom_id, kind }` value occurs at most once in the frame.

The second rule is deliberately equality of the complete ligand value, not merely of `atom_id` or
`StereoLigandKind`. Its consequences are:

- actual-atom ligands remain distinct by atom id;
- one implicit hydrogen and one lone pair borne by the same atom are distinct and may coexist;
- two implicit hydrogens borne by the same atom are invalid, as are two lone pairs borne by the
  same atom;
- the same virtual-ligand kind on opposite stereo-bond endpoints remains valid because the bearing
  atom ids differ; and
- explicit hydrogen atoms are ordinary atom ligands with distinct ids and therefore remain able to
  represent positions that must be distinguished.

These are tier-1 representation-integrity rules. They do not depend on a chemistry model, on whether
the site is stereogenic, or on whether the represented constraints are satisfiable.

The fixed maximum is a property of the current permutation representation, not a claim that no
coordination geometry can ever have a larger degree. Supporting larger frames later requires first
replacing or extending the permutation carrier; accepting an inoperable frame now is not forward
compatibility.

### Prochirality remains expressible without repeated virtual values

Prochirality, topicity, ligand symmetry, and fluxionality remain derived or assertable properties of
an integrity-valid stereo entity. The restriction removes only a representation in which two
indistinguishable virtual occurrences were treated as separate frame positions.

For implicit hydrogens, a use case that needs distinct positions may represent the hydrogens as
explicit atoms. Repeated lone-pair positions have no equivalent atom representation and are excluded;
that loss of an extreme edge case is accepted in exchange for one-valued frame transport throughout
the model.

Because legal frames are pairwise-distinct, there is no residual occurrence stabilizer. A stored
coset is a value in one ordered frame, not a representative of an orbit generated by exchanges of
equal frame entries. Ligand-symmetry, fluxionality, and topicity constraints keep their ordinary
kind, degree, and position-domain rules; they acquire no additional invariance or orbit-closure rule.

### `#T`, `#C`, TableIR, and other external sources

`#T` and `#C` assert stereo properties that perception or resolution must realize using an
integrity-valid frame. They do not create distinguishable occurrences of an equal virtual ligand.
If a site could be represented only by repeating an implicit hydrogen or lone pair, no stereo entity
is perceived from that completion; the operation's existing stereo failure policy governs that
absence. It is a semantic perception or resolution outcome, not permission to publish a malformed
frame.

The current TableIR tetrahedral and cis/trans raise paths already introduce at most one virtual
position in each relevant local ordering. Their final `Molecule::try_from_entries` boundary remains
authoritative. An explicit repeated or oversized frame from the graph-IR DSL, a future TableIR field,
Python, or another external source must fail at raise rather than be normalized, deduplicated, or
assigned an arbitrary coset.

### Integrity and semantic failure remain separate

The molecule error surface gains exact cases for:

```text
DuplicateStereoLigand { entity, ligand }
StereoFrameDegreeTooLarge { entity, degree, maximum }
```

The maximum check runs before kind-dependent arity or permutation work, so an oversized kindless
frame receives the bounded-frame diagnostic rather than reaching a panic or a later contradiction.
Duplicate checking is deterministic and reports the first repeated value in frame order.

Reaction integrity applies the same local checker to every added stereo entry and wraps or mirrors
the same entity-specific diagnosis; it does not implement a second version of the rule. A reaction's
lhs remains covered by molecule integrity. Reaction-span integrity remains projection-based, so both
sides pass through molecule integrity and receive side-specific errors.

The `StereoKind` carried by a molecule-level `Constraint::StereoAtom` or
`Constraint::StereoBond` is a site assertion. Its admissibility is therefore checked with the same
atom-site or bond-site table as an inline configuration. It need not equal an independently stored
configuration kind: disagreement between two individually admissible assertions is semantic
constraint mismatch or contradiction, not malformed representation.

An inability to perceive `#T` or `#C`, a frame-to-frame action outside a configured kind's group, a
failed pattern match, and an unsatisfiable combination of valid constraints remain semantic failures.
They must not be reported as molecule integrity errors.

### A reaction has one owning participant frame per overlay entity

A `Reaction` uses one owning participant coordinate frame for every overlay entity id it carries.
For an entity present on the lhs, that frame is the lhs entity's frame. For a created entity, it is
the frame in the unique `Add` delta. Every field or constraint delta for that id is stated in the
owning frame even though those delta variants do not repeat the participants.

Every overlay `Remove` repeats the source entity's site and structured participant incidence, but it
may state its recorded attributes in another ordering of those participants. The repeated sequence
is an explicit local frame. Because integrity prohibits repeated participants, equal structured
incidence determines one local removal-to-owner action. Removal matching transports the payload by
that action before applying its value relation; construction preserves the supplied frame and does
not normalize it eagerly.

The structured incidence comparison respects every entity kind's factors. Dative acceptors and
donors remain distinguished; a stereo-bond action permutes within endpoint blocks and may swap the
two complete blocks. Moving a ligand between blocks without that complete swap is
`ReactionIntegrityError::IncidenceMismatch`, as are a different site or participant multiset.
`ParticipantFrameMismatch` is removed because a compatible ordering is transportable and an
incompatible ordering is different structured incidence.

The single per-id owning action remains the aggregate witness. Let `q` be a removal's derived
local-to-owner action and `a` the owner-to-target action. `FrameTransport for Reaction` transports
that removal by the conjugate `q.compose(a).compose(q.inverse())`, keeping the same local-to-owner
relation while the owner moves under `a`. This is what gives raw locally framed removals the same
identity, inverse, and composition laws as the rest of `FrameTransport`. Normalization is a
different operation: it applies `q` to align the removal directly with its owner. Because `Reframe`
normalizes first, `q` is then identity and the conjugate reduces to `a`. Field and constraint deltas
consume the owner action directly. This does not require the deltas to materialize a reaction span
and does not change the older delta-normalization rule that a created entity subsequently removed is
an unobservable no-op. Changing an overlay's participants is still expressed by removing the old
entity and adding a new entity with a new id; no delta changes an existing entity's incidence.

### The three nested operations

Normalization, reframing, and canonicalization are the three prefixes of one semantic pipeline:

```text
normalize     = normal form
reframe       = normal form + participant frame
canonicalize  = normal form + participant frame + entity ids
```

Each outer operation includes every inner operation. `Reframe` therefore has `Normalize` as a
supertrait, and `Canonicalize` has `Reframe` as a supertrait once the aggregate implementations land.
Every canonicalized result is consequently also reframed and normalized. For integrity-valid inputs
and a fixed canonicalization context, their equality relations nest in the same order:

```text
==  implies  normalized_eq  implies  framed_eq  implies  canonical_eq
```

`normalized_eq` is a provided method on `Normalize`: it compares successful normal forms by
structural equality, counts two intrinsic contradictions as equal, and short-circuits structural
equality. The separate blanket `Equiv` trait and its ambiguous `equiv` method retire; there is no new
`NormalizedEq` trait. `framed_eq` is the corresponding provided method on `Reframe` and compares
the selected reframed values under the same contradiction rule.

The selected public shape separates receiving an action from owning enough frame information to
derive one:

```rust
pub trait FrameTransport: Sized {
    type Action;

    fn reframe_by(self, action: &Self::Action) -> Option<Self>;
}

pub trait Reframe: Normalize + FrameTransport {
    fn representative_action(&self) -> Self::Action;

    fn reframe_with_action(self) -> Result<(Self, Self::Action), Contradiction>;

    fn reframe(self) -> Result<Self, Contradiction>;

    fn framed_eq(&self, other: &Self) -> bool {
        self == other
            || self.clone().reframe() == other.clone().reframe()
    }
}
```

`FrameTransport` does not normalize and does not select a frame. It applies an independently
supplied compatible action, returning `None` when the action lacks required entity coverage or has
an incompatible domain, degree, or admissible subgroup observable by that receiver. A payload-only
receiver checks every condition encoded by its own value; the frame-owning aggregate additionally
checks the action against the actual participant sequence. `Reframe` is implemented only by a
carrier that owns participant frames.

For every integrity-valid carrier, `representative_action` is total and the action it returns
applies successfully to that same carrier. The action is derived from the input's structural frame
owners before semantic normalization, so an intrinsic contradiction cannot make action derivation
fallible. `reframe_with_action` obtains that input-domain action, normalizes the value, applies the
same action, and performs the final reduction required after position-bearing keys move. It returns
the action together with the result; action entries for values erased by normalization may remain
and are ignored. `reframe` produces the same reframed value but is an independently implemented
fused operation: it may derive, apply, and discard each local action as it visits the owning entry
and must not construct a complete aggregate action solely to discard it.

Both transformations consume their receiver, as `normalize` and `canonicalize` do. This lets an
owned copy-on-write aggregate reuse uniquely owned storage. The trait makes no exact allocation or
clone-count promise: sorting and transport may need temporary storage, and a downstream consumer may
require a selective internal lookup. The complete associated action is nevertheless a requested
result, not an obligatory intermediate representation. `framed_eq` may likewise be overridden by
an aggregate to compare through a fused path rather than materializing complete actions. Comparison
borrows and may clone in order to construct the representatives. Frame actions obey identity,
inverse, and composition, and applying a composite action is componentwise. The action returned by
`reframe_with_action` has the direction from the input's stored frame to the returned representative.

This prefix order is semantic, not a requirement that the canonicalization implementation mutate
storage in that literal order. Canonical search determines an entity-id witness using normalized,
frame-aware keys. Because the final participant presentation is expressed in the selected id space,
construction applies that id remapping and then runs aggregate `reframe` there. The returned value
nevertheless satisfies all three prefixes.

### Per-entity frame actions, with no occurrence search

A frame action always has the direction

```text
selected[i] = stored[action[i]]
```

For an ordinary relation this is a possibly unbounded `DynPermutation`: a checked
permutation of integer participant indices with identity, inverse, composition, and `between`
operations. `umol-perm` owns this storage-neutral action on a sequence; graph IR owns which
positions may move and how payload follows. Stereo uses the existing bounded `Permutation`, whose
`act` operation has the same direction. Pairwise-distinct legal participants make the action between
a stored and representative presentation unique. The graph-core storage-facing
`ParticipantPosition` type is not part of this action carrier.

Atoms and localized bonds do not bear graph-IR payload frames. A localized bond's endpoints are
unordered topology, and no current bond value is indexed by endpoint position. Reframing begins at
the six overlay kinds:

| Overlay kind | Representative frame | Action group | Values transported under the action |
| --- | --- | --- | --- |
| Dative bond | donors sorted by `AtomId`; singleton acceptor fixed | `S_n` on the donor factor, represented by `DynPermutation` | donor sequence only; current form and molecule-level dative constraints are frame-invariant |
| Aromatic system | participating atoms sorted by `AtomId` | `S_n`, represented by `DynPermutation` | participant electron counts; charge, unpaired electrons, and current constraints are invariant |
| Multicenter bond | participating atoms sorted by `AtomId` | `S_n`, represented by `DynPermutation` | participant electron counts; charge, unpaired electrons, and current constraints are invariant |
| Noncovalent bond | endpoints sorted by `AtomId` | `S_2`, represented by `DynPermutation` | endpoint sequence and every `NoncovalentBondEndsSatisfy` predicate pair referring to the entity; the entity form is invariant |
| Stereo atom | complete `StereoLigand` values sorted | `S_n` on the ligand factor, represented by bounded `Permutation` | configuration, ligand-symmetry and fluxionality permutations, topicity positions, inline constraints, and molecule-level `StereoAtom` leaves |
| Stereo bond | sort within each two-ligand endpoint block, then order the two blocks lexicographically | wreath product `S_2 wr S_2 = (S_2 x S_2) semidirect S_2`, represented by bounded `Permutation` | the corresponding stereo-bond configuration, inline constraints, and molecule-level `StereoBond` leaves |

Stereo-bond selection is structural rather than form-selected. Its storage is a birelation between
a singleton site factor and a ligand factor, but the latter has an additional two-block partition
not captured by the storage type. A legal action may reorder either endpoint block and may exchange
the two complete blocks; it may not move one ligand across the endpoint boundary. The resulting
block-preserving group is the wreath product above, also the `CisTrans` and `Axial` parent group.
Current atom-centered kinds all use the full symmetric parent group, so complete sorting is
admissible for every integrity-valid stereo atom assertion. Independently asserted admissible kinds
may be semantically inconsistent, but they do not make frame selection ambiguous.

Representative selection therefore belongs to the overlay aggregate, which knows the site and
frame-bearing factors. A form implements `FrameTransport` and can apply a local action, but it does
not implement `Reframe`. The existing stereo-only `FrameAction` trait and form-level `select_frame`
retire in favor of this common operation. `DynPermutation::between` and
`Permutation::between` derive the unique action between two legal presentations; representative
selection is the corresponding one-frame derivation. `find_reframed`, `Permutation::visit_between`,
and `Permutation::between_all` are removed after their callers migrate. `virtual_block_swaps`
and repeat-specific symmetry generators, fixtures, and property strategies are removed or rewritten
as integrity rejection cases.

### Frame actions are complete reusable witnesses

A correspondence is the witness for an entity-id action, and a frame action is the witness for a
participant-frame action. Both have identity, inverse, and composition laws. A frame action is not
merely assembly bookkeeping: the same action selected from an entity frame transports its form,
span sides, field and constraint delta values, and every frame-relative constraint that refers to
it. A removal stated in a compatible local frame consumes the conjugate of that owning action by
its local-to-owner alignment, so transport preserves the relation between the two frames.

`FrameTransport::Action` denotes the complete action on `Self`, not one row of a separately returned
vector. The carrier grows with the value being acted on:

| Carrier | Complete associated action |
| --- | --- |
| One ordinary entity form or `EntitySpan<form>` | One `DynPermutation` |
| One stereo entity form or `EntitySpan<form>` | One bounded `Permutation` admissible for that entity kind |
| One kind-specific overlay delta payload | One local action for that entity kind |
| One entity-kind aggregate or its `*Spans` peer | One local action per dense typed entity id |
| `Molecule` or `ReactionSpan` | One `OverlaysFrameAction` value covering every entry in all six overlay kinds |
| `Reaction` | One `OverlaysFrameAction` value covering every lhs and `Add`-owned overlay entity |
| `ConstraintDelta`, `Constraint`, or `ConstraintSpan` | The compatible `OverlaysFrameAction` value whose typed maps cover every frame-relative entity reference carried by the receiver |
| `Deltas` in isolation | No complete action: removals may use local frames, and only the containing `Reaction` supplies their owning frames |

`OverlaysFrameAction` is a typed six-component composite. Each field maps the corresponding entity
id to its local `DynPermutation` or bounded `Permutation`; atom and localized-bond entries are
absent because those entity kinds have no participant-frame action. Molecule and reaction-span
producers fill dense id domains, while a reaction adds sparse created-entity domains to the dense
domains of its lhs. A consumer requires every action needed by the frame-relative values it carries
and ignores entries for other entities. An independently supplied missing or incompatible action
makes `reframe_by` return `None`, not panic.

A complete action's domain is its typed id set together with the degree of every local action.
Identity is relative to that domain, inverse preserves it, and composition is defined only for equal
typed id sets with equal per-entry degrees; the operations then act componentwise. This exact-domain
algebra is distinct from consumer compatibility, which requires only that the action cover the
receiver's needed subset and permits irrelevant entries. Although operation-issued, a frame action
is not bound to one source object's identity: it may be reused with another structurally compatible
carrier, and the contextual consumer checks the coverage and compatibility its promised result
requires.

The carrier is operation-issued: `representative_action`, `reframe_with_action`, and frame-alignment
operations produce it with private fields establishing its internal invariants. No independently
assembled public constructor is required. Its presence in the public associated type is justified
by actual cross-carrier use rather than symmetry with `MoleculeCorrespondence`. This is a Rust
graph-IR action carrier; no Python action-map class is added in this scope because no supported
Python workflow independently consumes the witness.

For an entity-kind aggregate, the dense action sequence is indexed by its typed entity id and is
shared with its span-aggregate peer. `reframe_with_action` returns this complete witness. Ordinary
`reframe` remains the happy path and does not construct the complete witness.

### Aggregate normalization, reframing, and canonicalization

Entity forms do not implement `Reframe`: by definition they exclude their participants and cannot
know their current frame or select another one. They implement `FrameTransport` with the local
entity-kind action. The ordinary `reframe_to(from, to)` combination reduces to
`DynPermutation::between(from, to)` followed by `reframe_by`; stereo uses the bounded
counterpart. Compatibility is receiver-relative: a fixed-arity form or a determinate positional
payload checks the degree it exposes, while a variable-arity, dimensionless, or frame-invariant form
cannot reconstruct an absent owning frame and accepts every structurally valid action in its local
action group. The overlay aggregate always checks the actual frame degree and admissible action
group. Frame-invariant forms return themselves under every such compatible local action.

`EntitySpan<T>` is likewise payload-only rather than frame-owning. It has the generic
`FrameTransport` implementation whenever `T` does, applying one action to every present side through
`try_map`; it does not implement `Reframe`. This covers every kind-specific entity span without
introducing span-form newtypes. `Constraint`, `ConstraintSpan`, and kind-specific overlay delta
payloads also implement `FrameTransport` for a supplied local action, but none derives an action
independently. `Deltas` does not implement transport under owning `OverlaysFrameAction` in isolation:
only `Reaction` has the owning frames needed to conjugate each per-id action into a compatible
removal-local frame.

The six overlay aggregates and their span counterparts are the first frame-owning carriers.
They were introduced because each overlay entity kind, unlike a graph-core relation set, knows which
factor bears the frame, which is a site, how its payload moves, and which typed entity id keys the
action. Each aggregate implements reduction-only `Normalize` without changing relation ids or stored
participant sequences, implements `representative_action` as the dense per-id aggregate action, and
applies that complete action through `FrameTransport`. Its plain `reframe` instead derives and
immediately applies each local action without collecting the dense aggregate action. A span aggregate
uses the same aggregate action type as its molecule peer and applies each local action to every
carried side.

`Normalize for Molecule` reduces every carried form and constraint without changing ids or
participant frames. `representative_action` assembles the six dense aggregate actions into one
`OverlaysFrameAction`. `FrameTransport for Molecule` applies that complete action to all six overlay
aggregates and to the complete constraint tree; atom and localized-bond vectors are unchanged apart
from normalization performed by `Reframe`. The final reduction restores normalized constraint
ordering and deduplication after position-bearing keys change. `reframe_with_action` returns the same
complete action that was applied. Plain `reframe` does not call `reframe_with_action` on its overlay
aggregates and then apply and discard their complete actions. It derives each local action once and
uses it immediately for the owning entry and any frame-relative molecule constraint. Each
kind-specific constraint form exhaustively declares whether it uses participant positions and
defines its transport, so adding a constraint variant requires an explicit classification rather
than silently treating it as invariant. A single constraint-tree traversal builds a sparse action
domain over all six overlay kinds. Empty domains and action maps allocate no backing storage; for a
kind with no frame-relative consumer, plain `reframe` takes the ordinary fused aggregate path. For a
nonempty domain, it retains only the requested local actions and applies them to the constraint tree
after all overlays have been reframed.

`Molecule::framed_eq_under(other, correspondence)` is the explicit entity-id-witness comparison. It
is `try_remap(correspondence)` followed by `framed_eq`, returns `false` when the open correspondence
does not describe a complete dense source remapping, and remains an inherent molecule operation
rather than adding an id-witness type to `Reframe`. The old `equiv_under` name retires. Under the
identity correspondence it reduces to `framed_eq`. For integrity-valid inputs whose canonicalization
succeeds, `canonical_eq` holds exactly when some admissible total correspondence makes
`framed_eq_under` hold. Two intrinsic contradictions still compare equal under the equality
failure-totalization rule, but that convention does not promise a correspondence witness.

Mapped equality and pushout derive the one frame action for each matched entry from the supplied
entity correspondence and transport the complete assertion before comparing or meeting it.
Identical stored frame sequences imply identity; there is no hidden stabilizer action to enumerate.

Remapping relabels ids and preserves participant sequence. For the correspondence `c` returned by
molecule canonicalization, the required composition is:

```text
canonicalize_with_correspondence(x).0 == canonicalize(x)
reframe(x.remap(c)) == canonicalize(x)
```

Frame selection therefore occurs after id transport. The private complete-stereo candidate-product
search and the kindless position-order fallback are removed once the aggregate path transports the
complete constraint store. Graph automorphism and canonical-label search remain; they solve id-frame
selection, not repeated-ligand occurrence selection.

### Reaction and reaction-span transport

`Normalize for Reaction` reduces the lhs and deltas in their stored id and participant frames
without selecting another frame. Normalized `Deltas` fold entity chains and have a deterministic
semantic order; reframing does not promise to preserve incidental source order that normalization
has already erased. `Reframe for Reaction` obtains an existing entity's source frame from the raw
lhs and a created entity's source frame from its unique raw `Add` delta. Its complete input-domain
action covers every lhs overlay id and every overlay id introduced by an `Add`, including an entry
whose delta chain normalization later erases. `reframe_with_action` derives that total action before
normalization, applies it to the normalized lhs and deltas, and ignores entries no longer required by
the normalized receiver. Plain `reframe` derives and consumes the same actions incrementally and
retains only any sparse internal lookup needed when several deltas refer to the same entity.
`FrameTransport for Reaction` applies an independently supplied composite to the lhs and dispatches
the corresponding local action to every delta. Field and constraint deltas consume the owning
action directly. For a `Remove`, the reaction derives the local removal-to-owner alignment and
conjugates the owner action into the removal's coordinates before transporting:

- dative `Add` and `Remove` donor sequences;
- aromatic and multicenter `Add` and `Remove` forms and both sides of electron-count field changes;
- noncovalent `Add` and `Remove` endpoint sequences and ordered predicate pairs in
  `ConstraintDelta`;
- every stereo `Add`, `Remove`, configuration change, constraint change, and stereo leaf in
  `ConstraintDelta`.

`ConstraintDelta` delegates to `FrameTransport for Constraint`, which recursively transports only
the position-bearing leaves and accepts irrelevant entries in the composite action.
`FrameTransport for Reaction` returns the transported representation without reducing it;
`reframe` and `reframe_with_action` normalize after transport. `Delta::remap` relabels ids while
preserving the supplied participant sequence; it does not silently select a frame.

`Normalize for ReactionSpan` reduces every carried side and constraint span without selecting a
participant frame. Its six existing `*Spans` aggregates derive the same aggregate actions as their
molecule peers and assemble them into `OverlaysFrameAction`; atom and localized-bond spans remain
plain vectors because they have no frame action. `FrameTransport for ReactionSpan` applies the
composite to every span aggregate and every `ConstraintSpan` without reducing the result. Its
reframing operations perform the final reduction. Plain `reframe` derives each local action once and
transports every present side and affected constraint without materializing that complete composite.
A modified entity therefore cannot let its lhs and rhs choose different presentations.

Raw span construction preserves an explicitly supplied `Modified { lhs, rhs }` even when the two
values are `normalized_eq`; the tag has no additional semantic assertion beyond those side values.
`Normalize for ReactionSpan` collapses that no-op to `Unchanged`, so `Reframe` and `Canonicalize`
inherit the same collapse through their normalization prefix. `ReactionSpan::superimpose` is a
deriving producer rather than a faithful entry conversion: it may emit the standardized
`Unchanged(lhs)` directly for a `normalized_eq` pair, retaining exact lhs projection and semantic
rhs projection under the induced correspondence.

Reaction application performs the alignment counterpart of aggregate selection. For every matched
existing overlay it derives the unique action from the rule frame after atom mapping to the host's
stored frame and transports all values listed in the entity table. Added entities have no host
frame and retain the internally aligned participant sequence and payload they carry. The reframed
rule `old` is a pattern and is tested with `matches`; once it matches, lowering replaces it with the
concrete host value as the transaction's realized `old`. The same pass transports aromatic and
multicenter electron values, noncovalent ordered predicate pairs, stereo configuration changes,
inline `ModifyConstraint` values, and every frame-relative molecule-level `ConstraintDelta`. A
missing or inadmissible stereo action remains `StereoFrameMismatch`; it does not create alternative
products. No result or match-witness extension is needed.

## Publication and entry-point inventory

An **ingress** accepts independently assembled frame data. A **mutation boundary** changes a
published aggregate in place. A **publisher** returns an immutable aggregate after edits or
transformations. All three must uphold the same contract, but trusted transformations should prove
preservation rather than scatter duplicate validation logic.

| Surface | Current gate | Required disposition |
| --- | --- | --- |
| `StereoAtoms::new`, `StereoBonds::new`, `StereoAtomSpans::new`, `StereoBondSpans::new`, and raw relation-set `From` conversions | Public low-level construction can bypass aggregate checks. | Restrict raw assembly to graph IR; public consumers construct frames through checked/asserted aggregates, and no public unchecked route feeds frame operations. |
| `Molecule::try_from_entries` / `Molecule::from_entries` | Checked / asserted pair through `Molecule::check_integrity`. | Add the two frame rules and wrapper site-kind rule to the one authoritative check. |
| `Molecule::{stereo_atom_mut, stereo_bond_mut, modify_stereo_atoms, modify_stereo_bonds, constraints_mut}` and public overlay `attributes_mut` | Mutable references and unchecked replacement closures can install a kind, position domain, or reference that disagrees with the containing molecule after construction. | Restrict overlay mutation to graph IR and replace integrity-sensitive molecule mutation with checked replacement or transactional methods that validate before committing. No public mutation may leave a published `Molecule` outside its tier-1 contract. |
| `MoleculeBuilder::build`; `MoleculeEditor::{snapshot, try_build, build}`; `Molecule::apply` | Publish through editor integrity, with asserted and checked variants. | Retain the shared gate; add exact regressions for transient invalid frames rejected only at publication. |
| `Molecule::{try_remap, remap, extract, combine_all, combine, split}` and fragment, glue, pushout, canonicalization, resolution, and reaction-product paths | Trusted transforms, usually rebuilding through checked or asserted molecule construction. | Audit every publisher and add preservation properties; do not add redundant front-door checks to each operation. |
| `Reaction::new`, struct literals, `Reaction::check_integrity`, and `Reaction::from_sides` | `new` and public fields are unchecked; callers must remember a later check, and removal incidence is not established against the lhs or added entity. | Add `try_new`, make `new` the asserted sibling, make fields private with read/parts accessors, validate every added stereo frame and stereo constraint wrapper, and require each removal to repeat the source's structured incidence while accepting any uniquely alignable local ordering. |
| Reaction composition, derivation reversal, span conversion, canonicalization, and application | Trusted reaction publishers currently use unchecked `Reaction::new` or literals. | Route through the checked/asserted constructor and add preservation properties. |
| `ReactionSpan::{try_from_entries, from_entries, superimpose}` and reaction/span conversions | Checked / asserted entry pair; integrity projects both sides through `Molecule::try_from_entries`; construction currently normalizes an equivalent `Modified`. | Retain projection as authority and verify added, removed, modified, and unchanged stereo spans reject repeated or oversized frames on the affected side. Remove eager normalization from entry construction, preserve the raw tag there, and let `superimpose` emit `Unchanged(lhs)` for a `normalized_eq` pair. |
| Parsed molecule, reaction, and reaction-span DSL; `mol_dsl!` and typed `IntoIr` paths | Parsed molecule/span paths map aggregate integrity; reaction paths can construct unchecked IR. | Map checked failures to `ParseError::InvalidValue`; keep macros and typed conversion as asserted paths over trusted input. |
| TableIR raise and SMILES/CTfile ingestion | TableIR ends at `Molecule::try_from_entries`; current TableIR has no raw-frame field and its stereo raise emits no repeated frame. | Preserve the central gate, prove current `#T`/`#C` paths never publish an invalid frame, and document repeated-virtual completion as perception absence rather than deduplication. |
| Python `Molecule.from_entries`, parse/SMILES, and editor snapshot/build/apply | Molecule entry and editor publication are checked. | Preserve `ValueError` mapping and add exact duplicate/maximum cases. |
| Python `Reaction(...)`, parse, `from_sides`, reaction-SMILES, live component snapshots, and all operations consuming a live reaction | The facade can assemble or later expose unchecked Rust `Reaction` values. | Validate initial construction and every live snapshot consumed by an operation; map `ReactionIntegrityError` to `ValueError` without making transient Python edits impossible. |
| Python `ReactionSpan.from_entries`, parse, and conversions | Entry construction is checked; typed wrapping assumes valid Rust input. | Preserve the checked boundary and add both-side error cases. |

Leaf `StereoAtomForm`, `StereoBondForm`, and constraint constructors do not own a ligand frame and
cannot enforce frame uniqueness or length. They retain their own kind, coset, permutation-degree,
and position checks; the first aggregate that pairs them with a frame owns the frame invariant.

## Staged implementation plan

Tests and benchmarks begin with the first implementation stage. Exact fixtures use nonidentity
actions and nonuniform position-sensitive payloads; uniform values cannot prove that transport
happened. Property cases cover all six overlay kinds, every entity/span state, renumbering,
stereo-bond endpoint blocks, noncovalent ordered predicate pairs, constraints, and the
checked/asserted boundary.

The plan inherits fifteen semantic failures from the relation-storage migration. The repaired
rustdoc example is green before implementation begins and is not part of that ledger. Because a
stage boundary must be green, integrity closure, orbit removal, and the dependent aggregate
transport work form one S0 recovery stage. Its subitems may reduce the inherited ledger early but
must never add an unlisted failure. The expected latest closure points are:

| Checkpoint | Remaining failures | Closed at that checkpoint |
| --- | ---: | --- |
| S0 start through S0n | 15 | No expansion of the inherited ledger. |
| S0o | 12 | Aromatic glue, molecule pushout, and mapped molecule equality. |
| S0p | 10 | Constitution entity-kind minimum and participant-order canonicalization. |
| S0t | 8 | Two inherited stereo-application regressions, closed by generalized overlay transport. |
| S0u | 0 | Reaction-span integration plus seven reaction/span canonicalization properties. |
| S0v | 0 | Coherent frame-transport, reframing, and full-pipeline algebraic property suite. |

An earlier closure is welcome and is recorded at implementation time; a later closure or a new
failure is not. Each subitem runs its focused checks against this ledger, and S0 ends with the full
graph-IR normal and feature-gated property suites green.

The slow canonicalization integration and feature-gated property targets are not per-subitem gates.
Run them only at a checkpoint that owns a ledger closure, or at the S0 stage boundary; other
subitems use dependency-local tests plus the relevant compile, documentation, format, and lint
checks.

### S0 — Integrity, orbit removal, and frame-transport recovery

#### Integrity foundation and closed publication boundaries

- **S0a — current baseline and public degree bound** (`umol-perm/src/permutation.rs`, graph-IR
  focused suites, existing benchmarks): retain the exact pre-change test ledger above, record the
  pre-change benchmark baseline, and make the existing `MAX_DEGREE` a documented `pub const` rather
  than duplicating `6` in graph IR. Cover the constant through the checked permutation constructors.
  **Additive; inherited red ledger unchanged.** [dep: none] **Done.** `MAX_DEGREE` is public and
  documented as a representation limit. A downstream integration test imports it and checks that
  `Permutation::try_from` accepts degree `MAX_DEGREE` and returns the exact `ImageTooLong` error for
  `MAX_DEGREE + 1`.

  The first benchmark attempt exposed a stale `para_stereo_cascade` fixture that predated stereo-site
  incidence and kind-admissibility integrity. Its intended covalent incidence was restored and its
  atom-inadmissible `CisTrans` labels were replaced by distinct tetrahedral configurations. The full
  command
  `cargo bench -p umol-graph-ir --bench canonicalize -- --save-baseline doc214-s0a` then passed and
  saved the S0a baseline. Across the corpus, incidence construction spans approximately
  0.18-2.46 µs, remapping 2.00-13.48 µs, topology and constitution canonicalization 18-224 µs,
  structure and para-stereo structure canonicalization 71 µs-10.4 ms, and full canonicalization
  102 µs-15.7 ms. `cargo test -p umol-perm --all-features` passed 181 unit tests, the dedicated
  integration test, and 26 property tests; strict clippy passed for all `umol-perm` targets and the
  graph-IR canonicalization benchmark. The inherited fifteen-test semantic ledger remains the S0
  starting point.
- **S0b — molecule stereo-frame integrity and mutation gates**
  (`umol-graph-ir/src/ir/{molecule,molecule/integrity,view/stereo}.rs`, overlay modules): add
  `DuplicateStereoLigand` and `StereoFrameDegreeTooLarge`; check maximum then pairwise distinctness
  for stereo atoms and bonds before kind-dependent operations; make molecule-level stereo wrapper
  kinds use the existing atom/bond site-admissibility table. Restrict overlay `attributes_mut` to
  graph IR and replace public stereo-form and reference-bearing constraint mutation that can
  compromise integrity with checked replacement or transactional operations that validate before
  committing. Add exact atom, bond, kindless, opposite-endpoint, mixed-virtual, mutation rollback,
  and error-precedence tables plus generator laws that every published molecule satisfies the two
  frame invariants. **Breaking; accepted-domain, mutation-surface, and public-error changes, with the
  inherited red ledger unchanged after caller and fixture migration.** [dep: S0a] **Done.**
  `Molecule::check_integrity` now rejects frames longer than `MAX_DEGREE` before duplicate, arity,
  kind, or permutation checks, rejects the first repeated complete `StereoLigand`, and applies the
  atom/bond site-admissibility table to molecule-level stereo constraint wrappers. Reference
  validation remains earlier than both frame checks. Exact constructor tables cover atom and bond
  repeats of each virtual kind, mixed virtual kinds, the same virtual kind at opposite bond
  endpoints, kindless oversized frames, wrapper kinds, and diagnostic precedence.

  The raw stereo-form, overlay-attribute, and molecule-constraint mutation surfaces are now internal
  to graph IR. Public callers use five checked, rollback-safe molecule operations for one or all
  stereo atoms, one or all stereo bonds, and the molecule constraint tree; each operation has a
  dedicated success, failure, and rollback regression. The Python live views and setters use those
  same gates and preserve their published exception mapping. Graph, Python, reaction, transaction,
  incidence, canonicalization, and pushout fixtures were migrated to valid distinct frames; tests
  whose sole contract was the abandoned repeated-ligand orbit semantics were removed. The molecule
  generators now sample ligand frames without replacement, and a dedicated property checks both
  invariants on every generated stereo atom and bond.

  `cargo test -p umol-graph-ir --lib --no-fail-fast` ran 6,367 tests and retained exactly the four
  inherited library failures (`test_aromatic_systems_glue_differing_frames`, the two constitution
  canonicalization regressions, and `test_molecule_meet_pushout_overlays`). The focused molecule
  structure and edit property modules passed 9 and 17 tests. The full canonicalization integration
  target retained its one inherited failure, and the full feature-gated property target passed 322,
  failed exactly its ten inherited cases, and ignored one in 117.52 seconds. Graph-IR doctests, all
  974 graph library tests, 1,635 Python-binding Rust tests, and 1,324 Python tests passed. Workspace
  all-target compilation and strict clippy for graph IR, graph, and Python also passed.
- **S0c — reaction integrity and construction** (`umol-graph-ir/src/ir/reaction.rs`,
  `ir/reaction/integrity.rs`): add `Reaction::try_new`, make `Reaction::new` asserted, privatize
  `lhs`/`deltas` behind direct read and `into_parts` accessors, and apply the shared local stereo-entry
  checks to every `Add` plus site-kind checks to stereo constraint deltas. Record the source frame of
  every overlay id from the lhs or its unique `Add`; require every `Remove` to repeat the source site
  and exact participant sequence; and add `ParticipantFrameMismatch` for equal incidence in a
  different order while retaining `IncidenceMismatch` for a different site or participant multiset.
  Check invalid references, incidence, and frame equality in that order. Migrate struct literals and
  trusted producers in the same subitem. Exact cases cover invalid lhs, invalid additions, wrapper
  kinds on existing and added entities, removals for all six existing and created overlay kinds,
  diagnostic precedence, and valid explicit-H prochiral inputs. **Breaking; inherited red ledger
  unchanged after caller migration.** [dep: S0b] **Done.** `Reaction` now has private `lhs` and
  `deltas` fields, direct `lhs` / `deltas` borrows and `into_parts`, a checked `try_new`, and an
  asserted `new`. Reaction integrity reuses the molecule stereo-entry and site-kind checks for
  stereo additions, inline constraint changes, and molecule-level stereo constraint wrappers;
  these failures are reported as
  `ReactionIntegrityError::StereoIntegrityError(MoleculeIntegrityError)`. Exact addition tests cover
  repeated and oversized atom and bond frames, including maximum-degree rejection, while a valid
  explicit-H prochiral case confirms that distinct hydrogen atom ids remain admissible.

  Integrity records one private source frame for each lhs overlay and each unique raw overlay
  `Add`. Every overlay removal first checks references, then unordered incidence including stereo
  site and complete ligand values, then exact participant order. Equal incidence in another order
  is `ParticipantFrameMismatch`; changed incidence remains `IncidenceMismatch`. The exact table
  covers all six overlay kinds for both lhs-owned and `Add`-owned frames, and a separate table
  fixes reference/incidence/frame diagnostic precedence. Struct literals and trusted consumers were
  migrated to the accessors and checked/asserted construction pair; obsolete properties that
  depended on publicly constructing a malformed `Reaction` now exercise `try_new` directly.

  `cargo test -p umol-graph-ir --lib --no-fail-fast` ran 6,398 tests and retained exactly the four
  inherited library failures, with 6,390 passing and four ignored. The full feature-gated property
  target ran 332 tests and retained exactly the ten inherited failures, with 321 passing and one
  ignored. Graph-IR doctests passed, all 974 graph library tests, all 3,294 IO library tests, and
  1,635 Python-binding Rust tests passed with two ignored. The nightly graph-IR fuzz build, the
  Python-3.13 workspace all-target check, and strict clippy for graph IR, graph, IO, and Python all
  passed.

  **Superseded in part by doc 215.** The exact-order restriction and
  `ParticipantFrameMismatch` diagnostic implemented here are not part of the final integrity
  contract. Compatible removal-local frames are accepted and transported through their unique
  alignment with the owning lhs or `Add` frame; a non-compatible structured frame remains
  `IncidenceMismatch`. Doc 215 also removed `ReactionIntegrityError::Lhs` and the invalid-lhs case
  after molecule closure made that failure unreachable, then removed the defensive post-publication
  reaction checks and diagnostics that only those checks could produce. The retained S0c contract
  is private reaction storage, checked/asserted construction, local stereo validation, reference
  validation, and structured-incidence validation. Doc 215 completed this unwinding before S0e.
- **S0d — reaction-span and low-level carrier audit** (`umol-graph-ir/src/ir/stereo.rs`, overlay
  modules, `ir/reaction_span.rs`): restrict raw frame-bearing collection construction, conversions,
  and mutation to their actual graph-IR assembly role, and verify that checked/asserted span
  construction reports the lhs or rhs projection carrying an invalid frame. Retain no public
  unchecked frame publisher. **Breaking; inherited red ledger unchanged after carrier migration.**
  [dep: S0b] **Done.** Raw `new` / `into_entries` collection methods for all six overlay kinds
  and their reaction-span peers are now graph-IR-private, and the public raw relation-set `From`
  conversions were removed. The molecule editor's trusted publication path uses explicit
  crate-private `from_arc` assembly; public consumers receive read-only collections from checked
  aggregates, with only the intrinsically valid empty `Default` available independently.

  `ReactionSpanEntries` is documented as an open carrier, while `try_from_entries` and
  `from_entries` establish both projected molecules through the authoritative molecule integrity
  gate. A dedicated eight-case checked-constructor table covers atom and bond stereo spans in every
  added, removed, modified, and unchanged state, with repeated and oversized frames reported as the
  affected `Lhs` or `Rhs`. A separate asserted-constructor table verifies that the same side-specific
  error is retained in the panic message.

  The focused reaction-span run passed all 159 selected tests. The graph-IR library ran 6,408 tests,
  passing 6,400 and ignoring four while retaining exactly the four inherited failures. Graph-IR
  doctests, the Python-3.13 workspace all-target check, formatting, strict graph-IR clippy, and the
  public-surface audit passed. Per the checkpoint verification rule, the slow canonicalization
  integration and property targets remain deferred to their owning closure checkpoint.

At S0d, all independent Rust aggregate constructors are locally green and malformed repeated or
oversized frames cannot enter an immutable public aggregate through those constructors. The
inherited transport ledger remains.

The post-S0d admission audit found that the broader closed-container contract was not yet minimal or
fully closed. Doc [215](215-integrity-minimization-2026-08-28.md) completed the dative and
entity-count corrections, aromatic/multicenter mutation closure, defensive-check removal, and
reaction-span construction and derived-standardization semantics before S0e. S0e-S2 use that
corrected integrity domain and error surface.

#### External publication boundaries

- **S0e — graph-IR DSL raise** (`umol-graph-ir/src/dsl/{molecule,reaction,reaction_span,stereo}.rs`,
  DSL macros): route resolved reactions through `try_new`, map frame failures to
  `ParseError::InvalidValue`, and retain typed `IntoIr` and macros as asserted paths. Roundtrip valid
  frames and reject repeated/oversized atom, bond, reaction-add, and span frames plus reaction
  removal structured-incidence mismatches with exact parser diagnostics; accept compatible
  reordered removal frames without normalizing them during construction. **Breaking; parsed reaction
  failure propagation without expanding the inherited red ledger.** [dep: S0c, S0d] **Done.**
  Parsed reactions now publish through `Reaction::try_new` and translate its exact diagnostic to
  `ParseError::InvalidValue`; typed `IntoIr` and macro paths retain asserted publication. Focused
  molecule, reaction-add, and reaction-span tables cover repeated virtual ligands and frames above
  `MAX_DEGREE` for both stereo entity kinds, including side-specific span diagnostics. Existing
  valid-frame roundtrips remain green, and a typed reaction roundtrip proves that a compatible
  reordered aromatic removal keeps its recorded local frame. DSL removal syntax denotes its source
  entity and therefore reconstructs the source incidence; incompatible removal incidence is not a
  representable parsed surface value and remains covered by the checked Rust reaction constructor.
  The 4 molecule cases, 46 reaction-input cases, 10 span-error cases, and 8 typed reaction
  roundtrips passed.
- **S0f — TableIR and format raise** (`umol-io/src/table_ir/raise.rs`, SMILES and CTfile ingestion):
  retain `Molecule::try_from_entries` as the final gate, test that the current `#T`/`#C` paths never
  synthesize repeated frames, and add a semantic case where unavailable repeated-virtual completion
  follows stereo failure policy. Record that any future explicit TableIR frame field must continue
  through the same checked raise boundary. **Additive; inherited red ledger unchanged.** [dep: S0b]
  **Done.** TableIR raise still assembles `MoleculeEntries` and publishes only through
  `Molecule::try_from_entries`. Module documentation now records that a future explicit TableIR
  participant-frame field must remain behind that boundary. Focused SMILES and CTfile cases show
  that tetrahedral `#T` and cis/trans `#C` are currently raised as atom/bond constraints while both
  stereo-entity collections remain empty, so these paths cannot synthesize a repeated frame. The
  semantic `[C@H2](F)Cl` case, whose requested tetrahedral completion would require two virtual
  positions, returns the existing `TetrahedralLigandCount { count: 2 }` raise error rather than
  publishing or deduplicating a frame. All 18 focused TableIR conversion cases and both focused
  CTfile-to-IR stereo cases passed.
- **S0g — Python publication** (`umol-py/src/{molecule,reaction,reaction_span,transaction}.rs` and
  Python tests): preserve molecule/span `ValueError` mapping; validate `Reaction(...)`, parse,
  reaction-SMILES, and every operation snapshotting live lhs/deltas through `Reaction::try_new`.
  Cover construction and mutation-to-invalid-state, including invalid removal incidence and valid
  compatible removal reordering, without changing the live editing model. Build and test only under
  the repository Python 3.13 environment. **Breaking; snapshot methods become fallible where
  expanding the inherited red ledger.** [dep: S0c, S0d, S0e, S0f] **Done.** Python `Reaction`
  construction now publishes detached lhs and delta snapshots through `Reaction::try_new`, as do
  the post-lowering parse, side-comparison, and reaction-SMILES paths. The internal live-component
  snapshot is fallible and authoritative: rendering, equality, span conversion, reversal,
  composition, application, fingerprinting, canonicalization, and `Molecule.react` / `react_all`
  all propagate its `ValueError`. Component getters and setters retain the live editing model and
  deliberately permit a temporarily incoherent pair; the next interpreting operation reports the
  exact integrity failure.

  Dedicated Rust binding and public Python regressions reject an incidence-mismatched removal at
  construction, accept and retain a compatible reordered aromatic removal without eager
  reframing, and show that mutation into the same invalid state succeeds until `render` snapshots
  the reaction. A separate parse regression retains the public `ParseError` taxonomy for an
  invalid repeated-virtual stereo addition. The Python reaction Rust module passed 131 tests with
  its two inherited ignores, and the public reaction test module passed 100 tests with two skips,
  all under repository Python 3.13.
- **S0h — trusted publisher preservation** (`umol-graph-ir` remapping, editing, extraction,
  combination, pushout/glue, reaction/span conversion, composition and application; `umol-graph`
  resolution/transform modules): inventory every publisher in the table above, route construction
  through the checked/asserted pairs, and add preservation properties rather than repeating local
  checks. **Breaking; migrate every stale producer without expanding the inherited red ledger.**
  [dep: S0c, S0d] **Done.** The audit grouped every listed publisher by its terminal aggregate
  gate. Entry construction, raw carrier restriction, and checked mutation remain owned by S0b-S0d;
  DSL, format, and Python ingress remain owned by S0e-S0g. Molecule editors, extraction, splitting,
  fragments, reaction application, resolution, and chemistry transforms publish through
  `try_build` / `build`; combination and canonicalization publish through
  `try_from_entries` / `from_entries`; pushout/glue uses `try_from_entries`. Reaction derivation,
  composition, reversal, span conversion, and canonicalization use `Reaction::try_new` / `new`;
  span superimposition, conversion, remapping, and canonicalization use
  `ReactionSpan::try_from_entries` / `from_entries`.

  Dense molecule remapping was the sole stale path: it assembled a public molecule through the
  private raw-arc constructor. Arc assembly is now a checked `try_from_arcs` boundary shared by
  entry construction and editor publication, and trusted remapping asserts that the pure
  renumbering preserves integrity. No operation acquired a defensive input check.

  Five molecule properties republish remapping, extraction, combination, splitting, and successful
  edit application exactly through the checked editor gate. Reaction application and span
  properties republish product molecules, derived reactions, and both span projections; both
  composition properties republish every composite through `Reaction::try_new`. Four chemistry
  properties cover aromatization, charge delocalization, both kekulization algorithms, and full
  resolution. The seven focused graph-IR integrity properties, two composition properties, and four
  graph properties passed. The Python-3.13 workspace all-target check and strict all-target clippy
  for graph IR with properties, graph with properties, and Python passed. The slow inherited-ledger
  suites remain deferred to their designated closure checkpoints.

At S0h, Rust, DSL, TableIR/format, and Python boundaries enforce the same rejection, while transient
editors remain permitted and fail only when publishing. The inherited transport ledger remains.

#### Frame-action foundation and repeated-ligand orbit removal

- **S0i — action carriers and local transport** (`umol-perm/src/dynamic.rs`,
  `umol-graph-ir/src/ir/{traits,electrons,aromatic,multicenter,dative,noncovalent,stereo,delta}.rs`):
  add checked unbounded `DynPermutation` with identity, inverse, composition, and `between`;
  add `FrameTransport` with associated `Action` and consuming `reframe_by`; and implement local
  transport for the six entity forms, stereo configuration values, kind-specific overlay deltas,
  and generic `EntitySpan<T>` where `T: FrameTransport`. Keep the old helpers temporarily while
  callers migrate. Make compatibility receiver-relative: each leaf checks every exposed fixed degree,
  positional length, and stereo group condition, while dimensionless variable-arity payloads defer
  actual-frame checks to their owning aggregate. Assert the action laws, action direction,
  incompatible-domain absence, exact acceptance for leaves with and without an observable degree,
  stereo outside-group failure, and one-action transport of every present span side with nonuniform
  values. **Additive; inherited red ledger unchanged.** [dep: S0b] **Done.** `umol-perm` now
  exposes the checked, unbounded `DynPermutation` over integer indices with the stated direction and
  group operations. Graph IR exposes the consuming `FrameTransport` trait and implements local
  transport for all six overlay forms, stereo configurations and frame-relative inline constraints,
  all six kind-specific overlay deltas, and generic `EntitySpan<T>`. Compatibility is checked only
  where the receiver exposes it: positional values check length, fixed noncovalent and stereo-bond
  carriers check their action groups, kinded stereo values check their parent group even when the
  coset is undetermined, and dimensionless variable-arity values carry unchanged. The temporary
  `FrameAction`, `reframe_to`, and orbit-search helpers remain for S0j migration.

  Eighteen focused dynamic-action unit cases, two action-law properties, and forty-four focused
  graph-IR transport cases pass. Strict all-target Clippy with permutation and graph-IR property
  features passes. The inherited red ledger was not rerun or changed.
- **S0j — unique alignment at graph-IR consumers** (`ir/stereo.rs`, `ir/molecule/editor.rs`,
  `ir/molecule.rs`, `ir/substructure.rs`, `ir/view/stereo.rs`, `ir/reaction.rs`,
  `ir/reaction_span.rs`): replace `find_reframed` and candidate enumeration in lookup, glue, editor
  equality, mapped equality, matching, and span selection with one
  `DynPermutation::between` or `Permutation::between` followed by `FrameTransport`.
  Retire ordinary form `reframe_to` after its callers use the same two kernels. Preserve
  operation-specific `meet`, `matches`, or exact comparison semantics; do not collapse them into one
  generic callback. **Breaking; inherited red ledger unchanged after consumer migration.**
  [dep: S0h, S0i] **Done.** Every graph-IR alignment consumer now derives one checked action and
  transports through `FrameTransport`: overlay glue retains `meet`, editor and mapped equality retain
  normalized equivalence, substructure matching retains pattern-directed `matches`, reaction
  application retains its host-value precondition, and span superimposition retains exact side
  selection. Mapped molecule constraints consume the same unique per-entity stereo action rather
  than recursively searching a candidate product. Ordinary `reframe_to` methods and their duplicate
  tests are removed; `find_reframed` is unused but remains exported until the S0l legacy-trait
  removal, while the remaining `between_all` sites belong to S0k frame selection or S0p
  canonicalization.

  The migration exposed and corrected one S0i domain omission: kindless stereo-bond forms and
  positionless stereo-bond constraints now reject permutations outside the endpoint-block
  `S_2 wr S_2` parent just as kinded forms and stereo-bond deltas already did. Exact form,
  constraint, and editor tables cover within-endpoint swaps, complete endpoint-block swaps, and an
  illegal cross-endpoint movement. The obsolete active editor orbit test and ignored mapped-equality
  stabilizer test were removed because their repeated frames cannot be published under S0b.

  Focused graph-IR runs passed 46 editor cases, 10 mapped-equality cases, 45 substructure cases, 18
  local-removal cases, 12 span-superimposition cases, 12 stereo-pushout cases, 8 stereo-reaction
  application cases, 4 stereo-view coset cases, 42 local `reframe_by` cases, and the new five-case
  stereo-bond transport-domain tables. Graph-IR doctests and strict all-target Clippy with the
  property feature passed. The two directly affected inherited failures retain their prior
  frame-direction diffs (`test_aromatic_systems_glue_differing_frames` and
  `test_molecule_meet_pushout_overlays`); the untouched canonicalization failures and slow property
  ledger were not rerun.
- **S0k — frame selection and symmetry cleanup** (`ir/stereo.rs`, `ir/symmetry.rs`, property
  strategies, overlay modules): implement the per-kind representative rules in the six overlay
  aggregates. Stereo atoms sort the complete ligand frame under the full symmetric group; stereo
  bonds sort within their two endpoint blocks and then order the complete blocks under the wreath
  product `S_2 wr S_2`. Add dense complete aggregate-action carriers shared with each `*Spans` peer
  and the typed six-component `OverlaysFrameAction` composite. Make each complete action carry its typed id
  domain and local degrees: identity is domain-relative, inverse preserves the domain, composition
  requires equal domains and degrees, and consumers may accept a covering superset. Remove form-level
  `select_frame`,
  `virtual_block_swaps`, residual-stabilizer generation, orbit-representative minimization, and
  repeat-valid fixtures. Replace them with integrity rejection, atom/bond structural-domain cases,
  distinct-frame symmetry evidence, and complete-action identity/inverse/composition cases.
  **Breaking; inherited red ledger unchanged after fixture migration.** [dep: S0j] **Done.** The
  six overlay aggregates and their `*Spans` peers now return kind-specific complete action carriers
  keyed by their typed dense ids: four carry `DynPermutation`, and the stereo aggregates carry
  bounded `Permutation`. The operation-issued carriers expose typed lookup plus exact-domain
  identity, inverse, and composition; composition rejects unequal id domains or local degrees.
  `NoncovalentBondsFrameAction` admits only degree two and `StereoBondsFrameAction` admits only the
  endpoint-block `S_2 wr S_2` subgroup. The private-field `OverlaysFrameAction` composite provides
  the same algebra componentwise without an independently assembled public constructor.

  Representative selection is structural and owned by each overlay aggregate. Dative donors, aromatic and
  multicenter atoms, and noncovalent endpoints sort completely. Stereo atoms sort the complete
  distinct ligand frame under `S_n`; stereo bonds sort within endpoint blocks and then order the
  complete blocks under `S_2 wr S_2`, including a nonidentity block-exchange case. Form-level
  `select_frame`, span-side arbitration, candidate enumeration, `virtual_block_swaps`, and the
  repeat-valid unit/property generator surface are gone. Symmetry projection now obtains ligand
  permutations only from explicit graph automorphisms; the distinct explicit-ligand
  stereogenic/prochiral tests remain the evidence. The S0b integrity tables remain the rejection
  surface for repeated complete ligands. The temporary stereo `FrameAction` and `find_reframed`
  compatibility surface remains for S0l, while the only graph-IR `between_all` production
  uses left are the private canonicalization paths owned by S0p.

  Eight complete-action carrier cases, 85 focused reframing cases, eight framed-equality cases, 49
  symmetry-filter cases, seven focused stereo properties, and four duplicate-stereo-ligand
  integrity cases pass. The graph-IR all-target property build and strict Clippy pass, and
  `cargo fmt -p umol-graph-ir` completes. The workspace-wide formatter remains blocked by the
  unrelated pre-existing `iter:from_fn` parse error in `umol-py/src/reaction.rs`; no Python file was
  changed in this subitem. The inherited red ledger was not rerun or changed.
- **S0l — quotient and transport trait migration**
  (`umol-graph-ir/src/ir/{stereo,traits,constraint,delta}.rs`, overlay modules and re-exports;
  `umol-py` form bindings; `umol-perm/src/permutation.rs`): add `normalized_eq` directly to
  `Normalize`, implement reduction-only `Normalize` for every overlay and span aggregate,
  aggregate, and make consuming `Reframe: Normalize + FrameTransport` expose
  `representative_action`, `reframe_with_action`, `reframe`, and provided `framed_eq`. Its associated
  action is the complete `FrameTransport::Action`, not one row of a separately returned vector.
  Define `representative_action` on the input's structural frame owners before semantic
  normalization, keep it total for every integrity-valid value, and make `reframe_with_action`
  return that same input-domain action after normalization, transport, and final reduction.
  Keep `reframe` independently implementable rather than defining it through
  `reframe_with_action`; aggregate implementations fuse local derivation and transport without
  materializing a complete witness, and may override `framed_eq` for the same reason.
  Implement composite `FrameTransport` for `Constraint`, `Constraints`, `ConstraintSpan`,
  `ConstraintDelta`, and `Delta`; each requires actions for the frame-relative entity ids it carries
  and ignores irrelevant map entries. Do not implement it for standalone `Deltas`: different
  removal payloads may carry different local frames, and only the owning `Reaction` has the
  alignment needed to assemble one complete action. Remove `Equiv`, `Normalized`, the old
  stereo-only `FrameAction`,
  `find_reframed`, and `visit_between`. Migrate every Rust and Python caller according to whether it
  asks normalized equality, framed equality, frame transport, or pattern matching; retain the
  distinct inherent `Molecule::equiv` only until S0m and rename the Python form method to
  `normalized_eq`. Retain
  `Permutation::between_all` temporarily only for the private canonicalization path that S0p
  removes, and retain `Permutation::between` as the final alignment primitive. Build the Python
  caller migration under the repository Python 3.13 environment. **Breaking; inherited red ledger
  unchanged after all callers migrate.** [dep: S0j, S0k] **Done.** `Normalize` now owns
  `normalized_eq`; `Equiv` and `Normalized` are gone, and every form-level Rust and Python caller
  uses the new name. The six frame-owning overlay aggregates and all six span peers implement
  reduction-only `Normalize`, complete `FrameTransport`, and consuming `Reframe`; plain `reframe`
  fuses local action derivation and transport, while `reframe_with_action` returns the same dense
  typed witness used to reproduce the result. Recursive constraints, constraint stores, constraint
  spans and deltas, and individual `Delta` values transport through `OverlaysFrameAction` without
  assigning a false complete-action meaning to standalone `Deltas`.

  The stereo-only compatibility trait and search helper and `Permutation::visit_between` are
  removed. `Permutation::between_all` remains only in the private canonicalization path named
  above. Focused reframing tests (95), normalization tests (263), and permutation tests (48) pass;
  the graph-IR all-target property build and strict graph-IR/umol-perm Clippy pass. The Python
  caller compiles and its focused `normalized_eq` test passes under the repository Python 3.13
  environment. The inherited red ledger was not rerun or changed.

At S0l, no general graph-IR operation or public helper interprets an equal-value occurrence
permutation as observable stereo state. The isolated private canonicalization candidate product
remains only until S0p replaces it with complete aggregate transport. The inherited transport
ledger remains.

#### Complete molecule frame transport

- **S0m — reduction-only molecule normalization** (`umol-graph-ir/src/ir/molecule.rs` and overlay
  modules): implement `Normalize for Molecule` over every entity and constraint store without
  changing ids or participant frames, and remove the inherent `Molecule::equiv` in favor of the
  provided `normalized_eq`. Cover empty/nonempty overlay aggregates, shared storage, idempotence, contradiction
  propagation, and agreement of `normalized_eq` with structural equality on normalized values.
  **Breaking; inherited red ledger unchanged after molecule callers migrate.** [dep: S0l]
  **Done.** `Molecule` now implements consuming reduction over both inherent form vectors, all six
  overlay stores, and the molecule constraint store while leaving graph topology, ids, and stored
  participant frames unchanged. The inherent `equiv` method is removed and its callers use the
  provided `normalized_eq`. Exact unit coverage includes empty and populated aggregates, shared
  copy-on-write storage, idempotence, and contradictions from both inherent and overlay payloads;
  the three focused molecule equality properties agree with the settled laws. Focused unit and
  reaction-span caller tests and strict graph-IR Clippy pass. The comparison property target retains
  only its listed `test_molecule_equiv_under_reframed` failure, which remains assigned to S0o.
- **S0n — aggregate molecule reframing** (`ir/molecule.rs`, overlay modules): implement
  `FrameTransport` and consuming `Reframe for Molecule` from the normalized molecule. Assemble each
  dense aggregate action into one `OverlaysFrameAction`, apply that same complete action to all six
  overlay aggregates and the recursive constraint tree, and return it only on the
  `reframe_with_action` path. Implement plain `reframe` as a fused pass that derives each local
  action once and immediately transports the owning entry and affected constraints; use only a
  selective internal lookup when needed to avoid rescanning the recursive constraint tree.
  Reduce the transported molecule again so changed constraint keys sort and deduplicate. Prepare the
  aromatic/multicenter participant-order repair and benchmark `reframe`,
  `representative_action`, `reframe_with_action`, and `framed_eq` separately from this subitem
  onward. Include empty and few-overlay molecules, large aromatic systems, frame-relative
  constraints, and owned and shared copy-on-write inputs. Assert the per-entity-kind action table,
  complete-action identity/inverse/composition, action reuse against compatible carriers,
  missing/incompatible action absence,
  stereo-bond block preservation, nonuniform endpoint predicates, exact idempotence,
  `normalized_eq => framed_eq`, and no mutation of shared inputs. **Additive; inherited red ledger
  unchanged.** [dep: S0m]
  **Done.** `Molecule` now implements complete six-component `FrameTransport` and consuming
  `Reframe`. The witness-returning path assembles the operation-issued `OverlaysFrameAction`; plain
  `reframe` instead derives each local action once, transports the owning entry immediately, and
  retains only actions named by a sparse frame-relative constraint domain over all six overlay
  kinds. Each kind-specific constraint form exhaustively classifies its variants and implements
  transport; current frame-invariant variants require no action entry, while noncovalent endpoint
  predicates and position-bearing stereo constraints do. Both paths reduce after transport. Exact
  coverage exercises the classification and no-action paths, all six action components,
  nonuniform positional payloads and recursive constraints, compatible witness reuse,
  missing/incompatible witnesses, the empty representative, contradiction propagation,
  idempotence, and shared copy-on-write input preservation. Four molecule properties cover the
  complete action algebra, fused/witnessed agreement, replay, representative identity, pipeline
  implication, and framed equality. Separate Criterion groups cover `representative_action`,
  `reframe_with_action`, `reframe`, and `framed_eq` on empty, few-overlay, and 128-participant
  aromatic inputs. Focused reframing tests, four new properties, the benchmark build, and strict
  graph-IR Clippy pass. The comparison property target retains exactly its inherited
  `test_molecule_equiv_under_reframed` failure for S0o; its other six properties pass.
- **S0o — mapped equality and pushout** (`ir/molecule.rs`, glue/pushout consumers): replace the two
  hand-written `equiv_under` paths with inherent `framed_eq_under`, defined as checked remapping by
  the supplied correspondence followed by `framed_eq`; retire the old name. Use the same unique
  frame alignment and complete constraint transport in pushout, including noncovalent ordered
  endpoint predicates. Assert reduction to `framed_eq` under identity, inverse correspondence,
  renumbering, nonidentity ordinary and stereo actions, stereo-bond block swaps, contextual
  correspondence rejection, and nonuniform constraint cases. State and test the exact
  `canonical_eq` correspondence-witness characterization only where canonicalization succeeds;
  intrinsic-contradiction equality remains totalized without requiring a witness. Close
  `test_aromatic_systems_glue_differing_frames`, `test_molecule_meet_pushout_overlays`, and
  `molecule::comparison::test_molecule_equiv_under_reframed`. **Breaking; inherited red ledger
  decreases from fifteen to twelve.** [dep: S0n]
  **Done.** `Molecule::framed_eq_under` now performs checked entity-id remapping followed by
  `framed_eq`; the hand-written mapped-comparison paths and the `equiv_under` name are removed.
  Molecule pushout preserves the selected entity frame and transports every frame-relative
  constraint before entity-id remapping, including ordered noncovalent endpoint predicates and
  stereo atom and stereo bond constraints. The frame-invariant constraint path constructs no
  aggregate action map. Exact unit tests cover identity, renumbering, ordinary and stereo frame
  differences, stereo-bond endpoint-block swaps, partial and incompatible correspondences, and
  nonuniform pushout constraints. Properties cover identity reduction to `framed_eq`, inverse and
  composed correspondences, and participant-frame changes. The canonical-equality tests restrict
  the correspondence-witness characterization to successful canonicalization and separately
  preserve totalized equality for intrinsic contradictions. The three assigned regressions pass;
  the renamed comparison property is
  `molecule::comparison::test_molecule_framed_eq_under_participant_frame`. The all-core,
  feature-gated graph-IR checkpoint reports 6,878 passed, 12 inherited failures, and 4 skipped;
  those failures are exactly the two S0p canonicalization cases and the ten later reaction and
  reaction-span cases. Graph-IR doctests and strict all-target Clippy pass.
- **S0p — canonicalization composition** (`ir/canonicalize.rs` and focused benchmarks): after
  canonical id remapping, call the aggregate frame operation; remove complete-stereo candidate
  products and the above-`MAX_DEGREE` position-order fallback, then remove
  `Permutation::between_all` and its repeat-specific unit/property surface. Preserve
  graph-automorphism search, exact hash/equality behavior, bounded-exhaustive minimum checks,
  renumbering invariance, and the corrected post-remap law. Compare the unified path with the
  pre-change distinct-frame baseline. Close `test_canonicalize_constitution_entity_kind_minimum` and
  `test_canonicalize_constitution_participant_order`. **Breaking; inherited red ledger decreases
  from twelve to ten.** [dep: S0n, S0o]
  **Done.** Every molecule description-level path now applies the selected entity-id correspondence
  and then the aggregate `Reframe` operation. The complete-stereo candidate product, molecule-local
  frame helpers and duplicate normalizer, kindless stereo-atom position-order fallback, and
  `Permutation::between_all` repeat machinery are removed; reaction-span position-order transport
  remains for S0s. Unit and property coverage preserves exhaustive minima, renumbering invariance,
  canonical equality and hashing, and now asserts the exact witness law
  `canonical == source.remap(correspondence).reframe()`. The existing canonicalization operation
  benchmark compared against the S0a `doc214-s0a` baseline: affected structure cases improve by
  about 14–18%, para-stereo structure cases by about 13–21%, and full overlay/stereo cases by about
  22–29%; feature-free full cases improve by about 3–4%. Topology and constitution are mostly
  within noise, apart from small regressions and an 8% topology overlay regression. Both assigned
  unit regressions pass. Graph-IR's 6,505 unit cases, graph-IR and permutation strict all-target
  Clippy, graph-IR doctests, the focused canonicalization properties, and the permutation suite
  pass. The all-core feature-gated checkpoint reports 6,869 passed, 10 inherited failures, and 4
  skipped; the remaining failures are exactly the reaction and reaction-span cases assigned below.

At S0p, one complete molecule frame-transport path is shared by equality, pushout, and
canonicalization. Ten inherited reaction and reaction-span failures remain.

#### Reaction and reaction-span consumers

- **S0q — frame-preserving delta id transport** (`umol-graph-ir/src/ir/delta.rs`, coordinated with
  doc 212): make `Delta::remap` relabel entity and participant ids without sorting the participant
  frame or changing frame-relative payloads. Assert remap/inverse-remap and remap-then-reframe laws.
  **Additive; inherited red ledger unchanged.** [dep: S0l]
  **Done.** `delta::remap_delta` now performs only typed id substitution: dative donors, aromatic
  and multicenter atoms and electron counts, and noncovalent endpoints retain their supplied
  participant frame. The obsolete sorting helper and all four remap-time sorting paths are gone.
  Exact cases cover the former sorting paths, while separate roundtrip and frame-transport tables
  cover all six overlay kinds plus a frame-relative molecule constraint. Graph-IR's 6,519 unit
  cases, the focused reaction-reversal property, and strict all-target Clippy with the property
  feature pass. No inherited failure is assigned to S0q; the S0p ten-failure checkpoint remains the
  ledger baseline.
- **S0r — reaction reduction and reframing** (`ir/reaction.rs`): implement reduction-only
  `Normalize for Reaction` over the lhs and deltas in their stored frames, then `FrameTransport` and
  consuming `Reframe for Reaction`. Before normalization, derive one input-domain
  `OverlaysFrameAction`: existing entity frames come from every lhs overlay and created entity
  frames come from every unique raw `Add`, whose ids need not be dense. Retain entries for created
  entities whose delta chains normalization later erases. Apply it to the normalized lhs and use a
  reaction-owned contextual pass over the normalized deltas, carrying
  dative structural deltas, aromatic and multicenter electron values, noncovalent structural deltas
  and ordered predicate constraint deltas, and every stereo delta and stereo constraint delta under
  the owning entity's one action. Implement plain `reframe` as a fused pass, retaining only the
  sparse internal lookup needed for repeated references to one entity rather than the complete
  public witness. Return the composite from `reframe_with_action` and normalize the transported
  result; do not promise preservation of incidental pre-normalization delta order.
  Cover the trait and action laws, every delta arm and relation-valued entity kind, the complete lhs
  plus sparse created-id domain, compatible removal-local frames and owner actions conjugated into
  those local frames, created entities erased by normalization, intrinsically contradictory but
  integrity-valid reactions with total representative actions, multiple changes to one entity, and
  missing/incompatible map entries. **Additive; inherited red ledger unchanged.** [dep: S0n, S0q]
  **Done.** `Reaction` now implements reduction-only normalization, complete input-domain
  representative actions, contextual `FrameTransport`, and fused `Reframe`. Raw removals conjugate
  the owner action into their local coordinates, while normalization first aligns them with their
  owner; exact identity and noncommuting composition cases enforce the distinction. Focused cases
  cover all six removal kinds, complete lhs plus sparse erased-`Add` domains, repeated changes,
  frame-relative molecule constraints, both stereo constraint kinds, contradiction, and missing or
  incompatible actions. The 256-case reaction reframe properties cover the action laws, witnessed
  and fused agreement, idempotence, and framed equality over comprehensive reactions. A dedicated
  six-kind domain generates compatible removals in a distinct local frame and checks exact action
  identity, inverse/composition, and convergence with the owner-framed representation, including
  noncommuting stereo-bond actions. Graph IR's 6,541 unit cases (three ignored), doctests, and strict
  all-target Clippy with the property feature pass. No inherited failure is assigned to S0r; the S0p
  ten-failure checkpoint remains the ledger baseline.
- **S0s — reaction-span reduction and reframing** (`ir/reaction_span.rs`): implement reduction-only
  `Normalize for ReactionSpan`, then `FrameTransport` and consuming `Reframe for ReactionSpan`.
  Preserve equivalent `Modified` entries at checked and asserted construction, collapse them only
  in normalization and its reframing/canonicalization prefixes, and retain `superimpose` as a
  standardized producer that emits `Unchanged(lhs)` for `normalized_eq` paired values.
  Reuse the existing six `*Spans` packages. On `reframe_with_action`, assemble their aggregate actions
  into `OverlaysFrameAction` and apply it to every entity side and `ConstraintSpan`; on plain
  `reframe`, derive each local action once and immediately transport all present sides and affected
  constraints without constructing the complete composite. Atom and localized-bond spans remain
  plain vectors. Return the composite from `reframe_with_action`, then reduce the complete span.
  Assert the trait and action laws, projection agreement, stereo-bond block
  preservation, nonuniform endpoint predicates, equivalent-`Modified` constructor preservation and
  pipeline collapse, exact lhs and semantic rhs superimposition projections, roundtrips through
  reactions, and modified spans with nonuniform sides. **Additive; inherited red ledger unchanged.**
  [dep: S0n, S0q]
  **Done.** `ReactionSpan` now implements reduction-only `Normalize`, complete
  `OverlaysFrameAction` transport, and fused `Reframe` over its six typed span aggregates. Checked
  and asserted construction still preserve raw equivalent `Modified` entries; normalization and
  the reframe prefix collapse them, normalize every entity side and constraint span, and
  sort/deduplicate the constraint-span set without rebuilding or rechecking the closed container.
  Raw transport leaves atom and localized-bond spans untouched, transports both sides of every
  overlay span under one local action, and transports frame-relative constraint spans through the
  same complete witness. The fused path derives each local action once, records only actions named
  by the sparse constraint domain, whose empty case allocates nothing, and avoids constructing the
  complete witness.
  The reusable internal constraint action map now lives beside that domain and serves molecule
  reframe, molecule pushout, and reaction-span reframe. Canonicalization uses the public span
  normalization implementation instead of its former duplicate. Exact tests cover raw constructor
  preservation and pipeline collapse across all eight entity kinds, constraint-span reduction,
  nonuniform endpoint predicates, incomplete actions, and a nonuniform stereo-bond `Modified`
  entry under a block-swapping representative. The 256-case reaction-span properties cover the
  action laws, fused/witness agreement, the normalization prefix, framed equality, idempotence,
  both projections, and the reaction roundtrip. Graph IR's 6,544 unit cases (three ignored),
  doctests, the focused molecule pushout/reframe regressions, and strict all-target Clippy with the
  property feature pass. No inherited failure is assigned to S0s; the S0p ten-failure checkpoint
  remains the ledger baseline.
- **S0t — rule-to-host overlay-frame transport** (`ir/reaction.rs`, application fixtures): replace
  the stereo-only helper with one application-owned alignment pass over all matched overlay kinds.
  Derive the unique action from each atom-mapped rule owner frame to the host frame, then transport
  the already normalized deltas contextually. Normalization has aligned every compatible removal
  with its owner, so the removal consumes the owner-to-host action directly. If an application path
  instead retains a raw local removal frame, it must align that payload directly by the composition
  of its local-to-owner and owner-to-host actions; this host alignment is intentionally distinct from
  the conjugation used by generic `FrameTransport for Reaction`. This transports
  aromatic and multicenter electron values, noncovalent ordered endpoint predicates, every stereo
  delta arm, and every frame-relative molecule-level constraint delta. Test a reframed rule `old`
  with `matches` and install the actual host value as realized `old`; added entities retain their
  internally aligned frame because they have no host counterpart. Restore both existing stereo
  regressions and add identity/nonidentity ordinary and stereo actions, stereo-bond block swaps,
  nonuniform electron vectors and endpoint predicates, less-specific-old, missing/incompatible
  action, and mismatch cases. Benchmark application separately from match enumeration. **Breaking
  correctness change; inherited red ledger decreases from ten to eight.** [dep: S0j, S0q, S0r]
  **Done.** Reaction application now derives one sparse, typed action map from atom-mapped rule
  owner frames to their matched host frames and transports normalized deltas before lowering. The
  pass covers all six overlay kinds, including aromatic and multicenter electron vectors,
  noncovalent ordered endpoint predicates, every stereo delta arm, inline stereo constraints, and
  molecule-level frame-relative constraint deltas. Direct additions retain their existing frame
  without constructing an action; an identity action is created for an added owner only when a
  frame-relative constraint delta consumes it. Reframed rule `old` values are checked with
  `matches`, then replaced with the concrete host value used by the transaction. Missing ordinary
  correspondences retain `CorrespondenceMismatch`, while missing or inadmissible stereo transport
  retains `StereoFrameMismatch`.

  Frame sensitivity is owned exhaustively by each constraint form and entity delta.
  `ModifyConstraint` delegates to its constraint form; aggregate domain collection and application
  delegate to those per-kind decisions instead of maintaining application-local variant lists.
  Adding a field, inline constraint form, molecule constraint variant, or top-level delta variant
  therefore fails to compile until its frame-domain and transport behavior are specified.

  Exact application tables cover identity and nonidentity ordinary and stereo frames, dative
  donor reordering, aromatic and multicenter nonuniform electron vectors, noncovalent nonuniform
  endpoint predicates, stereo-atom permutation, stereo-bond endpoint-block swap, added-owner frame
  preservation, less-specific old values, and correspondence/frame/match failures. All 24 focused
  `apply_at` cases passed. All 20 reaction-application properties passed at 256 cases, closing the
  two inherited stereo-update failures assigned here; graph IR's 6,553 unit cases passed with three
  ignores, as did doctests and strict all-target Clippy with the property feature. A dedicated
  Criterion target now measures match enumeration and `apply_at` separately; its focused smoke run
  measured approximately 1.69 us and 4.31 us respectively on the implementation machine. The
  inherited red ledger is now eight failures, all assigned to S0u.
- **S0u — aggregate reaction canonicalization** (`ir/canonicalize.rs`): route reaction and
  reaction-span canonicalization through the settled normalization and frame operations, make
  `Canonicalize: Reframe`, and preserve the semantic normalization/reframing/canonicalization prefix
  while constructing the result by selected id transport followed by reframing in the target id
  space. Cover exact canonical equality, hash agreement, reaction/span conversion, and algorithm
  agreement. Close
  `reaction_span::test_reaction_span_canonicalize::case_2_constitution`, the four reaction
  canonicalization properties, and the three reaction-span canonicalization properties.
  **Breaking; inherited red ledger decreases from eight to zero.**
  [dep: S0p, S0r, S0s]
  **Done.** `Canonicalize` now has `Reframe` as its supertrait. Reaction-span search applies the
  aggregate representative action to every id-remapped candidate, and construction reframes the
  selected candidate in its target id space. The previous canonicalization-only stereo
  position-order machinery is removed; ordinary and stereo overlays now share the six-kind
  `ReactionSpan` implementation. Reduced-level comparison keys apply frame transport without
  normalizing excluded data, preserving their documented contradiction independence. The
  correspondence-returning cases now assert the settled witness law: entity-id remapping followed
  by reframing reconstructs the canonical value, while remapping alone need not.

  All 279 canonicalization unit cases and all 15 exact canonicalization fixtures pass. All 19
  molecule, reaction, and reaction-span canonicalization properties pass at 256 cases, including
  the four reaction and three reaction-span failures assigned here. Graph IR's 6,552 unit cases
  pass with three ignores, as do doctests and strict all-target Clippy with the property feature.
  The inherited red ledger is now empty.
- **S0v — frame and quotient algebra property suite**
  (`umol-graph-ir/tests/property/{frame,strategies,molecule/comparison,molecule/canonicalize,
  reaction/canonicalize,reaction/span/canonicalize}.rs`): review the existing stereo-transport,
  molecule-comparison, and molecule/reaction/span canonicalization properties and reorganize them
  into one coherent executable specification without duplicating operation-specific evidence.
  Exercise `FrameTransport` identity, inverse, and composition with independently generated
  compatible actions at the local, entity-kind aggregate, and root aggregate levels. Use domain-relative identity,
  require exact typed domains and degrees for composition, permit covering supersets at consumers,
  and assert exact `None` behavior for missing or incompatible actions. Exercise `Reframe`
  idempotence, representative identity, exact agreement of fused `reframe` with the value returned
  by `reframe_with_action`, and application of the pre-normalization input-domain witness to the
  normalized input followed by final normalization.
  Include entries erased by reaction normalization and assert total `representative_action` on
  integrity-valid intrinsic contradictions. Assert that a reframed value's
  `representative_action` is identity. For `N = normalize`, `R = reframe`, and
  `C = canonicalize`, assert the successful-result fixpoint and absorption laws
  `N(N(x)) = N(x)`, `R(R(x)) = R(x)`, `C(C(x)) = C(x)`, `R(N(x)) = R(x)`,
  `N(R(x)) = R(x)`, `C(N(x)) = C(x)`, `C(R(x)) = C(x)`, `N(C(x)) = C(x)`, and
  `R(C(x)) = C(x)`, together with the comparison implication ladder
  `== => normalized_eq => framed_eq => canonical_eq` and the documented reflexivity, symmetry, and
  transitivity domains of those relations for integrity-valid values under a fixed context. Retain
  the canonical correspondence law
  `reframe(remap(x, c)) == canonicalize(x)` and make every equality relation and precondition
  explicit; require a correspondence witness for `canonical_eq` only on successful canonicalization,
  not for the totalized intrinsic-contradiction class. Cover all six overlay kinds plus
  `Molecule`, `Reaction`, and `ReactionSpan`; separate raw satisfiable, normalized, reframed,
  canonical, intrinsically contradictory, and incompatible-action domains. Use nonidentity
  actions, nonuniform position-sensitive payloads, and independently generated relabelings so a
  production witness is not its own oracle. Preserve minimized regression cases, remove obsolete
  repeat-valid and `equiv`-named assertions, and document why retained overlapping properties have
  different operational domains or validation methods. **Additive; S0 remains green with the
  feature-gated property target passing.** [dep: S0u]
  **Done.** A dedicated frame-property module now checks local forms, all six entity-kind span
  aggregates, and molecule/reaction/reaction-span roots under identity, inverse, and composition.
  It separately covers missing and incompatible action domains, covering supersets, nonidentity
  actions, and position-sensitive payloads. Shared generated scenarios exercise all six overlay
  kinds, independently derived participant actions, entries erased by reaction normalization, and
  integrity-valid intrinsic contradictions.

  The molecule, reaction, and reaction-span suites now assert the complete nine-law
  normalization/reframe/canonicalize matrix, the successful comparison implication ladder and
  relation laws, the totalized contradictory canonical-equality class without a correspondence
  witness, and the correspondence law that remapping followed by reframing reconstructs the
  canonical result. Fused and witness-returning reframing agree, the returned input-domain action
  transports normalized input correctly, and representative actions remain total on intrinsic
  contradictions. Obsolete inverse-only stereo transport and duplicate standalone fixpoint
  properties were removed in favor of the systematic suite. All 370 feature-gated property tests
  pass at 256 cases; strict all-target Clippy with the property feature also passes.

S0 ends green with the pre-orbit frame transport gaps fixed for molecules, reactions, and spans. It
also has one systematic executable specification of the complete quotient pipeline and does not
require a new application-result witness in doc 204.

### S1 — Complete-only public canonicalization surface

- **S1a — Rust API retirement** (`umol-graph-ir`, `umol-graph`, Rust callers and benchmarks): remove
  public description-level selectors and level-parameterized aggregate canonicalization, retain
  private effective-level reduction, and expose only complete public canonicalization/equality/hash
  behavior. Migrate callers in the same subitem. **Breaking (red→green).** [dep: S0v]
  **Done.** `DescriptionLevel` is now private canonicalization machinery, and the duplicate private
  `CanonicalizeLevel` name is eliminated. `Molecule::description_level` and the graph-IR root export
  are removed. `Canonicalize` now exposes only `canonicalize`,
  `canonicalize_with_correspondence`, `canonical_hash`, and `canonical_eq`; its three aggregate
  implementations retain private effective-description-level selection for exact complete
  reduction. The reaction projection helpers and public transformation, equality, and hash paths
  that existed only for caller-selected levels are removed.

  Rust unit and property suites retain complete canonicalization, correspondence, remapping,
  contradiction, equality, hash, and quotient-pipeline laws. Forced-level cases remain only inside
  the canonicalization module as evidence for private effective-level reduction; the public
  description-level unit/property surface and its dedicated property module are removed. The
  benchmark now measures complete canonicalization with and without para-stereo refinement. All
  261 focused canonicalization unit cases and all 20 focused canonicalization properties at 256
  cases pass. All non-Python workspace targets compile, including the canonicalization benchmark,
  and strict all-target Clippy passes for `umol-graph-ir` with the property feature and for
  `umol-graph`. The Python bindings still reference the retired Rust surface and are assigned to
  S1b.
- **S1b — Python API retirement** (`umol-py`): remove Python description-level selection, migrate
  supported methods to the complete operation, and retain exact Rust/Python agreement. Use the
  Python 3.13 build gate. **Breaking (red→green).** [dep: S1a]
- **S1c — living documentation** (`docs/development/{data-types,nomenclature,property-tests,python-api}.md`,
  rustdoc, Python docs): align the public guides with fixed relation repeated-participant integrity,
  including the explicit prohibition on repeated complete stereo-ligand values, bounded frames,
  direct frame transport, `DynPermutation`, the `FrameTransport` / `Reframe` distinction,
  and the rule that plain aggregate `reframe` fuses local action derivation and transport without
  materializing a complete witness. Document the local,
  entity-kind aggregate, and `OverlaysFrameAction` carrier hierarchy, the per-entity-kind representative-action table,
  stereo-bond endpoint-block wreath-product actions, noncovalent ordered predicate transport, and
  delta behavior: id remapping preserves participant frames, reframing transports every
  frame-relative delta payload under the owning entity action, a reaction has one lhs- or
  `Add`-owned participant frame per overlay id, a compatible removal may carry another explicit
  local ordering whose action is composed with the owner action, structured incompatibility is
  `IncidenceMismatch`, raw span construction preserves an equivalent `Modified` tag while
  normalization and standardized producers collapse it, and reaction application installs the
  realized host value as `old` after pattern matching. Also cover input-domain representative
  actions, receiver-relative compatibility,
  exact-domain action algebra with subset-compatible consumers, and the `normalized_eq` /
  `framed_eq` / `framed_eq_under` ladder, operation-issued action and correspondence witnesses,
  the coherent action/reframing/pipeline property laws and their operational domains, `#T`/`#C`
  failure semantics, and the complete-only canonicalization surface. Historical discussion text
  remains in place behind its supersession notices. **Additive (green).**
  [dep: S0g, S0v, S1a, S1b]

S1 is deferrable after S0 if the immediate deliverable is integrity plus correct frame transport;
the core semantics do not depend on retiring public description levels. It remains required before
doc 214 can close as a replacement for all moved doc-209 work.

### S2 — Verification and closeout

- **S2a — cross-surface verification** (workspace): run format, graph-IR unit tests, feature-gated
  property/conformance suites, benchmark builds, clippy, workspace tests, and Python 3.13 extension
  plus pytest. Record exact commands and results in this document. **Additive (green).**
  [dep: S0v, S1c]
- **S2b — evidence review and lifecycle closeout** (docs 168, 186, 208, 209, 211, 212, 214 and
  `discussion/000-status.md`): verify that no living guide or API claims repeat-valid or unbounded
  stereo frames, move any genuinely independent follow-up rather than hiding it, and mark this
  document `Completed` only after every non-deferrable and retained doc-209 item is implemented and
  green. **Additive (green).** [dep: S2a]

## Dependency summary

The critical path is:

```text
S0a-S0h integrity and publication boundaries
  -> S0i-S0l action foundation, orbit removal, and direct transport
  -> S0m-S0p complete molecule transport
  -> S0q-S0u reaction and span transport, restoring the full green tree
  -> S0v coherent quotient-pipeline property suite
  -> S1 public-surface cleanup
  -> S2 verification and closeout
```

S1 is the only deferrable implementation stage. S0 is one deliberately indivisible red-to-green
recovery stage because the relation-storage handoff already left the tree red. Deleting orbit
machinery before closing publication boundaries would turn malformed input into panics, while
implementing aggregate or reaction transport before the deletion would preserve the candidate-set
complexity this decision removes. S0's subitems remain separately reviewable, but none is a stage
boundary and no implementation handoff may treat a still-red checkpoint as complete.

## Directions abandoned rather than moved

The following doc-209 and doc-211 directions are not part of this plan:

- treating a stored coset as an orbit representative under equal virtual-ligand occurrences;
- searching, minimizing, or accepting any of several occurrence permutations in equality,
  canonicalization, matching, pushout, or reaction application;
- requiring each frame-relative constraint to be invariant under a repeated-ligand stabilizer;
- allowing kindless stereo frames above `MAX_DEGREE` and carrying them through a separate general
  position-order fallback;
- requiring canonical-form construction to mutate storage in semantic prefix order: the returned
  value contains all three prefixes, while its selected id witness is applied before final
  reframing in that id space;
- requiring `Reframe` to borrow its receiver or promising a particular allocation count;
- absorbing `remap`, deleting `try_remap`, or manufacturing `Reaction::remap` through a materialized
  span merely as a consequence of complete-only canonicalization; and
- closing doc 209 as `Completed`. It remains `Superseded`; this document owns the moved work.

The historical discussions are intentionally retained. Their dated supersession notices identify
the exact repeated-ligand claims that no longer govern current design.
