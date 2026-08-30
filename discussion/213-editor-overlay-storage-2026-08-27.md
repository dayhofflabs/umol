# 213 — Collapsing the editor's overlay storage

Status: Proposed
Date: 2026-08-27
Relates: [211](211-relation-frames-and-api-2026-08-26.md),
[214](214-aggregate-frame-semantics-2026-08-28.md),
[data-type guide](../docs/development/data-types.md),
[nomenclature guide](../docs/development/nomenclature.md)

## Purpose

`MoleculeEditor` holds its six overlays through three copy-on-write wrappers — `FixedSetStorage`,
`VarSetStorage` and `FixedVarSetStorage` — carrying 31 methods between them. They exist only because
they are generic over the *storage shape*, and most of what they carry now duplicates reads the six
overlay types already own.

## The shape

One wrapper parameterised by the overlay type replaces the three:

```rust
enum OverlayEditor<O: Overlays> {
    Shared(O),
    Mutable(Vec<(O::Participants, O::Attributes)>),
}
```

`Shared` holds the published overlay, so reads go to it; `Mutable` accumulates entries and is
positional, which is why the editor indexes where the overlay does not.

**Do not fold the two states into the overlay type itself.** The shared state is a real
copy-on-write win — an edit touching one atom republishes the other five overlays' handles with no
work — and a mutable variant on `Molecule`'s field would give the published type a state it must
never hold. The accumulating state belongs to the editor.

## `Overlays`

**Open before implementation.** `Overlays` is not an approved trait name: the plural denotes all
overlays, while each proposed implementation represents one entity set. Whether the editor needs a
generic set trait at all, whether such a trait should be public, and what it should be called remain
unsettled. The surface below records the earlier design sketch rather than an implementation-ready
API.

The earlier sketch placed the shared surface of the six overlay storage types in `traits.rs` and
specified it as follows:

Public, carrying the members that were already `pub` and inherent on all six, plus three:

```rust
type Id: Copy + From<usize>;
type Participants: Clone;
type Attributes: Clone + Normalize + FrameTransport<Action = Self::LocalAction>;
type LocalAction;

fn into_entries(self) -> Vec<(Self::Participants, Self::Attributes)>;
fn entry(&self, id: Self::Id) -> (Self::Participants, Self::Attributes);         // new
fn count(&self) -> usize;
fn contains(&self, id: Self::Id) -> bool;
fn ids(&self) -> impl ExactSizeIterator<Item = Self::Id> + '_;
fn attributes(&self, id: Self::Id) -> &Self::Attributes;
fn attributes_mut(&mut self, id: Self::Id) -> &mut Self::Attributes;
fn incident_ids(&self, atom: AtomId) -> impl ExactSizeIterator<Item = Self::Id> + '_;
fn has_incident(&self, atom: AtomId) -> bool;
fn compact(&self, compaction: &GraphCompaction) -> (Self, Compaction<RelationId>);  // new

fn alignment_action(
    from: &Self::Participants,
    to: &Self::Participants,
) -> Option<Self::LocalAction>;                                                  // new
```

Members whose shape is the entity kind's own stay inherent: the participant accessors, `is_coincident`
and `coincident_id`, whose arguments are the participants themselves, `into_arc`, the node-level
accessors, `StereoBonds`'s bond-keyed incidence, and the crate-private `remap`, `glue` and
`attributes_iter_mut`. Nothing crate-private is published by the trait.

`ids` needs an explicit `+ '_`: in trait position it captures the `&self` lifetime where the
inherent `impl Trait` did not, and that propagates to the four `*Views::ids` methods.

## Participant alignment

The editor's six `*_equiv` old-state checks ask two questions of one entry — is this the entity you
mean, and does its value agree in the stored participant frame. The first cannot be a Boolean
participant comparison. It must retain the unique local frame action from the offered participants
to the stored participants so that the second can transport the offered attributes before comparing
their normal forms.

This operation has no home in the three storage-shaped editor wrappers. Their relation sets know
only factors and multisets; they do not know which factor bears the entity frame or, for a stereo
bond, that its ligand factor consists of two endpoint blocks. The six overlay aggregates do
know that structure, so `Overlays::alignment_action` owns the per-kind derivation. Its direction is

```text
to[i] = from[action[i]]
```

and it returns `None` unless `from` and `to` are two frames of the same entity under the
structured-participant semantics of that entity kind. Ordinary unordered factors use
`DynPermutation`; stereo factors use the bounded `Permutation`. Stereo-bond alignment permits
permutations within each endpoint block and exchange of the two complete blocks, but not movement
of one ligand across the endpoint boundary. Integrity-valid frames have distinct complete
participant values, so the admissible action is unique.

The participant and local-action types are:

| overlay | `Participants` | `LocalAction` |
| --- | --- | --- |
| aromatic system, multicenter bond | `Vec<AtomId>` | `DynPermutation` |
| noncovalent bond | `[AtomId; 2]` | `DynPermutation` in `S_2` |
| dative bond | `(Vec<AtomId>, AtomId)` — donors, acceptor | `DynPermutation` on donors |
| stereo atom | `(AtomId, Vec<StereoLigand>)` — site, ligands | `Permutation` on ligands |
| stereo bond | `(BondId, Vec<StereoLigand>)` — site, ligands | `Permutation` in `S_2 wr S_2` |

`OverlayEditor<O>` then owns one `entry_framed_eq` old-state check for every overlay kind. It
obtains the stored entry, calls `O::alignment_action(offered, stored)`, transports the offered
attributes with `FrameTransport`, and calls `normalized_eq` against the stored attributes. It
neither selects a representative frame nor implements `Reframe`; this is pairwise alignment
between two supplied local frames. Participant or transport incompatibility returns `false`, and
the transaction retains its existing `TransactionError::OldStateMismatch` boundary.

The current six `MoleculeEditor::*_equiv` methods therefore retire with the wrapper migration. Their
per-kind knowledge moves into the six `alignment_action` implementations rather than remaining as
six copies of the transaction operation.

Two consequences:

- Four constructors nest one level deeper — dative, noncovalent and both stereo kinds. `new` and
  `into_entries` are `pub`, but every caller is inside `umol-graph-ir` and most are in test modules;
  nothing in `umol-graph`, `umol-io` or `umol-py` constructs or destructures an overlay. The span
  types mirror the overlays and want the same `(Participants, Attributes)` entry shape or they drift
  from what they mirror.
- If `new` leaves the trait, `OverlayEditor`'s publish step still has to rebuild the overlay from
  entries, so construction must stay reachable generically — either a construction member under a
  name that is not `new`, or a
  `From<Vec<(Self::Participants, Self::Attributes)>>` bound on the wrapper.

## Editor views

Three of the six editor views hold `&'a [NodeId]` and convert on read, which worked while the
mutable state also held node ids. Overlay entries are in graph-IR ids, so the two states no longer
agree and no borrowed slice covers both.

The immutable views therefore own their participants, `Vec<AtomId>`. The `*ViewMut` variants are
always constructed after `materialize`, so they borrow `&'a [AtomId]` from the entry. This is one
allocation per editor view construction, on the editing path: the `*EditorView` types are built at
nine sites, all in `editor.rs`, and the published read path is `*View` on `Molecule`, untouched. The
editor already copies a whole entry list into the mutable state on the first write to an overlay.

The stereo views need no change — `StereoLigand` is the same type in both states — and the
noncovalent views already hold `[AtomId; 2]` owned.

## Evidence

Retain every editor and transaction assertion, including rollback. Assert that publishing an
unmodified overlay returns the same handle without rebuilding, and that a batch of pushes
materialises once rather than per operation.

Exercise `alignment_action` and the generic old-state check for every overlay kind: identity and a
legal reordering succeed, a participant mismatch fails, and stereo-bond cases distinguish
within-block and complete-block exchange from an illegal cross-block movement. Use a
position-sensitive aromatic or multicenter payload and a stereo payload so success demonstrates
transport rather than only frame-invariant comparison.

## Handoff

Nothing here is in the tree. `Overlays` and its six impls were written and reverted along with the
rest, which is no great loss: the members are lifted verbatim out of the six inherent blocks, and
the specification above is what took the thinking.

The whole of it is to be done: `Overlays` and its impls, the `Participants` and `LocalAction`
associated types with `alignment_action`, `OverlayEditor`, the generic framed old-state check, the
editor view change, and the migration of the editor call sites — around 66 of them, every one a
"trait not in scope" or an index-keyed read becoming an entry.

None of it is required by doc 211's remaining stages.
