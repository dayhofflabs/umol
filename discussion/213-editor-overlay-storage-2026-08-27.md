# 213 — Collapsing the editor's overlay storage

Status: Proposed
Date: 2026-08-27
Relates: [211](211-relation-frames-and-api-2026-08-26.md),
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
    Mutable(Vec<O::Entry>),
}
```

`Shared` holds the published overlay, so reads go to it; `Mutable` accumulates entries and is
positional, which is why the editor indexes where the overlay does not.

**Do not fold the two states into the overlay type itself.** The shared state is a real
copy-on-write win — an edit touching one atom republishes the other five overlays' handles with no
work — and a mutable variant on `Molecule`'s field would give the published type a state it must
never hold. The accumulating state belongs to the editor.

## `Overlays`

The shared surface of the six overlay storage types, in `traits.rs`. Settled, and specified here in
full because it has to be written along with everything else.

Public, carrying the members that were already `pub` and inherent on all six, plus two:

```rust
type Id: Copy + From<usize>;
type Entry: Clone;
type Attributes;

fn into_entries(self) -> Vec<Self::Entry>;
fn entry(&self, id: Self::Id) -> Self::Entry;                                    // new
fn count(&self) -> usize;
fn contains(&self, id: Self::Id) -> bool;
fn ids(&self) -> impl ExactSizeIterator<Item = Self::Id> + '_;
fn attributes(&self, id: Self::Id) -> &Self::Attributes;
fn attributes_mut(&mut self, id: Self::Id) -> &mut Self::Attributes;
fn incident_ids(&self, atom: AtomId) -> impl ExactSizeIterator<Item = Self::Id> + '_;
fn has_incident(&self, atom: AtomId) -> bool;
fn compact(&self, compaction: &GraphCompaction) -> (Self, Compaction<RelationId>);  // new
```

Members whose shape is the entity kind's own stay inherent: the participant accessors, `is_coincident`
and `coincident_id`, whose arguments are the participants themselves, `into_arc`, the node-level
accessors, `StereoBonds`'s bond-keyed incidence, and the crate-private `remap`, `glue` and
`attributes_iter_mut`. Nothing crate-private is published by the trait.

`ids` needs an explicit `+ '_`: in trait position it captures the `&self` lifetime where the
inherent `impl Trait` did not, and that propagates to the four `*Views::ids` methods.

## The open piece: `Participants`

The editor's six `*_equiv` old-state checks ask two questions of one entry — is this the entity you
mean, and does its value agree. The first is a participant comparison, each factor as a multiset,
with the attributes taking no part. It has no home today: the overlays' `is_coincident` takes an id
and loose per-kind participants, which the editor's accumulating state cannot supply.

Attempting it as a static comparison of two `Entry` values named the operation three times and got
it wrong three times, because the thing being compared has no type. It is not entry equality — the
caller builds an offered entry whose attributes are cloned and never read.

The fix is an associated type:

| overlay | `Participants` |
| --- | --- |
| aromatic system, multicenter bond | `Vec<AtomId>` |
| noncovalent bond | `[AtomId; 2]` |
| dative bond | `(Vec<AtomId>, AtomId)` — donors, acceptor |
| stereo atom, stereo bond | `(AtomId \| BondId, Vec<StereoLigand>)` — site, ligands |

With it, `Entry` stops being an arbitrary per-kind tuple and is `(Participants, Attributes)`
everywhere, and the comparison is on `Participants`, where it can be named for the question rather
than for whatever shape was at hand. The name is not settled: `coincident` and `is_coincident`
already differ by three characters while asking a search and a predicate respectively, so a third
member on that root compounds a collision rather than joining a family.

Two consequences:

- Four constructors nest one level deeper — dative, noncovalent and both stereo kinds. `new` and
  `into_entries` are `pub`, but every caller is inside `umol-graph-ir` and most are in test modules;
  nothing in `umol-graph`, `umol-io` or `umol-py` constructs or destructures an overlay. The span
  types mirror the overlays and want the same `Entry` shape or they drift from what they mirror.
- If `new` leaves the trait, `OverlayEditor`'s publish step still has to rebuild the overlay from
  entries, so construction must stay reachable generically — either a construction member under a
  name that is not `new`, or a `From<Vec<Self::Entry>>` bound on the wrapper.

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

## Handoff

Nothing here is in the tree. `Overlays` and its six impls were written and reverted along with the
rest, which is no great loss: the members are lifted verbatim out of the six inherent blocks, and
the specification above is what took the thinking.

The whole of it is to be done: `Overlays` and its impls, the `Participants` associated type with its
comparison member, `OverlayEditor`, the editor view change, and the migration of the editor call
sites — around 66 of them, every one a "trait not in scope" or an index-keyed read becoming an
entry.

None of it is required by doc 211's remaining stages.
