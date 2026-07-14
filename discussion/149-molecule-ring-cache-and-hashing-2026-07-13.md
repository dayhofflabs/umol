# 149 · MoleculeAst ring cache placement & hashing semantics

Status: Active
Date: 2026-07-13
Relates: 137 (Python-binding R1 audit — where this surfaced), 113/095 (lazy
canonicalization; structural vs canonical equality), 114 (interning; the future
immutable/finalized molecule form), 143 (vendored nauty — the canonical-numbering
oracle)

## Context

Doc 137's R1 fold-back ("make `MoleculeAst` hashable") surfaced that the lazily
populated ring cache is the *only* thing keeping `MoleculeAst` from being a plain
value type, and that "just add `Hash`" is not a local edit. Two separable
questions fell out, captured here:

- **A.** Should the ring cache live on `MoleculeAst` at all?
- **B.** How should `MoleculeAst` hash?

They are independent — removing the cache does not, by itself, make the type
hashable, and hashing does not depend on the cache. This doc records the current
state, the decisions, and the next steps.

## A. Ring cache — remove it from the AST

### Current state

- `MoleculeAst` carries `rings_cache: OnceLock<RingSet>`, filled lazily on the
  first call to `rings()` (Vismara relevant cycles, max ring size 22). Every other
  field is structural/value data (graph, per-atom/bond vectors, the overlay
  relation sets, molecule constraints).
- The cache is the sole reason `MoleculeAst` needs a hand-written `Clone` (which
  resets the cache to empty), a hand-written `PartialEq` (which skips the cache and
  compares the ten value fields), and interior mutability.
- Ring queries on the views come in two forms:
  - **cached:** `AtomView`/`BondView` `is_in_ring`, `ring_membership`,
    `ring_count`, `ring_size_count`, `ring_degree`, `ring_valence` — all route
    through `self.molecule.rings()` (the `OnceLock`).
  - **caller-supplied `RingSet` ("bring your own"):** `is_in_ring_from(&RingSet)`
    and `rings_from(&RingSet)` — the explicitly-intended model, but only two of the
    roughly six ring queries have a `_from` variant.
- The cache is never seeded from a precomputed set: every constructor initialises
  `OnceLock::new()`. `rings_with(family, max_size, filter)` is the uncached,
  owned-return enumerator.

### Decision

Introduce a `RingsView` and move all topology-derived ring queries onto it, then
delete the cache. This supersedes both the cache and the (only half-built)
caller-supplied `_from` surface.

`RingsView<'a>` owns a computed `RingSet` and borrows `&'a MoleculeAst`. Its shape
mirrors the entity namespace accessors (`AtomsView` / `BondsView`):

- ring-set level: `count()`, `ids() -> impl Iterator<Item = RingId>`,
  `iter() -> impl Iterator<Item = RingView>`, `get(id) -> Option<RingView>` — the
  existing `RingView` (one ring: `atoms` / `bonds` / `len` / relations) is reused;
- per-entity ring sub-views: `atom(id) -> RingAtomView`, `bond(id) -> RingBondView`.

`RingAtomView<'a>` / `RingBondView<'a>` hold `&RingSet` + `&MoleculeAst` + the id,
so the current view ring queries move onto them unchanged and stay **no-arg**:
`is_in_ring`, `ring_count`, `ring_size_count`, `ring_degree`, `ring_valence`,
`ring_membership`, plus a per-entity `rings()` iterator. `mol.rings()` is repointed
to compute and return a `RingsView` (owned `RingSet`, uncached);
`mol.rings_with(family, max, filter)` returns the same for a custom set.

Why this over the cache or the `_from` surface:

- **No cached data, no recomputation hazard.** The `RingSet` is computed once when
  the `RingsView` is created and held for its lifetime; there is no interior
  mutability and nothing to invalidate.
- **The extra hop is a signifier, not a wart.** `mol.rings().atom(id).ring_count()`
  makes visible that a ring enumeration is happening — a real, non-free operation —
  where `atom.ring_count()` hid it behind a lazy cache.
- **The validity contract becomes compile-enforced.** `RingsView` borrows `&mol`, so
  a ring set cannot be held across a structural mutation; the borrow checker enforces
  what the cache only assumed (topology invariant for the molecule's lifetime).
- **Consistent, read-only API.** The shape matches `AtomsView` / `BondsView`; rings
  are derived and never mutated, so no mutable accessor / `*ViewMut` is needed.
- **`MoleculeAst` becomes a plain value.** The cache was its only non-value field, so
  with it gone `Clone` / `PartialEq` / `Eq` derive.

### Next steps (staged, tree green after each stage)

- **S0 — add the `RingsView` type family** *(additive).* `RingsView<'a>`,
  `RingAtomView<'a>`, `RingBondView<'a>`, with the surface above (ring-set-level
  accessors returning `RingView`; the no-arg per-entity ring queries). Expose it via a
  temporary accessor (e.g. `rings_view()` / `rings_view_with(...)`) alongside the
  still-present cache. Tests for the new surface. The cache and the existing
  `AtomView` / `BondView` ring methods are untouched — green.
- **S1 — migrate consumers to `RingsView`** *(breaking → green).* Repoint every
  consumer of the cached ring queries — `AtomView` / `BondView`
  `is_in_ring` / `ring_count` / …, the two `_from` methods, `derive_constraints`, and
  any perception / aromatic-system step reaching `self.molecule.rings()` — to compute
  a `RingsView` once and query it. Per-atom / per-bond loops (e.g. constraint
  derivation) compute the `RingsView` up front and reuse it, never per element.
- **S2 — remove the cache and the old query surface** *(breaking → green).* Delete
  `rings_cache`, the cached `rings() -> &RingSet`, and the `AtomView` / `BondView`
  ring methods (including the `_from` pair, now subsumed by the sub-views). Rename the
  S0 accessor to `rings()` (now `-> RingsView`). Replace the hand-written `Clone` /
  `PartialEq` on `MoleculeAst` with derives. Payoff: `MoleculeAst` is a plain value.

Design note: `ring_valence` / `ring_degree` need the atom's incident bonds, so the
sub-views hold `&MoleculeAst` alongside `&RingSet`; pure membership / count queries
use the `RingSet` alone.

### API shapes

`RingsView<'a>` — borrows `&mol`, owns the `RingSet`. Ring-collection surface
mirrors `AtomsView`, plus per-entity sub-views and the ring-system queries:

    count() -> usize
    ids() -> impl Iterator<Item = RingId>
    iter() -> impl Iterator<Item = RingView<'_>>
    get(RingId) -> Option<RingView<'_>>
    contains(RingId) -> bool
    atom(AtomId) -> RingAtomView<'_>
    bond(BondId) -> RingBondView<'_>
    family() -> RingFamily · max_ring_size() -> usize

There is no plural `atoms()` / `bonds()` on `RingsView` — per-atom/bond questions go
through `atom(id)` / `bond(id)`. (No current caller wants "all atoms in some ring";
it can be added if one appears.)

`RingAtomView<'a>` / `RingBondView<'a>` — hold `&RingSet` + `&MoleculeAst` + the id;
carry exactly the ring queries that live on `AtomView` / `BondView` today, now no-arg:

    is_in_ring() -> bool
    ring_count() -> ValueAst
    ring_size_count(u8) -> ValueAst
    ring_membership(RingScope) -> ValueAst
    ring_degree() -> ValueAst                       // atom only
    ring_valence() -> ValueAst                      // atom only
    rings() -> impl Iterator<Item = RingView<'_>>   // the rings through this atom/bond
    smallest_ring_size() -> Option<usize>

`RingView<'a>` — gains a backref `&'a RingSet` and is brought into line with the
entity views:

    id: RingId
    len() -> usize
    is_empty() -> bool
    atom_ids() -> impl Iterator<Item = AtomId>      // was atoms() -> &[AtomId]
    atoms() -> impl Iterator<Item = AtomView<'_>>   // matches AromaticSystemView::atoms
    bond_ids() -> impl Iterator<Item = BondId>
    bonds() -> impl Iterator<Item = BondView<'_>>
    contains_atom(AtomId) -> bool
    contains_bond(BondId) -> bool
    // relations, now id-keyed via the backref (was: &RingView):
    shared_atoms(other: RingId) -> Vec<AtomId>
    shared_bonds(other: RingId) -> Vec<BondId>
    relation(other: RingId) -> RingRelation
    is_spiro(other: RingId) -> bool
    is_fused(other: RingId) -> bool
    is_bridged(other: RingId) -> bool
    spiro_neighbors() -> impl Iterator<Item = RingId>
    fused_neighbors() -> impl Iterator<Item = RingId>
    bridged_neighbors() -> impl Iterator<Item = RingId>

Settled calls:

- **`RingView::atoms()` / `bonds()` return iterators, not a nested collection view.**
  The convention for an entity's members is the `atom_ids()` + `atoms()` iterator pair
  (`AromaticSystemView` / `MulticenterBondView`); a ring's members follow it, with
  `len()` for count and `contains_atom` / `contains_bond` for membership. A nested
  collection view would be inconsistent with every other entity view and is out of
  scope; introducing one would be a separate, codebase-wide convention change.
- **`RingSet` keeps its id-pair primitives** — `shared_atoms(a, b)`, `relation(a, b)`,
  and the binary relation accessors `spiro_neighbors(i)` / `fused_neighbors(i)` /
  `bridged_neighbors(i)`; `RingView`'s relation methods delegate to them via the
  backref. `RingView::shared_atoms(&RingView)` / `shared_bonds(&RingView)` are removed.
- **Ring-system component analyses are dropped from the view surface.**
  `fused_components` (on `RingsView`) and `fused_component` (on `RingView`) are not
  exposed: they have no consumer today, and are cheap to reconstruct on demand from the
  binary relation accessors (`relation`, `*_neighbors`) plus union-find. If they become
  useful, the coherent form is likely a single `ring_systems()` over the combined
  fused + bridged connectivity rather than per-relation `X_components()` variants. The
  existing consumer-less `RingSet::fused_components` / `fused_component` can be removed
  in S2 or kept as internal utilities.

## B. Hashing — two tiers, mirroring equality

### Agreed semantics (lazy canonicalization)

Equality is already two-tiered and stays that way; hashing follows the same shape:

- **Structural tier (cheap, the default):** `==` is the field-by-field structural
  compare that exists today; `Hash` is the matching cheap structural hash. Both are
  numbering-dependent — two orderings of the same molecule compare unequal and hash
  differently. This is the default identity used for ordinary `HashSet`/`HashMap`
  keys and Python `set`/`dict`.
- **Canonical tier (opt-in):** `canonical_eq` (exists on the AST value types today)
  and a new `canonical_hash`. Isomorphism-invariant — the identity used for
  deduplicating molecules across a reaction network. Computed lazily on demand, not
  eagerly, so a mutate-heavy workload does not re-canonicalise on every edit.

This split is deliberate: the structural tier is cheap and numbering-sensitive; the
canonical tier is the explicit path when chemical identity is what's wanted. Whether
structural-default / canonical-opt-in is the right ergonomics (versus
canonical-by-default) is left for alpha feedback; it is not being changed ahead of
that, and lazy canonicalization is the better default for the mutate-heavy case
regardless.

### Implementation notes

- **Structural `Hash` must agree with structural `==`.** Since `==` compares the
  value fields, `Hash` must hash the same fields such that `a == b` implies equal
  hashes. This requires `Hash` on the field types — the graph, the relation-set
  containers, molecule `Constraints`, and the entity ASTs — each consistent with its
  own `==`. Mostly mechanical, with one real check per type: wherever a container's
  `==` is order-independent (logical membership regardless of internal storage
  order), its `Hash` must digest a canonical internal form, not raw storage.
  Once the ring cache is gone (Part A) and the fields are `Hash`, `MoleculeAst` can
  simply `derive(Hash)`.
- **The Python side is not blocked by immutability.** A hand-written
  `fn __hash__(&self) -> u64` works on the mutable pyclass; only the
  `#[pyclass(hash)]` auto-derive requires `frozen`. A manual `__hash__` (structural)
  mirrors the existing hand-written structural `__eq__`, and the mutable-dict-key
  caveat is inherent and already accepted for structural `==` on a mutable value.
- **`canonical_hash` is the relational extension of `canonical_eq`.** `canonical_eq`
  is settled on the AST value types, but a molecule's identity is dominated by its
  graph (topology plus numbering). Canonical molecule identity therefore needs graph
  canonicalization — a canonical numbering (via the coloring / the nauty oracle) plus
  remapping the per-entity data into that numbering — followed by a structural
  compare/hash of the canonical form. `canonical_hash` hashes that canonical form and
  is consistent with `canonical_eq` by construction. This graph-canonicalization +
  remapping step is the part that does not exist yet, and is the right tool for
  network-node deduplication.

### Relation to interning (114)

The interning plan introduces an immutable / finalized molecule form (a network
node, an interned handle). That form is the eventual natural home for both a stored
canonical identity and, if it returns, a ring set — topology and identity are stable
there, so caching and canonical hashing are unproblematic. Nothing in this doc
commits to that form; the structural tier and `canonical_hash` are defined on the
mutable `MoleculeAst` and carry over unchanged if a finalized form is later added.

## Sequencing across A and B

- Part A (ring-cache removal) is a self-contained value-cleanliness change,
  independent of hashing; it can proceed first.
- Part B's structural tier is a broad but mechanical `Hash`-on-the-field-types
  ripple through graph-core and the AST; it can land independently, and is cleaner
  once A is done.
- Part B's `canonical_hash` is a focused but deeper piece that depends on graph
  canonicalization + remapping; it lands when a real consumer (network dedup) needs
  it, and pairs with the nauty work.
