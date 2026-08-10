# 119 — umol-perm defensiveness and API-consistency review (2026-06-21)

## Trigger

The stage-E6 substructure-matching proptest (umol-ast) feeds `MoleculeAst`s with
structurally-invalid stereo overlays produced by `molecule_strategy`: a ligand
frame whose size differs from the stereo kind's degree (e.g. a Tetrahedral center
with two ligands), and a coset index `Lit(0..=6)` that can exceed the kind's coset
count. The stereo coset post-filter calls `CosetSpace::coset_for`/`reindex`, which
panics — `coset.rs:114` (the `parent.contains` assert) and, after guarding that,
`coset.rs:106` (`unindex` out of range). Resolver-validated molecules never hit
this, but it prompted a crate-wide look at where umol-perm should be defensive.

## Principle

Low-level libraries should not be reflexively defensive; input validation generally
belongs in higher layers. The exception is when the violated condition is much
cheaper or clearer to detect at the low level than at the call site. For
`CosetSpace`, "index out of range" / "permutation not in the parent group" is a
one-liner against `self.count()` / `self.parent`, whereas the umol-ast caller would
have to re-derive stereo kind → degree → ligand count to pre-check. Consistency is
the second rule: guarding one method of a family while its siblings keep panicking
is itself a defect.

## Findings

### A. `CosetSpace` index/permutation family is inconsistent (primary)

- `index(perm)` panics via a missing `HashMap` key when `perm` is outside the
  parent group (`coset.rs:100-102`); no validation. `reindex` *does* validate the
  analogous condition for its `relabeling` (`coset.rs:114`), but `index`,
  `merge_under` (`coset.rs:148`), and `observable_coset` do not.
- No `u32` index is range-checked: `unindex` (`coset.rs:106`), `reindex`,
  `enantiomer` (`coset.rs:123`), `observable_coset` (`coset.rs:160`) all bare-index
  `representatives` or the merge vector.
- `Coset::new(key, index)` (`class.rs:208`) stores an unvalidated index, deferring
  the panic to an arbitrary later `unindex`.

### B. `compose` degree check is debug-only (correctness hole)

`Permutation::compose` uses `debug_assert_eq!` on degree (`permutation.rs:74`); in
release a degree mismatch silently composes over the fixed `[u8; 6]`. It propagates
through `group`, `coset`, and `oriented` — the only invariant dropped in release.

### C. Mixed degree-validation in `permutation.rs`

`identity`/`from_image`/`from_cycles` assert (and `from_image` fully validates the
image); `compose` debug-asserts; `unrank` (rank ≥ degree! → `Vec::remove` panic,
`permutation.rs:158`; also no `degree ≤ MAX_DEGREE` entry assert) and `apply`
(i ≥ 6 panics; degree ≤ i < 6 silently returns a fixed point, `permutation.rs:64`)
do neither. `act` bare-indexes a too-short slice (`permutation.rs:69`).

### D. Efficiency (no blockers; a few per-call allocations)

- One-time construction dominates: `space()` interns one `CosetSpace` per
  `ClassKey` (`LazyLock<Mutex<HashMap>>` + `Box::leak`, `class.rs:185`); group
  generation enumerates ≤ 720 elements once.
- Per call: `index` is O(|R|) (`coset_rep` scans R, `coset.rs:168`);
  `OrientedPermutationGroup::elements`/`star_orbit_of` allocate the full element
  vector each call (`oriented.rs:160,190`); `OrientedPermutationGroup::generate`
  enumerates the proper subgroup twice (`oriented.rs:135`).
- All sizes are bounded by degree ≤ 6; nothing here is a bottleneck.

### F. `MAX_DEGREE` is private and underused

`MAX_DEGREE = 6` (`permutation.rs:13`) is the permutation-array bound and is used
consistently inside `permutation.rs`, but it is private and not referenced
elsewhere: `ClassKey::from_str` does not reject `degree > MAX_DEGREE`, and doc prose
hardcodes "≤ 6". A literature review suggests the bound may need to rise to 8
(7/8-coordinate geometries), so the bound must have a single named home.

Note the distinction: the `6`/`5`/`4` literals in `coset.rs`/`class.rs` are
**geometry degrees** (octahedral 6, trigonal-bipyramidal 5, square-planar 4), not
the max bound — octahedral stays degree 6 if `MAX_DEGREE` rises — so they are not
`MAX_DEGREE` and stay as-is (they may be named per geometry separately).

Plan: promote `MAX_DEGREE` to `pub(crate)`; use it in `from_str`'s degree check;
fix prose. Raising the bound then stays a one-line change.

### E. Smaller items

- Implicit MSRV: `usize::is_multiple_of` (`permutation.rs:108`) and
  `Option::is_none_or` (`oriented.rs:129`) need a recent stdlib though the edition
  is 2021; no `rust-version` is set.
- `ClassKey::from_str` accepts out-of-range degrees (e.g. `Symmetric(7)`), failing
  only later in `build`.
- The `space()` registry mutex poisons permanently if any `build()` panics
  (`class.rs:190`).
- TB/OH enantiomer pairing is flagged unverified against the OpenSMILES @/@@
  numbering (`class.rs:122`) — a pre-existing open item, out of scope here.

## Recommendations

### Guard low, return `Option`, consistently — the `CosetSpace` family (finding A)

Make the index/permutation-taking methods total:

| method | → | guard returns `None` |
|---|---|---|
| `unindex(index)` | `Option<Permutation>` | `index >= count` |
| `index(perm)` | `Option<u32>` | `perm.degree() != degree` or coset not numbered |
| `reindex(index, relabeling)` | `Option<u32>` | index out of range or `relabeling` ∉ parent |
| `enantiomer(index)` | `Option<u32>` | index out of range |
| `observable_coset(index, …)` | `Option<u32>` | index out of range (or a fluxional gen ∉ parent) |
| `coset_rep(perm)` | `Option<Permutation>` | `perm.degree() != degree` |

`count`/`degree`/`group`/`improper`/`is_chiral`/`merge_under` keep their signatures;
the latter two use the now-`Option` ops internally with `.expect()` since `0..count`
is valid by construction. `Coset::new` validates its index (open question on shape,
below). The umol-ast caller's `coset_for(...)?` and `coset_matches` (already `false`
when the meet is `None`) then handle malformed stereo with no re-validation.

### Tighten to a real precondition (panic, not `Option`) — findings B, C

These are programmer-error preconditions, not runtime-data conditions, so a clear
assert is the right tool (an `Option` would be too intrusive for the hot algebra):

- `compose`: `debug_assert_eq!` → `assert_eq!` on degree — one `u8` comparison,
  closes the release hole.
- `apply`: assert `i < degree` (or formally adopt the fixed-point-fill as the
  contract and document it).
- `act`: assert `items.len() >= degree`.
- `unrank`: assert `rank < degree!` and add the missing `degree ≤ MAX_DEGREE` entry
  assert.

### Leave to higher layers / note only

- MSRV: set `rust-version` to the real minimum or replace the two recent-stdlib
  calls. Hygiene; decide separately.
- `ClassKey::from_str` degree range, registry poisoning, TB/OH pairing: recorded,
  not part of this change.

## Proposed plan

1. `coset.rs`: the `Option` family + the perm-membership / index-range guards.
2. `permutation.rs`: `compose` debug_assert → assert; `apply`/`act`/`unrank`
   precondition asserts.
3. Ripple: umol-ast `stereo.rs` (`coset_apply_permutation`/`coset_meet`/`coset_join`
   propagate `Option`) and `symmetry.rs` (propagate or `.expect()` on resolved
   input); the substructure matcher needs no change.
4. `Coset::new` validation.
5. Tests: umol-perm unit cases for the new `None` paths; the substructure proptest
   then passes with no matcher-side checks.

## Resolved (2026-06-21)

- `Coset`: infallible **checked** `new` (asserts `index < count`) plus
  `new_unchecked`.
- `compose`/`apply`: full `assert` (and `act`/`unrank` likewise).
- MSRV: out of scope here. `is_multiple_of` stays (clippy's `manual_is_multiple_of`
  prefers it). A `rust-version` / MSRV policy should be set umol-wide, not patched
  in one corner — tracked separately.

## Implemented (2026-06-21)

All five plan steps landed. `coset.rs` `index`/`unindex`/`reindex`/`enantiomer`/
`observable_coset`/`coset_rep` return `Option` (cheap range / perm-membership
checks); `permutation.rs` `compose`/`apply`/`act`/`unrank` use `assert`; the umol-ast
ripple propagates `Option` through `coset_apply_permutation` (and `.expect()` on
resolver-validated input in `raise.rs`/`symmetry.rs`); `Coset::new` is checked with
`new_unchecked` alongside.

The substructure matcher needs no defensive code: with the `Option` family a
malformed generated stereo overlay yields no match instead of a panic. The planted
proptest dropped its `!is_empty` self-match assertion — `molecule_strategy` can
emit stored constraints inconsistent with the derived topology (e.g. a stored valence
differing from the bond-derived one), so a self-match is not guaranteed; the surviving
cross-strategy / cross-algorithm agreement invariant holds.

`MAX_DEGREE = 6` is now a named `pub(crate)` const in `permutation.rs`, used by
`unrank` and `ClassKey::from_str`'s degree-range check (may rise to 8 for 7/8-coordinate
geometries — single named home).

Verified: `cargo test --workspace --all-features` (all suites pass) and
`cargo clippy --workspace --all-features --all-targets` clean for the touched crates.
