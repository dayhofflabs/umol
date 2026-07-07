# 138 — Constraint container API: the map model (reversed)

Status: Reversed — map model implemented, then reverted. Flat model kept.
Date: 2026-07-06

## Outcome and lessons (2026-07-06)

The map model below (`RingMembership` holding `BTreeMap<RingScope, ValueAst>`) was built end
to end — container, delta, transact, DSL, binding — and **reverted**. It does not earn its
keep. The lessons, so this isn't re-attempted:

- **The map relocates the by-key API; it doesn't remove it.** Deltas, composition, and
  serialization operate per **(kind, sub-key)**, not per kind. The map stores per **kind**,
  so every boundary needs an *expand* (map → single-key units) on the way out and a
  *collect* (units → map) on the way in. The "clean keyed-unit API" that expand/collect
  needs (`unit_at`/`set_unit`/`remove_unit`/`key_units`) is the old `get_by_key`/
  `remove_by_key`/`iter` renamed. So you keep the by-key surface **and** pay a round-trip.

- **The flat single-key model was already uniform at the boundary.** With ring stored as N
  single-key entries (one per scope), every constraint is one keyed unit; diff / apply /
  transact / serialize / parse all go through `iter`/`get_by_key`/`remove_by_key`/`add` with
  **zero single-vs-map branching**. Every special-case the map forced — in `delta.rs`,
  `transact.rs`, the DSL, the binding — traced to one decision: making ring a map.

- **The map bought nothing the flat model lacks.** Canonical `Eq`/`Hash` come from a
  key-sorted `Vec`; a map-shaped read is a `ring_membership()` *view* over the flat entries;
  per-`(kind, sub-key)` uniqueness is the natural guarantee; flat serialization
  (`#f(..)'#f(..)=`, the preferred form) falls out of iterating flat entries — the map had
  to *expand* to produce it.

- **What the container simplification actually is** (flat model, small change, not this
  refactor): make `add` meet by full key (find same `key()`, meet-or-insert) → removes
  `is_unique` and the canonicalize same-scope dedup; expose typed accessors / a
  `ring_membership()` view → `get_all` goes private; keep `get_by_key`/`remove_by_key`/`iter`
  as the uniform boundary. Open follow-ups: how the by-key/by-kind split reads at the flat
  container level; where `is_unique` is used across the containers and how it's replaced;
  whether any `*Constraint::` variant is still special-cased in delta/transact.

The rest of this doc is the record of the map design and why it was chosen, kept for context.
Everything under "## Why this doc" onward describes the reverted approach.

## Flat cleanup plan (current direction, 2026-07-06)

Keep the flat single-key `Vec<Constraint>` storage. This cleanup kills *kind* as an
*addressing* input (the mirror image of the map model's failed attempt to kill *key*),
unifies `add` on meet, and gives the delta boundary one transactional primitive. Nothing
here is map-shaped, and every entry stays a single-key unit, so the boundary is uniform with
no per-variant special-casing.

### Scope — `AtomConstraints` only (S1–S4)

S1–S4 change **`AtomConstraints` alone**; the other containers stay on the current API and are
a **separate replication** (below). This is feasible because the delta/transact boundary is
*per-family hand-written* — `apply_constraint`, `apply_modify_*_constraint`, and
`constraint_key` are one method per family — and the only shared piece, the
`diff_field_ops!`-generated `diff_constraints`, uses only `constraints.iter()` and
`Constraint::key()`, neither of which changes. So atom moves to the new surface without
touching bond / dative / aromatic / multicenter / noncovalent / stereo. The `*Ast` rename
(S5) is the one cross-cutting step and waits until every family is converted.

The flat model needs **no `*Dsl` boundary type** (unlike the reverted map model, which needed
an `AtomConstraintsDsl` to collect/expand the map for flat serialization). The only DSL in
scope is the *lowering* functions in `dsl/atom.rs` that call the container's by-kind/`is_unique`
API — `raise_atom_constraints`, `lower_atom_constraints`, the dup-check (S3b/S4b). The
per-constraint `AtomConstraintDsl` (EDN) and the string `fmt_constraint`/parsers are unaffected.

### Target container API (`AtomConstraints`; peers get the same surface in replication)

Read — by full key only:
- `get(key) -> Option<&Constraint>`
- `contains(key) -> bool`
- `iter()`

Assert — monotone, **fallible**:
- `add(constraint) -> Result<(), Contradiction>` — meet the value at `constraint.key()` via
  `Lattice::meet`; `Err` on incompatible meet. Fallible because it must work for every
  payload and most have no ⊥ (`BooleanAst` is `true`/`false`/`undetermined`, no bottom) — the
  old empty-`LitSet` sentinel was `ValueAst`-only and is dropped.

Overwrite — non-monotone:
- `set(constraint)` — overwrite at `constraint.key()`; a vacuous constraint removes.
  Infallible. (Takes a whole `Constraint`, not `(key, value)`: the payload type varies by
  kind — `ValueAst` / `AromaticValenceAst` / `MulticenterValenceAst` / `TetrahedralStereoAst`
  — so `set(key, value)` would need a payload union whose tag duplicates the key and can
  mismatch it; a `Constraint` is that pair with the tag shared once.)
- `update(&other)` — bulk overwrite from another container (overlay; vacuous ⇒ remove).
  Infallible. Replaces the `remove`+`add` loop in `AtomAst::update`.
- `compare_and_set(old: Option<Constraint>, new: Option<Constraint>) -> Result<(), Contradiction>`
  — verify the value at the key equals `old` (`canonical_eq`; `None` = expect absent), then
  set/remove `new`. The single delta apply/undo primitive.

Remove:
- `remove(key) -> Option<Constraint>` — unconditional delete.

Bulk:
- `retain(pred)`

### Removed / renamed / kept
- **Removed:** `is_unique` (method + tests); by-kind `get`/`get_mut`/`contains`/`remove`
  *taking a `Kind`*; `get_all`, `remove_all`; the mutable accessors `get_by_key_mut` and
  `iter_mut` (both test-only — a `&mut Constraint` lets a caller mutate the sub-key in place and
  corrupt the key-sort; all mutation goes through `set`/`add`/`compare_and_set`); `remove_entry`
  (→ `compare_and_set`); `AtomConstraintKey::kind()` (+ its test).
- **Kept private, renamed `find_by_key`→`find`:** `find(key) -> Result<usize, usize>` — the
  binary-search-by-key that `add`/`set`/`get`/`remove`/`contains`/`compare_and_set` are all
  built on. It is the primitive they fold *into*, not something that folds away; the `_by_key`
  suffix is redundant now that key is the only search.
- **Renamed:** `get_by_key`→`get`, `remove_by_key`→`remove`, `contains_key`→`contains`. Key is
  the only addressing left, so the `_by_key` suffix is redundant; the names free up because
  the by-kind versions are removed.
- **Kept:** `Constraint::kind()` (DSL tags, default dispatch), the `AtomConstraintKind` type, the
  `AtomConstraintKey` type, all `*Dsl`. `Kind` lives only on the constraint side; `Key` lives
  only on the addressing side; they never cross.
- **Note on `kind()`:** `AtomConstraintKind` + `AtomConstraint::kind()` are the strum
  `EnumDiscriminants` on the enum (and drive `AtomConstraintKind::iter()`), so they stay.
  `AtomConstraintKey::kind()` was a *hand-maintained* 13-arm `match` duplicating that
  projection with zero real callers once by-kind addressing is gone — deleting it removes
  upkeep, not capability.

### Behavior change (pin with tests before it lands)
Meet-`add` makes duplicate assertions conjunctive: `#V4#V3` → contradiction (was
last-wins/`DuplicateAtomPredicate`); `#V*#V3` → `#V3` (vacuous meets away). Confirmed
acceptable — also lets repeated tags join.

### Call-site migration — atom slice (audited 2026-07-06, whole workspace)

Atom sites are the S1–S4 slice. Peer-family sites (bond/dative/aromatic/multicenter/stereo)
are the replication, shown only where they prove the pattern is identical.

**`get_all(Kind::RingMembership)` — 3 sites, all internal ring accessors:**
| site | now | after |
|---|---|---|
| `atom.rs:595` `ring_memberships()` | `get_all(Kind::RingMembership).filter_map(..)` | delete `ring_memberships`/`ring_membership_value`; `ring_count()` = `get(Key::RingMembership(All))`, `ring_size_count(s)` = `get(Key::RingMembership(Size(s)))` |
| `bond.rs:151`, `dative.rs:152` (replication) | same shape | same, their own `Key` |

**`remove_all` — 0 call sites. Delete.**

**`get(Kind)` / `get_mut(Kind)` (atom):**
- internal typed accessors — `atom.rs:531–631` (`valence`…`tetrahedral_stereo`): `self.get(Kind::X)` → `self.get(Key::X)` (single-valued: literal key ≡ literal kind).
- external — `umol-graph/ops/aromaticity.rs:288` `.get(Kind::AromaticValence)?` → `.aromatic_valence()` (existing accessor).

**`contains(Kind)` (atom):**
- DSL default-fill `raise_atom_constraints` (`dsl/atom.rs:755–796`): `!contains(kind)` → `!contains(Key::X)` per branch; the pattern-only arm (`TotalValence | Degree | … | RingMembership`) stays a no-op — **ring membership is never filled**, so no branch ever needs `contains(Key::RingMembership(?))`.
- DSL dup-check `dsl/atom.rs:548`: **deleted** (meet-`add` subsumes it).
- `Lattice` meet/join `atom.rs:902/903`: folded into the per-key rewrite below.
- external — `umol-py/constraint.rs:462–464` (tests): → `contains(Key::X)`.

**`remove(Kind)` (atom):** the only two callers are the `lower_atom_constraints` `remove(kind)`
(gone with the `retain` rewrite below) and external `umol-graph/ops/transform/kekulizer.rs:88`
`.remove(Kind::AromaticValence)` → `.remove(Key::AromaticValence)`.

**DSL default-elision rewrite (`lower_atom_constraints`, `dsl/atom.rs:864`):** replace the whole
`for kind in Kind::iter()`/`match` (which used `get(kind)`+`remove(kind)`) with
`constraints.retain(|c| !is_elidable_default(c, cfg))` + a new exhaustive
`is_elidable_default(&AtomConstraint, &AtomDefaults) -> bool` (`match c`, **no** `_`;
`RingMembership(_) => false`). The exhaustiveness guard moves from `match kind` to `match c`.
(Its counterpart `raise_atom_constraints` keeps the kind-loop — fill must see *absent* kinds.)

**`Lattice` meet/join restructure (atom, `atom.rs`≈770–935):** today per-kind (`contains(kind)`
+ typed accessors) plus a separate ring-scope loop. → one pass over `self.iter()` ∪
`other.iter()` **keys**, meet/join per key (absent = ⊤). Unifies scalar and ring (ring scopes
are just keys) and deletes the ring special-case.

**`is_unique` (atom):** `AtomConstraints::add` (`atom:641`) → meet-`add`; DSL dup-check
`dsl/atom:548` → deleted; `AtomConstraint::is_unique` (`atom:116`) +
`test_atom_constraint_is_unique` deleted.

**`remove_entry` → `compare_and_set` (atom):** atom `apply_constraint` (`delta.rs:1024`) +
`apply_modify_atom_constraint` (`transact.rs:1671`) → `compare_and_set`. `remove_entry` itself
stays until bond/dative also migrate (replication), then is deleted. Replication note: the
aromatic/multicenter/stereo `apply_constraint` currently `remove_by_key` **without** verifying
the old value — `compare_and_set` fixes that.

**`AtomConstraintKey::kind()`:** delete + `test_atom_constraint_key_kind`; only tests reference it.

### Subitems (atom slice; green after each subitem unless marked)

**S1 — additive primitives on `AtomConstraints` (all green, nothing removed):** **Done**
- **S1a** `constraint/atom.rs` — `set(&mut self, c: AtomConstraint)`. `find_by_key(c.key())`:
  `Ok(i)` replace `entries[i]`, `Err(i)` insert at `i`; if `c.is_undetermined()`, remove at that
  key instead (vacuous ⇒ absent). Infallible. Tests: set-into-empty, set-overwrite-same-key,
  set-vacuous-removes, set-new-key-sorts. [dep: —]
- **S1b** `constraint/atom.rs` — `update(&mut self, other: &AtomConstraints)` = `for c in
  other.iter() { self.set(c.clone()) }`. Tests: overlay overwrites shared, keeps disjoint,
  vacuous-in-other removes. [dep: S1a]
- **S1c** `constraint/atom.rs` — `compare_and_set(&mut self, old: Option<AtomConstraint>, new:
  Option<AtomConstraint>) -> Result<(), Contradiction>`. Key from `old`/`new` (both `Some` must
  agree on `key()` → else `Err`); verify `get_by_key(key)` `canonical_eq` `old` (both `None`
  ok) → else `Err`; then `set(new)` or `remove_by_key(key)`. (Error type: reuse `Contradiction`
  for the mismatch; transact maps it to `OldStateMismatch` — open to a dedicated variant if
  that conflation grates.) Tests: verified-modify, verified-remove (`new=None`), add-from-absent
  (`old=None`), old-mismatch → `Err`, key-mismatch → `Err`. [dep: S1a]
- **S1d** `ast/atom.rs` — `AtomAst::update` body → `constraints.update(&other.constraints)`
  (replaces the `remove_by_key`+`add` loop). Tests: existing `AtomAst::update` tests unchanged.
  [dep: S1b]

**S2 — atom delta/transact onto `compare_and_set` (breaking → green):**
- **S2a** `ast/delta.rs` `AtomDelta::apply_constraint` (≈1018) → `ast.constraints.compare_and_set
  (old, new)`. `remove_entry` left in place (bond/dative still call it). Tests: reaction property
  tests over atom constraint deltas; apply/undo round-trip. [dep: S1c]
- **S2b** `ast/molecule/transact.rs` `apply_modify_atom_constraint` (1671) → `compare_and_set`
  (kind-mismatch + old-verify fold into it; `map_err` to `OldStateMismatch`). Tests: transact
  add/modify/remove atom constraint; mismatch → `OldStateMismatch`. [dep: S1c]

**S3 — meet-`add` + `is_unique` removal (breaking → green):**
- **S3a** `constraint/atom.rs` — `add(&mut self, c) -> Result<(), Contradiction>`:
  `find(c.key())` → `Ok(i)` meet `entries[i]`'s payload with `c`'s: `existing.meet(&c).ok_or
  (Contradiction)?` (`Lattice::meet` returns `Option`; `None` = incompatible → `add` maps it to
  the AST `Contradiction`); `Err(i)` insert. Delete the `is_unique`/append branch. `Canonicalize`
  drops the same-scope ring dedup. No caller reads `add`'s old `Option<replaced>` return, so the
  switch is clean — and because `Result` is `#[must_use]` (unlike `Option`), the compiler flags
  every unhandled site. Routing rule: **`add` (fallible meet) only where *user* input can
  conflict; `set` (infallible, last-wins) for all construction of known/computed values.** Caller
  migration (audited 2026-07-06):
  - **`add` (fallible) → `ParseError` — the *only* fallible site:** the parse assembly
    `ast.constraints.add(c)` (`dsl/atom.rs:553`), where user predicates accumulate (`#V4#V3`
    meets). → `add(c).map_err(|_| ParseError::ContradictoryPredicate(tag))?` (see S3b; replaces
    `DuplicateAtomPredicate`'s constraint role). Today `Contradiction` never reaches the parse
    path, so this bridge is new.
  - **`set` (infallible) — everything else:** `with_constraint`/`with_constraints`
    (`ast/atom.rs:75/89`, stay fluent); perception umol-graph `aromaticity.rs:236/301`,
    `clar.rs:152`, `hmo.rs:322`, `hueckel_rule.rs:283`, `kekulizer.rs:192`, `validate.rs:285`,
    `validate/aromaticity.rs:132`; umol-io `table_ir/raise.rs:167/179`; `from_iter`/`from`;
    `AtomView` materialization (`view/atom.rs:476–521`); `raise_atom_constraints`
    (`dsl/atom.rs:756–799`); `Lattice::meet`/`join` build (`atom.rs:815–929`, rewritten in S4a).
    All infallible — no `Result`, no `.expect()`, no fallible builders. (`ast/atom.rs:102` is the
    `update` loop S1d deletes.)
  - **tests of `add` itself → `.unwrap()`**: `molecule/tests.rs`, `dsl/atom.rs:1610–1611`.
  Tests: `add(V(4))`→`add(V(3))` = `Err`; `add(V(*))`→`add(V(3))` = `V(3)`;
  `add(R(6,1))`→`add(R(6,2))` = `Err`; disjoint scopes coexist; insert-new. [dep: S1]
- **S3b** `constraint/atom.rs` + `dsl/atom.rs` + `dsl/error.rs` — delete `AtomConstraint::is_unique`
  (116) + `test_atom_constraint_is_unique`; delete the atom DSL dup-check (`dsl/atom.rs:548`); add
  a **new `ParseError::ContradictoryPredicate(String)` variant** (parallels
  `DuplicateAtomPredicate`, carrying the tag) and map the assembly `add` at
  `dsl/atom.rs:553` with `.map_err(|_| ParseError::ContradictoryPredicate(tag))?`. This is the sole new `Contradiction →
  ParseError` bridge and is what `#V4#V3` now produces (replacing `DuplicateAtomPredicate` for the
  constraint case; the field cases keep `DuplicateAtomPredicate` until its own follow-up). Tests:
  `C#V4#V3` → the new variant (not `DuplicateAtomPredicate`); `C#V*#V3` → `#V3`; duplicate ring
  scopes → meet/contradiction. [dep: S3a]

**S4 — kill kind-addressing on `AtomConstraints` (breaking → green):**
- **S4a** `constraint/atom.rs` — internal callers → by-key: typed accessors (531–631) `get(Kind::X)`
  → `get_by_key(Key::X)`; ring accessors → direct `get_by_key` (delete `ring_memberships`/
  `ring_membership_value`/the `get_all` call); rewrite `Lattice::meet`/`join` (≈770–935) to the
  per-key pass. Still calls `*_by_key`. Tests: existing accessor + lattice + ring_count/size. [dep: S3]
- **S4b** `dsl/atom.rs` — `raise_atom_constraints` (747): per-branch `contains(kind)` →
  `contains_key(Key::X)`, ring arm stays no-op, `add(..)?`; `lower_atom_constraints` (864):
  replace with `retain(|c| !is_elidable_default(c, cfg))` + new `is_elidable_default`. Tests:
  default fill/elision round-trips. [dep: S3; S4a-independent]
- **S4c** external — `umol-graph/ops/aromaticity.rs:288` → `.aromatic_valence()`;
  `umol-graph/ops/transform/kekulizer.rs:88` → `.remove_by_key(Key::AromaticValence)`;
  `umol-py/src/constraint.rs:462–464` (tests) → `contains_key(Key::X)`. [dep: S3]
- **S4d** `constraint/atom.rs` — teardown + rename (atomic): delete by-kind `get`/`get_mut`/
  `contains`/`remove` (taking `Kind`), `get_all`, `remove_all`, `get_by_key_mut`, `iter_mut`,
  `AtomConstraintKey::kind()` (+ tests); then rename `get_by_key`→`get`, `contains_key`→
  `contains`, `remove_by_key`→`remove`, and the private `find_by_key`→`find`; mechanically
  update the S4a/S4b/S4c call sites + S1–S3 internals to the freed names. Migrate the deleted-accessor
  tests (`test_atom_constraints_get_by_key_mut`, the by-kind `get_mut` tests) to `set` +
  `get`-assert. [dep: S4a, S4b, S4c]

Critical path within atom: **S1 → {S2, S3} → S4** (S2 and S3 are independent given S1; both
precede S4). Each subitem carries its tests and ends green (S3a briefly red while `add` callers
migrate, green by S3a's end).

### After the atom slice
- **Replication** — repeat S1–S4 for `BondConstraints`, `DativeBondConstraints`,
  `AromaticSystemConstraints`, `MulticenterBondConstraints`, `NoncovalentBondConstraints`
  (near-trivial: uninhabited enum), `StereoAtomConstraints`, `StereoBondConstraints`. Each is
  independent (per-family boundary); the shared `diff_field_ops!` macro is untouched throughout.
  Their `apply_constraint` (`remove_entry`/`remove_by_key`) → `compare_and_set`; delete
  `remove_entry` once the last family drops it. DSL peers (`dsl/bond.rs`, `dsl/dative.rs`,
  `dsl/stereo.rs`) get the same `raise`/`lower`/dup-check treatment.
- **S5 — `*Ast` rename (family-wide, mechanical).** Once all families are converted:
  `AtomConstraint`→`AtomConstraintAst`, `AtomConstraints`→`AtomConstraintsAst` and peers;
  `*ConstraintKind`/`*ConstraintKey`/`*Dsl` unchanged. Not atom-scoped — renaming atom alone
  would leave the family inconsistent.
- **Scheduled: remove `ParseError::DuplicateAtomPredicate`.** *Not* constraint-specific — of its
  ~10 raisers only `dsl/atom.rs:549` is the constraint dup-check (deleted in S3b); the rest guard
  duplicate atom *fields* (`#i`/`#c`/`#h`/`#n`, spin `#u`/`#s` via `apply_spin_pair`). Full
  removal requires extending meet-on-duplicate to those fields (`#c+#c-` → contradiction,
  `#c*#c-` → join, mirroring `#V4#V3`) — a change to the field parsers + `apply_spin_pair`, not
  the constraint container. Schedule after the constraint work; do not assume S3b removes the
  variant.

## Why this doc

The Python binding (doc 137) forced the question of what a good `AtomConstraints`
interface is. `AtomConstraints` is the largest of a family — the same container shape
recurs for bond, dative-bond, and stereo constraints — and building the Python surface
surfaced a cleaner structural model than the current `kind` + `is_unique` split.

This doc is about the container's *structure and API*. What constraints *represent*
(their role in resolution and atom typing) is doc 125's subject; the Rust use-inside-
the-atom changes attach there.

## Current state

- `AtomConstraints` = `SmallVec<[AtomConstraint; 2]>`, kept sorted by `key()`, unique
  by key.
- `AtomConstraintKind` — strum discriminant, 13 variants.
- `AtomConstraintKey` — kind plus an optional `RingScope`: `RingMembership(RingScope)`;
  every other kind is keyless.
- `RingScope` = `All | Size(u8)`.
- Lookup: `get(kind)` (first match) or `get_by_key(key)`. `is_unique()` is false only
  for `RingMembership` (several per atom, one per scope).

## Settled principle: per-key constraints are unique

Every constraint is unique by *key*. The `kind → key` relation is a map with unique
keys: a kind that carries a sub-discriminator (ring scope) maps to several keys, each
unique. The `is_unique`/non-unique distinction was a distraction — it conflated "one
key" with "one kind," and only `RingMembership` broke the latter.

Multiplicity semantics live *in the value*, never in repeated keys:

- **Disjunction** is already the `LitSet` value (`{a, b}` = "a or b").
- **Conjunction**, if a constraint ever needs it (none identified so far), would be an
  explicit value form — the counterpart to `LitSet` — not repeated keys.

So the container is a pure map from key to a single value.

## The reframe (leaning): ring membership as a map-valued constraint

Take the sub-discriminator one step further and make `RingScope` the *actual map key*
inside a single `RingMembership` constraint:

- `RingMembership(BTreeMap<RingScope, ValueAst>)` — one constraint holding the per-scope
  counts.
- Then every `AtomConstraintKind` is unique *in the container*, which becomes a pure
  `kind → payload` map and sheds `is_unique`/`get_all`/`get_by_key`. `AtomConstraintKey`
  itself is **retained** at the delta/DSL boundaries — the delta needs per-scope keying
  (see the delta-boundary note under Implementation plan). The idea was implicit in
  "unique by key" all along; this names it.

Rust:

- `RingScope` already derives `Ord` and `Hash` (it is a semantic `Option<u8>`, `All`
  first), so either `BTreeMap` or `HashMap` works with no new derive. `BTreeMap` gives
  canonical iteration order for free.
- Efficiency is at least the current linear scan over the `SmallVec`; a handful of
  scopes in practice.
- The map is owned, so mutation is direct: `.insert(scope, v)`, `.remove(scope)`,
  `.entry(RingScope::Size(6)).or_insert(1)`.
- Lookup: `constraints.ring_membership[RingScope::All]`, `[RingScope::Size(6)]`.

DSL surface: unchanged. `#R2#R(6)1` still parses — each `#R…` is a scope→count entry.
A canonical ordering `All, Size(n)` (n increasing) keeps ser/de deterministic. A block
form (`#R{all=>2, 6=>1}`) is possible but the per-entry form is simpler; deferred.

Naming: the current variant is `RingScope::Size(u8)`; `Sized` was raised as an
alternative (open).

This is the lead decision. **Settled**: adopt the map model in the container; the extra
`*Dsl` boundary (below) is justified by removing the container's by-key machinery
(`is_unique`/`get_all`/`get_by_key`). `AtomConstraintKey` is kept at the delta/DSL
boundaries.

## Type mechanics

- **`BTreeMap` is forced, not chosen.** `AtomConstraint` derives `Ord` and `Hash`;
  `HashMap` implements neither, `BTreeMap<K: Ord, V: Ord>` implements both (plus
  `Clone`/`Eq`/`Default`) over its sorted entries, so the map field must be `BTreeMap`.
  Its `Eq`/`Ord`/`Hash` are order-independent, giving canonical equality/hash/order for
  free (no `Vec` insertion-order hazard). `ValueAst` already derives `Ord`+`Hash`, so
  `BTreeMap<RingScope, ValueAst>` is fully derivable.
- **`AtomConstraintKey` leaves the container, not the codebase.** The container drops the
  by-key API (`is_unique`/`get_all`/`get_by_key`/`remove_by_key`) and keys lookups on the
  discriminant alone. But `AtomConstraintKey` is **retained** for the delta and DSL
  boundaries (per-scope identity), so `constraint_key`/`fold_preserved` stay unchanged. The
  original "delete `*Key` outright" aim was wrong — see the delta-boundary note under
  Implementation plan.
- **Newtype, not bare map.** Keep a named `RingMembershipAst(BTreeMap<…>)` (methods + the
  Python mirror); a blanket `impl … for BTreeMap` would not reach a newtype, so it
  delegates. The proxied surface is small: `get`, `Index` (explicit — `[]` does not
  forward through `Deref`), `iter`, `len`, `is_empty`, `contains_key`, `insert`, `remove`,
  `entry`, `clear`, plus `Default`/`FromIterator`. `entry` exposes std's
  `btree_map::Entry<'_, RingScope, ValueAst>` (mutable-borrow-for-the-chain) — no Entry
  newtype, and no `Deref`-to-map (which would leak the full `BTreeMap` surface).
- **Traits.** `Canonicalize` — cheap (map each value + drop undetermined, propagate
  `Contradiction`); a `BTreeMap` blanket needs no trait change, though the newtype
  delegates anyway. `AsLit` — a map has no single lit; do not impl it on the map (use
  per-value `ring_membership[All].as_lit()`), sidestepping a "None if any non-lit"
  decision. `Lattice` — hand-write it: per-key `meet`/`join`, missing key =
  `Undetermined` (vacuous); `is_undetermined` = empty ∨ all-undetermined; `is_ground` =
  non-empty ∧ all-ground. A `Lattice` blanket would require adding `top()`/`const TOP` to
  the trait, touching every existing impl (`ValueAst`, valence/stereo ASTs, `ElementAst`,
  …) — a trait-wide change to save a dozen lines; not worth it unless `top()` is wanted
  for its own sake.
- **Two-step canonicalization.** (1) the map canonicalizer canonicalizes each value and
  **drops** the entries that land on `Undetermined` — the exact analog of the container
  pruning a vacuous scalar constraint (an undetermined scope-value ≡ an absent scope, no
  information lost), so a canonical map has no undetermined values and `is_undetermined ⟺
  is_empty`; (2) the container then drops the whole constraint when `is_undetermined`.
- **Python mapping is standard.** A read-only `collections.abc.Mapping` as PyO3 dunders
  (`__getitem__`/`__len__`/`__iter__`/`__contains__` + `keys`/`values`/`items`/`get`),
  constructed from a Python dict; immutable, so no `__setitem__`, and `entry` does not
  cross. For `RingScope` to be a dict key the mirror needs `#[pyclass(eq, hash)]` (it has
  neither today), else key by `int | None`.
- **Empty-map edge.** An empty map (or an all-undetermined one, which the drop step
  reduces to empty) is vacuous and dropped, exactly as a lone `Undetermined` is today.

## Serialization (the second gate)

`#R2#R(6)1` *is* a serialized map — a sequence of scope→count entries — so the data
round-trips. The friction is that multiplicity moves from the **list level** into **one
field**:

- Today N `#R` units are N constraints: a flat reader parses each independently and adds
  it; the writer iterates the list.
- Under the map model N `#R` units are one constraint's map: the writer's `RingMembership`
  case must **expand** (map → N units) and the reader must **collect** (fold repeated
  same-tag units into one map field) rather than set a value. That collect-into-a-field is
  a real break from the uniform "one tagged unit → set one thing" reader — and it is not
  merely relocated from today's non-unique `add`, because those units never had to
  converge on a single field.

**Resolution: a boundary `*Dsl` type (the pattern already used throughout).**
`AtomConstraintsDsl` keeps the flat `(kind, subkey)→value` shape, so `FromEdn`/`ToEdn`
stay uniform — every `#R…` is an independent entry, no collect-mode in the framework. The
multiplicity reconciliation lives in the boundary conversion, as plain Rust:

- `to_ast` (Dsl → Ast) **collects** the flat entries into the `BTreeMap` — one fold over
  an already-parsed list.
- `from_ast` (Ast → Dsl) **expands** the map back to flat entries in canonical order.

The reader/writer never see the map, the flat surface is preserved, and the Ast is free to
be the clean map (data + traits + `Lattice`). One subtlety in the collect: **duplicate
keys must `meet`, not last-win** — two `#R(6)` entries are `RingMembership(6)` twice, which
today's `canonicalize` merges by value-`meet` (Err on contradiction), so `to_ast` folds
with `entry().and_modify(meet).or_insert()` (or a normalize pass) to preserve that.

**The block form (`#R{all=>2, 6=>1}`) is rejected.** It is verbose even for ring
membership, and it degrades badly across the rest of the map-valued constraint family:
the pattern recurs for stereo elements (multi-key fluxionality, ligand symmetry) whose
keys are permutations and values are glyphs, where a block reads as
`#f{(1,2)(3,4) => ', (1,3) => =}` — unconventional in all the wrong ways. Keeping the flat
repeated-tag surface and reconciling at the `*Dsl` boundary is the general answer for the
whole family, not just ring membership. The second gate closes.

## Container is a heterogeneous map (key ⇒ value type)

Under the map model the key determines the value *type*: `valence → ValueAst`,
`aromatic_valence → AromaticValenceAst`, `ring_membership → RingMembership` (itself a
scope map), `tetrahedral_stereo → TetrahedralStereoAst`. The container is closer to
sparse, optional atom properties than to a uniform collection, so `constraints[key]`
returns the **payload**, not the wrapping `AtomConstraint`.

## Python surface (assumes the map model)

**Q1 — key expression.** Rust keeps the `Kind` enum (exhaustive matching is the point);
the question is the Python facade.

- Leaning: **string keys = each kind's DSL keyword** (`constraints["valence"]`,
  `"valence" in constraints`, iteration over keys). The key doubles as S6 DSL
  consistency. Typed via `.pyi` `Literal` overloads, which also pin the per-key value
  type (the heterogeneous-value problem above).
- Alternative: enum keys (`constraints[Kind.Valence]`) — type-safe but verbose and
  exposes the enum. `StrEnum` (member *is* the string) would be ideal but is 3.11+ and
  the wheel is abi3-py39.

**Q2 — ring scope selector.** With the map model `ring_membership` is itself a map, so
lookup is indexing it: `constraints["ring_membership"][RingScope.All()]` /
`[RingScope.Size(6)]`, or a friendlier `int | None` selector (`[6]` = size, `None`/`.all`
= all — conflates `None` with `All`). Open which.

**Q3 — replace by kind.** Worth having, with the *same key-spec as lookup*. `AtomAst` is
immutable, so it is an immutable update: `constraints.set("valence", 4) → AtomConstraints`;
ring membership inserts/removes into its scope map.

**Q4 — bare lits.** The raw enum variant cannot coerce — `AtomConstraint.Valence(4)`
fails because the variant field is a fixed `Py<ValueAst>`. Bare lits enter through
`ValueArg` (int | mirror) arguments:

- the set/build path — `constraints.set("valence", 4)`, and
- factory functions — `AtomConstraint.valence(4)` over `ValueArg`, mirroring Rust's
  `AtomConstraint::valence(impl Into<ValueAst>)` — alongside the raw
  `AtomConstraint.Valence(ValueAst.Lit(4))`.

## Recurrence

Bond, dative-bond, and stereo constraint containers share the shape — each its own kind
set, with ring-membership-style parameterized entries collapsing to the same map model.
Whatever is settled here templates across all four.

The map-valued case is sharper for **stereo elements** than for ring membership:
fluxionality and ligand symmetry are multi-key, keyed by *permutations* with *glyph*
values. Their flat surface (`#f(1,2)(3,4)'#f(1,3)=` style) is far more legible than a
block map, which is why the flat-notation-plus-`*Dsl`-boundary answer is the one that
generalizes — the block form is essentially unusable here.

## Decisions

- **Settled** — per-key uniqueness; multiplicity carried in the value (disjunction =
  `LitSet`; conjunction = an explicit value form if ever needed).
- **Settled** — ring membership as `BTreeMap<RingScope, ValueAst>` in the *container* (the
  map model): the container sheds `is_unique`/`get_all`/`get_by_key`/`remove_by_key`. The
  extra `*Dsl` type is justified by that container tidy-up.
- **Settled (delta boundary)** — `AtomConstraintKey` is **retained** but relocated to the
  delta (`constraint_key`) and DSL boundaries; per-scope ring edits stay independently
  keyed, so `fold_preserved`/`ConstraintKey` are unchanged. `diff`/`apply`/`from_ast`
  expand the map into single-entry ring units, `add`/`to_ast` collect — the same pattern
  as the `*Dsl`. The reaction machinery ruled out deleting `*Key`: whole-map deltas collide
  disjoint-scope edits and cannot express a rule's substructure modification against a host
  carrying other scopes.
- **Settled** — serialization stays on the flat repeated-tag surface; the flat↔map
  reconciliation lives in an `AtomConstraintsDsl` boundary (`to_ast` collects, `from_ast`
  expands, `meet` on duplicate keys). The block-map form is rejected — verbose for rings,
  unusable for the stereo map-valued constraints (permutation keys, glyph values).
  Generalizes to the whole constraint family.
- **Settled (names)** — `AtomConstraints → AtomConstraintsAst`, `AtomConstraint →
  AtomConstraintAst`; `AtomConstraintKind` keeps `Kind`, no `Ast` (a discriminant, not a
  Lattice type); `AtomConstraintsDsl` is the boundary. `RingScope::Size` (not `Sized`). A
  family-wide pass (bond / dative-bond / stereo containers the same).
- **Open (Python surface, deferrable to the binding rework)** — string keys (DSL keyword)
  vs enum keys; ring-scope selector form (`RingScope` as a dict key needing `eq`/`hash`,
  or `int | None`).
- **Open (sequencing)** — the Rust structural rework coordinates with doc 125's in-flight
  constraint work; whether it is planned from here or folded there.

## Implementation plan

**Status: atom slice complete + green** (umol-ast 4278 tests, umol-py 84 pytest). Done ahead
of the wrapper's constraint surface so the binding is built once on the clean structure.
(Distinct from doc 125, which stays parked.)

The one mechanism beyond the written plan: `diff`'s expand is a uniform `key_units()` on each
constraint container — atom splits the ring map into single-scope units, the single-key
families pass their entries through. The macro-generated `diff_constraints` keys on
`key_units()` (not `iter()`), so `fold_preserved`/`ConstraintKey` stay untouched and the
other families are unaffected. `add` meets on same-scope ring (conjunctive), storing the
empty-`LitSet` bottom on an incompatible meet so `canonicalize` surfaces the contradiction;
`apply_constraint`/`transact` drop the verified ring scope before `add` so a modify *sets*
rather than meets. The umol-py `RingMembership` mirror stays single-scope (map-shaped mirror
is D1); its atom boundary bridges single-scope ↔ single-entry map and errors on a multi-scope
map rather than dropping scopes.

**Revised architecture — the delta boundary.** The earlier plan aimed to delete `*Key`
everywhere; the reaction machinery makes that wrong. `fold_preserved` composes
`ModifyConstraint`s keyed by `ConstraintKey`, so per-scope ring edits must stay
independently keyed: a whole-map constraint keyed on bare `Kind` would (a) collide two
disjoint-scope edits into one `HashMap` slot → contradiction, and (b) be unable to express
a rule's *partial* (substructure) modification against a host that carries other scopes —
a delta is a rule, not tied to a concrete old state. So the split is removed from the
**container only**: `AtomConstraintKey` is **retained** and relocated to the two boundaries
(`delta.rs::constraint_key`, the DSL), where per-scope identity is real and reads fine.
`fold_preserved` and `ConstraintKey` are **unchanged**; the container becomes the clean
kind-keyed map; `diff`/`apply`/`from_ast` **expand** the map into per-scope units and
`add`/`to_ast` **collect** — the same expand/collect the `*Dsl` already uses. A
`RingMembership` is a multi-entry map in the container, single-entry at the boundaries (the
one invariant; the same shape as the `*Dsl` flat entries). The complexity lands on the
reaction/serialization side, which is read less often than the container.

**Method.** Additive-first via a provisional `MapRingMembershipAst` alongside the retained
`RingMembershipAst` (bond/dative stay on the old type), so the tree is green after the atom
conversion. Other families and the mechanical rename are **out of this slice**.

### S1 — atom map container (green)

- **S1a** *(done)* — `ring.rs`: `MapRingMembershipAst`, a `BTreeMap<RingScope, ValueAst>`
  newtype + proxied surface (`get`/`Index`/`iter`/`len`/`is_empty`/`contains_key`/`insert`/
  `remove`/`entry` exposing `btree_map::Entry`/`clear`/`Default`/`FromIterator`) +
  `Canonicalize` (drop-undetermined) + hand-written `Lattice` (per-scope meet/join, missing
  = `Undetermined`; `is_undetermined` = empty ∨ all-undetermined; `is_ground` = non-empty ∧
  all-ground); no `AsLit`. The original `RingMembershipAst` is **retained** for the other
  families. Tests: map surface + lattice.
- **S1b** — `constraint/atom.rs` `AtomConstraint`: `RingMembership(MapRingMembershipAst)`;
  ring branches (`is_undetermined`/`as_undetermined`/`canonicalize`) via the map;
  `ring_membership(scope, count)` builds a **single-entry** map; `key()` retained (extracts
  the single scope — a boundary method, asserts single-entry); `AtomConstraintKey`
  unchanged. *Breaking.* [dep: S1a]
- **S1c** — `AtomConstraints` container, map-based. Add `ring_membership() ->
  Option<&MapRingMembershipAst>`; `ring_count`/`ring_size_count` wrap it. `add` merges a
  ring map; `remove_entry` reimplemented kind-based (no scan); `Canonicalize` sheds the
  same-scope ring-dedup (the map's own `Canonicalize` handles it); `Lattice`
  meet/join/matches via the map. **Remove** `is_unique`, `get_all`, `remove_all`,
  `contains_entry`, `find_by_key`, `contains_key`, `get_by_key`, `get_by_key_mut`,
  `remove_by_key`; migrate `ast/atom.rs::update` (`remove_by_key(c.key())` →
  `remove(c.kind())`). *Breaking.* [dep: S1b]

- Fix AtomConstraint ring membership constructors `ring_membership` replaces the map,
  `ring_count`, `ring_size_count` set individual values.
- Check   AtomPredicate::Constraint(AtomConstraint::ring_membership(rm.scope, rm.count)) in dsl/atom.rs
-    .expect("single-entry ring membership at the DSL boundary") in EDN shape??
- AtomConstraints::add needs reworking, cannot be special-casing RingMembership like this,
  needs `update` or smth like that on the inner AST type (replace for single value, merge for map-value).
- Method name `remove_entry` is confusing, `remove_if_matches`?
- Is `binary_search_by_key` still available?
- Canonicalize: Should be sorted, why sort again?
- Add Ord to AtomConstraintKind
- Fix path imports (transact.rs)
- All the stuff is manually unrolled in transact -> unacceptable, everything is special-cased.
- test_atom_constraints_ring_membership completely fucked up
- test_atom_constraints_remove_ring_membership - against skill
- key_units -> expands into single-kv maps

### S2 — atom delta + transact (the reaction boundary)

Expand/collect keyed by the **retained** `AtomConstraintKey`; `fold_preserved` and
`ConstraintKey` untouched.

- **S2a** — `delta.rs` `AtomDelta`: `ConstraintKey = AtomConstraintKey` (**unchanged**),
  `constraint_key(c) = c.key()`. `diff` for `RingMembership` **expands** the lhs/rhs
  map-diff into per-scope single-entry `ModifyConstraint` units; `apply_constraint` for a
  single-entry ring unit modifies the host's map entry (verify `host[scope] == old`,
  set/remove). *Breaking.* [dep: S1]
- **S2b** — `transact.rs::apply_modify_atom_constraint`: for ring, compare the single-entry
  delta against the host's map entry and modify. [dep: S1, S2a]

### S3 — atom DSL boundary (flat ↔ map)

- **S3a** — `dsl/`: `AtomConstraintsDsl` (new) — a flat per-scope list; `to_ast`
  **collects** ring entries into the map (`meet` on duplicate scope), `from_ast`
  **expands** the map to per-scope entries in `All, Size(n)` order. `AtomConstraintDsl`
  handles single-entry ring via the existing flat `RingMembershipDsl` shape; wire into
  `AtomDsl`. Tests: `#R2#R(6)1` roundtrip, dup-scope `meet`, canonical order. [dep: S1]

### Out of this slice

- Bond / dative-bond / stereo: the same treatment, one family at a time.
- **Mechanical rename** `MapRingMembershipAst → RingMembershipAst` (and the `*Ast` container
  rename) once every family is converted.
- Binding (umol-py): map-shaped `RingMembership` mirror + `AtomConstraints` dict surface
  (interim binding flagged in 137).

Critical path **S1 → S2 → S3**; green after each (bond/dative remain on the retained
`RingMembershipAst`). The reaction complexity lands in S2, deliberately — the container
(read more often) stays clean.

## Cross-references

- Doc 125 — what constraints represent (projection / view model), use inside the
  resolver and atom typing. The Rust use-side changes attach there.
- Doc 137 — the Python binding that drove this; its `AtomConstraints` surface
  (enum-`kind` `get`/`contains`) is interim, pending the decisions here.
