# 138 — Constraint container API

Status: Active
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

Storage: `SmallVec<[AtomConstraint; 2]>`, **sorted by `key()`, one entry per key** — both
structural invariants, maintained by every write (`set` overwrites the same key; `extend`/
`from_iter` are `set`-loops, last-wins; `update` is a `set`/`remove` loop). No same-key duplicate
ever reaches `canonicalize`.

Read — by full key:
- `get(key) -> Option<&AtomConstraint>`
- `contains(key) -> bool`
- `iter()`
- `find(key) -> Result<usize, usize>` — **private** binary-search insertion point; the primitive
  `get`/`set`/`remove`/`compare_and_set` fold into.

Write — **verbatim** (no meet, no vacuous trick):
- `set(c)` — insert-or-replace at `c.key()`, verbatim; a vacuous `c` is *stored*, not removed.
  Infallible. (Takes a whole `AtomConstraint`, not `(key, value)`: the payload type varies by kind
  — `ValueAst` / `AromaticValenceAst` / `MulticenterValenceAst` / `TetrahedralStereoAst` — so
  `set(key, value)` would need a payload union whose tag duplicates the key.)
- `update(&other)` — overlay: for each entry in `other`, a vacuous (Undetermined) one `remove`s
  that key, otherwise `set`. Infallible. Replaces the `remove`+`add` loop in `AtomAst::update`; the
  vacuous-remove is what lets a reaction's `#v*` rhs generate a remove-constraint delta (via
  `AtomAst::update` → `AtomDelta::diff`).
- `compare_and_set(old: Option<AtomConstraint>, new: Option<AtomConstraint>) -> Result<(), Contradiction>`
  — verify the value at the key `canonical_eq` `old` (`None` = expect absent), then `set(new)` /
  `remove(key)`. The delta apply/undo primitive.

Remove:
- `remove(key) -> Option<AtomConstraint>` — unconditional delete.
- `clear()`

Bulk / construction:
- `retain(pred)`
- `new()`, `FromIterator`/`extend` (`set` in a loop — last-wins, no `debug_assert`), `IntoIterator`
  (key order).

Lattice — sorted merge delegating the shared-key case to `AtomConstraint`'s value-op:
- `meet(&self, &other) -> Option<Self>` — shared key → `AtomConstraint::meet` (`None` aborts);
  A-only / B-only kept.
- `join(&self, &other) -> Result<Self, NoJoin>` — shared → `AtomConstraint::join` (`Err` ⇒ drop
  key); A/B-only dropped. Container has ⊤ (empty) ⇒ always `Ok`.
- `matches` / `is_compatible` — same merge shape.

**No `add`** — meet lives only inside the Lattice merge. **No `is_unique`** — the parse dup-check is
`contains(c.key())` (same-scope ring is a same-key duplicate; different scopes differ).

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

### Behavior change
None at parse — `DuplicateAtomPredicate` is retained (the reset dropped meet-`add`, so `#V4#V3`
and `#V*#V3` are both rejected as duplicates as before). The only shift is internal: `set` is now
verbatim (a vacuous constraint is stored and canonicalized away lazily) rather than eagerly
removing at the write.

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

**S2 — atom delta/transact onto `compare_and_set` (breaking → green):** **Done**
- **S2a** `ast/delta.rs` `AtomDelta::apply_constraint` (≈1018) → `ast.constraints.compare_and_set
  (old, new)`. `remove_entry` left in place (bond/dative still call it). Tests: reaction property
  tests over atom constraint deltas; apply/undo round-trip. [dep: S1c]
- **S2b** *(done)* `ast/molecule/transact.rs` `apply_modify_atom_constraint` (1671): **kept**
  the `KindMismatch` pre-check (behavior-preserving, and consistent with the 6 peer families
  that still raise it — the plan's "fold it in" would have made atom alone report
  `OldStateMismatch` for a key mismatch); routed verify+apply through `compare_and_set`
  (`map_err` `Contradiction` → `OldStateMismatch`). This also fixed a latent ring bug — the old
  body `add`-ed the new scope (append), leaving a duplicate that canonicalize meets, so a legit
  `1→2` ring modify contradicted; `set` overwrites the scope. old-mismatch → `OldStateMismatch`,
  key-mismatch → `KindMismatch`. [dep: S1c]

**S3 — verbatim `set`, drop `add`, Lattice merge, `is_unique` removal (breaking → green):** **Done (2026-07-07)**
- **S3a-1** `ast/traits.rs` + `ast/error.rs` + `umol-ast-macros` + **every `Lattice` impl** — make
  `join` fallible and give `AtomConstraint` a `Lattice`. *(Foundation; the one S3 subitem that is
  **not** atom-only — the `join` signature touches the whole `Lattice` surface. Breaking → green:
  the trait signature and all impls land together.)*
  - **`NoJoin`** error (`ast/error.rs`, beside `Contradiction`): the top-less join — no least
    upper bound. Distinct from `Contradiction` (meet's unrepresented ⊥); a failed join is not a
    contradiction (both operands are individually valid, they just have no common generalization).
  - `Lattice::join(&self, other) -> Result<Self, NoJoin>` (was `-> Self`). Bounded lattices return
    `Ok(..)`; meet-semilattices (`AtomConstraint`) return `Err(NoJoin)` for cross-key operands.
    `meet` stays `Option` (⊥). `widen_with` → `Result<bool, NoJoin>`. `matches`/`is_compatible`
    unchanged (meet-derived); **no** join-side predicate — `is_compatible` exists only as a
    cheap-override hook over `meet`, and join needs none.
  - `#[derive(Lattice)]` proc-macro (`umol-ast-macros/src/lib.rs:75`): derived `join` becomes
    `Ok(Self { field: <lattice>::join(&self.field, &other.field)?, .. })` (a field's `Err(NoJoin)`
    short-circuits) — regenerates for all 9 derive users below.
  - **New `impl Lattice for AtomConstraint`** (`constraint/atom.rs`) — disjoint-union /
    meet-semilattice keyed by `AtomConstraintKey`: `is_undetermined`/`is_ground` proxy to the
    payload; `meet` same-key → payload meet, diff-key → `None`; `join` same-key →
    `Ok(payload join)`, diff-key → `Err(NoJoin)`. Removes the inherent `AtomConstraint::
    is_undetermined` and `meet`; the container's Lattice merge (S3a-2) calls `Lattice::meet`.
  - **`StereoConfigurationAst::Undetermined` stays** — it is the `Default`, the DSL `*`, and the
    kind-unknown state (default stereo.rs:279; parse dsl/stereo.rs:169/187; render 832/885;
    `is_undetermined` checks reaction.rs:2576/2716), not merely a Lattice top. Its `join` stays
    total, wrapped `Ok`; only its top-crutch duty is retired.
  - **Full `Lattice` implementor list** (each `join` → `Ok`-wrap unless noted):
    - *Direct impls (15):* `ValueAst`, `BooleanAst`, `ElementAst`, `IsotopeMassAst`,
      `ElectronCountsAst`, `NoncovalentBondKindAst`, `StereoConfigurationAst`, `AromaticValenceAst`,
      `MulticenterValenceAst`, `AtomConstraints`, `BondConstraints`, `DativeBondConstraints`,
      `MulticenterBondConstraints`, `NoncovalentBondConstraints`, `AromaticSystemConstraints`.
    - *`#[derive(Lattice)]` field-wise (9):* `AtomAst`, `BondAst`, `DativeBondAst`,
      `MulticenterBondAst`, `NoncovalentBondAst`, `AromaticSystemAst`, `SpinStateAst`; plus
      `StereoAtomAst`, `StereoBondAst` (derived inside `stereo_element!`).
    - *`macro_rules!` manual impls (6):* `TetrahedralStereoAst`, `CisTransStereoAst`
      (`stereo_site!`); `StereoAtomConstraints`, `StereoBondConstraints` (`stereo_constraint!`);
      `TopicityRelationAst`, `StereogenicityAst` (`relation_ast!`).
    - *New (1):* `AtomConstraint` — the sole `Err(NoJoin)`-returning impl (cross-key join).
  - Tests: per-type `join` tests → `Ok(..)`; `AtomConstraint` same-key `join` → `Ok(payload join)`,
    cross-key → `Err(NoJoin)`; `widen_with` tests thread the `Result`. [dep: —] (foundation)
- **S3a-2** `constraint/atom.rs` — verbatim `set`, drop `add`, Lattice as sorted merge:
  - **`set` → verbatim.** Drop S1a's vacuous-remove branch: `set(c)` is insert-or-replace at
    `find(c.key())`, *storing* a vacuous `c`. S1a's `set-vacuous-removes` test → `set-vacuous-stores`.
  - **Delete `add`** (the meet-write) entirely — no caller keeps it; meet survives only in the
    merge below. Only `set` went verbatim; `update`/`AtomAst::update` keep their vacuous-remove (a
    vacuous entry in `other` `remove`s that key), so S1b's "vacuous-in-other removes" test stays.
  - **Rewrite `Lattice::meet`/`join`** (≈895–1015) as a **two-pointer sorted merge** over the two
    key-sorted slices, delegating the shared-key case to `AtomConstraint::meet`/`join` (S3a-1); an
    A-only/B-only key follows the meet (keep) / join (drop) rule; the result is built by ordered
    `push`. Deletes the hand-written 13-arm per-kind bodies and the `!is_undetermined()` guards.
    `matches`/`is_compatible` follow the same merge. `canonicalize` = value-canonicalize each +
    drop-vacuous (no same-key dedup — the container is always sorted-unique via last-wins `set`).
  - **`find` → private**; `from_iter`/`extend` are `set`-loops (last-wins, no `debug_assert`).
  - **Caller migration** of the former `add` sites is the re-walked call-site list below — each is
    `set` (verbatim), `collect` (fresh construction), or absorbed into the merge. No container write
    but `compare_and_set` produces a `Result`, so there is nothing to thread. [dep: S3a-1, S1]
- **S3b** `constraint/atom.rs` + `dsl/atom.rs` — delete `AtomConstraint::is_unique` (116) +
  `test_atom_constraint_is_unique`. The atom DSL dup-check (`dsl/atom.rs:548`) becomes
  `if constraints.contains(c.key())` → **`DuplicateAtomPredicate`** (same-scope ring is a same-key
  duplicate; different scopes are different keys); the assembly then does `set(c)`.
  **`DuplicateAtomPredicate` stays; no `ContradictoryPredicate`** — the permissive meet-merge that
  would accept `#V3#V*` is a later addition, if ever. Tests: `C#V4#V3` → `DuplicateAtomPredicate`;
  `C#R6#R6` → `DuplicateAtomPredicate`; `C#R6#R5` accepted (distinct scopes). [dep: S3a-2]

**Reset (2026-07-07): `add`-as-meet dropped — it was a footgun.** The parallel `add`(silent meet)/
`set`(silent vacuous-remove) system was a mistake. New model:
- Container write is a single **verbatim `set`** — insert-or-replace at `c.key()`, no vacuous
  special-case, no meet. `update`/`compare_and_set`/`find`/`get`/`remove`/`clear` stay; storage is
  a key-sorted vec; Lattice `meet`/`join`/`matches`/`is_compatible` are a sorted merge delegating
  the shared-key case to `AtomConstraint`'s value-`meet`; no kind-iteration.
- **Meet is used only by the container's Lattice `meet`/`join` merge** (shared key →
  `AtomConstraint` value-`meet`). **Parse does not meet** — it rejects duplicate predicates the old
  way, **`DuplicateAtomPredicate`** (any same-key duplicate → error), restored; **no
  `ContradictoryPredicate`**. That's a strict superset: `#V3#V4` is rejected either way, and
  `#V3#V*` (which a meet would accept as `V3`) is just rejected as a duplicate. A permissive
  meet-merge at parse is a possible later addition, if ever.

So each former `add` site resolves to `set` (verbatim), `collect` (fresh construction), or the
sorted merge. Re-walking the list:

1. **`Lattice::meet`/`join` builds** — ~15, `constraint/atom.rs:895–1015` — **merge, no `add`/`set`**:
   the sorted-merge rewrite builds the result by ordered `push` during the two-pointer walk, so the
   ~15 `result.add(..)` calls are *removed* (not converted). `meet`/`join`/`matches`/`is_compatible`
   delegate the shared-key case to `AtomConstraint`'s op; the `!is_undetermined()` guards go.
2. **`extend` / `from_iter`** — `constraint/atom.rs` — **both `set` in a loop** (last-wins, no
   `debug_assert`): `for c in iter { self.set(c) }` — `from_iter` into a fresh container, `extend`
   into an existing one. Last-wins `set` keeps the container sorted-*unique*, so no dup-key ever
   reaches `canonicalize`. Callers of `extend`: `AtomAst::with_constraints` (`ast/atom.rs:83`),
   `resolve_molecule_atom` (`atom_typing.rs:70`), `stereo.rs:71`. The singular `with_constraint`
   (`ast/atom.rs:74`) → one `set`. `into_iter` iterates; `remove` compacts.
3. **`AtomView` materialization** (`derive_constraints`) — ~10, `view/atom.rs:474–527` —
   **`from_iter`/collect**: a fresh copy from concrete perceived data — distinct keys, values
   never vacuous, so `add`/`set`/`collect` all coincide (plain insert). It's a can't-fail
   construction → `.collect()` (infallible). Minimal treatment — `derive_constraints` is itself
   slated for removal, so no further investment here.
4. **DSL `raise`/`lower` defaults** (`raise_atom_constraints`/`lower_atom_constraints`) —
   `dsl/atom.rs:747–…` — **`set`, restructured** (raise/lower are too complicated). Supersedes the
   S4b `add(..)?` and the `retain(...)` lowering below.
   - **Drop the global `retain(|c| !c.is_undetermined())`.** It strips *every* vacuous entry;
     we must only touch the kinds being defaulted — vacuous values of other kinds stay.
   - **No `for kind in Kind::iter()`.** One **explicit clause per relevant kind** — valence,
     donated, accepted, aromatic, multicenter, tetrahedral — written in **ascending key-sort
     order**. (These clauses will most likely go away later too.)
   - Each raise clause: if the key is **absent *or* undetermined** (vacuous ≈ absent for
     defaulting — confirmed ok), `set(default)`. Concrete user values are left untouched (the
     guard skips them); the default overwrites a vacuous entry of that same kind.
   - **`lower`:** likewise an explicit clause per kind, iterating in **reverse key-sort order**.
   - Any shared per-kind helper stays a **free fn in `dsl/atom.rs`** (raise: absent-or-undetermined
     → `set(default)`; lower: equals-default → `remove`), **not** a method on `AtomConstraints` —
     defaulting/elision is DSL-boundary logic and would pollute the container's primitive surface.
     Or just inline the ~6 guarded `set`s, since each is a one-liner over `get`/`set`.
5. **DSL parse accumulation** — 1, `dsl/atom.rs:548–553` — **`set` + `DuplicateAtomPredicate`**.
   Keep the dup-check (any same-key duplicate → `DuplicateAtomPredicate`); otherwise `set(c)`
   verbatim. No meet, no `ContradictoryPredicate`. (`#V3#V*` stays rejected as a duplicate; the
   permissive meet-merge is a later addition, if ever.)
6. **perception** — **1 production site + 5 test fixtures** — **`set`**. Only `aromaticity.rs:236`
   (`equalize_charges`, charge equalization for symmetric monoelement aromatic ions — Cp⁻,
   tropylium, `[S₄]²⁺`) writes a constraint in production; it pins `AromaticValence(Aromatic(Lit
   (k)))` alongside the `system.electrons` update on line 240 → a computed write, `set`. The other
   five `.add` sites are inside `mod tests` (aromaticity:301, clar:152, kekulizer:192, validate:285,
   validate/aromaticity:132) — fresh fixture builders, `set`/`add` coincide. (My earlier
   "perception overwrites, is_unique-replace policy" framing was wrong: `e`/`k` there are the
   system's electron counts, not the atom's constraint.)
7. **reaction partial rebuild** — 2, `dsl/reaction.rs:2144–2153` (the atom `ModifyConstraint` arm of
   `render_deltas`) — **`set` in both branches**: `Some(c)` → `set(c)`; `None` →
   `set(old.as_undetermined())`. Verbatim `set` *stores* the vacuous removal marker (renders `#v*`)
   rather than dropping it, which is what removes the earlier "set drops the marker" objection. The
   same block repeats for 6 peer families (bond 2191/2195, dative 2270/2274, aromatic 2344/2348,
   multicenter 2418/2422, stereo-atom 2560/2564, stereo-bond 2703/2707) — outside the atom slice
   (replication), keep `add`.
8. **IO raise** — 3, `umol-io/table_ir/raise.rs:67/167/179` — **`set`**. Fresh `AtomAst`
   construction from a table atom (distinct keys — tetrahedral stereo at 67, aromatic valence at
   167/179). Verbatim insert.
9. **molecule `inline_constraints`** — `molecule.rs:788` (atom arm) — **`set`**. Moves a
   molecule-level `Constraint::Atom` into the atom's inline store; the doc comment already says
   "last-wins per kind" — that's `set`. (Bond/dative arms are peers → keep `add`.)
10. **binding** — `umol-py/constraint.rs:314–319` — **`extend` / set-loop**. `for entry in entries {
    add }` building a fresh container is the set-loop; use `extend`/`collect` (last-wins per key),
    matching its docstring.
11. **tests + `replace` import** — `constraint/atom.rs` — test `add` sites → `set`; remove
    `use std::mem::replace` (only the old `add` used it).

**S4 — kill kind-addressing on `AtomConstraints` (breaking → green):** **Done (2026-07-07)** —
pure by-kind→by-key migration + teardown + rename, zero behavior change. **Deferred:** the item-4
`raise`/`lower` *restructure* (drop the global `retain`, explicit per-kind clauses in key-sort
order, `absent-or-undetermined → set`, lower in reverse) — a behavior change, separable from
killing by-kind addressing; S4b did only the minimal `contains(kind)`→`contains(Key::X)` swap and
kept the loop + global retain.
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
precede S4); within S3, **S3a-1 → S3a-2 → S3b**. Each subitem carries its tests and ends green
(S3a-1 red only while the `join` signature and its impls migrate together; S3a-2 briefly red while
`add` callers migrate, green by its end). S3a-1 is the only workspace-wide subitem — every other
S1–S4 item is `AtomConstraints`-scoped.

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
- **`ParseError::DuplicateAtomPredicate` stays** (the reset dropped `ContradictoryPredicate`).
  S3b keeps it for the constraint dup-check (`contains(c.key())` → duplicate), and its other ~10
  raisers guard duplicate atom *fields* (`#i`/`#c`/`#h`/`#n`, spin `#u`/`#s` via `apply_spin_pair`)
  unchanged. A permissive meet-on-duplicate (`#V3#V*` → `V3`, `#c+#c-` → contradiction) is a
  possible later addition; not scheduled.
- **Scheduled: remove `TransactionError::KindMismatch`.** Preserved through the atom slice (S2b)
  for behavior parity; used by all 7 `apply_modify_*_constraint`. Once every family routes
  verify+apply through `compare_and_set` (replication), drop the `KindMismatch` pre-checks and
  the variant — a key-mismatched modify then reports `OldStateMismatch` (or one merged error).
  Cross-cutting, so after replication, not in the atom slice.

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
