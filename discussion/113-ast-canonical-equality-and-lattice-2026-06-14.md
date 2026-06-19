# 113 · AST canonical equality and lattice algebra

Status: Active · 2026-06-14
Resolves: 095 open questions 1–3 (equality model, leak inventory, canonical form).

## General

- AST types are **open, transparent data carriers** (public variants, raw
  construction; no facade) — so canonical-by-construction is unenforceable and
  rejected; smart constructors are convenience.
- **`Canonicalize`**: `fn canonicalize(self) -> Result<Self, Contradiction>`,
  by-value, idempotent; the per-type step-lists (detail) are the normal form.
  `Err(Contradiction)` = ⊥ (same channel as `meet`'s `None`; no `⊥` variant). Pure
  field-wise types derive it; folding/validating leaves hand-write it. By-reference
  twin `fn canonical(&self) -> Result<Cow<'_, Self>, Contradiction>` (default = clone +
  `canonicalize`; overridden to **borrow** already-canonical `Lit`/`Undetermined`);
  `equiv`/`==`/`Hash`/`meet` build on `canonical`.
- Canonicalization is lazy: building or parsing a value does not normalize it — an AST node
  is exactly what was constructed.
- A value's normal form is computed only when it is compared, hashed, or combined. `equiv`
  is that comparison — it canonicalizes both sides and checks they match, so two spellings
  of the same value (`{0,1}` vs `{1,0}`, `#h2` vs `#h(2)`) count as equal — and `==`, `Hash`,
  and `Ord` use the same canonicalizing comparison. The dominant `Lit` and `Undetermined`
  values are already canonical, so they take a fast path that does no extra work.
- `Lattice` only on value/constraint atoms and index-free entities.
- **`matches` from `meet`** (default `matches(t) == (meet(t) == Some(canon t))`).
- **`join` stays `Self`** — every retained lattice has a faithful top.
- `meet`/`join` canonicalize (fast-path trivial; canonicalize complex operands and
  output); `meet: Option`, `join: Self`.
- **`PartialEq`/`Eq`/`Hash`/`Ord`** canonicalize, use fast-path for Lit/Undetermined
  variants

## Scope

Foundational: the new `Canonicalize` trait + derive, the `Lattice` trait
(`Lattice: Canonicalize`, default `matches`, `meet`/`join` = field-wise +
`canonicalize`) and its derive macro, the set-variant storage, the **structural**
`Eq`/`Hash`/`Ord` derives plus `equiv` / `Canonical<T>`, the (C) changes (`electrons`
→ `ElectronCountsAst`; `kind` → `StereoKindAst`; `pair` keeps `TopicityAst`
`matches`-only), the naming/boundary
cleanup (literal renames + `*Dsl` boundary types replacing the inline stereo-perm
serde; `perm` → `permutation`), and the construction sites that build sets. 095's
leak inventory is the re-audit surface: `MoleculeAst::PartialEq` (proptest
round-trip), `HashMap`/`HashSet` keys, constraint dedup, alias bijectivity
(`BiBTreeMap`). Expect an extended red period.

Out of scope here, **relocated atom → molecule** (not removed): dropping atom-level
`Bind`/`Ref` moves variable scope, cross-entity variable references (atom↔atom relations
such as neighbour-charge sums), and the storage of molecule-level joint/relational
constraints to the molecule level. **`JointDomain` is removed entirely** (it was a
field-tuple constraint masquerading as variables — a premature, parallel system; docs
097/098 obsoleted). Its capability (e.g. the Fe²⁺ high/low-spin split) is subsumed by the
general molecule-level **variable-constraint facility** (doc 115: typed vars + arithmetic
+ tuple/table constraints), designed later. `Var(String)` is a bare molecule-scoped
occurrence in the meantime.


## Implementation Plan

Built **bottom-up**; one extended red period (atomic in-place, no shims — broken tree
mid-refactor is acceptable). The per-type *designs* are in the detail + review tables
below; this section sequences the *work*. (`JointDomain` is already removed; `AtomAst` is
already field-wise.)

**What lazy buys us.** `Eq`/`Hash`/`Ord` stay **derived-structural** (no change) — so the
C4e.5 `lattice::` sweep (which compares with `==`) greens for a type the moment its
`meet`/`join` produce *canonical* output. The variant redesign, `Canonicalize`, and
`meet`/`join`-canonicalize therefore land **together per type** (the redesign is what
makes canonicalization clean). The broader semantic-equality adoption (`equiv` /
`Canonical<T>`, the 095 leak sites) and boundary cleanup come after — they don't gate the
sweep.

**Per-type recipe** (P1–P5): retype variants (detail) → impl `Canonicalize` (the step
list) → `canonical()` fast-path (borrow `Undetermined`/`Lit`, see the fast-path table) →
make `meet`/`join` canonicalize their output → faithful parser + canonical render
(`*Dsl`) → update the type's unit tests **and its `lattice::` proptest strategy**. Each
type greens its own `lattice::` test.

### P0 · Foundation (additive — compiles, no behavior change) **Done**
- `traits.rs`: `Canonicalize: Clone + PartialEq` (`canonicalize(self) -> Result<Self,
  Contradiction>`, `canonical(&self) -> Result<Cow<'_,Self>, Contradiction>` default
  clone+canonicalize, `equiv` default); `Canonical<T>` newtype (`new` canonicalizes once;
  derived structural `Eq`/`Hash`/`Ord`). Re-export from `ast.rs`.
- `umol-ast-macros`: add `#[derive(Canonicalize)]` (field-wise). **Do not** yet change
  `#[derive(Lattice)]` (that's P6) — keep P0 additive.

### P1 · Leaf value types **Done**
1. **`ValueAst`** **Done** - split `Expr`→`ValueTerm`(i64)/`ValuePredicate`(bool); `Bind`/`Ref`→
   `Var`; `Set`→`LitSet(BTreeSet)`; the full canonicalize step list (term fold, predicate
   NNF, lifts). Parser: build `Sum`/`Product`/`Div`/`Rem`, drop the `unary_expr` sign-XOR;
   delete `ValueAst::simplify`/`ValueExpr::simplify`. The biggest single change (all of
   `dsl/value.rs` + every `ValueAst` consumer).
2. **`ElementAst`** **Done** — `Undetermined|Lit|LitSet|NotSet|Var`; cardinality-canonical; drop `Not`.
3. **`IsotopeMassAst`** **Done** — `Undetermined|Natural|Lit|LitSet|Var`; positive-only; `u32`.
4. **`NoncovalentBondKindAst`**, **`StereoKindAst`**  **Done** — `Undetermined|Lit`; identity
   `canonicalize`; `canonical` always borrows.
5. **`ElectronCountsAst`** (new) **Done** — `Undetermined|Lit(Vec<u8>)`; identity (positional).

### P2 · Predicate / relation types **Done**
1. **`SpinStateAst`** **Done** — pure field-wise `#[derive(Canonicalize)]`; **no** cross-field
   parity gate (`unpaired`↔`multiplicity` parity is a tier-2 physical invariant, enforced
   at resolution, not the AST — invalid pairs like `(1,1)` are legal AST states); `From`
   stays infallible.
2. **`AromaticValenceAst`/`MulticenterValenceAst`** **Done** — delegate inner `ValueAst`; drop the
   hand `matches`.
3. **`TopicityRelationAst`/`StereogenicityAst`** **Done** — `relation_ast!` over `BTreeSet`; flatten
   stereogenicity (drop `StereogenicityRelationAst` + the wrapper).
4. **`TopicityAst`** — drop `impl Lattice` → matches-only `{pair, rel}`; remove its
   fixed-pair lattice proptest (covered by `test_topicity_relation_ast_lattice_laws`).

### P3 · Stereo configuration cluster **Done**

Restart from the AST and reason outward. Two facts drive the shape: a coset index
is meaningless without a `StereoKind`, and the kind reaches the two sides
differently — the **element** carries it as data (a pattern may leave the geometry
open), the **constraint** is split into one type *per kind* so the kind is the type's
identity. Canonicalization is **lazy**, so the DSL only places kind + coset next to
each other; no folding on the parse side.

Names and `AsLit` targets are settled: element type stays `StereoConfigurationAst` with
variants `Undetermined | Kind(StereoKind, StereoCosetAst)`; `AsLit` per the AST trait matrix.

#### P3a · Sub-part AST leaves (kind-relative — **no** `Canonicalize`/`Lattice`) **Done**
1. `StereoTerm = Var(Box<(String, Option<BTreeSet<u32>>)>) | Lit(u32) | LitSet(BTreeSet<u32>)
   | Swap(Box) | Mirror(Box) | Apply(Box, Permutation)`. Structural `Clone/Eq/Hash/Ord`
   only. **No `AsLit`** — a term is never a bare literal (even `Lit`/`LitSet` bases are
   pre-fold operands). (`Lit`/`LitSet` operands are approved: operator-on-literal parses
   faithfully and **folds at canonicalize**, not on parse; a canonical term is `Var`-rooted.)
2. `StereoCosetAst = Undetermined | Lit(u32) | LitSet(BTreeSet) | NotSet(BTreeSet) |
   Term(Box<StereoTerm>)`. Structural, **plus `impl AsLit<Lit = u32>`** (`Lit(i)`→`Some(i)`,
   else `None`) — extracting the literal index is kind-independent, so it lives on the
   sub-part; still **no `Canonicalize`/`Lattice`** (those need the kind). Replaces the
   `coset_as_lit` free fn.
3. Kind-aware coset algebra — free fns, each takes `kind: StereoKind`: `canon_coset`,
   `coset_meet`, `coset_join`, `coset_matches`, `coset_apply_permutation`
   (+ private `coset_to_set`/`coset_from_set`/`compose_term`). `NotSet` cardinality
   polarity per `ElementAst`. (Already implemented; keep.) Owners pass their kind — the
   element from its variant, each constraint type from its constant.

#### P3b · Element-side config — **remove `StereoKindAst`** **Done**
4. Delete `StereoKindAst` (only the config field + its own tests + the re-export use it).
5. `enum StereoConfigurationAst { Undetermined | Kind(StereoKind, StereoCosetAst) }`. The
   top-level `Undetermined` has **no coset** (kind unknown ⇒ coset meaningless); `Kind`
   binds a concrete geometry to a coset that may still be open. `*` (Undetermined) and
   `Th*` (`Kind(Tetrahedral, Undetermined)`) are **distinct**.
6. `StereoAtomAst`/`StereoBondAst` carry `configuration: StereoConfigurationAst`.

#### P3c · Constraint-side per-kind types — **replace `StereoSiteAst`** **Done**
7. `enum TetrahedralStereoAst { Undetermined | NotStereo | Stereo(StereoCosetAst) }`
   (kind ≡ `Tetrahedral` constant) and `enum CisTransStereoAst { … }` (kind ≡ `CisTrans`).
   No kind field — the type *is* the kind. (Tetrahedral/CisTrans have coset count 2, so
   their cosets are only `Undetermined`/`Lit`/`Term`; `LitSet`/`NotSet` never arise here.)
8. `AtomConstraint::TetrahedralStereo(TetrahedralStereoAst)`,
   `BondConstraint::CisTransStereo(CisTransStereoAst)`.

#### AST trait matrix
| Type | `Canonicalize` | `Lattice` | `AsLit` |
|---|---|---|---|
| `StereoKind` (plain enum) | — | — | — |
| `StereoTerm` | — | — | — |
| `StereoCosetAst` | — | — | `u32` (literal coset index; replaces `coset_as_lit`) |
| `StereoConfigurationAst` | ✓ `Undetermined`→id; `Kind(k,c)`→`Kind(k, canon_coset(c,k)?)` (never collapses `Kind`→`Undetermined`) | ✓ meet `(Und,x)`→`x`, `(Kind(k1,a),Kind(k2,b))`→ `k1≠k2?None:Kind(k1,coset_meet(a,b,k1)?)`; join same-kind→`Kind(k,coset_join)`, cross-kind→`Undetermined`; `is_undetermined`=top only (`Th*` is *not*); `is_ground`=`Kind(k,Lit)` | `StereoConfiguration{kind,coset}` (`Kind(k,Lit(i))`→`Some`) |
| `TetrahedralStereoAst` / `CisTransStereoAst` | ✓ `Stereo(c)`→`Stereo(canon_coset(c, KIND)?)`; others id | ✓ `Undetermined` wildcard; `NotStereo∧NotStereo`; `NotStereo∧Stereo`→`None`(meet)/`Undetermined`(join); `Stereo(a)∧Stereo(b)`→`Stereo(coset_meet/join(a,b,KIND))` | `StereoConfiguration{kind,coset}` (kind = type constant; `Stereo(Lit(i))`→`Some`, `NotStereo`/`Undetermined`→`None`) |

#### P3d · AST tests + property tests
9. Per-method `#[rstest]`/`#[case]` units for `StereoConfigurationAst`,
   `TetrahedralStereoAst`, `CisTransStereoAst`: `canonicalize`/`_identity`, `meet`, `join`,
   `matches`, `as_lit`, `is_undetermined`/`is_ground` (incl. the `*` vs `Th*` distinction).
   Fold the deleted `StereoKindAst` coverage into `StereoConfigurationAst`. Free-fn tests
   for the coset algebra (`canon_coset` incl. `NotSet` polarity + `Term` compose→priority
   `Mirror>Swap>Apply`, `coset_meet`/`coset_join`/`coset_matches`). **Done**
10. Property (lattice-law) tests for the three trait-bearing AST types under
    `--features proptest` (meet idempotent/commutative/associative, absorption,
    `canonicalize` idempotence).

#### P3e · DSL — element side **Done**
11. `StereoAtomDsl(StereoAtomAst)` / `StereoBondDsl(StereoBondAst)` — `Display`+`FromStr`
    (compact string), `FromEdn`+`ToEdn` (string/keyword in molecule EDN), `IntoAst`+`FromAst`.
12. `StereoConfigurationDsl(StereoConfigurationAst)` — `IntoAst`+`FromAst` only; **no
    `FromEdn`/`ToEdn`** (element side has no standalone EDN). String parse `*`→`Undetermined`,
    `<glyph><coset>`→`Kind`; render `Undetermined`→`*`, `Kind(k,c)`→glyph + coset.
13. Coset string round-trip via free fns `parse_stereo_coset(s, degree)` / `fmt_stereo_coset`.
    Coset EDN keeps the `StereoCosetDsl(StereoCosetAst)` boundary type
    (`FromEdn`/`ToEdn`/`IntoAst`/`FromAst`); degree fixed at 4 for `#T`/`#C`. The per-kind
    `TetrahedralStereoDsl`/`CisTransStereoDsl` and the streaming `read_stereo_coset_dsl` parse
    the `{:stereo <coset>}` value through it, then `into_ast` → `StereoCosetAst` wrapped as
    `Stereo(coset)`. No `StereoTermDsl` (term is inline in the coset grammar). `NotSet`↔`!{0,1}`.

#### P3f · DSL — constraint side **Done**
14. `TetrahedralStereoDsl(TetrahedralStereoAst)` / `CisTransStereoDsl(CisTransStereoAst)` —
    `FromEdn`+`ToEdn`, `IntoAst`+`FromAst`. The kind-free `FromEdn(edn)` signature **works**
    because the kind is the type's constant; parse from atom/bond DSL (`#T`/`#C` tag) and
    from EDN (`:tetrahedral-stereo`/`:cis-trans-stereo`). Render value
    `:undetermined`/`:not-stereo`/`{:stereo <coset>}`; kind emitted by the constraint key.
15. Remove the unauthorized `stereo_site_from_edn` fn and the `StereoSiteDsl` rename —
    superseded by the per-kind Dsl types' `FromEdn`.

#### DSL trait matrix
| Type | `FromEdn` | `ToEdn` | `IntoAst` | `FromAst` | string |
|---|---|---|---|---|---|
| `StereoAtomDsl` / `StereoBondDsl` | ✓ | ✓ | ✓ | ✓ | `Display`+`FromStr` |
| `TetrahedralStereoDsl` / `CisTransStereoDsl` | ✓ | ✓ | ✓ | ✓ | inline via `#T`/`#C` |

#### P3g · Consumers **Done**
16. Views `kind()`/`coset()` and `molecule.rs`/`transact.rs`/`symmetry.rs` field accesses
    (`.configuration.coset`) destructure the enum; the validator's `view.coset().as_lit()`
    resolves via `AsLit for StereoCosetAst` (drop the `coset_as_lit` free fn); `ast.rs`
    re-exports drop `StereoKindAst`
    /`StereoSiteAst`, add the enum `StereoConfigurationAst`, `TetrahedralStereoAst`,
    `CisTransStereoAst`, and the `*Dsl` types.

### P4 · Entities (pure field-wise)
1. `AtomAst`
2. `BondAst`
3. `NoncovalentBondAst` — `#[derive(Lattice, Canonicalize)]`; delete `simplify_values`.
4. Fix `NoncovalentBondAst` direct `Display`/`FromStr`/`FromEdn` (move to
  `NoncovalentBondDsl`).
5. `AromaticSystemAst`/`MulticenterBondAst` — `electrons: Vec<ValueAst>` → `ElectronCountsAst`
   (whole-vector `Undetermined | Lit(Vec<i64>)`). This narrows the field: per-cell partials,
   sets, terms, and mixed-undetermined are no longer representable. EDN: `:electrons [i j k]`
   → `Lit`, `:electrons :undetermined` (or the key omitted) → `Undetermined`. The `Canonicalize`
   derive on these two entities still waits on P5 (they hold constraint collections); this entry
   is the field retype only, keeping the existing hand `meet`/`join` and routing electrons through
   `ElectronCountsAst`.
   a. **AST (both types)** — retype the field; ctors become `new(ElectronCountsAst)`,
      `from_counts(Vec<i64>)` (→ `Lit`), `with_electrons(impl Into<ElectronCountsAst>)`
      (drop `from_electrons(Vec<u8>)` / `new(Vec<ValueAst>)` / `with_electrons(Vec<ValueAst>)`);
      hand `meet`/`join`/`is_ground`/`into_ground` route electrons via `ElectronCountsAst`
      `Lattice`/`AsLit` (whole-vector, not per-cell); the From-table lift builds `Lit` when every
      source cell is known else `Undetermined`; drop the per-cell electron step in
      `simplify_values`; retype the unit tests.
   b. **Views** — `AromaticSystemView::electrons` / `MulticenterBondView::electrons` return
      `&ElectronCountsAst` (not `&[ValueAst]`).
   c. **transact / edit** — `AromaticSystemFieldChange::Electrons` /
      `MulticenterBondFieldChange::Electrons` payload `Vec<ValueAst>` → `ElectronCountsAst`; the
      apply stays a whole-value set.
   d. **DSL** (`dsl/aromatic.rs`, `dsl/multicenter.rs`, `dsl/molecule.rs`) — `:electrons`
      parses/renders `[i j k]` → `Lit` and `:undetermined` → `Undetermined`; remove the per-cell
      set/term/undetermined grammar. Update `dsl/molecule/tests.rs`: drop `[[1 2] 1]`,
      `["?n + 1" 1]`, `[:undetermined :undetermined]`; add an `:electrons :undetermined` case.
   e. **umol-graph** — `ops/aromaticity.rs` aromatizer builds the whole `Vec<i64>` and assigns
      `electrons = Lit(v)` once (replacing the per-cell `electrons[i] = …` write); reads go through
      `as_lit`. `ops/validator/entity.rs` takes the length via `as_lit().map(Vec::len)`.
   f. **umol-io** — `table_ir/raise.rs` constructs `ElectronCountsAst::Undetermined` (replacing
      `vec![ValueAst::Undetermined; n]`) for aromatic systems and multicenter bonds.
   g. **doc comments** in `constraint/aromatic.rs` / `constraint/multicenter.rs` referencing the
      `electrons` field type.
6. `DativeBondAst` — derive after the birelation `acceptor_slot` drop.
7. Fix writing of stereo atom / bond constraints to use streaming deserializer instead of the tree-based impl.

### P5 · Constraint types — impl plan (bottom-up)

Design, lattice contracts, and notation: see *Per-type detail → Constraint types*. This
is the build order.

1. **Ring membership (model B — single-entry `RingMembershipAst` carrier, WET).** Replaces
   `RingCount` + `RingSize` on **three** enums (`AtomConstraint`, `BondConstraint`,
   `DativeBondConstraint`; atom-only `RingDegree`/`RingValence` stay). Each ring fact is **one
   constraint** carrying a `RingMembershipAst { scope: RingScope, count: ValueAst }` — a
   single-entry struct parallel to `TopicityAst { pair, rel }` (**not** a sub-map, **not** a
   shared keyed container). Several ring facts (one per `RingScope`) sit on an entity; each of
   the three collections hand-writes the per-scope handling inline (the former `RingSize` loop,
   now grouped by `RingScope`). **Done**

   a. **`RingScope` + `RingMembershipAst` + variant.** `constraint/ring.rs` (shared foundation
      module for the three sibling enums) holds `RingScope { All, Size(u8) }` (`Ord`: `All`
      first, then sizes ascending) **and** `RingMembershipAst { scope, count }` (pub fields +
      `new` ctor; carrier only — no `Lattice`, mirroring `TopicityAst`). On each enum replace
      `{RingCount(ValueAst), RingSize(ValueAst)}` with one variant
      `RingMembership(RingMembershipAst)` + a `ring_membership(scope, count)` ctor. `kind()`
      collapses to one `RingMembership` discriminant. `is_undetermined`/`is_ground`/`simplify`
      delegate to the inner `count`. `is_unique` returns **false** for `RingMembership` (the
      former `RingSize` non-unique case, renamed) — ring entries append.
   b. **Atom collection.** `AtomConstraints` (hand-written `SmallVec`): `add` **appends** ring
      entries (non-unique); per-scope dedup is **lazy** (at `meet`/`canonicalize`), consistent
      with the decided `add` policy. Accessors `ring_count()` (scope `All`), `ring_size_count(s)`
      (scope `Size(s)`), private `ring_memberships()` / `ring_membership_value(scope)`; the old
      `ring_count()`/`ring_sizes()` are gone. `meet`/`join`/`matches` run the per-scope product
      over `RingScope` (disjoint literal sizes ⇒ no subsumption, box-hull `join`), replacing the
      old single-`RingCount` path + `RingSize` loop. Tests.
   c. **Bond collection** — same for `BondConstraint` / `BondConstraints`.
   d. **Dative collection** — same for `DativeBondConstraint` / `DativeBondConstraints`. The
      `DativeBondConstraintDsl` gains the `RingMembership(RingMembershipAst)` variant; the DSL
      `apply_predicates` dup-check is gated on `is_unique` (matching bond), so multiple `#R`
      scopes on one dative bond are allowed.
   e. **DSL.** *String:* one glyph `#R<count>` (scope `All`) / `#R(s)<count>` (scope `Size(s)`),
      parsed by `predicates::ring_membership` → `RingMembershipAst`, rendered by
      `fmt_ring_membership`; `#r` is dropped. Count sugar: bare→`Lit(1)`, `+`→`var_at_least("r",
      1)` (≥1; fixed name `"r"` ok — nothing unifies yet, doc 115), `!`→`Lit(0)`, `n`→`Lit`,
      `{a,b}`→`LitSet`, `?v`→`Var`, `*`→`Undetermined`; canonical render `#R!` (0) / bare `#R`
      (1). *EDN:* a `RingMembershipDsl(RingMembershipAst)` boundary type (`FromEdn`/`ToEdn`)
      emitting `{:ring-membership {:size <int> :count <value>}}` (sized) /
      `{:ring-membership {:count <value>}}` (total), plus the streaming
      `read_ring_membership_dsl` in `dsl/constraint.rs` — the aromatic-valence pattern (a `*Dsl`
      type for the value path + a free streaming reader).
   f. **Consumers + tests** — `transact.rs` ring tests (ring is non-unique → the `Add`/`Remove`
      multi-valued path, not `Set`); molecule/dsl/property test suites; proptest generators and
      AST fuzzer seeds. Geometric ring refs (`ast/ring.rs` perception, umol-graph aromaticity,
      `views::{atom,bond}::ring_count`/`ring_size`) are graph-derived and **unchanged**. The
      `All = Σ_s Size(s)` cross-check is a **tier-2 validator** concern, **not** part of this
      AST/DSL work (deferred).
   g. Update umol-dsl-spec.md with the new syntax. **Done**
   h. `Canonicalize` (the item-3 requirement, applied to ring).** P5.1(b) already
      references per-scope dedup "at `meet`/`canonicalize`", but **no `Canonicalize` was added**:
      neither the per-enum `Canonicalize` (delegate to the inner `count`, replacing `simplify`)
      nor `impl Canonicalize for {Atom,Bond,Dative}Constraints` (group by `RingScope`, merge by
      value-`meet` → `Err` on contradiction, drop vacuous, order by scope). `meet` dedups lazily
      but there is **no container canonical form**, so equality/hashing see un-normalized
      collections. P5.1 is **not complete** until these land (identical shape to P5.2(c)); this
      also blocks P6 (`Lattice: Canonicalize`).
   i. `key()` infra.** The accepted design (a **dedicated key enum** + a **by-key
      API mirroring the by-kind one**) was conflated with "WET" and never landed — no `key()`
      method exists on any constraint. Each family needs a `<Enum>ConstraintKey` (kind +
      sub-key: `RingMembership(RingScope)`, every other kind unit) + `fn key(&self)` + a by-key
      collection API (`contains_key` / `get` / `remove` / update, mirroring the by-kind
      methods). "WET" means hand-written **per family**, **not** the absence of `key()` and
      **not** a generic shared keyed container. Same infra as P5.2(a); do together.
2. **Stereo constraints — `#g`/`#o`/`#f`/`#p`.** In `StereoAtomConstraint` /
   `StereoBondConstraint` (WET, macro-generated for atom + bond), several entries share a kind,
   distinguished by a sub-key; each entry is a `(key, value)` pair. **AST (a–c) Done 2026-06-19;
   DSL (d–e) pending — see "Stereo constraint EDN" below.**

   | glyph | kind | unique | sub-key | value | merge |
   |---|---|---|---|---|---|
   | `#g` | Stereogenicity | yes | — (unit) | `StereogenicityAst` | finite subset-lattice |
   | `#o` | Topicity | no | `StereoLigandPair` | `TopicityRelationAst` | finite subset-lattice |
   | `#p` | LigandSymmetry | no | `OrientedLigandPermutation` | `MemOp` | equal→keep, differ→`⊥` |
   | `#f` | Fluxionality | no | `LigandPermutation` | — (unit) | dedup (set presence) |

   a. **`key()` + key enum. Done.** Macro-generated **per collection**
      (`StereoAtomConstraintKey` / `StereoBondConstraintKey` — not one shared key, so `key().kind()`
      returns the collection's kind) with variants `{ LigandSymmetry(OrientedLigandPermutation),
      Fluxionality(LigandPermutation), Topicity(StereoLigandPair), Stereogenicity }` + `fn key()` +
      by-key API (`contains_key`/`get_by_key`/`get_by_key_mut`/`remove_by_key`).
      `OrientedLigandPermutation` derives `Ord` (added `Ord` to `umol_perm::Orientation`) so the key
      is orderable and the store key-sorted (binary-search by-key, parity with the other six).
      `mem` is the LigandSymmetry **value**, not part of its key — `(P,In)`/`(P,NotIn)` collide on
      key `P` and contradict.
   b. **Two-regime `Lattice`. Done.** `meet` rewritten so `#p` merges per perm (conflicting `mem`
      → `⊥`); `#o` per-pair `rel.meet`, `#g` `meet`, `#f` union. `join`/`matches` unchanged (already
      correct under mem-in-value). **No `MemOp` `Lattice` impl** — the 2-element merge is inline
      (equal→keep, differ→`⊥`), keeping `MemOp` a plain operator type.
   c. **`Canonicalize` (trait). Done.** `impl Canonicalize` for the enum (canonicalize the inner
      relation; `#f`/`#p` atomic) and the container (sort by `key()`, merge same-key by
      value-`meet` → `Err` on `⊥`, drop vacuous). `add` sort-inserts; `#o`/`#f`/`#p` lazy-append
      (was eager same-pair replace for `#o`).
   d. **DSL — EDN redesign so `kind` is positional, then stream. Done (B, 2-vector + true
      streaming).** The blocker was that
      `:ligand-symmetry`/`:fluxionality` values need `kind` (perm degree) to parse, but `kind` is a
      map key (`{:kind … <constraint-key> <value>}`) whose position EDN does not fix — so a
      streaming reader can't see it before the value, forcing the tree bridge. We will **not**
      force tree-parse, and **not** assume map-key order. Fix: mirror the AST tuple
      `Constraint::StereoAtom(id, kind, constraint)` in the EDN by giving `kind` a
      **container-fixed** position — see "Stereo constraint EDN" below for the exact form. Then
      `read_stereo_*_constraint_dsl` read `kind` positionally, then the value with kind known
      (fully incremental, no slice capture, no order assumption); drop the two
      `// TODO: FIX THIS TO USE streaming parser`.
      **Done:** `StereoAtomConstraintDsl`/`StereoBondConstraintDsl` `FromEdn`/`ToEdn` are the
      2-vector `[<kind> {<key> <value>}]` (all fixtures migrated). Streaming value readers
      (`read_member`, `read_perm_vov`, `read_relation_value` + `RelationValue`,
      `read_ligand_symmetry`/`read_topicity`/`read_stereogenicity`) live in `dsl/constraint.rs`
      with the other EDN constraint readers; `relation_serde!` gained a streaming `$from_parts`
      and `stereo_kind_from_edn` was split to share `stereo_kind_from_name` (both in `dsl/stereo.rs`,
      `pub(crate)`). `read_stereo_*_constraint_dsl` now read the 2-vector incrementally (kind first,
      then the single-key payload); the `read_value_slice → read_string → FromEdn` bridge and both
      `// TODO` are gone. (Module-size cleanup of `dsl/constraint.rs` deferred.)
   e. Replace the `_from_edn`/`_to_edn` helpers (including inside macros) by inlined
      `FromEdn`/`ToEdn` impls.
   f. **Tests** — DSL string + EDN streaming roundtrip; proptest generators; a stereo
      collection lattice-law property test (the other six have one; stereo currently has none).
3. **Per-entity enums + collections. Each per-entity enum gets
   `Canonicalize` = delegate to the inner value (replacing `simplify`); each collection gets
   `Canonicalize` (per-kind canonicalize + drop-vacuous) and keeps its **hand-written**
   `Lattice`. **No shared trait / macro / generic `KeyedConstraints`** — the collections
   already each hand-write their full surface; we leave that duplicated (WET) rather than
   abstract it. **Done**

   The only collection whose `meet`/`join`/`canonicalize`/`matches` gain new logic are the
   three ring-bearing ones (atom/bond/dative): they must **group ring entries by `RingScope`**
   and merge per scope (the former `RingSize` loop, now keyed). That is done inline, per
   collection, in P5.1(b–d). WET means **no generic / shared** key machinery (no cross-family
   trait, macro, or `KeyedConstraints`) — each family still hand-writes its **own**
   `<Enum>ConstraintKey` + `key()` + by-key API mirroring its by-kind API (P5.1(i)). Ring's
   keyed accessors are `ring_count()` / `ring_size_count(s)` / `ring_memberships()` on each of
   the three collections. Every other kind stays one-per-kind; its existing per-accessor
   `meet`/`join` is untouched.
4. **Relational** (`RelationalConstraint`) — `Canonicalize` (canonicalize inner values;
   refs unchanged); **not** a `Lattice`.
5. **Molecule** (`MoleculeConstraint`) — `Canonicalize` (canonicalize payloads; atom-sets
   sorted; `SubPattern` recurses into the inner `MoleculeAst`); **not** a `Lattice`.
6. **Logical** (`Constraint` `And`/`Or`/`Not` + the `Constraints` `Vec`) — `Canonicalize`:
   recurse, flatten nested same-combinator, sort + dedup children and the top `Vec` by the
   `Constraint` declaration order, drop empty `And`/`Or`; **not** a `Lattice`.
7. **Sweep + exit.** Replace remaining `simplify`/`simplify_*` with canonicalize calls;
   add the deferred P4 entity `Canonicalize` derives (`Atom/Bond/Aromatic/Multicenter/
   Dative`), now unblocked.

#### Stereo constraint EDN (P5.2 d)

Stereo constraints have two serialization surfaces:

- **Form A — entity inline string.** `StereoAtomDsl`/`StereoBondDsl` serialize the whole element as
  one EDN string (`:type "Th0#o=(0,1)#g/"`, or `:ccw`/`:cw` shorthand); the `#p`/`#f`/`#o`/`#g`
  predicates ride inside, after the `<Kind><coset>` head. Kind is read from the head before any
  predicate, so perm degree is known. **No change** — Form A is fine.
- **Form B — molecule-scope `Constraint`** (`Constraints(Vec<Constraint>)`). A single constraint
  detached from its element. The AST variant is `Constraint::StereoAtom(StereoAtomId, StereoKind,
  StereoAtomConstraint)` — **`kind` is a separate positional field**. The DSL boundary
  `ConstraintDsl::StereoAtom(StereoAtomRef, StereoAtomConstraintDsl)` serializes as the single-key
  map `{:stereo-atom [<ref> <constraint>]}` (the generic 2-element `[ref, constraint]` leaf), and
  `StereoAtomConstraintDsl` re-bundles kind back as a **map key**: `{:kind <k> <constraint-key>
  <value>}`.

Per-kind value forms (`<vov>` = vector-of-cycles, degree from kind; `<rel>` = `:undetermined` |
`:homotopic`… | `[:a :b]` | `[:x] :member :not-in`):

| constraint-key | value | needs kind |
|---|---|---|
| `:ligand-symmetry` | `{:perm <vov> [:orientation :improper] [:member :not-in]}` | yes |
| `:fluxionality` | `<vov>` | yes |
| `:topicity` | `{:pair [i j] :relation <rel> [:member :not-in]}` | no |
| `:stereogenicity` | `{:relation <rel> [:member :not-in]}` | no |

**Problem:** the DSL discards the AST's positional `kind` and makes it a map key; EDN map-key order
is not container-fixed, so the streaming reader can't read `kind` before a perm value → tree bridge.

**Redesign (chosen: B — nested 2-vector).** Mirror the AST tuple `StereoAtomConstraintDsl(StereoKind,
StereoAtomConstraint)` by making the *constraint's own* EDN a positional **2-vector**, leaving the
generic entity leaf untouched:

```
{:stereo-atom [<ref> [<kind> {<constraint-key> <value>}]]}
{:stereo-atom [a3 [:tetrahedral {:topicity {:pair [0 1] :relation :homotopic}}]]}
```

`<kind>` keyword (`:tetrahedral` | `:cis-trans` | `:axial` | `:square-planar` |
`:trigonal-bipyramidal` | `:octahedral`); the constraint payload is its own single-key map
`{<key> <value>}` (the four value forms above, **minus** the `:kind` key).
`StereoAtomConstraintDsl`'s `FromEdn`/`ToEdn` change from the `{:kind …}` map to the 2-vector
`[<kind> {<key> <value>}]` — it stays a **self-contained `C`**, so the generic `(ref, C)`
machinery (`parse_entity_leaf` / `read_entity_leaf` / leaf `to_edn`) is unchanged. The streaming
reader reads the 2-vector positionally: `kind` (elem 0) → degree known → stream the value (elem 1).

**Why not A (flat `[<ref> <kind> <constraint>]`).** The entity leaf is generic over `(Ref, C: FromEdn
+ ToEdn)` and hard-codes 2 elements on all three paths (`parse_entity_leaf` len-2; `read_entity_leaf`
`[`/ref/inner/`]`; leaf `to_edn` `[ref, c.to_edn()]`). A flat 3-element leaf forces 6 bespoke stereo
leaf sites and strips `StereoAtomConstraintDsl` of any self-contained codec (its `FromEdn`/`ToEdn` go
dead). The original kind-in-map design existed precisely to keep `StereoAtomConstraintDsl` a
self-contained `C`; B preserves that while making `kind` positional. Form A (inline string) and the
`#g`/`#o`/`#f`/`#p` grammar are untouched.

### P6 · `Lattice`-trait flip + macro (lands once P1–P5 all impl `Canonicalize`)
- `Lattice: Canonicalize`; `matches` becomes the `meet`-derived default; `join` stays
  `Self`. `#[derive(Lattice)]` generates `meet`/`join` = field-wise + `canonicalize`.
  Remove the now-redundant hand-written `matches` impls (keep only genuine cheaper
  overrides). Hand-written leaf `meet`/`join` already canonicalize from P1–P5.

### P7 · Semantic-equality adoption + boundary cleanup
1. 095 Q2 leak audit → route to `equiv`/`Canonical<T>` where semantic keys are needed:
   `MoleculeAst::PartialEq` (graph-canonical), `HashMap`/`HashSet` AST keys, alias
   `BiBTreeMap`, constraint dedup. Decide structural-vs-semantic per site.
2. Drop `matches_value` (replace by as_lit_matches where appropriate); delete `capture`/`evaluate`.
3. Literal renames + `*Dsl` per the boundary-type convention (one `*Dsl` per
   boundary-crossing type, owning its serde; literals drop `Ast`):
   `StereoConfigurationDsl`→`StereoSiteDsl`; add `TopicityDsl`/`StereogenicityDsl`.

### P8 · Verification (doc 111)
1. C4e.5(1): all retained `lattice::` tests green (raised `PROPTEST_CASES`); demoted types
   leave the sweep.
2. C4e.5(2): atom-DSL roundtrip; C4e.5(3): `canonicalize` idempotence beyond `ValueAst`.
3. `umol-ast` lib + `--features proptest` + workspace build + conformance all green.

## Per-type review

`Lattice` is supplied two ways: **`#[derive(Lattice)]`** (`umol-ast-macros`,
field-wise over named struct fields; all-field-top is the top) and
**hand-written** impls (enums, tuple structs, collections). Cross-field work lives
in **`Canonicalize::canonicalize`** (`fn canonicalize(self) -> Result<Self,
Contradiction>`): pure-field-wise types `#[derive(Canonicalize)]`; a type with a
cross-field step hand-writes `canonicalize` (field-wise + the step inline).
`meet`/`join` are a field-wise op + `canonicalize`, so there are no separate
`post_meet`/`post_join` hooks — and since `Lattice: Canonicalize`, `derive(Lattice)`
just delegates to `self.canonicalize()`, so **neither derive takes an attribute**.
Only the stereo elements hand-write a cross-field step (the coset/kind collapse);
every other entity — `AtomAst` included, now that `JointDomain` is removed — is pure
field-wise. Disposition codes:
`derive` / `hand` keep `Lattice`; `matches` = inherent `matches`, no trait; `—` =
neither (graph-mediated or pure key). The "cross-field step" column names the
hand-written step where present.

### Entity ASTs

| Type | Lattice | Index → per-index ASTs | `canonicalize` cross-field hook | On-construction simplification |
|---|---|---|---|---|
| `AtomAst` | derive | — | — | field-wise; each field canonical (no cross-field hook — `JointDomain` removed) |
| `BondAst` | derive | — | — | field-wise |
| `NoncovalentBondAst` | derive | — | — | field-wise |
| `DativeBondAst` | derive *(after birelation drops `acceptor_slot`)* | — | — | field-wise |
| `AromaticSystemAst` | derive | membership (external, not a field) | — | field-wise `(electrons: ElectronCountsAst, charge, spin, constraints)` |
| `MulticenterBondAst` | derive | membership (external) | — | field-wise (as above) |
| `StereoAtomAst`, `StereoBondAst` | derive | `configuration: StereoConfigurationAst` (kind+coset, self-contained) | — (kind↔coset hook is *inside* `StereoConfigurationAst`) | field-wise `(configuration, constraints)` |
| `MoleculeAst` | — (graph) | — | — | graph canonicalization (WL / canonical rank), not an AST normal form |
| `ReactionRuleAst` | — (graph) | — | — | graph-level |

### Leaf / predicate ASTs

| Type | Lattice | Index | On-construction simplification |
|---|---|---|---|
| `ValueAst` | hand (`Undetermined`/`Lit`/`LitSet`/`Term`/`Predicate`) | — | `LitSet` (`BTreeSet`) canonical; `Term` n-ary fold→`Lit`; `Predicate` NNF→⊤/⊥ (lift to `Undetermined`/`Err`); residual opaque |
| `ValueTerm` | *(in `ValueAst::Term`)* | — | n-ary `Sum`/`Product` flatten+sort, const-fold, identity/annihilator, `Neg` normal form |
| `ValuePredicate` | *(in `ValueAst::Predicate`)* | — | NNF; `And`/`Or` flatten+sort+dedup + ⊤/⊥ fold (private carrier) → lift at `ValueAst` |
| `ElementAst` | hand (`Undetermined`/`Lit`/`LitSet`/`NotSet`/`Var`) | — | cardinality-canonical (smaller side → positive/complement; tiebreak positive); sort/dedup; singleton→`Lit`; `Not` dropped; empty→`Err` |
| `IsotopeMassAst` | hand (`Undetermined`/`Natural`/`Lit`/`LitSet`/`Var`) | — | positive-only (no negation, no cardinality); sort/dedup; singleton→`Lit`; empty→`Err`; `Natural` ground |
| `NoncovalentBondKindAst` | hand (enum) | — | none (flat) |
| `ElectronCountsAst` | hand (`Undetermined`/`Lit(Vec<i64>)`) | — | none (atomic value; exact-match `meet`/`join`), i64 to mirror valence fields (ValueAst) |
| `StereoKindAst` | hand (`Undetermined`/`Lit`) | — | none (flat; set-lattice deferred like `NoncovalentBondKindAst`) |
| `StereoConfigurationAst` | hand (struct `{kind, coset}`, self-contained via `kind`) | — | kind-aware coset fold; `kind not Lit ⇒ coset Undetermined` (element-side joint AST) |
| `StereoSiteAst` | hand (tristate; `Stereo(StereoKind, coset)`) | — | coset folded under the arm's concrete `kind`; `NotStereo` ≠ `Stereo` (constraint-side, `#T`/`#C`) |
| `StereoCosetAst` | — (kind-relative; **no `Lattice`**) | — | bare coset value; normalized by the owning `kind` (set/operator fold) |
| `StereoTerm` | — (kind-relative sub-part) | — | `Var`-rooted; word composes to one net perm → `Var`/`Mirror`/`Swap`/`Apply` by priority Mirror>Swap>Apply (owner folds under `kind`) |
| `AromaticValenceAst` | hand (enum) | — | delegate `ValueAst`; `NotAromatic` ≠ `Aromatic(_)` |
| `MulticenterValenceAst` | hand (enum) | — | delegate `ValueAst`; `NotMulticenter` ≠ `Multicenter(_)` |
| `SpinStateAst` | derive (field-wise `unpaired`, `multiplicity`) | — | derive (field-wise; **no** cross-field parity gate — tier-2, not enforced at the AST) |
| `TopicityRelationAst` | hand (`relation_ast!`) | — | finite-domain `LitSet`/`NotSet` sort/dedup/collapse + negation polarity |
| `StereogenicityAst` | hand (`relation_ast!`) | — | as above; flattened — `relation_ast! { StereogenicityAst, Stereogenicity }`, no `*RelationAst`, no wrapper; macro also generates `AsLit` (`Lit`→`Some`) |
| `TopicityAst` | **matches** | `pair` (essential — topicity is per-pair) | per-pair `rel` (`TopicityRelationAst`) is the lattice |
| `LigandPermutation` | matches (= equality) | — | none; EDN via `LigandPermutationDsl` |
| `OrientedLigandPermutation` | matches | — | none; EDN via `OrientedLigandPermutationDsl`; field `permutation` |
| `LigandSymmetry` | matches | — | none; field `permutation` |
| `Fluxionality` | matches | — | none; field `permutation` |
| `LigandPair` | — (key) | *is* an index | normalize `first ≤ second`; EDN via `LigandPairDsl` |

### Constraint ASTs

The per-constraint **enums** (`AtomConstraint`, …) are not themselves `Lattice`;
the **collections** are. None need a `canonicalize` cross-field hook.

| Type | Lattice | On-construction simplification |
|---|---|---|
| `AtomConstraints` (+ `AtomConstraint`) | hand (collection) | fixed kind order; dedup by kind; inner values canonical |
| `BondConstraints` | hand | as above |
| `DativeBondConstraints` | hand | as above |
| `MulticenterBondConstraints` | hand | as above |
| `AromaticSystemConstraints` | hand | as above |
| `NoncovalentBondConstraints` | hand (trivial; inner enum uninhabited) | none |
| `StereoAtomConstraints` (+ `StereoAtomConstraint`) | hand (macro) | fixed kind order; inner (LigandSymmetry/Fluxionality/Topicity/Stereogenicity) canonical |
| `StereoBondConstraints` (+ `StereoBondConstraint`) | hand (macro) | as above |
| `Constraints` / `MoleculeConstraint` (molecule scope) | — (flat `Vec`: combinators + predicates) | — |
| `RelationalConstraint` | — | cross-entity references (dative donor/acceptor, aromatic-system atoms, …); molecule-level; not a lattice |

`StereogenicityAst` flattens to `relation_ast! { StereogenicityAst, Stereogenicity }`
directly — a per-site yes/no, no `pair`, so no wrapper and no `*RelationAst`.
`TopicityAst` does not flatten: topicity is intrinsically per-pair (no generic
topicity, none for non-binary kinds), so `pair` is essential and `TopicityAst {
pair, rel }` stays a `matches`-only carrier with `TopicityRelationAst` as the
per-pair lattice (decision 10).

### `canonical()` fast path (by type)

`canonical()` default = clone + `canonicalize`; a type overrides it only to **borrow**
cheaply-known-canonical values (`Cow::Borrowed`, no clone). Three classes:

| Class | Types | `canonical()` |
|---|---|---|
| identity `canonicalize` | `NoncovalentBondKindAst`, `StereoKindAst`, `ElectronCountsAst` | always `Cow::Borrowed` |
| pure field/variant-wise | `AtomAst`, `BondAst`, `NoncovalentBondAst`, `DativeBondAst`, `AromaticSystemAst`, `MulticenterBondAst`, `AromaticValenceAst`, `MulticenterValenceAst`, constraint enums | **derived** — borrow `self` iff every field/payload borrows, else assemble owned |
| folding / cross-field (hand-written) | below | borrow listed variants, else clone + fold |

Hand-written (borrow → otherwise fold):
- `ValueAst`: `Undetermined`, `Lit` — `LitSet`/`Term`/`Predicate` fold.
- `ElementAst`: `Undetermined`, `Lit` — `LitSet`/`NotSet`/`Var` fold.
- `IsotopeMassAst`: `Undetermined`, `Natural`, `Lit` — `LitSet`/`Var` fold.
- `TopicityRelationAst`, `StereogenicityAst`: `Undetermined`, `Lit` — `LitSet`/`NotSet` fold.
- `SpinStateAst`: borrow iff both fields are ground-canonical **and** `are_compatible`
  (parity); else fold.
- `StereoConfigurationAst`, `StereoSiteAst`: borrow `Undetermined` / `NotStereo` / ground
  `(kind: Lit, coset: Lit)` when `are_compatible(kind, coset)` (coset index ∈
  `[0, kind.count())`); else fold (kind-aware). `kind` not `Lit` + concrete coset →
  collapse `coset → Undetermined` (normalization); out-of-range when `kind = Lit` → `Err`.
- constraint collections (`AtomConstraints`, …): borrow if kind-sorted, deduped, no
  vacuous (`Undetermined`-payload) entry; else fold.

No `Canonicalize`, so no `canonical()`: `StereoCosetAst`, `StereoTerm` (kind-relative —
folded by their owner).

## Per-type detail (bottom-up)

### `ValueAst`

**Variants (resolved).**
- `ValueAst = Undetermined | Lit(i64) | LitSet(Box<BTreeSet<i64>>) | Term(Box<ValueTerm>) | Predicate(Box<ValuePredicate>)`. The old `Expr(ValueExpr)` splits **by sort** into a `Term` (i64-valued) and a `Predicate` (bool-valued field constraint). `Bind`/`Ref` are dropped — a variable is the identity instance `Var` (no privileged form; cross-field/cross-var relations go through the molecule-level variable-constraint facility, doc 115). `Set` → `LitSet` (`BTreeSet`).
- `ValueTerm = Lit(i64) | Var(String) | Neg(Box) | Sum(Vec) | Product(Vec) | Div(Box,Box) | Rem(Box,Box)` — n-ary `Sum`/`Product` (so `+`/`*` are associative + commutative by construction); `Sub` folds to `Sum([a, Neg(b)])`; `Div`/`Rem` stay binary.
- `ValuePredicate = Rel(ValueTerm, RelOp, ValueTerm) | Mem(ValueTerm, BTreeSet<i64>) | Not(Box) | And(Vec) | Or(Vec)` — **no ⊤/⊥ variant** (see canonicalization).
- Value-scoped, replacing `ValueExpr`; concrete over `i64` (the only domain with an algebra — not parameterized). No `Predicate` abbreviation anywhere.

**Use-sites.** `ValueAst` is the single integer-valued field/predicate type; there is
no per-quantity type (no separate "bond order" AST). Every site below accepts the full
grammar — `Lit`, `LitSet` (`{1,2}`), `Var` (`?x`), domain-restricted `Var`
(`?x :: {…}`), `Term`, `Predicate` — parsed by the shared `value` parser and rendered
by `fmt_value`, inline through the owning entity's `*Dsl`.

| Field | Entity | Meaning |
|---|---|---|
| `order` | `BondAst`, `DativeBondAst` | bond order |
| `charge` | `AtomAst`, `BondAst`, `AromaticSystemAst`, `MulticenterBondAst` | formal charge |
| `implicit_hydrogens` | `AtomAst` | H count |
| `lone_pairs` | `AtomAst` | lone-pair count |
| `unpaired`, `multiplicity` | `SpinStateAst` | spin |
| ring counts (`RingMembership` value, keyed by `RingScope`) | `AtomConstraint`, `BondConstraint` | ring membership |

(`edit.rs` `old`/`new` pairs are deltas over these same fields, not new semantics.)
**Bond order** carries no fractional/aromatic value: an aromatic bond is localized
`order: Lit(1)` plus the `#a` flag (`Bond` is the localized object). Order shares
`ValueAst`'s open `i64` domain — no positive-only or range constraint baked in, same as
charge.

**Lattice** (hand-written on `ValueAst`, `Lattice: Canonicalize`). Top
`Undetermined`; no bottom value (`meet → None`). `meet`: wildcard / `Lit==` /
`Lit ∈ LitSet` / `LitSet ∩ LitSet` / `Term==`, `Predicate==` (structural — sound
once canonical); incompatible symbolic vs other → `None`. `join`: `Undetermined`
absorbs / `Lit ∪` → `Lit | LitSet` / equal symbolic kept / else `Undetermined`.
`matches` is the `meet`-derived default. `AsLit = i64`. The earlier operand-order
and non-reflexive-`matches` bugs vanish — `LitSet` is canonical by type and
`matches` comes from `meet`.

**Boundary.** `ValueDsl(pub ValueAst)` — `FromEdn`/`ToEdn`/`Display` over the
surface (`*`, `n`, `{a,b}`, term/predicate syntax). The **parser is faithful** (no
folding); all normalization is in `canonicalize`, run **lazily** (at compare/`meet`
and at the `*Dsl` boundary for canonical rendering).

**Canonicalization.** `ValueTerm: Canonicalize` is closed (a term reduces to a
term; a ground term → `Lit`). `ValuePredicate` is canonicalized **by `ValueAst`** (a
predicate can leave the sort to ⊤/⊥), threading ⊤/⊥ through a **private**
`enum { Top, Bottom, Predicate(ValuePredicate) }` carrier — never a public variant.

Full primitive step list. Annotations: **[new]** = added to `canonicalize`;
**[move]** = currently done elsewhere (parser fold or the old `simplify`) and
moved here; the parallel removal is named.

*`ValueTerm`*
1. `Lit(n)`, n<0 → `Neg(Lit(-n))` (sign in `Neg`; `Lit` ≥ 0). **[new]** — invariant for programmatic construction; the parser already emits `Neg(Lit)` for `-n`.
2. `Neg(Neg e)` → `e`. **[move]** — remove the `unary_expr` sign-XOR collapse (`dsl/value.rs:438`); the parser instead emits one `Neg` per `-` (drops `+`), and this step collapses. Also subsumes the old `ValueExpr::simplify` `Neg(Neg)`.
3. `Neg(Lit 0)` → `Lit 0`. **[new]**
4. `Sum`/`Product`: flatten nested same-op. **[new]** — subsumes the old `ValueExpr::simplify` flatten; the parser preserves paren-nesting faithfully.
5. `Sum`/`Product`: sort operands (canonical order). **[new]**
6. `Sum`/`Product`: const-fold the literal operands to one `Lit`. **[new]**
7. Identity: drop `Lit 0` from `Sum`, `Lit 1` from `Product`. **[new]**
8. Annihilator: `Lit 0` ∈ `Product` → `Lit 0`. **[new]**
9. Unwrap single; empty `Sum`→`Lit 0`, empty `Product`→`Lit 1`. **[move]** — old `ValueExpr::simplify` single-child unwrap.
10. `Div`/`Rem`: const-fold `(Lit a) op (Lit b)` when `b ≠ 0`. **[new]**

*`ValuePredicate`* (⊤/⊥ via the private carrier)
11. `And`/`Or`: flatten nested. **[move]** — old `ValueExpr::simplify` flatten; parser keeps nesting faithful.
12. `And`/`Or`: sort + dedup operands. **[new]**
13. `And`: drop ⊤, any ⊥ → ⊥. `Or`: drop ⊥, any ⊤ → ⊤. **[new]**
14. Unwrap single; empty `And`→⊤, empty `Or`→⊥. **[move]** — old single-child unwrap.
15. `Not(Not e)` → `e`. **[new]** — parser emits nested `Not` faithfully (boolean `!` is *not* folded).
16. De Morgan: `Not(And xs)`→`Or(¬xs)`, `Not(Or xs)`→`And(¬xs)` (NNF). **[new]**
17. `Not(Rel(a,op,b))` → `Rel(a, flip(op), b)`; so `Not` survives only on `Mem`. **[new]**
18. `Rel(Lit a, op, Lit b)` → ⊤/⊥. **[new]**
19. `Rel` canonical orientation: `Ge`/`Gt` → `Le`/`Lt` with operands swapped; canonical operand order. **[new]**
20. `Mem` set sorted + deduped (`BTreeSet`). **[new]**
21. `Mem(e, {a})` → `Rel(e, Eq, Lit a)`. **[new]**
22. `Mem(e, {})` → ⊥. **[new]**

*Lift (in `ValueAst::canonicalize`)*
23. `Term` → `Lit(n)` becomes `ValueAst::Lit(n)`. **[move]** — old `ValueAst::simplify` `Expr(Lit)→Lit` / `Expr(Neg(Lit))→Lit`.
24. `Predicate` → ⊤ becomes `Undetermined`; → ⊥ becomes `Err(Contradiction)`. **[new]**
25. `LitSet([n])` → `Lit(n)`; empty `LitSet` → `Err(Contradiction)`. **[new]** (`BTreeSet` keeps it sorted/deduped by type).

**Parallel removals.** Delete the parser's `unary_expr` sign-XOR (step 2); switch
`add_expr`/`mult_expr` to build `Sum`/`Product` and `Sub`→`Sum([a, Neg b])`,
`Div`/`Rem` (removes `BinOp`); delete `ValueAst::simplify` and
`ValueExpr::simplify` entirely (subsumed by `canonicalize`). The old
`Expr(Var)→Ref` / `Expr(Mem(Var,set))→Bind` lifts vanish with `Bind`/`Ref`.

**Equality/hash.** Derived structural `==`/`Hash`/`Ord` ("same tree"); semantic for the
decidable fragment via `equiv` (canonicalize-on-compare): `Lit(5)` / a ground `Term` / a
singleton `LitSet` collapse to one form, and `Sum`/`And` ordering is canonical. Residual
symbolic forms (free vars, undistributed products) stay opaque — structurally
compared, sound but incomplete. `Ord` is the structural order (for `BiBTreeMap`).

### `ElementAst`

**Variants (resolved).**
```
ElementAst = Undetermined | Lit(Element) | LitSet(BTreeSet<Element>) | NotSet(BTreeSet<Element>)
           | Var(Box<(String, Option<(MemOp, BTreeSet<Element>)>)>)
```
- `Element` is a finite enum (~118 variants, `Ord` by atomic number); the universe `U` is its variants.
- `Set` → `LitSet` (set-typed). `Not(Element)` **dropped** — the negative singleton (`!H`) is rare and `NotSet({H})` covers it; removes a canonicalization rule and makes the set-lattice core exactly the `relation_ast!` shape (`Undetermined | Lit | LitSet | NotSet`).
- Variables are one `Var`: `domain: None` = free (field = x, any element); `Some((MemOp, set))` = local membership (x ∈/∉ set). (Was `Ref`/`Bind`; merged — the taxonomy is exactly "free, or one membership constraint," and a single variant avoids the `VarDomain`-vs-"a var's domain" name clash.) `Var` is `ElementAst`'s **sole non-simple branch** — the finite-domain analog of `ValueAst`'s `Term`/`Predicate`: the concrete values (`Undetermined`/`Lit`/`LitSet`/`NotSet`) are simple, and since there's no algebra, a variable (free, or with one membership) is the *entirety* of the non-simple case, so one `Var` covers it where `Value` needs the two-sort split. Cross-var relations (`x = y`, `x ≠ y`) are molecule-level (the variable-constraint facility, doc 115), not field-level; a variable has no definition/use distinction — it's instantiated on first occurrence, domain = intersection of its occurrences.

**No `And`/`Or` needed.** A finite domain is closed under ∪/∩/∁, so any boolean
combination of memberships collapses to a single set: `x∈A ∨ x∈B` = `LitSet(A∪B)`,
`x∈A ∧ x∈B` = `LitSet(A∩B)`, `x∈A ∧ x∉B` = `LitSet(A∖B)`, `¬(x∈A)` = `NotSet(A)`.
`LitSet`/`NotSet` *are* the disjunctive/conjunctive/complement normal form, fully
collapsed — which is why finite types lack the `And`/`Or`/`Rel` that `Value` needs
(integers are infinite, so a predicate like "even" isn't a finite set). The only
compound logic that doesn't collapse is cross-variable / cross-field
(`x∈A ∨ y∈B`, `x = y`) — relational, hence the molecule-level variable-constraint facility (doc 115), not field-level.

**Polarity — cardinality-canonical, universe-relative.** Store a semantic set on its smaller side: positive (`Lit` singleton / `LitSet`) iff `|S| ≤ ⌊|U|/2⌋`, else by complement `NotSet(U∖S)`; tiebreak `≤` → positive (a 59-of-118 set is `LitSet`; `NotSet` only for ≥60). Boundaries: empty positive → `Err`, full → `Undetermined`, empty complement → `Undetermined`, full complement → `Err`. The `Var` domain canonicalizes identically (`Some(In, S)`, `|S| > 59` → `Some(NotIn, U∖S)`; vacuous → `None`; empty positive → `Err`). Element constraints are therefore **universe-relative** (`LitSet(117) ≡ "not C"`) — acceptable, the element universe is fixed.

**Lattice** (hand-written, `Lattice: Canonicalize`). Set algebra over the finite
universe `U`. Top `Undetermined` (= `U`); no bottom value (`meet → None` = ∅ =
contradiction). `meet` = intersection (canonicalized; ∅ → `None`); `join` = union
(canonicalized; → `NotSet`/`Undetermined` when large). Concrete forms
(`Lit`/`LitSet`/`NotSet`) combine by set algebra across polarities
(`LitSet ∩ NotSet` = setminus, `NotSet ∩ NotSet` = complement of union, …). `Var`:
same var → narrow/widen its domain; var-vs-concrete → field stays a `Var` with a
narrowed domain; different vars (cross-var equality) → the variable-constraint facility (doc 115).
`matches` is the `meet`-derived default (target's admissible set ⊆ pattern's),
replacing the current hand-written `element_set_view` version. `AsLit = Element`.

**Boundary.** **No `*Dsl`** (per the pragmatic policy — an entity field component
that only renders inside the atom string). Element round-trips via `AtomDsl` + the
shared `element` parser / `fmt_element` (`dsl/atom.rs`). Surface: `*`, `C`, `{C,N}`,
`!H` / `!{F,Cl}`, `?x` / `?x :: {…}`. The parser is **faithful** (raw
unsorted/undeduped sets, no folds) — nothing to remove; it needs only variant
updates: `!lit` → `NotSet({lit})` (`Not` dropped), `Set` → `LitSet`, `Bind`/`Ref` →
`Var`.

**Canonicalization** (`ElementAst: Canonicalize`) — all **[new]** (the parser does
no folding; the current `canonicalize_set`/`canonicalize_not_set` helpers, which
neither sort nor apply cardinality, are **replaced** by this impl; no element
`simplify` exists):
1. `LitSet`/`NotSet` sorted + deduped (free via `BTreeSet`).
2. Singleton positive `LitSet({e})` → `Lit(e)`.
3. Cardinality polarity: `|S| ≤ ⌊|U|/2⌋` → positive (`Lit`/`LitSet`), else complement `NotSet(U∖S)`; tiebreak positive.
4. Boundaries: empty positive → `Err`; full positive (= `U`) → `Undetermined`; empty `NotSet` → `Undetermined`; full `NotSet` → `Err`.
5. `Var` domain: canonicalize its membership set by the same polarity (flip `MemOp` on cardinality); vacuous (`In U` / `NotIn ∅`) → `None` (free); empty positive (`In ∅`) → `Err`.

**Equality/hash.** Derived structural `==`/`Hash`/`Ord`; semantic for concrete sets via
`equiv` (canonicalize-on-compare): each semantic set has one rep (sorted `BTreeSet`,
cardinality-chosen polarity, singleton → `Lit`). `Var` compares structurally (same
`id` + domain) — sound; cross-var equality is the variable-constraint facility (doc 115), not `==`.

### `IsotopeMassAst`

**Variants (positive-only).**
```
IsotopeMassAst = Undetermined | Natural | Lit(u32) | LitSet(BTreeSet<u32>)
               | Var(Box<(String, Option<BTreeSet<u32>>)>)
```
- **No `Not`/`NotSet`, no `MemOp` on `Var`** — unlike `ElementAst`. "Not isotope X" isn't a useful constraint (it admits nonsensical masses); the meaningful form is positive enumeration (`13C or 14C` = `LitSet({13,14})`). A deliberate use-case divergence from `Element` (whose negation *is* useful: "not C/H" = heteroatom). This also moots the open/finite-domain complement question — there is no complement.
- `Natural` is a **separate top-level ground** = natural isotopic abundance, distinct from any specific mass. `u32` masses (finite-in-principle ~≤295, not enforced — same stance as not validating element↔isotope consistency).
- `Set`→`LitSet`, `Ref`/`Bind`→`Var` (positive domain only), `BTreeSet`.

**Lattice** (hand-written, `Lattice: Canonicalize`). Top `Undetermined`; two branches under it — the positive mass-set lattice and the isolated `Natural` ground.
- Mass-set: `meet` = ∩ (∅ → `None`); `join` = ∪ of finite sets. No complement, no cardinality.
- `Natural`: `meet(Natural, Natural)` = `Natural`; `meet(Natural, mass)` = `None` (⊥ — `Natural` ≠ any specific isotope, decision **(a)**); `join(Natural, mass)` = `Undetermined`; `meet(Undetermined, Natural)` = `Natural`.
- `matches` = `meet`-derived default. `AsLit = u32` for `Lit` (and ground singleton `LitSet`); **`Natural` is ground but `as_lit` = `None`** (it denotes no single mass), so the `is_ground == as_lit.is_some()` alignment that holds for `Value` does *not* hold here.

**Boundary.** **No `*Dsl`** (pragmatic policy — field component, inline). Round-trips
via `AtomDsl` + the shared isotope parser (`#i…`). Parser change: **drop the `#i!…`
negative arms** (`Not`/`NotSet` gone); `Set`→`LitSet`, `Ref`/`Bind`→`Var`. Otherwise
faithful (raw sets, no folds).

**Canonicalization** — all **[new]** (parser faithful apart from dropping the
`#i!…` arms; replaces `canonicalize_isotope_set`, removes `canonicalize_isotope_not_set`):
1. `LitSet`: `BTreeSet` sorted+deduped; empty → `Err`; singleton `{n}` → `Lit(n)`.
2. `Var{id, Some(domain)}`: domain sorted+deduped; empty → `Err`; otherwise kept (singleton domain stays a `Var` — pinned but still a name for cross-ref); `None` = free var.
3. `Natural` / `Lit` / `Undetermined`: atomic.

No cardinality, no complement, and **no "full domain → `Undetermined`"** collapse (the mass domain is open — a positive set is never the whole universe).

**Equality/hash.** Structural over canonical; positive sets canonical (sorted, singleton → `Lit`); `Natural` a distinct ground; `Var` structural.

### `NoncovalentBondKindAst`

**Variants.** `Undetermined | Lit(NoncovalentBondKind)`. Kept minimal — a set-lattice
(`LitSet` of kinds, e.g. "H-bond or halogen-bond") is imaginable but not needed yet;
add `LitSet` only if bond-kind sets become useful.

**Lattice** (hand, `Lattice: Canonicalize`). `Undetermined` top; `meet(Lit a, Lit b)`
= `a==b ? Lit : None`; `join` = `a==b ? Lit : Undetermined`; `matches` from `meet`.

**Boundary.** No separate `*Dsl` — the kind is parsed/rendered inline (`kind_expr`:
`*`, `Hbd`/`Xbd`/…) within **`NoncovalentBondDsl`**, the bond's boundary
(`FromStr`/`Display`/`FromEdn`/`ToEdn`). *Finding:* `NoncovalentBondAst` also impls
`Display`/`FromStr`/`FromEdn` directly (`dsl/noncovalent.rs:91/99/105`), duplicating
the `Dsl` — per the boundary rule the AST should be serde-free; remove those (and
audit other entity ASTs for the same).

**Canonicalization.** Trivial — identity (no sets, no folding); `canonicalize`
returns `self`.

**Equality/hash.** Structural = semantic (flat enum).

### `SpinStateAst`

**Variants.** Struct, two independent `ValueAst` fields `unpaired`, `multiplicity`. No
enum, no collapse: the physical relation `m ≤ u+1 ∧ m ≡ u+1 (mod 2)` admits multiple
`m` per `u` (open-shell singlet `(2,1)` vs triplet `(2,3)`), so both are needed.

**Lattice** (`#[derive(Lattice)]`, field-wise over the two `ValueAst` lattices; top is
`(Undetermined, Undetermined)`). `matches` from `meet`. `AsLit = SpinState` — `as_lit`
succeeds only when both fields are `Lit` **and** form a valid `SpinState`; an invalid pair
(e.g. `(1,1)`) is a legal AST state with no `SpinState` literal (`as_lit = None`), so
`is_ground ≠ as_lit.is_some()` there. The validity check is *projection* to the domain
type, not invariant enforcement.

**Boundary.** No `*Dsl` (pragmatic policy — field component of `AtomAst`/`BondAst`).
Round-trips inline via the shared `#u`/`#s` predicates (`apply_spin_pair` /
`fmt_spin_pair`, `dsl/predicates.rs`). Surface `#u2`, `#s3`, `#u2#s3`. (`SpinState`'s
own `Display`/`FromStr` in `umol-shared` is the *literal's* surface — a ground domain
type, not an `*Ast` — so not the serde-on-AST violation.)

**Canonicalization** (`#[derive(Canonicalize)]`, pure field-wise): canonicalize each
`ValueAst` field; **no cross-field parity gate.** The `unpaired`↔`multiplicity` parity
relation is a **tier-2 physical, model-independent invariant** — the AST/parsing layer
enforces only tier-1 syntactic invariants, so parity is checked later (resolution) and
syntactically-valid-but-physically-invalid pairs are allowed. `From<(u8,u8)>` /
`From<SpinState>` **stay infallible** (no `TryFrom`). `Contradiction` here is only for a
field's own tier-1 failure (e.g. an empty `LitSet`), never parity.

Parallel removals: `simplify_values` (`spin.rs:27`) is subsumed by `canonicalize` but kept
until P4 (entity `simplify_values` still call it). `is_plus_sugar` was already retyped in
P1.

**Equality/hash.** Derived structural, semantic once each field is canonical;
`(2,1) ≠ (0,1)` (distinct states) is correct.

### `AromaticValenceAst` / `MulticenterValenceAst`

Structurally parallel (reviewed together; kept as two concrete types):
- `AromaticValenceAst = Undetermined | NotAromatic | Aromatic(ValueAst)`
- `MulticenterValenceAst = Undetermined | NotMulticenter | Multicenter(ValueAst)`

**Variants.** Tristate: top (`Undetermined`), explicit-negative (`NotAromatic` /
`NotMulticenter`), and a positive variant carrying the count as a `ValueAst`. The three
are genuinely distinct — in particular **`Aromatic(Lit(0)) ≠ NotAromatic`**: `av=0` is
*aromatic with zero aromatic valence* (e.g. tropylium-derived B), which the aromaticity
criterion treats as aromatic, whereas `NotAromatic` is "not in an aromatic system." So
no `Aromatic(Lit(0)) → NotAromatic` collapse.

**Lattice** (hand, `Lattice: Canonicalize`). Top `Undetermined`. `meet`: wildcard;
`Not* ∧ Not* = Not*`; `Not* ∧ Positive = None` (contradiction); `Positive(a) ∧
Positive(b) = Positive(a.meet b)` (`None` if the inner contradicts). `join`:
`Undetermined` absorbs; `Not* ∨ Not* = Not*`; `Not* ∨ Positive = Undetermined` (faithful
top exists, so `join` stays `Self`); `Positive(a) ∨ Positive(b) = Positive(a.join b)`.
The hand-written `matches` is **exactly the `meet`-derived default** (verified case by
case) → drop it, take the trait default (decision B.4). `AsLit = i64` (`Not* → Some(0)`,
`Positive(Lit n) → Some(n)`).

**Boundary.** Constraint-*value* payloads → they **keep** a `*Dsl`
(`AromaticValenceDsl`/`MulticenterValenceDsl`, `dsl/atom.rs:904/978`) per the pragmatic
policy. Two surfaces, both faithful:
- EDN (inside the constraints map): `:undetermined` / `:not-aromatic` /
  `{:aromatic <value>}` — via the `*Dsl`.
- atom-string predicate (inline in `AtomDsl`, `aromatic_valence`/`multicenter_valence`
  parsers): `#a*`↔`Undetermined`, `#a!`↔`Not*`, `#a+`↔`Positive(Undetermined)`,
  `#a<n>`↔`Positive(Lit n)`, bare `#a`↔`Positive(Lit 1)`. These are fixed surface
  spellings of variants (a round-trip bijection), **not** value simplification — nothing
  to remove from the parser.

**Canonicalization.** Delegate to the inner `ValueAst` (`Positive(v) →
Positive(v.canonicalize())`); `Not*`/`Undetermined` identity; no cross-variant fold.
Parallel removal: delete `simplify` (`constraint/atom.rs:192/304`, subsumed).

**Equality/hash.** Derived structural, semantic once the inner `ValueAst` is canonical;
the three states stay distinct.

**`matches_value` / `AsLit` (not an inconsistency).** `matches_value(0)` is `true` for
both `NotAromatic` and `Aromatic(Lit(0))`, and both `AsLit` to `Some(0)` — **correct**: a
non-aromatic atom genuinely has aromatic-valence *count* 0, so `0` is its real value. What
distinguishes `NotAromatic` from `Aromatic(0)` is aromatic-system membership, a separate
axis (tested by the aromaticity criterion), not the count. This differs from `NotStereo →
None` for a real reason: `NotStereo` means *no coset exists* (not a stereocenter), whereas
`NotAromatic` means *count 0*. Both `AsLit` impls are correct as-is.

### `StereoKindAst` (new — the `kind` field of the joint config AST)

**Variants.** `Undetermined | Lit(StereoKind)`. It is the `kind` field of
`StereoConfigurationAst` (element side): a lattice top so a pattern can mean "a stereo
center of any geometry," with `coset` bound to it. No `LitSet`/`NotSet` — deferred
exactly like `NoncovalentBondKindAst` (add only if kind-sets become useful). (The
constraint-side `StereoSiteAst` instead carries a *concrete* `StereoKind` in its `Stereo`
arm — its `Undetermined` variant absorbs the unknown-geometry case.)

**Lattice** (hand, `Lattice: Canonicalize`). Trivial 2-state: `Undetermined` top;
`meet(Lit a, Lit b) = a==b ? Lit : None`; `join = a==b ? Lit : Undetermined`; `matches`
from `meet`. Same shape as `NoncovalentBondKindAst`.

**Boundary.** No `*Dsl` (field component, pragmatic policy). Kind round-trips inline via
`StereoAtomDsl`/`StereoBondDsl` (the `#T`/`#C`/… config tags). No `StereoKindDsl` exists
or is needed.

**Canonicalization.** Trivial identity (flat, no sets). The `kind`↔`coset` coupling
(coset range depends on `kind.count()`) is **not** here — it is the joint AST's hook
(`StereoConfigurationAst::canonicalize`: `coset → Undetermined` when `kind ≠ Lit`; Entity table,
decision C.9).

**Equality/hash.** Structural = semantic (flat).

### `TopicityRelationAst` / `StereogenicityAst` (relation lattices)

Two `relation_ast!` finite-domain set lattices, identical shape, over their 3-element
ground enums:
- `TopicityRelationAst` over `Topicity = Homotopic | Enantiotopic | Diastereotopic`
- `StereogenicityAst` over `Stereogenicity = Symmetric | Prochiral | Stereogenic`

**Flattening (stereogenicity only).** Stereogenicity is per-*site* (no ligand pair), so
the relation *is* the value: `StereogenicityAst` becomes the `relation_ast!` type
directly. Drop the current separate `StereogenicityRelationAst`
(`constraint/stereo.rs:184`) and the `StereogenicityAst(StereogenicityRelationAst)`
newtype (`:250`) — a redundant wrapper. Topicity does *not* flatten (keeps the
`TopicityAst { pair, rel }` carrier, below).

**Variants.** `Undetermined | Lit(x) | LitSet(BTreeSet) | NotSet(BTreeSet)` —
`Undetermined` top (any), `Lit` single, `LitSet` positive set, `NotSet` complement. Set
storage `Vec` → `BTreeSet` per the set-typed-storage decision (`NotSet(vec![…])` today).

**Lattice** (`relation_ast!`-generated, `Lattice: Canonicalize`). Finite-domain: meet =
intersection (∅ → `None`), join = union (full → `Undetermined`), `matches` from `meet`.
Macro also generates `AsLit` (`Lit → Some`).

**Boundary.** Constraint *values* (payloads of
`StereoAtomConstraint::Topicity`/`Stereogenicity`) → per the pragmatic policy they get a
`*Dsl`. **Finding:** today they are (de)serialized by hand-rolled free functions
(`topicity_to_edn`/`topicity_from_edn`, `dsl/stereo.rs:1034/1050`) and inline macro arms,
bypassing a boundary type — the `NoncovalentBondAst` anti-pattern. Fix:
`TopicityDsl`/`StereogenicityDsl`; `TopicityRelationAst` rides inside `TopicityDsl`
(field component → no `TopicityRelationDsl`).

**Canonicalization.** `relation_ast!` finite-domain: `LitSet` sort/dedup, singleton →
`Lit`, full set → `Undetermined`, empty → `Err(Contradiction)`; negation polarity (store
the smaller of positive/complement). Parser faithful.

**Equality/hash.** Structural = semantic once canonical (sorted sets + fixed polarity).

### `TopicityAst` (matches-only, per-pair carrier)

`TopicityAst { pair: LigandPairAst, rel: TopicityRelationAst }`. `pair` is **essential
identity** — topicity is intrinsically a relation between two specific ligands; there is
no generic, pair-free topicity (and none for non-binary kinds, which have many pairs).
So `pair` is neither removable nor liftable to a lattice dimension.

**Current → proposed.** Code currently `impl Lattice for TopicityAst`
(`constraint/stereo.rs:218`), combining `rel` only within a matching `pair`. Decision
C.10 makes it **matches-only**: different pairs are incomparable, so there is no global
top and it is not a single lattice. `matches` = `pair` equality ∧ `rel`-relation
`matches`. The lattice lives one level down (`TopicityRelationAst`); same-`pair` entries
in a `StereoAtomConstraints` collection are meet-combined there.

**Consequence (at implementation).** Remove `impl Lattice for TopicityAst` and its
fixed-`pair` lattice-law proptest (`property/lattice.rs:86`
`test_topicity_ast_lattice_laws`); the per-pair lattice is already covered by
`test_topicity_relation_ast_lattice_laws` (`:183`).

**Boundary.** `TopicityDsl` (carrying `pair` + `rel`), replacing the free-function serde
noted above.

### Stereo configuration cluster (`StereoConfigurationAst` / `StereoSiteAst` / `StereoCosetAst` / `StereoTerm`)

> **Superseded in part by P3 (revised 2026-06-17).** Current shapes: element config is
> `enum StereoConfigurationAst { Undetermined | Kind(StereoKind, StereoCosetAst) }`
> (`StereoKindAst` removed); the constraint side is two per-kind types
> `TetrahedralStereoAst`/`CisTransStereoAst` (not `StereoSiteAst`). The reasoning below
> about kind-relativity, the coset algebra, and `StereoTerm` still holds.

**Cosets are kind-relative.** A coset index/set/term is meaningless without its
`StereoKind` — and not only for evaluation: operator folding (`~Lit`), the identity
`~?x ≡ '?x` (true for chiral kinds, false for achiral), complement / cardinality-canonical
storage (the universe is `kind.count()`), the enantiomer map, and the meso pattern
(`?x` at one center, `'?x` at another) **all** need kind. So there is no honest kind-free
coset lattice; the crude structural approximation (opaque `Term`) can't fold, can't
identify `~?x` with `'?x`, can't do meso. **Kind therefore lives inside the AST**, on the
two self-contained owners; the bare coset/term carry no kind and are never operated on in
isolation.

**Types.**
```
StereoConfigurationAst { kind: StereoKindAst, coset: StereoCosetAst }   // in stereo elements
StereoSiteAst   = Undetermined | NotStereo | Stereo(StereoKind, StereoCosetAst)  // in #T / #C constraints
StereoCosetAst  = Undetermined | Lit(u32) | LitSet(BTreeSet<u32>) | NotSet(BTreeSet<u32>) | Term(Box<StereoTerm>)
StereoTerm      = Var(String, Option<BTreeSet<u32>>) | Swap(Box) | Mirror(Box) | Apply(Box, Permutation)
```
- `StereoConfigurationAst` (the joint AST, like `SpinStateAst`) is the **element**-side
  bundle. `kind` is `StereoKindAst` (`Undetermined | Lit`) so a pattern can mean "a stereo
  center of any geometry"; `coset` is bound to it.
- `StereoSiteAst` is the **constraint**-side tristate (`#T`/`#C` on an atom/bond). Its
  `Undetermined` absorbs the unknown-geometry case, so the `Stereo` arm carries a
  **concrete `StereoKind`** (a `#T` site is concretely tetrahedral), *not* `StereoKindAst`.
- `StereoCosetAst`, `StereoTerm` are kind-relative sub-parts: **no `Canonicalize`, no
  `Lattice`**; only derived structural `Clone`/`Eq`/`Hash`/`Ord`.

`Canonicalize` **and** `Lattice` are implemented on `StereoConfigurationAst` and
`StereoSiteAst` only — both self-contained via their own `kind`. The kind-aware coset
algebra (`canonicalize`/`meet`/`join`/`matches`/complement/operator-fold/enantiomer) are
helpers that take `kind` (on `StereoKind` / the coset space); the two owners call them
from their `self`-held kind.

**No `Predicate` sort.** The coset is finite with a group action but no order/arithmetic,
so (like `ElementAst`) memberships collapse to one set and there are no `Rel`/`And`/`Or` —
only the group-action `Term`. The set ops (complement, "store the smaller side") are
finite but **kind-relative** (universe = `kind.count()`), which is why they live in the
kind-aware algebra, not on the bare `StereoCosetAst`.

**No kind↔coset validity gate (tier-2).** The kind↔coset-index validity — when
`kind = Lit(k)`, coset indices lie in `[0, k.count())` — is a **tier-2 physical,
model-independent invariant**, so (like spin's u↔m parity) it is **not enforced at the
AST**: `canonicalize` does **not** range-check (an out-of-range index is a legal AST
state), and `meet`/`join` carry no gate. Validity is checked at resolution. The
kind-relative *syntactic* folds below (set sort/dedup, `Lit` collapse, operator compose)
are tier-1 and remain; they require a concrete `kind = Lit(k)` to run, so when
`kind = Undetermined` the coset keeps its kind-independent form (the exact handling is a
P3 detail).

**`StereoConfigurationAst::canonicalize` — full step list.**

Field-wise (the original draft's `i ≥ n` / `s ⊆ [0,n)` range-checks are **removed** —
tier-2 validity, not enforced here; see the no-gate note above):
1. canonicalize `kind` (`StereoKindAst`: identity).
2. `kind = Undetermined`: the kind-relative folds can't run; keep `coset`'s
   kind-independent form (exact handling a P3 detail — do **not** coerce or gate).

Otherwise `kind = Lit(k)`, `n = k.count()`; fold `coset` under `k` (no range-check):

`coset` = literal / set forms:
3. `Lit(i)`: identity.
4. `LitSet(s)`: `s = ∅` → `Err`; `|s| = 1` → `Lit`; `s` = full (`[0,n)`) → `Undetermined`;
   else sorted/deduped `BTreeSet`.
5. `NotSet(s)`: `s = ∅` → `Undetermined`; `s` = full → `Err`; cardinality polarity vs `n`
   — keep `NotSet(s)` iff `|s| < n − |s|`, else `LitSet([0,n) \ s)` (tiebreak positive),
   then re-apply step 4's collapse.
6. `Undetermined`: identity.

`coset` = `Term(t)` — operators are all permutation actions, so:
7. compose the operator word into **one net permutation `g_total`** over the inner `Var`
   (`Mirror` = μ_k, `Swap` = ι_k, `Apply(g)` = g; identity factors drop).
8. canonicalize the inner `Var(name, dom)`: if `dom = Some(s)`, sort/dedup, `s = ∅` →
   `Err`, `s` = full → `None` (no range-check — tier-2). `name` always kept; a `Var`-rooted
   term never folds to a `Lit`.
9. choose the representation by **priority Mirror > Swap > Apply**: `g_total = identity` →
   bare `Var`; `= μ_k` → `Mirror(Var)`; `= ι_k` → `Swap(Var)`; else → `Apply(Var,
   g_total)`. Consequence (since chiral `ι_k = μ_k`, achiral `μ_k = identity`): canonical
   `Mirror` is chiral-only, canonical `Swap` achiral-only, the tie only bites for chiral
   μ=ι (→ `Mirror`).

`Lattice`: `meet` = meet the two `StereoKindAst`s; if the result is `Lit(k)`, kind-aware
`coset` meet under `k`, else `coset = Undetermined`. `join` analogous; top
`(Undetermined, Undetermined)`. `matches` = meet-derived default. `AsLit =
StereoConfiguration { kind: StereoKind, coset: u32 }` — `Some` iff `kind = Lit(k)` ∧
`coset = Lit(i)`, else `None`.

**`StereoSiteAst`.** Tristate, structurally like `AromaticValenceAst` but the `Stereo`
arm pairs `(StereoKind, coset)`. `canonicalize` folds the `Stereo` arm's coset under its
concrete `kind` (steps 3–9). `meet`: wildcard / `NotStereo∧NotStereo` / `NotStereo∧Stereo
→ None` / `Stereo(k1,a)∧Stereo(k2,b)` = `k1≠k2 ? None : Stereo(k1, kind-aware coset
meet)`. `join` widens mismatches to `Undetermined`. `matches` meet-derived; `AsLit = u32`
(`NotStereo → None`).

**`StereoTerm`.** `Var` consolidates the old `Var`/`VarDomain` (optional finite domain);
`Swap`/`Mirror`/`Apply` rename `SwapOp`/`MirrorOp`/`ApplyOp`; no `Lit`/`LitSet` operand.
A canonical term is exactly one (or zero) operator layer over the `Var` — there is no
nesting, since the word composes to a single `g_total` (step 7) rendered per the priority
(step 9).

**Boundary.** No `*Dsl` for `StereoCosetAst`/`StereoTerm` (sub-parts; round-trip inline
via `StereoAtomDsl`/`StereoBondDsl`, config syntax `~ ' ^`). `StereoSiteAst` is the
constraint value → keeps a `*Dsl` (the renamed `StereoConfigurationDsl` → `StereoSiteDsl`).

**Parallel removals.** Delete `StereoCosetAst::simplify` / `StereoExpr::simplify`
(`ast/stereo.rs:205/315`) — the folding moves into the owners' `canonicalize` (kind from
`self`). The old `impl Lattice for StereoCosetAst` / `StereoConfigurationAst`
(`stereo.rs:247/145`) are removed; `Lattice` moves to the kind-carrying owners.

**Equality/hash.** `StereoCosetAst`/`StereoTerm` have only structural `==` (crude — `~?x`
vs `'?x` differ structurally even when a chiral kind identifies them). It becomes semantic
**after** the owner kind-canonicalizes the coset, which it does on `equiv`/`meet` (lazy) —
so `equiv` at the owner level is well-defined; a bare stored coset stays raw until then.
Cross-center meso (`?x`/`'?x`) is the relational
/ variable-constraint layer (doc 115), evaluating with each center's kind.

### `AtomAst`

**Fields.** `element: ElementAst`, `isotope_mass: IsotopeMassAst`,
`charge`/`implicit_hydrogens`/`lone_pairs: ValueAst`, `spin: SpinStateAst`,
`constraints: AtomConstraints` — every field a `Canonicalize`/`Lattice` type.

**Lattice/Canonicalize** — **pure field-wise**, `#[derive(Lattice)]` + `#[derive(Canonicalize)]`,
no attribute, **no cross-field hook** (with `JointDomain` removed, `AtomAst` has none).
`meet`/`join` = field-wise + `canonicalize`; `matches` = meet-derived default; top =
all-`Undetermined`, empty constraints. `canonicalize` = field-wise; any field `Err`
propagates (spin parity, empty `LitSet`, empty element set).

**Boundary.** `AtomDsl`; the field components (element/isotope/spin/valences) render
inline via the shared sub-parsers (pragmatic policy).

**Equality/hash.** Derived structural, semantic once fields are canonical.

**Parallel removals.** **`JointDomain` removed entirely** — `saturate_atom`,
`field_value_for_joint_var`, `narrow_joint_var_to_lit`, the `#[lattice(saturate)]`
attribute (and the macro's support), the `AtomConstraint::JointDomain` variant, and the
`Lattice::saturate` hook. Also `simplify_values` (`atom.rs:178`). `into_ground`/
`into_zeroed` are **grounding/defaulting** ops, *not* canonicalization — they stay.

### `ElectronCountsAst` + `AromaticSystemAst` / `MulticenterBondAst`

**`ElectronCountsAst` (new leaf).** `Undetermined | Lit(Vec<u8>)`, replacing both
entities' `electrons: Vec<ValueAst>`. Concrete counts only — no per-position `ValueAst`
pattern/var/expression; the whole per-atom vector is **one atomic lattice value**, with
the position-as-axis problem gone (it's compared whole, not cell-by-cell).
- **Lattice** (hand): `Undetermined` top; `meet(Undetermined, x) = x`; `meet(Lit a, Lit
  b) = a==b ? Lit(a) : None`; `join` unequal `Lit` → `Undetermined`; `matches` from
  `meet`. `AsLit = Vec<u8>`.
- **Canonicalize:** identity — the vector is positional (cell = member atom), so **no**
  sort/dedup; a `Lit(Vec<u8>)` is already canonical.
- **Boundary:** no `*Dsl` (field component; renders inline via
  `AromaticSystemDsl`/`MulticenterBondDsl`).

**`AromaticSystemAst` / `MulticenterBondAst`.** `{ electrons: ElectronCountsAst, charge:
ValueAst, spin: SpinStateAst, constraints }` — every field a lattice, so pure field-wise
`#[derive(Lattice)]` + `#[derive(Canonicalize)]`, **no cross-field hook**. Which atoms
belong to the system is external (molecule-level membership), not a struct field.
Constructors retype: `new(Vec<ValueAst>)`/`with_electrons(Vec<ValueAst>)` →
`ElectronCountsAst`; `from_electrons(Vec<u8>)` → `Lit`. Boundary: `AromaticSystemDsl` /
`MulticenterBondDsl` (entities). Parallel removal: `simplify_values`
(`aromatic.rs:92`, `multicenter.rs:92`).

### `BondAst` / `NoncovalentBondAst` / `DativeBondAst`

Pure field-wise, `#[derive(Lattice)]` + `#[derive(Canonicalize)]`, no cross-field hook.
`BondAst { order, charge, spin, constraints }` and `NoncovalentBondAst { kind,
constraints }` are unremarkable. `DativeBondAst` switches from its hand-written `Lattice`
to derive once the birelation promotion removes `acceptor_slot` (C.7). Parallel removals:
`simplify_values` on each (`bond.rs:105`, `noncovalent.rs:61`, `dative.rs:83`).
**Boundary finding:** `NoncovalentBondAst` impls `Display`/`FromStr`/`FromEdn` directly
(`dsl/noncovalent.rs:91/99/105`) — serde belongs on `NoncovalentBondDsl`; remove and
audit the other entities.

### Constraint types

**Per-entity enums** (`AtomConstraint`, `BondConstraint`, …) are **not** `Lattice` (a
fixed-kind constraint has no cross-kind meet). They impl `Canonicalize` = delegate to the
inner value's `canonicalize` (replacing `<Enum>::simplify`). Boundary-independent payloads
→ each keeps its `*Dsl`. `AtomConstraint::TetrahedralStereo`'s payload retypes
`StereoConfigurationAst` → `StereoSiteAst` (the rename); per-kind variants stay
(`TetrahedralStereo`, `CisTransStereo`): `#T` tetrahedral-specific, `#C` cis/trans-specific
(the duplicated concrete `StereoKind` in `StereoSiteAst::Stereo` is accepted). `JointDomain`
is removed (cross-field/cross-atom correlations → the molecule-level variable facility, doc
115).

**Per-entity collections** (`AtomConstraints`, `BondConstraints`, `DativeBondConstraints`,
`MulticenterBondConstraints`, `AromaticSystemConstraints`, `StereoAtomConstraints`,
`StereoBondConstraints`; `NoncovalentBondConstraints` trivial — inner enum uninhabited) are
`Lattice` + `Canonicalize`, and behave as **keyed maps** (key → value-lattice, unique per
key). The key is the constraint **kind**, plus a **sub-key** for the multi-entry kinds:
`RingScope` for ring, the pair/permutation for the stereo relations. Each collection is a
`Vec`-newtype (`AtomConstraints` a kind-sorted `SmallVec`) that **hand-writes** this surface
— WET, no shared macro/trait/generic (a future DRY pass could fold the duplication; deferred).
`simplify_each` removed everywhere; collections render inline via the owning entity's `*Dsl`.

**Keyed-collection lattice contract** (unique-per-key). The *product lattice over the
per-key value lattices*, each extended with a top; **absent key = top ≡ an
`Undetermined`-value entry**, which canonicalize drops. Each collection implements this by
hand over its key (for single-kind collections the key *is* the kind, so it reduces to the
existing per-accessor form; only ring/stereo collections add sub-key grouping). Laws hold
whenever every value is a proper lattice (`TopicityRelationAst`, `ValueAst`, the
valence/stereo value enums; `unit` trivially).
- **canonicalize** — group by key (one entry/key), merging same-key entries by `meet`-ing
  their values; drop vacuous (top/`Undetermined`-value ≡ absent); fixed key order. Same-key
  value contradiction → `Err(Contradiction)`. Unit-valued merge = dedup.
- **meet** (`a ∧ b`) — union of keys; key in both → value-`meet` (`None` if any
  contradicts); key in one → carried through; then canonical.
- **join** (`a ∨ b`) — keys present in **both** only; per shared key, value-`join`; keys on
  one side dropped (`join(top, x) = top`).
- **matches** — every key of `pattern` present in `target` with value-`matches`; ≡
  meet-derived `pattern.meet(target) == canonical(target)` (the form the P6 flip installs).
- **is_undetermined** = no entries; **is_ground** = every value ground.

**drop-vacuous moves into `canonicalize`.** Today a vacuous (`payload = Undetermined`)
constraint may sit in the AST, elided only at *render* time; `meet` already refuses to add
them. Canonical-by-construction *requires* dropping them in `canonicalize` so
`{Valence(Undetermined)}` and `{}` are structurally equal. A parsed `#v*` normalizes to no
entry — faithful (information-free top, like `1+1 → 2`), AST-level not just surface.
**Accepted** — roundtrip fidelity preserved.

**`add` semantics (decided).** `add` / `with_constraint` keep **replace** (set, last-wins,
infallible). Contradictions are *not* caught at add-time; they surface lazily at
canonicalize / `meet`, consistent with lazy canonicalization. So `add` (set) and the
canonical/lattice conjunction (per-key `meet`) intentionally diverge — `meet`/canonicalize
is the authority for the canonical form. An explicit conjoining verb, if ever wanted, is a
*separate* fallible `narrow` (= per-key `meet`), not an overload of `add`.

**Multi-entry kinds** (no subsumption multisets remain): ring membership, and the stereo
relations below. Each fact is **one constraint**; several may share a kind, distinguished by a
sub-key (`RingScope` for ring; the pair/permutation for stereo). They are **non-unique** —
`add` appends, and per-sub-key dedup is **lazy** (at `meet`/`canonicalize`); each collection
hand-writes the per-sub-key handling (WET — no shared container).

*Ring membership* (atom, bond, dative) replaces `RingCount` + `RingSize`. Membership is a
multiset `M : size → count` (`RingCount = |M|`, per-size = `M(s)`); we store bounds on its
projections as constraints `RingMembership(RingMembershipAst)`, where
`RingMembershipAst { scope: RingScope, count: ValueAst }` is a **single-entry carrier** (one
fact, parallel to `TopicityAst { pair, rel }`). `RingScope { All, Size(u8) }` (`All` = total
`|M|`, `Size(s)` = `M(s)`; `Ord` `All` first); the `count` is a `ValueAst`. One `RingMembership`
per scope on an entity; literal-int `Size` sub-keys are disjoint → no subsumption, clean
per-scope `join` (box-hull). The per-scope product lattice is hand-written in each collection's
`meet`/`join`/`matches`. Count sugar:
`n`→`Lit`, `+`→`var_at_least("r", 1)` (= `?r >= 1`; `ValueAst` has **no** `NotSet`, so `≥1`
is a `Predicate(Rel(Var, Ge, 1))` — the existing `dsl/predicates.rs` `+` sugar; the fixed
name `"r"` is fine here because nothing unifies variables yet — the anonymous-bound-vs-named-
var fix is deferred to doc 115), `0`→`Lit(0)` (canonical render `#R!`; bare `#R`→`Lit(1)`),
`{a,b}`→`LitSet`, `?v`→`Var`, `*`→`Undetermined`. Notation
`#R<count>` (`All`) / `#R(s)<count>` (`Size(s)`) — `#R+` (= SMARTS `R`), `#R2` (= `R2`),
`#R!` (acyclic, = `R0`), `#R(6)+` (= SMARTS `r6`), `#R(6)2` (naphthalene, beyond SMARTS),
`#R(6)!` (no six-ring). Cross-size disjunction ("5 or 6") is not a
map key (keys stay literal ints) — it goes through a variable size `#R(?r)…` (the variable
path); rings-as-entities would make it clean but needs graph-matching infra we
don't have. `All = Σ_s Size(s)` is a tier-2 validator check, not an AST-lattice invariant.

*Stereo `#o`/`#f`/`#p`* (stereo atoms and bonds) are keyed maps one level down; keys
canonical by construction, so dedup is exact-key grouping:
- `#o` Topicity — key `LigandPairAst` (unordered, normalized lower-first), value
  `TopicityRelationAst`. Per-pair dedup **meet**s the `rel`s (→ `Err` on contradiction) —
  stronger than `add`'s per-pair last-wins replace; canonicalize/meet is the authority.
- `#f` Fluxionality — key `PermutationAst`, **unit** value; dedup = drop identical.
- `#p` LigandSymmetry — key `(OrientedPermutationAst, MemOp)`, **unit** value; dedup = drop
  identical.
Only `#o` is per-pair (merges values); `#f`/`#p` are per-perm (drop duplicates).

**Molecule-level `Constraints` / `Constraint`.** `Constraint = Atom(AtomId, …) | Bond(…) |
… | Relational(…) | Molecule(…) | And(Vec) | Or(Vec) | Not(Box)` — a boolean combinator
tree over **ID-scoped** predicates; `Constraints` is a flat (conjunctive) `Vec<Constraint>`.
**Not a `Lattice`** (ID-bearing; the molecule order is graph subsumption, not algebraic).
`Canonicalize`: recurse into children, flatten nested same-combinator, sort + dedup, drop
empty `And`/`Or`; inner predicates canonical. **Canonical order = the `Constraint`
declaration order** (ABDAMNSS → `Relational` → `Molecule` → `And`/`Or`/`Not`); within a
variant by payload (ids / atom-sets / inner values, each canonicalized); combinator children
sorted recursively by the same total order. Equality is structural **with** the IDs.
`Constraint::simplify` → `Canonicalize`; boundary is `ConstraintsDsl`.

## Naming and boundary types

Convention clarified during the review:

- **`*Ast`** = in-memory type with *pattern* semantics — an `Undetermined` top and
  lattice / `matches` behaviour. Keep the suffix only for these (the value enums,
  `StereoAtomAst`/`StereoBondAst`, the relations, `StereoKindAst`, …).
- **Literals** (always concrete, no `Undetermined`) drop `Ast`: `PermutationAst`
  → `LigandPermutation`, `OrientedPermutationAst` → `OrientedLigandPermutation`,
  `LigandSymmetryAst` → `LigandSymmetry`, `FluxionalityAst` → `Fluxionality`,
  `LigandPairAst` → `LigandPair`.
- **`<Type>Dsl`** = the EDN/string **boundary** type implementing `FromEdn` /
  `ToEdn` (and `Display`/`FromStr` where applicable) — one per boundary-crossing
  type, including the literals: `LigandPairDsl`, `LigandPermutationDsl`,
  `OrientedLigandPermutationDsl`. **No manual DSL display/serde**, and nothing
  bypasses these: the current inline cycle-parsing in `dsl/stereo.rs`
  (`perm_cycles` → bare `umol_perm::Permutation`) is wrong and is replaced by the
  boundary types. The orphan rule (`FromEdn`/`ToEdn` foreign; `umol_perm` types
  foreign) forces these onto local newtypes — the standing reason the permutation
  newtype must exist; an inherent method on a foreign type is likewise impossible.
- **Fields**: no abbreviations — `perm` → `permutation`.

### Which types get a `*Dsl` (pragmatic policy)

A `*Dsl` exists **iff the type crosses the EDN/string boundary as an independent
payload** — a standalone entity, a constraint, or a constraint *value* — or is
broadly reused (`ValueAst`). That covers, today: the entity `*Dsl`s (`Atom`,
`Bond`, `Molecule`, `Dative`, `Multicenter`, `Noncovalent`, `AromaticSystem`,
`StereoAtom`, `StereoBond`); the constraint `*Dsl`s (`AtomConstraint`,
`BondConstraint`, `Constraints`, `SubPatternAnchor`); and the value/constraint-value
`*Dsl`s (`Value`, `AromaticValence`, `MulticenterValence`, `StereoConfiguration`,
`StereoCoset`).

**Entity field components that only render *inside* their entity's surface string
get no `*Dsl`** — they round-trip through the entity's `*Dsl` plus a shared
sub-grammar parser. These are exactly: **`ElementAst`, `IsotopeMassAst`,
`SpinStateAst`, `NoncovalentBondKindAst`**. *Why:* they never cross the boundary
standalone, so a struct would be pure boilerplate; the rule's real target — a
component **bypassing** the boundary into a bare, un-round-trippable AST (the
permutation case) — doesn't apply, since these round-trip via their entity's `*Dsl`.

## Equality strategy: precedents and the eager-vs-lazy trade-off

The per-type canonicalization logic (the step-lists above) is **strategy-independent**.
What is open is *when* a value acquires its normal form for the sake of `==`/`Hash`/`Ord`.
Real systems split on the **construct:compare ratio**, not on a universal best:

**Eager — canonicalize / intern on construction.** `==`/`Hash` become trivial (often
pointer equality). Chosen where equality/sharing *is* the workload:
- reduced-ordered **BDDs** — the unique table (hash-consing) makes canonical-by-construction
  the entire point;
- **rustc** interns `Ty`, consts, predicates in `TyCtxt` → pointer equality, because type
  equality is on the hot path;
- **CAS** — SymPy auto-evaluates on build (with `evaluate=False` as an escape *because*
  eager eval is a known cost); Mathematica rewrites to a fixpoint;
- **LLVM** `IRBuilder` constant-folds — but as a *convenience*; raw emission stays available.

**Plain tree + normalize-in-passes + equality-on-demand.** Build freely; normalize in
explicit passes; compare via value-numbering at the point of need. Chosen where you mostly
*traverse* and rarely compare:
- **rustc's surface AST** is a plain tree — not interned, not canonicalized (you walk it,
  you don't `==` it). The same project interns `Ty` — opposite choices in one codebase,
  driven purely by how often equality is invoked.
- **LLVM IR** — GVN/CSE are passes; equality is a value-numbering hash table, not a stored
  normal form.

**Quotient `==`/`Hash` (canonicalize-on-compare).** Normalize transiently inside
`==`/`Hash` (hash on the canonical form), store raw — doc 095's **"simplified-structural"**
model. Cost is one `canonicalize` per comparison; optionally **memoize** the normal form /
hash so it is paid once per value.

**Equivalence classes (no normal form).** e-graphs (egg) keep equality by union-find +
congruence closure incrementally — neither eager nor compare-time; "equal" = "same
e-class." Not applicable here (we have a normal form); listed for completeness.

**The trade-off for umol.** `Undetermined`/`Lit` dominate, and their `canonicalize` is an
O(1) no-op — so on the hot path eager and lazy are *both* trivial. They diverge only on the
rare complex values (sets, expressions), which chemistry keeps small (a charge, a short
valence expression, a few-element set), so the compare-time cost is negligible. Pattern
matching goes through `meet` (which canonicalizes regardless of this choice); standalone
`==`/`Hash` is mostly map keys and tests; molecule-scale dedup is **graph** canonical form
(WL/rank), separate from per-leaf AST equality. So the leaf-AST construct:compare ratio
points at compare-time.

And because the AST is a **transparent data carrier with no facade** (raw variant
construction stays public), eager canonical-by-construction is **not enforceable** — a
smart constructor is a convenience, never the sole entrypoint — so eager could never
*guarantee* "structural `==` is semantic" anyway. The honest model is then: *nothing*
assumes a stored value is canonical; the operations that need a normal form produce one on
demand.

**Status: open, leaning lazy (095 model 2 — "simplified-structural").** No decision
recorded. The choice affects only:
- Decision A's wording — eager "canonical-by-construction; structural `==` is semantic" vs
  lazy "raw construction; **custom** `Eq`/`Hash`/`Ord` canonicalize-then-compare";
- whether `Eq`/`Hash`/`Ord` are `#[derive]`d (eager) or hand-written (lazy, keying on the
  canonical form for `Eq`-consistency);
- whether `meet`/`join` may assume canonical inputs (eager) or must canonicalize inputs
  defensively (lazy).

The step-lists, the lattice algebra, the boundary types, and the per-type variants are the
same under either choice. Still-live from 095 before committing: its open Q2 — inventory of
where structural `==`/`Hash` leaks (`MoleculeAst::PartialEq`, `HashMap`/`HashSet` keys,
constraint dedup, any AST-keyed cache).

## Faithful parse and round-trip

The parser is **faithful**: it builds raw values preserving surface structure, with **no
folding** (no `simplify` in parsers). Canonicalization is **representation-normalization**,
not semantic `simplify` — it only chooses one normal form among *surface-equivalent*
inputs (`{2,0,1} ≡ {0,1,2}`, `#h2 ≡ #h(2)`, `~~x ≡ x`) and never changes the admissible
value set, so it is compatible with the no-`simplify`-in-parser rule whenever it runs.

Round-trip holds independent of the eager/lazy choice because **canonical rendering goes
through the boundary type**: the `*Dsl` emits the canonical surface. `#h(2)` and `#h2` both
render to canonical `#h2` and re-parse equal — under eager because the stored value is
already canonical, under lazy because the renderer canonicalizes at the boundary. A raw
stored value never escapes as a non-canonical *rendering*; whether it is stored canonical
is the strategy question above. This normalizes only surface spelling, not meaningful tags,
so it preserves roundtrip fidelity (no tag erasure).

## Why graph-shaped ⇒ not a lattice (order theory)

Order finite graphs by subgraph embedding (`G ≤ H` iff `G` is isomorphic to a
subgraph of `H`) — exactly substructure / subsumption matching. This is a poset
but **not a lattice**: meets and joins are not unique.

- No join: `A = 2K₂` (two disjoint edges) and `B = P₃` (2-edge path) have two
  incomparable minimal common supergraphs — `P₄` and `P₃ ⊔ K₂` — so `A ∨ B` does
  not exist.
- No meet: the maximum common subgraph is famously non-unique (multiple
  incomparable MCSs of equal size).

So molecule `matches`/`meet`/`join` are inherently search (alignment-dependent,
non-unique results), not algebraic `∧`/`∨`. (Aside: the *homomorphism* order on
graphs *is* a distributive lattice — categorical product / disjoint union — but
that is not subgraph embedding, so not what substructure matching uses.)

## Verification

The C4e.5 lattice sweep is the regression target: every `Lattice` type green
(including a raised `PROPTEST_CASES` pass), plus the existing umol-ast suites and
the molecule round-trip proptests. `AromaticSystemAst`/`MulticenterBondAst`
re-enter the sweep cleanly once `electrons` is retyped to `ElectronCountsAst`; the
concrete predicates (`matches`-only) get their own targeted `matches` tests.
