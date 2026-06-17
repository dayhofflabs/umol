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

### P1 · Leaf value types (bottom-up)
1. **`ValueAst`** **Done** - split `Expr`→`ValueTerm`(i64)/`ValuePredicate`(bool); `Bind`/`Ref`→
   `Var`; `Set`→`LitSet(BTreeSet)`; the full canonicalize step list (term fold, predicate
   NNF, lifts). Parser: build `Sum`/`Product`/`Div`/`Rem`, drop the `unary_expr` sign-XOR;
   delete `ValueAst::simplify`/`ValueExpr::simplify`. The biggest single change (all of
   `dsl/value.rs` + every `ValueAst` consumer).
2. **`ElementAst`** — `Undetermined|Lit|LitSet|NotSet|Var`; cardinality-canonical; drop `Not`.
3. **`IsotopeMassAst`** — `Undetermined|Natural|Lit|LitSet|Var`; positive-only; `u32`.
4. **`NoncovalentBondKindAst`**, **`StereoKindAst`** — `Undetermined|Lit`; identity
   `canonicalize`; `canonical` always borrows.
5. **`ElectronCountsAst`** (new) — `Undetermined|Lit(Vec<u8>)`; identity (positional).

### P2 · Predicate / relation types
- **`SpinStateAst`** — hand `canonicalize` with the `are_compatible` parity gate;
  `From<(u8,u8)>`/`From<SpinState>` → `TryFrom`.
- **`AromaticValenceAst`/`MulticenterValenceAst`** — delegate inner `ValueAst`; drop the
  hand `matches`.
- **`TopicityRelationAst`/`StereogenicityAst`** — `relation_ast!` over `BTreeSet`; flatten
  stereogenicity (drop `StereogenicityRelationAst` + the wrapper).
- **`TopicityAst`** — drop `impl Lattice` → matches-only `{pair, rel}`; remove its
  fixed-pair lattice proptest (covered by `test_topicity_relation_ast_lattice_laws`).

### P3 · Stereo configuration cluster
- `StereoExpr`→**`StereoTerm`** (`Var(name, opt domain)` | `Swap` | `Mirror` | `Apply`; no
  `Lit`/`LitSet`). `StereoCosetAst = Undetermined|Lit|LitSet|NotSet|Term` — **no**
  `Lattice`/`Canonicalize` (kind-relative).
- **`StereoConfigurationAst { kind: StereoKindAst, coset }`** (element side) +
  **`StereoSiteAst = Undetermined|NotStereo|Stereo(StereoKind, coset)`** (constraint side;
  renames the old tristate). Both: `Canonicalize`+`Lattice` reading `self`'s kind, the
  `are_compatible(kind, coset)` gate, the Term compose→priority(`Mirror>Swap>Apply`)
  normal form; `AsLit = StereoConfiguration { kind, coset }`. Update `dsl/stereo.rs`, the
  stereo views, and `StereoAtomAst`/`StereoBondAst` to carry the joint config.

### P4 · Entities (pure field-wise)
- `AtomAst`, `BondAst`, `NoncovalentBondAst` — `#[derive(Lattice, Canonicalize)]`; delete
  `simplify_values`. Fix `NoncovalentBondAst` direct `Display`/`FromStr`/`FromEdn` (move to
  `NoncovalentBondDsl`).
- `AromaticSystemAst`/`MulticenterBondAst` — `electrons: Vec<ValueAst>` → `ElectronCountsAst`;
  retype ctors; derive.
- `DativeBondAst` — derive after the birelation `acceptor_slot` drop.

### P5 · Constraint types
- `AtomConstraint` & siblings — `Canonicalize` = delegate inner (replaces `simplify`).
- Collections (`AtomConstraints` …) — kind-keyed `Lattice`+`Canonicalize`: fixed kind
  order, dedup, **drop vacuous (`Undetermined`-payload) entries**, multi-valued
  (`RingSize`) set-canonical. Delete `simplify_each`.
- Molecule-level `Constraints`/`Constraint` — `Canonicalize` (flatten/sort/dedup the
  `And`/`Or`/`Not` tree, ID-bearing); **not** a `Lattice`.

### P6 · `Lattice`-trait flip + macro (lands once P1–P5 all impl `Canonicalize`)
- `Lattice: Canonicalize`; `matches` becomes the `meet`-derived default; `join` stays
  `Self`. `#[derive(Lattice)]` generates `meet`/`join` = field-wise + `canonicalize`.
  Remove the now-redundant hand-written `matches` impls (keep only genuine cheaper
  overrides). Hand-written leaf `meet`/`join` already canonicalize from P1–P5.

### P7 · Semantic-equality adoption + boundary cleanup
- 095 Q2 leak audit → route to `equiv`/`Canonical<T>` where semantic keys are needed:
  `MoleculeAst::PartialEq` (graph-canonical), `HashMap`/`HashSet` AST keys, alias
  `BiBTreeMap`, constraint dedup. Decide structural-vs-semantic per site.
- Drop `matches_value` (≡ `matches(lift(v))`); relocate `capture`/`evaluate` to the
  resolver or delete.
- Literal renames + `*Dsl` per [[feedback_dsl_boundary_types]]
  (`StereoConfigurationDsl`→`StereoSiteDsl`; add `TopicityDsl`/`StereogenicityDsl`).

### P8 · Verification (doc 111)
- C4e.5(1): all retained `lattice::` tests green (raised `PROPTEST_CASES`); demoted types
  leave the sweep.
- C4e.5(2): atom-DSL roundtrip; C4e.5(3): `canonicalize` idempotence beyond `ValueAst`.
- `umol-ast` lib + `--features proptest` + workspace build + conformance all green.

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
| `SpinStateAst` | derive (field-wise `unpaired`, `multiplicity`) | — | cross-field parity: ground pair satisfies `are_compatible` else `Err` |
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
| `RingCount`, `RingSize` value | `BondConstraint` | ring constraints |

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
succeeds only when both fields are `Lit` (parity now guaranteed by canonicalize, so the
old parity re-check there is redundant).

**Boundary.** No `*Dsl` (pragmatic policy — field component of `AtomAst`/`BondAst`).
Round-trips inline via the shared `#u`/`#s` predicates (`apply_spin_pair` /
`fmt_spin_pair`, `dsl/predicates.rs`). Surface `#u2`, `#s3`, `#u2#s3`. (`SpinState`'s
own `Display`/`FromStr` in `umol-shared` is the *literal's* surface — a ground domain
type, not an `*Ast` — so not the serde-on-AST violation.)

**Canonicalization** (hand-written `Canonicalize`, like `AtomAst`):
1. field-wise canonicalize `unpaired`, `multiplicity`.
2. **cross-field parity:** both `Lit` and `!are_compatible(u, m)` → `Err(Contradiction)`
   (enforce the physical invariant on construction; matches the practice elsewhere). No
   cross-pruning of non-ground sets. **[new]**

Parallel removals: drop `simplify_values` (`spin.rs:27`, subsumed); update
`is_plus_sugar` (`predicates.rs:83`) — it still matches the old `ValueAst::Expr` /
`ValueExpr`, now `Term`/`Predicate`.

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

**Cross-field predicate `are_compatible(kind, coset)`** (parallel to
`SpinStateAst::are_compatible`): the kind↔coset-index validity — when `kind = Lit(k)`,
every coset index (in `Lit`/`LitSet`/`NotSet` and any `Var` domain) lies in
`[0, k.count())`. `canonicalize` gates on it (out-of-range → `Err`), and `canonical()`
borrows only when it holds. `StereoSiteAst::Stereo(kind, coset)` uses the same predicate.

**`StereoConfigurationAst::canonicalize` — full step list.**

Field-wise + kind gate:
1. canonicalize `kind` (`StereoKindAst`: identity).
2. **kind gate:** `kind = Undetermined` ⇒ set `coset = Undetermined`, return `Ok` (no
   concrete configuration without a concrete geometry — `SpinStateAst`-style hook).

Otherwise `kind = Lit(k)`, `n = k.count()`; fold `coset` under `k`:

`coset` = literal / set forms:
3. `Lit(i)`: `i ≥ n` → `Err(Contradiction)` (kind ↔ coset-index).
4. `LitSet(s)`: range-check `s ⊆ [0,n)` (else `Err`); `s = ∅` → `Err`; `|s| = 1` →
   `Lit`; `s` = full → `Undetermined`; else sorted/deduped `BTreeSet`.
5. `NotSet(s)`: range-check; `s = ∅` → `Undetermined`; `s` = full → `Err`; cardinality
   polarity vs `n` — keep `NotSet(s)` iff `|s| < n − |s|`, else `LitSet([0,n) \ s)`
   (tiebreak positive), then re-apply step 4's collapse.
6. `Undetermined`: identity.

`coset` = `Term(t)` — operators are all permutation actions, so:
7. compose the operator word into **one net permutation `g_total`** over the inner `Var`
   (`Mirror` = μ_k, `Swap` = ι_k, `Apply(g)` = g; identity factors drop).
8. canonicalize the inner `Var(name, dom)`: if `dom = Some(s)`, range-check `s ⊆ [0,n)`
   (else `Err`), sort/dedup, `s = ∅` → `Err`, `s` = full → `None`. `name` always kept;
   a `Var`-rooted term never folds to a `Lit`.
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
inner value's `canonicalize` (replacing `AtomConstraint::simplify`). Each is a
boundary-independent payload → keeps its `*Dsl` (`AtomConstraintDsl`, …).

**Per-entity collections** (`AtomConstraints`, `BondConstraints`,
`DativeBondConstraints`, `MulticenterBondConstraints`, `AromaticSystemConstraints`,
`StereoAtomConstraints`, `StereoBondConstraints`; `NoncovalentBondConstraints` trivial —
inner enum uninhabited) are `Lattice` + `Canonicalize`:
- **Lattice** (hand): kind-keyed field-wise. Each unique kind is read via a typed
  accessor that returns `Undetermined` when absent; `meet` per kind (`None` if any kind
  contradicts), `join` per kind. The multi-valued kind (`RingSize`) meets by
  set-union + dedup. Top = empty.
- **Canonicalize**: fixed kind order; dedup unique kinds (last-wins); **drop entries
  whose value is `Undetermined`**; canonicalize each inner value; multi-valued kinds
  sorted + deduped.

**Key step — drop-vacuous moves into `canonicalize`.** Today a vacuous
(`payload = Undetermined`) constraint may sit in the AST and is only elided at *render*
time (`dsl/predicates.rs` canonical-rendering note); `meet` already refuses to add them.
Canonical-by-construction *requires* dropping them in `canonicalize` so that
`{Valence(Undetermined)}` and `{}` are structurally equal. Consequence: a parsed
`#v*` (valence-undetermined) normalizes to no entry — faithful (same meaning, like
`1+1 → 2`), AST-level elision rather than just surface. **Accepted** — it only drops
information-free tops, so roundtrip fidelity is preserved.

(`JointDomain` is removed — see Scope; cross-field/cross-atom correlations move to the
molecule-level variable-constraint facility, doc 115.)
`AtomConstraint::TetrahedralStereo`'s payload retypes `StereoConfigurationAst` →
`StereoSiteAst` (the rename). **Per-kind constraint variants stay** (`TetrahedralStereo`,
`CisTransStereo`): `#T` is tetrahedral-specific, `#C` cis/trans-specific. The concrete
`StereoKind` duplicated in `StereoSiteAst::Stereo` is accepted.

**Boundary.** Collections render inline via the owning entity's `*Dsl`; the enums have
`*Dsl`. **Parallel removals:** `simplify_each` on every collection;
`<Enum>::simplify` → `Canonicalize`.

**Molecule-level `Constraints` / `Constraint`.** `Constraint = And(Vec) | Or(Vec) |
Not(Box) | Atom(AtomId, AtomConstraint) | Bond(BondId, BondConstraint) | <molecule-scope
predicates>` — a boolean combinator tree over **ID-scoped** predicates; `Constraints` is
a flat (conjunctive) `Vec<Constraint>`. **Not a `Lattice`** (ID-bearing; the molecule
order is graph subsumption, not algebraic). `Canonicalize` mirrors `ValuePredicate`:
recurse into children, flatten nested same-combinator, sort + dedup, drop empty
`And`/`Or`; inner predicates canonical. Equality is structural **with** the IDs.
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
