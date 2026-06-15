# 113 · AST canonical equality and lattice algebra

Status: Active · 2026-06-14
Resolves: 095 open questions 1–3 (equality model, leak inventory, canonical form).

## Decisions

### A. Equality: canonical-by-construction

AST value types move to a canonical representation so structural `==`/`Hash`
equal semantic equality:

1. **Set-typed storage** for disjunctive/negated variants (`Set`, `NotSet`) —
   `BTreeSet` not `Vec`. Ordering canonicalization disappears into the type;
   `meet` = intersection, `join` = union.
2. **Canonical-on-construction** for the decidable algebra (set order, dedup,
   singleton collapse, double negation). Smart constructors maintain the normal
   form; a non-canonical value is unconstructible.
3. **`Expr` lowered eagerly where the algebra is decidable** (coset `Expr` →
   apply the permutation → `Lit`, surviving only with a free `?var`; `ValueExpr`
   → fold its decidable fragment: const-fold, commutative operand sort, identity
   elimination). Residual symbolic `Expr` is **opaque**: structurally equal ⇒
   equal (sound); structurally different ⇒ conservatively incompatible. This is
   exact on the folded fragment, sound-but-incomplete beyond it (the trait's
   existing "modulo canonicalization" hedge).

`simplify` becomes the construction-time canonicalizer (an invariant), not a
compare-time step callers must remember.

### B. `Lattice` is for value/constraint atoms of information only

4. **`matches` derived from `meet`** — trait default `matches(t) == (meet(t) ==
   Some(t))`. Per-type hand-written `matches` is removed except where a genuine
   cheaper-but-equivalent override exists. This kills the class of bugs where a
   hand-written `matches` disagrees with `meet` (e.g. the `Expr` non-reflexivity).
5. **`join` stays `Self`** (not `Option`). Every type that remains a `Lattice`
   has a faithful top (`Undetermined` / empty collection / full relation-set),
   so `join` is total and sound. The `meet: Option`, `join: Self` asymmetry is
   correct for a bounded lattice with a representable top and an absent bottom.

### C. Identity is not a lattice axis; graphs are not lattices

The reason `join` looked partial (the `debug_assert` / unsound `Self::default()`
collapses) was a category error: a type's **structural identity** was modelled as
a lattice field. The fix is to take identity *out* of the lattice — not to drop
the lattice:

6. **`AromaticSystemAst` / `MulticenterBondAst` keep `Lattice`.** The per-atom
   `electrons` vector (whose length is the size index) is **electron-count
   information = a constraint**, not a lattice axis. Modelled in `constraints`,
   with the vector gone the struct is `(charge, spin, constraints)` and
   `#[derive(Lattice)]`s cleanly. (Sub-question: the exact constraint form —
   total `ElectronCount` vs per-atom-keyed.)
7. **`DativeBondAst`** keeps its attribute fields; the birelation promotion
   removes only `acceptor_slot`. With no index left it is an index-free
   named-field struct and switches from its hand-written impl to
   `#[derive(Lattice)]` (field-wise, like `BondAst`).
8. **`MoleculeAst` / `ReactionRuleAst` are never lattices.** Their order is
   substructure (subgraph embedding), which is intrinsically *not* a lattice
   (next section) — not fixable by relocating a field. Molecule-level
   `matches`/`meet`/`join` are correspondence-based graph algorithms (subgraph
   match; MCS / anti-unification — doc 108), delegating to per-node attribute
   `matches`.
9. **`StereoAtomAst` / `StereoBondAst` keep `Lattice`** via a kind set-lattice.
   `kind` becomes a finite-domain set-lattice `StereoKindAst` (`Undetermined` =
   all kinds | `Lit` | `Set`), idiomatically the `relation_ast!` shape. The
   element is then `{ kind, coset, constraints }`, field-wise `Lattice` with top
   `(Undetermined, Undetermined, ∅)` = "any stereo atom". This works because
   `StereoCosetAst`'s `meet`/`join` are kind-free (`kind` is used only by
   `simplify` / `apply_permutation` / `matches_value`). The `post_meet` and
   `post_join` hooks both **normalize `coset → Undetermined` whenever `kind` is not
   a single `Lit`**, so a concrete coset never coexists with a non-concrete kind
   (chosen over accepting strange `({TH,SP}, Lit(0))` states). Needing the
   normalization after `join` as well as `meet` is exactly why the single
   `saturate` is split into the two hooks. Chemistry rules (kind ↔ ligand count,
   kind ↔ element) are validator concerns, not lattice ones —
   AST-valid-but-molecule-invalid states are accepted, as `He` at a two-neighbour
   site is.
10. **`TopicityAst` is `matches`-only.** Topicity is intrinsically a relation
    *between a specific pair of ligands* — there is no generic, pair-free
    topicity, and none at all for non-binary stereo kinds (many ligand pairs). So
    `pair` is essential identity: not removable, and not a liftable set-lattice
    dimension the way `kind` is (a single relation over a *set* of pairs is
    meaningless). `TopicityAst { pair, rel }` therefore has no global top and is
    not a single lattice — it carries an inherent `matches`; the lattice is the
    per-pair `rel` (`TopicityRelationAst`), combined only between entries with the
    same `pair`. This is the structural difference from stereogenicity (a
    per-site yes/no, which flattens to one `StereogenicityAst` lattice).

Settled non-`Lattice` types: the graph-shaped ones (`MoleculeAst`,
`ReactionRuleAst`), the concrete predicates (`LigandPermutation`,
`OrientedLigandPermutation`, `LigandSymmetry`, `Fluxionality` — documented
"not a lattice; no top, no `join`"), the pure key `LigandPair`, and `TopicityAst`
(per-pair; point 10).

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

Entities are different — their non-lattice symptom was self-inflicted, not
intrinsic. An entity's position-indexed vector (`AromaticSystemAst.electrons`)
ties each cell to a specific member atom, so size/membership is **identity, not a
lattice axis**; a size-`Undetermined` top is definable (the count is an integer)
but semantically inert. The fix is to remove identity from the lattice
definition (C.6 — electron count becomes a constraint), after which the entity is
a clean field-wise lattice. A discrete identity can instead be lifted to a
finite-domain set-lattice *dimension* (C.9 — `kind` → `StereoKindAst`), kept
consistent by a `post_*` hook. Either route is *not* available to `MoleculeAst`:
its non-lattice-ness is the subgraph order itself, not a stray field.

## Empirical: who drives lattice ops today

Only the **atom-level valence resolver**: `atom.ast.meet(&derived)`
(`ops/invariants.rs`, `ops/valence/counts.rs`), `narrow_from`
(`ops/valence/atom_typing.rs`). `AtomAst` itself uses an **inherent** `meet` and
does not `impl Lattice`. Entity-level `meet`/`join` (aromatic / multicenter /
dative) and any molecule-level combination have **no production caller** — so the
(C) changes to these lattices touch speculative, production-untested code, not
live behaviour.

## Per-type review

`Lattice` is supplied two ways: **`#[derive(Lattice)]`** (`umol-ast-macros`,
field-wise over named struct fields; all-field-top is the top) and
**hand-written** impls (enums, tuple structs, collections). Two optional post-op
cross-field hooks replace the old single `saturate`: **`post_meet(&mut self) ->
Result<(), Contradiction>`** (fallible — `meet` narrows, can hit ⊥) and
**`post_join(&mut self)`** (infallible — `join` widens). The derive runs each
after its op; `#[lattice(post_meet = "fn", post_join = "fn")]` wires them.
`AtomAst` uses `post_meet` for JointDomain propagation; stereo elements use both
for the coset/kind collapse; all other types take the no-op defaults.
Disposition codes: `derive` / `hand` keep `Lattice`; `matches` = inherent
`matches`, no trait; `—` = neither (graph-mediated or pure key).

### Entity ASTs

| Type | Lattice | Index → per-index ASTs | post-op hook | On-construction simplification |
|---|---|---|---|---|
| `AtomAst` | derive | — | `post_meet` (JointDomain) | field-wise; each field canonical |
| `BondAst` | derive | — | — | field-wise |
| `SpinStateAst` | derive | — | — | field-wise (`unpaired`, `multiplicity`) |
| `NoncovalentBondAst` | derive | — | — | field-wise |
| `DativeBondAst` | derive *(after birelation drops `acceptor_slot`)* | — | — | field-wise |
| `AromaticSystemAst` | derive *(after `electrons`→constraint)* | membership (external, not a field) | — | electron count → constraint; `(charge, spin, constraints)` field-wise |
| `MulticenterBondAst` | derive *(after `electrons`→constraint)* | membership (external) | — | as above |
| `StereoAtomAst`, `StereoBondAst` | derive | `kind` is a set-lattice field (`StereoKindAst`); anchor = store key | `post_meet`/`post_join`: coset→`Undetermined` if kind ≠ `Lit` | field-wise `(kind, coset, constraints)` |
| `MoleculeAst` | — (graph) | — | — | graph canonicalization (WL / canonical rank), not an AST normal form |
| `ReactionRuleAst` | — (graph) | — | — | graph-level |

### Leaf / predicate ASTs

| Type | Lattice | Index | On-construction simplification |
|---|---|---|---|
| `ValueAst` | hand (enum) | — | `Set` sort/dedup/singleton→`Lit`; `Expr` fold decidable→`Lit`, opaque residual |
| `ValueExpr` | *(in `ValueAst`)* | — | const-fold, commutative-operand sort, identity elim, flatten assoc |
| `ElementAst` | hand (enum) | — | `Set`/`NotSet` sort/dedup/collapse; negation polarity (double-neg/De Morgan); ∅→⊥ |
| `IsotopeMassAst` | hand (enum) | — | as `ElementAst`; `Natural` separate channel |
| `NoncovalentBondKindAst` | hand (enum) | — | none (flat) |
| `StereoKindAst` | hand (`relation_ast!`) | — | finite-domain set: sort/dedup/collapse + negation polarity |
| `StereoConfigurationAst` | hand (enum) | — | delegate to coset; `NotStereo` ≠ `Stereo` |
| `StereoCosetAst` | hand (enum) | — | `Expr`→apply perm (decidable)→`Lit` unless free `?var`; `LitSet` sort/dedup/collapse |
| `StereoExpr` | *(in `StereoCosetAst`)* | — | operator normal form (`~` involution, `'` mirror, `^` apply)→`Lit`/`LitSet` |
| `AromaticValenceAst` | hand (enum) | — | delegate `ValueAst`; `NotAromatic` ≠ `Aromatic(_)` |
| `MulticenterValenceAst` | hand (enum) | — | delegate `ValueAst`; `NotMulticenter` ≠ `Multicenter(_)` |
| `JointDomainAst` | hand (enum) | — | domain set canonical |
| `TopicityRelationAst` | hand (`relation_ast!`) | — | finite-domain `Set`/`NotSet` sort/dedup/collapse + negation polarity |
| `StereogenicityAst` | hand (`relation_ast!`) | — | as above; flattened — `relation_ast! { StereogenicityAst, Stereogenicity }`, no `*RelationAst`, no wrapper; macro also generates `AsLit` (`Lit`→`Some`) |
| `TopicityAst` | **matches** | `pair` (essential — topicity is per-pair) | per-pair `rel` (`TopicityRelationAst`) is the lattice |
| `LigandPermutation` | matches (= equality) | — | none; EDN via `LigandPermutationDsl` |
| `OrientedLigandPermutation` | matches | — | none; EDN via `OrientedLigandPermutationDsl`; field `permutation` |
| `LigandSymmetry` | matches | — | none; field `permutation` |
| `Fluxionality` | matches | — | none; field `permutation` |
| `LigandPair` | — (key) | *is* an index | normalize `first ≤ second`; EDN via `LigandPairDsl` |

### Constraint ASTs

The per-constraint **enums** (`AtomConstraint`, …) are not themselves `Lattice`;
the **collections** are. None override the post-op hooks.

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
| `RelationalConstraint` (JointDomain) | — | relational; folded via `JointDomainAst` |

`StereogenicityAst` flattens to `relation_ast! { StereogenicityAst, Stereogenicity }`
directly — a per-site yes/no, no `pair`, so no wrapper and no `*RelationAst`.
`TopicityAst` does not flatten: topicity is intrinsically per-pair (no generic
topicity, none for non-binary kinds), so `pair` is essential and `TopicityAst {
pair, rel }` stays a `matches`-only carrier with `TopicityRelationAst` as the
per-pair lattice (decision 10).

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

## How canonical equality is typically achieved (survey)

- **Canonicalize-on-construction / hash-consing** — reduced-ordered BDDs, CAS
  automatic simplification, LLVM constant folding. `==`/`Hash` are structural and
  semantic; works for **decidable** domains, not open symbolic equivalence.
- **Faithful tree + separate `simplify`** — lossless parse, structural `==` (the
  current footgun).
- **Surface/core split** — faithful surface, canonical core via a lowering pass.
- **Quotient `==`/`Hash`** — normalize inside compare/hash; per-op cost.

Chosen: canonical core for the decidable classes; eager lowering + opaque
residual for the symbolic remainder.

## Faithful-parse reconciliation

The canonicalized forms carry **no meaningful surface information**: `{2,0,1}` and
`{0,1,2}` denote the same set; `#h2` and `#h(2)` the same value; double negation
is identity. Canonicalizing them is representation-normalization, not semantic
`simplify`, so it respects the no-`simplify`-in-parsers rule (the parser
constructs faithfully via the smart constructors, which only pick a normal form
for equivalent inputs; class-c lowering is a separate pass). Round-trip holds
because a non-canonical value never exists to be rendered.

## Scope

Foundational: the `Lattice` trait (default `matches`; `saturate` → `post_meet` /
`post_join`) and its derive macro, the set-variant storage, the
equality/`Hash`/`Ord` derives, the (C) changes (electron count → constraint;
`kind` → `StereoKindAst` set-lattice; `pair` unresolved), the naming/boundary
cleanup (literal renames + `*Dsl` boundary types replacing the inline stereo-perm
serde; `perm` → `permutation`), and the construction sites that build sets. 095's
leak inventory is the re-audit surface: `MoleculeAst::PartialEq` (proptest
round-trip), `HashMap`/`HashSet` keys, constraint dedup, alias bijectivity
(`BiBTreeMap`). Expect an extended red period.

## Open implementation questions

To resolve before an implementation plan (no plan/code yet):

1. **Enforcement of the construction invariant** — private-field newtypes
   (non-canonical unconstructible; every pattern site changes) vs smart
   constructors by convention (lighter, leakable).
2. **Negation canonical form** — pick a canonical polarity for `Not`/`NotSet`
   over the finite element universe (e.g. by cardinality) vs keep the written
   polarity. Settles double negation / De Morgan uniformly.
3. **`Expr` lowering site** — dedicated post-parse pass vs in `simplify` vs eager
   in smart constructors, under the no-simplify-in-parser constraint.
4. **Inherently non-ground patterns** (`Bind`, `Ref`, free-`?var` `Expr`) —
   confirm opaque structural `meet`/`matches` semantics.
5. **`matches` overrides** — which lattices keep a cheaper hand-written `matches`
   vs take the `meet`-derived default.
6. **Hash-consing / interning** as an optional dedup layer atop canonical values
   — in or out of scope.
7. **Migration sequencing** across the red period (storage → constructors →
   trait default → demotions → test-expectation canonicalization).
8. **Simplification framework.** Canonicalization is scattered today —
   `simplify` / `simplify_values` / `simplify_each` (`ValueAst`, `StereoExpr`,
   entities, collections), normalize-on-construction (`LigandPair`), the
   field-wise `#[derive(Lattice)]`, and the proposed smart constructors. Decide
   the single locus and mechanism (construction-time invariant via smart
   constructors vs a `simplify`/canonicalize trait vs a derive), how it composes
   for nested types (field-canonical ⇒ parent-canonical), and how it relates to
   the `Lattice` derive.
9. **`matches_value` consistency.** Present on `ValueAst` (`i64`),
   `Stereo{Configuration,Coset,Expr}` (`u32` + `kind`),
   `{Aromatic,Multicenter}ValenceAst` (`i64`); absent on `ElementAst` /
   `IsotopeMassAst` (which match via element sets). Decide a uniform
   "pattern admits a ground scalar" predicate, parameterized by the type's
   ground / `AsLit` type, vs keeping it ad-hoc; relate to deriving `matches`
   from `meet`.

## Verification

The C4e.5 lattice sweep is the regression target: every `Lattice` type green
(including a raised `PROPTEST_CASES` pass), plus the existing umol-ast suites and
the molecule round-trip proptests. `AromaticSystemAst`/`MulticenterBondAst`
re-enter the sweep cleanly once `electrons` moves to `constraints`; the
concrete predicates (`matches`-only) get their own targeted `matches` tests.
