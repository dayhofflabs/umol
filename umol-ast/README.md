# umol-ast

Semantic AST (abstract syntax tree) and surface DSL (domain-specific language) for molecule data in umol.

Two layers stacked over umol-graph-core:

| Layer | Module | Purpose |
|---|---|---|
| AST | `crate::ast` | semantic, `MoleculeAst` / `AtomAst` / `BondAst` / …; ground or pattern; algebraic |
| DSL | `crate::dsl` | surface, `AtomDsl` / `BondDsl` / …; EDN and string parsing; config-driven defaulting |

See `spec/umol-dsl-spec.md` for the normative DSL grammar. This README documents the AST invariants and the AST↔DSL boundary.

## Entity kinds

`MoleculeAst` holds six relation kinds over the atom set:

| Relation | Storage | Per-entity AST |
|---|---|---|
| atoms | positional `Vec<AtomAst>` | `AtomAst` |
| covalent bonds | `Graph` edges + parallel `Vec<BondAst>` | `BondAst` |
| dative bonds | `FixedRelationSet<DativeBondAst, 2>` | `DativeBondAst` |
| aromatic systems | `VarRelationSet<AromaticSystemAst>` | `AromaticSystemAst` |
| multicenter bonds | `VarRelationSet<MulticenterBondAst>` | `MulticenterBondAst` |
| noncovalent bonds | `FixedRelationSet<NoncovalentBondAst, 2>` | `NoncovalentBondAst` |

Index types (`AtomIdx`, `BondIdx`, `DativeBondIdx`, `AromaticSystemIdx`, `MulticenterBondIdx`, `NoncovalentBondIdx`) are disjoint newtype wrappers over `NodeId` / `EdgeId` / `RelationId`.

## AST is pure data

`MoleculeAst` and all per-entity AST types carry no caches, no derived views, no interior mutability. Equality and hashing are structural. Derived views (rings, distance matrix, match targets) live on runtime cells (`Molecule`, `Pattern`) outside this crate.

Storage is `Arc`-shared (`Arc<Vec<AtomAst>>`, etc.); cloning a `MoleculeAst` is an `Arc` bump. Structural edits go through `MoleculeBuilder` (`MoleculeAst::edit`); attribute edits go through per-entity mutation accessors (`atom_mut`, `bond_mut`, …) that use `Arc::make_mut` semantics.

## Grounding

Each AST form has a fixed set of **inherent fields** (the identity-bearing slots) and a `constraints: Vec<_>` tail. `is_ground` is `true` iff every inherent field holds a single concrete value (literal; not `Undetermined`, set, bind, ref, or symbolic state). Constraints — including `Undetermined`-valued ones — do not affect grounding.

| Form | Inherent fields |
|---|---|
| `AtomAst` | element, isotope mass, charge, implicit hydrogens, lone pairs, spin (unpaired, multiplicity) |
| `BondAst` | order, charge, spin |
| `AromaticSystemAst` | charge, spin, π-electron count |
| `MulticenterBondAst` | charge, spin, electron count |
| `DativeBondAst` | direction (order is fixed at two electrons by definition) |
| `NoncovalentBondAst` | interaction kind |

`MoleculeAst::is_ground` is the conjunction over every entity AST.

## `Undetermined` is a wildcard

Every inherent-field value type has an `Undetermined` state. In `matches(&self, target)`:

- `(Undetermined, _) => true` — wildcard in the pattern admits any target.
- `(_, Undetermined) => false` — a pattern with a concrete value does not match an undetermined target.

This rule is enforced in every `*Ast::matches` implementation; consumers of AST match semantics rely on it. `Undetermined` is never equality-matched against itself — it is only a wildcard in the pattern position.

## Dual constraint storage

Constraints are stored in **two places**, intentionally, and consumers read the union:

1. **Inline on each entity** — `AtomAst.constraints: Vec<AtomConstraint>`, `BondAst.constraints`, etc. Constraints that are inline-capable in the DSL (every `AtomConstraint` variant, every `BondConstraint` variant; the ring-membership variants on `DativeBondConstraint`) round-trip through the DSL's per-entity string form.

2. **Molecule-level `Constraints`** — per-scope `IndexMap<Idx, Vec<C>>` plus a flat `Vec<Constraint>` for molecule-scope and combinator (`And`/`Or`/`Not`) forms. Cross-entity constraints (`SubPattern`, `ChargeSum`, `Connected`, etc.) and non-inline-capable per-entity constraints (`DativeBondConstraint::Donor`, `AromaticSystemConstraint::Contains`, etc.) live here.

There is no invariant that the two stores are disjoint; consumers that need the full set of constraints on an entity read both.

## Canonical participant order

Every per-entity participant array is sorted ascending by `NodeId`:

| Entity | Storage | Guarantee |
|---|---|---|
| covalent bond | `Graph::edge_endpoints(eid)` | `[a, b]` with `a ≤ b` |
| dative bond | `DativeBondAst` + `FixedRelationSet` participants | participants `[a, b]` with `a ≤ b`; direction on `DativeBondAst::direction` |
| aromatic system | `VarRelationSet` participants | sorted ascending; duplicates preserved |
| multicenter bond | `VarRelationSet` participants | sorted ascending; duplicates preserved |
| noncovalent bond | `FixedRelationSet` participants | `[a, b]` with `a ≤ b` |

This is the invariant that makes `MoleculeAst: PartialEq + Hash` canonical — two molecules built from the same bonds in different construction orders compare equal and hash identically. Consumers MAY use `MoleculeAst` as a `HashSet` key, a fingerprint-cache key, or a structural-identity hash without a semantic-equality helper.

### Dative direction

Dative bonds are the only intrinsically directional relation. Direction lives in `DativeBondAst::direction: DativeDirection`, **not** in endpoint order:

- `Forward` — donor is `participants[0]`, acceptor is `participants[1]`.
- `Reverse` — donor is `participants[1]`, acceptor is `participants[0]`.

Participants are always sorted ascending (the `FixedRelationSet` invariant from §"Canonical participant order"); `direction` carries the semantic content that would otherwise be lost in the sort. `DativeBondAst::matches` requires direction equality — two dative bonds with the same participants but opposite direction are distinct identities. `DativeDirection::flip` swaps the two variants; it is a no-op on equality since the participants move with it at a higher level (there is no `flip` on the bond alone).

**Construction responsibility.** `MoleculeAst::new` and `MoleculeBuilder::add_dative_bond` take the authored `(donor, acceptor)` pair and set `direction` on the `DativeBondAst` to match: `Forward` when `donor.0 ≤ acceptor.0`, `Reverse` otherwise. Any `direction` value on the incoming `DativeBondAst` is overwritten. Callers MUST NOT hand-set `direction` before construction and expect it to survive; direction is computed from endpoint order at the chemistry-layer boundary.

**Accessors.** `DativeBondView::{donor, acceptor}` read `participants + direction` to return the semantic atoms in the order the caller authored. Consumers that need the authored direction MUST go through `DativeBondView`, not `FixedRelationSet::participants` directly.

**DSL surface.** The dative-string subgrammar (`spec §7.12`) carries no direction token — direction at the EDN surface is the `:donor` / `:acceptor` keys on `dative-bond-entry`. Raising / lowering translates between endpoint-plus-direction and the two EDN keys; the direction bit is invisible to DSL authors.

Covalent bonds carry no direction (symmetric in the AST; stereo wedges, when added, will need their own bit paralleling this pattern). Noncovalent bonds are symmetric and need no direction.

## AST ↔ DSL

The DSL layer (`crate::dsl`) owns surface representation; the AST layer owns semantics. Conversion uses the `FromAst` / `ToAst` traits with a configuration type per AST kind:

| Direction | Trait | Name | Behavior |
|---|---|---|---|
| DSL → AST | `ToAst<A>::to_ast(&self, cfg: &A::Config)` | **raising** | fill `Undetermined` inherent fields per `cfg` mode (Zero / Natural / Normal / Derived); fill cfg-implied default constraints if missing |
| AST → DSL | `FromAst<A>::from_ast(ast: &A, cfg: &A::Config)` | **lowering** | strip fields back to `Undetermined` where the current literal equals the cfg-implied default; drop default constraint entries |

Config types (`AtomAstConfig`, `BondAstConfig`, `AromaticSystemAstConfig`, `MulticenterBondAstConfig`, `DativeBondAstConfig`, `NoncovalentBondAstConfig`) carry per-field mode enums (`IsotopeMode`, `NumericMode`, `ImplicitHydrogenMode`, `UnpairedElectronsMode`, `MultiplicityMode`, `AromaticValenceMode`, `MulticenterValenceMode`). `zeroed()` and `open()` constructors give extreme policies; `with_overrides` applies partial overrides.

**Defaulting lives only in `FromAst` / `ToAst`.** There is no `coerce` or `release` function. The inverse relation (lowering then raising under the same cfg reproduces the original AST) holds for values ground or partial; under asymmetric cfg it is not expected.

## Three-tier invariants

`MoleculeAst::new` enforces **tier 1** (structural validity); `Solver::resolve` enforces **tier 2** (physics invariants); **tier 3** (chemistry conventions) is never enforced at construction.

| Tier | Examples | Enforced where |
|---|---|---|
| 1 — structural | index validity, relation arity, bond order ≥ 0, element/isotope domain | `MoleculeAst::new` |
| 2 — physics | per-atom electron-count balance, total-charge consistency, total-spin coupling, aromatic-system electron count | `Solver::resolve` (in umol-graph) |
| 3 — chemistry | octet rule, valence tables, Daylight aromaticity, drug-like charge bounds | opt-in validators only; **not** gated on construction |

The general-chemistry scope — organometallics, multicenter bonds, mixed-valence — is the reason for the tier-2 / tier-3 split: RDKit-style sanitization rejects valid compounds at tier 3; we deliberately don't. See `discussion/86-molecule-ast-api-2026-04-16.md` for the full rationale.

## What this crate does not do

- No solver. Resolution (propagation, valence perception, aromaticity) lives in umol-graph.
- No file-format parsers. SMILES / MOL / CTAB parsing lives in umol-graph; this crate only handles the DSL EDN surface and the inline per-entity string subgrammars.
- No coordinate handling. Coordinates live in umol-geometric as annotations on `Molecule`.
- No caching, no precomputation. `MoleculeAst` is pure data; derived caches live on `Molecule` / `Pattern` runtime wrappers outside this crate.
- No tier-3 validators. Opt-in; not part of construction.
