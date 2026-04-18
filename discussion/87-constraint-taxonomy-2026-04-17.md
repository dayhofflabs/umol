# 87 — Constraint taxonomy

Date: 2026-04-17

Prerequisite: docs 80 (unified constraint AST), 83 (unification architecture), 86 (molecule AST API).

## Terminology

- **Feature**: any structural element — atom, bond, aromatic system, multicenter bond, dative bond, noncovalent bond.
- **Topology**: the graph of atoms and bonds only. Dative, noncovalent, aromatic-system, and multicenter relations are overlays on topology.

## Primitive constraint categories

Eleven categories of primitive constraints. The primitive pattern can also appear as a subpattern.

| # | Name | Scope | Meaning |
|---|---|---|---|
| 1 | Atom predicate | single atom | Atom-typing registry match: element, charge, unpaired electrons, spin multiplicity, lone pairs, isotope. Not the valence predicates. |
| 2 | Atom/topology | atom + incident bonds | Σ bond\_electrons + unpaired + 2·lone\_pairs − charge = neutral\_valence(element). |
| 3 | Bond predicate | single bond | Bond order, charge, unpaired electrons, spin multiplicity. |
| 4 | ~~Bond/topology~~ | ~~bond + incident atoms~~ | Dissolved: perfect matching reduces to per-atom valence balance (category 2) with bond-order domains {1,2}. Maximum matching is COP (see below). |
| 5 | Topology | atom/bond subset | Connectedness, substructure topology, ring membership for atoms and bonds. |
| 6 | Atom/aromatic-system | atom + aromatic system | aromatic\_valence(atom) = # electrons contributed to the aromatic system. |
| 7 | Aromatic system | aromatic system | Total # of electrons (Hückel rule or Clar/HMO analog). |
| 8 | Atom/multicenter, atom/dative | atom + multicenter/dative bond | Valence contribution, parallel to #6. |
| 9 | Multi-atom | atom set | Σ atom\_charges = total charge; coupled spins = total spin. |
| 10 | Multi-bond | bond set | Σ bond\_orders = target (Johnson graphs, combinatorial enumeration). |
| 11 | Chirality | atom/bond + environment | Absolute and relative. Deferred. |

## Emitters and consumers

**Emitters** (produce constraints):

| Source | Categories emitted |
|---|---|
| SMILES | 1, 3, 5 (ring closures as bonds), 6 (aromatic hint on atom) |
| MOL / SDF | 1, 3, 5 |
| EDN DSL (molecule) | 1, 3, 5, 6, 7, 8, 9 |
| SMARTS | 1, 3, 5 (ring, degree, connectivity), 6, subpatterns recursively |
| SMIRKS | SMARTS + atom-mapping (ReactionRule level) |
| MOD | SMARTS superset + label constraints |
| Resolver (valence, aromaticity) | 2, 6, 7 |
| Resolver (atom typing) | 1 (narrows) |
| Operations (kekulize, aromatize) | 4 (if modeled), 7 (re-verify) |
| Manual / programmatic | any |

**Consumers** (evaluate constraints):

| Consumer | Operation |
|---|---|
| Solver refine | narrow features to satisfy 2, 6, 7, 9 |
| Solver validate | check all tier-2 invariants |
| Matcher | evaluate 1, 3, 5 inline; 2, 6, 7 as post-filters; subpatterns recursively |
| kekulize\_all | enumerate under 2 and 10 |
| tautomers | enumerate under 9, 10 with proton migration |
| E-graph saturation | hold 10 as a rewrite invariant |

## Design principles from CSP/SMT literature

### Three evaluation layers

Constraints divide into three layers by evaluation mode. These layers coexist — they are complementary, not competing.

**Attribute predicates.** Direct field comparison on a feature's AST. "This atom has element C and charge 0." Evaluated by comparing a pattern field against a target field. No narrowing, no topology dependence. This is the simplest layer.

**Derived predicates.** Quantities computed from topology and materialized views. "This atom is in exactly 2 SSSR rings." Pure functions over the target molecule, computed once and cached. SMARTS atom-local predicates (`D`, `X`, `H`, `R`, `r`, `v`) live here. Evaluated as lookups against cached views.

**Propagators.** Domain-narrowing engines with specialized algorithms. "The valence balance equation at this atom constrains (unpaired, lone\_pairs) to feasible set {(1,0), (0,2)}." A propagator consumes derived predicates internally — `ValenceBalance` calls `σ_bond_sum` — but its job is domain reduction, not observation.

| Layer | Purpose | When | Example |
|---|---|---|---|
| Attribute predicate | field comparison | match, validate | `atom.charge == 0` |
| Derived predicate | computed quantity lookup | match (SMARTS), validate | `ring_count(a)`, `degree(a)` |
| Propagator | domain narrowing | resolve, enumerate | `ValenceBalance`, `HückelCount` |

A SMARTS pattern uses layers 1 and 2. Resolution uses all three. Enumeration (Kekulization-all, tautomers) uses layer 3 in a loop.

### Why dedicated propagators

Dedicated propagators (as opposed to decomposing into generic arithmetic) exist for three reasons:

1. **Expressibility.** `PerfectMatching` and `Connected` are graph algorithms. They cannot be encoded as scalar arithmetic without bit-blasting.
2. **Propagation strength.** A global `sum(bonds) = N` propagator achieves bounds consistency in O(k). Decomposed into k individual bounds, inference is strictly weaker (arc consistency per pair only). The gap is real for `alldifferent`, `sum`, `gcc` in the CP literature.
3. **Algorithm dispatch.** Hückel verification is modular arithmetic + ring topology. A generic arithmetic engine would need to rediscover the ring structure every time.

### Fixed variants, not a term algebra

The constraint vocabulary is closed: ~9 derived predicates, ~6 propagators, and attribute predicates on feature fields. Additions will come from new applications (enumeration, sampling, optimization on molecular structures) but the number of new constraint kinds will be manageable. A term algebra (where constraints are symbolic expressions in a user-extensible signature) provides infinite extensibility at the cost of losing Rust's exhaustive `match` and compile-time sort-checking. Not justified for this vocabulary size.

### Reification

Constraints compose via Boolean combinators (`And`, `Or`, `Not`). These are operators in the constraint vocabulary, not adjacent to it. SMARTS uses all three: `[C,N]` = Or, `[C;R]` = And, `[!#7]` = Not. Reification enables meta-constraints ("at least 3 of these 5 hold") if needed in future.

### Lexical scope for bindings

Every atom carries its own variables (charge, element, etc.). A `?x` capture binds to the field of the atom (or bond) it is attached to. SubPattern inherits the parent's bindings (doc 80 D5). SMIRKS atom-mapping labels span LHS↔RHS but that is rule-level scope, not constraint-level. No dynamic scoping anywhere.

### Theory combination

Each theory (feature identity, scalar arithmetic, topology, rings/systems, element/registry) owns its own solver. Communication at boundaries follows the Nelson-Oppen pattern: shared equalities propagate across theories. Doc 80 D9's dispatch-by-variant is the realization of this.

### Enumeration semantics

The constraint representation must be stable under enumeration. A solution is a refinement (substitution) of the same constraint system, not an in-place edit. Kekulization-all and tautomer enumeration produce a stream of substitutions over a shared immutable constraint network.

### Materialized views are predicates

Derived views (ring set, distance matrix, biconnected components) are computed once per target molecule and cached. `InRing(a)` is a relation symbol whose extension is materialized on demand, semantically identical to a base relation. Datalog's EDB/IDB distinction applies directly.

### Stratification

Valence depends on bonds; aromaticity depends on valence; stereochemistry depends on aromaticity. A DAG of strata defines the solver's execution order. Pattern matching runs below the strata it consumes (a SMARTS aromatic-atom predicate is below the stratum that produced aromaticity).

## Extensional vs. intensional predicates

AtomAst fields divide into two kinds with different evaluation semantics:

- **Extensional (base) predicates.** Intrinsic to the atom, no external context needed: element, charge, spin, lone\_pairs, isotope\_mass, implicit\_hydrogens. A registry atom HAS `h=2` as a property. Evaluation is field comparison.
- **Intensional (derived) predicates.** Depend on the atom's relationships (topology, overlay relations). A registry atom does NOT have `v=2` as a property — it asserts that conforming instances must have bond-order sum = 2 when evaluated against topology. Evaluation requires molecule context.

In Datalog terms: base predicates are EDB (extensional database); derived predicates are IDB (intensional database, materialized via rules).

### Consequence: AtomAst carries only base fields

The five derived fields (`valence`, `aromatic_valence`, `multicenter_valence`, `donated_pairs`, `accepted_pairs`) move off `AtomAst`. Their values are:

| Field | Where it lives instead |
|---|---|
| `valence` (σ bond order sum) | computed from bonds: `bond_order_sum()` |
| `aromatic_valence` | attribute on the aromatic system relation (per-atom contribution) |
| `multicenter_valence` | attribute on the multicenter bond relation (per-atom contribution) |
| `donated_pairs` | computable from dative bond orders (donor side) |
| `accepted_pairs` | computable from dative bond orders (acceptor side) |

`AtomAst` retains six base fields: `element`, `isotope_mass`, `charge`, `implicit_hydrogens`, `lone_pairs`, `spin`.

### AtomPattern and BondPattern

Standalone containers pairing base AST with per-feature derived constraints. Used for registry entries, SMARTS atom/bond specs, and passing feature specifications independently of any molecule.

```rust
struct AtomPattern {
    ast: AtomAst,
    constraints: Vec<AtomConstraint>,
}

enum AtomConstraint {
    // Topology-derived (categories 2, 6, 8)
    ValenceSum(ValueAst),           // #v — σ bond order sum
    AromaticValence(ValueAst),      // #a — electrons contributed to aromatic system
    MulticenterValence(ValueAst),   // #m — electrons contributed to multicenter bond
    DonatedPairs(ValueAst),         // #d — pairs donated via dative bonds
    AcceptedPairs(ValueAst),        // #r — pairs accepted via dative bonds
    // SMARTS-parity (category 5)
    Degree(ValueAst),               // D — heavy-atom neighbor count
    Connectivity(ValueAst),         // X — total connections incl H
    TotalHCount(ValueAst),          // H — implicit + explicit H
    InRing,                         // R — in any SSSR ring
    RingCount(ValueAst),            // R<n> — # SSSR rings containing atom
    RingSize(ValueAst),             // r<n> — smallest containing ring
}

struct BondPattern {
    ast: BondAst,
    constraints: Vec<BondConstraint>,
}

enum BondConstraint {
    RingBond,                       // bond lies in some ring
}
```

`AtomPattern` and `BondPattern` are owning types (Clone, not Copy — `ValueAst` admits non-Copy `Expr` variants). They are not stored inside `MoleculeAst`; they are independent containers for registries, SMARTS atoms, and APIs.

On `Molecule`, the chemistry-facing read API is `AtomView<'a>` (borrowing), which exposes both base fields and derived values computed from the molecule's topology and cached views.

### Composition: feature constraints into molecule constraints

When an `AtomPattern` is placed at index `i` during `MoleculeAst` construction, its base AST goes into `MoleculeAst.atoms[i]` and its constraints are lifted into the molecule-level constraint vec by attaching the index:

```
AtomPattern {
    ast: AtomAst { element: Lit(C), h: Lit(2) },
    constraints: [ValenceSum(Lit(2)), InRing]
}
  → placed at index 3
  → MoleculeAst.atoms[3] = AtomAst { element: Lit(C), h: Lit(2) }
  → MoleculeAst.constraints += [
        AtomDerived(AtomIdx(3), ValenceSum(Lit(2))),
        AtomDerived(AtomIdx(3), InRing),
     ]
```

Same for `BondPattern` → `BondDerived(BondIdx, BondConstraint)`.

### MoleculeConstraint enum

The molecule-level constraint vec uses a flat enum. Atom/bond constraints are lifted from feature patterns; cross-feature propagators and combinators are molecule-level only.

```rust
enum MoleculeConstraint {
    // Lifted from per-feature patterns
    AtomDerived(AtomIdx, AtomConstraint),
    BondDerived(BondIdx, BondConstraint),

    // Cross-feature propagators (categories 2, 7, 9, 10)
    TotalCharge(ValueAst),
    TotalSpin(SpinStateAst),
    AromaticElectronCount(AromaticSystemIdx, ValueAst),
    MulticenterElectronCount(MulticenterBondIdx, ValueAst),
    BondOrderSum(Vec<BondIdx>, ValueAst),
    Connected(Vec<AtomIdx>),

    // Structural (recursive)
    SubPattern { anchor: AtomIdx, pattern: Box<MoleculeAst> },

    // Combinators
    And(Vec<MoleculeConstraint>),
    Or(Vec<MoleculeConstraint>),
    Not(Box<MoleculeConstraint>),
}
```

`ValenceBalance` is not a constraint-vec entry — it is a propagator invoked by the solver. The solver checks the electron invariant at each atom using base fields from `AtomAst` and derived values from the relation structures. The invariant holds for every ground `Molecule` by construction; listing it per atom in the constraint vec would be redundant.

### Constraint semantics on MoleculeAst

The `MoleculeAst::constraints` vec contains authored assertions — facts asserted about the molecule by its source (parser, user, programmatic construction). The solver verifies them but neither drains nor populates. Interpretation depends on context:

| Stage | Constraint vec entry means |
|---|---|
| Non-ground `MoleculeAst` | goal: solver must satisfy this |
| Ground `Molecule` | fact: known-true (solver verified) |
| `Pattern` | query: target must satisfy for a match |

Same representation, three readings. Homoiconicity at the AST level.

## Concrete mapping: categories to layers

### Attribute predicates (on feature AST fields)

Feature-local predicates. In patterns, the field carries `Undetermined` (matches anything), `Lit(v)` (exact match), or `Expr(...)` (expression match). No constraint-vec entry needed.

| Category | Feature | Fields |
|---|---|---|
| 1 | Atom | `element`, `charge`, `spin` (unpaired/multiplicity), `lone_pairs`, `isotope_mass` |
| 3 | Bond | `order`, `charge`, `spin` |

The atom-typing registry (category 1's "matches or doesn't") is a solver configuration that narrows these fields during resolution, not a constraint-vec entry.

### Derived predicates (in constraint vec via AtomDerived / BondDerived)

Computed from the target molecule's graph and views. Evaluated as lookups. Used in SMARTS patterns and validation. Emitted as `AtomConstraint` / `BondConstraint` on feature patterns; lifted into `MoleculeConstraint::AtomDerived` / `BondDerived` when placed in a molecule.

| Predicate | Sort | Computes | SMARTS | Backed by |
|---|---|---|---|---|
| `ValenceSum(ValueAst)` | Atom → Int | σ bond order sum | `v` (partial) | `bond_order_sum` |
| `AromaticValence(ValueAst)` | Atom → Int | electrons to aromatic system | — | aromatic system relation |
| `MulticenterValence(ValueAst)` | Atom → Int | electrons to MC bond | — | MC bond relation |
| `DonatedPairs(ValueAst)` | Atom → Int | pairs donated via dative | — | dative bond relation |
| `AcceptedPairs(ValueAst)` | Atom → Int | pairs accepted via dative | — | dative bond relation |
| `Degree(ValueAst)` | Atom → Int | heavy-atom neighbor count | `D<n>` | `Graph::degree` |
| `Connectivity(ValueAst)` | Atom → Int | total connections incl H | `X<n>` | degree + implicit\_hydrogens |
| `TotalHCount(ValueAst)` | Atom → Int | implicit + explicit H | `H<n>` | field + neighbor scan |
| `InRing` | Atom → Bool | in any SSSR ring | `R` | `RingSet::contains_atom` |
| `RingCount(ValueAst)` | Atom → Int | # SSSR rings containing atom | `R<n>` | `RingSet` |
| `RingSize(ValueAst)` | Atom → Int | smallest containing ring | `r<n>` | `RingSet::atom_smallest_ring_size` |
| `RingBond` | Bond → Bool | bond lies in some ring | `@` (bond) | `RingSet::contains_bond` |
| `Connected(Vec<AtomIdx>)` | {Atom} → Bool | all in same component | — | BFS / union-find |

### Propagators (in constraint vec or solver configuration)

Cross-feature consistency rules with specialized evaluators. Two modes: **goal** (solver narrows fields to satisfy) and **assertion** (discharged fact, verified on output).

| Propagator | Category | Scope | Equation / Algorithm |
|---|---|---|---|
| `ValenceBalance(atom)` | 2 | atom + all incident bond types | Full electron invariant (see below) |
| `AromaticElectronCount(system, ValueAst)` | 7 | aromatic system | Σ aromatic\_valence over system atoms = count |
| `MulticenterElectronCount(mc_bond, ValueAst)` | 7 parallel | multicenter bond | Σ multicenter\_valence over MC bond atoms = count |
| `TotalCharge(atoms, ValueAst)` | 9 | atom set (empty = all) | Σ atom\_charges = target |
| `TotalSpin(atoms, SpinStateAst)` | 9 | atom set (empty = all) | coupled spins = target |
| `BondOrderSum(bonds, ValueAst)` | 10 | bond set | Σ bond\_orders = target |

This list is not closed — new propagators can be added as applications require.

#### ValenceBalance electron invariant

The per-atom electron conservation equation includes all bond types (σ, aromatic, multicenter, dative). From the old `Atom::check_invariants`:

```
// Orbital side: how electrons are allocated
total_e_orbital = unpaired
    + 2·lone_pairs
    + 2·donated_pairs
    + 2·accepted_pairs
    + 2·implicit_hydrogens
    + 2·valence (σ bond order sum to heavy atoms)
    + aromatic_valence
    + aromatic_increment
    + multicenter_valence

// Source side: where electrons come from
total_e_source = valence_electrons(element) − charge
    + implicit_hydrogens
    + valence
    + aromatic_increment
    + 2·accepted_pairs

// Invariant: total_e_orbital = total_e_source
```

Dative bonds do not need their own propagator. `donated_pairs` and `accepted_pairs` are terms in `ValenceBalance`; the dative bond's order is a `BondAst` field. Cross-feature consistency (dative bond order matches what donor/acceptor claim) falls out of per-atom valence balance.

### Structural match (not in constraint vec)

Substructure topology — the bond graph of a pattern embedding in the target — is the MATCH phase (VF2). It operates on `MoleculeAst` topology directly. `SubPattern { anchor, pattern }` nests this recursively.

### Combinators

`And`, `Or`, `Not` wrap any constraint-vec entry (derived or propagator). SMARTS composition: `[C;R2,R3]` → `And(element=C, Or(RingCount(a,2), RingCount(a,3)))`. The attribute part (`element=C`) stays on `AtomAst`; the derived part goes in the constraint vec.

## Optimization (COP) extension path

Maximum matching and similar maximality constraints are COP (constraint optimization), not CSP. They require an objective function (`maximize Σ a_e`) in addition to feasibility constraints (`Σ_{e incident to v} a_e ≤ 1`).

The cleanest extension: keep the constraint vocabulary as-is, add an `Objective { maximize/minimize, expr: ValueAst }` at the solver or `MoleculeAst` level. The solver loop becomes: find a feasible solution, then iterate with a tightened bound (standard CP branch-and-bound). This is additive — the CSP machinery stays unchanged, the optimization wrapper is a separate layer.

Not needed now. When it is, it does not require redesigning the constraint enum.

## Operational nomenclature

Names for the engine layer. Settled in design; code rename deferred.

### Module layout

```
src/api/                     # chemistry-facing (moved from model/)
  molecule.rs                # Molecule, AtomView, BondView
  pattern.rs                 # MoleculePattern, AtomPattern, BondPattern
src/unify/
  chemistry.rs               # Chemistry
  resolve.rs                 # Resolver + ResolutionError + Progress
  validate.rs                # Validator + ValidationError
  enumerate.rs               # future — kekulize_all, tautomers
  valence.rs                 # ValenceTheory
  aromaticity.rs             # AromaticityTheory + submodules
  chirality.rs               # future
  config.rs                  # per-theory configs
src/matcher.rs               # Matcher — consumes &Chemistry
src/rewrite/                 # future — reactions, SMIRKS
```

### Type names

| Role | Name |
|---|---|
| Top-level theory combination | `Chemistry` |
| Valence theory | `ValenceTheory` (enum: `AtomTyping`, `Counts`) |
| Aromaticity theory | `AromaticityTheory` (enum: `Hueckel`, `Hmo`, `Clar`) |
| Chirality theory (future) | `ChiralityTheory` |
| Propagation state | `Progress` (`Advanced`, `Fixpoint`, `Contradictory`) |
| Errors | `ResolutionError`, `ValidationError` |
| Resolution engine | `Resolver` |
| Validation engine | `Validator` |
| Matching engine | `Matcher` |

### Call shape

Chemistry is passive data. Engines consume it. One shape per operation.

```rust
let chem = Chemistry::default();
let mol  = Resolver::new(&chem).resolve(ast)?;
Validator::new(&chem).validate(mol.ast())?;
let hits = Matcher::new(&chem).find(&pattern, &mol);
```

No free-function duplicate. No method on `Chemistry`. Pattern: `Engine::new(&chemistry).verb(...)`.

### Rationale

- **`unify` over `solver` / `csp` / `smt`.** Doc 83 frames every operation as unification. `csp` is explicitly rejected there (no SAT core); `smt` overclaims (no SAT core, no conflict learning); `solver` collides with "theory solver" used for the per-domain components.
- **`Chemistry` over `Profile` / `Model`.** `Profile` reads as vague configuration; `Model` collides with SMT "model" = satisfying assignment (what `Molecule` effectively is). `Chemistry` names what the container actually is — the chemical model (theories + configuration) under which unification runs.
- **Theory variants as enums, not sibling theory types.** Atom typing and counts solve the same problem with different algorithms. Tactics, not distinct theories. Doc 80 D9 dispatch-by-variant. Sibling theory naming is reserved for genuinely distinct theories (e.g., future `CentralChiralityTheory` vs `AxialChiralityTheory`).
- **Engines over methods on `Chemistry`.** Configuration should not act. `Resolver::new(&chem).resolve(ast)` reads resolver-uses-chemistry, not chemistry-resolves.
- **One shape.** Free functions alongside methods would duplicate the surface. Engines-with-methods chosen.

### Provisional: `Chemistry` name

Bold single-concept naming; may collide with future broader uses of the word. Revisit if usage becomes awkward. The rename is local to `unify/` and the engine constructors — low cost to change.

## Resolved

- **Category 4 (bond/topology).** Dissolved — perfect matching reduces to per-atom valence balance (category 2) with bond-order domains {1,2}. Maximum matching is COP, outside CSP scope; extension path via `Objective` at the solver level.
- **Hybridization.** Decomposes to element + `Degree`; no dedicated variant. C: 2 neighbors = sp, 3 = sp2, 4 = sp3.
- **Extensional vs intensional.** Split: `AtomAst` carries base (extensional) fields only. Derived (intensional) predicates go in constraint vec via `AtomConstraint` / `AtomDerived`.
- **Constraint semantics.** Authored assertions: solver verifies but neither drains nor populates. Same representation for goals, facts, and queries.
- **Solver emission.** Closed: no tier-2 witnesses in the constraint vec. `ValenceBalance` is implicit for every ground `Molecule`. Downstream consumers call `molecule.validate()` if they need verification. Additive if a future consumer requires explicit witnesses.
- **Multi-bond sum (category 10) EDN surface.** Closed: two emission paths. EDN DSL: constraint entry under `:constraints` (`{:bond-order-sum {:bonds [...] :equals N}}`), parser emits `BondOrderSum`. Kekulizer: constructs `BondOrderSum` programmatically per aromatic system. DSL key specified when EDN grammar is extended for constraints.
- **`Pattern` naming.** Renamed to `MoleculePattern` for consistency with `AtomPattern` / `BondPattern`.
- **Operational nomenclature.** Module `unify/`; container `Chemistry`; engines `Resolver`, `Validator`, `Matcher`; theories `ValenceTheory`, `AromaticityTheory`. Chemistry-facing types move to `api/`. See section above.

## Open points

- **`Chemistry` naming.** Provisional — may rename if usage becomes awkward.

## Implementation status (2026-04-17)

Structural refactor landed in phases P1–P5 of plan `sunny-seeking-fiddle`.

- [x] `AtomConstraint` and `BondConstraint` enums (`ast/constraint.rs`) — all 11 `AtomConstraint` variants and `BondConstraint::RingBond`.
- [x] `AtomPattern` and `BondPattern` containers (`api/pattern.rs`), with `coerce`/`release` helpers.
- [x] `MoleculeConstraint` reshape — flat enum with `AtomDerived` / `BondDerived`, cross-feature propagators, `SubPattern`, and `And`/`Or`/`Not` combinators.
- [x] `AtomAst` reduced to six base fields (element, isotope\_mass, charge, implicit\_hydrogens, lone\_pairs, spin). Five derived fields removed.
- [x] `AtomAstConfig` modes pruned; packed atom/bond DSL lifts constraints; serialization unlifts bare `AtomDerived` / `BondDerived` back into atom/bond sugar.
- [x] Registry accepts `AtomPattern` (sugar and explicit form). `ValenceTheory::AtomTyping` filter uses lifted constraints.
- [x] Solver writes lifted constraints during narrowing (`AromaticValence`, `ValenceSum`). Aromaticity stratum stores per-atom pi contribution on the constraint vec, not an atom field.
- [x] `ElectronInvariant` propagator (`unify/propagate.rs`): theory-independent per-atom electron-conservation check. Invoked by `Validator` and matcher post-filter. `ValenceTheory::validate` is the theory-specific feasibility check, kept separate.
- [x] Module layout, `Chemistry` container, engines `Resolver` / `Validator` / `Matcher`, theories `ValenceTheory` / `AromaticityTheory`.
- [x] 617 conformance tests + 3645 lib tests pass.

## Outstanding

- [ ] **Evaluators for SMARTS-parity `AtomConstraint` variants.** `Degree`, `Connectivity`, `TotalHCount`, `InRing`, `RingCount`, `RingSize`, and `BondConstraint::RingBond` — variants exist, no consumer. Requires a matcher-side evaluator that reads cached derived views (ring set, biconnected components, distance matrix) on the target.
- [ ] **Evaluators for cross-feature propagators.** `AromaticElectronCount`, `MulticenterElectronCount`, `TotalCharge`, `TotalSpin`, `BondOrderSum`, `Connected` — variants exist in `MoleculeConstraint`, no evaluator. Blocks tier-2 invariant verification in doc 86.
- [ ] **`SubPattern` matcher recursion.** Variant exists with (de)serialization; no matcher descends into a nested `MoleculeAst` anchored at a parent assignment. This is doc 80 step 9.
- [ ] **Combinator evaluation.** `And` / `Or` / `Not` have no dispatch. Required once any of the above evaluators land.
- [ ] **Derived view caches on `Molecule` / `ResolverCell`.** `OnceLock<RingSet>` is present; distance matrix, biconnected components, per-atom constraint index, packed pattern adjacency land when the first consumer arrives — "gated on actually needing them".
- [ ] **EDN DSL keys for cross-feature propagators** (`:total-charge`, `:total-spin`, `:bond-order-sum`, `:connected`, `:aromatic-electron-count`, `:multicenter-electron-count`). Parser currently emits `None` / ignores; programmatic construction works.
- [ ] **`Chemistry` naming.** Provisional.
- [ ] **Multicenter term in `ElectronInvariant`.** Omitted pending a `MoleculeAst` helper for per-atom multicenter-valence contribution. Covered separately by the `MulticenterElectronCount` propagator.
