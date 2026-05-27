# 83 — Constraint Unification Architecture

## Context

Step 6 of the unified AST migration (doc 80) ports valence resolution into the MoleculeAst framework. Scoping this step revealed that resolution, validation, and substructure matching share a common structure: constraint satisfaction via unification.

This doc defines the operational semantics. The pattern language design and constraint representation are in doc 79. The two are complementary: doc 79 defines what constraints look like (the `MoleculeConstraint` enum, `DerivedPred`, the `constraints` vec on `MoleculeAst`); this doc defines how constraints are processed (unification, three-valued outcomes, solver interface, resolution profile).

## MoleculeAst as a constraint system

The MoleculeAst is itself a set of constraints, not something constraints are applied to. Each field is either bound (`Lit(v)`) or free (`None`, `Wildcard`, `Var`). The topology (atoms, bonds, dative bonds, aromatic systems, multicenter bonds) is a structural constraint. A ground term has no free variables. A non-ground term is an underdetermined constraint system.

## Solver equations

Solvers contribute additional constraints to the system:

- **Atom typing**: a disjunction over the atom registry. The atom must unify with exactly one registry row (element, charge, bond types → atom type). Strict identity matching. Currently rejects unresolvable inputs (e.g., Fe2+ high/low spin).
- **Valence counts**: parametric arithmetic equations. Electron budget = valence_electrons(element) - charge; H_count = budget - sum(bond_orders). Parameterized by valence table (allowed valence states per element).
- **Aromaticity perception**: implicit rules. Hückel 4n+2 over ring π-electrons, Clar sextet rule, or HMO-based perception. Configurable: choose the perception algorithm.
- **Dative bond perception**: rules-based, specific to coordination compounds.
- **Chirality assignment**: split by chirality type (central, planar, axial).

Each solver has its own free variables and equations. Atom typing and valence counts are alternatives (not composed). Aromaticity, dative bonds, and chirality are independent and composable.

## Three operations, one structure

All three operations are unification over the constraint system. They differ in interpretation of the three possible outcomes:

| | Determined | Contradictory | Underdetermined |
|---|---|---|---|
| **Resolution** | success (ground term) | rejection | reject or return higher term |
| **Validation** | OK | rejection | OK |
| **Matching** | one assignment | no match | enumerate assignments |

### Resolution

Non-ground AST + solver equations → progressively assign free variables. If no free variables remain, resolution succeeds (ground term). Contradiction → rejection. Free variables remaining → underdetermined (reject or return higher-level term).

Formally: `MoleculeAst -(constraints)→ UnificationResult<MoleculeAst>`

Future extension: `MoleculeAst -(constraints)→ Vec<MoleculeAst>` (enumerate allowed structures) or e-graph encoding.

### Validation

AST + solver equations → unification, checking only for contradictions. Fully determined or underdetermined cases are both OK. Only contradictions lead to rejection.

### Matching

Pattern (free variables + equations) × target (usually ground) → unify, fill free variables on the pattern using target values. Contradiction → no match. Fully determined → exactly one match. Underdetermined → multiple matches. Values of free variables returned as assignments.

Formally: `MoleculeAst -(pattern)→ Vec<Assignment>`

The target does not have to be fully determined.

## Result type

```
UnificationResult<T> {
    Determined(T),
    Underdetermined(T),
    Contradictory,
}
```

Where T is MoleculeAst for resolution, Assignment for matching, () for validation.

## Common shape across all three operations

Each operation follows the same structure, regardless of code path:

1. Collect constraints (from MoleculeAst fields + active solvers)
2. Unify (propagate bindings, detect contradictions)
3. Classify outcome (determined / contradictory / underdetermined)
4. Interpret per operation semantics

The three operations do not have to share a code path. Specialized implementations are expected for performance. But this structure should be visible in each implementation.

## Configuration

Solvers are configurable:

- Atom typing vs valence counts: alternatives, user-selected
- Atom registry and valence table: configurable per domain (need not load the full table)
- Aromaticity perception: Hückel, Clar, or HMO — user-selected
- Dative bond perception: optional, for coordination compounds
- Chirality: optional, split by type

A **resolution profile** groups solver configuration into a coherent chemical model. The same profile is used for both resolution and matching, ensuring consistency.

## Solver interface

Each solver is a constraint source with two calling conventions:

- `refine(ast) → Refinement` — for resolution (propagate bindings, narrow AST)
- `validate(query, target, assignment) → bool` — for matching (check consistency)

Both are instances of unification. They share configuration (tables, registries) but have different input signatures. Two methods on one type, not two separate traits.

Dispatch is through the resolution profile, which holds solver instances and iterates them as a pipeline (resolution) or filter chain (matching/validation).

## Two sources of constraints

Unification draws constraints from two sources:

1. **Declarative constraints on MoleculeAst** — the `constraints: Vec<MoleculeConstraint>` field (doc 79). Sub-patterns, negation, derived predicates, matcher flags. These are part of the query/pattern and travel with the AST.

2. **Solver equations from the resolution profile** — atom typing disjunctions, valence arithmetic, aromaticity rules. These are procedural constraint sources external to the AST, configured per chemical model.

Both feed into the same unification process. The declarative constraints are evaluated via the two-pass matcher (doc 79: relational pass + constraint check). The solver equations are evaluated via the `refine`/`validate` interface. The resolution profile owns the solvers; the MoleculeAst owns the declarative constraints.

Derived views (ring membership, distance matrix, biconnected components) are cached on the target molecule and shared by both: `DerivedPred` constraints from the AST look them up, and solvers (e.g., aromaticity perception) compute into them.

## Formal analogues

The constraint system has substantive connections to four formal frameworks. The implementation is domain-specific (fixed schema, specialized dispatch), but the shape of each operation should be traceable to its formal analogue.

### Correspondences

| Concept | Datalog | SMT (DPLL(T)) | Property graph (Cypher) | Relational algebra |
|---|---|---|---|---|
| Base relations (atoms, bonds, ...) | EDB (extensional facts) | Ground literals | Node/edge store | Base tables |
| Derived predicates | IDB rules | Theory atoms | Computed properties | Materialized views |
| Resolution fixpoint | Semi-naïve evaluation | T-propagate loop | — | — |
| Solver `refine` | Rule firing | Theory propagation | — | — |
| Solver `validate` | Constraint check | Theory consistency check | — | Selection |
| Pattern matching | Query evaluation | — | Pattern-graph matching | Join + select |
| Three-valued result | — | SAT / UNSAT / UNKNOWN | — | — |
| Resolution profile | Stratification order | Theory combination | — | — |
| Constraints vec | Rule body | Formula | WHERE clause | Selection predicate |
| `And`/`Or`/`Not` | Horn clause combinators | Boolean connectives | Boolean operators | Set ops |
| Relation sorts | Relation schemas | Sorts | Node labels + edge types | Typed tables |

### What each framework contributes

**SMT / DPLL(T)** — solver architecture. Each solver in the resolution profile is a theory solver (atom typing theory, valence theory, aromaticity theory). They share constraint state and propagate independently. `refine` = T-propagate (narrow free variables, detect conflicts). `validate` = T-check (consistency without narrowing). Three-valued result maps directly: SAT → Determined, UNSAT → Contradictory, UNKNOWN → Underdetermined. The solver loop (propagate → decide → validate → repeat) is the DPLL(T) main loop minus the SAT core.

**Datalog** — derived predicates and the fixpoint. Base relations are EDB. Derived predicates (ring membership, aromatic electron count, valence sum) are IDB rules that materialize from the base. Resolution is semi-naïve evaluation: iterate constraint producers until no new bindings. Stratification is explicit — valence before aromaticity, aromaticity before stereochemistry. Cached views (ring set, distance matrix, biconnected components) are materialized IDB tables.

**Property graph / Cypher** — matcher shape. A MoleculeAst is a typed property graph: atoms are nodes with properties, bonds are edges with properties. Multiple edge types (localized, dative, noncovalent) with different schemas. Pattern matching is graph-pattern query: a query graph with labeled nodes, typed edges, and property predicates, matched against the target. Post-filters (dative direction, aromatic subset) are WHERE clauses on the match result.

**Relational algebra** — constraint evaluation pipeline. Each constraint evaluation is a selection or semi-join. The post-filter chain in the matcher is a pipeline of selections. `RelationRefs` scoping a `DerivedPred` over specific tuples is a projection + selection.

### Where the analogies break

- **No SAT core.** DPLL(T) has a propositional backbone that drives the search. Here the search is VF2 for matching and domain-specific heuristics (Hund's rule, canonical Kekulé) for resolution. D9's rejection of generic CSP (doc 80) is a rejection of the SAT core.
- **Demand-driven, not bottom-up.** Datalog materializes all derivable facts. Resolution is goal-directed: one molecule in, one ground molecule out.
- **Fixed, small schema.** Six relation sorts, known at compile time, each with a specific Rust type. The fixed schema is what makes specialized dispatch work and what makes a generic relational engine unnecessary.
- **Identity semantics.** Atoms are positionally indexed, not set elements. Atom 0 and atom 1 are distinct even if identical in all properties. Standard graph semantics, not relational-algebraic.

### Design implications

These connections constrain the method shapes even though the implementation is specialized:

- **Solver interface** (from SMT): `refine` should return not just the narrowed AST but whether progress was made (for fixpoint detection) and any new constraints discovered during propagation.
- **Resolution loop** (from Datalog): explicitly stratified. The resolution profile defines strata; within each stratum, solvers run to fixpoint. Cross-stratum dependencies flow in one direction.
- **Matcher** (from property graphs): a graph-pattern query engine. The assignment is a homomorphism from query graph to target graph. Constraints are post-conditions on the embedding.
- **Constraint evaluation** (from relational algebra): multi-relation checks (e.g., dative bond endpoints matching the atom assignment) are joins. The shape should be recognizably join-like even if implemented as direct lookups.

The heterogeneity — typed, enumerated relations rather than uniform tuples — justifies specialized dispatch over a generic engine. But the structure of each operation should remain traceable to its formal analogue.

## Operations on the constraint system

The formal analogues sharpen the three core operations and surface four additional operations on the same structure.

### Core operations

**Resolution** — model completion (SMT), forward chaining to fixpoint (Datalog). Non-ground AST + solver equations → progressively assign free variables. The solver loop follows the DPLL(T) shape: propagate (narrow variables via `refine`), decide (apply preference heuristics — Hund's rule, canonical Kekulé), check (validate consistency), repeat. Stratification follows Datalog: valence stratum completes before aromaticity stratum begins, aromaticity before stereochemistry.

Shape: `MoleculeAst × Profile → UnificationResult<MoleculeAst>`

**Validation** — satisfiability checking without model extraction (SMT), integrity constraint checking (Datalog). Same engine as resolution, different question: "is there a contradiction?" rather than "what is the ground term?" No narrowing, no decision — just propagate and check.

Shape: `MoleculeAst × Profile → UnificationResult<()>`

**Matching** — graph-pattern query (Cypher), query evaluation against a materialized database (Datalog). The assignment is a graph homomorphism from query to target, filtered by post-conditions (WHERE clauses in Cypher terms). Multi-relation checks (dative bond endpoints, aromatic system membership) are semi-joins in relational terms.

Shape: `MoleculeAst × MoleculeAst × Profile → Vec<Assignment>`

### Additional operations surfaced by the analogies

**Subsumption** — entailment between patterns. Does every molecule matching pattern A also match pattern B? This is `A ⊨ B` in SMT, query containment in Datalog, subtype ordering in type theory. A relation between two non-ground terms, not between a pattern and a ground term.

Use cases:
- Pattern ordering: "is `[#6]` more general than `[#6;X3]`?" Most-specific-match-wins dispatch in a pattern library.
- Redundancy elimination: if constraint A entails constraint B in a constraints vec, B is redundant.
- Atom typing hierarchy: is one atom type a specialization of another?

Shape: `(MoleculeAst, MoleculeAst) → bool` (partial order)

Status: useful, not urgent. Becomes necessary when pattern libraries or rule priority systems arrive.

**Enumeration** — given an underdetermined system, enumerate all solutions rather than deciding on one. All-SAT in SMT. Resolution returns the first/preferred solution (via decision heuristics); enumeration is the version that does not decide.

Use cases:
- Resonance structures: multiple valid Kekulé assignments for a single aromatic system.
- Tautomers: multiple valid proton/π-bond placements.
- Registry ambiguity: multiple atom types consistent with the constraints (high/low spin Fe²⁺).

Shape: `MoleculeAst × Profile → Iterator<MoleculeAst>`

Status: deferred. The current "decide" step collapses the enumeration to one. Enumeration is an additive extension — keep the resolution machinery, remove the decision step, iterate.

**Projection** — eliminate variables or relations from a constraint system. Variable elimination in SMT, relational π, logical forgetting. Produces a less specific (more general) term from a more specific one.

Use cases:
- "Forget stereochemistry" — widen chirality fields to None.
- "Project to carbon skeleton" — keep only C atoms and C-C bonds.
- The existing view system (MorganTarget, MatchTarget) is implicit projection: extract the relevant subset into a specialized representation.

Shape: `MoleculeAst × RelationMask → MoleculeAst`

Status: implicit in the view system. Making it explicit as an operation on MoleculeAst is straightforward (set fields to None, filter relation vecs) but not required until pattern algebra or molecule simplification becomes a use case.

**Bidirectional unification** — unification proper, as opposed to one-directional matching. Given two non-ground MoleculeAst values, compute the most general common refinement (most general unifier). Both sides have free variables; the result binds variables on both sides.

Use cases:
- Combining constraints from two partial sources (parser-produced partial molecule + user-supplied constraint pattern).
- Checking compatibility of two patterns without a ground target.

Shape: `(MoleculeAst, MoleculeAst) → UnificationResult<MoleculeAst>`

Status: not yet needed. Currently resolution is one-directional (partial molecule + solver equations) and matching is one-directional (pattern against ground target). The general case exists and should be recognized, but implementing it is deferred until a concrete use case demands it.

### Operation summary

| Operation | Formal source | Shape | Status |
|---|---|---|---|
| Resolution | SMT model completion, Datalog fixpoint | `Ast × Profile → UnificationResult<Ast>` | Step 6 |
| Validation | SMT satisfiability, Datalog integrity | `Ast × Profile → UnificationResult<()>` | Step 6 |
| Matching | Cypher pattern query, Datalog query eval | `Ast × Ast × Profile → Vec<Assignment>` | Done (step 4) |
| Subsumption | SMT entailment, Datalog containment | `Ast × Ast → bool` | Outstanding |
| Enumeration | SMT all-SAT | `Ast × Profile → Iterator<Ast>` | Outstanding |
| Projection | Relational π, logical forgetting | `Ast × Mask → Ast` | Outstanding |
| Bidirectional unification | Term unification | `Ast × Ast → UnificationResult<Ast>` | Outstanding |

Note: UnificationResult -> Solution, Profile -> Chemistry

## Absent vs undetermined: removing Option from AST fields

### The problem

`AtomAst` wraps most fields in `Option<>`. A field set to `None` is ambiguous: does it mean "absent from input" (the source did not mention it) or "undetermined" (the source mentioned it but left it free)? Both cases exist:

- SMILES `C` → charge was never mentioned → `charge: None`
- A pattern `[#6]` → charge is explicitly unconstrained → `charge: None`

The `Wildcard` variant on `ValueAst`, `ElementAst`, `IsotopeAst` already represents "explicitly unconstrained." So `None` and `Some(Wildcard)` collapse to the same semantic meaning during resolution and matching. The only place they differ is serialization: a round-tripping format might want to distinguish "field absent in source" from "field present but unconstrained."

This creates a concrete bug surface: `AtomAst::is_ground()` treats `None` as vacuously ground (correct for patterns where absence means "don't care"), but the solver's `refine` needs `None` to mean "needs narrowing" (undetermined). The workaround (`needs_narrowing()` as a separate predicate) patches the symptom without fixing the representation.

### Formal framework analysis

**Datalog.** Fixed schema. Every column in a relation is either bound or free. There is no "absent column" — if a column exists in the schema, it has a value or is a variable. No absent/undetermined distinction.

**SMT.** Variables exist in the formula or they don't. If a variable appears, it is either assigned (by the model) or free. The formula's signature is fixed. No absent/undetermined distinction.

**Property graphs (Cypher).** Cypher DOES distinguish `null` (property exists, value unknown) from property-not-present (key absent from the node). But this distinction serves schema validation and storage optimization, not constraint solving. During pattern matching, `WHERE n.charge IS NULL` and the absence of a `charge` property on a node produce the same match behavior for most practical queries.

**SQL.** `NULL` conflates "unknown" and "inapplicable" — widely recognized as a design mistake (Date, Codd's later work). The three-valued logic it forces (TRUE/FALSE/UNKNOWN) creates well-documented anomalies. SQL's experience argues against conflating distinct semantic states into one representation, but also against having a single null-like value carry two meanings.

### Conclusion

Absent vs undetermined is a schema/serialization concern, not a value-level concern. For all semantic operations (resolution, validation, matching, subsumption), there is no difference between "field not mentioned" and "field explicitly unconstrained." Both mean: this field places no constraint.

### Decision

1. **Remove `Option<>` from `AtomAst` fields.** All fields get a value. The unconstrained state is represented by the type's own variant.

2. **Rename `Wildcard` → `Undetermined`** across `ValueAst`, `ElementAst`, `IsotopeAst`. Rename `AromaticValenceAst::Unspecified` → `Undetermined`. `SpinStateAst::Wildcard` collapses into `Pair { unpaired: Undetermined, multiplicity: Undetermined }`. The new name reflects semantics (a value not yet determined) rather than syntax (a pattern wildcard `*`).

3. **`Undetermined` is the default.** Construction produces `Undetermined` for all fields. Progressive narrowing replaces `Undetermined` with `Lit(n)` during resolution. `is_ground()` returns true iff no field is `Undetermined`.

4. **No `Absent` variant.** The formal framework analysis shows no semantic operation needs the absent/undetermined distinction. If syntactic reshaping ever requires it (recording which fields were explicitly mentioned in the source), a separate `RawAtomAst` container with `Option` fields can serve that purpose. This is orthogonal to resolution, matching, and all other semantic operations.

5. **Roundtrip fidelity** is handled by `AtomAstConfig` modes (already in `ast/config.rs`), not by AST-level `Option`. The config determines how `Undetermined` is interpreted during serialization (as zero, as required, as derived, etc.).

## Scope for step 6

1. Define `UnificationResult<T>` with the three variants
2. Define the resolution profile as the configuration surface
3. Port atom typing engine: wrap as a solver with `refine` + `validate`
4. Port valence counts engine: wrap as a solver with `refine` + `validate`
5. Wire profile into resolution pipeline
6. Wire profile into matcher post-filter chain
7. Existing optimized implementations do not change; new code is the result type, profile, and dispatch glue

## Step 7 status: discovery half complete (2026-04-15)

Step 7 was previously suspended pending the graph representation question (doc 84). With `umol-graph-core` in place and `MoleculeAst` using CSR `Graph` for topology, the blocker is resolved. All four phases are complete for the discovery half.

### Phases

- **Phase 0**: `MoleculeAst` relation indexing — resolved by doc 84. `MoleculeAst` has CSR `Graph`, `neighbors()`, `bond_order_sum()`, `is_in_aromatic_system()`.
- **Phase 1**: `RingEnumerator::enumerate_ast` — uses `Graph::induced_subgraph` for aromatic atom filtering, BCC + cycle enumeration per component, maps back to original indices.
- **Phase 2**: `find_from_ast` on all three models (Hückel, HMO, Clar). `AromaticityModel::aromatic_systems_ast` dispatch. Temporary petgraph `AtomIndex` bridge via numeric conversion; goes away when GraphIR is removed (step 8).
- **Phase 3**: `AromaticityConfig` on `Solver`. Resolve loop: valence → aromaticity → re-valence. `Solver::resolve` returns `Result<Solution<()>, AromaticityError>`. `MoleculeAst::set_aromatic_systems` replaces the `VarRelationSet` wholesale.

