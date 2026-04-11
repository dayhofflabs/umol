# Pattern language design for MoleculeAst

Date: 2026-04-10

Context: while restructuring `umol_graph::ast` so that `AtomAst` / `BondAst` / `MoleculeAst` form a structural superset of ground molecules and pattern queries, the question came up of what "pattern" means at the molecule level. SMARTS is one answer; this doc surveys the design space and the theory it sits in.

## Graph-level pattern modes

Subgraph isomorphism is dominant but not the only graph-level pattern.

Worth representing:

- **Subgraph isomorphism** — the SMARTS model. The default.
- **Induced subgraph isomorphism** — same matching, but extra edges between matched atoms in the target are forbidden. Matcher flag, not a separate AST.
- **Recursive patterns** — atom/bond predicate references a sub-`MoleculeAst` (SMARTS `[$(...)]`). The only mode that forces `MoleculeAst` to be defined recursively.
- **Reaction patterns (SMIRKS)** — pair of `MoleculeAst`s plus an atom mapping. Distinct top-level AST, not a `MoleculeAst` variant.
- **Stereochemistry constraints** — chirality, E/Z. Orthogonal axis; absent from the current AST.

Not patterns — these are algorithms over a ground or pattern AST, no representation needed:

- Maximum common subgraph
- Similarity / fingerprint
- Path / distance / reachability queries
- Symmetry orbits

Already covered by atom/bond predicates, no graph mechanism needed:

- AND / OR / NOT, sets, ring membership, ring size, degree, aromaticity, formal charge ranges. The value-DSL handles boolean composition.

## Recursion: anchored vs. topological

Three distinct things get called "recursion":

1. **Atom-anchored sub-patterns** (SMARTS `[$(...)]`). Pattern hangs off a vertex; sub-pattern is graph-shaped but attachment is asymmetric. Bounded-radius local check. Cheap, finite, no fixpoint.
2. **Named sub-patterns / macro expansion.** Define `aromatic6 = ...`, reference it elsewhere. Purely syntactic, zero new expressivity.
3. **True graph recursion / fixpoint patterns.** Kleene-star paths (Cypher `-[*]->`), Datalog reachability, hyperedge-replacement grammars. Matches arbitrary-length structure: alkane chains, polymers, reaction networks. Genuine new expressive power.

The asymmetry in atom-anchored recursion is forced: any recursive query that "starts here and walks" needs a root. To get *symmetric* topological recursion you either define named patterns and reference them positionally, or move to a fixpoint language (Datalog / MSO / graph grammars) where recursion is on relations, not atoms.

A "recursive topology AST" *is* a coherent goal, but it's a different design point — Datalog-shaped, not SMARTS-shaped. SMARTS recursion is the weak local form because chemists rarely need unbounded patterns; biochemists doing polymers and reaction networks do.

## Ring size, valence, derived predicates

Both are derived vertex predicates: computable from the graph structure given a vertex, but their truth depends on global structure.

- **Valence** = Σ bond orders incident to v. Vertex-local once you have the edge set; can sit on the atom syntactically without losing anything.
- **Ring size** / in-ring = global computation, presented as a vertex predicate.

Putting them on the atom is convention, not necessity. The alternative is a separate `constraints` section in `MoleculeAst` carrying graph-level predicates indexed by atom (or bond, or atom-pair). That generalizes more cleanly to:

- Distance constraints (atoms i,j separated by ≥ k bonds)
- Connectivity constraints between specific named atoms
- Reaction-mapping constraints
- Symmetry / equivalence requirements

The atom-centric form is a degenerate case where every constraint happens to mention exactly one atom. The constraint-list form is uniform across arities and reads better once multi-atom predicates exist.

## Theory underneath

A molecular graph is a relational structure `(V, E, element, charge, ...)`. SMARTS is one specific point in a well-mapped landscape.

| Layer | Power | Tractability on molecules |
|---|---|---|
| First-order conjunctive queries (joins on relations) | Subgraph isomorphism | NP in general, fast in practice — molecule queries are tiny |
| + negation, disjunction | Full FO | Same |
| + recursion / fixpoint (Datalog) | Reachability, transitive closure, polymer patterns | PTIME data complexity |
| + Kleene paths (Cypher / regular path queries) | Path patterns of unbounded length | PTIME |
| Monadic second-order (MSO) | Connectedness, k-colorability, ring counts, planarity | Linear time on bounded-treewidth graphs (Courcelle) — and molecules have very low treewidth |
| Hyperedge-replacement grammars | Productive: generate molecule families | Used in reaction-network tools (e.g. Mod) |

SMARTS is essentially **conjunctive queries with negation, plus a few derived predicates (ring, valence, aromaticity) and atom-anchored sub-query references**. A fragment of well-studied logics, chosen for ergonomics rather than for expressivity ceilings.

Things to know:

- **Subgraph isomorphism vs. graph homomorphism.** SMARTS is injective (isomorphism). Homomorphism allows two query atoms to collapse onto the same target atom — different semantics, occasionally useful for symmetry queries.
- **Courcelle's theorem.** Any MSO property is linear-time checkable on bounded-treewidth graphs. Molecules are essentially always low-treewidth (treewidth bounded by largest fused ring system). In principle a very expressive declarative query language is available without paying for it.
- **Cypher / openCypher / GQL** is the closest production analogue and worth reading as prior art — gets graph queries right syntactically even without chemistry knowledge.
- **Datalog on graphs** is the canonical declarative recursion. "Find all atoms reachable from X through only single bonds" is Datalog, not SMARTS.

## Practical implication for the AST

The honest framing:

- The **chemistry attributes** (element, isotope, charge, aromatic, valence) are domain-specific and need their own treatment. That's the part SMARTS contributes that no general graph language has.
- The **query / pattern structure** (joins, negation, recursion, paths) is a generic graph language problem with decades of theory and three or four production query languages to crib from.

Designing this fresh and controlling the format, I'd separate them: a small chemistry attribute vocabulary plus a general graph-pattern language, instead of fusing them the way SMARTS does. That makes it natural to add Kleene paths, named sub-patterns, multi-atom constraints, and reaction patterns later without re-doing the AST every time. SMARTS is a cautionary example of what happens when you fuse the two and then can't extend either independently.

## The decision that gates everything else

**Conjunctive queries only, or fixpoint?** Everything else flows from that. SMARTS is the former with a few hacks toward the latter; Cypher / Datalog / MSO are honest commitments to the latter.

## Where the discussion moved

The CQ-vs-fixpoint framing above is the wrong binary. The interesting dimensions are schema generality, constraint language, and engineering budget — in that order.

## What doc 60 already commits to

Most of what looked like open theory in the opening of this doc is already closed in doc 60 (molecule-builder-dsl, 2026-03-21):

- **Homoiconicity.** Ground terms, partial structures, queries, and rules share one grammar distinguished by evaluation context. The Level 1–4 term algebra.
- **Multi-relational schema.** A molecule is five named relations (`:atoms`, `:bonds`, `:dative`, `:aromatic`, `:mc`, `:nc`). `MoleculeAst` is an instance of this schema.
- **Uniform feature matching (Q3).** Subgraph containment applies the same way across all five relations: predicates on each tuple, plus closure under references. This is CQ over a relation schema, not subgraph isomorphism on a single labeled graph.
- **Rules as LHS→RHS pairs with named-label mapping and implicit environment preservation.** This is DPO rewriting without the categorical vocabulary. Reaction networks are on solid ground.
- **EDN + Datalog-with-arithmetic** as the rule/query engine. Fixpoint is available at the engine level whenever a recursive rule is defined; it does not require new AST constructs.
- **Partial types and `Family<Sum>`/`Family<Product>`** (O25–O27) for genuine underdetermination. Already category-theoretic native in spirit.

The open part of 79 is narrower than its opening suggests: the grammar is designed, the matching semantics are committed, what is missing is the matcher implementation, path queries, MSO₂ forward-compatibility, and non-injective matching as an optional mode.

## Relation schemas as the organizing principle

The tableIR ancestry is explicit — two object types, atoms and bonds, which are unary and binary relations. Aromatic systems, stereocenters, and multicenter sets are already in the design as relations of higher arity. `MoleculeAst` is already a relational object; it is encoded as five specific struct fields rather than as a schema-parametric container.

The unifying move: make the AST parametric over a closed set of `RelationSchema`s. Each schema enumerates relation symbols, arities, attribute types, and invariants. The molecular-graph schema is one instance; a reaction schema is another (pair of molecular graphs with a correspondence relation); a crystal schema is a third. Schemas change slowly and are few in number — this is not the open-world type system of a generic database.

What the parametric form buys:

- New relations without AST rewrites — rings as a derived relation, reaction centers, distance restraints, symmetry orbits, stereochemistry as a relation over tuples of atoms
- A single matcher implementation that dispatches by arity and predicate type
- A common serialization surface — every schema round-trips through EDN the same way
- Reaction rules as a schema relating two molecule schemas, not a separate AST kind

Closed-world is the important part. The parametricity should be expressed as a trait with associated types or a top-level enum, not a runtime registry. No schema discovery, no plugin system, no dynamic loading. Adding a new schema means adding a Rust type; this matches the stance that schemas are few and slow to change.

## Constraint processing as the unifying frame

Most algorithms that feel distinct are instances of constraint processing.

- **Valence resolution** (`atom_typing::try_build_spec`, `counts_candidates`). Variables: unpaired electrons, lone pairs. Constraints: electron invariant, registry membership, Hund's rule as heuristic. Solution: one complete spec, or a candidate set when underdetermined.
- **Bond perception** (doc 70). Variables: integer bond orders. Constraints: per-atom valence sums. Objective: separable log-likelihood from distance model. Solvers: Lagrangian dual, LP relaxation, or belief propagation for distributions.
- **Aromaticity perception.** Variables: ring membership, π-electron assignment. Constraints: Hückel count, ring geometry, per-atom aromatic valence.
- **Pattern matching.** Variables: assignment of query atoms to target atoms. Constraints: element and attribute predicates, bond predicates, relation containment.
- **Stereochemistry assignment.** Variables: chirality labels. Constraints: CIP priority rules, ring-fusion geometry.

These share the same shape: variables over finite domains, predicates as relations over tuples of variables, a feasibility or optimality question, a solver. They differ in domain structure (integers, enums, atom indices), constraint arity, and in whether they want one, all, optimal, or marginal distributions over solutions.

The sharing is real for **expression** of constraints, less real for **solving**. The unified layer is the constraint language — Datalog over the relation schema, with arithmetic guards on tuple attributes. The engine evaluates joins, filters tuples by guards, and recurses on rules. Solvers for optimization and probabilistic inference plug in as backends that consume the same tuples; Datalog alone handles membership and containment queries.

This matches the doc 60 commitment directly: the rule engine is Datalog-with-arithmetic. Pattern matching, rule application, derived relations, and fixpoint closures all flow through it. Morgan fingerprints, SSSR, and connectivity have fast specialized implementations, but each can be expressed as a Datalog fragment for correctness-checking.

## Graphs versus relations: engineering budget

PostgreSQL exists because relational algebra on millions of rows needed a research program. On molecules with ≤ 200 atoms and ≤ 10 relations, it does not.

The relevant primitives:

- Joins over small tables: nested loops, no query planner
- One per-relation hash index keyed on the join column
- Semi-naive recursion on tiny input
- Arithmetic guards as direct Rust predicates on tuple fields

petgraph already provides the graph-algorithm side — traversals, connected components, biconnected components, shortest paths, cycles. These are not best expressed in relational algebra; they are best expressed as graph algorithms. The right design uses both:

- **Relational layer** is the source of truth. Canonical `MoleculeAst`, pattern queries, rules, Datalog all live here.
- **Graph view** is a computed cache — `petgraph::Graph<AtomAst, BondAst>` built once from the relation schema. Graph algorithms run on this.
- **Derived relations** are views in the database sense. "Is in ring" is a unary relation computed by petgraph's cycle basis and then queryable through the relation engine. Shortest-path results materialize as binary relations on demand.

The lift from graph algorithms to the relational view happens by exposing their **results** as relations, not by rewriting the algorithms in relational algebra. This is the database technique of materialized views over procedural computations, applied at molecule scale rather than database scale.

Engineering budget for the relational engine: small. The hard parts of SQL engines — query planning, index selection, transaction isolation, buffer management, concurrency — are not needed. What **is** needed: a join loop, a predicate evaluator, a recursion driver, and a small library of derived-relation computations (rings, biconnected components, distance matrix) that wrap petgraph.

## Principled core, chemist-facing surface

Non-negotiable. Input formats are MOL, SMILES, SMARTS, and the EDN DSL. Users are chemists. The theoretical scaffolding is not visible in normal use.

Architectural discipline that preserves both:

- The **surface** is the DSL of doc 60 plus the external format parsers. Chemists write what they already write.
- The **core** is the relation schema, the constraint engine, the Datalog evaluator, and the category-theoretic semantics that make the design coherent.
- **Parsers are the only boundary.** Every input format parses to `MoleculeAst` populated with tuples in the appropriate relations. SMARTS parses to a pattern `MoleculeAst`. Reaction SMARTS parses to a rule. Nothing in the core leaks upward unless the chemist explicitly drops into a more expressive layer.
- **Power users** (cheminformatics engineers, Claude in this project) can write direct queries against the relation schema, skipping the chemist surface. This is a library API, not a user-facing one.
- **Theoreticians** can write MSO formulas, fixpoint rules, or categorical constructions. This is even further from the chemist surface and is not required for any of the three-month deliverables.

The chemist surface breaks if the core requires annotations the surface does not (explicit schema tags, explicit quantifier scopes, explicit relation types on every tuple), or if the core is slow enough that chemist-scale workloads become unpleasant, or if the core leaks its vocabulary into error messages. These are engineering constraints on the implementation, not on the semantics.

## Specializations of general algorithms

The discipline: specialized algorithms must be specializations of theoretically-grounded general ones. The general form gives semantics and correctness; the specialization gives speed.

Illustrative pairs:

- **Morgan fingerprints** — general: fixpoint over atom hashes (each round's hash is a function of the previous round's hashes of neighbors), expressible as a Datalog rule set. Specialization: fast Rust loop over precomputed neighbor tables.
- **SSSR / cycle perception** — general: a specific MSO₂ property (smallest set of independent edge-cycles). Specialization: Horton's algorithm or equivalent. The general form is not evaluated at runtime but exists as the definitional ground truth.
- **Substructure search** — general: conjunctive query evaluation over the relation schema. Specialization: VF2 or a variant.
- **Bond perception** — general: maximum-likelihood integer programming under linear constraints on a bipartite incidence matrix. Specialization: the Lagrangian solver of doc 70.
- **Aromaticity perception** — general: a Datalog rule set over rings and electron counts. Specialization: a procedural pass that mutates `MoleculeBuilder`.

Workflow this enforces: for each new algorithm, the general form is written down first as documentation and test oracle. The specialized implementation is validated against the oracle on a small instance set. If the two disagree, the specialization has a bug.

## Delivery constraints and scope

Three-month window. Non-negotiable: SMILES round-tripping, fingerprints, reactions. Everything else is a nice-to-have.

In scope at the shovel depth described above:

1. **Relation schema as a Rust type** parametric over a closed set. The current five-relation AST becomes the molecular-graph instance. Roughly 1–2 weeks given the recent AST refactor.
2. **Matcher over the relation schema** for conjunctive queries with arithmetic guards. Covers SMARTS and DSL Level 2. Roughly 2–3 weeks; backtracking search with petgraph providing neighbor iteration.
3. **Morgan fingerprints** as a Rust loop, documented via Datalog fixpoint form. A few days of actual fingerprint code.
4. **Reaction rule application** as DPO on the relation schema: match LHS, substitute RHS, preserve environment. Roughly 1–2 weeks once the matcher lands.
5. **SMILES round-tripping.** In progress.

Deferred and documented as extension points, not foreclosed by the above:

- MSO₁ or MSO₂ query evaluator
- Fixpoint engine beyond what Morgan needs
- Path / Kleene queries
- Non-injective matching (graph homomorphism mode)
- Probabilistic modes (BP over bond perception, `Family` types)
- Tree decomposition infrastructure
- Distance constraints as first-class relations
- Reaction networks as a top-level schema

## The three decisions that gate the rest

1. **Is the relation schema parametric or hard-coded?** Parametric over a closed set of Rust types. The molecular-graph schema is the default; reaction, crystal, etc. are added as they land. No open-world registry. Cost: ~1 week refactor. Benefit: forward compatibility with essentially no semantic risk.
2. **Is the matcher built against the relation schema or against petgraph?** Against the relation schema, with petgraph used internally for neighbor iteration and graph-algorithm primitives. The relational view is the semantics; petgraph is an implementation detail of the matcher and of non-matching queries.
3. **Are graph algorithms implemented natively or routed through the constraint layer?** Natively, with each algorithm documented as a specialization of a general constraint-language form. The general form is test oracle and documentation; the specialization is production code.

These three answers together give: relation schema as the AST's organizing principle, matcher as a relational algorithm with petgraph for graph primitives, graph algorithms as specialized native code with documented theoretical lift. Principled core, chemist surface unchanged, delivery window feasible.

## The encoding for patterns

One homoiconic type. Per-tuple fields already pattern-capable. One new slot for cross-tuple constraints that ground terms leave empty.

**Aside on the DSL surface.** `umol-edn` provides native `FromEdn`/`ToEdn` derives, so the EDN wire format is a direct reflection of the Rust struct layout. Decisions below optimize for internal cleanness of `MoleculeAst`; the DSL shape falls out of the derive automatically. No separate parsing target, no ergonomics-driven contortions on the AST, and no obligation to preserve the doc-60 surface key-for-key — if the AST shape changes, the DSL changes with it.

### Core relations (unchanged from current)

The six relation vecs on `MoleculeAst` — atoms, bonds, dative, aromatic systems, multicenter bonds, noncovalent bonds. Each tuple's attribute fields already accept wildcards, bindings, and guards via `AtomAst` / `BondAst` and the value-DSL. Ground term = all fields concrete; pattern = any field may be a wildcard or bound variable. Atom identity inside the AST is still a `usize` index; DSL labels lower to indices at parse time.

### Charge and spin are constraints, not fields

Molecular charge and spin are not independent fields on `MoleculeAst`. They are derived from lower-level features:

- **Charge** = Σ atom.charge + Σ bond.charge + Σ aromatic_system.charge + Σ multicenter.charge. This is already the semantics described in doc 60 ("Charge and spin").
- **Spin** is not a sum — angular momentum coupling requires Clebsch–Gordan composition over the constituent spin states (`SpinState::is_constructible_from` in `umol-data`). Molecular spin is either asserted and validated against the achievable couplings, or computed as the set of consistent couplings.

The current `charge: Option<i64>` and `spin: Option<SpinState>` fields on `MoleculeAst` are *already* constraints in disguise. The `Option` means exactly "constrain the total if present, otherwise unconstrained." That is the definition of a constraint; they should fold into the constraints vec rather than exist as peer fields.

Fold both into `MoleculeConstraint::Derived`:

- `charge: Some(0)` → `Derived { predicate: TotalCharge(=, 0), atoms: all }`
- `spin: Some(singlet)` → `Derived { predicate: TotalMultiplicity(=, 1), atoms: all }`
- `charge: None` / `spin: None` → no entry

One place for everything that constrains a derived quantity, regardless of whether the quantity is a sum (charge), a coupling (spin), a structural property (in-ring), or a geometric one (distance).

### One new field on `MoleculeAst`

```rust
pub struct MoleculeAst {
    // ... existing six relation vecs ...
    pub constraints: Vec<MoleculeConstraint>,   // NEW, subsumes charge/spin
}
```

The current `charge: Option<i64>` and `spin: Option<SpinState>` fields are **removed**. `Vec::new()` is zero-allocation; serialization skips empty — ground terms with no total-quantity assertions are unaffected.

### Constraint enum (starting small)

```rust
pub enum MoleculeConstraint {
    // SMARTS [$(...)] — a sub-pattern hangs off an anchor atom
    SubPattern { anchor: usize, pattern: Box<MoleculeAst> },

    // Assert absence: "atom is not in any aromatic system", "atoms are not bonded"
    NotInRelation { relation: RelationSym, atoms: Vec<usize> },

    // Derived structural or numerical predicates from precomputed views
    // (in-ring, bridge, articulation, ring-size, total-charge, total-spin, distance, …)
    Derived { predicate: DerivedPred, atoms: Vec<usize> },

    // Kleene-style paths (deferred, slot exists)
    Path { from: usize, to: usize, expr: PathExpr },

    // Matcher-mode flag: non-injective, induced, …
    Matcher(MatcherFlag),
}
```

`DerivedPred` subsumes "is in ring", "distance between atoms ≥ k", "total charge = n", "total multiplicity = k", "ring of size s containing atoms", and anything else computable from the target molecule after assignment. All such predicates are backed by cached views on the target.

Sub-patterns *could* alternatively live as an extra field on `AtomAst` since SMARTS puts them inside the atom bracket. Either works; the constraints-vec location keeps `AtomAst` untouched and makes the matcher uniform. Pick one; don't do both.

### What goes where

| Query | Location |
|---|---|
| `[C;H3]`, element + attribute predicates | `AtomAst` fields |
| Bond-order / aromatic-hint predicates | `BondAst` fields |
| Connectivity between named atoms | `bonds` tuple in core relations |
| Aromatic-ring membership with size constraint | `aromatic_systems` tuple |
| 6-membered aromatic ring by specifying 6 member atoms | core relations, same shape as ground |
| SMARTS `[$(C=O)]` on an atom | `SubPattern` constraint |
| "Not in any ring" | `NotInRelation` or `Derived { NotInRing, [a] }` |
| "Total charge −1" | `Derived { TotalCharge(=, −1), all }` |
| "Singlet multiplicity" | `Derived { TotalMultiplicity(=, 1), all }` |
| "Distance ≥ 3 bonds apart" | `Derived { Distance(≥, 3), [a,b] }` later |
| "Any path of single bonds from A to B" | `Path` later |
| "Two query atoms may match one target atom" | `Matcher(NonInjective)` |

The contract: if your query can be expressed in core relations alone, it will be. Constraints are the escape hatch, not the default.

### Matcher: two passes

**Pass 1 — relational (the easy case).** Backtracking search assigning query atoms to target atoms. At each step, check the assigned atom's per-tuple predicates; prune using petgraph's neighbor iteration. When a query bond's endpoints are both assigned, check that a matching target bond exists. Same for higher-arity relations: when all member atoms of a query relation-tuple are assigned, check that some target tuple contains them with matching attributes.

For any pattern with `constraints: []`, this pass is the whole matcher. It reduces exactly to VF2-shaped subgraph isomorphism over the multi-relation schema.

**Pass 2 — constraint check (only if constraints present).** Evaluate each `MoleculeConstraint` against the fixed assignment.

- `SubPattern { anchor, pattern }`: recurse on the matcher with `anchor` pinned.
- `NotInRelation`: relational check against the target.
- `Derived`: look up in a precomputed view (is-in-ring bitset, distance matrix, biconnected components, total-charge sum, …) cached per target molecule.
- `Path`: run the path evaluator (when implemented).
- `Matcher`: change backtracker behavior in Pass 1 (these are read before Pass 1 starts, not after).

Derived views live on the target molecule, not the query. They are computed once per target on first use and cached. Ring membership, distance matrix, biconnected components, total charge, total spin — all computed from the resolved target via petgraph + arithmetic over relation tuples. Queryable as if they were relations but implemented procedurally.

### Why one type, not two

Doc 60's homoiconicity commitment says ground term = degenerate pattern. `MoleculeAst { constraints: [] }` with all concrete per-tuple fields *is* that degenerate case. Splitting into `MoleculeAst` + `MoleculePattern` forces a mode distinction at every API, breaks round-trip symmetry, and duplicates the parser and serializer. Keeping one type costs one almost-always-empty `Vec` on ground terms; it saves type proliferation everywhere else.

### Implementation sequence for the three-month window

1. **Relational matcher over core relations** (no constraints). Enough for SMARTS substructure search that doesn't use `[$(...)]`. ~2–3 weeks.
2. **Per-target derived view infrastructure**: ring membership, total charge, total spin at minimum. Needed for the `a!` / not-in-ring / charge / multiplicity common cases. Petgraph wrap plus arithmetic over relations. Days.
3. **`SubPattern` constraint** — SMARTS `[$(...)]`. Recursive matcher call. ~1 week.
4. **DPO rewriting** on matches: LHS match + RHS substitution + environment preservation. ~1–2 weeks after matcher lands.
5. **Morgan fingerprints** — independent of the matcher, direct Rust loop.

Deferred (slots exist, not implemented): `Path`, `Derived` beyond the starter set, `Matcher` non-injective mode, full MSO evaluator.

### Decision: constraints live on `MoleculeAst`

**Resolved 2026-04-10.** Homoiconicity takes precedence over keeping `MoleculeAst` field-minimal. The `constraints: Vec<MoleculeConstraint>` slot is added to `MoleculeAst` directly. Ground terms carry an empty vec; pattern terms use the slot as needed. No separate `MoleculePattern` wrapper type. One parser, one serializer, one API surface for both ground and pattern terms.

