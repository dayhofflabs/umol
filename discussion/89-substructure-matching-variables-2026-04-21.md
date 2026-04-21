# 89 — Substructure matching: variables and cross-type incidence

Date: 2026-04-21

Prerequisite: docs 80 (unified constraint AST), 83 (unification architecture), 86 (molecule AST API), 87 (constraint taxonomy).

## Context

`MoleculeAst` currently has a single pattern-match constraint variant:

```rust
MoleculeConstraint::SubPattern {
    target_anchor: AtomIdx,
    pattern_anchor: AtomIdx,
    pattern: Box<MoleculeAst>,
}
```

Six expressivity gaps identified:

1. No unanchored substructure match (current `SubPattern` is atom-to-atom anchored).
2. No way to assert a dative bond has a specific donor and/or acceptor atom.
3. Same for multicenter and noncovalent bond participants.
4. No way to assert an aromatic system is overlaid on a specific topology.
5. No way to assert a dative bond is overlaid on a specific localized bond.
6. No way to compose the above with logical operators.

Point 6 is free — `Constraint::{And, Or, Not}` already exists, just needs richer leaves.

## Survey of precedents

Six graph-pattern/reasoning substrates, contrasting variable mechanics and cross-type incidence.

| System | Variables | Incidence | Cross-type slots |
|---|---|---|---|
| Cypher | implicit, typed by slot (`(n)` vs `-[r]-`) | structural (pattern adjacency) | typed relationship paths |
| SPARQL | explicit `?x`, position-agnostic | triple pattern in any slot | one uniform term type |
| Datalog / ASP | implicit uppercase, untyped but position-typed | relational predicate of any arity | arbitrary-arity predicates |
| SMARTS | none — pattern IS binding | atom/bond in linear string | not first-class; aromaticity is a flag |
| Gremlin | `as('x')` + `select('x')` labels | pipeline steps | N/A — imperative |
| DPO (graph rewriting) | pattern nodes/edges ARE vars | homomorphism | typed graphs |

### What each system teaches

**Cypher.** `(a)-[r:TYPE]->(b)` — incidence is structural, no separate "incident" operator. The pattern's own adjacency tells you `a` is the source of `r` and `b` is the target. Variables are implicit on name-mention and typed by syntactic slot. Cross-type relations are just typed paths.

**SPARQL.** Variables are position-agnostic because RDF is uniform: every triple is `(term, term, term)`. `?x a ?y` puts a variable in any slot. Operators: BGP conjunction, `OPTIONAL`, `UNION`, `MINUS`, `FILTER`, `EXISTS`/`NOT EXISTS`.

**Datalog / ASP.** Variables are untyped symbols ranging over the Herbrand base. Typing is imposed by predicate position. Existential quantification is body-local; universal is head-universal. ASP adds choice rules and integrity constraints for NP-hard search. Most expressive for open-ended reasoning; most demanding mindset.

**SMARTS.** No variable layer. Each atom in the pattern string is a concrete pattern entity that unifies with a target atom at match. Atom maps `:1, :2` in SMIRKS are identity labels, not logic variables. Recursive SMARTS `[$(…)]` is a nested existential over atom environments. Aromaticity, rings, dative bonds are property flags on atoms/bonds, not first-class entities — so SMARTS cannot directly express "aromatic system X overlays atoms A, B, C".

**DPO.** Rule is a span L ← K → R; match is a graph homomorphism L → G. Every pattern node/edge is a variable by virtue of being in the pattern graph. Typed graphs and adhesive categories handle cross-type entities.

## Takeaways

1. **Cypher and DPO converge on the same answer**: pattern entities ARE the variables. There is no separate variable language. Incidence is expressed by structural co-occurrence in the pattern.
2. **SPARQL's position-agnosticism requires a uniform term type.** Chemistry graphs are heterogeneous (atoms, bonds, aromatic systems, …) — this approach doesn't transfer directly.
3. **Datalog/ASP offer maximum expressivity** — first-class structural variables — but add a whole substrate (named vars, unification, binding scope, grounding).
4. **SMARTS is the popularity benchmark for simple-things-simple.** Its expressivity ceiling (no first-class aromatic systems, dative bonds, multicenter bonds) is precisely what we want to exceed, but its usability lesson stands: the atom/bond string IS the pattern; nothing extra to learn.

## Decision

Adopt a Cypher-style model for now:

- **Pattern entities are implicit variables.** The pattern `MoleculeAst` carries its own `AtomIdx`, `BondIdx`, `DativeBondIdx`, `AromaticSystemIdx`, `MulticenterBondIdx`, `NoncovalentBondIdx`. At match time each pattern entity unifies with some target entity.
- **SubPattern anchor is a multi-correspondence struct (B3).** Homogeneous pairs for each of the six entity types:
    ```rust
    struct SubPatternAnchor {
        atoms: Vec<(AtomIdx, AtomIdx)>,           // (target, pattern)
        bonds: Vec<(BondIdx, BondIdx)>,
        dative_bonds: Vec<(DativeBondIdx, DativeBondIdx)>,
        aromatic_systems: Vec<(AromaticSystemIdx, AromaticSystemIdx)>,
        multicenter_bonds: Vec<(MulticenterBondIdx, MulticenterBondIdx)>,
        noncovalent_bonds: Vec<(NoncovalentBondIdx, NoncovalentBondIdx)>,
    }
    ```
- **Cross-type incidence emerges from pattern structure.** To assert "target aromatic A is incident on pattern atom X", include a pattern atom X that is a participant of a pattern aromatic system Y in the pattern itself, and anchor `(A ↔ Y)`. The incidence holds by structural correspondence under the homomorphism.
- **Per-scope entity-identity constraints expand**: add `Donor`, `Acceptor`, `Overlay` variants to `DativeBondConstraint`; `Contains`, `Atoms` to `MulticenterBondConstraint`; `Ends`, `Contains` to `NoncovalentBondConstraint`; `Atoms` to `AromaticSystemConstraint`. These cover points 2, 3, 4 (partial), 5 directly without pattern matching machinery.
- **Unanchored substructure**: an empty `SubPatternAnchor` = unanchored match.
- **Combinators (point 6)** continue to use `Constraint::{And, Or, Not}`.

## What this model cannot express

An unpinned target entity whose correspondence should constrain incidence. Example: "there exists a target atom, not fixed in the anchor, that participates in target aromatic system A and corresponds to some atom in the pattern satisfying property P". Under the Cypher-style model, one can express this by including the property-P atom as a pattern atom in the pattern graph, making it structurally a participant of a pattern aromatic system, and anchoring the aromatic systems. The cost: the pattern graph must carry every such atom as a structural node.

## Preserving the Datalog extension path

The EDN-based DSL is homoiconic — constraints are data. Extending to first-class structural variables later is a change to the DSL grammar, not to any host-language type:

- A pattern atom in an EDN pattern form could optionally carry a variable marker (e.g., `?a`), producing a variable binding in the parsed AST.
- At the AST layer, variable-carrying positions would be a union `AtomRef = Concrete(AtomIdx) | Var(String)`, introduced at anchor and constraint-argument slots only — not in `AtomAst` attribute fields.
- The matcher gains a unification phase; existing ground paths remain.

Nothing in the current design precludes this. The Cypher-style B3 anchor is a subset of what a future Datalog-style anchor would express; adding variables is strictly additive.

## Topology vs overlay: matching architecture implication

Aromatic systems and DMN (dative / multicenter / noncovalent) bonds differ in their relation to the localized bond graph:

- **Aromatic systems** are overlaid on topology. Their participants are topology atoms and their internal bonds are topology bonds. Once a topology match is established, an aromatic-system correspondence is a check on the overlay (participant set, electron count) against the already-fixed atom correspondence — not a separate structural match.
- **DMN bonds** are parallel relations. A dative bond's donor–acceptor pair is generally *not* a topology edge. DMN correspondence is an incidence query over the relation's participant set, independent of topology adjacency. Exceptions exist (borazine has dative bonds that also appear in the localized bond graph) and are expressed via explicit `Overlay` assertions.

Matching therefore splits:
1. **Structural match on topology** (atoms + localized bonds) — standard subgraph homomorphism.
2. **Incidence queries** on aromatic and DMN scopes — evaluated over participant sets with quantification (exact set, containment, existential/universal over participants).

## Per-scope incidence predicate additions

Expansions beyond the plain Cypher-style model to support the two-layer matching:

Electron count is an inherent property of aromatic systems and multicenter
bonds (already a `ValueAst` field on each entity AST), so it is not a
constraint variant.

```rust
AromaticSystemConstraint {
    Atoms(Vec<AtomIdx>),                    // exact participant set
    Contains(AtomIdx),                      // atom is a participant
    ContainsAll(Vec<AtomIdx>),              // participants ⊇ set
    AllAtoms(Box<AtomConstraint>),          // every participant satisfies
    AnyAtom(Box<AtomConstraint>),           // some participant satisfies
}

MulticenterBondConstraint {                 // same shape as AromaticSystemConstraint
    Atoms(Vec<AtomIdx>),
    Contains(AtomIdx),
    ContainsAll(Vec<AtomIdx>),
    AllAtoms(Box<AtomConstraint>),
    AnyAtom(Box<AtomConstraint>),
}

DativeBondConstraint {
    RingCount(ValueAst),                    // existing
    RingSize(ValueAst),                     // existing
    Donor(AtomIdx),
    Acceptor(AtomIdx),
    DonorSatisfies(Box<AtomConstraint>),
    AcceptorSatisfies(Box<AtomConstraint>),
    Parallels(BondIdx),                     // point 5: dative shares endpoints with a localized bond
}

NoncovalentBondConstraint {
    Ends([AtomIdx; 2]),
    Contains(AtomIdx),
    EndsSatisfy([Box<AtomConstraint>; 2]),  // order-free
}
```

### Quantifier payload width

`AllAtoms` and `AnyAtom` carry `Box<AtomConstraint>` — a single-scope atom predicate. Wide enough for current needs (element, valence, degree, etc.). If richer quantifier bodies become necessary (combinator trees, molecule predicates), the payload widens to `Box<Constraint>`. That is a non-breaking change: current call sites wrap their `AtomConstraint` in `Constraint::Atom(idx, c)` or a new per-scope `Constraint` leaf. Recorded as a deliberate choice, not an oversight.

### Dative exception: parallel to a localized bond

`DativeBondConstraint::Parallels(BondIdx)` is the point-5 cross-scope assertion. Its right-hand side is a topology bond index; evaluation checks that the dative's `[donor, acceptor]` equals that bond's endpoints (in either order, or in donor-order — decision deferred to matcher). Compositions like "dative bond runs parallel to a localized bond that is aromatic" use `Constraint::And([DativeBond(d, Parallels(b)), Bond(b, Aromatic)])`.

## Next steps

1. Expand per-scope constraint enums per the listing above (covers points 2, 3, 4 quantified, 5).
2. Replace current `SubPattern { target_anchor, pattern_anchor, pattern }` with `SubPattern { anchor: SubPatternAnchor, pattern: Box<MoleculeAst> }`.
3. Unanchored matches use `SubPatternAnchor::default()`.
4. Defer matching implementation to a later phase; this ticket is AST shape only.
