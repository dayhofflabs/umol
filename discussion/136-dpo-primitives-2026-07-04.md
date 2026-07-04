# 136 — DPO primitives — 2026-07-04

umol is a **double-pushout (DPO) attributed graph transformation** system. This note fixes the
vocabulary and maps each categorical primitive to its umol operation. Definitions cite Ehrig, Ehrig,
Prange & Taentzer, *Fundamentals of Algebraic Graph Transformation* (Springer, 2006) by number — Ch. 2–3
for the structural core, Ch. 8–9 for the attributed case.

## Setting

A molecule is an **attributed graph**: atoms are nodes, bonds are edges, each carrying attribute
values.

A reaction is a **production** (Def. 3.1 / Def. 9.1): a span

```
L ←l— K —r→ R
```

with `l`, `r` injective. `L` = reactant pattern, `R` = product, `K` = the **gluing graph** (Ehrig's
term) — the part preserved (neither deleted nor created); `L\K` is deleted, `R\K` is created. The
**inverse** production `p⁻¹ = (R ←r— K —l→ L)` (Def. 3.1) is umol's `reverse`. umol stores this as
`(lhs, deltas)` and can materialize the span as `ReactionSpanAst`.

In the attributed case (Def. 9.1) the three graphs share one attribute term-algebra with **variables**
`X`, and matching assigns those variables — a pattern's variable/wildcard matches any concrete value.
umol generalizes this term model to a **lattice** of attribute values (`⊑` = "more specific", `⊥` =
contradiction): sets and ranges, not just first-order terms. This is the symbolic-graph reading of
attributes (see references), and it is what MØD selects with `LabelType::Term` +
`LabelRelation::Specialisation`/`Unification`. Concretely: **umol's lattice order = `Specialisation`,
`Lattice::meet` = `Unification`.** The compatibility test `is_compatible()` (= `meet ≠ ⊥`, overridable)
is the label-relation check that drives matching and overlap enumeration.

A **match** is a morphism `m: L → G` (umol: `Correspondence` + `substructure_matches`).

## The four structural primitives

**1. Pushout — "glue along a shared part" (Def. 2.16).** Given `f: A→B`, `g: A→C`, the pushout
`D = B +_A C` is "the gluing of `B` and `C` via `A`" (Ehrig). In Sets/Graphs it is the disjoint union
identifying the two images of the shared part, `D = B ⊎ C / ≡` (Fact 2.17), built componentwise on
nodes and edges. `D` is **unique up to isomorphism** (Fact 2.20a).
*Example:* `1—2` glued with `2—3` over `{2}` → `1—2—3`.
*umol:* gluing `R_A` and `L_B` over their overlap; attributed → the meet-glue.

**2. Pushout complement — "delete a matched part, keep the context" (§3.2; Def. 9.8, Fact 9.9).** The
object `D` in the *left* DPO square: `G` with `L\K` removed, `K` still embedded. It exists **iff** the
**gluing condition** holds — identification points and dangling points are all gluing points,
`IP ∪ DP ⊆ GP` (Def. 9.8); when it exists it is **unique up to isomorphism** (Fact 9.9).
*umol:* the deletion half of `apply_at`; `ApplyError::Dangling` **is** the dangling condition failing.

**3. Pullback — "intersection of two maps into a common object" (Def. 2.22).** The categorical
**dual** of the pushout: given `f: C→D`, `g: B→D`, the pullback is `{(c,b) : f(c)=g(b)}`, "a
generalized intersection" (Ehrig); for inclusions, `B ∩ C` (Remark 2.24). A pushout with an injective
leg is also a pullback (Remark 2.25; attributed: Fact 8.16).
*umol:* two roles — (i) the **`Lattice::meet`** on attributes ("both hold" = intersection of admissible
values), and (ii) the composite rule's **interface** `K`, built as a pullback in the concurrency
construction (Def. 9.25, below).

**4. Double pushout = rule application (Def. 3.2 / Def. 9.2).** Applying `p = (L←K→R)` at `m: L→G`:

```
L ←l— K —r→ R
m│  (1)  │k  (2) │n
G ←f— D —g→ H
```

with (1), (2) pushouts: (1) deletes (`D` = pushout-complement, `G = L +_K D`), (2) adds
(`H = R +_K D`) — the two-step construction of Fact 9.11, unique up to isomorphism. This is `apply_at`;
reversing `G ⇒ H` gives `H ⇒ G` via `p⁻¹` (Fact 3.3 = `reverse`).

Because these constructions are **canonical (unique up to iso)** (Facts 2.20a, 9.9, 9.11), results can
be *verified without deriving them*: a passing `compose(A,B).apply(H) == B(A(H))` pins them down.

## The attributed twist — why `Lattice`

Two things happen when you glue attributed graphs, on two different axes:

- **Structure** pushes **out** — the graph pushout unions atoms/bonds (Ch. 2–3; the attributed pushout
  of Ch. 8 fixes the data algebra and glues the graph part, Def. 8.10 / Fact 8.12).
- **Attributes** pull **back** — reconciling the two labels on an identified element is a `meet`
  (their tightest common refinement). In Ehrig's term model this is **unification**: the composition
  construction factors the two matches through a common object using the **most-general unifier** on
  the term data (§9.4, pair factorization Def. 9.20, p. 198). umol replaces the term mgu with a
  general **lattice `meet`** (richer: sets, ranges), and `⊥` (no common refinement) makes the glue
  inadmissible.

So:

> **meet-glue = structural pushout (graph-core) + attribute meet/unification (umol-ast `Lattice`).**

That is why every entity AST implements `Lattice`: the meet is the attribute half of the composition,
and it is exactly MØD's label unification / the symbolic-graph constraint conjunction.

## Composition — the Concurrency Theorem (Def. 9.24–9.25, Thm. 9.26)

`compose(A, B)` is the **E-concurrent production** of the Concurrency Theorem. For each overlap of
`R_A` and `L_B`:

- The glue `E` where `R_A` and `L_B` are **jointly epimorphic** (`E` is covered by both) is an
  **E-dependency relation** (Def. 9.24), provided the two pushout-complements `C₁` (over `K_A→R_A→E`)
  and `C₂` (over `K_B→L_B→E`) exist. umol enumerates the underlying shared subgraphs with
  `enumerate_common_subgraphs`.
- The composite `p_A *_E p_B = (L ← K → R)` (Def. 9.25) is then:
  - `L = L_A +_{K_A} C₁` (pushout) — the composite reactant `L_c`,
  - `R = C₂ +_{K_B} R_B` (pushout) — the composite product `R_c`,
  - **`K = C₁ ×_E C₂` (pullback)** — the composite interface. *This pullback is doc-135's R2* ("the
    composite interface must be a meet, not A's `lhs` alone").
- **Thm. 9.26** gives *synthesis* (`G ⇒ H ⇒ G'` ⟹ `G ⇒ G'` via the composite), *analysis* (the
  converse), and — for epimorphic pairs — a **bijective correspondence up to isomorphism**. That
  bijection **is** the target property `compose(A,B).apply(H) == B(A(H))`: soundness *and*
  completeness.

Operationally (the spike, doc 135): glue `E = R_A ∪ L_B`, take `L_c = A⁻¹(E)` and `R_c = B(E)` via
`apply_at`, then `difference_to(L_c, R_c)`. The pullback interface and the meet fall out of the glue.

MØD calls this **rule composition**; its operators map onto umol's `CompositionScope`:

| MØD operator | meaning | umol |
|---|---|---|
| `rcCommon` | all common-subgraph overlaps | `CompositionScope::Full` (`enumerate_common_subgraphs`) |
| `rcParallel` | disjoint (empty overlap) | the "free rule-algebra sum" |
| `rcSuper` / `rcSub` | one side contained in the other's context | anchored / centre-touching subset |
| `rcId` | lift a graph to an identity rule | — |

## Layering

| primitive | over | home |
|---|---|---|
| match, pushout, pushout-complement, pullback | `Graph` + `Correspondence` | **umol-graph-core** (structural, reusable, testable in isolation) |
| meet-glue, rule application (`apply_at`), diff (`difference_to`) | `MoleculeAst` + `Lattice` | **umol-ast** (structural primitive + attribute meet) |

Today the structural core is embedded inside umol-ast's `transact`/`apply`. Extracting primitives 1–4
into graph-core makes the split explicit and yields graph-rewriting primitives reusable beyond
chemistry. It is optional for the compose work (compose can ride the existing umol-ast primitives), but
it is the principled home — MØD keeps exactly this split via `LabelSettings` (structural algorithms
orthogonal to the label relation), as `graph-core` is orthogonal to `Lattice`.

## References

- Ehrig, Ehrig, Prange, Taentzer, *Fundamentals of Algebraic Graph Transformation*, Springer 2006 —
  Defs 2.16 (pushout), 2.22 (pullback), 3.1/9.1 (production), 3.2/9.2 (DPO); Def. 9.8 + Fact 9.9
  (gluing condition, context uniqueness); Fact 9.11 (two-step application); Fact 8.16 (M-pushouts are
  pullbacks); Defs 9.24–9.25 + Thm. 9.26 (E-dependency, E-concurrent production, Concurrency Theorem).
  *The* primitive reference.
- Andersen, Flamm, Merkle, Stadler, *Rule Composition in Graph Transformation Models of Chemical
  Reactions*, MATCH Commun. Math. Comput. Chem. 80(3):661–704, 2018 — the reference for `compose`.
- Andersen, Flamm, Merkle, Stadler, *A Software Package for Chemically Inspired Graph Transformation*,
  ICGT 2016 (MØD) — the sibling implementation; `LabelSettings` / `LabelRelation`
  (`Specialisation`/`Unification`) = umol's lattice order/`meet`.
- Lack & Sobociński, *Adhesive Categories*, 2004 — the setting where pushouts / pushout-complements
  are unique and well-behaved (the guarantee behind testability).
- Orejas, *Symbolic graphs for attributed graph constraints*, 2011; Orejas & Lambers, *Delaying
  constraint solving in symbolic graph transformation*, 2010 — attributes as constraints with delayed
  solving; the closest match to umol's lattice-pattern attribute model (richer than Ehrig's term mgu).
- AGG (Runge et al. 2012) and Henshin (Arendt et al. 2010; Born et al. 2015) — attributed-graph DPO
  tools with critical-pair / conflict-and-dependency analysis; relevant if umol later analyses
  reaction conflicts and dependencies.
