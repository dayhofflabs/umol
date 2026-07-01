# 135 — Reaction-composition completeness — 2026-07-01

`ReactionAst::compose(A, B, scope)` builds the sequential composites of applying A then B, one
per admissible overlap of A's product `R_A` with B's reactant `L_B`. The target property
(concurrency): `⋃_C∈compose(A,B,Full) C.apply(H) = ⋃_H'∈A.apply(H) B.apply(H')` as a set of canonical
products, for every host `H`.

Doc 134's I3-prop split this into P1 (apply-equivalence), P2 (`RcAnchored` filter), P3 (well-formed),
P4 (determinism), P5 (empty overlap), P6 (correspondence). This doc is about the one still-open half:
**P1 completeness**, `seq ⊆ composed`.

## Current state

Everything compose emits is correct; the emitted set is not yet everything it should emit.

- **Sound** — `compose_sound` / `compose_sound_overlay`: every composite product is a sequential
  product (`composed ⊆ seq`). Green.
- **Well-formed** — `compose_well_formed(_overlay)`: every composite applied at its own `lhs`
  reproduces its `right()`. Green.
- **Dangling-free** — `compose_dangling_free` (tier-2 `DpoValidator`). Green.
- **Incomplete** — `compose_complete_overlay` (`seq ⊆ composed`) is `#[ignore]`d. Some sequential
  compositions have no composite.

So: every composed rule is a genuine A-then-B composite, but the *catalogue* of composites is
missing members. The enumeration split (`maximal_common_subgraphs` vs `enumerate_common_subgraphs`)
and compose's rewire to the complete enumeration are done and correct — they are a prerequisite,
not the whole story. Below are the three requirements that remain, each found by a minimal case the
completeness property shrank to, and how they interact.

## R1 — overlaps must be monomorphisms, not induced

**Symptom.** A = `F` → add `Cl` bonded to it (so `R_A` = F–Cl, **with** the bond). B = `[F, Cl]`
with **no** bond, modify Cl's charge. Sequentially: A builds F–Cl, then B matches `[F, Cl]` into it
and modifies Cl. But `compose(A, B, Full)` never produces that composite, so the sequential product
has no witness.

**Cause.** `apply` matches a reactant into a host by **monomorphism** (`substructure_matches`): the
host may carry edges/overlays the pattern does not mention. Compose enumerates overlaps as **common
induced subgraphs** — the modular product marks the two atom-pairs non-adjacent because `R_A` has the
F–Cl bond and `L_B` does not (`(Some, None)` is an induced disagreement). So the full overlap is not
a common *induced* subgraph and is dropped, even though B matches there monomorphically. The overlap
notion compose uses is stricter than the match notion `apply` uses.

**Approach.** Enumerate common subgraphs under the **monomorphism** (subgraph) notion, where an
overlap `E` embeds into both `R_A` and `L_B` but need not be induced in either. This is a different
**modular-product edge rule**: two node-pairs conflict only when *both* graphs have the edge and the
edges are incompatible; present-vs-absent (and absent-vs-absent) never conflict. Concretely the
adjacency becomes:

| R_A edge | L_B edge | induced (today) | subgraph (needed) |
|---|---|---|---|
| present | present, compatible | adjacent | adjacent |
| present | present, incompatible | not adjacent | not adjacent |
| present | absent | **not adjacent** | **adjacent** (E omits it; R_A's edge is context) |
| absent | present | not adjacent | **adjacent** |
| absent | absent | adjacent | adjacent |

Cliques of this subgraph modular product are the (maximal / all) common subgraphs under
monomorphism; `E`'s edges are the common compatible ones. Note this is not a stricter or looser
enumeration than induced — it is a *different set*.

### Does R1 need a new graph-core primitive?

No new clique-enumeration **algorithm** — the Bron–Kerbosch (maximal) and backtracking (complete)
clique walks are unchanged. What changes is the **modular product** that feeds them: a second edge
rule (induced vs subgraph). So R1 is a new *mode* on the common-subgraph enumeration surface, over
the existing `modular_product` / `subgraphs_from_cliques` machinery — not a new algorithm.

The maximal/complete axis and the induced/subgraph axis are orthogonal, a 2×2:

| | maximal | complete (all) |
|---|---|---|
| **induced** | `maximal_common_subgraphs` (BronKerbosch) | `enumerate_common_subgraphs` (Backtracking) |
| **subgraph (monomorphism)** | — | — (what compose needs) |

Open surface question: pass the edge rule as a parameter (an `induced: bool` / a `SubgraphNotion`
enum) to the two enumeration entry points, or add distinct entry points. The `mces`/`mcis` split
already models "one function per task"; the induced/subgraph choice is arguably a *parameter* of one
task rather than a separate task, since only the modular product differs. To settle.

### Is the complete induced enumeration superseded?

Partly, and it is worth being precise. The **machinery** — the backtracking all-cliques walk, the
`modular_product` helper, `subgraphs_from_cliques` — is exactly what the subgraph variant reuses (a
different edge rule feeding the same clique walk), so it is not wasted; it is the scaffolding R1
builds on. What is superseded is the **induced edge rule for compose's purpose**: compose (the only
intended consumer of the enumeration) needs the subgraph notion, so once R1 lands, neither induced
entry point (`maximal_common_subgraphs`, `enumerate_common_subgraphs`) has a consumer in this
codebase — the MCS operations use McGregor separately. Induced common-subgraph enumeration remains a
standard, legitimate graph operation, so the decision is a library-surface one: keep the induced
entries as a complete 2×2 capability, or trim to what has a consumer (the subgraph-complete variant)
and re-add induced if something needs it. Recommend deciding this together with the surface question
above.

## R2 — the composite interface must be a `meet`, not A's `lhs` alone

**Symptom / why R3 alone is unsound.** When R3's rebasing (below) was implemented on its own it
**broke `compose_sound`**: composites started producing products outside `seq`.

**Cause.** For an overlap entity, `lhs_c` currently carries **A's `lhs`** entity only. B's specificity
on that shared entity (its `lhs` requirements) is enforced *only* implicitly, by B's delta `old`
mismatching at apply-time — the composite happened to apply only where A's product exactly equalled
B's `lhs`. That is a fragile accident, but it kept compose sound. The moment R3 rebases B's `old`
onto A's product, that accidental guard is gone, and the composite matches hosts where A applies but
B would not — unsound.

**Approach.** The composite's overlap-entity interface must be the **pullback of B's `lhs` through
A's deltas**: `lhs_c` overlap entity, after A's deltas run, must match `L_B`. Field by field:

- fields A does **not** modify: `lhs_c` value = `meet(A-lhs, B-lhs)` (both patterns must hold, and the
  value is unchanged pre/post A);
- fields A **does** modify: A's *new* value must match `B-lhs`; if it does not, A produces a state B
  cannot match, so **the overlap is inadmissible** and the composite is skipped.

This makes the composite match exactly the hosts where A applies *and* B then applies — restoring
soundness while allowing R3 to broaden application correctly.

## R3 — B's overlap deltas must be rebased onto A's product

**Symptom.** A = atom `charge Und→0`; B = the same atom `Remove`. Sequentially the atom is removed
(empty product); compose drops the full-overlap composite. Compose accumulates A's + B's remapped
deltas and `Deltas::canonicalize()`s, skipping on error.

**Cause.** B's `Remove` carries B's `lhs` old-state (`charge Und`), but in the composite the atom is
at A's *product* state (`charge 0`) when B removes it. `fold_preserved` folds modify-then-remove by
reverting the field changes onto the removed ast, which requires that ast to be the post-modify
state — so the mismatch errors and the composite is lost. `remap_delta` remaps *ids*, not old-state
*values*.

**Approach.** For B deltas on overlap entities, reset the old-state to R_A's value before
accumulation: `Remove` ast, `ModifyField` `old`, `ModifyConstraint` `old` ← R_A's entity (compose
already has it via `r_a.atom(ru)` / the overlap-bond correspondence). Non-overlap B deltas pass
through. This is correct **only together with R2** (which keeps the interface specific enough to stay
sound).

## Interaction and ordering

The three are one change split three ways:

- **R1** decides *which* overlaps exist (enumerate monomorphism overlaps, incl. those where A's
  product carries context L_B lacks).
- **R2** builds the composite *interface* for each overlap (the `meet`/pullback), so the composite
  matches exactly where the sequential pair would.
- **R3** builds the composite *deltas* (A's, then B's rebased onto A's product).

R3 without R2 is unsound; R2 without R3 leaves the modify-then-remove/​modify folds erroring; both
without R1 still miss the context-edge overlaps. So they land together, most naturally: R1 first (a
graph-core edge-rule mode, self-contained and testable in isolation), then R2 + R3 as one compose
change (both touch the per-overlap composite build). `compose_complete_overlay` un-`#[ignore]`d only
when all three are in.

## A unifying alternative to weigh first

R2 + R3 are "compute the composite's overlap interface and its B-side deltas relative to A's
product." That is exactly what building the **composite span** and diffing would do: materialize
`L_A → R_A` and `L_B → R_B`, glue `R_A` and `L_B` over the overlap `E` (a pushout), read off the
composite `L_c → R_c`, then `diff` `L_c` vs `R_c` for the deltas. The gluing computes the interface
`meet` (R2) and the diff computes the rebased deltas (R3) uniformly, using the existing span /
`deltas_from_states` machinery, instead of hand-rewriting `old` fields per variant. This trades the
targeted per-field rebasing for a span construction over the overlap. Worth deciding between the
targeted approach (R2+R3 as delta surgery) and the span approach (build+diff) before implementing —
the span approach may be both cleaner and closer to the DPO concurrency construction.

## Open questions

- Enumeration surface for R1: parameter vs distinct entries; and whether to keep the induced entries
  (2×2) or trim to consumers.
- R2 pullback: is `Lattice::meet` on `AtomAst`/`BondAst` the right interface operation, and how do
  overlay entities (DAMN) participate?
- Targeted (R2+R3 delta surgery) vs span-based (glue + diff) composition.
- Whether the monomorphism enumeration's exponential blowup matters for real reaction pairs (compose
  overlaps are the small localized `R_A ∩ L_B` fragments, so likely not — but the subgraph edge rule
  admits *more* cliques than induced).
