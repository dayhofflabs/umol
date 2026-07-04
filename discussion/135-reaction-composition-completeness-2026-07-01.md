# 135 — Reaction-composition completeness — 2026-07-01

`ReactionAst::compose(A, B, scope)` builds the sequential composites of applying A then B, one
per admissible overlap of A's product `R_A` with B's reactant `L_B`. The target property
(concurrency): `⋃_C∈compose(A,B,Full) C.apply(H) = ⋃_H'∈A.apply(H) B.apply(H')` as a set of canonical
products, for every host `H`.

Doc 134's I3-prop split this into P1 (apply-equivalence), P2 (`RcAnchored` filter), P3 (well-formed),
P4 (determinism), P5 (empty overlap), P6 (correspondence). This doc is about the one still-open half:
**P1 completeness**, `seq ⊆ composed`.

**Set vs multiset comparison.** The target property compares the two product collections **as sets**
of canonical products. Multiset (multiplicity) equality — each product arising the same *number* of
ways on both sides — is the stronger statement worth attempting, but duplicate overlaps or symmetric
automorphisms can make one side produce a given product more than once, so exact multiplicity need
not hold even when the sets coincide. Land set equality first (`compose_complete_overlay`); if
multiplicity diverges, keep set equality as the master property and assert multiplicity separately
(likely narrower — e.g. modulo automorphism).

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

## Part C — core graph-core work: DPO primitives + R1 overlap enumeration (umol-graph-core)

The chemistry-agnostic pushout / pushout-complement / pullback over the adhesive category a molecule
lives in — `Graph` (atoms + bonds) and the relation-set families (`FixedRelationSet`,
`VarRelationSet`, `FixedVarBirelationSet`, …: typed attributed hyperedges, generic over an opaque data
`D`). The attribute/overlay asts are the `D`; the `meet` enters only through a caller-supplied
`combine: FnMut(&D,&D) -> Option<D>` (`None` = ⊥), so graph-core never inspects `D`. Every primitive
returns the object **plus its morphisms** — the bridge umol-ast reads to meet attributes and track
ids. All of Part C is **additive/green** (new methods; nothing existing changes signature).

- **C0 — graph-level primitives. Done (2026-07-04), green.** `graph.rs`: the `remove` split —
  `remove_cascading` (SqPO) + `try_remove` (DPO, `Option`, dangling check); `remove_node`/`remove_edge`
  → `_cascading`. `algorithms/rewriting.rs`: `pushout` / `pushout_complement` (on `try_remove`) /
  `pullback` over `Graph` + `GraphCorrespondence`, each returning object + morphisms. `[dep: —]`

- **C-a — `Correspondence<Id>::to_remapping()`.** `correspondence.rs`: a total-on-left correspondence
  (every left id mated) as a `Remapping` (left id → partner). Feeds a pushout's node morphism into the
  relation-set pushout, which relabels participants via `apply_remapping` (a `Remapping`). Additive.
  `[dep: C0]`

- **C-b — relation-set `pushout` on the overlay families.** `relation.rs`:
  `pushout(&self, right: &Self, combine: impl FnMut(&D,&D) -> Option<D>) ->
  Option<RelationPushout { object: Self, left: Correspondence<RelationId>, right: Correspondence<RelationId> }>`
  on `FixedRelationSet` (noncovalent), `VarRelationSet` (aromatic, multicenter), and
  `FixedVarBirelationSet` (dative, stereo-atom, stereo-bond). **Same-space contract (documented on the
  method): `self` and `right` must already be in the object participant id-space — the caller relabels
  each side with the existing `apply_remapping` and re-indexes its `D` to the new frame first.** This
  matters because `apply_remapping` re-canonicalizes participants (an `Unordered` factor re-sorts) and
  leaves `D` in the old frame — re-aligning per-participant data (aromatic/multicenter electron counts)
  to the new order is ast-specific, so it stays in umol-ast (the established `apply_remapping` +
  re-index pattern). The pushout is then the **same-space merge**: the overlap is implicit (two
  relations coincide iff their participants match — `find_by_participants`); union, `combine` the `D`
  of each coincidence (any `None` ⇒ inadmissible ⇒ `None`), append the rest. So `meet_glue`'s overlay
  half = `apply_remapping ×2 + umol-ast re-index + this pushout` — exactly parallel to the graph
  pushout (graph-core structure, umol-ast `meet`). Additive. `[dep: C-a]`

  Pushout-complement over relations **reuses the existing `apply_compaction`** — the SqPO overlay
  cascade (an overlay whose participant is deleted drops out); no new op. This is sound because
  **incidence is one-directional: an overlay references nodes/edges as participants, and nothing
  references an overlay.** So deletion only ever cascades *up* the hierarchy (node → its edges → the
  overlays on them), never back down. The DPO dangling discipline (`try_remove`) is therefore needed
  only in the node/edge layer, where an edge's survival constrains its endpoints; overlays, depended
  on by nothing, follow freely by `apply_compaction`.

- **C-c — complete the relation-set pushout surface.** `relation.rs`: the same `pushout` on
  `FixedFixedBirelationSet` and `VarVarBirelationSet` — the uniform library surface (no overlay uses
  these families yet, but the primitive is complete across the families). `[dep: C-b]`

- **C-d — relation-set `pullback`.** `relation.rs`: the shared-hyperedge intersection over the
  families, symmetric with the graph `pullback` (C0). Its named consumer is the composite interface
  and any future rule-conflict/dependency analysis; under the operational compose the interface `K`
  otherwise emerges from `A⁻¹`/`B` on the glue. `[dep: C-b]`

- **C-e — `EmbeddingKind` + the monomorphism edge rule in `modular_product` (R1, internal). (done)**
  `algorithms/common_subgraph.rs`: an `EmbeddingKind { Induced, Monomorphism }` enum + an `embedding`
  parameter on the private `modular_product`, branching its adjacency edge rule — the
  `(Some, None) | (None, Some)` arm becomes **allowed** under `Monomorphism` (`E` omits the edge; it is
  context in one graph) and stays **rejected** under `Induced` (edge-iff-edge, today's behaviour). The
  two entries pass `Induced` for now, so behaviour is unchanged. The Bron–Kerbosch / backtracking walks
  are untouched — a new *mode*, not a new algorithm, and **not** an `Algorithm` enum (no algorithm
  choice — it flips one match arm; §R1). This is the `modular_product` path, **not** `McsAlgorithmKind`
  / `pair_feasible` (that serves the separate *maximum*-search; compose uses
  `enumerate_common_subgraphs`). Additive/green. `[dep: —]`

- **C-f — thread `EmbeddingKind` through the two entries (R1, surface). (done)**
  `algorithms/common_subgraph.rs`: add the `EmbeddingKind` argument to `enumerate_common_subgraphs` and
  `maximal_common_subgraphs` (forwarded to `modular_product`). The {maximal, complete} × {induced,
  monomorphism} 2×2 is then **two methods × this one parameter** — the walk `match alg` untouched, no
  new entries (per no-YAGNI the parameter completes the surface). Migrate **all** callers to `Induced`
  (behaviour-preserving), `compose` (compose.rs:315) included — it stays `Induced` here: its overlap
  bond mapping asserts an overlap bond coincides in R_A (`.expect("an induced overlap bond exists in
  R_A")`, compose.rs:634), so it *panics* on a monomorphism overlap and cannot flip until Part D. Add
  the §R1 **F–Cl** graph-core unit test (the monomorphism overlap the induced rule drops). Breaking→green
  (entry-signature change + caller migration land together). `[dep: C-e]`

**Critical path** C0 (done) → C-a → C-b (the DPO primitives) and, independently, C-e → C-f (R1
enumeration) — **all of C0–C-f done**. Part D (`meet_glue` + the Concurrency-Theorem realization)
consumes **C0–C-b and C-f**, and is where `compose` flips to `EmbeddingKind::Monomorphism` — the
overlap bond mapping's meet-interface + delta-rebasing rewrite, the three fixes the ignored
`compose_complete_overlay` names, must land together. C-c and C-d complete the relation-set surface
(additive) and do not gate Part D.

## Part D — reaction-composition completeness: span-based compose (umol-ast)

Commits to the **span approach** (§172's targeted-vs-span, resolved: span — the C0–C-f DPO primitives
were built for glue+diff, and the spike confirmed the wiring). Instead of R2+R3 as per-field delta
surgery, each overlap builds the **composite span**: glue `R_A` and `L_B` over the overlap `E`, then
read off `L_c → R_c`. The glue computes the R2 interface `meet` and the side-diff computes the R3
rebased deltas **uniformly** — no hand-rewriting of `old` fields per variant. Consumes graph-core
C0–C-f. All work is in `umol-ast`; the span primitives already exist —
`ReactionSpanAst::reverse` (→ `ReactionAst`), `ReactionAst::apply_at` / `from_sides`,
`MoleculeCorrespondence::induce` / `reverse`, per-entity `Lattice::meet`, `difference_to` — so the one
new primitive is `meet_pushout`.

*Non-stereo core:*

- **D-a1 — `MoleculePushout` + `meet_pushout` node/edge layer. (additive)** `ast/compose.rs`:
  `MoleculePushout { object: MoleculeAst, left, right: MoleculeCorrespondence }` (fields mirror
  graph-core `Pushout`) + `meet_pushout(left, right, overlap: &MoleculeCorrespondence)
  -> Option<MoleculePushout>` over topology only — graph-core `pushout(left.raw_graph(),
  right.raw_graph(), overlap.atoms())`, each glued atom/bond datum copied from its origin side and, at a
  coincident entity, `Lattice::meet` (`⊥ → None`, the R2-inadmissible overlap). Tests: shared atom folds
  the two `AtomAst`s; conflicting shared atom `⊥ → None`. `[dep: graph-core C-b (done)]`

- **D-a2 — DAMN overlay gluing in `meet_pushout`. (additive)** The four/aromatic overlays via graph-core
  relation `pushout` with `combine = |x, y| x.meet(y)`; non-coinciding overlays kept as context
  (present-absent under monomorphism). Tests: an overlap carrying an overlay (met) and a context overlay
  (kept). `[dep: D-a1]`

- **D-b1 — `compose_overlay` (non-stereo). (additive)** `compose_overlay(span_a, a_inverse, b, overlap)
  -> Option<ReactionAst>`: `meet_pushout(span_a.rhs(), &b.lhs, overlap)?` → `L_c =
  a_inverse.apply_at(&glue.object, &glue.left)?.rhs()`, `R_c = b.apply_at(&glue.object,
  &glue.right)?.rhs()` → `ReactionAst::from_sides(L_c, R_c, corr)`, `corr` recovered from the two
  `ReactionDerivation` comaps (verify the accessor here). `a_inverse = span_a.reverse()?` hoisted once.
  R2 lives in the glue, R3 in `from_sides`' diff. Tests: the §R1 F–Cl pair. `[dep: D-a2]`

- **D-c1 — rewire `compose_all`, flip to `Monomorphism`, un-ignore. (breaking→green)** Enumerate with
  `EmbeddingKind::Monomorphism` (compose.rs:315), map each overlap through `compose_overlay`, collect;
  delete the superseded manual machinery (`created_atom_ids` / `db_atom` / `db_bond` / `lc_atoms` / the
  `ra_*` plumbing). Un-`#[ignore]` `compose_complete_overlay` (non-stereo generator). Milestone:
  **non-stereo completeness** — `sound` / `complete_overlay` / `dangling_free` / `well_formed` /
  `determinism` green. `[dep: D-b1]`

*Stereo (see "Part D stereo — design"):*

- **D-b0 — `ReactionSpanAst::reverse` stereo remappings. (additive)** reaction_span.rs:1468: build the
  two stereo `reversed_remapping`s (created↔removed swap) like the other families; deltas invert via I6a.
  Test: reverse of a stereo-carrying span. `[dep: —]`

- **D-a3 — stereo overlay gluing in `meet_pushout`. (additive)** Canonicalize each stereo overlay's
  ligand order (`transform_frame`, coset carried) → full-participant relation `pushout` over the two
  `FixedVarBirelationSet` families with `combine = StereoAtomAst` / `StereoBondAst::meet`. Test: two
  frames of one stereo center (reordered) glue + meet; contradictory cosets `⊥`. `[dep: D-a1]`

- **D-b2 — `compose_overlay` stereo. (additive)** Pre-apply `transform_frame` of `A⁻¹` / `B`'s stereo
  overlays into the glue's canonical frame before `apply_at` (rule-AST, per overlap). Test: a stereo
  compose pair (e.g. the cis/trans C=C carbon-swap overlap). `[dep: D-b1, D-a3, D-b0]`

- **D-c2 — drop the stereo bail + sample stereo. (breaking→green)** Remove the `has_stereo_*` +
  stereo-delta bail in `compose_all`; extend `overlay_reaction_strategy` to sample stereo overlays;
  `compose_complete_overlay` now covers stereo. Milestone: **stereo completeness** — full suite green.
  `[dep: D-c1, D-b2]`

**Stages** (green after each) — non-stereo phase: **S0** = {D-a1, D-a2}; **S1** = {D-b1}; **S2** =
{D-c1} → non-stereo completeness. Stereo phase: **S3** = {D-b0, D-a3}; **S4** = {D-b2}; **S5** = {D-c2}
→ stereo completeness. **Critical path** D-a1 → D-a2 → D-b1 → D-c1, then D-b0/D-a3 → D-b2 → D-c2; S3
needs only S0, so it can begin any time after S0 (parallel to S1/S2). Additive subitems carry a transient
`#[allow(dead_code)]` until wired.

**Stereo is in scope** — the `has_stereo_*` bail is *not* stale (I6 built single-reaction stereo but
deliberately deferred *compose*-stereo). Folding it into Part D adds `ReactionSpanAst::reverse` stereo
coverage, a frame-canonical stereo path in `meet_pushout`, and dropping the bail + sampling stereo in
the compose generator — see **Part D stereo — design** below for the frame-threading and the
differing-ligand-set decision.

**Opens** — #1 (glue name/shape) **resolved**: `MoleculePushout { object, left, right }` + `meet_pushout`
(placement: `compose.rs`; promote to `ast/glue.rs` only if it grows). #2 (interface op) **resolved**:
`meet` — the combined *value*, not `is_compatible`'s bool — per-entity `Lattice::meet` + the
relation-pushout `combine`; this is the R2 pullback (`R_A` is A's product state, overlap value =
`meet(R_A, L_B)`). #3: the `L_c`↔`R_c` correspondence recovery from each `ReactionDerivation`'s
glue↔result comap — verify in source at D-b. Remaining opens are the **stereo** decisions below.

## Part D stereo — design

Resolves the compose-stereo pieces I6 deferred, for the span approach.

**Where the frame lives.** `StereoAtomAst = { configuration: StereoConfigurationAst, constraints }`
(stereo.rs, `stereo_element!`; it *does* derive `Lattice`). The site atom and the **ligand frame** (the
`Ordered` `StereoLigand` list) are the overlay *participants*
(`FixedVarBirelationSet<NodeId, …, StereoLigand, …, EntitySpan<StereoAtomAst>>`), **not** the data. So
the coset in `configuration` is stated **relative to that participant frame**, and the derived
(site-less) `StereoAtomAst::meet` meets `configuration` + `constraints` field-wise — **correct only when
the two frames coincide** (cosets are frame-relative).

**The frame-threading problem.** Two sites need a common frame:
- *`meet_pushout` coincidence* — the relation pushout coincides overlays by *participant* equality, so
  the same stereo center stated in a different ligand *order* on the two sides has different
  participants: it won't coincide (lands as two stereo entities on one site), and a forced field-wise
  config-meet across differing frames is meaningless.
- *`apply_at` at the glue* — stereo apply is "same-frame, no reconciliation" (I6c / doc 104): it
  sets / relative-ops the host config in the *rule's* frame. Applying `A⁻¹` / `B` (stated in the
  `R_A` / `L_B` frames) at the glue (a third frame) breaks that.

**Approach — one canonical frame per shared stereo site.** Pick a deterministic ligand order per stereo
site in the glue (ascending glue atom id). Then all realignment localizes to `transform_frame`:
1. *Glue build (`meet_pushout`).* Relabel each side's stereo overlays into glue ids;
   `transform_frame(side_frame, canonical_frame)` each `configuration` so the coset follows the reorder
   (physical config fixed); now same-site / same-ligands specs share participants → the relation
   pushout coincides them → `combine = StereoAtomAst::meet` in the shared frame (`⊥ → None`, the same
   inadmissibility as atoms/bonds).
2. *Rule application (`apply_at`).* Pre-transform `A⁻¹` / `B`'s stereo overlays into the glue's
   canonical frame (the AST `transform_frame`, per overlap — the frame permutation comes from the
   overlap correspondence) so `apply_at`'s same-frame lowering holds. This keeps the realignment on the
   *rule AST* and needs **no** `Edit::TransformFrameStereo*` transact variant (that stays deferred).

**No differing-coordination case.** A stereo center is fully coordinated (virtual ligands — implicit-H,
lone pair — are site properties), and every explicit ligand is named in the rule's L/R, so there are
**no context ligands**. In an atom-admissible glue the shared site's explicit ligand *set* is therefore
identical on both sides — an extra / different explicit ligand would give the site a 5th neighbour →
**over-coordinated → atom-level `⊥`**, rejected before stereo is consulted. So the two overlays always
carry the *same* ligand atoms, differing only in *order* (the overlap's induced permutation), which
step 1's `transform_frame` canonicalization removes → **full-participant** coincidence (site **and**
ligands, exactly as `find_by_participants` keys it) → `meet`. The only stereo `⊥` is a genuine geometric
contradiction (the two cosets disagree in the shared frame). No union-frame-lift / admissibility choice.

**`reverse` stereo gap.** `ReactionSpanAst::reverse` (reaction_span.rs:1468) passes empty stereo
remappings. Fix: build the two stereo `reversed_remapping`s like the other families (created↔removed
swap); the span's stereo columns exist (I6c) and stereo deltas already invert (I6a `inverse`:
`Apply(σ) → Apply(σ⁻¹)`, `Swap` / `Mirror` self-inverse).

**New Part D pieces (fold in):** `meet_pushout` stereo path (per-site canonical frame + `transform_frame`
+ `StereoAtomAst`/`StereoBondAst::meet` over the two `FixedVarBirelationSet` families) → **D-a**;
pre-apply rule-stereo `transform_frame` → **D-b**; `ReactionSpanAst::reverse` stereo remappings → a
prereq of **D-b**; drop the `has_stereo_*` + stereo-delta bail and sample stereo in the compose
generator → **D-c**.

**No decision needed** — settled by this pass: canonical-frame threading → full-participant coincidence
→ `meet`; the `reverse` fix; no `Edit::TransformFrame` variant. Stereo folds into D-a/D-b/D-c as listed;
a differing explicit ligand is an atom-level `⊥`, never a stereo one.