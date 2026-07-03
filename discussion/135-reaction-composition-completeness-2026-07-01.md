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

**Scope (extended 2026-07-03).** Beyond P1 completeness, this doc also tracks **I5 — structural
entity refs**, reopened from 134: genuinely not built. See the final section.

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

## Structural entity refs (I5)

Reopened from 134 §3 — genuinely not built. Today every `<entity>-ref` in the reaction / constraint
surface is `int | keyword` (position or id): a bond or overlay with no `:id` can only be named by
position. **Want:** name a non-atom entity by its *constituents* — a bond by its endpoints, an
aromatic / multicenter system by its members, a dative bond by donors + acceptor, a stereo element by
site + ligands (atoms are the base; no structural form).

**Form** — a uniform structural-map variant, the §4-entry form minus `:type`/`:id`:
`<entity>-ref ::= int | keyword | <structural-map>`, where the map is `{:atoms [..]}` (bond,
noncovalent, aromatic, multicenter), `{:donors [..] :acceptor _}` (dative), or
`{:site _ [:ligands [..]]}` (stereo). Map form (not a bare vector) keeps it self-delimiting where refs
nest inside other vectors (anchor pairs, relational `[ref target]`).

**What exists.** The resolution kernel is done: `find_by_participants` (graph-core, S0a) / the
`<collection>.connecting(participants)` matchers, already driving `induce` and
`substructure::verify_overlays`. §4.1 uniqueness (no two same-constituent entries — extended to
noncovalent + multicenter, decided 2026-06-29) makes each structural match ≤1 hit.

**What remains** — the DSL surface + resolver. Extend the ref grammar with the structural-map variant
in one shared production so it reaches every non-atom ref site at once (reaction `:remove`/`:modify`,
entity + relational constraints, `:bond-order-sum :bonds`, anchor pairs, stereo-bond `:site`), and
resolve the structural variant per entity by its constituent payload (`[AtomRef; 2]` /
`Vec<AtomRef>` / donors + acceptor / site + ligands) through the kernel above. Not a `define_ref`
tweak — the structural variant carries a per-entity payload and a per-entity resolution, so the code
shape is the work.

Structural refs used as an *atom-map* input are tautological — a bond/overlay pair, endpoints being
unordered, only restates the atom bijection `induce` already derives — so `resolve` treats such a
pair as a consistency assertion against the induced correspondence (a contradicting one is an error),
never an override. The useful surface is naming an id-less entity by its parts.

### Resolution — the growing entity registry

Refs resolve during `*Input` → AST conversion, and at that point **there is no built `MoleculeAst`**:
`molecule.rs`'s `into_ast` collects entities into `Vec<(participants, ast)>` and calls `from_parts`
*last*, after constraints resolve; reaction deltas resolve against evolving state (lhs + deltas so
far), held as counts + metadata, not a queryable structure. So structural resolution can't call the
AST-level `find_by_participants` on a finished molecule — it resolves against the state built so far.

`EntityCounts` (per-kind running counts, already grown by the delta loop via `allocate_*`) reshapes
into an **`EntityRegistry`**: per kind a running count + name→id map + a **participant lookup**, grown
incrementally during molecule parsing (unifying with the delta loop). This also enables index-range
checks *as you parse* rather than only at the end. Structural resolution = resolve the inner
atom/bond refs → form the participant key → look up (≤1 hit).

Cost splits by kind, so the hot path stays cheap:

| kind | count | structural lookup | cost |
|---|---|---|---|
| **atom** | many | none — no structural form (base) | free, untouched |
| **bond** | many | `(min,max) → BondId` endpoint map, one insert per bond | O(1) insert, O(1) query |
| **overlays** (D/A/M/N/S) | few | `find_by_participants` over the small collection | O(few) |

The only numerous kind that takes a structural ref is the bond, and a bond is named by its endpoints —
an O(1) endpoint map, never a scan. Atoms have no structural form. Overlays are few, so they reuse the
`find_by_participants` kernel directly. Growing + querying is compatible because refs only ever point
**backward** (atoms before bonds before overlays before constraints; deltas at current state), so a
query always sees its target already registered; removal in the delta loop rides the existing
`IdCompaction`. The one honest asymmetry: bonds use a parse-time endpoint map (a bond is a graph edge,
not a relation set) while overlays call `find_by_participants`.

The resolution context unifies onto `&EntityRegistry`: `resolve(&registry)` replaces today's
`resolve(count, id_to_idx)` / `into_ast(count, metadata)` at every ref site, which is what makes all
sites light up from one change. `Structural` is **input-only** — the AST stores the resolved id with
no memory of structural authoring, so `ToEdn`/`from_ast` still render `Index`/`Id` (same lossiness as
writing index `3` for an entity that has an `:id`).

### Precondition — noncovalent uniqueness by endpoints alone

For a noncovalent structural ref to be unambiguous, noncovalent bonds must be disambiguated by their
**endpoints alone** — no two parallel noncovalent bonds of different kinds on the same pair (dropping
the current §4.1 allowance). The tier-1 entity-structure validator's `noncovalent_structure_check`
currently keys the parallel check on `(pair, kind)`; it must key on the unordered pair alone, and
`NoncovalentBondsParallel` drops its `kind` field. This is the doc-134 §3 decision (2026-06-29) and a
hard precondition. (`:electrons` is independent — structural refs read only participant keys, so the
electron-encoding relocation is an orthogonal cleanup, not a blocker.)

### Implementation plan

Modules: **ast** (precondition) → **dsl foundation** (registry, parsers) → **dsl surface** (refs).
Green after every stage; the sole breaking surfaces are S0a (validator) and S3 (resolve signature).

**S0 — precondition (ast)** — independent, land by S3b
- **S0a** `ast/validate/entity.rs`: `noncovalent_structure_check` keys on the unordered atom pair
  alone (drop `kind`); `NoncovalentBondsParallel` drops its `kind` field; update the §4.1 tier-1 note.
  **breaking (red→green)** — deliberate semantic change, migrate its `#[case]`s. `[dep: —]`

**S1 — shared participant parsers (dsl)** — additive
- **S1a** `dsl`: extract the participant-key readers from the entry parsers — `:atoms [..]`
  (bond/noncovalent/aromatic/multicenter), `:donors [..] :acceptor _` (dative), `:site _ :ligands [..]`
  (stereo) — into shared `read_*` fns; entry parsers delegate. **additive (green)**,
  behavior-preserving. `[dep: —]`

**S2 — entity registry (dsl)** — additive + internal restructure
- **S2a** `dsl`: `EntityRegistry` (reshape of `EntityCounts`) — per kind: running count, name→id map,
  participant lookup (bond `(min,max)→BondId`; overlays index their small collections for
  `find_by_participants`); `register_<entity>(id_name?, participants) -> Id`,
  `find_<entity>_by_participants(..) -> Option<Id>`, count/name accessors; counts preserved; mechanical
  `EntityCounts` rename folded in. **additive (green)** — queries unused yet. `[dep: —]`
- **S2b** `dsl`: grow the registry's participant data incrementally in `MoleculeInput::into_ast`
  (register each entity as parsed, so mid-build sites see it). Counts/results unchanged. **green.**
  `[dep: S2a]`
- **S2c** `dsl`: grow/shrink the registry in the reaction delta loop + reaction-span build (register on
  `Add`; unregister + compact on `Remove`, riding `IdCompaction`). **green.** `[dep: S2a]`

**S3 — structural refs (dsl)** — the breaking rewire that lights up every site
- **S3a** `dsl/refs.rs`: add `Structural(payload)` to the 7 non-atom refs (parametrize `define_ref!`
  with the per-entity payload + participant-resolution; `AtomRef` unchanged); `FromEdn` gains the
  `Edn::Map` arm reusing S1a's readers (rejecting `:type`/`:id`); `resolve` becomes
  `resolve(&EntityRegistry)` — `Index`/`Id` via count + name map, `Structural` resolves inner atom/bond
  refs then `find_*_by_participants` (`StereoBondRef` nests a `BondRef`, one level). `ToEdn`/`from_ast`
  unchanged. **breaking (resolve signature + enum).** `[dep: S1a, S2a]`
- **S3b** `dsl`: migrate every resolution site to `resolve(&registry)` and drop the per-loop
  `id_to_idx` maps — entity entries (stereo-bond `:site`), `constraint.rs`, `relational.rs` (18
  variants), `SubPatternAnchorDsl`, `:bond-order-sum :bonds`, `reaction.rs` deltas, `reaction_span.rs`.
  **red→green.** `[dep: S3a, S2b, S2c, S0a]`

**Critical path** S2a → S2b/S2c → S3a → S3b; S0 and S1 are independent foundations. **Deferrable within
S3b**: the stereo-bond *entry* `:site` structural form is the only mid-build site (the sole reason S2b
grows the registry incrementally rather than at end); if fiddly, ship the reference sites first and add
the entry-site form after. **Confirm before S3a**: stereo resolution keys on **site** (unique per the
validator), so `:ligands` is then an optional frame assertion vs required-match.
