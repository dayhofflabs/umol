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

### Resolution — the growing entity namespace

Refs resolve during `*Input` → AST conversion, and at that point **there is no built `MoleculeAst`**:
`molecule.rs`'s `into_ast` collects entities into `Vec<(participants, ast)>` and calls `from_parts`
*last*, after constraints resolve; reaction deltas resolve against evolving state (lhs + deltas so
far), held as counts + metadata, not a queryable structure. So structural resolution can't call the
AST-level `find_by_participants` on a finished molecule — it resolves against the state built so far.

`EntityCounts` (per-kind running counts, already grown by the delta loop via `allocate_*`) reshapes
into an **`MoleculeNamespace`**: per kind a running count + name→id map + a **participant lookup**, grown
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

The resolution context unifies onto `&MoleculeNamespace`: `resolve(&namespace)` replaces today's
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

### Ref-emission priority (roundtrip normalization) — to formalize

A use-site ref is not stored as authored; it is re-derived on render from its target (`render_atom_ref`
→ `Ref::from_ast` = the target's `:id` keyword if it has one, else its positional index). So ref
*form* is normalized on roundtrip: index-vs-keyword already collapses to keyword-if-named-else-index,
and mixed positional/keyword usage does not roundtrip to its authored mix. This is a designed
normalization, not a bug — but it is currently implicit and untested.

Formalize it as one rule: **a use site is emitted in the highest-priority form its target supports,
`keyword > positional > structural`** — descending specificity of what the ref denotes (a deliberate
label on *this* entity > *this* entity by slot > its participants, not the entity itself). Positional
is universal, so the emitted form is always keyword-else-positional; **structural, being lowest, is
never emitted — it is input-only** (this is *why* structural refs don't roundtrip, not a carve-out).
Consequence: I5 needs zero render/roundtrip work — structural is a parse-side-only input form.

**S-notation** (independent of the I5 dependency chain; pins existing behavior now, the structural row
lights up with S3):
- **Sn** spec + tests: state the use-site emission priority in `umol-ast/spec/umol-dsl-spec.md` (the
  ref-grammar section), and add a roundtrip test asserting the collapse — an entity referenced by a
  mix of positional + keyword renders to the keyword (extend with a structural ref once S3 lands). No
  `render_structural` path — it is dead by construction.

### Implementation plan

Modules: **ast** (precondition) → **dsl foundation** (namespace, parsers) → **dsl surface** (refs).
Green after every stage; the sole breaking surfaces are S0a (validator) and S3 (resolve signature).

**S0 — precondition (ast)** — independent, land by S3b **Done**
- **S0a** `ast/validate/entity.rs`: `noncovalent_structure_check` keys on the unordered atom pair
  alone (drop `kind`); `NoncovalentBondsParallel` drops its `kind` field; update the §4.1 tier-1 note.
  **breaking (red→green)** — deliberate semantic change, migrate its `#[case]`s. `[dep: —]`

**S1 — shared participant parsers (dsl)** — additive **Done**
- **S1a** `dsl`: extract the participant-key readers from the entry parsers — `:atoms [..]`
  (bond/noncovalent/aromatic/multicenter), `:donors [..] :acceptor _` (dative), `:site _ :ligands [..]`
  (stereo) — into shared `read_*` fns; entry parsers delegate. **additive (green)**,
  behavior-preserving. `[dep: —]`

**S2 — entity namespace (dsl)** — additive + internal restructure
- **S2a** `dsl`: `MoleculeNamespace` — a **new** type in `dsl/namespace.rs` (not an in-place `EntityCounts`
  rename, which would break the molecule end-of-parse struct literal). Per kind: running count, name→id
  map, participant lookup (bond `(min,max)→BondId`; overlays index their small collections), via
  private `NamedRegistry`/`KeyedRegistry`. `register_<entity>(name?, participants) -> Id`,
  `<entity>_by_participants(..)`, count/name accessors. Introduced alongside `EntityCounts`; the latter
  retires in S3 when resolution moves onto the namespace. **additive (green)** — module `#[allow(dead_code)]`
  until wired. `[dep: —]` **Done**
- **S2b** `dsl`: grow the namespace's participant data incrementally in `MoleculeInput::into_ast`
  (register each entity as parsed, so mid-build sites see it). Counts/results unchanged. **green.**
  `[dep: S2a]` **Done**
- **S2c** `dsl`: make the namespace the **source of truth for naming**, `MoleculeMetadata` a derived
  view (metadata ⊂ namespace — see the note below). Registry gains the **atom-alias table** (aliases and
  entity ids share one namespace, enforced by `check_id_disjoint`); add `MoleculeMetadata::from(&MoleculeNamespace)`
  (id→name by inverting `by_name`; aliases read directly). Change `MoleculeInput::into_ast →
  (MoleculeAst, MoleculeNamespace)` (an internal method, no trait constraint); `MoleculeDsl` formation
  derives the metadata via the projection. **breaking (into_ast return + MoleculeDsl formation),
  internal.** `[dep: S2b]` **Done**
- **S2d** `dsl`: reaction delta loop (`ReactionInput::into_ast`) — **phase A: grow `delta_namespace`,
  additively. Done.** `lhs_namespace` (returned by `into_ast`, S2c) stays as the lhs namespace,
  immutable through the deltas; build `delta_namespace = MoleculeNamespace::continuation(&lhs_namespace)` —
  a new ctor that copies the lhs per-kind counts (via `NamedRegistry`/`KeyedRegistry::with_count`) so
  `register_*` hands out **global** ids continuing the lhs id space, with empty name/participant/alias
  maps, so it holds only the delta-created entities. Grow it on each `Add` arm
  (`delta_namespace.register_<entity>(name, participants)`), monotonically — `Remove`/`Modify` never
  shrink it (the DSL delta pass never compacts — that's apply-time). **Purely additive/green** —
  `delta_namespace` is grown in parallel to everything else and read by nothing yet; `EntityCounts`,
  `from_ast`, the reaction's `namespace`, and the incremental `ReactionMetadata` are all left untouched.
  **Handoff to S3d + S3e (explicit, survives compaction):** S2d leaves in place, unchanged, the live
  `counts` (`EntityCounts`, still the id allocator + resolution bound), the reaction's `namespace` (a
  `MoleculeMetadata` clone grown on `Add`, still what resolution uses), and the **incremental
  `ReactionMetadata` building** (still the roundtrip artifact). `delta_namespace` is grown but read by
  nothing. **S3e** eliminates `EntityCounts` (below). **S3d** does one delta-loop rewrite: sole counter →
  `delta_namespace`, resolution → the two-namespace pair (drop `namespace`), `ReactionMetadata` → derived.
  The single-counter rewire is deliberately *not* done in S2d — it means reordering every `Add` arm to
  take its id from `register_*` and is the same delta-loop surgery S3d does for resolution, so folding it
  into S3d rewrites the loop once, not twice. `[dep: S2c]` **Done**
- **S2e** `dsl`: reaction-span build (`SpanInput::into_ast`) is molecule-shaped and **self-contained**
  (ids are molecule-shaped; it embeds no lhs molecule parse), so it builds its **own** namespace —
  **purely additive/green, done** (like S2d): grow it in parallel (register each atom/alias/bond/overlay
  as resolved), leaving the span's incremental `MoleculeMetadata` and its `EntityCounts` untouched. The
  namespace is grown but read by nothing yet. S3 resolves the span's structural refs (stereo-bond
  `:site`, constraints) against it, derives its metadata as its roundtrip projection
  (`MoleculeMetadata::from`), and drops its `EntityCounts` (S3e). `[dep: S2b]` **Done**

  *Note (metadata ⊂ namespace).* `MoleculeMetadata` is the roundtrip-relevant subset of the namespace:
  its eight `id→name` maps are the exact inverse of the namespace's `by_name`, and `atom_aliases` moves
  into the namespace. Everything else the namespace holds (name→id, participant indexes, counts) is
  parse-only and derivable from `ast + metadata`, so it never belongs in the persistent public type.
  Hence the namespace is the source and metadata a boundary projection — not a merged union, and not
  rebuilt from metadata.

**S3 — resolution & render on the `Namespace` / `Metadata` traits (dsl)** — the rewire that gives
molecule, reaction, and sub-pattern a *single* resolution and a *single* render, so the three can never
diverge. This supersedes the earlier "one grown namespace + stashed counts" sketch, which failed on the
reaction: it could not separate lhs from created at render without either a parallel structure or a
loose-counts boundary. The trait removes the tension entirely.

### The model — the parse/render asymmetry (settled 2026-07-03)

The DSL reads like a program: **a name is defined before it is used** — refs resolve in document order
against what already exists (no forward refs). Render is the inverse and is **order-agnostic** — it
substitutes an `id` by its keyword through the `id ↔ keyword` bijection, recording no provenance. Two
directions of the same bijection, so two matched traits:

| method | direction | trait | reads |
|---|---|---|---|
| `ref::resolve` | keyword / index / participants → id | **`Namespace`** | `<kind>_count`, `find_<kind>_by_keyword`, `find_<kind>_by_participants`, `contains_id`, `find_atom_alias` |
| `ref::from_ast` | id → keyword / index | **`Metadata`** | `<kind>_keyword(id)` (id→keyword) |

`from_ast` never emits `Structural`, so `Metadata` needs no participant index — that asymmetry *is* the
"metadata is a subset of the namespace" you noted. Both traits and the concrete `*Namespace` types are
**`pub`** — they are general lookup tools, not crate-private plumbing.

Three arrangements of the same two traits — resolution and render each written once, generic:
- **Molecule** — one `MoleculeNamespace` / `MoleculeMetadata`.
- **Reaction** — `ReactionNamespace { lhs, delta }` / `ReactionMetadata { lhs, created }`. `delta =
  MoleculeNamespace::continuation(&lhs)` (created-only, counts continuing lhs). Each trait method is
  delta-then-lhs (ids are unique across the reaction, so at most one hit); `<kind>_count` is `delta`'s
  (continuation carries the running total). **Provenance is intrinsic** — lhs entities live in `lhs`,
  created in `delta` — so the boundary needs no stashed counts and roundtrip needs no set-difference; the
  loop grows exactly one structure (`delta`).
- **Sub-pattern** — a molecule inside a constraint, so a *pair* of namespaces used per-side:
  `into_ast_pair(host: &impl Namespace, pattern: &impl Namespace)` resolves the target ref against the
  enclosing `host` (itself generic — `MoleculeNamespace` or `ReactionNamespace`) and the pattern ref
  against the pattern's. **Patterns are anonymous** (no `:id`, no `:atom-aliases`), so the pattern's
  namespace is derived on demand as `MoleculeNamespace::from_ast(&pattern_ast)` (counts + participants,
  no keyword map) and render is index-only against an empty `Metadata` — nothing to carry upward, no
  recursive metadata. Anonymity is a *stated rule*: the pattern parses through the molecule parser, then
  its namespace must have empty keyword maps and no aliases, else `ParseError::InvalidValue` ("a
  sub-pattern must not name entities (`:id`) or define `:atom-aliases`"). This turns today's silent drop
  of pattern `:id`s into a loud rejection; no new error variant.

### Subitems

- **S3a — done (2026-07-03).** `dsl/refs.rs`: `Structural(payload)` on the 7 non-atom refs via a single
  `define_ref!` arm with an optional `structural = <payload>, <parse>, <resolve>` tail (`AtomRef`
  unchanged); payloads mirror each entry's participant portion — `[AtomRef; 2]` (bond, noncovalent),
  `Vec<AtomRef>` (aromatic, multicenter), and the named `DativeBondParticipants` /
  `StereoAtomParticipants` / `StereoBondParticipants` (stereo ligands = `StereoLigandRef`, moved into
  `refs.rs` from molecule.rs). `FromEdn` gains the `Edn::Map` arm (reuses S1a `atoms_pair`/`atoms_vec`,
  rejects `:type`/`:id`). `resolve(&MoleculeNamespace)` implemented — `Index` via `<kind>_count`, `Id`
  via `find_<kind>_by_keyword`, `Structural` via a per-kind `resolve_<e>_structural` (resolve inner refs
  → `find_<kind>_by_participants`; `StereoBondRef` nests a `BondRef`). The molecule.rs entity loops
  migrate off `resolve(count, id_to_idx)`; the `id_to_idx` maps stay only for `check_id_disjoint`
  (folded away in S3h). `resolve` becomes trait-generic in S3b. Green. `[dep: S1a, S2a, S2b]`
  Also done alongside: `MoleculeNamespace` rename — `NamedRegistry`→`KeywordRegistry`,
  `KeyedRegistry`→`EntityRegistry` (flat, no wrap), `by_name`→`find_by_keyword`, `by_participants`→
  `find_by_participants`, `names`→`keywords`, `with_count`→`from_count`, `<kind>_by_name`→
  `find_<kind>_by_keyword`, `iter_atom_aliases`→`atom_aliases`; `MoleculeNamespace` moved to file top.

- **S3b — done (2026-07-03).** The two traits, green/transparent. `dsl/namespace.rs`,
  `dsl/molecule.rs`, `dsl/refs.rs`.
  - `pub trait Namespace` (namespace.rs): the 25-method query surface — per kind `<kind>_count`,
    `find_<kind>_by_keyword`, `find_<kind>_by_participants`, plus `contains_id(&str)` (id-uniqueness
    across all eight kinds + alias names) and `find_atom_alias(&str) -> Option<&AtomDsl>`.
    `impl Namespace for MoleculeNamespace` **delegates** to the existing inherent methods (inherent
    methods shadow trait methods, so `self.foo()` in the impl hits the inherent — no recursion); the two
    new members are direct over the maps. `MoleculeNamespace` is now `pub`.
  - `pub trait Metadata` (molecule.rs): `<kind>_id(&self, id) -> Option<&str>` per kind (the render
    surface). `impl Metadata for MoleculeMetadata` delegates.
  - `ref::resolve` → `pub fn resolve<N: Namespace>(self, &N)` (its structural resolvers `<N: Namespace>`
    too); `ref::from_ast` → `pub fn from_ast<M: Metadata>(id, &M)`. `into_ast` (metadata-scan) left
    concrete — deleted in S3f (it is the last caller — the `SubPattern` stopgap — that keeps it alive).
    Transparent: callers pass `MoleculeNamespace`/`MoleculeMetadata`.
  - **Deferred (not S3b):** renaming `MoleculeMetadata`'s `<kind>_id` accessors → `<kind>_keyword` and
    retiring `set_*_id`. Attempted; it cascades far past the trait work — `ReactionMetadata` carries
    parallel `<kind>_id` getters (delegating to `.lhs()`), and reaction/span render + the property tests
    call them pervasively, tangled with AST-view `.<kind>_id()` methods of the *same spelling*. Split out
    as **S3m** (below); the `Metadata` trait keeps the current `<kind>_id` names until then. `[dep: S3a]`

- **S3c — molecule-side resolution + render onto the traits (breaking → green).** `constraint.rs`,
  `relational.rs`, `molecule.rs`. Redo of the reverted S3b-a/b, now generic. Resolution methods
  (`MoleculeConstraintDsl`/`ConstraintDsl`/`ConstraintsDsl`/`RelationalConstraintDsl::into_ast`, the
  `atom_subset`/`bond_subset` helpers) drop `(counts, meta)` for a single `namespace: &impl Namespace`;
  leaf `ref.into_ast(count, meta)` → `ref.resolve(namespace)`. Render methods (`*::from_ast`,
  `*_subset_from_ast`) take `&impl Metadata`. `molecule.rs` resolves constraints against its
  `MoleculeNamespace`, renders against `MoleculeMetadata`; drop the mid-parse `EntityCounts` literal +
  metadata. `ref::into_ast` (metadata-scan resolution) is *not* deleted here — the `SubPattern`
  pattern-side stopgap (`into_ast_pair` resolving pattern refs via `into_ast`) keeps it alive; it goes in
  S3f. `constraint.rs`/`relational.rs` test fixtures build a `MoleculeNamespace`. Reaction/span callers
  break here, restored in S3d/e. `[dep: S3b]`

- **S3d — reaction (breaking → green).** `dsl/reaction.rs`. `ReactionNamespace { lhs: MoleculeNamespace,
  delta: MoleculeNamespace }` `impl Namespace` (delta-then-lhs, count = delta's); `ReactionMetadata { lhs:
  MoleculeMetadata, created: MoleculeMetadata }` `impl Metadata` (created-then-lhs). The delta loop: seed
  `delta = continuation(lhs_namespace)`; per `Add`, resolve participants against the `ReactionNamespace`
  then `delta.register_<kind>(...)` (the register return is the id — retire `counts.allocate_*`); dup-check
  via `ReactionNamespace::contains_id`. **Delete** the resolution `MoleculeMetadata` clone (the
  `metadata.set_*_id` grow) and `EntityCounts`. At the boundary: `ReactionMetadata.lhs =
  MoleculeMetadata::from(&lhs_namespace)`, `.created = MoleculeMetadata::from(&delta)` — both projected
  once, no incremental writes; reaction aliases are `delta`'s, lhs aliases render inside `.lhs`. `[dep:
  S3c]`

- **S3e — reaction-span (breaking → green).** `dsl/reaction_span.rs`. Same shape over its own
  `ReactionNamespace`/`ReactionMetadata` (or `MoleculeNamespace`/`MoleculeMetadata` where a span side is a
  plain molecule — settle when implementing). `[dep: S3d]`

- **S3f — sub-pattern (breaking → green).** `constraint.rs`, `dsl/namespace.rs`.
  `MoleculeNamespace::from_ast(&MoleculeAst)` — walk the AST entities, register each anonymously (counts +
  participants, empty keyword maps). `into_ast_pair(host: &impl Namespace, pattern: &impl Namespace)`
  (pattern = `from_ast(&pattern_ast)`); `from_ast_pair(host: &impl Metadata, …)` renders the pattern side
  index-only (empty `Metadata`). Add the anonymity check at the sub-pattern parse (pattern namespace must
  have empty keyword maps + no aliases → `InvalidValue`). Delete the S3c stopgap **and `ref::into_ast`
  itself (+ its tests)** — the stopgap was its last caller. In the same `define_ref!` edit, **rename the
  render leaf `ref::from_ast` → `ref::denote`** (the macro's `id → ref`, all eight refs — co-located with
  `into_ast`/`resolve`) and its render call sites. This is *not* `MoleculeNamespace::from_ast` (the
  pattern-namespace ctor above, a distinct AST→namespace constructor) nor/te value-DSL `FromAst`
  `from_ast` — those keep their names. `[dep: S3c]`

- **S3g — eliminate `EntityCounts` (cleanup, green).** Every remaining count reads the `Namespace`
  trait's `<kind>_count`; delete the struct + `from_ast` + `allocate_*` from `constraint.rs`. (Most users
  already gone in S3c/d/e/f.) `[dep: S3d, S3e, S3f]`

- **S3h — id-uniqueness on the namespace (cleanup, green).** Molecule build's scattered locals
  (`check_id_disjoint`, `entry_ids`, `atom_id_to_idx`/`bond_id_to_idx`) collapse onto
  `Namespace::contains_id` (or a register-time check returning `Err(DuplicateId)`); remove the
  `id_to_idx` maps (last use). Minor error-ordering change (check at register-time, after participant
  resolution). `[dep: S3c]`

- **S3i — proptest: structural refs resolve (feature `proptest`).** Off a generated molecule / reaction,
  pick each non-atom entity and form a *structural* ref to it (its constituent atom/bond refs) beside the
  positional ref; assert both resolve to the same id, and that a structural ref over the wrong constituent
  set fails. Cross-checks the `resolve_<e>_structural` path against positional resolution across all seven
  kinds (incl. the stereo `(site, ligand-multiset)` key). `[dep: S3c, S3d]`

- **S3j — proptest: `keyword > positional > structural` emission on roundtrip (feature `proptest`).** The
  render priority: a ref to a *named* entity re-emits as its keyword, to an unnamed entity as its index,
  and a *structural* ref is **never** re-emitted as structural (input-only — `from_ast` produces only
  `Index`/`Id`). Property: parse a DSL form with mixed positional / keyword / structural refs, roundtrip
  (parse → resolve → `from_ast` → render), and assert the rendered refs follow keyword-else-positional and
  carry no `Structural`. `[dep: S3c, S3d]`

- **S3k — fuzz seeds with structural refs.** Add corpus seeds exercising the `{:atoms […]}` /
  `{:donors … :acceptor …}` / `{:site … :ligands […]}` structural forms to the `umol-ast` targets whose
  grammar admits them — `fuzz_molecule`, `fuzz_reaction`, `fuzz_reaction_span`, `fuzz_constraints` — so the
  full parse→resolve path is fuzzed on the new arm. `[dep: S3c, S3d, S3e]`

- **S3m — rename the keyword-returning metadata accessors `<kind>_id` → `<kind>_keyword`.** Deferred out
  of S3b (it cascades past the trait work). Scope: the eight getters on `MoleculeMetadata` **and**
  `ReactionMetadata` that return `Option<&str>` (a keyword), their callers (molecule / reaction / span
  render fns, the property tests), the `Metadata` trait methods, and the `define_ref!` `$accessor`. The
  `set_<kind>_id` setters follow (either renamed `set_<kind>_keyword` or retired once `From`/boundary
  projection is the only builder). **Do not touch** the AST-view `.<kind>_id()` methods
  (`atom.aromatic_system_id()`, `neighbor.bond_id()`, `StereoLigand.atom_id()`, …) — those return actual
  ids and are correctly named; distinguish by the argument (metadata getters take an id, view methods
  don't). Green (pure rename). `[dep: S3d]`

- **S3l — update `umol-ast/spec/umol-dsl-spec.md`.** Document: the structural ref forms per kind and that
  they are accepted wherever a ref is (entries, entity/relational/molecule constraints, sub-pattern
  anchors, reaction deltas); the `keyword > positional > structural` emission rule with structural
  input-only; and anonymous sub-patterns (no `:id`, no `:atom-aliases`, rejected with `InvalidValue`).
  `[dep: S3c–S3f]`

**Critical path** S2a → S2b → {S2d, S2e} → S3a → S3b → S3c → {S3d, S3f, S3h} → S3e → S3g. S3b is a green,
transparent foundation; S3c is the first breaking cut (molecule green, reaction/span red until S3d/e).

**Stereo structural refs (settled 2026-07-03):** resolved by **(site, ligand multiset)** — both part of
the resolution, matching `connecting_id` (same site + same ligand multiset, frame order not matched,
repeats significant). The namespace keys stereo elements by `(site, Vec<StereoLigand>)` (sorted), and the
`:ligands` are required in the structural form, not an assertion tacked on after a site-only lookup.
