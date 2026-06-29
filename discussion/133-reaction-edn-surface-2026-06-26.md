# 133 — Reaction EDN surface syntax (design)

Design round for an EDN surface syntax for reactions. We have two reaction AST forms and a
lossless bidirectional conversion between them (doc 131/132): `ReactionAst` (operational —
`lhs` molecule + `Deltas`) and `ReactionSpanAst` (declarative — the superimposed `L ∪_K R`
graph, one `EntitySpan` per atom/bond). This doc designs the surface for **both**, even
though implementation may land one first; the conversion means an author can write either and
render the other.

The molecule parts are already expressible (`MoleculeDsl` is the boundary type, mirrored from the
`MoleculeAst` surface in `umol-dsl-spec.md`). So the weight here is the two genuinely new
encodings: **`Delta`** (work item 1) and **`EntitySpan`** (work item 2). Everything else is
reuse.

**Spec status.** The reaction notation currently in `umol-dsl-spec.md` (the §8.4 `:lhs`/`:rhs`
two-molecule sketch) predates the `ReactionAst` / `ReactionSpanAst` model and is explicitly
non-normative on this topic. This design supersedes it; the spec follows this doc, not the
reverse, and will be revised to whatever is decided here.

## Scope

In scope:
- `ReactionDsl` ↔ `ReactionAst`; `ReactionSpanDsl` ↔ `ReactionSpanAst`.
- The `Delta` and `EntitySpan` encodings, plus the `*FieldChange` leaf encoding.
- Reuse of `MoleculeDsl` (for `lhs`), `AtomDsl`/`BondDsl` (atom/bond strings), `ValueDsl`,
  `ConstraintDsl`, the int-or-keyword ref encoding, and the bond keyword shorthands.

Out of scope:
- SMIRKS / reaction-SMILES / GML import. Those are separate boundary types (the `TableIR`
  analogue), a later concern — but see the note on the two-molecule (SMIRKS-shaped) alternative.
- `:guards` on rules (spec L3 `:guards`), and any L4 compound form.

## What we are encoding (AST recap)

```
ReactionAst        = { lhs: MoleculeAst, deltas: Deltas }
Deltas             = Vec<Delta>
Delta              = Atom(AtomDelta) | Bond(BondDelta) | Constraint(ConstraintDelta)
AtomDelta          = Add{ id, ast } | Remove{ id, ast }
                   | ModifyField{ id, change: AtomFieldChange }
                   | ModifyConstraint{ id, old: Option<AtomConstraint>, new: Option<AtomConstraint> }
BondDelta          = Add{ id, atoms: [AtomId;2], ast } | Remove{ id, atoms, ast }
                   | ModifyField{ id, change: BondFieldChange } | ModifyConstraint{ id, old, new }
ConstraintDelta    = Add(Constraint) | Remove(Constraint)        // molecule-level
AtomFieldChange    = Element|IsotopeMass|Charge|ImplicitHydrogens|LonePairs|Spin  { old, new }
BondFieldChange    = Order|Charge|Spin  { old, new }

ReactionSpanAst    = { graph, atoms: Vec<EntitySpan<AtomAst>>, bonds: Vec<EntitySpan<BondAst>>,
                       constraints: Vec<ConstraintSpan> }
EntitySpan<T>      = Unchanged(T) | Modified{ left: T, right: T } | Added(T) | Removed(T)
ConstraintSpan     = Unchanged(Constraint) | Added(Constraint) | Removed(Constraint)  // molecule-level, no Modified
```

Two framing facts that drive the encodings:

- **Delta refs are in the `lhs` frame.** `ModifyField`/`ModifyConstraint`/`Remove` name an existing
  `lhs` entity by id; `Add` introduces a *new* id (≥ the lhs count) that later deltas reference.
- **Span entries are in the union frame.** The span's `atoms`/`bonds` vectors are index-parallel
  to the union graph: `lhs` entities in their original order (a `Removed` keeps its slot), created
  entities appended (`Added`). Bonds name endpoints by union index.

Running example (hydroxide + methyl bromide → methanol + bromide; a clean SN2 that exercises
bond add, bond remove, and two charge changes):

```
lhs:    [:c "C#h3"] [:br "Br"] [:nu "O#h1#c-1"]   with bond  c–br single
change: form c–nu single; break c–br; br charge 0→-1; nu charge -1→0
```

## Settled foundations (mirror `MoleculeDsl`)

These follow the existing conventions directly; not contested, listed for completeness.

- **Boundary types.** `ReactionDsl` and `ReactionSpanDsl`, each owning ser/de via `FromEdn`/`ToEdn`
  (+ `FromStr`/`Display`), and `FromAst`/`IntoAst` with a `Ctx`. Like `MoleculeDsl` they likely
  need a private-field struct carrying surface `Metadata` (the keyword↔id bindings) alongside the
  AST, rewrapped via `from_parts`/`into_parts` — created-entity ids and lhs ids are metadata.
- **`lhs` reuse.** `ReactionDsl`'s `:lhs` is a molecule map parsed by `MoleculeDsl` unchanged.
  `ReactionSpanDsl`'s `Unchanged`/`Modified`/`Added`/`Removed` values are `AtomDsl`/`BondDsl`
  strings, parsed unchanged.
- **Refs.** Reuse the existing int-or-keyword encoding (`int` = positional index, `keyword` = id)
  resolved against `Metadata`/`EntityCounts`, exactly as bond endpoints and constraint refs do
  today.
- **Leaf values are the existing compact DSLs — the reaction layer adds no leaf tokens.** Every
  value slot in either surface is one of the established leaf DSLs, reused verbatim with its own
  EDN form; what looks like a bare int or keyword is that DSL's rendering, not a reaction-surface
  primitive:
  - atom value → `AtomDsl` string (`"C#h3"`, `"O#h1#c-1"`);
  - bond value / type → `BondDsl` — the bond string (`"1"`, `"1#a"`, `"2#R"`, …) **or** its
    keyword shorthand (`:single`/`:double`/…). `:single` is *not* an independent keyword; it
    abbreviates the bond string `"1"`, and the string form carries any extra inline state
    (aromaticity `#a`, ring `#R`, charge, spin);
  - individual field values are **not** standalone slots — `:add` / `:modify` carry the whole entity
    DSL string (work item 1), so every field (element, `#i` isotope, `#u#s` spin, `#c` charge, …)
    rides inside it in its normal `#`-prefixed compact form; no bare-int-vs-string question arises at
    the delta surface. (Value literals *inside* a `ConstraintDsl` / molecule constraint follow the
    molecule surface's value-expr-string convention; refs in them are resolved when `ReactionInput`
    is converted, like every other ref.)
  - constraint slots are two distinct leaves by slot (Work item 1): the `ConstraintDelta` payload
    is a **`ConstraintDsl`** (molecule-level — `{:connected {…}}`, `{:and […]}`, or an entity leaf
    `{:atom [<ref> {:valence …}]}`); a `ModifyConstraint` `old`/`new` is a **per-entity constraint
    single-key map** (`AtomConstraintDsl` `{:valence …}`, `BondConstraintDsl` `:aromatic`).

  So a delta/span value is always a whole entity DSL string (`:single` is its `BondDsl` keyword
  shorthand); per-field values live inside that string, not as separate EDN slots.

## Work item 1 — `Delta` encoding  **[resolved — add / remove / modify]**

`ReactionDsl` shape:

```clojure
{:lhs    {:atoms [...] :bonds [...]}
 :deltas [<delta>*]}
```

Each delta is **entity-then-op** (B): `{<entity> {<op> <payload>}}`, where `<entity>` is one of the
eight kinds (`:atom :bond :aromatic-system :multicenter-bond :dative-bond :noncovalent-bond
:stereo-atom :stereo-bond`) and `<op>` is `:add` / `:remove` / `:modify`. Molecule-level constraints
(`ConstraintDelta`) use a separate top-level `:constraint` key with `:add` / `:remove`; per-entity
constraints are not separate ops — they live inside the entity DSL handled by `:modify`. The op set
collapses the earlier `:set-field` / `:set-constraint` pair into a single `:modify` (decision 5).

### Add — the molecule-map entry form for the entity

`:add` reuses the molecule-map entry syntax (spec §4) verbatim:

```
{:atom            {:add <added-atom>}}
{:bond            {:add <added-bond>}}
{:aromatic-system {:add <added-aromatic-system>}}
{:stereo-atom     {:add <added-stereo-atom>}}
{:dative-bond     {:add <added-dative-bond>}}        ; the asymmetric one; other overlays analogously

added-atom            ::= <atom-dsl> | [<id> <atom-dsl>]
added-bond            ::= [<atom-ref> <atom-ref> <bond-dsl>]
                        | {[:id <id>] :atoms [<atom-ref> <atom-ref>] :type <bond-dsl>}
added-aromatic-system ::= {[:id <id>] :atoms [<atom-ref>*] :type <aromatic-system-dsl>}
added-stereo-atom     ::= {[:id <id>] :site <atom-ref> :ligands [<atom-ref>*] :type <stereo-atom-dsl>}
added-dative-bond     ::= {[:id <id>] :donors [<atom-ref>*] :acceptor <atom-ref> :type <dative-bond-dsl>}
```

A created entity may carry an `<id>` that later deltas/bonds reference. (`:electrons` is
absent — it moves into the entity DSL; see below.) Dative is the asymmetric case: `:donors` (a
vector) + a single `:acceptor`, not an endpoint pair (the `:donor` → `:donors` rename is recorded
as a molecule-wide change below).

### Remove — a reference (index, id, or structural)

```
{:atom {:remove <atom-ref>}}   {:bond {:remove <bond-ref>}}   ; all entities analogously

atom-ref        ::= <index> | <id>
bond-ref        ::= <index> | <id> | [<atom-ref> <atom-ref>]
aromatic-ref    ::= <index> | <id> | [<atom-ref>*]
stereo-atom-ref ::= <index> | <id> | {:site <atom-ref> [:ligands [<atom-ref>*]]}
dative-ref      ::= <index> | <id> | {:donors [<atom-ref>*] :acceptor <atom-ref>}
```

A bond / overlay can be removed by index, id, or **structurally** (a bond by its endpoint pair, an
aromatic system by its member set, a dative bond by its donors + acceptor, a stereo atom by its site
\+ optional ligands). Boundary validation must check not only that the atom-refs exist, but that they
**cover an existing** bond / aromatic system / dative bond / stereo atom. `:ligands` is optional, but
a stereo atom needs the `{:site …}` map form (not a bare ref) so its id is not confused with the atom
id of its site/focus.

### Modify — ref + the RHS entity DSL (no 1:1 old/new)

```
{:atom            {:modify [<atom-ref> <atom-dsl>]}}
{:bond            {:modify [<bond-ref> <bond-dsl>]}}
{:aromatic-system {:modify [<aromatic-ref> <aromatic-system-dsl>]}}
{:stereo-atom     {:modify [<stereo-atom-ref> <stereo-atom-dsl>]}}   ; analogously
```

The surface deliberately does **not** mirror the AST's per-field `ModifyField` /
`ModifyConstraint` (which carry both old and new) 1:1 — that creates more problems than it solves (you would write both old
and new and keep them in agreement with the lhs, e.g. `"C#c?c"` in lhs, `"C#c?c"` as old,
`"C#c(?c+1)"` as new: much typing, easy to desync). The modify payload is just the **RHS** — the
entity DSL carrying only the parts being set. The old values come from the `lhs`; deltas are **not**
independently resolvable, and need not be. This mirrors the two-stage molecule boundary (raw
`MoleculeInput` → resolved `MoleculeDsl`, ids/aliases filled at the boundary): a reaction gets a raw
input form resolved against `lhs` into `ReactionAst`, with per-field old values filled from `lhs`.
(The span's `:modify` instead carries `[left right]` complete values — it has no `lhs` to recover
from; see work item 2.)

Semantics — uniform for fields and constraints (they look identical in the entity DSL, so one
coherent representation throughout, no field/constraint special-casing):
- a **field** present in the RHS replaces; absent = unchanged. `"C"` sets element; `"#c0"` sets
  charge 0; `"#c(?c+1)"` sets charge to `?c + 1` (`?c` bound in the lhs).
- a per-entity **constraint** with a concrete value (`"#v4"`) is set/added; the old need not exist —
  consistent with `ModifyConstraint` semantics (old may be absent).
- **constraint removal** exploits that an `Undetermined` constraint is vacuous: `"#v*"` parses as
  *remove the valence constraint*.
- the RHS is a **partial** entity DSL — only the changed parts, not a complete entity (`"#c0"` has
  no leading element). It need not satisfy the complete-entity grammar, so it must **not** be forced
  through the existing complete-entity parser; this wants a sparse/partial parse mode (R4).

Running example (the SN2):

```clojure
:deltas [{:bond {:add    [:c :nu :single]}}    ; form C–Nu
         {:bond {:remove [:c :br]}}            ; break C–Br (by endpoints)
         {:atom {:modify :br "#c-1"}}          ; Br → bromide (new charge only; old from lhs)
         {:atom {:modify :nu "#c0"}}]          ; Nu → neutral
```

### `:electrons` moves into the entity DSL

The per-atom electron-contribution vector currently sits as an `:electrons` map key on
aromatic-system / multicenter entries. It is **AST value data** (an electron-count AST), unlike the
other map keys (`:id`, `:atoms`, `:site`, `:ligands`) which are structural / refs — so it belongs in
the entity DSL string (the `:type` payload) with the rest of the value content, not in the map. The
original reason for keeping it out (avoiding vector parsing inside the compact DSL) no longer holds
now that the DSL carries other non-scalar content. Consequence: the aromatic-system / multicenter
DSL grammar gains an electrons clause, which ripples to the **molecule** surface and spec (those
entries appear in molecules too) — a change beyond reactions.

This needs an in-string **vector** form (the compact DSLs are otherwise scalar / predicate shaped).
The vector value syntax is **`[<nat>(,<nat>)*]`** — square brackets, at least one element,
comma-separated, whitespace ignored — paralleling `{1,2}` for sets and `(1,2,3)` for permutations;
the whole-vector undetermined case is `*` as usual. So a per-atom electrons clause reads e.g.
`[1,1,1,2]` or `*`, position `i` matching atom `i` of `:atoms` (carried over from the former
`:electrons` key). `[…]` thus becomes the DSL's vector bracket, distinct from `{…}` (sets) and `(…)`
(grouping / permutations), reusable for any future in-string vector. The per-atom vector takes
**no `#`-tag**: it is the **head** of the aromatic-system / multicenter DSL string — leading and
unprefixed, exactly as the element is the head of the atom-string. The scalar total remains the
separate `#e<n>` predicate, so there is no collision.

### `:donor` → `:donors` (molecule-wide)

A second molecule-surface change this design surfaces: the dative key is renamed `:donor` →
**`:donors`**, always a vector (a single donor is a one-element vector). Proper pluralization — a
vector-valued key is plural (`:atoms`, `:ligands`, `:electrons`, `:donors`), whereas the spec §4
`:donor` is a singular name holding ref-or-vector. Ripples to the molecule spec §4 dative-bond-entry
and `dsl/dative.rs`; a dative bond needs ≥1 donor.

## Work item 2 — `EntitySpan` encoding (span)  **[principle accepted]**

The span surface has the **same shape as a molecule map** (`:atoms`, `:bonds`, `:constraints`, …);
each entity is either a bare molecule entry — `Unchanged` — or that entry wrapped in a single-key
**verb** map `{:<op> <entity>}`, `<op>` ∈ `:add` / `:modify` / `:remove`. **Verbs, not participles**
— the span and the delta surface name the operations the same way.

**Homoiconicity.** A plain molecule map (every entry bare) is a valid reaction span — the degenerate
(identity) reaction. This is the inverse of the MOD convention, where `(∅, ∅, R)` embeds a molecule
in left+deltas semantics; here the molecule *is* the span.

Atoms / bonds:

```clojure
{:atoms ["C#h3"                              ; Unchanged — a bare molecule atom entry
         {:add "O#h1"}                        ; Added   (right only)
         {:remove "O"}                        ; Removed (left only)
         {:modify ["Br" "Br#c-1"]}]           ; Modified — [left right]
 :bonds [[0 1 :single]                        ; Unchanged
         {:add [0 2 :single]}                 ; Added
         {:remove [0 1 :single]}              ; Removed
         {:modify [0 2 [:single :double]]}]}  ; Modified — endpoints once, [left right] value
```

**`:modify` carries both sides** — `[left right]` (atoms), `[<ref> <ref> [left right]]` (bonds, since
endpoints are frame-invariant): the span is self-contained, with no `lhs` to recover the old value
from, so this is the one place the span `:modify` differs from the delta `:modify` (which carries
only the new value).

**Per-side consistency.** Each side must be internally ref-consistent on its own, like a standalone
molecule: the **left** (`Unchanged` ∪ `Removed` ∪ `Modified.left`) and the **right** (`Unchanged` ∪
`Added` ∪ `Modified.right`) must each form a valid molecule — every ref a side uses resolves within
that side. Verified when `SpanInput` is converted (S5).

**Molecule-level constraints** (`:constraints`): same shape — a bare constraint is `Unchanged`,
`{:add <constraint>}` / `{:remove <constraint>}` add or remove it. (Constraints are a multiset; no
`:modify`.) Backed by a dedicated `constraints: Vec<ConstraintSpan>` field on `ReactionSpanAst`,
where `ConstraintSpan` is a separate enum with only `Unchanged` / `Added` / `Removed` (no `Modified`).
`to_reaction_span` tags each: `Unchanged` from `lhs`, `Added` / `Removed` from the `ConstraintDelta`s;
`left()` / `right()` rebuild each side's constraint multiset (so the degenerate molecule→span→molecule
identity keeps its constraints). This is the one collection that is *not* a payload-lift — it has no
topology, so it gets its own field and enum rather than `EntitySpan<Constraint>`.

**Every valid entity / constraint form is acceptable** in the span. It reuses the *whole* molecule
entity/constraint grammar unchanged — inline ids (`[:id {:modify […]}]`), the bond map form, every
constraint form — and adds only the `{:<op> …}` wrapper. No new leaf or entity grammar.

**AST note.** Every `MoleculeAst` entity collection is parameterized by its AST payload type `R`:
atoms/bonds as `Vec<R>`, the six overlays as `*Relation` / `*Birelation` sets that split a
*topological* part (participants + incidence over `NodeId` / `EdgeId`, ordering) from a `data: Vec<R>`
payload column. `ReactionSpanAst` is the identical structure with that payload lifted
`R → EntitySpan<R>` in every collection — topology shared, only the per-hyperedge value becomes a
span; `EntitySpan<T>` itself is unchanged. The span does this for the two vecs today; lifting the six
overlay relation sets' `data` columns is the remaining work (doc 134, item 2). The *surface* above
already generalizes (every molecule section is op-wrappable). Molecule-level `Constraints` is the lone
exception — a flat multiset with no topology — so it gets its own `Vec<ConstraintSpan>` field and a
dedicated 3-variant enum (above), not a lifted payload.

## Cross-cutting

- **Frames.** Span = union frame (lhs slots preserved for `Removed`, created appended); deltas =
  lhs frame with fresh ids for `Add`. Both reuse the int-or-keyword ref + `Metadata`
  resolution. ids are optional sugar over positional indices in both surfaces, as in molecules.
- **Rendering is faithful.** Entity values render through the existing entity DSLs verbatim — no
  reaction-level grounding, zeroing, or default-dropping. Deltas are emitted in the AST's canonical
  (post-`canonicalize`) order so DSL→AST→DSL is stable.
- **Self-description (decision 7).** The two reaction surfaces and the molecule surface are
  distinguished by keys (`:lhs`/`:deltas` vs `:atoms`/`:bonds`), but a span and a molecule share
  keys — so the EDN is **not** self-describing; the caller picks the boundary type (matches
  today's model). Optional: a discriminating tag (`#reaction`, `#reaction-span`, `#molecule`) or a
  `:kind` key for self-describing streams. Default: no tag, caller-chosen, unless we want mixed
  streams.
- **Two-molecule (SMIRKS-shaped) alternative — non-goal.** A `lhs`/`rhs` pair with shared atom
  ids (the superseded §8.4 sketch; SMIRKS `reactant>>product`) is isomorphic to the span when ids correlate
  the sides, but it is a *third* boundary type and is lossy on the atom map without disciplined
  ids. It belongs with the SMIRKS import work, not here. The span's superimposed form is the
  faithful `ReactionSpanAst` boundary and strictly carries the map.

## Decisions

The **delta** (work item 1) and **span** (work item 2) surface designs are resolved. The remaining
work (the AST overlay / span-generalization / constraint-apply items) is deferred to doc 134.

1. **Delta key shape** — **resolved: B** (entity-then-op), for coherence with the entity-keyed
   `:constraints` records, a small reused keyspace, and a single `EntityDelta`-driven path.
2. **Field vs constraint ops** — **resolved**: unified into a single `:modify` (no `:set-field` /
   `:set-constraint` split); fields and constraints look the same in the entity DSL and are handled
   uniformly, constraint removal via the vacuous `"#v*"`.
3. **Carry vs recover** — **resolved: recover** via a raw **`ReactionInput`** stage (confirmed).
   `:modify` carries only the RHS; `ReactionInput` resolves against `lhs` → `ReactionAst`, filling
   old values. Ref / id / added-id resolution is **analogous to the existing `MoleculeInput`** —
   the same pattern, not new machinery.
4. **Value string shape** — **resolved**: delta/span values are whole entity DSL strings, so
   per-field values (incl. the inline-only `Element`/`IsotopeMass`/`Spin`) ride inside; no bare
   per-field slots. `ValueAst` literals *inside* a `ConstraintDsl` / molecule constraint follow the
   molecule surface's value-expr-string convention.
5. **`:modify` key name** — **resolved: `:modify`** (`:set` rejected — the op also adds/removes
   per-entity constraints; chosen over `:update` to share one word with the span's
   `EntitySpan::Modified` for a kept-but-changed entity).
6. **Modify RHS partial parse** — **resolved**: the `:modify` RHS is a *partial* entity DSL (only the
   changed parts), parsed sparsely (R4), not via the complete-entity parser.
7. **Electrons vector DSL** — **resolved**. Vector syntax `[<nat>(,<nat>)*]` (square brackets, ≥1,
   comma-separated, whitespace ignored; undetermined `*`), `[…]` paralleling `{…}`/`(…)`; position
   `i` matches atom `i` of `:atoms`. The vector is the **untagged head** of the aromatic-system /
   multicenter DSL (like element is the atom-string head); no `#`-tag, the scalar total stays the
   separate `#e<n>`. Ripples to the molecule DSL + spec.
8. **Remove cross-check** — **resolved**: a `:remove` ref must resolve to an existing entity, but
   that is the same ref-existence validation `MoleculeInput` already performs for atom-refs in bonds
   — reuse it, not new machinery.
9. **Span (`EntitySpan`) encoding** (work item 2) — **resolved**: molecule shape; each entity is
   bare (`Unchanged`) or wrapped `{:<op> <entity>}` with **verbs** `:add`/`:modify`/`:remove`; every
   valid molecule entity/constraint form is acceptable; a plain molecule is a degenerate span.
10. **Span `:modify` value form** — **resolved**: `[left right]` (the span is self-contained, so it
    carries both sides — unlike the delta `:modify`, which carries only the new value).
11. **Boundary vs constituent** — `Delta` is a constituent of `ReactionDsl` (recover-from-`lhs`
    rules out a standalone `DeltaDsl`); the span is self-contained, so `ReactionSpanDsl` is its own
    boundary type. **Resolved.**
12. **Self-describing tag** — **resolved: none** for now (caller picks the boundary type, as with
    molecules); revisit if mixed/heterogeneous streams arise.
13. **Dative bonds** — **resolved** (the asymmetric case: multi-donor + single acceptor, not an
    endpoint pair). `:add` = `{[:id <id>] :donors [<atom-ref>*] :acceptor <atom-ref> :type
    <dative-bond-dsl>}`; `:remove` ref = `<index> | <id> | {:donors [<atom-ref>*] :acceptor
    <atom-ref>}`. The `:donor` → `:donors` rename is a molecule-wide change, recorded separately.
14. **Span constraints** — **resolved (option a)**: `ReactionSpanAst` gets a
    `constraints: Vec<ConstraintSpan>` field, `ConstraintSpan` a dedicated enum `Unchanged` / `Added`
    / `Removed` (no `Modified` — molecule constraints are a by-value multiset). The lone collection
    that is not a payload-lift (no topology, so its own field + enum, not `EntitySpan<Constraint>`).
    `to_reaction_span` populates it (`Unchanged` from `lhs`, `Added`/`Removed` from the
    `ConstraintDelta`s) rather than dropping `Delta::Constraint`, preserving molecule→span
    homoiconicity for constraints. (Operational `apply_at` is the separate 134 item 1.)

**Deferred to doc 134 (implementation items, not surface-design questions):**
- **Overlay span generalization** — the six overlay relation sets' `data` columns need lifting to
  `EntitySpan<…>`; `EntitySpan<T>` is generic, so the work is the container, not the type (134 item 2).
- **Surface vs AST scope** — the surface spans all eight entities, but the AST `Delta` is
  `Atom`/`Bond`/`Constraint` only; overlay-entity deltas need increment-2 AST support (134 item 2).
- **Molecule-constraint reaction application** — `apply_at` currently drops `Delta::Constraint`;
  closing it needs a match-based constraint ref-remap (134 item 1).

## Implementation plan — delta surface (`ReactionDsl` ↔ `ReactionAst`)

Increment-1 entities (atom / bond / molecule-constraint); overlay-entity deltas deferred
(doc 134 item 2). Types and signatures below are approximate (shapes + API, no bodies).

**Module placement.** New module `umol-ast/src/dsl/reaction.rs`, registered in
`umol-ast/src/dsl.rs` (`pub(crate) mod reaction;` + `pub use reaction::{ReactionDsl,
ReactionMetadata};`). The partial-DSL parsers (R4) are free fns in the existing `dsl/atom.rs` /
`dsl/bond.rs`, next to `parse_atom` / `parse_bond`. No other new modules.

**D0 — Prerequisites** (molecule-DSL; no new types): **Done**
- 0.1 — relocate the per-atom `electrons: ElectronCountsAst` field from the `:electrons` map key to a
  **mandatory head** of the aromatic-system / multicenter-bond DSL string. No AST change —
  `ElectronCountsAst` stays `Lit(Vec<i64>)`; the head parses integers (the existing domain, not
  enforcing non-negative). Sub-points:
  - 0.1a — aromatic-string head parser (`dsl/aromatic.rs`): winnow `electron_head ->
    ElectronCountsAst` — mandatory, `*` → `Undetermined`, `[n(,n)*]` → `Lit` (whitespace ignored, ≥1
    element); wired into `aromatic_system` before the predicate loop; missing/malformed head → new
    `ParseError` variant(s). `raise`/`lower` keep passing `electrons` through (no defaulting).
  - 0.1b — aromatic-string `Display` (`dsl/aromatic.rs`): emit the head first — `*` for `Undetermined`,
    `[n,n,…]` for `Lit`.
  - 0.1c — multicenter-string parser + `Display` (`dsl/multicenter.rs`): same as 0.1a/b.
  - 0.1d — drop the `:electrons` map key (`dsl/molecule.rs`): remove the streaming + tree reads and the
    render for both entities, plus the now-unused `read_`/`parse_`/`render_electron_counts`; the `:type`
    DSL carries electrons now.
  - 0.1e — tests: string parse/`Display`/round-trip (concrete `[…]` heads, a dedicated `*` case, and
    invalid forms — empty `[]`, unmatched `[`/`]`, non-numeric, empty/trailing comma, missing head);
    `molecule/tests.rs` + fixtures move electrons into `:type` (length-matched 1s for synthetic, correct
    counts for real) and drop `:electrons`.
  - 0.1f — spec: drop `:electrons` from the §4 entries, rewrite ¶142, update §7.10/§7.11 (mandatory head
    grammar `[<n>(,<n>)*] | *`; `""` no longer valid for these strings; concrete examples).
- 0.2 — `:donor` → `:donors` (molecule-wide): rename the dative-bond-entry key and make it **always a
  vector** (a single donor is a 1-vector). Touches `dsl/molecule.rs` + spec only — *not* `dsl/dative.rs`
  (donor/acceptor live in the entry map, not the dative-string). Sub-steps:
  - 0.2a — streaming reader (`read_dative_bond_entry`): rename the `"donor"` arm to `"donors"`, parse
    the value as a vector only (`read_vec`, drop the `peek_byte` scalar/vector branch);
    `missing("donors", …)`.
  - 0.2b — tree reader (`parse_dative_bond_entry`): `required_key(m, "donors", …)`, parse via
    `parse_vec` only (drop the `Edn::Vector(_)`-vs-scalar match).
  - 0.2c — render (`render_dative`): emit `:donors` always as a vector (drop the `len() == 1`
    scalar-ref branch); key `"donors"`.
  - 0.2d — reject an empty `:donors` (a dative bond needs ≥1 donor); check whether
    `dative_structure_check` already covers it, else error at parse.
  - 0.2e — tests/fixtures: `:donor X` → `:donors [X]`, `:donor [X Y]` → `:donors [X Y]` across
    `molecule/tests.rs`, `validate.rs`, `macros.rs`, `benches/fixtures.rs`, examples, and other crates.
    Add an error test: a bare-int `:donors 1` (non-vector) is a parse error (streaming + tree).
  - 0.2f — spec §4: rename `:donor` → `:donors` in the `dative-bond-entry` production and the dative
    paragraph; state it is always a vector (≥1 donor).
- 0.3 — dative **donor relational constraints → multi-atom parity** (donors are a set, so mirror the
  aromatic/multicenter 5-variant set pattern; the **acceptor stays single**). Role stays in the name
  (`Donor`, not `Atom`). Sub-steps:
  - 0.3a — `ast/constraint/relational.rs` enum: keep `DativeBondDonor { bond, atom }` (semantics →
    **membership**, "`atom` is one of the donors" = `Contains` analogue; fix doc); add
    `DativeBondDonors { bond, atoms }` (exact set = `Atoms`), `DativeBondContainsAllDonors { bond,
    atoms }` (⊇ = `ContainsAll`); replace `DativeBondDonorSatisfies` with `DativeBondAllDonors { bond,
    predicate }` (∀ = `AllAtoms`) and `DativeBondAnyDonor { bond, predicate }` (∃ = `AnyAtom`). Acceptor
    variants unchanged.
  - 0.3b — `relational.rs` `remap`: arms for the new variants (vec-atom map for `Donors` /
    `ContainsAllDonors`, single-atom for `Donor`, predicate-only for `All`/`AnyDonors`), mirroring the
    aromatic arms.
  - 0.3c — `relational.rs` canonicalize + any exhaustive match (e.g. the `other @ (…)` arms): cover the
    new variants.
  - 0.3d — `dsl/relational.rs` (+ `dsl/constraint.rs` dispatch if needed): parse/render keys
    `:dative-bond-donors` (exact), `:dative-bond-donor` (membership), `:dative-bond-contains-all-donors`
    (⊇), `:dative-bond-all-donors` (∀), `:dative-bond-any-donor` (∃) — replacing
    `:dative-bond-donor-satisfies`; mirror the aromatic/multicenter key handling.
  - 0.3e — tests: `relational.rs`, `dsl/relational.rs`, `molecule/tests.rs` (replace the
    `:dative-bond-donor-satisfies` case, add cases for the new variants), `property/strategies.rs`.
  - 0.3f — spec §8: replace the single donor-satisfies form with the five donor forms, mirroring the
    aromatic/multicenter relational forms.

## Reaction build plan (R1–R12)

Bottom-up: each chunk is a coherent semantic context placed where its dependencies are satisfied,
ending at a compilable, tested checkpoint. The boundary-type traits (`FromAst`/`IntoAst`/`FromEdn`/
`ToEdn`/`FromStr`/`Display`) stay with the boundary type but split across R3/R8/R10 by dependence.
Within a chunk the entity order is atom → bond → constraint (bond resolution needs atom ids). The tree
may go red between work items inside a chunk; each chunk ends green and tested. Prereqs in brackets.

**R1 — Raw + metadata types** [—]: **Done**
- R1a — create `dsl/reaction.rs`; add `pub(crate) mod reaction;` to `dsl.rs` (between `noncovalent` and
  `refs`); imports mirror `dsl/molecule.rs`.
- R1b — `#[derive(Clone, Debug, Default, PartialEq, Eq)] struct ReactionMetadata { lhs: MoleculeMetadata,
  atom_ids: IndexMap<AtomId, String>, atom_aliases: BiBTreeMap<String, Box<AtomDsl>>, bond_ids:
  IndexMap<BondId, String> }` — created-entity id→name + atom-alias bindings (the reaction surface admits
  the `[:C "C#h3"]` alias notation), mirroring `MoleculeMetadata`; `pub use reaction::ReactionMetadata;`.
- R1c — `ReactionMetadata` accessors, mirroring `MoleculeMetadata`: `lhs()`, `atom_id(id)` /
  `bond_id(id)` (id→name), `set_atom_id` / `set_bond_id`, and the atom-alias accessors (`atom_alias_for`,
  `has_atom_alias`, `iter_atom_aliases`, …).
- R1d — `#[derive(Debug)] enum DeltaInput` (8 variants: `AtomAdd { id, value }`, `AtomRemove`,
  `AtomModify`, `BondAdd { id, atoms, value }`, `BondRemove`, `BondModify`, `ConstraintAdd`,
  `ConstraintRemove`; `id: Option<String>` is the created entity's `:id` name; refs unresolved,
  `:modify` RHS a partial AST; reuses `AtomRefDsl`/`BondRefDsl`/`ConstraintDsl`).
- R1e — `#[derive(Debug)] struct ReactionInput { lhs: MoleculeInput, deltas: Vec<DeltaInput> }`
  (`pub(crate)`; struct only, `into_ast` in R7).
- *Checkpoint:* compiles.

**R2 — `ReactionDsl` boundary shell** [R1]: **Done**
- R2a — `#[derive(Clone, Debug, Default)] struct ReactionDsl { ast: ReactionAst, metadata: ReactionMetadata }`
  (private fields) + doc; `pub use reaction::ReactionDsl;`.
- R2b — inherent API: `from_parts`, `into_parts`, `ast`, `metadata`.
- R2c — derive `PartialEq` + `Eq`.
- *Checkpoint:* type exists + exported.

**R3 — AST↔DSL conversion** [R2]: **Done**
- R3a — `FromAst<ReactionAst>` (`Ctx = ReactionDefaults` — atom/bond only): `lhs` via
  `MoleculeDsl::from_ast` (with a `MoleculeDefaults` carrying the reaction's atom/bond policy, no-op for
  the rest); `deltas` via `lower_delta` (per-delta `lower_atom`/`lower_bond` on `Add`/`Remove` payloads);
  `metadata = ReactionMetadata::default()`.
- R3b — `IntoAst<ReactionAst>`: inverse — `lhs` via `MoleculeDsl::into_ast`, `deltas` via `raise_delta`.
- R3c — test: DSL→AST→DSL round-trip under `ReactionDefaults::ground()` (wrap a `ReactionAst` as
  `ReactionDsl`, `into_ast` then `from_ast`, assert `.ast()` equal) — covers lhs + atom add + bond add +
  modify + constraint.
- *Checkpoint (tested):* programmatic round-trip, no EDN.

**R4 — Partial entity boundary types** [—] **Done**
  (`dsl/atom.rs`, `dsl/bond.rs`; reuse existing field parsers).  Each is a `*Dsl` boundary type owning
  its serde via traits — no free `read_edn_`/`render_`/`fmt_` partial fns:
- R4a — `PartialAtomDsl(pub AtomAst)`: `partial_atom` combinator (element optional ⇒ `Undetermined`,
  unspecified fields `Undetermined`) + `parse_partial_atom -> PartialAtomDsl`; `FromStr` (→
  `parse_partial_atom`), `Display` (complete render minus the element when `Undetermined`), `FromEdn`
  (`Edn::Str` → `s.parse`, recoded to `DeError`), `ToEdn` (`Edn::Str(to_string())`).
- R4b — `PartialBondDsl(pub BondAst)`: same shape, order optional.
- R4c — tests: parse + `FromEdn`/`ToEdn` round-trips (str + Edn), all-`Undetermined` and sparse cases,
  duplicate-predicate error.
- *Checkpoint (tested):* partial round-trips.

**R5 — Delta parsing** [R1, R4] (`dsl/reaction.rs`) **Done**
- R5a — dispatch: `read_delta_input` (streaming) / `parse_delta_input` (tree) route the outer
  entity keyword to per-entity `{read,parse}_delta_<entity>_input`, each owning its op-keyword
  routing (`:add` / `:remove` / `:modify`); arms are R5b–d. Content-independent EDN helpers
  (`single_key_map`, `read_vec`/`read_map`/`read_single_key_map_header`/`consume_single_key_map_close`,
  `parse_vec`, and the new tree `parse_single_key_map`) live in `dsl/edn_utils.rs`.
  `DeltaInput::Atom{Add,Remove,Modify}` etc. carry `AtomRef`/`BondRef` (index | keyword id, from
  `dsl/refs.rs`); `Add` carries the unresolved molecule entry, resolved in R7.
- R5b — atom arms (in `{read,parse}_delta_atom_input`): `:add` → `AtomAdd(AtomEntryInput)` (reuses
  `molecule::{read,parse}_atom_entry`; a bare atom, `[<id> <dsl>]`, or alias spec — aliases kept
  unresolved, not rejected); `:remove` → `AtomRemove(AtomRef)`; `:modify [<ref> <dsl>]` →
  `AtomModify(AtomRef, PartialAtomDsl)`.
- R5c — bond arms (in `{read,parse}_delta_bond_input`): `:add` → `BondAdd(BondEntryInput)` (reuses
  `molecule::{read,parse}_bond_entry`; vector or map+`:id` form); `:remove` → `BondRemove(BondRef)`;
  `:modify [<ref> <dsl>]` → `BondModify(BondRef, PartialBondDsl)`. The structural bond-ref is
  deferred to doc 134 item 3 (lands with the overlay entities); refs are index|id for now.
- R5d — constraint arms (in `{read,parse}_delta_constraint_input`): `:add` / `:remove` →
  `Constraint{Add,Remove}(ConstraintDsl)` value-based — molecule-level constraints are keyless (flat
  multiset), so removal names the whole constraint, not a ref.
- R5e — tests: `test_{parse,read}_delta_input_{atom,bond,constraint}` cover every `DeltaInput` shape,
  tree + streaming. (Landed per-arm in R5b–d.)
- *Checkpoint (tested):* every `DeltaInput` shape parses. **Done.**

**R6 — Top-level parse** [R1, R5] (`dsl/reaction.rs`) **Done**
- R6a — `read_reaction_input` / `parse_reaction_input`: `:lhs` via the molecule input parser
  (`read_/parse_molecule_input`, now `pub(super)`) — **mandatory** (`MissingField` if absent, like a
  molecule must give `:atoms`); `:deltas` looped via R5 (`read_vec`/`parse_vec` over the R5 dispatch);
  `:atom-aliases` (reaction-level, the molecule's flat name/atom-string vector form via
  `read_/parse_atom_aliases`) → `ReactionInput.atom_aliases`; unknown key → `UnknownField`.
- R6b — tests: `test_{parse,read}_reaction_input` — a full reaction map (lhs + atom/bond/constraint
  deltas) → raw `ReactionInput`.
- *Checkpoint (tested):* reaction map parses. **Done.** (The dispatch chain stays dead-code until R8
  wires `read_/parse_reaction_input` into `ReactionDsl`'s `FromEdn`/`from_edn_str` — the public root.)

**R7 — Resolution** [R1]: `ReactionInput::into_ast → (ReactionAst, ReactionMetadata)` **Done**
(`dsl/reaction.rs`). Work items in dependency order:
- R7a — **Done.** `into_ast(self) -> Result<(ReactionAst, ReactionMetadata), ParseError>` (no Ctx,
  mirrors `MoleculeInput::into_ast` — display-form AST; raise is the separate `IntoAst` step). Resolves
  `lhs`, seeds `ReactionMetadata.lhs`. Single in-order pass over `:deltas` (no forward refs).
- R7b — **Done.** Atom adds → fresh `AtomId` (lhs atom count + order), `:id` recorded via
  `set_atom_id`, `check_id_disjoint` against the id/alias namespace. `:add` spec resolves: `Bare` →
  `AtomDsl.0`; `Alias` → the **union** alias table (lhs aliases ∪ reaction `:atom-aliases`, bijective,
  collisions error), which is also seeded into `ReactionMetadata.atom_aliases`. Test:
  `test_reaction_input_into_ast`.
- R7c — **Done.** Atom namespace = lhs atom ids ∪ atom-add ids (`atom_id_to_idx` seeded from lhs,
  grown per add); atom refs resolve via `AtomRef::resolve` (unknown / out-of-range → error).
- R7d — **Done.** Bond adds → fresh `BondId` (lhs bond count + order), `:id` via `set_bond_id`,
  disjoint from the atom/alias/bond-id namespace; endpoint atom refs resolve against the running atom
  count (may name same-reaction added atoms). `BondDelta::Add { id, atoms: [a, b], ast }`.
- R7e — **Done.** Bond namespace = lhs bond ids ∪ bond-add ids (`bond_id_to_idx`); bond refs resolve
  via `BondRef::resolve`.
- R7f — **Done.** Removes (atom + bond) recover the entity AST from `lhs` by id (`lhs[id].clone()`;
  bond endpoints via `lhs.bond(id).atom_ids()`) → `Remove`. Removing an entity added in the same
  reaction is an error (recover-from-lhs cannot reach it).
- R7g — **Done.** Modifies (atom + bond): `lhs[id].update(&rhs)` overlays the partial RHS onto the lhs
  value (`AtomAst::update` / `BondAst::update`), then `EntityPatch::diff` emits the
  `ModifyField` / `ModifyConstraint` deltas (a possibly-symbolic `ValueAst` stored unevaluated; an
  `Undetermined` constraint = removal). Modifying a same-reaction-added entity is an error.
- R7h — **Done.** Constraint add/remove → `ConstraintDelta::Add` / `Remove`, each wrapping
  `ConstraintDsl::into_ast(&counts, &namespace)` (the existing full constraint resolver — all entity
  kinds + relational + molecule + And/Or/Not). Constraint `Remove` needs no lhs recovery; the DSL
  carries the whole constraint.
- **Unified resolution namespace** (replaces the per-entity `IndexMap`s across R7c–h): a running
  `EntityCounts` (`EntityCounts::from_ast(&lhs)` + `allocate_atom` / `allocate_bond` mutators) and a
  running `MoleculeMetadata` (`metadata.lhs().clone()`), both seeded from lhs and grown in delta order
  as entities are defined; every ref — entity and constraint — resolves against this pair via
  `Ref::into_ast(count, &namespace)` (linear id scan; O(count), parse-time). Disjointness uses
  `MoleculeMetadata::contains_id`. `metadata` (the output `ReactionMetadata`) keeps created ids
  separate for roundtrip. Generalizes to doc-134 overlays (each overlay add: `allocate_<kind>` +
  `set_<kind>_id` + `metadata.set_<kind>_id`, refs resolve against the same pair). Sub-linear is
  deferred: would make `MoleculeMetadata`'s id maps `BiBTreeMap` (O(log n) both ways), a separate
  shared-type change.
- R7i — **Done.** Assembles `Deltas` + `ReactionMetadata` (`ReactionAst { lhs, deltas }`).
- R7j — **Done.** `test_reaction_input_into_ast{,_alias_union}`,
  `..._atom_{remove,remove_error,modify}`, `..._bond_{add,remove,remove_error,modify}`,
  `..._constraint_{add,remove,added_atom_ref}` (the last exercises a constraint ref naming a
  same-reaction-added atom), plus `test_{atom,bond}_ast_update`.
- *Checkpoint (tested):* hand-built `ReactionInput` resolves (atoms, bonds, constraints).

**R8 — `FromEdn` path** [R2, R6, R7]: **Done**
- R8a — **Done.** `FromEdn` for `ReactionDsl`: `from_edn` → `parse_reaction_input` (tree) →
  `into_ast` → `from_parts`; `from_edn_str` → `read_reaction_input` (streaming) + `expect_eof` →
  `into_ast` → `from_parts`. `ParseError` mapped via `DeError::Custom`. Mirrors `MoleculeDsl`. Wires
  up the previously-unused reaction parse fns.
- R8b — **Done.** `FromStr for ReactionDsl` (`Err = ParseError`), delegating `from_edn_str` and
  mapping `EdnError` → `ParseError::EdnParse`. Mirrors `MoleculeDsl`.
- R8c — **Done.** `test_reaction_dsl_from_edn` (EDN → `ReactionDsl`, asserting the resolved
  `ReactionAst` for atom-modify and atom-add+bond-add), `test_reaction_dsl_from_edn_str_from_edn_parity`
  (tree vs streaming agree across modify / bond-modify+constraint / add cases), and
  `test_reaction_dsl_from_str` (`FromStr` ≡ `from_edn_str`).
- *Checkpoint (tested):* EDN → AST end to end.

**R9 — Render** [R1, R4] (`dsl/reaction.rs`): **Done**
- R9a — **Done.** `render_deltas` atom ops (`dsl/reaction.rs`): add → `{:atom {:add <entry>}}`
  (reaction-local `render_atom_entry`: created-id frame + lhs∪reaction alias), remove →
  `{:atom {:remove <ref>}}` (reaction-local `render_atom_ref`: lhs frame only — the two id spaces stay
  disjoint), modify → coalesce consecutive same-id `ModifyField`/`ModifyConstraint` into one partial
  `AtomAst` (fields set to `new`; constraint set adds `new`; constraint removal adds
  `old.as_undetermined()`), rendered `{:atom {:modify [<ref> <PartialAtomDsl>]}}`; `old` dropped.
  Constraint-removal `#v*`: the partial atom is the one place an undetermined constraint is **not**
  vacuous — `PartialAtomDsl`'s render now emits `#<tag>*` for an undetermined constraint (full
  `AtomDsl`/molecule still drops vacuous, unchanged). New AST primitive
  `AtomConstraint::as_undetermined()` (vacuous form, keeps `RingMembership` scope). Tests:
  `test_render_deltas`, `test_atom_constraint_as_undetermined`.
- R9b — **Done.** `render_deltas` bond ops: add → `{:bond {:add <entry>}}` (`render_bond_entry`:
  `[<a> <b> <type>]`, or `{:id :atoms :type}` when the bond has an id; `<type>` is full `BondDsl`
  so a plain order renders its keyword shorthand `:single`/`:double`); remove →
  `{:bond {:remove <ref>}}` (`render_bond_ref`, lhs frame); modify → coalesce same-id
  `ModifyField`(Order/Charge/Spin)/`ModifyConstraint` into a partial `BondAst` →
  `{:bond {:modify [<ref> <PartialBondDsl>]}}` (partial renders order as a bare string, `"2"`).
  Endpoints resolve against the **union** namespace (`render_atom_endpoint`: created ∪ lhs) — a bond
  may attach to a same-reaction atom, unlike a delta target ref. `BondConstraint::as_undetermined`
  added; `PartialBondDsl` renders undetermined constraints as `#<tag>*`. Tests:
  `test_render_deltas_bond`, `test_bond_constraint_as_undetermined`.
  - **Aromatic-as-lattice follow-up** (fixes the original "`Aromatic` is a flag, not removable"
    deficit): introduced `BooleanAst` (`Undetermined`/`Lit(bool)`, full lattice) + `BooleanDsl`
    (`+`/`` → true, `!` → false, `*` → undetermined; EDN `true`/`false`/`:undetermined`). Rewired
    `BondConstraint::Aromatic(BooleanAst)` and `DativeBondConstraint::Aromatic(BooleanAst)`: now
    value-aware in `is_undetermined`/`as_undetermined`/`canonicalize`/`Lattice` (`aromatic() ->
    BooleanAst`), so a bond/dative aromatic **is** removable via `#a*` like every other constraint.
    Surface: `#a`/`#a+` = true, `#a!` = not-aromatic, `#a*` = undetermined; the constraint EDN became
    `{:aromatic <bool>}` (the bare `:aromatic` keyword shorthand for the **constraint** is gone; the
    bond-keyword shorthand `:aromatic` for a whole bond stays). Made vacuous-`*` elision uniform in
    the per-constraint formatters: `#a*`, `#C*` (cis-trans) and `#T*` (tetrahedral) all elide on full
    render, while `#a+`/`#C+`/`#T+` (= `Stereo(Undetermined)`, "is a stereocenter, config unknown")
    are preserved. Cross-crate updates: `kekulizer`, `aromaticity`, `counts`, `table_ir/raise`.
- R9c — **Done.** `render_deltas` constraint ops: `ConstraintDelta::Add`/`Remove` →
  `{:constraint {:add|:remove <ConstraintDsl>}}` via `ConstraintDsl::from_ast(c, combined)`. Its
  entity refs span lhs ∪ created, so they resolve against `ReactionMetadata::combined_metadata()`
  (merges the lhs molecule metadata with the created-entity id bindings) — built once, lazily, only
  when a constraint delta is present. `render_deltas` renders deltas in their stored order (no
  ToEdn-side canonicalization). Test: `test_render_deltas_constraint` (molecule `:connected` +
  entity-leaf `{:atom [:o …]}` referencing a created atom + `:remove`).
- R9d — **Done.** `render_reaction_edn(&ReactionAst, &ReactionMetadata)` builds the `{:lhs …
  :deltas … :atom-aliases …}` map: `:lhs` via the molecule renderer (`render_molecule_edn`, made
  `pub(super)`), `:deltas` via `render_deltas`, then `:atom-aliases` (reaction-level, only when
  present) — aliases last, matching the molecule surface. Key order is cosmetic (both parse paths
  collect keys order-independently; alias resolution is deferred to `into_ast`).
  Test: `test_render_reaction_edn` (render → reparse preserves the `ReactionDsl`, covering a modify
  and a reaction-level alias used by an `:add`).
- R9e — **Done.** Per-op render table tests (`render_deltas`, shared `meta` fixture):
  `test_render_deltas_atom` (add / remove / modify-field / modify-set-constraint /
  modify-remove-constraint `#v*` / coalesced `#c-#v*`), `test_render_deltas_bond` (add / remove /
  modify-field / modify-constraint `#a`), `test_render_deltas_constraint` (molecule / entity-leaf via
  combined metadata / remove); plus `test_render_reaction_edn` (top-level render→reparse).
- *Checkpoint (tested):* AST renders to EDN.

**R10 — `ToEdn` path + AST routing** [R2, R8, R9]: **Done**
- R10a — **Done.** `impl ToEdn for ReactionDsl` delegates to `render_reaction_edn(&self.ast, &self.metadata)`;
  `test_reaction_dsl_from_edn_to_edn_roundtrip` exercises the trait path (from_edn → `to_edn` → from_edn).
- R10b — **Done.** `impl Display for ReactionDsl` → `write!("{}", self.to_edn())`.
- R10c — **Done.** `ReactionAst` routing through `ReactionDsl`: `FromEdn`/`FromStr` discard metadata
  (`dsl.into_parts().0`); `Display` → `to_edn`; `ToEdn` renders with `ReactionMetadata::default()`
  (positional lhs + refs, no ids/aliases), mirroring `MoleculeAst::to_edn`.
- R10d — **Done.** `test_reaction_dsl_from_edn_to_edn_roundtrip` (+molecule-`:connected` and entity-leaf
  `:atom` constraint deltas); `test_reaction_ast_to_edn` (metadata-free positional round-trip: lhs modify,
  created-atom-as-bond-endpoint, molecule constraint, entity-leaf constraint).
- *Checkpoint (tested):* full round-trip. **Done.**

**R11 — Spec** [R10]:
- R11a — **Done.** New normative **§8 Reaction map** in `umol-dsl-spec.md`: `reaction-map`/`delta`/
  `atom-delta`/`bond-delta`/`constraint-delta` grammar + normative paras (reference frames lhs vs
  lhs∪created, no-forward-refs, create-vs-edit, `:modify` partial w/ `#tag*` removal + coalescing,
  `:constraint` deltas, `:atom-aliases` resolved-after/emitted-last, serialization order/positional).
  Examples section renumbered §8→§9; §1 "EDN and rules" cross-refs §8.
- R11b — **Done.** §9.4 (was §8.4) rewritten from the obsolete `:lhs`/`:rhs` whole-molecule rewrite to
  the `:lhs` + `:deltas` surface: a `:modify` field edit, and a grow-and-constrain reaction (atom add →
  bond add → `:connected`). Both examples verified to parse + round-trip.
- R11c — **Done.** Bond `#a` (§7.5) and dative `#a` (§7.12) recast from no-payload flag to a boolean
  constraint: `#a`/`#a+`=true, `#a!`=false, `#a*`=undetermined (elided); EDN `:aromatic` → `{:aromatic
  <boolean>}`, `<boolean>` = `true | false | :undetermined`; added the `boolean` production to §7.9;
  bond/dative-constraint-form updated; §4 dative prose "flag"→"constraint". Atom `#a` (π-contribution,
  §5.4/§7.3) untouched — distinct namespace.
- *Checkpoint:* spec updated. **Done.**

**R12 — Property & fuzz tests** [R10] (`proptest`, feature-gated; `cargo-fuzz` targets): **Done**
Reuse the existing harness — `umol-ast/tests/property/{delta,reaction}.rs` and `umol-ast/fuzz/`.
- R12a — **Done.** Patch algebra (`EntityPatch`, now `pub`) for `AtomDelta` / `BondDelta` in
  `property/delta.rs`: `test_{atom,bond}_delta_diff_apply` (`apply(left, diff(left, right)) == right`)
  and `test_{atom,bond}_delta_diff_identity` (`diff(x, x)` empty; empty-diff apply is the identity),
  via local `apply_{atom,bond}_diff` helpers.
- R12b — **Done.** Partial entity AST↔DSL string round-trip in `property/entity.rs`:
  `test_partial_{atom,bond}_dsl_display_from_str_roundtrip` (`parse(display(x)) == x`).
- R12c — **Done.** Partial entity EDN round-trip in `property/entity.rs`:
  `test_partial_{atom,bond}_dsl_to_edn_from_edn_roundtrip` (`from_edn(to_edn(x)) == x`). The partial's
  EDN form is a string leaf; tree-vs-streaming parity is **not** asserted here — see the deferred note.
- R12d — **Done.** Delta EDN round-trip in `property/reaction.rs`:
  `test_reaction_ast_edn_roundtrip_stable` (render → parse reaches a fixpoint over `reaction_strategy`,
  covering atom/bond add / remove / modify-field). modify-constraint ops are covered by the R9e unit
  tests, not generated here.
- R12e — **Done.** `umol-ast/fuzz/fuzz_targets/fuzz_reaction.rs` (registered in `fuzz/Cargo.toml`):
  drives `ReactionDsl::from_edn_str` (streaming) and `read_string` + `from_edn` (tree); no panic, and
  parse-or-error parity (asserts equality when both succeed). Not executed here (needs `cargo-fuzz`).
- *Checkpoint (tested):* algebra and round-trip laws hold under generation (3854 lib + 95 property
  tests pass); fuzz target compiles into the harness.

**Deferred — streaming `from_edn_str` for the constituent `*Dsl` types.** Per dsl-serialization, every
EDN-shaped `*Dsl` should override `FromEdn::from_edn_str` to drive its hand-written streaming `read_*`
(never the default, which delegates to the tree `read_string` + `from_edn`). Today only the top-level /
full-entity types do (`MoleculeDsl`, `ReactionDsl`, `AtomDsl`, `BondDsl`, `DativeBondDsl`, … via
`read_subgrammar_all`); the constituent boundary types still inherit the tree-delegating default:
`ValueDsl`, `BooleanDsl`, `RingMembershipDsl`, `AromaticValenceDsl`, `MulticenterValenceDsl`,
`AtomConstraintDsl`, `BondConstraintDsl`, `DativeBondConstraintDsl`, `AromaticSystemConstraintDsl`,
`MulticenterBondConstraintDsl`, `NoncovalentBondConstraintDsl`, `SubPatternAnchorDsl`,
`MoleculeConstraintDsl`, `RelationalConstraintDsl`, `ConstraintDsl`, `ConstraintsDsl`,
`StereogenicityDsl`, `TopicityDsl`, `StereoCosetDsl`, `PartialAtomDsl`, `PartialBondDsl`, the stereo
constraint `*Dsl` macro, and the ref `*Dsl` macro. `fuzz_constraints` calls
`ConstraintDsl`/`ConstraintsDsl::from_edn_str` directly, so its "streaming" arm is currently the tree
path. Wiring notes: most have a single-arg standalone `read_<type>_dsl` to delegate to (trivial);
`read_topicity`/`read_stereogenicity`/the stereo constraint readers are private and would need
`pub(super)`; `read_molecule_constraint_dsl(de, key)` and `read_relational_constraint_dsl(de, key)`
have **no** standalone reader (the `{`+dispatch-key is consumed by the umbrella `read_constraint_dsl`),
so they need a thin standalone reader or a different structure — design decision, not mechanical.

## Implementation plan — span surface (`ReactionSpanDsl` ↔ `ReactionSpanAst`)

Increment-1 entities (atom / bond / molecule-constraint); overlay span columns deferred (doc 134
item 2). Span `:modify` carries complete `[left right]` values (the span is self-contained), so this
reuses the complete molecule entry parsers — no partial parser, unlike the delta surface.

**Module placement.** New module `umol-ast/src/dsl/reaction_span.rs`, registered in `dsl.rs`
(`pub(crate) mod reaction_span;` + `pub use reaction_span::ReactionSpanDsl;`). AST changes land in
`ast/delta.rs` and `ast/reaction_span.rs`.

## Reaction span build plan (C1–C12)

Bottom-up; mirrors the reaction R-plan. Differences: the span needs net-new AST work first (C1/C2),
there is no partial parser (span `:modify` carries complete `[left right]` values), and metadata reuses
`MoleculeMetadata` (no new type). The tree may go red between work items inside a chunk; each chunk ends
green and tested. Prereqs in brackets.

**C1 — AST: constraint span** [—] (`ast/delta.rs`, `ast/reaction_span.rs`): **Done**
- C1a — `#[derive(Clone, Debug, PartialEq, Eq)] enum ConstraintSpan { Unchanged(Constraint),
  Added(Constraint), Removed(Constraint) }` (no `Modified` — molecule constraints are a by-value
  multiset) + `left()` / `right()` accessors; re-export from `ast.rs` beside `EntitySpan`.
- C1b — `ReactionSpanAst.constraints: Vec<ConstraintSpan>` field + `constraints()` accessor.
- C1c — `ReactionSpanAst::from_parts(graph, atoms: Vec<EntitySpan<AtomAst>>, bonds:
  Vec<EntitySpan<BondAst>>, constraints: Vec<ConstraintSpan>)` (`pub(crate)`).
- *Checkpoint:* AST compiles.

**C2 — AST: carry constraints through conversion** [C1] (`ast/reaction_span.rs`): **Done**
- C2a — **Done.** `to_reaction_span` collects `ConstraintDelta::Add`/`Remove` from the delta loop, then
  builds the column by multiset subtraction: each lhs constraint matched by a `Remove` → `Removed`
  (consuming one match), the rest → `Unchanged`; each `Add` → `Added`.
- C2b — **Done.** `project()` takes a `constraint_side` selector (`left = Unchanged+Removed`,
  `right = Unchanged+Added`) and builds the molecule via `MoleculeAst::from_parts`. Per the decision,
  the side's constraints are **remapped through the projection's compaction** — an `IdRemapping` built
  from the removed-node/edge lists drives `Constraints::remap`, so refs to a removed atom/bond are
  dropped (consistent with how `project` already compacts bonds; this is removal-compaction, not the
  deferred arbitrary-match remap of doc 134 item 1).
- C2c — **Done.** `to_reaction` appends `Delta::Constraint(ConstraintDelta::Add/Remove)` per
  `Added`/`Removed` (skips `Unchanged`); refs stay in the span frame (lhs preserved via `left()`).
- C2d — **Done.** Three tests: `constraints_unchanged` (molecule → span, both projections carry it),
  `to_reaction_constraints` (`Added` round-trips to the operational reaction), and
  `project_drops_dangling_constraint` (a constraint naming a removed atom is kept on the left, dropped
  on the right by the remap).
- *Checkpoint (tested):* conversion carries constraints. **Done** (8 reaction_span tests, full suite
  3857 pass).

**C3 — DSL module + raw input types** [—] (`dsl/reaction_span.rs`): **Done**
- C3a — **Done.** Created `dsl/reaction_span.rs`; `pub(crate) mod reaction_span;` added to `dsl.rs`
  (before `refs`). Shell imports only what the input types reference (`ConstraintDsl`, `AtomRef`,
  `AtomAst`, `BondAst`, `EntitySpan`) to stay warning-clean; the rest grow with C6–C8.
- C3b — **Done.** `#[derive(Debug)] pub(crate) enum ConstraintSpanInput { Unchanged | Added |
  Removed }(ConstraintDsl)`.
- C3c — **Done.** `#[derive(Debug)] pub(crate) struct SpanInput { atoms: Vec<(Option<String>,
  EntitySpan<AtomAst>)>, bonds: Vec<(Option<String>, [AtomRef; 2], EntitySpan<BondAst>)>, constraints:
  Vec<ConstraintSpanInput> }`. Plan's `AtomRefDsl` corrected to **`AtomRef`** (the real ref type from
  `dsl/refs.rs`, as used for endpoints in `dsl/reaction.rs`). `into_ast` deferred to C8.
- *Checkpoint:* compiles. **Done** (two transient dead-code warnings — `ConstraintSpanInput` /
  `SpanInput` unused until C6 parsing / C8 `into_ast`).

**C4 — `ReactionSpanDsl` boundary shell** [C1, C3]: **Done**
- C4a — **Done.** `pub struct ReactionSpanDsl { ast: ReactionSpanAst, metadata: MoleculeMetadata }`
  (private fields) + doc; `pub use reaction_span::ReactionSpanDsl;` added to `dsl.rs`. No `Default`
  (`ReactionSpanAst` has none).
- C4b — **Done.** Inherent API: `from_parts`, `into_parts`, `ast`, `metadata` (mirrors `ReactionDsl`).
- C4c — **Done.** Derives `Clone, Debug, PartialEq, Eq` (`MoleculeMetadata` supplies all four).
- *Checkpoint:* type exists + exported. **Done** (build clean; `ReactionSpanDsl` warnings cleared by
  the `pub use` — only C3's `SpanInput`/`ConstraintSpanInput` remain, until C6/C8).

**C5 — AST↔DSL conversion** [C4]: **Done**
- C5a — **Done.** `FromAst<ReactionSpanAst>` (`Ctx = MoleculeDefaults`): rebuilds the span via
  `ReactionSpanAst::from_parts`, lowering each `EntitySpan` side with `AtomDsl/BondDsl::from_ast`
  (`cfg.atom` / `cfg.bond`); constraints pass through (as in `MoleculeDsl`); `metadata =
  MoleculeMetadata::default()`. A private `map_span` helper applies the per-side converter across the
  four `EntitySpan` variants (used 4×; no `EntitySpan::map` on the AST type).
- C5b — **Done.** `IntoAst<ReactionSpanAst>`: inverse raise via `AtomDsl/BondDsl(..).into_ast`, same
  `from_parts` rebuild (clones rather than `atoms_mut`, since `ReactionSpanAst` exposes only accessors).
- C5c — **Done.** `test_reaction_span_dsl_from_ast` — `#[case]` table asserting `from_ast → into_ast`
  identity: `modify` (Modified bond + Unchanged atoms + Unchanged constraint) and `add_remove`
  (Unchanged/Removed/Added atoms & bonds + Added constraint).
- *Checkpoint (tested):* programmatic `ReactionSpanAst` ↔ `ReactionSpanDsl`, no EDN. **Done** (full
  suite 3862 pass).

**C6 — Span entry parsers** [C3] (`dsl/reaction_span.rs`): **Done**
- C6a — **Dropped the `SpanOp` enum / separate `classify_span_op` pass** (it was redundant double
  classification — `EntitySpan`/`ConstraintSpanInput` *are* the classification). Replaced by a small
  `verb_wrapper(&Edn) -> Option<(&str, &Edn)>` discriminator (borrows the verb str, no clone; keys only
  on `add`/`modify`/`remove`, so a bare single-key constraint like `{:connected …}` is correctly
  `None`). Each entry parser matches it inline and builds the span variant directly.
- C6b — **Done.** `parse_atom_span_entry` → `(Option<String>, EntitySpan<AtomAst>)`: `split_span_entry`
  splits the optional outer `[<id> <body>]`; `verb_wrapper` then dispatches; values parsed to `AtomAst`
  via `AtomDsl::from_edn(..).0` (`Modify` = `[left right]`).
- C6c — **Done.** `parse_bond_span_entry` → `(Option<String>, [AtomRef; 2], EntitySpan<BondAst>)`:
  `Unchanged`/`Add`/`Remove` reuse `parse_bond_entry`; `Modify` uses `split_bond_frame` (handles
  `[a b [left right]]` and the `{:id :atoms :type [left right]}` map) + `pair`.
- C6d — **Done.** `parse_constraint_span_entry` → `ConstraintSpanInput`: `verb_wrapper`; `c` via
  `ConstraintDsl::from_edn`; `:modify` / other verbs error.
- C6e — **Done.** Three `#[case]` tables (atom 5, bond 6 incl. map + modify-map, constraint 3); 16
  reaction_span tests pass. `ConstraintSpanInput` gained `PartialEq` for the assertion.
- *Decisions:* span atom/bond **values are complete entities, not aliases** (`AtomDsl/BondDsl::from_edn`,
  matching `EntitySpan<AtomAst>`); the generic map helpers `two_atom_refs` / `required_key` /
  `optional_id` / `pair` moved from `molecule.rs` to `edn_utils.rs` (not molecule-specific). `AtomRefDsl`
  in the plan was the nonexistent name for `AtomRef`.
- *Checkpoint (tested):* entry parsers. **Done** (full suite 3876 pass).

**C7 — Top-level parse** [C3, C6] (`dsl/reaction_span.rs`): **Done**
- C7a — **Done.** `parse_span_input` (tree, via `parse_vec` + the C6 entry parsers) and
  `read_span_input` (streaming, via `read_map`/`read_vec`, buffering each section element with
  `read_value_slice` → `read_string` → the tree entry parser). Unknown top-level key errors; a plain
  molecule map (only `:atoms`/`:bonds`/`:constraints`, all entries bare) parses as an all-`Unchanged`
  span. (Per C7a the streaming reader intentionally buffers to a tree rather than duplicating the
  molecule entry grammar — a documented exception to streaming-never-delegates-to-tree.)
- C7b — **Done.** `test_parse_span_input` `#[case]` table (`full` span map, `plain_molecule`) asserts
  both `parse_span_input` and `read_span_input` produce the expected `SpanInput` (tree/streaming
  parity). `SpanInput` gained `PartialEq` for the assertion.
- *Checkpoint (tested):* span map parses. **Done** (18 reaction_span tests; full suite 3878 pass).

**Alias support + metadata decision (amends C3/C6/C7).** Spans support `:atom-aliases` exactly like a
molecule map — the earlier "span values are complete, no aliases" was an implementation gap, not a
design property. So: `SpanInput.atoms` holds `EntitySpan<AtomSpecInput>` (`Bare | Alias`, unresolved)
plus an `atom_aliases: Vec<(String, Box<AtomDsl>)>` field; C6 parses atom values via `parse_atom_entry`
(the molecule entry layer, alias-aware) not `AtomDsl::from_edn`; C7 parses the `:atom-aliases` section
(`parse_atom_aliases`/`read_atom_aliases`). Bonds/constraints unchanged (no aliases). This also settles
the metadata type: **`MoleculeMetadata` is the direct, hygienic correspondence** — the span's union
frame is molecule-shaped and `atom_aliases` is a legitimate span field, so no dedicated type / refactor
is needed. (Bond DSL keywords `:single`/`:double`/`:triple`/`:aromatic` flow through `BondDsl::from_edn`
— the same `:type` leaf parser `parse_bond_entry` uses — and are covered by tests.)

**C8 — Resolution** [C1, C3]: `SpanInput::into_ast → (ReactionSpanAst, MoleculeMetadata)`. Work items in
dependency order (no fresh-id assignment — union-frame; C8a also builds the bijective alias table and
populates `MoleculeMetadata.atom_aliases`, and atom `AtomSpecInput` sides resolve to `AtomAst` —
`Bare → .0`, `Alias → table lookup`, unknown → error):
- C8a — `into_ast` skeleton + union namespace = atom-entry positions ∪ inline `:id`s →
  `MoleculeMetadata`.
- C8b — resolve each bond's `[AtomRefDsl; 2]` against the namespace (unknown ref → error).
- C8c — build `Graph` from atom entries (nodes) + bond entries (edges) — every entry holds a union slot
  regardless of op.
- C8d — resolve `ConstraintSpanInput` refs → `ConstraintSpan` against the namespace.
- C8e — per-side ref consistency: the left projection (`Unchanged` ∪ `Removed` ∪ `Modified.left`) and
  the right projection (`Unchanged` ∪ `Added` ∪ `Modified.right`) must each be internally
  ref-consistent (run the molecule ref/consistency check on each side).
- C8f — assemble `Vec<EntitySpan<…>>` + `Vec<ConstraintSpan>` → `ReactionSpanAst::from_parts`.
- C8g — tests: hand-built `SpanInput` resolves; per-side inconsistency error; unknown-ref error.
- *Checkpoint (tested):* `SpanInput` resolves.

**C9 — `FromEdn` path** [C4, C7, C8]:
- C9a — `FromEdn` (`from_edn` + `from_edn_str`): parse (C7) → `into_ast` (C8) → `from_parts`.
- C9b — `FromStr` (`Err = ParseError`, delegates `from_edn_str`).
- C9c — tests: span EDN string → `ReactionSpanDsl` / `ReactionSpanAst`; a plain molecule map →
  all-`Unchanged`.
- *Checkpoint (tested):* EDN → AST end to end.

**C10 — Render** [C1] (`dsl/reaction_span.rs`; no partial parser — span values are complete):
- C10a — render atom span entries: bare for `Unchanged`, else `{:add v}` / `{:remove v}` /
  `{:modify [l r]}` via the molecule atom renderer.
- C10b — render bond span entries (same).
- C10c — render constraint span entries: `:constraints` bare / `{:add}` / `{:remove}`.
- C10d — `render_span_edn`: assemble `:atoms` / `:bonds` / `:constraints`.
- C10e — tests: `ReactionSpanAst` → EDN per op / entity.
- *Checkpoint (tested):* AST renders to EDN.

**C11 — `ToEdn` path + AST routing** [C4, C9, C10]:
- C11a — `ToEdn` (`render_span_edn`).
- C11b — `Display` (`write!("{}", self.to_edn())`).
- C11c — `ReactionSpanAst` routing: `FromEdn` / `ToEdn` / `FromStr` / `Display` through
  `ReactionSpanDsl` (discard `metadata`).
- C11d — tests: full span DSL → AST → DSL round-trip; homoiconicity; `:modify [left right]`; `:add` /
  `:remove`; `:constraints` add/remove survival through `left()` / `right()`.
- *Checkpoint (tested):* full round-trip.

**C12 — Spec** [C11]:
- C12a — add the span surface beside the delta surface in the reactions section of `umol-dsl-spec.md`.
- *Checkpoint:* spec updated.

**Increment-2 (doc 134):** overlay deltas, overlay span columns, molecule-constraint `apply_at`.
