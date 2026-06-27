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
  AST, rewrapped via `from_parts`/`into_parts` — created-entity handles and lhs ids are metadata.
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

A created entity may carry an `<id>` handle that later deltas/bonds reference. (`:electrons` is
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
{:atom            {:modify <atom-ref> <atom-dsl>}}
{:bond            {:modify <bond-ref> <bond-dsl>}}
{:aromatic-system {:modify <aromatic-ref> <aromatic-system-dsl>}}
{:stereo-atom     {:modify <stereo-atom-ref> <stereo-atom-dsl>}}   ; analogously
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
  through the existing complete-entity parser; this wants a sparse/partial parse mode (D2).

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
  lhs frame with fresh handles for `Add`. Both reuse the int-or-keyword ref + `Metadata`
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
   old values. Ref / id / added-handle resolution is **analogous to the existing `MoleculeInput`** —
   the same pattern, not new machinery.
4. **Value string shape** — **resolved**: delta/span values are whole entity DSL strings, so
   per-field values (incl. the inline-only `Element`/`IsotopeMass`/`Spin`) ride inside; no bare
   per-field slots. `ValueAst` literals *inside* a `ConstraintDsl` / molecule constraint follow the
   molecule surface's value-expr-string convention.
5. **`:modify` key name** — **resolved: `:modify`** (`:set` rejected — the op also adds/removes
   per-entity constraints; chosen over `:update` to share one word with the span's
   `EntitySpan::Modified` for a kept-but-changed entity).
6. **Modify RHS partial parse** — **resolved**: the `:modify` RHS is a *partial* entity DSL (only the
   changed parts), parsed sparsely (D2), not via the complete-entity parser.
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
ReactionMetadata};`). The partial-DSL parsers (D2) are free fns in the existing `dsl/atom.rs` /
`dsl/bond.rs`, next to `parse_atom` / `parse_bond`. No other new modules.

**D0 — Prerequisites** (molecule-DSL; no new types):
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
- 0.2 — `:donor` → `:donors`: `dsl/molecule.rs` dative-bond-entry parser + `dsl/dative.rs`, + spec §4.

**D1 — Types** (all in `dsl/reaction.rs`):
- 1.1 — `ReactionDsl` (boundary). Shape: `struct ReactionDsl { ast: ReactionAst, metadata:
  ReactionMetadata }` (private fields). API: `from_parts(ReactionAst, ReactionMetadata) -> Self`,
  `into_parts(self) -> (ReactionAst, ReactionMetadata)`, `ast(&self) -> &ReactionAst`, `metadata(&self)
  -> &ReactionMetadata`. Traits: `FromEdn<'de>`, `ToEdn`, `FromStr`, `Display`, `FromAst<ReactionAst>`,
  `IntoAst<ReactionAst>` (reuse `MoleculeDefaults` as `Ctx` — the reaction adds no fields or
  predicates, so a separate `ReactionDefaults` would be an empty duplicate; no new defaults type).
- 1.2 — `ReactionMetadata`. Shape: `struct ReactionMetadata { lhs: molecule::Metadata, atom_handles:
  BiBTreeMap<String, AtomId>, bond_handles: BiBTreeMap<String, BondId> }` (the lhs molecule metadata +
  created-entity handle ↔ id bindings). Traits: `Clone`, `Debug`, `Default`. API: accessors.
- 1.3 — `ReactionInput` (raw parse target, `pub(crate)`). Shape: `struct ReactionInput { lhs:
  MoleculeInput, deltas: Vec<DeltaInput> }`. API: `into_ast(self) -> Result<(ReactionAst,
  ReactionMetadata), DeError>` (resolution, D5). Traits: `Debug`.
- 1.4 — `DeltaInput` (raw delta, `pub(crate)`). Shape: enum (refs unresolved, `:modify` RHS a
  partial AST):
  ```
  enum DeltaInput {
      AtomAdd { handle: Option<String>, value: AtomAst },
      AtomRemove(AtomRefDsl),
      AtomModify(AtomRefDsl, AtomAst),            // partial: unspecified fields Undetermined
      BondAdd { handle: Option<String>, atoms: [AtomRefDsl; 2], value: BondAst },
      BondRemove(BondRefDsl),
      BondModify(BondRefDsl, BondAst),            // partial
      ConstraintAdd(ConstraintDsl),
      ConstraintRemove(ConstraintDsl),
  }
  ```
  Traits: `Debug`. (Reuses the existing `AtomRefDsl`/`BondRefDsl` and `ConstraintDsl`.)
- 1.5 — `ReactionAst` gets `FromEdn`/`ToEdn`/`FromStr`/`Display` routing through `ReactionDsl`
  (discarding `metadata`), mirroring `MoleculeAst` (in `dsl/reaction.rs`).

**D2 — Partial entity-DSL parsers** (free fns):
- 2.1 — `pub fn read_partial_atom(&str) -> Result<AtomAst, ParseError>` in `dsl/atom.rs`: sparse
  grammar — element optional, every unspecified field `Undetermined`; reuses the existing predicate
  sub-parsers. Calls `partial_atom` winnow parser.
- 2.2 — `pub fn read_partial_bond(&str) -> Result<BondAst, ParseError>` in `dsl/bond.rs`: order
  optional. Calls `partial_bond` winnow parser.
- 2.3. — `fmt_partial_atom` / `fmt_partial_bond` (render only the non-`Undetermined` fields),
  same files (for D6).
- 2.4 - `read_edn_partial_atom` / `read_edn_partial_bond` (used in `FromEdn` impl).
- 2.5 - `render_partial_atom` / `render_partial_bond` (used in `ToEdn` impl).

**D3 — Delta parser** (`dsl/reaction.rs`, free fns):
- 3.1 — `read_delta_input(&mut EdnStreamDeserializer) -> Result<DeltaInput, EdnError>` (streaming) and
  `parse_delta_input(&Edn) -> Result<DeltaInput, DeError>` (tree). Dispatch: entity keyword → op keyword →
  payload.
- 3.2 — payload parsers: `:add` reuses the molecule entry parser for the entity; `:remove` → an
  `AtomRefDsl`/`BondRefDsl` (incl. the structural form); `:modify` → ref + `parse_partial_*` (D2);
  `:constraint` → `ConstraintDsl::from_edn`.

**D4 — Top-level parser** (`dsl/reaction.rs`, free fns):
- 4.1 — `read_reaction_input(&mut …) -> Result<ReactionInput, EdnError>` and
  `parse_reaction_input(&Edn) -> Result<ReactionInput, DeError>`: `:lhs` via `read_molecule_input` /
  `parse_molecule_input`, `:deltas` via D3.
- 4.2 — `ReactionDsl::from_edn` / `from_edn_str` call 4.1 then `ReactionInput::into_ast` (D5).

**D5 — Resolution** (`ReactionInput::into_ast`, `dsl/reaction.rs`):
- 5.1 — resolve `lhs` (`MoleculeInput::into_ast`) → `(MoleculeAst, molecule::MoleculeMetadata)`.
- 5.2 — ref namespace = lhs ids ∪ `:add` handles; resolve each `AtomRefDsl`/`BondRefDsl` against it
  (reuse the molecule ref resolver); error on unknown / non-covering ref.
- 5.3 — `AtomAdd`/`BondAdd` → `Add` with a fresh id (lhs count + order) + register the handle.
- 5.4 — `AtomRemove`/`BondRemove` → recover the entity's `ast` from `lhs` → `Remove`.
- 5.5 — `AtomModify`/`BondModify` → for each field/constraint present in the partial RHS emit a
  `ModifyField` / `ModifyConstraint` whose `old` is the `lhs` entity's current value and whose `new`
  is the RHS value **as written** — a `ValueAst`, possibly a symbolic expression in `lhs`-bound vars
  (`?c+1`), stored unevaluated and resolved at apply, not at the boundary. An `Undetermined`
  constraint in the RHS = removal.
- 5.6 — `ConstraintAdd`/`ConstraintRemove` → `Delta::Constraint(ConstraintDelta::{Add,Remove})`
  (resolve constraint refs against the namespace).

**D6 — Render** (`dsl/reaction.rs`):
- 6.1 — `ReactionDsl::from_ast(&ReactionAst, &MoleculeDefaults) -> Self` and `ToEdn::to_edn`.
- 6.2 — `render_reaction_edn(&ReactionAst, &ReactionMetadata) -> Edn`: `:lhs` via the molecule
  renderer; `:deltas` via 6.3.
- 6.3 — `render_deltas`: group by entity; **coalesce** one entity's `ModifyField` + `ModifyConstraint`
  into a single `:modify <ref> <partial-DSL>` (via `fmt_partial_*`, dropping `old`, `#v*` for a
  removed constraint); emit `:add` / `:remove` / `:constraint`; canonical (post-`canonicalize`) order.

**D7 — Tests** (`#[cfg(test)]` in `dsl/reaction.rs`, `#[rstest]`): round-trip DSL→AST→DSL;
recover-from-`lhs` (a `?`-var `:modify` reading the lhs value); add-then-reference; structural remove;
molecule `:constraint`; partial-RHS parse.

**D8 — Spec** (`umol-ast/spec/umol-dsl-spec.md`)
- 8.1 - add normative reactions section; revise §8.4.
- 8.2 - rename :donor -> :donors, verify if single-int :donor is allowed, if so, remove.
- 8.3 - move aromatic system / multicenter bond :electrons to the string DSL head.

## Implementation plan — span surface (`ReactionSpanDsl` ↔ `ReactionSpanAst`)

Increment-1 entities (atom / bond / molecule-constraint); overlay span columns deferred (doc 134
item 2). Span `:modify` carries complete `[left right]` values (the span is self-contained), so this
reuses the complete molecule entry parsers — no partial parser, unlike the delta surface.

**Module placement.** New module `umol-ast/src/dsl/reaction_span.rs`, registered in `dsl.rs`
(`pub(crate) mod reaction_span;` + `pub use reaction_span::ReactionSpanDsl;`). AST changes land in
`ast/delta.rs` and `ast/reaction_span.rs`.

**S0 — AST: constraint span:**
- 0.1 — `ConstraintSpan` enum in `ast/delta.rs`, beside `EntitySpan`. Shape: `enum ConstraintSpan {
  Unchanged(Constraint), Added(Constraint), Removed(Constraint) }` (no `Modified` — molecule
  constraints are a by-value multiset). API: `left(&self) -> Option<&Constraint>` (`Unchanged`|
  `Removed`), `right(&self) -> Option<&Constraint>` (`Unchanged`|`Added`). Traits: `Clone, Debug,
  PartialEq, Eq`. Re-export from `ast.rs` beside `EntitySpan`.
- 0.2 — `ReactionSpanAst.constraints: Vec<ConstraintSpan>` field + `constraints(&self) ->
  &[ConstraintSpan]` accessor (`ast/reaction_span.rs`).
- 0.3 — `ReactionSpanAst::from_parts(graph: Graph, atoms: Vec<EntitySpan<AtomAst>>, bonds:
  Vec<EntitySpan<BondAst>>, constraints: Vec<ConstraintSpan>) -> Self` (`pub(crate)`) — the span is
  built only by `to_reaction_span` today; the DSL needs a constructor.

**S1 — AST: carry constraints through conversion** (`ast/reaction_span.rs`):
- 1.1 — `to_reaction_span`: replace `Delta::Constraint(_) => {}` — each surviving `lhs` constraint →
  `Unchanged`, each `ConstraintDelta::Add`/`Remove` → `Added`/`Removed`; populate `constraints`.
- 1.2 — `project()` (`left()`/`right()`): set the projected molecule's constraints (`left` =
  `Unchanged`+`Removed`, `right` = `Unchanged`+`Added`) via `MoleculeAst::from_parts` (empty relation
  vecs) instead of `from_atoms_and_bonds`, so molecule→span→molecule keeps constraints.
- 1.3 — `to_reaction`: append `Delta::Constraint(ConstraintDelta::Add/Remove)` for each `Added`/
  `Removed` span (skip `Unchanged`).

**S2 — DSL types** (`dsl/reaction_span.rs`):
- 2.1 — `ReactionSpanDsl` (boundary). Shape: `struct ReactionSpanDsl { ast: ReactionSpanAst, metadata:
  MoleculeMetadata }` (private fields; the span is union-frame molecule-shaped, so reuse
  `MoleculeMetadata` for keyword↔index bindings — no new metadata type). API: `from_parts`/
  `into_parts`/`ast`/`metadata`. Traits: `FromEdn<'de>`, `ToEdn`, `FromStr`, `Display`,
  `FromAst<ReactionSpanAst>`, `IntoAst<ReactionSpanAst>` (reuse `MoleculeDefaults` as `Ctx`, as
  `ReactionDsl` — no new defaults type).
- 2.2 — `SpanInput` (raw parse target, `pub(crate)`). Shape:
  ```
  struct SpanInput {
      atoms:       Vec<(Option<String>, EntitySpan<AtomAst>)>,                  // handle, op + complete value
      bonds:       Vec<(Option<String>, [AtomRefDsl; 2], EntitySpan<BondAst>)>, // endpoints once (frame-invariant)
      constraints: Vec<ConstraintSpanInput>,
  }
  ```
  API: `into_ast(self) -> Result<(ReactionSpanAst, MoleculeMetadata), DeError>` (S5). Traits: `Debug`.
- 2.3 — `ConstraintSpanInput` (`pub(crate)`). Shape: `enum ConstraintSpanInput {
  Unchanged(ConstraintDsl), Added(ConstraintDsl), Removed(ConstraintDsl) }` (refs unresolved). Traits:
  `Debug`.

**S3 — Entry parsers** (`dsl/reaction_span.rs`, free fns). A generic `read_span_entry<T>` does **not**
fit: a span entry splits into **slot-level** data (shared by both sides — the atom `:id` handle; the
bond handle + `[<ref> <ref>]` endpoints) and a **per-side** value (the `EntitySpan`), and the two
entities place that slot data differently. So a shared op-classifier handles the `bare | {:<op> …}`
wrapper and each per-entity parser owns its slot/value split:
- 3.1 — `classify_span_op(&Edn) -> (SpanOp, &Edn)`, `enum SpanOp { Unchanged, Add, Remove, Modify }`:
  a bare value ⇒ `Unchanged` (payload = the value); a single-key map `{:add|:remove|:modify <p>}` ⇒
  that op (payload = `<p>`). The only place the op wrapper is recognized; `SpanOp` is a parse-internal
  helper (no AST counterpart).
- 3.2 — `parse_atom_span_entry(&Edn) -> Result<(Option<String>, EntitySpan<AtomAst>), DeError>`: peel
  the optional **outer** `[<id> …]` handle (slot-level — one handle, not one per side), classify the
  op (3.1), then parse the value(s) with the molecule atom-entry parser (the one reused by delta
  D3.2): one complete atom for `Unchanged`/`Add`/`Remove`, the `[left right]` pair for `Modify`.
- 3.3 — `parse_bond_span_entry(&Edn) -> Result<(Option<String>, [AtomRefDsl; 2], EntitySpan<BondAst>), DeError>`:
  a bond carries its endpoints (and optional handle) **once, inside** the payload — `[<ref> <ref>
  <bond-dsl>]` for `Unchanged`/`Add`/`Remove`, `[<ref> <ref> [left right]]` for `Modify` (or the
  `{[:id …] :atoms […] :type …}` map form, which also supplies the handle). Classify the op (3.1),
  split off the shared `[<ref> <ref>]` endpoints + handle, wrap the per-side bond value(s) in
  `EntitySpan<BondAst>`.
- 3.4 — `parse_constraint_span_entry(&Edn) -> Result<ConstraintSpanInput, DeError>`: op wrapper
  `bare | {:add c} | {:remove c}` only (no `:modify`, no slot data); `c` via `ConstraintDsl::from_edn`.

**S4 — Top-level parser** (`dsl/reaction_span.rs`, free fns):
- 4.1 — `read_span_input(&mut EdnStreamDeserializer) -> Result<SpanInput, EdnError>` /
  `parse_span_input(&Edn) -> Result<SpanInput, DeError>`: `:atoms`/`:bonds`/`:constraints` via S3.
  Since S3's classifier needs the whole entry, the streaming reader buffers each section element to an
  `Edn` before dispatching to the tree-form entry parser (3.2/3.3/3.4).
- 4.2 — `ReactionSpanDsl::from_edn` / `from_edn_str` call 4.1 then `SpanInput::into_ast` (S5). A plain
  molecule map (all entries bare) parses as an all-`Unchanged` span (homoiconicity).

**S5 — Resolution** (`SpanInput::into_ast`, `dsl/reaction_span.rs`):
- 5.1 — union namespace = atom-entry positions ∪ inline `:id` handles → `MoleculeMetadata`.
- 5.2 — resolve each bond's `[AtomRefDsl; 2]` against the namespace (reuse the molecule ref resolver);
  error on unknown ref.
- 5.3 — build `Graph` from all atom entries (nodes) and all bond entries (edges) — every entry holds a
  union slot regardless of op.
- 5.4 — resolve `ConstraintSpanInput` refs → `ConstraintSpan` (reuse the molecule constraint resolver).
- 5.5 — validate per-side ref consistency: the left projection (`Unchanged` ∪ `Removed` ∪
  `Modified.left`) and the right projection (`Unchanged` ∪ `Added` ∪ `Modified.right`) must each be
  internally ref-consistent, like a standalone molecule (run the molecule ref/consistency check on
  each side).
- 5.6 — assemble `Vec<EntitySpan<…>>` + `Vec<ConstraintSpan>` → `ReactionSpanAst::from_parts`.

**S6 — Render** (`dsl/reaction_span.rs`):
- 6.1 — `ReactionSpanDsl::from_ast(&ReactionSpanAst, …) -> Self` + `ToEdn::to_edn`.
- 6.2 — `render_span_edn(&ReactionSpanAst, &MoleculeMetadata) -> Edn`: each entry bare for `Unchanged`,
  else `{:add v}` / `{:remove v}` / `{:modify [l r]}`; reuse the molecule entry renderers for the
  complete values; `:constraints` bare / `{:add}` / `{:remove}`.

**S7 — `ReactionSpanAst` ser/de routing** (`dsl/reaction_span.rs`): `FromEdn`/`ToEdn`/`FromStr`/
`Display` route through `ReactionSpanDsl` (discarding `metadata`), mirroring `MoleculeAst`.

**S8 — Tests** (`#[cfg(test)]`, `#[rstest]`): span DSL→AST→DSL round-trip; a plain molecule map parses
as an all-`Unchanged` span (homoiconicity); `:modify [left right]`; `:add` / `:remove`; `:constraints`
add/remove and their survival through `left()` / `right()`.

**S9 — Spec** (`umol-ast/spec/umol-dsl-spec.md`): add the span surface beside the delta surface in the
reactions section.

**Increment-2 (doc 134):** overlay deltas, overlay span columns, molecule-constraint `apply_at`.
