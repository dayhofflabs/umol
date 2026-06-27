# 133 — Reaction EDN surface syntax (design)

Design round for an EDN surface syntax for reactions. We have two reaction AST forms and a
lossless bidirectional conversion between them (doc 131/132): `ReactionAst` (operational —
`lhs` molecule + `Deltas`) and `ReactionSpanAst` (declarative — the superimposed `L ∪_K R`
graph, one `LeftRightState` per atom/bond). This doc designs the surface for **both**, even
though implementation may land one first; the conversion means an author can write either and
render the other.

The molecule parts are already expressible (`MoleculeDsl` is the boundary type, mirrored from the
`MoleculeAst` surface in `umol-dsl-spec.md`). So the weight here is the two genuinely new
encodings: **`Delta`** (work item 1) and **`LeftRightState`** (work item 2). Everything else is
reuse.

**Spec status.** The reaction notation currently in `umol-dsl-spec.md` (the §8.4 `:lhs`/`:rhs`
two-molecule sketch) predates the `ReactionAst` / `ReactionSpanAst` model and is explicitly
non-normative on this topic. This design supersedes it; the spec follows this doc, not the
reverse, and will be revised to whatever is decided here.

## Scope

In scope:
- `ReactionDsl` ↔ `ReactionAst`; `ReactionSpanDsl` ↔ `ReactionSpanAst`.
- The `Delta` and `LeftRightState` encodings, plus the `*FieldChange` leaf encoding.
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
                   | SetField{ id, change: AtomFieldChange }
                   | SetConstraint{ id, old: Option<AtomConstraint>, new: Option<AtomConstraint> }
BondDelta          = Add{ id, atoms: [AtomId;2], ast } | Remove{ id, atoms, ast }
                   | SetField{ id, change: BondFieldChange } | SetConstraint{ id, old, new }
ConstraintDelta    = Add(Constraint) | Remove(Constraint)        // molecule-level
AtomFieldChange    = Element|IsotopeMass|Charge|ImplicitHydrogens|LonePairs|Spin  { old, new }
BondFieldChange    = Order|Charge|Spin  { old, new }

ReactionSpanAst    = { graph, atoms: Vec<LeftRightState<AtomAst>>, bonds: Vec<LeftRightState<BondAst>> }
LeftRightState<T>  = Unchanged(T) | Modified{ left: T, right: T } | Added(T) | Removed(T)
```

Two framing facts that drive the encodings:

- **Delta refs are in the `lhs` frame.** `SetField`/`SetConstraint`/`Remove` name an existing
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
    the delta surface. (Value literals *inside* a `ConstraintDsl` / molecule constraint still follow
    the molecule surface's convention — the separate cross-surface reconciliation.)
  - constraint slots are two distinct leaves by slot (Work item 1): the `ConstraintDelta` payload
    is a **`ConstraintDsl`** (molecule-level — `{:connected {…}}`, `{:and […]}`, or an entity leaf
    `{:atom [<ref> {:valence …}]}`); a `SetConstraint` `old`/`new` is a **per-entity constraint
    single-key map** (`AtomConstraintDsl` `{:valence …}`, `BondConstraintDsl` `:aromatic`).

  So a delta/span value is always a whole entity DSL string (`:single` is its `BondDsl` keyword
  shorthand); per-field values live inside that string, not as separate EDN slots.
- **Constituent free fns.** `Delta`, `LeftRightState`, and `*FieldChange` are constituents of the
  two boundary types, so they get `read_`/`render_` (streaming) and `read_edn_`/`render_edn_`
  (tree) free functions, not `FromEdn`/`ToEdn` — unless we decide they are independently useful
  boundary types (see decision 6).
- **Ctx / defaults.** A `ReactionDefaults` aggregating `MoleculeDefaults` for `lhs` and the atom/
  bond defaults for delta/span values; `zeroed()` is the omission threshold for rendering (drop
  defaulted fields), as in molecules. Open whether `ReactionDefaults` adds anything beyond the
  reused `MoleculeDefaults` — likely not, decision deferred to implementation.

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
collapses the earlier `:set-field` / `:set-constraint` pair into a single `:modify`. The `:modify`
key name is **TBD** (`:set` / `:modify` / `:update`).

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

The surface deliberately does **not** mirror the AST's per-field `SetField{old,new}` /
`SetConstraint{old,new}` 1:1 — that creates more problems than it solves (you would write both old
and new and keep them in agreement with the lhs, e.g. `"C#c?c"` in lhs, `"C#c?c"` as old,
`"C#c(?c+1)"` as new: much typing, easy to desync). The modify payload is just the **RHS** — the
entity DSL carrying only the parts being set. The old values come from the `lhs`; deltas are **not**
independently resolvable, and need not be. This mirrors the two-stage molecule boundary (raw
`MoleculeInput` → resolved `MoleculeDsl`, ids/aliases filled at the boundary): a reaction gets a raw
input form resolved against `lhs` into `ReactionAst`, with per-field old values filled from `lhs`.

Semantics — uniform for fields and constraints (they look identical in the entity DSL, so one
coherent representation throughout, no field/constraint special-casing):
- a **field** present in the RHS replaces; absent = unchanged. `"C"` sets element; `"#c0"` sets
  charge 0; `"#c(?c+1)"` sets charge to `?c + 1` (`?c` bound in the lhs).
- a per-entity **constraint** with a concrete value (`"#v4"`) is set/added; the old need not exist —
  consistent with set-constraint semantics.
- **constraint removal** exploits that an `Undetermined` constraint is vacuous: `"#v*"` parses as
  *remove the valence constraint*.
- the RHS is a **partial** entity DSL — only the changed parts, not a complete entity (`"#c0"` has
  no leading element). It need not satisfy the complete-entity grammar, so it must **not** be forced
  through the existing complete-entity parser; this wants a sparse/partial parse mode. **Open.**

Running example (the SN2):

```clojure
:deltas [{:bond {:add    [:c :nu :single]}}    ; form C–Nu
         {:bond {:remove [:c :br]}}            ; break C–Br (by endpoints)
         {:atom {:modify :br "#c-1"}}          ; Br → bromide (new charge only; old from lhs)
         {:atom {:modify :nu "#c0"}}]          ; Nu → neutral
```

This **dissolves the inline-only serialization problem**: a modify (or add) value is a whole entity
DSL string, so element / isotope (`#i`) / spin (`#u#s`) ride inside it — no standalone leaf
serialization for `ElementAst` / `IsotopeMassAst` / `SpinStateAst` is needed.

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

## Work item 2 — `LeftRightState` encoding  **[decision]**

`ReactionSpanDsl` is the superimposed graph and reads like a molecule whose entries may carry a
transition. An all-`Unchanged` span **is** a molecule map (homoiconic — the same keys `:atoms`/
`:bonds`).

```clojure
{:atoms ["C#h3"
         {:modified ["Br" "Br#c-1"]}
         {:modified ["O#h1#c-1" "O#h1"]}]
 :bonds [{:removed [0 1 :single]}
         {:added   [0 2 :single]}]}
```

The contested axis is **how each of the four states is tagged** while keeping `Unchanged` (the
common case) as cheap as a plain molecule entry.

### Atom-state options

| state | Option I — word-tagged maps | Option II — sigil keywords |
|-------|------------------------------|-----------------------------|
| Unchanged | `"C"` (bare atom-string) | `"C"` |
| Modified | `{:modified ["C#h3" "C#h0"]}` | `{:>> ["C#h3" "C#h0"]}` |
| Added | `{:added "N"}` | `{:+ "N"}` |
| Removed | `{:removed "O"}` | `{:- "O"}` |

Option I is explicit and reads well; Option II is terser and visually CGR-like (`+`/`-`/`>>` are
legal EDN keyword characters) but cryptic. Either way, ids attach as today's inline-id 2-vector,
nesting the state: `[:br {:modified ["Br" "Br#c-1"]}]`, `[:oh {:added "O#h1"}]`.

A third option, **III — positional pair = Modified**: `"C"` = Unchanged, `["C#h3" "C#h0"]` =
Modified (a 2-vector of strings), with `Added`/`Removed` still tagged. Cheapest Modified, but the
2-vector collides with the inline-id form `[:id "C"]` (keyword-first vs string-first
disambiguates, yet mixing ids with Modified gets awkward). Mentioned for completeness; the
collision makes it the weakest.

### Bond-state options

Bonds always carry endpoints, so the state wraps a bond entry (`[<ref> <ref> <bond-spec>]`):

| state | wrapped form | type-on-the-entry form |
|-------|--------------|------------------------|
| Unchanged | `[0 1 :single]` | `[0 1 :single]` |
| Modified | `{:modified [0 1 [:single :double]]}` | `[0 1 {:modified [:single :double]}]` |
| Added | `{:added [0 1 :double]}` | `[0 1 {:added :double}]` |
| Removed | `{:removed [0 1 :single]}` | `[0 1 {:removed :single}]` |

Left column wraps the whole entry (uniform with the atom-state tagging); right column keeps the
endpoints bare and tags only the `:type` slot (endpoints are frame-invariant, only the bond value
changes — except for Added/Removed where the bond itself is present on one side only, which the
right column expresses as a one-sided type). Decision 3 picks the atom and bond encodings together
for consistency.

### Modified value form (decision 4)

`Modified` carries two full ASTs (`left`, `right`). The pair `["C#h3" "C#h0"]` is lossless and
needs no new grammar. Alternatives — an in-string arrow (`"C#h3>>C#h0"`, rejected: pollutes the
atom subgrammar) or a delta-style field diff (`{:charge ["0" "-1"]}`, rejected here: that is the
*delta* surface's job, and the span deliberately carries states not ops). The pair is the natural
span form; the only question is whether `Modified` ever wants to elide the unchanged half (it
cannot — both sides are full values by construction).

## Cross-cutting

- **Frames.** Span = union frame (lhs slots preserved for `Removed`, created appended); deltas =
  lhs frame with fresh handles for `Add`. Both reuse the int-or-keyword ref + `Metadata`
  resolution. ids are optional sugar over positional indices in both surfaces, as in molecules.
- **Canonical rendering.** Drop defaulted fields via the `zeroed()` Ctx; omit empty sections
  (`:deltas []`, an all-`Unchanged` span renders as a plain molecule); `ConstraintsDsl` already
  drops vacuous constraints. Render order follows the AST's canonical (post-`canonicalize`) delta
  order so DSL→AST→DSL is stable.
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

## Open decisions

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
   per-field slots. Residual: `ValueAst` literals *inside* a `ConstraintDsl` / molecule constraint
   follow the molecule surface's convention (the separate cross-surface reconciliation).
5. **`:modify` key name** — **resolved: `:modify`** (`:set` rejected — the op also adds/removes
   per-entity constraints; chosen over `:update` to share one word with the span's
   `LeftRightState::Modified` for a kept-but-changed entity).
6. **Modify RHS partial parse** — the RHS is a *partial* entity DSL (only changed parts); needs a
   sparse parse mode, not the existing complete-entity parser. **Open.**
7. **Electrons vector DSL** — **resolved**. Vector syntax `[<nat>(,<nat>)*]` (square brackets, ≥1,
   comma-separated, whitespace ignored; undetermined `*`), `[…]` paralleling `{…}`/`(…)`; position
   `i` matches atom `i` of `:atoms`. The vector is the **untagged head** of the aromatic-system /
   multicenter DSL (like element is the atom-string head); no `#`-tag, the scalar total stays the
   separate `#e<n>`. Ripples to the molecule DSL + spec.
8. **Remove cross-check** — **resolved**: a `:remove` ref must resolve to an existing entity, but
   that is the same ref-existence validation `MoleculeInput` already performs for atom-refs in bonds
   — reuse it, not new machinery.
9. **`LeftRightState` encoding** (work item 2) — atom-state (I word-tagged / II sigil / III
   positional) and bond-state (whole-entry wrap vs type-slot tag), chosen together. **Open.**
10. **`Modified` value form** (work item 2) — the `[left right]` string pair vs any alternative.
    **Open.**
11. **Boundary vs constituent** — `Delta` is now a constituent of `ReactionDsl` (recover-from-`lhs`
    rules out a standalone `DeltaDsl`); whether the span / `LeftRightState` is its own boundary type
    is **open** (the span is self-contained, so it could be). 
12. **Self-describing tag** — **resolved: none** for now (caller picks the boundary type, as with
    molecules); revisit if mixed/heterogeneous streams arise.
13. **Generalization debt** — `LeftRightState<T>` and the per-entity field-change split are
    atom/bond-shaped and will not survive the overlays; revisit when overlays land. **Deferred.**
14. **Dative bonds** — **resolved** (the asymmetric case: multi-donor + single acceptor, not an
    endpoint pair). `:add` = `{[:id <id>] :donors [<atom-ref>*] :acceptor <atom-ref> :type
    <dative-bond-dsl>}`; `:remove` ref = `<index> | <id> | {:donors [<atom-ref>*] :acceptor
    <atom-ref>}`. The `:donor` → `:donors` rename is a molecule-wide change, recorded separately.
15. **Surface vs AST scope** — the encoding spans all eight entities, but the AST `Delta` is
    `Atom` / `Bond` / `Constraint` only today; the overlay-entity deltas are forward-looking (need
    increment-2 AST support). **Note.**

## After decisions (implementation sketch)

`umol-ast/src/dsl/reaction.rs` (and `reaction_span.rs`): the two `*Dsl` boundary types and a raw
`ReactionInput` (resolved against `lhs`, mirroring `MoleculeInput`), with their
`FromEdn`/`ToEdn`/`FromStr`/`Display`/`FromAst`/`IntoAst`; `read_/render_` + `read_edn_/render_edn_`
free fns for `Delta`, `LeftRightState`, `*FieldChange`; `ReactionDefaults` in `dsl/config.rs` if
needed; round-trip tests (DSL→AST→DSL and the cross-form DSL→AST→span-AST→span-DSL); and the spec
updated once the encoding is fixed — a normative reactions section added and the stale §8.4 sketch
revised to match.
