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

## Work item 2 — `LeftRightState` encoding  **[open — leading proposal, not finalized]**

`ReactionSpanDsl` is the superimposed graph; it reads like a molecule whose entries may carry a
transition, and an all-`Unchanged` span **is** a plain molecule map (homoiconic — the same `:atoms`/
`:bonds` keys).

**Leading proposal (open — not finalized).** Per-element state is tagged with **participle
keywords** — `:added` / `:modified` / `:removed` — and a bare entry is `Unchanged`. This mirrors the `LeftRightState` variants and keeps
the established **verb (operations) vs participle (states)** split: the delta surface uses `:add` /
`:remove` / `:modify` (operations over `lhs`); the span uses `:added` / `:modified` / `:removed`
(per-element states). `Unchanged` stays as cheap as a plain molecule entry.

Atoms — a bare atom DSL is `Unchanged`; the tags wrap the value(s):

```clojure
{:atoms ["C#h3"                                  ; Unchanged
         {:modified ["Br" "Br#c-1"]}             ; Modified — [left right], two full atom DSLs
         {:added "O#h1"}                         ; Added   — right-only value
         {:removed "O"}                          ; Removed — left-only value
         [:nu {:modified ["O#h1#c-1" "O#h1"]}]]  ; with an id: [<id> <state>]
 …}
```

Bonds — endpoints are frame-invariant, so the state wraps the bond entry; `Modified` keeps the
endpoints once and pairs the bond value:

```clojure
 :bonds [[0 1 :single]                           ; Unchanged
         {:modified [0 2 [:single :double]]}     ; Modified — [<ref> <ref> [left right]]
         {:added   [0 2 :single]}                ; Added
         {:removed [0 1 :single]}]               ; Removed
```

`Modified` carries the two full values as a `[left right]` pair — lossless, no new grammar.
**Rejected:** sigil tags (`:+`/`:-`/`:>>`, cryptic); a positional `[a b]` for atom `Modified`
(collides with the inline-id `[:id "C"]` form); an in-string arrow (pollutes the atom subgrammar);
a delta-style field diff (that's the delta surface's job — the span carries states, not ops).

**Generalization debt.** `LeftRightState<T>` and the span's per-entity `atoms`/`bonds` vecs are
atom/bond-shaped and will not survive the six overlay entities; generalizing the span container is
tracked in doc 134 (gap 2). The encoding above is the localized-topology (atom/bond) surface.

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

## Decisions

The **delta-surface** decisions are resolved (work item 1). The **span** encoding (work item 2) is
still **open** — a leading proposal is recorded for context, not finalized. Other remaining items
are implementation gaps, tracked in doc 134.

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
9. **`LeftRightState` encoding** (work item 2) — **open**. Leading proposal: participle tags
   `:added`/`:modified`/`:removed`, bare = `Unchanged`; bond states wrap the entry; `Modified` =
   `[left right]`. Recorded for context, not finalized.
10. **`Modified` value form** — **open** (part of work item 2). Leading proposal: the `[left right]`
    full-value pair.
11. **Boundary vs constituent** — `Delta` is a constituent of `ReactionDsl` (recover-from-`lhs`
    rules out a standalone `DeltaDsl`); the span is self-contained, so `ReactionSpanDsl` is its own
    boundary type. **Resolved.**
12. **Self-describing tag** — **resolved: none** for now (caller picks the boundary type, as with
    molecules); revisit if mixed/heterogeneous streams arise.
13. **Dative bonds** — **resolved** (the asymmetric case: multi-donor + single acceptor, not an
    endpoint pair). `:add` = `{[:id <id>] :donors [<atom-ref>*] :acceptor <atom-ref> :type
    <dative-bond-dsl>}`; `:remove` ref = `<index> | <id> | {:donors [<atom-ref>*] :acceptor
    <atom-ref>}`. The `:donor` → `:donors` rename is a molecule-wide change, recorded separately.

**Deferred to doc 134 (implementation gaps, not surface-design questions):**
- **Generalization debt** — `LeftRightState<T>` and the per-entity span vecs are atom/bond-shaped
  and won't survive the six overlay entities (134 gap 2).
- **Surface vs AST scope** — the surface spans all eight entities, but the AST `Delta` is
  `Atom`/`Bond`/`Constraint` only; overlay-entity deltas need increment-2 AST support (134 gap 2).
- **Molecule-constraint reaction application** — `apply_at` currently drops `Delta::Constraint`;
  closing it needs a match-based constraint ref-remap (134 gap 1).

## Implementation plan — delta surface (`ReactionDsl` ↔ `ReactionAst`)

Increment-1 entities (atom / bond / molecule-constraint); overlay-entity deltas deferred
(doc 134 gap 2). Types and signatures below are approximate (shapes + API, no bodies).

**Module placement.** New module `umol-ast/src/dsl/reaction.rs`, registered in
`umol-ast/src/dsl.rs` (`pub(crate) mod reaction;` + `pub use reaction::{ReactionDsl,
ReactionMetadata};`). The partial-DSL parsers (W2) are free fns in the existing `dsl/atom.rs` /
`dsl/bond.rs`, next to `parse_atom` / `parse_bond`. No other new modules.

**W0 — Prerequisites** (molecule-DSL; no new types):
- 0.1 — `:electrons` → head of the aromatic-system / multicenter-bond DSL string, parsed as a
  `[<nat>(,<nat>)*]` vector: `dsl/aromatic.rs`, `dsl/multicenter.rs` (+ the entry parsers in
  `dsl/molecule.rs` that currently read the `:electrons` map key), + spec.
- 0.2 — `:donor` → `:donors`: `dsl/molecule.rs` dative-bond-entry parser + `dsl/dative.rs`, + spec §4.

**W1 — Types** (all in `dsl/reaction.rs`):
- 1.1 — `ReactionDsl` (boundary). Shape: `struct ReactionDsl { ast: ReactionAst, metadata:
  ReactionMetadata }` (private fields). API: `from_parts(ReactionAst, ReactionMetadata) -> Self`,
  `into_parts(self) -> (ReactionAst, ReactionMetadata)`, `ast(&self) -> &ReactionAst`, `metadata(&self)
  -> &ReactionMetadata`. Traits: `FromEdn<'de>`, `ToEdn`, `FromStr`, `Display`, `FromAst<ReactionAst>`,
  `IntoAst<ReactionAst>` (reusing the lhs `MoleculeDefaults` as `Ctx` — no new defaults type).
- 1.2 — `ReactionMetadata`. Shape: `struct ReactionMetadata { lhs: molecule::Metadata, atom_handles:
  BiBTreeMap<String, AtomId>, bond_handles: BiBTreeMap<String, BondId> }` (the lhs molecule metadata +
  created-entity handle ↔ id bindings). Traits: `Clone`, `Debug`, `Default`. API: accessors.
- 1.3 — `ReactionInput` (raw parse target, `pub(crate)`). Shape: `struct ReactionInput { lhs:
  MoleculeInput, deltas: Vec<DeltaInput> }`. API: `into_ast(self) -> Result<(ReactionAst,
  ReactionMetadata), DeError>` (resolution, W5). Traits: `Debug`.
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

**W2 — Partial entity-DSL parsers** (free fns):
- 2.1 — `pub fn parse_partial_atom(&str) -> Result<AtomAst, ParseError>` in `dsl/atom.rs`: sparse
  grammar — element optional, every unspecified field `Undetermined`; reuses the existing predicate
  sub-parsers.
- 2.2 — `pub fn parse_partial_bond(&str) -> Result<BondAst, ParseError>` in `dsl/bond.rs`: order
  optional.
- 2.3 — matching `fmt_partial_atom` / `fmt_partial_bond` (render only the non-`Undetermined` fields),
  same files (for W6).

**W3 — Delta parser** (`dsl/reaction.rs`, free fns):
- 3.1 — `read_delta(&mut EdnStreamDeserializer) -> Result<DeltaInput, EdnError>` (streaming) and
  `parse_delta(&Edn) -> Result<DeltaInput, DeError>` (tree). Dispatch: entity keyword → op keyword →
  payload.
- 3.2 — payload parsers: `:add` reuses the molecule entry parser for the entity; `:remove` → an
  `AtomRefDsl`/`BondRefDsl` (incl. the structural form); `:modify` → ref + `parse_partial_*` (W2);
  `:constraint` → `ConstraintDsl::from_edn`.

**W4 — Top-level parser** (`dsl/reaction.rs`, free fns):
- 4.1 — `read_reaction_input(&mut …) -> Result<ReactionInput, EdnError>` and
  `parse_reaction_input(&Edn) -> Result<ReactionInput, DeError>`: `:lhs` via `read_molecule_input` /
  `parse_molecule_input`, `:deltas` via W3.
- 4.2 — `ReactionDsl::from_edn` / `from_edn_str` call 4.1 then `ReactionInput::into_ast` (W5).

**W5 — Resolution** (`ReactionInput::into_ast`, `dsl/reaction.rs`):
- 5.1 — resolve `lhs` (`MoleculeInput::into_ast`) → `(MoleculeAst, molecule::Metadata)`.
- 5.2 — ref namespace = lhs ids ∪ `:add` handles; resolve each `AtomRefDsl`/`BondRefDsl` against it
  (reuse the molecule ref resolver); error on unknown / non-covering ref.
- 5.3 — `AtomAdd`/`BondAdd` → `Add` with a fresh id (lhs count + order) + register the handle.
- 5.4 — `AtomRemove`/`BondRemove` → recover the entity's `ast` from `lhs` → `Remove`.
- 5.5 — `AtomModify`/`BondModify` → diff the partial RHS against the `lhs` entity → `ModifyField` /
  `ModifyConstraint` deltas (old from `lhs`; an `Undetermined` constraint = removal).
- 5.6 — `ConstraintAdd`/`ConstraintRemove` → `Delta::Constraint(ConstraintDelta::{Add,Remove})`
  (resolve constraint refs against the namespace).

**W6 — Render** (`dsl/reaction.rs`):
- 6.1 — `ReactionDsl::from_ast(&ReactionAst, &MoleculeDefaults) -> Self` and `ToEdn::to_edn`.
- 6.2 — `render_reaction_edn(&ReactionAst, &ReactionMetadata) -> Edn`: `:lhs` via the molecule
  renderer; `:deltas` via 6.3.
- 6.3 — `render_deltas`: group by entity; **coalesce** one entity's `ModifyField` + `ModifyConstraint`
  into a single `:modify <ref> <partial-DSL>` (via `fmt_partial_*`, dropping `old`, `#v*` for a
  removed constraint); emit `:add` / `:remove` / `:constraint`; canonical (post-`canonicalize`) order.

**W7 — Tests** (`#[cfg(test)]` in `dsl/reaction.rs`, `#[rstest]`): round-trip DSL→AST→DSL;
recover-from-`lhs` (a `?`-var `:modify` reading the lhs value); add-then-reference; structural remove;
molecule `:constraint`; partial-RHS parse.

**W8 — Spec** (`umol-ast/spec/umol-dsl-spec.md`): a normative reactions section; revise §8.4.

**Span surface** (`ReactionSpanDsl` ↔ `ReactionSpanAst`): no plan yet — work item 2's encoding is
still open. **Increment-2 gaps** (overlay deltas, span generalization, molecule-constraint apply):
doc 134.
