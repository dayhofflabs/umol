---
description: Apply on ANY task that adds, extends, revises, or reviews DSL serialization/deserialization in a umol-workspace crate — `FromEdn`/`ToEdn` impls, `*Dsl` boundary types, compact-string DSL (`FromStr`/`Display`, winnow parsers), EDN readers/writers (tree or streaming), or the free functions that build/encode AST types. Trigger whenever serde code is created or edited, or whenever the request mentions EDN, the DSL, `FromEdn`/`ToEdn`, `*Dsl` types, `parse_`/`fmt_`/`read_`/`render_` functions, `_to_edn`/`_from_edn`, streaming deserialization, or roundtrip. Covers which types must own their serde via traits vs which may use free functions, the mandatory function-naming scheme, the same-prefix rule for shared helpers, and the ban on single-use helpers. Consult before writing or editing any DSL/EDN ser/de, and re-check naming and structure on every such edit.
---

# umol DSL serialization conventions

**Scope: all DSL/EDN serialization lives in the `umol-ast` crate (its `dsl` module). No serialization anywhere else.** Other crates (umol-graph, umol-io, …) must go through `umol-ast`'s DSL boundary types — they do not define their own `FromEdn`/`ToEdn`/`parse_`/`fmt_`/`read_`/`render_` for umol DSL or EDN. (Foreign *format* parsers — MOL/SDF/SMILES → AST in umol-io — are a separate concern, not this DSL/EDN serde.)

Two surfaces: the **compact string DSL** (`FromStr`/`Display`, winnow) and the **EDN** form (tree `FromEdn`/`ToEdn` + streaming readers). These rules govern how the work is structured and named. When a type can't serialize cleanly, **restructure it** — never declare the problem unsolvable and stop.

## Compact string DSL: combinators vs `parse_` wrappers

The winnow **combinator** `<type>()` (`fn(&mut &str) -> PResult<T>`, used with `.parse_next`) is for **composition only**. Higher-level code that parses a complete string **must not** call a combinator directly — it goes through the **top-level wrapper** `parse_<type>(input: &str) -> Result<T, ParseError>` (which does `<type>.parse(input)`).

- A `parse_<type>` wrapper exists **iff** the type is parsed at top level — i.e. there is (or should be) a `<type>.parse(…)` entry point for it. A combinator used only *inside* other combinators needs **no** wrapper.
- Any `<combinator>.parse(…)` call in non-combinator code is a violation — replace it with the `parse_<type>` wrapper. Do not add a wrapper that nothing calls just to "complete the set."

## Tree and streaming deserialization (EDN)

There is **deliberate** duplication of parsing between the tree-based approach (`FromEdn` over `&Edn`) and the streaming approach (`read_*` over the deserializer) — it is required, not a smell.

- Both **must** be available for types that appear in the **EDN syntax** — a streaming `read_<type>`/`render_<type>` **and** a tree `read_edn_<type>`/`render_edn_<type>`.
- **Only** streaming parsing must be provided for types that are encoded only inside the **entity string DSL**.
- Stream-based parsing **must not** delegate to the tree-based parsers (nor vice versa). They are independent and distinguished by name, not location.

## i. Top-level DSL types own their serde via traits

A top-level `*Dsl` boundary type **must** serialize via `FromEdn`/`ToEdn` (and, where it has a string form, `FromStr`/`Display`). **No exceptions.**

- If a `*Dsl` lacks the data to de/serialize itself, **restructure it to carry that data** — e.g. a stereo constraint needs `StereoKind` for the permutation degree, so the boundary type is `StereoAtomConstraintDsl(StereoKind, _)`. Do **not** fall back to a free helper that takes the missing data as a parameter.
- There must be **no `_to_edn` / `_from_edn` free helpers**. This includes the per-constraint `*Dsl` types.
- `*Dsl` types **should transparently wrap** their AST type (`StructDsl(pub StructAst)`); deviate only where absolutely necessary.

## ii. Constituent types may use free functions

Structs/enums that are used to build the AST types and do **not** have their own `*Dsl` type **may** have free-function ser/de, named:

| surface | deserialize | serialize |
|---|---|---|
| compact string DSL | `parse_<type>` | `fmt_<type>` |
| EDN — streaming (deserializer cursor) | `read_<type>` | `render_<type>` |
| EDN — tree (`&Edn` / `Edn`) | `read_edn_<type>` | `render_edn_<type>` |

`_to_edn` / `_from_edn` are **not** valid names. The tree vs streaming distinction is by **name** (`read_` vs `read_edn_`, `render_` vs `render_edn_`), **never by module** — code gets restructured/moved, so location is not a reliable discriminator.

A constituent free function **may** take extra parameters (e.g. `read_permutation(de, degree)`); needing a parameter is not a reason to skip the convention.

## iii. Shared helpers carry the served prefix

Additional free functions **may** factor out code shared by **multiple** group-(ii) functions. Such a helper **must** carry the **same prefix** as the functions it serves — e.g. `read_member` shared by `read_ligand_symmetry` / `read_topicity` / `read_stereogenicity`; `fmt_*` helpers serve `fmt_*` callers.

## iv. No single-use helpers

A free function used at exactly one call site **must not** exist — inline it at that site. (Macro-generated functions count per generated instance: if the macro emits a helper used once per expansion, inline it into the expansion.)

## Audit checklist

When touching `dsl/*`:
- Grep for `_to_edn` / `_from_edn` names → rename to `render_` / `read_` (or fold into a `*Dsl` trait impl).
- Every `*Dsl` type has `FromEdn` + `ToEdn` (+ `FromStr`/`Display` if it has a string form), self-sufficient — push any missing data up into the type.
- Every free ser/de fn matches `parse_`/`fmt_`/`read_`/`render_<type>`; shared helpers share the prefix.
- No single-use free helper survives; inline it.
- Tree and streaming readers for the same type stay consistent (same shape, same field/key handling) but remain **separate** — streaming must not call the tree path. Both exist for EDN-syntax types; entity-string-only types get streaming only.
