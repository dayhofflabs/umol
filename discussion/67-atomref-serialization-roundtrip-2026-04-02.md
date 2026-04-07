# AtomRef serialization and MoleculeAst ↔ EDN roundtrip fidelity

## Context

Removing the Clojure dialect eliminated digit-starting keywords (`:0`, `:1`).
Bond endpoints for indexed atoms changed from `:0 :1` (keywords) to `0 1` (integers).
This introduced `AtomRef` as a newtype for bond endpoint references, and exposed a
deeper problem: faithful roundtripping of the molecule DSL through EDN.

## The problem

EDN has three relevant scalar types: keywords (`:foo`), strings (`"foo"`), integers (`0`).
Serde's data model has two: strings and integers. Keywords do not exist in serde.

When the EDN deserializer encounters `:C` (keyword), it presents `visit_str("C")` to
serde. The keyword type is erased. On serialization, `serialize_str("C")` produces
`"C"` (string), not `:C` (keyword). The roundtrip breaks:

```
Input:   [:C :O :single]      (keywords)
Deser:   AtomRef("C"), AtomRef("O")
Ser:     ["C" "O" :single]    (strings — wrong)
```

The integer case roundtrips correctly because serde distinguishes integers from strings:

```
Input:   [0 1 :single]        (integers)
Deser:   AtomRef("0"), AtomRef("1")
Ser:     [0 1 :single]        (integers — correct, via parse::<usize>() heuristic)
```

## The `parse::<usize>()` heuristic

The initial AtomRef implementation stored a `String` internally and used a heuristic
on serialization: if the string parses as `usize`, emit an integer; otherwise emit a
string. This recovers the integer case but not the keyword case. It's a workaround
for tag erasure, not a principled encoding.

## Serialization strategies considered

### Strategy A: Via Edn tree

```
Write:  MoleculeAst → to_edn() → Edn tree → Display/Formatter → EDN text
Read:   EDN text → parser → Edn tree → EdnDeserializer → MoleculeAst
```

Full type fidelity. `to_edn()` produces `Edn::Keyword("C")`, `Edn::Int(0)`.
Allocates an intermediate tree on the write path.
The read path already uses the Edn tree (needed for alias resolution and pre-validation).
Correct by construction.

### Strategy B: Streaming serde with keyword sentinel

```
Write:  MoleculeAst → Serialize → EdnSerializer → EDN text
Read:   EDN text → EdnStreamDeserializer → Deserialize → MoleculeAst
```

A wrapper type (`EdnKeyword`) serializes via `serialize_newtype_struct("__edn_keyword", &s)`.
The `EdnSerializer` recognizes the sentinel name and emits `:keyword` syntax.
Non-EDN serializers see a transparent string. Same pattern as `serde_json::RawValue`.
No intermediate allocation. Private protocol between `umol-edn` types and `EdnSerializer`.

### Strategy C: Direct parser/formatter (bypass Edn and serde)

```
Write:  MoleculeAst → custom write_edn() → EDN text
Read:   EDN text → custom parse_molecule_edn() → MoleculeAst
```

Maximum performance. Full type fidelity. But every DSL type needs hand-written
parser + formatter. Duplicates EDN syntax knowledge. High maintenance cost.

### Strategy D: Hybrid (Edn tree for read, streaming serde for write)

```
Write:  MoleculeAst → Serialize (with keyword sentinel) → EdnSerializer → EDN text
Read:   EDN text → parser → Edn tree → EdnDeserializer → MoleculeAst
```

Read path stays unchanged (Edn tree required for alias resolution).
Write path is fast (streaming). Keyword fidelity via sentinel on write only.

## Decision

Strategy A chosen as the starting point. Rationale:

- Correct by construction — no heuristics, no sentinels, no private protocols.
- The write-path allocation cost is bounded and measurable for molecule-sized data.
- Optimization path to Strategy B or D is clear if profiling shows a bottleneck.
- JSON/TOML/YAML serialization is not a concern — EDN is the sole format for the molecule DSL (decided in prior discussions).

Strategy B rejected as "magic" — the `__edn_keyword` sentinel is a private protocol
that could silently break. Acceptable as a future optimization if needed, not as the
foundation.

## AtomRef design

With Strategy A, `AtomRef` should become an enum that preserves the EDN type distinction:

```rust
enum AtomRef {
    Index(usize),   // Serializes as Edn::Int
    Named(String),  // Serializes as Edn::Keyword
}
```

`to_edn()` produces the correct Edn variant. `from_edn()` (or the existing
EdnDeserializer path) constructs the correct variant from the Edn type.

## Atom label preservation

Currently `MoleculeBuilder` and `Molecule` discard atom labels during lowering.
The `to_ast()` method cannot reconstruct `Atoms::Named` — it always emits
`Atoms::Indexed`.

To support faithful roundtrip of named atoms:
- `MoleculeBuilder` stores `atom_labels: Option<IndexMap<AtomIndex, String>>`
- Populated during `from_ast()` when input is `Atoms::Named`; `None` for `Atoms::Indexed`
- `to_ast()` emits `Atoms::Named` when labels exist, `Atoms::Indexed` otherwise
- Labels transfer to `Molecule` during resolution

## Open question: mixed-key atom maps

The current syntax is strictly either named or indexed. A useful extension would allow
selected atoms to be tagged while the rest use auto-generated indices:

```edn
{:atoms {0 "C" 1 "C" :ring-O "O"} :bonds [[0 :ring-O :single]]}
```

This requires mixed integer/keyword keys in the atoms map. The `Atoms` enum and
deserialization would need a third variant or a unified representation. Validation
must prevent collisions between explicit integer keys and positional indices.

Not yet designed. Deferred until the core roundtrip machinery is in place.

---

## Implementation summary (2026-04-07)

The decisions above were revisited during implementation. Strategy B was adopted
instead of A, after examining `serde_json::RawValue` as prior art.

### Strategy B: `EdnKeyword` sentinel in `umol-edn`

- `EdnKeyword(String)` type in `umol-edn`, modeled after `serde_json::RawValue`.
- Serializes via `serialize_newtype_struct("$edn::keyword", &self.0)`.
- `EdnSerializer` recognizes the sentinel, sets a `keyword_mode` flag, emits `:keyword`.
- Non-EDN serializers see a transparent string. No coupling to molecule domain code.
- `EdnKeyword` is never stored in domain types — constructed ephemerally at the serde boundary.

### `AtomRef` enum

- `AtomRef { Index(usize), Tag(String) }` — `Tag` holds plain `String`, not `EdnKeyword`.
- Serialize: `Index(n)` → integer, `Tag(s)` → ephemeral `EdnKeyword(s)` → `:keyword`.
- Deserialize: `visit_i64/u64` → `Index`, `visit_str` → `Tag`.
- `parse::<usize>()` heuristic eliminated.

### `Atoms` struct (replaces enum)

- `Atoms { entries: Vec<AtomAst>, tags: IndexMap<String, usize>, aliases: BiMap<String, AtomAst> }`.
- Vector-only `:atoms` format. Map form dropped.
- Inline tags: `[:tag-name "atom-def"]` in the vector.
- Aliases roundtrippable via `BiMap` reverse lookup (`get_by_right`). `AtomAst` derives `Hash`.

### Streaming read path

- `MoleculeInput` adapter: stream-deserialized via `EdnStreamDeserializer`, no Edn tree.
- `AtomEntryInput` enum: `Str(String)` | `Tagged(String, Box<AtomEntryInput>)`.
- `into_ast()` pipeline: alias resolution via table lookup, tag extraction, validation.
- `:aliases` / `:atom-aliases` both accepted on read.

### `BondSpec` eliminated

- `BondSpec` enum removed. All bond fields are `BondAst` directly.
- Built-in bond aliases (`:single` → order 1, etc.) handled via `builtin_bond_aliases()` `BiMap`.
- `BondAst::Serialize` checks alias table, emits `EdnKeyword` for known aliases.
- `lower_bond_spec` match arms in `from_ast` collapsed to `BondPattern::from_ast(bond, cfg)`.

### `Metadata` struct

- `Metadata { atom_tags, atom_aliases }` — formatting metadata separate from `Molecule`.
- `from_ast_with_metadata()` returns `(MoleculeBuilder, Metadata)`.
- `to_ast_with_metadata(&meta)` reconstructs tags and aliases.
- `Molecule` stays clean — no DSL formatting concerns.

### Spec and validation

- `:atoms` vector-only format, inline tags, `atom-ref ::= int | keyword`.
- Keyword namespace disjointness: tags, alias names, bond IDs must be mutually distinct.
- Alias names must not shadow element symbols.
- Topology graph module removed (premature, unused).

### Mixed-key atom maps (resolved)

The open question about mixed integer/keyword keys was resolved differently:
inline tags `[:tag "def"]` within a vector, not mixed-key maps. Tags are optional
per-atom — no need to tag every atom.
