# A proper Rust EDN crate: motivation and approach

2026-03-31

## Context

We integrated `clojure-reader` 0.5.1 as our EDN reader/writer for the umol molecule
DSL. The crate's parser (`edn::read`) works correctly, but the serde layer required
~800 lines of workarounds in `dsl/edn_serde.rs` (see `discussion/62-*` for the full
issue list). Key problems: `deserialize_any` crashes on tagged literals, the serializer
ignores tag names, `Error` lacks `Clone`/`PartialEq`, `Display` delegates to `Debug`,
no string escaping, empty maps deserialize as unit.

## Decision: drop EDN tagged literals

The EDN spec reserves bare tags (`#atom`, `#bond`) for built-in use; user-defined tags
must be namespaced (`#umol/atom`). Rather than add 6 characters per atom/bond for
spec compliance, we dropped tags entirely. Atom and bond specs are plain strings —
context (`:atoms` map values, `:bond` fields, `:aliases` values) is sufficient for
disambiguation. This eliminated four of the ten workarounds.

## Strategic case for EDN

EDN has properties that JSON, YAML, and TOML lack, and that matter for where umol is
headed:

- **Homoiconicity.** Queries, rules, and data share the same format. A Datalog clause
  is just an EDN vector. You can store a query in the data it queries, compose rules
  programmatically, and serialize them without a separate query language.
- **Keywords as lightweight identifiers.** `:atom`, `:bond`, `:aromatic` — first-class
  values, not strings pretending to be enums. No quoting, no string-vs-enum impedance
  mismatch.
- **First-class sets.** Molecular properties, atom groups, ring membership. JSON arrays
  pretending to be sets are a constant source of bugs.
- **Tagged literals for domain types.** When needed (Datalog variables, rule references,
  domain-specific extensions), they're in the format, not hacked on top.
- **Extensibility without schema evolution.** A level-4 molecule can carry level-3 data
  without breaking parsers. Extra keys are ignored, not errors.

These properties are directly relevant to planned Datalog integration, level 3/4
buildout, and rule-based molecular transformations.

## Recommendation: write `umol-edn`

A clean-room EDN crate, not a subset — a proper implementation of the full spec.

### What exists today (reusable)

- `edn_serde.rs` is 70% of the serde integration (deserializer, serializer, seq/map
  access, error adapter). The main gap is that it wraps `clojure_reader::edn::Edn`
  instead of owning the type.
- The atom/bond/molecule DSL parsers use `nom`, which is the natural choice for an EDN
  parser too.

### What a proper crate needs

| Component | Estimate | Notes |
|---|---|---|
| `Edn` enum | Small | `Clone + PartialEq + Eq + Hash + Display`. Maps, vectors, lists, sets, keywords, symbols, strings, ints, floats, chars, bools, nil, tagged literals |
| Parser | ~400 lines | `nom`-based. Full grammar: all value types, `#_` discard, `;` comments, tagged literals, proper string escaping |
| Writer | ~150 lines | `Display for Edn` with proper string/char escaping |
| Error type | Small | Positions, structured variants, `Clone + PartialEq` |
| Serde deserializer | ~300 lines | Adapt from current `EdnDeserializer`. Handle tagged literals, enums, all numeric types |
| Serde serializer | ~250 lines | Adapt from current `EdnSerializer`. Proper `serialize_newtype_struct` for tags, string escaping |
| Tests | ~300 lines | Roundtrip, edge cases, conformance against the EDN spec |

Total: ~1500–2000 lines. Not a weekend, not a quarter.

### Parallel: PRs to `clojure-reader`

The fixes are straightforward and useful to others:
- `Display` that isn't `Debug`
- `Clone + PartialEq` on `Error` and `Code`
- `serialize_str` escaping
- `deserialize_any` handling `Tagged`
- `serialize_newtype_struct` emitting `#tag value`
- Empty map not deserializing as unit
- Bare tags in `deserialize_enum`

Submit these regardless. But don't block umol's timeline on upstream acceptance.

## Takeaways

- `clojure-reader`'s parser is solid; its serde layer is not production-ready.
- EDN tagged literals dropped from the molecule DSL (spec compliance, simpler code).
- Four simple struct types (`DativeBond`, `AromaticSystem`, `MulticenterBond`,
  `NoncovalentBond`) switched from hand-written visitors to `#[derive(Serialize,
  Deserialize)]` — custom visitors are only needed for polymorphic types.
- `recode_edn_error` extracts structured info from `clojure_reader::error::Error`
  instead of stringifying. `EdnParse(String)` remains as fallback.
- EDN is the right format for umol's direction (Datalog, rule systems, extensible
  molecular representations). Owning the implementation removes a fragile dependency
  and enables the features that make EDN worth using.
