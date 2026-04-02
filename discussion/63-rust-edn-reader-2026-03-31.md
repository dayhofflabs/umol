# A proper Rust EDN crate: motivation and approach

2026-03-31

## Context

We integrated `clojure-reader` 0.5.1 as our EDN reader/writer for the umol
molecule DSL. The crate's parser (`edn::read`) works correctly, but the serde
layer required ~800 lines of workarounds in `dsl/edn_serde.rs` (see
`discussion/62-*` for the full issue list). Key problems: `deserialize_any`
crashes on tagged literals, the serializer ignores tag names, `Error` lacks
`Clone`/`PartialEq`, `Display` delegates to `Debug`, no string escaping, empty
maps deserialize as unit.

## Decision: drop EDN tagged literals

The EDN spec reserves bare tags (`#atom`, `#bond`) for built-in use;
user-defined tags must be namespaced (`#umol/atom`). Rather than add 6
characters per atom/bond for spec compliance, we dropped tags entirely. Atom and
bond specs are plain strings — context (`:atoms` map values, `:bond` fields,
`:aliases` values) is sufficient for disambiguation. This eliminated four of the
ten workarounds.

## Strategic case for EDN

EDN has properties that JSON, YAML, and TOML lack, and that matter for where
umol is headed:

- **Homoiconicity.** Queries, rules, and data share the same format. A Datalog
clause is just an EDN vector. You can store a query in the data it queries,
compose rules programmatically, and serialize them without a separate query
language.
- **Keywords as lightweight identifiers.** `:atom`, `:bond`, `:aromatic` —
first-class values, not strings pretending to be enums. No quoting, no
string-vs-enum impedance mismatch.
- **First-class sets.** Molecular properties, atom groups, ring membership. JSON
arrays pretending to be sets are a constant source of bugs.
- **Tagged literals for domain types.** When needed (Datalog variables, rule
references, domain-specific extensions), they're in the format, not hacked on
top.
- **Extensibility without schema evolution.** A level-4 molecule can carry
level-3 data without breaking parsers. Extra keys are ignored, not errors.

These properties are directly relevant to planned Datalog integration, level 3/4
buildout, and rule-based molecular transformations.

## Recommendation: write `umol-edn`

A clean-room EDN crate, not a subset — a proper implementation of the full spec.

### What exists today (reusable)

- `edn_serde.rs` is 70% of the serde integration (deserializer, serializer,
seq/map access, error adapter). The main gap is that it wraps
`clojure_reader::edn::Edn` instead of owning the type.
- The atom/bond/molecule DSL parsers use `nom`, which is the natural choice for
an EDN parser too.

### What a proper crate needs


| Component          | Estimate   | Notes                                                                                                                                               |
| ------------------ | ---------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Edn` enum         | Small      | `Clone + PartialEq + Eq + Hash + Display`. Maps, vectors, lists, sets, keywords, symbols, strings, ints, floats, chars, bools, nil, tagged literals |
| Parser             | ~400 lines | `nom`-based. Full grammar: all value types, `#_` discard, `;` comments, tagged literals, proper string escaping                                     |
| Writer             | ~150 lines | `Display for Edn` with proper string/char escaping                                                                                                  |
| Error type         | Small      | Positions, structured variants, `Clone + PartialEq`                                                                                                 |
| Serde deserializer | ~300 lines | Adapt from current `EdnDeserializer`. Handle tagged literals, enums, all numeric types                                                              |
| Serde serializer   | ~250 lines | Adapt from current `EdnSerializer`. Proper `serialize_newtype_struct` for tags, string escaping                                                     |
| Tests              | ~300 lines | Roundtrip, edge cases, conformance against the EDN spec                                                                                             |


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

## EDN Spec

[GitHub](https://github.com/edn-format/edn)

## Takeaways

- `clojure-reader`'s parser is solid; its serde layer is not production-ready.
- EDN tagged literals dropped from the molecule DSL (spec compliance, simpler
code).
- Four simple struct types (`DativeBond`, `AromaticSystem`, `MulticenterBond`,
`NoncovalentBond`) switched from hand-written visitors to
`#[derive(Serialize, Deserialize)]` — custom visitors are only needed for
polymorphic types.
- `recode_edn_error` extracts structured info from
`clojure_reader::error::Error` instead of stringifying. `EdnParse(String)`
remains as fallback.
- EDN is the right format for umol's direction (Datalog, rule systems,
extensible molecular representations). Owning the implementation removes a
fragile dependency and enables the features that make EDN worth using.

## Additional concerns

- Which structures should hold map data? Should they be deterministic (even if
maps) are not required?
- Configurable formatting.

## Status of EDN outside Clojure

The honest assessment is that there's a lack of high-quality EDN implementations
for many languages, partly because the spec itself is not very formal or
complete (e.g., `clojure.edn` parses ratios as builtins despite the spec not
mentioning them). The ecosystem is thin compared to JSON/YAML/TOML, and most
libraries are hobby-grade or semi-abandoned.

**Most feature-complete implementations by language:**

**Rust:** The active choice is now `clojure-reader` (Grinkers), which aims to
match the behavior of Clojure's `tools.reader`. `edn-rs` (naomijub) is halted
and points to `clojure-reader` as successor. There's also
`bowbahdoe/edn-format`, which explicitly positions itself as more complete than
the other two crates, using `BTreeMap` for maps (ordered, deterministic).
`clojure-reader` supports serde, `no_std`, optional bigdecimal/num-bigint for
arbitrary precision, and ordered-float.

**Go:** `go-edn/edn` is complete, stable, and modeled on Go's `encoding/json` —
essentially a drop-in replacement for `encoding/json` patterns using struct
tags. Probably the most mature non-JVM implementation.

**Python:** `swaroopch/edn_format` implements all EDN features including custom
tagged elements, using PLY (lex/yacc). Last release was 0.7.5 in November 2020 —
effectively unmaintained. It returns `ImmutableDict` and `ImmutableList` by
default.

**Java (non-Clojure):** `bpsm/edn-java` is a pure Java parser and printer with
no external dependencies. This is actually one of the best-designed
implementations from an API perspective (more below).

**Haskell:** `hedn` is explicitly inspired by `Data.Aeson`, with
`ToEDN`/`FromEDN` typeclasses and a `Parser` monad. Well-typed but niche.

**JavaScript/TypeScript:** `edn-data` (jorinvo) works with plain JS data,
supports TypeScript, and has Node.js streaming support.

**Clojure-ecosystem fast path:** `edamame` (borkdude) is worth noting — highly
configurable, with location metadata, and used by prominent Clojure tooling
(malli, zprint, SCI, clerk). Not "outside" Clojure per se, but usable from
GraalVM native images.

## API Design Patterns

The interesting design decisions:

### 1. **bpsm/edn-java** — Builder-pattern configuration with `CollectionBuilder.Factory`

This is the most thoughtful API I've seen. The key insight: the parser
configuration is built via a `newParserConfigBuilder()` where you can override
the factory for each collection type — maps, sets, vectors, lists — via
`CollectionBuilder.Factory`. So you can swap in `TreeSet` for sets,
`LinkedHashMap` for maps, etc. The same pattern applies to tag handlers and
numeric type handlers (you can intercept `LONG_TAG`, `BIG_DECIMAL_TAG`, etc.).
The printer side mirrors this with a protocol-dispatch builder where you
register `Printer.Fn<T>` per Java class.

The crucial design: **collection construction is abstracted behind a factory +
builder interface**, not hardcoded. This cleanly separates "what the parser
emits" from "how it represents collections."

### 2. **edamame** — `:map` and `:set` constructor injection

Edamame lets you pass a `:map` option to inject an alternative map constructor —
e.g., `flatland.ordered.map/ordered-map`. Same for sets. Also has
`:auto-resolve` for namespace-qualified keywords, which accepts either a static
map or a function. This is the Clojure-idiomatic version of what edn-java does
with factories.

### 3. **edn-data (JS)** — Parse-time options for lossy vs. lossless

By default it returns a lossless representation (maps as
`{ map: [[key, value], ...] }`, keywords as `{ key: 'name' }`), but accepts
options like `{ mapAs: 'object', keywordAs: 'string' }` for a simplified
representation. This two-tier approach (faithful EDN AST vs.
convenient-but-lossy native structures) is pragmatic for languages without
keywords/symbols.

### 4. **hedn (Haskell)** — ADT + aeson-style typeclasses

The `Value` ADT has constructors for every EDN type
(`Nil | Boolean !Bool | String !Text | Character !Char | Symbol !ByteString !ByteString | Keyword !ByteString | Integer !Integer | ...`),
with a `Tagged a` wrapper (`NoTag !a | Tagged !a !ByteString !ByteString`).
Conversion uses `ToEDN`/`FromEDN` typeclasses with a `Parser` monad and
combinators like `(.:)` and `(.:?)` for map field access — directly mirroring
Aeson's API.

### 5. **Python edn_format** — Immutable wrappers + Keyword type

Returns `ImmutableDict` (implementing `collections.abc.Mapping`) and
`ImmutableList` (implementing `Sequence + Hashable`), with `Keyword` and
`Symbol` as distinct types that carry a `.name` attribute. The `loads`/`dumps`
API follows Python's `json` module convention. Honest but somewhat clunky — the
`ImmutableDict` is insertion-ordered but you have to do manual `isinstance`
dispatch to convert to native dicts.

## Key design axes for an EDN library

1. **Map representation:** `HashMap` (fast, unordered) vs.
  `BTreeMap`/`OrderedMap` (deterministic round-trip) vs. user-injectable
   factory. For spec compliance, insertion order shouldn't matter (EDN maps are
   unordered), but for human readability on round-trip, ordered maps are
   strongly preferred.
2. **Set representation:** Same trade-off. `HashSet` vs.
  `BTreeSet`/`OrderedSet`.
3. **Keyword/Symbol:** First-class newtype vs. prefixed string (`:foo` stored as
  string `":foo"`). The former is strictly better — you want
   `Keyword::namespaced("ns", "name")` with accessors, not string parsing at
   every use site.
4. **Tag handling:** Registry of `Tag → (Value → T)` handlers configured at
  parse time (edn-java, edamame) vs. returning a generic
   `TaggedValue(tag, inner)` for post-hoc dispatch. The registry approach is
   more ergonomic but couples parsing to application logic.
5. **Lossless AST vs. native projection:** Libraries either parse to a faithful
  AST enum (hedn's `Value`, edn-java's typed hierarchy) or directly into
   language-native containers (go-edn into structs, Python edn_format into
   dicts). The best APIs offer both: AST for tooling/round-tripping, native
   projection for convenience.
6. **Numeric tower:** EDN distinguishes int/float/bigint/bigdecimal/rational.
  Most implementations collapse these. The cleanest approach is feature-gated
   precision (like `clojure-reader` with optional `bigdecimal` + `num-bigint`
   crates).

If you're thinking about this for a Rust implementation (perhaps adjacent to
`umol`'s config or data interchange), the edn-java `CollectionBuilder.Factory`
pattern translates well to Rust traits, and `clojure-reader` is probably the
right starting point to evaluate or fork.

## Library API survey

Condensed view of the five implementations, focusing on what we'd want to steal.

### Core value types


| Library        | Value representation                                                       | Owns or borrows |
| -------------- | -------------------------------------------------------------------------- | --------------- |
| clojure-reader | `Edn<'e>` enum: 15 variants including `Tagged(&str, Box<Edn>)`             | Borrows input   |
| edamame        | Native Clojure data (no custom type)                                       | N/A             |
| hedn           | `Value` ADT (12 constructors), every node wrapped in `Tagged tag a`        | Owns (`Text`)   |
| go-edn         | `interface{}` + `Keyword`/`Symbol`/`Tag` newtypes                          | Owns            |
| edn_format     | Python natives + `Keyword`/`Symbol`/`Char`/`ImmutableDict`/`ImmutableList` | Owns            |


hedn's universal `Tagged` wrapper is the most principled — every value in every
collection carries an optional tag, so `#inst "2024-01-01"` inside a vector
doesn't need special cases. The tradeoff is verbosity: constructing untagged
values requires `NoTag(...)` everywhere. A middle ground: store tags in the
`Value` enum as a variant (`Tagged(tag, Box<Value>)`) like clojure-reader, but
ensure the serde layer handles it correctly (which clojure-reader does not).

### Collections


| Library        | Map                           | Set                    | Vector               | List            |
| -------------- | ----------------------------- | ---------------------- | -------------------- | --------------- |
| clojure-reader | `BTreeMap<Edn,Edn>`           | `BTreeSet<Edn>`        | `Vec<Edn>`           | `Vec<Edn>`      |
| hedn           | `Map TaggedValue TaggedValue` | `Set TaggedValue`      | `Vector TaggedValue` | `[TaggedValue]` |
| go-edn         | `map[interface{}]interface{}` | `map[interface{}]bool` | `[]interface{}`      | `[]interface{}` |
| edn_format     | `ImmutableDict`               | `frozenset`            | `ImmutableList`      | `tuple`         |


clojure-reader's `BTreeMap`/`BTreeSet` choice gives deterministic round-trip
ordering. This is the right default for snapshot testing and human readability.
go-edn uses hash maps (non-deterministic). Python uses `frozenset` for EDN sets
(correct semantics but unordered display). For umol-edn: `BTreeMap` and
`BTreeSet` as defaults, with the option to inject alternatives (via generic
parameter or builder, following edn-java's `CollectionBuilder.Factory` pattern).

### Keywords and symbols

All five store the full qualified name as a single string. hedn is the exception
for symbols: `Symbol Text Text` (namespace, name) as two fields. None split
keywords into namespace + name at the type level.

Recommendation: store as a single string internally (no allocation overhead for
the common unqualified case), expose `.namespace()` and `.name()` accessors that
split on `/`. A `Keyword` newtype and a `Symbol` newtype, both wrapping `String`
(or `&str` for borrowed parsing). Splitting at construction time wastes memory
for keywords that are only pattern-matched by full name.

### Tag handling

Three approaches in the wild:

1. **Data-only** (clojure-reader, hedn): tags are preserved as AST nodes,
  interpreted later. No dispatch during parsing.
2. **Registry dispatch** (go-edn, edn_format): `TagMap` / `add_tag()` maps tag
  names to handler functions, called during parsing.
3. **Configurable injection** (edamame): `:readers` option accepts a map from
  tag symbol to handler function.

For a serde-integrated Rust crate, (1) is the right base — parse to AST, let
serde dispatch on tags during deserialization. Offer an optional `TagRegistry`
for `edn::from_str`-style convenience where the user wants eager resolution.
This avoids coupling the parser to application types.

### Parse configuration


| Library        | Configuration                                                                                                    |
| -------------- | ---------------------------------------------------------------------------------------------------------------- |
| clojure-reader | None                                                                                                             |
| edamame        | Rich options map: `:readers`, `:auto-resolve`, `:map`/`:set` constructor injection, `:postprocess`, `:location?` |
| hedn           | Source name only                                                                                                 |
| go-edn         | `TagMap`, `MathContext`, `DisallowUnknownFields`                                                                 |
| edn_format     | `input_encoding` only                                                                                            |


edamame's configuration is the gold standard. The patterns directly useful for
umol-edn:

- `**:readers` → `TagRegistry`**: user-supplied tag → handler map.
- `**:auto-resolve` → `NamespaceResolver`**: for `::keyword` and
`::alias/keyword` resolution.
- `**:map`/`:set` constructor injection**: not needed if we hardcode `BTreeMap`
/ `BTreeSet`, but a generic parameter on the `Edn` type could allow
alternatives without runtime dispatch.

### Write / format configuration


| Library        | Formatting                                                              |
| -------------- | ----------------------------------------------------------------------- |
| clojure-reader | `Display` only, no options                                              |
| edamame        | No writer                                                               |
| hedn           | `renderText` only, no options                                           |
| go-edn         | `Marshal`, `MarshalIndent(prefix, indent)`, `MarshalPPrint(PPrintOpts)` |
| edn_format     | `dump(sort_keys, sort_sets, indent)`                                    |


go-edn has the richest printer: flat, indented, and pretty-printed modes.
edn_format adds `sort_keys` and `sort_sets` — directly useful for snapshot
stability. For umol-edn:

- `Display` for compact output (default).
- `EdnFormatter` with indent, sort_keys, sort_sets, line width. Analogous to
serde_json's `PrettyFormatter`.
- Deterministic output should be achievable without sorting if `BTreeMap` /
`BTreeSet` are the default collections.

### Serde integration

Only clojure-reader has serde. hedn has its own `ToEDN`/`FromEDN` typeclasses.
go-edn mirrors `encoding/json` with struct tags (`edn:"name,omitempty"`).

clojure-reader's serde gaps (the 10-issue list from discussion/62) define what
umol-edn must get right:

- `deserialize_any` must handle all variants including `Tagged`.
- `deserialize_enum` must support bare `:keyword` variants, not just
`#Type/Variant`.
- `serialize_str` must escape.
- `serialize_newtype_struct` must emit `#tag value`.
- Empty map `{}` must not deserialize as unit.
- Error type must derive `Clone + PartialEq`.

The serde layer is the hardest part of this crate and the most immediately
useful. It should be the first component brought to production quality after the
parser and `Edn` type.

## Parser infrastructure

### Requirements

- Composable combinators (not a grammar DSL).
- Zero-copy where possible (`&str` input, borrowed keywords/strings).
- Byte-position error reporting.
- Fuzzable (no panics on arbitrary input).
- Performance comparable to the nom-based MOL parser.
- ~15–20 productions (EDN values, strings, numbers, collections, comments,
`#_` discard, tagged literals).

### Framework comparison


|                | nom             | winnow           | chumsky       | pest         | hand-written |
| -------------- | --------------- | ---------------- | ------------- | ------------ | ------------ |
| Performance    | A               | A                | C+            | B            | A+           |
| Error messages | C               | B+               | A+            | B            | manual       |
| Combinator API | yes             | yes (cleaner)    | yes           | no (PEG DSL) | no           |
| Zero-copy      | yes             | yes              | partial       | partial      | yes          |
| Maintenance    | stable (frozen) | active (pre-1.0) | alpha rewrite | stable       | N/A          |
| EDN fit        | good            | best             | overkill      | poor         | good         |


**pest** — PEG DSL is the wrong paradigm for a library crate. Grammar in a
separate file, two-language friction, less control over allocation.

**chumsky** — error recovery is its selling point, but EDN doesn't need error
recovery (unlike a programming language). 2–5× slower than nom/winnow in
published benchmarks. Not worth the cost.

**nom** — proven in this project (MOL parser). But API is frozen at v7, error
handling requires significant boilerplate, and type signatures are painful.

**winnow** — nom's spiritual successor by epage (clap, toml_edit maintainer).
Same zero-copy architecture, same performance, but:

- `ContextError` with byte spans and `.context()` labels out of the box.
- `dispatch!` macro maps naturally to EDN's first-character dispatch (`(`→list,
`[`→vector, `{`→map, `#`→dispatch, `:`→keyword, `"`→string).
- `ErrMode::Cut` vs `ErrMode::Backtrack` gives committed parsing control.
- Already used by clojure-reader for the same grammar.

**hand-written** — justified when framework overhead is measurable (the SMILES
case: LALRPOP's table-driven dispatch dominated for a simple grammar). EDN's
grammar is simple enough that a combinator library adds negligible overhead, and
the error-reporting / composition benefits outweigh raw speed. If benchmarks
show winnow overhead is significant, a hand-written fallback is straightforward
for this grammar size.

### Recommendation: winnow

Same architecture as nom (no performance risk), better errors, actively
maintained, proven for this exact grammar (clojure-reader uses it). The pre-1.0
status is the only risk; mitigated by toml_edit's dependency ensuring API
stability. Fallback to nom v7 if winnow breaks compatibility.

### Validation strategy

- **Criterion benchmarks** from day one: parse a corpus of EDN files (reuse the
molecule DSL test fixtures), compare against `clojure_reader::edn::read`.
- **Fuzz with libfuzzer**: same setup as the existing SMILES fuzzer in
`umol-fuzz`. Feed arbitrary bytes, assert no panics.
- **Conformance suite**: hand-written EDN edge cases (nested collections, all
escape sequences, `#_` discard in every position, tagged literals, numeric
edge cases, comments in every position).

## Naming

The crate should be standalone and low-dependency. `umol-edn` ties it to umol,
which limits adoption and feels wrong for a general-purpose format library.
`edn` is taken (abandoned). Options:

- `**edn-rs`** — taken by naomijub (halted, points to clojure-reader).
- `**edn-format`** — mirrors the Python library name. Available on crates.io?
- `**reed`** — "read EDN". Short, memorable, available.
- **Just publish as `edn`** — the existing `edn` crate (v0.5, 2018) appears
abandoned. Could request ownership transfer.

Recommendation: develop as `umol-edn` internally, rename before publishing.
The name doesn't affect the implementation plan.

## Design decisions (2026-04-01)

### D1. Parser framework: winnow

Unless benchmarks show a performance wall. Fallback: nom v7 or hand-written.

### D2. Zero-copy borrowed strings

`Edn<'a>` borrows from input `&'a str`. Keywords, symbols, strings are
`Cow<'a, str>` — borrowed during parse-and-consume, clone-to-owned when the
value must outlive the input.

No string interning. The primary path is parse → serde → domain types — keywords
are consumed and discarded. If profiling shows keyword allocation pressure later
(e.g., Datalog engine), interning can be added behind a feature flag without
changing the `Cow` API.

### D2b. Numeric types in the Edn enum

The parser must choose concrete Rust types for the `Edn` enum independent of
serde. EDN has four numeric forms:


| EDN form         | Example | Default type                  | `bignum` feature                  |
| ---------------- | ------- | ----------------------------- | --------------------------------- |
| Integer          | `42`    | `i64` (error if out of range) | `i64`, fallback to `BigInt`       |
| Float            | `3.14`  | `f64`                         | `f64`                             |
| BigInt (`N`)     | `42N`   | parse error (unsupported)     | `BigInt`                          |
| BigDecimal (`M`) | `3.14M` | parse error (unsupported)     | `BigDecimal` (with `MathContext`) |


This matches clojure-reader's approach: `Int(i64)` by default,
`BigInt(num_bigint::BigInt)` and `BigDec(bigdecimal::BigDecimal)` behind
feature gates.

Rationals (`3/4`) are not in the formal EDN spec but are parsed by Clojure's
reader. Support as a feature flag: `Rational(i64, i64)`.

`i128` was considered for wider integer range without a dependency, but `i64`
matches Clojure's `Long` semantics and keeps the common case small.

**Integer overflow policy:** The EDN spec says "64-bit (signed integer)
precision is expected" but is silent on what to do when an unsuffixed literal
exceeds `i64` range. Three options:

1. **Parse error** — strict. If you want bigger, use `N`.
2. **Promote to `BigInt`** — lenient, matches Clojure's reader.
3. **Lossy truncation** — wrong, never do this.

Policy depends on feature flags, no configuration needed in the common case:

- **Without `bignum`**: overflow is always an error. No config option exposed.
`i64` is the only integer variant; there is nowhere to promote to.
- **With `bignum`**: overflow promotes to `BigInt` by default — if you opted
into arbitrary precision, you probably want it used. `overflow: Error`
available as opt-in for strict behavior even with bignum enabled.

**Narrowing to smaller types** (`u8`, `u32`, `u64`, etc.) is not a parser
concern. The `Edn` enum stores `i64`. Narrowing happens in the serde layer
(`deserialize_u8` checks range, errors if out of bounds) or via convenience
accessors (`.as_u8() -> Option<u8>`, `.as_u32() -> Option<u32>`, etc.) for
serde-free usage.

### D3. BTreeMap / BTreeSet by default

Deterministic round-trip ordering for snapshot testing and human readability.

Allow user-defined collection types via trait bounds or generic parameters on
the `Edn` type (following edn-java's `CollectionBuilder.Factory` and edamame's
`:map`/`:set` injection). All four collection types are injectable:

- **Maps, sets**: ordering and deduplication semantics vary (`BTreeMap` vs
`HashMap` vs `IndexMap`; `BTreeSet` vs `HashSet`).
- **Vectors, lists**: edamame doesn't inject these because Clojure's persistent
vector is the only sensible choice. In Rust, `SmallVec` and `ArrayVec` are
common alternatives — molecule DSL vectors (atom lists, bond tuples) are
typically small and benefit from stack allocation. This is a Rust-specific
concern with no analogue in the reference implementations.

### D4. Tagged as Value variant

`Tagged(String, Box<Edn>)` variant on the `Edn` enum. Not a universal wrapper
like hedn. Tags are uncommon in practice (dropped from the molecule DSL
entirely); wrapping every value in `Tagged`/`NoTag` is overhead for the common
case.

### D5. Keywords and symbols: single-string newtypes

`Keyword(Cow<'a, str>)` and `Symbol(Cow<'a, str>)` store the full qualified
name (`"ns/name"`). Accessors:

- `.name() -> &str` — everything after the last `/`, or the whole string.
- `.namespace() -> Option<&str>` — everything before the `/`, or `None`.

No allocation for the common unqualified case. Splitting at construction time
wastes memory for keywords that are only pattern-matched by full name.

### D6. Parsing configuration

Model on edamame's rich options. Initial API:


| Option           | Type                                                     | Purpose                                                                                                                                                                                                                                        |
| ---------------- | -------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `readers`        | `HashMap<String, Fn(Edn) -> Result<Edn>>`                | Tag dispatch during parsing (like edamame `:readers`)                                                                                                                                                                                          |
| `auto_resolve`   | `HashMap<String, String>`                                | Namespace aliases for `::alias/name` keywords                                                                                                                                                                                                  |
| `map_factory`    | trait / generic                                          | Custom map constructor (e.g., `IndexMap`)                                                                                                                                                                                                      |
| `set_factory`    | trait / generic                                          | Custom set constructor (e.g., `HashSet`)                                                                                                                                                                                                       |
| `vec_factory`    | trait / generic                                          | Custom vector constructor (e.g., `SmallVec`)                                                                                                                                                                                                   |
| `list_factory`   | trait / generic                                          | Custom list constructor (e.g., `ArrayVec`)                                                                                                                                                                                                     |
| `math_context`   | `MathContext { precision: u32, rounding: RoundingMode }` | Controls bigdecimal (`M` suffix) parsing precision and rounding. Only meaningful when the `bigdecimal` feature is enabled. Without it, `M` literals would silently lose precision by falling back to `f64`. Modeled on go-edn's `MathContext`. |
| `strict`         | `bool`                                                   | Strict EDN mode: reject Clojure extensions (extra escapes, `\formfeed`/`\backspace`, rationals, octal escapes). Default `false`.                                                                                                               |
| `duplicate_keys` | `DuplicateKeys::Error                                    | LastWins`                                                                                                                                                                                                                                      |


Deferred to tooling phase (see D10):


| Option        | Purpose                                          |
| ------------- | ------------------------------------------------ |
| `postprocess` | Per-form callback with `{obj, loc}`              |
| `location`    | Predicate filtering which forms get source spans |
| `uneval`      | Intercept `#_` discarded forms                   |


Not needed:


| Option                    | Reason                                                                                                                        |
| ------------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| `disallow_unknown_fields` | Serde already provides `#[serde(deny_unknown_fields)]` per target struct. Duplicating this at the parser level adds no value. |
| `input_encoding`          | Rust `&str` is UTF-8 by definition.                                                                                           |


### D7. EdnFormatter

Two-tier output, modeled on `serde_json::Formatter` trait and go-edn's
`MarshalPPrint`:

**Without serde** (operating on `Edn` values directly):

- `Display for Edn` — compact single-line output (default).
- `Edn::to_string_with(&EdnFormatter)` — configurable formatting.

**Through serde** (operating on `Serialize` types):

- `edn::to_string(v)` — compact.
- `edn::to_string_pretty(v)` — sensible default (2-space indent).
- `edn::Serializer::with_formatter(writer, EdnFormatter { .. })` — full
control. This is the serde_json pattern: serde itself has no formatting
config; formatting lives on the `Serializer` implementation via a `Formatter`
trait.

Both paths use the same `EdnFormatter` struct:


| Field        | Type            | Purpose                                                             |
| ------------ | --------------- | ------------------------------------------------------------------- |
| `indent`     | `&str`          | Indent string per nesting level (default: `" "`)                    |
| `line_width` | `Option<usize>` | Target line width for wrapping (None = no wrapping)                 |
| `sort_keys`  | `bool`          | Sort map keys (redundant with `BTreeMap`, needed for `HashMap`)     |
| `sort_sets`  | `bool`          | Sort set elements (redundant with `BTreeSet`, needed for `HashSet`) |


For reference — other format crates:

- **serde_json**: `Formatter` trait with `with_formatter()` on Serializer. Full
control. The gold standard.
- **toml_edit**: `to_string_pretty()` with fixed formatting. Real control
requires manipulating the lossless `Document` AST directly.
- **serde_yaml**: `to_string()` only, no formatting options exposed.

Deterministic output is free with `BTreeMap`/`BTreeSet` defaults. Formatting
options exist for the cases where users inject unordered collections.

Conformance suite snapshots should use a multi-line `EdnFormatter` for
diff-friendly `assert_snapshot!` output. This works identically whether the
value was constructed as `Edn` directly or serialized through serde.

### D8. Serde integration

Fix all 10 issues from discussion/62:

1. `deserialize_any` handles all variants including `Tagged`.
2. `deserialize_enum` supports bare `:keyword` variants.
3. `serialize_str` escapes properly.
4. `serialize_newtype_struct` emits `#tag value`.
5. Empty map `{}` stays a map.
6. Error type derives `Clone + PartialEq + Eq`.
7. `Display` is human-readable (not `Debug` delegation).
8. Sequences deserialize in correct order.
9. `from_str` rejects trailing content.
10. `serialize_unit_variant` emits `:keyword`.

### D9. Validation infrastructure

After basic parsing works:

- **Criterion benchmarks**: parse molecule DSL test fixtures, compare against
`clojure_reader::edn::read`.
- **Fuzz**: libfuzzer target in `umol-fuzz`, arbitrary bytes, assert no panics.
- **Conformance suite**: hand-written EDN edge cases (nested collections, all
escape sequences, `#_` discard in every position, tagged literals, numeric
edge cases, comments in every position, empty collections, `nil` in every
position).

### D10. Deferred: tooling API

Not needed for parse → serde → domain types. Add when building EDN-aware tooling
(linters, formatters, IDE support, source-mapping). These are well-understood
from edamame's design and can be added without breaking the core API.

`**postprocess`** — per-form callback. edamame signature:
`fn({obj, loc}) -> any`, called once for every parsed form (atoms, collections,
everything). The callback receives the parsed value and its source span
(`{row, col, end_row, end_col}`). Return value replaces the parsed object in the
tree. When set, edamame does NOT automatically attach location metadata —
the callback is responsible for wrapping or annotating. If `:source true` is
also set, the callback additionally receives the trimmed source string of that
form. Rust equivalent: `Fn(Edn, Span) -> Result<Edn>` or a trait with access
to source text.

`**location`** — predicate filtering which forms get source span metadata
attached. edamame signature: `fn(obj) -> bool`, called with the parsed object.
If truthy, location is attached via Clojure's `vary-meta`. Only called for
objects that support metadata (`IObj` — lists, vectors, maps, sets, symbols).
Without the predicate, all `IObj` forms get metadata. Rust equivalent: a filter
on `Edn` variant kind, or a `Fn(&Edn) -> bool` predicate. Requires the `Edn`
type to carry an optional `Span` field (or a parallel metadata map).

`**uneval`** — intercept `#_` discarded forms. edamame signature:
`fn({uneval, next}) -> any`. When set, edamame parses BOTH the discarded form
and the following form, then calls the callback with both. The return value
replaces both in the output. Without the callback, `#_form` is silently
skipped. Use case: formatters and linters that preserve `#_` comments in their
output. Rust equivalent: `Fn(Edn, Edn) -> Result<Edn>` where the first
argument is the discarded form and the second is the next form.

**Source spans in general** — the parser should track byte offsets internally
from the start (winnow provides this via `Located<&str>`). The question is
whether spans are stored on `Edn` nodes or returned separately. Storing them
on nodes adds a field to every value (16 bytes per node for start+end offsets).
Returning a parallel `SpanMap` avoids the per-node cost but complicates the
API. Decision deferred to when tooling needs materialize; the parser will track
positions regardless.

### D11. Spec conformance and extensions

**Strictness policy:** Clojure-compatible by default. Strict EDN as opt-in
(`strict: true` in parse config). The EDN spec is underspecified in several
areas, and real-world EDN files are written by Clojure programs that use
features beyond the spec text. Strict mode rejects anything not explicitly in
the spec (only 5 string escapes, no `\formfeed`/`\backspace` char literals, no
rationals, no octal escapes).

**Built-in tagged literals** (feature-gated, separate gates):


| Tag     | Feature  | Rust type                                      |
| ------- | -------- | ---------------------------------------------- |
| `#inst` | `chrono` | `chrono::DateTime<Utc>` or `chrono::NaiveDate` |
| `#uuid` | `uuid`   | `uuid::Uuid`                                   |


Without the feature, these parse as generic `Tagged("inst", Edn::Str(...))`.

**Character literals:**

Four named: `\newline`, `\return`, `\space`, `\tab`. Plus `\uNNNN` (4 hex
digits). Plus any single non-whitespace character (`\c`). All map to
`Edn::Char(char)`. `\`  (backslash + literal space) is invalid per spec.

Clojure adds `\formfeed` and `\backspace` — accept these for compatibility,
they map to `\u000C` and `\u0008`.

**String escaping:**


| Escape              | EDN spec | Clojure | umol-edn             |
| ------------------- | -------- | ------- | -------------------- |
| `\t \r \n \\ \"`    | yes      | yes     | yes                  |
| `\b \f`             | no       | yes     | yes (Clojure-compat) |
| `\uNNNN`            | no       | yes     | yes (Clojure-compat) |
| `\0`–`\377` (octal) | no       | yes     | yes (Clojure-compat) |
| `\'`                | no       | no      | no                   |


The EDN spec says "standard C/Java escape characters `\t, \r, \n, \\ and \"`
are supported" — only 5. But Clojure's reader silently extends this with `\b`,
`\f`, unicode, and octal escapes. Rejecting `\uNNNN` in strings would break
real-world EDN files. Accept the Clojure superset by default.

**Clojure reader extensions NOT supported** (strict EDN boundary):


| Feature                               | Rationale                        |
| ------------------------------------- | -------------------------------- |
| `#'var`                               | Clojure-specific, no EDN meaning |
| `@deref`                              | Clojure-specific                 |
| `#()` anonymous functions             | Clojure-specific                 |
| `'quote`, ``syntax-quote`,` ~unquote` | Clojure-specific                 |
| `#?` reader conditionals              | Clojure-specific                 |
| `#=` read-eval                        | Security risk, Clojure-specific  |


These could be added behind a `clojure-compat` feature flag later if there is
demand, following edamame's approach of per-extension opt-in.

**Rationals** (`3/4`): not in the EDN spec but parsed by Clojure's reader.
Behind `bignum` feature flag as `Rational(i64, i64)`.

### D12. Crate structure

**Feature flags:**


| Flag     | Adds                                                                                  | Default |
| -------- | ------------------------------------------------------------------------------------- | ------- |
| `serde`  | `Serialize`/`Deserialize` + `from_str`/`to_string` convenience                        | no      |
| `bignum` | `BigInt`, `BigDecimal` variants, integer overflow promotion, `MathContext`, rationals | no      |
| `chrono` | `#inst` → `DateTime` built-in handling                                                | no      |
| `uuid`   | `#uuid` → `Uuid` built-in handling                                                    | no      |


**Required dependency:** winnow only. Everything else feature-gated.

**Error type:** structured enum with byte-offset spans. `Clone + PartialEq + Eq`
(fixing the clojure-reader gap). Variants for: unexpected token, unexpected EOF,
integer overflow, invalid escape sequence, invalid character literal, duplicate
map key, trailing content. Byte offsets, not line/col (cheaper to track; line/col
computed on demand from the input).

### D13. Public API

**Input sources** (following serde_json's three-tier pattern):


| Source      | Type        | Borrowing                        | Notes                                                                                                                          |
| ----------- | ----------- | -------------------------------- | ------------------------------------------------------------------------------------------------------------------------------ |
| `&str`      | primary     | zero-copy, `Edn<'a>`             | Best performance.                                                                                                              |
| `&[u8]`     | convenience | zero-copy after UTF-8 validation | Useful for mmap'd files, network buffers.                                                                                      |
| `impl Read` | convenience | owned, `Edn<'static>`            | Reads into internal `String`, then parses. Can't borrow — all `Cow`s become owned. Same tradeoff as `serde_json::from_reader`. |


**Parsing entry points** (each available for all three input sources):


| Function              | Behavior                                                                         |
| --------------------- | -------------------------------------------------------------------------------- |
| `edn::read(s)`        | Parse one value, return `(Edn, &str)` remainder. EOF on empty input is an error. |
| `edn::read_string(s)` | Parse one value, reject trailing content.                                        |
| `edn::read_all(s)`    | Parse all values, return `Vec<Edn>`.                                             |
| `edn::Reader::new(s)` | Streaming iterator over `Result<Edn, EdnError>`.                                 |


With `serde` feature, additionally:


| Function                               | Behavior                                                                                 |
| -------------------------------------- | ---------------------------------------------------------------------------------------- |
| `edn::from_str::<T>(s)`                | Deserialize one value into `T`.                                                          |
| `edn::StreamDeserializer::<T>::new(s)` | Streaming iterator over `Result<T, EdnError>`. Mirrors `serde_json::StreamDeserializer`. |


`**Edn` accessors** — public enum, pattern matching is primary. Convenience
accessors for the common cases:

- Type checks: `.is_nil()`, `.is_keyword()`, `.is_map()`, etc.
- Narrowing: `.as_i64()`, `.as_f64()`, `.as_str()`, `.as_keyword()`, etc.
Return `Option<T>`.
- Numeric narrowing: `.as_u8()`, `.as_u32()`, `.as_u64()`, etc. Range-checked.
- Collection access: `.get("key")` for maps (looks up by keyword string),
`.iter()` for sequences.

`**edn!()` macro** — construct `Edn` values with EDN syntax at compile time:

```rust
let val = edn!({:name "ethanol" :atoms [:C :O :H]});
```

Proc macro in a separate `umol-edn-macros` crate, re-exported from `umol-edn`.
Parses EDN at compile time, emits `Edn` constructors. Avoids the verbosity of
manual `Edn::Map(BTreeMap::from([...]))` construction.

### D14. Round-tripping

Lossless round-trip (parse → print = identical text) is **not** a v1 goal. That
requires a CST preserving whitespace, comment positions, and formatting choices.
Semantic round-trip (parse → print → parse = identical `Edn` value) is a v1
requirement and should be covered by the conformance suite.

### D15. EDN spec clarification document

`spec/edn-spec.md` in the `umol-edn` crate. Not a new format — a precise
restatement of the existing EDN spec plus explicit decisions on every ambiguous
point. Same approach as the OpenSMILES spec in `umol-graph`. Each section links
to conformance test cases.

Ambiguities and gaps in the informal spec that need explicit decisions:


| Topic                     | Gap                                             | Decision                                                                                            |
| ------------------------- | ----------------------------------------------- | --------------------------------------------------------------------------------------------------- |
| Integer overflow          | "64-bit precision expected" — what on overflow? | D2b: error without `bignum`, promote with `bignum`                                                  |
| String escapes            | "Standard C/Java" — only lists 5                | D11: Clojure-compatible set, strict mode restricts to 5                                             |
| `\formfeed`, `\backspace` | Not in spec, accepted by Clojure                | D11: accept by default, reject in strict mode                                                       |
| Rationals                 | Not mentioned                                   | D2b: `bignum` feature flag                                                                          |
| Duplicate map keys        | Unspecified                                     | Configurable: error by default (safer), last-wins opt-in for Clojure-compat                         |
| `#_` nesting              | `#_ #_ a b` — discard two forms                 | Follow Clojure: each `#_` discards the next form, so `b` survives                                   |
| Whitespace                | Is `\f` whitespace? Is `,` whitespace?          | Follow Clojure: `,` is whitespace; `\f` is not                                                      |
| Tagged + `#_`             | `#_ #inst "2024"`                               | Follow Clojure: discard the entire tagged form                                                      |
| Numeric edge cases        | `+0`, `-0`, `00`, leading `+`                   | Follow spec grammar: `+` prefix allowed, `00` invalid (no leading zeros), `-0` valid and equals `0` |
| Empty symbol              | Is `/` alone valid?                             | Follow Clojure: `/` is the division symbol, a valid symbol                                          |
| Namespaced map            | `#:ns{:a 1}`                                    | Clojure extension, not in spec. Reject unless `clojure-compat` feature.                             |


### Additional quality checks

- llvm-cov
- clippy
- fuzz
- property testing -- generator for Edn objects -> test roundtripping

### Misc

- Remove incorrect comment:   // HashMap: iteration order is sorted, not insertion order.
- Remove fully qualified names from the code -> import types and traits
- stable or unstable sorting of maps and sets?
- use sorting in map / set tests, do not compare to specific ordering (or multiple orderings).
- how should the tag dispatch work?
- clarify that EDN uses UTF8.


### Spec definitions to verify

- [] The delimiters { } ( ) [ ] need not be separated from adjacent elements by whitespace.
- [] Tokens beginning with # are reserved. -> should error out
- [] # is not a delimiter. -> What does that actually mean? #a# is not allowed?
- [] Backslash cannot be followed by whitespace.
- [] Symbols begin with a non-numeric character and can contain alphanumeric characters and . * + ! - _ ? $ % & = < >. If -, + or . are the first character, the second character (if any) must be non-numeric. Additionally, : # are allowed as constituent characters in symbols other than as the first character. -> generate enough examples
- [] / has special meaning in symbols. It can be used once only in the middle of a symbol to separate the prefix (often a namespace) from the name, e.g. my-namespace/foo. / by itself is a legal symbol, but otherwise neither the prefix nor the name part can be empty when the symbol contains /. -> #foo/ , #/foo are illegal, #foo/bar and #/ are legal.
- [] If a symbol has a prefix and /, the following name component should follow the first-character restrictions for symbols as a whole. This is to avoid ambiguity in reading contexts where prefixes might be presumed as implicitly included namespaces and elided thereafter. -> #foo/1 is not allowed? #foo/_a is ok? 
- [] Per the symbol rules above, :/ and :/anything are not legal keywords. A keyword cannot begin with ::
- [] If the target platform supports some notion of interning, it is a further semantic of keywords that all instances of the same keyword yield the identical object. -> equality test?
- [] +num is allowed by the spec: "Integers consist of the digits 0 - 9, optionally prefixed by - to indicate a negative number, or (redundantly) by +. " -> fix the spec, add to suite
- [] No integer other than 0 may begin with 0 -> 01 is invalid, 0 is valid.
- [] -0 is a valid integer not distinct from 0.
- []  Integers -> generate examples

```
integer
  int
  int N
digit
  0-9
int
  digit
  1-9 digits
  + digit
  + 1-9 digits
  - digit
  - 1-9 digits
```

- [] Floats -> generate examples

```
floating-point-number
  int M
  int frac
  int exp
  int frac exp
digit
  0-9
int
  digit
  1-9 digits
  + digit
  + 1-9 digits
  - digit
  - 1-9 digits
frac
  . digits
exp
  ex digits
digits
  digit
  digit digits
ex
  e
  e+
  e-
  E
  E+
  E-
```

- [] In addition, a floating-point number may have the suffix M to indicate that exact precision is desired.
- [] A list is a sequence of values. Lists are represented by zero or more elements enclosed in parentheses (). Note that lists can be heterogeneous. -> Check (a b 1) list. -> Check list with ,
- [] A vector is a sequence of values that supports random access. Vectors are represented by zero or more elements enclosed in square brackets []. Note that vectors can be heterogeneous. -> Check [a b 1] -> Check commas
- [] A map is a collection of associations between keys and values. Maps are represented by zero or more key and value pairs enclosed in curly braces {}. Each key should appear at most once. No semantics should be associated with the order in which the pairs appear. -> Check heterogeneous maps. -> Check maps with different orderings, should give equal objects
- [] Note that keys and values can be elements of any type. The use of commas above is optional, as they are parsed as whitespace. -> Check keys of different types -> Check commas between pairs, inside of pairs
- [] A set is a collection of unique values. Sets are represented by zero or more elements enclosed in curly braces preceded by # #{}. No semantics should be associated with the order in which the elements appear. Note that sets can be heterogeneous. -> check heterogeneous sets. #{a b [1 2 3]} -> check different orderings give equal objects
- [] The semantics of a tag, and the type and interpretation of the tagged element are defined by the steward of the tag. #myapp/Person {:first "Fred" :last "Mertz"}.
- [] If a reader encounters a tag for which no handler is registered, the implementation can either report an error, call a designated 'unknown element' handler, or create a well-known generic representation that contains both the tag and the tagged element, as it sees fit. -> check error strategies. -> How can this system be implemented??
- [] Tag symbols without a prefix are reserved by **edn** for built-ins defined using the tag system. User tags ***must*** contain a prefix component, which must be owned by the user (e.g. trademark or domain) or known unique in the communication context. -> #tag "" is invalid. 
- [] A tag *may* specify more than one format for the tagged element, e.g. both a string and a vector representation.
- [] Tags themselves are not elements. It is an error to have a tag without a corresponding tagged element. -> "#tag" is invalid
- [] a `;` character is encountered outside of a string, that character and all subsequent characters to the next newline should be ignored. -> test
- [] `#` followed immediately by `_` is the discard sequence, indicating that the next element (whether separated from `#_` by whitespace or not) should be read and discarded. Note that the next element must still be a readable element. A reader should not call user-supplied tag handlers during the processing of the element to be discarded.
  `[a b #_foo 42] => [a b 42]`
- [] The discard sequence is not an element. It is an error to have a discard sequence without a following element.  -> This is EDN, not Clojure-dialect
  ```
- [] Sets and maps have requirements that their elements and keys respectively be unique, which requires a mechanism for determining when 2 values are not unique (i.e. are equal).
- [] nil, booleans, strings, characters, and symbols are equal to values of the same type with the same **edn** representation.
- [] integers and floating point numbers should be considered equal to values only of the same magnitude, *type, and precision*. Comingling numeric types and precision in map/set key/elements, or constituents therein, is not advised.
- [] sequences (lists and vectors) are equal to other sequences whose count of elements is the same, and for which each corresponding pair of elements (by ordinal) is equal.
- [] sets are equal if they have the same count of elements and, for every element in one set, an equal element is in the other.
- [] maps are equal if they have the same number of entries, and for every key/value entry in one map an equal key is present and mapped to an equal value in the other.
- [] tagged elements must define their own equality semantics. #uuid elements are equal if their canonic representations are equal. #inst elements are equal if their representation strings designate the same timestamp per ++[RFC-3339](http://www.ietf.org/rfc/rfc3339.txt)++.
- [] Exercise all format options.

