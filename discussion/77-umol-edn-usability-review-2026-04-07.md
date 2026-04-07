# umol-edn Usability Review

Date: 2026-04-07
Status: design notes; no implementation work scheduled.

Context: after merging the size-aware pretty-printing refactor (serde pretty
now routes through the `Edn` tree, sharing layout with `Edn::to_string_with`),
audit the wider public API surface for friction and inconsistencies.

## What works well

- **Three layers, clearly separated**: tree (`Edn` / `read_string`), serde
  compact (`from_str` / `to_string`), serde pretty (`to_string_pretty` /
  `to_string_with`). Each has a clear performance/ergonomics tradeoff.
- **`EdnFormatter` is the right abstraction**: every layout knob (width,
  sort, commas, indent, compact_*) lives in one config struct that flows
  through both pretty paths.
- **Single source of truth for layout**: after the recent refactor,
  `to_string_with` and `Edn::to_string_with` produce identical output,
  guarded by `test_pretty_path_matches_parse_then_format`.
- **`EdnKeyword` newtype** lets serde users round-trip keywords without
  leaking the lifetime-bound `Edn` enum.
- **`EdnKeyRef`** is the right shape for cross-lifetime map lookup without
  forcing a clone.
- **Streaming deserializer** exists (`StreamDeserializer`, `from_reader`) —
  important for large data.

## Real issues

### 1. The two parallel paths leak through the type system, not just performance

The serde path silently loses things the tree path preserves:
- arbitrary `#tag` values
- sets (`#{...}`)
- namespaced symbols
- `BigInt` / `BigDecimal`
- `char` distinct from one-character strings

A user who picks the serde API has no compile-time signal that those exist.
The crate doesn't document which path supports what — that asymmetry should
be a feature-matrix table at the crate root.

### 2. `Edn<'a>` lifetime parameter is viral

Most callers want owned `Edn<'static>`. The lifetime exists for zero-copy
parsing of `Cow::Borrowed("...")` strings, but it forces every downstream
type (`EdnMap<'a>`, `EdnSet<'a>`, `EdnKeyRef<'a>`) to thread `'a`.

Two cleaner alternatives:
- Drop the lifetime entirely; always own. Profile to see if `Cow` actually
  saves anything for typical inputs.
- Keep it but provide an `EdnOwned = Edn<'static>` alias and an
  `into_owned()` everywhere, prominently documented.

### 3. `to_value` requires constructing the full tree

This is the cost that motivates keeping a separate compact serde path. The
way out isn't to eliminate the tree; it's to make the tree path *cheap*:
- Arena-allocate `Edn` nodes (bump allocator) so building the tree is just
  pointer bumps.
- Make `Edn` itself a 16- or 24-byte enum (today it is likely larger because
  of `Cow`/`Box` payloads) so a tree of N nodes fits in cache.
- Then the "two paths" collapse to "one path with optional staging."

### 4. `TagReaders` uses `fn`, not `Box<dyn Fn>`

That blocks closures capturing context (e.g., a tag handler that needs a
registry of types). Either accept
`Box<dyn Fn(Edn) -> Result<Edn, EdnError> + Send + Sync>`, or expose both.

### 5. `DuplicateKeyPolicy` has only two variants (`Error`, `LastWins`)

Missing `FirstWins` and `Merge` (call a user fn). Not urgent, but the enum's
shape implies it is complete when it isn't.

### 6. `ParseConfig` and `EdnFormatter` are unrelated structs

Both configure I/O. Consider an `EdnConfig` that holds both, or at least
cross-link them in docs. Today a user has to discover them separately.

### 7. Error type opacity

`EdnError` is one big enum; users can't easily distinguish "syntax error at
line 5" from "tag handler failed" from "duplicate key" without matching
variants. A `kind()` method or splitting into `ParseError` / `SerError` /
`DeError` (like `serde_json`) would help. Connects to the crate-wide error
architecture work in `discussion/65-umol-error-handling-2026-03-31.md`.

### 8. `EdnMap` / `EdnSet` ordering is undocumented

Both are insertion-ordered via `Vec`-backed `HashMap`, but that is not
advertised in the public API. Clojure semantics say maps are unordered; if
we guarantee order, document it; if we don't, `sort_maps` in the formatter
is the only sensible default.

### 9. `edn!` macro uses `read_string` at runtime

A `proc-macro` variant that parses at compile time would catch typos and
avoid runtime cost. Optional.

### 10. No `Index` / `IndexMut` on `Edn`

Users will expect `edn["key"]` and `edn[0]` (serde_json has it). Missing
ergonomic accessors force tree-walking by hand.

### 11. `from_value` consumes `Edn`

There's no `from_value_ref(&Edn)`. Forces clones when you want to
deserialize the same value into multiple types.

### 12. Module visibility

`umol_edn::config::*`, `umol_edn::edn::*`, `umol_edn::error::*` are all
re-exported piecemeal at the crate root. The crate root should be the
*only* import path; submodules should be `pub(crate)` or sealed. Today
`use umol_edn::edn::Edn` and `use umol_edn::Edn` both work, which doubles
the API surface.

## Serde path limitations are structural, not implementation gaps

Most of the serde-path shortcomings listed in §1 cannot be fixed inside
`umol-edn` because they are limits of serde's data model, not of our
serializer. Worth recording so future contributors don't try to "fix"
them by reaching deeper into the serializer.

### Hard limits of serde's data model

Serde has no slot for:
- sets (only seq / map)
- arbitrary tags (only "newtype struct" and "enum variant" — both attach a
  *compile-time-known* name, not a runtime tag)
- distinct symbol vs keyword vs string (all collapse to `&str` on the wire)
- arbitrary-precision numbers (only fixed-width primitives)
- `char` distinct from a one-char string (serde has `char`, but most
  `Deserializer`s blur it)

A user's `Serialize` / `Deserialize` impl on an arbitrary Rust type cannot
signal "this is a set, not a vec" or "tag this with `#foo`" through the
standard serde traits — the traits don't have those slots. No work inside
`umol-edn` recovers these for arbitrary user types.

### Escape hatches that already exist

1. **Newtype tokens** — the `KEYWORD_TOKEN` trick `umol-edn` uses for
   `EdnKeyword` is generalizable. Define a wrapper (`EdnSet<T>`,
   `EdnTagged<T>`, `EdnBigInt`, …), give it a magic newtype name, and the
   serializer intercepts it. This recovers full fidelity *if the user opts
   into the wrapper type*. It is how `serde_json::Number`, `serde_bytes`,
   and `time::serde` all work.
2. **`Edn` as the universal type** — users who need full fidelity skip
   serde and work with `Edn` directly.

### What is fixable inside `umol-edn`

- More wrapper newtypes: `EdnSymbol`, `EdnSet<T>`, `EdnTagged<Tag, T>`,
  `EdnBigInt`, `EdnBigDecimal` — each with serde impls that round-trip via
  newtype tokens.
- Document the wrapper pattern so users know it exists.
- A derive macro (`#[derive(EdnSerialize)]`) that *augments* serde with
  EDN-specific attributes (`#[edn(tag = "Point")]`, `#[edn(set)]`) — the
  same approach `serde_with` takes for `serde_json` edge cases.

### What is not fixable

- Making a *plain* `HashSet<T>` serialize as an EDN set. Serde calls
  `serialize_seq` for it; the serializer cannot distinguish `HashSet` from
  `Vec` at that point.
- Making a user's enum variant emit as an arbitrary runtime `#tag`. Serde
  already owns enum encoding and reserves it for its own dispatch.
- Faithful round-trip of `Edn` itself through the serde traits without
  losing tag/set/bignum information — unless `Edn` becomes its own
  wrapper-token graph, which is exactly what `to_value` / `from_value`
  already are.

The gap is structural to serde, not to `umol-edn`. The wrapper-newtype
pattern can close most of it for users willing to opt in. The "transparent"
gap — serde-deriving an arbitrary Rust type and getting full EDN semantics
for free — cannot be closed without leaving serde's data model.

## Suggested priority order

1. Document the serde-vs-tree feature matrix (sets, tags, big numbers,
   namespaced symbols).
2. Decide on `Edn<'a>` vs owned-only — biggest source of friction.
3. `Box<dyn Fn>` for `TagReaders` — cheap fix, unlocks real use cases.
4. Seal submodules; export only via crate root.
5. Profile the tree-construction cost. If it is <2x compact serde, the
   "two paths" framing is overstated and the direct serde-compact path can
   be deprecated.
