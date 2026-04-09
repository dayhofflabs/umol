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

**Decision (2026-04-07): keep the lifetime.** Owned-only loses the
zero-copy substring story, which becomes load-bearing once `FromEdn::
from_edn_str` overrides start being written for hot types — that
optimization path needs to borrow from the source buffer, not clone out
of it. The native trait pair will be `FromEdn<'de>` / `ToEdn`, parallel
to toml-spanner's `Item<'de>`. Add an `EdnOwned = Edn<'static>` alias and
make `into_owned()` discoverable, but do not drop `'a`.

### 3. `to_value` requires constructing the full tree

This is the cost that motivates keeping a separate compact serde path. The
way out isn't to eliminate the tree; it's to make the tree path *cheap*:
- Arena-allocate `Edn` nodes (bump allocator) so building the tree is just
  pointer bumps.
- Make `Edn` itself a 16- or 24-byte enum (today it is likely larger because
  of `Cow`/`Box` payloads) so a tree of N nodes fits in cache.
- Then the "two paths" collapse to "one path with optional staging."

#### Benchmark reassessment (2026-04-07)

Measured `umol-edn/benches/parsing.rs` after the refactor:

| Operation                        | Tree path | Serde path | Ratio    |
|----------------------------------|-----------|------------|----------|
| Deserialize small → struct       | 582 ns    | 287 ns     | 2.03x    |
| Display small (MOLECULE_SMALL)   | 342 ns    | 142 ns     | 2.41x    |
| Display large (MOLECULE_LARGE)   | 3.01 µs   | 4.12 µs    | **0.73x** — tree faster |
| Pretty small (formatter)         | 491 ns    | (same)     | —        |
| Pretty large (formatter)         | 5.36 µs   | (same)     | —        |

Key finding: for serialization, **the tree path is faster on large inputs**.
Serde's per-field visitor dispatch overhead compounds with size; the tree
walker is a tight match loop over a prebuilt `Edn` enum, so the fixed tree
construction cost amortizes. The 2x gap only shows up on microbenchmark-
sized inputs where absolute cost is already hundreds of nanoseconds.

Revised conclusion: the second path is **not structurally necessary** for
throughput. The real reasons to keep it are narrower than originally stated:

1. **Streaming deserialization** — `StreamDeserializer` / `from_reader` can
   fail fast and skip allocation for never-read suffixes. Matters for huge
   files, not for molecule-sized inputs.
2. **Peak memory on very large inputs** — the tree path allocates an `Edn`
   node per parse node; direct serde allocates only the target struct. For
   a 100 MB input deserialized into a compact struct, tree path could use
   3–5x the peak memory.
3. **Small-value hot path** — direct serde wins ~300 ns on tiny values.
   Relevant only if parsing millions of tiny EDN blobs per second.

Reason 2 is real; 1 and 3 are weaker than the original review claimed.
For typical chemistry inputs (molecule-large and up) the tree path is
already competitive or faster, and within 2x on parse. Arena-allocating
`Edn` nodes would likely close the parse gap as well, leaving only the
streaming-prefix-parse use case as a genuine reason for two paths. At
that point the direct serde path becomes a specialized optimization for
`StreamDeserializer`, not a peer public API.

#### JSON-parity baseline (restored 2026-04-07)

Commit `63f3978b` (Apr 1) originally added a `serde_json` comparison to the
deserialize bench group; it was removed in `1ab0b940` (Apr 3). Restoring it
and re-running on current hardware confirms the original motivation for
keeping a direct streaming path:

| Path                                            | Time    |
|-------------------------------------------------|---------|
| `serde_json::from_str` (direct)                 | 252 ns  |
| `umol_edn::from_str` (direct streaming)         | **276 ns** |
| `serde_json::from_value` only                   | 312 ns  |
| `umol_edn::from_value` only                     | 338 ns  |
| `serde_json::from_str::<Value>` (tree only)     | 387 ns  |
| `umol_edn::read_string + from_value` (tree)     | 572 ns  |

Direct EDN streaming lands **24 ns behind serde_json direct** (+9.5%), on
input that is richer than the JSON equivalent (keywords, symbols, richer
literal set). This parity was a legitimate and hard-won design target —
it says "EDN can be a first-class peer to JSON for serde-speaking
consumers," which matters if the goal is ecosystem interoperability rather
than purely umol-internal use.

**What consolidation gives up.** A native `FromEdn` trait walking a
pre-built tree cannot match 276 ns. The minimum it can reach is roughly
`read_string` (~234 ns) + a direct trait walk (~100 ns estimated) ≈ 340 ns,
which is between the two current paths but short of JSON parity. Only
single-pass parser-deserializer *fusion* hits the 276 ns number, and that
is exactly what the direct streaming path is.

So the tradeoff is explicit, not hand-waved:

- **Keep direct streaming**: match serde_json on small inputs. Pay in code
  duplication, two-path maintenance, serde data-model limits on
  tags/sets/bignums, ~33K LLVM lines of per-type machinery in `umol-graph`.
- **Consolidate to tree + `FromEdn`**: lose ~50–100 ns on small inputs vs
  the JSON-parity baseline. Gain native EDN semantics, homoiconicity
  support, lower compile cost, single source of truth, and closed-form
  answers to §1's feature-gap list.

For umol-graph's actual workload (molecule DSL files parsed at load time,
not millions per second) the latency loss is immaterial. But if there is a
future use case where EDN parity with serde_json on tiny values matters —
parsing a stream of molecule fragments as a pipeline element, for example,
or positioning `umol-edn` as a general EDN library for serde consumers
outside umol — then the direct streaming path is the differentiator and
should be kept. That is a scope question, not a performance question.

#### Follow-up benchmarks (2026-04-07)

Two additional benches were added to quantify the memory and streaming
claims above, since the original microbenchmarks were time-only:

- `umol-edn/benches/memory.rs` — standalone binary with a tracking
  `GlobalAlloc` wrapper; measures peak bytes above baseline for each
  deserialization path across 10 KB, 100 KB, 1 MB, 10 MB inputs.
- `umol-edn/benches/parsing.rs::bench_stream_throughput` — criterion group
  comparing direct `StreamDeserializer` against `read_all` + per-element
  `from_value` on 1k and 10k concatenated records.

**Peak memory (Vec<Record> of concatenated MOLECULE_SMALL):**

| Input   | Direct `from_str` | Tree + `from_value` | Ratio |
|---------|-------------------|---------------------|-------|
|  10 KB  |      98 KB        |      223 KB         | 2.28x |
| 100 KB  |     1.10 MB       |     2.31 MB         | 2.10x |
|   1 MB  |    10.36 MB       |    23.26 MB         | 2.24x |
|  10 MB  |   100.48 MB       |   229.44 MB         | 2.28x |

The tree overhead is a flat **~2.2x, not the 3–5x speculated in reason 2**.
Per-record: direct ≈ 451 bytes/record (the `Record` struct itself), tree
≈ 1028 bytes/record. The ~580-byte gap is the `Edn` node + allocation
overhead per parse node. Arena allocation should compress this significantly
(bump-allocated `Edn` nodes don't carry per-allocation headers and share
one backing buffer), potentially bringing the ratio under 1.5x.

Also measured: `read_string` alone (no `from_value`) is indistinguishable
from `read_string + from_value` in peak. The deserialization step doesn't
add meaningful allocation beyond what the target `Vec<Record>` holds — the
tree is the whole cost.

**Stream throughput (records/sec via MiB/s):**

| Workload            | Direct stream | read_all + from_value | Ratio |
|---------------------|---------------|-----------------------|-------|
| 1k records (~39 KB) |  70.4 MiB/s   |       65.4 MiB/s      | 1.08x |
| 10k records (~390 KB) | 70.6 MiB/s  |       58.0 MiB/s      | 1.22x |

The direct stream path is **8–22% faster** in throughput — real but small.
It is *not* the order-of-magnitude advantage one might expect from "true
streaming vs eager batch." The `read_all + from_value` loop is within
striking distance on 1k records and only modestly behind on 10k. A lazy
iterator variant (`read_iter` returning `Iterator<Item = Edn>`) would
likely close most of the remaining gap by avoiding the intermediate
`Vec<Edn>` allocation in `read_all`.

**Reassessment summary:**

- Memory: real gap exists (~2.2x), smaller than feared, and most of it is
  per-node allocation overhead that an arena can compress.
- Throughput: direct streaming wins by 8–22%, not by an order of magnitude.
- The "keep two paths for peak memory and streaming" argument survives,
  but as a narrow optimization rather than a fundamental split. A lazy
  tree iterator + arena allocation would make the case for consolidation
  quantitatively strong.

**Still missing for a final decision:** an arena-allocated `Edn<'arena>`
prototype. The arena question cannot be benchmarked without writing it;
that is implementation work, not bench work. It is the load-bearing next
experiment.

#### MOL parity check (2026-04-07)

The relevant performance question for chemistry workloads is not "match
serde_json on tiny blobs" but "load a molecule DSL file at speed comparable
to MOL/CTfile." Measured against `umol-graph::io::ctfile`'s per-line
benchmarks:

| Bench                          | Time     |
|--------------------------------|----------|
| `mol_parsing/counts/valid`     | 39 ns    |
| `mol_parsing/atom/len69`       | 91 ns    |
| `mol_parsing/atom/len34`       | 50 ns    |
| `mol_parsing/bond/len21`       | 29 ns    |
| `mol_parsing/bond/len9`        | 16 ns    |
| `parse_collections/molecule_small` (EDN tree) | 378 ns |
| `parse_collections/molecule_large` (EDN tree) | 5.12 µs |

Estimated MOL parse for an equivalent 5-atom / 4-bond molecule:
`39 + 5×91 + 4×29 ≈ 610 ns` for line parsing, plus ~200–400 ns for struct
assembly and file boundary handling, total ~800 ns – 1 µs.

Estimated MOL parse for a `molecule_large`-equivalent (10 atoms, 11 bonds,
plus extended properties matching the EDN form's `config-overrides` and
`context` content): on the order of 2–4 µs once `M  CHG`, `M  RAD`, `M  CRS`,
`M  SDD` lines are parsed for each metadata field.

**Conclusion**: the EDN tree path is already at MOL parity for both small
and large molecule shapes. The 378 ns small-molecule parse is in the same
ballpark as a 5-atom MOL; the 5.12 µs large-molecule parse is within ~2x
of an equivalently-detailed MOL parse, with the gap accounted for by the
EDN form carrying *more* semantic content per byte (keywords vs fixed
columns, nested maps, aromatic hint vectors).

For a 10,000-molecule dataset load:

- MOL: ~10 ms total
- EDN tree: ~50 ms total
- EDN direct streaming (extrapolated): ~25 ms total

All three are imperceptible for interactive use and acceptable for batch
loading. Direct streaming buys real headroom but is not load-bearing for
chemistry workloads at the input shapes umol actually encounters.

#### Architecture this points to

"Easy things easy, hard things possible" lands on a single trait with a
default-implemented escape hatch:

```rust
trait FromEdn<'de>: Sized {
    /// Default path: walk a pre-built tree. Macro-derivable, MOL-parity speed.
    fn from_edn(edn: &Edn<'de>) -> Result<Self, EdnError>;

    /// Optional override for hot types: single-pass parser-deserializer fusion.
    /// Default impl materializes the tree first, then calls `from_edn`.
    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        let tree = parse_edn(input)?;
        Self::from_edn(&tree)
    }
}
```

- **Easy**: `#[derive(FromEdn)]` on `MoleculeAst`. Free, MOL-comparable,
  supports the full EDN data model (sets, tags, keywords, bignums).
- **Possible**: hand-write or specialize-derive `from_edn_str` for
  `MoleculeAst` if a benchmark ever shows tree construction is the
  bottleneck. Same architecture as the current direct streaming path,
  but EDN-native and per-type rather than parallel public API.
- **Possible**: tagged literals, sub-DSLs, isomer equilibria via tag
  readers and `Vec<Edn>` extension points (covered above).

The 24 ns JSON-parity number was real but it was bought at the cost of a
serde data model that fundamentally doesn't fit EDN — and that cost is
exactly the §1 limitations list. The native trait route preserves the
fast-path option (per-type, on demand) without requiring a parallel
public API to be maintained for every consumer.

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

## Reframing via toml-spanner (2026-04-07)

The `toml-spanner` crate (https://i64.dev/toml-spanner-no-we-have-serde-at-home/)
faced the same structural problem for TOML (datetimes, spans, inline-vs-block
style — none expressible through serde's data model) and chose to **bypass
serde entirely**, defining its own `FromToml<'de>` / `ToToml` traits plus a
derive macro. Three of its findings shift this review's conclusions:

1. **"You need an intermediate tree anyway."** toml-spanner's justification
   for a tree-first design: out-of-order parsing (and, for EDN, maps/sets)
   requires a staging tree regardless. Independent confirmation of §3's
   reassessment that the direct-serde path is not structurally necessary.

2. **Native traits close the transparent-derive gap.** §1 marked serde's
   limits as unrecoverable for users who want full fidelity from derive. A
   `FromEdn<'de>` / `ToEdn` trait pair takes `&Edn<'de>` directly and has
   its own derive macro — EDN-specific attributes (`#[edn(tag = "Point")]`,
   `#[edn(set)]`, `#[edn(keyword)]`, `#[edn(bignum)]`) become expressible
   without wrapper newtype gymnastics. This is a concrete, finite path, not
   a "documentation" fix. It also subsumes §11 — the traits naturally take
   `&Edn`, no clone needed for multi-target deserialization.

3. **The real payoff is compile time and binary size, not runtime.**
   toml-spanner reports 6x warm check, 7x warm build, ~10x less LLVM IR,
   3.8x smaller release binary size vs the `toml` crate. None of those
   wins come from SIMD, arenas, or interning — they come from removing
   generic indirection layers. Serde forces every consumer crate to
   re-monomorphize per-type deserializer machinery; a non-generic native
   trait is monomorphized once where defined.

### Measured on umol workspace (cargo llvm-lines)

- `umol-edn` lib standalone: ~63K LLVM lines, serde feature adds ~3K (~5%).
  Low because serde's tax is paid at *call sites*, not in the library
  defining the serializer.
- `umol-graph` lib: 785K LLVM lines total. Per-type streaming deserializer
  monomorphizations (`MapAccessor`, `EmptyMapAccessor`, `SeqAccessor` in
  `umol_edn::streaming`): **213 entries, ~12K lines**. Downstream visitor
  types from `umol_graph::dsl` (`MoleculeInputVisitor`, `AtomEntryInputVisitor`,
  `LocalizedBondVisitor`, etc.): **255 entries, ~16K lines**. All
  `deserialize_` / `visit_` entries combined: **522 entries, ~33K lines**
  (~4% of umol-graph's IR).
- `umol-params` lib: uses tree path (`from_value`), not serde streaming —
  zero streaming-accessor entries.

Smaller than toml-spanner's ratio because `umol-graph` is dominated by
SMILES/CTfile parsers, nom, and nalgebra, not by serde. But the 33K-line
tax is real, concentrated in exactly one crate, and exists *only* because
the direct streaming path is public API.

### Implication for the priority order

The consolidation case (collapse to tree path, add native traits) now rests
on three independent supports:

- **Runtime**: benchmarks show the tree path is competitive or faster on
  large inputs; memory gap is 2.2x, throughput gap 8–22% (§3).
- **Features**: native traits recover the transparent-derive gap for
  tags/sets/bignums that serde cannot express (§1).
- **Compile cost**: ~33K IR lines of per-type serde machinery in
  `umol-graph` alone would disappear. Check/build times improve; downstream
  crates pay no serde monomorphization tax.

The arena prototype answers a narrower question (parse throughput) that
only one of these three supports depends on, and only if the others fail.
It drops from "load-bearing next experiment" to "optimization to consider
after the trait layer exists."

## Homoiconicity: serde was never the right substrate

The preceding sections treat serde as something to be *replaced* reluctantly,
framing native traits as a workaround for limits. That framing is too
defensive. Serde was a convenience choice that leaks its core assumption;
the honest framing is that it was never the right substrate for an EDN-
backed homoiconic DSL in the first place.

### What serde actually offers

Exactly one thing of substance: **format-backend uniformity** — `Deserialize`
once, consume from JSON / TOML / YAML / bincode / MessagePack / ... without
rewriting the type. Everything else is replicable in finite work:

- **Derive automation** — `umol-edn-macros` already exists; `#[derive(FromEdn)]`
  is the same shape of work as `#[derive(Deserialize)]`.
- **Ecosystem gravity** — `clap`, `config`, `envy`, HTTP frameworks accept
  `Deserialize`. Real, but scoped to *config-shaped* types (CLI flags, env
  vars, settings files) — none of which are DSL types.
- **Zero-copy borrowed parsing** via `&'de str` — a native `FromEdn<'de>`
  gets this by construction.
- **Error plumbing conventions** — useful, not load-bearing.
- **Known attribute vocabulary** — users know `#[serde(rename = ...)]`; mild
  onboarding tax to learn `#[edn(...)]` equivalents.

None of these survive contact with the homoiconicity argument below.

### Why homoiconicity makes serde the wrong shape

Serde's data model *assumes* source format ≠ target type. Its premise is
that a `Deserializer` translates a foreign wire format into a distinct Rust
representation, and visitors exist to make that translation pluggable across
formats.

Homoiconic DSLs reject that premise. In EDN-backed DSLs:

- The surface syntax *is* the AST.
- `(molecule {:atoms [...] :bonds [...]})` is simultaneously program and data.
- `MoleculeAst` is not "the parsed form of the EDN" — it's a typed *view* of
  the same tree.
- Tagged literals (`#molecule`, `#isomer-equilibrium`, `#reaction`) are the
  extension mechanism for new DSL constructs; tag readers are the
  composition mechanism.

Serde has no concept of tagged literals, can't round-trip them, and forces
them to be encoded as sibling fields, wrapper structs, or `#[serde(untagged)]`
enum variants — none of which are homoiconic. You lose the ability to splice
partially-built AST fragments into larger forms, which is the core use case
when building higher-level DSLs on top of `MoleculeAst`.

### The wrapping scenario

Consider expressing isomer equilibria on top of the molecule DSL:

```edn
#isomer-equilibrium
{:species [#molecule {:atoms [C C ...] :bonds [...]}
           #molecule {:atoms [C C ...] :bonds [...]}]
 :temperature 298.15
 :populations [0.7 0.3]}
```

The natural representation is: `MoleculeAst` and `IsomerEquilibriumAst` both
implement `FromEdn`, both compose through tag dispatch, and each layer
validates its own slice. The key design choice at each layer is whether to
*eagerly* validate sub-forms (store `Vec<MoleculeAst>` and parse them at load
time) or to *defer* (store `Vec<Edn>` and validate on access). Deferred
validation is what makes the wrapper agnostic to which concrete sub-DSL is
spliced in — `#molecule`, `#conformer`, `#tautomer`, or a future construct
added by a downstream crate — without requiring the wrapper to recompile.

Note that this does **not** mean `MoleculeAst` itself becomes an `Edn`
newtype. It stays a typed struct with named fields (`atoms`, `bonds`,
`charge`, ...) and round-trips faithfully via `FromEdn` / `ToEdn`. The
homoiconic latitude shows up at two specific places: (1) forward-compatible
"unknown keys" fields that carry `Edn` values umol-graph does not recognize,
and (2) wrapper DSLs that choose to hold `Vec<Edn>` rather than
`Vec<MoleculeAst>` when they want deferred dispatch. The top-level struct
keeps static typing; the extension points keep tree fidelity.

Serde cannot do this cleanly because:

1. **No runtime tag dispatch.** `#molecule` would need a compile-time
   `Deserialize` impl, not a lookup in a tag-reader registry. New DSL
   constructs in a downstream crate would require modifying the enum.
2. **No deferred validation.** Serde's eager-parse model fights carrying
   unparsed sub-forms forward. You cannot say "this is a `MoleculeAst`
   *if and when* a molecule consumer asks for it."
3. **Closed enums for alternation.** `#[serde(untagged)]` doesn't compose
   across crates and doesn't support runtime extension.

Native `FromEdn` + tag readers sidesteps all three: the registry is runtime,
the intermediate representation is `Edn`, and new constructs are new
tag-reader registrations in whichever crate defines them.

### Consequence for the architecture

Keep serde support as a **compatibility shim** for config-shaped consumers
(a molecule file read by a general-purpose tool that expects `serde_json`-
style semantics), but treat `FromEdn` / `ToEdn` as the *primary* API for
DSL types. umol-graph's DSL types stop deriving `Deserialize` entirely and
derive `FromEdn` instead. Serde becomes a thin conversion layer on top of
the tree for external interop — not the way umol-graph's own types enter
or exit the EDN world.

This clarifies what "remove the direct serde path" actually means:

- **Gone**: `EdnStreamDeserializer`, per-type `visit_map::<MapAccessor>`
  monomorphizations, the `from_str` / `to_string` public API as a
  peer to the tree path, the dsl-type visitors in umol-graph.
- **Kept**: `Edn` tree + `FromEdn` / `ToEdn` as primary API; `from_value` /
  `to_value` as narrow interop with serde-ecosystem crates that still want
  the `serde::Deserialize` trait for their config types.

## Suggested priority order

1. **Prototype `FromEdn<'de>` / `ToEdn` traits + derive macro** for one
   representative type (e.g. `MoleculeAst`). Measure: round-trip fidelity
   for sets/tags/bignums, binary size delta, `cargo llvm-lines` delta on
   `umol-graph`. This is the load-bearing experiment — it validates both
   the feature-gap fix (§1) and the compile-cost payoff simultaneously.
2. Document the serde-vs-tree feature matrix (sets, tags, big numbers,
   namespaced symbols) — still valid, still cheap, useful regardless of
   whether (1) succeeds.
3. Decide on `Edn<'a>` vs owned-only. toml-spanner embraces the lifetime
   (`Item<'de>` is analogous); if native traits become the main path, the
   lifetime cost is paid by a narrower, more sophisticated audience and
   the "drop it entirely" argument weakens.
4. `Box<dyn Fn>` for `TagReaders` — cheap fix, unlocks real use cases.
5. Seal submodules; export only via crate root.
6. **Arena-allocated `Edn` nodes — demoted.** Only worthwhile if (1)
   succeeds and parse throughput becomes the dominant remaining cost. The
   benchmarks in §3 already suggest it isn't.

## Phase 2.5 measurement: fusion override closes the gap (2026-04-08)

The load-bearing experiment from §1 of the priority order has been run on
`MoleculeAst`, comparing three paths against the same five
chemistry-realistic inputs:

| Sample          | Serde stream | Native tree  | Native fused | Fused vs serde |
|-----------------|--------------|--------------|--------------|----------------|
| empty           | 179 ns       | 475 ns       | 103 ns       | −42%           |
| water           | 2.67 µs      | 3.07 µs      | 2.47 µs      | −7%            |
| tagged_ethanol  | 3.43 µs      | 3.94 µs      | 3.18 µs      | −7%            |
| aliased_benzene | 8.72 µs      | 9.63 µs      | 8.76 µs      | ±0%            |
| c20_chain       | 21.63 µs     | 23.29 µs     | 20.67 µs     | −4%            |

- **Native tree path** (the default `FromEdn` walking an `Edn` value) is
  7–17% slower than serde streaming on chemistry-realistic inputs and
  +167% on the trivial `empty` case. Cost is one allocation per `Edn`
  node plus a second walk to project to native types.
- **Native fused path** (`from_edn_str` override using
  `EdnStreamDeserializer`'s exposed primitives directly, no `Edn` tree)
  matches or beats serde on every input. The fusion is fully expressible
  with the public API added in this phase: `peek_byte`, `consume_byte`,
  `try_consume_byte`, `read_keyword_name`, `read_string`,
  `read_string_or_keyword`, `read_i64`, `read_skip_value`, `position`.

### Decision

Commit to the native `FromEdn` / `ToEdn` architecture. The fusion result
demonstrates that:

1. The default tree path is fast enough for cold types (within 17% of
   serde on realistic inputs) and ergonomic enough to author by hand or
   generate from a derive.
2. Hot types can opt into a hand-written or macro-generated `from_edn_str`
   override that beats serde streaming without depending on serde's data
   model.
3. The combination removes the structural objection to native traits:
   "fast paths require serde". They don't.

Whether `umol-edn` keeps the `serde` feature as an opt-in escape hatch
for external consumers is a secondary question, deferred until after
Phase 3 lands.

### Next steps (Phase 3)

1. Promote `parse_molecule_dsl_fused` to be the canonical
   `parse_molecule_dsl` entry point; retire the serde-coupled
   `MoleculeInput` deserialize visitors that are now redundant.
2. Design and prototype a `#[derive(FromEdn, ToEdn)]` macro that
   generates the equivalent of the hand-written fusion code from struct
   definitions. Target: have it produce code competitive with the
   hand-written `parse_molecule_dsl_fused`.
3. Extend the macro to cover the remaining DSL types (`AtomAst`,
   `BondAst`, `MoleculeAst` field-by-field) and remove the manual
   `Deserialize` impls in `umol-graph`.
4. Decide whether to drop the `serde` feature flag from `umol-edn`
   entirely, or keep it as a maintained escape hatch.

## Phase 3.5: serde compat layer redesign (2026-04-08)

Items 1–3 of the Phase 3 next-steps list are done. Item 4 — the fate of
the `serde` feature — is the topic of this section.

### Framing

The molecule DSL no longer goes through serde, so the serde feature has
no internal consumer. Its purpose is now exclusively external: let users
who bring serde-derived foreign types read/write EDN. The decision rule
that drove the redesign is **feature parity with the native path or
nothing**. Shipping a serde layer that lossily handles half the EDN data
model — keywords yes, sets and tagged literals no — is the
clojure-reader trap and is rejected.

### Architectural constraint: tree-walking only

The serde compat layer routes exclusively through
`Reader → Edn → EdnDeserializer`. The streaming serde
`Deserializer` impl on `EdnStreamDeserializer` is deleted entirely.
Consequences:

- One deserialize path to maintain instead of two. Every wrapper lands
  in `de.rs` only, not also in `streaming.rs`.
- `streaming.rs` becomes pure native primitives — every
  `#[cfg(feature = "serde")]` gate (27 sites, ~1500 LoC) goes away.
- `from_str_with(s, &ParseConfig)` keeps working because tag readers
  and parse config flow through `Reader`.
- Serde users pay one full pass over the input to build the `Edn` tree
  before deserialization. Phase 2.5 already showed native fused beats
  serde streaming on every input, so the casual-serde audience is not
  the audience that needs the fast path. Anyone who does have peak
  performance writes `FromEdn` directly.

### Required wrappers

Lossless feature parity requires one serde-reachable type per EDN
construct that does not have a native serde shape. The pattern is the
existing `EdnKeyword` newtype-token trick generalized to a dispatch
table over multiple tokens.

| EDN variant | Wrapper | State |
|---|---|---|
| Nil, Bool, Int, Float, Char, Str | native serde | done |
| Vector | `Vec<T>` (also accepts `Edn::List` for ergonomics) | done |
| Map | `HashMap` / `BTreeMap` / struct | done |
| Keyword | `Keyword` (current `EdnKeyword`, renamed) | exists |
| Symbol | `Symbol` | new |
| List | `List<T>` (strict opt-in for `Edn::List`) | new |
| Set | `Set<T>` | new |
| Tagged (enum-shaped) | `enum E { Variant(T) }` | exists |
| Tagged (dynamic) | `Tagged<T> { tag: String, value: T }` | new |
| BigInt | `BigInt` (cfg `bignum`) | new |
| BigDecimal | `BigDecimal` (cfg `bignum`) | new |
| Dynamic-typed lossless | `umol_edn::Value` | new |

Each wrapper claims a token name (`$edn::symbol`, `$edn::set`, …); the
serializer's `serialize_newtype_struct` and the deserializer's
`deserialize_newtype_struct` dispatch on the name and emit/consume the
correct `Edn` variant. Wrong-source-type errors at the boundary
(asking for `Symbol` from `Edn::Keyword` errors loudly).

### Implementation order

```
Phase 0 (delete streaming serde, reroute from_str)
   │
Phase 1 (generalize KEYWORD_TOKEN dispatch into a table)
   │
Phase 2 — wrappers, simplest first to validate the dispatch:
   ├─ 2.1 Symbol         (mirror of Keyword)
   ├─ 2.2 Set<T>         (first generic wrapper)
   ├─ 2.3 List<T>        (Vec<T> stays lenient on read)
   ├─ 2.4 Tagged<T>      (two-field tuple-struct token, the hard one)
   └─ 2.5 BigInt/Decimal (cfg-gated)
   │
Phase 3 (umol_edn::Value lossless mirror — mechanism TBD; may fall back
         to non-serde if serde-stable plumbing is too ugly)
   │
Phase 4 (parity test matrix: variant × wrapper × ser/de path)
   │
Phase 5 (re-exports under crate root, drop from_str_with, docs)
```

Phase 0 is reversible cleanup that should not change observable
behavior. Phase 1 generalizes the existing keyword dispatch. Phase 2
wrappers are mutually independent and can land as separate commits.
Phase 3 is the only step with an unproven mechanism (see "Open risks"
below).

### Net effect on the codebase

| Item | LoC delta |
|---|---|
| Delete streaming serde Deserializer in `streaming.rs` | −1500 |
| Six wrapper modules | +300 |
| Token dispatch in `de.rs` | +150 |
| Token dispatch in `ser.rs` (compact + tree) | +250 |
| `umol_edn::Value` mirror + impls | +300 |
| Parity test matrix | +300 |
| Docs | +50 |
| **Net** | **−150** |

End state: ~3500 LoC of serde compat (down from ~4250), single
tree-walking path, feature parity with native.

### Open risks

- **Phase 0 may surface latent test failures** (resolved 2026-04-08):
  no latent failures surfaced. All 790 tests pass after the tree-walking
  `from_value` path became the sole deserialize entry point.
- **`umol_edn::Value` plumbing** (resolved 2026-04-08): the
  carrier-pair mechanism. `Value::deserialize` uses `deserialize_any`.
  `EdnDeserializer::deserialize_any` routes standard variants (Nil,
  Bool, Int, Float, Char, Str, Vector, Map) through matching `visit_*`
  methods, and routes EDN-specific variants (Keyword, Symbol, List,
  Set, Tagged, BigInt, BigDecimal) through `visit_newtype_struct` with
  a purpose-built carrier deserializer that presents itself as a
  `(tag: &str, payload: Value)` tuple. `ValueVisitor::visit_newtype_struct`
  reads the pair and reconstructs the precise variant. Lossless over
  EDN; non-EDN deserializers degrade predictably (keyword/symbol →
  `Str`, list/set → `Vector`, tagged → tuple). Uses only standard
  serde machinery — no `Any` downcasting, no thread-local state, no
  byte smuggling.
- **`EdnTagged<T>` vs enum-variant tagged** (resolved 2026-04-08 by
  construction): `EdnTagged<T>` dispatches via
  `serialize_tuple_struct` / `deserialize_tuple_struct` keyed on the
  `TAGGED_TOKEN` name, while serde enum variants use
  `serialize_newtype_variant` / `deserialize_enum`. The paths never
  intersect. Phase 4 will still add a struct-holding-both spot test.
- **Wrapper composition with serde attributes** (`#[serde(default)]`,
  `#[serde(rename)]`, `#[serde(flatten)]`) must hold up. Phase 4
  spot-check.
