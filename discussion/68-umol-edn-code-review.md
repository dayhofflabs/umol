# umol-edn Code Review (Adversarial)

Date: 2026-04-02

## Scope

Reviewed:

- `umol-edn/src/edn.rs`
- `umol-edn/src/collections.rs`
- `umol-edn/src/de.rs`
- `umol-edn/src/error.rs`
- `umol-edn/src/tags.rs`
- `umol-edn/src/parser.rs`
- `umol-edn/src/display.rs`
- `umol-edn/spec/edn-spec.md`
- `umol-edn/README.md`
- `umol-edn/Cargo.toml`
- `umol-edn/tests/conformance.rs`
- `umol-edn/tests/macros.rs`
- `umol-edn-macros/src/lib.rs`

Review stance: maximally adversarial (correctness-first, then API semantics, then performance and maintainability).

## Executive Findings

Highest-risk issues:

1. `Eq`/`Hash` contract violation for `Edn::Float` with NaN payloads.
2. Tagged unit enum deserialization accepts and silently drops arbitrary payloads.
3. Public `bignum` claims in spec and Cargo features are not implemented in runtime behavior.
4. Spec and implementation disagree on tagged-literal behavior and collection ordering semantics.

These are not cosmetic issues; they can cause incorrect behavior, silent data loss, and broken user expectations.

## Fix Response to Original Findings (Follow-up)

Status legend:

- `Resolved`: concern addressed in code and tests.
- `Partially resolved`: some work done, but residual issue remains.
- `Not resolved`: original concern still present.
- `Superseded`: original concern removed, but replaced by a different risk.

| Original finding | Status | Follow-up comment |
| --- | --- | --- |
| P0-1 Float `Eq`/`Hash` mismatch | `Resolved` | `Edn` now distinguishes NaN payloads consistently in `cmp` + `hash`, with dedicated tests in `edn.rs`. |
| P0-2 Tagged unit payload dropped | `Resolved` | `unit_variant` now requires `nil` payload and errors otherwise; tests cover error cases. |
| P0-3 `bignum` claimed but unimplemented | `Resolved` | `bignum` paths exist in parser/streaming and pass under `--features bignum`. |
| P1-1 Built-in tag behavior mismatch | `Partially resolved` | Tree parser now permits bare built-ins as `Tagged`, but streaming path still rejects unqualified built-ins unless a reader is registered. |
| P1-2 BTree vs hash storage docs mismatch | `Resolved` | Spec now states hash-based storage and explains deterministic formatting separately. |
| P1-3 Macro/parser divergence on special floats/chars | `Not resolved` | Macro still supports values/parser forms that diverge from reader/serializer behavior. |
| P2-1 Struct key handling | `Resolved` | Behavior is now explicit and tested for both value-tree and streaming deserialization pathways. |
| P2-2 Avoidable `deserialize_any` allocations | `Partially resolved` | Significant improvements landed, but there are still hot-path allocation opportunities in streaming/value conversion edges. |
| P2-3 `get(&str)` keyword-only alloc | `Resolved` | API renamed to `get_keyword`; hidden allocation concern reduced by borrowed keyword lookup path. |
| P2-4 Expensive collection `Ord` | `Accepted debt` | Cost is now documented; implementation remains intentionally allocation-heavy for deterministic compare semantics. |
| P2-5 `expect(...)` panics in map access | `Resolved` | Panics were replaced with error returns. |
| P3-1 `unwrap_err` offset quality | `Partially resolved` | Offset handling improved, but `Incomplete` mapping still has ambiguity risk depending on parser context. |
| P3-2 `unreachable!()` in `cmp` | `Superseded` | Panic was removed; fallback now returns `Ordering::Equal`, which avoids panic but introduces a potential silent-ordering risk if invariants drift. |
| P3-3 README feature drift | `Partially resolved` | Feature table now includes `bignum`, but description still says "not yet implemented", which is stale. |

## Severity Legend

- `P0`: correctness/safety contract breach; can produce wrong behavior now.
- `P1`: spec/API semantic break or major reliability risk.
- `P2`: significant inefficiency/complexity or sharp-edge behavior.
- `P3`: maintainability debt, API roughness, or lower-impact issues.

## Findings

### P0-1 (`Resolved`): `Eq` and `Hash` are inconsistent for `Edn::Float`

**Evidence**

- `umol-edn/src/edn.rs`: `Ord` compares floats with `total_cmp` (`Edn::Float(a), Edn::Float(b) => a.total_cmp(b)`).
- `umol-edn/src/edn.rs`: `PartialEq` is implemented via `self.cmp(other) == Ordering::Equal`.
- `umol-edn/src/edn.rs`: `Hash` uses `f.to_bits().hash(state)` for floats.

**Why this is a bug**

With `total_cmp`, NaNs compare equal in ordering terms, but `to_bits()` differs across NaN payloads. That can produce equal keys with different hashes, violating Rust hash map/set key invariants.

**Impact**

- Incorrect behavior in `FxHashMap<Edn, _>` and `FxHashSet<Edn>` for float NaN keys.
- Hard-to-reproduce lookup failures and duplicate-key anomalies.

**Fix direction**

- Canonicalize NaN before hashing (single NaN bit pattern), or
- Redefine equality/hash semantics for floats so both operate on the same canonical representation.
- Add tests using `f64::from_bits` with different NaN payloads.

---

### P0-2 (`Resolved`): Tagged unit variant deserialization discards payload silently

**Evidence**

- `umol-edn/src/de.rs`: `EdnTaggedVariantAccess::unit_variant` returns `Ok(())` unconditionally.

**Why this is a bug**

Input like `#Variant 123` can deserialize into a unit enum variant without error, silently dropping `123`.

**Impact**

- Silent acceptance of malformed data.
- Data-loss behavior during deserialization.

**Fix direction**

- Reject non-empty/non-nil payloads for unit variants with explicit error.
- Add tests for `#Unit 1`, `#Unit {:a 1}` and assert error.

---

### P0-3 (`Resolved`): `bignum` is claimed but not implemented

**Evidence**

- `umol-edn/spec/edn-spec.md`: sections on `N`/`M` suffixes claim `BigInt`/decimal behavior with `bignum`.
- `umol-edn/Cargo.toml`: feature `bignum = ["dep:num-bigint", "dep:bigdecimal"]`.
- `umol-edn/src/parser.rs`: `N`/`M` suffix always returns `UnsupportedFeature { feature: "bignum" }`.
- No `num_bigint` / `bigdecimal` usage in `umol-edn/src`.

**Why this is a bug**

Public contract mismatch: users can enable `bignum` but cannot get bignum behavior.

**Impact**

- Broken expectations and misleading feature surface.
- Dead optional dependencies.

**Fix direction**

Pick one and enforce consistently:

1. Implement actual bignum parsing and AST representation, or
2. Remove `bignum` feature/dependencies and delete spec claims.

---

### P1-1 (`Partially resolved`): Spec contradiction for built-in unqualified tags without features

**Evidence**

- `umol-edn/spec/edn-spec.md`: states built-in tags parse as `Tagged(...)` without corresponding feature.
- `umol-edn/src/parser.rs`: rejects unqualified tags unless reader is registered.
- `umol-edn/src/config.rs`: default readers only register `inst`/`uuid` when features are enabled.

**Why this matters**

Feature-off behavior for `#inst`/`#uuid` is documented one way and implemented another.

**Impact**

- Behavioral surprises in production with feature toggles.
- Tests/docs confusion.

**Fix direction**

Either:

- Update parser behavior to allow these bare built-ins as plain `Tagged`, or
- Update spec/docs to state they are rejected unless feature-enabled reader exists.

---

### P1-2 (`Resolved`): Spec says BTree deterministic ordering; code uses hash collections

**Evidence**

- `umol-edn/spec/edn-spec.md`: claims `BTreeMap` / `BTreeSet` deterministic order.
- `umol-edn/src/collections.rs`: uses `FxHashMap` / `FxHashSet`.
- `umol-edn/src/display.rs` tests accept both map/set orderings.

**Why this matters**

Public documentation states deterministic iteration; implementation is explicitly unordered.

**Impact**

- Users may rely on stable formatting/iteration and get nondeterminism.
- Interop and snapshot-test friction.

**Fix direction**

Either:

- Switch internals to ordered collections, or
- Rewrite spec/docs to reflect unordered internals and define deterministic output only where explicitly implemented.

---

### P1-3 (`Not resolved`): Macro/parser behavior diverges on special floats and char names

**Evidence**

- `umol-edn-macros/src/lib.rs`: supports `##NaN`, `##Inf`, `##-Inf`; supports char names `formfeed`, `backspace`.
- `umol-edn/tests/macros.rs`: tests special floats succeed.
- `umol-edn/tests/conformance.rs`: parser rejects special floats; parser supports only `newline`, `return`, `space`, `tab`.
- `umol-edn/src/display.rs`: panics when asked to display NaN/Inf.

**Why this matters**

Two official construction paths (`read_string` vs `edn!`) produce different value spaces and round-trip properties.

**Impact**

- Non-obvious invariant break: some values can be created but cannot be serialized.
- Potential runtime panics from macro-created values.

**Fix direction**

- Align macro and parser semantics, or
- Explicitly document macro-only extensions and non-serializable values.

---

### P2-1 (`Resolved`): Struct deserialization silently ignores non-string-like map keys

**Evidence**

- `umol-edn/src/de.rs`: `EdnStructMapAccess` skips keys not in `Keyword|Symbol|Str`.
- `umol-edn/src/de.rs` test `test_deserialize_struct_ignores_non_string_keys` codifies behavior.

**Why this is risky**

Malformed config input may deserialize successfully with fields omitted, instead of failing loudly.

**Fix direction**

- Prefer strict mode by default (error on non-struct keys) and optionally offer lenient mode.

---

### P2-2 (`Partially resolved`): Avoidable allocations in `deserialize_any`

**Evidence**

- `umol-edn/src/de.rs`: keywords/symbols converted to owned `String` with `.to_string()`.
- `umol-edn/src/de.rs`: sets copied into `EdnSeq` via `collect` before visiting.

**Impact**

- Higher allocation cost for common paths.

**Fix direction**

- Use borrowed visitor paths when possible.
- Introduce dedicated set `SeqAccess` adapter to avoid intermediate `EdnSeq`.

---

### P2-3 (`Resolved`): `Edn::get(&str)` is keyword-only and allocates every call

**Evidence**

- `umol-edn/src/edn.rs`: `get(&str)` looks up only `Edn::Keyword(Keyword::owned(key.to_string()))`.

**Risk**

- Misleading API name (`get`) for a narrow key type.
- Repeated allocations in map-heavy reads.

**Fix direction**

- Rename to `get_keyword` (or equivalent), and add broader lookup API by `&Edn`.
- Consider borrowed lookup strategy.

---

### P2-4 (`Accepted debt`): Collection `Ord` implementations are expensive

**Evidence**

- `umol-edn/src/collections.rs`: `EdnMap::cmp` and `EdnSet::cmp` allocate vectors and sort on every comparison.

**Impact**

- `O(n log n)` with allocation per compare; poor behavior when sorting/comparing many nested EDN values.

**Fix direction**

- Document as intentionally expensive if only used for deterministic ordering logic.
- Revisit ordering strategy if `Ord` is required in hot paths.

---

### P2-5 (`Resolved`): Panics in serde map accessors (`expect(...)`)

**Evidence**

- `umol-edn/src/de.rs`: `next_value_seed` in `EdnMapAccess` and `EdnStructMapAccess` uses `expect(...)`.

**Risk**

- Panics if serde contract is violated by custom seeds/visitors.

**Fix direction**

- Replace `expect` with explicit `EdnError::Custom(...)` return.

---

### P3-1 (`Partially resolved`): `unwrap_err` loses offset context on incomplete input

**Evidence**

- `umol-edn/src/error.rs`: `ErrMode::Incomplete(_) => UnexpectedEof { offset: 0 }`.

**Impact**

- Poor diagnostics at exact failure location.

**Fix direction**

- Preserve best-known offset when mapping incomplete errors.

---

### P3-2 (`Superseded`): `unreachable!()` in `Edn::cmp` is a maintenance landmine

**Evidence**

- `umol-edn/src/edn.rs`: wildcard arm in `match (self, other)` after discriminant compare.

**Risk**

- Future variant changes can introduce panic if `variant_ord` and `cmp` updates drift.

**Fix direction**

- Keep matches structurally exhaustive without wildcard fallback; avoid runtime trap.

---

### P3-3 (`Partially resolved`): Documentation drift in README feature table

**Evidence**

- `umol-edn/README.md` feature table omits `bignum`, while `Cargo.toml` defines it.

**Impact**

- Confusing crate capability signaling.

**Fix direction**

- Sync README with real feature surface after resolving `bignum` implementation status.

## Test Coverage Gaps

Missing tests that should be added before further refactors:

1. `Eq`/`Hash` coherence tests for float edge cases:
   - Distinct NaN payloads via `from_bits`.
   - `+0.0` vs `-0.0` behavior for equality, ordering, hashing.
2. Tagged unit variant strictness tests (`#Variant` with non-empty payload must error).
3. Feature-matrix tests for `chrono`/`uuid` with and without features.
4. Set duplicate parse semantics tests (`#{1 1}` behavior explicitly asserted).
5. `bignum` positive-path tests (if feature remains), otherwise delete all bignum claims/tests.

## Technical Debt Inventory

### Correctness debt

- Float equality/hash mismatch.
- Unit variant payload dropping.

### Contract/documentation debt

- Spec contradicts code on collection internals and built-in tags.
- README feature drift.

### Complexity/performance debt

- Expensive map/set `Ord`.
- Avoidable allocs in deserializer.
- Narrow `get(&str)` helper with hidden allocation.

### Robustness debt

- `expect` panics in map accessors.
- Lossy error mapping for incomplete parsing.

## Recommended Remediation Order

1. Fix float `Eq`/`Hash` coherence.
2. Fix tagged unit variant payload validation.
3. Resolve `bignum` truth (implement or remove).
4. Align spec + README with actual parser behavior.
5. Address deserializer strictness/perf issues.
6. Add regression tests for all above before broader refactors.

## Residual Risk if Unchanged

- Hash-based containers can behave incorrectly with NaN keys.
- Malformed tagged inputs can be accepted silently.
- Users will continue to build against contradictory docs/spec behavior.
- Macro/parser divergence can continue to produce non-serializable values and panic paths.

## Fresh Adversarial Findings (Second Pass)

This section intentionally ignores prior findings and reports only current-code risks.

### P0

1. **Unsafe lifetime cast in map lookups (`collections.rs`)**
   - `EdnMap::get` and `contains_key` use raw pointer lifetime widening:
     - `unsafe { &*(key as *const Edn<'k> as *const Edn<'a>) }`
   - This is a soundness hazard: correctness depends on undocumented lifetime/layout invariants that the type system does not enforce.

2. **Discard parser swallows parse errors (`parser.rs`)**
   - `ws_and_comments` ignores all errors when parsing discarded form after `#_`.
   - This can desynchronize parser state and hide malformed discarded forms.

3. **Streaming discard skip-string escape bug (`streaming.rs`)**
   - `skip_string` always skips escape as backslash + one byte.
   - `\uXXXX` escapes are not skipped correctly in discard/ignored paths.

### P1

1. **Tree vs streaming mismatch for built-in bare tags**
   - Tree parser allows bare `inst`/`uuid` fallback to `Tagged`.
   - Streaming parser rejects unqualified tags unless a reader is registered.

2. **`DuplicateKeyPolicy` not applied in streaming map path**
   - Tree parse enforces duplicate key policy.
   - Streaming map accessor has no duplicate-key enforcement semantics.

3. **README bignum statement is stale**
   - README still says `bignum` is "not yet implemented" while implementation/tests exist.

### P2

1. **`Ord` fallback can silently equate distinct variants (`edn.rs`)**
   - `cmp` fallback now returns `Ordering::Equal` on invariant drift path.
   - This avoids panic but can silently violate ordering semantics if variant matching diverges.

2. **`parse_all` ignores whitespace/comment parser errors (`parser.rs`)**
   - Uses `.ok()` in loop, suppressing errors that should propagate.

### Evidence references for second pass

- `umol-edn/src/collections.rs` (`EdnMap::get`, `contains_key`)
- `umol-edn/src/parser.rs` (`ws_and_comments`, `edn_dispatch`, `parse_all`)
- `umol-edn/src/streaming.rs` (`skip_string`, `skip_tag_if_present`, map access path)
- `umol-edn/src/edn.rs` (`cmp` fallback)
- `umol-edn/README.md` (feature table entry for `bignum`)

## Claude Code Review (Third Pass)

Date: 2026-04-02

Adversarial review focused on correctness, performance, unnecessary complexity, and technical debt.

### Resolved

| # | Issue | Severity | Fix |
|---|-------|----------|-----|
| 1.1 | `Display` panics on NaN/Infinity floats | High | Returns `Err(fmt::Error)` instead of `assert!()`. `display.rs:51-53` |
| 1.3 | No recursion depth limit in winnow parser | High | `ParseCtx` with `Cell<u16>` depth counter, `MAX_DEPTH=128`. `enter_scope`/`leave_scope` in all four collection parsers. `ws_and_comments` returns `()` (zero overhead on atom parsing). `parser.rs:26-54` |
| 1.7 | Serde char serialization is lossy | High | `serialize_char` now handles `\newline`, `\return`, `\space`, `\tab`, `\uNNNN`. `ser.rs:105-119` |
| 1.6 | Streaming parser accepts leading-zero integers | Medium | `scan_number_str` rejects `007` etc. for non-decimal integers. `streaming.rs:191-196` |
| 2.3 | Serializer allocates temporary strings for integers/floats | Low | `serialize_i64`/`serialize_u64` use `itoa::Buffer`; `serialize_f64` uses `write!` directly. `ser.rs:72-96` |
| 2.4 | `EdnSeq::From<Vec>` round-trips through iterator for no reason | Low | `Self(v)` instead of `Self(v.into_iter().collect())`. `collections.rs:306` |

### Remaining: worth fixing

| # | Issue | Severity | Effort | Detail |
|---|-------|----------|--------|--------|
| 1.4 | Tagged `Display` can produce unparseable EDN | Medium | Small | `Edn::Tagged("my tag", ...)` renders as `#my tag nil`. No validation on tag name. Fix: validate tag name in `Display`, or change `Tagged` to use `Symbol` instead of `String`. `display.rs:45` |
| 2.1 | `EdnMap`/`EdnSet` equality allocates and sorts every call | Medium | Small | `PartialEq` delegates to `Ord`, which allocates two `Vec`s and sorts. Every `==` on a map or set pays `O(n log n)` + alloc. Fix: implement `PartialEq` directly via hash-map lookup (each key in `self` exists with equal value in `other`, and lengths match). `collections.rs:112-136` |
| 2.2 | `format_float` allocates via `format!()` | Low | Small | `let s = format!("{v}")` heap-allocates on every float display. Fix: use `ryu` crate for stack-based float formatting. `display.rs:54` |
| 4.2 | `Tagged` uses `String` not `Cow<'a, str>` | Low | Medium | Every other string-carrying `Edn` variant uses `Cow<'a, str>`. `Tagged` always allocates the tag name, even when parsing from borrowed input. Touches `Edn` enum, parser, display, serde. `edn.rs:262` |
| 4.3 | Unsafe lifetime coercion in `EdnMap::get`/`contains_key` missing SAFETY comment | Low | Trivial | The `unsafe` block widens `&Edn<'k>` to `&Edn<'a>` via pointer cast. Sound because `Hash`/`Eq` are content-based and lifetime-independent, but this invariant is undocumented. Adding any lifetime-dependent logic to `Hash` or `Eq` would silently make it UB. `collections.rs:40-44` |

### Remaining: not worth fixing

| # | Issue | Severity | Rationale for skipping |
|---|-------|----------|------------------------|
| 1.2 | Float `+0.0` and `-0.0` are unequal | Low | Intentional `total_cmp` semantics. Documented in tests (`edn.rs:625-629`). `Eq`/`Hash` agree. Changing would break the total-ordering property needed for map keys. |
| 1.5 | No duplicate detection in parsed sets | Low | Matches Clojure behavior (silent dedup). Maps have `DuplicateKeyPolicy`; adding an equivalent for sets is possible but low priority. |
| 2.5 | `compact_len` allocates for `BigInt`/`BigDecimal` length | Low | Feature-gated, rare path. `to_string()` called speculatively to measure length, then discarded. Only matters if bignum values appear in pretty-printed output. `formatter.rs:90,93` |
| 2.6 | XOR hashing for maps/sets | Low | Vulnerable to cancellation (`hash(a) ^ hash(a) = 0`), but sets enforce uniqueness so cancellation cannot occur. Weak mixing is a theoretical concern, not a practical one for typical EDN data. `collections.rs:142-148` |
| 3.1 | Two independent parser implementations (winnow + hand-rolled streaming) | Medium | 2308 lines of parser code with only 3 shared functions. Bug fixes must be applied twice. However, the two parsers serve fundamentally different purposes (value tree vs direct serde deserialization) and share little structure. Unifying them would require an abstraction layer that may not pay for itself. |
| 3.2 | `Keyword` and `Symbol` are structurally identical | Low | ~140 lines of verbatim duplication. A generic type or macro would reduce lines but add indirection. The types are conceptually distinct in EDN and benefit from separate `Display` impls. |
| 3.3 | `EdnFormatter` fields are all `pub` | Low | No invariants to protect. Public fields are simpler than a builder for a configuration struct in research code. |
| 3.4 | `TagReaders` uses linear search | Low | `Vec<(Box<str>, TagFn)>` with `iter().find()`. O(n) per lookup. Only 2 built-in tags; becomes a concern only with many custom readers. |
| 4.1 | `Edn::iter()` silently returns empty for non-collections | Low | Convenience method. Returning `Option<Iterator>` would force `unwrap` at every call site. Current behavior is consistent with Clojure's `seq` on non-sequentials returning nil. `edn.rs:487-493` |
| 4.4 | No `FromStr` impl for `Edn` | Low | Minor ergonomic gap. `read_string()` is the canonical entry point and takes a `&str`. Adding `FromStr` is trivial but adds a second way to do the same thing. |
