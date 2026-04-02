# umol-edn Code Review (Adversarial)

Date: 2026-03-26

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

## Severity Legend

- `P0`: correctness/safety contract breach; can produce wrong behavior now.
- `P1`: spec/API semantic break or major reliability risk.
- `P2`: significant inefficiency/complexity or sharp-edge behavior.
- `P3`: maintainability debt, API roughness, or lower-impact issues.

## Findings

### P0-1: `Eq` and `Hash` are inconsistent for `Edn::Float`

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

### P0-2: Tagged unit variant deserialization discards payload silently

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

### P0-3: `bignum` is claimed but not implemented

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

### P1-1: Spec contradiction for built-in unqualified tags without features

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

### P1-2: Spec says BTree deterministic ordering; code uses hash collections

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

### P1-3: Macro/parser behavior diverges on special floats and char names

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

### P2-1: Struct deserialization silently ignores non-string-like map keys

**Evidence**

- `umol-edn/src/de.rs`: `EdnStructMapAccess` skips keys not in `Keyword|Symbol|Str`.
- `umol-edn/src/de.rs` test `test_deserialize_struct_ignores_non_string_keys` codifies behavior.

**Why this is risky**

Malformed config input may deserialize successfully with fields omitted, instead of failing loudly.

**Fix direction**

- Prefer strict mode by default (error on non-struct keys) and optionally offer lenient mode.

---

### P2-2: Avoidable allocations in `deserialize_any`

**Evidence**

- `umol-edn/src/de.rs`: keywords/symbols converted to owned `String` with `.to_string()`.
- `umol-edn/src/de.rs`: sets copied into `EdnSeq` via `collect` before visiting.

**Impact**

- Higher allocation cost for common paths.

**Fix direction**

- Use borrowed visitor paths when possible.
- Introduce dedicated set `SeqAccess` adapter to avoid intermediate `EdnSeq`.

---

### P2-3: `Edn::get(&str)` is keyword-only and allocates every call

**Evidence**

- `umol-edn/src/edn.rs`: `get(&str)` looks up only `Edn::Keyword(Keyword::owned(key.to_string()))`.

**Risk**

- Misleading API name (`get`) for a narrow key type.
- Repeated allocations in map-heavy reads.

**Fix direction**

- Rename to `get_keyword` (or equivalent), and add broader lookup API by `&Edn`.
- Consider borrowed lookup strategy.

---

### P2-4: Collection `Ord` implementations are expensive

**Evidence**

- `umol-edn/src/collections.rs`: `EdnMap::cmp` and `EdnSet::cmp` allocate vectors and sort on every comparison.

**Impact**

- `O(n log n)` with allocation per compare; poor behavior when sorting/comparing many nested EDN values.

**Fix direction**

- Document as intentionally expensive if only used for deterministic ordering logic.
- Revisit ordering strategy if `Ord` is required in hot paths.

---

### P2-5: Panics in serde map accessors (`expect(...)`)

**Evidence**

- `umol-edn/src/de.rs`: `next_value_seed` in `EdnMapAccess` and `EdnStructMapAccess` uses `expect(...)`.

**Risk**

- Panics if serde contract is violated by custom seeds/visitors.

**Fix direction**

- Replace `expect` with explicit `EdnError::Custom(...)` return.

---

### P3-1: `unwrap_err` loses offset context on incomplete input

**Evidence**

- `umol-edn/src/error.rs`: `ErrMode::Incomplete(_) => UnexpectedEof { offset: 0 }`.

**Impact**

- Poor diagnostics at exact failure location.

**Fix direction**

- Preserve best-known offset when mapping incomplete errors.

---

### P3-2: `unreachable!()` in `Edn::cmp` is a maintenance landmine

**Evidence**

- `umol-edn/src/edn.rs`: wildcard arm in `match (self, other)` after discriminant compare.

**Risk**

- Future variant changes can introduce panic if `variant_ord` and `cmp` updates drift.

**Fix direction**

- Keep matches structurally exhaustive without wildcard fallback; avoid runtime trap.

---

### P3-3: Documentation drift in README feature table

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
