---
description: Use this skill whenever the task involves adding, modifying, restructuring, or reviewing unit tests inside any `#[cfg(test)] mod tests` block in a umol-workspace crate. Covers rstest table-test conventions, fixture usage, naming rules, ordering rules (parallel to module definitions), assertion style, identity-comparison patterns, circular-logic avoidance, and crate-specific patterns for macros, AST entities, and defaults configs. Apply before writing the first test in a session and consult on every subsequent test edit.
---

# umol test-writing conventions

## Framework

- `#[rstest]` on every test. No bare `#[test]`.
- `#[fixture]` for shared inputs that take no per-case variation.
- Table tests via `#[case]` rows for any function with more than one case.
- A single-case test is still written as `#[rstest]` with one `#[case]` row when the function has variation potential — keeps the surface uniform and makes it trivial to extend.

## Construction in cases

- Inline literal values directly in `#[case]` rows. Do not introduce helper constructor functions.
- Use the type's public API where it doesn't introduce circular logic (see below). Otherwise use struct literals.
- For collection containers with `From<T>` / `From<Vec<T>>` / `IntoIterator`, prefer the most direct path.

## Naming

- `test_<function>()` for free functions.
- `test_<struct>_<method>()` for inherent methods.
- `test_<struct>_<method>_error()` for the failure path of the same method.
- **Do not embed test behavior in the name.** No `test_X_returns_zero_when_empty`, no `test_X_handles_negative`. Behavior belongs in `#[case]` labels.
- Acceptable qualifiers: `_error` (failure path), `_identity` (when split out — see Identity comparisons), `_partial` (when distinguishing related tests on the same method).

## Assertions

- Assert specific return values or error variants. Avoid summary-stat-only assertions (`len()`, `is_some()`, presence-of-error) as the sole check.
- For methods returning `Self` (builder, transform), prefer full structural equality (`assert_eq!`) over per-field assertions when the expected can be constructed cleanly.

## Identity comparisons

When a transform method returns the input unchanged for some input shapes, split into two table tests:

- `test_X_method` — the transforming cases. `(input, expected)` pairs.
- `test_X_method_identity` — the identity cases. Single `#[case]` per row; assert `input.clone().method() == input`.

This avoids duplicating values in `#[case]` rows where input == expected.

## Avoid circular logic

A test for method `M` should not construct its input via method `N` whose correctness `M` relies on. Concretely:

- `test_X_is_ground` should construct inputs via struct literals, not via `X::into_ground()`. If `into_ground` has a bug that leaves a field `Undetermined`, the per-field test cases would fail ambiguously rather than precisely identifying the broken field.
- `test_X_into_ground` constructs the input however convenient, but the *expected* must be a struct literal (not built via `into_ground`).
- `test_X_into_zeroed` against `X_zeroed!()` macro output is a cross-path equivalence test, not circular — it verifies method-vs-macro convergence.

When the expected for a complex constructor is too verbose, prefer struct literals over delegating to other API methods.

## Ordering

**Tests within a `mod tests` block are ordered to parallel the module's type and method definitions.**

- Group by struct/enum: each type's tests are contiguous, in the order types appear in the source.
- Within a struct's tests, methods in declaration order: `new` / `from_*` first, then `with_*`, then transforms (`into_*`), then queries (`is_*`, `matches`), then mutators (`simplify_*`).
- `_error` variants follow their corresponding positive test.
- `From` impl tests precede the type's other method tests when the impl is declared before the impl block.

**Rearranging existing tests to fit this order is allowed and explicitly deviates from the "don't touch unrelated code" rule** — it's a permitted incidental cleanup when adding or modifying tests.

## What not to do

- No `// region:` / `// endregion:` markers or section divider comments.
- No helper constructor functions (`fn ground_atom() -> AtomAst { ... }`). Inline the construction in each `#[case]`.
- No "smoke tests" — either write a real test with specific assertions or do not test.
- No tautological tests (e.g., `assert_eq!(X::new(v), X::Variant(v))` with no other behavior to verify) when adding a `From` impl makes the constructor redundant — drop the `new`, drop the test.
- No behavior in test names. Use `#[case::descriptive_label]` rows instead.

## Crate-specific patterns

### Macros (`umol-ast::macros`)

- Each macro gets one positive test with `#[case]` rows for distinct shapes (e.g., empty + populated).
- Each fallible macro gets a separate `test_<macro>_macro_error` with `#[should_panic]`.
- For macros that wrap `FromStr` plus a config (e.g., `mol_ground!`), test against an explicitly-constructed expected via `from_atoms_and_bonds` + `with_*` chains + `into_ground()`. This verifies the macro applies the correct config.

### AST entity tests (`umol-ast::ast::*`)

- Method order: `new`/`from_*` → `with_*` → `into_*` → `is_ground` → `matches` → `simplify_*`.
- `with_*` methods consolidate into one `test_<entity>_with_methods` table.
- `into_ground` and `into_zeroed` get separate tests; preserve-existing-value cases belong in the same test as the from-empty case.

### Defaults configs (`umol-ast::dsl::config`)

- Test order: `MoleculeDefaults` first (composes), then `AtomDefaults`, then `BondDefaults`, etc., matching definition order.
- No-op default impls (e.g., `DativeBondDefaults`, `NoncovalentBondDefaults`) need no tests.
- One-line delegations (e.g., `BondDefaults::ground()` → `Self::zeroed()`) need no separate tests; `with_overrides` tests cover the surface.
