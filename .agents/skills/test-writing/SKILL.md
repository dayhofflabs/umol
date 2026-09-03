---
name: test-writing
description: MANDATORY — use whenever creating, editing, reviewing, renaming, or restructuring tests, cases, fixtures, or test assertions in a umol workspace crate, including test work inside a larger task. Enforces rstest structure, naming, construction, assertions, independence, ordering, and crate-specific patterns.
---

# Test writing

For property-test work, read `docs/development/property-tests.md` completely first.

## Scope

- Every touched test follows this skill. Do not copy a nonconforming neighboring pattern.
- Do not silently rewrite unrelated nonconforming tests; ask first. Reordering tests to match source
  definitions is the one permitted incidental cleanup.

## Structure

- Use `#[rstest]`, never bare `#[test]`.
- Use `#[fixture]` for shared inputs without per-case variation.
- Use `#[case]` table rows for multiple cases and for a single case with likely variation.
- Put literals directly in case rows; do not add test constructor helpers. Use public construction
  unless that makes the test circular; otherwise use struct literals. Prefer direct collection
  conversions (`From`, `From<Vec<_>>`, or `IntoIterator`).

## Names and assertions

- Free function: `test_<function>`; method: `test_<struct>_<method>`.
- An optional final noun may group a scenario, such as `_error`, `_identity`, `_roundtrip`, or
  `_partial`. Never encode behavior (`returns_zero`, `drops_x`) in the function name; put it in the
  `#[case]` label. Add a scenario suffix only when it simplifies assertions or separates a large
  surface.
- Assert exact values or error variants. A length, `is_some`, or error presence is insufficient by
  itself. For transforms returning `Self`, prefer full structural equality.
- Split unchanged transform cases into `test_X_method_identity`; each row supplies only the input and
  asserts `input.clone().method() == input`. Keep `(input, expected)` rows in `test_X_method`.

## Independence

Do not construct a test through behavior on which the subject depends. In particular:

- Test `is_ground` with struct literals, not `into_ground`.
- For `into_ground`, construct the expected value with a struct literal.
- Comparing `into_zeroed` with the corresponding macro is allowed as cross-path equivalence.

Prefer a verbose independent expected value to a concise circular one.

## Order

Parallel source definitions: keep each type's tests contiguous and types in source order. Within a
type, follow declaration order (`new`/`from_*`, `with_*`, transforms, queries, mutators); place an
error group after its positive group. If a `From` impl precedes the inherent impl, its tests precede
the method tests.

## Prohibitions

- No region markers, divider comments, helper constructors, smoke tests, or behavior in test names.
- No tautological constructor tests when `From` makes the constructor redundant.

## Crate patterns

### `umol-graph-ir::macros`

- One positive table per macro for distinct shapes; one separate
  `test_<macro>_macro_error` with `#[should_panic]` for fallible macros.
- For a macro wrapping `FromStr` plus config, compare with an explicit
  `from_atoms_and_bonds`/`with_*`/`into_ground` construction.

### `umol-graph-ir::ir::*`

- Order methods as `new`/`from_*`, `with_*`, `into_*`, `is_ground`, `matches`, `simplify_*`.
- Consolidate `with_*` methods in one `test_<entity>_with_methods` table.
- Test `into_ground` and `into_zeroed` separately; keep preservation and empty-input cases together.

### `umol-graph-ir::dsl::config`

- Follow definition order: `MoleculeDefaults`, `AtomDefaults`, `BondDefaults`, and so on.
- Do not test no-op defaults or one-line delegations already covered by `with_overrides`.
