# 183 — Parse and render entity updates in Python

Status: Completed
Date: 2026-08-03
Relates: [179](179-python-editing-and-transactions-2026-08-02.md)

Doc 179 exposed the entity `*Update` types needed to build `Edits`, but those values can currently be
constructed only through their field-based constructors. `*Update` is the all-optional analog of
`*Ast`, so its textual API should mirror the `*Ast` API directly. The required syntax is already
implemented by the Rust `*UpdateDsl` types.

## Verified current shape

The gap applies uniformly to all eight entity kinds.

| Entity AST | Entity DSL | Entity update | Update DSL |
| --- | --- | --- | --- |
| `AtomAst` | `AtomDsl` | `AtomUpdate` | `AtomUpdateDsl` |
| `BondAst` | `BondDsl` | `BondUpdate` | `BondUpdateDsl` |
| `DativeBondAst` | `DativeBondDsl` | `DativeBondUpdate` | `DativeBondUpdateDsl` |
| `AromaticSystemAst` | `AromaticSystemDsl` | `AromaticSystemUpdate` | `AromaticSystemUpdateDsl` |
| `MulticenterBondAst` | `MulticenterBondDsl` | `MulticenterBondUpdate` | `MulticenterBondUpdateDsl` |
| `NoncovalentBondAst` | `NoncovalentBondDsl` | `NoncovalentBondUpdate` | `NoncovalentBondUpdateDsl` |
| `StereoAtomAst` | `StereoAtomDsl` | `StereoAtomUpdate` | `StereoAtomUpdateDsl` |
| `StereoBondAst` | `StereoBondDsl` | `StereoBondUpdate` | `StereoBondUpdateDsl` |

Each `*UpdateDsl` is already a public transparent wrapper around the corresponding `*Update`. Each
already implements `FromStr`, `Display`, `FromEdn`, and `ToEdn`, and the public `parse_*_update`
function already parses its complete string subgrammar. The update syntax therefore does not need to
be designed or duplicated.

The existing `*Ast` pattern is:

```text
string -> *Dsl::from_str -> *Dsl::into_ast -> *Ast
*Ast -> *Dsl -> Display -> string
```

The `*Update` pattern should be identical:

```text
string -> *UpdateDsl::from_str -> *UpdateDsl::into_ast -> *Update
*Update -> *UpdateDsl -> Display -> string
```

This requires `FromAst<*Update>` and `IntoAst<*Update>` on every `*UpdateDsl`, followed by direct
`FromStr` and `Display` implementations on every `*Update`. The conversion context is `()` because
updates have no defaults: omission means “leave unchanged.” Conversion is otherwise only wrapping and
unwrapping the transparent DSL type.

The update types are semantically all-optional even though their constraint fields are containers
rather than `Option<ConstraintsAst>`. An empty or vacuous constraint container means that no
constraint update is present; individual non-vacuous entries express constraint changes or removals.

## Python surface

Each of the eight Python `*Update` classes should gain the same textual surface currently exposed by
the corresponding entity `*Ast` class:

```python
update = AtomUpdate.parse("#c-1#h*")
assert str(update) == "#c-#h*"
assert repr(update) == "AtomUpdate.parse('#c-#h*')"
```

`parse()` delegates to the direct Rust `FromStr` implementation and maps `ParseError` through the
existing Python error conversion. `str()` delegates to Rust `Display`. `repr()` remains evaluable and
uses the canonical rendered form, as the entity `*Ast` wrappers already do.

A separate `render()` method is not needed. This exactly mirrors the individual Python entity `*Ast`
classes, which expose `parse()`, `str()`, and `repr()`.

Parsing performs syntactic conversion and structural-integrity checks only. It does not validate the
chemical semantics of the resulting update.

The scope is the eight entity updates listed above. Component updates such as
`UnpairedElectronsUpdate` and `StereoConfigurationUpdate` have no corresponding entity `*Ast` /
`*Dsl` pair and are outside this parity change.

## Verification requirements

- For every Rust `*Update`, establish direct parse/render round trips through `FromStr` and `Display`,
  including the empty update, an omitted field, an explicit undetermined field, and a constraint
  removal where the grammar supports it.
- Establish that direct parsing reports the same `ParseError` as the existing `*UpdateDsl` parser.
- For every Python `*Update`, exercise `parse()`, `str()`, evaluable `repr()`, equality after
  round-trip, and representative invalid input.
- Keep the existing `*UpdateDsl` syntax tests as the detailed grammar tests; the new tests verify the
  direct Rust and Python surfaces rather than duplicating every grammar case.

## Staged implementation plan

All changes are additive. The tree remains green after every subitem and stage.

### S0 — Direct Rust textual APIs

Each subitem applies the same complete pattern in its entity DSL module: add the zero-copy
`*UpdateDsl::from_ref` projection, implement `FromAst<*Update>` and `IntoAst<*Update>` with `Ctx = ()`,
then implement `FromStr` and `Display` on `*Update`. Keep the detailed syntax cases on `*UpdateDsl`;
add a direct `*Update` display/from-string property using the existing update strategy and a focused
direct parse-error table.

- **S0a — Atom updates.** **Done.** In `umol-ast/src/dsl/atom.rs` and
  `umol-ast/tests/property/entity.rs`, add the direct textual API and tests for `AtomUpdate`.
  **Additive (green).** [dep: none]
- **S0b — Bond updates.** **Done.** In `umol-ast/src/dsl/bond.rs` and
  `umol-ast/tests/property/entity.rs`, add the direct textual API and tests for `BondUpdate`.
  **Additive (green).** [dep: none]
- **S0c — Dative-bond updates.** **Done.** In `umol-ast/src/dsl/dative.rs` and
  `umol-ast/tests/property/entity.rs`, add the direct textual API and tests for `DativeBondUpdate`.
  **Additive (green).** [dep: none]
- **S0d — Aromatic-system updates.** **Done.** In `umol-ast/src/dsl/aromatic.rs` and
  `umol-ast/tests/property/entity.rs`, add the direct textual API and tests for
  `AromaticSystemUpdate`. **Additive (green).** [dep: none]
- **S0e — Multicenter-bond updates.** **Done.** In `umol-ast/src/dsl/multicenter.rs` and
  `umol-ast/tests/property/entity.rs`, add the direct textual API and tests for
  `MulticenterBondUpdate`. **Additive (green).** [dep: none]
- **S0f — Noncovalent-bond updates.** **Done.** In `umol-ast/src/dsl/noncovalent.rs` and
  `umol-ast/tests/property/entity.rs`, add the direct textual API and tests for
  `NoncovalentBondUpdate`. **Additive (green).** [dep: none]
- **S0g — Stereo-entity updates.** **Done.** In `umol-ast/src/dsl/stereo.rs` and
  `umol-ast/tests/property/stereo/serialization.rs`, add the direct textual APIs and tests for the
  paired `StereoAtomUpdate` and `StereoBondUpdate` types. Keep the pair together because both
  implementations and their shared configuration grammar live in the same module. **Additive
  (green).** [dep: none]
- **S0h — Update-DSL consumer sweep.** **Done.** Migrate the reaction DSL streaming and tree
  readers, standalone edit DSL readers, and atom/bond update macros from direct tuple-field
  extraction to the public `IntoAst` or direct `FromStr` APIs. Audit all eight update wrappers so
  that tuple-field access remains only within their wrapper implementations. **Additive (green).**
  [dep: S0a, S0b, S0c, S0d, S0e, S0f, S0g]

Stage gate: the ordinary `umol-ast` suite and the generated entity-serialization properties pass.

### S1 — Mirrored Python APIs

For each entity module, add `*Update.parse()`, `__str__`, and evaluable `__repr__` methods that delegate
to the direct Rust API from S0. Extend the existing per-entity Python test module with table cases for
the empty update, a populated or explicitly undetermined field, a constraint removal, canonical
string output, evaluable representation, round-trip equality, and representative invalid input.

- **S1a — Atom updates.** **Done.** Update `umol-py/src/atom.rs` and `umol-py/tests/test_atom.py`.
  **Additive (green).** [dep: S0a]
- **S1b — Bond updates.** **Done.** Update `umol-py/src/bond.rs` and
  `umol-py/tests/test_bond.py`.
  **Additive (green).** [dep: S0b]
- **S1c — Dative-bond updates.** **Done.** Update `umol-py/src/dative.rs` and
  `umol-py/tests/test_dative.py`. **Additive (green).** [dep: S0c]
- **S1d — Aromatic-system updates.** **Done.** Update `umol-py/src/aromatic.rs` and
  `umol-py/tests/test_aromatic.py`. **Additive (green).** [dep: S0d]
- **S1e — Multicenter-bond updates.** **Done.** Update `umol-py/src/multicenter.rs` and
  `umol-py/tests/test_multicenter.py`. **Additive (green).** [dep: S0e]
- **S1f — Noncovalent-bond updates.** **Done.** Update `umol-py/src/noncovalent.rs` and
  `umol-py/tests/test_noncovalent.py`. **Additive (green).** [dep: S0f]
- **S1g — Stereo-entity updates.** **Done.** Update `umol-py/src/stereo.rs` and
  `umol-py/tests/test_stereo.py` for both `StereoAtomUpdate` and `StereoBondUpdate`. Keep the pair
  together because their bindings share the same module and method structure. **Additive (green).**
  [dep: S0g]

Stage gate: rebuild the extension in the Python 3.13 `umol-py/.venv`, then run all seven affected
Python test modules successfully.

### S2 — Integrated verification and closeout

- **S2a — Full verification.** **Done.** Format the workspace; run the complete `umol-ast` unit
  and property suites with the `proptest` feature, clippy for the affected Rust crates and targets
  with warnings denied, rebuild `umol-py` in its Python 3.13 virtual environment, and run the
  complete Python test suite. **Additive (green).** [dep: S0a, S0b, S0c, S0d, S0e, S0f, S0g, S1a,
  S1b, S1c, S1d, S1e, S1f, S1g]
- **S2b — Closeout.** **Done.** Mark this document `Completed` and update `000-status.md` after S2a
  passes. **Additive (green).** [dep: S2a]

The critical path is the Rust API for an entity update, its corresponding Python methods, and full
verification: `S0x -> S1x -> S2a -> S2b`. The S0 subitems are mutually independent, as are the S1
subitems once their matching Rust dependency is complete. No stage is deferrable: the work unit is
the uniform eight-entity API, not a partial binding for selected entity kinds.
