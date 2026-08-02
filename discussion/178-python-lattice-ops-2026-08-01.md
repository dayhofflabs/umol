# 178 — Expose the lattice operations on the Python surface

Status: Proposed
Date: 2026-08-01
Relates: [113](113-ast-canonical-equality-and-lattice-2026-06-14.md),
[177](177-nomenclature-guide-2026-07-31.md)

The `Lattice` and `Canonicalize` operations are Rust-only. Python callers can construct and inspect
AST values but cannot ask any of the four questions the algebra answers. This makes the central idea
of the model unreachable from the interface most users have, and it is the one part of the whitepaper
whose printed claims cannot be reproduced by a reader.

## Scope

**Every exported `Lattice` type, uniformly.** Every one is already exported to `umol-py`
(`add_class::<…>`), so no new classes are needed — only methods on existing ones. Verified
2026-08-01: the Rust property suite exercises 40 unique `Lattice` types, and all 40 are exported.
The trait bound, not this snapshot count, defines the scope.

Uniformity is the point. A scoped subset would leave callers to discover by trial which types happen
to answer `meet`, and "why does `AtomAst` have it but `IsotopeMassAst` not" has no good answer. The
per-type cost after the macro is written is close to nil.

`Canonicalize` has a broader scope. Its two operations belong on every exported type implementing
`Canonicalize`, including types that do not implement `Lattice` such as the molecule-constraint
expressions, `Deltas`, and `ReactionAst`.

## Methods

From `Lattice`:

| Rust | Python | Returns |
| --- | --- | --- |
| `is_undetermined` | `is_undetermined()` | `bool` |
| `is_ground` | `is_ground()` | `bool` |
| `meet` | `meet(other)` | the value, or `None` when no ground value satisfies both |
| `join` | `join(other)` | the value, or `None` when no join exists |
| `matches` | `matches(target)` | `bool` |
| `is_compatible` | `is_compatible(other)` | `bool` |

From `Canonicalize`, worth including in the same pass since they complete the picture and
`canonical_eq` is otherwise unreachable:

| Rust | Python | Returns |
| --- | --- | --- |
| `canonicalize` | `canonicalize()` | the canonical value; raises on an unsatisfiable value |
| `canonical_eq` | `canonical_eq(other)` | `bool` |

## Settled semantics

- `meet` and `join` both return `None` when the requested bound does not exist. Failure to relate two
  otherwise valid values is an ordinary result, not a Python exception. No `NoJoinError` is added.
- `canonicalize` returns the canonical value without mutating its receiver and raises the existing
  `ContradictionError` when the receiver is unsatisfiable. There is no useful result value in that
  case.
- `canonical_eq` follows Rust exactly, including treating two contradictory values as canonically
  equal. It does not raise merely because either operand is contradictory.
- `narrow_from` and `widen_with` are not exposed in this pass. Adding them later is not blocked by the
  non-mutating surface chosen here.
- `__eq__` stays structural. Semantic equality remains an explicit `canonical_eq` call, preserving
  the distinction in doc 113 and keeping Rust and Python aligned.
- Binary methods accept another instance of the same Python wrapper. PyO3 reports a `TypeError` for
  a different wrapper type.
- `meet`, `join`, and `canonicalize` return values and never mutate either input. Python object
  identity is not part of their contract, even when the result is structurally unchanged.
- Rust's borrowed `canonical()` fast path is not exposed. Python cannot preserve its `Cow` semantics,
  and `canonicalize()` is the complete value-level operation callers need.

## Verification

- Keep the algebraic property tests in `umol-ast`, where they already cover the Rust implementations.
  The binding tests instead verify complete method availability and representative cross-boundary
  results for bounded, fibered, contradictory, leaf, container, and entity AST values. This tests
  the binding without introducing a second property-testing stack or duplicating the Rust suite.
- The whitepaper's `tab:lattice-ops` (Section 8) was computed in Rust and is currently unreproducible
  by a reader. Once this lands, rerun it through Python and replace the table with a listing, putting
  Section 8 on the same footing as Sections 3–7 where everything printed executes.

## Staged implementation plan

### S0 — Binding kernel

- **S0a — Generate the uniform method surface.** **Done.** In `umol-py/Cargo.toml`, enable PyO3's
  `multiple-pymethods` support. Add `umol-py/src/lattice.rs` with complete macros for the
  `Canonicalize` pair and the six non-mutating `Lattice` methods; the invocations supply the existing
  per-wrapper `from_rust` and `to_rust` expressions rather than introducing another conversion trait.
  Wire the module from `lib.rs` and apply both surfaces to `BooleanAst` as the direct-conversion case.
  Add exact Python assertions for top, ground, compatible, incompatible, meet, join, matching,
  canonicalization, and canonical equality. **Additive (green).** [dep: none]
- **S0b — Cover Python-dependent conversions and contradictions.** **Done.** Apply the generated
  surfaces to `ValueAst` and the canonicalization-only surface to `ValueTerm`, exercising wrappers
  whose conversions require `Python` and can fail. Test structural-versus-canonical equality, `None` from
  absent bounds, `ContradictionError` from `canonicalize`, and Rust's
  two-contradictions-compare-equal rule.
  **Additive (green).** [dep: S0a]

### S1 — Leaf AST values

- **S1a — Atom and electron leaves.** **Done.** Apply the surfaces in `atom.rs`, `spin.rs`,
  `electrons.rs`, and `constraint/atom.rs` to `ElementAst`, `IsotopeMassAst`,
  `UnpairedElectronsAst`, `ElectronCountsAst`, `AromaticValenceAst`, and `MulticenterValenceAst`.
  Add table-driven Python
  assertions spanning undetermined, literal, set/range, compatible, and incompatible cases.
  **Additive (green).** [dep: S0b]
- **S1b — Stereo leaves.** **Done.** Apply the surfaces in `stereo.rs` and `constraint/stereo.rs` to
  `StereoConfigurationAst`, `TetrahedralStereoAst`, `CisTransStereoAst`, `StereogenicityAst`,
  `TopicityRelationAst`, `LigandSymmetryAst`, `FluxionalityAst`, and `TopicityAst`. Include a fibered
  no-join case so `join` returning `None` is verified independently of an incompatible meet.
  **Additive (green).** [dep: S0b]
- **S1c — Remaining leaves.** **Done.** Apply the surfaces in `noncovalent.rs` and `constraint/ring.rs` to
  `NoncovalentBondKindAst` and `RingMembershipAst`. Verify both operations through their public
  constructors and exact expected wrapper values. **Additive (green).** [dep: S0b]

### S2 — Entity and entity-constraint AST values

- **S2a — Atom and bond entities.** **Done.** Apply the lattice and canonicalization surfaces to `AtomAst` and
  `BondAst`. Test full-value results rather than reaching into individual fields, including one
  compatible refinement and one incompatible pair for each entity. **Additive (green).** [dep: S1a]
- **S2b — Overlay entities.** **Done.** Apply the surfaces to `DativeBondAst`, `AromaticSystemAst`,
  `MulticenterBondAst`, `NoncovalentBondAst`, `StereoAtomAst`, and `StereoBondAst` in their entity
  modules. The stereo entity wrappers are macro-generated and remain part of the inventory. Verify
  composite meet, join, matching direction, and canonical equality with exact expected objects.
  **Additive (green).** [dep: S1a, S1c]
- **S2c — Atom and bond constraint families.** **Done.** Apply the surfaces to the enum and
  container pairs in `constraint/atom.rs` and `constraint/bond.rs`: `AtomConstraintAst`,
  `AtomConstraintsAst`, `BondConstraintAst`, and `BondConstraintsAst`. Test same-kind bounds,
  different-fiber `join`, container canonicalization, and `matches` direction.
  **Additive (green).** [dep: S1a]
- **S2d — Overlay constraint families.** **Done.** Apply the surfaces to the enum and container
  pairs in `constraint/aromatic.rs`, `constraint/dative.rs`, `constraint/multicenter.rs`, and
  `constraint/noncovalent.rs`. Use table-driven Python cases covering every family and exact outputs
  for both enum and container operations. **Additive (green).** [dep: S1c, S2b]
- **S2e — Stereo constraint families.** Apply the surfaces to `StereoAtomConstraintAst`,
  `StereoAtomConstraintsAst`, `StereoBondConstraintAst`, and `StereoBondConstraintsAst`. Verify
  same-fiber and cross-fiber behavior and canonicalization of the corresponding containers.
  **Additive (green).** [dep: S1b]

### S3 — Canonicalization-only types and completeness

- **S3a — Molecule-constraint expressions.** Apply the canonicalization surface in
  `constraint/molecule.rs` to every exported canonicalizable non-lattice type there, including
  `RelationalConstraint`, `MoleculeConstraint`, `Constraint`, and `Constraints`. Test recursive
  canonicalization, contradiction propagation, structural inequality before canonicalization, and
  canonical equality afterward. **Additive (green).** [dep: S1a, S1b, S2c, S2d, S2e]
- **S3b — Delta and reaction containers.** Move the existing `Deltas.canonicalize` and
  `ReactionAst.canonicalize` bindings onto the uniform canonicalization implementation and add
  `canonical_eq` without changing their existing exception behavior. Test detached results and
  canonical equality for reordered deltas; for reactions, preserve the fixed LHS id-space semantics.
  **Additive (green).** [dep: S3a]
- **S3c — Close the exported-type inventory.** Audit every class registered by `lib.rs` against the
  Rust `Lattice` and `Canonicalize` implementations, add any canonicalizable non-lattice wrapper not
  covered by S3a/S3b, and add one explicit Python inventory test asserting the complete method set on
  every applicable exported class and its absence from unrelated classes. **Additive (green).**
  [dep: S2a, S2b, S2c, S2d, S2e, S3a, S3b]

### S4 — Reader-facing verification

- **S4a — Reproduce the whitepaper lattice example through Python.** Express the inputs and
  operations behind `tab:lattice-ops` solely with the public Python API, assert the printed values in
  the Python suite, and replace the static whitepaper table with the executable listing and its
  output. Run formatting, clippy, the `umol-py` Rust tests, rebuild the Python extension under the
  Python 3.13 virtual environment, and run the complete Python suite. **Additive (green).**
  [dep: S3c]

**Critical path:** S0a → S0b → S1a/S1b/S1c → S2 → S3 → S4a.

No implementation stage is deferrable: S3c establishes the promised uniform surface, and S4a closes
the reader-reproducibility requirement that motivated the work.

## Notes

- Method names are unaffected by the naming question in doc 176; only the classes they hang off would
  move. Landing this first is safe, at the cost of touching the binding surface twice.
- If mutating lattice operations are added later, `narrow_from` and `widen_with` remain the
  plain-English names for descending and ascending the order; see the *Narrow and widen* entry in doc
  177 before renaming anything.
