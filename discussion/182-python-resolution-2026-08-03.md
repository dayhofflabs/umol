# 182 — Expose resolution on the Python surface

Status: Completed
Date: 2026-08-03
Relates: [178](178-python-lattice-ops-2026-08-01.md),
[179](179-python-editing-and-transactions-2026-08-02.md)

A molecule built by editing cannot be resolved from Python. `MoleculeAst.from_smiles` accepts
`chemistry_model` and `resolve_config` because ingest invokes resolution; `MoleculeAst.parse` accepts
only `defaults` and does not resolve; there is no `resolve` method. The three `*ResolveConfig` classes
are exported and nothing outside reaction application consumes them.

## The Rust side, for reference (refreshed 2026-08-15 after the doc 194 S4 rework)

Fully worked out; the binding has a settled shape to mirror.

```rust
Resolver::new(&chemistry_model)                          // ResolveConfig::default()
Resolver::with_config(&chemistry_model, resolve_config)  // explicit; the resolver stores it
    .resolve(&mut molecule)
    // -> Result<Solution<ResolveReport, ResolveContradiction>, ResolveError>
```

The verdict carries the report in both non-error arms: `Determined(report)` (its
`tie_breaks` name the atoms where a configured policy decided), `Underdetermined(report)`
(its `unresolved` carries the per-atom candidate lists; nothing is committed — the
journal rolls back every failure path). Ingest constructs the same resolver any caller
would; nothing about it is privileged.

`Resolver` is not re-exported from `umol-graph/src/lib.rs`, so the path today is
`umol_graph::ops::resolve::Resolver`. Worth raising when the facade of doc
[180](180-umol-facade-crate-2026-08-02.md) is built, since this is a public operation reached through
a private-looking path.

## Justification

**Resolution is never automatic** (author, 2026-08-03), and that is the design. It follows that the
explicit operation has to be callable, otherwise the only way to reach it is to route a structure
back through SMILES.

**Two sections of the whitepaper are about it.** \Cref{sec:validity} is resolution and the chemistry
model; \Cref{sec:lattice} is the order that narrowing moves along. A reader can execute the lattice
operations after 178 and the edits after 179, and cannot execute the operation those two exist to
support.

**Section 9 needs it to state its own division of labour.** Mutation performs what was written;
resolution fills what was left open, when asked. Building methylamine has two routes — state the
hydrogen count explicitly, or leave `#h*` open and resolve — and only the first can be shown in a
listing today.

## Scope

**In:**

- Resolution on `MoleculeAst`, taking `chemistry_model` and `resolve_config`, both already exported.
- The verdict. See the open question below.
- Whatever is needed for the same operation on a reaction, if that is not already reachable through
  reaction application.

**Out:**

- The phase resolvers (`ValenceResolver`, `AromaticityResolver`, `StereoResolver`, `BondsResolver`,
  `MulticenterBondsResolver`). The composite operation is the interface; the phases are internal, and
  their ordering is a live design question (doc 174).
- `ResolverError` as a value. Rollback failure indicates a defect, not an outcome a caller plans for;
  an exception is right.

## Decided: verdict value (author, 2026-08-15)

The Rust verdict returns `Solution<ResolveReport, ResolveContradiction>`, and Python had
`ContradictionError`/`UnderdeterminedError` but no `Solution`. Decision: the solution and
the contradiction are integral to the model — the three-valued outcome is the
whitepaper's validity vocabulary — and are exposed faithfully as Python values. Ingest
keeps raising (a structure that will not resolve is a failure of that call); the
explicit operation returns the verdict; only `ResolveError` — the rollback-defect
class — remains an exception. This is the difference in kind the earlier analysis
anticipated, and it follows 178's rule that ordinary outcomes are values.

## Specified surface

Two new types plus one method; nothing else. (`MoleculeAst` is `Molecule` since doc
176; the older mentions above read accordingly.)

```python
class Solution:                      # pyo3 complex enum
    Determined(molecule: Molecule, report: ResolveReport)
    Underdetermined(report: ResolveReport)
    Contradictory(contradiction: ResolveContradiction)

class ResolveContradiction:          # value wrapper over the Rust enum
    __str__   # the Rust Display message
    __repr__  # ResolveContradiction("<message>")
    __eq__    # structural, via the Rust PartialEq

Molecule.resolve(
    self, *, chemistry_model=None, resolve_config=None
) -> Solution
```

- **Non-mutating**: the receiver is never touched; `Determined` carries the resolved
  molecule (the Rust in-place mutation runs on a clone). `Underdetermined` and
  `Contradictory` carry no molecule — nothing was committed, and the receiver already
  is the unchanged input.
- **Defaults**: `chemistry_model=None` means `ChemistryModel.default()` — presets are
  reader conventions applied at format boundaries, and a constructed molecule has no
  format; `resolve_config=None` means `ResolveConfig.default()`. The different ways of
  resolving are reached through these two arguments (preset models, tie-breaks,
  failure policies).
- `Determined.report` makes tie-break uses visible on success — the exception path
  never delivered that. `Underdetermined.report` is the same `ResolveReport` already
  bound for `e.report`.
- `ResolveContradiction` starts opaque-but-printable (message, equality); per-variant
  destructuring is additive if a consumer needs it.
- Reaction resolution stays out of the minimal set; reaction application already
  resolves internally, and an explicit reaction-side operation is additive later.
- **`Determined` owns its molecule** (settled 2026-08-15): the `Molecule` wrapper gains
  `Clone` (its inner is already `Clone`), so `Solution` derives structural equality like
  every binding type — no `Py<Molecule>` handle, no hand-rolled `__eq__`. Attribute
  access clones, as it already does for `ValenceCandidateSource`'s owned fields: value
  semantics on a frozen verdict value.
- **Resolve is not validate** (confirmed 2026-08-15): admission skips ground atoms, so a
  committed model-violating structure returns `Determined` unchecked. Resolution fills
  what is open; validating committed structure under a chemistry model is a separate
  operation outside this set (its Python surface is doc 166 territory). Stated in the
  method's docstring.

## Implementation plan

One stage, additive throughout; green after every subitem.

- S0a — `umol-py::molecule`: `Clone` derive on `Molecule`. Additive. [dep: none]
  **Done 2026-08-15:** `Clone` plus the explicit `from_py_object` opt-in pyo3 now
  requires on `Clone` pyclasses. Adjacent fact, recorded for doc 181: `ReactionSpan` is
  already `Clone`; `Reaction` is deliberately not — it holds live `Py<Molecule>`/
  `Py<Deltas>` handles and has no equality, the opposite ownership design from the
  owned-value wrappers.
- S0b — `umol-py::resolve`: `ResolveContradiction` pyclass (eq, frozen; `__str__` =
  Display message, `__repr__`, `from_rust`) plus repr/str/eq test rows. Additive.
  [dep: none]
- S0c — `umol-py::resolve`: `Solution` complex enum (frozen, derived `PartialEq`;
  keyword-only variant constructors per the spec block) plus `__repr__` and tests.
  Additive. [dep: S0a, S0b]
- S0d — `umol-py::molecule`: `Molecule.resolve` per the spec (clone inner, run
  `Resolver::with_config`, map the verdict; `ResolveError` → `RuntimeError`; docstring
  states non-mutation, the default model, and the resolve-is-not-validate boundary);
  `umol-py` gains the `umol-utils` dependency to name the Rust `Solution`. Rust-side
  tests: the three verdicts from constructed molecules, non-mutation of the receiver in
  all three arms, the default-model case. [dep: S0c]
- S0e — registration and Python tests: `lib.rs` exports and `add_class` for the two
  types, `__init__.py`, the import inventory; pytest rows for the three verdicts and
  non-mutation, including the `match`/`case` idiom — this is where pyo3's
  structural-pattern support is verified, and the documented idiom falls back to class
  patterns with attribute access if keyword patterns are unsupported. [dep: S0c, S0d]

Critical path S0a → S0c → S0d → S0e; S0b is parallel to S0a.

**S0b–S0e done 2026-08-15.** `ResolveContradiction` and `Solution` landed as specced
(frozen, derived equality; keyword-only constructors); `Molecule.resolve` maps the Rust
outcome with `ResolveError` → `RuntimeError`; `umol-py` depends on `umol-utils` directly
per the no-re-export policy. Rust-side tests pin the three arms, non-mutation of the
receiver in every arm, and the default-model case (a charge-open atom takes the
registry's charge-less lookup: nine carbon candidates); the pytest rows exercise the
same through `match`/`case` — pyo3 complex-enum keyword patterns verified working, so
the documented idiom stands without fallback. Vocabulary ruling applied while landing:
"verdict" does not appear in `umol-py` or `umol-graph` — the word is *solution*
(rustdoc, locals, and the `validate` helper renamed). Suites: pytest 1313, umol-py
1631, workspace clippy clean.

## Settled semantics

- Rust mutates the AST in place. Python should follow 178's precedent and **return the resolved
  molecule without mutating the receiver**, consistent with `canonicalize`.
- `chemistry_model` and `resolve_config` are optional and default as they do for `from_smiles`.
- Resolution is implemented as an edit transaction with a rollback journal, so a contradiction leaves
  the input unchanged. That property should hold at the Python boundary and be tested.

## Verification

Follow 178 and 179: algebraic properties stay in Rust, and the Python tests check availability and
representative cross-boundary results. At minimum, the three verdicts each reached from a constructed
molecule, and the rollback property — a contradictory resolution leaves the input untouched.

The reader-facing check is that \Cref{sec:validity} and \Cref{sec:mutation} listings execute against
the built module before those sections ship.

## Note

Naming is unaffected by doc [176](176-ast-naming-2026-07-31.md); only the class this hangs off would
move.
