# 137 — Python bindings and documentation strategy

Design entry for the Python wrapper and for the public-facing documentation that
grows alongside it. Records the settled infrastructure choices, the one forced
design axis (an owned facade over a borrowed/generic core), and the open decisions.
Documentation strategy is folded in here because the wrapper is the first artifact
written for external consumption and the two efforts share a surface.

## Scope

Two coupled deliverables:

1. A Python package over the existing Rust crates: functions to functions, types to
   types, with the Python surface used as pressure to refine the Rust APIs.
2. The start of a public **API Guide** (and later a **User Guide** with Getting
   Started), authored for the whitepaper's supporting information — a conceptual
   overview and background, not a symbol reference.

## Verified current state

- **Two disjoint domain roots.** `umol-graph` and `umol-geometric` share no
  dependency edge; a consumer of `umol-graph` compiles zero *msym* code. Domain
  isolation already exists at the crate level (doc 129).
- **Both domains carry a `-sys` C build** (corrected at S0a, 2026-07-05 — the
  earlier "graph = no C build" premise from doc 129 was wrong):
  - *graph* — `umol-graph-core` depends on `nauty-Traces-sys`, which builds via `cc`
    **and `bindgen` (libclang/LLVM)**. `umol-ast` depends on `umol-graph-core`, so
    *every* binding over the AST pulls it. Self-contained (vendored; no system lib or
    submodule), but needs a C + libclang toolchain to build from source.
  - *geometric* — `umol-msym-sys` (`cc` compiling vendored `libmsym/src/*.c` + a git
    submodule), reachable only through `umol-geometric`; no bindgen.
- **No umbrella crate.** Deferred in doc 129; the per-crate split already isolates
  the C build for Rust consumers, so no single entry point existed to hang features
  on.
- **Core type shapes** (the fact that drives the binding design):

  | Shape | Types | Can be a `#[pyclass]`? |
  | --- | --- | --- |
  | Owned, non-generic | `MoleculeAst`, `ReactionSpanAst` | Yes, directly |
  | Borrowed `<'a>` | `ReactionDerivation<'a>`, and the whole read API: `AtomView`, `BondView`, `GraphView`, `RingView`, `AromaticSystemView`, `DativeBondView`, `MulticenterBondView`, `NeighborView`, `StereoLigandView`, … | No |
  | Generic | `EntitySpan<T>`, `GraphSymmetryConfig<C>` | No |

  PyO3 `#[pyclass]` requires `'static`, non-generic, owned types. The owned AST
  roots cross directly; the borrowed view layer and the generic types do not.

## Binding stack (settled)

- **PyO3 + maturin.**
- **New crate** `umol-py` (confirmed), producing Python import module `umol`
  (crate name vs import name differ, as is normal — `py-polars`, `pydantic_core`).
  `crate-type = ["cdylib", "rlib"]` (the `rlib` lets the wrapper layer be unit
  tested in Rust). `pyo3` with `extension-module`, `abi3` for one wheel across
  CPython minors.
- **Mixed maturin layout.** A `python/umol/` package (pure-Python `__init__.py`,
  `.pyi` stubs, docstrings) over a thin compiled `umol._native`. The Python layer is
  a home for stubs/docstrings and eventual conveniences — but see the parallelism
  principle below: the interface tracks Rust 1:1 for now, so the Python layer stays
  thin and divergence is deferred to alpha-user iteration, not front-loaded.
- **Wheels (decided):** keep `nauty` (no feature gate — automorphisms/symmetry is
  core and the local drag is nil: ~5s one-time, cached, never run by the slice).
  Distribute **prebuilt abi3 wheels from CI** so colleagues `pip install` a binary and
  never need a toolchain. abi3 ⇒ one wheel per platform covering Python 3.9+, so the
  matrix is tiny: `macos-arm64`, `macos-x86_64`, `manylinux_2_28-x86_64` (+ `aarch64`
  only if needed). Tool: `PyO3/maturin-action` (`maturin generate-ci github`). The
  single extra step vs a pure-Rust crate: `bindgen` needs **libclang inside the
  manylinux container** — one line, `before-script-linux: yum install -y clang`;
  macOS runners already have clang. The geometric wheel adds the `libmsym` submodule
  on top later. The `automorphisms`/`nauty` feature gate stays available for the one
  case not in scope — a source build on a platform no wheel is produced for. A
  libclang-free fallback is feasible without porting nauty/bliss: pure-Rust
  individualization-refinement crates exist (`canonical-form`, `graph_symmetry`/CNAP,
  `graphica`) to evaluate, and worst case a small IR reusing graph-core's existing
  `refine.rs` — molecular graphs (small, sparse, strongly colored) are the easy case,
  and nauty serves as a differential-testing oracle. Deferred; not needed now.

## The facade boundary (lead open question)

The entire read API is borrowed `*View<'a>` and two of the reaction/config types are
generic. None can be a `#[pyclass]` directly. The binding therefore needs an owned
facade — forced by PyO3, not chosen. The resolution below is faithful to the Rust
view design (backref + index, no molecule ownership).

**Key point: a handle is not ownership.** `Py<Molecule>` is a refcounted shared
handle to the Python `Molecule` object — the direct analog of `&'a MoleculeAst`,
not a clone. A view pyclass holding `Py<Molecule>` pays a refcount bump, never a
molecule copy. This dissolves the earlier worry that a Python view would have to own
(clone) the molecule.

Per surface, resolved:

- **Owned roots** (`MoleculeAst`, `ReactionSpanAst`) wrap directly: a `#[pyclass]`
  holds the value by move. Parsing returns an owned `Molecule`.
- **Reads → handle-view objects.** Each borrowed view maps to a pyclass holding
  `Py<Molecule>` + index; the accessor rebuilds the transient Rust view per call:
  `AtomViews<'a>` → `Atoms { owner }`, `AtomView<'a>` → `Atom { owner, index }`, so
  `mol.atoms[i].charge` reads structurally with no copy. This preserves the
  organizing structure the views exist for — on the side where it is free. Flat
  accessors on `Molecule` would recreate the "everything on `MoleculeAst`" cramming
  the views were built to avoid.
- **Mutations → methods on the `Molecule` facade** (or a builder), never through
  child view objects. PyO3 borrows are checked at runtime (RefCell-style):
  `borrow_mut` while any borrow is live raises `PyBorrowMutError`. Shared reads
  compose; mutation through a live child object is the aliasing hazard, so mutators
  live on the owner where only one borrow is in play.
- **Borrowed derivations** (`ReactionDerivation<'a>`): materialize an owned result.
  Since it borrows strictly for efficiency, add a Rust method that produces an owned
  products type (`into_owned()` / an owned reaction result) and wrap that — a real
  Python-pull change to the Rust API, and the owned type may serve Rust callers too.
- **Generic types** monomorphize at the boundary:
  - `EntitySpan<T>` → one Python type per entity family (atoms, bonds, and the
    overlay relations): `AtomSpan`, `BondSpan`, … Bounded, mechanical.
  - `GraphSymmetryConfig<C>` → expose only the concrete coloring(s) actually
    surfaced; for v1 the default coloring → a single config type. No speculative
    monomorphization over the `MoleculeColoring` space.

**Scope for v1** (parse + fingerprint + substructure + reactions, whole-molecule
ops): build handle-views only for `atoms`, `bonds`, and substructure match mappings
(query→target indices); defer the rest of the view family (aromatic, dative,
multicenter, stereo) until a workflow drives it.

## Feature gating lands on `umol-py`

The umbrella that doc 129 deferred effectively arrives here: `umol-py` is the first
single artifact with a feature surface. The `geometric` feature still gates the
`libmsym` submodule build; note the `graph` build is *not* toolchain-free either (it
pulls `nauty-Traces-sys`, `cc` + `bindgen`; see Verified current state). The `graph`
/ `geometric` features from doc 129 §i belong on this crate:

```toml
[features]
default   = ["graph"]
graph     = ["dep:umol-graph"]
geometric = ["dep:umol-geometric"]   # sole path that pulls umol-msym-sys → the C build
```

This does not reopen the Rust-umbrella question. The Python crate is the umbrella
and stays Python-only; whether Rust consumers also get a `umol` facade remains the
deferred doc-129 decision.

## Error mapping

Mirrors the three-tier model of doc 065 (layering settled; that doc's crate names
predate the reorg — `umol-data` is gone, the `UmolError` trait now lives in
`umol-utils`, but the tier structure stands):

- **Base Python exception `umol.UmolError`** ↔ the `UmolError` trait.
- **One subclass per tier-2 dispatch enum** (e.g. `GraphError` ↔ `GraphIrError`),
  with finer subclasses per tier-1 concern where useful (`ParseError`,
  `ValidationError`), built with `create_exception!` on a parent chain.
- **Carry structured fields onto the exception object** (`err.field`, `err.value`,
  `err.min`, `err.max`) — the FFI analog of doc 065's no-string-flattening rule. A
  `ValidationError::OutOfRange { field, value, min, max }` must not degrade to a bare
  message at the boundary.

Python-pull on the Rust errors: the binding wants boundary functions to return a
**stable, exhaustively-matchable tier-2 enum** wherever they map to a typed
exception. Reserve `Box<dyn UmolError>` for the boundaries where Python is content
with the generic base exception. This is the constraint to weigh in the in-progress
error refactor.

## Naming / surface policy

Python expects `snake_case`, keyword arguments, `__repr__`/`__eq__`. The mixed
layout carries most of this: the native `_native` extension keeps terse Rust names;
the `python/umol/` layer presents the PEP8 surface. Explicit policy (to record as
each item is decided, then harvest into the ubook — not to fold the whole ubook in
here):

- **Parallelism first:** the Python surface mirrors the Rust API (names, structure,
  AST-typed returns). Rust methods are already `snake_case`, so mirroring gives PEP8
  for free. Divergences (Python-only conveniences) are deferred to alpha-user
  iteration, applied on both sides together, not front-loaded in the Python layer.
- **Idiomatic Python wins over fidelity when they conflict (S4b).** Fidelity is the
  default, but a Rust pattern that reads unidiomatically in Python is adapted, not
  slavishly mirrored — idiomatic style outranks fidelity. E.g. Rust's `atoms[id]` →
  `&AtomAst` (data) vs `atoms.get(id)`/`iter()` → `AtomView` (view) asymmetry is an
  `Index`-signature artifact (Index must return a borrow; a view can't be borrowed);
  Python has no such constraint, so `[]`/`get`/iteration return one consistent atom
  type rather than reproducing the split. Where Rust exposes overloaded generic
  `Index<T>` (`mol[AtomId]`, `mol[BondId]`), one Python `__getitem__` dispatches on
  the key type at runtime, with `@overload` stubs for static typing.
- **Builders → constructor kwargs + `replace`, not `with_*` chains (S4b, idiomatic
  rule).** Rust's consuming `with_*` chain moves `self` (zero copies); Python can't
  move out of a live pyclass, so each `with_*` clones the whole value — an N-step
  chain does N copies. Python's idiom for immutable-copy-with-changes is a single
  `replace(*, field=…)` (`datetime.replace`, `dataclasses.replace`,
  `namedtuple._replace`): keyword args on the constructor to build in one shot, one
  `replace` (single clone) to derive a modified copy. Applied to `AtomAst` at S4b —
  the six `with_*` methods dropped for kwargs on `new`/`from_element` + `replace`.
- **Reads:** mirror the Rust accessor shape (`atom.charge()` → `atom.charge`);
  computed/fallible stay methods (`mol.fingerprint()`).
- **Integer types (S2a):** a mirror's boundary integer type **matches the underlying
  Rust type** (fidelity), *not* the workspace "default to u32." Python `int` is
  arbitrary precision, so width is invisible on the Python side — but the Rust type
  sets the accepted input range (out-of-range → `OverflowError`) and should track its
  source: `Element.atomic_number` is `u8` (as `ChemElement`), `IsotopeMassAst.Lit` is
  `u32` (as the AST), `ValueAst.Lit` is `i64`. The "default u32" is for *new* fields;
  a `usize` *count* still casts to `u32` at the `.len()` boundary (`atom_count`).
- **Dunders (corrected S2a):** PyO3 has no automatic `Display`/`Debug` → dunder
  bridge. `__str__` comes from `Display` via `#[pyclass(str)]` (or a `str = "…"`
  format string); `__repr__` has **no** `Display`/`Debug` auto-option (a PyO3 todo) —
  hand-write it as the eval-able constructor form (`ValueTerm.Lit(5)`), except simple
  unit-only enums which auto-get `__repr__ = "Class.Variant"`. `__eq__`/`__hash__` via
  Rust `PartialEq`/`Hash` (the `eq`/`hash` pyclass options) where the type has them.
- **Mutation:** on the owner facade or a builder (see the facade section), never
  through child view objects.
- **Rust alias for a wrapped foreign type:** crate-suffix + the exact Rust type name —
  `umol_ast::ast::MoleculeAst` → `AstMoleculeAst`, `umol_chem::element::Element` →
  `ChemElement`. Duplication (`AstMoleculeAst`) is accepted; the wrapper keeps the
  bare Rust name (`MoleculeAst`, `Element`).

## Documentation strategy

Four artifacts describe aspects of the same system. They are **structurally
separate and coordinated, not formally linked**:

| Artifact | Role | Register |
| --- | --- | --- |
| `discussion/` docs | historical archive of design conversations | unpolished |
| DSL spec (`umol-ast/spec/umol-dsl-spec.md`) | normative surface definition | RFC-2119; needs tightening (became loose) |
| doc comments (rustdoc) | in-code reference | module-structured |
| public doc (API/User Guide) | user background + conceptual overview | polished, for the whitepaper SI |

No formal cross-linking or generation between them; they are kept consistent by
hand. In particular the API Guide is **not** generated from doc comments — rustdoc
is a type-system-organized reference, the Guide is task/concept-organized; auto-
extraction would degrade both.

### Doc-comment discipline

Terse, at the public boundary, mechanically enforced:

- `#![warn(missing_docs)]` on `umol-py` and on each domain crate with user-facing
  types — turns "is the public surface documented" into a lint rather than a
  judgment.
- One summary sentence per public item (first line = rustdoc summary + search
  text). No prose on self-evident fields.
- Prose reserved for the non-obvious: invariants, `# Panics`, `# Errors`, and the
  genuine surprises (Undetermined-as-wildcard, canonicalize-on-compare,
  frame-relative stereo cosets).
- `//!` module docs carry the conceptual glue (the few key types and how they
  relate), not per-function repetition.
- Doctests on the load-bearing entry points only — compiled and run by
  `cargo test`, so examples cannot rot. Skip on leaf accessors.
- Calibration references for tone: `std`, `regex`, `hashbrown`, `ndarray`.

### Doc pipeline (settled)

- **Whitepaper**: LaTeX, hand-authored.
- **API/User Guide (SI)**: Markdown → Pandoc → PDF as a standalone document;
  combine PDFs with the whitepaper after the fact. PDF output is the target, not a
  LaTeX fragment. Chemistry via `mhchem`, code via `minted`/`listings`, diagrams via
  TikZ or pre-rendered includes, all through Pandoc's LaTeX writer.
- **Typst rejected**: immature chemical-formula and diagram support; not worth it
  against 20 years of existing LaTeX practice.
- **mdbook**: not used now; reserved for a possible future HTML User Guide over the
  same Markdown.
- **External repo, symlinked in, gitignored** (like `materials/`; symlinked under
  `whitepaper/` at top level). The symlink serves internal cross-referencing only and
  is not committed. The DSL spec stays canonical in this repo (it is normative for
  the parser and versioned with the code); the doc repo references it, not the
  reverse.

## Toolchain (verified present)

- **Binding**: `maturin` 1.14, `uv`, `pytest`, `ruff`, `pyright`. `mypy` and
  `qpdf`/`pdftk` are absent but substituted (`pyright` covers stub checking,
  `pdfunite` covers PDF merge) — no need to install them.
- **Dev interpreter**: build/test against Python 3.13 (`/opt/homebrew/bin/python3.13`
  via `uv venv --python 3.13`); the default `python3` is 3.9 (EOL). Keep `abi3-py39`
  as the wheel floor if broad support is wanted, but 3.9 is not the dev target.
- **Graph C build**: `nauty-Traces-sys` needs `cc` + `bindgen`/libclang — present on
  the dev machine (S0a built clean). Required for *any* umol-py build, not just
  geometric.
- **Docs**: `pandoc`, `xelatex`/`lualatex`/`latexmk`, `pdfunite`, `pygmentize`,
  `dot`; LaTeX `mhchem`/`chemfig`/`minted`/`listings`/`pdfpages`/`tikz` all present.
- **Deferred**: `cibuildwheel` (not installed) for prebuilt multi-platform wheels —
  needed for *graph* wheels too, since source builds require the C+libclang toolchain;
  plus the `libmsym` submodule for the geometric wheel.

## First slice (de-risk)

Graph-only. Staged from a narrow AST slice, not the whole crate:

- **i.** `MoleculeAst`, `AtomViews`/`AtomView`, `AtomAst`, and the value types inside
  `AtomAst` (`ElementAst`, `IsotopeMassAst`, `ValueAst`, `SpinStateAst`,
  `AtomConstraints`).
- **ii.** atom DSL parsing.
- then build outward.

**Value-type representation (decided): full structural mirror.** Each AST enum is
mirrored as a Python algebraic type, not reduced to a DSL string — `ValueAst`
(`Lit`/`LitSet`/`RangeFrom`/`RangeTo`/`Term`/`Predicate`/`Undetermined`), `ElementAst`
(`Lit`/`LitSet`/`NotSet`/`Var`), and the `ValueTerm`/`ValuePredicate` trees. Python can
construct and inspect patterns directly, without string round-trips — the fit for
programmatic pattern / reaction-network building. The mirror follows the actual
(canonical) variants, e.g. `ElementAst.NotSet({…})`. Cost: a large, recursive
Python surface that exposes the enum structure; this is the accepted trade for a
homoiconic Python pattern API. The DSL-string form remains available via the `*Dsl`
boundary types (step ii), as a parse/print path, not the primary representation.

**`Element` (decided): a pyclass** wrapping `umol_chem::Element`, with `.symbol` and
`.atomic_number` (room for per-element data), returned from `ElementAst.Lit`.

**Interface parallels Rust 1:1 (decided).** The Python API mirrors the Rust API —
same names, same structure, same AST-typed returns (`atom.charge` yields a
`ValueAst`, not an `int`). The tension between AST-typed and primitive returns is
real and identical on both sides of the FFI; it is resolved by co-iterating the two
interfaces with alpha users, not by giving Python a divergent primitive surface.
Consequence: keep the binding thin and cheap to change.

**Read surface (decided): one Python class per variant, children as attributes,
`__match_args__` set.** This is the mainstream Python representation of a recursive
AST — the stdlib `ast` module is exactly this (one class per node under a base;
`_fields` as attributes; `__match_args__` since 3.10, so `case BinOp(l, op, r)`
works — AST processing is *the* canonical Python `match` use case), and polars'
Rust-backed `DataType` enum (`List(inner)`, `Struct(fields)`, …) is the same shape.
So `atom.charge` yields a `ValueAst` whose variant is a class (`Lit(0)`,
`RangeFrom(1)`, …), inspectable by attribute, `isinstance`, or `match`. The rejected
alternatives: per-field `as_<variant>()` accessors (no Python prior art — a Rust
idiom) and a tag-field/`type`-string discriminator (the pydantic-core / tree-sitter
serialization idiom, for external data, not a built-and-matched AST). Lesson stolen
in reverse from polars: **every variant is a consistent instance** — never some
variants classes and some instances (polars' `isinstance` broke on that).

**Mechanism (resolved by the S1a spike, 2026-07-05): native PyO3 complex enums
throughout.** The open question was whether a native complex enum can hold a
recursive `Py<Self>` / `Vec<Py<Self>>` field. The `ValueTerm` spike answered **yes**:
the enum compiles, and from Python it constructs (`ValueTerm.Lit(5)`,
`ValueTerm.Neg(...)`, `ValueTerm.Sum([...])`), reads fields (`_0`/`_1`), and `match`es
natively. So the pure-Python-over-opaque fallback is **not needed** — every value
mirror type is a native complex enum. AST recursion (`Box<Self>`/`Vec<Self>`) maps to
`Py<Self>`/`Vec<Py<Self>>`; the two are *distinct* types bridged by `pub(crate)`
`from_ast(py, &AstT) -> PyResult<T>` / `to_ast(&self, py) -> AstT` conversions (the
parallel-representation cost native enums carry — accepted, since native gives the
match ergonomics and the conversions are mechanical and roundtrip-tested).

Notes from the spike + S1b–d (native-enum mechanics that recur across every mirror):

- Tuple variants expose fields as `_0`/`_1` (positional `match` works:
  `case ValueTerm.Lit(n)`); `from_ast` allocates one Python object per node via
  `Py::new`; GIL-holding Rust tests use `Python::attach` (PyO3 0.29 renamed `with_gil`)
  under the `auto-initialize` dev-dependency (test-only, absent from the wheel build).
- **Unit variants are rejected** in a complex enum — use an empty-tuple variant
  (`Undetermined()`); Python constructs/matches `ValueAst.Undetermined()`. Recurs for
  `ElementAst` and `IsotopeMassAst` (`Natural`).
- **A simple enum used as a *field* of a complex enum** (operands like `RelOp`/`MemOp`)
  needs `#[derive(Clone)]` + `#[pyclass(…, from_py_object)]` — the getter clones it, the
  constructor extracts it. Opposite of a held leaf like `Element`, which drops `Clone`.
  Keys on: *is this type extracted from Python as a value?*
- `from_ast`/`to_ast` stay `#[allow(dead_code)]` until a live `#[pymethods]` caller
  reaches them (the `AtomAst` field getters at S3); the roundtrip tests keep them
  covered in the test build meanwhile. (Retired at S3.)
- **`Py::new` builds a *base*-type complex-enum instance** (S3): its Python-visible
  variant fields (`_0`, …) and `match` are absent. Wrap every nested `Py<…>` child in
  a `from_ast` with `IntoPyObject` instead — the `into_py_variant(py, value)` helper
  (`convert.rs`). Insidious because the Rust-side roundtrip tests (`to_ast`) never
  exercise Python-side variant access, so it stayed hidden until a getter (`atom.spin`)
  read a nested child; the fix carries a Python-side regression test.

Sequencing within the slice: `Element` → the scalar-field enums (`ElementAst`,
`IsotopeMassAst`, `ValueAst` + `ValueTerm`/`ValuePredicate`, `SpinStateAst`) →
`AtomConstraints` (13 variants + sub-ASTs) last.

### Wrapping strategies

Whether a wrapper carries AST-conversion methods is set by its wrapping *strategy*,
chosen per type — **not** by whether the Rust type is a `Lattice`/homoiconic AST, nor
by `ast`-module membership (`MoleculeAst` is both, and is *held*, not mirrored).

- **Hold-the-value** — `#[pyclass] struct W(RustValue)` stores the Rust value; reads
  are accessors over it. No parallel representation, so **no AST-conversion methods**;
  a getter exposing a value-algebra field runs *that field's* conversion on access.
  Used for owned roots (`MoleculeAst`, `AtomAst`) and leaf vocabulary (`Element` over
  `ChemElement`).
- **Structural mirror** — `#[pyclass] enum W { …variants… }`, a distinct native enum
  reproducing the Rust enum's variants (recursion re-expressed `Box<Self>` → `Py<Self>`).
  Being a separate type, it **carries conversions to/from the AST** (the `from_ast` /
  `to_ast` pair on `ValueTerm`). Used for the value/pattern algebra (`ValueAst`,
  `ValueTerm`, `ValuePredicate`, `ElementAst`, `IsotopeMassAst`, …).

Criterion: **does Python construct / `match` it by variant?** Yes → structural mirror
(it must be a umol-py-defined `#[pyclass] enum` — you can't `#[pyclass]` a foreign enum,
and the recursion becomes `Py<Self>`). No → hold the value. The mirrored set coincides
with the homoiconic `Lattice` algebra, but the *cause* is the destructure question, not
`Lattice`-ness.

That two-way cut is the coarse criterion; the concrete kinds observed through S4b:

| Kind | Examples | `#[pyclass]` form | field-hold / arg | AST bridge |
| --- | --- | --- | --- | --- |
| **Leaf newtype** | `Element` | `(eq, hash, frozen, from_py_object)` struct over a *vocab* type | by value | `From<Chem>` / `From<&Self>` — not an AST |
| **Id newtype** | `AtomId` | `(eq, hash, frozen, from_py_object)` struct over an id | by value (index arg) | `From<AstId>`; `#[new](u32)`, `.index` |
| **Simple enum** | `RelOp`, `MemOp` | `(eq, from_py_object)` enum, unit variants | by value | `from_ast`/`to_ast`, **py-free** |
| **Complex enum** | `ValueTerm`, `ValueAst`, `ValuePredicate`, `ElementAst`, `IsotopeMassAst` | `#[pyclass]` enum, data variants; unit → `()` | `Py<Self>` for recursion, by-value data else | `from_ast`/`to_ast` (**py** iff it holds `Py` children) |
| **Mirror struct** | `SpinStateAst` | `#[pyclass]` struct, `#[pyo3(get)]` + `#[new]` | `Py<mirror>` fields | `from_ast`/`to_ast` (**py**) |
| **Hold-the-value struct** | `MoleculeAst`, `AtomAst` | `(eq)` struct `W(AstT)` | holds `AstT`; getters mirror fields on read | **`inner`/`from_inner`** (borrow out / move in — no rebuild) |
| **Handle view** | `AtomView`, `AtomViews` | `#[pyclass]` struct `{ owner: Py<W>, id }` | — | none — rebuilds a transient Rust view via `owner.inner()` |
| **Iterator** | `AtomViewIter` | `#[pyclass]` struct, `__iter__`/`__next__` | holds `owner` + id iter | none — internal, unexported, not `add_class`'d |

The wrapping-infra methods that bridge each kind to the AST:

| Method | On | Role |
| --- | --- | --- |
| `from_ast(py?, &AstT) -> PyResult<T>` / `to_ast(&self, py?) -> AstT` | simple + complex enums, mirror struct | structural **convert** mirror ↔ AST (rebuilds; `py` when it allocates `Py` children) |
| `inner(&self) -> &AstT` / `from_inner(ast: AstT) -> Self` | hold-the-value structs | **wrap/unwrap** the held value (borrow out, move in — no rebuild, no `py`) |
| `into_py_variant(py, value) -> PyResult<Py<T>>` | every `from_ast` with `Py<…>` children | wrap a complex-enum child as its **variant** instance — `Py::new` yields a base instance (S3 bug) |
| `From<Chem>` / `From<&Self> for Chem` | leaf newtypes | vocab-type conversion (not an AST bridge) |

The two bridge *pairs* track the two strategies: `from_ast`/`to_ast` = structural
**conversion** (the mirror rebuilds a parallel structure), `inner`/`from_inner` =
**wrap/unwrap** of a held value (trivial). Which pair a type carries tells you its
strategy at a glance; a handle view carries neither — it reads through `owner.inner()`.

## Decisions

Settled:

- **Stack** — PyO3 + maturin, mixed layout, abi3, `cdylib`+`rlib`.
- **Feature gating** — `graph` / `geometric` on `umol-py`; default graph-only wheel;
  geometric wheel + `cibuildwheel` later.
- **Read pattern** — handle-view objects holding `Py<Molecule>` + index (faithful to
  the backref views, no molecule copy); mutations as methods on the owner facade /
  builder, never through child views (runtime-borrow aliasing).
- **Interface parallelism** — Python API mirrors the Rust API 1:1 (names, structure,
  AST-typed returns, not primitives); divergence deferred to alpha-user co-iteration.
- **Naming** — maximum Rust fidelity: verbatim type names incl. the `Ast` suffix
  (`MoleculeAst`/`AtomAst`/`AtomView`/`ValueAst`, no bare `Molecule`/`Atom`); `Id`
  types as Python newtypes, not bare `int`. Revisit on alpha-user feedback.
- **Value representation** — full structural mirror. Read surface is **one class per
  variant** (children as attributes, `__match_args__` set; `isinstance`/`match` both
  work) — the stdlib-`ast` / polars-`DataType` shape; not `as_*` accessors, not a
  tag-field discriminator. Every variant a consistent instance. `Element` is a
  pyclass (`.symbol`/`.atomic_number`).
- **Mirror mechanism** — **native PyO3 complex enums throughout** (S1a spike, 2026-07-05:
  a complex enum *does* hold recursive `Py<Self>`/`Vec<Py<Self>>` fields — construct,
  field-read, and `match` all work). Pure-Python fallback dropped. Each mirror type is
  a distinct native enum bridged to the AST by `pub(crate)` `from_ast`/`to_ast`.
- **Borrowed derivations** — `ReactionDerivation<'a>` gets an owned Rust result type
  the binding wraps (Python-pull change to the Rust API).
- **Generics** — monomorphize at the boundary: per-entity-family `EntitySpan<T>`
  types; only the surfaced coloring(s) for `GraphSymmetryConfig<C>`.
- **Errors** — Python exception hierarchy mirrors the three tiers; structured fields
  carried onto exception objects; typed exceptions want tier-2 enums, `Box<dyn>`
  only for the generic base.
- **Surface policy** — native layer keeps Rust names, Python layer is PEP8; getters
  for scalar fields, methods for computed/fallible; `__eq__`/`__hash__` from Rust
  `PartialEq`/`Hash`, `__str__` from `Display` via `#[pyclass(str)]`, `__repr__`
  hand-written (no `Display`/`Debug` auto-bridge). Recorded here, harvested into the
  ubook later.
- **`AtomAst` is an immutable value (S4b).** Construct via kwargs
  (`AtomAst(element, *, charge=…)`, `element` a concrete `Element` or an `ElementAst`),
  derive via one-copy `replace(*, field=…)` — **no field setters**, despite the Rust
  `pub` fields (which are storage, not a field-poke API — umol itself prefers
  accessors). Reasons: the surface is uniformly value-semantic, and a lone mutable
  type collides with `from_atoms`' clone-on-insert (`atom.charge = …` would silently
  not touch an already-inserted molecule); RDKit-style in-place editing is
  molecule-*owned* (`mol.atoms[i]`), so it belongs on the mutable-molecule facade (the
  deferred mutation story), not on a detached atom.
- **Literal coercion via `*Arg` unions (S5+).** Wherever the Rust builders take
  `impl Into<T>` (accepting a bare literal), the Python argument accepts the literal
  *or* the mirror, via a `#[derive(FromPyObject)]` union tried in order: `ElementArg`
  (`Element` | `ElementAst`), `ValueArg` (`int` → `ValueAst::Lit` | `ValueAst`),
  `IsotopeMassArg` (`int` → mass | `IsotopeMassAst`). So `AtomAst(Element("C"),
  charge=-1, isotope_mass=13)` works alongside the explicit mirror form. **Naming:**
  `*Arg` for these binding coercion inputs — `*Input` is reserved for the DSL side.
  Variants are `Ast` (the wrapper) / `Lit` (the literal). Same runtime-dispatch
  mechanism as the `__getitem__` overloads; input-only (getters still return mirrors).
  Applied to the atom value fields, `SpinStateAst(unpaired, multiplicity)`, and
  `RingMembershipAst.count`. The shared `ValueArg` lives in `value.rs` (with two coercions:
  `to_ast` → `AstValueAst` for `with_*` builders, `to_py` → `Py<ValueAst>` for mirror
  structs that store the field); `ElementArg`/`IsotopeMassArg` stay in `atom.rs`.
- **Doc separation** — four artifacts, coordinated by hand, no generation between
  them; Guide is not derived from doc comments.
- **Doc-comment discipline** — `missing_docs` lint, one-line summaries, prose only
  for the non-obvious, doctests on entry points, module `//!` for glue.
- **Doc pipeline** — Markdown + Pandoc → PDF SI; combine with the LaTeX whitepaper;
  external gitignored symlinked doc repo; spec stays canonical in-repo.
- **Toolchain** — verified present; dev/test on Python 3.13.

Open:

- **First AST slice** — `Molecule` (owned `MoleculeAst`), the atom read surface
  (`AtomViews`/`AtomView` handle + owned `AtomAst`), and the value-type mirror classes
  (`Element`, `ElementAst`, `IsotopeMassAst`, `ValueAst`/`ValueTerm`/`ValuePredicate`,
  `SpinStateAst`), then `AtomConstraints`; step ii adds atom DSL parsing. First plan
  task is the mechanism spike (recursive native complex enum on `ValueTerm`).
- **Owned reaction-result type** — its Rust name and shape (the `into_owned()`
  target for `ReactionDerivation`), when reactions are reached.

## Implementation plan (slice i + ii)

Scope: slice i (`MoleculeAst`, atom read surface, `AtomAst`, its value types) + ii
(atom DSL parsing). All work is **additive in the new `umol-py` crate** — existing
workspace crates are consumed read-only, so the workspace stays green throughout.
"Green" per stage = `cargo test -p umol-py` (wrapper/conversion unit tests) **and**
the pytest suite (via `maturin develop`) both pass. **Naming: maximum Rust fidelity**
(decided) — Python type names mirror Rust verbatim (`ValueAst`, `AtomAst`, `AtomView`,
`MoleculeAst`; no bare `Molecule`/`Atom`), and `Id` types are Python newtypes, not
bare `int`. Revisit on alpha-user feedback.

### S0 — Scaffold + leaves (green) **Done**

- **S0a** — crate `umol-py`: workspace member; `Cargo.toml` (`pyo3` extension-module
  + `abi3-py39`, `cdylib`+`rlib`, `umol-ast`/`umol-chem` behind a `graph` feature);
  `pyproject.toml` (maturin, mixed layout `python/umol/`); `#[pymodule] _native`;
  `python/umol/__init__.py`; `uv venv --python 3.13`; a pytest asserting `import
  umol`. *Additive.* [dep: —]
- **S0b** — `Element` pyclass (module `element`): `Element(symbol)` / from atomic
  number, `.symbol`, `.atomic_number`, `__repr__`/`__eq__`/`__hash__`. Rust
  conversion test + pytest roundtrip. *Additive.* [dep: S0a]
- **S0c** — skeletal `MoleculeAst` pyclass (module `molecule`): wrap the Rust value,
  `empty()`, `.atom_count`, `__eq__`/`__repr__` — proves the owned-root wrapper.
  *Additive.* [dep: S0a]

### S1 — Mechanism spike + value-expression core (green) **Done**

- **S1a** — **spike & `ValueTerm`** (module `value`): prototype `ValueTerm` as a PyO3
  native complex enum with recursive `Py<Self>` / `Vec<Py<Self>>` fields; confirm
  construct + `match` from Python. Record the outcome in this doc. If native fails →
  pure-Python variant classes over an opaque `#[pyclass] ValueTerm` core. Delivers
  `ValueTerm` + `__match_args__` + roundtrip tests either way. *Additive.* [dep: S0a]
- **S1b** — `RelOp`, `MemOp` simple pyclass enums (module `value`). *Additive.*
  [dep: S0a]
- **S1c** — `ValuePredicate` mirror (recursive; `Rel`/`Mem`/`Not`/`And`/`Or`).
  *Additive.* [dep: S1a, S1b]
- **S1d** — `ValueAst` mirror (7 variants incl. `Term`/`Predicate`), with the
  Rust-parallel constructors (`lit`/`lit_set`/`range_from`/`var`/`undetermined`).
  *Additive.* [dep: S1a, S1c]

### S2 — Atom-field value types (green) **Done**

- **S2a** — `ElementAst` (5 variants; `Lit` uses `Element`, `Var` uses `MemOp`).
  *Additive.* [dep: S0b, S1b]
- **S2b** — `IsotopeMassAst` (5 variants incl. `Natural`). *Additive.* [dep: S0a]
- **S2c** — `SpinStateAst` (`{unpaired, multiplicity}`; fields are `ValueAst`).
  *Additive.* [dep: S1d]

### S3 — `AtomAst` owned (green) **Done**

- **S3a** — `AtomAst` pyclass (module `atom`): mirror-typed fields (element /
  isotope_mass / charge / implicit_hydrogens / lone_pairs / spin); Rust-parallel
  constructors (`new`, `from_element`, `with_*`); getters; `__repr__`/`__eq__`.
  `constraints` deferred to S5c. *Additive.* [dep: S1d, S2a, S2b, S2c]

### S4 — Molecule atom read surface (green) **Done**

- **S4a** — `AtomId` newtype pyclass + `AtomViews` / `AtomView` handle pyclasses
  (module `atom`): views hold `Py<MoleculeAst>` + `AtomId`; `AtomView.id` is an
  `AtomId`, molecule indexing takes one; getters return the mirror types by reading a
  transient Rust `AtomView`. Bond / topology-derived methods (`valence`, `neighbors`,
  …) are out of scope (need the bond slice). *Additive.* [dep: S0c, S3a]
- **S4b** — enrich the `MoleculeAst` pyclass: `.atoms` → `AtomViews`, indexing /
  iteration → `AtomView`, and `from_atoms_and_bonds`-style construction taking Python
  `AtomAst`s. *Additive.* [dep: S4a, S3a]

### S5 — `AtomConstraints` (green) **Done**

- **S5a** *(done)* — constraint sub-ASTs (module `constraint`): `RingScope`,
  `RingMembershipAst`, `AromaticValenceAst`, `MulticenterValenceAst`. And **the full
  stereo sub-tree** (module `stereo`, decided over the concrete-only/defer options):
  `TetrahedralStereoAst` → `StereoCosetAst` → `StereoTerm` (recursive) → `Permutation`
  (hold-the-value over `umol_perm::Permutation`, a cross-crate dep added on `umol-py`).
  Finding: PyO3 maps `Vec<u8>` → `bytes`, so the permutation image surfaces as
  `Vec<u32>` (list of ints), consistent with the `u32` coset indices. *Additive.*
- **S5b** *(done)* — `AtomConstraint` (13 variants, `.kind`), `AtomConstraintKind`
  (simple enum, keyed lookup), `AtomConstraints` container (hold-the-value:
  `len`/iteration/`get(kind)`/`contains(kind)`, built from a list). Consuming the
  sub-ASTs here made every `from_ast`/`to_ast` live (cleared the S5a `dead_code`).
  *Additive.*
- **S5c** *(done)* — `constraints` getter on `AtomAst` and `AtomView`, plus a
  `constraints=` kwarg on `AtomAst(...)`/`replace` (wipe-and-set via the `pub`
  `constraints` field — replace-semantics, not `with_constraints`' add). `atom.rs` +
  `stereo.rs` now complete → module-level `#[allow(clippy::absolute_paths)]` for the
  `#[pyclass(hash)]` macro false-positives. **Whole atom slice green; clippy clean.**

### S6 (step ii) — atom DSL parsing (green) **Done**

- **S6a** *(done)* — `AtomAst.parse(str)` (staticmethod), `str(atom)`, and
  `repr(atom)` = `AtomAst.parse('<dsl>')`, wrapping `AtomAst`'s own `FromStr`/`Display`
  (the convenience shortcuts over the `FromAst`/`IntoAst` raise/lower system; the full
  raise/lower binding is deferred — Python todo 5). Roundtrip **canonicalizes** (`#c+1`→
  `#c+`, constraints reorder to key order), so it is `str`∘`parse` stable, not
  byte-verbatim. Parse errors → the dedicated, catchable `umol.ParseError`
  (`create_exception!`, `PyException` base); the three-tier hierarchy is a later slice.
  *Additive.* [dep: S3a, S5c]

### Critical path & deferrals

- **Critical path:** S0a → S1a → S1c → S1d → S2c → S3a → S4a → S4b → S5a → S5b → S5c
  → S6a. The spike (S1a) is the highest-risk item and sits first; everything
  downstream assumes its mechanism outcome.
- **Out of scope for this slice** (subsequent slices, building outward): `AtomView`
  bond & topology methods; molecule-level DSL/EDN parsing; fingerprints; substructure
  search; reactions; the full error/exception hierarchy. None are needed for the
  i + ii deliverable.
- **Naming (decided):** maximum Rust fidelity — verbatim `MoleculeAst`/`AtomAst`/
  `AtomView`/`ValueAst` (no bare `Molecule`/`Atom`), `Id` types as newtypes. Revisit
  on alpha-user feedback. *(Revised for `Id` types — see the mutation-facade note below:
  the Python-facing `AtomId` was retired in favor of bare ints.)*

### Post-slice mutation + ergonomics + API polish (2026-07-08)

Built out after i + ii; all umol-py, workspace green throughout. All exploratory — the
immutable-`AtomAst` premise was **not** treated as settled (Python leans mutable).

- **Mutable atom facade.** `AtomView` (what `mol.atoms[i]` returns) is a live mutable
  handle: settable field properties (`mol.atoms[i].charge = …`) route through `atom_mut`
  in place, and `.constraints` returns a live `AtomConstraintsView` whose mutators
  (`set`/`remove`/`update`, settable accessors) write straight through — reads borrow the
  item directly, no whole-container clone. The view backs onto **either** a molecule-atom
  or a standalone `AtomAst` (`ConstraintsBacking` enum; backref lives in the view, not the
  atom — deliberately unlike RDKit). `AtomAst` fields became settable and `replace` was
  retired; `mol.atoms[i] = atom` (`__setitem__`) is whole-atom replace. Model: the
  molecule is a mutable container of atom-*values* — you don't mutate an atom in place,
  you replace the value at a slot (or edit through the owner).
- **`AtomId` retired on the Python side** — bare ints (`mol.atoms[0]`, `atom.id → int`);
  the Rust `AstAtomId` newtype stays.
- **Literal ergonomics.** Settable per-kind accessors so `mol.atoms[i].constraints.
  aromatic_valence = 1` *is* `#a1`: `False → No*`, `int → positive` (aromatic/multicenter
  valence); `TetrahedralStereo.Cw|Ccw` enum (coset `Lit(1)`/`Lit(0)`); `ring_count` (All)
  property + `ring_size_count[size]` subscript proxy (`RingSizeCounts`); `as_lit()` on
  `ValueAst`/`ElementAst`/`IsotopeMassAst`; deeper `*Arg` int coercion. `E.H`/`E.As`
  element shorthand (`python/umol/` `__getattr__`) — Python todo 2, dynamic form.
- **API-review polish.** Value `__eq__`/`__hash__` (via `to_ast`) + eval-able `__repr__`
  across all mirrors — so `ValueAst.Lit(1) == ValueAst.Lit(1)` (was identity-false); dict
  protocol on the constraint containers (`in`/`[]`/`del`, alongside `get`/`remove`);
  dropped `is_empty`, `atom_count` (→ `len(mol.atoms)`), `AtomViews.get`; `AtomView.asdict`
  parity; `RingSizeCounts` len/iter/contains. Tests migrated `match`-as-equality → `==`.
  Still open: `append`/`extend` (todo 6), `__repr__` on a molecule DSL, `MoleculeAst`
  `Hash` (Rust fold-back below).

### Second API-review pass (2026-07-09)

A deeper sweep over the now-richer container surface; findings 1–4 implemented (umol-py
green: 147 Rust / 162 pytest / clippy clean), 5 and 6 deferred to the interning design.

- **1 — `AtomConstraintsAst` unhashable (bug).** The value container had a hand-written
  `__eq__` + `__hash__` but is mutable (`set`/`pop`/`update`) — an identity-inconsistent
  hash footgun. Fixed by converting to `#[pyclass(eq)]` + `derive(PartialEq)` (matches
  `AtomAst`/`MoleculeAst`), which makes it value-equal but **unhashable**. Note: merely
  deleting `__hash__` would *not* have worked — a hand-written `__eq__` leaves the default
  identity `__hash__` in place; only the `eq` macro nulls it (verified against PyO3 0.29).
- **2 — constraint containers are proper mappings.** `__iter__` now yields *keys*
  (`AtomConstraintKey`), matching `RingSizeCounts` and dict convention (it yielded values
  before, while `[]`/`in`/`del`/`get` were all keyed — a half-mapping). Added
  `keys()`/`values()`/`items()` on both `AtomConstraintsAst` and `AtomConstraintsView`.
- **3 — `atom.constraints` setter accepts a view.** Was value-only, so
  `dst.constraints = src.constraints` (RHS a live view) failed. New `ConstraintsArg`
  coercion snapshots either a value container or a view; both `AtomAst`/`AtomView` setters
  take it.
- **4 — `remove` → `pop`.** Remove-by-key-returning is `dict.pop`, not `list.remove`;
  renamed on both containers. `get(key, default=None)` gained the dict second arg (returns
  the object, else the default) — return type widened to `PyObject`.
- **5 — handle `__eq__` (pending, interning-entangled).** Giving `AtomView` /
  `AtomConstraintsView` an `__eq__` forces a value-vs-reference equality choice — the same
  identity axis as finding 6. Deferred to the interning discussion rather than guessed.
- **6 — `AtomViews.__contains__` (deferred).** `atom in mol.atoms` is a handle-membership
  question tied to interning; left for that design.

### Third API-review pass (2026-07-09)

A final sweep — the surface is largely clean by now. 1–3 implemented (umol-py green:
147 Rust / 165 pytest / clippy clean); 4 deferred.

- **1 — `AtomViews` negative indexing.** `mol.atoms[-1]` raised `OverflowError` (`usize`
  index). Both `__getitem__`/`__setitem__` now take `isize` and normalize `index < 0 →
  len + index` via a shared `resolve_atom_index` (bounds-check → `IndexError`). Slicing
  (`mol.atoms[1:3]`) is still unsupported — separate, heavier, left open.
- **2 — `RelOp` / `MemOp` / `TetrahedralStereo` made hashable.** They were
  `#[pyclass(eq, from_py_object)]` (macro `eq` nulls hash) while every other immutable
  value type hashes; now `#[pyclass(eq, hash, frozen, from_py_object)]` + `derive(Eq,
  Hash)`, matching `Element`. Their auto `__repr__` was already fine.
- **3 — `MoleculeAst.__repr__`** said `atom_count=N`, resurrecting the retired `atom_count`
  term; now `atoms=N` (consistent with `len(mol.atoms)`).
- **4 — `Permutation.image` is a method while `degree` is a property (deferred).** Both are
  no-arg pure reads; a list-returning property is borderline (cf. numpy `.tolist()`).
  Deferred to be taken up with the top-level stereo work (stereo atoms / stereo bonds),
  where the `Permutation` surface is revisited as a whole.

## Findings to fold back into Rust

Rust-side changes surfaced while wrapping (the Python-pull of the co-iteration
principle). Collected here as a running list; applied as a batch on the Rust side,
not mid-slice.

1. **`MoleculeAst` has no `Hash` impl** (S0c). `PartialEq`/`Eq` are hand-written
   (excluding `rings_cache`), but there is no `Hash` — yet the struct doc comment
   (`umol-ast/src/ast/molecule.rs:54`) claims the cache is "excluded from PartialEq /
   **Hash**." Consequence: the Python `MoleculeAst` is unhashable (`#[pyclass(eq)]`
   only). Fold-back: add a `Hash` impl mirroring the `PartialEq` field set (exclude
   `rings_cache`), then the wrapper takes `#[pyclass(eq, hash, frozen)]`; fix the
   stale comment.
2. **Atom-accessor surface: asymmetry + redundancy** (S4b — *to think about, not a
   settled fold-back*). Two smells in the by-id atom accessors:
   - *Asymmetry of return type.* `Index<AtomId>` → `&AtomAst` (intrinsic data) but
     `get(id)`/`iter()` → `AtomView` (contextual view). Partly forced — `Index` must
     return a borrow and a view can't be borrowed (it's built on demand) — but
     surprising even in Rust (indexing and `get` returning different types).
   - *Redundancy (TIMTOWTDI).* Four distinct by-id paths — `mol.atom(id)` → view,
     `mol[id]` → `&AtomAst`, `mol.atoms()[id]` → `&AtomAst`, `mol.atoms().get(id)` →
     view — and `mol[id]` fully duplicates `mol.atoms()[id]`. Does the `MoleculeAst`
     `Index` earn its keep over the `AtomViews` one?
   Open Rust-side question: consolidate toward one consistent return + fewer paths
   (e.g. a `view(id)` method, drop/rename the `Index` impls, one owner for indexing).
   Flagged for consideration; the binding meanwhile exposes a single lean surface —
   `mol.atoms[id]` / `.get(id)` / iteration, all returning the view — per the
   idiomatic-Python rule.
3. **`umol_perm::Permutation` image `u8` → `u32`** (S5a). `Permutation` stores
   `image: [u8; MAX_DEGREE]` (`MAX_DEGREE = 6`); the `u8` is a size optimization that's
   moot at that degree, and it deviates from the codebase's `u32` index convention
   (`AtomId`, coset `Lit`). Consequence: the binding round-trips `u8`↔`u32` at the
   boundary (image surfaced as `Vec<u32>`) with a silent `as u8` truncation on invalid
   input. Verified there's **no `u8`-specific logic** — only `as u8`/`as usize` casts —
   so the change is mechanical: `[u32; N]`, `from_image(&[u32])`, internal `Vec<u32>`;
   ripple is `from_image`'s signature + ~2 call sites in `umol-ast/symmetry.rs`, and the
   binding simplifies (drops both casts). Recommend doing it; batch with the other
   fold-backs (not mid-slice).
4. **Rename the exact-participant-set lookup `connecting` → `of`** (surfaced planning the
   relation bindings, doc 140; resolved 2026-07-09). `connecting` read well for a bond but
   not for the whole-set lookups (an aromatic system does not *connect* its members). The
   lookup is really "the entity keyed by this exact participant set" — a terse,
   converter-like selector, and the common case (vs the rarer `incident`/`induced`), so it
   earns the shortest name. Rust: `connecting` → `of`, `connecting_id` → `of_id`
   (`mol.bonds().of(first, second)`). Python mirrors it. Call convention: positional for a
   fixed symmetric pair (`bonds.of(0, 1)`), one iterable for a variable set
   (`aromatic_systems.of({…})`), keyword roles for a birelation (`dative_bonds.of(acceptor=…,
   donors=…)`, `stereo_atoms.of(site=…, ligands=[…])`). Note the lookup is multiset-matched,
   so `ligands` (which can repeat virtual ligands) must be an ordered sequence, not a set.
5. **`(a, b)` → `(first, second)` on participant-pair methods.** Ongoing normalization —
   the ast view `connecting` methods (bond, noncovalent) already use `first`/`second`;
   `umol-graph-core` `Graph::find_edge` / `add_edge` done (2026-07-09). Remaining:
   `find_bond_by_participants(a, b)` in the DSL layer (`dsl/namespace.rs` trait + impls,
   `dsl/reaction.rs`) + its call sites. Decided to leave as-is: the symmetric orbit
   predicates (`same_orbit` / `same_proper_orbit` / `same_star_orbit`) keep `a, b` (no
   first/second role in a symmetric relation); single-atom `a` params may become `atom`
   (non-urgent).
6. **Uniform `*ViewMut` across all eight families + retire the bulk `&mut` iterators**
   (surfaced planning the relation bindings, doc 140; interning-relevant, doc 114). Today
   the mutable accessors are asymmetric: `atom_mut`/`bond_mut` return a `*ViewMut` struct
   (`BondViewMut` usefully carries endpoints; `AtomViewMut` is `{ id, &mut ast }`), but the
   six relation `*_mut(id)` return a **bare `&mut XAst`**. Give the six relations a `*ViewMut`
   too (`id` + participants + `&mut ast`, mirroring their read `*View`). This is not only
   symmetry: 114's interning plan re-interns on the **view guard's `Drop`**, so a bare
   `&mut XAst` (no guard) is a latent interning blocker — uniform `*ViewMut` is the
   prerequisite. Same pass: retire the eight raw `&mut`-iterator `*s_mut()` (closure- or
   replace-based instead), per 114's discipline. The binding is unaffected (it edits only
   through per-id `atom_mut(id)`), so this is Rust-internal.

## Python-side todos

Binding-side work deferred and batched (distinct from the Rust fold-backs above —
these change only umol-py, not umol-ast).

1. **`__repr__` / `__str__` on the value mirrors + `AtomAst`** (S2a; pulled forward by
   the S4b `AtomAst.__repr__` request). `ValueTerm`, `ValueAst`, `ValuePredicate`,
   `ElementAst`, `IsotopeMassAst` (complex enums) and `SpinStateAst`, `AtomAst`
   (structs) currently have PyO3's default `<… object at 0x…>` repr (simple enums like
   `RelOp` auto-get `Class.Variant`; `Element`/`AtomId` hand-write it). **These are
   coupled**: a constructor-style `AtomAst.__repr__` (`AtomAst(ElementAst.Lit(…),
   charge=ValueAst.Lit(-1))`) embeds `repr()` of each field mirror, so it's only useful
   once the mirrors have reprs — do the set together. Each `__repr__` is the eval-able
   constructor form (recursing via Python `repr()` on children); separately decide
   whether `__str__` should be the compact **DSL string** (the AST value types impl
   `Display` as the value-expr, reachable via `to_ast`) — ties into the S6 DSL-string
   surface.
5. **`FromAst`/`IntoAst` (raise/lower) binding** — the primary, complete AST↔DSL
   conversion (perception boundary + structural raise/lower, parameterized by the
   `*Defaults` config); `Display`/`FromStr` — which S6a wraps for `str`/`repr`/`parse` —
   are convenience shortcuts layered on top, *not* a separate "raw" path. Exposing the
   full raise/lower surface to Python (ground a parsed pattern, lower with
   default-eliding, thread the defaults config) is a distinct, larger binding piece.
   Deferred; S6a delivers only the string shortcuts.
3. **Pickling** (`__reduce__` / `__getstate__`+`__setstate__`) — not automatic for
   pyclasses. Not needed now; may matter later for multiprocessing/distribution of the
   large reaction networks (atoms/molecules crossing process boundaries). The natural
   state is `to_ast`/`from_inner` (structs) or `to_ast` + variant tag (mirrors); the
   value types already round-trip, so `__reduce__` over that is mechanical.
4. **`AtomConstraints` container surface** *(done, 2026-07-08)* — reworked after doc 138
   settled: by-key `AtomConstraintKey` addressing (`get`/`contains`), Rust-mirroring named
   accessors (`valence()`… + `ring_count()`/`ring_size_count(n)`), and `asdict()` (canonical
   order; ring keyed `ring_count`/`ring_size_count_<n>`), with `AtomAst.asdict()` emitting
   the constraint dict. The type is now `AtomConstraintsAst` (the `*Ast` rename).
2. **Terse element notation `E.H` / `E.Cl`** (exploring; not needed now). A
   `python/umol/`-layer convenience mirroring the Rust `e!(H)` macro — a namespace
   object whose `__getattr__(symbol)` returns `Element(symbol)`; the only Python form
   that keeps the bare, unquoted symbol (`e("H")`/`E["H"]` need quotes; `e(H)` is a
   `NameError`). Works because every IUPAC symbol is a valid identifier. Tradeoff: the
   dynamic `__getattr__` has no IDE completion/type-checking without a `.pyi` stub
   (or generate the 118 as real constants for discoverability). No effect on the
   binding; decide the form with alpha users. *(done 2026-07-08, dynamic
   `__getattr__` form as `umol.E`; the `.pyi`-stub / 118-constants completion option
   remains open.)*
6. **Structural atom ops on `AtomViews` — `append` / `extend`** (possible; not settled).
   `mol.atoms.append(atom)` / `mol.atoms.extend(iterable_of_AtomAst)` add graph nodes,
   so they are *structural* — routed through `MoleculeBuilder` / `edit()` (one builder
   cycle; `extend` batches to a single cycle rather than one-per-atom), not the in-place
   `atom_mut` path. Appended atoms are isolated (bonding is a separate op). `extend`
   takes any iterable (Python idiom). `__delitem__` / remove is deferred — deleting an
   atom dangles its bonds (RDKit `RemoveAtom` semantics), needs the bond slice first.
   Together with `__setitem__` (whole-atom replace) this makes `AtomViews` a
   mutable-sequence facade.
