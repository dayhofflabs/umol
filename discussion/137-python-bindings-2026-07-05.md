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
- **Reads:** mirror the Rust accessor shape (`atom.charge()` → `atom.charge`);
  computed/fallible stay methods (`mol.fingerprint()`).
- **Dunders:** `__repr__` via Rust `Display`, `__eq__`/`__hash__` via Rust
  `PartialEq`/`Hash` where the type has them; unsupported comparisons raise, not
  silently return `NotImplemented`-by-omission.
- **Mutation:** on the owner facade or a builder (see the facade section), never
  through child view objects.

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

**Mechanism (spike-decided, not chosen in the abstract).** Two ways to realize the
per-variant-class surface, both yielding the identical Python API:

- **Native PyO3 complex enums** where they compile — PyO3 generates exactly the
  variant-class + `__match_args__` shape for free. The only blocker is the
  unverified recursion in `ValueTerm`/`ValuePredicate` (`Box`/`Vec` self-reference).
- **Pure-Python variant classes over an opaque Rust core** (polars-style) for
  whatever native enums can't take — the `python/umol/` layer defines the class,
  converting to/from the Rust value at the boundary.

The whole choice hinges on one unknown: **can a PyO3 native complex enum hold a
recursive `Py<Self>` / `Vec<Py<Self>>` field?** That is an early compile spike on
`ValueTerm` (below). If yes: native enums throughout (cleanest, most single-source).
If no: native enums for the non-recursive types, pure-Python variant classes over
the opaque value for the recursive trees only. Sequencing within the slice:
`Element` → the scalar-field enums (`ElementAst`, `IsotopeMassAst`, `ValueAst` +
`ValueTerm`/`ValuePredicate`, `SpinStateAst`) → `AtomConstraints` (13 variants +
sub-ASTs) last.

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
- **Mirror mechanism** — spike-decided: PyO3 native complex enums where they compile,
  pure-Python variant classes over the opaque Rust value where they can't. Hinges on
  one compile spike — whether a native complex enum can hold a recursive
  `Py<Self>`/`Vec<Py<Self>>` field (`ValueTerm`/`ValuePredicate`).
- **Borrowed derivations** — `ReactionDerivation<'a>` gets an owned Rust result type
  the binding wraps (Python-pull change to the Rust API).
- **Generics** — monomorphize at the boundary: per-entity-family `EntitySpan<T>`
  types; only the surfaced coloring(s) for `GraphSymmetryConfig<C>`.
- **Errors** — Python exception hierarchy mirrors the three tiers; structured fields
  carried onto exception objects; typed exceptions want tier-2 enums, `Box<dyn>`
  only for the generic base.
- **Surface policy** — native layer keeps Rust names, Python layer is PEP8; getters
  for scalar fields, methods for computed/fallible; dunders via Rust
  `Display`/`PartialEq`/`Hash`. Recorded here, harvested into the ubook later.
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

### S0 — Scaffold + leaves (green)

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

### S1 — Mechanism spike + value-expression core (green)

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

### S2 — Atom-field value types (green)

- **S2a** — `ElementAst` (5 variants; `Lit` uses `Element`, `Var` uses `MemOp`).
  *Additive.* [dep: S0b, S1b]
- **S2b** — `IsotopeMassAst` (5 variants incl. `Natural`). *Additive.* [dep: S0a]
- **S2c** — `SpinStateAst` (`{unpaired, multiplicity}`; fields are `ValueAst`).
  *Additive.* [dep: S1d]

### S3 — `AtomAst` owned (green)

- **S3a** — `AtomAst` pyclass (module `atom`): mirror-typed fields (element /
  isotope_mass / charge / implicit_hydrogens / lone_pairs / spin); Rust-parallel
  constructors (`new`, `from_element`, `with_*`); getters; `__repr__`/`__eq__`.
  `constraints` deferred to S5c. *Additive.* [dep: S1d, S2a, S2b, S2c]

### S4 — Molecule atom read surface (green)

- **S4a** — `AtomId` newtype pyclass + `AtomViews` / `AtomView` handle pyclasses
  (module `atom`): views hold `Py<MoleculeAst>` + `AtomId`; `AtomView.id` is an
  `AtomId`, molecule indexing takes one; getters return the mirror types by reading a
  transient Rust `AtomView`. Bond / topology-derived methods (`valence`, `neighbors`,
  …) are out of scope (need the bond slice). *Additive.* [dep: S0c, S3a]
- **S4b** — enrich the `MoleculeAst` pyclass: `.atoms` → `AtomViews`, indexing /
  iteration → `AtomView`, and `from_atoms_and_bonds`-style construction taking Python
  `AtomAst`s. *Additive.* [dep: S4a, S3a]

### S5 — `AtomConstraints` (green)

- **S5a** — constraint sub-ASTs (module `constraint`): `RingScope`,
  `RingMembershipAst`, `AromaticValenceAst`, `MulticenterValenceAst`,
  `TetrahedralStereoAst` (the last may pull one or two more small stereo types).
  *Additive.* [dep: S1d]
- **S5b** — `AtomConstraint` (13 variants) + `AtomConstraints` container (add / iter /
  get-by-key, mirroring the Rust API). *Additive.* [dep: S5a, S1d]
- **S5c** — wire a `constraints` getter onto `AtomAst` and `AtomView`. *Additive.*
  [dep: S5b, S3a, S4a]

### S6 (step ii) — atom DSL parsing (green)

- **S6a** — `AtomAst.parse(str)` / `str(atom)` via `AtomDsl` (`FromStr`/`Display`);
  full roundtrip incl. constraints. Parse errors → a minimal binding exception (the
  full three-tier hierarchy is a later slice). *Additive.* [dep: S3a, S5c]

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
  on alpha-user feedback.
