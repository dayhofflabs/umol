# Working in the umol repository

## Repository purpose

umol is a Rust workspace for explicit, algorithmically transparent molecular
representations and operations, with Python bindings in `umol-py`.

- `umol-graph-ir` defines the molecular and reaction DSLs and their semantic graph IR.
- `umol-graph-core` provides graph data structures and domain-independent graph
  algorithms.
- `umol-graph` provides chemistry-aware operations over the graph-IR model.
- `umol-io` owns external-format boundary representations and parsers.
- `umol-geometric*` provides geometric representations and their connection to
  the graph model.
- `umol-py` exposes the supported high-level surface to Python.

The graph IR, graph views, geometric, and external-format representations are distinct
models. Do not treat a format boundary type as the one true molecular model.

## Crate map

| Crate | Responsibility |
| --- | --- |
| `umol-coordgen-sys` | Feature-gated vendored CoordGen source and native 2D-coordinate boundary |
| `umol-graph-ir`, `umol-graph-ir-macros` | Molecular and reaction DSLs, graph IR, constraints, edits, deltas, and validation vocabulary |
| `umol-chem` | Elements, isotopes, spin, occupation, units, and other chemistry vocabulary |
| `umol-edn`, `umol-edn-macros` | EDN parsing, formatting, and macros |
| `umol-graph-core` | Graph storage, relations, rewriting primitives, and graph algorithms |
| `umol-graph` | Chemistry models, resolution, validation, transformations, matching, and fingerprints |
| `umol-io` | SMILES and CTfile boundary objects, TableIR, parsing, and format configuration |
| `umol-geometric-core` | Coordinates, orientations, and geometric primitives |
| `umol-geometric` | Geometric molecular models |
| `umol-geometric-graph` | Conversion between geometric and graph molecular models |
| `umol-perm` | Permutations, cosets, and related stereo algebra |
| `umol-params` | Model-dependent chemistry parameters |
| `umol-msym`, `umol-msym-sys` | Molecular-symmetry wrapper and native interface |
| `umol-nauty-sys` | Vendored native foundation for nauty integration |
| `umol-py` | Python bindings and Python package |
| `umol-utils` | Small cross-cutting infrastructure |

## Architectural policies

- At the `umol-graph-core` and `umol-io` algorithm layers, callers select the
  algorithm explicitly. Do not introduce silent algorithm defaults.
- Higher-level operations in `umol-graph`, `umol-io`, and `umol-py` may define
  defaults in operation-specific config objects. Algorithm choices are
  operational config, not chemistry-model parameters.
- Keep external formats behind explicit boundary types. Conversion from a
  boundary representation to the graph IR may be lossy or model-dependent.
- Rust-to-Python and Python-to-Rust boundary methods are named `from_rust` and
  `to_rust`. Rust imports in `umol-py` use the crate name without the `umol-`
  prefix when an import prefix is needed.
- Public names and visibility are API decisions. Do not hide unfinished design
  behind `pub(crate)` helpers, add indiscriminate re-exports, or proliferate
  public helper layers.
- Property tests document semantic laws as well as checking implementations.
  Preserve the stated property when changing generators or assertions.
- Algorithm work includes correctness fixtures, property tests, and benchmarks
  from the beginning. Benchmarks and external comparisons are evidence, not
  hidden runtime dependencies.
- Code describes current behavior, not the history of how it was reached.
  Discussion documents preserve design reasoning.

## Authority and status

- Current code and tests are authoritative for implemented behavior.
- `discussion/000-status.md` is authoritative for the status of discussion
  documents.
- Discussion-document filenames must be at most 55 characters so the status
  table remains readable in source form. Prefer concise area names over
  sentence-like summaries.
- A completed discussion document records the completed scope; it is not a
  substitute for inspecting the current API.
- A proposed document describes future work and must not be reported as
  implemented.
- `materials/` contains research inputs and reference implementations. It is not
  a location for runtime data or ordinary checked-in test fixtures.

## Session startup

1. Read this file and `CLAUDE.md`.
2. Run `git status --short` and preserve unrelated user changes.
3. Read the status vocabulary and relevant entries in
   `discussion/000-status.md`.
4. Open the discussion document named by the user or linked from the relevant
   status entry.
5. Inspect the current public API, implementation, tests, and benchmarks before
   making claims or edits.
6. Load every applicable repository skill before planning implementation work
   or editing tests.
7. For Python work, activate `umol-py/.venv` and confirm that `python` resolves
   to Python 3.13 before building the PyO3 crate or running pytest.

Do not read the discussion archive chronologically. Start from the status index
and follow only the references relevant to the task.

## Task routing

| Work area | Start with |
| --- | --- |
| Graph IR, DSL, constraints, entities, reactions | `umol-graph-ir`; docs 113, 131, 132, 164, 165 |
| Graph algorithms and rewriting primitives | `umol-graph-core`; docs 136, 157-162, 167 |
| Chemistry-aware graph operations | `umol-graph`; docs 145-149, 166 |
| SMILES, MOL, SDF, and TableIR | `umol-io`; docs 151-153, 155 |
| Python bindings and workflows | `umol-py`; docs 137, 140, 150, 151 |
| Stereo and permutation work | `umol-graph-ir`, `umol-perm`; docs 103, 110, 157 |
| Geometry and molecular symmetry | `umol-geometric*`, `umol-msym`; docs 69-76 |
| Property-testing policy | doc 161 and the crate's property-test modules |
| Release preparation | doc 163 |

## Verification

- Format: `cargo fmt --all`
- Workspace tests: `cargo test --workspace`
- Workspace lint: `cargo clippy --workspace --all-targets -- -D warnings`
- Python tests: activate `umol-py/.venv`, run `maturin develop`, then
  `pytest -q umol-py/tests`
- Feature-gated property and conformance suites must be run explicitly; inspect
  the relevant crate's `Cargo.toml` and discussion plan for the required command.

Use the narrowest relevant check while iterating, then run the verification
gate specified by the implementation plan. Do not assume the default workspace
test command covers feature-gated suites.
