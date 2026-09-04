# Working in umol

umol is a Rust workspace for explicit molecular representations and operations, with Python
bindings in `umol-py`. Graph IR, graph views, geometry, and external formats are distinct models;
a format boundary is not the molecular model.

## Map

- `umol-graph-ir` defines the molecular and reaction DSLs and their semantic graph IR.
- `umol-graph-core` provides graph data structures and domain-independent graph
  algorithms.
- `umol-graph` provides chemistry-aware operations over the graph-IR model.
- `umol-io` owns external-format boundary representations and parsers.
- `umol-coordgen-sys` provides the vendored native boundary for feature-gated
  2D coordinate generation.
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

## Repository rules

- At the `umol-graph-core` and `umol-io` algorithm layers, callers select algorithms explicitly.
  Higher layers may use operation-specific config defaults; algorithm choice is not a chemistry
  parameter.
- Keep external formats behind boundary types. Conversion to graph IR may be lossy or model-dependent.
- Rust/Python boundary methods are `from_rust` and `to_rust`; prefixed Rust imports in `umol-py`
  drop `umol-`.
- Public names and visibility are design decisions. Do not hide unfinished design in `pub(crate)`
  helpers, indiscriminately re-export, or proliferate public seams.
- Before changing an invariant-bearing public API, enumerate its intended public symbols and match
  each constructor, conversion, visibility, and failure boundary to a settled decision. Reconcile
  every changed symbol before completing a staged subitem; passing tests do not override the contract.
- Closed or producer-issued types gain no arbitrary-parts, bytes, handles, or independent-context
  constructor unless explicitly settled. Public transformations must preserve their invariant.
- Property tests state laws; preserve the stated law when changing generators or assertions.
  Algorithm work starts with correctness fixtures, properties, and benchmarks. External comparisons
  are evidence, not runtime dependencies.
- Code states current behavior; discussion documents preserve reasoning and history.

Normative guides are `docs/development/data-types.md` (construction and fallibility),
`integrity.md` (representation contracts), `nomenclature.md` (terms and public names), and
`property-tests.md` (property suites). Discussion documents are non-normative and must not be cited
from source comments or public rustdoc.

## Working rules

- Be direct and compact: answer first; omit preambles, restatements, closings, filler, flattery, and
  repeated context. Correct errors plainly and retain user corrections for the session.
- Inspect files and code before claims or edits. Say what is unknown; never invent paths, symbols,
  signatures, or behavior. Own mistakes instead of deflecting to pre-existing or out-of-scope work.
- Keep scope exact: no speculative features, adjacent refactors during a fix, or unnecessary files.
- Never perform mutating git operations unless the user explicitly authorizes them; read-only git is
  allowed.
- Prefer direct code for one-off work. Before load-bearing design changes, discuss simplicity,
  generality, and correctness. Expose structural problems rather than hiding them in stubs, shims,
  bridges, or helpers; ask when the principled design is unclear.
- Do not design for backward compatibility; make breaking changes when technically necessary. Do not
  infer real scale from guesses or treat current tests/benchmarks as representative during design.
  Prefer a maintained library to manual reimplementation.
- When asked only for options, do not append unsolicited recommendations.
- Prefer names to explanatory comments. No long comments, self-talk, implementation history,
  decorative dividers, or stock values such as 42. Put imports at module scope and use `module.rs`,
  not `module/mod.rs`.
- In prose, use bare type/trait/constant names and parent-qualified free functions.
- Apply `ir-literal-extraction` before graph-IR literal extraction in `umol-graph` or higher crates.

## Authority and synchronization

- `AGENTS.md` is the sole repository instruction file; `CLAUDE.md` only points here.
- `.agents/skills` and `.claude/skills` expose the same names. A shared skill has one canonical
  body and a redirect in the other catalog; keep redirect frontmatter aligned.
- Current code and tests define implemented behavior.
- `discussion/000-status.md` defines discussion status. Proposed work is not implemented; a
  completed record covers only its completed scope. Basenames are at most 55 characters.
- `materials/` holds research inputs and reference implementations, not runtime data or ordinary
  fixtures.

## Session start

1. Read this file.
2. Run `git status --short`; preserve unrelated changes.
3. Read the status vocabulary and relevant rows in `discussion/000-status.md`.
4. Open the named or linked discussion document.
5. Inspect the current API, implementation, tests, and benchmarks before claims or edits.
6. Load applicable repository skills before implementation planning or test edits.
7. Before PyO3 builds or Python tests, activate `umol-py/.venv` and confirm Python 3.13.

Start from the status index and follow relevant links; never read the archive chronologically.

## Verification

- Format: `cargo fmt --all`
- Test: `cargo test --workspace`
- Lint: `cargo clippy --workspace --all-targets -- -D warnings`
- Python: activate `umol-py/.venv`, run `maturin develop`, then `pytest -q umol-py/tests`

Use narrow checks while iterating, then the plan's gate. Run feature-gated property and conformance
suites explicitly; the default workspace test does not cover them.
