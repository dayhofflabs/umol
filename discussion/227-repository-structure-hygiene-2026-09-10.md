# 227 — Repository structure and hygiene at scale

Status: Proposed
Date: 2026-09-10
Relates: [117](117-entity-model-extensibility-2026-06-20.md),
[119](119-umol-perm-review-2026-06-21.md),
[code reviews guide](../docs/development/code-reviews.md),
[property tests guide](../docs/development/property-tests.md),
[development guides](../docs/development/README.md)

## Purpose

The workspace has grown past the size at which structure is carried by convention and by one
person's attention. This document records the current shape, what comparably sized Rust
workspaces do, and a set of proposed rules and measures for keeping the repository navigable,
searchable, and changeable. It is a scoping record. The rules are proposals until they are moved
into the living guides, and nothing here authorizes code or tooling changes.

Two concerns drive it. The practical one is readability: files over two thousand lines, mixed
types, traits, methods, and tests, and per-kind copies that must be changed together. The
process one is that every change currently passes through the same full gate, which makes the
maintainer the bottleneck and leaves agents loading normative text until their context is
exhausted.

## Current shape

Measured on the working tree at the document date. Line counts are `wc -l` over `.rs` files
outside `target/`.

| Measure | Value |
| --- | --- |
| Rust lines, workspace | 346,738 |
| `umol-graph-ir` (`src` plus `tests`) | 151,987 |
| `umol-py` | 51,946 |
| Files over 3,000 lines | 17 |
| Files over 2,000 lines | 33 |
| Files over 1,000 lines | 98 |
| Normative text loadable by an agent (`AGENTS.md`, `docs/development`, skills) | about 40,000 words |
| Discussion documents | 223 documents, about 945,000 words |

Two observations follow from the per-file split between code and inline test modules.

- Most of the largest files are large because of inline tests. `ir/reaction.rs` is 53 lines
  of code and 6,197 lines of tests; `umol-py/src/delta.rs` is 27 and 5,168. Of the ten largest
  files, seven are at least half tests.
- Four files have more than 2,500 lines of non-test code: `dsl/edit.rs`, `ir/delta.rs`,
  `umol-graph-core/src/relation.rs`, and `ir/reaction_span.rs`. These are the genuinely hard
  reads; the rest are a placement question.

`umol-graph-ir` is larger than every other crate combined, excluding the Python bindings. The
other eighteen crates are within ordinary Rust practice in size and shape. The internal layering
of `umol-graph-ir` between `ir`, `dsl`, `view`, `canonicalize`, and the edit, delta, and span
machinery is a directory convention, not a crate boundary, and nothing enforces it.

Doc 117 records the other structural cost: eight entity kinds with per-kind copies of storage,
views, constraints, DSL, deltas, edits, and bindings, and hand-modeled relations between kinds.
The cost of a fundamental change scales with the number of copies, not with line count. That
concern is addressed there, not here.

## What comparably sized Rust workspaces do

The lessons below are drawn from projects that are larger than umol and generally regarded as
navigable. Each is cited for one specific practice.

- **rust-analyzer.** Around forty crates layered as a strict dependency DAG, most under twenty
  thousand lines. Two short documents carry the structure: an architecture document of about
  two thousand words that says what each crate owns, which invariants cross crate boundaries,
  and where not to look; and a style document in which every rule has a one-paragraph rationale
  and a bad-versus-good example. Tests are mostly data: snapshot tests and test-data
  directories rather than inline test modules, which keeps source files readable.
- **Polars.** Data types and operations are separate crates: the frame and series types in one
  crate, the operations on them in another, the query plan and lazy layer in others. The split
  is by what a crate is for, not by size.
- **arrow-rs.** Crates are named by concern (`arrow-array`, `arrow-schema`, `arrow-cast`,
  `arrow-select`, `arrow-row`), each does one thing, and the name says which.
- **Cargo.** A single large crate kept navigable by module discipline alone (`core`, `ops`,
  `sources`, `util`) with all integration tests in one external test-suite directory. It shows
  that a single crate can work, given strict conventions and a full-time team.
- **rustc.** The counterexample. The crate graph is well layered and the files are enormous;
  the codebase is readable only to people who live in it. A 130,000-line crate drifts toward
  this shape.

The shared practices are: the crate is the unit of API boundary and the file is the unit of
reading; every crate and module has a short header saying what it owns; every rule that can be
mechanical is mechanical, so the prose guide covers only what tools cannot; and the navigational
document is short while reference material is looked up rather than loaded.

## Why prose rules decay under generated code

A language model conditions on the code it reads. Every file in the neighborhood of an edit is
evidence about what this codebase does, and a rule in a guide is one further piece of evidence
among many. When fifty files share a deviation, the deviation is the convention as far as the
model is concerned, and the guide loses. The repository's memory of corrections, which are
overwhelmingly corrections of copied deviations, is the record of this effect.

Two consequences shape the rules below. Structure has to live in the code the model reads, in
the form of correct neighbors, more than in the documents it might load. And any rule that can
be checked by a tool should be checked by a tool and removed from the agent's context, because a
build failure is evidence the model cannot argue with.

## Proposed rules

These are candidates for `docs/development/`, stated as rules so that each can be accepted,
rejected, or reworded on its own.

1. **The crate is the boundary unit; the file is the reading unit.** A file holds one concept:
   a type or a small family of types with their inherent impls and trait impls. Files over a
   fixed size, measured on non-test lines, are split by concept, never by line count.
2. **Every module opens with a header.** A `//!` block of a few lines stating what the module
   owns, what it deliberately does not own, and where to go instead. It is the first thing a
   reader or a model sees.
3. **Tests beyond a threshold live beside the code, not inside it.** An inline `#[cfg(test)]`
   module over the threshold moves to a sibling `tests.rs` or a `tests/` directory under the
   module. The threshold is a number to settle; the current split is three test submodules
   against eighty-four inline modules in `umol-graph-ir` alone.
4. **A rule that a tool can check is checked by the tool and leaves the prose.** Formatting,
   forbidden identifiers, file length, test placement, module headers, and dependency direction
   are all checkable. The guide retains only what needs judgment.
5. **The navigational document is short; reference documents are looked up.** One architecture
   document of about two thousand words is loaded at session start. The nomenclature glossary,
   the integrity inventory, and the data-type contracts are consulted by search on demand and
   are not loaded whole.
6. **Exemplars over rules.** For each recurring module shape there is one named exemplar module,
   and the rule for a new module of that shape is to match the exemplar. A guide that says
   "match `x.rs`" is more reliably followed than a guide that describes `x.rs`.
7. **Regenerate rather than patch.** Generated code is cheap to regenerate and expensive to
   patch. The durable assets are the design documents, the public type surfaces, and the tests
   and laws. A file that has accumulated several generations of patches is rewritten from those
   assets, not patched again.
8. **Changes pass a gate proportional to their kind.** Mechanical, tool-verified, test-preserving
   changes need no design document and no maintainer review. New public surfaces need a design
   document and sign-off on the type signatures only. Semantic changes keep the full current
   process: design, naming, plan approval, review cycle.

## Proposed measures

Each measure is independent of the others except where stated. None is planned; ordering is
noted only where one measure makes another cheaper.

### Test relocation

Move inline test modules over the threshold into sibling test files across the workspace. This
halves the largest files with no semantic change and is the natural first step because tests are
the asset every later step depends on. A CI check keeps it done.

### Mechanical checks

Candidate checks, all currently prose rules or memory items:

- non-test line count per file, with an explicit allowlist carrying a one-line reason per entry;
- forbidden identifier fragments (`n_`, `num_`, `slot`, `verdict`, and others recorded in the
  nomenclature guide's retired list);
- inline test module size;
- presence of a module header;
- `module.rs` rather than `module/mod.rs`;
- crate dependency direction.

Three ways to implement them, with the trade-offs:

- **Scripts under `scripts/` run from CI.** No toolchain coupling, trivially portable, cannot
  see the syntax tree, so identifier checks are textual and file checks are line-based.
- **Clippy configuration.** `disallowed-names`, `disallowed-methods`, and `disallowed-types` in
  `clippy.toml` cover identifier and API bans with real name resolution and no custom code. They
  do not cover file length, test placement, or headers.
- **dylint.** Custom lints with full syntax-tree and type access, so every check above is
  expressible precisely. A prior experiment found the setup costly; the toolchain pinning and the
  separate build are ongoing maintenance. It is the only option that can express a check such as
  "every trait in this family is implemented for every entity kind".

The options are not exclusive; a script layer and a clippy layer can carry most checks, with
dylint reserved for the checks the other two cannot express.

### Splitting `umol-graph-ir`

The existing directories are the candidate seams: the graph IR data model, the DSL, the edit,
delta, transaction, and span machinery, canonicalization, and stereo. A crate boundary enforces
layering that a directory only suggests, and per-crate test runs become fast again. The split
interacts with doc 117: if the per-kind copies collapse, the crates that hold them become
smaller before they are split. Which of the two happens first is open.

### A navigational architecture document

Either the `Map` and `Crate map` sections of `AGENTS.md` grow into it, or a new
`docs/development/architecture.md` is created and `AGENTS.md` points to it. It states, per
crate, the purpose, the key types, the invariants that cross its boundary, what depends on it,
and where its tests live. It is the one document loaded at session start.

### Reference documents as lookup

The nomenclature glossary at fourteen thousand words and the data-type contracts at ten
thousand are reference material. They stay where they are, but the skills and `AGENTS.md` stop
instructing agents to read them whole; the instruction becomes to search them for the term or
contract at hand. Closed discussion documents whose settled rules have been moved into the guides
are not read again.

### Snapshot tests

For DSL round trips and canonicalization, snapshot tests would shrink the test corpus by a large
factor, and rust-analyzer's experience is that they keep source files readable. This conflicts
with the repository rule that test cases use direct literals inline, and with the property-test
policy that tests state laws. Whether snapshot files are a third test form, a replacement for
some literal-heavy suites, or excluded, is open.

## Open questions

- The file-size threshold and the inline-test threshold, as numbers.
- Whether the tiered gate is recorded in `AGENTS.md`, in the code reviews guide, or in a new
  process guide.
- Which checks go to scripts, which to `clippy.toml`, and whether dylint is taken up at all.
- The order of test relocation, the doc 117 collapse, and the crate split.
- Whether snapshot tests are admitted, and under what rule.
- Names: the architecture document, the exemplar modules, and the check scripts are unnamed.
