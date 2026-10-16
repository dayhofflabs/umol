# 227 — Repository structure and hygiene at scale

Status: Proposed
Date: 2026-09-10
Relates: [117](117-entity-model-extensibility-2026-06-20.md),
[119](119-umol-perm-review-2026-06-21.md),
[166](166-molecule-ops-2026-07-27.md),
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

## Relation storage split experiment — 2026-09-15

A scratch experiment split the five relation-set implementations into private
child modules, retained shared vocabulary and implementation support in the
parent, and relocated the existing unit tests into separate files. The experiment
left the main workspace's production and test files unchanged. The proposed naming
and participant replacement work belongs to doc 166 and was not implemented here.

The experiment used relation.rs at revision
8a97660ed9e21238a86191de53ef23e41ef62e8d (SHA-256
020dd4ee3aceeb8e3b6e393c22d300098147c8c0b476510db91fef1992f8c5b6).
Its measurements, visibility findings, preservation results, and limitations are
recorded below. The scratch candidate, scripts, and logs were disposable; the
production implementation and its tests are retained in the repository.

### Sizes and dependencies

The original relation.rs is 5,585 lines: 2,819 before its inline test module and
2,766 in the test section. Candidate counts include comments and blank lines:

| Component | Production lines | Unit-test lines |
| --- | ---: | ---: |
| Shared relation parent | 353 | 65 support + 96 participant tests |
| FixedRelationSet | 421 | 471 |
| VarRelationSet | 450 | 580 |
| FixedFixedBirelationSet | 528 | 477 |
| FixedVarBirelationSet | 567 | 595 |
| VarVarBirelationSet | 566 | 499 |
| Total | 2,885 | 2,783 |

The total grows by 83 lines of module structure, imports, headers, and spacing.
No resulting file exceeds 595 lines. Each storage type retains its inherent and
trait implementations in one file. The five implementations use shared parent
helpers and existing graph/compaction/correspondence/remapping modules; none
calls another storage implementation. Four types refer to FixedRelationSet's
algebra documentation, which is a documentation dependency rather than a code
dependency. Explicit crate-root rustdoc targets preserve those links.

### Visibility

The experiment requires no new pub(crate), pub(super), or pub(in ...) markers
and no newly public helper. Private parent items are accessible to all child
modules. Existing public storage types are re-exported through the parent and
the existing crate-root API. Shared test support remains private in its own
parent module. The candidate also changes the existing crate-root declaration
from pub(crate) mod relation to mod relation; compilation confirms that the
scoped marker is unnecessary there.

Placing Incidence in a separate sibling module would instead require exposing
its type and methods to that sibling's consumers. Compiler probes confirm that
plain pub inside a private module works without scoped markers, but the shared
parent avoids even those additional public declarations. At 353 lines, it does
not need another split merely to reduce file size. This experiment supports
choosing the module hierarchy before adding visibility annotations.

### Ordering, documentation, and property evidence

The relocation groups tests by type but preserves existing method order and
test bodies. Before adopting these files as exemplars, settle a common method
order that distinguishes collection queries, per-relation operations, mutation,
id transport, and algebra. Keep checked/asserted and plain/tracked peers adjacent
and align test order with the methods. The current family differs in reader/
permutation placement, has scattered test order, and implements Hash for only
three of its five types. Those are explicit review items, not changes made by
the experiment.

Documentation also needs review independent of relocation. Existing comments
describe removed ordering markers, sorting during construction, canonical
entries, and uniqueness not established by constructors. Indexed-access panic
preconditions are mostly undocumented; permutation describes its panic in prose;
map/remap have Panics sections. Option-returning operations need clear absence
conditions and callback preconditions rather than invented Result errors.
Semantic-properties sections currently cover transport but do not collect all
the existing permutation, compaction, and algebra contracts coherently.

The external relation property module remains unchanged. Its ten properties
cover transport composition/inverse/remap agreement and exact-size iteration
across all five shapes. Its header should point back to the public properties
when that documentation is settled. Public operation-family documentation can
identify the property target; the executable target remains the inventory, as
required by the property-tests guide. File splitting is not evidence for a new
semantic law or complete property coverage.

### Validation and conclusion

The baseline and candidate each pass the same 220 expanded relation unit cases.
The candidate passes all ten relation properties, all-target graph-core
compilation with the proptest feature, and rustdoc with private items and
warnings denied. Preservation checks compare all 134 unit-test functions and
their attributes after formatting, and the five storage bodies apart from
whitespace and rustdoc link targets. No storage fields or implementation helpers
gain visibility.

This establishes a feasible module-and-test layout with ordinary public/private
visibility. It does not measure performance, approve the final method ordering,
or implement the relation API additions. A production adoption must reconcile
the documentation and naming decisions with doc 166; the scratch split is not
itself a completed repository cleanup.

### Agreed organization — 2026-09-15

Keep the five storage implementations in private child modules, with public
symbols explicitly re-exported through relation.rs. Extract the public
participant vocabulary into participant.rs: RelationParticipant, ParticipantRefs,
and ParticipantPosition remain pub and are re-exported by the parent. This
extraction needs no scoped visibility.

An incidence.rs module is appropriate if the forthcoming participant-mutation
machinery warrants a separate internal component. Its Incidence type and the
methods used by sibling storage modules may be pub(super), with fields and
implementation details private. This preserves the current access boundary:
a private item in relation.rs and a pub(super) item in its immediate child are
both accessible within the relation subtree. Keep small shared functions in
the parent until they form another coherent component; the parent need not
contain only re-exports.

Use visibility to state those roles consistently:

- pub for intended public API, explicitly re-exported;
- private for implementation owned by a module and its descendants;
- pub(super) for an internal component interface shared within the relation
  subsystem;
- pub(crate) only for an interface actually needed beyond that subsystem.

This is the accepted direction for the relation family, not a new workspace-wide
visibility rule. The experiment establishes that avoiding additional markers is
feasible; minimizing their count is not the final organization criterion. Its
measurements describe the tested parent-owned layout, not a subsequently
implemented participant/incidence extraction.

For production adoption, keep each unit-test file self-contained, with local
fixtures, test payload types, assertion helpers, and explicit imports of production
symbols and external test tools. Test modules must not import from one another.
Small duplication is acceptable; constructor shorthands such as n are removed
in favor of direct construction. The parent
tests.rs only declares child modules; no shared fixtures/utils module is planned.
This revises the experiment's shared test-support arrangement, not its recorded
measurements or results.

Doc 166's relation-storage S0–S4 plan sequences production adoption of this
layout, contract documentation, participant mutation, and the approved naming
migration. It records fresh size and visibility measurements after those
changes; the wider structural program here remains separately proposed.

Both arrangements have established precedents.
[Arrow's array module](https://github.com/apache/arrow-rs/blob/main/arrow-array/src/array/mod.rs)
defines its shared Array trait alongside concrete implementation modules and
re-exports. [Tokio's mpsc module](https://github.com/tokio-rs/tokio/blob/master/tokio/src/sync/mpsc/mod.rs)
uses sibling implementation modules and collected exports, with scoped
visibility on [shared channel machinery](https://github.com/tokio-rs/tokio/blob/master/tokio/src/sync/mpsc/chan.rs).
The [Rust visibility rules](https://doc.rust-lang.org/reference/visibility-and-privacy.html)
explain the access boundary preserved by a parent-private to child-pub(super)
extraction.

### Production adoption — completed 2026-09-17

Doc 166's S0–S4 relation-storage prerequisite is complete. Its
[closeout](166-molecule-ops-2026-07-27.md#s4--verify-and-record-the-delivered-prerequisite)
records the final sizes, API/delegation audit, dependencies, performance evidence,
and verification. This completes the relation pilot; the workspace-wide proposals
in this document remain Proposed.

The delivered parent is 200 lines, participant vocabulary 143, and Incidence 91.
The five storage implementations are 584–935 lines; their unit-test files are
858–1,655 lines, with 106 lines of participant tests and six declarations in
tests.rs. Production totals 4,283 lines and unit bodies/support 6,646. External
relation properties remain one 2,190-line module and benchmarks one 1,685-line
file. These counts include the complete mutation family and expanded contracts
and tests; the original experiment's counts describe its earlier snapshot.

The final implementation follows the agreed private-child organization and
explicit public re-exports. The only scoped visibility is six pub(super)
declarations for Incidence's type and five methods; fields and parent helpers
are private. Storage implementations share vocabulary/index/parent functions
without depending on sibling storage implementations. Their algebra rustdoc
is now self-contained too. Test files are independent and own their fixtures
and support; there is no shared test-support module or constructor shorthand.

Method order, failure contracts, semantic properties, and references to public-API
properties were reconciled during adoption. The original functions, exact cases,
and property laws are preserved apart from approved renames and relocation.
The final workspace gates passed, and 1,926 benchmark cases ran against the
delivered surface. The S1 comparison records both regressions and improvements;
the pilot does not establish unchanged performance or isolate their cause.

No workspace-wide threshold, exemplar policy, mechanical checker, tiered gate,
or further crate/module split is settled by this pilot. Those decisions and the
ordering of the wider program remain the next work here. Graph-IR/editor wiring
and hydrogen transformations remain follow-ups owned by doc 166.

## Open questions

- The file-size threshold and the inline-test threshold, as numbers.
- Whether the tiered gate is recorded in `AGENTS.md`, in the code reviews guide, or in a new
  process guide.
- Which checks go to scripts, which to `clippy.toml`, and whether dylint is taken up at all.
- The order of test relocation, the doc 117 collapse, and the crate split.
- Whether snapshot tests are admitted, and under what rule.
- Names: the architecture document, the exemplar modules, and the check scripts are unnamed.
