# 180 — The top-level `umol` facade crate

Status: Proposed
Date: 2026-08-02
Relates: [163](163-package-release-preparations-2026-07-26.md),
[176](176-ast-naming-2026-07-31.md)

The workspace has eighteen members and none of them is `umol`. A user who reads the whitepaper and
runs `cargo add umol` currently gets the `0.0.0` placeholder published on 2026-08-02, which contains
nothing. This document scopes the crate that has to replace it.

Distinct from 163, which covers versioning, CI, and the mechanics of publishing the crates that
already exist. This is about a crate that does not exist yet and about what its public surface should
be.

## Justification

**The placeholder created an obligation.** Reserving a name and leaving it empty is the practice this
project set out not to imitate. The reservation is defensible only if it is redeemed.

**Nothing currently defines the public API.** `umol-ast` has one `pub use`; `umol-graph` and
`umol-io` have none. What is public is whatever each member crate happens to expose, and the only
curated view of the library that exists anywhere is `umol-py`. A facade makes the question "what is
the public API" answerable in one file.

**It is the documentation entry point.** `docs.rs/umol` is where a reader arriving from the
whitepaper's Availability section lands. Today that would be a placeholder.

## Framing first: model classes, not a graph library with accessories

**Read [015](015-four-domains-semantic-model-2025-03-10.md) before this section.** Its 2026-07-03
note is the governing statement: model classes are represented by separate crates (`umol-graph`,
`umol-geometric`), model parameters by explicit objects (`ChemistryModel`), and conversion between
instances of different model classes by **bridge crates** (`umol-geometric-graph`), fallible and
one-to-many in either direction.

So `umol-geometric` is a **peer model class**, not a subsystem, and `umol-geometric-graph` is a
**bridge**, an instance of the Model Domain's conversion relation. That it has received less
attention than `umol-graph` is a fact about time, not about architecture (author, 2026-08-02).

Two readings of \umol\ are in tension, deliberately. Under the first, \umol\ is a cheminformatics
graph library with `umol-graph` at the centre and everything else supporting it. Under the second —
the original concept, and the one 015 describes — `umol-graph` is one model class among others,
alongside the BO-type model in `umol-geometric` and potential coarse-grained or electronic-structure
models later. The whitepaper gestures at the first, appropriately for 0.6, and its title says
"*with* the \umol\ library" precisely so that it makes no claim to be the whole of it.

**The facade must not foreclose the second reading.** Concretely: name features by model class rather
than by "core plus options", let module structure read as peers with bridges exposed as conversions,
and keep the shared substrate implicit rather than presenting it as the graph model's foundation.

**Settled (author, 2026-08-02): `umol-ast` is the IR for the graph model, no more and no less.** It
is not shared substrate. Consequences for the facade:

- `ast` is not a root-level peer of `graph`; it belongs inside the graph model's namespace.
- The genuinely shared layer is thinner than the crate graph suggests — `umol-chem`, `umol-edn`,
  `umol-perm`, `umol-utils`.
- **The prelude is model-scoped, not universal.** `Lattice`, `Canonicalize` and the entity ASTs are
  graph-model vocabulary. A model-neutral prelude would be nearly empty, so the natural shape is a
  prelude per model class with `umol::prelude` re-exporting the default model's.

**`umol-io` is shared** (author, 2026-08-02): MOL encodes graphs through atom and bond blocks but
also carries coordinates, making it the equivalent of XYZ and other formats. I/O is not graph I/O.

**The seam for this already exists.** Of the sixteen files in `umol-io/src/table_ir/`, fifteen
contain no reference to `umol_ast`; only `raise.rs` does. The table IR — atom and bond blocks,
coordinates, properties, sgroups, spans — is already format-level and model-neutral. `raise` is the
graph model's lift out of it.

**Direction of the dependency, settled (author, 2026-08-02).** `umol-graph` becomes the top-level
orchestrator: DSL parse and render via the AST, SMILES/MOL/SDF ingest and output, resolution,
validation, transformation. `umol-io` stays below it.

The decisive argument is XYZ. If `umol-io` sat *above* `umol-graph` it would also have to sit above
every other model class in order to serve coordinate formats, which is the monster-crate outcome.
Below everything, emitting table IR, one parser per format serves all model classes. The intuition
that "you do not read SMILES unless you want to do something with it" is about usage, and usage is a
weak reason to invert a dependency arrow.

**Cost of io being transitively in the default**, since `graph` is default and depends on it:
exactly one dependency. `umol-ast` and `umol-edn` use **winnow**; `umol-io` uses **nom** across 18
files; `regex` and `serde` are already pulled in by `umol-graph`. So io's marginal contribution is
`nom` alone — and the planned nom-to-winnow migration removes even that, after which io adds no new
dependency to the default build at all.

Alternative, if io in the default ever becomes a real objection: split the boundary type into
`umol-table-ir` so that `umol-io` holds the parsers and `umol-graph` depends only on the type. That
would make `io` a genuine feature. Judged not worth a further crate, and worth still less after the
winnow migration.

**Confirmed (author, 2026-08-02): `raise` does not touch the parser, and the table IR is the full
boundary type.** The move is a relocation, not a refactor.

**Recommended: move `raise` to the graph model class; leave `umol-io` unconditionally shared.**
Parsing yields table IR; each model class raises from it. The graph model lifts atom and bond blocks
into `umol-ast`; a geometric model lifts coordinates into a conformer. One CTfile parser, two lifts,
nothing duplicated and nothing straddling.

This puts `umol-io` *below* both model classes rather than above them, and it needs no feature gates
of its own. Convenience entry points such as `from_smiles` become compositions of parse-then-raise
and belong in the model class or in the facade.

Considered and rejected:

- **Splitting `umol-io` internally along feature lines.** Mechanically fine — `cfg(all(feature =
  "graph", feature = "io"))` works and features forward into dependencies — but it produces a crate
  whose public API changes shape with features, a combination matrix to test, and one docs.rs page
  that can only show one configuration.
- **Per-model I/O crates (`umol-graph-io`, `umol-geometric-io`).** Each would need the CTfile parser,
  so either it is duplicated or a shared core reappears — at which point this is the recommended
  option with extra crates.

**`umol-params` is shared** (author, 2026-08-02): covalent radii and PPP parameters, "a mixed bag,
generally shared." Consumed by `umol-graph` and `umol-geometric-graph`.

**Two distinct things are called parameters, and the facade must not conflate them.**

| | what | where |
| --- | --- | --- |
| empirical data | covalent radii, quantum/PPP parameters | `umol-params`, shared layer |
| model configuration | `ChemistryModel`, `ValenceTable`, `AtomTypeRegistry` | `umol-graph`, inside the graph model class |

The second is what [015](015-four-domains-semantic-model-2025-03-10.md) means by "the specific model
parameters by explicit model objects (ChemistryModel in umol-graph)". A geometric model class would
have its own, and they are not the same kind of object. Expose the data in the shared layer and the
model objects under `umol::graph`.

Note for the whitepaper: \Cref{sec:primer} imports `ChemistryModel`, `ValenceModel` and
`ValenceTable` from a flat `umol` namespace, which is correct for Python but becomes graph-scoped in
Rust.

`umol-params/quantum` (PPP) sits beside the existing `AromaticityModel::Hmo` and points toward an
electronic-structure model class. Under the four-domains reading it is less a mixed bag than an early
seed of a second model class.

Minor, confirmed stale (author, 2026-08-02): `umol-geometric` declares `umol-params` but no source
file references `umol_params::`. The covalent radii were used there for bond perception, which since
moved to `umol-geometric-graph` — correctly, since perceiving bonds from coordinates is a conversion
between model-class instances and therefore bridge work. The dependency line is residue; drop it.

## What the dependency graph decides

**These are observations about the code as it stands, not constraints on the design** (author,
2026-08-02: current behavior can be changed and must not drive it). Recorded so the implementer knows
the starting point. Verified 2026-08-02.

**Nauty stays, and this is already decided elsewhere.** `umol-ast` depends on `umol-graph-core`,
which depends on `umol-nauty-sys`. Do not reopen it here:

- [137](137-python-bindings-2026-07-05.md) decided to keep nauty with no feature gate — automorphism
  and symmetry are core, and the local cost is a cached one-time ~5s build. Distribution is solved by
  prebuilt abi3 wheels from CI rather than by removing the dependency. The gate remains available for
  the one case out of scope: a source build on a platform no wheel covers.
- [143](143-vendored-nauty-integration-2026-07-11.md) already reduced the burden substantially, moving
  off `nauty-Traces-sys` (bindgen, LLVM/libclang) to vendored C sources compiled with `cc` behind a
  handwritten FFI shim. The requirement is now a C compiler and nothing else.
- [121](121-linux-profiling-setup-2026-06-21.md) records the residual portability question as tracked
  future work.

Two facts worth carrying into the facade. The coupling is **shallow**: nauty appears in one file,
`umol-graph-core/src/algorithms/automorphism.rs`, with one import, and `umol-ast`'s 51 files of
graph-core usage take only types and traits (`NodeId`, `RelationData`, `ParticipantPosition`). And
the **blast radius is narrower than it looks**: nauty is reached only at runtime through
`AutomorphismAlgorithm::Nauty` and never by matching, so a hypothetical no-nauty build keeps
substructure search and loses canonical labeling. `AutomorphismAlgorithm` currently has that one
variant; 137 notes pure-Rust individualization-refinement crates (`canonical-form`,
`graph_symmetry`/CNAP) as a feasible fallback should one ever be wanted.

**I/O is separable, but only after ingest moves.** `umol-graph` depends on `umol-io` because
`umol-graph/src/ingest.rs` and `parse.rs` hold the SMILES and CTfile entry points. The edge points the
wrong way: the algebra should not depend on file formats. Relocating ingest — as `umol-ingest`, or
into `umol-io` with io depending on graph — makes the layering monotone and incidentally makes `io`
gateable, taking `umol-geometric-core` with it since `umol-io` is its only consumer.

The feature payoff alone is modest: the beneficiaries are DSL-only users, the EDN path survives
through `umol-ast` and `umol-edn` regardless, and `io` would be default-on. **Do it for the layering,
and treat the flag as a side effect** (author, 2026-08-02: relocating ingest is acceptable).

**One model class is separable today, without any restructuring.** `umol-geometric`,
`umol-geometric-graph`, `umol-msym` and `umol-msym-sys` are reachable only through `umol-geometric`.
This reflects how far each model class has been developed, not a claim that one is subordinate.

**Both `-sys` crates vendor their C sources** and build through `cc`. There is no system library, no
`pkg-config`, no package manager step. `cargo add umol` needs a C compiler and nothing else. This is
worth protecting — it is the difference between a library people try and one they bounce off.

**Caution on a name.** `umol-geometric-core` is *not* the core of `umol-geometric`. It is a leaf with
no dependencies, used only by `umol-io` for coordinates, whereas `umol-geometric` handles symmetry
via `umol-msym`. Nothing depends on both. The naming implies a relationship the graph does not have;
worth an entry in [177](177-nomenclature-guide-2026-07-31.md) or a rename.

## Scope

### Feature flags

**Two, named for model classes.** Settled 2026-08-02, after the placement questions above.

| feature | default | contents |
| --- | --- | --- |
| `graph` | **on** | the graph model class: `umol-graph-core`, `umol-nauty-sys`, `umol-ast`, `umol-ast-macros`, `umol-graph`, and `raise` once relocated from `umol-io` |
| `geometric` | off | the geometric model class: `umol-geometric`, `umol-msym`, `umol-msym-sys` |

Everything else is shared and ungated: `umol-chem`, `umol-edn`, `umol-edn-macros`, `umol-perm`,
`umol-utils`, `umol-params`, `umol-io`, `umol-geometric-core`.

**No `io` feature.** Once `raise` moves out, I/O is model-neutral and small. Gating it would break
every whitepaper listing for a compile-time saving that does not justify it. This supersedes the
earlier suggestion in this document that relocating ingest would make `io` gateable — relocating it
is still worth doing, but for the layering, not for a flag.

**No nauty feature.** Settled in [137](137-python-bindings-2026-07-05.md); see above.

#### What the default set does and does not assert

Recorded 2026-08-02, because the question is easy to reopen badly.

**Measured:** a clean build of the whole geometric stack including its C sources —
`cargo clean -p umol-msym-sys -p umol-msym -p umol-geometric` then `cargo build -p umol-geometric` —
takes **5.43s**, the same order as the ~5s that [137](137-python-bindings-2026-07-05.md) quotes for
nauty. So compile time is not an argument in this decision, in either direction.

**Two commitments pull opposite ways.** Model classes are peers ([015](015-four-domains-semantic-model-2025-03-10.md)),
which argues for enabling both by default. And nothing should be promised that cannot be delivered,
which argues against putting a rudimentary geometric surface in front of everyone who runs
`cargo add umol`. Only the second is legible to a user: `default = ["graph"]` reads as "geometric is
opt-in", not as "geometric is subordinate", to anyone who has not read 015. A thin API enabled by
default draws a conclusion immediately and needs no interpretive frame.

**The decisive reframing, and the real question** (author, 2026-08-02: the question is not whether
`geometric` is excluded but whether `graph` should be included). **Peerhood is a property of the
feature graph, not of the default set.** What makes two model classes peers is that either can be
switched off and the remainder still coheres. What would demote `geometric` is `graph` being
*mandatory*, not `graph` being convenient. So the test is:

```toml
umol = { version = "0.6", default-features = false, features = ["geometric"] }
```

**Today this would not compile into anything coherent**, because `umol-io`'s `raise` ties I/O to
`umol-ast`: there is no I/O without the graph model's IR. After the relocation recommended above,
it does — parsing yields table IR, `umol-chem` supplies elements, `umol-params` supplies radii and
PPP parameters, and a geometric model raises coordinates from the same table IR. Thin, but coherent.

**Therefore: relocating `raise` is what makes `graph` genuinely optional, and so what makes the
peer claim true rather than asserted.** That is a stronger reason to do it than the layering
argument, and it should be treated as a prerequisite for the equal-model reading rather than as
cleanup.

**Decision:** `default = ["graph"]` as a convenience for who shows up today, *conditional on* `graph`
being genuinely optional. The default set states what is common; `default-features = false` states
what \umol\ is. Document that distinction in the crate docs, since the manifest cannot express
"peer but early".

Reversibility, for whenever this is reopened: adding a feature to `default` is additive and
non-breaking, removing one is breaking. Moving `geometric` into the defaults later is therefore safe;
starting with both and retreating is not.

#### Bridges, and the one place a feature conjunction is principled

`umol-geometric-graph` is meaningful only when both model classes are present. That is not a
packaging workaround but what a bridge *is* under [015](015-four-domains-semantic-model-2025-03-10.md):
a conversion between instances of two model classes. So the two-feature condition encodes the
architecture rather than working around Cargo.

Cargo cannot declaratively activate a dependency on a conjunction of features, so pick a spelling:

- **(a) `geometric` implies `graph`.** `geometric = ["graph", "dep:umol-geometric",
  "dep:umol-geometric-graph", ...]`. Simple, and honest while most of `umol-geometric`'s current value
  is reached through the bridge. Cost: the geometric model cannot stand alone.
- **(b) Separate the bridge.** `geometric` for the model class, `graph-geometric` for the conversion.
  Pure under four domains — model classes independent, conversions explicit — and what is wanted the
  day `umol-geometric` is useful by itself.

**Take (a) now, record (b) as the destination.** One line to switch, and today the distinction has no
practical consequence. Within the facade, bridge re-exports are gated with
`#[cfg(all(feature = "graph", feature = "geometric"))]` under either spelling.

**Not `default = []`, but for reasons that survive either framing.** The empty-default pattern is
legitimate — `hyper` 1.0 and tokio use it. Three arguments against it here hold regardless of how
\umol\ is read: every whitepaper listing assumes `cargo add umol` yields a working library, so an
empty default breaks reproducibility at line one; a library without an established reputation cannot
afford a first experience of "nothing is here"; and feature unification blunts the compile-time
benefit in any real dependency graph.

A fourth argument — "there is nothing to opt into" — is **contingent and should not be relied on.** It
was true only under the cheminformatics-first reading. Model classes are the orthogonal axes that
make opt-in coherent, and under the 015 reading there will be more of them. Default to the model
class that is the entry point; do not default to nothing.

### Re-exports

`umol-py` registers **189 classes**, which settles the shape: a flat root namespace is unusable and
hand-curating 189 names at the root is worse. Use modules mirroring the whitepaper's own taxonomy,
which readers arriving from the paper already carry:

- `umol::graph` — the graph model class: its IR (the entity ASTs, constraints, deltas), resolution,
  validation, matching, transformation, and its own `prelude`
- `umol::geometric` — the geometric model class, behind its feature, with its own `prelude`
- `umol::chem` — shared: elements and the chemistry model
- `umol::io` — placement open, see above
- `umol::prelude` — re-exports the default model class's prelude

Plus a small root-level set — the dozen types a trivial program needs — so that simple use does not
require knowing the module layout.

Re-exports should be **curated rather than globbed**. `pub use umol_ast::*` makes the facade a
mirror; a chosen surface makes it the definition. The latter is the artifact worth having.

### Prelude

**Include one, and the traits are the reason.** `Lattice` and `Canonicalize` are traits, so `meet`,
`join`, `matches` and `canonical_eq` are unavailable unless they are in scope. A user who cannot call
`meet` on an `AtomAst` concludes the method does not exist. That is a stronger argument than the
usual convenience case.

Roughly a dozen items: the traits from `umol-ast/src/ast/traits.rs` that a caller actually needs
(`Lattice`, `Canonicalize`, the conversion pair), the constantly used types (`MoleculeAst`,
`ReactionAst`, `AtomAst`, `BondAst`, `Element`), `ChemistryModel`, and the error types. A prelude
that imports 189 names is a glob with extra steps.

## The feature matrix must be built in CI

**The feature matrix is the untyped part of a Rust API.** Everything else in the facade is
type-checked, but `#[cfg(all(feature = "graph", feature = "geometric"))]` code is compiled only when
something actually builds that combination. A configuration nobody builds is a configuration nobody
type-checks, and it fails the way an untyped API fails: silently, in the combination no one tried.

Two features give four combinations, and each is a real claim this document makes:

| features | what it must prove |
| --- | --- |
| default (`graph`) | the ordinary path; every whitepaper listing compiles |
| `--no-default-features --features geometric` | **the peer claim.** If this does not build, model classes are not peers whatever the manifest says |
| `--no-default-features --features graph` | `graph` is genuinely optional rather than merely listed |
| `--all-features` | the bridge — the only place `cfg(all(...))` code is reachable |

The second row is the one that matters most and the one most likely to be skipped, because nobody
uses it day to day. It is the executable form of the argument in this document that peerhood is a
property of the feature graph rather than of the default set. If it is not in CI, that argument is
an assertion.

Add all four to the CI matrix set up under
[163](163-package-release-preparations-2026-07-26.md), and add `--no-default-features` alone if a
model-class-free build is meant to be coherent.

## The completeness test

`umol-py` depends on six member crates directly — `umol-ast`, `umol-chem`, `umol-graph`,
`umol-graph-core`, `umol-io`, `umol-perm` — and exposes 189 classes.

**If the facade is complete, `umol-py`'s manifest reduces to `umol` alone.** Anything the bindings
must reach past the facade to obtain is a gap in the facade. This turns surface design into gap
analysis against something that already ships, rather than a first-principles exercise, and the 189
names are the checklist.

Converting `umol-py` to depend on the facade is optional as a shipped change; using it as a test is
not.

## Out of scope

- Relocating ingest out of `umol-graph`. Worth doing for the layering; wants its own document and
  should not block this crate.
- Anything touching the nauty dependency. Settled in 137, 143 and 121; not a facade question.
- Changing any member crate's public API. The facade selects and re-exports; it does not redesign.
- The naming question in [176](176-ast-naming-2026-07-31.md). The facade is where that decision
  becomes most visible, since it is the curated public surface, but re-exports are one line each and
  do not force it. Build the facade against whatever names exist; apply 176 when it lands.

## Notes

- Version with the workspace at 0.6.0 per 163. `0.0.0 -> 0.6.0` is a valid increase; the placeholder
  version cannot be reused or deleted.
- 163's publish list has eighteen crates and does not include `umol`. It needs a nineteenth entry,
  published last, after everything it re-exports.
- The implementer should choose the staging. The natural shape is: manifest and feature wiring, then
  one module per subsystem, then the prelude, then the `umol-py` completeness check.
