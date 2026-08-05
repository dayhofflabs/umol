# 185 — Expose the reaction span form on the Python surface

Status: Proposed
Date: 2026-08-04
Relates: [179](179-python-editing-and-transactions-2026-08-02.md),
[182](182-python-resolution-2026-08-03.md),
[184](184-deltas-and-edits-2026-08-04.md),
[168](168-api-hygiene-2026-07-27.md)

`ReactionSpanAst` is not exported to `umol-py`. The whitepaper's reactions section names the span form
as the route by which an existing corpus of rules reaches \umol, so the route should be reachable
from the interface most users have.

## Justification

The field uses both mapped side-pair forms, such as SMIRKS and reaction SMARTS, and explicitly
superimposed forms, such as the condensed graph of reaction. `ReactionAst::from_sides` bridges the
former; direct `ReactionSpanAst` construction bridges the latter without first splitting the native
span into two complete molecules. \umol\ can express both forms, and they are interconvertible,
which is the practical answer to "must I re-encode my rules?"

The paper states this in prose and shows no listing, deliberately: demonstrating a second
representation would ask the reader to learn one, and the section already carries anchoring, deltas,
completeness, and composition. **So this work does not block the paper.** It matters because the paper
names a path, and a named path that cannot be walked from Python is the same gap docs 178 and 179
closed for lattice operations and editing.

## Scope

**In:**

- `ReactionSpanAst`, with `parse`, `parse_with_metadata`, `render`, `render_with_metadata`, and
  `str`. The span form is notation first; a caller who has one written down needs to read it in.
- Direct asserted and checked construction of `ReactionSpanAst` from its superimposed topology and
  per-entity before/after states in Rust, with checked construction at the Python boundary. A format
  adapter must not have to generate umol DSL text or split its native span into two complete
  molecules before it can use the span model.
- Rename the existing public molecule construction vocabulary from `MoleculeParts` /
  `MoleculeAst::from_parts` to `MoleculeEntries` / `MoleculeAst::from_entries` in Rust and from
  `MoleculeAst.from_parts` to `MoleculeAst.from_entries` in Python. Introduce the reaction-span
  construction API directly under the same `Entries` / `from_entries` convention. Add paired
  `try_from_entries` methods for checked runtime input; ordinary Rust construction remains
  infallible.
- The conversion in both directions: `ReactionAst.to_reaction_span()` and
  `ReactionSpanAst.to_reaction()`. These exist in Rust at `umol-ast/src/ast/reaction_span.rs:900` and
  `:546`. Both directions, since the claim being supported is that the representations are equivalent.
- A public checked `Correspondence` constructor in Rust and Python. The same type returned by
  substructure matching, molecule combination, and splitting should be accepted by subsequent
  operations; raw pair lists must not form a second construction channel.
- `ReactionAst.from_sides(lhs, rhs, atom_correspondence)` — construction from two molecules and an
  atom correspondence, which is how a rule arrives from an atom-mapped reaction.
- Rename the atom-level correspondence accessor on `ReactionDerivation` from `atom_map` to
  `atom_correspondence` in Rust and Python. Reserve atom mapping terminology for external-format atom
  labels.
- `lhs()`, `rhs()`, and `correspondence()` on the span. Together they project the two molecule sides
  and their existing `MoleculeCorrespondence` from the superimposed representation.

**Out:**

- Direct raw-pair input to `ReactionAst.from_sides`.
- `superimpose` and `difference_to` on the Python surface. Their Rust correspondence-precondition
  contracts belong to the broader public-consumer audit in doc 168.
- The `EntitySpan` accessors — `atoms()`, `bonds()`, `constraints()`. The paper describes the span
  form; it does not inspect one entry by entry. Python construction uses explicit before/after
  values and does not require another eight-class wrapper family.
- `ReactionSpanAst.reverse()`.
- Preservation of source-format spelling, annotations, agents, or other boundary metadata. A format
  adapter owns data that has no representation in the reaction ASTs.
- A representation of non-bijective source identities. One reaction span encodes a partial
  bijection; an adapter must reject a source identity relation that cannot be projected into that
  model or apply a separately specified projection policy.
- A general redesign of `GraphCorrespondence`, `MoleculeCorrespondence`, or every operation that
  consumes them. Doc 168 owns that review.

## Direct span construction

Parsing the umol span DSL is a textual construction path, not a sufficient adapter API. A native
span-shaped format already has a union topology and before/after values; forcing its adapter to emit
an intermediate DSL string would make the serialization layer part of ordinary in-memory
construction. Splitting the input into two molecules and calling `ReactionAst::from_sides` is
semantically valid, but bypasses `ReactionSpanAst` rather than making it the format bridge.

### Rust

The public carrier and constructor vocabulary describes semantic entries, not raw storage parts:

```rust
MoleculeEntries
MoleculeAst::from_entries(MoleculeEntries)
MoleculeAst::try_from_entries(MoleculeEntries)

ReactionSpanEntries
ReactionSpanAst::from_entries(ReactionSpanEntries)
ReactionSpanAst::try_from_entries(ReactionSpanEntries)
```

Rename the existing public `MoleculeParts` and `MoleculeAst::from_parts` API accordingly, then add a
public `ReactionSpanEntries` and the paired reaction-span constructors. The entries contain:

- `atoms: Vec<EntitySpan<AtomAst>>`;
- `bonds: Vec<(AtomId, AtomId, EntitySpan<BondAst>)>`;
- dative, aromatic, multicenter, noncovalent, stereo-atom, and stereo-bond entries in the same flat
  participant-plus-value shapes used by `MoleculeEntries`, with `EntitySpan<T>` as each value;
- `constraints: Vec<ConstraintSpan>`.

The flat form is the format-adapter boundary. Callers supply semantic atom, bond, and relation ids;
they do not assemble `Graph`, `RelationSet`, or `BiRelationSet` storage.

`try_from_entries` validates structural integrity before materializing the internal graph and
relation containers:

- the topology references existing union-frame atoms, bonds, and stereo ligands;
- a bond or overlay present on one side refers only to participants present on that side;
- a stereo entity present on one side refers to a site and ligand frame present on that side;
- molecule-level constraints refer to entities present on each side on which the constraint exists;
- an entry cannot be absent from both sides.

These are representation checks, not chemical validation. Neither construction path invokes
chemistry models, validators, or resolvers.

`from_entries` is the asserted path for entries whose structural integrity is established by their
producer. It uses the same checks and panics on a violated construction contract rather than
returning `Result`. `MoleculeAst::try_from_entries` returns `MoleculeEntriesError`, and
`ReactionSpanAst::try_from_entries` returns `ReactionSpanEntriesError`. Python maps both errors to
`ValueError` because the model value has not yet been constructed.

Apply the same split to the renamed molecule API. `MoleculeAst::from_entries` preserves the existing
infallible Rust construction shape. `MoleculeAst::try_from_entries` checks bond endpoints, overlay
participants, stereo sites and ligands, and constraint references before constructing graph and
relation storage. Python `MoleculeAst.from_entries` uses the checked Rust path and reports invalid
runtime input as `ValueError`.

The current crate-private reaction-span storage constructor is removed. The DSL parser already resolves
entries and performs several of the side-presence checks above; those checks move into or are shared
with `try_from_entries`. Runtime-data boundaries use the checked path; AST/DSL conversion and
internally generated spans use `from_entries` because their inputs establish the invariant by
construction. There is no second crate-private constructor with weaker semantics. Private
`from_parts` functions that genuinely assemble internal storage or pair ASTs with metadata are not
part of the public rename.

Both reaction-span construction paths normalize a `Modified` entity whose lhs and rhs values are
canonically equal to `Unchanged`. DSL parsing, Python before/after pairs, direct Rust construction,
and `superimpose` therefore produce the same span state.

### Python

Rename `MoleculeAst.from_parts` to `MoleculeAst.from_entries`, then expose
`ReactionSpanAst.from_entries` with the same call shape: `atoms` is
the required positional argument, while the remaining entity-family collections are keyword-only
and default to empty. Python does not need generic `EntitySpan<T>` wrappers. Each value is supplied
as an explicit `(lhs, rhs)` pair whose members are the corresponding AST type or `None`:

```python
span = ReactionSpanAst.from_entries(
    atoms=[(carbon, carbon), (None, oxygen)],
    bonds=[(0, 1, (None, single_bond))],
)
```

The pair maps to the Rust span state as follows:

- canonically equal values on both sides become `Unchanged`;
- different values on both sides become `Modified`;
- lhs only becomes `Removed`;
- rhs only becomes `Added`;
- `(None, None)` is invalid.

Bond and overlay entries carry their union-frame participants before the value pair. Stereo entries
reuse the existing Python `StereoLigand` value. Constraint pairs support unchanged, removed, and
added constraints; unequal values present on both sides are represented as separate removal and
addition entries because `ConstraintSpan` has no `Modified` state.

This is an input representation only. It does not add the entity-column inspection surface or expose
the graph-core storage containers. Both Python methods call their checked Rust
`try_from_entries` paths; Python does not expose a redundant `try_from_entries` spelling because
exceptions are the ordinary constructor-failure mechanism.

The textual surface also follows `MoleculeAst`. `parse` and `parse_with_metadata` accept an optional
keyword-only `MoleculeDefaults`; the latter returns `(ReactionSpanAst, MoleculeMetadata)`. `render`
and `render_with_metadata` use the same defaults and metadata types, and `str(span)` uses positional
rendering without metadata. The span has one union-frame namespace, so `MoleculeMetadata`, not
`ReactionMetadata`, is the corresponding persistent metadata type.

## Correspondence construction

Rust's `from_sides` currently takes a `Correspondence<NodeId>`, and `MoleculeCorrespondence` stores
the same graph-core id type for its atom family. These are molecular APIs: migrate that family to
`Correspondence<AtomId>` throughout `umol-ast`, converting to or from `NodeId` only where a
graph-core operation is actually called. `GraphCorrespondence` remains node- and edge-typed.

Python currently exposes `Correspondence` read-only and compensates by accepting raw pairs in
`from_sides`. This duplicates construction and validation at one consumer instead of making the
shared value constructible.

`Correspondence::new(matched_pairs, left_count, right_count)` becomes fallible. It establishes the
partial-bijection invariant by checking that:

- every left id is below `left_count` and every right id is below `right_count`;
- no left id occurs more than once;
- no right id occurs more than once;
- the stored pairs are ordered by left id.

A correspondence is intentionally partial. Unmatched ids express deletion and creation, so totality
is not a construction requirement.

Python exposes the same constructor as
`Correspondence(matched_pairs, left_count, right_count)`. Python ids remain ordinary non-negative
integers; `NodeId` and the entity-specific Rust id types do not cross the boundary. Invalid input is
reported as `ValueError` through the Rust `CorrespondenceError` rather than reimplemented in a
Python-specific helper.

`CorrespondenceError` also reports declared left- or right-space sizes that disagree with the
context supplied by an operation such as `ReactionAst::from_sides`. This keeps intrinsic pair
validation and contextual id-space validation under the same correspondence error without adding a
one-operation error type.

Its public variants are `LeftIdOutOfRange`, `RightIdOutOfRange`, `DuplicateLeftId`,
`DuplicateRightId`, `LeftCountMismatch`, and `RightCountMismatch`.

`Correspondence::from_images(images, right_count)` remains infallible. It is the dense-left form
used for algorithm-produced embeddings and remappings: left id `i` is paired with `images[i]`.
Its contract requires unique, in-range images and must be asserted rather than silently permitting an
invalid correspondence. This follows the established distinction between the asserted
`Permutation::from_image` path and fallible construction from untrusted runtime data.

The aggregating `GraphCorrespondence::new` and `MoleculeCorrespondence::new` constructors remain
infallible. They receive already-valid correspondences and have no graph or molecule against which
to check contextual dimensions.

## `from_sides`

Rust `ReactionAst::from_sides` takes `Correspondence<AtomId>` and becomes fallible. Python accepts
the constructible `Correspondence` object directly:

```python
atom_correspondence = Correspondence([(0, 0), (1, 1)], 2, 3)
reaction = ReactionAst.from_sides(lhs, rhs, atom_correspondence)
```

The operation checks that the correspondence's declared left and right spaces equal the atom counts
of `lhs` and `rhs`. This is contextual coherence, distinct from the partial-bijection invariant
established by `Correspondence::new`. Once the dimensions agree, the existing induced entity
correspondence and molecule-difference construction remain the operation.

No raw iterable of pairs is accepted by `from_sides`. Callers construct the reusable value once and
can inspect, reverse, or compose it through the existing `Correspondence` API.

## Settled semantics

- `to_reaction_span` returns `Contradiction` where Rust does; a reaction whose sides cannot be
  superimposed is an ordinary failure, not a defect.
- `from_sides` accepts a partial atom correspondence; unmatched lhs and rhs atoms become removals and
  additions respectively.
- Correspondence construction validates structural integrity only. It does not validate chemistry or
  require that matched atoms or bonds have equal attributes.
- Conversion does not mutate its receiver, consistent with 178.
- Round-tripping a reaction through the span form and back yields a canonically equal reaction. This
  is the property the paper's prose asserts, so it is the property the tests should state.

## Verification

Per 178 and 179: algebraic properties stay in Rust; the Python tests check availability and
representative cross-boundary results.

Rust verification covers:

- migration of the existing molecule construction surface to `MoleculeEntries` and paired
  `from_entries` / `try_from_entries` methods;
- exact molecule-entry reference failures through the checked path and the asserted contract of the
  infallible path;
- direct construction of every entity family and each span state;
- exact structural failures for missing union-frame participants and side-presence mismatches;
- normalization of canonically equal `Modified` entries to `Unchanged` through every construction
  path;
- equivalence between direct construction, DSL parsing, and superimposition for the same span;
- exact `Correspondence::new` failures for duplicate and out-of-range ids;
- ordering and unmatched-id behavior for valid partial correspondences;
- the asserted valid-image contract of `from_images` without changing its result type;
- `from_sides` dimension mismatches as errors rather than indexing failures;
- reaction/span round trips under partial correspondences.

Python verification covers direct span construction from before/after pairs, construction and
inspection of a partial `Correspondence`, error mapping to `ValueError`, typed `from_sides`,
parse/render of `ReactionSpanAst`, both conversion directions, and the `atom_correspondence`
accessor, and projection of the span's `lhs`, `rhs`, and `correspondence`. At least one round trip
creates an atom, since a partial correspondence is the case where the two forms differ most visibly.

## Staged implementation plan

### S0 — Establish checked and domain-typed correspondence construction

- **S0a — Make `Correspondence::new` establish the partial-bijection invariant.** In
  `umol-graph-core/src/correspondence.rs`, add the public `CorrespondenceError` with
  `LeftIdOutOfRange`, `RightIdOutOfRange`, `DuplicateLeftId`, `DuplicateRightId`,
  `LeftCountMismatch`, and `RightCountMismatch`; export it from the crate root. Change `new` to
  return `Result`, reject each invalid pair set, and retain sorting by left id for valid input.
  Migrate every workspace caller in the same subitem: propagate the error at runtime-data
  boundaries and explicitly assert producer invariants where an algorithm or already-validated AST
  constructs the correspondence. Do not add a second unchecked sparse-pair constructor. Add exact
  unit tables for all four constructor error variants, valid unsorted and partial inputs, and empty
  id spaces. The two contextual count-mismatch variants are exercised when introduced in S3a.
  The molecule property generator must test duplicate entity incidence directly rather than use an
  invalid induced correspondence as a uniqueness probe. `MoleculeCorrespondence::induce` retains
  the relation containers' incidence lookup and only drops later matches to an already-used right
  entity: duplicate-incidence input is validator-invalid, and the infallible operation guarantees a
  non-panicking partial correspondence rather than correctness for that input. Do not build a
  second incidence index or assign duplicate entities an arbitrary stable pairing.
  **Breaking (red→green).** [dep: none] **Done.**
- **S0b — Assert and verify dense-image construction.** Keep `from_images` infallible, but assert
  that every image is in the declared right space and occurs only once. Add exact asserted-contract
  cases and properties over generated valid partial bijections proving sorted storage, unique
  columns, range validity, and that matched plus unmatched ids partition each declared id space.
  Preserve the existing composition, reverse, and `compose_all` properties under the checked
  constructor. **Additive validation and tests (green).** [dep: S0a]
- **S0c — Move molecular atom correspondences from `NodeId` to `AtomId`.** Change the atom family
  stored by `MoleculeCorrespondence`, its `new`, `induce`, and `atoms` methods, and every molecular
  correspondence helper to `Correspondence<AtomId>`. Carry the type through molecule
  split/combine and remapping, substructure results, reaction construction and derivation,
  reaction-span superimposition and recovery, reaction composition, graph-layer reaction ingestion,
  and the Python Rust-conversion boundary. Convert `NodeId` only at calls into graph-core or the raw
  graph escape hatch; keep `GraphCorrespondence` unchanged. Migrate all fixtures, examples,
  benchmarks, unit tests, and property generators, and assert atom pairs with `AtomId` so the
  domain type is visible in the contract. Verify by source inventory that no public `umol-ast`
  signature exposes `Correspondence<NodeId>`. **Breaking (red→green).** [dep: S0a]

### S1 — Establish checked molecule-entry construction

- **S1a — Rename the Rust molecule construction vocabulary.** Rename `MoleculeParts` to
  `MoleculeEntries` and `MoleculeAst::from_parts` to `from_entries`, including crate-root exports,
  rustdoc links, macros, production callers, DSL conversions, fixtures, unit and property tests,
  fuzz targets, examples, and benchmarks throughout the workspace. Preserve the current asserted
  construction semantics in this subitem and verify representative empty, topology-only, overlay,
  stereo, and constrained molecules by full structural equality. Do not retain aliases for the old
  names. The Python-visible method may remain `from_parts` until S1c, but its Rust implementation
  uses the renamed API. **Breaking (red→green).** [dep: none]
- **S1b — Add checked molecule-entry construction.** In `ast/molecule.rs`, add
  `MoleculeEntriesError` and `MoleculeAst::try_from_entries`. Validate bond endpoints; every atom
  participant of dative, aromatic, multicenter, and noncovalent entries; stereo-atom sites and
  ligands; stereo-bond sites and ligands; and every entity reference in molecule constraints before
  constructing graph or relation storage. Keep these as representation-integrity checks only:
  self-loops, parallel entities, chemistry, and constraint satisfiability remain validator concerns.
  Make `from_entries` use the same validation and panic on a violated asserted-construction
  contract. Add exact error tables for every reference family and properties showing that valid
  entry sets produce the same molecule through both constructors. **Additive API with strengthened
  asserted contract (green).** [dep: S1a]
- **S1c — Move the Python molecule constructor to the checked entry path.** Rename
  `MoleculeAst.from_parts` to `MoleculeAst.from_entries`, change it to return `PyResult`, and map
  `MoleculeEntriesError` to `ValueError`. Migrate all Python callers and tests without retaining the
  old spelling. Test the full keyword-only entity-family surface and representative invalid atom,
  bond-site, ligand, and constraint references with exact exception messages. **Breaking
  (red→green).** [dep: S1b]

### S2 — Establish the public reaction-span entry model

- **S2a — Add `ReactionSpanEntries` and its paired constructors.** In
  `ast/reaction_span.rs`, add the public flat `ReactionSpanEntries`,
  `ReactionSpanEntriesError`, `ReactionSpanAst::from_entries`, and
  `ReactionSpanAst::try_from_entries`. Validate union-frame references, side-specific participant
  presence, stereo sites and ligand frames, constraint references on each present side, and the
  prohibition on an entity absent from both sides before constructing storage. Normalize every
  canonically equal `Modified` value to `Unchanged` in both paths. Test all four entity span states,
  all eight entity families, all three constraint states, every structural failure category, and
  exact normalization. **Additive (green).** [dep: S1a]
- **S2b — Route generated spans through entries.** Migrate `ReactionSpanAst::superimpose`,
  `ReactionAst::to_reaction_span`, and other constructors inside
  `ast/reaction_span.rs` to build `ReactionSpanEntries` and use the asserted constructor. Preserve
  union ordering, participant remapping, constraint ordering, and side compaction. Test exact
  generated spans containing created and removed atoms, bonds, every overlay family, stereo frames,
  and constraints; compare whole span values rather than family counts. **Internal rewire (green).**
  [dep: S2a]
- **S2c — Route the reaction-span DSL through entries.** Refactor
  `dsl/reaction_span.rs` so parsed input resolves into `ReactionSpanEntries` and calls the checked
  constructor, mapping structural failures into the existing parse-error surface. Move duplicated
  side-presence checks to the entry validator. Make `IntoAst` and `FromAst` use the asserted entry
  path after their type-directed conversions. Preserve `MoleculeMetadata`, defaults, keyword, and
  alias behavior. Add exact parsing errors plus direct-construction/DSL and
  DSL/superimposition equivalence cases, including canonically equal input sides normalizing to
  `Unchanged`. **Internal rewire (green).** [dep: S2a, S2b]
- **S2d — Remove the raw-storage span constructor.** Remove the crate-private
  `ReactionSpanAst::from_parts` and migrate any remaining production, fixture, property, or fuzz
  construction to `from_entries` or `try_from_entries` according to provenance. Verify by source
  inventory that no caller can bypass the flat entry contract, then run the complete `umol-ast`
  unit and property suites. **Breaking cleanup (red→green).** [dep: S2b, S2c]

### S3 — Complete the Rust reaction bridge

- **S3a — Make `ReactionAst::from_sides` check correspondence dimensions.** Change the method to
  accept `Correspondence<AtomId>` and return `Result<ReactionAst, CorrespondenceError>`, require the
  correspondence's declared left and right counts to equal the two molecule atom counts, and only
  then induce the remaining entity correspondences and derive the reaction. Migrate every Rust
  caller according to provenance. Add exact left/right count-mismatch tables and
  partial-correspondence cases containing both removals and additions. **Breaking (red→green).**
  [dep: S0c]
- **S3b — Align reaction-derivation terminology.** Rename
  `ReactionDerivation::atom_map` to `atom_correspondence` in Rust and migrate all callers, rustdoc,
  repr expectations, and tests. Do not retain the old accessor. **Breaking (red→green).**
  [dep: S0c]
- **S3c — State the bridge laws as Rust properties.** Extend the reaction property suite to verify
  that `from_sides` followed by span conversion preserves both projected sides under a partial atom
  correspondence; reaction → span → reaction is canonically equal; and direct entries, DSL parsing,
  and superimposition agree for generated structurally valid spans. Include unmatched atoms and
  dependent bonds or overlays rather than restricting the generator to total correspondences.
  **Additive tests (green).** [dep: S2d, S3a]

### S4 — Expose reusable correspondence construction in Python

- **S4a — Make `Correspondence` constructible.** In `umol-py/src/correspondence.rs`, add
  `Correspondence(matched_pairs, left_count, right_count)` backed exclusively by the checked Rust
  constructor, add the internal typed conversion needed by Rust consumers, and map every
  `CorrespondenceError` variant to `ValueError`. Retain the existing immutable accessors,
  composition, reversal, equality, and repr. Test empty, unsorted, partial, duplicate, and
  out-of-range construction with exact values and exceptions. **Additive (green).** [dep: S0a]
- **S4b — Require the correspondence value in Python `from_sides`.** Change
  `ReactionAst.from_sides(lhs, rhs, atom_correspondence)` to accept `Correspondence` directly,
  remove the raw-pair parser and its duplicate `HashSet` validation, and map the Rust method's count
  mismatches to `ValueError`. Test reuse of a constructed correspondence, unmatched ids, both count
  mismatches, and rejection of raw pair lists by the typed signature. **Breaking (red→green).**
  [dep: S3a, S4a]
- **S4c — Align the Python derivation accessor.** Rename
  `ReactionDerivation.atom_map` to `atom_correspondence`, migrate Python tests and repr expectations,
  and verify that it returns the same constructible `Correspondence` value used by `from_sides`.
  **Breaking (red→green).** [dep: S3b, S4a]

### S5 — Expose `ReactionSpanAst` in Python

- **S5a — Bind direct reaction-span construction.** Add and register the Python
  `ReactionSpanAst` value, with `from_entries` taking required atom pairs and keyword-only bond,
  overlay, stereo, and constraint entries. Convert each `(lhs, rhs)` pair to the corresponding span
  state; reject `(None, None)`; split unequal constraint pairs into removal and addition; and call
  `ReactionSpanAst::try_from_entries` so all structural failures become `ValueError`. Reuse the
  existing AST, `StereoLigand`, and `Constraint` wrappers rather than exposing `EntitySpan` or
  graph-core storage. Test every entity family, all four entity states, constraint splitting,
  canonical normalization, and representative union- and side-reference failures. **Additive
  (green).** [dep: S2d]
- **S5b — Bind the textual reaction-span surface.** Add `parse`, `parse_with_metadata`, `render`,
  `render_with_metadata`, and `__str__` using `ReactionSpanDsl`, keyword-only
  `MoleculeDefaults`, and `MoleculeMetadata`. Match the established `MoleculeAst` parse/render
  behavior and error mapping. Test positional rendering, metadata-preserving rendering, defaults,
  keyword and alias retention, and parse/render round trips. **Additive (green).** [dep: S5a]
- **S5c — Bind span projections and conversions.** Add `lhs()`, `rhs()`, and
  `correspondence()` to `ReactionSpanAst`; add `ReactionSpanAst.to_reaction()` and
  `ReactionAst.to_reaction_span()`, mapping `Contradiction` through the existing Python error
  surface. Test both conversion directions, projected side values, correspondence contents, and a
  round trip with an unmatched rhs atom so the bridge exercises creation rather than only relabeling.
  **Additive (green).** [dep: S3c, S4a, S5a]

### S6 — Verify and close the public bridge

- **S6a — Run cross-boundary verification and close the inventories.** Run formatting, complete
  workspace tests, the enabled `umol-graph-core` and `umol-ast` property suites, clippy, rebuild the
  Python extension under the Python 3.13 `umol-py` virtual environment, and run the full Python
  suite. Verify by source inventory that `MoleculeParts`, public `from_parts`, raw-pair Python
  `from_sides`, `atom_map`, and the raw-storage reaction-span constructor are absent, while every
  declared Python method is registered and importable. Update this document and `000-status.md`
  only after all checks pass. **Verification and documentation (green).**
  [dep: S1c, S2d, S3c, S4c, S5b, S5c]

## Plan properties

- **Critical path:** S0a → S0c → S3a → S3c → S5c → S6a, converging with the span-construction
  branch S1a → S1b → S2a → S2b → S2c → S2d and the Python branch S0a → S4a.
- **Independent work:** S1a has no dependency on checked correspondence construction; S4a can
  proceed once S0a is complete, and S5a can proceed once S2d is complete.
- **Green boundaries:** every stage ends with the affected crate and all workspace dependents
  compiling and passing their focused tests; S0a, S1a, S1c, S3a, S3b, S4b, and S4c contain their
  complete caller migrations because they change public names or signatures.
- **Deferrability:** none. S0–S3 establish the safe Rust bridge required by the Python API; S4–S5
  make that bridge usable from Python; S6 establishes the advertised public contract.
