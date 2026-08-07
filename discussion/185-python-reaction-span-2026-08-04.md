# 185 — Expose the reaction span form on the Python surface

Status: **In Progress**
Date: 2026-08-04
Relates: [179](179-python-editing-and-transactions-2026-08-02.md),
[182](182-python-resolution-2026-08-03.md),
[184](184-deltas-and-edits-2026-08-04.md),
[168](168-api-hygiene-2026-07-27.md),
[data type guide](../docs/development/data-types.md)

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
  `ReactionSpanAst.to_reaction()`. The former remains fallible because a permissive `ReactionAst`
  may contain deltas whose projected side is not a structurally valid molecule; the latter is
  infallible because a constructed span guarantees both sides.
- A public checked `Correspondence` constructor in Rust and Python. The same type returned by
  substructure matching, molecule combination, and splitting should be accepted by subsequent
  operations; raw pair lists must not form a second construction channel.
- The correspondence API corrections required by this bridge: intrinsic constructor errors,
  directional totality, fallible dense remapping, contextual graph and molecule induction, and
  contextual span superimposition and molecule differencing. Doc 168 retains the repository-wide
  audit beyond these identified consumers.
- `ReactionAst.from_sides(lhs, rhs, atom_correspondence)` — construction from two molecules and an
  atom correspondence, which is how a rule arrives from an atom-mapped reaction.
- Rename the atom-level correspondence accessor on `ReactionDerivation` from `atom_map` to
  `atom_correspondence` in Rust and Python. Reserve atom mapping terminology for external-format atom
  labels.
- Infallible `lhs()`, `rhs()`, and `correspondence()` on the span. A `ReactionSpanAst` represents an
  actual span of two referentially intact molecules, so checked construction rejects entries that
  cannot form either side before any projection is exposed.

**Out:**

- Direct raw-pair input to `ReactionAst.from_sides`.
- `superimpose` and `difference_to` on the Python surface. Their Rust correspondence-precondition
  contracts are corrected in this work, but no Python methods are added.
- The `EntitySpan` accessors — `atoms()`, `bonds()`, `constraints()`. The paper describes the span
  form; it does not inspect one entry by entry. Python construction uses explicit before/after
  values and does not require another eight-class wrapper family.
- `ReactionSpanAst.reverse()`.
- Preservation of source-format spelling, annotations, agents, or other boundary metadata. A format
  adapter owns data that has no representation in the reaction ASTs.
- A representation of non-bijective source identities. One reaction span encodes a partial
  bijection; an adapter must reject a source identity relation that cannot be projected into that
  model or apply a separately specified projection policy.
- Correspondence consumers beyond the graph induction, molecule induction, remapping,
  superimposition, and differencing paths identified here. Doc 168 owns that repository-wide
  review.
- Public error-surface changes to `ReactionAst::reverse`, reaction composition, or reaction
  fingerprints. Their current use of a materialized span is an implementation choice;
  strengthening span construction does not change their operation-specific contracts.

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
- every bond, overlay, stereo, and constraint reference resolves in the union-frame namespace;
- the entries selected on each side independently form a referentially intact `MoleculeAst`;
- an entry cannot be absent from both sides.

These are representation checks, not chemical or model-level validation. `ReactionSpanAst` means a
span of two molecules, so a bond, overlay, stereo entity, or constraint cannot be present on a side
unless every reference it requires is present on that side. An arbitrary union graph annotated with
left/right values but lacking either projected graph is not a reaction span. Neither construction
path invokes chemistry models, resolvers, canonicalization, repair, or implicit closure. The
general boundary is defined in the [data type guide](../docs/development/data-types.md).

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

The current crate-private reaction-span storage constructor is removed. The DSL parser already
resolves entries; union-frame and per-side reference checks move into or are shared with
`try_from_entries`. Runtime-data boundaries use the checked path; AST/DSL conversion and
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

The pair maps to the Rust span state as follows. Unlike `ReactionSpanEntries`, the Python pair does
not carry an explicit `EntitySpan` variant, so equality is the only available distinction between
unchanged and modified input:

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

Rust's `from_sides` takes a `Correspondence<AtomId>`, and `MoleculeCorrespondence` stores the same
domain id type for its atom family. Conversion to or from `NodeId` occurs only where a graph-core
operation is called. `GraphCorrespondence` remains node- and edge-typed.

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

`CorrespondenceError` is limited to failures of the carrier's own partial-bijection invariant. Its
public variants are `LeftIdOutOfRange`, `RightIdOutOfRange`, `DuplicateLeftId`, and
`DuplicateRightId`. A declared id-space size that disagrees with a separately supplied graph or
molecule is a contextual application failure and does not belong to the constructor error.

`Correspondence::from_images(images, right_count)` remains infallible. It is the dense-left form
used for algorithm-produced embeddings and remappings: left id `i` is paired with `images[i]`.
Its contract requires unique, in-range images and must be asserted rather than silently permitting an
invalid correspondence. This follows the established distinction between the asserted
`Permutation::from_image` path and fallible construction from untrusted runtime data.

The aggregating `GraphCorrespondence::new` and `MoleculeCorrespondence::new` constructors remain
infallible. They receive already-valid correspondences and have no graph or molecule against which
to check contextual dimensions.

`Correspondence` adds `is_total_on_left` and `is_total_on_right`; `is_total` continues to mean both.
`GraphCorrespondence` and `MoleculeCorrespondence` aggregate the same three predicates over their
entity families. Dense remapping conversion returns `Option`: `to_remapping` returns `None` unless
every family is total on the left. This is a normal absence condition because a valid
correspondence may be partial even when it was produced for the exact object pair.

Graph- and molecule-level induction accept an open correspondence alongside independent objects.
They are therefore provisionally `Option`-returning pending the repository-wide error review. They
return `None` when declared carrier sizes disagree with the supplied objects or when the remaining
entity correspondence is not uniquely inducible. They do not silently retain the first of several
possible induced matches. Internal callers whose producer establishes the contextual invariant may
assert the result; public consumers of independently supplied objects propagate absence. The
provenance rule and the containment of fallibility are specified in
[data type guide](../docs/development/data-types.md).

## `from_sides`

Rust `ReactionAst::from_sides` takes `Correspondence<AtomId>` and provisionally returns `Option`
pending the error review. Python accepts the constructible `Correspondence` object directly:

```python
atom_correspondence = Correspondence([(0, 0), (1, 1)], 2, 3)
reaction = ReactionAst.from_sides(lhs, rhs, atom_correspondence)
```

The operation checks that the correspondence's declared left and right spaces equal the atom counts
of `lhs` and `rhs`, and that the remaining entity families can be induced uniquely. This is
contextual coherence, distinct from the partial-bijection invariant established by
`Correspondence::new`. Python maps `None` to `ValueError`; a later error review may replace the
provisional absence with a named operation error without changing the successful call shape.

No raw iterable of pairs is accepted by `from_sides`. Callers construct the reusable value once and
can inspect, reverse, or compose it through the existing `Correspondence` API.

## Settled semantics

- `ReactionAst` remains a permissive lhs-plus-deltas carrier and may contain a DPO-invalid rule;
- `ReactionAst::to_reaction_span` retains its existing `Contradiction` result and rejects a rule
  whose deltas cannot form two referentially intact projected molecules;
- `ReactionSpanAst::try_from_entries` validates union-frame references and the structural integrity
  of both projected sides;
- `ReactionSpanAst::lhs`, `rhs`, and `to_reaction` are infallible because their required structural
  invariants are established at span construction;
- `DpoValidator::validate_reaction_span` is removed as redundant. Its current rule-level dangling
  check is implied by the stronger span construction invariant; `DpoValidator::validate_reaction`
  and match-dependent application checks remain;
- `from_sides` accepts a partial atom correspondence; unmatched lhs and rhs atoms become removals and
  additions respectively.
- Correspondence construction validates structural integrity only. It does not validate chemistry or
  require that matched atoms or bonds have equal attributes.
- Conversion does not mutate its receiver, consistent with 178.
- Round-tripping a reaction through the span form and back may replace relative deltas with absolute
  updates. Materializing the recovered reaction must reproduce the same span; equality of delta
  syntax is not the semantic contract.

## Verification

Per 178 and 179: algebraic properties stay in Rust; the Python tests check availability and
representative cross-boundary results.

Rust verification covers:

- migration of the existing molecule construction surface to `MoleculeEntries` and paired
  `from_entries` / `try_from_entries` methods;
- exact molecule-entry reference failures through the checked path and the asserted contract of the
  infallible path;
- direct construction of every entity family and each span state;
- exact structural failures for missing union-frame references;
- rejection of entries that are union-valid but cannot form the lhs or rhs, including exact
  failures for bond, overlay, stereo, and constraint references absent from the selected side;
- infallible, faithful projection of both normalized sides for every constructed span, with exact
  lhs preservation and rhs equivalence under the induced total reaction-frame correspondence;
- removal of the redundant reaction-span DPO-validator entry point while retaining reaction and
  match-dependent DPO checks;
- normalization of canonically equal `Modified` entries to `Unchanged` through every construction
  path;
- equivalence between direct construction, DSL parsing, and superimposition for the same span;
- exact `Correspondence::new` failures for duplicate and out-of-range ids;
- ordering and unmatched-id behavior for valid partial correspondences;
- the asserted valid-image contract of `from_images` without changing its result type;
- directional totality laws and `None` from dense remapping of a correspondence that is not total on
  the left;
- contextual graph and molecule induction failures for carrier-size mismatch and non-unique
  incidence, while producer-established paths remain asserted;
- contextual superimposition and differencing failures for incompatible full correspondences,
  including acceptance of a coherent correspondence narrower than maximal induction;
- `from_sides` dimension mismatches as contextual failure rather than indexing failures;
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
  `LeftIdOutOfRange`, `RightIdOutOfRange`, `DuplicateLeftId`, and `DuplicateRightId`; export it from
  the crate root. Change `new` to
  return `Result`, reject each invalid pair set, and retain sorting by left id for valid input.
  Migrate every workspace caller in the same subitem: propagate the error at runtime-data
  boundaries and explicitly assert producer invariants where an algorithm or already-validated AST
  constructs the correspondence. Do not add a second unchecked sparse-pair constructor. Add exact
  unit tables for all four constructor error variants, valid unsorted and partial inputs, and empty
  id spaces. The molecule property generator must test duplicate entity incidence directly rather
  than use an invalid induced correspondence as a uniqueness probe.
  **Breaking (red→green).** [dep: none] **Done.**
- **S0b — Assert and verify dense-image construction.** Keep `from_images` infallible, but assert
  that every image is in the declared right space and occurs only once. Add exact asserted-contract
  cases and properties over generated valid partial bijections proving sorted storage, unique
  columns, range validity, and that matched plus unmatched ids partition each declared id space.
  Preserve the existing composition, reverse, and `compose_all` properties under the checked
  constructor. **Additive validation and tests (green).** [dep: S0a] **Done.**
- **S0c — Move molecular atom correspondences from `NodeId` to `AtomId`.** Change the atom family
  stored by `MoleculeCorrespondence`, its `new`, `induce`, and `atoms` methods, and every molecular
  correspondence helper to `Correspondence<AtomId>`. Carry the type through molecule
  split/combine and remapping, substructure results, reaction construction and derivation,
  reaction-span superimposition and recovery, reaction composition, graph-layer reaction ingestion,
  and the Python Rust-conversion boundary. Convert `NodeId` only at calls into graph-core or the raw
  graph escape hatch; keep `GraphCorrespondence` unchanged. Migrate all fixtures, examples,
  benchmarks, unit tests, and property generators, and assert atom pairs with `AtomId` so the
  domain type is visible in the contract. Verify by source inventory that no public `umol-ast`
  signature exposes `Correspondence<NodeId>`. **Breaking (red→green).** [dep: S0a] **Done.**
- **S0d — Separate intrinsic construction errors from contextual failures.** Remove
  `LeftCountMismatch` and `RightCountMismatch` from `CorrespondenceError`, its formatting, and the
  Python conversion match. Keep the four pair/range variants exhaustive and limit the type's
  rustdoc to failures of the carrier's own partial-bijection invariant. Update this document's
  inventories and focused constructor tables so no contextual object-size check is attributed to
  `Correspondence::new`. **Breaking cleanup (red→green).** [dep: S0a] **Done.**
- **S0e — Make graph-correspondence context and directional totality explicit.** In
  `umol-graph-core/src/correspondence.rs`, add `Correspondence::is_total_on_left` and
  `is_total_on_right`, retain `is_total` as their conjunction, and aggregate the predicates on
  `GraphCorrespondence`. Change `GraphCorrespondence::to_remapping` to return `Option<Remapping>`.
  Make `Correspondence<NodeId>::edge_matched_pairs`, `shared_edge_count`, and
  `GraphCorrespondence::induced` return `None` when the node carrier does not describe the supplied
  graph pair or when parallel-edge incidence prevents a unique induced edge correspondence. Migrate
  graph rewriting and matching callers according to provenance: algorithm-produced pairs assert
  their contract, while public paths propagate absence. Add exact tables for left/right dimension
  mismatch, partial and total-left remapping, and ambiguous parallel edges. Add properties comparing
  induced edge matching against an exhaustive reference relation over generated multigraphs and
  relating aggregate graph totality and remapping to the component correspondences. **Breaking
  (red→green).** [dep: S0d] **Done.**
- **S0f — Make molecule-correspondence context and remapping fallible.** In
  `umol-ast/src/ast/correspondence.rs`, aggregate `is_total_on_left` and `is_total_on_right` over all
  eight families, retain `is_total` as totality on both sides, and change `to_remapping` to return
  `Option<IdRemapping>`. Change `MoleculeCorrespondence::induce` to return `Option<Self>`: require
  the atom carrier's declared sizes to equal the two molecule atom counts and return `None` when
  bond, overlay, or stereo incidence does not induce a unique right partner. Remove
  `retain_unique_rights` and its first-match semantics. Migrate all workspace callers by provenance;
  this includes changing `ReactionAst::from_sides` to return `Option<ReactionAst>`, while callers
  using correspondences produced for the same unchanged object pair may assert the invariant. Add
  exact tables for dimension mismatch, partial and total-left remapping, and non-unique incidence in
  the affected entity families, plus generated properties for directional totality and successful
  induction over structurally valid molecule pairs. The `Option` surface is provisional pending the
  repository-wide error review. **Breaking (red→green).** [dep: S0c, S0e] **Done.**

### S1 — Establish checked molecule-entry construction

- **S1a — Rename the Rust molecule construction vocabulary.** Rename `MoleculeParts` to
  `MoleculeEntries` and `MoleculeAst::from_parts` to `from_entries`, including crate-root exports,
  rustdoc links, macros, production callers, DSL conversions, fixtures, unit and property tests,
  fuzz targets, examples, and benchmarks throughout the workspace. Preserve the current asserted
  construction semantics in this subitem and verify representative empty, topology-only, overlay,
  stereo, and constrained molecules by full structural equality. Do not retain aliases for the old
  names. The Python-visible method may remain `from_parts` until S1c, but its Rust implementation
  uses the renamed API. **Breaking (red→green).** [dep: none] **Done.**
- **S1b — Add checked molecule-entry construction.** In `ast/molecule.rs`, add
  `MoleculeEntriesError` and `MoleculeAst::try_from_entries`. Validate bond endpoints; every atom
  participant of dative, aromatic, multicenter, and noncovalent entries; stereo-atom sites and
  ligands; stereo-bond sites and ligands; and every entity reference in molecule constraints before
  constructing graph or relation storage. Keep these as representation-integrity checks only:
  self-loops, parallel entities, chemistry, and constraint satisfiability remain validator concerns.
  Make `from_entries` use the same validation and panic on a violated asserted-construction
  contract. Add exact error tables for every reference family and properties showing that valid
  entry sets produce the same molecule through both constructors. The constructor-routing audit
  keeps DSL conversion, molecule combination, pushout, reaction-span projection, and perceived
  molecule conversion on the asserted path because they establish references by construction;
  TableIR raising uses the checked path, and the Python entry constructor moves to it in S1c.
  **Additive API with strengthened asserted contract (green).** [dep: S1a] **Done.**
- **S1c — Move the Python molecule constructor to the checked entry path.** Rename
  `MoleculeAst.from_parts` to `MoleculeAst.from_entries`, change it to return `PyResult`, and map
  `MoleculeEntriesError` to `ValueError`. Migrate all Python callers and tests without retaining the
  old spelling. Test the full keyword-only entity-family surface and representative invalid atom,
  bond-site, ligand, and constraint references with exact exception messages. **Breaking
  (red→green).** [dep: S1b] **Done.**

### S2 — Establish the public reaction-span entry model

- **S2a — Add `ReactionSpanEntries` and its paired constructors.** In
  `ast/reaction_span.rs`, add the public flat `ReactionSpanEntries`,
  `ReactionSpanEntriesError`, `ReactionSpanAst::from_entries`, and
  `ReactionSpanAst::try_from_entries`. Validate union-frame participant, site, ligand, and constraint
  references and the prohibition on an entity absent from both sides before constructing storage;
  do not enforce side presence or DPO semantics. Normalize every canonically equal `Modified` value
  to `Unchanged` in both paths. Test all four entity span states,
  all eight entity families, all three constraint states, every structural failure category, and
  exact normalization. Before adding another exhaustive constraint-reference walk, review the
  overlap with molecule-entry and reaction-integrity validation and use a common internal traversal
  where that reduces duplication without expanding the public API. Also avoid retaining DSL-side
  checks that merely repeat checks owned by the checked entry constructor. **Additive (green).**
  [dep: S1a] **Done. The initially implemented side checks were removed in S2b; the revised
  actual-span contract restores side structural checks in S3b.**
- **S2b — Correct the construction boundary and route generated spans through entries.** Remove
  side-presence validation from `ReactionSpanAst::try_from_entries`, retaining only union-reference
  integrity, and add exact construction cases for DPO-invalid but representable spans. Preserve the
  explicit `DpoValidator` path rather than duplicating it in construction. Migrate
  `ReactionSpanAst::superimpose`,
  `ReactionAst::to_reaction_span`, and other constructors inside
  `ast/reaction_span.rs` to build `ReactionSpanEntries`. Use the asserted constructor for
  `superimpose`, whose inputs already establish the entry invariants; use the checked constructor
  for `to_reaction_span` and retain its existing `Contradiction` contract. Preserve union ordering,
  participant remapping, and constraint ordering. Test exact
  generated spans containing created and removed atoms, bonds, every overlay family, stereo frames,
  and constraints; compare whole span values rather than family counts. **Internal rewire (green).**
  [dep: S2a] **Done. Its union-only construction decision is superseded by S3b; the entry routing
  and generated-span migration remain in force.**
- **S2c — Route the reaction-span DSL through entries.** Refactor
  `dsl/reaction_span.rs` so parsed input resolves into `ReactionSpanEntries` and calls the checked
  constructor, mapping structural failures into the existing parse-error surface. Remove semantic
  side-presence checks rather than moving them into construction. Make `IntoAst` and `FromAst` use
  the asserted entry path after their type-directed conversions. Preserve `MoleculeMetadata`,
  defaults, keyword, and alias behavior. Add exact parsing errors plus direct-construction/DSL and
  DSL/superimposition equivalence cases, including canonically equal input sides normalizing to
  `Unchanged`. **Internal rewire (green).** [dep: S2a, S2b] **Done. Its temporary removal of
  side checks is superseded by S3b; routing through the checked constructor remains in force.**
- **S2d — Remove the raw-storage span constructor.** Remove the crate-private
  `ReactionSpanAst::from_parts` and migrate any remaining production, fixture, property, or fuzz
  construction to `from_entries` or `try_from_entries` according to provenance. Verify by source
  inventory that no caller can bypass the flat entry contract, then run the complete `umol-ast`
  unit and property suites. **Breaking cleanup (red→green).** [dep: S2b, S2c] **Done.**
- **S2e — Make span superimposition contextually fallible.** Change
  `ReactionSpanAst::superimpose` and `MoleculeAst::difference_to` to return `Option`, and migrate
  every caller in the same subitem. Validate that each family declares the supplied molecule
  counts and that every supplied matched bond, overlay, or stereo pair has compatible incidence
  under the atom family. Permit a coherent correspondence to be narrower than the maximally induced
  correspondence; unmatched entities continue to encode removal and addition. Return `None` for an
  independently supplied incompatible full correspondence, while trusted producer paths assert
  their established context. Add exact cases for count mismatch, incompatible bond, overlay, and
  stereo incidence, a coherent narrower correspondence, and successful whole-span equality. The
  `Option` surface is provisional pending the error review. **Breaking (red→green).**
  [dep: S0f, S2d] **Done.**

### S3 — Complete the Rust reaction bridge

- **S3a — Align reaction-derivation terminology.** Rename
  `ReactionDerivation::atom_map` to `atom_correspondence` in Rust and migrate all callers, rustdoc,
  repr expectations, and tests. Do not retain the old accessor. **Breaking (red→green).**
  [dep: S0c] **Done.**
- **S3b — Establish the actual-span invariant and keep projection infallible.** Strengthen
  `ReactionSpanAst::try_from_entries` so that, after checking the union namespace, it verifies that
  the entries selected on each side satisfy the same referential-integrity contract as
  `MoleculeAst::try_from_entries`. Reuse the molecule-entry reference traversal rather than adding
  a divergent validator. Keep `ReactionSpanAst::lhs`, `rhs`, and `to_reaction` infallible; make the
  shared projection retain and remap every selected entity and constraint, relying on the
  constructor invariant instead of silently dropping a selected entry with an absent participant.
  Retain `ReactionAst::to_reaction_span`'s existing `Contradiction` surface and reject reactions
  whose deltas cannot form an actual two-sided span. Remove
  `DpoValidator::validate_reaction_span` and migrate its structural cases to checked span
  construction; retain `validate_reaction` and the match-dependent application checks. Add exact
  `#[rstest]` cases for union-valid entries rejected because the lhs or rhs lacks a required bond,
  overlay, stereo, or constraint reference, and exact whole-value tests for both projections of
  every accepted span. **Breaking (red→green).** [dep: S2d] **Done.**
- **S3c — State the bridge laws as Rust properties.** Extend the reaction property suite to verify
  that `from_sides` followed by span conversion preserves the lhs exactly and preserves the rhs up
  to the total correspondence induced by lhs-anchored reaction-frame reindexing; exercise crossing
  partial atom correspondences explicitly rather than assuming order preservation. Verify that
  every reaction that converts to a span round-trips to the same materialized span, and that direct
  entries, DSL parsing, and superimposition agree for generated structurally valid, lhs-anchored
  spans.
  Include unmatched atoms and dependent bonds or overlays rather than restricting the generator to
  total correspondences.
  **Additive properties (green).** [dep: S0f, S2e, S3b] **Done.**
- **S3d — Add checked point lookup to `Remapping`.** In `umol-graph-core`, keep
  `Remapping::new` infallible: its image vectors define a total function over their own source
  ranges, and injectivity or bijectivity is a precondition only for operations that require it. Add
  `try_map_node` and `try_map_edge`, returning `None` for an id outside the corresponding source
  range. Keep `map_node` and `map_edge` as the asserted point-lookup routes and implement them over
  the same stored images. Document their panic condition and test exact covered and uncovered node
  and edge lookups. **Additive checked surface (green).** [dep: S3c] **Done.**
- **S3e — Make relation-set remapping own positional data transport.** In `umol-graph-core`, change
  all five relation-set and birelation-set `apply_remapping` operations to return the remapped set
  without exposing participant permutations. Map each factor's participants and route the entries
  through the set's existing constructor, so factor canonicalization and `RelationData` or
  `BiRelationData::on_permutation` happen exactly once at the abstraction that owns them; do not
  pre-permute payloads or change `RelationParticipant::remap`.

  Add `try_apply_remapping` to all five set types. Determine coverage from each participant's
  `ParticipantRefs` using `try_map_node` and `try_map_edge`; return `None` if either referenced id is
  outside the remapping's source ranges, otherwise delegate to the same transport implementation as
  `apply_remapping`. The asserted route is for a remapping known to cover every participant and
  documents that a mismatch panics; the checked route handles independently supplied values.

  State and test the concrete positional-data law: after relabeling and factor canonicalization,
  each payload position remains attached to the same logical participant. Use tagged
  position-sensitive payloads across all five storage shapes, covering ordered and unordered
  factors and both factors of birelation data. Add exact checked-application cases for uncovered
  node and edge participants. No production caller currently uses these methods, so migrate their
  graph-core tests with the signature change and keep the workspace green.
  **Breaking API correction with additive checked path (red→green).** [dep: S3d] **Done.**
- **S3f — Route reaction-span superimposition through relation remapping.** Derive the graph-core
  participant remapping and the typed `IdRemapping` from the same lhs-anchored correspondence in
  `ReactionSpanAst::superimpose`. Use the corrected relation-set operation for every overlay and
  remove local participant sorting and payload permutation. Retain the exact lhs law and the rhs
  `MoleculeAst::equiv_under` property, including crossing correspondences and position-sensitive
  payloads. **Internal rewire (green).** [dep: S3e] **Done.**
- **S3g — Route molecule pushout through relation remapping.** Replace the corresponding manual
  participant canonicalization and payload permutation in molecule pushout with the corrected
  relation-set operation. Preserve the existing pushout result and correspondence contracts and
  add exact position-sensitive coverage where the shared graph-core properties do not exercise the
  molecule-level assembly. **Internal rewire (green).** [dep: S3e] **Done.**

### S4 — Expose reusable correspondence construction in Python

- **S4a — Make `Correspondence` constructible.** In `umol-py/src/correspondence.rs`, add
  `Correspondence(matched_pairs, left_count, right_count)` backed exclusively by the checked Rust
  constructor, add the internal typed conversion needed by Rust consumers, and map every
  `CorrespondenceError` variant to `ValueError`. Retain the existing immutable accessors,
  composition, reversal, equality, and repr. Test empty, unsorted, partial, duplicate, and
  out-of-range construction with exact values and exceptions. **Additive (green).** [dep: S0a]
  **Done.**
- **S4b — Require the correspondence value in Python `from_sides`.** Change
  `ReactionAst.from_sides(lhs, rhs, atom_correspondence)` to accept `Correspondence` directly,
  remove the raw-pair parser and its duplicate `HashSet` validation, and map the Rust method's
  provisional `None` result to `ValueError`. Test reuse of a constructed correspondence, unmatched
  ids, and contextual incompatibility.
  **Breaking (red→green).** [dep: S0f, S4a] **Done.**
- **S4c — Align the Python derivation accessor.** Rename
  `ReactionDerivation.atom_map` to `atom_correspondence`, migrate Python tests and repr expectations,
  and verify that it returns the same constructible `Correspondence` value used by `from_sides`.
  **Breaking (red→green).** [dep: S3a, S4a] **Done.**

### S5 — Expose `ReactionSpanAst` in Python

- **S5a — Bind direct reaction-span construction.** Add and register the Python
  `ReactionSpanAst` value, with `from_entries` taking required atom pairs and keyword-only bond,
  overlay, stereo, and constraint entries. Convert each `(lhs, rhs)` pair to the corresponding span
  state; reject `(None, None)`; split unequal constraint pairs into removal and addition; and call
  `ReactionSpanAst::try_from_entries` so all structural failures become `ValueError`. Reuse the
  existing AST, `StereoLigand`, and `Constraint` wrappers rather than exposing `EntitySpan` or
  graph-core storage. Test every entity family, all four entity states, constraint splitting,
  canonical normalization, representative union-reference failures, and rejection of entries
  whose lhs or rhs is not referentially intact. **Additive
  (green).** [dep: S3b]
- **S5b — Bind the textual reaction-span surface.** Add `parse`, `parse_with_metadata`, `render`,
  `render_with_metadata`, and `__str__` using `ReactionSpanDsl`, keyword-only
  `MoleculeDefaults`, and `MoleculeMetadata`. Match the established `MoleculeAst` parse/render
  behavior and error mapping. Test positional rendering, metadata-preserving rendering, defaults,
  keyword and alias retention, and parse/render round trips. **Additive (green).** [dep: S5a]
- **S5c — Bind span projections and conversions.** Add infallible `lhs()` and `rhs()`, plus
  `correspondence()` to `ReactionSpanAst`; add `ReactionSpanAst.to_reaction()` and
  `ReactionAst.to_reaction_span()`. Map the latter's contradictory or structurally unrepresentable
  deltas through `ContradictionError`; direct checked construction already reports invalid span
  entries as `ValueError`. Test both conversion directions, projected side values, correspondence
  contents, and a round trip with an unmatched rhs atom so the bridge exercises creation rather
  than only relabeling.
  **Additive (green).** [dep: S3b, S3c, S4a, S5a]

### S6 — Verify and close the public bridge

- **S6a — Run cross-boundary verification and close the inventories.** Run formatting, complete
  workspace tests, the enabled `umol-graph-core` and `umol-ast` property suites, clippy, rebuild the
  Python extension under the Python 3.13 `umol-py` virtual environment, and run the full Python
  suite. Verify by source inventory that `MoleculeParts`, public `from_parts`, raw-pair Python
  `from_sides`, `atom_map`, and the raw-storage reaction-span constructor are absent, while every
  declared Python method is registered and importable. Update this document and `000-status.md`
  only after all checks pass. **Verification and documentation (green).**
  [dep: S0f, S1c, S2e, S3b, S4c, S5b, S5c]

## Plan properties

- **Critical path:** S0a → S0c → S0d → S0e → S0f → S2e → S3b → S5c → S6a, converging with the
  span-construction branch
  S1a → S1b → S2a → S2b → S2c → S2d, and the Python branch S0a → S4a.
- **Independent work:** S1a has no dependency on checked correspondence construction; S4a can
  proceed once S0a is complete, and S5a can proceed once S2d is complete.
- **Green boundaries:** every stage ends with the affected crate and all workspace dependents
  compiling and passing their focused tests; S0a, S0d, S0e, S0f, S1a, S1c, S2e, S3a, S3b, S4b,
  and S4c contain their complete caller migrations because they change public names or signatures.
- **Deferrability:** none. S0–S3 establish the safe Rust bridge required by the Python API; S4–S5
  make that bridge usable from Python; S6 establishes the advertised public contract.
