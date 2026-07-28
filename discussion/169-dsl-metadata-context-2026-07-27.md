# DSL metadata and context

Status: **Proposed**
Date: 2026-07-27
Relates: [133](133-reaction-edn-surface-2026-06-26.md),
[134](134-reaction-application-overlays-2026-06-26.md),
[151](151-python-molecule-workflows-2026-07-13.md),
[164](164-dsl-edn-worklist-2026-07-27.md)

## Scope

This document defines the relationship between persistent DSL metadata and the
temporary contexts used to resolve DSL references. The Rust API must be
reconciled before metadata-preserving molecule and reaction parsing and
rendering are exposed in Python.

The design applies uniformly to the eight entity kinds:

1. atoms;
2. bonds;
3. dative bonds;
4. aromatic systems;
5. multicenter bonds;
6. noncovalent bonds;
7. stereo atoms;
8. stereo bonds.

## Concepts and lifetimes

`MoleculeMetadata` and `ReactionMetadata` are persistent surface objects.
They currently retain the keyword bindings and atom aliases required to render
the DSL after it has been resolved into numerical AST identifiers. The
`Metadata` name is deliberately broader than `Namespace`: future DSL
extensions may require fragments, ports, or metadata for other graph
primitives.

The current `MoleculeNamespace` and `ReactionNamespace` are temporary
resolution objects. They additionally hold entity counts and participant
indexes used to resolve positional, keyword, and structural references while
parsing. Rename them:

- `MoleculeNamespace` to `MoleculeContext`;
- `ReactionNamespace` to `ReactionContext`.

The generic `Namespace` trait may retain its name. It describes the
reference-resolution query surface implemented by the two concrete contexts.
Its overlap with the `Ctx` associated type used by `IntoAst` and `FromAst` is
minor and does not require changing those traits.

The two layers remain distinct:

| Layer | Lifetime | Persistent contents | Additional temporary contents |
| --- | --- | --- | --- |
| `MoleculeMetadata` | parse/render roundtrip | entity keywords and atom aliases | none |
| `MoleculeContext` | reference resolution | the same keyword bindings and aliases | counts and participant indexes |
| `ReactionMetadata` | parse/render roundtrip | lhs metadata, delta keywords, reaction aliases | none |
| `ReactionContext` | reference resolution | the same scoped bindings and aliases | lhs/delta counts and participant indexes |

Counts and participant indexes do not belong in metadata. Conversely, keyword
bindings should not be stored twice in opposite directions.

## Shared keyword representation

Reuse [`Entity`](../umol-ast/src/ast/entity.rs), the existing typed reference
to any of the eight molecule entity kinds, directly as the target of the
bidirectional entity-keyword maps:

```rust
BiBTreeMap<Entity, String>
```

Add derived `PartialOrd` and `Ord` implementations to `Entity`, complementing
its existing `PartialEq`, `Eq`, and `Hash` implementations. Ordering by entity
variant and then numerical id is well-defined and generally useful. The single
table makes both entity uniqueness and cross-kind keyword uniqueness
structural: an `Entity` has at most one keyword and a keyword names at most one
`Entity`. Typed operations such as `atom_id(keyword)` wrap and match the
corresponding `Entity` variant.

Atom aliases remain separate. They name atom DSL templates rather than AST
entities and retain their `BiBTreeMap<String, Box<AtomDsl>>` representation.
The containing metadata type checks keyword disjointness between the entity
table and the alias table.

```rust
pub struct MoleculeMetadata {
    keywords: BiBTreeMap<Entity, String>,
    atom_aliases: BiBTreeMap<String, Box<AtomDsl>>,
}

pub struct ReactionMetadata {
    lhs: MoleculeMetadata,
    delta_keywords: BiBTreeMap<Entity, String>,
    atom_aliases: BiBTreeMap<String, Box<AtomDsl>>,
}
```

A separate `EntityKeywords` newtype would add no representational invariant
beyond `BiBTreeMap`. Its collision policy would still be incomplete because
keyword/alias disjointness and reaction parent-scope checks require the
containing metadata object. Keep those operations on `MoleculeMetadata` and
`ReactionMetadata`, whose private fields already prevent unchecked external
mutation.

`keywords` is preferable to `entity_keywords` on `MoleculeMetadata`: the
receiver already identifies the metadata as molecule-shaped.
`delta_keywords` distinguishes reaction keywords defined in deltas from those
defined in the lhs.

`MoleculeContext` owns or builds the same `MoleculeMetadata` bindings alongside
its counts and participant indexes. Parsing moves the persistent portion out
of the context instead of reconstructing it by inverting eight maps.
`ReactionContext` does the same for the lhs and delta scopes.

## Lookup and mutation semantics

Metadata supports both directions:

```rust
metadata.atom_keyword(AtomId(3)) // -> Option<&str>
metadata.atom_id("carbonyl")     // -> Option<AtomId>
```

Here `AtomId` is the numerical AST identifier, while `"carbonyl"` is a DSL
keyword. The same pair of operations is required for every entity kind.

Insertion must preserve the DSL namespace invariants. It is fallible when it
would violate:

- bijectivity within one entity kind;
- keyword uniqueness across all eight entity kinds;
- disjointness between entity keywords and atom aliases;
- bijectivity of atom aliases.

The current metadata setters and alias insertion silently replace colliding
entries, whereas parsing rejects the same collisions. The reconciled API must
not retain this discrepancy. Context registration validates and inserts a
keyword atomically into the shared bindings.

`ReactionMetadata` provides layered lookup over delta and lhs bindings and
implements `Metadata` directly. Consumers that render constraints or overlay
references can accept that view without constructing a combined
`MoleculeMetadata`. The current `combined_metadata()` operation and its
render-time cloning should disappear.

## Metadata remapping

`MoleculeMetadata` follows a molecule through a
`MoleculeCorrespondence`. Add entity-level queries to the correspondence:

```rust
correspondence.right_of(Entity::Atom(AtomId(3)))
    // -> Option<Entity>
```

and the symmetric `left_of`. These operations dispatch to the existing typed
correspondence for the `Entity` variant.

`MoleculeMetadata::remap` maps the left entity namespace to the right:

- a keyword attached to a matched entity is attached to its right-hand entity;
- a keyword attached to an unmatched left entity is omitted;
- an unmatched right entity acquires no keyword;
- atom aliases are carried through unchanged because they name reusable atom
  DSL templates and contain no entity ids.

The operation is infallible. A `MoleculeCorrespondence` is a partial
bijection, so remapping a valid entity-keyword bimap cannot create a target
collision. It follows the existing consuming `remap(self, ...) -> Self`
convention in `umol-ast`; callers that need both values clone explicitly.

For a total correspondence, remapping through the correspondence and then its
reverse recovers the original metadata. Remapping respects correspondence
composition. For a partial correspondence, it is the restriction of metadata
to the matched left entities followed by renaming into the right id space.

This covers molecule reindexing and the correspondences returned by
combination, extraction, and splitting. A split correspondence maps the new
component to the source molecule, so source metadata is transported to the
component through its reverse.

`ReactionMetadata` does not accept a `MoleculeCorrespondence`: it spans an lhs
namespace and a delta-added namespace, which are not jointly described by one
existing molecule correspondence. Reaction metadata remapping should be added
only together with a corresponding reaction reindexing operation and mapping
type, rather than assigning molecule-correspondence semantics to it.

## Python correspondence operations

The current Python `Correspondence` and `MoleculeCorrespondence` values are
inspectable but not operational. `Correspondence` exposes matched pairs,
id-space sizes, and unmatched ids; `MoleculeCorrespondence` exposes the eight
family correspondences. Neither supports lookup, reversal, composition, or
totality.

Keep both types return-only, but add:

```python
len(correspondence)
correspondence.right_of(left)
correspondence.left_of(right)
correspondence.is_total()
correspondence.reverse()
correspondence.compose(other)

molecule_correspondence.is_total()
molecule_correspondence.reverse()
molecule_correspondence.compose(other)
```

`Correspondence` should wrap `umol_graph_core::Correspondence<usize>`
internally instead of copying its three fields into a parallel Python-only
representation. Expose the properties as `matched_pairs`, `left_count`,
`right_count`, `left_unmatched`, and `right_unmatched`. A separate iterator is
unnecessary because `matched_pairs` already supplies the natural iterable.

Python composition validates that the intermediate id-space sizes agree and
raises `ValueError` when they do not. For `MoleculeCorrespondence`, this check
applies to all eight families before composition. This is boundary validation
of the Rust operation's documented precondition. Construction remains
unexposed: no current Python operation needs it, and accepting arbitrary
matched-pair lists or eight independently constructed family correspondences
would require additional coherence validation.

These operations make metadata remapping usable for all existing molecule
operations. In particular, `split()` returns component-to-source
correspondences, so source metadata is mapped to a component with
`metadata.remap(correspondence.reverse())`.

The correspondence vocabulary should be made suitable for the broader
chemistry-facing API before these operations are exposed:

- `mates` becomes `matched_pairs`;
- `mate_count` becomes `matched_pair_count`;
- `left_exposed` becomes `left_unmatched`;
- `right_exposed` becomes `right_unmatched`;
- `edge_mates` becomes `edge_matched_pairs`.

This rename applies to `Correspondence` and its consumers across Rust and
Python. It does not apply to `Matching::mate(node)` or algorithm-local mate
arrays in the matching implementations, where `mate` is the established
graph-theory term for the matching partner of one vertex.

## Metadata-bearing DSL wrappers

Only three DSL wrappers pair an AST with persistent metadata:

- `MoleculeDsl` with `MoleculeMetadata`;
- `ReactionDsl` with `ReactionMetadata`;
- `ReactionSpanDsl` with `MoleculeMetadata`.

Entity DSL wrappers carry no metadata and are unaffected by this constructor
policy.

The wrappers' private fields preserve coherence after construction, but their
current public, non-validating `from_parts` constructors allow arbitrary
AST/metadata pairs to bypass that invariant. Keep `from_parts` as the private
low-level constructor and add a checked public `new`:

```rust
impl MoleculeDsl {
    pub fn new(
        ast: MoleculeAst,
        metadata: MoleculeMetadata,
    ) -> Result<Self, MetadataError>;

    fn from_parts(ast: MoleculeAst, metadata: MoleculeMetadata) -> Self;
}
```

Apply the same split to `ReactionDsl` and `ReactionSpanDsl`. A Rust `new`
constructor may return `Result`; `try_new`, `new_checked`, and `checked_new`
are unnecessary when there is no public infallible or unchecked `new`.

Parsing uses private `from_parts` because the context establishes coherence.
`FromAst` remains infallible because it constructs empty metadata. External
callers use `new`, which validates both the metadata's intrinsic namespace
invariants and its coherence with the supplied AST:

- molecule keywords must name existing entities of the corresponding kind;
- reaction lhs keywords must name lhs entities;
- reaction delta keywords must name entities introduced by the corresponding
  add deltas;
- keyword and alias disjointness and alias bijectivity must hold across the
  complete molecule or reaction scope.

`into_parts` remains public because it consumes an already coherent wrapper.

## DSL terminology

Within `umol-ast`, `Id` means a numerical AST identifier such as `AtomId` or
`BondId`. A symbolic EDN value such as `:carbonyl` is a keyword, even when it
appears as the value of the DSL's literal `:id` field.

Apply the following terminology changes throughout `umol-ast/src/dsl`:

| Current | Replacement |
| --- | --- |
| `AtomRef::Id(String)` and the other entity-ref variants | `AtomRef::Keyword(String)` and corresponding variants |
| `contains_id()` | `contains_keyword()` |
| `ParseError::DuplicateId` | `ParseError::DuplicateKeyword` |
| metadata fields such as `atom_ids` | the entity-kind fields inside `keywords` or `delta_keywords` |
| `id: Option<String>` in parsed entry structs | `keyword: Option<String>` |
| string locals named `id` | `keyword` |
| “symbolic id” or “id keyword” | “keyword reference” or “keyword” |
| `value::id`, which parses variable names | `variable_name` |

The helper that reads an optional literal `:id` field should be named
`optional_id_keyword`: it retains the connection to the surface key while
making the returned value explicit. It populates an entry's `keyword` field.

The following uses of `id` remain correct:

- `AtomId`, `BondId`, and the other numerical AST identifier types;
- variables and fields whose values have one of those types;
- delta fields such as `AtomDelta::Add { id, .. }`;
- numerical entity “ID space” terminology;
- the literal DSL map key `:id`;
- parser and renderer code matching or emitting `Edn::keyword("id")`.

Test case names should follow the same distinction. A case specifically
testing the literal `:id` map field may retain `id`; a case testing an entity
reference or lookup uses `keyword`.

## Python implications

Do not bind the current metadata containers before this reconciliation.
Python should receive the stable `MoleculeMetadata` and `ReactionMetadata`
objects, with mapping operations rather than exposure of their physical Rust
fields. This permits new metadata categories without freezing the internal
layout.

The ordinary AST-only operations remain:

```python
MoleculeAst.parse(text, *, defaults=None) -> MoleculeAst
molecule.render(*, defaults=None) -> str
```

Metadata-preserving operations are parallel and explicit:

```python
MoleculeAst.parse_with_metadata(
    text,
    *,
    defaults=None,
) -> tuple[MoleculeAst, MoleculeMetadata]

molecule.render_with_metadata(
    metadata,
    *,
    defaults=None,
) -> str
```

Reaction parsing and rendering use the same signatures with
`ReactionDefaults` and `ReactionMetadata`. Tuple unpacking is ordinary,
long-standing Python syntax:

```python
molecule, metadata = MoleculeAst.parse_with_metadata(text)
```

`str(molecule)` is equivalent to `molecule.render()` and therefore renders
without retained metadata. It uses positional references where no keyword is
available.

Detached metadata may become incompatible with a modified AST.
`render_with_metadata` constructs the relevant DSL wrapper through its checked
`new` operation rather than implementing a separate renderer-specific
validation path.

## Completion criteria

- Contexts and metadata share one bidirectional keyword representation.
- `MoleculeContext` and `ReactionContext` replace the current concrete
  namespace names.
- Metadata lookup works efficiently in both directions for all eight entity
  kinds.
- `MoleculeMetadata` remaps entity keywords through a
  `MoleculeCorrespondence` while preserving aliases.
- Python correspondence values support the lookup, reversal, composition, and
  totality operations required to transport metadata.
- Metadata and context construction enforce the same keyword and alias
  invariants.
- The three metadata-bearing DSL wrappers have private non-validating
  `from_parts` constructors and public checked `new` constructors.
- `ReactionMetadata` renders the lhs/delta union without building a combined
  clone.
- DSL code uses `id` only for numerical AST identifiers or the literal `:id`
  surface field.
- Python metadata-preserving parse and render APIs are defined against the
  reconciled Rust types.

## Staged implementation plan

`MetadataError` is extended incrementally. Each subitem introduces only the
variants required by the invariant implemented in that subitem; the error enum
is not populated in advance with hypothetical failure modes.

### S0 — Bidirectional keyword foundation

**S0a — `Entity` ordering.** In `umol-ast/src/ast/entity.rs`, derive
`PartialOrd` and `Ord` for `Entity` and test ordering across variants and ids.
This is additive and stays green. **Implemented (green).** [dep: none]

**S0b — Public correspondence vocabulary.** In
`umol-graph-core/src/correspondence.rs`, rename `mates` to `matched_pairs`,
`mate_count` to `matched_pair_count`, `left_exposed` to `left_unmatched`, and
`right_exposed` to `right_unmatched`. Rename `edge_mates` to
`edge_matched_pairs`. Migrate all graph-core, AST, graph, and Python consumers,
tests, documentation, fields, parameters, and local variables whose values are
correspondence pairs. Do not rename
`Matching::mate(node)`, its mate array, or matching-algorithm locals that
represent one vertex's matching partner. This is breaking and goes red→green
across the workspace. **Implemented (green).** [dep: none]

**S0c — Entity-level correspondence lookup.** In
`umol-ast/src/ast/correspondence.rs`, add
`MoleculeCorrespondence::right_of(Entity)` and `left_of(Entity)`, dispatching
to the existing per-family correspondences and preserving the input entity
kind. Test both directions for all eight variants and unmatched entities on
both sides. This is additive and stays green. **Implemented (green).**
[dep: S0a, S0b]

### S1 — DSL keyword terminology

**S1a — Entity reference variants.** In `umol-ast/src/dsl/refs.rs`, rename the
string-bearing `Id` variant of each of the eight entity-reference enums to
`Keyword`, then migrate parsing, resolution, denotation, rendering, fixtures,
and property strategies in the same change. Numerical id variants and
structural reference variants are unchanged. This is breaking and goes
red→green while all enum consumers are migrated. [dep: none]

**S1b — Parsed-entry fields and the literal `:id` helper.** Across the molecule,
reaction, reaction-span, overlay, constraint, and delta DSL input structs,
rename `Option<String>` fields and local variables from `id` to `keyword`.
Rename the helper that reads the literal `:id` field to
`optional_id_keyword`; the surface key remains `:id`. Update focused parser
tests for absent, valid, and malformed literal `:id` values. This is breaking
and goes red→green with all constructors and pattern matches migrated.
[dep: S1a]

**S1c — Namespace and parser vocabulary.** Rename `contains_id` to
`contains_keyword`, `ParseError::DuplicateId` to `DuplicateKeyword`, and
`value::id` to `variable_name`. Update prose, diagnostics, fixtures, and test
case names so `id` remains only for numerical AST ids, numerical id spaces, or
the literal `:id` field. Preserve the existing parse-error behavior apart from
the corrected terminology. This is breaking and goes red→green with all
callers migrated. [dep: S1b]

### S2 — Persistent metadata representation

**S2a — Shared render query surface.** Move the `Metadata` trait from
`dsl/molecule.rs` into `dsl/metadata.rs` and update private molecule rendering
helpers that only denote entity references to accept `&impl Metadata`.
Alias-specific whole-molecule rendering remains on `MoleculeMetadata`.
Re-export the public trait from the same DSL surface as today and test that
the generic reference renderers produce the existing positional and keyword
forms. This is an internal rewire and stays green after its caller migration.
[dep: S0a, S1c]

**S2b — `MoleculeMetadata`.** Replace the eight one-way maps in
`MoleculeMetadata` with
`keywords: BiBTreeMap<Entity, String>`, retaining the atom-alias bimap.
Introduce and re-export `MetadataError` with only the variants needed by this
subitem. Add efficient `*_keyword(id)` and `*_id(keyword)` methods for all
eight entity kinds. Make the existing setter and builder mutation surfaces
fallible and atomic, enforce keyword/alias disjointness and alias bijectivity,
and test all eight `Entity` variants, idempotence, rebinding, cross-kind
collisions, collision rollback, and empty/default construction. Migrate unit
tests, fixtures, macros, and property strategies to the fallible API. This is
breaking and goes red→green with every Rust caller migrated. [dep: S0a, S2a]

**S2c — Molecule metadata remapping.** In `dsl/metadata.rs` and
`dsl/molecule.rs`, add public consuming
`MoleculeMetadata::remap(self, &MoleculeCorrespondence) -> Self`. Rebuild the
entity-keyword bimap with matched targets mapped left-to-right, omit unmatched
left entity bindings, and move the atom-alias bimap through unchanged.
Test identity, a nontrivial permutation, partial restriction, total reverse
roundtrip, composition, every entity kind, and alias invariance. This is
additive and stays green. [dep: S0c, S2b]

**S2d — `ReactionMetadata`.** Replace the eight delta maps with
`delta_keywords: BiBTreeMap<Entity, String>`, retain
`lhs: MoleculeMetadata` and the reaction alias bimap, and implement `Metadata`
directly with unambiguous delta-then-lhs lookup. Provide explicit lhs/delta
mutation operations rather than exposing either physical field as the public
mapping API. Enforce keyword and alias uniqueness over the complete reaction
scope and test all eight entity kinds in both scopes, including collisions
across the scope boundary. Extend `MetadataError` only for newly introduced
reaction-scope failures. This is breaking and goes red→green with all callers
migrated. [dep: S2b]

**S2e — Clone-free reaction rendering.** In `dsl/reaction.rs`, pass
`&ReactionMetadata` through the generic reference-rendering helpers and remove
`combined_metadata()` together with the lazy combined clone in
`render_reaction_edn`. Preserve exact reaction rendering for positional,
lhs-keyword, delta-keyword, constraint, and overlay references. Add regression
tests that exercise each reference scope. This is a breaking removal of an
obsolete public method and goes red→green with all callers migrated.
[dep: S2d]

### S3 — Parse-time contexts

**S3a — `MoleculeContext`.** In `dsl/namespace.rs`, rename
`MoleculeNamespace` to `MoleculeContext` and remove the duplicate
keyword-to-id storage from its per-kind registries. Let the context own/build
`MoleculeMetadata` alongside counts and participant indexes; keyword and alias
registration updates metadata atomically. Replace the inversion-based
`From<&MoleculeNamespace> for MoleculeMetadata` projection with extraction of
the already-built metadata. Migrate molecule and subpattern parsing, public
exports, tests, and documentation in the same change. Test count allocation,
participant lookup, keyword lookup, collision rollback, and parsed metadata
equivalence between tree and streaming parsers. This is breaking and goes
red→green across the workspace. [dep: S2b]

**S3b — `ReactionContext`.** In `dsl/reaction.rs`, rename
`ReactionNamespace` to `ReactionContext`. Build `ReactionMetadata` directly
from its lhs and delta contexts, preserve continuation id allocation, and
route layered keyword, alias, participant, and entity-scope queries through
the shared storage. Remove projection/inversion code and migrate reaction
parsing, tests, and public exports. Test lhs/delta lookup, continuation
indices, cross-scope collisions, added-entity references, and tree/streaming
metadata parity. This is breaking and goes red→green across the workspace.
[dep: S3a, S2d]

**S3c — Reaction-span parsing.** Rewire `dsl/reaction_span.rs` to use
`MoleculeContext` and to move the context's `MoleculeMetadata` into the
resulting wrapper. Update its unit and roundtrip tests for all eight entity
kinds and aliases. This is an internal rewire and stays green after caller
migration. [dep: S3a]

### S4 — Checked metadata-bearing wrappers

**S4a — `MoleculeDsl::new`.** Add the public fallible
`MoleculeDsl::new(ast, metadata)` constructor. Validate that every metadata
binding names an existing entity of the corresponding kind; add only the
required out-of-range `MetadataError` representation. Keep `from_parts`
temporarily available so the additive constructor can land green. Test every
entity kind, boundary indices, empty metadata, and parsed coherent parts.
This is additive and stays green. [dep: S2b]

**S4b — `ReactionDsl::new`.** Add the public fallible
`ReactionDsl::new(ast, metadata)` constructor. Validate lhs metadata against
the lhs and each delta keyword against an `Add` delta of the corresponding
kind, while reusing the intrinsic scope checks already enforced by
`ReactionMetadata`. Add error variants only for failure modes not expressible
after S4a. Test all eight add-delta kinds, wrong-kind and absent additions,
lhs/delta scope errors, and coherent parsed parts. This is additive and stays
green. [dep: S2d, S4a]

**S4c — `ReactionSpanDsl::new`.** Add the public fallible
`ReactionSpanDsl::new(ast, metadata)` constructor and validate each keyword
against the corresponding span entity collection. Reuse the S4a error
representation where possible. Test all eight kinds, empty metadata, and
coherent parsed parts. This is additive and stays green. [dep: S2b, S4a]

**S4d — Privatize unchecked construction.** Make all three `from_parts`
constructors private. Keep parser paths on private `from_parts`, because their
contexts establish coherence, and keep `FromAst` infallible with empty
metadata. Migrate macros, ordinary Rust callers, fixtures, and property
strategies to checked `new`; keep `into_parts` public. Property generators that
intend to produce wrappers must generate coherent metadata, while dedicated
invalid-metadata strategies assert the exact `MetadataError`. This is breaking
and goes red→green with all external callers migrated in the same subitem.
[dep: S4a, S4b, S4c, S3b, S3c]

### S5 — Python correspondence and metadata types

**S5a — Operational correspondence values.** In
`umol-py/src/correspondence.rs`, replace the copied-field representation of
Python `Correspondence` with a wrapper over
`umol_graph_core::Correspondence<usize>`. Preserve the existing properties and
add `__len__`, `right_of`, `left_of`, `is_total`, `reverse`, and `compose`.
Add `is_total`, `reverse`, and `compose` to `MoleculeCorrespondence`.
Composition checks intermediate counts and raises `ValueError` on mismatch;
the molecule-level operation validates all eight families before composing.
Keep both classes frozen and return-only. Test empty, partial, total, reversed,
composed, and dimension-mismatched correspondences, plus every molecule entity
family. This is additive apart from the internal representation change and
stays green. [dep: S0c]

**S5b — Metadata wrappers and exception mapping.** Add
`umol-py/src/metadata.rs` with Python `MoleculeMetadata` and
`ReactionMetadata` wrappers over the reconciled Rust types, and register them
in `_native`. Expose per-kind id-to-keyword and keyword-to-id operations and
collision-aware mutation methods without exposing Rust fields. Preserve atom
aliases inside the wrappers for roundtripping; do not add a separate,
half-designed alias-authoring surface in this stage. Add a Python
`MetadataError` exception and map Rust `MetadataError` values to it. Test all
eight lookup pairs, lhs/delta reaction scope, mutation collisions, exception
types/messages, repr, and molecule metadata remapping through the existing
Python `MoleculeCorrespondence`. This is additive and stays green.
[dep: S2c, S4d, S5a]

### S6 — Python DSL parse and render operations

**S6a — Molecule operations.** In `umol-py/src/molecule.rs`, retain
`MoleculeAst.parse(text, *, defaults=None)`, add
`parse_with_metadata`, `render`, and `render_with_metadata`, and make
`__str__` delegate to `render()` with default `MoleculeDefaults`.
`parse_with_metadata` returns the AST after applying defaults together with
the parsed metadata. `render_with_metadata` first lowers the AST under the
requested defaults, then constructs `MoleculeDsl` through checked `new`.
Test keyword- and alias-bearing roundtrips, explicit defaults, positional
metadata-free rendering, detached incompatible metadata, keyword-only
arguments, and `str(molecule) == molecule.render()`. This is additive and
stays green. [dep: S5b]

**S6b — Reaction operations.** In `umol-py/src/reaction.rs`, add the parallel
`ReactionAst.parse_with_metadata`, `render`, and `render_with_metadata`
operations using `ReactionDefaults`, and change `ReactionAst.__str__` to
delegate to `render()`. Route metadata incompatibility through the same checked
`ReactionDsl::new` and Python `MetadataError` mapping. Test lhs and delta
keywords across all entity families, reaction aliases, explicit defaults,
metadata-free positional rendering, incompatible metadata, keyword-only
arguments, and `str(reaction) == reaction.render()`. This is additive except
for the corrected `__str__` implementation; it goes red→green with Python
expectations updated in the same subitem. [dep: S5b, S6a]

### S7 — Cross-cutting verification and documentation

**S7a — Metadata properties.** In the `umol-ast` property suite, add explicit
properties for bidirectional lookup, atomic failed insertion, namespace
disjointness, context-to-metadata preservation, checked-constructor
acceptance of parsed parts, rejection of incoherent generated parts, and
parse/render/reparse preservation for molecule, reaction, and reaction-span
DSLs. Add metadata-remapping properties for identity, composition, total
reverse roundtrip, partial restriction, and alias invariance. Keep intentional
overlap with unit tests documented beside the property modules. This is
additive and stays green. [dep: S2c, S4d]

**S7b — Terminology and API audit.** Update the DSL module documentation,
`umol-ast/spec/umol-dsl-spec.md`, Python docstrings, and related discussion
forward references. Audit the workspace for stale concrete `*Namespace`
names, symbolic `Id` variants, `contains_id`, `DuplicateId`, one-way metadata
maps, public `from_parts` constructors, and `combined_metadata`. Run Rust
formatting, focused and workspace tests, property tests, clippy, the Python
3.13 binding build, and Python tests before marking this document complete.
This is documentation plus verification and stays green. [dep: S1c, S2e,
S3b, S3c, S4d, S6b, S7a]

The critical path is
`S0 → S1 → S2 → S3 → S4 → S5 → S6 → S7`.
No stage is deferrable for the metadata-preserving Python API: Rust metadata
must be coherent and stable before it is bound, and the final property and
documentation audit is part of the contract rather than optional cleanup.
