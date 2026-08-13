# AST API worklist

Status: **Proposed**
Date: 2026-07-27
Relates: [086](086-molecule-ast-api-2026-04-16.md),
[104](104-stereochemistry-implementation-plan-2026-05-31.md),
[105](105-dsl-fixes-2026-06-06.md),
[123](123-ast-allocation-survey-2026-06-21.md),
[149](149-molecule-ring-cache-and-hashing-2026-07-13.md)

## Scope

This document tracks unresolved work on the semantic AST and its public API. It
does not revive the chemist-facing `Molecule`/`Pattern` wrapper hierarchy from
doc 086. `MoleculeAst` and `ReactionAst` remain the public semantic models.

Variable representation, entity-model expansion, allocation optimization, and
hashing remain in docs 115, 117, 123, and 149 respectively.

## Public API organization

- Break long `MoleculeAst` and entity implementation blocks into modules grouped
  by operation: construction, views, comparison, remapping, substructure, and
  reaction operations. Do not split solely to reduce file length.
- Audit mutable collection iterators. Prefer mutation APIs that preserve
  container invariants; determine separately whether individual entity
  `*_mut` accessors remain sound.
- Document iteration order for public `ids`, `iter`, incidence, neighbor, and
  relation queries where callers can observe it.
- Keep external view APIs on semantic accessors. Reaching through `raw_graph`
  is reserved for graph operations that have no AST-level equivalent.
- Audit constructor argument types such as bond order. Public arguments should
  express the accepted domain and convert once at the boundary.
- Inventory public operations duplicated across entity kinds and determine
  which are better expressed by one `Entity`-based operation. Use `Entity`
  where the contract and return type are uniform across all eight kinds; retain
  kind-specific operations where their typed arguments or results prevent
  invalid calls.
- Review imports inside the `ast` and `dsl` trees after the module split:
  `super` for a direct parent and `crate` for non-local dependencies.

## Naming and value semantics

- Decide whether the `ValueAst` name still describes the numeric predicate AST
  or whether a `NumAst` migration is warranted. This is a cross-cutting rename
  and must not begin without a complete public-name inventory.
- Rename `aromatic_increment` to `aromatic_covalence` if the latter accurately
  describes every caller and formula.
- Replace remaining special-purpose builder equivalence methods with the
  established comparison vocabulary rather than adding another equivalence
  relation.
- Audit `matches_value` and other legacy matching names against `Lattice`,
  `AsLit`, `matches`, and literal extraction.
- Ensure reverse-remapping APIs return the semantic remapping type rather than
  an untyped tuple or raw vector.

## Relation and view integrity

- Replace `Option<Vec<AtomId>>` and `Option<Vec<BondId>>` in molecule-scope
  constraints with semantic scope enums. `All` and `Set` have different
  behavior under structural growth and should be explicit in the type.
- Verify that sorting aromatic-system and multicenter-bond participants applies
  the same permutation to their per-participant electron counts. The stored
  participant/electron association is semantic.
- Audit atom and bond reverse-incidence accessors for stereo atoms and stereo
  bonds. Add missing predicates or iterators only where the underlying
  cardinality and role are well defined.
- Review whether molecule-level constraints need direct entity- and kind-based
  query methods. Add them only if validator and matching consumers otherwise
  repeat scans.

## Ring views

Doc 149 removed the ring cache and introduced `RingViews`, but its proposed
`RingView` cleanup remains unfinished. The current `RingView` still exposes
identifier slices and accepts another `RingView` in relation methods.

The follow-up is:

- give `RingView` the context needed to return atom and bond views;
- provide `atom_ids`/`atoms` and `bond_ids`/`bonds` with the same conventions as
  other entity views;
- make binary ring relations id-keyed through the owning ring set;
- retain `RingSet` primitives for algorithms that operate directly on ids;
- remove view-as-argument APIs rather than duplicating them.

The last two items are folded into doc 194 (S0g): the id-keyed `RingSet::shared_atoms`/
`shared_bonds` already exist, and the view-as-argument `RingView::shared_atoms`/`shared_bonds`
are deleted there. The accessor-convention items above remain here.

This is an AST API change, not part of hashing.

## Substructure matching

Molecule-level pattern constraints (relational and molecule leaves, combinators) are silently
ignored by the matcher; doc 195 owns that design.

### Stereo constraints

Stereo overlays are frame-aligned during substructure matching. A pattern that
expresses stereo only through atom `#T` or bond `#C` constraints is not yet
handled equivalently. Implement the constraint path without grounding the
pattern:

- derive or locate the corresponding pattern stereo element;
- map its ligand frame through the candidate atom correspondence;
- reindex the host coset into the pattern frame;
- apply the constraint lattice;
- cover tetrahedral and cis/trans cases, frame permutations, achiral targets,
  and enantiomer rejection.

### Reference fixtures

Doc 104 planned committed RDKit and ArcMatch reference fixtures but only the
internal algorithm cross-validation and development benchmark scripts landed.
Decide whether the external fixtures add an independent semantic check. If they
do, capture their outputs as repository fixtures; the external libraries must
not become test dependencies. Otherwise remove the requirement explicitly.

## Completion criteria

- The public AST surface is organized by semantic responsibility.
- Relation payloads remain aligned under every construction and remapping path.
- Ring views follow the same id/view conventions as other entity views.
- `#T` and `#C` pattern constraints have the same frame-correct matching
  semantics as stereo overlays.
