# AST comparison and property-suite organization

## Scope

This document refines the comparison API and property-test portions of the
validated reaction-application follow-up. It records the relation between the
comparison operations that already exist and the molecule-level operation
needed by application tests, then gives the implementation sequence for
comparison, panic removal, and property-suite cleanup.

## Comparison vocabulary

The current APIs are not interchangeable:

| Operation | Domain | Meaning |
| --- | --- | --- |
| `==` | every AST value, including `MoleculeAst` | Derived representation equality. The stored tree, ordering, and presence/absence details must agree exactly. It is intentionally cheap and is not the semantic comparison used by canonicalization. |
| `canonical_eq` | values implementing `Canonicalize` | Semantic equality of one value type after canonicalization. Equivalent encodings of the same entity compare equal; unsatisfiable values follow the existing canonicalization contract. It does not compare molecule topology or remap relation participants. |
| `equiv` / `equiv_under` | relation data implementing `Equiv`/`BiEquiv` | Relation-data comparison. The value axis is `canonical_eq`; `equiv_under` additionally reindexes participant positions into a supplied frame. It is about relation data, not whole-molecule identity. |
| `equiv` | `MoleculeAst` | Complete semantic equality in the existing ID and participant frame: identical topology/incidence, with entity ASTs and constraints compared canonically. |
| `equiv_under` | `MoleculeAst` plus a total `MoleculeCorrespondence` | Complete semantic equality after translating the left molecule into the right molecule's ID and participant frame. |

There is no molecule-level `canonical_eq` in this round. Full canonical
equality remains dependent on full molecule canonical labeling. `==` remains
exact stored-representation equality. Molecule `equiv` and `equiv_under` extend
the existing relation-data vocabulary: `equiv` compares in the current frame;
`equiv_under` compares after applying the supplied frame correspondence.

## Staged implementation plan

### S0 — Framed semantic equality for `MoleculeAst`

- **S0a — complete-correspondence operations**
  (`umol-ast/src/ast/correspondence.rs`, `remap.rs`): add
  `MoleculeCorrespondence::is_total`, requiring all eight constituent
  correspondences to be total, and `MoleculeCorrespondence::to_remapping`,
  producing the existing complete `IdRemapping` from a correspondence that is
  total on the left. Follow the existing `Correspondence::is_total` and
  `GraphCorrespondence::to_remapping` contracts. Table tests cover one exposed
  entity in each family and a complete non-identity correspondence; assertions
  use exact booleans and remappings. `IdRemapping` is now part of the public
  `ast` surface so the public conversion has a nameable return type.
  **Implemented (green).** `[dep: —]`
- **S0b — same-frame `MoleculeAst::equiv`**
  (`umol-ast/src/ast/molecule.rs`): compare the ordinary graph and the stored
  participants/sites of all six overlay/stereo families in their current ID
  frame; compare atom, bond, overlay, and stereo ASTs with `canonical_eq`,
  position-dependent relation data with its existing `equiv`, and molecule
  constraints with `canonical_eq`. Return `false` for any count, topology,
  participant, site, entity-AST, or constraint difference. Do not construct a
  correspondence or derive deltas. Table tests vary one dimension at a time
  across all entity families; property tests cover reflexivity, symmetry, and
  agreement with `==` for already-canonical molecules.
  **Implemented (green).**
  `[dep: —]`
- **S0c — correspondence-aware `MoleculeAst::equiv_under`**
  (`umol-ast/src/ast/molecule.rs`, `correspondence.rs`): require S0a's total
  correspondence, verify that every mapped bond endpoint, overlay participant,
  stereo site, and stereo ligand agrees with the right molecule, and compare
  each mated entity AST canonically. Reindex position-dependent relation data
  with the participant permutations already used by overlay matching; remap
  the left molecule constraints through `IdRemapping` before canonical
  comparison. Return `false` for a partial or structurally inconsistent
  correspondence. Tests cover non-identity atom numbering, each overlay family,
  ordered dative/stereo frames, constraint references, partial
  correspondences, and deliberately inconsistent correspondences. Property
  tests cover symmetry under `correspondence.reverse()` and reduction to S0b
  under the identity frame. **Additive (green).** `[dep: S0a, S0b]`
- **S0d — comparison contract and application-property migration**
  (`umol-ast/tests/property/reaction.rs`): record the comparison matrix for
  `==`, entity `canonical_eq`, relation
  `equiv`/`equiv_under`, and molecule `equiv`/`equiv_under`. Replace complete
  result inspection in the eight focused host-refinement application properties
  with `MoleculeAst::equiv`. Run the focused properties and the 4,096-case
  reaction-application soak.
  **Additive (green).** `[dep: S0c]`

### S1 — Replace application-path panics with typed outcomes

- **S1a — panic inventory and error taxonomy** (`umol-ast/src/ast/reaction.rs`,
  `error.rs`, application validators): enumerate every `expect`/panic reachable from reaction
  preflight, matching, lowering, and derivation construction; classify each as
  invalid input/precondition or internal failure and assign a stable
  `ApplyPreconditionError`/`ApplyError` variant. Add one regression test per
  classification. **Additive (green).**
- **S1b — reference and incidence preflight** (`reaction.rs`, application
  validators): validate all
  delta IDs, entity references, endpoints, overlay participants, stereo sites,
  and constraint references before entering the iterator. Return the typed
  precondition error; never use a lookup `expect` for caller-controlled data.
  **Additive (green).** [dep: S1a]
- **S1c — checked lowering and iterator termination** (`reaction.rs`): replace
  remaining lowering/reframing assumptions with checked lookups. Convert
  impossible-after-preflight failures into one fatal `ApplyError`, emit it once,
  and terminate permanently; preserve local match rejection behavior.
  **Breaking internal migration (red→green).** [dep: S1b]
- **S1d — malformed-input property and fuzz coverage**
  (`tests/property/malformed.rs`, `fuzz/fuzz_targets/fuzz_reaction.rs`): generate invalid IDs, missing incidence, incompatible
  overlay references, and malformed update combinations; assert errors rather
  than panics and assert post-fatal iterator termination. **Additive (green).**
  [dep: S1c]

### S2 — Property-suite structure and purpose

- **S2a — invariant inventory** (`umol-ast/tests/property`): write a table
  mapping each property to its invariant, generator domain, oracle, and failure
  class. Explicitly distinguish identity/canonicalization, update/difference,
  delta composition, reaction application, malformed-input safety, and
  serialization/span properties. **Additive (green).**
- **S2b — split the large reaction property module**
  (`umol-ast/tests/property.rs`, `tests/property/reaction.rs`): retain the
  existing `strategies.rs`, `value.rs`, `lattice.rs`, `entity.rs`, `stereo.rs`,
  `molecule.rs`, `edit.rs`, `delta.rs`, and `substructure.rs` modules. Replace
  `reaction.rs` with `reaction_application.rs`, `reaction_composition.rs`,
  `reaction_span.rs`, and `reaction_serialization.rs`; register
  `malformed.rs` from S1d alongside them. Move each existing reaction property
  according to the operation named in the test and change no generators or
  assertions in this subitem. Keep the single `property` test target.
  **Breaking file migration (red→green).** [dep: S2a]
- **S2c — remove accidental overlap while retaining deliberate overlap**:
  preserve one minimal identity property and one canonical-equivalence
  property per entity family; retain the eight focused host-refinement
  application properties because they exercise distinct delta variants; remove
  duplicate generators/assertions whose only difference is incidental data
  shape. Record the purpose of any retained overlap beside the property.
  **Additive (green).** [dep: S2b]
- **S2d — migrate semantic assertions to framed equivalence**: use molecule
  `equiv` for complete same-frame results, molecule `equiv_under` where a
  correspondence is part of the property, entity `canonical_eq` for
  entity-focused results, and `==` only where representation identity is the
  invariant. Keep direct topology assertions only for properties explicitly
  testing topology. No property may reach into private molecule storage merely
  to obtain a general-purpose equality check. **Additive (green).** [dep: S0d, S2b]
- **S2e — suite and soak gate**: run the default property suite, targeted
  reaction application properties at the enlarged case count, malformed-input
  properties, clippy, formatting, and `git diff --check`; retain minimized
  regressions only when they exercise a documented invariant. **Additive
  (green).** [dep: S1d, S2c, S2d]

## Dependencies and deferrals

The comparison critical path is `S0a → S0c → S0d → S2d`, with `S0b → S0c`
joining it. Panic removal proceeds
independently as `S1a → S1b → S1c → S1d`. The property-suite restructuring
follows `S2a → S2b → {S2c, S2d} → S2e` and joins the panic path at S2e.

Full molecule `canonical_eq`, graph-plus-overlays canonical labeling, scalable
`umol-perm`, and canonical-image search are deferred. They are not dependencies
of the framed comparison or property-suite work. Module splitting and duplicate
cleanup are organizational and can be deferred after S1 if necessary, but the
final suite gate remains required.
