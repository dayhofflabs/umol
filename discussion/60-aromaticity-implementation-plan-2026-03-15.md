# Aromaticity Implementation Plan

Based on discussion/59, the user's TODOs, and current code review.

## Core principle

**Resolution preserves input representation.** Kekule stays Kekule. Aromatic-hinted
input gets aromatic treatment. Transformations (perceive/kekulize) are explicit,
separate operations on resolved Molecules.

## Current state

- `AtomBuilder` carries `aromatic_hint: Option<bool>` from SMILES/MOL parsing
- `BondBuilder` carries `aromatic_hint: Option<bool>` (MOL aromatic bond type, SMILES `:`)
- `AtomTypeQuery::from_builder_atom` handles `Some(false)` → constrain to non-aromatic,
  but `Some(true)` and `None` both map to unconstrained (broken)
- Bond aromatic hints are not consulted during atom hint inference
- `ValenceValidator::Counts` hardcodes `AromaticValence::None` — cannot represent aromatic atoms
- Aromaticity phase adds `AromaticSystem` objects but does not narrow atom candidates
- No bond hint validation at build time
- No transform module exists (`kekule.rs` exists but not as a Molecule→Builder transform)

## Phases

### Phase 1: AromaticConstraint and valence query fix

**Goal**: Aromatic SMILES (`c1ccccc1`) resolves correctly under AtomTyping strategy.

1a. Add `AromaticConstraint` enum to `atom_type.rs`:

```rust
pub enum AromaticConstraint {
    None,       // must be non-aromatic (a absent)
    Any,        // must be aromatic (a present, any value)
    Valence(u8) // must be aromatic with specific valence
}
```

1b. Change `AtomTypeQuery.aromatic_valence: Option<AromaticValence>` →
    `aromatic_constraint: Option<AromaticConstraint>`.
    Update `matches()`:
    - `None` → unconstrained (match anything)
    - `Some(AromaticConstraint::None)` → match only `AromaticValence::None`
    - `Some(AromaticConstraint::Any)` → match only `AromaticValence::Valence(_)`
    - `Some(AromaticConstraint::Valence(n))` → match only `AromaticValence::Valence(n)`

1c. Bond→atom hint propagation in `AtomTypeQuery::from_builder_atom`:
    If `aromatic_hint == Some(true)` OR any incident bond has `aromatic_hint == Some(true)` →
    `aromatic_constraint = Some(AromaticConstraint::Any)`.
    If `aromatic_hint == Some(false)` → `Some(AromaticConstraint::None)`.
    Otherwise → `None`.
    Requires access to incident bonds in `from_builder_atom` (already takes `builder` + `atom_index`).

1d. Tests: benzene from aromatic SMILES resolves with `{Cv2a1H}` candidates.
    Pyridine, pyrrole, furan from aromatic SMILES.
    Mixed: `c1ccccc1CC` — aromatic ring + aliphatic tail.

### Phase 2: Aromaticity phase candidate narrowing

**Goal**: After aromaticity detection, atom candidates are narrowed to be consistent
with the detected AromaticSystem assignments.

2a. In `resolve_aromaticity_with`, after `model.aromatic_systems()` returns:
    For each atom in an `AromaticSystem`, narrow candidates to those whose
    `aromatic_valence` matches the system's contribution for that atom.
    For each atom with `aromatic_hint == Some(true)` NOT in any system → error
    (per hint policy).

2b. Add `AromaticityHintPolicy` to config: `Strict` (error) / `Lenient` (warn) / `Ignore`.

2c. Bond hint validation: for each bond with `aromatic_hint == Some(true)`,
    verify both endpoints are in some `AromaticSystem`. Apply policy.

### Phase 3: Counts strategy aromatic support

**Goal**: `ValenceStrategy::Counts` can represent aromatic atoms.

3a. When `aromatic_hint == Some(true)` (or inferred from bond hints):
    Infer aromatic valence from element + charge + sigma-skeleton bond order sum.
    The registry isn't available but the inference is usually unambiguous
    (C: a=1, N-pyridine: a=1, N-pyrrole: a=2, O-furan: a=2, S: a=2).

3b. Adjust implicit H calculation to account for aromatic pi-electron budget.

3c. Emit `AtomTypeSpec` with correct `AromaticValence::Valence(n)` instead of `::None`.

### Phase 4: Separate aromaticity primitives from resolver

**Goal**: Aromaticity detection code works on both MoleculeBuilder (resolution)
and Molecule (transformation).

4a. Move `resolver/aromaticity/` → `graph_ir/aromaticity/` (or similar).
    The models (hueckel_rule, hmo, clar) are not resolver-specific.
    `AromaticityModel` becomes a standalone dispatch.

4b. `RingEnumerator` enum lives alongside rings.rs:
    - `new(strategy: &RingStrategy)` → constructs from config
    - `enumerate_builder(&self, builder) → MoleculeRings`
    - `enumerate_molecule(&self, molecule) → MoleculeRings`
    Both use the same underlying ring enumeration; PiSubgraph filtering
    uses aromatic candidates (builder) or resolved aromatic types (molecule).

4c. Resolver's aromaticity phase becomes a thin caller of these primitives.

### Phase 5: Transform module

**Goal**: Explicit Molecule → MoleculeBuilder transformations.

5a. New module `graph_ir/transform/` with:
    - `aromatize(mol: &Molecule, config) → Result<MoleculeBuilder, _>`
      Perceives aromaticity on a Kekule molecule. Runs ring enumeration +
      aromaticity model on the resolved structure. Produces builder with
      AromaticSystem objects attached.
    - `kekulize(mol: &Molecule, config) → Result<MoleculeBuilder, _>`
      Assigns alternating single/double bonds to aromatic systems.
      Produces builder with explicit bond orders, no AromaticSystem.

5b. No-op guard: if `aromatize` is called on a molecule that already has
    AromaticSystem objects, return early or error (configurable).

5c. Round-trip tests: aromatic → kekulize → re-resolve → aromatize → compare.

### Phase 6: Config cleanup

6a. Discuss naming: replace `enabled` flag with behavior-based config
    (the aromaticity phase always runs; behavior depends on whether hints
    are present and what the policy says).

6b. `AromaticityHintPolicy` (from Phase 2).

6c. Align discussion docs with actual API.

## Dependencies

```
Phase 1 ← standalone (highest priority, unblocks aromatic SMILES)
Phase 2 ← Phase 1 (needs aromatic candidates to narrow)
Phase 3 ← Phase 1 (uses same AromaticConstraint)
Phase 4 ← Phase 2 (needs working aromaticity to refactor)
Phase 5 ← Phase 4 (uses extracted primitives)
Phase 6 ← can happen alongside any phase
```

## Open questions

- Phase 4 module location: `graph_ir/aromaticity/` vs keep under `resolver/aromaticity/`
  but with public interface?
- RingEnumerator: PiSubgraph on Molecule — what predicate? Resolved `aromatic_valence`?
  Alternating bond pattern? Both?
- Phase 6 naming: what replaces `enabled` + `input_mode`?
