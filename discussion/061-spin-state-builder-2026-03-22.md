# SpinStateBuilder Design and Implementation Plan

Revised 2026-03-22 after design review.

## Scope and status

This document locks design decisions and provides an implementation plan for:

1. Unified partial-spin handling via `SpinStateBuilder`
2. Bond spin invariants with full `SpinState` in resolved `Bond`
3. Fallible `from_table_atom` / `from_table_bond`
4. Aromatic systems carrying explicit `SpinState`
5. TableIR atom spin fields split into `unpaired_electrons` and `multiplicity`
6. No high-spin guessing for molecular spin

This is implementation-plan ready.

## Design decisions

### D1. Partial spin state type (to be implemented in `src/spin.rs`)

`SpinStateBuilder` remains the pre-resolution type for partial knowledge:

```rust
pub struct SpinStateBuilder {
    unpaired_electrons: Option<u8>,
    multiplicity: Option<SpinMultiplicity>,
}
```

Invariant: when both fields are `Some`, they must satisfy
`is_valid_spin_state(unpaired_electrons, multiplicity)`.

API:

```rust
impl SpinStateBuilder {
    pub fn unknown() -> Self;
    pub fn with_unpaired(n: u8) -> Self;
    pub fn with_multiplicity(m: SpinMultiplicity) -> Self;
    pub fn try_new(
        unpaired_electrons: Option<u8>,
        multiplicity: Option<SpinMultiplicity>,
    ) -> Result<Self, ResolutionError>;

    pub fn set_unpaired(&mut self, n: u8) -> Result<&mut Self, ResolutionError>;
    pub fn set_multiplicity(&mut self, m: SpinMultiplicity) -> Result<&mut Self, ResolutionError>;

    pub fn unpaired_electrons(&self) -> Option<u8>;
    pub fn multiplicity(&self) -> Option<SpinMultiplicity>;
    pub fn can_build(&self) -> bool;
    pub fn build(&self) -> Result<SpinState, ResolutionError>; // exact only, no guessing
}
```

No default Hund path in normal resolution. If max-multiplicity behavior is needed, it must
be opt-in and named at call sites.

### D1a. Builder API alignment (Atom/Bond/SpinState)

Builder method conventions are aligned across `AtomBuilder`, `BondBuilder`,
`SpinStateBuilder`:

1. `can_build(&self) -> bool`:
   - fast readiness predicate for call sites
   - equivalent to `self.build().is_ok()` semantically
2. `build(&self) -> Result<T, ResolutionError>`:
   - canonical checked conversion to resolved value
   - returns structured error on incompleteness or inconsistency
3. No `try_build` method:
   - redundant with `build`
   - avoids ambiguity around `Option` vs `Result` semantics

### D2. Resolved `Bond` stores full `SpinState`

Resolved bond shape:

```rust
pub struct Bond {
    order: u8,
    charge: i8,
    spin: SpinState,
}
```

Bond invariants (enforced in constructor/build path):

1. `electron_count = 2 * order - charge`
2. `electron_count >= 0`
3. `spin.unpaired_electrons() <= electron_count`
4. `(electron_count - spin.unpaired_electrons()) % 2 == 0` (remaining electrons can pair)
5. `SpinState` validity (`m <= n + 1`, parity) is already enforced by `SpinState`

`BondBuilder` changes from `multiplicity: Option<SpinMultiplicity>` to
`spin: SpinStateBuilder`.

### D3. `from_table_*` constructors are fallible

`AtomBuilder::from_table_atom` and `BondBuilder::from_table_bond` become fallible:

```rust
pub fn from_table_atom(atom: &TableAtom) -> Result<Self, ResolutionError>;
pub fn from_table_bond(bond: &TableBond) -> Result<Self, ResolutionError>;
```

They do consistency checks immediately and return `ResolutionError` for impossible input
combinations, while keeping TableIR as a faithful raw container.

### D4. Aromatic systems have explicit resolved spin

Resolved aromatic system carries:

```rust
pub struct AromaticSystem {
    contributions: Vec<AromaticContribution>,
    charge: i8,
    spin: SpinState,
    rings: Vec<Ring>,
}
```

No `Option<SpinState>` on the resolved type.

Builder-time aromatic spin is partial (`SpinStateBuilder`) and must resolve by build time.

### D5. Molecular spin: infer only if unique

Molecular multiplicity is never guessed as high-spin by default.

If molecular spin is explicit in input:

1. Verify unpaired sum and coupling compatibility
2. Error on mismatch

If molecular spin is not explicit:

1. Compute all compatible molecular multiplicities from feature-level `SpinState` values
2. If exactly one, infer it
3. If more than one, return ambiguity error (resolved `Molecule` cannot carry unknown spin)

### D6. TableIR atom spin fields are split

Replace `UnpairedElectrons` on TableIR atom structs with separate fields:

```rust
unpaired_electrons: Option<u8>,
multiplicity: Option<SpinMultiplicity>,
```

Applies to both basic and extended atom structs.

Parser semantics:

1. SMILES radical / MOL RAD set `unpaired_electrons`
2. CX `^` codes set:
   - unpaired only for codes `1,2,5`
   - unpaired plus multiplicity for codes `3,4,6,7`

`UnpairedElectrons` is removed from shared model types.

## Integration targets

### AtomBuilder

`AtomBuilder` moves from:

```rust
unpaired_electrons: Option<UnpairedElectrons>
```

to:

```rust
spin: SpinStateBuilder
```

### BondBuilder

`BondBuilder` moves from:

```rust
multiplicity: Option<SpinMultiplicity>
```

to:

```rust
spin: SpinStateBuilder
```

Resolved `Bond` now exposes:

```rust
pub fn spin(&self) -> SpinState;
pub fn multiplicity(&self) -> SpinMultiplicity; // convenience passthrough
pub fn unpaired_electrons(&self) -> u8;         // convenience passthrough
```

### Aromatic systems in `MoleculeBuilder`

`MoleculeBuilder` should carry aromatic-system builders (partial spin state) and produce
resolved `AromaticSystem` with definite `SpinState` at build.

### `Molecule::charge()` and `Molecule::spin()`

Implement `todo!()` methods in resolved molecule:

1. `charge()` = atom + bond + aromatic + multicenter contributions
2. `spin()` = resolved molecular spin guaranteed by build-time validation rules above

## Error model updates

Add explicit spin-related `ResolutionError` variants (names can change, semantics fixed):

1. Invalid partial spin state combination (`unpaired`, `multiplicity`)
2. Incomplete feature spin at build (`Atom`, `Bond`, `AromaticSystem`, `Multicenter`)
3. Bond electron/spin inconsistency (`order`, `charge`, `spin`)
4. Molecular spin ambiguity (multiple compatible multiplicities)
5. Molecular spin incompatibility with explicit annotation

## Implementation plan

### Phase 1. Introduce `SpinStateBuilder`

Files:

1. `umol-models-graph/src/graph_ir/spin_state_builder.rs` (new)
2. `umol-models-graph/src/graph_ir/mod.rs` (exports)
3. `umol-models-graph/src/graph_ir/error.rs` (new spin errors)

Tasks:

1. Implement type + API + invariants
2. Add unit tests for constructor/setter validation and `build()`

Acceptance:

1. Invalid `(n,m)` combos fail at mutation time
2. `can_build()` reflects readiness correctly
3. `build()` returns a structured incompleteness error when fields are missing

### Phase 2. TableIR atom spin-field migration

Files:

1. `umol-models-graph/src/table_ir/atom.rs`
2. `umol-models-graph/src/io/ctfile/parser/atom.rs`
3. `umol-models-graph/src/io/ctfile/parser/accumulator.rs`
4. `umol-models-graph/src/io/smiles/parser/cx.rs`
5. `umol-models-graph/src/io/smiles/parser/tests/cx.rs`
6. `umol-models-graph/src/io/ctfile/parser/accumulator/tests.rs`
7. `umol-models-graph/src/io/ctfile/parser/convert.rs` tests as needed

Tasks:

1. Replace `Option<UnpairedElectrons>` with separate `Option<u8>` and
   `Option<SpinMultiplicity>`
2. Update CX radical decoding to produce split values
3. Update CT/MOL paths to set unpaired count and leave multiplicity unset unless explicitly known

Acceptance:

1. Parser tests pass with unchanged semantics
2. CX code mapping remains exact (`1..7`)

### Phase 3. AtomBuilder migration and fallible TableIR conversion

Files:

1. `umol-models-graph/src/graph_ir/atom.rs`
2. `umol-models-graph/src/graph_ir/resolution.rs`
3. `umol-models-graph/src/graph_ir/symmetry.rs`
4. `umol-models-graph/src/graph_ir/valence.rs`

Tasks:

1. Replace atom builder spin representation with `SpinStateBuilder`
2. Change `from_table_atom` to `Result`
3. Update call sites in resolution pipeline
4. Add `can_build()` on `AtomBuilder` (matching builder convention)
5. Keep valence-facing getter surface stable (`unpaired_electrons()`, `multiplicity()`)

Acceptance:

1. No behavior regression in valence candidate generation except new early spin consistency errors
2. `resolve_builder` propagates `from_table_atom` errors via `?`
3. `AtomBuilder::can_build()` and `build()` are semantically consistent

### Phase 4. Bond refactor to full `SpinState`

Files:

1. `umol-models-graph/src/graph_ir/bond.rs`
2. `umol-models-graph/src/table_ir/bond.rs` (if helper methods added)
3. `umol-models-graph/src/graph_ir/resolution.rs`
4. Tests under `graph_ir` and parser integration

Tasks:

1. Update resolved `Bond` to store `SpinState`
2. Update `BondBuilder` to carry `SpinStateBuilder`
3. Make `from_table_bond` fallible
4. Make `BondBuilder::build()` fallible (`Result<Bond, ResolutionError>`)
5. Add `BondBuilder::can_build()` and enforce electron-count/spin invariants in build path

Acceptance:

1. Invalid table bond combinations fail deterministically
2. Resolved bond always has internally consistent `(order, charge, spin)`
3. `BondBuilder` follows the shared builder contract (`can_build` + fallible `build`)

### Phase 5. Aromatic system spin integration

Files:

1. `umol-models-graph/src/graph_ir/aromaticity.rs`
2. `umol-models-graph/src/graph_ir/molecule/builder.rs`
3. `umol-models-graph/src/graph_ir/molecule.rs`

Tasks:

1. Add resolved `spin: SpinState` to aromatic systems
2. Add builder-time aromatic spin representation
3. Ensure aromaticity perception/model layer sets spin explicitly
4. Add aromatic spin/electron consistency checks in build

Acceptance:

1. No resolved aromatic system without a definite spin
2. Existing aromatic tests are updated to set or derive spin explicitly

### Phase 6. Molecular spin and charge finalization

Files:

1. `umol-models-graph/src/graph_ir/molecule/builder.rs`
2. `umol-models-graph/src/graph_ir/molecule.rs`
3. `umol-models-graph/src/graph_ir/error.rs`

Tasks:

1. Implement `Molecule::charge()` (remove `todo!()`)
2. Implement molecular spin resolution rules:
   - explicit annotation: validate
   - no annotation: infer only when unique
   - ambiguity: error
3. Remove current implicit high-spin fallback

Acceptance:

1. `Molecule::spin()` is deterministic for resolved molecules
2. Ambiguous cases are rejected with dedicated diagnostics

### Phase 7. Remove `UnpairedElectrons` and finish migration

Files:

1. `umol-models-graph/src/atom.rs`
2. All remaining imports/usages found by `rg "UnpairedElectrons"`

Tasks:

1. Delete `UnpairedElectrons` type
2. Replace all remaining uses with split fields or `SpinStateBuilder`
3. Update docs/tests/snapshots

Acceptance:

1. No `UnpairedElectrons` references remain
2. Full workspace builds and tests pass

## Test plan

Add/adjust tests for:

1. `SpinStateBuilder` invariant enforcement
2. Fallible `from_table_atom` / `from_table_bond`
3. Bond electron-count/spin invariant failures
4. Aromatic-system spin presence and consistency
5. Molecular spin ambiguity rejection (no default high-spin)
6. CX radical code split-field mapping
7. `can_build()` <-> `build()` consistency on all three builder types

## Non-goals in this rollout

1. Full multicenter spin model (keep current behavior, design separately)
2. New external syntax for top-level molecular multiplicity annotations
3. Performance optimization of spin-coupling enumeration beyond correctness and clarity
