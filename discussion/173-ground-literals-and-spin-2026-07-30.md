# 173 — Ground literals and spin-state boundaries

Status: **In Progress**
Date: 2026-07-30
Relates: [061](061-spin-state-builder-2026-03-22.md),
[172](172-ast-literal-extraction-2026-07-30.md),
[175](175-ground-ast-api-2026-07-31.md)

## Scope

This document defines the contract needed before adding a checked `Ground<T>` wrapper. It separates
three operations that are currently mixed in some `AsLit` implementations:

1. exact projection of a structurally ground AST value;
2. operation-specific reduction of that value, such as its contribution to electron counting;
3. validation or conversion into a physically constrained chemistry type.

The general extraction policy remains in doc 172. This document covers the underlying literal
contract, the exceptional leaf types, and the relation among the existing spin types and validation
stages. The concrete `Ground<T>` API is tracked separately in doc 175; the staged plan below
covers this document's exact-literal and spin migration.

The implementation is a repository-wide semantic migration, not only a type rename. Every stage
must sweep the affected DSL specification, unit tests, property-test strategies and properties,
conformance suites, and fuzz targets and seed corpora. Fixtures, examples, benchmarks, snapshots,
and Python tests that construct or inspect the renamed values are included in the same sweep. The
later staged implementation plan must assign these consumers to the subitem that changes their
contract rather than collecting them as deferred cleanup.

## Terms

### Structurally ground

`Lattice::is_ground` describes the AST lattice. A value is ground when it is resolved to a bottom
element of that lattice. For a product such as `UnpairedElectronsAst`, this means that every
component is ground.

Groundness does not imply that the value satisfies chemistry invariants or that the surrounding
entity structure is internally consistent. Those are separate validation layers.

### Exact literal projection

`AsLit` should project a structurally ground AST value into a non-lattice representation without:

- canonicalizing the input;
- applying defaults;
- validating chemistry invariants;
- collapsing distinct ground AST states for the convenience of an operation.

For every type implementing both `Lattice` and `AsLit`, the required totality law is:

```rust
value.is_ground() == value.as_lit().is_some()
```

The projection should also be faithful on ground values: two distinct canonical ground values must
not become indistinguishable solely because they have the same downstream numerical effect.

These laws make an infallible `Ground<T>::lit()` sound. They do not make all operations on a ground
molecule infallible: topology-derived values may additionally depend on entity-integrity
preconditions.

### Operation-specific reduction

An operation may intentionally identify distinct ground values. This is not literal projection.

`AromaticValenceAst::NotAromatic` and `AromaticValenceAst::Aromatic(ValueAst::Lit(0))` have distinct
structural meanings but both supply zero to the relevant electron-counting calculations.
`MulticenterValenceAst::NotMulticenter` and
`MulticenterValenceAst::Multicenter(ValueAst::Lit(0))` have the analogous relationship.

These types therefore need both:

- an exact `AsLit` projection that preserves the structural distinction;
- a separately named total method on the exact carrier with the calculation semantics: explicit
  absence maps to zero and a present value maps to that value.

The calculation method is not a second general extraction trait. It belongs on each domain type
because its meaning is specific to that value. Both methods are named `valence_count()`. The
AST-level calculation first performs exact projection with `as_lit()` and then maps the carrier's
total method, so a non-literal state remains `None`. The
existing `aromatic_increment` calculation is separately renamed `aromatic_covalence`, reflecting
the Langmuir covalence supplied by aromatic bonding rather than a generic increment.

## Leaf-type corrections

Most existing leaf implementations already make `is_ground` and `as_lit` agree. The following need
correction:

| AST type | Current mismatch | Required separation |
| --- | --- | --- |
| `IsotopeMassAst` | `Natural` is ground but projects to `None` | Exact carrier distinguishes natural composition from an exact mass number |
| `TetrahedralStereoAst` | `NotStereo` is ground but projects to `None` | Exact carrier distinguishes absence from a literal tetrahedral configuration |
| `CisTransStereoAst` | `NotStereo` is ground but projects to `None` | Exact carrier distinguishes absence from a literal cis/trans configuration |
| `AromaticValenceAst` | Ground forms project, but absence and present zero collapse | Exact carrier plus a separate numerical calculation method |
| `MulticenterValenceAst` | Ground forms project, but absence and present zero collapse | Exact carrier plus a separate numerical calculation method |

The approved exact non-lattice carriers are:

```rust
pub enum IsotopeMass {
    Natural,
    MassNumber(u32),
}

pub enum AromaticValence {
    NotAromatic,
    Aromatic(i64),
}

pub enum MulticenterValence {
    NotMulticenter,
    Multicenter(i64),
}

pub enum TetrahedralStereo {
    NotStereo,
    Stereo(u32),
}

pub enum CisTransStereo {
    NotStereo,
    Stereo(u32),
}
```

Their `AsLit` projections preserve the distinction between the two ground variants. The exact
valence carriers expose total `valence_count() -> i64` methods; AST callers obtain the optional
calculation result with `as_lit().map(...)`. The per-kind stereo carrier determines the
kind, so its `Stereo` variant stores the structurally literal coset index directly and retains
out-of-range values for later validation.

The current Python configuration shorthands are renamed and their semantic definitions move to
`umol-ast`:

```rust
pub enum TetrahedralConfiguration {
    Ccw,
    Cw,
}

pub enum CisTransConfiguration {
    Z,
    E,
}
```

`umol-py` retains thin PyO3 wrappers with these names. This leaves `TetrahedralStereo` and
`CisTransStereo` available for the exact Rust carriers rather than using the same names for the
narrower configuration shorthands.

Expose all five exact carriers in Python under the same names:

- `IsotopeMass`;
- `AromaticValence`;
- `MulticenterValence`;
- `TetrahedralStereo`;
- `CisTransStereo`.

The corresponding AST `as_lit()` methods return the exact carrier or `None`, where `None` means
only that the AST value is not structurally ground. Python primitives cannot preserve these
distinctions uniformly: `None` is already the non-ground result, and Python equality identifies
`False` with integer zero. `AromaticValence` and `MulticenterValence` separately expose total
`valence_count() -> int` methods for the calculation-specific projection.

## Unpaired-electron and spin types

### `SpinMultiplicity`

Replace the fixed `Singlet` through `Decet` enum with a validated numeric newtype:

```rust
pub struct SpinMultiplicity(NonZeroU8);
```

It represents multiplicities `1..=255`, removing the artificial decet ceiling while keeping a
one-byte value. Its checked constructor follows the underlying nonzero integer:

```rust
pub const fn new(value: u8) -> Option<Self>;
```

It returns `None` only for zero. The conventional singlet through decet names remain available as
uppercase associated constants and explicit optional lookup:

```rust
pub const SINGLET: Self;
pub const DOUBLET: Self;
pub const TRIPLET: Self;
pub const QUARTET: Self;
pub const QUINTET: Self;
pub const SEXTET: Self;
pub const SEPTET: Self;
pub const OCTET: Self;
pub const NONET: Self;
pub const DECET: Self;

pub const fn name(self) -> Option<&'static str>;
pub fn from_name(name: &str) -> Option<Self>;
```

The forward lookup is ASCII-case-insensitive, and multiplicities above ten have no conventional
name. Canonical display and serialization use the numeric multiplicity. Remove the current
word-based `FromStr`: numeric construction belongs to `new`, while the optional chemical
nomenclature remains explicit rather than defining the type's domain.

An `u8` unpaired-electron count admits the theoretical maximum-spin boundary
`#u255#s256`, which this representation cannot express. This single extreme boundary does not
justify a wider multiplicity representation for the intended chemical domain. Compatibility
arithmetic must nevertheless widen before evaluating `unpaired_electrons + 1`.

### `UnpairedElectrons`

`UnpairedElectrons` is the exact, non-lattice representation of the two values recorded
structurally:

```rust
pub struct UnpairedElectrons {
    pub count: i64,
    pub multiplicity: i64,
}
```

The integer widths preserve every literal accepted by `ValueAst`, including negative and
out-of-range values. Constructing this type asserts only structural completeness, not physical
validity. It belongs in `umol-chem::spin` next to the validated conversion target rather than in
`umol-ast` as a one-off projection type.

### `SpinState`

`SpinState` in `umol-chem` is a validated physical pair:

- the unpaired-electron count is a bounded `u8`;
- multiplicity is a `SpinMultiplicity`;
- multiplicity cannot exceed `unpaired_electrons + 1`;
- the parity of multiplicity and `unpaired_electrons + 1` must agree.

Its constructor and parser enforce these invariants. Its parser may also derive a missing component:
`#u` alone chooses maximum multiplicity, while `#s` alone chooses the minimum compatible unpaired
count.

The fields and accessors remain descriptive at this level:

```rust
pub struct SpinState {
    unpaired_electrons: u8,
    multiplicity: SpinMultiplicity,
}
```

`SpinState` is produced by a checked conversion from `UnpairedElectrons`. The conversion must use
checked integer conversions before applying the compatibility rule; `as u8` would incorrectly wrap
negative or oversized literals.

The approved construction and conversion API is:

```rust
impl SpinState {
    pub fn new(
        unpaired_electrons: u8,
        multiplicity: SpinMultiplicity,
    ) -> Result<Self, SpinStateError>;
}

impl TryFrom<UnpairedElectrons> for SpinState {
    type Error = SpinStateError;
}

impl From<SpinState> for UnpairedElectrons;
```

The physical error variants are:

```rust
SpinStateError::UnpairedElectronsOutOfRange { count: i64 }
SpinStateError::MultiplicityOutOfRange { multiplicity: i64 }
SpinStateError::Incompatible {
    unpaired_electrons: u8,
    multiplicity: SpinMultiplicity,
}
```

`TryFrom<UnpairedElectrons>` checks both integer ranges and delegates compatibility validation to
`SpinState::new`. The existing syntax-error variants remain available to `SpinState::from_str`.

### `UnpairedElectronsAst`

`UnpairedElectronsAst` in `umol-ast` is a product of two independent `ValueAst` fields:

```rust
pub struct UnpairedElectronsAst {
    pub count: ValueAst,
    pub multiplicity: ValueAst,
}
```

It serves several purposes that `SpinState` cannot:

- partial knowledge of either component;
- value patterns and expressions;
- componentwise lattice operations;
- componentwise updates, including removal to `Undetermined`;
- delayed application of DSL defaults.

Its canonicalization and lattice operations are intentionally structural. A literal pair that fails
spin parity is retained and is still structurally ground; physics belongs to a later validator.

`UnpairedElectronsUpdate` is the componentwise update type:

```rust
pub struct UnpairedElectronsUpdate {
    pub count: Option<ValueAst>,
    pub multiplicity: Option<ValueAst>,
}
```

Atoms, localized bonds, aromatic systems, and multicenter bonds store this value in an
`unpaired_electrons` field. Their views expose `unpaired_electrons()`, and their builders, updates,
edits, deltas, dictionaries, and Python properties use the same terminology. The short field names
remain unambiguous under the containing value:

```rust
atom.unpaired_electrons.count
atom.unpaired_electrons.multiplicity
```

`AsLit` projects a structurally ground `UnpairedElectronsAst` to `UnpairedElectrons`. Physical
conversion to `SpinState` remains a separate fallible operation:

```text
UnpairedElectronsAst
    -> Option<UnpairedElectrons>
    -> Result<SpinState, SpinStateError>
```

### Unpaired-electron coupling constraint

The molecule-level `SpinSum` name incorrectly suggests arithmetic addition even though spin
multiplicities are coupled. Replace it with:

```rust
MoleculeConstraint::UnpairedElectronCoupling {
    atoms: Option<Vec<AtomId>>,
    unpaired_electrons: UnpairedElectronsAst,
}
```

The EDN form uses the same component names as `UnpairedElectronsAst`:

```edn
{:unpaired-electron-coupling
 {:atoms [0]
  :unpaired-electrons {:count 2 :multiplicity 3}}}
```

Omitting `:atoms` retains the existing all-atoms scope. A constraint with both component values
`Undetermined` remains vacuous. The compact entity syntax continues to use `#u` and `#s`.

The DSL grammar names the reusable value directly:

```text
unpaired-electrons ::=
    { :count value-expr :multiplicity value-expr }

molecule-constraint ::=
    { :unpaired-electron-coupling
      { [:atoms [atom-ref+]]?
        :unpaired-electrons unpaired-electrons } }
```

### Python exposure

Python currently exposes `SpinStateAst` and entity properties named `spin`. Rename these to
`UnpairedElectronsAst` and `unpaired_electrons`, including constructor arguments and dictionary
keys.

Expose `UnpairedElectrons` as an immutable Python value type with integer `count` and
`multiplicity` properties. It is the exact return type of `UnpairedElectronsAst.as_lit()` and does
not impose physical compatibility.

Expose `SpinState` as the physically validated Python value type. Its constructor accepts
keyword-only integer `unpaired_electrons` and `multiplicity` arguments, rejects invalid combinations
with `ValueError`, and provides read-only properties with the same names. Do not expose
`SpinMultiplicity` as a separate Python class: it is Rust representation and validation machinery,
while an ordinary integer is the natural Python representation. Construction of `SpinState` is the
explicit physical conversion; no additional conversion helpers are needed.

## Current validation timing

The current code applies spin-related rules at several boundaries:

1. Entity-string and EDN parsing currently build `SpinStateAst`, planned as
   `UnpairedElectronsAst`, structurally and accept partial or physics-invalid literal combinations.
2. DSL raising may derive a missing component according to `UnpairedElectronsDefault` and
   `MultiplicityDefault`; lowering removes derivable defaults for faithful rendering.
3. Valence and invariant enumeration filter candidate multiplicities using `SpinMultiplicity` and
   `SpinState` compatibility rules.
4. `SpinStateAst::as_lit` currently performs validated conversion to `SpinState`. After the rename,
   `UnpairedElectronsAst::as_lit` must instead perform exact structural projection without physical
   validation.
5. `SpinInvariantsValidator` is intended to validate literal spin pairs, but is currently a stub
   that always returns `Determined`. Its module comment names atoms, aromatic systems, and
   multicenter bonds, omitting localized bonds and molecule-level `SpinSum` constraints. The latter
   are planned as `UnpairedElectronCoupling` constraints carrying `UnpairedElectronsAst`.

The harmonized API must preserve the useful differences among these stages instead of validating
everything at parse or lattice-operation time.

### Spin-invariant validation

`SpinInvariantsValidator` validates the `UnpairedElectronsAst` values stored on atoms, localized
bonds, aromatic systems, and multicenter bonds. A value with either non-literal component
contributes `Underdetermined`. A complete literal pair is converted with
`SpinState::try_from`; conversion failure is a contradiction. An underdetermined value does not
short-circuit traversal because a later entity may still be contradictory.

The contradiction enum identifies each entity kind explicitly:

```rust
pub enum SpinInvariantsContradiction {
    MoleculeAtom {
        atom: AtomId,
        error: SpinStateError,
    },
    Bond {
        bond: BondId,
        error: SpinStateError,
    },
    AromaticSystem {
        system: AromaticSystemId,
        error: SpinStateError,
    },
    MulticenterBond {
        bond: MulticenterBondId,
        error: SpinStateError,
    },
    Atom {
        error: SpinStateError,
    },
    UnpairedElectronCoupling {
        constraint_index: usize,
        error: SpinStateError,
    },
}
```

`Atom` is used only by the context-free `validate_atom(&AtomAst)` entry point. `MoleculeAtom`
identifies an atom within a molecule. This avoids inventing an `AtomId` for an atom that does not
belong to a molecule.

The validator recursively inspects `UnpairedElectronCoupling` leaves inside each top-level
constraint tree. The top-level index identifies the containing tree, consistently with the
existing constraint validators. An invalid literal coupling target is contradictory. Until the
angular-momentum operation exists, every non-vacuous coupling constraint remains underdetermined
after this local target check because its actual satisfiability has not been evaluated.

`SpinInvariantsError` remains empty: these checks have no setup or operational failure mode.

## `SpinState` API cleanup

`SpinState` should validate one unpaired-electron count and multiplicity. It should not also own
angular-momentum coupling operations.

The local API cleanup is:

- make `SpinState::new` the fallible constructor and remove `try_new`;
- remove the panicking constructor path;
- remove public `are_compatible`, which duplicates the constructor's validation semantics;
- remove public `max_multiplicity`; the parser can derive its explicitly requested high-spin
  default internally;
- remove unused `is_compatible_with`;
- rename the `unpaired` field, accessor, arguments, locals, and error fields consistently;
- parse numerical tags without accumulating into `u8`, then perform checked conversion;
- widen error payloads so invalid `i64` AST literals and oversized parsed values are preserved;
- make `SpinMultiplicity` serialization agree with its canonical textual representation;
- correct the stale `#m` rustdoc examples to the actual `#s` syntax;
- remove unused `HIGHEST_SPIN_MULTIPLICITY` and the unused `FromRepr` conversion path.

The current arithmetic must not evaluate `unpaired_electrons + 1` until the supported range has
been established.

## Angular-momentum coupling boundary

Remove both coupling methods from `SpinState`:

```rust
SpinState::high_spin_combine
SpinState::is_constructible_from
```

Their current implementations are not a suitable basis for a general coupling API.
`high_spin_combine` derives multiplicity from the sum of unpaired-electron counts and can therefore
produce a result that `is_constructible_from` rejects when a component is itself in a lower-spin
state.

Angular-momentum coupling belongs in a separate, generic crate. The proposed dependency direction
is:

```text
umol-chem
  SpinMultiplicity, UnpairedElectrons, SpinState
  local physical validity only

umol-angular-momentum
  exact angular-momentum values
  admissible coupling calculations and validation

umol-graph
  SpinState/angular-momentum conversion
  molecular coupling operations
  UnpairedElectronCoupling constraint validation
```

`umol-chem` remains independent of the coupling crate. `umol-graph` depends on both and owns the
chemical interpretation. The generic crate may also support later `umol-msym` work.

For the current molecular operation, the generic crate needs exact doubled-integer angular momentum
(`2j`), pair and multiple-resultant enumeration, admissibility checks, maximum-coupling selection,
and checked arithmetic. Clebsch-Gordan coefficients, Wigner symbols, explicit coupled-basis trees,
and recoupling coefficients are separate extensions rather than prerequisites for validating final
resultants.

Existing implementations are useful references but not suitable direct dependencies for this
narrow core:

- [`wigner-symbols`](https://docs.rs/crate/wigner-symbols/latest) supports half-integer angular
  momenta and exact coefficients using doubled integers, but brings the substantially larger
  coefficient implementation and a `rug` dependency;
- [`wigners`](https://docs.rs/crate/wigners/latest) is pure Rust but currently supports only integer
  angular momenta;
- [SymPy coupled spin](https://docs.sympy.org/latest/modules/physics/quantum/spin.html) provides a
  suitable development-time reference for captured validation cases.

The generic admissibility crate is a small-to-moderate implementation unit; integration with
`umol-graph` and the molecule constraint is a separate unit. Exact coefficient calculation would be
materially larger.

The immediately actionable result in this document is only removal of the two misplaced
`SpinState` coupling methods. The generic crate and higher-level graph operation require their own
API design and implementation plan.

## Settled spin design

The following decisions are settled:

- add `UnpairedElectrons` as the exact non-lattice carrier;
- rename `SpinStateAst` to `UnpairedElectronsAst`;
- rename `SpinStateUpdate` to `UnpairedElectronsUpdate`;
- use `count` and `multiplicity` inside the `UnpairedElectrons` family;
- use `unpaired_electrons` on surrounding entities and on `SpinState`;
- keep `UnpairedElectronsAst` structurally permissive and componentwise;
- keep `SpinState` physically validated;
- keep `SpinState` limited to local validity and remove its two coupling methods;
- use one fallible `SpinState::new` construction path and remove public `are_compatible`;
- keep compact `#u` and `#s` syntax;
- replace `SpinSum` with `UnpairedElectronCoupling` and use the EDN form above;
- do not redefine lattice groundness as physical validity;
- make physical conversion explicit rather than hiding it in exact literal projection;
- expose `UnpairedElectrons` and `SpinState` in Python, but keep `SpinMultiplicity` internal to the
  Rust boundary;
- implement the tier-2 validator with entity-specific contradictions rather than relying on
  incidental validation in consumers.

## Ground-wrapper boundary

The intended wrapper remains:

```rust
pub struct Ground<T>(T);
```

The wrapper owns the value passed to it, not necessarily the underlying AST. `T` may itself be a
borrowed handle:

```rust
Ground<&AtomAst>
Ground<AtomView<'a>>
```

Owning the handle lets `Ground` wrap view values directly, including temporary views returned by
lookup operations, while the view continues to carry the borrow of its molecule. A
`Ground<'a, T>(&'a T)` wrapper would instead add another reference layer and could not own such a
view. The private field and checked construction establish groundness; ownership of `T` does not
imply ownership or cloning of the underlying molecular structure.

`Ground` is evidence of structural groundness only. For an exact leaf projection:

```rust
impl<T: AsLit> Ground<T> {
    pub fn lit(&self) -> T::Lit;
}
```

Aggregate ground views may return ground wrappers around stored fields. They may expose a concrete
derived value only where the derivation is proven to preserve groundness. Entity integrity and
chemistry validity must not be silently added to the meaning of `Ground`.

Define the wrapper in a new top-level `umol-ast/src/ast/ground.rs` module and re-export it from
`umol_ast::ast`. Its approved construction and access surface is:

```rust
Ground::new(value) -> Option<Self>
ground.as_ref() -> &T
ground.into_inner() -> T
```

`new` is implemented for the supported concrete wrapped forms and checks structural groundness.
`Ground<T>` implements `AsRef<T>` and does not implement `Deref`, keeping the checked surface
visible rather than silently falling back to the ordinary AST API. No public `Groundable` trait is
introduced.

The concrete useful API is broader than this generic wrapper contract: it must define groundness
recursively across every AST family, ground-preserving molecule/entity navigation, optional
constraint projections, and the boundary between stored and topology-derived values. That design
and its implementation are tracked in [doc 175](175-ground-ast-api-2026-07-31.md). This document
does not specify the supported concrete `Ground<T>` forms or their entity accessors.

## Staged implementation plan

Every subitem includes its directly affected tests under the test-writing conventions. Property,
conformance, fuzz, fixture, benchmark, snapshot, specification, and Python changes are named
explicitly where their contracts change; the final verification stage does not defer those
migrations.

### S0 — Establish the physical spin types

- **S0a — `umol-chem/src/spin.rs`: add `UnpairedElectrons`.** **Done.** Add the exact `i64` pair
  carrier with the approved field names and ordinary value-type trait implementations. Add unit
  tests for construction, equality, ordering, and serialization. **Additive (green).** [dep: none]
- **S0b — `umol-chem`: replace `SpinMultiplicity` and clean `SpinState`.** **Done.** Replace the
  ten-variant enum with the `NonZeroU8` newtype, associated conventional constants and lookups,
  numeric display/serialization, and checked `new`. Replace the `SpinState` construction and
  conversion surface, widen parse arithmetic and error payloads, rename `unpaired` throughout,
  and remove `are_compatible`, `max_multiplicity`, `is_compatible_with`, `high_spin_combine`, and
  `is_constructible_from`. Remove the decet-derived maximum-unpaired-electron assertion and the
  unused `strum` dependency. Test the full numeric range, conventional-name roundtrips, parsing,
  checked conversions, compatibility boundaries, serialization, and the `UnpairedElectrons`
  conversion roundtrip. **Breaking (red→green within S0).** [dep: S0a]
- **S0c — `umol-io`, `umol-geometric`, `umol-geometric-graph`: migrate physical-spin
  consumers.** **Done.** Update TableIR, CXSMILES and CTfile radical handling, geometric molecule
  conversion, and their fixtures to use the numeric multiplicity API and approved field names.
  Preserve external-format behavior in the unit and parsing-conformance suites; sweep the
  OpenSMILES fuzz target and seeds that exercise radicals. **Breaking caller migration
  (red→green within S0).** [dep: S0b]
- **S0d — `umol-ast`, `umol-graph`: migrate direct `SpinState` consumers.** **Done.** Update current
  construction, comparison, resolution, validation, and test code that consumes the physical
  type. Remove tests for the deleted coupling methods rather than replacing them with a local
  approximation. Run the affected graph conformance tests. **Breaking caller migration
  (red→green within S0).** [dep: S0b]

S0 ends with the workspace green and no fixed-decet or local angular-momentum-coupling API
remaining.

### S1 — Correct exact projection for the exceptional leaves

- **S1a — `umol-ast/src/ast/{atom,constraint/atom,stereo}.rs`: add exact carriers and
  configuration enums.** **Done.** Add `IsotopeMass`, `AromaticValence`, `MulticenterValence`,
  `TetrahedralStereo`, and `CisTransStereo`, plus `TetrahedralConfiguration` and
  `CisTransConfiguration`. Test every variant and its canonical value semantics.
  **Additive (green).** [dep: none]
- **S1b — the same AST modules and `traits.rs`: make `AsLit` exact.** **Done.** Change the five `AsLit`
  implementations to return their exact carriers, add total `valence_count()` methods to the two
  exact valence carriers, and update Rustdoc so `AsLit` states the `is_ground` totality and
  faithfulness contract.
  Update all Rust callers and unit tests in the same subitem. Add property tests over all affected
  variants for `value.is_ground() == value.as_lit().is_some()` and preservation of distinct ground
  forms; update the relevant property strategies and inspect the constraint/entity-string fuzz
  targets and seeds for both absent and zero-valued forms. **Breaking (red→green).** [dep: S1a]
- **S1c — `umol-ast/src/ast/view/atom.rs` and `umol-graph`: clarify aromatic calculations.**
  **Done.** Rename `aromatic_increment` to `aromatic_covalence`, retain the calculation-specific
  collapsing behavior through `valence_count()`, and migrate electron-counting, aromaticity, and
  valence callers. Add focused unit and property cases distinguishing structural absence from
  present zero while preserving equal numerical contribution. **Breaking (red→green).** [dep: S1b]
- **S1d — `umol-py`: expose the exact carriers.**
  **Done.** Rename the existing `TetrahedralStereo` and `CisTransStereo` configuration shorthands
  to `TetrahedralConfiguration` and
  `CisTransConfiguration`, bind all five exact carriers, and make each Python `as_lit()` return
  the carrier or `None`. Expose both carrier `valence_count()` methods. Update exports, import
  tests, unit tests, and Python fixtures without using Python primitive sentinels.
  **Breaking (red→green).**
  [dep: S1a, S1b]

S1 ends with exact literal projection having the same structural meaning in Rust and Python.

### S2 — Migrate the AST unpaired-electron vocabulary

- **S2a — `umol-ast/src/ast/spin.rs`, `ast.rs`: rename the AST and update types.** **Done.** Replace
  `SpinStateAst`/`SpinStateUpdate` with `UnpairedElectronsAst`/`UnpairedElectronsUpdate`, rename
  the component field to `count`, keep componentwise update and lattice behavior, and change
  `AsLit` to exact projection into `umol_chem::spin::UnpairedElectrons`. Add unit and property
  tests for partial values, groundness, canonicalization, lattice operations,
  `difference_to`/`update`, and exact projection. **Breaking (red→green within S2).** [dep: S0a,
  S0b]
- **S2b — `umol-ast/src/ast`: migrate entity storage and operations.** **Done.** Rename the `spin` field to
  `unpaired_electrons` on atoms, bonds, aromatic systems, and multicenter bonds; migrate builders,
  views, coloring, edits, deltas, transactions, reactions, dictionaries, and all corresponding
  field-change variants and tests. Extend property coverage for update/difference and
  edit/undo/delta roundtrips under partial component updates. **Breaking caller migration
  (red→green within S2).** [dep: S2a]
- **S2c — `umol-ast/src/dsl`, `umol-ast/spec/umol-dsl-spec.md`: migrate parsing and
  rendering.** **Done.** Rename Rust DSL fields and internal types while retaining `#u`/`#s`; update
  defaults, raising/lowering, EDN-shaped forms, parsing benchmarks, fixtures, and serialization
  expectations. Update parser roundtrip properties and sweep `fuzz_entity_strings`,
  `fuzz_molecule`, `fuzz_reaction`, and their spin-bearing seeds. **Breaking caller migration
  (red→green within S2).** [dep: S2a, S2b]
- **S2d — `umol-graph`, `umol-io`, `umol-geometric*`: migrate AST consumers.** **Done.** Update
  resolution, valence, aromaticity, transformation, validation, TableIR raising, and geometric
  conversion to the renamed AST and entity fields. Bring unit and conformance suites, examples,
  fixtures, snapshots, and benchmarks into conformance. **Breaking caller migration
  (red→green within S2).** [dep: S2b, S2c]
- **S2e — `umol-py`: migrate the AST bindings.** Replace `SpinStateAst` with
  `UnpairedElectronsAst`, rename entity constructor keywords, properties, update/delta variants,
  dictionary keys, reprs, and exports to `unpaired_electrons`, `count`, and `multiplicity`.
  Update the complete affected Python fixture and test surface. **Breaking caller migration
  (red→green within S2).** [dep: S2a, S2b]

S2 is intentionally one red→green stage because the public type and field rename crosses the
workspace. It ends with no compatibility aliases or stale code/specification names.

### S3 — Rename and reshape the coupling constraint

- **S3a — `umol-ast/src/ast/constraint/molecule.rs`: replace `SpinSum`.** Add
  `MoleculeConstraint::UnpairedElectronCoupling` with `atoms` and `unpaired_electrons`, and migrate
  vacuity, canonicalization, compaction, remapping, equality, reaction integrity, and constraint
  traversal. Update unit and property tests, including all-atoms and explicit-subset forms.
  **Breaking (red→green within S3).** [dep: S2a, S2b]
- **S3b — `umol-ast/src/dsl/constraint.rs`, DSL specification: replace the serialized form.**
  Implement the approved `:unpaired-electron-coupling` and nested
  `:unpaired-electrons {:count ... :multiplicity ...}` grammar, raising, lowering, rendering, and
  parse/render roundtrips. Replace the `mol_spin_sum` fuzz seed and add partial and complete
  coupling seeds; sweep molecule and reaction fuzz corpora for the old form. **Breaking caller
  migration (red→green within S3).** [dep: S3a]
- **S3c — downstream Rust consumers:** migrate geometric conversion and every remaining AST,
  graph, IO, fixture, example, benchmark, property, and conformance reference to the new
  constraint. Preserve the existing all-atoms scope semantics. **Breaking caller migration
  (red→green within S3).** [dep: S3a, S3b]
- **S3d — `umol-py`: migrate the constraint binding.** Rename the Python variant, constructor
  keywords, dictionary form, repr, and tests to `UnpairedElectronCoupling` and
  `unpaired_electrons`. **Breaking caller migration (red→green within S3).** [dep: S3a, S3b]

S3 ends with the old arithmetic `SpinSum` terminology absent from active code, specifications, and
test/fuzz inputs.

### S4 — Implement spin-invariant validation

- **S4a — `umol-graph/src/ops/validate/spin.rs`: validate an `AtomAst`.** Add the shared exact
  pair check and implement `validate_atom`, returning `Underdetermined` for a partial pair,
  `Determined` for a valid literal pair, and `SpinInvariantsContradiction::Atom` for a physically
  invalid literal pair. Test range, parity, and partial cases. **Additive behavior (green).**
  [dep: S0b, S2a]
- **S4b — the same module: validate molecule entities.** Traverse atoms, bonds, aromatic systems,
  and multicenter bonds, reporting `MoleculeAtom`, `Bond`, `AromaticSystem`, and
  `MulticenterBond` with their concrete IDs. Accumulate underdetermination without allowing it to
  mask a later contradiction. Add table tests for every entity kind and property tests over
  generated mixed entity states. **Additive behavior (green).** [dep: S4a]
- **S4c — the same module and constraint traversal: inspect coupling constraints.** Recursively
  inspect `UnpairedElectronCoupling` leaves in `And`, `Or`, and `Not` trees, report the containing
  top-level `constraint_index` for invalid targets, and leave every non-vacuous coupling
  underdetermined until the angular-momentum operation exists. Test nested trees, vacuous
  constraints, invalid targets, and valid-but-not-yet-evaluated targets. **Additive behavior
  (green).** [dep: S3a, S4a]
- **S4d — composite validation and consumers:** verify `Validator::validate_atom`,
  `validate_invariants`, and transformation callers propagate the new contradictions and
  underdetermination exactly. Add integration/property cases and run the resolution,
  aromaticity, and kekulization conformance suites affected by the formerly inert validator.
  **Additive behavior (green).** [dep: S4b, S4c]

### S5 — Complete the Python spin value surface

- **S5a — `umol-py/src/spin.rs`: bind `UnpairedElectrons`.** Add the immutable exact value type
  with read-only integer `count` and `multiplicity`, equality/hash/repr consistent with other
  Python value wrappers, Rust conversion, exports, and import/unit tests. **Additive (green).**
  [dep: S0a]
- **S5b — `umol-py/src/spin.rs`: bind `SpinState`.** Add the immutable validated value type with
  keyword-only integer `unpaired_electrons` and `multiplicity`, read-only properties, and
  `ValueError` conversion for range or compatibility failures. Do not expose
  `SpinMultiplicity`. Test valid construction, each physical error class, equality/hash/repr, and
  Rust roundtrips. **Additive (green).** [dep: S0b]
- **S5c — `umol-py`: connect exact projection.** Expose
  `UnpairedElectronsAst.as_lit() -> Optional[UnpairedElectrons]`, update annotations, exports,
  fixtures, and tests so partial and physics-invalid ground pairs remain distinguishable.
  **Additive (green).** [dep: S2a, S2e, S5a]

### S6 — Verify and close the migration

- **S6a — active specification and source audit:** search active Rust, Python, DSL specification,
  test, benchmark, fixture, snapshot, property strategy, fuzz target, and named seed files for
  stale `SpinStateAst`, `SpinStateUpdate`, entity `spin`, `SpinSum`, old stereo shorthand names,
  `aromatic_increment`, and fixed-enum multiplicity assumptions. Correct findings in their owning
  earlier subitem rather than treating S6 as a cleanup bucket. **Verification (green).** [dep:
  S0d, S1d, S2e, S3d, S4d, S5c]
- **S6b — full validation gate:** run formatting and workspace clippy/tests; the extended
  `umol-ast` and relevant `umol-graph-core` property suites; IO and graph conformance suites; the
  Python 3.13 binding build and pytest suite; fuzz-target builds and replay of all affected named
  seeds; and benchmark compilation for changed benches. Update doc 173 and status tracking only
  after all gates pass. **Verification (green).** [dep: S6a]

The critical path is **S0a → S0b → S1a/S2a → S2b → S2c → S2d/S2e → S3a → S3b →
S3c/S3d → S4c → S4d → S6**. S1 can otherwise proceed after the spin foundation, and S5 can
proceed once its Rust types are stable. No stage in this plan is deferrable: exact projection, the
full vocabulary migration, validation, Python parity, and every specified
test/specification/fuzz sweep are part of the deliverable. The concrete ground AST API is the
separate doc 175 work unit; the angular-momentum crate also remains outside this plan.

## Follow-up work

Design the generic angular-momentum crate and the higher-level `umol-graph` coupling operation in a
separate discussion document. This does not block the exact-carrier, spin-type, validator, or
generic `Ground<T>` contract defined here.
