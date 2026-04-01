# Error Handling Architecture

## Problem

The current error types in umol are a mix of catch-all enums (`AtomError` with 18+ variants spanning parse/validation/query concerns), `#[from]`/`#[transparent]` wrapping creating false IS-A hierarchies, and `.map_err(|e| SomeError::Variant(e.to_string()))` that discards structured information.

These are not domain-crossing boundaries — parsing, lowering, validation, resolution, and transformation are distinct concerns within the same domain. They need different error types, but those types must compose without nesting or flattening.

The top-level `umol::error::Error` is a god-enum with 9 sub-enums, all `#[from]`/`#[transparent]`. `DataError` alone has 30+ variants. Most variants are `SomeThing(String)`.

## Survey of Existing Approaches

| Project | Public error shape | Cross-concern composition |
|---|---|---|
| Diesel | Per-concern enums | Explicit variant wrapping + `Box<dyn DatabaseErrorInformation>` for DB details |
| SQLx | Single flat enum (~20 variants) | `Box<dyn DatabaseError>` + `BoxDynError` for inner causes |
| Hyper | Opaque struct, private `Kind` enum | Inspector methods, no matching. `Box<dyn StdError>` cause |
| Axum | Newtype over `Box<dyn StdError>` | Full erasure at boundary |

Common thread: `Box<dyn StdError + Send + Sync>` at concern boundaries. Variation is how much structure surrounds the box.

## Proposed Design: Three-Tier Error Architecture

### Tier 1: Sub-concern enums (direct `?`)

Focused enums for each concern within a module. 5-15 variants each. Flat — no `#[from]` wrapping of other error types. Cross-cutting errors from external crates (e.g. `SpinStateError` from `umol-data`) are mapped via explicit `From` impls that destructure into native variants.

Internal functions return these directly.

### Tier 2: Module dispatch enums (`#[from]` wrapping)

One dispatch enum per top-level crate module. Maps to first-level modules in the crate. Each dispatch enum derives `thiserror::Error` and implements `UmolError`.

### Tier 3: Cross-module boundary (`Box<dyn UmolError>`)

`UmolError` trait in the umol base crate. Functions that combine concerns across top-level modules return `Box<dyn UmolError>`. Type erasure at this level. Callers can downcast if needed.

### Flow

```
ValidationError ──┐
                   ├──#[from]──► GraphIrError ──┐
ResolutionError ──┘                              ├──box──► Box<dyn UmolError>
                                                 │
DslError ────────────────────────────────────────┘
```

## Rules

- **Public surface = module error type.** `Atom::from_str` returns `GraphIrError`, not `ValidationError`.
- **Internal functions use narrow types.** `check_invariants()` returns `ValidationError`; `?` promotes via `#[from]` on `GraphIrError`.
- **No `#[from]` across tiers.** `GraphIrError` does not impl `From<DslError>`. Cross-module = box.
- **No string flattening.** `SomeError::Variant(e.to_string())` is not acceptable — it discards structured information.
- **No wrapped error variants.** Sub-concern enums do not wrap other error enums. Use explicit `From` impls that destructure source errors into native variants.
- **thiserror for both tiers.** snafu context chains solve a different problem (annotations, not organization) and can be layered on later.

## Cleanup: `umol::error::Error` and `umol-models`

- **Remove `umol-models` crate entirely.** It was exploratory.
- **Remove `umol::error::Error` god-enum** and all its sub-enums (`ModelError`, `PropertyError`, `ParseError`, `ConversionError`, `OperationError`, `EntityError`, `ValidationError`, `SerializationError`, `DataError`).
- **Add flat `DataError` in `umol-data`**, derives `thiserror::Error` and `UmolError`. No sub-concern split — not worth the trouble for a data crate. Replaces the current dependency inversion where `umol-data` imports `DataError` from `umol`.

## GraphIR Error Design

### `ParseError` (DSL syntax errors)

Variants:
- `UnexpectedToken(...)` — nom-level errors converted to syntactic form
- `InvalidTag(...)` — unknown predicate tags
- `DuplicateTag(...)` — repeated predicates
- `OutOfRange(...)` — syntactic range violations (e.g. "value 999 doesn't fit u8")

Nom's low-level errors are mapped via explicit `From` into these syntactic variants. No raw nom types leak.

### `ValidationError` (invariant violations)

Variants:
- `OutOfRange { field, value, min, max }` — semantic range violations
- `ChargeOutOfBounds { element, charge, min_charge, max_charge }`
- `ElectronInvariantMismatch { element, orbital_invariant, electron_invariant }`
- `InvalidMultiplicity(...)` — bad spin multiplicity value
- `SpinUnderdetermined { ... }` — from `SpinStateError::Underdetermined`
- `SpinIncompatible { ... }` — from `SpinStateError::Incompatible`
- `NonGround { field }` — existing variant

`SpinStateError` (from `umol-data`) is mapped via explicit `From<SpinStateError> for ValidationError` that destructures into `SpinUnderdetermined` / `SpinIncompatible`. No wrapping.

### `GraphIrError` (dispatch enum)

```rust
#[derive(Debug, Error)]
enum GraphIrError {
    #[error(transparent)]
    Validation(#[from] ValidationError),
    #[error(transparent)]
    Resolution(#[from] ResolutionError),
    #[error(transparent)]
    Parse(#[from] ParseError),
    #[error(transparent)]
    Aromaticity(#[from] AromaticityError),
    #[error(transparent)]
    Kekulization(#[from] KekulizationError),
    #[error(transparent)]
    Transform(#[from] TransformError),
    #[error(transparent)]
    TopologyExport(#[from] TopologyExportError),
}
```

## Migration: Removing `AtomError` and `BondError`

### `AtomError` decomposition

| Current variant | Disposition |
|---|---|
| `InvalidQueryFormat`, `EmptyQuery`, `InvalidElement` | Dead code (legacy query parser). Delete. |
| `InvalidTag`, `DuplicateTag`, `UnexpectedTag` | → `ParseError` |
| `SpinState(SpinStateError)` | Unwrap. Map `SpinStateError` → `ValidationError` via `From`. |
| `ChargeOutOfBounds`, `OutOfRange`, `ElectronInvariantMismatch` | → `ValidationError` |
| `InvalidMultiplicity` | → `ValidationError` |
| `InvalidCharge`, `InvalidImplicitHydrogens`, `InvalidLonePairs`, `InvalidUnpairedElectrons`, `InvalidValence`, `InvalidDonatedPairs`, `InvalidAcceptedPairs`, `InvalidAromaticValence`, `InvalidMulticenterValence` | Dead code (only used in removed `AtomTypeQuery` parser). Delete. |

### `BondError` decomposition

Same approach: move variants to `ParseError` or `ValidationError` based on concern, delete dead code.

### `Atom::from_str` and other public functions

All public functions in `graph_ir` return `Result<_, GraphIrError>`. Internal functions (`check_invariants`, lowering helpers) return narrow sub-concern errors; `?` promotes via `#[from]` on `GraphIrError`.

## Migration Order

1. Remove `umol-models` crate.
2. Add `DataError` to `umol-data`, remove `umol-data`'s dependency on `umol::error::DataError`.
3. Define `UmolError` trait in umol base crate.
4. Consolidate GraphIR error types into `graph_ir::error`:
   - Define `ParseError` with syntactic variants, add `From` for nom errors.
   - Extend `ValidationError` with invariant variants from `AtomError`, add `From<SpinStateError>`.
   - Define `GraphIrError` dispatch enum.
5. Remove `AtomError` and `BondError`, update all call sites.
6. Change all public GraphIR functions to return `GraphIrError`.
7. Remove `umol::error::Error` god-enum and remaining sub-enums.
