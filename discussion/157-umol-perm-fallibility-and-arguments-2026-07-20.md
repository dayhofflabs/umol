# umol-perm fallibility and argument design

Status: **Active**
Date: 2026-07-20
Relates: [109](109-permutation-infrastructure-2026-06-09.md),
[119](119-umol-perm-review-2026-06-21.md),
[156](156-ast-comparison-and-property-suite-2026-07-20.md)

## Scope

Reaction application needs stereo reframing to reject incompatible ligand
frames without allowing `Permutation::between` to panic. The same audit found
several other public operations whose failure behavior or arguments are
unclear.

The solution is not to make every operation whose domain is narrower than its
Rust arguments return `Result`. That would turn ordinary small-degree algebra
into a chain of error propagation without providing useful recovery. The API
distinguishes:

- construction from runtime structural data;
- queries or derivations for which absence or incompatibility is an ordinary
  outcome;
- algebraic operations whose receiver and documented argument contract define
  their domain.

The first two categories may be fallible. The third remains infallible with
asserted programmer contracts. Crate-private mathematical invariants may also
continue to assert.

Construction and parsing use separate errors:

- `PermutationError` reports malformed runtime permutation images and cycles;
- `ParseClassKeyError` reports invalid `ClassKey` text.

Both errors are ordinary public Rust error types with `Display` and
`std::error::Error` implementations. `umol-perm` remains dependency-free; the
migration does not add an error-derive dependency.

## Fixed maximum degree

The maximum permutation degree is fixed by the largest supported stereo
`ClassKey`. It defines the domain of the compact representation and is not
runtime policy.

Exceeding that degree remains an assertion for ordinary constructors and
algebraic operations. A conversion whose stated purpose is to accept arbitrary
external image data rejects an image that cannot fit the representation; this
does not make the maximum configurable runtime policy. Neither a validated
degree wrapper nor a const-generic permutation is introduced.

## Existing fallible surface

Eight public operations currently return `Option` or `Result`:

| Operation | Return | Meaning |
| --- | --- | --- |
| `ClassKey::from_str` | `Result<ClassKey, String>` | Invalid spelling or degree. |
| `CosetSpace::coset_rep` | `Option<Permutation>` | The permutation has the wrong degree. |
| `CosetSpace::index` | `Option<u32>` | The permutation has the wrong degree or is not numbered by the space. |
| `CosetSpace::unindex` | `Option<Permutation>` | The index is out of range. |
| `CosetSpace::reindex` | `Option<u32>` | The index or relabeling is invalid for the space. |
| `CosetSpace::enantiomer` | `Option<u32>` | The index is out of range. |
| `CosetSpace::observable_coset` | `Option<u32>` | The index is out of range. Its current implementation can still panic through `merge_under` for an invalid fluxional generator. |
| `OrientedPermutationGroup::improper_rep` | `Option<OrientedPermutation>` | The group has no improper component. This is optional structure, not operation failure. |

Seven existing operations therefore signal that the requested result cannot be
produced. `improper_rep` merely exposes optional data. The established `rep`
abbreviation is retained in these names, consistent with APIs such as
`class_reps`.

The proposed surface has 11 failure-signalling operations: the seven existing
ones plus fallible image conversion, `from_cycles`, `between`, and
`orbit_reps` (currently `merge_under`). `from_image` and `unrank` remain
infallible. Changing how `observable_coset` obtains its `Option` does not add
another operation to the count.

## Return-shape rule

- Use `Option` for bounded lookup where the useful distinction is present or
  absent.
- Use `Result` for malformed runtime structure or incompatible arguments where
  the reason is useful to diagnostics and Python conversion.
- Use `false` for membership in an incompatible degree or domain.
- Use an assertion for violations of the fixed representation domain or of an
  ordinary algebraic operation's documented contract.

## `Permutation`

| Operation | Proposed API | Reason |
| --- | --- | --- |
| `identity(degree)` | Unchanged and infallible. | The fixed maximum degree is an asserted representation invariant. |
| `from_image(image)` | Remove the redundant degree argument; accept `&[usize]`; remain infallible with an asserted validity contract. | Degree is exactly `image.len()`. Public points and cycle notation already use `usize`; the packed internal integer width should not leak into the API. This is the ordinary construction path for fixed class definitions and already-established images. |
| fallible image conversion | Implement standard fallible conversion from `&[usize]`. | Python and other runtime-data boundaries need image range, bijectivity, and maximum-length failures as values. Keeping this separate avoids adding assertions to every static and internally derived valid image. |
| `degree` | Unchanged. | Total observation. |
| `apply(point)` | Unchanged and infallible. | `point < degree` is the ordinary domain contract of applying the permutation. |
| `act(items)` | Remain infallible; require `items.len() == degree` and assert otherwise; weaken `T: Copy` to `T: Clone`. | The current acceptance and truncation of longer slices obscures a mismatch. Acting on owned values requires cloning, not the stronger `Copy` bound. |
| `compose(other)` | Unchanged and infallible. | Equal degree is the ordinary composition contract. Returning `Result` would infect every algebraic call site. |
| `inverse`, `sign`, `rank`, `cycles` | Unchanged. | Total on a valid receiver. |
| `between(from, to)` | Return `Option<Self>`. | Two runtime frames may differ in length, contain repeated members, or have unequal membership. No caller distinguishes those causes; the useful result is whether a relabeling exists. Exceeding the fixed maximum remains an assertion. |
| `unrank(degree, rank)` | Remain infallible with an asserted rank contract. | Every call site constructs `rank` in `0..degree!`: symmetric/alternating group enumeration, property generators, or `rank` round trips. Both arguments remain necessary. |
| `from_cycles(degree, cycles)` | Keep both arguments and return `Result<Self, _>`. | Degree cannot be inferred because unmentioned points are fixed. Runtime cycles may contain out-of-range or repeated positions. Exceeding the fixed maximum remains an assertion. |

## Groups and oriented permutations

| Operation family | Proposed behavior | Reason |
| --- | --- | --- |
| `PermutationGroup::generate` | Remain infallible; assert that every generator has the requested degree. Retain the explicit degree. | The degree is necessary for an empty generator set. Degree agreement is the algebraic construction contract. Splitting out a separate trivial-group path would make dynamic generator lists less convenient. |
| `symmetric`, `alternating`, `cyclic`, `dihedral` | Remain infallible with documented degree contracts. | They construct named groups inside the fixed representation domain. Degree-zero behavior for cyclic and dihedral groups must be documented rather than discovered through arithmetic failure. |
| `PermutationGroup::contains` | Keep `bool`; return `false` for a wrong-degree permutation. | Membership is a predicate. |
| `OrientedPermutation::identity`, `apply`, `compose` | Remain infallible with the same contracts as their underlying permutation operations. | The oriented wrapper does not introduce a new runtime failure class. |
| `OrientedPermutationGroup::generate` | Remain infallible; assert generator degree agreement. | Same group-construction contract as `PermutationGroup::generate`. |
| `OrientedPermutationGroup::contains` | Keep `bool`; return `false` for a wrong-degree operation. | The current improper branch can panic while trying to compose different degrees; membership should instead be false. |
| `proper_orbit_of`, `star_orbit_of` | Remain infallible; assert that the point lies in the group domain. | Orbit calculation is ordinary algebra over a documented point domain. |
| Other observations and transformations | Unchanged. | Total on valid receivers. |

## Classes, cosets, and coset spaces

| Operation | Proposed behavior | Reason and argument assessment |
| --- | --- | --- |
| `ClassKey::from_str` | Return `Result<Self, ParseClassKeyError>`. | Text parsing is an external-data boundary. The current `String` error was the direct formatting used by the initial implementation, not a deliberate public error contract. Class-key syntax errors are separate from permutation construction errors. |
| `space(key)` | Move to the infallible `ClassKey::space(self)`. | The class key is the natural receiver. Class construction stays inside the fixed representation domain. |
| `Coset::new(key, index)` | Remain an infallible checked constructor that asserts `index < count`. | `Coset` otherwise has no ordinary infallible construction path, and current external stereo representations do not construct it directly from untrusted indices. |
| `Coset::new_unchecked` | Unchanged. | It remains the explicit bypass for callers that already establish the index invariant. |
| `CosetSpace::coset_rep`, `index`, `unindex`, `reindex`, `enantiomer` | Keep their current `Option` returns. | These are bounded lookups. |
| `CosetSpace::merge_under(generators)` | Rename to `orbit_reps(generators)` and return `Option<Vec<u32>>`. | The result contains the canonical orbit representative for every coset; `orbit_reps` follows the established `class_reps` abbreviation. Runtime generators may not belong to the parent group. No current caller distinguishes causes, so absence is sufficient. A generator slice is the correct argument because the operation computes their generated closure; requiring a materialized group would do extra work without removing the parent-membership check. |
| `CosetSpace::observable_coset(index, fluxional)` | Keep `Option<u32>` and propagate absence from `orbit_reps`. | An absent index or generator outside the parent group means no observable coset. Current graph-symmetry callers already use that meaning. |
| Other observations | Unchanged. | Total on valid receivers. |

## Construction and acquisition audit

No public type is left without either an infallible constructor or its ordinary
infallible acquisition path:

| Type | Infallible construction or acquisition |
| --- | --- |
| `Permutation` | `identity`; additionally obtained from class spaces and other valid algebraic operations. |
| `PermutationGroup` | `generate`, `symmetric`, `alternating`, `cyclic`, `dihedral`. |
| `Orientation` | Enum variants. |
| `OrientedPermutation` | `new`, `proper`, `improper`, `identity`. |
| `OrientedPermutationGroup` | `generate`. |
| `ClassKey` | Enum variants. |
| `Coset` | Checked `new`; `new_unchecked` for an already-established index invariant. |
| `CosetSpace` | `ClassKey::space`. |

Fallible image conversion therefore supplements rather than replaces
`from_image`, the ordinary construction path for `Permutation`. `Coset::new`
does not become fallible merely to eliminate its assertion.

## Error vocabulary

`PermutationError` has the following variants:

- `ImageTooLong { length, maximum }`;
- `ImageValueOutOfRange { position, value, degree }`;
- `DuplicateImageValue { value }`;
- `CyclePointOutOfRange { cycle, position, point, degree }`;
- `DuplicateCyclePoint { point }`.

`ParseClassKeyError` has the following variants:

- `UnknownClassKey { input }` for an unknown fixed key or family;
- `InvalidDegree { input }` for a missing or malformed degree;
- `DegreeTooLarge { degree, maximum }` for a parsed degree beyond the fixed
  representation width.

The parser reads a family degree as `usize` before checking the fixed maximum,
so a numeric value larger than `u8` is still classified as `DegreeTooLarge`
rather than as malformed text.

`between` and `orbit_reps` use `Option` and therefore need no public error
variants. Their callers add operation-specific meaning such as
`StereoFrameMismatch` where necessary.

## Call-site migration

### Image construction

| Call sites | Handling |
| --- | --- |
| `umol-perm` group rotation/reflection, class definitions, and coset representatives | Continue through infallible `from_image`; remove the degree argument and use `usize` images. These formulas and fixed tables establish bijectivity. |
| Fixed permutations in `umol-ast`, validators, tests, and Python conversion fixtures | Continue through infallible `from_image`; remove the degree argument. |
| `umol-ast` property strategies | Continue through infallible `from_image`; the shuffled complete range establishes a valid image. |
| `project_onto_ligands` | Use fallible image conversion and propagate `None`; the constructed projection may fail to be injective for malformed ligand data. |
| Python `Permutation(image)` | Use fallible image conversion and map construction failure to `ValueError`. |

### Cycle and rank construction

| Call sites | Handling |
| --- | --- |
| Text DSL, tree EDN, and streaming EDN cycle readers | Propagate `from_cycles` failure through their existing syntax/deserialization errors. Negative integer conversion remains local to the readers. |
| Symmetry-generated transpositions in `observable_descriptor` and `virtual_block_swaps` | Use an explicit invariant assertion: loop bounds and adjacent position selection establish valid disjoint cycles. |
| DSL expected values that currently call `from_cycles` | Construct expected permutations with infallible `from_image`; this avoids using the parser's construction path as its own oracle. |
| Symmetric and alternating group enumeration | Keep infallible `unrank`; `0..degree!` establishes the rank contract. |
| Permutation property strategies and rank round trips | Keep infallible `unrank`; their generated or recovered ranks establish the contract. |

### Frame relabeling

| Call sites | Handling |
| --- | --- |
| `StereoAtomAst::transform_frame` and `StereoBondAst::transform_frame` | Return `Option<Self>` and propagate absence from `between`. |
| Molecule pushout | Make stereo glue collection fallible and propagate `None` through the existing optional pushout result. |
| Reaction reframing | Return `Result<(), ApplyError>` and map `between`/`transform_frame` absence to `ApplyError::StereoFrameMismatch` for the affected stereo entity. |
| `permutation_for_ligands` and symmetry `reexpress` | Propagate `between` directly as `Option`; retain only their independent coset-index handling. |
| TableIR tetrahedral and cis/trans raising | Retain explicit invariant assertions. The source and target orders are derived from the same ligand set after geometry validation; failure would be an implementation defect rather than malformed syntax. |
| `umol-perm` and `umol-ast` tests | Positive cases compare with `Some(expected)`; incompatible frames use exact `None` cases. |

### Coset merging

| Call sites | Handling |
| --- | --- |
| `CosetSpace::observable_coset` | Propagate `None` from `orbit_reps`, then perform the existing bounded index lookup. Its public signature remains `Option<u32>`. |
| `MoleculeAst::observable_descriptor` | Propagate `observable_coset` directly; malformed coset or generator data yields no descriptor. |
| `StereoSymmetry::is_stereogenic` | Return `false` when merging fails or the stored coset index is out of range, consistent with malformed data not establishing stereogenicity. |
| Unit and property tests | Compare valid merges with `Some(expected)` and add exact `None` cases for generators outside the parent group. |

## Redundant validation and retained assertions

The lower-level fallible operations make these caller-side checks redundant:

- the explicit range/disjointness loop in the text DSL cycle parser;
- the `seen` vectors and range/repetition filters in both EDN cycle readers;
- the length, uniqueness, and set-equality checks in
  `permutation_for_ligands`;
- the length, uniqueness, and set-equality checks in symmetry `reexpress`;
- the separate injectivity tracking in `project_onto_ligands`, replaced by
  fallible image conversion;
- the explicit stored-coset range check in `observable_descriptor`, now covered
  by `observable_coset`;
- the caller-controlled `expect("generator stays in the parent group")` in
  the current `merge_under`, replaced by `None` in `orbit_reps`;
- the `#[should_panic]` cycle-construction test, replaced by exact construction
  error assertions.

These assertions and validations remain:

- the fixed maximum-degree assertions;
- the point, slice-length, and degree-agreement contracts of `apply`, `act`, and
  `compose`;
- group-generator degree assertions;
- internal `0..count` coset-loop assertions;
- TableIR tetrahedral/cis-trans geometry validation and the subsequent
  same-ligand-set `between` assertions;
- reaction correspondence checks, which S1c handles independently;
- the asserted invalid-image contract of infallible `from_image`, alongside
  exact error tests for the fallible image conversion.

## Boundary consequences

- `ReactionAst` stereo reframing calls optional `between` and maps frame
  incompatibility to `ApplyError::StereoFrameMismatch`.
- DSL parsing propagates fallible cycle construction through its own syntax
  errors.
- Python maps malformed runtime permutation data to `ValueError`; it checks the
  fallible image conversion rather than calling asserted `from_image`.
- Internal literal permutations and algebraic operations may use explicit
  assertions where their inputs are statically evident or already validated.

The fallible `between` migration is a prerequisite of S1c in document 156. The
remaining API changes form a separate crate-wide migration.

## Staged implementation plan

Every changed test uses the workspace test conventions: `rstest` tables for
input classes, exact values or error variants rather than `is_err`/`is_some`,
and expected permutations constructed independently of the operation under
test. Unit tests stay ordered with the methods they cover; algebraic laws remain
in the property suites.

### S0 — additive error and lookup foundation

- **S0a — `PermutationError`** (`umol-perm/src/error.rs`, `src/lib.rs`): add
  the five settled variants, manual `Display` and `std::error::Error`
  implementations, and the public re-export without adding a dependency.
  Table tests pin each diagnostic and its stored fields.
  **Implemented (green).**
- **S0b — `ParseClassKeyError`** (`umol-perm/src/error.rs`, `src/lib.rs`): add
  the three settled variants and manual standard-error implementations. Table
  tests distinguish unknown keys, invalid degrees, and excessive degrees.
  **Implemented (green).**
- **S0c — typed `ClassKey` parsing** (`umol-perm/src/class.rs`): change
  `FromStr::Err` from `String` to `ParseClassKeyError`, parse family degrees as
  `usize`, and return the exact variant for every invalid syntax class. Replace
  presence-only error assertions with exact expected errors and retain display /
  parse round trips for every fixed and parameterized family.
  **Implemented (green).** [dep: S0b]
- **S0d — fallible image conversion** (`umol-perm/src/permutation.rs`,
  `tests/property.rs`): implement `TryFrom<&[usize]> for Permutation` through one
  checked image-construction path. Unit tables cover all five relevant image
  outcomes (valid, too long, out of range, duplicate, identity); property tests
  cover valid shuffled images and image/rank round trips.
  **Implemented (green).** [dep: S0a]
- **S0e — `ClassKey::space`** (`umol-perm/src/class.rs`): add the inherent
  method while retaining the free `space` function temporarily. Move the
  interning implementation behind the method and test counts, pointer identity,
  and acquisition through `Coset` without using the free function as the test
  oracle.
  **Implemented (green).**
- **S0f — `CosetSpace::orbit_reps`** (`umol-perm/src/coset.rs`,
  `tests/property.rs`): add the optional operation with the current merge
  algorithm, returning `None` before composition when any generator lies
  outside the parent group. Keep `merge_under` temporarily as a compatibility
  wrapper over the valid-input path. Unit tables assert complete representative
  vectors and exact `None`; properties cover identity generators and closure
  under valid generated actions.
  **Implemented (green).**

S0 verification: `cargo test -p umol-perm --features proptest` and
`cargo clippy -p umol-perm --all-targets --features proptest -- -D warnings`
pass, as does `cargo check --workspace --lib --all-features`. The workspace
all-target check is currently blocked in an existing `umol-ast` property test
that still calls the old one-argument `DpoValidator::validate_reaction` API;
that failure is outside the S0 changes.

### S1 — image API migration

- **S1a — infallible `from_image` contract**
  (`umol-perm/src/permutation.rs`, `group.rs`, `class.rs`, `coset.rs`,
  `oriented.rs`): change the signature to `from_image(image: &[usize])`, infer
  degree from the slice, and share construction with the checked image path
  while preserving an asserted contract. Migrate all `umol-perm` formulas,
  fixed class tables, examples, unit tests, and property generators; retain an
  exact panic-contract test independently of the fallible-conversion tests.
  **Implemented (green).** [dep: S0d]
- **S1b — Rust consumer migration** (`umol-ast`, `umol-graph`, `umol-io`):
  remove degree arguments and packed-width casts from fixed images and property
  strategies. Change `project_onto_ligands` to `Permutation::try_from`, remove
  its separate injectivity tracking, and propagate `None`; add exact projection
  cases for a valid image and a collision.
  **Implemented (green).** [dep: S1a]
- **S1c — Python boundary migration** (`umol-py/src/stereo.rs`, `delta.rs`,
  Python tests): accept Python image values as `usize`, construct through
  `TryFrom<&[usize]>`, and map every `PermutationError` to `ValueError`.
  Migrate Rust-side Python fixtures to the new infallible signature and test
  exact valid images plus each invalid Python input class. This restores the
  whole workspace to green.
  **Implemented (green).** [dep: S1a, S1b]

S1 verification: `cargo test -p umol-perm --features proptest`, the focused
`umol-ast` projection tests, Rust test compilation for `umol-ast`, `umol-graph`,
`umol-io`, and `umol-py` pass. Clippy also passes for all workspace libraries
with all features and warnings denied. After rebuilding the extension in the
Python 3.13 virtual environment, all 900 Python tests pass. The complete
`umol-ast` property target remains blocked by the unrelated pre-existing
`DpoValidator::validate_reaction` call noted after S0.

### S2 — algebraic contracts

- **S2a — sequence action** (`umol-perm/src/permutation.rs`): change `act` to
  `T: Clone`, require `items.len() == degree`, and clone selected items. Test a
  non-`Copy` element type, ordinary action, and both short and long contract
  violations.
  **Implemented (green).** [dep: S1a]
- **S2b — permutation-group contracts** (`umol-perm/src/group.rs`): make
  `generate` assert generator-degree agreement before closure, make
  wrong-degree `contains` explicitly return `false`, and give `cyclic` and
  `dihedral` an explicit nonzero-degree assertion instead of allowing modulo by
  zero to fail incidentally. Add exact membership and contract tables; preserve
  order tests for valid named groups.
  **Implemented (green).** [dep: S1a]
- **S2c — oriented-group contracts** (`umol-perm/src/oriented.rs`): assert
  generator degree agreement in `generate`, return `false` from `contains` for
  either orientation at the wrong degree, and assert the point domain directly
  in `proper_orbit_of` and `star_orbit_of`. Test the proper and improper
  membership branches and both orbit families.
  **Implemented (green).** [dep: S2b]

S2 verification: the complete `umol-perm` unit/property suite and all 5,085
`umol-ast` unit tests pass. Clippy passes for every `umol-perm` target with the
property feature and for all workspace libraries with all features, with
warnings denied.

### S3 — cycle construction and readers

- **S3a — fallible `from_cycles`** (`umol-perm/src/permutation.rs`,
  `tests/property.rs`): return `Result<Self, PermutationError>`, detect
  out-of-range and duplicate points before writing the image, and retain the
  asserted fixed-maximum contract. Construct expected values with
  `from_image`, replace the old panic test with exact error tables, and preserve
  cycle-decomposition round-trip properties.
  **Implemented (green).** [dep: S0a, S1a]
- **S3b — text DSL cycles** (`umol-ast/src/dsl/stereo.rs`): propagate
  `PermutationError` through the existing `ParseError::InvalidValue` boundary,
  remove the duplicate range/disjointness validation, and replace parser
  expected values built with `from_cycles` by fixed images. Parsing tables cover
  valid, overlapping, repeated, and out-of-range cycle notation with exact
  errors.
  **Implemented (green).** [dep: S3a]
- **S3c — EDN cycle readers** (`umol-ast/src/dsl/stereo.rs`,
  `dsl/constraint.rs`): make both tree and streaming readers delegate structural
  validation to `from_cycles`, retaining only EDN type and negative-integer
  checks, and map construction errors through `DeError::Custom`. Add matching
  exact malformed-cycle tables to both reader paths and use fixed images as
  successful expected values.
  **Implemented (green).** [dep: S3a]
- **S3d — generated internal cycles** (`umol-ast/src/ast/symmetry.rs`): make
  the loop-derived transpositions in `observable_descriptor` and
  `virtual_block_swaps` use explicit invariant assertions on the fallible
  result. Tests cover the generated virtual-block images directly and the
  observable effect of a class-pair transposition. The property suite retains
  its coverage of observable descriptors and virtual ligand blocks. This
  restores the workspace to green.
  **Implemented (green).** [dep: S3a, S3b, S3c]

S3 verification: all 150 `umol-perm` unit tests, its 17 property tests, and all
5,107 `umol-ast` unit tests pass. The workspace library check and Clippy for all
`umol-perm` targets and the `umol-ast` library pass with all features and
warnings denied. The complete `umol-ast` property target compiles, and the
composition dangling-invariant property passes.

### S4 — class-space API migration

- **S4a — inherent-space callers** (`umol-perm/src/class.rs`, `coset.rs`,
  `umol-ast/src/ast/stereo.rs`, `symmetry.rs`, `umol-io/src/table_ir/raise.rs`):
  replace every `space(key)` call with `key.space()` while both entry points
  still exist. Update affected tests to acquire the same interned space through
  the receiver API.
  **Implemented (green).** [dep: S0e]
- **S4b — retire free `space`** (`umol-perm/src/class.rs`, `src/lib.rs`): remove
  the free function and its re-export after the workspace has no callers;
  retain the registry and its poisoning assertion behind `ClassKey::space`.
  Run the class/coset unit and property suites to pin counts, interning, and
  `Coset` identity.
  The property audit added generated checks for repeated `ClassKey::space`
  pointer identity and for `Coset::space` returning that same interned value.
  **Implemented (green).** [dep: S4a]

S4 verification: all 150 `umol-perm` unit tests and 19 property tests pass, as
do all 5,107 `umol-ast` and 3,288 `umol-io` unit tests. The workspace library
check and Clippy pass with all features and warnings denied.

### S5 — coset-orbit API migration

- **S5a — observable cosets** (`umol-perm/src/coset.rs`,
  `tests/property.rs`): implement `observable_coset` through `orbit_reps` and
  propagate `None` before the bounded index lookup. Replace valid merge
  expectations by `Some(expected)` and add wrong-parent generators as exact
  `None` cases.
  Property tests cover generated valid-generator delegation and distinguish
  invalid indices from wrong-degree generators.
  **Implemented (green).** [dep: S0f]
- **S5b — symmetry consumers** (`umol-ast/src/ast/symmetry.rs`): replace
  `merge_under` with `orbit_reps`; remove the stored-coset precheck now covered
  by `observable_coset`; make `StereoSymmetry::is_stereogenic` return `false`
  when the orbit representatives or stored index are invalid. Unit and property
  tests distinguish valid singleton orbits, merged orbits, invalid generators,
  and invalid coset indices.
  **Implemented (green).** [dep: S5a]
- **S5c — retire `merge_under`** (`umol-perm/src/coset.rs`): remove the
  compatibility method after all callers use `orbit_reps`; rerun the complete
  coset unit/property suite against the retained `coset_rep` and
  `improper_rep` names.
  **Implemented (green).** [dep: S5b]

S5 verification: all 152 `umol-perm` unit tests and 21 property tests pass, as
do all 5,112 `umol-ast` unit tests and its complete 161-test property suite. The
workspace library check and Clippy pass across all targets and features with
warnings denied.

### S6 — frame relabeling and application

- **S6a — optional `Permutation::between`** (`umol-perm/src/permutation.rs`,
  `class.rs`, `tests/property.rs`): return `Option<Self>`,
  reject unequal lengths, repeated members, and unequal membership without
  panicking, and retain the fixed-maximum assertion. Exact unit tables cover
  every incompatibility; properties assert that successful relabeling acts the
  source frame into the target and round-trips through the reverse frame.
  **Breaking (red).** [dep: S1a]
- **S6b — stereo AST and view relabeling** (`umol-ast/src/ast/stereo.rs`,
  `ast/view/stereo.rs`, `ast/symmetry.rs`): return `Option<Self>` from both
  `transform_frame` implementations; let `permutation_for_ligands` delegate
  compatibility to `between`; inline the one-use `reexpress` helper and retain
  only its independent coset lookup. Tests compare successful values with
  `Some(expected)` and incompatible frames with exact `None`.
  **Breaking caller migration (red).** [dep: S6a, S4b]
- **S6c — molecule pushout** (`umol-ast/src/ast/molecule/pushout.rs`): make
  stereo glue-entry construction fallible and propagate incompatible frames
  through the existing optional pushout result. Add atom- and bond-stereo
  pushout tables for reordered compatible frames and changed ligand sets.
  **Breaking caller migration (red).** [dep: S6b]
- **S6d — reaction reframing** (`umol-ast/src/ast/reaction.rs`, `ast/error.rs`):
  make `reframe_stereo` return `Result<(), ApplyError>`, map failed `between` or
  `transform_frame` to `StereoFrameMismatch` for the affected entity, and
  propagate the fatal error once through the application iterator. Regression
  tables cover atom and bond field changes/removals and assert both the exact
  entity-bearing error and permanent iterator termination.
  **Breaking caller migration (red).** [dep: S6b]
- **S6e — TableIR invariant callers** (`umol-io/src/table_ir/raise.rs`): keep
  tetrahedral and cis/trans geometry validation, then explicitly assert that
  `between` succeeds because both frames were derived from that validated
  ligand set. Conformance tests retain fixed successful references and add no
  user-facing resolution error for an internal invariant. This restores the
  workspace to green.
  **Breaking caller migration (red→green).** [dep: S6a]

### S7 — property tests

- **S7a — class-key text round trips** (`umol-perm/tests/property.rs`): extend
  the class-key strategies across every fixed class and every representable
  degree of the symmetric, alternating, cyclic, and dihedral families. Assert
  `ClassKey::from_str(key.to_string()) == Ok(key)` for the full generated
  domain; retain the exact malformed-input tables as unit tests.
  **Additive (green).** [dep: S0c]
- **S7b — independent image round trips** (`umol-perm/tests/property.rs`):
  generate shuffled images directly for every supported degree, construct with
  `Permutation::try_from`, and assert both the inferred degree and the complete
  recovered image. Replace the rank-derived image property rather than keeping
  two tests of the same invariant; retain the separate Lehmer-rank round trip.
  **Additive (green).** [dep: S0d, S1a]
- **S7c — canonical coset orbits** (`umol-perm/tests/property.rs`): generate
  arbitrary valid multi-generator sets for each fixed class space and compare
  every `orbit_reps` entry with the minimum index reached by an independent
  traversal under those generators. This covers multi-generator closure,
  representative membership, and the canonical-minimum contract; retain the
  exact invalid-generator unit tables.
  **Additive (green).** [dep: S0f]
- **S7d — sequence-action composition** (`umol-perm/tests/property.rs`): for
  equal-degree permutations and generated sequences, assert
  `a.compose(b).act(items) == b.act(&a.act(items))`, matching the documented
  composition and right-action conventions. Keep the non-`Copy` and slice-size
  contracts in the unit tests; do not add a redundant inverse-action property.
  **Additive (green).** [dep: S2a]
- **S7e — generated permutation groups** (`umol-perm/tests/property.rs`): retain
  arbitrary valid generator lists alongside each generated `PermutationGroup`.
  Assert that the identity and every supplied generator are members and that
  the returned elements are closed under inverse and composition. Keep
  wrong-degree inputs and named-group orders in the exact unit tables.
  **Additive (green).** [dep: S2b]
- **S7f — generated oriented groups** (`umol-perm/tests/property.rs`): retain
  the input generators in the existing oriented-group strategy and assert that
  every supplied proper or improper generator belongs to the result. Preserve
  the existing identity, inverse, and composition closure properties without
  duplicating them; keep wrong-degree membership and orbit-domain contracts in
  the unit tests.
  **Additive (green).** [dep: S2c]

The critical consumer path joins **S0d → S1** and **S0e → S4** before
**S6a → S6b → S6d**; S6 is the prerequisite for document 156 S1c. The cycle
path joins **S0a** and **S1a** at S3, and the coset path is **S0f → S5**. S2 is
independent once S1 lands. S7 is additive and may run at any point after S2;
it is part of the migration's completion rather than a dependency of S3–S6.
No stage is deferrable within this API migration; each stage ends green,
although the explicitly breaking subitems may be red until their remaining
subitems restore all callers.

Final verification runs formatting and `git diff --check`, the complete
`umol-perm` unit/property suites, affected `umol-ast`, `umol-graph`, and
`umol-io` suites with all relevant features, Python tests from the activated
`umol-py/.venv`, and workspace Clippy across all targets.
