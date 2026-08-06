# 189 — Property-test organization

Status: **Proposed**
Date: 2026-08-06
Relates: [156](156-ast-comparison-and-property-suite-2026-07-20.md),
[161](161-property-tests-as-specs-2026-07-25.md),
[168](168-api-hygiene-2026-07-27.md),
[property-test guide](../docs/development/property-tests.md)

## Purpose

Restructure the property-test suites so that their semantic properties, generated operational
domains, and validation methods are easier to locate and review. This is separate from production
API hygiene: it changes test organization and test-support boundaries, not public Rust or Python
APIs.

The cleanup must preserve the existing properties and their evidence while making later semantic
reviews less dependent on one maintainer remembering where a generator or reference comparison is
hidden. Passing more cases is not the objective; making the executable specification legible is.

## Current structure

Several existing choices should be retained:

- `umol-ast/tests/property.rs` groups the larger molecule, reaction, and stereo subjects
  hierarchically while leaving uniform lattice, entity, delta, and edit law families flat.
- Molecule, reaction, and stereo parent modules own stable regression-file paths so operation-module
  splits do not orphan minimized failures.
- `umol-graph-core` groups properties by graph operation. Its cycle subtree distinguishes
  definition-level exhaustive checks, captured results, and literature fixtures rather than
  presenting them as interchangeable extra cases.
- The test targets and `cargo test --test property -- --list` are the authoritative executable
  inventories. A property-suite README was intentionally removed because it duplicated source and
  drifted.

The primary structural defect is `umol-ast/tests/property/strategies.rs`, currently above 5,000
lines. It combines unrelated generated domains, operation-specific scenarios, reference-adjacent
logic, large transaction fixtures, and a wildcard prelude that re-exports production symbols. Most
property modules therefore import `crate::strategies::*`, obscuring both the public API under test
and the particular operational domain generated for it.

Other large property modules require content review but not automatic splitting. A cohesive family
of laws may remain in one file; line count is evidence that review is needed, not the desired module
count.

## Organization model

Apply the three-axis model from the permanent property-test guide:

| Axis | Organization |
| --- | --- |
| Semantic property | Subject and operation modules |
| Operational domain | Domain-specific strategy modules |
| Validation method | Operation-local reference support, exhaustive support, captured results, or literature fixtures |

The axes must not be collapsed. A generator is not the property definition, and a reference
implementation is not merely another generator helper.

### Property modules

- Retain the present subject-first hierarchy where it already exposes the public operation clearly.
- Split a flat subject by operation only when independent property families are otherwise difficult
  to locate or their operational domains are materially different.
- Keep uniform laws together across entity families when the shared law is the reason the property
  exists.
- State the public properties, generated domains, validation method, and reason for deliberate
  overlap in module documentation. Do not duplicate a test-by-test inventory.
- Keep exact failure properties alongside the public operation whose failure contract they specify,
  not in a generic malformed-input collection unless the malformed domain itself is the subject.

### Strategy modules

Replace the monolithic strategy file with a hierarchy organized by generated domain. An illustrative
shape is:

```text
tests/property/strategy.rs
tests/property/strategy/value.rs
tests/property/strategy/entity.rs
tests/property/strategy/molecule.rs
tests/property/strategy/reaction.rs
tests/property/strategy/edit.rs
tests/property/strategy/metadata.rs
```

The exact child boundaries follow the dependency review; the layout is not a requirement to create
one file for every listed noun.

- Remove production-symbol re-exports from strategy support. Every property module imports the
  public API it tests directly.
- Replace wildcard strategy imports with the specific generators and scenario values used by the
  property module.
- Put a general generator at the lowest domain that owns it. Keep operation-specific scenario
  generation beside that operation unless multiple property families genuinely share it.
- Keep raw, canonical, ground, structurally valid, semantically valid, and deliberately malformed
  strategies distinct and named accordingly.
- Prefer generating independent primitive choices and deterministically assembling a structurally
  valid aggregate. Review broad `prop_filter` use for hidden domain restrictions and poor shrinking.
- For exact failure properties, generate a valid base and introduce one named defect. Carry expected
  information independently when deriving it through the production operation would make the test
  circular.

The strategy surface is internal to the integration-test crate. Use only the visibility needed by
the selected hierarchy and do not widen production visibility to accommodate test organization.

### Reference support and fixtures

- Separate definition-level reference implementations from both generation and production code.
- Locate reference support under the operation it validates rather than creating a global reference
  dumping ground.
- Keep checked-in corpora and captured outputs under the existing test-data facilities, with their
  provenance and regeneration procedures. External tools remain one-off producers, not test
  dependencies.
- Retain regression files at stable subject roots and verify that every child continues to use the
  intended file after movement.
- Do not create `common`, `shared`, or `utils` modules. A support module is named for its generated
  domain, reference operation, fixture family, or other concrete role.

## Review scope

### `umol-ast`

1. Inventory every exported strategy, scenario value, local reference helper, wildcard import, and
   consumer.
2. Classify each item by generated domain or owning operation before proposing the module tree.
3. Split the monolithic strategy file in dependency order, starting with leaf value/entity domains
   and ending with molecules, edits, and reactions.
4. Move operation-specific transaction, composition, application, serialization, and malformed
   scenarios to their owning property areas when they are not genuinely shared.
5. Expand imports in every consumer and remove the strategy prelude.
6. Review the remaining large property modules for cohesive ownership, without forcing a split by
   size alone.

### `umol-graph-core`

1. Preserve the current operation hierarchy and the cycle suite's distinction between exhaustive,
   captured, and literature evidence.
2. Review whether its small common strategy module remains genuinely shared and whether reference
   support is consistently operation-local.
3. Apply the same explicit-import and module-documentation rules; do not reorganize a clear module
   merely for symmetry with `umol-ast`.

### Other crates

Inventory property targets in the remaining workspace crates after the two primary suites establish
the pattern. Apply the pattern where a suite has enough independent domains to benefit; do not build
hierarchy around a handful of cohesive properties.

## Semantic-change boundary

This work preserves property semantics, admitted domains, comparison relations, algorithm selectors,
case configuration, and regression inputs unless a change is separately reviewed. Moving a strategy
must not quietly broaden or narrow it.

If the inventory exposes a mistaken property, duplicated property with no distinct evidentiary role,
missing operational domain, circular expectation, or production semantic defect:

- record the finding explicitly;
- identify the owning public contract and permanent documentation;
- consult before changing the property or production behavior;
- apply an approved correction throughout its relevant implementation, tests, and documentation
  rather than hiding it inside the structural move.

## Verification

- Capture the property inventories and regression-file ownership before movement.
- After each subject migration, run its focused property subtree with an elevated case count and then
  the complete crate property target.
- Compare the before/after inventory by property identity, accounting for intentionally changed
  module paths; do not use a raw count as the sole parity check.
- Run formatting, clippy with property features, and the complete workspace property targets after
  the restructuring.
- Verify by source inventory that the monolithic strategy prelude, production-symbol re-exports, and
  wildcard strategy imports are absent.
- No external program or generated capture becomes a runtime or CI prerequisite.

## Non-goals

- Changing production APIs, module visibility, or crate boundaries.
- Adding a new property-testing framework or replacing proptest.
- Increasing case counts as a substitute for defining properties.
- Creating a property-suite README or hand-maintained coverage matrix.
- Folding missing semantic feature work into a module move without review.
- Splitting every large property file or creating one strategy module per public type.

## Completion criteria

- A property's module path identifies its public subject and operation.
- Each property module imports the production symbols and strategies it actually uses.
- Generated domains are separated and named precisely enough to expose their preconditions.
- Reference implementations, corpora, and literature fixtures have explicit and distinct roles.
- Regression files remain stable and every pre-existing property is accounted for.
- Large remaining files have a recorded cohesive owner rather than being omitted from review.
- The permanent property-test guide and actual suite organization agree.
