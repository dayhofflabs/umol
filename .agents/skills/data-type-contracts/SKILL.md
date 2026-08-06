---
name: data-type-contracts
description: MANDATORY — load and apply before designing, adding, changing, reviewing, or binding a public data type whose construction or use carries non-trivial invariants. Trigger for public constructors; checked/asserted constructor pairs; conversions between representations; validators and transformations; open carriers containing ids, handles, correspondences, entries, or references; operation-issued provenance-bound values; contextual consumers combining independently supplied objects; changes to Option/Result/panic behavior; and Rust/Python construction parity. Establishes the semantic contract before implementation and derives rustdoc and tests from it.
---

# Data type contracts

Make the semantic obligations of a public value explicit before changing its implementation or
surface. Do not introduce new runtime layers merely to describe the contract.

Read `docs/development/data-types.md` and the relevant entries in
`docs/development/nomenclature.md` completely before taking task actions. These living guides are
normative; dated discussion documents are not.

## Inspect the complete surface

Read the type definition and inventory:

- every public constructor and conversion;
- every named validator and transformation;
- every public consumer that combines the value with an independently supplied object or id space;
- algorithmic and internal producers that establish stronger properties by construction;
- Rust and Python bindings;
- rustdoc, exact tests, property tests, examples, and specifications.

Search the whole workspace. Do not infer the contract from one constructor or its nearest tests.

## Classify each property

- **Representation integrity:** the value can be stored and all references in its owning namespace
  resolve. Public construction establishes this.
- **Contextual validity:** the value agrees with a separately supplied graph, molecule, id space,
  model, or operation history. The first public consumer requiring the relationship establishes it.
- **Semantic validity:** a named model-independent or model-dependent predicate. A validator checks
  it by explicit request.
- **Transformation:** canonicalization, repair, stripping, resolution, cascading removal, or another
  change of state. Construction and faithful conversion do not perform it implicitly.
- **Provenance:** an open carrier may be assembled independently and must be checked when combined
  with context. An operation-issued value may rely on the issuing operation's documented contract;
  misuse must not panic but need not produce a semantically correct result.

## Write the contract sheet

Before implementation, state:

```text
Type and role:
Open carrier or operation-issued value:
Intrinsic representation invariants:
Contextual properties and supplied context:
Semantic predicates and validators:
Public constructors:
Conversions and preserved information:
Explicit transformations:
First public consumer requiring each contextual property:
Failure, absence, and panic behavior:
Algebraic, preservation, or roundtrip properties:
Rust/Python boundary:
```

Resolve semantic gaps and public names with the user before editing. A staged implementation plan
sequences a settled contract; it does not settle one.

## Shape the API directly

- Constructors establish intrinsic representation integrity, not every property that can be checked.
- Pair `from_x` for asserted producers with `try_from_x` for independently assembled input when both
  are useful. Both establish the same invariant and differ only in failure reporting.
- Put contextual fallibility on the first public operation whose promised result requires the
  context. Do not propagate it through unrelated operations because of an internal call path.
- Use `Option` for one ordinary absence condition with no useful cause; use an operation-specific
  `Result` when callers can act on distinct causes.
- Preserve every source state representable by the target. Conversion does not silently drop,
  canonicalize, repair, resolve, or select an interpretation.
- Keep asserted internal paths only where a producer establishes the property by construction, and
  document that producer contract.
- Do not add `Validated<T>`, prepared/witness types, typestate, internal flags, or parallel checked
  APIs unless the state is a distinct reusable concept that materially clarifies ordinary use.

## Document the contract in rustdoc

- Start with the standard one-line summary and explanatory prose. State intrinsic invariants and
  important contextual limitations there; headings are not required merely to label them.
- Use `# Errors`, `# Panics`, and `# Safety` exactly when applicable. Document ordinary `None`
  conditions in the main prose.
- Use `# Examples` where an example explains why or how the API is used.
- Use `# Semantic properties` for non-trivial laws, invariants, soundness/completeness, or roundtrip
  guarantees. State preconditions and the equality or equivalence relation.
- Rustdoc accepts ordinary Markdown headings; there is no closed vocabulary. Keep custom headings
  rare and repository-wide rather than inventing a local taxonomy for one type.
- Keep public rustdoc self-contained. Do not cite dated discussion documents from source code.

## Derive verification from the contract

- Add exact unit cases for each meaningful success and failure boundary.
- Add properties for algebraic laws, preservation, roundtrips, and contextual soundness.
- Prefer a definition-level reference implementation, exhaustive bounded domain, independent
  implementation, or metamorphic relation over repeating production logic.
- State why overlapping properties use different operational domains or validation methods.
- Verify trusted producer paths remain infallible and externally reachable mismatches cannot become
  indexing panics.
- Migrate all callers, bindings, tests, specifications, examples, benchmarks, and fuzz targets
  affected by the settled contract.

Passing tests are evidence for the stated contract, not a substitute for stating it.
