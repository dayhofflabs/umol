# DSL and EDN worklist

Status: **Proposed**
Date: 2026-07-27
Relates: [105](105-dsl-fixes-2026-06-06.md),
[115](115-variable-facility-2026-06-16.md),
[153](153-format-parsing-outstanding-tasks-2026-07-18.md),
[169](169-dsl-metadata-context-2026-07-27.md)

## Scope

This document tracks the remaining work on the native umol DSL and its EDN
representation. It replaces the DSL and tooling lists in doc 105 after removing
completed items and routing work owned by other subsystems to their focused
documents.

The three relevant layers are distinct:

1. entity-string grammars such as atom, bond, value, and stereo strings;
2. the tree-shaped EDN representation used by `FromEdn` and `ToEdn`;
3. streaming or source-preserving parsing used for diagnostics and large inputs.

Changes must preserve this separation. A streaming parser is another execution
path for the same language, not a second DSL.

Persistent metadata, parse-time contexts, and metadata-preserving Python
parse/render operations are owned by doc 169 rather than this remaining-work
inventory.

## Language and specification

### Construction forms

- Design graph construction shorthand for rings, chains, and fragments. The
  shorthand must have a defined expansion into the ordinary molecule DSL before
  syntax is selected. Traversal notation and construction notation must not be
  conflated.
- Revisit stereo-constraint shorthand such as free rotation, Berry
  pseudorotation, and ring flips only after the corresponding constraint
  semantics are settled. The shorthand must lower to ordinary AST constraints.
- Decide whether stereo membership constraints need a more readable surface
  than the current membership operators. This is a syntax decision, not a new
  constraint model.

### Specification conformance

- Sweep the specification and examples for residual `ref` and `bind`
  terminology. Variables use the `var` vocabulary; historic names must not
  survive in prose, grammar productions, methods, or test names.
- Update the normative-keyword declaration from RFC 2119 to RFC 8174 and verify
  that every `MUST` and `MUST NOT` expresses an implemented, testable rule.
- Add or identify property tests for normative grammar laws. Table tests remain
  appropriate for individual examples; property tests should cover laws such
  as whitespace insertion, map-key permutation, and parse/render stability.
- Keep element casing, isotope grammar, numerical ranges, and stereo grammar in
  explicit subsections. Validate that the grammar and the integrity validators
  agree on their accepted ranges.
- Decide and document the EDN representation of cycles before adding a public
  parser or formatter for them.

Variable scope and cross-object variable references are owned by doc 115 and
are not redesigned here.

## Conversion and entry-point APIs

- Audit DSL-to-AST and AST-to-DSL crossings. Use `IntoAst` and `FromAst` for
  semantic conversion instead of reaching through tuple-wrapper `.0` fields.
  Direct field access inside a wrapper's own parser or formatter is not itself
  a conversion defect.
- Inventory the current string and EDN entry points and settle one naming
  convention for:
  - parsing an entity-string grammar;
  - reading a tree-shaped EDN value;
  - reading EDN source text;
  - streaming multiple forms.
- Separate parser and formatter modules where the tree and streaming paths have
  different state, while sharing the grammar and conversion logic.
- Add streaming entry points for constraint forms once the common streaming
  boundary is settled.
- Complete structured parse errors for missing mandatory entity fields,
  especially missing localized- and dative-bond order.

External SMILES, MOL, and SDF entry points remain in umol-io and are tracked by
docs 112 and 153.

## Source-preserving EDN

Editor diagnostics require source locations without changing the existing
runtime EDN value:

```rust
struct Span {
    start: usize,
    end: usize,
}

struct Spanned<T> {
    value: T,
    span: Span,
}
```

The intended split is:

```text
source text
├── parse_spanned() -> SpannedEdn
└── read_string()   -> Edn
```

Before implementation:

- define whether spans are byte offsets, character offsets, or both;
- specify how map keys, collection delimiters, and synthetic values are
  represented;
- ensure ordinary parsing does not pay to retain spans;
- assess whether schema validation adds useful diagnostics beyond the typed
  `FromEdn` conversions. Malli is a design reference, not a dependency
  decision.

## Module organization

Split long DSL modules by grammar or representation boundary, not by line count
alone. Parsing, rendering, EDN conversion, and AST conversion may occupy
separate modules when they have independent responsibilities. Avoid collections
of one-use helper functions.

## Completion criteria

- Every live DSL task has a settled owner and no task is duplicated with docs
  112, 115, or 153.
- Public entry points use one consistent naming scheme.
- Normative specification clauses are covered by a discoverable test or are
  explicitly marked as design-only.
- Tree, streaming, and source-preserving parsing share one language semantics.
