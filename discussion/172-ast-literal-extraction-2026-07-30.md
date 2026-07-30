# 172 — AST literal extraction

Status: **Informational**
Date: 2026-07-30
Relates: [168](168-api-hygiene-2026-07-27.md),
[171](171-aromaticity-inconsistency-policy-2026-07-29.md)

## Scope

Operations in `umol-graph` consume AST values whose fields may be literals, patterns, or
undetermined lattice values. As more chemistry operations move into this layer, repeatedly checking
and extracting literals becomes both a readability cost and a coordination problem. A local
convenience chosen without regard for the operation's contract can silently turn underdetermination
into an error, a skipped candidate, a fallback, or a panic.

This document records the layered policy for literal extraction and the coordination surfaces used
to keep it consistent. It is not an implementation plan. The Rustdoc pass and any concrete ground-view
design are separate work.

## Existing primitive

`AsLit` is the leaf-level abstraction:

```rust
pub trait AsLit {
    type Lit;

    fn as_lit(&self) -> Option<Self::Lit>;
}
```

The trait performs one AST-specific operation: projection to a concrete literal. `Option` supplies
defaults, error conversion, assertions, predicates, and control flow. The former derived methods
`as_lit_ok_or`, `as_lit_ok_or_else`, `as_lit_or`, `as_lit_or_else`, and `as_lit_expect` merely
renamed standard `Option` combinators and were removed.

`as_lit_matches` was also removed. It meant exact equality with an already-literal AST, had no
production users, and conflicted with two established meanings: `Lattice::matches` tests refinement,
while type-specific `matches_value` operations test whether a possibly non-literal constraint admits
a literal. Exact literal equality remains directly expressible as `ast.as_lit() == Some(value)`.

`AsLit::Lit` is returned by value. This is cheap for the common scalar fields, but aggregate
implementations may clone. Any higher-level extraction abstraction must preserve the option to borrow
aggregate data.

## The semantic decision comes first

The same missing literal has different meanings in different operations:

| Meaning | Normal representation |
| --- | --- |
| Candidate cannot be derived yet or operation is inapplicable | `Option` |
| Several fields trigger the same alternative path | tuple `let ... else` |
| Input violates a public operation precondition | domain-specific `Result` |
| Groundness or another invariant was already validated | internal `as_lit().expect` or a checked view |
| Partial state must be retained for later processing | explicit accumulation or derivation state |

Extraction syntax must follow this classification. It must not decide the classification.

## Layered policy

### One-off extraction

Use ordinary Rust control flow:

- `?` in an `Option`-returning derivation;
- tuple `let (Some(...), ...) = (...) else { ... };` for a shared fallback;
- `as_lit().ok_or_else` at a typed error boundary;
- iterator collection into `Option<C>` for homogeneous collections.

Direct enum matching remains appropriate when it is clearer than `AsLit`. A single site does not
justify a helper, wrapper, or projection type.

### Operation-input projections

When the same coherent subset of literal fields is extracted repeatedly, define a named input type
for that operation family. Its constructor performs the extraction once and returns `Option` or a
domain-specific `Result`, according to the operation's contract.

An operation-input projection:

- expresses only the fields the operation requires;
- keeps unrelated AST fields free to remain non-ground;
- provides concrete scalar values to the algorithm;
- borrows aggregate values where practical;
- may include derived values or consistency checks that belong to candidate preparation.

It must not become a generic `PartialGround*` family parameterized over every possible field subset.
Names and public visibility require the ordinary API design review.

### Ground views

A fully ground view is appropriate only where complete groundness is the real precondition. Fingerprint
calculation is the clearest current example: the operation validates groundness and its inner loops
then repeatedly assert `"ground atom"`.

A future ground view should:

- validate once during construction;
- borrow the original AST rather than clone it merely to change the access surface;
- return concrete literal values from entity views;
- make downstream literal assertions unnecessary;
- prevent unchecked construction.

Ground views do not replace operation-input projections. Stereo, aromaticity, and resolution commonly
require only selected literal fields and must continue to operate on partially determined molecules.

### Macros and closure blocks

An immediately invoked closure can create a local `Option` or `Result` scope for `?`, but it is useful
only when it is clearer than `let ... else`.

A macro can remove syntax but cannot determine whether a non-literal value means absence, error,
fallback, skip, or retained partial state. Macros also obscure owned extraction and cloning. No
declarative macro, proc macro, custom carrier, or generic partial-ground abstraction should be added
until repeated call sites show that the remaining duplication is purely syntactic.

Custom implementations of `Try` and `try` blocks are not available on stable Rust at the time of this
decision. Their eventual stabilization would improve syntax but would not remove the semantic
classification above.

## Coordination hierarchy

The policy is exposed through several surfaces because no single surface reaches agents, contributors,
and API consumers reliably:

1. `.agents/skills/ast-literal-extraction/SKILL.md` and its `.claude` counterpart contain the
   normative implementation procedure and trigger on relevant work.
2. `CLAUDE.md` routes AST-consuming implementation work to the skill.
3. This discussion document preserves the design rationale, rejected universal mechanisms, and
   boundaries between the layers.
4. A later Rustdoc pass will document the leaf-level contract at `AsLit` and the checked contracts of
   any projection or ground-view types.
5. Concrete types should enforce validated boundaries wherever possible; documentation is not a
   substitute for a checked constructor.

Memory is not an authoritative policy surface. It is agent-specific, may be stale, and is not visible
to human contributors.

Routine extraction work should load the skill but need not reread this document. Work that introduces
a new abstraction or changes the policy must consult this document and update both the decision record
and the skill.

## Review questions

For new or changed AST-consuming code:

1. What does a non-literal value mean for this operation?
2. Does the public return type express that meaning?
3. Is extraction local, repeated for an operation-specific subset, or guarded by complete groundness?
4. Is `as_lit().expect` dominated by an actual validation boundary?
5. Does owned extraction clone aggregate data?
6. Is a proposed helper removing semantic duplication or only hiding control flow?
