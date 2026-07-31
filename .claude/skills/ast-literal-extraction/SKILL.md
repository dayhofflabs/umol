---
name: ast-literal-extraction
description: MANDATORY — load and apply before creating or editing code in umol-graph or higher-level crates that extracts literal values from AST types, calls AsLit methods, adds a ground view, or introduces an operation-specific literal input type. Also apply when reviewing repeated as_lit calls, literal extraction followed by expect, non-literal failure handling, or a proposed macro/helper for AST extraction. Classifies what non-literal values mean and selects the corresponding Rust control flow without hiding operation semantics.
---

# AST literal extraction

Treat `AsLit` as the leaf-level primitive. Decide what a non-literal value means for the operation
before selecting extraction syntax.

## Procedure

1. Read the operation and its public return contract.
2. Classify a non-literal value as one of:
   - normal absence or an inapplicable candidate;
   - a shared fallback path;
   - a domain error at a fallible boundary;
   - impossible after a validated precondition;
   - state that must be retained or accumulated rather than returned immediately.
3. Select the narrowest pattern below.
4. Check whether extraction clones an aggregate `AsLit::Lit`.
5. If the same literal subset recurs, stop duplicating extraction and consider a named
   operation-input type.
6. If complete groundness is the actual operation precondition, consider a validated borrowed
   ground view.

## Pattern selection

| Contract | Pattern |
| --- | --- |
| Function already returns `Option` | Use `as_lit()?`. |
| Several independent values share one fallback | Use tuple `let (Some(...), ...) = (...) else { ... };`. |
| Public operation rejects non-literal input | Use `as_lit().ok_or_else` with the operation's domain error. |
| Validation already established the invariant | Use `as_lit().expect` only within the dominated internal code. |
| Homogeneous collection must be literal | Collect an iterator of `Option<T>` into `Option<Vec<T>>` or the required collection. |
| Repeated stable subset of fields | Introduce a named operation-input projection with `Option` or `Result` construction. |
| Entire structure must be ground | Validate once and operate through a borrowed ground view returning concrete values. |

Use ordinary pattern matching when it is clearer than a trait method. Do not turn a single extraction
site into a helper or projection type.

## Operation-input projections

Use a projection only when a coherent field set recurs within an operation family. It should:

- express that operation's actual precondition, not complete AST groundness;
- extract and validate once at the algorithm boundary;
- carry concrete scalar values and borrow aggregate values where practical;
- have `Option` or a domain-specific `Result` construction contract;
- avoid a generic family of every possible partially ground entity.

Do not invent a public name without the normal API naming review.

## Ground views

A ground view is evidence of a checked precondition, not an unchecked marker:

- validate during construction;
- borrow the original AST rather than clone it merely to change the access API;
- expose concrete literal accessors;
- remove downstream `expect` calls where the view already proves groundness.

Do not require complete groundness for an operation that needs only a subset of literal fields.

## Restrictions

- Do not use `as_lit().expect` to enforce an unchecked caller precondition.
- Do not add `AsLit` methods that duplicate `Option` combinators.
- Do not add increasingly specialized `AsLit` helpers to avoid modeling a repeated operation input.
- Do not hide `Option`, fallback, error, skip, or accumulation semantics in a macro.
- Do not introduce a proc macro, declarative extraction macro, custom carrier, or generic partial-ground
  abstraction without a repo-wide design review.
- Do not assume `AsLit` is free: its associated literal is returned by value and aggregate
  implementations may clone.

If work proposes a new abstraction or changes this policy, read
[`discussion/172-ast-literal-extraction-2026-07-30.md`](../../../discussion/172-ast-literal-extraction-2026-07-30.md)
before proceeding. Routine call-site work does not require rereading the discussion history.
