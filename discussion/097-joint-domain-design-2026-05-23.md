# JointDomain constraint design

Status: AST + Lattice + saturate **Done** (doc 96 substeps 2f–2j). DSL surface and EDN roundtrip **Deferred** to a later pass — see doc 96 step 2 status table and doc 98 (bind scope).

## Why this exists

The valence resolver needs to represent "this atom is one of (charge=0, h=1) or (charge=1, h=2), pick one." Per-field domains alone can't say that: setting `charge ∈ {0,1}` and `implicit_hydrogens ∈ {1,2}` independently admits the unwanted (0,2) and (1,1) combinations. The constraint is a *joint* domain over a tuple of fields — a finite set of admissible tuples that cuts the cross-product down to the legal cases.

Concretely from the iron(II) example: `Fe²⁺` can be (lone-pairs=3, unpaired=0) or (lone-pairs=1, unpaired=4) but not the off-diagonal combinations. JointDomain encodes that.

This is a relational / CSP-style constraint. The umol lattice machinery handles per-variable narrowing natively; JointDomain adds the inter-variable layer.

## DSL notation (deferred)

```
Fe#c+2#n?n#u?u#E(?n,?u) :: {(3,0), (1,4)}
```

- `?n`, `?u` — named binds declared on per-field predicates (`#n?n` declares `?n` on lone-pairs).
- `#E(vars) :: {tuples}` — the joint-domain predicate. Set semantics (unordered tuple list); `::` matches the existing bind-domain membership glyph (`?h :: {1,2,3}`).
- Inner tuples `(v1, v2)`. Outer set `{(...), (...)}`. Comma at both levels; `()` brackets each tuple so the outer comma is unambiguous.
- Tuple values are always literals (no expressions, no nested sets) — neither delimiter can be confused with `value-expr`.

The parser + formatter are not implemented yet. Their work is blocked on bind-scope design (doc 98); see "Deferred work" below.

## Type shape (done)

```rust
pub enum JointDomainAst {
    Undetermined,                            // top of the lattice
    Domain(DomainState),                     // proper constraint
}

pub struct DomainState {
    vars: Vec<JointVar>,                     // private fields
    tuples: Vec<Vec<JointValue>>,
}
```

`Undetermined` is the lattice top, matching the convention every other lattice type in `umol-ast` follows. `DomainState` has module-private fields — all construction goes through `JointDomainAst::from_ints` (or sibling constructors for non-numeric types when those land), so stored invariants hold by construction. External code reads via `JointDomainAst::vars()` and `tuples()` accessors, or via the `Lattice` trait methods.

```rust
pub enum JointVar {
    Charge, ImplicitHydrogens, LonePairs,
    UnpairedElectrons, Multiplicity,
    Valence, DonatedPairs, AcceptedPairs,
}

pub enum JointValue {
    Int(i64),
}
```

Both `#[non_exhaustive]`. `JointVar` variants spell the atom-AST field name in full (no `ImplicitH` shortenings). Spin's `unpaired` field surfaces as `UnpairedElectrons` for read clarity. `JointValue` is currently integer-only; will extend to `Element`/`Isotope` when those become joint vars.

Files: `umol-ast/src/ast/constraint/joint_domain.rs`.

## Invariants (done)

Constructor `JointDomainAst::from_ints(vars, tuples)` rejects:
- `vars.len() < 1` (zero is degenerate)
- `tuples.len() < 1` (empty is bottom; signaled via `Lattice::meet -> None`, not stored)
- `tuples[i].len() != vars.len()` for any `i` (arity mismatch)
- duplicate vars

Canonicalization on success:
- Sort `vars` (and permute each tuple to match)
- Sort `tuples` lexicographically
- Dedup `tuples`

A `Domain { tuples: [t] }` (single tuple) is a *ground* state, not a transient invariant violation. A `Domain { vars: [single] }` is redundant with a per-field constraint on `single` but harmless — saturate may normalize it.

## Lattice behavior (done)

Every operation is total against `Undetermined`:

- `is_undetermined`: matches `Undetermined`
- `is_ground`: matches `Domain` with exactly one tuple
- `meet(Undetermined, x) = Some(x)`
- `meet(Domain, Domain)` is the **relational meet** (natural join): cartesian product on disjoint vars, equijoin on shared, intersection on identical. Returns `None` only when the joined tuple set is empty (genuine contradiction).
- `join(Undetermined, _) = Undetermined`
- `join(Domain, Domain)` projects both to shared vars, unions, dedups, wraps as `Domain`. If shared vars are empty, returns `Undetermined`.
- `matches(Undetermined, _) = true`
- `Domain` pattern matches `Domain` target iff `pattern.vars ⊆ target.vars` and every target tuple projected to `pattern.vars` is in `pattern.tuples`.

## Simplify vs saturate (done)

Two responsibilities, two operations:

- **`JointDomainAst::simplify`** — infallible per-type canonical-form normalization. Re-sorts and dedups `tuples`. Idempotent. Safety net for code paths that build a `Domain` via internal constructors (e.g., the relational meet path).
- **`Lattice::saturate`** — fallible cross-field propagation. For `AtomAst`, this is `saturate_atom`: walks every JointDomain in `atom.constraints`, projects its tuples against current field values (forward propagation only), and either prunes the JD, narrows fields to a surviving single-tuple's literals and drops the JD, or returns `Err(Contradiction)`. Wired into `AtomAst::meet` via `#[lattice(saturate = "saturate_atom")]`.

Forward-only is the current scope. Backward propagation (narrowing field domains from JD tuples) is a possible extension.

## Settled decisions

- **Set semantics (unordered) for the tuple list.** Meet is intersection — commutative, associative, idempotent, monotone. Constructor + meet canonicalize.
- **DSL syntax**: `#E(?v1,…,?vn) :: {(l1,…), …}`, `::` glyph, comma at both levels.
- **Naming**: `JointDomainAst`, `JointVar`, `JointValue`. Variants follow atom AST field names. `UnpairedElectrons` is the explicit rename of `Unpaired`.
- **Encapsulation**: `DomainState` has module-private fields. All construction routes through `from_ints`.
- **JointVar separation**: `JointVar` is constraint-internal; the bind-scope work (doc 98) introduces `AtomBindTarget` separately. They may merge later.
- **No CSP-vs-lattice framing battle**: lattice and CSP vocabularies are equivalent for the algorithm. Code stays lattice; references to CSP literature live in docs where helpful.

## Alternatives considered

- **Ordered (ranked) tuple list** — first-match preference for solver ranking. Breaks commutativity of meet; would force the lattice into asymmetric behavior. Rejected; if ranking is ever needed, attach it at the search-driver level over candidate sets, not inside the constraint type.
- **Weighted CSP (soft constraints)** — each tuple carries a cost; meet adds costs; narrowing shrinks to min-cost tuples; solving becomes optimization. Principled and well-studied. Not adopted in this pass because it changes propagator/solver semantics across the codebase. Open option for a future pass if quantitative ranking has to compose under meet.

## Deferred work

These are blocked on adjacent design decisions and not in scope for the current pass.

### DSL parser + formatter for `#E` (doc 96 step 2k)

The compact-DSL parser/formatter for `#E` requires resolving bind names like `?n` to a `JointVar` (parser direction) and back (formatter direction). That resolution depends on having a bind scope per atom — see doc 98.

Until doc 98 lands, `umol-ast/src/dsl/atom.rs::fmt_constraint` panics via `todo!()` on the `JointDomain` arm. Any test or call site that needs to format a JD-bearing atom must either:
- avoid creating JD on atoms that get formatted, or
- mark the failing path with `#[should_panic]` and revisit when 2k lands.

Programmatic construction via `JointDomainAst::from_ints` works fine and does not exercise the formatter.

### EDN roundtrip (doc 96 step 2l)

Same dependency chain: EDN serialization needs to express the JD, which needs bind-name representation, which needs the scope. Deferred to the same pass as 2k.

### Bind-name scoping

Lifted out into doc 98. Originally an open question in this doc's predecessor; the work to give atoms an explicit variable scope is large enough to be its own design doc. The JointDomain AST already commits to `JointVar` (not strings), so the scope question is about the DSL surface and validation paths, not the AST shape.

### JointVar extensions

Reserved (not yet implemented): `Element`, `Isotope`, `AromaticValence`, `MulticenterValence`, `HapticValence`. Each requires a corresponding `JointValue` variant and lattice handling for the projection step. The `#[non_exhaustive]` markers on `JointVar` / `JointValue` make adding them backward-compatible.

Concrete deferred use case: `[Fe; Co] #E(?el,?u) :: {(Fe,4), (Co,3)}` — metalloprotein active sites where metal identity couples to spin state.

### `#E` as one or many constraint kinds

Product-domain (the current `JointDomainAst`) is one relational shape. Others — equality (`?n = f(?u)`), linear (`?n + ?u ≤ k`), single-variable membership — could be sibling constraint variants or could compose through `ValueAst::Expr`. Not pursued in this pass.

### Saturate: backward propagation

Current saturate is forward-only (field state prunes JD tuples). Backward (JD tuples narrow field domains) would give stronger arc consistency. Optional for a future pass; not needed for the resolver's current "fail if narrowing produces non-singleton" gate.

## Implementation status table

| Substep | Description                                                             | Status   |
| ------- | ----------------------------------------------------------------------- | -------- |
| 2f      | `JointDomainAst` type + `from_ints` constructor                         | Done     |
| 2g      | `AtomConstraint::JointDomain` variant wiring + container ops            | Done     |
| 2h      | Hand-rolled `Lattice` impl on `JointDomainAst`                          | Done     |
| 2i      | `#[derive(Lattice)]` proc-macro (applies to `AtomAst`)                  | Done     |
| 2j      | `Lattice::saturate` + `saturate_atom`                                   | Done     |
| 2k      | DSL syntax + parser/formatter for `#E`                                  | Deferred |
| 2l      | EDN serialization roundtrip                                             | Deferred |

## CSP precedent (for reference)

JointDomain + per-variable domains + arc-consistency narrowing + (eventually) labeling search is classical finite-domain CSP. The lattice + `narrow_from` machinery already in place is arc consistency in lattice clothing. Relevant literature:

- **Mackworth's arc consistency** (AC-3, AC-4) — original algorithms.
- **CLP(FD)** — Prolog dialects with `library(clpfd)`; the "name your variables, post relational constraints, solve" surface.
- **MiniZinc** — declarative constraint-modeling language; `table` predicate matches our `#E`.
- **Saraswat's Concurrent Constraint Programming** — lattice-theoretic formulation; maps to umol's `Lattice` + `narrow_from`.
- **Apt, *Principles of Constraint Programming*** (Cambridge, 2003).
- **SMARTS/SMIRKS** matching, **RDKit**'s VF2-derived substructure matcher — chemistry-internal precedent.

Vocabulary coexistence is fine: lattice and CSP are equivalent views of the algorithm. Switching to an external CSP solver (`copper`, `pumpkin`, MiniZinc port) is a possible future move; not on the table right now because the inference problems here stay atom-scale.

## CSP-to-umol mapping

| CSP concept                     | Current umol equivalent                                                                   |
| ------------------------------- | ----------------------------------------------------------------------------------------- |
| variable                        | AST field (e.g., `atom.charge`, `atom.implicit_hydrogens`)                                |
| domain                          | `Set` / `Lit` / `Undetermined` on field AST                                               |
| propagator                      | `Lattice::narrow_from`                                                                    |
| arc consistency                 | meet-driven narrowing in resolver passes                                                  |
| named variable                  | `ElementAst::Bind` / `Ref`, `ValueAst::Bind` / `Ref`, `IsotopeAst::Bind` / `Ref`          |
| cross-variable constraint       | `JointDomainAst` (per-atom, value-typed slots only today)                                 |
| labeling / search               | not yet — resolver enumerates candidates but has no general search driver                 |
| fixed-point engine              | `Lattice::saturate` loops per-entity to fixpoint; no molecule-level engine yet            |

## Related docs

- **96**: valence resolution plan (top-level). Step 2 = JointDomain (this doc).
- **98**: bind scope. Blocks 2k/2l in doc 96 and the DSL surface here.
- **80**: unified constraint AST (older context).
- **83**: constraint unification architecture (older context).
