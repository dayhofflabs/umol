# 115 · Variable facility (prospecting, deferred)

Status: Deferred (forward-looking; not built in the 113 restructure)
Date: 2026-06-16
Relates: 097 (JointDomain), 098 (bind scope), 113 (AST restructure), 114 (interning)

## Purpose

Decide what — if anything — must be done **now** to let equation solving / constraint
satisfaction be added cleanly later (same calculus as interning, 114): do the minimum
unless deferring makes it more work. This is prospecting, nothing committed.

## A variable has three attributes

- **Type** — the domain it ranges over: `Element`, a `Value` slot (charge/H/lone-pairs/
  order/…, all i64), `Coset` (per stereo kind), maybe `Isotope`/`Spin`. Mirrors the leaf
  `AsLit` types. *Today implicit:* a `Var` in `ValueTerm`/`ElementAst`/`StereoTerm` is
  typed by which leaf holds it (a variable ranges over exactly one domain).
- **Identity** — a name (`String` today; could be an interned `Sym`).
- **Scope** — the container level whose environment owns the binding:
  **atom < molecule < reaction-pair < network**. *Today implicit:* no environment exists;
  same-name occurrences co-refer within whatever region a solver walks.

## Where each piece lives, by level

| Piece | atom | molecule | reaction | network |
|---|---|---|---|---|
| variable declarations / domains | — | `VarEnv`? | `VarEnv`? | `VarEnv`? |
| correlations (joint / relational) | `JointDomainAst` (today) | molecule joint constraint | cross-side conservation / mapping | parametrization |
| solver | `saturate` (today) | molecule engine | reaction engine | network engine |

A correlation lives at the **lowest level containing all its variables**; an env lives at
the level across which a name is meant to co-refer.

## Two reasonable shapes

**A — occurrence-based (current trajectory, minimal).** Vars stay bare names in leaf
ASTs; no env; type inferred from location; scope = the region the solver collects
occurrences over. Higher levels add joint-constraint variants + a per-level solver.
- Pros: nothing to build now; purely additive.
- Cons: type/scope implicit (re-inferred each pass); same-name-different-type collisions
  surface only at solve; heterogeneous correlations (`#E(?el,?u)` — element × spin)
  reconstruct each var's type from its occurrence.

**B — explicit typed, scoped vars.** `Var(VarId)` resolving in a `VarEnv`
(name → `VarType` + level + optional domain) owned by the container; joint/relational
constraints reference `VarId`s.
- Pros: types explicit and checkable; collisions caught; cross-type/cross-level
  correlation and network parametrization clean.
- Cons: an env type + plumbing; constraints carry `VarId`s.

## Retrofit calculus

113 already did the one enabling move: removing atom-level `Bind`/`Ref` made vars
**scope-agnostic bare names** and pushed solving out of the AST. Consequently:

- **Scope / env / higher-level correlations are additive** — a `VarEnv` is a
  container-level layer keyed on the stable var name and consulted by the solver (like the
  intern pool, 114); molecule/reaction/network joint constraints are new variants on
  containers that don't exist yet. None forces an AST change now. **No `atoms_mut`-style
  hazard** — vars are referenced by a stable name a later env can key on.
- **Per-leaf typed vars are correct, not a limitation** — a variable ranges over one
  domain, so nothing multi-level needs to break that.
- **The only borderline "do-now":** the `Var` *representation* itself, which 113 edits
  anyway (`Bind`/`Ref` → `Var`). Enriching it now (`Var(String)` → typed/interned
  `VarId`) is incremental while in there; doing it later is a *contained* change (per-leaf
  `Var` variant + parser + constraint refs) that does **not** ripple into
  `canonicalize`/`equiv` (vars stay opaque). Later cost ≈ now cost → no "later is more
  work" pressure; it's a preference (carry type from day one for heterogeneous
  correlations), not a necessity.

## Read

Same shape as interning: **minimal now suffices.** The facility is additive on top of
113's scope-agnostic vars; no retrofit trap. Envs, multi-level joint constraints, and the
solver are built when the molecule/reaction/network levels are.

## Subsumes (and replaces) JointDomain

`JointDomain` (doc 097) is **removed** (113, 2026-06-16): it correlated atom *fields*
(`JointVar` tags), not variables — a parallel premature system that only expressed the
bare-var case and forced a privileged-bare-var position, and whose serde "blocker" (name ↔
field resolution) was an artifact of the field-tag representation. This facility replaces
it: a relation/table constraint over **named variables** stores the names directly (no
field resolution, no serde blocker), expresses the general case (vars in arithmetic
expressions, cross-atom), and the Fe²⁺ high/low-spin split is just
`(?lp, ?u) :: {(3,0),(1,4)}` over two vars. So the table-shaped constraint comes back as
part of this facility, var-based — not as the removed `JointDomain`.
