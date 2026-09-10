# 115 · Variable facility (prospecting, deferred)

Status: Proposed
Date: 2026-06-16
Relates: [097](097-joint-domain-design-2026-05-23.md),
[098](098-bind-scope-2026-05-23.md),
[113](113-ast-canonical-equality-and-lattice-2026-06-14.md),
[114](114-atom-bond-interning-2026-06-16.md),
[223](223-pattern-semantics-spike-2026-09-05.md)

2026-09-05: The [nested-query scope decision](#nested-query-scope-decision-2026-09-05)
settles lexical visibility and an initial evaluation restriction. It supersedes the earlier
suggestion that scope may simply follow the region a solver walks. The representation options
and historical implementation observations below remain prospecting; no implementation is committed.

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

## Anonymous bounds vs named variables (the `+` problem)

113's bound sugars introduce variables *implicitly*: `+` (≥1) lowers to
`var_at_least("r", 1)` = `Predicate(Rel(Var("r"), Ge, Lit(1)))`, because `ValueAst` has no
bare field lower-bound — a relation needs a term subject, and the only subject available is
a named variable. That conflates two things this facility must keep apart:

- an **anonymous, field-local bound** — "this count is ≥1," no name, no reuse, each
  occurrence independent;
- a **molecule-scoped named variable** — `?n`, same name = same variable, *deliberately*
  correlating (how `(?lp,?u) :: {…}` ties fields together).

Because `+` uses a *fixed* name (`"r"`), once this facility unifies same-named variables,
multiple bounds correlate spuriously — `#R(5)+#R(6)+` would force `M(5) == M(6)`. **No
correctness bug today** (113 has no unification; `var_at_least` just matches
`value ≥ threshold` per occurrence), which is why 113 ships `+` as `var_at_least("r", 1)`
unchanged. The facility must resolve it on arrival. Options:

1. **A `Self`/`Field` term** — a distinguished `ValueTerm` meaning "the value being
   constrained." `+` → `Rel(Field, Ge, 1)`: non-variable (nothing to unify), structurally
   canonical (identical across occurrences), composes by construction. Cleanly splits
   `Field` (anonymous/local) from `Var` (named/scoped). Cost: the term + field-relative
   handling in `matches`/`meet`/canonicalize (bound intersection). Fixes *every* bound
   sugar, not just `+`.
2. **A bound/range `ValueAst` form** (`AtLeast(n)` / interval) — `+` → `AtLeast(1)`, no
   predicate, no variable; reverses 113's "bounds via `Predicate`" choice; new lattice.
3. **Reserved anonymous name + alpha-normalization** — `+` mints a reserved,
   non-user-writable name the facility excludes from unification; canonicalize alpha-renames
   gensyms. Keeps the predicate encoding but needs the local/scoped flag *and*
   alpha-equivalence.

All three need the same core: **distinguish anonymous local bounds from molecule-scoped
variables, treat the anonymous ones as non-correlating**, keep them structurally canonical.
(1) is the most contained principled answer. Until then, 113 uses the fixed name `"r"` —
the one place to revisit when this facility lands.

## Nested-query scope decision (2026-09-05)

The SMARTS/SMIRKS exploration in [doc 223](223-pattern-semantics-spike-2026-09-05.md)
adopts lexical scope with correlated existential queries for the variable facility. This is
a semantic decision, not a new SMARTS spelling or a choice of variable-environment API.

### Precedents

Lexical binding gives each variable an owning scope: consistently renaming a local binding
and its references, without capture, preserves meaning. Nested existential queries supply
the corresponding visibility and witness rules:

- [SQL correlated EXISTS](https://www.postgresql.org/docs/current/functions-subquery.html#FUNCTIONS-SUBQUERY-EXISTS)
  permits references to enclosing values, which act as constants during an individual
  subquery evaluation. Its result reports existence rather than exporting a witness.
- [SPARQL EXISTS and NOT EXISTS](https://www.w3.org/TR/sparql11-query/#negation)
  similarly test a nested graph pattern using enclosing bindings without generating
  additional bindings for the surrounding query.
- [Soufflé's Datalog grounding rules](https://souffle-lang.github.io/rules) require
  variables to obtain values from positive body predicates and require grounded arguments
  for arithmetic and string functors. This supplies a precedent for restricting evaluation
  without requiring arbitrary equation solving.

These precedents supply reusable rules; adopting them does not require SQL null semantics,
the full SPARQL algebra, or a recursive Datalog evaluator.

### Agreed initial rules

1. Every variable has a type and an owning scope. References resolve structurally,
   independently of matching order. Initially, shadowing an enclosing variable name is
   rejected. Equal spellings in independent sibling local scopes do not establish identity.
2. Nested predicates may read enclosing bindings but cannot reassign them. An inner
   condition may reject an outer candidate, so shared values still impose correlations.
3. Variables introduced within a nested predicate are existential and local. Neither their
   values nor their entity witnesses escape into the surrounding pattern or sibling
   predicates. Locality does not impose disjointness: independently chosen witnesses may
   overlap, subject to each nested pattern's own structural constraints.
4. Negation means that no satisfying local witness exists for the supplied enclosing
   bindings. It exports no bindings and cannot obtain a value from a failed match.
5. The initial evaluator requires referenced values to be supplied by parameters or positive
   structural matches. Conditions compare supplied values; they need not solve equations
   to discover them. This is an evaluation capability restriction, not a change to the
   denotation of symbolic patterns or set-valued results.

For example, an outer match supplies an atom and its charge q. A nested predicate asks
whether that atom has a neighbor with charge q. The neighbor is a local witness; q is an
enclosing input. Negating the predicate asks whether no such neighbor exists. Ordinary
rooted recursive SMARTS already has the analogous fixed input in its root atom; named
attribute variables extend that correlation mechanism.

Evaluation can use an environment of supplied values and nested matching with fixed inputs.
These rules do not require mutable bindings, general unification, or fixed-point evaluation.
Binding syntax, environment representation, and richer symbolic evaluation remain deferred.
Explicit witness export or richer solving can be considered later without changing the
meaning of patterns admitted by this initial fragment.
