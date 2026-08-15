# 128 — Substructure matching of derived topological predicates

Status: Superseded — the matcher on the views' keyed readings
([194](194-constraint-assertion-semantics-2026-08-10.md) S1) evaluates derived predicates
during matching; ring-scope queries (`#R(6)`) verified end to end.
Date: 2026-06-23

## Problem

Substructure queries that constrain a *derived topological* predicate — ring
membership (`#R`, `[R]`, ring-bond `@`), degree (`#D`), connectivity (`#X`), ring
connectivity (`#x`), total hydrogens (`#H`) — do not match any concrete target,
even when the target plainly satisfies them. Surfaced while building the RDKit
`PatternFingerprint` replica (doc 126, slice 5 direction B): the four ring-junction
templates (`[R](@[R])(@[R])~[R]~…`) match nothing, so fused-ring fingerprints are
incomplete. Naphthalene pins it — umol reproduces 117 of RDKit's 126 bits; the 9
missing are exactly the junction features. Nothing is wrong, only absent.

Per the DSL spec §6.1 these predicates are **derived**: "a topological query
evaluated against the target graph once an embedding is proposed." They filter
matches and do not affect grounding. The expectation is that `[R]` matches any
target atom the graph shows to be in a ring.

## Two implementation gaps

`MoleculeAst::substructure_matches` builds host match-targets via
`host_match_targets` → `derive_constraints(true)`, folding derived facts into each
host atom/bond, then matches the pattern field-wise.

**Gap 1 — derived ring/degree facts are never folded in.** `AtomView::
derive_constraints` emits localized valence, dative donated/accepted pairs, aromatic
valence, multicenter valence, and tetrahedral stereo. It does **not** emit ring
membership, degree, total degree, ring degree, ring valence, or total hydrogens.
`BondView::derive_constraints` emits only cis/trans stereo — not ring membership.
The matcher for `AtomConstraint::RingMembership` reads
`target.ring_membership_value(scope)`, which looks up a **stored** constraint on the
target and returns `None` when absent, falling back to `Undetermined`. A pattern's
`RingMembership(All, Lit(k))` against `Undetermined` is false — only an
`Undetermined` pattern passes. So every non-trivial ring/degree pattern fails. The
data exists: `AtomView::ring_membership(scope)`, `degree()`, `total_hydrogens()`,
`ring_degree()`, and `BondView::ring_membership(scope)` all compute the value from
the perceived graph; it is simply never folded into the match target.

**Gap 2 — the `+` ("at least one") form is a predicate the matcher does not
evaluate.** Spec §7.3: `#R+` (and `[R]`) is sugar for `?r >= 1`, a `var_at_least`
predicate. Spec §6.2: a bool-expr pattern matches a literal target iff the
expression holds on the literal. But `ValueAst::matches` evaluates a
predicate-vs-literal pair through `meet`, and the meet of a free-variable predicate
with `Lit(n)` is `None`, so the match is false regardless of `n`. Even once Gap 1
supplies the host's concrete ring count, `#R+` would not match it. This affects
every `>=` / bool-expr query slot (e.g. `#D(?d >= 2)`), not only ring membership.

Both gaps must close for `[R]` / `@` / `#R+` to match per spec.

## Design axes

### A. Where the derived fact is supplied (Gap 1)

1. **Fold into `derive_constraints`** — consistent with how valence/aromatic facts
   already reach the matcher. `AtomView`/`BondView::derive_constraints` additionally
   emit the topological predicates computed from the view; the matcher is unchanged.
   Open question: ring **size** scopes — `#R(6)` asks `RingScope::Size(6)`. Emit
   `All` plus one entry per ring size present in the molecule. A query for a size
   that is absent then finds no stored entry (`Undetermined`), which correctly
   matches `#R(s)*` but not `#R(s)!` ("no ring of size s") — a gap for negative
   sized queries. The `All` scope (what `[R]` / `@` / the junctions need) has no
   such issue.

2. **Compute in the matcher on demand** — the constraint matcher evaluates ring
   membership/degree against the target *view* at match time rather than a folded
   constraint. Handles arbitrary scopes uniformly, so negative and sized queries
   work. Cost: the matcher must hold the target view (graph access), not just the
   target AST constraint set — a larger change to the matching signature than A1.

### B. Predicate evaluation (Gap 2)

1. **Evaluate bool-expr / `var_at_least` against concrete targets** in
   `ValueAst::matches`, per spec §6.2 (the expression holds on the literal).
   General — fixes every query predicate, not only ring `+`. Scope: the
   value-matching core, used everywhere.

2. **Avoid the predicate form for these templates** — encode `[R]` / `@` as an
   exact/set form the current matcher already handles (e.g. a finite `LitSet` of
   admissible ring counts). Not faithful (an atom in more rings than enumerated is
   missed) and specific to this one caller; recorded only for completeness.

## Scope of the change

- **Minimal (unblock junctions):** `All`-scope ring membership on atoms and bonds
  (Gap 1 via A1 for the `All` scope) + predicate evaluation for `>= n` (Gap 2 via
  B1). Enables `[R]`, `@`, `#R+`.
- **Spec parity:** the full derived set — ring membership (all scopes), degree,
  total degree, ring degree, ring valence, total hydrogens, connectivity, ring
  connectivity — folded/computed, plus general bool-expr evaluation. Larger, but
  closes the derived-predicate conformance gap once.

## Decisions needed

1. Gap 1 mechanism: fold into `derive_constraints` (A1) vs compute-in-matcher (A2).
2. Gap 2: evaluate predicates in `ValueAst::matches` per spec (B1) — confirm this is
   the intended fix and not deferred.
3. Scope: minimal (`All`-scope ring + `>=`) vs full derived-predicate parity.
4. Sized-ring negative queries (`#R(s)!`) under A1 — accept the gap, or require A2.

## Downstream

Closing this unblocks doc 126 slice 5 B junction templates (and naphthalene
bit-exactness), and makes ring/degree/`@` substructure queries usable generally
(BRIDGIT reactive sites, a future SMARTS-level query layer).
