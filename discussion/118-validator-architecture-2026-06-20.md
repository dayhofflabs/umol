# Validator architecture — validation tiers and model-carrying wiring

Status: **Active / design + impl plan.** Pseudo code agreed; open names and one layering
question must close before coding. No code authorized yet.
Date: 2026-06-20.
Trigger: wiring the stereo validator (doc 104 Phase C / doc 111 C4c) surfaced that the composite
`Validator` excludes the model-carrying validators, nothing invokes validation in production, and
there is no tier-3 valence conformance validator (the valence counterpart of the aromaticity/stereo
validators). This doc formulates the validator architecture rather than extending doc 104.

Note: the "three tiers" in doc 065 are an **error-type** architecture (sub-concern → module dispatch →
cross-module box). The tiers here are a **validation** taxonomy — unrelated; do not conflate.

## Validation tiers

The boundary between tiers is **representability**: *in which type can the broken state exist?* That
dictates where the check runs; phase (raise vs construction) is a consequence, not the criterion.

| tier | concept | violation representable in | enforced at | a validator? |
|---|---|---|---|---|
| 0 | reference resolution | `MoleculeDsl` only (refs `AtomRefDsl`, aliases in `Metadata`) | raise (DSL→AST) | **no** |
| 1 | integrity | `MoleculeAst` | construction (`MoleculeAst::new`) + standalone | yes |
| 2 | invariants | `MoleculeAst` | standalone | yes |
| 3 | conformance | `MoleculeAst` + chemistry model | standalone | yes |

tier-0 is **not** a peer validator: refs are erased on raise (the AST carries no metadata), so a
dangling ref is unrepresentable in `MoleculeAst` and cannot recur — it is a translation precondition
owned by raise, not an AST check. Roundtripping is preserved at the `MoleculeDsl` level (AST +
`Metadata`), independent of this. The tier-0/tier-1 line *is* the `MoleculeDsl`↔`MoleculeAst` boundary,
already enforced by the type system.

tier-1 through tier-3 all operate on a constructed `MoleculeAst`; they differ by what knowledge a
violation requires: integrity needs only the AST, invariants need physics, conformance needs a chemistry
model.

## Tier membership and ownership

| tier | members | model-carrying | ownership |
|---|---|---|---|
| 1 integrity | `EntityStructureValidator`, `ConstraintValidator` (currently a stub) | no | own |
| 2 invariants | `ValenceInvariantsValidator`, `SpinInvariantsValidator` | no | own |
| 3 conformance | `ValenceConformanceValidator` (new), `AromaticityConformanceValidator`, `StereoConformanceValidator` | yes | valence **borrows**; aromaticity/stereo **own** |

Ownership mirrors the resolvers and is driven by data size: the valence table and atom-typing registry
are large, so the conformance validator borrows `&'a ValenceModel` (like `ValenceResolver<'a>`); the
aromaticity and stereo models are small and owned (like their resolvers). The composite therefore gains
a lifetime, `Validator<'a>`, mirroring `Resolver<'a>`.

## Shared model-carrying structs

A shared struct is justified only where the resolver and validator perform the **same** derivation.
That holds for valence and aromaticity but **not** stereo:

| domain | shared piece | resolver use | validator use |
|---|---|---|---|
| valence | `CountsValence<'a>` / `AtomTypingValence<'a>` | `resolve_atom` (assign) | add read-only `conforms_atom` (check) |
| aromaticity | `AromaticityPerception` | `find_systems` + `add_systems` | `find_systems` + compare to stored |
| stereo | none needed | transcribes the coset from the `#T`/`#C` marker — **computes no symmetry** | derives `graph_symmetry` / `stereo_*_symmetry` (umol-ast) + checks stored |

Consequences:
- **`AromaticityPerception` keeps its name.** Both ops call `find_systems` — they (re-)perceive aromatic
  systems; "perception" is accurate.
- **No stereo counterpart struct.** The stereo resolver derives nothing (it copies the marker coset
  verbatim), so there is no computation shared with the validator to factor out. The validator calls
  umol-ast's symmetry methods directly — keeping doc 111 C4's original decision. "Perception" would also
  be wrong here: in stereo it is reserved for the deferred geometry→`#T`/`#C` path.
- Valence's `CountsValence`/`AtomTypingValence` gain a read-only `conforms_atom`; the validator reuses
  the engines' derive/lookup (no duplicated table/registry logic).

## Composite `Validator` API

Separate methods, one per tier; an all-tiers `validate`. No exposed tier enum (tiers are internal
jargon), no intermediate sub-composite layer.

```
struct Validator<'a> {
    // tier 1 — integrity
    entity_structure: EntityStructureValidator,
    constraint:       ConstraintValidator,
    // tier 2 — invariants
    valence_invariants: ValenceInvariantsValidator,
    spin_invariants:    SpinInvariantsValidator,
    // tier 3 — conformance
    valence_conformance: ValenceConformanceValidator<'a>,  // borrows &'a ValenceModel; dispatch AtomTyping|Counts
    aromaticity:         AromaticityConformanceValidator,    // owns AromaticityModel; wraps AromaticityPerception
    stereo:              StereoConformanceValidator,         // owns StereoModel; calls umol-ast symmetry directly
}
impl<'a> Validator<'a> {
    fn new(model: &'a ChemistryModel) -> Self
    fn validate_integrity   (&self, ast) -> Result<Solution<…>>  // entity_structure, constraint
    fn validate_invariants  (&self, ast) -> Result<Solution<…>>  // valence_invariants, spin_invariants
    fn validate_conformance (&self, ast) -> Result<Solution<…>>  // valence_conformance, aromaticity, stereo (stereo last)
    fn validate             (&self, ast) -> Result<Solution<…>>  // integrity → invariants → conformance
    fn validate_atom        (&self, atom) -> Result<Solution<…>> // unchanged: invariants per-atom subset
}
```

`ValidatorContradiction` / `ValidatorError` gain `#[from]` variants for the three tier-3 sub-validators
(aromaticity, stereo, and valence conformance contradictions; aromaticity's error remains
`AromaticityError`). Contradiction/error type names are deferred to the doc 065 umol-graph error survey.

Scope this round is **composer-only**: make `Validator` compose and order all tiers. No new production
run site (no `resolve_and_validate`, harness still drives `Resolver` directly); callers invoke the
tier methods they want.

## tier-3 valence conformance validator

Read-only twin of `ValenceResolver`, dispatching on `ValenceModel`:

```
enum ValenceConformanceValidator<'a> { AtomTyping(AtomTypingValence<'a>), Counts(CountsValence<'a>) }
impl<'a> ValenceConformanceValidator<'a> {
    fn new(model: &'a ValenceModel) -> Self        // mirror ValenceResolver::new
    fn validate(&self, ast) -> Result<Solution<(), …Contradiction>, …Error>  // fold conforms_atom over atoms
}
```

Conformance reuses the engines' derive/lookup — no duplicated table/registry logic:

```
impl CountsValence<'_> {
    fn conforms_atom(&self, atom) -> Solution<(), CountsMismatch> {
        // ground atom must be admitted by the table-derived fields:
        //   entry = table.entry(element.shift(2·accepted_pairs − charge)); some target_covalence ≥ valence;
        //   if aromatic, aromatic valence ∈ aromatic_valences.
        // re-derive via derive_fields; Determined iff atom.matches(derived); Contradictory on NoMatch.
        // non-ground / undetermined element|charge → Underdetermined.
    }
}
impl AtomTypingValence<'_> {
    fn conforms_atom(&self, atom) -> Solution<(), AtomTypingMismatch> {
        // patterns = registry.lookup(element, Some(charge)) (+ element-only fallback);
        // empty → Contradictory(NoType); any pattern admits atom → Determined; else Contradictory(NoMatch).
    }
}
```

## Vars and interning

Current free-string `ValueTerm::Var(String)` is self-contained: nothing to dangle, so tier-1 has no
var check yet (vacuously a member). Interning at raise (doc 114-style) would bifurcate vars along the
same `MoleculeDsl`↔`MoleculeAst` boundary as atom refs:

- name → `VarId` resolution and type-clash unification happen at raise → tier-0 (erased afterward).
- `VarId` ↔ interned-var-table correspondence becomes a real tier-1 invariant (a `from_parts` AST could
  dangle a `VarId`), provided the var table lives in the AST (it must, for a pattern AST to be
  self-describing).

So tier-1's var content sharpens from vacuous to a ref-correspondence check if/when interning lands
(doc 114 / doc 115). Not wired now.

## Settled naming
- tier-1 = **integrity** (`validate_integrity`); members `EntityStructureValidator`, `ConstraintValidator`.
- tier-2 = **invariants** (`validate_invariants`); `ValenceInvariantsValidator`, `SpinInvariantsValidator`
  (rename of `SpinCouplingValidator`).
- tier-3 = **conformance** (`validate_conformance`); `ValenceConformanceValidator` (new),
  `AromaticityConformanceValidator` / `StereoConformanceValidator` (renames of `AromaticityValidator` /
  `StereoValidator`).
- `AromaticityPerception` kept; no stereo counterpart struct.
- Contradiction/error type names deferred to the doc 065 umol-graph error survey.

## Module naming

A module is named for what it exports: a **result-object noun** if it has one, otherwise the **verb**
for the action with the agent type inside (std `fmt::Formatter` pattern; also kills the
`resolver::Resolver` stutter). Renames this work entails:
- `ops/resolver.rs` + `resolver/` → `resolve.rs` + `resolve/`
- `ops/validator.rs` + `validator/` → `validate.rs` + `validate/` (umol-ast tier-1 validators live in a
  umol-ast `validate` module; umol-graph `validate` holds tier-2/3 + the composite)
- `ops/transformer.rs` + `transformer/` → `transform.rs` + `transform/`
- `ops/invariants.rs` → `invariant.rs` (singular)

Unchanged: `parse`/`edit` (already verbs); domain/data modules stay nouns (`model`, `valence`,
`aromaticity`, `value`, `ring`, `symmetry`, `reaction`); `coloring`/`embedding`/`matching` stay nouns —
they export `Coloring`/`Embedding`/`Matching`, so the data-module rule applies (not a gerund exception).

## Layering — tier-1 in umol-ast (resolved)

All three `MoleculeAst` construction paths must enforce tier-1 uniformly: direct construction and DSL
raise (umol-ast), and TableIR raise (umol-io → depends on umol-ast). umol-ast is the common floor, so
`EntityStructureValidator` and `ConstraintValidator` (with their contradiction/error types) **move from
umol-graph to umol-ast**.

This pulls one shared type down. `Solution<T, C>` (currently `umol-graph/src/ops/solution.rs`) has **no
lattice content**: a generic three-valued outcome (`Determined`/`Underdetermined`/`Contradictory`) with
process combinators (`map`, `map_contradiction`, `into_observation`, `into_decisive`), generic over both
payloads. The lattice reading (Determined≈ground, Contradictory≈⊥) is interpretive, applied by callers.
It is the outcome floor for fallible/iterative constraint-satisfaction passes — parallel to `Error` — so
it **relocates to umol-shared** (umol-ast already depends on umol-shared; umol-graph re-imports from
there). The lattice `Contradiction` (`Canonicalize::canonicalize`'s `Err`, a meet/⊥ failure) is genuine
lattice algebra and **stays in umol-ast** — it is one possible `C`, not part of `Solution`.

Construction stays **open**: `from_parts` builds any AST without enforcement (open data type); the raise
paths run tier-1 after building, and tier-1 is callable standalone. The umol-graph composite `Validator`
composes the umol-ast tier-1 validators alongside its tier-2/3 members.

## Implementation plan

Sequenced bottom-up by crate to minimize rework; red periods between phases are acceptable. Paths use the
post-rename names once a rename has happened in an earlier phase.

1. **umol-shared — relocate `Solution`.**
   - `umol-shared/src/solution.rs` (new) ← move `Solution<T, C>` + impls + tests verbatim from
     `umol-graph/src/ops/solution.rs`; declare `pub mod solution;` in `umol-shared/src/lib.rs`. No
     top-level re-export — accessed by submodule path, like `umol_shared::element::Element`.
   - umol-graph: repoint every `crate::ops::solution::Solution` to `umol_shared::solution::Solution`;
     delete `ops/solution.rs`. Mechanical sweep; umol-graph stays green.

2. **umol-ast — tier-1 validators move in.**
   - `umol-ast/src/ast/validate.rs` (new) ← move `EntityStructureValidator`, `ConstraintValidator` and
     their contradiction/error types from `umol-graph/src/ops/validator/{entity,constraint}.rs`; return
     `umol_shared::solution::Solution`. Register `mod validate`; re-export the two validators.
   - DSL raise entry (umol-ast): run tier-1 after building the AST; `from_parts` stays open (no
     enforcement); validators remain callable standalone. [confirm exact raise fn at code time]
   - umol-graph: delete the two moved files from `ops/validator/`; the composite temporarily drops tier-1
     (rewired in step 6).

3. **umol-io — TableIR raise.**
   - TableIR→AST raise path: call the umol-ast tier-1 validators after building, mirroring the DSL path, so
     all construction routes enforce tier-1 uniformly.

4. **umol-graph — renames (mechanical).**
   - modules: `ops/resolver.rs`+`resolver/`→`resolve`; `ops/validator.rs`+`validator/`→`validate`;
     `ops/transformer.rs`+`transformer/`→`transform`; `ops/invariants.rs`→`invariant.rs`. Update `mod`/`use`
     paths and `ops` re-exports.
   - types: `SpinCouplingValidator`→`SpinInvariantsValidator` (+ its contradiction/error);
     `AromaticityValidator`→`AromaticityConformanceValidator`; `StereoValidator`→`StereoConformanceValidator`.

5. **umol-graph — tier-3 valence conformance.**
   - `ops/valence/counts.rs`, `ops/valence/atom_typing.rs`: add read-only `conforms_atom`
     (+ `CountsMismatch` / `AtomTypingMismatch`) reusing `derive_fields` / `registry.lookup`.
   - `ops/validate/valence.rs` (new): `ValenceConformanceValidator<'a>` enum (AtomTyping | Counts),
     `new(&'a ValenceModel)` mirroring `ValenceResolver::new`, `validate` folding `conforms_atom` over atoms.

6. **umol-graph — compose `Validator`.**
   - `ops/validate.rs`: `Validator<'a>` with seven members (tier-1 from umol-ast, tier-2 invariants, tier-3
     conformance); `new(&'a ChemistryModel)`; `validate_integrity` / `validate_invariants` /
     `validate_conformance` / `validate`; `validate_atom` unchanged. Add `#[from]` variants to
     `ValidatorContradiction` / `ValidatorError` for the three tier-3 sub-validators.

7. **tests.**
   - new: `conforms_atom` (counts + atom-typing; determined / contradictory / underdetermined);
     `ValenceConformanceValidator`; composite per-tier (`validate_integrity` / `_invariants` /
     `_conformance`) and full `validate`; tier-1 enforced on the DSL and TableIR raise paths while
     `from_parts` stays open.
   - update: `Validator::new()` → `new(&ChemistryModel::default())`; references to renamed types/modules.

## Remaining
- umol-graph error/contradiction type survey (doc 065) — deferred; not a blocker for this work.

## Out of scope
- tier-0 validator (raise owns reference resolution; no AST to validate).
- Linter op category (needs broader cross-entity infra, not validator-specific).
- A production resolve→validate run site (composer-only this round).
