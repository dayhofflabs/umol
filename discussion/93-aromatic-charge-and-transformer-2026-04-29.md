# Aromatic charge equalization, and aromatize/kekulize transformer

## 1. Charge-delocalized aromatic perception

### What the perception code already does

Hückel rule (`hueckel_rule.rs::ring_electron_count`) and HMO (`hmo.rs::build_calculator`) both compute total π for the candidate system as `sum(aromatic_pi_contribution(atom))`. The per-atom `#a` value already encodes local charge effects through the registry (Cp⁻ atom `C #c- #v3 #a2` contributes 2; tropylium atom `C #c+ #v3 #a0` contributes 0; pyridinium N `N #c+ #v2 #a #h` contributes 1). Summing these gives the correct total π in every case. No subtraction of any charge field is involved at the perception step, because no system exists yet — charge lives on the atoms.

### What the resolver does after perception (and where it stops)

After the algorithm reports a system, `aromaticity.rs::resolve` (lines 110–151):

1. Sums per-atom charges → writes that into `system.charge`.
2. Sums per-atom unpaired electrons → derives `system.spin`.
3. Resets every system atom's `charge = 0` and `spin = closed_shell`.
4. Leaves per-atom `electrons` (the `AromaticSystemAst.electrons` vector) at the localized values.

Step 4 is the gap. After charge has been moved up to the system, per-atom `electrons` still encode the localized resonance form. Conservation of total π then has to read both fields together with the charge subtracted back out:

```
total_pi  =  sum(electrons)  -  system.charge
```

For Cp⁻ as currently rendered (`electrons = [2 1 1 1 1]`, `system.charge = -1`), this gives `6 - (-1) = 7` — wrong. The "right" interpretation is "ignore `system.charge` and just sum `electrons`", but then `system.charge` is redundant bookkeeping that contradicts itself.

### What changes: pattern-match equalization

Per atom in the system, examine the pair `(atom.charge, electrons[i])` and rewrite only the two patterns that correspond to a π-electron-displaced atom:

| Match | Rewrite |
|-------|---------|
| `(+1, 0)` (tropylium-style C⁺) | `(0, 1)`, `system.charge += +1` |
| `(-1, 2)` (Cp⁻-style C⁻)       | `(0, 1)`, `system.charge += -1` |
| anything else                   | leave atom and `electrons[i]` untouched |

Pyridinium / pyrylium-type heteroatoms with `(+1, 1)` or `(-1, 1)` fall through unchanged: their per-atom contribution already matches the canonical 1, so the charge is σ-localized and stays on the atom.

Each rewrite preserves `electrons[i] - atom.charge`, so `sum(electrons) - system.charge` stays invariant per atom and therefore over the whole system. Total π is conserved.

| Case        | electrons (post) | system.charge | per-atom charges (post) | total π |
|-------------|------------------|---------------|-------------------------|---------|
| Cp⁻         | `[1 1 1 1 1]`    | -1 | all 0       | 5 - (-1) = 6 ✓ |
| Tropylium   | `[1 1 1 1 1 1 1]`| +1 | all 0       | 7 - 1 = 6 ✓ |
| Benzene     | `[1 1 1 1 1 1]`  | 0  | all 0       | 6 - 0 = 6 ✓ |
| Pyridinium  | `[1 1 1 1 1 1]`  | 0  | N=+1, C=0   | 6 - 0 = 6 ✓ |
| Pyrylium    | `[1 1 1 1 1 1]`  | 0  | O=+1, C=0   | 6 - 0 = 6 ✓ |

Two behavioral changes vs current code:
1. Charge: the resolver currently zeroes every aromatic atom's charge unconditionally, transferring everything into `system.charge`. Under the new rule only the matched atoms are zeroed, so pyridinium-like σ-charges stay on the heteroatom.
2. Spin: the resolver currently sums all atom unpaired electrons into `system.spin` and resets every atom to closed-shell. Drop this entirely. Per-atom spin stays where the input put it; `system.spin` keeps its construction default (closed-shell). No analogous pattern rule is introduced for spin.

### Equalization step

The transformation that satisfies the new rule and conserves total π:

```
for atom i in system.atoms:
    new_electrons[i]  =  old_electrons[i] + atom[i].charge
    atom[i].charge    =  0
system.charge  =  sum(old atom[i].charge)
```

Conservation check: `sum(new_electrons) - system.charge = sum(old_electrons) + sum(charges) - sum(charges) = sum(old_electrons)`. The right-hand side equals the original total π, so the count is preserved.

Applied to the cases above:

| Case        | new electrons         | system.charge | total |
|-------------|-----------------------|---------------|-------|
| Cp⁻         | `[1 1 1 1 1]`         | -1 | 5 - (-1) = 6 ✓ |
| Tropylium   | `[1 1 1 1 1 1 1]`     | +1 | 7 - 1 = 6 ✓ |
| Benzene     | `[1 1 1 1 1 1]`       | 0  | 6 ✓ |
| Pyridinium  | `[2 1 1 1 1 1]` (N→2) | +1 | 7 - 1 = 6 ✓ |

The first three rows match the user-prescribed output. The fourth row is the pyridinium edge case.

### Cascade of changes

1. `aromaticity.rs::resolve` — replace the current unconditional charge/spin aggregation with the pattern-match rule: walk system atoms, apply the two charge patterns, leave everything else (including all spin handling) untouched. The per-atom-spin reset and `system.spin` aggregation block goes away.
2. Test fixtures `cyclopentadienyl_anion` and `tropylium` in `hueckel_rule.rs::tests` — current `cyclopentadienyl_anion` puts `(C, charge=-1, #a2)` on every atom (sum 10, trivially `4n+2` for `n=2`); not representative of any real input. Replace with the realistic inputs the user described.
3. Conformance snapshots `aromatic_cyclopentadienyl-anion.snap`, `aromatic_tropylium.snap`, and any other charged-aromatic snapshot — expected diff: per-atom strings collapse to the canonical `C#h#v2#a`, `:electrons` flatten to all-1s.
4. `AromaticSystemConstraint::ElectronCount` (the `#e<n>` predicate at system level) — the validation rule it must enforce is `ElectronCount == sum(electrons) - system.charge`. Status of the codebase against this rule:

   - The docstring at `umol-ast/src/ast/constraint/aromatic.rs:12` claims cross-checking by a `ConsistencyValidator` against `sum(electrons)`. Inaccurate on two counts: no `ConsistencyValidator` exists in the workspace; and the formula omits `- system.charge`.
   - The validator that took its place per discussion 92, `ConstraintValidator` (`umol-graph/src/ops/validator.rs:391–398`), has a stub body returning `Solution::Determined(())` — no `ElectronCount` (or any other constraint) is consumed.
   - Result: no object currently obeys, or contradicts, the rule. Fixing it means (a) updating the docstring to state the rule including `- system.charge`, and (b) wiring the actual check into `ConstraintValidator`. The same applies to `MulticenterBondConstraint`'s analogous `#e` (the multicenter docstring at `constraint/multicenter.rs:12` carries the same broken `ConsistencyValidator` reference; needs the same correction once we know the multicenter-side rule).
5. Atom-typing registry — Cp⁻ and tropylium localized-form patterns (`C #c- #v3 #a2`, `C #c+ #v3 #a0`, etc.) remain the patterns matched on *input* atoms before perception runs. After perception their atoms switch to the canonical `C #v3 #a` / `C #v2 #a #h` patterns. No registry edit is required for the perception fix alone; whether to keep the localized patterns as legitimate input forms is a separate curation question.

### What I want to confirm before coding

- Nothing on §1 — confirmed: charge by pattern match, spin untouched.

---

## 2. Aromatize / kekulize transformer

### Scope

Two operations on a fully resolved `MoleculeAst`:

- **Aromatize**: collapse a Kekulé form (alternating single/double bonds) into an aromatic system + aromatic bonds. Inverse of kekulize.
- **Kekulize**: pick one Kekulé structure for each aromatic system; remove the system entry; assign bond orders 1/2 to the system's bonds.

Both operate on resolved input and produce resolved output. They are not part of the resolution pipeline.

User-facing shape (given): `Transformer::new(&model).transform(&ast) -> ast`.

Open design questions are below. Each carries options, no recommendation.

### Q1. One transformer trait, or two distinct operations?

a. **Single `Transformer` trait** with `Aromatize` and `Kekulize` enum variants (parallel to `AromaticityResolver`'s current shape). One entry point, dispatch by config.

b. **Two separate types**: `Aromatizer::new(&AromatizationModel).transform(...)`, `Kekulizer::new(&KekulizationModel).transform(...)`. No shared trait.

c. **`Transformer` trait** parameterized by an enum of all possible transforms (current and future: tautomerize, normalize, etc.).

The current `AromaticityResolver` follows shape (a). `Resolver` is a struct enum dispatching by model. Pattern continuity suggests (a); composability with future transforms suggests (c); minimum surface suggests (b).

### Q2. In-place vs return-new

a. `fn transform(&self, ast: &mut MoleculeAst) -> Result<Solution<...>, ...>` — mirrors the resolver shape.

b. `fn transform(&self, ast: &MoleculeAst) -> Result<MoleculeAst, ...>` — value return, matches the user's stated shape.

(b) is what the user wrote. (a) is what the codebase already uses for resolvers. The two differ in how a partial-failure result is reported (Solution vs Result) and in cost (no clone vs always clone).

### Q3. Aromatize: what model parameters?

The aromatize step needs the same chemistry policy as the perception phase of resolution. Two options:

a. Reuse `AromaticityModel` — aromatize is "just run perception on a Kekulé input". The `c1ccc[cH-]1` vs `C1=CC=C[CH-]1` distinction collapses: both run through perception, and the user simply doesn't run aromatize on the second form. This makes aromatize a thin wrapper around the existing resolver.

b. Distinct `AromatizationModel` — aromatize accepts a Kekulé form and only that. Adds input validation: reject input that already has aromatic systems / aromatic bonds. Conceptually cleaner separation, but duplicates most of perception's logic.

### Q4. Kekulize: which Kekulé structure?

Discussion 85 already covers the matching algorithms. The transformer-level question is which one to surface:

a. **Single canonical structure** (RDKit-style greedy DFS in canonical atom order). Deterministic, one output. Discussion 85 §"Phase 1".

b. **All Kekulé structures** (Uno enumeration). Returns `Vec<MoleculeAst>` instead of one. Different return type.

c. **One configurable** — `KekulizationModel::SingleCanonical` vs `KekulizationModel::All`, dispatched at the transformer level.

Phase 1 of discussion 85 is single-canonical. (c) leaves the door open without committing to (b) now; (a) commits to single-output and forces a separate API later for enumeration.

### Q5. Where does the transformer module live?

a. `umol-graph/src/ops/transform/aromatize.rs`, `umol-graph/src/ops/transform/kekulize.rs` — sibling to `ops/aromaticity.rs`, `ops/valence.rs`, etc.

b. `umol-graph/src/transform.rs` (top-level) — distinct from `ops` (which is resolution-pipeline territory).

c. `umol-graph/src/ops/aromaticity/` — both aromatize and kekulize are aromatic-system transformations; live next to perception.

`ops/` currently contains the resolver pipeline. If transformers are explicitly *not* part of resolution, (b) keeps that boundary visible; (a) keeps everything in `ops/` for discoverability; (c) groups by chemistry topic.

### Q6. SMILES input policy (user's point 2a)

The user stated that both `c1ccc[cH-]1` (lowercase, aromatic) and `C1=CC=C[CH-]1` (uppercase, Kekulé) must parse to valid `MoleculeAst`s, with the first triggering perception and the second not.

The SMILES parser already distinguishes lowercase aromatic atoms from uppercase. Confirming the current behavior and writing it down explicitly is a separate small task; the perception-fix work in §1 doesn't depend on it.

### What I want to confirm before coding

- Q1, Q2, Q3, Q4, Q5 each with one of the listed options (or another).
- Whether the transformer work blocks on §1 (perception fix) or can proceed in parallel.
