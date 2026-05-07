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

- **Aromatize**: collapse a Kekulé form (alternating single/double bonds) into an aromatic system + aromatic bonds.
- **Kekulize**: pick a Kekulé structure for each aromatic system; remove the system entry; assign bond orders 1/2 to the system's bonds.

Aromatize and Kekulize are *not* inverse pairs in general (kekulize loses the aromatic-system identity that aromatize would re-derive from topology + chemistry; aromatize commits a single canonical aromatic form when multiple were possible at the Kekulé level). They share a trait but have independent semantics.

Both operate on resolved input and produce resolved output. They are not part of the resolution pipeline.

### Trait

```rust
pub trait Transformer {
    fn transform_into(&self, ast: &mut MoleculeAst) -> Result<…>;

    fn transform(&self, ast: &MoleculeAst) -> Result<MoleculeAst, …> {
        let mut out = ast.clone();
        self.transform_into(&mut out)?;
        Ok(out)
    }

    fn generate_all<'a>(
        &'a self,
        ast: &'a MoleculeAst,
    ) -> Box<dyn Iterator<Item = MoleculeAst> + 'a>;
}
```

- `transform_into` — required, in-place. Picks the canonical result for transformers that have multiple (Kekulize: deterministic single matching in canonical atom order, per discussion 85 §"Phase 1").
- `transform` — value return, default impl via clone.
- `generate_all` — required, yields all results. Boxed iterator preserves dyn-compat for `Vec<Box<dyn Transformer>>`. For Aromatize the iterator yields one element; for Kekulize it enumerates Kekulé structures (Uno's algorithm, doc 85 §"All Kekulé structures").

### Concrete types

`AromaticityModel` is reused for aromatize — aromatize is perception applied to a Kekulé-form input. The `c1ccc[cH-]1` vs `C1=CC=C[CH-]1` distinction collapses: both parse into valid `MoleculeAst`s; the second form runs through aromatize to acquire aromatic systems if the user wants them.

```rust
pub struct Aromatize { model: AromaticityModel }
pub struct Kekulize { model: KekulizationModel }
```

### Module location

`umol-graph/src/ops/transform/aromatize.rs`, `umol-graph/src/ops/transform/kekulize.rs`. The `Transformer` trait and shared error types live in `ops/transform.rs` (parent module).

### SMILES input policy (point 2a)

Both `c1ccc[cH-]1` (lowercase, aromatic) and `C1=CC=C[CH-]1` (uppercase, Kekulé) parse to valid `MoleculeAst`s. The lowercase form already carries an aromatic system after parsing; the uppercase form does not, and the user runs Aromatize if they want one.

The SMILES parser already distinguishes lowercase aromatic atoms from uppercase. Verifying the current behavior end-to-end is a separate small task; the perception-fix work in §1 doesn't depend on it.

### Open: ordering vs §1

§2 reuses `AromaticityModel` and the perception path (Q3). The §1 pattern-match fix in `aromaticity.rs::resolve` lands first; §2 builds on the corrected behavior. Whether to overlap implementation in parallel branches or sequence them is a workflow question, not a design one.

---

## 3. Equalization rule: generality and limits (2026-05-07)

§1 introduced a pattern-match equalization with two patterns, `(+1, 0)` and `(-1, 2)` → `(0, 1)`. Reframing in terms of a per-atom invariant clarifies what the rule covers and what it cannot reach.

### The k_i invariant

For each atom in an aromatic system, define the local canonical occupancy `k_i`: the value `atom.charge + electrons[i]` would take if the atom were neutral and contributing canonically to the ring. Equalization rewrites `(q, e)` → `(0, k_i)` when `q + e == k_i` and `e` is offset by ±1 from `k_i` (the offset lives in π). The shifted charge accumulates into `system.charge`.

The current implementation hard-codes `k = 1`: it triggers only on `(+1, 0)` and `(-1, 2)`. The two `(q, e)` patterns above are the complete `k = 1` instance; nothing more is needed for `k = 1` alone.

Correctly handled:

- Atoms with canonical π = 1: sp² C, pyridine-like N (sp² N with the lone pair in σ), and analogous P, As at the same kind of position. `(charge, π)` triggers are `(+1, 0)` and `(-1, 2)`.
- ±n totals via per-atom accumulation: COT²⁻ has two atoms each at `(-1, 2)`; cyclooctatetraene dication would have two atoms each at `(+1, 0)`. `accumulated += c` rolls them into `system.charge`. No "n=2" branch is needed. (Benzene dication itself is anti-aromatic — 4 π electrons — so Hückel rejects it before equalization runs; not a useful test.)

Correctly left alone:

- σ-side charges on pyridine-like heteroatoms: pyridinium N, pyrylium O at `(+1, 1)` — `q+e=2`, no match.
- Canonical π=2 atoms: pyrrole N, furan O, thiophene S at `(0, 2)` — `q+e=2`, no match.

The `q+e=2` coincidence collapses both "leave alone" categories into the same skip path.

### Heterocyclic and ring-size coverage

Concrete examples within the `k = 1` rule, all carbocyclic:

- Cyclopropenium `[C₃H₃]⁺` (n=0, 2 π electrons): one C at `(+1, 0)` → `(0, 1)`.
- Tropylium `[C₇H₇]⁺` and Cp⁻ `[C₅H₅]⁻` (n=1, 6 π electrons): a single `(+1, 0)` or `(-1, 2)` atom respectively.
- Cyclononatetraenide `[C₉H₉]⁻` (n=2, 10 π electrons): one C at `(-1, 2)` → `(0, 1)`.
- COT²⁻ `[C₈H₈]²⁻` (n=2, 10 π electrons): two C at `(-1, 2)` → `(0, 1)` each, accumulating into `system.charge = -2`.

Heteroaromatic atoms with canonical π=1 (pyridine-like sp² N, P, As) work the same way in principle, but `(±1, 0)` / `(±1, 2)` configurations on these atoms are exotic — typical heteroaromatics keep the heteroatom at canonical `(0, 1)` or `(0, 2)`. Boron in the borabenzene/boratabenzene family does NOT fit the `k = 1` rule: B in pyridine-analog has canonical π = 0, not 1, so the boratabenzene anion's B is at `(-1, 1)` (q+e=0) — outside the rule's domain.

Untouched (rule correctly skips):

- pyridinium, pyrylium, thiopyrylium — heteroatom at `(+1, 1)`, σ-side charge.
- pyrrole, furan, thiophene — heteroatom at `(0, 2)`, k=2 canonical.
- boratabenzene anion — B at `(-1, 1)`, q+e=0, not a `k = 1` configuration at all.

### k = 2 case: S₄²⁺

Square-planar S₄²⁺ is the canonical `k = 2` system. 4 S, total charge +2, 6 π electrons (4n+2, n=1). Each S is a lone-pair donor (canonical π = 2 like furan O). Electron accounting: 22 valence e⁻ (4×6 − 2); σ skeleton 16 e⁻ (4 σ bonds + 4 σ in-plane lone pairs); π = 6. The fully equalized form has `(q=0, π=2)` on every atom with `system.charge = +2`; total π = 4×2 − 2 = 6.

A natural input description has two S at `(+1, 1)` (lost a lone-pair electron from π) and two S at `(0, 2)`. All atoms have `q + e = 2`; the `k = 1` rule does nothing. The molecule remains in the mixed input form — not wrong (total π is correct) but not uniform.

Reaching the equalized `(0, 2)`-uniform form requires per-atom `k_i`, which is element + σ-environment dependent and effectively another atom-typing concern. Deferred until concrete demand exists. The `k = 1` rule is documented as such; S₄²⁺ is recorded as a known gap.

### Test coverage to add

Tests will live alongside `equalize_charges`.

Equalizes (`k = 1`):

- cyclopropenium `[C₃H₃]⁺` — single `(+1, 0)` atom.
- cyclooctatetraenide dianion `[C₈H₈]²⁻` — two `(-1, 2)` atoms; demonstrates ±n accumulation.
- (tropylium and Cp⁻ are already covered indirectly by the conformance suite.)

Leaves alone:

- pyridinium, pyrylium.
- pyrrole, furan, thiophene.

Out of scope (`k = 2` known gap):

- S₄²⁺ — assert input passes through unchanged.

### Notation work, sequenced first

Writing these inputs against `MoleculeAst::new(vec![...], vec![...], ...)` with manual atom/bond construction is the friction the EDN DSL was supposed to remove. Before adding tests, decide:

1. **Inline EDN packaging** — per-module `fn mol(edn: &str) -> MoleculeAst` helper, or shared `umol_ast::test_support::mol`, or a `mol!` macro.
2. **Output assertion** — direct accessors (`ast.atom(i).data.charge` etc.), or EDN roundtrip equality, or insta snapshots.
3. **Aromatic-system input form** — the existing EDN `:aromatic [{:atoms [...] :type "#e2" :charge n}]` block already supports declaring a pre-existing aromatic system, which is what these equalization tests need (input has the system inserted, run equalization, check resulting fields).

Picks 1 and 2 are independent.
