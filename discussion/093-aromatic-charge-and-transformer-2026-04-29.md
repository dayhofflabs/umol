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
| `(+1, 0)` (tropylium-class C⁺) | `(0, 1)`, `system.charge += +1` |
| `(-1, 2)` (Cp⁻-class C⁻)       | `(0, 1)`, `system.charge += -1` |
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

## 3. Equalization rule: delocalization criterion (2026-05-07)

§1 introduced a pattern-match equalization. After working through the test cases, the principled framing is **delocalization**: equalization captures the case where the system's charge has no chemically meaningful localization on any specific atom, because the atoms are π-MO-equivalent. Operationally that's the **monoelement** condition — every atom in the aromatic system has the same element. In a heterogeneous system the heteroatom is the natural locus of the molecule's charge (pyridinium's `+1` on N, pyrylium's on O, boratabenzene's `−1` on B), and equalization would erase that.

### Rule

For each aromatic system:

1. If every member atom has the same `Element::Lit` value:
   - For each atom, compute `K = V − σ − 2·LP` from closed-shell electron accounting, with `V = Element::valence_electrons()`, `σ` the σ-bond count (graph degree + implicit H), `LP` the σ-lone-pair count.
   - If `electrons[i] ≠ K`, rewrite the atom to `(q = 0, π = K)` and accumulate the per-atom `q` into `system.charge`.
2. Otherwise (heterogeneous): skip the system. Charges stay on whichever atom carried them.

The K formula is element-uniform — `V` is read from `umol_shared::Element`, the rest is atom state. No element-specific table or scope list is buried in code; the only datum that varies per element is V, which is intrinsic chemistry and already lives in shared element data.

### Concrete coverage

Equalizes (monoelement):

- Cyclopropenium `[C₃H₃]⁺` (3 π in 3-ring): one C at `(+1, 0)` → `(0, 1)`. K_C = 4 − 3 − 0 = 1.
- Tropylium `[C₇H₇]⁺` and Cp⁻ `[C₅H₅]⁻` (6 π): single `(+1, 0)` or `(-1, 2)` atom equalizes.
- Cyclononatetraenide `[C₉H₉]⁻` (10 π): one C at `(-1, 2)` equalizes.
- COT²⁻ `[C₈H₈]²⁻` (10 π): two C at `(-1, 2)`, both equalize, system.charge = −2.
- `[S₄]²⁺` (square-planar, 6 π in 4-ring): K_S = 6 − 2 − 2 = 2; the two `(+1, 1)` atoms equalize to `(0, 2)`; system.charge = +2.

Skips (heterogeneous):

- Pyridinium `[C₅H₅NH]⁺`, pyrylium `[C₅H₅O]⁺`, thiopyrylium — N/O/S keeps its charge.
- Boratabenzene anion `[C₅H₅BR]⁻`, borabenzene cation — B keeps its charge.
- Pyrrole, furan, thiophene, borepin, borazine, 1,2-azaborine, boroxine — no charge to equalize anyway, but the heterogeneous check makes the no-op explicit.

±n accumulation (e.g., COT²⁻'s `−2`, S₄²⁺'s `+2`) falls out of per-atom `accumulated += q`. No special handling for n.

### Why this answers the σ/π-origin question

The atom state `(q, e, σ, LP)` doesn't disambiguate a σ-side charge (pyridinium = pyridine + H⁺) from a π-side charge in the same per-atom shape (hypothetical "oxidized pyrrole-N" at the same `(q, π, σ, LP)`). Tying the rule to the monoelement condition sidesteps the ambiguity: whenever the rule fires, every atom in the ring is the same element, so per-atom localization of the system's charge can only be a symmetry-breaking artifact (the atoms are equivalent under the ring's automorphisms), not a chemical fact about which atom the charge "actually" sits on. In heterogeneous systems the heteroatom *is* a chemical fact about the charge's location, and skipping is correct.

### Tests

Implemented in `umol-graph/src/ops/aromaticity.rs::tests`:

Equalizes:
- `test_equalize_charges_cyclopropenium_cation`
- `test_equalize_charges_cot_dianion_accumulates_to_minus_two`
- `test_equalize_charges_s4_dication`

Skips (heterogeneous):
- `test_equalize_charges_boratabenzene_anion_leaves_b_charged`
- `test_equalize_charges_borepin_neutral_aromatic_no_op`
- `test_equalize_charges_pyridinium_leaves_n_charged`
- `test_equalize_charges_pyrylium_leaves_o_charged`
- `test_equalize_charges_pyrrole_leaves_n_at_pi_two`
- `test_equalize_charges_furan_thiophene_leave_heteroatom_alone` (rstest, two cases)

All inputs use the `mol!` macro from doc 94's notation work; assertions read AST via direct accessors.

### Open

- AtomView lacks a `sigma_bond_count()` accessor. The rule currently inlines `view.neighbors().count() + implicit_h` with a `TODO` to replace once the AtomView API is fleshed out. Tracked separately as part of the umol-ast view-API completion.
- The discussion of σ/π-origin distinction is folded into the monoelement framing here; if the σ/π origin ever needs to be preserved in the AST (e.g., to support equalization within heteroaromatics), it's a separate AST-level decision (extra bits per atom recording charge origin), not a refinement of this rule.
