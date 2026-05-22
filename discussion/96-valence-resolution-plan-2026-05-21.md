# Valence resolution plan

Design and case analysis for the valence resolvers (counts, atom-typing). Covers (i) the per-atom electron-conservation invariants from [doc 52 §10.1.3](52-graph-ir-2026-02-11.md#1013-invariants), (ii) the case-by-case behaviour of `implicit_h`, `lone_pairs`, `unpaired`, `charge` resolution under partially-pinned inputs, and (iii) the proposed architectural refactor that pulls the conservation equations onto `ValenceModel` so resolver and validator share one source of truth. Inputs are read under the AST's `Default::new()` open-defaults convention — absent tags mean `Undetermined`, not zero. SMILES bare-atom semantic (`#h* #u0`) and the MOL implicit-H flag live at the parser boundary; this doc is about what the AST-level resolver sees after parsing.

## Setup

The per-atom electron-conservation equation, from the orbital-count and electron-count invariants in [doc 52 §10.1.3](52-graph-ir-2026-02-11.md#1013-invariants):

```
Z(element) − charge = valence + implicit_h
                    + aromatic_valence + multicenter_valence
                    + 2·lone_pairs + 2·donated_pairs + unpaired
```

with side constraints:

- each term ≥ 0
- parity `unpaired ↔ multiplicity`
- `valence` ∈ `allowed_valences(element)`
- `aromatic_valence` ∈ {0, 1, 2} (per doc 52 fn 8)
- `multicenter_valence` ≥ 0

Accepted pairs (`a`, dative bond acceptor) do not appear in this equation — they're balanced by the borrowed-electron increment `+2a` on the electron-count side and cancel out. Donated pairs (`d`, dative donor) contribute `2d` because the donor keeps the electrons in its orbital while sharing them.

The cases below isolate one slice at a time. Isolated atoms — `valence = aromatic_valence = multicenter_valence = donated_pairs = 0` — live in §C / §N / §B, where the equation reduces to `Z − q = h + 2n + u`. Bonded (`v > 0`), aromatic (`aromatic_valence > 0`), dative (`d > 0`), and multicenter (`mc_v > 0`) tables to come.

DSL surface: `#u` parses as `unpaired = Lit(1)`; `#u0` as `Lit(0)`; `#u2` as `Lit(2)`; `#u?u > 0` as expression-constrained. Same shape for `#h`, `#c`, `#n`, `#v`, `#a`.

## Conventional `lp_conv` per element (reference)

Lone-pair count of the neutral closed-shell ground state:

| element | Z | nv(neutral) | lp_conv(neutral) |
|---|---|---|---|
| B  | 3 | 3 | 0 |
| C  | 4 | 4 | 0 |
| N  | 5 | 3 | 1 |
| O  | 6 | 2 | 2 |
| F  | 7 | 1 | 3 |
| Ne |10 | 0 | 4 |

The previous shift-by-charge stopgap (look up `lp_conv` of the isoelectronic neutral via `element.shift(-charge)`) and the `NormalImplicitHydrogensTable` together implemented the now-removed `Normal` sentinel — they picked `lp = lp_conv` and derived `h` from conservation. With `Normal` removed, `lone_pairs` (`n`) is a free conservation variable on the same footing as `h`, `c`, `u`; this table is reference data, not the spec for resolution.

## C cases (Z = 4, isolated atom: v = aromatic_valence = mc_v = d = 0)

Conservation reduces to `4 − q = h + 2n + u`. Inputs under `Default::new()` open-defaults — absent tags are `Undetermined`.

| input | parsed `(h, c, u, n)` | conservation | outcome |
|---|---|---|---|
| `C#h4#c0#u0` | (4, 0, 0, ?) | `4 = 4 + 2n` → n=0 | commit (4, 0, 0, 0) — CH₄ |
| `C#h5#c0#u0` | (5, 0, 0, ?) | `4 = 5 + 2n` → 2n=−1 | infeasible |
| `C#h3#c0#u0` | (3, 0, 0, ?) | `4 = 3 + 2n` → 2n=1 | parity violation |
| `C#h*#c0#u0#n0` | (?, 0, 0, 0) | `4 = h` | commit (4, 0, 0, 0) — CH₄ |
| `C#h*#c0#u0#n1` | (?, 0, 0, 1) | `4 = h + 2` → h=2 | commit (2, 0, 0, 1) — singlet :CH₂ |
| `C#h*#c0#u0#n2` | (?, 0, 0, 2) | `4 = h + 4` → h=0 | commit (0, 0, 0, 2) — :C: |
| `C#h*#c+#u0#n0` | (?, +1, 0, 0) | `3 = h` | commit (3, +1, 0, 0) — CH₃⁺ |
| `C#h*#c+5#u0#n0` | (?, +5, 0, 0) | `−1 = h` | infeasible |
| `C#h*#c-#u0#n1` | (?, −1, 0, 1) | `5 = h + 2` → h=3 | commit (3, −1, 0, 1) — CH₃⁻ |
| `C#h*#c-2#u0#n2` | (?, −2, 0, 2) | `6 = h + 4` → h=2 | commit (2, −2, 0, 2) — CH₂²⁻ |
| `C#h*#c0#u#n0` | (?, 0, 1, 0) | `4 = h + 1` → h=3 | commit (3, 0, 1, 0) — ·CH₃ |
| `C#h*#c0#u3#n0` | (?, 0, 3, 0) | `4 = h + 3` → h=1 | commit (1, 0, 3, 0) — ·CH quartet |
| `C#h*#c0#u4#n0` | (?, 0, 4, 0) | `4 = h + 4` → h=0 | commit (0, 0, 4, 0) — C atom quintet |
| `C#h*#c0#u5#n0` | (?, 0, 5, 0) | `4 = h + 5` → h=−1 | infeasible |
| `C#h*#c0#u0` | (?, 0, 0, ?) | `4 = h + 2n` | underdetermined: `(h, n) ∈ {(4,0), (2,1), (0,2)}` |
| `C#h*#c+#u` | (?, +1, 1, ?) | `2 = h + 2n` | underdetermined: `(h, n) ∈ {(2,0), (0,1)}` |
| `C#h*#c0` | (?, 0, ?, ?) | `4 = h + 2n + u` | underdetermined (3 free vars) |
| `C#h*` | (?, ?, ?, ?) | `4 − q = h + 2n + u` | underdetermined (4 free vars) |

## N cases (Z = 5, isolated atom: v = aromatic_valence = mc_v = d = 0)

Conservation reduces to `5 − q = h + 2n + u`.

| input | parsed `(h, c, u, n)` | conservation | outcome |
|---|---|---|---|
| `N#h4#c0#u0` | (4, 0, 0, ?) | `5 = 4 + 2n` → 2n=1 | parity violation |
| `N#h5#c0#u0` | (5, 0, 0, ?) | `5 = 5 + 2n` → n=0 | commit (5, 0, 0, 0); validator rejects as hypervalent (`allowed_valences(N) = [3]`) |
| `N#h3#c0#u0` | (3, 0, 0, ?) | `5 = 3 + 2n` → n=1 | commit (3, 0, 0, 1) — NH₃ |
| `N#h*#c0#u0#n1` | (?, 0, 0, 1) | `5 = h + 2` → h=3 | commit (3, 0, 0, 1) — NH₃ |
| `N#h*#c+#u0#n0` | (?, +1, 0, 0) | `4 = h` | commit (4, +1, 0, 0) — NH₄⁺ |
| `N#h*#c-#u0#n2` | (?, −1, 0, 2) | `6 = h + 4` → h=2 | commit (2, −1, 0, 2) — NH₂⁻ |
| `N#h*#c-2#u0#n3` | (?, −2, 0, 3) | `7 = h + 6` → h=1 | commit (1, −2, 0, 3) — NH²⁻ |
| `N#h*#c+#u#n0` | (?, +1, 1, 0) | `4 = h + 1` → h=3 | commit (3, +1, 1, 0) — ·NH₃⁺ |
| `N#h*#c-#u#n1` | (?, −1, 1, 1) | `6 = h + 1 + 2` → h=3 | commit (3, −1, 1, 1) — ·NH₃⁻ |
| `N#h*#c0#u4` | (?, 0, 4, ?) | `5 = h + 2n + 4` → h+2n=1 | conservation pins (h, n) = (1, 0); commit (1, 0, 4, 0) — quintet ·NH |
| `N#h*#c0#u5` | (?, 0, 5, ?) | `5 = h + 2n + 5` → h+2n=0 | conservation pins (0, 0); commit (0, 0, 5, 0) — quintet N atom |
| `N#h*#c0#u0` | (?, 0, 0, ?) | `5 = h + 2n` | underdetermined: `{(5,0), (3,1), (1,2)}` |
| `N#h*#c0` | (?, 0, ?, ?) | `5 = h + 2n + u` | underdetermined |
| `N#h*` | (?, ?, ?, ?) | `5 − q = h + 2n + u` | underdetermined |

## B cases (Z = 3, isolated atom: v = aromatic_valence = mc_v = d = 0)

Conservation reduces to `3 − q = h + 2n + u`.

| input | parsed `(h, c, u, n)` | conservation | outcome |
|---|---|---|---|
| `B#h4#c0#u0` | (4, 0, 0, ?) | `3 = 4 + 2n` → 2n=−1 | infeasible |
| `B#h3#c0#u0` | (3, 0, 0, ?) | `3 = 3 + 2n` → n=0 | commit (3, 0, 0, 0) — BH₃ |
| `B#h*#c0#u0#n0` | (?, 0, 0, 0) | `3 = h` | commit (3, 0, 0, 0) — BH₃ |
| `B#h*#c+#u0#n0` | (?, +1, 0, 0) | `2 = h` | commit (2, +1, 0, 0) — BH₂⁺ |
| `B#h*#c-#u0#n0` | (?, −1, 0, 0) | `4 = h` | commit (4, −1, 0, 0) — BH₄⁻ |
| `B#h*#c-#u0#n1` | (?, −1, 0, 1) | `4 = h + 2` → h=2 | commit (2, −1, 0, 1) — BH₂⁻ |
| `B#h*#c0#u3#n0` | (?, 0, 3, 0) | `3 = h + 3` → h=0 | commit (0, 0, 3, 0) — B atom quartet |
| `B#h*#c0#u4#n0` | (?, 0, 4, 0) | `3 = h + 4` → h=−1 | infeasible |
| `B#h*#c0#u0` | (?, 0, 0, ?) | `3 = h + 2n` | underdetermined: `{(3,0), (1,1)}` |
| `B#h*` | (?, ?, ?, ?) | `3 − q = h + 2n + u` | underdetermined |

## Bonded-atom cases (`?v > 0`) — TBD

For atoms with localized valence > 0 but no aromatic, multicenter, or dative participation, conservation becomes `Z − q − v = h + 2n + u`. Tables to enumerate per element × `v ∈ allowed_valences(element)`. Examples to fill in:

- Methyl C (`v=1` in CH₃–X): `C#v1#c0#u0#n0` → `3 = h` → CH₃ commit.
- Methylene C (`v=2` in =CH₂): `C#v2#c0#u0#n0` → `2 = h` → CH₂ commit.
- Carbonyl C (`v=3` in chains like RC(=O)R): `C#v3#c0#u0#n0` → `1 = h` → CH commit; `C#v4#c0#u0#n0` → `0 = h` → no implicit H.
- NH in chain (`v=2`): `N#v2#c0#u0#n1` → `1 = h` → NH commit.
- Tetracoordinated N⁺ (`v=4`): `N#v4#c+#u0#n0` → `0 = h` → ammonium commit.

Pattern same as isolated: pinned `(c, u, n, v)` → conservation determines `h`; multiple unknowns → underdetermined.

## Aromatic-atom cases (`?aromatic_valence > 0`) — TBD

For atoms in an aromatic system (`aromatic_valence ∈ {1, 2}`), no multicenter or dative participation, conservation becomes `Z − q − v − aromatic_valence = h + 2n + u`. Tables to enumerate per element × `(v, aromatic_valence)` allowed combinations. Examples to fill in:

- Benzene C (`v=2`, `a=1`): `C#v2#a1#c0#u0#n0` → `1 = h` → aromatic CH commit.
- Pyridine N (`v=2`, `a=1`): `N#v2#a1#c0#u0#n1` → `0 = h` → no implicit H commit; lp=1.
- Pyrrole N (`v=2`, `a=2`): `N#v2#a2#c0#u0#n0` → `1 = h` → NH commit; lp=0 (donated to system).
- Furan O (`v=2`, `a=2`): `O#v2#a2#c0#u0#n1` → `0 = h` → commit; lp=1 (one of two donated).
- Pyrrolyl-like anion N (`v=2`, `a=2`, `c=−1`): `N#v2#a2#c-#u0#n2` → `0 = h` → commit; lp=2.

## Observations

1. **Conservation alone leaves most partial inputs underdetermined.** Without `n` pinned, even `(h, c, u)` fully pinned leaves `n` to be derived from conservation — one unknown is fine. But the typical SMILES-derived AST input pins `(c, u)` (parser fills in `c=0`, `u=0` defaults at the surface) and leaves `(h, n)` open → two unknowns → underdetermined unless additional inputs narrow.
2. **Pinned `h` (e.g., `#h4`) ⇒ resolver leaves `h` alone**; conservation pins `n` if `c, u` are also pinned. Hypervalent rejection (e.g., NH₅) is the validator's job — resolver doesn't enforce `allowed_valences`.
3. **Conservation can pin without all-pinned input**: when the residual after constraint substitution is small (≤ 0), the only feasible `(h, n, u)` may be unique. `N#h*#c0#u4` has `h + 2n = 1` → forces `(h=1, n=0)`. Resolver commits.
4. **Parity violations** (odd residual where lp doubles) → reject as infeasible. Distinct from constraint conflict — pure arithmetic.
5. **`#u?u > 0`** (positive lower bound, specific count undetermined) is non-Lit; conservation can't substitute it. Same handling as any other free variable — underdetermined until a downstream pass pins.

## Open design questions

1. **Underdetermined output shape.** The central question after dropping `Normal`: when multiple `(h, n, u, q)` tuples satisfy conservation + inputs, what does the resolver commit? Options enumerated in the design draft below.
2. **`lp_conv` table use.** The reference table above is still meaningful for the closed-shell neutral default — a chemist asking "how many H does C have" expects 4, not "underdetermined". The question is whether to fold this default into the resolver (defeats the `Default::new()` open-defaults principle) or to surface it only as a separate "applied normal-defaults" transformation that callers opt into.
3. **Aromatic-specific lp accounting.** Pyrrole-N has `lp=0` (lone pair donated into the ring); pyridine-N has `lp=1`. The aromatic table column previously called `aromatic_normal_valence` encoded this. Same `lp_conv` question (point 2) applies to aromatic atoms.

## Resolver design draft (partial narrowing, fixed-point with aromaticity)

The valence resolvers (counts and atom-typing) internally enumerate `(h, n, u, q, v, a)` candidates against conservation + the element's `allowed_valences` and `allowed_aromatic_valences`. The candidate set has three sizes:

- **0** — input is infeasible (conservation impossible) or violates `allowed_valences`. Resolver returns an error. Current behavior; keep.
- **1** — input is feasible and conservation + constraints pin a unique candidate. Resolver narrows the atom. Current behavior; keep.
- **> 1** — input is underdetermined. **This is the case the design draft covers.**

### Three options for size > 1

a. **Error.** Current behavior in `counts.rs::candidates_for` and `atom_typing.rs::resolve`. Conservative; rejects realistic inputs (a neutral `C` with no `n` constraint has 3 candidates: methane, singlet :CH₂, atomic carbon).

b. **Partial narrowing (meet over candidates).** Write the meet of candidate field values into the atom. Fields where all candidates agree become `Lit`; fields where they differ stay `Undetermined` (or become `LitSet` if we want to record the alternatives). Atom remains non-ground; downstream passes may close the gap.

c. **Conservation-only narrowing.** Narrow only what conservation forces *regardless of which candidate is chosen*. E.g., if `c` and `u` are pinned and the only free vars are `h` and `n`, conservation pins `h + 2n`; if every candidate happens to share `h` (different `n`), commit `h`. Otherwise no-op.

(b) and (c) often coincide; (c) is the safer subset.

### Fixed-point interaction with aromaticity

When candidate sets are size > 1 after the first sweep, aromaticity perception can prune them on the second sweep:

- A neutral C with candidate `(h, n, u) ∈ {(4,0,0), (2,1,0), (0,2,0)}` after first counts pass reduces to `{(2,0,0,a=1)}` once Hückel assigns it to an aromatic ring (a=1, forced n=0).
- Charge equalization in aromatic systems (`equalize_charges` in `ops/aromaticity.rs`) writes per-atom `c` after the system lands. A counts pass that left `h` underdetermined because `c` was free can now narrow `h`.

The fixed-point driver wraps `Resolver::resolve` in `while any_resolver_made_progress { ... }`. Phase 14's `narrow_from` is monotone on the meet-semilattice, so the chain of refinements is finite and terminates.

Implicit-H itself doesn't *need* fixed-point — it's a derived value from conservation, not a separate resolver. But under options (b)/(c), the implicit-H decision may legitimately defer until another resolver pass pins the inputs that close conservation.

### Enumeration as method, not phase

The enumeration of candidates lives inside the existing `counts.rs::candidates_for` and `atom_typing.rs::prepare_atom` — they already do this for the size-1 case. The design extension is the output shape for size > 1 (option a/b/c) and the entry/exit conditions for the fixed-point loop driving them.

(Needs concrete examples per element, per (v, aromatic_valence) combination, and a worked-through fixed-point trace before turning into code. Aromatic charge equalization → counts re-pass is the canonical motivating example.)

## Proposed architecture: invariants on `ValenceModel`

The conservation equations and the enumeration solver they drive belong on `ValenceModel` (in `umol-graph/src/ops/config.rs`). Resolver and validator both call through; neither inlines the equations. Goal: a single source of truth so resolver and validator never diverge.

### Architectural principle: resolver / validator are model-agnostic shells

**The resolver and validator never reach into model-specific data.** They iterate atoms and call `ValenceModel` API methods; everything else — `AtomTypeRegistry` lookups, `ValenceTable` allowed-valence lists, `valence_capacity` element bounds, the conservation equations themselves — lives behind that API.

Consequence: switching `ValenceModel::AtomTyping` to `ValenceModel::Counts` (or to a future variant) requires *zero* changes in the resolver / validator loops. Adding a new constraint kind to atom-typing (e.g., per-element ring-membership priors) is a model-internal change; the resolver/validator are oblivious.

Things that move behind the model API:

- `counts.rs::try_build_candidate` and `resolve_unpaired_lone_pairs` — the conservation arithmetic — collapse into `model.resolve()` for the Counts variant.
- `atom_typing.rs::prepare_atom` — the registry-lookup setup — collapses into `model.resolve()` for the AtomTyping variant.
- `validator/invariant.rs::validate{,_atom}` — the inline orbital/electron-count equations — become a single `model.invariants().check()` call.
- Any future per-atom chemistry priors land as new private state on a `ValenceModel` variant, not as new logic in the resolver/validator.

This is the load-bearing contract for the refactor. If a `pub` method on `ValenceModel` would surface implementation detail (e.g., `model.registry() -> &AtomTypeRegistry`), it's a code smell — fold the operation into a higher-level method instead.

### Two-method contract per model: resolve + validate

Each `ValenceModel` variant exposes a matched pair of operations on the same set of atoms:

| variant | `resolve(partial)` | `validate(full)` |
|---|---|---|
| `AtomTyping` | filter registry patterns by `prepared.matches(pat)`; narrow / output JointDomain per match count | atom must match at least one registry pattern for its `(element, charge)` |
| `Counts` | enumerate via `invariants.solve(view)` — conservation + element bounds | atom satisfies `invariants.check(view)` (orbital count = electron count + aromatic-bound check) |

Resolve and validate are the same chemistry seen from two directions: resolve fills in missing information against the model's allowed states; validate checks that pinned information is in an allowed state. Both must agree — if `resolve` would have admitted the atom as a candidate, `validate` must pass; if `resolve` would have rejected, `validate` must fail.

The conservation invariant (`Invariants::check`) is universal physics — it applies regardless of which `ValenceModel` is active. `Counts::validate` IS the conservation check. `AtomTyping::validate` is the registry-membership check; conservation is *additionally* enforced by the standalone `ElectronInvariantValidator` (which also goes through `model.invariants().check()`).

### Current state (what to fix)

Equations and chemistry data are reached into from outside the model, in violation of the above:

| site | violation |
|---|---|
| `validator/invariant.rs::validate` | inlines orbital + electron count equations (~30 lines) |
| `validator/invariant.rs::validate_atom` | duplicates the equations (~30 more lines) |
| `counts.rs::try_build_candidate` + `resolve_unpaired_lone_pairs` | inlines simplified conservation; missing `2·d`, `2·a`, `mc_v`, `I(ar_v)` |
| `counts.rs::candidates_for` | reaches into `element.valence_electrons()`, `element.shift(-charge)`, `self.table.entry(element).allowed_aromatic_valences` |
| `atom_typing.rs::prepare_atom` | reaches into `self.normal_valence.implicit_hydrogens_for(...)`, `view.is_in_aromatic_system()` for aromatic-branch synthesis |
| `atom_typing.rs::resolve` | reaches into `self.registry.lookup(element, charge_key)` directly |

The resolver's simplified form happens to match the full form for the closed-shell main-group cases the conformance suite exercises (where `d = a = mc_v = 0` and `ar_v ∈ {0, 1}` with the table pre-pinning the answer). Off-table inputs — dative donors, multicenter participants, pyrrole-like (`ar_v = 2`) atoms outside the table — diverge silently.

### API surface

```rust
impl ValenceModel {
    /// Resolve candidates consistent with the model + the atom's pinned fields.
    /// AtomTyping: filter registry by matches. Counts: enumerate conservation candidates.
    pub fn resolve(&self, view: &AtomView<'_>) -> Vec<AtomAst>;
    pub fn resolve_atom(&self, atom: &AtomAst) -> Vec<AtomAst>;

    /// Validate a fully-pinned atom against the model.
    /// AtomTyping: atom matches at least one registry pattern.
    /// Counts: conservation equation holds (orbital_count == electron_count).
    pub fn validate(&self, view: &AtomView<'_>) -> Solution<(), ValenceValidationError>;
    pub fn validate_atom(&self, atom: &AtomAst) -> Solution<(), ValenceValidationError>;

    /// Universal physics — orbital/electron count equations, conservation invariant.
    /// Model-independent; both AtomTyping and Counts return the same `Invariants`.
    /// Exposed so the standalone ElectronInvariantValidator can call into it
    /// without going through a specific model's `validate`.
    pub fn invariants(&self) -> &Invariants;
}

pub struct Invariants;  // unit; equations are pure functions of inputs

impl Invariants {
    /// Orbital occupancy from the atom's pinned fields. None if any
    /// required field is non-Lit.
    pub fn orbital_count(&self, view: &AtomView<'_>) -> Option<i64>;
    pub fn orbital_count_atom(&self, atom: &AtomAst) -> Option<i64>;

    /// Electron count from the atom's pinned fields.
    pub fn electron_count(&self, view: &AtomView<'_>) -> Option<i64>;
    pub fn electron_count_atom(&self, atom: &AtomAst) -> Option<i64>;

    /// Conservation invariant: `orbital_count == electron_count` plus the
    /// aromatic_valence range check (∈ {0, 1, 2}). Universal physics.
    pub fn check(&self, view: &AtomView<'_>) -> Solution<(), Mismatch>;
    pub fn check_atom(&self, atom: &AtomAst) -> Solution<(), Mismatch>;

    /// Enumerate candidate atoms consistent with conservation +
    /// element-level bounds + the atom's pinned fields. Used by
    /// `ValenceModel::Counts::resolve` and (as a sanity gate) by
    /// `ValenceModel::AtomTyping::resolve` to invariant-check registry
    /// candidates.
    ///
    /// Output:
    /// - empty → infeasible
    /// - one candidate → unique
    /// - many candidates → underdetermined (caller writes JointDomain)
    pub fn solve(&self, view: &AtomView<'_>) -> Vec<AtomAst>;
    pub fn solve_atom(&self, atom: &AtomAst) -> Vec<AtomAst>;
}
```

The resolver code becomes a one-liner per atom:

```rust
// CountsValenceResolver::resolve and AtomTypingValenceResolver::resolve both become:
for atom in ast.atoms_mut() {
    if atom.is_ground() { continue; }
    match model.resolve(&view).len() {
        0 => return Err(...),
        1 => atom.narrow_from(&candidate),
        _ => attach JointDomain via field-wise meet,
    }
}
```

Same shape for both variants. The two resolver structs (`CountsValenceResolver`, `AtomTypingValenceResolver`) can probably collapse into a single `ValenceResolver` holding only `model: ValenceModel` — variant-specific behavior lives entirely behind `model.resolve()`.

Atom-only vs view-based methods: keep both pairs on both `ValenceModel` and `Invariants`. View-based pulls topology-derived counts (donated/accepted pairs from incident dative bonds, valence from incident bonds); atom-only relies on constraints alone. Same equation core called from both.

`solve` returns the candidate set, not a "uniquely solved" / "underdetermined" enum — the caller (resolver) compares `len()` to decide its action. Pure data return; no embedded policy.

### Element-level bounds tighten the enumeration

Conservation alone often leaves the search space infinite-in-principle (any non-negative integer is a candidate for `h`, `n`, `u`). The element-level bounds from `umol-shared/src/element.rs` close this:

| bound | source | role |
|---|---|---|
| `max_valence(element)` | `element.rs:687` | upper bound on `valence` |
| `charge_bounds(element) -> (min, max)` | `element.rs:693` | `min ≤ charge ≤ max` |
| `max_unpaired_electrons(element)` | `element.rs:699` | upper bound on `unpaired` |
| `max_implicit_hydrogens(element)` | `element.rs:705` | upper bound on `implicit_h` |
| `valence_electrons(element)` | `element.rs:681` | `Z` in the conservation equation |
| `valence_capacity(element)` | **not yet in element.rs** — needs adding | upper bound on orbital occupancy via the shell rules (2 for H/He, 8 for octet rows, 18 for d-block, 26 for f-block) |

`solve` uses these bounds in addition to:

- the atom's pinned field literals (filter candidates to those that match);
- `allowed_valences(element)` and `allowed_aromatic_valences(element)` from the resolver's `ValenceTable` (when in `ValenceModel::Counts`);
- the conservation equation itself (drops candidates that fail).

This makes the enumeration finite by construction and small in practice (the bounds are tight — `max_implicit_hydrogens` ≤ 4 for first-row, `max_unpaired_electrons` ≤ 7).

#### Shell-capacity bound (`total_e ≤ valence_capacity(element)`)

From inv o (doc 52 §10.1.3): `total_e = u + 2n + 2d + 2a + 2h + 2v + ar_v + I(ar_v) + mc_v`. The shell-capacity bound says this sum cannot exceed the element's valence-shell electron capacity:

- 1st row (H, He): 2 (1s only)
- 2nd row (Li–Ne), 3rd row (Na–Ar) under strict octet: 8
- d-block transition metals: 18 (s + 5d + 3p; "18-electron rule")
- f-block (Ce–Lu, Th–Lr): 26 (s + d + f sub-shells; 6p/7p don't bond meaningfully). Admits [Gd(H₂O)₉]³⁺ (orbital count 25 from 7 unpaired f-electrons + 9 dative bonds) and similar high-coordinate hydrates.

For standard main-group atoms with `d = a = mc_v = 0` and conventional valences, **the shell bound is redundant with `max_valence`** — `max_valence` is calibrated to satisfy the octet (e.g., C max_v=4, N max_v=3, O max_v=2 all leave space for `2n` to fill the octet exactly). The conservation equation then pins `n` and the shell bound is automatically satisfied.

The shell bound becomes **load-bearing** in three cases:

1. **Terms not separately bounded.** `lone_pairs`, `donated_pairs`, `accepted_pairs`, `aromatic_valence`, `multicenter_valence` have no per-element max accessors today. A run-away `n`, `d`, or `a` is only caught by conservation (which limits `2n + 2d + 2a` from below by Z − q minus the rest) and by the shell-capacity bound (which limits from above by `shell − 2v − 2h − ...`). For atoms with `q ≪ 0` (deep anions), conservation by itself allows large `n`; shell capacity clamps it.

2. **Transition metals.** `max_valence` for d-block elements is typically loose (Fe: 6, Os: 8, etc.); the 18-electron rule is the real constraint. For metal carbonyl / dative-bond chemistry, accepted_pairs from ligands push orbital occupancy toward 18, and the bound rejects 18+. Without the shell bound the resolver would silently commit invalid states for organometallic inputs.

3. **Hypervalent main-group.** P, S, Cl can carry max_valence ≥ 5 (PCl₅, SF₆, ClO₄⁻). The octet bound is violated by design in these cases — for them the third-row "hypervalent capacity" (12 for S, 14 for Cl) replaces the strict octet. Per-element shell capacity needs to carry the actual capacity, not a blanket "8 for row 2/3".

**Per-element bound table** (reasonably permissive — admits all standard cheminformatics inputs including high-oxidation-state Kekulés like XeO₄; the only edge case it rejects is ClO₄⁻ in the all-double-bond representation, which is a borderline pathological structure):

| group | elements | bound |
|---|---|---|
| Row 1 | H, He | 2 |
| Row 2 main-group | Li–Ne | 8 (strict octet) |
| Row 3 non-hypervalent | Na, Mg, Al, Ar | 8 |
| Row 3 hypervalent | Si, P, S | 12 |
| Row 3 halogen | Cl | 14 |
| Row 4+ groups 13-15 | Ga, Ge, As, In, Sn, Sb, Tl, Pb, Bi | 12 |
| Se (row 4 chalcogen) | Se | 12 (no real compound needs 14) |
| Te, Po (row 5+ chalcogen) | Te, Po | 14 (Te: [TeF₆]²⁻; Po by extrapolation, sparse chemistry) |
| Br (row 4 halogen) | Br | 14 |
| I, At (row 5+ halogen) | I, At | 16 (I: [IF₈]⁻; At by extrapolation) |
| Noble gas (heavier than Ar) | Kr, Xe, Rn | 18 (admits XeO₄, KrF₄) |
| d-block transition metals | Sc–Zn, Y–Cd, La–Hg, Ac–Cn | 18 (18-electron rule) |
| f-block (lanthanides, actinides) | Ce–Lu, Th–Lr | 26 (s + d + f sub-shells) |

Calibration examples:

| compound | center | orbital count | bound used | admitted? |
|---|---|---|---|---|
| PF₆⁻ | P | 12 | 12 | ✓ |
| SF₆ | S | 12 | 12 | ✓ |
| [SeF₆] | Se | 12 | 12 | ✓ |
| **[TeF₆]²⁻** (K₂TeF₆) | **Te** | **14** | **14** | **✓ (per-element bump; S, Se not bumped — no real compound forces it)** |
| IF₇ | I | 14 | 16 | ✓ |
| **[IF₈]⁻** (e.g., [NMe₄]⁺[IF₈]⁻) | **I** | **16** | **16** | **✓ (drives I=16, not 14)** |
| [IO₄]⁻ (all double) | I | 16 | 16 | ✓ |
| XeF₆ | Xe | 14 | 18 | ✓ |
| XeO₄ (all double) | Xe | 16 | 18 | ✓ |
| ClO₄⁻ (3 double + 1 single) | Cl | 14 | 14 | ✓ |
| ClO₄⁻ (all double `O=[Cl⁻](=O)(=O)(=O)`) | Cl | 16 | 14 | ✗ (borderline; rejected by design) |
| [SiF₆]⁴⁻ (hypothetical) | Si | 16 | 12 | ✗ (no known compound) |
| [SnF₆]⁴⁻ (hypothetical) | Sn | 16 | 12 | ✗ (no known compound) |

Per-element-evidence note: bounds are bumped above the row-default only when a real compound forces it. Group uniformity is not enforced — O at 8 (octet) already breaks the chalcogen group regardless, so partial uniformity within a row+group cell is not worth defending at the cost of tightness elsewhere. Po and At are extrapolated from their lighter heavier-row neighbors (Te, I) — their chemistry is sparse enough that strict per-evidence isn't practical, so they inherit the bound of the closest characterized neighbor.

**Required code change**: add `valence_capacity(element) -> u8` to `umol-shared/src/element.rs`. New per-element column in `ELEMENT_DATA`; const accessor. Data table populated from the per-group bounds above.

### Solver strategy: hand-roll + JointDomain

Decision: **hand-rolled per-atom enumeration solver** (nested loops over bounded ranges) inside `Invariants::solve` and `ValenceModel::resolve`, with `AtomConstraint::JointDomain` as the output representation when more than one candidate survives. No external CSP library; no general CSP framework. Rationale:

- **Scale**: per-atom problem has 7–9 bounded integer variables (h, n, u, q, v, ar_v, mc_v, optionally d, a, future haptic_v); product of element-level bounds gives ≤ ~10⁴ pre-filter, ~1–15 post-conservation. Trivial for nested loops; setup cost of an external solver exceeds the solve cost.
- **Stability**: variables won't grow much. Forecast addition is `haptic_valence` (splitting dative bonds into two-center vs multicenter); one more nested range inside `Invariants::solve`, not a new constraint class.
- **Underdetermined output**: candidate sets with >1 element are captured as `AtomConstraint::JointDomain { vars, tuples }` written onto the atom. Preserves the inter-field correlations (e.g., `h + 2n = 4` ties h, n together) that a naive field-wise meet would lose. Downstream resolvers narrow individual fields → JointDomain propagation prunes its tuple list → fixed-point converges.
- **What this is NOT**: a general CSP solver. JointDomain is a table-constraint AST variant. The per-atom solver is hand-rolled. Items 4–5 of doc 96's "Progressive adoption" (cross-atom binds, search/labeling driver) stay out of scope until template/coupled-atom chemistry actually needs them.

`ValenceModel::resolve(view)` output contract (same for both variants):

| candidate set size | `ValenceResolver` action |
|---|---|
| 0 | error (infeasible input) |
| 1 | narrow the atom to the single candidate; commit all fields |
| > 1 | take field-wise meet for shared values; attach `AtomConstraint::JointDomain { vars: <fields with disagreement>, tuples: <enumerated values> }` for the rest |

The `ValenceResolver` shell does not know whether `model.resolve` came from registry filtering (AtomTyping) or conservation enumeration (Counts) — same output shape either way.

Propagation: handled inside `Lattice::meet` via the `saturate` hook (see step 2). Any pass that calls `atom.narrow_from(&new_info)` triggers field-wise meet → saturate → JointDomain pruning automatically. No separate scan or trigger needed; propagation is part of the lattice contract. Resolved JointDomains (tuple count drops to 1) expand back into field commits within the same meet call; contradictions (tuple count drops to 0) surface as `None` from meet. Cascading propagation across multiple JointDomain entries on the same atom runs to fixpoint inside `saturate`.

#### Chemistry-dependent variable ranges live in `ValenceModel::resolve`, not in `Invariants::solve`

`Invariants::solve` iterates each variable over whatever range is already on the atom — `Lit(n)` → fixed, `LitSet([a,b,c])` → iterate, `Undetermined` → fall back to the element-level bound (`max_valence`, `valence_capacity`, etc.). It has no chemistry knowledge of which valences or aromatic-π counts are *allowed* per element.

`ValenceModel::resolve` (both variants) is where that chemistry lives. Before invoking variant-specific resolution (registry filter for AtomTyping, `Invariants::solve` for Counts), it narrows the atom's fields by writing the allowed ranges from `ValenceTable` as standard `LitSet` constraints.

#### Variant-discriminator field handling (aromatic, multicenter)

Some fields are tagged enums whose tag itself carries chemistry meaning: `AromaticValenceAst::{Aromatic(_), NotAromatic, Undetermined}`, `MulticenterValenceAst::{Multicenter(_), NotMulticenter, Undetermined}`. The pre-narrow step uses a three-way decision shared between AtomTyping and Counts (mirror of the existing `is_aromatic` disjunction in `counts.rs:96` and `atom_typing.rs:118`, with the Undetermined case added):

For `aromatic_valence`:

| input atom state | branch decision |
|---|---|
| `view.is_in_aromatic_system()` (idempotency arm — membership from a prior sweep) | activate aromatic — narrow to `Aromatic(LitSet(allowed_aromatic_valences))` |
| `AromaticValenceAst::aromatic(Undetermined).matches(&atom.aromatic_valence)` (first-pass arm — declared `Aromatic(_)` from the parser) | activate aromatic — same narrow |
| `atom.aromatic_valence == NotAromatic` (explicit non-membership) | run non-aromatic only — atom stays with `NotAromatic` |
| `atom.aromatic_valence == Undetermined` (no commitment either way) | **run both branches; union the candidate sets** |

The Undetermined case is the key one: strict narrowing doesn't allow the resolver to *assume* non-aromatic just because membership isn't already established. Both possibilities must remain candidates until something downstream pins.

```rust
fn aromatic_branch_atoms(view: &AtomView<'_>, allowed: &[u8]) -> Vec<AtomAst> {
    let atom = view.ast;
    let allowed_lits: Vec<i64> = allowed.iter().map(|&v| v as i64).collect();

    let activated = view.is_in_aromatic_system()
        || AromaticValenceAst::aromatic(ValueAst::Undetermined)
            .matches(&atom.aromatic_valence);
    let explicit_not = matches!(atom.aromatic_valence, AromaticValenceAst::NotAromatic);

    let aromatic_pinned = || {
        let mut a = atom.clone();
        narrow(&mut a.aromatic_valence,
               AromaticValenceAst::Aromatic(ValueAst::LitSet(allowed_lits.clone())));
        a
    };
    let not_aromatic_pinned = || {
        let mut a = atom.clone();
        narrow(&mut a.aromatic_valence, AromaticValenceAst::NotAromatic);
        a
    };

    match (activated, explicit_not) {
        (true, _) => vec![aromatic_pinned()],
        (_, true) => vec![not_aromatic_pinned()],
        _ => vec![aromatic_pinned(), not_aromatic_pinned()],  // Undetermined: both
    }
}
```

`MulticenterValenceAst` gets the same three-way treatment. Pattern: any tagged-enum field with `{Active(_), NotActive, Undetermined}` shape and chemistry-dependent membership is pre-narrowed by this scheme.

#### Plain-`ValueAst` fields

For fields without a variant discriminator (`valence`, `implicit_hydrogens`, `lone_pairs`, `unpaired`, `charge`, `donated_pairs`, `accepted_pairs`), `ValenceModel::resolve` pre-narrows via plain `LitSet` when there's a chemistry-derived range:

```rust
narrow(&mut atom.valence,
       ValueAst::LitSet(entry.allowed_valences.iter().map(|&v| v as i64).collect()));
```

For variables with no chemistry-derived bound beyond element-level (`unpaired`, `charge`, `lone_pairs`, etc.), the pre-narrow does nothing — `Invariants::solve` falls back to the element-level bound (`max_unpaired_electrons`, `charge_bounds`, `valence_capacity / 2`) when the field is `Undetermined`.

#### Counts and AtomTyping share the pre-narrow step

Both variants call a private `ValenceModel` method that takes a view and returns the list of pre-narrowed atoms (each representing a branch). Per-variant resolution runs against each pre-narrowed atom, results are unioned:

```rust
impl ValenceModel {
    /// Pre-narrow chemistry-dependent ranges. Returns one or more pre-narrowed
    /// atoms (multiple when a variant-discriminator field is Undetermined).
    fn pre_narrow(&self, view: &AtomView<'_>) -> Vec<AtomAst> { ... }

    pub fn resolve(&self, view: &AtomView<'_>) -> Vec<AtomAst> {
        let pre = self.pre_narrow(view);
        match self {
            Self::AtomTyping { registry, .. } => pre.iter()
                .flat_map(|a| filter_registry(registry, a))
                .collect(),
            Self::Counts { .. } => pre.iter()
                .flat_map(|a| self.invariants().solve_atom(a))
                .collect(),
        }
    }
}
```

This keeps the universal physics in `Invariants`, the chemistry-bound setting in `ValenceModel::pre_narrow`, and the variant-specific resolution mechanism (registry filter vs conservation enumeration) in the per-variant dispatch.

### Migration plan

0. **Add `valence_capacity` to `umol-shared/src/element.rs`** — new per-element constant column in `ELEMENT_DATA`, new `valence_capacity(&self) -> u8` const accessor. Per-element values per the table above. Prerequisite for `solve`'s shell-capacity bound. **Done**

1. **Add `Bind` / `Ref` to `ValueAst`** in `umol-ast/src/ast/value.rs`. Mirrors the existing `ElementAst::Bind { id, set } / Ref(id)` shape. Pure additive — existing matches keep working via wildcard arms. Prerequisite for `JointDomain` to reference named field variables. **Done (partial — surface design revised in 2a)**

2a. **AST + DSL cleanup for ElementAst, IsotopeAst, ValueAst.** **Done**
    Before adding `JointDomain` in step 2, the three ASTs need a consistent and principled surface for negation and binds. Three ground rules from the DSL audit:

   - **(i) Same syntax, same AST.** `?h` and `(?h)` must produce identical AST. Parens are transparent at any depth, in all positions where they're optional.
   - **(ii) Isotope cannot route through Value.** Mass numbers are tagged-enum-like, not arithmetic; Isotope gets its own dedicated parser and its own variant set.
   - **(iii) Parens never required.** Anywhere a parenthesized form is accepted, the unparenthesized form parses identically.

   Plus the design conclusions:

   - **Negation as first-class for Element + Isotope.** Their domains are finite Boolean algebras with two compact encodings: positive (Lit/LitSet) and negative (Not/NotSet). All boolean operations close on this representation. Value has no top-level negation; Expr covers it.
   - **Bind with negation.** A bind's admissible domain can be positive or negative. Use a `Polarity` tag rather than a separate enum or per-case variants — keeps the bind struct flat and the lattice ops compact.
   - **Bind/Ref at top level for Value too.** Programmatic construction of `Expr(Var(id))` / `Expr(Mem(Var(id), set))` is canonicalized to `Ref(id)` / `Bind { id, set }` via simplify. Parser produces the canonical forms directly.

   #### Polarity enum

   ```rust
   // umol-ast/src/ast/atom.rs (alongside ElementAst and IsotopeAst)
   #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
   pub enum Polarity {
       Include,  // set is the admissible values
       Exclude,  // set is the excluded values; domain is complement
   }
   ```

   Used by `ElementAst::Bind` and `IsotopeAst::Bind`. ValueAst's `Bind` stays positive-only (`Polarity` not present); negation routes through Expr. If a third user appears later, lift to a shared module then; no value in pre-promoting it.

   #### ElementAst (revised)

   ```rust
   pub enum ElementAst {
       Undetermined,
       Lit(Element),
       LitSet(Vec<Element>),
       Not(Element),                                                      // NEW
       NotSet(Vec<Element>),                                              // NEW
       Bind { id: String, set: Vec<Element>, polarity: Polarity },        // polarity NEW
       Ref(String),
   }
   ```

   Lattice impl extensions:
   - `meet`: 6×6 table over the 5 concrete cases (Lit, LitSet, Not, NotSet, Bind) plus Undetermined/Ref. Mechanical. Examples: `Lit(a) & Not(b)` = `Lit(a)` if `a ≠ b` else None; `Not(a) & Not(b)` = `Not(a)` if `a == b` else `NotSet([a, b])`; `NotSet(s) & NotSet(t)` = `NotSet(union(s, t))`.
   - `join`: dual. `Lit(a) | Not(b)` = Undetermined if `a == b` else `Not(b)`; `NotSet(s) | NotSet(t)` = if `intersection(s, t)` empty then Undetermined else `NotSet(intersection)`.
   - `matches`: superset semantics extended. `Not(x).matches(Lit(y))` = `y ≠ x`; `Not(x).matches(Not(y))` = `x == y`; `NotSet(s).matches(Lit(y))` = `y ∉ s`; etc.
   - `is_ground` / `is_undetermined`: false for Not/NotSet/Bind/Ref.

   Canonicalization in constructors / simplify:
   - `NotSet([])` → `Undetermined` (excluding nothing = full domain).
   - `NotSet([x])` → `Not(x)`.
   - `LitSet([x])` → `Lit(x)`.
   - `LitSet([])` → reject (empty positive = bottom; never stored).
   - Set dedup, first-occurrence preserving.

   #### IsotopeAst (revised — own parser, no Expr)

   ```rust
   pub enum IsotopeAst {
       Undetermined,
       Natural,                                                           // existing
       Lit(i64),
       LitSet(Vec<i64>),
       Not(i64),                                                          // NEW
       NotSet(Vec<i64>),                                                  // NEW
       Bind { id: String, set: Vec<i64>, polarity: Polarity },            // NEW (replaces missing variants)
       Ref(String),                                                       // NEW
   }
   ```

   `Expr` variant **removed**. Isotope mass arithmetic / predicates aren't meaningful in any chemistry use case; routing through Value's Expr was historical accident.

   Lattice impl: same shape as ElementAst (modulo `Natural` handling — meets `Natural & Natural` = `Natural`; `Natural & Lit/anything-else` = None unless we say `Natural` is a separate "channel" from the integer domain).

   `Natural` semantics question: is `Natural` a *value* (specific isotope, the most-abundant one) or a *flag* (use the natural distribution)? Today it's effectively a flag and doesn't interact with literal isotopes. Keep that semantic — `Natural & Lit(n)` = None (different meanings).

   Drop `From<ValueAst> for IsotopeAst` — no longer needed since Isotope doesn't route through Value's parser.

   #### ValueAst (cleanup)

   Variants stay as added in step 1: `Undetermined, Lit, LitSet, Expr, Bind { id, set }, Ref`. **No Polarity on Bind** (negation goes through Expr).

   Parser cleanup:
   - Drop the forced-parens `value_bind` / `value_ref` arms.
   - Add bare-form parsers: `?h` → `Ref("h")`, `?h :: {1,2,3}` → `Bind { id: "h", set: [1,2,3] }`.
   - All forms paren-transparent at any nesting depth.
   - Compound expressions involving `?h` (e.g., `?h + 1`, `?h == 0`) still parse to `Expr(...)` with `Var(id)` inside, unchanged.

   Simplify additions:
   - `Expr(Var(id))` → `Ref(id)`.
   - `Expr(Mem(Var(id), set))` → `Bind { id, set: Box::new(set) }`.
   - Conservative: only when the Expr is exactly that shape; nested Var inside BinOp/Rel stays as Expr.

   Drop the `Bind { id, set }` → `Expr(Mem(Var(id), set))` lowering currently in `From<ValueAst> for IsotopeAst` (the impl goes away entirely with Isotope's dedicated parser).

   #### DSL grammar (shared shape across all three)

   The user-facing syntax aligns:

   | input | Element | Isotope | Value |
   |---|---|---|---|
   | `*` | Undetermined | Undetermined | Undetermined |
   | `H` / `12` / `42` | Lit | Lit | Lit |
   | `{H, F}` / `{12, 13}` / `{1, 2}` | LitSet | LitSet | LitSet |
   | `!H` / `!12` | Not | Not | — (Value uses Expr) |
   | `!{F, Cl}` / `!{12, 13}` | NotSet | NotSet | — |
   | `?e` / `?m` / `?h` | Ref | Ref | Ref |
   | `?e :: H` / `?m :: 12` / `?h :: 2` | Bind (Include, singleton) | Bind (Include) | Bind (positive only) |
   | `?e :: {F, Cl}` / `?m :: {12, 13}` / `?h :: {0, 1}` | Bind (Include) | Bind (Include) | Bind |
   | `?e :: !H` / `?m :: !12` | Bind (Exclude, singleton) | Bind (Exclude) | — |
   | `?e :: !{F, Cl}` / `?m :: !{12, 13}` | Bind (Exclude) | Bind (Exclude) | — |
   | (any of above) wrapped in parens | identical | identical | identical |
   | Isotope-only: `=` | — | Natural | — |
   | Value-only: `?h + 1`, `?h == 0`, etc. | — | — | Expr |

   #### Parser implementation pattern

   For each AST's parser, the structure is:

   ```rust
   fn ast(i: &mut &str) -> PResult<...> {
       alt((
           // Recursive paren wrapper: (X) parses same as X
           delimited('(', ast, ')'),
           // ...all unparenthesized forms
       )).parse_next(i)
   }
   ```

   For the `?id` / `?id :: ...` ambiguity in ValueAst (where Expr also has `?id`-as-Var), the ref/bind arms succeed only when followed by a terminator or `::` (not by an arithmetic/relational operator). Compound expressions fall through to `bool_expr`.

   #### Tests to update

   - Existing isotope tests `case_5_bind` and `case_6_ref_` (in `dsl/atom.rs::test_isotope`): currently expect `Expr(Mem(...))` / `Expr(Var(...))`. Update to expect `Bind { ... }` / `Ref(...)` directly.
   - Existing atom test `case_19_h_bind` (`#h(?h)`): if ImplicitHydrogensAst still routes through Value at this stage, update to expect `Ref` via the From-conversion (or just accept the breakage since ImplicitH is going away).
   - All element_bind / element_ref tests: drop the paren requirement; bare `?e` should also parse to Ref.
   - New round-trip tests across all three ASTs covering: bare/paren forms; Include/Exclude bind polarity; compound Expr forms still produce Expr.

   #### Order of work

   a. Add `Polarity` enum (single trait-adjacent file). **Done** (in `ast/atom.rs`).
   b. ElementAst: add `Not`/`NotSet`, change `Bind` shape, update Lattice, update parser, update tests. **Done** — AST variants, `Polarity` on `Bind`, Lattice impl (polarity-aware `element_set_view` helper), parser arms for `!H`/`!{F, Cl}`/`?e :: !H`/`?e :: !{F, Cl}` (folded in via 2e), and unit tests in place.
   c. IsotopeAst: full rewrite of variants (add `Not`/`NotSet`/`Bind`/`Ref`, drop `Expr`), write dedicated parser, remove `From<ValueAst>`, update lattice + tests. **Done.** `Natural` semantics chosen: own channel — `Natural & Lit = None`, `join(Natural, x) = Undetermined` for `x ≠ Natural`. Dedicated `isotope` parser produces variants directly with no `value`-routing. `IsotopeAst::simplify` deleted (no Expr to lift). `From<ValueAst> for IsotopeAst` deleted. `AtomAst::simplify_values` no longer touches `isotope_mass`.
   d. ValueAst: add simplify normalization rules. **Done** (per 2d). Parser-arm cleanup (drop forced-parens `value_bind` / `value_ref`) handled by 2e.

   Each substep keeps tests green at its end. Substep (c) is the largest — Isotope basically gets rewritten.

   #### Completion criteria

   - `?h` and `(?h)` produce identical AST in all three contexts (Element, Isotope, Value).
   - `?e :: !{F, Cl}` parses to `ElementAst::Bind { id: "e", set: [F, Cl], polarity: Exclude }`.
   - IsotopeAst has its own parser; no `value` routing.
   - `Expr(Var(id))` canonicalizes to `Ref(id)` via simplify.
   - Full DSL round-trip for all forms across all three ASTs.

   #### Dependencies & risk

   - **Dependencies**: none (step 1 partially supersedes itself; this is the cleanup).
   - **Risk**: medium. Touches three AST types and their parsers, plus existing tests. Substantial but self-contained. The `Natural` semantics on IsotopeAst is the one open chemistry call; default to "Natural is its own channel, doesn't meet with Lit" unless someone surfaces a use case.

2b. **Full conversion `ImplicitHydrogensAst` → `ValueAst`.** **Done**
    The doc 96 takeaway "`ImplicitHydrogensAst` can be replaced by `ValueAst` in the future" — that future is now. With the AST cleanup in 2a done, `ValueAst` covers everything `ImplicitHydrogensAst` does, plus the bind/ref machinery for joint constraints, minus the obsolete `Normal` sentinel.

   #### Type-level change

   - Change `AtomAst::implicit_hydrogens` field type from `ImplicitHydrogensAst` to `ValueAst`.
   - Delete `ImplicitHydrogensAst` entirely from `umol-ast/src/ast/atom.rs` along with its impls (Lattice, AsLit, From/Into conversions, `simplify`, etc.).
   - Delete `From<ValueAst> for ImplicitHydrogensAst` and `From<ImplicitHydrogensAst> for ValueAst` — both go.

   #### Defaults

   - Delete `ImplicitHydrogensDefaults` (the per-field defaults config).
   - In `AtomDefaults`, replace the `implicit_hydrogens: ImplicitHydrogensDefaults` field with `implicit_hydrogens: ValueDefaults` (reuse the existing `ValueDefaults` used by `charge`, `lone_pairs`, etc.).
   - All callers / macro expansions / tests that constructed `ImplicitHydrogensDefaults` instances move to `ValueDefaults`.

   #### Parser changes

   - Remove `implicit_hydrogens` parser in `dsl/atom.rs` (currently routes through `value.map(ImplicitHydrogensAst::from)` plus the `=` → Normal arm).
   - Replace the `#h` predicate handler with a direct call to the `value` parser — same as `#c`, `#n`, etc.
   - Drop the `=` → Normal special case. `#h=` is no longer a valid syntax. Use `#h*` (Undetermined) in its place where the intent was "unspecified".
   - The `empty` → `Lit(1)` sugar (`#h` alone meaning one H) stays. That's a Value-side convention now.

   #### Tests

   - Remove `ImplicitHydrogensAst`-specific tests (parser tests, lattice tests, conversion tests).
   - Add `#h` test cases to the value parser test suite. Specifically, the cases currently in `dsl::atom::test_parse_atom` for `h_count`, `h_undetermined`, `h_bind`, `h_set`, `h_expr`, `h_omit` migrate to demonstrate that `#h` now produces `ValueAst` directly.
   - Migrate the existing `case_19_h_bind` test (`C#h(?h)`) to expect `ValueAst::Ref("h")` (per 2a) rather than `Expr(Var)`.
   - The `#h=` → Normal case (currently `case::h_normal` in atom tests) deletes entirely. Replace any inline usage of `#h=` in tests with `#h*` (or with a specific literal if the test intended a specific count).

   #### Migration of `#h=` in existing artifacts

   - **Parsing examples / docs**: review every `#h=` occurrence in `discussion/`, `umol-ast/spec/`, examples, and inline rustdoc. Replace with `#h*` unless redundant (a `#h=` adjacent to e.g. `#u0#c0` was relying on `Normal` semantics; with `Normal` gone, the explicit value-resolver pathway computes the H count, so `#h*` is the correct "let the resolver figure it out" intent).
   - **Conformance test inputs / snapshots**: review `umol-graph/tests/.../*.edn` for `#h=` occurrences. Replace with `#h*` and re-run conformance to confirm the resolved output matches (it should — the semantics for the cases where Normal was working are the cases where `#h*` resolves to the same H count via valence resolution).
   - Snapshots that show resolved H counts shouldn't change (the resolver produces the same result). Snapshots that show pre-resolution form change (`#h=` → `#h*`).

   #### Completion criteria

   - `ImplicitHydrogensAst` is no longer defined anywhere in the workspace.
   - `AtomAst::implicit_hydrogens` is `ValueAst`.
   - `#h` predicate parses to `ValueAst`; `#h=` is rejected by the parser; `#h*` parses to `Undetermined`.
   - Workspace tests pass; conformance suite passes with snapshots updated for the `#h=` → `#h*` surface change.

   #### Dependencies & risk

   - **Dependencies**: 2a (ValueAst cleanup for Bind/Ref); also requires the implicit-H resolution work to already be folded into valence resolution (per the doc 96 takeaway — already the design intent).
   - **Risk**: medium. Touches many call sites because `implicit_hydrogens` is referenced throughout the resolver / validator / DSL / table_ir / etc. But each change is mechanical (`ImplicitHydrogensAst::Lit(n)` → `ValueAst::Lit(n)`, etc.). The `#h=` → `#h*` text substitution in conformance inputs is the largest set of file changes; it's a sed-then-snapshot-update.

2d. **Ensure `ValueAst::simplify` canonicalizes `Expr(Var(id))` → `Ref(id)` and `Expr(Mem(Var(id), set))` → `Bind { id, set }`.** **Done**
    With Bind/Ref as first-class top-level variants, the simplify path is the safety net for programmatic AST construction that produces equivalent Expr-wrapped forms. After simplify, no top-level `Expr(Var)` or `Expr(Mem(Var, set))` should remain — they collapse to their canonical Bind/Ref representations.

   Conservative: the rule fires only when the Expr is *exactly* `Var(id)` or `Mem(Var(id), <literal set>)`. `Expr(BinOp(Var, Add, Lit))` doesn't lift because the outer is BinOp; `Expr(Mem(BinOp(Var, ...), [..]))` doesn't lift because the Mem's first arg isn't a bare Var.

   Add tests covering:
   - `Expr(Var("h"))` simplifies to `Ref("h")`.
   - `Expr(Mem(Var("h"), [1, 2, 3]))` simplifies to `Bind { id: "h", set: [1, 2, 3] }`.
   - `Expr(BinOp(Var("h"), Add, Lit(1)))` is left alone (not a simple shape).
   - `Expr(Mem(BinOp(Var("h"), Add, Lit(1)), [2, 3]))` is left alone.
   - Idempotency: simplify twice on each of the above gives the same result.

   No change needed for `ElementAst` or `IsotopeAst`: neither has an `Expr` variant, so there's no Expr-form-to-canonical-Ref/Bind path. Set canonicalization (`NotSet([x])` → `Not(x)`, etc.) is already handled inside meet/join via the `canonicalize_*` helpers; doesn't need a public `simplify` method.

2c. **Update `umol-ast/spec/umol-dsl-spec.md`** to reflect the surface changes from 2a and 2b. **Done**
    The spec is **normative** for new language constructs — when adding a new form to the DSL, the spec is the binding source of truth for what the parser must accept and emit.

   #### Read the spec first

   Before writing, read `umol-ast/spec/umol-dsl-spec.md` end-to-end. Match the existing conventions for:
   - Section ordering (which AST types appear in what order; how predicates are grouped).
   - EBNF-style vs prose grammar (whichever the spec uses).
   - Worked-example density and shape.
   - Cross-references to other specs (`edn-spec.md`, `opensmiles-spec.md`, `atom-type-spec-query.md`).

   #### What to add

   From 2a:
   - **Negation forms** at the top level for ElementAst and IsotopeAst: `!H`, `!{F, Cl}`, `!12`, `!{12, 13}`. Document the SMARTS-style reading.
   - **Bind with negation** for ElementAst and IsotopeAst: `?e :: !H`, `?e :: !{F, Cl}`, etc. Polarity tag at AST level.
   - **Paren-transparency rule** as a normative grammar rule: anywhere a value-position form is accepted, an arbitrarily-deeply-parenthesized form of it parses identically. Make this an explicit invariant readers can rely on.
   - **Bare bind/ref syntax**: `?id`, `?id :: <domain>` without surrounding parens. Document that `(?id)` and `?id` are equivalent.
   - **ValueAst::Bind / Ref** as first-class top-level forms (parallel to ElementAst).

   From 2b:
   - **Remove `#h*`** from the spec. Remove all references to "Normal" implicit-hydrogens semantics.
   - **Document `#h*`** as the canonical "unspecified, let the resolver compute" form.
   - **`#h` predicate** is now a `ValueAst` predicate, same shape as `#c`, `#n`, etc. — document under the unified value-predicate section if the spec has one.

   #### What to update

   - Any worked examples currently showing `#h*` get rewritten to `#h*`.
   - Any worked examples showing `(?id)`-only parens get a sibling example showing the bare form `?id`.
   - The "implicit hydrogens" section either deletes entirely or merges into the generic value-predicate section.

   #### Cross-reference check

   - `edn-spec.md`: confirm Bind/Ref EDN serialization shape (per 2a, currently EDN strings). Add spec entry if missing.
   - `atom-type-spec-query.md`: check if it references implicit-H forms; update.
   - `opensmiles-spec.md`: SMILES doesn't have binds; nothing to change here unless the doc cross-references the DSL form for implicit H.

   #### Completion criteria

   - Spec compiles (no broken cross-refs); all examples are valid DSL under the post-2a/2b parser.
   - `#h*` does not appear in the spec.
   - Negation, bind, and ref forms are documented for all three AST types where they apply.
   - Paren-transparency is an explicit normative rule.

   #### Dependencies & risk

   - **Dependencies**: 2a and 2b complete.
   - **Risk**: low — documentation update. The risk is missing a cross-reference or example; mitigated by reading the full spec end-to-end before editing.

2. **Add `AtomConstraint::JointDomain(JointDomain)`** in `umol-ast/src/ast/constraint/atom.rs`, with `JointDomain` defined in a new file `umol-ast/src/ast/joint_domain.rs`.

   Decomposed into substeps 2f–2l (see migration plan). Schema (locked — see joint-domain design section below for the full type-shape rationale):

   ```rust
   pub enum JointDomain {
       Undetermined,
       Domain { vars: Vec<JointVar>, tuples: Vec<Vec<JointValue>> },
   }

   #[non_exhaustive]
   pub enum JointVar {
       Charge, ImplicitHydrogens, LonePairs,
       UnpairedElectrons, Multiplicity,             // spin AST fields (Unpaired → UnpairedElectrons rename for clarity)
       Valence, DonatedPairs, AcceptedPairs,
       // Reserved: Element, Isotope, AromaticValence, MulticenterValence, HapticValence
   }

   #[non_exhaustive]
   pub enum JointValue {
       Int(i64),
       // Reserved: Element(Element), Isotope(u32), AromaticValence(AromaticValenceAst), ...
   }
   ```

   **Naming rule**: `JointVar` variants spell the underlying AST field in full (matches the `_count`/`_index` naming convention; no `ImplicitH` shortenings). Spin's `unpaired` field surfaces as `UnpairedElectrons` for read clarity; `Multiplicity` keeps its natural name.

   **Field projection**:
   - `Charge` → `atom.charge`
   - `ImplicitHydrogens` → `atom.implicit_hydrogens`
   - `LonePairs` → `atom.lone_pairs`
   - `UnpairedElectrons` → `atom.spin.unpaired`
   - `Multiplicity` → `atom.spin.multiplicity`
   - `Valence` / `DonatedPairs` / `AcceptedPairs` → `atom.constraints` entries of the corresponding kind

   **Constructor**: `JointDomain::from_ints(vars, tuples: Vec<Vec<i64>>) -> Result<Self, JointDomainError>`. Rejects degenerate inputs:
   - `vars.len() < 1` (zero vars is degenerate)
   - `tuples.len() < 1` (empty tuples is bottom; signaled via `Lattice::meet -> None`, not stored)
   - `tuples[i].len() != vars.len()` for any `i`
   - duplicate vars

   Stored `Domain` values satisfy `vars.len() ≥ 1` and `tuples.len() ≥ 1`, vars sorted and unique, tuples sorted and dedup'd. Internal: wraps each `i64` as `JointValue::Int`. Sibling constructors (`from_mixed`, etc.) added when non-numeric variants land.

   **Lattice impl for `JointDomain`**: standard lattice operations against the `Undetermined` top — see "Lattice behavior" in the joint-domain design section. `meet` is the relational meet (natural join: cartesian / equijoin / intersection); `join` projects to shared vars (or returns `Undetermined` if shared is empty); `matches` is "pattern vars are a subset of target vars and every projected target tuple is in pattern tuples". `is_undetermined` is `true` only for the `Undetermined` variant; `is_ground` is `true` only for single-tuple `Domain` values.

   **Add `Lattice::saturate(&mut self) -> Result<(), Contradiction>`** to the trait with a no-op default impl. New `#[derive(Lattice)]` proc-macro (no derive exists today — net new) generates field-wise `meet`/`join`/`matches`/`is_undetermined`/`is_ground` for struct types, calling `result.saturate()?` at the end of `meet`. Existing hand-rolled struct impls (`AtomAst`, `SpinStateAst`, `BondAst`, `DativeBondAst`, `MulticenterBondAst`, `NoncovalentBondAst`, `AromaticSystemAst`, `NoncovalentBondConstraints`) all migrate to the derive — they're mechanical field-wise compositions today. Enum impls (`ValueAst`, `ElementAst`, `IsotopeAst`) stay hand-rolled (variant-specific rules, not mechanical).

   Types without relational constraints inherit the no-op default `saturate` and pay nothing. Types with relational constraints opt in via an attribute on the derive:

   ```rust
   #[derive(Lattice)]
   #[lattice(saturate = "saturate_atom")]
   pub struct AtomAst { /* fields */ }

   fn saturate_atom(atom: &mut AtomAst) -> Result<(), Contradiction> {
       propagate_joint_domains(atom)
   }
   ```

   The derive's generated `meet` calls `saturate_atom(&mut result)?` at the end. (Rationale: Rust trait specialization isn't stable, so we can't simply override a single trait method on a derived impl. The attribute-on-derive pattern matches `serde` / `derive_builder` conventions.)

   `saturate_atom` walks the constraint container's `JointDomain` entries, prunes tuples against the current field values, and either (i) returns `Err(Contradiction)` if any tuple list becomes empty, (ii) extracts the single remaining tuple's values into the atom's fields (and drops the resolved `JointDomain`) if a tuple list collapses to size 1, or (iii) keeps the pruned `JointDomain` if size stays ≥ 2. Loops to fixpoint — a resolved JointDomain may pin a field that another JointDomain references, cascading.

   **Bottom convention**: empty tuple set is bottom, signaled via `Lattice::meet -> None` per the existing trait convention. Constructor rejects empty `tuples` at construction; meet returns `None` when intersection / projection yields empty. No "transient bottom" workaround needed — the type follows the same `Option`-via-`meet` pattern every other lattice type uses.

   **Substep plan**:

   - **2f**. `JointDomain` type (enum with `Undetermined` and `Domain` variants), `JointVar` / `JointValue` enums, `from_ints` constructor with all rejections + canonicalization. Unit tests for constructor invariants and equality.
   - **2g**. `AtomConstraint::JointDomain(JointDomain)` variant wiring (kind, simplify routing, container ops, `is_unique = false` so multiple JointDomains can coexist).
   - **2h**. Hand-rolled `Lattice` impl on `JointDomain` itself: relational `meet` (natural join), `join` (projection-to-shared-or-top), `matches` (projection rule). Total operations against the `Undetermined` top — no special-case fallbacks.
   - **2i**. `#[derive(Lattice)]` proc-macro + migration of the eight struct impls listed above. No saturate hook yet (derive emits `meet` body without the saturate call). Separate design discussion before coding.
   - **2j**. `Lattice::saturate` trait method (default no-op) + `saturate_atom` implementation (cross-field propagation; optional 1-tuple/1-var canonicalizations) + derive macro hook (`#[lattice(saturate = "…")]`).
   - **2k**. DSL syntax + parser for `#E(?v1,…,?vn) :: {(l1,…), …}` per the "Notation" section above. Includes `JointDomain::from_ints` integration so parse errors surface constructor failures.
   - **2l**. EDN serialization roundtrip.

2e. **Remove forced-parens requirement on bind/ref in DSL parsers.** Per the audit ground rules, the parser must emit `Bind` / `Ref` directly for bare and parenthesized surface forms, never generating `Expr(Var(_))` or `Expr(Mem(Var(_), _))` shapes for these. Parens become transparent at any nesting depth: `?h`, `(?h)`, `((?h))` all produce identical AST. Compound expressions involving `?h` (e.g., `?h + 1`, `?h == 0`) continue to route through `Expr` because they genuinely need expression structure.

   This step also subsumes the **parser-side work from 2a(b)** for ElementAst: bare-form `?e` / `?e :: {set}` plus negation (`!H`, `!{F, Cl}`, `?e :: !H`, `?e :: !{F, Cl}`).

   #### Disambiguation rule

   Inside the `value` parser's `alt(...)`, the `value_bind` / `value_ref` arms succeed only when followed by a `terminator` (end-of-input or `#`). This is the existing pattern used for `signed_int` — extending it lets `(?id :: {set})` be intercepted as `Bind` when it's the entire predicate body, but fall through to `bool_expr` when it appears as a parenthesized operand inside a larger expression (e.g., `(?a :: {0}) & 0 <= 0`). The narrow terminator fix has already been applied to `value_bind` / `value_ref` (existing fns) but the bare-form variants remain to be added.

   #### What changes

   - **ValueAst parser**: add bare-form arms for `?id` (Ref) and `?id :: {set}` (Bind) that emit those variants directly, with the same terminator check. Drop the forced-parens-only restriction.
   - **ElementAst parser**: add bare-form arms for `?e` (Ref), `?e :: {set}` (Bind, Include), `!H` (Not), `!{F, Cl}` (NotSet), `?e :: !H` (Bind, Exclude, singleton), `?e :: !{F, Cl}` (Bind, Exclude, multi). Existing `element_bind`/`element_ref` paren-only arms become a subset of the new bare-form parsers.
   - **IsotopeAst parser**: deferred to 2a(c) when Isotope gets its dedicated parser.
   - **Display impls**: render in the canonical bare form (no outer parens) — `?id`, `?id :: {set}`, `!H`, etc. The parser accepts paren-wrapped equivalents but Display picks one canonical shape.

   #### Test migrations (pre-planned, in-scope)

   These are explicit migrations that are part of this step's plan, not "test edits to fix code":

   - `dsl::value::tests::test_value`: cases that currently expect `Expr(Var("h"))` for `?h` migrate to expect `ValueAst::Ref("h")`. Cases that expect `Expr(Mem(Var("h"), set))` for `?h :: {set}` migrate to `ValueAst::Bind { id: "h", set }`. The `?h + 1`, `?h == 0` etc. cases stay as Expr (no migration).
   - `dsl::atom::tests::test_parse_atom::case::h_set` migrates accordingly (`N#h?h :: {2,3}` → `ValueAst::Bind { id: "h", set: [2,3] }`).
   - `dsl::atom::tests::test_element` for ElementAst: add cases for `!H` / `!{F, Cl}` / `?e :: !H` / `?e :: !{F, Cl}`; add bare-form cases for `?e` / `?e :: {C, N}` paralleling the existing paren-wrapped ones.
   - Add paren-transparency tests at one and two levels of nesting for each AST type.

   #### Completion criteria

   - `?h` (bare) and `(?h)` and `((?h))` all parse to `ValueAst::Ref("h")` in standalone value position.
   - `?h :: {set}` (bare) and `(?h :: {set})` and `((?h :: {set}))` all parse to `ValueAst::Bind { ... }`.
   - `(?a :: {0}) & 0 <= 0` parses to `Expr(And([Mem(Var("a"), [0]), Rel(...)]))` (the parens are interpreted as Expr grouping, not as a Bind intercept).
   - Element parser accepts `!H` and `!{F, Cl}` at top level and `?e :: !H` / `?e :: !{F, Cl}` inside Bind.
   - All workspace + property tests pass post-migration.

   #### Dependencies & risk

   - **Dependencies**: 2a (Polarity + AST variants), 2d (simplify rules — though with parser emitting directly, simplify becomes the safety net for programmatic construction only).
   - **Risk**: medium. The disambiguation rule is subtle (terminator-after-paren); some Expr-shaped tests will need pre-planned migration to the new canonical Bind/Ref expectations. Care needed to avoid the conflation between "I'm changing the test because code changed" (forbidden) and "this test migration is in the plan" (in-scope). Each test edit in this step must point back to this section.

3. **Add `Invariants` module** (`umol-graph/src/ops/valence/invariants.rs`) with `orbital_count`, `electron_count`, `check{,_atom}`, `solve{,_atom}`. Equations transcribed verbatim from `validator/invariant.rs::validate` (matches doc 52 §10.1.3). `solve` is the hand-rolled nested-loop enumerator using `valence_capacity`, `max_valence`, `charge_bounds`, `max_unpaired_electrons`, `max_implicit_hydrogens` as range bounds.

4. **`ValenceModel` API methods** (`resolve`, `resolve_atom`, `validate`, `validate_atom`, `invariants`) on the enum in `umol-graph/src/ops/config.rs`. Per-variant dispatch:
   - `AtomTyping::resolve` — move `prepare_atom` here (constraint synthesis, registry filter, optional invariant sanity gate). No `normal_valence` pre-fill.
   - `AtomTyping::validate` — registry-membership check: `prepared.matches(pat)` for at least one `pat` in `registry.lookup(element, charge)`.
   - `Counts::resolve` — call `self.invariants().solve(view)`. Bounded by `allowed_valences` / `allowed_aromatic_valences` from the embedded `ValenceTable`.
   - `Counts::validate` — call `self.invariants().check(view)`.
   - `invariants()` returns the same `Invariants` value for both variants — universal physics.

   All access to `AtomTypeRegistry`, `ValenceTable`, `element.valence_electrons()`, `element.shift(...)`, etc. happens *inside* these methods. No `pub` accessor exposes raw registry / table state — the model API is the boundary.

5. **Collapse resolvers into one thin shell.** `CountsValenceResolver` and `AtomTypingValenceResolver` both become `ValenceResolver { model: ValenceModel }`. Per-atom body:
   ```rust
   for atom in ast.atoms_mut() {
       if atom.is_ground() { continue; }
       match self.model.resolve(&view).as_slice() {
           [] => return Err(ValenceError::NoValidState { atom: id, ... }),
           [cand] => atom.narrow_from(cand),
           candidates => attach_joint_domain(atom, candidates),
       }
   }
   ```
   No model-specific code lives in this loop. `try_build_candidate`, `resolve_unpaired_lone_pairs`, `prepare_atom`, `candidates_for`, `build_aromatic_candidates` all move into `ValenceModel` (per step 4). Both resolver files (`counts.rs`, `atom_typing.rs`) are deleted; `resolver/valence.rs` becomes the single entry point.

6. **Collapse validators similarly.** `validator/invariant.rs::validate{,_atom}` becomes a thin caller of `model.invariants().check{,_atom}(view)` — universal-physics validator, runs regardless of `ValenceModel` variant. A new `validator/valence.rs` calls `model.validate(view)` — model-specific (registry membership for AtomTyping, conservation for Counts). Equations deleted from validator code; both files are pure iteration shells.

7. **Cleanup: remove `NormalValenceTable` and its config plumbing.** Counts no longer consults it (step 4 routes counts through `Invariants::solve`); atom-typing no longer pre-fills `h` from it (registry patterns carry pinned `h`). Drop the type, the TOML, and the `ValenceModel::AtomTyping::normal_valence` / `ValenceModel::Counts::normal_valence` fields.

   `ValenceTable` stays — it provides `allowed_valences` and `allowed_aromatic_valences` per element, which `Counts::resolve` still needs to bound the enumeration. `AtomTypeRegistry` stays — it's the atom-typing source of truth.

Each step lands independently with tests green. Steps 1–2 are pure additions; steps 3–4 add new infrastructure; steps 5–8 swap call sites one at a time.

### Tracked cleanups (small, opportunistic)

These are small follow-ups noticed during the in-progress work. Not blocking any other step; pick up alongside the nearest related substep.

- **`umol-graph/src/ops/valence/counts.rs:181`** — unreachable-pattern warning. The old `ImplicitHydrogensAst::Normal` arm became unreachable when 2b merged it with `Undetermined`. One-line tidy.
- **`umol-graph/src/table_ir/lift.rs:135`** — `Some(TableImplicitH::Normal) => ValueAst::Undetermined`. With `Normal` retired at the AST level, the table-IR variant `TableImplicitH::Normal` is now load-bearing only for the SMILES/MOL parser default. Review whether the variant should also be retired (and parsers map directly to `Undetermined` at the table layer) or if the distinction still carries useful information for the lift step.
- **`NormalValenceTable` last consumer** — `atom_typing.rs::prepare_atom:124` still pre-narrows `implicit_hydrogens` via `normal_valence.implicit_hydrogens_for(...)`. Step 7 plans removal; flagging that this is the final consumer to retire.

### Open design choices (pick before implementing)

1. **`Invariants` as unit type vs owned per-model.** The equations are pure functions; a unit type works. If a future model carries equation-altering policy (e.g., "ignore multicenter contributions"), the unit becomes a config-bearing struct. Default: unit type now, expand if needed.
2. **`solve` cost when many free vars.** A fully open atom (`Default::new()` only) has 4–5 free vars × per-element bound products — at most low hundreds, but worth a perf check on the conformance suite. If hot, add fast paths (single-unknown closed-form) inside `solve`.
3. **`Mismatch` payload.** Carry both `orbital_count` and `electron_count` (matches today's `ElectronInvariantContradiction::AtomInvariantMismatch`) plus the aromatic-bounds case as a separate variant; or one variant with a free-form `reason: String`. Prefer typed variants.
4. **Aromatic-valence range check.** Include in `check` (per doc 52 fn 8) or split into `check_aromatic_bounds`? Suggest include — it's part of the same per-atom invariant; splitting fragments error reporting.
5. **`JointDomain` tuple-list ordering.** **Decided: set semantics (unordered).** See "Things to nail down" above for the WCSP alternative noted but not adopted.
6. **`JointDomain` round-trip EDN/DSL syntax.** **Decided**: `#E(?v1,…,?vn) :: {(l1,…), …}` — `{}` for the set (matches existing LitSet/element-set), `()` per tuple, comma at both levels, `::` for membership (matches bind-domain). See "Notation" in the joint-domain design section. Parser + EDN roundtrip land in substeps 2k–2l.

## Takeaways

- **`Normal` is not useful for implicit hydrogens.** The sentinel goes away.
- **No separate implicit-H perception phase.** Subsumed into valence resolution (atom typing / counts).
- **SMILES bare-atom semantic**: `#h* #u0`. The "normal valence" rule was effectively the sum formula with `u` left out — equivalent to `u = 0` by default. Defensible because all `u > 0` states require bracket annotation with explicit H (e.g., `[NH2]` = `N #h2 #c0`; charge is always explicit in bracket atoms). Everything else can be filled in.
- **MOL has no bare-atom distinction.** Add an "add hydrogens" parsing flag. Set → emit `N #c0 #h*` (or whatever explicit atom properties are present). Unset → emit `N #c0 #h0`.
- **Implicit ↔ explicit H transformations are still needed.** Mechanical, model-independent: `{:atoms ["C#h4"] :bonds []}` ↔ `{:atoms ["C" "H" "H" "H" "H"] :bonds [[0 1 :single] [0 2 :single] [0 3 :single] [0 4 :single]]}`.
- **Valence resolution runs before aromaticity.** Uses atom typing where needed. Counts needs to include charge + unpaired electrons when computing implicit H. May not resolve some valence states without specific evidence; that's acceptable.
- **Fixed-point iteration is fundamentally useful** (e.g., aromaticity charge equalization). Validation can catch the affected cases for now.
- **True disjunctions in AtomAst are missing.** Cases like Fe²⁺ — `Fe#c+2#n3#u0 | Fe#c+2#n1#u4` — aren't easily expressible today. The previous smallvec-of-candidates was a disjunction in disguise; the AST rewrite removed it. Four shapes were considered (top-level DNF; per-field domains + relational constraints; atom-type-id indirection into an inlined registry; hierarchical ADT with stable-prefix / varying-suffix). The latter two have structural problems: atom-type indirection makes resolution a different operation (promote `Set<TypeId>` to `Single`) rather than uniform field narrowing, and the hierarchy approach requires picking a global outer/inner split that doesn't generalize (different element classes vary along different axes — `u/n` for radicals and TM spin, `h/v` for hypervalent main-group, `v/a` for aromatic donation alternatives, `h/lp` for tautomers).
- **Per-field domains + relational constraints (CSP-style) is the most reasonable tradeoff.** Each field keeps its existing AST shape and may carry a domain (`LitSet`, `Set`, etc.); inter-field correlations live as new constraint variants. Rough DSL shape for the Fe²⁺ joint-domain example: `Fe#c+2#n?n#u?u#E(?n,?u)=[(3,0);(1,4)]`, where `#E` introduces a domain over a tuple of named variables (`?n`, `?u`) and lists the allowed combinations. Syntax details TBD. Open question: whether the joint-domain operator applies only to numerical fields or also to the element field and other non-numeric AST components.
- **`ImplicitHydrogensAst` can be replaced by `ValueAst` in the future.** For now, treat `Normal` operationally identical to `Undetermined`.

## Joint-domain constraint design (per-field + `#E`)

Notes on the per-field domains + relational constraints (CSP-style) approach.

### Notation

```
Fe#c+2#n?n#u?u#E(?n,?u) :: {(3,0), (1,4)}
```

- `?n`, `?u` — named field-level binds. Extend the existing `ElementAst::Bind { id, set }` / `Ref(id)` machinery to `ValueAst` so any field can declare or reference a named variable.
- `#E(vars) :: {tuples}` — joint-domain constraint: the tuple of named binds must be drawn from the listed set. Set semantics (unordered); same `::` glyph as bind-domain membership (`?h :: {1,2,3}`).
- Inner tuples: `(v1, v2)`. Outer set: `{(...), (...)}`. Comma at both levels — `()` brackets each tuple, so the outer comma is unambiguous.
- Tuple values are always lits (no expressions, no nested sets), so neither delimiter can be confused with `value-expr`.

### Type shape

`JointDomain` is an enum with the standard top variant — matching the convention every other lattice type in `umol-ast` follows (`AromaticValenceAst::Undetermined`, `MulticenterValenceAst::Undetermined`, `ValueAst::Undetermined`, …):

```rust
pub enum JointDomain {
    Undetermined,                                                  // top
    Domain { vars: Vec<JointVar>, tuples: Vec<Vec<JointValue>> },  // proper constraint
}
```

### Invariants

The `Domain` variant satisfies:

- `vars.len() ≥ 1` (zero is genuinely degenerate)
- `tuples.len() ≥ 1` (empty is bottom — signaled via the existing `Lattice::meet -> Option<Self>` convention, not stored)
- `tuples[i].len() == vars.len()` for all `i`
- `vars` is sorted and unique; `tuples` is sorted and dedup'd

`from_ints` rejects degenerate inputs (zero vars, zero tuples, arity mismatch, duplicate vars), canonicalizes (sort vars + permute tuples + sort + dedup), and returns `Err` on any violation — the parser surfaces these as parse errors.

A `Domain { tuples: [t] }` is a ground state (every var pinned by the single tuple). It is *valid storage*, not a transient invariant violation. A `Domain { vars: [single] }` is redundant with a per-field constraint on `single` but harmless — saturate may normalize it to the per-field form (an optional canonicalization), but the type system doesn't require it.

### Lattice behavior

Every operation is total and clean — no special-case workarounds:

- `is_undetermined`: `matches!(self, Self::Undetermined)`
- `is_ground`: `matches!(self, Self::Domain { tuples, .. } if tuples.len() == 1)`
- `meet(Undetermined, x) = Some(x)`; `meet(Domain, Domain)` is the relational meet (natural join: cartesian product on disjoint vars, equijoin on shared, intersection on identical) — returns `None` only when the result tuple set is empty (genuine contradiction).
- `join(Undetermined, _) = Undetermined`; `join(Domain, Domain)` returns `Undetermined` cleanly when shared vars are empty; otherwise projects to shared vars, unions, dedups, wraps as `Domain`.
- `matches(Undetermined, _) = true`; `Domain` pattern matches `Domain` target iff `pattern.vars ⊆ target.vars` and every target tuple projected to `pattern.vars` is in `pattern.tuples`.

### Simplify vs saturate

Two responsibilities, two operations:

- **`simplify`** is infallible per-type canonical-form normalization — idempotent. For `JointDomain`: re-sort and dedup `tuples`; vars are already canonical from the constructor.
- **`saturate`** (new `Lattice::saturate(&mut self) -> Result<(), Contradiction>` in step 2j) is the *cross-field propagation* step — projects per-field constraints into JointDomain tuples and vice versa, propagates pinned values across the constraint container. Returns `Err(Contradiction)` if propagation reveals incompatibility; that propagates to `None` at the `AtomAst::meet` boundary.

Optional saturate canonicalization (not load-bearing for correctness):

- Single-tuple `Domain` → pin each var to its tuple value, drop the JointDomain. Equivalent representation; saves storage.
- Single-var `Domain` → fold into the corresponding per-field constraint as a `Set`. Equivalent representation.

These are optimizations, not invariants. The lattice is well-defined either way.

### Things to nail down

- **Order semantics of the tuple list.** **Decided: set semantics (unordered).** Meet is intersection — commutative, associative, idempotent, monotone. Tuples normalized (sort + dedup) inside the constructor and after each meet so equal sets compare equal. *Alternative considered: weighted CSP (soft constraints)* — each tuple carries a cost; meet adds costs; narrowing shrinks to min-cost tuples; solving becomes optimization. Principled but materially heavier — propagators turn cost-aware, the solver becomes an optimizer. Not adopted: out of scope for the current narrowing pipeline. Revisit if quantitative ranking ever has to compose under meet (vs. being applied at the search-driver level over candidate sets, where the lattice doesn't need to know).
- **`#E` as one or many joint-constraint kinds.** Product-domain is one shape. Others: equality (`?n = f(?u)`), linear (`?n + ?u ≤ k`), single-variable set-membership (`?n ∈ {1,3}`). `ValueAst::Expr` already covers some of this. Decide whether `#E` is its own kind or whether broader Bind+Expr machinery covers product-domain via conjunction (`?n ∈ {1,3} & ?u ∈ {0,4} & (?n=3 ↔ ?u=0)`). Product-domain form is more compact and reads chemistry-naturally; Expr form is more uniform with existing AST.
- **Bind-name scoping.** Atom-local (atom-id-qualified) or molecule-scope (cross-atom referenceable)? Existing `ElementAst::Ref` supports cross-atom for elements. Generalizing uniformly to values is consistent; whether the chemistry needs cross-atom for values (e.g., coupled diradical pairs) is a separate question.

### Applicability beyond numeric fields

The natural targets are numeric fields (`h`, `q`, `u`, `n`/`lp`, `v`) plus `element`. Element is already half-supported via `ElementAst::Bind/Ref` — extending the joint-constraint machinery to include it costs nothing extra. Concrete use cases:

- **Element + numeric joint** (`[Fe; Co]#E(?el,?u)=[(Fe,4);(Co,3)]`): metalloprotein active sites where metal identity and spin state are coupled.
- **Charge + numeric joint** (as in the Fe²⁺ example).
- **Isotope**: niche (isotope-labeling experiments); probably not worth designing for upfront.
- **Bond-order fields on incident bonds**: cross-entity, different scope. Not in this design.

The mechanism doesn't need to be restricted to numeric, but the practical targets are numeric + element.

### Implementation surface

- Constraint lives in the atom's `constraints` container as a new variant carrying `binds: Vec<BindId>` and `tuples: Vec<Vec<Value>>` (or a more general expression form).
- Resolver gains an arc-consistency propagation step: narrowing one bind's domain via any source (R1 pinning, conservation, validator-driven contradiction) re-projects the joint constraint onto the others. Standard CSP propagation. Cheap on tables this small.
- Compose-ability: multiple joint constraints on the same atom AND-conjoin; their meet narrows the joint domain. Single-variable domain constraints (`LitSet` on `ValueAst`) compose with joint constraints uniformly.
- Round-trip: the named binds must round-trip through EDN. Bind-name uniqueness and scoping decisions feed back into the parser/serializer.

### Precedent — this is CSP

Largely, yes. Per-field domains + relational constraints + arc-consistency narrowing + (eventually) labeling/search-with-ranking is exactly classical constraint-satisfaction. The lattice + `narrow_from` machinery already in place is arc consistency in lattice clothing. Naming the relevant literature so we don't re-derive:

- **Mackworth's arc consistency** (AC-3, AC-4) — the original algorithms for finite-domain CSP propagation.
- **Constraint Logic Programming over finite domains (CLP(FD))** — Prolog dialects (SICStus, SWI-Prolog) with `library(clpfd)`; ECLiPSe-CLP. The exact "name your variables, declare their domains, post relational constraints, solve" surface.
- **MiniZinc** — modern declarative constraint-modeling language; clean syntax for product-domain constraints (`table([?n, ?u], [|3, 0|, 1, 4|])`).
- **Saraswat's Concurrent Constraint Programming** — the lattice-theoretic formulation; constraints as monotone narrowing on a meet-semilattice. Maps directly to umol's `Lattice` trait + `narrow_from`.
- **Apt's *Principles of Constraint Programming*** (Cambridge, 2003) — textbook coverage.
- **Cheminformatics-internal precedent**: SMARTS/SMIRKS pattern matching is implicit CSP; RDKit's substructure matcher (VF2-derived) is a specialized CSP solver. The atom-typing problem itself is finite-domain CSP. Worth mapping umol's vocabulary onto the CSP one explicitly so the cheminformatics literature is reachable rather than re-derived.

The relevant question isn't "are we rediscovering CSP" (we are) but "do we want to embrace the CSP framing explicitly?" Embracing it means:

- Adopting the vocabulary (variable, domain, constraint, propagator, labeling, search) consistently in code and docs.
- Possibly using or porting an existing finite-domain solver instead of building one (Rust ecosystem candidates exist but maturity varies — `copper`, `pumpkin`; explicit evaluation needed if this route is taken).
- Getting termination / correctness guarantees and standard algorithms (AC-3, GAC-Schema, MAC, conflict-directed backjumping) for free from the literature.

Not embracing it means continuing the lattice-and-narrowing-from-first-principles approach. Algorithmically equivalent; vocabulary divergent.

### Progressive adoption

Progressive addition is possible — the current code already implements the algorithmic core of CSP under lattice vocabulary; adding the missing features is additive, no restructure required.

What's already in place that maps onto CSP:

| CSP concept | Current umol equivalent |
|---|---|
| variable | AST field (e.g., `atom.charge`, `atom.implicit_hydrogens`) |
| domain | `LitSet` for ValueAst, `Set` for ElementAst, `Lit(n)` for singleton, `Undetermined` for full domain |
| propagator | `Lattice::narrow_from` |
| arc consistency | meet-driven narrowing already in resolver passes |
| named variable / first-order var | `ElementAst::Bind { id, set }` + `Ref(id)` |
| cross-variable constraint | per-entity `*Constraints` containers (single-field constraints only today) |
| labeling / search | not yet — resolver enumerates candidates in counts/atom-typing but doesn't have a general search driver |
| fixed-point engine | not yet — resolver is a single topological sweep |

What's missing and additive (each item independent, no caller breakage, no renames):

1. **`Bind` / `Ref` on `ValueAst`** (and other value-bearing AST types). Mirrors what `ElementAst` already has.
2. **Joint-domain constraint variant** (e.g., `AtomConstraint::JointDomain { binds, tuples }`). New variant on the existing enum; existing `meet`/`matches` for the constraint container handles it field-wise.
3. **Propagation step** in the resolver. When one bind narrows, project the joint constraint onto the others' domains. Small loop inside the existing resolver pipeline.
4. **Cross-atom bind scoping** (if needed for templates with coupled variables). `ElementAst::Ref` already covers this for elements; generalize uniformly.
5. **Optional later**: search/labeling driver for cases where propagation alone doesn't determine unique values. Standalone resolver pass.

What would require restructure and is **not** necessary for the progressive path:

- **Codebase-wide CSP vocabulary rename** (`Lattice` → something CSP-ish, `narrow_from` → `propagate`, `LitSet` → `Domain`, etc.). Pure cosmetic, large diff, no algorithmic gain. The lattice vocabulary is mathematically equivalent and already established.
- **Switching to an external CSP solver** (`copper`, `pumpkin`, MiniZinc/CHR port). Big API surface change at the resolver boundary; loses chemistry-specific lattice integration; only worth considering if inference problems get large enough that a generic solver outperforms the bespoke one — they won't for atom-scale problems.

Vocabulary coexistence is fine: the lattice viewpoint and the CSP viewpoint are two vocabularies for the same algorithm. Code stays lattice; docs can map to CSP when useful for connecting to literature.

Pacing: each of the five additive pieces can land independently, in any order, with tests passing throughout. There's no half-done state where the codebase breaks — the existing lattice machinery works for field-level cases today; new pieces extend it without replacing it.
