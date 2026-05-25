# Counts model on top of invariants

Status: Spin model revised. The original min-u / max-n rule (decisions 1, 5)
is **superseded** — it is anti-Hund and gives wrong open-shell ground states.
Current model: `covalence_set` fixes `h`, Hund's first rule fixes `u`. See
"Spin resolution" below.

## Decisions

| #   | Decision                                                                                              |
| --- | ----------------------------------------------------------------------------------------------------- |
| 1   | ~~Counts = invariant enumeration + min-u + max-n.~~ **Superseded** — see decision 7.                  |
| 2   | `ValenceTable` field renames only: `allowed_valences` → `covalence_set`, `allowed_aromatic_valences` → `aromatic_valence_set`. |
| 3   | `NormalValenceTable` removed. Doc 96 step 7 folds into this redesign.                                  |
| 4   | Counts and AtomTyping coexist. AtomTyping covers multicenter / dative / specific Lewis modes Counts cannot express. |
| 5   | ~~Min-u and max-n apply per atom.~~ **Superseded** — see decision 7.                                  |
| 6   | `ValenceInvariants::solve_atom` returns `Vec<AtomAst>` (ground candidates); the caller selects.        |
| 7   | **Spin model**: `covalence_set` (charge-shifted) fixes `h`; **Hund's first rule** over the octet open shell fixes `u`; `n` follows from conservation. No min-u / max-n. |
| 8   | **Element scope**: main group only. d/f-block `u` is left **undetermined** (free-atom multiplicity does not determine molecular spin — that is oxidation-state / ligand-field dependent). |
| 9   | Hund open-shell capacity is the **octet** (`s`+`p` = 4 orbitals), decoupled from `valence_capacity` (which is bonding capacity). |
| 10  | Charge enters `covalence_set` lookup as an **isoelectronic element shift** (`element.shift(−charge)`), not an arithmetic adjustment — charge is sign-unconstrained and cannot be a slack variable. |
| 11  | **Aromatic branch is not a separate branch.** Aromatic iff `is_in_aromatic_system()` **or** `Aromatic(Undetermined)` matches the `aromatic_valence` constraint (see "Aromaticity criterion"). For aromatic atoms `av` enumerates over `aromatic_valence_set`; `covalence_set` constrains `v + h + ai` (`ai = av==1`), shared with non-aromatic. The reserved π orbital is keyed on that criterion (any `av`, including 0), never on membership alone or `av>0`. |
| 12  | **Covalence** (Langmuir, *J. Am. Chem. Soc.* **1919**, *41*, 868; def. in doc 52 §10.1.4) = ordinary covalent bonds formed from the atom's own electrons = `v + h + ai`. `#d` (donated), `av=2`, and `#m` (multicenter) are not covalence (the atom donates/delocalizes, not 1:1 sharing). Charge **and** accepted dative pairs are also not covalence — they enter as **isoelectronic shifts** of the `covalence_set` lookup (decision 10 extended; shift `= 2·accepted − charge`: `−1` charge → `+1`, an accepted pair → `+2`). So `covalence ∈ covalence_set(element.shift(2·accepted − charge))` uniformly. Exposed as `AtomView::covalence()` (= the checked quantity; no separate "gained" form). `#d0 #t0` in the counts resolver is a temporary scope limit, not part of the definition. Distinct from `total_valence` (`v+h+av+mv`). |
| 13  | Charged aromatic ions resolve **directly**, no equalization: both `covalence_set` and `aromatic_valence_set` are charge-shifted, so Cp⁻ `C⁻ ≅ N` (`av=2` donor) and tropylium `C⁺ ≅ B` (`av=0` acceptor). Charge equalization is output canonicalization, not a resolvability requirement. |

## Spin resolution

Non-aromatic branch assumes `#a! #m! #d0 #t0` (aromatic / multicenter / dative
contributions zero), so `AtomView::total_valence` collapses to `v + h`
(topology bond order + implicit H) and the conservation equation is

```
e − c = v + h + u + 2·n
```

`e` = `valence_electrons`, `c` = charge, `v` = topology valence, `h` = implicit
hydrogens, `u` = unpaired electrons, `n` = lone pairs.

Two independent jobs:

### h — valence saturation (`covalence_set`)

If `h` is free (implicit), fill it so `v + h ∈ covalence_set(element.shift(−charge))`.
Charge is an isoelectronic shift: `O⁻` reads `F`'s `covalence_set = [1]` (so a
carboxylate oxygen saturates at `v+h = 1`, not the neutral `2`); `N⁺` reads
`C`'s `[4]`. If `h` is pinned (explicit, e.g. a radical with a fixed H count),
`covalence_set` does not fill — the given count stands. `covalence_set` only
*narrows*; an empty entry means "element default, no narrowing."

This is why a separate `NormalValenceTable` is unnecessary: the charge-shifted
`covalence_set` lookup is the entire valence-saturation rule.

### u — spin (Hund's first rule)

`u` is **not** minimized. It is Hund's first rule over the octet open shell
(`s` + `p` = 4 orbitals; `valence_capacity` is *not* used — it is bonding
capacity, and would wrongly give e.g. `PH₃` three open `p` orbitals). Let
`R = e − c − (v+h)` be the non-bonding electrons (`R = u + 2n`). Bonds consume
orbitals `s`-first:

- **bonded** (`v+h ≥ 1`): the `s` orbital is in a bond/hybrid, so the
  remaining `K = 4 − (v+h)` orbitals are degenerate `p` → `u = Hund(R, K)`
  where `Hund(x, k) = x` if `x ≤ k`, else `2k − x`.
- **bare** (`v+h = 0`): the `s` orbital is non-bonding and lowest, so `s²`
  pairs first, then Hund over the 3 `p` orbitals: `u = Hund(R − 2, 3)` for
  `R ≥ 2`, else `u = R`.

`n = (R − u) / 2`.

### Element scope

Applies to main-group elements. For d/f-block, `u` is left undetermined: a
bare metal's Hund ground state (Fe ⁵D, Cr ⁷S) has no bearing on its spin once
bonded, which is set by oxidation state and ligand field — information the
counts model does not have. A bare or ionic transition metal resolves as
*underdetermined* in counts (to be supplied explicitly or by AtomTyping),
never forced low-spin (the old min-u bug) nor to a meaningless free-atom value.

### Why min-u / max-n was wrong

`min-u` is anti-Hund: it pairs electrons maximally, destroying the open-shell
ground state. It only ever produced correct atomic spin for groups 15–17, and
only because `covalence_set` coincidentally equalled the ground-state unpaired
count there (`N→3, O→2, F→1`). It fails wherever valence is reached by
promotion (groups 2/13/14: `min-u` and the `covalence_set` coincidence both
break) and across the empty-`covalence_set` d/f-block (collapses to wrong
low-spin). The earlier "effective valence `v+h+u ∈ covalence_set`" slack
formulation has the same defect — it forces `u = valence` for a bare atom
(`C → u=4`, not the ³P `u=2`).

### Verification

Worked through main-group species spanning s-block, groups 13–17, ions,
radicals, hypervalent, and implicit-H saturation. All resolve to the correct
ground state except two genuinely ambiguous cases the model handles honestly:

- `CH₂` / `[BH₂]⁻` (`v+h=2`, `K=2`, `R=2`): Hund gives **triplet** (`u=2`),
  the methylene ground state and the right default. The singlet (`u=0, n=1`)
  is a distinct electronic state electron-counting cannot rule out; it needs
  an explicit `#u0`. `[CH₃]⁻` (`K=1`) is *not* ambiguous — one orbital forces
  pairing.
- `NO₂`: a delocalized radical; the per-atom model places the unpaired
  electron per the input Kekulé form (on the singly-bonded O).

Representative results: `O atom → u2 n2` (³P), `H₂O → u0 n2`, `OH• → u1 n2`,
`CH₄ → u0`, `CH₃• → u1`, `[CH₃]⁻ → u0 n1`, `NH₃ → u0 n1`, `NH₂• → u1 n1`,
`PH₃ → u0 n1` (octet `K=1`, not hypervalent), `Mg → u0 n1`, `F atom → u1 n3`,
`F⁻ → u0 n4`, carboxylate `O⁻ → u0 n3` (via `O⁻→F` shift).

## Aromatic branch

Same scheme with `av ≠ 0` — no separate code branch. Conservation gains the
`av` term:

```
e − c = v + h + av + u + 2·n
```

### Aromaticity criterion

An atom is aromatic iff

```
is_in_aromatic_system()  ||  Aromatic(Undetermined).matches(aromatic_valence)
```

The first arm is idempotency — a prior pass inserted an aromatic system. The
second is the resolution input: SMILES lowercase `c` carries
`Aromatic(Undetermined)`, and SMILES has **no aromatic-system notion**, so on
the first pass there is no membership — the constraint match is the only
signal. `NotAromatic` is distinct from `Aromatic(Lit(0))`: an aromatic atom may
contribute `av=0` (acceptor, e.g. borazine B) and is still aromatic. So `av=0`
does **not** mean non-aromatic, and the criterion is never "membership" or
"`av>0`."

### Resolution

- For an aromatic atom, `av` is enumerated over `aromatic_valence_set`
  (charge-shifted, exactly like `covalence_set`). For a non-aromatic atom
  `av = 0`, `ai = 0`, no π orbital reserved — the non-aromatic formulas. Same
  arithmetic; only `av`'s domain differs.
- `aromatic_increment ai = (av == 1)`. `covalence_set` constrains **`v + h + ai`**
  — the covalent valence (σ bonds + implicit H + the `av=1` π half-bond). An
  `av=2` *donated lone pair* does **not** count (`ai=0`). This unifies with the
  non-aromatic branch (`ai=0` → `v+h`) and shares one `covalence_set`: benzene C
  and CH₄ both hit 4; pyridine N, pyrrole N, NH₃ all hit 3; furan O and H₂O
  both hit 2. No aromatic-specific valence numbers.
- `u` via Hund on residual `R = e − c − v − h − av`, over
  `K = 4 − (v+h) − [aromatic]` (the π orbital is reserved whenever the atom is
  aromatic by the criterion above — any `av`, including the `av=0` acceptor —
  not by membership alone, not by `av>0`). Aromatic atoms are ≈ always
  closed-shell (`u=0`).

Table records two per-element sets, both doubling as element-scope markers:

- `covalence_set` — admissible `v + h + ai` (shared σ/localized valences).
- `aromatic_valence_set` — admissible `av` (C→`[1]`, O→`[2]`, N→`[1,2]`,
  B→`[0]`); **non-empty = aromatic-capable** (the aromatic element scope).

`aromatic_valence_set` bounds `av` because `covalence_set` alone cannot: benzene
C admits both `(av=1, h=1)` and `(av=2, h=2)` under `[4]`; only
`aromatic_valence_set(C)=[1]` rules out the spurious `av=2`. The chain is
`aromatic_valence_set → av → ai → covalence_set check`.

### Aromatic verification

| species | atom | v | h | av | `v+h+ai` | `covalence_set` | u | n |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| benzene | C | 2 | 1 | 1 | 4 | [4] | 0 | 0 |
| pyridine | N | 2 | 0 | 1 | 3 | [3] | 0 | 1 |
| pyrrole | N | 2 | 1 | 2 | 3 | [3] | 0 | 0 |
| furan | O | 2 | 0 | 2 | 2 | [2] | 0 | 1 |
| naphthalene | peripheral CH | 2 | 1 | 1 | 4 | [4] | 0 | 0 |
| naphthalene | ring-fusion C | 3 | 0 | 1 | 4 | [4] | 0 | 0 |
| borazine | N | 2 | 1 | 2 | 3 | [3] | 0 | 0 |
| borazine | B | 2 | 1 | 0 | 3 | [3] | 0 | 0 |
| Cp⁻ | ring CH | 2 | 1 | 1 | 4 | [4] | 0 | 0 |
| Cp⁻ | `[cH-]` (C⁻≅N) | 2 | 1 | 2 | 3 | [3] shifted | 0 | 0 |
| tropylium | ring CH | 2 | 1 | 1 | 4 | [4] | 0 | 0 |
| tropylium | `[cH+]` (C⁺≅B) | 2 | 1 | 0 | 3 | [3] shifted | 0 | 0 |

### Charged aromatic ions resolve directly (no equalization)

The charge-shift applies to **both** `covalence_set` and `aromatic_valence_set`,
so a charged aromatic carbon is handled as its isoelectronic neighbor:

- Cp⁻ `[cH-]`: `C⁻ → N`. `aromatic_valence_set(N)=[1,2]` admits `av=2`;
  `covalence_set(N)=[3]` selects it (`v+h+ai = 3`). The carbanion behaves like
  pyrrole N — a lone-pair π donor — charge stays on the atom. 6 π e = 4·1 + 2.
- tropylium `[cH+]`: `C⁺ → B`. `aromatic_valence_set(B)=[0]` → `av=0`;
  `covalence_set(B)=[3]` (`v+h+ai = 3`). The carbocation behaves like borazine B —
  an empty-p acceptor. 6 π e = 6·1 + 0.

Charge equalization (monoelement → move charge to `system.charge`) is a
*charge-representation canonicalization* (doc 93), independent of and not
required for resolvability. The earlier "equalization needed" reading was an
error of looking up `covalence_set(C)=[4]` without the charge shift.

### `covalence` vs `total_valence` (distinct accessors)

`AtomView::covalence()` = `v + h + ai` (ordinary covalent bonds from the atom's
own electrons; decision 12). `AtomView::total_valence()` (atom.rs:340) =
`v + h + av + mv` (full electron-sharing sum; diverges from SMARTS `v` for
lone-pair donors). They differ on `ai` vs full `av` (the `av=2` donated pair)
and `mv` — so both exist as separate accessors. The counts `covalence_set`
check is `covalence() ∈ covalence_set(element.shift(2·accepted − charge))`;
charge and accepted dative perturb the lookup element, not `covalence()`
itself.

## Comparison harness (temporary)

`umol-graph/src/ops/valence/counts_new.rs` (`CountsNewResolver`, `ValenceScheme`)
and `umol-graph/tests/counts_comparison.rs` run the new path against a
reference over the resolution corpus, comparing chemically meaningful per-atom
fields (charge, implicit H, lone pairs, spin), ignoring constraint promotion.
Both are marked for deletion once the model lands.

Findings that drove decisions 7–10:

- Adding the charge **and** unpaired adjustment to the `covalence_set` bound
  fixed every carboxylate (the original `O⁻` divergence) — dropped new-vs-old
  errors from 153 to ~34/45.
- Comparing against **atom-typing** (the better reference for atomic spin)
  showed the residual divergences are bare d/f-block atoms where old counts'
  `min-u` gives wrong low-spin and atom-typing's registry gives correct
  high-spin — i.e. the comparison metric was penalizing the *fix*. This is
  what motivated scoping d/f-block `u` out (decision 8).

## Background (superseded analysis)

The initial restatement modelled disambiguation as min-u then max-n over the
invariant enumeration, with `covalence_set` / `aromatic_valence_set` as
narrowing ranges and a three-table shortcut (`covalence_set`,
`aromatic_valence_set`, `NormalValenceTable`). The pyridine/pyrrole traces
worked under that rule because aromatic π enumeration plus an input `#h` pin
masked the spin defect. The spin defect only surfaced on bare atoms and
open-shell main-group species, leading to the Hund-based model above. The
aromatic branch folds into the same scheme via `av`/`ai` — see "Aromatic
branch".

## Doc 96 sequencing impact

- Step 4 (`ValenceModel` API methods): `Counts::resolve` per the spin model
  above; `Counts::validate` calls `ValenceInvariants::check`.
- Steps 5–6 (collapse resolvers / validators): unchanged.
- Step 7 (`NormalValenceTable` removal): folds in — the charge-shifted
  `covalence_set` lookup replaces it.

## Critical files

- `umol-graph/src/ops/valence/table.rs` — `covalence_set` / `aromatic_valence_set`; charge-shift lookup; `compute_implicit_hydrogens` removal.
- `umol-graph/src/ops/valence/normal_valence.rs` — removal.
- `umol-graph/src/ops/model.rs` — `ValenceModel::Counts` spin model (h via covalence_set, u via Hund).
- `umol-shared/src/element.rs` — Hund needs `group_base8` / block to derive the open subshell; octet capacity (4) is constant for main group.
- `umol-graph/src/ops/valence/counts_new.rs`, `umol-graph/tests/counts_comparison.rs` — temporary comparison scaffolding, delete after.
- `umol-graph/config/default-valence-table.toml` — renamed keys.
