# Counts valence resolver

Status: **Implemented** in `umol-graph/src/ops/valence/counts.rs`. `CountsValence`
derives per-atom `#h`, `#a`, `#n`, `#u`, `#s` from a `ValenceTable`, topology
`#v`, charge, and accepted dative pairs. `ValenceInvariants` remains a
separate check/enumeration module; Counts does not call it.

## Decisions (landed)

| # | Decision |
| --- | --- |
| 1 | Counts enumerates admissible `(h, a)` assignments, then fills `#n`/`#u` from the electron budget; one candidate is chosen by a fixed preference order. |
| 2 | `NormalValenceTable` removed; saturation comes from `target_covalences` on a charge-shifted element. |
| 3 | `ValenceTable` fields: `target_covalences`, `aromatic_valences` (TOML + `ValenceEntry`). |
| 4 | Counts and AtomTyping coexist. AtomTyping covers registry Lewis types Counts cannot pin; Counts covers electron bookkeeping from the table. |
| 5 | Shared tie-break: `compare::compare_valence_preference` — max `#h`, then max `#n`, then min `#u`. Used by Counts and AtomTyping multi-match paths. |
| 6 | Charge and accepted dative pairs enter table lookup as an **isoelectronic shift**: `element.shift(2·accepted − charge)`. |
| 7 | **Out-of-table** (no row after shift): electron bookkeeping only — no `target_covalences`. `#h*` / `#a*` → `0`; `#h` / `#a` literals honored; `#a+` (`Aromatic(Undetermined)`) → error. In-table `#a+` still enumerates from `aromatic_valences`. |
| 8 | Aromatic context: `is_in_aromatic_system()` **or** aromatic bond on a neighbor **or** `aromatic_valence()` constraint is `Aromatic(_)`. |
| 9 | `aromatic_increment` = `1` iff resolved `av == 1`, else `0`; headroom check `h + ai ≤ bonding_budget` when `#h` is free. |

## Conservation and inputs

Per atom (main-group table path):

```
nonbonding = valence_electrons(element) − charge − v − av − h
```

- `v` — topology valence (`AtomView::valence()`, localized σ).
- `h` — implicit hydrogens (field).
- `av` — aromatic valence count (constraint/field; enumerated when aromatic).
- Accepted dative pairs affect **which table row** is read, not this sum directly.

`#d`, `#m`, and multicenter contributions are outside Counts today (implicit
`#d0 #t0` scope in the invariant checker only).

## Algorithm (`derive_fields`)

1. **Table row** — `entry = table.entry(element.shift(2·accepted − charge))`.
2. **Out-of-table guard** — if `entry` is `None` and constraint is `#a+`, return
   `UndeterminedAromaticValence`.
3. **Bonding budget** — when `entry` is `Some` and `#h` is not a literal:
   `budget = first target_covalence ≥ v` minus `v` (targets sorted ascending).
   Literal `#h` skips the budget (no saturation from the table).
4. **Candidate `h`** — literal → `[h]`; in-table free `#h` → `0..=budget`;
   out-of-table free `#h*` → `[0]`; other non-literal shapes → `NoMatch`.
5. **Candidate `av`** — literal `#a` → `[a]`; in-table aromatic → values from
   `aromatic_valences`; otherwise → `[0]`.
6. **Filter** — drop pairs with `h + ai > budget` (`ai` from `aromatic_increment`);
   drop negative `nonbonding`; satisfy pinned `#n` / `#u` when present.
7. **Free `#n` / `#u`** — when both free: `u = nonbonding % 2`,
   `n = (nonbonding − u) / 2`, capped by `valence_capacity/2` for `n`.
   Pinned `#n` or `#u` infers the other from `nonbonding`.
8. **Select** — `max_by(compare_valence_preference)` over survivors; `meet` onto
   the input atom.

Multiplicity: `2u + 1` when `#s` is free and compatible with pinned `#u`.

## Aromatic behavior

- `#a+` on an in-table atom: `av` runs over `aromatic_valences` for the
  shifted element (e.g. C → `[1]`, N → `[1,2]`, O → `[2]`, B → `[0]`).
- `#a0`, `#a1`, `#a2` literals constrain enumeration via `matches_value`.
- Preference (max `#h`, max `#n`, min `#u`) disambiguates when several
  `(h, av)` pairs satisfy the table — e.g. benzene C with free `#h*` picks
  `h = 1, av = 1` over `h = 0, av = 2` when both fit the budget.
- Heteroatom discrimination in conformance often needs explicit `#h0` or `#h`
  on the input (pyridine `N#h0`, pyrrole `N#h`) so AtomTyping registry
  matching is not ambiguous; Counts alone does not replace that pin.

## Out-of-table elements

Elements absent from `default-valence-table.toml` (most transition metals,
lanthanides, etc.) after the isoelectronic shift:

| Input | Behavior |
| --- | --- |
| `#h` literal | used as-is |
| `#h*` | `h = 0` |
| `#a` literal | used as-is |
| `#a*` | `av = 0` |
| `#a+` | **error** (`UndeterminedAromaticValence`) |
| free `#n` / `#u` | same residual split as in-table |

Example: neutral Fe, `v = h = av = 0` → `nonbonding = 8` → `n = 4, u = 0`
(`Fe#n4`), without inventing a Lewis `target_covalence`.

No molecule-wide upfront rejection for out-of-table symbols; each atom is
resolved on its own row lookup.

## ValenceInvariants (separate)

`invariants.rs` implements the orbital–electron equation for **checking** and
for `enumerate_atom` (full `(h, av, n, u, …)` enumeration on a grounded
molecule). Counts reimplements a narrower derivation loop tuned to the table
and does not delegate to `enumerate_atom`.

## AtomTyping

Registry patterns supply concrete Lewis types (including d/f and exotic
modes). Counts does not replace that when the table has no row or when
registry tie-break semantics differ. Resolution conformance runs both
pipelines separately on the same EDN inputs.

## Conformance inputs

Resolution EDN uses `MoleculeDefaults` (`#h` defaults to required/undetermined
until resolved). Explicit `#h0` / `#h` pins on heteroatoms are load-bearing
for AtomTyping; bare element symbols rely on Counts/table saturation.

## Critical files

- `umol-graph/src/ops/valence/counts.rs` — `CountsValence::derive_fields`.
- `umol-graph/src/ops/valence/table.rs` — `ValenceTable`, `ValenceEntry`.
- `umol-graph/src/ops/valence/compare.rs` — `compare_valence_preference`.
- `umol-graph/src/ops/valence/atom_typing.rs` — registry resolver (shared preference).
- `umol-graph/src/ops/valence/invariants.rs` — equation check / `enumerate_atom`.
- `umol-graph/config/default-valence-table.toml` — per-element targets.
- `umol-graph/src/ops/resolver/valence.rs` — `ValenceModel::Counts` dispatch.
