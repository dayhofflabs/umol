# Table IR → AST raise (SMILES / MOL bridge)

Status: **Active** — implementation plan. Renamed module: `umol-graph/src/table_ir/raise.rs`
(was `lift.rs`). Depends on valence resolvers in [99-counts-as-invariants-2026-05-23.md](99-counts-as-invariants-2026-05-23.md).

Related: [86-molecule-ast-api-2026-04-16.md](86-molecule-ast-api-2026-04-16.md) (coordinates undecided),
[94-dsl-ast-io-ergonomics-2026-05-07.md](94-dsl-ast-io-ergonomics-2026-05-07.md) (`TryIntoAst` for Table IR),
[92-umol-graph-engines-restructure-2026-04-25.md](92-umol-graph-engines-restructure-2026-04-25.md) (`Molecule` = AST + positions).

## Problem

`parse_smiles` / `parse_mol` build `TableMolecule`, then `try_into_ast` produces
`MoleculeAst` for `Resolver`. The bridge today:

- Maps fields ad hoc in `raise.rs` (no shared `raise_atom`).
- Uses `ImplicitHydrogens::Normal` in Table IR (redundant with neutral + undetermined `#h`).
- Sets atom `#a+` only when `aromatic == Some(true)`; never `#a!`.
- Maps CTAB `vvv` to `#v` without semantic review.
- Leaves MOL atoms with `aromatic: None` despite aromatic bonds.

Valence resolution (Counts) now works on EDN inputs; IO-raised AST must supply
compatible ground pins (`#c0`, `#h*`, `#a+` / `#a!`, …) without over-constraining
`#n` or spurious `#v`.

## Terminology

| Term | Meaning |
| --- | --- |
| **Table IR** | `table_ir::Molecule` — faithful parse of the input string |
| **Raise** | Table IR → `MoleculeAst` via `TryIntoAst` + defaults (not “lift”) |
| **`raise_atom`** | DSL-only (`umol-ast`); IO uses separate **`raise_ground_atom`** in `table_ir/raise.rs` until unified |
| **`lift_constraints`** | Unrelated — moves inline constraints to molecule scope on `MoleculeAst` |

`SmilesIoConfig` / `CtfileIoConfig` are **parse** configuration only. The
Table IR → AST contract is **`GROUND_ATOM_DEFAULTS` hardcoded in `raise.rs`**
(not user-configurable).

## Design principles

1. **Table IR = input truth** — parsers record what the format encodes (`Option`
   where absent). No `Normal` marker for implicit H.
2. **Raise = ground interpretation** — `AtomDefaults` in `raise.rs` turns
   neutral/absent fields into resolution-ready AST (e.g. `#c0`, `#u0`, `#h*`).
3. **Special cases at the boundary** — aromatic membership (`#a+` / `#a!`) and
   bond `1#a` are allowed as raise-time policy; SMILES/MOL do not carry π counts.
4. **Coordinates stay out of AST** — `TableMolecule.positions` remains a sidecar;
   raise does not copy them (see doc 86).

## Implicit hydrogens

### MOL (`hhh` + `M  HYD`)

Corpus (~113k V2000 atoms, scifinder/zinc/wwpdb): **~95%** infer (`hhh = 0` or
absent), **~2.5%** explicit (`hhh ≥ 1`), **0%** `M  HYD` in tree scanned.

| CTAB | Table IR `implicit_hydrogens` | Raised AST |
| --- | --- | --- |
| `hhh = 0` or field absent | `None` | `#h*` (`implicit_hydrogens: Required`) |
| `hhh = k ≥ 1` | `Some(k − 1)` | `#h Lit(k − 1)` (`k = 1` → `#h0`) |
| `M  HYD n` | `Some(n)` (overrides line) | `#h Lit(n)` |

No MOL-specific “infer vs zero” flag. **`hhh = 0` is not `#h0`.**

Retire `ImplicitHydrogens::Normal`; use `Option<u8>` only.

### SMILES

| Source | Table IR | Raised AST |
| --- | --- | --- |
| Organic | `None` (stop using `Normal` in builder) | `#h*` |
| Bracket `[CH4]` | `Some(4)` (`H` omitted → 0 in parser) | `#h4` |
| Bracket `[C]` | `Some(0)` | `#h0` |

## Charge, spin, lone pairs

| Field | Table IR (parser) | Raise (`GROUND_ATOM_DEFAULTS`) |
| --- | --- | --- |
| `charge` | `None` unless bracket `+`/`-`, CTAB `ccc` / `M CHG`; SMILES organic still `Some(0)` in v1 | `NumericDefault::Zero` → `#c0` when `Undetermined` |
| `unpaired_electrons` | `None` unless radical | `UnpairedElectronsDefault::Zero` |
| `multiplicity` | never in SMILES/MOL Table IR | `MultiplicityDefault::Required` (not `Derived` — nothing to derive from) |
| `lone_pairs` | `None` unless set | **`NumericDefault::Required`** — never `Zero` (would break Counts) |

**v1:** leave SMILES `charge_opt.or(Some(0))` and `on_atom_fast` `charge: Some(0)`.
Raised neutral SMILES matches MOL `None` → raise. Cleanup deferred.

Radical/charge exclusivity on MOL property application unchanged (`M RAD` clears
charge, etc.).

## Aromatic valence

SMILES/MOL record **membership**, not aromatic π count.

| Table IR `aromatic` | Pre-seed before `raise_atom` | Notes |
| --- | --- | --- |
| `Some(true)` | `#a+` (`Aromatic(Undetermined)`) | SMILES lowercase / flag |
| `Some(false)` | `#a!` (`NotAromatic`) | all aliphatic ground MOL atoms |
| `None` | (none) | queries only |

**MOL v1:** after parse, set `aromatic: Some(false)` on every atom, then
`Some(true)` on atoms incident to a bond with `BondOrder::Aromatic` (type 4).

**Bonds:** `BondOrder::Aromatic` → `order = Lit(1)` + `BondConstraint::Aromatic` (`1#a`).

Raise `AtomDefaults` uses `aromatic_valence: Required` so `raise_atom_constraints`
does not blanket-insert `#a!` over pre-seeded `#a+`.

## Valence (`vvv`) — deferred mapping

| | CTAB `vvv` | AST `#v` |
| --- | --- | --- |
| Meaning | MDL connection-table valence pin | Localized σ valence (sum of localized bond orders; see DSL spec) |
| `0` | unspecified | omit |
| `1..14`, `15→0` | rare explicit pin | **do not map to `#v` in v1** |

Table IR may keep `valence: Option<u8>` when pinned. Remove current
`raise.rs` post-pass that adds `Valence(Lit(v))` until MDL ↔ AST equivalence is
documented.

## Coordinates

`MoleculeAst` has no position container (doc 86: likely wrapper/sidecar).
`TableMolecule.positions` and CXSMILES `|(...)|` stay on Table IR; **raise does
not copy them.**

## Raise step (IO)

After mechanical copy `TableAtom` → `AtomAst`:

```text
copy_table_atom(ast)              // Lit / Undetermined from Table IR
preseed_aromatic_constraints()    // aromatic: Some(bool) → #a+ / #a!
raise_ground_atom(ast)            // IO-specific; see below
```

**Do not export** `raise_atom` from `umol-ast` (`dsl/atom.rs` stays private). v1
implements a **separate** `raise_ground_atom` (and bond analogue if needed) in
`umol-graph/src/table_ir/raise.rs`, mirroring the same `AtomDefaults` rules as
DSL `raise_atom` / `raise_atom_constraints` / `raise_spin`.

**After tests pass:** refactor to share one implementation (move to `umol-ast`
or a small shared module) without changing raised AST behavior.

Reference behavior (DSL, not called from IO):

- Fill `Undetermined` struct fields per `AtomDefaults`.
- `raise_spin` for `unpaired` / `multiplicity`.
- Constraint slots only where cfg demands and not already pre-seeded.

### Proposed `GROUND_ATOM_DEFAULTS`

| Slot | Mode |
| --- | --- |
| `isotope` | `Natural` |
| `charge` | `Zero` |
| `implicit_hydrogens` | `Required` |
| `lone_pairs` | `Required` |
| `unpaired_electrons` | `Zero` |
| `multiplicity` | `Required` |
| `valence` | `Required` |
| `donated_pairs` | `Required` |
| `accepted_pairs` | `Required` |
| `aromatic_valence` | `Required` (constraints pre-seeded) |
| `multicenter_valence` | `Required` |

## Implementation plan

### Phase 1 — Table IR model

- Replace `ImplicitHydrogens` with `implicit_hydrogens: Option<u8>` on `Atom`.
- CTAB: `hhh = 0` or absent → `None`; `hhh = k ≥ 1` → `Some(k − 1)`;
  `convert_atom_hydrogen_count_code` drops `Normal` arm.
- SMILES builder: organic/aromatic paths use `None` instead of `Normal`;
  bracket keeps `Some(n)` for `H` count.
- Update `table_ir` / CTAB / SMILES tests and `hydrogens_to_count` helpers.

### Phase 2 — MOL aromatic flags

- Post-parse (or during molecule build): every atom `aromatic: Some(false)`,
  then `Some(true)` if any incident bond is `BondOrder::Aromatic`.
- Tests: benzene/allene patterns in mol fixtures.

### Phase 3 — Raise integration

- Module-local `GROUND_ATOM_DEFAULTS` (may use `AtomDefaults` type from `umol_ast`).
- Per-atom: mechanical field copy from Table IR.
- Pre-seed `#a+` / `#a!` from `aromatic` bool.
- `raise_ground_atom` in `raise.rs` (duplicate of DSL raise logic, not a call into DSL).
- Remove aromatic-only post-loop and `vvv` → `#v` constraint injection.
- Bond/dative/multicenter copy unchanged from current `raise.rs`.

### Phase 4 — Verification

- Unit tests in `table_ir/raise/tests`.
- Spot-check: methane/benzene MOL → `#h*`, `#c0`, `#a+`/`#a!`; `hhh=1` → `#h0`.
- `parse_smiles_to_ast` / `parse_mol_to_ast` + Counts or resolution conformance
  where applicable.

## Deferred

| Item | Notes |
| --- | --- |
| SMILES charge coercion cleanup | Remove `charge_opt.or(Some(0))` / fast-path `Some(0)` so Table IR uses `None` for neutral; raised AST should be unchanged |
| `vvv` → `#v` | Semantic audit of MDL valence vs AST `#v`; rare in corpus |
| Query MOL `hhh ≥ 2` | CTAB minimum-H semantics → pattern, not `Lit(n)` |
| `M  HYD` in corpus | Parser support exists; no test pressure yet |
| CXSMILES coordination / multicenter | `#m+`, coordination bonds |
| Coordinates on `Molecule` wrapper | Thread through DSL / chemist API (doc 86) |
| `BondDefaults` raise pass | Only if bond-level ground defaults needed beyond copy |
| ChemDoodle / misaligned CTAB lines | Invalid or separate normalization pass |
| Re-export `RaiseError` variants | When strict raise checks are added |
| Unify IO raise with DSL `raise_atom` | After tests pass; deduplicate `raise_ground_atom` vs `dsl/atom.rs` without exporting private DSL API |

## Out of scope (this doc)

- Changing `Resolver` / Counts algorithms.
- EDN/DSL raise defaults (`AtomDefaults::default()` vs `zeroed()` remain DSL concerns).
- Full Table IR → DSL round-trip.
