# Table IR → AST raise (SMILES / MOL bridge)

Status: **Implemented** (v1). Module: `umol-graph/src/table_ir/raise.rs` (was `lift.rs`).
Depends on valence resolvers in [99-counts-as-invariants-2026-05-23.md](99-counts-as-invariants-2026-05-23.md).

Related: [86-molecule-ast-api-2026-04-16.md](86-molecule-ast-api-2026-04-16.md) (coordinates undecided),
[94-dsl-ast-io-ergonomics-2026-05-07.md](94-dsl-ast-io-ergonomics-2026-05-07.md) (`TryIntoAst` for Table IR),
[92-umol-graph-engines-restructure-2026-04-25.md](92-umol-graph-engines-restructure-2026-04-25.md) (`Molecule` = AST + positions).

## Problem

`parse_smiles` / `parse_mol` build `TableMolecule`, then `try_into_ast` produces
`MoleculeAst` for `Resolver`. Remaining gaps:

- IO raise duplicates DSL `raise_atom` (`raise_ground_atom` in `raise.rs`; unify later).
- CTAB `vvv` → `#v` not mapped (deferred; post-pass removed).

Valence resolution (Counts) works on EDN and IO-raised AST; neutral ground fields
come from raise defaults (`#c0`, `#h*`, …), not from parse-time coercion.

## Terminology

| Term | Meaning |
| --- | --- |
| **Table IR** | `table_ir::Molecule` — faithful parse of the input string |
| **Raise** | Table IR → `MoleculeAst` via `TryIntoAst` + `GROUND_ATOM_DEFAULTS` |
| **`raise_atom`** | DSL-only (`umol-ast`); IO uses **`raise_ground_atom`** until unified |
| **`lift_constraints`** | Unrelated — moves inline constraints to molecule scope on `MoleculeAst` |

`SmilesIoConfig` / `CtfileIoConfig` are **parse** configuration only. Raise defaults
are **hardcoded in `raise.rs`** (not user-configurable).

## Design principles

1. **Table IR = input truth** — parsers record what the format encodes (`Option`
   where absent). No implicit-H `Normal` marker; no neutral `charge: Some(0)` on SMILES.
2. **Raise = ground interpretation** — fixed rules fill neutral fields (`#c0`,
   `#u0`); `#h*`, `#n*`, … stay undetermined for Counts.
3. **Aromatic membership** — explicit on atoms only when the input encodes it
   (SMILES `aromatic: Some(bool)`). MOL ground atoms stay `aromatic: None`; CTAB
   type-4 bonds carry aromaticity on the bond. **Counts** may treat an atom as
   aromatic when a neighbor bond has `BondConstraint::Aromatic` (see doc 99 §8).
4. **Coordinates stay out of AST** — `TableMolecule.positions` remains a sidecar.

## Implicit hydrogens

### MOL (`hhh` + `M  HYD`)

| CTAB | Table IR `implicit_hydrogens` | Raised AST |
| --- | --- | --- |
| `hhh = 0` or field absent | `None` | `#h*` (`Required`) |
| `hhh = k ≥ 1` | `Some(k − 1)` | `#h Lit(k − 1)` (`k = 1` → `#h0`) |
| `M  HYD n` | `Some(n)` (overrides line) | `#h Lit(n)` |

**`hhh = 0` is not `#h0`.** `ImplicitHydrogens::Normal` retired; `Option<u8>` only.

### SMILES

| Source | Table IR | Raised AST |
| --- | --- | --- |
| Organic / aromatic organic | `None` | `#h*` |
| Bracket `[CH4]` | `Some(4)` | `#h4` |
| Bracket `[C]` | `Some(0)` | `#h0` |

## Charge, spin, lone pairs

| Field | Table IR (parser) | Raise |
| --- | --- | --- |
| `charge` | `None` unless bracket `+`/`-` (incl. `[C+0]`), CTAB `ccc` / `M CHG` | undetermined → `#c0` |
| `unpaired_electrons` | `None` unless radical | `Zero` → `#u0` |
| `multiplicity` | not in SMILES/MOL Table IR | `Required` (undetermined) |
| `lone_pairs` | `None` unless set | **`Required`** — never `Zero` at raise |

SMILES: no `charge_opt.or(Some(0))` / no `on_atom_fast` `charge: Some(0)`; neutral
organic and bracket atoms without charge token use `charge: None` in Table IR.

Radical/charge exclusivity on MOL property application unchanged (`M RAD` clears
charge, etc.).

## Aromatic valence

SMILES/MOL record **membership** on atoms when the format says so; π counts come
from resolution, not Table IR.

### Table IR (parse)

| Source | `Atom::aromatic` | Bonds |
| --- | --- | --- |
| SMILES organic lowercase | `Some(true)` | usually plain single/double (not `1#a`) |
| SMILES organic uppercase | `Some(false)` | — |
| SMILES bracket | not set from bracket alone | — |
| MOL / CTAB ground | **`None`** (atom line does not encode membership) | `BondOrder::Aromatic` (type 4) when present |

**Do not** infer `aromatic: Some(bool)` on MOL atoms from incident aromatic bonds at
parse time.

### Raise (`raise.rs`)

| Table IR `aromatic` | Atom constraint after raise |
| --- | --- |
| `Some(true)` | `#a+` (`Aromatic(Undetermined)`) |
| `Some(false)` | `#a!` (`NotAromatic`) |
| `None` | (none) |

`raise_ground_atom` does not add aromatic constraints; those come only from
`aromatic: Some(bool)` on Table IR.

### Bonds (raise)

`BondOrder::Aromatic` → `order = Lit(1)` + `BondConstraint::Aromatic` (`1#a`).

### Resolver (Counts)

For `CountsValence::resolve` on a molecule, aromatic **context** for enumerating
`#a` includes (doc 99): `is_in_aromatic_system()`, **neighbor bond
`BondConstraint::Aromatic`**, or atom `Aromatic(_)`. So MOL benzene with only
`1#a` bonds and `aromatic: None` on atoms still resolves after raise. Unit test:
`counts.rs` `test_counts_valence_resolve_molecule_atom::benzene_ring`.

Hueckel / aromaticity perception runs after valence and reads **`Aromatic(Lit(n))`**
written by Counts, not undetermined atom flags alone.

## Valence (`vvv`) — deferred mapping

| | CTAB `vvv` | AST `#v` |
| --- | --- | --- |
| `0` | unspecified | omit |
| `1..14`, `15→0` | rare explicit pin | **do not map in v1** |

Table IR may keep `valence: Option<u8>` when pinned. No `vvv` → `#v` post-pass in raise.

## Coordinates

`MoleculeAst` has no position container. `TableMolecule.positions` and CXSMILES
`|(...)|` stay on Table IR; raise does not copy them.

## Raise pipeline (IO)

```text
copy TableAtom → AtomAst     // Lit / Undetermined; constraints empty
copy bonds                   // Aromatic → Lit(1) + BondConstraint::Aromatic
if aromatic == Some(true)  → #a+
if aromatic == Some(false) → #a!
raise_ground_atom            // fixed IO semantics (see below)
```

`raise_ground_atom` encodes IO semantics directly (no `AtomDefaults` / config table).
DSL `raise_atom` stays private; dedupe later if worthwhile.

### `raise_ground_atom` (landed)

| Field / slot | If still `Undetermined` after Table IR copy |
| --- | --- |
| `isotope_mass` | `Natural` |
| `charge` | `Lit(0)` |
| `implicit_hydrogens`, `lone_pairs`, `multiplicity` | unchanged |
| `spin.unpaired` | `Lit(0)` |
| constraints | drop undetermined; `#a+`/`#a!` only from `aromatic: Some(bool)` above |

## Implementation status

| Phase | Status |
| --- | --- |
| 1 — `implicit_hydrogens: Option<u8>`, CTAB/SMILES parsers | **Done** |
| 2 — MOL aromatic flags from bonds at parse | **Rejected** — bonds only in Table IR |
| 3 — `raise_ground_atom` (direct semantics), aromatic from `Some(bool)`, bond `1#a` | **Done** |
| 4 — `raise.rs` tests, MOL/SMILES → AST + Counts spot checks | **Done** |

## Deferred

| Item | Notes |
| --- | --- |
| `vvv` → `#v` | Semantic audit of MDL valence vs AST `#v` |
| Query MOL `hhh ≥ 2` | CTAB minimum-H semantics → pattern, not `Lit(n)` |
| `M  HYD` in corpus | Parser support exists; little test pressure |
| CXSMILES coordination / multicenter | `#m+`, coordination bonds |
| Coordinates on `Molecule` wrapper | doc 86 |
| `BondDefaults` raise pass | Only if bond-level ground defaults needed |
| ChemDoodle / misaligned CTAB lines | Separate normalization |
| Re-export `RaiseError` variants | When strict raise checks are added |
| Share `raise_spin` / constraint retain with DSL | Optional; IO path is ~15 lines, no `AtomDefaults` |
| Resolution EDN fixtures with bond-only `#a` | Conformance inputs still use `#a+` on atoms; optional alignment |

## Out of scope (this doc)

- Changing `Resolver` / Counts algorithms (except documenting IO interaction).
- EDN/DSL raise defaults (`AtomDefaults::default()` vs `zeroed()` remain DSL concerns).
- Full Table IR → DSL round-trip.
