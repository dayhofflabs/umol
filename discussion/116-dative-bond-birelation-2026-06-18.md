# 116 — Dative bond → fixed-var birelation

Status: Active · 2026-06-18

## Problem

Dative bonds are stored as a single variable-arity relation:

```rust
dative_bonds: Arc<VarRelationSet<NodeId, Unordered, DativeBondAst>>
```

Donors and the acceptor are flattened into one sorted participant list, and the
acceptor is identified by an index *into that list* carried inside the data:

```rust
pub struct DativeBondAst {
    pub acceptor_slot: u8,   // index of the acceptor in the sorted participants
    pub order: ValueAst,
    pub constraints: DativeBondConstraints,
}
```

`from_parts` / `add_dative_bond` concatenate `donors + acceptor`, `sort_unstable`,
then search for the acceptor's position to set `acceptor_slot`. The view reads
`atoms[acceptor_slot]` to recover the acceptor and filters it out to recover the
donors.

This is the index-into-storage pattern the lattice redesign is removing
elsewhere: a structural role (which participant is the acceptor) lives as a
numeric offset inside the value, meaningless on a standalone `DativeBondAst`,
and a `meet`/`join`/`matches` anchor that has to special-case `acceptor_slot`
equality. It also blocks a clean `Canonicalize` for `DativeBondAst` (deferred
"P4.6" in doc 113): the slot is a tier-1 structural field tangled into the data.

## Target

A two-factor birelation, with the acceptor as a fixed single-element factor and
the donors as the variable factor:

```rust
dative_bonds:
    Arc<FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondAst>>
```

This is the stereo overlay structure (`stereo_atoms`):

```rust
FixedVarBirelationSet<NodeId, Ordered, 1, StereoLigand, Ordered, StereoAtomAst>
```

with two differences:
- factor 2 is plain `NodeId`, not `StereoLigand` (no ligand wrapper);
- factor 2 is `Unordered`, not `Ordered` — donors are a set with no sequence or
  parity, so the factor canonicalizes by sorting (the same observable order the
  old `Unordered` `VarRelationSet` produced). Stereo uses `Ordered` only because
  ligand sequence encodes the permutation/coset.

`DativeBondAst` loses `acceptor_slot` and becomes a pure value type
`{ order, constraints }`; the acceptor role is structural in the birelation, not
data. `meet`/`join`/`matches` drop the slot anchor and reduce to the value
fields, exactly like `StereoAtomAst` (whose anchor is the birelation, not the
data).

## Scope boundary

`acceptor_slot` / `AcceptorSlot` occur only in umol-ast (`dative.rs`,
`molecule.rs`, `views/dative.rs`, `builder.rs`, `transact.rs`, `edit.rs`, their
tests, the property strategy, and one `dsl/dative.rs` test). umol-graph and
umol-io never see the slot: they build dative bonds through the high-level
`(donors: Vec<AtomId>, acceptor: AtomId, DativeBondAst)` input (the `from_parts`
tuple and `add_dative_bond`). Those signatures stay; only the internal storage
and the slot change. The DSL render/parse go through the view API
(`donor_ids()` / `acceptor_id()`), so they need no structural change either.

`FixedVarBirelationSet` (graph-core) and its mutable builder wrapper
`FixedVarSetStorage` already exist and are exercised by stereo, so no new
infrastructure is required.

## Combined accessors removed

`DativeBondView::atoms()` / `atom_ids()` returned the union of donors and
acceptor as one sorted list — a shape that only made sense while storage held a
single flat participant list. They have no current consumers (only their own
view tests), so they are dropped for now rather than reconstructed from the two
factors. Callers use `donor_ids()` / `donors()` and `acceptor_id()` /
`acceptor()`; `connecting_id` / `induced_ids` and the builder's undo records
compute the participant union internally from the two factors. If a combined
accessor is needed later, it returns over `participants_1 ∪ participants_2` with
an order chosen then.

## Plan

| Phase | File(s) | Change |
|---|---|---|
| 1 | umol-graph-core | none — `FixedVarBirelationSet` already present |
| 2 | `ast/dative.rs` | drop `acceptor_slot` field, `with_acceptor_slot`, and slot handling in `meet`/`join`/`matches`; struct → `{ order, constraints }`; retype tests; land `Canonicalize` (the deferred P4.6) as a plain value-type derive |
| 3 | `ast/molecule.rs` | field + `from_arcs` type → `FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondAst>`; `from_parts` maps `(donors, acceptor, d)` → `([acceptor.into()], donors→NodeId, d)` with no sort/slot; update clone/eq |
| 4 | `ast/views/dative.rs` | views hold the birelation; `acceptor_id` ← `participants_1(rid)[0]`, donors ← `participants_2(rid)`; **remove `atoms()` / `atom_ids()`** (no consumers); drop `acceptor_slot()` accessor; `connecting_id`/`induced_ids` union both factors internally; remove the `atom_ids`/`atoms` view tests, retype the rest |
| 5 | `ast/molecule/builder.rs` | storage → `FixedVarSetStorage<NodeId, Ordered, 1, NodeId, Unordered, DativeBondAst>`; `add_dative_bond(donors, acceptor, bond)` pushes `([acceptor], donors)` — no slot; builder views read the factor split; remap/restore via the FixedVar paths (mirror stereo) |
| 6 | `ast/edit.rs`, `ast/molecule/transact.rs` | remove `DativeBondFieldChange::AcceptorSlot` (variant, apply arm, inverse); keep `Order`; retype affected transact tests |
| 7 | `dsl/dative.rs` (tests), `dsl/molecule.rs` | no structural change (render/parse use the view API); drop the `acceptor_slot` struct field from the `dsl/dative.rs` test |
| 8 | `tests/property/strategies.rs` | drop `acceptor_slot` from the dative strategy; molecule-level strategy emits `(donors, acceptor)` |

## Consequence

`DativeBondAst` joins the other entity value types as a clean
`{ order, constraints }` lattice value with a derivable `Canonicalize`, closing
the P4.6 item from doc 113. The acceptor/donor distinction is expressed once,
structurally, in the birelation — no offset markers, no standalone-meaningless
fields, and no slot special-casing in the lattice ops.
