# 116 — Dative bond → fixed-var birelation

Status: Completed · 2026-06-18

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

## Combined accessors kept, unsorted

`DativeBondView::atoms()` / `atom_ids()` / `atom_count()` return the union of
donors and acceptor. They have real consumers — incidence-graph construction
(`incidence.rs`), host-subgraph filtering (`molecule.rs`), and reaction rewrite
(`molecule/rewrite.rs`) — every one of which uses the union as an *unordered
set* (membership tests, edge building). So they stay, and the ordering question
is moot: `atom_ids()` yields donors (factor 2, already sorted by the `Unordered`
canonicalization) followed by the acceptor (factor 1), with no extra sort.
`connecting_id` / `induced_ids` build the same union internally from the two
factors.

## Plan

| Phase | File(s) | Change |
|---|---|---|
| 1 | umol-graph-core | none — `FixedVarBirelationSet` already present |
| 2 | `ast/dative.rs` | drop `acceptor_slot` field, `with_acceptor_slot`, and slot handling in `meet`/`join`/`matches`; struct → `{ order, constraints }`; retype tests. **`Canonicalize` is *not* landed here** — `DativeBondConstraints` doesn't impl it yet, so no entity has `Canonicalize`; dropping `acceptor_slot` clears one of P4.6's two blockers, the constraint-collection blocker remains (P5) |
| 3 | `ast/molecule.rs` | field + `from_arcs` type → `FixedVarBirelationSet<NodeId, Ordered, 1, NodeId, Unordered, DativeBondAst>`; `from_parts` maps `(donors, acceptor, d)` → `([acceptor.into()], donors→NodeId, d)` with no sort/slot; update clone/eq |
| 4 | `ast/views/dative.rs` | views hold the birelation; `acceptor_id` ← `participants_1(rid)[0]`, donors ← `participants_2(rid)`; drop `acceptor_slot()` accessor; `atom_ids`/`atoms`/`atom_count` kept, yielding donors-then-acceptor unsorted; `connecting_id`/`induced_ids` union both factors internally; retype tests |
| 5 | `ast/molecule/builder.rs` | storage → `FixedVarSetStorage<NodeId, Ordered, 1, NodeId, Unordered, DativeBondAst>`; `add_dative_bond(donors, acceptor, bond)` pushes `([acceptor], donors)` — no slot; builder views read the factor split; build-remove uses `birelation_removed` + FixedVar default; restore: current entries via `restore_birelation_participants`, removed records via `atoms.split_last()` → `([acceptor], donors)`. `Added/RemovedDativeBond` keep a flat `atoms: Vec<AtomId>` (acceptor last, matching `atom_ids()` and the Edit layer's `split_last`) rather than splitting into `acceptor`/`donors` like the stereo records — keeps the change off the `edit.rs`/`transact.rs` record construction |
| 6 | `ast/edit.rs`, `ast/molecule/transact.rs` | remove `DativeBondFieldChange::AcceptorSlot` (variant, apply arm, inverse); keep `Order`; retype affected transact tests |
| 7 | `dsl/dative.rs` (tests), `dsl/molecule.rs` | no structural change (render/parse use the view API); drop the `acceptor_slot` struct field from the `dsl/dative.rs` test |
| 8 | `tests/property/strategies.rs` | drop `acceptor_slot` from the dative strategy; molecule-level strategy emits `(donors, acceptor)` |

## Consequence

`DativeBondAst` becomes a clean `{ order, constraints }` lattice value. P4.6
(`Canonicalize`) is no longer blocked by `acceptor_slot`, but still waits on P5
(constraint collections implementing `Canonicalize`) like every other entity —
it is not landed by this migration. The acceptor/donor distinction is expressed
once, structurally, in the birelation — no offset markers, no
standalone-meaningless fields, and no slot special-casing in the lattice ops.
