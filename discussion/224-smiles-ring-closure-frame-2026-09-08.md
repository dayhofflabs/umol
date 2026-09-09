# 224 — SMILES ring-closure stereo frame

Status: Proposed
Date: 2026-09-08
Relates: [104](104-stereochemistry-implementation-plan-2026-05-31.md),
[153](153-format-parsing-outstanding-tasks-2026-07-18.md),
[217](217-rhea-participant-failures-2026-08-30.md)

## Purpose

The SMILES tetrahedral reading is inverted at every stereocenter that carries a ring-closing
digit. This document records the finding, the rule, why the bond list cannot carry the frame, the
settled design, and the staged plan.

## Finding

OpenSMILES §3.8.2 defines the frame of a tetrahedral center as the order of its bonds in the
string: the bond to the preceding atom first; an implicit hydrogen next, or first when the center
opens the SMILES; then, left to right after the bracket, each ring-closure digit at its position,
each branch, and the continuing chain. `@` lists the neighbors after the first anticlockwise when
viewed from the first toward the center, `@@` clockwise. "One exception to the atom order is when
these atoms are bonded to the chiral center via a ring bond. In these cases, it is the order of
the bonds to these atoms that should be considered."

umol takes the order from the TableIR bond list (`first_neighbor_toward_ordering` through
`AtomNeighbors`). The SMILES builder reserves a ring-closure bond's slot at its opening digit and
fills it at the closing (doc 104, 2026-06-07). At the opening atom the list order is the written
order. At the closing atom the ring bond precedes the bond to the preceding atom, one
transposition away from the rule, so the coset is inverted there.

Minimal pair, same isomer by the rule: `C[C@H]1CCCCO1` has the frame (CH3, H, O, CH2) with `@`;
`O1CCCC[C@@H]1C` has (CH2, H, O, CH3) with `@@`, and swapping CH3 and CH2 turns `@@` into `@`.
umol reads the second in the frame (O, H, CH2, CH3), yielding the mirror image. Asked directly,
umol reports the pair as canonically unequal and reports `C[C@H]1CCCCO1` equal to its mirror
image `O1CCCC[C@H]1C`. Controls whose digits open rings at the center, such as
`C[C@]1([H])CCCCO1` against `[H][C@@]1(C)CCCCO1`, read correctly.

The census surfaced it through CHEBI:40, (+)-pinoresinol, whose Rhea SMILES closes rings at
two of its four centers; the MOL wedge reading of the same record agrees with an independent 3D
reading at all four centers (doc 217, "Paired MOL and SMILES input for CHEBI:40"; scripts
`tetrahedral_arbitration.py` and the `wedge-pair-probe` bin under `scratch/rhea-census/census/`).
Three resolution fixtures were lowered from affected SMILES and encode the wrong isomers:
`stereo_tetrahedral/cis-1-2-dichlorocyclohexane.edn` (from `Cl[C@H]1CCCC[C@H]1Cl`) holds the
chiral trans isomer, `trans-1-2-dichlorocyclohexane.edn` the meso cis isomer, and
`alpha-d-glucopyranose.edn` the C4 epimer, alpha-D-galactopyranose. Any SMILES with a
stereocenter at a ring-closing digit is affected, including the usual writing of pyranoses.

## No bond order carries the frame

For `X A1(B) ... Y Z1 D` with the ring bond between the opening atom A and the closing atom Z,
the rule needs the ring bond before A–B at A and after Y–Z at Z, while A–B is written before
Y–Z. The three inequalities form a cycle, so no single position in one bond list satisfies both
ends whenever the opening atom has any bond written after its digit. Open order serves the
opening atom, close order the closing atom. The frame therefore has to be recorded by the parser,
the only place that sees it; the bond list order is free.

## Design

TableIR is a boundary record and reflects the input; it does not acquire graph semantics.

### TableIR ring closures

```rust
/// A SMILES ring-closure bond: its number as written and the position of the bond among the
/// bonds of each end atom, in the order the atom's bonds are written.
pub struct RingClosure {
    pub bond: u32,
    pub number: u32,
    pub opening: RingClosureEnd,
    pub closing: RingClosureEnd,
}

pub struct RingClosureEnd {
    pub atom: u32,
    /// Rank of the ring bond among this atom's bonds in written order; an implicit hydrogen is
    /// not a bond and does not count.
    pub position: u32,
}
```

`Molecule` and `ExtendedMolecule` gain `pub ring_closures: Vec<RingClosure>` after `bonds`;
`empty()` yields an empty vector, the record conversions copy it, and CTfile input leaves it
empty. `Bond::ring` and `ExtendedBond::ring`, never set and never read, are removed. For
`X A1(B)C ... Y Z1 D` the opening end is (A, 1) and the closing end is (Z, 1). Recording both
ends makes the record independent of the bond list order.

### SMILES builder

Both editors append the ring bond when the closing digit is seen, so the bond list is in close
order and no slot is reserved. Each editor keeps a per-atom count of bonds written so far,
incremented for the bond to the preceding atom, for each ring digit at the atom whether it opens
or closes, and for each branch or chain bond; a digit's `position` is the count before its own
increment. Ring closures are recorded for every parse; the `store_rings` flag and the
`(close_rank, open_index)` tuples go. Reaction SMILES record them per component molecule.
Direction and donation flips at the closing end, the mismatch errors, and the reuse of digits
are unchanged.

### CXSMILES

CXSMILES bond indices count ring-closure bonds at their closing digit. With the list in close
order they equal umol's bond indices; `BondIndexMap` and `remap_cx_bond_indices` are removed and
the CX entries apply their indices directly.

### Raise

`first_neighbor_toward_ordering(neighbors, ring_closures, atom_idx)` takes the atom's neighbors
in bond-list order and moves each ring-closure end at this atom to its `position`, then inserts
the virtual hydrogen at slot 0 when the atom opens the SMILES and at slot 1 otherwise, as today.
The `ChiralityFrame::FirstNeighborToward` documentation states that the frame is the written
order, with ring-closure bonds placed by `ring_closures`.

### Rejected

- Transporting the coset into the raise's frame in the parser: TableIR would carry umol's frame
  instead of the input's.
- A per-bond record (`Bond::ring: Option<RingClosure>`): bonds outnumber ring closures, and the
  record would be pair-relative.
- Any single bond list order: impossible, above.

### Consequences

Bond indices of parsed SMILES move from open order to close order; the parser `test_ring` rows in
`smiles/parser/tests/basic.rs` state the new order, the SMILES conformance snapshots record only
counts and are unchanged, and no stability guarantee exists before 1.0. The three fixtures are
regenerated from their SMILES with `verify_stereo write` and their snapshots recorded. Doc 104's
note that the bond list equals the written order is superseded at the closing atom by this
document; doc 104 gains a link. Doc 153 T4 (`ChiralityFrame`) and T9 (arrangement numbering)
remain open and unaffected in scope. Python exposes no TableIR, so no binding changes.

## Implementation plan

Each stage ends green; breaking subitems carry their caller migration.

S0 — TableIR ring closures

- S0a `table_ir::molecule` (or a `ring_closure` submodule): `RingClosure`, `RingClosureEnd`,
  the `ring_closures` field on both records, `empty()`, and the conversions. Tests: conversion
  rows carrying a ring closure; `empty()` rows. Breaking for struct-literal construction
  [dep: none].
- S0b `table_ir::bond`: remove `ring` from `Bond` and `ExtendedBond` and from the conversions.
  Tests: the conversion tables drop the field. Breaking [dep: none].

Gate: `cargo test -p umol-io --features conformance,proptest`.

S1 — Builder in close order with positions

- S1a `smiles::parser::builder`: both editors append ring bonds at the closing digit, keep the
  per-atom written-bond counts, and emit `ring_closures`; `store_rings`, `ring_bonds`,
  `closed_bonds`, `take_ring_bonds` are removed; `parse_smiles_inner` and the reaction parser
  stop threading ring-bond tuples. Tests: ring closures for a digit that opens after the
  preceding atom, after a branch, and at the first atom; two digits at one atom; `%nn`; a digit
  with a bond symbol; a reused digit; positions at both ends; reaction components. The
  `test_ring` and extended `test_ring_invalid_topology` rows state close order. Breaking
  [dep: S0a].
- S1b `smiles::parser::cx`: remove `BondIndexMap` and `remap_cx_bond_indices`; CX bond-indexed
  entries apply directly. Tests: the CX bond-index rows for ring-containing inputs keep their
  expected bonds under the new numbering; the map's unit tests are deleted. Migration
  [dep: S1a].

Gate: `cargo test -p umol-io --features conformance,proptest`.

S2 — Raise frame

- S2a `table_ir::raise::utils`, `table_ir::raise`: `first_neighbor_toward_ordering` takes
  `ring_closures` and places ring-closure ends at their positions; `raise_tetrahedral_stereo`
  passes `mol.ring_closures`; the `ChiralityFrame::FirstNeighborToward` doc comment. Tests:
  `C[C@H]1CCCCO1` and `O1CCCC[C@@H]1C` raise to cosets that are equal after transport to a common
  frame, `O1CCCC[C@H]1C` to the opposite; the explicit-hydrogen forms; a center opening two
  rings; a center closing one ring and opening another; the doc 104 pair unchanged. Breaking
  [dep: S1a].
- S2b `umol-graph` ingest tests: `Molecule` from the minimal pair is `canonical_eq`, from the
  mirror pair is not; the same for `Cl[C@H]1CCCC[C@H]1Cl` against its mirror image (meso) and
  `Cl[C@H]1CCCC[C@@H]1Cl` against its mirror image (chiral). Additive [dep: S2a].

Gate: `cargo test -p umol-io --features conformance,proptest` and
`cargo test -p umol-graph --features conformance`.

S3 — Fixtures and records

- S3a `umol-graph/tests/resolution`: regenerate `cis-1-2-dichlorocyclohexane.edn`,
  `trans-1-2-dichlorocyclohexane.edn`, `alpha-d-glucopyranose.edn` from their SMILES and record
  their snapshots; verify the cis fixture is equal to its mirror image and the glucose fixture is
  the alpha-D-gluco configuration. [dep: S2a].
- S3b Docs: doc 104 link, doc 217 open item 7 pointing here, this document's status. [dep: S3a].

Critical path: S0a → S1a → S1b, S2a → S2b, S3a → S3b.
