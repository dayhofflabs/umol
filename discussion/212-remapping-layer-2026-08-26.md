# 212 — The remapping id-transport layer

Status: Proposed
Date: 2026-08-26
Relates: [211](211-relation-frames-and-api-2026-08-26.md),
[214](214-aggregate-frame-semantics-2026-08-28.md),
[nomenclature guide](../docs/development/nomenclature.md)

## Purpose

Doc [211](211-relation-frames-and-api-2026-08-26.md) completes the compaction row of the
id-transport types because relation-set `compact` returns a relation compaction. That work exposed
the same missing layer in the remapping row. Remapping is not on 211's critical path and its
extraction turns on a representation question that compaction does not raise, so it is recorded
separately.

This document records the findings and the open decision. No implementation plan is appropriate
until the decision below is settled.

## Finding

There are three id-transport concepts. Each should exist at three layers: one id space, the graph's
node and edge spaces, and the molecule's eight entity families. Only correspondence has all three.

| concept | single id space | graph, node and edge | molecule, eight families |
| --- | --- | --- | --- |
| partial bijection | `Correspondence<Id>` `correspondence.rs:50` | `GraphCorrespondence` `correspondence.rs:296` | `MoleculeCorrespondence` `ir/correspondence.rs:25` |
| removal and dense shift | absent, added by doc 211 | `Compaction` `graph.rs:566` | `IdCompaction` `ir/remap.rs:23` |
| total relabel | absent | `Remapping` `graph.rs:643` | `IdRemapping` `ir/remap.rs:235` |

The third layer is also irregularly named. `MoleculeCorrespondence` names the id space it covers;
`IdCompaction` and `IdRemapping` do not.

## How remapping differs from compaction

Doc 211 extracts the compaction layer because the arithmetic is duplicated: `compact_relation`,
`uncompact_dense`, and `normalize_removed` in `ir/remap.rs` are line-for-line reimplementations of
`Compaction::compact_node`, `Compaction::uncompact_node`, and the sort-and-dedup in
`Compaction::new`. That extraction is mechanical and deletes code.

Remapping has no such duplication. The two layers use different representations for the same
concept.

```rust
pub struct Remapping {                  // graph.rs:643
    nodes: Vec<NodeId>,                 // dense, indexed by source id
    edges: Vec<EdgeId>,
}
// try_map_node(old) = self.nodes.get(old.0 as usize).copied()

pub struct IdRemapping {                // ir/remap.rs:235
    atom: HashMap<AtomId, AtomId>,      // sparse
    bond: HashMap<BondId, BondId>,
    dative: HashMap<DativeBondId, DativeBondId>,
    // five more
}
```

Extracting a shared `Remapping<Id>` therefore requires choosing a representation first. That is a
design decision, not an extraction.

## Evidence bearing on the representation

The two `to_remapping` operations are the same operation at two layers, and both are built from a
correspondence that is total on the left over a dense left id space.

`GraphCorrespondence::to_remapping` (`correspondence.rs:365`) checks `is_total_on_left`, then
collects only the right id of each matched pair into a positional vector. That is valid precisely
because the left space is dense and `matched_pairs` is in left order.

`MoleculeCorrespondence::to_remapping` (`ir/correspondence.rs:246`) performs the identical check and
then collects the whole pair into a hash map, discarding the density it just established.

So on this construction path the sparse representation is not required. Whether it is required at
all turns on the other construction sites, which have not been audited: `ir/edit.rs:1350` and
`:1385`, `ir/reaction_span.rs:1176`, `:1500`, `:1827`, and `:2507`, `ir/molecule.rs:1925`, `:2098`,
and `:2622`, `ir/molecule/pushout.rs:233`, and `ir/constraint/molecule.rs:826`.

The declared contract points the other way. `IdRemapping` is documented as mapping "every referenced
atom / bond / overlay id to its image in the target id space", used to move `Delta` values between
id spaces for `reverse` and `compose`, with the requirement that "every id a moved delta references
must be present". A delta may reference an arbitrary subset of a molecule's ids, which is consistent
with a sparse map but does not by itself exclude a dense one, since the molecule id spaces are
themselves dense.

Resolving this requires reading the eleven construction sites above and establishing whether any
supplies a source id set that is genuinely sparse or non-dense, rather than merely partial over a
dense space.

## Open decision

Select the representation for the extracted single-id-space layer:

1. **Dense positional**, as `Remapping` uses today. `IdRemapping` becomes eight positional vectors.
   Requires that every construction site supplies a dense source space.
2. **Sparse map**, as `IdRemapping` uses today. `Remapping` becomes two hash maps, which changes the
   graph-level lookup from an indexed read to a hash lookup on a hot path.
3. **Generic over the representation**, keeping the dense form at the graph layer and the sparse form
   at the molecule layer behind one interface. Retains both costs and adds a parameter.

The audit described above must precede the choice. Option 1 is only available if the audit supports
it, and option 2 has a measurable cost that should be established rather than assumed.

## Naming

If the layer is extracted, the row should be regularised alongside the other two:
`Remapping<Id>`, `GraphRemapping`, `MoleculeRemapping`. These are proposed and subject to the
nomenclature guide. Doc 211 makes the corresponding compaction change, so a decision here that
leaves the remapping row unchanged would leave the three concepts named inconsistently.

## Scope boundary

This document owns the remapping row only. The compaction row, the relation-set surface, frame
transport, and the `Reframe` operation belong to doc
[211](211-relation-frames-and-api-2026-08-26.md) and are not reopened here. Nothing in doc 211
depends on this document, and doc 211 may be implemented and closed while this remains proposed.

Doc 211 should nonetheless land its compaction row first, so the two aggregates are renamed against a
settled precedent rather than concurrently. The audit of the eleven construction sites above must
precede the representation decision, and the staged plan follows that.
