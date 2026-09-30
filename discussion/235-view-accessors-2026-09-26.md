# 235 — View names, collection access, and algorithm placement

Status: Proposed
Date: 2026-09-26
Relates: [213](213-editor-overlay-storage-2026-08-27.md),
[nomenclature guide](../docs/development/nomenclature.md),
[data-type guide](../docs/development/data-types.md),
[Python API guide](../docs/development/python-api.md)

## Scope and status

Review graph-IR views and the underlying accessors they call for misleading id names,
unnecessary iterator interfaces, hidden allocations, and graph algorithms implemented
at the view layer. Include entity sets and corresponding reaction-span accessors where
they expose the same operations. Trace Rust and Python consumers before proposing changes.

This is a separate review, not a prerequisite for doc 213 S2j. S2i remains complete.
Molecule/editor ownership and mutation contracts are unchanged. No implementation is
authorized here; the examples below seed the review rather than constitute a complete audit.

## Required review output

For each finding, record the current signature and implementation, proposed signature
or delegation, reason, allocation/work differences, and affected callers. Distinguish
confirmed behavior from proposed changes. Keep straightforward renames separate from
return-type or algorithm changes so each group can be reviewed independently.

- Identify methods returning ids or positions without naming them accordingly. Stereo views use
  ligand_ids() for the ordered compound StereoLigand identifiers (atom id and kind), and
  atom_ligand_ids() for the actual-atom subset as AtomId values.
- Compare borrowed slices, fixed arrays, owned vectors, and iterators according to what
  the operation actually produces. Include iterator inputs that are immediately collected.
- Trace allocation through delegated calls. An iterator return type does not establish
  allocation-free or fully lazy execution. Identify materialization followed by iteration
  or recollection, without assuming every such sequence is unnecessary.
- Retain useful lazy conversion/filtering. Do not replace iterators with allocated vectors
  merely for uniformity, or introduce representation casts to obtain borrowed slices.
- Keep graph traversal and selection algorithms in graph-core. Graph IR supplies entity
  ids and molecular meaning; views expose access to those operations.
- Preserve ordering, multiplicity, borrowing, and failure contracts unless a specific
  change is proposed and approved. Account for Python ownership and materialization explicitly.

## Initial evidence and candidates

### Id names

The stereo-view rename ligand_frame() → ligand_ids() is settled and implemented separately
from this review. It retains the borrowed &[StereoLigand] return type and stored order.

[RingView](../umol-graph-ir/src/ir/view/ring.rs) exposes atoms() -> &[AtomId] and
bonds() -> &[BondId]. Proposed names are atom_ids() and bond_ids(), preserving return
types and behavior. Underlying entity sets also need review: DativeBonds::donors and
acceptor return atom ids; AromaticSystems::atoms, MulticenterBonds::atoms, and
NoncovalentBonds::atoms likewise return ids. Their span containers have corresponding
accessors. These are concrete inventory entries, not approval for a blanket rename.

### Collection shapes

| Current interface | Candidate | Reason and limit |
| --- | --- | --- |
| BondView::atoms() -> impl ExactSizeIterator<Item = AtomView> | Return [AtomView; 2] | Matches NoncovalentBondView::atoms(); both avoid heap allocation. Constructs both views immediately and permits direct indexing/destructuring. Review iterator-style callers. |
| AromaticSystemViews / MulticenterBondViews::of_id(atoms: impl IntoIterator<Item = AtomId>) | Accept &[AtomId], with matching of() signature | Both collect a Vec solely to pass a slice to coincident_id. Removes that intermediate collection for callers already holding a sequence; arbitrary iterator callers would collect explicitly. |
| DativeBondView::donor_ids() -> impl ExactSizeIterator<Item = AtomId> | Retain pending review | Storage contains NodeId; lazy conversion produces AtomId without allocation. A borrowed AtomId slice is not available directly. |

Sources: [bond views](../umol-graph-ir/src/ir/view/bond.rs),
[noncovalent views](../umol-graph-ir/src/ir/view/noncovalent.rs),
[aromatic views](../umol-graph-ir/src/ir/view/aromatic.rs),
[multicenter views](../umol-graph-ir/src/ir/view/multicenter.rs),
[dative storage](../umol-graph-ir/src/ir/dative.rs).

### Induced-edge membership allocation

[Graph::induced_edges](../umol-graph-core/src/graph.rs) constructs a HashSet<NodeId>
from the supplied slice on every call, then returns a lazy filter over all graph edges.
AromaticSystemView::bond_ids() and bonds() delegate to it. The membership construction
is eager even if the result iterator is only partially consumed.

Review both the membership allocation and the full-edge scan. Compare direct slice
membership, sorted membership where ordering is established, and adjacency-based selection
where useful; do not choose a replacement without checking costs and contracts. Preserve
the existing independence from input order/repetitions and retention of loops and parallel
edges. Do not require callers to construct a hash set to compensate for the implementation.
Also inspect induced selection in overlay views, which constructs its own membership sets.

### Noncovalent connectivity at the wrong layer

In [noncovalent views](../umol-graph-ir/src/ir/view/noncovalent.rs),
noncovalent_bond_derived_constraint computes intramolecularity through
same_bond_component(molecule, a, b). That function implements traversal itself, allocating
a visited vector and a work vector and walking molecule.neighbors(). Its comment says
breadth-first, but Vec::pop makes the work vector a stack.

The required correction is delegation to graph-core. The molecular meaning stays here:
the two endpoints belong to the same localized-bond component; overlay relations do not
connect components for this predicate. Graph-core already exposes visit_breadth_first and
visit_connected_components. Review whether an existing operation efficiently answers the
pair query or a core connectivity predicate is warranted. Avoid enumerating every component
or materializing full component membership just to obtain one boolean. The specific core
interface and its algorithm selection remain to be settled in this review.

## Next action

Complete the accessor/caller inventory and present bounded findings with exact signatures.
Settle names and algorithm placement before writing implementation subitems. No general
replacement of iterators, new cache, or new view hierarchy is implied by this review.
