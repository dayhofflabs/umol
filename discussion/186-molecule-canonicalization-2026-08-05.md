# 186 - Molecule AST canonicalization

Status: Proposed
Date: 2026-08-05
Relates: [156](156-ast-comparison-and-property-suite-2026-07-20.md),
[113](113-ast-canonical-equality-and-lattice-2026-06-14.md),
[185](185-python-reaction-span-2026-08-04.md)

`MoleculeAst` does not implement Canonicalize and does not provide a canonical representation with respect to atom renumbering. canonical_key() is available from umol-graph-core for the nodes and edges (algorithms/automorphism). Incidence (Levi) graph construction is available umol-ast (ast/incidence.rs). Relevant considerations are:
- Need to compute a canonical frame and renumber accordingly. Consult 156-ast-comparison-and-property-suite-2026-07-20.md for the relation between equiv (relation under frame changes) and canonical equality.
- Canonicalizaton of disconnected molecules: does it require a splitting and recombining the connected components? Canonical ordering does not group connected components in general.
- Is the graph and overlays scheme sound for canonicalization? Consider meso compounds as examples. Any other complications? Very likely need to iterate to a fixpoint. StereoModel contains a parameter for parastereo.
- Aromatic and delocalized compounds: need to define canonical Kekule structure (in kekulizer); need to identify an aromatic system (in aromaticity perception) and define a canonical Kekule structure under the aromaticity model without running aromatizer. Should be formally equivalent to aromatization, followed by canonical kekulization.
- To consider, not required in this scope. Is cost of the subdivided graph construction for automorphism calculation a blocker? This complication is particular to the nauty implementation because the latter cannot use edge colors. If the cost becomes determining factor, consider adding an alternative implementation that can directly include edge colors.
`ReactionAst` and `ReactionSpanAst` also do not implement canonicalization. Deltas can be canonicalized, reaction AST canonicalization needs molecule canonicalization plus remapping of the deltas. Reaction span can generate its own node and edge color scheme based on entity spans.

## Required molecule-remapping operation

Canonicalization requires a general end-to-end `MoleculeAst` remapping operation. This is a public
AST transformation rather than a canonicalization-specific helper: canonicalization derives a
canonical labeling, represents it as a `MoleculeCorrespondence`, and applies the same operation that
other callers can use to transport a standalone molecule between dense id spaces. The public method
name remains to be approved before implementation.

The operation has the following contract:

- the correspondence source counts equal the molecule's counts for all eight entity families;
- every component correspondence is total on both sides and therefore defines a bijection onto a
  dense target id space;
- topology, relation participants, position-sensitive relation data, stereo frames, entity ASTs,
  and all references in constraints are transported together;
- it performs no chemistry validation, resolution, attribute canonicalization, repair, compaction,
  or entity removal; and
- failure of the correspondence to describe such a dense renumbering is ordinary absence and is
  reported with `Option`.

The implementation coordinates `umol_graph_core::Remapping`, which owns topology and relation
participant transport, with `IdRemapping`, which owns typed references across all eight entity
families. It must not reimplement relation participant sorting or payload permutation at the
`MoleculeAst` layer. The graph-core relation-remapping correction and its immediate reaction-span
consumer remain in doc 185, S3c; this work consumes that corrected facility.

Required properties are:

- semantic preservation, stated directly as
  `source.equiv_under(&remapped, &correspondence)`;
- exact identity remapping;
- inverse round-tripping;
- agreement between sequential remapping and correspondence composition; and
- preservation of referential integrity for every entity family and constraint reference.

The generated cases must exercise crossing permutations, all eight entity families,
position-sensitive relation data and stereo frames, and constraints that reference remapped
entities. These are properties of molecule remapping itself, independent of its use by
canonicalization. Canonicalization properties additionally verify that applying the remapping
derived from the canonical labeling produces the canonical representative.

This operation does not replace embedding into an ambient union namespace. A reaction-side mapping
may target a sparse or larger union id space and can transport entries into that space, but it cannot
produce a standalone dense `MoleculeAst` without a separate dense reindexing.
