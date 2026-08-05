# 186 - Molecule AST canonicalization

Status: Proposed
Date: 2026-08-05
Relates: [156](156-ast-comparison-and-property-suite-2026-07-20.md),
[113](113-ast-canonical-equality-and-lattice-2026-06-14.md)

`MoleculeAst` does not implement Canonicalize and does not provide a canonical representation with respect to atom renumbering. canonical_key() is available from umol-graph-core for the nodes and edges (algorithms/automorphism). Incidence (Levi) graph construction is available umol-ast (ast/incidence.rs). Relevant considerations are:
- Should use IdRemapping internally? The IdRemapping API is quite sparse, does it need a redesign to support the current use cases?
- Need to be able to generate renumbered molecule, and canonical equality.
- Need to compute a canonical frame and renumber accordingly. Consult 156-ast-comparison-and-property-suite-2026-07-20.md for the relation between equiv (relation under frame changes) and canonical equality.
- Canonicalizaton of disconnected molecules: does it require a splitting and recombining the connected components? Canonical ordering does not group connected components in general.
- Is the graph and overlays scheme sound for canonicalization? Consider meso compounds as examples. Any other complications? Very likely need to iterate to a fixpoint. StereoModel contains a parameter for parastereo.
- Aromatic and delocalized compounds: need to define canonical Kekule structure (in kekulizer); need to identify an aromatic system (in aromaticity perception) and define a canonical Kekule structure under the aromaticity model without running aromatizer. Should be formally equivalent to aromatization, followed by canonical kekulization.
- To consider, not required in this scope. Is cost of the subdivided graph construction for automorphism calculation a blocker? This complication is particular to the nauty implementation because the latter cannot use edge colors. If the cost becomes determining factor, consider adding an alternative implementation that can directly include edge colors.
`ReactionAst` and `ReactionSpanAst` also do not implement canonicalization. Deltas can be canonicalized, reaction AST canonicalization needs molecule canonicalization plus remapping of the deltas. Reaction span can generate its own node and edge color scheme based on entity spans.
