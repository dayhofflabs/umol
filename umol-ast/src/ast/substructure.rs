//! Substructure matching: occurrences of a pattern `MoleculeAst` within a host
//! `MoleculeAst`, as [`MoleculeEmbedding`]s. The receiver is the pattern, so it
//! parallels `pattern.matches(target)`: `pattern.substructure_matches(host, ...)`.
//!
//! Two strategies compose over the chosen subgraph-isomorphism algorithm:
//! `GraphAndOverlays` matches the localized atom-bond skeleton then post-verifies
//! overlays; `Incidence` matches the incidence (Levi) graph for hyperedge-only
//! connectivity.

use umol_graph_core::SubgraphIsomorphismAlgorithm;

use super::atom::AtomAst;
use super::bond::BondAst;
use super::embedding::MoleculeEmbedding;
use super::molecule::MoleculeAst;
use super::traits::Lattice;

/// Strategy for `substructure_matches`, each composing over a
/// [`SubgraphIsomorphismAlgorithm`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SubstructureMatchAlgorithm {
    /// Match the localized atom-bond skeleton, then post-verify the N-ary / special
    /// overlays against the atom correspondence.
    GraphAndOverlays,
    /// Match the incidence (Levi) graph — true hypergraph matching for connectivity
    /// carried solely by overlays.
    Incidence,
}

impl MoleculeAst {
    /// Occurrences of `self` (the pattern) within `host`, one [`MoleculeEmbedding`]
    /// per occurrence. Pattern predicates are evaluated as `pattern.matches(host)`
    /// against the host atom/bond augmented with its derived topological
    /// constraints.
    pub fn substructure_matches<'h>(
        &self,
        host: &'h MoleculeAst,
        strategy: SubstructureMatchAlgorithm,
        subiso: SubgraphIsomorphismAlgorithm,
    ) -> Vec<MoleculeEmbedding<'h>> {
        match strategy {
            SubstructureMatchAlgorithm::GraphAndOverlays => {
                self.substructure_matches_graph_and_overlays(host, subiso)
            }
            SubstructureMatchAlgorithm::Incidence => {
                unimplemented!("Incidence substructure strategy")
            }
        }
    }

    fn substructure_matches_graph_and_overlays<'h>(
        &self,
        host: &'h MoleculeAst,
        subiso: SubgraphIsomorphismAlgorithm,
    ) -> Vec<MoleculeEmbedding<'h>> {
        let pattern = self;
        if pattern.atoms().count() > host.atoms().count() {
            return Vec::new();
        }
        // Host match-targets: stored constraints with derived topological ones folded
        // in (last-wins). Precomputed once, not per candidate.
        let host_atoms: Vec<AtomAst> = host
            .atoms()
            .iter()
            .map(|a| a.ast.clone().with_constraints(a.derive_constraints()))
            .collect();
        let host_bonds: Vec<BondAst> = host
            .bonds()
            .iter()
            .map(|b| b.ast.clone().with_constraints(b.derive_constraints()))
            .collect();

        host.graph()
            .subgraph_isomorphisms(
                &pattern.graph(),
                &mut |query_atom, host_atom| {
                    pattern.atom(query_atom).ast.matches(&host_atoms[host_atom.index()])
                },
                &mut |query_bond, host_bond| {
                    pattern.bond(query_bond).ast.matches(&host_bonds[host_bond.index()])
                },
                subiso,
            )
            .into_iter()
            .map(|atom_map| MoleculeEmbedding::from_correspondence(host, pattern, atom_map))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_core::SubgraphIsomorphismAlgorithm::{
        ArcMatch, RayKirsch, Ri, Ullmann, Vf2, Vf2Rdkit,
    };
    use umol_graph_core::ARCMATCH_DEFAULT_PATH_LENGTH;
    use umol_shared::element::Element;

    use super::super::atom::AtomAst;
    use super::super::bond::BondAst;
    use super::super::id::AtomId;
    use super::super::molecule::MoleculeAst;
    use super::SubstructureMatchAlgorithm::GraphAndOverlays;

    // Pattern C-C in host C-C-O: matches the C-C edge in both orientations, never the
    // C-O edge. The same occurrence set under every subgraph-isomorphism algorithm.
    #[rstest]
    fn test_molecule_ast_substructure_matches() {
        let host = MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::O),
            ],
            vec![
                (AtomId(0), AtomId(1), BondAst::default()),
                (AtomId(1), AtomId(2), BondAst::default()),
            ],
        );
        let pattern = MoleculeAst::from_atoms_and_bonds(
            vec![
                AtomAst::from_element(Element::C),
                AtomAst::from_element(Element::C),
            ],
            vec![(AtomId(0), AtomId(1), BondAst::default())],
        );
        for subiso in [
            Vf2,
            Ullmann,
            Ri,
            ArcMatch { path_length: ARCMATCH_DEFAULT_PATH_LENGTH },
            Vf2Rdkit,
            RayKirsch,
        ] {
            let mut maps: Vec<Vec<AtomId>> = pattern
                .substructure_matches(&host, GraphAndOverlays, subiso)
                .iter()
                .map(|e| e.host_atoms().to_vec())
                .collect();
            maps.sort();
            assert_eq!(
                maps,
                vec![vec![AtomId(0), AtomId(1)], vec![AtomId(1), AtomId(0)]],
                "subiso {subiso:?}"
            );
        }
    }
}
