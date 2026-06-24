//! Typed matching result wrapper over `umol_graph_core::Matching`.

use umol_graph_core::{Matching, NodeId};

use super::id::{AtomId, BondId};

/// Bond-level wrapper over `umol_graph_core::Matching`. Exposes matched
/// bonds and matched-atom membership in terms of `AtomId` / `BondId`.
#[derive(Clone, Debug)]
pub struct BondMatching(pub(crate) Matching);

impl BondMatching {
    pub fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
        self.0.edges().iter().map(|&e| BondId::from(e))
    }

    pub fn size(&self) -> usize {
        self.0.size()
    }

    pub fn is_perfect(&self, atom_count: usize) -> bool {
        self.0.is_perfect(atom_count)
    }

    pub fn mate(&self, atom: AtomId) -> Option<AtomId> {
        self.0.mate(NodeId::from(atom)).map(AtomId::from)
    }

    pub fn is_matched(&self, atom: AtomId) -> bool {
        self.0.is_matched(NodeId::from(atom))
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::MaxMatchingAlgorithm;
    use umol_chem::element::Element;

    use super::*;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::molecule::MoleculeAst;

    fn chain(n: usize) -> MoleculeAst {
        let atoms = vec![AtomAst::from_element(Element::C); n];
        let bonds: Vec<_> = (0..n.saturating_sub(1))
            .map(|i| {
                (
                    AtomId(i as u32),
                    AtomId((i + 1) as u32),
                    BondAst::from_order(1),
                )
            })
            .collect();
        MoleculeAst::from_atoms_and_bonds(atoms, bonds)
    }

    fn ring(n: usize) -> MoleculeAst {
        let atoms = vec![AtomAst::from_element(Element::C); n];
        let bonds: Vec<_> = (0..n)
            .map(|i| {
                (
                    AtomId(i as u32),
                    AtomId(((i + 1) % n) as u32),
                    BondAst::from_order(1),
                )
            })
            .collect();
        MoleculeAst::from_atoms_and_bonds(atoms, bonds)
    }

    #[fixture]
    fn chain_4_matching() -> BondMatching {
        chain(4)
            .graph()
            .maximum_matching(MaxMatchingAlgorithm::Edmonds)
    }

    #[fixture]
    fn ring_6_matching() -> BondMatching {
        ring(6)
            .graph()
            .maximum_matching(MaxMatchingAlgorithm::Edmonds)
    }

    #[rstest]
    #[case::chain_4(chain_4_matching(), 2)]
    #[case::ring_6(ring_6_matching(), 3)]
    fn test_bond_matching_size(#[case] matching: BondMatching, #[case] expected: usize) {
        assert_eq!(matching.size(), expected);
    }

    #[rstest]
    fn test_bond_matching_bonds_enumerates_matched_edges(chain_4_matching: BondMatching) {
        let bonds: Vec<BondId> = chain_4_matching.bonds().collect();
        assert_eq!(bonds.len(), 2);
        // All bond indices are < bond count of chain(4) == 3.
        for b in bonds {
            assert!(b.0 < 3);
        }
    }

    #[rstest]
    #[case::chain_4_is_perfect(chain_4_matching(), 4, true)]
    #[case::chain_4_wrong_atom_count(chain_4_matching(), 5, false)]
    #[case::ring_6_is_perfect(ring_6_matching(), 6, true)]
    fn test_bond_matching_is_perfect(
        #[case] matching: BondMatching,
        #[case] atom_count: usize,
        #[case] expected: bool,
    ) {
        assert_eq!(matching.is_perfect(atom_count), expected);
    }

    #[rstest]
    fn test_bond_matching_mate_and_is_matched(chain_4_matching: BondMatching) {
        for i in 0..4 {
            assert!(chain_4_matching.is_matched(AtomId(i)));
            let mate = chain_4_matching.mate(AtomId(i)).unwrap();
            assert_ne!(mate, AtomId(i));
            // Matching is symmetric.
            assert_eq!(chain_4_matching.mate(mate), Some(AtomId(i)));
        }
    }

    #[rstest]
    fn test_bond_matching_mate_unmatched_atom() {
        // Single atom is not matched; a "matching" on just {0} is size 0 and
        // atom 0 has no mate.
        let ast = chain(1);
        let m = ast.graph().maximum_matching(MaxMatchingAlgorithm::Edmonds);
        assert!(!m.is_matched(AtomId(0)));
        assert_eq!(m.mate(AtomId(0)), None);
    }
}
