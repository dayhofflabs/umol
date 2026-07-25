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
    use umol_chem::element::Element;
    use umol_graph_core::GeneralMaximumMatchingAlgorithm;

    use super::*;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::molecule::{MoleculeAst, MoleculeParts};

    #[fixture]
    fn chain_4_matching() -> BondMatching {
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C); 4],
            bonds: vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
                (AtomId(2), AtomId(3), BondAst::from_order(1)),
            ],
            ..Default::default()
        })
        .graph()
        .general_maximum_matching(
            &[AtomId(0), AtomId(1), AtomId(2), AtomId(3)],
            GeneralMaximumMatchingAlgorithm::Edmonds,
        )
    }

    #[fixture]
    fn ring_6_matching() -> BondMatching {
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C); 6],
            bonds: vec![
                (AtomId(0), AtomId(1), BondAst::from_order(1)),
                (AtomId(1), AtomId(2), BondAst::from_order(1)),
                (AtomId(2), AtomId(3), BondAst::from_order(1)),
                (AtomId(3), AtomId(4), BondAst::from_order(1)),
                (AtomId(4), AtomId(5), BondAst::from_order(1)),
                (AtomId(5), AtomId(0), BondAst::from_order(1)),
            ],
            ..Default::default()
        })
        .graph()
        .general_maximum_matching(
            &[
                AtomId(0),
                AtomId(1),
                AtomId(2),
                AtomId(3),
                AtomId(4),
                AtomId(5),
            ],
            GeneralMaximumMatchingAlgorithm::Edmonds,
        )
    }

    #[fixture]
    fn singleton_matching() -> BondMatching {
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst::from_element(Element::C)],
            ..Default::default()
        })
        .graph()
        .general_maximum_matching(&[AtomId(0)], GeneralMaximumMatchingAlgorithm::Edmonds)
    }

    #[rstest]
    #[case::chain_4(chain_4_matching(), 2)]
    #[case::ring_6(ring_6_matching(), 3)]
    fn test_bond_matching_size(#[case] matching: BondMatching, #[case] expected: usize) {
        assert_eq!(matching.size(), expected);
    }

    #[rstest]
    fn test_bond_matching_bonds(chain_4_matching: BondMatching) {
        assert_eq!(
            chain_4_matching.bonds().collect::<Vec<_>>(),
            vec![BondId(0), BondId(2)],
        );
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
    #[case::first(chain_4_matching(), AtomId(0), Some(AtomId(1)))]
    #[case::second(chain_4_matching(), AtomId(1), Some(AtomId(0)))]
    #[case::third(chain_4_matching(), AtomId(2), Some(AtomId(3)))]
    #[case::fourth(chain_4_matching(), AtomId(3), Some(AtomId(2)))]
    #[case::unmatched(singleton_matching(), AtomId(0), None)]
    fn test_bond_matching_mate(
        #[case] matching: BondMatching,
        #[case] atom: AtomId,
        #[case] expected: Option<AtomId>,
    ) {
        assert_eq!(matching.mate(atom), expected);
    }

    #[rstest]
    #[case::matched(chain_4_matching(), AtomId(0), true)]
    #[case::unmatched(singleton_matching(), AtomId(0), false)]
    fn test_bond_matching_is_matched(
        #[case] matching: BondMatching,
        #[case] atom: AtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(matching.is_matched(atom), expected);
    }
}
