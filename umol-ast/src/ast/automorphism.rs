//! Typed automorphism result wrapper over `umol_graph_core::Automorphism`.

use umol_graph_core::{AutoGroupOrder, Automorphism, NodeId};

use super::idx::AtomId;

/// Atom-level wrapper over `umol_graph_core::Automorphism`. Indexes the
/// permutation in terms of `AtomId` rather than raw `NodeId`.
#[derive(Clone, Debug)]
pub struct AtomAutomorphism(pub(crate) Automorphism);

impl AtomAutomorphism {
    pub fn atom_count(&self) -> usize {
        self.0.node_count()
    }

    pub fn orbit_count(&self) -> usize {
        self.0.orbit_count()
    }

    pub fn orbit_of(&self, atom: AtomId) -> AtomId {
        AtomId::from(self.0.orbit_of(NodeId::from(atom)))
    }

    pub fn same_orbit(&self, a: AtomId, b: AtomId) -> bool {
        self.0.same_orbit(NodeId::from(a), NodeId::from(b))
    }

    pub fn canonical_labeling(&self) -> Vec<AtomId> {
        self.0
            .canonical_labeling()
            .iter()
            .map(|&n| AtomId::from(n))
            .collect()
    }

    pub fn auto_group_order(&self) -> AutoGroupOrder {
        self.0.auto_group_order()
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_graph_core::AutomorphismAlgorithm;
    use umol_shared::element::Element;

    use super::*;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    
    use crate::ast::idx::AtomId;
    use crate::ast::molecule::MoleculeAst;

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
        MoleculeAst::from_atoms_and_bonds(
            atoms,
            bonds,
        )
    }

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
        MoleculeAst::from_atoms_and_bonds(
            atoms,
            bonds,
        )
    }

    #[fixture]
    fn hexagon() -> AtomAutomorphism {
        ring(6).graph().automorphisms(|_| 0u8, AutomorphismAlgorithm::Nauty)
    }

    #[fixture]
    fn chain_3() -> AtomAutomorphism {
        chain(3).graph().automorphisms(|_| 0u8, AutomorphismAlgorithm::Nauty)
    }

    #[rstest]
    #[case::hexagon_wraps_six(hexagon(), 6, 1)]
    #[case::chain_3_wraps_three(chain_3(), 3, 2)]
    fn test_atom_automorphism_counts(
        #[case] auto: AtomAutomorphism,
        #[case] expected_atom_count: usize,
        #[case] expected_orbit_count: usize,
    ) {
        assert_eq!(auto.atom_count(), expected_atom_count);
        assert_eq!(auto.orbit_count(), expected_orbit_count);
    }

    #[rstest]
    fn test_atom_automorphism_orbit_of_agrees_within_orbit(hexagon: AtomAutomorphism) {
        // In a hexagon every atom is in the same orbit, so `orbit_of` must
        // return the same representative for every atom.
        let rep = hexagon.orbit_of(AtomId(0));
        for i in 1..6 {
            assert_eq!(hexagon.orbit_of(AtomId(i)), rep);
        }
    }

    #[rstest]
    fn test_atom_automorphism_orbit_of_distinguishes_non_equivalent(chain_3: AtomAutomorphism) {
        // In C-C-C the center atom sits in a different orbit from the endpoints.
        assert_ne!(chain_3.orbit_of(AtomId(0)), chain_3.orbit_of(AtomId(1)));
        assert_eq!(chain_3.orbit_of(AtomId(0)), chain_3.orbit_of(AtomId(2)));
    }

    #[rstest]
    #[case::hexagon_equivalent(hexagon(), AtomId(0), AtomId(3), true)]
    #[case::chain_endpoints_equivalent(chain_3(), AtomId(0), AtomId(2), true)]
    #[case::chain_endpoint_vs_middle(chain_3(), AtomId(0), AtomId(1), false)]
    fn test_atom_automorphism_same_orbit(
        #[case] auto: AtomAutomorphism,
        #[case] a: AtomId,
        #[case] b: AtomId,
        #[case] expected: bool,
    ) {
        assert_eq!(auto.same_orbit(a, b), expected);
    }

    #[rstest]
    fn test_atom_automorphism_canonical_labeling_is_permutation(hexagon: AtomAutomorphism) {
        let labeling = hexagon.canonical_labeling();
        assert_eq!(labeling.len(), 6);
        let mut sorted = labeling.clone();
        sorted.sort_unstable();
        assert_eq!(
            sorted,
            (0..6).map(|i| AtomId(i as u32)).collect::<Vec<_>>()
        );
    }

    #[rstest]
    fn test_atom_automorphism_auto_group_order_hexagon(hexagon: AtomAutomorphism) {
        // D6 has 12 automorphisms on the vertex set.
        assert_eq!(hexagon.auto_group_order(), AutoGroupOrder::Exact(12));
    }

    #[rstest]
    fn test_atom_automorphism_auto_group_order_chain(chain_3: AtomAutomorphism) {
        // Linear C-C-C: only identity and the endpoint-swap reflection.
        assert_eq!(chain_3.auto_group_order(), AutoGroupOrder::Exact(2));
    }
}
