//! Typed automorphism result wrapper over `umol_graph_core::Automorphism`.

use umol_graph_core::{AutoGroupOrder, Automorphism, NodeId};

use super::idx::AtomIdx;

/// Atom-level wrapper over `umol_graph_core::Automorphism`. Indexes the
/// permutation in terms of `AtomIdx` rather than raw `NodeId`.
#[derive(Clone, Debug)]
pub struct AtomAutomorphism(pub(crate) Automorphism);

impl AtomAutomorphism {
    pub fn atom_count(&self) -> usize {
        self.0.node_count()
    }

    pub fn num_orbits(&self) -> usize {
        self.0.num_orbits()
    }

    pub fn orbit_of(&self, atom: AtomIdx) -> AtomIdx {
        AtomIdx::from(self.0.orbit_of(NodeId::from(atom)))
    }

    pub fn same_orbit(&self, a: AtomIdx, b: AtomIdx) -> bool {
        self.0.same_orbit(NodeId::from(a), NodeId::from(b))
    }

    pub fn canonical_labeling(&self) -> Vec<AtomIdx> {
        self.0
            .canonical_labeling()
            .iter()
            .map(|&n| AtomIdx::from(n))
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
    
    use crate::ast::idx::AtomIdx;
    use crate::ast::molecule::MoleculeAst;

    fn ring(n: usize) -> MoleculeAst {
        let atoms = vec![AtomAst::from_element(Element::C); n];
        let bonds: Vec<_> = (0..n)
            .map(|i| {
                (
                    AtomIdx(i as u32),
                    AtomIdx(((i + 1) % n) as u32),
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
                    AtomIdx(i as u32),
                    AtomIdx((i + 1) as u32),
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
        ring(6).automorphisms(|_| 0u8, AutomorphismAlgorithm::Nauty)
    }

    #[fixture]
    fn chain_3() -> AtomAutomorphism {
        chain(3).automorphisms(|_| 0u8, AutomorphismAlgorithm::Nauty)
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
        assert_eq!(auto.num_orbits(), expected_orbit_count);
    }

    #[rstest]
    fn test_atom_automorphism_orbit_of_agrees_within_orbit(hexagon: AtomAutomorphism) {
        // In a hexagon every atom is in the same orbit, so `orbit_of` must
        // return the same representative for every atom.
        let rep = hexagon.orbit_of(AtomIdx(0));
        for i in 1..6 {
            assert_eq!(hexagon.orbit_of(AtomIdx(i)), rep);
        }
    }

    #[rstest]
    fn test_atom_automorphism_orbit_of_distinguishes_non_equivalent(chain_3: AtomAutomorphism) {
        // In C-C-C the center atom sits in a different orbit from the endpoints.
        assert_ne!(chain_3.orbit_of(AtomIdx(0)), chain_3.orbit_of(AtomIdx(1)));
        assert_eq!(chain_3.orbit_of(AtomIdx(0)), chain_3.orbit_of(AtomIdx(2)));
    }

    #[rstest]
    #[case::hexagon_equivalent(hexagon(), AtomIdx(0), AtomIdx(3), true)]
    #[case::chain_endpoints_equivalent(chain_3(), AtomIdx(0), AtomIdx(2), true)]
    #[case::chain_endpoint_vs_middle(chain_3(), AtomIdx(0), AtomIdx(1), false)]
    fn test_atom_automorphism_same_orbit(
        #[case] auto: AtomAutomorphism,
        #[case] a: AtomIdx,
        #[case] b: AtomIdx,
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
            (0..6).map(|i| AtomIdx(i as u32)).collect::<Vec<_>>()
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
