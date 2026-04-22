//! Typed automorphism result wrapper over `umol_graph_core::Automorphism`.

use umol_graph_core::{AutoGroupOrder, Automorphism, NodeId};

use super::idx::AtomIdx;

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
