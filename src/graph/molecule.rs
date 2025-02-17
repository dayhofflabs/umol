// GraphMolecule implementation

use super::types::{AtomIndex, BondIndex};
use crate::atom::Element;
use crate::error::MoleculeError;
use petgraph::stable_graph::StableGraph;
use std::collections::HashMap;

use super::{bond::BondOrder, GraphAtom, GraphBond};

pub struct GraphMolecule {
    graph: StableGraph<GraphAtom, GraphBond>,
    properties: HashMap<String, String>,
}

impl GraphMolecule {
    pub fn new() -> Self {
        Self {
            graph: StableGraph::new(),
            properties: HashMap::new(),
        }
    }

    pub fn validate(&self) -> Result<(), MoleculeError> {
        Ok(())
    }

    pub fn add_atom(&mut self, element: Element) -> AtomIndex {
        self.graph.add_node(GraphAtom::new(element))
    }

    pub fn add_bond(
        &mut self,
        from: AtomIndex,
        to: AtomIndex,
        order: BondOrder,
    ) -> Result<BondIndex, MoleculeError> {
        if !self.graph.contains_node(from) || !self.graph.contains_node(to) {
            return Err(MoleculeError::InvalidAtomIndex);
        }

        let bond = GraphBond::new(from, to, order);
        Ok(self.graph.add_edge(from, to, bond))
    }
}
