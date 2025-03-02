// GraphMolecule implementation

use super::types::{AtomIndex, BondIndex};
use super::{bond::BondOrder, GraphAtom, GraphBond};
use crate::atom::AtomSite;
use crate::error::MoleculeError;
use petgraph::stable_graph::StableGraph;
use std::collections::HashMap;
use std::fmt::{self, Display};

#[derive(Debug, Clone)]
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

    // pub fn with_atom_charge(self, idx: AtomIndex, charge: i8) -> Result<Self, ValidationError> {
    //     let mut graph = self.graph;
    //     let atom = graph.node_weight_mut(idx).unwrap().with_charge(charge);
    //     // Validation logic
    //     Ok(Self {
    //         graph,
    //         properties: self.properties,
    //     })
    // }

    // Other mutation functions

    // Basic accessors
    pub fn atom_count(&self) -> usize {
        self.graph.node_count()
    }

    pub fn bond_count(&self) -> usize {
        self.graph.edge_count()
    }

    // Atom access
    pub fn atoms(&self) -> impl Iterator<Item = (AtomIndex, &GraphAtom)> {
        self.graph
            .node_indices()
            .map(move |idx| (idx.into(), self.graph.node_weight(idx).unwrap()))
    }

    pub fn atom(&self, idx: AtomIndex) -> Option<&GraphAtom> {
        self.graph.node_weight(idx.into())
    }

    // Bond access
    pub fn bonds(&self) -> impl Iterator<Item = (BondIndex, &GraphBond)> {
        self.graph
            .edge_indices()
            .map(move |idx| (idx.into(), self.graph.edge_weight(idx).unwrap()))
    }

    pub fn bond(&self, idx: BondIndex) -> Option<&GraphBond> {
        self.graph.edge_weight(idx.into())
    }

    // Connectivity
    pub fn neighbors(&self, idx: AtomIndex) -> Vec<(AtomIndex, BondIndex)> {
        self.graph
            .neighbors(idx.into())
            .map(|neighbor_idx| {
                let edge_idx = self.graph.find_edge(idx.into(), neighbor_idx).unwrap();
                (neighbor_idx.into(), edge_idx.into())
            })
            .collect()
    }

    pub fn atoms_connected_by_bond(&self, bond_idx: BondIndex) -> Option<(AtomIndex, AtomIndex)> {
        self.graph
            .edge_endpoints(bond_idx.into())
            .map(|(a, b)| (a.into(), b.into()))
    }

    // Temporary placeholder for formula calculation
    pub fn formula(&self) -> String {
        // This is a simplified implementation
        // A real implementation would count atoms by element and format properly
        let mut elements = HashMap::new();

        for (_, atom) in self.atoms() {
            if let Some(element) = atom.element() {
                *elements.entry(element.symbol()).or_insert(0) += 1;
            }
        }

        // Sort elements with C and H first, then alphabetically
        let mut element_counts: Vec<_> = elements.into_iter().collect();
        element_counts.sort_by(|a, b| {
            if a.0 == "C" {
                return std::cmp::Ordering::Less;
            }
            if b.0 == "C" {
                return std::cmp::Ordering::Greater;
            }
            if a.0 == "H" {
                return std::cmp::Ordering::Less;
            }
            if b.0 == "H" {
                return std::cmp::Ordering::Greater;
            }
            a.0.cmp(b.0)
        });

        // Format the formula
        element_counts
            .into_iter()
            .map(|(symbol, count)| {
                if count == 1 {
                    symbol.to_string()
                } else {
                    format!("{}{}", symbol, count)
                }
            })
            .collect::<Vec<_>>()
            .join("")
    }
}

impl Display for GraphMolecule {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Header with formula and atom count
        writeln!(
            f,
            "Molecule: {} ({} atoms)",
            self.formula(),
            self.atom_count()
        )?;

        // Connection table
        writeln!(f, "\nAtoms:")?;
        for (idx, atom) in self.atoms() {
            let neighbors = self.neighbors(idx);
            let connections: Vec<String> = neighbors
                .iter()
                .map(|(neighbor_idx, bond_idx)| {
                    let bond = self.bond(*bond_idx).unwrap();
                    let neighbor = self.atom(*neighbor_idx).unwrap();
                    format!(
                        "{}{}",
                        neighbor.element().unwrap().symbol(),
                        match bond.order() {
                            BondOrder::Single => "-",
                            BondOrder::Double => "=",
                            BondOrder::Triple => "#",
                            BondOrder::Quadruple => "$",
                        }
                    )
                })
                .collect();

            writeln!(
                f,
                "{:3}: {} connected to: {}",
                idx,
                atom,
                connections.join(", ")
            )?;
        }

        Ok(())
    }
}

pub struct GraphMoleculeBuilder {
    graph: StableGraph<GraphAtom, GraphBond>,
}

impl GraphMoleculeBuilder {
    pub fn new() -> Self {
        Self {
            graph: StableGraph::new(),
        }
    }

    pub fn add_atom<T: Into<GraphAtom>>(&mut self, atom: T) -> AtomIndex {
        self.graph.add_node(atom.into()).into()
    }

    pub fn add_bond(
        &mut self,
        from: AtomIndex,
        to: AtomIndex,
        order: BondOrder,
    ) -> Result<BondIndex, MoleculeError> {
        if !self.graph.contains_node(from.into()) {
            return Err(MoleculeError::InvalidAtomIndex(from));
        }
        if !self.graph.contains_node(to.into()) {
            return Err(MoleculeError::InvalidAtomIndex(to));
        }

        let bond = GraphBond::new(order);
        Ok(self
            .graph
            .add_edge(from.into(), to.into(), bond.into())
            .into())
    }

    pub fn build(self) -> GraphMolecule {
        GraphMolecule {
            graph: self.graph,
            properties: HashMap::new(),
        }
    }
}
