// GraphMolecule implementation

use super::types::{AtomIndex, BondIndex};
use super::{bond::BondOrder, GraphAtom, GraphBond};
use crate::{AtomSite, Error};
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

    pub fn add_bond<T: Into<GraphBond>>(
        &mut self,
        from: AtomIndex,
        to: AtomIndex,
        order: T,
    ) -> Result<BondIndex, Error> {
        if !self.graph.contains_node(from.into()) {
            return Err(Error::InvalidAtomIndex(from));
        }
        if !self.graph.contains_node(to.into()) {
            return Err(Error::InvalidAtomIndex(to));
        }

        let bond = order.into();
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

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::atom::Element;
//     use crate::graph::atom::GraphAtom;
//     use crate::graph::bond::{BondOrder, GraphBond};

//     #[test]
//     fn test_new_molecule() {
//         let mol = GraphMolecule::new();
//         assert_eq!(mol.atom_count(), 0);
//         assert_eq!(mol.bond_count(), 0);
//     }

//     #[test]
//     fn test_molecule_builder() {
//         // Test building a simple molecule with the builder
//         let mol = GraphMoleculeBuilder::new()
//             .add_atom(GraphAtom::new(Element::C))
//             .add_atom(GraphAtom::new(Element::O))
//             .add_bond(0.into(), 1.into(), GraphBond::new(BondOrder::Double))
//             .unwrap()
//             .build();

//         assert_eq!(mol.atom_count(), 2);
//         assert_eq!(mol.bond_count(), 1);

//         // Verify atom types
//         let atoms: Vec<_> = mol.atoms().collect();
//         for (idx, atom) in &atoms {
//             if idx.index() == 0 {
//                 assert_eq!(atom.element(), Some(Element::C));
//             } else if idx.index() == 1 {
//                 assert_eq!(atom.element(), Some(Element::O));
//             }
//         }

//         // Verify bond
//         let bonds: Vec<_> = mol.bonds().collect();
//         assert_eq!(bonds.len(), 1);
//         assert_eq!(bonds[0].1.order(), BondOrder::Double);
//     }

//     #[test]
//     fn test_atom_access() {
//         // Build a test molecule
//         let mol = GraphMoleculeBuilder::new()
//             .add_atom(GraphAtom::new(Element::C))
//             .add_atom(GraphAtom::new(Element::O))
//             .build();

//         // Test atom retrieval by index
//         let c_idx = AtomIndex::new(0);
//         let o_idx = AtomIndex::new(1);

//         let c_atom = mol.atom(c_idx).unwrap();
//         let o_atom = mol.atom(o_idx).unwrap();

//         assert_eq!(c_atom.element(), Some(Element::C));
//         assert_eq!(o_atom.element(), Some(Element::O));

//         // Test invalid atom index
//         let invalid_idx = AtomIndex::new(999);
//         assert!(mol.atom(invalid_idx).is_none());
//     }

//     #[test]
//     fn test_bond_access() {
//         // Build a test molecule with a bond
//         let mut builder = GraphMoleculeBuilder::new();
//         let c_idx = builder.add_atom(GraphAtom::new(Element::C));
//         let o_idx = builder.add_atom(GraphAtom::new(Element::O));
//         let bond_idx = builder.add_bond(c_idx, o_idx, GraphBond::new(BondOrder::Double)).unwrap();
//         let mol = builder.build();

//         // Test bond retrieval
//         let bond = mol.bond(bond_idx).unwrap();
//         assert_eq!(bond.order(), BondOrder::Double);

//         // Test atoms connected by bond
//         let (atom1, atom2) = mol.atoms_connected_by_bond(bond_idx).unwrap();
//         assert!(
//             (atom1 == c_idx && atom2 == o_idx) ||
//             (atom1 == o_idx && atom2 == c_idx)
//         );

//         // Test invalid bond index
//         let invalid_bond_idx = BondIndex::new(999);
//         assert!(mol.bond(invalid_bond_idx).is_none());
//         assert!(mol.atoms_connected_by_bond(invalid_bond_idx).is_none());
//     }

//     #[test]
//     fn test_neighbors() {
//         // Build a test molecule with a chain C-O-N
//         let mut builder = GraphMoleculeBuilder::new();
//         let c_idx = builder.add_atom(GraphAtom::new(Element::C));
//         let o_idx = builder.add_atom(GraphAtom::new(Element::O));
//         let n_idx = builder.add_atom(GraphAtom::new(Element::N));

//         builder.add_bond(c_idx, o_idx, GraphBond::new(BondOrder::Single)).unwrap();
//         builder.add_bond(o_idx, n_idx, GraphBond::new(BondOrder::Single)).unwrap();

//         let mol = builder.build();

//         // Test neighbors
//         let c_neighbors = mol.neighbors(c_idx);
//         assert_eq!(c_neighbors.len(), 1);
//         assert_eq!(c_neighbors[0].0, o_idx);

//         let o_neighbors = mol.neighbors(o_idx);
//         assert_eq!(o_neighbors.len(), 2);
//         assert!(o_neighbors.iter().any(|(idx, _)| *idx == c_idx));
//         assert!(o_neighbors.iter().any(|(idx, _)| *idx == n_idx));

//         // Test neighbors of invalid atom
//         let invalid_idx = AtomIndex::new(999);
//         assert!(mol.neighbors(invalid_idx).is_empty());
//     }

//     #[test]
//     fn test_atom_iteration() {
//         // Build a test molecule
//         let mut builder = GraphMoleculeBuilder::new();
//         let c_idx = builder.add_atom(GraphAtom::new(Element::C));
//         let o_idx = builder.add_atom(GraphAtom::new(Element::O));
//         let n_idx = builder.add_atom(GraphAtom::new(Element::N));
//         let mol = builder.build();

//         // Test atom iteration
//         let atoms: Vec<_> = mol.atoms().collect();
//         assert_eq!(atoms.len(), 3);

//         // Check that all atoms are present
//         let indices: Vec<_> = atoms.iter().map(|(idx, _)| *idx).collect();
//         assert!(indices.contains(&c_idx));
//         assert!(indices.contains(&o_idx));
//         assert!(indices.contains(&n_idx));

//         // Check that atoms match their indices
//         for (idx, atom) in atoms {
//             if idx == c_idx {
//                 assert_eq!(atom.element(), Some(Element::C));
//             } else if idx == o_idx {
//                 assert_eq!(atom.element(), Some(Element::O));
//             } else if idx == n_idx {
//                 assert_eq!(atom.element(), Some(Element::N));
//             }
//         }
//     }

//     #[test]
//     fn test_bond_iteration() {
//         // Build a test molecule with bonds
//         let mut builder = GraphMoleculeBuilder::new();
//         let c_idx = builder.add_atom(GraphAtom::new(Element::C));
//         let o_idx = builder.add_atom(GraphAtom::new(Element::O));
//         let n_idx = builder.add_atom(GraphAtom::new(Element::N));

//         let co_bond = builder.add_bond(c_idx, o_idx, GraphBond::new(BondOrder::Single)).unwrap();
//         let on_bond = builder.add_bond(o_idx, n_idx, GraphBond::new(BondOrder::Double)).unwrap();

//         let mol = builder.build();

//         // Test bond iteration
//         let bonds: Vec<_> = mol.bonds().collect();
//         assert_eq!(bonds.len(), 2);

//         // Check that all bonds are present with correct order
//         let bond_indices: Vec<_> = bonds.iter().map(|(idx, _)| *idx).collect();
//         assert!(bond_indices.contains(&co_bond));
//         assert!(bond_indices.contains(&on_bond));

//         for (idx, bond) in bonds {
//             if idx == co_bond {
//                 assert_eq!(bond.order(), BondOrder::Single);
//             } else if idx == on_bond {
//                 assert_eq!(bond.order(), BondOrder::Double);
//             }
//         }
//     }

//     #[test]
//     fn test_formula() {
//         // Test formula generation for a simple molecule
//         let mol = GraphMoleculeBuilder::new()
//             .add_atom(GraphAtom::new(Element::C))
//             .add_atom(GraphAtom::new(Element::O))
//             .add_atom(GraphAtom::new(Element::O))
//             .add_atom(GraphAtom::new(Element::H))
//             .add_atom(GraphAtom::new(Element::H))
//             .build();

//         assert_eq!(mol.formula(), "CH2O2");
//     }

//     #[test]
//     fn test_display() {
//         // Test string representation of a molecule
//         let mut builder = GraphMoleculeBuilder::new();
//         let c_idx = builder.add_atom(GraphAtom::new(Element::C));
//         let o_idx = builder.add_atom(GraphAtom::new(Element::O));
//         builder.add_bond(c_idx, o_idx, GraphBond::new(BondOrder::Double)).unwrap();

//         let mol = builder.build();
//         let display_str = format!("{}", mol);

//         // Basic checks that the display output contains expected information
//         assert!(display_str.contains("CO"));
//         assert!(display_str.contains("2 atoms"));
//         assert!(display_str.contains("C"));
//         assert!(display_str.contains("O"));
//     }
// }
