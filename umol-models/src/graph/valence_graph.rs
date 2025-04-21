//! Valence graph model
//!
//! Graph model of atoms and bonds, with valence constraints

use crate::graph::{ValenceAtom, ValenceBond};
use petgraph::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fmt;
use umol::error::DataError;
use umol::{Error, Result};

/// The type used for original indices in the input vectors. Must be sortable, hashable, and copyable.
pub type IndexType = u16;
/// The type used for internal atom indices.
pub type AtomIndex = NodeIndex<IndexType>;
/// The type used for internal bond indices.
pub type BondIndex = EdgeIndex<IndexType>;

/// Graph model of atoms and bonds, with valence constraints
#[derive(Debug, Clone)]
pub struct ValenceGraph {
    data: StableGraph<ValenceAtom, ValenceBond, Undirected, IndexType>,
}

impl ValenceGraph {
    /// Creates a ValenceGraph from vectors of atoms and bonds with arbitrary original indices.
    ///
    /// Atom and bond lists are processed to map potentially non-sequential, 1-based,
    /// or otherwise arbitrary indices (`I`) to dense, 0-based internal indices used
    /// by the underlying graph structure (`AtomIndex`, `BondIndex`).
    ///
    /// # Arguments
    /// * \`atoms\`: A vector of tuples \`(original_index, atom_data)\`.
    /// * \`bonds\`: A vector of tuples \`(original_bond_index, original_atom_index1, original_atom_index2, bond_data)\`.
    ///
    /// # Type Parameters
    /// * \`I\`: The type used for original indices in the input vectors. Must be sortable, hashable, and copyable.
    ///
    /// # Returns
    /// A \`Result\` containing the new \`ValenceGraph\` on success, or a \`GraphBuildError\`
    /// if issues like duplicate or missing atom indices are found.
    ///
    /// # Behavior
    /// 1. Checks for duplicate original atom indices.
    /// 2. Sorts atoms by their original index.
    /// 3. Adds atoms to the graph, creating a mapping from original index (\`I\`) to internal \`AtomIndex\`.
    /// 4. Sorts bonds by their original index.
    /// 5. Adds bonds to the graph using the mapped internal \`AtomIndex\` values.
    pub fn from_atoms_bonds(
        mut atoms: Vec<(IndexType, ValenceAtom)>,
        mut bonds: Vec<(IndexType, IndexType, IndexType, ValenceBond)>,
    ) -> Self {
        let mut graph = StableGraph::with_capacity(atoms.len(), bonds.len());
        let mut idx_map = HashMap::with_capacity(atoms.len());
        let mut seen_indices = HashSet::with_capacity(atoms.len());

        // 1. Check for duplicate original atom indices before sorting
        for (old_idx, _) in atoms.iter() {
            if !seen_indices.insert(*old_idx) {
                panic!("Duplicate atom index: {}", old_idx);
            }
        }

        // 2. Sort atoms by original index
        atoms.sort_unstable_by_key(|(k, _)| *k);

        // 3. Add atoms and build the index map
        for (old_idx, atom) in atoms {
            let new_idx = graph.add_node(atom);
            idx_map.insert(old_idx, new_idx);
        }

        // 4. Sort bonds by original index
        bonds.sort_unstable_by_key(|(k, _, _, _)| *k);

        // 5. Add bonds using the mapped indices
        for (_old_idx, old_idx1, old_idx2, bond) in bonds {
            // Look up the internal indices corresponding to the original atom indices
            let new_idx1 = idx_map
                .get(&old_idx1)
                .unwrap_or_else(|| panic!("Missing atom index: {}", old_idx1));
            let new_idx2 = idx_map
                .get(&old_idx2)
                .unwrap_or_else(|| panic!("Missing atom index: {}", old_idx2));

            // Add the edge using the *internal* indices
            graph.add_edge(*new_idx1, *new_idx2, bond);
            // We ignore the EdgeIndex returned by add_edge here, as petgraph manages it.
        }

        Self { data: graph }
    }

    /// Tries to create a ValenceGraph from vectors of atoms and bonds with arbitrary original indices.
    ///
    /// Atom and bond lists are processed to map potentially non-sequential, 1-based,
    /// or otherwise arbitrary indices (\`IndexType\`) to dense, 0-based internal indices used
    /// by the underlying graph structure (\`AtomIndex\`, \`BondIndex\`).
    ///
    /// # Arguments
    /// * \`atoms\`: A vector of tuples \`(original_index, atom_data)\`.
    /// * \`bonds\`: A vector of tuples \`(original_bond_index, original_atom_index1, original_atom_index2, bond_data)\`.
    ///
    /// # Returns
    /// A \`Result\` containing the new \`ValenceGraph\` on success, or a \`umol::Error\`
    /// wrapping a \`DataError\` if issues like duplicate or missing atom indices are found.
    ///
    /// # Behavior
    /// 1. Checks for duplicate original atom indices.
    /// 2. Sorts atoms by their original index.
    /// 3. Adds atoms to the graph, creating a mapping from original index (\`IndexType\`) to internal \`AtomIndex\`.
    /// 4. Sorts bonds by their original index.
    /// 5. Adds bonds to the graph using the mapped internal \`AtomIndex\` values.
    pub fn try_from_atoms_bonds(
        mut atoms: Vec<(IndexType, ValenceAtom)>,
        mut bonds: Vec<(IndexType, IndexType, IndexType, ValenceBond)>,
    ) -> Result<Self> {
        let mut graph = StableGraph::with_capacity(atoms.len(), bonds.len());
        let mut idx_map = HashMap::with_capacity(atoms.len());
        let mut seen_indices = HashSet::with_capacity(atoms.len());

        // 1. Check for duplicate original atom indices before sorting
        for (old_idx, _) in atoms.iter() {
            if !seen_indices.insert(*old_idx) {
                // Return error on duplicate
                return Err(DataError::DuplicateAtomIndex(*old_idx as usize).into());
            }
        }

        // 2. Sort atoms by original index
        atoms.sort_unstable_by_key(|(k, _)| *k);

        // 3. Add atoms and build the index map
        for (old_idx, atom) in atoms {
            let new_idx = graph.add_node(atom);
            idx_map.insert(old_idx, new_idx);
        }

        // 4. Sort bonds by original index
        bonds.sort_unstable_by_key(|(k, _, _, _)| *k);

        // 5. Add bonds using the mapped indices
        for (_old_idx, old_idx1, old_idx2, bond) in bonds {
            // Look up the internal indices corresponding to the original atom indices
            let new_idx1 = idx_map
                .get(&old_idx1)
                // Return error if missing
                .ok_or::<Error>(DataError::MissingAtomIndex(old_idx1 as usize).into())?;
            let new_idx2 = idx_map
                .get(&old_idx2)
                // Return error if missing
                .ok_or::<Error>(DataError::MissingAtomIndex(old_idx2 as usize).into())?;

            // Add the edge using the *internal* indices
            graph.add_edge(*new_idx1, *new_idx2, bond);
            // We ignore the EdgeIndex returned by add_edge here, as petgraph manages it.
        }

        Ok(Self { data: graph })
    }

    pub fn atom_count(&self) -> usize {
        self.data.node_count()
    }

    pub fn atoms<'graph>(&'graph self) -> impl Iterator<Item = &'graph ValenceAtom> + 'graph {
        self.data.node_weights()
    }

    pub fn atom<'graph>(&'graph self, index: AtomIndex) -> Option<&'graph ValenceAtom> {
        self.data.node_weight(index)
    }

    pub fn bond_count(&self) -> usize {
        self.data.edge_count()
    }

    pub fn bonds<'graph>(&'graph self) -> impl Iterator<Item = &'graph ValenceBond> + 'graph {
        self.data.edge_weights()
    }

    pub fn bond<'graph>(&'graph self, index: BondIndex) -> Option<&'graph ValenceBond> {
        self.data.edge_weight(index)
    }

    pub fn bond_atoms<'graph>(
        &'graph self,
        index: BondIndex,
    ) -> Option<(&'graph ValenceAtom, &'graph ValenceAtom)> {
        self.data.edge_endpoints(index).map(|(a, b)| {
            (
                self.data.node_weight(a).unwrap(),
                self.data.node_weight(b).unwrap(),
            )
        })
    }

    pub fn bond_atom_indices<'graph>(
        &'graph self,
        index: BondIndex,
    ) -> Option<(AtomIndex, AtomIndex)> {
        self.data.edge_endpoints(index)
    }

    pub fn atom_bonds<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = &'graph ValenceBond> + 'graph {
        self.data.edges(index).map(|e| e.weight())
    }

    pub fn atom_bond_indices<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = BondIndex> + 'graph {
        self.data.edges(index).map(|e| e.id())
    }

    pub fn atom_neighbors<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = &'graph ValenceAtom> + 'graph {
        self.data
            .neighbors(index)
            .map(|n| self.data.node_weight(n).unwrap())
    }

    pub fn atom_neighbor_indices<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = AtomIndex> + 'graph {
        self.data.neighbors(index)
    }
}

impl fmt::Display for ValenceGraph {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(
            f,
            "Molecule with {} atoms and {} bonds:",
            self.atom_count(),
            self.bond_count()
        )?;
        for (i, atom) in self.atoms().enumerate() {
            writeln!(f, "  Atom {}: {}", i, atom)?;
        }
        for (i, bond) in self.bonds().enumerate() {
            if let Some((a, b)) = self.bond_atoms(BondIndex::new(i)) {
                writeln!(f, "  Bond {}: {} between atoms {} and {}", i, bond, a, b)?;
            }
        }
        Ok(())
    }
}

/// Builder type for ValenceGraphs
pub struct ValenceGraphBuilder {
    graph: StableGraph<ValenceAtom, ValenceBond, Undirected, IndexType>,
    active_atom: Option<AtomIndex>,
    ops: Vec<ValenceGraphOp>,
    validations: ValidationSet,
}

enum ValenceGraphOp {
    AddAtom {
        index: IndexType,
        atom: ValenceAtom,
    },
    AddBond {
        index: IndexType,
        atom1: AtomIndex,
        atom2: AtomIndex,
        bond: ValenceBond,
    },
    ModifyAtom {
        index: AtomIndex,
        property: String,
    },
    ModifyBond {
        index: BondIndex,
        property: String,
    },
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{ValenceAtom, ValenceBond};
    use umol_data::{BondOrder, Element};

    #[test]
    fn test_valence_graph_from_atoms_bonds() {
        let atoms = vec![
            (1, ValenceAtom::new(Element::H)),
            (2, ValenceAtom::new(Element::H)),
        ];
        let bonds = vec![(1, 1, 2, ValenceBond::new(BondOrder::Single))];
        let graph = ValenceGraph::from_atoms_bonds(atoms, bonds);
        assert_eq!(graph.atom_count(), 2);
        assert_eq!(graph.bond_count(), 1);
    }

    #[test]
    #[should_panic]
    fn test_valence_graph_from_atoms_bonds_duplicate_atom_index() {
        let atoms = vec![
            (1, ValenceAtom::new(Element::H)),
            (1, ValenceAtom::new(Element::H)),
        ];
        let bonds = vec![(1, 1, 2, ValenceBond::new(BondOrder::Single))];
        ValenceGraph::from_atoms_bonds(atoms, bonds);
    }

    #[test]
    #[should_panic]
    fn test_valence_graph_from_atoms_bonds_missing_atom_index() {
        let atoms = vec![
            (1, ValenceAtom::new(Element::H)),
            (2, ValenceAtom::new(Element::H)),
        ];
        let bonds = vec![(1, 1, 3, ValenceBond::new(BondOrder::Single))];
        ValenceGraph::from_atoms_bonds(atoms, bonds);
    }

    #[test]
    fn test_valence_graph_try_from_atoms_bonds() {
        let atoms = vec![
            (7, ValenceAtom::new(Element::H)),
            (12, ValenceAtom::new(Element::O)),
            (1, ValenceAtom::new(Element::H)),
        ];
        let bonds = vec![
            (2, 7, 12, ValenceBond::new(BondOrder::Single)),
            (0, 1, 7, ValenceBond::new(BondOrder::Single)),
        ];
        let graph = ValenceGraph::try_from_atoms_bonds(atoms, bonds);
        assert!(graph.is_ok());
        let graph = graph.unwrap();
        assert_eq!(graph.atom_count(), 3);
        assert_eq!(graph.bond_count(), 2);
    }

    #[test]
    fn test_valence_graph_try_from_atoms_bonds_duplicate_atom_index() {
        let atoms = vec![
            (1, ValenceAtom::new(Element::H)),
            (1, ValenceAtom::new(Element::H)),
        ];
        let bonds = vec![(1, 1, 2, ValenceBond::new(BondOrder::Single))];
        let graph = ValenceGraph::try_from_atoms_bonds(atoms, bonds);
        assert!(graph.is_err());
        assert!(matches!(
            graph,
            Err(Error::Data(DataError::DuplicateAtomIndex(1)))
        ));
    }
    #[test]
    fn test_valence_graph_try_from_atoms_bonds_missing_atom_index() {
        let atoms = vec![
            (1, ValenceAtom::new(Element::H)),
            (2, ValenceAtom::new(Element::H)),
        ];
        let bonds = vec![(1, 1, 3, ValenceBond::new(BondOrder::Single))];
        let graph = ValenceGraph::try_from_atoms_bonds(atoms, bonds);
        assert!(graph.is_err());
        assert!(matches!(
            graph,
            Err(Error::Data(DataError::MissingAtomIndex(3)))
        ));
    }

    #[test]
    fn test_valence_graph_atoms() {
        let atoms = vec![
            (1, ValenceAtom::new(Element::H)),
            (2, ValenceAtom::new(Element::H)),
        ];
        let bonds = vec![(1, 1, 2, ValenceBond::new(BondOrder::Single))];
        let graph = ValenceGraph::from_atoms_bonds(atoms, bonds);
        assert_eq!(graph.atom_count(), 2);
        assert_eq!(graph.atom(AtomIndex::new(0)).unwrap().element(), Element::H);
        assert_eq!(graph.atom(AtomIndex::new(1)).unwrap().element(), Element::H);
    }
}
