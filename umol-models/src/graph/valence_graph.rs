//! Valence graph model
//!
//! Graph model of atoms and bonds, with valence constraints

use crate::graph::{ValenceAtom, ValenceBond, ValenceAtomBuilder, ValenceBondBuilder};
use petgraph::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fmt;
use umol::error::DataError;
use umol::{Error, Result};
use umol_data::{BondOrder, BondDonation, ValenceState, Element};
// use crate::graph::find_matching_states;

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

// impl ValenceGraph {
//     /// Creates a ValenceGraph from vectors of atoms and bonds with arbitrary original indices.
//     ///
//     /// Atom and bond lists are processed to map potentially non-sequential, 1-based,
//     /// or otherwise arbitrary indices (`I`) to dense, 0-based internal indices used
//     /// by the underlying graph structure (`AtomIndex`, `BondIndex`).
//     ///
//     /// # Arguments
//     /// * \`atoms\`: A vector of tuples \`(original_index, atom_data)\`.
//     /// * \`bonds\`: A vector of tuples \`(original_bond_index, original_atom_index1, original_atom_index2, bond_data)\`.
//     ///
//     /// # Type Parameters
//     /// * \`I\`: The type used for original indices in the input vectors. Must be sortable, hashable, and copyable.
//     ///
//     /// # Returns
//     /// A \`Result\` containing the new \`ValenceGraph\` on success, or a \`GraphBuildError\`
//     /// if issues like duplicate or missing atom indices are found.
//     ///
//     /// # Behavior
//     /// 1. Checks for duplicate original atom indices.
//     /// 2. Sorts atoms by their original index.
//     /// 3. Adds atoms to the graph, creating a mapping from original index (\`I\`) to internal \`AtomIndex\`.
//     /// 4. Sorts bonds by their original index.
//     /// 5. Adds bonds to the graph using the mapped internal \`AtomIndex\` values.
//     pub fn from_atoms_bonds(
//         mut atoms: Vec<(IndexType, ValenceAtom)>,
//         mut bonds: Vec<(IndexType, IndexType, IndexType, ValenceBond)>,
//     ) -> Self {
//         let mut graph = StableGraph::with_capacity(atoms.len(), bonds.len());
//         let mut idx_map = HashMap::with_capacity(atoms.len());
//         let mut seen_indices = HashSet::with_capacity(atoms.len());

//         // 1. Check for duplicate original atom indices before sorting
//         for (old_idx, _) in atoms.iter() {
//             if !seen_indices.insert(*old_idx) {
//                 panic!("Duplicate atom index: {}", old_idx);
//             }
//         }

//         // 2. Sort atoms by original index
//         atoms.sort_unstable_by_key(|(k, _)| *k);

//         // 3. Add atoms and build the index map
//         for (old_idx, atom) in atoms {
//             let new_idx = graph.add_node(atom);
//             idx_map.insert(old_idx, new_idx);
//         }

//         // 4. Sort bonds by original index
//         bonds.sort_unstable_by_key(|(k, _, _, _)| *k);

//         // 5. Add bonds using the mapped indices
//         for (_old_idx, old_idx1, old_idx2, bond) in bonds {
//             // Look up the internal indices corresponding to the original atom indices
//             let new_idx1 = idx_map
//                 .get(&old_idx1)
//                 .unwrap_or_else(|| panic!("Missing atom index: {}", old_idx1));
//             let new_idx2 = idx_map
//                 .get(&old_idx2)
//                 .unwrap_or_else(|| panic!("Missing atom index: {}", old_idx2));

//             // Add the edge using the *internal* indices
//             graph.add_edge(*new_idx1, *new_idx2, bond);
//             // We ignore the EdgeIndex returned by add_edge here, as petgraph manages it.
//         }

//         Self { data: graph }
//     }

//     /// Tries to create a ValenceGraph from vectors of atoms and bonds with arbitrary original indices.
//     ///
//     /// Atom and bond lists are processed to map potentially non-sequential, 1-based,
//     /// or otherwise arbitrary indices (\`IndexType\`) to dense, 0-based internal indices used
//     /// by the underlying graph structure (\`AtomIndex\`, \`BondIndex\`).
//     ///
//     /// # Arguments
//     /// * \`atoms\`: A vector of tuples \`(original_index, atom_data)\`.
//     /// * \`bonds\`: A vector of tuples \`(original_bond_index, original_atom_index1, original_atom_index2, bond_data)\`.
//     ///
//     /// # Returns
//     /// A \`Result\` containing the new \`ValenceGraph\` on success, or a \`umol::Error\`
//     /// wrapping a \`DataError\` if issues like duplicate or missing atom indices are found.
//     ///
//     /// # Behavior
//     /// 1. Checks for duplicate original atom indices.
//     /// 2. Sorts atoms by their original index.
//     /// 3. Adds atoms to the graph, creating a mapping from original index (\`IndexType\`) to internal \`AtomIndex\`.
//     /// 4. Sorts bonds by their original index.
//     /// 5. Adds bonds to the graph using the mapped internal \`AtomIndex\` values.
//     pub fn try_from_atoms_bonds(
//         mut atoms: Vec<(IndexType, ValenceAtom)>,
//         mut bonds: Vec<(IndexType, IndexType, IndexType, ValenceBond)>,
//     ) -> Result<Self> {
//         let mut graph = StableGraph::with_capacity(atoms.len(), bonds.len());
//         let mut idx_map = HashMap::with_capacity(atoms.len());
//         let mut seen_indices = HashSet::with_capacity(atoms.len());

//         // 1. Check for duplicate original atom indices before sorting
//         for (old_idx, _) in atoms.iter() {
//             if !seen_indices.insert(*old_idx) {
//                 // Return error on duplicate
//                 return Err(DataError::DuplicateAtomIndex(*old_idx as usize).into());
//             }
//         }

//         // 2. Sort atoms by original index
//         atoms.sort_unstable_by_key(|(k, _)| *k);

//         // 3. Add atoms and build the index map
//         for (old_idx, atom) in atoms {
//             let new_idx = graph.add_node(atom);
//             idx_map.insert(old_idx, new_idx);
//         }

//         // 4. Sort bonds by original index
//         bonds.sort_unstable_by_key(|(k, _, _, _)| *k);

//         // 5. Add bonds using the mapped indices
//         for (_old_idx, old_idx1, old_idx2, bond) in bonds {
//             // Look up the internal indices corresponding to the original atom indices
//             let new_idx1 = idx_map
//                 .get(&old_idx1)
//                 // Return error if missing
//                 .ok_or::<Error>(DataError::MissingAtomIndex(old_idx1 as usize).into())?;
//             let new_idx2 = idx_map
//                 .get(&old_idx2)
//                 // Return error if missing
//                 .ok_or::<Error>(DataError::MissingAtomIndex(old_idx2 as usize).into())?;

//             // Add the edge using the *internal* indices
//             graph.add_edge(*new_idx1, *new_idx2, bond);
//             // We ignore the EdgeIndex returned by add_edge here, as petgraph manages it.
//         }

//         Ok(Self { data: graph })
//     }

//     pub fn atom_count(&self) -> usize {
//         self.data.node_count()
//     }

//     pub fn atoms<'graph>(&'graph self) -> impl Iterator<Item = &'graph ValenceAtom> + 'graph {
//         self.data.node_weights()
//     }

//     pub fn atom<'graph>(&'graph self, index: AtomIndex) -> Option<&'graph ValenceAtom> {
//         self.data.node_weight(index)
//     }

//     pub fn bond_count(&self) -> usize {
//         self.data.edge_count()
//     }

//     pub fn bonds<'graph>(&'graph self) -> impl Iterator<Item = &'graph ValenceBond> + 'graph {
//         self.data.edge_weights()
//     }

//     pub fn bond<'graph>(&'graph self, index: BondIndex) -> Option<&'graph ValenceBond> {
//         self.data.edge_weight(index)
//     }

//     pub fn bond_atoms<'graph>(
//         &'graph self,
//         index: BondIndex,
//     ) -> Option<(&'graph ValenceAtom, &'graph ValenceAtom)> {
//         self.data.edge_endpoints(index).map(|(a, b)| {
//             (
//                 self.data.node_weight(a).unwrap(),
//                 self.data.node_weight(b).unwrap(),
//             )
//         })
//     }

//     pub fn bond_atom_indices<'graph>(
//         &'graph self,
//         index: BondIndex,
//     ) -> Option<(AtomIndex, AtomIndex)> {
//         self.data.edge_endpoints(index)
//     }

//     pub fn atom_bonds<'graph>(
//         &'graph self,
//         index: AtomIndex,
//     ) -> impl Iterator<Item = &'graph ValenceBond> + 'graph {
//         self.data.edges(index).map(|e| e.weight())
//     }

//     pub fn atom_bond_indices<'graph>(
//         &'graph self,
//         index: AtomIndex,
//     ) -> impl Iterator<Item = BondIndex> + 'graph {
//         self.data.edges(index).map(|e| e.id())
//     }

//     pub fn atom_neighbors<'graph>(
//         &'graph self,
//         index: AtomIndex,
//     ) -> impl Iterator<Item = &'graph ValenceAtom> + 'graph {
//         self.data
//             .neighbors(index)
//             .map(|n| self.data.node_weight(n).unwrap())
//     }

//     pub fn atom_neighbor_indices<'graph>(
//         &'graph self,
//         index: AtomIndex,
//     ) -> impl Iterator<Item = AtomIndex> + 'graph {
//         self.data.neighbors(index)
//     }
// }

// impl fmt::Display for ValenceGraph {
//     fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//         writeln!(
//             f,
//             "Molecule with {} atoms and {} bonds:",
//             self.atom_count(),
//             self.bond_count()
//         )?;
//         for (i, atom) in self.atoms().enumerate() {
//             writeln!(f, "  Atom {}: {:?}", i, atom)?;
//         }
//         for (i, bond) in self.bonds().enumerate() {
//             if let Some((a, b)) = self.bond_atoms(BondIndex::new(i)) {
//                 writeln!(f, "  Bond {}: {} between atoms {:?} and {:?}", i, bond, a, b)?;
//             }
//         }
//         Ok(())
//     }
// }

/// Builder type for ValenceGraphs, allowing incremental construction and validation.
pub struct ValenceGraphBuilder {
    atom_builders: HashMap<IndexType, ValenceAtomBuilder>,
    bond_builders: Vec<(IndexType, IndexType, IndexType, ValenceBondBuilder)>,
}

impl ValenceGraphBuilder {
    /// Creates a new, empty ValenceGraphBuilder.
    pub fn new() -> Self {
        Self {
            atom_builders: HashMap::new(),
            bond_builders: Vec::new(),
        }
    }

    /// Gets a mutable reference to the ValenceAtomBuilder for the given original index,
    /// inserting a default builder if one doesn't exist.
    fn atom_builder_mut(&mut self, idx: IndexType) -> Result<&mut ValenceAtomBuilder> {
        self.atom_builders.get_mut(&idx)
            .ok_or(DataError::MissingAtomIndex(idx as usize).into())
    }

    /// Sets the element for the atom with the given original index.
    /// Returns an error if the element was already set to a different value.
    pub fn set_atom_element(&mut self, idx: IndexType, element: Element) -> Result<()> {
        let builder = self.atom_builder_mut(idx)?;
        builder.element(element);
        Ok(())
    }

    /// Sets the formal charge for the atom with the given original index.
    /// Returns an error if the charge was already set to a different value.
    pub fn set_atom_charge(&mut self, orig_idx: IndexType, charge: i8) -> Result<()> {
        let builder = self.atom_builder_mut(orig_idx);
        match builder.charge {
            Some(existing) if existing != charge => {
                 Err(DataError::AtomPropertyAlreadySet(
                    orig_idx as usize,
                    "charge".to_string(),
                ).into())
            }
             _ => {
                builder.charge = Some(charge);
                 Ok(())
             }
        }
    }

    /// Sets the number of lone pairs for the atom with the given original index.
    /// Returns an error if the lone pairs were already set to a different value.
    pub fn set_atom_lone_pairs(&mut self, orig_idx: IndexType, lp: u8) -> Result<()> {
         let builder = self.atom_builder_mut(orig_idx);
         match builder.lone_pairs {
             Some(existing) if existing != lp => {
                 Err(DataError::AtomPropertyAlreadySet(
                     orig_idx as usize,
                     "lone_pairs".to_string(),
                 ).into())
             }
             _ => {
                 builder.lone_pairs = Some(lp);
                 Ok(())
             }
         }
    }

    /// Sets the number of unpaired electrons for the atom with the given original index.
    /// Returns an error if the unpaired electrons were already set to a different value.
    pub fn set_atom_unpaired_electrons(&mut self, orig_idx: IndexType, unpaired: u8) -> Result<()> {
         let builder = self.atom_builder_mut(orig_idx);
         match builder.unpaired_electrons {
             Some(existing) if existing != unpaired => {
                 Err(DataError::AtomPropertyAlreadySet(
                     orig_idx as usize,
                     "unpaired_electrons".to_string(),
                 ).into())
             }
             _ => {
                 builder.unpaired_electrons = Some(unpaired);
                 Ok(())
             }
         }
    }

    /// Sets the multiplicity for the atom with the given original index.
    /// Returns an error if the multiplicity was already set to a different value.
    pub fn set_atom_multiplicity(&mut self, orig_idx: IndexType, mult: u8) -> Result<()> {
        let builder = self.atom_builder_mut(orig_idx);
        match builder.multiplicity {
            Some(existing) if existing != mult => {
                Err(DataError::AtomPropertyAlreadySet(
                    orig_idx as usize,
                    "multiplicity".to_string(),
                ).into())
            }
            _ => {
                builder.multiplicity = Some(mult);
                Ok(())
            }
        }
    }

    /// Sets the explicit number of implicit hydrogens for the atom with the given original index.
    /// Returns an error if the explicit count was already set to a different value.
    pub fn set_atom_explicit_implicit_h(&mut self, orig_idx: IndexType, count: u8) -> Result<()> {
        let builder = self.atom_builder_mut(orig_idx);
        match builder.explicit_implicit_h {
            Some(existing) if existing != count => {
                 Err(DataError::AtomPropertyAlreadySet(
                    orig_idx as usize,
                    "explicit_implicit_h".to_string(),
                 ).into())
            }
            _ => {
                builder.explicit_implicit_h = Some(count);
                 Ok(())
             }
        }
    }

    /// Adds a bond request to the builder.
    /// Performs basic validation like checking for self-loops.
    pub fn add_bond(
        &mut self,
        orig_bond_idx: IndexType,
        orig_idx1: IndexType,
        orig_idx2: IndexType,
        order: BondOrder,
        donation: BondDonation,
    ) -> Result<()> {
        // Basic validation: prevent self-loops
        if orig_idx1 == orig_idx2 {
            return Err(DataError::SelfLoopBond(orig_idx1 as usize).into());
        }

        // Create the bond builder
        let bond_builder = ValenceBondBuilder::new(order).donation(donation);

        // Store the request
        self.bond_builders.push((orig_bond_idx, orig_idx1, orig_idx2, bond_builder));

        Ok(())
    }

    /// Builds the ValenceGraph after all atoms and bonds have been added.
    /// This performs the final validation against valence states and connectivity.
    pub fn build(self) -> Result<ValenceGraph> {
        // 1. Create mutable copy to update bond sums
        let mut atom_builders = self.atom_builders;

        // 2. Iterate through bond requests to calculate bond sums
        for (_orig_bond_idx, orig_idx1, orig_idx2, ref bond_builder) in &self.bond_builders {
            let bond_order_value = bond_builder.order().value();

            // Update first atom's bond sum
            let builder1 = atom_builders.get_mut(orig_idx1)
                .ok_or_else(|| Error::Data(DataError::MissingAtomIndex(*orig_idx1 as usize)))?;
            builder1.increment_bond_sum(bond_order_value);

            // Update second atom's bond sum
            let builder2 = atom_builders.get_mut(orig_idx2)
                .ok_or_else(|| Error::Data(DataError::MissingAtomIndex(*orig_idx2 as usize)))?;
            builder2.increment_bond_sum(bond_order_value);
        }

        // --- Step 5: Atom Validation & Finalization --- 
        let mut final_atoms: Vec<(IndexType, ValenceAtom)> = Vec::with_capacity(atom_builders.len());
        for (orig_idx, builder) in atom_builders {
            let element = builder.element().ok_or_else(|| Error::Data(DataError::MissingAtomProperty(orig_idx, "element".to_string())))?;
            let charge = builder.charge();
            let lone_pairs = builder.lone_pairs();
            let unpaired_electrons = builder.unpaired_electrons();
            let multiplicity = builder.multiplicity();
            let explicit_implicit_h = builder.explicit_implicit_h();
            let bond_sum = builder.bond_sum();
            
            let candidate_states = find_matching_states(
                element, 
                charge, 
                lone_pairs, 
                unpaired_electrons, 
                multiplicity
            )?;

            let mut valid_assignments = Vec::new();

            for state in candidate_states {
                let required_valence = state.valence();

                if required_valence < bond_sum {
                    continue;
                }

                let calculated_implicit_h = required_valence - bond_sum;
                
                if let Some(explicit_h) = explicit_implicit_h {
                    if explicit_h != calculated_implicit_h {
                        continue;
                    }
                }
                
                valid_assignments.push((state, calculated_implicit_h));
            }

            if valid_assignments.len() == 1 {
                let (final_state, final_implicit_h) = valid_assignments[0];
                
                let final_atom = ValenceAtom {
                    element: final_state.element(),
                    charge: final_state.charge(),
                    lone_pairs: final_state.lone_pairs(),
                    unpaired_electrons: final_state.unpaired_electrons(),
                    multiplicity: final_state.multiplicity(),
                    implicit_hydrogens: final_implicit_h,
                    bond_sum: bond_sum,
                };
                final_atoms.push((orig_idx, final_atom));
            } else {
                return Err(Error::Data(DataError::ValenceCheckFailed(orig_idx, valid_assignments.len())));
            }
        }

        // --- Step 6: Finalize Bonds & Construct Graph --- 
        let mut final_bonds: Vec<(IndexType, IndexType, IndexType, ValenceBond)> = Vec::with_capacity(self.bond_builders.len());
        for (orig_bond_idx, orig_idx1, orig_idx2, bond_builder) in self.bond_builders {
            let final_bond = bond_builder.build()?;
            final_bonds.push((orig_bond_idx, orig_idx1, orig_idx2, final_bond));
        }
        
        ValenceGraph::try_from_atoms_bonds(final_atoms, final_bonds)
    }
}

// Temporary placeholder for the assumed function
// This would normally live in umol-data::valence_state
// or be imported if umol-models depends on umol-data directly.
pub fn find_matching_states(
    element: Element,
    charge: Option<i8>,
    lone_pairs: Option<u8>,
    unpaired_electrons: Option<u8>,
    multiplicity: Option<u8>,
) -> Result<Vec<ValenceState>> {
    eprintln!(
        "Warning: Using placeholder find_matching_states({:?}, {:?}, {:?}, {:?}, {:?})",
        element, charge, lone_pairs, unpaired_electrons, multiplicity
    );
    if element == umol_data::Element::C {
        Ok(vec![ValenceState::new(element, 0, 0, 0, 1, 4)])
    } else {
        Ok(vec![])
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::{ValenceAtom, ValenceBond};
//     use umol_data::{BondOrder, Element};

//     #[test]
//     fn test_valence_graph_from_atoms_bonds() {
//         let atoms = vec![
//             (1, ValenceAtom::new(Element::H)),
//             (2, ValenceAtom::new(Element::H)),
//         ];
//         let bonds = vec![(1, 1, 2, ValenceBond::new(BondOrder::Single))];
//         let graph = ValenceGraph::from_atoms_bonds(atoms, bonds);
//         assert_eq!(graph.atom_count(), 2);
//         assert_eq!(graph.bond_count(), 1);
//     }

//     #[test]
//     #[should_panic]
//     fn test_valence_graph_from_atoms_bonds_duplicate_atom_index() {
//         let atoms = vec![
//             (1, ValenceAtom::new(Element::H)),
//             (1, ValenceAtom::new(Element::H)),
//         ];
//         let bonds = vec![(1, 1, 2, ValenceBond::new(BondOrder::Single))];
//         ValenceGraph::from_atoms_bonds(atoms, bonds);
//     }

//     #[test]
//     #[should_panic]
//     fn test_valence_graph_from_atoms_bonds_missing_atom_index() {
//         let atoms = vec![
//             (1, ValenceAtom::new(Element::H)),
//             (2, ValenceAtom::new(Element::H)),
//         ];
//         let bonds = vec![(1, 1, 3, ValenceBond::new(BondOrder::Single))];
//         ValenceGraph::from_atoms_bonds(atoms, bonds);
//     }

//     #[test]
//     fn test_valence_graph_try_from_atoms_bonds() {
//         let atoms = vec![
//             (7, ValenceAtom::new(Element::H)),
//             (12, ValenceAtom::new(Element::O)),
//             (1, ValenceAtom::new(Element::H)),
//         ];
//         let bonds = vec![
//             (2, 7, 12, ValenceBond::new(BondOrder::Single)),
//             (0, 1, 7, ValenceBond::new(BondOrder::Single)),
//         ];
//         let graph = ValenceGraph::try_from_atoms_bonds(atoms, bonds);
//         assert!(graph.is_ok());
//         let graph = graph.unwrap();
//         assert_eq!(graph.atom_count(), 3);
//         assert_eq!(graph.bond_count(), 2);
//     }

//     #[test]
//     fn test_valence_graph_try_from_atoms_bonds_duplicate_atom_index() {
//         let atoms = vec![
//             (1, ValenceAtom::new(Element::H)),
//             (1, ValenceAtom::new(Element::H)),
//         ];
//         let bonds = vec![(1, 1, 2, ValenceBond::new(BondOrder::Single))];
//         let graph = ValenceGraph::try_from_atoms_bonds(atoms, bonds);
//         assert!(graph.is_err());
//         assert!(matches!(
//             graph,
//             Err(Error::Data(DataError::DuplicateAtomIndex(1)))
//         ));
//     }
//     #[test]
//     fn test_valence_graph_try_from_atoms_bonds_missing_atom_index() {
//         let atoms = vec![
//             (1, ValenceAtom::new(Element::H)),
//             (2, ValenceAtom::new(Element::H)),
//         ];
//         let bonds = vec![(1, 1, 3, ValenceBond::new(BondOrder::Single))];
//         let graph = ValenceGraph::try_from_atoms_bonds(atoms, bonds);
//         assert!(graph.is_err());
//         assert!(matches!(
//             graph,
//             Err(Error::Data(DataError::MissingAtomIndex(3)))
//         ));
//     }

//     #[test]
//     fn test_valence_graph_atoms() {
//         let atoms = vec![
//             (1, ValenceAtom::new(Element::H)),
//             (2, ValenceAtom::new(Element::H)),
//         ];
//         let bonds = vec![(1, 1, 2, ValenceBond::new(BondOrder::Single))];
//         let graph = ValenceGraph::from_atoms_bonds(atoms, bonds);
//         assert_eq!(graph.atom_count(), 2);
//         assert_eq!(graph.atom(AtomIndex::new(0)).unwrap().element(), Element::H);
//         assert_eq!(graph.atom(AtomIndex::new(1)).unwrap().element(), Element::H);
//     }
// }
