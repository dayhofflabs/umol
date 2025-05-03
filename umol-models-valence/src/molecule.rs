//! Molecular representation as valence graph

use super::{Atom, Bond};
use petgraph::prelude::*;
use std::collections::{HashMap, HashSet};
use std::{fmt, Display};
use umol::error::DataError;
use umol::{Error, Result};
use umol_data::{BondOrder, BondDonation, ValenceState, Element};

/// The type used for internal atom indices.
pub type AtomIndex = NodeIndex<usize>;
/// The type used for internal bond indices.
pub type BondIndex = EdgeIndex<usize>;

/// Graph model of atoms and bonds, with valence constraints
#[derive(Debug, Clone)]
pub struct Molecule {
    data: StableGraph<Atom, Bond, Undirected, usize>,
}

impl Molecule {
    pub fn atom_count(&self) -> usize {
        self.data.node_count()
    }

    pub fn atoms<'graph>(&'graph self) -> impl Iterator<Item = &'graph Atom> + 'graph {
        self.data.node_weights()
    }

    pub fn atom<'graph>(&'graph self, index: AtomIndex) -> Option<&'graph Atom> {
        self.data.node_weight(index)
    }

    pub fn bond_count(&self) -> usize {
        self.data.edge_count()
    }

    pub fn bonds<'graph>(&'graph self) -> impl Iterator<Item = &'graph Bond> + 'graph {
        self.data.edge_weights()
    }

    pub fn bond<'graph>(&'graph self, index: BondIndex) -> Option<&'graph Bond> {
        self.data.edge_weight(index)
    }

    // TODO: Review naming
    pub fn bond_atoms<'graph>(
        &'graph self,
        index: BondIndex,
    ) -> Option<(&'graph Atom, &'graph Atom)> {
        self.data.edge_endpoints(index).map(|(a, b)| {
            (
                self.data.node_weight(a).unwrap(),
                self.data.node_weight(b).unwrap(),
            )
        })
    }

    // TODO: Review naming
    pub fn bond_atom_indices<'graph>(
        &'graph self,
        index: BondIndex,
    ) -> Option<(AtomIndex, AtomIndex)> {
        self.data.edge_endpoints(index)
    }

    // TODO: Review naming
    pub fn atom_bonds<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = &'graph Bond> + 'graph {
        self.data.edges(index).map(|e| e.weight())
    }

    // TODO: Review naming
    pub fn atom_bond_indices<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = BondIndex> + 'graph {
        self.data.edges(index).map(|e| e.id())
    }

    pub fn atom_neighbors<'graph>(
        &'graph self,
        index: AtomIndex,
    ) -> impl Iterator<Item = &'graph Atom> + 'graph {
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
    // TODO: Add methods for converting some atoms/bonds to builders
}

impl Display for Molecule {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(
            f,
            "Molecule with {} atoms and {} bonds:",
            self.atom_count(),
            self.bond_count()
        )?;
        for (i, atom) in self.atoms().enumerate() {
            writeln!(f, "  Atom {}: {:?}", i, atom)?;
        }
        for (i, bond) in self.bonds().enumerate() {
            if let Some((a, b)) = self.bond_atoms(BondIndex::new(i)) {
                writeln!(f, "  Bond {}: {} between atoms {:?} and {:?}", i, bond, a, b)?;
            }
        }
        Ok(())
    }
}

/// Builder type for ValenceGraphs, allowing incremental construction and validation.
pub struct MoleculeBuilder {
    atom_builders: Vec<(usize, AtomBuilder)>,
    bond_builders: Vec<(usize, usize, usize, BondBuilder)>,
}

impl MoleculeBuilder {
    pub fn new() -> Self {
        Self {
            atom_builders: Vec::new(),
            bond_builders: Vec::new(),
        }
    }

    pub fn with_capacity(atom_capacity: usize, bond_capacity: usize) -> Self {
        Self {
            atom_builders: Vec::with_capacity(atom_capacity),
            bond_builders: Vec::with_capacity(bond_capacity),
        }
    }

    pub fn add_atom(
        &mut self,
        atom: Into<Atom>,
    ) -> Result<usize> {
        let builder = AtomBuilder::from(atom);
        let idx = self.atom_builders.len();
        self.atom_builders.push((idx, builder));
        Ok(idx)
    }

    pub fn add_atoms(&mut self, atoms: impl IntoIterator<Item = Atom>) -> Result<Vec<usize>> {
        let builders = atoms.into_iter().map(|atom| AtomBuilder::from(atom));
        let atom_count = self.atom_builders.len();
        let mut indices = Vec::with_capacity(atom_count);   
        for (idx, builder) in (atom_count..).zip(builders) {
            self.atom_builders.push((idx, builder));
            indices.push(idx);
        }
        Ok(indices)
    }

    pub fn add_bond(
        &mut self,
        idx1: usize,
        idx2: usize,
        bond: Into<Bond>,
    ) -> Result<usize> {
        if idx1 == idx2 {
            return Err(DataError::LoopBond(idx1).into());
        }
        if idx1 >= self.atom_builders.len() || idx2 >= self.atom_builders.len() {
            return Err(DataError::MissingAtomIndex(idx1).into());
        }
        let builder = BondBuilder::from(bond);
        let idx = self.bond_builders.len();
        self.bond_builders.push((idx, idx1, idx2, builder));
        Ok(idx)
    }

    pub fn add_bonds(&mut self, bonds: impl IntoIterator<Item = (usize, usize, Bond)>) -> Result<Vec<usize>> {
        let builders = bonds.into_iter().map(|(idx1, idx2, bond)| (idx1, idx2, BondBuilder::from(bond)));
        let bond_count = self.bond_builders.len();
        let mut indices = Vec::with_capacity(bond_count);
        for (idx, (idx1, idx2, builder)) in (bond_count..).zip(builders) {
            if idx1 >= self.atom_builders.len() || idx2 >= self.atom_builders.len() {
                return Err(DataError::MissingAtomIndex(idx1).into());
            }   
            self.bond_builders.push((idx, idx1, idx2, builder));
            indices.push(idx);
        }
        Ok(indices)
    }

    // TODO: Review naming
    fn atom_builder_mut(&mut self, idx: usize) -> Result<&mut ValenceAtomBuilder> {
        self.atom_builders.get_mut(&idx)
            .ok_or(DataError::MissingAtomIndex(idx).into())
    }

    // TODO: Review naming
    fn bond_builder_mut(&mut self, idx: usize) -> Result<&mut ValenceBondBuilder> {
        self.bond_builders.get_mut(&idx)
            .ok_or(DataError::MissingBondIndex(idx).into())
    }

    pub fn set_atom_element(&mut self, idx: usize, element: Element) -> Result<&mut ValenceAtomBuilder> {
        let builder = self.atom_builder_mut(idx)?;
        builder.set_element(element);
        Ok(builder)
    }

    pub fn set_atom_charge(&mut self, idx: usize, charge: i8) -> Result<&mut ValenceAtomBuilder> {
        let builder = self.atom_builder_mut(idx)?;
        builder.set_charge(charge);
        Ok(builder)
    }

    pub fn set_atom_lone_pairs(&mut self, idx: usize, lp: u8) -> Result<&mut ValenceAtomBuilder> {
         let builder = self.atom_builder_mut(idx)?;
         builder.set_lone_pairs(lp);
         Ok(builder)
    }

    pub fn set_atom_unpaired_electrons(&mut self, idx: usize, unpaired: u8) -> Result<&mut ValenceAtomBuilder> {
         let builder = self.atom_builder_mut(idx)?;
         builder.set_unpaired_electrons(unpaired);
         Ok(builder)
    }

    pub fn set_atom_multiplicity(&mut self, idx: usize, mult: u8) -> Result<&mut ValenceAtomBuilder> {
         let builder = self.atom_builder_mut(idx)?;
         builder.set_multiplicity(mult);
         Ok(builder)
    }

    pub fn set_atom_implicit_hydrogens(&mut self, idx: usize, count: u8) -> Result<&mut ValenceAtomBuilder> {
         let builder = self.atom_builder_mut(idx)?;
         builder.set_implicit_hydrogens(count);
         Ok(builder)
    }

    pub fn set_atom_bond_sum(&mut self, idx: usize, sum: u8) -> Result<&mut ValenceAtomBuilder> {
         let builder = self.atom_builder_mut(idx)?;
         builder.with_bond_sum(sum);
         Ok(builder)
    }

    pub fn set_bond_order(&mut self, idx: usize, order: BondOrder) -> Result<&mut ValenceBondBuilder> {
        let builder = self.bond_builder_mut(idx)?;
        builder.with_order(order);
        Ok(builder)
    }

    pub fn set_bond_donation(&mut self, idx: usize, donation: BondDonation) -> Result<&mut ValenceBondBuilder> {
        let builder = self.bond_builder_mut(idx)?;
        builder.with_donation(donation);
        Ok(builder)
    }

    pub fn build(self) -> Result<Molecule> {
        let mut atom_builders = self.atom_builders;

        // Update bond sums and lone pairs
        for (idx, idx1, idx2, ref bond_builder) in &self.bond_builders {
            let bond_order = bond_builder.order().ok_or_else(
                || Error::Data(DataError::MissingBondProperty(*idx, "order".to_string())))?;
            let bond_donation = bond_builder.donation().unwrap_or(BondDonation::Shared).value();

            let builder1 = atom_builders.get_mut(idx1)
                .ok_or_else(|| Error::Data(DataError::MissingAtomIndex(*idx1)))?;
            builder1.update_bond_sum(|sum| sum + bond_order);
            if bond_donation != 0 {
                builder1.update_lone_pairs(|lp| lp + bond_donation * bond_order);
            }

            let builder2 = atom_builders.get_mut(idx2)
                .ok_or_else(|| Error::Data(DataError::MissingAtomIndex(*idx2)))?;
            builder2.update_bond_sum(|sum| sum + bond_order);
            if bond_donation != 0 {
                builder2.update_lone_pairs(|lp| lp + bond_donation * bond_order);
            }
        }

        // Validate and finalize atoms
        let mut atoms: Vec<(usize, Atom)> = Vec::with_capacity(atom_builders.len());
        for (idx, builder) in atom_builders {
            let element = builder.element().ok_or_else(
                || Error::Data(DataError::MissingAtomProperty(idx, "element".to_string())))?;
            let charge = builder.charge();
            let lone_pairs = builder.lone_pairs();
            let unpaired_electrons = builder.unpaired_electrons();
            let multiplicity = builder.multiplicity();
            let explicit_implicit_h = builder.explicit_implicit_h();
            let bond_sum = builder.bond_sum();
            
            let candidate_states = infer_types(
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
                
                let final_atom = Atom {
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
                return Err(Error::Data(DataError::InvalidValenceAtom(
                    format!("Invalid valence atom {}: {:?}", orig_idx, valid_assignments)
                ).into());
            }
        }

        // --- Step 6: Finalize Bonds & Construct Graph --- 
        let mut final_bonds: Vec<(usize, usize, usize, Bond)> = Vec::with_capacity(self.bond_builders.len());
        for (orig_bond_idx, orig_idx1, orig_idx2, bond_builder) in self.bond_builders {
            let final_bond = bond_builder.build()?;
            final_bonds.push((orig_bond_idx, orig_idx1, orig_idx2, final_bond));
        }
        
        todo!()
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
