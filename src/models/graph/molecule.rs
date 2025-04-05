// Graph molecule implementation

use std::collections::HashMap;
use std::fmt;

use crate::core::types::{AtomIndex, BondIndex};
use crate::core::Model;
use super::atom::Atom;
use super::bond::Bond;
use crate::error::Error;

/// A molecule represented as a graph of atoms and bonds
#[derive(Debug, Clone)]
pub struct Molecule {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    atom_bonds: HashMap<AtomIndex, Vec<BondIndex>>,
    bond_atoms: HashMap<BondIndex, (AtomIndex, AtomIndex)>,
}

impl Model for Molecule {
    type Atom = Atom;
    type Bond = Bond;
    type Molecule = Self;

    fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    fn num_bonds(&self) -> usize {
        self.bonds.len()
    }

    fn atom(&self, index: AtomIndex) -> Option<&Atom> {
        self.atoms.get(index.0)
    }

    fn bond(&self, index: BondIndex) -> Option<&Bond> {
        self.bonds.get(index.0)
    }

    fn bond_atoms(&self, index: BondIndex) -> Option<(AtomIndex, AtomIndex)> {
        self.bond_atoms.get(&index).copied()
    }

    fn atom_bonds(&self, index: AtomIndex) -> Vec<BondIndex> {
        self.atom_bonds.get(&index).cloned().unwrap_or_default()
    }

    fn atom_neighbors(&self, index: AtomIndex) -> Vec<AtomIndex> {
        self.atom_bonds(index)
            .iter()
            .filter_map(|&bond_idx| {
                self.bond_atoms(bond_idx).map(|(a, b)| {
                    if a == index { b } else { a }
                })
            })
            .collect()
    }
}

impl Molecule {
    /// Create a new empty molecule
    pub fn new() -> Self {
        Self {
            atoms: Vec::new(),
            bonds: Vec::new(),
            atom_bonds: HashMap::new(),
            bond_atoms: HashMap::new(),
        }
    }

    /// Add an atom to the molecule
    pub fn add_atom(&mut self, atom: Atom) -> AtomIndex {
        let index = AtomIndex(self.atoms.len());
        self.atoms.push(atom);
        self.atom_bonds.insert(index, Vec::new());
        index
    }

    /// Add a bond between two atoms
    pub fn add_bond(&mut self, a: AtomIndex, b: AtomIndex, bond: Bond) -> BondIndex {
        let index = BondIndex(self.bonds.len());
        self.bonds.push(bond);
        self.bond_atoms.insert(index, (a, b));
        self.atom_bonds.entry(a).or_default().push(index);
        self.atom_bonds.entry(b).or_default().push(index);
        index
    }

    /// Get all atoms in the molecule
    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }

    /// Get all bonds in the molecule
    pub fn bonds(&self) -> &[Bond] {
        &self.bonds
    }
}

impl fmt::Display for Molecule {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "Molecule with {} atoms and {} bonds:", self.num_atoms(), self.num_bonds())?;
        for (i, atom) in self.atoms.iter().enumerate() {
            writeln!(f, "  Atom {}: {}", i, atom)?;
        }
        for (i, bond) in self.bonds.iter().enumerate() {
            if let Some((a, b)) = self.bond_atoms.get(&BondIndex(i)) {
                writeln!(f, "  Bond {}: {} between atoms {} and {}", i, bond, a.0, b.0)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Element;

    #[test]
    fn test_molecule_new() {
        let mol = Molecule::new();
        assert_eq!(mol.num_atoms(), 0);
        assert_eq!(mol.num_bonds(), 0);
    }

    #[test]
    fn test_molecule_add_atom() {
        let mut mol = Molecule::new();
        let atom = Atom::new(Element::C);
        let index = mol.add_atom(atom);
        assert_eq!(mol.num_atoms(), 1);
        assert_eq!(mol.atom(index), Some(&Atom::new(Element::C)));
    }

    #[test]
    fn test_molecule_add_bond() {
        let mut mol = Molecule::new();
        let c1 = mol.add_atom(Atom::new(Element::C));
        let c2 = mol.add_atom(Atom::new(Element::C));
        let bond = Bond::single();
        let bond_index = mol.add_bond(c1, c2, bond);
        
        assert_eq!(mol.num_bonds(), 1);
        assert_eq!(mol.bond(bond_index), Some(&Bond::single()));
        assert_eq!(mol.bond_atoms(bond_index), Some((c1, c2)));
        assert_eq!(mol.atom_bonds(c1), vec![bond_index]);
        assert_eq!(mol.atom_bonds(c2), vec![bond_index]);
        assert_eq!(mol.atom_neighbors(c1), vec![c2]);
        assert_eq!(mol.atom_neighbors(c2), vec![c1]);
    }

    #[test]
    fn test_molecule_invalid_bond() {
        let mut mol = Molecule::new();
        let c1 = mol.add_atom(Atom::new(Element::C));
        let invalid_index = AtomIndex(1);
        let bond = Bond::single();
        
        assert!(mol.add_bond(c1, invalid_index, bond).is_err());
    }
}
