//! MoleculeBuilder: IR molecule builder for parsers.

use crate::io::ir::{AtomSymbol, SourceFormat};
use std::mem;

use super::{Atom, Bond, BondDir, BondOrder, BondSymbol, Chirality, Molecule};
use umol_data::Element;

pub struct AtomData {
    pub element: Element,
    pub isotope: Option<u32>,
    pub charge: Option<i32>,
    pub hydrogen_count: Option<u32>,
    pub class: Option<u32>,
    pub aromatic: bool,
    pub implicit_h: bool,
    pub chirality: Option<Chirality>,
    pub unknown_symbol: bool,
}

pub struct BondData {
    pub order: BondOrder,
    pub dir: Option<BondDir>,
}

pub struct MoleculeBuilder {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    molecules: Vec<Molecule>,
}

impl MoleculeBuilder {
    pub fn with_capacity(approx_atoms: usize, approx_bonds: usize) -> Self {
        Self {
            atoms: Vec::with_capacity(approx_atoms),
            bonds: Vec::with_capacity(approx_bonds),
            molecules: Vec::new(),
        }
    }

    pub fn clear_reuse(&mut self) {
        self.atoms.clear();
        self.bonds.clear();
        self.molecules.clear();
    }

    #[inline]
    pub fn on_atom(&mut self, a: AtomData) -> u32 {
        let atom = if a.unknown_symbol {
            Atom {
                symbol: AtomSymbol::Unknown,
                position: None,
                charge: a.charge,
                isotope: a.isotope,
                radical: None,
                hydrogen_count: a.hydrogen_count,
                implicit_h: a.implicit_h,
                aromatic: Some(a.aromatic),
                chirality: a.chirality,
                class: a.class,
                source_format: SourceFormat::SMILES,
            }
        } else {
            Atom {
                symbol: AtomSymbol::Element(a.element),
                position: None,
                isotope: a.isotope,
                radical: None,
                charge: a.charge,
                hydrogen_count: a.hydrogen_count,
                implicit_h: a.implicit_h,
                aromatic: Some(a.aromatic),
                chirality: a.chirality,
                class: a.class,
                source_format: SourceFormat::SMILES,
            }
        };

        self.atoms.push(atom);
        self.atoms.len() as u32
    }

    #[inline]
    pub fn on_bond(&mut self, start: u32, end: u32, b: BondData) {
        let bond = Bond {
            start_atom: start,
            end_atom: end,
            symbol: BondSymbol::Bond(b.order),
            direction: b.dir,
            ring: None,
            stereo: None,
            source_format: SourceFormat::SMILES,
        };
        self.bonds.push(bond);
    }

    // Fast-path constructors for hot parser loops
    #[inline]
    pub fn on_atom_fast(&mut self, element: Element, implicit_h: bool, aromatic: bool) -> u32 {
        let atom = Atom {
            symbol: AtomSymbol::Element(element),
            position: None,
            charge: None,
            isotope: None,
            radical: None,
            hydrogen_count: None,
            chirality: None,
            class: None,
            implicit_h,
            aromatic: Some(aromatic),
            source_format: SourceFormat::SMILES,
        };
        self.atoms.push(atom);
        self.atoms.len() as u32
    }

    #[inline]
    pub fn on_bond_single_fast(&mut self, start: u32, end: u32) {
        let bond = Bond {
            start_atom: start,
            end_atom: end,
            symbol: BondSymbol::Bond(BondOrder::Single),
            direction: None,
            ring: None,
            stereo: None,
            source_format: SourceFormat::SMILES,
        };
        self.bonds.push(bond);
    }

    pub fn on_component_end(&mut self) {
        if self.atoms.is_empty() {
            return;
        }
        let mut mol = Molecule::default();
        mol.atoms = mem::take(&mut self.atoms);
        mol.bonds = mem::take(&mut self.bonds);
        self.molecules.push(mol);
    }

    pub fn finish(&mut self) -> Vec<Molecule> {
        if !self.atoms.is_empty() || !self.bonds.is_empty() {
            self.on_component_end();
        }
        mem::take(&mut self.molecules)
    }
}
