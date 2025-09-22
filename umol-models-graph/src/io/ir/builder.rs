//! MoleculeBuilder: fast, remap-free IR construction for parsers.

use super::{Atom as IRAtom, Bond as IRBond, BondDir, BondOrder, Molecule as IRMolecule, SourceFormat};
use umol_data::Element;

pub struct AtomData {
    pub element: Element,
    pub isotope: Option<u32>,
    pub charge: Option<i32>,
    pub hydrogen_count: Option<u8>,
    pub aromatic: bool,
    pub implicit_h: bool,
    pub chirality: Option<super::Chirality>,
}

pub struct BondData {
    pub order: BondOrder,
    pub dir: Option<BondDir>,
}

pub struct MoleculeBuilder {
    atoms: Vec<IRAtom>,
    bonds: Vec<IRBond>,
    molecules: Vec<IRMolecule>,
}

impl MoleculeBuilder {
    pub fn with_capacity(approx_atoms: usize, approx_bonds: usize) -> Self {
        Self { atoms: Vec::with_capacity(approx_atoms), bonds: Vec::with_capacity(approx_bonds), molecules: Vec::new() }
    }

    pub fn clear_reuse(&mut self) {
        self.atoms.clear();
        self.bonds.clear();
        self.molecules.clear();
    }

    pub fn on_atom(&mut self, a: AtomData) -> u32 {
        let mut atom = IRAtom::from_element(a.element);
        atom.isotope = a.isotope;
        atom.charge = a.charge;
        atom.hydrogen_count = a.hydrogen_count.map(|v| v as u32);
        atom.aromatic = Some(a.aromatic);
        atom.implicit_h = a.implicit_h;
        atom.chirality = a.chirality;
        atom.source_format = SourceFormat::SMILES;
        let idx = self.atoms.len() as u32;
        // Keep index for compatibility now
        atom.index = Some(idx);
        self.atoms.push(atom);
        idx
    }

    pub fn on_bond(&mut self, start: u32, end: u32, b: BondData) {
        let mut bond = IRBond::from_order(b.order);
        bond.start_atom = Some(start);
        bond.end_atom = Some(end);
        bond.direction = b.dir;
        bond.source_format = SourceFormat::SMILES;
        self.bonds.push(bond);
    }

    pub fn on_component_end(&mut self) {
        if self.atoms.is_empty() { return; }
        let mut mol = IRMolecule::default();
        mol.source_format = SourceFormat::SMILES;
        mol.atoms = std::mem::take(&mut self.atoms);
        mol.bonds = std::mem::take(&mut self.bonds);
        self.molecules.push(mol);
    }

    pub fn finish(&mut self) -> Vec<IRMolecule> {
        if !self.atoms.is_empty() || !self.bonds.is_empty() { self.on_component_end(); }
        std::mem::take(&mut self.molecules)
    }
}


