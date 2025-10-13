//! MoleculeBuilder: IR molecule builder for parsers.

use crate::io::ir::{AtomSymbol, SourceFormat};
use std::mem;

use super::{Atom, Bond, BondDir, BondOrder, BondSymbol, Chirality, Molecule, Ring};
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
    pub span_start: Option<u32>,
}

pub struct BondData {
    pub order: BondOrder,
    pub dir: Option<BondDir>,
    pub span_start: Option<u32>,
}

pub struct MoleculeBuilder {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    rings: Vec<Ring>,
    molecules: Vec<Molecule>,
}

impl MoleculeBuilder {
    pub fn with_capacity(approx_atoms: usize, approx_bonds: usize) -> Self {
        Self {
            atoms: Vec::with_capacity(approx_atoms),
            bonds: Vec::with_capacity(approx_bonds),
            rings: Vec::new(),
            molecules: Vec::new(),
        }
    }

    pub fn clear_reuse(&mut self) {
        self.atoms.clear();
        self.bonds.clear();
        self.rings.clear();
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
                span_start: a.span_start,
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
                span_start: a.span_start,
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
            span_start: b.span_start,
            source_format: SourceFormat::SMILES,
        };
        self.bonds.push(bond);
    }

    // Fast-path constructors for hot parser loops
    #[inline]
    pub fn on_atom_fast(&mut self, element: Element, implicit_h: bool, aromatic: bool, span_start: Option<u32>) -> u32 {
        let atom = Atom {
            symbol: AtomSymbol::Element(element),
            position: None,
            span_start,
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
    pub fn on_bond_single_fast(&mut self, start: u32, end: u32, span_start: Option<u32>) {
        let bond = Bond {
            start_atom: start,
            end_atom: end,
            symbol: BondSymbol::Bond(BondOrder::Single),
            direction: None,
            ring: None,
            stereo: None,
            span_start,
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
        mol.ring_events = mem::take(&mut self.rings);
        self.molecules.push(mol);
    }

    pub fn finish(&mut self) -> Vec<Molecule> {
        if !self.atoms.is_empty() || !self.bonds.is_empty() {
            self.on_component_end();
        }
        mem::take(&mut self.molecules)
    }

    #[inline]
    pub fn annotate_last_atom_span(&mut self, start: u32) {
        if let Some(a) = self.atoms.last_mut() {
            a.span_start = Some(start);
        }
    }

    #[inline]
    pub fn annotate_last_bond_span(&mut self, start: u32) {
        if let Some(b) = self.bonds.last_mut() {
            b.span_start = Some(start);
        }
    }

    #[inline]
    pub fn on_ring_open(&mut self, ring_idx: u32, pos: u32, atom_id: u32) {
        self.rings.push(Ring {
            ring_idx,
            open_pos: Some(pos),
            close_pos: None,
            atom_a: Some(atom_id),
            atom_b: None,
        });
    }

    #[inline]
    pub fn on_ring_close(&mut self, ring_idx: u32, pos: u32, atom_id: u32) {
        // find last open event for this index without close_pos
        for ev in self.rings.iter_mut().rev() {
            if ev.ring_idx == ring_idx && ev.close_pos.is_none() {
                ev.close_pos = Some(pos);
                ev.atom_b = Some(atom_id);
                return;
            }
        }
        // if none found, record a close-only event
        self.rings.push(Ring {
            ring_idx,
            open_pos: None,
            close_pos: Some(pos),
            atom_a: None,
            atom_b: Some(atom_id),
        });
    }
}
