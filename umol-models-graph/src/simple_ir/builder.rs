//! MoleculeBuilder: IR molecule builder for parsers.

use std::mem;

use umol_data::Element;

use crate::simple_ir::{
    Atom, AtomSymbol, Bond, BondDir, BondOrder, BondSymbol, Chirality, Molecule, Ring,
};
use crate::span::Span;

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
    pub span: Option<Span>,
}

pub struct BondData {
    pub order: BondOrder,
    pub dir: Option<BondDir>,
    pub span: Option<Span>,
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
        let span = a.span;
        let atom = if a.unknown_symbol {
            Atom {
                symbol: AtomSymbol::Unknown,
                position: None,
                charge: a.charge,
                isotope: a.isotope,
                unpaired_e: None,
                hydrogens: a.hydrogen_count,
                implicit_h: a.implicit_h,
                aromatic: Some(a.aromatic),
                chirality: a.chirality,
                class: a.class,
                span,
            }
        } else {
            Atom {
                symbol: AtomSymbol::Element(a.element),
                position: None,
                isotope: a.isotope,
                unpaired_e: None,
                charge: a.charge,
                hydrogens: a.hydrogen_count,
                implicit_h: a.implicit_h,
                aromatic: Some(a.aromatic),
                chirality: a.chirality,
                class: a.class,
                span,
            }
        };

        self.atoms.push(atom);
        self.atoms.len() as u32
    }

    #[inline]
    pub fn on_bond(&mut self, start: u32, end: u32, b: BondData) {
        let span = b.span;
        let bond = Bond {
            start_atom: start,
            end_atom: end,
            symbol: BondSymbol::Bond(b.order),
            direction: b.dir,
            ring: None,
            stereo: None,
            span,
        };
        self.bonds.push(bond);
    }

    // Fast-path constructors for hot parser loops
    #[inline]
    pub fn on_atom_fast(
        &mut self,
        element: Element,
        aromatic: bool,
        span_start: Option<u32>,
        span_end: Option<u32>,
    ) -> u32 {
        let span = Span::from_bytes_opt(span_start, span_end);
        let atom = if aromatic {
            let mut a = Atom::from_aromatic_atom(element);
            a.span = span;
            a
        } else {
            let mut a = Atom::from_aliphatic_atom(element);
            a.span = span;
            a
        };
        self.atoms.push(atom);
        self.atoms.len() as u32
    }

    #[inline]
    pub fn on_bond_single_fast(
        &mut self,
        start_atom: u32,
        end_atom: u32,
        span_start: Option<u32>,
        span_end: Option<u32>,
    ) {
        let span = Span::from_bytes_opt(span_start, span_end);
        let mut bond = Bond::from_order(start_atom, end_atom, BondOrder::Single);
        bond.span = span;
        self.bonds.push(bond);
    }

    pub fn on_component_end(&mut self) {
        if self.atoms.is_empty() {
            return;
        }
        let mut mol = Molecule::default();
        mol.atoms = mem::take(&mut self.atoms);
        mol.bonds = mem::take(&mut self.bonds);
        mol.rings = mem::take(&mut self.rings);
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
            a.span = Some(match a.span {
                Some(span) => span.with_start(start),
                None => Span::bytes(start, start),
            });
        }
    }

    #[inline]
    pub fn annotate_last_bond_span(&mut self, start: u32) {
        if let Some(b) = self.bonds.last_mut() {
            b.span = Some(match b.span {
                Some(span) => span.with_start(start),
                None => Span::bytes(start, start),
            });
        }
    }

    #[inline]
    pub fn on_ring_open(
        &mut self,
        ring_idx: u32,
        start: Option<u32>,
        end: Option<u32>,
        atom_idx: Option<u32>,
    ) {
        self.rings.push(Ring {
            ring_idx,
            open_start: start,
            close_start: None,
            start_atom: atom_idx,
            end_atom: None,
            open_end: end,
            close_end: None,
        });
    }

    #[inline]
    pub fn on_ring_close(
        &mut self,
        ring_idx: u32,
        start: Option<u32>,
        end: Option<u32>,
        atom_idx: Option<u32>,
    ) {
        // find last open event for this index without close_start
        for ev in self.rings.iter_mut().rev() {
            if ev.ring_idx == ring_idx && ev.close_start.is_none() {
                ev.close_start = start;
                ev.end_atom = atom_idx;
                ev.close_end = end;
                return;
            }
        }
        // if none found, record a close-only event
        self.rings.push(Ring {
            ring_idx,
            open_start: None,
            close_start: start,
            start_atom: None,
            end_atom: atom_idx,
            open_end: None,
            close_end: end,
        });
    }
}
