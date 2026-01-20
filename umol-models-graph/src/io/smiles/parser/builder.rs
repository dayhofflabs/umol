//! Molecule builder for SMILES parser

use std::mem;

use umol_data::Element;

use crate::span::Span;
use crate::table_ir::{Atom, AtomPair, Bond, BondDirection, BondOrder, Chirality, Molecule, Ring};

/// Atom event data
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct AtomData {
    pub element: Element,
    pub isotope: Option<u32>,
    pub charge: Option<i8>,
    pub hydrogen_count: Option<u8>,
    pub class: Option<u32>,
    pub aromatic: bool,
    pub implicit_h: bool,
    pub chirality: Option<Chirality>,
    pub span: Option<Span>,
}

/// Bond event data
pub(super) struct BondData {
    pub order: BondOrder,
    pub direction: Option<BondDirection>,
    pub span: Option<Span>,
}

/// Molecule builder
pub(super) struct MoleculeBuilder {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    rings: Vec<Ring>,
    molecules: Vec<Molecule>,
}

impl MoleculeBuilder {
    pub(crate) fn with_capacity(approx_atoms: usize, approx_bonds: usize) -> Self {
        Self {
            atoms: Vec::with_capacity(approx_atoms),
            bonds: Vec::with_capacity(approx_bonds),
            rings: Vec::new(),
            molecules: Vec::new(),
        }
    }

    #[allow(dead_code)]
    pub(crate) fn clear_reuse(&mut self) {
        self.atoms.clear();
        self.bonds.clear();
        self.rings.clear();
        self.molecules.clear();
    }

    #[inline]
    pub(crate) fn on_atom(&mut self, a: AtomData) -> u32 {
        let span = a.span;
        let atom = Atom {
            element: a.element,
            charge: a.charge,
            isotope_mass: a.isotope,
            hydrogens: a.hydrogen_count,
            implicit_h: a.implicit_h,
            valence: None,
            unpaired_e: None,
            aromatic: Some(a.aromatic),
            chirality: a.chirality,
            class: a.class,
            label: None,
            value: None,
            span,
        };

        self.atoms.push(atom);
        self.atoms.len() as u32
    }

    #[inline]
    pub(crate) fn on_bond(&mut self, start: u32, end: u32, b: BondData) {
        let span = b.span;
        let bond = Bond {
            atoms: AtomPair::new(start, end),
            order: b.order,
            direction: b.direction,
            ring: None,
            stereo: None,
            span,
        };
        self.bonds.push(bond);
    }

    #[inline]
    pub(crate) fn on_atom_fast(
        &mut self,
        element: Element,
        aromatic: bool,
        span_start: Option<u32>,
        span_end: Option<u32>,
    ) -> u32 {
        let span = Span::from_bytes_opt(span_start, span_end);
        let atom = if aromatic {
            let mut a = Atom::aromatic_atom(element);
            a.span = span;
            a
        } else {
            let mut a = Atom::aliphatic_atom(element);
            a.span = span;
            a
        };
        self.atoms.push(atom);
        self.atoms.len() as u32
    }

    #[inline]
    pub(crate) fn on_bond_single_fast(
        &mut self,
        start_atom: u32,
        end_atom: u32,
        span_start: Option<u32>,
        span_end: Option<u32>,
    ) {
        let span = Span::from_bytes_opt(span_start, span_end);
        let mut bond = Bond::new(start_atom, end_atom, BondOrder::Single);
        bond.span = span;
        self.bonds.push(bond);
    }

    pub(crate) fn on_component_end(&mut self) {
        if self.atoms.is_empty() {
            return;
        }
        let mut mol = Molecule::empty();
        mol.atoms = mem::take(&mut self.atoms);
        mol.bonds = mem::take(&mut self.bonds);
        mol.rings = mem::take(&mut self.rings);
        self.molecules.push(mol);
    }

    pub(crate) fn finish(&mut self) -> Vec<Molecule> {
        if !self.atoms.is_empty() || !self.bonds.is_empty() {
            self.on_component_end();
        }
        mem::take(&mut self.molecules)
    }

    #[inline]
    pub(crate) fn on_ring_open(
        &mut self,
        ring_idx: u32,
        start: Option<u32>,
        end: Option<u32>,
        atom_idx: Option<u32>,
    ) {
        self.rings.push(Ring {
            ring_idx,
            open_span: Span::from_bytes_opt(start, end),
            close_span: None,
            start_atom: atom_idx,
            end_atom: None,
        });
    }

    #[inline]
    pub(crate) fn on_ring_close(
        &mut self,
        ring_idx: u32,
        start: Option<u32>,
        end: Option<u32>,
        atom_idx: Option<u32>,
    ) {
        for ev in self.rings.iter_mut().rev() {
            if ev.ring_idx == ring_idx && ev.close_span.is_none() {
                ev.close_span = Span::from_bytes_opt(start, end);
                ev.end_atom = atom_idx;
                return;
            }
        }

        self.rings.push(Ring {
            ring_idx,
            open_span: None,
            close_span: Span::from_bytes_opt(start, end),
            start_atom: None,
            end_atom: atom_idx,
        });
    }
}
