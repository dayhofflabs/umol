//! Molecule builders for SMILES parser

use std::mem;

use umol_shared::element::Element;

use crate::span::Span;
use crate::table_ir::atom::Chirality;
use crate::table_ir::{
    Atom, AtomPair, AtomSymbol, Bond, BondDonation, BondOrder, BondWedge, ExtendedAtom,
    ExtendedBond, ExtendedMolecule, Molecule, Ring, WildcardAtom,
};

/// Atom event data
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct AtomData {
    pub element: Element,
    pub isotope: Option<u32>,
    pub charge: Option<i8>,
    pub implicit_hydrogens: Option<u8>,
    pub class: Option<u32>,
    pub aromatic: bool,
    pub chirality: Option<Chirality>,
    pub span: Option<Span>,
}

/// Bond event data
pub(super) struct BondData {
    pub order: BondOrder,
    pub wedge: Option<BondWedge>,
    pub donation: Option<BondDonation>,
    pub span: Option<Span>,
}

/// Molecule builder
/// TODO: Consider removing Vec allocation and always return a single Molecule object
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

    #[inline]
    pub(crate) fn on_atom(&mut self, a: AtomData) -> u32 {
        let span = a.span;
        let atom = Atom {
            element: a.element,
            charge: a.charge,
            isotope_mass: a.isotope,
            implicit_hydrogens: a.implicit_hydrogens,
            valence: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic: Some(a.aromatic),
            chirality: a.chirality,
            class: a.class,
            label: None,
            value: None,
            span,
        };

        self.atoms.push(atom);
        (self.atoms.len() - 1) as u32
    }

    #[inline]
    pub(crate) fn on_bond(&mut self, start: u32, end: u32, b: BondData) {
        let span = b.span;
        // Adjust donation for AtomPair normalization (swap flips donation)
        let donation = if start > end {
            b.donation.map(|d| d.flip())
        } else {
            b.donation
        };
        let bond = Bond {
            atoms: AtomPair::new(start, end),
            order: b.order,
            charge: None,
            unpaired_electrons: None,
            multiplicity: None,
            wedge: b.wedge,
            donation,
            noncovalent: None,
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
        let mut atom = if aromatic {
            Atom::aromatic_atom(element)
        } else {
            Atom::aliphatic_atom(element)
        };
        atom.span = span;
        self.atoms.push(atom);
        (self.atoms.len() - 1) as u32
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

pub(super) struct ExtendedAtomData {
    pub symbol: AtomSymbol,
    pub isotope: Option<u32>,
    pub charge: Option<i8>,
    pub implicit_hydrogens: Option<u8>,
    pub class: Option<u32>,
    pub aromatic: bool,
    pub chirality: Option<Chirality>,
    pub span: Option<Span>,
}

/// TODO: Consider removing Vec allocation and always return a single Molecule object
pub(super) struct ExtendedMoleculeBuilder {
    atoms: Vec<ExtendedAtom>,
    bonds: Vec<ExtendedBond>,
    rings: Vec<Ring>,
    molecules: Vec<ExtendedMolecule>,
}

impl ExtendedMoleculeBuilder {
    pub(crate) fn with_capacity(approx_atoms: usize, approx_bonds: usize) -> Self {
        Self {
            atoms: Vec::with_capacity(approx_atoms),
            bonds: Vec::with_capacity(approx_bonds),
            rings: Vec::new(),
            molecules: Vec::new(),
        }
    }

    #[inline]
    pub(crate) fn on_atom(&mut self, a: ExtendedAtomData) -> u32 {
        let atom = ExtendedAtom {
            symbol: a.symbol,
            charge: a.charge,
            isotope_mass: a.isotope,
            implicit_hydrogens: a.implicit_hydrogens,
            valence: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic: Some(a.aromatic),
            chirality: a.chirality,
            class: a.class,
            label: None,
            value: None,
            pattern: None,
            stereo_care: None,
            inversion_retention: None,
            exact_change: None,
            attachment_point: None,
            attachment_order: None,
            ligand_order: None,
            ring_bond_count: None,
            substitution_count: None,
            unsaturated: None,
            link_atom: None,
            properties: Default::default(),
            span: a.span,
        };

        self.atoms.push(atom);
        (self.atoms.len() - 1) as u32
    }

    #[inline]
    pub(crate) fn on_bond(&mut self, start: u32, end: u32, b: BondData) {
        let mut bond = ExtendedBond::new(start, end, b.order);
        bond.wedge = b.wedge;
        // Adjust donation for AtomPair normalization (swap flips donation)
        bond.donation = if start > end {
            b.donation.map(|d| d.flip())
        } else {
            b.donation
        };
        bond.span = b.span;
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
        let atom = ExtendedAtom {
            symbol: AtomSymbol::Element(element),
            charge: None,
            isotope_mass: None,
            implicit_hydrogens: None,
            valence: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic: Some(aromatic),
            chirality: None,
            class: None,
            label: None,
            value: None,
            pattern: None,
            stereo_care: None,
            inversion_retention: None,
            exact_change: None,
            attachment_point: None,
            attachment_order: None,
            ligand_order: None,
            ring_bond_count: None,
            substitution_count: None,
            unsaturated: None,
            link_atom: None,
            properties: Default::default(),
            span,
        };
        self.atoms.push(atom);
        (self.atoms.len() - 1) as u32
    }

    #[inline]
    pub(crate) fn on_wildcard(
        &mut self,
        wildcard: WildcardAtom,
        class: Option<u32>,
        span_start: Option<u32>,
        span_end: Option<u32>,
    ) -> u32 {
        let span = Span::from_bytes_opt(span_start, span_end);
        let atom = ExtendedAtom {
            symbol: AtomSymbol::WildcardAtom(wildcard),
            charge: None,
            isotope_mass: None,
            implicit_hydrogens: None,
            valence: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic: None,
            chirality: None,
            class,
            label: None,
            value: None,
            pattern: None,
            stereo_care: None,
            inversion_retention: None,
            exact_change: None,
            attachment_point: None,
            attachment_order: None,
            ligand_order: None,
            ring_bond_count: None,
            substitution_count: None,
            unsaturated: None,
            link_atom: None,
            properties: Default::default(),
            span,
        };
        self.atoms.push(atom);
        (self.atoms.len() - 1) as u32
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
        let mut bond = ExtendedBond::new(start_atom, end_atom, BondOrder::Single);
        bond.span = span;
        self.bonds.push(bond);
    }

    pub(crate) fn on_component_end(&mut self) {
        if self.atoms.is_empty() {
            return;
        }
        let mut mol = ExtendedMolecule::empty();
        mol.atoms = mem::take(&mut self.atoms);
        mol.bonds = mem::take(&mut self.bonds);
        mol.rings = mem::take(&mut self.rings);
        self.molecules.push(mol);
    }

    pub(crate) fn finish(&mut self) -> Vec<ExtendedMolecule> {
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
