//! Molecule builders for SMILES parser

use std::mem;

use umol_shared::element::Element;

use super::super::error::ParseError;
use super::utils::{make_bond, make_extended_bond};
use crate::table_ir::{
    Atom, AtomSymbol, Bond, BondDonation, BondOrder, BondDirection, Chirality, ExtendedAtom,
    ExtendedBond, ExtendedMolecule, Molecule, Span, WildcardAtom,
};

/// Open ring-closure bond awaiting its matching digit. `bond_idx` is the
/// reserved entry in `bond_table`, filled in when the ring closes.
#[derive(Debug, Clone, Copy)]
struct OpenRing {
    atom_idx: usize,
    bond_idx: usize,
    order: Option<BondOrder>,
    direction: Option<BondDirection>,
    donation: Option<BondDonation>,
    open_pos: usize,
    open_end: usize,
}

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
    pub direction: Option<BondDirection>,
    pub donation: Option<BondDonation>,
    pub span: Option<Span>,
}

/// Molecule builder
pub(super) struct MoleculeBuilder {
    atoms: Vec<Atom>,
    bond_table: Vec<Option<Bond>>,
    ring_table: Vec<Option<OpenRing>>,
    /// (close_rank, open_index) per ring-closure bond, for CX bond-index remapping.
    ring_bonds: Vec<(usize, usize)>,
    /// Count of completed bonds; a bond's completion order is its CX close index.
    closed_bonds: usize,
    /// Whether to record ring closures (set only when a CX block is present).
    store_rings: bool,
    molecules: Vec<Molecule>,
}

impl MoleculeBuilder {
    pub(crate) fn with_capacity(
        approx_atoms: usize,
        approx_bonds: usize,
        store_rings: bool,
    ) -> Self {
        Self {
            atoms: Vec::with_capacity(approx_atoms),
            bond_table: Vec::with_capacity(approx_bonds),
            ring_table: Vec::new(),
            ring_bonds: Vec::new(),
            closed_bonds: 0,
            store_rings,
            molecules: Vec::new(),
        }
    }

    #[inline]
    pub(crate) fn on_atom(&mut self, a: AtomData) -> usize {
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
        self.atoms.len() - 1
    }

    #[inline]
    pub(crate) fn on_bond(&mut self, start: usize, end: usize, b: BondData) {
        self.bond_table.push(Some(make_bond(start, end, b)));
        self.closed_bonds += 1;
    }

    /// Reserve an empty bond entry at the current write position, returning its
    /// index. Completed later by [`on_ring_bond_close`] so that ring-closure
    /// bonds are recorded at their opening position, not their closing position.
    #[inline]
    fn on_ring_bond_open(&mut self) -> usize {
        self.bond_table.push(None);
        self.bond_table.len() - 1
    }

    #[inline]
    fn on_ring_bond_close(&mut self, bond_idx: usize, start: usize, end: usize, b: BondData) {
        self.bond_table[bond_idx] = Some(make_bond(start, end, b));
        if self.store_rings {
            self.ring_bonds.push((self.closed_bonds, bond_idx));
        }
        self.closed_bonds += 1;
    }

    /// Whether the atom at `atom_idx` is aromatic (false for wildcards / unset).
    pub(crate) fn is_aromatic(&self, atom_idx: usize) -> bool {
        self.atoms[atom_idx].aromatic == Some(true)
    }

    /// Process a ring-bond digit `ring_idx`: open it (reserve a bond entry) on first
    /// sight, close it (fill the reserved entry) on the matching second sight.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn on_ring_bond(
        &mut self,
        last_atom_idx: usize,
        ring_idx: usize,
        order_opt: Option<BondOrder>,
        direction_opt: Option<BondDirection>,
        donation_opt: Option<BondDonation>,
        pos: usize,
        token_end: usize,
        offset: usize,
    ) -> Result<(), ParseError> {
        if self.ring_table.len() <= ring_idx {
            self.ring_table.resize_with(ring_idx + 1, || None);
        }
        match self.ring_table[ring_idx].take() {
            None => {
                let bond_idx = self.on_ring_bond_open();
                self.ring_table[ring_idx] = Some(OpenRing {
                    atom_idx: last_atom_idx,
                    bond_idx,
                    order: order_opt,
                    direction: direction_opt,
                    donation: donation_opt,
                    open_pos: pos,
                    open_end: token_end,
                });
            }
            Some(open) => {
                // Once the close end's view is flipped (below), a consistent both-ends spec has
                // opposite raw symbols; equal raw symbols conflict.
                if let (Some(d1), Some(d2)) = (open.direction, direction_opt) {
                    if d1 == d2 {
                        return Err(ParseError::MismatchedRingBondDirections {
                            pos: offset + pos,
                            open_pos: offset + open.open_pos,
                        });
                    }
                }
                // Same donation on both ends = conflict (both donating or both receiving)
                if let (Some(don1), Some(don2)) = (open.donation, donation_opt) {
                    if don1 == don2 {
                        return Err(ParseError::MismatchedRingBondDonations {
                            pos: offset + pos,
                            open_pos: offset + open.open_pos,
                        });
                    }
                }
                if let (Some(o1), Some(o2)) = (open.order, order_opt) {
                    if o1 != o2 {
                        return Err(ParseError::MismatchedRingBondOrders {
                            pos: offset + pos,
                            open_pos: offset + open.open_pos,
                        });
                    }
                }
                let mut final_order = match (open.order, order_opt) {
                    (Some(o1), Some(o2)) => {
                        if o1 == o2 {
                            o1
                        } else {
                            o2
                        }
                    }
                    (Some(o), None) | (None, Some(o)) => o,
                    (None, None) => BondOrder::Single,
                };
                // Direction and donation: use the opening atom's perspective; if only the close
                // specifies it, flip it (it is from the closing atom's perspective).
                let final_direction = match (open.direction, direction_opt) {
                    (Some(d), _) => Some(d),
                    (None, Some(d)) => Some(d.flip()),
                    (None, None) => None,
                };
                let final_donation = match (open.donation, donation_opt) {
                    (Some(d), _) => Some(d),
                    (None, Some(d)) => Some(d.flip()),
                    (None, None) => None,
                };
                let a = open.atom_idx;
                let b = last_atom_idx;
                // Promote to aromatic only when the ring bond is implicit (no explicit bond token);
                // an explicit order or a directional /,\ keeps the bond as written.
                if open.order.is_none()
                    && order_opt.is_none()
                    && open.direction.is_none()
                    && direction_opt.is_none()
                    && self.is_aromatic(open.atom_idx)
                    && self.is_aromatic(b)
                {
                    final_order = BondOrder::Aromatic;
                }
                self.on_ring_bond_close(
                    open.bond_idx,
                    a,
                    b,
                    BondData {
                        order: final_order,
                        direction: final_direction,
                        donation: final_donation,
                        span: Span::from_bytes_opt(
                            Some(open.open_pos as u32),
                            Some(open.open_end as u32),
                        ),
                    },
                );
            }
        }
        Ok(())
    }

    /// Byte position of the latest still-open ring, if any (for reporting an
    /// unbalanced ring index once parsing is complete).
    pub(crate) fn unclosed_ring_pos(&self) -> Option<usize> {
        self.ring_table.iter().flatten().map(|o| o.open_pos).max()
    }

    #[inline]
    pub(crate) fn on_atom_fast(
        &mut self,
        element: Element,
        aromatic: bool,
        span_start: Option<u32>,
        span_end: Option<u32>,
    ) -> usize {
        let span = Span::from_bytes_opt(span_start, span_end);
        let mut atom = if aromatic {
            Atom::aromatic_atom(element)
        } else {
            Atom::aliphatic_atom(element)
        };
        atom.span = span;
        self.atoms.push(atom);
        self.atoms.len() - 1
    }

    #[inline]
    pub(crate) fn on_bond_single_fast(
        &mut self,
        start_atom: usize,
        end_atom: usize,
        span_start: Option<u32>,
        span_end: Option<u32>,
    ) {
        let span = Span::from_bytes_opt(span_start, span_end);
        let mut bond = Bond::new(start_atom as u32, end_atom as u32, BondOrder::Single);
        bond.span = span;
        self.bond_table.push(Some(bond));
        self.closed_bonds += 1;
    }

    pub(crate) fn on_component_end(&mut self) {
        if self.atoms.is_empty() {
            return;
        }
        let mut mol = Molecule::empty();
        mol.atoms = mem::take(&mut self.atoms);
        mol.bonds = self.bond_table.drain(..).flatten().collect();
        self.molecules.push(mol);
    }

    pub(crate) fn finish(&mut self) -> Vec<Molecule> {
        if !self.atoms.is_empty() || !self.bond_table.is_empty() {
            self.on_component_end();
        }
        mem::take(&mut self.molecules)
    }

    /// Take the recorded ring-closure bond records for CX bond-index remapping.
    pub(crate) fn take_ring_bonds(&mut self) -> Vec<(usize, usize)> {
        mem::take(&mut self.ring_bonds)
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

pub(super) struct ExtendedMoleculeBuilder {
    atoms: Vec<ExtendedAtom>,
    bond_table: Vec<Option<ExtendedBond>>,
    ring_table: Vec<Option<OpenRing>>,
    /// (close_rank, open_index) per ring-closure bond, for CX bond-index remapping.
    ring_bonds: Vec<(usize, usize)>,
    /// Count of completed bonds; a bond's completion order is its CX close index.
    closed_bonds: usize,
    /// Whether to record ring closures (set only when a CX block is present).
    store_rings: bool,
    molecules: Vec<ExtendedMolecule>,
}

impl ExtendedMoleculeBuilder {
    pub(crate) fn with_capacity(
        approx_atoms: usize,
        approx_bonds: usize,
        store_rings: bool,
    ) -> Self {
        Self {
            atoms: Vec::with_capacity(approx_atoms),
            bond_table: Vec::with_capacity(approx_bonds),
            ring_table: Vec::new(),
            ring_bonds: Vec::new(),
            closed_bonds: 0,
            store_rings,
            molecules: Vec::new(),
        }
    }

    #[inline]
    pub(crate) fn on_atom(&mut self, a: ExtendedAtomData) -> usize {
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
        self.atoms.len() - 1
    }

    #[inline]
    pub(crate) fn on_bond(&mut self, start: usize, end: usize, b: BondData) {
        self.bond_table
            .push(Some(make_extended_bond(start, end, b)));
        self.closed_bonds += 1;
    }

    /// Reserve an empty bond entry at the current write position, returning its
    /// index. Completed later by [`on_ring_bond_close`] so that ring-closure
    /// bonds are recorded at their opening position, not their closing position.
    #[inline]
    fn on_ring_bond_open(&mut self) -> usize {
        self.bond_table.push(None);
        self.bond_table.len() - 1
    }

    #[inline]
    fn on_ring_bond_close(&mut self, bond_idx: usize, start: usize, end: usize, b: BondData) {
        self.bond_table[bond_idx] = Some(make_extended_bond(start, end, b));
        if self.store_rings {
            self.ring_bonds.push((self.closed_bonds, bond_idx));
        }
        self.closed_bonds += 1;
    }

    /// Whether the atom at `atom_idx` is aromatic (false for wildcards / unset).
    pub(crate) fn is_aromatic(&self, atom_idx: usize) -> bool {
        self.atoms[atom_idx].aromatic == Some(true)
    }

    /// Process a ring-bond digit `ring_idx`: open it (reserve a bond entry) on first
    /// sight, close it (fill the reserved entry) on the matching second sight.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn on_ring_bond(
        &mut self,
        last_atom_idx: usize,
        ring_idx: usize,
        order_opt: Option<BondOrder>,
        direction_opt: Option<BondDirection>,
        donation_opt: Option<BondDonation>,
        pos: usize,
        token_end: usize,
        offset: usize,
    ) -> Result<(), ParseError> {
        if self.ring_table.len() <= ring_idx {
            self.ring_table.resize_with(ring_idx + 1, || None);
        }
        match self.ring_table[ring_idx].take() {
            None => {
                let bond_idx = self.on_ring_bond_open();
                self.ring_table[ring_idx] = Some(OpenRing {
                    atom_idx: last_atom_idx,
                    bond_idx,
                    order: order_opt,
                    direction: direction_opt,
                    donation: donation_opt,
                    open_pos: pos,
                    open_end: token_end,
                });
            }
            Some(open) => {
                // Once the close end's view is flipped (below), a consistent both-ends spec has
                // opposite raw symbols; equal raw symbols conflict.
                if let (Some(d1), Some(d2)) = (open.direction, direction_opt) {
                    if d1 == d2 {
                        return Err(ParseError::MismatchedRingBondDirections {
                            pos: offset + pos,
                            open_pos: offset + open.open_pos,
                        });
                    }
                }
                // Same donation on both ends = conflict (both donating or both receiving)
                if let (Some(don1), Some(don2)) = (open.donation, donation_opt) {
                    if don1 == don2 {
                        return Err(ParseError::MismatchedRingBondDonations {
                            pos: offset + pos,
                            open_pos: offset + open.open_pos,
                        });
                    }
                }
                if let (Some(o1), Some(o2)) = (open.order, order_opt) {
                    if o1 != o2 {
                        return Err(ParseError::MismatchedRingBondOrders {
                            pos: offset + pos,
                            open_pos: offset + open.open_pos,
                        });
                    }
                }
                let mut final_order = match (open.order, order_opt) {
                    (Some(o1), Some(o2)) => {
                        if o1 == o2 {
                            o1
                        } else {
                            o2
                        }
                    }
                    (Some(o), None) | (None, Some(o)) => o,
                    (None, None) => BondOrder::Single,
                };
                // Direction and donation: use the opening atom's perspective; if only the close
                // specifies it, flip it (it is from the closing atom's perspective).
                let final_direction = match (open.direction, direction_opt) {
                    (Some(d), _) => Some(d),
                    (None, Some(d)) => Some(d.flip()),
                    (None, None) => None,
                };
                let final_donation = match (open.donation, donation_opt) {
                    (Some(d), _) => Some(d),
                    (None, Some(d)) => Some(d.flip()),
                    (None, None) => None,
                };
                let a = open.atom_idx;
                let b = last_atom_idx;
                // Promote to aromatic only when the ring bond is implicit (no explicit bond token);
                // an explicit order or a directional /,\ keeps the bond as written.
                if open.order.is_none()
                    && order_opt.is_none()
                    && open.direction.is_none()
                    && direction_opt.is_none()
                    && self.is_aromatic(open.atom_idx)
                    && self.is_aromatic(b)
                {
                    final_order = BondOrder::Aromatic;
                }
                self.on_ring_bond_close(
                    open.bond_idx,
                    a,
                    b,
                    BondData {
                        order: final_order,
                        direction: final_direction,
                        donation: final_donation,
                        span: Span::from_bytes_opt(
                            Some(open.open_pos as u32),
                            Some(open.open_end as u32),
                        ),
                    },
                );
            }
        }
        Ok(())
    }

    /// Byte position of the latest still-open ring, if any (for reporting an
    /// unbalanced ring index once parsing is complete).
    pub(crate) fn unclosed_ring_pos(&self) -> Option<usize> {
        self.ring_table.iter().flatten().map(|o| o.open_pos).max()
    }

    #[inline]
    pub(crate) fn on_atom_fast(
        &mut self,
        element: Element,
        aromatic: bool,
        span_start: Option<u32>,
        span_end: Option<u32>,
    ) -> usize {
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
        self.atoms.len() - 1
    }

    #[inline]
    pub(crate) fn on_wildcard(
        &mut self,
        wildcard: WildcardAtom,
        class: Option<u32>,
        span_start: Option<u32>,
        span_end: Option<u32>,
    ) -> usize {
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
        self.atoms.len() - 1
    }

    #[inline]
    pub(crate) fn on_bond_single_fast(
        &mut self,
        start_atom: usize,
        end_atom: usize,
        span_start: Option<u32>,
        span_end: Option<u32>,
    ) {
        let span = Span::from_bytes_opt(span_start, span_end);
        let mut bond = ExtendedBond::new(start_atom as u32, end_atom as u32, BondOrder::Single);
        bond.span = span;
        self.bond_table.push(Some(bond));
        self.closed_bonds += 1;
    }

    pub(crate) fn on_component_end(&mut self) {
        if self.atoms.is_empty() {
            return;
        }
        let mut mol = ExtendedMolecule::empty();
        mol.atoms = mem::take(&mut self.atoms);
        mol.bonds = self.bond_table.drain(..).flatten().collect();
        self.molecules.push(mol);
    }

    pub(crate) fn finish(&mut self) -> Vec<ExtendedMolecule> {
        if !self.atoms.is_empty() || !self.bond_table.is_empty() {
            self.on_component_end();
        }
        mem::take(&mut self.molecules)
    }

    /// Take the recorded ring-closure bond records for CX bond-index remapping.
    pub(crate) fn take_ring_bonds(&mut self) -> Vec<(usize, usize)> {
        mem::take(&mut self.ring_bonds)
    }
}
