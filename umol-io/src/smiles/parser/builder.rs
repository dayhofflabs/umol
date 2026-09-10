//! Molecule builder for SMILES parser

use std::collections::BTreeMap;
use std::iter;

use smallvec::SmallVec;
use umol_chem::element::Element;

use super::super::error::ParseError;
use super::utils::{invalid_ring_context, make_bond, make_extended_bond, Frame};
use crate::table_ir::{
    Atom, AtomSymbol, Bond, BondDirection, BondDonation, BondOrder, Chirality, ExtendedAtom,
    ExtendedBond, ExtendedMolecule, Molecule, SourceFormat, Span, StereoAtom, StereoLigand,
    WildcardAtom, Winding,
};

#[derive(Debug, Clone, Copy)]
pub(super) struct Cursor {
    atom: usize,
    aromatic: bool,
    stereo: Option<usize>,
}

type BondSpec = (
    BondOrder,
    Option<BondDirection>,
    Option<BondDonation>,
    usize,
);
pub(super) type AtomMapping = BTreeMap<u32, (Vec<u32>, Vec<u32>)>;

struct PendingStereo {
    atom: u32,
    winding: Winding,
    is_root: bool,
    bonds: SmallVec<[usize; 4]>,
}

#[derive(Debug, Clone, Copy)]
struct OpenRing {
    atom_idx: usize,
    aromatic: bool,
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
    pub element: Option<Element>,
    pub isotope: Option<u32>,
    pub charge: Option<i8>,
    pub implicit_hydrogens: Option<u8>,
    pub class: Option<u32>,
    pub aromatic: Option<bool>,
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

pub(super) struct ExtendedAtomData {
    pub symbol: AtomSymbol,
    pub isotope: Option<u32>,
    pub charge: Option<i8>,
    pub implicit_hydrogens: Option<u8>,
    pub class: Option<u32>,
    pub aromatic: Option<bool>,
    pub chirality: Option<Chirality>,
    pub span: Option<Span>,
}

pub(super) trait Target {
    type Atom;
    type Bond;
    type Molecule;
    type Input;
    const EXTENDED: bool;
    fn bracket(
        a: Self::Input,
    ) -> (
        Self::Atom,
        bool,
        Option<Chirality>,
        Option<u32>,
        Option<Span>,
    );
    fn organic(element: Element, aromatic: bool, span: Option<Span>) -> Self::Atom;
    fn wildcard(span: Option<Span>) -> Self::Atom;
    fn bond(start: usize, end: usize, data: BondData) -> Self::Bond;
    fn bond_atoms(bond: &Self::Bond) -> (u32, u32);
    fn hydrogens(atom: &Self::Atom) -> u8;
    fn molecule(
        atoms: Vec<Self::Atom>,
        bonds: Vec<Self::Bond>,
        stereo: Vec<StereoAtom>,
    ) -> Self::Molecule;
}

pub(super) struct Basic;
pub(super) struct Extended;
pub(super) type MoleculeEditor<'a> = Builder<'a, Basic>;
pub(super) type ExtendedMoleculeBuilder<'a> = Builder<'a, Extended>;

impl Target for Basic {
    type Atom = Atom;
    type Bond = Bond;
    type Molecule = Molecule;
    type Input = AtomData;
    const EXTENDED: bool = false;

    fn bracket(
        a: Self::Input,
    ) -> (
        Self::Atom,
        bool,
        Option<Chirality>,
        Option<u32>,
        Option<Span>,
    ) {
        let atom = Atom {
            element: a.element,
            charge: a.charge,
            isotope_mass: a.isotope,
            implicit_hydrogens: a.implicit_hydrogens,
            valence: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic: a.aromatic,
            chirality: a.chirality,
            class: a.class,
            label: None,
            value: None,
            span: a.span,
        };

        (atom, a.aromatic == Some(true), a.chirality, a.class, a.span)
    }

    fn organic(element: Element, aromatic: bool, span: Option<Span>) -> Self::Atom {
        let mut atom = Atom::from_element(element);
        atom.aromatic = Some(aromatic);
        atom.span = span;
        atom
    }

    fn wildcard(span: Option<Span>) -> Self::Atom {
        let mut atom = Atom::wildcard();
        atom.span = span;
        atom
    }

    fn bond(start: usize, end: usize, data: BondData) -> Self::Bond {
        make_bond(start, end, data)
    }

    fn bond_atoms(bond: &Self::Bond) -> (u32, u32) {
        (bond.atoms.first(), bond.atoms.second())
    }

    fn hydrogens(atom: &Self::Atom) -> u8 {
        atom.implicit_hydrogens.unwrap_or(0)
    }

    fn molecule(
        atoms: Vec<Self::Atom>,
        bonds: Vec<Self::Bond>,
        stereo: Vec<StereoAtom>,
    ) -> Self::Molecule {
        let mut mol = Molecule::empty();
        if !atoms.is_empty() {
            mol.source_format = SourceFormat::SMILES;
        }
        mol.atoms = atoms;
        mol.bonds = bonds;
        mol.stereo_atoms = stereo;
        mol
    }
}

impl Target for Extended {
    type Atom = ExtendedAtom;
    type Bond = ExtendedBond;
    type Molecule = ExtendedMolecule;
    type Input = ExtendedAtomData;
    const EXTENDED: bool = true;

    fn bracket(
        a: Self::Input,
    ) -> (
        Self::Atom,
        bool,
        Option<Chirality>,
        Option<u32>,
        Option<Span>,
    ) {
        let atom = ExtendedAtom {
            symbol: a.symbol,
            charge: a.charge,
            isotope_mass: a.isotope,
            implicit_hydrogens: a.implicit_hydrogens,
            valence: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic: a.aromatic,
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

        (atom, a.aromatic == Some(true), a.chirality, a.class, a.span)
    }

    fn organic(element: Element, aromatic: bool, span: Option<Span>) -> Self::Atom {
        let mut atom = ExtendedAtom::from_element(element);
        atom.aromatic = Some(aromatic);
        atom.span = span;
        atom
    }

    fn wildcard(span: Option<Span>) -> Self::Atom {
        let mut atom = ExtendedAtom::from_atom_symbol(AtomSymbol::WildcardAtom(WildcardAtom::Any));
        atom.span = span;
        atom
    }

    fn bond(start: usize, end: usize, data: BondData) -> Self::Bond {
        make_extended_bond(start, end, data)
    }

    fn bond_atoms(bond: &Self::Bond) -> (u32, u32) {
        (bond.atoms.first(), bond.atoms.second())
    }

    fn hydrogens(atom: &Self::Atom) -> u8 {
        atom.implicit_hydrogens.unwrap_or(0)
    }

    fn molecule(
        atoms: Vec<Self::Atom>,
        bonds: Vec<Self::Bond>,
        stereo: Vec<StereoAtom>,
    ) -> Self::Molecule {
        let mut mol = ExtendedMolecule::empty();
        if !atoms.is_empty() {
            mol.source_format = SourceFormat::SMILES;
        }
        mol.atoms = atoms;
        mol.bonds = bonds;
        mol.stereo_atoms = stereo;
        mol
    }
}

pub(super) struct Builder<'a, T: Target> {
    atoms: Vec<T::Atom>,
    bond_table: Vec<Option<T::Bond>>,
    ring_table: Vec<Option<OpenRing>>,
    ring_bonds: Vec<(usize, usize)>,
    open_rings: usize,
    store_rings: bool,
    stereo: Vec<PendingStereo>,
    current: Option<Cursor>,
    branches: Vec<Frame>,
    pending_bond: Option<BondSpec>,
    after_closed_group: bool,
    mapping: Option<(&'a mut AtomMapping, bool)>,
}

impl<'a, T: Target> Builder<'a, T> {
    pub(super) fn with_capacity(
        approx_atoms: usize,
        approx_bonds: usize,
        store_rings: bool,
        mapping: Option<(&'a mut AtomMapping, bool)>,
    ) -> Self {
        Self {
            atoms: Vec::with_capacity(approx_atoms),
            bond_table: Vec::with_capacity(approx_bonds),
            ring_table: Vec::new(),
            ring_bonds: Vec::new(),
            open_rings: 0,
            store_rings,
            stereo: Vec::new(),
            current: None,
            branches: Vec::new(),
            pending_bond: None,
            after_closed_group: false,
            mapping,
        }
    }

    #[inline]
    pub(super) fn token(&mut self, byte: u8) {
        if byte != b'(' {
            self.after_closed_group = false;
        }
    }

    pub(super) fn open_branch(&mut self, pos: usize, offset: usize) -> Result<(), ParseError> {
        if let Some((_, _, _, pos)) = self.pending_bond {
            return Err(ParseError::TrailingBond { pos: offset + pos });
        }
        if self.after_closed_group {
            self.current = None;
            self.branches.push(Frame::Group {
                had_atom: false,
                open_pos: pos,
            });
            self.after_closed_group = false;
        } else {
            self.branches.push(match self.current {
                Some(base) => Frame::Branch {
                    base,
                    had_atom: false,
                    open_pos: pos,
                },
                None => Frame::Group {
                    had_atom: false,
                    open_pos: pos,
                },
            });
        }
        Ok(())
    }

    pub(super) fn close_branch(
        &mut self,
        input: &[u8],
        pos: usize,
        offset: usize,
    ) -> Result<(), ParseError> {
        if let Some((_, _, _, pos)) = self.pending_bond {
            return Err(ParseError::TrailingBond { pos: offset + pos });
        }
        let Some(frame) = self.branches.pop() else {
            return Err(ParseError::UnbalancedCloseParen { pos: offset + pos });
        };
        match frame {
            Frame::Branch { base, had_atom, .. } => {
                if !had_atom {
                    return Err(ParseError::EmptyBranch { pos: offset + pos });
                }
                self.current = Some(base);
            }
            Frame::Group { had_atom, .. } => {
                if !had_atom {
                    return Err(ParseError::EmptyGroup { pos: offset + pos });
                }
                self.after_closed_group = true;
                self.mark_branch();
                if self.branches.is_empty() && pos + 1 != input.len() && input[pos + 1] != b'.' {
                    return Err(ParseError::NonfinalGroup { pos: offset + pos });
                }
            }
        }
        Ok(())
    }

    pub(super) fn dot(
        &mut self,
        input: &[u8],
        pos: usize,
        offset: usize,
        as_reaction: bool,
    ) -> Result<(), ParseError> {
        let local_offset = if T::EXTENDED { 0 } else { offset };
        if let Some((_, _, _, pos)) = self.pending_bond {
            return Err(ParseError::TrailingBond {
                pos: local_offset + pos,
            });
        }
        if pos == 0
            || matches!(
                self.branches.last(),
                Some(Frame::Group {
                    had_atom: false,
                    ..
                })
            )
        {
            return Err(ParseError::LeadingDot {
                pos: local_offset + pos,
            });
        }
        if pos + 1 == input.len() {
            return Err(ParseError::TrailingDot {
                pos: local_offset + pos,
            });
        }
        if as_reaction && input[pos + 1] == b'>' {
            return Err(ParseError::TrailingDot { pos: offset + pos });
        }
        if input[pos + 1] == b'.' {
            return Err(ParseError::ConsecutiveDots { pos: offset + pos });
        }
        if input[pos + 1].is_ascii_digit() || input[pos + 1] == b'%' {
            return Err(ParseError::DotBeforeRing { pos: offset + pos });
        }
        self.current = None;
        Ok(())
    }

    pub(super) fn on_bond(
        &mut self,
        order: BondOrder,
        direction: Option<BondDirection>,
        donation: Option<BondDonation>,
        pos: usize,
        offset: usize,
    ) -> Result<(), ParseError> {
        if self.pending_bond.is_some() {
            return Err(ParseError::ConsecutiveBonds { pos: offset + pos });
        }
        if self.current.is_none() {
            return Err(ParseError::LeadingBond { pos: offset + pos });
        }
        self.pending_bond = Some((order, direction, donation, pos));
        Ok(())
    }

    #[inline]
    fn mark_branch(&mut self) {
        if let Some(Frame::Branch { had_atom, .. } | Frame::Group { had_atom, .. }) =
            self.branches.last_mut()
        {
            *had_atom = true;
        }
    }

    #[inline]
    pub(super) fn on_atom(&mut self, data: T::Input) {
        let (atom, aromatic, chirality, class, span) = T::bracket(data);
        let index = self.atoms.len();
        self.atoms.push(atom);
        if let (Some(class), Some((mapping, is_product))) = (class, self.mapping.as_mut()) {
            let entry = mapping.entry(class).or_default();
            if *is_product {
                entry.1.push(index as u32);
            } else {
                entry.0.push(index as u32);
            }
        }
        let winding = match chirality {
            Some(Chirality::CounterClockwise | Chirality::Tetrahedral { arr: 1 }) => {
                Some(Winding::CounterClockwise)
            }
            Some(Chirality::Clockwise | Chirality::Tetrahedral { arr: 2 }) => {
                Some(Winding::Clockwise)
            }
            _ => None,
        };
        let stereo = winding.map(|winding| {
            let index = self.stereo.len();
            self.stereo.push(PendingStereo {
                atom: (self.atoms.len() - 1) as u32,
                winding,
                is_root: self.current.is_none(),
                bonds: SmallVec::new(),
            });
            index
        });
        self.attach(
            Cursor {
                atom: index,
                aromatic,
                stereo,
            },
            span,
        );
    }

    #[inline]
    pub(super) fn on_atom_fast(
        &mut self,
        element: Element,
        aromatic: bool,
        span_start: Option<u32>,
        span_end: Option<u32>,
    ) {
        let span = Span::from_bytes_opt(span_start, span_end);
        let index = self.atoms.len();
        self.atoms.push(T::organic(element, aromatic, span));
        self.attach(
            Cursor {
                atom: index,
                aromatic,
                stereo: None,
            },
            span,
        );
    }

    #[inline]
    pub(super) fn on_wildcard(&mut self, span_start: Option<u32>, span_end: Option<u32>) {
        let span = Span::from_bytes_opt(span_start, span_end);
        let index = self.atoms.len();
        self.atoms.push(T::wildcard(span));
        self.attach(
            Cursor {
                atom: index,
                aromatic: false,
                stereo: None,
            },
            span,
        );
    }

    #[inline]
    fn incidence(&mut self, cursor: Cursor, bond: usize) {
        if let Some(stereo) = cursor.stereo {
            self.stereo[stereo].bonds.push(bond);
        }
    }

    #[inline]
    fn attach(&mut self, current: Cursor, span: Option<Span>) {
        if let Some(previous) = self.current {
            let data = match self.pending_bond.take() {
                Some((order, direction, donation, pos)) => BondData {
                    order,
                    direction,
                    donation,
                    span: Span::from_bytes_opt(Some(pos as u32), Some(pos as u32 + 1)),
                },
                None => BondData {
                    order: if previous.aromatic && current.aromatic {
                        BondOrder::Aromatic
                    } else {
                        BondOrder::Single
                    },
                    direction: None,
                    donation: None,
                    span,
                },
            };
            let bond = self.append_bond(T::bond(previous.atom, current.atom, data));
            self.incidence(previous, bond);
            self.incidence(current, bond);
        }
        self.current = Some(current);
        self.mark_branch();
    }

    #[inline]
    fn append_bond(&mut self, bond: T::Bond) -> usize {
        let index = self.bond_table.len();
        self.bond_table.push(Some(bond));
        index
    }

    #[inline]
    fn reserve_bond(&mut self) -> usize {
        let index = self.bond_table.len();
        self.bond_table.push(None);
        index
    }

    #[inline]
    fn complete_bond(&mut self, index: usize, bond: T::Bond) {
        self.bond_table[index] = Some(bond);
    }

    pub(super) fn on_ring_bond(
        &mut self,
        ring_idx: usize,
        pos: usize,
        token_end: usize,
        offset: usize,
    ) -> Result<(), ParseError> {
        let Some(current) = self.current else {
            return Err(ParseError::LeadingRing { pos: offset + pos });
        };
        if invalid_ring_context(&self.branches) {
            return Err(ParseError::LeadingRing { pos: offset });
        }
        let (order_opt, direction_opt, donation_opt) = self
            .pending_bond
            .take()
            .map_or((None, None, None), |(order, direction, donation, _)| {
                (Some(order), direction, donation)
            });
        if self.ring_table.len() <= ring_idx {
            self.ring_table.resize_with(ring_idx + 1, || None);
        }
        match self.ring_table[ring_idx].take() {
            None => {
                let bond_idx = self.reserve_bond();
                self.open_rings += 1;
                self.incidence(current, bond_idx);
                self.ring_table[ring_idx] = Some(OpenRing {
                    atom_idx: current.atom,
                    aromatic: current.aromatic,
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
                let b = current.atom;
                // Promote to aromatic only when the ring bond is implicit (no explicit bond token);
                // an explicit order or a directional /,\ keeps the bond as written.
                if open.order.is_none()
                    && order_opt.is_none()
                    && open.direction.is_none()
                    && direction_opt.is_none()
                    && open.aromatic
                    && current.aromatic
                {
                    final_order = BondOrder::Aromatic;
                }
                self.complete_bond(
                    open.bond_idx,
                    T::bond(
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
                    ),
                );
                self.open_rings -= 1;
                self.incidence(current, open.bond_idx);
                if self.store_rings {
                    self.ring_bonds
                        .push((self.bond_table.len() - self.open_rings - 1, open.bond_idx));
                }
            }
        }
        Ok(())
    }

    #[allow(clippy::type_complexity)]
    pub(super) fn finish(
        self,
        offset: usize,
    ) -> Result<(T::Molecule, Vec<(usize, usize)>), ParseError> {
        if let Some((_, _, _, pos)) = self.pending_bond {
            return Err(ParseError::TrailingBond { pos: offset + pos });
        }
        if let Some(Frame::Branch { open_pos, .. } | Frame::Group { open_pos, .. }) =
            self.branches.last()
        {
            return Err(ParseError::UnbalancedOpenParen {
                pos: offset + open_pos,
            });
        }
        if let Some(open_pos) = self
            .ring_table
            .iter()
            .flatten()
            .map(|ring| ring.open_pos)
            .max()
        {
            return Err(ParseError::UnbalancedRingIndex {
                open_pos: offset + open_pos,
            });
        }
        let bonds: Vec<_> = self
            .bond_table
            .into_iter()
            .map(|bond| bond.expect("all ring slots completed"))
            .collect();
        let stereo = self
            .stereo
            .into_iter()
            .map(|pending| {
                let hydrogens = T::hydrogens(&self.atoms[pending.atom as usize]) as usize;
                let lone_pair = hydrogens == 0 && pending.bonds.len() == 3;
                let virtual_count = hydrogens + usize::from(lone_pair);
                let position = usize::from(!pending.is_root);
                let actual_count = pending.bonds.len();
                let mut ligands = Vec::with_capacity(actual_count + virtual_count);
                for (index, bond) in pending.bonds.into_iter().enumerate() {
                    if index == position {
                        ligands.extend(iter::repeat_n(StereoLigand::ImplicitHydrogen, hydrogens));
                        if lone_pair {
                            ligands.push(StereoLigand::LonePair);
                        }
                    }
                    let (a, b) = T::bond_atoms(&bonds[bond]);
                    ligands.push(StereoLigand::Atom(if a == pending.atom { b } else { a }));
                }
                if position == actual_count {
                    ligands.extend(iter::repeat_n(StereoLigand::ImplicitHydrogen, hydrogens));
                }
                StereoAtom {
                    atom: pending.atom,
                    ligands,
                    winding: pending.winding,
                }
            })
            .collect();
        Ok((T::molecule(self.atoms, bonds, stereo), self.ring_bonds))
    }
}

#[cfg(test)]
mod tests;
