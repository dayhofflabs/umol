//! Dative (coordinate) bond representation for GraphIR.

use crate::graph_ir::molecule::AtomIndex;
use crate::table_ir::bond::{Bond as TableBond, BondDonation};

/// A dative (coordinate) bond in GraphIR. Carries the donor and acceptor atom
/// indices and the bond order (typically 1). Unlike shared bonds, dative bonds
/// are not stored in the main connectivity graph and do not contribute to
/// `atom_bond_order_sum`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct DativeBond {
    donor: AtomIndex,
    acceptor: AtomIndex,
    order: u8,
}

impl DativeBond {
    pub fn new(donor: AtomIndex, acceptor: AtomIndex, order: u8) -> Self {
        Self {
            donor,
            acceptor,
            order,
        }
    }

    /// Construct from a TableIR bond and the resolved node index map.
    ///
    /// `bond.donation` must be `Some(Donating)` or `Some(Accepting)`.
    /// Direction relative to `bond.atoms` (which normalizes to `first <= second`):
    /// - `Donating`: first atom donates → first is donor
    /// - `Accepting`: first atom accepts → second is donor
    pub fn from_table_bond(bond: &TableBond, node_indices: &[AtomIndex]) -> Self {
        let first = bond.atoms.first() as usize;
        let second = bond.atoms.second() as usize;
        let order = bond.order.value().unwrap_or(1);
        let (donor, acceptor) = match bond.donation {
            Some(BondDonation::Donating) | Some(BondDonation::Shared) => {
                (node_indices[first], node_indices[second])
            }
            Some(BondDonation::Accepting) => (node_indices[second], node_indices[first]),
            None => panic!("DativeBond::from_table_bond called on non-dative bond"),
        };
        Self {
            donor,
            acceptor,
            order,
        }
    }

    pub fn donor(&self) -> AtomIndex {
        self.donor
    }

    pub fn acceptor(&self) -> AtomIndex {
        self.acceptor
    }

    pub fn order(&self) -> u8 {
        self.order
    }

    pub fn contains_atom(&self, atom: AtomIndex) -> bool {
        self.donor == atom || self.acceptor == atom
    }
}
