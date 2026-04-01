//! Non-covalent bond representation for GraphIR.

use serde::{Deserialize, Serialize};

use super::molecule::AtomIndex;
use crate::bond::BondNoncovalent;
use crate::table_ir::bond::Bond as TableBond;

/// A non-covalent interaction (hydrogen bond, halogen bond, etc.) in GraphIR.
/// Non-covalent bonds are not stored in the main connectivity graph and do not
/// contribute to valence calculations.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct NoncovalentBond {
    a: AtomIndex,
    b: AtomIndex,
    kind: BondNoncovalent,
}

impl NoncovalentBond {
    pub fn new(a: AtomIndex, b: AtomIndex, kind: BondNoncovalent) -> Self {
        Self { a, b, kind }
    }

    /// Construct from a TableIR bond and the resolved node index map.
    ///
    /// `bond.noncovalent` must be `Some(_)`.
    pub fn from_table_bond(bond: &TableBond, node_indices: &[AtomIndex]) -> Self {
        let kind = bond
            .noncovalent
            .expect("NoncovalentBond::from_table_bond called on non-noncovalent bond");
        Self {
            a: node_indices[bond.atoms.first() as usize],
            b: node_indices[bond.atoms.second() as usize],
            kind,
        }
    }

    pub fn a(&self) -> AtomIndex {
        self.a
    }

    pub fn b(&self) -> AtomIndex {
        self.b
    }

    pub fn kind(&self) -> BondNoncovalent {
        self.kind
    }

    pub fn contains_atom(&self, atom: AtomIndex) -> bool {
        self.a == atom || self.b == atom
    }
}
