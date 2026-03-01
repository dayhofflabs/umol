//! Bond type for GraphIR.

use super::error::ResolutionError;
use crate::table_ir::bond::{Bond as TableBond, BondOrder};

/// Resolved shared (covalent) bond in GraphIR. Order is the localized (σ-skeleton)
/// bond order. Dative and non-covalent bonds are stored separately.
#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    order: u8,
}

impl Bond {
    pub fn new(order: u8) -> Self {
        Self { order }
    }

    pub fn order(&self) -> u8 {
        self.order
    }

    /// Construct from a shared TableIR bond.
    ///
    /// Note: Aromaticity is handled separately; aromatic bonds are mapped to
    /// single bonds. Query bonds are not allowed. Dative and non-covalent bonds
    /// must be routed to their respective constructors before calling this.
    pub fn from_table_bond(bond: &TableBond) -> Result<Bond, ResolutionError> {
        let order = match bond.order {
            BondOrder::Aromatic => 1,
            o if o.is_query() => return Err(ResolutionError::InvalidBondOrder(bond.order)),
            o => o
                .value()
                .ok_or(ResolutionError::InvalidBondOrder(bond.order))?,
        };
        Ok(Bond::new(order))
    }
}
