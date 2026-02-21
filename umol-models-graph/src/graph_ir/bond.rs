//! Bond type for GraphIR.

use super::error::ResolutionError;
use crate::bond::BondDonation;
use crate::table_ir::bond::{Bond as TableBond, BondOrder};

/// Resolved bond in GraphIR. Order is the localized (σ-skeleton) bond order.
#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    order: u8,
    donation: Option<BondDonation>,
}

impl Bond {
    pub fn new(order: u8) -> Self {
        Self {
            order,
            donation: None,
        }
    }

    pub fn new_dative(order: u8, donation: BondDonation) -> Self {
        Self {
            order,
            donation: Some(donation),
        }
    }

    pub fn order(&self) -> u8 {
        self.order
    }

    pub fn donation(&self) -> Option<BondDonation> {
        self.donation
    }

    /// Construct from table bond
    /// Note: Aromaticity is handled separately, bond order here refers *only*
    /// to the sigma skeleton. Aromatic bonds are mapped to single bonds.
    /// Query bonds are not allowed.
    pub fn from_table_bond(bond: &TableBond) -> Result<Bond, ResolutionError> {
        let sigma_order = match bond.order {
            BondOrder::Aromatic => 1,
            o if o.is_query() => return Err(ResolutionError::InvalidBondOrder(bond.order)),
            o => o
                .value()
                .ok_or(ResolutionError::InvalidBondOrder(bond.order))?,
        };
        Ok(match bond.donation {
            Some(d) => Bond::new_dative(sigma_order, d),
            None => Bond::new(sigma_order),
        })
    }
}
