//! Bond types for GraphIR.

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

    pub fn to_builder(&self) -> BondBuilder {
        BondBuilder::new(self.order, None)
    }
}

/// Mutable bond representation used during resolution phases.
/// Carries an aromaticity hint that is consumed by Kekulization;
/// `build()` produces the final `Bond` with a definite order.
#[derive(Debug, Clone, PartialEq)]
pub struct BondBuilder {
    order: u8,
    aromatic_hint: Option<bool>,
}

impl BondBuilder {
    pub fn new(order: u8, aromatic_hint: Option<bool>) -> Self {
        Self {
            order,
            aromatic_hint,
        }
    }

    pub fn order(&self) -> u8 {
        self.order
    }

    pub fn set_order(&mut self, order: u8) {
        self.order = order;
    }

    pub fn aromatic_hint(&self) -> Option<bool> {
        self.aromatic_hint
    }

    pub fn set_aromatic_hint(&mut self, aromatic: Option<bool>) {
        self.aromatic_hint = aromatic;
    }

    pub fn from_table_bond(bond: &TableBond) -> Result<Self, ResolutionError> {
        match bond.order {
            BondOrder::Aromatic => Ok(Self::new(1, Some(true))),
            o if o.is_query() => Err(ResolutionError::InvalidBondOrder(bond.order)),
            o => {
                let order = o
                    .value()
                    .ok_or(ResolutionError::InvalidBondOrder(bond.order))?;
                Ok(Self::new(order, Some(false)))
            }
        }
    }

    pub fn build(&self) -> Bond {
        Bond::new(self.order)
    }
}
