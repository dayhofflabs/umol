//! Bond type for GraphIR.

use crate::bond::BondDonation;

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
}
