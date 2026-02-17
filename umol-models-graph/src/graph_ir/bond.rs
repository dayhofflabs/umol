//! Bond type for GraphIR.

use crate::bond::{BondDonation};

/// Basic bond IR
#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    order: u8,
    donation: Option<BondDonation>,
}

impl Bond {
    pub fn order(&self) -> u8 {
        self.order
    }

    pub fn donation(&self) -> Option<BondDonation> {
        self.donation
    }
}
