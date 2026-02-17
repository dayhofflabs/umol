//! Bond type for GraphIR.

use crate::bond::BondDonation;

/// Valence bond type including strict typing.
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

impl From<Bond> for BondBuilder {
    fn from(bond: Bond) -> Self {
        Self {
            order: bond.order,
            donation: bond.donation,
        }
    }
}

/// Builder type for creating and mutating `Bond` types including strict typing.
#[derive(Debug)]
pub struct BondBuilder {
    order: u8,
    donation: Option<BondDonation>,
}

impl BondBuilder {
    pub fn new(order: u8) -> Self {
        Self {
            order,
            donation: None,
        }
    }

    pub fn order(&self) -> u8 {
        self.order
    }

    pub fn donation(&self) -> Option<BondDonation> {
        self.donation
    }

    pub fn set_order(&mut self, order: u8) -> &mut Self {
        self.order = order;
        self
    }

    pub fn set_donation(&mut self, donation: Option<BondDonation>) -> &mut Self {
        self.donation = donation;
        self
    }

    pub fn update_order(&mut self, f: impl FnOnce(u8) -> u8) -> &mut Self {
        self.order = f(self.order);
        self
    }

    pub fn update_donation(&mut self, f: impl FnOnce(BondDonation) -> BondDonation) -> &mut Self {
        self.donation = Some(f(self.donation.unwrap_or(BondDonation::Shared)));
        self
    }

    //     pub fn build(self) -> Result<Bond, ResolutionError> {
    //         self.build_with(&DEFAULT_BOND_MATCHER)
    //     }

    //     pub fn build_with(self, matcher: &BondMatcher) -> Result<Bond, ResolutionError> {
    //         let bond_specs = matcher.find(&self)?;
    //         if bond_specs.is_empty() {
    //             return Err(ResolutionError::InvalidBondSpec(format!("{:?}", self)));
    //         } else if bond_specs.len() > 1 {
    //             return Err(ResolutionError::InvalidBondSpec(format!("{:?}", self)));
    //         }
    //         let bond_spec = bond_specs.first().unwrap();
    //         Ok(Bond {
    //             order: bond_spec.order(),
    //             donation: bond_spec.donation(),
    //             wedge: self.wedge,
    //             stereo: self.stereo,
    //             ring: self.ring,
    //             span: self.span,
    //         })
    //     }
}
