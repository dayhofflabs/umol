//! Valence bond implementation
//!
//! Valence bond is the edge type of valence graphs and is
//! defined by its bond order and bond donation.

use serde::{Deserialize, Serialize};
use std::fmt::{self, Display};
use umol::Result;
use umol_data::{BondDonation, BondOrder, CovalentBond};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ValenceBond {
    order: BondOrder,
    donation: BondDonation,
}

impl ValenceBond {
    pub fn order(&self) -> BondOrder {
        self.order
    }

    pub fn donation(&self) -> BondDonation {
        self.donation
    }
}

impl Display for ValenceBond {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        CovalentBond::new(self.order, self.donation).fmt(f)
    }
}

#[derive(Debug)]
pub struct ValenceBondBuilder {
    order: BondOrder,
    donation: BondDonation,
}

impl ValenceBondBuilder {
    pub fn new(order: BondOrder) -> Self {
        Self {
            order,
            donation: BondDonation::Shared,
        }
    }

    pub fn order(&mut self, order: BondOrder) -> &mut Self {
        self.order = order;
        self
    }

    pub fn donation(&mut self, donation: BondDonation) -> &mut Self {
        self.donation = donation;
        self
    }

    pub fn build(self) -> Result<ValenceBond> {
        self.validate()?;
        self.infer_type()
    }

    fn validate(&self) -> Result<()> {
        Ok(())
    }

    fn infer_type(self) -> Result<ValenceBond> {
        Ok(ValenceBond {
            order: self.order,
            donation: self.donation,
        })
    }
}

impl From<CovalentBond> for ValenceBondBuilder {
    fn from(bond: CovalentBond) -> Self {
        ValenceBondBuilder {
            order: bond.order(),
            donation: bond.donation(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_valence_bond_display() {
        let bond = ValenceBondBuilder::new(BondOrder::Single).build().unwrap();
        assert_eq!(bond.to_string(), "-");
    }

    #[test]
    fn test_valence_bond_builder() {
        let bond = ValenceBondBuilder::new(BondOrder::Single).build().unwrap();
        assert_eq!(bond.order(), BondOrder::Single);
        assert_eq!(bond.donation(), BondDonation::Shared);

        let mut builder = ValenceBondBuilder::new(BondOrder::Single);
        builder
            .order(BondOrder::Double)
            .donation(BondDonation::Donating);
        let bond = builder.build().unwrap();
        assert_eq!(bond.order(), BondOrder::Double);
        assert_eq!(bond.donation(), BondDonation::Donating);
    }

    #[test]
    fn test_valence_bond_builder_validation() {
        let builder = ValenceBondBuilder::new(BondOrder::Quadruple);
        println!("{:?}", builder);
        let result = builder.build();
        println!("{:?}", result);
        assert!(result.is_ok());
    }
}
