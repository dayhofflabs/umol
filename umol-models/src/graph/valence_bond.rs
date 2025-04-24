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

    pub fn to_builder(self) -> ValenceBondBuilder {
        ValenceBondBuilder::from(self)
    }
}

impl Display for ValenceBond {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        CovalentBond::new(self.order, self.donation).fmt(f)
    }
}

#[derive(Debug)]
pub struct ValenceBondBuilder {
    order: Option<BondOrder>,
    donation: Option<BondDonation>,
}

impl ValenceBondBuilder {
    pub fn new(order: BondOrder) -> Self {
        Self {
            order: Some(order),
            donation: Some(BondDonation::Shared),
        }
    }

    pub fn from_covalent_bond(bond: CovalentBond) -> Self {
        Self {
            order: Some(bond.order()),
            donation: Some(bond.donation()),
        }
    }

    pub fn order(&self) -> Option<BondOrder> {
        self.order
    }

    pub fn donation(&self) -> Option<BondDonation> {
        self.donation
    }

    pub fn set_order(&mut self, order: BondOrder) -> &mut Self {
        self.order = Some(order);
        self
    }

    pub fn set_donation(&mut self, donation: BondDonation) -> &mut Self {
        self.donation = Some(donation);
        self
    }

    pub fn update_order(&mut self, f: impl FnOnce(BondOrder) -> BondOrder) -> &mut Self {
        self.order = Some(f(self.order.unwrap()));
        self
    }

    pub fn update_donation(&mut self, f: impl FnOnce(BondDonation) -> BondDonation) -> &mut Self {
        self.donation = Some(f(self.donation.unwrap_or(BondDonation::Shared)));
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
            order: self.order.unwrap(),
            donation: self.donation.unwrap(),
        })
    }
}

impl From<CovalentBond> for ValenceBondBuilder {
    fn from(bond: CovalentBond) -> Self {
        ValenceBondBuilder {
            order: Some(bond.order()),
            donation: Some(bond.donation()),
        }
    }
}

impl From<ValenceBond> for ValenceBondBuilder {
    fn from(bond: ValenceBond) -> Self {
        ValenceBondBuilder::from_covalent_bond(CovalentBond::new(bond.order(), bond.donation()))
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
    fn test_valence_bond_serialize() {
        let bond = ValenceBondBuilder::new(BondOrder::Single).build().unwrap();
        let serialized = serde_json::to_string(&bond).unwrap();
        assert_eq!(serialized, "{\"order\":\"Single\",\"donation\":\"Shared\"}");
    }

    #[test]
    fn test_valence_bond_to_builder() {
        let bond = ValenceBondBuilder::new(BondOrder::Single).build().unwrap();
        let builder = bond.to_builder();
        assert_eq!(builder.order(), Some(BondOrder::Single));
        assert_eq!(builder.donation(), Some(BondDonation::Shared));
    }

    #[test]
    fn test_valence_bond_builder_new() {
        let bond = ValenceBondBuilder::new(BondOrder::Single).build().unwrap();
        assert_eq!(bond.order(), BondOrder::Single);
        assert_eq!(bond.donation(), BondDonation::Shared);
    }

    #[test]
    fn test_valence_bond_builder_from_covalent_bond() {
        let bond = CovalentBond::new(BondOrder::Single, BondDonation::Shared);
        let builder = ValenceBondBuilder::from_covalent_bond(bond);
        assert_eq!(builder.order(), Some(BondOrder::Single));
        assert_eq!(builder.donation(), Some(BondDonation::Shared));
    }

    #[test]
    fn test_valence_bond_builder_properties() {
        let builder = ValenceBondBuilder::new(BondOrder::Single);
        assert_eq!(builder.order(), Some(BondOrder::Single));
        assert_eq!(builder.donation(), Some(BondDonation::Shared));
    }

    #[test]
    fn test_valence_bond_builder_set() {
        let mut builder = ValenceBondBuilder::new(BondOrder::Single);
        builder.set_order(BondOrder::Double);
        builder.set_donation(BondDonation::Donating);

        let bond = builder.build().unwrap();
        assert_eq!(bond.order(), BondOrder::Double);
        assert_eq!(bond.donation(), BondDonation::Donating);
    }

    #[test]
    fn test_valence_bond_builder_update() {
        let mut builder = ValenceBondBuilder::new(BondOrder::Single);
        builder
            .set_donation(BondDonation::Accepting)
            .update_order(|x| x.increase().unwrap())
            .update_donation(|x| x.reverse().unwrap());

        let bond = builder.build().unwrap();
        assert_eq!(bond.order(), BondOrder::Double);
        assert_eq!(bond.donation(), BondDonation::Donating);
    }

    #[test]
    fn test_valence_bond_builder_build() {
        let mut builder = ValenceBondBuilder::new(BondOrder::Single);
        builder
            .set_order(BondOrder::Double)
            .set_donation(BondDonation::Donating);

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

    #[test]
    fn test_covalent_bond_into_valence_bond_builder() {
        let bond = CovalentBond::new(BondOrder::Single, BondDonation::Shared);
        let builder: ValenceBondBuilder = bond.into();
        assert_eq!(builder.order(), Some(BondOrder::Single));
        assert_eq!(builder.donation(), Some(BondDonation::Shared));
    }

    #[test]
    fn test_valence_bond_into_valence_bond_builder() {
        let bond = ValenceBondBuilder::new(BondOrder::Single).build().unwrap();
        let builder: ValenceBondBuilder = bond.into();
        assert_eq!(builder.order(), Some(BondOrder::Single));
        assert_eq!(builder.donation(), Some(BondDonation::Shared));
    }
}
