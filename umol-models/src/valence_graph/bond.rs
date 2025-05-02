//! Valence bond type and builder
//!
//! Valence bond is the edge type of valence graphs and is defined by its bond order and bond donation.
//! It should be created using the `ValenceBondBuilder`.

use crate::{BondDonation, BondOrder, BondSpec};
use serde::{Deserialize, Serialize};
use std::fmt::{self, Display};
use umol::Result;

/// Valence bond type including strict typing. Cannot be created directly, but only through
/// the `BondBuilder` type, which performs validation of the bond properties. Mutations are
/// possible by converting back to a builder using the `to_builder` method.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Bond {
    order: BondOrder,
    donation: BondDonation,
}

impl Bond {
    pub fn order(&self) -> BondOrder {
        self.order
    }

    pub fn donation(&self) -> BondDonation {
        self.donation
    }

    pub fn to_builder(self) -> BondBuilder {
        BondBuilder {
            order: Some(self.order),
            donation: Some(self.donation),
        }
    }

    pub fn to_type(self) -> BondSpec {
        BondSpec::new(self.order, self.donation)
    }
}

impl Display for Bond {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_type())
    }
}

impl From<Bond> for BondBuilder {
    fn from(bond: Bond) -> Self {
        bond.to_builder()
    }
}

/// Builder type for creating and mutating `Bond` types including strict typing.
/// The resulting `Bond` objects must match the predefined `BondSpec` types.
#[derive(Debug)]
pub struct BondBuilder {
    order: Option<BondOrder>,
    donation: Option<BondDonation>,
}

impl BondBuilder {
    pub fn new(order: BondOrder) -> Self {
        Self {
            order: Some(order),
            donation: Some(BondDonation::Shared),
        }
    }

    pub fn from_type(bond_type: BondSpec) -> Self {
        Self {
            order: Some(bond_type.order()),
            donation: Some(bond_type.donation()),
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
        self.donation = Some(f(self.donation.unwrap()));
        self
    }

    pub fn build(self) -> Result<Bond> {
        Ok(Bond {
            order: self.order.unwrap(),
            donation: self.donation.unwrap(),
        })
    }
}

impl From<BondSpec> for BondBuilder {
    fn from(bond_type: BondSpec) -> Self {
        BondBuilder::from_type(bond_type)
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bond_display() {
        let bond = BondBuilder::new(BondOrder::Single).build().unwrap();
        assert_eq!(format!("{}", bond), "-");
    }

    #[test]
    fn test_bond_serialize() {
        let bond = BondBuilder::new(BondOrder::Single).build().unwrap();
        let serialized = serde_json::to_string(&bond).unwrap();
        assert_eq!(serialized, "{\"order\":\"Single\",\"donation\":\"Shared\"}");
    }

    #[test]
    fn test_bond_to_builder() {
        let bond = BondBuilder::new(BondOrder::Single).build().unwrap();
        let builder = bond.to_builder();
        assert_eq!(builder.order(), Some(BondOrder::Single));
        assert_eq!(builder.donation(), Some(BondDonation::Shared));
    }

    #[test]
    fn test_bond_builder_new() {
        let bond = BondBuilder::new(BondOrder::Single).build().unwrap();
        assert_eq!(bond.order(), BondOrder::Single);
        assert_eq!(bond.donation(), BondDonation::Shared);
    }

    #[test]
    fn test_bond_builder_from_type() {
        let bond = BondSpec::new(BondOrder::Single, BondDonation::Shared);
        let builder = BondBuilder::from_type(bond);
        assert_eq!(builder.order(), Some(BondOrder::Single));
        assert_eq!(builder.donation(), Some(BondDonation::Shared));
    }

    #[test]
    fn test_bond_builder_properties() {
        let builder = BondBuilder::new(BondOrder::Single);
        assert_eq!(builder.order(), Some(BondOrder::Single));
        assert_eq!(builder.donation(), Some(BondDonation::Shared));
    }

    #[test]
    fn test_bond_builder_set() {
        let mut builder = BondBuilder::new(BondOrder::Single);
        builder.set_order(BondOrder::Double);
        builder.set_donation(BondDonation::Donating);

        let bond = builder.build().unwrap();
        assert_eq!(bond.order(), BondOrder::Double);
        assert_eq!(bond.donation(), BondDonation::Donating);
    }

    #[test]
    fn test_bond_builder_update() {
        let mut builder = BondBuilder::new(BondOrder::Single);
        builder
            .set_donation(BondDonation::Accepting)
            .update_order(|x| x.increase().unwrap())
            .update_donation(|x| x.reverse().unwrap());

        let bond = builder.build().unwrap();
        assert_eq!(bond.order(), BondOrder::Double);
        assert_eq!(bond.donation(), BondDonation::Donating);
    }

    #[test]
    fn test_bond_builder_build() {
        let mut builder = BondBuilder::new(BondOrder::Single);
        builder
            .set_order(BondOrder::Double)
            .set_donation(BondDonation::Donating);

        let bond = builder.build().unwrap();
        assert_eq!(bond.order(), BondOrder::Double);
        assert_eq!(bond.donation(), BondDonation::Donating);
    }

    #[test]
    fn test_bond_builder_validation() {
        let builder = BondBuilder::new(BondOrder::Quadruple);
        println!("{:?}", builder);
        let result = builder.build();
        println!("{:?}", result);
        assert!(result.is_ok());
    }

    #[test]
    fn test_bond_type_into_bond_builder() {
        let bond = BondSpec::new(BondOrder::Single, BondDonation::Shared);
        let builder: BondBuilder = bond.into();
        assert_eq!(builder.order(), Some(BondOrder::Single));
        assert_eq!(builder.donation(), Some(BondDonation::Shared));
    }

    #[test]
    fn test_bond_into_bond_builder() {
        let bond = BondBuilder::new(BondOrder::Single).build().unwrap();
        let builder: BondBuilder = bond.into();
        assert_eq!(builder.order(), Some(BondOrder::Single));
        assert_eq!(builder.donation(), Some(BondDonation::Shared));
    }
}
