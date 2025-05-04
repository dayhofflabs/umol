//! Valence bond type and builder
//!
//! Valence bond is the edge type of valence graphs and is defined by its bond order and bond donation.
//! It should be created using the `ValenceBondBuilder`.

use crate::{BondDonation, BondMatcher, BondOrder, BondSpec, DEFAULT_BOND_MATCHER};
use serde::{Deserialize, Serialize};
use std::fmt::{self, Display};
use umol::{error::DataError, Result};

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
            order: self.order,
            donation: Some(self.donation),
        }
    }

    pub fn to_spec(self) -> BondSpec {
        BondSpec::new(self.order, self.donation)
    }
}

impl Display for Bond {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_spec())
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
    order: BondOrder,
    donation: Option<BondDonation>,
}

impl BondBuilder {
    pub fn new(order: BondOrder) -> Self {
        Self {
            order,
            donation: None,
        }
    }

    pub fn from_spec(bond_spec: BondSpec) -> Self {
        Self {
            order: bond_spec.order(),
            donation: Some(bond_spec.donation()),
        }
    }

    pub fn order(&self) -> BondOrder {
        self.order
    }

    pub fn donation(&self) -> Option<BondDonation> {
        self.donation
    }

    pub fn set_order(&mut self, order: BondOrder) -> &mut Self {
        self.order = order;
        self
    }

    pub fn set_donation(&mut self, donation: BondDonation) -> &mut Self {
        self.donation = Some(donation);
        self
    }

    pub fn update_order(&mut self, f: impl FnOnce(BondOrder) -> BondOrder) -> &mut Self {
        self.order = f(self.order);
        self
    }

    pub fn update_donation(&mut self, f: impl FnOnce(BondDonation) -> BondDonation) -> &mut Self {
        self.donation = Some(f(self.donation.unwrap_or(BondDonation::Shared)));
        self
    }

    pub fn build(self) -> Result<Bond> {
        self.build_with(&DEFAULT_BOND_MATCHER)
    }

    pub fn build_with(self, matcher: &BondMatcher) -> Result<Bond> {
        let bond_specs = matcher.find(&self)?;
        if bond_specs.is_empty() {
            return Err(DataError::NoBondSpec(format!("{:?}", self)).into());
        } else if bond_specs.len() > 1 {
            return Err(DataError::MultipleBondSpecs(format!("{:?}", self)).into());
        }
        let bond_spec = bond_specs.first().unwrap();
        Ok(Bond {
            order: bond_spec.order(),
            donation: bond_spec.donation(),
        })
    }
}

impl From<BondSpec> for BondBuilder {
    fn from(bond_spec: BondSpec) -> Self {
        BondBuilder::from_spec(bond_spec)
    }
}

impl From<BondOrder> for BondBuilder {
    fn from(order: BondOrder) -> Self {
        BondBuilder::new(order)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{b, ALWAYS_BOND_MATCHER};

    #[test]
    fn test_bond_display() {
        let mut builder = BondBuilder::new(BondOrder::Single);
        builder.set_donation(BondDonation::Shared);
        let bond = builder.build().unwrap();
        assert_eq!(format!("{}", bond), "-");
    }

    #[test]
    fn test_bond_serialize() {
        let mut builder = BondBuilder::new(BondOrder::Single);
        builder.set_donation(BondDonation::Shared);
        let bond = builder.build().unwrap();
        let serialized = serde_json::to_string(&bond).unwrap();
        assert_eq!(serialized, "{\"order\":\"Single\",\"donation\":\"Shared\"}");
    }

    #[test]
    fn test_bond_to_builder() {
        let mut builder = BondBuilder::new(BondOrder::Single);
        builder.set_donation(BondDonation::Shared);
        let bond = builder.build().unwrap();
        let builder = bond.to_builder();
        assert_eq!(builder.order(), BondOrder::Single);
        assert_eq!(builder.donation(), Some(BondDonation::Shared));
    }

    #[test]
    fn test_bond_builder_new() {
        let builder = BondBuilder::new(BondOrder::Single);
        assert_eq!(builder.order(), BondOrder::Single);
        assert_eq!(builder.donation(), None);
    }

    #[test]
    fn test_bond_builder_from_spec() {
        let builder = BondBuilder::from_spec(b!("->"));
        assert_eq!(builder.order(), BondOrder::Single);
        assert_eq!(builder.donation(), Some(BondDonation::Donating));
    }

    #[test]
    fn test_bond_builder_properties() {
        let builder = BondBuilder::new(BondOrder::Single);
        assert_eq!(builder.order(), BondOrder::Single);
        assert_eq!(builder.donation(), None);
    }

    #[test]
    fn test_bond_builder_set() {
        let mut builder = BondBuilder::new(BondOrder::Single);
        builder.set_order(BondOrder::Double);
        builder.set_donation(BondDonation::Donating);

        assert_eq!(builder.order(), BondOrder::Double);
        assert_eq!(builder.donation(), Some(BondDonation::Donating));
    }

    #[test]
    fn test_bond_builder_update() {
        let mut builder = BondBuilder::new(BondOrder::Single);
        builder
            .set_donation(BondDonation::Accepting)
            .update_order(|x| x.increase().unwrap())
            .update_donation(|x| x.reverse().unwrap());

        assert_eq!(builder.order(), BondOrder::Double);
        assert_eq!(builder.donation(), Some(BondDonation::Donating));
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
    fn test_bond_builder_build_with() {
        let builder = BondBuilder::new(BondOrder::Quadruple);

        let bond = builder.build_with(&ALWAYS_BOND_MATCHER).unwrap();
        assert_eq!(bond.order(), BondOrder::Quadruple);
        assert_eq!(bond.donation(), BondDonation::Shared);
    }

    #[test]
    fn test_bond_spec_into_bond_builder() {
        let bond_spec = BondSpec::new(BondOrder::Single, BondDonation::Shared);
        let builder: BondBuilder = bond_spec.into();
        assert_eq!(builder.order(), BondOrder::Single);
        assert_eq!(builder.donation(), Some(BondDonation::Shared));
    }

    #[test]
    fn test_bond_into_bond_builder() {
        let mut builder = BondBuilder::new(BondOrder::Single);
        builder.set_donation(BondDonation::Accepting);

        let bond = builder.build().unwrap();
        let builder: BondBuilder = bond.into();
        assert_eq!(builder.order(), BondOrder::Single);
        assert_eq!(builder.donation(), Some(BondDonation::Accepting));
    }

    #[test]
    fn test_bond_order_into_bond_builder() {
        let builder: BondBuilder = BondOrder::Single.into();
        assert_eq!(builder.order(), BondOrder::Single);
        assert_eq!(builder.donation(), None);
    }
}
