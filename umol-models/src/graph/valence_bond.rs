//! Valence bond implementation
//!
//! Valence bond is the edge type of valence graphs and is
//! defined by its bond order and bond donation.

use serde::{Deserialize, Serialize};
use std::fmt::{self, Display};
use std::str::FromStr;
use umol::error::DataError;
use umol::{Error, Result};
use umol_data::{BondDonation, BondOrder, CovalentBond};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ValenceBond {
    order: BondOrder,
    donation: BondDonation,
}

impl ValenceBond {
    pub fn new(order: BondOrder, donation: BondDonation) -> Self {
        ValenceBond { order, donation }
    }

    pub fn single() -> Self {
        Self::new(BondOrder::Single, BondDonation::Shared)
    }

    pub fn double() -> Self {
        Self::new(BondOrder::Double, BondDonation::Shared)
    }

    pub fn triple() -> Self {
        Self::new(BondOrder::Triple, BondDonation::Shared)
    }

    pub fn quadruple() -> Self {
        Self::new(BondOrder::Quadruple, BondDonation::Shared)
    }

    pub fn donating() -> Self {
        Self::new(BondOrder::Single, BondDonation::Donating)
    }

    pub fn accepting() -> Self {
        Self::new(BondOrder::Single, BondDonation::Accepting)
    }

    pub fn order(&self) -> BondOrder {
        self.order
    }

    pub fn donation(&self) -> BondDonation {
        self.donation
    }
}

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

    pub fn order(mut self, order: BondOrder) -> Self {
        self.order = order;
        self
    }

    pub fn donation(mut self, donation: BondDonation) -> Self {
        self.donation = donation;
        self
    }

    pub fn build(self) -> Result<ValenceBond> {
        Ok(ValenceBond {
            order: self.order,
            donation: self.donation,
        })
    }
}

impl From<BondOrder> for ValenceBond {
    fn from(order: BondOrder) -> Self {
        ValenceBond::new(order, BondDonation::Shared)
    }
}

impl From<BondDonation> for ValenceBond {
    fn from(donation: BondDonation) -> Self {
        ValenceBond::new(BondOrder::Single, donation)
    }
}

impl From<CovalentBond> for ValenceBond {
    fn from(bond: CovalentBond) -> Self {
        ValenceBond::new(bond.order(), bond.donation())
    }
}

impl TryFrom<&str> for ValenceBond {
    type Error = Error;

    fn try_from(s: &str) -> Result<Self> {
        CovalentBond::from_str(s).map_or(
            Err(DataError::InvalidValenceBond(s.to_string()).into()),
            |bond| Ok(bond.into()),
        )
    }
}

impl FromStr for ValenceBond {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        ValenceBond::try_from(s)
    }
}

impl Display for ValenceBond {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        CovalentBond::new(self.order, self.donation).fmt(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;

    #[test]
    fn test_valence_bond_new() {
        let bond = ValenceBond::single();
        assert_eq!(bond.order(), BondOrder::Single);
        assert_eq!(bond.donation(), BondDonation::Shared);

        let bond = ValenceBond::double();
        assert_eq!(bond.order(), BondOrder::Double);
        assert_eq!(bond.donation(), BondDonation::Shared);

        let bond = ValenceBond::triple();
        assert_eq!(bond.order(), BondOrder::Triple);
        assert_eq!(bond.donation(), BondDonation::Shared);

        let bond = ValenceBond::quadruple();
        assert_eq!(bond.order(), BondOrder::Quadruple);
        assert_eq!(bond.donation(), BondDonation::Shared);

        let bond = ValenceBond::donating();
        assert_eq!(bond.order(), BondOrder::Single);
        assert_eq!(bond.donation(), BondDonation::Donating);

        let bond = ValenceBond::accepting();
        assert_eq!(bond.order(), BondOrder::Single);
        assert_eq!(bond.donation(), BondDonation::Accepting);
    }

    #[test]
    fn test_valence_bond_from_bond_order() {
        let bond = ValenceBond::from(BondOrder::Single);
        assert_eq!(bond.order(), BondOrder::Single);
        assert_eq!(bond.donation(), BondDonation::Shared);
    }

    #[test]
    fn test_valence_bond_from_bond_donation() {
        let bond = ValenceBond::from(BondDonation::Donating);
        assert_eq!(bond.order(), BondOrder::Single);
        assert_eq!(bond.donation(), BondDonation::Donating);

        let bond = ValenceBond::from(BondDonation::Accepting);
        assert_eq!(bond.order(), BondOrder::Single);
        assert_eq!(bond.donation(), BondDonation::Accepting);
    }

    #[test]
    fn test_valence_bond_builder() {
        let bond = ValenceBondBuilder::new(BondOrder::Single).build().unwrap();
        assert_eq!(bond.order(), BondOrder::Single);
        assert_eq!(bond.donation(), BondDonation::Shared);

        let bond = ValenceBondBuilder::new(BondOrder::Single)
            .order(BondOrder::Double)
            .donation(BondDonation::Donating)
            .build()
            .unwrap();
        assert_eq!(bond.order(), BondOrder::Double);
        assert_eq!(bond.donation(), BondDonation::Donating);
    }

    #[test]
    fn test_valence_bond_builder_validation() {
        let result = ValenceBondBuilder::new(BondOrder::Quadruple).build();
        assert!(result.is_ok());
    }

    #[rstest]
    #[case("-", ValenceBond { order: BondOrder::Single, donation: BondDonation::Shared })]
    #[case("=", ValenceBond { order: BondOrder::Double, donation: BondDonation::Shared })]
    fn test_valence_bond_try_from_str(#[case] s: &str, #[case] bond: ValenceBond) {
        assert_eq!(ValenceBond::try_from(s).unwrap(), bond);
    }

    #[test]
    fn test_valence_bond_try_from_str_invalid() {
        assert!(ValenceBond::try_from("nobond").is_err());
    }

    #[rstest]
    #[case("-", ValenceBond { order: BondOrder::Single, donation: BondDonation::Shared })]
    #[case("-|", ValenceBond { order: BondOrder::Single, donation: BondDonation::Shared })]
    #[case("->", ValenceBond { order: BondOrder::Single, donation: BondDonation::Donating })]
    #[case("-<", ValenceBond { order: BondOrder::Single, donation: BondDonation::Accepting })]
    #[case("=", ValenceBond { order: BondOrder::Double, donation: BondDonation::Shared })]
    #[case("=|", ValenceBond { order: BondOrder::Double, donation: BondDonation::Shared })]
    #[case("=>", ValenceBond { order: BondOrder::Double, donation: BondDonation::Donating })]
    #[case("=<", ValenceBond { order: BondOrder::Double, donation: BondDonation::Accepting })]
    #[case("#", ValenceBond { order: BondOrder::Triple, donation: BondDonation::Shared })]
    #[case("#|", ValenceBond { order: BondOrder::Triple, donation: BondDonation::Shared })]
    #[case("#>", ValenceBond { order: BondOrder::Triple, donation: BondDonation::Donating })]
    #[case("#<", ValenceBond { order: BondOrder::Triple, donation: BondDonation::Accepting })]
    #[case("$", ValenceBond { order: BondOrder::Quadruple, donation: BondDonation::Shared })]
    #[case("$|", ValenceBond { order: BondOrder::Quadruple, donation: BondDonation::Shared })]
    #[case("$>", ValenceBond { order: BondOrder::Quadruple, donation: BondDonation::Donating })]
    #[case("$<", ValenceBond { order: BondOrder::Quadruple, donation: BondDonation::Accepting })]
    fn test_covalent_bond_from_str(#[case] s: &str, #[case] bond: ValenceBond) {
        assert_eq!(ValenceBond::from_str(s).unwrap(), bond);
    }

    #[rstest]
    #[case(ValenceBond { order: BondOrder::Single, donation: BondDonation::Shared }, "-")]
    #[case(ValenceBond { order: BondOrder::Single, donation: BondDonation::Donating }, "->")]
    #[case(ValenceBond { order: BondOrder::Single, donation: BondDonation::Accepting }, "-<")]
    #[case(ValenceBond { order: BondOrder::Double, donation: BondDonation::Shared }, "=")]
    #[case(ValenceBond { order: BondOrder::Double, donation: BondDonation::Donating }, "=>")]
    #[case(ValenceBond { order: BondOrder::Double, donation: BondDonation::Accepting }, "=<")]
    #[case(ValenceBond { order: BondOrder::Triple, donation: BondDonation::Shared }, "#")]
    #[case(ValenceBond { order: BondOrder::Triple, donation: BondDonation::Donating }, "#>")]
    #[case(ValenceBond { order: BondOrder::Triple, donation: BondDonation::Accepting }, "#<")]
    #[case(ValenceBond { order: BondOrder::Quadruple, donation: BondDonation::Shared }, "$")]
    #[case(ValenceBond { order: BondOrder::Quadruple, donation: BondDonation::Donating }, "$>")]
    #[case(ValenceBond { order: BondOrder::Quadruple, donation: BondDonation::Accepting }, "$<")]
    fn test_valence_bond_display(#[case] bond: ValenceBond, #[case] symbol: &str) {
        assert_eq!(bond.to_string(), symbol);
    }
}
