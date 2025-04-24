//! Covalent bond data and validation
//!
//! Uses a slightly extended SMILES notation for bond representation
//! Bond order: "-": single, "=": double, "#": triple, "$": quadruple
//! Donation: "|": shared, ">": donating, "<": accepting
//! If bond donation is not specified, it is assumed to be shared

use map_macro::hash_map;
use once_cell::sync::Lazy;
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt::{self, Display};
use std::str::FromStr;
use umol::error::DataError;
use umol::{Error, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, PartialOrd, Ord)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Quadruple,
}

pub const MAX_BOND_ORDER: u8 = 4;

// Static data for bond orders
static BOND_DATA: Lazy<HashMap<BondOrder, (u8, &'static [&'static str; 5])>> = Lazy::new(|| {
    hash_map! {
        BondOrder::Single => (1, &["-", "–", "1", "S", "single"]),
        BondOrder::Double => (2, &["=", "⹀", "2", "D", "double"]),
        BondOrder::Triple => (3, &["#", "≡", "3", "T", "triple"]),
        BondOrder::Quadruple => (4, &["$", "⩸", "4", "Q", "quadruple"]),
    }
});

// Map from symbol to bond order
static SYMBOL_TO_BOND: Lazy<HashMap<&'static str, BondOrder>> = Lazy::new(|| {
    BOND_DATA
        .iter()
        .flat_map(|(order, (_, symbols))| symbols.iter().map(|symbol| (*symbol, *order)))
        .collect()
});

// Map from value to bond order
static VALUE_TO_BOND: Lazy<HashMap<u8, BondOrder>> = Lazy::new(|| {
    BOND_DATA
        .iter()
        .map(|(order, (value, _))| (*value, *order))
        .collect()
});

impl BondOrder {
    pub fn from_value(value: u8) -> Option<Self> {
        VALUE_TO_BOND.get(&value).copied()
    }

    pub fn from_symbol(symbol: &str) -> Option<Self> {
        SYMBOL_TO_BOND.get(symbol).copied()
    }

    pub fn value(&self) -> u8 {
        BOND_DATA.get(self).unwrap().0
    }

    pub fn symbol(&self) -> &'static str {
        // Return the first (canonical) symbol for this bond order
        BOND_DATA.get(self).unwrap().1[0]
    }

    pub fn increase(self) -> Option<Self> {
        match self {
            BondOrder::Single => Some(BondOrder::Double),
            BondOrder::Double => Some(BondOrder::Triple),
            BondOrder::Triple => Some(BondOrder::Quadruple),
            BondOrder::Quadruple => None,
        }
    }

    pub fn decrease(self) -> Option<Self> {
        match self {
            BondOrder::Single => None,
            BondOrder::Double => Some(BondOrder::Single),
            BondOrder::Triple => Some(BondOrder::Double),
            BondOrder::Quadruple => Some(BondOrder::Triple),
        }
    }
}

impl TryFrom<&str> for BondOrder {
    type Error = Error;

    fn try_from(symbol: &str) -> Result<Self> {
        Self::from_symbol(symbol)
            .ok_or_else(|| DataError::InvalidBondOrder(symbol.to_string()).into())
    }
}

impl TryFrom<u8> for BondOrder {
    type Error = Error;

    fn try_from(value: u8) -> Result<Self> {
        Self::from_value(value).ok_or_else(|| DataError::InvalidBondOrder(value.to_string()).into())
    }
}

impl FromStr for BondOrder {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        Self::try_from(s)
    }
}

impl Display for BondOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BondOrder::Single => write!(f, "-"),
            BondOrder::Double => write!(f, "="),
            BondOrder::Triple => write!(f, "#"),
            BondOrder::Quadruple => write!(f, "$"),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondDonation {
    Shared,
    Donating,  // atom1 (donor) -> atom2 (acceptor)
    Accepting, // atom1 (acceptor) -> atom2 (donor)
}

impl BondDonation {
    pub fn from_symbol(symbol: &str) -> Option<Self> {
        match symbol {
            "|" => Some(BondDonation::Shared),
            ">" => Some(BondDonation::Donating),
            "<" => Some(BondDonation::Accepting),
            _ => None,
        }
    }

    pub fn symbol(&self) -> &'static str {
        match self {
            BondDonation::Shared => "|",
            BondDonation::Donating => ">",
            BondDonation::Accepting => "<",
        }
    }

    pub fn reverse(self) -> Option<Self> {
        match self {
            BondDonation::Shared => None,
            BondDonation::Donating => Some(BondDonation::Accepting),
            BondDonation::Accepting => Some(BondDonation::Donating),
        }
    }
}

impl Display for BondDonation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BondDonation::Shared => write!(f, "|"),
            BondDonation::Donating => write!(f, ">"),
            BondDonation::Accepting => write!(f, "<"),
        }
    }
}

impl TryFrom<&str> for BondDonation {
    type Error = Error;

    fn try_from(symbol: &str) -> Result<Self> {
        Self::from_symbol(symbol)
            .ok_or_else(|| DataError::InvalidBondDonation(symbol.to_string()).into())
    }
}

impl FromStr for BondDonation {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        Self::try_from(s)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CovalentBond {
    order: BondOrder,
    donation: BondDonation,
}

impl CovalentBond {
    pub fn new(order: BondOrder, donation: BondDonation) -> Self {
        Self { order, donation }
    }

    pub fn from_symbol(symbol: &str) -> Option<Self> {
        if symbol.is_empty() {
            return None;
        }

        let bond_pattern = Regex::new(r"^([-=:#$])([<>|])?$").unwrap();
        if !bond_pattern.is_match(symbol) {
            return None;
        }

        let caps = bond_pattern.captures(symbol).unwrap();
        let order = caps[1].parse::<BondOrder>().unwrap();
        let donation = caps
            .get(2)
            .map(|m| BondDonation::from_symbol(m.as_str()).unwrap());
        Some(CovalentBond::new(
            order,
            donation.unwrap_or(BondDonation::Shared),
        ))
    }

    pub fn order(&self) -> BondOrder {
        self.order
    }

    pub fn donation(&self) -> BondDonation {
        self.donation
    }

    pub fn increase(self) -> Option<Self> {
        if let Some(order) = self.order.increase() {
            Some(CovalentBond::new(order, self.donation))
        } else {
            None
        }
    }

    pub fn decrease(self) -> Option<Self> {
        if let Some(order) = self.order.decrease() {
            Some(CovalentBond::new(order, self.donation))
        } else {
            None
        }
    }

    pub fn reverse(self) -> Option<Self> {
        if let Some(donation) = self.donation.reverse() {
            Some(CovalentBond::new(self.order, donation))
        } else {
            None
        }
    }
}

impl Display for CovalentBond {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}{}",
            self.order().symbol(),
            if self.donation() == BondDonation::Shared {
                ""
            } else {
                self.donation().symbol()
            }
        )
    }
}

impl TryFrom<&str> for CovalentBond {
    type Error = Error;

    fn try_from(s: &str) -> Result<Self> {
        Self::from_symbol(s).ok_or_else(|| DataError::InvalidCovalentBond(s.to_string()).into())
    }
}

impl FromStr for CovalentBond {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        Self::try_from(s)
    }
}

impl From<BondOrder> for CovalentBond {
    fn from(order: BondOrder) -> Self {
        Self {
            order,
            donation: BondDonation::Shared,
        }
    }
}

impl From<BondDonation> for CovalentBond {
    fn from(donation: BondDonation) -> Self {
        Self {
            order: BondOrder::Single,
            donation,
        }
    }
}

impl From<CovalentBond> for BondOrder {
    fn from(bond: CovalentBond) -> Self {
        bond.order
    }
}

impl From<CovalentBond> for BondDonation {
    fn from(bond: CovalentBond) -> Self {
        bond.donation
    }
}

/// Shorthand macro for covalent bond access
#[macro_export]
macro_rules! b {
    ($bond:expr) => {
        $bond.parse::<CovalentBond>().unwrap()
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;
    use serde_json;

    #[rstest]
    #[case(BondOrder::Single, "single")]
    #[case(BondOrder::Double, "double")]
    #[case(BondOrder::Triple, "triple")]
    #[case(BondOrder::Quadruple, "quadruple")]
    fn test_bond_order_from_symbol(#[case] bond_order: BondOrder, #[case] symbol: &str) {
        assert_eq!(BondOrder::from_symbol(symbol).unwrap(), bond_order);
    }

    #[rstest]
    #[case(BondOrder::Single, 1)]
    #[case(BondOrder::Double, 2)]
    #[case(BondOrder::Triple, 3)]
    #[case(BondOrder::Quadruple, 4)]
    fn test_bond_order_from_value(#[case] bond_order: BondOrder, #[case] value: u8) {
        assert_eq!(BondOrder::from_value(value).unwrap(), bond_order);
    }

    #[rstest]
    #[case(BondOrder::Single, "single")]
    #[case(BondOrder::Double, "double")]
    #[case(BondOrder::Triple, "triple")]
    #[case(BondOrder::Quadruple, "quadruple")]
    fn test_bond_order_from_str(#[case] bond_order: BondOrder, #[case] symbol: &str) {
        assert_eq!(BondOrder::from_str(symbol).unwrap(), bond_order);
    }

    #[rstest]
    #[case(BondOrder::Single, "-")]
    #[case(BondOrder::Double, "=")]
    #[case(BondOrder::Triple, "#")]
    #[case(BondOrder::Quadruple, "$")]
    fn test_bond_order_display(#[case] bond_order: BondOrder, #[case] symbol: &str) {
        assert_eq!(bond_order.to_string(), symbol);
    }

    #[rstest]
    #[case(BondOrder::Single, "-")]
    #[case(BondOrder::Double, "=")]
    #[case(BondOrder::Triple, "#")]
    #[case(BondOrder::Quadruple, "$")]
    fn test_bond_order_symbol(#[case] bond_order: BondOrder, #[case] symbol: &str) {
        assert_eq!(bond_order.symbol(), symbol);
    }

    #[rstest]
    #[case(BondOrder::Single, 1)]
    #[case(BondOrder::Double, 2)]
    #[case(BondOrder::Triple, 3)]
    #[case(BondOrder::Quadruple, 4)]
    fn test_bond_order_value(#[case] bond_order: BondOrder, #[case] value: u8) {
        assert_eq!(bond_order.value(), value);
    }

    #[test]
    fn test_bond_order_ordering() {
        assert!(BondOrder::Single < BondOrder::Double);
        assert!(BondOrder::Double < BondOrder::Triple);
        assert!(BondOrder::Triple < BondOrder::Quadruple);
    }

    #[test]
    fn test_bond_order_serialization() {
        let bonds = vec![
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Triple,
            BondOrder::Quadruple,
        ];
        let serialized = serde_json::to_string(&bonds).unwrap();
        assert_eq!(serialized, r#"["Single","Double","Triple","Quadruple"]"#);
    }

    #[test]
    fn test_bond_order_deserialization() {
        let serialized = r#"["Single","Double","Triple","Quadruple"]"#;
        let bonds: Vec<BondOrder> = serde_json::from_str(serialized).unwrap();
        assert_eq!(
            bonds,
            vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Triple,
                BondOrder::Quadruple,
            ]
        );
    }

    #[rstest]
    #[case(BondOrder::Single, Some(BondOrder::Double))]
    #[case(BondOrder::Double, Some(BondOrder::Triple))]
    #[case(BondOrder::Triple, Some(BondOrder::Quadruple))]
    #[case(BondOrder::Quadruple, None)]
    fn test_bond_order_increase(
        #[case] bond_order: BondOrder,
        #[case] expected: Option<BondOrder>,
    ) {
        assert_eq!(bond_order.increase(), expected);
    }

    #[rstest]
    #[case(BondOrder::Single, None)]
    #[case(BondOrder::Double, Some(BondOrder::Single))]
    #[case(BondOrder::Triple, Some(BondOrder::Double))]
    #[case(BondOrder::Quadruple, Some(BondOrder::Triple))]
    fn test_bond_order_decrease(
        #[case] bond_order: BondOrder,
        #[case] expected: Option<BondOrder>,
    ) {
        assert_eq!(bond_order.decrease(), expected);
    }

    #[rstest]
    #[case(BondDonation::Donating, ">")]
    #[case(BondDonation::Accepting, "<")]
    #[case(BondDonation::Shared, "|")]
    fn test_bond_donation_from_symbol(#[case] donation: BondDonation, #[case] symbol: &str) {
        assert_eq!(BondDonation::from_symbol(symbol).unwrap(), donation);
    }

    #[test]
    fn test_bond_donation_symbol() {
        assert_eq!(BondDonation::Donating.symbol(), ">");
        assert_eq!(BondDonation::Accepting.symbol(), "<");
        assert_eq!(BondDonation::Shared.symbol(), "|");
    }

    #[rstest]
    #[case(BondDonation::Donating, ">")]
    #[case(BondDonation::Accepting, "<")]
    #[case(BondDonation::Shared, "|")]
    fn test_bond_donation_display(#[case] donation: BondDonation, #[case] symbol: &str) {
        assert_eq!(donation.to_string(), symbol);
    }

    #[rstest]
    #[case(BondDonation::Donating, ">")]
    #[case(BondDonation::Accepting, "<")]
    #[case(BondDonation::Shared, "|")]
    fn test_bond_donation_from_str(#[case] donation: BondDonation, #[case] symbol: &str) {
        assert_eq!(BondDonation::from_str(symbol).unwrap(), donation);
    }

    #[test]
    fn test_bond_donation_serialization() {
        let donations = vec![
            BondDonation::Donating,
            BondDonation::Accepting,
            BondDonation::Shared,
        ];
        let serialized = serde_json::to_string(&donations).unwrap();
        assert_eq!(serialized, r#"["Donating","Accepting","Shared"]"#);
    }

    #[test]
    fn test_bond_donation_deserialization() {
        let serialized = r#"["Donating","Accepting","Shared"]"#;
        let donations: Vec<BondDonation> = serde_json::from_str(serialized).unwrap();
        assert_eq!(
            donations,
            vec![
                BondDonation::Donating,
                BondDonation::Accepting,
                BondDonation::Shared
            ]
        );
    }

    #[rstest]
    #[case(BondDonation::Donating, Some(BondDonation::Accepting))]
    #[case(BondDonation::Accepting, Some(BondDonation::Donating))]
    #[case(BondDonation::Shared, None)]
    fn test_bond_donation_reverse(
        #[case] donation: BondDonation,
        #[case] expected: Option<BondDonation>,
    ) {
        assert_eq!(donation.reverse(), expected);
    }

    #[rstest]
    #[case("-", CovalentBond { order: BondOrder::Single, donation: BondDonation::Shared })]
    #[case("-|", CovalentBond { order: BondOrder::Single, donation: BondDonation::Shared })]
    #[case("->", CovalentBond { order: BondOrder::Single, donation: BondDonation::Donating })]
    #[case("-<", CovalentBond { order: BondOrder::Single, donation: BondDonation::Accepting })]
    #[case("=", CovalentBond { order: BondOrder::Double, donation: BondDonation::Shared })]
    #[case("=|", CovalentBond { order: BondOrder::Double, donation: BondDonation::Shared })]
    #[case("=>", CovalentBond { order: BondOrder::Double, donation: BondDonation::Donating })]
    #[case("=<", CovalentBond { order: BondOrder::Double, donation: BondDonation::Accepting })]
    #[case("#", CovalentBond { order: BondOrder::Triple, donation: BondDonation::Shared })]
    #[case("#|", CovalentBond { order: BondOrder::Triple, donation: BondDonation::Shared })]
    #[case("#>", CovalentBond { order: BondOrder::Triple, donation: BondDonation::Donating })]
    #[case("#<", CovalentBond { order: BondOrder::Triple, donation: BondDonation::Accepting })]
    #[case("$", CovalentBond { order: BondOrder::Quadruple, donation: BondDonation::Shared })]
    #[case("$|", CovalentBond { order: BondOrder::Quadruple, donation: BondDonation::Shared })]
    #[case("$>", CovalentBond { order: BondOrder::Quadruple, donation: BondDonation::Donating })]
    #[case("$<", CovalentBond { order: BondOrder::Quadruple, donation: BondDonation::Accepting })]
    fn test_covalent_bond_from_str(#[case] s: &str, #[case] bond: CovalentBond) {
        assert_eq!(CovalentBond::from_str(s).unwrap(), bond);
    }

    #[test]
    fn test_covalent_bond_from_str_invalid() {
        assert!(CovalentBond::from_str("invalid").is_err());
    }

    #[rstest]
    #[case(CovalentBond { order: BondOrder::Single, donation: BondDonation::Shared }, "-")]
    #[case(CovalentBond { order: BondOrder::Single, donation: BondDonation::Donating }, "->")]
    #[case(CovalentBond { order: BondOrder::Single, donation: BondDonation::Accepting }, "-<")]
    #[case(CovalentBond { order: BondOrder::Double, donation: BondDonation::Shared }, "=")]
    #[case(CovalentBond { order: BondOrder::Double, donation: BondDonation::Donating }, "=>")]
    #[case(CovalentBond { order: BondOrder::Double, donation: BondDonation::Accepting }, "=<")]
    #[case(CovalentBond { order: BondOrder::Triple, donation: BondDonation::Shared }, "#")]
    #[case(CovalentBond { order: BondOrder::Triple, donation: BondDonation::Donating }, "#>")]
    #[case(CovalentBond { order: BondOrder::Triple, donation: BondDonation::Accepting }, "#<")]
    #[case(CovalentBond { order: BondOrder::Quadruple, donation: BondDonation::Shared }, "$")]
    #[case(CovalentBond { order: BondOrder::Quadruple, donation: BondDonation::Donating }, "$>")]
    #[case(CovalentBond { order: BondOrder::Quadruple, donation: BondDonation::Accepting }, "$<")]
    fn test_covalent_bond_display(#[case] bond: CovalentBond, #[case] symbol: &str) {
        assert_eq!(bond.to_string(), symbol);
    }

    #[test]
    fn test_covalent_bond_serialization() {
        let bonds = vec![
            CovalentBond::new(BondOrder::Single, BondDonation::Shared),
            CovalentBond::new(BondOrder::Double, BondDonation::Donating),
            CovalentBond::new(BondOrder::Triple, BondDonation::Accepting),
            CovalentBond::new(BondOrder::Quadruple, BondDonation::Shared),
        ];
        let serialized = serde_json::to_string(&bonds).unwrap();
        let expected = concat!(
            r#"[{"order":"Single","donation":"Shared"},{"order":"Double","donation":"Donating"},"#,
            r#"{"order":"Triple","donation":"Accepting"},{"order":"Quadruple","donation":"Shared"}]"#
        );
        assert_eq!(serialized, expected);
    }

    #[test]
    fn test_covalent_bond_deserialization() {
        let serialized = concat!(
            r#"[{"order":"Single","donation":"Shared"},{"order":"Double","donation":"Donating"},"#,
            r#"{"order":"Triple","donation":"Accepting"},{"order":"Quadruple","donation":"Shared"}]"#
        );
        let bonds: Vec<CovalentBond> = serde_json::from_str(serialized).unwrap();
        assert_eq!(
            bonds,
            vec![
                CovalentBond::new(BondOrder::Single, BondDonation::Shared),
                CovalentBond::new(BondOrder::Double, BondDonation::Donating),
                CovalentBond::new(BondOrder::Triple, BondDonation::Accepting),
                CovalentBond::new(BondOrder::Quadruple, BondDonation::Shared),
            ]
        );
    }

    #[rstest]
    #[case(CovalentBond { order: BondOrder::Single, donation: BondDonation::Shared },
        Some(CovalentBond { order: BondOrder::Double, donation: BondDonation::Shared }))]
    fn test_covalent_bond_increase(
        #[case] bond: CovalentBond,
        #[case] expected: Option<CovalentBond>,
    ) {
        assert_eq!(bond.increase(), expected);
    }

    #[rstest]
    #[case(CovalentBond { order: BondOrder::Double, donation: BondDonation::Shared },
        Some(CovalentBond { order: BondOrder::Single, donation: BondDonation::Shared }))]
    fn test_covalent_bond_decrease(
        #[case] bond: CovalentBond,
        #[case] expected: Option<CovalentBond>,
    ) {
        assert_eq!(bond.decrease(), expected);
    }

    #[rstest]
    #[case(CovalentBond { order: BondOrder::Single, donation: BondDonation::Shared }, None)]
    #[case(CovalentBond { order: BondOrder::Single, donation: BondDonation::Donating },
        Some(CovalentBond { order: BondOrder::Single, donation: BondDonation::Accepting }))]
    #[case(CovalentBond { order: BondOrder::Single, donation: BondDonation::Accepting },
        Some(CovalentBond { order: BondOrder::Single, donation: BondDonation::Donating }))]
    fn test_covalent_bond_reverse(
        #[case] bond: CovalentBond,
        #[case] expected: Option<CovalentBond>,
    ) {
        assert_eq!(bond.reverse(), expected);
    }

    #[rstest]
    #[case(b!("-"), CovalentBond::new(BondOrder::Single, BondDonation::Shared))]
    #[case(b!("="), CovalentBond::new(BondOrder::Double, BondDonation::Shared))]
    #[case(b!("#"), CovalentBond::new(BondOrder::Triple, BondDonation::Shared))]
    #[case(b!("$"), CovalentBond::new(BondOrder::Quadruple, BondDonation::Shared))]
    #[case(b!("->"), CovalentBond::new(BondOrder::Single, BondDonation::Donating))]
    #[case(b!("-<"), CovalentBond::new(BondOrder::Single, BondDonation::Accepting))]
    #[case(b!("=>"), CovalentBond::new(BondOrder::Double, BondDonation::Donating))]
    #[case(b!("=<"), CovalentBond::new(BondOrder::Double, BondDonation::Accepting))]
    #[case(b!("#<"), CovalentBond::new(BondOrder::Triple, BondDonation::Accepting))]
    #[case(b!("#>"), CovalentBond::new(BondOrder::Triple, BondDonation::Donating))]
    #[case(b!("$<"), CovalentBond::new(BondOrder::Quadruple, BondDonation::Accepting))]
    #[case(b!("$>"), CovalentBond::new(BondOrder::Quadruple, BondDonation::Donating))]
    fn test_covalent_bond_macro(#[case] bond: CovalentBond, #[case] expected: CovalentBond) {
        assert_eq!(bond, expected);
    }
}
