//! Bond spec
//!
//! Defines bond specs for strictly typed molecular valence graphs and an internal
//! string notation, primarily intended for easy definition within data files or
//! code, not for general exchange.
//!
//! ## Internal Notation Format for Bond Specs
//!
//! The notation resembles SMARTS bond primitives with an additional notation for
//! bond donation (for dative bonds).
//!
//! The format is: `[BondOrder][Donation]`
//!
//! - `Bond order`: ".": zero, "-": single, "=": double, "#": triple, "$": quadruple
//! - `Donation`: "|": shared, ">": donating, "<": accepting
//!
//! If bond donation is not specified, it is assumed to be shared.
//!
//! ### Examples
//! - `.` -> zero (no) bond
//! - `-` -> single bond
//! - `=` -> double bond
//! - `#` -> triple bond
//! - `$` -> quadruple bond
//! - `->` -> single bond with donating atom
//! - `-<` -> single bond with accepting atom

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
    Zero,
    Single,
    Double,
    Triple,
    Quadruple,
}

pub const MAX_BOND_ORDER: u8 = 4;

// Static data for bond orders
static BOND_DATA: Lazy<HashMap<BondOrder, (u8, &'static [&'static str; 5])>> = Lazy::new(|| {
    hash_map! {
        BondOrder::Zero => (0, &[".", "·", "0", "N", "none"]),
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
            BondOrder::Zero => Some(BondOrder::Single),
            BondOrder::Single => Some(BondOrder::Double),
            BondOrder::Double => Some(BondOrder::Triple),
            BondOrder::Triple => Some(BondOrder::Quadruple),
            BondOrder::Quadruple => None,
        }
    }

    pub fn decrease(self) -> Option<Self> {
        match self {
            BondOrder::Zero => None,
            BondOrder::Single => Some(BondOrder::Zero),
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
            BondOrder::Zero => write!(f, "."),
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
    pub fn from_value(value: i8) -> Option<Self> {
        match value {
            0 => Some(BondDonation::Shared),
            -1 => Some(BondDonation::Donating),
            1 => Some(BondDonation::Accepting),
            _ => None,
        }
    }
    pub fn from_symbol(symbol: &str) -> Option<Self> {
        match symbol {
            "|" => Some(BondDonation::Shared),
            ">" => Some(BondDonation::Donating),
            "<" => Some(BondDonation::Accepting),
            _ => None,
        }
    }

    pub fn value(&self) -> i8 {
        match self {
            BondDonation::Shared => 0,
            BondDonation::Donating => -1,
            BondDonation::Accepting => 1,
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
pub struct BondSpec {
    order: BondOrder,
    donation: BondDonation,
}

impl BondSpec {
    pub fn new(order: BondOrder, donation: BondDonation) -> Self {
        Self { order, donation }
    }

    pub fn from_symbol(symbol: &str) -> Option<Self> {
        if symbol.is_empty() {
            return None;
        }

        let bond_pattern = Regex::new(r"^([-.=:#$])([<>|])?$").unwrap();
        if !bond_pattern.is_match(symbol) {
            return None;
        }

        let caps = bond_pattern.captures(symbol).unwrap();
        let order = caps[1].parse::<BondOrder>().unwrap();
        let donation = caps
            .get(2)
            .map(|m| BondDonation::from_symbol(m.as_str()).unwrap());
        Some(BondSpec::new(
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
        self.order.increase().map(|order| Self::new(order, self.donation))
    }

    pub fn decrease(self) -> Option<Self> {
        self.order.decrease().map(|order| Self::new(order, self.donation))
    }

    pub fn reverse(self) -> Option<Self> {
        self.donation.reverse().map(|donation| Self::new(self.order, donation))
    }
}

impl Display for BondSpec {
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

impl TryFrom<&str> for BondSpec {
    type Error = Error;

    fn try_from(s: &str) -> Result<Self> {
        Self::from_symbol(s).ok_or_else(|| DataError::InvalidBondSpec(s.to_string()).into())
    }
}

impl FromStr for BondSpec {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        Self::try_from(s)
    }
}

impl From<BondOrder> for BondSpec {
    fn from(order: BondOrder) -> Self {
        Self {
            order,
            donation: BondDonation::Shared,
        }
    }
}

impl From<BondDonation> for BondSpec {
    fn from(donation: BondDonation) -> Self {
        Self {
            order: BondOrder::Single,
            donation,
        }
    }
}

impl From<BondSpec> for BondOrder {
    fn from(bond: BondSpec) -> Self {
        bond.order
    }
}

impl From<BondSpec> for BondDonation {
    fn from(bond: BondSpec) -> Self {
        bond.donation
    }
}

/// Shorthand macro for bond spec parsing
#[macro_export]
macro_rules! b {
    ($s:expr) => {
        $s.parse::<BondSpec>().unwrap()
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;
    use serde_json;

    #[rstest]
    #[case(BondOrder::Zero, "none")]
    #[case(BondOrder::Single, "single")]
    #[case(BondOrder::Double, "double")]
    #[case(BondOrder::Triple, "triple")]
    #[case(BondOrder::Quadruple, "quadruple")]
    fn test_bond_order_from_symbol(#[case] bond_order: BondOrder, #[case] symbol: &str) {
        assert_eq!(BondOrder::from_symbol(symbol).unwrap(), bond_order);
    }

    #[rstest]
    #[case(BondOrder::Zero, 0)]
    #[case(BondOrder::Single, 1)]
    #[case(BondOrder::Double, 2)]
    #[case(BondOrder::Triple, 3)]
    #[case(BondOrder::Quadruple, 4)]
    fn test_bond_order_from_value(#[case] bond_order: BondOrder, #[case] value: u8) {
        assert_eq!(BondOrder::from_value(value).unwrap(), bond_order);
    }

    #[rstest]
    #[case(BondOrder::Zero, "none")]
    #[case(BondOrder::Single, "single")]
    #[case(BondOrder::Double, "double")]
    #[case(BondOrder::Triple, "triple")]
    #[case(BondOrder::Quadruple, "quadruple")]
    fn test_bond_order_from_str(#[case] bond_order: BondOrder, #[case] symbol: &str) {
        assert_eq!(BondOrder::from_str(symbol).unwrap(), bond_order);
    }

    #[rstest]
    #[case(BondOrder::Zero, ".")]
    #[case(BondOrder::Single, "-")]
    #[case(BondOrder::Double, "=")]
    #[case(BondOrder::Triple, "#")]
    #[case(BondOrder::Quadruple, "$")]
    fn test_bond_order_display(#[case] bond_order: BondOrder, #[case] symbol: &str) {
        assert_eq!(bond_order.to_string(), symbol);
    }

    #[rstest]
    #[case(BondOrder::Zero, ".")]
    #[case(BondOrder::Single, "-")]
    #[case(BondOrder::Double, "=")]
    #[case(BondOrder::Triple, "#")]
    #[case(BondOrder::Quadruple, "$")]
    fn test_bond_order_symbol(#[case] bond_order: BondOrder, #[case] symbol: &str) {
        assert_eq!(bond_order.symbol(), symbol);
    }

    #[rstest]
    #[case(BondOrder::Zero, 0)]
    #[case(BondOrder::Single, 1)]
    #[case(BondOrder::Double, 2)]
    #[case(BondOrder::Triple, 3)]
    #[case(BondOrder::Quadruple, 4)]
    fn test_bond_order_value(#[case] bond_order: BondOrder, #[case] value: u8) {
        assert_eq!(bond_order.value(), value);
    }

    #[test]
    fn test_bond_order_ordering() {
        assert!(BondOrder::Zero < BondOrder::Single);
        assert!(BondOrder::Single < BondOrder::Double);
        assert!(BondOrder::Double < BondOrder::Triple);
        assert!(BondOrder::Triple < BondOrder::Quadruple);
    }

    #[test]
    fn test_bond_order_serialization() {
        let bonds = vec![
            BondOrder::Zero,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Triple,
            BondOrder::Quadruple,
        ];
        let serialized = serde_json::to_string(&bonds).unwrap();
        assert_eq!(
            serialized,
            r#"["Zero","Single","Double","Triple","Quadruple"]"#
        );
    }

    #[test]
    fn test_bond_order_deserialization() {
        let serialized = r#"["Zero","Single","Double","Triple","Quadruple"]"#;
        let bonds: Vec<BondOrder> = serde_json::from_str(serialized).unwrap();
        assert_eq!(
            bonds,
            vec![
                BondOrder::Zero,
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Triple,
                BondOrder::Quadruple,
            ]
        );
    }

    #[rstest]
    #[case(BondOrder::Zero, Some(BondOrder::Single))]
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
    #[case(BondOrder::Zero, None)]
    #[case(BondOrder::Single, Some(BondOrder::Zero))]
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
    #[case(BondDonation::Donating, -1)]
    #[case(BondDonation::Accepting, 1)]
    #[case(BondDonation::Shared, 0)]
    fn test_bond_donation_value(#[case] donation: BondDonation, #[case] value: i8) {
        assert_eq!(donation.value(), value);
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
    #[case(".", BondSpec { order: BondOrder::Zero, donation: BondDonation::Shared })]
    #[case("-", BondSpec { order: BondOrder::Single, donation: BondDonation::Shared })]
    #[case("-|", BondSpec { order: BondOrder::Single, donation: BondDonation::Shared })]
    #[case("->", BondSpec { order: BondOrder::Single, donation: BondDonation::Donating })]
    #[case("-<", BondSpec { order: BondOrder::Single, donation: BondDonation::Accepting })]
    #[case("=", BondSpec { order: BondOrder::Double, donation: BondDonation::Shared })]
    #[case("=|", BondSpec { order: BondOrder::Double, donation: BondDonation::Shared })]
    #[case("=>", BondSpec { order: BondOrder::Double, donation: BondDonation::Donating })]
    #[case("=<", BondSpec { order: BondOrder::Double, donation: BondDonation::Accepting })]
    #[case("#", BondSpec { order: BondOrder::Triple, donation: BondDonation::Shared })]
    #[case("#|", BondSpec { order: BondOrder::Triple, donation: BondDonation::Shared })]
    #[case("#>", BondSpec { order: BondOrder::Triple, donation: BondDonation::Donating })]
    #[case("#<", BondSpec { order: BondOrder::Triple, donation: BondDonation::Accepting })]
    #[case("$", BondSpec { order: BondOrder::Quadruple, donation: BondDonation::Shared })]
    #[case("$|", BondSpec { order: BondOrder::Quadruple, donation: BondDonation::Shared })]
    #[case("$>", BondSpec { order: BondOrder::Quadruple, donation: BondDonation::Donating })]
    #[case("$<", BondSpec { order: BondOrder::Quadruple, donation: BondDonation::Accepting })]
    fn test_bond_spec_from_str(#[case] s: &str, #[case] bond: BondSpec) {
        assert_eq!(BondSpec::from_str(s).unwrap(), bond);
    }

    #[test]
    fn test_bond_spec_from_str_invalid() {
        assert!(BondSpec::from_str("invalid").is_err());
    }

    #[rstest]
    #[case(BondSpec::new(BondOrder::Zero, BondDonation::Shared), ".")]
    #[case(BondSpec::new(BondOrder::Single, BondDonation::Shared), "-")]
    #[case(BondSpec::new(BondOrder::Single, BondDonation::Donating), "->")]
    #[case(BondSpec::new(BondOrder::Single, BondDonation::Accepting), "-<")]
    #[case(BondSpec::new(BondOrder::Double, BondDonation::Shared), "=")]
    #[case(BondSpec::new(BondOrder::Double, BondDonation::Donating), "=>")]
    #[case(BondSpec::new(BondOrder::Double, BondDonation::Accepting), "=<")]
    #[case(BondSpec::new(BondOrder::Triple, BondDonation::Shared), "#")]
    #[case(BondSpec::new(BondOrder::Triple, BondDonation::Donating), "#>")]
    #[case(BondSpec::new(BondOrder::Triple, BondDonation::Accepting), "#<")]
    #[case(BondSpec::new(BondOrder::Quadruple, BondDonation::Shared), "$")]
    #[case(BondSpec::new(BondOrder::Quadruple, BondDonation::Donating), "$>")]
    #[case(BondSpec::new(BondOrder::Quadruple, BondDonation::Accepting), "$<")]
    fn test_bond_spec_display(#[case] bond: BondSpec, #[case] symbol: &str) {
        assert_eq!(bond.to_string(), symbol);
    }

    #[test]
    fn test_bond_spec_serialization() {
        let bonds = vec![
            BondSpec::new(BondOrder::Zero, BondDonation::Shared),
            BondSpec::new(BondOrder::Single, BondDonation::Shared),
            BondSpec::new(BondOrder::Double, BondDonation::Donating),
            BondSpec::new(BondOrder::Triple, BondDonation::Accepting),
            BondSpec::new(BondOrder::Quadruple, BondDonation::Shared),
        ];
        let serialized = serde_json::to_string(&bonds).unwrap();
        let expected = concat!(
            r#"[{"order":"Zero","donation":"Shared"},{"order":"Single","donation":"Shared"},"#,
            r#"{"order":"Double","donation":"Donating"},{"order":"Triple","donation":"Accepting"},"#,
            r#"{"order":"Quadruple","donation":"Shared"}]"#
        );
        assert_eq!(serialized, expected);
    }

    #[test]
    fn test_bond_spec_deserialization() {
        let serialized = concat!(
            r#"[{"order":"Zero","donation":"Shared"},{"order":"Single","donation":"Shared"},"#,
            r#"{"order":"Double","donation":"Donating"},{"order":"Triple","donation":"Accepting"},"#,
            r#"{"order":"Quadruple","donation":"Shared"}]"#
        );
        let bonds: Vec<BondSpec> = serde_json::from_str(serialized).unwrap();
        assert_eq!(
            bonds,
            vec![
                BondSpec::new(BondOrder::Zero, BondDonation::Shared),
                BondSpec::new(BondOrder::Single, BondDonation::Shared),
                BondSpec::new(BondOrder::Double, BondDonation::Donating),
                BondSpec::new(BondOrder::Triple, BondDonation::Accepting),
                BondSpec::new(BondOrder::Quadruple, BondDonation::Shared),
            ]
        );
    }

    #[rstest]
    #[case(
        BondSpec::new(BondOrder::Single, BondDonation::Shared),
        Some(BondSpec::new(BondOrder::Double, BondDonation::Shared))
    )]
    fn test_bond_spec_increase(#[case] bond: BondSpec, #[case] expected: Option<BondSpec>) {
        assert_eq!(bond.increase(), expected);
    }

    #[rstest]
    #[case(
        BondSpec::new(BondOrder::Double, BondDonation::Shared),
        Some(BondSpec::new(BondOrder::Single, BondDonation::Shared))
    )]
    fn test_bond_spec_decrease(#[case] bond: BondSpec, #[case] expected: Option<BondSpec>) {
        assert_eq!(bond.decrease(), expected);
    }

    #[rstest]
    #[case(BondSpec::new(BondOrder::Single, BondDonation::Shared), None)]
    #[case(
        BondSpec::new(BondOrder::Single, BondDonation::Donating),
        Some(BondSpec::new(BondOrder::Single, BondDonation::Accepting))
    )]
    #[case(
        BondSpec::new(BondOrder::Single, BondDonation::Accepting),
        Some(BondSpec::new(BondOrder::Single, BondDonation::Donating))
    )]
    fn test_bond_spec_reverse(#[case] bond: BondSpec, #[case] expected: Option<BondSpec>) {
        assert_eq!(bond.reverse(), expected);
    }

    #[rstest]
    #[case(b!("."), BondSpec::new(BondOrder::Zero, BondDonation::Shared))]
    #[case(b!("-"), BondSpec::new(BondOrder::Single, BondDonation::Shared))]
    #[case(b!("="), BondSpec::new(BondOrder::Double, BondDonation::Shared))]
    #[case(b!("#"), BondSpec::new(BondOrder::Triple, BondDonation::Shared))]
    #[case(b!("$"), BondSpec::new(BondOrder::Quadruple, BondDonation::Shared))]
    #[case(b!("->"), BondSpec::new(BondOrder::Single, BondDonation::Donating))]
    #[case(b!("-<"), BondSpec::new(BondOrder::Single, BondDonation::Accepting))]
    #[case(b!("=>"), BondSpec::new(BondOrder::Double, BondDonation::Donating))]
    #[case(b!("=<"), BondSpec::new(BondOrder::Double, BondDonation::Accepting))]
    #[case(b!("#<"), BondSpec::new(BondOrder::Triple, BondDonation::Accepting))]
    #[case(b!("#>"), BondSpec::new(BondOrder::Triple, BondDonation::Donating))]
    #[case(b!("$<"), BondSpec::new(BondOrder::Quadruple, BondDonation::Accepting))]
    #[case(b!("$>"), BondSpec::new(BondOrder::Quadruple, BondDonation::Donating))]
    fn test_bond_spec_macro(#[case] bond: BondSpec, #[case] expected: BondSpec) {
        assert_eq!(bond, expected);
    }
}
