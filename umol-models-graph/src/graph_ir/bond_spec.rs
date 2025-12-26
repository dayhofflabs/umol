//! Bond spec definitions for GraphIR.

use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;

use map_macro::hash_map;
use once_cell::sync::Lazy;
use regex::Regex;
use serde::{Deserialize, Serialize};
use crate::graph_ir::error::GraphError;

type Result<T> = std::result::Result<T, GraphError>;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, PartialOrd, Ord)]
pub enum BondOrder {
    Zero,
    Single,
    Double,
    Triple,
    Quadruple,
}

pub const MAX_BOND_ORDER: u8 = 4;

static BOND_DATA: Lazy<HashMap<BondOrder, (u8, &'static [&'static str; 5])>> = Lazy::new(|| {
    hash_map! {
        BondOrder::Zero => (0, &[".", "·", "0", "N", "none"]),
        BondOrder::Single => (1, &["-", "–", "1", "S", "single"]),
        BondOrder::Double => (2, &["=", "⹀", "2", "D", "double"]),
        BondOrder::Triple => (3, &["#", "≡", "3", "T", "triple"]),
        BondOrder::Quadruple => (4, &["$", "⩸", "4", "Q", "quadruple"]),
    }
});

static SYMBOL_TO_BOND: Lazy<HashMap<&'static str, BondOrder>> = Lazy::new(|| {
    BOND_DATA
        .iter()
        .flat_map(|(order, (_, symbols))| symbols.iter().map(|symbol| (*symbol, *order)))
        .collect()
});

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
    type Error = GraphError;

    fn try_from(symbol: &str) -> Result<Self> {
        Self::from_symbol(symbol)
            .ok_or_else(|| GraphError::InvalidBondSpec(format!("Invalid bond order symbol: {}", symbol)))
    }
}

impl TryFrom<u8> for BondOrder {
    type Error = GraphError;

    fn try_from(value: u8) -> Result<Self> {
        Self::from_value(value)
            .ok_or_else(|| GraphError::InvalidBondSpec(format!("Invalid bond order value: {}", value)))
    }
}

impl FromStr for BondOrder {
    type Err = GraphError;

    fn from_str(s: &str) -> Result<Self> {
        Self::try_from(s)
    }
}

impl fmt::Display for BondOrder {
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
    Donating,
    Accepting,
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

impl fmt::Display for BondDonation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BondDonation::Shared => write!(f, "|"),
            BondDonation::Donating => write!(f, ">"),
            BondDonation::Accepting => write!(f, "<"),
        }
    }
}

impl TryFrom<&str> for BondDonation {
    type Error = GraphError;

    fn try_from(symbol: &str) -> Result<Self> {
        Self::from_symbol(symbol)
            .ok_or_else(|| GraphError::InvalidBondSpec(format!("Invalid bond donation symbol: {}", symbol)))
    }
}

impl FromStr for BondDonation {
    type Err = GraphError;

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
        self.order
            .increase()
            .map(|order| Self::new(order, self.donation))
    }

    pub fn decrease(self) -> Option<Self> {
        self.order
            .decrease()
            .map(|order| Self::new(order, self.donation))
    }

    pub fn reverse(self) -> Option<Self> {
        self.donation
            .reverse()
            .map(|donation| Self::new(self.order, donation))
    }
}

impl fmt::Display for BondSpec {
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
    type Error = GraphError;

    fn try_from(s: &str) -> Result<Self> {
        Self::from_symbol(s).ok_or_else(|| GraphError::InvalidBondSpec(format!("Invalid bond spec symbol: {}", s)))
    }
}

impl FromStr for BondSpec {
    type Err = GraphError;

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

#[macro_export]
macro_rules! b {
    ($s:expr) => {
        $s.parse::<$crate::graph_ir::BondSpec>().unwrap()
    };
}
