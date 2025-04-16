// Bond order data and validation

use once_cell::sync::Lazy;
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

// Static data for bond orders
static BOND_DATA: Lazy<HashMap<BondOrder, (u8, &'static [&'static str; 5])>> = Lazy::new(|| {
    let mut m = HashMap::new();
    m.insert(BondOrder::Single, (1, &["-", "–", "1", "S", "single"]));
    m.insert(BondOrder::Double, (2, &["=", "⹀", "2", "D", "double"]));
    m.insert(BondOrder::Triple, (3, &["#", "≡", "3", "T", "triple"]));
    m.insert(
        BondOrder::Quadruple,
        (4, &["$", "⩸", "4", "Q", "quadruple"]),
    );
    m
});

// Map from symbol to bond order
static SYMBOL_TO_BOND: Lazy<HashMap<&'static str, BondOrder>> = Lazy::new(|| {
    let mut m = HashMap::new();
    for (order, (_, symbols)) in BOND_DATA.iter() {
        for &symbol in *symbols {
            m.insert(symbol, *order);
        }
    }
    m
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

/// Shorthand macro for bond order access
#[macro_export]
macro_rules! b {
    ($bond:expr) => {
        BondOrder::from_symbol($bond).unwrap()
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
    fn test_bond_order_from_symbol(#[case] bond: BondOrder, #[case] symbol: &str) {
        assert_eq!(BondOrder::from_symbol(symbol).unwrap(), bond);
    }

    #[rstest]
    #[case(BondOrder::Single, 1)]
    #[case(BondOrder::Double, 2)]
    #[case(BondOrder::Triple, 3)]
    #[case(BondOrder::Quadruple, 4)]
    fn test_bond_order_from_value(#[case] bond: BondOrder, #[case] value: u8) {
        assert_eq!(BondOrder::from_value(value).unwrap(), bond);
    }

    #[rstest]
    #[case(BondOrder::Single, "single")]
    #[case(BondOrder::Double, "double")]
    #[case(BondOrder::Triple, "triple")]
    #[case(BondOrder::Quadruple, "quadruple")]
    fn test_bond_order_from_str(#[case] bond: BondOrder, #[case] symbol: &str) {
        assert_eq!(BondOrder::from_str(symbol).unwrap(), bond);
    }

    #[rstest]
    #[case(BondOrder::Single, "-")]
    #[case(BondOrder::Double, "=")]
    #[case(BondOrder::Triple, "#")]
    #[case(BondOrder::Quadruple, "$")]
    fn test_bond_order_display(#[case] bond: BondOrder, #[case] symbol: &str) {
        assert_eq!(bond.to_string(), symbol);
    }

    #[rstest]
    #[case(BondOrder::Single, "-")]
    #[case(BondOrder::Double, "=")]
    #[case(BondOrder::Triple, "#")]
    #[case(BondOrder::Quadruple, "$")]
    fn test_bond_order_symbol(#[case] bond: BondOrder, #[case] symbol: &str) {
        assert_eq!(bond.symbol(), symbol);
    }

    #[rstest]
    #[case(BondOrder::Single, 1)]
    #[case(BondOrder::Double, 2)]
    #[case(BondOrder::Triple, 3)]
    #[case(BondOrder::Quadruple, 4)]
    fn test_bond_order_value(#[case] bond: BondOrder, #[case] value: u8) {
        assert_eq!(bond.value(), value);
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

    #[test]
    fn test_bond_order_macro() {
        assert_eq!(b!("-"), BondOrder::Single);
        assert_eq!(b!("="), BondOrder::Double);
        assert_eq!(b!("#"), BondOrder::Triple);
        assert_eq!(b!("$"), BondOrder::Quadruple);
    }
}
