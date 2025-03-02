// GraphBond implementation

use super::atom::GraphAtom;
use super::types::AtomIndex;
use crate::error::MoleculeError;
use crate::link::AtomLink;
use once_cell::sync::Lazy;
use std::collections::HashMap;
use std::fmt::{self, Display};
use std::str::FromStr;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
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

    pub fn value(&self) -> u8 {
        BOND_DATA.get(self).unwrap().0
    }

    pub fn from_symbol(symbol: &str) -> Option<Self> {
        SYMBOL_TO_BOND.get(symbol).copied()
    }

    pub fn symbol(&self) -> &'static str {
        // Return the first (canonical) symbol for this bond order
        BOND_DATA.get(self).unwrap().1[0]
    }
}

impl TryFrom<&str> for BondOrder {
    type Error = MoleculeError;

    fn try_from(symbol: &str) -> Result<Self, Self::Error> {
        Self::from_symbol(symbol).ok_or_else(|| MoleculeError::InvalidBondOrder(symbol.to_string()))
    }
}

impl TryFrom<u8> for BondOrder {
    type Error = MoleculeError;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        Self::from_value(value).ok_or_else(|| MoleculeError::InvalidBondOrder(value.to_string()))
    }
}

impl FromStr for BondOrder {
    type Err = MoleculeError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct GraphBond {
    order: BondOrder,
    // Add a properties field if you want to support properties
}

impl GraphBond {
    pub fn new(order: BondOrder) -> Self {
        GraphBond { order }
    }

    pub fn single() -> Self {
        Self::new(BondOrder::Single)
    }

    pub fn double() -> Self {
        Self::new(BondOrder::Double)
    }

    pub fn triple() -> Self {
        Self::new(BondOrder::Triple)
    }

    pub fn quadruple() -> Self {
        Self::new(BondOrder::Quadruple)
    }

    pub fn order(&self) -> BondOrder {
        self.order
    }

    pub fn with_order(mut self, order: BondOrder) -> Self {
        self.order = order;
        self
    }

    // Add property methods if you want to support properties
}

impl Display for GraphBond {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.order)
    }
}

impl AtomLink<GraphAtom> for GraphBond {
    type SiteRef = AtomIndex;
}

impl From<BondOrder> for GraphBond {
    fn from(order: BondOrder) -> Self {
        GraphBond::new(order)
    }
}

impl TryFrom<&str> for GraphBond {
    type Error = MoleculeError;

    fn try_from(symbol: &str) -> Result<Self, Self::Error> {
        BondOrder::try_from(symbol).map(|order| order.into())
    }
}

// Add GraphBondBuilder if you want to support the builder pattern

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    #[test]
    fn test_bond_new() {
        let bond = GraphBond::new(BondOrder::Single);
        assert_eq!(bond.order(), BondOrder::Single);

        let bond = GraphBond::new(BondOrder::Double);
        assert_eq!(bond.order(), BondOrder::Double);

        let bond = GraphBond::new(BondOrder::Triple);
        assert_eq!(bond.order(), BondOrder::Triple);
    }

    #[test]
    fn test_bond_with_order() {
        let bond = GraphBond::new(BondOrder::Single);
        let bond = bond.with_order(BondOrder::Double);
        assert_eq!(bond.order(), BondOrder::Double);

        let bond = bond.with_order(BondOrder::Triple);
        assert_eq!(bond.order(), BondOrder::Triple);
    }

    #[test]
    fn test_bond_order_values() {
        assert_eq!(BondOrder::Single.value(), 1);
        assert_eq!(BondOrder::Double.value(), 2);
        assert_eq!(BondOrder::Triple.value(), 3);
        assert_eq!(BondOrder::Quadruple.value(), 4);
    }

    #[test]
    fn test_bond_order_from_value() {
        assert_eq!(BondOrder::from_value(1), Some(BondOrder::Single));
        assert_eq!(BondOrder::from_value(2), Some(BondOrder::Double));
        assert_eq!(BondOrder::from_value(3), Some(BondOrder::Triple));
        assert_eq!(BondOrder::from_value(4), Some(BondOrder::Quadruple));
        assert_eq!(BondOrder::from_value(5), None);
        assert_eq!(BondOrder::from_value(0), None);
    }

    #[rstest]
    #[case(BondOrder::Single, "-")]
    #[case(BondOrder::Double, "=")]
    #[case(BondOrder::Triple, "#")]
    #[case(BondOrder::Quadruple, "$")]
    fn test_bond_order_display(#[case] order: BondOrder, #[case] expected: &str) {
        assert_eq!(format!("{}", order), expected);
    }

    #[test]
    fn test_bond_display() {
        let bond = GraphBond::new(BondOrder::Single);
        assert_eq!(format!("{}", bond), "-");

        let bond = GraphBond::new(BondOrder::Double);
        assert_eq!(format!("{}", bond), "=");
    }

    #[test]
    fn test_bond_from_str() {
        assert_eq!("-".parse::<BondOrder>(), Ok(BondOrder::Single));
        assert_eq!("=".parse::<BondOrder>(), Ok(BondOrder::Double));
        assert_eq!("#".parse::<BondOrder>(), Ok(BondOrder::Triple));
        assert_eq!("$".parse::<BondOrder>(), Ok(BondOrder::Quadruple));
        assert!("x".parse::<BondOrder>().is_err());
    }

    #[test]
    fn test_bond_try_from() {
        assert_eq!(BondOrder::try_from(1), Ok(BondOrder::Single));
        assert_eq!(BondOrder::try_from(2), Ok(BondOrder::Double));
        assert_eq!(BondOrder::try_from(3), Ok(BondOrder::Triple));
        assert_eq!(BondOrder::try_from(4), Ok(BondOrder::Quadruple));
        assert!(BondOrder::try_from(5).is_err());
    }

    #[test]
    fn test_bond_from_bond_order() {
        let bond: GraphBond = BondOrder::Single.into();
        assert_eq!(bond.order(), BondOrder::Single);

        let bond: GraphBond = BondOrder::Double.into();
        assert_eq!(bond.order(), BondOrder::Double);
    }

    // Add property tests if you implement properties
}
