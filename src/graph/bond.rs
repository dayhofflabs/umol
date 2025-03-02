// GraphBond implementation

use super::atom::GraphAtom;
use super::types::AtomIndex;
use crate::error::MoleculeError;
use crate::link::AtomLink;
use std::fmt::{self, Display};

#[derive(Debug, Clone, Copy)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Quadruple,
}

impl TryFrom<&str> for BondOrder {
    type Error = MoleculeError;

    fn try_from(symbol: &str) -> Result<Self, Self::Error> {
        match symbol {
            "-" | "–" | "1" | "S" | "single" => Ok(BondOrder::Single),
            "=" | "⹀" | "2" | "D" | "double" => Ok(BondOrder::Double),
            "#" | "≡" | "3" | "T" | "triple" => Ok(BondOrder::Triple),
            "$" | "4" | "Q" | "quadruple" => Ok(BondOrder::Quadruple),
            _ => Err(MoleculeError::InvalidBondOrder(symbol.to_string())),
        }
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

#[derive(Debug, Clone)]
pub struct GraphBond {
    order: BondOrder,
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
        match BondOrder::try_from(symbol) {
            Ok(order) => Ok(order.into()),
            Err(e) => Err(e),
        }
    }
}

