//! Bond types for GraphIR.

use std::fmt::{self, Display};
use std::str::FromStr;

use serde::de::{Deserializer, Error as SerdeError};
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};
use thiserror::Error;
use umol_data::{SpinMultiplicity, SpinState, SpinStateError};

use crate::dsl::ast::{FromAst, ToAst};
use crate::dsl::bond::{parse_bond_dsl, BondAst, BondLowerConfig};
use crate::dsl::error::LoweringError;
use crate::dsl::value::ValueAst;
use crate::graph_ir::bond_pattern::BondPattern;
use crate::graph_ir::error::{ResolutionError, ValidationError};

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum BondError {
    #[error("invalid bond order value: {0}")]
    InvalidOrder(String),
    #[error("invalid charge value: {0}")]
    InvalidCharge(String),
    #[error("invalid multiplicity value: {0}")]
    InvalidMultiplicity(String),
    #[error(transparent)]
    SpinState(#[from] SpinStateError),
    #[error("invalid bond state: {0}")]
    InvalidState(String),
}

impl From<BondError> for ResolutionError {
    fn from(value: BondError) -> Self {
        ResolutionError::InvalidBond(value.to_string())
    }
}

/// Resolved shared (covalent) bond in GraphIR. Order is the localized (σ-skeleton)
/// bond order. Dative and non-covalent bonds are stored separately.
#[derive(Debug, Clone, PartialEq)]
pub struct Bond {
    order: u8,
    charge: i8,
    spin: SpinState,
}

impl Bond {
    pub fn new(order: u8) -> Self {
        Self {
            order,
            charge: 0,
            spin: SpinState::closed_shell(),
        }
    }

    pub(crate) fn from_parts(order: u8, charge: i8, spin: SpinState) -> Self {
        Self {
            order,
            charge,
            spin,
        }
    }

    pub fn charge(&self) -> i8 {
        self.charge
    }

    pub fn spin(&self) -> SpinState {
        self.spin
    }

    pub fn multiplicity(&self) -> SpinMultiplicity {
        self.spin.multiplicity()
    }

    pub fn unpaired_electrons(&self) -> u8 {
        self.spin.unpaired_electrons()
    }

    pub fn order(&self) -> u8 {
        self.order
    }
}

impl FromAst<BondAst> for Bond {
    fn from_ast(ast: BondAst, cfg: &BondLowerConfig) -> Result<Self, LoweringError> {
        let pattern = BondPattern::from_ast(ast, cfg)?;
        pattern.to_bond().map_err(|e| match e {
            ValidationError::NonGround { field } => LoweringError::NonGround { field },
            ValidationError::InvalidMultiplicity(n) => LoweringError::InvalidMultiplicity(n),
            ValidationError::Bond(be) => match be {
                BondError::SpinState(se) => LoweringError::SpinState(se),
                other => LoweringError::Atom(other.to_string()),
            },
            other => LoweringError::Atom(other.to_string()),
        })
    }
}

impl ToAst<BondAst> for Bond {
    fn to_ast(&self) -> BondAst {
        BondAst {
            order: ValueAst::Lit(self.order() as i32),
            charge: Some(ValueAst::Lit(self.charge() as i32)),
            unpaired_electrons: Some(ValueAst::Lit(self.unpaired_electrons() as i32)),
            multiplicity: Some(ValueAst::Lit(self.multiplicity().multiplicity() as i32)),
        }
    }
}

impl FromStr for Bond {
    type Err = BondError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let ast = parse_bond_dsl(s).map_err(|e| BondError::InvalidState(e.to_string()))?;
        Bond::from_ast(ast, &BondLowerConfig::default())
            .map_err(|e| BondError::InvalidState(e.to_string()))
    }
}

impl Display for Bond {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.to_ast().fmt(f)
    }
}

impl Serialize for Bond {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for Bond {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let s = String::deserialize(deserializer)?;
        s.parse().map_err(SerdeError::custom)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_data::SpinMultiplicity;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::defaults(BondAst { order: ValueAst::Lit(1), charge: None, unpaired_electrons: None, multiplicity: None }, Bond::new(1))]
    #[case::charged(BondAst { order: ValueAst::Lit(2), charge: Some(ValueAst::Lit(1)), unpaired_electrons: Some(ValueAst::Lit(1)), multiplicity: None },
        Bond::from_parts(2, 1, SpinState::new(1, SpinMultiplicity::Doublet)))]
    #[case::full(BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(0)), unpaired_electrons: Some(ValueAst::Lit(2)), multiplicity: Some(ValueAst::Lit(1)), },
        Bond::from_parts(1, 0, SpinState::new(2, SpinMultiplicity::Singlet)))]
    fn test_bond_from_ast(#[case] ast: BondAst, #[case] expected: Bond) {
        assert_eq!(Bond::from_ast(ast, &BondLowerConfig::default()).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(Bond::new(1), BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(0)), unpaired_electrons: Some(ValueAst::Lit(0)), multiplicity: Some(ValueAst::Lit(1)) })]
    #[case::charged_doublet( Bond::from_parts(2, 1, SpinState::new(1, SpinMultiplicity::Doublet)), BondAst { order: ValueAst::Lit(2), charge: Some(ValueAst::Lit(1)),
        unpaired_electrons: Some(ValueAst::Lit(1)), multiplicity: Some(ValueAst::Lit(2)) })]
    fn test_bond_to_ast(#[case] bond: Bond, #[case] expected: BondAst) {
        assert_eq!(bond.to_ast(), expected);
    }

    #[rstest]
    #[case::single("1", 1, 0, 0, SpinMultiplicity::Singlet)]
    #[case::double_high_spin("2#c+#u1", 2, 1, 1, SpinMultiplicity::Doublet)]
    fn test_bond_from_str(
        #[case] input: &str,
        #[case] expected_order: u8,
        #[case] expected_charge: i8,
        #[case] expected_unpaired: u8,
        #[case] expected_multiplicity: SpinMultiplicity,
    ) {
        let bond: Bond = input.parse().expect("expected parse success");
        assert_eq!(bond.order(), expected_order);
        assert_eq!(bond.charge(), expected_charge);
        assert_eq!(bond.unpaired_electrons(), expected_unpaired);
        assert_eq!(bond.multiplicity(), expected_multiplicity);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(Bond::new(1), "1")]
    #[case::charged_doublet(Bond::from_parts(2, 1, SpinState::new(1, SpinMultiplicity::Doublet)), "2#c+#u")]
    #[case::triplet(Bond::from_parts(2, 0, SpinState::new(2, SpinMultiplicity::Triplet)), "2#u2")]
    fn test_bond_display(#[case] bond: Bond, #[case] expected: &str) {
        assert_eq!(bond.to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(Bond::new(1), r#""1""#)]
    #[case::charged_doublet(Bond::from_parts(2, 1, SpinState::new(1, SpinMultiplicity::Doublet)), r#""2#c+#u""#)]
    fn test_bond_serialize(#[case] bond: Bond, #[case] expected: &str) {
        let json = serde_json::to_string(&bond).unwrap();
        assert_eq!(json, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(r#""1""#, 1, 0, 0, SpinMultiplicity::Singlet)]
    #[case::charged_doublet(r#""2#c+#u""#, 2, 1, 1, SpinMultiplicity::Doublet)]
    fn test_bond_deserialize(
        #[case] input: &str,
        #[case] expected_order: u8,
        #[case] expected_charge: i8,
        #[case] expected_unpaired: u8,
        #[case] expected_multiplicity: SpinMultiplicity,
    ) {
        let bond: Bond = serde_json::from_str(input).unwrap();
        assert_eq!(bond.order(), expected_order);
        assert_eq!(bond.charge(), expected_charge);
        assert_eq!(bond.unpaired_electrons(), expected_unpaired);
        assert_eq!(bond.multiplicity(), expected_multiplicity);
    }
}
