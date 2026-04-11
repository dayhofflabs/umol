//! Bond types for GraphIR.

use std::fmt::{self, Display};
use std::str::FromStr;

use umol_data::{SpinMultiplicity, SpinState, SpinStateError};
use umol_edn::{DeError, Edn, FromEdn, ToEdn};

use super::ast_utils::raise_spin_ground;
use super::bond_pattern::BondPattern;
use super::error::ValidationError;
use crate::ast::bond::BondAst;
use crate::ast::config::{BondAstConfig, NumericMode};
use crate::ast::error::LoweringError;
use crate::ast::value::ValueAst;
use crate::ast::{FromAst, ToAst};
use crate::dsl::bond::parse_bond_dsl;

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
    fn from_ast(ast: &BondAst, cfg: &BondAstConfig) -> Result<Self, LoweringError> {
        let pattern = BondPattern::from_ast(ast, cfg)?;
        pattern.to_bond().map_err(|e| match e {
            ValidationError::NonGround { field } => LoweringError::NonGround { field },
            ValidationError::InvalidMultiplicity(n) => LoweringError::InvalidMultiplicity(n),
            ValidationError::SpinUnderdetermined => {
                LoweringError::SpinState(SpinStateError::Underdetermined)
            }
            ValidationError::SpinIncompatible {
                unpaired_electrons,
                multiplicity,
            } => LoweringError::SpinState(SpinStateError::Incompatible {
                unpaired_electrons,
                multiplicity,
            }),
            other => LoweringError::Atom(other.to_string()),
        })
    }
}

impl ToAst<BondAst> for Bond {
    fn to_ast(&self, cfg: &BondAstConfig) -> BondAst {
        let (spin_u, spin_m) = raise_spin_ground(
            self.unpaired_electrons(),
            self.multiplicity(),
            &cfg.unpaired_electrons_mode,
            &cfg.multiplicity_mode,
        );
        BondAst {
            order: ValueAst::Lit(self.order() as i32),
            charge: match (self.charge(), &cfg.charge_mode) {
                (0, NumericMode::Zero) => None,
                (n, _) => Some(ValueAst::Lit(n as i32)),
            },
            unpaired_electrons: spin_u,
            multiplicity: spin_m,
        }
    }
}

impl FromStr for Bond {
    type Err = LoweringError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let ast = parse_bond_dsl(s).map_err(|e| LoweringError::Atom(e.to_string()))?;
        Bond::from_ast(&ast, &BondAstConfig::zeroed())
    }
}

impl Display for Bond {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.to_ast(&BondAstConfig::zeroed()).fmt(f)
    }
}

impl<'de> FromEdn<'de> for Bond {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let ast = BondAst::from_edn(edn)?;
        Self::from_ast(&ast, &BondAstConfig::zeroed())
            .map_err(|e| DeError::subgrammar("bond", e))
    }
}

impl ToEdn for Bond {
    fn to_edn(&self) -> Edn<'static> {
        self.to_ast(&BondAstConfig::zeroed()).to_edn()
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
        assert_eq!(Bond::from_ast(&ast, &BondAstConfig::zeroed()).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(Bond::new(1), BondAst { order: ValueAst::Lit(1), charge: None, unpaired_electrons: None, multiplicity: None })]
    #[case::charged_doublet(Bond::from_parts(2, 1, SpinState::new(1, SpinMultiplicity::Doublet)), BondAst { order: ValueAst::Lit(2), charge: Some(ValueAst::Lit(1)),
        unpaired_electrons: Some(ValueAst::Lit(1)), multiplicity: None })]
    fn test_bond_to_ast(#[case] bond: Bond, #[case] expected: BondAst) {
        assert_eq!(bond.to_ast(&BondAstConfig::zeroed()), expected);
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
    #[case::single(Bond::new(1), ":single")]
    #[case::charged_doublet(Bond::from_parts(2, 1, SpinState::new(1, SpinMultiplicity::Doublet)), r#""2#c+#u""#)]
    fn test_bond_to_edn(#[case] bond: Bond, #[case] expected: &str) {
        assert_eq!(bond.to_edn().to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(":single", 1, 0, 0, SpinMultiplicity::Singlet)]
    #[case::charged_doublet(r#""2#c+#u""#, 2, 1, 1, SpinMultiplicity::Doublet)]
    fn test_bond_from_edn(
        #[case] input: &str,
        #[case] expected_order: u8,
        #[case] expected_charge: i8,
        #[case] expected_unpaired: u8,
        #[case] expected_multiplicity: SpinMultiplicity,
    ) {
        let bond = Bond::from_edn_str(input).unwrap();
        assert_eq!(bond.order(), expected_order);
        assert_eq!(bond.charge(), expected_charge);
        assert_eq!(bond.unpaired_electrons(), expected_unpaired);
        assert_eq!(bond.multiplicity(), expected_multiplicity);
    }
}
