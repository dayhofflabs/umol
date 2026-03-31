//! Bond pattern types for GraphIR.

use std::fmt::{self, Display};
use std::str::FromStr;

use serde::de::{Deserializer, Error as SerdeError};
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};
use umol_data::{SpinMultiplicity, SpinState};

use crate::dsl::ast::{FromAst, ToAst};
use crate::dsl::bond::{parse_bond_dsl, BondAst, BondLowerConfig, ChargeMode};
use crate::dsl::error::LoweringError;
use crate::dsl::value::ValueAst;
use crate::graph_ir::atom_pattern::Pattern;
use crate::graph_ir::bond::{Bond, BondError};
use crate::graph_ir::error::ValidationError;
use crate::table_ir::bond::{Bond as TableBond, BondOrder};

/// Bond pattern: carries order/charge/spin fields as `Pattern<T>` through the
/// resolution pipeline, grounding to `Bond` via `to_bond()`.
///
/// Fields are public to allow direct mutation in transform/kekulize code.
/// Aromaticity hints are stored separately on `MoleculeBuilder`, not here.
#[derive(Debug, Clone, PartialEq)]
pub struct BondPattern {
    pub order: Pattern<u8>,
    pub charge: Pattern<i8>,
    pub unpaired_electrons: Pattern<u8>,
    pub multiplicity: Pattern<SpinMultiplicity>,
}

impl BondPattern {
    /// Minimal pattern: concrete order, all other fields unconstrained.
    pub fn new(order: u8) -> Self {
        Self {
            order: Pattern::Is(order),
            charge: Pattern::Any,
            unpaired_electrons: Pattern::Any,
            multiplicity: Pattern::Any,
        }
    }

    /// Create a pattern from a concrete ground bond.
    pub fn from_bond(bond: &Bond) -> Self {
        Self {
            order: Pattern::Is(bond.order()),
            charge: Pattern::Is(bond.charge()),
            unpaired_electrons: Pattern::Is(bond.unpaired_electrons()),
            multiplicity: Pattern::Is(bond.multiplicity()),
        }
    }

    /// Create a pattern from a table IR bond.
    ///
    /// Aromatic order is normalized to 1; the caller is responsible for
    /// recording the aromatic hint on `MoleculeBuilder` via
    /// `set_bond_aromatic_hint`.
    pub fn from_table_bond(bond: &TableBond) -> Self {
        debug_assert!(
            !bond.order.is_query(),
            "query bond orders must be resolved before conversion to BondPattern"
        );
        let order = match bond.order {
            BondOrder::Aromatic => Pattern::Is(1),
            o => Pattern::Is(
                o.value()
                    .expect("non-query, non-aromatic bond order must have a value"),
            ),
        };
        Self {
            order,
            charge: bond.charge.map_or(Pattern::Any, Pattern::Is),
            unpaired_electrons: bond.unpaired_electrons.map_or(Pattern::Any, Pattern::Is),
            multiplicity: bond.multiplicity.map_or(Pattern::Any, Pattern::Is),
        }
    }

    /// Return the concrete order. Panics if `order` is `Any`.
    pub fn order(&self) -> u8 {
        match self.order {
            Pattern::Is(o) => o,
            Pattern::Any => panic!("bond order is not ground"),
        }
    }

    /// Ground this pattern into a concrete bond.
    ///
    /// `order` must be `Is`; `Any` on other fields defaults to 0 / closed-shell.
    pub fn to_bond(&self) -> Result<Bond, ValidationError> {
        let order = match self.order {
            Pattern::Is(o) => o,
            Pattern::Any => return Err(ValidationError::NonGround { field: "order" }),
        };

        let charge = match self.charge {
            Pattern::Any => 0,
            Pattern::Is(c) => c,
        };

        let electrons: i16 = 2 * order as i16 - charge as i16;
        if electrons < 0 {
            return Err(ValidationError::Bond(BondError::InvalidState(format!(
                "bond electron count is negative: order={order}, charge={charge}"
            ))));
        }

        let unpaired = self.unpaired_electrons.into_option();
        let mult = self.multiplicity.into_option();
        let spin = match (unpaired, mult) {
            (Some(u), Some(m)) => SpinState::try_new(u, m).map_err(BondError::SpinState)?,
            (Some(u), None) => SpinState::max_multiplicity(u).ok_or_else(|| {
                BondError::InvalidState(format!("unpaired electrons {u} out of range"))
            })?,
            (None, Some(m)) => {
                SpinState::try_new(m.multiplicity() - 1, m).map_err(BondError::SpinState)?
            }
            (None, None) => SpinState::closed_shell(),
        };

        if !spin.is_compatible_with(electrons as u8) {
            return Err(ValidationError::Bond(BondError::InvalidState(format!(
                "bond spin is not compatible with electron count: order={order}, charge={charge}, \
                 unpaired_electrons={unpaired:?}, multiplicity={mult:?}"
            ))));
        }

        Ok(Bond::from_parts(order, charge, spin))
    }

    pub fn matches_bond(&self, bond: &Bond) -> bool {
        self.order.matches(bond.order())
            && self.charge.matches(bond.charge())
            && self.unpaired_electrons.matches(bond.unpaired_electrons())
            && self.multiplicity.matches(bond.multiplicity())
    }
}

impl FromAst<BondAst> for BondPattern {
    fn from_ast(ast: BondAst, cfg: &BondLowerConfig) -> Result<Self, LoweringError> {
        let order = match ast.order {
            ValueAst::Lit(n) => {
                Pattern::Is(u8::try_from(n).map_err(|_| LoweringError::OutOfRange {
                    field: "order",
                    value: n as i64,
                })?)
            }
            ValueAst::Wildcard => Pattern::Any,
            _ => return Err(LoweringError::NonGround { field: "order" }),
        };

        let charge = match ast.charge.or_else(|| match cfg.charge_mode {
            ChargeMode::Zero => Some(ValueAst::Lit(0)),
            ChargeMode::Provided => None,
        }) {
            None => Pattern::Any,
            Some(ValueAst::Wildcard) => Pattern::Any,
            Some(ValueAst::Lit(n)) => Pattern::Is(
                i8::try_from(n).map_err(|_| LoweringError::NonGround { field: "charge" })?,
            ),
            Some(_) => return Err(LoweringError::NonGround { field: "charge" }),
        };

        let lower_u8_opt =
            |v: Option<ValueAst>, field: &'static str| -> Result<Pattern<u8>, LoweringError> {
                match v {
                    None => Ok(Pattern::Any),
                    Some(ValueAst::Wildcard) => Ok(Pattern::Any),
                    Some(ValueAst::Lit(n)) => u8::try_from(n)
                        .map(Pattern::Is)
                        .map_err(|_| LoweringError::NonGround { field }),
                    Some(_) => Err(LoweringError::NonGround { field }),
                }
            };

        let multiplicity = match ast.multiplicity {
            None => Pattern::Any,
            Some(ValueAst::Wildcard) => Pattern::Any,
            Some(ValueAst::Lit(n)) => {
                let m = u8::try_from(n).map_err(|_| LoweringError::NonGround {
                    field: "multiplicity",
                })?;
                Pattern::Is(
                    SpinMultiplicity::from_multiplicity(m)
                        .ok_or(LoweringError::InvalidMultiplicity(m))?,
                )
            }
            Some(_) => {
                return Err(LoweringError::NonGround {
                    field: "multiplicity",
                })
            }
        };

        Ok(BondPattern {
            order,
            charge,
            unpaired_electrons: lower_u8_opt(ast.unpaired_electrons, "unpaired_electrons")?,
            multiplicity,
        })
    }
}

impl ToAst<BondAst> for BondPattern {
    fn to_ast(&self) -> BondAst {
        BondAst {
            order: match self.order {
                Pattern::Any => ValueAst::Wildcard,
                Pattern::Is(n) => ValueAst::Lit(n as i32),
            },
            charge: match self.charge {
                Pattern::Any => None,
                Pattern::Is(n) => Some(ValueAst::Lit(n as i32)),
            },
            unpaired_electrons: match self.unpaired_electrons {
                Pattern::Any => None,
                Pattern::Is(n) => Some(ValueAst::Lit(n as i32)),
            },
            multiplicity: match self.multiplicity {
                Pattern::Any => None,
                Pattern::Is(m) => Some(ValueAst::Lit(m.multiplicity() as i32)),
            },
        }
    }
}

impl FromStr for BondPattern {
    type Err = LoweringError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let ast = parse_bond_dsl(s).map_err(|e| LoweringError::Atom(e.to_string()))?;
        Self::from_ast(ast, &BondLowerConfig::default())
    }
}

impl Display for BondPattern {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.to_ast().fmt(f)
    }
}

impl Serialize for BondPattern {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for BondPattern {
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
    use crate::table_ir::bond::BondOrder;

    #[rustfmt::skip]
    #[rstest]
    #[case::closed_shell(BondPattern::new(1), 0, 0, SpinMultiplicity::Singlet)]
    #[case::high_spin(BondPattern { charge: Pattern::Is(1), unpaired_electrons: Pattern::Is(1), ..BondPattern::new(1) }, 1, 1, SpinMultiplicity::Doublet)]
    #[case::from_multiplicity(BondPattern { charge: Pattern::Is(0), multiplicity: Pattern::Is(SpinMultiplicity::Triplet), ..BondPattern::new(2) }, 0, 2, SpinMultiplicity::Triplet)]
    #[case::complete(BondPattern { charge: Pattern::Is(0), unpaired_electrons: Pattern::Is(2), multiplicity: Pattern::Is(SpinMultiplicity::Singlet), ..BondPattern::new(1) }, 0, 2, SpinMultiplicity::Singlet)]
    fn test_bond_pattern_to_bond(
        #[case] pattern: BondPattern,
        #[case] expected_charge: i8,
        #[case] expected_unpaired: u8,
        #[case] expected_multiplicity: SpinMultiplicity,
    ) {
        let bond = pattern.to_bond().expect("expected ground success");
        assert_eq!(bond.order(), pattern.order());
        assert_eq!(bond.charge(), expected_charge);
        assert_eq!(bond.unpaired_electrons(), expected_unpaired);
        assert_eq!(bond.multiplicity(), expected_multiplicity);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::negative_electrons(BondPattern { charge: Pattern::Is(1), ..BondPattern::new(0) })]
    #[case::incompatible_spin_pair(BondPattern { unpaired_electrons: Pattern::Is(0), multiplicity: Pattern::Is(SpinMultiplicity::Triplet), ..BondPattern::new(1) })]
    #[case::electron_parity_mismatch(BondPattern { charge: Pattern::Is(0), unpaired_electrons: Pattern::Is(1), ..BondPattern::new(1) })]
    #[case::max_unpaired_exceeded(BondPattern { charge: Pattern::Is(0), unpaired_electrons: Pattern::Is(10), ..BondPattern::new(1) })]
    fn test_bond_pattern_to_bond_error(#[case] pattern: BondPattern) {
        assert!(pattern.to_bond().is_err(), "{:?} should have failed", pattern);
    }

    #[rstest]
    #[case::aromatic(1, BondOrder::Aromatic, None, None, None)]
    #[case::single(1, BondOrder::Single, None, None, None)]
    #[case::double(2, BondOrder::Double, Some(-1), Some(1), Some(SpinMultiplicity::Doublet))]
    fn test_bond_pattern_from_table_bond(
        #[case] expected_order: u8,
        #[case] order: BondOrder,
        #[case] charge: Option<i8>,
        #[case] unpaired_electrons: Option<u8>,
        #[case] multiplicity: Option<SpinMultiplicity>,
    ) {
        let mut bond = TableBond::new(0, 1, order);
        bond.charge = charge;
        bond.unpaired_electrons = unpaired_electrons;
        bond.multiplicity = multiplicity;

        let pattern = BondPattern::from_table_bond(&bond);
        assert_eq!(pattern.order, Pattern::Is(expected_order));
        assert_eq!(pattern.charge, charge.map_or(Pattern::Any, Pattern::Is));
        assert_eq!(
            pattern.unpaired_electrons,
            unpaired_electrons.map_or(Pattern::Any, Pattern::Is)
        );
        assert_eq!(
            pattern.multiplicity,
            multiplicity.map_or(Pattern::Any, Pattern::Is)
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::defaults(BondAst { order: ValueAst::Lit(1), charge: None, unpaired_electrons: None, multiplicity: None },
        BondPattern::new(1))]
    #[case::charged(BondAst { order: ValueAst::Lit(2), charge: Some(ValueAst::Lit(1)), unpaired_electrons: Some(ValueAst::Lit(1)), multiplicity: None },
        BondPattern { charge: Pattern::Is(1), unpaired_electrons: Pattern::Is(1), ..BondPattern::new(2) })]
    #[case::wildcard_charge(BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Wildcard), unpaired_electrons: None, multiplicity: None },
        BondPattern::new(1))]
    #[case::wildcard_order(BondAst { order: ValueAst::Wildcard, charge: None, unpaired_electrons: None, multiplicity: None },
        BondPattern { order: Pattern::Any, ..BondPattern::new(1) })]
    #[case::full(BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(0)), unpaired_electrons: Some(ValueAst::Lit(2)), multiplicity: Some(ValueAst::Lit(1)) },
        BondPattern { charge: Pattern::Is(0), unpaired_electrons: Pattern::Is(2), multiplicity: Pattern::Is(SpinMultiplicity::Singlet), ..BondPattern::new(1) })]
    fn test_bond_pattern_from_ast(#[case] ast: BondAst, #[case] expected: BondPattern) {
        assert_eq!(BondPattern::from_ast(ast, &BondLowerConfig::default()).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::minimal(BondPattern::new(1),
        BondAst { order: ValueAst::Lit(1), charge: None, unpaired_electrons: None, multiplicity: None })]
    #[case::charged(BondPattern { charge: Pattern::Is(1), unpaired_electrons: Pattern::Is(1), ..BondPattern::new(2) },
        BondAst { order: ValueAst::Lit(2), charge: Some(ValueAst::Lit(1)), unpaired_electrons: Some(ValueAst::Lit(1)), multiplicity: None })]
    #[case::wildcard_order(BondPattern { order: Pattern::Any, ..BondPattern::new(1) },
        BondAst { order: ValueAst::Wildcard, charge: None, unpaired_electrons: None, multiplicity: None })]
    #[case::full(BondPattern { charge: Pattern::Is(0), unpaired_electrons: Pattern::Is(2), multiplicity: Pattern::Is(SpinMultiplicity::Triplet), ..BondPattern::new(2) },
        BondAst { order: ValueAst::Lit(2), charge: Some(ValueAst::Lit(0)), unpaired_electrons: Some(ValueAst::Lit(2)), multiplicity: Some(ValueAst::Lit(3)) })]
    fn test_bond_pattern_to_ast(#[case] pattern: BondPattern, #[case] expected: BondAst) {
        assert_eq!(pattern.to_ast(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single("1", Pattern::Is(1), Pattern::Any, Pattern::Any, Pattern::Any)]
    #[case::double_charged("2#c+", Pattern::Is(2), Pattern::Is(1), Pattern::Any, Pattern::Any)]
    #[case::full("1#c0#u1#s2", Pattern::Is(1), Pattern::Is(0), Pattern::Is(1), Pattern::Is(SpinMultiplicity::Doublet))]
    fn test_bond_pattern_from_str(
        #[case] input: &str,
        #[case] expected_order: Pattern<u8>,
        #[case] expected_charge: Pattern<i8>,
        #[case] expected_unpaired: Pattern<u8>,
        #[case] expected_multiplicity: Pattern<SpinMultiplicity>,
    ) {
        let pattern: BondPattern = input.parse().expect("expected parse success");
        assert_eq!(pattern.order, expected_order);
        assert_eq!(pattern.charge, expected_charge);
        assert_eq!(pattern.unpaired_electrons, expected_unpaired);
        assert_eq!(pattern.multiplicity, expected_multiplicity);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(BondPattern::new(1), "1")]
    #[case::charged(BondPattern { charge: Pattern::Is(1), unpaired_electrons: Pattern::Is(1), ..BondPattern::new(2) }, "2#c+#u")]
    #[case::wildcard_order(BondPattern { order: Pattern::Any, ..BondPattern::new(1) }, "*")]
    #[case::triplet(BondPattern { unpaired_electrons: Pattern::Is(2), ..BondPattern::new(2) }, "2#u2")]
    fn test_bond_pattern_display(#[case] pattern: BondPattern, #[case] expected: &str) {
        assert_eq!(pattern.to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(BondPattern::new(1), r#""1""#)]
    #[case::charged(BondPattern { charge: Pattern::Is(1), unpaired_electrons: Pattern::Is(1), ..BondPattern::new(2) }, r#""2#c+#u""#)]
    #[case::wildcard_order(BondPattern { order: Pattern::Any, ..BondPattern::new(1) }, r#""*""#)]
    fn test_bond_pattern_serialize(#[case] pattern: BondPattern, #[case] expected: &str) {
        let json = serde_json::to_string(&pattern).unwrap();
        assert_eq!(json, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(r#""1""#, Pattern::Is(1), Pattern::Any, Pattern::Any, Pattern::Any)]
    #[case::charged(r#""2#c+#u""#, Pattern::Is(2), Pattern::Is(1), Pattern::Is(1), Pattern::Any)]
    #[case::full(r#""1#c0#u1#s2""#, Pattern::Is(1), Pattern::Is(0), Pattern::Is(1), Pattern::Is(SpinMultiplicity::Doublet))]
    fn test_bond_pattern_deserialize(
        #[case] input: &str,
        #[case] expected_order: Pattern<u8>,
        #[case] expected_charge: Pattern<i8>,
        #[case] expected_unpaired: Pattern<u8>,
        #[case] expected_multiplicity: Pattern<SpinMultiplicity>,
    ) {
        let pattern: BondPattern = serde_json::from_str(input).unwrap();
        assert_eq!(pattern.order, expected_order);
        assert_eq!(pattern.charge, expected_charge);
        assert_eq!(pattern.unpaired_electrons, expected_unpaired);
        assert_eq!(pattern.multiplicity, expected_multiplicity);
    }
}
