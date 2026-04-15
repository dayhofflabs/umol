//! Bond-string DSL: parser, AST, and display

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use nom::branch::alt;
use nom::bytes::complete::{tag, take};
use nom::character::complete::multispace0;
use nom::combinator::{all_consuming, map, success, value};
use nom::multi::many0;
use nom::sequence::{delimited, pair, preceded, terminated};
use nom::{Err, IResult, Parser};
use umol_edn::{DeError, Edn, EdnKeyword, FromEdn, ToEdn};
use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;

use crate::ast::bond::BondAst;
use crate::dsl::error::BondDslError;
use crate::dsl::value::value_dsl;

impl FromStr for BondAst {
    type Err = BondDslError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_bond_dsl(s)
    }
}

impl Display for BondAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.order {
            ValueAst::Lit(n) => write!(f, "{}", n)?,
            ValueAst::Undetermined => write!(f, "*")?,
            v => {
                write!(f, "{{")?;
                if let ValueAst::LitSet(s) = v {
                    for (i, n) in s.iter().enumerate() {
                        if i > 0 {
                            write!(f, ",")?;
                        }
                        write!(f, "{}", n)?;
                    }
                }
                write!(f, "}}")?;
            }
        }

        match &self.charge {
            ValueAst::Undetermined | ValueAst::Lit(0) => {}
            ValueAst::Lit(1) => write!(f, "#c+")?,
            ValueAst::Lit(-1) => write!(f, "#c-")?,
            ValueAst::Lit(n) if *n > 0 => write!(f, "#c+{}", n)?,
            ValueAst::Lit(n) => write!(f, "#c{}", n)?,
            v => {
                write!(f, "#c")?;
                fmt_bond_value(f, v)?;
            }
        }

        let (u_field, m_field) = self.spin.to_pair();
        if matches!((&u_field, &m_field), (ValueAst::Undetermined, ValueAst::Undetermined)) {
            return Ok(());
        }

        match &u_field {
            ValueAst::Undetermined | ValueAst::Lit(0) => {}
            ValueAst::Lit(1) => write!(f, "#u")?,
            ValueAst::Lit(n) => write!(f, "#u{}", n)?,
            v => {
                write!(f, "#u")?;
                fmt_bond_value(f, v)?;
            }
        }

        let m = match &m_field {
            ValueAst::Undetermined => return Ok(()),
            ValueAst::Lit(m) => *m,
            v => {
                write!(f, "#s")?;
                return fmt_bond_value(f, v);
            }
        };
        let u: i64 = match &u_field {
            ValueAst::Lit(u) => *u,
            ValueAst::Undetermined => 0,
            _ => -1,
        };
        if m != u + 1 {
            if m == 1 {
                write!(f, "#s")?;
            } else {
                write!(f, "#s{}", m)?;
            }
        }

        Ok(())
    }
}

/// Built-in bond keyword aliases (:single, :double, :triple, :quadruple).
pub fn builtin_bond_aliases() -> bimap::BiMap<String, BondAst> {
    bimap::BiMap::from_iter([
        ("single".into(), BondAst::from_order(1)),
        ("double".into(), BondAst::from_order(2)),
        ("triple".into(), BondAst::from_order(3)),
        ("quadruple".into(), BondAst::from_order(4)),
    ])
}

impl<'de> FromEdn<'de> for BondAst {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let s: &str = match edn {
            Edn::Str(s) => s,
            Edn::Keyword(k) => k.as_str(),
            other => {
                return Err(DeError::TypeMismatch {
                    expected: "string or keyword",
                    got: other.kind(),
                    path: Vec::new(),
                });
            }
        };
        let aliases = builtin_bond_aliases();
        if let Some(ast) = aliases.get_by_left(s) {
            return Ok(ast.clone());
        }
        parse_bond_dsl(s).map_err(|e| umol_edn::DeError::subgrammar("bond", e))
    }
}

impl ToEdn for BondAst {
    fn to_edn(&self) -> Edn<'static> {
        let aliases = builtin_bond_aliases();
        if let Some(name) = aliases.get_by_right(self) {
            Edn::Keyword(EdnKeyword::owned(name.clone()))
        } else {
            Edn::Str(Cow::Owned(self.to_string()))
        }
    }
}

/// Parse a bond subgrammar string
pub fn parse_bond_dsl(input: &str) -> Result<BondAst, BondDslError> {
    all_consuming(bond_dsl)
        .parse(input)
        .map(|(_, result)| result)
        .map_err(|e| match e {
            Err::Error(e) | Err::Failure(e) => e,
            Err::Incomplete(_) => BondDslError::Incomplete,
        })
}

/// Bond subgrammar parser
pub fn bond_dsl(i: &str) -> IResult<&str, BondAst, BondDslError> {
    let (remaining, (order, preds)) = pair(
        delimited(multispace0, bond_order, multispace0),
        many0(terminated(bond_predicate, multispace0)),
    )
    .parse(i)?;

    let mut ast = BondAst::new(order);
    update_bond_ast(&mut ast, preds).map_err(Err::Error)?;
    Ok((remaining, ast))
}

/// Merge a list of bond predicates into a `BondAst`
fn update_bond_ast(ast: &mut BondAst, preds: Vec<BondPredicate>) -> Result<(), BondDslError> {
    for pred in preds {
        match pred {
            BondPredicate::Charge(v) => {
                if !matches!(ast.charge, ValueAst::Undetermined) {
                    return Err(BondDslError::DuplicateBondPredicate("#c".to_string()));
                }
                ast.charge = v;
            }
            BondPredicate::UnpairedElectrons(v) => {
                let SpinStateAst::Pair { unpaired, .. } = &mut ast.spin else {
                    unreachable!("default is Pair")
                };
                if !matches!(unpaired, ValueAst::Undetermined) {
                    return Err(BondDslError::DuplicateBondPredicate("#u".to_string()));
                }
                *unpaired = v;
            }
            BondPredicate::Multiplicity(v) => {
                let SpinStateAst::Pair { multiplicity, .. } = &mut ast.spin else {
                    unreachable!("default is Pair")
                };
                if !matches!(multiplicity, ValueAst::Undetermined) {
                    return Err(BondDslError::DuplicateBondPredicate("#s".to_string()));
                }
                *multiplicity = v;
            }
        }
    }
    Ok(())
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BondPredicate {
    Charge(ValueAst),
    UnpairedElectrons(ValueAst),
    Multiplicity(ValueAst),
}

fn bond_predicate(i: &str) -> IResult<&str, BondPredicate, BondDslError> {
    let (remaining, prefix) = take(2usize)(i)?;
    match prefix {
        "#c" => map(charge_value, BondPredicate::Charge).parse(remaining),
        "#u" => map(optional_value, BondPredicate::UnpairedElectrons).parse(remaining),
        "#s" => map(optional_value, BondPredicate::Multiplicity).parse(remaining),
        p if p.starts_with("#") => Err(Err::Failure(BondDslError::UnknownBondPredicate(
            p.to_string(),
        ))),
        _ => Err(Err::Failure(BondDslError::TrailingInput(i.to_string()))),
    }
}

fn bond_order(i: &str) -> IResult<&str, ValueAst, BondDslError> {
    value_dsl
        .parse(i)
        .map_err(|_| Err::Failure(BondDslError::InvalidBondOrder(i.to_string())))
}

fn charge_value(i: &str) -> IResult<&str, ValueAst, BondDslError> {
    preceded(
        multispace0,
        alt((
            value_dsl,
            value(ValueAst::Lit(1), tag("+")),
            value(ValueAst::Lit(-1), tag("-")),
        )),
    )
    .parse(i)
    .map_err(|_| Err::Failure(BondDslError::InvalidCharge(i.to_string())))
}

fn optional_value(i: &str) -> IResult<&str, ValueAst, BondDslError> {
    preceded(multispace0, alt((value_dsl, success(ValueAst::Lit(1)))))
        .parse(i)
        .map_err(|_| Err::Failure(BondDslError::InvalidValue(i.to_string())))
}

fn fmt_bond_value(f: &mut fmt::Formatter<'_>, v: &ValueAst) -> fmt::Result {
    match v {
        ValueAst::Undetermined => write!(f, "*"),
        ValueAst::Lit(n) => write!(f, "{}", n),
        ValueAst::LitSet(s) => {
            write!(f, "{{")?;
            for (i, n) in s.iter().enumerate() {
                if i > 0 {
                    write!(f, ",")?;
                }
                write!(f, "{}", n)?;
            }
            write!(f, "}}")
        }
        ValueAst::Expr(_) => write!(f, "<expr>"),
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::single("1", BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default() })]
    #[case::double("2", BondAst { order: ValueAst::Lit(2), charge: ValueAst::Undetermined, spin: SpinStateAst::default() })]
    #[case::triple("3", BondAst { order: ValueAst::Lit(3), charge: ValueAst::Undetermined, spin: SpinStateAst::default() })]
    #[case::quadruple("4", BondAst { order: ValueAst::Lit(4), charge: ValueAst::Undetermined, spin: SpinStateAst::default() })]
    #[case::single_whitespace("  1  ", BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default() })]
    #[case::single_pos_charge("1#c+2", BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(2), spin: SpinStateAst::default() })]
    #[case::single_neg_charge("1#c-2", BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(-2), spin: SpinStateAst::default() })]
    #[case::single_zero_charge("1#c0", BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(0), spin: SpinStateAst::default() })]
    #[case::single_plus_only("1#c+", BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(1), spin: SpinStateAst::default() })]
    #[case::single_minus_only("1#c-",  BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(-1), spin: SpinStateAst::default() })]
    #[case::single_plus_whitespace("1#c +", BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(1), spin: SpinStateAst::default() })]
    #[case::single_minus_whitespace("1#c -", BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(-1), spin: SpinStateAst::default() })]
    #[case::single_pos_charge_whitespace("1#c +2", BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(2), spin: SpinStateAst::default() })]
    #[case::double_unpaired("2#u3", BondAst { order: ValueAst::Lit(2), charge: ValueAst::Undetermined, spin: SpinStateAst::Pair { unpaired: ValueAst::Lit(3), multiplicity: ValueAst::Undetermined } })]
    #[case::single_u_only("1#u", BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::Pair { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Undetermined } })]
    #[case::single_mult("1#s2", BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::Pair { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(2) } })]
    #[case::single_s_only("1#s", BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::Pair { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(1) } })]
    #[case::double_charge_unpaired("2#c+#u2", BondAst { order: ValueAst::Lit(2), charge: ValueAst::Lit(1), spin: SpinStateAst::Pair { unpaired: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined } })]
    #[case::double_charge_mult("2#c-1#s3", BondAst { order: ValueAst::Lit(2), charge: ValueAst::Lit(-1), spin: SpinStateAst::Pair { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(3) } })]
    #[case::double_charge_unpaired_mult("1#c0#u1#s1", BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(0), spin: SpinStateAst::Pair { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Lit(1) } })]
    #[case::double_plus_only_unpaired("1 #c+ #u2", BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(1), spin: SpinStateAst::Pair { unpaired: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined } })]
    fn test_parse_bond_dsl(#[case] input: &str, #[case] expected: BondAst) {
        let result = bond_dsl(input);
        assert!(
            result.is_ok(),
            "{:?} should succeed, got {:?}",
            input,
            result.unwrap_err()
        );
        let (remaining, ast) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "{input:?} should have consumed all input, remaining: {remaining:?}"
        );
        assert_eq!(ast, expected);
    }

    #[rstest]
    #[case::empty("", BondDslError::InvalidBondOrder("".to_string()))]
    #[case::tag_whitespace("1# c", BondDslError::UnknownBondPredicate("# ".to_string()))]
    #[case::invalid_tag("1#x", BondDslError::UnknownBondPredicate("#x".to_string()))]
    #[case::trailing("1#c+ foo", BondDslError::TrailingInput("foo".to_string()))]
    #[case::dup_charge("1#c+#c-", BondDslError::DuplicateBondPredicate("#c".to_string()))]
    #[case::dup_unpaired("1#u2#u3", BondDslError::DuplicateBondPredicate("#u".to_string()))]
    #[case::dup_multiplicity("1#s1#s2", BondDslError::DuplicateBondPredicate("#s".to_string()))]
    fn test_parse_bond_dsl_invalid(#[case] input: &str, #[case] expected: BondDslError) {
        let result = bond_dsl(input);
        assert!(
            result.is_err(),
            "{:?} should fail, got {:?}",
            input,
            result.unwrap_err()
        );
        let err = match result.unwrap_err() {
            Err::Error(e) | Err::Failure(e) => e,
            Err::Incomplete(_) => BondDslError::Incomplete,
        };
        assert_eq!(
            err, expected,
            "{:?} should fail with {:?}, got {:?}",
            input, expected, err
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge_pos("#c+2", BondPredicate::Charge(ValueAst::Lit(2)))]
    #[case::charge_neg("#c-2", BondPredicate::Charge(ValueAst::Lit(-2)))]
    #[case::charge_plus("#c+", BondPredicate::Charge(ValueAst::Lit(1)))]
    #[case::charge_minus("#c-", BondPredicate::Charge(ValueAst::Lit(-1)))]
    #[case::charge_zero("#c0", BondPredicate::Charge(ValueAst::Lit(0)))]
    #[case::charge_undetermined("#c*", BondPredicate::Charge(ValueAst::Undetermined))]
    #[case::unpaired("#u2", BondPredicate::UnpairedElectrons(ValueAst::Lit(2)))]
    #[case::unpaired_omit("#u", BondPredicate::UnpairedElectrons(ValueAst::Lit(1)))]
    #[case::unpaired_undetermined("#u*", BondPredicate::UnpairedElectrons(ValueAst::Undetermined))]
    #[case::multiplicity("#s3", BondPredicate::Multiplicity(ValueAst::Lit(3)))]
    #[case::multiplicity_omit("#s", BondPredicate::Multiplicity(ValueAst::Lit(1)))]
    #[case::multiplicity_undetermined("#s*", BondPredicate::Multiplicity(ValueAst::Undetermined))]
    fn test_bond_predicate(#[case] input: &str, #[case] expected: BondPredicate) {
        let result = bond_predicate(input);
        assert!(result.is_ok(), "{input:?} should succeed, got {:?}", result.unwrap_err());
        let (_, pred) = result.unwrap();
        assert_eq!(pred, expected);
    }

    #[rstest]
    #[case::unknown("#x", BondDslError::UnknownBondPredicate("#x".to_string()))]
    #[case::unknown_tag("#z", BondDslError::UnknownBondPredicate("#z".to_string()))]
    #[case::trailing_no_hash("fo", BondDslError::TrailingInput("fo".to_string()))]
    fn test_bond_predicate_error(#[case] input: &str, #[case] expected: BondDslError) {
        let result = bond_predicate(input);
        assert!(result.is_err(), "{input:?} should fail, got {:?}", result.unwrap());
        let err = match result.unwrap_err() {
            Err::Error(e) | Err::Failure(e) => e,
            Err::Incomplete(_) => BondDslError::Incomplete,
        };
        assert_eq!(err, expected);
    }
}
