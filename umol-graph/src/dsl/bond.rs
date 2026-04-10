//! Bond-string DSL: parser, AST, and display

use std::fmt::{self, Display};
use std::str::FromStr;

use nom::character::complete::multispace0;
use nom::combinator::all_consuming;
use nom::multi::many0;
use nom::sequence::{delimited, pair, terminated};
use nom::{Err, IResult, Parser};
use super::ast::DslAst;
use super::config::BondDslConfig;
use super::error::ParseError;
use super::predicates::{bond_order, bond_predicate, BondPredicate};
use super::value::ValueAst;

/// Parsed bond-string AST
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct BondAst {
    pub order: ValueAst,
    pub charge: Option<ValueAst>,
    pub unpaired_electrons: Option<ValueAst>,
    pub multiplicity: Option<ValueAst>,
}

impl BondAst {
    pub fn new(order: ValueAst) -> Self {
        Self {
            order,
            charge: None,
            unpaired_electrons: None,
            multiplicity: None,
        }
    }

    pub fn from_order(order: u8) -> Self {
        Self {
            order: ValueAst::Lit(order as i32),
            charge: None,
            unpaired_electrons: None,
            multiplicity: None,
        }
    }
}

impl DslAst for BondAst {
    type Config = BondDslConfig;
}

impl FromStr for BondAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_bond_dsl(s)
    }
}

impl Display for BondAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.order {
            ValueAst::Lit(n) => write!(f, "{}", n)?,
            ValueAst::Wildcard => write!(f, "*")?,
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
            None | Some(ValueAst::Lit(0)) => {}
            Some(ValueAst::Lit(1)) => write!(f, "#c+")?,
            Some(ValueAst::Lit(-1)) => write!(f, "#c-")?,
            Some(ValueAst::Lit(n)) if *n > 0 => write!(f, "#c+{}", n)?,
            Some(ValueAst::Lit(n)) => write!(f, "#c{}", n)?,
            Some(ValueAst::Wildcard) => write!(f, "#c*")?,
            Some(v) => {
                write!(f, "#c")?;
                fmt_bond_value(f, v)?;
            }
        }

        match &self.unpaired_electrons {
            None | Some(ValueAst::Lit(0)) => {}
            Some(ValueAst::Lit(1)) => write!(f, "#u")?,
            Some(ValueAst::Lit(n)) => write!(f, "#u{}", n)?,
            Some(ValueAst::Wildcard) => write!(f, "#u*")?,
            Some(v) => {
                write!(f, "#u")?;
                fmt_bond_value(f, v)?;
            }
        }

        let m = match &self.multiplicity {
            None => return Ok(()),
            Some(ValueAst::Lit(m)) => *m,
            Some(ValueAst::Wildcard) => return write!(f, "#s*"),
            Some(v) => {
                write!(f, "#s")?;
                return fmt_bond_value(f, v);
            }
        };
        let u: i32 = match &self.unpaired_electrons {
            Some(ValueAst::Lit(u)) => *u,
            None => 0,
            _ => -1,
        };
        if m as i32 != u + 1 {
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

impl<'de> umol_edn::FromEdn<'de> for BondAst {
    fn from_edn(edn: &umol_edn::Edn<'de>) -> Result<Self, umol_edn::DeError> {
        let s: &str = match edn {
            umol_edn::Edn::Str(s) => s,
            umol_edn::Edn::Keyword(k) => k.as_str(),
            other => {
                return Err(umol_edn::DeError::TypeMismatch {
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
        parse_bond_dsl(s).map_err(|e| umol_edn::DeError::Custom(e.to_string()))
    }
}

impl umol_edn::ToEdn for BondAst {
    fn to_edn(&self) -> umol_edn::Edn<'static> {
        let aliases = builtin_bond_aliases();
        if let Some(name) = aliases.get_by_right(self) {
            umol_edn::Edn::Keyword(umol_edn::EdnKeyword::owned(name.clone()))
        } else {
            umol_edn::Edn::Str(std::borrow::Cow::Owned(self.to_string()))
        }
    }
}

/// Parse a bond subgrammar string
pub fn parse_bond_dsl(input: &str) -> Result<BondAst, ParseError> {
    all_consuming(bond_dsl)
        .parse(input)
        .map(|(_, result)| result)
        .map_err(|e| match e {
            Err::Error(e) | Err::Failure(e) => e,
            Err::Incomplete(_) => ParseError::Incomplete,
        })
}

/// Bond subgrammar parser
pub fn bond_dsl(i: &str) -> IResult<&str, BondAst, ParseError> {
    let (remaining, (order, preds)) = pair(
        delimited(multispace0, bond_order, multispace0),
        many0(terminated(bond_predicate, multispace0)),
    )
    .parse(i)?;

    let mut ast = BondAst {
        order,
        charge: None,
        unpaired_electrons: None,
        multiplicity: None,
    };
    update_bond_ast(&mut ast, preds).map_err(Err::Error)?;
    Ok((remaining, ast))
}

/// Merge a list of bond predicates into a `BondAst`
fn update_bond_ast(ast: &mut BondAst, preds: Vec<BondPredicate>) -> Result<(), ParseError> {
    for pred in preds {
        match pred {
            BondPredicate::Charge(v) => {
                if ast.charge.is_some() {
                    return Err(ParseError::DuplicateBondPredicate("#c".to_string()));
                }
                ast.charge = Some(v);
            }
            BondPredicate::UnpairedElectrons(v) => {
                if ast.unpaired_electrons.is_some() {
                    return Err(ParseError::DuplicateBondPredicate("#u".to_string()));
                }
                ast.unpaired_electrons = Some(v);
            }
            BondPredicate::Multiplicity(v) => {
                if ast.multiplicity.is_some() {
                    return Err(ParseError::DuplicateBondPredicate("#s".to_string()));
                }
                ast.multiplicity = Some(v);
            }
        }
    }
    Ok(())
}

fn fmt_bond_value(f: &mut fmt::Formatter<'_>, v: &ValueAst) -> fmt::Result {
    match v {
        ValueAst::Wildcard => write!(f, "*"),
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
    #[case::single("1", BondAst { order: ValueAst::Lit(1), charge: None, unpaired_electrons:None, multiplicity: None })]
    #[case::double("2", BondAst { order: ValueAst::Lit(2), charge: None, unpaired_electrons:None, multiplicity: None })]
    #[case::triple("3", BondAst { order: ValueAst::Lit(3), charge: None, unpaired_electrons:None, multiplicity: None })]
    #[case::quadruple("4", BondAst { order: ValueAst::Lit(4), charge: None, unpaired_electrons:None, multiplicity: None })]
    #[case::single_whitespace("  1  ", BondAst { order: ValueAst::Lit(1), charge: None, unpaired_electrons:None, multiplicity: None })]
    #[case::single_pos_charge("1#c+2", BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(2)), unpaired_electrons:None, multiplicity: None })]
    #[case::single_neg_charge("1#c-2", BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(-2)), unpaired_electrons:None, multiplicity: None })]
    #[case::single_zero_charge("1#c0", BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(0)), unpaired_electrons:None, multiplicity: None })]
    #[case::single_plus_only("1#c+", BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(1)), unpaired_electrons:None, multiplicity: None })]
    #[case::single_minus_only("1#c-",  BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(-1)), unpaired_electrons:None, multiplicity: None })]
    #[case::single_plus_whitespace("1#c +", BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(1)), unpaired_electrons:None, multiplicity: None })]
    #[case::single_minus_whitespace("1#c -", BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(-1)), unpaired_electrons:None, multiplicity: None })]
    #[case::single_pos_charge_whitespace("1#c +2", BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(2)), unpaired_electrons:None, multiplicity: None })]
    #[case::double_unpaired("2#u3", BondAst { order: ValueAst::Lit(2), charge: None, unpaired_electrons:Some(ValueAst::Lit(3)), multiplicity: None })]
    #[case::single_u_only("1#u", BondAst { order: ValueAst::Lit(1), charge: None, unpaired_electrons:Some(ValueAst::Lit(1)), multiplicity: None })]
    #[case::single_mult("1#s2", BondAst { order: ValueAst::Lit(1), charge: None, unpaired_electrons:None, multiplicity: Some(ValueAst::Lit(2)) })]
    #[case::single_s_only("1#s", BondAst { order: ValueAst::Lit(1), charge: None, unpaired_electrons:None, multiplicity: Some(ValueAst::Lit(1)) })]
    #[case::double_charge_unpaired("2#c+#u2", BondAst { order: ValueAst::Lit(2), charge: Some(ValueAst::Lit(1)), unpaired_electrons:Some(ValueAst::Lit(2)), multiplicity: None })]
    #[case::double_charge_mult("2#c-1#s3", BondAst { order: ValueAst::Lit(2), charge: Some(ValueAst::Lit(-1)), unpaired_electrons:None, multiplicity: Some(ValueAst::Lit(3)) })]
    #[case::double_charge_unpaired_mult("1#c0#u1#s1", BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(0)), unpaired_electrons:Some(ValueAst::Lit(1)), multiplicity: Some(ValueAst::Lit(1)) })]
    #[case::double_plus_only_unpaired("1 #c+ #u2", BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(1)), unpaired_electrons:Some(ValueAst::Lit(2)), multiplicity: None })]
    fn test_parse_bond_dsl(#[case] input: &str, #[case] expected: BondAst) {
        let result = bond_dsl(input);
        assert!(
            result.is_ok(),
            "{input:?} should succeed, got {:?}",
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
    #[case::empty("", ParseError::InvalidBondOrder("".to_string()))]
    #[case::tag_whitespace("1# c", ParseError::UnknownBondPredicate("# ".to_string()))]
    #[case::invalid_tag("1#x", ParseError::UnknownBondPredicate("#x".to_string()))]
    #[case::trailing("1#c+ foo", ParseError::TrailingInput("foo".to_string()))]
    #[case::dup_charge("1#c+#c-", ParseError::DuplicateBondPredicate("#c".to_string()))]
    #[case::dup_unpaired("1#u2#u3", ParseError::DuplicateBondPredicate("#u".to_string()))]
    #[case::dup_multiplicity("1#s1#s2", ParseError::DuplicateBondPredicate("#s".to_string()))]
    fn test_parse_bond_dsl_invalid(#[case] input: &str, #[case] expected: ParseError) {
        let result = bond_dsl(input);
        assert!(
            result.is_err(),
            "{input:?} should fail, got {:?}",
            result.unwrap_err()
        );
        let err = match result.unwrap_err() {
            Err::Error(e) | Err::Failure(e) => e,
            Err::Incomplete(_) => ParseError::Incomplete,
        };
        assert_eq!(
            err, expected,
            "{:?} should fail with {:?}, got {:?}",
            input, expected, err
        );
    }
}
