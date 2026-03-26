//! Bond-string DSL parser — `spec/umol-dsl-spec.md` §7.5.

use nom::character::complete::multispace0;
use nom::multi::many0;
use nom::sequence::{delimited, pair, terminated};
use nom::{Err as NomErr, IResult, Parser};

use super::error::ParseError;
use super::predicates::{bond_predicate, BondPredicate};
use super::value::{value_dsl, ValueAst};

/// Parsed bond-string AST.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BondAst {
    pub order: ValueAst<u8>,
    pub charge: Option<ValueAst<i8>>,
    pub unpaired: Option<ValueAst<u8>>,
    pub multiplicity: Option<ValueAst<u8>>,
}

/// Parse the bond order prefix (`value_dsl::<u8>`).
pub fn bond_order(i: &str) -> IResult<&str, ValueAst<u8>, ParseError> {
    value_dsl::<u8>(i)
}

/// Merge a list of bond predicates into a `BondAst`.
///
/// Returns `Err(ParseError::DuplicateBondPredicate(tag))` on duplicate.
fn update_bond_ast(ast: &mut BondAst, preds: Vec<BondPredicate>) -> Result<(), ParseError> {
    for pred in preds {
        match pred {
            BondPredicate::Charge(v) => {
                if ast.charge.is_some() {
                    return Err(ParseError::DuplicateBondPredicate("#c".to_string()));
                }
                ast.charge = Some(v);
            }
            BondPredicate::Unpaired(v) => {
                if ast.unpaired.is_some() {
                    return Err(ParseError::DuplicateBondPredicate("#u".to_string()));
                }
                ast.unpaired = Some(v);
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

/// Combinator parser for a full bond-string (without `all_consuming`).
pub fn bond_dsl(i: &str) -> IResult<&str, BondAst, ParseError> {
    let (remaining, (order, preds)) = pair(
        delimited(multispace0, bond_order, multispace0),
        many0(terminated(bond_predicate, multispace0)),
    )
    .parse(i)?;

    let mut ast = BondAst {
        order,
        charge: None,
        unpaired: None,
        multiplicity: None,
    };
    update_bond_ast(&mut ast, preds).map_err(|e| NomErr::Error(e))?;
    Ok((remaining, ast))
}

/// Top-level entry point: parse a complete bond-string.
///
/// Errors are domain-meaningful; `ParseError::NomError` never escapes this function.
pub fn parse_bond_dsl(input: &str) -> Result<BondAst, ParseError> {
    match bond_dsl(input) {
        Ok((remaining, ast)) => {
            let rest = remaining.trim_start_matches(|c: char| c.is_ascii_whitespace());
            if rest.is_empty() {
                Ok(ast)
            } else if rest.starts_with('#') {
                Err(ParseError::UnknownBondPredicate)
            } else {
                Err(ParseError::TrailingContent)
            }
        }
        Err(NomErr::Incomplete(_)) => Err(ParseError::Incomplete),
        Err(NomErr::Error(e) | NomErr::Failure(e)) => match e {
            dup @ ParseError::DuplicateBondPredicate(_) => Err(dup),
            _ => Err(ParseError::InvalidBondOrder),
        },
    }
}

#[cfg(test)]
mod tests {
    use nom::Parser;
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    use super::*;

    fn lit_u8(n: u8) -> ValueAst<u8> {
        ValueAst::Lit(n)
    }

    fn lit_i8(n: i8) -> ValueAst<i8> {
        ValueAst::Lit(n)
    }

    fn ast(
        order: u8,
        charge: Option<i8>,
        unpaired: Option<u8>,
        multiplicity: Option<u8>,
    ) -> BondAst {
        BondAst {
            order: lit_u8(order),
            charge: charge.map(lit_i8),
            unpaired: unpaired.map(lit_u8),
            multiplicity: multiplicity.map(lit_u8),
        }
    }

    #[rstest]
    #[case("1", ast(1, None, None, None))]
    #[case("2", ast(2, None, None, None))]
    #[case("3", ast(3, None, None, None))]
    #[case("  1  ", ast(1, None, None, None))]
    #[case("1#c+2", ast(1, Some(2), None, None))]
    #[case("1#c-2", ast(1, Some(-2), None, None))]
    #[case("1#c0", ast(1, Some(0), None, None))]
    #[case("1#c+", ast(1, Some(1), None, None))]
    #[case("1#c-",  ast(1, Some(-1), None, None))]
    #[case("1#c +", ast(1, Some(1), None, None))]
    #[case("1#c -", ast(1, Some(-1), None, None))]
    #[case("1#c +2", ast(1, Some(2), None, None))]
    #[case("1#u3", ast(1, None, Some(3), None))]
    #[case("1#u", ast(1, None, Some(1), None))]
    #[case("1#s2", ast(1, None, None, Some(2)))]
    #[case("1#s", ast(1, None, None, Some(1)))]
    #[case("2#c+#u2", ast(2, Some(1), Some(2), None))]
    #[case("2#c-1#s3", ast(2, Some(-1), None, Some(3)))]
    #[case("1#c0#u1#s1", ast(1, Some(0), Some(1), Some(1)))]
    #[case("1 #c+ #u2", ast(1, Some(1), Some(2), None))]
    fn test_parse_ok(#[case] input: &str, #[case] expected: BondAst) {
        let result = parse_bond_dsl(input);
        assert!(
            result.is_ok(),
            "{input:?} should succeed, got {:?}",
            result.unwrap_err()
        );
        assert_eq!(result.unwrap(), expected);
    }

    #[rstest]
    #[case("1# c", ParseError::UnknownBondPredicate)]
    #[case("1#x", ParseError::UnknownBondPredicate)]
    #[case("1#c+ foo", ParseError::TrailingContent)]
    #[case("", ParseError::InvalidBondOrder)]
    fn test_parse_err(#[case] input: &str, #[case] expected: ParseError) {
        let result = parse_bond_dsl(input);
        assert_eq!(result, Err(expected), "{input:?}");
    }

    #[rstest]
    #[case("1#c+#c-", "#c")]
    #[case("1#u2#u3", "#u")]
    #[case("1#s1#s2", "#s")]
    fn test_parse_duplicate(#[case] input: &str, #[case] tag: &str) {
        let result = parse_bond_dsl(input);
        assert_eq!(
            result,
            Err(ParseError::DuplicateBondPredicate(tag.to_string())),
            "{input:?}"
        );
    }

    #[test]
    fn test_bond_dsl_partial_remainder() {
        let result = bond_dsl.parse("1#c+ rest");
        assert!(result.is_ok(), "expected Ok, got {:?}", result.unwrap_err());
        let (remaining, ast) = result.unwrap();
        assert_eq!(remaining, "rest");
        assert_eq!(
            ast,
            BondAst {
                order: lit_u8(1),
                charge: Some(lit_i8(1)),
                unpaired: None,
                multiplicity: None
            }
        );
    }

    #[test]
    fn test_bond_dsl_stops_before_hash_unknown() {
        // bond_dsl stops when it encounters an unknown predicate tag, leaving it unconsumed
        let result = bond_dsl.parse("2#u #x remaining");
        assert!(result.is_ok(), "expected Ok, got {:?}", result.unwrap_err());
        let (remaining, ast) = result.unwrap();
        assert_eq!(remaining, "#x remaining");
        assert_eq!(
            ast,
            BondAst {
                order: lit_u8(2),
                charge: None,
                unpaired: Some(lit_u8(1)),
                multiplicity: None
            }
        );
    }
}
