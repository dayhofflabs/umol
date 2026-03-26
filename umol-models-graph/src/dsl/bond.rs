//! Bond-string DSL parser — `spec/umol-dsl-spec.md` §7.5.

use nom::character::complete::multispace0;
use nom::combinator::all_consuming;
use nom::multi::many0;
use nom::sequence::{delimited, pair, terminated};
use nom::{Err, IResult, Parser};

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
        unpaired: None,
        multiplicity: None,
    };
    update_bond_ast(&mut ast, preds).map_err(|e| Err::Error(e))?;
    Ok((remaining, ast))
}

/// Parse the bond order prefix (`value_dsl::<u8>`).
pub fn bond_order(i: &str) -> IResult<&str, ValueAst<u8>, ParseError> {
    value_dsl::<u8>(i).map_err(|_| Err::Failure(ParseError::InvalidBondOrder(i.to_string())))
}

/// Merge a list of bond predicates into a `BondAst`.
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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::single("1", BondAst { order: ValueAst::Lit(1), charge: None, unpaired: None, multiplicity: None })]
    #[case::double("2", BondAst { order: ValueAst::Lit(2), charge: None, unpaired: None, multiplicity: None })]
    #[case::triple("3", BondAst { order: ValueAst::Lit(3), charge: None, unpaired: None, multiplicity: None })]
    #[case::quadruple("4", BondAst { order: ValueAst::Lit(4), charge: None, unpaired: None, multiplicity: None })]
    #[case::single_whitespace("  1  ", BondAst { order: ValueAst::Lit(1), charge: None, unpaired: None, multiplicity: None })]
    #[case::single_pos_charge("1#c+2", BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(2)), unpaired: None, multiplicity: None })]
    #[case::single_neg_charge("1#c-2", BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(-2)), unpaired: None, multiplicity: None })]
    #[case::single_zero_charge("1#c0", BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(0)), unpaired: None, multiplicity: None })]
    #[case::single_plus_only("1#c+", BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(1)), unpaired: None, multiplicity: None })]
    #[case::single_minus_only("1#c-",  BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(-1)), unpaired: None, multiplicity: None })]
    #[case::single_plus_whitespace("1#c +", BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(1)), unpaired: None, multiplicity: None })]
    #[case::single_minus_whitespace("1#c -", BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(-1)), unpaired: None, multiplicity: None })]
    #[case::single_pos_charge_whitespace("1#c +2", BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(2)), unpaired: None, multiplicity: None })]
    #[case::double_unpaired("2#u3", BondAst { order: ValueAst::Lit(2), charge: None, unpaired: Some(ValueAst::Lit(3)), multiplicity: None })]
    #[case::single_u_only("1#u", BondAst { order: ValueAst::Lit(1), charge: None, unpaired: Some(ValueAst::Lit(1)), multiplicity: None })]
    #[case::single_mult("1#s2", BondAst { order: ValueAst::Lit(1), charge: None, unpaired: None, multiplicity: Some(ValueAst::Lit(2)) })]
    #[case::single_s_only("1#s", BondAst { order: ValueAst::Lit(1), charge: None, unpaired: None, multiplicity: Some(ValueAst::Lit(1)) })]
    #[case::double_charge_unpaired("2#c+#u2", BondAst { order: ValueAst::Lit(2), charge: Some(ValueAst::Lit(1)), unpaired: Some(ValueAst::Lit(2)), multiplicity: None })]
    #[case::double_charge_mult("2#c-1#s3", BondAst { order: ValueAst::Lit(2), charge: Some(ValueAst::Lit(-1)), unpaired: None, multiplicity: Some(ValueAst::Lit(3)) })]
    #[case::double_charge_unpaired_mult("1#c0#u1#s1", BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(0)), unpaired: Some(ValueAst::Lit(1)), multiplicity: Some(ValueAst::Lit(1)) })]
    #[case::double_plus_only_unpaired("1 #c+ #u2", BondAst { order: ValueAst::Lit(1), charge: Some(ValueAst::Lit(1)), unpaired: Some(ValueAst::Lit(2)), multiplicity: None })]
    fn test_bond_dsl(#[case] input: &str, #[case] expected: BondAst) {
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
    fn test_parse_invalid(#[case] input: &str, #[case] expected: ParseError) {
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
