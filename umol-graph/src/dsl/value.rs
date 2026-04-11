//! `value-dsl` — `spec/umol-dsl-spec.md` §5

use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::{char, i64 as nom_i64, multispace0, satisfy};
use nom::combinator::{all_consuming, map, opt, recognize, value};
use nom::error::{Error as NomError, ErrorKind};
use nom::multi::{many0, separated_list1};
use nom::sequence::{delimited, pair, preceded, terminated};
use nom::{Err, IResult, Parser};

use umol_shared::value_ast::{ArithOp, Expr, RelOp, ValueAst};

pub fn parse_value_dsl(input: &str) -> Result<ValueAst, Err<NomError<&str>>> {
    all_consuming(value_dsl).parse(input).map(|(_, v)| v)
}

pub fn value_dsl(i: &str) -> IResult<&str, ValueAst, NomError<&str>> {
    alt((
        map(
            terminated(nom_i64, (multispace0, terminator)),
            ValueAst::Lit,
        ),
        value(ValueAst::Wildcard, tag("*")),
        map(lit_set, ValueAst::LitSet),
        map(bool_expr, ValueAst::Expr),
    ))
    .parse(i)
}

fn terminator(i: &str) -> IResult<&str, (), NomError<&str>> {
    if i.is_empty() || i.starts_with('#') {
        Ok((i, ()))
    } else {
        Err(Err::Error(NomError::new(i, ErrorKind::Tag)))
    }
}

fn bool_expr(i: &str) -> IResult<&str, Expr, NomError<&str>> {
    map(
        pair(and_expr, many0(preceded(op_char('|'), and_expr))),
        |(first, rest)| {
            if rest.is_empty() {
                first
            } else {
                let mut disjuncts = vec![first];
                disjuncts.extend(rest);
                Expr::Or(disjuncts)
            }
        },
    )
    .parse(i)
}

fn and_expr(i: &str) -> IResult<&str, Expr, NomError<&str>> {
    map(
        pair(not_expr, many0(preceded(op_char('&'), not_expr))),
        |(first, rest)| {
            if rest.is_empty() {
                first
            } else {
                let mut conjuncts = vec![first];
                conjuncts.extend(rest);
                Expr::And(conjuncts)
            }
        },
    )
    .parse(i)
}

fn not_expr(i: &str) -> IResult<&str, Expr, NomError<&str>> {
    alt((
        map(preceded((char('!'), multispace0), not_expr), |n| {
            Expr::Not(Box::new(n))
        }),
        rel_expr,
        map(
            delimited(
                char('('),
                delimited(multispace0, bool_expr, multispace0),
                char(')'),
            ),
            |b| b,
        ),
    ))
    .parse(i)
}

fn rel_expr(i: &str) -> IResult<&str, Expr, NomError<&str>> {
    map(
        pair(
            mem_expr,
            opt(preceded(
                multispace0,
                pair(rel_op, preceded(multispace0, mem_expr)),
            )),
        ),
        |(left, right)| match right {
            None => left,
            Some((op, r)) => Expr::Rel(Box::new(left), op, Box::new(r)),
        },
    )
    .parse(i)
}

fn mem_expr(i: &str) -> IResult<&str, Expr, NomError<&str>> {
    map(
        pair(
            add_expr,
            opt(preceded(
                multispace0,
                preceded(map(tag("::"), |_| ()), preceded(multispace0, lit_set)),
            )),
        ),
        |(expr, set)| match set {
            None => expr,
            Some(s) => Expr::Mem(Box::new(expr), s),
        },
    )
    .parse(i)
}

pub(crate) fn lit_set(i: &str) -> IResult<&str, Vec<i64>, NomError<&str>> {
    delimited(
        char('{'),
        delimited(
            multispace0,
            separated_list1(op_char(','), nom_i64),
            multispace0,
        ),
        char('}'),
    )
    .parse(i)
}

fn rel_op(i: &str) -> IResult<&str, RelOp, NomError<&str>> {
    alt((
        value(RelOp::Le, tag("<=")),
        value(RelOp::Ge, tag(">=")),
        value(RelOp::Eq, tag("==")),
        value(RelOp::Lt, char('<')),
        value(RelOp::Gt, char('>')),
    ))
    .parse(i)
}

fn add_expr(i: &str) -> IResult<&str, Expr, NomError<&str>> {
    map(
        pair(
            mult_expr,
            many0(pair(delimited(multispace0, add_op, multispace0), mult_expr)),
        ),
        |(head, tail)| {
            tail.into_iter().fold(head, |acc, (op, rhs)| {
                Expr::BinOp(Box::new(acc), op, Box::new(rhs))
            })
        },
    )
    .parse(i)
}

fn add_op(i: &str) -> IResult<&str, ArithOp, NomError<&str>> {
    alt((
        value(ArithOp::Add, char('+')),
        value(ArithOp::Sub, char('-')),
    ))
    .parse(i)
}

fn mult_expr(i: &str) -> IResult<&str, Expr, NomError<&str>> {
    map(
        pair(
            unary_expr,
            many0(pair(
                delimited(multispace0, mult_op, multispace0),
                unary_expr,
            )),
        ),
        |(head, tail)| {
            tail.into_iter().fold(head, |acc, (op, rhs)| {
                Expr::BinOp(Box::new(acc), op, Box::new(rhs))
            })
        },
    )
    .parse(i)
}

fn mult_op(i: &str) -> IResult<&str, ArithOp, NomError<&str>> {
    alt((
        value(ArithOp::Mul, char('*')),
        value(ArithOp::Div, char('/')),
        value(ArithOp::Rem, char('%')),
    ))
    .parse(i)
}

fn unary_expr(i: &str) -> IResult<&str, Expr, NomError<&str>> {
    map(
        pair(
            map(
                many0(alt((value(true, char('-')), value(false, char('+'))))),
                |marks: Vec<bool>| marks.into_iter().fold(false, |acc, m| acc ^ m),
            ),
            base_expr,
        ),
        |(negate, base)| {
            if negate {
                Expr::Neg(Box::new(base))
            } else {
                base
            }
        },
    )
    .parse(i)
}

fn base_expr(i: &str) -> IResult<&str, Expr, NomError<&str>> {
    alt((
        map(nom_i64, Expr::Lit),
        map(preceded(char('?'), parse_id), Expr::Var),
        map(
            delimited(
                char('('),
                delimited(multispace0, add_expr, multispace0),
                char(')'),
            ),
            |a| a,
        ),
    ))
    .parse(i)
}

pub(crate) fn parse_id(i: &str) -> IResult<&str, String, NomError<&str>> {
    map(
        recognize(pair(
            satisfy(|c: char| c.is_ascii_alphabetic()),
            many0(satisfy(|c: char| c.is_ascii_alphanumeric() || c == '_')),
        )),
        |s: &str| s.to_string(),
    )
    .parse(i)
}

pub(crate) fn op_char<'a>(
    c: char,
) -> impl Parser<&'a str, Output = char, Error = NomError<&'a str>> {
    delimited(multispace0, char(c), multispace0)
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::star("*", ValueAst::Wildcard)]
    #[case::num("0", ValueAst::Lit(0))]
    #[case::set("{0,1,2}", ValueAst::LitSet(vec![0, 1, 2]))]
    #[case::set_spaced("{ 0, 1 ,2}", ValueAst::LitSet(vec![0, 1, 2]))]
    #[case::sum("1+2", ValueAst::Expr(Expr::BinOp(Box::new(Expr::Lit(1)), ArithOp::Add, Box::new(Expr::Lit(2)))))]
    #[case::sum_spaced("1 + 2", ValueAst::Expr(Expr::BinOp(Box::new(Expr::Lit(1)), ArithOp::Add, Box::new(Expr::Lit(2)))))]
    #[case::diff("1-2", ValueAst::Expr(Expr::BinOp(Box::new(Expr::Lit(1)), ArithOp::Sub, Box::new(Expr::Lit(2)))))]
    #[case::diff_spaced("1 - 2", ValueAst::Expr(Expr::BinOp(Box::new(Expr::Lit(1)), ArithOp::Sub, Box::new(Expr::Lit(2)))))]
    #[case::mult("1*2", ValueAst::Expr(Expr::BinOp(Box::new(Expr::Lit(1)), ArithOp::Mul, Box::new(Expr::Lit(2)))))]
    #[case::mult_spaced("1 * 2", ValueAst::Expr(Expr::BinOp(Box::new(Expr::Lit(1)), ArithOp::Mul, Box::new(Expr::Lit(2)))))]
    #[case::div("2/2", ValueAst::Expr(Expr::BinOp(Box::new(Expr::Lit(2)), ArithOp::Div, Box::new(Expr::Lit(2)))))]
    #[case::div_spaced("2 / 2", ValueAst::Expr(Expr::BinOp(Box::new(Expr::Lit(2)), ArithOp::Div, Box::new(Expr::Lit(2)))))]
    #[case::var("?h", ValueAst::Expr(Expr::Var("h".to_string())))]
    #[case::var_2char("?ha", ValueAst::Expr(Expr::Var("ha".to_string())))]
    #[case::var_number("?h1", ValueAst::Expr(Expr::Var("h1".to_string())))]
    #[case::var_underscore("?h_", ValueAst::Expr(Expr::Var("h_".to_string())))]
    #[case::membership("?h + 0 :: {0,1}", ValueAst::Expr(Expr::Mem(Box::new(Expr::BinOp(Box::new(Expr::Var("h".to_string())), ArithOp::Add, Box::new(Expr::Lit(0)))), vec![0, 1])))]
    #[case::double_neg("--0", ValueAst::Expr(Expr::Lit(0)))]
    #[case::not_and_precedence("! ?h == 0 & ?v == 1", ValueAst::Expr(Expr::And(vec![
        Expr::Not(Box::new(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Eq, Box::new(Expr::Lit(0))))),
        Expr::Rel(Box::new(Expr::Var("v".to_string())), RelOp::Eq, Box::new(Expr::Lit(1))),
    ])))]
    #[case::paren_arith("(0 + 1) * 1", ValueAst::Expr(Expr::BinOp(Box::new(Expr::BinOp(Box::new(Expr::Lit(0)), ArithOp::Add, Box::new(Expr::Lit(1)))), ArithOp::Mul, Box::new(Expr::Lit(1)))))]
    fn test_value_dsl(#[case] input: &str, #[case] expected: ValueAst) {
        let result = value_dsl.parse(input);
        assert!(result.is_ok(), "{:?} should have succeeded, error: {:?}", input, result.clone().unwrap_err());
        let (remaining, value) = result.unwrap();
        assert!(remaining.is_empty(), "{:?} should have consumed all input, remaining: {:?}", input, remaining);
        assert_eq!(value, expected);
    }

    #[rstest]
    #[case::empty("")]
    #[case::invalid_char("[]")]
    #[case::bare_open_paren("(")]
    #[case::bare_close_paren(")")]
    #[case::whitespace_id("? h")]
    #[case::adjacent_ops("a + * 3")]
    #[case::bare_plus("+")]
    #[case::bare_minus("-")]
    #[case::bare_equal("=")]
    #[case::bare_lt("<")]
    #[case::bare_gt(">")]
    #[case::leading_op("/ 3")]
    #[case::missing_id("? ")]
    #[case::invalid_id_1("?&x ")]
    #[case::invalid_id_2("?_x ")]
    #[case::triple_q("???")]
    #[case::empty_set("{}")]
    #[case::unclosed_paren_add("(0 + 1")]
    fn test_value_dsl_error(#[case] input: &str) {
        let res = value_dsl.parse(input);
        assert!(res.is_err(), "{input:?} should fail, got {:?}", res.unwrap());
        assert!(
            matches!(&res, Err(Err::Error(e)) if e.code == ErrorKind::Char),
            "{input:?} should fail with ErrorKind::Char, got {:?}",
            res.unwrap_err().map(|e| e.code)
        );
    }
}
