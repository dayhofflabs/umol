//! `value-dsl` — `spec/umol-dsl-spec.md` §5.

use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::{char, multispace0, satisfy};
use nom::combinator::{all_consuming, map, opt, recognize, value};
use nom::multi::{many0, separated_list1};
use nom::sequence::{delimited, pair, preceded, terminated};
use nom::{Err, IResult, Parser};

use super::error::ParseError;
use super::utils::IntParser;

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ValueAst<T: IntParser> {
    Wildcard,
    LitSet(Vec<T>),
    Lit(T),
    Expr(Expr<T>),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Expr<T: IntParser> {
    Lit(T),
    Var(String),
    Neg(Box<Expr<T>>),
    BinOp(Box<Expr<T>>, ArithOp, Box<Expr<T>>),
    Mem(Box<Expr<T>>, Vec<T>),
    Rel(Box<Expr<T>>, RelOp, Box<Expr<T>>),
    Not(Box<Expr<T>>),
    And(Vec<Expr<T>>),
    Or(Vec<Expr<T>>),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ArithOp {
    Add,
    Sub,
    Mul,
    Div,
    Rem,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RelOp {
    Le,
    Ge,
    Eq,
    Lt,
    Gt,
}

pub fn parse_value_dsl<T: IntParser>(input: &str) -> Result<ValueAst<T>, ParseError> {
    all_consuming(value_dsl::<T>)
        .parse(input)
        .map(|(_, v)| v)
        .map_err(|e| match e {
            Err::Error(e) | Err::Failure(e) => e,
            Err::Incomplete(_) => ParseError::Incomplete,
        })
}

pub fn value_dsl<T: IntParser>(i: &str) -> IResult<&str, ValueAst<T>, ParseError> {
    alt((
        map(
            terminated(T::nom_parser(), (multispace0, terminator)),
            ValueAst::Lit,
        ),
        value(ValueAst::Wildcard, tag("*")),
        map(lit_set::<T>, ValueAst::LitSet),
        map(bool_expr::<T>, ValueAst::Expr),
    ))
    .parse(i)
    .map_err(|_| Err::Error(ParseError::InvalidValue(i.to_string())))
}

fn terminator(i: &str) -> IResult<&str, (), ParseError> {
    if i.is_empty() || i.starts_with('#') {
        Ok((i, ()))
    } else {
        Err(Err::Error(ParseError::InvalidValue(i.to_string())))
    }
}

fn bool_expr<T: IntParser>(i: &str) -> IResult<&str, Expr<T>, ParseError> {
    map(
        pair(and_expr::<T>, many0(preceded(op_char('|'), and_expr::<T>))),
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

fn and_expr<T: IntParser>(i: &str) -> IResult<&str, Expr<T>, ParseError> {
    map(
        pair(not_expr::<T>, many0(preceded(op_char('&'), not_expr::<T>))),
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

fn not_expr<T: IntParser>(i: &str) -> IResult<&str, Expr<T>, ParseError> {
    alt((
        map(preceded((char('!'), multispace0), not_expr), |n| {
            Expr::Not(Box::new(n))
        }),
        rel_expr,
        map(
            delimited(
                char('('),
                delimited(multispace0, bool_expr::<T>, multispace0),
                char(')'),
            ),
            |b| b,
        ),
    ))
    .parse(i)
}

fn rel_expr<T: IntParser>(i: &str) -> IResult<&str, Expr<T>, ParseError> {
    map(
        pair(
            mem_expr::<T>,
            opt(preceded(
                multispace0,
                pair(rel_op, preceded(multispace0, mem_expr::<T>)),
            )),
        ),
        |(left, right)| match right {
            None => left,
            Some((op, r)) => Expr::Rel(Box::new(left), op, Box::new(r)),
        },
    )
    .parse(i)
}

fn mem_expr<T: IntParser>(i: &str) -> IResult<&str, Expr<T>, ParseError> {
    map(
        pair(
            add_expr::<T>,
            opt(preceded(
                multispace0,
                preceded(map(tag("::"), |_| ()), preceded(multispace0, lit_set::<T>)),
            )),
        ),
        |(expr, set)| match set {
            None => expr,
            Some(s) => Expr::Mem(Box::new(expr), s),
        },
    )
    .parse(i)
}

pub(crate) fn lit_set<T: IntParser>(i: &str) -> IResult<&str, Vec<T>, ParseError> {
    delimited(
        char('{'),
        delimited(
            multispace0,
            separated_list1(op_char(','), T::nom_parser()),
            multispace0,
        ),
        char('}'),
    )
    .parse(i)
}

fn rel_op(i: &str) -> IResult<&str, RelOp, ParseError> {
    alt((
        value(RelOp::Le, tag("<=")),
        value(RelOp::Ge, tag(">=")),
        value(RelOp::Eq, tag("==")),
        value(RelOp::Lt, char('<')),
        value(RelOp::Gt, char('>')),
    ))
    .parse(i)
}

fn add_expr<T: IntParser>(i: &str) -> IResult<&str, Expr<T>, ParseError> {
    map(
        pair(
            mult_expr::<T>,
            many0(pair(
                delimited(multispace0, add_op, multispace0),
                mult_expr::<T>,
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

fn add_op(i: &str) -> IResult<&str, ArithOp, ParseError> {
    alt((
        value(ArithOp::Add, char('+')),
        value(ArithOp::Sub, char('-')),
    ))
    .parse(i)
}

fn mult_expr<T: IntParser>(i: &str) -> IResult<&str, Expr<T>, ParseError> {
    map(
        pair(
            unary_expr::<T>,
            many0(pair(
                delimited(multispace0, mult_op, multispace0),
                unary_expr::<T>,
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

fn mult_op(i: &str) -> IResult<&str, ArithOp, ParseError> {
    alt((
        value(ArithOp::Mul, char('*')),
        value(ArithOp::Div, char('/')),
        value(ArithOp::Rem, char('%')),
    ))
    .parse(i)
}

fn unary_expr<T: IntParser>(i: &str) -> IResult<&str, Expr<T>, ParseError> {
    map(
        pair(
            map(
                many0(alt((value(true, char('-')), value(false, char('+'))))),
                |marks: Vec<bool>| marks.into_iter().fold(false, |acc, m| acc ^ m),
            ),
            base_expr::<T>,
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

fn base_expr<T: IntParser>(i: &str) -> IResult<&str, Expr<T>, ParseError> {
    alt((
        map(T::nom_parser(), Expr::Lit),
        map(preceded(char('?'), parse_id), Expr::Var),
        map(
            delimited(
                char('('),
                delimited(multispace0, add_expr::<T>, multispace0),
                char(')'),
            ),
            |a| a,
        ),
    ))
    .parse(i)
}

pub(crate) fn parse_id(i: &str) -> IResult<&str, String, ParseError> {
    map(
        recognize(pair(
            satisfy(|c: char| c.is_ascii_alphabetic()),
            many0(satisfy(|c: char| c.is_ascii_alphanumeric() || c == '_')),
        )),
        |s: &str| s.to_string(),
    )
    .parse(i)
}

pub(crate) fn op_char<'a>(c: char) -> impl Parser<&'a str, Output = char, Error = ParseError> {
    delimited(multispace0, char(c), multispace0)
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    fn elit(n: u8) -> Expr<u8> {
        Expr::Lit(n)
    }

    fn evar(s: &str) -> Expr<u8> {
        Expr::Var(s.into())
    }

    fn enot(e: Expr<u8>) -> Expr<u8> {
        Expr::Not(Box::new(e))
    }

    fn eadd(l: Expr<u8>, r: Expr<u8>) -> Expr<u8> {
        Expr::BinOp(Box::new(l), ArithOp::Add, Box::new(r))
    }

    fn esub(l: Expr<u8>, r: Expr<u8>) -> Expr<u8> {
        Expr::BinOp(Box::new(l), ArithOp::Sub, Box::new(r))
    }

    fn emul(l: Expr<u8>, r: Expr<u8>) -> Expr<u8> {
        Expr::BinOp(Box::new(l), ArithOp::Mul, Box::new(r))
    }

    fn ediv(l: Expr<u8>, r: Expr<u8>) -> Expr<u8> {
        Expr::BinOp(Box::new(l), ArithOp::Div, Box::new(r))
    }

    fn emem(e: Expr<u8>, set: Vec<u8>) -> Expr<u8> {
        Expr::Mem(Box::new(e), set)
    }

    fn erel(l: Expr<u8>, op: RelOp, r: Expr<u8>) -> Expr<u8> {
        Expr::Rel(Box::new(l), op, Box::new(r))
    }

    fn vexpr(e: Expr<u8>) -> ValueAst<u8> {
        ValueAst::Expr(e)
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::star("*", ValueAst::Wildcard)]
    #[case::num("0", ValueAst::Lit(0))]
    #[case::set("{0,1,2}", ValueAst::LitSet(vec![0, 1, 2]))]
    #[case::set_spaced("{ 0, 1 ,2}", ValueAst::LitSet(vec![0, 1, 2]))]
    #[case::sum("1+2", vexpr(eadd(elit(1), elit(2))))]
    #[case::sum_spaced("1 + 2", vexpr(eadd(elit(1), elit(2))))]
    #[case::diff("1-2", vexpr(esub(elit(1), elit(2))))]
    #[case::diff_spaced("1 - 2", vexpr(esub(elit(1), elit(2))))]
    #[case::mult("1*2", vexpr(emul(elit(1), elit(2))))]
    #[case::mult_spaced("1 * 2", vexpr(emul(elit(1), elit(2))))]
    #[case::div("2/2", vexpr(ediv(elit(2), elit(2))))]
    #[case::div_spaced("2 / 2", vexpr(ediv(elit(2), elit(2))))]
    #[case::var("?h", vexpr(evar("h")))]
    #[case::var_2char("?ha", vexpr(evar("ha")))]
    #[case::var_number("?h1", vexpr(evar("h1")))]
    #[case::var_underscore("?h_", vexpr(evar("h_")))]
    #[case::membership("?h + 0 :: {0,1}", vexpr(emem(eadd(evar("h"), elit(0)), vec![0, 1])))]
    #[case::double_neg("--0", vexpr(elit(0)))]
    #[case::not_and_precedence("! ?h == 0 & ?v == 1", vexpr(Expr::And(vec![
        enot(erel(evar("h"), RelOp::Eq, elit(0))),
        erel(evar("v"), RelOp::Eq, elit(1)),
    ])))]
    #[case::paren_arith("(0 + 1) * 1", vexpr(emul(eadd(elit(0), elit(1)), elit(1))))]
    fn test_value_dsl(#[case] input: &str, #[case] expected: ValueAst<u8>) {
        let result = value_dsl(input);
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
        let result = value_dsl::<u8>(input);
        assert!(result.is_err(), "{input:?} should fail, got {result:?}");
        let err = match result.unwrap_err() {
            Err::Error(e) | Err::Failure(e) => e,
            Err::Incomplete(_) => ParseError::Incomplete,
        };
        assert_eq!(
            err,
            ParseError::InvalidValue(input.to_string()),
            "{:?} should fail with InvalidValue, got {:?}",
            input,
            err
        );
    }
}
