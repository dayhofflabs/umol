//! `value-expr` — `spec/umol-dsl-spec.md` §5.

use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::{char, multispace0, satisfy};
use nom::combinator::{all_consuming, map, opt, recognize, value};
use nom::error::{Error as NomError, ErrorKind};
use nom::multi::{many0, separated_list1};
use nom::sequence::{delimited, pair, preceded};
use nom::{Err, IResult, Parser};

use super::error::ParseError;
use super::utils::IntParser;

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum ValueExpr<T: IntParser> {
    Wildcard,
    NatSet(Vec<T>),
    Lit(T),
    Bool(BoolExpr<T>),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct BoolExpr<T: IntParser> {
    pub disjuncts: Vec<AndExpr<T>>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct AndExpr<T: IntParser> {
    pub conjuncts: Vec<NotExpr<T>>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum NotExpr<T: IntParser> {
    Not(Box<NotExpr<T>>),
    Block(Box<BoolExpr<T>>),
    Rel(RelExpr<T>),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RelExpr<T: IntParser> {
    pub left: MemExpr<T>,
    pub right: Option<(RelOp, MemExpr<T>)>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct MemExpr<T: IntParser> {
    pub expr: AddExpr<T>,
    pub set: Option<Vec<T>>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct AddExpr<T: IntParser> {
    pub head: MultExpr<T>,
    pub tail: Vec<(AddOp, MultExpr<T>)>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct MultExpr<T: IntParser> {
    pub head: UnaryExpr<T>,
    pub tail: Vec<(MultOp, UnaryExpr<T>)>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct UnaryExpr<T: IntParser> {
    pub negate: bool,
    pub base: BaseExpr<T>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum BaseExpr<T: IntParser> {
    Lit(T),
    Var(String),
    Paren(Box<AddExpr<T>>),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum RelOp {
    Le,
    Ge,
    Eq,
    Lt,
    Gt,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum AddOp {
    Add,
    Sub,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum MultOp {
    Mul,
    Div,
    Rem,
}

pub(crate) fn parse_value_expr<T: IntParser>(input: &str) -> Result<ValueExpr<T>, ParseError> {
    all_consuming(value_expr::<T>)
        .parse(input)
        .map(|(_, v)| v)
        .map_err(ParseError::value_from_nom)
}

pub(crate) fn value_expr<T: IntParser>(i: &str) -> IResult<&str, ValueExpr<T>> {
    alt((
        value(ValueExpr::Wildcard, tag("*")),
        map(nat_set::<T>, ValueExpr::NatSet),
        map(
            pair(T::nom_parser(), preceded(multispace0, terminator)),
            |(n, _)| ValueExpr::Lit(n),
        ),
        map(bool_expr::<T>, ValueExpr::Bool),
    ))
    .parse(i)
}

fn terminator(i: &str) -> IResult<&str, ()> {
    if i.is_empty() || i.starts_with('#') {
        Ok((i, ()))
    } else {
        Err(Err::Error(NomError::new(i, ErrorKind::Verify)))
    }
}

fn bool_expr<T: IntParser>(i: &str) -> IResult<&str, BoolExpr<T>> {
    map(separated_list1(ws_sym('|'), and_expr), |disjuncts| {
        BoolExpr { disjuncts }
    })
    .parse(i)
}

fn and_expr<T: IntParser>(i: &str) -> IResult<&str, AndExpr<T>> {
    map(separated_list1(ws_sym('&'), not_expr), |conjuncts| {
        AndExpr { conjuncts }
    })
    .parse(i)
}

fn not_expr<T: IntParser>(i: &str) -> IResult<&str, NotExpr<T>> {
    alt((
        map(preceded((char('!'), multispace0), not_expr), |n| {
            NotExpr::Not(Box::new(n))
        }),
        map(rel_expr, NotExpr::Rel),
        map(
            delimited(
                char('('),
                delimited(multispace0, bool_expr::<T>, multispace0),
                char(')'),
            ),
            |b| NotExpr::Block(Box::new(b)),
        ),
    ))
    .parse(i)
}

fn rel_expr<T: IntParser>(i: &str) -> IResult<&str, RelExpr<T>> {
    map(
        pair(
            mem_expr,
            opt(preceded(
                multispace0,
                pair(rel_op, preceded(multispace0, mem_expr)),
            )),
        ),
        |(left, right)| RelExpr { left, right },
    )
    .parse(i)
}

fn mem_expr<T: IntParser>(i: &str) -> IResult<&str, MemExpr<T>> {
    map(
        pair(
            add_expr,
            opt(preceded(
                multispace0,
                preceded(map(tag("::"), |_| ()), preceded(multispace0, nat_set)),
            )),
        ),
        |(expr, set)| MemExpr { expr, set },
    )
    .parse(i)
}

fn nat_set<T: IntParser>(i: &str) -> IResult<&str, Vec<T>> {
    delimited(
        char('{'),
        delimited(
            multispace0,
            separated_list1(ws_sym(','), T::nom_parser()),
            multispace0,
        ),
        char('}'),
    )
    .parse(i)
}

fn rel_op(i: &str) -> IResult<&str, RelOp> {
    alt((
        value(RelOp::Le, tag("<=")),
        value(RelOp::Ge, tag(">=")),
        value(RelOp::Eq, tag("==")),
        value(RelOp::Lt, char('<')),
        value(RelOp::Gt, char('>')),
    ))
    .parse(i)
}

fn add_expr<T: IntParser>(i: &str) -> IResult<&str, AddExpr<T>> {
    map(
        pair(
            mult_expr::<T>,
            many0(pair(
                delimited(multispace0, add_op, multispace0),
                mult_expr::<T>,
            )),
        ),
        |(head, tail)| AddExpr { head, tail },
    )
    .parse(i)
}

fn add_op(i: &str) -> IResult<&str, AddOp> {
    alt((value(AddOp::Add, char('+')), value(AddOp::Sub, char('-')))).parse(i)
}

fn mult_expr<T: IntParser>(i: &str) -> IResult<&str, MultExpr<T>> {
    map(
        pair(
            unary_expr::<T>,
            many0(pair(
                delimited(multispace0, mult_op, multispace0),
                unary_expr::<T>,
            )),
        ),
        |(head, tail)| MultExpr { head, tail },
    )
    .parse(i)
}

fn mult_op(i: &str) -> IResult<&str, MultOp> {
    alt((
        value(MultOp::Mul, char('*')),
        value(MultOp::Div, char('/')),
        value(MultOp::Rem, char('%')),
    ))
    .parse(i)
}

fn unary_expr<T: IntParser>(i: &str) -> IResult<&str, UnaryExpr<T>> {
    map(
        pair(
            map(
                many0(alt((value(true, char('-')), value(false, char('+'))))),
                |marks: Vec<bool>| marks.into_iter().fold(false, |acc, m| acc ^ m),
            ),
            base_expr::<T>,
        ),
        |(negate, base)| UnaryExpr { negate, base },
    )
    .parse(i)
}

fn base_expr<T: IntParser>(i: &str) -> IResult<&str, BaseExpr<T>> {
    alt((
        map(T::nom_parser(), BaseExpr::Lit),
        map(preceded(char('?'), parse_id), BaseExpr::Var),
        map(
            delimited(
                char('('),
                delimited(multispace0, add_expr::<T>, multispace0),
                char(')'),
            ),
            |a| BaseExpr::Paren(Box::new(a)),
        ),
    ))
    .parse(i)
}

fn parse_id(i: &str) -> IResult<&str, String> {
    map(
        recognize(pair(
            satisfy(|c: char| c.is_ascii_alphabetic()),
            many0(satisfy(|c: char| c.is_ascii_alphanumeric() || c == '_')),
        )),
        |s: &str| s.to_string(),
    )
    .parse(i)
}

fn ws_sym<'a>(c: char) -> impl Parser<&'a str, Output = char, Error = NomError<&'a str>> {
    delimited(multispace0, char(c), multispace0)
}

#[cfg(test)]
mod tests {
    use nom::error::ErrorKind as NomErrorKind;
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    fn lit(n: u8) -> UnaryExpr<u8> {
        UnaryExpr {
            negate: false,
            base: BaseExpr::Lit(n),
        }
    }

    fn var(name: &str) -> UnaryExpr<u8> {
        UnaryExpr {
            negate: false,
            base: BaseExpr::Var(name.into()),
        }
    }

    fn mult1(head: UnaryExpr<u8>, tail: Vec<(MultOp, UnaryExpr<u8>)>) -> MultExpr<u8> {
        MultExpr { head, tail }
    }

    fn add1(head: MultExpr<u8>, tail: Vec<(AddOp, MultExpr<u8>)>) -> AddExpr<u8> {
        AddExpr { head, tail }
    }

    fn mem(expr: AddExpr<u8>, set: Option<Vec<u8>>) -> MemExpr<u8> {
        MemExpr { expr, set }
    }

    fn ax_nat(n: u8) -> AddExpr<u8> {
        add1(mult1(lit(n), vec![]), vec![])
    }

    fn ax_var(s: &str) -> AddExpr<u8> {
        add1(mult1(var(s), vec![]), vec![])
    }

    fn ax_var_add(s: &str, n: u8) -> AddExpr<u8> {
        add1(
            mult1(var(s), vec![]),
            vec![(AddOp::Add, mult1(lit(n), vec![]))],
        )
    }

    fn ax_lit_add(a: u8, b: u8) -> AddExpr<u8> {
        add1(
            mult1(lit(a), vec![]),
            vec![(AddOp::Add, mult1(lit(b), vec![]))],
        )
    }

    fn ax_lit_sub(a: u8, b: u8) -> AddExpr<u8> {
        add1(
            mult1(lit(a), vec![]),
            vec![(AddOp::Sub, mult1(lit(b), vec![]))],
        )
    }

    fn ax_lit_mul(a: u8, b: u8) -> AddExpr<u8> {
        add1(mult1(lit(a), vec![(MultOp::Mul, lit(b))]), vec![])
    }

    fn ax_lit_div(a: u8, b: u8) -> AddExpr<u8> {
        add1(mult1(lit(a), vec![(MultOp::Div, lit(b))]), vec![])
    }

    fn ax_paren_times(inner: AddExpr<u8>, k: u8) -> AddExpr<u8> {
        add1(
            mult1(
                UnaryExpr {
                    negate: false,
                    base: BaseExpr::Paren(Box::new(inner)),
                },
                vec![(MultOp::Mul, lit(k))],
            ),
            vec![],
        )
    }

    fn mem_nat(n: u8) -> MemExpr<u8> {
        mem(ax_nat(n), None)
    }

    fn mem_var(s: &str) -> MemExpr<u8> {
        mem(ax_var(s), None)
    }

    fn rel_only(left: MemExpr<u8>) -> RelExpr<u8> {
        RelExpr { left, right: None }
    }

    fn rel_eq(left: MemExpr<u8>, r: MemExpr<u8>) -> RelExpr<u8> {
        RelExpr {
            left,
            right: Some((RelOp::Eq, r)),
        }
    }

    fn vx_rel_single(rel: RelExpr<u8>) -> ValueExpr<u8> {
        ValueExpr::Bool(BoolExpr {
            disjuncts: vec![AndExpr {
                conjuncts: vec![NotExpr::Rel(rel)],
            }],
        })
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::star("*", ValueExpr::Wildcard)]
    #[case::num("0", ValueExpr::Lit(0))]
    #[case::set("{0,1,2}", ValueExpr::NatSet(vec![0, 1, 2]))]
    #[case::set_spaced("{ 0, 1 ,2}", ValueExpr::NatSet(vec![0, 1, 2]))]
    #[case::sum("1+2", vx_rel_single(rel_only(mem(ax_lit_add(1, 2), None))))]
    #[case::sum_spaced("1 + 2", vx_rel_single(rel_only(mem(ax_lit_add(1, 2), None))))]
    #[case::diff("1-2", vx_rel_single(rel_only(mem(ax_lit_sub(1, 2), None))))]
    #[case::diff_spaced("1 - 2", vx_rel_single(rel_only(mem(ax_lit_sub(1, 2), None))))]
    #[case::mult("1*2", vx_rel_single(rel_only(mem(ax_lit_mul(1, 2), None))))]
    #[case::mult_spaced("1 * 2", vx_rel_single(rel_only(mem(ax_lit_mul(1, 2), None))))]
    #[case::div("2/2", vx_rel_single(rel_only(mem(ax_lit_div(2, 2), None))))]
    #[case::div_spaced("2 / 2", vx_rel_single(rel_only(mem(ax_lit_div(2, 2), None))))]
    #[case::var("?h", vx_rel_single(rel_only(mem_var("h"))))]
    #[case::var_2char("?ha", vx_rel_single(rel_only(mem_var("ha"))))]
    #[case::var_number("?h1", vx_rel_single(rel_only(mem_var("h1"))))]
    #[case::var_underscore("?h_", vx_rel_single(rel_only(mem_var("h_"))))]
    #[case::membership("?h + 0 :: {0,1}", vx_rel_single(rel_only(mem(ax_var_add("h", 0), Some(vec![0, 1])))))]
    #[case::double_neg("--0", vx_rel_single(rel_only(mem_nat(0))))]
    #[case::not_and_precedence("! ?h == 0 & ?v == 1", ValueExpr::Bool(BoolExpr { disjuncts: vec![AndExpr { conjuncts:
        vec![ NotExpr::Not(Box::new(NotExpr::Rel(rel_eq(mem_var("h"), mem_nat(0))))), NotExpr::Rel(rel_eq(mem_var("v"), mem_nat(1))) ] }] }))]
    #[case::paren_arith("(0 + 1) * 1", vx_rel_single(rel_only(mem(ax_paren_times(ax_lit_add(0, 1), 1), None))))]
    fn test_parse_value_expr(#[case] input: &str, #[case] expected: ValueExpr<u8>) {
        let result = value_expr(input);
        assert!(result.is_ok(), "{:?} should have succeeded, error: {:?}", input, result.clone().unwrap_err());
        let (remaining, value) = result.unwrap(); 
        assert!(remaining.is_empty(), "{:?} should have consumed all input, remaining: {:?}", input, remaining);
        assert_eq!(value, expected);
    }

    #[rstest]
    #[case::literal("1#v3", ValueExpr::Lit(1), "#v3")]
    fn test_value_expr_literal(
        #[case] input: &str,
        #[case] expected_value: ValueExpr<u8>,
        #[case] expected_remaining: &str,
    ) {
        let result = value_expr::<u8>.parse("1#v3");
        assert!(
            result.is_ok(),
            "{:?} should have succeeded, error: {:?}",
            input,
            result.clone().unwrap_err()
        );
        let (remaining, value) = result.unwrap();
        assert_eq!(value, expected_value);
        assert_eq!(remaining, expected_remaining);
    }

    #[rstest]
    #[case::empty("", NomErrorKind::Char)]
    #[case::invalid_char("[]", NomErrorKind::Char)]
    #[case::bare_open_paren("(", NomErrorKind::Char)]
    #[case::bare_close_paren(")", NomErrorKind::Char)]
    #[case::whitespace_id("? h", NomErrorKind::Char)]
    #[case::adjacent_ops("a + * 3", NomErrorKind::Char)]
    #[case::bare_plus("+", NomErrorKind::Char)]
    #[case::bare_minus("-", NomErrorKind::Char)]
    #[case::bare_equal("=", NomErrorKind::Char)]
    #[case::bare_lt("<", NomErrorKind::Char)]
    #[case::bare_gt("<", NomErrorKind::Char)]
    #[case::leading_op("/ 3", NomErrorKind::Char)]
    #[case::missing_id("? ", NomErrorKind::Char)]
    #[case::invalid_id_1("?&x ", NomErrorKind::Char)]
    #[case::invalid_id_2("?_x ", NomErrorKind::Char)]
    #[case::triple_q("???", NomErrorKind::Char)]
    #[case::empty_set("{}", NomErrorKind::Char)]
    #[case::unclosed_paren_add("(0 + 1", NomErrorKind::Char)]
    fn test_parse_value_expr_invalid(#[case] input: &str, #[case] expected_kind: NomErrorKind) {
        let result = value_expr::<u8>(input);
        assert!(
            result.is_err(),
            "{:?} should have failed, output: {:?}",
            input,
            result.clone().unwrap()
        );
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "{:?} should have failed with error kind {:?}, got {:?}",
            input,
            expected_kind,
            result.clone().unwrap_err().map(|e| e.code)
        );
    }
}
