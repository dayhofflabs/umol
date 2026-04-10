//! `value-dsl` — `spec/umol-dsl-spec.md` §5

use std::collections::HashMap;

use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::{char, i32 as nom_i32, multispace0, satisfy};
use nom::combinator::{all_consuming, map, opt, recognize, value};
use nom::error::{Error as NomError, ErrorKind};
use nom::multi::{many0, separated_list1};
use nom::sequence::{delimited, pair, preceded, terminated};
use nom::{Err, IResult, Parser};

use super::error::EvaluationError;

/// Variable bindings used by [`Expr::evaluate`] and [`Expr::evaluate_bool`]
pub type Bindings = HashMap<String, i32>;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum ValueAst {
    Wildcard,
    LitSet(Vec<i32>),
    Lit(i32),
    Expr(Expr),
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Expr {
    Lit(i32),
    Var(String),
    Neg(Box<Expr>),
    BinOp(Box<Expr>, ArithOp, Box<Expr>),
    Mem(Box<Expr>, Vec<i32>),
    Rel(Box<Expr>, RelOp, Box<Expr>),
    Not(Box<Expr>),
    And(Vec<Expr>),
    Or(Vec<Expr>),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ArithOp {
    Add,
    Sub,
    Mul,
    Div,
    Rem,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RelOp {
    Le,
    Ge,
    Eq,
    Lt,
    Gt,
}

impl ValueAst {
    /// Match a concrete integer value against this pattern
    pub fn matches(&self, value: i32) -> bool {
        self.capture(value).is_some()
    }

    /// Match a concrete integer value against this pattern, returning variable bindings
    ///
    /// Variables in the pattern are bound to `value`. For boolean expressions the
    /// predicate is evaluated with those bindings; for arithmetic expressions the
    /// result is compared to `value`
    pub fn capture(&self, value: i32) -> Option<Bindings> {
        match self {
            ValueAst::Wildcard => Some(Bindings::new()),
            ValueAst::Lit(n) => {
                if *n == value {
                    Some(Bindings::new())
                } else {
                    None
                }
            }
            ValueAst::LitSet(s) => {
                if s.contains(&value) {
                    Some(Bindings::new())
                } else {
                    None
                }
            }
            ValueAst::Expr(e) => {
                let mut bindings = Bindings::new();
                collect_bindings(e, value, &mut bindings);
                if e.is_arithmetic() {
                    match e.evaluate(&bindings) {
                        Ok(v) if v == value => Some(bindings),
                        _ => None,
                    }
                } else {
                    match e.evaluate_bool(&bindings) {
                        Ok(true) => Some(bindings),
                        _ => None,
                    }
                }
            }
        }
    }
}

/// Recursively bind every variable in `expr` to `value`
fn collect_bindings(expr: &Expr, value: i32, bindings: &mut Bindings) {
    match expr {
        Expr::Var(name) => {
            bindings.insert(name.clone(), value);
        }
        Expr::Neg(e) => collect_bindings(e, value, bindings),
        Expr::BinOp(l, _, r) => {
            collect_bindings(l, value, bindings);
            collect_bindings(r, value, bindings);
        }
        Expr::Mem(e, _) => collect_bindings(e, value, bindings),
        Expr::Rel(l, _, r) => {
            collect_bindings(l, value, bindings);
            collect_bindings(r, value, bindings);
        }
        Expr::Not(e) => collect_bindings(e, value, bindings),
        Expr::And(exprs) | Expr::Or(exprs) => {
            for e in exprs {
                collect_bindings(e, value, bindings);
            }
        }
        Expr::Lit(_) => {}
    }
}

impl Expr {
    fn is_arithmetic(&self) -> bool {
        matches!(
            self,
            Expr::Lit(..) | Expr::Var(..) | Expr::Neg(..) | Expr::BinOp(..)
        )
    }

    /// Evaluate an arithmetic expression to an `i32`
    ///
    /// Returns [`EvaluationError::TypeMismatch`] if called on a boolean-domain
    /// expression (`Rel`, `Mem`, `Not`, `And`, `Or`)
    pub fn evaluate(&self, vars: &Bindings) -> Result<i32, EvaluationError> {
        match self {
            Expr::Lit(n) => Ok(*n),
            Expr::Var(name) => vars
                .get(name)
                .copied()
                .ok_or_else(|| EvaluationError::UnboundVariable(name.clone())),
            Expr::Neg(e) => Ok(-e.evaluate(vars)?),
            Expr::BinOp(l, op, r) => {
                let l = l.evaluate(vars)?;
                let r = r.evaluate(vars)?;
                match op {
                    ArithOp::Add => Ok(l + r),
                    ArithOp::Sub => Ok(l - r),
                    ArithOp::Mul => Ok(l * r),
                    ArithOp::Div => {
                        if r == 0 {
                            Err(EvaluationError::DivisionByZero)
                        } else {
                            Ok(l / r)
                        }
                    }
                    ArithOp::Rem => {
                        if r == 0 {
                            Err(EvaluationError::DivisionByZero)
                        } else {
                            Ok(l % r)
                        }
                    }
                }
            }
            Expr::Rel(..) | Expr::Mem(..) | Expr::Not(..) | Expr::And(..) | Expr::Or(..) => {
                Err(EvaluationError::TypeMismatch)
            }
        }
    }

    /// Evaluate a boolean expression to a `bool`
    ///
    /// Returns [`EvaluationError::TypeMismatch`] if called on an arithmetic-domain
    /// expression (`Lit`, `Var`, `Neg`, `BinOp`)
    pub fn evaluate_bool(&self, vars: &Bindings) -> Result<bool, EvaluationError> {
        match self {
            Expr::Rel(l, op, r) => {
                let l = l.evaluate(vars)?;
                let r = r.evaluate(vars)?;
                Ok(match op {
                    RelOp::Le => l <= r,
                    RelOp::Ge => l >= r,
                    RelOp::Eq => l == r,
                    RelOp::Lt => l < r,
                    RelOp::Gt => l > r,
                })
            }
            Expr::Mem(e, set) => Ok(set.contains(&e.evaluate(vars)?)),
            Expr::Not(e) => Ok(!e.evaluate_bool(vars)?),
            Expr::And(exprs) => {
                for e in exprs {
                    if !e.evaluate_bool(vars)? {
                        return Ok(false);
                    }
                }
                Ok(true)
            }
            Expr::Or(exprs) => {
                for e in exprs {
                    if e.evaluate_bool(vars)? {
                        return Ok(true);
                    }
                }
                Ok(false)
            }
            Expr::Lit(..) | Expr::Var(..) | Expr::Neg(..) | Expr::BinOp(..) => {
                Err(EvaluationError::TypeMismatch)
            }
        }
    }
}

pub fn parse_value_dsl(input: &str) -> Result<ValueAst, Err<NomError<&str>>> {
    all_consuming(value_dsl).parse(input).map(|(_, v)| v)
}

pub fn value_dsl(i: &str) -> IResult<&str, ValueAst, NomError<&str>> {
    alt((
        map(
            terminated(nom_i32, (multispace0, terminator)),
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

pub(crate) fn lit_set(i: &str) -> IResult<&str, Vec<i32>, NomError<&str>> {
    delimited(
        char('{'),
        delimited(
            multispace0,
            separated_list1(op_char(','), nom_i32),
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
        map(nom_i32, Expr::Lit),
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
    use crate::dsl::error::EvaluationError;

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

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(Expr::Lit(5), Bindings::new(), 5)]
    #[case::var_bound(Expr::Var("x".to_string()), Bindings::from([("x".to_string(), 3)]), 3)]
    #[case::neg(Expr::Neg(Box::new(Expr::Lit(3))), Bindings::new(), -3)]
    #[case::add(Expr::BinOp(Box::new(Expr::Lit(2)), ArithOp::Add, Box::new(Expr::Lit(3))), Bindings::new(), 5)]
    #[case::sub(Expr::BinOp(Box::new(Expr::Lit(5)), ArithOp::Sub, Box::new(Expr::Lit(3))), Bindings::new(), 2)]
    #[case::mul(Expr::BinOp(Box::new(Expr::Lit(3)), ArithOp::Mul, Box::new(Expr::Lit(4))), Bindings::new(), 12)]
    #[case::div(Expr::BinOp(Box::new(Expr::Lit(10)), ArithOp::Div, Box::new(Expr::Lit(3))), Bindings::new(), 3)]
    #[case::rem(Expr::BinOp(Box::new(Expr::Lit(10)), ArithOp::Rem, Box::new(Expr::Lit(3))), Bindings::new(), 1)]
    fn test_evaluate(#[case] expr: Expr, #[case] vars: Bindings, #[case] expected: i32) {
        let result = expr.evaluate(&vars);
        assert!(result.is_ok(), "{:?} should have succeeded, error: {:?}", expr, result.clone().unwrap_err());
        assert_eq!(result.unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::var_unbound(Expr::Var("x".to_string()), Bindings::new(), EvaluationError::UnboundVariable("x".to_string()))]
    #[case::div_zero(Expr::BinOp(Box::new(Expr::Lit(10)), ArithOp::Div, Box::new(Expr::Lit(0))), Bindings::new(), EvaluationError::DivisionByZero)]
    #[case::rem_zero(Expr::BinOp(Box::new(Expr::Lit(10)), ArithOp::Rem, Box::new(Expr::Lit(0))), Bindings::new(), EvaluationError::DivisionByZero)]
    #[case::type_mismatch(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Eq, Box::new(Expr::Lit(1))), Bindings::new(), EvaluationError::TypeMismatch)]
    fn test_evaluate_invalid(#[case] expr: Expr, #[case] vars: Bindings, #[case] expected: EvaluationError) {
        let result = expr.evaluate(&vars);
        assert!(result.is_err(), "{:?} should have failed, error: {:?}", expr, result.clone().unwrap_err());
        assert_eq!(result.unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::rel_eq_true(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Eq, Box::new(Expr::Lit(1))), Bindings::new(), true)]
    #[case::rel_eq_false(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Eq, Box::new(Expr::Lit(2))), Bindings::new(), false)]
    #[case::rel_lt_true(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Lt, Box::new(Expr::Lit(2))), Bindings::new(), true)]
    #[case::rel_lt_false(Expr::Rel(Box::new(Expr::Lit(2)), RelOp::Lt, Box::new(Expr::Lit(1))), Bindings::new(), false)]
    #[case::rel_le_true(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Le, Box::new(Expr::Lit(1))), Bindings::new(), true)]
    #[case::rel_le_false(Expr::Rel(Box::new(Expr::Lit(2)), RelOp::Le, Box::new(Expr::Lit(1))), Bindings::new(), false)]
    #[case::rel_gt_true(Expr::Rel(Box::new(Expr::Lit(2)), RelOp::Gt, Box::new(Expr::Lit(1))), Bindings::new(), true)]
    #[case::rel_gt_false(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Gt, Box::new(Expr::Lit(2))), Bindings::new(), false)]
    #[case::rel_ge_true(Expr::Rel(Box::new(Expr::Lit(2)), RelOp::Ge, Box::new(Expr::Lit(2))), Bindings::new(), true)]
    #[case::rel_ge_false(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Ge, Box::new(Expr::Lit(2))), Bindings::new(), false)]
    #[case::mem_true(Expr::Mem(Box::new(Expr::Lit(2)), vec![1, 2, 3]), Bindings::new(), true)]
    #[case::mem_false(Expr::Mem(Box::new(Expr::Lit(4)), vec![1, 2, 3]), Bindings::new(), false)]
    #[case::not_true(Expr::Not(Box::new(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Eq, Box::new(Expr::Lit(2))))), Bindings::new(), true)]
    #[case::not_false(Expr::Not(Box::new(Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Eq, Box::new(Expr::Lit(1))))), Bindings::new(), false)]
    #[case::and_true(Expr::And(vec![Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Lt, Box::new(Expr::Lit(2))), Expr::Rel(Box::new(Expr::Lit(3)), RelOp::Gt, Box::new(Expr::Lit(2)))]), Bindings::new(), true)]
    #[case::and_false(Expr::And(vec![Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Lt, Box::new(Expr::Lit(2))), Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Gt, Box::new(Expr::Lit(2)))]), Bindings::new(), false)]
    #[case::or_true(Expr::Or(vec![Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Eq, Box::new(Expr::Lit(2))), Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Lt, Box::new(Expr::Lit(2)))]), Bindings::new(), true)]
    #[case::or_false(Expr::Or(vec![Expr::Rel(Box::new(Expr::Lit(1)), RelOp::Eq, Box::new(Expr::Lit(2))), Expr::Rel(Box::new(Expr::Lit(3)), RelOp::Lt, Box::new(Expr::Lit(2)))]), Bindings::new(), false)]
    #[case::var_in_rel(Expr::Rel(Box::new(Expr::Var("x".to_string())), RelOp::Gt, Box::new(Expr::Lit(0))), Bindings::from([("x".to_string(), 5)]), true)]
    fn test_evaluate_bool(#[case] expr: Expr, #[case] vars: Bindings, #[case] expected: bool) {
        let result = expr.evaluate_bool(&vars);
        assert!(result.is_ok(), "{:?} should have succeeded, error: {:?}", expr, result.clone().unwrap_err());
        assert_eq!(result.unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::unbound_in_rel(Expr::Rel(Box::new(Expr::Var("x".to_string())), RelOp::Gt, Box::new(Expr::Lit(0))), Bindings::new(), EvaluationError::UnboundVariable("x".to_string()))]
    #[case::type_mismatch(Expr::Lit(1), Bindings::new(), EvaluationError::TypeMismatch)]
    fn test_evaluate_bool_invalid(#[case] expr: Expr, #[case] vars: Bindings, #[case] expected: EvaluationError) {
        let result = expr.evaluate_bool(&vars);
        assert!(result.is_err(), "{:?} should have failed, error: {:?}", expr, result.clone().unwrap_err());
        assert_eq!(result.unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(ValueAst::Wildcard, 3, true)]
    #[case::lit_match(ValueAst::Lit(3), 3,  true)]
    #[case::lit_set_match(ValueAst::LitSet(vec![1, 2, 3]), 2, true)]
    #[case::expr_var(ValueAst::Expr(Expr::Var("h".to_string())), 5, true)]
    #[case::expr_lit_match(ValueAst::Expr(Expr::Lit(3)), 3, true)]
    #[case::expr_rel_match(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))), 3, true)]
    #[case::expr_mem_match(ValueAst::Expr(Expr::Mem(Box::new(Expr::Var("h".to_string())), vec![0, 1])), 1, true)]
    #[case::lit_no_match(ValueAst::Lit(3), 4, false)]
    #[case::expr_lit_no_match(ValueAst::Expr(Expr::Lit(3)), 4, false)]
    #[case::expr_rel_no_match(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))), 0, false)]
    fn test_matches(#[case] pattern: ValueAst, #[case] value: i32, #[case] expected: bool) {
        assert_eq!(pattern.matches(value), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(ValueAst::Wildcard, 3, Bindings::new())]
    #[case::lit_match(ValueAst::Lit(3), 3, Bindings::new())]
    #[case::lit_set_match(ValueAst::LitSet(vec![1, 2, 3]), 2, Bindings::new())]
    #[case::expr_var(ValueAst::Expr(Expr::Var("h".to_string())), 5, Bindings::from([("h".to_string(), 5)]))]
    #[case::expr_lit_match(ValueAst::Expr(Expr::Lit(3)), 3, Bindings::new())]
    #[case::expr_rel_match(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))), 3, Bindings::from([("h".to_string(), 3)]))]
    #[case::expr_mem_match(ValueAst::Expr(Expr::Mem(Box::new(Expr::Var("h".to_string())), vec![0, 1])), 1, Bindings::from([("h".to_string(), 1)]))]
    fn test_capture(#[case] pattern: ValueAst, #[case] value: i32, #[case] expected: Bindings) {
        assert_eq!(pattern.capture(value), Some(expected));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit_no_match(ValueAst::Lit(3), 4)]
    #[case::expr_lit_no_match(ValueAst::Expr(Expr::Lit(3)), 4)]
    #[case::expr_rel_no_match(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))), 0)]
    fn test_capture_no_match(#[case] pattern: ValueAst, #[case] value: i32) {
        assert_eq!(pattern.capture(value), None);
    }
}
