//! Value DSL: parser and display helpers.

use winnow::ascii::{dec_int, multispace0};
use winnow::combinator::{alt, delimited, opt, preceded, repeat, separated, terminated};
use winnow::error::ErrMode;
use winnow::token::one_of;
use winnow::Parser;

use super::error::{PResult, ParseError};
use crate::ast::value::{ArithOp, Expr, RelOp, ValueAst};

pub fn parse_value(input: &str) -> Result<ValueAst, ParseError> {
    value.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn value(i: &mut &str) -> PResult<ValueAst> {
    alt((
        terminated(dec_int, (multispace0, terminator)).map(ValueAst::Lit),
        "*".value(ValueAst::Undetermined),
        lit_set.map(ValueAst::LitSet),
        bool_expr.map(ValueAst::Expr),
    ))
    .parse_next(i)
}

pub(crate) fn id(i: &mut &str) -> PResult<String> {
    (
        one_of(|c: char| c.is_ascii_alphabetic()),
        repeat::<_, _, (), _, _>(0.., one_of(|c: char| c.is_ascii_alphanumeric() || c == '_')),
    )
        .take()
        .map(|s: &str| s.to_string())
        .parse_next(i)
}

fn terminator(i: &mut &str) -> PResult<()> {
    if i.is_empty() || i.starts_with('#') {
        Ok(())
    } else {
        Err(ErrMode::Backtrack(ParseError::Syntax))
    }
}

fn bool_expr(i: &mut &str) -> PResult<Expr> {
    let first = and_expr.parse_next(i)?;
    let rest: Vec<Expr> = repeat(
        0..,
        preceded(delimited(multispace0, '|', multispace0), and_expr),
    )
    .parse_next(i)?;
    Ok(if rest.is_empty() {
        first
    } else {
        let mut disjuncts = vec![first];
        disjuncts.extend(rest);
        Expr::Or(disjuncts)
    })
}

fn and_expr(i: &mut &str) -> PResult<Expr> {
    let first = not_expr.parse_next(i)?;
    let rest: Vec<Expr> = repeat(
        0..,
        preceded(delimited(multispace0, '&', multispace0), not_expr),
    )
    .parse_next(i)?;
    Ok(if rest.is_empty() {
        first
    } else {
        let mut conjuncts = vec![first];
        conjuncts.extend(rest);
        Expr::And(conjuncts)
    })
}

fn not_expr(i: &mut &str) -> PResult<Expr> {
    alt((
        preceded(('!', multispace0), not_expr).map(|n| Expr::Not(Box::new(n))),
        rel_expr,
        delimited('(', delimited(multispace0, bool_expr, multispace0), ')'),
    ))
    .parse_next(i)
}

fn rel_expr(i: &mut &str) -> PResult<Expr> {
    let left = mem_expr.parse_next(i)?;
    let right = opt(preceded(
        multispace0,
        (rel_op, preceded(multispace0, mem_expr)),
    ))
    .parse_next(i)?;
    Ok(match right {
        None => left,
        Some((op, r)) => Expr::Rel(Box::new(left), op, Box::new(r)),
    })
}

fn mem_expr(i: &mut &str) -> PResult<Expr> {
    let expr = add_expr.parse_next(i)?;
    let set = opt(preceded(
        multispace0,
        preceded("::", preceded(multispace0, lit_set)),
    ))
    .parse_next(i)?;
    Ok(match set {
        None => expr,
        Some(s) => Expr::Mem(Box::new(expr), s),
    })
}

pub(crate) fn lit_set(i: &mut &str) -> PResult<Vec<i64>> {
    delimited(
        '{',
        delimited(
            multispace0,
            separated(
                1..,
                dec_int::<_, i64, _>,
                delimited(multispace0, ',', multispace0),
            ),
            multispace0,
        ),
        '}',
    )
    .parse_next(i)
}

fn rel_op(i: &mut &str) -> PResult<RelOp> {
    alt((
        "<=".value(RelOp::Le),
        ">=".value(RelOp::Ge),
        "==".value(RelOp::Eq),
        '<'.value(RelOp::Lt),
        '>'.value(RelOp::Gt),
    ))
    .parse_next(i)
}

fn add_expr(i: &mut &str) -> PResult<Expr> {
    let head = mult_expr.parse_next(i)?;
    let tail: Vec<(ArithOp, Expr)> = repeat(
        0..,
        (delimited(multispace0, add_op, multispace0), mult_expr),
    )
    .parse_next(i)?;
    Ok(tail.into_iter().fold(head, |acc, (op, rhs)| {
        Expr::BinOp(Box::new(acc), op, Box::new(rhs))
    }))
}

fn add_op(i: &mut &str) -> PResult<ArithOp> {
    alt(('+'.value(ArithOp::Add), '-'.value(ArithOp::Sub))).parse_next(i)
}

fn mult_expr(i: &mut &str) -> PResult<Expr> {
    let head = unary_expr.parse_next(i)?;
    let tail: Vec<(ArithOp, Expr)> = repeat(
        0..,
        (delimited(multispace0, mult_op, multispace0), unary_expr),
    )
    .parse_next(i)?;
    Ok(tail.into_iter().fold(head, |acc, (op, rhs)| {
        Expr::BinOp(Box::new(acc), op, Box::new(rhs))
    }))
}

fn mult_op(i: &mut &str) -> PResult<ArithOp> {
    alt((
        '*'.value(ArithOp::Mul),
        '/'.value(ArithOp::Div),
        '%'.value(ArithOp::Rem),
    ))
    .parse_next(i)
}

fn unary_expr(i: &mut &str) -> PResult<Expr> {
    let marks: Vec<bool> = repeat(0.., alt(('-'.value(true), '+'.value(false)))).parse_next(i)?;
    let negate = marks.into_iter().fold(false, |acc, m| acc ^ m);
    let base = base_expr.parse_next(i)?;
    Ok(if negate {
        Expr::Neg(Box::new(base))
    } else {
        base
    })
}

fn base_expr(i: &mut &str) -> PResult<Expr> {
    alt((
        dec_int::<_, i64, _>.map(Expr::Lit),
        preceded('?', id).map(Expr::Var),
        delimited('(', delimited(multispace0, add_expr, multispace0), ')'),
    ))
    .parse_next(i)
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::star("*", ValueAst::Undetermined)]
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
    fn test_value(#[case] input: &str, #[case] expected: ValueAst) {
        let result = value.parse(input);
        assert!(result.is_ok(), "{:?} should have succeeded, error: {:?}", input, result.clone().unwrap_err());
        let value = result.unwrap();
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
    fn test_value_error(#[case] input: &str) {
        let res = value.parse(input);
        assert!(
            res.is_err(),
            "{input:?} should fail, got {:?}",
            res.unwrap()
        );
    }
}
