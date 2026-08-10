//! Numeric DSL: parser, `Display`, and EDN boundary.

use std::borrow::Cow;
use std::collections::BTreeSet;
use std::fmt::{self, Display, Write};

use umol_edn::{DeError, Edn, EdnKeyword, FromEdn, ToEdn};
use winnow::ascii::{dec_int, digit1, multispace0};
use winnow::combinator::{alt, delimited, opt, preceded, repeat, separated, terminated};
use winnow::error::ErrMode;
use winnow::token::one_of;
use winnow::Parser;

use super::error::{PResult, ParseError};
use super::operators::{mem_op, mem_op_str, rel_op, rel_op_str};
use crate::ir::num::{ArithExpr, NumForm, PredExpr};
use crate::ir::traits::{FromIr, IntoIr};

/// Surface DSL wrapper around `NumForm`. EDN form is hybrid: `Lit` → `Int`,
/// `Undetermined` → `:undetermined`, `LitSet` → vector of ints, `ArithExpr`/
/// `PredExpr` → string via the numeric subgrammar (EDN has no native form for
/// the arithmetic/boolean grammar and round-trip fidelity is mandatory).
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct NumDsl(pub NumForm);

impl FromIr<NumForm> for NumDsl {
    type Ctx = ();

    fn from_ir(form: &NumForm, _ctx: &Self::Ctx) -> Self {
        Self(form.clone())
    }
}

impl IntoIr<NumForm> for NumDsl {
    type Ctx = ();

    fn into_ir(self, _ctx: &Self::Ctx) -> NumForm {
        self.0
    }
}

impl Display for NumDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_num(f, &self.0)
    }
}

impl<'de> FromEdn<'de> for NumDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let v = match edn {
            Edn::Int(n) => NumForm::Lit(*n),
            Edn::Keyword(k) if k.name() == "undetermined" => NumForm::Undetermined,
            Edn::Vector(xs) => {
                let mut out = BTreeSet::new();
                for e in xs.iter() {
                    let Edn::Int(n) = e else {
                        return Err(DeError::TypeMismatch {
                            expected: "int (numeric-set element)",
                            got: e.kind(),
                            path: Vec::new(),
                        });
                    };
                    out.insert(*n);
                }
                NumForm::from(out)
            }
            Edn::Str(s) => parse_num(s).map_err(|e| DeError::subgrammar("num", e))?,
            other => {
                return Err(DeError::TypeMismatch {
                    expected: "numeric form (int, :undetermined, vector, or string)",
                    got: other.kind(),
                    path: Vec::new(),
                });
            }
        };
        Ok(Self(v))
    }
}

impl ToEdn for NumDsl {
    fn to_edn(&self) -> Edn<'static> {
        match &self.0 {
            NumForm::Lit(n) => Edn::Int(*n),
            NumForm::Undetermined => Edn::Keyword(EdnKeyword::owned("undetermined".to_string())),
            NumForm::LitSet(xs) => {
                Edn::Vector(xs.iter().map(|n| Edn::Int(*n)).collect::<Vec<_>>().into())
            }
            NumForm::RangeFrom(_)
            | NumForm::RangeTo(_)
            | NumForm::ArithExpr(_)
            | NumForm::PredExpr(_) => Edn::Str(Cow::Owned(self.to_string())),
        }
    }
}

/// Precedence levels (lowest-binding `Or` to highest-binding atom), matching the
/// parser layering. `fmt_*` wraps a child in parens when its level would reparse
/// incorrectly under the parent.
const PREC_OR: u8 = 0;
const PREC_AND: u8 = 1;
const PREC_NOT: u8 = 2;
const PREC_REL: u8 = 3;
const PREC_MEM: u8 = 4;
const PREC_SUM: u8 = 5;
const PREC_PRODUCT: u8 = 6;
const PREC_NEG: u8 = 7;
const PREC_ATOM: u8 = 8;

pub(crate) fn fmt_num(f: &mut fmt::Formatter<'_>, v: &NumForm) -> fmt::Result {
    match v {
        NumForm::Undetermined => f.write_char('*'),
        NumForm::Lit(n) => write!(f, "{}", n),
        NumForm::LitSet(s) => fmt_set(f, s.iter().copied()),
        NumForm::RangeFrom(n) => write!(f, "({n}..)"),
        NumForm::RangeTo(n) => write!(f, "(..{n})"),
        NumForm::ArithExpr(t) => fmt_arith_expr(f, t, PREC_OR),
        NumForm::PredExpr(p) => fmt_pred_expr(f, p, PREC_OR),
    }
}

pub(crate) fn fmt_set<T: fmt::Display>(
    f: &mut fmt::Formatter<'_>,
    values: impl Iterator<Item = T>,
) -> fmt::Result {
    f.write_char('{')?;
    for (i, n) in values.enumerate() {
        if i > 0 {
            f.write_char(',')?;
        }
        write!(f, "{}", n)?;
    }
    f.write_char('}')
}

fn arith_expr_prec(t: &ArithExpr) -> u8 {
    match t {
        ArithExpr::Lit(_) | ArithExpr::Var(_) => PREC_ATOM,
        ArithExpr::Neg(_) => PREC_NEG,
        ArithExpr::Product(_) | ArithExpr::Div(..) | ArithExpr::Rem(..) => PREC_PRODUCT,
        ArithExpr::Sum(_) => PREC_SUM,
    }
}

fn fmt_arith_expr(f: &mut fmt::Formatter<'_>, t: &ArithExpr, parent: u8) -> fmt::Result {
    let prec = arith_expr_prec(t);
    let wrap = prec < parent;
    if wrap {
        f.write_char('(')?;
    }
    match t {
        ArithExpr::Lit(n) => write!(f, "{}", n)?,
        ArithExpr::Var(name) => write!(f, "?{}", name)?,
        ArithExpr::Neg(inner) => {
            f.write_char('-')?;
            fmt_arith_expr(f, inner, PREC_NEG)?;
        }
        ArithExpr::Sum(operands) => {
            for (i, operand) in operands.iter().enumerate() {
                match (i, operand) {
                    (0, _) => fmt_arith_expr(f, operand, PREC_PRODUCT)?,
                    (_, ArithExpr::Neg(inner)) => {
                        f.write_str(" - ")?;
                        fmt_arith_expr(f, inner, PREC_PRODUCT)?;
                    }
                    (_, _) => {
                        f.write_str(" + ")?;
                        fmt_arith_expr(f, operand, PREC_PRODUCT)?;
                    }
                }
            }
        }
        ArithExpr::Product(operands) => {
            for (i, operand) in operands.iter().enumerate() {
                if i > 0 {
                    f.write_str(" * ")?;
                }
                fmt_arith_expr(f, operand, if i == 0 { PREC_PRODUCT } else { PREC_NEG })?;
            }
        }
        ArithExpr::Div(a, b) => {
            fmt_arith_expr(f, a, PREC_PRODUCT)?;
            f.write_str(" / ")?;
            fmt_arith_expr(f, b, PREC_NEG)?;
        }
        ArithExpr::Rem(a, b) => {
            fmt_arith_expr(f, a, PREC_PRODUCT)?;
            f.write_str(" % ")?;
            fmt_arith_expr(f, b, PREC_NEG)?;
        }
    }
    if wrap {
        f.write_char(')')?;
    }
    Ok(())
}

fn pred_expr_prec(p: &PredExpr) -> u8 {
    match p {
        PredExpr::Or(_) => PREC_OR,
        PredExpr::And(_) => PREC_AND,
        PredExpr::Not(_) => PREC_NOT,
        PredExpr::Rel(..) => PREC_REL,
        PredExpr::Mem(..) => PREC_MEM,
    }
}

fn fmt_pred_expr(f: &mut fmt::Formatter<'_>, p: &PredExpr, parent: u8) -> fmt::Result {
    let prec = pred_expr_prec(p);
    let wrap = prec < parent;
    if wrap {
        f.write_char('(')?;
    }
    match p {
        PredExpr::Or(operands) => {
            for (i, operand) in operands.iter().enumerate() {
                if i > 0 {
                    f.write_str(" | ")?;
                }
                fmt_pred_expr(f, operand, PREC_AND)?;
            }
        }
        PredExpr::And(operands) => {
            for (i, operand) in operands.iter().enumerate() {
                if i > 0 {
                    f.write_str(" & ")?;
                }
                fmt_pred_expr(f, operand, PREC_NOT)?;
            }
        }
        PredExpr::Not(inner) => {
            f.write_char('!')?;
            fmt_pred_expr(f, inner, PREC_NOT)?;
        }
        PredExpr::Rel(l, op, r) => {
            fmt_arith_expr(f, l, PREC_OR)?;
            write!(f, " {} ", rel_op_str(*op))?;
            fmt_arith_expr(f, r, PREC_OR)?;
        }
        PredExpr::Mem(e, op, set) => {
            fmt_arith_expr(f, e, PREC_OR)?;
            write!(f, " {} ", mem_op_str(*op))?;
            fmt_set(f, set.iter().copied())?;
        }
    }
    if wrap {
        f.write_char(')')?;
    }
    Ok(())
}

/// Parse a complete numeric DSL string into a `NumForm`.
pub fn parse_num(input: &str) -> Result<NumForm, ParseError> {
    num.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn num(i: &mut &str) -> PResult<NumForm> {
    alt((
        terminated(signed_int, (multispace0, terminator)).map(NumForm::Lit),
        "*".value(NumForm::Undetermined),
        set.map(|v: Vec<i64>| NumForm::lit_set(v)),
        range,
        or_expr.map(parsed_to_num),
    ))
    .parse_next(i)
}

/// A half-open range: `(i..)` → `RangeFrom(i)`, `(..j)` → `RangeTo(j)`. The
/// both-bounded form is intentionally absent (it is a finite set; use `{…}`).
fn range(i: &mut &str) -> PResult<NumForm> {
    delimited(
        ('(', multispace0),
        alt((
            (signed_int, multispace0, "..").map(|(lo, _, _)| NumForm::RangeFrom(lo)),
            ("..", multispace0, signed_int).map(|(_, _, hi)| NumForm::RangeTo(hi)),
        )),
        (multispace0, ')'),
    )
    .parse_next(i)
}

/// A parsed expression carries its sort.
enum Parsed {
    ArithExpr(ArithExpr),
    PredExpr(PredExpr),
}

fn parsed_to_num(parsed: Parsed) -> NumForm {
    match parsed {
        Parsed::ArithExpr(expression) => NumForm::arith_expr(expression),
        Parsed::PredExpr(expression) => NumForm::pred_expr(expression),
    }
}

fn require_arith_expr(parsed: Parsed) -> PResult<ArithExpr> {
    match parsed {
        Parsed::ArithExpr(expression) => Ok(expression),
        Parsed::PredExpr(_) => Err(ErrMode::Backtrack(ParseError::Syntax)),
    }
}

fn require_pred_expr(parsed: Parsed) -> PResult<PredExpr> {
    match parsed {
        Parsed::PredExpr(expression) => Ok(expression),
        Parsed::ArithExpr(_) => Err(ErrMode::Backtrack(ParseError::Syntax)),
    }
}

/// Parse a signed decimal integer matching `[-+]?\d+`, accepting redundant
/// signed-zero spellings and explicit `+`.
pub(crate) fn signed_int(i: &mut &str) -> PResult<i64> {
    let span: &str = (opt(one_of(['-', '+'])), digit1).take().parse_next(i)?;
    span.parse::<i64>()
        .map_err(|_| ErrMode::Backtrack(ParseError::Syntax))
}

fn uint(i: &mut &str) -> PResult<i64> {
    let span: &str = digit1.parse_next(i)?;
    span.parse::<i64>()
        .map_err(|_| ErrMode::Backtrack(ParseError::Syntax))
}

pub(crate) fn variable_name(i: &mut &str) -> PResult<String> {
    (
        one_of(|c: char| c.is_ascii_alphabetic()),
        repeat::<_, _, (), _, _>(0.., one_of(|c: char| c.is_ascii_alphanumeric() || c == '_')),
    )
        .take()
        .map(|s: &str| s.to_string())
        .parse_next(i)
}

pub(crate) fn terminator(i: &mut &str) -> PResult<()> {
    if i.is_empty() || i.starts_with('#') {
        Ok(())
    } else {
        Err(ErrMode::Backtrack(ParseError::Syntax))
    }
}

pub(crate) fn set(i: &mut &str) -> PResult<Vec<i64>> {
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

fn or_expr(i: &mut &str) -> PResult<Parsed> {
    let head = and_expr.parse_next(i)?;
    let rest: Vec<Parsed> = repeat(
        0..,
        preceded(delimited(multispace0, '|', multispace0), and_expr),
    )
    .parse_next(i)?;
    if rest.is_empty() {
        return Ok(head);
    }
    let mut operands = vec![require_pred_expr(head)?];
    for p in rest {
        operands.push(require_pred_expr(p)?);
    }
    Ok(Parsed::PredExpr(PredExpr::Or(operands)))
}

fn and_expr(i: &mut &str) -> PResult<Parsed> {
    let head = not_expr.parse_next(i)?;
    let rest: Vec<Parsed> = repeat(
        0..,
        preceded(delimited(multispace0, '&', multispace0), not_expr),
    )
    .parse_next(i)?;
    if rest.is_empty() {
        return Ok(head);
    }
    let mut operands = vec![require_pred_expr(head)?];
    for p in rest {
        operands.push(require_pred_expr(p)?);
    }
    Ok(Parsed::PredExpr(PredExpr::And(operands)))
}

fn not_expr(i: &mut &str) -> PResult<Parsed> {
    if opt(terminated('!', multispace0)).parse_next(i)?.is_some() {
        let inner = require_pred_expr(not_expr.parse_next(i)?)?;
        Ok(Parsed::PredExpr(PredExpr::Not(Box::new(inner))))
    } else {
        rel_expr.parse_next(i)
    }
}

fn rel_expr(i: &mut &str) -> PResult<Parsed> {
    let lhs = mem_expr.parse_next(i)?;
    let rhs = opt(preceded(
        multispace0,
        (rel_op, preceded(multispace0, mem_expr)),
    ))
    .parse_next(i)?;
    match rhs {
        None => Ok(lhs),
        Some((op, r)) => {
            let l = require_arith_expr(lhs)?;
            let r = require_arith_expr(r)?;
            Ok(Parsed::PredExpr(PredExpr::Rel(l, op, r)))
        }
    }
}

fn mem_expr(i: &mut &str) -> PResult<Parsed> {
    let head = add_expr.parse_next(i)?;
    let membership =
        opt(preceded(multispace0, (mem_op, preceded(multispace0, set)))).parse_next(i)?;
    match membership {
        None => Ok(head),
        Some((op, values)) => {
            let term = require_arith_expr(head)?;
            Ok(Parsed::PredExpr(PredExpr::Mem(
                term,
                op,
                values.into_iter().collect(),
            )))
        }
    }
}

fn add_expr(i: &mut &str) -> PResult<Parsed> {
    let head = mult_expr.parse_next(i)?;
    let tail: Vec<(char, Parsed)> = repeat(
        0..,
        (delimited(multispace0, add_op, multispace0), mult_expr),
    )
    .parse_next(i)?;
    if tail.is_empty() {
        return Ok(head);
    }
    let mut operands = vec![require_arith_expr(head)?];
    for (op, rhs) in tail {
        let rhs = require_arith_expr(rhs)?;
        operands.push(if op == '-' {
            ArithExpr::Neg(Box::new(rhs))
        } else {
            rhs
        });
    }
    Ok(Parsed::ArithExpr(ArithExpr::Sum(operands)))
}

fn add_op(i: &mut &str) -> PResult<char> {
    alt(('+'.value('+'), '-'.value('-'))).parse_next(i)
}

fn mult_expr(i: &mut &str) -> PResult<Parsed> {
    let head = unary_expr.parse_next(i)?;
    let tail: Vec<(char, Parsed)> = repeat(
        0..,
        (delimited(multispace0, mult_op, multispace0), unary_expr),
    )
    .parse_next(i)?;
    if tail.is_empty() {
        return Ok(head);
    }
    let mut acc = require_arith_expr(head)?;
    for (op, rhs) in tail {
        let rhs = require_arith_expr(rhs)?;
        acc = match op {
            '*' => match acc {
                ArithExpr::Product(mut factors) => {
                    factors.push(rhs);
                    ArithExpr::Product(factors)
                }
                other => ArithExpr::Product(vec![other, rhs]),
            },
            '/' => ArithExpr::Div(Box::new(acc), Box::new(rhs)),
            _ => ArithExpr::Rem(Box::new(acc), Box::new(rhs)),
        };
    }
    Ok(Parsed::ArithExpr(acc))
}

fn mult_op(i: &mut &str) -> PResult<char> {
    alt(('*'.value('*'), '/'.value('/'), '%'.value('%'))).parse_next(i)
}

fn unary_expr(i: &mut &str) -> PResult<Parsed> {
    let signs: Vec<char> = repeat(0.., one_of(['-', '+'])).parse_next(i)?;
    let negations = signs.iter().filter(|&&c| c == '-').count();
    let base = base_expr.parse_next(i)?;
    if negations == 0 {
        return Ok(base);
    }
    let mut term = require_arith_expr(base)?;
    for _ in 0..negations {
        term = ArithExpr::Neg(Box::new(term));
    }
    Ok(Parsed::ArithExpr(term))
}

fn base_expr(i: &mut &str) -> PResult<Parsed> {
    alt((
        uint.map(|n| Parsed::ArithExpr(ArithExpr::Lit(n))),
        preceded('?', variable_name).map(|name| Parsed::ArithExpr(ArithExpr::Var(name))),
        delimited('(', delimited(multispace0, or_expr, multispace0), ')'),
    ))
    .parse_next(i)
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::{FromEdn, ToEdn};

    use super::*;
    use crate::ir::operators::{MemOp, RelOp};

    #[rustfmt::skip]
    #[rstest]
    #[case::star("*", NumForm::Undetermined)]
    #[case::num("0", NumForm::Lit(0))]
    #[case::num_neg("-1", NumForm::Lit(-1))]
    #[case::num_pos("+1", NumForm::Lit(1))]
    #[case::num_i64_min("-9223372036854775808", NumForm::Lit(i64::MIN))]
    #[case::set("{0,1,2}", NumForm::lit_set([0, 1, 2]))]
    #[case::set_spaced("{ 0, 1 ,2}", NumForm::lit_set([0, 1, 2]))]
    #[case::range_from("(1..)", NumForm::RangeFrom(1))]
    #[case::range_to("(..3)", NumForm::RangeTo(3))]
    #[case::range_from_neg("(-2..)", NumForm::RangeFrom(-2))]
    #[case::var("?h", NumForm::arith_expr(ArithExpr::Var("h".to_string())))]
    #[case::sum("1 + 2", NumForm::arith_expr(ArithExpr::Sum(vec![ArithExpr::Lit(1), ArithExpr::Lit(2)])))]
    #[case::diff("1 - 2", NumForm::arith_expr(ArithExpr::Sum(vec![ArithExpr::Lit(1), ArithExpr::Neg(Box::new(ArithExpr::Lit(2)))])))]
    #[case::mult("3 * ?h", NumForm::arith_expr(ArithExpr::Product(vec![ArithExpr::Lit(3), ArithExpr::Var("h".to_string())])))]
    #[case::div("10 / 3", NumForm::arith_expr(ArithExpr::Div(Box::new(ArithExpr::Lit(10)), Box::new(ArithExpr::Lit(3)))))]
    #[case::rem("10 % 3", NumForm::arith_expr(ArithExpr::Rem(Box::new(ArithExpr::Lit(10)), Box::new(ArithExpr::Lit(3)))))]
    #[case::neg_var("-?x", NumForm::arith_expr(ArithExpr::Neg(Box::new(ArithExpr::Var("x".to_string())))))]
    #[case::rel_eq("?h == 0", NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Eq, ArithExpr::Lit(0))))]
    #[case::rel_ne("?h != 0", NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Ne, ArithExpr::Lit(0))))]
    #[case::rel_ge("?h >= 1", NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Ge, ArithExpr::Lit(1))))]
    #[case::rel_le("?h <= 1", NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Le, ArithExpr::Lit(1))))]
    #[case::rel_lt("?h < 0", NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Lt, ArithExpr::Lit(0))))]
    #[case::rel_gt("?h > 0", NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Gt, ArithExpr::Lit(0))))]
    #[case::mem_in("?h :: {0,1}", NumForm::pred_expr(PredExpr::Mem(ArithExpr::Var("h".to_string()), MemOp::In, BTreeSet::from([0, 1]))))]
    #[case::mem_notin("?h !: {0,1}", NumForm::pred_expr(PredExpr::Mem(ArithExpr::Var("h".to_string()), MemOp::NotIn, BTreeSet::from([0, 1]))))]
    #[case::not("!?h == 0", NumForm::pred_expr(PredExpr::Not(Box::new(PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Eq, ArithExpr::Lit(0))))))]
    #[case::and("?h == 0 & ?v == 1", NumForm::pred_expr(PredExpr::And(vec![
        PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Eq, ArithExpr::Lit(0)),
        PredExpr::Rel(ArithExpr::Var("v".to_string()), RelOp::Eq, ArithExpr::Lit(1)),
    ])))]
    #[case::or("?h == 0 | ?v == 1", NumForm::pred_expr(PredExpr::Or(vec![
        PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Eq, ArithExpr::Lit(0)),
        PredExpr::Rel(ArithExpr::Var("v".to_string()), RelOp::Eq, ArithExpr::Lit(1)),
    ])))]
    #[case::paren_arith("(0 + 1) * 1", NumForm::arith_expr(ArithExpr::Product(vec![
        ArithExpr::Sum(vec![ArithExpr::Lit(0), ArithExpr::Lit(1)]),
        ArithExpr::Lit(1),
    ])))]
    fn test_parse_num(#[case] input: &str, #[case] expected: NumForm) {
        let result = num.parse(input);
        assert!(result.is_ok(), "{:?} error: {:?}", input, result.clone().unwrap_err());
        assert_eq!(result.unwrap(), expected);
    }

    #[rstest]
    #[case::empty("")]
    #[case::invalid_char("[]")]
    #[case::bare_open_paren("(")]
    #[case::whitespace_variable_name("? h")]
    #[case::spaced_le("?h < = 1")]
    #[case::spaced_ge("?h > = 1")]
    #[case::spaced_eq("?h = = 0")]
    #[case::spaced_mem("?h : : {0,1}")]
    #[case::bare_plus("+")]
    #[case::bare_equal("=")]
    #[case::empty_set("{}")]
    #[case::open_range("(..)")]
    #[case::finite_range("(1..3)")]
    #[case::unclosed_paren_add("(0 + 1")]
    #[case::not_term("!?h")]
    fn test_parse_num_error(#[case] input: &str) {
        let res = num.parse(input);
        assert!(
            res.is_err(),
            "{input:?} should fail, got {:?}",
            res.unwrap()
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(NumForm::Undetermined, "*")]
    #[case::lit_neg(NumForm::Lit(-3), "-3")]
    #[case::set(NumForm::lit_set([0, 1, 2]), "{0,1,2}")]
    #[case::range_from(NumForm::RangeFrom(1), "(1..)")]
    #[case::range_to(NumForm::RangeTo(3), "(..3)")]
    #[case::term_var(NumForm::arith_expr(ArithExpr::Var("h".to_string())), "?h")]
    #[case::term_neg(NumForm::arith_expr(ArithExpr::Neg(Box::new(ArithExpr::Var("x".to_string())))), "-?x")]
    #[case::term_sum(NumForm::arith_expr(ArithExpr::Sum(vec![ArithExpr::Lit(1), ArithExpr::Lit(2)])), "1 + 2")]
    #[case::term_mul(NumForm::arith_expr(ArithExpr::Product(vec![ArithExpr::Lit(3), ArithExpr::Var("h".to_string())])), "3 * ?h")]
    #[case::term_div(NumForm::arith_expr(ArithExpr::Div(Box::new(ArithExpr::Lit(10)), Box::new(ArithExpr::Lit(3)))), "10 / 3")]
    #[case::term_rem(NumForm::arith_expr(ArithExpr::Rem(Box::new(ArithExpr::Lit(10)), Box::new(ArithExpr::Lit(3)))), "10 % 3")]
    #[case::pred_rel(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Eq, ArithExpr::Lit(0))), "?h == 0")]
    #[case::pred_ne(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Ne, ArithExpr::Lit(0))), "?h != 0")]
    #[case::pred_lt(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Lt, ArithExpr::Lit(0))), "?h < 0")]
    #[case::pred_le(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Le, ArithExpr::Lit(1))), "?h <= 1")]
    #[case::pred_gt(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Gt, ArithExpr::Lit(0))), "?h > 0")]
    #[case::pred_ge(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Ge, ArithExpr::Lit(1))), "?h >= 1")]
    #[case::pred_mem(NumForm::pred_expr(PredExpr::Mem(ArithExpr::Var("h".to_string()), MemOp::In, BTreeSet::from([0, 1, 2]))), "?h :: {0,1,2}")]
    #[case::pred_mem_notin(NumForm::pred_expr(PredExpr::Mem(ArithExpr::Var("h".to_string()), MemOp::NotIn, BTreeSet::from([0, 1]))), "?h !: {0,1}")]
    #[case::pred_and_of_or(NumForm::pred_expr(PredExpr::And(vec![
        PredExpr::Or(vec![
            PredExpr::Rel(ArithExpr::Var("a".to_string()), RelOp::Eq, ArithExpr::Lit(0)),
            PredExpr::Rel(ArithExpr::Var("b".to_string()), RelOp::Eq, ArithExpr::Lit(0)),
        ]),
        PredExpr::Rel(ArithExpr::Var("c".to_string()), RelOp::Eq, ArithExpr::Lit(0)),
    ])), "(?a == 0 | ?b == 0) & ?c == 0")]
    fn test_num_dsl_display(#[case] input: NumForm, #[case] expected: &str) {
        assert_eq!(NumDsl::from_ir(&input, &()).to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined("*")]
    #[case::lit("2")]
    #[case::set("{0,1,2}")]
    #[case::var("?h")]
    #[case::add("1 + 2")]
    #[case::sub("1 + -2")]
    #[case::mul_of_add("(0 + 1) * 1")]
    #[case::div("10 / 3")]
    #[case::rem("10 % 3")]
    #[case::rel("?h == 0")]
    #[case::ne("?h != 0")]
    #[case::lt("?h < 0")]
    #[case::le("?h <= 1")]
    #[case::gt("?h > 0")]
    #[case::ge("?h >= 1")]
    #[case::not("!?h == 0")]
    #[case::and("?h == 0 & ?v == 1")]
    #[case::or("?h == 0 | ?v == 1")]
    #[case::and_of_or("(?a == 0 | ?b == 0) & ?c == 0")]
    #[case::mem("?h :: {0,1,2}")]
    #[case::mem_notin("?h !: {0,1}")]
    #[case::range_from("(1..)")]
    #[case::range_to("(..3)")]
    fn test_num_dsl_display_roundtrip(#[case] input: &str) {
        let parsed = num.parse(input).unwrap();
        let rendered = NumDsl::from_ir(&parsed, &()).to_string();
        let reparsed = num.parse(&rendered).unwrap();
        assert_eq!(parsed, reparsed, "input={input:?} rendered={rendered:?}");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(NumForm::Lit(4), Edn::Int(4))]
    #[case::undetermined(NumForm::Undetermined, Edn::Keyword(EdnKeyword::owned("undetermined".into())))]
    #[case::set(NumForm::lit_set([1, 2, 3]), Edn::Vector(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)].into()))]
    #[case::term_var(NumForm::arith_expr(ArithExpr::Var("h".to_string())), Edn::Str(Cow::Borrowed("?h")))]
    #[case::range_from(NumForm::RangeFrom(1), Edn::Str(Cow::Borrowed("(1..)")))]
    #[case::pred_rel(
        NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Eq, ArithExpr::Lit(0))),
        Edn::Str(Cow::Borrowed("?h == 0")),
    )]
    fn test_num_dsl_to_edn(#[case] v: NumForm, #[case] expected: Edn<'static>) {
        assert_eq!(NumDsl::from_ir(&v, &()).to_edn(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::int(Edn::Int(5), NumForm::Lit(5))]
    #[case::keyword(Edn::Keyword(EdnKeyword::owned("undetermined".into())), NumForm::Undetermined)]
    #[case::vector(Edn::Vector(vec![Edn::Int(0), Edn::Int(2)].into()), NumForm::lit_set([0, 2]))]
    #[case::str_int(Edn::Str(Cow::Borrowed("4")), NumForm::Lit(4))]
    #[case::str_set(Edn::Str(Cow::Borrowed("{1,2}")), NumForm::lit_set([1, 2]))]
    #[case::str_range(Edn::Str(Cow::Borrowed("(1..)")), NumForm::RangeFrom(1))]
    #[case::str_pred(
        Edn::Str(Cow::Borrowed("?h == 0")),
        NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Eq, ArithExpr::Lit(0))),
    )]
    fn test_num_dsl_from_edn(#[case] input: Edn<'static>, #[case] expected: NumForm) {
        assert_eq!(NumDsl::from_edn(&input).unwrap().into_ir(&()), expected);
    }

    #[rstest]
    fn test_num_dsl_from_edn_error() {
        let err = NumDsl::from_edn(&Edn::Nil).unwrap_err();
        assert!(matches!(err, DeError::TypeMismatch { .. }));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(NumForm::Lit(3))]
    #[case::undetermined(NumForm::Undetermined)]
    #[case::set(NumForm::lit_set([1, 2, 3]))]
    #[case::range_from(NumForm::RangeFrom(1))]
    #[case::range_to(NumForm::RangeTo(3))]
    #[case::pred(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Ge, ArithExpr::Lit(1))))]
    fn test_num_dsl_edn_roundtrip(#[case] v: NumForm) {
        let edn = NumDsl::from_ir(&v, &()).to_edn();
        let back = NumDsl::from_edn(&edn).unwrap().into_ir(&());
        assert_eq!(back, v);
    }
}
