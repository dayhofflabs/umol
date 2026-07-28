//! Value DSL: parser, `Display`, EDN boundary. The canonicalization is run lazily.

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
use crate::ast::traits::{FromAst, IntoAst};
use crate::ast::value::{ValueAst, ValuePredicate, ValueTerm};

/// Surface DSL wrapper around `ValueAst`. EDN form is hybrid: `Lit` → `Int`,
/// `Undetermined` → `:undetermined`, `LitSet` → vector of ints, `Term`/
/// `Predicate` → string via the value subgrammar (EDN has no native form for
/// the arithmetic/boolean grammar and round-trip fidelity is mandatory).
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct ValueDsl(pub ValueAst);

impl FromAst<ValueAst> for ValueDsl {
    type Ctx = ();

    fn from_ast(ast: &ValueAst, _ctx: &Self::Ctx) -> Self {
        Self(ast.clone())
    }
}

impl IntoAst<ValueAst> for ValueDsl {
    type Ctx = ();

    fn into_ast(self, _ctx: &Self::Ctx) -> ValueAst {
        self.0
    }
}

impl Display for ValueDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_value(f, &self.0)
    }
}

impl<'de> FromEdn<'de> for ValueDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let v = match edn {
            Edn::Int(n) => ValueAst::Lit(*n),
            Edn::Keyword(k) if k.name() == "undetermined" => ValueAst::Undetermined,
            Edn::Vector(xs) => {
                let mut out = BTreeSet::new();
                for e in xs.iter() {
                    let Edn::Int(n) = e else {
                        return Err(DeError::TypeMismatch {
                            expected: "int (value-set element)",
                            got: e.kind(),
                            path: Vec::new(),
                        });
                    };
                    out.insert(*n);
                }
                ValueAst::from(out)
            }
            Edn::Str(s) => parse_value(s).map_err(|e| DeError::subgrammar("value", e))?,
            other => {
                return Err(DeError::TypeMismatch {
                    expected: "value (int, :undetermined, vector, or string)",
                    got: other.kind(),
                    path: Vec::new(),
                });
            }
        };
        Ok(Self(v))
    }
}

impl ToEdn for ValueDsl {
    fn to_edn(&self) -> Edn<'static> {
        match &self.0 {
            ValueAst::Lit(n) => Edn::Int(*n),
            ValueAst::Undetermined => Edn::Keyword(EdnKeyword::owned("undetermined".to_string())),
            ValueAst::LitSet(xs) => {
                Edn::Vector(xs.iter().map(|n| Edn::Int(*n)).collect::<Vec<_>>().into())
            }
            ValueAst::RangeFrom(_)
            | ValueAst::RangeTo(_)
            | ValueAst::Term(_)
            | ValueAst::Predicate(_) => Edn::Str(Cow::Owned(self.to_string())),
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

pub(crate) fn fmt_value(f: &mut fmt::Formatter<'_>, v: &ValueAst) -> fmt::Result {
    match v {
        ValueAst::Undetermined => f.write_char('*'),
        ValueAst::Lit(n) => write!(f, "{}", n),
        ValueAst::LitSet(s) => fmt_set(f, s.iter().copied()),
        ValueAst::RangeFrom(n) => write!(f, "({n}..)"),
        ValueAst::RangeTo(n) => write!(f, "(..{n})"),
        ValueAst::Term(t) => fmt_term(f, t, PREC_OR),
        ValueAst::Predicate(p) => fmt_predicate(f, p, PREC_OR),
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

fn term_prec(t: &ValueTerm) -> u8 {
    match t {
        ValueTerm::Lit(_) | ValueTerm::Var(_) => PREC_ATOM,
        ValueTerm::Neg(_) => PREC_NEG,
        ValueTerm::Product(_) | ValueTerm::Div(..) | ValueTerm::Rem(..) => PREC_PRODUCT,
        ValueTerm::Sum(_) => PREC_SUM,
    }
}

fn fmt_term(f: &mut fmt::Formatter<'_>, t: &ValueTerm, parent: u8) -> fmt::Result {
    let prec = term_prec(t);
    let wrap = prec < parent;
    if wrap {
        f.write_char('(')?;
    }
    match t {
        ValueTerm::Lit(n) => write!(f, "{}", n)?,
        ValueTerm::Var(name) => write!(f, "?{}", name)?,
        ValueTerm::Neg(inner) => {
            f.write_char('-')?;
            fmt_term(f, inner, PREC_NEG)?;
        }
        ValueTerm::Sum(operands) => {
            for (i, operand) in operands.iter().enumerate() {
                match (i, operand) {
                    (0, _) => fmt_term(f, operand, PREC_PRODUCT)?,
                    (_, ValueTerm::Neg(inner)) => {
                        f.write_str(" - ")?;
                        fmt_term(f, inner, PREC_PRODUCT)?;
                    }
                    (_, _) => {
                        f.write_str(" + ")?;
                        fmt_term(f, operand, PREC_PRODUCT)?;
                    }
                }
            }
        }
        ValueTerm::Product(operands) => {
            for (i, operand) in operands.iter().enumerate() {
                if i > 0 {
                    f.write_str(" * ")?;
                }
                fmt_term(f, operand, if i == 0 { PREC_PRODUCT } else { PREC_NEG })?;
            }
        }
        ValueTerm::Div(a, b) => {
            fmt_term(f, a, PREC_PRODUCT)?;
            f.write_str(" / ")?;
            fmt_term(f, b, PREC_NEG)?;
        }
        ValueTerm::Rem(a, b) => {
            fmt_term(f, a, PREC_PRODUCT)?;
            f.write_str(" % ")?;
            fmt_term(f, b, PREC_NEG)?;
        }
    }
    if wrap {
        f.write_char(')')?;
    }
    Ok(())
}

fn pred_prec(p: &ValuePredicate) -> u8 {
    match p {
        ValuePredicate::Or(_) => PREC_OR,
        ValuePredicate::And(_) => PREC_AND,
        ValuePredicate::Not(_) => PREC_NOT,
        ValuePredicate::Rel(..) => PREC_REL,
        ValuePredicate::Mem(..) => PREC_MEM,
    }
}

fn fmt_predicate(f: &mut fmt::Formatter<'_>, p: &ValuePredicate, parent: u8) -> fmt::Result {
    let prec = pred_prec(p);
    let wrap = prec < parent;
    if wrap {
        f.write_char('(')?;
    }
    match p {
        ValuePredicate::Or(operands) => {
            for (i, operand) in operands.iter().enumerate() {
                if i > 0 {
                    f.write_str(" | ")?;
                }
                fmt_predicate(f, operand, PREC_AND)?;
            }
        }
        ValuePredicate::And(operands) => {
            for (i, operand) in operands.iter().enumerate() {
                if i > 0 {
                    f.write_str(" & ")?;
                }
                fmt_predicate(f, operand, PREC_NOT)?;
            }
        }
        ValuePredicate::Not(inner) => {
            f.write_char('!')?;
            fmt_predicate(f, inner, PREC_NOT)?;
        }
        ValuePredicate::Rel(l, op, r) => {
            fmt_term(f, l, PREC_OR)?;
            write!(f, " {} ", rel_op_str(*op))?;
            fmt_term(f, r, PREC_OR)?;
        }
        ValuePredicate::Mem(e, op, set) => {
            fmt_term(f, e, PREC_OR)?;
            write!(f, " {} ", mem_op_str(*op))?;
            fmt_set(f, set.iter().copied())?;
        }
    }
    if wrap {
        f.write_char(')')?;
    }
    Ok(())
}

/// Parse a complete value-string into a `ValueAst`.
pub fn parse_value(input: &str) -> Result<ValueAst, ParseError> {
    value.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn value(i: &mut &str) -> PResult<ValueAst> {
    alt((
        terminated(signed_int, (multispace0, terminator)).map(ValueAst::Lit),
        "*".value(ValueAst::Undetermined),
        set.map(|v: Vec<i64>| ValueAst::lit_set(v)),
        range,
        or_expr.map(parsed_to_value),
    ))
    .parse_next(i)
}

/// A half-open range: `(i..)` → `RangeFrom(i)`, `(..j)` → `RangeTo(j)`. The
/// both-bounded form is intentionally absent (it is a finite set; use `{…}`).
fn range(i: &mut &str) -> PResult<ValueAst> {
    delimited(
        ('(', multispace0),
        alt((
            (signed_int, multispace0, "..").map(|(lo, _, _)| ValueAst::RangeFrom(lo)),
            ("..", multispace0, signed_int).map(|(_, _, hi)| ValueAst::RangeTo(hi)),
        )),
        (multispace0, ')'),
    )
    .parse_next(i)
}

/// A parsed expression carries its sort: the arithmetic operators build a
/// `Term`, the relational/membership/boolean operators build a `Predicate`.
enum Parsed {
    Term(ValueTerm),
    Predicate(ValuePredicate),
}

fn parsed_to_value(parsed: Parsed) -> ValueAst {
    match parsed {
        Parsed::Term(t) => ValueAst::term(t),
        Parsed::Predicate(p) => ValueAst::predicate(p),
    }
}

fn require_term(parsed: Parsed) -> PResult<ValueTerm> {
    match parsed {
        Parsed::Term(t) => Ok(t),
        Parsed::Predicate(_) => Err(ErrMode::Backtrack(ParseError::Syntax)),
    }
}

fn require_predicate(parsed: Parsed) -> PResult<ValuePredicate> {
    match parsed {
        Parsed::Predicate(p) => Ok(p),
        Parsed::Term(_) => Err(ErrMode::Backtrack(ParseError::Syntax)),
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
    let mut operands = vec![require_predicate(head)?];
    for p in rest {
        operands.push(require_predicate(p)?);
    }
    Ok(Parsed::Predicate(ValuePredicate::Or(operands)))
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
    let mut operands = vec![require_predicate(head)?];
    for p in rest {
        operands.push(require_predicate(p)?);
    }
    Ok(Parsed::Predicate(ValuePredicate::And(operands)))
}

fn not_expr(i: &mut &str) -> PResult<Parsed> {
    if opt(terminated('!', multispace0)).parse_next(i)?.is_some() {
        let inner = require_predicate(not_expr.parse_next(i)?)?;
        Ok(Parsed::Predicate(ValuePredicate::Not(Box::new(inner))))
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
            let l = require_term(lhs)?;
            let r = require_term(r)?;
            Ok(Parsed::Predicate(ValuePredicate::Rel(l, op, r)))
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
            let term = require_term(head)?;
            Ok(Parsed::Predicate(ValuePredicate::Mem(
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
    let mut operands = vec![require_term(head)?];
    for (op, rhs) in tail {
        let rhs = require_term(rhs)?;
        operands.push(if op == '-' {
            ValueTerm::Neg(Box::new(rhs))
        } else {
            rhs
        });
    }
    Ok(Parsed::Term(ValueTerm::Sum(operands)))
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
    let mut acc = require_term(head)?;
    for (op, rhs) in tail {
        let rhs = require_term(rhs)?;
        acc = match op {
            '*' => match acc {
                ValueTerm::Product(mut factors) => {
                    factors.push(rhs);
                    ValueTerm::Product(factors)
                }
                other => ValueTerm::Product(vec![other, rhs]),
            },
            '/' => ValueTerm::Div(Box::new(acc), Box::new(rhs)),
            _ => ValueTerm::Rem(Box::new(acc), Box::new(rhs)),
        };
    }
    Ok(Parsed::Term(acc))
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
    let mut term = require_term(base)?;
    for _ in 0..negations {
        term = ValueTerm::Neg(Box::new(term));
    }
    Ok(Parsed::Term(term))
}

fn base_expr(i: &mut &str) -> PResult<Parsed> {
    alt((
        uint.map(|n| Parsed::Term(ValueTerm::Lit(n))),
        preceded('?', variable_name).map(|name| Parsed::Term(ValueTerm::Var(name))),
        delimited('(', delimited(multispace0, or_expr, multispace0), ')'),
    ))
    .parse_next(i)
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::operators::{MemOp, RelOp};

    #[rustfmt::skip]
    #[rstest]
    #[case::star("*", ValueAst::Undetermined)]
    #[case::num("0", ValueAst::Lit(0))]
    #[case::num_neg("-1", ValueAst::Lit(-1))]
    #[case::num_pos("+1", ValueAst::Lit(1))]
    #[case::num_i64_min("-9223372036854775808", ValueAst::Lit(i64::MIN))]
    #[case::set("{0,1,2}", ValueAst::lit_set([0, 1, 2]))]
    #[case::set_spaced("{ 0, 1 ,2}", ValueAst::lit_set([0, 1, 2]))]
    #[case::range_from("(1..)", ValueAst::RangeFrom(1))]
    #[case::range_to("(..3)", ValueAst::RangeTo(3))]
    #[case::range_from_neg("(-2..)", ValueAst::RangeFrom(-2))]
    #[case::var("?h", ValueAst::term(ValueTerm::Var("h".to_string())))]
    #[case::sum("1 + 2", ValueAst::term(ValueTerm::Sum(vec![ValueTerm::Lit(1), ValueTerm::Lit(2)])))]
    #[case::diff("1 - 2", ValueAst::term(ValueTerm::Sum(vec![ValueTerm::Lit(1), ValueTerm::Neg(Box::new(ValueTerm::Lit(2)))])))]
    #[case::mult("3 * ?h", ValueAst::term(ValueTerm::Product(vec![ValueTerm::Lit(3), ValueTerm::Var("h".to_string())])))]
    #[case::div("10 / 3", ValueAst::term(ValueTerm::Div(Box::new(ValueTerm::Lit(10)), Box::new(ValueTerm::Lit(3)))))]
    #[case::rem("10 % 3", ValueAst::term(ValueTerm::Rem(Box::new(ValueTerm::Lit(10)), Box::new(ValueTerm::Lit(3)))))]
    #[case::neg_var("-?x", ValueAst::term(ValueTerm::Neg(Box::new(ValueTerm::Var("x".to_string())))))]
    #[case::rel_eq("?h == 0", ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Eq, ValueTerm::Lit(0))))]
    #[case::rel_ne("?h != 0", ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Ne, ValueTerm::Lit(0))))]
    #[case::rel_ge("?h >= 1", ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Ge, ValueTerm::Lit(1))))]
    #[case::rel_le("?h <= 1", ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Le, ValueTerm::Lit(1))))]
    #[case::rel_lt("?h < 0", ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Lt, ValueTerm::Lit(0))))]
    #[case::rel_gt("?h > 0", ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Gt, ValueTerm::Lit(0))))]
    #[case::mem_in("?h :: {0,1}", ValueAst::predicate(ValuePredicate::Mem(ValueTerm::Var("h".to_string()), MemOp::In, BTreeSet::from([0, 1]))))]
    #[case::mem_notin("?h !: {0,1}", ValueAst::predicate(ValuePredicate::Mem(ValueTerm::Var("h".to_string()), MemOp::NotIn, BTreeSet::from([0, 1]))))]
    #[case::not("!?h == 0", ValueAst::predicate(ValuePredicate::Not(Box::new(ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Eq, ValueTerm::Lit(0))))))]
    #[case::and("?h == 0 & ?v == 1", ValueAst::predicate(ValuePredicate::And(vec![
        ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Eq, ValueTerm::Lit(0)),
        ValuePredicate::Rel(ValueTerm::Var("v".to_string()), RelOp::Eq, ValueTerm::Lit(1)),
    ])))]
    #[case::or("?h == 0 | ?v == 1", ValueAst::predicate(ValuePredicate::Or(vec![
        ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Eq, ValueTerm::Lit(0)),
        ValuePredicate::Rel(ValueTerm::Var("v".to_string()), RelOp::Eq, ValueTerm::Lit(1)),
    ])))]
    #[case::paren_arith("(0 + 1) * 1", ValueAst::term(ValueTerm::Product(vec![
        ValueTerm::Sum(vec![ValueTerm::Lit(0), ValueTerm::Lit(1)]),
        ValueTerm::Lit(1),
    ])))]
    fn test_value(#[case] input: &str, #[case] expected: ValueAst) {
        let result = value.parse(input);
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
    fn test_value_error(#[case] input: &str) {
        let res = value.parse(input);
        assert!(
            res.is_err(),
            "{input:?} should fail, got {:?}",
            res.unwrap()
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(ValueAst::Undetermined, "*")]
    #[case::lit_neg(ValueAst::Lit(-3), "-3")]
    #[case::set(ValueAst::lit_set([0, 1, 2]), "{0,1,2}")]
    #[case::range_from(ValueAst::RangeFrom(1), "(1..)")]
    #[case::range_to(ValueAst::RangeTo(3), "(..3)")]
    #[case::term_var(ValueAst::term(ValueTerm::Var("h".to_string())), "?h")]
    #[case::term_neg(ValueAst::term(ValueTerm::Neg(Box::new(ValueTerm::Var("x".to_string())))), "-?x")]
    #[case::term_sum(ValueAst::term(ValueTerm::Sum(vec![ValueTerm::Lit(1), ValueTerm::Lit(2)])), "1 + 2")]
    #[case::term_mul(ValueAst::term(ValueTerm::Product(vec![ValueTerm::Lit(3), ValueTerm::Var("h".to_string())])), "3 * ?h")]
    #[case::term_div(ValueAst::term(ValueTerm::Div(Box::new(ValueTerm::Lit(10)), Box::new(ValueTerm::Lit(3)))), "10 / 3")]
    #[case::term_rem(ValueAst::term(ValueTerm::Rem(Box::new(ValueTerm::Lit(10)), Box::new(ValueTerm::Lit(3)))), "10 % 3")]
    #[case::pred_rel(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Eq, ValueTerm::Lit(0))), "?h == 0")]
    #[case::pred_ne(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Ne, ValueTerm::Lit(0))), "?h != 0")]
    #[case::pred_lt(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Lt, ValueTerm::Lit(0))), "?h < 0")]
    #[case::pred_le(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Le, ValueTerm::Lit(1))), "?h <= 1")]
    #[case::pred_gt(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Gt, ValueTerm::Lit(0))), "?h > 0")]
    #[case::pred_ge(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Ge, ValueTerm::Lit(1))), "?h >= 1")]
    #[case::pred_mem(ValueAst::predicate(ValuePredicate::Mem(ValueTerm::Var("h".to_string()), MemOp::In, BTreeSet::from([0, 1, 2]))), "?h :: {0,1,2}")]
    #[case::pred_mem_notin(ValueAst::predicate(ValuePredicate::Mem(ValueTerm::Var("h".to_string()), MemOp::NotIn, BTreeSet::from([0, 1]))), "?h !: {0,1}")]
    #[case::pred_and_of_or(ValueAst::predicate(ValuePredicate::And(vec![
        ValuePredicate::Or(vec![
            ValuePredicate::Rel(ValueTerm::Var("a".to_string()), RelOp::Eq, ValueTerm::Lit(0)),
            ValuePredicate::Rel(ValueTerm::Var("b".to_string()), RelOp::Eq, ValueTerm::Lit(0)),
        ]),
        ValuePredicate::Rel(ValueTerm::Var("c".to_string()), RelOp::Eq, ValueTerm::Lit(0)),
    ])), "(?a == 0 | ?b == 0) & ?c == 0")]
    fn test_value_display(#[case] input: ValueAst, #[case] expected: &str) {
        assert_eq!(ValueDsl::from_ast(&input, &()).to_string(), expected);
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
    fn test_value_display_roundtrip(#[case] input: &str) {
        let parsed = value.parse(input).unwrap();
        let rendered = ValueDsl::from_ast(&parsed, &()).to_string();
        let reparsed = value.parse(&rendered).unwrap();
        assert_eq!(parsed, reparsed, "input={input:?} rendered={rendered:?}");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(ValueAst::Lit(4), Edn::Int(4))]
    #[case::undetermined(ValueAst::Undetermined, Edn::Keyword(EdnKeyword::owned("undetermined".into())))]
    #[case::set(ValueAst::lit_set([1, 2, 3]), Edn::Vector(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)].into()))]
    #[case::term_var(ValueAst::term(ValueTerm::Var("h".to_string())), Edn::Str(Cow::Borrowed("?h")))]
    #[case::range_from(ValueAst::RangeFrom(1), Edn::Str(Cow::Borrowed("(1..)")))]
    #[case::pred_rel(
        ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Eq, ValueTerm::Lit(0))),
        Edn::Str(Cow::Borrowed("?h == 0")),
    )]
    fn test_value_dsl_to_edn(#[case] v: ValueAst, #[case] expected: Edn<'static>) {
        use umol_edn::ToEdn;
        assert_eq!(ValueDsl::from_ast(&v, &()).to_edn(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::int(Edn::Int(5), ValueAst::Lit(5))]
    #[case::keyword(Edn::Keyword(EdnKeyword::owned("undetermined".into())), ValueAst::Undetermined)]
    #[case::vector(Edn::Vector(vec![Edn::Int(0), Edn::Int(2)].into()), ValueAst::lit_set([0, 2]))]
    #[case::str_int(Edn::Str(Cow::Borrowed("4")), ValueAst::Lit(4))]
    #[case::str_set(Edn::Str(Cow::Borrowed("{1,2}")), ValueAst::lit_set([1, 2]))]
    #[case::str_range(Edn::Str(Cow::Borrowed("(1..)")), ValueAst::RangeFrom(1))]
    #[case::str_pred(
        Edn::Str(Cow::Borrowed("?h == 0")),
        ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Eq, ValueTerm::Lit(0))),
    )]
    fn test_value_dsl_from_edn(#[case] input: Edn<'static>, #[case] expected: ValueAst) {
        use umol_edn::FromEdn;
        assert_eq!(ValueDsl::from_edn(&input).unwrap().into_ast(&()), expected);
    }

    #[rstest]
    fn test_value_dsl_from_edn_error() {
        use umol_edn::FromEdn;
        let err = ValueDsl::from_edn(&Edn::Nil).unwrap_err();
        assert!(matches!(err, DeError::TypeMismatch { .. }));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(ValueAst::Lit(3))]
    #[case::undetermined(ValueAst::Undetermined)]
    #[case::set(ValueAst::lit_set([1, 2, 3]))]
    #[case::range_from(ValueAst::RangeFrom(1))]
    #[case::range_to(ValueAst::RangeTo(3))]
    #[case::pred(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Ge, ValueTerm::Lit(1))))]
    fn test_value_dsl_edn_roundtrip(#[case] v: ValueAst) {
        use umol_edn::{FromEdn, ToEdn};
        let edn = ValueDsl::from_ast(&v, &()).to_edn();
        let back = ValueDsl::from_edn(&edn).unwrap().into_ast(&());
        assert_eq!(back, v);
    }
}
