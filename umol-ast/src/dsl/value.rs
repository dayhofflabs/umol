//! Value DSL: parser, `Display`, EDN boundary. The string and EDN forms of a
//! `ValueAst` are DSL concerns, so they all live on `ValueDsl` here. The AST
//! type itself has no `Display` impl.

use std::borrow::Cow;
use std::fmt::{self, Display, Write};

use umol_edn::{DeError, Edn, EdnKeyword, FromEdn, ToEdn};
use winnow::ascii::{dec_int, digit1, multispace0};
use winnow::combinator::{alt, delimited, opt, preceded, repeat, separated, terminated};
use winnow::error::ErrMode;
use winnow::token::one_of;
use winnow::Parser;

use super::error::{PResult, ParseError};
use crate::ast::traits::{FromAst, IntoAst};
use crate::ast::value::{ArithOp, Expr, RelOp, ValueAst};

/// Surface DSL wrapper around `ValueAst`. EDN form is hybrid: `Lit` → `Int`,
/// `Undetermined` → `:undetermined`, `LitSet` → vector of ints, `Expr` →
/// string via the value subgrammar. `Expr` is string-encoded because EDN has
/// no native representation for the boolean/arithmetic grammar and round-trip
/// fidelity is mandatory.
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
                let mut out = Vec::with_capacity(xs.len());
                for e in xs.iter() {
                    let Edn::Int(n) = e else {
                        return Err(DeError::TypeMismatch {
                            expected: "int (value-set element)",
                            got: e.kind(),
                            path: Vec::new(),
                        });
                    };
                    out.push(*n);
                }
                ValueAst::LitSet(Box::new(out))
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
            ValueAst::Expr(_) => Edn::Str(Cow::Owned(self.to_string())),
        }
    }
}

// region: Format
//
// `Display for ValueDsl` delegates to the helpers below, which know how to
// render any `ValueAst` or `Expr` node in the surface subgrammar. They are
// module-private: DSL callers go through `ValueDsl`; predicates/atom/bond
// strings use `fmt_value` directly (same input type, same output form).

pub(crate) fn fmt_value(f: &mut fmt::Formatter<'_>, v: &ValueAst) -> fmt::Result {
    match v {
        ValueAst::Undetermined => f.write_char('*'),
        ValueAst::Lit(n) => write!(f, "{}", n),
        ValueAst::LitSet(s) => {
            f.write_char('{')?;
            for (i, n) in s.iter().enumerate() {
                if i > 0 {
                    f.write_char(',')?;
                }
                write!(f, "{}", n)?;
            }
            f.write_char('}')
        }
        ValueAst::Expr(e) => fmt_expr(f, e),
    }
}

fn arith_op_str(op: ArithOp) -> &'static str {
    match op {
        ArithOp::Add => "+",
        ArithOp::Sub => "-",
        ArithOp::Mul => "*",
        ArithOp::Div => "/",
        ArithOp::Rem => "%",
    }
}

fn rel_op_str(op: RelOp) -> &'static str {
    match op {
        RelOp::Le => "<=",
        RelOp::Ge => ">=",
        RelOp::Eq => "==",
        RelOp::Lt => "<",
        RelOp::Gt => ">",
    }
}

/// Precedence level for the `Expr` grammar: lowest-binding (`Or`) to
/// highest-binding (atom). Matches the parser's recursive-descent layering.
/// Used to decide where `fmt_expr` wraps a child in parens to reparse to the
/// same tree.
fn expr_prec(e: &Expr) -> u8 {
    match e {
        Expr::Or(_) => 0,
        Expr::And(_) => 1,
        Expr::Not(_) => 2,
        Expr::Rel(..) => 3,
        Expr::Mem(..) => 4,
        Expr::BinOp(_, ArithOp::Add | ArithOp::Sub, _) => 5,
        Expr::BinOp(_, ArithOp::Mul | ArithOp::Div | ArithOp::Rem, _) => 6,
        Expr::Neg(_) => 7,
        Expr::Lit(_) | Expr::Var(_) => 8,
    }
}

fn fmt_paren(f: &mut fmt::Formatter<'_>, e: &Expr, paren: bool) -> fmt::Result {
    if paren {
        f.write_char('(')?;
        fmt_expr(f, e)?;
        f.write_char(')')
    } else {
        fmt_expr(f, e)
    }
}

fn fmt_expr(f: &mut fmt::Formatter<'_>, e: &Expr) -> fmt::Result {
    let parent = expr_prec(e);
    match e {
        Expr::Or(xs) => {
            for (i, x) in xs.iter().enumerate() {
                if i > 0 {
                    f.write_str(" | ")?;
                }
                fmt_paren(f, x, expr_prec(x) < parent)?;
            }
            Ok(())
        }
        Expr::And(xs) => {
            for (i, x) in xs.iter().enumerate() {
                if i > 0 {
                    f.write_str(" & ")?;
                }
                fmt_paren(f, x, expr_prec(x) < parent)?;
            }
            Ok(())
        }
        Expr::Not(x) => {
            f.write_char('!')?;
            fmt_paren(f, x, expr_prec(x) < parent)
        }
        Expr::Rel(l, op, r) => {
            fmt_paren(f, l, expr_prec(l) <= parent)?;
            write!(f, " {} ", rel_op_str(*op))?;
            fmt_paren(f, r, expr_prec(r) <= parent)
        }
        Expr::Mem(inner, set) => {
            fmt_paren(f, inner, expr_prec(inner) < parent)?;
            f.write_str(" :: {")?;
            for (i, n) in set.iter().enumerate() {
                if i > 0 {
                    f.write_char(',')?;
                }
                write!(f, "{}", n)?;
            }
            f.write_char('}')
        }
        Expr::BinOp(l, op, r) => {
            fmt_paren(f, l, expr_prec(l) < parent)?;
            write!(f, " {} ", arith_op_str(*op))?;
            fmt_paren(f, r, expr_prec(r) <= parent)
        }
        Expr::Neg(x) => {
            f.write_char('-')?;
            fmt_paren(f, x, expr_prec(x) < parent)
        }
        Expr::Lit(n) => write!(f, "{}", n),
        Expr::Var(name) => write!(f, "?{}", name),
    }
}

// endregion: Format

// region: Parse

/// Parse a complete value-string into a `ValueAst` (literal, set, or expression).
pub fn parse_value(input: &str) -> Result<ValueAst, ParseError> {
    value.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn value(i: &mut &str) -> PResult<ValueAst> {
    alt((
        terminated(signed_int, (multispace0, terminator)).map(ValueAst::Lit),
        "*".value(ValueAst::Undetermined),
        lit_set.map(ValueAst::lit_set),
        bool_expr.map(ValueAst::expr),
    ))
    .parse_next(i)
}

/// Parse a signed decimal integer matching `[-+]?\d+`. Unlike winnow's
/// `dec_int`, this accepts redundant signed-zero spellings (`-0`, `+0`,
/// `-00`) and explicit `+` sign on positive values, all of which the
/// top-level `value` parser should treat as ground integer literals.
fn signed_int(i: &mut &str) -> PResult<i64> {
    let span: &str = (opt(one_of(['-', '+'])), digit1).take().parse_next(i)?;
    span.parse::<i64>()
        .map_err(|_| ErrMode::Backtrack(ParseError::Syntax))
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

// endregion: Parse

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::star("*", ValueAst::Undetermined)]
    #[case::num("0", ValueAst::Lit(0))]
    #[case::num_neg("-1", ValueAst::Lit(-1))]
    #[case::num_pos("+1", ValueAst::Lit(1))]
    #[case::num_neg_zero("-0", ValueAst::Lit(0))]
    #[case::num_pos_zero("+0", ValueAst::Lit(0))]
    #[case::num_pos_multi("+42", ValueAst::Lit(42))]
    #[case::num_neg_multi("-42", ValueAst::Lit(-42))]
    #[case::num_i64_min("-9223372036854775808", ValueAst::Lit(i64::MIN))]
    #[case::set("{0,1,2}", ValueAst::LitSet(Box::new(vec![0, 1, 2])))]
    #[case::set_spaced("{ 0, 1 ,2}", ValueAst::LitSet(Box::new(vec![0, 1, 2])))]
    #[case::sum("1+2", ValueAst::Expr(Box::new(Expr::BinOp(Box::new(Expr::Lit(1)), ArithOp::Add, Box::new(Expr::Lit(2))))))]
    #[case::sum_spaced("1 + 2", ValueAst::Expr(Box::new(Expr::BinOp(Box::new(Expr::Lit(1)), ArithOp::Add, Box::new(Expr::Lit(2))))))]
    #[case::diff("1-2", ValueAst::Expr(Box::new(Expr::BinOp(Box::new(Expr::Lit(1)), ArithOp::Sub, Box::new(Expr::Lit(2))))))]
    #[case::diff_spaced("1 - 2", ValueAst::Expr(Box::new(Expr::BinOp(Box::new(Expr::Lit(1)), ArithOp::Sub, Box::new(Expr::Lit(2))))))]
    #[case::mult("1*2", ValueAst::Expr(Box::new(Expr::BinOp(Box::new(Expr::Lit(1)), ArithOp::Mul, Box::new(Expr::Lit(2))))))]
    #[case::mult_spaced("1 * 2", ValueAst::Expr(Box::new(Expr::BinOp(Box::new(Expr::Lit(1)), ArithOp::Mul, Box::new(Expr::Lit(2))))))]
    #[case::div("2/2", ValueAst::Expr(Box::new(Expr::BinOp(Box::new(Expr::Lit(2)), ArithOp::Div, Box::new(Expr::Lit(2))))))]
    #[case::div_spaced("2 / 2", ValueAst::Expr(Box::new(Expr::BinOp(Box::new(Expr::Lit(2)), ArithOp::Div, Box::new(Expr::Lit(2))))))]
    #[case::var("?h", ValueAst::Expr(Box::new(Expr::Var("h".to_string()))))]
    #[case::var_2char("?ha", ValueAst::Expr(Box::new(Expr::Var("ha".to_string()))))]
    #[case::var_number("?h1", ValueAst::Expr(Box::new(Expr::Var("h1".to_string()))))]
    #[case::var_underscore("?h_", ValueAst::Expr(Box::new(Expr::Var("h_".to_string()))))]
    #[case::membership("?h + 0 :: {0,1}", ValueAst::Expr(Box::new(Expr::Mem(Box::new(Expr::BinOp(Box::new(Expr::Var("h".to_string())), ArithOp::Add, Box::new(Expr::Lit(0)))), vec![0, 1]))))]
    #[case::double_neg("--0", ValueAst::Expr(Box::new(Expr::Lit(0))))]
    #[case::not_and_precedence("! ?h == 0 & ?v == 1", ValueAst::Expr(Box::new(Expr::And(vec![
        Expr::Not(Box::new(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Eq, Box::new(Expr::Lit(0))))),
        Expr::Rel(Box::new(Expr::Var("v".to_string())), RelOp::Eq, Box::new(Expr::Lit(1))),
    ]))))]
    #[case::paren_arith("(0 + 1) * 1", ValueAst::Expr(Box::new(Expr::BinOp(Box::new(Expr::BinOp(Box::new(Expr::Lit(0)), ArithOp::Add, Box::new(Expr::Lit(1)))), ArithOp::Mul, Box::new(Expr::Lit(1))))))]
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

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(ValueAst::Undetermined, "*")]
    #[case::lit_zero(ValueAst::Lit(0), "0")]
    #[case::lit_neg(ValueAst::Lit(-3), "-3")]
    #[case::set(ValueAst::LitSet(Box::new(vec![0, 1, 2])), "{0,1,2}")]
    #[case::set_single(ValueAst::LitSet(Box::new(vec![5])), "{5}")]
    #[case::expr_lit(ValueAst::Expr(Box::new(Expr::Lit(7))), "7")]
    #[case::expr_var(ValueAst::Expr(Box::new(Expr::Var("h".into()))), "?h")]
    #[case::expr_neg(ValueAst::Expr(Box::new(Expr::Neg(Box::new(Expr::Var("x".into()))))), "-?x")]
    #[case::expr_add(ValueAst::Expr(Box::new(Expr::BinOp(Box::new(Expr::Lit(1)), ArithOp::Add, Box::new(Expr::Lit(2))))), "1 + 2")]
    #[case::expr_mul(ValueAst::Expr(Box::new(Expr::BinOp(Box::new(Expr::Lit(3)), ArithOp::Mul, Box::new(Expr::Var("h".into()))))), "3 * ?h")]
    #[case::expr_rel(ValueAst::Expr(Box::new(Expr::Rel(Box::new(Expr::Var("h".into())), RelOp::Eq, Box::new(Expr::Lit(0))))), "?h == 0")]
    #[case::expr_mem(ValueAst::Expr(Box::new(Expr::Mem(Box::new(Expr::Var("h".into())), vec![0, 1, 2]))), "?h :: {0,1,2}")]
    #[case::expr_not(ValueAst::Expr(Box::new(Expr::Not(Box::new(Expr::Rel(Box::new(Expr::Var("h".into())), RelOp::Eq, Box::new(Expr::Lit(0))))))), "!?h == 0")]
    #[case::expr_and(ValueAst::Expr(Box::new(Expr::And(vec![
        Expr::Rel(Box::new(Expr::Var("h".into())), RelOp::Eq, Box::new(Expr::Lit(0))),
        Expr::Rel(Box::new(Expr::Var("v".into())), RelOp::Eq, Box::new(Expr::Lit(1))),
    ]))), "?h == 0 & ?v == 1")]
    #[case::expr_or(ValueAst::Expr(Box::new(Expr::Or(vec![
        Expr::Rel(Box::new(Expr::Var("h".into())), RelOp::Eq, Box::new(Expr::Lit(0))),
        Expr::Rel(Box::new(Expr::Var("v".into())), RelOp::Eq, Box::new(Expr::Lit(1))),
    ]))), "?h == 0 | ?v == 1")]
    // And of Or: Or children need parens because And binds tighter.
    #[case::and_of_or(ValueAst::Expr(Box::new(Expr::And(vec![
        Expr::Or(vec![Expr::Var("a".into()), Expr::Var("b".into())]),
        Expr::Var("c".into()),
    ]))), "(?a | ?b) & ?c")]
    // Subtraction is non-right-associative at same precedence.
    #[case::sub_right_nests(ValueAst::Expr(Box::new(Expr::BinOp(
        Box::new(Expr::Lit(1)),
        ArithOp::Sub,
        Box::new(Expr::BinOp(Box::new(Expr::Lit(2)), ArithOp::Sub, Box::new(Expr::Lit(3)))),
    ))), "1 - (2 - 3)")]
    // (0 + 1) * 1 — left child of Mul is Add (lower prec) → parens needed.
    #[case::mul_of_add(ValueAst::Expr(Box::new(Expr::BinOp(
        Box::new(Expr::BinOp(Box::new(Expr::Lit(0)), ArithOp::Add, Box::new(Expr::Lit(1)))),
        ArithOp::Mul,
        Box::new(Expr::Lit(1)),
    ))), "(0 + 1) * 1")]
    fn test_value_display(#[case] input: ValueAst, #[case] expected: &str) {
        assert_eq!(ValueDsl::from_ast(&input, &()).to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined("*")]
    #[case::lit("2")]
    #[case::lit_neg("-3")]
    #[case::set("{0,1,2}")]
    #[case::var("?h")]
    #[case::add("1 + 2")]
    #[case::mul_of_add("(0 + 1) * 1")]
    #[case::sub_right("1 - (2 - 3)")]
    #[case::rel("?h == 0")]
    #[case::not("!?h == 0")]
    #[case::and("?h == 0 & ?v == 1")]
    #[case::or("?h == 0 | ?v == 1")]
    #[case::and_of_or("(?a | ?b) & ?c")]
    #[case::mem("?h :: {0,1,2}")]
    #[case::chained_and_or("?a & ?b | ?c & ?d")]
    fn test_value_display_roundtrip(#[case] input: &str) {
        let parsed = value.parse(input).unwrap();
        let rendered = ValueDsl::from_ast(&parsed, &()).to_string();
        let reparsed = value.parse(&rendered).unwrap();
        assert_eq!(parsed, reparsed, "input={input:?} rendered={rendered:?}");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(ValueAst::Lit(4), Edn::Int(4))]
    #[case::lit_neg(ValueAst::Lit(-2), Edn::Int(-2))]
    #[case::undetermined(ValueAst::Undetermined, Edn::Keyword(EdnKeyword::owned("undetermined".into())))]
    #[case::set(ValueAst::LitSet(Box::new(vec![1, 2, 3])), Edn::Vector(vec![Edn::Int(1), Edn::Int(2), Edn::Int(3)].into()))]
    #[case::expr_var(ValueAst::Expr(Box::new(Expr::Var("h".into()))), Edn::Str(Cow::Borrowed("?h")))]
    #[case::expr_rel(
        ValueAst::Expr(Box::new(Expr::Rel(Box::new(Expr::Var("h".into())), RelOp::Eq, Box::new(Expr::Lit(0))))),
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
    #[case::vector(Edn::Vector(vec![Edn::Int(0), Edn::Int(2)].into()), ValueAst::LitSet(Box::new(vec![0, 2])))]
    #[case::str_int(Edn::Str(Cow::Borrowed("4")), ValueAst::Lit(4))]
    #[case::str_undetermined(Edn::Str(Cow::Borrowed("*")), ValueAst::Undetermined)]
    #[case::str_set(Edn::Str(Cow::Borrowed("{1,2}")), ValueAst::LitSet(Box::new(vec![1, 2])))]
    #[case::str_expr(
        Edn::Str(Cow::Borrowed("?h == 0")),
        ValueAst::Expr(Box::new(Expr::Rel(Box::new(Expr::Var("h".into())), RelOp::Eq, Box::new(Expr::Lit(0))))),
    )]
    fn test_value_dsl_from_edn(#[case] input: Edn<'static>, #[case] expected: ValueAst) {
        use umol_edn::FromEdn;
        assert_eq!(ValueDsl::from_edn(&input).unwrap().into_ast(&()), expected);
    }

    #[rstest]
    fn test_value_dsl_from_edn_rejects_wrong_kind() {
        use umol_edn::FromEdn;
        let err = ValueDsl::from_edn(&Edn::Nil).unwrap_err();
        assert!(matches!(err, DeError::TypeMismatch { .. }));
    }

    #[rstest]
    fn test_value_dsl_from_edn_rejects_non_int_in_vector() {
        use umol_edn::FromEdn;
        let err = ValueDsl::from_edn(&Edn::Vector(vec![Edn::Int(1), Edn::Nil].into())).unwrap_err();
        assert!(matches!(err, DeError::TypeMismatch { .. }));
    }

    #[rstest]
    fn test_value_dsl_from_edn_rejects_invalid_string() {
        use umol_edn::FromEdn;
        let err = ValueDsl::from_edn(&Edn::Str(Cow::Borrowed("???"))).unwrap_err();
        assert!(matches!(
            err,
            DeError::Subgrammar {
                grammar: "value",
                ..
            }
        ));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(ValueAst::Lit(3))]
    #[case::undetermined(ValueAst::Undetermined)]
    #[case::set(ValueAst::LitSet(Box::new(vec![1, 2, 3])))]
    #[case::expr(ValueAst::Expr(Box::new(Expr::Rel(Box::new(Expr::Var("h".into())), RelOp::Ge, Box::new(Expr::Lit(1))))))]
    fn test_value_dsl_edn_roundtrip(#[case] v: ValueAst) {
        use umol_edn::{FromEdn, ToEdn};
        let edn = ValueDsl::from_ast(&v, &()).to_edn();
        let back = ValueDsl::from_edn(&edn).unwrap().into_ast(&());
        assert_eq!(back, v);
    }
}
