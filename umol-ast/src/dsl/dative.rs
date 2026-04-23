//! Dative-bond-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, FromEdn, ToEdn};
use winnow::ascii::multispace0;
use winnow::combinator::{repeat, terminated};
use winnow::error::ErrMode;
use winnow::token::take;
use winnow::Parser;

use super::error::{PResult, ParseError};
use super::predicates::{fmt_ring_count, fmt_value, optional_value, ring_count};
use crate::ast::config::DativeBondAstConfig;
use crate::ast::constraint::DativeBondConstraint;
use crate::ast::dative::DativeBondAst;
use crate::ast::traits::{FromAst, ToAst};
use crate::ast::value::ValueAst;

/// Surface DSL wrapper around `DativeBondAst`. No leading token; the string
/// form is a sequence of `#…` predicates. Inline-capable constraints from
/// `DativeBondConstraint` are `RingCount` (`#R`) and `RingSize` (`#r`); the
/// remaining variants reference other entities and stay in the molecule
/// constraints container.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DativeBondDsl(pub DativeBondAst);

impl FromStr for DativeBondDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_dative(s)
    }
}

impl Display for DativeBondDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for c in self.0.constraints.iter() {
            fmt_constraint(f, c)?;
        }
        Ok(())
    }
}

impl<'de> FromEdn<'de> for DativeBondDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s.parse().map_err(|e| DeError::subgrammar("dative", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for DativeBondDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

impl FromAst<DativeBondAst> for DativeBondDsl {
    type Error = ParseError;

    fn from_ast(ast: &DativeBondAst, _cfg: &DativeBondAstConfig) -> Result<Self, ParseError> {
        Ok(DativeBondDsl(ast.clone()))
    }
}

impl ToAst<DativeBondAst> for DativeBondDsl {
    type Error = ParseError;

    fn to_ast(&self, _cfg: &DativeBondAstConfig) -> Result<DativeBondAst, ParseError> {
        Ok(self.0.clone())
    }
}

// -- Parse --------------------

pub fn parse_dative(input: &str) -> Result<DativeBondDsl, ParseError> {
    dative.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn dative(i: &mut &str) -> PResult<DativeBondDsl> {
    multispace0.parse_next(i)?;
    let preds: Vec<DativePredicate> =
        repeat(0.., terminated(dative_predicate, multispace0)).parse_next(i)?;
    let mut form = DativeBondDsl::default();
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

fn inline_constraint_tag(c: &DativeBondConstraint) -> Option<&'static str> {
    match c {
        DativeBondConstraint::RingCount(_) => Some("#R"),
        DativeBondConstraint::RingSize(_) => Some("#r"),
        _ => None,
    }
}

fn constraint_tag(c: &DativeBondConstraint) -> &'static str {
    inline_constraint_tag(c).expect("non-inline-capable dative constraint produced by parser")
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DativePredicate {
    Constraint(DativeBondConstraint),
}

fn dative_predicate(i: &mut &str) -> PResult<DativePredicate> {
    let start = *i;
    let prefix: &str = take(2usize).parse_next(i)?;
    match prefix {
        "#R" => ring_count
            .map(|v| DativePredicate::Constraint(DativeBondConstraint::RingCount(v)))
            .parse_next(i),
        "#r" => optional_value
            .map(|v| DativePredicate::Constraint(DativeBondConstraint::RingSize(v)))
            .parse_next(i),
        p if p.starts_with('#') => Err(ErrMode::Cut(ParseError::UnknownDativePredicate(
            p.to_string(),
        ))),
        _ => Err(ErrMode::Cut(ParseError::TrailingInput(start.to_string()))),
    }
}

fn apply_predicates(
    form: &mut DativeBondDsl,
    preds: Vec<DativePredicate>,
) -> Result<(), ParseError> {
    let ast = &mut form.0;
    for pred in preds {
        let DativePredicate::Constraint(c) = pred;
        let tag = constraint_tag(&c);
        if ast
            .constraints
            .iter()
            .any(|existing| inline_constraint_tag(existing) == Some(tag))
        {
            return Err(ParseError::DuplicateDativePredicate(tag.to_string()));
        }
        ast.constraints.add(c);
    }
    Ok(())
}

// -- Format --------------------

fn fmt_constraint(f: &mut fmt::Formatter<'_>, c: &DativeBondConstraint) -> fmt::Result {
    match c {
        DativeBondConstraint::RingCount(v) => fmt_ring_count(f, v),
        DativeBondConstraint::RingSize(v) => match v {
            ValueAst::Undetermined => Ok(()),
            ValueAst::Lit(1) => write!(f, "#r"),
            ValueAst::Lit(n) => write!(f, "#r{}", n),
            v => {
                write!(f, "#r")?;
                fmt_value(f, v)
            }
        },
        _ => Ok(()),
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::constraint::DativeBondConstraints;
    use crate::ast::dative::DativeDirection;
    use crate::ast::value::{Expr, RelOp};

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", DativeBondDsl(DativeBondAst::default()))]
    #[case::whitespace("   ", DativeBondDsl(DativeBondAst::default()))]
    #[case::ring_count("#R2", DativeBondDsl(DativeBondAst { direction: DativeDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Lit(2))]) }))]
    #[case::ring_bare("#R", DativeBondDsl(DativeBondAst { direction: DativeDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Lit(1))]) }))]
    #[case::ring_plus("#R+", DativeBondDsl(DativeBondAst { direction: DativeDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("r".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))))]) }))]
    #[case::ring_undetermined("#R*", DativeBondDsl(DativeBondAst { direction: DativeDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Undetermined)]) }))]
    #[case::ring_size("#r6", DativeBondDsl(DativeBondAst { direction: DativeDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingSize(ValueAst::Lit(6))]) }))]
    #[case::ring_size_bare("#r", DativeBondDsl(DativeBondAst { direction: DativeDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingSize(ValueAst::Lit(1))]) }))]
    #[case::ring_count_and_size("#R2#r6", DativeBondDsl(DativeBondAst { direction: DativeDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Lit(2)), DativeBondConstraint::RingSize(ValueAst::Lit(6))]) }))]
    fn test_parse_dative(#[case] input: &str, #[case] expected: DativeBondDsl) {
        let result = dative.parse(input);
        assert!(result.is_ok(), "{:?} should succeed, got {:?}", input, result.clone().unwrap_err());
        let form = result.unwrap();
        assert_eq!(form, expected);
    }

    #[rstest]
    #[case::unknown("#x", ParseError::UnknownDativePredicate("#x".to_string()))]
    #[case::unknown_c("#c", ParseError::UnknownDativePredicate("#c".to_string()))]
    #[case::dup_ring("#R1#R2", ParseError::DuplicateDativePredicate("#R".to_string()))]
    #[case::dup_ring_size("#r6#r5", ParseError::DuplicateDativePredicate("#r".to_string()))]
    #[case::trailing("#R2 foo", ParseError::TrailingInput("foo".to_string()))]
    fn test_parse_dative_error(#[case] input: &str, #[case] expected: ParseError) {
        let result = dative.parse(input);
        assert!(result.is_err(), "{:?} should fail", input);
        let err = result.unwrap_err().into_inner();
        assert_eq!(err, expected);
    }

    #[rstest]
    #[case::empty("")]
    #[case::ring_count("#R2")]
    #[case::ring_size("#r6")]
    #[case::both("#R2#r6")]
    fn test_dative_roundtrip(#[case] input: &str) {
        let form: DativeBondDsl = input.parse().unwrap();
        let rendered = form.to_string();
        let reparsed: DativeBondDsl = rendered.parse().unwrap();
        assert_eq!(form, reparsed);
    }

    #[rstest]
    fn test_dative_dsl_to_ast_passthrough() {
        let dsl = DativeBondDsl(DativeBondAst {
            direction: DativeDirection::Forward,
            constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(
                ValueAst::Lit(2),
            )]),
        });
        let cfg = DativeBondAstConfig::zeroed();
        let ast = dsl.to_ast(&cfg).unwrap();
        assert_eq!(
            ast.constraints,
            DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Lit(2))])
        );
    }
}
