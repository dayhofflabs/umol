//! Dative-bond-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, EdnError, EdnStreamDeserializer, FromEdn, ToEdn};
use winnow::ascii::multispace0;
use winnow::combinator::{repeat, terminated};
use winnow::error::ErrMode;
use winnow::token::take;
use winnow::Parser;

use super::config::DativeBondDefaults;
use super::error::{PResult, ParseError};
use super::predicates::{fmt_ring_count, optional_value, ring_count};
use super::value::{fmt_value, ValueDsl};
use crate::ast::constraint::DativeBondConstraint;
use crate::ast::dative::DativeBondAst;
use crate::ast::traits::{FromAst, IntoAst};
use crate::ast::value::ValueAst;

/// Surface DSL wrapper around `DativeBondAst`. No leading token; the string
/// form is a sequence of `#…` predicates. Inline-capable constraints from
/// `DativeBondConstraint` are `RingCount` (`#R`) and `RingSize` (`#r`); the
/// remaining variants reference other entities and stay in the molecule
/// constraints container.
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DativeBondDsl(pub DativeBondAst);

impl DativeBondDsl {
    /// Zero-cost reference cast from `&DativeBondAst`. Relies on `repr(transparent)`.
    pub fn from_ref(ast: &DativeBondAst) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(ast as *const DativeBondAst as *const Self) }
    }
}

impl FromStr for DativeBondDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_dative_bond(s)
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

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        EdnStreamDeserializer::new(input).read_subgrammar_all("dative")
    }
}

impl ToEdn for DativeBondDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

impl FromAst<DativeBondAst> for DativeBondDsl {
    type Ctx = DativeBondDefaults;
    type Error = ParseError;

    fn from_ast(ast: &DativeBondAst, _cfg: &Self::Ctx) -> Result<Self, ParseError> {
        Ok(DativeBondDsl(ast.clone()))
    }
}

impl IntoAst<DativeBondAst> for DativeBondDsl {
    type Ctx = DativeBondDefaults;
    type Error = ParseError;

    fn into_ast(self, _cfg: &Self::Ctx) -> Result<DativeBondAst, ParseError> {
        Ok(self.0)
    }
}

// -- Parse --------------------

/// Parse a complete dative-bond-string into a `DativeBondDsl`.
pub fn parse_dative_bond(input: &str) -> Result<DativeBondDsl, ParseError> {
    dative_bond.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn dative_bond(i: &mut &str) -> PResult<DativeBondDsl> {
    multispace0.parse_next(i)?;
    let preds: Vec<DativeBondPredicate> =
        repeat(0.., terminated(dative_bond_predicate, multispace0)).parse_next(i)?;
    let mut form = DativeBondDsl::default();
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

fn constraint_tag(c: &DativeBondConstraint) -> &'static str {
    match c {
        DativeBondConstraint::RingCount(_) => "#R",
        DativeBondConstraint::RingSize(_) => "#r",
    }
}

/// One predicate from a dative-bond-string; the parser yields a `Vec` of
/// these and the applier folds them into the `DativeBondAst`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DativeBondPredicate {
    Constraint(DativeBondConstraint),
}

fn dative_bond_predicate(i: &mut &str) -> PResult<DativeBondPredicate> {
    let start = *i;
    let prefix: &str = take(2usize).parse_next(i)?;
    match prefix {
        "#R" => ring_count
            .map(|v| DativeBondPredicate::Constraint(DativeBondConstraint::RingCount(v)))
            .parse_next(i),
        "#r" => optional_value
            .map(|v| DativeBondPredicate::Constraint(DativeBondConstraint::RingSize(v)))
            .parse_next(i),
        p if p.starts_with('#') => Err(ErrMode::Cut(ParseError::UnknownDativeBondPredicate(
            p.to_string(),
        ))),
        _ => Err(ErrMode::Cut(ParseError::TrailingInput(start.to_string()))),
    }
}

fn apply_predicates(
    form: &mut DativeBondDsl,
    preds: Vec<DativeBondPredicate>,
) -> Result<(), ParseError> {
    let ast = &mut form.0;
    for pred in preds {
        let DativeBondPredicate::Constraint(c) = pred;
        let tag = constraint_tag(&c);
        if ast
            .constraints
            .iter()
            .any(|existing| constraint_tag(existing) == tag)
        {
            return Err(ParseError::DuplicateDativeBondPredicate(tag.to_string()));
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
    }
}

// -- Constraint DSL -------------------

/// Surface DSL wrapper around the narrow `DativeBondConstraint`. EDN form is
/// a single-key map keyed by the constraint kind: `{:ring-count <value>}` or
/// `{:ring-size <value>}`. Ref-bearing variants moved to
/// [`super::relational::RelationalConstraintDsl`].
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum DativeBondConstraintDsl {
    RingCount(ValueAst),
    RingSize(ValueAst),
}

impl<'de> FromEdn<'de> for DativeBondConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let Edn::Map(m) = edn else {
            return Err(DeError::TypeMismatch {
                expected: "dative-bond-constraint single-key map",
                got: edn.kind(),
                path: Vec::new(),
            });
        };
        if m.len() != 1 {
            return Err(DeError::Custom(format!(
                "dative-bond-constraint must have exactly one key, got {}",
                m.len()
            )));
        }
        let (k, v) = m.iter().next().unwrap();
        let Edn::Keyword(key) = k else {
            return Err(DeError::TypeMismatch {
                expected: "keyword key",
                got: k.kind(),
                path: vec!["dative-bond-constraint".into()],
            });
        };
        Ok(match key.name() {
            "ring-count" => Self::RingCount(ValueDsl::from_edn(v)?.into_ast(&()).unwrap()),
            "ring-size" => Self::RingSize(ValueDsl::from_edn(v)?.into_ast(&()).unwrap()),
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["dative-bond-constraint".into()],
                });
            }
        })
    }
}

impl ToEdn for DativeBondConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        let (key, value) = match self {
            Self::RingCount(v) => ("ring-count", ValueDsl::from_ast(v, &()).unwrap().to_edn()),
            Self::RingSize(v) => ("ring-size", ValueDsl::from_ast(v, &()).unwrap().to_edn()),
        };
        let mut m = umol_edn::EdnMap::with_capacity(1);
        m.insert(Edn::Keyword(umol_edn::EdnKeyword::owned(key.into())), value);
        Edn::Map(m)
    }
}

impl DativeBondConstraintDsl {
    /// Build from the narrow inline AST form.
    pub(crate) fn from_ast(c: &DativeBondConstraint) -> Self {
        match c {
            DativeBondConstraint::RingCount(v) => Self::RingCount(v.clone()),
            DativeBondConstraint::RingSize(v) => Self::RingSize(v.clone()),
        }
    }

    /// Convert into the narrow inline AST form.
    pub(crate) fn into_ast(self) -> DativeBondConstraint {
        match self {
            Self::RingCount(v) => DativeBondConstraint::RingCount(v),
            Self::RingSize(v) => DativeBondConstraint::RingSize(v),
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::constraint::DativeBondConstraints;
    use crate::ast::dative::DativeBondDirection;
    use crate::ast::value::{Expr, RelOp};

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", DativeBondDsl(DativeBondAst::default()))]
    #[case::whitespace("   ", DativeBondDsl(DativeBondAst::default()))]
    #[case::ring_count("#R2", DativeBondDsl(DativeBondAst { direction: DativeBondDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Lit(2))]) }))]
    #[case::ring_bare("#R", DativeBondDsl(DativeBondAst { direction: DativeBondDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Lit(1))]) }))]
    #[case::ring_plus("#R+", DativeBondDsl(DativeBondAst { direction: DativeBondDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("r".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))))]) }))]
    #[case::ring_undetermined("#R*", DativeBondDsl(DativeBondAst { direction: DativeBondDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Undetermined)]) }))]
    #[case::ring_size("#r6", DativeBondDsl(DativeBondAst { direction: DativeBondDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingSize(ValueAst::Lit(6))]) }))]
    #[case::ring_size_bare("#r", DativeBondDsl(DativeBondAst { direction: DativeBondDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingSize(ValueAst::Lit(1))]) }))]
    #[case::ring_count_and_size("#R2#r6", DativeBondDsl(DativeBondAst { direction: DativeBondDirection::Forward, constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Lit(2)), DativeBondConstraint::RingSize(ValueAst::Lit(6))]) }))]
    fn test_parse_dative(#[case] input: &str, #[case] expected: DativeBondDsl) {
        let result = dative_bond.parse(input);
        assert!(result.is_ok(), "{:?} should succeed, got {:?}", input, result.clone().unwrap_err());
        let form = result.unwrap();
        assert_eq!(form, expected);
    }

    #[rstest]
    #[case::unknown("#x", ParseError::UnknownDativeBondPredicate("#x".to_string()))]
    #[case::unknown_c("#c", ParseError::UnknownDativeBondPredicate("#c".to_string()))]
    #[case::dup_ring("#R1#R2", ParseError::DuplicateDativeBondPredicate("#R".to_string()))]
    #[case::dup_ring_size("#r6#r5", ParseError::DuplicateDativeBondPredicate("#r".to_string()))]
    #[case::trailing("#R2 foo", ParseError::TrailingInput("foo".to_string()))]
    fn test_parse_dative_error(#[case] input: &str, #[case] expected: ParseError) {
        let result = dative_bond.parse(input);
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
            direction: DativeBondDirection::Forward,
            constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(
                ValueAst::Lit(2),
            )]),
        });
        let cfg = DativeBondDefaults::zeroed();
        let ast = dsl.into_ast(&cfg).unwrap();
        assert_eq!(
            ast.constraints,
            DativeBondConstraints::from_iter([DativeBondConstraint::RingCount(ValueAst::Lit(2))])
        );
    }

    #[rstest]
    #[case::empty(r##""""##)]
    #[case::ring_count(r##""#R2""##)]
    #[case::ring_count_and_size(r##""#R2#r6""##)]
    fn test_dative_dsl_from_edn_str_matches_from_edn(#[case] input: &str) {
        let via_stream = DativeBondDsl::from_edn_str(input).unwrap();
        let tree = umol_edn::read_string(input).unwrap();
        let via_tree = DativeBondDsl::from_edn(&tree).unwrap();
        assert_eq!(via_stream, via_tree);
    }

    // -- DativeBondConstraintDsl ----------------

    #[rustfmt::skip]
    #[rstest]
    #[case::ring_count(DativeBondConstraint::RingCount(ValueAst::Lit(2)), "{:ring-count 2}")]
    #[case::ring_size(DativeBondConstraint::RingSize(ValueAst::Lit(6)), "{:ring-size 6}")]
    fn test_dative_bond_constraint_dsl_roundtrip(
        #[case] input: DativeBondConstraint,
        #[case] edn_source: &str,
    ) {
        let dsl = DativeBondConstraintDsl::from_ast(&input);
        let edn = dsl.clone().to_edn();
        let expected = umol_edn::read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = DativeBondConstraintDsl::from_edn(&edn).unwrap();
        let back = parsed.into_ast();
        assert_eq!(back, input, "parse-back mismatch");
    }

    #[rstest]
    fn test_dative_bond_constraint_dsl_rejects_wrong_shape() {
        let err = DativeBondConstraintDsl::from_edn(&Edn::Int(3)).unwrap_err();
        assert!(matches!(err, DeError::TypeMismatch { .. }));
    }

    #[rstest]
    fn test_dative_bond_constraint_dsl_rejects_unknown_key() {
        let edn = umol_edn::read_string("{:bogus 1}").unwrap();
        let err = DativeBondConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::UnknownField { .. }));
    }
}
