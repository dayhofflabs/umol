//! Aromatic-system-string DSL.

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
use super::predicates::{
    apply_spin_pair, charge, fmt_charge, fmt_spin_pair, fmt_value, optional_value, SpinPredicate,
};
use crate::ast::aromatic::AromaticSystemAst;
use crate::ast::constraint::AromaticSystemConstraint;
use crate::ast::value::ValueAst;

/// AST-layer bundle pairing an `AromaticSystemAst` with the system-level
/// constraints that round-trip through the aromatic-string DSL.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AromaticSystemTypeAst {
    pub ast: AromaticSystemAst,
    pub constraints: Vec<AromaticSystemConstraint>,
}

impl AromaticSystemTypeAst {
    pub fn new(ast: AromaticSystemAst) -> Self {
        Self {
            ast,
            constraints: Vec::new(),
        }
    }

    pub fn with_constraints(
        ast: AromaticSystemAst,
        constraints: Vec<AromaticSystemConstraint>,
    ) -> Self {
        Self { ast, constraints }
    }
}

impl FromStr for AromaticSystemTypeAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_aromatic(s)
    }
}

impl Display for AromaticSystemTypeAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_aromatic_ast(f, &self.ast)?;
        for c in &self.constraints {
            fmt_constraint(f, c)?;
        }
        Ok(())
    }
}

impl<'de> FromEdn<'de> for AromaticSystemTypeAst {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s.parse().map_err(|e| DeError::subgrammar("aromatic", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for AromaticSystemTypeAst {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

pub fn parse_aromatic(input: &str) -> Result<AromaticSystemTypeAst, ParseError> {
    aromatic.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn aromatic(i: &mut &str) -> PResult<AromaticSystemTypeAst> {
    multispace0.parse_next(i)?;
    let preds: Vec<AromaticPredicate> =
        repeat(0.., terminated(aromatic_predicate, multispace0)).parse_next(i)?;
    let mut form = AromaticSystemTypeAst::default();
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

fn apply_predicates(
    form: &mut AromaticSystemTypeAst,
    preds: Vec<AromaticPredicate>,
) -> Result<(), ParseError> {
    let ast = &mut form.ast;
    for pred in preds {
        match pred {
            AromaticPredicate::Charge(v) => {
                if !matches!(ast.charge, ValueAst::Undetermined) {
                    return Err(ParseError::DuplicateAromaticPredicate("#c".to_string()));
                }
                ast.charge = v;
            }
            AromaticPredicate::Spin(sp) => {
                apply_spin_pair(&mut ast.spin, sp, ParseError::DuplicateAromaticPredicate)?;
            }
            AromaticPredicate::Electrons(v) => {
                if !matches!(ast.electrons, ValueAst::Undetermined) {
                    return Err(ParseError::DuplicateAromaticPredicate("#e".to_string()));
                }
                ast.electrons = v;
            }
        }
    }
    Ok(())
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AromaticPredicate {
    Charge(ValueAst),
    Spin(SpinPredicate),
    Electrons(ValueAst),
}

fn aromatic_predicate(i: &mut &str) -> PResult<AromaticPredicate> {
    let start = *i;
    let prefix: &str = take(2usize).parse_next(i)?;
    match prefix {
        "#c" => charge.map(AromaticPredicate::Charge).parse_next(i),
        "#u" => optional_value
            .map(|v| AromaticPredicate::Spin(SpinPredicate::Unpaired(v)))
            .parse_next(i),
        "#s" => optional_value
            .map(|v| AromaticPredicate::Spin(SpinPredicate::Multiplicity(v)))
            .parse_next(i),
        "#e" => optional_value
            .map(AromaticPredicate::Electrons)
            .parse_next(i),
        p if p.starts_with('#') => Err(ErrMode::Cut(ParseError::UnknownAromaticPredicate(
            p.to_string(),
        ))),
        _ => Err(ErrMode::Cut(ParseError::TrailingInput(start.to_string()))),
    }
}

fn fmt_aromatic_ast(f: &mut fmt::Formatter<'_>, ast: &AromaticSystemAst) -> fmt::Result {
    fmt_charge(f, &ast.charge)?;
    fmt_spin_pair(f, &ast.spin)?;
    fmt_electrons(f, &ast.electrons)
}

fn fmt_electrons(f: &mut fmt::Formatter<'_>, v: &ValueAst) -> fmt::Result {
    match v {
        ValueAst::Undetermined => Ok(()),
        ValueAst::Lit(1) => write!(f, "#e"),
        ValueAst::Lit(n) => write!(f, "#e{}", n),
        v => {
            write!(f, "#e")?;
            fmt_value(f, v)
        }
    }
}

fn fmt_constraint(_f: &mut fmt::Formatter<'_>, c: &AromaticSystemConstraint) -> fmt::Result {
    match c {
        AromaticSystemConstraint::Atoms(_)
        | AromaticSystemConstraint::Contains(_)
        | AromaticSystemConstraint::ContainsAll(_)
        | AromaticSystemConstraint::AllAtoms(_)
        | AromaticSystemConstraint::AnyAtom(_) => {
            unreachable!("molecule-scope aromatic constraint in aromatic-system DSL")
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::spin::SpinStateAst;

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", AromaticSystemTypeAst::new(AromaticSystemAst::default()))]
    #[case::whitespace("   ", AromaticSystemTypeAst::new(AromaticSystemAst::default()))]
    #[case::charge_pos("#c+1", AromaticSystemTypeAst::new(AromaticSystemAst { charge: ValueAst::Lit(1), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: Vec::new() }))]
    #[case::charge_neg("#c-2", AromaticSystemTypeAst::new(AromaticSystemAst { charge: ValueAst::Lit(-2), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: Vec::new() }))]
    #[case::charge_plus_only("#c+", AromaticSystemTypeAst::new(AromaticSystemAst { charge: ValueAst::Lit(1), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: Vec::new() }))]
    #[case::charge_minus_only("#c-", AromaticSystemTypeAst::new(AromaticSystemAst { charge: ValueAst::Lit(-1), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: Vec::new() }))]
    #[case::electrons("#e6", AromaticSystemTypeAst::new(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Lit(6), constraints: Vec::new() }))]
    #[case::electrons_bare("#e", AromaticSystemTypeAst::new(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Lit(1), constraints: Vec::new() }))]
    #[case::electrons_wild("#e*", AromaticSystemTypeAst::new(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: Vec::new() }))]
    #[case::unpaired("#u1", AromaticSystemTypeAst::new(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Undetermined }, electrons: ValueAst::Undetermined, constraints: Vec::new() }))]
    #[case::mult("#s2", AromaticSystemTypeAst::new(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(2) }, electrons: ValueAst::Undetermined, constraints: Vec::new() }))]
    #[case::charge_electrons("#c+#e6", AromaticSystemTypeAst::new(AromaticSystemAst { charge: ValueAst::Lit(1), spin: SpinStateAst::default(), electrons: ValueAst::Lit(6), constraints: Vec::new() }))]
    #[case::full("#c0#u0#s1#e6", AromaticSystemTypeAst::new(AromaticSystemAst { charge: ValueAst::Lit(0), spin: SpinStateAst::new(0, 1), electrons: ValueAst::Lit(6), constraints: Vec::new() }))]
    fn test_parse_aromatic(#[case] input: &str, #[case] expected: AromaticSystemTypeAst) {
        let result = aromatic.parse(input);
        assert!(result.is_ok(), "{:?} should succeed, got {:?}", input, result.clone().unwrap_err());
        let form = result.unwrap();
        assert_eq!(form, expected);
    }

    #[rstest]
    #[case::unknown("#x", ParseError::UnknownAromaticPredicate("#x".to_string()))]
    #[case::unknown_a("#a", ParseError::UnknownAromaticPredicate("#a".to_string()))]
    #[case::dup_charge("#c+#c-", ParseError::DuplicateAromaticPredicate("#c".to_string()))]
    #[case::dup_electrons("#e6#e4", ParseError::DuplicateAromaticPredicate("#e".to_string()))]
    #[case::dup_unpaired("#u1#u2", ParseError::DuplicateAromaticPredicate("#u".to_string()))]
    #[case::dup_multiplicity("#s1#s2", ParseError::DuplicateAromaticPredicate("#s".to_string()))]
    #[case::trailing("#c+ foo", ParseError::TrailingInput("foo".to_string()))]
    fn test_parse_aromatic_error(#[case] input: &str, #[case] expected: ParseError) {
        let result = aromatic.parse(input);
        assert!(
            result.is_err(),
            "{:?} should fail, got {:?}",
            input,
            result.unwrap()
        );
        let err = result.unwrap_err().into_inner();
        assert_eq!(err, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(AromaticSystemTypeAst::default(), "")]
    #[case::charge_one(AromaticSystemTypeAst::new(AromaticSystemAst { charge: ValueAst::Lit(1), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: Vec::new() }), "#c+")]
    #[case::charge_neg_two(AromaticSystemTypeAst::new(AromaticSystemAst { charge: ValueAst::Lit(-2), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: Vec::new() }), "#c-2")]
    #[case::electrons_six(AromaticSystemTypeAst::new(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Lit(6), constraints: Vec::new() }), "#e6")]
    #[case::electrons_one(AromaticSystemTypeAst::new(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Lit(1), constraints: Vec::new() }), "#e")]
    #[case::full(AromaticSystemTypeAst::new(AromaticSystemAst { charge: ValueAst::Lit(0), spin: SpinStateAst::new(0, 1), electrons: ValueAst::Lit(6), constraints: Vec::new() }), "#e6")]
    fn test_display_aromatic(#[case] form: AromaticSystemTypeAst, #[case] expected: &str) {
        assert_eq!(form.to_string(), expected);
    }

    #[rstest]
    #[case::empty("")]
    #[case::charge("#c+1")]
    #[case::electrons("#e6")]
    #[case::unpaired("#u2")]
    #[case::explicit_mult("#s2")]
    fn test_aromatic_roundtrip(#[case] input: &str) {
        let form: AromaticSystemTypeAst = input.parse().unwrap();
        let rendered = form.to_string();
        let reparsed: AromaticSystemTypeAst = rendered.parse().unwrap();
        assert_eq!(form, reparsed);
    }
}
