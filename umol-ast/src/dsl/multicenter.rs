//! Multicenter-bond-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, EdnError, EdnStreamDeserializer, FromEdn, ToEdn};
use winnow::ascii::multispace0;
use winnow::combinator::{repeat, terminated};
use winnow::error::ErrMode;
use winnow::token::take;
use winnow::Parser;

use super::error::{PResult, ParseError};
use super::predicates::{
    apply_spin_pair, charge, fmt_charge, fmt_spin_pair, lower_spin, optional_value,
    raise_spin, SpinPredicate,
};
use super::value::fmt_value;
use crate::ast::config::{MulticenterBondAstConfig, NumericMode};
use crate::ast::multicenter::MulticenterBondAst;
use crate::ast::traits::{FromAst, ToAst};
use crate::ast::value::ValueAst;

/// Surface DSL wrapper around `MulticenterBondAst`. Parses and renders the
/// multicenter-bond-string form. All `MulticenterBondConstraint` variants are
/// molecule-scope, so nothing from the constraint vec serializes inline.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MulticenterBondDsl(pub MulticenterBondAst);

impl FromStr for MulticenterBondDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_multicenter(s)
    }
}

impl Display for MulticenterBondDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_multicenter_ast(f, &self.0)
    }
}

impl<'de> FromEdn<'de> for MulticenterBondDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s.parse().map_err(|e| DeError::subgrammar("multicenter", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        EdnStreamDeserializer::new(input).read_subgrammar_all("multicenter")
    }
}

impl ToEdn for MulticenterBondDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

impl FromAst<MulticenterBondAst> for MulticenterBondDsl {
    type Error = ParseError;

    fn from_ast(
        ast: &MulticenterBondAst,
        cfg: &MulticenterBondAstConfig,
    ) -> Result<Self, ParseError> {
        let mut out = ast.clone();
        lower_multicenter(&mut out, cfg);
        Ok(MulticenterBondDsl(out))
    }
}

impl ToAst<MulticenterBondAst> for MulticenterBondDsl {
    type Error = ParseError;

    fn to_ast(&self, cfg: &MulticenterBondAstConfig) -> Result<MulticenterBondAst, ParseError> {
        let mut out = self.0.clone();
        raise_multicenter(&mut out, cfg);
        Ok(out)
    }
}

// -- Parse --------------------

pub fn parse_multicenter(input: &str) -> Result<MulticenterBondDsl, ParseError> {
    multicenter.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn multicenter(i: &mut &str) -> PResult<MulticenterBondDsl> {
    multispace0.parse_next(i)?;
    let preds: Vec<MulticenterPredicate> =
        repeat(0.., terminated(multicenter_predicate, multispace0)).parse_next(i)?;
    let mut form = MulticenterBondDsl::default();
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MulticenterPredicate {
    Charge(ValueAst),
    Spin(SpinPredicate),
    Electrons(ValueAst),
}

fn multicenter_predicate(i: &mut &str) -> PResult<MulticenterPredicate> {
    let start = *i;
    let prefix: &str = take(2usize).parse_next(i)?;
    match prefix {
        "#c" => charge.map(MulticenterPredicate::Charge).parse_next(i),
        "#u" => optional_value
            .map(|v| MulticenterPredicate::Spin(SpinPredicate::Unpaired(v)))
            .parse_next(i),
        "#s" => optional_value
            .map(|v| MulticenterPredicate::Spin(SpinPredicate::Multiplicity(v)))
            .parse_next(i),
        "#e" => optional_value
            .map(MulticenterPredicate::Electrons)
            .parse_next(i),
        p if p.starts_with('#') => Err(ErrMode::Cut(ParseError::UnknownMulticenterPredicate(
            p.to_string(),
        ))),
        _ => Err(ErrMode::Cut(ParseError::TrailingInput(start.to_string()))),
    }
}

fn apply_predicates(
    form: &mut MulticenterBondDsl,
    preds: Vec<MulticenterPredicate>,
) -> Result<(), ParseError> {
    let ast = &mut form.0;
    for pred in preds {
        match pred {
            MulticenterPredicate::Charge(v) => {
                if !matches!(ast.charge, ValueAst::Undetermined) {
                    return Err(ParseError::DuplicateMulticenterPredicate("#c".to_string()));
                }
                ast.charge = v;
            }
            MulticenterPredicate::Spin(sp) => {
                apply_spin_pair(&mut ast.spin, sp, ParseError::DuplicateMulticenterPredicate)?;
            }
            MulticenterPredicate::Electrons(v) => {
                if !matches!(ast.electrons, ValueAst::Undetermined) {
                    return Err(ParseError::DuplicateMulticenterPredicate("#e".to_string()));
                }
                ast.electrons = v;
            }
        }
    }
    Ok(())
}

// -- Format --------------------

fn fmt_multicenter_ast(f: &mut fmt::Formatter<'_>, ast: &MulticenterBondAst) -> fmt::Result {
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

// -- Raise --------------------

fn raise_multicenter(ast: &mut MulticenterBondAst, cfg: &MulticenterBondAstConfig) {
    // Exhaustive destructure: adding a new MulticenterBondAst field is a
    // compile error here, forcing the author to decide how raising should
    // handle it.
    let MulticenterBondAst {
        charge,
        spin,
        electrons,
        constraints: _,
    } = ast;

    if matches!(*charge, ValueAst::Undetermined) {
        *charge = match cfg.charge_mode {
            NumericMode::Zero => ValueAst::Lit(0),
            NumericMode::Required => ValueAst::Undetermined,
        };
    }
    if matches!(*electrons, ValueAst::Undetermined) {
        *electrons = match cfg.electrons_mode {
            NumericMode::Zero => ValueAst::Lit(0),
            NumericMode::Required => ValueAst::Undetermined,
        };
    }
    raise_spin(spin, cfg.unpaired_electrons_mode, cfg.multiplicity_mode);
}

// -- Format --------------------

fn lower_multicenter(ast: &mut MulticenterBondAst, cfg: &MulticenterBondAstConfig) {
    // Exhaustive destructure: adding a new MulticenterBondAst field is a
    // compile error here, forcing the author to decide how lowering should
    // handle it.
    let MulticenterBondAst {
        charge,
        spin,
        electrons,
        constraints: _,
    } = ast;

    if matches!(
        (&cfg.charge_mode, &*charge),
        (NumericMode::Zero, ValueAst::Lit(0))
    ) {
        *charge = ValueAst::Undetermined;
    }
    if matches!(
        (&cfg.electrons_mode, &*electrons),
        (NumericMode::Zero, ValueAst::Lit(0))
    ) {
        *electrons = ValueAst::Undetermined;
    }
    lower_spin(spin, cfg.unpaired_electrons_mode, cfg.multiplicity_mode);
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::constraint::MulticenterBondConstraints;
    use crate::ast::spin::SpinStateAst;

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", MulticenterBondDsl(MulticenterBondAst::default()))]
    #[case::whitespace("   ", MulticenterBondDsl(MulticenterBondAst::default()))]
    #[case::charge_pos("#c+1", MulticenterBondDsl(MulticenterBondAst { charge: ValueAst::Lit(1), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: MulticenterBondConstraints::new() }))]
    #[case::charge_neg("#c-2", MulticenterBondDsl(MulticenterBondAst { charge: ValueAst::Lit(-2), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: MulticenterBondConstraints::new() }))]
    #[case::electrons("#e6", MulticenterBondDsl(MulticenterBondAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Lit(6), constraints: MulticenterBondConstraints::new() }))]
    #[case::electrons_bare("#e", MulticenterBondDsl(MulticenterBondAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Lit(1), constraints: MulticenterBondConstraints::new() }))]
    #[case::unpaired("#u1", MulticenterBondDsl(MulticenterBondAst { charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Undetermined }, electrons: ValueAst::Undetermined, constraints: MulticenterBondConstraints::new() }))]
    #[case::mult("#s2", MulticenterBondDsl(MulticenterBondAst { charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(2) }, electrons: ValueAst::Undetermined, constraints: MulticenterBondConstraints::new() }))]
    #[case::charge_electrons("#c+#e2", MulticenterBondDsl(MulticenterBondAst { charge: ValueAst::Lit(1), spin: SpinStateAst::default(), electrons: ValueAst::Lit(2), constraints: MulticenterBondConstraints::new() }))]
    #[case::full("#c0#u0#s1#e2", MulticenterBondDsl(MulticenterBondAst { charge: ValueAst::Lit(0), spin: SpinStateAst::new(0, 1), electrons: ValueAst::Lit(2), constraints: MulticenterBondConstraints::new() }))]
    fn test_parse_multicenter(#[case] input: &str, #[case] expected: MulticenterBondDsl) {
        let result = multicenter.parse(input);
        assert!(result.is_ok(), "{:?} should succeed, got {:?}", input, result.clone().unwrap_err());
        let form = result.unwrap();
        assert_eq!(form, expected);
    }

    #[rstest]
    #[case::unknown("#x", ParseError::UnknownMulticenterPredicate("#x".to_string()))]
    #[case::unknown_a("#a", ParseError::UnknownMulticenterPredicate("#a".to_string()))]
    #[case::dup_charge("#c+#c-", ParseError::DuplicateMulticenterPredicate("#c".to_string()))]
    #[case::dup_electrons("#e2#e4", ParseError::DuplicateMulticenterPredicate("#e".to_string()))]
    #[case::dup_unpaired("#u1#u2", ParseError::DuplicateMulticenterPredicate("#u".to_string()))]
    #[case::dup_multiplicity("#s1#s2", ParseError::DuplicateMulticenterPredicate("#s".to_string()))]
    #[case::trailing("#c+ foo", ParseError::TrailingInput("foo".to_string()))]
    fn test_parse_multicenter_error(#[case] input: &str, #[case] expected: ParseError) {
        let result = multicenter.parse(input);
        assert!(result.is_err(), "{:?} should fail", input);
        let err = result.unwrap_err().into_inner();
        assert_eq!(err, expected);
    }

    #[rstest]
    #[case::empty("")]
    #[case::charge("#c+1")]
    #[case::electrons("#e6")]
    #[case::unpaired("#u2")]
    #[case::explicit_mult("#s2")]
    fn test_multicenter_roundtrip(#[case] input: &str) {
        let form: MulticenterBondDsl = input.parse().unwrap();
        let rendered = form.to_string();
        let reparsed: MulticenterBondDsl = rendered.parse().unwrap();
        assert_eq!(form, reparsed);
    }

    #[rstest]
    fn test_multicenter_dsl_to_ast_fills_zero_defaults() {
        let dsl = MulticenterBondDsl::default();
        let cfg = MulticenterBondAstConfig::zeroed();
        let ast = dsl.to_ast(&cfg).unwrap();
        assert_eq!(ast.charge, ValueAst::Lit(0));
        assert_eq!(ast.electrons, ValueAst::Lit(0));
        assert_eq!(ast.spin, SpinStateAst::new(0, 1));
    }

    #[rstest]
    fn test_multicenter_dsl_from_ast_strips_zero_defaults() {
        let ast = MulticenterBondAst {
            charge: ValueAst::Lit(0),
            spin: SpinStateAst::new(0, 1),
            electrons: ValueAst::Lit(0),
            constraints: MulticenterBondConstraints::new(),
        };
        let cfg = MulticenterBondAstConfig::zeroed();
        let dsl = MulticenterBondDsl::from_ast(&ast, &cfg).unwrap();
        assert_eq!(dsl.0.charge, ValueAst::Undetermined);
        assert_eq!(dsl.0.electrons, ValueAst::Undetermined);
        assert_eq!(dsl.0.spin, SpinStateAst::default());
    }

    #[rstest]
    #[case::empty(r##""""##)]
    #[case::charge(r##""#c+""##)]
    #[case::full(r##""#c0#u0#s1#e2""##)]
    fn test_multicenter_dsl_from_edn_str_matches_from_edn(#[case] input: &str) {
        let via_stream = MulticenterBondDsl::from_edn_str(input).unwrap();
        let tree = umol_edn::read_string(input).unwrap();
        let via_tree = MulticenterBondDsl::from_edn(&tree).unwrap();
        assert_eq!(via_stream, via_tree);
    }
}
