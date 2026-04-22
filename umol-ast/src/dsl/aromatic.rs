//! Aromatic-system-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::mem;
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
use crate::ast::config::{
    AromaticSystemAstConfig, MultiplicityMode, NumericMode, UnpairedElectronsMode,
};
use crate::ast::spin::SpinStateAst;
use crate::ast::traits::{FromAst, ToAst};
use crate::ast::value::ValueAst;

/// Surface DSL wrapper around `AromaticSystemAst`. Parses and renders the
/// aromatic-system-string form. All `AromaticSystemConstraint` variants are
/// molecule-scope, so nothing from the constraint vec serializes inline.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AromaticSystemDsl(pub AromaticSystemAst);

impl FromStr for AromaticSystemDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_aromatic(s)
    }
}

impl Display for AromaticSystemDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_aromatic_ast(f, &self.0)
    }
}

impl<'de> FromEdn<'de> for AromaticSystemDsl {
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

impl ToEdn for AromaticSystemDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

impl FromAst<AromaticSystemAst> for AromaticSystemDsl {
    type Error = ParseError;

    fn from_ast(
        ast: &AromaticSystemAst,
        cfg: &AromaticSystemAstConfig,
    ) -> Result<Self, ParseError> {
        let mut out = ast.clone();
        lower_aromatic(&mut out, cfg);
        Ok(AromaticSystemDsl(out))
    }
}

impl ToAst<AromaticSystemAst> for AromaticSystemDsl {
    type Error = ParseError;

    fn to_ast(&self, cfg: &AromaticSystemAstConfig) -> Result<AromaticSystemAst, ParseError> {
        let mut out = self.0.clone();
        raise_aromatic(&mut out, cfg);
        Ok(out)
    }
}

fn raise_aromatic(ast: &mut AromaticSystemAst, cfg: &AromaticSystemAstConfig) {
    if matches!(ast.charge, ValueAst::Undetermined) {
        ast.charge = match cfg.charge_mode {
            NumericMode::Zero => ValueAst::Lit(0),
            NumericMode::Required => ValueAst::Undetermined,
        };
    }
    if matches!(ast.electrons, ValueAst::Undetermined) {
        ast.electrons = match cfg.electrons_mode {
            NumericMode::Zero => ValueAst::Lit(0),
            NumericMode::Required => ValueAst::Undetermined,
        };
    }
    raise_spin(&mut ast.spin, cfg);
}

fn raise_spin(spin: &mut SpinStateAst, cfg: &AromaticSystemAstConfig) {
    let u = mem::replace(&mut spin.unpaired, ValueAst::Undetermined);
    let m = mem::replace(&mut spin.multiplicity, ValueAst::Undetermined);
    let resolved_u = if matches!(u, ValueAst::Undetermined) {
        match cfg.unpaired_electrons_mode {
            UnpairedElectronsMode::Zero => ValueAst::Lit(0),
            UnpairedElectronsMode::Required => ValueAst::Undetermined,
            UnpairedElectronsMode::Derived => match &m {
                ValueAst::Lit(mm) => ValueAst::Lit(mm - 1),
                _ => ValueAst::Undetermined,
            },
        }
    } else {
        u
    };
    let resolved_m = if matches!(m, ValueAst::Undetermined) {
        match cfg.multiplicity_mode {
            MultiplicityMode::Required => ValueAst::Undetermined,
            MultiplicityMode::Derived => match &resolved_u {
                ValueAst::Lit(uu) => ValueAst::Lit(uu + 1),
                _ => ValueAst::Undetermined,
            },
        }
    } else {
        m
    };
    spin.unpaired = resolved_u;
    spin.multiplicity = resolved_m;
}

fn lower_aromatic(ast: &mut AromaticSystemAst, cfg: &AromaticSystemAstConfig) {
    if matches!(
        (&cfg.charge_mode, &ast.charge),
        (NumericMode::Zero, ValueAst::Lit(0))
    ) {
        ast.charge = ValueAst::Undetermined;
    }
    if matches!(
        (&cfg.electrons_mode, &ast.electrons),
        (NumericMode::Zero, ValueAst::Lit(0))
    ) {
        ast.electrons = ValueAst::Undetermined;
    }
    lower_spin(&mut ast.spin, cfg);
}

fn lower_spin(spin: &mut SpinStateAst, cfg: &AromaticSystemAstConfig) {
    if let (ValueAst::Lit(uu), ValueAst::Lit(mm)) = (&spin.unpaired, &spin.multiplicity) {
        let derived = *mm == uu + 1;
        let strip_u = match cfg.unpaired_electrons_mode {
            UnpairedElectronsMode::Zero => *uu == 0,
            UnpairedElectronsMode::Derived => {
                derived && matches!(cfg.multiplicity_mode, MultiplicityMode::Derived)
            }
            UnpairedElectronsMode::Required => false,
        };
        let strip_m = matches!(cfg.multiplicity_mode, MultiplicityMode::Derived) && derived;
        if strip_u {
            spin.unpaired = ValueAst::Undetermined;
        }
        if strip_m {
            spin.multiplicity = ValueAst::Undetermined;
        }
    }
}

pub fn parse_aromatic(input: &str) -> Result<AromaticSystemDsl, ParseError> {
    aromatic.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn aromatic(i: &mut &str) -> PResult<AromaticSystemDsl> {
    multispace0.parse_next(i)?;
    let preds: Vec<AromaticPredicate> =
        repeat(0.., terminated(aromatic_predicate, multispace0)).parse_next(i)?;
    let mut form = AromaticSystemDsl::default();
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

fn apply_predicates(
    form: &mut AromaticSystemDsl,
    preds: Vec<AromaticPredicate>,
) -> Result<(), ParseError> {
    let ast = &mut form.0;
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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::spin::SpinStateAst;

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", AromaticSystemDsl(AromaticSystemAst::default()))]
    #[case::whitespace("   ", AromaticSystemDsl(AromaticSystemAst::default()))]
    #[case::charge_pos("#c+1", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Lit(1), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: Vec::new() }))]
    #[case::charge_neg("#c-2", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Lit(-2), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: Vec::new() }))]
    #[case::charge_plus_only("#c+", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Lit(1), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: Vec::new() }))]
    #[case::charge_minus_only("#c-", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Lit(-1), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: Vec::new() }))]
    #[case::electrons("#e6", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Lit(6), constraints: Vec::new() }))]
    #[case::electrons_bare("#e", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Lit(1), constraints: Vec::new() }))]
    #[case::electrons_wild("#e*", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: Vec::new() }))]
    #[case::unpaired("#u1", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Undetermined }, electrons: ValueAst::Undetermined, constraints: Vec::new() }))]
    #[case::mult("#s2", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(2) }, electrons: ValueAst::Undetermined, constraints: Vec::new() }))]
    #[case::charge_electrons("#c+#e6", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Lit(1), spin: SpinStateAst::default(), electrons: ValueAst::Lit(6), constraints: Vec::new() }))]
    #[case::full("#c0#u0#s1#e6", AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Lit(0), spin: SpinStateAst::new(0, 1), electrons: ValueAst::Lit(6), constraints: Vec::new() }))]
    fn test_parse_aromatic(#[case] input: &str, #[case] expected: AromaticSystemDsl) {
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
    #[case::empty(AromaticSystemDsl::default(), "")]
    #[case::charge_one(AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Lit(1), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: Vec::new() }), "#c+")]
    #[case::charge_neg_two(AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Lit(-2), spin: SpinStateAst::default(), electrons: ValueAst::Undetermined, constraints: Vec::new() }), "#c-2")]
    #[case::electrons_six(AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Lit(6), constraints: Vec::new() }), "#e6")]
    #[case::electrons_one(AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Undetermined, spin: SpinStateAst::default(), electrons: ValueAst::Lit(1), constraints: Vec::new() }), "#e")]
    #[case::full(AromaticSystemDsl(AromaticSystemAst { charge: ValueAst::Lit(0), spin: SpinStateAst::new(0, 1), electrons: ValueAst::Lit(6), constraints: Vec::new() }), "#c0#u0#s#e6")]
    fn test_display_aromatic(#[case] form: AromaticSystemDsl, #[case] expected: &str) {
        assert_eq!(form.to_string(), expected);
    }

    #[rstest]
    #[case::empty("")]
    #[case::charge("#c+1")]
    #[case::electrons("#e6")]
    #[case::unpaired("#u2")]
    #[case::explicit_mult("#s2")]
    fn test_aromatic_roundtrip(#[case] input: &str) {
        let form: AromaticSystemDsl = input.parse().unwrap();
        let rendered = form.to_string();
        let reparsed: AromaticSystemDsl = rendered.parse().unwrap();
        assert_eq!(form, reparsed);
    }

    #[rstest]
    fn test_aromatic_dsl_to_ast_fills_zero_defaults() {
        let dsl = AromaticSystemDsl::default();
        let cfg = AromaticSystemAstConfig::zeroed();
        let ast = dsl.to_ast(&cfg).unwrap();
        assert_eq!(ast.charge, ValueAst::Lit(0));
        assert_eq!(ast.electrons, ValueAst::Lit(0));
        assert_eq!(ast.spin, SpinStateAst::new(0, 1));
    }

    #[rstest]
    fn test_aromatic_dsl_from_ast_strips_zero_defaults() {
        let ast = AromaticSystemAst {
            charge: ValueAst::Lit(0),
            spin: SpinStateAst::new(0, 1),
            electrons: ValueAst::Lit(0),
            constraints: Vec::new(),
        };
        let cfg = AromaticSystemAstConfig::zeroed();
        let dsl = AromaticSystemDsl::from_ast(&ast, &cfg).unwrap();
        assert_eq!(dsl.0.charge, ValueAst::Undetermined);
        assert_eq!(dsl.0.electrons, ValueAst::Undetermined);
        assert_eq!(dsl.0.spin, SpinStateAst::default());
    }

    #[rstest]
    fn test_aromatic_dsl_roundtrip_zeroed() {
        let input = AromaticSystemDsl::default();
        let cfg = AromaticSystemAstConfig::zeroed();
        let raised = input.to_ast(&cfg).unwrap();
        let lowered = AromaticSystemDsl::from_ast(&raised, &cfg).unwrap();
        assert_eq!(input, lowered);
    }
}
