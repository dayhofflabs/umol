//! Aromatic-system-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn};
use winnow::ascii::multispace0;
use winnow::combinator::{repeat, terminated};
use winnow::error::ErrMode;
use winnow::token::take;
use winnow::Parser;

use super::config::{AromaticSystemDefaults, NumericDefault};
use super::electrons::{electron_counts, fmt_electron_counts};
use super::error::{PResult, ParseError};
use super::predicates::{
    apply_spin_pair, charge, fmt_charge, fmt_spin_pair, lower_spin, optional_value, raise_spin,
    SpinPredicate,
};
use super::value::{fmt_value, ValueDsl};
use crate::ast::aromatic::AromaticSystemAst;
use crate::ast::constraint::AromaticSystemConstraint;
use crate::ast::traits::{FromAst, IntoAst};
use crate::ast::value::ValueAst;

/// Surface DSL wrapper around `AromaticSystemAst`.
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AromaticSystemDsl(pub AromaticSystemAst);

impl AromaticSystemDsl {
    /// Zero-cost reference cast from `&AromaticSystemAst`. Relies on `repr(transparent)`.
    pub fn from_ref(ast: &AromaticSystemAst) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(ast as *const AromaticSystemAst as *const Self) }
    }
}

impl From<AromaticSystemAst> for AromaticSystemDsl {
    fn from(ast: AromaticSystemAst) -> Self {
        Self(ast)
    }
}

impl FromStr for AromaticSystemDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_aromatic_system(s)
    }
}

impl Display for AromaticSystemDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_aromatic_system_ast(f, &self.0)
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

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        EdnStreamDeserializer::new(input).read_subgrammar_all("aromatic")
    }
}

impl ToEdn for AromaticSystemDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

impl FromAst<AromaticSystemAst> for AromaticSystemDsl {
    type Ctx = AromaticSystemDefaults;

    fn from_ast(ast: &AromaticSystemAst, cfg: &Self::Ctx) -> Self {
        let mut out = ast.clone();
        lower_aromatic_system(&mut out, cfg);
        AromaticSystemDsl(out)
    }
}

impl IntoAst<AromaticSystemAst> for AromaticSystemDsl {
    type Ctx = AromaticSystemDefaults;

    fn into_ast(mut self, cfg: &Self::Ctx) -> AromaticSystemAst {
        raise_aromatic_system(&mut self.0, cfg);
        self.0
    }
}

impl FromStr for AromaticSystemAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(AromaticSystemDsl::from_str(s)?.into_ast(&AromaticSystemDefaults::default()))
    }
}

impl Display for AromaticSystemAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        AromaticSystemDsl::from_ref(self).fmt(f)
    }
}

impl<'de> FromEdn<'de> for AromaticSystemAst {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        Ok(AromaticSystemDsl::from_edn(edn)?.into_ast(&AromaticSystemDefaults::default()))
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        Ok(AromaticSystemDsl::from_edn_str(input)?.into_ast(&AromaticSystemDefaults::default()))
    }
}

impl ToEdn for AromaticSystemAst {
    fn to_edn(&self) -> Edn<'static> {
        AromaticSystemDsl::from_ref(self).to_edn()
    }
}

pub fn parse_aromatic_system(input: &str) -> Result<AromaticSystemDsl, ParseError> {
    aromatic_system.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn aromatic_system(i: &mut &str) -> PResult<AromaticSystemDsl> {
    multispace0.parse_next(i)?;
    let electrons = electron_counts(i)?;
    multispace0.parse_next(i)?;
    let preds: Vec<AromaticSystemPredicate> =
        repeat(0.., terminated(aromatic_system_predicate, multispace0)).parse_next(i)?;
    let mut form = AromaticSystemDsl(AromaticSystemAst::new(electrons));
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AromaticSystemPredicate {
    Charge(ValueAst),
    Spin(SpinPredicate),
    Electrons(ValueAst),
}

fn aromatic_system_predicate(i: &mut &str) -> PResult<AromaticSystemPredicate> {
    let start = *i;
    let prefix: &str = take(2usize).parse_next(i)?;
    match prefix {
        "#c" => charge.map(AromaticSystemPredicate::Charge).parse_next(i),
        "#u" => optional_value
            .map(|v| AromaticSystemPredicate::Spin(SpinPredicate::Unpaired(v)))
            .parse_next(i),
        "#s" => optional_value
            .map(|v| AromaticSystemPredicate::Spin(SpinPredicate::Multiplicity(v)))
            .parse_next(i),
        "#e" => optional_value
            .map(AromaticSystemPredicate::Electrons)
            .parse_next(i),
        p if p.starts_with('#') => Err(ErrMode::Cut(ParseError::UnknownAromaticSystemPredicate(
            p.to_string(),
        ))),
        _ => Err(ErrMode::Cut(ParseError::TrailingInput(start.to_string()))),
    }
}

fn apply_predicates(
    form: &mut AromaticSystemDsl,
    preds: Vec<AromaticSystemPredicate>,
) -> Result<(), ParseError> {
    let ast = &mut form.0;
    for pred in preds {
        match pred {
            AromaticSystemPredicate::Charge(v) => {
                if !matches!(ast.charge, ValueAst::Undetermined) {
                    return Err(ParseError::DuplicateAromaticSystemPredicate(
                        "#c".to_string(),
                    ));
                }
                ast.charge = v;
            }
            AromaticSystemPredicate::Spin(sp) => {
                apply_spin_pair(
                    &mut ast.spin,
                    sp,
                    ParseError::DuplicateAromaticSystemPredicate,
                )?;
            }
            AromaticSystemPredicate::Electrons(v) => {
                if has_electron_count(ast) {
                    return Err(ParseError::DuplicateAromaticSystemPredicate(
                        "#e".to_string(),
                    ));
                }
                ast.constraints
                    .add(AromaticSystemConstraint::ElectronCount(v));
            }
        }
    }
    Ok(())
}

fn has_electron_count(ast: &AromaticSystemAst) -> bool {
    ast.constraints
        .iter()
        .any(|c| matches!(c, AromaticSystemConstraint::ElectronCount(_)))
}

fn electron_count_value(ast: &AromaticSystemAst) -> Option<&ValueAst> {
    ast.constraints
        .iter()
        .map(|AromaticSystemConstraint::ElectronCount(v)| v)
        .next()
}

fn fmt_aromatic_system_ast(f: &mut fmt::Formatter<'_>, ast: &AromaticSystemAst) -> fmt::Result {
    fmt_electron_counts(f, &ast.electrons)?;
    fmt_charge(f, &ast.charge)?;
    fmt_spin_pair(f, &ast.spin)?;
    if let Some(v) = electron_count_value(ast) {
        fmt_electrons(f, v)?;
    }
    Ok(())
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

/// Partial aromatic system for a reaction `:modify` payload: reuses the aromatic parser but renders
/// an undetermined `#e` electron-count constraint explicitly (`#e*`) so its removal round-trips.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PartialAromaticSystemDsl(pub AromaticSystemAst);

impl FromStr for PartialAromaticSystemDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Self(parse_aromatic_system(s)?.0))
    }
}

impl Display for PartialAromaticSystemDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_electron_counts(f, &self.0.electrons)?;
        fmt_charge(f, &self.0.charge)?;
        fmt_spin_pair(f, &self.0.spin)?;
        match electron_count_value(&self.0) {
            Some(ValueAst::Undetermined) => write!(f, "#e*"),
            Some(v) => fmt_electrons(f, v),
            None => Ok(()),
        }
    }
}

impl<'de> FromEdn<'de> for PartialAromaticSystemDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s.parse().map_err(|e| DeError::subgrammar("aromatic-system", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for PartialAromaticSystemDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

fn raise_aromatic_system(ast: &mut AromaticSystemAst, cfg: &AromaticSystemDefaults) {
    let AromaticSystemAst {
        charge,
        spin,
        electrons: _,
        constraints: _,
    } = ast;

    if matches!(*charge, ValueAst::Undetermined) {
        *charge = match cfg.charge {
            NumericDefault::Zero => ValueAst::Lit(0),
            NumericDefault::Required => ValueAst::Undetermined,
        };
    }
    raise_spin(spin, cfg.unpaired_electrons, cfg.multiplicity);
}

fn lower_aromatic_system(ast: &mut AromaticSystemAst, cfg: &AromaticSystemDefaults) {
    let AromaticSystemAst {
        charge,
        spin,
        electrons: _,
        constraints: _,
    } = ast;

    if matches!(
        (&cfg.charge, &*charge),
        (NumericDefault::Zero, ValueAst::Lit(0))
    ) {
        *charge = ValueAst::Undetermined;
    }
    lower_spin(spin, cfg.unpaired_electrons, cfg.multiplicity);
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AromaticSystemConstraintDsl {
    ElectronCount(ValueAst),
}

impl AromaticSystemConstraintDsl {
    pub(crate) fn from_ast(c: &AromaticSystemConstraint) -> Self {
        match c {
            AromaticSystemConstraint::ElectronCount(v) => Self::ElectronCount(v.clone()),
        }
    }

    pub(crate) fn into_ast(self) -> AromaticSystemConstraint {
        match self {
            Self::ElectronCount(v) => AromaticSystemConstraint::ElectronCount(v),
        }
    }
}

impl<'de> FromEdn<'de> for AromaticSystemConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let map = match edn {
            Edn::Map(m) => m,
            other => {
                return Err(DeError::TypeMismatch {
                    expected: "map",
                    got: other.kind(),
                    path: vec!["aromatic-system-constraint".into()],
                });
            }
        };
        let mut entries = map.iter();
        let (key, value) = entries.next().ok_or_else(|| {
            DeError::Custom("expected single-key map for aromatic-system constraint".to_string())
        })?;
        if entries.next().is_some() {
            return Err(DeError::Custom(
                "aromatic-system constraint map has multiple keys".to_string(),
            ));
        }
        let kw = match key {
            Edn::Keyword(k) => k.name(),
            other => {
                return Err(DeError::TypeMismatch {
                    expected: "keyword",
                    got: other.kind(),
                    path: vec!["aromatic-system-constraint".into()],
                });
            }
        };
        match kw {
            "electron-count" => {
                let v = ValueDsl::from_edn(value)?;
                Ok(Self::ElectronCount(v.0))
            }
            other => Err(DeError::Custom(format!(
                "unknown aromatic-system constraint keyword :{}",
                other,
            ))),
        }
    }
}

impl ToEdn for AromaticSystemConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        match self {
            Self::ElectronCount(v) => {
                let value_edn = ValueDsl(v.clone()).to_edn();
                let mut map = EdnMap::with_capacity(1);
                map.insert(
                    Edn::Keyword(EdnKeyword::owned("electron-count".to_string())),
                    value_edn,
                );
                Edn::Map(map)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::read_string;

    use super::*;
    use crate::ast::constraint::AromaticSystemConstraints;
    use crate::ast::electrons::ElectronCountsAst;
    use crate::ast::spin::SpinStateAst;

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined("*", AromaticSystemDsl(AromaticSystemAst::default()))]
    #[case::whitespace("  [1,1,1]  #c+1  ", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: AromaticSystemConstraints::new() }))]
    #[case::charge_pos("[1,1,1]#c+1", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: AromaticSystemConstraints::new() }))]
    #[case::charge_neg("[1,1,1]#c-2", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(-2), spin: SpinStateAst::default(), constraints: AromaticSystemConstraints::new() }))]
    #[case::charge_plus_only("[1,1,1]#c+", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: AromaticSystemConstraints::new() }))]
    #[case::charge_minus_only("[1,1,1]#c-", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(-1), spin: SpinStateAst::default(), constraints: AromaticSystemConstraints::new() }))]
    #[case::electron_count("[1,1,1]#e6", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: AromaticSystemConstraints::from_iter([AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6))]) }))]
    #[case::electron_count_bare("[1,1,1]#e", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: AromaticSystemConstraints::from_iter([AromaticSystemConstraint::ElectronCount(ValueAst::Lit(1))]) }))]
    #[case::electron_count_wild("[1,1,1]#e*", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: AromaticSystemConstraints::from_iter([AromaticSystemConstraint::ElectronCount(ValueAst::Undetermined)]) }))]
    #[case::unpaired("[1,1,1]#u1", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Undetermined }, constraints: AromaticSystemConstraints::new() }))]
    #[case::mult("[1,1,1]#s2", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(2) }, constraints: AromaticSystemConstraints::new() }))]
    #[case::charge_electron_count("[1,1,1]#c+#e6", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: AromaticSystemConstraints::from_iter([AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6))]) }))]
    #[case::full("[1,1,1]#c0#u0#s1#e6", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)), constraints: AromaticSystemConstraints::from_iter([AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6))]) }))]
    fn test_parse_aromatic(#[case] input: &str, #[case] expected: AromaticSystemDsl) {
        let result = aromatic_system.parse(input);
        assert!(result.is_ok(), "{:?} should succeed, got {:?}", input, result.clone().unwrap_err());
        let form = result.unwrap();
        assert_eq!(form, expected);
    }

    #[rstest]
    #[case::missing_head("#c+1", ParseError::ExpectedElectronCounts)]
    #[case::unknown("[1,1,1]#x", ParseError::UnknownAromaticSystemPredicate("#x".to_string()))]
    #[case::unknown_a("[1,1,1]#a", ParseError::UnknownAromaticSystemPredicate("#a".to_string()))]
    #[case::dup_charge("[1,1,1]#c+#c-", ParseError::DuplicateAromaticSystemPredicate("#c".to_string()))]
    #[case::dup_electron_count("[1,1,1]#e6#e4", ParseError::DuplicateAromaticSystemPredicate("#e".to_string()))]
    #[case::dup_unpaired("[1,1,1]#u1#u2", ParseError::DuplicateAromaticSystemPredicate("#u".to_string()))]
    #[case::dup_multiplicity("[1,1,1]#s1#s2", ParseError::DuplicateAromaticSystemPredicate("#s".to_string()))]
    #[case::trailing("[1,1,1]#c+ foo", ParseError::TrailingInput("foo".to_string()))]
    fn test_parse_aromatic_error(#[case] input: &str, #[case] expected: ParseError) {
        let result = aromatic_system.parse(input);
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
    #[case::undetermined(AromaticSystemDsl::default(), "*")]
    #[case::charge_one(AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: AromaticSystemConstraints::new() }), "[1,1,1]#c+")]
    #[case::charge_neg_two(AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(-2), spin: SpinStateAst::default(), constraints: AromaticSystemConstraints::new() }), "[1,1,1]#c-2")]
    #[case::electron_count_six(AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: AromaticSystemConstraints::from_iter([AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6))]) }), "[1,1,1]#e6")]
    #[case::electron_count_one(AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: AromaticSystemConstraints::from_iter([AromaticSystemConstraint::ElectronCount(ValueAst::Lit(1))]) }), "[1,1,1]#e")]
    #[case::full(AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)), constraints: AromaticSystemConstraints::from_iter([AromaticSystemConstraint::ElectronCount(ValueAst::Lit(6))]) }), "[1,1,1]#c0#u0#s#e6")]
    fn test_display_aromatic(#[case] form: AromaticSystemDsl, #[case] expected: &str) {
        assert_eq!(form.to_string(), expected);
    }

    #[rstest]
    #[case::undetermined("*")]
    #[case::charge("[1,1,1]#c+1")]
    #[case::electron_count("[1,1,1]#e6")]
    #[case::unpaired("[1,1,1]#u2")]
    #[case::explicit_mult("[1,1,1]#s2")]
    fn test_aromatic_roundtrip(#[case] input: &str) {
        let form: AromaticSystemDsl = input.parse().unwrap();
        let rendered = form.to_string();
        let reparsed: AromaticSystemDsl = rendered.parse().unwrap();
        assert_eq!(form, reparsed);
    }

    #[rstest]
    fn test_aromatic_dsl_to_ast_fills_zero_defaults() {
        let dsl = AromaticSystemDsl::default();
        let cfg = AromaticSystemDefaults::zeroed();
        let ast = dsl.into_ast(&cfg);
        assert_eq!(ast.charge, ValueAst::Lit(0));
        assert_eq!(ast.spin, SpinStateAst::from((0_u8, 1_u8)));
        assert_eq!(ast.electrons, ElectronCountsAst::Undetermined);
        assert!(ast.constraints.is_empty());
    }

    #[rstest]
    fn test_aromatic_dsl_from_ast_strips_zero_defaults() {
        let ast = AromaticSystemAst {
            charge: ValueAst::Lit(0),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            electrons: ElectronCountsAst::Undetermined,
            constraints: AromaticSystemConstraints::new(),
        };
        let cfg = AromaticSystemDefaults::zeroed();
        let dsl = AromaticSystemDsl::from_ast(&ast, &cfg);
        assert_eq!(dsl.0.charge, ValueAst::Undetermined);
        assert_eq!(dsl.0.spin, SpinStateAst::default());
        assert_eq!(dsl.0.electrons, ElectronCountsAst::Undetermined);
        assert!(dsl.0.constraints.is_empty());
    }

    #[rstest]
    fn test_aromatic_dsl_roundtrip_zeroed() {
        let input = AromaticSystemDsl::default();
        let cfg = AromaticSystemDefaults::zeroed();
        let raised = input.clone().into_ast(&cfg);
        let lowered = AromaticSystemDsl::from_ast(&raised, &cfg);
        assert_eq!(input, lowered);
    }

    #[rstest]
    #[case::undetermined(r##""*""##)]
    #[case::charge(r##""[1,1,1]#c+""##)]
    #[case::full(r##""[1,1,1]#c0#u0#s1#e6""##)]
    fn test_aromatic_dsl_from_edn_str_matches_from_edn(#[case] input: &str) {
        let via_stream = AromaticSystemDsl::from_edn_str(input).unwrap();
        let tree = read_string(input).unwrap();
        let via_tree = AromaticSystemDsl::from_edn(&tree).unwrap();
        assert_eq!(via_stream, via_tree);
    }

    #[rstest]
    fn test_aromatic_system_constraint_dsl_from_edn_errors() {
        let edn = read_string("{:contains 1}").unwrap();
        let err = AromaticSystemConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    /// Vacuous aromatic-system inline constraint elides on rendering.
    /// `#e*` parses to `ElectronCount(Undetermined)` but the canonical
    /// surface form drops it.
    #[rstest]
    fn test_aromatic_render_elides_vacuous_electron_count() {
        let parsed: AromaticSystemDsl = aromatic_system.parse("[1,1,1]#e*").unwrap();
        assert_eq!(parsed.to_string(), "[1,1,1]");
        let reparsed: AromaticSystemDsl = aromatic_system.parse(&parsed.to_string()).unwrap();
        assert!(
            reparsed.0.constraints.is_empty(),
            "vacuous ElectronCount should be absent after render → reparse, got {:?}",
            reparsed.0.constraints,
        );
    }

    #[rstest]
    #[case::undetermined("*")]
    #[case::charged("[1,1,1]#c+")]
    fn test_aromatic_system_ast_from_str_to_string_roundtrip(#[case] s: &str) {
        let ast: AromaticSystemAst = s.parse().unwrap();
        assert_eq!(ast.to_string(), s);
    }
}
