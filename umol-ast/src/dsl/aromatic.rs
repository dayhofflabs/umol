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
use super::predicate::{
    apply_spin_pair, charge, fmt_charge, fmt_spin_pair, lower_spin, optional_value, raise_spin,
    SpinPredicate,
};
use super::value::{fmt_value, ValueDsl};
use crate::ast::aromatic::{AromaticSystemAst, AromaticSystemUpdate};
use crate::ast::constraint::AromaticSystemConstraintAst;
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

/// Parse a complete aromatic-system update string.
pub fn parse_aromatic_system_update(input: &str) -> Result<AromaticSystemUpdateDsl, ParseError> {
    aromatic_system_update
        .parse(input)
        .map_err(|e| e.into_inner())
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
                let c = AromaticSystemConstraintAst::ElectronCount(v);
                if ast.constraints.contains(c.key()) {
                    return Err(ParseError::DuplicateAromaticSystemPredicate(
                        "#e".to_string(),
                    ));
                }
                ast.constraints.set(c);
            }
        }
    }
    Ok(())
}

fn electron_count_value(ast: &AromaticSystemAst) -> Option<&ValueAst> {
    ast.constraints
        .iter()
        .map(|AromaticSystemConstraintAst::ElectronCount(v)| v)
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

/// Surface DSL wrapper around an [`AromaticSystemUpdate`].
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AromaticSystemUpdateDsl(pub AromaticSystemUpdate);

impl FromStr for AromaticSystemUpdateDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_aromatic_system_update(s)
    }
}

impl Display for AromaticSystemUpdateDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(electrons) = &self.0.electrons {
            fmt_electron_counts(f, electrons)?;
        }
        if let Some(charge) = &self.0.charge {
            if matches!(charge, ValueAst::Undetermined) {
                write!(f, "#c*")?;
            } else {
                fmt_charge(f, charge)?;
            }
        }
        if let Some(unpaired) = &self.0.spin.unpaired {
            fmt_update_value_field(f, "#u", unpaired)?;
        }
        if let Some(multiplicity) = &self.0.spin.multiplicity {
            fmt_update_value_field(f, "#s", multiplicity)?;
        }
        for constraint in self.0.constraints.iter() {
            match constraint {
                AromaticSystemConstraintAst::ElectronCount(ValueAst::Undetermined) => {
                    write!(f, "#e*")?;
                }
                AromaticSystemConstraintAst::ElectronCount(value) => fmt_electrons(f, value)?,
            }
        }
        Ok(())
    }
}

impl<'de> FromEdn<'de> for AromaticSystemUpdateDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s
                .parse()
                .map_err(|e| DeError::subgrammar("aromatic-system-update", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for AromaticSystemUpdateDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

fn aromatic_system_update(i: &mut &str) -> PResult<AromaticSystemUpdateDsl> {
    multispace0.parse_next(i)?;
    let electrons = if i.starts_with('*') || i.starts_with('[') {
        Some(electron_counts(i)?)
    } else {
        None
    };
    multispace0.parse_next(i)?;
    let preds: Vec<AromaticSystemPredicate> =
        repeat(0.., terminated(aromatic_system_predicate, multispace0)).parse_next(i)?;
    let mut update = AromaticSystemUpdate {
        electrons,
        ..Default::default()
    };
    apply_update_predicates(&mut update, preds).map_err(ErrMode::Cut)?;
    Ok(AromaticSystemUpdateDsl(update))
}

fn apply_update_predicates(
    update: &mut AromaticSystemUpdate,
    preds: Vec<AromaticSystemPredicate>,
) -> Result<(), ParseError> {
    for pred in preds {
        match pred {
            AromaticSystemPredicate::Charge(value) => {
                if update.charge.replace(value).is_some() {
                    return Err(ParseError::DuplicateAromaticSystemPredicate(
                        "#c".to_string(),
                    ));
                }
            }
            AromaticSystemPredicate::Spin(SpinPredicate::Unpaired(value)) => {
                if update.spin.unpaired.replace(value).is_some() {
                    return Err(ParseError::DuplicateAromaticSystemPredicate(
                        "#u".to_string(),
                    ));
                }
            }
            AromaticSystemPredicate::Spin(SpinPredicate::Multiplicity(value)) => {
                if update.spin.multiplicity.replace(value).is_some() {
                    return Err(ParseError::DuplicateAromaticSystemPredicate(
                        "#s".to_string(),
                    ));
                }
            }
            AromaticSystemPredicate::Electrons(value) => {
                let constraint = AromaticSystemConstraintAst::ElectronCount(value);
                if update.constraints.contains(constraint.key()) {
                    return Err(ParseError::DuplicateAromaticSystemPredicate(
                        "#e".to_string(),
                    ));
                }
                update.constraints.set(constraint);
            }
        }
    }
    Ok(())
}

fn fmt_update_value_field(f: &mut fmt::Formatter<'_>, prefix: &str, v: &ValueAst) -> fmt::Result {
    match v {
        ValueAst::Undetermined => write!(f, "{}*", prefix),
        ValueAst::Lit(1) => write!(f, "{}", prefix),
        ValueAst::Lit(n) => write!(f, "{}{}", prefix, n),
        value => {
            write!(f, "{}", prefix)?;
            fmt_value(f, value)
        }
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
    pub(crate) fn from_ast(c: &AromaticSystemConstraintAst) -> Self {
        match c {
            AromaticSystemConstraintAst::ElectronCount(v) => Self::ElectronCount(v.clone()),
        }
    }

    pub(crate) fn into_ast(self) -> AromaticSystemConstraintAst {
        match self {
            Self::ElectronCount(v) => AromaticSystemConstraintAst::ElectronCount(v),
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
    use crate::ast::constraint::AromaticSystemConstraintsAst;
    use crate::ast::electrons::ElectronCountsAst;
    use crate::ast::spin::{SpinStateAst, SpinStateUpdate};

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined("*", AromaticSystemDsl(AromaticSystemAst::default()))]
    #[case::whitespace("  [1,1,1]  #c+1  ", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: AromaticSystemConstraintsAst::new() }))]
    #[case::charge_pos("[1,1,1]#c+1", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: AromaticSystemConstraintsAst::new() }))]
    #[case::charge_neg("[1,1,1]#c-2", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(-2), spin: SpinStateAst::default(), constraints: AromaticSystemConstraintsAst::new() }))]
    #[case::charge_plus_only("[1,1,1]#c+", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: AromaticSystemConstraintsAst::new() }))]
    #[case::charge_minus_only("[1,1,1]#c-", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(-1), spin: SpinStateAst::default(), constraints: AromaticSystemConstraintsAst::new() }))]
    #[case::electron_count("[1,1,1]#e6", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::ElectronCount(ValueAst::Lit(6))) }))]
    #[case::electron_count_bare("[1,1,1]#e", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::ElectronCount(ValueAst::Lit(1))) }))]
    #[case::electron_count_wild("[1,1,1]#e*", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::ElectronCount(ValueAst::Undetermined)) }))]
    #[case::unpaired("[1,1,1]#u1", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Undetermined }, constraints: AromaticSystemConstraintsAst::new() }))]
    #[case::mult("[1,1,1]#s2", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(2) }, constraints: AromaticSystemConstraintsAst::new() }))]
    #[case::charge_electron_count("[1,1,1]#c+#e6", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::ElectronCount(ValueAst::Lit(6))) }))]
    #[case::full("[1,1,1]#c0#u0#s1#e6", AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)), constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::ElectronCount(ValueAst::Lit(6))) }))]
    fn test_parse_aromatic_system(#[case] input: &str, #[case] expected: AromaticSystemDsl) {
        assert_eq!(parse_aromatic_system(input).unwrap(), expected);
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
    fn test_parse_aromatic_system_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(parse_aromatic_system(input).unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(AromaticSystemDsl::default(), "*")]
    #[case::charge_one(AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: AromaticSystemConstraintsAst::new() }), "[1,1,1]#c+")]
    #[case::charge_neg_two(AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(-2), spin: SpinStateAst::default(), constraints: AromaticSystemConstraintsAst::new() }), "[1,1,1]#c-2")]
    #[case::electron_count_six(AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::ElectronCount(ValueAst::Lit(6))) }), "[1,1,1]#e6")]
    #[case::electron_count_one(AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::ElectronCount(ValueAst::Lit(1))) }), "[1,1,1]#e")]
    #[case::full(AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)), constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::ElectronCount(ValueAst::Lit(6))) }), "[1,1,1]#c0#u0#s#e6")]
    fn test_aromatic_system_dsl_display(
        #[case] input: AromaticSystemDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_string(), expected);
    }

    #[rstest]
    #[case::undetermined("*")]
    #[case::charge("[1,1,1]#c+1")]
    #[case::electron_count("[1,1,1]#e6")]
    #[case::unpaired("[1,1,1]#u2")]
    #[case::explicit_mult("[1,1,1]#s2")]
    fn test_aromatic_system_dsl_display_roundtrip(#[case] input: &str) {
        let dsl = parse_aromatic_system(input).unwrap();
        assert_eq!(parse_aromatic_system(&dsl.to_string()).unwrap(), dsl);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count(
        AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(ValueAst::Undetermined)) }),
        AromaticSystemDsl(AromaticSystemAst::from_electrons(vec![1, 1, 1])),
    )]
    fn test_aromatic_system_dsl_display_vacuous_constraints(
        #[case] input: AromaticSystemDsl,
        #[case] expected: AromaticSystemDsl,
    ) {
        assert_eq!(parse_aromatic_system(&input.to_string()).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(r##""*""##, AromaticSystemDsl(AromaticSystemAst::default()))]
    #[case::charge(r##""[1,1,1]#c+""##, AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: AromaticSystemConstraintsAst::new() }))]
    #[case::full(r##""[1,1,1]#c0#u0#s1#e6""##, AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)), constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6_i64)) }))]
    fn test_aromatic_system_dsl_from_edn(
        #[case] input: &str,
        #[case] expected: AromaticSystemDsl,
    ) {
        assert_eq!(
            AromaticSystemDsl::from_edn(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rstest]
    #[case::wrong_type("1", DeError::TypeMismatch { expected: "string", got: "int", path: Vec::new() })]
    fn test_aromatic_system_dsl_from_edn_error(#[case] input: &str, #[case] expected: DeError) {
        assert_eq!(
            AromaticSystemDsl::from_edn(&read_string(input).unwrap()).unwrap_err(),
            expected,
        );
    }

    #[rstest]
    #[case::undetermined(r##""*""##)]
    #[case::charge(r##""[1,1,1]#c+""##)]
    #[case::full(r##""[1,1,1]#c0#u0#s1#e6""##)]
    fn test_aromatic_system_dsl_from_edn_parity(#[case] input: &str) {
        assert_eq!(
            AromaticSystemDsl::from_edn_str(input).unwrap(),
            AromaticSystemDsl::from_edn(&read_string(input).unwrap()).unwrap(),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(AromaticSystemDsl(AromaticSystemAst::default()), r##""*""##)]
    #[case::charge(AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: AromaticSystemConstraintsAst::new() }), r##""[1,1,1]#c+""##)]
    #[case::full(AromaticSystemDsl(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)), constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6_i64)) }), r##""[1,1,1]#c0#u0#s#e6""##)]
    fn test_aromatic_system_dsl_to_edn(
        #[case] input: AromaticSystemDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::zeroed(
        AromaticSystemAst { electrons: ElectronCountsAst::Undetermined, charge: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)), constraints: AromaticSystemConstraintsAst::new() },
        AromaticSystemDsl(AromaticSystemAst::default()),
    )]
    fn test_aromatic_system_dsl_from_ast(
        #[case] input: AromaticSystemAst,
        #[case] expected: AromaticSystemDsl,
    ) {
        assert_eq!(
            AromaticSystemDsl::from_ast(&input, &AromaticSystemDefaults::zeroed()),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::zeroed(
        AromaticSystemDsl(AromaticSystemAst::default()),
        AromaticSystemAst { electrons: ElectronCountsAst::Undetermined, charge: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)), constraints: AromaticSystemConstraintsAst::new() },
    )]
    fn test_aromatic_system_dsl_into_ast(
        #[case] input: AromaticSystemDsl,
        #[case] expected: AromaticSystemAst,
    ) {
        assert_eq!(input.into_ast(&AromaticSystemDefaults::zeroed()), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined("*", AromaticSystemAst::default())]
    #[case::charged("[1,1,1]#c+", AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: AromaticSystemConstraintsAst::new() })]
    fn test_aromatic_system_ast_from_str(
        #[case] input: &str,
        #[case] expected: AromaticSystemAst,
    ) {
        assert_eq!(input.parse::<AromaticSystemAst>().unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(AromaticSystemAst::default(), "*")]
    #[case::charged(AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: AromaticSystemConstraintsAst::new() }, "[1,1,1]#c+")]
    fn test_aromatic_system_ast_display(
        #[case] input: AromaticSystemAst,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", AromaticSystemUpdateDsl(AromaticSystemUpdate::default()))]
    #[case::electrons("[2,2,2]", AromaticSystemUpdateDsl(AromaticSystemUpdate { electrons: Some(ElectronCountsAst::Lit(vec![2, 2, 2])), ..Default::default() }))]
    #[case::electrons_undetermined("*", AromaticSystemUpdateDsl(AromaticSystemUpdate { electrons: Some(ElectronCountsAst::Undetermined), ..Default::default() }))]
    #[case::charge("#c-1", AromaticSystemUpdateDsl(AromaticSystemUpdate { charge: Some(ValueAst::Lit(-1)), ..Default::default() }))]
    #[case::spin_unpaired("#u2", AromaticSystemUpdateDsl(AromaticSystemUpdate { spin: SpinStateUpdate { unpaired: Some(ValueAst::Lit(2)), multiplicity: None }, ..Default::default() }))]
    #[case::spin_multiplicity("#s1", AromaticSystemUpdateDsl(AromaticSystemUpdate { spin: SpinStateUpdate { unpaired: None, multiplicity: Some(ValueAst::Lit(1)) }, ..Default::default() }))]
    #[case::explicit_undetermined("*#c*#u*#s*#e*", AromaticSystemUpdateDsl(AromaticSystemUpdate { electrons: Some(ElectronCountsAst::Undetermined), charge: Some(ValueAst::Undetermined), spin: SpinStateUpdate { unpaired: Some(ValueAst::Undetermined), multiplicity: Some(ValueAst::Undetermined) }, constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(ValueAst::Undetermined)) }))]
    #[case::constraint_removal("#e*", AromaticSystemUpdateDsl(AromaticSystemUpdate { constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(ValueAst::Undetermined)), ..Default::default() }))]
    fn test_parse_aromatic_system_update(
        #[case] input: &str,
        #[case] expected: AromaticSystemUpdateDsl,
    ) {
        assert_eq!(parse_aromatic_system_update(input).unwrap(), expected);
    }

    #[rstest]
    #[case::duplicate_charge("#c+#c-", ParseError::DuplicateAromaticSystemPredicate("#c".to_string()))]
    #[case::duplicate_unpaired("#u1#u2", ParseError::DuplicateAromaticSystemPredicate("#u".to_string()))]
    #[case::duplicate_multiplicity("#s1#s2", ParseError::DuplicateAromaticSystemPredicate("#s".to_string()))]
    #[case::duplicate_electron_count("#e6#e4", ParseError::DuplicateAromaticSystemPredicate("#e".to_string()))]
    fn test_parse_aromatic_system_update_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(parse_aromatic_system_update(input).unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(AromaticSystemUpdateDsl(AromaticSystemUpdate::default()), "")]
    #[case::fields(AromaticSystemUpdateDsl(AromaticSystemUpdate { electrons: Some(ElectronCountsAst::Lit(vec![2, 2, 2])), charge: Some(ValueAst::Lit(-1)), spin: SpinStateUpdate { unpaired: None, multiplicity: Some(ValueAst::Lit(1)) }, ..Default::default() }), "[2,2,2]#c-#s")]
    #[case::explicit_undetermined(AromaticSystemUpdateDsl(AromaticSystemUpdate { electrons: Some(ElectronCountsAst::Undetermined), charge: Some(ValueAst::Undetermined), spin: SpinStateUpdate { unpaired: Some(ValueAst::Undetermined), multiplicity: Some(ValueAst::Undetermined) }, constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(ValueAst::Undetermined)) }), "*#c*#u*#s*#e*")]
    fn test_aromatic_system_update_dsl_display(
        #[case] input: AromaticSystemUpdateDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fields(r##""[2,2,2]#c-#s""##, AromaticSystemUpdateDsl(AromaticSystemUpdate { electrons: Some(ElectronCountsAst::Lit(vec![2, 2, 2])), charge: Some(ValueAst::Lit(-1)), spin: SpinStateUpdate { unpaired: None, multiplicity: Some(ValueAst::Lit(1)) }, ..Default::default() }))]
    #[case::constraint_removal(r##""#e*""##, AromaticSystemUpdateDsl(AromaticSystemUpdate { constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(ValueAst::Undetermined)), ..Default::default() }))]
    fn test_aromatic_system_update_dsl_from_edn(
        #[case] input: &str,
        #[case] expected: AromaticSystemUpdateDsl,
    ) {
        assert_eq!(
            AromaticSystemUpdateDsl::from_edn(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rstest]
    #[case::wrong_type("1", DeError::TypeMismatch { expected: "string", got: "int", path: Vec::new() })]
    fn test_aromatic_system_update_dsl_from_edn_error(
        #[case] input: &str,
        #[case] expected: DeError,
    ) {
        assert_eq!(
            AromaticSystemUpdateDsl::from_edn(&read_string(input).unwrap()).unwrap_err(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(AromaticSystemUpdateDsl(AromaticSystemUpdate::default()), r##""""##)]
    #[case::fields(AromaticSystemUpdateDsl(AromaticSystemUpdate { electrons: Some(ElectronCountsAst::Lit(vec![2, 2, 2])), charge: Some(ValueAst::Lit(-1)), spin: SpinStateUpdate { unpaired: None, multiplicity: Some(ValueAst::Lit(1)) }, ..Default::default() }), r##""[2,2,2]#c-#s""##)]
    #[case::explicit_undetermined(AromaticSystemUpdateDsl(AromaticSystemUpdate { electrons: Some(ElectronCountsAst::Undetermined), charge: Some(ValueAst::Undetermined), spin: SpinStateUpdate { unpaired: Some(ValueAst::Undetermined), multiplicity: Some(ValueAst::Undetermined) }, constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(ValueAst::Undetermined)) }), r##""*#c*#u*#s*#e*""##)]
    fn test_aromatic_system_update_dsl_to_edn(
        #[case] input: AromaticSystemUpdateDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count(AromaticSystemConstraintAst::electron_count(6_i64), AromaticSystemConstraintDsl::ElectronCount(ValueAst::Lit(6)))]
    fn test_aromatic_system_constraint_dsl_from_ast(
        #[case] input: AromaticSystemConstraintAst,
        #[case] expected: AromaticSystemConstraintDsl,
    ) {
        assert_eq!(AromaticSystemConstraintDsl::from_ast(&input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count(AromaticSystemConstraintDsl::ElectronCount(ValueAst::Lit(6)), AromaticSystemConstraintAst::electron_count(6_i64))]
    fn test_aromatic_system_constraint_dsl_into_ast(
        #[case] input: AromaticSystemConstraintDsl,
        #[case] expected: AromaticSystemConstraintAst,
    ) {
        assert_eq!(input.into_ast(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count("{:electron-count 6}", AromaticSystemConstraintDsl::ElectronCount(ValueAst::Lit(6)))]
    fn test_aromatic_system_constraint_dsl_from_edn(
        #[case] input: &str,
        #[case] expected: AromaticSystemConstraintDsl,
    ) {
        assert_eq!(
            AromaticSystemConstraintDsl::from_edn(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rstest]
    #[case::unknown_key(
        "{:contains 1}",
        DeError::Custom("unknown aromatic-system constraint keyword :contains".to_string()),
    )]
    fn test_aromatic_system_constraint_dsl_from_edn_error(
        #[case] input: &str,
        #[case] expected: DeError,
    ) {
        assert_eq!(
            AromaticSystemConstraintDsl::from_edn(&read_string(input).unwrap()).unwrap_err(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count(AromaticSystemConstraintDsl::ElectronCount(ValueAst::Lit(6)), "{:electron-count 6}")]
    fn test_aromatic_system_constraint_dsl_to_edn(
        #[case] input: AromaticSystemConstraintDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_edn(), read_string(expected).unwrap());
    }
}
