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
    apply_unpaired_electrons_predicate, charge, fmt_charge, fmt_unpaired_electrons,
    lower_unpaired_electrons, optional_value, raise_unpaired_electrons, UnpairedElectronsPredicate,
};
use super::value::{fmt_value, ValueDsl};
use crate::ir::aromatic::{AromaticSystemForm, AromaticSystemUpdate};
use crate::ir::constraint::AromaticSystemConstraintForm;
use crate::ir::traits::{FromIr, IntoIr};
use crate::ir::value::NumForm;

/// Surface DSL wrapper around `AromaticSystemForm`.
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AromaticSystemDsl(pub AromaticSystemForm);

impl AromaticSystemDsl {
    /// Zero-cost reference cast from `&AromaticSystemForm`. Relies on `repr(transparent)`.
    pub fn from_ref(ast: &AromaticSystemForm) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(ast as *const AromaticSystemForm as *const Self) }
    }
}

impl From<AromaticSystemForm> for AromaticSystemDsl {
    fn from(ast: AromaticSystemForm) -> Self {
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
        fmt_aromatic_system_form(f, &self.0)
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

impl FromIr<AromaticSystemForm> for AromaticSystemDsl {
    type Ctx = AromaticSystemDefaults;

    fn from_ir(ast: &AromaticSystemForm, cfg: &Self::Ctx) -> Self {
        let mut out = ast.clone();
        lower_aromatic_system(&mut out, cfg);
        AromaticSystemDsl(out)
    }
}

impl IntoIr<AromaticSystemForm> for AromaticSystemDsl {
    type Ctx = AromaticSystemDefaults;

    fn into_ir(mut self, cfg: &Self::Ctx) -> AromaticSystemForm {
        raise_aromatic_system(&mut self.0, cfg);
        self.0
    }
}

impl FromStr for AromaticSystemForm {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(AromaticSystemDsl::from_str(s)?.into_ir(&AromaticSystemDefaults::default()))
    }
}

impl Display for AromaticSystemForm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        AromaticSystemDsl::from_ref(self).fmt(f)
    }
}

impl<'de> FromEdn<'de> for AromaticSystemForm {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        Ok(AromaticSystemDsl::from_edn(edn)?.into_ir(&AromaticSystemDefaults::default()))
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        Ok(AromaticSystemDsl::from_edn_str(input)?.into_ir(&AromaticSystemDefaults::default()))
    }
}

impl ToEdn for AromaticSystemForm {
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
    let mut form = AromaticSystemDsl(AromaticSystemForm::new(electrons));
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AromaticSystemPredicate {
    Charge(NumForm),
    UnpairedElectrons(UnpairedElectronsPredicate),
    Electrons(NumForm),
}

fn aromatic_system_predicate(i: &mut &str) -> PResult<AromaticSystemPredicate> {
    let start = *i;
    let prefix: &str = take(2usize).parse_next(i)?;
    match prefix {
        "#c" => charge.map(AromaticSystemPredicate::Charge).parse_next(i),
        "#u" => optional_value
            .map(|v| {
                AromaticSystemPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Count(v))
            })
            .parse_next(i),
        "#s" => optional_value
            .map(|v| {
                AromaticSystemPredicate::UnpairedElectrons(
                    UnpairedElectronsPredicate::Multiplicity(v),
                )
            })
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
                if !matches!(ast.charge, NumForm::Undetermined) {
                    return Err(ParseError::DuplicateAromaticSystemPredicate(
                        "#c".to_string(),
                    ));
                }
                ast.charge = v;
            }
            AromaticSystemPredicate::UnpairedElectrons(predicate) => {
                apply_unpaired_electrons_predicate(
                    &mut ast.unpaired_electrons,
                    predicate,
                    ParseError::DuplicateAromaticSystemPredicate,
                )?;
            }
            AromaticSystemPredicate::Electrons(v) => {
                let c = AromaticSystemConstraintForm::ElectronCount(v);
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

fn electron_count_value(ast: &AromaticSystemForm) -> Option<&NumForm> {
    ast.constraints
        .iter()
        .map(|AromaticSystemConstraintForm::ElectronCount(v)| v)
        .next()
}

fn fmt_aromatic_system_form(f: &mut fmt::Formatter<'_>, ast: &AromaticSystemForm) -> fmt::Result {
    fmt_electron_counts(f, &ast.electrons)?;
    fmt_charge(f, &ast.charge)?;
    fmt_unpaired_electrons(f, &ast.unpaired_electrons)?;
    if let Some(v) = electron_count_value(ast) {
        fmt_electrons(f, v)?;
    }
    Ok(())
}

fn fmt_electrons(f: &mut fmt::Formatter<'_>, v: &NumForm) -> fmt::Result {
    match v {
        NumForm::Undetermined => Ok(()),
        NumForm::Lit(1) => write!(f, "#e"),
        NumForm::Lit(n) => write!(f, "#e{}", n),
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

impl AromaticSystemUpdateDsl {
    /// Zero-cost reference cast from `&AromaticSystemUpdate`. Relies on `repr(transparent)`.
    pub fn from_ref(update: &AromaticSystemUpdate) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(update as *const AromaticSystemUpdate as *const Self) }
    }
}

impl FromIr<AromaticSystemUpdate> for AromaticSystemUpdateDsl {
    type Ctx = ();

    fn from_ir(update: &AromaticSystemUpdate, _ctx: &Self::Ctx) -> Self {
        Self(update.clone())
    }
}

impl IntoIr<AromaticSystemUpdate> for AromaticSystemUpdateDsl {
    type Ctx = ();

    fn into_ir(self, _ctx: &Self::Ctx) -> AromaticSystemUpdate {
        self.0
    }
}

impl FromStr for AromaticSystemUpdateDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_aromatic_system_update(s)
    }
}

impl FromStr for AromaticSystemUpdate {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(AromaticSystemUpdateDsl::from_str(s)?.into_ir(&()))
    }
}

impl Display for AromaticSystemUpdate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        AromaticSystemUpdateDsl::from_ref(self).fmt(f)
    }
}

impl Display for AromaticSystemUpdateDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(electrons) = &self.0.electrons {
            fmt_electron_counts(f, electrons)?;
        }
        if let Some(charge) = &self.0.charge {
            if matches!(charge, NumForm::Undetermined) {
                write!(f, "#c*")?;
            } else {
                fmt_charge(f, charge)?;
            }
        }
        if let Some(unpaired_electrons) = &self.0.unpaired_electrons.count {
            fmt_update_value_field(f, "#u", unpaired_electrons)?;
        }
        if let Some(multiplicity) = &self.0.unpaired_electrons.multiplicity {
            fmt_update_value_field(f, "#s", multiplicity)?;
        }
        for constraint in self.0.constraints.iter() {
            match constraint {
                AromaticSystemConstraintForm::ElectronCount(NumForm::Undetermined) => {
                    write!(f, "#e*")?;
                }
                AromaticSystemConstraintForm::ElectronCount(value) => fmt_electrons(f, value)?,
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
            AromaticSystemPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Count(
                value,
            )) => {
                if update.unpaired_electrons.count.replace(value).is_some() {
                    return Err(ParseError::DuplicateAromaticSystemPredicate(
                        "#u".to_string(),
                    ));
                }
            }
            AromaticSystemPredicate::UnpairedElectrons(
                UnpairedElectronsPredicate::Multiplicity(value),
            ) => {
                if update
                    .unpaired_electrons
                    .multiplicity
                    .replace(value)
                    .is_some()
                {
                    return Err(ParseError::DuplicateAromaticSystemPredicate(
                        "#s".to_string(),
                    ));
                }
            }
            AromaticSystemPredicate::Electrons(value) => {
                let constraint = AromaticSystemConstraintForm::ElectronCount(value);
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

fn fmt_update_value_field(f: &mut fmt::Formatter<'_>, prefix: &str, v: &NumForm) -> fmt::Result {
    match v {
        NumForm::Undetermined => write!(f, "{}*", prefix),
        NumForm::Lit(1) => write!(f, "{}", prefix),
        NumForm::Lit(n) => write!(f, "{}{}", prefix, n),
        value => {
            write!(f, "{}", prefix)?;
            fmt_value(f, value)
        }
    }
}

fn raise_aromatic_system(ast: &mut AromaticSystemForm, cfg: &AromaticSystemDefaults) {
    let AromaticSystemForm {
        charge,
        unpaired_electrons,
        electrons: _,
        constraints: _,
    } = ast;

    if matches!(*charge, NumForm::Undetermined) {
        *charge = match cfg.charge {
            NumericDefault::Zero => NumForm::Lit(0),
            NumericDefault::Required => NumForm::Undetermined,
        };
    }
    raise_unpaired_electrons(unpaired_electrons, cfg.unpaired_electrons, cfg.multiplicity);
}

fn lower_aromatic_system(ast: &mut AromaticSystemForm, cfg: &AromaticSystemDefaults) {
    let AromaticSystemForm {
        charge,
        unpaired_electrons,
        electrons: _,
        constraints: _,
    } = ast;

    if matches!(
        (&cfg.charge, &*charge),
        (NumericDefault::Zero, NumForm::Lit(0))
    ) {
        *charge = NumForm::Undetermined;
    }
    lower_unpaired_electrons(unpaired_electrons, cfg.unpaired_electrons, cfg.multiplicity);
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum AromaticSystemConstraintDsl {
    ElectronCount(NumForm),
}

impl AromaticSystemConstraintDsl {
    pub(crate) fn from_ir(c: &AromaticSystemConstraintForm) -> Self {
        match c {
            AromaticSystemConstraintForm::ElectronCount(v) => Self::ElectronCount(v.clone()),
        }
    }

    pub(crate) fn into_ir(self) -> AromaticSystemConstraintForm {
        match self {
            Self::ElectronCount(v) => AromaticSystemConstraintForm::ElectronCount(v),
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
    use crate::ir::constraint::AromaticSystemConstraintsForm;
    use crate::ir::electrons::ElectronCountsForm;
    use crate::ir::spin::{UnpairedElectronsForm, UnpairedElectronsUpdate};

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined("*", AromaticSystemDsl(AromaticSystemForm::default()))]
    #[case::whitespace("  [1,1,1]  #c+1  ", AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: AromaticSystemConstraintsForm::new() }))]
    #[case::charge_pos("[1,1,1]#c+1", AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: AromaticSystemConstraintsForm::new() }))]
    #[case::charge_neg("[1,1,1]#c-2", AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Lit(-2), unpaired_electrons: UnpairedElectronsForm::default(), constraints: AromaticSystemConstraintsForm::new() }))]
    #[case::charge_plus_only("[1,1,1]#c+", AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: AromaticSystemConstraintsForm::new() }))]
    #[case::charge_minus_only("[1,1,1]#c-", AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Lit(-1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: AromaticSystemConstraintsForm::new() }))]
    #[case::electron_count("[1,1,1]#e6", AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::ElectronCount(NumForm::Lit(6))) }))]
    #[case::electron_count_bare("[1,1,1]#e", AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::ElectronCount(NumForm::Lit(1))) }))]
    #[case::electron_count_wild("[1,1,1]#e*", AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::ElectronCount(NumForm::Undetermined)) }))]
    #[case::unpaired_electrons("[1,1,1]#u1", AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Undetermined }, constraints: AromaticSystemConstraintsForm::new() }))]
    #[case::mult("[1,1,1]#s2", AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(2) }, constraints: AromaticSystemConstraintsForm::new() }))]
    #[case::charge_electron_count("[1,1,1]#c+#e6", AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::ElectronCount(NumForm::Lit(6))) }))]
    #[case::full("[1,1,1]#c0#u0#s1#e6", AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::ElectronCount(NumForm::Lit(6))) }))]
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
    #[case::charge_one(AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: AromaticSystemConstraintsForm::new() }), "[1,1,1]#c+")]
    #[case::charge_neg_two(AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Lit(-2), unpaired_electrons: UnpairedElectronsForm::default(), constraints: AromaticSystemConstraintsForm::new() }), "[1,1,1]#c-2")]
    #[case::electron_count_six(AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::ElectronCount(NumForm::Lit(6))) }), "[1,1,1]#e6")]
    #[case::electron_count_one(AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::ElectronCount(NumForm::Lit(1))) }), "[1,1,1]#e")]
    #[case::full(AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::ElectronCount(NumForm::Lit(6))) }), "[1,1,1]#c0#u0#s#e6")]
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
    #[case::unpaired_electrons("[1,1,1]#u2")]
    #[case::explicit_mult("[1,1,1]#s2")]
    fn test_aromatic_system_dsl_display_roundtrip(#[case] input: &str) {
        let dsl = parse_aromatic_system(input).unwrap();
        assert_eq!(parse_aromatic_system(&dsl.to_string()).unwrap(), dsl);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count(
        AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(NumForm::Undetermined)) }),
        AromaticSystemDsl(AromaticSystemForm::from_electrons(vec![1, 1, 1])),
    )]
    fn test_aromatic_system_dsl_display_vacuous_constraints(
        #[case] input: AromaticSystemDsl,
        #[case] expected: AromaticSystemDsl,
    ) {
        assert_eq!(parse_aromatic_system(&input.to_string()).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(r##""*""##, AromaticSystemDsl(AromaticSystemForm::default()))]
    #[case::charge(r##""[1,1,1]#c+""##, AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: AromaticSystemConstraintsForm::new() }))]
    #[case::full(r##""[1,1,1]#c0#u0#s1#e6""##, AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6_i64)) }))]
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
    #[case::undetermined(AromaticSystemDsl(AromaticSystemForm::default()), r##""*""##)]
    #[case::charge(AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: AromaticSystemConstraintsForm::new() }), r##""[1,1,1]#c+""##)]
    #[case::full(AromaticSystemDsl(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(6_i64)) }), r##""[1,1,1]#c0#u0#s#e6""##)]
    fn test_aromatic_system_dsl_to_edn(
        #[case] input: AromaticSystemDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::zeroed(
        AromaticSystemForm { electrons: ElectronCountsForm::Undetermined, charge: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), constraints: AromaticSystemConstraintsForm::new() },
        AromaticSystemDsl(AromaticSystemForm::default()),
    )]
    fn test_aromatic_system_dsl_from_ast(
        #[case] input: AromaticSystemForm,
        #[case] expected: AromaticSystemDsl,
    ) {
        assert_eq!(
            AromaticSystemDsl::from_ir(&input, &AromaticSystemDefaults::zeroed()),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::zeroed(
        AromaticSystemDsl(AromaticSystemForm::default()),
        AromaticSystemForm { electrons: ElectronCountsForm::Undetermined, charge: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)), constraints: AromaticSystemConstraintsForm::new() },
    )]
    fn test_aromatic_system_dsl_into_ast(
        #[case] input: AromaticSystemDsl,
        #[case] expected: AromaticSystemForm,
    ) {
        assert_eq!(input.into_ir(&AromaticSystemDefaults::zeroed()), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined("*", AromaticSystemForm::default())]
    #[case::charged("[1,1,1]#c+", AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: AromaticSystemConstraintsForm::new() })]
    fn test_aromatic_system_form_from_str(
        #[case] input: &str,
        #[case] expected: AromaticSystemForm,
    ) {
        assert_eq!(input.parse::<AromaticSystemForm>().unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(AromaticSystemForm::default(), "*")]
    #[case::charged(AromaticSystemForm { electrons: ElectronCountsForm::Lit(vec![1, 1, 1]), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: AromaticSystemConstraintsForm::new() }, "[1,1,1]#c+")]
    fn test_aromatic_system_form_display(
        #[case] input: AromaticSystemForm,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", AromaticSystemUpdateDsl(AromaticSystemUpdate::default()))]
    #[case::electrons("[2,2,2]", AromaticSystemUpdateDsl(AromaticSystemUpdate { electrons: Some(ElectronCountsForm::Lit(vec![2, 2, 2])), ..Default::default() }))]
    #[case::electrons_undetermined("*", AromaticSystemUpdateDsl(AromaticSystemUpdate { electrons: Some(ElectronCountsForm::Undetermined), ..Default::default() }))]
    #[case::charge("#c-1", AromaticSystemUpdateDsl(AromaticSystemUpdate { charge: Some(NumForm::Lit(-1)), ..Default::default() }))]
    #[case::unpaired_electrons_unpaired("#u2", AromaticSystemUpdateDsl(AromaticSystemUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Lit(2)), multiplicity: None }, ..Default::default() }))]
    #[case::unpaired_electrons_multiplicity("#s1", AromaticSystemUpdateDsl(AromaticSystemUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) }, ..Default::default() }))]
    #[case::explicit_undetermined("*#c*#u*#s*#e*", AromaticSystemUpdateDsl(AromaticSystemUpdate { electrons: Some(ElectronCountsForm::Undetermined), charge: Some(NumForm::Undetermined), unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Undetermined), multiplicity: Some(NumForm::Undetermined) }, constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(NumForm::Undetermined)) }))]
    #[case::constraint_removal("#e*", AromaticSystemUpdateDsl(AromaticSystemUpdate { constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(NumForm::Undetermined)), ..Default::default() }))]
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

    #[rstest]
    #[case::duplicate_charge(
        "#c+#c-",
        ParseError::DuplicateAromaticSystemPredicate("#c".to_string())
    )]
    fn test_aromatic_system_update_from_str_error(
        #[case] input: &str,
        #[case] expected: ParseError,
    ) {
        assert_eq!(input.parse::<AromaticSystemUpdate>().unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(AromaticSystemUpdateDsl(AromaticSystemUpdate::default()), "")]
    #[case::fields(AromaticSystemUpdateDsl(AromaticSystemUpdate { electrons: Some(ElectronCountsForm::Lit(vec![2, 2, 2])), charge: Some(NumForm::Lit(-1)), unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) }, ..Default::default() }), "[2,2,2]#c-#s")]
    #[case::explicit_undetermined(AromaticSystemUpdateDsl(AromaticSystemUpdate { electrons: Some(ElectronCountsForm::Undetermined), charge: Some(NumForm::Undetermined), unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Undetermined), multiplicity: Some(NumForm::Undetermined) }, constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(NumForm::Undetermined)) }), "*#c*#u*#s*#e*")]
    fn test_aromatic_system_update_dsl_display(
        #[case] input: AromaticSystemUpdateDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fields(r##""[2,2,2]#c-#s""##, AromaticSystemUpdateDsl(AromaticSystemUpdate { electrons: Some(ElectronCountsForm::Lit(vec![2, 2, 2])), charge: Some(NumForm::Lit(-1)), unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) }, ..Default::default() }))]
    #[case::constraint_removal(r##""#e*""##, AromaticSystemUpdateDsl(AromaticSystemUpdate { constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(NumForm::Undetermined)), ..Default::default() }))]
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
    #[case::fields(AromaticSystemUpdateDsl(AromaticSystemUpdate { electrons: Some(ElectronCountsForm::Lit(vec![2, 2, 2])), charge: Some(NumForm::Lit(-1)), unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) }, ..Default::default() }), r##""[2,2,2]#c-#s""##)]
    #[case::explicit_undetermined(AromaticSystemUpdateDsl(AromaticSystemUpdate { electrons: Some(ElectronCountsForm::Undetermined), charge: Some(NumForm::Undetermined), unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Undetermined), multiplicity: Some(NumForm::Undetermined) }, constraints: AromaticSystemConstraintsForm::from(AromaticSystemConstraintForm::electron_count(NumForm::Undetermined)) }), r##""*#c*#u*#s*#e*""##)]
    fn test_aromatic_system_update_dsl_to_edn(
        #[case] input: AromaticSystemUpdateDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count(AromaticSystemConstraintForm::electron_count(6_i64), AromaticSystemConstraintDsl::ElectronCount(NumForm::Lit(6)))]
    fn test_aromatic_system_constraint_dsl_from_ast(
        #[case] input: AromaticSystemConstraintForm,
        #[case] expected: AromaticSystemConstraintDsl,
    ) {
        assert_eq!(AromaticSystemConstraintDsl::from_ir(&input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count(AromaticSystemConstraintDsl::ElectronCount(NumForm::Lit(6)), AromaticSystemConstraintForm::electron_count(6_i64))]
    fn test_aromatic_system_constraint_dsl_into_ast(
        #[case] input: AromaticSystemConstraintDsl,
        #[case] expected: AromaticSystemConstraintForm,
    ) {
        assert_eq!(input.into_ir(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count("{:electron-count 6}", AromaticSystemConstraintDsl::ElectronCount(NumForm::Lit(6)))]
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
    #[case::electron_count(AromaticSystemConstraintDsl::ElectronCount(NumForm::Lit(6)), "{:electron-count 6}")]
    fn test_aromatic_system_constraint_dsl_to_edn(
        #[case] input: AromaticSystemConstraintDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_edn(), read_string(expected).unwrap());
    }
}
