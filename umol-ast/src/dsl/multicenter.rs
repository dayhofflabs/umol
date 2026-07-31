//! Multicenter-bond-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn};
use winnow::ascii::multispace0;
use winnow::combinator::{repeat, terminated};
use winnow::error::ErrMode;
use winnow::token::take;
use winnow::Parser;

use super::config::{MulticenterBondDefaults, NumericDefault};
use super::electrons::{electron_counts, fmt_electron_counts};
use super::error::{PResult, ParseError};
use super::predicate::{
    apply_unpaired_electrons_predicate, charge, fmt_charge, fmt_unpaired_electrons,
    lower_unpaired_electrons, optional_value, raise_unpaired_electrons, UnpairedElectronsPredicate,
};
use super::value::{fmt_value, ValueDsl};
use crate::ast::constraint::MulticenterBondConstraintAst;
use crate::ast::multicenter::{MulticenterBondAst, MulticenterBondUpdate};
use crate::ast::traits::{FromAst, IntoAst};
use crate::ast::value::ValueAst;

/// Surface DSL wrapper around `MulticenterBondAst`. The `electrons` field
/// (per-atom contributions) is serialized at the molecule level. The
/// `ElectronCount` constraint round-trips here as `#e<n>`.
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MulticenterBondDsl(pub MulticenterBondAst);

impl MulticenterBondDsl {
    /// Zero-cost reference cast from `&MulticenterBondAst`. Relies on `repr(transparent)`.
    pub fn from_ref(ast: &MulticenterBondAst) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(ast as *const MulticenterBondAst as *const Self) }
    }
}

impl From<MulticenterBondAst> for MulticenterBondDsl {
    fn from(ast: MulticenterBondAst) -> Self {
        Self(ast)
    }
}

impl FromStr for MulticenterBondDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_multicenter_bond(s)
    }
}

impl Display for MulticenterBondDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_multicenter_bond_ast(f, &self.0)
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
    type Ctx = MulticenterBondDefaults;

    fn from_ast(ast: &MulticenterBondAst, cfg: &Self::Ctx) -> Self {
        let mut out = ast.clone();
        lower_multicenter_bond(&mut out, cfg);
        MulticenterBondDsl(out)
    }
}

impl IntoAst<MulticenterBondAst> for MulticenterBondDsl {
    type Ctx = MulticenterBondDefaults;

    fn into_ast(mut self, cfg: &Self::Ctx) -> MulticenterBondAst {
        raise_multicenter_bond(&mut self.0, cfg);
        self.0
    }
}

impl FromStr for MulticenterBondAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(MulticenterBondDsl::from_str(s)?.into_ast(&MulticenterBondDefaults::default()))
    }
}

impl Display for MulticenterBondAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        MulticenterBondDsl::from_ref(self).fmt(f)
    }
}

impl<'de> FromEdn<'de> for MulticenterBondAst {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        Ok(MulticenterBondDsl::from_edn(edn)?.into_ast(&MulticenterBondDefaults::default()))
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        Ok(MulticenterBondDsl::from_edn_str(input)?.into_ast(&MulticenterBondDefaults::default()))
    }
}

impl ToEdn for MulticenterBondAst {
    fn to_edn(&self) -> Edn<'static> {
        MulticenterBondDsl::from_ref(self).to_edn()
    }
}

pub fn parse_multicenter_bond(input: &str) -> Result<MulticenterBondDsl, ParseError> {
    multicenter_bond.parse(input).map_err(|e| e.into_inner())
}

/// Parse a complete multicenter-bond update string.
pub fn parse_multicenter_bond_update(input: &str) -> Result<MulticenterBondUpdateDsl, ParseError> {
    multicenter_bond_update
        .parse(input)
        .map_err(|e| e.into_inner())
}

pub(crate) fn multicenter_bond(i: &mut &str) -> PResult<MulticenterBondDsl> {
    multispace0.parse_next(i)?;
    let electrons = electron_counts(i)?;
    multispace0.parse_next(i)?;
    let preds: Vec<MulticenterBondPredicate> =
        repeat(0.., terminated(multicenter_bond_predicate, multispace0)).parse_next(i)?;
    let mut form = MulticenterBondDsl(MulticenterBondAst::new(electrons));
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MulticenterBondPredicate {
    Charge(ValueAst),
    UnpairedElectrons(UnpairedElectronsPredicate),
    Electrons(ValueAst),
}

fn multicenter_bond_predicate(i: &mut &str) -> PResult<MulticenterBondPredicate> {
    let start = *i;
    let prefix: &str = take(2usize).parse_next(i)?;
    match prefix {
        "#c" => charge.map(MulticenterBondPredicate::Charge).parse_next(i),
        "#u" => optional_value
            .map(|v| {
                MulticenterBondPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Count(v))
            })
            .parse_next(i),
        "#s" => optional_value
            .map(|v| {
                MulticenterBondPredicate::UnpairedElectrons(
                    UnpairedElectronsPredicate::Multiplicity(v),
                )
            })
            .parse_next(i),
        "#e" => optional_value
            .map(MulticenterBondPredicate::Electrons)
            .parse_next(i),
        p if p.starts_with('#') => Err(ErrMode::Cut(ParseError::UnknownMulticenterBondPredicate(
            p.to_string(),
        ))),
        _ => Err(ErrMode::Cut(ParseError::TrailingInput(start.to_string()))),
    }
}

fn apply_predicates(
    form: &mut MulticenterBondDsl,
    preds: Vec<MulticenterBondPredicate>,
) -> Result<(), ParseError> {
    let ast = &mut form.0;
    for pred in preds {
        match pred {
            MulticenterBondPredicate::Charge(v) => {
                if !matches!(ast.charge, ValueAst::Undetermined) {
                    return Err(ParseError::DuplicateMulticenterBondPredicate(
                        "#c".to_string(),
                    ));
                }
                ast.charge = v;
            }
            MulticenterBondPredicate::UnpairedElectrons(predicate) => {
                apply_unpaired_electrons_predicate(
                    &mut ast.unpaired_electrons,
                    predicate,
                    ParseError::DuplicateMulticenterBondPredicate,
                )?;
            }
            MulticenterBondPredicate::Electrons(v) => {
                let c = MulticenterBondConstraintAst::ElectronCount(v);
                if ast.constraints.contains(c.key()) {
                    return Err(ParseError::DuplicateMulticenterBondPredicate(
                        "#e".to_string(),
                    ));
                }
                ast.constraints.set(c);
            }
        }
    }
    Ok(())
}

fn electron_count_value(ast: &MulticenterBondAst) -> Option<&ValueAst> {
    ast.constraints
        .iter()
        .map(|MulticenterBondConstraintAst::ElectronCount(v)| v)
        .next()
}

fn fmt_multicenter_bond_ast(f: &mut fmt::Formatter<'_>, ast: &MulticenterBondAst) -> fmt::Result {
    fmt_electron_counts(f, &ast.electrons)?;
    fmt_charge(f, &ast.charge)?;
    fmt_unpaired_electrons(f, &ast.unpaired_electrons)?;
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

/// Surface DSL wrapper around a [`MulticenterBondUpdate`].
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct MulticenterBondUpdateDsl(pub MulticenterBondUpdate);

impl FromStr for MulticenterBondUpdateDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_multicenter_bond_update(s)
    }
}

impl Display for MulticenterBondUpdateDsl {
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
        if let Some(unpaired) = &self.0.unpaired_electrons.count {
            fmt_update_value_field(f, "#u", unpaired)?;
        }
        if let Some(multiplicity) = &self.0.unpaired_electrons.multiplicity {
            fmt_update_value_field(f, "#s", multiplicity)?;
        }
        for constraint in self.0.constraints.iter() {
            match constraint {
                MulticenterBondConstraintAst::ElectronCount(ValueAst::Undetermined) => {
                    write!(f, "#e*")?;
                }
                MulticenterBondConstraintAst::ElectronCount(value) => fmt_electrons(f, value)?,
            }
        }
        Ok(())
    }
}

impl<'de> FromEdn<'de> for MulticenterBondUpdateDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s
                .parse()
                .map_err(|e| DeError::subgrammar("multicenter-bond-update", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for MulticenterBondUpdateDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

fn multicenter_bond_update(i: &mut &str) -> PResult<MulticenterBondUpdateDsl> {
    multispace0.parse_next(i)?;
    let electrons = if i.starts_with('*') || i.starts_with('[') {
        Some(electron_counts(i)?)
    } else {
        None
    };
    multispace0.parse_next(i)?;
    let preds: Vec<MulticenterBondPredicate> =
        repeat(0.., terminated(multicenter_bond_predicate, multispace0)).parse_next(i)?;
    let mut update = MulticenterBondUpdate {
        electrons,
        ..Default::default()
    };
    apply_update_predicates(&mut update, preds).map_err(ErrMode::Cut)?;
    Ok(MulticenterBondUpdateDsl(update))
}

fn apply_update_predicates(
    update: &mut MulticenterBondUpdate,
    preds: Vec<MulticenterBondPredicate>,
) -> Result<(), ParseError> {
    for pred in preds {
        match pred {
            MulticenterBondPredicate::Charge(value) => {
                if update.charge.replace(value).is_some() {
                    return Err(ParseError::DuplicateMulticenterBondPredicate(
                        "#c".to_string(),
                    ));
                }
            }
            MulticenterBondPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Count(
                value,
            )) => {
                if update.unpaired_electrons.count.replace(value).is_some() {
                    return Err(ParseError::DuplicateMulticenterBondPredicate(
                        "#u".to_string(),
                    ));
                }
            }
            MulticenterBondPredicate::UnpairedElectrons(
                UnpairedElectronsPredicate::Multiplicity(value),
            ) => {
                if update
                    .unpaired_electrons
                    .multiplicity
                    .replace(value)
                    .is_some()
                {
                    return Err(ParseError::DuplicateMulticenterBondPredicate(
                        "#s".to_string(),
                    ));
                }
            }
            MulticenterBondPredicate::Electrons(value) => {
                let constraint = MulticenterBondConstraintAst::ElectronCount(value);
                if update.constraints.contains(constraint.key()) {
                    return Err(ParseError::DuplicateMulticenterBondPredicate(
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

fn raise_multicenter_bond(ast: &mut MulticenterBondAst, cfg: &MulticenterBondDefaults) {
    let MulticenterBondAst {
        charge,
        unpaired_electrons,
        electrons: _,
        constraints: _,
    } = ast;

    if matches!(*charge, ValueAst::Undetermined) {
        *charge = match cfg.charge {
            NumericDefault::Zero => ValueAst::Lit(0),
            NumericDefault::Required => ValueAst::Undetermined,
        };
    }
    raise_unpaired_electrons(unpaired_electrons, cfg.unpaired_electrons, cfg.multiplicity);
}

fn lower_multicenter_bond(ast: &mut MulticenterBondAst, cfg: &MulticenterBondDefaults) {
    let MulticenterBondAst {
        charge,
        unpaired_electrons,
        electrons: _,
        constraints: _,
    } = ast;

    if matches!(
        (&cfg.charge, &*charge),
        (NumericDefault::Zero, ValueAst::Lit(0))
    ) {
        *charge = ValueAst::Undetermined;
    }
    lower_unpaired_electrons(unpaired_electrons, cfg.unpaired_electrons, cfg.multiplicity);
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum MulticenterBondConstraintDsl {
    ElectronCount(ValueAst),
}

impl MulticenterBondConstraintDsl {
    pub(crate) fn from_ast(c: &MulticenterBondConstraintAst) -> Self {
        match c {
            MulticenterBondConstraintAst::ElectronCount(v) => Self::ElectronCount(v.clone()),
        }
    }

    pub(crate) fn into_ast(self) -> MulticenterBondConstraintAst {
        match self {
            Self::ElectronCount(v) => MulticenterBondConstraintAst::ElectronCount(v),
        }
    }
}

impl<'de> FromEdn<'de> for MulticenterBondConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let map = match edn {
            Edn::Map(m) => m,
            other => {
                return Err(DeError::TypeMismatch {
                    expected: "map",
                    got: other.kind(),
                    path: vec!["multicenter-bond-constraint".into()],
                });
            }
        };
        let mut entries = map.iter();
        let (key, value) = entries.next().ok_or_else(|| {
            DeError::Custom("expected single-key map for multicenter-bond constraint".to_string())
        })?;
        if entries.next().is_some() {
            return Err(DeError::Custom(
                "multicenter-bond constraint map has multiple keys".to_string(),
            ));
        }
        let kw = match key {
            Edn::Keyword(k) => k.name(),
            other => {
                return Err(DeError::TypeMismatch {
                    expected: "keyword",
                    got: other.kind(),
                    path: vec!["multicenter-bond-constraint".into()],
                });
            }
        };
        match kw {
            "electron-count" => {
                let v = ValueDsl::from_edn(value)?;
                Ok(Self::ElectronCount(v.0))
            }
            other => Err(DeError::Custom(format!(
                "unknown multicenter-bond constraint keyword :{}",
                other,
            ))),
        }
    }
}

impl ToEdn for MulticenterBondConstraintDsl {
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
    use crate::ast::constraint::MulticenterBondConstraintsAst;
    use crate::ast::electrons::ElectronCountsAst;
    use crate::ast::spin::{UnpairedElectronsAst, UnpairedElectronsUpdate};

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined("*", MulticenterBondDsl(MulticenterBondAst::default()))]
    #[case::whitespace("  [1,1,1]  #c+1  ", MulticenterBondDsl(MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), unpaired_electrons: UnpairedElectronsAst::default(), constraints: MulticenterBondConstraintsAst::new() }))]
    #[case::charge_pos("[1,1,1]#c+1", MulticenterBondDsl(MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), unpaired_electrons: UnpairedElectronsAst::default(), constraints: MulticenterBondConstraintsAst::new() }))]
    #[case::charge_neg("[1,1,1]#c-2", MulticenterBondDsl(MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(-2), unpaired_electrons: UnpairedElectronsAst::default(), constraints: MulticenterBondConstraintsAst::new() }))]
    #[case::electron_count("[1,1,1]#e6", MulticenterBondDsl(MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, unpaired_electrons: UnpairedElectronsAst::default(), constraints: MulticenterBondConstraintsAst::from_iter([MulticenterBondConstraintAst::ElectronCount(ValueAst::Lit(6))]) }))]
    #[case::electron_count_bare("[1,1,1]#e", MulticenterBondDsl(MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, unpaired_electrons: UnpairedElectronsAst::default(), constraints: MulticenterBondConstraintsAst::from_iter([MulticenterBondConstraintAst::ElectronCount(ValueAst::Lit(1))]) }))]
    #[case::unpaired("[1,1,1]#u1", MulticenterBondDsl(MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, unpaired_electrons: UnpairedElectronsAst { count: ValueAst::Lit(1), multiplicity: ValueAst::Undetermined }, constraints: MulticenterBondConstraintsAst::new() }))]
    #[case::mult("[1,1,1]#s2", MulticenterBondDsl(MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, unpaired_electrons: UnpairedElectronsAst { count: ValueAst::Undetermined, multiplicity: ValueAst::Lit(2) }, constraints: MulticenterBondConstraintsAst::new() }))]
    #[case::charge_electron_count("[1,1,1]#c+#e2", MulticenterBondDsl(MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), unpaired_electrons: UnpairedElectronsAst::default(), constraints: MulticenterBondConstraintsAst::from_iter([MulticenterBondConstraintAst::ElectronCount(ValueAst::Lit(2))]) }))]
    #[case::full("[1,1,1]#c0#u0#s1#e2", MulticenterBondDsl(MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(0), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)), constraints: MulticenterBondConstraintsAst::from_iter([MulticenterBondConstraintAst::ElectronCount(ValueAst::Lit(2))]) }))]
    fn test_parse_multicenter_bond(#[case] input: &str, #[case] expected: MulticenterBondDsl) {
        assert_eq!(parse_multicenter_bond(input).unwrap(), expected);
    }

    #[rstest]
    #[case::missing_head("#c+1", ParseError::ExpectedElectronCounts)]
    #[case::unknown("[1,1,1]#x", ParseError::UnknownMulticenterBondPredicate("#x".to_string()))]
    #[case::unknown_a("[1,1,1]#a", ParseError::UnknownMulticenterBondPredicate("#a".to_string()))]
    #[case::dup_charge("[1,1,1]#c+#c-", ParseError::DuplicateMulticenterBondPredicate("#c".to_string()))]
    #[case::dup_electron_count("[1,1,1]#e2#e4", ParseError::DuplicateMulticenterBondPredicate("#e".to_string()))]
    #[case::dup_unpaired("[1,1,1]#u1#u2", ParseError::DuplicateMulticenterBondPredicate("#u".to_string()))]
    #[case::dup_multiplicity("[1,1,1]#s1#s2", ParseError::DuplicateMulticenterBondPredicate("#s".to_string()))]
    #[case::trailing("[1,1,1]#c+ foo", ParseError::TrailingInput("foo".to_string()))]
    fn test_parse_multicenter_bond_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(parse_multicenter_bond(input).unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(MulticenterBondDsl::default(), "*")]
    #[case::charge_one(MulticenterBondDsl(MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), unpaired_electrons: UnpairedElectronsAst::default(), constraints: MulticenterBondConstraintsAst::new() }), "[1,1,1]#c+")]
    #[case::full(MulticenterBondDsl(MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(0), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)), constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(2_i64)) }), "[1,1,1]#c0#u0#s#e2")]
    fn test_multicenter_bond_dsl_display(
        #[case] input: MulticenterBondDsl,
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
    fn test_multicenter_bond_dsl_display_roundtrip(#[case] input: &str) {
        let dsl = parse_multicenter_bond(input).unwrap();
        assert_eq!(parse_multicenter_bond(&dsl.to_string()).unwrap(), dsl);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count(
        MulticenterBondDsl(MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Undetermined, unpaired_electrons: UnpairedElectronsAst::default(), constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(ValueAst::Undetermined)) }),
        MulticenterBondDsl(MulticenterBondAst::from_electrons(vec![1, 1, 1])),
    )]
    fn test_multicenter_bond_dsl_display_vacuous_constraints(
        #[case] input: MulticenterBondDsl,
        #[case] expected: MulticenterBondDsl,
    ) {
        assert_eq!(parse_multicenter_bond(&input.to_string()).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(r##""*""##, MulticenterBondDsl(MulticenterBondAst::default()))]
    #[case::full(r##""[1,1,1]#c0#u0#s1#e2""##, MulticenterBondDsl(MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(0), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)), constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(2_i64)) }))]
    fn test_multicenter_bond_dsl_from_edn(
        #[case] input: &str,
        #[case] expected: MulticenterBondDsl,
    ) {
        assert_eq!(
            MulticenterBondDsl::from_edn(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rstest]
    #[case::wrong_type("1", DeError::TypeMismatch { expected: "string", got: "int", path: Vec::new() })]
    fn test_multicenter_bond_dsl_from_edn_error(#[case] input: &str, #[case] expected: DeError) {
        assert_eq!(
            MulticenterBondDsl::from_edn(&read_string(input).unwrap()).unwrap_err(),
            expected,
        );
    }

    #[rstest]
    #[case::undetermined(r##""*""##)]
    #[case::charge(r##""[1,1,1]#c+""##)]
    #[case::full(r##""[1,1,1]#c0#u0#s1#e2""##)]
    fn test_multicenter_bond_dsl_from_edn_parity(#[case] input: &str) {
        assert_eq!(
            MulticenterBondDsl::from_edn_str(input).unwrap(),
            MulticenterBondDsl::from_edn(&read_string(input).unwrap()).unwrap(),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(MulticenterBondDsl(MulticenterBondAst::default()), r##""*""##)]
    #[case::full(MulticenterBondDsl(MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(0), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)), constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(2_i64)) }), r##""[1,1,1]#c0#u0#s#e2""##)]
    fn test_multicenter_bond_dsl_to_edn(
        #[case] input: MulticenterBondDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::zeroed(
        MulticenterBondAst { electrons: ElectronCountsAst::Undetermined, charge: ValueAst::Lit(0), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)), constraints: MulticenterBondConstraintsAst::new() },
        MulticenterBondDsl(MulticenterBondAst::default()),
    )]
    fn test_multicenter_bond_dsl_from_ast(
        #[case] input: MulticenterBondAst,
        #[case] expected: MulticenterBondDsl,
    ) {
        assert_eq!(
            MulticenterBondDsl::from_ast(&input, &MulticenterBondDefaults::zeroed()),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::zeroed(
        MulticenterBondDsl(MulticenterBondAst::default()),
        MulticenterBondAst { electrons: ElectronCountsAst::Undetermined, charge: ValueAst::Lit(0), unpaired_electrons: UnpairedElectronsAst::from((0_u8, 1_u8)), constraints: MulticenterBondConstraintsAst::new() },
    )]
    fn test_multicenter_bond_dsl_into_ast(
        #[case] input: MulticenterBondDsl,
        #[case] expected: MulticenterBondAst,
    ) {
        assert_eq!(input.into_ast(&MulticenterBondDefaults::zeroed()), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined("*", MulticenterBondAst::default())]
    #[case::charged("[1,1,1]#c+", MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), unpaired_electrons: UnpairedElectronsAst::default(), constraints: MulticenterBondConstraintsAst::new() })]
    fn test_multicenter_bond_ast_from_str(
        #[case] input: &str,
        #[case] expected: MulticenterBondAst,
    ) {
        assert_eq!(input.parse::<MulticenterBondAst>().unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(MulticenterBondAst::default(), "*")]
    #[case::charged(MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]), charge: ValueAst::Lit(1), unpaired_electrons: UnpairedElectronsAst::default(), constraints: MulticenterBondConstraintsAst::new() }, "[1,1,1]#c+")]
    fn test_multicenter_bond_ast_display(
        #[case] input: MulticenterBondAst,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", MulticenterBondUpdateDsl(MulticenterBondUpdate::default()))]
    #[case::electrons("[2,2,2]", MulticenterBondUpdateDsl(MulticenterBondUpdate { electrons: Some(ElectronCountsAst::Lit(vec![2, 2, 2])), ..Default::default() }))]
    #[case::electrons_undetermined("*", MulticenterBondUpdateDsl(MulticenterBondUpdate { electrons: Some(ElectronCountsAst::Undetermined), ..Default::default() }))]
    #[case::charge("#c-1", MulticenterBondUpdateDsl(MulticenterBondUpdate { charge: Some(ValueAst::Lit(-1)), ..Default::default() }))]
    #[case::unpaired_electrons_unpaired("#u2", MulticenterBondUpdateDsl(MulticenterBondUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: Some(ValueAst::Lit(2)), multiplicity: None }, ..Default::default() }))]
    #[case::unpaired_electrons_multiplicity("#s1", MulticenterBondUpdateDsl(MulticenterBondUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(ValueAst::Lit(1)) }, ..Default::default() }))]
    #[case::explicit_undetermined("*#c*#u*#s*#e*", MulticenterBondUpdateDsl(MulticenterBondUpdate { electrons: Some(ElectronCountsAst::Undetermined), charge: Some(ValueAst::Undetermined), unpaired_electrons: UnpairedElectronsUpdate { count: Some(ValueAst::Undetermined), multiplicity: Some(ValueAst::Undetermined) }, constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(ValueAst::Undetermined)) }))]
    #[case::constraint_removal("#e*", MulticenterBondUpdateDsl(MulticenterBondUpdate { constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(ValueAst::Undetermined)), ..Default::default() }))]
    fn test_parse_multicenter_bond_update(
        #[case] input: &str,
        #[case] expected: MulticenterBondUpdateDsl,
    ) {
        assert_eq!(parse_multicenter_bond_update(input).unwrap(), expected);
    }

    #[rstest]
    #[case::duplicate_charge("#c+#c-", ParseError::DuplicateMulticenterBondPredicate("#c".to_string()))]
    #[case::duplicate_unpaired("#u1#u2", ParseError::DuplicateMulticenterBondPredicate("#u".to_string()))]
    #[case::duplicate_multiplicity("#s1#s2", ParseError::DuplicateMulticenterBondPredicate("#s".to_string()))]
    #[case::duplicate_electron_count("#e6#e4", ParseError::DuplicateMulticenterBondPredicate("#e".to_string()))]
    fn test_parse_multicenter_bond_update_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(parse_multicenter_bond_update(input).unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(MulticenterBondUpdateDsl(MulticenterBondUpdate::default()), "")]
    #[case::fields(MulticenterBondUpdateDsl(MulticenterBondUpdate { electrons: Some(ElectronCountsAst::Lit(vec![2, 2, 2])), charge: Some(ValueAst::Lit(-1)), unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(ValueAst::Lit(1)) }, ..Default::default() }), "[2,2,2]#c-#s")]
    #[case::explicit_undetermined(MulticenterBondUpdateDsl(MulticenterBondUpdate { electrons: Some(ElectronCountsAst::Undetermined), charge: Some(ValueAst::Undetermined), unpaired_electrons: UnpairedElectronsUpdate { count: Some(ValueAst::Undetermined), multiplicity: Some(ValueAst::Undetermined) }, constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(ValueAst::Undetermined)) }), "*#c*#u*#s*#e*")]
    fn test_multicenter_bond_update_dsl_display(
        #[case] input: MulticenterBondUpdateDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fields(r##""[2,2,2]#c-#s""##, MulticenterBondUpdateDsl(MulticenterBondUpdate { electrons: Some(ElectronCountsAst::Lit(vec![2, 2, 2])), charge: Some(ValueAst::Lit(-1)), unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(ValueAst::Lit(1)) }, ..Default::default() }))]
    #[case::constraint_removal(r##""#e*""##, MulticenterBondUpdateDsl(MulticenterBondUpdate { constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(ValueAst::Undetermined)), ..Default::default() }))]
    fn test_multicenter_bond_update_dsl_from_edn(
        #[case] input: &str,
        #[case] expected: MulticenterBondUpdateDsl,
    ) {
        assert_eq!(
            MulticenterBondUpdateDsl::from_edn(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rstest]
    #[case::wrong_type("1", DeError::TypeMismatch { expected: "string", got: "int", path: Vec::new() })]
    fn test_multicenter_bond_update_dsl_from_edn_error(
        #[case] input: &str,
        #[case] expected: DeError,
    ) {
        assert_eq!(
            MulticenterBondUpdateDsl::from_edn(&read_string(input).unwrap()).unwrap_err(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(MulticenterBondUpdateDsl(MulticenterBondUpdate::default()), r##""""##)]
    #[case::fields(MulticenterBondUpdateDsl(MulticenterBondUpdate { electrons: Some(ElectronCountsAst::Lit(vec![2, 2, 2])), charge: Some(ValueAst::Lit(-1)), unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(ValueAst::Lit(1)) }, ..Default::default() }), r##""[2,2,2]#c-#s""##)]
    #[case::explicit_undetermined(MulticenterBondUpdateDsl(MulticenterBondUpdate { electrons: Some(ElectronCountsAst::Undetermined), charge: Some(ValueAst::Undetermined), unpaired_electrons: UnpairedElectronsUpdate { count: Some(ValueAst::Undetermined), multiplicity: Some(ValueAst::Undetermined) }, constraints: MulticenterBondConstraintsAst::from(MulticenterBondConstraintAst::electron_count(ValueAst::Undetermined)) }), r##""*#c*#u*#s*#e*""##)]
    fn test_multicenter_bond_update_dsl_to_edn(
        #[case] input: MulticenterBondUpdateDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count(MulticenterBondConstraintAst::electron_count(6_i64), MulticenterBondConstraintDsl::ElectronCount(ValueAst::Lit(6)))]
    fn test_multicenter_bond_constraint_dsl_from_ast(
        #[case] input: MulticenterBondConstraintAst,
        #[case] expected: MulticenterBondConstraintDsl,
    ) {
        assert_eq!(MulticenterBondConstraintDsl::from_ast(&input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count(MulticenterBondConstraintDsl::ElectronCount(ValueAst::Lit(6)), MulticenterBondConstraintAst::electron_count(6_i64))]
    fn test_multicenter_bond_constraint_dsl_into_ast(
        #[case] input: MulticenterBondConstraintDsl,
        #[case] expected: MulticenterBondConstraintAst,
    ) {
        assert_eq!(input.into_ast(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count("{:electron-count 6}", MulticenterBondConstraintDsl::ElectronCount(ValueAst::Lit(6)))]
    fn test_multicenter_bond_constraint_dsl_from_edn(
        #[case] input: &str,
        #[case] expected: MulticenterBondConstraintDsl,
    ) {
        assert_eq!(
            MulticenterBondConstraintDsl::from_edn(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rstest]
    #[case::unknown_key(
        "{:contains 1}",
        DeError::Custom("unknown multicenter-bond constraint keyword :contains".to_string()),
    )]
    fn test_multicenter_bond_constraint_dsl_from_edn_error(
        #[case] input: &str,
        #[case] expected: DeError,
    ) {
        assert_eq!(
            MulticenterBondConstraintDsl::from_edn(&read_string(input).unwrap()).unwrap_err(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electron_count(MulticenterBondConstraintDsl::ElectronCount(ValueAst::Lit(6)), "{:electron-count 6}")]
    fn test_multicenter_bond_constraint_dsl_to_edn(
        #[case] input: MulticenterBondConstraintDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_edn(), read_string(expected).unwrap());
    }
}
