//! Bond-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnStreamDeserializer, FromEdn, ToEdn};
use winnow::ascii::multispace0;
use winnow::combinator::{delimited, opt, repeat, terminated};
use winnow::error::ErrMode;
use winnow::token::take;
use winnow::Parser;

use super::boolean::{boolean, BooleanDsl};
use super::config::{BondDefaults, NumericDefault, StereoDefault};
use super::constraint::RingMembershipDsl;
use super::edn_utils::single_key_map;
use super::error::{PResult, ParseError};
use super::predicate::{
    apply_unpaired_electrons_predicate, charge, fmt_charge, fmt_ring_membership,
    fmt_unpaired_electrons, lower_unpaired_electrons, optional_value, raise_unpaired_electrons,
    ring_membership, UnpairedElectronsPredicate,
};
use super::stereo::{cis_trans_stereo_config, fmt_cis_trans_stereo_config, CisTransStereoDsl};
use super::value::{fmt_value, value};
use crate::ir::bond::{BondForm, BondUpdate};
use crate::ir::boolean::BooleanForm;
use crate::ir::constraint::{BondConstraintAst, BondConstraintKey, BondConstraintsAst, RingScope};
use crate::ir::spin::UnpairedElectronsForm;
use crate::ir::stereo::CisTransStereoForm;
use crate::ir::traits::{FromIr, IntoIr, Lattice};
use crate::ir::value::NumForm;

/// Surface DSL wrapper around `BondForm`.
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BondDsl(pub BondForm);

impl BondDsl {
    /// Zero-cost reference cast from `&BondForm`. Relies on `repr(transparent)`.
    pub fn from_ref(ast: &BondForm) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(ast as *const BondForm as *const Self) }
    }
}

impl From<BondForm> for BondDsl {
    fn from(ast: BondForm) -> Self {
        Self(ast)
    }
}

impl FromStr for BondDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_bond(s)
    }
}

impl Display for BondDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_bond_form(f, &self.0)?;
        for c in self.0.constraints.iter() {
            fmt_constraint(f, c)?;
        }
        Ok(())
    }
}

impl<'de> FromEdn<'de> for BondDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s.parse().map_err(|e| DeError::subgrammar("bond", e)),
            Edn::Keyword(k) => {
                let s = expand_bond_keyword(k.name()).ok_or_else(|| {
                    DeError::Custom(format!("unknown bond keyword :{}", k.name()))
                })?;
                s.parse().map_err(|e| DeError::subgrammar("bond", e))
            }
            other => Err(DeError::TypeMismatch {
                expected: "string or bond-keyword",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        EdnStreamDeserializer::new(input).read_subgrammar_all("bond")
    }
}

/// Expand a bond-entry keyword shorthand to its equivalent bond-string
/// payload. The five recognized keywords mirror the spec §7.6 table:
///
/// - `:single` → `"1"`
/// - `:double` → `"2"`
/// - `:triple` → `"3"`
/// - `:quadruple` → `"4"`
/// - `:aromatic` → `"1#a"` (single-order localized bond participating in
///   an aromatic system; `#a` is the bond-string aromatic predicate)
///
/// Returns `None` for unrecognized keywords.
pub(crate) fn expand_bond_keyword(name: &str) -> Option<&'static str> {
    match name {
        "single" => Some("1"),
        "double" => Some("2"),
        "triple" => Some("3"),
        "quadruple" => Some("4"),
        "aromatic" => Some("1#a"),
        _ => None,
    }
}

impl ToEdn for BondDsl {
    fn to_edn(&self) -> Edn<'static> {
        match bond_keyword_for(&self.0) {
            Some(kw) => Edn::Keyword(EdnKeyword::owned(kw.to_string())),
            None => Edn::Str(Cow::Owned(self.to_string())),
        }
    }
}

/// Return the bond-keyword shorthand for canonical bond shapes, or `None`
/// when the bond requires the full bond-string form. Inverse of
/// [`expand_bond_keyword`]: every shape this returns must round-trip.
///
/// Canonical means: charge/unpaired electrons at their defaults (Undetermined / default
/// pair) and either no constraints (orders 1–4) or exactly the `Aromatic`
/// flag (order 1 → `:aromatic`).
fn bond_keyword_for(ast: &BondForm) -> Option<&'static str> {
    if !matches!(ast.charge, NumForm::Undetermined)
        || ast.unpaired_electrons != UnpairedElectronsForm::default()
    {
        return None;
    }
    let constraints: Vec<&BondConstraintAst> = ast.constraints.iter().collect();
    match (&ast.order, constraints.as_slice()) {
        (NumForm::Lit(1), []) => Some("single"),
        (NumForm::Lit(2), []) => Some("double"),
        (NumForm::Lit(3), []) => Some("triple"),
        (NumForm::Lit(4), []) => Some("quadruple"),
        (NumForm::Lit(1), [BondConstraintAst::Aromatic(BooleanForm::Lit(true))]) => {
            Some("aromatic")
        }
        _ => None,
    }
}

impl FromIr<BondForm> for BondDsl {
    type Ctx = BondDefaults;

    fn from_ir(ast: &BondForm, cfg: &Self::Ctx) -> Self {
        let mut out = ast.clone();
        lower_bond(&mut out, cfg);
        BondDsl(out)
    }
}

impl IntoIr<BondForm> for BondDsl {
    type Ctx = BondDefaults;

    fn into_ir(mut self, cfg: &Self::Ctx) -> BondForm {
        raise_bond(&mut self.0, cfg);
        self.0
    }
}

impl FromStr for BondForm {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(BondDsl::from_str(s)?.into_ir(&BondDefaults::default()))
    }
}

impl Display for BondForm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        BondDsl::from_ref(self).fmt(f)
    }
}

impl<'de> FromEdn<'de> for BondForm {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        Ok(BondDsl::from_edn(edn)?.into_ir(&BondDefaults::default()))
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        Ok(BondDsl::from_edn_str(input)?.into_ir(&BondDefaults::default()))
    }
}

impl ToEdn for BondForm {
    fn to_edn(&self) -> Edn<'static> {
        BondDsl::from_ref(self).to_edn()
    }
}

/// Parse bond string into a `BondDsl`.
pub fn parse_bond(input: &str) -> Result<BondDsl, ParseError> {
    bond.parse(input).map_err(|e| e.into_inner())
}

/// Bond-string parser (does not require consuming all input).
fn bond(i: &mut &str) -> PResult<BondDsl> {
    let order = delimited(multispace0, value, multispace0).parse_next(i)?;
    let preds: Vec<BondPredicate> =
        repeat(0.., terminated(bond_predicate, multispace0)).parse_next(i)?;
    let mut form = BondDsl(BondForm::new(order));
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

/// Surface DSL wrapper around a [`BondUpdate`].
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BondUpdateDsl(pub BondUpdate);

impl BondUpdateDsl {
    /// Zero-cost reference cast from `&BondUpdate`. Relies on `repr(transparent)`.
    pub fn from_ref(update: &BondUpdate) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(update as *const BondUpdate as *const Self) }
    }
}

impl FromIr<BondUpdate> for BondUpdateDsl {
    type Ctx = ();

    fn from_ir(update: &BondUpdate, _ctx: &Self::Ctx) -> Self {
        Self(update.clone())
    }
}

impl IntoIr<BondUpdate> for BondUpdateDsl {
    type Ctx = ();

    fn into_ir(self, _ctx: &Self::Ctx) -> BondUpdate {
        self.0
    }
}

impl FromStr for BondUpdateDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_bond_update(s)
    }
}

impl FromStr for BondUpdate {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(BondUpdateDsl::from_str(s)?.into_ir(&()))
    }
}

impl Display for BondUpdate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        BondUpdateDsl::from_ref(self).fmt(f)
    }
}

impl Display for BondUpdateDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let update = &self.0;
        if let Some(order) = &update.order {
            match order {
                NumForm::Undetermined => write!(f, "*")?,
                NumForm::Lit(n) => write!(f, "{}", n)?,
                value => fmt_value(f, value)?,
            }
        }
        if let Some(charge) = &update.charge {
            if charge.is_undetermined() {
                write!(f, "#c*")?;
            } else {
                fmt_charge(f, charge)?;
            }
        }
        if let Some(unpaired_electrons) = &update.unpaired_electrons.count {
            fmt_update_value_field(f, "#u", unpaired_electrons)?;
        }
        if let Some(multiplicity) = &update.unpaired_electrons.multiplicity {
            fmt_update_value_field(f, "#s", multiplicity)?;
        }
        for c in update.constraints.iter() {
            if c.is_undetermined() {
                fmt_undetermined_constraint(f, c)?;
            } else {
                fmt_constraint(f, c)?;
            }
        }
        Ok(())
    }
}

impl<'de> FromEdn<'de> for BondUpdateDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s.parse().map_err(|e| DeError::subgrammar("bond-update", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for BondUpdateDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

pub fn parse_bond_update(input: &str) -> Result<BondUpdateDsl, ParseError> {
    bond_update.parse(input).map_err(|e| e.into_inner())
}

fn bond_update(i: &mut &str) -> PResult<BondUpdateDsl> {
    let order = delimited(multispace0, opt(value), multispace0).parse_next(i)?;
    let preds: Vec<BondPredicate> =
        repeat(0.., terminated(bond_predicate, multispace0)).parse_next(i)?;
    let mut update = BondUpdate {
        order,
        ..Default::default()
    };
    apply_update_predicates(&mut update, preds).map_err(ErrMode::Cut)?;
    Ok(BondUpdateDsl(update))
}

fn apply_update_predicates(
    update: &mut BondUpdate,
    preds: Vec<BondPredicate>,
) -> Result<(), ParseError> {
    for pred in preds {
        match pred {
            BondPredicate::Charge(value) => {
                if update.charge.replace(value).is_some() {
                    return Err(ParseError::DuplicateBondPredicate("#c".to_string()));
                }
            }
            BondPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Count(value)) => {
                if update.unpaired_electrons.count.replace(value).is_some() {
                    return Err(ParseError::DuplicateBondPredicate("#u".to_string()));
                }
            }
            BondPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Multiplicity(value)) => {
                if update
                    .unpaired_electrons
                    .multiplicity
                    .replace(value)
                    .is_some()
                {
                    return Err(ParseError::DuplicateBondPredicate("#s".to_string()));
                }
            }
            BondPredicate::Constraint(constraint) => {
                if update.constraints.contains(constraint.key()) {
                    return Err(ParseError::DuplicateBondPredicate(
                        constraint_tag(&constraint).to_string(),
                    ));
                }
                update.constraints.set(constraint);
            }
        }
    }
    Ok(())
}

fn fmt_undetermined_constraint(f: &mut fmt::Formatter<'_>, c: &BondConstraintAst) -> fmt::Result {
    match c {
        BondConstraintAst::RingMembership(membership) => match membership.scope {
            RingScope::All => write!(f, "#R*"),
            RingScope::Size(size) => write!(f, "#R({})*", size),
        },
        _ => write!(f, "{}*", constraint_tag(c)),
    }
}

fn fmt_update_value_field(f: &mut fmt::Formatter<'_>, prefix: &str, v: &NumForm) -> fmt::Result {
    if v.is_undetermined() {
        write!(f, "{}*", prefix)
    } else {
        match v {
            NumForm::Lit(1) => write!(f, "{}", prefix),
            NumForm::Lit(n) => write!(f, "{}{}", prefix, n),
            value => {
                write!(f, "{}", prefix)?;
                fmt_value(f, value)
            }
        }
    }
}

fn constraint_tag(c: &BondConstraintAst) -> &'static str {
    match c {
        BondConstraintAst::Aromatic(_) => "#a",
        BondConstraintAst::RingMembership(..) => "#R",
        BondConstraintAst::CisTransStereo(_) => "#C",
    }
}

/// One predicate from a bond-string; the parser yields a `Vec` of these
/// and the applier folds them into the `BondForm`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BondPredicate {
    Charge(NumForm),
    UnpairedElectrons(UnpairedElectronsPredicate),
    Constraint(BondConstraintAst),
}

fn bond_predicate(i: &mut &str) -> PResult<BondPredicate> {
    let start = *i;
    let prefix: &str = take(2usize).parse_next(i)?;
    match prefix {
        "#c" => charge.map(BondPredicate::Charge).parse_next(i),
        "#u" => optional_value
            .map(|v| BondPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Count(v)))
            .parse_next(i),
        "#s" => optional_value
            .map(|v| BondPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Multiplicity(v)))
            .parse_next(i),
        "#a" => boolean
            .map(|b| BondPredicate::Constraint(BondConstraintAst::Aromatic(b.0)))
            .parse_next(i),
        "#R" => ring_membership
            .map(|m| BondPredicate::Constraint(BondConstraintAst::RingMembership(m)))
            .parse_next(i),
        "#C" => (|i: &mut &str| cis_trans_stereo_config(i))
            .map(|c| BondPredicate::Constraint(BondConstraintAst::CisTransStereo(c)))
            .parse_next(i),
        p if p.starts_with('#') => Err(ErrMode::Cut(ParseError::UnknownBondPredicate(
            p.to_string(),
        ))),
        _ => Err(ErrMode::Cut(ParseError::TrailingInput(start.to_string()))),
    }
}

fn apply_predicates(form: &mut BondDsl, preds: Vec<BondPredicate>) -> Result<(), ParseError> {
    let ast = &mut form.0;
    for pred in preds {
        match pred {
            BondPredicate::Charge(v) => {
                if !matches!(ast.charge, NumForm::Undetermined) {
                    return Err(ParseError::DuplicateBondPredicate("#c".to_string()));
                }
                ast.charge = v;
            }
            BondPredicate::UnpairedElectrons(predicate) => {
                apply_unpaired_electrons_predicate(
                    &mut ast.unpaired_electrons,
                    predicate,
                    ParseError::DuplicateBondPredicate,
                )?;
            }
            BondPredicate::Constraint(c) => {
                if ast.constraints.contains(c.key()) {
                    return Err(ParseError::DuplicateBondPredicate(
                        constraint_tag(&c).to_string(),
                    ));
                }
                ast.constraints.set(c);
            }
        }
    }
    Ok(())
}

fn fmt_bond_form(f: &mut fmt::Formatter<'_>, ast: &BondForm) -> fmt::Result {
    match &ast.order {
        NumForm::Lit(n) => write!(f, "{}", n)?,
        NumForm::Undetermined => write!(f, "*")?,
        v => fmt_value(f, v)?,
    }

    fmt_charge(f, &ast.charge)?;
    fmt_unpaired_electrons(f, &ast.unpaired_electrons)
}

fn fmt_constraint(f: &mut fmt::Formatter<'_>, c: &BondConstraintAst) -> fmt::Result {
    match c {
        BondConstraintAst::Aromatic(BooleanForm::Lit(true)) => write!(f, "#a"),
        BondConstraintAst::Aromatic(BooleanForm::Lit(false)) => write!(f, "#a!"),
        BondConstraintAst::Aromatic(BooleanForm::Undetermined) => Ok(()),
        BondConstraintAst::RingMembership(m) => fmt_ring_membership(f, m),
        BondConstraintAst::CisTransStereo(CisTransStereoForm::Undetermined) => Ok(()),
        BondConstraintAst::CisTransStereo(c) => {
            write!(f, "#C")?;
            fmt_cis_trans_stereo_config(f, c)
        }
    }
}

pub(crate) fn lower_bond(ast: &mut BondForm, cfg: &BondDefaults) {
    // Exhaustive destructure: adding a new BondForm field is a compile error
    // here, forcing the author to decide how lowering should handle it.
    let BondForm {
        order: _,
        charge,
        unpaired_electrons,
        constraints,
    } = ast;

    if matches!(
        (&cfg.charge, &*charge),
        (NumericDefault::Zero, NumForm::Lit(0))
    ) {
        *charge = NumForm::Undetermined;
    }
    lower_unpaired_electrons(unpaired_electrons, cfg.unpaired_electrons, cfg.multiplicity);
    lower_bond_constraints(constraints, cfg);
}

pub(crate) fn raise_bond(ast: &mut BondForm, cfg: &BondDefaults) {
    // Exhaustive destructure: adding a new BondForm field is a compile error
    // here, forcing the author to decide how raising should handle it.
    let BondForm {
        order: _,
        charge,
        unpaired_electrons,
        constraints,
    } = ast;

    if matches!(*charge, NumForm::Undetermined) {
        *charge = match cfg.charge {
            NumericDefault::Zero => NumForm::Lit(0),
            NumericDefault::Required => NumForm::Undetermined,
        };
    }
    raise_unpaired_electrons(unpaired_electrons, cfg.unpaired_electrons, cfg.multiplicity);
    raise_bond_constraints(constraints, cfg);
}

fn raise_bond_constraints(constraints: &mut BondConstraintsAst, cfg: &BondDefaults) {
    // CisTransStereo is the only defaulted bond constraint; Aromatic/RingMembership are pattern-only.
    if matches!(cfg.cis_trans_stereo, StereoDefault::NotStereo)
        && constraints
            .get(BondConstraintKey::CisTransStereo)
            .is_none_or(|c| c.is_undetermined())
    {
        constraints.set(BondConstraintAst::CisTransStereo(
            CisTransStereoForm::NotStereo,
        ));
    }
}

fn lower_bond_constraints(constraints: &mut BondConstraintsAst, cfg: &BondDefaults) {
    // Elide the default CisTransStereo (NotStereo); Aromatic/RingMembership are pattern-only.
    if matches!(cfg.cis_trans_stereo, StereoDefault::NotStereo)
        && constraints.get(BondConstraintKey::CisTransStereo)
            == Some(&BondConstraintAst::CisTransStereo(
                CisTransStereoForm::NotStereo,
            ))
    {
        constraints.remove(BondConstraintKey::CisTransStereo);
    }
}

/// Surface DSL wrapper around `BondConstraintAst`. EDN form: the keyword
/// `:aromatic` (flag variant, no value) or a single-key map
/// `{:ring-membership {:size? <int> :count <value>}}` / `{:cis-trans-stereo …}`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct BondConstraintDsl(pub BondConstraintAst);

impl FromIr<BondConstraintAst> for BondConstraintDsl {
    type Ctx = ();

    fn from_ir(ast: &BondConstraintAst, _ctx: &Self::Ctx) -> Self {
        Self(ast.clone())
    }
}

impl IntoIr<BondConstraintAst> for BondConstraintDsl {
    type Ctx = ();

    fn into_ir(self, _ctx: &Self::Ctx) -> BondConstraintAst {
        self.0
    }
}

impl<'de> FromEdn<'de> for BondConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Map(m) if m.len() == 1 => {
                let (k, v) = m.iter().next().unwrap();
                let Edn::Keyword(key) = k else {
                    return Err(DeError::TypeMismatch {
                        expected: "keyword key",
                        got: k.kind(),
                        path: vec!["bond-constraint".into()],
                    });
                };
                let c = match key.name() {
                    "aromatic" => BondConstraintAst::Aromatic(BooleanDsl::from_edn(v)?.0),
                    "ring-membership" => {
                        BondConstraintAst::RingMembership(RingMembershipDsl::from_edn(v)?.0)
                    }
                    "cis-trans-stereo" => BondConstraintAst::CisTransStereo(
                        CisTransStereoDsl::from_edn(v)?.into_ir(&()),
                    ),
                    other => {
                        return Err(DeError::UnknownField {
                            key: other.to_string(),
                            path: vec!["bond-constraint".into()],
                        });
                    }
                };
                Ok(Self(c))
            }
            other => Err(DeError::TypeMismatch {
                expected: "{:aromatic …} / {:ring-membership …} / {:cis-trans-stereo …}",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for BondConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        match &self.0 {
            BondConstraintAst::Aromatic(b) => single_key_map("aromatic", BooleanDsl(*b).to_edn()),
            BondConstraintAst::RingMembership(m) => {
                single_key_map("ring-membership", RingMembershipDsl(m.clone()).to_edn())
            }
            BondConstraintAst::CisTransStereo(c) => single_key_map(
                "cis-trans-stereo",
                CisTransStereoDsl::from_ir(c, &()).to_edn(),
            ),
        }
    }
}

#[cfg(test)]
mod tests {

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::read_string;

    use super::*;
    use crate::bond_dsl;
    use crate::ir::constraint::{BondConstraintsAst, RingScope};
    use crate::ir::spin::{UnpairedElectronsForm, UnpairedElectronsUpdate};
    use crate::ir::stereo::StereoCoset;

    #[rustfmt::skip]
    #[rstest]
    #[case::single("1", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::new() }))]
    #[case::double("2", BondDsl(BondForm { order: NumForm::Lit(2), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::new() }))]
    #[case::triple("3", BondDsl(BondForm { order: NumForm::Lit(3), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::new() }))]
    #[case::quadruple("4", BondDsl(BondForm { order: NumForm::Lit(4), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::new() }))]
    #[case::single_whitespace("  1  ", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::new() }))]
    #[case::single_pos_charge("1#c+2", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(2), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::new() }))]
    #[case::single_neg_charge("1#c-2", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(-2), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::new() }))]
    #[case::single_zero_charge("1#c0", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::new() }))]
    #[case::single_plus_only("1#c+", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::new() }))]
    #[case::single_minus_only("1#c-", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(-1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::new() }))]
    #[case::double_unpaired("2#u3", BondDsl(BondForm { order: NumForm::Lit(2), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(3), multiplicity: NumForm::Undetermined }, constraints: BondConstraintsAst::new() }))]
    #[case::single_u_only("1#u", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Undetermined }, constraints: BondConstraintsAst::new() }))]
    #[case::single_mult("1#s2", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(2) }, constraints: BondConstraintsAst::new() }))]
    #[case::single_s_only("1#s", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(1) }, constraints: BondConstraintsAst::new() }))]
    #[case::double_charge_unpaired("2#c+#u2", BondDsl(BondForm { order: NumForm::Lit(2), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Undetermined }, constraints: BondConstraintsAst::new() }))]
    #[case::double_charge_mult("2#c-1#s3", BondDsl(BondForm { order: NumForm::Lit(2), charge: NumForm::Lit(-1), unpaired_electrons: UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(3) }, constraints: BondConstraintsAst::new() }))]
    #[case::aromatic("1#a", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanForm::Lit(true))]) }))]
    #[case::aromatic_plus("1#a+", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanForm::Lit(true))]) }))]
    #[case::aromatic_undetermined("1#a*", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanForm::Undetermined)]) }))]
    #[case::not_aromatic("1#a!", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanForm::Lit(false))]) }))]
    #[case::charged_aromatic("1#c+#a", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanForm::Lit(true))]) }))]
    #[case::ring_membership_all("1#R2", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::All, NumForm::Lit(2))]) }))]
    #[case::ring_membership_all_bare("1#R", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::All, NumForm::Lit(1))]) }))]
    #[case::ring_membership_all_plus("1#R+", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::All, NumForm::RangeFrom(1))]) }))]
    #[case::ring_bang("1#R!", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::All, NumForm::Lit(0))]) }))]
    #[case::ring_zero("1#R0", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::All, NumForm::Lit(0))]) }))]
    #[case::ring_membership_all_star("1#R*", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::All, NumForm::Undetermined)]) }))]
    #[case::ring_membership_size("1#R(6)", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::Size(6), 1)]) }))]
    #[case::ring_membership_size_conj("1#R(5)#R(6)", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::from_iter([BondConstraintAst::ring_membership(RingScope::Size(5), 1), BondConstraintAst::ring_membership(RingScope::Size(6), 1)]) }))]
    #[case::whitespace_before_predicate("2 #c+", BondDsl(BondForm { order: NumForm::Lit(2), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::new() }))]
    #[case::whitespace_between_predicates("2#c+ #a", BondDsl(BondForm { order: NumForm::Lit(2), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanForm::Lit(true))]) }))]
    #[case::whitespace_surrounding_predicates("  2  #c+  #a  ", BondDsl(BondForm { order: NumForm::Lit(2), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanForm::Lit(true))]) }))]
    fn test_parse_bond(#[case] input: &str, #[case] expected: BondDsl) {
        let result = bond.parse(input);
        assert!(result.is_ok(), "{:?} should succeed, got {:?}", input, result.clone().unwrap_err());
        let form = result.unwrap();
        assert_eq!(form, expected);
    }

    #[rstest]
    #[case::empty("", ParseError::Syntax)]
    #[case::unknown_pred("1#x", ParseError::UnknownBondPredicate("#x".to_string()))]
    #[case::dup_charge("1#c+#c-", ParseError::DuplicateBondPredicate("#c".to_string()))]
    #[case::dup_unpaired("1#u2#u3", ParseError::DuplicateBondPredicate("#u".to_string()))]
    #[case::dup_multiplicity("1#s1#s2", ParseError::DuplicateBondPredicate("#s".to_string()))]
    #[case::dup_aromatic("1#a#a", ParseError::DuplicateBondPredicate("#a".to_string()))]
    #[case::trailing("1#c+ foo", ParseError::TrailingInput("foo".to_string()))]
    fn test_parse_bond_error(#[case] input: &str, #[case] expected: ParseError) {
        let result = bond.parse(input);
        assert!(
            result.is_err(),
            "{:?} should fail, got {:?}",
            input,
            result.unwrap()
        );
        let err = result.unwrap_err().into_inner();
        assert_eq!(
            err, expected,
            "{:?} should fail with {:?}, got {:?}",
            input, expected, err
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", BondUpdateDsl(BondUpdate::default()))]
    #[case::order_only("2", BondUpdateDsl(BondUpdate { order: Some(NumForm::Lit(2)), ..Default::default() }))]
    #[case::order_undetermined("*", BondUpdateDsl(BondUpdate { order: Some(NumForm::Undetermined), ..Default::default() }))]
    #[case::charge_only("#c-1", BondUpdateDsl(BondUpdate { charge: Some(NumForm::Lit(-1)), ..Default::default() }))]
    #[case::unpaired_electrons_unpaired("#u2", BondUpdateDsl(BondUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Lit(2)), multiplicity: None }, ..Default::default() }))]
    #[case::unpaired_electrons_multiplicity("#s1", BondUpdateDsl(BondUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) }, ..Default::default() }))]
    #[case::explicit_undetermined("*#c*#u*#s*", BondUpdateDsl(BondUpdate { order: Some(NumForm::Undetermined), charge: Some(NumForm::Undetermined), unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Undetermined), multiplicity: Some(NumForm::Undetermined) }, constraints: Default::default() }))]
    #[case::order_and_pred("1#a", BondUpdateDsl(BondUpdate { order: Some(NumForm::Lit(1)), constraints: BondConstraintsAst::from(BondConstraintAst::Aromatic(BooleanForm::Lit(true))), ..Default::default() }))]
    #[case::constraint_removal("#R(6)*", BondUpdateDsl(BondUpdate { constraints: BondConstraintsAst::from(BondConstraintAst::ring_membership(RingScope::Size(6), NumForm::Undetermined)), ..Default::default() }))]
    fn test_parse_bond_update(#[case] input: &str, #[case] expected: BondUpdateDsl) {
        assert_eq!(parse_bond_update(input).unwrap(), expected);
    }

    #[rstest]
    #[case::dup_charge("#c+#c-", ParseError::DuplicateBondPredicate("#c".to_string()))]
    #[case::dup_undetermined_charge("#c*#c-", ParseError::DuplicateBondPredicate("#c".to_string()))]
    #[case::dup_undetermined_unpaired("#u*#u2", ParseError::DuplicateBondPredicate("#u".to_string()))]
    #[case::unknown_pred("1#x", ParseError::UnknownBondPredicate("#x".to_string()))]
    fn test_parse_bond_update_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(parse_bond_update(input).unwrap_err(), expected);
    }

    #[rstest]
    #[case::duplicate_charge("#c+#c-", ParseError::DuplicateBondPredicate("#c".to_string()))]
    fn test_bond_update_from_str_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(input.parse::<BondUpdate>().unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge_only(r##""#c-1""##, BondUpdateDsl(BondUpdate { charge: Some(NumForm::Lit(-1)), ..Default::default() }))]
    fn test_bond_update_dsl_from_edn(#[case] input: &str, #[case] expected: BondUpdateDsl) {
        assert_eq!(
            BondUpdateDsl::from_edn(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::non_string("1")]
    fn test_bond_update_dsl_from_edn_error(#[case] input: &str) {
        assert!(matches!(
            BondUpdateDsl::from_edn(&read_string(input).unwrap()),
            Err(DeError::TypeMismatch {
                expected: "string",
                ..
            })
        ));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::order(BondUpdateDsl(BondUpdate { order: Some(NumForm::Lit(1)), ..Default::default() }), r##""1""##)]
    #[case::order_undetermined(BondUpdateDsl(BondUpdate { order: Some(NumForm::Undetermined), ..Default::default() }), r##""*""##)]
    #[case::charge_only(BondUpdateDsl(BondUpdate { charge: Some(NumForm::Lit(-1)), ..Default::default() }), r##""#c-""##)]
    #[case::unpaired_electrons_multiplicity(BondUpdateDsl(BondUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) }, ..Default::default() }), r##""#s""##)]
    #[case::explicit_undetermined(BondUpdateDsl(BondUpdate { order: Some(NumForm::Undetermined), charge: Some(NumForm::Undetermined), unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Undetermined), multiplicity: Some(NumForm::Undetermined) }, constraints: Default::default() }), r##""*#c*#u*#s*""##)]
    #[case::aromatic(BondUpdateDsl(BondUpdate { constraints: BondConstraintsAst::from(BondConstraintAst::Aromatic(BooleanForm::Lit(true))), ..Default::default() }), r##""#a""##)]
    #[case::aromatic_undetermined(BondUpdateDsl(BondUpdate { constraints: BondConstraintsAst::from(BondConstraintAst::Aromatic(BooleanForm::Undetermined)), ..Default::default() }), r##""#a*""##)]
    #[case::ring_size_removal(BondUpdateDsl(BondUpdate { constraints: BondConstraintsAst::from(BondConstraintAst::ring_membership(RingScope::Size(6), NumForm::Undetermined)), ..Default::default() }), r##""#R(6)*""##)]
    fn test_bond_update_dsl_to_edn(#[case] input: BondUpdateDsl, #[case] expected: &str) {
        assert_eq!(input.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::charge_pos("#c+2", BondPredicate::Charge(NumForm::Lit(2)))]
    #[case::charge_neg("#c-2", BondPredicate::Charge(NumForm::Lit(-2)))]
    #[case::charge_plus("#c+", BondPredicate::Charge(NumForm::Lit(1)))]
    #[case::charge_minus("#c-", BondPredicate::Charge(NumForm::Lit(-1)))]
    #[case::charge_zero("#c0", BondPredicate::Charge(NumForm::Lit(0)))]
    #[case::charge_undetermined("#c*", BondPredicate::Charge(NumForm::Undetermined))]
    #[case::unpaired_electrons("#u2", BondPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Count(NumForm::Lit(2))))]
    #[case::unpaired_omit("#u", BondPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Count(NumForm::Lit(1))))]
    #[case::unpaired_undetermined("#u*", BondPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Count(NumForm::Undetermined)))]
    #[case::multiplicity("#s3", BondPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Multiplicity(NumForm::Lit(3))))]
    #[case::multiplicity_omit("#s", BondPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Multiplicity(NumForm::Lit(1))))]
    #[case::multiplicity_undetermined("#s*", BondPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Multiplicity(NumForm::Undetermined)))]
    #[case::aromatic("#a", BondPredicate::Constraint(BondConstraintAst::Aromatic(BooleanForm::Lit(true))))]
    #[case::aromatic_plus("#a+", BondPredicate::Constraint(BondConstraintAst::Aromatic(BooleanForm::Lit(true))))]
    #[case::aromatic_false("#a!", BondPredicate::Constraint(BondConstraintAst::Aromatic(BooleanForm::Lit(false))))]
    #[case::aromatic_undetermined("#a*", BondPredicate::Constraint(BondConstraintAst::Aromatic(BooleanForm::Undetermined)))]
    #[case::ring_membership_all("#R2", BondPredicate::Constraint(BondConstraintAst::ring_membership(RingScope::All, NumForm::Lit(2))))]
    #[case::ring_membership_all_plus("#R+", BondPredicate::Constraint(BondConstraintAst::ring_membership(RingScope::All, NumForm::RangeFrom(1))))]
    #[case::ring_membership_zero("#R0", BondPredicate::Constraint(BondConstraintAst::ring_membership(RingScope::All, NumForm::Lit(0))))]
    #[case::ring_membership_all_undetermined("#R*", BondPredicate::Constraint(BondConstraintAst::ring_membership(RingScope::All, NumForm::Undetermined)))]
    #[case::ring_membership_size("#R(6)", BondPredicate::Constraint(BondConstraintAst::ring_membership(RingScope::Size(6), 1)))]
    #[case::ring_membership_size_undetermined("#R(6)*", BondPredicate::Constraint(BondConstraintAst::ring_membership(RingScope::Size(6), NumForm::Undetermined)))]
    #[case::cis_trans_stereo_undetermined("#C*", BondPredicate::Constraint(BondConstraintAst::CisTransStereo(CisTransStereoForm::Undetermined)))]
    #[case::cis_trans_stereo_plus("#C+", BondPredicate::Constraint(BondConstraintAst::CisTransStereo(CisTransStereoForm::Stereo(StereoCoset::Undetermined))))]
    #[case::cis_trans_stereo_not_stereo("#C!", BondPredicate::Constraint(BondConstraintAst::CisTransStereo(CisTransStereoForm::NotStereo)))]
    #[case::cis_trans_stereo("#C1", BondPredicate::Constraint(BondConstraintAst::CisTransStereo(CisTransStereoForm::Stereo(StereoCoset::Lit(1)))))]
    fn test_bond_predicate(#[case] input: &str, #[case] expected: BondPredicate) {
        let result = bond_predicate.parse(input);
        assert!(result.is_ok(), "{input:?} should succeed, got {:?}", result.unwrap_err());
        let pred = result.unwrap();
        assert_eq!(pred, expected);
    }

    #[rstest]
    #[case::unknown("#x", ParseError::UnknownBondPredicate("#x".to_string()))]
    #[case::unknown_tag("#z", ParseError::UnknownBondPredicate("#z".to_string()))]
    #[case::trailing_no_tag("fo", ParseError::TrailingInput("fo".to_string()))]
    fn test_bond_predicate_error(#[case] input: &str, #[case] expected: ParseError) {
        let result = bond_predicate.parse(input);
        assert!(
            result.is_err(),
            "{input:?} should fail, got {:?}",
            result.unwrap()
        );
        let err = result.unwrap_err().into_inner();
        assert_eq!(err, expected);
    }

    #[rstest]
    fn test_bond_dsl_from_ast() {
        let mut ast = BondForm::new(NumForm::Lit(1));
        ast.charge = NumForm::Lit(0);
        ast.unpaired_electrons = UnpairedElectronsForm::from((0_u8, 1_u8));
        let cfg = BondDefaults::zeroed();
        let dsl = BondDsl::from_ir(&ast, &cfg);
        assert_eq!(dsl.0.charge, NumForm::Undetermined);
        assert_eq!(dsl.0.unpaired_electrons, UnpairedElectronsForm::default());
    }

    #[rstest]
    fn test_bond_dsl_into_ast() {
        let dsl = BondDsl(BondForm::new(NumForm::Lit(1)));
        let cfg = BondDefaults::zeroed();
        let ast = dsl.into_ir(&cfg);
        assert_eq!(ast.charge, NumForm::Lit(0));
        assert_eq!(
            ast.unpaired_electrons,
            UnpairedElectronsForm::from((0_u8, 1_u8))
        );
    }

    #[rstest]
    fn test_bond_dsl_roundtrip_zeroed() {
        let input = BondDsl(BondForm::new(NumForm::Lit(2)));
        let cfg = BondDefaults::zeroed();
        let raised = input.clone().into_ir(&cfg);
        let lowered = BondDsl::from_ir(&raised, &cfg);
        assert_eq!(input, lowered);
    }

    #[rstest]
    #[case::single(r##""1""##)]
    #[case::aromatic(r##""1#a""##)]
    #[case::ring_membership_all(r##""2#R+""##)]
    #[case::with_cis_trans_stereo(r##""2#C1""##)]
    fn test_bond_dsl_from_edn_str_matches_from_edn(#[case] input: &str) {
        let via_stream = BondDsl::from_edn_str(input).unwrap();
        let tree = read_string(input).unwrap();
        let via_tree = BondDsl::from_edn(&tree).unwrap();
        assert_eq!(via_stream, via_tree);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(":single", BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::new() })]
    #[case::double(":double", BondForm { order: NumForm::Lit(2), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::new() })]
    #[case::triple(":triple", BondForm { order: NumForm::Lit(3), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::new() })]
    #[case::quadruple(":quadruple", BondForm { order: NumForm::Lit(4), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::new() })]
    #[case::aromatic(":aromatic", BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsAst::from_iter([BondConstraintAst::Aromatic(BooleanForm::Lit(true))]) })]
    fn test_bond_dsl_keyword_shorthand(#[case] input: &str, #[case] expected: BondForm) {
        let edn = read_string(input).unwrap();
        let dsl = BondDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.0, expected);
    }

    #[rstest]
    fn test_bond_dsl_keyword_shorthand_error() {
        let edn = read_string(":bogus").unwrap();
        let err = BondDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraintAst::Aromatic(BooleanForm::Lit(true)), "{:aromatic true}")]
    #[case::aromatic_false(BondConstraintAst::Aromatic(BooleanForm::Lit(false)), "{:aromatic false}")]
    #[case::ring_membership_all(BondConstraintAst::ring_membership(RingScope::All, NumForm::Lit(1)), "{:ring-membership {:count 1}}")]
    #[case::ring_membership_all_undetermined(BondConstraintAst::ring_membership(RingScope::All, NumForm::Undetermined), "{:ring-membership {:count :undetermined}}")]
    #[case::ring_membership_size(BondConstraintAst::ring_membership(RingScope::Size(6), 1), "{:ring-membership {:size 6 :count 1}}")]
    #[case::ring_membership_size_count_set(BondConstraintAst::ring_membership(RingScope::Size(6), NumForm::lit_set([5, 6])), "{:ring-membership {:size 6 :count [5 6]}}")]
    #[case::cis_trans_stereo_undetermined(BondConstraintAst::CisTransStereo(CisTransStereoForm::Undetermined), "{:cis-trans-stereo :undetermined}")]
    #[case::cis_trans_stereo_not_stereo(BondConstraintAst::CisTransStereo(CisTransStereoForm::NotStereo), "{:cis-trans-stereo :not-stereo}")]
    #[case::cis_trans_stereo_lit(BondConstraintAst::CisTransStereo(CisTransStereoForm::Stereo(StereoCoset::Lit(1))), "{:cis-trans-stereo {:stereo 1}}")]
    #[case::cis_trans_stereo_coset_undetermined(BondConstraintAst::CisTransStereo(CisTransStereoForm::Stereo(StereoCoset::Undetermined)), "{:cis-trans-stereo {:stereo :undetermined}}")]
    #[case::cis_trans_stereo_set(BondConstraintAst::CisTransStereo(CisTransStereoForm::Stereo(StereoCoset::lit_set([1, 2]))), "{:cis-trans-stereo {:stereo [1 2]}}")]
    fn test_bond_constraint_dsl_roundtrip(
        #[case] input: BondConstraintAst,
        #[case] edn_source: &str,
    ) {
        let dsl = BondConstraintDsl::from_ir(&input, &());
        let edn = dsl.to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = BondConstraintDsl::from_edn(&edn).unwrap();
        assert_eq!(parsed.into_ir(&()), input, "parse-back mismatch");
    }

    #[rstest]
    #[case::ring_membership_all("1#R*", "1")]
    #[case::ring_membership_size("1#R(6)*", "1")]
    #[case::aromatic("1#a*", "1")]
    #[case::cis_trans_stereo("2#C*", "2")]
    fn test_bond_render_vacuous_constraints(#[case] input: &str, #[case] expected_canonical: &str) {
        let parsed: BondDsl = bond.parse(input).unwrap();
        assert_eq!(parsed.to_string(), expected_canonical);
        let reparsed: BondDsl = bond.parse(&parsed.to_string()).unwrap();
        assert!(
            reparsed.0.constraints.is_empty(),
            "vacuous constraint should be absent after render → reparse, got {:?}",
            reparsed.0.constraints,
        );
    }

    #[rstest]
    #[case::charge_before_multiplicity("2#s3#c+", "2#c+#s3")]
    #[case::aromatic_before_ring("2#R2#a", "2#a#R2")]
    #[case::stereo_before_ring("1#R2#C1", "1#C1#R2")]
    fn test_bond_render_canonical_order(#[case] input: &str, #[case] expected_canonical: &str) {
        let parsed: BondDsl = bond.parse(input).unwrap();
        assert_eq!(parsed.to_string(), expected_canonical);
    }

    #[rstest]
    #[case::wrong_shape(Edn::Int(3), DeError::TypeMismatch { expected: "{:aromatic …} / {:ring-membership …} / {:cis-trans-stereo …}", got: "int", path: vec![] })]
    #[case::unknown_key("{:bogus 1}", DeError::UnknownField { key: "bogus".to_string(), path: vec!["bond-constraint".into()] })]
    fn test_bond_constraint_dsl_error(#[case] input: Edn<'static>, #[case] expected: DeError) {
        let err = BondConstraintDsl::from_edn(&input).unwrap_err();
        assert_eq!(err, expected);
    }

    #[rstest]
    #[case::single("1")]
    #[case::double("2")]
    #[case::aromatic("1#a")]
    #[case::ring_membership_all("1#R2")]
    #[case::ring_membership_size("1#R(6)")]
    #[case::cis_trans_stereo("2#C1")]
    fn test_bond_form_from_str_to_string_roundtrip(#[case] s: &str) {
        let ast: BondForm = s.parse().unwrap();
        assert_eq!(ast.to_string(), s);
    }

    #[rstest]
    #[case::double(bond_dsl!("2"))]
    #[case::aromatic(bond_dsl!("1#a"))]
    #[case::ring_membership_all(bond_dsl!("1#R+"))]
    #[case::ring_membership_size(bond_dsl!("1#R(6)+"))]
    #[case::cis_trans_stereo(bond_dsl!("1#C1"))]
    #[case::cis_trans_stereo_coset_undetermined(bond_dsl!("1#C+"))]
    #[case::cis_trans_stereo_not_stereo(bond_dsl!("1#C!"))]
    #[case::cis_trans_stereo_set(bond_dsl!("1#C{1,2}"))]
    fn test_bond_form_to_edn_from_edn_roundtrip(#[case] input: BondForm) {
        let edn = input.to_edn();
        let parsed = BondForm::from_edn(&edn).unwrap();
        assert_eq!(parsed, input);
    }
}
