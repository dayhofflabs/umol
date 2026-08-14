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
use super::config::{BondDefaults, NumDefault, StereoDefault};
use super::constraint::RingMembershipDsl;
use super::edn_utils::single_key_map;
use super::error::{PResult, ParseError};
use super::num::{fmt_num, num};
use super::predicate::{
    apply_unpaired_electrons_predicate, charge, fmt_charge, fmt_ring_membership,
    fmt_unpaired_electrons, lower_unpaired_electrons, optional_value, raise_unpaired_electrons,
    ring_membership, UnpairedElectronsPredicate,
};
use super::stereo::{cis_trans_stereo_config, fmt_cis_trans_stereo_config, CisTransStereoDsl};
use crate::ir::bond::{BondForm, BondUpdate};
use crate::ir::boolean::BooleanForm;
use crate::ir::constraint::{
    BondConstraintForm, BondConstraintKey, BondConstraintsForm, RingScope,
};
use crate::ir::num::NumForm;
use crate::ir::spin::UnpairedElectronsForm;
use crate::ir::stereo::CisTransStereoForm;
use crate::ir::traits::{FromIr, IntoIr, Lattice};

/// Surface DSL wrapper around `BondForm`.
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BondDsl(pub BondForm);

impl BondDsl {
    /// Zero-cost reference cast from `&BondForm`. Relies on `repr(transparent)`.
    pub fn from_ref(form: &BondForm) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(form as *const BondForm as *const Self) }
    }
}

impl From<BondForm> for BondDsl {
    fn from(form: BondForm) -> Self {
        Self(form)
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
fn bond_keyword_for(form: &BondForm) -> Option<&'static str> {
    if !matches!(form.charge, NumForm::Undetermined)
        || form.unpaired_electrons != UnpairedElectronsForm::default()
    {
        return None;
    }
    let constraints: Vec<&BondConstraintForm> = form.constraints.iter().collect();
    match (&form.order, constraints.as_slice()) {
        (NumForm::Lit(1), []) => Some("single"),
        (NumForm::Lit(2), []) => Some("double"),
        (NumForm::Lit(3), []) => Some("triple"),
        (NumForm::Lit(4), []) => Some("quadruple"),
        (NumForm::Lit(1), [BondConstraintForm::Aromatic(BooleanForm::Lit(true))]) => {
            Some("aromatic")
        }
        _ => None,
    }
}

impl FromIr<BondForm> for BondDsl {
    type Context = BondDefaults;

    fn from_ir(form: &BondForm, context: &Self::Context) -> Self {
        let mut out = form.clone();
        lower_bond(&mut out, context);
        BondDsl(out)
    }
}

impl IntoIr<BondForm> for BondDsl {
    type Context = BondDefaults;

    fn into_ir(mut self, context: &Self::Context) -> BondForm {
        raise_bond(&mut self.0, context);
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
    let order = delimited(multispace0, num, multispace0).parse_next(i)?;
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
    type Context = ();

    fn from_ir(update: &BondUpdate, _context: &Self::Context) -> Self {
        Self(update.clone())
    }
}

impl IntoIr<BondUpdate> for BondUpdateDsl {
    type Context = ();

    fn into_ir(self, _context: &Self::Context) -> BondUpdate {
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
                value => fmt_num(f, value)?,
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
    let order = delimited(multispace0, opt(num), multispace0).parse_next(i)?;
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

fn fmt_undetermined_constraint(f: &mut fmt::Formatter<'_>, c: &BondConstraintForm) -> fmt::Result {
    match c {
        BondConstraintForm::RingMembership(membership) => match membership.scope {
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
                fmt_num(f, value)
            }
        }
    }
}

fn constraint_tag(c: &BondConstraintForm) -> &'static str {
    match c {
        BondConstraintForm::Aromatic(_) => "#a",
        BondConstraintForm::RingMembership(..) => "#R",
        BondConstraintForm::CisTransStereo(_) => "#C",
    }
}

/// One predicate from a bond-string; the parser yields a `Vec` of these
/// and the applier folds them into the `BondForm`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BondPredicate {
    Charge(NumForm),
    UnpairedElectrons(UnpairedElectronsPredicate),
    Constraint(BondConstraintForm),
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
            .map(|b| BondPredicate::Constraint(BondConstraintForm::Aromatic(b.0)))
            .parse_next(i),
        "#R" => ring_membership
            .map(|m| BondPredicate::Constraint(BondConstraintForm::RingMembership(m)))
            .parse_next(i),
        "#C" => (|i: &mut &str| cis_trans_stereo_config(i))
            .map(|c| BondPredicate::Constraint(BondConstraintForm::CisTransStereo(c)))
            .parse_next(i),
        p if p.starts_with('#') => Err(ErrMode::Cut(ParseError::UnknownBondPredicate(
            p.to_string(),
        ))),
        _ => Err(ErrMode::Cut(ParseError::TrailingInput(start.to_string()))),
    }
}

fn apply_predicates(dsl: &mut BondDsl, preds: Vec<BondPredicate>) -> Result<(), ParseError> {
    let bond = &mut dsl.0;
    for pred in preds {
        match pred {
            BondPredicate::Charge(v) => {
                if !matches!(bond.charge, NumForm::Undetermined) {
                    return Err(ParseError::DuplicateBondPredicate("#c".to_string()));
                }
                bond.charge = v;
            }
            BondPredicate::UnpairedElectrons(predicate) => {
                apply_unpaired_electrons_predicate(
                    &mut bond.unpaired_electrons,
                    predicate,
                    ParseError::DuplicateBondPredicate,
                )?;
            }
            BondPredicate::Constraint(c) => {
                if bond.constraints.contains(c.key()) {
                    return Err(ParseError::DuplicateBondPredicate(
                        constraint_tag(&c).to_string(),
                    ));
                }
                bond.constraints.set(c);
            }
        }
    }
    Ok(())
}

fn fmt_bond_form(f: &mut fmt::Formatter<'_>, form: &BondForm) -> fmt::Result {
    match &form.order {
        NumForm::Lit(n) => write!(f, "{}", n)?,
        NumForm::Undetermined => write!(f, "*")?,
        v => fmt_num(f, v)?,
    }

    fmt_charge(f, &form.charge)?;
    fmt_unpaired_electrons(f, &form.unpaired_electrons)
}

fn fmt_constraint(f: &mut fmt::Formatter<'_>, c: &BondConstraintForm) -> fmt::Result {
    match c {
        BondConstraintForm::Aromatic(BooleanForm::Lit(true)) => write!(f, "#a"),
        BondConstraintForm::Aromatic(BooleanForm::Lit(false)) => write!(f, "#a!"),
        BondConstraintForm::Aromatic(BooleanForm::Undetermined) => Ok(()),
        BondConstraintForm::RingMembership(m) => fmt_ring_membership(f, m),
        BondConstraintForm::CisTransStereo(CisTransStereoForm::Undetermined) => Ok(()),
        BondConstraintForm::CisTransStereo(c) => {
            write!(f, "#C")?;
            fmt_cis_trans_stereo_config(f, c)
        }
    }
}

pub(crate) fn lower_bond(bond: &mut BondForm, cfg: &BondDefaults) {
    // Exhaustive destructure: adding a new BondForm field is a compile error
    // here, forcing the author to decide how lowering should handle it.
    let BondForm {
        order: _,
        charge,
        unpaired_electrons,
        constraints,
    } = bond;

    if matches!((&cfg.charge, &*charge), (NumDefault::Zero, NumForm::Lit(0))) {
        *charge = NumForm::Undetermined;
    }
    lower_unpaired_electrons(unpaired_electrons, cfg.unpaired_electrons, cfg.multiplicity);
    lower_bond_constraints(constraints, cfg);
}

pub(crate) fn raise_bond(bond: &mut BondForm, cfg: &BondDefaults) {
    // Exhaustive destructure: adding a new BondForm field is a compile error
    // here, forcing the author to decide how raising should handle it.
    let BondForm {
        order: _,
        charge,
        unpaired_electrons,
        constraints,
    } = bond;

    if matches!(*charge, NumForm::Undetermined) {
        *charge = match cfg.charge {
            NumDefault::Zero => NumForm::Lit(0),
            NumDefault::Required => NumForm::Undetermined,
        };
    }
    raise_unpaired_electrons(unpaired_electrons, cfg.unpaired_electrons, cfg.multiplicity);
    raise_bond_constraints(constraints, cfg);
}

fn raise_bond_constraints(constraints: &mut BondConstraintsForm, cfg: &BondDefaults) {
    // CisTransStereo is the only defaulted bond constraint; Aromatic/RingMembership are pattern-only.
    if matches!(cfg.cis_trans_stereo, StereoDefault::NotStereo)
        && constraints
            .get(BondConstraintKey::CisTransStereo)
            .is_none_or(|c| c.is_undetermined())
    {
        constraints.set(BondConstraintForm::CisTransStereo(
            CisTransStereoForm::NotStereo,
        ));
    }
}

fn lower_bond_constraints(constraints: &mut BondConstraintsForm, cfg: &BondDefaults) {
    // Elide the default CisTransStereo (NotStereo); Aromatic/RingMembership are pattern-only.
    if matches!(cfg.cis_trans_stereo, StereoDefault::NotStereo)
        && constraints.get(BondConstraintKey::CisTransStereo)
            == Some(&BondConstraintForm::CisTransStereo(
                CisTransStereoForm::NotStereo,
            ))
    {
        constraints.remove(BondConstraintKey::CisTransStereo);
    }
}

/// Surface DSL wrapper around `BondConstraintForm`. EDN form: the keyword
/// `:aromatic` (flag variant, no value) or a single-key map
/// `{:ring-membership {:size? <int> :count <value>}}` / `{:cis-trans-stereo …}`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct BondConstraintDsl(pub BondConstraintForm);

impl FromIr<BondConstraintForm> for BondConstraintDsl {
    type Context = ();

    fn from_ir(form: &BondConstraintForm, _context: &Self::Context) -> Self {
        Self(form.clone())
    }
}

impl IntoIr<BondConstraintForm> for BondConstraintDsl {
    type Context = ();

    fn into_ir(self, _context: &Self::Context) -> BondConstraintForm {
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
                    "aromatic" => BondConstraintForm::Aromatic(BooleanDsl::from_edn(v)?.0),
                    "ring-membership" => {
                        BondConstraintForm::RingMembership(RingMembershipDsl::from_edn(v)?.0)
                    }
                    "cis-trans-stereo" => BondConstraintForm::CisTransStereo(
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
            BondConstraintForm::Aromatic(b) => single_key_map("aromatic", BooleanDsl(*b).to_edn()),
            BondConstraintForm::RingMembership(m) => {
                single_key_map("ring-membership", RingMembershipDsl(m.clone()).to_edn())
            }
            BondConstraintForm::CisTransStereo(c) => single_key_map(
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
    use crate::ir::constraint::{BondConstraintsForm, RingScope};
    use crate::ir::spin::{UnpairedElectronsForm, UnpairedElectronsUpdate};
    use crate::ir::stereo::StereoCoset;

    #[rustfmt::skip]
    #[rstest]
    #[case::single("1", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() }))]
    #[case::double("2", BondDsl(BondForm { order: NumForm::Lit(2), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() }))]
    #[case::triple("3", BondDsl(BondForm { order: NumForm::Lit(3), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() }))]
    #[case::quadruple("4", BondDsl(BondForm { order: NumForm::Lit(4), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() }))]
    #[case::single_whitespace("  1  ", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() }))]
    #[case::single_pos_charge("1#c+2", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(2), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() }))]
    #[case::single_neg_charge("1#c-2", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(-2), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() }))]
    #[case::single_zero_charge("1#c0", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() }))]
    #[case::single_plus_only("1#c+", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() }))]
    #[case::single_minus_only("1#c-", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(-1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() }))]
    #[case::double_unpaired("2#u3", BondDsl(BondForm { order: NumForm::Lit(2), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(3), multiplicity: NumForm::Undetermined }, constraints: BondConstraintsForm::new() }))]
    #[case::single_u_only("1#u", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Undetermined }, constraints: BondConstraintsForm::new() }))]
    #[case::single_mult("1#s2", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(2) }, constraints: BondConstraintsForm::new() }))]
    #[case::single_s_only("1#s", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(1) }, constraints: BondConstraintsForm::new() }))]
    #[case::double_charge_unpaired("2#c+#u2", BondDsl(BondForm { order: NumForm::Lit(2), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Undetermined }, constraints: BondConstraintsForm::new() }))]
    #[case::double_charge_mult("2#c-1#s3", BondDsl(BondForm { order: NumForm::Lit(2), charge: NumForm::Lit(-1), unpaired_electrons: UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(3) }, constraints: BondConstraintsForm::new() }))]
    #[case::aromatic("1#a", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]) }))]
    #[case::aromatic_plus("1#a+", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]) }))]
    #[case::aromatic_undetermined("1#a*", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Undetermined)]) }))]
    #[case::not_aromatic("1#a!", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(false))]) }))]
    #[case::charged_aromatic("1#c+#a", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]) }))]
    #[case::ring_membership_all("1#R2", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::All, NumForm::Lit(2))]) }))]
    #[case::ring_membership_all_bare("1#R", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::All, NumForm::Lit(1))]) }))]
    #[case::ring_membership_all_plus("1#R+", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::All, NumForm::RangeFrom(1))]) }))]
    #[case::ring_bang("1#R!", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::All, NumForm::Lit(0))]) }))]
    #[case::ring_zero("1#R0", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::All, NumForm::Lit(0))]) }))]
    #[case::ring_membership_all_star("1#R*", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::All, NumForm::Undetermined)]) }))]
    #[case::ring_membership_size("1#R(6)", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::Size(6), 1)]) }))]
    #[case::ring_membership_size_conj("1#R(5)#R(6)", BondDsl(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::from_iter([BondConstraintForm::ring_membership(RingScope::Size(5), 1), BondConstraintForm::ring_membership(RingScope::Size(6), 1)]) }))]
    #[case::whitespace_before_predicate("2 #c+", BondDsl(BondForm { order: NumForm::Lit(2), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() }))]
    #[case::whitespace_between_predicates("2#c+ #a", BondDsl(BondForm { order: NumForm::Lit(2), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]) }))]
    #[case::whitespace_surrounding_predicates("  2  #c+  #a  ", BondDsl(BondForm { order: NumForm::Lit(2), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]) }))]
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
    #[case::order_and_pred("1#a", BondUpdateDsl(BondUpdate { order: Some(NumForm::Lit(1)), constraints: BondConstraintsForm::from(BondConstraintForm::Aromatic(BooleanForm::Lit(true))), ..Default::default() }))]
    #[case::constraint_removal("#R(6)*", BondUpdateDsl(BondUpdate { constraints: BondConstraintsForm::from(BondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined)), ..Default::default() }))]
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
    #[case::aromatic(BondUpdateDsl(BondUpdate { constraints: BondConstraintsForm::from(BondConstraintForm::Aromatic(BooleanForm::Lit(true))), ..Default::default() }), r##""#a""##)]
    #[case::aromatic_undetermined(BondUpdateDsl(BondUpdate { constraints: BondConstraintsForm::from(BondConstraintForm::Aromatic(BooleanForm::Undetermined)), ..Default::default() }), r##""#a*""##)]
    #[case::ring_size_removal(BondUpdateDsl(BondUpdate { constraints: BondConstraintsForm::from(BondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined)), ..Default::default() }), r##""#R(6)*""##)]
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
    #[case::aromatic("#a", BondPredicate::Constraint(BondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    #[case::aromatic_plus("#a+", BondPredicate::Constraint(BondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    #[case::aromatic_false("#a!", BondPredicate::Constraint(BondConstraintForm::Aromatic(BooleanForm::Lit(false))))]
    #[case::aromatic_undetermined("#a*", BondPredicate::Constraint(BondConstraintForm::Aromatic(BooleanForm::Undetermined)))]
    #[case::ring_membership_all("#R2", BondPredicate::Constraint(BondConstraintForm::ring_membership(RingScope::All, NumForm::Lit(2))))]
    #[case::ring_membership_all_plus("#R+", BondPredicate::Constraint(BondConstraintForm::ring_membership(RingScope::All, NumForm::RangeFrom(1))))]
    #[case::ring_membership_zero("#R0", BondPredicate::Constraint(BondConstraintForm::ring_membership(RingScope::All, NumForm::Lit(0))))]
    #[case::ring_membership_all_undetermined("#R*", BondPredicate::Constraint(BondConstraintForm::ring_membership(RingScope::All, NumForm::Undetermined)))]
    #[case::ring_membership_size("#R(6)", BondPredicate::Constraint(BondConstraintForm::ring_membership(RingScope::Size(6), 1)))]
    #[case::ring_membership_size_undetermined("#R(6)*", BondPredicate::Constraint(BondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined)))]
    #[case::cis_trans_stereo_undetermined("#C*", BondPredicate::Constraint(BondConstraintForm::CisTransStereo(CisTransStereoForm::Undetermined)))]
    #[case::cis_trans_stereo_plus("#C+", BondPredicate::Constraint(BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(StereoCoset::Undetermined))))]
    #[case::cis_trans_stereo_not_stereo("#C!", BondPredicate::Constraint(BondConstraintForm::CisTransStereo(CisTransStereoForm::NotStereo)))]
    #[case::cis_trans_stereo("#C1", BondPredicate::Constraint(BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(StereoCoset::Lit(1)))))]
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
    fn test_bond_dsl_from_ir() {
        let mut form = BondForm::new(NumForm::Lit(1));
        form.charge = NumForm::Lit(0);
        form.unpaired_electrons = UnpairedElectronsForm::from((0_u8, 1_u8));
        let cfg = BondDefaults::zeroed();
        let dsl = BondDsl::from_ir(&form, &cfg);
        assert_eq!(dsl.0.charge, NumForm::Undetermined);
        assert_eq!(dsl.0.unpaired_electrons, UnpairedElectronsForm::default());
    }

    #[rstest]
    fn test_bond_dsl_into_ir() {
        let dsl = BondDsl(BondForm::new(NumForm::Lit(1)));
        let cfg = BondDefaults::zeroed();
        let form = dsl.into_ir(&cfg);
        assert_eq!(form.charge, NumForm::Lit(0));
        assert_eq!(
            form.unpaired_electrons,
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
    #[case::single(":single", BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() })]
    #[case::double(":double", BondForm { order: NumForm::Lit(2), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() })]
    #[case::triple(":triple", BondForm { order: NumForm::Lit(3), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() })]
    #[case::quadruple(":quadruple", BondForm { order: NumForm::Lit(4), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() })]
    #[case::aromatic(":aromatic", BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::from_iter([BondConstraintForm::Aromatic(BooleanForm::Lit(true))]) })]
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
    #[case::aromatic(BondConstraintForm::Aromatic(BooleanForm::Lit(true)), "{:aromatic true}")]
    #[case::aromatic_false(BondConstraintForm::Aromatic(BooleanForm::Lit(false)), "{:aromatic false}")]
    #[case::ring_membership_all(BondConstraintForm::ring_membership(RingScope::All, NumForm::Lit(1)), "{:ring-membership {:count 1}}")]
    #[case::ring_membership_all_undetermined(BondConstraintForm::ring_membership(RingScope::All, NumForm::Undetermined), "{:ring-membership {:count :undetermined}}")]
    #[case::ring_membership_size(BondConstraintForm::ring_membership(RingScope::Size(6), 1), "{:ring-membership {:size 6 :count 1}}")]
    #[case::ring_membership_size_count_set(BondConstraintForm::ring_membership(RingScope::Size(6), NumForm::lit_set([5, 6])), "{:ring-membership {:size 6 :count [5 6]}}")]
    #[case::cis_trans_stereo_undetermined(BondConstraintForm::CisTransStereo(CisTransStereoForm::Undetermined), "{:cis-trans-stereo :undetermined}")]
    #[case::cis_trans_stereo_not_stereo(BondConstraintForm::CisTransStereo(CisTransStereoForm::NotStereo), "{:cis-trans-stereo :not-stereo}")]
    #[case::cis_trans_stereo_lit(BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(StereoCoset::Lit(1))), "{:cis-trans-stereo {:stereo 1}}")]
    #[case::cis_trans_stereo_coset_undetermined(BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(StereoCoset::Undetermined)), "{:cis-trans-stereo {:stereo :undetermined}}")]
    #[case::cis_trans_stereo_set(BondConstraintForm::CisTransStereo(CisTransStereoForm::Stereo(StereoCoset::lit_set([1, 2]))), "{:cis-trans-stereo {:stereo [1 2]}}")]
    fn test_bond_constraint_dsl_roundtrip(
        #[case] input: BondConstraintForm,
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
        let form: BondForm = s.parse().unwrap();
        assert_eq!(form.to_string(), s);
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
