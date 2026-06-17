//! Bond-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use strum::IntoEnumIterator;
use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn};
use winnow::ascii::multispace0;
use winnow::combinator::{delimited, repeat, terminated};
use winnow::error::ErrMode;
use winnow::token::take;
use winnow::Parser;

use super::atom::single_key_map;
use super::config::{BondDefaults, NumericDefault, StereoDefault};
use super::error::{PResult, ParseError};
use super::predicates::{
    apply_spin_pair, charge, fmt_charge, fmt_ring_count, fmt_spin_pair, lower_spin, optional_value,
    raise_spin, ring_count, SpinPredicate,
};
use super::stereo::{fmt_stereo_config, stereo_config, StereoConfigurationDsl};
use super::value::{fmt_value, value, ValueDsl};
use crate::ast::bond::BondAst;
use crate::ast::constraint::{BondConstraint, BondConstraintKind, BondConstraints};
use crate::ast::spin::SpinStateAst;
use crate::ast::stereo::{StereoConfigurationAst, StereoKind};
use crate::ast::traits::{FromAst, IntoAst};
use crate::ast::value::ValueAst;

/// Surface DSL wrapper around `BondAst`. Parses and renders the bond-string
/// form (order plus `#…` predicates); inline constraints land in
/// `self.0.constraints`.
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct BondDsl(pub BondAst);

impl BondDsl {
    /// Zero-cost reference cast from `&BondAst`. Relies on `repr(transparent)`.
    pub fn from_ref(ast: &BondAst) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(ast as *const BondAst as *const Self) }
    }
}

impl From<BondAst> for BondDsl {
    fn from(ast: BondAst) -> Self {
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
        fmt_bond_ast(f, &self.0)?;
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
/// payload. The five recognized keywords mirror the spec §7.7 table:
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
/// Canonical means: charge/spin at their defaults (Undetermined / default
/// pair) and either no constraints (orders 1–4) or exactly the `Aromatic`
/// flag (order 1 → `:aromatic`).
fn bond_keyword_for(ast: &BondAst) -> Option<&'static str> {
    if !matches!(ast.charge, ValueAst::Undetermined) || ast.spin != SpinStateAst::default() {
        return None;
    }
    let constraints: Vec<&BondConstraint> = ast.constraints.iter().collect();
    match (&ast.order, constraints.as_slice()) {
        (ValueAst::Lit(1), []) => Some("single"),
        (ValueAst::Lit(2), []) => Some("double"),
        (ValueAst::Lit(3), []) => Some("triple"),
        (ValueAst::Lit(4), []) => Some("quadruple"),
        (ValueAst::Lit(1), [BondConstraint::Aromatic]) => Some("aromatic"),
        _ => None,
    }
}

impl FromAst<BondAst> for BondDsl {
    type Ctx = BondDefaults;

    fn from_ast(ast: &BondAst, cfg: &Self::Ctx) -> Self {
        let mut out = ast.clone();
        lower_bond(&mut out, cfg);
        BondDsl(out)
    }
}

impl IntoAst<BondAst> for BondDsl {
    type Ctx = BondDefaults;

    fn into_ast(mut self, cfg: &Self::Ctx) -> BondAst {
        raise_bond(&mut self.0, cfg);
        self.0
    }
}

impl FromStr for BondAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(BondDsl::from_str(s)?.into_ast(&BondDefaults::default()))
    }
}

impl Display for BondAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        BondDsl::from_ref(self).fmt(f)
    }
}

impl<'de> FromEdn<'de> for BondAst {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        Ok(BondDsl::from_edn(edn)?.into_ast(&BondDefaults::default()))
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        Ok(BondDsl::from_edn_str(input)?.into_ast(&BondDefaults::default()))
    }
}

impl ToEdn for BondAst {
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
    let mut form = BondDsl(BondAst::new(order));
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

fn constraint_tag(c: &BondConstraint) -> &'static str {
    match c {
        BondConstraint::Aromatic => "#a",
        BondConstraint::RingCount(_) => "#R",
        BondConstraint::RingSize(_) => "#r",
        BondConstraint::CisTransStereo(_) => "#C",
    }
}

/// One predicate from a bond-string; the parser yields a `Vec` of these
/// and the applier folds them into the `BondAst`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BondPredicate {
    Charge(ValueAst),
    Spin(SpinPredicate),
    Constraint(BondConstraint),
}

fn bond_predicate(i: &mut &str) -> PResult<BondPredicate> {
    let start = *i;
    let prefix: &str = take(2usize).parse_next(i)?;
    match prefix {
        "#c" => charge.map(BondPredicate::Charge).parse_next(i),
        "#u" => optional_value
            .map(|v| BondPredicate::Spin(SpinPredicate::Unpaired(v)))
            .parse_next(i),
        "#s" => optional_value
            .map(|v| BondPredicate::Spin(SpinPredicate::Multiplicity(v)))
            .parse_next(i),
        "#a" => Ok(BondPredicate::Constraint(BondConstraint::Aromatic)),
        "#R" => ring_count
            .map(|v| BondPredicate::Constraint(BondConstraint::RingCount(v)))
            .parse_next(i),
        "#r" => optional_value
            .map(|v| BondPredicate::Constraint(BondConstraint::RingSize(v)))
            .parse_next(i),
        "#C" => (|i: &mut &str| stereo_config(i, StereoKind::CisTrans.degree()))
            .map(|c| BondPredicate::Constraint(BondConstraint::CisTransStereo(c)))
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
                if !matches!(ast.charge, ValueAst::Undetermined) {
                    return Err(ParseError::DuplicateBondPredicate("#c".to_string()));
                }
                ast.charge = v;
            }
            BondPredicate::Spin(sp) => {
                apply_spin_pair(&mut ast.spin, sp, ParseError::DuplicateBondPredicate)?;
            }
            BondPredicate::Constraint(c) => {
                let tag = constraint_tag(&c);
                if c.is_unique()
                    && ast
                        .constraints
                        .iter()
                        .any(|existing| constraint_tag(existing) == tag)
                {
                    return Err(ParseError::DuplicateBondPredicate(tag.to_string()));
                }
                ast.constraints.add(c);
            }
        }
    }
    Ok(())
}

fn fmt_bond_ast(f: &mut fmt::Formatter<'_>, ast: &BondAst) -> fmt::Result {
    match &ast.order {
        ValueAst::Lit(n) => write!(f, "{}", n)?,
        ValueAst::Undetermined => write!(f, "*")?,
        v => fmt_value(f, v)?,
    }

    fmt_charge(f, &ast.charge)?;
    fmt_spin_pair(f, &ast.spin)
}

fn fmt_constraint(f: &mut fmt::Formatter<'_>, c: &BondConstraint) -> fmt::Result {
    match c {
        BondConstraint::Aromatic => write!(f, "#a"),
        BondConstraint::RingCount(v) => fmt_ring_count(f, v),
        BondConstraint::RingSize(v) => match v {
            ValueAst::Undetermined => Ok(()),
            ValueAst::Lit(1) => write!(f, "#r"),
            ValueAst::Lit(n) => write!(f, "#r{}", n),
            v => {
                write!(f, "#r")?;
                fmt_value(f, v)
            }
        },
        BondConstraint::CisTransStereo(c) => {
            write!(f, "#C")?;
            fmt_stereo_config(f, c)
        }
    }
}

fn lower_bond(ast: &mut BondAst, cfg: &BondDefaults) {
    // Exhaustive destructure: adding a new BondAst field is a compile error
    // here, forcing the author to decide how lowering should handle it.
    let BondAst {
        order: _,
        charge,
        spin,
        constraints,
    } = ast;

    if matches!(
        (&cfg.charge, &*charge),
        (NumericDefault::Zero, ValueAst::Lit(0))
    ) {
        *charge = ValueAst::Undetermined;
    }
    lower_spin(spin, cfg.unpaired_electrons, cfg.multiplicity);
    lower_bond_constraints(constraints, cfg);
}

fn raise_bond(ast: &mut BondAst, cfg: &BondDefaults) {
    // Exhaustive destructure: adding a new BondAst field is a compile error
    // here, forcing the author to decide how raising should handle it.
    let BondAst {
        order: _,
        charge,
        spin,
        constraints,
    } = ast;

    if matches!(*charge, ValueAst::Undetermined) {
        *charge = match cfg.charge {
            NumericDefault::Zero => ValueAst::Lit(0),
            NumericDefault::Required => ValueAst::Undetermined,
        };
    }
    raise_spin(spin, cfg.unpaired_electrons, cfg.multiplicity);
    raise_bond_constraints(constraints, cfg);
}

fn raise_bond_constraints(constraints: &mut BondConstraints, cfg: &BondDefaults) {
    for kind in BondConstraintKind::iter() {
        match kind {
            BondConstraintKind::CisTransStereo => {
                if !constraints.contains(kind) {
                    match cfg.cis_trans_stereo {
                        StereoDefault::NotStereo => {
                            constraints.add(BondConstraint::CisTransStereo(
                                StereoConfigurationAst::NotStereo,
                            ));
                        }
                        StereoDefault::Required => {}
                    }
                }
            }
            BondConstraintKind::Aromatic
            | BondConstraintKind::RingCount
            | BondConstraintKind::RingSize => {
                // Pattern-only constraint: no defaulting mode in BondDefaults.
            }
        }
    }
}

fn lower_bond_constraints(constraints: &mut BondConstraints, cfg: &BondDefaults) {
    for kind in BondConstraintKind::iter() {
        match kind {
            BondConstraintKind::CisTransStereo => match cfg.cis_trans_stereo {
                StereoDefault::NotStereo => {
                    if matches!(
                        constraints.get(kind),
                        Some(BondConstraint::CisTransStereo(
                            StereoConfigurationAst::NotStereo
                        ))
                    ) {
                        constraints.remove(kind);
                    }
                }
                StereoDefault::Required => {}
            },
            BondConstraintKind::Aromatic
            | BondConstraintKind::RingCount
            | BondConstraintKind::RingSize => {
                // Pattern-only constraint: no defaulting mode in BondDefaults.
            }
        }
    }
}

/// Surface DSL wrapper around `BondConstraint`. EDN form: the keyword
/// `:aromatic` (flag variant, no value) or a single-key map
/// `{:ring-count <value>}` / `{:ring-size <value>}`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct BondConstraintDsl(pub BondConstraint);

impl FromAst<BondConstraint> for BondConstraintDsl {
    type Ctx = ();

    fn from_ast(ast: &BondConstraint, _ctx: &Self::Ctx) -> Self {
        Self(ast.clone())
    }
}

impl IntoAst<BondConstraint> for BondConstraintDsl {
    type Ctx = ();

    fn into_ast(self, _ctx: &Self::Ctx) -> BondConstraint {
        self.0
    }
}

impl<'de> FromEdn<'de> for BondConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Keyword(k) if k.name() == "aromatic" => Ok(Self(BondConstraint::Aromatic)),
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
                    "ring-count" => BondConstraint::RingCount(ValueDsl::from_edn(v)?.into_ast(&())),
                    "ring-size" => BondConstraint::RingSize(ValueDsl::from_edn(v)?.into_ast(&())),
                    "cis-trans-stereo" => BondConstraint::CisTransStereo(
                        StereoConfigurationDsl::from_edn(v)?.into_ast(&()),
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
                expected: ":aromatic / {:ring-count <value>} / {:ring-size <value>}",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for BondConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        match &self.0 {
            BondConstraint::Aromatic => Edn::Keyword(EdnKeyword::owned("aromatic".into())),
            BondConstraint::RingCount(v) => bond_constraint_single_key("ring-count", v),
            BondConstraint::RingSize(v) => bond_constraint_single_key("ring-size", v),
            BondConstraint::CisTransStereo(c) => single_key_map(
                "cis-trans-stereo",
                StereoConfigurationDsl::from_ast(c, &()).to_edn(),
            ),
        }
    }
}

fn bond_constraint_single_key(key: &str, v: &ValueAst) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(1);
    m.insert(
        Edn::Keyword(EdnKeyword::owned(key.into())),
        ValueDsl::from_ast(v, &()).to_edn(),
    );
    Edn::Map(m)
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::constraint::BondConstraints;
    use std::collections::BTreeSet;

    use crate::ast::operators::MemOp;
    use crate::ast::spin::SpinStateAst;
    use crate::ast::stereo::{StereoCosetAst, StereoExpr};
    use crate::ast::value::{ValuePredicate, ValueTerm};

    #[rustfmt::skip]
    #[rstest]
    #[case::single("1", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::new() }))]
    #[case::double("2", BondDsl(BondAst { order: ValueAst::Lit(2), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::new() }))]
    #[case::triple("3", BondDsl(BondAst { order: ValueAst::Lit(3), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::new() }))]
    #[case::quadruple("4", BondDsl(BondAst { order: ValueAst::Lit(4), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::new() }))]
    #[case::single_whitespace("  1  ", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::new() }))]
    #[case::single_pos_charge("1#c+2", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(2), spin: SpinStateAst::default(), constraints: BondConstraints::new() }))]
    #[case::single_neg_charge("1#c-2", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(-2), spin: SpinStateAst::default(), constraints: BondConstraints::new() }))]
    #[case::single_zero_charge("1#c0", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(0), spin: SpinStateAst::default(), constraints: BondConstraints::new() }))]
    #[case::single_plus_only("1#c+", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: BondConstraints::new() }))]
    #[case::single_minus_only("1#c-", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(-1), spin: SpinStateAst::default(), constraints: BondConstraints::new() }))]
    #[case::double_unpaired("2#u3", BondDsl(BondAst { order: ValueAst::Lit(2), charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Lit(3), multiplicity: ValueAst::Undetermined }, constraints: BondConstraints::new() }))]
    #[case::single_u_only("1#u", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Undetermined }, constraints: BondConstraints::new() }))]
    #[case::single_mult("1#s2", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(2) }, constraints: BondConstraints::new() }))]
    #[case::single_s_only("1#s", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(1) }, constraints: BondConstraints::new() }))]
    #[case::double_charge_unpaired("2#c+#u2", BondDsl(BondAst { order: ValueAst::Lit(2), charge: ValueAst::Lit(1), spin: SpinStateAst { unpaired: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined }, constraints: BondConstraints::new() }))]
    #[case::double_charge_mult("2#c-1#s3", BondDsl(BondAst { order: ValueAst::Lit(2), charge: ValueAst::Lit(-1), spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(3) }, constraints: BondConstraints::new() }))]
    #[case::aromatic("1#a", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::from_iter([BondConstraint::Aromatic]) }))]
    #[case::charged_aromatic("1#c+#a", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: BondConstraints::from_iter([BondConstraint::Aromatic]) }))]
    #[case::ring_count("1#R2", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::from_iter([BondConstraint::RingCount(ValueAst::Lit(2))]) }))]
    #[case::ring_bare("1#R", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::from_iter([BondConstraint::RingCount(ValueAst::Lit(1))]) }))]
    #[case::ring_plus("1#R+", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::from_iter([BondConstraint::RingCount(ValueAst::var_at_least("r", 1))]) }))]
    #[case::ring_bang("1#R!", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::from_iter([BondConstraint::RingCount(ValueAst::Lit(0))]) }))]
    #[case::ring_zero("1#R0", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::from_iter([BondConstraint::RingCount(ValueAst::Lit(0))]) }))]
    #[case::ring_undetermined("1#R*", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::from_iter([BondConstraint::RingCount(ValueAst::Undetermined)]) }))]
    #[case::ring_size("1#r6", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::from_iter([BondConstraint::RingSize(ValueAst::Lit(6))]) }))]
    #[case::ring_size_conj("1#r5#r6", BondDsl(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::from_iter([BondConstraint::RingSize(ValueAst::Lit(5)), BondConstraint::RingSize(ValueAst::Lit(6))]) }))]
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
    #[case::charge_pos("#c+2", BondPredicate::Charge(ValueAst::Lit(2)))]
    #[case::charge_neg("#c-2", BondPredicate::Charge(ValueAst::Lit(-2)))]
    #[case::charge_plus("#c+", BondPredicate::Charge(ValueAst::Lit(1)))]
    #[case::charge_minus("#c-", BondPredicate::Charge(ValueAst::Lit(-1)))]
    #[case::charge_zero("#c0", BondPredicate::Charge(ValueAst::Lit(0)))]
    #[case::charge_undetermined("#c*", BondPredicate::Charge(ValueAst::Undetermined))]
    #[case::unpaired("#u2", BondPredicate::Spin(SpinPredicate::Unpaired(ValueAst::Lit(2))))]
    #[case::unpaired_omit("#u", BondPredicate::Spin(SpinPredicate::Unpaired(ValueAst::Lit(1))))]
    #[case::unpaired_undetermined("#u*", BondPredicate::Spin(SpinPredicate::Unpaired(ValueAst::Undetermined)))]
    #[case::multiplicity("#s3", BondPredicate::Spin(SpinPredicate::Multiplicity(ValueAst::Lit(3))))]
    #[case::multiplicity_omit("#s", BondPredicate::Spin(SpinPredicate::Multiplicity(ValueAst::Lit(1))))]
    #[case::multiplicity_undetermined("#s*", BondPredicate::Spin(SpinPredicate::Multiplicity(ValueAst::Undetermined)))]
    #[case::aromatic("#a", BondPredicate::Constraint(BondConstraint::Aromatic))]
    fn test_bond_predicate(#[case] input: &str, #[case] expected: BondPredicate) {
        let result = bond_predicate.parse(input);
        assert!(result.is_ok(), "{input:?} should succeed, got {:?}", result.unwrap_err());
        let pred = result.unwrap();
        assert_eq!(pred, expected);
    }

    #[rstest]
    #[case::unknown("#x", ParseError::UnknownBondPredicate("#x".to_string()))]
    #[case::unknown_tag("#z", ParseError::UnknownBondPredicate("#z".to_string()))]
    #[case::trailing_no_hash("fo", ParseError::TrailingInput("fo".to_string()))]
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
        let mut ast = BondAst::new(ValueAst::Lit(1));
        ast.charge = ValueAst::Lit(0);
        ast.spin = SpinStateAst::from((0_u8, 1_u8));
        let cfg = BondDefaults::zeroed();
        let dsl = BondDsl::from_ast(&ast, &cfg);
        assert_eq!(dsl.0.charge, ValueAst::Undetermined);
        assert_eq!(dsl.0.spin, SpinStateAst::default());
    }

    #[rstest]
    fn test_bond_dsl_into_ast() {
        let dsl = BondDsl(BondAst::new(ValueAst::Lit(1)));
        let cfg = BondDefaults::zeroed();
        let ast = dsl.into_ast(&cfg);
        assert_eq!(ast.charge, ValueAst::Lit(0));
        assert_eq!(ast.spin, SpinStateAst::from((0_u8, 1_u8)));
    }

    #[rstest]
    fn test_bond_dsl_roundtrip_zeroed() {
        let input = BondDsl(BondAst::new(ValueAst::Lit(2)));
        let cfg = BondDefaults::zeroed();
        let raised = input.clone().into_ast(&cfg);
        let lowered = BondDsl::from_ast(&raised, &cfg);
        assert_eq!(input, lowered);
    }

    #[rstest]
    #[case::single(r##""1""##)]
    #[case::aromatic(r##""1#a""##)]
    #[case::ring_count(r##""2#R+""##)]
    #[case::with_cis_trans_stereo(r##""2#C1""##)]
    fn test_bond_dsl_from_edn_str_matches_from_edn(#[case] input: &str) {
        let via_stream = BondDsl::from_edn_str(input).unwrap();
        let tree = umol_edn::read_string(input).unwrap();
        let via_tree = BondDsl::from_edn(&tree).unwrap();
        assert_eq!(via_stream, via_tree);
    }

    /// Bond-keyword shorthands expand to their bond-string equivalents
    /// in the tree-parse path (`BondDsl::from_edn`). Mirrors the spec
    /// §7.7 expansion table.
    #[rustfmt::skip]
    #[rstest]
    #[case::single(":single", BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::new() })]
    #[case::double(":double", BondAst { order: ValueAst::Lit(2), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::new() })]
    #[case::triple(":triple", BondAst { order: ValueAst::Lit(3), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::new() })]
    #[case::quadruple(":quadruple", BondAst { order: ValueAst::Lit(4), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::new() })]
    #[case::aromatic(":aromatic", BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::from_iter([BondConstraint::Aromatic]) })]
    fn test_bond_dsl_keyword_shorthand_tree(#[case] input: &str, #[case] expected: BondAst) {
        let edn = umol_edn::read_string(input).unwrap();
        let dsl = BondDsl::from_edn(&edn).unwrap();
        assert_eq!(dsl.0, expected);
    }

    #[rstest]
    fn test_bond_dsl_keyword_shorthand_unknown_rejected() {
        let edn = umol_edn::read_string(":bogus").unwrap();
        let err = BondDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(BondConstraint::Aromatic, ":aromatic")]
    #[case::ring_count(BondConstraint::RingCount(ValueAst::Lit(1)), "{:ring-count 1}")]
    #[case::ring_count_undetermined(BondConstraint::RingCount(ValueAst::Undetermined), "{:ring-count :undetermined}")]
    #[case::ring_size(BondConstraint::RingSize(ValueAst::Lit(6)), "{:ring-size 6}")]
    #[case::ring_size_set(BondConstraint::RingSize(ValueAst::lit_set([5, 6])), "{:ring-size [5 6]}")]
    #[case::cis_trans_stereo_undetermined(BondConstraint::CisTransStereo(StereoConfigurationAst::Undetermined), "{:cis-trans-stereo :undetermined}")]
    #[case::cis_trans_stereo_not_stereo(BondConstraint::CisTransStereo(StereoConfigurationAst::NotStereo), "{:cis-trans-stereo :not-stereo}")]
    #[case::cis_trans_stereo_lit(BondConstraint::CisTransStereo(StereoConfigurationAst::Stereo(StereoCosetAst::Lit(1))), "{:cis-trans-stereo {:stereo 1}}")]
    #[case::cis_trans_stereo_coset_undetermined(BondConstraint::CisTransStereo(StereoConfigurationAst::Stereo(StereoCosetAst::Undetermined)), "{:cis-trans-stereo {:stereo :undetermined}}")]
    #[case::cis_trans_stereo_set(BondConstraint::CisTransStereo(StereoConfigurationAst::Stereo(StereoCosetAst::Expr(Box::new(StereoExpr::LitSet(vec![1, 2]))))), "{:cis-trans-stereo {:stereo [1 2]}}")]
    fn test_bond_constraint_dsl_roundtrip(
        #[case] input: BondConstraint,
        #[case] edn_source: &str,
    ) {
        let dsl = BondConstraintDsl::from_ast(&input, &());
        let edn = dsl.to_edn();
        let expected = umol_edn::read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = BondConstraintDsl::from_edn(&edn).unwrap();
        assert_eq!(parsed.into_ast(&()), input, "parse-back mismatch");
    }

    /// Vacuous bond constraints elide on rendering (canonical-rendering
    /// rule, see `dsl::predicates`). `#R*` and `#r*` are admitted on
    /// parse but the rendered surface drops them.
    #[rstest]
    #[case::ring_count("1#R*", "1")]
    #[case::ring_size("1#r*", "1")]
    fn test_bond_render_elides_vacuous_constraints(
        #[case] input: &str,
        #[case] expected_canonical: &str,
    ) {
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
    fn test_bond_constraint_dsl_rejects_wrong_shape() {
        let err = BondConstraintDsl::from_edn(&Edn::Int(3)).unwrap_err();
        assert!(matches!(err, DeError::TypeMismatch { .. }));
    }

    #[rstest]
    fn test_bond_constraint_dsl_rejects_unknown_key() {
        let edn = umol_edn::read_string("{:bogus 1}").unwrap();
        let err = BondConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::UnknownField { .. }));
    }

    #[rstest]
    fn test_bond_constraint_dsl_accepts_value_as_string_subgrammar() {
        let edn = umol_edn::read_string(r##"{:ring-size "?n :: {5,6}"}"##).unwrap();
        let parsed = BondConstraintDsl::from_edn(&edn).unwrap();
        assert_eq!(
            parsed.into_ast(&()),
            BondConstraint::RingSize(ValueAst::predicate(ValuePredicate::Mem(
                ValueTerm::Var("n".to_string()),
                MemOp::In,
                BTreeSet::from([5, 6]),
            ))),
        );
    }

    #[rstest]
    #[case::single("1")]
    #[case::double("2")]
    #[case::aromatic("1#a")]
    fn test_bond_ast_from_str_to_string_roundtrip(#[case] s: &str) {
        let ast: BondAst = s.parse().unwrap();
        assert_eq!(ast.to_string(), s);
    }

    #[rstest]
    fn test_bond_ast_to_edn_roundtrip() {
        use umol_edn::ToEdn;
        let ast: BondAst = "2".parse().unwrap();
        let edn = ast.to_edn();
        let back = BondAst::from_edn(&edn).unwrap();
        assert_eq!(back, ast);
    }
}
