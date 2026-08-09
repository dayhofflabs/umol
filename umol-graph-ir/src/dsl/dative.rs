//! Dative-bond-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnStreamDeserializer, FromEdn, ToEdn};
use winnow::ascii::multispace0;
use winnow::combinator::{opt, preceded, repeat, terminated};
use winnow::error::ErrMode;
use winnow::token::take;
use winnow::Parser;

use super::boolean::{boolean, BooleanDsl};
use super::config::DativeBondDefaults;
use super::constraint::RingMembershipDsl;
use super::edn_utils::single_key_map;
use super::error::{PResult, ParseError};
use super::predicate::{fmt_ring_membership, ring_membership};
use super::value::{fmt_value, value};
use crate::ir::boolean::BooleanAst;
use crate::ir::constraint::{DativeBondConstraintAst, RingMembershipAst, RingScope};
use crate::ir::dative::{DativeBondAst, DativeBondUpdate};
use crate::ir::traits::{FromIr, IntoIr, Lattice};
use crate::ir::value::NumForm;

/// Surface DSL wrapper around `DativeBondAst`. The string form is the order
/// (number of donated electron pairs) followed by `#…` predicates,
/// paralleling `BondDsl`. Both `DativeBondConstraintAst` variants are inline:
/// `Aromatic` (`#a`) and `RingMembership` (`#R`, bare = total ring count,
/// `(s)` = count of size-`s` rings).
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DativeBondDsl(pub DativeBondAst);

impl DativeBondDsl {
    /// Zero-cost reference cast from `&DativeBondAst`. Relies on `repr(transparent)`.
    pub fn from_ref(ast: &DativeBondAst) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(ast as *const DativeBondAst as *const Self) }
    }
}

impl From<DativeBondAst> for DativeBondDsl {
    fn from(ast: DativeBondAst) -> Self {
        Self(ast)
    }
}

impl FromStr for DativeBondDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_dative_bond(s)
    }
}

impl Display for DativeBondDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_order(f, &self.0.order)?;
        for c in self.0.constraints.iter() {
            fmt_constraint(f, c)?;
        }
        Ok(())
    }
}

impl<'de> FromEdn<'de> for DativeBondDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s.parse().map_err(|e| DeError::subgrammar("dative", e)),
            Edn::Keyword(k) => {
                let s = expand_dative_keyword(k.name()).ok_or_else(|| {
                    DeError::Custom(format!("unknown dative keyword :{}", k.name()))
                })?;
                s.parse().map_err(|e| DeError::subgrammar("dative", e))
            }
            other => Err(DeError::TypeMismatch {
                expected: "string or dative-keyword",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        EdnStreamDeserializer::new(input).read_subgrammar_all("dative")
    }
}

/// Expand a dative-entry keyword shorthand to its equivalent dative-string
/// payload. Parallels [`super::bond::expand_bond_keyword`]:
///
/// - `:single` → `"1"`
/// - `:double` → `"2"`
/// - `:triple` → `"3"`
/// - `:quadruple` → `"4"`
///
/// Returns `None` for unrecognized keywords. Input sugar only — the AST
/// renders back to dative-string form.
pub(crate) fn expand_dative_keyword(name: &str) -> Option<&'static str> {
    match name {
        "single" => Some("1"),
        "double" => Some("2"),
        "triple" => Some("3"),
        "quadruple" => Some("4"),
        _ => None,
    }
}

impl ToEdn for DativeBondDsl {
    fn to_edn(&self) -> Edn<'static> {
        match dative_keyword_for(&self.0) {
            Some(kw) => Edn::Keyword(EdnKeyword::owned(kw.to_string())),
            None => Edn::Str(Cow::Owned(self.to_string())),
        }
    }
}

/// Return the dative-keyword shorthand for canonical dative shapes, or
/// `None` when the bond requires the full dative-string form. Inverse of
/// [`expand_dative_keyword`]: every shape this returns must round-trip.
///
/// Canonical means: no constraints and an integer order in 1..=4.
fn dative_keyword_for(ast: &DativeBondAst) -> Option<&'static str> {
    if !ast.constraints.is_empty() {
        return None;
    }
    match &ast.order {
        NumForm::Lit(1) => Some("single"),
        NumForm::Lit(2) => Some("double"),
        NumForm::Lit(3) => Some("triple"),
        NumForm::Lit(4) => Some("quadruple"),
        _ => None,
    }
}

impl FromIr<DativeBondAst> for DativeBondDsl {
    type Ctx = DativeBondDefaults;

    fn from_ir(ast: &DativeBondAst, _cfg: &Self::Ctx) -> Self {
        DativeBondDsl(ast.clone())
    }
}

impl IntoIr<DativeBondAst> for DativeBondDsl {
    type Ctx = DativeBondDefaults;

    fn into_ir(self, _cfg: &Self::Ctx) -> DativeBondAst {
        self.0
    }
}

/// Surface DSL wrapper around a [`DativeBondUpdate`].
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DativeBondUpdateDsl(pub DativeBondUpdate);

impl DativeBondUpdateDsl {
    /// Zero-cost reference cast from `&DativeBondUpdate`. Relies on `repr(transparent)`.
    pub fn from_ref(update: &DativeBondUpdate) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(update as *const DativeBondUpdate as *const Self) }
    }
}

impl FromIr<DativeBondUpdate> for DativeBondUpdateDsl {
    type Ctx = ();

    fn from_ir(update: &DativeBondUpdate, _ctx: &Self::Ctx) -> Self {
        Self(update.clone())
    }
}

impl IntoIr<DativeBondUpdate> for DativeBondUpdateDsl {
    type Ctx = ();

    fn into_ir(self, _ctx: &Self::Ctx) -> DativeBondUpdate {
        self.0
    }
}

impl FromStr for DativeBondUpdateDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_dative_bond_update(s)
    }
}

impl FromStr for DativeBondUpdate {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(DativeBondUpdateDsl::from_str(s)?.into_ir(&()))
    }
}

impl Display for DativeBondUpdate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        DativeBondUpdateDsl::from_ref(self).fmt(f)
    }
}

impl Display for DativeBondUpdateDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(order) = &self.0.order {
            fmt_order(f, order)?;
        }
        for c in self.0.constraints.iter() {
            if c.is_undetermined() {
                fmt_undetermined_constraint(f, c)?;
            } else {
                fmt_constraint(f, c)?;
            }
        }
        Ok(())
    }
}

impl<'de> FromEdn<'de> for DativeBondUpdateDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s
                .parse()
                .map_err(|e| DeError::subgrammar("dative-bond-update", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for DativeBondUpdateDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

impl FromStr for DativeBondAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(DativeBondDsl::from_str(s)?.into_ir(&DativeBondDefaults::default()))
    }
}

impl Display for DativeBondAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        DativeBondDsl::from_ref(self).fmt(f)
    }
}

impl<'de> FromEdn<'de> for DativeBondAst {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        Ok(DativeBondDsl::from_edn(edn)?.into_ir(&DativeBondDefaults::default()))
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        Ok(DativeBondDsl::from_edn_str(input)?.into_ir(&DativeBondDefaults::default()))
    }
}

impl ToEdn for DativeBondAst {
    fn to_edn(&self) -> Edn<'static> {
        DativeBondDsl::from_ref(self).to_edn()
    }
}

/// Parse a complete dative-bond-string into a `DativeBondDsl`.
pub fn parse_dative_bond(input: &str) -> Result<DativeBondDsl, ParseError> {
    dative_bond.parse(input).map_err(|e| e.into_inner())
}

/// Parse a complete dative-bond update string.
pub fn parse_dative_bond_update(input: &str) -> Result<DativeBondUpdateDsl, ParseError> {
    dative_bond_update.parse(input).map_err(|e| e.into_inner())
}

fn dative_bond_update(i: &mut &str) -> PResult<DativeBondUpdateDsl> {
    let order = preceded(multispace0, terminated(opt(value), multispace0)).parse_next(i)?;
    let preds: Vec<DativeBondPredicate> =
        repeat(0.., terminated(dative_bond_predicate, multispace0)).parse_next(i)?;
    let mut update = DativeBondUpdate {
        order,
        ..Default::default()
    };
    apply_update_predicates(&mut update, preds).map_err(ErrMode::Cut)?;
    Ok(DativeBondUpdateDsl(update))
}

pub(crate) fn dative_bond(i: &mut &str) -> PResult<DativeBondDsl> {
    let order = preceded(multispace0, terminated(value, multispace0)).parse_next(i)?;
    let preds: Vec<DativeBondPredicate> =
        repeat(0.., terminated(dative_bond_predicate, multispace0)).parse_next(i)?;
    let mut form = DativeBondDsl(DativeBondAst::new(order));
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

fn constraint_tag(c: &DativeBondConstraintAst) -> &'static str {
    match c {
        DativeBondConstraintAst::Aromatic(_) => "#a",
        DativeBondConstraintAst::RingMembership(..) => "#R",
    }
}

/// One predicate from a dative-bond-string; the parser yields a `Vec` of
/// these and the applier folds them into the `DativeBondAst`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DativeBondPredicate {
    Constraint(DativeBondConstraintAst),
}

fn dative_bond_predicate(i: &mut &str) -> PResult<DativeBondPredicate> {
    let start = *i;
    let prefix: &str = take(2usize).parse_next(i)?;
    match prefix {
        "#a" => boolean
            .map(|b| DativeBondPredicate::Constraint(DativeBondConstraintAst::Aromatic(b.0)))
            .parse_next(i),
        "#R" => ring_membership
            .map(|m| DativeBondPredicate::Constraint(DativeBondConstraintAst::RingMembership(m)))
            .parse_next(i),
        p if p.starts_with('#') => Err(ErrMode::Cut(ParseError::UnknownDativeBondPredicate(
            p.to_string(),
        ))),
        _ => Err(ErrMode::Cut(ParseError::TrailingInput(start.to_string()))),
    }
}

fn apply_predicates(
    form: &mut DativeBondDsl,
    preds: Vec<DativeBondPredicate>,
) -> Result<(), ParseError> {
    let ast = &mut form.0;
    for pred in preds {
        let DativeBondPredicate::Constraint(c) = pred;
        if ast.constraints.contains(c.key()) {
            return Err(ParseError::DuplicateDativeBondPredicate(
                constraint_tag(&c).to_string(),
            ));
        }
        ast.constraints.set(c);
    }
    Ok(())
}

fn apply_update_predicates(
    update: &mut DativeBondUpdate,
    preds: Vec<DativeBondPredicate>,
) -> Result<(), ParseError> {
    for pred in preds {
        let DativeBondPredicate::Constraint(constraint) = pred;
        if update.constraints.contains(constraint.key()) {
            return Err(ParseError::DuplicateDativeBondPredicate(
                constraint_tag(&constraint).to_string(),
            ));
        }
        update.constraints.set(constraint);
    }
    Ok(())
}

fn fmt_order(f: &mut fmt::Formatter<'_>, order: &NumForm) -> fmt::Result {
    match order {
        NumForm::Lit(n) => write!(f, "{}", n),
        NumForm::Undetermined => write!(f, "*"),
        v => fmt_value(f, v),
    }
}

fn fmt_constraint(f: &mut fmt::Formatter<'_>, c: &DativeBondConstraintAst) -> fmt::Result {
    match c {
        DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)) => write!(f, "#a"),
        DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false)) => write!(f, "#a!"),
        DativeBondConstraintAst::Aromatic(BooleanAst::Undetermined) => Ok(()),
        DativeBondConstraintAst::RingMembership(m) => fmt_ring_membership(f, m),
    }
}

fn fmt_undetermined_constraint(
    f: &mut fmt::Formatter<'_>,
    constraint: &DativeBondConstraintAst,
) -> fmt::Result {
    match constraint {
        DativeBondConstraintAst::RingMembership(membership) => match membership.scope {
            RingScope::All => write!(f, "#R*"),
            RingScope::Size(size) => write!(f, "#R({})*", size),
        },
        _ => write!(f, "{}*", constraint_tag(constraint)),
    }
}

/// Surface DSL wrapper around the narrow `DativeBondConstraintAst`. EDN form is a single-key map
/// `{:aromatic <bool>}` or `{:ring-membership {:size? <int> :count <value>}}`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum DativeBondConstraintDsl {
    Aromatic(BooleanAst),
    RingMembership(RingMembershipAst),
}

impl<'de> FromEdn<'de> for DativeBondConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Map(m) => {
                if m.len() != 1 {
                    return Err(DeError::Custom(format!(
                        "dative-bond-constraint must have exactly one key, got {}",
                        m.len()
                    )));
                }
                let (k, v) = m.iter().next().unwrap();
                let Edn::Keyword(key) = k else {
                    return Err(DeError::TypeMismatch {
                        expected: "keyword key",
                        got: k.kind(),
                        path: vec!["dative-bond-constraint".into()],
                    });
                };
                Ok(match key.name() {
                    "aromatic" => Self::Aromatic(BooleanDsl::from_edn(v)?.0),
                    "ring-membership" => Self::RingMembership(RingMembershipDsl::from_edn(v)?.0),
                    other => {
                        return Err(DeError::UnknownField {
                            key: other.to_string(),
                            path: vec!["dative-bond-constraint".into()],
                        });
                    }
                })
            }
            other => Err(DeError::TypeMismatch {
                expected: "{:aromatic …} / {:ring-membership …}",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for DativeBondConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        match self {
            Self::Aromatic(b) => single_key_map("aromatic", BooleanDsl(*b).to_edn()),
            Self::RingMembership(m) => {
                single_key_map("ring-membership", RingMembershipDsl(m.clone()).to_edn())
            }
        }
    }
}

impl DativeBondConstraintDsl {
    /// Build from the narrow inline AST form.
    pub(crate) fn from_ir(c: &DativeBondConstraintAst) -> Self {
        match c {
            DativeBondConstraintAst::Aromatic(b) => Self::Aromatic(*b),
            DativeBondConstraintAst::RingMembership(m) => Self::RingMembership(m.clone()),
        }
    }

    /// Convert into the narrow inline AST form.
    pub(crate) fn into_ir(self) -> DativeBondConstraintAst {
        match self {
            Self::Aromatic(b) => DativeBondConstraintAst::Aromatic(b),
            Self::RingMembership(m) => DativeBondConstraintAst::RingMembership(m),
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::read_string;

    use super::*;
    use crate::ir::constraint::{DativeBondConstraintsAst, RingScope};

    #[rustfmt::skip]
    #[rstest]
    #[case::single("1", DativeBondDsl(DativeBondAst::from_order(1)))]
    #[case::triple("3", DativeBondDsl(DativeBondAst::from_order(3)))]
    #[case::single_whitespace("  1  ", DativeBondDsl(DativeBondAst::from_order(1)))]
    #[case::undetermined_order("*", DativeBondDsl(DativeBondAst::default()))]
    #[case::aromatic("1#a", DativeBondDsl(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))) }))]
    #[case::aromatic_plus("1#a+", DativeBondDsl(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))) }))]
    #[case::aromatic_false("1#a!", DativeBondDsl(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false))) }))]
    #[case::aromatic_undetermined("1#a*", DativeBondDsl(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::Aromatic(BooleanAst::Undetermined)) }))]
    #[case::aromatic_with_ring("1#a#R(6)", DativeBondDsl(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1_i64)]) }))]
    #[case::ring_membership_all("1#R2", DativeBondDsl(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::All, NumForm::Lit(2))) }))]
    #[case::ring_membership_all_bare("1#R", DativeBondDsl(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::All, NumForm::Lit(1))) }))]
    #[case::ring_membership_all_plus("1#R+", DativeBondDsl(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::All, NumForm::RangeFrom(1))) }))]
    #[case::ring_membership_all_undetermined("1#R*", DativeBondDsl(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::All, NumForm::Undetermined)) }))]
    #[case::ring_membership_size("1#R(6)", DativeBondDsl(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1_i64)) }))]
    #[case::ring_membership_size_one("1#R(1)", DativeBondDsl(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::Size(1), 1_i64)) }))]
    #[case::ring_membership_all_and_size("1#R2#R(6)", DativeBondDsl(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from_iter([DativeBondConstraintAst::ring_membership(RingScope::All, NumForm::Lit(2)), DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1_i64)]) }))]
    #[case::triple_with_constraint("3#R+", DativeBondDsl(DativeBondAst { order: NumForm::Lit(3), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::All, NumForm::RangeFrom(1))) }))]
    fn test_parse_dative_bond(#[case] input: &str, #[case] expected: DativeBondDsl) {
        assert_eq!(parse_dative_bond(input).unwrap(), expected);
    }

    #[rstest]
    #[case::unknown("1#x", ParseError::UnknownDativeBondPredicate("#x".to_string()))]
    #[case::unknown_charge("1#c", ParseError::UnknownDativeBondPredicate("#c".to_string()))]
    #[case::duplicate_aromatic("1#a#a", ParseError::DuplicateDativeBondPredicate("#a".to_string()))]
    #[case::trailing("1#R2 foo", ParseError::TrailingInput("foo".to_string()))]
    fn test_parse_dative_bond_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(parse_dative_bond(input).unwrap_err(), expected);
    }

    #[rstest]
    #[case::single("1")]
    #[case::triple("3")]
    #[case::undetermined("*")]
    #[case::ring_membership_all("1#R2")]
    #[case::ring_membership_size("1#R(6)")]
    #[case::ring_membership_scopes("1#R2#R(6)")]
    #[case::aromatic("1#a")]
    #[case::aromatic_false("1#a!")]
    #[case::aromatic_with_ring("1#a#R(6)")]
    fn test_dative_bond_dsl_display_roundtrip(#[case] input: &str) {
        let dsl = parse_dative_bond(input).unwrap();
        assert_eq!(parse_dative_bond(&dsl.to_string()).unwrap(), dsl);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(
        DativeBondDsl(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::Aromatic(BooleanAst::Undetermined)) }),
        DativeBondDsl(DativeBondAst::from_order(1)),
    )]
    #[case::ring_membership_all(
        DativeBondDsl(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::All, NumForm::Undetermined)) }),
        DativeBondDsl(DativeBondAst::from_order(1)),
    )]
    #[case::ring_membership_size(
        DativeBondDsl(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::Size(6), NumForm::Undetermined)) }),
        DativeBondDsl(DativeBondAst::from_order(1)),
    )]
    fn test_dative_bond_dsl_display_vacuous_constraints(
        #[case] input: DativeBondDsl,
        #[case] expected: DativeBondDsl,
    ) {
        assert_eq!(parse_dative_bond(&input.to_string()).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::string(r##""1#a""##, DativeBondDsl(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))) }))]
    #[case::single_keyword(":single", DativeBondDsl(DativeBondAst::from_order(1)))]
    #[case::double_keyword(":double", DativeBondDsl(DativeBondAst::from_order(2)))]
    #[case::triple_keyword(":triple", DativeBondDsl(DativeBondAst::from_order(3)))]
    #[case::quadruple_keyword(":quadruple", DativeBondDsl(DativeBondAst::from_order(4)))]
    fn test_dative_bond_dsl_from_edn(#[case] input: &str, #[case] expected: DativeBondDsl) {
        assert_eq!(
            DativeBondDsl::from_edn(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rstest]
    #[case::unknown_keyword(":bogus", DeError::Custom("unknown dative keyword :bogus".to_string()))]
    #[case::wrong_type("3", DeError::TypeMismatch { expected: "string or dative-keyword", got: "int", path: Vec::new() })]
    fn test_dative_bond_dsl_from_edn_error(#[case] input: &str, #[case] expected: DeError) {
        assert_eq!(
            DativeBondDsl::from_edn(&read_string(input).unwrap()).unwrap_err(),
            expected,
        );
    }

    #[rstest]
    #[case::single(r##""1""##)]
    #[case::aromatic(r##""1#a""##)]
    #[case::ring_membership_all(r##""1#R2""##)]
    #[case::ring_membership_scopes(r##""1#R2#R(6)""##)]
    fn test_dative_bond_dsl_from_edn_parity(#[case] input: &str) {
        assert_eq!(
            DativeBondDsl::from_edn_str(input).unwrap(),
            DativeBondDsl::from_edn(&read_string(input).unwrap()).unwrap(),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(DativeBondDsl(DativeBondAst::from_order(1)), ":single")]
    #[case::double(DativeBondDsl(DativeBondAst::from_order(2)), ":double")]
    #[case::triple(DativeBondDsl(DativeBondAst::from_order(3)), ":triple")]
    #[case::quadruple(DativeBondDsl(DativeBondAst::from_order(4)), ":quadruple")]
    #[case::undetermined(DativeBondDsl(DativeBondAst::default()), r##""*""##)]
    #[case::aromatic(DativeBondDsl(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))) }), r##""1#a""##)]
    fn test_dative_bond_dsl_to_edn(#[case] input: DativeBondDsl, #[case] expected: &str) {
        assert_eq!(input.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::constrained(
        DativeBondDsl(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::All, NumForm::Lit(2))) }),
        DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::All, NumForm::Lit(2))) },
    )]
    fn test_dative_bond_dsl_into_ast(
        #[case] input: DativeBondDsl,
        #[case] expected: DativeBondAst,
    ) {
        assert_eq!(input.into_ir(&DativeBondDefaults::zeroed()), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", DativeBondUpdateDsl(DativeBondUpdate::default()))]
    #[case::order("2", DativeBondUpdateDsl(DativeBondUpdate { order: Some(NumForm::Lit(2)), ..Default::default() }))]
    #[case::order_undetermined("*", DativeBondUpdateDsl(DativeBondUpdate { order: Some(NumForm::Undetermined), ..Default::default() }))]
    #[case::aromatic("#a", DativeBondUpdateDsl(DativeBondUpdate { constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))), ..Default::default() }))]
    #[case::aromatic_undetermined("#a*", DativeBondUpdateDsl(DativeBondUpdate { constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::Aromatic(BooleanAst::Undetermined)), ..Default::default() }))]
    #[case::ring_size_removal("#R(6)*", DativeBondUpdateDsl(DativeBondUpdate { constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::Size(6), NumForm::Undetermined)), ..Default::default() }))]
    #[case::order_and_constraint("2#R", DativeBondUpdateDsl(DativeBondUpdate { order: Some(NumForm::Lit(2)), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::All, NumForm::Lit(1))) }))]
    fn test_parse_dative_bond_update(
        #[case] input: &str,
        #[case] expected: DativeBondUpdateDsl,
    ) {
        assert_eq!(parse_dative_bond_update(input).unwrap(), expected);
    }

    #[rstest]
    #[case::unknown("#x", ParseError::UnknownDativeBondPredicate("#x".to_string()))]
    #[case::duplicate("#a#a", ParseError::DuplicateDativeBondPredicate("#a".to_string()))]
    fn test_parse_dative_bond_update_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(parse_dative_bond_update(input).unwrap_err(), expected);
    }

    #[rstest]
    #[case::duplicate_aromatic("#a#a", ParseError::DuplicateDativeBondPredicate("#a".to_string()))]
    fn test_dative_bond_update_from_str_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(input.parse::<DativeBondUpdate>().unwrap_err(), expected);
    }

    #[rstest]
    #[case::string(r##""#R(6)*""##, DativeBondUpdateDsl(DativeBondUpdate { constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::Size(6), NumForm::Undetermined)), ..Default::default() }))]
    fn test_dative_bond_update_dsl_from_edn(
        #[case] input: &str,
        #[case] expected: DativeBondUpdateDsl,
    ) {
        assert_eq!(
            DativeBondUpdateDsl::from_edn(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rstest]
    #[case::wrong_type("3", DeError::TypeMismatch { expected: "string", got: "int", path: Vec::new() })]
    fn test_dative_bond_update_dsl_from_edn_error(#[case] input: &str, #[case] expected: DeError) {
        assert_eq!(
            DativeBondUpdateDsl::from_edn(&read_string(input).unwrap()).unwrap_err(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(DativeBondUpdateDsl(DativeBondUpdate::default()), r##""""##)]
    #[case::order(DativeBondUpdateDsl(DativeBondUpdate { order: Some(NumForm::Lit(2)), ..Default::default() }), r##""2""##)]
    #[case::order_undetermined(DativeBondUpdateDsl(DativeBondUpdate { order: Some(NumForm::Undetermined), ..Default::default() }), r##""*""##)]
    #[case::aromatic_removal(DativeBondUpdateDsl(DativeBondUpdate { constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::Aromatic(BooleanAst::Undetermined)), ..Default::default() }), r##""#a*""##)]
    #[case::ring_size_removal(DativeBondUpdateDsl(DativeBondUpdate { constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::Size(6), NumForm::Undetermined)), ..Default::default() }), r##""#R(6)*""##)]
    fn test_dative_bond_update_dsl_to_edn(
        #[case] update: DativeBondUpdateDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(update.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single("1", DativeBondAst::from_order(1))]
    #[case::triple("3", DativeBondAst::from_order(3))]
    #[case::aromatic("1#a", DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))) })]
    #[case::aromatic_false("1#a!", DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false))) })]
    #[case::ring_membership_all("1#R2", DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::All, 2_i64)) })]
    #[case::ring_membership_size("1#R(6)", DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1_i64)) })]
    fn test_dative_bond_ast_from_str(#[case] input: &str, #[case] expected: DativeBondAst) {
        assert_eq!(input.parse::<DativeBondAst>().unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(DativeBondAst::from_order(1), "1")]
    #[case::triple(DativeBondAst::from_order(3), "3")]
    #[case::aromatic(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))) }, "1#a")]
    #[case::aromatic_false(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false))) }, "1#a!")]
    #[case::ring_membership_all(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::All, 2_i64)) }, "1#R2")]
    #[case::ring_membership_size(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1_i64)) }, "1#R(6)")]
    fn test_dative_bond_ast_display(#[case] input: DativeBondAst, #[case] expected: &str) {
        assert_eq!(input.to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::string(r##""1#a""##, DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))) })]
    #[case::keyword(":double", DativeBondAst::from_order(2))]
    fn test_dative_bond_ast_from_edn(#[case] input: &str, #[case] expected: DativeBondAst) {
        assert_eq!(
            DativeBondAst::from_edn(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::single(DativeBondAst::from_order(1), ":single")]
    #[case::aromatic(DativeBondAst { order: NumForm::Lit(1), constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))) }, r##""1#a""##)]
    fn test_dative_bond_ast_to_edn(#[case] input: DativeBondAst, #[case] expected: &str) {
        assert_eq!(input.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic("{:aromatic true}", DativeBondConstraintDsl::Aromatic(BooleanAst::Lit(true)))]
    #[case::aromatic_false("{:aromatic false}", DativeBondConstraintDsl::Aromatic(BooleanAst::Lit(false)))]
    #[case::aromatic_undetermined("{:aromatic :undetermined}", DativeBondConstraintDsl::Aromatic(BooleanAst::Undetermined))]
    #[case::ring_membership_all("{:ring-membership {:count 2}}", DativeBondConstraintDsl::RingMembership(RingMembershipAst { scope: RingScope::All, count: NumForm::Lit(2) }))]
    #[case::ring_membership_size("{:ring-membership {:size 6 :count 1}}", DativeBondConstraintDsl::RingMembership(RingMembershipAst { scope: RingScope::Size(6), count: NumForm::Lit(1) }))]
    fn test_dative_bond_constraint_dsl_from_edn(
        #[case] input: &str,
        #[case] expected: DativeBondConstraintDsl,
    ) {
        assert_eq!(
            DativeBondConstraintDsl::from_edn(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rstest]
    #[case::wrong_shape(Edn::Int(3), DeError::TypeMismatch { expected: "{:aromatic …} / {:ring-membership …}", got: "int", path: vec![] })]
    #[case::unknown_key("{:bogus 1}", DeError::UnknownField { key: "bogus".to_string(), path: vec!["dative-bond-constraint".into()] })]
    fn test_dative_bond_constraint_dsl_from_edn_error(
        #[case] input: Edn<'static>,
        #[case] expected: DeError,
    ) {
        assert_eq!(
            DativeBondConstraintDsl::from_edn(&input).unwrap_err(),
            expected
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintDsl::Aromatic(BooleanAst::Lit(true)), "{:aromatic true}")]
    #[case::aromatic_false(DativeBondConstraintDsl::Aromatic(BooleanAst::Lit(false)), "{:aromatic false}")]
    #[case::aromatic_undetermined(DativeBondConstraintDsl::Aromatic(BooleanAst::Undetermined), "{:aromatic :undetermined}")]
    #[case::ring_membership_all(DativeBondConstraintDsl::RingMembership(RingMembershipAst { scope: RingScope::All, count: NumForm::Lit(2) }), "{:ring-membership {:count 2}}")]
    #[case::ring_membership_size(DativeBondConstraintDsl::RingMembership(RingMembershipAst { scope: RingScope::Size(6), count: NumForm::Lit(1) }), "{:ring-membership {:size 6 :count 1}}")]
    fn test_dative_bond_constraint_dsl_to_edn(
        #[case] input: DativeBondConstraintDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintDsl::Aromatic(BooleanAst::Lit(true)))]
    #[case::aromatic_false(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false)), DativeBondConstraintDsl::Aromatic(BooleanAst::Lit(false)))]
    #[case::aromatic_undetermined(DativeBondConstraintAst::Aromatic(BooleanAst::Undetermined), DativeBondConstraintDsl::Aromatic(BooleanAst::Undetermined))]
    #[case::ring_membership(DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1_i64), DativeBondConstraintDsl::RingMembership(RingMembershipAst { scope: RingScope::Size(6), count: NumForm::Lit(1) }))]
    fn test_dative_bond_constraint_dsl_from_ast(
        #[case] input: DativeBondConstraintAst,
        #[case] expected: DativeBondConstraintDsl,
    ) {
        assert_eq!(DativeBondConstraintDsl::from_ir(&input), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::aromatic(DativeBondConstraintDsl::Aromatic(BooleanAst::Lit(true)), DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)))]
    #[case::aromatic_false(DativeBondConstraintDsl::Aromatic(BooleanAst::Lit(false)), DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false)))]
    #[case::aromatic_undetermined(DativeBondConstraintDsl::Aromatic(BooleanAst::Undetermined), DativeBondConstraintAst::Aromatic(BooleanAst::Undetermined))]
    #[case::ring_membership(DativeBondConstraintDsl::RingMembership(RingMembershipAst { scope: RingScope::Size(6), count: NumForm::Lit(1) }), DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1_i64))]
    fn test_dative_bond_constraint_dsl_into_ast(
        #[case] input: DativeBondConstraintDsl,
        #[case] expected: DativeBondConstraintAst,
    ) {
        assert_eq!(input.into_ir(), expected);
    }
}
