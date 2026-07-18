//! Noncovalent-bond-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, EdnError, EdnStreamDeserializer, FromEdn, ToEdn};
use winnow::ascii::multispace0;
use winnow::combinator::{alt, preceded, repeat, terminated};
use winnow::error::{ErrMode, ParserError};
use winnow::token::{one_of, take};
use winnow::Parser;

use super::boolean::{boolean, BooleanDsl};
use super::config::NoncovalentBondDefaults;
use super::edn_utils::single_key_map;
use super::error::{PResult, ParseError};
use crate::ast::boolean::BooleanAst;
use crate::ast::constraint::NoncovalentBondConstraintAst;
use crate::ast::noncovalent::{
    NoncovalentBondAst, NoncovalentBondKind, NoncovalentBondKindAst, NoncovalentBondUpdate,
};
use crate::ast::traits::{FromAst, IntoAst, Lattice};

/// Surface DSL wrapper around `NoncovalentBondAst`.
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NoncovalentBondDsl(pub NoncovalentBondAst);

impl NoncovalentBondDsl {
    /// Zero-cost reference cast from `&NoncovalentBondAst`. Relies on `repr(transparent)`.
    pub fn from_ref(ast: &NoncovalentBondAst) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(ast as *const NoncovalentBondAst as *const Self) }
    }
}

impl From<NoncovalentBondAst> for NoncovalentBondDsl {
    fn from(ast: NoncovalentBondAst) -> Self {
        Self(ast)
    }
}

impl FromStr for NoncovalentBondDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_noncovalent_bond(s)
    }
}

impl Display for NoncovalentBondDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_noncovalent_bond_ast(f, &self.0)
    }
}

impl<'de> FromEdn<'de> for NoncovalentBondDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s.parse().map_err(|e| DeError::subgrammar("noncovalent", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        EdnStreamDeserializer::new(input).read_subgrammar_all("noncovalent")
    }
}

impl ToEdn for NoncovalentBondDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

impl FromAst<NoncovalentBondAst> for NoncovalentBondDsl {
    type Ctx = NoncovalentBondDefaults;

    fn from_ast(ast: &NoncovalentBondAst, _ctx: &Self::Ctx) -> Self {
        NoncovalentBondDsl(ast.clone())
    }
}

impl IntoAst<NoncovalentBondAst> for NoncovalentBondDsl {
    type Ctx = NoncovalentBondDefaults;

    fn into_ast(self, _ctx: &Self::Ctx) -> NoncovalentBondAst {
        self.0
    }
}

impl FromStr for NoncovalentBondAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(NoncovalentBondDsl::from_str(s)?.into_ast(&NoncovalentBondDefaults::default()))
    }
}

impl Display for NoncovalentBondAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        NoncovalentBondDsl::from_ref(self).fmt(f)
    }
}

impl<'de> FromEdn<'de> for NoncovalentBondAst {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        Ok(NoncovalentBondDsl::from_edn(edn)?.into_ast(&NoncovalentBondDefaults::default()))
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        Ok(NoncovalentBondDsl::from_edn_str(input)?.into_ast(&NoncovalentBondDefaults::default()))
    }
}

impl ToEdn for NoncovalentBondAst {
    fn to_edn(&self) -> Edn<'static> {
        NoncovalentBondDsl::from_ref(self).to_edn()
    }
}

/// Parse a complete noncovalent-bond-string into a `NoncovalentBondDsl`.
pub fn parse_noncovalent_bond(input: &str) -> Result<NoncovalentBondDsl, ParseError> {
    noncovalent_bond.parse(input).map_err(|e| e.into_inner())
}

/// Parse a complete noncovalent-bond update string.
pub fn parse_noncovalent_bond_update(input: &str) -> Result<NoncovalentBondUpdateDsl, ParseError> {
    noncovalent_bond_update
        .parse(input)
        .map_err(|e| e.into_inner())
}

pub(crate) fn noncovalent_bond(i: &mut &str) -> PResult<NoncovalentBondDsl> {
    let kind = preceded(multispace0, terminated(kind_expr, multispace0)).parse_next(i)?;
    let preds: Vec<NoncovalentBondPredicate> =
        repeat(0.., terminated(noncovalent_bond_predicate, multispace0)).parse_next(i)?;
    let mut form = NoncovalentBondDsl(NoncovalentBondAst::new(kind));
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

/// One predicate from a noncovalent-bond-string; the parser yields a `Vec` of
/// these and the applier folds them into the `NoncovalentBondAst`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NoncovalentBondPredicate {
    Constraint(NoncovalentBondConstraintAst),
}

fn noncovalent_bond_predicate(i: &mut &str) -> PResult<NoncovalentBondPredicate> {
    let start = *i;
    let prefix: &str = take(2usize).parse_next(i)?;
    match prefix {
        "#I" => boolean
            .map(|b| {
                NoncovalentBondPredicate::Constraint(NoncovalentBondConstraintAst::Intramolecular(
                    b.0,
                ))
            })
            .parse_next(i),
        p if p.starts_with('#') => Err(ErrMode::Cut(ParseError::UnknownNoncovalentBondPredicate(
            p.to_string(),
        ))),
        _ => Err(ErrMode::Cut(ParseError::TrailingInput(start.to_string()))),
    }
}

fn apply_predicates(
    form: &mut NoncovalentBondDsl,
    preds: Vec<NoncovalentBondPredicate>,
) -> Result<(), ParseError> {
    let ast = &mut form.0;
    for pred in preds {
        let NoncovalentBondPredicate::Constraint(c) = pred;
        if ast.constraints.contains(c.key()) {
            return Err(ParseError::DuplicateNoncovalentBondPredicate(
                constraint_tag(&c).to_string(),
            ));
        }
        ast.constraints.set(c);
    }
    Ok(())
}

fn constraint_tag(c: &NoncovalentBondConstraintAst) -> &'static str {
    match c {
        NoncovalentBondConstraintAst::Intramolecular(_) => "#I",
    }
}

fn kind_expr(i: &mut &str) -> PResult<NoncovalentBondKindAst> {
    alt((
        '*'.value(NoncovalentBondKindAst::Undetermined),
        kind_literal.map(NoncovalentBondKindAst::Lit),
    ))
    .parse_next(i)
    .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedNoncovalentBondKind))
}

fn kind_literal(i: &mut &str) -> PResult<NoncovalentBondKind> {
    let sym: &str = (
        one_of(|c: char| c.is_ascii_uppercase()),
        one_of(|c: char| c.is_ascii_lowercase()),
        one_of(|c: char| c.is_ascii_lowercase()),
    )
        .take()
        .parse_next(i)?;
    match kind_from_symbol(sym) {
        Some(k) => Ok(k),
        None => Err(ErrMode::Backtrack(ParseError::from_input(i))),
    }
}

fn kind_from_symbol(sym: &str) -> Option<NoncovalentBondKind> {
    match sym {
        "Hbd" => Some(NoncovalentBondKind::HydrogenBond),
        "Xbd" => Some(NoncovalentBondKind::HalogenBond),
        "Ybd" => Some(NoncovalentBondKind::ChalcogenBond),
        "Ion" => Some(NoncovalentBondKind::Ionic),
        "Vdw" => Some(NoncovalentBondKind::VanDerWaals),
        _ => None,
    }
}

fn kind_symbol(k: NoncovalentBondKind) -> &'static str {
    match k {
        NoncovalentBondKind::HydrogenBond => "Hbd",
        NoncovalentBondKind::HalogenBond => "Xbd",
        NoncovalentBondKind::ChalcogenBond => "Ybd",
        NoncovalentBondKind::Ionic => "Ion",
        NoncovalentBondKind::VanDerWaals => "Vdw",
    }
}

fn fmt_noncovalent_bond_ast(f: &mut fmt::Formatter<'_>, ast: &NoncovalentBondAst) -> fmt::Result {
    fmt_kind(f, &ast.kind)?;
    for c in ast.constraints.iter() {
        fmt_constraint(f, c)?;
    }
    Ok(())
}

fn fmt_constraint(f: &mut fmt::Formatter<'_>, c: &NoncovalentBondConstraintAst) -> fmt::Result {
    match c {
        NoncovalentBondConstraintAst::Intramolecular(BooleanAst::Lit(true)) => write!(f, "#I"),
        NoncovalentBondConstraintAst::Intramolecular(BooleanAst::Lit(false)) => write!(f, "#I!"),
        NoncovalentBondConstraintAst::Intramolecular(BooleanAst::Undetermined) => Ok(()),
    }
}

fn fmt_kind(f: &mut fmt::Formatter<'_>, kind: &NoncovalentBondKindAst) -> fmt::Result {
    match kind {
        NoncovalentBondKindAst::Lit(k) => write!(f, "{}", kind_symbol(*k)),
        NoncovalentBondKindAst::Undetermined => write!(f, "*"),
    }
}

/// Surface DSL wrapper around a [`NoncovalentBondUpdate`].
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct NoncovalentBondUpdateDsl(pub NoncovalentBondUpdate);

impl FromStr for NoncovalentBondUpdateDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_noncovalent_bond_update(s)
    }
}

impl Display for NoncovalentBondUpdateDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if let Some(kind) = &self.0.kind {
            fmt_kind(f, kind)?;
        }
        for c in self.0.constraints.iter() {
            if c.is_undetermined() {
                write!(f, "{}*", constraint_tag(c))?;
            } else {
                fmt_constraint(f, c)?;
            }
        }
        Ok(())
    }
}

impl<'de> FromEdn<'de> for NoncovalentBondUpdateDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s
                .parse()
                .map_err(|e| DeError::subgrammar("noncovalent-bond-update", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for NoncovalentBondUpdateDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

fn noncovalent_bond_update(i: &mut &str) -> PResult<NoncovalentBondUpdateDsl> {
    multispace0.parse_next(i)?;
    let kind = if i.starts_with('*') || i.as_bytes().first().is_some_and(u8::is_ascii_uppercase) {
        Some(terminated(kind_expr, multispace0).parse_next(i)?)
    } else {
        None
    };
    let preds: Vec<NoncovalentBondPredicate> =
        repeat(0.., terminated(noncovalent_bond_predicate, multispace0)).parse_next(i)?;
    let mut update = NoncovalentBondUpdate {
        kind,
        ..Default::default()
    };
    apply_update_predicates(&mut update, preds).map_err(ErrMode::Cut)?;
    Ok(NoncovalentBondUpdateDsl(update))
}

fn apply_update_predicates(
    update: &mut NoncovalentBondUpdate,
    preds: Vec<NoncovalentBondPredicate>,
) -> Result<(), ParseError> {
    for pred in preds {
        let NoncovalentBondPredicate::Constraint(constraint) = pred;
        if update.constraints.contains(constraint.key()) {
            return Err(ParseError::DuplicateNoncovalentBondPredicate(
                constraint_tag(&constraint).to_string(),
            ));
        }
        update.constraints.set(constraint);
    }
    Ok(())
}

/// Surface DSL wrapper around the narrow `NoncovalentBondConstraintAst`. EDN form is a single-key
/// map `{:intramolecular <bool>}`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum NoncovalentBondConstraintDsl {
    Intramolecular(BooleanAst),
}

impl<'de> FromEdn<'de> for NoncovalentBondConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Map(m) => {
                if m.len() != 1 {
                    return Err(DeError::Custom(format!(
                        "noncovalent-bond-constraint must have exactly one key, got {}",
                        m.len()
                    )));
                }
                let (k, v) = m.iter().next().unwrap();
                let Edn::Keyword(key) = k else {
                    return Err(DeError::TypeMismatch {
                        expected: "keyword key",
                        got: k.kind(),
                        path: vec!["noncovalent-bond-constraint".into()],
                    });
                };
                Ok(match key.name() {
                    "intramolecular" => Self::Intramolecular(BooleanDsl::from_edn(v)?.0),
                    other => {
                        return Err(DeError::UnknownField {
                            key: other.to_string(),
                            path: vec!["noncovalent-bond-constraint".into()],
                        });
                    }
                })
            }
            other => Err(DeError::TypeMismatch {
                expected: "{:intramolecular …}",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for NoncovalentBondConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        match self {
            Self::Intramolecular(b) => single_key_map("intramolecular", BooleanDsl(*b).to_edn()),
        }
    }
}

impl NoncovalentBondConstraintDsl {
    /// Build from the narrow inline AST form.
    pub(crate) fn from_ast(c: &NoncovalentBondConstraintAst) -> Self {
        match c {
            NoncovalentBondConstraintAst::Intramolecular(b) => Self::Intramolecular(*b),
        }
    }

    /// Convert into the narrow inline AST form.
    pub(crate) fn into_ast(self) -> NoncovalentBondConstraintAst {
        match self {
            Self::Intramolecular(b) => NoncovalentBondConstraintAst::Intramolecular(b),
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::read_string;

    use super::*;
    use crate::ast::constraint::NoncovalentBondConstraintsAst;

    #[rustfmt::skip]
    #[rstest]
    #[case::hbond("Hbd", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond)))]
    #[case::xbond("Xbd", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::HalogenBond)))]
    #[case::ybond("Ybd", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::ChalcogenBond)))]
    #[case::ion("Ion", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::Ionic)))]
    #[case::vdw("Vdw", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::VanDerWaals)))]
    #[case::whitespace("  Hbd  ", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond)))]
    #[case::undetermined("*", NoncovalentBondDsl(NoncovalentBondAst::new(NoncovalentBondKindAst::Undetermined)))]
    #[case::intramolecular("Hbd#I", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintAst::intramolecular(true))))]
    #[case::intramolecular_plus("Hbd#I+", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintAst::intramolecular(true))))]
    #[case::intermolecular("Hbd#I!", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintAst::intramolecular(false))))]
    #[case::intramolecular_undetermined("Hbd#I*", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintAst::Intramolecular(BooleanAst::Undetermined))))]
    #[case::undetermined_kind_with_pred("*#I", NoncovalentBondDsl(NoncovalentBondAst::new(NoncovalentBondKindAst::Undetermined).with_constraint(NoncovalentBondConstraintAst::intramolecular(true))))]
    fn test_parse_noncovalent_bond(#[case] input: &str, #[case] expected: NoncovalentBondDsl) {
        assert_eq!(parse_noncovalent_bond(input).unwrap(), expected);
    }

    #[rstest]
    #[case::empty("", ParseError::ExpectedNoncovalentBondKind)]
    #[case::unknown_literal("Abc", ParseError::ExpectedNoncovalentBondKind)]
    #[case::two_letter("Hb", ParseError::ExpectedNoncovalentBondKind)]
    #[case::bare_paren("(", ParseError::ExpectedNoncovalentBondKind)]
    #[case::unknown_predicate("Hbd#z", ParseError::UnknownNoncovalentBondPredicate("#z".into()))]
    #[case::duplicate("Hbd#I#I", ParseError::DuplicateNoncovalentBondPredicate("#I".into()))]
    fn test_parse_noncovalent_bond_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(parse_noncovalent_bond(input).unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::hydrogen_bond(NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond)), "Hbd")]
    #[case::intermolecular(NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintAst::intramolecular(false))), "Hbd#I!")]
    #[case::undetermined_constraint(NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintAst::Intramolecular(BooleanAst::Undetermined))), "Hbd")]
    fn test_noncovalent_bond_dsl_display(
        #[case] input: NoncovalentBondDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_string(), expected);
    }

    #[rstest]
    #[case::hbond("Hbd")]
    #[case::ion("Ion")]
    #[case::undetermined("*")]
    #[case::intramolecular("Hbd#I")]
    #[case::intermolecular("Hbd#I!")]
    fn test_noncovalent_bond_dsl_display_roundtrip(#[case] input: &str) {
        let dsl = parse_noncovalent_bond(input).unwrap();
        assert_eq!(parse_noncovalent_bond(&dsl.to_string()).unwrap(), dsl);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::hydrogen_bond(r##""Hbd""##, NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond)))]
    #[case::undetermined(r##""*""##, NoncovalentBondDsl(NoncovalentBondAst::default()))]
    fn test_noncovalent_bond_dsl_from_edn(
        #[case] input: &str,
        #[case] expected: NoncovalentBondDsl,
    ) {
        assert_eq!(
            NoncovalentBondDsl::from_edn(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rstest]
    #[case::wrong_type("1", DeError::TypeMismatch { expected: "string", got: "int", path: Vec::new() })]
    fn test_noncovalent_bond_dsl_from_edn_error(#[case] input: &str, #[case] expected: DeError) {
        assert_eq!(
            NoncovalentBondDsl::from_edn(&read_string(input).unwrap()).unwrap_err(),
            expected,
        );
    }

    #[rstest]
    #[case::hydrogen_bond(r##""Hbd""##)]
    #[case::undetermined(r##""*""##)]
    fn test_noncovalent_bond_dsl_from_edn_parity(#[case] input: &str) {
        assert_eq!(
            NoncovalentBondDsl::from_edn_str(input).unwrap(),
            NoncovalentBondDsl::from_edn(&read_string(input).unwrap()).unwrap(),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::hydrogen_bond(NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond)), r##""Hbd""##)]
    #[case::undetermined(NoncovalentBondDsl(NoncovalentBondAst::default()), r##""*""##)]
    fn test_noncovalent_bond_dsl_to_edn(
        #[case] input: NoncovalentBondDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::hydrogen_bond(
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
        NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond)),
    )]
    fn test_noncovalent_bond_dsl_from_ast(
        #[case] input: NoncovalentBondAst,
        #[case] expected: NoncovalentBondDsl,
    ) {
        assert_eq!(
            NoncovalentBondDsl::from_ast(&input, &NoncovalentBondDefaults::zeroed()),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::hydrogen_bond(
        NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond)),
        NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond),
    )]
    fn test_noncovalent_bond_dsl_into_ast(
        #[case] input: NoncovalentBondDsl,
        #[case] expected: NoncovalentBondAst,
    ) {
        assert_eq!(input.into_ast(&NoncovalentBondDefaults::zeroed()), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::hbond("Hbd", NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond))]
    #[case::xbond("Xbd", NoncovalentBondAst::from_kind(NoncovalentBondKind::HalogenBond))]
    #[case::ybond("Ybd", NoncovalentBondAst::from_kind(NoncovalentBondKind::ChalcogenBond))]
    fn test_noncovalent_bond_ast_from_str(
        #[case] input: &str,
        #[case] expected: NoncovalentBondAst,
    ) {
        assert_eq!(input.parse::<NoncovalentBondAst>().unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::intramolecular(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintAst::intramolecular(true)), "Hbd#I")]
    #[case::intermolecular(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintAst::intramolecular(false)), "Hbd#I!")]
    #[case::undetermined_constraint(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(NoncovalentBondConstraintAst::Intramolecular(BooleanAst::Undetermined)), "Hbd")]
    fn test_noncovalent_bond_ast_display(
        #[case] input: NoncovalentBondAst,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", NoncovalentBondUpdateDsl(NoncovalentBondUpdate::default()))]
    #[case::kind("Hbd", NoncovalentBondUpdateDsl(NoncovalentBondUpdate { kind: Some(NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond)), ..Default::default() }))]
    #[case::kind_undetermined("*", NoncovalentBondUpdateDsl(NoncovalentBondUpdate { kind: Some(NoncovalentBondKindAst::Undetermined), ..Default::default() }))]
    #[case::constraint("#I", NoncovalentBondUpdateDsl(NoncovalentBondUpdate { constraints: NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(true)), ..Default::default() }))]
    #[case::constraint_removal("#I*", NoncovalentBondUpdateDsl(NoncovalentBondUpdate { constraints: NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::Intramolecular(BooleanAst::Undetermined)), ..Default::default() }))]
    #[case::kind_and_constraint("Ion#I!", NoncovalentBondUpdateDsl(NoncovalentBondUpdate { kind: Some(NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic)), constraints: NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(false)) }))]
    fn test_parse_noncovalent_bond_update(
        #[case] input: &str,
        #[case] expected: NoncovalentBondUpdateDsl,
    ) {
        assert_eq!(parse_noncovalent_bond_update(input).unwrap(), expected);
    }

    #[rstest]
    #[case::unknown_predicate("#z", ParseError::UnknownNoncovalentBondPredicate("#z".into()))]
    #[case::duplicate("#I#I", ParseError::DuplicateNoncovalentBondPredicate("#I".into()))]
    fn test_parse_noncovalent_bond_update_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(parse_noncovalent_bond_update(input).unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(NoncovalentBondUpdateDsl(NoncovalentBondUpdate::default()), "")]
    #[case::kind_undetermined(NoncovalentBondUpdateDsl(NoncovalentBondUpdate { kind: Some(NoncovalentBondKindAst::Undetermined), ..Default::default() }), "*")]
    #[case::constraint_removal(NoncovalentBondUpdateDsl(NoncovalentBondUpdate { constraints: NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::Intramolecular(BooleanAst::Undetermined)), ..Default::default() }), "#I*")]
    #[case::kind_and_constraint(NoncovalentBondUpdateDsl(NoncovalentBondUpdate { kind: Some(NoncovalentBondKindAst::Lit(NoncovalentBondKind::Ionic)), constraints: NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::intramolecular(false)) }), "Ion#I!")]
    fn test_noncovalent_bond_update_dsl_display(
        #[case] input: NoncovalentBondUpdateDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(r##""""##, NoncovalentBondUpdateDsl(NoncovalentBondUpdate::default()))]
    #[case::kind_undetermined(r##""*""##, NoncovalentBondUpdateDsl(NoncovalentBondUpdate { kind: Some(NoncovalentBondKindAst::Undetermined), ..Default::default() }))]
    #[case::constraint_removal(r##""#I*""##, NoncovalentBondUpdateDsl(NoncovalentBondUpdate { constraints: NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::Intramolecular(BooleanAst::Undetermined)), ..Default::default() }))]
    fn test_noncovalent_bond_update_dsl_from_edn(
        #[case] input: &str,
        #[case] expected: NoncovalentBondUpdateDsl,
    ) {
        assert_eq!(
            NoncovalentBondUpdateDsl::from_edn(&read_string(input).unwrap()).unwrap(),
            expected,
        );
    }

    #[rstest]
    #[case::wrong_type("1", DeError::TypeMismatch { expected: "string", got: "int", path: Vec::new() })]
    fn test_noncovalent_bond_update_dsl_from_edn_error(
        #[case] input: &str,
        #[case] expected: DeError,
    ) {
        assert_eq!(
            NoncovalentBondUpdateDsl::from_edn(&read_string(input).unwrap()).unwrap_err(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty(NoncovalentBondUpdateDsl(NoncovalentBondUpdate::default()), r##""""##)]
    #[case::kind_undetermined(NoncovalentBondUpdateDsl(NoncovalentBondUpdate { kind: Some(NoncovalentBondKindAst::Undetermined), ..Default::default() }), r##""*""##)]
    #[case::constraint_removal(NoncovalentBondUpdateDsl(NoncovalentBondUpdate { constraints: NoncovalentBondConstraintsAst::from(NoncovalentBondConstraintAst::Intramolecular(BooleanAst::Undetermined)), ..Default::default() }), r##""#I*""##)]
    fn test_noncovalent_bond_update_dsl_to_edn(
        #[case] input: NoncovalentBondUpdateDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(input.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::true_("{:intramolecular true}", NoncovalentBondConstraintDsl::Intramolecular(BooleanAst::Lit(true)))]
    #[case::false_("{:intramolecular false}", NoncovalentBondConstraintDsl::Intramolecular(BooleanAst::Lit(false)))]
    #[case::undetermined("{:intramolecular :undetermined}", NoncovalentBondConstraintDsl::Intramolecular(BooleanAst::Undetermined))]
    fn test_noncovalent_bond_constraint_dsl_from_edn(
        #[case] input: &str,
        #[case] expected: NoncovalentBondConstraintDsl,
    ) {
        let edn = read_string(input).unwrap();
        assert_eq!(NoncovalentBondConstraintDsl::from_edn(&edn).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::unknown_key("{:contains 1}", DeError::UnknownField { key: "contains".to_string(), path: vec!["noncovalent-bond-constraint".into()] })]
    #[case::two_keys("{:intramolecular true, :contains 1}", DeError::Custom("noncovalent-bond-constraint must have exactly one key, got 2".to_string()))]
    #[case::not_a_map("42", DeError::TypeMismatch { expected: "{:intramolecular …}", got: "int", path: Vec::new() })]
    fn test_noncovalent_bond_constraint_dsl_from_edn_error(
        #[case] input: &str,
        #[case] expected: DeError,
    ) {
        assert_eq!(
            NoncovalentBondConstraintDsl::from_edn(&read_string(input).unwrap()).unwrap_err(),
            expected,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::true_(NoncovalentBondConstraintDsl::Intramolecular(BooleanAst::Lit(true)), "{:intramolecular true}")]
    #[case::false_(NoncovalentBondConstraintDsl::Intramolecular(BooleanAst::Lit(false)), "{:intramolecular false}")]
    #[case::undetermined(NoncovalentBondConstraintDsl::Intramolecular(BooleanAst::Undetermined), "{:intramolecular :undetermined}")]
    fn test_noncovalent_bond_constraint_dsl_to_edn(
        #[case] dsl: NoncovalentBondConstraintDsl,
        #[case] expected: &str,
    ) {
        assert_eq!(dsl.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::intramolecular(NoncovalentBondConstraintAst::intramolecular(true), NoncovalentBondConstraintDsl::Intramolecular(BooleanAst::Lit(true)))]
    #[case::undetermined(NoncovalentBondConstraintAst::Intramolecular(BooleanAst::Undetermined), NoncovalentBondConstraintDsl::Intramolecular(BooleanAst::Undetermined))]
    fn test_noncovalent_bond_constraint_dsl_from_ast(
        #[case] ast: NoncovalentBondConstraintAst,
        #[case] expected: NoncovalentBondConstraintDsl,
    ) {
        assert_eq!(NoncovalentBondConstraintDsl::from_ast(&ast), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::intramolecular(NoncovalentBondConstraintDsl::Intramolecular(BooleanAst::Lit(false)), NoncovalentBondConstraintAst::intramolecular(false))]
    fn test_noncovalent_bond_constraint_dsl_into_ast(
        #[case] dsl: NoncovalentBondConstraintDsl,
        #[case] expected: NoncovalentBondConstraintAst,
    ) {
        assert_eq!(dsl.into_ast(), expected);
    }
}
