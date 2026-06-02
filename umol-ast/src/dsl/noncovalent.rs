//! Noncovalent-bond-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, EdnError, EdnStreamDeserializer, FromEdn, ToEdn};
use winnow::ascii::multispace0;
use winnow::combinator::{alt, delimited};
use winnow::error::{ErrMode, ParserError};
use winnow::token::one_of;
use winnow::Parser;

use super::config::NoncovalentBondDefaults;
use super::error::{PResult, ParseError};
use crate::ast::noncovalent::{NoncovalentBondAst, NoncovalentBondKind, NoncovalentBondKindAst};
use crate::ast::traits::{FromAst, IntoAst};

/// Surface DSL wrapper around `NoncovalentBondAst`. String form is the
/// noncovalent-kind expression (three-letter literal, set, bind, ref, or `*`).
/// All `NoncovalentBondConstraint` variants are molecule-scope.
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

    fn from_ast(ast: &NoncovalentBondAst, _cfg: &Self::Ctx) -> Self {
        NoncovalentBondDsl(ast.clone())
    }
}

impl IntoAst<NoncovalentBondAst> for NoncovalentBondDsl {
    type Ctx = NoncovalentBondDefaults;

    fn into_ast(self, _cfg: &Self::Ctx) -> NoncovalentBondAst {
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

pub(crate) fn noncovalent_bond(i: &mut &str) -> PResult<NoncovalentBondDsl> {
    let kind = delimited(multispace0, kind_expr, multispace0).parse_next(i)?;
    Ok(NoncovalentBondDsl(NoncovalentBondAst::new(kind)))
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
    fmt_kind(f, &ast.kind)
}

fn fmt_kind(f: &mut fmt::Formatter<'_>, kind: &NoncovalentBondKindAst) -> fmt::Result {
    match kind {
        NoncovalentBondKindAst::Lit(k) => write!(f, "{}", kind_symbol(*k)),
        NoncovalentBondKindAst::Undetermined => write!(f, "*"),
    }
}

/// Surface DSL wrapper around the narrow `NoncovalentBondConstraint`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum NoncovalentBondConstraintDsl {}

impl<'de> FromEdn<'de> for NoncovalentBondConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        Err(DeError::TypeMismatch {
            expected: "no value-only noncovalent-bond constraints exist yet",
            got: edn.kind(),
            path: vec!["noncovalent-bond-constraint".into()],
        })
    }
}

impl ToEdn for NoncovalentBondConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        match *self {}
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::read_string;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::hbond("Hbd", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond)))]
    #[case::xbond("Xbd", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::HalogenBond)))]
    #[case::ybond("Ybd", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::ChalcogenBond)))]
    #[case::ion("Ion", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::Ionic)))]
    #[case::vdw("Vdw", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::VanDerWaals)))]
    #[case::whitespace("  Hbd  ", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond)))]
    #[case::undetermined("*", NoncovalentBondDsl(NoncovalentBondAst::new(NoncovalentBondKindAst::Undetermined)))]
    fn test_parse_noncovalent(#[case] input: &str, #[case] expected: NoncovalentBondDsl) {
        let result = noncovalent_bond.parse(input);
        assert!(result.is_ok(), "{:?} should succeed, got {:?}", input, result.clone().unwrap_err());
        let form = result.unwrap();
        assert_eq!(form, expected);
    }

    #[rstest]
    #[case::empty("")]
    #[case::unknown_literal("Abc")]
    #[case::two_letter("Hb")]
    #[case::bare_paren("(")]
    fn test_parse_noncovalent_invalid(#[case] input: &str) {
        let result = noncovalent_bond.parse(input);
        assert!(result.is_err(), "{:?} should fail", input);
    }

    #[rstest]
    #[case::hbond("Hbd")]
    #[case::ion("Ion")]
    #[case::undetermined("*")]
    fn test_noncovalent_roundtrip(#[case] input: &str) {
        let form: NoncovalentBondDsl = input.parse().unwrap();
        let rendered = form.to_string();
        let reparsed: NoncovalentBondDsl = rendered.parse().unwrap();
        assert_eq!(form, reparsed);
    }

    #[rstest]
    fn test_noncovalent_dsl_to_ast_passthrough() {
        let dsl = NoncovalentBondDsl(NoncovalentBondAst::from_kind(
            NoncovalentBondKind::HydrogenBond,
        ));
        let cfg = NoncovalentBondDefaults::zeroed();
        let ast = dsl.into_ast(&cfg);
        assert_eq!(
            ast.kind,
            NoncovalentBondKindAst::Lit(NoncovalentBondKind::HydrogenBond)
        );
    }

    #[rstest]
    #[case::single(r##""Hbd""##)]
    #[case::undetermined(r##""*""##)]
    fn test_noncovalent_dsl_from_edn_str_matches_from_edn(#[case] input: &str) {
        let via_stream = NoncovalentBondDsl::from_edn_str(input).unwrap();
        let tree = read_string(input).unwrap();
        let via_tree = NoncovalentBondDsl::from_edn(&tree).unwrap();
        assert_eq!(via_stream, via_tree);
    }

    #[rstest]
    fn test_noncovalent_bond_constraint_dsl_from_edn_errors() {
        let edn = read_string("{:contains 1}").unwrap();
        let err = NoncovalentBondConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::TypeMismatch { .. }));
    }

    #[rstest]
    #[case::hbond("Hbd")]
    #[case::xbond("Xbd")]
    #[case::ybond("Ybd")]
    fn test_noncovalent_bond_ast_from_str_to_string_roundtrip(#[case] s: &str) {
        let ast: NoncovalentBondAst = s.parse().unwrap();
        assert_eq!(ast.to_string(), s);
    }
}
