//! Noncovalent-bond-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, FromEdn, ToEdn};
use winnow::ascii::multispace0;
use winnow::combinator::{alt, delimited, preceded, separated, terminated};
use winnow::error::{ErrMode, ParserError};
use winnow::token::one_of;
use winnow::Parser;

use super::error::{PResult, ParseError};
use super::value::id;
use crate::ast::config::NoncovalentBondAstConfig;
use crate::ast::noncovalent::{NoncovalentBondAst, NoncovalentKind, NoncovalentKindAst};
use crate::ast::traits::{FromAst, ToAst};

/// Surface DSL wrapper around `NoncovalentBondAst`. String form is the
/// noncovalent-kind expression (three-letter literal, set, bind, ref, or `*`).
/// All `NoncovalentBondConstraint` variants are molecule-scope.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct NoncovalentBondDsl(pub NoncovalentBondAst);

impl FromStr for NoncovalentBondDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_noncovalent(s)
    }
}

impl Display for NoncovalentBondDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_noncovalent_ast(f, &self.0)
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
}

impl ToEdn for NoncovalentBondDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

impl FromAst<NoncovalentBondAst> for NoncovalentBondDsl {
    type Error = ParseError;

    fn from_ast(
        ast: &NoncovalentBondAst,
        _cfg: &NoncovalentBondAstConfig,
    ) -> Result<Self, ParseError> {
        Ok(NoncovalentBondDsl(ast.clone()))
    }
}

impl ToAst<NoncovalentBondAst> for NoncovalentBondDsl {
    type Error = ParseError;

    fn to_ast(&self, _cfg: &NoncovalentBondAstConfig) -> Result<NoncovalentBondAst, ParseError> {
        Ok(self.0.clone())
    }
}

// -- Parse --------------------

pub fn parse_noncovalent(input: &str) -> Result<NoncovalentBondDsl, ParseError> {
    noncovalent.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn noncovalent(i: &mut &str) -> PResult<NoncovalentBondDsl> {
    let kind = delimited(multispace0, kind_expr, multispace0).parse_next(i)?;
    Ok(NoncovalentBondDsl(NoncovalentBondAst::new(kind)))
}

fn kind_expr(i: &mut &str) -> PResult<NoncovalentKindAst> {
    alt((
        '*'.value(NoncovalentKindAst::Undetermined),
        kind_set.map(NoncovalentKindAst::Set),
        kind_bind.map(|(id, set)| NoncovalentKindAst::Bind { id, set }),
        kind_ref.map(NoncovalentKindAst::Ref),
        kind_literal.map(NoncovalentKindAst::Lit),
    ))
    .parse_next(i)
    .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedNoncovalentKind))
}

fn kind_literal(i: &mut &str) -> PResult<NoncovalentKind> {
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

fn kind_set(i: &mut &str) -> PResult<Vec<NoncovalentKind>> {
    delimited(
        '{',
        delimited(
            multispace0,
            separated(1.., kind_literal, delimited(multispace0, ',', multispace0)),
            multispace0,
        ),
        '}',
    )
    .parse_next(i)
}

fn kind_bind(i: &mut &str) -> PResult<(String, Vec<NoncovalentKind>)> {
    delimited(
        '(',
        (
            delimited(multispace0, preceded('?', id), multispace0),
            preceded(("::", multispace0), terminated(kind_set, multispace0)),
        ),
        ')',
    )
    .parse_next(i)
}

fn kind_ref(i: &mut &str) -> PResult<String> {
    delimited(
        '(',
        delimited(multispace0, preceded('?', id), multispace0),
        ')',
    )
    .parse_next(i)
}

fn kind_from_symbol(sym: &str) -> Option<NoncovalentKind> {
    match sym {
        "Hbd" => Some(NoncovalentKind::HydrogenBond),
        "Xbd" => Some(NoncovalentKind::HalogenBond),
        "Ybd" => Some(NoncovalentKind::ChalcogenBond),
        "Ion" => Some(NoncovalentKind::Ionic),
        "Vdw" => Some(NoncovalentKind::VanDerWaals),
        _ => None,
    }
}

fn kind_symbol(k: NoncovalentKind) -> &'static str {
    match k {
        NoncovalentKind::HydrogenBond => "Hbd",
        NoncovalentKind::HalogenBond => "Xbd",
        NoncovalentKind::ChalcogenBond => "Ybd",
        NoncovalentKind::Ionic => "Ion",
        NoncovalentKind::VanDerWaals => "Vdw",
    }
}

// -- Format --------------------

fn fmt_noncovalent_ast(f: &mut fmt::Formatter<'_>, ast: &NoncovalentBondAst) -> fmt::Result {
    fmt_kind(f, &ast.kind)
}

fn fmt_kind(f: &mut fmt::Formatter<'_>, kind: &NoncovalentKindAst) -> fmt::Result {
    match kind {
        NoncovalentKindAst::Lit(k) => write!(f, "{}", kind_symbol(*k)),
        NoncovalentKindAst::Undetermined => write!(f, "*"),
        NoncovalentKindAst::Set(ks) => {
            write!(f, "{{")?;
            for (i, k) in ks.iter().enumerate() {
                if i > 0 {
                    write!(f, ",")?;
                }
                write!(f, "{}", kind_symbol(*k))?;
            }
            write!(f, "}}")
        }
        NoncovalentKindAst::Bind { id, set } => {
            write!(f, "(?{} :: {{", id)?;
            for (i, k) in set.iter().enumerate() {
                if i > 0 {
                    write!(f, ",")?;
                }
                write!(f, "{}", kind_symbol(*k))?;
            }
            write!(f, "}})")
        }
        NoncovalentKindAst::Ref(id) => write!(f, "(?{})", id),
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::hbond("Hbd", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentKind::HydrogenBond)))]
    #[case::xbond("Xbd", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentKind::HalogenBond)))]
    #[case::ybond("Ybd", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentKind::ChalcogenBond)))]
    #[case::ion("Ion", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentKind::Ionic)))]
    #[case::vdw("Vdw", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentKind::VanDerWaals)))]
    #[case::whitespace("  Hbd  ", NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentKind::HydrogenBond)))]
    #[case::undetermined("*", NoncovalentBondDsl(NoncovalentBondAst::new(NoncovalentKindAst::Undetermined)))]
    #[case::set("{Hbd,Ion}", NoncovalentBondDsl(NoncovalentBondAst::new(NoncovalentKindAst::Set(vec![NoncovalentKind::HydrogenBond, NoncovalentKind::Ionic]))))]
    #[case::set_spaced("{ Hbd, Vdw }", NoncovalentBondDsl(NoncovalentBondAst::new(NoncovalentKindAst::Set(vec![NoncovalentKind::HydrogenBond, NoncovalentKind::VanDerWaals]))))]
    #[case::bind("(?k :: {Hbd,Ion})", NoncovalentBondDsl(NoncovalentBondAst::new(NoncovalentKindAst::Bind { id: "k".to_string(), set: vec![NoncovalentKind::HydrogenBond, NoncovalentKind::Ionic] })))]
    #[case::ref_("(?k)", NoncovalentBondDsl(NoncovalentBondAst::new(NoncovalentKindAst::Ref("k".to_string()))))]
    fn test_parse_noncovalent(#[case] input: &str, #[case] expected: NoncovalentBondDsl) {
        let result = noncovalent.parse(input);
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
        let result = noncovalent.parse(input);
        assert!(result.is_err(), "{:?} should fail", input);
    }

    #[rstest]
    #[case::hbond("Hbd")]
    #[case::ion("Ion")]
    #[case::undetermined("*")]
    #[case::set("{Hbd,Ion}")]
    #[case::bind("(?k :: {Hbd,Ion})")]
    #[case::ref_("(?k)")]
    fn test_noncovalent_roundtrip(#[case] input: &str) {
        let form: NoncovalentBondDsl = input.parse().unwrap();
        let rendered = form.to_string();
        let reparsed: NoncovalentBondDsl = rendered.parse().unwrap();
        assert_eq!(form, reparsed);
    }

    #[rstest]
    fn test_noncovalent_dsl_to_ast_passthrough() {
        let dsl = NoncovalentBondDsl(NoncovalentBondAst::from_kind(NoncovalentKind::HydrogenBond));
        let cfg = NoncovalentBondAstConfig::zeroed();
        let ast = dsl.to_ast(&cfg).unwrap();
        assert_eq!(
            ast.kind,
            NoncovalentKindAst::Lit(NoncovalentKind::HydrogenBond)
        );
    }
}
