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
use crate::ast::noncovalent::{NoncovalentBondAst, NoncovalentBondKind, NoncovalentBondKindAst};
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

/// Partial noncovalent bond for a reaction `:modify` payload. The kind (`*` = unchanged) followed
/// by any `#I` intramolecular constraint overlay.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PartialNoncovalentBondDsl(pub NoncovalentBondAst);

impl FromStr for PartialNoncovalentBondDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Self(parse_noncovalent_bond(s)?.0))
    }
}

impl Display for PartialNoncovalentBondDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_kind(f, &self.0.kind)?;
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

impl<'de> FromEdn<'de> for PartialNoncovalentBondDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s
                .parse()
                .map_err(|e| DeError::subgrammar("noncovalent-bond", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for PartialNoncovalentBondDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
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
    #[case::unknown_predicate("Hbd#z", ParseError::UnknownNoncovalentBondPredicate("#z".into()))]
    #[case::duplicate("Hbd#I#I", ParseError::DuplicateNoncovalentBondPredicate("#I".into()))]
    fn test_parse_noncovalent_predicate_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(parse_noncovalent_bond(input).unwrap_err(), expected);
    }

    #[rstest]
    #[case::hbond("Hbd")]
    #[case::ion("Ion")]
    #[case::undetermined("*")]
    #[case::intramolecular("Hbd#I")]
    #[case::intermolecular("Hbd#I!")]
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

    #[rustfmt::skip]
    #[rstest]
    #[case::intramolecular(NoncovalentBondConstraintAst::intramolecular(true), "Hbd#I")]
    #[case::intermolecular(NoncovalentBondConstraintAst::intramolecular(false), "Hbd#I!")]
    #[case::undetermined_elides(NoncovalentBondConstraintAst::Intramolecular(BooleanAst::Undetermined), "Hbd")]
    fn test_fmt_noncovalent_bond_ast(#[case] c: NoncovalentBondConstraintAst, #[case] expected: &str) {
        let bond = NoncovalentBondAst::from_kind(NoncovalentBondKind::HydrogenBond).with_constraint(c);
        assert_eq!(bond.to_string(), expected);
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
    #[case::unknown_key("{:contains 1}", |e: &DeError| matches!(e, DeError::UnknownField { .. }))]
    #[case::two_keys("{:intramolecular true, :contains 1}", |e: &DeError| matches!(e, DeError::Custom(_)))]
    #[case::not_a_map("42", |e: &DeError| matches!(e, DeError::TypeMismatch { .. }))]
    fn test_noncovalent_bond_constraint_dsl_from_edn_error(
        #[case] input: &str,
        #[case] is_expected: impl Fn(&DeError) -> bool,
    ) {
        let edn = read_string(input).unwrap();
        let err = NoncovalentBondConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(is_expected(&err), "unexpected error: {err:?}");
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

    #[rstest]
    #[case::hbond("Hbd")]
    #[case::xbond("Xbd")]
    #[case::ybond("Ybd")]
    fn test_noncovalent_bond_ast_from_str_to_string_roundtrip(#[case] s: &str) {
        let ast: NoncovalentBondAst = s.parse().unwrap();
        assert_eq!(ast.to_string(), s);
    }
}
