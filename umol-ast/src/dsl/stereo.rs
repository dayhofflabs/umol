//! Stereo config-string DSL

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnStreamDeserializer, FromEdn, ToEdn};
use umol_perm::Permutation;
use winnow::ascii::{digit1, multispace0};
use winnow::combinator::{alt, delimited, opt, preceded, separated, terminated};
use winnow::error::ErrMode;
use winnow::Parser;

use super::error::{PResult, ParseError};
use super::value::{id, terminator};
use crate::ast::stereo::{
    StereoAtomAst, StereoBondAst, StereoConfigurationAst, StereoExpr, StereoIndexAst, StereoKind,
};

/// Surface DSL wrapper for `StereoAtomAst`
#[repr(transparent)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoAtomDsl(pub StereoAtomAst);

impl StereoAtomDsl {
    /// Zero-cost reference cast from `&StereoAtomAst`. Relies on `repr(transparent)`.
    pub fn from_ref(ast: &StereoAtomAst) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(ast as *const StereoAtomAst as *const Self) }
    }
}

impl From<StereoAtomAst> for StereoAtomDsl {
    fn from(ast: StereoAtomAst) -> Self {
        Self(ast)
    }
}

impl FromStr for StereoAtomDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_stereo_atom(s)
    }
}

impl Display for StereoAtomDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_stereo_atom(f, &self.0)?;
        Ok(())
    }
}

impl<'de> FromEdn<'de> for StereoAtomDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s.parse().map_err(|e| DeError::subgrammar("stereo atom", e)),
            Edn::Keyword(k) => {
                let s = expand_stereo_atom_keyword(k.name()).ok_or_else(|| {
                    DeError::Custom(format!("unknown stereo atom keyword :{}", k.name()))
                })?;
                s.parse().map_err(|e| DeError::subgrammar("stereo atom", e))
            }
            other => Err(DeError::TypeMismatch {
                expected: "string or stereo atom keyword",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        EdnStreamDeserializer::new(input).read_subgrammar_all("stereo atom")
    }
}

/// Expand a stereo atom keyword shorthand to its equivalent stereo atom string payload.
/// The recognized keywords are as follows:
///
/// - `:ccw` -> `"Th1"`
/// - `:cw` -> `"Th2"`
///
/// Returns `None` for unrecorgnized keywords.
pub(crate) fn expand_stereo_atom_keyword(name: &str) -> Option<&'static str> {
    match name {
        "ccw" => Some("Th1"),
        "cw" => Some("Th2"),
        _ => None,
    }
}

impl ToEdn for StereoAtomDsl {
    fn to_edn(&self) -> Edn<'static> {
        match stereo_atom_keyword_for(&self.0) {
            Some(kw) => Edn::Keyword(EdnKeyword::owned(kw.to_string())),
            None => Edn::Str(Cow::Owned(self.to_string())),
        }
    }
}

/// Return the stereo atom keyword for canonical stereo atom shapes, or `None`
/// when the full definition is required. Inverse of `expand_stereo_atom_keyword`.
fn stereo_atom_keyword_for(ast: &StereoAtomAst) -> Option<&'static str> {
    match (ast.kind, &ast.configuration) {
        (StereoKind::Tetrahedral, &StereoConfigurationAst::Stereo(StereoIndexAst::Lit(1))) => {
            Some("ccw")
        }
        (StereoKind::Tetrahedral, &StereoConfigurationAst::Stereo(StereoIndexAst::Lit(2))) => {
            Some("cw")
        }
        _ => None,
    }
}

pub fn parse_stereo_atom(input: &str) -> Result<StereoAtomDsl, ParseError> {
    stereo_atom.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn stereo_atom(i: &mut &str) -> PResult<StereoAtomDsl> {
    let kind = delimited(multispace0, stereo_kind, multispace0).parse_next(i)?;
    let config = terminated(stereo_config, multispace0).parse_next(i)?;
    Ok(StereoAtomDsl(StereoAtomAst::new(kind, config)))
}

/// Surface DSL wrapper for `StereoBondAst`
#[repr(transparent)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoBondDsl(pub StereoBondAst);

impl StereoBondDsl {
    /// Zero-cost reference cast from `&StereoBondAst`. Relies on `repr(transparent)`.
    pub fn from_ref(ast: &StereoBondAst) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(ast as *const StereoBondAst as *const Self) }
    }
}

impl From<StereoBondAst> for StereoBondDsl {
    fn from(ast: StereoBondAst) -> Self {
        Self(ast)
    }
}

impl FromStr for StereoBondDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_stereo_bond(s)
    }
}

impl Display for StereoBondDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_stereo_bond(f, &self.0)?;
        Ok(())
    }
}

impl<'de> FromEdn<'de> for StereoBondDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s.parse().map_err(|e| DeError::subgrammar("stereo bond", e)),
            Edn::Keyword(k) => {
                let s = expand_stereo_bond_keyword(k.name()).ok_or_else(|| {
                    DeError::Custom(format!("unknown stereo bond keyword :{}", k.name()))
                })?;
                s.parse().map_err(|e| DeError::subgrammar("stereo bond", e))
            }
            other => Err(DeError::TypeMismatch {
                expected: "string or stereo bond keyword",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        EdnStreamDeserializer::new(input).read_subgrammar_all("stereo bond")
    }
}

/// Expand a stereo bond keyword shorthand to its equivalent bond string payload.
/// The recognized keywords are as follows:
///
/// - `:z` -> `"Ct1"`
/// - `:e` -> `"Ct2"`
///
/// Returns `None` for unrecorgnized keywords.
pub(crate) fn expand_stereo_bond_keyword(name: &str) -> Option<&'static str> {
    match name {
        "z" => Some("Ct1"),
        "e" => Some("Ct2"),
        _ => None,
    }
}

impl ToEdn for StereoBondDsl {
    fn to_edn(&self) -> Edn<'static> {
        match stereo_bond_keyword_for(&self.0) {
            Some(kw) => Edn::Keyword(EdnKeyword::owned(kw.to_string())),
            None => Edn::Str(Cow::Owned(self.to_string())),
        }
    }
}

/// Return the stereo bond keyword for canonical stereo bond shapes, or `None`
/// when the full definition is required. Inverse of `expand_stereo_bond_keyword`.
fn stereo_bond_keyword_for(ast: &StereoBondAst) -> Option<&'static str> {
    match (ast.kind, &ast.configuration) {
        (StereoKind::CisTrans, &StereoConfigurationAst::Stereo(StereoIndexAst::Lit(1))) => {
            Some("z")
        }
        (StereoKind::CisTrans, &StereoConfigurationAst::Stereo(StereoIndexAst::Lit(2))) => {
            Some("e")
        }
        _ => None,
    }
}

pub fn parse_stereo_bond(input: &str) -> Result<StereoBondDsl, ParseError> {
    stereo_bond.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn stereo_bond(i: &mut &str) -> PResult<StereoBondDsl> {
    let kind = delimited(multispace0, stereo_kind, multispace0).parse_next(i)?;
    let config = terminated(stereo_config, multispace0).parse_next(i)?;
    Ok(StereoBondDsl(StereoBondAst::new(kind, config)))
}

/// `Th` / `Ct` / `Sp` / `Tb` / `Oh` symbol → `StereoKind`.
pub(crate) fn stereo_kind(i: &mut &str) -> PResult<StereoKind> {
    alt((
        "Th".value(StereoKind::Tetrahedral),
        "Ct".value(StereoKind::CisTrans),
        "Sp".value(StereoKind::SquarePlanar),
        "Tb".value(StereoKind::TrigonalBipyramidal),
        "Oh".value(StereoKind::Octahedral),
    ))
    .parse_next(i)
}

/// Parse the `config` grammar into a `StereoConfigurationAst`.
pub(crate) fn stereo_config(i: &mut &str) -> PResult<StereoConfigurationAst> {
    alt((
        '*'.value(StereoConfigurationAst::Undetermined),
        '!'.value(StereoConfigurationAst::NotStereo),
        '+'.value(StereoConfigurationAst::Stereo(StereoIndexAst::Undetermined)),
        terminated(stereo_lit, (multispace0, terminator))
            .map(|n| StereoConfigurationAst::Stereo(StereoIndexAst::Lit(n))),
        stereo_expr.map(|e| StereoConfigurationAst::Stereo(StereoIndexAst::Expr(Box::new(e)))),
    ))
    .parse_next(i)
}

/// `stereo-expr`: a `~`-prefixed base carrying zero or more left-associative
/// `^image` postfixes.
fn stereo_expr(i: &mut &str) -> PResult<StereoExpr> {
    let mut e = stereo_swap(i)?;
    loop {
        multispace0.parse_next(i)?;
        if opt('^').parse_next(i)?.is_some() {
            e = StereoExpr::apply(e, perm_image(i)?);
        } else {
            return Ok(e);
        }
    }
}

/// `'~' prefix-term | base` — unary `~` binds tighter than `^`.
fn stereo_swap(i: &mut &str) -> PResult<StereoExpr> {
    multispace0.parse_next(i)?;
    if opt('~').parse_next(i)?.is_some() {
        Ok(StereoExpr::swap(stereo_swap(i)?))
    } else {
        stereo_base(i)
    }
}

/// `nat | '?' id ('::' set)? | set`.
fn stereo_base(i: &mut &str) -> PResult<StereoExpr> {
    preceded(
        multispace0,
        alt((
            stereo_lit_set.map(StereoExpr::LitSet),
            stereo_var,
            stereo_lit.map(StereoExpr::Lit),
        )),
    )
    .parse_next(i)
}

fn stereo_var(i: &mut &str) -> PResult<StereoExpr> {
    let name = preceded('?', id).parse_next(i)?;
    let domain = opt(preceded((multispace0, "::", multispace0), stereo_lit_set)).parse_next(i)?;
    Ok(match domain {
        Some(set) => StereoExpr::VarDomain(name, set),
        None => StereoExpr::Var(name),
    })
}

fn stereo_lit(i: &mut &str) -> PResult<u32> {
    let s: &str = digit1.parse_next(i)?;
    s.parse::<u32>()
        .map_err(|_| ErrMode::Backtrack(ParseError::Syntax))
}

/// `'{' nat (',' nat)* '}'`, whitespace-insensitive.
fn stereo_lit_set(i: &mut &str) -> PResult<Vec<u32>> {
    preceded(
        ('{', multispace0),
        terminated(
            separated(1.., stereo_lit, (multispace0, ',', multispace0)),
            (multispace0, '}'),
        ),
    )
    .parse_next(i)
}

/// The 1-indexed one-line image read as digits → `Permutation`. Kept literal
/// (never canonicalized); validated as a bijection of `1..=degree`.
fn perm_image(i: &mut &str) -> PResult<Permutation> {
    let digits: &str = digit1.parse_next(i)?;
    let img: Vec<u8> = digits
        .bytes()
        .map(|b| b.checked_sub(b'1'))
        .collect::<Option<Vec<u8>>>()
        .ok_or(ErrMode::Cut(ParseError::Syntax))?;
    let degree = img.len();
    let mut seen = [false; 6];
    for &x in &img {
        let x = x as usize;
        if x >= degree || x >= 6 || seen[x] {
            return Err(ErrMode::Cut(ParseError::Syntax));
        }
        seen[x] = true;
    }
    Ok(Permutation::from_image(degree, &img))
}

/// Write the stereo atom DSL
pub(crate) fn fmt_stereo_atom(f: &mut fmt::Formatter<'_>, atom: &StereoAtomAst) -> fmt::Result {
    fmt_stereo_kind(f, atom.kind)?;
    fmt_stereo_config(f, &atom.configuration)?;
    Ok(())
}

/// Write the stereo bond DSL
pub(crate) fn fmt_stereo_bond(f: &mut fmt::Formatter<'_>, bond: &StereoBondAst) -> fmt::Result {
    fmt_stereo_kind(f, bond.kind)?;
    fmt_stereo_config(f, &bond.configuration)?;
    Ok(())
}

/// Write the stereo kind for `kind`.
pub(crate) fn fmt_stereo_kind(f: &mut fmt::Formatter<'_>, kind: StereoKind) -> fmt::Result {
    match kind {
        StereoKind::Tetrahedral => f.write_str("Th"),
        StereoKind::CisTrans => f.write_str("Ct"),
        StereoKind::SquarePlanar => f.write_str("Sp"),
        StereoKind::TrigonalBipyramidal => f.write_str("Tb"),
        StereoKind::Octahedral => f.write_str("Oh"),
    }
}

/// Write the `config` for `c` (no class head — the caller writes `#T`/`#C` or
/// the `Th`/`Ct` head).
pub(crate) fn fmt_stereo_config(
    f: &mut fmt::Formatter<'_>,
    config: &StereoConfigurationAst,
) -> fmt::Result {
    match config {
        StereoConfigurationAst::Undetermined => write!(f, "*"),
        StereoConfigurationAst::NotStereo => write!(f, "!"),
        StereoConfigurationAst::Stereo(StereoIndexAst::Undetermined) => write!(f, "+"),
        StereoConfigurationAst::Stereo(StereoIndexAst::Lit(n)) => write!(f, "{n}"),
        StereoConfigurationAst::Stereo(StereoIndexAst::Expr(e)) => fmt_stereo_expr(f, e),
    }
}

fn fmt_stereo_expr(f: &mut fmt::Formatter<'_>, e: &StereoExpr) -> fmt::Result {
    match e {
        StereoExpr::Lit(n) => write!(f, "{n}"),
        StereoExpr::Var(name) => write!(f, "?{name}"),
        StereoExpr::SwapOp(inner) => {
            write!(f, "~")?;
            fmt_stereo_expr(f, inner)
        }
        StereoExpr::ApplyOp(inner, perm) => {
            fmt_stereo_expr(f, inner)?;
            write!(f, "^")?;
            fmt_perm_image(f, *perm)
        }
        StereoExpr::LitSet(set) => fmt_stereo_lit_set(f, set),
        StereoExpr::VarDomain(name, set) => {
            write!(f, "?{name} :: ")?;
            fmt_stereo_lit_set(f, set)
        }
    }
}

fn fmt_stereo_lit_set(f: &mut fmt::Formatter<'_>, set: &[u32]) -> fmt::Result {
    write!(f, "{{")?;
    for (i, n) in set.iter().enumerate() {
        if i > 0 {
            write!(f, ",")?;
        }
        write!(f, "{n}")?;
    }
    write!(f, "}}")
}

fn fmt_perm_image(f: &mut fmt::Formatter<'_>, perm: Permutation) -> fmt::Result {
    for i in 0..perm.degree() {
        write!(f, "{}", perm.apply(i) + 1)?;
    }
    Ok(())
}

/// DSL for `StereoAtomConstraint`
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum StereoAtomConstraintDsl {}

impl<'de> FromEdn<'de> for StereoAtomConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        Err(DeError::TypeMismatch {
            expected: "no value-only stereo atom constraints exist yet",
            got: edn.kind(),
            path: vec!["stereo-atom-constraint".into()],
        })
    }
}

impl ToEdn for StereoAtomConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        match *self {}
    }
}

/// DSL for `StereoBondConstraint`
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum StereoBondConstraintDsl {}

impl<'de> FromEdn<'de> for StereoBondConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        Err(DeError::TypeMismatch {
            expected: "no value-only stereo bond constraints exist yet",
            got: edn.kind(),
            path: vec!["stereo-bond-constraint".into()],
        })
    }
}

impl ToEdn for StereoBondConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        match *self {}
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined("*", StereoConfigurationAst::Undetermined)]
    #[case::not_stereo("!", StereoConfigurationAst::NotStereo)]
    #[case::stereogenic("+", StereoConfigurationAst::Stereo(StereoIndexAst::Undetermined))]
    #[case::lit("1", StereoConfigurationAst::from(1_u32))]
    #[case::var("?o", StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::Var("o".to_string()))))]
    #[case::lit_set("{1,2}", StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::LitSet(vec![1, 2]))))]
    #[case::var_domain("?o :: {1,2}", StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::VarDomain("o".to_string(), vec![1, 2]))))]
    #[case::swap("~1", StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::swap(StereoExpr::Lit(1)))))]
    #[case::apply("1^2134", StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::apply(StereoExpr::Lit(1), Permutation::from_image(4, &[1, 0, 2, 3])))))]
    #[case::swap_binds_tighter_than_apply("~1^2134", StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::apply(StereoExpr::swap(StereoExpr::Lit(1)), Permutation::from_image(4, &[1, 0, 2, 3])))))]
    #[case::whitespace_ignored("  ?o :: { 1 , 2 }", StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::VarDomain("o".to_string(), vec![1, 2]))))]
    fn test_stereo_config(#[case] input: &str, #[case] expected: StereoConfigurationAst) {
        assert_eq!(stereo_config.parse(input).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(StereoConfigurationAst::Undetermined, "*")]
    #[case::not_stereo(StereoConfigurationAst::NotStereo, "!")]
    #[case::stereogenic(StereoConfigurationAst::Stereo(StereoIndexAst::Undetermined), "+")]
    #[case::lit(StereoConfigurationAst::from(1_u32), "1")]
    #[case::var(StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::Var("o".to_string()))), "?o")]
    #[case::lit_set(StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::LitSet(vec![1, 2]))), "{1,2}")]
    #[case::var_domain(StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::VarDomain("o".to_string(), vec![1, 2]))), "?o :: {1,2}")]
    #[case::swap(StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::swap(StereoExpr::Lit(1)))), "~1")]
    #[case::apply(StereoConfigurationAst::Stereo(StereoIndexAst::expr(StereoExpr::apply(StereoExpr::Lit(1), Permutation::from_image(4, &[1, 0, 2, 3])))), "1^2134")]
    fn test_fmt_stereo_config(#[case] c: StereoConfigurationAst, #[case] expected: &str) {
        struct W<'a>(&'a StereoConfigurationAst);
        impl fmt::Display for W<'_> {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                fmt_stereo_config(f, self.0)
            }
        }
        assert_eq!(W(&c).to_string(), expected);
    }

    #[rstest]
    #[case::tetrahedral("Th", StereoKind::Tetrahedral)]
    #[case::cis_trans("Ct", StereoKind::CisTrans)]
    #[case::square_planar("Sp", StereoKind::SquarePlanar)]
    #[case::trigonal_bipyramidal("Tb", StereoKind::TrigonalBipyramidal)]
    #[case::octahedral("Oh", StereoKind::Octahedral)]
    fn test_class(#[case] input: &str, #[case] expected: StereoKind) {
        assert_eq!(stereo_kind.parse(input).unwrap(), expected);
    }
}
