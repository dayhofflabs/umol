//! Stereo config-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnStreamDeserializer, FromEdn, ToEdn};
use umol_perm::Permutation;
use winnow::ascii::{digit1, multispace0};
use winnow::combinator::{alt, delimited, opt, preceded, separated, terminated};
use winnow::error::ErrMode;
use winnow::Parser;

use super::atom::single_key_map;
use super::error::{PResult, ParseError};
use super::value::{id, terminator};
use crate::ast::constraint::{StereoAtomConstraint, StereoBondConstraint};
use crate::ast::stereo::{
    StereoAtomAst, StereoBondAst, StereoConfigurationAst, StereoCosetAst, StereoExpr, StereoKind,
};
use crate::ast::traits::{FromAst, IntoAst};
use crate::dsl::config::{StereoAtomDefaults, StereoBondDefaults};

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

impl FromStr for StereoAtomAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(StereoAtomDsl::from_str(s)?.into_ast(&StereoAtomDefaults::default()))
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
    match (ast.kind, &ast.coset) {
        (StereoKind::Tetrahedral, &StereoCosetAst::Lit(1)) => Some("ccw"),
        (StereoKind::Tetrahedral, &StereoCosetAst::Lit(2)) => Some("cw"),
        _ => None,
    }
}

impl FromAst<StereoAtomAst> for StereoAtomDsl {
    type Ctx = StereoAtomDefaults;

    fn from_ast(ast: &StereoAtomAst, _ctx: &Self::Ctx) -> Self {
        let ast = ast.clone();
        StereoAtomDsl(ast.clone())
    }
}

impl IntoAst<StereoAtomAst> for StereoAtomDsl {
    type Ctx = StereoAtomDefaults;

    fn into_ast(self, _ctx: &Self::Ctx) -> StereoAtomAst {
        self.0
    }
}

pub fn parse_stereo_atom(input: &str) -> Result<StereoAtomDsl, ParseError> {
    stereo_atom.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn stereo_atom(i: &mut &str) -> PResult<StereoAtomDsl> {
    let kind = delimited(multispace0, stereo_kind, multispace0).parse_next(i)?;
    let coset = terminated(stereo_coset, multispace0).parse_next(i)?;
    Ok(StereoAtomDsl(StereoAtomAst::new(kind, coset)))
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

impl FromStr for StereoBondAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(StereoBondDsl::from_str(s)?.into_ast(&StereoBondDefaults::default()))
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
    match (ast.kind, &ast.coset) {
        (StereoKind::CisTrans, &StereoCosetAst::Lit(1)) => Some("z"),
        (StereoKind::CisTrans, &StereoCosetAst::Lit(2)) => Some("e"),
        _ => None,
    }
}

impl FromAst<StereoBondAst> for StereoBondDsl {
    type Ctx = StereoBondDefaults;

    fn from_ast(ast: &StereoBondAst, _ctx: &Self::Ctx) -> Self {
        StereoBondDsl(ast.clone())
    }
}

impl IntoAst<StereoBondAst> for StereoBondDsl {
    type Ctx = StereoBondDefaults;

    fn into_ast(self, _ctx: &Self::Ctx) -> StereoBondAst {
        self.0
    }
}

pub fn parse_stereo_bond(input: &str) -> Result<StereoBondDsl, ParseError> {
    stereo_bond.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn stereo_bond(i: &mut &str) -> PResult<StereoBondDsl> {
    let kind = delimited(multispace0, stereo_kind, multispace0).parse_next(i)?;
    let coset = terminated(stereo_coset, multispace0).parse_next(i)?;
    Ok(StereoBondDsl(StereoBondAst::new(kind, coset)))
}

/// `Th` / `Ct` / `Sp` / `Tb` / `Oh` symbol → `StereoKind`.
pub(crate) fn stereo_kind(i: &mut &str) -> PResult<StereoKind> {
    alt((
        "Th".value(StereoKind::Tetrahedral),
        "Ct".value(StereoKind::CisTrans),
        "Ax".value(StereoKind::Axial),
        "Sp".value(StereoKind::SquarePlanar),
        "Tb".value(StereoKind::TrigonalBipyramidal),
        "Oh".value(StereoKind::Octahedral),
    ))
    .parse_next(i)
}

/// Parse the `config` grammar into `StereoConfigurationAst` for constraints.
pub(crate) fn stereo_config(i: &mut &str) -> PResult<StereoConfigurationAst> {
    alt((
        '*'.value(StereoConfigurationAst::Undetermined),
        '!'.value(StereoConfigurationAst::NotStereo),
        '+'.value(StereoConfigurationAst::Stereo(StereoCosetAst::Undetermined)),
        terminated(stereo_lit, (multispace0, terminator))
            .map(|n| StereoConfigurationAst::Stereo(StereoCosetAst::Lit(n))),
        stereo_expr.map(|e| StereoConfigurationAst::Stereo(StereoCosetAst::Expr(Box::new(e)))),
    ))
    .parse_next(i)
}

pub(crate) fn parse_stereo_coset(input: &str) -> Result<StereoCosetAst, ParseError> {
    stereo_coset.parse(input).map_err(|e| e.into_inner())
}

/// Parse the `coset` grammar into `StereoCosetAst` (used in stereo elements).
fn stereo_coset(i: &mut &str) -> PResult<StereoCosetAst> {
    alt((
        '*'.value(StereoCosetAst::Undetermined),
        terminated(stereo_lit, (multispace0, terminator)).map(StereoCosetAst::Lit),
        stereo_expr.map(|e| StereoCosetAst::Expr(Box::new(e))),
    ))
    .parse_next(i)
}

/// `stereo-expr`: a `~`-prefixed base carrying zero or more left-associative
/// `^image` postfixes.
fn stereo_expr(i: &mut &str) -> PResult<StereoExpr> {
    let mut e = stereo_prefix(i)?;
    loop {
        multispace0.parse_next(i)?;
        if opt('^').parse_next(i)?.is_some() {
            e = StereoExpr::apply(e, perm_image(i)?);
        } else {
            return Ok(e);
        }
    }
}

/// `('~' | '\'') prefix-term | base` — unary `~` (swap) and `'` (mirror) bind
/// tighter than `^`.
fn stereo_prefix(i: &mut &str) -> PResult<StereoExpr> {
    multispace0.parse_next(i)?;
    if opt('~').parse_next(i)?.is_some() {
        Ok(StereoExpr::swap(stereo_prefix(i)?))
    } else if opt('\'').parse_next(i)?.is_some() {
        Ok(StereoExpr::mirror(stereo_prefix(i)?))
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
    fmt_stereo_coset(f, &atom.coset)
}

/// Write the stereo bond DSL
pub(crate) fn fmt_stereo_bond(f: &mut fmt::Formatter<'_>, bond: &StereoBondAst) -> fmt::Result {
    fmt_stereo_kind(f, bond.kind)?;
    fmt_stereo_coset(f, &bond.coset)
}

/// Write the stereo kind for `kind`.
pub(crate) fn fmt_stereo_kind(f: &mut fmt::Formatter<'_>, kind: StereoKind) -> fmt::Result {
    match kind {
        StereoKind::Tetrahedral => f.write_str("Th"),
        StereoKind::CisTrans => f.write_str("Ct"),
        StereoKind::Axial => f.write_str("Ax"),
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
        StereoConfigurationAst::Stereo(StereoCosetAst::Undetermined) => write!(f, "+"),
        StereoConfigurationAst::Stereo(StereoCosetAst::Lit(n)) => write!(f, "{n}"),
        StereoConfigurationAst::Stereo(StereoCosetAst::Expr(e)) => fmt_stereo_expr(f, e),
    }
}

/// Write a `StereoCosetAst` for the element `:type` body: `*` (open coset), a
/// literal, or an operator expression. `fmt_stereo_config` reuses the literal
/// and expression arms but writes its own `+` for `Stereo(Undetermined)`.
fn fmt_stereo_coset(f: &mut fmt::Formatter<'_>, coset: &StereoCosetAst) -> fmt::Result {
    match coset {
        StereoCosetAst::Undetermined => write!(f, "*"),
        StereoCosetAst::Lit(n) => write!(f, "{n}"),
        StereoCosetAst::Expr(e) => fmt_stereo_expr(f, e),
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
        StereoExpr::MirrorOp(inner) => {
            write!(f, "'")?;
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

impl FromAst<StereoAtomConstraint> for StereoAtomConstraintDsl {
    type Ctx = ();

    fn from_ast(ast: &StereoAtomConstraint, _ctx: &Self::Ctx) -> Self {
        match *ast {}
    }
}

impl IntoAst<StereoAtomConstraint> for StereoAtomConstraintDsl {
    type Ctx = ();

    fn into_ast(self, _ctx: &Self::Ctx) -> StereoAtomConstraint {
        match self {}
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

impl FromAst<StereoBondConstraint> for StereoBondConstraintDsl {
    type Ctx = ();

    fn from_ast(ast: &StereoBondConstraint, _ctx: &Self::Ctx) -> Self {
        match *ast {}
    }
}

impl IntoAst<StereoBondConstraint> for StereoBondConstraintDsl {
    type Ctx = ();

    fn into_ast(self, _ctx: &Self::Ctx) -> StereoBondConstraint {
        match self {}
    }
}

pub(crate) fn coset_lit(n: i64) -> Result<u32, DeError> {
    u32::try_from(n).map_err(|_| DeError::OutOfRange {
        value: n.to_string(),
        target: "u32",
        path: Vec::new(),
    })
}

/// Surface DSL wrapper around `StereoCosetAst` — the coset value under
/// `:stereo`. EDN form: int (`Lit`), `:undetermined`, a vector of ints
/// (`Expr(LitSet)`), or a string carrying the operator-expression subgrammar.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoCosetDsl(pub StereoCosetAst);

impl FromAst<StereoCosetAst> for StereoCosetDsl {
    type Ctx = ();

    fn from_ast(ast: &StereoCosetAst, _ctx: &Self::Ctx) -> Self {
        Self(ast.clone())
    }
}

impl IntoAst<StereoCosetAst> for StereoCosetDsl {
    type Ctx = ();

    fn into_ast(self, _ctx: &Self::Ctx) -> StereoCosetAst {
        self.0
    }
}

impl Display for StereoCosetDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_stereo_coset(f, &self.0)
    }
}

impl<'de> FromEdn<'de> for StereoCosetDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let coset = match edn {
            Edn::Int(n) => StereoCosetAst::Lit(coset_lit(*n)?),
            Edn::Keyword(k) if k.name() == "undetermined" => StereoCosetAst::Undetermined,
            Edn::Vector(xs) => {
                let mut set = Vec::with_capacity(xs.len());
                for e in xs.iter() {
                    let Edn::Int(n) = e else {
                        return Err(DeError::TypeMismatch {
                            expected: "int (coset-set element)",
                            got: e.kind(),
                            path: Vec::new(),
                        });
                    };
                    set.push(coset_lit(*n)?);
                }
                StereoCosetAst::Expr(Box::new(StereoExpr::LitSet(set)))
            }
            Edn::Str(s) => {
                parse_stereo_coset(s).map_err(|e| DeError::subgrammar("stereo coset", e))?
            }
            other => {
                return Err(DeError::TypeMismatch {
                    expected: "coset (int, :undetermined, vector, or string)",
                    got: other.kind(),
                    path: Vec::new(),
                });
            }
        };
        Ok(Self(coset))
    }
}

impl ToEdn for StereoCosetDsl {
    fn to_edn(&self) -> Edn<'static> {
        match &self.0 {
            StereoCosetAst::Lit(n) => Edn::Int(*n as i64),
            StereoCosetAst::Undetermined => {
                Edn::Keyword(EdnKeyword::owned("undetermined".to_string()))
            }
            StereoCosetAst::Expr(e) => match e.as_ref() {
                StereoExpr::LitSet(set) => Edn::Vector(
                    set.iter()
                        .map(|n| Edn::Int(*n as i64))
                        .collect::<Vec<_>>()
                        .into(),
                ),
                _ => Edn::Str(Cow::Owned(self.to_string())),
            },
        }
    }
}

/// Surface DSL wrapper around `StereoConfigurationAst` — the value under a
/// `:tetrahedral-stereo` / `:cis-trans-stereo` key. EDN form: `:undetermined`,
/// `:not-stereo`, or `{:stereo <coset>}`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoConfigurationDsl(pub StereoConfigurationAst);

impl FromAst<StereoConfigurationAst> for StereoConfigurationDsl {
    type Ctx = ();

    fn from_ast(ast: &StereoConfigurationAst, _ctx: &Self::Ctx) -> Self {
        Self(ast.clone())
    }
}

impl IntoAst<StereoConfigurationAst> for StereoConfigurationDsl {
    type Ctx = ();

    fn into_ast(self, _ctx: &Self::Ctx) -> StereoConfigurationAst {
        self.0
    }
}

impl<'de> FromEdn<'de> for StereoConfigurationDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Keyword(k) if k.name() == "undetermined" => {
                Ok(Self(StereoConfigurationAst::Undetermined))
            }
            Edn::Keyword(k) if k.name() == "not-stereo" => {
                Ok(Self(StereoConfigurationAst::NotStereo))
            }
            Edn::Map(m) if m.len() == 1 => {
                let (k, v) = m.iter().next().unwrap();
                let Edn::Keyword(key) = k else {
                    return Err(DeError::TypeMismatch {
                        expected: "keyword key",
                        got: k.kind(),
                        path: Vec::new(),
                    });
                };
                match key.name() {
                    "stereo" => Ok(Self(StereoConfigurationAst::Stereo(
                        StereoCosetDsl::from_edn(v)?.into_ast(&()),
                    ))),
                    other => Err(DeError::UnknownField {
                        key: other.to_string(),
                        path: vec!["stereo-configuration".into()],
                    }),
                }
            }
            other => Err(DeError::TypeMismatch {
                expected: ":undetermined / :not-stereo / {:stereo <coset>}",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for StereoConfigurationDsl {
    fn to_edn(&self) -> Edn<'static> {
        match &self.0 {
            StereoConfigurationAst::Undetermined => {
                Edn::Keyword(EdnKeyword::owned("undetermined".to_string()))
            }
            StereoConfigurationAst::NotStereo => {
                Edn::Keyword(EdnKeyword::owned("not-stereo".to_string()))
            }
            StereoConfigurationAst::Stereo(coset) => {
                single_key_map("stereo", StereoCosetDsl::from_ast(coset, &()).to_edn())
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

    #[rustfmt::skip]
    #[rstest]
    #[case::tetrahedral_ccw("Th1", StereoAtomDsl(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(1))))]
    #[case::tetrahedral_cw("Th2", StereoAtomDsl(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(2))))]
    #[case::open("Th*", StereoAtomDsl(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined)))]
    #[case::square_planar("Sp3", StereoAtomDsl(StereoAtomAst::new(StereoKind::SquarePlanar, StereoCosetAst::Lit(3))))]
    #[case::octahedral("Oh6", StereoAtomDsl(StereoAtomAst::new(StereoKind::Octahedral, StereoCosetAst::Lit(6))))]
    fn test_parse_stereo_atom(#[case] input: &str, #[case] expected: StereoAtomDsl) {
        assert_eq!(parse_stereo_atom(input).unwrap(), expected);
    }

    #[rstest]
    #[case::not_stereo("Th!")]
    fn test_parse_stereo_atom_error(#[case] input: &str) {
        assert_eq!(parse_stereo_atom(input).unwrap_err(), ParseError::Syntax);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::tetrahedral_ccw(StereoAtomDsl(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(1))), "Th1")]
    #[case::open(StereoAtomDsl(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined)), "Th*")]
    #[case::square_planar(StereoAtomDsl(StereoAtomAst::new(StereoKind::SquarePlanar, StereoCosetAst::Lit(3))), "Sp3")]
    #[case::octahedral(StereoAtomDsl(StereoAtomAst::new(StereoKind::Octahedral, StereoCosetAst::Lit(6))), "Oh6")]
    fn test_fmt_stereo_atom(#[case] form: StereoAtomDsl, #[case] expected: &str) {
        assert_eq!(form.to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::string("\"Th1\"", StereoAtomDsl(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(1))))]
    #[case::keyword_ccw(":ccw", StereoAtomDsl(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(1))))]
    #[case::keyword_cw(":cw", StereoAtomDsl(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(2))))]
    #[case::string_square_planar("\"Sp3\"", StereoAtomDsl(StereoAtomAst::new(StereoKind::SquarePlanar, StereoCosetAst::Lit(3))))]
    fn test_stereo_atom_dsl_from_edn(#[case] input: &str, #[case] expected: StereoAtomDsl) {
        assert_eq!(StereoAtomDsl::from_edn(&read_string(input).unwrap()).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::unknown_keyword(":xyz", DeError::Custom("unknown stereo atom keyword :xyz".to_string()))]
    #[case::wrong_type("5", DeError::TypeMismatch { expected: "string or stereo atom keyword", got: "int", path: Vec::new() })]
    fn test_stereo_atom_dsl_from_edn_error(#[case] input: &str, #[case] expected: DeError) {
        assert_eq!(StereoAtomDsl::from_edn(&read_string(input).unwrap()).unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::canonical_ccw(StereoAtomDsl(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(1))), ":ccw")]
    #[case::canonical_cw(StereoAtomDsl(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(2))), ":cw")]
    #[case::open_string(StereoAtomDsl(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined)), "\"Th*\"")]
    #[case::non_tetrahedral_string(StereoAtomDsl(StereoAtomAst::new(StereoKind::SquarePlanar, StereoCosetAst::Lit(1))), "\"Sp1\"")]
    #[case::tetrahedral_three_string(StereoAtomDsl(StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Lit(3))), "\"Th3\"")]
    fn test_stereo_atom_dsl_to_edn(#[case] form: StereoAtomDsl, #[case] expected: &str) {
        assert_eq!(form.to_edn(), read_string(expected).unwrap());
    }

    #[rstest]
    fn test_stereo_atom_dsl_into_ast() {
        let ast = StereoAtomAst::new(StereoKind::Tetrahedral, StereoCosetAst::Undetermined);
        assert_eq!(
            StereoAtomDsl(ast.clone()).into_ast(&StereoAtomDefaults::default()),
            ast
        );
        assert_eq!(
            StereoAtomDsl::from_ast(&ast, &StereoAtomDefaults::default()),
            StereoAtomDsl(ast)
        );
    }

    #[rstest]
    #[case::ccw("ccw", Some("Th1"))]
    #[case::cw("cw", Some("Th2"))]
    #[case::unknown("xyz", None)]
    fn test_expand_stereo_atom_keyword(#[case] name: &str, #[case] expected: Option<&str>) {
        assert_eq!(expand_stereo_atom_keyword(name), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::cis_trans_z("Ct1", StereoBondDsl(StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(1))))]
    #[case::cis_trans_e("Ct2", StereoBondDsl(StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(2))))]
    #[case::open("Ct*", StereoBondDsl(StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Undetermined)))]
    fn test_parse_stereo_bond(#[case] input: &str, #[case] expected: StereoBondDsl) {
        assert_eq!(parse_stereo_bond(input).unwrap(), expected);
    }

    #[rstest]
    #[case::not_stereo("Ct!")]
    fn test_parse_stereo_bond_error(#[case] input: &str) {
        assert_eq!(parse_stereo_bond(input).unwrap_err(), ParseError::Syntax);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::cis_trans_z(StereoBondDsl(StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(1))), "Ct1")]
    #[case::open(StereoBondDsl(StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Undetermined)), "Ct*")]
    fn test_fmt_stereo_bond(#[case] form: StereoBondDsl, #[case] expected: &str) {
        assert_eq!(form.to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::string("\"Ct1\"", StereoBondDsl(StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(1))))]
    #[case::keyword_z(":z", StereoBondDsl(StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(1))))]
    #[case::keyword_e(":e", StereoBondDsl(StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(2))))]
    fn test_stereo_bond_dsl_from_edn(#[case] input: &str, #[case] expected: StereoBondDsl) {
        assert_eq!(StereoBondDsl::from_edn(&read_string(input).unwrap()).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::unknown_keyword(":xyz", DeError::Custom("unknown stereo bond keyword :xyz".to_string()))]
    #[case::wrong_type("5", DeError::TypeMismatch { expected: "string or stereo bond keyword", got: "int", path: Vec::new() })]
    fn test_stereo_bond_dsl_from_edn_error(#[case] input: &str, #[case] expected: DeError) {
        assert_eq!(StereoBondDsl::from_edn(&read_string(input).unwrap()).unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::canonical_z(StereoBondDsl(StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(1))), ":z")]
    #[case::canonical_e(StereoBondDsl(StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Lit(2))), ":e")]
    #[case::open_string(StereoBondDsl(StereoBondAst::new(StereoKind::CisTrans, StereoCosetAst::Undetermined)), "\"Ct*\"")]
    fn test_stereo_bond_dsl_to_edn(#[case] form: StereoBondDsl, #[case] expected: &str) {
        assert_eq!(form.to_edn(), read_string(expected).unwrap());
    }

    #[rstest]
    #[case::z("z", Some("Ct1"))]
    #[case::e("e", Some("Ct2"))]
    #[case::unknown("xyz", None)]
    fn test_expand_stereo_bond_keyword(#[case] name: &str, #[case] expected: Option<&str>) {
        assert_eq!(expand_stereo_bond_keyword(name), expected);
    }

    #[rstest]
    #[case::tetrahedral("Th", StereoKind::Tetrahedral)]
    #[case::cis_trans("Ct", StereoKind::CisTrans)]
    #[case::axial("Ax", StereoKind::Axial)]
    #[case::square_planar("Sp", StereoKind::SquarePlanar)]
    #[case::trigonal_bipyramidal("Tb", StereoKind::TrigonalBipyramidal)]
    #[case::octahedral("Oh", StereoKind::Octahedral)]
    fn test_stereo_kind(#[case] input: &str, #[case] expected: StereoKind) {
        assert_eq!(stereo_kind.parse(input).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined("*", StereoConfigurationAst::Undetermined)]
    #[case::not_stereo("!", StereoConfigurationAst::NotStereo)]
    #[case::stereogenic("+", StereoConfigurationAst::Stereo(StereoCosetAst::Undetermined))]
    #[case::lit("1", StereoConfigurationAst::from(1_u32))]
    #[case::var("?o", StereoConfigurationAst::Stereo(StereoCosetAst::expr(StereoExpr::Var("o".to_string()))))]
    #[case::lit_set("{1,2}", StereoConfigurationAst::Stereo(StereoCosetAst::expr(StereoExpr::LitSet(vec![1, 2]))))]
    #[case::var_domain("?o :: {1,2}", StereoConfigurationAst::Stereo(StereoCosetAst::expr(StereoExpr::VarDomain("o".to_string(), vec![1, 2]))))]
    #[case::swap("~1", StereoConfigurationAst::Stereo(StereoCosetAst::expr(StereoExpr::swap(StereoExpr::Lit(1)))))]
    #[case::mirror("'1", StereoConfigurationAst::Stereo(StereoCosetAst::expr(StereoExpr::mirror(StereoExpr::Lit(1)))))]
    #[case::apply("1^2134", StereoConfigurationAst::Stereo(StereoCosetAst::expr(StereoExpr::apply(StereoExpr::Lit(1), Permutation::from_image(4, &[1, 0, 2, 3])))))]
    #[case::swap_binds_tighter_than_apply("~1^2134", StereoConfigurationAst::Stereo(StereoCosetAst::expr(StereoExpr::apply(StereoExpr::swap(StereoExpr::Lit(1)), Permutation::from_image(4, &[1, 0, 2, 3])))))]
    #[case::mirror_binds_tighter_than_apply("'1^2134", StereoConfigurationAst::Stereo(StereoCosetAst::expr(StereoExpr::apply(StereoExpr::mirror(StereoExpr::Lit(1)), Permutation::from_image(4, &[1, 0, 2, 3])))))]
    #[case::whitespace_ignored("  ?o :: { 1 , 2 }", StereoConfigurationAst::Stereo(StereoCosetAst::expr(StereoExpr::VarDomain("o".to_string(), vec![1, 2]))))]
    fn test_stereo_config(#[case] input: &str, #[case] expected: StereoConfigurationAst) {
        assert_eq!(stereo_config.parse(input).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(StereoConfigurationAst::Undetermined, "*")]
    #[case::not_stereo(StereoConfigurationAst::NotStereo, "!")]
    #[case::stereogenic(StereoConfigurationAst::Stereo(StereoCosetAst::Undetermined), "+")]
    #[case::lit(StereoConfigurationAst::from(1_u32), "1")]
    #[case::var(StereoConfigurationAst::Stereo(StereoCosetAst::expr(StereoExpr::Var("o".to_string()))), "?o")]
    #[case::lit_set(StereoConfigurationAst::Stereo(StereoCosetAst::expr(StereoExpr::LitSet(vec![1, 2]))), "{1,2}")]
    #[case::var_domain(StereoConfigurationAst::Stereo(StereoCosetAst::expr(StereoExpr::VarDomain("o".to_string(), vec![1, 2]))), "?o :: {1,2}")]
    #[case::swap(StereoConfigurationAst::Stereo(StereoCosetAst::expr(StereoExpr::swap(StereoExpr::Lit(1)))), "~1")]
    #[case::apply(StereoConfigurationAst::Stereo(StereoCosetAst::expr(StereoExpr::apply(StereoExpr::Lit(1), Permutation::from_image(4, &[1, 0, 2, 3])))), "1^2134")]
    fn test_fmt_stereo_config(#[case] c: StereoConfigurationAst, #[case] expected: &str) {
        struct W<'a>(&'a StereoConfigurationAst);
        impl fmt::Display for W<'_> {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                fmt_stereo_config(f, self.0)
            }
        }
        assert_eq!(W(&c).to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(StereoCosetAst::Lit(2))]
    #[case::undetermined(StereoCosetAst::Undetermined)]
    #[case::expr_lit_set(StereoCosetAst::Expr(Box::new(StereoExpr::LitSet(vec![1, 2]))))]
    #[case::expr_swap(StereoCosetAst::Expr(Box::new(StereoExpr::swap(StereoExpr::Lit(1)))))]
    fn test_stereo_coset_dsl_into_ast(#[case] ast: StereoCosetAst) {
        assert_eq!(StereoCosetDsl(ast.clone()).into_ast(&()), ast);
        assert_eq!(StereoCosetDsl::from_ast(&ast, &()), StereoCosetDsl(ast));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::int("2", StereoCosetDsl(StereoCosetAst::Lit(2)))]
    #[case::undetermined(":undetermined", StereoCosetDsl(StereoCosetAst::Undetermined))]
    #[case::vector("[1 2]", StereoCosetDsl(StereoCosetAst::Expr(Box::new(StereoExpr::LitSet(vec![1, 2])))))]
    #[case::string_lit("\"3\"", StereoCosetDsl(StereoCosetAst::Lit(3)))]
    #[case::string_expr("\"~1\"", StereoCosetDsl(StereoCosetAst::Expr(Box::new(StereoExpr::swap(StereoExpr::Lit(1))))))]
    fn test_stereo_coset_dsl_from_edn(#[case] input: &str, #[case] expected: StereoCosetDsl) {
        assert_eq!(StereoCosetDsl::from_edn(&read_string(input).unwrap()).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wrong_type("nil", DeError::TypeMismatch { expected: "coset (int, :undetermined, vector, or string)", got: "nil", path: Vec::new() })]
    #[case::non_int_in_vector("[1 nil]", DeError::TypeMismatch { expected: "int (coset-set element)", got: "nil", path: Vec::new() })]
    #[case::negative("-1", DeError::OutOfRange { value: "-1".to_string(), target: "u32", path: Vec::new() })]
    fn test_stereo_coset_dsl_from_edn_error(#[case] input: &str, #[case] expected: DeError) {
        assert_eq!(StereoCosetDsl::from_edn(&read_string(input).unwrap()).unwrap_err(), expected);
    }

    #[rstest]
    fn test_stereo_coset_dsl_from_edn_rejects_invalid_string() {
        let err = StereoCosetDsl::from_edn(&read_string("\"???\"").unwrap()).unwrap_err();
        assert!(matches!(err, DeError::Subgrammar { grammar: "stereo coset", .. }));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(StereoCosetDsl(StereoCosetAst::Lit(2)), "2")]
    #[case::undetermined(StereoCosetDsl(StereoCosetAst::Undetermined), ":undetermined")]
    #[case::expr_lit_set(StereoCosetDsl(StereoCosetAst::Expr(Box::new(StereoExpr::LitSet(vec![1, 2])))), "[1 2]")]
    #[case::expr_swap(StereoCosetDsl(StereoCosetAst::Expr(Box::new(StereoExpr::swap(StereoExpr::Lit(1))))), "\"~1\"")]
    fn test_stereo_coset_dsl_to_edn(#[case] form: StereoCosetDsl, #[case] expected: &str) {
        assert_eq!(form.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(StereoConfigurationAst::Undetermined)]
    #[case::not_stereo(StereoConfigurationAst::NotStereo)]
    #[case::stereo_lit(StereoConfigurationAst::Stereo(StereoCosetAst::Lit(1)))]
    fn test_stereo_configuration_dsl_into_ast(#[case] ast: StereoConfigurationAst) {
        assert_eq!(StereoConfigurationDsl(ast.clone()).into_ast(&()), ast);
        assert_eq!(StereoConfigurationDsl::from_ast(&ast, &()), StereoConfigurationDsl(ast));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(":undetermined", StereoConfigurationDsl(StereoConfigurationAst::Undetermined))]
    #[case::not_stereo(":not-stereo", StereoConfigurationDsl(StereoConfigurationAst::NotStereo))]
    #[case::stereo_lit("{:stereo 1}", StereoConfigurationDsl(StereoConfigurationAst::Stereo(StereoCosetAst::Lit(1))))]
    #[case::stereo_undetermined("{:stereo :undetermined}", StereoConfigurationDsl(StereoConfigurationAst::Stereo(StereoCosetAst::Undetermined)))]
    #[case::stereo_set("{:stereo [1 2]}", StereoConfigurationDsl(StereoConfigurationAst::Stereo(StereoCosetAst::Expr(Box::new(StereoExpr::LitSet(vec![1, 2]))))))]
    fn test_stereo_configuration_dsl_from_edn(#[case] input: &str, #[case] expected: StereoConfigurationDsl) {
        assert_eq!(StereoConfigurationDsl::from_edn(&read_string(input).unwrap()).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::unknown_keyword(":bogus", DeError::TypeMismatch { expected: ":undetermined / :not-stereo / {:stereo <coset>}", got: "keyword", path: Vec::new() })]
    #[case::unknown_key("{:bogus 1}", DeError::UnknownField { key: "bogus".to_string(), path: vec!["stereo-configuration".into()] })]
    #[case::wrong_type("1", DeError::TypeMismatch { expected: ":undetermined / :not-stereo / {:stereo <coset>}", got: "int", path: Vec::new() })]
    fn test_stereo_configuration_dsl_from_edn_error(#[case] input: &str, #[case] expected: DeError) {
        assert_eq!(StereoConfigurationDsl::from_edn(&read_string(input).unwrap()).unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(StereoConfigurationDsl(StereoConfigurationAst::Undetermined), ":undetermined")]
    #[case::not_stereo(StereoConfigurationDsl(StereoConfigurationAst::NotStereo), ":not-stereo")]
    #[case::stereo_lit(StereoConfigurationDsl(StereoConfigurationAst::Stereo(StereoCosetAst::Lit(1))), "{:stereo 1}")]
    #[case::stereo_undetermined(StereoConfigurationDsl(StereoConfigurationAst::Stereo(StereoCosetAst::Undetermined)), "{:stereo :undetermined}")]
    fn test_stereo_configuration_dsl_to_edn(#[case] form: StereoConfigurationDsl, #[case] expected: &str) {
        assert_eq!(form.to_edn(), read_string(expected).unwrap());
    }
}
