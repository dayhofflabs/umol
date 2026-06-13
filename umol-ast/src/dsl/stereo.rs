//! Stereo config-string DSL.

use std::borrow::Cow;
use std::collections::BTreeSet;
use std::fmt::{self, Display};
use std::str::FromStr;

use strum::VariantArray;
use umol_edn::{
    DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnSet, EdnStreamDeserializer, FromEdn, ToEdn,
};
use umol_perm::{Orientation, Permutation};
use winnow::ascii::{digit1, multispace0};
use winnow::combinator::{
    alt, delimited, opt, preceded, repeat, separated, separated_pair, terminated,
};
use winnow::error::ErrMode;
use winnow::token::one_of;
use winnow::Parser;

use super::atom::single_key_map;
use super::config::{StereoAtomDefaults, StereoBondDefaults};
use super::error::{PResult, ParseError};
use super::value::{id, terminator};
use crate::ast::constraint::{
    FluxionalityAst, LigandPairAst, LigandSymmetryAst, OrientedPermutationAst, PermutationAst,
    StereoAtomConstraint, StereoBondConstraint, StereogenicityAst, StereogenicityRelationAst,
    TopicityAst, TopicityRelationAst,
};
use crate::ast::ids::StereoLigandId;
use crate::ast::operators::MemOp;
use crate::ast::stereo::{
    StereoAtomAst, StereoBondAst, StereoConfigurationAst, StereoCosetAst, StereoExpr, StereoKind,
    Stereogenicity, Topicity,
};
use crate::ast::traits::{FromAst, IntoAst};

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
    if !ast.constraints.is_empty() {
        return None;
    }
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
    let constraints: Vec<StereoAtomConstraint> =
        repeat(0.., move |i: &mut &str| stereo_atom_predicate(i, kind)).parse_next(i)?;
    Ok(StereoAtomDsl(
        StereoAtomAst::new(kind, coset).with_constraints(constraints),
    ))
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
    if !ast.constraints.is_empty() {
        return None;
    }
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
    let constraints: Vec<StereoBondConstraint> =
        repeat(0.., move |i: &mut &str| stereo_bond_predicate(i, kind)).parse_next(i)?;
    Ok(StereoBondDsl(
        StereoBondAst::new(kind, coset).with_constraints(constraints),
    ))
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

fn cycle_point(i: &mut &str) -> PResult<usize> {
    let s: &str = digit1.parse_next(i)?;
    s.parse::<usize>()
        .map_err(|_| ErrMode::Backtrack(ParseError::Syntax))
}

/// One cycle `(p0,p1,…)` mapping `p0→p1→…→p0`; the empty `()` is the identity.
fn cycle(i: &mut &str) -> PResult<Vec<usize>> {
    delimited('(', separated(0.., cycle_point, ','), ')').parse_next(i)
}

/// A permutation in disjoint-cycle notation → `Permutation` of `degree`
/// (0-indexed); validated as in-range and disjoint, since `from_cycles` panics
/// on a non-bijection.
fn perm_cycles(i: &mut &str, degree: usize) -> PResult<Permutation> {
    let cycles: Vec<Vec<usize>> = repeat(1.., cycle).parse_next(i)?;
    let mut seen = vec![false; degree];
    for cycle in &cycles {
        for &p in cycle {
            if p >= degree || seen[p] {
                return Err(ErrMode::Cut(ParseError::Syntax));
            }
            seen[p] = true;
        }
    }
    Ok(Permutation::from_cycles(degree, &cycles))
}

/// `~` (the kind involution, eager) or an explicit permutation in cycle notation.
fn stereo_perm(i: &mut &str, kind: StereoKind) -> PResult<Permutation> {
    if opt('~').parse_next(i)?.is_some() {
        Ok(kind.involution())
    } else {
        perm_cycles(i, kind.degree())
    }
}

fn ligand_pair(i: &mut &str) -> PResult<LigandPairAst> {
    let (a, b) =
        delimited('(', separated_pair(cycle_point, ',', cycle_point), ')').parse_next(i)?;
    Ok(LigandPairAst::new(
        StereoLigandId(a as u8),
        StereoLigandId(b as u8),
    ))
}

fn topicity_glyph(i: &mut &str) -> PResult<Topicity> {
    alt((
        '='.value(Topicity::Homotopic),
        '\''.value(Topicity::Enantiotopic),
        '/'.value(Topicity::Diastereotopic),
    ))
    .parse_next(i)
}

fn stereogenicity_glyph(i: &mut &str) -> PResult<Stereogenicity> {
    alt((
        '='.value(Stereogenicity::Symmetric),
        '\''.value(Stereogenicity::Prochiral),
        '/'.value(Stereogenicity::Stereogenic),
    ))
    .parse_next(i)
}

fn topicity_relation_inline(i: &mut &str) -> PResult<TopicityRelationAst> {
    if opt('*').parse_next(i)?.is_some() {
        return Ok(TopicityRelationAst::Undetermined);
    }
    let neg = opt('!').parse_next(i)?.is_some();
    let v = topicity_glyph(i)?;
    Ok(if neg {
        TopicityRelationAst::NotSet(vec![v])
    } else {
        TopicityRelationAst::Lit(v)
    })
}

fn stereogenicity_relation_inline(i: &mut &str) -> PResult<StereogenicityRelationAst> {
    if opt('*').parse_next(i)?.is_some() {
        return Ok(StereogenicityRelationAst::Undetermined);
    }
    let neg = opt('!').parse_next(i)?.is_some();
    let v = stereogenicity_glyph(i)?;
    Ok(if neg {
        StereogenicityRelationAst::NotSet(vec![v])
    } else {
        StereogenicityRelationAst::Lit(v)
    })
}

/// Parse one `#p`/`#f`/`#o`/`#g` predicate for `kind`. Atom and bond share the
/// grammar; the macro emits a parser per constraint type.
macro_rules! stereo_predicate_parser {
    ($name:ident, $constraint:ident) => {
        fn $name(i: &mut &str, kind: StereoKind) -> PResult<$constraint> {
            '#'.parse_next(i)?;
            let tag = one_of(['p', 'f', 'o', 'g']).parse_next(i)?;
            match tag {
                'p' => {
                    let mem = if opt('!').parse_next(i)?.is_some() {
                        MemOp::NotIn
                    } else {
                        MemOp::In
                    };
                    let (perm, orientation) = if opt('~').parse_next(i)?.is_some() {
                        let orientation = if kind.is_chiral_class() {
                            Orientation::Improper
                        } else {
                            Orientation::Proper
                        };
                        (kind.involution(), orientation)
                    } else {
                        let orientation = if opt('\'').parse_next(i)?.is_some() {
                            Orientation::Improper
                        } else {
                            Orientation::Proper
                        };
                        (perm_cycles(i, kind.degree())?, orientation)
                    };
                    Ok($constraint::LigandSymmetry(LigandSymmetryAst {
                        perm: OrientedPermutationAst {
                            perm: PermutationAst(perm),
                            orientation,
                        },
                        mem,
                    }))
                }
                'f' => Ok($constraint::Fluxionality(FluxionalityAst {
                    perm: PermutationAst(stereo_perm(i, kind)?),
                })),
                'o' => {
                    let rel = topicity_relation_inline(i)?;
                    let pair = ligand_pair(i)?;
                    Ok($constraint::Topicity(TopicityAst { pair, rel }))
                }
                'g' => Ok($constraint::Stereogenicity(StereogenicityAst(
                    stereogenicity_relation_inline(i)?,
                ))),
                _ => unreachable!("one_of restricts the tag"),
            }
        }
    };
}

stereo_predicate_parser! { stereo_atom_predicate, StereoAtomConstraint }
stereo_predicate_parser! { stereo_bond_predicate, StereoBondConstraint }

/// `~` (when the perm is the kind involution) or an explicit permutation in
/// cycle notation (`Permutation`'s `Display`).
fn fmt_stereo_perm(f: &mut fmt::Formatter<'_>, perm: Permutation, kind: StereoKind) -> fmt::Result {
    if perm == kind.involution() {
        f.write_str("~")
    } else {
        write!(f, "{perm}")
    }
}

fn topicity_char(t: Topicity) -> char {
    match t {
        Topicity::Homotopic => '=',
        Topicity::Enantiotopic => '\'',
        Topicity::Diastereotopic => '/',
    }
}

fn stereogenicity_char(s: Stereogenicity) -> char {
    match s {
        Stereogenicity::Symmetric => '=',
        Stereogenicity::Prochiral => '\'',
        Stereogenicity::Stereogenic => '/',
    }
}

fn fmt_topicity_relation(f: &mut fmt::Formatter<'_>, rel: &TopicityRelationAst) -> fmt::Result {
    let members = rel.to_set();
    match members.len() {
        1 => write!(f, "{}", topicity_char(members.into_iter().next().unwrap())),
        2 => {
            let complement = Topicity::VARIANTS
                .iter()
                .copied()
                .find(|v| !members.contains(v))
                .unwrap();
            write!(f, "!{}", topicity_char(complement))
        }
        _ => f.write_str("*"),
    }
}

fn fmt_stereogenicity_relation(
    f: &mut fmt::Formatter<'_>,
    rel: &StereogenicityRelationAst,
) -> fmt::Result {
    let members = rel.to_set();
    match members.len() {
        1 => write!(
            f,
            "{}",
            stereogenicity_char(members.into_iter().next().unwrap())
        ),
        2 => {
            let complement = Stereogenicity::VARIANTS
                .iter()
                .copied()
                .find(|v| !members.contains(v))
                .unwrap();
            write!(f, "!{}", stereogenicity_char(complement))
        }
        _ => f.write_str("*"),
    }
}

/// Render one stereo constraint as its inline predicate. `~` is emitted for a
/// `#p`/`#f` perm equal to the kind involution (matching the involution's
/// orientation for `#p`); otherwise explicit cycles.
macro_rules! stereo_constraint_fmt {
    ($name:ident, $constraint:ident) => {
        fn $name(f: &mut fmt::Formatter<'_>, c: &$constraint, kind: StereoKind) -> fmt::Result {
            match c {
                $constraint::LigandSymmetry(ls) => {
                    f.write_str("#p")?;
                    if ls.mem == MemOp::NotIn {
                        f.write_str("!")?;
                    }
                    let involution_orientation = if kind.is_chiral_class() {
                        Orientation::Improper
                    } else {
                        Orientation::Proper
                    };
                    if ls.perm.perm.0 == kind.involution()
                        && ls.perm.orientation == involution_orientation
                    {
                        f.write_str("~")
                    } else {
                        if ls.perm.orientation == Orientation::Improper {
                            f.write_str("'")?;
                        }
                        write!(f, "{}", ls.perm.perm.0)
                    }
                }
                $constraint::Fluxionality(fx) => {
                    f.write_str("#f")?;
                    fmt_stereo_perm(f, fx.perm.0, kind)
                }
                $constraint::Topicity(t) => {
                    f.write_str("#o")?;
                    fmt_topicity_relation(f, &t.rel)?;
                    write!(f, "({},{})", t.pair.first().0, t.pair.second().0)
                }
                $constraint::Stereogenicity(g) => {
                    f.write_str("#g")?;
                    fmt_stereogenicity_relation(f, &g.0)
                }
            }
        }
    };
}

stereo_constraint_fmt! { fmt_stereo_atom_constraint, StereoAtomConstraint }
stereo_constraint_fmt! { fmt_stereo_bond_constraint, StereoBondConstraint }

/// Write the stereo atom DSL
pub(crate) fn fmt_stereo_atom(f: &mut fmt::Formatter<'_>, atom: &StereoAtomAst) -> fmt::Result {
    fmt_stereo_kind(f, atom.kind)?;
    fmt_stereo_coset(f, &atom.coset)?;
    // Vacuous (Undetermined) predicates are elided on canonical render, as for
    // the atom-string (§7.1); they remain admissible on parse.
    for c in atom.constraints.iter().filter(|c| !c.is_undetermined()) {
        fmt_stereo_atom_constraint(f, c, atom.kind)?;
    }
    Ok(())
}

/// Write the stereo bond DSL
pub(crate) fn fmt_stereo_bond(f: &mut fmt::Formatter<'_>, bond: &StereoBondAst) -> fmt::Result {
    fmt_stereo_kind(f, bond.kind)?;
    fmt_stereo_coset(f, &bond.coset)?;
    for c in bond.constraints.iter().filter(|c| !c.is_undetermined()) {
        fmt_stereo_bond_constraint(f, c, bond.kind)?;
    }
    Ok(())
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

/// A permutation as a vector of disjoint cycles (`[[0 1 2] [3 4]]`; identity `[]`).
/// The degree is not encoded — fixed points drop out — so the reader supplies it
/// from the stereo kind.
fn perm_to_vov(perm: Permutation) -> Edn<'static> {
    let cycles: Vec<Edn<'static>> = perm
        .cycles()
        .into_iter()
        .map(|cycle| {
            Edn::Vector(
                cycle
                    .into_iter()
                    .map(|p| Edn::Int(p as i64))
                    .collect::<Vec<_>>()
                    .into(),
            )
        })
        .collect();
    Edn::Vector(cycles.into())
}

fn perm_from_vov(edn: &Edn, degree: usize) -> Result<Permutation, DeError> {
    let Edn::Vector(cycles_edn) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "vector of cycles",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    let mut seen = vec![false; degree];
    let mut cycles: Vec<Vec<usize>> = Vec::with_capacity(cycles_edn.len());
    for cycle_edn in cycles_edn.iter() {
        let Edn::Vector(points) = cycle_edn else {
            return Err(DeError::TypeMismatch {
                expected: "cycle (vector of ints)",
                got: cycle_edn.kind(),
                path: Vec::new(),
            });
        };
        let mut cycle = Vec::with_capacity(points.len());
        for p in points.iter() {
            let Edn::Int(n) = p else {
                return Err(DeError::TypeMismatch {
                    expected: "int (cycle point)",
                    got: p.kind(),
                    path: Vec::new(),
                });
            };
            let point = usize::try_from(*n)
                .ok()
                .filter(|&x| x < degree && !seen[x])
                .ok_or_else(|| DeError::OutOfRange {
                    value: n.to_string(),
                    target: "ligand position",
                    path: Vec::new(),
                })?;
            seen[point] = true;
            cycle.push(point);
        }
        cycles.push(cycle);
    }
    Ok(Permutation::from_cycles(degree, &cycles))
}

/// Generates the structured-EDN serialization/deserialization for a relation type:
/// `:undetermined` (the full domain), a single keyword (singleton), or a keyword set `#{…}`.
/// Keyword names map 1:1 to the domain variants per the table; the AST itself carries none.
macro_rules! relation_serde {
    ($to:ident, $from:ident, $relation:ident, $domain:ty, $($variant:path => $kw:literal),+ $(,)?) => {
        fn $to(rel: &$relation) -> Edn<'static> {
            let members = rel.to_set();
            if members.len() == <$domain as VariantArray>::VARIANTS.len() {
                Edn::Keyword(EdnKeyword::owned("undetermined".to_string()))
            } else if members.len() == 1 {
                let kw = match members.into_iter().next().unwrap() {
                    $($variant => $kw,)+
                };
                Edn::Keyword(EdnKeyword::owned(kw.to_string()))
            } else {
                let set: EdnSet<'static> = members
                    .into_iter()
                    .map(|m| {
                        let kw = match m {
                            $($variant => $kw,)+
                        };
                        Edn::Keyword(EdnKeyword::owned(kw.to_string()))
                    })
                    .collect();
                Edn::Set(set)
            }
        }

        fn $from(edn: &Edn) -> Result<$relation, DeError> {
            fn keyword_member(name: &str) -> Option<$domain> {
                match name {
                    $($kw => Some($variant),)+
                    _ => None,
                }
            }
            match edn {
                Edn::Keyword(k) if k.name() == "undetermined" => Ok($relation::Undetermined),
                Edn::Keyword(k) => {
                    let m = keyword_member(k.name()).ok_or_else(|| DeError::TypeMismatch {
                        expected: concat!(stringify!($relation), " keyword"),
                        got: edn.kind(),
                        path: Vec::new(),
                    })?;
                    Ok($relation::from_set(BTreeSet::from([m])).unwrap())
                }
                Edn::Set(s) => {
                    let mut members = BTreeSet::new();
                    for e in s.iter() {
                        let Edn::Keyword(k) = e else {
                            return Err(DeError::TypeMismatch {
                                expected: "relation keyword",
                                got: e.kind(),
                                path: Vec::new(),
                            });
                        };
                        let m = keyword_member(k.name()).ok_or_else(|| DeError::TypeMismatch {
                            expected: "relation keyword",
                            got: e.kind(),
                            path: Vec::new(),
                        })?;
                        members.insert(m);
                    }
                    $relation::from_set(members).ok_or_else(|| DeError::TypeMismatch {
                        expected: "non-empty relation set",
                        got: edn.kind(),
                        path: Vec::new(),
                    })
                }
                other => Err(DeError::TypeMismatch {
                    expected: "relation (keyword, set, or :undetermined)",
                    got: other.kind(),
                    path: Vec::new(),
                }),
            }
        }
    };
}

relation_serde! {
    topicity_relation_to_edn, topicity_relation_from_edn, TopicityRelationAst, Topicity,
    Topicity::Homotopic => "homotopic",
    Topicity::Enantiotopic => "enantiotopic",
    Topicity::Diastereotopic => "diastereotopic",
}

relation_serde! {
    stereogenicity_relation_to_edn, stereogenicity_relation_from_edn,
    StereogenicityRelationAst, Stereogenicity,
    Stereogenicity::Symmetric => "symmetric",
    Stereogenicity::Prochiral => "prochiral",
    Stereogenicity::Stereogenic => "stereogenic",
}

/// `StereoKind` ↔ kebab keyword (`:tetrahedral`, `:cis-trans`, …).
pub(crate) fn stereo_kind_to_edn(kind: StereoKind) -> Edn<'static> {
    let name = match kind {
        StereoKind::Tetrahedral => "tetrahedral",
        StereoKind::CisTrans => "cis-trans",
        StereoKind::Axial => "axial",
        StereoKind::SquarePlanar => "square-planar",
        StereoKind::TrigonalBipyramidal => "trigonal-bipyramidal",
        StereoKind::Octahedral => "octahedral",
    };
    Edn::keyword(name)
}

pub(crate) fn stereo_kind_from_edn(edn: &Edn) -> Result<StereoKind, DeError> {
    let Edn::Keyword(k) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "stereo-kind keyword",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    match k.name() {
        "tetrahedral" => Ok(StereoKind::Tetrahedral),
        "cis-trans" => Ok(StereoKind::CisTrans),
        "axial" => Ok(StereoKind::Axial),
        "square-planar" => Ok(StereoKind::SquarePlanar),
        "trigonal-bipyramidal" => Ok(StereoKind::TrigonalBipyramidal),
        "octahedral" => Ok(StereoKind::Octahedral),
        other => Err(DeError::Custom(format!("unknown stereo kind :{other}"))),
    }
}

fn ligand_position(edn: &Edn) -> Result<StereoLigandId, DeError> {
    let Edn::Int(n) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "int (ligand position)",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    let v = u8::try_from(*n).map_err(|_| DeError::OutOfRange {
        value: n.to_string(),
        target: "ligand position",
        path: Vec::new(),
    })?;
    Ok(StereoLigandId(v))
}

fn ligand_symmetry_to_edn(ls: &LigandSymmetryAst) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(3);
    m.insert(Edn::keyword("perm"), perm_to_vov(ls.perm.perm.0));
    if ls.perm.orientation == Orientation::Improper {
        m.insert(Edn::keyword("orientation"), Edn::keyword("improper"));
    }
    if ls.mem == MemOp::NotIn {
        m.insert(Edn::keyword("member"), Edn::keyword("not-in"));
    }
    Edn::Map(m)
}

fn ligand_symmetry_from_edn(edn: &Edn, kind: StereoKind) -> Result<LigandSymmetryAst, DeError> {
    let Edn::Map(m) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "ligand-symmetry map",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    let perm_edn = m.get_keyword("perm").ok_or_else(|| DeError::MissingField {
        key: "perm".into(),
        path: vec!["ligand-symmetry".into()],
    })?;
    let perm = PermutationAst(perm_from_vov(perm_edn, kind.degree())?);
    let orientation = match m.get_keyword("orientation") {
        None => Orientation::Proper,
        Some(Edn::Keyword(k)) if k.name() == "proper" => Orientation::Proper,
        Some(Edn::Keyword(k)) if k.name() == "improper" => Orientation::Improper,
        Some(other) => {
            return Err(DeError::TypeMismatch {
                expected: ":proper | :improper",
                got: other.kind(),
                path: vec!["ligand-symmetry".into()],
            })
        }
    };
    let mem = match m.get_keyword("member") {
        None => MemOp::In,
        Some(Edn::Keyword(k)) if k.name() == "in" => MemOp::In,
        Some(Edn::Keyword(k)) if k.name() == "not-in" => MemOp::NotIn,
        Some(other) => {
            return Err(DeError::TypeMismatch {
                expected: ":in | :not-in",
                got: other.kind(),
                path: vec!["ligand-symmetry".into()],
            })
        }
    };
    Ok(LigandSymmetryAst {
        perm: OrientedPermutationAst { perm, orientation },
        mem,
    })
}

fn topicity_to_edn(t: &TopicityAst) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(2);
    m.insert(
        Edn::keyword("pair"),
        Edn::Vector(
            vec![
                Edn::Int(t.pair.first().0 as i64),
                Edn::Int(t.pair.second().0 as i64),
            ]
            .into(),
        ),
    );
    m.insert(Edn::keyword("relation"), topicity_relation_to_edn(&t.rel));
    Edn::Map(m)
}

fn topicity_from_edn(edn: &Edn) -> Result<TopicityAst, DeError> {
    let Edn::Map(m) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "topicity map",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    let pair_edn = m.get_keyword("pair").ok_or_else(|| DeError::MissingField {
        key: "pair".into(),
        path: vec!["topicity".into()],
    })?;
    let Edn::Vector(p) = pair_edn else {
        return Err(DeError::TypeMismatch {
            expected: "[i j] (ligand pair)",
            got: pair_edn.kind(),
            path: vec!["topicity".into()],
        });
    };
    if p.len() != 2 {
        return Err(DeError::Custom(
            "topicity pair must have 2 positions".into(),
        ));
    }
    let pair = LigandPairAst::new(ligand_position(&p[0])?, ligand_position(&p[1])?);
    let rel_edn = m
        .get_keyword("relation")
        .ok_or_else(|| DeError::MissingField {
            key: "relation".into(),
            path: vec!["topicity".into()],
        })?;
    Ok(TopicityAst {
        pair,
        rel: topicity_relation_from_edn(rel_edn)?,
    })
}

/// Render/parse a stereo constraint as its keyword tag plus value (the single
/// constraint entry inside the kind-bearing map). Atom and bond share the four
/// variants, so the codec is generated for each.
macro_rules! stereo_constraint_entry_codec {
    ($to:ident, $from:ident, $constraint:ident) => {
        fn $to(c: &$constraint) -> (&'static str, Edn<'static>) {
            match c {
                $constraint::LigandSymmetry(ls) => ("ligand-symmetry", ligand_symmetry_to_edn(ls)),
                $constraint::Fluxionality(f) => ("fluxionality", perm_to_vov(f.perm.0)),
                $constraint::Topicity(t) => ("topicity", topicity_to_edn(t)),
                $constraint::Stereogenicity(g) => {
                    ("stereogenicity", stereogenicity_relation_to_edn(&g.0))
                }
            }
        }

        fn $from(key: &str, value: &Edn, kind: StereoKind) -> Result<$constraint, DeError> {
            match key {
                "ligand-symmetry" => Ok($constraint::LigandSymmetry(ligand_symmetry_from_edn(
                    value, kind,
                )?)),
                "fluxionality" => Ok($constraint::Fluxionality(FluxionalityAst {
                    perm: PermutationAst(perm_from_vov(value, kind.degree())?),
                })),
                "topicity" => Ok($constraint::Topicity(topicity_from_edn(value)?)),
                "stereogenicity" => Ok($constraint::Stereogenicity(StereogenicityAst(
                    stereogenicity_relation_from_edn(value)?,
                ))),
                other => Err(DeError::Custom(format!(
                    "unknown stereo constraint keyword :{other}"
                ))),
            }
        }
    };
}

stereo_constraint_entry_codec! {
    stereo_atom_constraint_entry, stereo_atom_constraint_from_entry, StereoAtomConstraint
}
stereo_constraint_entry_codec! {
    stereo_bond_constraint_entry, stereo_bond_constraint_from_entry, StereoBondConstraint
}

/// Molecule-scope DSL wrapper for a stereo constraint. It carries the element
/// kind (the stereo subtype) so the permutation degree is known when parsing —
/// the EDN is a single map `{:kind <kind> <constraint-key> <value>}`, self-
/// contained, so the generic 2-field entity-leaf machinery applies.
macro_rules! stereo_constraint_dsl {
    ($dsl:ident, $constraint:ident, $entry:ident, $from_entry:ident, $context:literal) => {
        #[derive(Clone, Debug, PartialEq, Eq, Hash)]
        pub struct $dsl(pub StereoKind, pub $constraint);

        impl ToEdn for $dsl {
            fn to_edn(&self) -> Edn<'static> {
                let mut m = EdnMap::with_capacity(2);
                m.insert(Edn::keyword("kind"), stereo_kind_to_edn(self.0));
                let (key, value) = $entry(&self.1);
                m.insert(Edn::keyword(key), value);
                Edn::Map(m)
            }
        }

        impl<'de> FromEdn<'de> for $dsl {
            fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
                let Edn::Map(m) = edn else {
                    return Err(DeError::TypeMismatch {
                        expected: "stereo constraint map",
                        got: edn.kind(),
                        path: vec![$context.into()],
                    });
                };
                let kind = stereo_kind_from_edn(m.get_keyword("kind").ok_or_else(|| {
                    DeError::MissingField {
                        key: "kind".into(),
                        path: vec![$context.into()],
                    }
                })?)?;
                let mut entry = None;
                for (k, v) in m.iter() {
                    let Edn::Keyword(key) = k else {
                        return Err(DeError::TypeMismatch {
                            expected: "keyword key",
                            got: k.kind(),
                            path: vec![$context.into()],
                        });
                    };
                    if key.name() == "kind" {
                        continue;
                    }
                    if entry.is_some() {
                        return Err(DeError::Custom(format!(
                            "{} map has multiple constraint keys",
                            $context
                        )));
                    }
                    entry = Some((key.name().to_string(), v));
                }
                let (key, value) = entry.ok_or_else(|| {
                    DeError::Custom(format!("{} map missing the constraint key", $context))
                })?;
                Ok($dsl(kind, $from_entry(&key, value, kind)?))
            }
        }
    };
}

stereo_constraint_dsl! {
    StereoAtomConstraintDsl, StereoAtomConstraint,
    stereo_atom_constraint_entry, stereo_atom_constraint_from_entry, "stereo-atom-constraint"
}
stereo_constraint_dsl! {
    StereoBondConstraintDsl, StereoBondConstraint,
    stereo_bond_constraint_entry, stereo_bond_constraint_from_entry, "stereo-bond-constraint"
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
    #[case::fluxionality("Th2#f(0,1,2)")]
    #[case::ligand_symmetry_involution("Th2#p~")]
    #[case::ligand_symmetry_not_in("Th2#p!~")]
    #[case::ligand_symmetry_explicit("Th2#p(0,1,2)")]
    #[case::topicity("Th2#o=(0,1)")]
    #[case::topicity_negated("Th2#o!'(0,1)")]
    #[case::stereogenicity("Th2#g/")]
    #[case::multiple("Th2#f(0,1,2)#o=(0,1)#g/")]
    fn test_stereo_atom_inline_render_identity(#[case] s: &str) {
        assert_eq!(parse_stereo_atom(s).unwrap().to_string(), s);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::topicity_open("Th2#o*(0,1)", "Th2")]
    #[case::stereogenicity_open("Th2#g*", "Th2")]
    fn test_stereo_atom_inline_render(#[case] input: &str, #[case] canonical: &str) {
        assert_eq!(parse_stereo_atom(input).unwrap().to_string(), canonical);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fluxionality("Th2#f(0,1,2)",
        StereoAtomConstraint::Fluxionality(FluxionalityAst { perm: PermutationAst(Permutation::from_cycles(4, &[vec![0, 1, 2]])) }))]
    #[case::ligand_symmetry("Th2#p(0,1,2)",
        StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst {
            perm: OrientedPermutationAst { perm: PermutationAst(Permutation::from_cycles(4, &[vec![0, 1, 2]])), orientation: Orientation::Proper },
            mem: MemOp::In }))]
    #[case::topicity_negated("Th2#o!'(0,1)",
        StereoAtomConstraint::Topicity(TopicityAst { pair: LigandPairAst::new(StereoLigandId(0), StereoLigandId(1)), rel: TopicityRelationAst::NotSet(vec![Topicity::Enantiotopic]) }))]
    #[case::topicity_open("Th2#o*(0,1)",
        StereoAtomConstraint::Topicity(TopicityAst { pair: LigandPairAst::new(StereoLigandId(0), StereoLigandId(1)), rel: TopicityRelationAst::Undetermined }))]
    #[case::stereogenicity("Th2#g/",
        StereoAtomConstraint::Stereogenicity(StereogenicityAst(StereogenicityRelationAst::Lit(Stereogenicity::Stereogenic))))]
    fn test_stereo_atom_predicate(#[case] input: &str, #[case] expected: StereoAtomConstraint) {
        let dsl = parse_stereo_atom(input).unwrap();
        assert_eq!(dsl.0.constraints.iter().cloned().collect::<Vec<_>>(), vec![expected]);
    }

    #[rstest]
    fn test_stereo_atom_predicate_involution() {
        let dsl = parse_stereo_atom("Th2#p~").unwrap();
        let expected = StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst {
            perm: OrientedPermutationAst {
                perm: PermutationAst(StereoKind::Tetrahedral.involution()),
                orientation: Orientation::Improper,
            },
            mem: MemOp::In,
        });
        assert_eq!(
            dsl.0.constraints.iter().cloned().collect::<Vec<_>>(),
            vec![expected],
        );
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
    #[case::fluxionality_involution("Ct1#f~")]
    #[case::topicity("Ct1#o=(0,1)")]
    #[case::stereogenicity("Ct1#g/")]
    fn test_stereo_bond_inline_render_identity(#[case] s: &str) {
        assert_eq!(parse_stereo_bond(s).unwrap().to_string(), s);
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
        assert!(matches!(
            err,
            DeError::Subgrammar {
                grammar: "stereo coset",
                ..
            }
        ));
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
