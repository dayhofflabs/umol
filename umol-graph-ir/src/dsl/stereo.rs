//! Stereo config-string DSL.

use std::borrow::Cow;
use std::collections::BTreeSet;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnMap, EdnStreamDeserializer, FromEdn, ToEdn};
use umol_perm::{Orientation, Permutation};
use winnow::ascii::{digit1, multispace0};
use winnow::combinator::{
    alt, delimited, opt, preceded, repeat, separated, separated_pair, terminated,
};
use winnow::error::ErrMode;
use winnow::token::any;
use winnow::Parser;

use super::boolean::{boolean, fmt_boolean, BooleanDsl};
use super::config::{StereoAtomDefaults, StereoBondDefaults};
use super::edn_utils::single_key_map;
use super::error::{PResult, ParseError};
use super::value::variable_name;
use crate::ir::boolean::BooleanForm;
use crate::ir::constraint::{
    FluxionalityForm, LigandPermutation, LigandSymmetryForm, OrientedLigandPermutation,
    StereoAtomConstraintForm, StereoAtomConstraintsForm, StereoBondConstraintForm,
    StereoBondConstraintsForm, StereoLigandPair, StereogenicityForm, TopicityForm,
    TopicityRelationForm,
};
use crate::ir::id::StereoLigandPosition;
use crate::ir::stereo::{
    CisTransStereoForm, StereoAtomForm, StereoAtomUpdate, StereoBondForm, StereoBondUpdate,
    StereoConfigurationForm, StereoConfigurationUpdate, StereoCoset, StereoKind, StereoTerm,
    Stereogenicity, TetrahedralStereoForm, Topicity,
};
use crate::ir::traits::{FromIr, IntoIr, Lattice};

/// Surface DSL wrapper for `StereoAtomForm`
#[repr(transparent)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoAtomDsl(pub StereoAtomForm);

impl StereoAtomDsl {
    /// Zero-cost reference cast from `&StereoAtomForm`. Relies on `repr(transparent)`.
    pub fn from_ref(form: &StereoAtomForm) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(form as *const StereoAtomForm as *const Self) }
    }
}

impl From<StereoAtomForm> for StereoAtomDsl {
    fn from(form: StereoAtomForm) -> Self {
        Self(form)
    }
}

impl FromStr for StereoAtomDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_stereo_atom(s)
    }
}

impl FromStr for StereoAtomForm {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(StereoAtomDsl::from_str(s)?.into_ir(&StereoAtomDefaults::default()))
    }
}

impl Display for StereoAtomForm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        StereoAtomDsl::from_ref(self).fmt(f)
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
/// - `:ccw` -> `"Th0"`
/// - `:cw` -> `"Th1"`
///
/// Returns `None` for unrecorgnized keywords.
pub(crate) fn expand_stereo_atom_keyword(name: &str) -> Option<&'static str> {
    match name {
        "ccw" => Some("Th0"),
        "cw" => Some("Th1"),
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
fn stereo_atom_keyword_for(form: &StereoAtomForm) -> Option<&'static str> {
    if !form.constraints.is_empty() {
        return None;
    }
    match &form.configuration {
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(0)) => {
            Some("ccw")
        }
        StereoConfigurationForm::Kinded(StereoKind::Tetrahedral, StereoCoset::Lit(1)) => Some("cw"),
        _ => None,
    }
}

impl FromIr<StereoAtomForm> for StereoAtomDsl {
    type Ctx = StereoAtomDefaults;

    fn from_ir(form: &StereoAtomForm, _ctx: &Self::Ctx) -> Self {
        let form = form.clone();
        StereoAtomDsl(form.clone())
    }
}

impl IntoIr<StereoAtomForm> for StereoAtomDsl {
    type Ctx = StereoAtomDefaults;

    fn into_ir(self, _ctx: &Self::Ctx) -> StereoAtomForm {
        self.0
    }
}

pub fn parse_stereo_atom(input: &str) -> Result<StereoAtomDsl, ParseError> {
    stereo_atom.parse(input).map_err(|e| e.into_inner())
}

/// Parse the configuration head of a stereo element: `*` (undetermined geometry,
/// no coset) or `<kind><coset>` (a concrete kind with a coset). `*` and `Th*` are distinct.
fn stereo_configuration(i: &mut &str) -> PResult<StereoConfigurationForm> {
    multispace0.parse_next(i)?;
    if opt('*').parse_next(i)?.is_some() {
        multispace0.parse_next(i)?;
        return Ok(StereoConfigurationForm::Undetermined);
    }
    let kind = delimited(multispace0, stereo_kind, multispace0).parse_next(i)?;
    let coset = terminated(
        move |i: &mut &str| stereo_coset(i, kind.degree()),
        multispace0,
    )
    .parse_next(i)?;
    Ok(StereoConfigurationForm::kinded(kind, coset))
}

/// Partial-modify variant of `stereo_configuration`: the coset is optional (omitted =
/// `Undetermined`, "coset unchanged"). So `*` alone, `<kind>`, or `<kind><coset>`; the kind is
/// mandatory once anything past `*` appears.
pub(crate) fn stereo_atom(i: &mut &str) -> PResult<StereoAtomDsl> {
    let configuration = stereo_configuration(i)?;
    stereo_atom_tail(i, configuration)
}

/// Parse the kind-scoped constraint predicates and trailing-input check given an already-parsed
/// configuration.
fn stereo_atom_tail(
    i: &mut &str,
    configuration: StereoConfigurationForm,
) -> PResult<StereoAtomDsl> {
    let mut constraints = StereoAtomConstraintsForm::new();
    if let StereoConfigurationForm::Kinded(kind, _) = &configuration {
        let kind = *kind;
        loop {
            let before = *i;
            match stereo_atom_predicate(i, kind) {
                Ok(c) => {
                    if constraints.contains(c.key()) {
                        return Err(ErrMode::Cut(ParseError::DuplicateStereoPredicate(
                            before[..2].to_string(),
                        )));
                    }
                    constraints.set(c);
                }
                Err(ErrMode::Backtrack(_)) => {
                    *i = before;
                    break;
                }
                Err(e) => return Err(e),
            }
        }
    }
    multispace0.parse_next(i)?;
    if !i.is_empty() {
        return Err(ErrMode::Cut(ParseError::TrailingInput((*i).to_string())));
    }
    Ok(StereoAtomDsl(StereoAtomForm {
        configuration,
        constraints,
    }))
}

/// Surface DSL wrapper for `StereoBondForm`
#[repr(transparent)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoBondDsl(pub StereoBondForm);

impl StereoBondDsl {
    /// Zero-cost reference cast from `&StereoBondForm`. Relies on `repr(transparent)`.
    pub fn from_ref(form: &StereoBondForm) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(form as *const StereoBondForm as *const Self) }
    }
}

impl From<StereoBondForm> for StereoBondDsl {
    fn from(form: StereoBondForm) -> Self {
        Self(form)
    }
}

impl FromStr for StereoBondDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_stereo_bond(s)
    }
}

impl FromStr for StereoBondForm {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(StereoBondDsl::from_str(s)?.into_ir(&StereoBondDefaults::default()))
    }
}

impl Display for StereoBondForm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        StereoBondDsl::from_ref(self).fmt(f)
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
/// - `:z` -> `"Ct0"`
/// - `:e` -> `"Ct1"`
///
/// Returns `None` for unrecorgnized keywords.
pub(crate) fn expand_stereo_bond_keyword(name: &str) -> Option<&'static str> {
    match name {
        "z" => Some("Ct0"),
        "e" => Some("Ct1"),
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
fn stereo_bond_keyword_for(form: &StereoBondForm) -> Option<&'static str> {
    if !form.constraints.is_empty() {
        return None;
    }
    match &form.configuration {
        StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(0)) => Some("z"),
        StereoConfigurationForm::Kinded(StereoKind::CisTrans, StereoCoset::Lit(1)) => Some("e"),
        _ => None,
    }
}

impl FromIr<StereoBondForm> for StereoBondDsl {
    type Ctx = StereoBondDefaults;

    fn from_ir(form: &StereoBondForm, _ctx: &Self::Ctx) -> Self {
        StereoBondDsl(form.clone())
    }
}

impl IntoIr<StereoBondForm> for StereoBondDsl {
    type Ctx = StereoBondDefaults;

    fn into_ir(self, _ctx: &Self::Ctx) -> StereoBondForm {
        self.0
    }
}

pub fn parse_stereo_bond(input: &str) -> Result<StereoBondDsl, ParseError> {
    stereo_bond.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn stereo_bond(i: &mut &str) -> PResult<StereoBondDsl> {
    let configuration = stereo_configuration(i)?;
    stereo_bond_tail(i, configuration)
}

/// Parse the kind-scoped constraint predicates and trailing-input check given an already-parsed
/// configuration.
fn stereo_bond_tail(
    i: &mut &str,
    configuration: StereoConfigurationForm,
) -> PResult<StereoBondDsl> {
    let mut constraints = StereoBondConstraintsForm::new();
    if let StereoConfigurationForm::Kinded(kind, _) = &configuration {
        let kind = *kind;
        loop {
            let before = *i;
            match stereo_bond_predicate(i, kind) {
                Ok(c) => {
                    if constraints.contains(c.key()) {
                        return Err(ErrMode::Cut(ParseError::DuplicateStereoPredicate(
                            before[..2].to_string(),
                        )));
                    }
                    constraints.set(c);
                }
                Err(ErrMode::Backtrack(_)) => {
                    *i = before;
                    break;
                }
                Err(e) => return Err(e),
            }
        }
    }
    multispace0.parse_next(i)?;
    if !i.is_empty() {
        return Err(ErrMode::Cut(ParseError::TrailingInput((*i).to_string())));
    }
    Ok(StereoBondDsl(StereoBondForm {
        configuration,
        constraints,
    }))
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

pub(crate) fn parse_stereo_coset(input: &str, degree: usize) -> Result<StereoCoset, ParseError> {
    (move |i: &mut &str| stereo_coset(i, degree))
        .parse(input)
        .map_err(|e| e.into_inner())
}

/// Parse a permutation in compact cycle notation (`(0,1)(2,3)`) over a `degree`-position
/// ligand frame. Points are validated in-range and disjoint.
pub(crate) fn parse_permutation(input: &str, degree: usize) -> Result<Permutation, ParseError> {
    (move |i: &mut &str| perm_cycles(i, degree))
        .parse(input)
        .map_err(|e| e.into_inner())
}

/// Parse the `coset` grammar into `StereoCoset` over a `degree`-position
/// ligand frame (a coset-expression's `^` permutation acts on those positions).
fn stereo_coset(i: &mut &str, degree: usize) -> PResult<StereoCoset> {
    alt((
        '*'.value(StereoCoset::Undetermined),
        (|i: &mut &str| stereo_term(i, degree)).map(|t| match t {
            StereoTerm::Lit(n) => StereoCoset::Lit(n),
            StereoTerm::LitSet(s) => StereoCoset::LitSet(s),
            other => StereoCoset::term(other),
        }),
    ))
    .parse_next(i)
}

/// `stereo-expr`: a `~`-prefixed base carrying zero or more left-associative
/// `^cycles` postfixes (the permutation in 0-indexed disjoint-cycle notation).
fn stereo_term(i: &mut &str, degree: usize) -> PResult<StereoTerm> {
    let mut e = stereo_prefix(i)?;
    loop {
        multispace0.parse_next(i)?;
        if opt('^').parse_next(i)?.is_some() {
            e = StereoTerm::apply(e, perm_cycles(i, degree)?);
        } else {
            return Ok(e);
        }
    }
}

/// `('~' | '\'') prefix-term | base` — unary `~` (swap) and `'` (mirror) bind
/// tighter than `^`.
fn stereo_prefix(i: &mut &str) -> PResult<StereoTerm> {
    multispace0.parse_next(i)?;
    if opt('~').parse_next(i)?.is_some() {
        Ok(StereoTerm::swap(stereo_prefix(i)?))
    } else if opt('\'').parse_next(i)?.is_some() {
        Ok(StereoTerm::mirror(stereo_prefix(i)?))
    } else {
        stereo_base(i)
    }
}

/// `nat | '?' variable-name ('::' set)? | set`.
fn stereo_base(i: &mut &str) -> PResult<StereoTerm> {
    preceded(
        multispace0,
        alt((
            stereo_lit_set.map(StereoTerm::lit_set),
            stereo_var,
            stereo_lit.map(StereoTerm::Lit),
        )),
    )
    .parse_next(i)
}

fn stereo_var(i: &mut &str) -> PResult<StereoTerm> {
    let name = preceded('?', variable_name).parse_next(i)?;
    let domain = opt(preceded((multispace0, "::", multispace0), stereo_lit_set)).parse_next(i)?;
    Ok(match domain {
        Some(set) => StereoTerm::var_in(name, set),
        None => StereoTerm::var(name),
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
/// (0-indexed).
fn perm_cycles(i: &mut &str, degree: usize) -> PResult<Permutation> {
    let cycles: Vec<Vec<usize>> = repeat(1.., cycle).parse_next(i)?;
    Permutation::from_cycles(degree, &cycles)
        .map_err(|error| ErrMode::Cut(ParseError::InvalidValue(error.to_string())))
}

/// `~` (the kind involution, eager) or an explicit permutation in cycle notation.
fn stereo_permutation(i: &mut &str, kind: StereoKind) -> PResult<Permutation> {
    if opt('~').parse_next(i)?.is_some() {
        Ok(kind.involution())
    } else {
        perm_cycles(i, kind.degree())
    }
}

fn ligand_pair(i: &mut &str) -> PResult<StereoLigandPair> {
    let (a, b) =
        delimited('(', separated_pair(cycle_point, ',', cycle_point), ')').parse_next(i)?;
    Ok(StereoLigandPair::new(
        StereoLigandPosition(a as u32),
        StereoLigandPosition(b as u32),
    ))
}

fn topicity_char(i: &mut &str) -> PResult<Topicity> {
    alt((
        '='.value(Topicity::Homotopic),
        '\''.value(Topicity::Enantiotopic),
        '/'.value(Topicity::Diastereotopic),
    ))
    .parse_next(i)
}

fn stereogenicity_char(i: &mut &str) -> PResult<Stereogenicity> {
    alt((
        '='.value(Stereogenicity::Symmetric),
        '\''.value(Stereogenicity::Prochiral),
        '/'.value(Stereogenicity::Stereogenic),
    ))
    .parse_next(i)
}

/// Parse a `{ glyph (',' glyph)* }` glyph set — the `LitSet` / multi-`NotSet` body of a relation.
fn char_set<T: Ord>(
    i: &mut &str,
    mut ch: impl FnMut(&mut &str) -> PResult<T>,
) -> PResult<BTreeSet<T>> {
    '{'.parse_next(i)?;
    let mut set = BTreeSet::from([ch(i)?]);
    while opt(',').parse_next(i)?.is_some() {
        set.insert(ch(i)?);
    }
    '}'.parse_next(i)?;
    Ok(set)
}

fn topicity_relation_inline(i: &mut &str) -> PResult<TopicityRelationForm> {
    if opt('*').parse_next(i)?.is_some() {
        return Ok(TopicityRelationForm::Undetermined);
    }
    let neg = opt('!').parse_next(i)?.is_some();
    if i.starts_with('{') {
        let set = char_set(i, topicity_char)?;
        Ok(if neg {
            TopicityRelationForm::NotSet(set)
        } else {
            TopicityRelationForm::LitSet(set)
        })
    } else {
        let v = topicity_char(i)?;
        Ok(if neg {
            TopicityRelationForm::NotSet(BTreeSet::from([v]))
        } else {
            TopicityRelationForm::Lit(v)
        })
    }
}

fn stereogenicity_relation_inline(i: &mut &str) -> PResult<StereogenicityForm> {
    if opt('*').parse_next(i)?.is_some() {
        return Ok(StereogenicityForm::Undetermined);
    }
    let neg = opt('!').parse_next(i)?.is_some();
    if i.starts_with('{') {
        let set = char_set(i, stereogenicity_char)?;
        Ok(if neg {
            StereogenicityForm::NotSet(set)
        } else {
            StereogenicityForm::LitSet(set)
        })
    } else {
        let v = stereogenicity_char(i)?;
        Ok(if neg {
            StereogenicityForm::NotSet(BTreeSet::from([v]))
        } else {
            StereogenicityForm::Lit(v)
        })
    }
}

/// Parse one `#p`/`#f`/`#o`/`#g` predicate for `kind`. Atom and bond share the
/// grammar; the macro emits a parser per constraint type.
macro_rules! stereo_predicate_parser {
    ($name:ident, $constraint:ident) => {
        fn $name(i: &mut &str, kind: StereoKind) -> PResult<$constraint> {
            '#'.parse_next(i)?;
            let tag = any.parse_next(i)?;
            match tag {
                'p' => {
                    let (permutation, orientation) = if opt('~').parse_next(i)?.is_some() {
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
                    let invariant = boolean(i)?.0;
                    Ok($constraint::LigandSymmetry(LigandSymmetryForm {
                        permutation: OrientedLigandPermutation {
                            permutation: LigandPermutation(permutation),
                            orientation,
                        },
                        invariant,
                    }))
                }
                'f' => {
                    let permutation = LigandPermutation(stereo_permutation(i, kind)?);
                    let active = boolean(i)?.0;
                    Ok($constraint::Fluxionality(FluxionalityForm {
                        permutation,
                        active,
                    }))
                }
                'o' => {
                    let pair = ligand_pair(i)?;
                    let relation = topicity_relation_inline(i)?;
                    Ok($constraint::Topicity(TopicityForm { pair, relation }))
                }
                'g' => Ok($constraint::Stereogenicity(stereogenicity_relation_inline(
                    i,
                )?)),
                other => Err(ErrMode::Cut(ParseError::UnknownStereoPredicate(format!(
                    "#{other}"
                )))),
            }
        }
    };
}

stereo_predicate_parser! { stereo_atom_predicate, StereoAtomConstraintForm }
stereo_predicate_parser! { stereo_bond_predicate, StereoBondConstraintForm }

/// `~` (when the permutation is the kind involution) or an explicit permutation in
/// cycle notation (`Permutation`'s `Display`).
fn fmt_stereo_permutation(
    f: &mut fmt::Formatter<'_>,
    permutation: Permutation,
    kind: StereoKind,
) -> fmt::Result {
    if permutation == kind.involution() {
        f.write_str("~")
    } else {
        write!(f, "{permutation}")
    }
}

fn write_topicity_char(t: Topicity) -> char {
    match t {
        Topicity::Homotopic => '=',
        Topicity::Enantiotopic => '\'',
        Topicity::Diastereotopic => '/',
    }
}

fn write_stereogenicity_char(s: Stereogenicity) -> char {
    match s {
        Stereogenicity::Symmetric => '=',
        Stereogenicity::Prochiral => '\'',
        Stereogenicity::Stereogenic => '/',
    }
}

/// Write a `{ glyph,glyph,… }` glyph set (a `LitSet` body, or a multi-element `NotSet` after `!`).
fn fmt_char_set<T: Copy>(
    f: &mut fmt::Formatter<'_>,
    set: &BTreeSet<T>,
    ch: impl Fn(T) -> char,
) -> fmt::Result {
    write!(f, "{{")?;
    for (idx, v) in set.iter().enumerate() {
        if idx > 0 {
            write!(f, ",")?;
        }
        write!(f, "{}", ch(*v))?;
    }
    write!(f, "}}")
}

fn fmt_topicity_relation(f: &mut fmt::Formatter<'_>, rel: &TopicityRelationForm) -> fmt::Result {
    match rel {
        TopicityRelationForm::Undetermined => f.write_str("*"),
        TopicityRelationForm::Lit(t) => write!(f, "{}", write_topicity_char(*t)),
        TopicityRelationForm::LitSet(s) => fmt_char_set(f, s, write_topicity_char),
        TopicityRelationForm::NotSet(s) => {
            write!(f, "!")?;
            match s.iter().next() {
                Some(&t) if s.len() == 1 => write!(f, "{}", write_topicity_char(t)),
                _ => fmt_char_set(f, s, write_topicity_char),
            }
        }
    }
}

fn fmt_stereogenicity_relation(
    f: &mut fmt::Formatter<'_>,
    rel: &StereogenicityForm,
) -> fmt::Result {
    match rel {
        StereogenicityForm::Undetermined => f.write_str("*"),
        StereogenicityForm::Lit(s) => write!(f, "{}", write_stereogenicity_char(*s)),
        StereogenicityForm::LitSet(set) => fmt_char_set(f, set, write_stereogenicity_char),
        StereogenicityForm::NotSet(set) => {
            write!(f, "!")?;
            match set.iter().next() {
                Some(&s) if set.len() == 1 => write!(f, "{}", write_stereogenicity_char(s)),
                _ => fmt_char_set(f, set, write_stereogenicity_char),
            }
        }
    }
}

/// Render one stereo constraint as its inline predicate. `~` is emitted for a
/// `#p`/`#f` permutation equal to the kind involution (matching the involution's
/// orientation for `#p`); otherwise explicit cycles.
macro_rules! stereo_constraint_fmt {
    ($name:ident, $constraint:ident) => {
        fn $name(f: &mut fmt::Formatter<'_>, c: &$constraint, kind: StereoKind) -> fmt::Result {
            match c {
                $constraint::LigandSymmetry(ls) => {
                    f.write_str("#p")?;
                    let involution_orientation = if kind.is_chiral_class() {
                        Orientation::Improper
                    } else {
                        Orientation::Proper
                    };
                    if ls.permutation.permutation.0 == kind.involution()
                        && ls.permutation.orientation == involution_orientation
                    {
                        f.write_str("~")?;
                    } else {
                        if ls.permutation.orientation == Orientation::Improper {
                            f.write_str("'")?;
                        }
                        write!(f, "{}", ls.permutation.permutation.0)?;
                    }
                    fmt_boolean(f, &ls.invariant)
                }
                $constraint::Fluxionality(fx) => {
                    f.write_str("#f")?;
                    fmt_stereo_permutation(f, fx.permutation.0, kind)?;
                    fmt_boolean(f, &fx.active)
                }
                $constraint::Topicity(t) => {
                    f.write_str("#o")?;
                    write!(f, "({},{})", t.pair.first().0, t.pair.second().0)?;
                    fmt_topicity_relation(f, &t.relation)
                }
                $constraint::Stereogenicity(g) => {
                    f.write_str("#g")?;
                    fmt_stereogenicity_relation(f, &g)
                }
            }
        }
    };
}

stereo_constraint_fmt! { fmt_stereo_atom_constraint, StereoAtomConstraintForm }
stereo_constraint_fmt! { fmt_stereo_bond_constraint, StereoBondConstraintForm }

/// Write the stereo atom DSL
/// Render the configuration head of a stereo element: `*` or `<kind><coset>`.
fn fmt_stereo_configuration(
    f: &mut fmt::Formatter<'_>,
    configuration: &StereoConfigurationForm,
) -> fmt::Result {
    match configuration {
        StereoConfigurationForm::Undetermined => write!(f, "*"),
        StereoConfigurationForm::Kinded(kind, coset) => {
            fmt_stereo_kind(f, *kind)?;
            fmt_stereo_coset(f, coset)
        }
    }
}

pub(crate) fn fmt_stereo_atom(f: &mut fmt::Formatter<'_>, atom: &StereoAtomForm) -> fmt::Result {
    fmt_stereo_configuration(f, &atom.configuration)?;
    // Vacuous (Undetermined) predicates are elided on canonical render, as for
    // the atom-string; they remain admissible on parse.
    if let StereoConfigurationForm::Kinded(kind, _) = &atom.configuration {
        for c in atom.constraints.iter().filter(|c| !c.is_undetermined()) {
            fmt_stereo_atom_constraint(f, c, *kind)?;
        }
    }
    Ok(())
}

/// Write the stereo bond DSL
pub(crate) fn fmt_stereo_bond(f: &mut fmt::Formatter<'_>, bond: &StereoBondForm) -> fmt::Result {
    fmt_stereo_configuration(f, &bond.configuration)?;
    if let StereoConfigurationForm::Kinded(kind, _) = &bond.configuration {
        for c in bond.constraints.iter().filter(|c| !c.is_undetermined()) {
            fmt_stereo_bond_constraint(f, c, *kind)?;
        }
    }
    Ok(())
}

pub fn parse_stereo_atom_update(input: &str) -> Result<StereoAtomUpdateDsl, ParseError> {
    stereo_atom_update.parse(input).map_err(|e| e.into_inner())
}

pub fn parse_stereo_bond_update(input: &str) -> Result<StereoBondUpdateDsl, ParseError> {
    stereo_bond_update.parse(input).map_err(|e| e.into_inner())
}

/// Surface DSL wrapper around a [`StereoAtomUpdate`].
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StereoAtomUpdateDsl(pub StereoAtomUpdate);

impl StereoAtomUpdateDsl {
    /// Zero-cost reference cast from `&StereoAtomUpdate`. Relies on `repr(transparent)`.
    pub fn from_ref(update: &StereoAtomUpdate) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(update as *const StereoAtomUpdate as *const Self) }
    }
}

impl FromIr<StereoAtomUpdate> for StereoAtomUpdateDsl {
    type Ctx = ();

    fn from_ir(update: &StereoAtomUpdate, _ctx: &Self::Ctx) -> Self {
        Self(update.clone())
    }
}

impl IntoIr<StereoAtomUpdate> for StereoAtomUpdateDsl {
    type Ctx = ();

    fn into_ir(self, _ctx: &Self::Ctx) -> StereoAtomUpdate {
        self.0
    }
}

impl FromStr for StereoAtomUpdateDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_stereo_atom_update(s)
    }
}

impl FromStr for StereoAtomUpdate {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(StereoAtomUpdateDsl::from_str(s)?.into_ir(&()))
    }
}

impl Display for StereoAtomUpdate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        StereoAtomUpdateDsl::from_ref(self).fmt(f)
    }
}

impl Display for StereoAtomUpdateDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let kind = match &self.0.configuration {
            StereoConfigurationUpdate::Unchanged => None,
            StereoConfigurationUpdate::Undetermined => {
                f.write_str("*")?;
                None
            }
            StereoConfigurationUpdate::Kinded { kind, coset } => {
                fmt_stereo_kind(f, *kind)?;
                if let Some(coset) = coset {
                    fmt_stereo_coset(f, coset)?;
                }
                Some(*kind)
            }
        };
        if let Some(kind) = kind {
            for constraint in self.0.constraints.iter() {
                fmt_stereo_atom_constraint(f, constraint, kind)?;
            }
        } else if !self.0.constraints.is_empty() {
            return Err(fmt::Error);
        }
        Ok(())
    }
}

impl<'de> FromEdn<'de> for StereoAtomUpdateDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s
                .parse()
                .map_err(|e| DeError::subgrammar("stereo-atom-update", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for StereoAtomUpdateDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

fn stereo_atom_update(i: &mut &str) -> PResult<StereoAtomUpdateDsl> {
    multispace0.parse_next(i)?;
    if i.is_empty() {
        return Ok(StereoAtomUpdateDsl::default());
    }
    let configuration = if opt('*').parse_next(i)?.is_some() {
        StereoConfigurationUpdate::Undetermined
    } else {
        let kind = stereo_kind(i)?;
        let coset = opt(move |i: &mut &str| stereo_coset(i, kind.degree())).parse_next(i)?;
        StereoConfigurationUpdate::Kinded { kind, coset }
    };
    multispace0.parse_next(i)?;
    let mut constraints = StereoAtomConstraintsForm::new();
    if let StereoConfigurationUpdate::Kinded { kind, .. } = &configuration {
        let kind = *kind;
        loop {
            let before = *i;
            match stereo_atom_predicate(i, kind) {
                Ok(constraint) => {
                    if constraints.contains(constraint.key()) {
                        return Err(ErrMode::Cut(ParseError::DuplicateStereoPredicate(
                            before[..2].to_string(),
                        )));
                    }
                    constraints.set(constraint);
                }
                Err(ErrMode::Backtrack(_)) => {
                    *i = before;
                    break;
                }
                Err(error) => return Err(error),
            }
        }
    }
    multispace0.parse_next(i)?;
    if !i.is_empty() {
        return Err(ErrMode::Cut(ParseError::TrailingInput((*i).to_string())));
    }
    Ok(StereoAtomUpdateDsl(StereoAtomUpdate {
        configuration,
        constraints,
    }))
}

/// Surface DSL wrapper around a [`StereoBondUpdate`].
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct StereoBondUpdateDsl(pub StereoBondUpdate);

impl StereoBondUpdateDsl {
    /// Zero-cost reference cast from `&StereoBondUpdate`. Relies on `repr(transparent)`.
    pub fn from_ref(update: &StereoBondUpdate) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(update as *const StereoBondUpdate as *const Self) }
    }
}

impl FromIr<StereoBondUpdate> for StereoBondUpdateDsl {
    type Ctx = ();

    fn from_ir(update: &StereoBondUpdate, _ctx: &Self::Ctx) -> Self {
        Self(update.clone())
    }
}

impl IntoIr<StereoBondUpdate> for StereoBondUpdateDsl {
    type Ctx = ();

    fn into_ir(self, _ctx: &Self::Ctx) -> StereoBondUpdate {
        self.0
    }
}

impl FromStr for StereoBondUpdateDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_stereo_bond_update(s)
    }
}

impl FromStr for StereoBondUpdate {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(StereoBondUpdateDsl::from_str(s)?.into_ir(&()))
    }
}

impl Display for StereoBondUpdate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        StereoBondUpdateDsl::from_ref(self).fmt(f)
    }
}

impl Display for StereoBondUpdateDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let kind = match &self.0.configuration {
            StereoConfigurationUpdate::Unchanged => None,
            StereoConfigurationUpdate::Undetermined => {
                f.write_str("*")?;
                None
            }
            StereoConfigurationUpdate::Kinded { kind, coset } => {
                fmt_stereo_kind(f, *kind)?;
                if let Some(coset) = coset {
                    fmt_stereo_coset(f, coset)?;
                }
                Some(*kind)
            }
        };
        if let Some(kind) = kind {
            for constraint in self.0.constraints.iter() {
                fmt_stereo_bond_constraint(f, constraint, kind)?;
            }
        } else if !self.0.constraints.is_empty() {
            return Err(fmt::Error);
        }
        Ok(())
    }
}

impl<'de> FromEdn<'de> for StereoBondUpdateDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s
                .parse()
                .map_err(|e| DeError::subgrammar("stereo-bond-update", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for StereoBondUpdateDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

fn stereo_bond_update(i: &mut &str) -> PResult<StereoBondUpdateDsl> {
    multispace0.parse_next(i)?;
    if i.is_empty() {
        return Ok(StereoBondUpdateDsl::default());
    }
    let configuration = if opt('*').parse_next(i)?.is_some() {
        StereoConfigurationUpdate::Undetermined
    } else {
        let kind = stereo_kind(i)?;
        let coset = opt(move |i: &mut &str| stereo_coset(i, kind.degree())).parse_next(i)?;
        StereoConfigurationUpdate::Kinded { kind, coset }
    };
    multispace0.parse_next(i)?;
    let mut constraints = StereoBondConstraintsForm::new();
    if let StereoConfigurationUpdate::Kinded { kind, .. } = &configuration {
        let kind = *kind;
        loop {
            let before = *i;
            match stereo_bond_predicate(i, kind) {
                Ok(constraint) => {
                    if constraints.contains(constraint.key()) {
                        return Err(ErrMode::Cut(ParseError::DuplicateStereoPredicate(
                            before[..2].to_string(),
                        )));
                    }
                    constraints.set(constraint);
                }
                Err(ErrMode::Backtrack(_)) => {
                    *i = before;
                    break;
                }
                Err(error) => return Err(error),
            }
        }
    }
    multispace0.parse_next(i)?;
    if !i.is_empty() {
        return Err(ErrMode::Cut(ParseError::TrailingInput((*i).to_string())));
    }
    Ok(StereoBondUpdateDsl(StereoBondUpdate {
        configuration,
        constraints,
    }))
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

/// Write a `StereoCoset` for the element `:type` body: `*` (open coset), a
/// literal, or an operator expression. `fmt_stereo_config` reuses the literal
/// and expression arms but writes its own `+` for `Stereo(Undetermined)`.
fn fmt_stereo_coset(f: &mut fmt::Formatter<'_>, coset: &StereoCoset) -> fmt::Result {
    match coset {
        StereoCoset::Undetermined => write!(f, "*"),
        StereoCoset::Lit(n) => write!(f, "{n}"),
        StereoCoset::LitSet(s) => fmt_stereo_lit_set(f, s),
        StereoCoset::Term(t) => fmt_stereo_term(f, t),
    }
}

fn fmt_stereo_term(f: &mut fmt::Formatter<'_>, t: &StereoTerm) -> fmt::Result {
    match t {
        StereoTerm::Var(v) => match &v.1 {
            None => write!(f, "?{}", v.0),
            Some(set) => {
                write!(f, "?{} :: ", v.0)?;
                fmt_stereo_lit_set(f, set)
            }
        },
        StereoTerm::Lit(n) => write!(f, "{n}"),
        StereoTerm::LitSet(set) => fmt_stereo_lit_set(f, set),
        StereoTerm::Swap(inner) => {
            write!(f, "~")?;
            fmt_stereo_term(f, inner)
        }
        StereoTerm::Mirror(inner) => {
            write!(f, "'")?;
            fmt_stereo_term(f, inner)
        }
        StereoTerm::Apply(inner, permutation) => {
            fmt_stereo_term(f, inner)?;
            write!(f, "^{permutation}")
        }
    }
}

fn fmt_stereo_lit_set(f: &mut fmt::Formatter<'_>, set: &BTreeSet<u32>) -> fmt::Result {
    write!(f, "{{")?;
    for (i, n) in set.iter().enumerate() {
        if i > 0 {
            write!(f, ",")?;
        }
        write!(f, "{n}")?;
    }
    write!(f, "}}")
}

/// A permutation as a vector of disjoint cycles (`[[0 1 2] [3 4]]`; identity `[]`).
/// The degree is not encoded — fixed points drop out — so the reader supplies it
/// from the stereo kind.
fn render_edn_permutation(permutation: Permutation) -> Edn<'static> {
    let cycles: Vec<Edn<'static>> = permutation
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

fn read_edn_permutation(edn: &Edn, degree: usize) -> Result<Permutation, DeError> {
    let Edn::Vector(cycles_edn) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "vector of cycles",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
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
            let point = usize::try_from(*n).map_err(|_| DeError::OutOfRange {
                value: n.to_string(),
                target: "ligand position",
                path: Vec::new(),
            })?;
            cycle.push(point);
        }
        cycles.push(cycle);
    }
    Permutation::from_cycles(degree, &cycles).map_err(|error| DeError::Custom(error.to_string()))
}

/// Streaming-read intermediate for a `:relation` value (`:undetermined`, a single
/// keyword, a keyword vector, or a `{:not-in [members]}` complement map) before the
/// domain mapping is applied by `relation_serde!`'s `$from_parts`.
#[derive(Clone, Debug)]
pub(crate) enum RelationValue {
    Undetermined,
    One(String),
    Many(Vec<String>),
    NotIn(Vec<String>),
}

/// Generates the structured-EDN serialization/deserialization for a relation type
/// into/out of the owning constraint map: `:relation` is `:undetermined`, a single
/// keyword (`Lit`), a member vector (`LitSet`), or a nested `{:not-in [members]}`
/// complement map (`NotSet`). Mechanical (no folding); keyword names map 1:1 to the
/// domain variants per the table. `$from_parts` is the streaming twin of `$from`,
/// mapping a `RelationValue` to the relation.
macro_rules! relation_serde {
    ($to:ident, $from:ident, $from_parts:ident, $relation:ident, $domain:ty, $($variant:path => $kw:literal),+ $(,)?) => {
        pub(crate) fn $from_parts(value: RelationValue) -> Result<$relation, DeError> {
            fn keyword_member(name: &str) -> Option<$domain> {
                match name {
                    $($kw => Some($variant),)+
                    _ => None,
                }
            }
            let member = |k: &str| {
                keyword_member(k).ok_or_else(|| {
                    DeError::Custom(format!(
                        concat!("unknown ", stringify!($relation), " keyword :{}"),
                        k
                    ))
                })
            };
            let member_set = |ks: &[String]| -> Result<BTreeSet<$domain>, DeError> {
                let mut set = BTreeSet::new();
                for k in ks {
                    set.insert(member(k)?);
                }
                Ok(set)
            };
            match value {
                RelationValue::Undetermined => Ok($relation::Undetermined),
                RelationValue::One(k) => Ok($relation::lit(member(&k)?)),
                RelationValue::Many(ks) => Ok($relation::lit_set(member_set(&ks)?)),
                RelationValue::NotIn(ks) => Ok($relation::not_set(member_set(&ks)?)),
            }
        }

        fn $to(rel: &$relation, m: &mut EdnMap<'static>) {
            fn member_kw(v: $domain) -> Edn<'static> {
                let name = match v { $($variant => $kw,)+ };
                Edn::Keyword(EdnKeyword::owned(name.to_string()))
            }
            fn member_vec(s: &BTreeSet<$domain>) -> Edn<'static> {
                Edn::Vector(s.iter().map(|v| member_kw(*v)).collect::<Vec<_>>().into())
            }
            let value = match rel {
                $relation::Undetermined => Edn::keyword("undetermined"),
                $relation::Lit(v) => member_kw(*v),
                $relation::LitSet(s) => member_vec(s),
                $relation::NotSet(s) => single_key_map("not-in", member_vec(s)),
            };
            m.insert(Edn::keyword("relation"), value);
        }

        fn $from(m: &EdnMap, path: &'static str) -> Result<$relation, DeError> {
            fn keyword_member(name: &str) -> Option<$domain> {
                match name {
                    $($kw => Some($variant),)+
                    _ => None,
                }
            }
            let member_set = |xs: &[Edn]| -> Result<BTreeSet<$domain>, DeError> {
                let mut set = BTreeSet::new();
                for e in xs {
                    let Edn::Keyword(k) = e else {
                        return Err(DeError::TypeMismatch {
                            expected: "relation keyword",
                            got: e.kind(),
                            path: vec![path.into()],
                        });
                    };
                    set.insert(keyword_member(k.name()).ok_or_else(|| DeError::TypeMismatch {
                        expected: "relation keyword",
                        got: e.kind(),
                        path: vec![path.into()],
                    })?);
                }
                Ok(set)
            };
            let value = m.get_keyword("relation").ok_or_else(|| DeError::MissingField {
                key: "relation".into(),
                path: vec![path.into()],
            })?;
            match value {
                Edn::Keyword(k) if k.name() == "undetermined" => Ok($relation::Undetermined),
                Edn::Keyword(k) => {
                    let v = keyword_member(k.name()).ok_or_else(|| DeError::TypeMismatch {
                        expected: concat!(stringify!($relation), " keyword"),
                        got: value.kind(),
                        path: vec![path.into()],
                    })?;
                    Ok($relation::lit(v))
                }
                Edn::Vector(xs) => Ok($relation::lit_set(member_set(xs)?)),
                Edn::Map(inner) => {
                    let complement =
                        inner.get_keyword("not-in").ok_or_else(|| DeError::MissingField {
                            key: "not-in".into(),
                            path: vec![path.into()],
                        })?;
                    let Edn::Vector(xs) = complement else {
                        return Err(DeError::TypeMismatch {
                            expected: "[members]",
                            got: complement.kind(),
                            path: vec![path.into()],
                        });
                    };
                    Ok($relation::not_set(member_set(xs)?))
                }
                other => Err(DeError::TypeMismatch {
                    expected: ":undetermined | keyword | [members] | {:not-in [members]}",
                    got: other.kind(),
                    path: vec![path.into()],
                }),
            }
        }
    };
}

relation_serde! {
    render_edn_topicity_relation, read_edn_topicity_relation, read_topicity_relation,
    TopicityRelationForm, Topicity,
    Topicity::Homotopic => "homotopic",
    Topicity::Enantiotopic => "enantiotopic",
    Topicity::Diastereotopic => "diastereotopic",
}

relation_serde! {
    render_edn_stereogenicity_relation, read_edn_stereogenicity_relation,
    read_stereogenicity_relation, StereogenicityForm, Stereogenicity,
    Stereogenicity::Symmetric => "symmetric",
    Stereogenicity::Prochiral => "prochiral",
    Stereogenicity::Stereogenic => "stereogenic",
}

/// EDN boundary for the `#g` stereogenicity constraint value. `StereogenicityForm`
/// is itself the relation, so this is its (de)serialization at the constraints map.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct StereogenicityDsl(pub StereogenicityForm);

impl ToEdn for StereogenicityDsl {
    fn to_edn(&self) -> Edn<'static> {
        let mut m = EdnMap::with_capacity(2);
        render_edn_stereogenicity_relation(&self.0, &mut m);
        Edn::Map(m)
    }
}

impl<'de> FromEdn<'de> for StereogenicityDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let Edn::Map(m) = edn else {
            return Err(DeError::TypeMismatch {
                expected: "stereogenicity map",
                got: edn.kind(),
                path: Vec::new(),
            });
        };
        Ok(Self(read_edn_stereogenicity_relation(m, "stereogenicity")?))
    }
}

impl FromIr<StereogenicityForm> for StereogenicityDsl {
    type Ctx = ();

    fn from_ir(form: &StereogenicityForm, _ctx: &Self::Ctx) -> Self {
        Self(form.clone())
    }
}

impl IntoIr<StereogenicityForm> for StereogenicityDsl {
    type Ctx = ();

    fn into_ir(self, _ctx: &Self::Ctx) -> StereogenicityForm {
        self.0
    }
}

/// EDN boundary for the `#o` topicity constraint value: `{:pair [i j] :relation
/// <rel>}`. The `TopicityRelationForm` rides inside (no separate relation DSL).
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct TopicityDsl(pub TopicityForm);

impl ToEdn for TopicityDsl {
    fn to_edn(&self) -> Edn<'static> {
        let t = &self.0;
        let mut m = EdnMap::with_capacity(3);
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
        render_edn_topicity_relation(&t.relation, &mut m);
        Edn::Map(m)
    }
}

impl<'de> FromEdn<'de> for TopicityDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
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
        let position = |e: &Edn| -> Result<StereoLigandPosition, DeError> {
            let Edn::Int(n) = e else {
                return Err(DeError::TypeMismatch {
                    expected: "int (ligand position)",
                    got: e.kind(),
                    path: vec!["topicity".into()],
                });
            };
            let v = u32::try_from(*n).map_err(|_| DeError::OutOfRange {
                value: n.to_string(),
                target: "ligand position",
                path: Vec::new(),
            })?;
            Ok(StereoLigandPosition(v))
        };
        let pair = StereoLigandPair::new(position(&p[0])?, position(&p[1])?);
        Ok(Self(TopicityForm {
            pair,
            relation: read_edn_topicity_relation(m, "topicity")?,
        }))
    }
}

impl FromIr<TopicityForm> for TopicityDsl {
    type Ctx = ();

    fn from_ir(form: &TopicityForm, _ctx: &Self::Ctx) -> Self {
        Self(form.clone())
    }
}

impl IntoIr<TopicityForm> for TopicityDsl {
    type Ctx = ();

    fn into_ir(self, _ctx: &Self::Ctx) -> TopicityForm {
        self.0
    }
}

/// `StereoKind` ↔ kebab keyword (`:tetrahedral`, `:cis-trans`, …).
pub(crate) fn render_edn_stereo_kind(kind: StereoKind) -> Edn<'static> {
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

pub(crate) fn stereo_kind_from_name(name: &str) -> Result<StereoKind, DeError> {
    match name {
        "tetrahedral" => Ok(StereoKind::Tetrahedral),
        "cis-trans" => Ok(StereoKind::CisTrans),
        "axial" => Ok(StereoKind::Axial),
        "square-planar" => Ok(StereoKind::SquarePlanar),
        "trigonal-bipyramidal" => Ok(StereoKind::TrigonalBipyramidal),
        "octahedral" => Ok(StereoKind::Octahedral),
        other => Err(DeError::Custom(format!("unknown stereo kind :{other}"))),
    }
}

pub(crate) fn read_edn_stereo_kind(edn: &Edn) -> Result<StereoKind, DeError> {
    let Edn::Keyword(k) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "stereo-kind keyword",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    stereo_kind_from_name(k.name())
}

fn render_edn_ligand_symmetry(ls: &LigandSymmetryForm) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(3);
    m.insert(
        Edn::keyword("permutation"),
        render_edn_permutation(ls.permutation.permutation.0),
    );
    if ls.permutation.orientation == Orientation::Improper {
        m.insert(Edn::keyword("orientation"), Edn::keyword("improper"));
    }
    if ls.invariant != BooleanForm::Lit(true) {
        m.insert(Edn::keyword("invariant"), BooleanDsl(ls.invariant).to_edn());
    }
    Edn::Map(m)
}

fn read_edn_ligand_symmetry(edn: &Edn, kind: StereoKind) -> Result<LigandSymmetryForm, DeError> {
    let Edn::Map(m) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "ligand-symmetry map",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    let permutation_edn = m
        .get_keyword("permutation")
        .ok_or_else(|| DeError::MissingField {
            key: "permutation".into(),
            path: vec!["ligand-symmetry".into()],
        })?;
    let permutation = LigandPermutation(read_edn_permutation(permutation_edn, kind.degree())?);
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
    let invariant = match m.get_keyword("invariant") {
        None => BooleanForm::Lit(true),
        Some(edn) => BooleanDsl::from_edn(edn)?.0,
    };
    Ok(LigandSymmetryForm {
        permutation: OrientedLigandPermutation {
            permutation,
            orientation,
        },
        invariant,
    })
}

fn render_edn_fluxionality(f: &FluxionalityForm) -> Edn<'static> {
    let mut m = EdnMap::with_capacity(2);
    m.insert(
        Edn::keyword("permutation"),
        render_edn_permutation(f.permutation.0),
    );
    if f.active != BooleanForm::Lit(true) {
        m.insert(Edn::keyword("active"), BooleanDsl(f.active).to_edn());
    }
    Edn::Map(m)
}

fn read_edn_fluxionality(edn: &Edn, kind: StereoKind) -> Result<FluxionalityForm, DeError> {
    let Edn::Map(m) = edn else {
        return Err(DeError::TypeMismatch {
            expected: "fluxionality map",
            got: edn.kind(),
            path: Vec::new(),
        });
    };
    let permutation_edn = m
        .get_keyword("permutation")
        .ok_or_else(|| DeError::MissingField {
            key: "permutation".into(),
            path: vec!["fluxionality".into()],
        })?;
    let permutation = LigandPermutation(read_edn_permutation(permutation_edn, kind.degree())?);
    let active = match m.get_keyword("active") {
        None => BooleanForm::Lit(true),
        Some(edn) => BooleanDsl::from_edn(edn)?.0,
    };
    Ok(FluxionalityForm {
        permutation,
        active,
    })
}

/// Molecule-scope DSL wrapper for a stereo constraint. It carries the element
/// kind (the stereo subtype) so the permutation degree is known when parsing.
/// The EDN is a positional 2-vector `[<kind> {<constraint-key> <value>}]`,
/// mirroring the `(StereoKind, _)` tuple: kind first (container-fixed position),
/// then the single-key constraint payload. Self-contained, so the generic
/// entity-leaf machinery applies, and `kind` is readable before the value.
macro_rules! stereo_constraint_dsl {
    ($dsl:ident, $constraint:ident, $context:literal) => {
        #[derive(Clone, Debug, PartialEq, Eq, Hash)]
        pub struct $dsl(pub StereoKind, pub $constraint);

        impl ToEdn for $dsl {
            fn to_edn(&self) -> Edn<'static> {
                let (key, value) = match &self.1 {
                    $constraint::LigandSymmetry(ls) => {
                        ("ligand-symmetry", render_edn_ligand_symmetry(ls))
                    }
                    $constraint::Fluxionality(f) => ("fluxionality", render_edn_fluxionality(f)),
                    $constraint::Topicity(t) => ("topicity", TopicityDsl(t.clone()).to_edn()),
                    $constraint::Stereogenicity(g) => {
                        ("stereogenicity", StereogenicityDsl(g.clone()).to_edn())
                    }
                };
                Edn::Vector(vec![render_edn_stereo_kind(self.0), single_key_map(key, value)].into())
            }
        }

        impl<'de> FromEdn<'de> for $dsl {
            fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
                let Edn::Vector(v) = edn else {
                    return Err(DeError::TypeMismatch {
                        expected: "stereo constraint [kind {key value}]",
                        got: edn.kind(),
                        path: vec![$context.into()],
                    });
                };
                if v.len() != 2 {
                    return Err(DeError::Custom(format!(
                        "{} must be [kind {{key value}}], got {}-element vector",
                        $context,
                        v.len()
                    )));
                }
                let kind = read_edn_stereo_kind(&v[0])?;
                let Edn::Map(m) = &v[1] else {
                    return Err(DeError::TypeMismatch {
                        expected: "single-key constraint map",
                        got: v[1].kind(),
                        path: vec![$context.into()],
                    });
                };
                if m.len() != 1 {
                    return Err(DeError::Custom(format!(
                        "{} payload must have one key, got {}",
                        $context,
                        m.len()
                    )));
                }
                let (k, value) = m.iter().next().unwrap();
                let Edn::Keyword(key) = k else {
                    return Err(DeError::TypeMismatch {
                        expected: "keyword key",
                        got: k.kind(),
                        path: vec![$context.into()],
                    });
                };
                let constraint = match key.name() {
                    "ligand-symmetry" => {
                        $constraint::LigandSymmetry(read_edn_ligand_symmetry(value, kind)?)
                    }
                    "fluxionality" => {
                        $constraint::Fluxionality(read_edn_fluxionality(value, kind)?)
                    }
                    "topicity" => $constraint::Topicity(TopicityDsl::from_edn(value)?.0),
                    "stereogenicity" => {
                        $constraint::Stereogenicity(StereogenicityDsl::from_edn(value)?.0)
                    }
                    other => {
                        return Err(DeError::Custom(format!(
                            "unknown stereo constraint keyword :{other}"
                        )))
                    }
                };
                Ok($dsl(kind, constraint))
            }
        }

        impl FromIr<$constraint> for $dsl {
            type Ctx = StereoKind;

            fn from_ir(form: &$constraint, ctx: &Self::Ctx) -> Self {
                Self(*ctx, form.clone())
            }
        }

        impl IntoIr<$constraint> for $dsl {
            type Ctx = ();

            fn into_ir(self, _ctx: &Self::Ctx) -> $constraint {
                self.1
            }
        }
    };
}

stereo_constraint_dsl! {
    StereoAtomConstraintDsl, StereoAtomConstraintForm, "stereo-atom-constraint"
}
stereo_constraint_dsl! {
    StereoBondConstraintDsl, StereoBondConstraintForm, "stereo-bond-constraint"
}

pub(crate) fn coset_lit(n: i64) -> Result<u32, DeError> {
    u32::try_from(n).map_err(|_| DeError::OutOfRange {
        value: n.to_string(),
        target: "u32",
        path: Vec::new(),
    })
}

/// Surface DSL wrapper around `StereoCoset` — the coset value under
/// `:stereo`. EDN form: int (`Lit`), `:undetermined`, a vector of ints
/// (`Expr(LitSet)`), or a string carrying the operator-expression subgrammar.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoCosetDsl(pub StereoCoset);

impl FromIr<StereoCoset> for StereoCosetDsl {
    type Ctx = ();

    fn from_ir(coset: &StereoCoset, _ctx: &Self::Ctx) -> Self {
        Self(coset.clone())
    }
}

impl IntoIr<StereoCoset> for StereoCosetDsl {
    type Ctx = ();

    fn into_ir(self, _ctx: &Self::Ctx) -> StereoCoset {
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
            Edn::Int(n) => StereoCoset::Lit(coset_lit(*n)?),
            Edn::Keyword(k) if k.name() == "undetermined" => StereoCoset::Undetermined,
            Edn::Vector(xs) => {
                let mut set = BTreeSet::new();
                for e in xs.iter() {
                    let Edn::Int(n) = e else {
                        return Err(DeError::TypeMismatch {
                            expected: "int (coset-set element)",
                            got: e.kind(),
                            path: Vec::new(),
                        });
                    };
                    set.insert(coset_lit(*n)?);
                }
                StereoCoset::LitSet(set)
            }
            Edn::Str(s) => {
                // The coset-form is the payload of a `#T` / `#C` config, both degree 4.
                parse_stereo_coset(s, 4).map_err(|e| DeError::subgrammar("stereo coset", e))?
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
            StereoCoset::Lit(n) => Edn::Int(*n as i64),
            StereoCoset::Undetermined => {
                Edn::Keyword(EdnKeyword::owned("undetermined".to_string()))
            }
            StereoCoset::LitSet(set) => Edn::Vector(
                set.iter()
                    .map(|n| Edn::Int(*n as i64))
                    .collect::<Vec<_>>()
                    .into(),
            ),
            StereoCoset::Term(_) => Edn::Str(Cow::Owned(self.to_string())),
        }
    }
}

/// Generates the constraint-side DSL for a fixed-kind stereo site (`#T`/`#C`):
/// the surface config-string parser/formatter and the EDN boundary type
/// (`:undetermined`, `:not-stereo`, or `{:stereo <coset>}`). `$kind` fixes the
/// coset degree; the per-kind type's `Stereo` arm carries the coset.
macro_rules! stereo_site_dsl {
    ($dsl:ident, $form:ident, $kind:expr, $parse:ident, $fmt:ident) => {
        pub(crate) fn $parse(i: &mut &str) -> PResult<$form> {
            alt((
                '*'.value($form::Undetermined),
                '!'.value($form::NotStereo),
                '+'.value($form::Stereo(StereoCoset::Undetermined)),
                (|i: &mut &str| stereo_coset(i, $kind.degree())).map($form::Stereo),
            ))
            .parse_next(i)
        }

        pub(crate) fn $fmt(f: &mut fmt::Formatter<'_>, config: &$form) -> fmt::Result {
            match config {
                $form::Undetermined => write!(f, "*"),
                $form::NotStereo => write!(f, "!"),
                $form::Stereo(StereoCoset::Undetermined) => write!(f, "+"),
                $form::Stereo(coset) => fmt_stereo_coset(f, coset),
            }
        }

        #[derive(Clone, Debug, PartialEq, Eq)]
        pub struct $dsl(pub $form);

        impl FromIr<$form> for $dsl {
            type Ctx = ();

            fn from_ir(form: &$form, _ctx: &Self::Ctx) -> Self {
                Self(form.clone())
            }
        }

        impl IntoIr<$form> for $dsl {
            type Ctx = ();

            fn into_ir(self, _ctx: &Self::Ctx) -> $form {
                self.0
            }
        }

        impl<'de> FromEdn<'de> for $dsl {
            fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
                match edn {
                    Edn::Keyword(k) if k.name() == "undetermined" => Ok(Self($form::Undetermined)),
                    Edn::Keyword(k) if k.name() == "not-stereo" => Ok(Self($form::NotStereo)),
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
                            "stereo" => Ok(Self($form::Stereo(
                                StereoCosetDsl::from_edn(v)?.into_ir(&()),
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

        impl ToEdn for $dsl {
            fn to_edn(&self) -> Edn<'static> {
                match &self.0 {
                    $form::Undetermined => {
                        Edn::Keyword(EdnKeyword::owned("undetermined".to_string()))
                    }
                    $form::NotStereo => Edn::Keyword(EdnKeyword::owned("not-stereo".to_string())),
                    $form::Stereo(coset) => {
                        single_key_map("stereo", StereoCosetDsl::from_ir(coset, &()).to_edn())
                    }
                }
            }
        }
    };
}

stereo_site_dsl! {
    TetrahedralStereoDsl,
    TetrahedralStereoForm,
    StereoKind::Tetrahedral,
    tetrahedral_stereo_config,
    fmt_tetrahedral_stereo_config
}
stereo_site_dsl! {
    CisTransStereoDsl,
    CisTransStereoForm,
    StereoKind::CisTrans,
    cis_trans_stereo_config,
    fmt_cis_trans_stereo_config
}

#[cfg(test)]
mod tests {
    use std::fs;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_edn::read_string;

    use super::*;

    /// Every `fuzz_entity_strings` stereo seed must parse with its stereo parser — guards the seed
    /// corpus against notation rot.
    #[rstest]
    fn test_fuzz_entity_strings_stereo_seeds_valid() {
        let dir = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/fuzz/seeds/fuzz_entity_strings"
        );
        let mut failures: Vec<String> = Vec::new();
        for entry in fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            let name = path.file_name().unwrap().to_string_lossy().into_owned();
            let data = fs::read_to_string(&path).unwrap();
            let result = if name.starts_with("stereo_atom_") {
                parse_stereo_atom(&data).map(|_| ())
            } else if name.starts_with("stereo_bond_") {
                parse_stereo_bond(&data).map(|_| ())
            } else {
                continue;
            };
            if let Err(e) = result {
                failures.push(format!("{name}: {e:?}"));
            }
        }
        assert!(
            failures.is_empty(),
            "invalid stereo seeds:\n{}",
            failures.join("\n")
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::tetrahedral_ccw("Th0", StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))))]
    #[case::tetrahedral_cw("Th1", StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1))))]
    #[case::open("Th*", StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined)))]
    #[case::square_planar("Sp2", StereoAtomDsl(StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Lit(2))))]
    #[case::octahedral("Oh6", StereoAtomDsl(StereoAtomForm::new(StereoKind::Octahedral, StereoCoset::Lit(6))))]
    #[case::undetermined("*", StereoAtomDsl(StereoAtomForm::default()))]
    #[case::no_canonicalization("Th~1", StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(1))))))]
    fn test_parse_stereo_atom(#[case] input: &str, #[case] expected: StereoAtomDsl) {
        assert_eq!(parse_stereo_atom(input).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::not_stereo("Th!", ParseError::Syntax)]
    #[case::trailing_after_coset("Th1x", ParseError::TrailingInput("x".to_string()))]
    #[case::unknown_predicate("Th1#x", ParseError::UnknownStereoPredicate("#x".to_string()))]
    #[case::duplicate("Th1#g/#g=", ParseError::DuplicateStereoPredicate("#g".to_string()))]
    #[case::undetermined_rejects_constraint("*#g/", ParseError::TrailingInput("#g/".to_string()))]
    fn test_parse_stereo_atom_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(parse_stereo_atom(input).unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", StereoAtomUpdateDsl::default())]
    #[case::undetermined("*", StereoAtomUpdateDsl(StereoAtomUpdate { configuration: StereoConfigurationUpdate::Undetermined, ..Default::default() }))]
    #[case::absolute("Th1", StereoAtomUpdateDsl(StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: Some(StereoCoset::Lit(1)) }, ..Default::default() }))]
    #[case::relative("Th", StereoAtomUpdateDsl(StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: None }, ..Default::default() }))]
    #[case::explicit_open("Th*", StereoAtomUpdateDsl(StereoAtomUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::Tetrahedral, coset: Some(StereoCoset::Undetermined) }, ..Default::default() }))]
    fn test_parse_stereo_atom_update(#[case] input: &str, #[case] expected: StereoAtomUpdateDsl) {
        assert_eq!(parse_stereo_atom_update(input).unwrap(), expected);
    }

    #[rstest]
    #[case::undetermined_kind_with_constraint("*#o(0,1)=", ParseError::TrailingInput("#o(0,1)=".to_string()))]
    fn test_parse_stereo_atom_update_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(parse_stereo_atom_update(input).unwrap_err(), expected);
    }

    #[rstest]
    #[case::undetermined_kind_with_constraint(
        "*#o(0,1)=",
        ParseError::TrailingInput("#o(0,1)=".to_string())
    )]
    fn test_stereo_atom_update_from_str_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(input.parse::<StereoAtomUpdate>().unwrap_err(), expected);
    }

    #[rstest]
    #[case::empty("")]
    #[case::undetermined("*")]
    #[case::absolute("Th1")]
    #[case::relative("Th")]
    #[case::explicit_open("Th*")]
    #[case::topicity_removal("Th#o(0,1)*")]
    #[case::topicity_change("Th#o(0,1)/")]
    #[case::ligand_symmetry_removal("Th#p(0,1)*")]
    #[case::ligand_symmetry("Th#p(0,1)")]
    fn test_stereo_atom_update_dsl_display_roundtrip(#[case] input: &str) {
        let dsl = parse_stereo_atom_update(input).unwrap();
        assert_eq!(dsl.to_string(), input);
        assert_eq!(parse_stereo_atom_update(&dsl.to_string()).unwrap(), dsl);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", StereoBondUpdateDsl::default())]
    #[case::undetermined("*", StereoBondUpdateDsl(StereoBondUpdate { configuration: StereoConfigurationUpdate::Undetermined, ..Default::default() }))]
    #[case::absolute("Ct1", StereoBondUpdateDsl(StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: Some(StereoCoset::Lit(1)) }, ..Default::default() }))]
    #[case::relative("Ct", StereoBondUpdateDsl(StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: None }, ..Default::default() }))]
    #[case::explicit_open("Ct*", StereoBondUpdateDsl(StereoBondUpdate { configuration: StereoConfigurationUpdate::Kinded { kind: StereoKind::CisTrans, coset: Some(StereoCoset::Undetermined) }, ..Default::default() }))]
    fn test_parse_stereo_bond_update(#[case] input: &str, #[case] expected: StereoBondUpdateDsl) {
        assert_eq!(parse_stereo_bond_update(input).unwrap(), expected);
    }

    #[rstest]
    #[case::undetermined_kind_with_constraint("*#g/", ParseError::TrailingInput("#g/".to_string()))]
    fn test_parse_stereo_bond_update_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(parse_stereo_bond_update(input).unwrap_err(), expected);
    }

    #[rstest]
    #[case::undetermined_kind_with_constraint(
        "*#g/",
        ParseError::TrailingInput("#g/".to_string())
    )]
    fn test_stereo_bond_update_from_str_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(input.parse::<StereoBondUpdate>().unwrap_err(), expected);
    }

    #[rstest]
    #[case::empty_cycle("()", Permutation::identity(4))]
    #[case::single_cycle("(0,1,2)", Permutation::from_image(&[1, 2, 0, 3]))]
    #[case::disjoint_cycles("(0,1)(2,3)", Permutation::from_image(&[1, 0, 3, 2]))]
    fn test_parse_permutation(#[case] input: &str, #[case] expected: Permutation) {
        assert_eq!(parse_permutation(input, 4), Ok(expected));
    }

    #[rstest]
    #[case::overlap(
        "(0,1)(1,2)",
        ParseError::InvalidValue("cycle point 1 occurs more than once".to_string()),
    )]
    #[case::repeated(
        "(0,1,0)",
        ParseError::InvalidValue("cycle point 0 occurs more than once".to_string()),
    )]
    #[case::out_of_range(
        "(0,4)",
        ParseError::InvalidValue(
            "cycle point 4 at cycle 0, position 1 is outside 0..4".to_string(),
        ),
    )]
    fn test_parse_permutation_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(parse_permutation(input, 4), Err(expected));
    }

    #[rstest]
    #[case::empty("")]
    #[case::undetermined("*")]
    #[case::absolute("Ct1")]
    #[case::relative("Ct")]
    #[case::explicit_open("Ct*")]
    #[case::topicity_removal("Ct#o(0,1)*")]
    #[case::stereogenicity_change("Ct#g/")]
    fn test_stereo_bond_update_dsl_display_roundtrip(#[case] input: &str) {
        let dsl = parse_stereo_bond_update(input).unwrap();
        assert_eq!(dsl.to_string(), input);
        assert_eq!(parse_stereo_bond_update(&dsl.to_string()).unwrap(), dsl);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::tetrahedral_ccw(StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))), "Th0")]
    #[case::open(StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined)), "Th*")]
    #[case::square_planar(StereoAtomDsl(StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Lit(2))), "Sp2")]
    #[case::octahedral(StereoAtomDsl(StereoAtomForm::new(StereoKind::Octahedral, StereoCoset::Lit(6))), "Oh6")]
    fn test_stereo_atom_dsl_to_string(#[case] form: StereoAtomDsl, #[case] expected: &str) {
        assert_eq!(form.to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fluxionality("Th1#f(0,1,2)")]
    #[case::ligand_symmetry_involution("Th1#p~")]
    #[case::ligand_symmetry_not_present("Th1#p~!")]
    #[case::ligand_symmetry_explicit("Th1#p(0,1,2)")]
    #[case::topicity("Th1#o(0,1)=")]
    #[case::topicity_negated("Th1#o(0,1)!'")]
    #[case::topicity_lit_set("Th1#o(0,1){=,'}")]
    #[case::topicity_not_set("Th1#o(0,1)!{=,'}")]
    #[case::stereogenicity("Th1#g/")]
    #[case::stereogenicity_lit_set("Th1#g{=,/}")]
    #[case::multiple("Th1#f(0,1,2)#o(0,1)=#g/")]
    fn test_stereo_atom_inline_render_identity(#[case] s: &str) {
        assert_eq!(parse_stereo_atom(s).unwrap().to_string(), s);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::topicity_open("Th1#o(0,1)*", "Th1")]
    #[case::stereogenicity_open("Th1#g*", "Th1")]
    fn test_stereo_atom_inline_render(#[case] input: &str, #[case] canonical: &str) {
        assert_eq!(parse_stereo_atom(input).unwrap().to_string(), canonical);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fluxionality("Th1#f(0,1,2)",
        StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 2, 0, 3])), active: BooleanForm::Lit(true) }))]
    #[case::ligand_symmetry("Th1#p(0,1,2)",
        StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
            permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 2, 0, 3])), orientation: Orientation::Proper },
            invariant: BooleanForm::Lit(true) }))]
    #[case::ligand_symmetry_absent("Th1#p(0,1,2)!",
        StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
            permutation: OrientedLigandPermutation { permutation: LigandPermutation(Permutation::from_image(&[1, 2, 0, 3])), orientation: Orientation::Proper },
            invariant: BooleanForm::Lit(false) }))]
    #[case::fluxionality_absent("Th1#f(0,1,2)!",
        StereoAtomConstraintForm::Fluxionality(FluxionalityForm { permutation: LigandPermutation(Permutation::from_image(&[1, 2, 0, 3])), active: BooleanForm::Lit(false) }))]
    #[case::topicity_negated("Th1#o(0,1)!'",
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Enantiotopic])) }))]
    #[case::topicity_lit_set("Th1#o(0,1){=,'}",
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::LitSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])) }))]
    #[case::topicity_not_set("Th1#o(0,1)!{=,'}",
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])) }))]
    #[case::topicity_open("Th1#o(0,1)*",
        StereoAtomConstraintForm::Topicity(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined }))]
    #[case::stereogenicity("Th1#g/",
        StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)))]
    fn test_stereo_atom_predicate(#[case] input: &str, #[case] expected: StereoAtomConstraintForm) {
        let dsl = parse_stereo_atom(input).unwrap();
        assert_eq!(dsl.0.constraints.iter().cloned().collect::<Vec<_>>(), vec![expected]);
    }

    #[rstest]
    fn test_stereo_atom_predicate_involution() {
        let dsl = parse_stereo_atom("Th1#p~").unwrap();
        let expected = StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
            permutation: OrientedLigandPermutation {
                permutation: LigandPermutation(StereoKind::Tetrahedral.involution()),
                orientation: Orientation::Improper,
            },
            invariant: BooleanForm::Lit(true),
        });
        assert_eq!(
            dsl.0.constraints.iter().cloned().collect::<Vec<_>>(),
            vec![expected],
        );
    }

    #[rstest]
    #[case::empty_cycle("[[]]", Permutation::identity(4))]
    #[case::single_cycle("[[0 1 2]]", Permutation::from_image(&[1, 2, 0, 3]))]
    #[case::disjoint_cycles("[[0 1] [2 3]]", Permutation::from_image(&[1, 0, 3, 2]))]
    fn test_read_edn_permutation(#[case] input: &str, #[case] expected: Permutation) {
        assert_eq!(
            read_edn_permutation(&read_string(input).unwrap(), 4),
            Ok(expected),
        );
    }

    #[rstest]
    #[case::overlap(
        "[[0 1] [1 2]]",
        DeError::Custom("cycle point 1 occurs more than once".to_string()),
    )]
    #[case::repeated(
        "[[0 1 0]]",
        DeError::Custom("cycle point 0 occurs more than once".to_string()),
    )]
    #[case::out_of_range(
        "[[0 4]]",
        DeError::Custom(
            "cycle point 4 at cycle 0, position 1 is outside 0..4".to_string(),
        ),
    )]
    #[case::negative(
        "[[0 -1]]",
        DeError::OutOfRange {
            value: "-1".to_string(),
            target: "ligand position",
            path: Vec::new(),
        },
    )]
    fn test_read_edn_permutation_error(#[case] input: &str, #[case] expected: DeError) {
        assert_eq!(
            read_edn_permutation(&read_string(input).unwrap(), 4),
            Err(expected),
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::string("\"Th1\"", StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1))))]
    #[case::keyword_ccw(":ccw", StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))))]
    #[case::keyword_cw(":cw", StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1))))]
    #[case::string_square_planar("\"Sp2\"", StereoAtomDsl(StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Lit(2))))]
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
    #[case::canonical_ccw(StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))), ":ccw")]
    #[case::canonical_cw(StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1))), ":cw")]
    #[case::open_string(StereoAtomDsl(StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined)), "\"Th*\"")]
    #[case::non_tetrahedral_string(StereoAtomDsl(StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Lit(1))), "\"Sp1\"")]
    fn test_stereo_atom_dsl_to_edn(#[case] form: StereoAtomDsl, #[case] expected: &str) {
        assert_eq!(form.to_edn(), read_string(expected).unwrap());
    }

    #[rstest]
    fn test_stereo_atom_dsl_into_ast() {
        let form = StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Undetermined);
        assert_eq!(
            StereoAtomDsl(form.clone()).into_ir(&StereoAtomDefaults::default()),
            form
        );
        assert_eq!(
            StereoAtomDsl::from_ir(&form, &StereoAtomDefaults::default()),
            StereoAtomDsl(form)
        );
    }

    #[rstest]
    #[case::ccw("ccw", Some("Th0"))]
    #[case::cw("cw", Some("Th1"))]
    #[case::unknown("xyz", None)]
    fn test_expand_stereo_atom_keyword(#[case] name: &str, #[case] expected: Option<&str>) {
        assert_eq!(expand_stereo_atom_keyword(name), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::cis_trans_z("Ct0", StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0))))]
    #[case::cis_trans_e("Ct1", StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1))))]
    #[case::open("Ct*", StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined)))]
    #[case::undetermined("*", StereoBondDsl(StereoBondForm::default()))]
    #[case::no_canonicalization("Ct~1", StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(1))))))]
    fn test_parse_stereo_bond(#[case] input: &str, #[case] expected: StereoBondDsl) {
        assert_eq!(parse_stereo_bond(input).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::not_stereo("Ct!", ParseError::Syntax)]
    #[case::trailing_after_coset("Ct1x", ParseError::TrailingInput("x".to_string()))]
    #[case::unknown_predicate("Ct1#x", ParseError::UnknownStereoPredicate("#x".to_string()))]
    #[case::duplicate("Ct1#g/#g=", ParseError::DuplicateStereoPredicate("#g".to_string()))]
    #[case::undetermined_rejects_constraint("*#g/", ParseError::TrailingInput("#g/".to_string()))]
    fn test_parse_stereo_bond_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(parse_stereo_bond(input).unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::cis_trans_z(StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0))), "Ct0")]
    #[case::open(StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined)), "Ct*")]
    fn test_stereo_bond_dsl_to_string(#[case] form: StereoBondDsl, #[case] expected: &str) {
        assert_eq!(form.to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fluxionality_involution("Ct1#f~")]
    #[case::topicity("Ct1#o(0,1)=")]
    #[case::stereogenicity("Ct1#g/")]
    fn test_stereo_bond_inline_render_identity(#[case] s: &str) {
        assert_eq!(parse_stereo_bond(s).unwrap().to_string(), s);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::string("\"Ct1\"", StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1))))]
    #[case::keyword_z(":z", StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0))))]
    #[case::keyword_e(":e", StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1))))]
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
    #[case::canonical_z(StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0))), ":z")]
    #[case::canonical_e(StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1))), ":e")]
    #[case::open_string(StereoBondDsl(StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Undetermined)), "\"Ct*\"")]
    fn test_stereo_bond_dsl_to_edn(#[case] form: StereoBondDsl, #[case] expected: &str) {
        assert_eq!(form.to_edn(), read_string(expected).unwrap());
    }

    #[rstest]
    #[case::z("z", Some("Ct0"))]
    #[case::e("e", Some("Ct1"))]
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
    #[case::undetermined("*", TetrahedralStereoForm::Undetermined)]
    #[case::not_stereo("!", TetrahedralStereoForm::NotStereo)]
    #[case::stereo_undetermined("+", TetrahedralStereoForm::Stereo(StereoCoset::Undetermined))]
    #[case::lit("1", TetrahedralStereoForm::Stereo(StereoCoset::Lit(1)))]
    #[case::var("?o", TetrahedralStereoForm::Stereo(StereoCoset::term(StereoTerm::var("o"))))]
    #[case::lit_set("{1,2}", TetrahedralStereoForm::Stereo(StereoCoset::lit_set([1, 2])))]
    #[case::var_domain("?o :: {1,2}", TetrahedralStereoForm::Stereo(StereoCoset::term(StereoTerm::var_in("o", [1, 2]))))]
    #[case::swap("~1", TetrahedralStereoForm::Stereo(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(1)))))]
    #[case::mirror("'1", TetrahedralStereoForm::Stereo(StereoCoset::term(StereoTerm::mirror(StereoTerm::Lit(1)))))]
    #[case::apply("1^(0,1)", TetrahedralStereoForm::Stereo(StereoCoset::term(StereoTerm::apply(StereoTerm::Lit(1), Permutation::from_image(&[1, 0, 2, 3])))))]
    #[case::swap_binds_tighter_than_apply("~1^(0,1)", TetrahedralStereoForm::Stereo(StereoCoset::term(StereoTerm::apply(StereoTerm::swap(StereoTerm::Lit(1)), Permutation::from_image(&[1, 0, 2, 3])))))]
    #[case::mirror_binds_tighter_than_apply("'1^(0,1)", TetrahedralStereoForm::Stereo(StereoCoset::term(StereoTerm::apply(StereoTerm::mirror(StereoTerm::Lit(1)), Permutation::from_image(&[1, 0, 2, 3])))))]
    #[case::whitespace_ignored("  ?o :: { 1 , 2 }", TetrahedralStereoForm::Stereo(StereoCoset::term(StereoTerm::var_in("o", [1, 2]))))]
    #[case::no_canonicalization("{1}", TetrahedralStereoForm::Stereo(StereoCoset::lit_set([1])))]
    fn test_tetrahedral_stereo_config(#[case] input: &str, #[case] expected: TetrahedralStereoForm) {
        assert_eq!(
            (|i: &mut &str| tetrahedral_stereo_config(i)).parse(input).unwrap(),
            expected
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined("*", CisTransStereoForm::Undetermined)]
    #[case::not_stereo("!", CisTransStereoForm::NotStereo)]
    #[case::stereo_undetermined("+", CisTransStereoForm::Stereo(StereoCoset::Undetermined))]
    #[case::lit("1", CisTransStereoForm::Stereo(StereoCoset::Lit(1)))]
    #[case::no_canonicalization("~1", CisTransStereoForm::Stereo(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(1)))))]
    fn test_cis_trans_stereo_config(#[case] input: &str, #[case] expected: CisTransStereoForm) {
        assert_eq!(
            (|i: &mut &str| cis_trans_stereo_config(i)).parse(input).unwrap(),
            expected
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(TetrahedralStereoForm::Undetermined, "*")]
    #[case::not_stereo(TetrahedralStereoForm::NotStereo, "!")]
    #[case::stereogenic(TetrahedralStereoForm::Stereo(StereoCoset::Undetermined), "+")]
    #[case::lit(TetrahedralStereoForm::Stereo(StereoCoset::Lit(1)), "1")]
    #[case::var(TetrahedralStereoForm::Stereo(StereoCoset::term(StereoTerm::var("o"))), "?o")]
    #[case::lit_set(TetrahedralStereoForm::Stereo(StereoCoset::lit_set([1, 2])), "{1,2}")]
    #[case::var_domain(TetrahedralStereoForm::Stereo(StereoCoset::term(StereoTerm::var_in("o", [1, 2]))), "?o :: {1,2}")]
    #[case::swap(TetrahedralStereoForm::Stereo(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(1)))), "~1")]
    #[case::mirror(TetrahedralStereoForm::Stereo(StereoCoset::term(StereoTerm::mirror(StereoTerm::Lit(1)))), "'1")]
    #[case::apply(TetrahedralStereoForm::Stereo(StereoCoset::term(StereoTerm::apply(StereoTerm::Lit(1), Permutation::from_image(&[1, 0, 2, 3])))), "1^(0,1)")]
    fn test_fmt_tetrahedral_stereo_config(#[case] c: TetrahedralStereoForm, #[case] expected: &str) {
        struct W(TetrahedralStereoForm);
        impl fmt::Display for W {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                fmt_tetrahedral_stereo_config(f, &self.0)
            }
        }
        assert_eq!(W(c).to_string(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(StereoCoset::Lit(2))]
    #[case::undetermined(StereoCoset::Undetermined)]
    #[case::lit_set(StereoCoset::lit_set([1, 2]))]
    #[case::term_swap(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(1))))]
    fn test_stereo_coset_dsl_into_ast(#[case] coset: StereoCoset) {
        assert_eq!(StereoCosetDsl(coset.clone()).into_ir(&()), coset);
        assert_eq!(StereoCosetDsl::from_ir(&coset, &()), StereoCosetDsl(coset));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::int("2", StereoCosetDsl(StereoCoset::Lit(2)))]
    #[case::undetermined(":undetermined", StereoCosetDsl(StereoCoset::Undetermined))]
    #[case::vector("[1 2]", StereoCosetDsl(StereoCoset::lit_set([1, 2])))]
    #[case::string_lit("\"3\"", StereoCosetDsl(StereoCoset::Lit(3)))]
    #[case::string_term("\"~1\"", StereoCosetDsl(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(1)))))]
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
    #[case::lit(StereoCosetDsl(StereoCoset::Lit(2)), "2")]
    #[case::undetermined(StereoCosetDsl(StereoCoset::Undetermined), ":undetermined")]
    #[case::lit_set(StereoCosetDsl(StereoCoset::lit_set([1, 2])), "[1 2]")]
    #[case::term_swap(StereoCosetDsl(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(1)))), "\"~1\"")]
    fn test_stereo_coset_dsl_to_edn(#[case] form: StereoCosetDsl, #[case] expected: &str) {
        assert_eq!(form.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(TetrahedralStereoForm::Undetermined)]
    #[case::not_stereo(TetrahedralStereoForm::NotStereo)]
    #[case::stereo_lit(TetrahedralStereoForm::Stereo(StereoCoset::Lit(1)))]
    fn test_tetrahedral_stereo_dsl_into_ast(#[case] form: TetrahedralStereoForm) {
        assert_eq!(TetrahedralStereoDsl(form.clone()).into_ir(&()), form);
        assert_eq!(TetrahedralStereoDsl::from_ir(&form, &()), TetrahedralStereoDsl(form));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(":undetermined", TetrahedralStereoDsl(TetrahedralStereoForm::Undetermined))]
    #[case::not_stereo(":not-stereo", TetrahedralStereoDsl(TetrahedralStereoForm::NotStereo))]
    #[case::stereo_lit("{:stereo 1}", TetrahedralStereoDsl(TetrahedralStereoForm::Stereo(StereoCoset::Lit(1))))]
    #[case::stereo_undetermined("{:stereo :undetermined}", TetrahedralStereoDsl(TetrahedralStereoForm::Stereo(StereoCoset::Undetermined)))]
    #[case::stereo_set("{:stereo [1 2]}", TetrahedralStereoDsl(TetrahedralStereoForm::Stereo(StereoCoset::lit_set([1, 2]))))]
    fn test_tetrahedral_stereo_dsl_from_edn(#[case] input: &str, #[case] expected: TetrahedralStereoDsl) {
        assert_eq!(TetrahedralStereoDsl::from_edn(&read_string(input).unwrap()).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::unknown_keyword(":bogus", DeError::TypeMismatch { expected: ":undetermined / :not-stereo / {:stereo <coset>}", got: "keyword", path: Vec::new() })]
    #[case::unknown_key("{:bogus 1}", DeError::UnknownField { key: "bogus".to_string(), path: vec!["stereo-configuration".into()] })]
    #[case::wrong_type("1", DeError::TypeMismatch { expected: ":undetermined / :not-stereo / {:stereo <coset>}", got: "int", path: Vec::new() })]
    fn test_tetrahedral_stereo_dsl_from_edn_error(#[case] input: &str, #[case] expected: DeError) {
        assert_eq!(TetrahedralStereoDsl::from_edn(&read_string(input).unwrap()).unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(TetrahedralStereoDsl(TetrahedralStereoForm::Undetermined), ":undetermined")]
    #[case::not_stereo(TetrahedralStereoDsl(TetrahedralStereoForm::NotStereo), ":not-stereo")]
    #[case::stereo_lit(TetrahedralStereoDsl(TetrahedralStereoForm::Stereo(StereoCoset::Lit(1))), "{:stereo 1}")]
    #[case::stereo_undetermined(TetrahedralStereoDsl(TetrahedralStereoForm::Stereo(StereoCoset::Undetermined)), "{:stereo :undetermined}")]
    fn test_tetrahedral_stereo_dsl_to_edn(#[case] form: TetrahedralStereoDsl, #[case] expected: &str) {
        assert_eq!(form.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(TopicityDsl(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }), "{:pair [0 1] :relation :homotopic}")]
    #[case::undetermined(TopicityDsl(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined }), "{:pair [0 1] :relation :undetermined}")]
    #[case::lit_set(TopicityDsl(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(1), StereoLigandPosition(2)), relation: TopicityRelationForm::LitSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])) }), "{:pair [1 2] :relation [:homotopic :enantiotopic]}")]
    #[case::not_set(TopicityDsl(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(1), StereoLigandPosition(2)), relation: TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])) }), "{:pair [1 2] :relation {:not-in [:diastereotopic]}}")]
    fn test_topicity_dsl_to_edn(#[case] form: TopicityDsl, #[case] expected: &str) {
        assert_eq!(form.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit("{:pair [0 1] :relation :homotopic}", TopicityDsl(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Lit(Topicity::Homotopic) }))]
    #[case::undetermined("{:pair [0 1] :relation :undetermined}", TopicityDsl(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)), relation: TopicityRelationForm::Undetermined }))]
    #[case::lit_set("{:pair [1 2] :relation [:homotopic :enantiotopic]}", TopicityDsl(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(1), StereoLigandPosition(2)), relation: TopicityRelationForm::LitSet(BTreeSet::from([Topicity::Homotopic, Topicity::Enantiotopic])) }))]
    #[case::not_set("{:pair [1 2] :relation {:not-in [:diastereotopic]}}", TopicityDsl(TopicityForm { pair: StereoLigandPair::new(StereoLigandPosition(1), StereoLigandPosition(2)), relation: TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])) }))]
    fn test_topicity_dsl_from_edn(#[case] input: &str, #[case] expected: TopicityDsl) {
        assert_eq!(TopicityDsl::from_edn(&read_string(input).unwrap()).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wrong_type("nil", DeError::TypeMismatch { expected: "topicity map", got: "nil", path: Vec::new() })]
    #[case::missing_pair("{:relation :homotopic}", DeError::MissingField { key: "pair".into(), path: vec!["topicity".into()] })]
    #[case::missing_relation("{:pair [0 1]}", DeError::MissingField { key: "relation".into(), path: vec!["topicity".into()] })]
    fn test_topicity_dsl_from_edn_error(#[case] input: &str, #[case] expected: DeError) {
        assert_eq!(TopicityDsl::from_edn(&read_string(input).unwrap()).unwrap_err(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit(StereogenicityDsl(StereogenicityForm::Lit(Stereogenicity::Stereogenic)), "{:relation :stereogenic}")]
    #[case::undetermined(StereogenicityDsl(StereogenicityForm::Undetermined), "{:relation :undetermined}")]
    #[case::not_set(StereogenicityDsl(StereogenicityForm::NotSet(BTreeSet::from([Stereogenicity::Symmetric]))), "{:relation {:not-in [:symmetric]}}")]
    fn test_stereogenicity_dsl_to_edn(#[case] form: StereogenicityDsl, #[case] expected: &str) {
        assert_eq!(form.to_edn(), read_string(expected).unwrap());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::lit("{:relation :stereogenic}", StereogenicityDsl(StereogenicityForm::Lit(Stereogenicity::Stereogenic)))]
    #[case::undetermined("{:relation :undetermined}", StereogenicityDsl(StereogenicityForm::Undetermined))]
    #[case::lit_set("{:relation [:prochiral :stereogenic]}", StereogenicityDsl(StereogenicityForm::LitSet(BTreeSet::from([Stereogenicity::Prochiral, Stereogenicity::Stereogenic]))))]
    #[case::not_set("{:relation {:not-in [:symmetric]}}", StereogenicityDsl(StereogenicityForm::NotSet(BTreeSet::from([Stereogenicity::Symmetric]))))]
    fn test_stereogenicity_dsl_from_edn(#[case] input: &str, #[case] expected: StereogenicityDsl) {
        assert_eq!(StereogenicityDsl::from_edn(&read_string(input).unwrap()).unwrap(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wrong_type("nil", DeError::TypeMismatch { expected: "stereogenicity map", got: "nil", path: Vec::new() })]
    #[case::unknown_keyword("{:relation :bogus}", DeError::TypeMismatch { expected: "StereogenicityForm keyword", got: "keyword", path: vec!["stereogenicity".into()] })]
    fn test_stereogenicity_dsl_from_edn_error(#[case] input: &str, #[case] expected: DeError) {
        assert_eq!(StereogenicityDsl::from_edn(&read_string(input).unwrap()).unwrap_err(), expected);
    }
}
