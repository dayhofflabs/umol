//! Atom-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_chem::element::Element;
use umol_edn::{DeError, Edn, EdnError, EdnKeyword, EdnStreamDeserializer, FromEdn, ToEdn};
use winnow::ascii::{dec_uint, multispace0};
use winnow::combinator::{alt, delimited, empty, opt, preceded, repeat, separated, terminated};
use winnow::error::{ErrMode, ParserError};
use winnow::token::{one_of, take};
use winnow::Parser;

use super::config::{
    AromaticValenceDefault, AtomDefaults, IsotopeDefault, MulticenterValenceDefault,
    NumericDefault, StereoDefault,
};
use super::constraint::RingMembershipDsl;
use super::edn_utils::single_key_map;
use super::error::{PResult, ParseError};
use super::predicate::{
    apply_spin_pair, charge, fmt_charge, fmt_ring_membership, fmt_spin_pair, lower_spin,
    optional_value, raise_spin, ring_membership, SpinPredicate,
};
use super::stereo::{
    fmt_tetrahedral_stereo_config, tetrahedral_stereo_config, TetrahedralStereoDsl,
};
use super::value::{fmt_set, fmt_value, id, terminator, value, ValueDsl};
use crate::ast::atom::{AtomAst, ElementAst, IsotopeMassAst};
use crate::ast::constraint::{
    AromaticValenceAst, AtomConstraintAst, AtomConstraintKey, AtomConstraintsAst,
    MulticenterValenceAst,
};
use crate::ast::operators::MemOp;
use crate::ast::stereo::TetrahedralStereoAst;
use crate::ast::traits::{FromAst, IntoAst, Lattice};
use crate::ast::value::ValueAst;

/// Surface DSL wrapper around `AtomAst`. Parses and renders the atom-string form
/// (element plus `#…` predicates); inline-capable constraints land in
/// `self.0.constraints`.
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AtomDsl(pub AtomAst);

impl AtomDsl {
    /// Zero-cost reference cast from `&AtomAst`. Relies on `repr(transparent)`.
    pub fn from_ref(ast: &AtomAst) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(ast as *const AtomAst as *const Self) }
    }
}

impl From<AtomAst> for AtomDsl {
    fn from(ast: AtomAst) -> Self {
        Self(ast)
    }
}

impl FromStr for AtomDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_atom(s)
    }
}

impl Display for AtomDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_atom_ast(f, &self.0)?;
        for c in self.0.constraints.iter() {
            fmt_constraint(f, c)?;
        }
        Ok(())
    }
}

impl<'de> FromEdn<'de> for AtomDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s.parse().map_err(|e| DeError::subgrammar("atom", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        EdnStreamDeserializer::new(input).read_subgrammar_all("atom")
    }
}

impl ToEdn for AtomDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

impl FromAst<AtomAst> for AtomDsl {
    type Ctx = AtomDefaults;

    fn from_ast(ast: &AtomAst, cfg: &Self::Ctx) -> Self {
        let mut out = ast.clone();
        lower_atom(&mut out, cfg);
        AtomDsl(out)
    }
}

impl IntoAst<AtomAst> for AtomDsl {
    type Ctx = AtomDefaults;

    fn into_ast(mut self, cfg: &Self::Ctx) -> AtomAst {
        raise_atom(&mut self.0, cfg);
        self.0
    }
}

impl FromStr for AtomAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(AtomDsl::from_str(s)?.into_ast(&AtomDefaults::default()))
    }
}

impl Display for AtomAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        AtomDsl::from_ref(self).fmt(f)
    }
}

impl<'de> FromEdn<'de> for AtomAst {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        Ok(AtomDsl::from_edn(edn)?.into_ast(&AtomDefaults::default()))
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        Ok(AtomDsl::from_edn_str(input)?.into_ast(&AtomDefaults::default()))
    }
}

impl ToEdn for AtomAst {
    fn to_edn(&self) -> Edn<'static> {
        AtomDsl::from_ref(self).to_edn()
    }
}

/// Parse a complete atom-string into an `AtomDsl`.
pub fn parse_atom(input: &str) -> Result<AtomDsl, ParseError> {
    if let Some(dsl) = parse_bare_element(input) {
        return Ok(dsl);
    }
    atom.parse(input).map_err(|e| e.into_inner())
}

/// Fast path for element-only atom string. Currently requires element names in title case.
fn parse_bare_element(input: &str) -> Option<AtomDsl> {
    let bytes = input.as_bytes();
    if bytes.len() != 1 && bytes.len() != 2 {
        return None;
    }
    if !bytes[0].is_ascii_uppercase() {
        return None;
    }
    if bytes.len() == 2 && !bytes[1].is_ascii_lowercase() {
        return None;
    }
    let el = Element::from_symbol_bytes(input.as_bytes())?;
    Some(AtomDsl(AtomAst::from_element(el)))
}

/// Atom-string parser (does not require consuming all input).
pub(crate) fn atom(i: &mut &str) -> PResult<AtomDsl> {
    let el = delimited(multispace0, element, multispace0).parse_next(i)?;
    let preds: Vec<AtomPredicate> =
        repeat(0.., terminated(atom_predicate, multispace0)).parse_next(i)?;
    let mut form = AtomDsl(AtomAst::new(el));
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

/// Partial atom with all fields, including element, optional.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PartialAtomDsl(pub AtomAst);

impl FromStr for PartialAtomDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_partial_atom(s)
    }
}

/// Render a partial atom-string (omits the element when `Undetermined`).
/// Undetermined constraints are rendered explicitly as `#<tag>*`, unlike
/// AtomDsl which elides them.
impl Display for PartialAtomDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let ast = &self.0;
        if !matches!(ast.element, ElementAst::Undetermined) {
            fmt_element(f, &ast.element)?;
        }
        fmt_isotope_mass(f, &ast.isotope_mass)?;
        fmt_charge(f, &ast.charge)?;
        fmt_value_field(f, "#h", &ast.implicit_hydrogens)?;
        fmt_value_field(f, "#n", &ast.lone_pairs)?;
        fmt_spin_pair(f, &ast.spin)?;
        for c in ast.constraints.iter() {
            if c.is_undetermined() {
                write!(f, "{}*", constraint_tag(c.key()))?;
            } else {
                fmt_constraint(f, c)?;
            }
        }
        Ok(())
    }
}

impl<'de> FromEdn<'de> for PartialAtomDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s.parse().map_err(|e| DeError::subgrammar("atom", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for PartialAtomDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

pub fn parse_partial_atom(input: &str) -> Result<PartialAtomDsl, ParseError> {
    partial_atom.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn partial_atom(i: &mut &str) -> PResult<PartialAtomDsl> {
    let el = delimited(multispace0, opt(element), multispace0)
        .parse_next(i)?
        .unwrap_or(ElementAst::Undetermined);
    let preds: Vec<AtomPredicate> =
        repeat(0.., terminated(atom_predicate, multispace0)).parse_next(i)?;
    let mut form = AtomDsl(AtomAst::new(el));
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(PartialAtomDsl(form.0))
}

fn constraint_tag(key: AtomConstraintKey) -> &'static str {
    match key {
        AtomConstraintKey::Valence => "#v",
        AtomConstraintKey::DonatedPairs => "#d",
        AtomConstraintKey::AcceptedPairs => "#t",
        AtomConstraintKey::AromaticValence => "#a",
        AtomConstraintKey::MulticenterValence => "#m",
        AtomConstraintKey::TetrahedralStereo => "#T",
        AtomConstraintKey::Degree => "#D",
        AtomConstraintKey::TotalDegree => "#X",
        AtomConstraintKey::TotalValence => "#V",
        AtomConstraintKey::RingDegree => "#x",
        AtomConstraintKey::RingValence => "#y",
        AtomConstraintKey::TotalHydrogens => "#H",
        AtomConstraintKey::RingMembership(_) => "#R",
    }
}

/// One predicate from an atom-string; the parser yields a `Vec` of these
/// and the applier folds them into the `AtomAst`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AtomPredicate {
    IsotopeMass(IsotopeMassAst),
    Charge(ValueAst),
    ImplicitHydrogens(ValueAst),
    LonePairs(ValueAst),
    Spin(SpinPredicate),
    Constraint(AtomConstraintAst),
}

fn atom_predicate(i: &mut &str) -> PResult<AtomPredicate> {
    let start = *i;
    if i.is_empty() {
        return Err(ErrMode::Backtrack(ParseError::Syntax));
    }
    if i.len() < 2 {
        return Err(ErrMode::Cut(ParseError::TrailingInput(start.to_string())));
    }
    let prefix: &str = take(2usize).parse_next(i)?;
    match prefix {
        "#i" => isotope.map(AtomPredicate::IsotopeMass).parse_next(i),
        "#c" => charge.map(AtomPredicate::Charge).parse_next(i),
        "#h" => optional_value
            .map(AtomPredicate::ImplicitHydrogens)
            .parse_next(i),
        "#n" => optional_value.map(AtomPredicate::LonePairs).parse_next(i),
        "#u" => optional_value
            .map(|v| AtomPredicate::Spin(SpinPredicate::Unpaired(v)))
            .parse_next(i),
        "#s" => optional_value
            .map(|v| AtomPredicate::Spin(SpinPredicate::Multiplicity(v)))
            .parse_next(i),
        "#v" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraintAst::Valence(v)))
            .parse_next(i),
        "#V" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraintAst::TotalValence(v)))
            .parse_next(i),
        "#d" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraintAst::DonatedPairs(v)))
            .parse_next(i),
        "#t" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraintAst::AcceptedPairs(v)))
            .parse_next(i),
        "#a" => aromatic_valence
            .map(|c| AtomPredicate::Constraint(AtomConstraintAst::AromaticValence(c)))
            .parse_next(i),
        "#m" => multicenter_valence
            .map(|c| AtomPredicate::Constraint(AtomConstraintAst::MulticenterValence(c)))
            .parse_next(i),
        "#D" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraintAst::Degree(v)))
            .parse_next(i),
        "#X" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraintAst::TotalDegree(v)))
            .parse_next(i),
        "#x" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraintAst::RingDegree(v)))
            .parse_next(i),
        "#y" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraintAst::RingValence(v)))
            .parse_next(i),
        "#H" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraintAst::TotalHydrogens(v)))
            .parse_next(i),
        "#R" => ring_membership
            .map(|m| AtomPredicate::Constraint(AtomConstraintAst::RingMembership(m)))
            .parse_next(i),
        "#T" => (|i: &mut &str| tetrahedral_stereo_config(i))
            .map(|c| AtomPredicate::Constraint(AtomConstraintAst::TetrahedralStereo(c)))
            .parse_next(i),
        p if p.starts_with('#') => Err(ErrMode::Cut(ParseError::UnknownAtomPredicate(
            p.to_string(),
        ))),
        _ => Err(ErrMode::Cut(ParseError::TrailingInput(start.to_string()))),
    }
}

fn element(i: &mut &str) -> PResult<ElementAst> {
    alt((
        '*'.value(ElementAst::Undetermined),
        preceded('!', element_set).map(ElementAst::not_set),
        preceded('!', element_literal).map(ElementAst::not),
        element_set.map(ElementAst::lit_set),
        element_bind.map(|(id, set, polarity)| match polarity {
            MemOp::In => ElementAst::var_in(id, set),
            MemOp::NotIn => ElementAst::var_not_in(id, set),
        }),
        element_ref.map(ElementAst::var),
        element_literal.map(ElementAst::Lit),
    ))
    .parse_next(i)
    .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedElement))
}

fn element_literal(i: &mut &str) -> PResult<Element> {
    let sym: &str = (
        one_of(|c: char| c.is_ascii_uppercase()),
        repeat::<_, _, (), _, _>(0.., one_of(|c: char| c.is_ascii_lowercase())),
    )
        .take()
        .parse_next(i)?;
    match Element::from_symbol(sym) {
        Some(el) => Ok(el),
        None => Err(ErrMode::Backtrack(ParseError::from_input(i))),
    }
}

fn element_set(i: &mut &str) -> PResult<Vec<Element>> {
    delimited(
        '{',
        delimited(
            multispace0,
            separated(
                1..,
                element_literal,
                delimited(multispace0, ',', multispace0),
            ),
            multispace0,
        ),
        '}',
    )
    .parse_next(i)
}

fn element_bind(i: &mut &str) -> PResult<(String, Vec<Element>, MemOp)> {
    alt((
        delimited('(', delimited(multispace0, element_bind, multispace0), ')'),
        (
            preceded('?', id),
            preceded(
                delimited(multispace0, "::", multispace0),
                element_bind_domain,
            ),
        )
            .map(|(id, (set, polarity))| (id, set, polarity)),
    ))
    .parse_next(i)
}

fn element_bind_domain(i: &mut &str) -> PResult<(Vec<Element>, MemOp)> {
    alt((
        preceded('!', element_set).map(|s| (s, MemOp::NotIn)),
        preceded('!', element_literal).map(|e| (vec![e], MemOp::NotIn)),
        element_set.map(|s| (s, MemOp::In)),
    ))
    .parse_next(i)
}

fn element_ref(i: &mut &str) -> PResult<String> {
    alt((
        delimited('(', delimited(multispace0, element_ref, multispace0), ')'),
        preceded('?', id),
    ))
    .parse_next(i)
}

fn isotope(i: &mut &str) -> PResult<IsotopeMassAst> {
    preceded(
        multispace0,
        alt((
            '='.value(IsotopeMassAst::Natural),
            '*'.value(IsotopeMassAst::Undetermined),
            isotope_set.map(IsotopeMassAst::lit_set),
            terminated(isotope_var, (multispace0, terminator)).map(|(id, domain)| match domain {
                Some(set) => IsotopeMassAst::var_in(id, set),
                None => IsotopeMassAst::var(id),
            }),
            terminated(dec_uint::<_, u32, _>, (multispace0, terminator)).map(IsotopeMassAst::Lit),
        )),
    )
    .parse_next(i)
    .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedPredicateBody))
}

fn isotope_set(i: &mut &str) -> PResult<Vec<u32>> {
    delimited(
        '{',
        delimited(
            multispace0,
            separated(
                1..,
                dec_uint::<_, u32, _>,
                delimited(multispace0, ',', multispace0),
            ),
            multispace0,
        ),
        '}',
    )
    .parse_next(i)
}

fn isotope_var(i: &mut &str) -> PResult<(String, Option<Vec<u32>>)> {
    alt((
        delimited('(', delimited(multispace0, isotope_var, multispace0), ')'),
        (
            preceded('?', id),
            opt(preceded(
                delimited(multispace0, "::", multispace0),
                isotope_set,
            )),
        ),
    ))
    .parse_next(i)
}

fn aromatic_valence(i: &mut &str) -> PResult<AromaticValenceAst> {
    preceded(
        multispace0,
        alt((
            "*".value(AromaticValenceAst::Undetermined),
            "!".value(AromaticValenceAst::NotAromatic),
            // `#a+` encodes "aromatic, count unspecified" — structurally
            // distinct from the outer Undetermined; canonical form is
            // Aromatic(Undetermined).
            "+".value(AromaticValenceAst::Aromatic(ValueAst::Undetermined)),
            value.map(AromaticValenceAst::Aromatic),
            empty.value(AromaticValenceAst::Aromatic(ValueAst::Lit(1))),
        )),
    )
    .parse_next(i)
    .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedPredicateBody))
}

fn multicenter_valence(i: &mut &str) -> PResult<MulticenterValenceAst> {
    preceded(
        multispace0,
        alt((
            "*".value(MulticenterValenceAst::Undetermined),
            "!".value(MulticenterValenceAst::NotMulticenter),
            // `#m+` mirrors `#a+` — "multicenter, count unspecified".
            "+".value(MulticenterValenceAst::Multicenter(ValueAst::Undetermined)),
            value.map(MulticenterValenceAst::Multicenter),
            empty.value(MulticenterValenceAst::Multicenter(ValueAst::Lit(1))),
        )),
    )
    .parse_next(i)
    .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedPredicateBody))
}

fn apply_predicates(form: &mut AtomDsl, preds: Vec<AtomPredicate>) -> Result<(), ParseError> {
    let ast = &mut form.0;
    for pred in preds {
        match pred {
            AtomPredicate::IsotopeMass(v) => {
                if !matches!(ast.isotope_mass, IsotopeMassAst::Undetermined) {
                    return Err(ParseError::DuplicateAtomPredicate("#i".to_string()));
                }
                ast.isotope_mass = v;
            }
            AtomPredicate::Charge(v) => {
                if !matches!(ast.charge, ValueAst::Undetermined) {
                    return Err(ParseError::DuplicateAtomPredicate("#c".to_string()));
                }
                ast.charge = v;
            }
            AtomPredicate::ImplicitHydrogens(v) => {
                if !matches!(ast.implicit_hydrogens, ValueAst::Undetermined) {
                    return Err(ParseError::DuplicateAtomPredicate("#h".to_string()));
                }
                ast.implicit_hydrogens = v;
            }
            AtomPredicate::LonePairs(v) => {
                if !matches!(ast.lone_pairs, ValueAst::Undetermined) {
                    return Err(ParseError::DuplicateAtomPredicate("#n".to_string()));
                }
                ast.lone_pairs = v;
            }
            AtomPredicate::Spin(sp) => {
                apply_spin_pair(&mut ast.spin, sp, ParseError::DuplicateAtomPredicate)?;
            }
            AtomPredicate::Constraint(c) => {
                if ast.constraints.contains(c.key()) {
                    return Err(ParseError::DuplicateAtomPredicate(
                        constraint_tag(c.key()).to_string(),
                    ));
                }
                ast.constraints.set(c);
            }
        }
    }
    Ok(())
}

fn fmt_atom_ast(f: &mut fmt::Formatter<'_>, ast: &AtomAst) -> fmt::Result {
    fmt_element(f, &ast.element)?;
    fmt_isotope_mass(f, &ast.isotope_mass)?;
    fmt_charge(f, &ast.charge)?;
    fmt_value_field(f, "#h", &ast.implicit_hydrogens)?;
    fmt_value_field(f, "#n", &ast.lone_pairs)?;
    fmt_spin_pair(f, &ast.spin)
}

fn fmt_element(f: &mut fmt::Formatter<'_>, expr: &ElementAst) -> fmt::Result {
    match expr {
        ElementAst::Lit(e) => write!(f, "{}", e),
        ElementAst::Undetermined => write!(f, "*"),
        ElementAst::LitSet(es) => fmt_element_set(f, es.iter().copied()),
        ElementAst::NotSet(es) if es.len() == 1 => {
            // A singleton complement renders as `!e` (no braces), e.g. `!H`.
            write!(f, "!{}", es.iter().next().unwrap())
        }
        ElementAst::NotSet(es) => {
            write!(f, "!")?;
            fmt_element_set(f, es.iter().copied())
        }
        ElementAst::Var(v) => {
            let (id, domain) = &**v;
            write!(f, "?{}", id)?;
            if let Some((op, set)) = domain {
                write!(f, " :: ")?;
                if matches!(op, MemOp::NotIn) {
                    write!(f, "!")?;
                }
                fmt_element_set(f, set.iter().copied())?;
            }
            Ok(())
        }
    }
}

fn fmt_element_set(
    f: &mut fmt::Formatter<'_>,
    es: impl IntoIterator<Item = Element>,
) -> fmt::Result {
    write!(f, "{{")?;
    for (i, e) in es.into_iter().enumerate() {
        if i > 0 {
            write!(f, ",")?;
        }
        write!(f, "{}", e)?;
    }
    write!(f, "}}")
}

fn fmt_isotope_mass(f: &mut fmt::Formatter<'_>, iso: &IsotopeMassAst) -> fmt::Result {
    match iso {
        IsotopeMassAst::Undetermined => Ok(()),
        IsotopeMassAst::Natural => write!(f, "#i="),
        IsotopeMassAst::Lit(n) => write!(f, "#i{}", n),
        IsotopeMassAst::LitSet(s) => {
            write!(f, "#i")?;
            fmt_set(f, s.iter().copied())
        }
        IsotopeMassAst::Var(v) => {
            let (id, domain) = &**v;
            write!(f, "#i?{}", id)?;
            match domain {
                Some(set) => {
                    write!(f, " :: ")?;
                    fmt_set(f, set.iter().copied())
                }
                None => Ok(()),
            }
        }
    }
}

/// Format a value field with `Lit(1)` sugared as the bare prefix. Only
/// `Undetermined` elides; every literal (including `Lit(0)`) must render so
/// parsing recovers it.
fn fmt_value_field(f: &mut fmt::Formatter<'_>, prefix: &str, v: &ValueAst) -> fmt::Result {
    match v {
        ValueAst::Undetermined => Ok(()),
        ValueAst::Lit(1) => write!(f, "{}", prefix),
        ValueAst::Lit(n) => write!(f, "{}{}", prefix, n),
        v => {
            write!(f, "{}", prefix)?;
            fmt_value(f, v)
        }
    }
}

/// Format an inline-constraint value field. Per the canonical-rendering
/// rules in `dsl::predicates`, vacuous constraints (`Undetermined`) elide.
/// `Lit(0)` is a meaningful constraint and renders.
fn fmt_value_field_required(f: &mut fmt::Formatter<'_>, prefix: &str, v: &ValueAst) -> fmt::Result {
    match v {
        ValueAst::Undetermined => Ok(()),
        ValueAst::Lit(1) => write!(f, "{}", prefix),
        ValueAst::Lit(n) => write!(f, "{}{}", prefix, n),
        v => {
            write!(f, "{}", prefix)?;
            fmt_value(f, v)
        }
    }
}

fn fmt_constraint(f: &mut fmt::Formatter<'_>, c: &AtomConstraintAst) -> fmt::Result {
    match c {
        AtomConstraintAst::Valence(v) => fmt_value_field_required(f, "#v", v),
        AtomConstraintAst::DonatedPairs(v) => fmt_value_field_required(f, "#d", v),
        AtomConstraintAst::AcceptedPairs(v) => fmt_value_field_required(f, "#t", v),
        AtomConstraintAst::MulticenterValence(c) => match c {
            MulticenterValenceAst::Undetermined => Ok(()),
            MulticenterValenceAst::NotMulticenter => write!(f, "#m!"),
            MulticenterValenceAst::Multicenter(ValueAst::Undetermined) => write!(f, "#m+"),
            MulticenterValenceAst::Multicenter(ValueAst::Lit(1)) => write!(f, "#m"),
            MulticenterValenceAst::Multicenter(ValueAst::Lit(n)) => write!(f, "#m{}", n),
            MulticenterValenceAst::Multicenter(v) => {
                write!(f, "#m")?;
                fmt_value(f, v)
            }
        },
        AtomConstraintAst::AromaticValence(c) => match c {
            AromaticValenceAst::Undetermined => Ok(()),
            AromaticValenceAst::NotAromatic => write!(f, "#a!"),
            AromaticValenceAst::Aromatic(ValueAst::Undetermined) => write!(f, "#a+"),
            AromaticValenceAst::Aromatic(ValueAst::Lit(1)) => write!(f, "#a"),
            AromaticValenceAst::Aromatic(ValueAst::Lit(n)) => write!(f, "#a{}", n),
            AromaticValenceAst::Aromatic(v) => {
                write!(f, "#a")?;
                fmt_value(f, v)
            }
        },
        AtomConstraintAst::Degree(v) => fmt_value_field_required(f, "#D", v),
        AtomConstraintAst::TotalDegree(v) => fmt_value_field_required(f, "#X", v),
        AtomConstraintAst::RingDegree(v) => fmt_value_field_required(f, "#x", v),
        AtomConstraintAst::RingValence(v) => fmt_value_field_required(f, "#y", v),
        AtomConstraintAst::TotalValence(v) => fmt_value_field_required(f, "#V", v),
        AtomConstraintAst::TotalHydrogens(v) => fmt_value_field_required(f, "#H", v),
        AtomConstraintAst::RingMembership(m) => fmt_ring_membership(f, m),
        AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::Undetermined) => Ok(()),
        AtomConstraintAst::TetrahedralStereo(c) => {
            write!(f, "#T")?;
            fmt_tetrahedral_stereo_config(f, c)
        }
    }
}

pub(crate) fn raise_atom(ast: &mut AtomAst, cfg: &AtomDefaults) {
    // Exhaustive destructure: adding a new AtomAst field is a compile error
    // here, forcing the author to decide how raising should handle it.
    let AtomAst {
        element: _,
        isotope_mass,
        charge,
        implicit_hydrogens,
        lone_pairs,
        spin,
        constraints,
    } = ast;

    if matches!(*isotope_mass, IsotopeMassAst::Undetermined) {
        *isotope_mass = match cfg.isotope {
            IsotopeDefault::Natural => IsotopeMassAst::Natural,
            IsotopeDefault::Required => IsotopeMassAst::Undetermined,
        };
    }
    if matches!(*charge, ValueAst::Undetermined) {
        *charge = match cfg.charge {
            NumericDefault::Zero => ValueAst::Lit(0),
            NumericDefault::Required => ValueAst::Undetermined,
        };
    }
    if matches!(*implicit_hydrogens, ValueAst::Undetermined) {
        *implicit_hydrogens = match cfg.implicit_hydrogens {
            NumericDefault::Zero => ValueAst::Lit(0),
            NumericDefault::Required => ValueAst::Undetermined,
        };
    }
    if matches!(*lone_pairs, ValueAst::Undetermined) {
        *lone_pairs = match cfg.lone_pairs {
            NumericDefault::Zero => ValueAst::Lit(0),
            NumericDefault::Required => ValueAst::Undetermined,
        };
    }
    raise_spin(spin, cfg.unpaired_electrons, cfg.multiplicity);
    raise_atom_constraints(constraints, cfg);
}

/// A defaulted key wants filling iff it is absent or holds a vacuous (`Undetermined`) value.
/// A concrete user value is left alone.
fn is_unset_or_vacuous(constraints: &AtomConstraintsAst, key: AtomConstraintKey) -> bool {
    constraints.get(key).is_none_or(|c| c.is_undetermined())
}

fn raise_atom_constraints(constraints: &mut AtomConstraintsAst, cfg: &AtomDefaults) {
    // One explicit clause per defaulted kind, in ascending key-sort order. No global vacuous
    // strip: a defaulted kind fills its own absent/vacuous entry; vacuous entries of other kinds
    // are left for lazy canonicalization.
    if matches!(cfg.valence, NumericDefault::Zero)
        && is_unset_or_vacuous(constraints, AtomConstraintKey::Valence)
    {
        constraints.set(AtomConstraintAst::Valence(ValueAst::Lit(0)));
    }
    if matches!(cfg.aromatic_valence, AromaticValenceDefault::NotAromatic)
        && is_unset_or_vacuous(constraints, AtomConstraintKey::AromaticValence)
    {
        constraints.set(AtomConstraintAst::AromaticValence(
            AromaticValenceAst::NotAromatic,
        ));
    }
    if matches!(
        cfg.multicenter_valence,
        MulticenterValenceDefault::NotMulticenter
    ) && is_unset_or_vacuous(constraints, AtomConstraintKey::MulticenterValence)
    {
        constraints.set(AtomConstraintAst::MulticenterValence(
            MulticenterValenceAst::NotMulticenter,
        ));
    }
    if matches!(cfg.donated_pairs, NumericDefault::Zero)
        && is_unset_or_vacuous(constraints, AtomConstraintKey::DonatedPairs)
    {
        constraints.set(AtomConstraintAst::DonatedPairs(ValueAst::Lit(0)));
    }
    if matches!(cfg.accepted_pairs, NumericDefault::Zero)
        && is_unset_or_vacuous(constraints, AtomConstraintKey::AcceptedPairs)
    {
        constraints.set(AtomConstraintAst::AcceptedPairs(ValueAst::Lit(0)));
    }
    if matches!(cfg.tetrahedral_stereo, StereoDefault::NotStereo)
        && is_unset_or_vacuous(constraints, AtomConstraintKey::TetrahedralStereo)
    {
        constraints.set(AtomConstraintAst::TetrahedralStereo(
            TetrahedralStereoAst::NotStereo,
        ));
    }
}

pub(crate) fn lower_atom(ast: &mut AtomAst, cfg: &AtomDefaults) {
    // Exhaustive destructure: adding a new AtomAst field is a compile error
    // here, forcing the author to decide how lowering should handle it.
    let AtomAst {
        element: _,
        isotope_mass,
        charge,
        implicit_hydrogens,
        lone_pairs,
        spin,
        constraints,
    } = ast;

    if matches!(
        (&cfg.isotope, &*isotope_mass),
        (IsotopeDefault::Natural, IsotopeMassAst::Natural)
    ) {
        *isotope_mass = IsotopeMassAst::Undetermined;
    }
    if matches!(
        (&cfg.charge, &*charge),
        (NumericDefault::Zero, ValueAst::Lit(0))
    ) {
        *charge = ValueAst::Undetermined;
    }
    match (&cfg.implicit_hydrogens, &*implicit_hydrogens) {
        (NumericDefault::Required, ValueAst::Undetermined) => {
            *implicit_hydrogens = ValueAst::Undetermined;
        }
        (NumericDefault::Zero, ValueAst::Lit(0)) => {
            *implicit_hydrogens = ValueAst::Undetermined;
        }
        _ => {}
    }
    if matches!(
        (&cfg.lone_pairs, &*lone_pairs),
        (NumericDefault::Zero, ValueAst::Lit(0))
    ) {
        *lone_pairs = ValueAst::Undetermined;
    }
    lower_spin(spin, cfg.unpaired_electrons, cfg.multiplicity);
    lower_atom_constraints(constraints, cfg);
}

/// Elide a defaulted entry: if its key holds exactly `default`, remove it (raise would refill it
/// from the same `cfg`). No-op otherwise.
fn remove_if_default(constraints: &mut AtomConstraintsAst, default: AtomConstraintAst) {
    let key = default.key();
    if constraints.get(key) == Some(&default) {
        constraints.remove(key);
    }
}

fn lower_atom_constraints(constraints: &mut AtomConstraintsAst, cfg: &AtomDefaults) {
    // Elide each defaulted entry equal to its default, in reverse key-sort order (mirror of raise).
    if matches!(cfg.tetrahedral_stereo, StereoDefault::NotStereo) {
        remove_if_default(
            constraints,
            AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::NotStereo),
        );
    }
    if matches!(cfg.accepted_pairs, NumericDefault::Zero) {
        remove_if_default(
            constraints,
            AtomConstraintAst::AcceptedPairs(ValueAst::Lit(0)),
        );
    }
    if matches!(cfg.donated_pairs, NumericDefault::Zero) {
        remove_if_default(
            constraints,
            AtomConstraintAst::DonatedPairs(ValueAst::Lit(0)),
        );
    }
    if matches!(
        cfg.multicenter_valence,
        MulticenterValenceDefault::NotMulticenter
    ) {
        remove_if_default(
            constraints,
            AtomConstraintAst::MulticenterValence(MulticenterValenceAst::NotMulticenter),
        );
    }
    if matches!(cfg.aromatic_valence, AromaticValenceDefault::NotAromatic) {
        remove_if_default(
            constraints,
            AtomConstraintAst::AromaticValence(AromaticValenceAst::NotAromatic),
        );
    }
    if matches!(cfg.valence, NumericDefault::Zero) {
        remove_if_default(constraints, AtomConstraintAst::Valence(ValueAst::Lit(0)));
    }
}

/// Surface DSL wrapper around `AromaticValenceAst`. EDN form: `:undetermined`,
/// `:not-aromatic`, or `{:aromatic <value>}`.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct AromaticValenceDsl(pub AromaticValenceAst);

impl<'de> FromEdn<'de> for AromaticValenceDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Keyword(k) if k.name() == "undetermined" => {
                Ok(Self(AromaticValenceAst::Undetermined))
            }
            Edn::Keyword(k) if k.name() == "not-aromatic" => {
                Ok(Self(AromaticValenceAst::NotAromatic))
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
                    "aromatic" => Ok(Self(AromaticValenceAst::Aromatic(
                        ValueDsl::from_edn(v)?.into_ast(&()),
                    ))),
                    other => Err(DeError::UnknownField {
                        key: other.to_string(),
                        path: vec!["aromatic-valence".into()],
                    }),
                }
            }
            other => Err(DeError::TypeMismatch {
                expected: ":undetermined / :not-aromatic / {:aromatic <value>}",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for AromaticValenceDsl {
    fn to_edn(&self) -> Edn<'static> {
        match &self.0 {
            AromaticValenceAst::Undetermined => {
                Edn::Keyword(EdnKeyword::owned("undetermined".into()))
            }
            AromaticValenceAst::NotAromatic => {
                Edn::Keyword(EdnKeyword::owned("not-aromatic".into()))
            }
            AromaticValenceAst::Aromatic(v) => {
                single_key_map("aromatic", ValueDsl::from_ast(v, &()).to_edn())
            }
        }
    }
}

impl FromAst<AromaticValenceAst> for AromaticValenceDsl {
    type Ctx = ();

    fn from_ast(ast: &AromaticValenceAst, _ctx: &Self::Ctx) -> Self {
        Self(ast.clone())
    }
}

impl IntoAst<AromaticValenceAst> for AromaticValenceDsl {
    type Ctx = ();

    fn into_ast(self, _ctx: &Self::Ctx) -> AromaticValenceAst {
        self.0
    }
}

/// Surface DSL wrapper around `MulticenterValenceAst`. EDN form:
/// `:undetermined`, `:not-multicenter`, or `{:multicenter <value>}`.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct MulticenterValenceDsl(pub MulticenterValenceAst);

impl<'de> FromEdn<'de> for MulticenterValenceDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Keyword(k) if k.name() == "undetermined" => {
                Ok(Self(MulticenterValenceAst::Undetermined))
            }
            Edn::Keyword(k) if k.name() == "not-multicenter" => {
                Ok(Self(MulticenterValenceAst::NotMulticenter))
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
                    "multicenter" => Ok(Self(MulticenterValenceAst::Multicenter(
                        ValueDsl::from_edn(v)?.into_ast(&()),
                    ))),
                    other => Err(DeError::UnknownField {
                        key: other.to_string(),
                        path: vec!["multicenter-valence".into()],
                    }),
                }
            }
            other => Err(DeError::TypeMismatch {
                expected: ":undetermined / :not-multicenter / {:multicenter <value>}",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for MulticenterValenceDsl {
    fn to_edn(&self) -> Edn<'static> {
        match &self.0 {
            MulticenterValenceAst::Undetermined => {
                Edn::Keyword(EdnKeyword::owned("undetermined".into()))
            }
            MulticenterValenceAst::NotMulticenter => {
                Edn::Keyword(EdnKeyword::owned("not-multicenter".into()))
            }
            MulticenterValenceAst::Multicenter(v) => {
                single_key_map("multicenter", ValueDsl::from_ast(v, &()).to_edn())
            }
        }
    }
}

impl FromAst<MulticenterValenceAst> for MulticenterValenceDsl {
    type Ctx = ();

    fn from_ast(ast: &MulticenterValenceAst, _ctx: &Self::Ctx) -> Self {
        Self(ast.clone())
    }
}

impl IntoAst<MulticenterValenceAst> for MulticenterValenceDsl {
    type Ctx = ();

    fn into_ast(self, _ctx: &Self::Ctx) -> MulticenterValenceAst {
        self.0
    }
}

/// Surface DSL wrapper around `AtomConstraintAst`. EDN form is a single-key map
/// keyed by the constraint kind: e.g. `{:valence 4}`, `{:degree *}`,
/// `{:aromatic-valence :not-aromatic}`, `{:total-hydrogens {?h :: {0,1,2}}}`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct AtomConstraintDsl(pub AtomConstraintAst);

impl<'de> FromEdn<'de> for AtomConstraintDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        let Edn::Map(m) = edn else {
            return Err(DeError::TypeMismatch {
                expected: "atom-constraint single-key map",
                got: edn.kind(),
                path: Vec::new(),
            });
        };
        if m.len() != 1 {
            return Err(DeError::Custom(format!(
                "atom-constraint must have exactly one key, got {}",
                m.len()
            )));
        }
        let (k, v) = m.iter().next().unwrap();
        let Edn::Keyword(key) = k else {
            return Err(DeError::TypeMismatch {
                expected: "keyword key",
                got: k.kind(),
                path: vec!["atom-constraint".into()],
            });
        };
        let c = match key.name() {
            "valence" => AtomConstraintAst::Valence(ValueDsl::from_edn(v)?.into_ast(&())),
            "total-valence" => {
                AtomConstraintAst::TotalValence(ValueDsl::from_edn(v)?.into_ast(&()))
            }
            "aromatic-valence" => {
                AtomConstraintAst::AromaticValence(AromaticValenceDsl::from_edn(v)?.into_ast(&()))
            }
            "multicenter-valence" => AtomConstraintAst::MulticenterValence(
                MulticenterValenceDsl::from_edn(v)?.into_ast(&()),
            ),
            "donated-pairs" => {
                AtomConstraintAst::DonatedPairs(ValueDsl::from_edn(v)?.into_ast(&()))
            }
            "accepted-pairs" => {
                AtomConstraintAst::AcceptedPairs(ValueDsl::from_edn(v)?.into_ast(&()))
            }
            "degree" => AtomConstraintAst::Degree(ValueDsl::from_edn(v)?.into_ast(&())),
            "total-degree" => AtomConstraintAst::TotalDegree(ValueDsl::from_edn(v)?.into_ast(&())),
            "ring-degree" => AtomConstraintAst::RingDegree(ValueDsl::from_edn(v)?.into_ast(&())),
            "ring-valence" => AtomConstraintAst::RingValence(ValueDsl::from_edn(v)?.into_ast(&())),
            "total-hydrogens" => {
                AtomConstraintAst::TotalHydrogens(ValueDsl::from_edn(v)?.into_ast(&()))
            }
            "ring-membership" => {
                AtomConstraintAst::RingMembership(RingMembershipDsl::from_edn(v)?.0)
            }
            "tetrahedral-stereo" => AtomConstraintAst::TetrahedralStereo(
                TetrahedralStereoDsl::from_edn(v)?.into_ast(&()),
            ),
            other => {
                return Err(DeError::UnknownField {
                    key: other.to_string(),
                    path: vec!["atom-constraint".into()],
                });
            }
        };
        Ok(Self(c))
    }
}

impl ToEdn for AtomConstraintDsl {
    fn to_edn(&self) -> Edn<'static> {
        match &self.0 {
            AtomConstraintAst::Valence(v) => {
                single_key_map("valence", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraintAst::TotalValence(v) => {
                single_key_map("total-valence", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraintAst::AromaticValence(c) => single_key_map(
                "aromatic-valence",
                AromaticValenceDsl::from_ast(c, &()).to_edn(),
            ),
            AtomConstraintAst::MulticenterValence(c) => single_key_map(
                "multicenter-valence",
                MulticenterValenceDsl::from_ast(c, &()).to_edn(),
            ),
            AtomConstraintAst::DonatedPairs(v) => {
                single_key_map("donated-pairs", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraintAst::AcceptedPairs(v) => {
                single_key_map("accepted-pairs", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraintAst::Degree(v) => {
                single_key_map("degree", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraintAst::TotalDegree(v) => {
                single_key_map("total-degree", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraintAst::RingDegree(v) => {
                single_key_map("ring-degree", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraintAst::RingValence(v) => {
                single_key_map("ring-valence", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraintAst::TotalHydrogens(v) => {
                single_key_map("total-hydrogens", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraintAst::RingMembership(m) => {
                single_key_map("ring-membership", RingMembershipDsl(m.clone()).to_edn())
            }
            AtomConstraintAst::TetrahedralStereo(c) => single_key_map(
                "tetrahedral-stereo",
                TetrahedralStereoDsl::from_ast(c, &()).to_edn(),
            ),
        }
    }
}

impl FromAst<AtomConstraintAst> for AtomConstraintDsl {
    type Ctx = ();

    fn from_ast(ast: &AtomConstraintAst, _ctx: &Self::Ctx) -> Self {
        Self(ast.clone())
    }
}

impl IntoAst<AtomConstraintAst> for AtomConstraintDsl {
    type Ctx = ();

    fn into_ast(self, _ctx: &Self::Ctx) -> AtomConstraintAst {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_chem::element::Element;
    use umol_edn::read_string;

    use super::*;
    use crate::ast::constraint::RingScope;
    use crate::ast::operators::{MemOp, RelOp};
    use crate::ast::spin::SpinStateAst;
    use crate::ast::stereo::{StereoCosetAst, StereoTerm};
    use crate::ast::value::{ValuePredicate, ValueTerm};

    #[rstest]
    #[case::single("C", AtomDsl(AtomAst::new(ElementAst::Lit(Element::C))))]
    #[case::double("Cl", AtomDsl(AtomAst::new(ElementAst::Lit(Element::Cl))))]
    #[case::transuranic("Og", AtomDsl(AtomAst::new(ElementAst::Lit(Element::Og))))]
    fn test_parse_bare_element(#[case] input: &str, #[case] expected: AtomDsl) {
        assert_eq!(parse_bare_element(input), Some(expected));
    }

    #[rstest]
    #[case::empty("")]
    #[case::lowercase_first("cl")]
    #[case::unknown_symbol("Zx")]
    #[case::three_bytes("Abc")]
    #[case::leading_whitespace(" C")]
    #[case::trailing_whitespace("C ")]
    #[case::wildcard("*")]
    #[case::set("{C,N}")]
    #[case::var("(?e)")]
    #[case::h_count("C#h3")]
    #[case::charge_plus("N#c+")]
    #[case::full("C#c+1#R+#v4")]
    fn test_parse_bare_element_error(#[case] input: &str) {
        assert_eq!(parse_bare_element(input), None);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::carbon("C", AtomDsl(AtomAst::new(ElementAst::Lit(Element::C))))]
    #[case::iron("Fe", AtomDsl(AtomAst::new(ElementAst::Lit(Element::Fe))))]
    #[case::chlorine("Cl", AtomDsl(AtomAst::new(ElementAst::Lit(Element::Cl))))]
    #[case::whitespace("  C  ", AtomDsl(AtomAst::new(ElementAst::Lit(Element::C))))]
    #[case::undetermined("*", AtomDsl(AtomAst::new(ElementAst::Undetermined)))]
    #[case::element_set("{C,N,O}", AtomDsl(AtomAst::new(ElementAst::lit_set(vec![Element::C, Element::N, Element::O]))))]
    #[case::element_set_singleton("{C}", AtomDsl(AtomAst::new(ElementAst::lit_set(vec![Element::C]))))]
    #[case::element_var_mem_in("(?e :: {C,N})", AtomDsl(AtomAst::new(ElementAst::var_in("e", vec![Element::C, Element::N]))))]
    #[case::element_var("?e", AtomDsl(AtomAst::new(ElementAst::var("e".to_string()))))]
    #[case::isotope("C#i12", AtomDsl(AtomAst { isotope_mass: IsotopeMassAst::Lit(12), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::isotope_natural("C#i=", AtomDsl(AtomAst { isotope_mass: IsotopeMassAst::Natural, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_pos("C#c+2", AtomDsl(AtomAst { charge: ValueAst::Lit(2), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_neg("C#c-2", AtomDsl(AtomAst { charge: ValueAst::Lit(-2), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_plus("C#c+", AtomDsl(AtomAst { charge: ValueAst::Lit(1), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_minus("C#c-", AtomDsl(AtomAst { charge: ValueAst::Lit(-1), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_zero("C#c0", AtomDsl(AtomAst { charge: ValueAst::Lit(0), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_count("C#h3", AtomDsl(AtomAst { implicit_hydrogens: ValueAst::Lit(3), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_undetermined("C#h*", AtomDsl(AtomAst { implicit_hydrogens: ValueAst::Undetermined, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_bind("C#h(?h)", AtomDsl(AtomAst { implicit_hydrogens: ValueAst::var("h"), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_set("N#h?h :: {2,3}", AtomDsl(AtomAst { implicit_hydrogens: ValueAst::predicate(ValuePredicate::Mem(ValueTerm::Var("h".to_string()), MemOp::In, BTreeSet::from([2, 3]))), ..AtomAst::new(ElementAst::Lit(Element::N)) }))]
    #[case::h_expr("C#h?h >= 1", AtomDsl(AtomAst { implicit_hydrogens: ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Ge, ValueTerm::Lit(1))), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_omit("C#h", AtomDsl(AtomAst { implicit_hydrogens: ValueAst::Lit(1), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::lone_pairs("O#n2", AtomDsl(AtomAst { lone_pairs: ValueAst::Lit(2), ..AtomAst::new(ElementAst::Lit(Element::O)) }))]
    #[case::lone_pairs_omit("O#n", AtomDsl(AtomAst { lone_pairs: ValueAst::Lit(1), ..AtomAst::new(ElementAst::Lit(Element::O)) }))]
    #[case::unpaired("C#u2", AtomDsl(AtomAst { spin: SpinStateAst { unpaired: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined }, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::unpaired_no_canonicalization("C#u{1}", AtomDsl(AtomAst { spin: SpinStateAst { unpaired: ValueAst::lit_set([1]), multiplicity: ValueAst::Undetermined }, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::unpaired_omit("C#u", AtomDsl(AtomAst { spin: SpinStateAst { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Undetermined }, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multiplicity("C#s3", AtomDsl(AtomAst { spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(3) }, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multiplicity_omit("C#s", AtomDsl(AtomAst { spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(1) }, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::valence("C#v4", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::Valence(ValueAst::Lit(4))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::donated_pairs("N#d1", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::DonatedPairs(ValueAst::Lit(1))]), ..AtomAst::new(ElementAst::Lit(Element::N)) }))]
    #[case::accepted_pairs("B#t1", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::AcceptedPairs(ValueAst::Lit(1))]), ..AtomAst::new(ElementAst::Lit(Element::B)) }))]
    #[case::ring_membership_size("C#R(6)", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::ring_membership(RingScope::Size(6), 1)]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::arom_not_aromatic("C#a!", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::AromaticValence(AromaticValenceAst::NotAromatic)]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::arom_undetermined("C#a*", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::AromaticValence(AromaticValenceAst::Undetermined)]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::arom_plus("C#a+", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Undetermined))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::arom_zero("C#a0", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(0)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::arom_one("C#a1", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(1)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::arom_omit("C#a", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(1)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multicenter_not("C#m!", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::MulticenterValence(MulticenterValenceAst::NotMulticenter)]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multicenter_undetermined("C#m*", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::MulticenterValence(MulticenterValenceAst::Undetermined)]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multicenter_plus("C#m+", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Undetermined))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multicenter_zero("C#m0", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(0)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multicenter_one("C#m", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(1)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multicenter("C#m2", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(2)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::degree("C#D2", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::Degree(ValueAst::Lit(2))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::total_degree("C#X3", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::TotalDegree(ValueAst::Lit(3))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_degree("C#x2", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::RingDegree(ValueAst::Lit(2))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_valence("C#y3", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::RingValence(ValueAst::Lit(3))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::total_valence("C#V5", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::TotalValence(ValueAst::Lit(5))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::total_hydrogens("C#H1", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::TotalHydrogens(ValueAst::Lit(1))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_membership_all_bare("C#R", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::ring_membership(RingScope::All, ValueAst::Lit(1))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_membership_all_star("C#R*", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::ring_membership(RingScope::All, ValueAst::Undetermined)]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_membership_all_plus("C#R+", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::ring_membership(RingScope::All, ValueAst::RangeFrom(1))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_not_in_ring("C#R!", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::ring_membership(RingScope::All, ValueAst::Lit(0))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_zero("C#R0", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::ring_membership(RingScope::All, ValueAst::Lit(0))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_membership_all("C#R2", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::ring_membership(RingScope::All, ValueAst::Lit(2))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_membership_size_conj("C#R(5)#R(6)", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::ring_membership(RingScope::Size(5), 1), AtomConstraintAst::ring_membership(RingScope::Size(6), 1)]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::tetrahedral_stereo_stereo("C#T+", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::Stereo(StereoCosetAst::Undetermined))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::tetrahedral_stereo_lit("C#T1", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(1)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::tetrahedral_stereo_not_stereo("C#T!", AtomDsl(AtomAst { constraints: AtomConstraintsAst::from_iter([AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::NotStereo)]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::whitespace_before_predicate("C #h3", AtomDsl(AtomAst { implicit_hydrogens: ValueAst::Lit(3), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::whitespace_between_predicates("C#c+ #h3", AtomDsl(AtomAst { charge: ValueAst::Lit(1), implicit_hydrogens: ValueAst::Lit(3), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::whitespace_surrounding_predicates("  C  #c+  #h3  ", AtomDsl(AtomAst { charge: ValueAst::Lit(1), implicit_hydrogens: ValueAst::Lit(3), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    fn test_parse_atom(#[case] input: &str, #[case] expected: AtomDsl) {
        let result = atom.parse(input);
        assert!(result.is_ok(), "{:?} should succeed, got {:?}", input, result.clone().unwrap_err());
        let form = result.unwrap();
        assert_eq!(form, expected);
    }

    #[rstest]
    #[case::empty("", ParseError::ExpectedElement)]
    #[case::no_element("#h3", ParseError::ExpectedElement)]
    #[case::unknown_pred("C#z", ParseError::UnknownAtomPredicate("#z".to_string()))]
    #[case::whitespace_after_hash("C# h3", ParseError::UnknownAtomPredicate("# ".to_string()))]
    #[case::dup_implicit_hydrogens("C#h3#h2", ParseError::DuplicateAtomPredicate("#h".to_string()))]
    #[case::dup_charge("C#c+#c-", ParseError::DuplicateAtomPredicate("#c".to_string()))]
    #[case::dup_valence("C#v3#v4", ParseError::DuplicateAtomPredicate("#v".to_string()))]
    #[case::dup_ring_same_scope("C#R(6)#R(6)", ParseError::DuplicateAtomPredicate("#R".to_string()))]
    #[case::invalid_special_slash("C#h/", ParseError::TrailingInput("/".to_string()))]
    #[case::invalid_special_minus("C#h-", ParseError::TrailingInput("-".to_string()))]
    #[case::invalid_special_equal("C#h=", ParseError::TrailingInput("=".to_string()))]
    #[case::trailing("C#h3 foo", ParseError::TrailingInput("foo".to_string()))]
    fn test_parse_atom_error(#[case] input: &str, #[case] expected: ParseError) {
        let result = atom.parse(input);
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

    #[rstest]
    #[case::empty("", PartialAtomDsl(AtomAst::new(ElementAst::Undetermined)))]
    #[case::element_only("O", PartialAtomDsl(AtomAst::new(ElementAst::Lit(Element::O))))]
    #[case::charge_only("#c-1", PartialAtomDsl(AtomAst::new(ElementAst::Undetermined).with_charge(ValueAst::Lit(-1))))]
    #[case::element_and_pred("O#h1", PartialAtomDsl(AtomAst::new(ElementAst::Lit(Element::O)).with_implicit_hydrogens(ValueAst::Lit(1))))]
    fn test_parse_partial_atom(#[case] input: &str, #[case] expected: PartialAtomDsl) {
        assert_eq!(parse_partial_atom(input).unwrap(), expected);
    }

    #[rstest]
    #[case::dup_hydrogens("#h1#h2", ParseError::DuplicateAtomPredicate("#h".to_string()))]
    fn test_parse_partial_atom_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(parse_partial_atom(input).unwrap_err(), expected);
    }

    #[rstest]
    #[case::charge_only(r##""#c-1""##, PartialAtomDsl(AtomAst::new(ElementAst::Undetermined).with_charge(ValueAst::Lit(-1))))]
    fn test_partial_atom_dsl_from_edn(#[case] input: &str, #[case] expected: PartialAtomDsl) {
        assert_eq!(
            PartialAtomDsl::from_edn(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::non_string("1")]
    fn test_partial_atom_dsl_from_edn_error(#[case] input: &str) {
        assert!(matches!(
            PartialAtomDsl::from_edn(&read_string(input).unwrap()),
            Err(DeError::TypeMismatch {
                expected: "string",
                ..
            })
        ));
    }

    #[rstest]
    #[case::charge_only(PartialAtomDsl(AtomAst::new(ElementAst::Undetermined).with_charge(ValueAst::Lit(-1))), r##""#c-""##)]
    fn test_partial_atom_dsl_to_edn(#[case] input: PartialAtomDsl, #[case] expected: &str) {
        assert_eq!(input.to_edn(), read_string(expected).unwrap());
    }

    #[rstest]
    #[case::arom_not_aromatic("C#a!")]
    #[case::arom_plus("C#a+")]
    #[case::arom_zero("C#a0")]
    #[case::arom_omit("C#a")]
    #[case::multicenter_not("C#m!")]
    #[case::multicenter_plus("C#m+")]
    #[case::multicenter_zero("C#m0")]
    #[case::multicenter_omit("C#m")]
    #[case::ring_membership_all_bare("C#R")]
    #[case::ring_membership_all_plus("C#R+")]
    #[case::ring_bang("C#R!")]
    #[case::ring_membership_all("C#R2")]
    #[case::ring_membership_size_conj("C#R(5)#R(6)")]
    #[case::element_not_lit("!H")]
    #[case::element_not_set("!{F,Cl}")]
    #[case::element_var_bare("?e")]
    #[case::element_var_domain("?e :: {C,N}")]
    #[case::element_var_domain_not_in("?e :: !{F,Cl}")]
    fn test_atom_display_roundtrip(#[case] input: &str) {
        let parsed = atom.parse(input).unwrap();
        assert_eq!(parsed.to_string(), input);
    }

    #[rstest]
    #[case::valence("C#v*", "C")]
    #[case::donated("C#d*", "C")]
    #[case::accepted("C#t*", "C")]
    #[case::degree("C#D*", "C")]
    #[case::total_degree("C#X*", "C")]
    #[case::ring_degree("C#x*", "C")]
    #[case::ring_valence("C#y*", "C")]
    #[case::total_valence("C#V*", "C")]
    #[case::total_h("C#H*", "C")]
    #[case::ring_membership_all("C#R*", "C")]
    #[case::ring_membership_size("C#R(6)*", "C")]
    #[case::aromatic_undetermined("C#a*", "C")]
    #[case::multicenter_undetermined("C#m*", "C")]
    #[case::tetrahedral_undetermined("C#T*", "C")]
    fn test_atom_render_vacuous_constraints(#[case] input: &str, #[case] expected_canonical: &str) {
        let parsed: AtomDsl = atom.parse(input).unwrap();
        assert_eq!(parsed.to_string(), expected_canonical);
        let reparsed: AtomDsl = atom.parse(&parsed.to_string()).unwrap();
        assert!(
            reparsed.0.constraints.is_empty(),
            "vacuous constraint should be absent after render → reparse, got {:?}",
            reparsed.0.constraints,
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::carbon("C", ElementAst::Lit(Element::C))]
    #[case::iron("Fe", ElementAst::Lit(Element::Fe))]
    #[case::chlorine("Cl", ElementAst::Lit(Element::Cl))]
    #[case::undetermined("*", ElementAst::Undetermined)]
    #[case::set("{C,N,O}", ElementAst::lit_set(vec![Element::C, Element::N, Element::O]))]
    #[case::set_spaced("{ C, N}", ElementAst::lit_set(vec![Element::C, Element::N]))]
    #[case::var_domain_paren("(?e :: {C,N})", ElementAst::var_in("e", vec![Element::C, Element::N]))]
    #[case::var_domain_bare("?e :: {C,N}", ElementAst::var_in("e", vec![Element::C, Element::N]))]
    #[case::var_domain_paren_paren("((?e :: {C,N}))", ElementAst::var_in("e", vec![Element::C, Element::N]))]
    #[case::var_domain_not_in_lit("?e :: !H", ElementAst::var_not_in("e", vec![Element::H]))]
    #[case::var_domain_not_in_set("?e :: !{F,Cl}", ElementAst::var_not_in("e", vec![Element::F, Element::Cl]))]
    #[case::var_domain_not_in_paren("(?e :: !{F,Cl})", ElementAst::var_not_in("e", vec![Element::F, Element::Cl]))]
    #[case::var_bare("?e", ElementAst::var("e".to_string()))]
    #[case::var_paren("(?e)", ElementAst::var("e".to_string()))]
    #[case::var_paren_paren("((?e))", ElementAst::var("e".to_string()))]
    #[case::not_lit("!H", ElementAst::not(Element::H))]
    #[case::not_set("!{F,Cl}", ElementAst::not_set(vec![Element::F, Element::Cl]))]
    fn test_element(#[case] input: &str, #[case] expected: ElementAst) {
        let result = element.parse(input);
        assert!(result.is_ok(), "{input:?} should succeed, got {:?}", result.unwrap_err());
        let expr = result.unwrap();
        assert_eq!(expr, expected);
    }

    #[rstest]
    #[case::empty("")]
    #[case::lowercase("c")]
    #[case::invalid("123")]
    #[case::unknown_element("Xx")]
    fn test_element_error(#[case] input: &str) {
        let result = element.parse(input);
        assert!(
            result.is_err(),
            "{input:?} should fail, got {:?}",
            result.unwrap()
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::natural("=", IsotopeMassAst::Natural)]
    #[case::lit("12", IsotopeMassAst::Lit(12))]
    #[case::undetermined("*", IsotopeMassAst::Undetermined)]
    #[case::set("{12,13,14}", IsotopeMassAst::lit_set([12, 13, 14]))]
    #[case::set_spaced("{ 12, 13 }", IsotopeMassAst::lit_set([12, 13]))]
    #[case::set_singleton("{12}", IsotopeMassAst::lit_set([12]))]
    #[case::var_in_paren("(?m :: {12,13})", IsotopeMassAst::var_in("m", [12, 13]))]
    #[case::var_in_bare("?m :: {12,13}", IsotopeMassAst::var_in("m", [12, 13]))]
    #[case::var_in_paren_paren("((?m :: {12,13}))", IsotopeMassAst::var_in("m", [12, 13]))]
    #[case::var_paren("(?m)", IsotopeMassAst::var("m"))]
    #[case::var_bare("?m", IsotopeMassAst::var("m"))]
    #[case::var_paren_paren("((?m))", IsotopeMassAst::var("m"))]
    fn test_isotope(#[case] input: &str, #[case] expected: IsotopeMassAst) {
        let result = isotope.parse(input);
        assert!(result.is_ok(), "{input:?} should succeed, got {:?}", result.unwrap_err());
        let expr = result.unwrap();
        assert_eq!(expr, expected);
    }

    #[rstest]
    #[case::empty("")]
    #[case::negative("-14")]
    #[case::not_lit("!14")]
    #[case::not_set("!{12,13}")]
    #[case::var_not_in("?m :: !14")]
    fn test_isotope_error(#[case] input: &str) {
        assert!(isotope.parse(input).is_err(), "{input:?} should fail");
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::isotope_lit("#i12", AtomPredicate::IsotopeMass(IsotopeMassAst::Lit(12)))]
    #[case::isotope_natural("#i=", AtomPredicate::IsotopeMass(IsotopeMassAst::Natural))]
    #[case::isotope_undetermined("#i*", AtomPredicate::IsotopeMass(IsotopeMassAst::Undetermined))]
    #[case::charge_pos("#c+2", AtomPredicate::Charge(ValueAst::Lit(2)))]
    #[case::charge_neg("#c-2", AtomPredicate::Charge(ValueAst::Lit(-2)))]
    #[case::charge_plus("#c+", AtomPredicate::Charge(ValueAst::Lit(1)))]
    #[case::charge_minus("#c-", AtomPredicate::Charge(ValueAst::Lit(-1)))]
    #[case::charge_zero("#c0", AtomPredicate::Charge(ValueAst::Lit(0)))]
    #[case::charge_undetermined("#c*", AtomPredicate::Charge(ValueAst::Undetermined))]
    #[case::h_count("#h3", AtomPredicate::ImplicitHydrogens(ValueAst::Lit(3)))]
    #[case::h_undetermined("#h*", AtomPredicate::ImplicitHydrogens(ValueAst::Undetermined))]
    #[case::h_omit("#h", AtomPredicate::ImplicitHydrogens(ValueAst::Lit(1)))]
    #[case::lone_pairs("#n2", AtomPredicate::LonePairs(ValueAst::Lit(2)))]
    #[case::lone_pairs_omit("#n", AtomPredicate::LonePairs(ValueAst::Lit(1)))]
    #[case::unpaired("#u2", AtomPredicate::Spin(SpinPredicate::Unpaired(ValueAst::Lit(2))))]
    #[case::unpaired_omit("#u", AtomPredicate::Spin(SpinPredicate::Unpaired(ValueAst::Lit(1))))]
    #[case::multiplicity("#s3", AtomPredicate::Spin(SpinPredicate::Multiplicity(ValueAst::Lit(3))))]
    #[case::multiplicity_omit("#s", AtomPredicate::Spin(SpinPredicate::Multiplicity(ValueAst::Lit(1))))]
    #[case::valence("#v4", AtomPredicate::Constraint(AtomConstraintAst::Valence(ValueAst::Lit(4))))]
    #[case::total_valence("#V5", AtomPredicate::Constraint(AtomConstraintAst::TotalValence(ValueAst::Lit(5))))]
    #[case::total_valence_omit("#V", AtomPredicate::Constraint(AtomConstraintAst::TotalValence(ValueAst::Lit(1))))]
    #[case::ring_valence("#y2", AtomPredicate::Constraint(AtomConstraintAst::RingValence(ValueAst::Lit(2))))]
    #[case::ring_valence_omit("#y", AtomPredicate::Constraint(AtomConstraintAst::RingValence(ValueAst::Lit(1))))]
    #[case::donated_pairs("#d1", AtomPredicate::Constraint(AtomConstraintAst::DonatedPairs(ValueAst::Lit(1))))]
    #[case::accepted_pairs("#t1", AtomPredicate::Constraint(AtomConstraintAst::AcceptedPairs(ValueAst::Lit(1))))]
    #[case::ring_membership_size("#R(6)", AtomPredicate::Constraint(AtomConstraintAst::ring_membership(RingScope::Size(6), 1)))]
    #[case::arom_not_aromatic("#a!", AtomPredicate::Constraint(AtomConstraintAst::AromaticValence(AromaticValenceAst::NotAromatic)))]
    #[case::arom_undetermined("#a*", AtomPredicate::Constraint(AtomConstraintAst::AromaticValence(AromaticValenceAst::Undetermined)))]
    #[case::arom_plus("#a+", AtomPredicate::Constraint(AtomConstraintAst::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Undetermined))))]
    #[case::arom_lit("#a2", AtomPredicate::Constraint(AtomConstraintAst::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(2)))))]
    #[case::arom_omit("#a", AtomPredicate::Constraint(AtomConstraintAst::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(1)))))]
    #[case::multicenter_not("#m!", AtomPredicate::Constraint(AtomConstraintAst::MulticenterValence(MulticenterValenceAst::NotMulticenter)))]
    #[case::multicenter_undetermined("#m*", AtomPredicate::Constraint(AtomConstraintAst::MulticenterValence(MulticenterValenceAst::Undetermined)))]
    #[case::multicenter_plus("#m+", AtomPredicate::Constraint(AtomConstraintAst::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Undetermined))))]
    #[case::multicenter_zero("#m0", AtomPredicate::Constraint(AtomConstraintAst::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(0)))))]
    #[case::multicenter_omit("#m", AtomPredicate::Constraint(AtomConstraintAst::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(1)))))]
    #[case::multicenter("#m2", AtomPredicate::Constraint(AtomConstraintAst::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(2)))))]
    #[case::degree("#D2", AtomPredicate::Constraint(AtomConstraintAst::Degree(ValueAst::Lit(2))))]
    #[case::degree_omit("#D", AtomPredicate::Constraint(AtomConstraintAst::Degree(ValueAst::Lit(1))))]
    #[case::total_degree("#X3", AtomPredicate::Constraint(AtomConstraintAst::TotalDegree(ValueAst::Lit(3))))]
    #[case::ring_degree("#x2", AtomPredicate::Constraint(AtomConstraintAst::RingDegree(ValueAst::Lit(2))))]
    #[case::ring_degree_omit("#x", AtomPredicate::Constraint(AtomConstraintAst::RingDegree(ValueAst::Lit(1))))]
    #[case::total_hydrogens("#H1", AtomPredicate::Constraint(AtomConstraintAst::TotalHydrogens(ValueAst::Lit(1))))]
    #[case::ring_membership_all_bare("#R", AtomPredicate::Constraint(AtomConstraintAst::ring_membership(RingScope::All, ValueAst::Lit(1))))]
    #[case::ring_membership_all_star("#R*", AtomPredicate::Constraint(AtomConstraintAst::ring_membership(RingScope::All, ValueAst::Undetermined)))]
    #[case::ring_membership_all_plus("#R+", AtomPredicate::Constraint(AtomConstraintAst::ring_membership(RingScope::All, ValueAst::RangeFrom(1))))]
    #[case::ring_membership_all("#R2", AtomPredicate::Constraint(AtomConstraintAst::ring_membership(RingScope::All, ValueAst::Lit(2))))]
    fn test_atom_predicate(#[case] input: &str, #[case] expected: AtomPredicate) {
        let result = atom_predicate.parse(input);
        assert!(result.is_ok(), "{input:?} should succeed, got {:?}", result.unwrap_err());
        let pred = result.unwrap();
        assert_eq!(pred, expected);
    }

    #[rstest]
    #[case::unknown_tag("#z", ParseError::UnknownAtomPredicate("#z".to_string()))]
    #[case::trailing_no_hash("fo", ParseError::TrailingInput("fo".to_string()))]
    fn test_atom_predicate_error(#[case] input: &str, #[case] expected: ParseError) {
        let result = atom_predicate.parse(input);
        assert!(
            result.is_err(),
            "{input:?} should fail, got {:?}",
            result.unwrap()
        );
        let err = result.unwrap_err().into_inner();
        assert_eq!(err, expected);
    }

    #[rstest]
    fn test_atom_dsl_from_ast() {
        let mut ast = AtomAst::new(ElementAst::Lit(Element::C));
        ast.charge = ValueAst::Lit(0);
        ast.lone_pairs = ValueAst::Lit(0);
        ast.implicit_hydrogens = ValueAst::Lit(0);
        ast.isotope_mass = IsotopeMassAst::Natural;
        ast.spin = SpinStateAst::from((0_u8, 1_u8));
        ast.constraints
            .set(AtomConstraintAst::Valence(ValueAst::Lit(0)));
        ast.constraints.set(AtomConstraintAst::AromaticValence(
            AromaticValenceAst::NotAromatic,
        ));
        let cfg = AtomDefaults::zeroed();
        let dsl = AtomDsl::from_ast(&ast, &cfg);
        assert_eq!(dsl.0.charge, ValueAst::Undetermined);
        assert_eq!(dsl.0.lone_pairs, ValueAst::Undetermined);
        assert_eq!(dsl.0.implicit_hydrogens, ValueAst::Undetermined);
        assert_eq!(dsl.0.isotope_mass, IsotopeMassAst::Undetermined);
        assert_eq!(dsl.0.spin, SpinStateAst::default());
        assert!(dsl.0.constraints.is_empty());
    }

    #[rstest]
    fn test_atom_dsl_into_ast() {
        let dsl = AtomDsl(AtomAst::new(ElementAst::Lit(Element::C)));
        let cfg = AtomDefaults::zeroed();
        let ast = dsl.into_ast(&cfg);
        assert_eq!(ast.charge, ValueAst::Lit(0));
        assert_eq!(ast.lone_pairs, ValueAst::Lit(0));
        assert_eq!(ast.implicit_hydrogens, ValueAst::Lit(0));
        assert_eq!(ast.isotope_mass, IsotopeMassAst::Natural);
        assert_eq!(ast.spin, SpinStateAst::from((0_u8, 1_u8)));
        assert_eq!(
            ast.constraints.get(AtomConstraintKey::Valence),
            Some(&AtomConstraintAst::Valence(ValueAst::Lit(0)))
        );
        assert_eq!(
            ast.constraints.get(AtomConstraintKey::AromaticValence),
            Some(&AtomConstraintAst::AromaticValence(
                AromaticValenceAst::NotAromatic
            ))
        );
    }

    #[rstest]
    fn test_raise_atom_constraints() {
        // A vacuous defaulted kind is overwritten with its default; a vacuous NON-defaulted kind
        // survives (no global vacuous strip); a concrete value is left alone.
        let mut constraints = AtomConstraintsAst::from_iter([
            AtomConstraintAst::Valence(ValueAst::Undetermined),
            AtomConstraintAst::TotalValence(ValueAst::Undetermined),
            AtomConstraintAst::degree(3),
        ]);
        raise_atom_constraints(&mut constraints, &AtomDefaults::zeroed());
        assert_eq!(
            constraints.get(AtomConstraintKey::Valence),
            Some(&AtomConstraintAst::Valence(ValueAst::Lit(0)))
        );
        assert_eq!(
            constraints.get(AtomConstraintKey::TotalValence),
            Some(&AtomConstraintAst::TotalValence(ValueAst::Undetermined))
        );
        assert_eq!(
            constraints.get(AtomConstraintKey::Degree),
            Some(&AtomConstraintAst::degree(3))
        );
    }

    #[rstest]
    fn test_lower_atom_constraints() {
        // A defaulted entry equal to its default is elided; a non-default value is kept.
        let mut constraints = AtomConstraintsAst::from_iter([
            AtomConstraintAst::Valence(ValueAst::Lit(0)),
            AtomConstraintAst::AromaticValence(AromaticValenceAst::NotAromatic),
            AtomConstraintAst::degree(3),
        ]);
        lower_atom_constraints(&mut constraints, &AtomDefaults::zeroed());
        assert_eq!(
            constraints.iter().cloned().collect::<Vec<_>>(),
            vec![AtomConstraintAst::degree(3)]
        );
    }

    #[rstest]
    fn test_atom_dsl_roundtrip_zeroed() {
        let input = AtomDsl(AtomAst::new(ElementAst::Lit(Element::C)));
        let cfg = AtomDefaults::zeroed();
        let raised = input.clone().into_ast(&cfg);
        let lowered = AtomDsl::from_ast(&raised, &cfg);
        assert_eq!(input, lowered);
    }

    #[rstest]
    #[case::simple(r##""C""##)]
    #[case::with_charge(r##""C#c+""##)]
    #[case::with_constraint(r##""N#v3#a""##)]
    #[case::with_tetrahedral_stereo(r##""C#T1""##)]
    fn test_atom_dsl_from_edn_str_matches_from_edn(#[case] input: &str) {
        let via_stream = AtomDsl::from_edn_str(input).unwrap();
        let tree = read_string(input).unwrap();
        let via_tree = AtomDsl::from_edn(&tree).unwrap();
        assert_eq!(via_stream, via_tree);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(AromaticValenceAst::Undetermined, ":undetermined")]
    #[case::not_aromatic(AromaticValenceAst::NotAromatic, ":not-aromatic")]
    #[case::aromatic_lit(AromaticValenceAst::Aromatic(ValueAst::Lit(6)), "{:aromatic 6}")]
    #[case::aromatic_undetermined(AromaticValenceAst::Aromatic(ValueAst::Undetermined), "{:aromatic :undetermined}")]
    #[case::aromatic_no_canonicalization(AromaticValenceAst::Aromatic(ValueAst::lit_set([5])), "{:aromatic [5]}")]
    fn test_aromatic_valence_dsl_roundtrip(
        #[case] input: AromaticValenceAst,
        #[case] edn_source: &str,
    ) {
        let dsl = AromaticValenceDsl::from_ast(&input, &());
        let edn = dsl.to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected);
        let parsed = AromaticValenceDsl::from_edn(&edn).unwrap();
        assert_eq!(parsed.into_ast(&()), input);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(MulticenterValenceAst::Undetermined, ":undetermined")]
    #[case::not_multicenter(MulticenterValenceAst::NotMulticenter, ":not-multicenter")]
    #[case::multicenter_lit(MulticenterValenceAst::Multicenter(ValueAst::Lit(3)), "{:multicenter 3}")]
    #[case::multicenter_no_canonicalization(MulticenterValenceAst::Multicenter(ValueAst::lit_set([5])), "{:multicenter [5]}")]
    fn test_multicenter_valence_dsl_roundtrip(
        #[case] input: MulticenterValenceAst,
        #[case] edn_source: &str,
    ) {
        let dsl = MulticenterValenceDsl::from_ast(&input, &());
        let edn = dsl.to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected);
        let parsed = MulticenterValenceDsl::from_edn(&edn).unwrap();
        assert_eq!(parsed.into_ast(&()), input);
    }

    #[rstest]
    fn test_aromatic_valence_dsl_rejects_unknown_key() {
        let edn = read_string("{:bogus 1}").unwrap();
        let err = AromaticValenceDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::UnknownField { .. }));
    }

    #[rstest]
    fn test_multicenter_valence_dsl_rejects_wrong_shape() {
        let err = MulticenterValenceDsl::from_edn(&Edn::Int(3)).unwrap_err();
        assert!(matches!(err, DeError::TypeMismatch { .. }));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(AtomConstraintAst::Valence(ValueAst::Lit(4)), "{:valence 4}")]
    #[case::valence_undetermined(AtomConstraintAst::Valence(ValueAst::Undetermined), "{:valence :undetermined}")]
    #[case::valence_set(AtomConstraintAst::Valence(ValueAst::lit_set([3, 4])), "{:valence [3 4]}")]
    #[case::degree(AtomConstraintAst::Degree(ValueAst::Lit(3)), "{:degree 3}")]
    #[case::total_degree(AtomConstraintAst::TotalDegree(ValueAst::Lit(4)), "{:total-degree 4}")]
    #[case::ring_degree(AtomConstraintAst::RingDegree(ValueAst::Lit(2)), "{:ring-degree 2}")]
    #[case::ring_valence(AtomConstraintAst::RingValence(ValueAst::Lit(3)), "{:ring-valence 3}")]
    #[case::total_valence(AtomConstraintAst::TotalValence(ValueAst::Lit(5)), "{:total-valence 5}")]
    #[case::total_h(AtomConstraintAst::TotalHydrogens(ValueAst::Lit(3)), "{:total-hydrogens 3}")]
    #[case::ring_membership_all(AtomConstraintAst::ring_membership(RingScope::All, ValueAst::Lit(1)), "{:ring-membership {:count 1}}")]
    #[case::ring_membership_size(AtomConstraintAst::ring_membership(RingScope::Size(6), 1), "{:ring-membership {:size 6 :count 1}}")]
    #[case::donated(AtomConstraintAst::DonatedPairs(ValueAst::Lit(1)), "{:donated-pairs 1}")]
    #[case::accepted(AtomConstraintAst::AcceptedPairs(ValueAst::Lit(2)), "{:accepted-pairs 2}")]
    #[case::aromatic_not(AtomConstraintAst::AromaticValence(AromaticValenceAst::NotAromatic), "{:aromatic-valence :not-aromatic}")]
    #[case::aromatic_value(AtomConstraintAst::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(6))), "{:aromatic-valence {:aromatic 6}}")]
    #[case::multicenter_not(AtomConstraintAst::MulticenterValence(MulticenterValenceAst::NotMulticenter), "{:multicenter-valence :not-multicenter}")]
    #[case::multicenter_value(AtomConstraintAst::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(3))), "{:multicenter-valence {:multicenter 3}}")]
    #[case::valence_expr(AtomConstraintAst::Valence(ValueAst::predicate(ValuePredicate::Rel(ValueTerm::Var("h".to_string()), RelOp::Ge, ValueTerm::Lit(1)))), "{:valence \"?h >= 1\"}")]
    #[case::ring_membership_size_count_set(AtomConstraintAst::ring_membership(RingScope::Size(6), ValueAst::lit_set([5, 6])), "{:ring-membership {:size 6 :count [5 6]}}")]
    #[case::tetrahedral_stereo_undetermined(AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::Undetermined), "{:tetrahedral-stereo :undetermined}")]
    #[case::tetrahedral_stereo_not_stereo(AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::NotStereo), "{:tetrahedral-stereo :not-stereo}")]
    #[case::tetrahedral_stereo_lit(AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(1))), "{:tetrahedral-stereo {:stereo 1}}")]
    #[case::tetrahedral_stereo_coset_undetermined(AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::Stereo(StereoCosetAst::Undetermined)), "{:tetrahedral-stereo {:stereo :undetermined}}")]
    #[case::tetrahedral_stereo_set(AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::Stereo(StereoCosetAst::lit_set([1, 2]))), "{:tetrahedral-stereo {:stereo [1 2]}}")]
    #[case::tetrahedral_stereo_term(AtomConstraintAst::TetrahedralStereo(TetrahedralStereoAst::Stereo(StereoCosetAst::term(StereoTerm::swap(StereoTerm::Lit(1))))), "{:tetrahedral-stereo {:stereo \"~1\"}}")]
    fn test_atom_constraint_dsl_roundtrip(#[case] input: AtomConstraintAst, #[case] edn_source: &str) {

        let dsl = AtomConstraintDsl::from_ast(&input, &());
        let edn = dsl.to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = AtomConstraintDsl::from_edn(&edn).unwrap();
        assert_eq!(parsed.into_ast(&()), input, "parse-back mismatch");
    }

    #[rstest]
    #[case::non_map(Edn::Int(3), DeError::TypeMismatch { expected: "atom-constraint single-key map", got: "int", path: vec![] })]
    #[case::multiple_keys(read_string("{:valence 4 :degree 3}").unwrap(), DeError::Custom("atom-constraint must have exactly one key, got 2".to_string()))]
    #[case::unknown_key(read_string("{:bogus 1}").unwrap(), DeError::UnknownField { key: "bogus".to_string(), path: vec!["atom-constraint".into()] })]
    fn test_atom_constraint_dsl_error(#[case] input: Edn<'static>, #[case] expected: DeError) {
        let err = AtomConstraintDsl::from_edn(&input).unwrap_err();
        assert_eq!(err, expected);
    }

    #[rstest]
    #[case::bare("C")]
    #[case::charged("N#c+")]
    #[case::aromatic("C#h3#a+")]
    fn test_atom_ast_from_str_to_string_roundtrip(#[case] s: &str) {
        let ast: AtomAst = s.parse().unwrap();
        assert_eq!(ast.to_string(), s);
    }

    #[rstest]
    fn test_atom_ast_from_str_carbon_element() {
        let ast: AtomAst = "C".parse().unwrap();
        assert_eq!(ast.element, ElementAst::Lit(Element::C));
        assert_eq!(ast.charge, ValueAst::Undetermined);
    }

    #[rstest]
    fn test_atom_ast_to_edn_roundtrip() {
        use umol_edn::ToEdn;
        let ast: AtomAst = "C#c+#h3".parse().unwrap();
        let edn = ast.to_edn();
        let back = AtomAst::from_edn(&edn).unwrap();
        assert_eq!(back, ast);
    }
}
