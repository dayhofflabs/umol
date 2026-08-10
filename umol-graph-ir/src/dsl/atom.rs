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
use super::num::{fmt_num, fmt_set, num, terminator, variable_name, NumDsl};
use super::operators::{mem_op, mem_op_str};
use super::predicate::{
    apply_unpaired_electrons_predicate, charge, fmt_charge, fmt_ring_membership,
    fmt_unpaired_electrons, lower_unpaired_electrons, optional_value, raise_unpaired_electrons,
    ring_membership, UnpairedElectronsPredicate,
};
use super::stereo::{
    fmt_tetrahedral_stereo_config, tetrahedral_stereo_config, TetrahedralStereoDsl,
};
use crate::ir::atom::{AtomForm, AtomUpdate, ElementForm, IsotopeMassForm};
use crate::ir::constraint::{
    AromaticValenceForm, AtomConstraintForm, AtomConstraintKey, AtomConstraintsForm,
    MulticenterValenceForm, RingScope,
};
use crate::ir::num::NumForm;
use crate::ir::operators::MemOp;
use crate::ir::stereo::TetrahedralStereoForm;
use crate::ir::traits::{FromIr, IntoIr, Lattice};

/// Surface DSL wrapper around `AtomForm`. Parses and renders the atom-string form
/// (element plus `#…` predicates); inline-capable constraints land in
/// `self.0.constraints`.
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AtomDsl(pub AtomForm);

impl AtomDsl {
    /// Zero-cost reference cast from `&AtomForm`. Relies on `repr(transparent)`.
    pub fn from_ref(form: &AtomForm) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(form as *const AtomForm as *const Self) }
    }
}

impl From<AtomForm> for AtomDsl {
    fn from(form: AtomForm) -> Self {
        Self(form)
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
        fmt_atom_form(f, &self.0)?;
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

impl FromIr<AtomForm> for AtomDsl {
    type Ctx = AtomDefaults;

    fn from_ir(form: &AtomForm, cfg: &Self::Ctx) -> Self {
        let mut out = form.clone();
        lower_atom(&mut out, cfg);
        AtomDsl(out)
    }
}

impl IntoIr<AtomForm> for AtomDsl {
    type Ctx = AtomDefaults;

    fn into_ir(mut self, cfg: &Self::Ctx) -> AtomForm {
        raise_atom(&mut self.0, cfg);
        self.0
    }
}

impl FromStr for AtomForm {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(AtomDsl::from_str(s)?.into_ir(&AtomDefaults::default()))
    }
}

impl Display for AtomForm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        AtomDsl::from_ref(self).fmt(f)
    }
}

impl<'de> FromEdn<'de> for AtomForm {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        Ok(AtomDsl::from_edn(edn)?.into_ir(&AtomDefaults::default()))
    }

    fn from_edn_str(input: &'de str) -> Result<Self, EdnError> {
        Ok(AtomDsl::from_edn_str(input)?.into_ir(&AtomDefaults::default()))
    }
}

impl ToEdn for AtomForm {
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
    Some(AtomDsl(AtomForm::from_element(el)))
}

/// Atom-string parser (does not require consuming all input).
pub(crate) fn atom(i: &mut &str) -> PResult<AtomDsl> {
    let el = delimited(multispace0, element, multispace0).parse_next(i)?;
    let preds: Vec<AtomPredicate> =
        repeat(0.., terminated(atom_predicate, multispace0)).parse_next(i)?;
    let mut form = AtomDsl(AtomForm::new(el));
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

/// Surface DSL wrapper around an [`AtomUpdate`].
#[repr(transparent)]
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AtomUpdateDsl(pub AtomUpdate);

impl AtomUpdateDsl {
    /// Zero-cost reference cast from `&AtomUpdate`. Relies on `repr(transparent)`.
    pub fn from_ref(update: &AtomUpdate) -> &Self {
        // SAFETY: `#[repr(transparent)]` guarantees identical layout.
        unsafe { &*(update as *const AtomUpdate as *const Self) }
    }
}

impl FromIr<AtomUpdate> for AtomUpdateDsl {
    type Ctx = ();

    fn from_ir(update: &AtomUpdate, _ctx: &Self::Ctx) -> Self {
        Self(update.clone())
    }
}

impl IntoIr<AtomUpdate> for AtomUpdateDsl {
    type Ctx = ();

    fn into_ir(self, _ctx: &Self::Ctx) -> AtomUpdate {
        self.0
    }
}

impl FromStr for AtomUpdateDsl {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_atom_update(s)
    }
}

impl FromStr for AtomUpdate {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(AtomUpdateDsl::from_str(s)?.into_ir(&()))
    }
}

impl Display for AtomUpdate {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        AtomUpdateDsl::from_ref(self).fmt(f)
    }
}

impl Display for AtomUpdateDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let update = &self.0;
        if let Some(element) = &update.element {
            fmt_element(f, element)?;
        }
        if let Some(isotope_mass) = &update.isotope_mass {
            if isotope_mass.is_undetermined() {
                write!(f, "#i*")?;
            } else {
                fmt_isotope_mass(f, isotope_mass)?;
            }
        }
        if let Some(charge) = &update.charge {
            if charge.is_undetermined() {
                write!(f, "#c*")?;
            } else {
                fmt_charge(f, charge)?;
            }
        }
        if let Some(implicit_hydrogens) = &update.implicit_hydrogens {
            fmt_update_value_field(f, "#h", implicit_hydrogens)?;
        }
        if let Some(lone_pairs) = &update.lone_pairs {
            fmt_update_value_field(f, "#n", lone_pairs)?;
        }
        if let Some(unpaired_electrons) = &update.unpaired_electrons.count {
            fmt_update_value_field(f, "#u", unpaired_electrons)?;
        }
        if let Some(multiplicity) = &update.unpaired_electrons.multiplicity {
            fmt_update_value_field(f, "#s", multiplicity)?;
        }
        for c in update.constraints.iter() {
            if c.is_undetermined() {
                fmt_undetermined_constraint(f, c)?;
            } else {
                fmt_constraint(f, c)?;
            }
        }
        Ok(())
    }
}

impl<'de> FromEdn<'de> for AtomUpdateDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => s.parse().map_err(|e| DeError::subgrammar("atom-update", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for AtomUpdateDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

pub fn parse_atom_update(input: &str) -> Result<AtomUpdateDsl, ParseError> {
    atom_update.parse(input).map_err(|e| e.into_inner())
}

pub(crate) fn atom_update(i: &mut &str) -> PResult<AtomUpdateDsl> {
    let element = delimited(multispace0, opt(element), multispace0).parse_next(i)?;
    let preds: Vec<AtomPredicate> =
        repeat(0.., terminated(atom_predicate, multispace0)).parse_next(i)?;
    let mut update = AtomUpdate {
        element,
        ..Default::default()
    };
    apply_update_predicates(&mut update, preds).map_err(ErrMode::Cut)?;
    Ok(AtomUpdateDsl(update))
}

fn apply_update_predicates(
    update: &mut AtomUpdate,
    preds: Vec<AtomPredicate>,
) -> Result<(), ParseError> {
    for pred in preds {
        match pred {
            AtomPredicate::IsotopeMass(value) => {
                if update.isotope_mass.replace(value).is_some() {
                    return Err(ParseError::DuplicateAtomPredicate("#i".to_string()));
                }
            }
            AtomPredicate::Charge(value) => {
                if update.charge.replace(value).is_some() {
                    return Err(ParseError::DuplicateAtomPredicate("#c".to_string()));
                }
            }
            AtomPredicate::ImplicitHydrogens(value) => {
                if update.implicit_hydrogens.replace(value).is_some() {
                    return Err(ParseError::DuplicateAtomPredicate("#h".to_string()));
                }
            }
            AtomPredicate::LonePairs(value) => {
                if update.lone_pairs.replace(value).is_some() {
                    return Err(ParseError::DuplicateAtomPredicate("#n".to_string()));
                }
            }
            AtomPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Count(value)) => {
                if update.unpaired_electrons.count.replace(value).is_some() {
                    return Err(ParseError::DuplicateAtomPredicate("#u".to_string()));
                }
            }
            AtomPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Multiplicity(value)) => {
                if update
                    .unpaired_electrons
                    .multiplicity
                    .replace(value)
                    .is_some()
                {
                    return Err(ParseError::DuplicateAtomPredicate("#s".to_string()));
                }
            }
            AtomPredicate::Constraint(constraint) => {
                if update.constraints.contains(constraint.key()) {
                    return Err(ParseError::DuplicateAtomPredicate(
                        constraint_tag(constraint.key()).to_string(),
                    ));
                }
                update.constraints.set(constraint);
            }
        }
    }
    Ok(())
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

fn fmt_undetermined_constraint(f: &mut fmt::Formatter<'_>, c: &AtomConstraintForm) -> fmt::Result {
    match c {
        AtomConstraintForm::RingMembership(membership) => match membership.scope {
            RingScope::All => write!(f, "#R*"),
            RingScope::Size(size) => write!(f, "#R({})*", size),
        },
        _ => write!(f, "{}*", constraint_tag(c.key())),
    }
}

/// One predicate from an atom-string; the parser yields a `Vec` of these
/// and the applier folds them into the `AtomForm`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AtomPredicate {
    IsotopeMass(IsotopeMassForm),
    Charge(NumForm),
    ImplicitHydrogens(NumForm),
    LonePairs(NumForm),
    UnpairedElectrons(UnpairedElectronsPredicate),
    Constraint(AtomConstraintForm),
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
            .map(|v| AtomPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Count(v)))
            .parse_next(i),
        "#s" => optional_value
            .map(|v| AtomPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Multiplicity(v)))
            .parse_next(i),
        "#v" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraintForm::Valence(v)))
            .parse_next(i),
        "#V" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraintForm::TotalValence(v)))
            .parse_next(i),
        "#d" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraintForm::DonatedPairs(v)))
            .parse_next(i),
        "#t" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraintForm::AcceptedPairs(v)))
            .parse_next(i),
        "#a" => aromatic_valence
            .map(|c| AtomPredicate::Constraint(AtomConstraintForm::AromaticValence(c)))
            .parse_next(i),
        "#m" => multicenter_valence
            .map(|c| AtomPredicate::Constraint(AtomConstraintForm::MulticenterValence(c)))
            .parse_next(i),
        "#D" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraintForm::Degree(v)))
            .parse_next(i),
        "#X" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraintForm::TotalDegree(v)))
            .parse_next(i),
        "#x" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraintForm::RingDegree(v)))
            .parse_next(i),
        "#y" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraintForm::RingValence(v)))
            .parse_next(i),
        "#H" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraintForm::TotalHydrogens(v)))
            .parse_next(i),
        "#R" => ring_membership
            .map(|m| AtomPredicate::Constraint(AtomConstraintForm::RingMembership(m)))
            .parse_next(i),
        "#T" => (|i: &mut &str| tetrahedral_stereo_config(i))
            .map(|c| AtomPredicate::Constraint(AtomConstraintForm::TetrahedralStereo(c)))
            .parse_next(i),
        p if p.starts_with('#') => Err(ErrMode::Cut(ParseError::UnknownAtomPredicate(
            p.to_string(),
        ))),
        _ => Err(ErrMode::Cut(ParseError::TrailingInput(start.to_string()))),
    }
}

fn element(i: &mut &str) -> PResult<ElementForm> {
    alt((
        '*'.value(ElementForm::Undetermined),
        preceded('!', element_set).map(ElementForm::not_set),
        preceded('!', element_literal).map(ElementForm::not),
        element_set.map(ElementForm::lit_set),
        element_bind.map(|(name, set, mem_op)| match mem_op {
            MemOp::In => ElementForm::var_in(name, set),
            MemOp::NotIn => ElementForm::var_not_in(name, set),
        }),
        element_ref.map(ElementForm::var),
        element_literal.map(ElementForm::Lit),
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
            preceded('?', variable_name),
            delimited(multispace0, mem_op, multispace0),
            element_bind_domain,
        )
            .map(|(name, op, set)| (name, set, op)),
    ))
    .parse_next(i)
}

fn element_bind_domain(i: &mut &str) -> PResult<Vec<Element>> {
    alt((element_set, element_literal.map(|e| vec![e]))).parse_next(i)
}

fn element_ref(i: &mut &str) -> PResult<String> {
    alt((
        delimited('(', delimited(multispace0, element_ref, multispace0), ')'),
        preceded('?', variable_name),
    ))
    .parse_next(i)
}

fn isotope(i: &mut &str) -> PResult<IsotopeMassForm> {
    preceded(
        multispace0,
        alt((
            '='.value(IsotopeMassForm::Natural),
            '*'.value(IsotopeMassForm::Undetermined),
            isotope_set.map(IsotopeMassForm::lit_set),
            terminated(isotope_var, (multispace0, terminator)).map(|(name, domain)| match domain {
                Some(set) => IsotopeMassForm::var_in(name, set),
                None => IsotopeMassForm::var(name),
            }),
            terminated(dec_uint::<_, u32, _>, (multispace0, terminator)).map(IsotopeMassForm::Lit),
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
            preceded('?', variable_name),
            opt(preceded(
                delimited(multispace0, "::", multispace0),
                isotope_set,
            )),
        ),
    ))
    .parse_next(i)
}

fn aromatic_valence(i: &mut &str) -> PResult<AromaticValenceForm> {
    preceded(
        multispace0,
        alt((
            "*".value(AromaticValenceForm::Undetermined),
            "!".value(AromaticValenceForm::NotAromatic),
            // `#a+` encodes "aromatic, count unspecified" — structurally
            // distinct from the outer Undetermined; canonical form is
            // Aromatic(Undetermined).
            "+".value(AromaticValenceForm::Aromatic(NumForm::Undetermined)),
            num.map(AromaticValenceForm::Aromatic),
            empty.value(AromaticValenceForm::Aromatic(NumForm::Lit(1))),
        )),
    )
    .parse_next(i)
    .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedPredicateBody))
}

fn multicenter_valence(i: &mut &str) -> PResult<MulticenterValenceForm> {
    preceded(
        multispace0,
        alt((
            "*".value(MulticenterValenceForm::Undetermined),
            "!".value(MulticenterValenceForm::NotMulticenter),
            // `#m+` mirrors `#a+` — "multicenter, count unspecified".
            "+".value(MulticenterValenceForm::Multicenter(NumForm::Undetermined)),
            num.map(MulticenterValenceForm::Multicenter),
            empty.value(MulticenterValenceForm::Multicenter(NumForm::Lit(1))),
        )),
    )
    .parse_next(i)
    .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedPredicateBody))
}

fn apply_predicates(dsl: &mut AtomDsl, preds: Vec<AtomPredicate>) -> Result<(), ParseError> {
    let atom = &mut dsl.0;
    for pred in preds {
        match pred {
            AtomPredicate::IsotopeMass(v) => {
                if !matches!(atom.isotope_mass, IsotopeMassForm::Undetermined) {
                    return Err(ParseError::DuplicateAtomPredicate("#i".to_string()));
                }
                atom.isotope_mass = v;
            }
            AtomPredicate::Charge(v) => {
                if !matches!(atom.charge, NumForm::Undetermined) {
                    return Err(ParseError::DuplicateAtomPredicate("#c".to_string()));
                }
                atom.charge = v;
            }
            AtomPredicate::ImplicitHydrogens(v) => {
                if !matches!(atom.implicit_hydrogens, NumForm::Undetermined) {
                    return Err(ParseError::DuplicateAtomPredicate("#h".to_string()));
                }
                atom.implicit_hydrogens = v;
            }
            AtomPredicate::LonePairs(v) => {
                if !matches!(atom.lone_pairs, NumForm::Undetermined) {
                    return Err(ParseError::DuplicateAtomPredicate("#n".to_string()));
                }
                atom.lone_pairs = v;
            }
            AtomPredicate::UnpairedElectrons(predicate) => {
                apply_unpaired_electrons_predicate(
                    &mut atom.unpaired_electrons,
                    predicate,
                    ParseError::DuplicateAtomPredicate,
                )?;
            }
            AtomPredicate::Constraint(c) => {
                if atom.constraints.contains(c.key()) {
                    return Err(ParseError::DuplicateAtomPredicate(
                        constraint_tag(c.key()).to_string(),
                    ));
                }
                atom.constraints.set(c);
            }
        }
    }
    Ok(())
}

fn fmt_atom_form(f: &mut fmt::Formatter<'_>, form: &AtomForm) -> fmt::Result {
    fmt_element(f, &form.element)?;
    fmt_isotope_mass(f, &form.isotope_mass)?;
    fmt_charge(f, &form.charge)?;
    fmt_num_field(f, "#h", &form.implicit_hydrogens)?;
    fmt_num_field(f, "#n", &form.lone_pairs)?;
    fmt_unpaired_electrons(f, &form.unpaired_electrons)
}

fn fmt_element(f: &mut fmt::Formatter<'_>, expr: &ElementForm) -> fmt::Result {
    match expr {
        ElementForm::Lit(e) => write!(f, "{}", e),
        ElementForm::Undetermined => write!(f, "*"),
        ElementForm::LitSet(es) => fmt_element_set(f, es.iter().copied()),
        ElementForm::NotSet(es) if es.len() == 1 => {
            // A singleton complement renders as `!e` (no braces), e.g. `!H`.
            write!(f, "!{}", es.iter().next().unwrap())
        }
        ElementForm::NotSet(es) => {
            write!(f, "!")?;
            fmt_element_set(f, es.iter().copied())
        }
        ElementForm::Var(v) => {
            let (name, domain) = &**v;
            write!(f, "?{}", name)?;
            if let Some((op, set)) = domain {
                write!(f, " {} ", mem_op_str(*op))?;
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

fn fmt_isotope_mass(f: &mut fmt::Formatter<'_>, iso: &IsotopeMassForm) -> fmt::Result {
    match iso {
        IsotopeMassForm::Undetermined => Ok(()),
        IsotopeMassForm::Natural => write!(f, "#i="),
        IsotopeMassForm::Lit(n) => write!(f, "#i{}", n),
        IsotopeMassForm::LitSet(s) => {
            write!(f, "#i")?;
            fmt_set(f, s.iter().copied())
        }
        IsotopeMassForm::Var(v) => {
            let (name, domain) = &**v;
            write!(f, "#i?{}", name)?;
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
fn fmt_num_field(f: &mut fmt::Formatter<'_>, prefix: &str, v: &NumForm) -> fmt::Result {
    match v {
        NumForm::Undetermined => Ok(()),
        NumForm::Lit(1) => write!(f, "{}", prefix),
        NumForm::Lit(n) => write!(f, "{}{}", prefix, n),
        v => {
            write!(f, "{}", prefix)?;
            fmt_num(f, v)
        }
    }
}

fn fmt_update_value_field(f: &mut fmt::Formatter<'_>, prefix: &str, v: &NumForm) -> fmt::Result {
    if v.is_undetermined() {
        write!(f, "{}*", prefix)
    } else {
        fmt_num_field(f, prefix, v)
    }
}

/// Format an inline-constraint value field. Per the canonical-rendering
/// rules in `dsl::predicates`, vacuous constraints (`Undetermined`) elide.
/// `Lit(0)` is a meaningful constraint and renders.
fn fmt_num_field_required(f: &mut fmt::Formatter<'_>, prefix: &str, v: &NumForm) -> fmt::Result {
    match v {
        NumForm::Undetermined => Ok(()),
        NumForm::Lit(1) => write!(f, "{}", prefix),
        NumForm::Lit(n) => write!(f, "{}{}", prefix, n),
        v => {
            write!(f, "{}", prefix)?;
            fmt_num(f, v)
        }
    }
}

fn fmt_constraint(f: &mut fmt::Formatter<'_>, c: &AtomConstraintForm) -> fmt::Result {
    match c {
        AtomConstraintForm::Valence(v) => fmt_num_field_required(f, "#v", v),
        AtomConstraintForm::DonatedPairs(v) => fmt_num_field_required(f, "#d", v),
        AtomConstraintForm::AcceptedPairs(v) => fmt_num_field_required(f, "#t", v),
        AtomConstraintForm::MulticenterValence(c) => match c {
            MulticenterValenceForm::Undetermined => Ok(()),
            MulticenterValenceForm::NotMulticenter => write!(f, "#m!"),
            MulticenterValenceForm::Multicenter(NumForm::Undetermined) => write!(f, "#m+"),
            MulticenterValenceForm::Multicenter(NumForm::Lit(1)) => write!(f, "#m"),
            MulticenterValenceForm::Multicenter(NumForm::Lit(n)) => write!(f, "#m{}", n),
            MulticenterValenceForm::Multicenter(v) => {
                write!(f, "#m")?;
                fmt_num(f, v)
            }
        },
        AtomConstraintForm::AromaticValence(c) => match c {
            AromaticValenceForm::Undetermined => Ok(()),
            AromaticValenceForm::NotAromatic => write!(f, "#a!"),
            AromaticValenceForm::Aromatic(NumForm::Undetermined) => write!(f, "#a+"),
            AromaticValenceForm::Aromatic(NumForm::Lit(1)) => write!(f, "#a"),
            AromaticValenceForm::Aromatic(NumForm::Lit(n)) => write!(f, "#a{}", n),
            AromaticValenceForm::Aromatic(v) => {
                write!(f, "#a")?;
                fmt_num(f, v)
            }
        },
        AtomConstraintForm::Degree(v) => fmt_num_field_required(f, "#D", v),
        AtomConstraintForm::TotalDegree(v) => fmt_num_field_required(f, "#X", v),
        AtomConstraintForm::RingDegree(v) => fmt_num_field_required(f, "#x", v),
        AtomConstraintForm::RingValence(v) => fmt_num_field_required(f, "#y", v),
        AtomConstraintForm::TotalValence(v) => fmt_num_field_required(f, "#V", v),
        AtomConstraintForm::TotalHydrogens(v) => fmt_num_field_required(f, "#H", v),
        AtomConstraintForm::RingMembership(m) => fmt_ring_membership(f, m),
        AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Undetermined) => Ok(()),
        AtomConstraintForm::TetrahedralStereo(c) => {
            write!(f, "#T")?;
            fmt_tetrahedral_stereo_config(f, c)
        }
    }
}

pub(crate) fn raise_atom(atom: &mut AtomForm, cfg: &AtomDefaults) {
    // Exhaustive destructure: adding a new AtomForm field is a compile error
    // here, forcing the author to decide how raising should handle it.
    let AtomForm {
        element: _,
        isotope_mass,
        charge,
        implicit_hydrogens,
        lone_pairs,
        unpaired_electrons,
        constraints,
    } = atom;

    if matches!(*isotope_mass, IsotopeMassForm::Undetermined) {
        *isotope_mass = match cfg.isotope {
            IsotopeDefault::Natural => IsotopeMassForm::Natural,
            IsotopeDefault::Required => IsotopeMassForm::Undetermined,
        };
    }
    if matches!(*charge, NumForm::Undetermined) {
        *charge = match cfg.charge {
            NumericDefault::Zero => NumForm::Lit(0),
            NumericDefault::Required => NumForm::Undetermined,
        };
    }
    if matches!(*implicit_hydrogens, NumForm::Undetermined) {
        *implicit_hydrogens = match cfg.implicit_hydrogens {
            NumericDefault::Zero => NumForm::Lit(0),
            NumericDefault::Required => NumForm::Undetermined,
        };
    }
    if matches!(*lone_pairs, NumForm::Undetermined) {
        *lone_pairs = match cfg.lone_pairs {
            NumericDefault::Zero => NumForm::Lit(0),
            NumericDefault::Required => NumForm::Undetermined,
        };
    }
    raise_unpaired_electrons(unpaired_electrons, cfg.unpaired_electrons, cfg.multiplicity);
    raise_atom_constraints(constraints, cfg);
}

/// A defaulted key wants filling iff it is absent or holds a vacuous (`Undetermined`) value.
/// A concrete user value is left alone.
fn is_unset_or_vacuous(constraints: &AtomConstraintsForm, key: AtomConstraintKey) -> bool {
    constraints.get(key).is_none_or(|c| c.is_undetermined())
}

fn raise_atom_constraints(constraints: &mut AtomConstraintsForm, cfg: &AtomDefaults) {
    // One explicit clause per defaulted kind, in ascending key-sort order. No global vacuous
    // strip: a defaulted kind fills its own absent/vacuous entry; vacuous entries of other kinds
    // are left for lazy canonicalization.
    if matches!(cfg.valence, NumericDefault::Zero)
        && is_unset_or_vacuous(constraints, AtomConstraintKey::Valence)
    {
        constraints.set(AtomConstraintForm::Valence(NumForm::Lit(0)));
    }
    if matches!(cfg.aromatic_valence, AromaticValenceDefault::NotAromatic)
        && is_unset_or_vacuous(constraints, AtomConstraintKey::AromaticValence)
    {
        constraints.set(AtomConstraintForm::AromaticValence(
            AromaticValenceForm::NotAromatic,
        ));
    }
    if matches!(
        cfg.multicenter_valence,
        MulticenterValenceDefault::NotMulticenter
    ) && is_unset_or_vacuous(constraints, AtomConstraintKey::MulticenterValence)
    {
        constraints.set(AtomConstraintForm::MulticenterValence(
            MulticenterValenceForm::NotMulticenter,
        ));
    }
    if matches!(cfg.donated_pairs, NumericDefault::Zero)
        && is_unset_or_vacuous(constraints, AtomConstraintKey::DonatedPairs)
    {
        constraints.set(AtomConstraintForm::DonatedPairs(NumForm::Lit(0)));
    }
    if matches!(cfg.accepted_pairs, NumericDefault::Zero)
        && is_unset_or_vacuous(constraints, AtomConstraintKey::AcceptedPairs)
    {
        constraints.set(AtomConstraintForm::AcceptedPairs(NumForm::Lit(0)));
    }
    if matches!(cfg.tetrahedral_stereo, StereoDefault::NotStereo)
        && is_unset_or_vacuous(constraints, AtomConstraintKey::TetrahedralStereo)
    {
        constraints.set(AtomConstraintForm::TetrahedralStereo(
            TetrahedralStereoForm::NotStereo,
        ));
    }
}

pub(crate) fn lower_atom(atom: &mut AtomForm, cfg: &AtomDefaults) {
    // Exhaustive destructure: adding a new AtomForm field is a compile error
    // here, forcing the author to decide how lowering should handle it.
    let AtomForm {
        element: _,
        isotope_mass,
        charge,
        implicit_hydrogens,
        lone_pairs,
        unpaired_electrons,
        constraints,
    } = atom;

    if matches!(
        (&cfg.isotope, &*isotope_mass),
        (IsotopeDefault::Natural, IsotopeMassForm::Natural)
    ) {
        *isotope_mass = IsotopeMassForm::Undetermined;
    }
    if matches!(
        (&cfg.charge, &*charge),
        (NumericDefault::Zero, NumForm::Lit(0))
    ) {
        *charge = NumForm::Undetermined;
    }
    match (&cfg.implicit_hydrogens, &*implicit_hydrogens) {
        (NumericDefault::Required, NumForm::Undetermined) => {
            *implicit_hydrogens = NumForm::Undetermined;
        }
        (NumericDefault::Zero, NumForm::Lit(0)) => {
            *implicit_hydrogens = NumForm::Undetermined;
        }
        _ => {}
    }
    if matches!(
        (&cfg.lone_pairs, &*lone_pairs),
        (NumericDefault::Zero, NumForm::Lit(0))
    ) {
        *lone_pairs = NumForm::Undetermined;
    }
    lower_unpaired_electrons(unpaired_electrons, cfg.unpaired_electrons, cfg.multiplicity);
    lower_atom_constraints(constraints, cfg);
}

/// Elide a defaulted entry: if its key holds exactly `default`, remove it (raise would refill it
/// from the same `cfg`). No-op otherwise.
fn remove_if_default(constraints: &mut AtomConstraintsForm, default: AtomConstraintForm) {
    let key = default.key();
    if constraints.get(key) == Some(&default) {
        constraints.remove(key);
    }
}

fn lower_atom_constraints(constraints: &mut AtomConstraintsForm, cfg: &AtomDefaults) {
    // Elide each defaulted entry equal to its default, in reverse key-sort order (mirror of raise).
    if matches!(cfg.tetrahedral_stereo, StereoDefault::NotStereo) {
        remove_if_default(
            constraints,
            AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::NotStereo),
        );
    }
    if matches!(cfg.accepted_pairs, NumericDefault::Zero) {
        remove_if_default(
            constraints,
            AtomConstraintForm::AcceptedPairs(NumForm::Lit(0)),
        );
    }
    if matches!(cfg.donated_pairs, NumericDefault::Zero) {
        remove_if_default(
            constraints,
            AtomConstraintForm::DonatedPairs(NumForm::Lit(0)),
        );
    }
    if matches!(
        cfg.multicenter_valence,
        MulticenterValenceDefault::NotMulticenter
    ) {
        remove_if_default(
            constraints,
            AtomConstraintForm::MulticenterValence(MulticenterValenceForm::NotMulticenter),
        );
    }
    if matches!(cfg.aromatic_valence, AromaticValenceDefault::NotAromatic) {
        remove_if_default(
            constraints,
            AtomConstraintForm::AromaticValence(AromaticValenceForm::NotAromatic),
        );
    }
    if matches!(cfg.valence, NumericDefault::Zero) {
        remove_if_default(constraints, AtomConstraintForm::Valence(NumForm::Lit(0)));
    }
}

/// Surface DSL wrapper around `AromaticValenceForm`. EDN form: `:undetermined`,
/// `:not-aromatic`, or `{:aromatic <value>}`.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct AromaticValenceDsl(pub AromaticValenceForm);

impl<'de> FromEdn<'de> for AromaticValenceDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Keyword(k) if k.name() == "undetermined" => {
                Ok(Self(AromaticValenceForm::Undetermined))
            }
            Edn::Keyword(k) if k.name() == "not-aromatic" => {
                Ok(Self(AromaticValenceForm::NotAromatic))
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
                    "aromatic" => Ok(Self(AromaticValenceForm::Aromatic(
                        NumDsl::from_edn(v)?.into_ir(&()),
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
            AromaticValenceForm::Undetermined => {
                Edn::Keyword(EdnKeyword::owned("undetermined".into()))
            }
            AromaticValenceForm::NotAromatic => {
                Edn::Keyword(EdnKeyword::owned("not-aromatic".into()))
            }
            AromaticValenceForm::Aromatic(v) => {
                single_key_map("aromatic", NumDsl::from_ir(v, &()).to_edn())
            }
        }
    }
}

impl FromIr<AromaticValenceForm> for AromaticValenceDsl {
    type Ctx = ();

    fn from_ir(form: &AromaticValenceForm, _ctx: &Self::Ctx) -> Self {
        Self(form.clone())
    }
}

impl IntoIr<AromaticValenceForm> for AromaticValenceDsl {
    type Ctx = ();

    fn into_ir(self, _ctx: &Self::Ctx) -> AromaticValenceForm {
        self.0
    }
}

/// Surface DSL wrapper around `MulticenterValenceForm`. EDN form:
/// `:undetermined`, `:not-multicenter`, or `{:multicenter <value>}`.
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct MulticenterValenceDsl(pub MulticenterValenceForm);

impl<'de> FromEdn<'de> for MulticenterValenceDsl {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Keyword(k) if k.name() == "undetermined" => {
                Ok(Self(MulticenterValenceForm::Undetermined))
            }
            Edn::Keyword(k) if k.name() == "not-multicenter" => {
                Ok(Self(MulticenterValenceForm::NotMulticenter))
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
                    "multicenter" => Ok(Self(MulticenterValenceForm::Multicenter(
                        NumDsl::from_edn(v)?.into_ir(&()),
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
            MulticenterValenceForm::Undetermined => {
                Edn::Keyword(EdnKeyword::owned("undetermined".into()))
            }
            MulticenterValenceForm::NotMulticenter => {
                Edn::Keyword(EdnKeyword::owned("not-multicenter".into()))
            }
            MulticenterValenceForm::Multicenter(v) => {
                single_key_map("multicenter", NumDsl::from_ir(v, &()).to_edn())
            }
        }
    }
}

impl FromIr<MulticenterValenceForm> for MulticenterValenceDsl {
    type Ctx = ();

    fn from_ir(form: &MulticenterValenceForm, _ctx: &Self::Ctx) -> Self {
        Self(form.clone())
    }
}

impl IntoIr<MulticenterValenceForm> for MulticenterValenceDsl {
    type Ctx = ();

    fn into_ir(self, _ctx: &Self::Ctx) -> MulticenterValenceForm {
        self.0
    }
}

/// Surface DSL wrapper around `AtomConstraintForm`. EDN form is a single-key map
/// keyed by the constraint kind: e.g. `{:valence 4}`, `{:degree *}`,
/// `{:aromatic-valence :not-aromatic}`, `{:total-hydrogens {?h :: {0,1,2}}}`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct AtomConstraintDsl(pub AtomConstraintForm);

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
            "valence" => AtomConstraintForm::Valence(NumDsl::from_edn(v)?.into_ir(&())),
            "total-valence" => AtomConstraintForm::TotalValence(NumDsl::from_edn(v)?.into_ir(&())),
            "aromatic-valence" => {
                AtomConstraintForm::AromaticValence(AromaticValenceDsl::from_edn(v)?.into_ir(&()))
            }
            "multicenter-valence" => AtomConstraintForm::MulticenterValence(
                MulticenterValenceDsl::from_edn(v)?.into_ir(&()),
            ),
            "donated-pairs" => AtomConstraintForm::DonatedPairs(NumDsl::from_edn(v)?.into_ir(&())),
            "accepted-pairs" => {
                AtomConstraintForm::AcceptedPairs(NumDsl::from_edn(v)?.into_ir(&()))
            }
            "degree" => AtomConstraintForm::Degree(NumDsl::from_edn(v)?.into_ir(&())),
            "total-degree" => AtomConstraintForm::TotalDegree(NumDsl::from_edn(v)?.into_ir(&())),
            "ring-degree" => AtomConstraintForm::RingDegree(NumDsl::from_edn(v)?.into_ir(&())),
            "ring-valence" => AtomConstraintForm::RingValence(NumDsl::from_edn(v)?.into_ir(&())),
            "total-hydrogens" => {
                AtomConstraintForm::TotalHydrogens(NumDsl::from_edn(v)?.into_ir(&()))
            }
            "ring-membership" => {
                AtomConstraintForm::RingMembership(RingMembershipDsl::from_edn(v)?.0)
            }
            "tetrahedral-stereo" => AtomConstraintForm::TetrahedralStereo(
                TetrahedralStereoDsl::from_edn(v)?.into_ir(&()),
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
            AtomConstraintForm::Valence(v) => {
                single_key_map("valence", NumDsl::from_ir(v, &()).to_edn())
            }
            AtomConstraintForm::TotalValence(v) => {
                single_key_map("total-valence", NumDsl::from_ir(v, &()).to_edn())
            }
            AtomConstraintForm::AromaticValence(c) => single_key_map(
                "aromatic-valence",
                AromaticValenceDsl::from_ir(c, &()).to_edn(),
            ),
            AtomConstraintForm::MulticenterValence(c) => single_key_map(
                "multicenter-valence",
                MulticenterValenceDsl::from_ir(c, &()).to_edn(),
            ),
            AtomConstraintForm::DonatedPairs(v) => {
                single_key_map("donated-pairs", NumDsl::from_ir(v, &()).to_edn())
            }
            AtomConstraintForm::AcceptedPairs(v) => {
                single_key_map("accepted-pairs", NumDsl::from_ir(v, &()).to_edn())
            }
            AtomConstraintForm::Degree(v) => {
                single_key_map("degree", NumDsl::from_ir(v, &()).to_edn())
            }
            AtomConstraintForm::TotalDegree(v) => {
                single_key_map("total-degree", NumDsl::from_ir(v, &()).to_edn())
            }
            AtomConstraintForm::RingDegree(v) => {
                single_key_map("ring-degree", NumDsl::from_ir(v, &()).to_edn())
            }
            AtomConstraintForm::RingValence(v) => {
                single_key_map("ring-valence", NumDsl::from_ir(v, &()).to_edn())
            }
            AtomConstraintForm::TotalHydrogens(v) => {
                single_key_map("total-hydrogens", NumDsl::from_ir(v, &()).to_edn())
            }
            AtomConstraintForm::RingMembership(m) => {
                single_key_map("ring-membership", RingMembershipDsl(m.clone()).to_edn())
            }
            AtomConstraintForm::TetrahedralStereo(c) => single_key_map(
                "tetrahedral-stereo",
                TetrahedralStereoDsl::from_ir(c, &()).to_edn(),
            ),
        }
    }
}

impl FromIr<AtomConstraintForm> for AtomConstraintDsl {
    type Ctx = ();

    fn from_ir(form: &AtomConstraintForm, _ctx: &Self::Ctx) -> Self {
        Self(form.clone())
    }
}

impl IntoIr<AtomConstraintForm> for AtomConstraintDsl {
    type Ctx = ();

    fn into_ir(self, _ctx: &Self::Ctx) -> AtomConstraintForm {
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
    use crate::ir::constraint::RingScope;
    use crate::ir::num::{ArithExpr, PredExpr};
    use crate::ir::operators::{MemOp, RelOp};
    use crate::ir::spin::{UnpairedElectronsForm, UnpairedElectronsUpdate};
    use crate::ir::stereo::{StereoCoset, StereoTerm};

    #[rstest]
    #[case::single("C", AtomDsl(AtomForm::new(ElementForm::Lit(Element::C))))]
    #[case::double("Cl", AtomDsl(AtomForm::new(ElementForm::Lit(Element::Cl))))]
    #[case::transuranic("Og", AtomDsl(AtomForm::new(ElementForm::Lit(Element::Og))))]
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
    #[case::carbon("C", AtomDsl(AtomForm::new(ElementForm::Lit(Element::C))))]
    #[case::iron("Fe", AtomDsl(AtomForm::new(ElementForm::Lit(Element::Fe))))]
    #[case::chlorine("Cl", AtomDsl(AtomForm::new(ElementForm::Lit(Element::Cl))))]
    #[case::whitespace("  C  ", AtomDsl(AtomForm::new(ElementForm::Lit(Element::C))))]
    #[case::undetermined("*", AtomDsl(AtomForm::new(ElementForm::Undetermined)))]
    #[case::element_set("{C,N,O}", AtomDsl(AtomForm::new(ElementForm::lit_set(vec![Element::C, Element::N, Element::O]))))]
    #[case::element_set_singleton("{C}", AtomDsl(AtomForm::new(ElementForm::lit_set(vec![Element::C]))))]
    #[case::element_var_mem_in("(?e :: {C,N})", AtomDsl(AtomForm::new(ElementForm::var_in("e", vec![Element::C, Element::N]))))]
    #[case::element_var("?e", AtomDsl(AtomForm::new(ElementForm::var("e".to_string()))))]
    #[case::isotope("C#i12", AtomDsl(AtomForm { isotope_mass: IsotopeMassForm::Lit(12), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::isotope_natural("C#i=", AtomDsl(AtomForm { isotope_mass: IsotopeMassForm::Natural, ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::charge_pos("C#c+2", AtomDsl(AtomForm { charge: NumForm::Lit(2), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::charge_neg("C#c-2", AtomDsl(AtomForm { charge: NumForm::Lit(-2), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::charge_plus("C#c+", AtomDsl(AtomForm { charge: NumForm::Lit(1), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::charge_minus("C#c-", AtomDsl(AtomForm { charge: NumForm::Lit(-1), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::charge_zero("C#c0", AtomDsl(AtomForm { charge: NumForm::Lit(0), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::h_count("C#h3", AtomDsl(AtomForm { implicit_hydrogens: NumForm::Lit(3), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::h_undetermined("C#h*", AtomDsl(AtomForm { implicit_hydrogens: NumForm::Undetermined, ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::h_bind("C#h(?h)", AtomDsl(AtomForm { implicit_hydrogens: NumForm::var("h"), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::h_set("N#h?h :: {2,3}", AtomDsl(AtomForm { implicit_hydrogens: NumForm::pred_expr(PredExpr::Mem(ArithExpr::Var("h".to_string()), MemOp::In, BTreeSet::from([2, 3]))), ..AtomForm::new(ElementForm::Lit(Element::N)) }))]
    #[case::h_expr("C#h?h >= 1", AtomDsl(AtomForm { implicit_hydrogens: NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Ge, ArithExpr::Lit(1))), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::h_omit("C#h", AtomDsl(AtomForm { implicit_hydrogens: NumForm::Lit(1), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::lone_pairs("O#n2", AtomDsl(AtomForm { lone_pairs: NumForm::Lit(2), ..AtomForm::new(ElementForm::Lit(Element::O)) }))]
    #[case::lone_pairs_omit("O#n", AtomDsl(AtomForm { lone_pairs: NumForm::Lit(1), ..AtomForm::new(ElementForm::Lit(Element::O)) }))]
    #[case::unpaired_electrons("C#u2", AtomDsl(AtomForm { unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Undetermined }, ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::unpaired_no_canonicalization("C#u{1}", AtomDsl(AtomForm { unpaired_electrons: UnpairedElectronsForm { count: NumForm::lit_set([1]), multiplicity: NumForm::Undetermined }, ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::unpaired_omit("C#u", AtomDsl(AtomForm { unpaired_electrons: UnpairedElectronsForm { count: NumForm::Lit(1), multiplicity: NumForm::Undetermined }, ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::multiplicity("C#s3", AtomDsl(AtomForm { unpaired_electrons: UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(3) }, ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::multiplicity_omit("C#s", AtomDsl(AtomForm { unpaired_electrons: UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(1) }, ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::valence("C#v4", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::Valence(NumForm::Lit(4))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::donated_pairs("N#d1", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::DonatedPairs(NumForm::Lit(1))]), ..AtomForm::new(ElementForm::Lit(Element::N)) }))]
    #[case::accepted_pairs("B#t1", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::AcceptedPairs(NumForm::Lit(1))]), ..AtomForm::new(ElementForm::Lit(Element::B)) }))]
    #[case::ring_membership_size("C#R(6)", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(RingScope::Size(6), 1)]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::arom_not_aromatic("C#a!", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::AromaticValence(AromaticValenceForm::NotAromatic)]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::arom_undetermined("C#a*", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::AromaticValence(AromaticValenceForm::Undetermined)]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::arom_plus("C#a+", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::AromaticValence(AromaticValenceForm::Aromatic(NumForm::Undetermined))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::arom_zero("C#a0", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::AromaticValence(AromaticValenceForm::Aromatic(NumForm::Lit(0)))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::arom_one("C#a1", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::AromaticValence(AromaticValenceForm::Aromatic(NumForm::Lit(1)))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::arom_omit("C#a", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::AromaticValence(AromaticValenceForm::Aromatic(NumForm::Lit(1)))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::multicenter_not("C#m!", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::MulticenterValence(MulticenterValenceForm::NotMulticenter)]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::multicenter_undetermined("C#m*", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::MulticenterValence(MulticenterValenceForm::Undetermined)]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::multicenter_plus("C#m+", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::MulticenterValence(MulticenterValenceForm::Multicenter(NumForm::Undetermined))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::multicenter_zero("C#m0", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::MulticenterValence(MulticenterValenceForm::Multicenter(NumForm::Lit(0)))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::multicenter_one("C#m", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::MulticenterValence(MulticenterValenceForm::Multicenter(NumForm::Lit(1)))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::multicenter("C#m2", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::MulticenterValence(MulticenterValenceForm::Multicenter(NumForm::Lit(2)))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::degree("C#D2", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::Degree(NumForm::Lit(2))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::total_degree("C#X3", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::TotalDegree(NumForm::Lit(3))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::ring_degree("C#x2", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::RingDegree(NumForm::Lit(2))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::ring_valence("C#y3", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::RingValence(NumForm::Lit(3))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::total_valence("C#V5", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::TotalValence(NumForm::Lit(5))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::total_hydrogens("C#H1", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::TotalHydrogens(NumForm::Lit(1))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::ring_membership_all_bare("C#R", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(RingScope::All, NumForm::Lit(1))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::ring_membership_all_star("C#R*", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(RingScope::All, NumForm::Undetermined)]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::ring_membership_all_plus("C#R+", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(RingScope::All, NumForm::RangeFrom(1))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::ring_not_in_ring("C#R!", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(RingScope::All, NumForm::Lit(0))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::ring_zero("C#R0", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(RingScope::All, NumForm::Lit(0))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::ring_membership_all("C#R2", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(RingScope::All, NumForm::Lit(2))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::ring_membership_size_conj("C#R(5)#R(6)", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::ring_membership(RingScope::Size(5), 1), AtomConstraintForm::ring_membership(RingScope::Size(6), 1)]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::tetrahedral_stereo_stereo("C#T+", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(StereoCoset::Undetermined))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::tetrahedral_stereo_lit("C#T1", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(StereoCoset::Lit(1)))]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::tetrahedral_stereo_not_stereo("C#T!", AtomDsl(AtomForm { constraints: AtomConstraintsForm::from_iter([AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::NotStereo)]), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::whitespace_before_predicate("C #h3", AtomDsl(AtomForm { implicit_hydrogens: NumForm::Lit(3), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::whitespace_between_predicates("C#c+ #h3", AtomDsl(AtomForm { charge: NumForm::Lit(1), implicit_hydrogens: NumForm::Lit(3), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
    #[case::whitespace_surrounding_predicates("  C  #c+  #h3  ", AtomDsl(AtomForm { charge: NumForm::Lit(1), implicit_hydrogens: NumForm::Lit(3), ..AtomForm::new(ElementForm::Lit(Element::C)) }))]
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

    #[rustfmt::skip]
    #[rstest]
    #[case::empty("", AtomUpdateDsl(AtomUpdate::default()))]
    #[case::element_only("O", AtomUpdateDsl(AtomUpdate { element: Some(ElementForm::Lit(Element::O)), ..Default::default() }))]
    #[case::charge_only("#c-1", AtomUpdateDsl(AtomUpdate { charge: Some(NumForm::Lit(-1)), ..Default::default() }))]
    #[case::element_and_pred("O#h1", AtomUpdateDsl(AtomUpdate { element: Some(ElementForm::Lit(Element::O)), implicit_hydrogens: Some(NumForm::Lit(1)), ..Default::default() }))]
    #[case::unpaired_electrons_unpaired("#u2", AtomUpdateDsl(AtomUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Lit(2)), multiplicity: None }, ..Default::default() }))]
    #[case::unpaired_electrons_multiplicity("#s1", AtomUpdateDsl(AtomUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) }, ..Default::default() }))]
    #[case::explicit_undetermined("*#i*#c*#h*#n*#u*#s*", AtomUpdateDsl(AtomUpdate { element: Some(ElementForm::Undetermined), isotope_mass: Some(IsotopeMassForm::Undetermined), charge: Some(NumForm::Undetermined), implicit_hydrogens: Some(NumForm::Undetermined), lone_pairs: Some(NumForm::Undetermined), unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Undetermined), multiplicity: Some(NumForm::Undetermined) }, constraints: Default::default() }))]
    #[case::constraint_removal("#v*", AtomUpdateDsl(AtomUpdate { constraints: AtomConstraintsForm::from(AtomConstraintForm::valence(NumForm::Undetermined)), ..Default::default() }))]
    fn test_parse_atom_update(#[case] input: &str, #[case] expected: AtomUpdateDsl) {
        assert_eq!(parse_atom_update(input).unwrap(), expected);
    }

    #[rstest]
    #[case::dup_hydrogens("#h1#h2", ParseError::DuplicateAtomPredicate("#h".to_string()))]
    #[case::dup_undetermined_charge("#c*#c1", ParseError::DuplicateAtomPredicate("#c".to_string()))]
    #[case::duplicate_undetermined_count("#u*#u1", ParseError::DuplicateAtomPredicate("#u".to_string()))]
    fn test_parse_atom_update_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(parse_atom_update(input).unwrap_err(), expected);
    }

    #[rstest]
    #[case::duplicate_hydrogens("#h1#h2", ParseError::DuplicateAtomPredicate("#h".to_string()))]
    fn test_atom_update_from_str_error(#[case] input: &str, #[case] expected: ParseError) {
        assert_eq!(input.parse::<AtomUpdate>().unwrap_err(), expected);
    }

    #[rstest]
    #[case::charge_only(r##""#c-1""##, AtomUpdateDsl(AtomUpdate { charge: Some(NumForm::Lit(-1)), ..Default::default() }))]
    fn test_atom_update_dsl_from_edn(#[case] input: &str, #[case] expected: AtomUpdateDsl) {
        assert_eq!(
            AtomUpdateDsl::from_edn(&read_string(input).unwrap()).unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::non_string("1")]
    fn test_atom_update_dsl_from_edn_error(#[case] input: &str) {
        assert!(matches!(
            AtomUpdateDsl::from_edn(&read_string(input).unwrap()),
            Err(DeError::TypeMismatch {
                expected: "string",
                ..
            })
        ));
    }

    #[rstest]
    #[case::charge_only(AtomUpdateDsl(AtomUpdate { charge: Some(NumForm::Lit(-1)), ..Default::default() }), r##""#c-""##)]
    #[case::unpaired_electrons_multiplicity(AtomUpdateDsl(AtomUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) }, ..Default::default() }), r##""#s""##)]
    #[case::explicit_undetermined(AtomUpdateDsl(AtomUpdate { element: Some(ElementForm::Undetermined), isotope_mass: Some(IsotopeMassForm::Undetermined), charge: Some(NumForm::Undetermined), implicit_hydrogens: Some(NumForm::Undetermined), lone_pairs: Some(NumForm::Undetermined), unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Undetermined), multiplicity: Some(NumForm::Undetermined) }, constraints: Default::default() }), r##""*#i*#c*#h*#n*#u*#s*""##)]
    #[case::ring_size_removal(AtomUpdateDsl(AtomUpdate { constraints: AtomConstraintsForm::from(AtomConstraintForm::ring_membership(RingScope::Size(3), NumForm::Undetermined)), ..Default::default() }), r##""#R(3)*""##)]
    fn test_atom_update_dsl_to_edn(#[case] input: AtomUpdateDsl, #[case] expected: &str) {
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
    #[case::element_var_domain_not_in("?e !: {F,Cl}")]
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

    #[rstest]
    #[case::charge_before_h("C#h3#c+", "C#c+#h3")]
    #[case::aromatic_before_ring("C#R2#a2", "C#a2#R2")]
    #[case::multicenter_before_stereo("C#T1#m2", "C#m2#T1")]
    #[case::stereo_before_ring("C#R2#T1", "C#T1#R2")]
    fn test_atom_render_canonical_order(#[case] input: &str, #[case] expected_canonical: &str) {
        let parsed: AtomDsl = atom.parse(input).unwrap();
        assert_eq!(parsed.to_string(), expected_canonical);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::carbon("C", ElementForm::Lit(Element::C))]
    #[case::iron("Fe", ElementForm::Lit(Element::Fe))]
    #[case::chlorine("Cl", ElementForm::Lit(Element::Cl))]
    #[case::undetermined("*", ElementForm::Undetermined)]
    #[case::set("{C,N,O}", ElementForm::lit_set(vec![Element::C, Element::N, Element::O]))]
    #[case::set_spaced("{ C, N}", ElementForm::lit_set(vec![Element::C, Element::N]))]
    #[case::var_domain_paren("(?e :: {C,N})", ElementForm::var_in("e", vec![Element::C, Element::N]))]
    #[case::var_domain_bare("?e :: {C,N}", ElementForm::var_in("e", vec![Element::C, Element::N]))]
    #[case::var_domain_paren_paren("((?e :: {C,N}))", ElementForm::var_in("e", vec![Element::C, Element::N]))]
    #[case::var_domain_not_in_lit("?e !: H", ElementForm::var_not_in("e", vec![Element::H]))]
    #[case::var_domain_not_in_set("?e !: {F,Cl}", ElementForm::var_not_in("e", vec![Element::F, Element::Cl]))]
    #[case::var_domain_not_in_paren("(?e !: {F,Cl})", ElementForm::var_not_in("e", vec![Element::F, Element::Cl]))]
    #[case::var_bare("?e", ElementForm::var("e".to_string()))]
    #[case::var_paren("(?e)", ElementForm::var("e".to_string()))]
    #[case::var_paren_paren("((?e))", ElementForm::var("e".to_string()))]
    #[case::not_lit("!H", ElementForm::not(Element::H))]
    #[case::not_set("!{F,Cl}", ElementForm::not_set(vec![Element::F, Element::Cl]))]
    fn test_element(#[case] input: &str, #[case] expected: ElementForm) {
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
    #[case::natural("=", IsotopeMassForm::Natural)]
    #[case::lit("12", IsotopeMassForm::Lit(12))]
    #[case::undetermined("*", IsotopeMassForm::Undetermined)]
    #[case::set("{12,13,14}", IsotopeMassForm::lit_set([12, 13, 14]))]
    #[case::set_spaced("{ 12, 13 }", IsotopeMassForm::lit_set([12, 13]))]
    #[case::set_singleton("{12}", IsotopeMassForm::lit_set([12]))]
    #[case::var_in_paren("(?m :: {12,13})", IsotopeMassForm::var_in("m", [12, 13]))]
    #[case::var_in_bare("?m :: {12,13}", IsotopeMassForm::var_in("m", [12, 13]))]
    #[case::var_in_paren_paren("((?m :: {12,13}))", IsotopeMassForm::var_in("m", [12, 13]))]
    #[case::var_paren("(?m)", IsotopeMassForm::var("m"))]
    #[case::var_bare("?m", IsotopeMassForm::var("m"))]
    #[case::var_paren_paren("((?m))", IsotopeMassForm::var("m"))]
    fn test_isotope(#[case] input: &str, #[case] expected: IsotopeMassForm) {
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
    #[case::isotope_lit("#i12", AtomPredicate::IsotopeMass(IsotopeMassForm::Lit(12)))]
    #[case::isotope_natural("#i=", AtomPredicate::IsotopeMass(IsotopeMassForm::Natural))]
    #[case::isotope_undetermined("#i*", AtomPredicate::IsotopeMass(IsotopeMassForm::Undetermined))]
    #[case::charge_pos("#c+2", AtomPredicate::Charge(NumForm::Lit(2)))]
    #[case::charge_neg("#c-2", AtomPredicate::Charge(NumForm::Lit(-2)))]
    #[case::charge_plus("#c+", AtomPredicate::Charge(NumForm::Lit(1)))]
    #[case::charge_minus("#c-", AtomPredicate::Charge(NumForm::Lit(-1)))]
    #[case::charge_zero("#c0", AtomPredicate::Charge(NumForm::Lit(0)))]
    #[case::charge_undetermined("#c*", AtomPredicate::Charge(NumForm::Undetermined))]
    #[case::h_count("#h3", AtomPredicate::ImplicitHydrogens(NumForm::Lit(3)))]
    #[case::h_undetermined("#h*", AtomPredicate::ImplicitHydrogens(NumForm::Undetermined))]
    #[case::h_omit("#h", AtomPredicate::ImplicitHydrogens(NumForm::Lit(1)))]
    #[case::lone_pairs("#n2", AtomPredicate::LonePairs(NumForm::Lit(2)))]
    #[case::lone_pairs_omit("#n", AtomPredicate::LonePairs(NumForm::Lit(1)))]
    #[case::unpaired_electrons("#u2", AtomPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Count(NumForm::Lit(2))))]
    #[case::unpaired_omit("#u", AtomPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Count(NumForm::Lit(1))))]
    #[case::multiplicity("#s3", AtomPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Multiplicity(NumForm::Lit(3))))]
    #[case::multiplicity_omit("#s", AtomPredicate::UnpairedElectrons(UnpairedElectronsPredicate::Multiplicity(NumForm::Lit(1))))]
    #[case::valence("#v4", AtomPredicate::Constraint(AtomConstraintForm::Valence(NumForm::Lit(4))))]
    #[case::total_valence("#V5", AtomPredicate::Constraint(AtomConstraintForm::TotalValence(NumForm::Lit(5))))]
    #[case::total_valence_omit("#V", AtomPredicate::Constraint(AtomConstraintForm::TotalValence(NumForm::Lit(1))))]
    #[case::ring_valence("#y2", AtomPredicate::Constraint(AtomConstraintForm::RingValence(NumForm::Lit(2))))]
    #[case::ring_valence_omit("#y", AtomPredicate::Constraint(AtomConstraintForm::RingValence(NumForm::Lit(1))))]
    #[case::donated_pairs("#d1", AtomPredicate::Constraint(AtomConstraintForm::DonatedPairs(NumForm::Lit(1))))]
    #[case::accepted_pairs("#t1", AtomPredicate::Constraint(AtomConstraintForm::AcceptedPairs(NumForm::Lit(1))))]
    #[case::ring_membership_size("#R(6)", AtomPredicate::Constraint(AtomConstraintForm::ring_membership(RingScope::Size(6), 1)))]
    #[case::arom_not_aromatic("#a!", AtomPredicate::Constraint(AtomConstraintForm::AromaticValence(AromaticValenceForm::NotAromatic)))]
    #[case::arom_undetermined("#a*", AtomPredicate::Constraint(AtomConstraintForm::AromaticValence(AromaticValenceForm::Undetermined)))]
    #[case::arom_plus("#a+", AtomPredicate::Constraint(AtomConstraintForm::AromaticValence(AromaticValenceForm::Aromatic(NumForm::Undetermined))))]
    #[case::arom_lit("#a2", AtomPredicate::Constraint(AtomConstraintForm::AromaticValence(AromaticValenceForm::Aromatic(NumForm::Lit(2)))))]
    #[case::arom_omit("#a", AtomPredicate::Constraint(AtomConstraintForm::AromaticValence(AromaticValenceForm::Aromatic(NumForm::Lit(1)))))]
    #[case::multicenter_not("#m!", AtomPredicate::Constraint(AtomConstraintForm::MulticenterValence(MulticenterValenceForm::NotMulticenter)))]
    #[case::multicenter_undetermined("#m*", AtomPredicate::Constraint(AtomConstraintForm::MulticenterValence(MulticenterValenceForm::Undetermined)))]
    #[case::multicenter_plus("#m+", AtomPredicate::Constraint(AtomConstraintForm::MulticenterValence(MulticenterValenceForm::Multicenter(NumForm::Undetermined))))]
    #[case::multicenter_zero("#m0", AtomPredicate::Constraint(AtomConstraintForm::MulticenterValence(MulticenterValenceForm::Multicenter(NumForm::Lit(0)))))]
    #[case::multicenter_omit("#m", AtomPredicate::Constraint(AtomConstraintForm::MulticenterValence(MulticenterValenceForm::Multicenter(NumForm::Lit(1)))))]
    #[case::multicenter("#m2", AtomPredicate::Constraint(AtomConstraintForm::MulticenterValence(MulticenterValenceForm::Multicenter(NumForm::Lit(2)))))]
    #[case::degree("#D2", AtomPredicate::Constraint(AtomConstraintForm::Degree(NumForm::Lit(2))))]
    #[case::degree_omit("#D", AtomPredicate::Constraint(AtomConstraintForm::Degree(NumForm::Lit(1))))]
    #[case::total_degree("#X3", AtomPredicate::Constraint(AtomConstraintForm::TotalDegree(NumForm::Lit(3))))]
    #[case::ring_degree("#x2", AtomPredicate::Constraint(AtomConstraintForm::RingDegree(NumForm::Lit(2))))]
    #[case::ring_degree_omit("#x", AtomPredicate::Constraint(AtomConstraintForm::RingDegree(NumForm::Lit(1))))]
    #[case::total_hydrogens("#H1", AtomPredicate::Constraint(AtomConstraintForm::TotalHydrogens(NumForm::Lit(1))))]
    #[case::ring_membership_all_bare("#R", AtomPredicate::Constraint(AtomConstraintForm::ring_membership(RingScope::All, NumForm::Lit(1))))]
    #[case::ring_membership_all_star("#R*", AtomPredicate::Constraint(AtomConstraintForm::ring_membership(RingScope::All, NumForm::Undetermined)))]
    #[case::ring_membership_all_plus("#R+", AtomPredicate::Constraint(AtomConstraintForm::ring_membership(RingScope::All, NumForm::RangeFrom(1))))]
    #[case::ring_membership_all("#R2", AtomPredicate::Constraint(AtomConstraintForm::ring_membership(RingScope::All, NumForm::Lit(2))))]
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
        let mut form = AtomForm::new(ElementForm::Lit(Element::C));
        form.charge = NumForm::Lit(0);
        form.lone_pairs = NumForm::Lit(0);
        form.implicit_hydrogens = NumForm::Lit(0);
        form.isotope_mass = IsotopeMassForm::Natural;
        form.unpaired_electrons = UnpairedElectronsForm::from((0_u8, 1_u8));
        form.constraints
            .set(AtomConstraintForm::Valence(NumForm::Lit(0)));
        form.constraints.set(AtomConstraintForm::AromaticValence(
            AromaticValenceForm::NotAromatic,
        ));
        let cfg = AtomDefaults::zeroed();
        let dsl = AtomDsl::from_ir(&form, &cfg);
        assert_eq!(dsl.0.charge, NumForm::Undetermined);
        assert_eq!(dsl.0.lone_pairs, NumForm::Undetermined);
        assert_eq!(dsl.0.implicit_hydrogens, NumForm::Undetermined);
        assert_eq!(dsl.0.isotope_mass, IsotopeMassForm::Undetermined);
        assert_eq!(dsl.0.unpaired_electrons, UnpairedElectronsForm::default());
        assert!(dsl.0.constraints.is_empty());
    }

    #[rstest]
    fn test_atom_dsl_into_ast() {
        let dsl = AtomDsl(AtomForm::new(ElementForm::Lit(Element::C)));
        let cfg = AtomDefaults::zeroed();
        let form = dsl.into_ir(&cfg);
        assert_eq!(form.charge, NumForm::Lit(0));
        assert_eq!(form.lone_pairs, NumForm::Lit(0));
        assert_eq!(form.implicit_hydrogens, NumForm::Lit(0));
        assert_eq!(form.isotope_mass, IsotopeMassForm::Natural);
        assert_eq!(
            form.unpaired_electrons,
            UnpairedElectronsForm::from((0_u8, 1_u8))
        );
        assert_eq!(
            form.constraints.get(AtomConstraintKey::Valence),
            Some(&AtomConstraintForm::Valence(NumForm::Lit(0)))
        );
        assert_eq!(
            form.constraints.get(AtomConstraintKey::AromaticValence),
            Some(&AtomConstraintForm::AromaticValence(
                AromaticValenceForm::NotAromatic
            ))
        );
    }

    #[rstest]
    fn test_raise_atom_constraints() {
        // A vacuous defaulted kind is overwritten with its default; a vacuous NON-defaulted kind
        // survives (no global vacuous strip); a concrete value is left alone.
        let mut constraints = AtomConstraintsForm::from_iter([
            AtomConstraintForm::Valence(NumForm::Undetermined),
            AtomConstraintForm::TotalValence(NumForm::Undetermined),
            AtomConstraintForm::degree(3),
        ]);
        raise_atom_constraints(&mut constraints, &AtomDefaults::zeroed());
        assert_eq!(
            constraints.get(AtomConstraintKey::Valence),
            Some(&AtomConstraintForm::Valence(NumForm::Lit(0)))
        );
        assert_eq!(
            constraints.get(AtomConstraintKey::TotalValence),
            Some(&AtomConstraintForm::TotalValence(NumForm::Undetermined))
        );
        assert_eq!(
            constraints.get(AtomConstraintKey::Degree),
            Some(&AtomConstraintForm::degree(3))
        );
    }

    #[rstest]
    fn test_lower_atom_constraints() {
        // A defaulted entry equal to its default is elided; a non-default value is kept.
        let mut constraints = AtomConstraintsForm::from_iter([
            AtomConstraintForm::Valence(NumForm::Lit(0)),
            AtomConstraintForm::AromaticValence(AromaticValenceForm::NotAromatic),
            AtomConstraintForm::degree(3),
        ]);
        lower_atom_constraints(&mut constraints, &AtomDefaults::zeroed());
        assert_eq!(
            constraints.iter().cloned().collect::<Vec<_>>(),
            vec![AtomConstraintForm::degree(3)]
        );
    }

    #[rstest]
    fn test_atom_dsl_roundtrip_zeroed() {
        let input = AtomDsl(AtomForm::new(ElementForm::Lit(Element::C)));
        let cfg = AtomDefaults::zeroed();
        let raised = input.clone().into_ir(&cfg);
        let lowered = AtomDsl::from_ir(&raised, &cfg);
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
    #[case::undetermined(AromaticValenceForm::Undetermined, ":undetermined")]
    #[case::not_aromatic(AromaticValenceForm::NotAromatic, ":not-aromatic")]
    #[case::aromatic_lit(AromaticValenceForm::Aromatic(NumForm::Lit(6)), "{:aromatic 6}")]
    #[case::aromatic_undetermined(AromaticValenceForm::Aromatic(NumForm::Undetermined), "{:aromatic :undetermined}")]
    #[case::aromatic_no_canonicalization(AromaticValenceForm::Aromatic(NumForm::lit_set([5])), "{:aromatic [5]}")]
    fn test_aromatic_valence_dsl_roundtrip(
        #[case] input: AromaticValenceForm,
        #[case] edn_source: &str,
    ) {
        let dsl = AromaticValenceDsl::from_ir(&input, &());
        let edn = dsl.to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected);
        let parsed = AromaticValenceDsl::from_edn(&edn).unwrap();
        assert_eq!(parsed.into_ir(&()), input);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(MulticenterValenceForm::Undetermined, ":undetermined")]
    #[case::not_multicenter(MulticenterValenceForm::NotMulticenter, ":not-multicenter")]
    #[case::multicenter_lit(MulticenterValenceForm::Multicenter(NumForm::Lit(3)), "{:multicenter 3}")]
    #[case::multicenter_no_canonicalization(MulticenterValenceForm::Multicenter(NumForm::lit_set([5])), "{:multicenter [5]}")]
    fn test_multicenter_valence_dsl_roundtrip(
        #[case] input: MulticenterValenceForm,
        #[case] edn_source: &str,
    ) {
        let dsl = MulticenterValenceDsl::from_ir(&input, &());
        let edn = dsl.to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected);
        let parsed = MulticenterValenceDsl::from_edn(&edn).unwrap();
        assert_eq!(parsed.into_ir(&()), input);
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
    #[case::valence(AtomConstraintForm::Valence(NumForm::Lit(4)), "{:valence 4}")]
    #[case::valence_undetermined(AtomConstraintForm::Valence(NumForm::Undetermined), "{:valence :undetermined}")]
    #[case::valence_set(AtomConstraintForm::Valence(NumForm::lit_set([3, 4])), "{:valence [3 4]}")]
    #[case::degree(AtomConstraintForm::Degree(NumForm::Lit(3)), "{:degree 3}")]
    #[case::total_degree(AtomConstraintForm::TotalDegree(NumForm::Lit(4)), "{:total-degree 4}")]
    #[case::ring_degree(AtomConstraintForm::RingDegree(NumForm::Lit(2)), "{:ring-degree 2}")]
    #[case::ring_valence(AtomConstraintForm::RingValence(NumForm::Lit(3)), "{:ring-valence 3}")]
    #[case::total_valence(AtomConstraintForm::TotalValence(NumForm::Lit(5)), "{:total-valence 5}")]
    #[case::total_h(AtomConstraintForm::TotalHydrogens(NumForm::Lit(3)), "{:total-hydrogens 3}")]
    #[case::ring_membership_all(AtomConstraintForm::ring_membership(RingScope::All, NumForm::Lit(1)), "{:ring-membership {:count 1}}")]
    #[case::ring_membership_size(AtomConstraintForm::ring_membership(RingScope::Size(6), 1), "{:ring-membership {:size 6 :count 1}}")]
    #[case::donated(AtomConstraintForm::DonatedPairs(NumForm::Lit(1)), "{:donated-pairs 1}")]
    #[case::accepted(AtomConstraintForm::AcceptedPairs(NumForm::Lit(2)), "{:accepted-pairs 2}")]
    #[case::aromatic_not(AtomConstraintForm::AromaticValence(AromaticValenceForm::NotAromatic), "{:aromatic-valence :not-aromatic}")]
    #[case::aromatic_value(AtomConstraintForm::AromaticValence(AromaticValenceForm::Aromatic(NumForm::Lit(6))), "{:aromatic-valence {:aromatic 6}}")]
    #[case::multicenter_not(AtomConstraintForm::MulticenterValence(MulticenterValenceForm::NotMulticenter), "{:multicenter-valence :not-multicenter}")]
    #[case::multicenter_value(AtomConstraintForm::MulticenterValence(MulticenterValenceForm::Multicenter(NumForm::Lit(3))), "{:multicenter-valence {:multicenter 3}}")]
    #[case::valence_expr(AtomConstraintForm::Valence(NumForm::pred_expr(PredExpr::Rel(ArithExpr::Var("h".to_string()), RelOp::Ge, ArithExpr::Lit(1)))), "{:valence \"?h >= 1\"}")]
    #[case::ring_membership_size_count_set(AtomConstraintForm::ring_membership(RingScope::Size(6), NumForm::lit_set([5, 6])), "{:ring-membership {:size 6 :count [5 6]}}")]
    #[case::tetrahedral_stereo_undetermined(AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Undetermined), "{:tetrahedral-stereo :undetermined}")]
    #[case::tetrahedral_stereo_not_stereo(AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::NotStereo), "{:tetrahedral-stereo :not-stereo}")]
    #[case::tetrahedral_stereo_lit(AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(StereoCoset::Lit(1))), "{:tetrahedral-stereo {:stereo 1}}")]
    #[case::tetrahedral_stereo_coset_undetermined(AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(StereoCoset::Undetermined)), "{:tetrahedral-stereo {:stereo :undetermined}}")]
    #[case::tetrahedral_stereo_set(AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(StereoCoset::lit_set([1, 2]))), "{:tetrahedral-stereo {:stereo [1 2]}}")]
    #[case::tetrahedral_stereo_term(AtomConstraintForm::TetrahedralStereo(TetrahedralStereoForm::Stereo(StereoCoset::term(StereoTerm::swap(StereoTerm::Lit(1))))), "{:tetrahedral-stereo {:stereo \"~1\"}}")]
    fn test_atom_constraint_dsl_roundtrip(#[case] input: AtomConstraintForm, #[case] edn_source: &str) {

        let dsl = AtomConstraintDsl::from_ir(&input, &());
        let edn = dsl.to_edn();
        let expected = read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = AtomConstraintDsl::from_edn(&edn).unwrap();
        assert_eq!(parsed.into_ir(&()), input, "parse-back mismatch");
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
    fn test_atom_form_from_str_to_string_roundtrip(#[case] s: &str) {
        let form: AtomForm = s.parse().unwrap();
        assert_eq!(form.to_string(), s);
    }

    #[rstest]
    fn test_atom_form_from_str_carbon_element() {
        let form: AtomForm = "C".parse().unwrap();
        assert_eq!(form.element, ElementForm::Lit(Element::C));
        assert_eq!(form.charge, NumForm::Undetermined);
    }

    #[rstest]
    fn test_atom_form_to_edn_roundtrip() {
        use umol_edn::ToEdn;
        let form: AtomForm = "C#c+#h3".parse().unwrap();
        let edn = form.to_edn();
        let back = AtomForm::from_edn(&edn).unwrap();
        assert_eq!(back, form);
    }
}
