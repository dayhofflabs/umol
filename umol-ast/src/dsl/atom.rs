//! Atom-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use strum::IntoEnumIterator;
use umol_edn::{DeError, Edn, EdnError, EdnStreamDeserializer, FromEdn, ToEdn};
use umol_shared::element::Element;
use winnow::ascii::multispace0;
use winnow::combinator::{alt, delimited, empty, preceded, repeat, separated, terminated};
use winnow::error::{ErrMode, ParserError};
use winnow::token::{one_of, take};
use winnow::Parser;

use super::config::{
    AromaticValenceDefault, AtomDefaults, IsotopeDefault, MulticenterValenceDefault, NumericDefault,
};
use super::error::{PResult, ParseError};
use super::predicates::{
    apply_spin_pair, charge, fmt_charge, fmt_ring_count, fmt_spin_pair, lower_spin, optional_value,
    raise_spin, ring_count, SpinPredicate,
};
use super::value::{fmt_value, id, value, ValueDsl};
use crate::ast::atom::{AtomAst, ElementAst, IsotopeAst, Polarity};
use crate::ast::constraint::{
    AromaticValenceAst, AtomConstraint, AtomConstraintKind, AtomConstraints, MulticenterValenceAst,
};
use crate::ast::traits::{FromAst, IntoAst};
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

// region: Parse

/// Parse a complete atom-string into an `AtomDsl`.
pub fn parse_atom(input: &str) -> Result<AtomDsl, ParseError> {
    if let Some(dsl) = parse_bare_element(input) {
        return Ok(dsl);
    }
    atom.parse(input).map_err(|e| e.into_inner())
}

/// Fast path for the overwhelmingly common case of an element-only
/// atom-string ("C", "N", "Cl", "Og", …). Returns `None` if the input
/// has anything beyond a single 1- or 2-byte ASCII element symbol, in
/// which case the caller falls through to the winnow parser.
fn parse_bare_element(input: &str) -> Option<AtomDsl> {
    let bytes = input.as_bytes();
    let is_first = |b: u8| b.is_ascii_uppercase();
    let is_rest = |b: u8| b.is_ascii_lowercase();
    let ok = match bytes {
        [a] => is_first(*a),
        [a, b] => is_first(*a) && is_rest(*b),
        _ => false,
    };
    if !ok {
        return None;
    }
    let el = Element::from_symbol_bytes(bytes)?;
    Some(AtomDsl(AtomAst::new(ElementAst::Lit(el))))
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

fn is_set(v: &ValueAst) -> bool {
    !matches!(v, ValueAst::Undetermined)
}

fn constraint_tag(kind: AtomConstraintKind) -> &'static str {
    match kind {
        AtomConstraintKind::Valence => "#v",
        AtomConstraintKind::TotalValence => "#V",
        AtomConstraintKind::DonatedPairs => "#d",
        AtomConstraintKind::AcceptedPairs => "#t",
        AtomConstraintKind::AromaticValence => "#a",
        AtomConstraintKind::MulticenterValence => "#m",
        AtomConstraintKind::Degree => "#D",
        AtomConstraintKind::TotalDegree => "#X",
        AtomConstraintKind::RingDegree => "#x",
        AtomConstraintKind::RingValence => "#y",
        AtomConstraintKind::TotalHydrogens => "#H",
        AtomConstraintKind::RingCount => "#R",
        AtomConstraintKind::RingSize => "#r",
    }
}

/// One predicate from an atom-string; the parser yields a `Vec` of these
/// and the applier folds them into the `AtomAst`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AtomPredicate {
    IsotopeMass(IsotopeAst),
    Charge(ValueAst),
    ImplicitHydrogens(ValueAst),
    LonePairs(ValueAst),
    Spin(SpinPredicate),
    Constraint(AtomConstraint),
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
            .map(|v| AtomPredicate::Constraint(AtomConstraint::Valence(v)))
            .parse_next(i),
        "#V" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraint::TotalValence(v)))
            .parse_next(i),
        "#d" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraint::DonatedPairs(v)))
            .parse_next(i),
        "#t" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraint::AcceptedPairs(v)))
            .parse_next(i),
        "#a" => aromatic_valence
            .map(|c| AtomPredicate::Constraint(AtomConstraint::AromaticValence(c)))
            .parse_next(i),
        "#m" => multicenter_valence
            .map(|c| AtomPredicate::Constraint(AtomConstraint::MulticenterValence(c)))
            .parse_next(i),
        "#D" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraint::Degree(v)))
            .parse_next(i),
        "#X" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraint::TotalDegree(v)))
            .parse_next(i),
        "#x" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraint::RingDegree(v)))
            .parse_next(i),
        "#y" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraint::RingValence(v)))
            .parse_next(i),
        "#H" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraint::TotalHydrogens(v)))
            .parse_next(i),
        "#R" => ring_count
            .map(|v| AtomPredicate::Constraint(AtomConstraint::RingCount(v)))
            .parse_next(i),
        "#r" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraint::RingSize(v)))
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
        preceded('!', element_set).map(ElementAst::NotSet),
        preceded('!', element_literal).map(ElementAst::Not),
        element_set.map(ElementAst::LitSet),
        element_bind.map(|(id, set, polarity)| ElementAst::Bind { id, set, polarity }),
        element_ref.map(ElementAst::Ref),
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

fn element_bind(i: &mut &str) -> PResult<(String, Vec<Element>, Polarity)> {
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

fn element_bind_domain(i: &mut &str) -> PResult<(Vec<Element>, Polarity)> {
    alt((
        preceded('!', element_set).map(|s| (s, Polarity::Exclude)),
        preceded('!', element_literal).map(|e| (vec![e], Polarity::Exclude)),
        element_set.map(|s| (s, Polarity::Include)),
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

fn isotope(i: &mut &str) -> PResult<IsotopeAst> {
    preceded(
        multispace0,
        alt(('='.value(IsotopeAst::Natural), value.map(IsotopeAst::from))),
    )
    .parse_next(i)
    .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedPredicateBody))
}

fn aromatic_valence(i: &mut &str) -> PResult<AromaticValenceAst> {
    preceded(
        multispace0,
        alt((
            "*".value(AromaticValenceAst::Undetermined),
            "!".value(AromaticValenceAst::NotAromatic),
            // `#a+` encodes "aromatic, count unspecified" — semantically the
            // same as the (a >= 0) Expr form but structurally distinct from
            // the outer Undetermined. Older sugar used Expr(Rel(Var("a"),
            // Ge, Lit(0))); the canonical form is now Aromatic(Undetermined).
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
                if !matches!(ast.isotope_mass, IsotopeAst::Undetermined) {
                    return Err(ParseError::DuplicateAtomPredicate("#i".to_string()));
                }
                ast.isotope_mass = v;
            }
            AtomPredicate::Charge(v) => {
                if is_set(&ast.charge) {
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
                if is_set(&ast.lone_pairs) {
                    return Err(ParseError::DuplicateAtomPredicate("#n".to_string()));
                }
                ast.lone_pairs = v;
            }
            AtomPredicate::Spin(sp) => {
                apply_spin_pair(&mut ast.spin, sp, ParseError::DuplicateAtomPredicate)?;
            }
            AtomPredicate::Constraint(c) => {
                if c.is_unique() && ast.constraints.contains(c.kind()) {
                    return Err(ParseError::DuplicateAtomPredicate(
                        constraint_tag(c.kind()).to_string(),
                    ));
                }
                ast.constraints.add(c);
            }
        }
    }
    Ok(())
}

// endregion: Parse

// region: Format

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
        ElementAst::LitSet(es) => fmt_element_set(f, es),
        ElementAst::Not(e) => write!(f, "!{}", e),
        ElementAst::NotSet(es) => {
            write!(f, "!")?;
            fmt_element_set(f, es)
        }
        ElementAst::Bind { id, set, polarity } => {
            write!(f, "?{} :: ", id)?;
            if matches!(polarity, Polarity::Exclude) {
                write!(f, "!")?;
            }
            fmt_element_set(f, set)
        }
        ElementAst::Ref(id) => write!(f, "?{}", id),
    }
}

fn fmt_element_set(f: &mut fmt::Formatter<'_>, es: &[Element]) -> fmt::Result {
    write!(f, "{{")?;
    for (i, e) in es.iter().enumerate() {
        if i > 0 {
            write!(f, ",")?;
        }
        write!(f, "{}", e)?;
    }
    write!(f, "}}")
}

fn fmt_isotope_mass(f: &mut fmt::Formatter<'_>, iso: &IsotopeAst) -> fmt::Result {
    match iso {
        IsotopeAst::Undetermined => Ok(()),
        IsotopeAst::Natural => write!(f, "#i="),
        IsotopeAst::Lit(n) => write!(f, "#i{}", n),
        IsotopeAst::LitSet(s) => {
            write!(f, "#i")?;
            fmt_value(f, &ValueAst::LitSet(s.clone()))
        }
        IsotopeAst::Expr(e) => {
            write!(f, "#i")?;
            fmt_value(f, &ValueAst::Expr(e.clone()))
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

fn fmt_constraint(f: &mut fmt::Formatter<'_>, c: &AtomConstraint) -> fmt::Result {
    match c {
        AtomConstraint::Valence(v) => fmt_value_field_required(f, "#v", v),
        AtomConstraint::DonatedPairs(v) => fmt_value_field_required(f, "#d", v),
        AtomConstraint::AcceptedPairs(v) => fmt_value_field_required(f, "#t", v),
        AtomConstraint::MulticenterValence(c) => match c {
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
        AtomConstraint::AromaticValence(c) => match c {
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
        AtomConstraint::Degree(v) => fmt_value_field_required(f, "#D", v),
        AtomConstraint::TotalDegree(v) => fmt_value_field_required(f, "#X", v),
        AtomConstraint::RingDegree(v) => fmt_value_field_required(f, "#x", v),
        AtomConstraint::RingValence(v) => fmt_value_field_required(f, "#y", v),
        AtomConstraint::TotalValence(v) => fmt_value_field_required(f, "#V", v),
        AtomConstraint::TotalHydrogens(v) => fmt_value_field_required(f, "#H", v),
        AtomConstraint::RingCount(v) => fmt_ring_count(f, v),
        AtomConstraint::RingSize(v) => fmt_value_field_required(f, "#r", v),
    }
}

// endregion: Format

// region: Raise

fn raise_atom(ast: &mut AtomAst, cfg: &AtomDefaults) {
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

    if matches!(*isotope_mass, IsotopeAst::Undetermined) {
        *isotope_mass = match cfg.isotope {
            IsotopeDefault::Natural => IsotopeAst::Natural,
            IsotopeDefault::Required => IsotopeAst::Undetermined,
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

fn raise_atom_constraints(constraints: &mut AtomConstraints, cfg: &AtomDefaults) {
    constraints.retain(|c| !c.is_undetermined());

    // Exhaustive dispatch over every kind: a new AtomConstraintKind variant
    // fails to build here until it has an explicit branch.
    for kind in AtomConstraintKind::iter() {
        match kind {
            AtomConstraintKind::Valence => {
                if matches!(cfg.valence, NumericDefault::Zero) && !constraints.contains(kind) {
                    constraints.add(AtomConstraint::Valence(ValueAst::Lit(0)));
                }
            }
            AtomConstraintKind::DonatedPairs => {
                if matches!(cfg.donated_pairs, NumericDefault::Zero) && !constraints.contains(kind)
                {
                    constraints.add(AtomConstraint::DonatedPairs(ValueAst::Lit(0)));
                }
            }
            AtomConstraintKind::AcceptedPairs => {
                if matches!(cfg.accepted_pairs, NumericDefault::Zero) && !constraints.contains(kind)
                {
                    constraints.add(AtomConstraint::AcceptedPairs(ValueAst::Lit(0)));
                }
            }
            AtomConstraintKind::AromaticValence => {
                if !constraints.contains(kind) {
                    match cfg.aromatic_valence {
                        AromaticValenceDefault::NotAromatic => {
                            constraints.add(AtomConstraint::AromaticValence(
                                AromaticValenceAst::NotAromatic,
                            ));
                        }
                        AromaticValenceDefault::Aromatic => {
                            constraints.add(AtomConstraint::AromaticValence(
                                AromaticValenceAst::Aromatic(ValueAst::Undetermined),
                            ));
                        }
                        AromaticValenceDefault::Required => {}
                    }
                }
            }
            AtomConstraintKind::MulticenterValence => {
                if !constraints.contains(kind) {
                    match cfg.multicenter_valence {
                        MulticenterValenceDefault::NotMulticenter => {
                            constraints.add(AtomConstraint::MulticenterValence(
                                MulticenterValenceAst::NotMulticenter,
                            ));
                        }
                        MulticenterValenceDefault::Multicenter => {
                            constraints.add(AtomConstraint::MulticenterValence(
                                MulticenterValenceAst::Multicenter(ValueAst::Undetermined),
                            ));
                        }
                        MulticenterValenceDefault::Required => {}
                    }
                }
            }
            AtomConstraintKind::TotalValence
            | AtomConstraintKind::Degree
            | AtomConstraintKind::TotalDegree
            | AtomConstraintKind::RingDegree
            | AtomConstraintKind::RingValence
            | AtomConstraintKind::TotalHydrogens
            | AtomConstraintKind::RingCount
            | AtomConstraintKind::RingSize => {
                // Pattern-only constraint: no defaulting mode in AtomDefaults.
            }
        }
    }
}

// endregion: Raise

// region: Lower

fn lower_atom(ast: &mut AtomAst, cfg: &AtomDefaults) {
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
        (IsotopeDefault::Natural, IsotopeAst::Natural)
    ) {
        *isotope_mass = IsotopeAst::Undetermined;
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

fn lower_atom_constraints(constraints: &mut AtomConstraints, cfg: &AtomDefaults) {
    // Exhaustive dispatch over every kind: a new AtomConstraintKind variant
    // fails to build here until it has an explicit branch.
    for kind in AtomConstraintKind::iter() {
        match kind {
            AtomConstraintKind::Valence => {
                if matches!(cfg.valence, NumericDefault::Zero)
                    && matches!(
                        constraints.get(kind),
                        Some(AtomConstraint::Valence(ValueAst::Lit(0)))
                    )
                {
                    constraints.remove(kind);
                }
            }
            AtomConstraintKind::DonatedPairs => {
                if matches!(cfg.donated_pairs, NumericDefault::Zero)
                    && matches!(
                        constraints.get(kind),
                        Some(AtomConstraint::DonatedPairs(ValueAst::Lit(0)))
                    )
                {
                    constraints.remove(kind);
                }
            }
            AtomConstraintKind::AcceptedPairs => {
                if matches!(cfg.accepted_pairs, NumericDefault::Zero)
                    && matches!(
                        constraints.get(kind),
                        Some(AtomConstraint::AcceptedPairs(ValueAst::Lit(0)))
                    )
                {
                    constraints.remove(kind);
                }
            }
            AtomConstraintKind::MulticenterValence => match cfg.multicenter_valence {
                MulticenterValenceDefault::NotMulticenter => {
                    if matches!(
                        constraints.get(kind),
                        Some(AtomConstraint::MulticenterValence(
                            MulticenterValenceAst::NotMulticenter
                        ))
                    ) {
                        constraints.remove(kind);
                    }
                }
                MulticenterValenceDefault::Multicenter => {
                    if matches!(
                        constraints.get(kind),
                        Some(AtomConstraint::MulticenterValence(
                            MulticenterValenceAst::Multicenter(ValueAst::Undetermined)
                        ))
                    ) {
                        constraints.remove(kind);
                    }
                }
                MulticenterValenceDefault::Required => {}
            },
            AtomConstraintKind::AromaticValence => match cfg.aromatic_valence {
                AromaticValenceDefault::NotAromatic => {
                    if matches!(
                        constraints.get(kind),
                        Some(AtomConstraint::AromaticValence(
                            AromaticValenceAst::NotAromatic
                        ))
                    ) {
                        constraints.remove(kind);
                    }
                }
                AromaticValenceDefault::Aromatic => {
                    if matches!(
                        constraints.get(kind),
                        Some(AtomConstraint::AromaticValence(
                            AromaticValenceAst::Aromatic(ValueAst::Undetermined)
                        ))
                    ) {
                        constraints.remove(kind);
                    }
                }
                AromaticValenceDefault::Required => {}
            },
            AtomConstraintKind::TotalValence
            | AtomConstraintKind::Degree
            | AtomConstraintKind::TotalDegree
            | AtomConstraintKind::RingDegree
            | AtomConstraintKind::RingValence
            | AtomConstraintKind::TotalHydrogens
            | AtomConstraintKind::RingCount
            | AtomConstraintKind::RingSize => {
                // Pattern-only constraint: no defaulting mode in AtomDefaults.
            }
        }
    }
}

// endregion: Lower

// region: Constraint DSLs

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
                Edn::Keyword(umol_edn::EdnKeyword::owned("undetermined".into()))
            }
            AromaticValenceAst::NotAromatic => {
                Edn::Keyword(umol_edn::EdnKeyword::owned("not-aromatic".into()))
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
                Edn::Keyword(umol_edn::EdnKeyword::owned("undetermined".into()))
            }
            MulticenterValenceAst::NotMulticenter => {
                Edn::Keyword(umol_edn::EdnKeyword::owned("not-multicenter".into()))
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

/// Surface DSL wrapper around `AtomConstraint`. EDN form is a single-key map
/// keyed by the constraint kind: e.g. `{:valence 4}`, `{:degree *}`,
/// `{:aromatic-valence :not-aromatic}`, `{:total-hydrogens {?h :: {0,1,2}}}`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct AtomConstraintDsl(pub AtomConstraint);

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
            "valence" => AtomConstraint::Valence(ValueDsl::from_edn(v)?.into_ast(&())),
            "total-valence" => AtomConstraint::TotalValence(ValueDsl::from_edn(v)?.into_ast(&())),
            "aromatic-valence" => {
                AtomConstraint::AromaticValence(AromaticValenceDsl::from_edn(v)?.into_ast(&()))
            }
            "multicenter-valence" => AtomConstraint::MulticenterValence(
                MulticenterValenceDsl::from_edn(v)?.into_ast(&()),
            ),
            "donated-pairs" => AtomConstraint::DonatedPairs(ValueDsl::from_edn(v)?.into_ast(&())),
            "accepted-pairs" => AtomConstraint::AcceptedPairs(ValueDsl::from_edn(v)?.into_ast(&())),
            "degree" => AtomConstraint::Degree(ValueDsl::from_edn(v)?.into_ast(&())),
            "total-degree" => AtomConstraint::TotalDegree(ValueDsl::from_edn(v)?.into_ast(&())),
            "ring-degree" => AtomConstraint::RingDegree(ValueDsl::from_edn(v)?.into_ast(&())),
            "ring-valence" => AtomConstraint::RingValence(ValueDsl::from_edn(v)?.into_ast(&())),
            "total-hydrogens" => {
                AtomConstraint::TotalHydrogens(ValueDsl::from_edn(v)?.into_ast(&()))
            }
            "ring-count" => AtomConstraint::RingCount(ValueDsl::from_edn(v)?.into_ast(&())),
            "ring-size" => AtomConstraint::RingSize(ValueDsl::from_edn(v)?.into_ast(&())),
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
            AtomConstraint::Valence(v) => {
                single_key_map("valence", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraint::TotalValence(v) => {
                single_key_map("total-valence", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraint::AromaticValence(c) => single_key_map(
                "aromatic-valence",
                AromaticValenceDsl::from_ast(c, &()).to_edn(),
            ),
            AtomConstraint::MulticenterValence(c) => single_key_map(
                "multicenter-valence",
                MulticenterValenceDsl::from_ast(c, &()).to_edn(),
            ),
            AtomConstraint::DonatedPairs(v) => {
                single_key_map("donated-pairs", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraint::AcceptedPairs(v) => {
                single_key_map("accepted-pairs", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraint::Degree(v) => {
                single_key_map("degree", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraint::TotalDegree(v) => {
                single_key_map("total-degree", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraint::RingDegree(v) => {
                single_key_map("ring-degree", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraint::RingValence(v) => {
                single_key_map("ring-valence", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraint::TotalHydrogens(v) => {
                single_key_map("total-hydrogens", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraint::RingCount(v) => {
                single_key_map("ring-count", ValueDsl::from_ast(v, &()).to_edn())
            }
            AtomConstraint::RingSize(v) => {
                single_key_map("ring-size", ValueDsl::from_ast(v, &()).to_edn())
            }
        }
    }
}

impl FromAst<AtomConstraint> for AtomConstraintDsl {
    type Ctx = ();

    fn from_ast(ast: &AtomConstraint, _ctx: &Self::Ctx) -> Self {
        Self(ast.clone())
    }
}

impl IntoAst<AtomConstraint> for AtomConstraintDsl {
    type Ctx = ();

    fn into_ast(self, _ctx: &Self::Ctx) -> AtomConstraint {
        self.0
    }
}

fn single_key_map(key: &str, value: Edn<'static>) -> Edn<'static> {
    let mut m = umol_edn::EdnMap::with_capacity(1);
    m.insert(Edn::Keyword(umol_edn::EdnKeyword::owned(key.into())), value);
    Edn::Map(m)
}

// endregion: Constraint DSLs

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::element::Element;

    use super::*;
    use crate::ast::spin::SpinStateAst;
    use crate::ast::value::{ArithOp, Expr, RelOp};

    #[rstest]
    #[case::single("C", AtomDsl(AtomAst::new(ElementAst::Lit(Element::C))))]
    #[case::double("Cl", AtomDsl(AtomAst::new(ElementAst::Lit(Element::Cl))))]
    #[case::transuranic("Og", AtomDsl(AtomAst::new(ElementAst::Lit(Element::Og))))]
    fn test_parse_bare_element_accepts(#[case] input: &str, #[case] expected: AtomDsl) {
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
    #[case::bind("(?e :: {C,N})")]
    #[case::ref_("(?e)")]
    #[case::h_count("C#h3")]
    #[case::charge_plus("N#c+")]
    #[case::full("C#c+1#R+#v4")]
    fn test_parse_bare_element_rejects(#[case] input: &str) {
        assert_eq!(parse_bare_element(input), None);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::carbon("C", AtomDsl(AtomAst::new(ElementAst::Lit(Element::C))))]
    #[case::iron("Fe", AtomDsl(AtomAst::new(ElementAst::Lit(Element::Fe))))]
    #[case::chlorine("Cl", AtomDsl(AtomAst::new(ElementAst::Lit(Element::Cl))))]
    #[case::whitespace("  C  ", AtomDsl(AtomAst::new(ElementAst::Lit(Element::C))))]
    #[case::undetermined("*", AtomDsl(AtomAst::new(ElementAst::Undetermined)))]
    #[case::element_set("{C,N,O}", AtomDsl(AtomAst::new(ElementAst::LitSet(vec![Element::C, Element::N, Element::O]))))]
    #[case::element_bind("(?e :: {C,N})", AtomDsl(AtomAst::new(ElementAst::Bind { id: "e".to_string(), set: vec![Element::C, Element::N], polarity: Polarity::Include })))]
    #[case::element_ref("(?e)", AtomDsl(AtomAst::new(ElementAst::Ref("e".to_string()))))]
    #[case::isotope("C#i12", AtomDsl(AtomAst { isotope_mass: IsotopeAst::Lit(12), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::isotope_natural("C#i=", AtomDsl(AtomAst { isotope_mass: IsotopeAst::Natural, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_pos("C#c+2", AtomDsl(AtomAst { charge: ValueAst::Lit(2), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_neg("C#c-2", AtomDsl(AtomAst { charge: ValueAst::Lit(-2), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_plus("C#c+", AtomDsl(AtomAst { charge: ValueAst::Lit(1), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_minus("C#c-", AtomDsl(AtomAst { charge: ValueAst::Lit(-1), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_zero("C#c0", AtomDsl(AtomAst { charge: ValueAst::Lit(0), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_count("C#h3", AtomDsl(AtomAst { implicit_hydrogens: ValueAst::Lit(3), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_undetermined("C#h*", AtomDsl(AtomAst { implicit_hydrogens: ValueAst::Undetermined, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_bind("C#h(?h)", AtomDsl(AtomAst { implicit_hydrogens: ValueAst::Ref("h".to_string()), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_set("N#h?h :: {2,3}", AtomDsl(AtomAst { implicit_hydrogens: ValueAst::bind("h", vec![2, 3]), ..AtomAst::new(ElementAst::Lit(Element::N)) }))]
    #[case::h_expr("C#h?h >= 1", AtomDsl(AtomAst { implicit_hydrogens: ValueAst::Expr(Box::new(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1))))), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_omit("C#h", AtomDsl(AtomAst { implicit_hydrogens: ValueAst::Lit(1), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::lone_pairs("O#n2", AtomDsl(AtomAst { lone_pairs: ValueAst::Lit(2), ..AtomAst::new(ElementAst::Lit(Element::O)) }))]
    #[case::lone_pairs_omit("O#n", AtomDsl(AtomAst { lone_pairs: ValueAst::Lit(1), ..AtomAst::new(ElementAst::Lit(Element::O)) }))]
    #[case::unpaired("C#u2", AtomDsl(AtomAst { spin: SpinStateAst { unpaired: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined }, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::unpaired_omit("C#u", AtomDsl(AtomAst { spin: SpinStateAst { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Undetermined }, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multiplicity("C#s3", AtomDsl(AtomAst { spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(3) }, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multiplicity_omit("C#s", AtomDsl(AtomAst { spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(1) }, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::valence("C#v4", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::Valence(ValueAst::Lit(4))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::donated_pairs("N#d1", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::DonatedPairs(ValueAst::Lit(1))]), ..AtomAst::new(ElementAst::Lit(Element::N)) }))]
    #[case::accepted_pairs("B#t1", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::AcceptedPairs(ValueAst::Lit(1))]), ..AtomAst::new(ElementAst::Lit(Element::B)) }))]
    #[case::ring_size("C#r6", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::RingSize(ValueAst::Lit(6))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::arom_not_aromatic("C#a!", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::AromaticValence(AromaticValenceAst::NotAromatic)]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::arom_undetermined("C#a*", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::AromaticValence(AromaticValenceAst::Undetermined)]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::arom_plus("C#a+", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Undetermined))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::arom_zero("C#a0", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(0)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::arom_one("C#a1", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(1)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::arom_omit("C#a", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(1)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multicenter_not("C#m!", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::MulticenterValence(MulticenterValenceAst::NotMulticenter)]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multicenter_undetermined("C#m*", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::MulticenterValence(MulticenterValenceAst::Undetermined)]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multicenter_plus("C#m+", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Undetermined))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multicenter_zero("C#m0", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(0)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multicenter_one("C#m", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(1)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multicenter("C#m2", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(2)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::degree("C#D2", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::Degree(ValueAst::Lit(2))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::total_degree("C#X3", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::TotalDegree(ValueAst::Lit(3))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_degree("C#x2", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::RingDegree(ValueAst::Lit(2))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_valence("C#y3", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::RingValence(ValueAst::Lit(3))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::total_valence("C#V5", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::TotalValence(ValueAst::Lit(5))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::total_hydrogens("C#H1", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::TotalHydrogens(ValueAst::Lit(1))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_bare("C#R", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::RingCount(ValueAst::Lit(1))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_undetermined("C#R*", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::RingCount(ValueAst::Undetermined)]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_plus("C#R+", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::RingCount(ValueAst::Expr(Box::new(Expr::Rel(Box::new(Expr::Var("r".to_string())), RelOp::Ge, Box::new(Expr::Lit(1))))))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_bang("C#R!", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::RingCount(ValueAst::Lit(0))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_zero("C#R0", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::RingCount(ValueAst::Lit(0))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_count("C#R2", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::RingCount(ValueAst::Lit(2))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_size_conj("C#r5#r6", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::RingSize(ValueAst::Lit(5)), AtomConstraint::RingSize(ValueAst::Lit(6))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
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
    #[case::dup_h("C#h3#h2", ParseError::DuplicateAtomPredicate("#h".to_string()))]
    #[case::dup_charge("C#c+#c-", ParseError::DuplicateAtomPredicate("#c".to_string()))]
    #[case::dup_valence("C#v3#v4", ParseError::DuplicateAtomPredicate("#v".to_string()))]
    #[case::invalid_special_slash("C#h/", ParseError::TrailingInput("/".to_string()))]
    #[case::invalid_special_minus("C#h-", ParseError::TrailingInput("-".to_string()))]
    #[case::invalid_special_equal("C#h=", ParseError::TrailingInput("=".to_string()))]
    #[case::trailing("C#h3 foo", ParseError::TrailingInput("foo".to_string()))]
    fn test_parse_atom_invalid(#[case] input: &str, #[case] expected: ParseError) {
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
    #[case::arom_not_aromatic("C#a!")]
    #[case::arom_plus("C#a+")]
    #[case::arom_zero("C#a0")]
    #[case::arom_omit("C#a")]
    #[case::multicenter_not("C#m!")]
    #[case::multicenter_plus("C#m+")]
    #[case::multicenter_zero("C#m0")]
    #[case::multicenter_omit("C#m")]
    #[case::ring_bare("C#R")]
    #[case::ring_plus("C#R+")]
    #[case::ring_bang("C#R!")]
    #[case::ring_count("C#R2")]
    #[case::ring_size_conj("C#r5#r6")]
    fn test_atom_display_roundtrip(#[case] input: &str) {
        let parsed = atom.parse(input).unwrap();
        assert_eq!(parsed.to_string(), input);
    }

    /// Vacuous constraints (those with `Undetermined` payload) elide on
    /// rendering per the canonical-rendering rule (see `dsl::predicates`).
    /// `#v*`, `#R*`, `#m*`, `#a*` etc. are admitted on parse but the
    /// rendered surface drops them entirely, so the constraint is gone
    /// after a render → reparse cycle.
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
    #[case::ring_count("C#R*", "C")]
    #[case::ring_size("C#r*", "C")]
    #[case::aromatic_undetermined("C#a*", "C")]
    #[case::multicenter_undetermined("C#m*", "C")]
    fn test_atom_render_elides_vacuous_constraints(
        #[case] input: &str,
        #[case] expected_canonical: &str,
    ) {
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
    #[case::set("{C,N,O}", ElementAst::LitSet(vec![Element::C, Element::N, Element::O]))]
    #[case::set_spaced("{ C, N}", ElementAst::LitSet(vec![Element::C, Element::N]))]
    #[case::bind("(?e :: {C,N})", ElementAst::Bind { id: "e".to_string(), set: vec![Element::C, Element::N], polarity: Polarity::Include })]
    #[case::ref_("(?e)", ElementAst::Ref("e".to_string()))]
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
    fn test_element_invalid(#[case] input: &str) {
        let result = element.parse(input);
        assert!(
            result.is_err(),
            "{input:?} should fail, got {:?}",
            result.unwrap()
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::natural("=", IsotopeAst::Natural)]
    #[case::lit("12", IsotopeAst::Lit(12))]
    #[case::undetermined("*", IsotopeAst::Undetermined)]
    #[case::set("{12,13,14}", IsotopeAst::LitSet(Box::new(vec![12, 13, 14])))]
    #[case::bind("(?m :: {12,13})", IsotopeAst::Expr(Box::new(Expr::Mem(Box::new(Expr::Var("m".to_string())), vec![12, 13]))))]
    #[case::ref_("(?m)", IsotopeAst::Expr(Box::new(Expr::Var("m".to_string()))))]
    fn test_isotope(#[case] input: &str, #[case] expected: IsotopeAst) {
        let result = isotope.parse(input);
        assert!(result.is_ok(), "{input:?} should succeed, got {:?}", result.unwrap_err());
        let expr = result.unwrap();
        assert_eq!(expr, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::isotope_lit("#i12", AtomPredicate::IsotopeMass(IsotopeAst::Lit(12)))]
    #[case::isotope_natural("#i=", AtomPredicate::IsotopeMass(IsotopeAst::Natural))]
    #[case::isotope_undetermined("#i*", AtomPredicate::IsotopeMass(IsotopeAst::Undetermined))]
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
    #[case::valence("#v4", AtomPredicate::Constraint(AtomConstraint::Valence(ValueAst::Lit(4))))]
    #[case::total_valence("#V5", AtomPredicate::Constraint(AtomConstraint::TotalValence(ValueAst::Lit(5))))]
    #[case::total_valence_omit("#V", AtomPredicate::Constraint(AtomConstraint::TotalValence(ValueAst::Lit(1))))]
    #[case::ring_valence("#y2", AtomPredicate::Constraint(AtomConstraint::RingValence(ValueAst::Lit(2))))]
    #[case::ring_valence_omit("#y", AtomPredicate::Constraint(AtomConstraint::RingValence(ValueAst::Lit(1))))]
    #[case::donated_pairs("#d1", AtomPredicate::Constraint(AtomConstraint::DonatedPairs(ValueAst::Lit(1))))]
    #[case::accepted_pairs("#t1", AtomPredicate::Constraint(AtomConstraint::AcceptedPairs(ValueAst::Lit(1))))]
    #[case::ring_size("#r6", AtomPredicate::Constraint(AtomConstraint::RingSize(ValueAst::Lit(6))))]
    #[case::arom_not_aromatic("#a!", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceAst::NotAromatic)))]
    #[case::arom_undetermined("#a*", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceAst::Undetermined)))]
    #[case::arom_plus("#a+", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Undetermined))))]
    #[case::arom_lit("#a2", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(2)))))]
    #[case::arom_omit("#a", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(1)))))]
    #[case::multicenter_not("#m!", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(MulticenterValenceAst::NotMulticenter)))]
    #[case::multicenter_undetermined("#m*", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(MulticenterValenceAst::Undetermined)))]
    #[case::multicenter_plus("#m+", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Undetermined))))]
    #[case::multicenter_zero("#m0", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(0)))))]
    #[case::multicenter_omit("#m", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(1)))))]
    #[case::multicenter("#m2", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(2)))))]
    #[case::degree("#D2", AtomPredicate::Constraint(AtomConstraint::Degree(ValueAst::Lit(2))))]
    #[case::degree_omit("#D", AtomPredicate::Constraint(AtomConstraint::Degree(ValueAst::Lit(1))))]
    #[case::total_degree("#X3", AtomPredicate::Constraint(AtomConstraint::TotalDegree(ValueAst::Lit(3))))]
    #[case::ring_degree("#x2", AtomPredicate::Constraint(AtomConstraint::RingDegree(ValueAst::Lit(2))))]
    #[case::ring_degree_omit("#x", AtomPredicate::Constraint(AtomConstraint::RingDegree(ValueAst::Lit(1))))]
    #[case::total_hydrogens("#H1", AtomPredicate::Constraint(AtomConstraint::TotalHydrogens(ValueAst::Lit(1))))]
    #[case::ring_bare("#R", AtomPredicate::Constraint(AtomConstraint::RingCount(ValueAst::Lit(1))))]
    #[case::ring_undetermined("#R*", AtomPredicate::Constraint(AtomConstraint::RingCount(ValueAst::Undetermined)))]
    #[case::ring_plus("#R+", AtomPredicate::Constraint(AtomConstraint::RingCount(ValueAst::Expr(Box::new(Expr::Rel(Box::new(Expr::Var("r".to_string())), RelOp::Ge, Box::new(Expr::Lit(1))))))))]
    #[case::ring_count("#R2", AtomPredicate::Constraint(AtomConstraint::RingCount(ValueAst::Lit(2))))]
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
    fn test_atom_dsl_to_ast_fills_zero_defaults() {
        let dsl = AtomDsl(AtomAst::new(ElementAst::Lit(Element::C)));
        let cfg = AtomDefaults::zeroed();
        let ast = dsl.into_ast(&cfg);
        assert_eq!(ast.charge, ValueAst::Lit(0));
        assert_eq!(ast.lone_pairs, ValueAst::Lit(0));
        assert_eq!(ast.implicit_hydrogens, ValueAst::Lit(0));
        assert_eq!(ast.isotope_mass, IsotopeAst::Natural);
        assert_eq!(ast.spin, SpinStateAst::from((0_u8, 1_u8)));
        assert_eq!(
            ast.constraints.get(AtomConstraintKind::Valence),
            Some(&AtomConstraint::Valence(ValueAst::Lit(0)))
        );
        assert_eq!(
            ast.constraints.get(AtomConstraintKind::AromaticValence),
            Some(&AtomConstraint::AromaticValence(
                AromaticValenceAst::NotAromatic
            ))
        );
    }

    #[rstest]
    fn test_atom_dsl_from_ast_strips_zero_defaults() {
        let mut ast = AtomAst::new(ElementAst::Lit(Element::C));
        ast.charge = ValueAst::Lit(0);
        ast.lone_pairs = ValueAst::Lit(0);
        ast.implicit_hydrogens = ValueAst::Lit(0);
        ast.isotope_mass = IsotopeAst::Natural;
        ast.spin = SpinStateAst::from((0_u8, 1_u8));
        ast.constraints
            .add(AtomConstraint::Valence(ValueAst::Lit(0)));
        ast.constraints.add(AtomConstraint::AromaticValence(
            AromaticValenceAst::NotAromatic,
        ));
        let cfg = AtomDefaults::zeroed();
        let dsl = AtomDsl::from_ast(&ast, &cfg);
        assert_eq!(dsl.0.charge, ValueAst::Undetermined);
        assert_eq!(dsl.0.lone_pairs, ValueAst::Undetermined);
        assert_eq!(dsl.0.implicit_hydrogens, ValueAst::Undetermined);
        assert_eq!(dsl.0.isotope_mass, IsotopeAst::Undetermined);
        assert_eq!(dsl.0.spin, SpinStateAst::default());
        assert!(dsl.0.constraints.is_empty());
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
    fn test_atom_dsl_from_edn_str_matches_from_edn(#[case] input: &str) {
        let via_stream = AtomDsl::from_edn_str(input).unwrap();
        let tree = umol_edn::read_string(input).unwrap();
        let via_tree = AtomDsl::from_edn(&tree).unwrap();
        assert_eq!(via_stream, via_tree);
    }

    // region: AromaticValenceDsl / MulticenterValenceDsl

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(AromaticValenceAst::Undetermined, ":undetermined")]
    #[case::not_aromatic(AromaticValenceAst::NotAromatic, ":not-aromatic")]
    #[case::aromatic_lit(AromaticValenceAst::Aromatic(ValueAst::Lit(6)), "{:aromatic 6}")]
    #[case::aromatic_undetermined(AromaticValenceAst::Aromatic(ValueAst::Undetermined), "{:aromatic :undetermined}")]
    fn test_aromatic_valence_dsl_roundtrip(
        #[case] input: AromaticValenceAst,
        #[case] edn_source: &str,
    ) {
        let dsl = AromaticValenceDsl::from_ast(&input, &());
        let edn = dsl.to_edn();
        let expected = umol_edn::read_string(edn_source).unwrap();
        assert_eq!(edn, expected);
        let parsed = AromaticValenceDsl::from_edn(&edn).unwrap();
        assert_eq!(parsed.into_ast(&()), input);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::undetermined(MulticenterValenceAst::Undetermined, ":undetermined")]
    #[case::not_multicenter(MulticenterValenceAst::NotMulticenter, ":not-multicenter")]
    #[case::multicenter_lit(MulticenterValenceAst::Multicenter(ValueAst::Lit(3)), "{:multicenter 3}")]
    fn test_multicenter_valence_dsl_roundtrip(
        #[case] input: MulticenterValenceAst,
        #[case] edn_source: &str,
    ) {
        let dsl = MulticenterValenceDsl::from_ast(&input, &());
        let edn = dsl.to_edn();
        let expected = umol_edn::read_string(edn_source).unwrap();
        assert_eq!(edn, expected);
        let parsed = MulticenterValenceDsl::from_edn(&edn).unwrap();
        assert_eq!(parsed.into_ast(&()), input);
    }

    #[rstest]
    fn test_aromatic_valence_dsl_rejects_unknown_key() {
        let edn = umol_edn::read_string("{:bogus 1}").unwrap();
        let err = AromaticValenceDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::UnknownField { .. }));
    }

    #[rstest]
    fn test_multicenter_valence_dsl_rejects_wrong_shape() {
        let err = MulticenterValenceDsl::from_edn(&Edn::Int(3)).unwrap_err();
        assert!(matches!(err, DeError::TypeMismatch { .. }));
    }

    // endregion: AromaticValenceDsl / MulticenterValenceDsl

    // region: AtomConstraintDsl

    #[rustfmt::skip]
    #[rstest]
    #[case::valence(AtomConstraint::Valence(ValueAst::Lit(4)), "{:valence 4}")]
    #[case::valence_undetermined(AtomConstraint::Valence(ValueAst::Undetermined), "{:valence :undetermined}")]
    #[case::valence_set(AtomConstraint::Valence(ValueAst::LitSet(Box::new(vec![3, 4]))), "{:valence [3 4]}")]
    #[case::degree(AtomConstraint::Degree(ValueAst::Lit(3)), "{:degree 3}")]
    #[case::total_degree(AtomConstraint::TotalDegree(ValueAst::Lit(4)), "{:total-degree 4}")]
    #[case::ring_degree(AtomConstraint::RingDegree(ValueAst::Lit(2)), "{:ring-degree 2}")]
    #[case::ring_valence(AtomConstraint::RingValence(ValueAst::Lit(3)), "{:ring-valence 3}")]
    #[case::total_valence(AtomConstraint::TotalValence(ValueAst::Lit(5)), "{:total-valence 5}")]
    #[case::total_h(AtomConstraint::TotalHydrogens(ValueAst::Lit(3)), "{:total-hydrogens 3}")]
    #[case::ring_count(AtomConstraint::RingCount(ValueAst::Lit(1)), "{:ring-count 1}")]
    #[case::ring_size(AtomConstraint::RingSize(ValueAst::Lit(6)), "{:ring-size 6}")]
    #[case::donated(AtomConstraint::DonatedPairs(ValueAst::Lit(1)), "{:donated-pairs 1}")]
    #[case::accepted(AtomConstraint::AcceptedPairs(ValueAst::Lit(2)), "{:accepted-pairs 2}")]
    #[case::aromatic_not(AtomConstraint::AromaticValence(AromaticValenceAst::NotAromatic), "{:aromatic-valence :not-aromatic}")]
    #[case::aromatic_value(AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(6))), "{:aromatic-valence {:aromatic 6}}")]
    #[case::multicenter_not(AtomConstraint::MulticenterValence(MulticenterValenceAst::NotMulticenter), "{:multicenter-valence :not-multicenter}")]
    #[case::multicenter_value(AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(3))), "{:multicenter-valence {:multicenter 3}}")]
    #[case::valence_expr(AtomConstraint::Valence(ValueAst::Expr(Box::new(Expr::Rel(Box::new(Expr::Var("h".into())), RelOp::Ge, Box::new(Expr::Lit(1)))))), "{:valence \"?h >= 1\"}")]
    #[case::ring_size_litset(AtomConstraint::RingSize(ValueAst::LitSet(Box::new(vec![5, 6]))), "{:ring-size [5 6]}")]
    fn test_atom_constraint_dsl_roundtrip(#[case] input: AtomConstraint, #[case] edn_source: &str) {
        let dsl = AtomConstraintDsl::from_ast(&input, &());
        let edn = dsl.to_edn();
        let expected = umol_edn::read_string(edn_source).unwrap();
        assert_eq!(edn, expected, "render mismatch");
        let parsed = AtomConstraintDsl::from_edn(&edn).unwrap();
        assert_eq!(parsed.into_ast(&()), input, "parse-back mismatch");
    }

    #[rstest]
    fn test_atom_constraint_dsl_rejects_non_map() {
        let err = AtomConstraintDsl::from_edn(&Edn::Int(3)).unwrap_err();
        assert!(matches!(err, DeError::TypeMismatch { .. }));
    }

    #[rstest]
    fn test_atom_constraint_dsl_rejects_multiple_keys() {
        let edn = umol_edn::read_string("{:valence 4 :degree 3}").unwrap();
        let err = AtomConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::Custom(_)));
    }

    #[rstest]
    fn test_atom_constraint_dsl_rejects_unknown_key() {
        let edn = umol_edn::read_string("{:bogus 1}").unwrap();
        let err = AtomConstraintDsl::from_edn(&edn).unwrap_err();
        assert!(matches!(err, DeError::UnknownField { .. }));
    }

    #[rstest]
    fn test_atom_constraint_dsl_accepts_value_as_string_subgrammar() {
        let edn = umol_edn::read_string(r##"{:valence "?h + 1"}"##).unwrap();
        let parsed = AtomConstraintDsl::from_edn(&edn).unwrap();
        assert_eq!(
            parsed.into_ast(&()),
            AtomConstraint::Valence(ValueAst::Expr(Box::new(Expr::BinOp(
                Box::new(Expr::Var("h".into())),
                ArithOp::Add,
                Box::new(Expr::Lit(1)),
            ))))
        );
    }
    // endregion: AtomConstraintDsl

    // region: AtomAst symmetric I/O

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

    // endregion: AtomAst symmetric I/O
}
