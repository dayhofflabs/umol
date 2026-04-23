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

use super::error::{PResult, ParseError};
use super::predicates::{
    apply_spin_pair, charge, fmt_charge, fmt_ring_count, fmt_spin_pair, fmt_value, is_plus_sugar,
    lower_spin, optional_value, raise_spin, ring_count, SpinPredicate,
};
use super::value::{id, value};
use crate::ast::atom::{AtomAst, ElementAst, ImplicitHydrogensAst, IsotopeAst};
use crate::ast::config::{
    AromaticValenceMode, AtomAstConfig, ImplicitHydrogenMode, IsotopeMode, MulticenterValenceMode,
    NumericMode,
};
use crate::ast::constraint::{
    AromaticValenceAst, AtomConstraint, AtomConstraintKind, AtomConstraints, MulticenterValenceAst,
};
use crate::ast::traits::{FromAst, ToAst};
use crate::ast::value::{Expr, RelOp, ValueAst};

/// Surface DSL wrapper around `AtomAst`. Parses and renders the atom-string form
/// (element plus `#…` predicates); inline-capable constraints land in
/// `self.0.constraints`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AtomDsl(pub AtomAst);

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
    type Error = ParseError;

    fn from_ast(ast: &AtomAst, cfg: &AtomAstConfig) -> Result<Self, ParseError> {
        let mut out = ast.clone();
        lower_atom(&mut out, cfg);
        Ok(AtomDsl(out))
    }
}

impl ToAst<AtomAst> for AtomDsl {
    type Error = ParseError;

    fn to_ast(&self, cfg: &AtomAstConfig) -> Result<AtomAst, ParseError> {
        let mut out = self.0.clone();
        raise_atom(&mut out, cfg);
        Ok(out)
    }
}

// -- Parse -------------------------------------------------------

/// Parse a complete atom-string into an `AtomDsl`.
pub fn parse_atom(input: &str) -> Result<AtomDsl, ParseError> {
    atom.parse(input).map_err(|e| e.into_inner())
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
        AtomConstraintKind::DonatedPairs => "#d",
        AtomConstraintKind::AcceptedPairs => "#t",
        AtomConstraintKind::AromaticValence => "#a",
        AtomConstraintKind::MulticenterValence => "#m",
        AtomConstraintKind::Degree => "#D",
        AtomConstraintKind::Connectivity => "#X",
        AtomConstraintKind::RingConnectivity => "#x",
        AtomConstraintKind::TotalHydrogens => "#H",
        AtomConstraintKind::RingCount => "#R",
        AtomConstraintKind::RingSize => "#r",
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AtomPredicate {
    IsotopeMass(IsotopeAst),
    Charge(ValueAst),
    ImplicitHydrogens(ImplicitHydrogensAst),
    LonePairs(ValueAst),
    Spin(SpinPredicate),
    Constraint(AtomConstraint),
}

fn atom_predicate(i: &mut &str) -> PResult<AtomPredicate> {
    let start = *i;
    let prefix: &str = take(2usize).parse_next(i)?;
    match prefix {
        "#i" => isotope.map(AtomPredicate::IsotopeMass).parse_next(i),
        "#c" => charge.map(AtomPredicate::Charge).parse_next(i),
        "#h" => implicit_hydrogens
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
            .map(|v| AtomPredicate::Constraint(AtomConstraint::Connectivity(v)))
            .parse_next(i),
        "#x" => optional_value
            .map(|v| AtomPredicate::Constraint(AtomConstraint::RingConnectivity(v)))
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
        element_set.map(ElementAst::Set),
        element_bind.map(|(id, set)| ElementAst::Bind { id, set }),
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

fn element_bind(i: &mut &str) -> PResult<(String, Vec<Element>)> {
    delimited(
        '(',
        (
            delimited(multispace0, preceded('?', id), multispace0),
            preceded(("::", multispace0), terminated(element_set, multispace0)),
        ),
        ')',
    )
    .parse_next(i)
}

fn element_ref(i: &mut &str) -> PResult<String> {
    delimited(
        '(',
        delimited(multispace0, preceded('?', id), multispace0),
        ')',
    )
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

fn implicit_hydrogens(i: &mut &str) -> PResult<ImplicitHydrogensAst> {
    preceded(
        multispace0,
        alt((
            "=".value(ImplicitHydrogensAst::Normal),
            value.map(ImplicitHydrogensAst::from),
            empty.value(ImplicitHydrogensAst::Lit(1)),
        )),
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
            "+".value(AromaticValenceAst::Aromatic(ValueAst::Expr(Expr::Rel(
                Box::new(Expr::Var("a".to_string())),
                RelOp::Ge,
                Box::new(Expr::Lit(0)),
            )))),
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
            "+".value(MulticenterValenceAst::Multicenter(ValueAst::Expr(
                Expr::Rel(
                    Box::new(Expr::Var("m".to_string())),
                    RelOp::Ge,
                    Box::new(Expr::Lit(0)),
                ),
            ))),
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
                if !matches!(ast.implicit_hydrogens, ImplicitHydrogensAst::Undetermined) {
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
                let kind = c.kind();
                if ast.constraints.contains(kind) {
                    return Err(ParseError::DuplicateAtomPredicate(
                        constraint_tag(kind).to_string(),
                    ));
                }
                ast.constraints.add(c);
            }
        }
    }
    Ok(())
}

// -- Format -------------------------------------------------------

fn fmt_atom_ast(f: &mut fmt::Formatter<'_>, ast: &AtomAst) -> fmt::Result {
    fmt_element(f, &ast.element)?;

    match &ast.isotope_mass {
        IsotopeAst::Undetermined => {}
        IsotopeAst::Natural => write!(f, "#i=")?,
        IsotopeAst::Lit(n) => write!(f, "#i{}", n)?,
        IsotopeAst::LitSet(s) => {
            write!(f, "#i")?;
            fmt_value(f, &ValueAst::LitSet(s.clone()))?;
        }
        IsotopeAst::Expr(e) => {
            write!(f, "#i")?;
            fmt_value(f, &ValueAst::Expr(e.clone()))?;
        }
    }

    fmt_charge(f, &ast.charge)?;

    match &ast.implicit_hydrogens {
        ImplicitHydrogensAst::Undetermined => {}
        ImplicitHydrogensAst::Normal => write!(f, "#h=")?,
        ImplicitHydrogensAst::Lit(1) => write!(f, "#h")?,
        ImplicitHydrogensAst::Lit(n) => write!(f, "#h{}", n)?,
        ImplicitHydrogensAst::LitSet(s) => {
            write!(f, "#h")?;
            fmt_value(f, &ValueAst::LitSet(s.clone()))?;
        }
        ImplicitHydrogensAst::Expr(e) => {
            write!(f, "#h")?;
            fmt_value(f, &ValueAst::Expr(e.clone()))?;
        }
    }

    fmt_value_field(f, "#n", &ast.lone_pairs)?;
    fmt_spin_pair(f, &ast.spin)
}

fn fmt_element(f: &mut fmt::Formatter<'_>, expr: &ElementAst) -> fmt::Result {
    match expr {
        ElementAst::Lit(e) => write!(f, "{}", e),
        ElementAst::Undetermined => write!(f, "*"),
        ElementAst::Set(es) => {
            write!(f, "{{")?;
            for (i, e) in es.iter().enumerate() {
                if i > 0 {
                    write!(f, ",")?;
                }
                write!(f, "{}", e)?;
            }
            write!(f, "}}")
        }
        ElementAst::Bind { id, set } => {
            write!(f, "(?{} :: {{", id)?;
            for (i, e) in set.iter().enumerate() {
                if i > 0 {
                    write!(f, ",")?;
                }
                write!(f, "{}", e)?;
            }
            write!(f, "}})")
        }
        ElementAst::Ref(id) => write!(f, "(?{})", id),
    }
}

/// Format a value field that suppresses zero (DSL convention for AST fields with
/// implicit-zero defaults like lone_pairs).
fn fmt_value_field(f: &mut fmt::Formatter<'_>, prefix: &str, v: &ValueAst) -> fmt::Result {
    match v {
        ValueAst::Undetermined | ValueAst::Lit(0) => Ok(()),
        ValueAst::Lit(1) => write!(f, "{}", prefix),
        ValueAst::Lit(n) => write!(f, "{}{}", prefix, n),
        v => {
            write!(f, "{}", prefix)?;
            fmt_value(f, v)
        }
    }
}

/// Format a value field that always emits (constraint sugar — zero is meaningful).
fn fmt_value_field_required(f: &mut fmt::Formatter<'_>, prefix: &str, v: &ValueAst) -> fmt::Result {
    match v {
        ValueAst::Undetermined => write!(f, "{}*", prefix),
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
            MulticenterValenceAst::Undetermined => write!(f, "#m*"),
            MulticenterValenceAst::NotMulticenter => write!(f, "#m!"),
            MulticenterValenceAst::Multicenter(v) if is_plus_sugar(v, "m", 0) => {
                write!(f, "#m+")
            }
            MulticenterValenceAst::Multicenter(ValueAst::Lit(1)) => write!(f, "#m"),
            MulticenterValenceAst::Multicenter(ValueAst::Lit(n)) => write!(f, "#m{}", n),
            MulticenterValenceAst::Multicenter(v) => {
                write!(f, "#m")?;
                fmt_value(f, v)
            }
        },
        AtomConstraint::AromaticValence(c) => match c {
            AromaticValenceAst::Undetermined => write!(f, "#a*"),
            AromaticValenceAst::NotAromatic => write!(f, "#a!"),
            AromaticValenceAst::Aromatic(v) if is_plus_sugar(v, "a", 0) => {
                write!(f, "#a+")
            }
            AromaticValenceAst::Aromatic(ValueAst::Lit(1)) => write!(f, "#a"),
            AromaticValenceAst::Aromatic(ValueAst::Lit(n)) => write!(f, "#a{}", n),
            AromaticValenceAst::Aromatic(v) => {
                write!(f, "#a")?;
                fmt_value(f, v)
            }
        },
        AtomConstraint::Degree(v) => fmt_value_field_required(f, "#D", v),
        AtomConstraint::Connectivity(v) => fmt_value_field_required(f, "#X", v),
        AtomConstraint::RingConnectivity(v) => fmt_value_field_required(f, "#x", v),
        AtomConstraint::TotalHydrogens(v) => fmt_value_field_required(f, "#H", v),
        AtomConstraint::RingCount(v) => fmt_ring_count(f, v),
        AtomConstraint::RingSize(v) => fmt_value_field_required(f, "#r", v),
    }
}

// -- Raise -------------------------------------------------------

fn raise_atom(ast: &mut AtomAst, cfg: &AtomAstConfig) {
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
        *isotope_mass = match cfg.isotope_mode {
            IsotopeMode::Natural => IsotopeAst::Natural,
            IsotopeMode::Required => IsotopeAst::Undetermined,
        };
    }
    if matches!(*charge, ValueAst::Undetermined) {
        *charge = match cfg.charge_mode {
            NumericMode::Zero => ValueAst::Lit(0),
            NumericMode::Required => ValueAst::Undetermined,
        };
    }
    if matches!(*implicit_hydrogens, ImplicitHydrogensAst::Undetermined) {
        *implicit_hydrogens = match cfg.implicit_h_mode {
            ImplicitHydrogenMode::Normal => ImplicitHydrogensAst::Normal,
            ImplicitHydrogenMode::Zero => ImplicitHydrogensAst::Lit(0),
            ImplicitHydrogenMode::Required => ImplicitHydrogensAst::Undetermined,
        };
    }
    if matches!(*lone_pairs, ValueAst::Undetermined) {
        *lone_pairs = match cfg.lone_pairs_mode {
            NumericMode::Zero => ValueAst::Lit(0),
            NumericMode::Required => ValueAst::Undetermined,
        };
    }
    raise_spin(spin, cfg.unpaired_electrons_mode, cfg.multiplicity_mode);
    raise_atom_constraints(constraints, cfg);
}

fn raise_atom_constraints(constraints: &mut AtomConstraints, cfg: &AtomAstConfig) {
    constraints.retain(|c| !c.is_undetermined());

    // Exhaustive dispatch over every kind: a new AtomConstraintKind variant
    // fails to build here until it has an explicit branch.
    for kind in AtomConstraintKind::iter() {
        match kind {
            AtomConstraintKind::Valence => {
                if matches!(cfg.valence_mode, NumericMode::Zero) && !constraints.contains(kind) {
                    constraints.add(AtomConstraint::Valence(ValueAst::Lit(0)));
                }
            }
            AtomConstraintKind::DonatedPairs => {
                if matches!(cfg.donated_pairs_mode, NumericMode::Zero)
                    && !constraints.contains(kind)
                {
                    constraints.add(AtomConstraint::DonatedPairs(ValueAst::Lit(0)));
                }
            }
            AtomConstraintKind::AcceptedPairs => {
                if matches!(cfg.accepted_pairs_mode, NumericMode::Zero)
                    && !constraints.contains(kind)
                {
                    constraints.add(AtomConstraint::AcceptedPairs(ValueAst::Lit(0)));
                }
            }
            AtomConstraintKind::AromaticValence => {
                if !constraints.contains(kind) {
                    match cfg.aromatic_valence_mode {
                        AromaticValenceMode::NotAromatic => {
                            constraints.add(AtomConstraint::AromaticValence(
                                AromaticValenceAst::NotAromatic,
                            ));
                        }
                        AromaticValenceMode::Aromatic => {
                            constraints.add(AtomConstraint::AromaticValence(
                                AromaticValenceAst::Aromatic(ValueAst::Undetermined),
                            ));
                        }
                        AromaticValenceMode::Required => {}
                    }
                }
            }
            AtomConstraintKind::MulticenterValence => {
                if !constraints.contains(kind) {
                    match cfg.multicenter_valence_mode {
                        MulticenterValenceMode::NotMulticenter => {
                            constraints.add(AtomConstraint::MulticenterValence(
                                MulticenterValenceAst::NotMulticenter,
                            ));
                        }
                        MulticenterValenceMode::Multicenter => {
                            constraints.add(AtomConstraint::MulticenterValence(
                                MulticenterValenceAst::Multicenter(ValueAst::Undetermined),
                            ));
                        }
                        MulticenterValenceMode::Required => {}
                    }
                }
            }
            AtomConstraintKind::Degree
            | AtomConstraintKind::Connectivity
            | AtomConstraintKind::RingConnectivity
            | AtomConstraintKind::TotalHydrogens
            | AtomConstraintKind::RingCount
            | AtomConstraintKind::RingSize => {
                // Pattern-only constraint: no defaulting mode in AtomAstConfig.
            }
        }
    }
}

// -- Lower -------------------------------------------------------

fn lower_atom(ast: &mut AtomAst, cfg: &AtomAstConfig) {
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
        (&cfg.isotope_mode, &*isotope_mass),
        (IsotopeMode::Natural, IsotopeAst::Natural)
    ) {
        *isotope_mass = IsotopeAst::Undetermined;
    }
    if matches!(
        (&cfg.charge_mode, &*charge),
        (NumericMode::Zero, ValueAst::Lit(0))
    ) {
        *charge = ValueAst::Undetermined;
    }
    match (&cfg.implicit_h_mode, &*implicit_hydrogens) {
        (ImplicitHydrogenMode::Normal, ImplicitHydrogensAst::Normal) => {
            *implicit_hydrogens = ImplicitHydrogensAst::Undetermined;
        }
        (ImplicitHydrogenMode::Zero, ImplicitHydrogensAst::Lit(0)) => {
            *implicit_hydrogens = ImplicitHydrogensAst::Undetermined;
        }
        _ => {}
    }
    if matches!(
        (&cfg.lone_pairs_mode, &*lone_pairs),
        (NumericMode::Zero, ValueAst::Lit(0))
    ) {
        *lone_pairs = ValueAst::Undetermined;
    }
    lower_spin(spin, cfg.unpaired_electrons_mode, cfg.multiplicity_mode);
    lower_atom_constraints(constraints, cfg);
}

fn lower_atom_constraints(constraints: &mut AtomConstraints, cfg: &AtomAstConfig) {
    // Exhaustive dispatch over every kind: a new AtomConstraintKind variant
    // fails to build here until it has an explicit branch.
    for kind in AtomConstraintKind::iter() {
        match kind {
            AtomConstraintKind::Valence => {
                if matches!(cfg.valence_mode, NumericMode::Zero)
                    && matches!(
                        constraints.get(kind),
                        Some(AtomConstraint::Valence(ValueAst::Lit(0)))
                    )
                {
                    constraints.remove(kind);
                }
            }
            AtomConstraintKind::DonatedPairs => {
                if matches!(cfg.donated_pairs_mode, NumericMode::Zero)
                    && matches!(
                        constraints.get(kind),
                        Some(AtomConstraint::DonatedPairs(ValueAst::Lit(0)))
                    )
                {
                    constraints.remove(kind);
                }
            }
            AtomConstraintKind::AcceptedPairs => {
                if matches!(cfg.accepted_pairs_mode, NumericMode::Zero)
                    && matches!(
                        constraints.get(kind),
                        Some(AtomConstraint::AcceptedPairs(ValueAst::Lit(0)))
                    )
                {
                    constraints.remove(kind);
                }
            }
            AtomConstraintKind::MulticenterValence => match cfg.multicenter_valence_mode {
                MulticenterValenceMode::NotMulticenter => {
                    if matches!(
                        constraints.get(kind),
                        Some(AtomConstraint::MulticenterValence(
                            MulticenterValenceAst::NotMulticenter
                        ))
                    ) {
                        constraints.remove(kind);
                    }
                }
                MulticenterValenceMode::Multicenter => {
                    if let Some(AtomConstraint::MulticenterValence(
                        MulticenterValenceAst::Multicenter(v),
                    )) = constraints.get(kind)
                    {
                        if is_plus_sugar(v, "m", 0) {
                            constraints.remove(kind);
                        }
                    }
                }
                MulticenterValenceMode::Required => {}
            },
            AtomConstraintKind::AromaticValence => match cfg.aromatic_valence_mode {
                AromaticValenceMode::NotAromatic => {
                    if matches!(
                        constraints.get(kind),
                        Some(AtomConstraint::AromaticValence(
                            AromaticValenceAst::NotAromatic
                        ))
                    ) {
                        constraints.remove(kind);
                    }
                }
                AromaticValenceMode::Aromatic => {
                    if let Some(AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(v))) =
                        constraints.get(kind)
                    {
                        if is_plus_sugar(v, "a", 0) {
                            constraints.remove(kind);
                        }
                    }
                }
                AromaticValenceMode::Required => {}
            },
            AtomConstraintKind::Degree
            | AtomConstraintKind::Connectivity
            | AtomConstraintKind::RingConnectivity
            | AtomConstraintKind::TotalHydrogens
            | AtomConstraintKind::RingCount
            | AtomConstraintKind::RingSize => {
                // Pattern-only constraint: no defaulting mode in AtomAstConfig.
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::element::Element;

    use super::*;
    use crate::ast::spin::SpinStateAst;
    use crate::ast::value::{Expr, RelOp};

    #[rustfmt::skip]
    #[rstest]
    #[case::carbon("C", AtomDsl(AtomAst::new(ElementAst::Lit(Element::C))))]
    #[case::iron("Fe", AtomDsl(AtomAst::new(ElementAst::Lit(Element::Fe))))]
    #[case::chlorine("Cl", AtomDsl(AtomAst::new(ElementAst::Lit(Element::Cl))))]
    #[case::whitespace("  C  ", AtomDsl(AtomAst::new(ElementAst::Lit(Element::C))))]
    #[case::undetermined("*", AtomDsl(AtomAst::new(ElementAst::Undetermined)))]
    #[case::element_set("{C,N,O}", AtomDsl(AtomAst::new(ElementAst::Set(vec![Element::C, Element::N, Element::O]))))]
    #[case::element_bind("(?e :: {C,N})", AtomDsl(AtomAst::new(ElementAst::Bind { id: "e".to_string(), set: vec![Element::C, Element::N] })))]
    #[case::element_ref("(?e)", AtomDsl(AtomAst::new(ElementAst::Ref("e".to_string()))))]
    #[case::isotope("C#i12", AtomDsl(AtomAst { isotope_mass: IsotopeAst::Lit(12), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::isotope_natural("C#i=", AtomDsl(AtomAst { isotope_mass: IsotopeAst::Natural, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_pos("C#c+2", AtomDsl(AtomAst { charge: ValueAst::Lit(2), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_neg("C#c-2", AtomDsl(AtomAst { charge: ValueAst::Lit(-2), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_plus("C#c+", AtomDsl(AtomAst { charge: ValueAst::Lit(1), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_minus("C#c-", AtomDsl(AtomAst { charge: ValueAst::Lit(-1), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_zero("C#c0", AtomDsl(AtomAst { charge: ValueAst::Lit(0), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_count("C#h3", AtomDsl(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Lit(3), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_normal("C#h=", AtomDsl(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Normal, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_undetermined("C#h*", AtomDsl(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Undetermined, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_bind("C#h(?h)", AtomDsl(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Expr(Expr::Var("h".to_string())), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_set("N#h?h :: {2,3}", AtomDsl(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Expr(Expr::Mem(Box::new(Expr::Var("h".to_string())), vec![2, 3])), ..AtomAst::new(ElementAst::Lit(Element::N)) }))]
    #[case::h_expr("C#h?h >= 1", AtomDsl(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Expr(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_omit("C#h", AtomDsl(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Lit(1), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
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
    #[case::arom_plus("C#a+", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("a".to_string())), RelOp::Ge, Box::new(Expr::Lit(0))))))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::arom_zero("C#a0", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(0)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::arom_one("C#a1", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(1)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::arom_omit("C#a", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(1)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multicenter_not("C#m!", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::MulticenterValence(MulticenterValenceAst::NotMulticenter)]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multicenter_undetermined("C#m*", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::MulticenterValence(MulticenterValenceAst::Undetermined)]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multicenter_plus("C#m+", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("m".to_string())), RelOp::Ge, Box::new(Expr::Lit(0))))))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multicenter_zero("C#m0", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(0)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multicenter_one("C#m", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(1)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multicenter("C#m2", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(2)))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::degree("C#D2", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::Degree(ValueAst::Lit(2))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::connectivity("C#X3", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::Connectivity(ValueAst::Lit(3))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_connectivity("C#x2", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::RingConnectivity(ValueAst::Lit(2))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::total_hydrogens("C#H1", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::TotalHydrogens(ValueAst::Lit(1))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_bare("C#R", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::RingCount(ValueAst::Lit(1))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_undetermined("C#R*", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::RingCount(ValueAst::Undetermined)]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_plus("C#R+", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::RingCount(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("r".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::ring_count("C#R2", AtomDsl(AtomAst { constraints: AtomConstraints::from_iter([AtomConstraint::RingCount(ValueAst::Lit(2))]), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    fn test_parse_atom(#[case] input: &str, #[case] expected: AtomDsl) {
        let result = atom.parse(input);
        assert!(result.is_ok(), "{:?} should succeed, got {:?}", input, result.clone().unwrap_err());
        let form = result.unwrap();
        assert_eq!(form, expected);
    }

    #[rstest]
    #[case::empty("", ParseError::ExpectedElement)]
    #[case::no_element("#h3", ParseError::ExpectedElement)]
    #[case::unknown_pred("C#y", ParseError::UnknownAtomPredicate("#y".to_string()))]
    #[case::dup_h("C#h3#h2", ParseError::DuplicateAtomPredicate("#h".to_string()))]
    #[case::dup_charge("C#c+#c-", ParseError::DuplicateAtomPredicate("#c".to_string()))]
    #[case::dup_valence("C#v3#v4", ParseError::DuplicateAtomPredicate("#v".to_string()))]
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
    #[case::arom_undetermined("C#a*")]
    #[case::multicenter_not("C#m!")]
    #[case::multicenter_plus("C#m+")]
    #[case::multicenter_zero("C#m0")]
    #[case::multicenter_omit("C#m")]
    #[case::multicenter_undetermined("C#m*")]
    #[case::ring_bare("C#R")]
    #[case::ring_plus("C#R+")]
    #[case::ring_count("C#R2")]
    #[case::ring_undetermined("C#R*")]
    fn test_atom_display_roundtrip(#[case] input: &str) {
        let parsed = atom.parse(input).unwrap();
        assert_eq!(parsed.to_string(), input);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::carbon("C", ElementAst::Lit(Element::C))]
    #[case::iron("Fe", ElementAst::Lit(Element::Fe))]
    #[case::chlorine("Cl", ElementAst::Lit(Element::Cl))]
    #[case::undetermined("*", ElementAst::Undetermined)]
    #[case::set("{C,N,O}", ElementAst::Set(vec![Element::C, Element::N, Element::O]))]
    #[case::set_spaced("{ C, N}", ElementAst::Set(vec![Element::C, Element::N]))]
    #[case::bind("(?e :: {C,N})", ElementAst::Bind { id: "e".to_string(), set: vec![Element::C, Element::N] })]
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
    #[case::set("{12,13,14}", IsotopeAst::LitSet(vec![12, 13, 14]))]
    #[case::bind("(?m :: {12,13})", IsotopeAst::Expr(Expr::Mem(Box::new(Expr::Var("m".to_string())), vec![12, 13])))]
    #[case::ref_("(?m)", IsotopeAst::Expr(Expr::Var("m".to_string())))]
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
    #[case::h_count("#h3", AtomPredicate::ImplicitHydrogens(ImplicitHydrogensAst::Lit(3)))]
    #[case::h_normal("#h=", AtomPredicate::ImplicitHydrogens(ImplicitHydrogensAst::Normal))]
    #[case::h_undetermined("#h*", AtomPredicate::ImplicitHydrogens(ImplicitHydrogensAst::Undetermined))]
    #[case::h_omit("#h", AtomPredicate::ImplicitHydrogens(ImplicitHydrogensAst::Lit(1)))]
    #[case::lone_pairs("#n2", AtomPredicate::LonePairs(ValueAst::Lit(2)))]
    #[case::lone_pairs_omit("#n", AtomPredicate::LonePairs(ValueAst::Lit(1)))]
    #[case::unpaired("#u2", AtomPredicate::Spin(SpinPredicate::Unpaired(ValueAst::Lit(2))))]
    #[case::unpaired_omit("#u", AtomPredicate::Spin(SpinPredicate::Unpaired(ValueAst::Lit(1))))]
    #[case::multiplicity("#s3", AtomPredicate::Spin(SpinPredicate::Multiplicity(ValueAst::Lit(3))))]
    #[case::multiplicity_omit("#s", AtomPredicate::Spin(SpinPredicate::Multiplicity(ValueAst::Lit(1))))]
    #[case::valence("#v4", AtomPredicate::Constraint(AtomConstraint::Valence(ValueAst::Lit(4))))]
    #[case::donated_pairs("#d1", AtomPredicate::Constraint(AtomConstraint::DonatedPairs(ValueAst::Lit(1))))]
    #[case::accepted_pairs("#t1", AtomPredicate::Constraint(AtomConstraint::AcceptedPairs(ValueAst::Lit(1))))]
    #[case::ring_size("#r6", AtomPredicate::Constraint(AtomConstraint::RingSize(ValueAst::Lit(6))))]
    #[case::arom_not_aromatic("#a!", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceAst::NotAromatic)))]
    #[case::arom_undetermined("#a*", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceAst::Undetermined)))]
    #[case::arom_plus("#a+", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("a".to_string())), RelOp::Ge, Box::new(Expr::Lit(0))))))))]
    #[case::arom_lit("#a2", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(2)))))]
    #[case::arom_omit("#a", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceAst::Aromatic(ValueAst::Lit(1)))))]
    #[case::multicenter_not("#m!", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(MulticenterValenceAst::NotMulticenter)))]
    #[case::multicenter_undetermined("#m*", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(MulticenterValenceAst::Undetermined)))]
    #[case::multicenter_plus("#m+", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("m".to_string())), RelOp::Ge, Box::new(Expr::Lit(0))))))))]
    #[case::multicenter_zero("#m0", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(0)))))]
    #[case::multicenter_omit("#m", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(1)))))]
    #[case::multicenter("#m2", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(MulticenterValenceAst::Multicenter(ValueAst::Lit(2)))))]
    #[case::degree("#D2", AtomPredicate::Constraint(AtomConstraint::Degree(ValueAst::Lit(2))))]
    #[case::degree_omit("#D", AtomPredicate::Constraint(AtomConstraint::Degree(ValueAst::Lit(1))))]
    #[case::connectivity("#X3", AtomPredicate::Constraint(AtomConstraint::Connectivity(ValueAst::Lit(3))))]
    #[case::ring_connectivity("#x2", AtomPredicate::Constraint(AtomConstraint::RingConnectivity(ValueAst::Lit(2))))]
    #[case::ring_connectivity_omit("#x", AtomPredicate::Constraint(AtomConstraint::RingConnectivity(ValueAst::Lit(1))))]
    #[case::total_hydrogens("#H1", AtomPredicate::Constraint(AtomConstraint::TotalHydrogens(ValueAst::Lit(1))))]
    #[case::ring_bare("#R", AtomPredicate::Constraint(AtomConstraint::RingCount(ValueAst::Lit(1))))]
    #[case::ring_undetermined("#R*", AtomPredicate::Constraint(AtomConstraint::RingCount(ValueAst::Undetermined)))]
    #[case::ring_plus("#R+", AtomPredicate::Constraint(AtomConstraint::RingCount(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("r".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))))))]
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
        let cfg = AtomAstConfig::zeroed();
        let ast = dsl.to_ast(&cfg).unwrap();
        assert_eq!(ast.charge, ValueAst::Lit(0));
        assert_eq!(ast.lone_pairs, ValueAst::Lit(0));
        assert_eq!(ast.implicit_hydrogens, ImplicitHydrogensAst::Lit(0));
        assert_eq!(ast.isotope_mass, IsotopeAst::Natural);
        assert_eq!(ast.spin, SpinStateAst::new(0, 1));
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
        ast.implicit_hydrogens = ImplicitHydrogensAst::Lit(0);
        ast.isotope_mass = IsotopeAst::Natural;
        ast.spin = SpinStateAst::new(0, 1);
        ast.constraints
            .add(AtomConstraint::Valence(ValueAst::Lit(0)));
        ast.constraints.add(AtomConstraint::AromaticValence(
            AromaticValenceAst::NotAromatic,
        ));
        let cfg = AtomAstConfig::zeroed();
        let dsl = AtomDsl::from_ast(&ast, &cfg).unwrap();
        assert_eq!(dsl.0.charge, ValueAst::Undetermined);
        assert_eq!(dsl.0.lone_pairs, ValueAst::Undetermined);
        assert_eq!(dsl.0.implicit_hydrogens, ImplicitHydrogensAst::Undetermined);
        assert_eq!(dsl.0.isotope_mass, IsotopeAst::Undetermined);
        assert_eq!(dsl.0.spin, SpinStateAst::default());
        assert!(dsl.0.constraints.is_empty());
    }

    #[rstest]
    fn test_atom_dsl_roundtrip_zeroed() {
        let input = AtomDsl(AtomAst::new(ElementAst::Lit(Element::C)));
        let cfg = AtomAstConfig::zeroed();
        let raised = input.to_ast(&cfg).unwrap();
        let lowered = AtomDsl::from_ast(&raised, &cfg).unwrap();
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
}
