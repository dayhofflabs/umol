//! Atom-string DSL.

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use umol_edn::{DeError, Edn, FromEdn, ToEdn};
use umol_shared::element::Element;
use winnow::ascii::multispace0;
use winnow::combinator::{alt, delimited, empty, preceded, repeat, separated, terminated};
use winnow::error::{ErrMode, ParserError};
use winnow::token::{one_of, take};
use winnow::Parser;

use crate::ast::atom::{AtomAst, ElementAst, ImplicitHydrogensAst, IsotopeAst};
use crate::ast::constraint::{AromaticValenceConstraint, AtomConstraint, MulticenterValenceConstraint};
use crate::ast::value::{Expr, RelOp, ValueAst};
use super::error::{PResult, ParseError};
use super::predicates::{
    apply_spin_pair, charge, fmt_charge, fmt_ring_count, fmt_spin_pair, fmt_value, is_plus_sugar,
    optional_value, ring_count, SpinPredicate,
};
use super::value::{id, value};

/// `AtomAst` combined with the atom-level constraints.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AtomTypeAst {
    pub ast: AtomAst,
    pub constraints: Vec<AtomConstraint>,
}

impl AtomTypeAst {
    pub fn new(ast: AtomAst) -> Self {
        Self {
            ast,
            constraints: Vec::new(),
        }
    }

    pub fn with_constraints(ast: AtomAst, constraints: Vec<AtomConstraint>) -> Self {
        Self { ast, constraints }
    }
}

impl FromStr for AtomTypeAst {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_atom(s)
    }
}

impl Display for AtomTypeAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_atom_ast(f, &self.ast)?;
        for c in &self.constraints {
            fmt_constraint(f, c)?;
        }
        Ok(())
    }
}

impl<'de> FromEdn<'de> for AtomTypeAst {
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

impl ToEdn for AtomTypeAst {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

/// Parse a complete atom-string into an `AtomTypeAst` (base AST + lifted constraints).
pub fn parse_atom(input: &str) -> Result<AtomTypeAst, ParseError> {
    atom.parse(input).map_err(|e| e.into_inner())
}

/// Atom-string parser (does not require consuming all input).
pub(crate) fn atom(i: &mut &str) -> PResult<AtomTypeAst> {
    let el = delimited(multispace0, element, multispace0).parse_next(i)?;
    let preds: Vec<AtomPredicate> =
        repeat(0.., terminated(atom_predicate, multispace0)).parse_next(i)?;
    let mut form = AtomTypeAst::new(AtomAst::new(el));
    apply_predicates(&mut form, preds).map_err(ErrMode::Cut)?;
    Ok(form)
}

fn is_set(v: &ValueAst) -> bool {
    !matches!(v, ValueAst::Undetermined)
}

fn apply_predicates(form: &mut AtomTypeAst, preds: Vec<AtomPredicate>) -> Result<(), ParseError> {
    let ast = &mut form.ast;
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
                let tag = constraint_tag(&c);
                if form
                    .constraints
                    .iter()
                    .any(|existing| constraint_tag(existing) == tag)
                {
                    return Err(ParseError::DuplicateAtomPredicate(tag.to_string()));
                }
                form.constraints.push(c);
            }
        }
    }
    Ok(())
}

fn constraint_tag(c: &AtomConstraint) -> &'static str {
    match c {
        AtomConstraint::Valence(_) => "#v",
        AtomConstraint::DonatedPairs(_) => "#d",
        AtomConstraint::AcceptedPairs(_) => "#t",
        AtomConstraint::AromaticValence(_) => "#a",
        AtomConstraint::MulticenterValence(_) => "#m",
        AtomConstraint::Degree(_) => "#D",
        AtomConstraint::Connectivity(_) => "#X",
        AtomConstraint::RingConnectivity(_) => "#x",
        AtomConstraint::TotalHydrogens(_) => "#H",
        AtomConstraint::RingCount(_) => "#R",
        AtomConstraint::RingSize(_) => "#r",
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
        alt(('='.value(IsotopeAst::Natural), value.map(IsotopeAst::Value))),
    )
    .parse_next(i)
    .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedPredicateBody))
}

fn implicit_hydrogens(i: &mut &str) -> PResult<ImplicitHydrogensAst> {
    preceded(
        multispace0,
        alt((
            "=".value(ImplicitHydrogensAst::Normal),
            value.map(ImplicitHydrogensAst::Value),
            empty.value(ImplicitHydrogensAst::Value(ValueAst::Lit(1))),
        )),
    )
    .parse_next(i)
    .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedPredicateBody))
}

fn aromatic_valence(i: &mut &str) -> PResult<AromaticValenceConstraint> {
    preceded(
        multispace0,
        alt((
            "!".value(AromaticValenceConstraint::NotAromatic),
            value.map(AromaticValenceConstraint::Value),
            "+".value(AromaticValenceConstraint::Value(ValueAst::Expr(Expr::Rel(
                Box::new(Expr::Var("a".to_string())),
                RelOp::Ge,
                Box::new(Expr::Lit(0)),
            )))),
            empty.value(AromaticValenceConstraint::Value(ValueAst::Lit(1))),
        )),
    )
    .parse_next(i)
    .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedPredicateBody))
}

fn multicenter_valence(i: &mut &str) -> PResult<MulticenterValenceConstraint> {
    preceded(
        multispace0,
        alt((
            "!".value(MulticenterValenceConstraint::NotMulticenter),
            value.map(MulticenterValenceConstraint::Value),
            "+".value(MulticenterValenceConstraint::Value(ValueAst::Expr(Expr::Rel(
                Box::new(Expr::Var("m".to_string())),
                RelOp::Ge,
                Box::new(Expr::Lit(0)),
            )))),
            empty.value(MulticenterValenceConstraint::Value(ValueAst::Lit(1))),
        )),
    )
    .parse_next(i)
    .map_err(|_: ErrMode<ParseError>| ErrMode::Backtrack(ParseError::ExpectedPredicateBody))
}

fn fmt_atom_ast(f: &mut fmt::Formatter<'_>, ast: &AtomAst) -> fmt::Result {
    fmt_element(f, &ast.element)?;

    match &ast.isotope_mass {
        IsotopeAst::Undetermined => {}
        IsotopeAst::Natural => write!(f, "#i=")?,
        IsotopeAst::Value(v) => {
            write!(f, "#i")?;
            fmt_value(f, v)?;
        }
    }

    fmt_charge(f, &ast.charge)?;

    match &ast.implicit_hydrogens {
        ImplicitHydrogensAst::Undetermined | ImplicitHydrogensAst::Value(ValueAst::Lit(0)) => {}
        ImplicitHydrogensAst::Normal => write!(f, "#h=")?,
        ImplicitHydrogensAst::Value(ValueAst::Lit(1)) => write!(f, "#h")?,
        ImplicitHydrogensAst::Value(ValueAst::Lit(n)) => write!(f, "#h{}", n)?,
        ImplicitHydrogensAst::Value(ValueAst::Undetermined) => write!(f, "#h*")?,
        ImplicitHydrogensAst::Value(v) => {
            write!(f, "#h")?;
            fmt_value(f, v)?;
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
            MulticenterValenceConstraint::NotMulticenter => write!(f, "#m!"),
            MulticenterValenceConstraint::Value(v) if is_plus_sugar(v, "m", 0) => write!(f, "#m+"),
            MulticenterValenceConstraint::Value(ValueAst::Lit(1)) => write!(f, "#m"),
            MulticenterValenceConstraint::Value(ValueAst::Lit(n)) => write!(f, "#m{}", n),
            MulticenterValenceConstraint::Value(v) => {
                write!(f, "#m")?;
                fmt_value(f, v)
            }
        },
        AtomConstraint::AromaticValence(c) => match c {
            AromaticValenceConstraint::NotAromatic => write!(f, "#a!"),
            AromaticValenceConstraint::Value(v) if is_plus_sugar(v, "a", 0) => write!(f, "#a+"),
            AromaticValenceConstraint::Value(ValueAst::Lit(1)) => write!(f, "#a"),
            AromaticValenceConstraint::Value(ValueAst::Lit(n)) => write!(f, "#a{}", n),
            AromaticValenceConstraint::Value(v) => {
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
    #[case::carbon("C", AtomTypeAst::new(AtomAst::new(ElementAst::Lit(Element::C))))]
    #[case::iron("Fe", AtomTypeAst::new(AtomAst::new(ElementAst::Lit(Element::Fe))))]
    #[case::chlorine("Cl", AtomTypeAst::new(AtomAst::new(ElementAst::Lit(Element::Cl))))]
    #[case::whitespace("  C  ", AtomTypeAst::new(AtomAst::new(ElementAst::Lit(Element::C))))]
    #[case::undetermined("*", AtomTypeAst::new(AtomAst::new(ElementAst::Undetermined)))]
    #[case::element_set("{C,N,O}", AtomTypeAst::new(AtomAst::new(ElementAst::Set(vec![Element::C, Element::N, Element::O]))))]
    #[case::element_bind("(?e :: {C,N})", AtomTypeAst::new(AtomAst::new(ElementAst::Bind { id: "e".to_string(), set: vec![Element::C, Element::N] })))]
    #[case::element_ref("(?e)", AtomTypeAst::new(AtomAst::new(ElementAst::Ref("e".to_string()))))]
    #[case::isotope("C#i12", AtomTypeAst::new(AtomAst { isotope_mass: IsotopeAst::Value(ValueAst::Lit(12)), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::isotope_natural("C#i=", AtomTypeAst::new(AtomAst { isotope_mass: IsotopeAst::Natural, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_pos("C#c+2", AtomTypeAst::new(AtomAst { charge: ValueAst::Lit(2), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_neg("C#c-2", AtomTypeAst::new(AtomAst { charge: ValueAst::Lit(-2), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_plus("C#c+", AtomTypeAst::new(AtomAst { charge: ValueAst::Lit(1), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_minus("C#c-", AtomTypeAst::new(AtomAst { charge: ValueAst::Lit(-1), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_zero("C#c0", AtomTypeAst::new(AtomAst { charge: ValueAst::Lit(0), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_count("C#h3", AtomTypeAst::new(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Value(ValueAst::Lit(3)), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_normal("C#h=", AtomTypeAst::new(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Normal, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_undetermined("C#h*", AtomTypeAst::new(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Value(ValueAst::Undetermined), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_bind("C#h(?h)", AtomTypeAst::new(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Value(ValueAst::Expr(Expr::Var("h".to_string()))), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_set("N#h?h :: {2,3}", AtomTypeAst::new(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Value(ValueAst::Expr(Expr::Mem(Box::new(Expr::Var("h".to_string())), vec![2, 3]))), ..AtomAst::new(ElementAst::Lit(Element::N)) }))]
    #[case::h_expr("C#h?h >= 1", AtomTypeAst::new(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Value(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1))))), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_omit("C#h", AtomTypeAst::new(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Value(ValueAst::Lit(1)), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::lone_pairs("O#n2", AtomTypeAst::new(AtomAst { lone_pairs: ValueAst::Lit(2), ..AtomAst::new(ElementAst::Lit(Element::O)) }))]
    #[case::lone_pairs_omit("O#n", AtomTypeAst::new(AtomAst { lone_pairs: ValueAst::Lit(1), ..AtomAst::new(ElementAst::Lit(Element::O)) }))]
    #[case::unpaired("C#u2", AtomTypeAst::new(AtomAst { spin: SpinStateAst { unpaired: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined }, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::unpaired_omit("C#u", AtomTypeAst::new(AtomAst { spin: SpinStateAst { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Undetermined }, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multiplicity("C#s3", AtomTypeAst::new(AtomAst { spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(3) }, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multiplicity_omit("C#s", AtomTypeAst::new(AtomAst { spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(1) }, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::valence("C#v4", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::Valence(ValueAst::Lit(4))]))]
    #[case::donated_pairs("N#d1", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::N)), vec![AtomConstraint::DonatedPairs(ValueAst::Lit(1))]))]
    #[case::accepted_pairs("B#t1", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::B)), vec![AtomConstraint::AcceptedPairs(ValueAst::Lit(1))]))]
    #[case::ring_size("C#r6", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::RingSize(ValueAst::Lit(6))]))]
    #[case::arom_not_aromatic("C#a!", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::AromaticValence(AromaticValenceConstraint::NotAromatic)]))]
    #[case::arom_undetermined("C#a*", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(ValueAst::Undetermined))]))]
    #[case::arom_plus("C#a+", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("a".to_string())), RelOp::Ge, Box::new(Expr::Lit(0))))))]))]
    #[case::arom_zero("C#a0", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(ValueAst::Lit(0)))]))]
    #[case::arom_one("C#a1", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(ValueAst::Lit(1)))]))]
    #[case::arom_omit("C#a", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(ValueAst::Lit(1)))]))]
    #[case::multicenter_not("C#m!", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::MulticenterValence(MulticenterValenceConstraint::NotMulticenter)]))]
    #[case::multicenter_undetermined("C#m*", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::MulticenterValence(MulticenterValenceConstraint::Value(ValueAst::Undetermined))]))]
    #[case::multicenter_plus("C#m+", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::MulticenterValence(MulticenterValenceConstraint::Value(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("m".to_string())), RelOp::Ge, Box::new(Expr::Lit(0))))))]))]
    #[case::multicenter_zero("C#m0", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::MulticenterValence(MulticenterValenceConstraint::Value(ValueAst::Lit(0)))]))]
    #[case::multicenter_one("C#m", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::MulticenterValence(MulticenterValenceConstraint::Value(ValueAst::Lit(1)))]))]
    #[case::multicenter("C#m2", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::MulticenterValence(MulticenterValenceConstraint::Value(ValueAst::Lit(2)))]))]
    #[case::degree("C#D2", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::Degree(ValueAst::Lit(2))]))]
    #[case::connectivity("C#X3", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::Connectivity(ValueAst::Lit(3))]))]
    #[case::ring_connectivity("C#x2", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::RingConnectivity(ValueAst::Lit(2))]))]
    #[case::total_hydrogens("C#H1", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::TotalHydrogens(ValueAst::Lit(1))]))]
    #[case::ring_bare("C#R", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::RingCount(ValueAst::Lit(1))]))]
    #[case::ring_undetermined("C#R*", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::RingCount(ValueAst::Undetermined)]))]
    #[case::ring_plus("C#R+", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::RingCount(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("r".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))))]))]
    #[case::ring_count("C#R2", AtomTypeAst::with_constraints(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::RingCount(ValueAst::Lit(2))]))]
    fn test_parse_atom(#[case] input: &str, #[case] expected: AtomTypeAst) {
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
    #[case::lit("12", IsotopeAst::Value(ValueAst::Lit(12)))]
    #[case::undetermined("*", IsotopeAst::Value(ValueAst::Undetermined))]
    #[case::set("{12,13,14}", IsotopeAst::Value(ValueAst::LitSet(vec![12, 13, 14])))]
    #[case::bind("(?m :: {12,13})", IsotopeAst::Value(ValueAst::Expr(Expr::Mem(Box::new(Expr::Var("m".to_string())), vec![12, 13]))))]
    #[case::ref_("(?m)", IsotopeAst::Value(ValueAst::Expr(Expr::Var("m".to_string()))))]
    fn test_isotope(#[case] input: &str, #[case] expected: IsotopeAst) {
        let result = isotope.parse(input);
        assert!(result.is_ok(), "{input:?} should succeed, got {:?}", result.unwrap_err());
        let expr = result.unwrap();
        assert_eq!(expr, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::isotope_lit("#i12", AtomPredicate::IsotopeMass(IsotopeAst::Value(ValueAst::Lit(12))))]
    #[case::isotope_natural("#i=", AtomPredicate::IsotopeMass(IsotopeAst::Natural))]
    #[case::isotope_undetermined("#i*", AtomPredicate::IsotopeMass(IsotopeAst::Value(ValueAst::Undetermined)))]
    #[case::charge_pos("#c+2", AtomPredicate::Charge(ValueAst::Lit(2)))]
    #[case::charge_neg("#c-2", AtomPredicate::Charge(ValueAst::Lit(-2)))]
    #[case::charge_plus("#c+", AtomPredicate::Charge(ValueAst::Lit(1)))]
    #[case::charge_minus("#c-", AtomPredicate::Charge(ValueAst::Lit(-1)))]
    #[case::charge_zero("#c0", AtomPredicate::Charge(ValueAst::Lit(0)))]
    #[case::charge_undetermined("#c*", AtomPredicate::Charge(ValueAst::Undetermined))]
    #[case::h_count("#h3", AtomPredicate::ImplicitHydrogens(ImplicitHydrogensAst::Value(ValueAst::Lit(3))))]
    #[case::h_normal("#h=", AtomPredicate::ImplicitHydrogens(ImplicitHydrogensAst::Normal))]
    #[case::h_undetermined("#h*", AtomPredicate::ImplicitHydrogens(ImplicitHydrogensAst::Value(ValueAst::Undetermined)))]
    #[case::h_omit("#h", AtomPredicate::ImplicitHydrogens(ImplicitHydrogensAst::Value(ValueAst::Lit(1))))]
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
    #[case::arom_not_aromatic("#a!", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceConstraint::NotAromatic)))]
    #[case::arom_undetermined("#a*", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(ValueAst::Undetermined))))]
    #[case::arom_plus("#a+", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("a".to_string())), RelOp::Ge, Box::new(Expr::Lit(0))))))))]
    #[case::arom_lit("#a2", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(ValueAst::Lit(2)))))]
    #[case::arom_omit("#a", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(ValueAst::Lit(1)))))]
    #[case::multicenter_not("#m!", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(MulticenterValenceConstraint::NotMulticenter)))]
    #[case::multicenter_undetermined("#m*", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(MulticenterValenceConstraint::Value(ValueAst::Undetermined))))]
    #[case::multicenter_plus("#m+", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(MulticenterValenceConstraint::Value(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("m".to_string())), RelOp::Ge, Box::new(Expr::Lit(0))))))))]
    #[case::multicenter_zero("#m0", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(MulticenterValenceConstraint::Value(ValueAst::Lit(0)))))]
    #[case::multicenter_omit("#m", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(MulticenterValenceConstraint::Value(ValueAst::Lit(1)))))]
    #[case::multicenter("#m2", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(MulticenterValenceConstraint::Value(ValueAst::Lit(2)))))]
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
    #[case::unknown("#y", ParseError::UnknownAtomPredicate("#y".to_string()))]
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
}
