//! Atom-string DSL: parser, AST, and display

use std::borrow::Cow;
use std::fmt::{self, Display};
use std::str::FromStr;

use nom::branch::alt;
use nom::bytes::complete::{tag, take};
use nom::character::complete::{char, multispace0, satisfy, u32 as nom_u32};
use nom::combinator::{all_consuming, map, recognize, success, value};
use nom::error::{Error as NomError, ErrorKind};
use nom::multi::{many0, separated_list1};
use nom::sequence::{delimited, pair, preceded, terminated};
use nom::{Err, IResult, Parser};
use umol_ast::ast::atom::{ElementAst, ImplicitHydrogensAst, IsotopeAst};
use umol_shared::element::Element;
use umol_ast::ast::spin::SpinStateAst;
use umol_ast::ast::value::ValueAst;
use umol_edn::{DeError, Edn, FromEdn, ToEdn};

use super::error::AtomDslError;
use super::value::{op_char, parse_id, value_dsl};
use crate::api::pattern::AtomPattern;
use crate::ast::atom::AtomAst;
use crate::ast::config::AtomAstConfig;
use crate::ast::constraint::{AromaticValenceConstraint, AtomConstraint};
use crate::ast::error::LoweringError;
use crate::ast::{FromAst, ToAst};

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct AtomDsl(String);

impl AtomDsl {
    pub fn from_pattern(pattern: AtomPattern) -> Self {
        Self(pattern.to_string())
    }

    pub fn from_parts(ast: AtomAst, constraints: Vec<AtomConstraint>) -> Self {
        Self(AtomPattern::with_constraints(ast, constraints).to_string())
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }

    pub fn into_pattern(self) -> Result<AtomPattern, AtomDslError> {
        parse_atom_dsl(&self.0)
    }

    pub fn lower_parts(self) -> Result<(AtomAst, Vec<AtomConstraint>), AtomDslError> {
        let pattern = self.into_pattern()?;
        Ok((pattern.ast, pattern.constraints))
    }

    pub fn into_ast(self) -> Result<AtomAst, AtomDslError> {
        let (ast, constraints) = self.lower_parts()?;
        if !constraints.is_empty() {
            return Err(AtomDslError::ConstraintsNotAllowed);
        }
        Ok(ast)
    }

    pub fn is_inline_constraint(constraint: &AtomConstraint) -> bool {
        matches!(
            constraint,
            AtomConstraint::Valence(_)
                | AtomConstraint::DonatedPairs(_)
                | AtomConstraint::AcceptedPairs(_)
                | AtomConstraint::MulticenterValence(_)
                | AtomConstraint::AromaticValence(_)
                | AtomConstraint::Degree(_)
                | AtomConstraint::Connectivity(_)
                | AtomConstraint::TotalHCount(_)
                | AtomConstraint::InRing
                | AtomConstraint::RingCount(_)
        )
    }
}

impl Display for AtomDsl {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl FromStr for AtomDsl {
    type Err = AtomDslError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_atom_dsl(s)?;
        Ok(Self(s.to_string()))
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
}

impl ToEdn for AtomDsl {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.0.clone()))
    }
}

impl FromAst<AtomAst> for AtomDsl {
    fn from_ast(ast: &AtomAst, cfg: &AtomAstConfig) -> Result<Self, LoweringError> {
        let mut ast = ast.clone();
        ast.release(cfg);
        Ok(Self(ast.to_string()))
    }
}

impl ToAst<AtomAst> for AtomDsl {
    fn to_ast(&self, _cfg: &AtomAstConfig) -> Result<AtomAst, LoweringError> {
        self.clone()
            .into_ast()
            .map_err(|e| LoweringError::Custom(e.to_string()))
    }
}

impl FromStr for AtomPattern {
    type Err = AtomDslError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_atom_dsl(s)
    }
}

impl FromStr for AtomAst {
    type Err = AtomDslError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let pat = parse_atom_dsl(s)?;
        if !pat.constraints.is_empty() {
            return Err(AtomDslError::ConstraintsNotAllowed);
        }
        Ok(pat.ast)
    }
}

impl Display for AtomAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_element(f, &self.element)?;

        match &self.isotope_mass {
            IsotopeAst::Undetermined => {}
            IsotopeAst::Natural => write!(f, "#i=")?,
            IsotopeAst::Lit(n) => write!(f, "#i{}", n)?,
            IsotopeAst::Set(ns) => {
                write!(f, "#i{{")?;
                for (i, n) in ns.iter().enumerate() {
                    if i > 0 {
                        write!(f, ",")?;
                    }
                    write!(f, "{}", n)?;
                }
                write!(f, "}}")?;
            }
            IsotopeAst::Bind { id, set } => {
                write!(f, "#i(?{} :: {{", id)?;
                for (i, n) in set.iter().enumerate() {
                    if i > 0 {
                        write!(f, ",")?;
                    }
                    write!(f, "{}", n)?;
                }
                write!(f, "}})")?;
            }
            IsotopeAst::Ref(id) => write!(f, "#i(?{})", id)?,
        }

        match &self.charge {
            ValueAst::Undetermined | ValueAst::Lit(0) => {}
            ValueAst::Lit(1) => write!(f, "#c+")?,
            ValueAst::Lit(-1) => write!(f, "#c-")?,
            ValueAst::Lit(n) if *n > 0 => write!(f, "#c+{}", n)?,
            ValueAst::Lit(n) => write!(f, "#c{}", n)?,
            v => {
                write!(f, "#c")?;
                fmt_value(f, v)?;
            }
        }

        match &self.implicit_hydrogens {
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

        fmt_value_field(f, "#n", &self.lone_pairs)?;
        let (u_field, m_field) = self.spin.to_pair();
        fmt_value_field(f, "#u", &u_field)?;
        fmt_multiplicity(f, &m_field, &u_field)?;

        Ok(())
    }
}

impl Display for AtomPattern {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        Display::fmt(&self.ast, f)?;
        for c in &self.constraints {
            fmt_constraint(f, c)?;
        }
        Ok(())
    }
}

fn fmt_constraint(f: &mut fmt::Formatter<'_>, c: &AtomConstraint) -> fmt::Result {
    match c {
        AtomConstraint::Valence(v) => fmt_value_field_required(f, "#v", v),
        AtomConstraint::DonatedPairs(v) => fmt_value_field_required(f, "#d", v),
        AtomConstraint::AcceptedPairs(v) => fmt_value_field_required(f, "#r", v),
        AtomConstraint::MulticenterValence(v) => fmt_value_field_required(f, "#m", v),
        AtomConstraint::AromaticValence(c) => match c {
            AromaticValenceConstraint::NotAromatic => write!(f, "#a!"),
            AromaticValenceConstraint::Value(ValueAst::Lit(1)) => write!(f, "#a"),
            AromaticValenceConstraint::Value(ValueAst::Lit(n)) => write!(f, "#a{}", n),
            AromaticValenceConstraint::Value(v) => {
                write!(f, "#a")?;
                fmt_value(f, v)
            }
        },
        AtomConstraint::Degree(v) => fmt_value_field_required(f, "#D", v),
        AtomConstraint::Connectivity(v) => fmt_value_field_required(f, "#X", v),
        AtomConstraint::TotalHCount(v) => fmt_value_field_required(f, "#H", v),
        AtomConstraint::InRing => write!(f, "#R"),
        AtomConstraint::RingCount(v) => fmt_ring_count(f, v),
        AtomConstraint::RingSize(_) => Ok(()),
    }
}

impl<'de> FromEdn<'de> for AtomAst {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => AtomAst::from_str(s).map_err(|e| DeError::subgrammar("atom", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for AtomAst {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

impl<'de> FromEdn<'de> for AtomPattern {
    fn from_edn(edn: &Edn<'de>) -> Result<Self, DeError> {
        match edn {
            Edn::Str(s) => parse_atom_dsl(s).map_err(|e| DeError::subgrammar("atom", e)),
            other => Err(DeError::TypeMismatch {
                expected: "string",
                got: other.kind(),
                path: Vec::new(),
            }),
        }
    }
}

impl ToEdn for AtomPattern {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

/// Parse a complete atom-string into an `AtomPattern` (base AST + lifted constraints).
pub fn parse_atom_dsl(input: &str) -> Result<AtomPattern, AtomDslError> {
    all_consuming(atom_dsl)
        .parse(input)
        .map(|(_, r)| r)
        .map_err(|e| match e {
            Err::Error(e) | Err::Failure(e) => e,
            Err::Incomplete(_) => AtomDslError::Incomplete,
        })
}

/// Atom-string parser (does not require consuming all input)
pub fn atom_dsl(i: &str) -> IResult<&str, AtomPattern, AtomDslError> {
    let (remaining, (element, preds)) = pair(
        delimited(multispace0, element_expr, multispace0),
        many0(terminated(atom_predicate, multispace0)),
    )
    .parse(i)?;

    let mut pattern = AtomPattern::new(AtomAst::new(element));
    apply_predicates(&mut pattern, preds).map_err(Err::Error)?;
    Ok((remaining, pattern))
}

fn is_set(v: &ValueAst) -> bool {
    !matches!(v, ValueAst::Undetermined)
}

fn apply_predicates(
    pattern: &mut AtomPattern,
    preds: Vec<AtomPredicate>,
) -> Result<(), AtomDslError> {
    let ast = &mut pattern.ast;
    for pred in preds {
        match pred {
            AtomPredicate::IsotopeMass(v) => {
                if !matches!(ast.isotope_mass, IsotopeAst::Undetermined) {
                    return Err(AtomDslError::DuplicateAtomPredicate("#i".to_string()));
                }
                ast.isotope_mass = v;
            }
            AtomPredicate::Charge(v) => {
                if is_set(&ast.charge) {
                    return Err(AtomDslError::DuplicateAtomPredicate("#c".to_string()));
                }
                ast.charge = v;
            }
            AtomPredicate::ImplicitHydrogens(v) => {
                if !matches!(ast.implicit_hydrogens, ImplicitHydrogensAst::Undetermined) {
                    return Err(AtomDslError::DuplicateAtomPredicate("#h".to_string()));
                }
                ast.implicit_hydrogens = v;
            }
            AtomPredicate::LonePairs(v) => {
                if is_set(&ast.lone_pairs) {
                    return Err(AtomDslError::DuplicateAtomPredicate("#n".to_string()));
                }
                ast.lone_pairs = v;
            }
            AtomPredicate::UnpairedElectrons(v) => {
                let SpinStateAst { unpaired, .. } = &mut ast.spin else {
                    unreachable!("default is Pair")
                };
                if !matches!(unpaired, ValueAst::Undetermined) {
                    return Err(AtomDslError::DuplicateAtomPredicate("#u".to_string()));
                }
                *unpaired = v;
            }
            AtomPredicate::Multiplicity(v) => {
                let SpinStateAst { multiplicity, .. } = &mut ast.spin else {
                    unreachable!("default is Pair")
                };
                if !matches!(multiplicity, ValueAst::Undetermined) {
                    return Err(AtomDslError::DuplicateAtomPredicate("#s".to_string()));
                }
                *multiplicity = v;
            }
            AtomPredicate::Constraint(c) => {
                let tag = constraint_tag(&c);
                if pattern.constraints.iter().any(|existing| constraint_tag(existing) == tag) {
                    return Err(AtomDslError::DuplicateAtomPredicate(tag.to_string()));
                }
                pattern.constraints.push(c);
            }
        }
    }
    Ok(())
}

fn constraint_tag(c: &AtomConstraint) -> &'static str {
    match c {
        AtomConstraint::Valence(_) => "#v",
        AtomConstraint::DonatedPairs(_) => "#d",
        AtomConstraint::AcceptedPairs(_) => "#r",
        AtomConstraint::AromaticValence(_) => "#a",
        AtomConstraint::MulticenterValence(_) => "#m",
        AtomConstraint::Degree(_) => "#D",
        AtomConstraint::Connectivity(_) => "#X",
        AtomConstraint::TotalHCount(_) => "#H",
        AtomConstraint::InRing => "#R",
        AtomConstraint::RingCount(_) => "#R",
        AtomConstraint::RingSize(_) => "#rs",
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AtomPredicate {
    IsotopeMass(IsotopeAst),
    Charge(ValueAst),
    ImplicitHydrogens(ImplicitHydrogensAst),
    LonePairs(ValueAst),
    UnpairedElectrons(ValueAst),
    Multiplicity(ValueAst),
    Constraint(AtomConstraint),
}

fn atom_predicate(i: &str) -> IResult<&str, AtomPredicate, AtomDslError> {
    let (remaining, prefix) = take(2usize)(i)?;
    match prefix {
        "#i" => map(isotope_expr, AtomPredicate::IsotopeMass).parse(remaining),
        "#c" => map(charge_value, AtomPredicate::Charge).parse(remaining),
        "#h" => map(hydrogen_expr, AtomPredicate::ImplicitHydrogens).parse(remaining),
        "#n" => map(optional_value, AtomPredicate::LonePairs).parse(remaining),
        "#u" => map(optional_value, AtomPredicate::UnpairedElectrons).parse(remaining),
        "#s" => map(optional_value, AtomPredicate::Multiplicity).parse(remaining),
        "#v" => map(
            optional_value,
            |v| AtomPredicate::Constraint(AtomConstraint::Valence(v)),
        )
        .parse(remaining),
        "#d" => map(
            optional_value,
            |v| AtomPredicate::Constraint(AtomConstraint::DonatedPairs(v)),
        )
        .parse(remaining),
        "#r" => map(
            optional_value,
            |v| AtomPredicate::Constraint(AtomConstraint::AcceptedPairs(v)),
        )
        .parse(remaining),
        "#a" => map(aromatic_valence_expr, |c| {
            AtomPredicate::Constraint(AtomConstraint::AromaticValence(c))
        })
        .parse(remaining),
        "#m" => map(
            optional_value,
            |v| AtomPredicate::Constraint(AtomConstraint::MulticenterValence(v)),
        )
        .parse(remaining),
        "#D" => map(
            optional_value,
            |v| AtomPredicate::Constraint(AtomConstraint::Degree(v)),
        )
        .parse(remaining),
        "#X" => map(
            optional_value,
            |v| AtomPredicate::Constraint(AtomConstraint::Connectivity(v)),
        )
        .parse(remaining),
        "#H" => map(
            optional_value,
            |v| AtomPredicate::Constraint(AtomConstraint::TotalHCount(v)),
        )
        .parse(remaining),
        "#R" => ring_predicate(remaining),
        p if p.starts_with("#") => Err(Err::Failure(AtomDslError::UnknownAtomPredicate(
            p.to_string(),
        ))),
        _ => Err(Err::Failure(AtomDslError::TrailingInput(i.to_string()))),
    }
}

fn ring_predicate(i: &str) -> IResult<&str, AtomPredicate, AtomDslError> {
    preceded(
        multispace0,
        alt((
            map(value_dsl, |v| {
                AtomPredicate::Constraint(AtomConstraint::RingCount(v))
            }),
            success(AtomPredicate::Constraint(AtomConstraint::InRing)),
        )),
    )
    .parse(i)
    .map_err(|_| Err::Error(AtomDslError::InvalidValue(i.to_string())))
}

fn element_expr(i: &str) -> IResult<&str, ElementAst, AtomDslError> {
    alt((
        value(ElementAst::Undetermined, char('*')),
        map(element_set, ElementAst::Set),
        map(element_bind, |(id, set)| ElementAst::Bind { id, set }),
        map(element_ref, ElementAst::Ref),
        map(element_literal, ElementAst::Lit),
    ))
    .parse(i)
    .map_err(|_| Err::Error(AtomDslError::InvalidElement(i.to_string())))
}

fn element_literal(i: &str) -> IResult<&str, Element, NomError<&str>> {
    let (rest, sym) = recognize(pair(
        satisfy(|c: char| c.is_ascii_uppercase()),
        many0(satisfy(|c: char| c.is_ascii_lowercase())),
    ))
    .parse(i)?;
    match Element::from_symbol(sym) {
        Some(el) => Ok((rest, el)),
        None => Err(Err::Error(NomError::new(sym, ErrorKind::Verify))),
    }
}

fn element_set(i: &str) -> IResult<&str, Vec<Element>, NomError<&str>> {
    delimited(
        char('{'),
        delimited(
            multispace0,
            separated_list1(op_char(','), element_literal),
            multispace0,
        ),
        char('}'),
    )
    .parse(i)
}

fn element_bind(i: &str) -> IResult<&str, (String, Vec<Element>), NomError<&str>> {
    delimited(
        char('('),
        pair(
            delimited(multispace0, preceded(char('?'), parse_id), multispace0),
            preceded(
                pair(tag("::"), multispace0),
                terminated(element_set, multispace0),
            ),
        ),
        char(')'),
    )
    .parse(i)
}

fn element_ref(i: &str) -> IResult<&str, String, NomError<&str>> {
    delimited(
        char('('),
        delimited(multispace0, preceded(char('?'), parse_id), multispace0),
        char(')'),
    )
    .parse(i)
}

fn isotope_expr(i: &str) -> IResult<&str, IsotopeAst, AtomDslError> {
    preceded(
        multispace0,
        alt((
            value(IsotopeAst::Natural, char('=')),
            value(IsotopeAst::Undetermined, char('*')),
            map(isotope_set, IsotopeAst::Set),
            map(isotope_bind, |(id, set)| IsotopeAst::Bind { id, set }),
            map(isotope_ref, IsotopeAst::Ref),
            map(nom_u32, IsotopeAst::Lit),
        )),
    )
    .parse(i)
    .map_err(|_| Err::Error(AtomDslError::InvalidIsotope(i.to_string())))
}

fn isotope_set(i: &str) -> IResult<&str, Vec<u32>, NomError<&str>> {
    delimited(
        char('{'),
        delimited(
            multispace0,
            separated_list1(op_char(','), nom_u32),
            multispace0,
        ),
        char('}'),
    )
    .parse(i)
}

fn isotope_bind(i: &str) -> IResult<&str, (String, Vec<u32>), NomError<&str>> {
    delimited(
        char('('),
        pair(
            delimited(multispace0, preceded(char('?'), parse_id), multispace0),
            preceded(
                pair(tag("::"), multispace0),
                terminated(isotope_set, multispace0),
            ),
        ),
        char(')'),
    )
    .parse(i)
}

fn isotope_ref(i: &str) -> IResult<&str, String, NomError<&str>> {
    delimited(
        char('('),
        delimited(multispace0, preceded(char('?'), parse_id), multispace0),
        char(')'),
    )
    .parse(i)
}

fn charge_value(i: &str) -> IResult<&str, ValueAst, AtomDslError> {
    preceded(
        multispace0,
        alt((
            value_dsl,
            value(ValueAst::Lit(1), tag("+")),
            value(ValueAst::Lit(-1), tag("-")),
        )),
    )
    .parse(i)
    .map_err(|_| Err::Error(AtomDslError::InvalidCharge(i.to_string())))
}

fn hydrogen_expr(i: &str) -> IResult<&str, ImplicitHydrogensAst, AtomDslError> {
    preceded(
        multispace0,
        alt((
            value(ImplicitHydrogensAst::Normal, tag("=")),
            map(value_dsl, ImplicitHydrogensAst::Value),
            success(ImplicitHydrogensAst::Value(ValueAst::Lit(1))),
        )),
    )
    .parse(i)
    .map_err(|_| Err::Error(AtomDslError::InvalidImplicitHydrogens(i.to_string())))
}

fn aromatic_valence_expr(i: &str) -> IResult<&str, AromaticValenceConstraint, AtomDslError> {
    preceded(
        multispace0,
        alt((
            value(AromaticValenceConstraint::NotAromatic, tag("!")),
            value(
                AromaticValenceConstraint::Value(ValueAst::Undetermined),
                tag("?"),
            ),
            map(value_dsl, AromaticValenceConstraint::Value),
            success(AromaticValenceConstraint::Value(ValueAst::Lit(1))),
        )),
    )
    .parse(i)
    .map_err(|_| Err::Error(AtomDslError::InvalidValue(i.to_string())))
}

fn optional_value(i: &str) -> IResult<&str, ValueAst, AtomDslError> {
    preceded(multispace0, alt((value_dsl, success(ValueAst::Lit(1)))))
        .parse(i)
        .map_err(|_| Err::Error(AtomDslError::InvalidValue(i.to_string())))
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

/// Format RingCount: always emits the digit, since bare `#R` parses as `InRing`.
fn fmt_ring_count(f: &mut fmt::Formatter<'_>, v: &ValueAst) -> fmt::Result {
    match v {
        ValueAst::Undetermined => write!(f, "#R*"),
        ValueAst::Lit(n) => write!(f, "#R{}", n),
        v => {
            write!(f, "#R")?;
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

/// Suppress multiplicity when it equals unpaired_electrons + 1 (derivable default).
fn fmt_multiplicity(
    f: &mut fmt::Formatter<'_>,
    multiplicity: &ValueAst,
    unpaired: &ValueAst,
) -> fmt::Result {
    let m = match multiplicity {
        ValueAst::Undetermined => return Ok(()),
        ValueAst::Lit(m) => *m,
        v => {
            write!(f, "#s")?;
            return fmt_value(f, v);
        }
    };
    let u: i64 = match unpaired {
        ValueAst::Lit(u) => *u,
        ValueAst::Undetermined => 0,
        _ => -1,
    };
    if m == u + 1 {
        Ok(())
    } else if m == 1 {
        write!(f, "#s")
    } else {
        write!(f, "#s{}", m)
    }
}

fn fmt_value(f: &mut fmt::Formatter<'_>, v: &ValueAst) -> fmt::Result {
    match v {
        ValueAst::Undetermined => write!(f, "*"),
        ValueAst::Lit(n) => write!(f, "{}", n),
        ValueAst::LitSet(s) => {
            write!(f, "{{")?;
            for (i, n) in s.iter().enumerate() {
                if i > 0 {
                    write!(f, ",")?;
                }
                write!(f, "{}", n)?;
            }
            write!(f, "}}")
        }
        ValueAst::Expr(_) => write!(f, "<expr>"),
    }
}

#[cfg(test)]
mod tests {
    use nom::Err;
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::element::Element;
    use umol_ast::ast::value::{Expr, RelOp, ValueAst};

    use super::*;

    fn pat(ast: AtomAst) -> AtomPattern {
        AtomPattern::new(ast)
    }

    fn pat_with(ast: AtomAst, constraints: Vec<AtomConstraint>) -> AtomPattern {
        AtomPattern::with_constraints(ast, constraints)
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::carbon("C", pat(AtomAst::new(ElementAst::Lit(Element::C))))]
    #[case::iron("Fe", pat(AtomAst::new(ElementAst::Lit(Element::Fe))))]
    #[case::chlorine("Cl", pat(AtomAst::new(ElementAst::Lit(Element::Cl))))]
    #[case::whitespace("  C  ", pat(AtomAst::new(ElementAst::Lit(Element::C))))]
    #[case::undetermined("*", pat(AtomAst::new(ElementAst::Undetermined)))]
    #[case::element_set("{C,N,O}", pat(AtomAst::new(ElementAst::Set(vec![Element::C, Element::N, Element::O]))))]
    #[case::element_bind("(?e :: {C,N})", pat(AtomAst::new(ElementAst::Bind { id: "e".to_string(), set: vec![Element::C, Element::N] })))]
    #[case::element_ref("(?e)", pat(AtomAst::new(ElementAst::Ref("e".to_string()))))]
    #[case::isotope("C#i12", pat(AtomAst { isotope_mass: IsotopeAst::Lit(12), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::isotope_natural("C#i=", pat(AtomAst { isotope_mass: IsotopeAst::Natural, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_pos("C#c+2", pat(AtomAst { charge: ValueAst::Lit(2), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_neg("C#c-2", pat(AtomAst { charge: ValueAst::Lit(-2), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_plus("C#c+", pat(AtomAst { charge: ValueAst::Lit(1), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_minus("C#c-", pat(AtomAst { charge: ValueAst::Lit(-1), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::charge_zero("C#c0", pat(AtomAst { charge: ValueAst::Lit(0), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_count("C#h3", pat(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Value(ValueAst::Lit(3)), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_normal("C#h=", pat(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Normal, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_undetermined("C#h*", pat(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Value(ValueAst::Undetermined), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_bind("C#h(?h)", pat(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Value(ValueAst::Expr(Expr::Var("h".to_string()))), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_set("N#h?h :: {2,3}", pat(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Value(ValueAst::Expr(Expr::Mem(Box::new(Expr::Var("h".to_string())), vec![2, 3]))), ..AtomAst::new(ElementAst::Lit(Element::N)) }))]
    #[case::h_expr("C#h?h >= 1", pat(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Value(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1))))), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::h_omit("C#h", pat(AtomAst { implicit_hydrogens: ImplicitHydrogensAst::Value(ValueAst::Lit(1)), ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::lone_pairs("O#n2", pat(AtomAst { lone_pairs: ValueAst::Lit(2), ..AtomAst::new(ElementAst::Lit(Element::O)) }))]
    #[case::lone_pairs_omit("O#n", pat(AtomAst { lone_pairs: ValueAst::Lit(1), ..AtomAst::new(ElementAst::Lit(Element::O)) }))]
    #[case::unpaired("C#u2", pat(AtomAst { spin: SpinStateAst { unpaired: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined }, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::unpaired_omit("C#u", pat(AtomAst { spin: SpinStateAst { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Undetermined }, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multiplicity("C#s3", pat(AtomAst { spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(3) }, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::multiplicity_omit("C#s", pat(AtomAst { spin: SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(1) }, ..AtomAst::new(ElementAst::Lit(Element::C)) }))]
    #[case::valence("C#v4", pat_with(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::Valence(ValueAst::Lit(4))]))]
    #[case::donated_pairs("N#d1", pat_with(AtomAst::new(ElementAst::Lit(Element::N)), vec![AtomConstraint::DonatedPairs(ValueAst::Lit(1))]))]
    #[case::accepted_pairs("B#r1", pat_with(AtomAst::new(ElementAst::Lit(Element::B)), vec![AtomConstraint::AcceptedPairs(ValueAst::Lit(1))]))]
    #[case::arom_not_aromatic("C#a!", pat_with(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::AromaticValence(AromaticValenceConstraint::NotAromatic)]))]
    #[case::arom_undetermined("C#a*", pat_with(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(ValueAst::Undetermined))]))]
    #[case::arom_zero("C#a0", pat_with(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(ValueAst::Lit(0)))]))]
    #[case::arom_one("C#a1", pat_with(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(ValueAst::Lit(1)))]))]
    #[case::arom_omit("C#a", pat_with(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(ValueAst::Lit(1)))]))]
    #[case::multicenter("C#m2", pat_with(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::MulticenterValence(ValueAst::Lit(2))]))]
    #[case::degree("C#D2", pat_with(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::Degree(ValueAst::Lit(2))]))]
    #[case::connectivity("C#X3", pat_with(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::Connectivity(ValueAst::Lit(3))]))]
    #[case::total_h_count("C#H1", pat_with(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::TotalHCount(ValueAst::Lit(1))]))]
    #[case::in_ring_bare("C#R", pat_with(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::InRing]))]
    #[case::ring_count("C#R2", pat_with(AtomAst::new(ElementAst::Lit(Element::C)), vec![AtomConstraint::RingCount(ValueAst::Lit(2))]))]
    fn test_parse_atom_dsl(#[case] input: &str, #[case] expected: AtomPattern) {
        let result = atom_dsl(input);
        assert!(result.is_ok(), "{:?} should succeed, got {:?}", input, result.unwrap_err());
        let (remaining, ast) = result.unwrap();
        assert!(remaining.is_empty(), "{:?} should consume all input, remaining: {:?}", input, remaining);
        assert_eq!(ast, expected);
    }

    #[rstest]
    #[case::empty("", AtomDslError::InvalidElement("".to_string()))]
    #[case::no_element("#h3", AtomDslError::InvalidElement("#h3".to_string()))]
    #[case::unknown_pred("C#x", AtomDslError::UnknownAtomPredicate("#x".to_string()))]
    #[case::dup_h("C#h3#h2", AtomDslError::DuplicateAtomPredicate("#h".to_string()))]
    #[case::dup_charge("C#c+#c-", AtomDslError::DuplicateAtomPredicate("#c".to_string()))]
    #[case::dup_valence("C#v3#v4", AtomDslError::DuplicateAtomPredicate("#v".to_string()))]
    #[case::trailing("C#h3 foo", AtomDslError::TrailingInput("foo".to_string()))]
    fn test_parse_atom_dsl_invalid(#[case] input: &str, #[case] expected: AtomDslError) {
        let result = atom_dsl(input);
        assert!(
            result.is_err(),
            "{:?} should fail, got {:?}",
            input,
            result.unwrap()
        );
        let err = match result.unwrap_err() {
            Err::Error(e) | Err::Failure(e) => e,
            Err::Incomplete(_) => AtomDslError::Incomplete,
        };
        assert_eq!(
            err, expected,
            "{:?} should fail with {:?}, got {:?}",
            input, expected, err
        );
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
    fn test_element_expr(#[case] input: &str, #[case] expected: ElementAst) {
        let result = element_expr(input);
        assert!(result.is_ok(), "{input:?} should succeed, got {:?}", result.unwrap_err());
        let (remaining, expr) = result.unwrap();
        assert!(remaining.is_empty(), "{input:?} should consume all input, remaining: {remaining:?}");
        assert_eq!(expr, expected);
    }

    #[rstest]
    #[case::empty("")]
    #[case::lowercase("c")]
    #[case::invalid("123")]
    #[case::unknown_element("Xx")]
    fn test_element_expr_invalid(#[case] input: &str) {
        let result = element_expr(input);
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
    #[case::set("{12,13,14}", IsotopeAst::Set(vec![12, 13, 14]))]
    #[case::bind("(?m :: {12,13})", IsotopeAst::Bind { id: "m".to_string(), set: vec![12, 13] })]
    #[case::ref_("(?m)", IsotopeAst::Ref("m".to_string()))]
    fn test_isotope_expr(#[case] input: &str, #[case] expected: IsotopeAst) {
        let result = isotope_expr(input);
        assert!(result.is_ok(), "{input:?} should succeed, got {:?}", result.unwrap_err());
        let (remaining, expr) = result.unwrap();
        assert!(remaining.is_empty(), "{input:?} should consume all input, remaining: {remaining:?}");
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
    #[case::h_count("#h3", AtomPredicate::ImplicitHydrogens(ImplicitHydrogensAst::Value(ValueAst::Lit(3))))]
    #[case::h_normal("#h=", AtomPredicate::ImplicitHydrogens(ImplicitHydrogensAst::Normal))]
    #[case::h_undetermined("#h*", AtomPredicate::ImplicitHydrogens(ImplicitHydrogensAst::Value(ValueAst::Undetermined)))]
    #[case::h_omit("#h", AtomPredicate::ImplicitHydrogens(ImplicitHydrogensAst::Value(ValueAst::Lit(1))))]
    #[case::lone_pairs("#n2", AtomPredicate::LonePairs(ValueAst::Lit(2)))]
    #[case::lone_pairs_omit("#n", AtomPredicate::LonePairs(ValueAst::Lit(1)))]
    #[case::unpaired("#u2", AtomPredicate::UnpairedElectrons(ValueAst::Lit(2)))]
    #[case::unpaired_omit("#u", AtomPredicate::UnpairedElectrons(ValueAst::Lit(1)))]
    #[case::multiplicity("#s3", AtomPredicate::Multiplicity(ValueAst::Lit(3)))]
    #[case::multiplicity_omit("#s", AtomPredicate::Multiplicity(ValueAst::Lit(1)))]
    #[case::valence("#v4", AtomPredicate::Constraint(AtomConstraint::Valence(ValueAst::Lit(4))))]
    #[case::donated_pairs("#d1", AtomPredicate::Constraint(AtomConstraint::DonatedPairs(ValueAst::Lit(1))))]
    #[case::accepted_pairs("#r1", AtomPredicate::Constraint(AtomConstraint::AcceptedPairs(ValueAst::Lit(1))))]
    #[case::arom_not_aromatic("#a!", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceConstraint::NotAromatic)))]
    #[case::arom_undetermined("#a*", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(ValueAst::Undetermined))))]
    #[case::arom_lit("#a2", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(ValueAst::Lit(2)))))]
    #[case::arom_omit("#a", AtomPredicate::Constraint(AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(ValueAst::Lit(1)))))]
    #[case::multicenter("#m2", AtomPredicate::Constraint(AtomConstraint::MulticenterValence(ValueAst::Lit(2))))]
    #[case::degree("#D2", AtomPredicate::Constraint(AtomConstraint::Degree(ValueAst::Lit(2))))]
    #[case::degree_omit("#D", AtomPredicate::Constraint(AtomConstraint::Degree(ValueAst::Lit(1))))]
    #[case::connectivity("#X3", AtomPredicate::Constraint(AtomConstraint::Connectivity(ValueAst::Lit(3))))]
    #[case::total_h_count("#H1", AtomPredicate::Constraint(AtomConstraint::TotalHCount(ValueAst::Lit(1))))]
    #[case::in_ring_bare("#R", AtomPredicate::Constraint(AtomConstraint::InRing))]
    #[case::ring_count("#R2", AtomPredicate::Constraint(AtomConstraint::RingCount(ValueAst::Lit(2))))]
    fn test_atom_predicate(#[case] input: &str, #[case] expected: AtomPredicate) {
        let result = atom_predicate(input);
        assert!(result.is_ok(), "{input:?} should succeed, got {:?}", result.unwrap_err());
        let (_, pred) = result.unwrap();
        assert_eq!(pred, expected);
    }

    #[rstest]
    #[case::unknown("#x", AtomDslError::UnknownAtomPredicate("#x".to_string()))]
    #[case::unknown_tag("#z", AtomDslError::UnknownAtomPredicate("#z".to_string()))]
    #[case::trailing_no_hash("fo", AtomDslError::TrailingInput("fo".to_string()))]
    fn test_atom_predicate_error(#[case] input: &str, #[case] expected: AtomDslError) {
        let result = atom_predicate(input);
        assert!(result.is_err(), "{input:?} should fail, got {:?}", result.unwrap());
        let err = match result.unwrap_err() {
            Err::Error(e) | Err::Failure(e) => e,
            Err::Incomplete(_) => AtomDslError::Incomplete,
        };
        assert_eq!(err, expected);
    }
}
