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
use umol_shared::atom_ast::{AromaticAst, ElementAst, HydrogenAst, IsotopeAst};
use umol_shared::element::Element;
use umol_shared::spin_ast::SpinStateAst;
use umol_shared::value_ast::ValueAst;
use umol_edn::{DeError, Edn, FromEdn, ToEdn};

use super::error::AtomDslError;
use super::value::{op_char, parse_id, value_dsl};
use crate::ast::atom::AtomAst;

impl FromStr for AtomAst {
    type Err = AtomDslError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        parse_atom_dsl(s)
    }
}

impl Display for AtomAst {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fmt_element(f, &self.element)?;

        match &self.isotope_mass {
            None => {}
            Some(IsotopeAst::Natural) => write!(f, "#i=")?,
            Some(IsotopeAst::Lit(n)) => write!(f, "#i{}", n)?,
            Some(IsotopeAst::Wildcard) => write!(f, "#i*")?,
            Some(IsotopeAst::Set(ns)) => {
                write!(f, "#i{{")?;
                for (i, n) in ns.iter().enumerate() {
                    if i > 0 {
                        write!(f, ",")?;
                    }
                    write!(f, "{}", n)?;
                }
                write!(f, "}}")?;
            }
            Some(IsotopeAst::Bind { id, set }) => {
                write!(f, "#i(?{} :: {{", id)?;
                for (i, n) in set.iter().enumerate() {
                    if i > 0 {
                        write!(f, ",")?;
                    }
                    write!(f, "{}", n)?;
                }
                write!(f, "}})")?;
            }
            Some(IsotopeAst::Ref(id)) => write!(f, "#i(?{})", id)?,
        }

        match &self.charge {
            None | Some(ValueAst::Lit(0)) => {}
            Some(ValueAst::Lit(1)) => write!(f, "#c+")?,
            Some(ValueAst::Lit(-1)) => write!(f, "#c-")?,
            Some(ValueAst::Lit(n)) if *n > 0 => write!(f, "#c+{}", n)?,
            Some(ValueAst::Lit(n)) => write!(f, "#c{}", n)?,
            Some(ValueAst::Wildcard) => write!(f, "#c*")?,
            Some(v) => {
                write!(f, "#c")?;
                fmt_value(f, v)?;
            }
        }

        match &self.implicit_hydrogens {
            None | Some(HydrogenAst::Value(ValueAst::Lit(0))) => {}
            Some(HydrogenAst::Normal) => write!(f, "#h=")?,
            Some(HydrogenAst::Value(ValueAst::Lit(1))) => write!(f, "#h")?,
            Some(HydrogenAst::Value(ValueAst::Lit(n))) => write!(f, "#h{}", n)?,
            Some(HydrogenAst::Value(ValueAst::Wildcard)) => write!(f, "#h*")?,
            Some(HydrogenAst::Value(v)) => {
                write!(f, "#h")?;
                fmt_value(f, v)?;
            }
        }

        fmt_unsigned_field(f, "#n", &self.lone_pairs)?;
        let (u_field, m_field) = self
            .spin
            .as_ref()
            .map(SpinStateAst::to_pair)
            .unwrap_or((None, None));
        fmt_unsigned_field(f, "#u", &u_field)?;
        fmt_multiplicity(f, &m_field, &u_field)?;
        fmt_unsigned_field(f, "#v", &self.valence)?;
        fmt_unsigned_field(f, "#d", &self.donated_pairs)?;
        fmt_unsigned_field(f, "#r", &self.accepted_pairs)?;

        match &self.aromatic_valence {
            None => {}
            Some(AromaticAst::Unspecified) => write!(f, "#a?")?,
            Some(AromaticAst::NotAromatic) => write!(f, "#a!")?,
            Some(AromaticAst::Value(ValueAst::Lit(1))) => write!(f, "#a")?,
            Some(AromaticAst::Value(ValueAst::Lit(n))) => write!(f, "#a{}", n)?,
            Some(AromaticAst::Value(v)) => {
                write!(f, "#a")?;
                fmt_value(f, v)?;
            }
        }

        fmt_unsigned_field(f, "#m", &self.multicenter_valence)?;

        Ok(())
    }
}

impl<'de> FromEdn<'de> for AtomAst {
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

impl ToEdn for AtomAst {
    fn to_edn(&self) -> Edn<'static> {
        Edn::Str(Cow::Owned(self.to_string()))
    }
}

/// Parse a complete atom-string
pub fn parse_atom_dsl(input: &str) -> Result<AtomAst, AtomDslError> {
    all_consuming(atom_dsl)
        .parse(input)
        .map(|(_, r)| r)
        .map_err(|e| match e {
            Err::Error(e) | Err::Failure(e) => e,
            Err::Incomplete(_) => AtomDslError::Incomplete,
        })
}

/// Atom-string parser (does not require consuming all input)
pub fn atom_dsl(i: &str) -> IResult<&str, AtomAst, AtomDslError> {
    let (remaining, (element, preds)) = pair(
        delimited(multispace0, element_expr, multispace0),
        many0(terminated(atom_predicate, multispace0)),
    )
    .parse(i)?;

    let mut ast = AtomAst::new(element);
    update_atom_ast(&mut ast, preds).map_err(Err::Error)?;
    Ok((remaining, ast))
}

fn update_atom_ast(ast: &mut AtomAst, preds: Vec<AtomPredicate>) -> Result<(), AtomDslError> {
    for pred in preds {
        match pred {
            AtomPredicate::IsotopeMass(v) => {
                if ast.isotope_mass.is_some() {
                    return Err(AtomDslError::DuplicateAtomPredicate("#i".to_string()));
                }
                ast.isotope_mass = Some(v);
            }
            AtomPredicate::Charge(v) => {
                if ast.charge.is_some() {
                    return Err(AtomDslError::DuplicateAtomPredicate("#c".to_string()));
                }
                ast.charge = Some(v);
            }
            AtomPredicate::ImplicitHydrogens(v) => {
                if ast.implicit_hydrogens.is_some() {
                    return Err(AtomDslError::DuplicateAtomPredicate("#h".to_string()));
                }
                ast.implicit_hydrogens = Some(v);
            }
            AtomPredicate::LonePairs(v) => {
                if ast.lone_pairs.is_some() {
                    return Err(AtomDslError::DuplicateAtomPredicate("#n".to_string()));
                }
                ast.lone_pairs = Some(v);
            }
            AtomPredicate::UnpairedElectrons(v) => {
                let pair = ast.spin.get_or_insert(SpinStateAst::Pair {
                    unpaired: None,
                    multiplicity: None,
                });
                let SpinStateAst::Pair { unpaired, .. } = pair else {
                    unreachable!("parser only constructs Pair")
                };
                if unpaired.is_some() {
                    return Err(AtomDslError::DuplicateAtomPredicate("#u".to_string()));
                }
                *unpaired = Some(v);
            }
            AtomPredicate::Multiplicity(v) => {
                let pair = ast.spin.get_or_insert(SpinStateAst::Pair {
                    unpaired: None,
                    multiplicity: None,
                });
                let SpinStateAst::Pair { multiplicity, .. } = pair else {
                    unreachable!("parser only constructs Pair")
                };
                if multiplicity.is_some() {
                    return Err(AtomDslError::DuplicateAtomPredicate("#s".to_string()));
                }
                *multiplicity = Some(v);
            }
            AtomPredicate::Valence(v) => {
                if ast.valence.is_some() {
                    return Err(AtomDslError::DuplicateAtomPredicate("#v".to_string()));
                }
                ast.valence = Some(v);
            }
            AtomPredicate::DonatedPairs(v) => {
                if ast.donated_pairs.is_some() {
                    return Err(AtomDslError::DuplicateAtomPredicate("#d".to_string()));
                }
                ast.donated_pairs = Some(v);
            }
            AtomPredicate::AcceptedPairs(v) => {
                if ast.accepted_pairs.is_some() {
                    return Err(AtomDslError::DuplicateAtomPredicate("#r".to_string()));
                }
                ast.accepted_pairs = Some(v);
            }
            AtomPredicate::AromaticValence(v) => {
                if ast.aromatic_valence.is_some() {
                    return Err(AtomDslError::DuplicateAtomPredicate("#a".to_string()));
                }
                ast.aromatic_valence = Some(v);
            }
            AtomPredicate::MulticenterValence(v) => {
                if ast.multicenter_valence.is_some() {
                    return Err(AtomDslError::DuplicateAtomPredicate("#m".to_string()));
                }
                ast.multicenter_valence = Some(v);
            }
        }
    }
    Ok(())
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AtomPredicate {
    IsotopeMass(IsotopeAst),
    Charge(ValueAst),
    ImplicitHydrogens(HydrogenAst),
    LonePairs(ValueAst),
    UnpairedElectrons(ValueAst),
    Multiplicity(ValueAst),
    Valence(ValueAst),
    DonatedPairs(ValueAst),
    AcceptedPairs(ValueAst),
    AromaticValence(AromaticAst),
    MulticenterValence(ValueAst),
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
        "#v" => map(optional_value, AtomPredicate::Valence).parse(remaining),
        "#d" => map(optional_value, AtomPredicate::DonatedPairs).parse(remaining),
        "#r" => map(optional_value, AtomPredicate::AcceptedPairs).parse(remaining),
        "#a" => map(aromatic_valence_expr, AtomPredicate::AromaticValence).parse(remaining),
        "#m" => map(optional_value, AtomPredicate::MulticenterValence).parse(remaining),
        p if p.starts_with("#") => Err(Err::Failure(AtomDslError::UnknownAtomPredicate(
            p.to_string(),
        ))),
        _ => Err(Err::Failure(AtomDslError::TrailingInput(i.to_string()))),
    }
}

fn element_expr(i: &str) -> IResult<&str, ElementAst, AtomDslError> {
    alt((
        value(ElementAst::Wildcard, char('*')),
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
            value(IsotopeAst::Wildcard, char('*')),
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

fn hydrogen_expr(i: &str) -> IResult<&str, HydrogenAst, AtomDslError> {
    preceded(
        multispace0,
        alt((
            value(HydrogenAst::Normal, tag("=")),
            map(value_dsl, HydrogenAst::Value),
            success(HydrogenAst::Value(ValueAst::Lit(1))),
        )),
    )
    .parse(i)
    .map_err(|_| Err::Error(AtomDslError::InvalidImplicitHydrogens(i.to_string())))
}

fn aromatic_valence_expr(i: &str) -> IResult<&str, AromaticAst, AtomDslError> {
    preceded(
        multispace0,
        alt((
            value(AromaticAst::NotAromatic, tag("!")),
            value(AromaticAst::Unspecified, tag("?")),
            map(value_dsl, AromaticAst::Value),
            success(AromaticAst::Value(ValueAst::Lit(1))),
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
        ElementAst::Wildcard => write!(f, "*"),
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

/// Suppress None and Lit(0); abbreviate Lit(1) to just the prefix.
fn fmt_unsigned_field(
    f: &mut fmt::Formatter<'_>,
    prefix: &str,
    v: &Option<ValueAst>,
) -> fmt::Result {
    match v {
        None | Some(ValueAst::Lit(0)) => Ok(()),
        Some(ValueAst::Lit(1)) => write!(f, "{}", prefix),
        Some(ValueAst::Lit(n)) => write!(f, "{}{}", prefix, n),
        Some(ValueAst::Wildcard) => write!(f, "{}*", prefix),
        Some(v) => {
            write!(f, "{}", prefix)?;
            fmt_value(f, v)
        }
    }
}

/// Suppress multiplicity when it equals unpaired_electrons + 1 (derivable default).
fn fmt_multiplicity(
    f: &mut fmt::Formatter<'_>,
    multiplicity: &Option<ValueAst>,
    unpaired: &Option<ValueAst>,
) -> fmt::Result {
    let m = match multiplicity {
        None => return Ok(()),
        Some(ValueAst::Lit(m)) => *m,
        Some(ValueAst::Wildcard) => return write!(f, "#s*"),
        Some(v) => {
            write!(f, "#s")?;
            return fmt_value(f, v);
        }
    };
    let u: i64 = match unpaired {
        Some(ValueAst::Lit(u)) => *u,
        None => 0,
        _ => -1, // non-literal: can't determine derivability, always print
    };
    if m as i64 == u + 1 {
        Ok(()) // derivable from unpaired, suppress
    } else if m == 1 {
        write!(f, "#s")
    } else {
        write!(f, "#s{}", m)
    }
}

fn fmt_value(f: &mut fmt::Formatter<'_>, v: &ValueAst) -> fmt::Result {
    match v {
        ValueAst::Wildcard => write!(f, "*"),
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
    use umol_shared::value_ast::{Expr, RelOp, ValueAst};

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::carbon("C", AtomAst::new(ElementAst::Lit(Element::C)))]
    #[case::iron("Fe", AtomAst::new(ElementAst::Lit(Element::Fe)))]
    #[case::chlorine("Cl", AtomAst::new(ElementAst::Lit(Element::Cl)))]
    #[case::whitespace("  C  ", AtomAst::new(ElementAst::Lit(Element::C)))]
    #[case::wildcard("*", AtomAst::new(ElementAst::Wildcard))]
    #[case::element_set("{C,N,O}", AtomAst::new(ElementAst::Set(vec![Element::C, Element::N, Element::O])))]
    #[case::element_bind("(?e :: {C,N})", AtomAst::new(ElementAst::Bind { id: "e".to_string(), set: vec![Element::C, Element::N] }))]
    #[case::element_ref("(?e)", AtomAst::new(ElementAst::Ref("e".to_string())))]
    #[case::isotope("C#i12", AtomAst { isotope_mass: Some(IsotopeAst::Lit(12)), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::isotope_natural("C#i=", AtomAst { isotope_mass: Some(IsotopeAst::Natural), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::charge_pos("C#c+2", AtomAst { charge: Some(ValueAst::Lit(2)), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::charge_neg("C#c-2", AtomAst { charge: Some(ValueAst::Lit(-2)), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::charge_plus("C#c+", AtomAst { charge: Some(ValueAst::Lit(1)), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::charge_minus("C#c-", AtomAst { charge: Some(ValueAst::Lit(-1)), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::charge_zero("C#c0", AtomAst { charge: Some(ValueAst::Lit(0)), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::h_count("C#h3", AtomAst { implicit_hydrogens: Some(HydrogenAst::Value(ValueAst::Lit(3))), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::h_normal("C#h=", AtomAst { implicit_hydrogens: Some(HydrogenAst::Normal), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::h_wildcard("C#h*", AtomAst { implicit_hydrogens: Some(HydrogenAst::Value(ValueAst::Wildcard)), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::h_bind("C#h(?h)", AtomAst { implicit_hydrogens: Some(HydrogenAst::Value(ValueAst::Expr(Expr::Var("h".to_string())))), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::h_set("N#h?h :: {2,3}", AtomAst { implicit_hydrogens: Some(HydrogenAst::Value(ValueAst::Expr(Expr::Mem(Box::new(Expr::Var("h".to_string())), vec![2, 3])))), ..AtomAst::new(ElementAst::Lit(Element::N)) })]
    #[case::h_expr("C#h?h >= 1", AtomAst { implicit_hydrogens: Some(HydrogenAst::Value(ValueAst::Expr(Expr::Rel(Box::new(Expr::Var("h".to_string())), RelOp::Ge, Box::new(Expr::Lit(1)))))), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::h_omit("C#h", AtomAst { implicit_hydrogens: Some(HydrogenAst::Value(ValueAst::Lit(1))), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::lone_pairs("O#n2", AtomAst { lone_pairs: Some(ValueAst::Lit(2)), ..AtomAst::new(ElementAst::Lit(Element::O)) })]
    #[case::lone_pairs_omit("O#n", AtomAst { lone_pairs: Some(ValueAst::Lit(1)), ..AtomAst::new(ElementAst::Lit(Element::O)) })]
    #[case::unpaired("C#u2", AtomAst { spin: Some(SpinStateAst::Pair { unpaired: Some(ValueAst::Lit(2)), multiplicity: None }), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::unpaired_omit("C#u", AtomAst { spin: Some(SpinStateAst::Pair { unpaired: Some(ValueAst::Lit(1)), multiplicity: None }), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::multiplicity("C#s3", AtomAst { spin: Some(SpinStateAst::Pair { unpaired: None, multiplicity: Some(ValueAst::Lit(3)) }), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::multiplicity_omit("C#s", AtomAst { spin: Some(SpinStateAst::Pair { unpaired: None, multiplicity: Some(ValueAst::Lit(1)) }), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::valence("C#v4", AtomAst { valence: Some(ValueAst::Lit(4)), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::donated_pairs("N#d1", AtomAst { donated_pairs: Some(ValueAst::Lit(1)), ..AtomAst::new(ElementAst::Lit(Element::N)) })]
    #[case::accepted_pairs("B#r1", AtomAst { accepted_pairs: Some(ValueAst::Lit(1)), ..AtomAst::new(ElementAst::Lit(Element::B)) })]
    #[case::arom_unspecified("C#a?", AtomAst { aromatic_valence: Some(AromaticAst::Unspecified), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::arom_not_aromatic("C#a!", AtomAst { aromatic_valence: Some(AromaticAst::NotAromatic), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::arom_aromatic("C#a*", AtomAst { aromatic_valence: Some(AromaticAst::Value(ValueAst::Wildcard)), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::arom_zero("C#a0", AtomAst { aromatic_valence: Some(AromaticAst::Value(ValueAst::Lit(0))), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::arom_one("C#a1", AtomAst { aromatic_valence: Some(AromaticAst::Value(ValueAst::Lit(1))), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::arom_omit("C#a", AtomAst { aromatic_valence: Some(AromaticAst::Value(ValueAst::Lit(1))), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    #[case::multicenter("C#m2", AtomAst { multicenter_valence: Some(ValueAst::Lit(2)), ..AtomAst::new(ElementAst::Lit(Element::C)) })]
    fn test_parse_atom_dsl(#[case] input: &str, #[case] expected: AtomAst) {
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
    #[case::wildcard("*", ElementAst::Wildcard)]
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
    #[case::wildcard("*", IsotopeAst::Wildcard)]
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
    #[case::isotope_wildcard("#i*", AtomPredicate::IsotopeMass(IsotopeAst::Wildcard))]
    #[case::charge_pos("#c+2", AtomPredicate::Charge(ValueAst::Lit(2)))]
    #[case::charge_neg("#c-2", AtomPredicate::Charge(ValueAst::Lit(-2)))]
    #[case::charge_plus("#c+", AtomPredicate::Charge(ValueAst::Lit(1)))]
    #[case::charge_minus("#c-", AtomPredicate::Charge(ValueAst::Lit(-1)))]
    #[case::charge_zero("#c0", AtomPredicate::Charge(ValueAst::Lit(0)))]
    #[case::charge_wildcard("#c*", AtomPredicate::Charge(ValueAst::Wildcard))]
    #[case::h_count("#h3", AtomPredicate::ImplicitHydrogens(HydrogenAst::Value(ValueAst::Lit(3))))]
    #[case::h_normal("#h=", AtomPredicate::ImplicitHydrogens(HydrogenAst::Normal))]
    #[case::h_wildcard("#h*", AtomPredicate::ImplicitHydrogens(HydrogenAst::Value(ValueAst::Wildcard)))]
    #[case::h_omit("#h", AtomPredicate::ImplicitHydrogens(HydrogenAst::Value(ValueAst::Lit(1))))]
    #[case::lone_pairs("#n2", AtomPredicate::LonePairs(ValueAst::Lit(2)))]
    #[case::lone_pairs_omit("#n", AtomPredicate::LonePairs(ValueAst::Lit(1)))]
    #[case::unpaired("#u2", AtomPredicate::UnpairedElectrons(ValueAst::Lit(2)))]
    #[case::unpaired_omit("#u", AtomPredicate::UnpairedElectrons(ValueAst::Lit(1)))]
    #[case::multiplicity("#s3", AtomPredicate::Multiplicity(ValueAst::Lit(3)))]
    #[case::multiplicity_omit("#s", AtomPredicate::Multiplicity(ValueAst::Lit(1)))]
    #[case::valence("#v4", AtomPredicate::Valence(ValueAst::Lit(4)))]
    #[case::donated_pairs("#d1", AtomPredicate::DonatedPairs(ValueAst::Lit(1)))]
    #[case::accepted_pairs("#r1", AtomPredicate::AcceptedPairs(ValueAst::Lit(1)))]
    #[case::arom_unspecified("#a?", AtomPredicate::AromaticValence(AromaticAst::Unspecified))]
    #[case::arom_not_aromatic("#a!", AtomPredicate::AromaticValence(AromaticAst::NotAromatic))]
    #[case::arom_wildcard("#a*", AtomPredicate::AromaticValence(AromaticAst::Value(ValueAst::Wildcard)))]
    #[case::arom_lit("#a2", AtomPredicate::AromaticValence(AromaticAst::Value(ValueAst::Lit(2))))]
    #[case::arom_omit("#a", AtomPredicate::AromaticValence(AromaticAst::Value(ValueAst::Lit(1))))]
    #[case::multicenter("#m2", AtomPredicate::MulticenterValence(ValueAst::Lit(2)))]
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
