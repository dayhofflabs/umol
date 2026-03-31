//! Bond and atom predicate parsers

use nom::branch::alt;
use nom::bytes::complete::{tag, take};
use nom::character::complete::{char, multispace0, satisfy, u32 as nom_u32};
use nom::combinator::{map, recognize, success, value};
use nom::multi::{many0, separated_list1};
use nom::sequence::{delimited, pair, preceded, terminated};
use nom::{Err, IResult, Parser};
use serde::{Deserialize, Serialize};
use umol_data::Element;

use super::error::ParseError;
use super::value::{op_char, parse_id, value_dsl, ValueAst};

/// Element expressions
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum ElementExpr {
    Lit(Element),
    Wildcard,
    Set(Vec<Element>),
    Bind { id: String, set: Vec<Element> },
    Ref(String),
}

impl ElementExpr {
    pub fn new(element: Element) -> Self {
        Self::Lit(element)
    }
}

/// Isotope-mass expressions (Natural = #i=)
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum IsotopeExpr {
    Natural,
    Lit(u32),
    Wildcard,
    Set(Vec<u32>),
    Bind { id: String, set: Vec<u32> },
    Ref(String),
}

/// Implicit hydrogen expressions (Normal = #h=)
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum HydrogenExpr {
    Normal,
    Value(ValueAst),
}

impl HydrogenExpr {
    pub fn from_value(value: ValueAst) -> Self {
        Self::Value(value)
    }
}

/// Aromatic valence expressions
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum AromaticExpr {
    Unspecified,
    NotAromatic,
    Value(ValueAst),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AtomPredicate {
    IsotopeMass(IsotopeExpr),
    Charge(ValueAst),
    ImplicitHydrogens(HydrogenExpr),
    LonePairs(ValueAst),
    UnpairedElectrons(ValueAst),
    Multiplicity(ValueAst),
    Valence(ValueAst),
    DonatedPairs(ValueAst),
    AcceptedPairs(ValueAst),
    AromaticValence(AromaticExpr),
    MulticenterValence(ValueAst),
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BondPredicate {
    Charge(ValueAst),
    UnpairedElectrons(ValueAst),
    Multiplicity(ValueAst),
}

pub(crate) fn element_expr(i: &str) -> IResult<&str, ElementExpr, ParseError> {
    alt((
        value(ElementExpr::Wildcard, char('*')),
        map(element_set, ElementExpr::Set),
        map(element_bind, |(id, set)| ElementExpr::Bind { id, set }),
        map(element_ref, ElementExpr::Ref),
        map(element_literal, ElementExpr::Lit),
    ))
    .parse(i)
    .map_err(|_| Err::Error(ParseError::InvalidElement(i.to_string())))
}

fn element_literal(i: &str) -> IResult<&str, Element, ParseError> {
    let (rest, sym) = recognize(pair(
        satisfy(|c: char| c.is_ascii_uppercase()),
        many0(satisfy(|c: char| c.is_ascii_lowercase())),
    ))
    .parse(i)?;
    match Element::from_symbol(sym) {
        Some(el) => Ok((rest, el)),
        None => Err(Err::Error(ParseError::InvalidElement(sym.to_string()))),
    }
}

fn element_set(i: &str) -> IResult<&str, Vec<Element>, ParseError> {
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

fn element_bind(i: &str) -> IResult<&str, (String, Vec<Element>), ParseError> {
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

fn element_ref(i: &str) -> IResult<&str, String, ParseError> {
    delimited(
        char('('),
        delimited(multispace0, preceded(char('?'), parse_id), multispace0),
        char(')'),
    )
    .parse(i)
}

pub(crate) fn atom_predicate(i: &str) -> IResult<&str, AtomPredicate, ParseError> {
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
        p if p.starts_with("#") => Err(Err::Failure(ParseError::UnknownAtomPredicate(
            p.to_string(),
        ))),
        _ => Err(Err::Failure(ParseError::TrailingInput(i.to_string()))),
    }
}

fn isotope_set(i: &str) -> IResult<&str, Vec<u32>, ParseError> {
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

fn isotope_bind(i: &str) -> IResult<&str, (String, Vec<u32>), ParseError> {
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

fn isotope_ref(i: &str) -> IResult<&str, String, ParseError> {
    delimited(
        char('('),
        delimited(multispace0, preceded(char('?'), parse_id), multispace0),
        char(')'),
    )
    .parse(i)
}

fn isotope_expr(i: &str) -> IResult<&str, IsotopeExpr, ParseError> {
    preceded(
        multispace0,
        alt((
            value(IsotopeExpr::Natural, char('=')),
            value(IsotopeExpr::Wildcard, char('*')),
            map(isotope_set, IsotopeExpr::Set),
            map(isotope_bind, |(id, set)| IsotopeExpr::Bind { id, set }),
            map(isotope_ref, IsotopeExpr::Ref),
            map(nom_u32, IsotopeExpr::Lit),
        )),
    )
    .parse(i)
}

fn charge_value(i: &str) -> IResult<&str, ValueAst, ParseError> {
    preceded(
        multispace0,
        alt((
            value_dsl,
            value(ValueAst::Lit(1), tag("+")),
            value(ValueAst::Lit(-1), tag("-")),
        )),
    )
    .parse(i)
}

fn hydrogen_expr(i: &str) -> IResult<&str, HydrogenExpr, ParseError> {
    preceded(
        multispace0,
        alt((
            value(HydrogenExpr::Normal, tag("=")),
            map(value_dsl, |v| HydrogenExpr::Value(v)),
            success(HydrogenExpr::Value(ValueAst::Lit(1))),
        )),
    )
    .parse(i)
}

fn aromatic_valence_expr(i: &str) -> IResult<&str, AromaticExpr, ParseError> {
    preceded(
        multispace0,
        alt((
            value(AromaticExpr::NotAromatic, tag("!")),
            value(AromaticExpr::Unspecified, tag("?")),
            map(value_dsl, |v| AromaticExpr::Value(v)),
            success(AromaticExpr::Value(ValueAst::Lit(1))),
        )),
    )
    .parse(i)
}

fn optional_value(i: &str) -> IResult<&str, ValueAst, ParseError> {
    preceded(multispace0, alt((value_dsl, success(ValueAst::Lit(1))))).parse(i)
}

pub(crate) fn bond_order(i: &str) -> IResult<&str, ValueAst, ParseError> {
    value_dsl(i).map_err(|_| Err::Failure(ParseError::InvalidBondOrder(i.to_string())))
}

pub(crate) fn bond_predicate(i: &str) -> IResult<&str, BondPredicate, ParseError> {
    let (remaining, prefix) = take(2usize)(i)?;
    match prefix {
        "#c" => map(charge_value, BondPredicate::Charge).parse(remaining),
        "#u" => map(optional_value, BondPredicate::UnpairedElectrons).parse(remaining),
        "#s" => map(optional_value, BondPredicate::Multiplicity).parse(remaining),
        p if p.starts_with("#") => Err(Err::Failure(ParseError::UnknownBondPredicate(
            p.to_string(),
        ))),
        _ => Err(Err::Failure(ParseError::TrailingInput(i.to_string()))),
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_data::Element;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::carbon("C", ElementExpr::Lit(Element::C))]
    #[case::iron("Fe", ElementExpr::Lit(Element::Fe))]
    #[case::chlorine("Cl", ElementExpr::Lit(Element::Cl))]
    #[case::wildcard("*", ElementExpr::Wildcard)]
    #[case::set("{C,N,O}", ElementExpr::Set(vec![Element::C, Element::N, Element::O]))]
    #[case::set_spaced("{ C, N}", ElementExpr::Set(vec![Element::C, Element::N]))]
    #[case::bind("(?e :: {C,N})", ElementExpr::Bind { id: "e".to_string(), set: vec![Element::C, Element::N] })]
    #[case::ref_("(?e)", ElementExpr::Ref("e".to_string()))]
    fn test_element_expr(#[case] input: &str, #[case] expected: ElementExpr) {
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
    #[case::natural("=", IsotopeExpr::Natural)]
    #[case::lit("12", IsotopeExpr::Lit(12))]
    #[case::wildcard("*", IsotopeExpr::Wildcard)]
    #[case::set("{12,13,14}", IsotopeExpr::Set(vec![12, 13, 14]))]
    #[case::bind("(?m :: {12,13})", IsotopeExpr::Bind { id: "m".to_string(), set: vec![12, 13] })]
    #[case::ref_("(?m)", IsotopeExpr::Ref("m".to_string()))]
    fn test_isotope_expr(#[case] input: &str, #[case] expected: IsotopeExpr) {
        let result = isotope_expr(input);
        assert!(result.is_ok(), "{input:?} should succeed, got {:?}", result.unwrap_err());
        let (remaining, expr) = result.unwrap();
        assert!(remaining.is_empty(), "{input:?} should consume all input, remaining: {remaining:?}");
        assert_eq!(expr, expected);
    }
}
