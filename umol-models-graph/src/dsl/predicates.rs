//! Bond and atom predicate parsers.

use nom::branch::alt;
use nom::bytes::complete::{tag, take};
use nom::character::complete::multispace0;
use nom::combinator::{map, success, value};
use nom::sequence::preceded;
use nom::{Err, IResult, Parser};

use super::atom::{AromaticExpr, HydrogenExpr, IsotopeExpr, isotope_expr};
use super::error::ParseError;
use super::value::{value_dsl, ValueAst};

pub(crate) fn optional_value(i: &str) -> IResult<&str, ValueAst, ParseError> {
    preceded(multispace0, alt((value_dsl, success(ValueAst::Lit(1))))).parse(i)
}

pub(crate) fn charge_value(i: &str) -> IResult<&str, ValueAst, ParseError> {
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

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BondPredicate {
    Charge(ValueAst),
    UnpairedElectrons(ValueAst),
    Multiplicity(ValueAst),
}

pub fn bond_predicate(i: &str) -> IResult<&str, BondPredicate, ParseError> {
    let (remaining, prefix) = take(2usize)(i)?;
    match prefix {
        "#c" => map(charge_value, BondPredicate::Charge).parse(remaining),
        "#u" => map(optional_value, BondPredicate::UnpairedElectrons).parse(remaining),
        "#s" => map(optional_value, BondPredicate::Multiplicity).parse(remaining),
        p if p.starts_with("#") => Err(Err::Failure(ParseError::UnknownBondPredicate(p.to_string()))),
        _ => Err(Err::Failure(ParseError::TrailingInput(i.to_string()))),
    }
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

pub fn atom_predicate(i: &str) -> IResult<&str, AtomPredicate, ParseError> {
    let (remaining, prefix) = take(2usize)(i)?;
    match prefix {
        "#i" => preceded(multispace0, map(isotope_expr, AtomPredicate::IsotopeMass)).parse(remaining),
        "#c" => map(charge_value, AtomPredicate::Charge).parse(remaining),
        "#h" => preceded(
            multispace0,
            alt((
                value(AtomPredicate::ImplicitHydrogens(HydrogenExpr::Normal), tag("=")),
                map(value_dsl, |v| AtomPredicate::ImplicitHydrogens(HydrogenExpr::Value(v))),
                success(AtomPredicate::ImplicitHydrogens(HydrogenExpr::Value(ValueAst::Lit(1)))),
            )),
        )
        .parse(remaining),
        "#n" => map(optional_value, AtomPredicate::LonePairs).parse(remaining),
        "#u" => map(optional_value, AtomPredicate::UnpairedElectrons).parse(remaining),
        "#s" => map(optional_value, AtomPredicate::Multiplicity).parse(remaining),
        "#v" => map(optional_value, AtomPredicate::Valence).parse(remaining),
        "#d" => map(optional_value, AtomPredicate::DonatedPairs).parse(remaining),
        "#r" => map(optional_value, AtomPredicate::AcceptedPairs).parse(remaining),
        "#a" => preceded(
            multispace0,
            alt((
                value(AtomPredicate::AromaticValence(AromaticExpr::None), tag("!")),
                map(value_dsl, |v| AtomPredicate::AromaticValence(AromaticExpr::Value(v))),
                success(AtomPredicate::AromaticValence(AromaticExpr::Value(ValueAst::Lit(1)))),
            )),
        )
        .parse(remaining),
        "#m" => map(optional_value, AtomPredicate::MulticenterValence).parse(remaining),
        p if p.starts_with("#") => Err(Err::Failure(ParseError::UnknownAtomPredicate(p.to_string()))),
        _ => Err(Err::Failure(ParseError::TrailingInput(i.to_string()))),
    }
}
