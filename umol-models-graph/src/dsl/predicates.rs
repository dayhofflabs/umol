//! Bond predicate parsers — `#c`, `#u`, `#s` dispatch.

use nom::branch::alt;
use nom::bytes::complete::{tag, take};
use nom::character::complete::multispace0;
use nom::combinator::{map, success, value};
use nom::sequence::preceded;
use nom::{Err, IResult, Parser};

use super::error::ParseError;
use super::value::{value_dsl, ValueAst};

/// A single parsed bond predicate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BondPredicate {
    Charge(ValueAst),
    Unpaired(ValueAst),
    Multiplicity(ValueAst),
}

/// Parse one bond predicate (`#c…`, `#u…`, or `#s…`).
///
/// Expects input starting at `#`. Consumes the 2-char prefix via `take(2)`,
/// dispatches on the prefix string, then calls the appropriate body parser
/// on the remainder.
pub fn bond_predicate(i: &str) -> IResult<&str, BondPredicate, ParseError> {
    let (remaining, prefix) = take(2usize)(i)?;
    match prefix {
        "#c" => bond_charge_predicate(remaining),
        "#u" => bond_unpaired_predicate(remaining),
        "#s" => bond_multiplicity_predicate(remaining),
        p if p.starts_with("#") => Err(Err::Failure(ParseError::UnknownBondPredicate(
            p.to_string(),
        ))),
        _ => Err(Err::Failure(ParseError::TrailingInput(i.to_string()))),
    }
}

/// Body parser for `#c` — formal bond charge.
///
/// Tries full `value_dsl` first; falls back to bare `+` (+1) or `-` (-1).
/// Optional whitespace between the tag letter and the payload is allowed.
pub fn bond_charge_predicate(i: &str) -> IResult<&str, BondPredicate, ParseError> {
    preceded(
        multispace0,
        alt((
            map(value_dsl, BondPredicate::Charge),
            value(BondPredicate::Charge(ValueAst::Lit(1)), tag("+")),
            value(BondPredicate::Charge(ValueAst::Lit(-1)), tag("-")),
        )),
    )
    .parse(i)
}

/// Body parser for `#u` — unpaired electrons; omitted payload = 1.
///
/// Optional whitespace between the tag letter and the payload is allowed.
pub fn bond_unpaired_predicate(i: &str) -> IResult<&str, BondPredicate, ParseError> {
    preceded(
        multispace0,
        alt((
            map(value_dsl, BondPredicate::Unpaired),
            success(BondPredicate::Unpaired(ValueAst::Lit(1))),
        )),
    )
    .parse(i)
}

/// Body parser for `#s` — spin multiplicity; omitted payload = 1.
///
/// Optional whitespace between the tag letter and the payload is allowed.
pub fn bond_multiplicity_predicate(i: &str) -> IResult<&str, BondPredicate, ParseError> {
    preceded(
        multispace0,
        alt((
            map(value_dsl, BondPredicate::Multiplicity),
            success(BondPredicate::Multiplicity(ValueAst::Lit(1))),
        )),
    )
    .parse(i)
}
