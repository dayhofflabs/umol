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
    Charge(ValueAst<i8>),
    Unpaired(ValueAst<u8>),
    Multiplicity(ValueAst<u8>),
}

/// Parse one bond predicate (`#c…`, `#u…`, or `#s…`).
///
/// Expects input starting at `#`. Consumes the 2-char prefix via `take(2)`,
/// dispatches on the prefix string, then calls the appropriate body parser
/// on the remainder.
pub fn bond_predicate(i: &str) -> IResult<&str, BondPredicate, ParseError> {
    let (remaining, prefix) = take(2usize)(i)?;
    match prefix {
        "#c" => charge_predicate(remaining),
        "#u" => unpaired_predicate(remaining),
        "#s" => multiplicity_predicate(remaining),
        _ => Err(Err::Error(ParseError::UnknownBondPredicate)),
    }
}

/// Body parser for `#c` — formal bond charge (`i8`).
///
/// Tries full `value_dsl::<i8>` first; falls back to bare `+` (+1) or `-` (-1).
/// Optional whitespace between the tag letter and the payload is allowed.
pub fn charge_predicate(i: &str) -> IResult<&str, BondPredicate, ParseError> {
    preceded(
        multispace0,
        alt((
            map(value_dsl::<i8>, BondPredicate::Charge),
            value(BondPredicate::Charge(ValueAst::Lit(1)), tag("+")),
            value(BondPredicate::Charge(ValueAst::Lit(-1)), tag("-")),
        )),
    )
    .parse(i)
}

/// Body parser for `#u` — unpaired electrons (`u8`); omitted payload = 1.
///
/// Optional whitespace between the tag letter and the payload is allowed.
pub fn unpaired_predicate(i: &str) -> IResult<&str, BondPredicate, ParseError> {
    preceded(
        multispace0,
        alt((
            map(value_dsl::<u8>, BondPredicate::Unpaired),
            success(BondPredicate::Unpaired(ValueAst::Lit(1))),
        )),
    )
    .parse(i)
}

/// Body parser for `#s` — spin multiplicity (`u8`); omitted payload = 1.
///
/// Optional whitespace between the tag letter and the payload is allowed.
pub fn multiplicity_predicate(i: &str) -> IResult<&str, BondPredicate, ParseError> {
    preceded(
        multispace0,
        alt((
            map(value_dsl::<u8>, BondPredicate::Multiplicity),
            success(BondPredicate::Multiplicity(ValueAst::Lit(1))),
        )),
    )
    .parse(i)
}
