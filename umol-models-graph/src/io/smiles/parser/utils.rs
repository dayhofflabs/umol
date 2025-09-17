//! Utilities for SMILES parser.

use crate::io::ir::Chirality;
use nom::branch::alt;
use nom::bytes::complete::tag;
use nom::character::complete::{digit1, one_of};
use nom::combinator::{map, map_opt, opt, recognize, value};
use nom::error::Error;
use nom::multi::fold_many0;
use nom::sequence::{pair, preceded};
use nom::Parser;
use smallvec::SmallVec;
use umol_data::Element;

const DIGITS: &str = "0123456789";

#[derive(Clone, Debug)]
pub enum BracketField {
    Chiral(Chirality),
    HydrogenCount(u32),
    Charge(i32),
    Class(u32),
}

/// Zero-allocation bracket field parser used by both parser and linter.
/// Returns (element, isotope, tails) where tails are bracket fields in appearance order.
pub fn parse_bracket(inner: &str) -> (Option<Element>, Option<u32>, SmallVec<[BracketField; 4]>) {
    fn isotope<'a>() -> impl Parser<&'a str, Output = u32, Error = Error<&'a str>> {
        map(digit1, |digits: &str| digits.parse::<u32>().unwrap_or(0))
    }

    fn element_symbol<'a>() -> impl Parser<&'a str, Output = Option<Element>, Error = Error<&'a str>> {
        alt((
            value(None, tag("*")),
            map(
                map_opt(
                    recognize(pair(one_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ"), opt(one_of("abcdefghijklmnopqrstuvwxyz")))),
                    |s: &str| Element::from_symbol(s),
                ),
                Some,
            ),
        ))
    }

    fn chiral<'a>() -> impl Parser<&'a str, Output = BracketField, Error = Error<&'a str>> {
        alt((
            value(BracketField::Chiral(Chirality::CounterClockwise), tag("@@")),
            map(preceded(tag("@TH"), digit1), |d: &str| {
                let n = d.as_bytes()[0] - b'0';
                BracketField::Chiral(Chirality::Tetrahedral { arr: n as u32 })
            }),
            map(preceded(tag("@AL"), one_of("12")), |c: char| {
                BracketField::Chiral(Chirality::Allenal { arr: c.to_digit(10).unwrap() as u32 })
            }),
            map(preceded(tag("@SP"), one_of("123")), |c: char| {
                BracketField::Chiral(Chirality::SquarePlanar { arr: c.to_digit(10).unwrap() as u32 })
            }),
            map(preceded(tag("@TB"), digit1), |d: &str| {
                let n = d.as_bytes()[0] - b'0';
                BracketField::Chiral(Chirality::TrigonalBipyramidal { arr: n as u32 })
            }),
            map(preceded(tag("@OH"), digit1), |d: &str| {
                let n = d.as_bytes()[0] - b'0';
                BracketField::Chiral(Chirality::Octahedral { arr: n as u32 })
            }),
            value(BracketField::Chiral(Chirality::Clockwise), tag("@")),
        ))
    }

    fn d1<'a>() -> impl Parser<&'a str, Output = u32, Error = Error<&'a str>> {
        map(one_of(DIGITS), |c| c.to_digit(10).unwrap())
    }

    fn d1_to_2<'a>() -> impl Parser<&'a str, Output = u32, Error = Error<&'a str>> {
        map(pair(one_of(DIGITS), opt(one_of(DIGITS))), |(d1, d2)| {
            let mut v = d1.to_digit(10).unwrap();
            if let Some(c2) = d2 {
                v = v * 10 + c2.to_digit(10).unwrap();
            }
            v
        })
    }

    fn hcount<'a>() -> impl Parser<&'a str, Output = BracketField, Error = Error<&'a str>> {
        map(pair(tag("H"), opt(d1())), |(_, d): (&str, Option<u32>)| BracketField::HydrogenCount(d.unwrap_or(1)))
    }

    fn charge<'a>() -> impl Parser<&'a str, Output = BracketField, Error = Error<&'a str>> {
        alt((
            value(BracketField::Charge(2), tag("++")),
            value(BracketField::Charge(-2), tag("--")),
            map(pair(tag("+"), opt(d1_to_2())), |(_, n): (&str, Option<u32>)| BracketField::Charge(n.unwrap_or(1) as i32)),
            map(pair(tag("-"), opt(d1_to_2())), |(_, n): (&str, Option<u32>)| BracketField::Charge(-(n.unwrap_or(1) as i32))),
        ))
    }

    fn class<'a>() -> impl Parser<&'a str, Output = BracketField, Error = Error<&'a str>> {
        map(preceded(tag(":"), digit1), |digits: &str| BracketField::Class(digits.parse::<u32>().unwrap_or(0)))
    }

    fn field_tail<'a>() -> impl Parser<&'a str, Output = BracketField, Error = Error<&'a str>> {
        alt((chiral(), hcount(), charge(), class()))
    }

    fn fields<'a>() -> impl Parser<&'a str, Output = (Option<Element>, Option<u32>, SmallVec<[BracketField; 4]>), Error = Error<&'a str>> {
        map(
            (
                opt(isotope()),
                element_symbol(),
                fold_many0::<&str, Error<&str>, _, _, _, SmallVec<[BracketField; 4]>>(
                    field_tail(),
                    || SmallVec::<[BracketField; 4]>::new(),
                    |mut acc, f| { acc.push(f); acc }
                ),
            ),
            |(iso_opt, elem_opt, tails): (Option<u32>, Option<Element>, SmallVec<[BracketField; 4]>)| (elem_opt, iso_opt, tails),
        )
    }

    fields().parse(inner).map(|(_, f)| f).unwrap_or_else(|_| (None, None, SmallVec::new()))
}
