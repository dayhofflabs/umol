//! Utilities for SMILES parser.

use atoi;
use nom::branch::alt;
use nom::bytes::complete::{tag, take_while_m_n};
use nom::character::complete::{digit1, one_of};
use nom::combinator::{map, map_res, opt, value};
use nom::error::Error;
use nom::multi::fold_many0;
use nom::sequence::{pair, preceded};
use nom::Parser;
use smallvec::SmallVec;
use umol_data::Element;

use crate::io::ir::Chirality;

#[derive(Clone, Debug, PartialEq)]
pub enum BracketField {
    Chiral(Chirality),
    HydrogenCount(u32),
    Charge(i32),
    Class(u32),
}

#[derive(Default, Debug, Clone, Copy)]
pub struct BracketFields {
    pub element: Option<Element>,
    pub isotope: Option<u32>,
    pub hcount: Option<u32>,
    pub charge: Option<i32>,
    pub class: Option<u32>,
}

/// Zero-allocation bracket field parser used by both parser and linter.
/// Returns (element, isotope, fields) where fields are bracket fields in appearance order.
pub fn parse_bracket(inner: &str) -> (Option<Element>, Option<u32>, SmallVec<[BracketField; 4]>) {
    fields()
        .parse(inner)
        .map(|(_, f)| f)
        .unwrap_or_else(|_| (None, None, SmallVec::new()))
}

/// Tail bracket fields can be in any order, return in order of appearance
fn fields<'a>() -> impl Parser<
    &'a str,
    Output = (Option<Element>, Option<u32>, SmallVec<[BracketField; 4]>),
    Error = Error<&'a str>,
> {
    map(
        (
            opt(isotope()),
            element_symbol(),
            fold_many0::<&str, Error<&str>, _, _, _, SmallVec<[BracketField; 4]>>(
                field(),
                || SmallVec::<[BracketField; 4]>::new(),
                |mut acc, f| {
                    acc.push(f);
                    acc
                },
            ),
        ),
        |(iso_opt, elem_opt, tails): (
            Option<u32>,
            Option<Element>,
            SmallVec<[BracketField; 4]>,
        )| (elem_opt, iso_opt, tails),
    )
}

/// Parse tail bracket field: chiral, hydrogen count, charge, class
fn field<'a>() -> impl Parser<&'a str, Output = BracketField, Error = Error<&'a str>> {
    alt((chiral(), hcount(), charge(), class()))
}

/// Parse isotope field: d, must come first in bracket
fn isotope<'a>() -> impl Parser<&'a str, Output = u32, Error = Error<&'a str>> {
    map(digit1, |digits: &str| digits.parse::<u32>().unwrap_or(0))
}

/// Parse element symbol, follows isotope field
fn element_symbol<'a>() -> impl Parser<&'a str, Output = Option<Element>, Error = Error<&'a str>> {
    alt((
        // Try two-letter element first (e.g., Cl, Br). If not a valid element, backtrack.
        map_res(
            take_while_m_n(2, 2, |c: char| c.is_ascii_alphabetic()),
            |s: &str| Element::from_symbol(s).ok_or("Invalid element symbol").map(Some),
        ),
        // Then try one-letter element (e.g., C, N, O)
        map_res(
            take_while_m_n(1, 1, |c: char| c.is_ascii_alphabetic()),
            |s: &str| Element::from_symbol(s).ok_or("Invalid element symbol").map(Some),
        ),
        // Or wildcard '*'
        value(None, tag("*")),
    ))
}

/// Parse chiral field: @, @@, @THn, n = 1..2, @ALn, n = 1..2, @SPn, n = 1..3,
/// @TBn, n = 1..20, @OHn, n = 1..30
fn chiral<'a>() -> impl Parser<&'a str, Output = BracketField, Error = Error<&'a str>> {
    alt((
        value(BracketField::Chiral(Chirality::CounterClockwise), tag("@@")),
        map(preceded(tag("@TH"), digit1), |d: &str| {
            let n = d.as_bytes()[0] - b'0';
            BracketField::Chiral(Chirality::Tetrahedral { arr: n as u32 })
        }),
        map(preceded(tag("@AL"), one_of("12")), |c| {
            BracketField::Chiral(Chirality::Allenal {
                arr: c.to_digit(10).unwrap() as u32,
            })
        }),
        map(preceded(tag("@SP"), one_of("123")), |c| {
            BracketField::Chiral(Chirality::SquarePlanar {
                arr: c.to_digit(10).unwrap() as u32,
            })
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

/// Parse hydrogen count field: H = H1 or Hd, d = 0..9
fn hcount<'a>() -> impl Parser<&'a str, Output = BracketField, Error = Error<&'a str>> {
    map(
        preceded(
            tag("H"),
            map(
                take_while_m_n(0, 1, |c: char| c.is_ascii_digit()),
                |s: &str| atoi::atoi::<u32>(s.as_bytes()).unwrap_or(1),
            ),
        ),
        BracketField::HydrogenCount,
    )
}

/// Parse charge field: + = +1, - = -1, ++ = 2, -- = -2, +d = d, -d = -d, d = 0..99
fn charge<'a>() -> impl Parser<&'a str, Output = BracketField, Error = Error<&'a str>> {
    alt((
        value(BracketField::Charge(2), tag("++")),
        value(BracketField::Charge(-2), tag("--")),
        map(
            pair(
                alt((value(1, tag("+")), value(-1, tag("-")))),
                map(
                    take_while_m_n(0, 2, |c: char| c.is_ascii_digit()),
                    |s: &str| atoi::atoi::<i32>(s.as_bytes()).unwrap_or(1),
                ),
            ),
            |(s, n)| BracketField::Charge(n * s),
        ),
    ))
}

/// Parse class field: :d, d = 0..999
fn class<'a>() -> impl Parser<&'a str, Output = BracketField, Error = Error<&'a str>> {
    map(preceded(tag(":"), digit1), |digits: &str| {
        BracketField::Class(atoi::atoi::<u32>(digits.as_bytes()).unwrap_or(0))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;

    #[rstest]
    #[case("C", Some(Element::C), None, vec![])]
    #[case("*", None, None, vec![])]
    #[case("13C", Some(Element::C), Some(13), vec![])]
    #[case("13*", None, Some(13), vec![])]
    #[case("*H", None, None, vec![BracketField::HydrogenCount(1)])]
    #[case("*H0", None, None, vec![BracketField::HydrogenCount(0)])]
    #[case("C+", Some(Element::C), None, vec![BracketField::Charge(1)])]
    #[case("C-0", Some(Element::C), None, vec![BracketField::Charge(0)])]
    #[case("C++", Some(Element::C), None, vec![BracketField::Charge(2)])]
    #[case("C--", Some(Element::C), None, vec![BracketField::Charge(-2)])]
    #[case("C:12", Some(Element::C), None, vec![BracketField::Class(12)])]
    #[case("C@", Some(Element::C), None, vec![BracketField::Chiral(Chirality::Clockwise)])]
    #[case("C@@", Some(Element::C), None, vec![BracketField::Chiral(Chirality::CounterClockwise)])]
    #[case("C@TH2", Some(Element::C), None, vec![BracketField::Chiral(Chirality::Tetrahedral { arr: 2 })])]
    #[case("C@AL1", Some(Element::C), None, vec![BracketField::Chiral(Chirality::Allenal { arr: 1 })])]
    #[case("C@SP3", Some(Element::C), None, vec![BracketField::Chiral(Chirality::SquarePlanar { arr: 3 })])]
    #[case("C@TB5", Some(Element::C), None, vec![BracketField::Chiral(Chirality::TrigonalBipyramidal { arr: 5 })])]
    #[case("C@OH7", Some(Element::C), None, vec![BracketField::Chiral(Chirality::Octahedral { arr: 7 })])]
    #[case("X", None, None, vec![])]
    fn test_parse_bracket(
        #[case] input: &str,
        #[case] elem: Option<Element>,
        #[case] iso: Option<u32>,
        #[case] tails: Vec<BracketField>,
    ) {
        let (e, i, t) = parse_bracket(input);
        assert_eq!(e, elem);
        assert_eq!(i, iso);
        assert_eq!(t.as_slice(), tails.as_slice());
    }

    #[rstest]
    #[case(
        "C@H2+3:12",
        vec![
            BracketField::Chiral(Chirality::Clockwise),
            BracketField::HydrogenCount(2),
            BracketField::Charge(3),
            BracketField::Class(12),
        ]
    )]
    #[case(
        "C:12H2-2@@",
        vec![
            BracketField::Class(12),
            BracketField::HydrogenCount(2),
            BracketField::Charge(-2),
            BracketField::Chiral(Chirality::CounterClockwise),
        ]
    )]
    #[case(
        "*+2H",
        vec![
            BracketField::Charge(2),
            BracketField::HydrogenCount(1),
        ]
    )]
    fn test_parse_bracket_tail_fields(#[case] input: &str, #[case] tails: Vec<BracketField>) {
        let (_e, _i, t) = parse_bracket(input);
        assert_eq!(t.as_slice(), tails.as_slice());
    }
}
