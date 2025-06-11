//! Properties block parser for CTab files.

use super::utils::{fixed_width_int, fixed_width_int_minus1};
use nom::{
    branch::alt,
    bytes::complete::tag,
    combinator::{all_consuming, map},
    error,
    multi::length_count,
    sequence::{preceded, tuple},
    Parser,
};

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum PropertyLine {
    ChargeLine { entries: Vec<(usize, i8)> },
}

// fn charge_line<'a>(
// ) -> impl Parser<&'a [u8], Output = PropertyLine, Error = error::Error<&'a [u8]>> {
//     let count_parser = fixed_width_int_in_range::<u8>(3, 1..=8);
//     let item_parser = tuple((
//         fixed_width_usize_minus1(3),
//         fixed_width_int_in_range::<i8>(3, -15..=15),
//     ));
//     map(length_count(count_parser, item_parser), |entries| {
//         PropertyLine::ChargeLine { entries }
//     })
// }

// pub(crate) fn property_line<'a>(
// ) -> impl Parser<&'a [u8], Output = PropertyLine, Error = error::Error<&'a [u8]>> {
//     all_consuming(alt((preceded(tag("M  CHG"), charge_line()),)))
// }

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use rstest::rstest;

//     #[rstest]
//     // A single entry
//     #[case(b"M  CHG  1  1 -1", PropertyLine::ChargeLine { entries: vec![(0, -1)] })]
//     // Multiple entries
//     #[case(b"M  CHG  2  1 -1  4  1", PropertyLine::ChargeLine { entries: vec![(0, -1), (3, 1)] })]
//     // Max entries (8)
//     #[case(b"M  CHG  8  1  1  2  2  3  3  4  4  5  5  6  6  7  7  8  8",
//         PropertyLine::ChargeLine { entries: vec![(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8)] }
//     )]
//     // Different spacing
//     #[case(b"M  CHG  1  25 15", PropertyLine::ChargeLine { entries: vec![(24, 15)] })]
//     fn test_charge_line_parser(#[case] input: &[u8], #[case] expected: PropertyLine) {
//         let (remaining, prop) = property_line().parse(input).unwrap();
//         assert_eq!(remaining, b"");
//         assert_eq!(prop, expected);
//     }

//     #[rstest]
//     #[case(b"M  CHG  1  1 -1  ", "trailing space")] // all_consuming should fail
//     #[case(b"M  CHG  2  1 -1", "count does not match item list")] // length_count should fail
//     #[case(b"M  CHG  1  1 -1  4  1", "item list longer than count")] // all_consuming should fail
//     #[case(b"M  CHG  0  1 -1", "count is zero")] // verify on count should fail
//     #[case(b"M  CHG  9  1 -1", "count is > 8")] // verify on count should fail
//     #[case(b"M  XXX  1  1 -1", "invalid property tag")] // alt should fail
//     #[case(b"X  CHG  1  1 -1", "invalid prefix")] // outer tag should fail
//     fn test_charge_line_parser_invalid(#[case] input: &[u8], #[case] message: &str) {
//         let res = property_line().parse(input);
//         assert!(res.is_err(), "Test failed: '{}'. Expected error, got Ok.", message);
//     }
// }
