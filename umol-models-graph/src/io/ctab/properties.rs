// //! Properties block parser for CTab files.

// use super::utils::{fixed_width_int_in_range, fixed_width_int_minus1};
// use nom::{
//     branch::alt,
//     bytes::complete::tag,
//     combinator::{all_consuming, map},
//     error,
//     multi::length_count,
//     sequence::{preceded, tuple},
//     Parser,
// };

// #[derive(Debug, Clone, PartialEq)]
// pub(crate) enum PropertyLine {
//     ChargeLine { entries: Vec<(usize, i8)> },
// }

// /// Parses the data part of a `CHG` line, which consists of a count
// /// followed by that many (atom_index, charge) pairs.
// fn charge_data<'a>(
// ) -> impl Parser<&'a [u8], Output = Vec<(usize, i8)>, Error = error::Error<&'a [u8]>> {
//     length_count(
//         // `nn`: Number of entries (1-8)
//         fixed_width_int_in_range::<u8, _>(3, 1..=8),
//         // A pair of (atom_index, charge_value)
//         tuple((
//             // `aaa`: Atom index (1-based)
//             fixed_width_int_minus1::<usize>(3),
//             // `vvv`: Charge value (-15 to 15)
//             fixed_width_int_in_range::<i8, _>(3, -15..=15),
//         )),
//     )
// }

// /// Parses a complete `M ...` property line.
// pub(crate) fn property_line<'a>(
// ) -> impl Parser<&'a [u8], Output = PropertyLine, Error = error::Error<&'a [u8]>> {
//     all_consuming(alt((
//         // M  CHG line
//         map(
//             preceded(tag("M  CHG"), charge_data()),
//             |entries| PropertyLine::ChargeLine { entries },
//         ),
//         // TODO: Add other properties like RAD, ISO etc. here
//     )))
// }

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use nom::Parser;
//     use rstest::rstest;

//     #[rstest]
//     // A single entry
//     #[case(b"M  CHG  1   1  -1", PropertyLine::ChargeLine { entries: vec![(0, -1)] })]
//     // Multiple entries
//     #[case(b"M  CHG  2   1  -1   4   1", PropertyLine::ChargeLine { entries: vec![(0, -1), (3, 1)] })]
//     // Max entries (8)
//     #[case(b"M  CHG  8   1   1   2   2   3   3   4   4   5   5   6   6   7   7   8   8",
//         PropertyLine::ChargeLine { entries: vec![(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8)] }
//     )]
//     // Different spacing for values
//     #[case(b"M  CHG  1  25  15", PropertyLine::ChargeLine { entries: vec![(24, 15)] })]
//     fn test_charge_line_parser(#[case] input: &[u8], #[case] expected: PropertyLine) {
//         let (remaining, prop) = property_line().parse(input).unwrap();
//         assert_eq!(remaining, b"");
//         assert_eq!(prop, expected);
//     }

//     #[rstest]
//     #[case(b"M  CHG  1   1  -1  ", "trailing space")]
//     #[case(b"M  CHG  2   1  -1", "count does not match item list")]
//     #[case(b"M  CHG  1   1  -1   4   1", "item list longer than count")]
//     #[case(b"M  CHG  0   1  -1", "count is zero")]
//     #[case(b"M  CHG  9   1  -1", "count is > 8")]
//     #[case(b"M  CHG  1   1 -16", "charge out of range")]
//     #[case(b"M  CHG  1   1  16", "charge out of range")]
//     #[case(b"M  XXX  1   1  -1", "invalid property tag")]
//     #[case(b"X  CHG  1   1  -1", "invalid prefix")]
//     fn test_charge_line_parser_invalid(#[case] input: &[u8], #[case] message: &str) {
//         let res = property_line().parse(input);
//         assert!(res.is_err(), "Test failed: '{}'. Expected error, got Ok.", message);
//     }
// }
