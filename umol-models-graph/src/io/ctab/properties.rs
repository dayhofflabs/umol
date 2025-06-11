//! Properties block parser for CTab files.

use nom::{
    branch::alt,
    bytes::complete::tag,
    character::complete::space0,
    combinator::{all_consuming, map, peek, verify},
    error,
    multi::length_count,
    sequence::{delimited, preceded},
    Parser,
};

use super::utils::{fixed_width_int, fixed_width_int_in_range, fixed_width_int_minus1};

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ChargeLine {
    entries: Vec<(usize, i8)>,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RadicalLine {
    entries: Vec<(usize, i8)>,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct IsotopeLine {
    entries: Vec<(usize, u32)>,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum PropertyLine {
    ChargeLine(ChargeLine),
    RadicalLine(RadicalLine),
    IsotopeLine(IsotopeLine),
}

/// Parse charge property line (*[Generic]*).
/// M  CHGnn8 aaa vvv ...
/// vvv: -15..= 15.
fn charge_line<'a>() -> impl Parser<&'a [u8], Output = PropertyLine, Error = error::Error<&'a [u8]>>
{
    all_consuming(map(
        delimited(
            tag("M  CHG"),
            length_count(
                fixed_width_int_in_range::<u8, _>(3, 1..=8),
                (
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                    preceded(tag(" "), fixed_width_int_in_range::<i8, _>(3, -15..=15)),
                ),
            ),
            space0,
        ),
        |entries| PropertyLine::ChargeLine(ChargeLine { entries }),
    ))
}

/// Parse radical property line (*[Generic]*).
/// M  RADnn8 aaa vvv ...
/// vvv: 0..= 3: 0 = no radical, 1 = singlet (:), 2 = doublet (. or ^), 3 = triplet (^^).
fn radical_line<'a>() -> impl Parser<&'a [u8], Output = PropertyLine, Error = error::Error<&'a [u8]>>
{
    all_consuming(map(
        delimited(
            tag("M  RAD"),
            length_count(
                fixed_width_int_in_range::<u8, _>(3, 1..=8),
                (
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                    preceded(tag(" "), fixed_width_int_in_range::<i8, _>(3, 0..=3)),
                ),
            ),
            space0,
        ),
        |entries| PropertyLine::RadicalLine(RadicalLine { entries }),
    ))
}

/// Parse isotope property line (*[Generic]*).
/// M  ISOnn8 aaa vvv ...
/// vvv: isotope mass number (not difference)
/// Difference between the isotope mass number and reference isotope mass number
/// should be in the range -18..=12.
fn isotope_line<'a>() -> impl Parser<&'a [u8], Output = PropertyLine, Error = error::Error<&'a [u8]>>
{
    all_consuming(map(
        delimited(
            tag("M  ISO"),
            length_count(
                fixed_width_int_in_range::<u8, _>(3, 1..=8),
                (
                    preceded(tag(" "), fixed_width_int_minus1::<usize>(3)),
                    preceded(tag(" "), fixed_width_int::<u32>(3)),
                ),
            ),
            space0,
        ),
        |entries| PropertyLine::IsotopeLine(IsotopeLine { entries }),
    ))
}

/// Parse property line
pub(crate) fn property_line<'a>(
) -> impl Parser<&'a [u8], Output = PropertyLine, Error = error::Error<&'a [u8]>> {
    alt((
        preceded(peek(tag("M  CHG")), charge_line()),
        preceded(peek(tag("M  RAD")), radical_line()),
        preceded(peek(tag("M  ISO")), isotope_line()),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nom::{error::ErrorKind, Err, Parser};
    use rstest::rstest;

    #[rstest]
    #[case(b"M  CHG  1   1  -1", PropertyLine::ChargeLine(ChargeLine { entries: vec![(0, -1)] }))]
    #[case(b"M  CHG  2   1  -1   4   1", PropertyLine::ChargeLine(ChargeLine { entries: vec![(0, -1), (3, 1)] }))]
    #[case(b"M  CHG  8   1   1   2   2   3   3   4   4   5   5   6   6   7   7   8   8",
        PropertyLine::ChargeLine(ChargeLine { entries: vec![(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8)] }
    ))]
    #[case(b"M  CHG  1  25  15", PropertyLine::ChargeLine(ChargeLine { entries: vec![(24, 15)] }))]
    fn test_charge_line(#[case] input: &[u8], #[case] expected: PropertyLine) {
        let (remaining, result) = charge_line().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"M  CHG  1   1  -1  a", "trailing chars", ErrorKind::Eof)]
    #[case(b"M  CHG  2   1  -1", "count does not match item list", ErrorKind::Tag)]
    #[case(
        b"M  CHG  1   1  -1   4   1",
        "item list longer than count",
        ErrorKind::Eof
    )]
    #[case(b"M  CHG  0", "count is zero", ErrorKind::Verify)]
    #[case(b"M  CHG  1   0 -10", "atom index is zero", ErrorKind::Verify)]
    #[case(b"M  XXX  1   1  -1", "invalid property tag", ErrorKind::Tag)]
    #[case(b"X  CHG  1   1  -1", "invalid prefix", ErrorKind::Tag)]
    fn test_charge_line_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let result = charge_line().parse(input);
        assert!(result.is_err(), "{} should have failed", desc);
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {}",
            desc,
            expected_kind,
            result.clone().unwrap_err().map(|e| e.code),
        );
    }

    #[rstest]
    #[case(b"M  RAD  1   1   2", PropertyLine::RadicalLine(RadicalLine { entries: vec![(0, 2)] }))]
    #[case(b"M  RAD  2   1   1   4   3", PropertyLine::RadicalLine(RadicalLine { entries: vec![(0, 1), (3, 3)] }))]
    fn test_radical_line(#[case] input: &[u8], #[case] expected: PropertyLine) {
        let (remaining, result) = radical_line().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"M  RAD  1   1   4", "value out of range", ErrorKind::Verify)]
    #[case(b"M  RAD  1   1  -1", "value out of range", ErrorKind::Verify)]
    fn test_radical_line_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let result = radical_line().parse(input);
        assert!(result.is_err(), "{} should have failed", desc);
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {}",
            desc,
            expected_kind,
            result.clone().unwrap_err().map(|e| e.code),
        );
    }

    #[rstest]
    #[case(b"M  ISO  1   1  13", PropertyLine::IsotopeLine(IsotopeLine { entries: vec![(0, 13)] }))]
    #[case(b"M  ISO  2  12   2  15  14", PropertyLine::IsotopeLine(IsotopeLine { entries: vec![(11, 2), (14, 14)] }))]
    fn test_isotope_line(#[case] input: &[u8], #[case] expected: PropertyLine) {
        let (remaining, result) = isotope_line().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"M  CHG  1   1  -1", PropertyLine::ChargeLine(ChargeLine { entries: vec![(0, -1)] }))]
    #[case(b"M  RAD  1   1   2", PropertyLine::RadicalLine(RadicalLine { entries: vec![(0, 2)] }))]
    #[case(b"M  ISO  1   1  13", PropertyLine::IsotopeLine(IsotopeLine { entries: vec![(0, 13)] }))]
    fn test_property_line_dispatch(#[case] input: &[u8], #[case] expected: PropertyLine) {
        let (remaining, result) = property_line().parse(input).unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }
}
