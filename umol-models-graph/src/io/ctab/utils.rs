//! Parsing utilities for CTab files.

use std::{
    fmt::Debug,
    ops::{Range, RangeInclusive},
};

use fast_float::FastFloat;
use nom::{
    branch::alt,
    bytes::{complete::tag, take, take_while_m_n},
    character::complete::{
        digit0, digit1, i32 as nom_i32, i8 as nom_i8, space0, u8 as nom_u8, usize as nom_usize,
    },
    combinator::{all_consuming, complete, map, map_parser, opt, recognize, value, verify},
    error,
    number::complete::{double, float},
    sequence::delimited,
    Parser,
};
use num::{Float, Integer};

pub(crate) trait Contains<T: PartialOrd> {
    fn contains(&self, value: &T) -> bool;
}

impl<T: PartialOrd> Contains<T> for Range<T> {
    fn contains(&self, value: &T) -> bool {
        Range::contains(self, value)
    }
}

impl<T: PartialOrd> Contains<T> for RangeInclusive<T> {
    fn contains(&self, value: &T) -> bool {
        RangeInclusive::contains(self, value)
    }
}

pub(crate) trait IntParser: Sized + Copy + PartialOrd + Debug + Default + Integer {
    fn nom_parser<'a>() -> impl Parser<&'a [u8], Output = Self, Error = error::Error<&'a [u8]>>;
}

impl IntParser for i8 {
    fn nom_parser<'a>() -> impl Parser<&'a [u8], Output = Self, Error = error::Error<&'a [u8]>> {
        nom_i8
    }
}

impl IntParser for i32 {
    fn nom_parser<'a>() -> impl Parser<&'a [u8], Output = Self, Error = error::Error<&'a [u8]>> {
        nom_i32
    }
}

impl IntParser for u8 {
    fn nom_parser<'a>() -> impl Parser<&'a [u8], Output = Self, Error = error::Error<&'a [u8]>> {
        nom_u8
    }
}

impl IntParser for usize {
    fn nom_parser<'a>() -> impl Parser<&'a [u8], Output = Self, Error = error::Error<&'a [u8]>> {
        nom_usize
    }
}

pub(crate) trait FloatParser: Sized + Copy + Debug + Default + Float {
    fn nom_parser<'a>() -> impl Parser<&'a [u8], Output = Self, Error = error::Error<&'a [u8]>>;
}

impl FloatParser for f32 {
    fn nom_parser<'a>() -> impl Parser<&'a [u8], Output = Self, Error = error::Error<&'a [u8]>> {
        float
    }
}

impl FloatParser for f64 {
    fn nom_parser<'a>() -> impl Parser<&'a [u8], Output = Self, Error = error::Error<&'a [u8]>> {
        double
    }
}

/// Parse a fixed-width field as an integer type. Interprets empty/whitespace field as default.
pub(crate) fn fixed_width_int<'a, T>(
    width: usize,
) -> impl Parser<&'a [u8], Output = T, Error = error::Error<&'a [u8]>>
where
    T: IntParser,
{
    complete(alt((
        map_parser(
            take(width),
            all_consuming(delimited(space0, T::nom_parser(), space0)),
        ),
        value(T::zero(), take_while_m_n(width, width, |c| c == b' ')),
    )))
}

/// Parse a fixed-width field as an integer type, applying range bounds.
pub(crate) fn fixed_width_int_in_range<'a, T, R>(
    width: usize,
    range: R,
) -> impl Parser<&'a [u8], Output = T, Error = error::Error<&'a [u8]>>
where
    T: IntParser,
    R: Contains<T> + Clone,
{
    complete(map_parser(
        take(width),
        verify(
            all_consuming(delimited(space0, T::nom_parser(), space0)),
            move |val: &T| range.contains(val),
        ),
    ))
}

/// Parse a fixed-width field as an integer type, subtracting one.
pub(crate) fn fixed_width_int_minus1<'a, T>(
    width: usize,
) -> impl Parser<&'a [u8], Output = T, Error = error::Error<&'a [u8]>>
where
    T: IntParser,
{
    map(fixed_width_int(width), |x: T| x - T::one())
}

/// Parse a fixed-width field as an integer type with a range check, subtracting one.
pub(crate) fn fixed_width_int_in_range_minus1<'a, T, R>(
    width: usize,
    range: R,
) -> impl Parser<&'a [u8], Output = T, Error = error::Error<&'a [u8]>>
where
    T: IntParser,
    R: Contains<T> + Clone,
{
    map(fixed_width_int_in_range(width, range), |x: T| x - T::one())
}

/// Parse a fixed-width field as a float type with Fortran semantics (Fw.d).
pub(crate) fn fixed_width_float<'a, T>(
    width: usize,
    precision: usize,
) -> impl Parser<&'a [u8], Output = T, Error = error::Error<&'a [u8]>>
where
    T: FloatParser + FastFloat,
{
    map_parser(
        take(width),
        alt((
            map(
                all_consuming(delimited(
                    space0,
                    recognize((opt(tag(&b"-"[..])), digit1, tag(&b"."[..]), digit0)),
                    space0,
                )),
                |s| fast_float::parse::<T, _>(s).unwrap(),
            ),
            map(
                all_consuming(delimited(
                    space0,
                    recognize((opt(tag(&b"-"[..])), digit1)),
                    space0,
                )),
                move |s| {
                    fast_float::parse::<T, _>(s).unwrap()
                        / T::from(10.0).unwrap().powi(precision as i32)
                },
            ),
            value(T::zero(), all_consuming(space0)),
        )),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use nom::error::ErrorKind;
    use nom::Err;
    use rstest::rstest;

    #[rstest]
    #[case(b"123", 123i32)]
    #[case(b"-98", -98i32)]
    #[case(b"  8", 8i32)]
    #[case(b"   ", 0i32)]
    fn test_fixed_width_int(#[case] input: &[u8], #[case] expected: i32) {
        let mut parser = all_consuming(fixed_width_int::<i32>(3));
        let result = parser.parse(input);
        assert!(result.is_ok(), "Test for '{}' should have succeeded", String::from_utf8_lossy(input));
        let (remaining, result) = result.unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"1234", "too many characters", ErrorKind::Eof)]
    #[case(b"12", "too few characters", ErrorKind::TakeWhileMN)]
    #[case(b"abc", "non-numeric input", ErrorKind::TakeWhileMN)]
    #[case(b"1a ", "trailing characters", ErrorKind::TakeWhileMN)]
    fn test_fixed_width_int_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let mut parser = all_consuming(fixed_width_int::<i32>(3));
        let result = parser.parse(input);
        assert!(result.is_err(), "Test for should have failed for {}", desc);
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {}",
            desc,
            expected_kind,
            result.clone().unwrap_err().map(|e| e.code),
        );
    }

    #[rstest]
    #[case(b"100", 100i8)]
    #[case(b" -9", -9i8)]
    #[case(b"8  ", 8i8)]
    #[case(b" 1 ", 1i8)]
    fn test_fixed_width_int_in_range(#[case] input: &[u8], #[case] expected: i8) {
        let mut parser = all_consuming(fixed_width_int_in_range::<i8, _>(3, -10i8..=110i8));
        let result = parser.parse(input);
        assert!(result.is_ok(), "Test for '{}' should have succeeded", String::from_utf8_lossy(input));
        let (remaining, result) = result.unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"1234", "too many characters", ErrorKind::Verify)]
    #[case(b"8", "too few characters", ErrorKind::Eof)]
    #[case(b"abc", "non-numeric input", ErrorKind::Digit)]
    #[case(b"1a ", "trailing characters", ErrorKind::Eof)]
    fn test_fixed_width_int_in_range_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let mut parser = all_consuming(fixed_width_int_in_range::<i8, _>(3, -10i8..=10i8));
        let result = parser.parse(input);
        assert!(result.is_err(), "Test for should have failed for {}", desc);
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {}",
            desc,
            expected_kind,
            result.clone().unwrap_err().map(|e| e.code),
        );
    }

    #[rstest]
    #[case(b"100", 100u8)]
    #[case(b"  9", 9u8)]
    #[case(b"8  ", 8u8)]
    #[case(b" 1 ", 1u8)]
    fn test_fixed_width_int_in_range_inclusive(#[case] input: &[u8], #[case] expected: u8) {
        let mut parser = all_consuming(fixed_width_int_in_range::<u8, _>(3, 0u8..=100u8));
        let result = parser.parse(input);
        assert!(result.is_ok(), "Test for '{}' should have succeeded", String::from_utf8_lossy(input));
        let (remaining, result) = result.unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"  1", 0usize)]
    #[case(b"123", 122usize)]
    fn test_fixed_width_int_minus1(#[case] input: &[u8], #[case] expected: usize) {
        let mut parser = all_consuming(fixed_width_int_minus1::<usize>(3));
        let result = parser.parse(input);
        assert!(result.is_ok(), "Test for '{}' should have succeeded", String::from_utf8_lossy(input));
        let (remaining, result) = result.unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"  2", 1usize)]
    #[case(b"100", 99usize)]
    fn test_fixed_width_int_in_range_minus1(#[case] input: &[u8], #[case] expected: usize) {
        let mut parser = all_consuming(fixed_width_int_in_range_minus1::<usize, _>(3, 1..=100));
        let result = parser.parse(input);
        assert!(result.is_ok(), "Test for '{}' should have succeeded", String::from_utf8_lossy(input));
        let (remaining, result) = result.unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"101", "out of range", ErrorKind::Verify)]
    #[case(b"  0", "out of range", ErrorKind::Verify)]
    fn test_fixed_width_int_in_range_minus1_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let mut parser = all_consuming(fixed_width_int_in_range_minus1::<usize, _>(3, 1..=100));
        let result = parser.parse(input);
        assert!(result.is_err(), "Test for '{}' should have failed", desc);
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {}",
            desc,
            expected_kind,
            result.clone().unwrap_err().map(|e| e.code),
        );
    }

    #[rstest]
    #[case(b"  1.2345  ", 1.2345)]
    #[case(b"    -1.234", -1.234)]
    #[case(b"1.0       ", 1.0)]
    #[case(b"1.        ", 1.0)]
    #[case(b"1.23456   ", 1.23456)]
    #[case(b"   1234567", 123.4567)]
    #[case(b"  -1234567", -123.4567)]
    #[case(b"       123", 0.0123)]
    #[case(b"          ", 0.0)]
    fn test_fixed_width_float(#[case] input: &[u8], #[case] expected: f64) {
        let mut parser = all_consuming(fixed_width_float::<f64>(10, 4)); // precision is ignored here
        let result = parser.parse(input);
        let (_, parsed_val) = result.unwrap();
        assert!((parsed_val - expected).abs() < 1e-9);
    }

    #[rstest]
    #[case(b"1.23a     ", "trailing characters", ErrorKind::Eof)]
    #[case(b"1.2.3     ", "invalid decimal point", ErrorKind::Eof)]
    #[case(b"          a", "trailing characters", ErrorKind::Eof)]
    fn test_fixed_width_float_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let mut parser = all_consuming(fixed_width_float::<f64>(10, 4));
        let result = parser.parse(input);
        assert!(result.is_err(), "Test for '{}' should have failed", desc);
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {}",
            desc,
            expected_kind,
            result.clone().unwrap_err().map(|e| e.code),
        );
    }
}
