//! Parsing utilities for CTab files.

use std::{
    fmt::Debug,
    ops::{Range, RangeInclusive},
};

use fast_float::FastFloat;
use nom::{
    branch::alt,
    bytes::{complete::tag, take},
    character::complete::{
        digit0, digit1, i16 as nom_i16, i32 as nom_i32, i8 as nom_i8, space0, u32 as nom_u32,
        u8 as nom_u8, usize as nom_usize,
    },
    combinator::{complete, map, opt, recognize, verify},
    error,
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

impl IntParser for i16 {
    fn nom_parser<'a>() -> impl Parser<&'a [u8], Output = Self, Error = error::Error<&'a [u8]>> {
        nom_i16
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

impl IntParser for u32 {
    fn nom_parser<'a>() -> impl Parser<&'a [u8], Output = Self, Error = error::Error<&'a [u8]>> {
        nom_u32
    }
}

impl IntParser for usize {
    fn nom_parser<'a>() -> impl Parser<&'a [u8], Output = Self, Error = error::Error<&'a [u8]>> {
        nom_usize
    }
}

/// Parse a fixed-width field, making it optional.
///
/// This is the foundational parser for fixed-width fields. It handles two cases of "optionality":
/// 1. The line is truncated, and the field is not present at all.
/// 2. The field is present but consists only of whitespace.
///
/// In both cases, it succeeds with `None`. If the field contains non-whitespace data, it runs
/// the `inner_parser`. If `inner_parser` fails on that data, it's a fatal error, which is the
/// desired behavior for malformed data.
pub(crate) fn fixed_width_opt<'a, O, P>(
    width: usize,
    mut inner: P,
) -> impl Parser<&'a [u8], Output = Option<O>, Error = error::Error<&'a [u8]>>
where
    P: Parser<&'a [u8], Output = O, Error = error::Error<&'a [u8]>>,
{
    move |input: &'a [u8]| {
        let n_to_take = width.min(input.len());
        let (remaining, field_slice) = take(n_to_take).parse(input)?;

        // If the slice is shorter than the expected width, it's only valid if it's all whitespace.
        if field_slice.len() < width && !field_slice.iter().all(|&b| b == b' ') {
            return Err(nom::Err::Error(error::Error::new(
                input,
                nom::error::ErrorKind::Eof,
            )));
        }

        if field_slice.iter().all(|&b| b == b' ') {
            return Ok((remaining, None));
        }

        match inner.parse(field_slice) {
            Ok((unconsumed, val)) => {
                if unconsumed.is_empty() {
                    Ok((remaining, Some(val)))
                } else {
                    Err(nom::Err::Error(error::Error::new(
                        input,
                        nom::error::ErrorKind::Eof,
                    )))
                }
            }
            Err(nom::Err::Error(e)) => Err(nom::Err::Error(error::Error::new(input, e.code))),
            Err(nom::Err::Failure(e)) => Err(nom::Err::Failure(error::Error::new(input, e.code))),
            Err(nom::Err::Incomplete(needed)) => Err(nom::Err::Incomplete(needed)),
        }
    }
}

/// Parse a fixed-width field as an integer type. Interprets empty/whitespace field as default.
pub(crate) fn fixed_width_int<'a, T>(
    width: usize,
) -> impl Parser<&'a [u8], Output = T, Error = error::Error<&'a [u8]>>
where
    T: IntParser,
{
    complete(map(
        fixed_width_opt(width, delimited(space0, T::nom_parser(), space0)),
        |opt| opt.unwrap_or_else(T::zero),
    ))
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
    complete(map(
        verify(
            fixed_width_opt(width, delimited(space0, T::nom_parser(), space0)),
            move |opt: &Option<T>| opt.map_or(true, |val| range.contains(&val)),
        ),
        |opt| opt.unwrap_or_else(T::zero),
    ))
}

/// Parse a fixed-width field as an integer type, subtracting one.
pub(crate) fn fixed_width_int_minus1<'a, T>(
    width: usize,
) -> impl Parser<&'a [u8], Output = T, Error = error::Error<&'a [u8]>>
where
    T: IntParser,
{
    map(
        verify(fixed_width_int(width), |val: &T| *val >= T::one()),
        |x: T| x - T::one(),
    )
}

/// Parse a fixed-width field as an integer type with a range check, subtracting one.
pub(crate) fn fixed_width_int_in_range_minus1<'a, T, R>(
    width: usize,
    range: R,
) -> impl Parser<&'a [u8], Output = T, Error = error::Error<&'a [u8]>>
where
    T: IntParser,
    R: Contains<T> + Clone + Debug,
{
    map(
        verify(fixed_width_int_in_range(width, range), |val: &T| {
            *val >= T::one()
        }),
        |x: T| x - T::one(),
    )
}

/// Parse a fixed-width field as an optional integer type. If range check fails, return None.
pub(crate) fn fixed_width_int_in_range_opt<'a, T, R>(
    width: usize,
    range: R,
) -> impl Parser<&'a [u8], Output = Option<T>, Error = error::Error<&'a [u8]>>
where
    T: IntParser,
    R: Contains<T> + Clone,
{
    map(
        fixed_width_opt(width, delimited(space0, T::nom_parser(), space0)),
        move |opt| opt.filter(|val| range.contains(val)),
    )
}

/// Parse a fixed-width field as a float type with Fortran semantics (Fw.d).
pub(crate) fn fixed_width_float<'a, T>(
    width: usize,
    precision: usize,
) -> impl Parser<&'a [u8], Output = T, Error = error::Error<&'a [u8]>>
where
    T: Float + FastFloat,
{
    let value_parser = alt((
        map(
            recognize((opt(tag(&b"-"[..])), digit1, tag(&b"."[..]), digit0)),
            |s| fast_float::parse::<T, _>(s).unwrap(),
        ),
        map(recognize((opt(tag(&b"-"[..])), digit1)), move |s| {
            fast_float::parse::<T, _>(s).unwrap() / T::from(10.0).unwrap().powi(precision as i32)
        }),
    ));
    complete(map(
        fixed_width_opt(width, delimited(space0, value_parser, space0)),
        |opt| opt.unwrap_or_else(T::zero),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nom::combinator::all_consuming;
    use nom::error::ErrorKind;
    use nom::Err;
    use rstest::rstest;

    #[rstest]
    // Field absent (line is empty, a subset of too short)
    #[case(b"", None, b"")]
    // Field absent (line is too short but all whitespace)
    #[case(b"  ", None, b"")]
    // Field present but blank
    #[case(b"   ", None, b"")]
    // Field present and parsable
    #[case(b" 42", Some(42), b"")]
    // Field present and parsable with remaining data
    #[case(b" 42 leftover", Some(42), b" leftover")]
    fn test_fixed_width_opt(
        #[case] input: &[u8],
        #[case] expected_val: Option<i32>,
        #[case] expected_rem: &[u8],
    ) {
        let mut parser = fixed_width_opt(3, delimited(space0, nom_i32, space0));
        let result = parser.parse(input);
        assert!(
            result.is_ok(),
            "Test for '{}' should have succeeded but failed with {:?}",
            String::from_utf8_lossy(input),
            result
        );
        let (remaining, value) = result.unwrap();
        assert_eq!(
            value,
            expected_val,
            "Mismatched value for '{}'",
            String::from_utf8_lossy(input)
        );
        assert_eq!(
            remaining,
            expected_rem,
            "Mismatched remaining for '{}'",
            String::from_utf8_lossy(input)
        );
    }

    #[rstest]
    // Field present but not parsable
    #[case(b" abc ", "abc", ErrorKind::Digit)]
    // Field partially present (too short) and not all whitespace -> ambiguous
    #[case(b" 1", "1", ErrorKind::Eof)]
    // Field partially present (one char) and not all whitespace -> ambiguous
    #[case(b"1", "1", ErrorKind::Eof)]
    fn test_fixed_width_opt_invalid(
        #[case] input: &[u8],
        #[case] _desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let mut parser = fixed_width_opt(5, delimited(space0, nom_i32, space0));
        let result = parser.parse(input);
        assert!(
            result.is_err(),
            "Test for '{}' should have failed",
            String::from_utf8_lossy(input)
        );
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for '{}', expected {:?}, got {:?}",
            String::from_utf8_lossy(input),
            expected_kind,
            result
        );
    }

    #[rstest]
    #[case(b"123", 123i32)]
    #[case(b"-98", -98i32)]
    #[case(b"  8", 8i32)]
    #[case(b"   ", 0i32)]
    fn test_fixed_width_int(#[case] input: &[u8], #[case] expected: i32) {
        let mut parser = all_consuming(fixed_width_int::<i32>(3));
        let result = parser.parse(input);
        assert!(
            result.is_ok(),
            "Test for '{}' should have succeeded",
            String::from_utf8_lossy(input)
        );
        let (remaining, result) = result.unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"1234", "too many characters", ErrorKind::Eof)]
    #[case(b"12", "too few characters", ErrorKind::Eof)]
    #[case(b"abc", "non-numeric input", ErrorKind::Digit)]
    #[case(b"1a ", "trailing characters", ErrorKind::Eof)]
    fn test_fixed_width_int_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let mut parser = all_consuming(fixed_width_int::<i32>(3));
        let result = parser.parse(input);
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
    #[case(b"100", 100i8)]
    #[case(b" -9", -9i8)]
    #[case(b"8  ", 8i8)]
    #[case(b" 1 ", 1i8)]
    fn test_fixed_width_int_in_range(#[case] input: &[u8], #[case] expected: i8) {
        let mut parser = all_consuming(fixed_width_int_in_range::<i8, _>(3, -10i8..=110i8));
        let result = parser.parse(input);
        assert!(
            result.is_ok(),
            "{} should have succeeded",
            String::from_utf8_lossy(input)
        );
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
    #[case(b"100", 100u8)]
    #[case(b"  9", 9u8)]
    #[case(b"8  ", 8u8)]
    #[case(b" 1 ", 1u8)]
    fn test_fixed_width_int_in_range_inclusive(#[case] input: &[u8], #[case] expected: u8) {
        let mut parser = all_consuming(fixed_width_int_in_range::<u8, _>(3, 0u8..=100u8));
        let result = parser.parse(input);
        assert!(
            result.is_ok(),
            "Test for '{}' should have succeeded",
            String::from_utf8_lossy(input)
        );
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
        assert!(
            result.is_ok(),
            "Test for '{}' should have succeeded",
            String::from_utf8_lossy(input)
        );
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
        assert!(
            result.is_ok(),
            "Test for '{}' should have succeeded",
            String::from_utf8_lossy(input)
        );
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
    #[case(b"  5", Some(5i32))]
    #[case(b" 10", Some(10i32))]
    #[case(b"  0", Some(0i32))]
    #[case(b" 11", None)]
    #[case(b" -1", None)]
    #[case(b"   ", None)]
    #[case(b"  ", None)]
    #[case(b"", None)]
    fn test_fixed_width_int_in_range_opt(#[case] input: &[u8], #[case] expected: Option<i32>) {
        let mut parser = all_consuming(fixed_width_int_in_range_opt::<i32, _>(3, 0..=10));
        let result = parser.parse(input);
        assert!(
            result.is_ok(),
            "Test for '{}' should have succeeded but failed with {:?}",
            String::from_utf8_lossy(input),
            result
        );
        let (remaining, value) = result.unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(value, expected);
    }

    #[rstest]
    #[case(b"abc", "non-numeric input", ErrorKind::Digit)]
    #[case(b"1a ", "trailing characters", ErrorKind::Eof)]
    fn test_fixed_width_int_in_range_opt_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let mut parser = all_consuming(fixed_width_int_in_range_opt::<i32, _>(3, 0..=10));
        let result = parser.parse(input);
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
    #[case(b"  0", "value too small", ErrorKind::Verify)]
    fn test_fixed_width_int_minus1_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: ErrorKind,
    ) {
        let mut parser = all_consuming(fixed_width_int_minus1::<usize>(3));
        let result = parser.parse(input);
        assert!(result.is_err(), "{} should have failed", desc);
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {}",
            desc,
            expected_kind,
            result.clone().unwrap_err().map(|e| e.code),
        );
    }
}
