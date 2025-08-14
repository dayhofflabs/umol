//! Parsing utilities for CTab files.

use std::{
    fmt::Debug,
    ops::{Range, RangeInclusive},
};

use umol_data::Element;

use fast_float::FastFloat;
use nom::branch::alt;
use nom::bytes::complete::{tag, take};
use nom::character::complete::{
    alpha1, digit0, digit1, i16 as nom_i16, i32 as nom_i32, i8 as nom_i8, space0, u32 as nom_u32,
    u8 as nom_u8, usize as nom_usize,
};
use nom::combinator::{complete, map, map_opt, map_res, opt, recognize, verify};
use nom::multi::{fold_many_m_n, separated_list1};
use nom::sequence::delimited;
use nom::{error, Err, Input, Parser};
use num::{Float, Integer};
use smallvec::{Array, SmallVec};

use crate::io::ctab::rgroup::RGroupOccurrence;

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
/// See `fixed_width_opt_partial` for more details.
pub(crate) fn fixed_width_opt<'a, O, P>(
    width: usize,
    inner: P,
) -> impl Parser<&'a [u8], Output = Option<O>, Error = error::Error<&'a [u8]>>
where
    P: Parser<&'a [u8], Output = O, Error = error::Error<&'a [u8]>>,
{
    fixed_width_partial(width, inner, false)
}

/// Parse a fixed-width field, making it optional.
///
/// This parser handles two cases of "optionality":
/// 1. The line is truncated, and the field is not present at all.
/// 2. The field is present but consists only of whitespace.
///
/// In both cases, it succeeds with `None`. If the field contains non-whitespace data, it runs
/// the `inner_parser`. If `inner_parser` fails on that data, it's a fatal error.
///
/// If `partial_ok` is true, the parser will succeed with `None` if the field is present but
/// consists only of whitespace. Otherwise, it will return an error.
pub(crate) fn fixed_width_partial<'a, O, P>(
    width: usize,
    mut inner: P,
    partial_ok: bool,
) -> impl Parser<&'a [u8], Output = Option<O>, Error = error::Error<&'a [u8]>>
where
    P: Parser<&'a [u8], Output = O, Error = error::Error<&'a [u8]>>,
{
    move |input: &'a [u8]| {
        let n_to_take = width.min(input.len());
        let (remaining, field) = take(n_to_take).parse(input)?;

        // If the slice is shorter than the expected width, it's only valid if it's all whitespace
        // or if `partial_ok` is true.
        if !partial_ok && field.len() < width && !field.iter().all(|&b| b == b' ') {
            return Err(Err::Error(error::Error::new(input, error::ErrorKind::Eof)));
        }

        if field.iter().all(|&b| b == b' ') {
            return Ok((remaining, None));
        }

        match inner.parse(field) {
            Ok((unconsumed, val)) => {
                if unconsumed.is_empty() {
                    Ok((remaining, Some(val)))
                } else {
                    Err(Err::Error(error::Error::new(input, error::ErrorKind::Eof)))
                }
            }
            Err(Err::Error(e)) => Err(Err::Error(error::Error::new(input, e.code))),
            Err(Err::Failure(e)) => Err(Err::Failure(error::Error::new(input, e.code))),
            Err(Err::Incomplete(needed)) => Err(Err::Incomplete(needed)),
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
    complete(verify(
        map(
            fixed_width_opt(width, delimited(space0, T::nom_parser(), space0)),
            |opt| opt.unwrap_or_else(T::zero),
        ),
        move |val: &T| range.contains(val),
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

/// Parse a fixed-width field as optional integer type. If range check fails, return None.
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

/// Parse a fixed-width field as integer, allow partial fields
pub(crate) fn fixed_width_int_partial<'a, T>(
    width: usize,
) -> impl Parser<&'a [u8], Output = T, Error = error::Error<&'a [u8]>>
where
    T: IntParser,
{
    complete(map(
        fixed_width_partial(width, delimited(space0, T::nom_parser(), space0), true),
        |opt| opt.unwrap_or_else(T::zero),
    ))
}

/// Parse a fixed-width field as float with Fortran semantics (Fw.d).
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

/// Parse a fixed-width field as element symbol, allow partial fields
pub(crate) fn fixed_width_element_partial<'a>(
    width: usize,
) -> impl Parser<&'a [u8], Output = Element, Error = error::Error<&'a [u8]>> {
    map_res(
        fixed_width_partial(
            width,
            delimited(space0, map_opt(alpha1, Element::from_symbol_bytes), space0),
            true,
        ),
        |opt| opt.ok_or_else(|| error::Error::new(&b""[..], error::ErrorKind::Verify)),
    )
}

/// Verify that a slice contains only blanks or zeros
pub(crate) fn is_blanks_or_zeros(input: &[u8]) -> bool {
    // Trim leading spaces.
    let mut start = 0;
    while start < input.len() && input[start] == b' ' {
        start += 1;
    }

    // If the slice is all spaces, it's valid.
    if start == input.len() {
        return true;
    }

    // Trim trailing spaces.
    let mut end = input.len();
    while end > start && input[end - 1] == b' ' {
        end -= 1;
    }

    // Check if the remaining (non-empty) slice is all zeros.
    let num_slice = &input[start..end];
    num_slice.iter().all(|&b| b == b'0')
}

/// Remove leading and trailing whitespace from a 2-character byte slice.
pub(crate) fn trim_whitespace_2char(s: &[u8]) -> &[u8] {
    debug_assert_eq!(s.len(), 2, "Input must be 2 characters");
    if s[0] == b' ' {
        if s[1] == b' ' {
            &b""[..]
        } else {
            &s[1..2]
        }
    } else if s[1] == b' ' {
        &s[0..1]
    } else {
        s
    }
}

/// Remove leading and trailing whitespace from a 3-character byte slice.
pub(crate) fn trim_whitespace_3char(s: &[u8]) -> &[u8] {
    debug_assert_eq!(s.len(), 3, "Input must be 3 characters");

    // Find start (skip leading whitespace)
    let start = if s[0] == b' ' {
        if s[1] == b' ' {
            if s[2] == b' ' {
                3 // All whitespace
            } else {
                2 // First two are whitespace
            }
        } else {
            1 // Only first is whitespace
        }
    } else {
        0 // No leading whitespace
    };

    if start == 3 {
        return &s[3..3]; // Empty slice - all whitespace
    }

    // Find end (skip trailing whitespace)
    let end = if s[2] == b' ' {
        if s[1] == b' ' {
            1 // Last two are whitespace
        } else {
            2 // Only last is whitespace
        }
    } else {
        3 // No trailing whitespace
    };

    // Ensure end is never less than start
    let end = end.max(start);

    &s[start..end]
}

/// Remove leading and trailing whitespace
#[allow(dead_code)]
pub(crate) fn trim_whitespace(input: &[u8]) -> &[u8] {
    match input.len() {
        0 => input,
        1 => {
            if input[0] == b' ' {
                &input[1..1]
            } else {
                input
            }
        }
        2 => trim_whitespace_2char(input),
        3 => trim_whitespace_3char(input),
        _ => {
            let mut start = 0;
            while start < input.len() && input[start] == b' ' {
                start += 1;
            }
            let mut end = input.len();
            while end > start && input[end - 1] == b' ' {
                end -= 1;
            }
            &input[start..end]
        }
    }
}

/// Apply the parser `p` exactly `n` times, discarding the results.
pub(crate) fn repeat<I, O, P>(
    n: usize,
    p: P,
) -> impl Parser<I, Output = (), Error = error::Error<I>>
where
    I: Input,
    P: Parser<I, Output = O, Error = error::Error<I>>,
{
    fold_many_m_n(n, n, p, || (), |_, _| ())
}

/// SmallVec-based parser combinator for length_count expressions.
#[allow(dead_code)]
pub(crate) fn small_length_count<I, A, C, F>(
    mut count: C,
    mut f: F,
) -> impl Parser<I, Output = SmallVec<A>, Error = error::Error<I>>
where
    I: nom::Input,
    A: Array,
    C: Parser<I, Output = usize, Error = error::Error<I>>,
    F: Parser<I, Output = <A as Array>::Item, Error = error::Error<I>>,
{
    move |input: I| {
        let (remaining, count) = count.parse(input)?;
        let mut v = SmallVec::new();

        let mut input = remaining;
        for _ in 0..count {
            let (remaining, val) = f.parse(input)?;
            v.push(val);
            input = remaining;
        }
        Ok((input, v))
    }
}

/// Parse a single RGroup occurrence.
pub(crate) fn rgroup_occurrence<'a>(
) -> impl Parser<&'a [u8], Output = RGroupOccurrence, Error = error::Error<&'a [u8]>> {
    alt((
        map((nom_u8, tag("-"), nom_u8), |(n, _, m)| {
            RGroupOccurrence::Range(n, m)
        }),
        map(nom_u8, RGroupOccurrence::Exactly),
        map((tag(">"), nom_u8), |(_, n)| {
            RGroupOccurrence::GreaterThan(n)
        }),
        map((tag("<"), nom_u8), |(_, n)| RGroupOccurrence::FewerThan(n)),
    ))
}

/// Parse a comma-separated list of RGroup occurrences.
pub(crate) fn rgroup_occurrences<'a>(
) -> impl Parser<&'a [u8], Output = Vec<RGroupOccurrence>, Error = error::Error<&'a [u8]>> {
    alt((
        map(
            delimited(
                space0,
                separated_list1(tag(","), rgroup_occurrence()),
                space0,
            ),
            |occurrences| occurrences,
        ),
        map(tag(""), |_| vec![RGroupOccurrence::GreaterThan(0)]),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nom::bytes::complete::take;
    use nom::combinator::{all_consuming, map_parser};
    use nom::{error, Err};
    use pretty_assertions::assert_eq;
    use rstest::*;
    use smallvec::smallvec;

    #[rstest]
    #[case(b"", None)]
    #[case(b"  ", None)]
    #[case(b"   ", None)]
    #[case(b"42", Some(42))]
    #[case(b" 42", Some(42))]
    #[case(b"42 ", Some(42))]
    #[case(b"042", Some(42))]
    fn test_fixed_width_partial(#[case] input: &[u8], #[case] expected_val: Option<i32>) {
        let mut parser = fixed_width_partial(3, delimited(space0, nom_i32, space0), true);

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
        assert!(remaining.is_empty(), "remaining should be empty");
    }

    #[rstest]
    #[case(b"abc", "non-numeric input", error::ErrorKind::Digit)]
    #[case(b"1a ", "trailing characters", error::ErrorKind::Eof)]
    fn test_fixed_width_partial_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
    ) {
        let mut parser = fixed_width_partial(3, delimited(space0, nom_i32, space0), true);
        let result = parser.parse(input);
        assert!(result.is_err(), "{} should have failed", desc);
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for '{}', expected {:?}, got {:?}",
            String::from_utf8_lossy(input),
            expected_kind,
            result.clone().unwrap_err().map(|e| e.code)
        );
    }

    #[rstest]
    #[case(b"", None)]
    #[case(b"  ", None)]
    #[case(b"   ", None)]
    #[case(b" 42", Some(42))]
    fn test_fixed_width_opt(#[case] input: &[u8], #[case] expected_val: Option<i32>) {
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
        assert!(remaining.is_empty(), "remaining should be empty");
    }

    #[rstest]
    #[case(b" abc ", "non-numeric input", error::ErrorKind::Digit)]
    #[case(b" 1", "too few characters", error::ErrorKind::Eof)]
    #[case(b"1", "too few characters", error::ErrorKind::Eof)]
    fn test_fixed_width_opt_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
    ) {
        let mut parser = fixed_width_opt(5, delimited(space0, nom_i32, space0));
        let result = parser.parse(input);
        assert!(
            result.is_err(),
            "{} should have failed with {:?}",
            desc,
            result.clone().unwrap_err().map(|e| e.code),
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
    #[case(b"1234", "too many characters", error::ErrorKind::Eof)]
    #[case(b"12", "too few characters", error::ErrorKind::Eof)]
    #[case(b"abc", "non-numeric input", error::ErrorKind::Digit)]
    #[case(b"1a ", "trailing characters", error::ErrorKind::Eof)]
    fn test_fixed_width_int_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
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
    #[case(b"   ", "blank field not in range", error::ErrorKind::Verify)]
    #[case(b"11 ", "value is out of range", error::ErrorKind::Verify)]
    #[case(b"1234", "too many characters", error::ErrorKind::Verify)]
    #[case(b"8", "too few characters", error::ErrorKind::Eof)]
    #[case(b"abc", "non-numeric input", error::ErrorKind::Digit)]
    #[case(b"1a ", "trailing characters", error::ErrorKind::Eof)]
    fn test_fixed_width_int_in_range_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
    ) {
        let mut parser = all_consuming(fixed_width_int_in_range::<i8, _>(3, 1i8..=10i8));
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
    #[case(b"abc", "non-numeric input", error::ErrorKind::Digit)]
    #[case(b"1a ", "trailing characters", error::ErrorKind::Eof)]
    fn test_fixed_width_int_in_range_opt_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
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
    #[case(b"123", 123i32)]
    #[case(b"12", 12i32)]
    #[case(b"1", 1i32)]
    #[case(b"", 0i32)]
    #[case(b"12 ", 12i32)]
    #[case(b"1  ", 1i32)]
    #[case(b" 12", 12i32)]
    #[case(b"  1", 1i32)]
    #[case(b" 1 ", 1i32)]
    #[case(b"   ", 0i32)]
    #[case(b"  ", 0i32)]
    #[case(b" ", 0i32)]
    #[case(b" -1", -1i32)]
    fn test_fixed_width_int_partial(#[case] input: &[u8], #[case] expected: i32) {
        let mut parser = all_consuming(fixed_width_int_partial::<i32>(3));
        let result = parser.parse(input);
        assert!(
            result.is_ok(),
            "Test for '{}' should have succeeded but failed with {:?}",
            String::from_utf8_lossy(input),
            result
        );
        let (remaining, result) = result.unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(result, expected);
    }

    #[rstest]
    #[case(b"1234", "too many characters", error::ErrorKind::Eof)]
    #[case(b"abc", "non-numeric input", error::ErrorKind::Digit)]
    #[case(b"1a ", "trailing characters", error::ErrorKind::Eof)]
    fn test_fixed_width_int_partial_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
    ) {
        let mut parser = all_consuming(fixed_width_int_partial::<i32>(3));
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
    #[case(b"C", Element::C)]
    #[case(b" C", Element::C)]
    #[case(b"Cu", Element::Cu)]
    #[case(b" Cu", Element::Cu)]
    #[case(b"  Cu", Element::Cu)]
    #[case(b" Cu ", Element::Cu)]
    #[case(b"Cu  ", Element::Cu)]
    #[case(b"   C", Element::C)]
    #[case(b"  C ", Element::C)]
    #[case(b" C  ", Element::C)]
    #[case(b"C   ", Element::C)]
    fn test_fixed_width_element_partial(#[case] input: &[u8], #[case] expected: Element) {
        let mut parser = all_consuming(fixed_width_element_partial(4));
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
    #[case(b"Cu   ", "trailing characters", error::ErrorKind::Eof)]
    #[case(b" X  ", "invalid element symbol", error::ErrorKind::MapOpt)]
    fn test_fixed_width_element_partial_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
    ) {
        let mut parser = all_consuming(fixed_width_element_partial(4));
        let result = parser.parse(input);
        assert!(
            result.is_err(),
            "{} should have failed",
            String::from_utf8_lossy(input)
        );
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {:?}",
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
    #[case(b"1.23a     ", "trailing characters", error::ErrorKind::Eof)]
    #[case(b"1.2.3     ", "invalid decimal point", error::ErrorKind::Eof)]
    #[case(b"          a", "trailing characters", error::ErrorKind::Eof)]
    fn test_fixed_width_float_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
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
    #[case(b"  0", "value too small", error::ErrorKind::Verify)]
    fn test_fixed_width_int_minus1_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
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

    #[rstest]
    #[case(b"", true)]
    #[case(b"   ", true)]
    #[case(b"  0", true)]
    #[case(b" 0 ", true)]
    #[case(b"0  ", true)]
    #[case(b" 00", true)]
    #[case(b"00 ", true)]
    #[case(b"000", true)]
    #[case(b"0", true)]
    #[case(b"00", true)]
    #[case(b"0 0", false)]
    #[case(b"  1", false)]
    fn test_is_blanks_or_zeros(#[case] input: &[u8], #[case] expected: bool) {
        assert_eq!(is_blanks_or_zeros(input), expected);
    }

    #[rstest]
    #[case(b"  ", &b""[..])]
    #[case(b" 0", &b"0"[..])]
    #[case(b"0 ", &b"0"[..])]
    #[case(b"00", &b"00"[..])]
    fn test_trim_whitespace_2char(#[case] input: &[u8], #[case] expected: &[u8]) {
        assert_eq!(trim_whitespace_2char(input), expected);
    }

    #[rstest]
    #[case(b"   ", &b""[..])]
    #[case(b"  0", &b"0"[..])]
    #[case(b" 0 ", &b"0"[..])]
    #[case(b"0  ", &b"0"[..])]
    #[case(b" 00", &b"00"[..])]
    #[case(b"00 ", &b"00"[..])]
    #[case(b"000", &b"000"[..])]
    fn test_trim_whitespace_3char(#[case] input: &[u8], #[case] expected: &[u8]) {
        assert_eq!(trim_whitespace_3char(input), expected);
    }

    #[rstest]
    #[case(b"", &b""[..])]
    #[case(b" ", &b""[..])]
    #[case(b"0", &b"0"[..])]
    #[case(b"  ", &b""[..])]
    #[case(b" 0", &b"0"[..])]
    #[case(b"   ", &b""[..])]
    #[case(b"  0", &b"0"[..])]
    #[case(b"  00", &b"00"[..])]
    #[case(b"  000", &b"000"[..])]
    fn test_trim_whitespace(#[case] input: &[u8], #[case] expected: &[u8]) {
        assert_eq!(trim_whitespace(input), expected);
    }

    #[rstest]
    #[case(b"abcabcabc")]
    fn test_repeat(#[case] input: &[u8]) {
        let mut parser = all_consuming(repeat::<_, _, _>(3, tag(&b"abc"[..])));
        let result = parser.parse(input);
        assert!(
            result.is_ok(),
            "{:?}: should have succeeded",
            String::from_utf8_lossy(input)
        );
        let (remaining, _) = result.unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
    }

    #[rstest]
    #[case(b"abcabc", error::ErrorKind::Tag)]
    #[case(b"abc", error::ErrorKind::Tag)]
    #[case(b"", error::ErrorKind::Tag)]
    fn test_repeat_invalid(#[case] input: &[u8], #[case] expected_kind: error::ErrorKind) {
        let mut parser = all_consuming(repeat::<_, _, _>(3, tag(&b"abc"[..])));
        let result = parser.parse(input);
        assert!(result.is_err(), "{:?}: should have failed", input);
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {:?}, expected {:?}, got {:?}",
            input,
            expected_kind,
            result.clone().unwrap_err().map(|e| e.code),
        );
    }

    #[rstest]
    #[case(b"3123", smallvec![1, 2, 3], "three items")]
    #[case(b"0", smallvec![], "zero count")]
    fn test_small_length_count(
        #[case] input: &[u8],
        #[case] expected_val: SmallVec<[u8; 8]>,
        #[case] desc: &str,
    ) {
        let count = map_parser(take(1u8), nom_usize);
        let item = map_parser(take(1u8), nom_u8);
        let mut parser = small_length_count::<_, [u8; 8], _, _>(count, item);
        let result = parser.parse(input);

        assert!(result.is_ok(), "{}: should have succeeded", desc);
        let (remaining, val) = result.unwrap();
        assert_eq!(val, expected_val, "{}: value mismatch", desc);
        assert!(remaining.is_empty(), "remaining should be empty");
    }

    #[rstest]
    #[case(b"312", error::ErrorKind::Eof, "incomplete items")]
    #[case(b"x12", error::ErrorKind::Digit, "invalid count character")]
    fn test_small_length_count_invalid(
        #[case] input: &[u8],
        #[case] expected_kind: error::ErrorKind,
        #[case] desc: &str,
    ) {
        let count = map_parser(take(1u8), nom_usize);
        let item = map_parser(take(1u8), nom_u8);
        let mut parser = small_length_count::<_, [u8; 8], _, _>(count, item);
        let result = parser.parse(input);

        assert!(result.is_err(), "{}: should have failed", desc);
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {}",
            desc,
            expected_kind,
            result.clone().unwrap_err().map(|e| e.code),
        );
    }

    #[rstest]
    #[case(b"", vec![RGroupOccurrence::GreaterThan(0)])]
    #[case(b"1", vec![RGroupOccurrence::Exactly(1)])]
    #[case(b"1,2", vec![RGroupOccurrence::Exactly(1), RGroupOccurrence::Exactly(2)])]
    #[case(b">1", vec![RGroupOccurrence::GreaterThan(1)])]
    #[case(b"<2", vec![RGroupOccurrence::FewerThan(2)])]
    #[case(b"1-3", vec![RGroupOccurrence::Range(1, 3)])]
    #[case(b"0,>0", vec![RGroupOccurrence::Exactly(0), RGroupOccurrence::GreaterThan(0)])]
    fn test_rgroup_occurrences(#[case] input: &[u8], #[case] expected: Vec<RGroupOccurrence>) {
        let mut parser = all_consuming(rgroup_occurrences());
        let result = parser.parse(input);
        assert!(
            result.is_ok(),
            "{}: should have succeeded",
            String::from_utf8_lossy(input)
        );
        let (remaining, val) = result.unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
        assert_eq!(
            val,
            expected,
            "{}: value mismatch",
            String::from_utf8_lossy(input)
        );
    }

    #[rstest]
    #[case(b"a", "invalid character", error::ErrorKind::Eof)]
    #[case(b"-3", "negative value", error::ErrorKind::Eof)]
    fn test_rgroup_occurrences_invalid(
        #[case] input: &[u8],
        #[case] desc: &str,
        #[case] expected_kind: error::ErrorKind,
    ) {
        let mut parser = all_consuming(rgroup_occurrences());
        let result = parser.parse(input);
        assert!(
            result.is_err(),
            "{}: should have failed",
            String::from_utf8_lossy(input)
        );
        assert!(
            matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
            "Mismatched error kind for {}, expected {:?}, got {}",
            desc,
            expected_kind,
            result.clone().unwrap_err().map(|e| e.code),
        );
    }
}
