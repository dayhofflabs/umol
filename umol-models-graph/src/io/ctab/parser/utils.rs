//! Parsing utilities for CTab files.

use std::fmt::Debug;
use std::ops::{Range, RangeInclusive};

use bstr::ByteSlice;
use umol_data::Element;

use fast_float::FastFloat;
use nom::branch::alt;
use nom::bytes::complete::{tag, take};
use nom::character::complete::{
    alpha1, digit0, digit1, i16 as nom_i16, i32 as nom_i32, i8 as nom_i8, u32 as nom_u32,
    u8 as nom_u8, usize as nom_usize,
};
use nom::combinator::{complete, map, map_opt, opt, recognize, rest, verify};
use nom::multi::{count as nom_count, separated_list1};
use nom::{error, Err, Parser};
use num::{Float, Integer};

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

/// Trim whitespace
pub(crate) fn trim_whitespace(input: &[u8], allow_unicode: bool) -> &[u8] {
    if allow_unicode {
        input.trim()
    } else {
        input.trim_ascii()
    }
}

/// Check if all bytes in slice are whitespace characters
pub(crate) fn is_all_whitespace(input: &[u8], allow_unicode: bool) -> bool {
    if allow_unicode {
        use bstr::ByteSlice;
        input.chars().all(|c| c.is_whitespace())
    } else {
        input.iter().all(|&b| matches!(b, b' ' | b'\t'))
    }
}

/// Verify that a slice contains only whitespace or zeroes
pub(crate) fn is_all_whitespace_or_zeroes(input: &[u8], allow_unicode: bool) -> bool {
    let input = if allow_unicode {
        input.trim()
    } else {
        input.trim_ascii()
    };
    input.find_not_byteset(b"0").is_none()
}

/// Convert byte slice to string, trimming leading and trailing whitespace
pub(crate) fn to_string(bytes: &[u8], allow_unicode: bool) -> Result<String, error::Error<&[u8]>> {
    let s = trim_whitespace(bytes, allow_unicode);
    if !allow_unicode && s.find_non_ascii_byte().is_some() {
        Err(error::Error::new(bytes, error::ErrorKind::Verify))
    } else {
        Ok(s.to_str_lossy().into_owned())
    }
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
    allow_unicode: bool,
) -> impl Parser<&'a [u8], Output = Option<O>, Error = error::Error<&'a [u8]>>
where
    P: Parser<&'a [u8], Output = O, Error = error::Error<&'a [u8]>>,
{
    move |input: &'a [u8]| {
        let min_width = width.min(input.len());
        let (remaining, field) = take(min_width).parse(input)?;

        // If the slice is shorter than the expected width, it's only valid if it's all whitespace
        // or if `partial_ok` is true.
        if field.len() < width && !partial_ok && !is_all_whitespace(field, allow_unicode) {
            return Err(Err::Error(error::Error::new(input, error::ErrorKind::Eof)));
        }

        if is_all_whitespace(field, allow_unicode) {
            return Ok((remaining, None));
        }

        match inner.parse(field) {
            Ok((remaining_inner, val)) => {
                if remaining_inner.is_empty() {
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

/// Parse an optional fixed-width field. If the field is present but consists only of whitespace,
/// it succeeds with `None`. Otherwise, it runs the `inner` parser. Partial fields are not allowed.
pub(crate) fn fixed_width_opt<'a, O, P>(
    width: usize,
    inner: P,
    allow_unicode: bool,
) -> impl Parser<&'a [u8], Output = Option<O>, Error = error::Error<&'a [u8]>>
where
    P: Parser<&'a [u8], Output = O, Error = error::Error<&'a [u8]>>,
{
    fixed_width_partial(width, inner, false, allow_unicode)
}

/// Parse a fixed-width field as an integer type. Interprets empty/whitespace field as default.
pub(crate) fn fixed_width_int<'a, T>(
    width: usize,
    allow_unicode: bool,
) -> impl Parser<&'a [u8], Output = T, Error = error::Error<&'a [u8]>>
where
    T: IntParser,
{
    complete(map(
        fixed_width_opt(
            width,
            move |input| {
                let s = trim_whitespace(input, allow_unicode);
                T::nom_parser().parse(s)
            },
            allow_unicode,
        ),
        |opt| opt.unwrap_or_else(T::zero),
    ))
}

/// Parse a fixed-width field as an integer type, applying range bounds.
pub(crate) fn fixed_width_int_in_range<'a, T, R>(
    width: usize,
    range: R,
    allow_unicode: bool,
) -> impl Parser<&'a [u8], Output = T, Error = error::Error<&'a [u8]>>
where
    T: IntParser,
    R: Contains<T> + Clone,
{
    complete(verify(
        map(
            fixed_width_opt(
                width,
                move |input| {
                    let s = trim_whitespace(input, allow_unicode);
                    T::nom_parser().parse(s)
                },
                allow_unicode,
            ),
            |opt| opt.unwrap_or_else(T::zero),
        ),
        move |val: &T| range.contains(val),
    ))
}

/// Parse a fixed-width field as optional integer type. If range check fails, return None.
pub(crate) fn fixed_width_int_in_range_opt<'a, T, R>(
    width: usize,
    range: R,
    allow_unicode: bool,
) -> impl Parser<&'a [u8], Output = Option<T>, Error = error::Error<&'a [u8]>>
where
    T: IntParser,
    R: Contains<T> + Clone,
{
    map(
        fixed_width_opt(
            width,
            move |input| {
                let s = trim_whitespace(input, allow_unicode);
                T::nom_parser().parse(s)
            },
            allow_unicode,
        ),
        move |opt| opt.filter(|val| range.contains(val)),
    )
}

/// Parse a fixed-width field as an integer type, subtracting one.
pub(crate) fn fixed_width_int_minus1<'a, T>(
    width: usize,
    allow_unicode: bool,
) -> impl Parser<&'a [u8], Output = T, Error = error::Error<&'a [u8]>>
where
    T: IntParser,
{
    map(
        verify(fixed_width_int(width, allow_unicode), |val: &T| {
            *val >= T::one()
        }),
        |x: T| x - T::one(),
    )
}

/// Parse a fixed-width field as integer, allow partial fields
pub(crate) fn fixed_width_int_partial<'a, T>(
    width: usize,
    allow_unicode: bool,
) -> impl Parser<&'a [u8], Output = T, Error = error::Error<&'a [u8]>>
where
    T: IntParser,
{
    complete(map(
        fixed_width_partial(
            width,
            move |input| {
                let s = trim_whitespace(input, allow_unicode);
                T::nom_parser().parse(s)
            },
            true,
            allow_unicode,
        ),
        |opt| opt.unwrap_or_else(T::zero),
    ))
}

/// Parse a fixed-width field as float with Fortran semantics (Fw.d).
pub(crate) fn fixed_width_float<'a, T>(
    width: usize,
    precision: usize,
    allow_unicode: bool,
) -> impl Parser<&'a [u8], Output = T, Error = error::Error<&'a [u8]>>
where
    T: Float + FastFloat,
{
    complete(map(
        fixed_width_opt(
            width,
            move |input| {
                let s = trim_whitespace(input, allow_unicode);
                alt((
                    map(
                        recognize((opt(tag(&b"-"[..])), digit1, tag(&b"."[..]), digit0)),
                        |s| fast_float::parse::<T, _>(s).unwrap(),
                    ),
                    map(recognize((opt(tag(&b"-"[..])), digit1)), move |s| {
                        fast_float::parse::<T, _>(s).unwrap()
                            / T::from(10.0).unwrap().powi(precision as i32)
                    }),
                ))
                .parse(s)
            },
            allow_unicode,
        ),
        |opt| opt.unwrap_or_else(T::zero),
    ))
}

/// Parse a fixed-width field as element symbol
pub(crate) fn fixed_width_element_partial<'a>(
    width: usize,
    allow_unicode: bool,
) -> impl Parser<&'a [u8], Output = Option<Element>, Error = error::Error<&'a [u8]>> {
    fixed_width_partial(
        width,
        move |s| {
            let s = trim_whitespace(s, allow_unicode);
            map_opt(alpha1, Element::from_symbol_bytes).parse(s)
        },
        true,
        allow_unicode,
    )
}

/// Padding field of fixed width `width`
/// If `allow_unicode` is true, allow unicode whitespace.
/// If `strict_padding` is true, require strict padding (only whitespace or zeroes).
pub(crate) fn fixed_width_padding<'a>(
    width: usize,
    allow_unicode: bool,
    strict_padding: bool,
) -> impl Parser<&'a [u8], Output = (), Error = error::Error<&'a [u8]>> {
    move |input: &'a [u8]| {
        let (remaining, padding) = take(width).parse(input)?;
        if strict_padding && width > 0 && !is_all_whitespace_or_zeroes(padding, allow_unicode) {
            Err(Err::Error(error::Error::new(
                input,
                error::ErrorKind::Verify,
            )))
        } else {
            Ok((remaining, ()))
        }
    }
}

/// Multiple fixed-width padding fields of width `width`
/// If `allow_unicode` is true, allow unicode whitespace.
/// If `strict_padding` is true, require strict padding (only whitespace or zeroes).
pub(crate) fn fixed_width_padding_n<'a>(
    count: usize,
    width: usize,
    allow_unicode: bool,
    strict_padding: bool,
) -> impl Parser<&'a [u8], Output = (), Error = error::Error<&'a [u8]>> {
    move |input: &'a [u8]| {
        let (remaining, padding) = take(count * width).parse(input)?;
        if count > 0 && width > 0 && strict_padding {
            nom_count(
                fixed_width_padding(width, allow_unicode, strict_padding),
                count,
            )
            .parse(padding)
            .map(|(_, _)| (remaining, ()))
        } else {
            Ok((remaining, ()))
        }
    }
}

/// Parse a fixed-width field as a string, allow partial fields
pub(crate) fn fixed_width_str_partial<'a>(
    width: usize,
    allow_unicode: bool,
) -> impl Parser<&'a [u8], Output = Option<String>, Error = error::Error<&'a [u8]>> {
    map(
        fixed_width_partial(width, rest, true, allow_unicode),
        move |opt| opt.and_then(|s| to_string(s, allow_unicode).ok()),
    )
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
    allow_unicode: bool,
) -> impl Parser<&'a [u8], Output = Vec<RGroupOccurrence>, Error = error::Error<&'a [u8]>> {
    alt((
        move |input: &'a [u8]| {
            let s = trim_whitespace(input, allow_unicode);
            separated_list1(tag(","), rgroup_occurrence()).parse(s)
        },
        map(tag(""), |_| vec![RGroupOccurrence::GreaterThan(0)]),
    ))
}

#[cfg(test)]
mod tests;
