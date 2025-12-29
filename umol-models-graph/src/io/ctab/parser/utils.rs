//! Parsing utilities for CTab files.

use std::fmt::Debug;
use std::ops::{Range, RangeInclusive};

use bstr::ByteSlice;
use fast_float::FastFloat;
use nom::branch::alt;
use nom::bytes::complete::{tag, take};
use nom::character::complete::{
    alpha1, digit0, digit1, i16 as nom_i16, i32 as nom_i32, i8 as nom_i8, space0, u32 as nom_u32,
    u8 as nom_u8, usize as nom_usize,
};
use nom::combinator::{map, map_opt, opt, recognize, rest, success, verify};
use nom::error::{Error as NomError, ErrorKind as NomErrorKind};
use nom::multi::{count as nom_count, separated_list1};
use nom::sequence::delimited;
use nom::{Err, Parser};
use num::{Float, Integer};
use umol_data::{Element, NamedIsotope};

use crate::position::Point3D;
use crate::table_ir::RGroupOccurrence;

pub(super) trait Contains<T: PartialOrd> {
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

pub(super) trait IntParser: Sized + Copy + PartialOrd + Debug + Default + Integer {
    fn nom_parser<'a>() -> impl Parser<&'a [u8], Output = Self, Error = NomError<&'a [u8]>>;
}

impl IntParser for i8 {
    fn nom_parser<'a>() -> impl Parser<&'a [u8], Output = Self, Error = NomError<&'a [u8]>> {
        nom_i8
    }
}

impl IntParser for i16 {
    fn nom_parser<'a>() -> impl Parser<&'a [u8], Output = Self, Error = NomError<&'a [u8]>> {
        nom_i16
    }
}

impl IntParser for i32 {
    fn nom_parser<'a>() -> impl Parser<&'a [u8], Output = Self, Error = NomError<&'a [u8]>> {
        nom_i32
    }
}

impl IntParser for u8 {
    fn nom_parser<'a>() -> impl Parser<&'a [u8], Output = Self, Error = NomError<&'a [u8]>> {
        nom_u8
    }
}

impl IntParser for u32 {
    fn nom_parser<'a>() -> impl Parser<&'a [u8], Output = Self, Error = NomError<&'a [u8]>> {
        nom_u32
    }
}

impl IntParser for usize {
    fn nom_parser<'a>() -> impl Parser<&'a [u8], Output = Self, Error = NomError<&'a [u8]>> {
        nom_usize
    }
}

/// Verify that a slice contains only whitespace or zeroes
pub(super) fn is_all_whitespace_or_zeroes(input: &[u8]) -> bool {
    input.trim_ascii().find_not_byteset(b"0").is_none()
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
pub(super) fn fixed_width_partial<'a, O, P>(
    width: usize,
    mut inner: P,
    partial_ok: bool,
) -> impl Parser<&'a [u8], Output = Option<O>, Error = NomError<&'a [u8]>>
where
    P: Parser<&'a [u8], Output = O, Error = NomError<&'a [u8]>>,
{
    move |input: &'a [u8]| {
        let min_width = width.min(input.len());
        let (remaining, field) = take(min_width).parse(input)?;

        // If the slice is shorter than the expected width, it's only valid if it's all whitespace
        // or if `partial_ok` is true.
        if field.len() < width && !partial_ok && field.find_not_byteset(b"  \t").is_some() {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Eof)));
        }

        if field.find_not_byteset(b"  \t").is_none() {
            return Ok((remaining, None));
        }

        match inner.parse(field) {
            Ok((remaining_inner, val)) => {
                if remaining_inner.is_empty() {
                    Ok((remaining, Some(val)))
                } else {
                    Err(Err::Error(NomError::new(input, NomErrorKind::Eof)))
                }
            }
            Err(Err::Error(e)) => Err(Err::Error(NomError::new(input, e.code))),
            Err(Err::Failure(e)) => Err(Err::Failure(NomError::new(input, e.code))),
            Err(Err::Incomplete(needed)) => Err(Err::Incomplete(needed)),
        }
    }
}

/// Parse an optional fixed-width field. If the field is present but consists only of whitespace,
/// it succeeds with `None`. Otherwise, it runs the `inner` parser. Partial fields are not allowed.
pub(super) fn fixed_width_opt<'a, O, P>(
    width: usize,
    inner: P,
) -> impl Parser<&'a [u8], Output = Option<O>, Error = NomError<&'a [u8]>>
where
    P: Parser<&'a [u8], Output = O, Error = NomError<&'a [u8]>>,
{
    fixed_width_partial(width, inner, false)
}

/// Parse a fixed-width field as an integer type. Interprets empty/whitespace field as default.
pub(super) fn fixed_width_int<'a, T>(
    width: usize,
) -> impl Parser<&'a [u8], Output = T, Error = NomError<&'a [u8]>>
where
    T: IntParser,
{
    map(
        fixed_width_opt(width, delimited(space0, T::nom_parser(), space0)),
        |opt| opt.unwrap_or_else(T::zero),
    )
}

/// Parse a fixed-width field as an integer type, applying range bounds.
pub(super) fn fixed_width_int_in_range<'a, T, R>(
    width: usize,
    range: R,
) -> impl Parser<&'a [u8], Output = T, Error = NomError<&'a [u8]>>
where
    T: IntParser,
    R: Contains<T> + Clone,
{
    verify(
        map(
            fixed_width_opt(width, delimited(space0, T::nom_parser(), space0)),
            |opt| opt.unwrap_or_else(T::zero),
        ),
        move |val: &T| range.contains(val),
    )
}

/// Parse a fixed-width field as optional integer type. If range check fails, return None.
pub(super) fn fixed_width_int_in_range_opt<'a, T, R>(
    width: usize,
    range: R,
) -> impl Parser<&'a [u8], Output = Option<T>, Error = NomError<&'a [u8]>>
where
    T: IntParser,
    R: Contains<T> + Clone,
{
    map(
        fixed_width_opt(width, delimited(space0, T::nom_parser(), space0)),
        move |opt| opt.filter(|val| range.contains(val)),
    )
}

/// Parse a fixed-width field as an integer type, subtracting one.
pub(super) fn fixed_width_int_minus1<'a, T>(
    width: usize,
) -> impl Parser<&'a [u8], Output = T, Error = NomError<&'a [u8]>>
where
    T: IntParser,
{
    map(
        verify(fixed_width_int(width), |val: &T| *val >= T::one()),
        |x: T| x - T::one(),
    )
}

/// Parse a fixed-width field as integer, allow partial fields
pub(super) fn fixed_width_int_partial<'a, T>(
    width: usize,
) -> impl Parser<&'a [u8], Output = T, Error = NomError<&'a [u8]>>
where
    T: IntParser,
{
    map(
        fixed_width_partial(width, delimited(space0, T::nom_parser(), space0), true),
        |opt| opt.unwrap_or_else(T::zero),
    )
}

/// Parse a fixed-width field as float with Fortran semantics (Fw.d).
pub(super) fn fixed_width_float<'a, T>(
    width: usize,
    precision: usize,
) -> impl Parser<&'a [u8], Output = T, Error = NomError<&'a [u8]>>
where
    T: Float + FastFloat,
{
    map(
        fixed_width_opt(
            width,
            delimited(
                space0,
                alt((
                    map(
                        recognize((opt(tag(&b"-"[..])), digit1, tag(&b"."[..]), digit0)),
                        |s| fast_float::parse::<T, _>(s).unwrap(),
                    ),
                    map(
                        recognize((opt(tag(&b"-"[..])), digit0, tag(&b"."[..]), digit1)),
                        |s| fast_float::parse::<T, _>(s).unwrap(),
                    ),
                    map(recognize((opt(tag(&b"-"[..])), digit1)), move |s| {
                        fast_float::parse::<T, _>(s).unwrap()
                            / T::from(10.0).unwrap().powi(precision as i32)
                    }),
                )),
                space0,
            ),
        ),
        |opt| opt.unwrap_or_else(T::zero),
    )
}

/// Parse a fixed-width field as element symbol
pub(super) fn fixed_width_element_partial<'a>(
    width: usize,
) -> impl Parser<&'a [u8], Output = Option<Element>, Error = NomError<&'a [u8]>> {
    fixed_width_partial(
        width,
        delimited(space0, map_opt(alpha1, Element::from_symbol_bytes), space0),
        true,
    )
}

/// Padding field of fixed width `width`
/// Only validate padding if `skip_padding` is false.
pub(super) fn fixed_width_padding<'a>(
    width: usize,
    skip_padding: bool,
) -> impl Parser<&'a [u8], Output = (), Error = NomError<&'a [u8]>> {
    move |input: &'a [u8]| {
        let (remaining, padding) = take(width).parse(input)?;
        if !skip_padding && width > 0 && !is_all_whitespace_or_zeroes(padding) {
            Err(Err::Error(NomError::new(input, NomErrorKind::Verify)))
        } else {
            Ok((remaining, ()))
        }
    }
}

/// Multiple fixed-width padding fields of width `width`
/// Only validate padding if `skip_padding` is false.
pub(super) fn fixed_width_padding_n<'a>(
    count: usize,
    width: usize,
    skip_padding: bool,
) -> impl Parser<&'a [u8], Output = (), Error = NomError<&'a [u8]>> {
    move |input: &'a [u8]| {
        let (remaining, padding) = take(count * width).parse(input)?;
        if !skip_padding && count > 0 && width > 0 {
            nom_count(fixed_width_padding(width, skip_padding), count)
                .parse(padding)
                .map(|(_, _)| (remaining, ()))
        } else {
            Ok((remaining, ()))
        }
    }
}

/// Parse a fixed-width field as a string, allow partial fields
pub(super) fn fixed_width_str_partial<'a>(
    width: usize,
) -> impl Parser<&'a [u8], Output = Option<String>, Error = NomError<&'a [u8]>> {
    map(fixed_width_partial(width, rest, true), move |opt| {
        opt.and_then(|s| Some(s.trim_ascii().to_str_lossy().into_owned()))
    })
}

/// Check if a symbol is a reserved atom symbol that requires a specific flag.
///
/// Returns `true` if the symbol is reserved and should be rejected when the corresponding
/// flag is not set. This prevents reserved symbols from being incorrectly parsed as pseudoatoms.
pub(super) fn is_reserved_atom_symbol(
    s: &[u8],
    allow_named_isotopes: bool,
    allow_wildcards: bool,
    allow_chemaxon_wildcards: bool,
    allow_electrons: bool,
    allow_rgroups: bool,
) -> bool {
    // Check for named isotopes (D, T)
    if !allow_named_isotopes && NamedIsotope::is_named_isotope_bytes(s) {
        return true;
    }

    // Check for wildcard atoms (A, Q, *, X, M) and atom lists (L)
    if !allow_wildcards {
        match s {
            b"A" | b"Q" | b"*" | b"X" | b"M" | b"L" => return true,
            _ => {}
        }
    }

    // Check for ChemAxon wildcard atoms (AH, QH, XH, MH)
    if !allow_chemaxon_wildcards {
        match s {
            b"AH" | b"QH" | b"XH" | b"MH" => return true,
            _ => {}
        }
    }

    // Check for lone pairs (LP)
    if !allow_electrons && s == b"LP" {
        return true;
    }

    // Check for R-groups (R, R#, R0, R1, R2, etc.)
    if !allow_rgroups {
        if s.starts_with(b"R") {
            // Check if it's a valid R-group pattern: "R", "R#", or "R" followed by digits
            if s.len() == 1 {
                // "R"
                return true;
            } else if s == b"R#" {
                // "R#"
                return true;
            } else if s.len() > 1 && s[1..].iter().all(|&b| b.is_ascii_digit()) {
                // "R0", "R1", "R12", etc.
                return true;
            }
        }
    }

    false
}

/// Parse position data from 3f10.4 format
pub(super) fn position30<'a>(
    ignore_positions: bool,
) -> impl Parser<&'a [u8], Output = Point3D, Error = NomError<&'a [u8]>> {
    move |input: &'a [u8]| {
        if ignore_positions {
            let (remaining, _) = take(30usize).parse(input)?;
            Ok((remaining, Point3D::zero()))
        } else {
            let x = fixed_width_float::<f64>(10, 4);
            let y = fixed_width_float::<f64>(10, 4);
            let z = fixed_width_float::<f64>(10, 4);
            map((x, y, z), |(x, y, z)| Point3D::new(x, y, z)).parse(input)
        }
    }
}

/// Parse a single RGroup occurrence.
pub(super) fn rgroup_occurrence<'a>(
) -> impl Parser<&'a [u8], Output = RGroupOccurrence, Error = NomError<&'a [u8]>> {
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
pub(super) fn rgroup_occurrences<'a>(
) -> impl Parser<&'a [u8], Output = Vec<RGroupOccurrence>, Error = NomError<&'a [u8]>> {
    delimited(
        space0,
        separated_list1(tag(","), rgroup_occurrence()),
        space0,
    )
    .or(success(vec![RGroupOccurrence::GreaterThan(0)]))
}

#[cfg(test)]
mod tests;
