//! Parsing utilities for CTab files.

use std::fmt::Debug;
use std::ops::{Range, RangeInclusive};

use bstr::ByteSlice;
use fast_float2::FastFloat;
use nom::branch::alt;
use nom::bytes::complete::{tag, take};
use nom::character::complete::{
    alpha1, i16 as nom_i16, i32 as nom_i32, i8 as nom_i8, space0, u32 as nom_u32, u8 as nom_u8,
    usize as nom_usize,
};
use nom::combinator::{map, map_opt, rest, success, verify};
use nom::error::{Error as NomError, ErrorKind as NomErrorKind};
use nom::multi::separated_list1;
use nom::sequence::delimited;
use nom::{Err, Parser};
use num::{Float, Integer};
use umol_data::{Element, NamedIsotope};

use crate::table_ir::RGroupOccurrence;

/// Iterator over lines that yields each line without terminator and its byte length including terminator.
pub(super) struct LinesWithOffset<'inp> {
    inner: bstr::LinesWithTerminator<'inp>,
}

impl<'inp> LinesWithOffset<'inp> {
    pub(super) fn new(input: &'inp [u8]) -> Self {
        Self {
            inner: input.lines_with_terminator(),
        }
    }
}

impl<'inp> Iterator for LinesWithOffset<'inp> {
    type Item = (&'inp [u8], usize);

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|line_with_term| {
            let byte_len = line_with_term.len();
            let line = line_with_term.trim_end_with(|c| c == '\r' || c == '\n');
            (line, byte_len)
        })
    }
}

pub(super) trait LinesWithOffsetExt<'inp> {
    fn lines_with_offset(self) -> LinesWithOffset<'inp>;
}

impl<'inp> LinesWithOffsetExt<'inp> for &'inp [u8] {
    fn lines_with_offset(self) -> LinesWithOffset<'inp> {
        LinesWithOffset::new(self)
    }
}

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
    fn nom_parser<'inp>() -> impl Parser<&'inp [u8], Output = Self, Error = NomError<&'inp [u8]>>;
}

impl IntParser for i8 {
    fn nom_parser<'inp>() -> impl Parser<&'inp [u8], Output = Self, Error = NomError<&'inp [u8]>> {
        nom_i8
    }
}

impl IntParser for i16 {
    fn nom_parser<'inp>() -> impl Parser<&'inp [u8], Output = Self, Error = NomError<&'inp [u8]>> {
        nom_i16
    }
}

impl IntParser for i32 {
    fn nom_parser<'inp>() -> impl Parser<&'inp [u8], Output = Self, Error = NomError<&'inp [u8]>> {
        nom_i32
    }
}

impl IntParser for u8 {
    fn nom_parser<'inp>() -> impl Parser<&'inp [u8], Output = Self, Error = NomError<&'inp [u8]>> {
        nom_u8
    }
}

impl IntParser for u32 {
    fn nom_parser<'inp>() -> impl Parser<&'inp [u8], Output = Self, Error = NomError<&'inp [u8]>> {
        nom_u32
    }
}

impl IntParser for usize {
    fn nom_parser<'inp>() -> impl Parser<&'inp [u8], Output = Self, Error = NomError<&'inp [u8]>> {
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
pub(super) fn fixed_width_partial<'inp, O, P>(
    width: usize,
    mut inner: P,
    partial_ok: bool,
) -> impl Parser<&'inp [u8], Output = Option<O>, Error = NomError<&'inp [u8]>>
where
    P: Parser<&'inp [u8], Output = O, Error = NomError<&'inp [u8]>>,
{
    move |input: &'inp [u8]| {
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
pub(super) fn fixed_width_opt<'inp, O, P>(
    width: usize,
    inner: P,
) -> impl Parser<&'inp [u8], Output = Option<O>, Error = NomError<&'inp [u8]>>
where
    P: Parser<&'inp [u8], Output = O, Error = NomError<&'inp [u8]>>,
{
    fixed_width_partial(width, inner, false)
}

/// Parse a fixed-width field as an integer type. Interprets empty/whitespace field as default.
pub(super) fn fixed_width_int<'inp, T>(
    width: usize,
) -> impl Parser<&'inp [u8], Output = T, Error = NomError<&'inp [u8]>>
where
    T: IntParser,
{
    map(
        fixed_width_opt(width, delimited(space0, T::nom_parser(), space0)),
        |opt| opt.unwrap_or_else(T::zero),
    )
}

/// Parse a fixed-width field as an integer type, applying range bounds.
pub(super) fn fixed_width_int_in_range<'inp, T, R>(
    width: usize,
    range: R,
) -> impl Parser<&'inp [u8], Output = T, Error = NomError<&'inp [u8]>>
where
    T: IntParser,
    R: Contains<T> + Clone,
{
    verify(fixed_width_int::<T>(width), move |val: &T| {
        range.contains(val)
    })
}

/// Parse a fixed-width field as an integer type, subtracting one.
pub(super) fn fixed_width_int_minus1<'inp, T>(
    width: usize,
) -> impl Parser<&'inp [u8], Output = T, Error = NomError<&'inp [u8]>>
where
    T: IntParser,
{
    map(
        verify(fixed_width_int(width), |val: &T| *val >= T::one()),
        |x: T| x - T::one(),
    )
}

/// Parse a fixed-width field as integer, allow partial fields
pub(super) fn fixed_width_int_partial<'inp, T>(
    width: usize,
) -> impl Parser<&'inp [u8], Output = T, Error = NomError<&'inp [u8]>>
where
    T: IntParser,
{
    map(
        fixed_width_partial(width, delimited(space0, T::nom_parser(), space0), true),
        |opt| opt.unwrap_or_else(T::zero),
    )
}

/// Parse a fixed-width field as float with Fortran semantics (Fw.d).
/// Parser combinator version for use with nom combinators.
#[inline(always)]
pub(super) fn fixed_width_float_f10_4<'inp, T>(
) -> impl Parser<&'inp [u8], Output = T, Error = NomError<&'inp [u8]>>
where
    T: Float + FastFloat,
{
    move |input: &'inp [u8]| {
        let min_width = 10.min(input.len());
        let (remaining, field) = take(min_width).parse(input)?;
        let trimmed = field.trim_ascii();
        if trimmed.is_empty() {
            return Ok((remaining, T::zero()));
        }
        if trimmed.find_not_byteset(b"0123456789+-.").is_some() {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Digit)));
        }
        let val = if trimmed.find_byte(b'.').is_some() {
            fast_float2::parse::<T, _>(trimmed)
                .map_err(|_| Err::Error(NomError::new(input, NomErrorKind::Digit)))?
        } else {
            let val = fast_float2::parse::<T, _>(trimmed)
                .map_err(|_| Err::Error(NomError::new(input, NomErrorKind::Digit)))?;
            val / T::from(10.0).unwrap().powi(4)
        };
        Ok((remaining, val))
    }
}

/// Parse a fixed-width field as element symbol
pub(super) fn fixed_width_element_partial<'inp>(
    width: usize,
) -> impl Parser<&'inp [u8], Output = Option<Element>, Error = NomError<&'inp [u8]>> {
    fixed_width_partial(
        width,
        delimited(space0, map_opt(alpha1, Element::from_symbol_bytes), space0),
        true,
    )
}

/// Unused field of fixed width `width`
/// Only validate unused field if `skip_unused_fields` is false.
pub(super) fn fixed_width_unused<'inp>(
    width: usize,
    skip_unused_fields: bool,
) -> impl Parser<&'inp [u8], Output = (), Error = NomError<&'inp [u8]>> {
    move |input: &'inp [u8]| {
        let (remaining, unused) = take(width).parse(input)?;
        if skip_unused_fields {
            return Ok((remaining, ()));
        }
        if !is_all_whitespace_or_zeroes(&unused) {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Verify)));
        }
        Ok((remaining, ()))
    }
}

/// Parse a fixed-width field as a string, allow partial fields
pub(super) fn fixed_width_str_partial<'inp>(
    width: usize,
) -> impl Parser<&'inp [u8], Output = Option<String>, Error = NomError<&'inp [u8]>> {
    map(fixed_width_partial(width, rest, true), move |opt| {
        opt.and_then(|s| Some(s.trim_ascii().to_str_lossy().into_owned()))
    })
}

/// Parse a single RGroup occurrence.
pub(super) fn rgroup_occurrence<'inp>(
) -> impl Parser<&'inp [u8], Output = RGroupOccurrence, Error = NomError<&'inp [u8]>> {
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
pub(super) fn rgroup_occurrences<'inp>(
) -> impl Parser<&'inp [u8], Output = Vec<RGroupOccurrence>, Error = NomError<&'inp [u8]>> {
    delimited(
        space0,
        separated_list1(tag(","), rgroup_occurrence()),
        space0,
    )
    .or(success(vec![RGroupOccurrence::GreaterThan(0)]))
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

/// Parse an optional integer field. Returns None if field is whitespace only.
#[inline(always)]
pub(super) fn parse_int_opt<'inp, T: IntParser>(
    input: &'inp [u8],
    field: &[u8],
) -> Result<Option<T>, Err<NomError<&'inp [u8]>>> {
    let trimmed = field.trim_ascii();
    if trimmed.is_empty() {
        return Ok(None);
    }
    match T::nom_parser().parse(trimmed) {
        Ok((remaining, val)) if remaining.is_empty() => Ok(Some(val)),
        Ok(_) => Err(Err::Error(NomError::new(input, NomErrorKind::Eof))), // trailing garbage
        Err(Err::Error(_)) => Err(Err::Error(NomError::new(input, NomErrorKind::Digit))),
        Err(Err::Failure(e)) => Err(Err::Failure(NomError::new(input, e.code))),
        Err(Err::Incomplete(n)) => Err(Err::Incomplete(n)),
    }
}

/// Parse a float field with F10.4 format. Returns 0.0 if field is whitespace only.
#[inline(always)]
pub(super) fn parse_float_f10_4<'inp>(
    input: &'inp [u8],
    field: &[u8],
) -> Result<f64, Err<NomError<&'inp [u8]>>> {
    let trimmed = field.trim_ascii();
    if trimmed.is_empty() {
        return Ok(0.0);
    }
    if trimmed.find_not_byteset(b"0123456789+-.").is_some() {
        return Err(Err::Error(NomError::new(input, NomErrorKind::Digit)));
    }
    if trimmed.find_byte(b'.').is_some() {
        fast_float2::parse::<f64, _>(trimmed)
            .map_err(|_| Err::Error(NomError::new(input, NomErrorKind::Digit)))
    } else {
        let val = fast_float2::parse::<f64, _>(trimmed)
            .map_err(|_| Err::Error(NomError::new(input, NomErrorKind::Digit)))?;
        Ok(val / 10_000.0)
    }
}

/// Validate that `count` consecutive fields of `width` bytes each contain only whitespace or zeros.
/// If `skip_unused_fields` is true, validation is skipped and always succeeds.
#[inline(always)]
pub(super) fn validate_unused_n<'inp>(
    input: &'inp [u8],
    field: &[u8],
    count: usize,
    width: usize,
    skip_unused_fields: bool,
) -> Result<(), Err<NomError<&'inp [u8]>>> {
    if field.len() < count * width {
        return Err(Err::Error(NomError::new(input, NomErrorKind::Eof)));
    }
    if skip_unused_fields || count == 0 || width == 0 {
        return Ok(());
    }
    for i in 0..count {
        let start = i * width;
        let end = start + width;
        if end > field.len() {
            break;
        }
        if !is_all_whitespace_or_zeroes(&field[start..end]) {
            return Err(Err::Error(NomError::new(input, NomErrorKind::Verify)));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests;
