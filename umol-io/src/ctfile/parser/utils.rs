//! Parsing utilities for CTfile inputs.

use std::fmt::Debug;
use std::ops::{Range, RangeInclusive};

use atoi::{FromRadix10, FromRadix10Signed};
use bstr::ByteSlice;
use fast_float2::FastFloat;
use num::{Float, Integer};
use umol_chem::element::Element;
use umol_chem::isotope::NamedIsotope;
use winnow::ascii::{alpha1, space0};
use winnow::combinator::{alt, delimited, empty, separated};
use winnow::error::{ErrMode, ParserError};
use winnow::stream::{Location, Stream};
use winnow::token::{rest, take};
use winnow::{LocatingSlice, ModalResult, Parser};

use crate::table_ir::RGroupOccurrence;

pub(super) type Input<'inp> = LocatingSlice<&'inp [u8]>;
pub(super) type PResult<T> = ModalResult<T, InputError>;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) struct InputError {
    pub(super) column: u32,
}

impl InputError {
    fn at(input: &Input<'_>) -> Self {
        Self::at_column(input.current_token_start())
    }

    pub(super) fn at_column(column: usize) -> Self {
        Self {
            column: column.min(u32::MAX as usize) as u32,
        }
    }
}

impl ParserError<Input<'_>> for InputError {
    type Inner = Self;

    fn from_input(input: &Input<'_>) -> Self {
        Self::at(input)
    }

    fn into_inner(self) -> Result<Self::Inner, Self> {
        Ok(self)
    }
}

pub(super) fn next_line<'inp>(input: &mut &'inp [u8]) -> ModalResult<Input<'inp>, ()> {
    if input.is_empty() {
        return Err(ErrMode::Backtrack(()));
    }

    let line_with_terminator_len = input
        .iter()
        .position(|byte| *byte == b'\n')
        .map_or(input.len(), |index| index + 1);
    let line_with_terminator: &[u8] = take(line_with_terminator_len).parse_next(input)?;
    let line = line_with_terminator.trim_end_with(|byte| byte == '\r' || byte == '\n');

    Ok(Input::new(line))
}

pub(super) fn finish_line(input: &mut Input<'_>) -> PResult<()> {
    let _: &[u8] = space0.parse_next(input)?;
    if input.is_empty() {
        Ok(())
    } else {
        Err(ErrMode::Backtrack(InputError::at(input)))
    }
}

pub(super) fn input_error_column(error: ErrMode<InputError>, input: &Input<'_>) -> u32 {
    error
        .into_inner()
        .map_or_else(|_| InputError::at(input).column, |error| error.column)
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
    fn parse_ascii(input: &[u8]) -> (Self, usize);

    fn parse<'inp>(input: &mut Input<'inp>) -> PResult<Self> {
        let start = input.checkpoint();
        let column = input.current_token_start();
        let bytes: &[u8] = input.as_ref();
        let (value, consumed) = Self::parse_ascii(bytes);
        if consumed == 0 {
            return Err(field_error(input, &start, column));
        }
        let _: &[u8] = take(consumed).parse_next(input)?;
        Ok(value)
    }
}

macro_rules! impl_signed_int_parser {
    ($($ty:ty),+ $(,)?) => {
        $(
            impl IntParser for $ty {
                fn parse_ascii(input: &[u8]) -> (Self, usize) {
                    Self::from_radix_10_signed(input)
                }
            }
        )+
    };
}

macro_rules! impl_unsigned_int_parser {
    ($($ty:ty),+ $(,)?) => {
        $(
            impl IntParser for $ty {
                fn parse_ascii(input: &[u8]) -> (Self, usize) {
                    Self::from_radix_10(input)
                }
            }
        )+
    };
}

impl_signed_int_parser!(i8, i16, i32);
impl_unsigned_int_parser!(u8, u32, usize);

pub(super) fn is_all_whitespace_or_zeroes(input: &[u8]) -> bool {
    input.trim_ascii().find_not_byteset(b"0").is_none()
}

fn at_field_start(mode: ErrMode<InputError>, column: usize) -> ErrMode<InputError> {
    let error = InputError::at_column(column);
    match mode {
        ErrMode::Backtrack(_) | ErrMode::Incomplete(_) => ErrMode::Backtrack(error),
        ErrMode::Cut(_) => ErrMode::Cut(error),
    }
}

fn field_error<'inp>(
    input: &mut Input<'inp>,
    checkpoint: &<Input<'inp> as Stream>::Checkpoint,
    column: usize,
) -> ErrMode<InputError> {
    input.reset(checkpoint);
    ErrMode::Backtrack(InputError::at_column(column))
}

pub(super) fn fixed_width_partial<'inp, O, P>(
    width: usize,
    mut inner: P,
    partial_ok: bool,
) -> impl Parser<Input<'inp>, Option<O>, ErrMode<InputError>>
where
    P: Parser<Input<'inp>, O, ErrMode<InputError>>,
{
    move |input: &mut Input<'inp>| {
        let start = input.checkpoint();
        let column = input.current_token_start();
        let available = width.min(input.len());
        let field: &[u8] = take(available).parse_next(input)?;

        if field.len() < width && !partial_ok && field.find_not_byteset(b"  \t").is_some() {
            return Err(field_error(input, &start, column));
        }

        if field.find_not_byteset(b"  \t").is_none() {
            return Ok(None);
        }

        let mut field_input = Input::new(field);
        match inner.parse_next(&mut field_input) {
            Ok(value) if field_input.is_empty() => Ok(Some(value)),
            Ok(_) => {
                input.reset(&start);
                Err(ErrMode::Backtrack(InputError::at_column(column)))
            }
            Err(error) => {
                input.reset(&start);
                Err(at_field_start(error, column))
            }
        }
    }
}

pub(super) fn fixed_width_opt<'inp, O, P>(
    width: usize,
    inner: P,
) -> impl Parser<Input<'inp>, Option<O>, ErrMode<InputError>>
where
    P: Parser<Input<'inp>, O, ErrMode<InputError>>,
{
    fixed_width_partial(width, inner, false)
}

fn int<'inp, T: IntParser>(input: &mut Input<'inp>) -> PResult<T> {
    T::parse(input)
}

pub(super) fn fixed_width_int<'inp, T>(
    width: usize,
) -> impl Parser<Input<'inp>, T, ErrMode<InputError>>
where
    T: IntParser,
{
    fixed_width_opt(width, delimited(space0, int::<T>, space0))
        .map(|value| value.unwrap_or_else(T::zero))
}

pub(super) fn fixed_width_int_in_range<'inp, T, R>(
    width: usize,
    range: R,
) -> impl Parser<Input<'inp>, T, ErrMode<InputError>>
where
    T: IntParser,
    R: Contains<T> + Clone,
{
    move |input: &mut Input<'inp>| {
        let start = input.checkpoint();
        let column = input.current_token_start();
        let value = fixed_width_int::<T>(width).parse_next(input)?;
        if range.contains(&value) {
            Ok(value)
        } else {
            input.reset(&start);
            Err(ErrMode::Backtrack(InputError::at_column(column)))
        }
    }
}

pub(super) fn fixed_width_int_minus1<'inp, T>(
    width: usize,
) -> impl Parser<Input<'inp>, T, ErrMode<InputError>>
where
    T: IntParser,
{
    move |input: &mut Input<'inp>| {
        let start = input.checkpoint();
        let column = input.current_token_start();
        let value = fixed_width_int::<T>(width).parse_next(input)?;
        if value >= T::one() {
            Ok(value - T::one())
        } else {
            input.reset(&start);
            Err(ErrMode::Backtrack(InputError::at_column(column)))
        }
    }
}

#[inline(always)]
pub(super) fn fixed_width_float_f10_4<'inp, T>() -> impl Parser<Input<'inp>, T, ErrMode<InputError>>
where
    T: Float + FastFloat,
{
    move |input: &mut Input<'inp>| {
        let start = input.checkpoint();
        let column = input.current_token_start();
        let available = 10.min(input.len());
        let field: &[u8] = take(available).parse_next(input)?;
        let trimmed = field.trim_ascii();
        if trimmed.is_empty() {
            return Ok(T::zero());
        }
        if trimmed.find_not_byteset(b"0123456789+-.").is_some() {
            input.reset(&start);
            return Err(ErrMode::Backtrack(InputError::at_column(column)));
        }
        let value = match fast_float2::parse::<T, _>(trimmed) {
            Ok(value) => value,
            Err(_) => {
                input.reset(&start);
                return Err(ErrMode::Backtrack(InputError::at_column(column)));
            }
        };
        if trimmed.find_byte(b'.').is_some() {
            Ok(value)
        } else {
            Ok(value / T::from(10.0).unwrap().powi(4))
        }
    }
}

pub(super) fn fixed_width_element_partial<'inp>(
    width: usize,
) -> impl Parser<Input<'inp>, Option<Element>, ErrMode<InputError>> {
    fixed_width_partial(
        width,
        delimited(
            space0,
            alpha1.verify_map(Element::from_symbol_bytes),
            space0,
        ),
        true,
    )
}

pub(super) fn fixed_width_unused<'inp>(
    width: usize,
    skip_unused_fields: bool,
) -> impl Parser<Input<'inp>, (), ErrMode<InputError>> {
    move |input: &mut Input<'inp>| {
        let start = input.checkpoint();
        let column = input.current_token_start();
        let unused: &[u8] = take(width).parse_next(input)?;
        if skip_unused_fields || is_all_whitespace_or_zeroes(unused) {
            Ok(())
        } else {
            input.reset(&start);
            Err(ErrMode::Backtrack(InputError::at_column(column)))
        }
    }
}

fn string_field<'inp>(input: &mut Input<'inp>) -> PResult<String> {
    let value: &[u8] = rest.parse_next(input)?;
    Ok(value.trim_ascii().to_str_lossy().into_owned())
}

pub(super) fn fixed_width_str_partial<'inp>(
    width: usize,
) -> impl Parser<Input<'inp>, Option<String>, ErrMode<InputError>> {
    fixed_width_partial(width, string_field, true)
}

fn u8_value<'inp>(input: &mut Input<'inp>) -> PResult<u8> {
    <u8 as IntParser>::parse(input)
}

pub(super) fn rgroup_occurrence<'inp>(
) -> impl Parser<Input<'inp>, RGroupOccurrence, ErrMode<InputError>> {
    alt((
        (u8_value, b'-', u8_value).map(|(lower, _, upper)| RGroupOccurrence::Range(lower, upper)),
        u8_value.map(RGroupOccurrence::Exactly),
        (b'>', u8_value).map(|(_, value)| RGroupOccurrence::GreaterThan(value)),
        (b'<', u8_value).map(|(_, value)| RGroupOccurrence::FewerThan(value)),
    ))
}

pub(super) fn rgroup_occurrences<'inp>(
) -> impl Parser<Input<'inp>, Vec<RGroupOccurrence>, ErrMode<InputError>> {
    alt((
        delimited(space0, separated(1.., rgroup_occurrence(), b','), space0),
        empty.value(vec![RGroupOccurrence::GreaterThan(0)]),
    ))
}

pub(super) fn is_reserved_atom_symbol(
    symbol: &[u8],
    allow_named_isotopes: bool,
    allow_wildcards: bool,
    allow_chemaxon_wildcards: bool,
    allow_electrons: bool,
    allow_rgroups: bool,
) -> bool {
    if !allow_named_isotopes && NamedIsotope::is_named_isotope_bytes(symbol) {
        return true;
    }

    if !allow_wildcards && matches!(symbol, b"A" | b"Q" | b"*" | b"X" | b"M" | b"L") {
        return true;
    }

    if !allow_chemaxon_wildcards && matches!(symbol, b"AH" | b"QH" | b"XH" | b"MH") {
        return true;
    }

    if !allow_electrons && symbol == b"LP" {
        return true;
    }

    if !allow_rgroups && symbol.starts_with(b"R") {
        return symbol.len() == 1
            || symbol == b"R#"
            || symbol[1..].iter().all(|byte| byte.is_ascii_digit());
    }

    false
}

#[inline(always)]
pub(super) fn parse_int_opt<T: IntParser>(field: &[u8], column: u32) -> PResult<Option<T>> {
    let trimmed = field.trim_ascii();
    if trimmed.is_empty() {
        return Ok(None);
    }

    let mut input = Input::new(trimmed);
    match T::parse(&mut input) {
        Ok(value) if input.is_empty() => Ok(Some(value)),
        Ok(_) | Err(_) => Err(ErrMode::Backtrack(InputError { column })),
    }
}

#[inline(always)]
pub(super) fn parse_float_f10_4(field: &[u8], column: u32) -> PResult<f64> {
    let trimmed = field.trim_ascii();
    if trimmed.is_empty() {
        return Ok(0.0);
    }
    if trimmed.find_not_byteset(b"0123456789+-.").is_some() {
        return Err(ErrMode::Backtrack(InputError { column }));
    }

    let value = fast_float2::parse::<f64, _>(trimmed)
        .map_err(|_| ErrMode::Backtrack(InputError { column }))?;
    if trimmed.find_byte(b'.').is_some() {
        Ok(value)
    } else {
        Ok(value / 10_000.0)
    }
}

#[inline(always)]
pub(super) fn validate_unused_n(
    field: &[u8],
    count: usize,
    width: usize,
    skip_unused_fields: bool,
    column: u32,
) -> PResult<()> {
    let required = count.saturating_mul(width);
    if field.len() < required {
        return Err(ErrMode::Backtrack(InputError { column }));
    }
    if skip_unused_fields || count == 0 || width == 0 {
        return Ok(());
    }
    if field[..required]
        .chunks_exact(width)
        .all(is_all_whitespace_or_zeroes)
    {
        Ok(())
    } else {
        Err(ErrMode::Backtrack(InputError { column }))
    }
}

#[cfg(test)]
mod tests;
