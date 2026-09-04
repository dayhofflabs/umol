use float_cmp::approx_eq;
use pretty_assertions::assert_eq;
use rstest::rstest;
use winnow::ascii::space0;
use winnow::combinator::delimited;
use winnow::error::ErrMode;
use winnow::stream::Location;
use winnow::token::take;
use winnow::Parser;

use super::*;

#[rstest]
#[case::single_line_lf(b"abc\n", b"abc", b"")]
#[case::single_line_crlf(b"abc\r\n", b"abc", b"")]
#[case::single_line_no_term(b"abc", b"abc", b"")]
#[case::first_of_two_lf(b"abc\ndef\n", b"abc", b"def\n")]
#[case::first_of_two_crlf(b"abc\r\ndef\r\n", b"abc", b"def\r\n")]
#[case::interior_cr(b"abc\rdef\n", b"abc\rdef", b"")]
#[case::trailing_cr(b"abc\r", b"abc", b"")]
#[case::empty_line(b"\nnext", b"", b"next")]
fn test_next_line(
    #[case] input: &[u8],
    #[case] expected_line: &[u8],
    #[case] expected_remaining: &[u8],
) {
    let mut remaining = input;
    let line = next_line(&mut remaining).unwrap();
    assert_eq!(*line.as_ref(), expected_line);
    assert_eq!(remaining, expected_remaining);
}

#[rstest]
#[case::empty(b"")]
fn test_next_line_error(#[case] input: &[u8]) {
    let mut remaining = input;
    assert_eq!(next_line(&mut remaining), Err(ErrMode::Backtrack(())));
}

#[rstest]
#[case::empty(b"")]
#[case::spaces(b"   ")]
#[case::tabs(b"\t\t")]
fn test_finish_line(#[case] input: &[u8]) {
    let mut input = Input::new(input);
    assert_eq!(finish_line(&mut input), Ok(()));
}

#[rstest]
#[case::text(b"x", 0)]
#[case::after_spaces(b"  x", 2)]
fn test_finish_line_error(#[case] input: &[u8], #[case] column: u32) {
    let mut input = Input::new(input);
    assert_eq!(
        finish_line(&mut input),
        Err(ErrMode::Backtrack(InputError { column }))
    );
}

#[rstest]
#[case::empty(b"", true)]
#[case::whitespace(b"   ", true)]
#[case::zero(b"  0", true)]
#[case::zero_pos2(b" 0 ", true)]
#[case::zero_pos1(b"0  ", true)]
#[case::two_zeros(b" 00", true)]
#[case::two_zeros_pos1(b"00 ", true)]
#[case::three_zeros(b"000", true)]
#[case::zero_width1(b"0", true)]
#[case::two_zeros_width1(b"00", true)]
#[case::two_zeros_separated(b"0 0", false)]
#[case::one(b"  1", false)]
fn test_is_all_whitespace_or_zeroes(#[case] input: &[u8], #[case] expected: bool) {
    assert_eq!(is_all_whitespace_or_zeroes(input), expected);
}

#[rstest]
#[case::empty(b"", None)]
#[case::blank_width2(b"  ", None)]
#[case::blank_width3(b"   ", None)]
#[case::number(b"12", Some(12))]
#[case::number_pos2(b" 12", Some(12))]
#[case::number_pos1(b"12 ", Some(12))]
#[case::number_zero_padded(b"012", Some(12))]
fn test_fixed_width_partial(#[case] input: &[u8], #[case] expected: Option<i32>) {
    let mut parser = fixed_width_partial(3, delimited(space0, int::<i32>, space0), true);
    assert_eq!(parser.parse(Input::new(input)), Ok(expected));
}

#[rstest]
#[case::non_numeric(b"abc", 0)]
#[case::trailing_characters(b"1a ", 0)]
fn test_fixed_width_partial_error(#[case] input: &[u8], #[case] column: u32) {
    let mut parser = fixed_width_partial(3, delimited(space0, int::<i32>, space0), true);
    let error = parser.parse(Input::new(input)).unwrap_err().into_inner();
    assert_eq!(error, InputError { column });
}

#[rstest]
#[case::empty(b"", None)]
#[case::blank_width2(b"  ", None)]
#[case::blank_width3(b"   ", None)]
#[case::number(b" 12", Some(12))]
fn test_fixed_width_opt(#[case] input: &[u8], #[case] expected: Option<i32>) {
    let mut parser = fixed_width_opt(3, delimited(space0, int::<i32>, space0));
    assert_eq!(parser.parse(Input::new(input)), Ok(expected));
}

#[rstest]
#[case::non_numeric(b" abc ", 0)]
#[case::two_characters(b" 1", 0)]
#[case::one_character(b"1", 0)]
fn test_fixed_width_opt_error(#[case] input: &[u8], #[case] column: u32) {
    let mut parser = fixed_width_opt(5, delimited(space0, int::<i32>, space0));
    let error = parser.parse(Input::new(input)).unwrap_err().into_inner();
    assert_eq!(error, InputError { column });
}

#[rstest]
#[case::positive(b"123", 123)]
#[case::negative(b"-98", -98)]
#[case::padded(b"  8", 8)]
#[case::blank(b"   ", 0)]
#[case::blank_width2(b"  ", 0)]
fn test_fixed_width_int(#[case] input: &[u8], #[case] expected: i32) {
    let mut parser = fixed_width_int::<i32>(3);
    assert_eq!(parser.parse(Input::new(input)), Ok(expected));
}

#[rstest]
#[case::four_characters(b"1234", 3)]
#[case::two_characters(b"12", 0)]
#[case::non_numeric(b"abc", 0)]
#[case::trailing_characters(b"1a ", 0)]
fn test_fixed_width_int_error(#[case] input: &[u8], #[case] column: u32) {
    let mut parser = fixed_width_int::<i32>(3);
    let error = parser.parse(Input::new(input)).unwrap_err().into_inner();
    assert_eq!(error, InputError { column });
}

#[rstest]
#[case::after_prefix(b"xxabc", 2, 2)]
fn test_fixed_width_int_location_error(
    #[case] input: &[u8],
    #[case] prefix: usize,
    #[case] column: u32,
) {
    let mut input = Input::new(input);
    let parsed: PResult<&[u8]> = take(prefix).parse_next(&mut input);
    parsed.unwrap();
    let error = fixed_width_int::<i32>(3)
        .parse_next(&mut input)
        .unwrap_err();
    assert_eq!(error, ErrMode::Backtrack(InputError { column }));
    assert_eq!(input.current_token_start(), prefix);
}

#[rstest]
#[case::positive(b"100", 100)]
#[case::negative(b" -9", -9)]
#[case::padded_right(b"8  ", 8)]
#[case::padded_both_sides(b" 1 ", 1)]
fn test_fixed_width_int_in_range(#[case] input: &[u8], #[case] expected: i8) {
    let mut parser = fixed_width_int_in_range::<i8, _>(3, -10..=110);
    assert_eq!(parser.parse(Input::new(input)), Ok(expected));
}

#[rstest]
#[case::blank(b"   ", 0)]
#[case::out_of_range(b"11 ", 0)]
#[case::four_characters(b"1234", 0)]
#[case::one_character(b"8", 0)]
#[case::non_numeric(b"abc", 0)]
#[case::trailing_characters(b"1a ", 0)]
fn test_fixed_width_int_in_range_error(#[case] input: &[u8], #[case] column: u32) {
    let mut parser = fixed_width_int_in_range::<i8, _>(3, 1..=10);
    let error = parser.parse(Input::new(input)).unwrap_err().into_inner();
    assert_eq!(error, InputError { column });
}

#[rstest]
#[case::positive(b"100", 100)]
#[case::padded_left(b"  9", 9)]
#[case::padded_right(b"8  ", 8)]
#[case::padded_both_sides(b" 1 ", 1)]
fn test_fixed_width_int_in_range_inclusive(#[case] input: &[u8], #[case] expected: u8) {
    let mut parser = fixed_width_int_in_range::<u8, _>(3, 0..=100);
    assert_eq!(parser.parse(Input::new(input)), Ok(expected));
}

#[rstest]
#[case::padded_left(b"  1", 0)]
#[case::three_digits(b"123", 122)]
fn test_fixed_width_int_minus1(#[case] input: &[u8], #[case] expected: usize) {
    let mut parser = fixed_width_int_minus1::<usize>(3);
    assert_eq!(parser.parse(Input::new(input)), Ok(expected));
}

#[rstest]
#[case::zero(b"  0", 0)]
fn test_fixed_width_int_minus1_error(#[case] input: &[u8], #[case] column: u32) {
    let mut parser = fixed_width_int_minus1::<usize>(3);
    let error = parser.parse(Input::new(input)).unwrap_err().into_inner();
    assert_eq!(error, InputError { column });
}

#[rstest]
#[case::three_digits(b"123", 123)]
#[case::two_digits(b"12", 12)]
#[case::one_digit(b"1", 1)]
#[case::empty(b"", 0)]
#[case::two_digits_padded_right(b"12 ", 12)]
#[case::one_digit_padded_right(b"1  ", 1)]
#[case::two_digits_padded_left(b" 12", 12)]
#[case::one_digit_padded_left(b"  1", 1)]
#[case::two_digits_padded_both_sides(b" 1 ", 1)]
#[case::blank_width3(b"   ", 0)]
#[case::blank_width2(b"  ", 0)]
#[case::blank_width1(b" ", 0)]
#[case::negative(b" -1", -1)]
fn test_fixed_width_int_partial(#[case] input: &[u8], #[case] expected: i32) {
    let mut parser = fixed_width_int_partial::<i32>(3);
    assert_eq!(parser.parse(Input::new(input)), Ok(expected));
}

#[rstest]
#[case::too_many_characters(b"1234", 3)]
#[case::non_numeric(b"abc", 0)]
#[case::trailing_characters(b"1a ", 0)]
fn test_fixed_width_int_partial_error(#[case] input: &[u8], #[case] column: u32) {
    let mut parser = fixed_width_int_partial::<i32>(3);
    let error = parser.parse(Input::new(input)).unwrap_err().into_inner();
    assert_eq!(error, InputError { column });
}

#[rstest]
#[case::padded_both_sides(b"  1.2345  ", 1.2345)]
#[case::negative(b"    -1.234", -1.234)]
#[case::dot_zero_padded_right(b"1.0       ", 1.0)]
#[case::no_fractional_part_padded_right(b"1.        ", 1.0)]
#[case::no_integer_part_padded_both_sides(b" .1       ", 0.1)]
#[case::padded_right(b"1.23456   ", 1.23456)]
#[case::integer_padded_left(b"   1234567", 123.4567)]
#[case::negative_padded_left(b"  -1234567", -123.4567)]
#[case::blank(b"          ", 0.0)]
#[case::blank_width9(b"         ", 0.0)]
fn test_fixed_width_float_f10_4(#[case] input: &[u8], #[case] expected: f64) {
    let mut parser = fixed_width_float_f10_4::<f64>();
    let value = parser.parse(Input::new(input)).unwrap();
    assert!(approx_eq!(f64, value, expected, ulps = 4));
}

#[rstest]
#[case::trailing_characters(b"1.23a     ", 0)]
#[case::invalid_decimal_point(b"1.2.3     ", 0)]
#[case::trailing_data(b"          a", 10)]
fn test_fixed_width_float_f10_4_error(#[case] input: &[u8], #[case] column: u32) {
    let mut parser = fixed_width_float_f10_4::<f64>();
    let error = parser.parse(Input::new(input)).unwrap_err().into_inner();
    assert_eq!(error, InputError { column });
}

#[rstest]
#[case::one_character(b"   C", Some(Element::C))]
#[case::one_character_pos3(b"  C ", Some(Element::C))]
#[case::one_character_pos2(b" C  ", Some(Element::C))]
#[case::one_character_pos1(b"C   ", Some(Element::C))]
#[case::two_characters_pos1(b"Cu  ", Some(Element::Cu))]
#[case::two_characters_pos2(b" Cu ", Some(Element::Cu))]
#[case::two_characters_pos3(b"  Cu", Some(Element::Cu))]
#[case::blank(b"   ", None)]
#[case::blank_width2(b"  ", None)]
fn test_fixed_width_element_partial(#[case] input: &[u8], #[case] expected: Option<Element>) {
    let mut parser = fixed_width_element_partial(4);
    assert_eq!(parser.parse(Input::new(input)), Ok(expected));
}

#[rstest]
#[case::trailing_characters(b"Cu   ", 4)]
#[case::invalid_element_symbol(b" X  ", 0)]
fn test_fixed_width_element_partial_error(#[case] input: &[u8], #[case] column: u32) {
    let mut parser = fixed_width_element_partial(4);
    let error = parser.parse(Input::new(input)).unwrap_err().into_inner();
    assert_eq!(error, InputError { column });
}

#[rstest]
#[case::blank(b"   ", 3, false)]
#[case::zero(b"  0", 3, false)]
#[case::zero_pos2(b" 0 ", 3, false)]
#[case::zero_pos1(b"0  ", 3, false)]
#[case::two_zeros(b" 00", 3, false)]
#[case::two_zeros_pos1(b"00 ", 3, false)]
#[case::three_zeros(b"000", 3, false)]
#[case::three_zeros_skipped(b"000", 3, true)]
#[case::nonzero_skipped(b"  1", 3, true)]
fn test_fixed_width_unused(
    #[case] input: &[u8],
    #[case] width: usize,
    #[case] skip_unused_fields: bool,
) {
    let mut parser = fixed_width_unused(width, skip_unused_fields);
    assert_eq!(parser.parse(Input::new(input)), Ok(()));
}

#[rstest]
#[case::empty(b"", 3, false, 0)]
#[case::width1(b"0", 3, false, 0)]
#[case::width2(b"00", 3, false, 0)]
#[case::zeros_separated(b"0 0", 3, false, 0)]
#[case::nonzero(b"  1", 3, false, 0)]
#[case::non_numeric(b" a ", 3, false, 0)]
#[case::short_skipped(b"  ", 3, true, 0)]
fn test_fixed_width_unused_error(
    #[case] input: &[u8],
    #[case] width: usize,
    #[case] skip_unused_fields: bool,
    #[case] column: u32,
) {
    let mut parser = fixed_width_unused(width, skip_unused_fields);
    let error = parser.parse(Input::new(input)).unwrap_err().into_inner();
    assert_eq!(error, InputError { column });
}

#[rstest]
#[case::full_field(b"abcd", 4, Some("abcd".to_string()))]
#[case::padded_right(b"abc ", 4, Some("abc".to_string()))]
#[case::too_short(b"abc", 4, Some("abc".to_string()))]
#[case::padded_left(b" abc", 4, Some("abc".to_string()))]
#[case::padded_both_sides(b" ab ", 4, Some("ab".to_string()))]
#[case::empty(b"", 4, None)]
#[case::blank(b"   ", 4, None)]
fn test_fixed_width_str_partial(
    #[case] input: &[u8],
    #[case] width: usize,
    #[case] expected: Option<String>,
) {
    let mut parser = fixed_width_str_partial(width);
    assert_eq!(parser.parse(Input::new(input)), Ok(expected));
}

#[rstest]
#[case::empty(b"", vec![RGroupOccurrence::GreaterThan(0)])]
#[case::one(b"1", vec![RGroupOccurrence::Exactly(1)])]
#[case::zero_padded(b"01", vec![RGroupOccurrence::Exactly(1)])]
#[case::two(b"1,2", vec![RGroupOccurrence::Exactly(1), RGroupOccurrence::Exactly(2)])]
#[case::greater_than(b">1", vec![RGroupOccurrence::GreaterThan(1)])]
#[case::fewer_than(b"<2", vec![RGroupOccurrence::FewerThan(2)])]
#[case::range(b"1-3", vec![RGroupOccurrence::Range(1, 3)])]
#[case::mixed(b"0,>0", vec![RGroupOccurrence::Exactly(0), RGroupOccurrence::GreaterThan(0)])]
fn test_rgroup_occurrences(#[case] input: &[u8], #[case] expected: Vec<RGroupOccurrence>) {
    let mut parser = rgroup_occurrences();
    assert_eq!(parser.parse(Input::new(input)), Ok(expected));
}

#[rstest]
#[case::invalid_character(b"a", 0)]
#[case::negative_value(b"-3", 0)]
fn test_rgroup_occurrences_error(#[case] input: &[u8], #[case] column: u32) {
    let mut parser = rgroup_occurrences();
    let error = parser.parse(Input::new(input)).unwrap_err().into_inner();
    assert_eq!(error, InputError { column });
}

#[rstest]
#[case::named_isotope_disallowed(b"D", false, false, false, false, false, true)]
#[case::named_isotope_allowed(b"D", true, false, false, false, false, false)]
#[case::wildcard_disallowed(b"A", false, false, false, false, false, true)]
#[case::wildcard_allowed(b"A", false, true, false, false, false, false)]
#[case::electrons_disallowed(b"LP", false, false, false, false, false, true)]
#[case::electrons_allowed(b"LP", false, false, false, true, false, false)]
#[case::rgroup_disallowed(b"R1", false, false, false, false, false, true)]
#[case::rgroup_allowed(b"R1", false, false, false, false, true, false)]
#[case::element_symbol(b"C", false, false, false, false, false, false)]
fn test_is_reserved_atom_symbol(
    #[case] symbol: &[u8],
    #[case] allow_named_isotopes: bool,
    #[case] allow_wildcards: bool,
    #[case] allow_chemaxon_wildcards: bool,
    #[case] allow_electrons: bool,
    #[case] allow_rgroups: bool,
    #[case] expected: bool,
) {
    assert_eq!(
        is_reserved_atom_symbol(
            symbol,
            allow_named_isotopes,
            allow_wildcards,
            allow_chemaxon_wildcards,
            allow_electrons,
            allow_rgroups,
        ),
        expected
    );
}

#[rstest]
#[case::whitespace(b"   ", None::<i32>)]
#[case::zero(b"  0", Some(0))]
#[case::positive(b" 17", Some(17))]
#[case::negative(b" -5", Some(-5))]
#[case::padded(b"123", Some(123))]
fn test_parse_int_opt(#[case] field: &[u8], #[case] expected: Option<i32>) {
    assert_eq!(parse_int_opt::<i32>(field, 7), Ok(expected));
}

#[rstest]
#[case::whitespace(b"   ", None::<u8>)]
#[case::zero(b"  0", Some(0))]
#[case::positive(b" 17", Some(17))]
#[case::maximum(b"255", Some(255))]
fn test_parse_int_opt_unsigned(#[case] field: &[u8], #[case] expected: Option<u8>) {
    assert_eq!(parse_int_opt::<u8>(field, 7), Ok(expected));
}

#[rstest]
#[case::non_digit(b"abc", 7)]
#[case::mixed(b"12a", 7)]
#[case::interstitial_space(b"12 3", 7)]
fn test_parse_int_opt_error(#[case] field: &[u8], #[case] column: u32) {
    assert_eq!(
        parse_int_opt::<i32>(field, column),
        Err(ErrMode::Backtrack(InputError { column }))
    );
}

#[rstest]
#[case::padded_both_sides(b"  1.2345  ", 1.2345)]
#[case::negative(b"    -1.234", -1.234)]
#[case::dot_zero_padded_right(b"1.0       ", 1.0)]
#[case::no_fractional_part_padded_right(b"1.        ", 1.0)]
#[case::no_integer_part_padded_both_sides(b" .1       ", 0.1)]
#[case::padded_right(b"1.23456   ", 1.23456)]
#[case::integer_padded_left(b"   1234567", 123.4567)]
#[case::negative_padded_left(b"  -1234567", -123.4567)]
#[case::blank(b"          ", 0.0)]
#[case::blank_width9(b"         ", 0.0)]
fn test_parse_float_f10_4(#[case] field: &[u8], #[case] expected: f64) {
    assert_eq!(parse_float_f10_4(field, 6), Ok(expected));
}

#[rstest]
#[case::trailing_characters(b"1.23a     ", 6)]
#[case::invalid_decimal_point(b"1.2.3     ", 6)]
#[case::trailing_data(b"          a", 6)]
fn test_parse_float_f10_4_error(#[case] field: &[u8], #[case] column: u32) {
    assert_eq!(
        parse_float_f10_4(field, column),
        Err(ErrMode::Backtrack(InputError { column }))
    );
}

#[rstest]
#[case::empty(b"", 0, 0, false)]
#[case::blank(b"   ", 1, 3, false)]
#[case::zero(b"  0", 1, 3, false)]
#[case::two_fields(b"000000", 2, 3, false)]
#[case::blank_then_zero(b"   000", 2, 3, false)]
#[case::nonzero_skipped(b"  1", 1, 3, true)]
fn test_validate_unused_n(
    #[case] field: &[u8],
    #[case] count: usize,
    #[case] width: usize,
    #[case] skip_unused_fields: bool,
) {
    assert_eq!(
        validate_unused_n(field, count, width, skip_unused_fields, 5),
        Ok(())
    );
}

#[rstest]
#[case::zeros_separated(b"0 0000", 2, 3, false, 5)]
#[case::too_short(b"000", 2, 3, false, 5)]
#[case::nonzero(b"  1", 1, 3, false, 5)]
#[case::too_short_skipped(b"000", 2, 3, true, 5)]
#[case::blank_too_short(b"     ", 2, 3, false, 5)]
fn test_validate_unused_n_error(
    #[case] field: &[u8],
    #[case] count: usize,
    #[case] width: usize,
    #[case] skip_unused_fields: bool,
    #[case] column: u32,
) {
    assert_eq!(
        validate_unused_n(field, count, width, skip_unused_fields, column),
        Err(ErrMode::Backtrack(InputError { column }))
    );
}
