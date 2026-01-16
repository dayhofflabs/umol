use nom::character::complete::space0;
use nom::combinator::all_consuming;
use nom::error::ErrorKind as NomErrorKind;
use nom::sequence::delimited;
use nom::Err;
use pretty_assertions::assert_eq;
use rstest::*;

use super::*;

#[rstest]
#[case::empty(b"", vec![])]
#[case::single_line_lf(b"abc\n", vec![(b"abc".as_slice(), 4)])]
#[case::single_line_crlf(b"abc\r\n", vec![(b"abc".as_slice(), 5)])]
#[case::single_line_no_term(b"abc", vec![(b"abc".as_slice(), 3)])]
#[case::two_lines_lf(b"abc\ndef\n", vec![(b"abc".as_slice(), 4), (b"def".as_slice(), 4)])]
#[case::two_lines_crlf(b"abc\r\ndef\r\n", vec![(b"abc".as_slice(), 5), (b"def".as_slice(), 5)])]
#[case::two_lines_mixed(b"abc\r\ndef\n", vec![(b"abc".as_slice(), 5), (b"def".as_slice(), 4)])]
#[case::two_lines_no_final_term(b"abc\ndef", vec![(b"abc".as_slice(), 4), (b"def".as_slice(), 3)])]
#[case::empty_line_lf(b"\n", vec![(b"".as_slice(), 1)])]
#[case::empty_line_crlf(b"\r\n", vec![(b"".as_slice(), 2)])]
#[case::two_empty_lines(b"\n\n", vec![(b"".as_slice(), 1), (b"".as_slice(), 1)])]
#[case::blank_line_between(b"a\n\nb\n", vec![(b"a".as_slice(), 2), (b"".as_slice(), 1), (b"b".as_slice(), 2)])]
fn test_lines_with_offset(#[case] input: &[u8], #[case] expected: Vec<(&[u8], usize)>) {
    let result: Vec<_> = input.lines_with_offset().collect();
    assert_eq!(result, expected);
}

#[test]
fn test_lines_with_offset_offset_sum() {
    let input = b"line1\r\nline2\nline3";
    let mut total = 0;
    for (_, len) in input.lines_with_offset() {
        total += len;
    }
    assert_eq!(total, input.len());
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
fn test_is_whitespace_or_zeroes(#[case] input: &[u8], #[case] expected: bool) {
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
#[case::non_numeric_input(b"abc", NomErrorKind::Digit)]
#[case::trailing_characters(b"1a ", NomErrorKind::Eof)]
fn test_fixed_width_partial_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let mut parser = fixed_width_partial(3, delimited(space0, nom_i32, space0), true);
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case::empty(b"", None)]
#[case::blank_width2(b"  ", None)]
#[case::blank_width3(b"   ", None)]
#[case::number(b" 12", Some(12))]
fn test_fixed_width_opt(#[case] input: &[u8], #[case] expected_val: Option<i32>) {
    let mut parser = fixed_width_opt(3, delimited(space0, nom_i32, space0));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(
        result.is_ok(),
        "{:?} should have succeeded but failed with {:?}",
        input_str,
        result
    );
    let (remaining, value) = result.unwrap();
    assert_eq!(
        value, expected_val,
        "{:?} has returned value {:?}, expected {:?}",
        input_str, value, expected_val
    );
    assert!(remaining.is_empty(), "remaining should be empty");
}

#[rstest]
#[case::non_numeric_input(b" abc ", NomErrorKind::Digit)]
#[case::two_characters(b" 1", NomErrorKind::Eof)]
#[case::one_character(b"1", NomErrorKind::Eof)]
fn test_fixed_width_opt_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let mut parser = fixed_width_opt(5, delimited(space0, nom_i32, space0));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(
        result.is_err(),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case::positive(b"123", 123i32)]
#[case::negative(b"-98", -98i32)]
#[case::padded(b"  8", 8i32)]
#[case::blank(b"   ", 0i32)]
#[case::blank_width2(b"  ", 0i32)]
fn test_fixed_width_int(#[case] input: &[u8], #[case] expected: i32) {
    let mut parser = all_consuming(fixed_width_int::<i32>(3));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(
        result.is_ok(),
        "{:?} should have succeeded but failed with {:?}",
        input_str,
        result
    );
    let (remaining, value) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(
        value, expected,
        "{:?} has returned value {:?}, expected {:?}",
        input_str, value, expected
    );
}

#[rstest]
#[case::four_characters(b"1234", NomErrorKind::Eof)]
#[case::two_characters(b"12", NomErrorKind::Eof)]
#[case::non_numeric_input(b"abc", NomErrorKind::Digit)]
#[case::trailing_characters(b"1a ", NomErrorKind::Eof)]
fn test_fixed_width_int_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let mut parser = all_consuming(fixed_width_int::<i32>(3));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case::positive(b"100", 100i8)]
#[case::negative(b" -9", -9i8)]
#[case::padded_right(b"8  ", 8i8)]
#[case::padded_both_sides(b" 1 ", 1i8)]
fn test_fixed_width_int_in_range(#[case] input: &[u8], #[case] expected: i8) {
    let mut parser = all_consuming(fixed_width_int_in_range::<i8, _>(3, -10i8..=110i8));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(
        result.is_ok(),
        "{:?} should have succeeded but failed with {:?}",
        input_str,
        result
    );
    let (remaining, value) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(
        value, expected,
        "{:?} has returned value {:?}, expected {:?}",
        input_str, value, expected
    );
}

#[rstest]
#[case::blank_not_in_range(b"   ", NomErrorKind::Verify)]
#[case::out_of_range(b"11 ", NomErrorKind::Verify)]
#[case::four_characters(b"1234", NomErrorKind::Verify)]
#[case::one_character(b"8", NomErrorKind::Eof)]
#[case::non_numeric_input(b"abc", NomErrorKind::Digit)]
#[case::trailing_characters(b"1a ", NomErrorKind::Eof)]
fn test_fixed_width_int_in_range_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let mut parser = all_consuming(fixed_width_int_in_range::<i8, _>(3, 1i8..=10i8));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case::positive(b"100", 100u8)]
#[case::padded_left(b"  9", 9u8)]
#[case::padded_right(b"8  ", 8u8)]
#[case::padded_both_sides(b" 1 ", 1u8)]
fn test_fixed_width_int_in_range_inclusive(#[case] input: &[u8], #[case] expected: u8) {
    let mut parser = all_consuming(fixed_width_int_in_range::<u8, _>(3, 0u8..=100u8));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(
        result.is_ok(),
        "{:?} should have succeeded but failed with {:?}",
        input_str,
        result
    );
    let (remaining, value) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(
        value, expected,
        "{:?} has returned value {:?}, expected {:?}",
        input_str, value, expected
    );
}

#[rstest]
#[case::padded_left_one_digit(b"  5", Some(5i32))]
#[case::padded_left_two_digits(b" 10", Some(10i32))]
#[case::padded_left_zero(b"  0", Some(0i32))]
#[case::two_digits_out_of_range(b" 11", None)]
#[case::negative_out_of_range(b" -1", None)]
#[case::blank(b"   ", None)]
#[case::blank_width2(b"  ", None)]
fn test_fixed_width_int_in_range_opt(#[case] input: &[u8], #[case] expected: Option<i32>) {
    let mut parser = all_consuming(fixed_width_int_in_range_opt::<i32, _>(3, 0..=10));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(
        result.is_ok(),
        "{:?} should have succeeded but failed with {:?}",
        input_str,
        result
    );
    let (remaining, value) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(
        value, expected,
        "{:?} has returned value {:?}, expected {:?}",
        input_str, value, expected
    );
}

#[rstest]
#[case::non_numeric_input(b"abc", NomErrorKind::Digit)]
#[case::trailing_characters(b"1a ", NomErrorKind::Eof)]
fn test_fixed_width_int_in_range_opt_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let mut parser = all_consuming(fixed_width_int_in_range_opt::<i32, _>(3, 0..=10));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case::padded_left(b"  1", 0usize)]
#[case::three_digits(b"123", 122usize)]
fn test_fixed_width_int_minus1(#[case] input: &[u8], #[case] expected: usize) {
    let mut parser = all_consuming(fixed_width_int_minus1::<usize>(3));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(
        result.is_ok(),
        "{:?} should have succeeded but failed with {:?}",
        input_str,
        result
    );
    let (remaining, value) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(
        value, expected,
        "{:?} has returned value {:?}, expected {:?}",
        input_str, value, expected
    );
}

#[rstest]
#[case::value_too_small(b"  0", NomErrorKind::Verify)]
fn test_fixed_width_int_minus1_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let mut parser = all_consuming(fixed_width_int_minus1::<usize>(3));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case::three_digits(b"123", 123i32)]
#[case::two_digits(b"12", 12i32)]
#[case::one_digit(b"1", 1i32)]
#[case::empty(b"", 0i32)]
#[case::two_digits_padded_right(b"12 ", 12i32)]
#[case::one_digit_padded_right(b"1  ", 1i32)]
#[case::two_digits_padded_left(b" 12", 12i32)]
#[case::one_digit_padded_left(b"  1", 1i32)]
#[case::two_digits_padded_both_sides(b" 1 ", 1i32)]
#[case::blank_width3(b"   ", 0i32)]
#[case::blank_width2(b"  ", 0i32)]
#[case::blank_width1(b" ", 0i32)]
#[case::negative(b" -1", -1i32)]
fn test_fixed_width_int_partial(#[case] input: &[u8], #[case] expected: i32) {
    let mut parser = all_consuming(fixed_width_int_partial::<i32>(3));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(
        result.is_ok(),
        "{:?} should have succeeded but failed with {:?}",
        input_str,
        result
    );
    let (remaining, value) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(
        value, expected,
        "{:?} has returned value {:?}, expected {:?}",
        input_str, value, expected
    );
}

#[rstest]
#[case::too_many_characters(b"1234", NomErrorKind::Eof)]
#[case::non_numeric_input(b"abc", NomErrorKind::Digit)]
#[case::trailing_characters(b"1a ", NomErrorKind::Eof)]
fn test_fixed_width_int_partial_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let mut parser = all_consuming(fixed_width_int_partial::<i32>(3));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
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
fn test_fixed_width_float(#[case] input: &[u8], #[case] expected: f64) {
    let mut parser = all_consuming(fixed_width_float::<f64>(10, 4));
    let result = parser.parse(input);
    let (_, parsed_val) = result.unwrap();
    assert!((parsed_val - expected).abs() < 1e-9);
}

#[rstest]
#[case::trailing_characters(b"1.23a     ", NomErrorKind::Eof)]
#[case::invalid_decimal_point(b"1.2.3     ", NomErrorKind::Eof)]
#[case::trailing_characters_after_blank(b"          a", NomErrorKind::Eof)]
fn test_fixed_width_float_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let mut parser = all_consuming(fixed_width_float::<f64>(10, 4));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
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
    let mut parser = all_consuming(fixed_width_element_partial(4));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, value) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(
        value, expected,
        "{:?} has returned value {:?}, expected {:?}",
        input_str, value, expected
    );
}

#[rstest]
#[case::trailing_characters(b"Cu   ", NomErrorKind::Eof)]
#[case::invalid_element_symbol(b" X  ", NomErrorKind::MapOpt)]
fn test_fixed_width_element_partial_invalid(
    #[case] input: &[u8],
    #[case] expected_kind: NomErrorKind,
) {
    let mut parser = all_consuming(fixed_width_element_partial(4));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case::empty(b"", 3, false, false)]
#[case::blank(b"   ", 3, false, true)]
#[case::zero(b"  0", 3, false, true)]
#[case::zero_pos2(b" 0 ", 3, false, true)]
#[case::zero_pos1(b"0  ", 3, false, true)]
#[case::two_zeros(b" 00", 3, false, true)]
#[case::two_zeros_pos1(b"00 ", 3, false, true)]
#[case::three_zeros(b"000", 3, false, true)]
#[case::zero_width1(b"0", 3, false, false)]
#[case::two_zeros_width2(b"00", 3, false, false)]
#[case::two_zeros_separated(b"0 0", 3, false, false)]
#[case::one(b"  1", 3, false, false)]
#[case::three_zeros_skip_unused_fields(b"000", 3, true, true)]
#[case::two_zeros_separated_skip_unused_fields(b"0 0", 3, true, true)]
#[case::one_skip_unused_fields(b"  1", 3, true, true)]
#[case::non_numeric(b" a ", 3, false, false)]
#[case::blank_width2(b"  ", 3, false, false)]
fn test_fixed_width_unused(
    #[case] input: &[u8],
    #[case] width: usize,
    #[case] skip_unused_fields: bool,
    #[case] expected_success: bool,
) {
    let mut parser = all_consuming(fixed_width_unused(width, skip_unused_fields));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert_eq!(
        result.is_ok(),
        expected_success,
        "{:?}: should have {}",
        input_str,
        if expected_success {
            "succeeded"
        } else {
            "failed"
        }
    );
    if expected_success {
        let (remaining, _) = result.unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
    }
}

#[rstest]
#[case::empty(b"", 0, 0, false, true)]
#[case::blank(b"   ", 1, 3, false, true)]
#[case::zero(b"  0", 1, 3, false, true)]
#[case::zero_pos2(b" 0 ", 1, 3, false, true)]
#[case::zero_pos1(b"0  ", 1, 3, false, true)]
#[case::two_zeros(b" 00", 1, 3, false, true)]
#[case::two_zeros_pos1(b"00 ", 1, 3, false, true)]
#[case::three_zeros(b"000", 1, 3, false, true)]
#[case::two_fields_zeros(b"000000", 2, 3, false, true)]
#[case::blank_field_zero(b"   000", 2, 3, false, true)]
#[case::zeros_separated(b"0 0000", 2, 3, false, false)]
#[case::two_fields_too_short(b"000", 2, 3, false, false)]
#[case::one(b"  1", 1, 3, false, false)]
#[case::zeros_separated_skip_unused_fields(b"0 0000", 2, 3, true, true)]
#[case::one_skip_unused_fields(b"  1", 1, 3, true, true)]
#[case::one_zero_padded_left(b"000001", 2, 3, false, false)]
#[case::two_fields_too_short_skip_unused_fields(b"000", 2, 3, true, false)]
#[case::two_fields_blank_width5(b"     ", 2, 3, false, false)]
fn test_fixed_width_unused_n(
    #[case] input: &[u8],
    #[case] count: usize,
    #[case] width: usize,
    #[case] skip_unused_fields: bool,
    #[case] expected_success: bool,
) {
    let mut parser = all_consuming(fixed_width_unused_n(count, width, skip_unused_fields));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert_eq!(
        result.is_ok(),
        expected_success,
        "{:?}: should have {}",
        input_str,
        if expected_success {
            "succeeded"
        } else {
            "failed"
        }
    );
    if expected_success {
        let (remaining, _) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "remaining should be empty for {:?}",
            input_str
        );
    }
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
    let mut parser = all_consuming(fixed_width_str_partial(width));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_ok(), "{:?} should have succeeded", input_str);
    let (remaining, result) = result.unwrap();
    assert!(remaining.is_empty(), "remaining should be empty");
    assert_eq!(result, expected);
}

#[rstest]
#[case::zero(b"    0.0000    0.0000    0.0000", false, false, Point3D::zero())]
#[case::nonzero(b"    1.2345   -2.3456    3.4567", false, false, Point3D::new(1.2345, -2.3456, 3.4567))]
#[case::tens(b"   10.0000  -20.0000   30.0000", false, false, Point3D::new(10.0, -20.0, 30.0))]
#[case::hundreds(b"  123.4567 -234.5678  345.6789", false, false, Point3D::new(123.4567, -234.5678, 345.6789))]
#[case::short_zero(b"       0.0       0.0       0.0", false, false, Point3D::zero())]
#[case::no_integer_part(b"        .0        .0        .0", false, false, Point3D::zero())]
#[case::zero_ignored(b"    0.0000    0.0000    0.0000", true, false, Point3D::zero())]
#[case::nonzero_ignored(b"    1.2345   -2.3456    3.4567", true, false, Point3D::zero())]
#[case::zero_skipped(b"    0.0000    0.0000    0.0000", true, true, Point3D::zero())]
#[case::invalid_skipped(b"    x.0000    0.0000    0.0000", true, true, Point3D::zero())]
fn test_fixed_width_position(
    #[case] input: &[u8],
    #[case] ignore_positions: bool,
    #[case] skip_unused_fields: bool,
    #[case] expected: Point3D,
) {
    let mut parser = all_consuming(fixed_width_position(ignore_positions, skip_unused_fields));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(
        result.is_ok(),
        "{:?} should have succeeded, error: {:?}",
        input_str,
        result.clone().unwrap_err()
    );
    let (remaining, pos) = result.unwrap();
    assert!(
        remaining.is_empty(),
        "{:?} should have consumed all input, remaining: {:?}",
        input_str,
        remaining
    );
    assert_eq!(pos, expected);
}

#[rstest]
#[case::too_short(b"    0.0000    0.0000    0.000", false, false, NomErrorKind::Eof)]
#[case::too_long(b"    0.0000    0.0000    0.00000", false, false, NomErrorKind::Eof)]
#[case::invalid_x(b"    x.0000    0.0000    0.0000", false, false, NomErrorKind::Digit)]
#[case::invalid_y(b"    0.0000    y.0000    0.0000", false, false, NomErrorKind::Digit)]
#[case::invalid_z(b"    0.0000    0.0000    z.0000", false, false, NomErrorKind::Digit)]
#[case::invalid_x_unused(b"    x.0000    0.0000    0.0000", true, false, NomErrorKind::Digit)]
#[case::invalid_y_unused(b"    0.0000    y.0000    0.0000", true, false, NomErrorKind::Digit)]
#[case::invalid_z_unused(b"    0.0000    0.0000    z.0000", true, false, NomErrorKind::Digit)]
fn test_fixed_width_position_invalid(
    #[case] input: &[u8],
    #[case] ignore_positions: bool,
    #[case] skip_unused_fields: bool,
    #[case] expected_kind: NomErrorKind,
) {
    let mut parser = all_consuming(fixed_width_position(ignore_positions, skip_unused_fields));
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(
        result.is_err(),
        "{:?} should have failed, output: {:?}",
        input_str,
        result.clone().unwrap()
    );
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code)
    );
}

#[rstest]
#[case::empty(b"", vec![RGroupOccurrence::GreaterThan(0)])]
#[case::one(b"1", vec![RGroupOccurrence::Exactly(1)])]
#[case::two(b"1,2", vec![RGroupOccurrence::Exactly(1), RGroupOccurrence::Exactly(2)])]
#[case::greater_than(b">1", vec![RGroupOccurrence::GreaterThan(1)])]
#[case::fewer_than(b"<2", vec![RGroupOccurrence::FewerThan(2)])]
#[case::range(b"1-3", vec![RGroupOccurrence::Range(1, 3)])]
#[case::exactly_greater_than(b"0,>0", vec![RGroupOccurrence::Exactly(0), RGroupOccurrence::GreaterThan(0)])]
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
            allow_rgroups
        ),
        expected
    );
}

#[rstest]
#[case::invalid_character(b"a", NomErrorKind::Eof)]
#[case::negative_value(b"-3", NomErrorKind::Eof)]
fn test_rgroup_occurrences_invalid(#[case] input: &[u8], #[case] expected_kind: NomErrorKind) {
    let mut parser = all_consuming(rgroup_occurrences());
    let result = parser.parse(input);
    let input_str = input.to_str_lossy();
    assert!(result.is_err(), "{:?} should have failed", input_str);
    assert!(
        matches!(result.clone(), Err(Err::Error(e)) if e.code == expected_kind),
        "{:?} should have failed with error kind {:?}, got {:?}",
        input_str,
        expected_kind,
        result.clone().unwrap_err().map(|e| e.code),
    );
}
