use nom::character::complete::space0;
use nom::combinator::all_consuming;
use nom::sequence::delimited;
use nom::{error, Err};
use pretty_assertions::assert_eq;
use rstest::*;

use super::*;

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
fn test_is_whitespace_or_zeroes(#[case] input: &[u8], #[case] expected: bool) {
    assert_eq!(is_all_whitespace_or_zeroes(input), expected);
}

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
#[case(b"  1.2345  ", 1.2345)]
#[case(b"    -1.234", -1.234)]
#[case(b"1.0       ", 1.0)]
#[case(b"1.        ", 1.0)]
#[case(b" .1       ", 0.1)]
#[case(b"1.23456   ", 1.23456)]
#[case(b"   1234567", 123.4567)]
#[case(b"  -1234567", -123.4567)]
#[case(b"       123", 0.0123)]
#[case(b"          ", 0.0)]
fn test_fixed_width_float(#[case] input: &[u8], #[case] expected: f64) {
    let mut parser = all_consuming(fixed_width_float::<f64>(10, 4));
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
#[case(b"   C", Some(Element::C))]
#[case(b"  C ", Some(Element::C))]
#[case(b" C  ", Some(Element::C))]
#[case(b"C   ", Some(Element::C))]
#[case(b"Cu  ", Some(Element::Cu))]
#[case(b" Cu ", Some(Element::Cu))]
#[case(b"  Cu", Some(Element::Cu))]
#[case(b" Cu ", Some(Element::Cu))]
#[case(b"Cu  ", Some(Element::Cu))]
fn test_fixed_width_element_partial(#[case] input: &[u8], #[case] expected: Option<Element>) {
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
#[case(b"", 3, true, false)]
#[case(b"   ", 3, true, true)]
#[case(b"  0", 3, true, true)]
#[case(b" 0 ", 3, true, true)]
#[case(b"0  ", 3, true, true)]
#[case(b" 00", 3, true, true)]
#[case(b"00 ", 3, true, true)]
#[case(b"000", 3, true, true)]
#[case(b"0", 3, true, false)]
#[case(b"00", 3, true, false)]
#[case(b"0 0", 3, true, false)]
#[case(b"  1", 3, true, false)]
#[case(b"000", 3, false, true)]
#[case(b"0 0", 3, false, true)]
#[case(b"  1", 3, false, true)]
fn test_fixed_width_padding(
    #[case] input: &[u8],
    #[case] width: usize,
    #[case] strict_padding: bool,
    #[case] expected_success: bool,
) {
    let mut parser = all_consuming(fixed_width_padding(width, strict_padding));
    let result = parser.parse(input);
    assert_eq!(
        result.is_ok(),
        expected_success,
        "{:?}: should have {}",
        input,
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
#[case(b"", 0, 0, true, true)]
#[case(b"   ", 1, 3, true, true)]
#[case(b"  0", 1, 3, true, true)]
#[case(b" 0 ", 1, 3, true, true)]
#[case(b"0  ", 1, 3, true, true)]
#[case(b" 00", 1, 3, true, true)]
#[case(b"00 ", 1, 3, true, true)]
#[case(b"000", 1, 3, true, true)]
#[case(b"000000", 2, 3, true, true)]
#[case(b"   000", 2, 3, true, true)]
#[case(b"0 0000", 2, 3, true, false)]
#[case(b"  1", 1, 3, true, false)]
#[case(b"0 0000", 2, 3, false, true)]
#[case(b"  1", 1, 3, false, true)]
#[case(b"000001", 2, 3, true, false)]
fn test_fixed_width_padding_n(
    #[case] input: &[u8],
    #[case] count: usize,
    #[case] width: usize,
    #[case] strict_padding: bool,
    #[case] expected_success: bool,
) {
    let mut parser = all_consuming(fixed_width_padding_n(count, width, strict_padding));
    let result = parser.parse(input);
    assert_eq!(
        result.is_ok(),
        expected_success,
        "{:?}: should have succeeded",
        input
    );
    if expected_success {
        let (remaining, _) = result.unwrap();
        assert!(remaining.is_empty(), "remaining should be empty");
    }
}

#[rstest]
#[case(b"abcd", 4, Some("abcd".to_string()))]
#[case(b"abc ", 4, Some("abc".to_string()))]
#[case(b"abc", 4, Some("abc".to_string()))]
#[case(b" abc", 4, Some("abc".to_string()))]
#[case(b" ab ", 4, Some("ab".to_string()))]
#[case(b"", 4, None)]
#[case(b"   ", 4, None)]
fn test_fixed_width_str_partial(
    #[case] input: &[u8],
    #[case] width: usize,
    #[case] expected: Option<String>,
) {
    let mut parser = all_consuming(fixed_width_str_partial(width));
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
