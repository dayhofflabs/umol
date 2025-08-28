use super::*;
use nom::character::complete::space0;
use nom::combinator::all_consuming;
use nom::sequence::delimited;
use nom::{error, Err};
use pretty_assertions::assert_eq;
use rstest::*;

#[rstest]
#[case(b"  0", false, b"0")]
#[case(b" 0 ", false, b"0")]
#[case(b"0  ", false, b"0")]
#[case("\u{00A0}\u{00A0}0".as_bytes(), false, "\u{00A0}\u{00A0}0".as_bytes())]
#[case("0\u{00A0}\u{00A0}".as_bytes(), false, "0\u{00A0}\u{00A0}".as_bytes())]
#[case("\u{00A0}0\u{00A0}".as_bytes(), false, "\u{00A0}0\u{00A0}".as_bytes())]
#[case(b"  0", true, b"0")]
#[case(b" 0 ", true, b"0")]
#[case(b"0  ", true, b"0")]
#[case("\u{00A0}\u{00A0}0".as_bytes(), true, b"0")]
#[case("0\u{00A0}\u{00A0}".as_bytes(), true, b"0")]
#[case("\u{00A0}0\u{00A0}".as_bytes(), true, b"0")]
fn test_trim_whitespace(
    #[case] input: &[u8],
    #[case] allow_unicode: bool,
    #[case] expected: &[u8],
) {
    assert_eq!(trim_whitespace(input, allow_unicode), expected);
}

#[rstest]
#[case(b"", false, true)]
#[case(b"   ", false, true)]
#[case(b"  0", false, true)]
#[case(b" 0 ", false, true)]
#[case(b"0  ", false, true)]
#[case(b" 00", false, true)]
#[case(b"00 ", false, true)]
#[case(b"000", false, true)]
#[case(b"0", false, true)]
#[case(b"00", false, true)]
#[case(b"0 0", false, false)]
#[case(b"  1", false, false)]
#[case(b"", true, true)]
#[case(b"   ", true, true)]
#[case(b"  0", true, true)]
#[case(b" 0 ", true, true)]
#[case(b"0  ", true, true)]
#[case(b" 00", true, true)]
#[case(b"00 ", true, true)]
#[case(b"000", true, true)]
#[case(b"0", true, true)]
#[case(b"00", true, true)]
#[case(b"0 0", true, false)]
#[case(b"  1", true, false)]
#[case("\u{00A0}\u{00A0}0".as_bytes(), true, true)]
fn test_is_whitespace_or_zeroes(
    #[case] input: &[u8],
    #[case] allow_unicode: bool,
    #[case] expected: bool,
) {
    assert_eq!(is_all_whitespace_or_zeroes(input, allow_unicode), expected);
}

#[rstest]
#[case(b"", false, None)]
#[case(b"  ", false, None)]
#[case(b"   ", false, None)]
#[case(b"42", false, Some(42))]
#[case(b" 42", false, Some(42))]
#[case(b"42 ", false, Some(42))]
#[case(b"042", false, Some(42))]
#[case("2\u{00A0}".as_bytes(), true, Some(2))]
fn test_fixed_width_partial(
    #[case] input: &[u8],
    #[case] allow_unicode: bool,
    #[case] expected_val: Option<i32>,
) {
    let mut parser = fixed_width_partial(
        3,
        move |input| {
            let s = trim_whitespace(input, allow_unicode);
            nom_i32.parse(s)
        },
        true,
        allow_unicode,
    );

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
#[case("2\u{00A0}".as_bytes(), "unicode whitespace", error::ErrorKind::Eof)]
fn test_fixed_width_partial_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let mut parser = fixed_width_partial(
        3,
        move |input| {
            let s = trim_whitespace(input, false);
            nom_i32.parse(s)
        },
        true,
        false,
    );
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
#[case(b"", false, None)]
#[case(b"  ", false, None)]
#[case(b"   ", false, None)]
#[case(b" 42", false, Some(42))]
#[case("2\u{00A0}".as_bytes(), true, Some(2))]
fn test_fixed_width_opt(
    #[case] input: &[u8],
    #[case] allow_unicode: bool,
    #[case] expected_val: Option<i32>,
) {
    let mut parser = fixed_width_opt(3,
        move |input| {
            let s = trim_whitespace(input, allow_unicode);
            nom_i32.parse(s)
        },
        allow_unicode,
    );
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
#[case("2\u{00A0}".as_bytes(), "unicode whitespace", error::ErrorKind::Eof)]
fn test_fixed_width_opt_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let mut parser = fixed_width_opt(5, delimited(space0, nom_i32, space0), false);
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
#[case(b"123", false, 123i32)]
#[case(b"-98", false, -98i32)]
#[case(b"  8", false, 8i32)]
#[case(b"   ", false, 0i32)]
#[case("2\u{00A0}".as_bytes(), true, 2i32)]
fn test_fixed_width_int(#[case] input: &[u8], #[case] allow_unicode: bool, #[case] expected: i32) {
    let mut parser = all_consuming(fixed_width_int::<i32>(3, allow_unicode));
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
#[case("2\u{00A0}".as_bytes(), "unicode whitespace", error::ErrorKind::Eof)]
fn test_fixed_width_int_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let mut parser = all_consuming(fixed_width_int::<i32>(3, false));
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
#[case(b"100", false, 100i8)]
#[case(b" -9", false, -9i8)]
#[case(b"8  ", false, 8i8)]
#[case(b" 1 ", false, 1i8)]
#[case("2\u{00A0}".as_bytes(), true, 2i8)]
fn test_fixed_width_int_in_range(
    #[case] input: &[u8],
    #[case] allow_unicode: bool,
    #[case] expected: i8,
) {
    let mut parser = all_consuming(fixed_width_int_in_range::<i8, _>(
        3,
        -10i8..=110i8,
        allow_unicode,
    ));
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
#[case("2\u{00A0}".as_bytes(), "unicode whitespace", error::ErrorKind::Eof)]
fn test_fixed_width_int_in_range_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let mut parser = all_consuming(fixed_width_int_in_range::<i8, _>(3, 1i8..=10i8, false));
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
#[case(b"100", false, 100u8)]
#[case(b"  9", false, 9u8)]
#[case(b"8  ", false, 8u8)]
#[case(b" 1 ", false, 1u8)]
#[case("2\u{00A0}".as_bytes(), true, 2u8)]
fn test_fixed_width_int_in_range_inclusive(
    #[case] input: &[u8],
    #[case] allow_unicode: bool,
    #[case] expected: u8,
) {
    let mut parser = all_consuming(fixed_width_int_in_range::<u8, _>(
        3,
        0u8..=100u8,
        allow_unicode,
    ));
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
#[case(b"  5", false, Some(5i32))]
#[case(b" 10", false, Some(10i32))]
#[case(b"  0", false, Some(0i32))]
#[case(b" 11", false, None)]
#[case(b" -1", false, None)]
#[case(b"   ", false, None)]
#[case(b"  ", false, None)]
#[case(b"", false, None)]
#[case("2\u{00A0}".as_bytes(), true, Some(2i32))]
fn test_fixed_width_int_in_range_opt(
    #[case] input: &[u8],
    #[case] allow_unicode: bool,
    #[case] expected: Option<i32>,
) {
    let mut parser = all_consuming(fixed_width_int_in_range_opt::<i32, _>(
        3,
        0..=10,
        allow_unicode,
    ));
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
#[case("2\u{00A0}".as_bytes(), "unicode whitespace", error::ErrorKind::Eof)]
fn test_fixed_width_int_in_range_opt_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let mut parser = all_consuming(fixed_width_int_in_range_opt::<i32, _>(3, 0..=10, false));
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
#[case(b"  1", false, 0usize)]
#[case(b"123", false, 122usize)]
#[case("2\u{00A0}".as_bytes(), true, 1usize)]
fn test_fixed_width_int_minus1(
    #[case] input: &[u8],
    #[case] allow_unicode: bool,
    #[case] expected: usize,
) {
    let mut parser = all_consuming(fixed_width_int_minus1::<usize>(3, allow_unicode));
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
#[case("2\u{00A0}".as_bytes(), "unicode whitespace", error::ErrorKind::Eof)]
fn test_fixed_width_int_minus1_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let mut parser = all_consuming(fixed_width_int_minus1::<usize>(3, false));
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
#[case(b"123", false, 123i32)]
#[case(b"12", false, 12i32)]
#[case(b"1", false, 1i32)]
#[case(b"", false, 0i32)]
#[case(b"12 ", false, 12i32)]
#[case(b"1  ", false, 1i32)]
#[case(b" 12", false, 12i32)]
#[case(b"  1", false, 1i32)]
#[case(b" 1 ", false, 1i32)]
#[case(b"   ", false, 0i32)]
#[case(b"  ", false, 0i32)]
#[case(b" ", false, 0i32)]
#[case(b" -1", false, -1i32)]
#[case("2\u{00A0}".as_bytes(), true, 2i32)]
fn test_fixed_width_int_partial(
    #[case] input: &[u8],
    #[case] allow_unicode: bool,
    #[case] expected: i32,
) {
    let mut parser = all_consuming(fixed_width_int_partial::<i32>(3, allow_unicode));
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
#[case("2\u{00A0}".as_bytes(), "unicode whitespace", error::ErrorKind::Eof)]
fn test_fixed_width_int_partial_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let mut parser = all_consuming(fixed_width_int_partial::<i32>(3, false));
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
#[case(b"  1.2345  ", false, 1.2345)]
#[case(b"    -1.234", false, -1.234)]
#[case(b"1.0       ", false, 1.0)]
#[case(b"1.        ", false, 1.0)]
#[case(b"1.23456   ", false, 1.23456)]
#[case(b"   1234567", false, 123.4567)]
#[case(b"  -1234567", false, -123.4567)]
#[case(b"       123", false, 0.0123)]
#[case(b"          ", false, 0.0)]
#[case("1.234567\u{00A0}".as_bytes(), true, 1.234567)]
fn test_fixed_width_float(
    #[case] input: &[u8],
    #[case] allow_unicode: bool,
    #[case] expected: f64,
) {
    let mut parser = all_consuming(fixed_width_float::<f64>(10, 4, allow_unicode)); // precision is ignored here
    let result = parser.parse(input);
    let (_, parsed_val) = result.unwrap();
    assert!((parsed_val - expected).abs() < 1e-9);
}

#[rstest]
#[case(b"1.23a     ", "trailing characters", error::ErrorKind::Eof)]
#[case(b"1.2.3     ", "invalid decimal point", error::ErrorKind::Eof)]
#[case(b"          a", "trailing characters", error::ErrorKind::Eof)]
#[case("1.234567\u{00A0}".as_bytes(), "unicode whitespace", error::ErrorKind::Eof)]
fn test_fixed_width_float_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let mut parser = all_consuming(fixed_width_float::<f64>(10, 4, false));
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
#[case(b"   C", false, Element::C)]
#[case(b"  C ", false, Element::C)]
#[case(b" C  ", false, Element::C)]
#[case(b"C   ", false, Element::C)]
#[case(b"Cu  ", false, Element::Cu)]
#[case(b" Cu ", false, Element::Cu)]
#[case(b"  Cu", false, Element::Cu)]
#[case(b" Cu ", false, Element::Cu)]
#[case(b"Cu  ", false, Element::Cu)]
#[case("Cu\u{00A0}".as_bytes(), true, Element::Cu)]
fn test_fixed_width_element(
    #[case] input: &[u8],
    #[case] allow_unicode: bool,
    #[case] expected: Element,
) {
    let mut parser = all_consuming(fixed_width_element(4, allow_unicode));
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
#[case("Cu\u{00A0}".as_bytes(), "unicode whitespace", error::ErrorKind::MapOpt)]
fn test_fixed_width_element_invalid(
    #[case] input: &[u8],
    #[case] desc: &str,
    #[case] expected_kind: error::ErrorKind,
) {
    let mut parser = all_consuming(fixed_width_element(4, false));
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
#[case(b"", 3, false, true, false)]
#[case(b"   ", 3, false, true, true)]
#[case(b"  0", 3, false, true, true)]
#[case(b" 0 ", 3, false, true, true)]
#[case(b"0  ", 3, false, true, true)]
#[case(b" 00", 3, false, true, true)]
#[case(b"00 ", 3, false, true, true)]
#[case(b"000", 3, false, true, true)]
#[case(b"0", 3, false, true, false)]
#[case(b"00", 3, false, true, false)]
#[case(b"0 0", 3, false, true, false)]
#[case(b"  1", 3, false, true, false)]
#[case(b"000", 3, false, false, true)]
#[case(b"0 0", 3, false, false, true)]
#[case(b"  1", 3, false, false, true)]
#[case("\u{00A0}0".as_bytes(), 3, true, false, true)]
fn test_fixed_width_padding(
    #[case] input: &[u8],
    #[case] width: usize,
    #[case] allow_unicode: bool,
    #[case] strict_padding: bool,
    #[case] expected_success: bool,
) {
    let mut parser = all_consuming(fixed_width_padding(width, allow_unicode, strict_padding));
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
#[case(b"", 0, 0, false, true, true)]
#[case(b"   ", 1, 3, false, true, true)]
#[case(b"  0", 1, 3, false, true, true)]
#[case(b" 0 ", 1, 3, false, true, true)]
#[case(b"0  ", 1, 3, false, true, true)]
#[case(b" 00", 1, 3, false, true, true)]
#[case(b"00 ", 1, 3, false, true, true)]
#[case(b"000", 1, 3, false, true, true)]
#[case(b"000000", 2, 3, false, true, true)]
#[case(b"   000", 2, 3, false, true, true)]
#[case(b"0 0000", 2, 3, false, true, false)]
#[case(b"  1", 1, 3, false, true, false)]
#[case(b"0 0000", 2, 3, false, false, true)]
#[case(b"  1", 1, 3, false, false, true)]
#[case(b"000001", 2, 3, false, true, false)]
#[case("\u{00A0}0\u{00A0}0".as_bytes(), 2, 3, true, false, true)]
fn test_fixed_width_padding_n(
    #[case] input: &[u8],
    #[case] count: usize,
    #[case] width: usize,
    #[case] allow_unicode: bool,
    #[case] strict_padding: bool,
    #[case] expected_success: bool,
) {
    let mut parser = all_consuming(fixed_width_padding_n(
        count,
        width,
        allow_unicode,
        strict_padding,
    ));
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
#[case(b"abcd", 4, false, Some("abcd".to_string()))]
#[case(b"abc ", 4, false, Some("abc".to_string()))]
#[case(b"abc", 4, false, Some("abc".to_string()))]
#[case(b" abc", 4, false, Some("abc".to_string()))]
#[case(b" ab ", 4, false, Some("ab".to_string()))]
#[case(b"", 4, false, None)]
#[case(b"   ", 4, false, None)]
#[case("ab\u{00A0}".as_bytes(), 4, true, Some("ab".to_string()))]
fn test_fixed_width_str_partial(
    #[case] input: &[u8],
    #[case] width: usize,
    #[case] allow_unicode: bool,
    #[case] expected: Option<String>,
) {
    let mut parser = all_consuming(fixed_width_str_partial(width, allow_unicode));
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
    let mut parser = all_consuming(rgroup_occurrences(false));
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
    let mut parser = all_consuming(rgroup_occurrences(false));
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
