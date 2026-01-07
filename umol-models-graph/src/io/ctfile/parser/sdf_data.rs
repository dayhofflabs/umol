//! SDF data field parsers

use bstr::ByteSlice;
use indexmap::IndexMap;
use nom::branch::alt;
use nom::bytes::complete::{tag, take_until};
use nom::character::complete::{line_ending, not_line_ending, space1};
use nom::combinator::{eof, map, opt, peek, value};
use nom::error::Error as NomError;
use nom::multi::many_till;
use nom::sequence::{delimited, terminated};
use nom::Parser;

/// Parse data field header: `> <Field Name>`
pub(super) fn sdf_data_header<'inp>(
) -> impl Parser<&'inp [u8], Output = String, Error = NomError<&'inp [u8]>> {
    map(
        delimited(
            (tag(">"), space1),
            delimited(tag("<"), take_until(">"), tag(">")),
            not_line_ending,
        ),
        |field_name: &[u8]| field_name.to_str_lossy().into_owned(),
    )
}

/// Parse multi-line data value until blank line
pub(super) fn sdf_data_value<'inp>(
) -> impl Parser<&'inp [u8], Output = String, Error = NomError<&'inp [u8]>> {
    map(
        many_till(
            terminated(not_line_ending, line_ending),
            // TODO: Why is this a blank line?
            alt((
                peek(line_ending), // blank line
                eof,
            )),
        ),
        // TODO: Understand why this is necessary.
        |(lines, _): (Vec<&[u8]>, _)| {
            lines
                .iter()
                .map(|line| line.to_str_lossy())
                .collect::<Vec<_>>()
                .join("\n")
        },
    )
}

/// Parse complete data field (header + value)
pub(super) fn sdf_data_field<'inp>(
) -> impl Parser<&'inp [u8], Output = (String, String), Error = NomError<&'inp [u8]>> {
    (terminated(sdf_data_header(), line_ending), sdf_data_value())
}

/// Parse SDF record delimiter
pub(super) fn sdf_delimiter<'inp>(
) -> impl Parser<&'inp [u8], Output = (), Error = NomError<&'inp [u8]>> {
    value((), (tag("$$$$"), opt(line_ending)))
}

/// Parse multiple data fields
pub(super) fn sdf_data_block<'inp>(
) -> impl Parser<&'inp [u8], Output = IndexMap<String, String>, Error = NomError<&'inp [u8]>> {
    map(
        many_till(
            alt((
                sdf_data_field(),
                // Skip non-data lines
                value(
                    (String::new(), String::new()),
                    terminated(not_line_ending, line_ending),
                ),
            )),
            peek(alt((tag("$$$$"), eof))),
        ),
        |(fields, _): (Vec<(String, String)>, _)| {
            fields
                .into_iter()
                .filter(|(name, _)| !name.is_empty())
                .collect::<IndexMap<_, _>>()
        },
    )
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;

    use super::*;

    #[test]
    fn test_sdf_data_header() {
        let input = b"> <MELTING.POINT>\n";
        let result = sdf_data_header().parse(input);
        assert!(result.is_ok());
        let (_, field_name) = result.unwrap();
        assert_eq!(field_name, "MELTING.POINT");
    }

    #[test]
    fn test_sdf_data_header_multiple_spaces() {
        let input = b">  <BOILING.POINT>\n";
        let result = sdf_data_header().parse(input);
        assert!(result.is_ok());
        let (_, field_name) = result.unwrap();
        assert_eq!(field_name, "BOILING.POINT");
    }

    #[test]
    fn test_sdf_data_value() {
        let input = b"100.5\n\n";
        let result = sdf_data_value().parse(input);
        assert!(result.is_ok());
        let (_, value) = result.unwrap();
        assert_eq!(value, "100.5");
    }

    #[test]
    fn test_sdf_data_value_multiline() {
        let input = b"Line 1\nLine 2\nLine 3\n\n";
        let result = sdf_data_value().parse(input);
        assert!(result.is_ok());
        let (_, value) = result.unwrap();
        assert_eq!(value, "Line 1\nLine 2\nLine 3");
    }

    #[test]
    fn test_sdf_data_field() {
        let input = b"> <TEST_FIELD>\n100.5\n\n";
        let result = sdf_data_field().parse(input);
        assert!(result.is_ok());
        let (_, (field_name, field_value)) = result.unwrap();
        assert_eq!(field_name, "TEST_FIELD");
        assert_eq!(field_value, "100.5");
    }

    #[test]
    fn test_sdf_delimiter() {
        let input = b"$$$$\n";
        let result = sdf_delimiter().parse(input);
        assert!(result.is_ok());
    }

    #[test]
    fn test_sdf_delimiter_no_newline() {
        let input = b"$$$$";
        let result = sdf_delimiter().parse(input);
        assert!(result.is_ok());
    }

    #[test]
    fn test_sdf_data_block() {
        let input = b"> <FIELD1>\nValue1\n\n> <FIELD2>\nValue2\n\n$$$$\n";
        let result = sdf_data_block().parse(input);
        assert!(result.is_ok());
        let (remaining, fields): (&[u8], IndexMap<String, String>) = result.unwrap();
        assert_eq!(fields.len(), 2);
        assert_eq!(fields.get("FIELD1"), Some(&"Value1".to_string()));
        assert_eq!(fields.get("FIELD2"), Some(&"Value2".to_string()));
        assert!(remaining.starts_with(b"$$$$"));
    }
}
