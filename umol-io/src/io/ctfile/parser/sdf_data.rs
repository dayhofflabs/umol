//! SDF data field parsers

use bstr::{join, ByteSlice};
use indexmap::IndexMap;
use nom::bytes::complete::{is_not, tag, take_until1};
use nom::character::complete::multispace0;
use nom::combinator::{map, value};
use nom::sequence::{delimited, preceded};
use nom::{Err, Parser};

use super::utils::LinesWithOffsetExt;
use crate::io::ctfile::error::ParseError;

/// Parse data field header: `> <Field Name>`
pub(super) fn sdf_data_header<'inp>(
    line_offset: u32,
) -> impl Parser<&'inp [u8], Output = (String, u32), Error = ParseError> + use<'inp> {
    move |input: &'inp [u8]| {
        let (line, byte_len) =
            input
                .lines_with_offset()
                .next()
                .ok_or(Err::Error(ParseError::UnexpectedEof {
                    line: line_offset,
                    block: "sdf data",
                }))?;
        let mut name_input = map(
            preceded(
                (tag(">"), is_not("<")),
                delimited(tag("<"), take_until1(">"), tag(">")),
            ),
            move |name: &[u8]| name.trim().to_str_lossy().into_owned(),
        );
        let (_, name) = name_input.parse(line)?;
        let remaining = &input[byte_len..];
        Ok((remaining, (name, line_offset + 1)))
    }
}

/// Parse multi-line data value until blank line
pub(super) fn sdf_data_value<'inp>(
    line_offset: u32,
) -> impl Parser<&'inp [u8], Output = (String, u32), Error = ParseError> + use<'inp> {
    move |input: &'inp [u8]| {
        let mut byte_offset = 0;
        let mut line_index = 0;
        let mut value_lines = Vec::new();

        for (line, byte_len) in input.lines_with_offset() {
            byte_offset += byte_len;
            line_index += 1;
            if line.trim().is_empty() {
                break;
            }

            value_lines.push(line.trim())
        }
        let value = join(",", value_lines).to_str_lossy().to_string();
        let remaining = &input[byte_offset..];
        Ok((remaining, (value, line_offset + line_index)))
    }
}

/// Parse complete data field (header + value)
pub(super) fn sdf_data_field<'inp>(
    line_offset: u32,
) -> impl Parser<&'inp [u8], Output = ((String, String), u32), Error = ParseError> + use<'inp> {
    move |input: &'inp [u8]| {
        let (remaining, (name, line_offset)) = sdf_data_header(line_offset).parse(input)?;
        let (remaining, (data, line_offset)) = sdf_data_value(line_offset).parse(remaining)?;
        Ok((remaining, ((name, data), line_offset)))
    }
}

/// Parse SDF record delimiter
pub(super) fn sdf_delimiter<'inp>(
    line_offset: u32,
) -> impl Parser<&'inp [u8], Output = ((), u32), Error = ParseError> {
    value(((), line_offset + 1), (tag("$$$$"), multispace0))
}

/// Parse multiple data fields
pub(super) fn sdf_data_block<'inp>(
    line_offset: u32,
) -> impl Parser<&'inp [u8], Output = (IndexMap<String, String>, u32), Error = ParseError> + use<'inp>
{
    move |input: &'inp [u8]| {
        let mut remaining = input;
        let mut remaining_offset = line_offset;
        let mut data = IndexMap::new();
        loop {
            if remaining.is_empty() {
                return Err(Err::Error(ParseError::UnexpectedEof {
                    line: line_offset,
                    block: "sdf data",
                }));
            }
            if let Ok((new_remaining, ((name, value), new_line_offset))) =
                sdf_data_field(remaining_offset).parse(remaining)
            {
                data.insert(name, value);
                remaining = new_remaining;
                remaining_offset = new_line_offset;
            } else if let Ok((new_remaining, (_, new_line_offset))) =
                sdf_delimiter(line_offset).parse(remaining)
            {
                return Ok((new_remaining, (data, new_line_offset)));
            } else {
                return Err(Err::Error(ParseError::InvalidSdfDataHeader {
                    line: line_offset + 1,
                }));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use bstr::ByteSlice;
    use indexmap::indexmap;
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rstest]
    #[case::dotted(b"> <MELTING.POINT>\n", "MELTING.POINT".to_string())]
    #[case::whitespace(b"> <CAS NR>\n", "CAS NR".to_string())]
    #[case::multiple_space(b">  <BOILING.POINT>\n", "BOILING.POINT".to_string())]
    #[case::interstitial_data(b"> (MD-0894) <CAS NR>\n", "CAS NR".to_string())]
    #[case::trailing_data(b"> <CAS NR> DT12\n", "CAS NR".to_string())]
    #[case::surrounding_data(b"> (MD-0894) <BOILING.POINT> FROM ARCHIVES\n", "BOILING.POINT".to_string())]
    fn test_sdf_data_header(#[case] input: &[u8], #[case] expected: String) {
        let result = sdf_data_header(0).parse(input);
        let input_str = input.to_str_lossy();
        assert!(
            result.is_ok(),
            "{:?} should have succeeded, result: {:?}",
            input_str,
            result
        );
        let (remaining, (name, _)) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "{:?} leaved unparsed data: {:?}",
            input_str,
            remaining
        );
        assert_eq!(
            name, expected,
            "{:?}: name {:?} != expected {:?}",
            input_str, name, expected
        );
    }

    #[rstest]
    #[case::single_line(b"100.5\n\n", "100.5".to_string())]
    #[case::whitespace(b" 100.5 \n\n", "100.5".to_string())]
    #[case::multiple_lines(b"benzene\nBenzol\n\n", "benzene,Benzol".to_string())]
    fn test_sdf_data_value(#[case] input: &[u8], #[case] expected: String) {
        let result = sdf_data_value(0).parse(input);
        let input_str = input.to_str_lossy();
        assert!(
            result.is_ok(),
            "{:?} should have succeeded, result: {:?}",
            input_str,
            result
        );
        let (remaining, (value, _)) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "{:?} leaves unparsed data: {:?}",
            input_str,
            remaining
        );
        assert_eq!(
            value, expected,
            "{:?}: value {:?} != expected {:?}",
            input_str, value, expected
        );
    }

    #[rstest]
    #[case::single_line(b"> <BOILING.POINT>\n100.5\n\n", "BOILING.POINT".to_string(), "100.5".to_string())]
    #[case::multiple_line(b"> <NAMES>\nbenzene\nBenzol\n\n", "NAMES".to_string(), "benzene,Benzol".to_string())]
    fn test_sdf_data_field(
        #[case] input: &[u8],
        #[case] expected_name: String,
        #[case] expected_value: String,
    ) {
        let result = sdf_data_field(0).parse(input);
        let input_str = input.to_str_lossy();
        assert!(
            result.is_ok(),
            "{:?} should have succeeded, result: {:?}",
            input_str,
            result
        );
        let (remaining, ((name, value), _)) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "{:?} leaves unparsed data: {:?}",
            input_str,
            remaining
        );
        assert_eq!(
            name, expected_name,
            "{:?}: name {:?} != expected {:?}",
            input_str, name, expected_name
        );
        assert_eq!(
            value, expected_value,
            "{:?}: value {:?} != expected {:?}",
            input_str, value, expected_value
        );
    }

    #[rstest]
    #[case::terminated(b"$$$$\n")]
    #[case::no_newline(b"$$$$")]
    fn test_sdf_delimiter(#[case] input: &[u8]) {
        let input_str = input.to_str_lossy();
        let result = sdf_delimiter(0).parse(input);
        assert!(
            result.is_ok(),
            "{:?} should have succeeded, result: {:?}",
            input_str,
            result
        );
        let (remaining, _) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "{:?} leaves unparsed data: {:?}",
            input_str,
            remaining
        );
    }

    #[rstest]
    #[case::single_entry(b"> <NAMES>\nbenzene\nBenzol\n\n$$$$\n", indexmap! {"NAMES".to_string() => "benzene,Benzol".to_string()})]
    #[case::two_entries(b"> <BOILING.POINT>\n100.5\n\n> <CAS NR>\n110-82-7\n12217-02-6\n\n$$$$\n",
        indexmap! {"BOILING.POINT".to_string() => "100.5".to_string(), "CAS NR".to_string() => "110-82-7,12217-02-6".to_string()})]
    fn test_sdf_data_block(#[case] input: &[u8], #[case] expected: IndexMap<String, String>) {
        let result = sdf_data_block(0).parse(input);
        let input_str = input.to_str_lossy().into_owned();
        assert!(
            result.is_ok(),
            "{:?} should have succeeded, result: {:?}",
            input_str,
            result
        );
        let (remaining, (data, _)) = result.unwrap();
        assert!(
            remaining.is_empty(),
            "{:?} leaves unparsed data: {:?}",
            input_str,
            remaining
        );
        assert_eq!(
            data, expected,
            "{:?}: data {:?} != expected {:?}",
            input_str, data, expected
        );
    }
}
