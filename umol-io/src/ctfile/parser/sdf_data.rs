//! SDF data field parsers

use bstr::{join, ByteSlice};
use indexmap::IndexMap;
use winnow::error::ErrMode;
use winnow::ModalResult;

use super::utils::next_line;
use crate::ctfile::error::ParseError;

/// Parse data field header: `> <Field Name>`
fn sdf_data_header(input: &mut &[u8], line_offset: u32) -> ModalResult<(String, u32), ParseError> {
    let line = next_line(input).map_err(|_| {
        ErrMode::Cut(ParseError::UnexpectedEof {
            line: line_offset,
            block: "sdf data",
        })
    })?;
    let bytes: &[u8] = line.as_ref();
    if bytes.first() != Some(&b'>') {
        return Err(ErrMode::Cut(ParseError::InvalidSdfDataHeader {
            line: line_offset,
            col: 0,
        }));
    }

    let Some(open) = bytes[1..]
        .iter()
        .position(|byte| *byte == b'<')
        .map(|position| position + 1)
    else {
        return Err(ErrMode::Cut(ParseError::InvalidSdfDataHeader {
            line: line_offset,
            col: bytes.len() as u32,
        }));
    };
    if open == 1 {
        return Err(ErrMode::Cut(ParseError::InvalidSdfDataHeader {
            line: line_offset,
            col: 1,
        }));
    }

    let name_start = open + 1;
    let Some(close) = bytes[name_start..]
        .iter()
        .position(|byte| *byte == b'>')
        .map(|position| position + name_start)
    else {
        return Err(ErrMode::Cut(ParseError::InvalidSdfDataHeader {
            line: line_offset,
            col: bytes.len() as u32,
        }));
    };
    if close == name_start {
        return Err(ErrMode::Cut(ParseError::InvalidSdfDataHeader {
            line: line_offset,
            col: name_start as u32,
        }));
    }

    let name = bytes[name_start..close]
        .trim_ascii()
        .to_str_lossy()
        .into_owned();
    Ok((name, line_offset + 1))
}

/// Parse multi-line data value until blank line
fn sdf_data_value(input: &mut &[u8], line_offset: u32) -> ModalResult<(String, u32), ParseError> {
    let mut line_index = 0;
    let mut value_lines = Vec::new();

    while !input.is_empty() {
        let line = next_line(input).expect("non-empty input contains a physical line");
        line_index += 1;
        let bytes: &[u8] = line.as_ref();
        if bytes.trim_ascii().is_empty() {
            break;
        }

        value_lines.push(bytes.trim_ascii());
    }

    let value = join(",", value_lines).to_str_lossy().into_owned();
    Ok((value, line_offset + line_index))
}

/// Parse complete data field (header + value)
fn sdf_data_field(
    input: &mut &[u8],
    line_offset: u32,
) -> ModalResult<((String, String), u32), ParseError> {
    let (name, line_offset) = sdf_data_header(input, line_offset)?;
    let (data, line_offset) = sdf_data_value(input, line_offset)?;
    Ok(((name, data), line_offset))
}

/// Parse SDF record delimiter
fn sdf_delimiter(input: &mut &[u8], line_offset: u32) -> ModalResult<u32, ParseError> {
    let line = next_line(input)
        .map_err(|_| ErrMode::Cut(ParseError::MissingDelimiter { line: line_offset }))?;
    let bytes: &[u8] = line.as_ref();
    if !bytes.starts_with(b"$$$$") || !bytes[4..].trim_ascii().is_empty() {
        return Err(ErrMode::Cut(ParseError::MissingDelimiter {
            line: line_offset,
        }));
    }

    let whitespace_len = input
        .iter()
        .take_while(|byte| byte.is_ascii_whitespace())
        .count();
    let skipped = &input[..whitespace_len];
    let skipped_lines = skipped.iter().filter(|byte| **byte == b'\n').count() as u32;
    *input = &input[whitespace_len..];

    Ok(line_offset + 1 + skipped_lines)
}

/// Parse multiple data fields
pub(super) fn sdf_data_block(
    input: &mut &[u8],
    line_offset: u32,
) -> ModalResult<(IndexMap<String, String>, u32), ParseError> {
    let mut remaining_offset = line_offset;
    let mut data = IndexMap::new();

    loop {
        if input.is_empty() {
            return Err(ErrMode::Cut(ParseError::MissingDelimiter {
                line: remaining_offset,
            }));
        }

        if input.starts_with(b">") {
            let ((name, value), new_line_offset) = sdf_data_field(input, remaining_offset)?;
            data.insert(name, value);
            remaining_offset = new_line_offset;
        } else {
            let new_line_offset = sdf_delimiter(input, remaining_offset)?;
            return Ok((data, new_line_offset));
        }
    }
}

#[cfg(test)]
mod tests {
    use indexmap::indexmap;
    use pretty_assertions::assert_eq;
    use rstest::rstest;
    use winnow::error::ErrMode;

    use super::*;

    #[rstest]
    #[case::dotted(b"> <MELTING.POINT>\n", "MELTING.POINT".to_string())]
    #[case::whitespace(b"> <CAS NR>\n", "CAS NR".to_string())]
    #[case::multiple_space(b">  <BOILING.POINT>\n", "BOILING.POINT".to_string())]
    #[case::interstitial_data(b"> (MD-0894) <CAS NR>\n", "CAS NR".to_string())]
    #[case::trailing_data(b"> <CAS NR> DT12\n", "CAS NR".to_string())]
    #[case::surrounding_data(b"> (MD-0894) <BOILING.POINT> FROM ARCHIVES\n", "BOILING.POINT".to_string())]
    #[case::crlf(b"> <NAME>\r\n", "NAME".to_string())]
    fn test_sdf_data_header(#[case] input: &[u8], #[case] expected: String) {
        let mut remaining = input;
        assert_eq!(sdf_data_header(&mut remaining, 7), Ok((expected, 8)));
        assert!(remaining.is_empty());
    }

    #[rstest]
    #[case::missing_open(b"> NAME\n", 6)]
    #[case::missing_prefix_space(b"><NAME>\n", 1)]
    #[case::empty_name(b"> <>\n", 3)]
    #[case::missing_close(b"> <NAME\n", 7)]
    fn test_sdf_data_header_error(#[case] input: &[u8], #[case] col: u32) {
        let mut remaining = input;
        assert_eq!(
            sdf_data_header(&mut remaining, 7),
            Err(ErrMode::Cut(ParseError::InvalidSdfDataHeader {
                line: 7,
                col,
            }))
        );
    }

    #[rstest]
    #[case::single_line(b"100.5\n\n", "100.5".to_string())]
    #[case::whitespace(b" 100.5 \n\n", "100.5".to_string())]
    #[case::multiple_lines(b"benzene\nBenzol\n\n", "benzene,Benzol".to_string())]
    #[case::crlf(b"benzene\r\nBenzol\r\n\r\n", "benzene,Benzol".to_string())]
    fn test_sdf_data_value(#[case] input: &[u8], #[case] expected: String) {
        let line_count = input.iter().filter(|byte| **byte == b'\n').count() as u32;
        let mut remaining = input;
        assert_eq!(
            sdf_data_value(&mut remaining, 4),
            Ok((expected, 4 + line_count))
        );
        assert!(remaining.is_empty());
    }

    #[rstest]
    #[case::single_line(b"> <BOILING.POINT>\n100.5\n\n", "BOILING.POINT".to_string(), "100.5".to_string())]
    #[case::multiple_line(b"> <NAMES>\nbenzene\nBenzol\n\n", "NAMES".to_string(), "benzene,Benzol".to_string())]
    fn test_sdf_data_field(
        #[case] input: &[u8],
        #[case] expected_name: String,
        #[case] expected_value: String,
    ) {
        let line_count = input.iter().filter(|byte| **byte == b'\n').count() as u32;
        let mut remaining = input;
        assert_eq!(
            sdf_data_field(&mut remaining, 3),
            Ok(((expected_name, expected_value), 3 + line_count))
        );
        assert!(remaining.is_empty());
    }

    #[rstest]
    #[case::terminated(b"$$$$\n", b"", 5)]
    #[case::no_newline(b"$$$$", b"", 5)]
    #[case::trailing_space(b"$$$$  \n", b"", 5)]
    #[case::inter_record_whitespace(b"$$$$\r\n\nNext", b"Next", 6)]
    fn test_sdf_delimiter(
        #[case] input: &[u8],
        #[case] expected_remaining: &[u8],
        #[case] expected_offset: u32,
    ) {
        let mut remaining = input;
        assert_eq!(sdf_delimiter(&mut remaining, 4), Ok(expected_offset));
        assert_eq!(remaining, expected_remaining);
    }

    #[rstest]
    #[case::short(b"$$$\n")]
    #[case::trailing_data(b"$$$$x\n")]
    #[case::other_record(b"next\n")]
    fn test_sdf_delimiter_error(#[case] input: &[u8]) {
        let mut remaining = input;
        assert_eq!(
            sdf_delimiter(&mut remaining, 4),
            Err(ErrMode::Cut(ParseError::MissingDelimiter { line: 4 }))
        );
    }

    #[rstest]
    #[case::empty(b"$$$$\n", indexmap! {}, 1)]
    #[case::single_entry(b"> <NAMES>\nbenzene\nBenzol\n\n$$$$\n", indexmap! {"NAMES".to_string() => "benzene,Benzol".to_string()}, 5)]
    #[case::two_entries(b"> <BOILING.POINT>\n100.5\n\n> <CAS NR>\n110-82-7\n12217-02-6\n\n$$$$\n",
        indexmap! {"BOILING.POINT".to_string() => "100.5".to_string(), "CAS NR".to_string() => "110-82-7,12217-02-6".to_string()}, 8)]
    fn test_sdf_data_block(
        #[case] input: &[u8],
        #[case] expected: IndexMap<String, String>,
        #[case] expected_offset: u32,
    ) {
        let mut remaining = input;
        assert_eq!(
            sdf_data_block(&mut remaining, 0),
            Ok((expected, expected_offset))
        );
        assert!(remaining.is_empty());
    }

    #[rstest]
    #[case::missing_delimiter(b"", ParseError::MissingDelimiter { line: 5 })]
    #[case::missing_delimiter_after_field(
        b"> <NAME>\nvalue\n\n",
        ParseError::MissingDelimiter { line: 8 },
    )]
    #[case::malformed_header(
        b"> NAME\n",
        ParseError::InvalidSdfDataHeader { line: 5, col: 6 },
    )]
    #[case::malformed_delimiter(b"$$$\n", ParseError::MissingDelimiter { line: 5 })]
    fn test_sdf_data_block_error(#[case] input: &[u8], #[case] expected: ParseError) {
        let mut remaining = input;
        assert_eq!(
            sdf_data_block(&mut remaining, 5),
            Err(ErrMode::Cut(expected))
        );
    }
}
