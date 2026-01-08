//! SDF data field parsers

use bstr::{join, ByteSlice};
use indexmap::IndexMap;
use nom::bytes::complete::{is_not, tag, take_until1};
use nom::character::complete::{line_ending, multispace0, not_line_ending, space0};
use nom::combinator::{map, value};
use nom::error::Error as NomError;
use nom::multi::many_till;
use nom::sequence::{delimited, terminated};
use nom::Parser;

/// Parse data field header: `> <Field Name>`
pub(super) fn sdf_data_header<'inp>(
) -> impl Parser<&'inp [u8], Output = String, Error = NomError<&'inp [u8]>> {
    map(
        delimited(
            (tag(">"), is_not("<")),
            delimited(tag("<"), take_until1(">"), tag(">")),
            (not_line_ending, line_ending),
        ),
        |field_name: &[u8]| field_name.trim().to_str_lossy().into_owned(),
    )
}

/// Parse multi-line data value until blank line
pub(super) fn sdf_data_value<'inp>(
) -> impl Parser<&'inp [u8], Output = String, Error = NomError<&'inp [u8]>> {
    map(
        many_till(
            terminated(not_line_ending, line_ending),
            value((), (space0, line_ending)),
        ),
        |(lines, _)| {
            join(",", lines.iter().map(|line: &&[u8]| line.trim()))
                .to_str_lossy()
                .to_string()
        },
    )
}

/// Parse complete data field (header + value)
pub(super) fn sdf_data_field<'inp>(
) -> impl Parser<&'inp [u8], Output = (String, String), Error = NomError<&'inp [u8]>> {
    (sdf_data_header(), sdf_data_value())
}

/// Parse SDF record delimiter
pub(super) fn sdf_delimiter<'inp>(
) -> impl Parser<&'inp [u8], Output = (), Error = NomError<&'inp [u8]>> {
    value((), (tag("$$$$"), multispace0))
}

/// Parse multiple data fields
pub(super) fn sdf_data_block<'inp>(
) -> impl Parser<&'inp [u8], Output = IndexMap<String, String>, Error = NomError<&'inp [u8]>> {
    map(
        many_till(sdf_data_field(), sdf_delimiter()),
        |(fields, _)| fields.into_iter().collect(),
    )
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
        let result = sdf_data_header().parse(input);
        let input_str = input.to_str_lossy();
        assert!(
            result.is_ok(),
            "{:?} should have succeeded, result: {:?}",
            input_str,
            result
        );
        let (remaining, name) = result.unwrap();
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
        let result = sdf_data_value().parse(input);
        let input_str = input.to_str_lossy();
        assert!(
            result.is_ok(),
            "{:?} should have succeeded, result: {:?}",
            input_str,
            result
        );
        let (remaining, value) = result.unwrap();
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
        let result = sdf_data_field().parse(input);
        let input_str = input.to_str_lossy();
        assert!(
            result.is_ok(),
            "{:?} should have succeeded, result: {:?}",
            input_str,
            result
        );
        let (remaining, (name, value)) = result.unwrap();
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
        let result = sdf_delimiter().parse(input);
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
        let result = sdf_data_block().parse(input);
        let input_str = input.to_str_lossy().to_owned();
        assert!(
            result.is_ok(),
            "{:?} should have succeeded, result: {:?}",
            input_str,
            result
        );
        let (remaining, data) = result.unwrap();
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
