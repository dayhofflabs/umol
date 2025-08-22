//! SDF (Structure Data File) format parsing

use bstr::ByteSlice;
use nom::{
    bytes::complete::{tag, take_until},
    character::complete::{line_ending, multispace0, not_line_ending},
    combinator::{all_consuming, complete, map, opt, value},
    error,
    multi::many1,
    sequence::{delimited, terminated},
    Parser,
};
use serde::{Deserialize, Serialize};

use crate::io::mol::parser::{parse_mol_file, MolFile};
use umol::error::DataError;
use umol::Result;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SdfFile {
    pub compounds: Vec<SdfCompound>,
}

impl SdfFile {
    pub fn new(compounds: Vec<SdfCompound>) -> Self {
        Self { compounds }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SdfCompound {
    pub mol_file: MolFile,
    pub data_fields: Vec<(String, String)>, // Preserves insertion order
}

impl SdfCompound {
    pub fn new(mol_file: MolFile, data_fields: Vec<(String, String)>) -> Self {
        Self {
            mol_file,
            data_fields,
        }
    }
}

/// Parse data field header: `> <Field Name>`
fn data_header<'a>() -> impl Parser<&'a [u8], Output = String, Error = error::Error<&'a [u8]>> {
    map(
        delimited(
            (tag(">"), many1(tag(" "))), // Allow multiple spaces after >
            delimited(tag("<"), take_until(">"), tag(">")),
            not_line_ending, // Ignore any trailing content
        ),
        |field_name: &[u8]| field_name.to_str_lossy().into_owned(),
    )
}

/// Parse multi-line data value until blank line
fn data_value<'a>() -> impl Parser<&'a [u8], Output = String, Error = error::Error<&'a [u8]>> {
    move |input: &'a [u8]| {
        let mut lines = Vec::new();
        let mut consumed = 0;

        for line in input.lines_with_terminator() {
            if line.trim_ascii().is_empty() {
                consumed += line.len();
                break; // Blank line terminates data
            }
            lines.push(line.trim_end_with(|c| c == '\r' || c == '\n'));
            consumed += line.len();
        }

        let value = lines
            .iter()
            .map(|l| l.to_str_lossy())
            .collect::<Vec<_>>()
            .join("\n");

        Ok((&input[consumed..], value))
    }
}

/// Parse complete data field (header + value)
fn data_field<'a>(
) -> impl Parser<&'a [u8], Output = (String, String), Error = error::Error<&'a [u8]>> {
    (terminated(data_header(), line_ending), data_value())
}

/// Parse SDF record delimiter
fn sdf_delimiter<'a>() -> impl Parser<&'a [u8], Output = (), Error = error::Error<&'a [u8]>> {
    value((), (tag("$$$$"), opt(line_ending)))
}

/// Parse multiple data fields
fn data_block<'a>(
) -> impl Parser<&'a [u8], Output = Vec<(String, String)>, Error = error::Error<&'a [u8]>> {
    move |input: &'a [u8]| {
        let mut fields = Vec::new();
        let mut remaining = input;

        while !remaining.is_empty() {
            if remaining.starts_with(b"$$$$") {
                break;
            } else if remaining.starts_with(b">") {
                match data_field().parse(remaining) {
                    Ok((next_remaining, (field_name, field_value))) => {
                        fields.push((field_name, field_value));
                        remaining = next_remaining;
                    }
                    Err(_) => {
                        // Skip invalid data field line
                        if let Some(next_line_start) = remaining.lines_with_terminator().next() {
                            remaining = &remaining[next_line_start.len()..];
                        } else {
                            break;
                        }
                    }
                }
            } else {
                // Skip non-data lines
                if let Some(next_line) = remaining.lines_with_terminator().next() {
                    remaining = &remaining[next_line.len()..];
                } else {
                    break;
                }
            }
        }

        Ok((remaining, fields))
    }
}

/// Parse single SDF compound (MOL + data + $$$$)
fn sdf_compound<'a>() -> impl Parser<&'a [u8], Output = SdfCompound, Error = error::Error<&'a [u8]>>
{
    move |input: &'a [u8]| {
        let mut mol_end = input.len();
        let mut offset = 0;

        // Find where MOL data ends and SDF data begins
        for line in input.lines_with_terminator() {
            if line.starts_with(b">") || line.starts_with(b"$$$$") {
                mol_end = offset;
                break;
            }
            offset += line.len();
        }

        // Parse MOL file
        let mol_input = &input[..mol_end];
        let mol_file = match parse_mol_file(mol_input) {
            Ok(mol_file) => mol_file,
            Err(_) => {
                return Err(nom::Err::Error(error::Error::new(
                    input,
                    error::ErrorKind::Verify,
                )))
            }
        };

        // Parse data fields
        let (remaining, data_fields) = data_block().parse(&input[mol_end..])?;

        // Parse delimiter
        let (remaining, _) = sdf_delimiter().parse(remaining)?;

        Ok((remaining, SdfCompound::new(mol_file, data_fields)))
    }
}

/// Parse complete SDF file (multiple compounds)
fn sdf_file<'a>() -> impl Parser<&'a [u8], Output = SdfFile, Error = error::Error<&'a [u8]>> {
    map(many1(sdf_compound()), |compounds| SdfFile::new(compounds))
}

/// Public API function to parse SDF files
pub fn parse_sdf(input: &[u8]) -> Result<SdfFile> {
    all_consuming(complete(terminated(sdf_file(), multispace0)))
        .parse(input)
        .map(|(_, sdf_file)| sdf_file)
        .map_err(|e| DataError::InvalidSdfFormat(format!("SDF parsing failed: {:?}", e)).into())
}

#[cfg(test)]
mod tests;
