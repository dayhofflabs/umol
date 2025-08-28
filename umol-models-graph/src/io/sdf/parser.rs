//! SDF (Structure Data File) format parsing

use crate::io::mol::parser::{mol_file_moleculelike, MolFileLike};
use bstr::ByteSlice;
use indexmap::IndexMap;
use nom::bytes::complete::{tag, take_until};
use nom::character::complete::{line_ending, multispace0, not_line_ending};
use nom::combinator::{all_consuming, complete, eof, map, opt, peek, value};
use nom::multi::{many1, many_till};
use nom::sequence::{delimited, terminated};
use nom::{branch::alt, error, Parser};
use serde::{Deserialize, Serialize};
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
    pub mol_file: MolFileLike,
    pub data_fields: IndexMap<String, String>,
}

impl SdfCompound {
    pub fn new(mol_file: MolFileLike, data_fields: IndexMap<String, String>) -> Self {
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
    map(
        many_till(
            terminated(not_line_ending::<&[u8], error::Error<&[u8]>>, line_ending),
            alt((
                peek(line_ending), // blank line
                eof,
            )),
        ),
        |(lines, _)| {
            lines
                .iter()
                .map(|line| line.to_str_lossy())
                .collect::<Vec<_>>()
                .join("\n")
        },
    )
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
) -> impl Parser<&'a [u8], Output = IndexMap<String, String>, Error = error::Error<&'a [u8]>> {
    map(
        many_till(
            alt((
                data_field(),
                // Skip non-data lines
                value(
                    (String::new(), String::new()),
                    terminated(not_line_ending, line_ending),
                ),
            )),
            peek(alt((tag("$$$$"), eof))),
        ),
        |(fields, _)| {
            fields
                .into_iter()
                .filter(|(name, _)| !name.is_empty())
                .collect::<IndexMap<_, _>>()
        },
    )
}

/// Parse single SDF compound (MOL + data + $$$$)
fn sdf_compound<'a>() -> impl Parser<&'a [u8], Output = SdfCompound, Error = error::Error<&'a [u8]>>
{
    move |input: &'a [u8]| {
        let (remaining, mol_input) = alt((
            take_until(">"),
            take_until("$$$$"),
            // If no data fields, take everything
            |input: &'a [u8]| Ok((&[][..], input)),
        ))
        .parse(input)?;

        let (_, mol_file) = mol_file_moleculelike()
            .parse(mol_input)
            .map_err(|_| nom::Err::Error(error::Error::new(input, error::ErrorKind::Verify)))?;

        let (remaining, data_fields) = data_block().parse(remaining)?;
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
