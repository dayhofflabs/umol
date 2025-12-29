//! SDF (Structure Data File) format parsing

use bstr::ByteSlice;
use indexmap::IndexMap;
use nom::branch::alt;
use nom::bytes::complete::{tag, take_until};
use nom::character::complete::{line_ending, not_line_ending};
use nom::combinator::{eof, map, opt, peek, value};
use nom::multi::{many1, many_till};
use nom::sequence::{delimited, terminated};
use nom::Parser;
use crate::io::ctab::config::MolIoConfig;
use crate::io::ctab::parser::extended_ctab_block;
use crate::io::ctfile::error::ParseError;
use crate::io::mol::parser::{header, ExtendedMolFile};

#[derive(Debug, Clone)]
pub struct SdfFile {
    pub compounds: Vec<SdfCompound>,
}

impl SdfFile {
    pub fn new(compounds: Vec<SdfCompound>) -> Self {
        Self { compounds }
    }
}

#[derive(Debug, Clone)]
pub struct SdfCompound {
    pub mol_file: ExtendedMolFile,
    pub data_fields: IndexMap<String, String>,
}

impl SdfCompound {
    pub fn new(mol_file: ExtendedMolFile, data_fields: IndexMap<String, String>) -> Self {
        Self {
            mol_file,
            data_fields,
        }
    }
}

/// Parse data field header: `> <Field Name>`
fn data_header<'inp>(
) -> impl Parser<&'inp [u8], Output = String, Error = nom::error::Error<&'inp [u8]>> {
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
fn data_value<'inp>(
) -> impl Parser<&'inp [u8], Output = String, Error = nom::error::Error<&'inp [u8]>> {
    map(
        many_till(
            terminated(not_line_ending::<&[u8], nom::error::Error<&[u8]>>, line_ending),
            alt((
                peek(line_ending), // blank line
                eof,
            )),
        ),
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
fn data_field<'inp>(
) -> impl Parser<&'inp [u8], Output = (String, String), Error = nom::error::Error<&'inp [u8]>> {
    (terminated(data_header(), line_ending), data_value())
}

/// Parse SDF record delimiter
fn sdf_delimiter<'inp>(
) -> impl Parser<&'inp [u8], Output = (), Error = nom::error::Error<&'inp [u8]>> {
    value((), (tag("$$$$"), opt(line_ending)))
}

/// Parse multiple data fields
fn data_block<'inp>(
) -> impl Parser<&'inp [u8], Output = IndexMap<String, String>, Error = nom::error::Error<&'inp [u8]>> {
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
        |(fields, _): (Vec<(String, String)>, _)| {
            fields
                .into_iter()
                .filter(|(name, _)| !name.is_empty())
                .collect::<IndexMap<_, _>>()
        },
    )
}

/// Parse single SDF compound
pub fn parse_sdf_compound<'inp>(
    input: &'inp [u8],
    mut current_line: u32,
) -> Result<(&'inp [u8], SdfCompound), ParseError> {
    let config = MolIoConfig::lenient();
    let flags = config.parse_flags;

    let (remaining, header) = header::header()
        .parse(input)
        .map_err(|e| ParseError::header_from_nom(e, current_line))?;
    current_line += 3;

    let (remaining, molecule) = extended_ctab_block(current_line, &flags).parse(remaining)?;

    // Recalculate current_line after ctab_block
    let consumed_by_mol = input.len() - remaining.len();
    current_line += input[..consumed_by_mol].lines_with_terminator().count() as u32 - 3;

    let (remaining, data_fields) = data_block()
        .parse(remaining)
        .map_err(|e| ParseError::sdf_data_from_nom(e, current_line))?;

    let (remaining, _) = opt(sdf_delimiter())
        .parse(remaining)
        .map_err(|e| ParseError::delimiter_from_nom(e, current_line))?;

    Ok((
        remaining,
        SdfCompound::new(ExtendedMolFile::new(header, molecule), data_fields),
    ))
}

/// Parse SDF from bytes
pub fn parse_sdf_bytes(input: &[u8]) -> Result<SdfFile, ParseError> {
    let mut compounds = Vec::new();
    let mut current_line = 0;
    let mut remaining = input;

    while !remaining.trim_ascii().is_empty() {
        let (rem, compound) = parse_sdf_compound(remaining, current_line)?;
        compounds.push(compound);

        let consumed = remaining.len() - rem.len();
        current_line += remaining[..consumed].lines_with_terminator().count() as u32;
        remaining = rem;
    }

    Ok(SdfFile::new(compounds))
}

/// Parse SDF from string
pub fn parse_sdf(input: &str) -> Result<SdfFile, ParseError> {
    parse_sdf_bytes(input.as_bytes())
}

#[cfg(test)]
mod tests;
