// umol-models-graph/src/io/mol.rs
//
// This file will contain the top-level parser for the MOL V2000 format.
// The design follows a line-oriented state machine approach to be robust
// against common format variations found in real-world files.

use crate::atom::AtomStandard;
use crate::bond::BondStandard;
use crate::conformer::Point3D;
use crate::io::ctab::{
    atom, bond,
    counts::{self, Counts},
    properties::{self, ChargeEntry, IsotopeEntry, PropertyEntries, RadicalEntry},
};
use nom::Parser;
use std::io::BufRead;

// --- Final Data Structures ---

/// Represents a fully parsed MOL file.
#[derive(Debug, Default)]
pub struct MolFile {
    header: Vec<String>,
    counts: Counts,
    atoms: Vec<(AtomStandard, Point3D)>,
    bonds: Vec<BondStandard>,
    charges: Vec<ChargeEntry>,
    radicals: Vec<RadicalEntry>,
    isotopes: Vec<IsotopeEntry>,
}

impl MolFile {
    pub(crate) fn header(&self) -> &[String] {
        &self.header
    }

    pub(crate) fn counts(&self) -> &Counts {
        &self.counts
    }

    pub(crate) fn atoms(&self) -> &[(AtomStandard, Point3D)] {
        &self.atoms
    }

    pub(crate) fn bonds(&self) -> &[BondStandard] {
        &self.bonds
    }

    pub(crate) fn charges(&self) -> &[ChargeEntry] {
        &self.charges
    }

    pub(crate) fn radicals(&self) -> &[RadicalEntry] {
        &self.radicals
    }

    pub(crate) fn isotopes(&self) -> &[IsotopeEntry] {
        &self.isotopes
    }
}

#[derive(Debug, PartialEq)]
enum MolParserState {
    ParseHeader,
    ParseCounts,
    ParseAtomBlock { remaining: usize },
    ParseBondBlock { remaining: usize },
    ParseProperties,
    Done,
}

#[derive(Debug, Default)]
struct MolFileBuilder {
    header_lines: Vec<String>,
    counts_line: Option<Counts>,
    atom_lines: Vec<(AtomStandard, Point3D)>,
    bond_lines: Vec<BondStandard>,
    charges: Vec<ChargeEntry>,
    radicals: Vec<RadicalEntry>,
    isotopes: Vec<IsotopeEntry>,
}

impl MolFileBuilder {
    fn build(self) -> Result<MolFile, MolFileError> {
        Ok(MolFile {
            header: self.header_lines,
            counts: self.counts_line.unwrap(),
            atoms: self.atom_lines,
            bonds: self.bond_lines,
            charges: self.charges,
            radicals: self.radicals,
            isotopes: self.isotopes,
        })
    }
}

#[derive(Debug)]
pub enum MolFileError {
    Io(std::io::Error),
    ParseError {
        line_number: usize,
        line: String,
        message: String,
    },
    Incomplete,
    Invalid(String),
}

impl From<std::io::Error> for MolFileError {
    fn from(err: std::io::Error) -> Self {
        MolFileError::Io(err)
    }
}

/// Parse MOL file from streaming source.
pub fn parse_mol_stream<R: BufRead>(mut source: R) -> Result<MolFile, MolFileError> {
    let mut state = MolParserState::ParseHeader;
    let mut builder = MolFileBuilder::default();
    let mut buffer = String::new();
    let mut line_number = 0;

    while source.read_line(&mut buffer).map_err(MolFileError::Io)? > 0 {
        line_number += 1;
        let line = buffer.trim_end();
        println!("L{}: {}", line_number, line);

        if line == "M  END" {
            state = MolParserState::Done;
            break;
        }

        println!("about to process L{}", line_number);
        process_line(&mut state, &mut builder, line, line_number)?;
        println!("processed L{}", line_number);

        buffer.clear();
    }

    if !matches!(
        state,
        MolParserState::ParseProperties | MolParserState::Done
    ) {
        return Err(MolFileError::Incomplete);
    }

    builder.build()
}

/// Process single line of input, updating parser state and builder.
fn process_line(
    state: &mut MolParserState,
    builder: &mut MolFileBuilder,
    line: &str,
    line_number: usize,
) -> Result<(), MolFileError> {
    // Helper to create a contextualized parse error
    let to_parse_error = |e: nom::Err<nom::error::Error<&[u8]>>| MolFileError::ParseError {
        line_number,
        line: line.to_string(),
        message: e.to_string(),
    };

    println!("processing L{}: {} [state: {:?}]", line_number, line, state);
    match *state {
        MolParserState::ParseHeader => {
            builder.header_lines.push(line.to_string());
            if builder.header_lines.len() >= 3 {
                *state = MolParserState::ParseCounts;
            }
        }
        MolParserState::ParseCounts => {
            let counts = counts::counts_input()
                .parse(line.as_bytes())
                .map_err(to_parse_error)?
                .1;

            let atom_count = counts.atoms() as usize;
            let bond_count = counts.bonds() as usize;

            builder.atom_lines = Vec::with_capacity(atom_count);
            builder.bond_lines = Vec::with_capacity(bond_count);
            builder.counts_line = Some(counts);
            *state = MolParserState::ParseAtomBlock {
                remaining: atom_count,
            };
        }
        MolParserState::ParseAtomBlock { ref mut remaining } => {
            if *remaining > 0 {
                let atom = atom::atom_input_standard()
                    .parse(line.as_bytes())
                    .map_err(to_parse_error)?
                    .1;
                builder.atom_lines.push(atom);
                *remaining -= 1;
            }
            if *remaining == 0 {
                *state = MolParserState::ParseBondBlock {
                    remaining: builder.bond_lines.capacity(),
                };
            }
        }
        MolParserState::ParseBondBlock { ref mut remaining } => {
            if *remaining > 0 {
                let bond = bond::bond_input_standard()
                    .parse(line.as_bytes())
                    .map_err(to_parse_error)?
                    .1
                     .2;
                builder.bond_lines.push(bond);
                *remaining -= 1;
            }
            if *remaining == 0 {
                *state = MolParserState::ParseProperties;
            }
        }
        MolParserState::ParseProperties => {
            if let Ok((_, prop)) = properties::property_input_standard(line.as_bytes()) {
                match prop {
                    PropertyEntries::ChargeEntries(entries) => builder.charges.extend(entries),
                    PropertyEntries::RadicalEntries(entries) => builder.radicals.extend(entries),
                    PropertyEntries::IsotopeEntries(entries) => builder.isotopes.extend(entries),
                }
            }
        }
        MolParserState::Done => {}
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_simple_mol_file() {
        let mol_data = [
            "benzene",
            "  -ISIS-  07031514222D",
            "",
            "  6  6  0  0  0  0  0  0  0  0999 V2000",
            "    2.2282   -0.5133    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "    3.0949    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "    3.0949    1.0267    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "    2.2282    1.5400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "    1.3616    1.0267    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "    1.3616    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "  1  2  2  0  0  0  0",
            "  2  3  1  0  0  0  0",
            "  3  4  2  0  0  0  0",
            "  4  5  1  0  0  0  0",
            "  5  6  2  0  0  0  0",
            "  6  1  1  0  0  0  0",
            "M  END",
        ]
        .join("\n");

        let cursor = Cursor::new(mol_data);
        let result = parse_mol_stream(cursor);

        assert!(result.is_ok());
        let mol_file = result.unwrap();

        assert_eq!(mol_file.header().len(), 3);
        assert_eq!(mol_file.counts().atoms(), 6);
        assert_eq!(mol_file.counts().bonds(), 6);
        assert_eq!(mol_file.atoms().len(), 6);
        assert_eq!(mol_file.bonds().len(), 6);
        assert!(mol_file.charges().is_empty());
    }

    #[test]
    fn test_parse_with_properties() {
        let mol_data = [
            "prop_test",
            "  -UMOL-",
            "",
            "  2  1  0  0  0  0  0  0  0  0999 V2000",
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0",
            "    1.5000    0.0000    0.0000 O   0  0  0  0  0  0",
            "  1  2  1  0",
            "M  CHG  2   2  -1   1   1",
            "M  ISO  1   1  13",
            "M  END",
        ]
        .join("\n");

        let cursor = Cursor::new(mol_data);
        let result = parse_mol_stream(cursor);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            MolFileError::ParseError { line_number: 5, .. }
        ));
    }

    #[test]
    fn test_missing_m_end() {
        let mol_data = [
            "benzene",
            "  -ISIS-  07031514222D",
            "",
            "  6  1  0  0  0  0  0  0  0  0999 V2000",
            "    2.2282   -0.5133    0.0000 C   0  0  0",
            "    3.0949    0.0000    0.0000 C   0  0  0",
            "    3.0949    1.0267    0.0000 C   0  0  0",
            "    2.2282    1.5400    0.0000 C   0  0  0",
            "    1.3616    1.0267    0.0000 C   0  0  0",
            "    1.3616    0.0000    0.0000 C   0  0  0",
            "  1  2  1",
        ]
        .join("\n");

        let cursor = Cursor::new(mol_data);
        let result = parse_mol_stream(cursor);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            MolFileError::ParseError { line_number: 5, .. }
        ));
    }

    #[test]
    fn test_malformed_atom_line() {
        let mol_data = [
            "malformed atom",
            "  -UMOL-",
            "",
            "  2  1  0  0  0  0  0  0  0  0999 V2000",
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0", // Malformed
            "    1.5000    0.0000    0.0000 Xx  0  0  0  0  0  0",
            "  1  2  1  0",
            "M  END",
        ]
        .join("\n");

        let cursor = Cursor::new(mol_data);
        let result = parse_mol_stream(cursor);
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            MolFileError::ParseError { line_number: 5, .. }
        ));
    }

    #[test]
    fn test_incomplete() {
        let mol_data = ["header1", "header2", "header3"].join("\n");
        let cursor = Cursor::new(mol_data);
        let result = parse_mol_stream(cursor);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), MolFileError::Incomplete));
    }
}
