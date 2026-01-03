//! MOL file parsing conformance tests.
//!
//! This module runs parsing tests against a collection of MOL files from various sources.
//! Files are organized by parse result category:
//! - molecule: passes all 4 parsers
//! - molecule_lenient: needs lenient flags for basic parser
//! - extended_molecule: needs extended parser
//! - extended_molecule_lenient: needs lenient flags for extended parser
//! - invalid: fails all parsers

use std::path::{Path, PathBuf};

use insta::{assert_yaml_snapshot, Settings};
use rstest::rstest;
use serde::Serialize;
use umol_models_graph::io::ctfile::config::{CtabParseFlags, CtfileIoConfig};
use umol_models_graph::io::ctfile::parser::{
    parse_extended_mol_bytes_with, parse_mol_bytes_with,
};
use umol_models_graph::table_ir::{ExtendedMolecule, Molecule};

/// Summary of molecule statistics for snapshot comparison
#[derive(Serialize)]
struct MoleculeSummary {
    sum_formula: String,
    atom_count: usize,
    bond_count: usize,
    property_count: usize,
}

/// Extended summary including extended features
#[derive(Serialize)]
struct ExtendedMoleculeSummary {
    sum_formula: String,
    atom_count: usize,
    bond_count: usize,
    extended_atoms: usize,
    extended_bonds: usize,
    property_count: usize,
    rgroup_count: usize,
    sgroup_count: usize,
}

/// Error summary for failed parses
#[derive(Serialize)]
struct ParseErrorSummary {
    error_type: String,
    message: String,
}

/// Result of a single parser invocation
#[derive(Serialize)]
struct ParseResult {
    success: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    summary: Option<MoleculeSummary>,
    #[serde(skip_serializing_if = "Option::is_none")]
    extended_summary: Option<ExtendedMoleculeSummary>,
    #[serde(skip_serializing_if = "Option::is_none")]
    error: Option<ParseErrorSummary>,
}

/// Combined results from all 4 parsers for a single file
#[derive(Serialize)]
struct FileParseResults {
    molecule: ParseResult,
    molecule_lenient: ParseResult,
    extended_molecule: ParseResult,
    extended_molecule_lenient: ParseResult,
}

impl From<&Molecule> for MoleculeSummary {
    fn from(mol: &Molecule) -> Self {
        Self {
            sum_formula: mol.sum_formula(),
            atom_count: mol.atom_count(),
            bond_count: mol.bond_count(),
            property_count: mol.property_count(),
        }
    }
}

impl From<&ExtendedMolecule> for ExtendedMoleculeSummary {
    fn from(mol: &ExtendedMolecule) -> Self {
        Self {
            sum_formula: mol.sum_formula(),
            atom_count: mol.atom_count(),
            bond_count: mol.bond_count(),
            extended_atoms: mol.extended_atom_count(),
            extended_bonds: mol.extended_bond_count(),
            property_count: mol.property_count(),
            rgroup_count: mol.rgroup_count(),
            sgroup_count: mol.sgroup_count(),
        }
    }
}

fn parse_with_basic_flags(bytes: &[u8]) -> ParseResult {
    let config = CtfileIoConfig::with_parse_flags(CtabParseFlags::BASIC);
    match parse_mol_bytes_with(bytes, &config) {
        Ok(mol) => ParseResult {
            success: true,
            summary: Some(MoleculeSummary::from(&mol)),
            extended_summary: None,
            error: None,
        },
        Err(e) => ParseResult {
            success: false,
            summary: None,
            extended_summary: None,
            error: Some(ParseErrorSummary {
                error_type: format!("{:?}", std::mem::discriminant(&e)),
                message: e.to_string(),
            }),
        },
    }
}

fn parse_with_lenient_flags(bytes: &[u8]) -> ParseResult {
    let config = CtfileIoConfig::with_parse_flags(CtabParseFlags::BASIC_MAX);
    match parse_mol_bytes_with(bytes, &config) {
        Ok(mol) => ParseResult {
            success: true,
            summary: Some(MoleculeSummary::from(&mol)),
            extended_summary: None,
            error: None,
        },
        Err(e) => ParseResult {
            success: false,
            summary: None,
            extended_summary: None,
            error: Some(ParseErrorSummary {
                error_type: format!("{:?}", std::mem::discriminant(&e)),
                message: e.to_string(),
            }),
        },
    }
}

fn parse_extended_with_extended_flags(bytes: &[u8]) -> ParseResult {
    let config = CtfileIoConfig::with_parse_flags(CtabParseFlags::EXTENDED);
    match parse_extended_mol_bytes_with(bytes, &config) {
        Ok(mol) => ParseResult {
            success: true,
            summary: None,
            extended_summary: Some(ExtendedMoleculeSummary::from(&mol)),
            error: None,
        },
        Err(e) => ParseResult {
            success: false,
            summary: None,
            extended_summary: None,
            error: Some(ParseErrorSummary {
                error_type: format!("{:?}", std::mem::discriminant(&e)),
                message: e.to_string(),
            }),
        },
    }
}

fn parse_extended_with_lenient_flags(bytes: &[u8]) -> ParseResult {
    let config = CtfileIoConfig::with_parse_flags(CtabParseFlags::LENIENT);
    match parse_extended_mol_bytes_with(bytes, &config) {
        Ok(mol) => ParseResult {
            success: true,
            summary: None,
            extended_summary: Some(ExtendedMoleculeSummary::from(&mol)),
            error: None,
        },
        Err(e) => ParseResult {
            success: false,
            summary: None,
            extended_summary: None,
            error: Some(ParseErrorSummary {
                error_type: format!("{:?}", std::mem::discriminant(&e)),
                message: e.to_string(),
            }),
        },
    }
}

fn parse_file(path: &Path) -> FileParseResults {
    let bytes = std::fs::read(path).expect("Failed to read file");

    FileParseResults {
        molecule: parse_with_basic_flags(&bytes),
        molecule_lenient: parse_with_lenient_flags(&bytes),
        extended_molecule: parse_extended_with_extended_flags(&bytes),
        extended_molecule_lenient: parse_extended_with_lenient_flags(&bytes),
    }
}

fn run_conformance_test(file_path: &PathBuf) {
    let source = file_path
        .parent()
        .unwrap()
        .file_name()
        .unwrap()
        .to_str()
        .unwrap();
    let filename = file_path.file_name().unwrap().to_str().unwrap();

    let results = parse_file(file_path);

    let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("mol_parsing");

    let mut settings = Settings::clone_current();
    settings.set_snapshot_path(base.join("snapshots"));
    settings.set_snapshot_suffix(format!("{}_{}", source, filename));
    settings.bind(|| {
        assert_yaml_snapshot!(results);
    });
}

#[rstest]
fn test_conformance(#[files("tests/mol_parsing/data/**/*.mol")] file_path: PathBuf) {
    run_conformance_test(&file_path);
}
