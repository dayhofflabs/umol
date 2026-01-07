//! Tests for CTAB block parsing

use nom::Parser;
use pretty_assertions::assert_eq;
use rstest::*;

use crate::io::ctfile::config::CtabParseFlags;
use crate::io::ctfile::parser::{ctab_block, extended_ctab_block};
use crate::table_ir::AtomSymbol;

#[fixture]
fn header_atoms_only() -> &'static [u8] {
    b"  2  0  0  0  0  0  0  0  0  0999 V2000\n"
}

#[fixture]
fn header_atoms_bonds() -> &'static [u8] {
    b"  2  1  0  0  0  0  0  0  0  0999 V2000\n"
}

#[fixture]
fn header_atoms_bonds_legacy() -> &'static [u8] {
    b"  2  1  1  0  0  0  0  0  0  0999 V2000\n"
}

#[fixture]
fn atoms() -> &'static [u8] {
    b"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5400    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
}

#[fixture]
fn bonds() -> &'static [u8] {
    b"  1  2  1  0  0  0  0\n"
}

#[fixture]
fn legacy_atom_list() -> &'static [u8] {
    b"  1 F    3   9   7   8  \n"
}

#[fixture]
fn properties() -> &'static [u8] {
    b"M  CHG  1   2  -1\n"
}

#[fixture]
fn m_end() -> &'static [u8] {
    b"M  END\n"
}

#[fixture]
fn m_end_no_newline() -> &'static [u8] {
    b"M  END"
}

// FIX: Use fixtures
#[test]
fn test_ctab_block_error_extended_features() {
    let ctab_data = b"  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  ALS   1  2 F Cl  Br
M  END
";
    let flags = CtabParseFlags::BASIC;
    let result = ctab_block(0, flags).parse(ctab_data);
    assert!(
        result.is_err(),
        "ctab_block parser should fail on query features"
    );
}

// FIX: Use fixtures
#[test]
fn test_extended_ctab_block_truncated_lines() {
    let ctab_data = b"  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
    1.5400    0.0000    0.0000 C   0  0
  1  2  1
M  END
";
    let flags = CtabParseFlags::LENIENT;
    let result = extended_ctab_block(0, flags).parse(ctab_data);
    assert!(
        result.is_ok(),
        "CTAB block with truncated lines should parse successfully"
    );

    let (remaining, (molecule, _)) = result.unwrap();
    assert!(remaining.is_empty(), "All input should be consumed");
    assert_eq!(molecule.atom_count(), 2, "Should have 2 atoms");
    assert_eq!(molecule.bond_count(), 1, "Should have 1 bond");
}

// FIX: Use fixtures
#[test]
fn test_extended_ctab_block_insufficient_atoms() {
    let ctab_data = b"  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
";
    let flags = CtabParseFlags::EXTENDED;
    let result = extended_ctab_block(0, flags).parse(ctab_data);
    assert!(result.is_err(), "Should fail with insufficient atoms");
}

// FIX: Use fixtures
#[test]
fn test_extended_ctab_block_insufficient_bonds() {
    let ctab_data = b"  2  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
";
    let flags = CtabParseFlags::EXTENDED;
    let result = extended_ctab_block(0, flags).parse(ctab_data);
    assert!(result.is_err(), "Should fail with insufficient bonds");
}

// FIX: Use fixtures
#[test]
fn test_extended_ctab_block() {
    let ctab_data = b"  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  ALS   1  2 F Cl  Br
M  END
";
    let flags = CtabParseFlags::EXTENDED;
    let result = extended_ctab_block(0, flags).parse(ctab_data);
    assert!(
        result.is_ok(),
        "CTAB block with query features should parse successfully"
    );

    let (_, (molecule, _)) = result.unwrap();
    let atom1 = &molecule.atoms[0];
    assert!(matches!(atom1.symbol, AtomSymbol::AtomList(_)));
}

#[rstest]
fn test_extended_ctab_block_termination(
    header_atoms_bonds: &[u8],
    atoms: &[u8],
    bonds: &[u8],
    m_end: &[u8],
) {
    let mut data = Vec::new();
    data.extend_from_slice(header_atoms_bonds);
    data.extend_from_slice(atoms);
    data.extend_from_slice(bonds);
    data.extend_from_slice(m_end);

    let flags = CtabParseFlags::EXTENDED;
    let result = extended_ctab_block(0, flags).parse(&data);
    assert!(result.is_ok(), "Should parse successfully");

    let (remaining, (molecule, _)) = result.unwrap();
    assert!(remaining.is_empty(), "All input should be consumed");
    assert_eq!(molecule.atom_count(), 2, "Should have 2 atoms");
    assert_eq!(molecule.bond_count(), 1, "Should have 1 bond");
}

#[rstest]
fn test_extended_ctab_block_termination_with_properties(
    header_atoms_bonds: &[u8],
    atoms: &[u8],
    bonds: &[u8],
    properties: &[u8],
    m_end: &[u8],
) {
    let mut data = Vec::new();
    data.extend_from_slice(header_atoms_bonds);
    data.extend_from_slice(atoms);
    data.extend_from_slice(bonds);
    data.extend_from_slice(properties);
    data.extend_from_slice(m_end);

    let flags = CtabParseFlags::EXTENDED;
    let result = extended_ctab_block(0, flags).parse(&data);
    assert!(result.is_ok(), "Should parse successfully");

    let (remaining, (molecule, _)) = result.unwrap();
    assert!(remaining.is_empty(), "All input should be consumed");
    assert_eq!(molecule.atom_count(), 2, "Should have 2 atoms");
    assert_eq!(molecule.bond_count(), 1, "Should have 1 bond");

    let atom2 = &molecule.atoms[1];
    assert_eq!(
        atom2.charge,
        Some(-1),
        "Oxygen should have -1 charge from M CHG"
    );
}

#[rstest]
fn test_extended_ctab_block_missing_m_end(header_atoms_only: &[u8], atoms: &[u8]) {
    let mut data = Vec::new();
    data.extend_from_slice(header_atoms_only);
    data.extend_from_slice(atoms);

    let flags = CtabParseFlags::LENIENT;
    let result = extended_ctab_block(0, flags).parse(&data);
    assert!(result.is_ok(), "Should parse successfully");

    let (remaining, (molecule, _)) = result.unwrap();
    assert!(remaining.is_empty(), "All input should be consumed");
    assert_eq!(molecule.atom_count(), 2, "Should have 2 atoms");
}

#[rstest]
fn test_extended_ctab_block_missing_m_end_no_newline(header_atoms_only: &[u8], atoms: &[u8]) {
    let mut data = Vec::new();
    data.extend_from_slice(header_atoms_only);
    data.extend_from_slice(atoms);

    // Remove trailing newline
    if data.ends_with(b"\n") {
        data.pop();
    }

    let flags = CtabParseFlags::LENIENT;
    let result = extended_ctab_block(0, flags).parse(&data);
    assert!(result.is_ok(), "Should parse successfully");

    let (remaining, (molecule, _)) = result.unwrap();
    assert!(remaining.is_empty(), "All input should be consumed");
    assert_eq!(molecule.atom_count(), 2, "Should have 2 atoms");
}

#[rstest]
#[case::legacy_list(vec![atoms(), bonds(), legacy_atom_list()])]
#[case::properties_with_legacy(vec![atoms(), bonds(), legacy_atom_list(), properties()])]
fn test_extended_ctab_block_missing_m_end_legacy(
    #[case] blocks: Vec<&[u8]>,
    header_atoms_bonds_legacy: &[u8],
) {
    let mut data = Vec::new();
    data.extend_from_slice(header_atoms_bonds_legacy);

    for block in blocks {
        data.extend_from_slice(block);
    }

    let result = extended_ctab_block(0, CtabParseFlags::LENIENT).parse(&data);
    assert!(result.is_ok(), "Should parse successfully");

    let (remaining, (molecule, _)) = result.unwrap();
    assert!(remaining.is_empty(), "All input should be consumed");
    assert_eq!(molecule.atom_count(), 2, "Should have 2 atoms");
}

#[rstest]
#[case::legacy_list(vec![atoms(), bonds(), legacy_atom_list()])]
#[case::properties_with_legacy(vec![atoms(), bonds(), legacy_atom_list(), properties()])]
fn test_extended_ctab_block_missing_m_end_legacy_no_newline(
    #[case] blocks: Vec<&[u8]>,
    header_atoms_bonds_legacy: &[u8],
) {
    let mut data = Vec::new();
    data.extend_from_slice(header_atoms_bonds_legacy);

    for block in blocks {
        data.extend_from_slice(block);
    }

    // Remove trailing newline
    if data.ends_with(b"\n") {
        data.pop();
    }

    let result = extended_ctab_block(0, CtabParseFlags::LENIENT).parse(&data);
    assert!(result.is_ok(), "Should parse successfully");

    let (remaining, (molecule, _)) = result.unwrap();
    assert!(remaining.is_empty(), "All input should be consumed");
    assert_eq!(molecule.atom_count(), 2, "Should have 2 atoms");
}

#[rstest]
fn test_extended_ctab_block_m_end_no_newline(
    header_atoms_bonds: &[u8],
    atoms: &[u8],
    bonds: &[u8],
    properties: &[u8],
    m_end_no_newline: &[u8],
) {
    let mut data = Vec::new();
    data.extend_from_slice(header_atoms_bonds);
    data.extend_from_slice(atoms);
    data.extend_from_slice(bonds);
    data.extend_from_slice(properties);
    data.extend_from_slice(m_end_no_newline);

    let flags = CtabParseFlags::EXTENDED;
    let result = extended_ctab_block(0, flags).parse(&data);
    assert!(
        result.is_ok(),
        "M END without newline should parse successfully"
    );

    let (remaining, (molecule, _)) = result.unwrap();
    assert!(remaining.is_empty(), "All input should be consumed");
    assert_eq!(molecule.atom_count(), 2, "Should have 2 atoms");
    assert_eq!(molecule.bond_count(), 1, "Should have 1 bond");

    let atom2 = &molecule.atoms[1];
    assert_eq!(
        atom2.charge,
        Some(-1),
        "Oxygen should have -1 charge from M CHG"
    );
}
