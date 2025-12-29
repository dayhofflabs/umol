//! Tests for CTAB block parsing

use rstest::{fixture, rstest};
use umol_data::Element;

use super::*;
use crate::table_ir::{AtomSymbol, BondOrder};

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

#[test]
fn test_atom_block() {
    let atom_data = b"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
";
    let flags = CtabParseFlags::BASIC;
    let result = atom_block(2, 0, flags).parse(atom_data);
    assert!(result.is_ok(), "Atom block should parse successfully");

    let (remaining, (atoms, _, _)) = result.unwrap();
    assert_eq!(remaining, b"", "All input should be consumed");
    assert_eq!(atoms.len(), 2, "Should have 2 atoms");
    assert_eq!(atoms[0].element, Element::C);
    assert_eq!(atoms[1].element, Element::O);
}

#[test]
fn test_extended_atom_block() {
    let atom_data = b"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
";
    let flags = CtabParseFlags::EXTENDED;
    let result = extended_atom_block(2, 0, flags).parse(atom_data);
    assert!(result.is_ok(), "Atom block should parse successfully");

    let (remaining, (atoms, _, _)) = result.unwrap();
    assert_eq!(remaining, b"", "All input should be consumed");
    assert_eq!(atoms.len(), 2, "Should have 2 atoms");
    assert_eq!(atoms[0].symbol, AtomSymbol::Element(Element::C));
    assert_eq!(atoms[1].symbol, AtomSymbol::Element(Element::O));
}

#[test]
fn test_bond_block() {
    let bond_data = b"  1  2  1  0  0  0  0\n  1  3  2  0  0  0  0\n";
    let flags = CtabParseFlags::BASIC;
    let result = bond_block(2, 0, flags).parse(bond_data);
    assert!(result.is_ok(), "Bond block should parse successfully");

    let (remaining, (bonds, _)) = result.unwrap();
    assert_eq!(remaining, b"", "All input should be consumed");
    assert_eq!(bonds.len(), 2, "Should have 2 bonds");
    assert_eq!(bonds[0].2.order, BondOrder::Single);
    assert_eq!(bonds[1].2.order, BondOrder::Double);
}

#[test]
fn test_extended_bond_block() {
    let bond_data = b"  1  2  1  0  0  0  0\n  1  3  2  0  0  0  0\n";
    let flags = CtabParseFlags::EXTENDED;
    let result = extended_bond_block(2, 0, flags).parse(bond_data);
    assert!(result.is_ok(), "Bond block should parse successfully");

    let (remaining, (bonds, _)) = result.unwrap();
    assert_eq!(remaining, b"", "All input should be consumed");
    assert_eq!(bonds.len(), 2, "Should have 2 bonds");
    assert_eq!(bonds[0].2.order, BondOrder::Single);
    assert_eq!(bonds[1].2.order, BondOrder::Double);
}

#[test]
fn test_legacy_atom_list_block() {
    let atom_list_data = b"  1 F    3   9   7   8  ";
    let flags = CtabParseFlags::LENIENT;
    let result = legacy_atom_list_block(1, 0, flags).parse(atom_list_data);
    assert!(
        result.is_ok(),
        "Legacy atom list block should parse successfully"
    );

    let (remaining, (atom_list, _)) = result.unwrap();
    assert_eq!(remaining, b"", "All input should be consumed");
    assert_eq!(atom_list.len(), 1, "Should have 1 atom list");
    assert!(matches!(atom_list[0], PropertyEntries::AtomListEntry(_)));
}

#[test]
fn test_properties_block() {
    let ctab_data = b"M  CHG  1   2  -1\nM  END";
    let flags = CtabParseFlags::BASIC;
    let result = properties_block(0, flags).parse(ctab_data);
    assert!(
        result.is_ok(),
        "Properties block should parse successfully in BASIC mode"
    );

    let (remaining, (property_entries, _)) = result.unwrap();
    assert_eq!(remaining, b"M  END", "All input should be consumed");
    assert_eq!(property_entries.len(), 1, "Should have 1 property entry");
    assert!(matches!(
        property_entries[0],
        PropertyEntries::ChargeEntries(_)
    ));
}

#[test]
fn test_extended_properties_block() {
    let ctab_data = b"M  ALS   1  2 F Cl  Br\nM  END";
    let flags = CtabParseFlags::EXTENDED;
    let result = extended_properties_block(0, flags).parse(ctab_data);
    assert!(
        result.is_ok(),
        "Extended properties block should parse successfully in EXTENDED mode"
    );
}

#[test]
fn test_properties_block_missing_newline() {
    let ctab_data = b"M  END";
    let flags = CtabParseFlags::BASIC;
    let result = properties_block(0, flags).parse(ctab_data);
    assert!(
        result.is_ok(),
        "Properties block should parse with terminating newline in LENIENT mode"
    );

    let (remaining, (property_entries, _)) = result.unwrap();
    assert_eq!(remaining, b"M  END", "Should leave M  END for next parser");
    assert_eq!(property_entries.len(), 0, "Should have 0 property entries");
}

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

    let (remaining, molecule) = result.unwrap();
    assert!(remaining.is_empty(), "All input should be consumed");
    assert_eq!(molecule.atom_count(), 2, "Should have 2 atoms");
    assert_eq!(molecule.bond_count(), 1, "Should have 1 bond");
}

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

    let (_, molecule) = result.unwrap();
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

    let (remaining, molecule) = result.unwrap();
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

    let (remaining, molecule) = result.unwrap();
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

    let (remaining, molecule) = result.unwrap();
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

    let (remaining, molecule) = result.unwrap();
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

    let (remaining, molecule) = result.unwrap();
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

    let (remaining, molecule) = result.unwrap();
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

    let (remaining, molecule) = result.unwrap();
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
