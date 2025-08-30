//! Tests for CTAB block parsing

use super::*;
use crate::io::ctab::atom::AtomSymbol;
use crate::io::ctab::bond::BondType;
use rstest::{fixture, rstest};
use umol_data::Element;

#[fixture]
fn header_atoms_only() -> &'static [u8] {
    b"  2  0  0  0  0  0  0  0  0  0999 V2000\n"
}

#[fixture]
fn header_atoms_bonds() -> &'static [u8] {
    b"  2  1  0  0  0  0  0  0  0  0999 V2000\n"
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
    let atom_count = 2;
    let flags = ParseFlags::LENIENT;
    let result = atom_block(atom_count, flags).parse(atom_data);
    assert!(result.is_ok(), "Atom block should parse successfully");

    let (remaining, atoms) = result.unwrap();
    assert_eq!(remaining, b"", "All input should be consumed");
    assert_eq!(atoms.len(), 2, "Should have 2 atoms");
    assert_eq!(atoms[0].element, Element::C);
    assert_eq!(atoms[1].element, Element::O);
}

#[test]
fn test_atomlike_block() {
    let atom_data = b"    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
";
    let atom_count = 2;
    let flags = ParseFlags::LENIENT;
    let result = atomlike_block(atom_count, flags).parse(atom_data);
    assert!(result.is_ok(), "Atom block should parse successfully");

    let (remaining, atoms) = result.unwrap();
    assert_eq!(remaining, b"", "All input should be consumed");
    assert_eq!(atoms.len(), 2, "Should have 2 atoms");
    assert_eq!(atoms[0].symbol, AtomSymbol::Element(Element::C));
    assert_eq!(atoms[1].symbol, AtomSymbol::Element(Element::O));
}

#[test]
fn test_bond_block() {
    let bond_data = b"  1  2  1  0  0  0  0\n  1  3  2  0  0  0  0\n";
    let bond_count = 2;
    let flags = ParseFlags::LENIENT;
    let result = bond_block(bond_count, flags).parse(bond_data);
    assert!(result.is_ok(), "Bond block should parse successfully");

    let (remaining, bonds) = result.unwrap();
    assert_eq!(remaining, b"", "All input should be consumed");
    assert_eq!(bonds.len(), 2, "Should have 2 bonds");
    assert_eq!(bonds[0].2.bond_type, BondType::Single);
    assert_eq!(bonds[1].2.bond_type, BondType::Double);
}

#[test]
fn test_bondlike_block() {
    let bond_data = b"  1  2  1  0  0  0  0\n  1  3  2  0  0  0  0\n";
    let bond_count = 2;
    let flags = ParseFlags::LENIENT;
    let result = bondlike_block(bond_count, flags).parse(bond_data);
    assert!(result.is_ok(), "Bond block should parse successfully");

    let (remaining, bonds) = result.unwrap();
    assert_eq!(remaining, b"", "All input should be consumed");
    assert_eq!(bonds.len(), 2, "Should have 2 bonds");
    assert_eq!(bonds[0].2.bond_type, BondType::Single);
    assert_eq!(bonds[1].2.bond_type, BondType::Double);
}

#[test]
fn test_legacy_atom_list_block() {
    let atom_list_data = b"  1 F    3   9   7   8  ";
    let flags = ParseFlags::LENIENT;
    let result = legacy_atom_list_block(flags).parse(atom_list_data);
    assert!(
        result.is_ok(),
        "Legacy atom list block should parse successfully"
    );

    let (remaining, atom_list) = result.unwrap();
    assert_eq!(remaining, b"", "All input should be consumed");
    assert_eq!(atom_list.len(), 1, "Should have 1 atom list");
    assert!(matches!(atom_list[0], PropertyEntries::AtomListEntry(_)));
}

#[test]
fn test_basic_properties_block_missing_newline() {
    let ctab_data = b"M  END";
    let flags = ParseFlags::LENIENT;
    let result = basic_properties_block(flags).parse(ctab_data);
    assert!(
        result.is_ok(),
        "CTAB block without terminating newline should parse successfully"
    );

    let (remaining, property_entries) = result.unwrap();
    assert_eq!(
        remaining, b"M  END",
        "Should leave M  END for next parser"
    );
    assert_eq!(property_entries.len(), 0, "Should have 0 property entries");
}

#[test]
fn test_basic_ctab_block_fails_on_query() {
    let ctab_data = b"  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  ALS   1  2 F Cl  Br
M  END
";
    let flags = ParseFlags::BASIC;
    let result = basic_ctab_block(flags).parse(ctab_data);
    assert!(
        result.is_err(),
        "Basic parser should fail on query features"
    );
}

#[test]
fn test_ctab_block_truncated_lines() {
    let ctab_data = b"  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
    1.5400    0.0000    0.0000 C   0  0
  1  2  1
M  END
";
    let flags = ParseFlags::LENIENT;
    let result = ctab_block(flags).parse(ctab_data);
    assert!(
        result.is_ok(),
        "CTAB block with truncated lines should parse successfully"
    );

    let (remaining, molecule) = result.unwrap();
    assert!(remaining.is_empty(), "All input should be consumed");
    assert_eq!(molecule.graph.node_count(), 2, "Should have 2 atoms");
    assert_eq!(molecule.graph.edge_count(), 1, "Should have 1 bond");
}

#[test]
fn test_ctab_block_insufficient_atoms() {
    let ctab_data = b"  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
";
    let flags = ParseFlags::LENIENT;
    let result = ctab_block(flags).parse(ctab_data);
    assert!(result.is_err(), "Should fail with insufficient atoms");
}

#[test]
fn test_ctab_block_insufficient_bonds() {
    let ctab_data = b"  2  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
";
    let flags = ParseFlags::LENIENT;
    let result = ctab_block(flags).parse(ctab_data);
    assert!(result.is_err(), "Should fail with insufficient bonds");
}

#[test]
fn test_ctab_block_query_features() {
    let ctab_data = b"  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  ALS   1  2 F Cl  Br
M  END
";
    let flags = ParseFlags::LENIENT;
    let result = ctab_block(flags).parse(ctab_data);
    assert!(
        result.is_ok(),
        "CTAB block with query features should parse successfully"
    );

    let (_, molecule) = result.unwrap();
    let atom1 = molecule.graph.node_weight(0.into()).unwrap();
    assert!(matches!(atom1.symbol, AtomSymbol::AtomList(_)));
}

#[rstest]
#[case("basic", vec![], m_end())]
#[case("no_m_end", vec![], b"")]
#[case("with_properties", vec![properties()], m_end())]
fn test_ctab_block_termination(
    #[case] name: &str,
    #[case] extra_blocks: Vec<&[u8]>,
    #[case] ending: &[u8],
    header_atoms_bonds: &[u8],
    atoms: &[u8],
    bonds: &[u8],
) {
    let mut data = Vec::new();
    data.extend_from_slice(header_atoms_bonds);
    data.extend_from_slice(atoms);
    data.extend_from_slice(bonds);

    for block in extra_blocks {
        data.extend_from_slice(block);
    }

    data.extend_from_slice(ending);

    let flags = ParseFlags::LENIENT;
    let result = ctab_block(flags).parse(&data);
    assert!(result.is_ok(), "{} case should parse successfully", name);

    let (remaining, molecule) = result.unwrap();
    assert!(remaining.is_empty(), "All input should be consumed");
    assert_eq!(molecule.graph.node_count(), 2, "Should have 2 atoms");
    assert_eq!(molecule.graph.edge_count(), 1, "Should have 1 bond");

    // Check charge property is applied correctly
    if name == "with_properties" {
        let atom2 = molecule.graph.node_weight(1.into()).unwrap();
        assert_eq!(atom2.charge, -1, "Oxygen should have -1 charge from M CHG");
    }
}

#[rstest]
#[case("atoms_only_newline", header_atoms_only(), vec![atoms()], b"")]
#[case("atoms_only_no_newline", header_atoms_only(), vec![atoms()], b"")]
#[case("legacy_list_newline", header_atoms_bonds(), vec![atoms(), bonds(), legacy_atom_list()], b"")]
#[case("legacy_list_no_newline", header_atoms_bonds(), vec![atoms(), bonds(), legacy_atom_list()], b"")]
#[case("properties_with_legacy_newline", header_atoms_bonds(), vec![atoms(), bonds(), legacy_atom_list(), properties()], b"")]
#[case("properties_with_legacy_no_newline", header_atoms_bonds(), vec![atoms(), bonds(), legacy_atom_list(), properties()], b"")]
fn test_ctab_block_missing_m_end(
    #[case] name: &str,
    #[case] header: &[u8],
    #[case] blocks: Vec<&[u8]>,
    #[case] ending: &[u8],
) {
    let mut data = Vec::new();
    data.extend_from_slice(header);

    for block in blocks {
        data.extend_from_slice(block);
    }

    // Remove final newline for "no_newline" cases
    if name.ends_with("no_newline") && data.ends_with(b"\n") {
        data.pop();
    }

    data.extend_from_slice(ending);

    let flags = ParseFlags::LENIENT;
    let result = ctab_block(flags).parse(&data);
    assert!(result.is_ok(), "{} case should parse successfully", name);

    let (remaining, molecule) = result.unwrap();
    assert!(remaining.is_empty(), "All input should be consumed");
    assert_eq!(molecule.graph.node_count(), 2, "Should have 2 atoms");
}

#[rstest]
fn test_ctab_block_m_end_no_newline(
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

    let flags = ParseFlags::LENIENT;
    let result = ctab_block(flags).parse(&data);
    assert!(
        result.is_ok(),
        "M END without newline should parse successfully"
    );

    let (remaining, molecule) = result.unwrap();
    assert!(remaining.is_empty(), "All input should be consumed");
    assert_eq!(molecule.graph.node_count(), 2, "Should have 2 atoms");
    assert_eq!(molecule.graph.edge_count(), 1, "Should have 1 bond");

    let atom2 = molecule.graph.node_weight(1.into()).unwrap();
    assert_eq!(atom2.charge, -1, "Oxygen should have -1 charge from M CHG");
}
