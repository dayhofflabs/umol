//! Tests for CTAB block parsing

use super::*;
use crate::io::ctab::atom::AtomSymbol;
use crate::io::ctab::bond::BondType;
use umol_data::Element;

#[test]
fn test_properties_block() {
    let input = b"M  CHG  1   2  -1\nM  END\n";
    let result = properties_block().parse(input);
    assert!(result.is_ok(), "Should parse properties block");
    let (remaining, properties) = result.unwrap();
    assert_eq!(properties.len(), 1, "Should have 1 property");
    assert!(remaining.is_empty(), "All input should be consumed");
}

#[test]
fn test_bond_block() {
    // Test the bond block parser directly
    let bond_data = b"  1  2  1  0  0  0  0\n  1  3  2  0  0  0  0\n";
    let bond_count = 2;
    let result = bond_block(bond_count).parse(bond_data);
    assert!(result.is_ok(), "Bond block should parse successfully");

    let (remaining, bonds) = result.unwrap();
    assert_eq!(bonds.len(), 2, "Should have 2 bonds");
    assert_eq!(remaining, b"", "All input should be consumed");
}

#[test]
fn test_ctab_block() {
    let ctab_data = b"  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  END
";
    let result = ctab_block().parse(ctab_data);
    assert!(result.is_ok(), "CTAB block should parse successfully");

    let (remaining, molecule) = result.unwrap();
    assert!(remaining.is_empty(), "All input should be consumed");

    // Check molecule structure
    assert_eq!(molecule.graph.node_count(), 2, "Should have 2 atoms");
    assert_eq!(molecule.graph.edge_count(), 1, "Should have 1 bond");

    // Check atoms
    let atom1 = molecule.graph.node_weight(0.into()).unwrap();
    assert!(matches!(atom1.symbol, AtomSymbol::Element(Element::C)));

    let atom2 = molecule.graph.node_weight(1.into()).unwrap();
    assert!(matches!(atom2.symbol, AtomSymbol::Element(Element::C)));

    // Check bond
    let edge = molecule.graph.edge_indices().next().unwrap();
    let bond = molecule.graph.edge_weight(edge).unwrap();
    assert_eq!(bond.bond_type, BondType::Single);
}

#[test]
fn test_ctab_block_with_properties() {
    let ctab_data = b"  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  CHG  1   2  -1
M  END
";
    let result = ctab_block().parse(ctab_data);
    assert!(
        result.is_ok(),
        "CTAB block with properties should parse successfully"
    );

    let (remaining, molecule) = result.unwrap();
    assert!(remaining.is_empty(), "All input should be consumed");

    // Check that we have the expected structure
    assert_eq!(molecule.graph.node_count(), 2, "Should have 2 atoms");
    assert_eq!(molecule.graph.edge_count(), 1, "Should have 1 bond");

    // Check charge was applied to the oxygen atom (index 1)
    let atom2 = molecule.graph.node_weight(1.into()).unwrap();
    assert_eq!(atom2.charge, -1, "Oxygen should have -1 charge from M CHG");
}

#[test]
fn test_ctab_block_truncated_lines() {
    // Test with truncated atom lines (missing trailing fields)
    let ctab_data = b"  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
    1.5400    0.0000    0.0000 C   0  0
  1  2  1
M  END
";
    let result = ctab_block().parse(ctab_data);
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
fn test_ctab_block_missing_m_end() {
    // M END is optional according to spec
    let ctab_data = b"  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
";
    let result = ctab_block().parse(ctab_data);
    assert!(
        result.is_ok(),
        "CTAB block without M END should parse successfully"
    );

    let (_, molecule) = result.unwrap();
    assert_eq!(molecule.graph.node_count(), 2, "Should have 2 atoms");
}

#[test]
fn test_ctab_block_insufficient_atoms() {
    let ctab_data = b"  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
";
    let result = ctab_block().parse(ctab_data);
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
    let result = ctab_block().parse(ctab_data);
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
    let result = ctab_block().parse(ctab_data);
    assert!(
        result.is_ok(),
        "CTAB block with query features should parse successfully"
    );

    let (_, molecule) = result.unwrap();
    let atom1 = molecule.graph.node_weight(0.into()).unwrap();
    assert!(matches!(atom1.symbol, AtomSymbol::AtomList(_)));
}

#[test]
fn test_ctab_block_standard_fails_on_query() {
    let ctab_data = b"  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 L   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  ALS   1  2 F Cl  Br
M  END
";
    let result = ctab_block_standard().parse(ctab_data);
    assert!(
        result.is_err(),
        "Standard parser should fail on query features"
    );
}
