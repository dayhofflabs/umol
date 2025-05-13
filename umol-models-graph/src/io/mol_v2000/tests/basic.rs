use super::super::*;
use float_cmp::*;
use std::io::Cursor;
use umol_data::Element;

#[test]
fn test_read_mol_minimal() {
    let mol_str = r#"
  Methane

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
"#;

    let cursor = Cursor::new(mol_str);
    let result = read_mol_v2000(cursor);

    assert!(result.is_ok(), "Parsing failed: {:?}", result.err());
    let molecule = result.unwrap();

    assert_eq!(molecule.atom_count(), 1, "Incorrect number of atoms");
    assert_eq!(molecule.bond_count(), 0, "Incorrect number of bonds");
    assert!(molecule.sgroups.is_empty(), "SGroups should be empty");

    // Check atom properties
    let atom = molecule.atom(0).unwrap();
    assert_eq!(atom.element, Element::C, "Atom element should be Carbon");
    assert_eq!(atom.charge, 0, "Atom charge should be 0");
    assert_eq!(atom.radical, None, "Atom radical should be None");
    assert_eq!(
        atom.isotope_mass, None,
        "Atom mass difference should be None"
    );
    assert_eq!(atom.stereo_parity, None, "Atom stereo should be None");
    assert_eq!(atom.valence, None, "Atom valence should be None"); // Default valence code 0 -> None
    assert_eq!(atom.hydrogen_count, None, "Atom H count should be None");
    assert_eq!(atom.atom_map_num, None, "Atom map number should be None");

    // Check conformer
    assert_eq!(molecule.conformers.len(), 1, "Should have one conformer");
    let conformer = &molecule.conformers[0]; // Access via indexing
    assert!(!conformer.is_3d, "Conformer should be marked as 2D");
    assert_eq!(
        conformer.positions.len(),
        1,
        "Conformer should have one position"
    );
    let pos = conformer.get_position(0).unwrap();
    assert_eq!(pos.x, 0.0, "Position x should be 0.0");
    assert_eq!(pos.y, 0.0, "Position y should be 0.0");
    assert_eq!(pos.z, 0.0, "Position z should be 0.0");
}


#[test]
fn test_read_mol_atom_props() {
    let mol_str = r#"
  -ISIS-  05110910502D

  4  3  0  0  0  0  0  0  0  0999 V2000
   -3.0000   -7.8750    0.0000 C   0  0  0  0  0  0  0  0  0  1  0  0
   -2.2855   -7.4625    0.0000 C   0  0  0  0  0  0  0  0  0  2  0  0
   -2.2855   -6.6375    0.0000 O   0  0  0  0  0  0  0  0  0  3  0  0
   -1.5711   -7.8750    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
  2  3  2  0  0  0  0
  1  2  1  0  0  0  0
  2  4  1  0  0  0  0
V    2 acidchloride
M  END
"#;

    let cursor = Cursor::new(mol_str);
    let result = read_mol_v2000(cursor);

    assert!(result.is_ok(), "Parsing failed: {:?}", result.err());
    let _molecule = result.unwrap();
}


#[test]
fn test_read_mol_radical() {
    let mol_str = r#"
  Methyl radical

  4  3  0     0  0  0  0  0  0999 V2000
    2.5369    0.1550    0.0000 C   0  4  0  0  0  0  0  0  0  0  0  0
    3.0739    0.4650    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    0.4650    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.5369   -0.4650    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
M  RAD  1   1   2
M  END
"#;

    let cursor = Cursor::new(mol_str);
    let result = read_mol_v2000(cursor);

    assert!(result.is_ok(), "Parsing failed: {:?}", result.err());
    let molecule = result.unwrap();

    assert_eq!(molecule.atom_count(), 4, "Incorrect number of atoms");
    assert_eq!(molecule.bond_count(), 3, "Incorrect number of bonds");

    // Check atom properties
    let atom = molecule.atom(0).unwrap();
    assert_eq!(atom.element, Element::C, "Atom element should be Carbon");
    assert_eq!(atom.charge, 0, "Atom charge should be 0");
    assert_eq!(atom.radical, Some(2), "Atom radical should be 2");
    assert_eq!(
        atom.isotope_mass, None,
        "Atom mass difference should be None"
    );
    assert_eq!(atom.stereo_parity, None, "Atom stereo should be None");

    // Check conformer
    assert_eq!(molecule.conformers.len(), 1, "Should have one conformer");
    let conformer = &molecule.conformers[0];

    // Check conformer position
    let pos = conformer.get_position(0).unwrap();
    assert!(approx_eq!(f64, 2.5369, pos.x, ulps = 4));
    assert!(approx_eq!(f64, 0.1550, pos.y, ulps = 4));
    assert!(approx_eq!(f64, 0.0, pos.z, ulps = 4));
}