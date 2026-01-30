use indexmap::IndexMap;
use pretty_assertions::assert_eq;
use umol_data::Element;

use super::*;
use crate::table_ir::{
    Atom, AtomStereoCare, Bond, BondOrder, ConversionError, ExtendedAtom, ExtendedBond, RGroup,
    RGroupOccurrence, SGroup, SGroupType, SourceFormat,
};

#[test]
fn test_molecule_empty() {
    let mol = Molecule::empty();
    assert!(mol.atoms.is_empty());
    assert!(mol.bonds.is_empty());
    assert!(mol.rings.is_empty());
    assert_eq!(mol.source_format, SourceFormat::UNKNOWN);
}

#[test]
fn test_molecule_with_atoms_and_bonds() {
    let mol = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::O),
        ],
        bonds: vec![Bond::new(0, 1, BondOrder::Single)],
        rings: vec![],
        positions: None,
        comments: vec![],
        properties: IndexMap::new(),
        stereo_interpretation: None,
        source_format: SourceFormat::SMILES,
    };
    assert_eq!(mol.atom_count(), 2);
    assert_eq!(mol.bond_count(), 1);
    assert_eq!(mol.bonds[0].start_atom(), 0);
    assert_eq!(mol.bonds[0].end_atom(), 1);
}

#[test]
fn test_molecule_atom_count() {
    let mol = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::C),
            Atom::aromatic_atom(Element::N),
            Atom::aliphatic_atom(Element::O),
        ],
        bonds: vec![],
        rings: vec![],
        positions: None,
        comments: vec![],
        stereo_interpretation: None,
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    assert_eq!(mol.atom_count(), 3);
}

#[test]
fn test_molecule_bond_count() {
    let mol = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::C),
        ],
        bonds: vec![
            Bond::new(0, 1, BondOrder::Single),
            Bond::new(1, 2, BondOrder::Double),
        ],
        rings: vec![],
        positions: None,
        comments: vec![],
        stereo_interpretation: None,
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    assert_eq!(mol.bond_count(), 2);
}

#[test]
fn test_molecule_sum_formula_simple() {
    let mol = Molecule {
        atoms: vec![Atom::aliphatic_atom(Element::C)],
        bonds: vec![],
        rings: vec![],
        positions: None,
        comments: vec![],
        stereo_interpretation: None,
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    assert_eq!(mol.sum_formula(), "C");
}

#[test]
fn test_molecule_sum_formula_with_hydrogen() {
    let mol = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::H),
            Atom::aliphatic_atom(Element::H),
            Atom::aliphatic_atom(Element::H),
            Atom::aliphatic_atom(Element::H),
        ],
        bonds: vec![],
        rings: vec![],
        positions: None,
        comments: vec![],
        stereo_interpretation: None,
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    assert_eq!(mol.sum_formula(), "CH4");
}

#[test]
fn test_molecule_sum_formula_with_charge() {
    let mut atom = Atom::aliphatic_atom(Element::C);
    atom.charge = Some(1);
    let mol = Molecule {
        atoms: vec![atom],
        bonds: vec![],
        rings: vec![],
        positions: None,
        comments: vec![],
        stereo_interpretation: None,
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    let formula = mol.sum_formula();
    assert_eq!(formula, "C+");
}

#[test]
fn test_molecule_sum_formula_multiple_elements() {
    let mol = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::C),
            Atom::aliphatic_atom(Element::O),
            Atom::aliphatic_atom(Element::N),
        ],
        bonds: vec![],
        rings: vec![],
        positions: None,
        comments: vec![],
        stereo_interpretation: None,
        properties: IndexMap::new(),
        source_format: SourceFormat::SMILES,
    };
    let formula = mol.sum_formula();
    assert_eq!(formula, "C2NO");
}

#[test]
fn test_extended_molecule_empty() {
    let ext = ExtendedMolecule::empty();
    assert!(ext.atoms.is_empty());
    assert!(ext.bonds.is_empty());
    assert!(ext.rings.is_empty());
    assert!(ext.fragments.is_empty());
    assert!(ext.links.is_empty());
    assert!(ext.properties.is_empty());
    assert!(ext.comments.is_empty());
    assert_eq!(ext.source_format, SourceFormat::UNKNOWN);
    assert!(ext.ctfile_data.is_none());
    assert!(ext.electrons.is_none());
}

#[test]
fn test_extended_molecule_direct_construction() {
    let ext = ExtendedMolecule {
        atoms: vec![ExtendedAtom::from_element(Element::C)],
        bonds: vec![ExtendedBond::new(0, 0, BondOrder::Single)],
        rings: vec![],
        positions: None,
        fragments: vec![],
        links: vec![],
        electrons: Some(0),
        comments: vec!["test".to_string()],
        stereo_interpretation: None,
        properties: IndexMap::new(),
        ctfile_data: None,
        cx_data: None,
        source_format: SourceFormat::SMILES,
    };

    assert_eq!(ext.atom_count(), 1);
    assert_eq!(ext.bond_count(), 1);
    assert_eq!(ext.electrons, Some(0));
    assert_eq!(ext.comments.len(), 1);
}

#[test]
fn test_extended_molecule_from_molecule() {
    let mol = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::C),
            Atom::aromatic_atom(Element::N),
        ],
        bonds: vec![Bond::new(0, 1, BondOrder::Double)],
        rings: vec![],
        positions: None,
        comments: vec![],
        stereo_interpretation: None,
        properties: IndexMap::new(),
        source_format: SourceFormat::MOL,
    };

    let ext = ExtendedMolecule::from(mol);
    assert_eq!(ext.atoms.len(), 2);
    assert_eq!(ext.bonds.len(), 1);
    assert_eq!(ext.source_format, SourceFormat::MOL);
    assert!(ext.sgroups().is_empty());
    assert!(ext.rgroups().is_empty());
    assert!(ext.ctfile_data.is_none());
}

#[test]
fn test_extended_molecule_atom_count() {
    let mut ext = ExtendedMolecule::empty();
    ext.atoms.push(ExtendedAtom::from_element(Element::C));
    ext.atoms.push(ExtendedAtom::from_element(Element::N));
    assert_eq!(ext.atom_count(), 2);
}

#[test]
fn test_extended_molecule_bond_count() {
    let mut ext = ExtendedMolecule::empty();
    ext.bonds.push(ExtendedBond::new(0, 1, BondOrder::Single));
    ext.bonds.push(ExtendedBond::new(1, 2, BondOrder::Double));
    assert_eq!(ext.bond_count(), 2);
}

#[test]
fn test_extended_molecule_sum_formula() {
    let mut ext = ExtendedMolecule::empty();
    ext.atoms.push(ExtendedAtom::from_element(Element::C));
    ext.atoms.push(ExtendedAtom::from_element(Element::H));
    ext.atoms.push(ExtendedAtom::from_element(Element::H));
    ext.atoms.push(ExtendedAtom::from_element(Element::H));
    ext.atoms.push(ExtendedAtom::from_element(Element::H));
    assert_eq!(ext.sum_formula(), "CH4");
}

#[test]
fn test_extended_molecule_sum_formula_with_wildcards() {
    use super::super::atom::{AtomSymbol, WildcardAtom};

    let mut ext = ExtendedMolecule::empty();
    // *C(=O)OCC - ethyl ester with wildcard
    ext.atoms
        .push(ExtendedAtom::from_atom_symbol(AtomSymbol::WildcardAtom(
            WildcardAtom::Any,
        )));
    ext.atoms.push(ExtendedAtom::from_element(Element::C));
    ext.atoms.push(ExtendedAtom::from_element(Element::O));
    ext.atoms.push(ExtendedAtom::from_element(Element::O));
    ext.atoms.push(ExtendedAtom::from_element(Element::C));
    ext.atoms.push(ExtendedAtom::from_element(Element::C));
    assert_eq!(ext.sum_formula(), "C3O2*");

    // Two wildcards
    let mut ext2 = ExtendedMolecule::empty();
    ext2.atoms
        .push(ExtendedAtom::from_atom_symbol(AtomSymbol::WildcardAtom(
            WildcardAtom::Any,
        )));
    ext2.atoms.push(ExtendedAtom::from_element(Element::C));
    ext2.atoms
        .push(ExtendedAtom::from_atom_symbol(AtomSymbol::WildcardAtom(
            WildcardAtom::Any,
        )));
    assert_eq!(ext2.sum_formula(), "C*2");
}

#[test]
fn test_extended_molecule_sgroups() {
    let mut ext = ExtendedMolecule::empty();
    assert!(ext.sgroups().is_empty());

    ext.sgroups_mut()
        .insert(0, SGroup::new(SGroupType::Superatom));
    assert_eq!(ext.sgroups().len(), 1);
    assert!(ext.sgroups().contains_key(&0));
}

#[test]
fn test_extended_molecule_rgroups() {
    let mut ext = ExtendedMolecule::empty();
    assert!(ext.rgroups().is_empty());

    let rgroup = RGroup {
        label: Some(1),
        dependent_label: None,
        rgroup_or_h: false,
        occurrence: vec![RGroupOccurrence::Range(1, 3)],
    };
    ext.rgroups_mut().insert(1, rgroup);
    assert_eq!(ext.rgroups().len(), 1);
    assert!(ext.rgroups().contains_key(&1));
}

#[test]
fn test_extended_molecule_to_molecule() {
    let mut ext = ExtendedMolecule::empty();
    ext.atoms.push(ExtendedAtom::from_element(Element::C));
    ext.atoms.push(ExtendedAtom::from_element(Element::O));
    ext.bonds.push(ExtendedBond::new(0, 1, BondOrder::Single));
    ext.source_format = SourceFormat::SMILES;

    let mol = ext.to_molecule().unwrap();
    assert_eq!(mol.atoms.len(), 2);
    assert_eq!(mol.bonds.len(), 1);
    assert_eq!(mol.source_format, SourceFormat::SMILES);
}

#[test]
fn test_extended_molecule_to_molecule_error() {
    let mut ext = ExtendedMolecule::empty();
    let mut atom = ExtendedAtom::from_element(Element::C);
    atom.stereo_care = Some(AtomStereoCare::Care);
    ext.atoms.push(atom);

    let result = ext.to_molecule();
    assert!(result.is_err());
    assert_eq!(result.unwrap_err(), ConversionError::HasExtendedFeatures);
}
