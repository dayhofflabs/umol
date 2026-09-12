use rstest::rstest;
use umol_chem::element::Element;
use umol_chem::spin::SpinMultiplicity;

use super::{atom_winding, render, RenderError};
use crate::smiles::{Smiles, SmilesIoConfig};
use crate::table_ir::{
    Atom, Bond, BondConfiguration, BondDonation, BondOrder, Molecule, StereoAtom, StereoBond,
    StereoLigand, Winding,
};

#[rstest]
#[case::empty("", "")]
#[case::chain("CCC", "CCC")]
#[case::branches("CC(C)(O)N", "CC(C)(O)N")]
#[case::nested("CC(C(C)O)N", "CC(C(C)O)N")]
#[case::components("C.O.N", "C.O.N")]
#[case::branch_before_component("CC(C)O.N", "CC(C)O.N")]
#[case::ring("C%12CCCCC%12", "C1CCCCC1")]
#[case::aromatic("c1ccccc1", "c1ccccc1")]
#[case::aromatic_link("c1ccccc1-c2ccccc2", "c1ccccc1-c1ccccc1")]
#[case::explicit_aromatic("c1:c:c:c:c:c:1", "c1ccccc1")]
#[case::explicit_single("C-C", "CC")]
#[case::fixed_hydrogens("[CH3][CH2][OH]", "[CH3][CH2][OH]")]
#[case::zero_hydrogens("[CH0]", "[C]")]
#[case::explicit_hydrogen("[H]C", "[H]C")]
#[case::stereo_explicit_hydrogen("[H][C@](F)(Cl)Br", "[H][C@](F)(Cl)Br")]
#[case::isotope_charge_class("[13CH2-:17]", "[13CH2-:17]")]
#[case::multiple_charge("[Fe++]", "[Fe+2]")]
#[case::zero_charge("[C+0]", "[C]")]
#[case::wildcard("*", "*")]
#[case::bracket_wildcard("[*:3]", "[*:3]")]
#[case::selenium("[seH]", "[seH]")]
#[case::arsenic("[as]", "[as]")]
#[case::double_triple("CC#CC=C", "CC#CC=C")]
#[case::quadruple("[Mo]$[Mo]", "[Mo]$[Mo]")]
#[case::tetrahedral("F[C@](Cl)(Br)I", "F[C@](Cl)(Br)I")]
#[case::tetrahedral_class("F[C@TH2](Cl)(Br)I", "F[C@@](Cl)(Br)I")]
#[case::root_hydrogen("[C@H](F)(Cl)Br", "[C@H](F)(Cl)Br")]
#[case::nonroot_hydrogen("F[C@@H](Cl)Br", "F[C@@H](Cl)Br")]
#[case::root_lone_pair("[S@](=O)(C)CC", "[S@](=O)(C)CC")]
#[case::nonroot_lone_pair("C[S@@](=O)CC", "C[S@@](=O)CC")]
#[case::ring_hydrogen("C[C@H]1CCCCO1", "C[C@H]1CCCCO1")]
#[case::ring_root("[C@H]1(F)CCCCO1", "[C@H]1(F)CCCCO1")]
#[case::two_stereo_closures("N[C@]12CCCC1CCC2", "N[C@]12CCCC1CCC2")]
#[case::trans("F/C=C/F", "F/C=C/F")]
#[case::cis("F/C=C\\F", "F/C=C\\F")]
#[case::global_reversal("F\\C=C\\F", "F/C=C/F")]
#[case::partial("F/C=CF", "FC=CF")]
#[case::four_substituents("F/C(Cl)=C(Br)/I", "FC(/Cl)=C(Br)/I")]
#[case::shared_chain("C/C=C/C=C/C", "C/C=C/C=C/C")]
#[case::ring_markers("C1/C=C/C=C/CCC1", "C1/C=C/C=C/CCC1")]
fn test_render(#[case] input: &str, #[case] expected: &str) {
    let molecule = Smiles::parse(input).unwrap().into_table_ir();
    let output = render(&molecule, &SmilesIoConfig::opensmiles()).unwrap();
    assert_eq!(output, expected);
    Smiles::parse(&output).unwrap();
}

#[rstest]
#[case::donating("N->[Cu]", "N->[Cu]")]
#[case::accepting("N<-[Cu]", "N<-[Cu]")]
#[case::ring_donation("N1CC[Cu]->1", "N<-1CC[Cu]1")]
#[case::any("C~C", "C~C")]
#[case::tellurium("[te]", "[te]")]
#[case::silicon("[siH]", "[siH]")]
fn test_render_extensions(#[case] input: &str, #[case] expected: &str) {
    let config = SmilesIoConfig::lenient();
    let molecule = Smiles::parse_with(input, &config).unwrap().into_table_ir();
    let output = render(&molecule, &config).unwrap();
    assert_eq!(output, expected);
    Smiles::parse_with(&output, &config).unwrap();
}

#[rstest]
#[case::unchanged(vec![StereoLigand::Atom(0),StereoLigand::Atom(2),StereoLigand::Atom(3),StereoLigand::Atom(4)], Some(0), vec![0,2,3,4], true, Ok(Winding::CounterClockwise))]
#[case::odd(vec![StereoLigand::Atom(0),StereoLigand::Atom(2),StereoLigand::Atom(3),StereoLigand::Atom(4)], Some(0), vec![0,3,2,4], true, Ok(Winding::Clockwise))]
#[case::even(vec![StereoLigand::Atom(0),StereoLigand::Atom(2),StereoLigand::Atom(3),StereoLigand::Atom(4)], Some(0), vec![2,0,4,3], true, Ok(Winding::CounterClockwise))]
#[case::root_hydrogen(vec![StereoLigand::Atom(0),StereoLigand::ImplicitHydrogen,StereoLigand::Atom(2),StereoLigand::Atom(3)], Some(1), vec![0,2,3], false, Ok(Winding::Clockwise))]
#[case::root_lone_pair(vec![StereoLigand::Atom(0),StereoLigand::LonePair,StereoLigand::Atom(2),StereoLigand::Atom(3)], Some(0), vec![0,2,3], false, Ok(Winding::Clockwise))]
#[case::wrong_ligand(vec![StereoLigand::Atom(0),StereoLigand::Atom(2),StereoLigand::Atom(3),StereoLigand::Atom(5)], Some(0), vec![0,2,3,4], true, Err(RenderError::InvalidStereoAtom { atom: 1 }))]
#[case::missing_hydrogen(vec![StereoLigand::Atom(0),StereoLigand::ImplicitHydrogen,StereoLigand::Atom(2),StereoLigand::Atom(3)], Some(0), vec![0,2,3], true, Err(RenderError::InvalidStereoAtom { atom: 1 }))]
#[case::unknown_hydrogens(vec![StereoLigand::Atom(0),StereoLigand::Atom(2),StereoLigand::Atom(3),StereoLigand::Atom(4)], None, vec![0,2,3,4], true, Err(RenderError::InvalidStereoAtom { atom: 1 }))]
#[case::duplicate_ligand(vec![StereoLigand::Atom(0),StereoLigand::Atom(2),StereoLigand::Atom(3),StereoLigand::Atom(3)], Some(0), vec![0,2,3,4], true, Err(RenderError::InvalidStereoAtom { atom: 1 }))]
fn test_atom_winding(
    #[case] ligands: Vec<StereoLigand>,
    #[case] hydrogens: Option<u8>,
    #[case] neighbors: Vec<u32>,
    #[case] parent: bool,
    #[case] expected: Result<Winding, RenderError>,
) {
    let atom = Atom {
        implicit_hydrogens: hydrogens,
        ..Atom::aliphatic_atom(Element::C)
    };
    let frame = StereoAtom {
        atom: 1,
        ligands,
        winding: Winding::CounterClockwise,
    };
    assert_eq!(atom_winding(&frame, &atom, &neighbors, parent), expected);
}

#[rstest]
#[case::missing_atom(vec![Bond::new(0,2,BondOrder::Single)], RenderError::AtomIndexOutOfBounds { atom: 2 })]
#[case::self_bond(vec![Bond::new(0,0,BondOrder::Single)], RenderError::SelfBond { bond: 0 })]
#[case::parallel(vec![Bond::new(0,1,BondOrder::Single),Bond::new(0,1,BondOrder::Single)], RenderError::DuplicateBond { bond: 0 })]
#[case::zero(vec![Bond::new(0,1,BondOrder::Zero)], RenderError::UnsupportedBond { bond: 0, field: "order" })]
#[case::charge(vec![Bond { charge: Some(1), ..Bond::new(0,1,BondOrder::Single) }], RenderError::UnsupportedBond { bond: 0, field: "charge" })]
#[case::spin(vec![Bond { multiplicity: Some(SpinMultiplicity::TRIPLET), ..Bond::new(0,1,BondOrder::Single) }], RenderError::UnsupportedBond { bond: 0, field: "multiplicity" })]
#[case::donation(vec![Bond::new_dative(0,1,BondOrder::Single,BondDonation::Donating)], RenderError::UnsupportedBond { bond: 0, field: "donation" })]
fn test_render_bond_error(#[case] bonds: Vec<Bond>, #[case] expected: RenderError) {
    let molecule = Molecule {
        atoms: vec![Atom::aliphatic_atom(Element::C); 2],
        bonds,
        ..Molecule::empty()
    };
    assert_eq!(
        render(&molecule, &SmilesIoConfig::opensmiles()),
        Err(expected)
    );
}

#[rstest]
#[case::isotope_inferred_h(Atom { isotope_mass: Some(13), ..Atom::aliphatic_atom(Element::C) }, RenderError::InferredHydrogens { atom: 0 })]
#[case::h_count(Atom { implicit_hydrogens: Some(10), ..Atom::aliphatic_atom(Element::C) }, RenderError::UnsupportedAtom { atom: 0, field: "implicit_hydrogens" })]
#[case::charge(Atom { charge: Some(-128), ..Atom::aliphatic_atom(Element::C) }, RenderError::UnsupportedAtom { atom: 0, field: "charge" })]
#[case::electron_count(Atom { lone_pairs: Some(1), ..Atom::aliphatic_atom(Element::C) }, RenderError::UnsupportedAtom { atom: 0, field: "lone_pairs" })]
#[case::aromatic(Atom::aromatic_atom(Element::Fe), RenderError::UnsupportedAtom { atom: 0, field: "aromatic" })]
fn test_render_atom_error(#[case] atom: Atom, #[case] expected: RenderError) {
    let molecule = Molecule {
        atoms: vec![atom],
        ..Molecule::empty()
    };
    assert_eq!(
        render(&molecule, &SmilesIoConfig::opensmiles()),
        Err(expected)
    );
}

#[rstest]
fn test_render_either() {
    let mut molecule = Smiles::parse("FC=CF").unwrap().into_table_ir();
    molecule.stereo_bonds.push(StereoBond {
        bond: 1,
        configuration: BondConfiguration::Either,
    });
    assert_eq!(
        render(&molecule, &SmilesIoConfig::opensmiles()),
        Err(RenderError::UnsupportedBond {
            bond: 1,
            field: "Either"
        })
    );
}

#[rstest]
#[case::percent(12, Ok("C123456789%10CC1C2C3C4C5C6C7C8C9C%10"))]
#[case::limit(102, Err(RenderError::RingLabel { label: 100 }))]
fn test_render_ring_labels(#[case] count: usize, #[case] expected: Result<&str, RenderError>) {
    let molecule = Molecule {
        atoms: vec![Atom::aliphatic_atom(Element::C); count],
        bonds: (1..count)
            .map(|atom| Bond::new(atom as u32 - 1, atom as u32, BondOrder::Single))
            .chain((2..count).map(|atom| Bond::new(0, atom as u32, BondOrder::Single)))
            .collect(),
        ..Molecule::empty()
    };
    let output = render(&molecule, &SmilesIoConfig::opensmiles());
    assert_eq!(output.as_deref(), expected.as_deref());
    if let Ok(output) = output {
        Smiles::parse(&output).unwrap();
    }
}

#[rstest]
#[case::donating(BondDonation::Donating, "NC<-[Cu]")]
#[case::accepting(BondDonation::Accepting, "NC->[Cu]")]
fn test_render_donation(#[case] donation: BondDonation, #[case] expected: &str) {
    let molecule = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::N),
            Atom {
                implicit_hydrogens: Some(0),
                ..Atom::aliphatic_atom(Element::Cu)
            },
            Atom::aliphatic_atom(Element::C),
        ],
        bonds: vec![
            Bond::new_dative(1, 2, BondOrder::Single, donation),
            Bond::new(0, 2, BondOrder::Single),
        ],
        ..Molecule::empty()
    };
    assert_eq!(
        render(&molecule, &SmilesIoConfig::lenient()),
        Ok(expected.into())
    );
}

#[rstest]
#[case::positions(Molecule { positions: Some(vec![]), ..Molecule::empty() }, "positions")]
#[case::comments(Molecule { comments: vec!["sample".into()], ..Molecule::empty() }, "comments")]
fn test_render_molecule_error(#[case] molecule: Molecule, #[case] field: &'static str) {
    assert_eq!(
        render(&molecule, &SmilesIoConfig::opensmiles()),
        Err(RenderError::UnsupportedMolecule { field })
    );
}

#[rstest]
#[case::missing(5, RenderError::AtomIndexOutOfBounds { atom: 5 })]
#[case::duplicate(1, RenderError::DuplicateStereoAtom { atom: 1 })]
fn test_render_stereo_atom_error(#[case] atom: u32, #[case] expected: RenderError) {
    let mut molecule = Smiles::parse("F[C@](Cl)(Br)I").unwrap().into_table_ir();
    let mut frame = molecule.stereo_atoms[0].clone();
    frame.atom = atom;
    molecule.stereo_atoms.push(frame);
    assert_eq!(
        render(&molecule, &SmilesIoConfig::opensmiles()),
        Err(expected)
    );
}

#[rstest]
#[case::opening(BondDonation::Donating, "N->1C[Cu]1")]
#[case::reverse(BondDonation::Accepting, "N<-1C[Cu]1")]
fn test_render_ring_donation(#[case] donation: BondDonation, #[case] expected: &str) {
    let molecule = Molecule {
        atoms: vec![
            Atom::aliphatic_atom(Element::N),
            Atom::aliphatic_atom(Element::C),
            Atom {
                implicit_hydrogens: Some(0),
                ..Atom::aliphatic_atom(Element::Cu)
            },
        ],
        bonds: vec![
            Bond::new(0, 1, BondOrder::Single),
            Bond::new(1, 2, BondOrder::Single),
            Bond::new_dative(0, 2, BondOrder::Single, donation),
        ],
        ..Molecule::empty()
    };
    let output = render(&molecule, &SmilesIoConfig::lenient()).unwrap();
    assert_eq!(output, expected);
    Smiles::parse_with(&output, &SmilesIoConfig::lenient()).unwrap();
}
