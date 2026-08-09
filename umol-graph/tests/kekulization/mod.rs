use std::collections::HashSet;

use rstest::rstest;
use umol_chem::element::Element;
use umol_graph::ops::invariant::{ValenceInvariants, ValenceMismatch};
use umol_graph::ops::transform::{KekulizationConfig, Kekulizer, KekulizerError, Transformer};
use umol_graph_ir::dsl::{MoleculeDefaults, MoleculeDsl};
use umol_graph_ir::ir::{
    AromaticSystemId, AtomConstraintKey, AtomId, BondConstraintKey, BondId, ElectronCountsAst,
    ElementAst, IntoIr, NumForm, UnpairedElectronsAst,
};
use umol_utils::solution::Solution;

#[derive(Debug, PartialEq, Eq)]
struct KekulizationFixture {
    participants: Vec<AtomId>,
    electrons: Vec<i64>,
    charge: i64,
    unpaired_electrons: UnpairedElectronsAst,
    elements: Vec<ElementAst>,
    nonzero_atom_charges: Vec<(AtomId, i64)>,
    nonzero_lone_pairs: Vec<(AtomId, i64)>,
}

#[rstest]
#[case::benzene(
    include_str!("data/benzene_aromatic_input.edn"),
    KekulizationFixture {
        participants: (0..6).map(AtomId).collect(),
        electrons: vec![1, 1, 1, 1, 1, 1],
        charge: 0,
        unpaired_electrons: UnpairedElectronsAst::from((0, 1)),
        elements: vec![ElementAst::Lit(Element::C); 6],
        nonzero_atom_charges: vec![],
        nonzero_lone_pairs: vec![],
    }
)]
#[case::pyridine(
    include_str!("data/pyridine_aromatic_input.edn"),
    KekulizationFixture {
        participants: (0..6).map(AtomId).collect(),
        electrons: vec![1, 1, 1, 1, 1, 1],
        charge: 0,
        unpaired_electrons: UnpairedElectronsAst::from((0, 1)),
        elements: vec![
            ElementAst::Lit(Element::N),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
        ],
        nonzero_atom_charges: vec![],
        nonzero_lone_pairs: vec![(AtomId(0), 1)],
    }
)]
#[case::pyrrole(
    include_str!("data/pyrrole_aromatic_input.edn"),
    KekulizationFixture {
        participants: (0..5).map(AtomId).collect(),
        electrons: vec![2, 1, 1, 1, 1],
        charge: 0,
        unpaired_electrons: UnpairedElectronsAst::from((0, 1)),
        elements: vec![
            ElementAst::Lit(Element::N),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
        ],
        nonzero_atom_charges: vec![],
        nonzero_lone_pairs: vec![(AtomId(0), 1)],
    }
)]
#[case::furan(
    include_str!("data/furan_aromatic_input.edn"),
    KekulizationFixture {
        participants: (0..5).map(AtomId).collect(),
        electrons: vec![2, 1, 1, 1, 1],
        charge: 0,
        unpaired_electrons: UnpairedElectronsAst::from((0, 1)),
        elements: vec![
            ElementAst::Lit(Element::O),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
        ],
        nonzero_atom_charges: vec![],
        nonzero_lone_pairs: vec![(AtomId(0), 2)],
    }
)]
#[case::thiophene(
    include_str!("data/thiophene_aromatic_input.edn"),
    KekulizationFixture {
        participants: (0..5).map(AtomId).collect(),
        electrons: vec![2, 1, 1, 1, 1],
        charge: 0,
        unpaired_electrons: UnpairedElectronsAst::from((0, 1)),
        elements: vec![
            ElementAst::Lit(Element::S),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
        ],
        nonzero_atom_charges: vec![],
        nonzero_lone_pairs: vec![(AtomId(0), 2)],
    }
)]
#[case::borepin(
    include_str!("data/borepin_aromatic_input.edn"),
    KekulizationFixture {
        participants: (0..7).map(AtomId).collect(),
        electrons: vec![0, 1, 1, 1, 1, 1, 1],
        charge: 0,
        unpaired_electrons: UnpairedElectronsAst::from((0, 1)),
        elements: vec![
            ElementAst::Lit(Element::B),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
        ],
        nonzero_atom_charges: vec![],
        nonzero_lone_pairs: vec![],
    }
)]
#[case::boratabenzene(
    include_str!("data/boratabenzene_aromatic_input.edn"),
    KekulizationFixture {
        participants: (0..6).map(AtomId).collect(),
        electrons: vec![1, 1, 1, 1, 1, 1],
        charge: 0,
        unpaired_electrons: UnpairedElectronsAst::from((0, 1)),
        elements: vec![
            ElementAst::Lit(Element::B),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
        ],
        nonzero_atom_charges: vec![(AtomId(0), -1)],
        nonzero_lone_pairs: vec![],
    }
)]
#[case::cyclopentadienyl_anion(
    include_str!("data/cyclopentadienyl_anion_aromatic_input.edn"),
    KekulizationFixture {
        participants: (0..5).map(AtomId).collect(),
        electrons: vec![1, 1, 1, 1, 1],
        charge: -1,
        unpaired_electrons: UnpairedElectronsAst::from((0, 1)),
        elements: vec![ElementAst::Lit(Element::C); 5],
        nonzero_atom_charges: vec![],
        nonzero_lone_pairs: vec![],
    }
)]
#[case::tropylium(
    include_str!("data/tropylium_aromatic_input.edn"),
    KekulizationFixture {
        participants: (0..7).map(AtomId).collect(),
        electrons: vec![1, 1, 1, 1, 1, 1, 1],
        charge: 1,
        unpaired_electrons: UnpairedElectronsAst::from((0, 1)),
        elements: vec![ElementAst::Lit(Element::C); 7],
        nonzero_atom_charges: vec![],
        nonzero_lone_pairs: vec![],
    }
)]
#[case::azulene(
    include_str!("data/azulene_aromatic_input.edn"),
    KekulizationFixture {
        participants: (0..10).map(AtomId).collect(),
        electrons: vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        charge: 0,
        unpaired_electrons: UnpairedElectronsAst::from((0, 1)),
        elements: vec![ElementAst::Lit(Element::C); 10],
        nonzero_atom_charges: vec![],
        nonzero_lone_pairs: vec![],
    }
)]
#[case::indole_prescribed_donor(
    include_str!("data/indole_prescribed_donor_aromatic_input.edn"),
    KekulizationFixture {
        participants: (0..9).map(AtomId).collect(),
        electrons: vec![2, 1, 1, 1, 1, 1, 1, 1, 1],
        charge: 0,
        unpaired_electrons: UnpairedElectronsAst::from((0, 1)),
        elements: vec![
            ElementAst::Lit(Element::N),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
            ElementAst::Lit(Element::C),
        ],
        nonzero_atom_charges: vec![],
        nonzero_lone_pairs: vec![(AtomId(0), 1)],
    }
)]
fn test_kekulization_fixture(#[case] source: &str, #[case] expected: KekulizationFixture) {
    let dsl: MoleculeDsl = source.parse().unwrap();
    let molecule = dsl.into_ir(&MoleculeDefaults::ground());

    assert!(molecule.is_ground());
    assert_eq!(molecule.aromatic_systems().count(), 1);

    let system = molecule.aromatic_systems().iter().next().unwrap();
    let ElectronCountsAst::Lit(electrons) = system.electrons() else {
        panic!("fixture aromatic electron contributions are undetermined");
    };
    let NumForm::Lit(charge) = system.charge() else {
        panic!("fixture aromatic charge is undetermined");
    };
    let actual = KekulizationFixture {
        participants: system.atom_ids().collect(),
        electrons: electrons.clone(),
        charge: *charge,
        unpaired_electrons: system.unpaired_electrons().clone(),
        elements: molecule
            .atoms()
            .iter()
            .map(|atom| atom.ast.element.clone())
            .collect(),
        nonzero_atom_charges: molecule
            .atoms()
            .iter()
            .filter_map(|atom| match atom.ast.charge {
                NumForm::Lit(0) => None,
                NumForm::Lit(charge) => Some((atom.id, charge)),
                _ => panic!("fixture atom charge is undetermined"),
            })
            .collect(),
        nonzero_lone_pairs: molecule
            .atoms()
            .iter()
            .filter_map(|atom| match atom.ast.lone_pairs {
                NumForm::Lit(0) => None,
                NumForm::Lit(lone_pairs) => Some((atom.id, lone_pairs)),
                _ => panic!("fixture atom lone-pair count is undetermined"),
            })
            .collect(),
    };
    assert_eq!(actual, expected);
}

#[rstest]
#[case::benzene(
    include_str!("data/benzene_aromatic_input.edn"),
    include_str!("data/benzene_kekulized_expected.edn"),
    vec![BondId(0), BondId(2), BondId(4)],
    None
)]
#[case::pyridine(
    include_str!("data/pyridine_aromatic_input.edn"),
    include_str!("data/pyridine_kekulized_expected.edn"),
    vec![BondId(0), BondId(2), BondId(4)],
    None
)]
#[case::pyrrole(
    include_str!("data/pyrrole_aromatic_input.edn"),
    include_str!("data/pyrrole_kekulized_expected.edn"),
    vec![BondId(1), BondId(3)],
    Some(AtomId(0))
)]
#[case::furan(
    include_str!("data/furan_aromatic_input.edn"),
    include_str!("data/furan_kekulized_expected.edn"),
    vec![BondId(1), BondId(3)],
    Some(AtomId(0))
)]
#[case::thiophene(
    include_str!("data/thiophene_aromatic_input.edn"),
    include_str!("data/thiophene_kekulized_expected.edn"),
    vec![BondId(1), BondId(3)],
    Some(AtomId(0))
)]
#[case::borepin(
    include_str!("data/borepin_aromatic_input.edn"),
    include_str!("data/borepin_kekulized_expected.edn"),
    vec![BondId(1), BondId(3), BondId(5)],
    Some(AtomId(0))
)]
#[case::boratabenzene(
    include_str!("data/boratabenzene_aromatic_input.edn"),
    include_str!("data/boratabenzene_kekulized_expected.edn"),
    vec![BondId(0), BondId(2), BondId(4)],
    None
)]
#[case::cyclopentadienyl_anion(
    include_str!("data/cyclopentadienyl_anion_aromatic_input.edn"),
    include_str!("data/cyclopentadienyl_anion_kekulized_expected.edn"),
    vec![BondId(0), BondId(2)],
    Some(AtomId(4))
)]
#[case::tropylium(
    include_str!("data/tropylium_aromatic_input.edn"),
    include_str!("data/tropylium_kekulized_expected.edn"),
    vec![BondId(0), BondId(2), BondId(4)],
    Some(AtomId(6))
)]
#[case::azulene(
    include_str!("data/azulene_aromatic_input.edn"),
    include_str!("data/azulene_kekulized_expected.edn"),
    vec![BondId(0), BondId(2), BondId(5), BondId(7), BondId(9)],
    None
)]
#[case::indole_prescribed_donor(
    include_str!("data/indole_prescribed_donor_aromatic_input.edn"),
    include_str!("data/indole_prescribed_donor_kekulized_expected.edn"),
    vec![BondId(1), BondId(3), BondId(6), BondId(8)],
    Some(AtomId(0))
)]
fn test_kekulization_fixture_output(
    #[case] source: &str,
    #[case] expected_source: &str,
    #[case] expected_double_bonds: Vec<BondId>,
    #[case] expected_exposed_atom: Option<AtomId>,
) {
    let input_dsl: MoleculeDsl = source.parse().unwrap();
    let input = input_dsl.into_ir(&MoleculeDefaults::ground());
    let expected_dsl: MoleculeDsl = expected_source.parse().unwrap();
    let expected = expected_dsl.into_ir(&MoleculeDefaults::ground());
    let node_order: Vec<AtomId> = input.atoms().iter().map(|atom| atom.id).collect();
    let kekulizer = Kekulizer::new(KekulizationConfig::default(), node_order);

    let first = kekulizer.transform(&input).unwrap();
    let second = kekulizer.transform(&input).unwrap();
    let double_bonds: Vec<BondId> = first
        .bonds()
        .iter()
        .filter(|bond| bond.ast.order == NumForm::Lit(2))
        .map(|bond| bond.id)
        .collect();
    let covered_atoms: HashSet<AtomId> = first
        .bonds()
        .iter()
        .filter(|bond| bond.ast.order == NumForm::Lit(2))
        .flat_map(|bond| bond.atom_ids())
        .collect();
    let expected_covered_atoms: HashSet<AtomId> = first
        .atoms()
        .ids()
        .filter(|atom| Some(*atom) != expected_exposed_atom)
        .collect();
    let input_total_charge: i64 = input
        .atoms()
        .iter()
        .map(|atom| match atom.ast.charge {
            NumForm::Lit(charge) => charge,
            _ => panic!("input atom charge is undetermined"),
        })
        .chain(
            input
                .aromatic_systems()
                .iter()
                .map(|system| match system.ast.charge {
                    NumForm::Lit(charge) => charge,
                    _ => panic!("input aromatic-system charge is undetermined"),
                }),
        )
        .sum();
    let output_total_charge: i64 = first
        .atoms()
        .iter()
        .map(|atom| match atom.ast.charge {
            NumForm::Lit(charge) => charge,
            _ => panic!("output atom charge is undetermined"),
        })
        .sum();
    let system_charge = match input.aromatic_systems().iter().next().unwrap().ast.charge {
        NumForm::Lit(charge) => charge,
        _ => panic!("input aromatic-system charge is undetermined"),
    };

    assert_eq!(first, expected);
    assert_eq!(second, expected);
    assert_eq!(double_bonds, expected_double_bonds);
    assert_eq!(covered_atoms, expected_covered_atoms);
    assert_eq!(output_total_charge, input_total_charge);
    assert_eq!(ValenceInvariants::check(&first), Solution::Determined(()));
    assert_eq!(first.aromatic_systems().count(), 0);
    assert!(first.atoms().iter().all(|atom| !atom
        .ast
        .constraints
        .contains(AtomConstraintKey::AromaticValence)));
    assert!(first
        .bonds()
        .iter()
        .all(|bond| !bond.ast.constraints.contains(BondConstraintKey::Aromatic)));

    if let Some(exposed) = expected_exposed_atom {
        let before = input.atom(exposed).ast;
        let after = first.atom(exposed).ast;
        let NumForm::Lit(before_charge) = before.charge else {
            panic!("input exposed-atom charge is undetermined");
        };
        assert_eq!(after.charge, NumForm::Lit(before_charge + system_charge));
        let NumForm::Lit(before_lone_pairs) = before.lone_pairs else {
            panic!("input exposed-atom lone pairs are undetermined");
        };
        let expected_lone_pairs = before_lone_pairs + i64::from(system_charge == -1);
        assert_eq!(after.lone_pairs, NumForm::Lit(expected_lone_pairs));
        assert_eq!(after.implicit_hydrogens, before.implicit_hydrogens);
        assert_eq!(after.unpaired_electrons, before.unpaired_electrons);
        if system_charge != 0 {
            assert_eq!(after.unpaired_electrons, UnpairedElectronsAst::from((0, 1)));
        }
    }
}

#[rstest]
#[case::localization(
    NumForm::Undetermined,
    NumForm::Lit(0),
    KekulizerError::UndeterminedExposedAtomCharge {
        system: AromaticSystemId(0),
        atom: AtomId(4),
    }
)]
#[case::lone_pair_localization(
    NumForm::Lit(0),
    NumForm::Undetermined,
    KekulizerError::UndeterminedExposedAtomLonePairs {
        system: AromaticSystemId(0),
        atom: AtomId(4),
    }
)]
#[case::valence_invariant(
    NumForm::Lit(0),
    NumForm::Lit(1),
    KekulizerError::PostLocalizationValenceInvariant(ValenceMismatch::OrbitalCount {
        atom_id: AtomId(4),
        orbital_count: 10,
        electron_count: 8,
    })
)]
fn test_kekulization_fixture_output_error(
    #[case] exposed_charge: NumForm,
    #[case] exposed_lone_pairs: NumForm,
    #[case] expected: KekulizerError,
) {
    let dsl: MoleculeDsl = include_str!("data/cyclopentadienyl_anion_aromatic_input.edn")
        .parse()
        .unwrap();
    let mut input = dsl.into_ir(&MoleculeDefaults::ground());
    input.atom_mut(AtomId(4)).ast.charge = exposed_charge;
    input.atom_mut(AtomId(4)).ast.lone_pairs = exposed_lone_pairs;
    let original = input.clone();
    let node_order: Vec<AtomId> = input.atoms().ids().collect();

    let result =
        Kekulizer::new(KekulizationConfig::default(), node_order).transform_into(&mut input);

    assert_eq!(result, Err(expected));
    assert_eq!(input, original);
}
