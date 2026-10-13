use std::borrow::Cow;
use std::collections::HashSet;
use std::fs;
use std::path::PathBuf;

use rstest::rstest;
use umol_chem::element::Element;
use umol_edn::FromEdn;
use umol_graph::ops::aromaticity::AromaticityPerceiver;
use umol_graph::ops::invariant::ValenceMismatch;
use umol_graph::ops::model::{AromaticityModel, ChemistryModel, ValenceModel, ValenceTieBreak};
use umol_graph::ops::resolve::{IsotopePolicy, ResolveConfig, Resolver};
use umol_graph::ops::transform::{KekulizeConfig, KekulizeError, Kekulizer, Transformer};
use umol_graph::ops::valence::ValenceTable;
use umol_graph::ops::validate::ValenceInvariantsValidator;
use umol_graph::parse::parse_mol;
use umol_graph_ir::dsl::{MoleculeDefaults, MoleculeDsl, MoleculeOverrides};
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticSystemId, AtomConstraintKey, AtomId, BondConstraintKey, BondId,
    ElectronCountsForm, ElementForm, IntoIr, NumForm, UnpairedElectronsForm,
};
use umol_utils::solution::Solution;

#[derive(Debug, PartialEq, Eq)]
struct KekulizationFixture {
    participants: Vec<AtomId>,
    electrons: Vec<i64>,
    charge: i64,
    unpaired_electrons: UnpairedElectronsForm,
    elements: Vec<ElementForm>,
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
        unpaired_electrons: UnpairedElectronsForm::from((0, 1)),
        elements: vec![ElementForm::Lit(Element::C); 6],
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
        unpaired_electrons: UnpairedElectronsForm::from((0, 1)),
        elements: vec![
            ElementForm::Lit(Element::N),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
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
        unpaired_electrons: UnpairedElectronsForm::from((0, 1)),
        elements: vec![
            ElementForm::Lit(Element::N),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
        ],
        nonzero_atom_charges: vec![],
        nonzero_lone_pairs: vec![],
    }
)]
#[case::furan(
    include_str!("data/furan_aromatic_input.edn"),
    KekulizationFixture {
        participants: (0..5).map(AtomId).collect(),
        electrons: vec![2, 1, 1, 1, 1],
        charge: 0,
        unpaired_electrons: UnpairedElectronsForm::from((0, 1)),
        elements: vec![
            ElementForm::Lit(Element::O),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
        ],
        nonzero_atom_charges: vec![],
        nonzero_lone_pairs: vec![(AtomId(0), 1)],
    }
)]
#[case::thiophene(
    include_str!("data/thiophene_aromatic_input.edn"),
    KekulizationFixture {
        participants: (0..5).map(AtomId).collect(),
        electrons: vec![2, 1, 1, 1, 1],
        charge: 0,
        unpaired_electrons: UnpairedElectronsForm::from((0, 1)),
        elements: vec![
            ElementForm::Lit(Element::S),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
        ],
        nonzero_atom_charges: vec![],
        nonzero_lone_pairs: vec![(AtomId(0), 1)],
    }
)]
#[case::borepin(
    include_str!("data/borepin_aromatic_input.edn"),
    KekulizationFixture {
        participants: (0..7).map(AtomId).collect(),
        electrons: vec![0, 1, 1, 1, 1, 1, 1],
        charge: 0,
        unpaired_electrons: UnpairedElectronsForm::from((0, 1)),
        elements: vec![
            ElementForm::Lit(Element::B),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
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
        unpaired_electrons: UnpairedElectronsForm::from((0, 1)),
        elements: vec![
            ElementForm::Lit(Element::B),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
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
        unpaired_electrons: UnpairedElectronsForm::from((0, 1)),
        elements: vec![ElementForm::Lit(Element::C); 5],
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
        unpaired_electrons: UnpairedElectronsForm::from((0, 1)),
        elements: vec![ElementForm::Lit(Element::C); 7],
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
        unpaired_electrons: UnpairedElectronsForm::from((0, 1)),
        elements: vec![ElementForm::Lit(Element::C); 10],
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
        unpaired_electrons: UnpairedElectronsForm::from((0, 1)),
        elements: vec![
            ElementForm::Lit(Element::N),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
            ElementForm::Lit(Element::C),
        ],
        nonzero_atom_charges: vec![],
        nonzero_lone_pairs: vec![],
    }
)]
fn test_kekulizer_transform_input(#[case] source: &str, #[case] expected: KekulizationFixture) {
    let dsl: MoleculeDsl = source.parse().unwrap();
    let molecule = dsl.into_ir(&MoleculeDefaults::concrete());

    assert!(molecule.is_concrete());
    assert_eq!(molecule.aromatic_systems().count(), 1);

    let system = molecule.aromatic_systems().iter().next().unwrap();
    let ElectronCountsForm::Lit(electrons) = system.electrons() else {
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
            .map(|atom| atom.attributes().element.clone())
            .collect(),
        nonzero_atom_charges: molecule
            .atoms()
            .iter()
            .filter_map(|atom| match atom.attributes().charge {
                NumForm::Lit(0) => None,
                NumForm::Lit(charge) => Some((atom.id(), charge)),
                _ => panic!("fixture atom charge is undetermined"),
            })
            .collect(),
        nonzero_lone_pairs: molecule
            .atoms()
            .iter()
            .filter_map(|atom| match atom.attributes().lone_pairs {
                NumForm::Lit(0) => None,
                NumForm::Lit(lone_pairs) => Some((atom.id(), lone_pairs)),
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
fn test_kekulizer_transform(
    #[case] source: &str,
    #[case] expected_source: &str,
    #[case] expected_double_bonds: Vec<BondId>,
    #[case] expected_exposed_atom: Option<AtomId>,
) {
    let input_dsl: MoleculeDsl = source.parse().unwrap();
    let input = input_dsl.into_ir(&MoleculeDefaults::concrete());
    let expected_dsl: MoleculeDsl = expected_source.parse().unwrap();
    let expected = expected_dsl.into_ir(&MoleculeDefaults::concrete());
    let node_order: Vec<AtomId> = input.atoms().iter().map(|atom| atom.id()).collect();
    let kekulizer = Kekulizer::new(KekulizeConfig::default(), node_order);

    assert_eq!(
        ValenceInvariantsValidator.validate(&input),
        Ok(Solution::Determined(()))
    );
    let first = kekulizer.transform(input.clone()).unwrap();
    let second = kekulizer.transform(input.clone()).unwrap();
    let double_bonds: Vec<BondId> = first
        .bonds()
        .iter()
        .filter(|bond| bond.attributes().order == NumForm::Lit(2))
        .map(|bond| bond.id())
        .collect();
    let covered_atoms: HashSet<AtomId> = first
        .bonds()
        .iter()
        .filter(|bond| bond.attributes().order == NumForm::Lit(2))
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
        .map(|atom| match atom.attributes().charge {
            NumForm::Lit(charge) => charge,
            _ => panic!("input atom charge is undetermined"),
        })
        .chain(
            input
                .aromatic_systems()
                .iter()
                .map(|system| match system.attributes().charge {
                    NumForm::Lit(charge) => charge,
                    _ => panic!("input aromatic-system charge is undetermined"),
                }),
        )
        .sum();
    let output_total_charge: i64 = first
        .atoms()
        .iter()
        .map(|atom| match atom.attributes().charge {
            NumForm::Lit(charge) => charge,
            _ => panic!("output atom charge is undetermined"),
        })
        .sum();
    let system_charge = match input
        .aromatic_systems()
        .iter()
        .next()
        .unwrap()
        .attributes()
        .charge
    {
        NumForm::Lit(charge) => charge,
        _ => panic!("input aromatic-system charge is undetermined"),
    };

    assert_eq!(first, expected);
    assert_eq!(second, expected);
    assert_eq!(double_bonds, expected_double_bonds);
    assert_eq!(covered_atoms, expected_covered_atoms);
    assert_eq!(output_total_charge, input_total_charge);
    assert_eq!(
        ValenceInvariantsValidator.validate(&first),
        Ok(Solution::Determined(()))
    );
    assert_eq!(first.aromatic_systems().count(), 0);
    assert!(first.atoms().iter().all(|atom| !atom
        .attributes()
        .constraints
        .contains(AtomConstraintKey::AromaticValence)));
    assert!(first.bonds().iter().all(|bond| !bond
        .attributes()
        .constraints
        .contains(BondConstraintKey::Aromatic)));

    if let Some(exposed) = expected_exposed_atom {
        let before = input.atom(exposed).attributes();
        let after = first.atom(exposed).attributes();
        let NumForm::Lit(before_charge) = before.charge else {
            panic!("input exposed-atom charge is undetermined");
        };
        assert_eq!(after.charge, NumForm::Lit(before_charge + system_charge));
        assert_eq!(after.lone_pairs, *expected.atom(exposed).lone_pairs());
        assert_eq!(after.implicit_hydrogens, before.implicit_hydrogens);
        assert_eq!(after.unpaired_electrons, before.unpaired_electrons);
        if system_charge != 0 {
            assert_eq!(
                after.unpaired_electrons,
                UnpairedElectronsForm::from((0, 1))
            );
        }
    }
}

#[rstest]
#[case::c60(include_str!(
    "../../../umol-io/tests/mol_parsing/data/molecule/scifinder/99685-96-8.mol"
))]
fn test_kekulization_fixture_output_cage(#[case] source: &str) {
    let mut input = parse_mol(source).unwrap();
    let system_atoms: Vec<AtomId> = input.atoms().ids().collect();
    let system_bonds: Vec<BondId> = input.bonds().ids().collect();
    for bond in system_bonds {
        input.bond_mut(bond).attributes_mut().order = NumForm::Lit(1);
    }
    AromaticityPerceiver::new(&AromaticityModel::daylight()).add_systems(
        &mut input,
        vec![(
            system_atoms.clone(),
            AromaticSystemForm {
                electrons: ElectronCountsForm::Lit(vec![1; 60]),
                charge: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm::closed_shell(),
                ..Default::default()
            },
        )],
    );
    let original = input.clone();
    let node_order: Vec<AtomId> = input.atoms().ids().collect();
    let kekulizer = Kekulizer::new(KekulizeConfig::default(), node_order);

    let first = kekulizer.transform(input.clone()).unwrap();
    let second = kekulizer.transform(input.clone()).unwrap();
    let double_bonds: Vec<BondId> = first
        .bonds()
        .iter()
        .filter(|bond| bond.attributes().order == NumForm::Lit(2))
        .map(|bond| bond.id())
        .collect();
    let covered_atoms: HashSet<AtomId> = first
        .bonds()
        .iter()
        .filter(|bond| bond.attributes().order == NumForm::Lit(2))
        .flat_map(|bond| bond.atom_ids())
        .collect();
    let expected_covered_atoms: HashSet<AtomId> = system_atoms.into_iter().collect();

    assert_eq!(input, original);
    assert_eq!(first, second);
    assert_eq!(double_bonds.len(), 30);
    assert_eq!(covered_atoms, expected_covered_atoms);
    assert_eq!(
        ValenceInvariantsValidator.validate(&first),
        Ok(Solution::Determined(()))
    );
    assert_eq!(first.aromatic_systems().count(), 0);
    assert!(first.bonds().iter().all(|bond| {
        matches!(bond.attributes().order, NumForm::Lit(1 | 2))
            && !bond
                .attributes()
                .constraints
                .contains(BondConstraintKey::Aromatic)
    }));
}

#[rstest]
#[case::localization(
    NumForm::Undetermined,
    NumForm::Lit(0),
    KekulizeError::UndeterminedExposedAtomCharge {
        system: AromaticSystemId(0),
        atom: AtomId(4),
    }
)]
#[case::lone_pair_localization(
    NumForm::Lit(0),
    NumForm::Undetermined,
    KekulizeError::UndeterminedExposedAtomLonePairs {
        system: AromaticSystemId(0),
        atom: AtomId(4),
    }
)]
#[case::valence_invariant(
    NumForm::Lit(0),
    NumForm::Lit(1),
    KekulizeError::PostLocalizationValenceInvariant(ValenceMismatch::OrbitalCount {
        atom_id: AtomId(4),
        orbital_count: 10,
        electron_count: 8,
    })
)]
fn test_kekulization_fixture_output_error(
    #[case] exposed_charge: NumForm,
    #[case] exposed_lone_pairs: NumForm,
    #[case] expected: KekulizeError,
) {
    let dsl: MoleculeDsl = include_str!("data/cyclopentadienyl_anion_aromatic_input.edn")
        .parse()
        .unwrap();
    let mut input = dsl.into_ir(&MoleculeDefaults::concrete());
    input.atom_mut(AtomId(4)).attributes_mut().charge = exposed_charge;
    input.atom_mut(AtomId(4)).attributes_mut().lone_pairs = exposed_lone_pairs;
    let original = input.clone();
    let node_order: Vec<AtomId> = input.atoms().ids().collect();

    let result = Kekulizer::new(KekulizeConfig::default(), node_order).transform_into(&mut input);

    assert_eq!(result, Err(expected));
    assert_eq!(input, original);
}

#[derive(FromEdn)]
struct ResolutionFixture {
    input: MoleculeDsl,
    #[edn(default)]
    config_overrides: MoleculeOverrides,
}

#[rstest]
fn test_kekulizer_transform_corpus(
    #[files("tests/resolution/data/aromatic/*.edn")] path: PathBuf,
    #[values(false, true)] counts: bool,
) {
    let fixture = ResolutionFixture::from_edn_str(&fs::read_to_string(&path).unwrap()).unwrap();
    let defaults = MoleculeDefaults::default().with_overrides(fixture.config_overrides);
    let mut model = ChemistryModel {
        valence: if counts {
            ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table()))
        } else {
            ValenceModel::default()
        },
        ..Default::default()
    };
    model.valence.tie_break = ValenceTieBreak::MostSaturated;
    let resolver = Resolver::with_config(
        &model,
        ResolveConfig {
            isotope: IsotopePolicy::Natural,
            ..Default::default()
        },
    );
    let solution = resolver.resolve(fixture.input.into_ir(&defaults)).unwrap();
    let Solution::Determined(input) = solution else {
        panic!("corpus resolution failed: {solution:?}");
    };
    assert_eq!(
        ValenceInvariantsValidator.validate(&input),
        Ok(Solution::Determined(()))
    );
    let kekulizer = Kekulizer::new(KekulizeConfig::default(), input.atoms().ids().collect());
    let output = kekulizer.transform(input.clone());
    let mut in_place = input.clone();
    let result = kekulizer.transform_into(&mut in_place);
    if path.file_stem().unwrap() == "f420" {
        assert_eq!(output, Err(KekulizeError::NoMatching(AromaticSystemId(0))));
        assert_eq!(result, Err(KekulizeError::NoMatching(AromaticSystemId(0))));
        assert_eq!(in_place, input);
        assert_eq!(kekulizer.transform_iter(&input).collect::<Vec<_>>(), vec![]);
    } else {
        let output = output.unwrap();
        assert_eq!(result, Ok(()));
        assert_eq!(in_place, output);
        assert_eq!(
            ValenceInvariantsValidator.validate(&output),
            Ok(Solution::Determined(()))
        );
        assert_eq!(output.aromatic_systems().count(), 0);
        assert_eq!(output.atoms().count(), input.atoms().count());
        assert_eq!(output.bonds().count(), input.bonds().count());
        for atom in input.atoms().iter() {
            let localized = output.atom(atom.id());
            assert_eq!(localized.element(), atom.element());
            assert_eq!(localized.charge(), atom.charge());
            assert_eq!(localized.implicit_hydrogens(), atom.implicit_hydrogens());
            assert_eq!(localized.unpaired_electrons(), atom.unpaired_electrons());
        }
    }
}

#[rstest]
#[case::f420(include_str!("../resolution/data/aromatic/f420.edn"))]
fn test_kekulizer_transform_hydrogens(#[case] source: &str, #[values(false, true)] counts: bool) {
    let fixture = ResolutionFixture::from_edn_str(source).unwrap();
    let defaults = MoleculeDefaults::default().with_overrides(fixture.config_overrides);
    let mut input = fixture.input.into_ir(&defaults);
    input
        .atom_mut(AtomId(17))
        .attributes_mut()
        .implicit_hydrogens = NumForm::Lit(0);
    input
        .atom_mut(AtomId(20))
        .attributes_mut()
        .implicit_hydrogens = NumForm::Lit(1);
    let mut model = ChemistryModel {
        valence: if counts {
            ValenceModel::counts(Cow::Borrowed(ValenceTable::default_table()))
        } else {
            ValenceModel::default()
        },
        ..Default::default()
    };
    model.valence.tie_break = ValenceTieBreak::MostSaturated;
    let resolver = Resolver::with_config(
        &model,
        ResolveConfig {
            isotope: IsotopePolicy::Natural,
            ..Default::default()
        },
    );
    let solution = resolver.resolve(input).unwrap();
    let Solution::Determined(input) = solution else {
        panic!("corpus resolution failed: {solution:?}");
    };
    assert_eq!(
        ValenceInvariantsValidator.validate(&input),
        Ok(Solution::Determined(()))
    );
    assert_eq!(
        input.aromatic_system(AromaticSystemId(0)).electrons(),
        &ElectronCountsForm::Lit(vec![2, 1, 1, 0, 2, 0, 1, 1, 1, 1, 1, 1, 1, 1])
    );
    let kekulizer = Kekulizer::new(KekulizeConfig::default(), input.atoms().ids().collect());
    let output = kekulizer.transform(input.clone()).unwrap();
    let mut in_place = input;
    assert_eq!(kekulizer.transform_into(&mut in_place), Ok(()));
    assert_eq!(in_place, output);
    assert_eq!(
        ValenceInvariantsValidator.validate(&output),
        Ok(Solution::Determined(()))
    );
    assert_eq!(output.aromatic_systems().count(), 0);
    for atom in [AtomId(15), AtomId(17), AtomId(20)] {
        assert_eq!(output.atom(atom).lone_pairs(), &NumForm::Lit(1));
    }
    for bond in [BondId(18), BondId(25), BondId(28), BondId(30), BondId(33)] {
        assert_eq!(output.bond(bond).order(), &NumForm::Lit(2));
    }
}

#[rstest]
#[case::pyrrole(include_str!("data/pyrrole_aromatic_input.edn"))]
#[case::furan(include_str!("data/furan_aromatic_input.edn"))]
fn test_kekulizer_transform_donor_error(#[case] source: &str) {
    let mut input = source
        .parse::<MoleculeDsl>()
        .unwrap()
        .into_ir(&MoleculeDefaults::concrete());
    input.atom_mut(AtomId(0)).attributes_mut().lone_pairs = NumForm::Undetermined;
    let expected = KekulizeError::UndeterminedExposedAtomLonePairs {
        system: AromaticSystemId(0),
        atom: AtomId(0),
    };
    let kekulizer = Kekulizer::new(KekulizeConfig::default(), input.atoms().ids().collect());
    assert_eq!(kekulizer.transform(input.clone()), Err(expected.clone()));
    let original = input.clone();
    assert_eq!(kekulizer.transform_into(&mut input), Err(expected));
    assert_eq!(input, original);
}
