use rstest::rstest;
use umol_ast::ast::{AtomId, ElectronCountsAst, ElementAst, IntoAst, SpinStateAst, ValueAst};
use umol_ast::dsl::{MoleculeDefaults, MoleculeDsl};
use umol_chem::element::Element;

#[derive(Debug, PartialEq, Eq)]
struct FixtureDemand {
    participants: Vec<AtomId>,
    electrons: Vec<i64>,
    charge: i64,
    spin: SpinStateAst,
    elements: Vec<ElementAst>,
    nonzero_atom_charges: Vec<(AtomId, i64)>,
    nonzero_lone_pairs: Vec<(AtomId, i64)>,
}

#[rstest]
#[case::benzene(
    include_str!("data/benzene_aromatic_input.edn"),
    FixtureDemand {
        participants: (0..6).map(AtomId).collect(),
        electrons: vec![1, 1, 1, 1, 1, 1],
        charge: 0,
        spin: SpinStateAst::from((0, 1)),
        elements: vec![ElementAst::Lit(Element::C); 6],
        nonzero_atom_charges: vec![],
        nonzero_lone_pairs: vec![],
    }
)]
#[case::pyridine(
    include_str!("data/pyridine_aromatic_input.edn"),
    FixtureDemand {
        participants: (0..6).map(AtomId).collect(),
        electrons: vec![1, 1, 1, 1, 1, 1],
        charge: 0,
        spin: SpinStateAst::from((0, 1)),
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
    FixtureDemand {
        participants: (0..5).map(AtomId).collect(),
        electrons: vec![2, 1, 1, 1, 1],
        charge: 0,
        spin: SpinStateAst::from((0, 1)),
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
    FixtureDemand {
        participants: (0..5).map(AtomId).collect(),
        electrons: vec![2, 1, 1, 1, 1],
        charge: 0,
        spin: SpinStateAst::from((0, 1)),
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
    FixtureDemand {
        participants: (0..5).map(AtomId).collect(),
        electrons: vec![2, 1, 1, 1, 1],
        charge: 0,
        spin: SpinStateAst::from((0, 1)),
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
    FixtureDemand {
        participants: (0..7).map(AtomId).collect(),
        electrons: vec![0, 1, 1, 1, 1, 1, 1],
        charge: 0,
        spin: SpinStateAst::from((0, 1)),
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
    FixtureDemand {
        participants: (0..6).map(AtomId).collect(),
        electrons: vec![1, 1, 1, 1, 1, 1],
        charge: 0,
        spin: SpinStateAst::from((0, 1)),
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
    FixtureDemand {
        participants: (0..5).map(AtomId).collect(),
        electrons: vec![1, 1, 1, 1, 1],
        charge: -1,
        spin: SpinStateAst::from((0, 1)),
        elements: vec![ElementAst::Lit(Element::C); 5],
        nonzero_atom_charges: vec![],
        nonzero_lone_pairs: vec![],
    }
)]
#[case::tropylium(
    include_str!("data/tropylium_aromatic_input.edn"),
    FixtureDemand {
        participants: (0..7).map(AtomId).collect(),
        electrons: vec![1, 1, 1, 1, 1, 1, 1],
        charge: 1,
        spin: SpinStateAst::from((0, 1)),
        elements: vec![ElementAst::Lit(Element::C); 7],
        nonzero_atom_charges: vec![],
        nonzero_lone_pairs: vec![],
    }
)]
#[case::azulene(
    include_str!("data/azulene_aromatic_input.edn"),
    FixtureDemand {
        participants: (0..10).map(AtomId).collect(),
        electrons: vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        charge: 0,
        spin: SpinStateAst::from((0, 1)),
        elements: vec![ElementAst::Lit(Element::C); 10],
        nonzero_atom_charges: vec![],
        nonzero_lone_pairs: vec![],
    }
)]
#[case::indole_prescribed_donor(
    include_str!("data/indole_prescribed_donor_aromatic_input.edn"),
    FixtureDemand {
        participants: (0..9).map(AtomId).collect(),
        electrons: vec![2, 1, 1, 1, 1, 1, 1, 1, 1],
        charge: 0,
        spin: SpinStateAst::from((0, 1)),
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
fn test_kekulization_fixture(#[case] source: &str, #[case] expected: FixtureDemand) {
    let dsl: MoleculeDsl = source.parse().unwrap();
    let molecule = dsl.into_ast(&MoleculeDefaults::ground());

    assert!(molecule.is_ground());
    assert_eq!(molecule.aromatic_systems().count(), 1);

    let system = molecule.aromatic_systems().iter().next().unwrap();
    let ElectronCountsAst::Lit(electrons) = system.electrons() else {
        panic!("fixture aromatic electron contributions are undetermined");
    };
    let ValueAst::Lit(charge) = system.charge() else {
        panic!("fixture aromatic charge is undetermined");
    };
    let actual = FixtureDemand {
        participants: system.atom_ids().collect(),
        electrons: electrons.clone(),
        charge: *charge,
        spin: system.spin().clone(),
        elements: molecule
            .atoms()
            .iter()
            .map(|atom| atom.ast.element.clone())
            .collect(),
        nonzero_atom_charges: molecule
            .atoms()
            .iter()
            .filter_map(|atom| match atom.ast.charge {
                ValueAst::Lit(0) => None,
                ValueAst::Lit(charge) => Some((atom.id, charge)),
                _ => panic!("fixture atom charge is undetermined"),
            })
            .collect(),
        nonzero_lone_pairs: molecule
            .atoms()
            .iter()
            .filter_map(|atom| match atom.ast.lone_pairs {
                ValueAst::Lit(0) => None,
                ValueAst::Lit(lone_pairs) => Some((atom.id, lone_pairs)),
                _ => panic!("fixture atom lone-pair count is undetermined"),
            })
            .collect(),
    };
    assert_eq!(actual, expected);
}
