//! Aggregate canonicalization and participant-frame benchmarks.
//!
//! Criterion ids include the measured graph's node and edge counts.

use std::fmt::Debug;
use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use umol_chem::element::Element;
use umol_edn::FromEdn;
use umol_graph_core::{AutomorphismAlgorithm, EdgeId, GraphRemapping, NodeId, Remapping};
use umol_graph_ir::dsl::MoleculeDsl;
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticSystemId, AtomConstraintForm, AtomDelta, AtomFieldChange, AtomForm,
    AtomId, BondForm, BondId, Canonicalize, CanonicalizeContext, Constraint, DativeBondForm,
    DativeBondId, Delta, Deltas, IncidenceLevel, Molecule, MoleculeConstraint, MoleculeEntries,
    MoleculeRemapping, MulticenterBondForm, MulticenterBondId, NoncovalentBondForm,
    NoncovalentBondId, NoncovalentBondKind, NumForm, Reaction, Reframe, StereoAtomConstraintForm,
    StereoAtomDelta, StereoAtomFieldChange, StereoAtomForm, StereoAtomId, StereoBondDelta,
    StereoBondFieldChange, StereoBondForm, StereoBondId, StereoConfigurationForm, StereoCoset,
    StereoKind, StereoLigand, StereoLigandKind, StereoLigandPair, Topicity, TopicityForm,
    TopicityRelationForm,
};

const ALGORITHM: AutomorphismAlgorithm = AutomorphismAlgorithm::Nauty;

struct CorpusCase {
    name: &'static str,
    molecule: Molecule,
}

// Extended C/H propane network, seed flask 0.
const FEATURE_FREE_CONNECTED: &str = r#"
{:atoms ["H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "C#i=#c0#h0#n0#u0#s"
         "C#i=#c0#h0#n0#u0#s"
         "C#i=#c0#h0#n0#u0#s"]
 :bonds [[0 8 "1#c0#u0#s"]
         [1 8 "1#c0#u0#s"]
         [2 8 "1#c0#u0#s"]
         [3 9 "1#c0#u0#s"]
         [4 9 "1#c0#u0#s"]
         [5 9 "1#c0#u0#s"]
         [6 10 "1#c0#u0#s"]
         [7 10 "1#c0#u0#s"]
         [8 10 "1#c0#u0#s"]
         [9 10 "1#c0#u0#s"]]}
"#;

// Extended C/H propane network, product flask 72.
const FEATURE_FREE_DISCONNECTED: &str = r#"
{:atoms ["H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "H#i=#c0#h0#n0#u0#s"
         "C#i=#c0#h0#n0#u0#s"
         "C#i=#c0#h0#n0#u#s2"
         "C#i=#c0#h0#n0#u#s2"]
 :bonds [[0 8 "1#c0#u0#s"]
         [1 8 "1#c0#u0#s"]
         [2 8 "1#c0#u0#s"]
         [3 8 "1#c0#u0#s"]
         [4 9 "1#c0#u0#s"]
         [5 9 "1#c0#u0#s"]
         [6 10 "1#c0#u0#s"]
         [7 10 "1#c0#u0#s"]
         [9 10 "1#c0#u0#s"]]}
"#;

// Extended C/H ethane network, product flask 99.
const SYMMETRY_HEAVY_RADICALS: &str = r#"
{:atoms ["H#i=#c0#h0#n0#u#s2"
         "H#i=#c0#h0#n0#u#s2"
         "H#i=#c0#h0#n0#u#s2"
         "H#i=#c0#h0#n0#u#s2"
         "H#i=#c0#h0#n0#u#s2"
         "H#i=#c0#h0#n0#u#s2"
         "C#i=#c0#h0#n0#u#s2"
         "C#i=#c0#h0#n0#u#s2"]
 :bonds [[6 7 "3#c0#u0#s"]]}
"#;

const SCALING_CASES: [(&str, &str, [usize; 8]); 3] = [
    (
        "feature_free_connected",
        FEATURE_FREE_CONNECTED,
        [11, 10, 0, 0, 0, 0, 0, 0],
    ),
    (
        "feature_free_disconnected",
        FEATURE_FREE_DISCONNECTED,
        [11, 9, 0, 0, 0, 0, 0, 0],
    ),
    (
        "symmetry_heavy_radicals",
        SYMMETRY_HEAVY_RADICALS,
        [8, 1, 0, 0, 0, 0, 0, 0],
    ),
];

fn atom(element: Element) -> AtomForm {
    AtomForm::from_element(element)
}

fn bond(first: u32, second: u32, order: u8) -> (AtomId, AtomId, BondForm) {
    (AtomId(first), AtomId(second), BondForm::from_order(order))
}

fn ligand(atom: u32) -> StereoLigand {
    StereoLigand::new(AtomId(atom), StereoLigandKind::Atom)
}

fn implicit_hydrogen(site: u32) -> StereoLigand {
    StereoLigand::new(AtomId(site), StereoLigandKind::ImplicitHydrogen)
}

fn carbon_graph(atom_count: usize, edges: &[(u32, u32)]) -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![atom(Element::C); atom_count],
        bonds: edges
            .iter()
            .map(|&(first, second)| bond(first, second, 1))
            .collect(),
        ..Default::default()
    })
}

fn ordinary_naphthalene() -> Molecule {
    carbon_graph(
        10,
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (5, 6),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 4),
        ],
    )
}

fn disconnected_rings() -> Molecule {
    carbon_graph(
        12,
        &[
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 0),
            (6, 7),
            (7, 8),
            (8, 9),
            (9, 10),
            (10, 11),
            (11, 6),
        ],
    )
}

fn carbon_path(atom_count: usize) -> Molecule {
    carbon_graph(
        atom_count,
        &(0..atom_count as u32 - 1)
            .map(|first| (first, first + 1))
            .collect::<Vec<_>>(),
    )
}

fn overlay_heavy() -> Molecule {
    let mut atoms = [
        Element::C,
        Element::C,
        Element::C,
        Element::C,
        Element::N,
        Element::O,
        Element::F,
        Element::Cl,
    ]
    .into_iter()
    .map(atom)
    .collect::<Vec<_>>();
    atoms[0].constraints = AtomConstraintForm::valence(4).into();

    Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds: vec![
            bond(0, 1, 1),
            bond(1, 2, 2),
            bond(2, 3, 1),
            bond(3, 0, 1),
            bond(1, 4, 1),
            bond(1, 5, 1),
            bond(2, 6, 1),
            bond(2, 7, 1),
        ],
        dative: vec![(
            vec![AtomId(5), AtomId(4)],
            AtomId(3),
            DativeBondForm::from_order(1),
        )],
        aromatic: vec![(
            vec![AtomId(3), AtomId(2), AtomId(1), AtomId(0)],
            AromaticSystemForm::default(),
        )],
        multicenter: vec![(
            vec![AtomId(5), AtomId(4), AtomId(0)],
            MulticenterBondForm::default(),
        )],
        noncovalent: vec![(
            [AtomId(7), AtomId(6)],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        stereo_atoms: vec![(
            AtomId(1),
            vec![ligand(5), ligand(4), ligand(2), ligand(0)],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        stereo_bonds: vec![(
            BondId(1),
            vec![ligand(6), ligand(3), ligand(4), ligand(0)],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
        )],
        constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3)]),
            sum: NumForm::Lit(0),
        })
        .into(),
    })
}

fn overlay_unconstrained() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::N,
            Element::O,
            Element::F,
            Element::Cl,
        ]
        .into_iter()
        .map(atom)
        .collect(),
        bonds: vec![
            bond(0, 1, 1),
            bond(1, 2, 2),
            bond(2, 3, 1),
            bond(3, 0, 1),
            bond(1, 4, 1),
            bond(1, 5, 1),
            bond(2, 6, 1),
            bond(2, 7, 1),
        ],
        dative: vec![(
            vec![AtomId(5), AtomId(4)],
            AtomId(3),
            DativeBondForm::from_order(1),
        )],
        aromatic: vec![(
            vec![AtomId(3), AtomId(2), AtomId(1), AtomId(0)],
            AromaticSystemForm::default(),
        )],
        multicenter: vec![(
            vec![AtomId(5), AtomId(4), AtomId(0)],
            MulticenterBondForm::default(),
        )],
        noncovalent: vec![(
            [AtomId(7), AtomId(6)],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        ..Default::default()
    })
}

fn aromatic_electron_counts() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![atom(Element::C); 5],
        bonds: vec![
            bond(0, 1, 1),
            bond(1, 2, 1),
            bond(2, 3, 1),
            bond(3, 4, 1),
            bond(4, 0, 1),
        ],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4)],
            AromaticSystemForm::from_electrons(vec![2, 1, 1, 1, 1]),
        )],
        ..Default::default()
    })
}

fn charge_sum_ring(atoms: Option<Vec<AtomId>>) -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![atom(Element::C); 6],
        bonds: vec![
            bond(0, 1, 1),
            bond(1, 2, 1),
            bond(2, 3, 1),
            bond(3, 4, 1),
            bond(4, 5, 1),
            bond(5, 0, 1),
        ],
        constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms,
            sum: NumForm::Lit(0),
        })
        .into(),
        ..Default::default()
    })
}

fn symmetric_charge_sum() -> Molecule {
    charge_sum_ring(None)
}

fn asymmetric_charge_sum() -> Molecule {
    charge_sum_ring(Some(vec![AtomId(0), AtomId(1)]))
}

fn large_aromatic() -> Molecule {
    const ATOM_COUNT: usize = 128;

    Molecule::from_entries(MoleculeEntries {
        atoms: vec![atom(Element::C); ATOM_COUNT],
        aromatic: vec![(
            (0..ATOM_COUNT).rev().map(AtomId::from).collect::<Vec<_>>(),
            AromaticSystemForm::from_electrons((0..ATOM_COUNT).map(|value| value as i64).collect()),
        )],
        ..Default::default()
    })
}

fn reframing_corpus() -> Vec<CorpusCase> {
    vec![
        CorpusCase {
            name: "empty",
            molecule: Molecule::default(),
        },
        CorpusCase {
            name: "overlay_heavy",
            molecule: overlay_heavy(),
        },
        CorpusCase {
            name: "large_aromatic",
            molecule: large_aromatic(),
        },
    ]
}

fn tetrahedral_stereo() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .into_iter()
            .map(atom)
            .collect(),
        bonds: vec![bond(0, 1, 1), bond(0, 2, 1), bond(0, 3, 1), bond(0, 4, 1)],
        stereo_atoms: vec![(
            AtomId(0),
            vec![ligand(1), ligand(2), ligand(3), ligand(4)],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    })
}

fn cis_trans_stereo_bond() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [
            Element::C,
            Element::C,
            Element::F,
            Element::Cl,
            Element::Br,
            Element::I,
        ]
        .into_iter()
        .map(atom)
        .collect(),
        bonds: vec![
            bond(0, 1, 2),
            bond(0, 2, 1),
            bond(0, 3, 1),
            bond(1, 4, 1),
            bond(1, 5, 1),
        ],
        stereo_bonds: vec![(
            BondId(0),
            vec![ligand(3), ligand(2), ligand(5), ligand(4)],
            StereoBondForm::new(StereoKind::CisTrans, 1u32),
        )],
        ..Default::default()
    })
}

fn mixed_atom_and_bond_stereo() -> Molecule {
    let stereo_atom = tetrahedral_stereo();
    let stereo_bond = cis_trans_stereo_bond();
    Molecule::combine_all([&stereo_atom, &stereo_bond])
}

fn frame_relative_stereo_constraint() -> Molecule {
    let constraint = StereoAtomConstraintForm::Topicity(TopicityForm {
        pair: StereoLigandPair::new(0usize.into(), 2usize.into()),
        relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
    });
    Molecule::from_entries(MoleculeEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .into_iter()
            .map(atom)
            .collect(),
        bonds: vec![bond(0, 1, 1), bond(0, 2, 1), bond(0, 3, 1), bond(0, 4, 1)],
        stereo_atoms: vec![(
            AtomId(0),
            vec![ligand(1), ligand(2), ligand(3), ligand(4)],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))
                .with_constraint(constraint.clone()),
        )],
        constraints: Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral, constraint)
            .into(),
        ..Default::default()
    })
}

fn meso_dichlorobutane() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::Cl,
            Element::Cl,
        ]
        .into_iter()
        .map(atom)
        .collect(),
        bonds: vec![
            bond(0, 1, 1),
            bond(0, 2, 1),
            bond(1, 3, 1),
            bond(0, 4, 1),
            bond(1, 5, 1),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                vec![ligand(1), ligand(2), ligand(4), implicit_hydrogen(0)],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(1),
                vec![ligand(0), ligand(3), ligand(5), implicit_hydrogen(1)],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            ),
        ],
        ..Default::default()
    })
}

fn para_stereo_trichloropentane() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::Cl,
            Element::Cl,
            Element::Cl,
        ]
        .into_iter()
        .map(atom)
        .collect(),
        bonds: vec![
            bond(0, 1, 1),
            bond(1, 2, 1),
            bond(2, 3, 1),
            bond(3, 4, 1),
            bond(1, 5, 1),
            bond(2, 6, 1),
            bond(3, 7, 1),
        ],
        stereo_atoms: vec![
            (
                AtomId(1),
                vec![ligand(0), ligand(2), ligand(5), implicit_hydrogen(1)],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(2),
                vec![ligand(1), ligand(3), ligand(6), implicit_hydrogen(2)],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(3),
                vec![ligand(2), ligand(4), ligand(7), implicit_hydrogen(3)],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            ),
        ],
        ..Default::default()
    })
}

fn para_stereo_cascade() -> Molecule {
    let outer_ligands = vec![ligand(10), ligand(11), ligand(12), ligand(13)];

    Molecule::from_entries(MoleculeEntries {
        atoms: [
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::F,
            Element::Cl,
            Element::Br,
            Element::I,
        ]
        .into_iter()
        .map(atom)
        .collect(),
        bonds: [
            (0, [2, 3, 4, 5]),
            (1, [6, 8, 7, 9]),
            (2, [10, 11, 12, 13]),
            (3, [10, 11, 12, 13]),
            (4, [10, 11, 12, 13]),
            (5, [10, 11, 12, 13]),
            (6, [10, 11, 12, 13]),
            (7, [10, 11, 12, 13]),
            (8, [10, 11, 12, 13]),
            (9, [10, 11, 12, 13]),
        ]
        .into_iter()
        .flat_map(|(site, ligands)| ligands.into_iter().map(move |ligand| bond(site, ligand, 1)))
        .collect(),
        stereo_atoms: vec![
            (
                AtomId(0),
                [2, 3, 4, 5].map(ligand).into(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(1),
                [6, 8, 7, 9].map(ligand).into(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(2),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(3),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
            ),
            (
                AtomId(4),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Axial, 0u32),
            ),
            (
                AtomId(5),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::SquarePlanar, 0u32),
            ),
            (
                AtomId(6),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(7),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
            ),
            (
                AtomId(8),
                outer_ligands.clone(),
                StereoAtomForm::new(StereoKind::Axial, 0u32),
            ),
            (
                AtomId(9),
                outer_ligands,
                StereoAtomForm::new(StereoKind::SquarePlanar, 0u32),
            ),
        ],
        ..Default::default()
    })
}

/// Triisopropanolamine with three equally configured stereocentres, as a heavy-atom skeleton.
///
/// Each arm is a methylene carbon, the stereocentre, its oxygen, and its methyl carbon. The three
/// arms are equivalent, so the stereo-blind graph has six automorphisms and each preserves the
/// stereo atoms.
fn symmetric_stereo_triisopropanolamine() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [Element::N]
            .into_iter()
            .chain([Element::C, Element::C, Element::O, Element::C].repeat(3))
            .map(atom)
            .collect(),
        bonds: [1, 5, 9]
            .into_iter()
            .flat_map(|methylene| {
                [
                    bond(0, methylene, 1),
                    bond(methylene, methylene + 1, 1),
                    bond(methylene + 1, methylene + 2, 1),
                    bond(methylene + 1, methylene + 3, 1),
                ]
            })
            .collect(),
        stereo_atoms: [2, 6, 10]
            .into_iter()
            .map(|site| {
                (
                    AtomId(site),
                    vec![
                        ligand(site - 1),
                        ligand(site + 1),
                        ligand(site + 2),
                        implicit_hydrogen(site),
                    ],
                    StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
                )
            })
            .collect(),
        ..Default::default()
    })
}

/// Butan-2-ol with explicit hydrogens and its stereocentre configured.
///
/// Atoms 0 to 3 are the carbon chain, atom 4 is the oxygen, and atoms 5 to 14 are the hydrogens.
/// The hydrogens of the two methyl groups and of the methylene group are exchangeable, so the
/// stereo-blind graph has 6 * 2 * 6 = 72 automorphisms and each preserves the stereo atom.
fn symmetric_stereo_butanol() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [Element::C, Element::C, Element::C, Element::C, Element::O]
            .into_iter()
            .chain([Element::H; 10])
            .map(atom)
            .collect(),
        bonds: vec![
            bond(0, 1, 1),
            bond(1, 2, 1),
            bond(2, 3, 1),
            bond(1, 4, 1),
            bond(0, 5, 1),
            bond(0, 6, 1),
            bond(0, 7, 1),
            bond(1, 8, 1),
            bond(2, 9, 1),
            bond(2, 10, 1),
            bond(3, 11, 1),
            bond(3, 12, 1),
            bond(3, 13, 1),
            bond(4, 14, 1),
        ],
        stereo_atoms: vec![(
            AtomId(1),
            vec![ligand(0), ligand(2), ligand(4), ligand(8)],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    })
}

/// 3,3-Dimethylbutan-2-ol with explicit hydrogens and its stereocentre configured.
///
/// Atom 0 is the methyl carbon, atom 1 the stereocentre, atom 2 the quaternary carbon, atoms 3 to 5
/// the carbons of the tert-butyl group, atom 6 the oxygen, and atoms 7 to 20 the hydrogens. The
/// three tert-butyl carbons and the hydrogens of each of the four methyl groups are exchangeable,
/// so the stereo-blind graph has 6 * 6 * 6^3 = 7776 automorphisms and each preserves the stereo
/// atom.
fn symmetric_stereo_dimethylbutanol() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [Element::C; 6]
            .into_iter()
            .chain([Element::O])
            .chain([Element::H; 14])
            .map(atom)
            .collect(),
        bonds: vec![
            bond(0, 1, 1),
            bond(1, 2, 1),
            bond(2, 3, 1),
            bond(2, 4, 1),
            bond(2, 5, 1),
            bond(1, 6, 1),
            bond(0, 7, 1),
            bond(0, 8, 1),
            bond(0, 9, 1),
            bond(1, 10, 1),
            bond(3, 11, 1),
            bond(3, 12, 1),
            bond(3, 13, 1),
            bond(4, 14, 1),
            bond(4, 15, 1),
            bond(4, 16, 1),
            bond(5, 17, 1),
            bond(5, 18, 1),
            bond(5, 19, 1),
            bond(6, 20, 1),
        ],
        stereo_atoms: vec![(
            AtomId(1),
            vec![ligand(0), ligand(2), ligand(6), ligand(10)],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    })
}

fn corpus() -> Vec<CorpusCase> {
    let mut corpus = vec![
        CorpusCase {
            name: "ordinary_naphthalene",
            molecule: ordinary_naphthalene(),
        },
        CorpusCase {
            name: "disconnected_rings",
            molecule: disconnected_rings(),
        },
        CorpusCase {
            name: "overlay_heavy",
            molecule: overlay_heavy(),
        },
        CorpusCase {
            name: "overlay_unconstrained",
            molecule: overlay_unconstrained(),
        },
        CorpusCase {
            name: "aromatic_electron_counts",
            molecule: aromatic_electron_counts(),
        },
        CorpusCase {
            name: "symmetric_charge_sum",
            molecule: symmetric_charge_sum(),
        },
        CorpusCase {
            name: "asymmetric_charge_sum",
            molecule: asymmetric_charge_sum(),
        },
        CorpusCase {
            name: "tetrahedral_stereo",
            molecule: tetrahedral_stereo(),
        },
        CorpusCase {
            name: "cis_trans_stereo_bond",
            molecule: cis_trans_stereo_bond(),
        },
        CorpusCase {
            name: "mixed_atom_and_bond_stereo",
            molecule: mixed_atom_and_bond_stereo(),
        },
        CorpusCase {
            name: "frame_relative_stereo_constraint",
            molecule: frame_relative_stereo_constraint(),
        },
        CorpusCase {
            name: "meso_dichlorobutane",
            molecule: meso_dichlorobutane(),
        },
        CorpusCase {
            name: "para_stereo_trichloropentane",
            molecule: para_stereo_trichloropentane(),
        },
        CorpusCase {
            name: "para_stereo_cascade",
            molecule: para_stereo_cascade(),
        },
        CorpusCase {
            name: "symmetric_stereo_triisopropanolamine",
            molecule: symmetric_stereo_triisopropanolamine(),
        },
        CorpusCase {
            name: "symmetric_stereo_butanol",
            molecule: symmetric_stereo_butanol(),
        },
        CorpusCase {
            name: "symmetric_stereo_dimethylbutanol",
            molecule: symmetric_stereo_dimethylbutanol(),
        },
    ];
    corpus.extend(retained_scaling_corpus());
    corpus
}

fn retained_scaling_corpus() -> Vec<CorpusCase> {
    SCALING_CASES
        .map(|(name, source, expected_counts)| {
            let molecule = MoleculeDsl::from_edn_str(source)
                .unwrap_or_else(|error| panic!("benchmark case {name} must parse: {error}"))
                .into_parts()
                .0;
            assert_eq!(
                molecule_counts(&molecule),
                expected_counts,
                "benchmark case {name} changed shape"
            );
            assert!(
                !molecule.has_constraints()
                    && molecule
                        .atoms()
                        .iter()
                        .all(|atom| atom.attributes().constraints.is_empty())
                    && molecule
                        .bonds()
                        .iter()
                        .all(|bond| bond.attributes().constraints.is_empty()),
                "benchmark case {name} must remain feature-free"
            );
            CorpusCase { name, molecule }
        })
        .into()
}

fn inline_constraint_corpus() -> Vec<CorpusCase> {
    SCALING_CASES
        .map(|(name, source, _)| {
            let constrained = source.replace("\"H#", "\"H#v1#");
            let molecule = MoleculeDsl::from_edn_str(&constrained)
                .unwrap_or_else(|error| panic!("benchmark case {name} must parse: {error}"))
                .into_parts()
                .0;
            assert_eq!(
                molecule
                    .atoms()
                    .iter()
                    .filter(|atom| !atom.attributes().constraints.is_empty())
                    .count(),
                constrained.matches("\"H#v1#").count(),
                "benchmark case {name} must constrain every hydrogen"
            );
            CorpusCase { name, molecule }
        })
        .into()
}

const LEVELS: [IncidenceLevel; 3] = [
    IncidenceLevel::Topology,
    IncidenceLevel::Constitution,
    IncidenceLevel::Full,
];

fn level_name(level: IncidenceLevel) -> &'static str {
    match level {
        IncidenceLevel::Topology => "topology",
        IncidenceLevel::Constitution => "constitution",
        IncidenceLevel::Full => "full",
    }
}

fn graph_size(nodes: usize, edges: usize) -> String {
    format!("n{nodes}_e{edges}")
}

fn reverse_remapping(molecule: &Molecule) -> MoleculeRemapping {
    fn reverse<Id>(count: usize) -> Remapping<Id>
    where
        Id: Copy + Debug + Into<usize> + From<usize>,
    {
        let images = (0..count).rev().map(Id::from).collect::<Vec<_>>();
        Remapping::new(images).unwrap()
    }

    MoleculeRemapping::new(
        GraphRemapping::new(
            reverse::<NodeId>(molecule.atoms().count()),
            reverse::<EdgeId>(molecule.bonds().count()),
        ),
        reverse::<DativeBondId>(molecule.dative_bonds().count()),
        reverse::<AromaticSystemId>(molecule.aromatic_systems().count()),
        reverse::<MulticenterBondId>(molecule.multicenter_bonds().count()),
        reverse::<NoncovalentBondId>(molecule.noncovalent_bonds().count()),
        reverse::<StereoAtomId>(molecule.stereo_atoms().count()),
        reverse::<StereoBondId>(molecule.stereo_bonds().count()),
    )
}

fn atom_addition(lhs: Molecule, element: Element) -> Reaction {
    let id = AtomId(lhs.atoms().count() as u32);
    Reaction::new(
        lhs,
        Deltas::from_iter([Delta::Atom(AtomDelta::Add {
            id,
            attributes: AtomForm::from_element(element),
        })]),
    )
}

fn modification(lhs: Molecule) -> Reaction {
    let changed = |configuration: &StereoConfigurationForm| {
        let StereoConfigurationForm::Kinded(kind, StereoCoset::Lit(coset)) = configuration else {
            panic!("benchmark stereo configurations are literal");
        };
        StereoConfigurationForm::kinded(*kind, u32::from(*coset == 0))
    };
    let mut deltas = vec![Delta::Atom(AtomDelta::ModifyField {
        id: AtomId(0),
        change: AtomFieldChange::Charge {
            old: lhs.atom(AtomId(0)).attributes().charge.clone(),
            new: NumForm::Lit(1),
        },
    })];
    if let Some(stereo) = lhs.stereo_atoms().iter().next() {
        let old = stereo.attributes().configuration.clone();
        deltas.push(Delta::StereoAtom(StereoAtomDelta::ModifyField {
            id: stereo.id(),
            change: StereoAtomFieldChange::Configuration {
                new: changed(&old),
                old,
            },
        }));
    } else if let Some(stereo) = lhs.stereo_bonds().iter().next() {
        let old = stereo.attributes().configuration.clone();
        deltas.push(Delta::StereoBond(StereoBondDelta::ModifyField {
            id: stereo.id(),
            change: StereoBondFieldChange::Configuration {
                new: changed(&old),
                old,
            },
        }));
    }
    Reaction::new(lhs, Deltas::from_iter(deltas))
}

fn bench_incidence_construction(c: &mut Criterion) {
    let corpus = corpus();

    for level in LEVELS {
        let mut group = c.benchmark_group(format!(
            "canonicalize/incidence_construction/{}",
            level_name(level)
        ));
        for case in &corpus {
            let incidence = case.molecule.incidence_graph(level);
            let size = graph_size(
                incidence.graph().node_count(),
                incidence.graph().edge_count(),
            );
            group.bench_function(BenchmarkId::new(case.name, size), |b| {
                b.iter(|| black_box(&case.molecule).incidence_graph(level))
            });
        }
        group.finish();
    }
}

fn bench_remapping(c: &mut Criterion) {
    let corpus = corpus();
    let mut group = c.benchmark_group("canonicalize/remapping");

    for case in &corpus {
        let correspondence = reverse_remapping(&case.molecule);
        let counts = molecule_counts(&case.molecule)
            .into_iter()
            .map(|count| count.to_string())
            .collect::<Vec<_>>()
            .join("_");
        group.bench_function(BenchmarkId::new(case.name, counts), |b| {
            b.iter_batched(
                || case.molecule.clone(),
                |molecule| black_box(molecule.remap(black_box(&correspondence))),
                BatchSize::SmallInput,
            )
        });
    }

    group.finish();
}

fn bench_canonicalize(c: &mut Criterion) {
    let corpus = corpus();

    let context = CanonicalizeContext {
        automorphism_algorithm: ALGORITHM,
    };
    let mut group = c.benchmark_group("canonicalize/complete/canonicalize");
    for case in &corpus {
        let counts = molecule_counts(&case.molecule)
            .into_iter()
            .map(|count| count.to_string())
            .collect::<Vec<_>>()
            .join("_");
        group.bench_function(BenchmarkId::new(case.name, counts), |b| {
            b.iter_batched(
                || case.molecule.clone(),
                |molecule| {
                    black_box(
                        molecule
                            .canonicalize(&context)
                            .expect("benchmark corpus canonicalizes"),
                    )
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();
}

fn bench_retained_scaling_cases(c: &mut Criterion) {
    let context = CanonicalizeContext {
        automorphism_algorithm: ALGORITHM,
    };
    let cases = retained_scaling_corpus()
        .into_iter()
        .map(|case| {
            let renumbered = case
                .molecule
                .clone()
                .remap(&reverse_remapping(&case.molecule));
            let additional_atom = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::O)],
                ..Default::default()
            });
            let unequal = Molecule::combine_all([&case.molecule, &additional_atom]);
            assert!(case.molecule.canonical_eq(&renumbered, &context));
            assert!(!case.molecule.canonical_eq(&unequal, &context));
            (case, renumbered, unequal)
        })
        .collect::<Vec<_>>();

    let mut group = c.benchmark_group("canonicalize/scaling/canonicalize");
    for (case, _, _) in &cases {
        group.bench_function(case.name, |b| {
            b.iter_batched(
                || case.molecule.clone(),
                |molecule| {
                    black_box(
                        molecule
                            .canonicalize(&context)
                            .expect("retained scaling case canonicalizes"),
                    )
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();

    let mut group = c.benchmark_group("canonicalize/scaling/tracked_canonicalize");
    for (case, _, _) in &cases {
        group.bench_function(case.name, |b| {
            b.iter_batched(
                || case.molecule.clone(),
                |molecule| {
                    black_box(
                        molecule
                            .tracked_canonicalize(&context)
                            .expect("retained scaling case canonicalizes"),
                    )
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();

    let mut group = c.benchmark_group("canonicalize/scaling/canonical_hash");
    for (case, _, _) in &cases {
        group.bench_function(case.name, |b| {
            b.iter_batched(
                || case.molecule.clone(),
                |molecule| {
                    black_box(
                        molecule
                            .canonical_hash(&context)
                            .expect("retained scaling case canonicalizes"),
                    )
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();

    let mut group = c.benchmark_group("canonicalize/scaling/canonical_eq_equal");
    for (case, renumbered, _) in &cases {
        group.bench_function(case.name, |b| {
            b.iter(|| black_box(&case.molecule).canonical_eq(black_box(renumbered), &context))
        });
    }
    group.finish();

    let mut group = c.benchmark_group("canonicalize/scaling/canonical_eq_unequal");
    for (case, _, unequal) in &cases {
        group.bench_function(case.name, |b| {
            b.iter(|| black_box(&case.molecule).canonical_eq(black_box(unequal), &context))
        });
    }
    group.finish();
}

fn bench_topology_path_scaling(c: &mut Criterion) {
    let context = CanonicalizeContext {
        automorphism_algorithm: ALGORITHM,
    };
    let cases = [8, 16, 32, 64, 77, 128].map(|atom_count| (atom_count, carbon_path(atom_count)));
    let mut group = c.benchmark_group("canonicalize/scaling/topology_path");

    for (atom_count, molecule) in &cases {
        group.bench_function(BenchmarkId::from_parameter(atom_count), |b| {
            b.iter_batched(
                || molecule.clone(),
                |molecule| {
                    black_box(
                        molecule
                            .canonicalize(&context)
                            .expect("topology path canonicalizes"),
                    )
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();
}

fn bench_reframe(c: &mut Criterion) {
    let corpus = reframing_corpus();

    let mut group = c.benchmark_group("reframe/representative_action");
    for case in &corpus {
        group.bench_function(case.name, |b| {
            b.iter(|| black_box(&case.molecule).representative_action())
        });
    }
    group.finish();

    let mut group = c.benchmark_group("reframe/tracked_reframe");
    for case in &corpus {
        group.bench_function(case.name, |b| {
            b.iter_batched(
                || case.molecule.clone(),
                |molecule| {
                    black_box(
                        molecule
                            .tracked_reframe()
                            .expect("benchmark corpus reframes"),
                    )
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();

    let mut group = c.benchmark_group("reframe/reframe");
    for case in &corpus {
        group.bench_function(case.name, |b| {
            b.iter_batched(
                || case.molecule.clone(),
                |molecule| black_box(molecule.reframe().expect("benchmark corpus reframes")),
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();

    let representatives = corpus
        .iter()
        .map(|case| {
            case.molecule
                .clone()
                .reframe()
                .expect("benchmark corpus reframes")
        })
        .collect::<Vec<_>>();
    let mut group = c.benchmark_group("reframe/framed_eq");
    for (case, representative) in corpus.iter().zip(&representatives) {
        group.bench_function(case.name, |b| {
            b.iter(|| black_box(&case.molecule).framed_eq(black_box(representative)))
        });
    }
    group.finish();
}

fn bench_reaction_canonical_eq(c: &mut Criterion) {
    let context = CanonicalizeContext {
        automorphism_algorithm: ALGORITHM,
    };
    let reaction = |lhs: Molecule, element: Element| {
        let id = AtomId(lhs.atoms().count() as u32);
        Reaction::new(
            lhs,
            Deltas::from_iter([Delta::Atom(AtomDelta::Add {
                id,
                attributes: AtomForm::from_element(element),
            })]),
        )
    };
    let cases = corpus()
        .into_iter()
        .map(|case| {
            let source = reaction(case.molecule.clone(), Element::N);
            let renumbered = reaction(
                case.molecule
                    .clone()
                    .remap(&reverse_remapping(&case.molecule)),
                Element::N,
            );
            let unequal = reaction(case.molecule.clone(), Element::O);
            let spans = (
                source.to_reaction_span().unwrap(),
                renumbered.to_reaction_span().unwrap(),
            );
            (case.name, source, renumbered, unequal, spans)
        })
        .collect::<Vec<_>>();

    let mut group = c.benchmark_group("canonicalize/reaction/canonical_eq_equal");
    for (name, source, renumbered, _, _) in &cases {
        group.bench_function(*name, |b| {
            b.iter(|| black_box(source).canonical_eq(black_box(renumbered), &context))
        });
    }
    group.finish();

    let mut group = c.benchmark_group("canonicalize/reaction/canonical_eq_unequal");
    for (name, source, _, unequal, _) in &cases {
        group.bench_function(*name, |b| {
            b.iter(|| black_box(source).canonical_eq(black_box(unequal), &context))
        });
    }
    group.finish();

    let mut group = c.benchmark_group("canonicalize/reaction_span/canonical_eq_equal");
    for (name, _, _, _, (left, right)) in &cases {
        group.bench_function(*name, |b| {
            b.iter(|| black_box(left).canonical_eq(black_box(right), &context))
        });
    }
    group.finish();
}

fn bench_complete_comparison(c: &mut Criterion) {
    let context = CanonicalizeContext {
        automorphism_algorithm: ALGORITHM,
    };
    let cases = corpus()
        .into_iter()
        .map(|case| {
            let renumbered = case
                .molecule
                .clone()
                .remap(&reverse_remapping(&case.molecule));
            assert!(case.molecule.canonical_eq(&renumbered, &context));
            (case, renumbered)
        })
        .collect::<Vec<_>>();

    let mut group = c.benchmark_group("canonicalize/complete/canonical_hash");
    for (case, _) in &cases {
        group.bench_function(case.name, |b| {
            b.iter(|| {
                black_box(&case.molecule)
                    .canonical_hash(&context)
                    .expect("benchmark corpus canonicalizes")
            })
        });
    }
    group.finish();

    let mut group = c.benchmark_group("canonicalize/complete/canonical_eq_equal");
    for (case, renumbered) in &cases {
        group.bench_function(case.name, |b| {
            b.iter(|| black_box(&case.molecule).canonical_eq(black_box(renumbered), &context))
        });
    }
    group.finish();
}

fn bench_reaction_span(c: &mut Criterion) {
    let context = CanonicalizeContext {
        automorphism_algorithm: ALGORITHM,
    };
    let cases = corpus()
        .into_iter()
        .map(|case| {
            let remapping = reverse_remapping(&case.molecule);
            let added = atom_addition(case.molecule.clone(), Element::N)
                .to_reaction_span()
                .expect("benchmark addition materializes");
            let modified = modification(case.molecule)
                .to_reaction_span()
                .expect("benchmark modification materializes");
            let renumbered = modified.clone().remap(&remapping);
            assert!(modified.canonical_eq(&renumbered, &context));
            (case.name, added, modified, renumbered)
        })
        .collect::<Vec<_>>();

    let mut group = c.benchmark_group("canonicalize/reaction_span/canonicalize");
    for (name, added, _, _) in &cases {
        group.bench_function(*name, |b| {
            b.iter_batched(
                || added.clone(),
                |span| {
                    black_box(
                        span.canonicalize(&context)
                            .expect("benchmark span canonicalizes"),
                    )
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();

    let mut group = c.benchmark_group("canonicalize/reaction_span/canonical_hash");
    for (name, added, _, _) in &cases {
        group.bench_function(*name, |b| {
            b.iter(|| {
                black_box(added)
                    .canonical_hash(&context)
                    .expect("benchmark span canonicalizes")
            })
        });
    }
    group.finish();

    let mut group = c.benchmark_group("canonicalize/reaction_span_modified/canonicalize");
    for (name, _, modified, _) in &cases {
        group.bench_function(*name, |b| {
            b.iter_batched(
                || modified.clone(),
                |span| {
                    black_box(
                        span.canonicalize(&context)
                            .expect("benchmark span canonicalizes"),
                    )
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();

    let mut group = c.benchmark_group("canonicalize/reaction_span_modified/canonical_eq_equal");
    for (name, _, modified, renumbered) in &cases {
        group.bench_function(*name, |b| {
            b.iter(|| black_box(modified).canonical_eq(black_box(renumbered), &context))
        });
    }
    group.finish();
}

fn bench_inline_constraint(c: &mut Criterion) {
    let context = CanonicalizeContext {
        automorphism_algorithm: ALGORITHM,
    };
    let cases = inline_constraint_corpus()
        .into_iter()
        .map(|case| {
            let renumbered = case
                .molecule
                .clone()
                .remap(&reverse_remapping(&case.molecule));
            let span = atom_addition(case.molecule.clone(), Element::N)
                .to_reaction_span()
                .expect("benchmark addition materializes");
            let renumbered_span = atom_addition(renumbered.clone(), Element::N)
                .to_reaction_span()
                .expect("benchmark addition materializes");
            assert!(case.molecule.canonical_eq(&renumbered, &context));
            assert!(span.canonical_eq(&renumbered_span, &context));
            (case, renumbered, span, renumbered_span)
        })
        .collect::<Vec<_>>();

    let mut group = c.benchmark_group("canonicalize/inline_constraint/molecule/canonicalize");
    for (case, _, _, _) in &cases {
        group.bench_function(case.name, |b| {
            b.iter_batched(
                || case.molecule.clone(),
                |molecule| {
                    black_box(
                        molecule
                            .canonicalize(&context)
                            .expect("benchmark molecule canonicalizes"),
                    )
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();

    let mut group = c.benchmark_group("canonicalize/inline_constraint/molecule/canonical_eq_equal");
    for (case, renumbered, _, _) in &cases {
        group.bench_function(case.name, |b| {
            b.iter(|| black_box(&case.molecule).canonical_eq(black_box(renumbered), &context))
        });
    }
    group.finish();

    let mut group = c.benchmark_group("canonicalize/inline_constraint/reaction_span/canonicalize");
    for (case, _, span, _) in &cases {
        group.bench_function(case.name, |b| {
            b.iter_batched(
                || span.clone(),
                |span| {
                    black_box(
                        span.canonicalize(&context)
                            .expect("benchmark span canonicalizes"),
                    )
                },
                BatchSize::SmallInput,
            )
        });
    }
    group.finish();

    let mut group =
        c.benchmark_group("canonicalize/inline_constraint/reaction_span/canonical_eq_equal");
    for (case, _, span, renumbered_span) in &cases {
        group.bench_function(case.name, |b| {
            b.iter(|| black_box(span).canonical_eq(black_box(renumbered_span), &context))
        });
    }
    group.finish();
}

fn molecule_counts(molecule: &Molecule) -> [usize; 8] {
    [
        molecule.atoms().count(),
        molecule.bonds().count(),
        molecule.dative_bonds().count(),
        molecule.aromatic_systems().count(),
        molecule.multicenter_bonds().count(),
        molecule.noncovalent_bonds().count(),
        molecule.stereo_atoms().count(),
        molecule.stereo_bonds().count(),
    ]
}

criterion_group!(
    canonicalize,
    bench_incidence_construction,
    bench_remapping,
    bench_canonicalize,
    bench_retained_scaling_cases,
    bench_topology_path_scaling,
    bench_reframe,
    bench_reaction_canonical_eq,
    bench_complete_comparison,
    bench_reaction_span,
    bench_inline_constraint,
);
criterion_main!(canonicalize);
