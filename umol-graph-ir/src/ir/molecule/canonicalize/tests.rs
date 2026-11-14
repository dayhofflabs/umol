//! Unit and bounded-exhaustive tests for molecule canonicalization.

use rstest::{fixture, rstest};
use umol_chem::element::Element;
use umol_edn::FromEdn;
use umol_graph_core::AutomorphismAlgorithm;
use umol_perm::Permutation;

use super::super::symmetry::tests::initial_color_molecule;
use super::*;
use crate::ir::canonicalize::tests::{
    canonicalize_context, explicitly_dense_equivalent, permutations,
};
use crate::ir::symmetry::{same_cells, AutomorphismAdapterOutput, CanonicalSearchStats};
use crate::ir::{
    AromaticSystemForm, AromaticSystemId, AtomConstraintForm, AtomForm, AtomId, BondForm, BondId,
    Constraint, Contradiction, DativeBondForm, DativeBondId, ElementForm, Entity, Incidence,
    IncidenceLevel, MoleculeConstraint, MoleculeCorrespondence, MoleculeEntries,
    MulticenterBondForm, MulticenterBondId, NoncovalentBondForm, NoncovalentBondId,
    NoncovalentBondKind, NoncovalentBondKindForm, NumForm, StereoAtomConstraintForm,
    StereoAtomForm, StereoAtomId, StereoBondForm, StereoBondId, StereoConfigurationForm,
    StereoCoset, StereoKind, StereoLigand, StereoLigandKind, StereoLigandPair, Topicity,
    TopicityForm, TopicityRelationForm,
};
use crate::mol_dsl;

#[rstest]
#[case::element_order(
    vec![AtomForm::from_element(Element::N), AtomForm::from_element(Element::C)],
    vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)]
)]
#[case::raw_element(
    vec![AtomForm::from_element(Element::N), AtomForm::new(ElementForm::lit_set([Element::C]))],
    vec![AtomForm::from_element(Element::C), AtomForm::from_element(Element::N)]
)]
#[case::raw_constraint(
    vec![AtomForm::from_element(Element::C)
        .with_constraint(AtomConstraintForm::Valence(NumForm::lit_set([4])))],
    vec![AtomForm::from_element(Element::C)
        .with_constraint(AtomConstraintForm::Valence(NumForm::Lit(4)))]
)]
#[case::constraint_order(
    vec![
        AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4)),
        AtomForm::from_element(Element::C),
    ],
    vec![
        AtomForm::from_element(Element::C),
        AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4)),
    ]
)]
fn test_molecule_canonicalize_colors(
    canonicalize_context: CanonicalizeContext,
    #[case] atoms: Vec<AtomForm>,
    #[case] expected_atoms: Vec<AtomForm>,
) {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms,
        ..Default::default()
    });
    let expected = Molecule::from_entries(MoleculeEntries {
        atoms: expected_atoms,
        ..Default::default()
    });
    assert_eq!(
        molecule.clone().canonicalize(&canonicalize_context),
        Ok(expected.clone())
    );
    assert!(molecule.canonical_eq(&expected, &canonicalize_context));
}

#[fixture]
fn stereo_atom_canonicalization_molecule() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        )],
        ..Default::default()
    })
}

#[fixture]
fn stereo_bond_canonicalization_molecule() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(1), AtomId(4), BondForm::from_order(1)),
            (AtomId(1), AtomId(5), BondForm::from_order(1)),
        ],
        stereo_bonds: vec![(
            BondId(0),
            vec![
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            StereoBondForm::new(StereoKind::CisTrans, 1u32),
        )],
        ..Default::default()
    })
}

#[fixture]
fn stereo_constraint_canonicalization_molecule() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::H),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        )],
        ..Default::default()
    })
}

/// Eight outer stereo atoms in two mirrored groups of four, each site bearing a distinct stereo
/// descriptor so partition refinement can separate them. Only three stereo kinds are admissible on
/// an atom at degree 4, so the fourth descriptor is a second `Tetrahedral` coset rather than a
/// fourth kind.
#[fixture]
fn para_stereo_canonicalization_molecule() -> Molecule {
    let outer_ligands = vec![
        StereoLigand::new(AtomId(10), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(11), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(12), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(13), StereoLigandKind::Atom),
    ];

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
        .map(AtomForm::from_element)
        .collect(),
        bonds: (2..=5)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
            .chain((6..=9).map(|ligand| (AtomId(1), AtomId(ligand), BondForm::from_order(1))))
            .chain((2..=9).flat_map(|site| {
                (10..=13).map(move |ligand| (AtomId(site), AtomId(ligand), BondForm::from_order(1)))
            }))
            .collect(),
        stereo_atoms: vec![
            (
                AtomId(0),
                [2, 3, 4, 5]
                    .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                    .into(),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(1),
                [6, 8, 7, 9]
                    .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                    .into(),
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

fn reverse_remapping(molecule: &Molecule) -> MoleculeRemapping {
    let images = molecule
        .entity_counts()
        .map(|count| (0..count).rev().collect::<Vec<_>>());
    molecule_remapping(&images)
}

impl AutomorphismAdapter {
    pub(super) fn automorphisms(
        &self,
        algorithm: AutomorphismAlgorithm,
    ) -> AutomorphismAdapterOutput {
        let output = self
            .graph()
            .automorphisms(|node| self.color(node), algorithm);
        self.project_automorphisms(&output)
    }
}

#[fixture]
fn symmetric_stereo_canonicalization_molecule() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 5],
        bonds: (1..=4)
            .map(|id| (AtomId(0), AtomId(id), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
        )],
        ..Default::default()
    })
}

#[fixture]
fn stereo_site_exchange_canonicalization_molecule() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: [
            Element::C,
            Element::F,
            Element::Cl,
            Element::Br,
            Element::I,
            Element::C,
            Element::F,
            Element::Cl,
            Element::Br,
            Element::I,
        ]
        .into_iter()
        .map(AtomForm::from_element)
        .collect(),
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
            (AtomId(5), AtomId(6), BondForm::from_order(1)),
            (AtomId(5), AtomId(7), BondForm::from_order(1)),
            (AtomId(5), AtomId(8), BondForm::from_order(1)),
            (AtomId(5), AtomId(9), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                (1..=4)
                    .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(5),
                (6..=9)
                    .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
        ],
        ..Default::default()
    })
}

#[fixture]
fn meso_canonicalization_molecule() -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Cl),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
            (AtomId(1), AtomId(3), BondForm::from_order(1)),
            (AtomId(1), AtomId(5), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
            ),
        ],
        ..Default::default()
    })
}

#[rstest]
#[case::topology(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::C),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    }),
    DescriptionLevel::Topology,
)]
#[case::constitution(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::B),
        ],
        dative: vec![(vec![AtomId(0)], AtomId(1), DativeBondForm::from_order(1))],
        ..Default::default()
    }),
    DescriptionLevel::Constitution,
)]
#[case::structure(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::H),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
        ],
        bonds: (1..=4)
            .map(|id| (AtomId(0), AtomId(id), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 0_u32),
        )],
        ..Default::default()
    }),
    DescriptionLevel::Structure,
)]
#[case::full(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: None,
            sum: NumForm::Lit(0),
        })
        .into(),
        ..Default::default()
    }),
    DescriptionLevel::Full,
)]
#[case::atom_constraint(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_constraint(AtomConstraintForm::valence(4)),
            AtomForm::from_element(Element::C),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Topology,
)]
#[case::aromatic_system(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::default(),
        )],
        ..Default::default()
    }),
    DescriptionLevel::Constitution,
)]
#[case::multicenter_bond(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::B); 3],
        multicenter: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondForm::default(),
        )],
        ..Default::default()
    }),
    DescriptionLevel::Constitution,
)]
#[case::noncovalent_bond(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::H),
        ],
        noncovalent: vec![(
            [AtomId(0), AtomId(1)],
            NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
        )],
        ..Default::default()
    }),
    DescriptionLevel::Constitution,
)]
#[case::stereo_bond(stereo_bond_canonicalization_molecule(), DescriptionLevel::Structure)]
fn test_molecule_canonicalize_level(
    #[case] molecule: Molecule,
    #[case] expected: DescriptionLevel,
) {
    assert_eq!(molecule_canonicalize_level(&molecule), expected);
}

#[rstest]
fn test_canonicalize_structure_para_stereo(
    para_stereo_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let atom_images = [1, 0, 6, 8, 7, 9, 2, 4, 3, 5, 13, 11, 10, 12].map(AtomId);
    let mut images = para_stereo_canonicalization_molecule
        .entity_counts()
        .map(|count| (0..count).collect::<Vec<_>>());
    images[0] = atom_images.iter().map(|id| id.index()).collect();
    let remapping = molecule_remapping(&images);
    let renumbered = para_stereo_canonicalization_molecule
        .clone()
        .remap(&remapping);
    let context = canonicalize_context;
    let (canonical, remapping) =
        canonicalize_structure(para_stereo_canonicalization_molecule.clone(), &context)
            .expect("fixed molecule canonicalizes");

    assert_eq!(
        para_stereo_canonicalization_molecule
            .remap(&remapping)
            .reframe(),
        Ok(canonical.clone()),
    );
    assert_eq!(
        renumbered
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &context)
            .map(|(canonical, _)| canonical),
        Ok(canonical.clone()),
    );
    assert_eq!(
        canonical
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &context)
            .map(|(canonical, _)| canonical),
        Ok(canonical)
    );
}

#[rstest]
fn test_canonicalize_structure_stereo_atom(
    canonicalize_context: CanonicalizeContext,
    stereo_atom_canonicalization_molecule: Molecule,
) {
    let renumbering = molecule_remapping(&[
        vec![4, 2, 0, 3, 1],
        vec![2, 0, 3, 1],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        vec![0],
        Vec::new(),
    ]);
    let renumbered = stereo_atom_canonicalization_molecule
        .clone()
        .remap(&renumbering);
    let canonical = stereo_atom_canonicalization_molecule
        .clone()
        .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
        .map(|(canonical, _)| canonical)
        .expect("fixed molecule canonicalizes");

    assert_eq!(
        renumbered
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        Ok(canonical.clone()),
    );
    assert_eq!(
        canonical
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        Ok(canonical.clone()),
    );
    assert_eq!(canonical.check_integrity(), Ok(()));
    assert_eq!(canonical.stereo_atoms().count(), 1);
}

#[rstest]
fn test_canonicalize_structure_configuration(
    canonicalize_context: CanonicalizeContext,
    stereo_atom_canonicalization_molecule: Molecule,
) {
    let mut opposite = stereo_atom_canonicalization_molecule.clone();
    opposite
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32);

    assert_ne!(
        stereo_atom_canonicalization_molecule
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        opposite
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
    );
}

#[rstest]
fn test_canonicalize_structure_stereo_bond(
    canonicalize_context: CanonicalizeContext,
    stereo_bond_canonicalization_molecule: Molecule,
) {
    let renumbering = molecule_remapping(&[
        vec![5, 3, 1, 4, 2, 0],
        vec![4, 2, 0, 3, 1],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        vec![0],
    ]);
    let renumbered = stereo_bond_canonicalization_molecule
        .clone()
        .remap(&renumbering);
    let canonical = stereo_bond_canonicalization_molecule
        .clone()
        .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
        .map(|(canonical, _)| canonical)
        .expect("fixed molecule canonicalizes");

    assert_eq!(
        renumbered
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        Ok(canonical.clone()),
    );
    assert_eq!(
        canonical
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        Ok(canonical.clone()),
    );
    assert_eq!(canonical.check_integrity(), Ok(()));
    assert_eq!(canonical.stereo_bonds().count(), 1);
}

#[rstest]
#[case::stereo_atom(stereo_atom_canonicalization_molecule(), 1)]
#[case::stereo_bond(stereo_bond_canonicalization_molecule(), 1)]
#[case::symmetric_stereo(symmetric_stereo_canonicalization_molecule(), 24)]
#[case::para_stereo(para_stereo_canonicalization_molecule(), 8)]
fn test_canonicalize_structure_visited_leaves(
    canonicalize_context: CanonicalizeContext,
    #[case] molecule: Molecule,
    #[case] expected_visited_leaves: usize,
) {
    let molecule = molecule.normalize().expect("fixed molecule normalizes");
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(
        &molecule,
        &incidence_graph,
        DescriptionLevel::Structure,
        true,
    )
    .expect("fixed molecule has colors");
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let leaf_candidate = |order: &[NodeId]| {
        let entity_order = entity_order(
            &molecule,
            &incidence_graph,
            DescriptionLevel::Structure,
            order,
        );
        CanonicalCandidate {
            key: color_key(&molecule, &incidence_graph, &colors, &entity_order),
            entity_order,
        }
    };
    let filter_generators = |generators: &mut Vec<Vec<NodeId>>| {
        retain_stereo_preserving_generators(&molecule, &incidence_graph, generators);
    };
    let search = canonical_search_with_generator_filter(
        &adapter,
        &adapter.colors,
        canonicalize_context.automorphism_algorithm,
        &leaf_candidate,
        &filter_generators,
    );

    assert_eq!(search.stats.visited_leaves, expected_visited_leaves);
}

#[rstest]
fn test_canonicalize_structure_renumbering(
    symmetric_stereo_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let canonical = symmetric_stereo_canonicalization_molecule
        .clone()
        .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
        .map(|(canonical, _)| canonical)
        .expect("fixed molecule canonicalizes");

    for rank in 0..(1..=5).product() {
        let permutation = Permutation::unrank(5, rank);
        let atom_images = (0..5)
            .map(|index| AtomId(permutation.apply(index) as u32))
            .collect::<Vec<_>>();
        let remapping = molecule_remapping(&[
            atom_images.iter().map(|id| id.index()).collect(),
            vec![0, 1, 2, 3],
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            vec![0],
            Vec::new(),
        ]);
        let renumbered = symmetric_stereo_canonicalization_molecule
            .clone()
            .remap(&remapping);
        let (actual, canonical_remapping) =
            canonicalize_structure(renumbered.clone(), &canonicalize_context)
                .expect("renumbered molecule canonicalizes");

        assert_eq!(actual, canonical, "rank {rank}");
        assert!(
            renumbered.framed_eq_under(&actual, &canonical_remapping),
            "rank {rank}",
        );
        assert_eq!(
            renumbered.remap(&canonical_remapping).reframe(),
            Ok(actual),
            "rank {rank}",
        );
    }
}

#[rstest]
#[case::nauty(AutomorphismAlgorithm::Nauty)]
fn test_canonicalize_structure_orbit_pruning(
    stereo_site_exchange_canonicalization_molecule: Molecule,
    #[case] algorithm: AutomorphismAlgorithm,
) {
    let molecule = stereo_site_exchange_canonicalization_molecule
        .normalize()
        .expect("fixed molecule normalizes");
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(
        &molecule,
        &incidence_graph,
        DescriptionLevel::Structure,
        true,
    )
    .expect("fixed molecule has colors");
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let leaf_candidate = |order: &[NodeId]| {
        let entity_order = entity_order(
            &molecule,
            &incidence_graph,
            DescriptionLevel::Structure,
            order,
        );
        CanonicalCandidate {
            key: color_key(&molecule, &incidence_graph, &colors, &entity_order),
            entity_order,
        }
    };
    let filter_generators = |generators: &mut Vec<Vec<NodeId>>| {
        retain_stereo_preserving_generators(&molecule, &incidence_graph, generators);
    };
    let actual = canonical_search_with_generator_filter(
        &adapter,
        &adapter.colors,
        algorithm,
        &leaf_candidate,
        &filter_generators,
    );
    let remapping =
        remapping_from_order(&molecule, &incidence_graph, &actual.candidate.entity_order);
    let canonical = molecule
        .clone()
        .remap(&remapping)
        .reframe()
        .expect("selected leaf preserves compatible frames");

    assert!(molecule.framed_eq_under(&canonical, &remapping));
    assert_eq!(
        actual.stats,
        CanonicalSearchStats {
            initial_residual_cell_sizes: vec![2; 10],
            refinement_calls: 2,
            branch_order_calls: 1,
            backend_calls: 1,
            visited_leaves: 1,
            leaf_comparisons: 0,
            orbit_pruned_branches: 1,
        },
    );
}

#[rstest]
fn test_canonicalize_structure_meso(
    meso_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let remapping = molecule_remapping(&[
        vec![1, 0, 3, 2, 5, 4],
        vec![0, 3, 4, 1, 2],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        vec![1, 0],
        Vec::new(),
    ]);
    let renumbered = meso_canonicalization_molecule.clone().remap(&remapping);

    assert_eq!(
        renumbered
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        meso_canonicalization_molecule
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
    );
}

#[rstest]
#[case::kinded(StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32))]
#[case::undetermined(StereoConfigurationForm::Undetermined)]
fn test_canonicalize_structure_stereo_frame(
    stereo_constraint_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
    #[case] configuration: StereoConfigurationForm,
) {
    let mut source = stereo_constraint_canonicalization_molecule;
    source
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = configuration;
    let reframed = source
        .clone()
        .reframe()
        .expect("integrity-valid molecule reframes");

    assert_eq!(
        reframed
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        source
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Structure, &canonicalize_context)
            .map(|(canonical, _)| canonical),
    );
}

#[rstest]
#[case::atom_constraint(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(1)), false)]
#[case::atom_element(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(2)), false)]
#[case::atom_charge(Entity::Atom(AtomId(0)), Entity::Atom(AtomId(3)), false)]
#[case::normalized_bond(Entity::Bond(BondId(0)), Entity::Bond(BondId(1)), true)]
#[case::bond_order(Entity::Bond(BondId(0)), Entity::Bond(BondId(2)), false)]
#[case::normalized_dative(
    Entity::DativeBond(DativeBondId(0)),
    Entity::DativeBond(DativeBondId(1)),
    true
)]
#[case::dative_order(
    Entity::DativeBond(DativeBondId(0)),
    Entity::DativeBond(DativeBondId(2)),
    false
)]
#[case::aromatic_electrons(
    Entity::AromaticSystem(AromaticSystemId(0)),
    Entity::AromaticSystem(AromaticSystemId(1)),
    true
)]
#[case::aromatic_charge(
    Entity::AromaticSystem(AromaticSystemId(0)),
    Entity::AromaticSystem(AromaticSystemId(2)),
    false
)]
#[case::multicenter_electrons(
    Entity::MulticenterBond(MulticenterBondId(0)),
    Entity::MulticenterBond(MulticenterBondId(1)),
    true
)]
#[case::noncovalent_kind_equal(
    Entity::NoncovalentBond(NoncovalentBondId(0)),
    Entity::NoncovalentBond(NoncovalentBondId(1)),
    true
)]
#[case::noncovalent_kind_distinct(
    Entity::NoncovalentBond(NoncovalentBondId(0)),
    Entity::NoncovalentBond(NoncovalentBondId(2)),
    false
)]
#[case::stereo_atom_configuration(
    Entity::StereoAtom(StereoAtomId(0)),
    Entity::StereoAtom(StereoAtomId(1)),
    true
)]
#[case::stereo_atom_kind(
    Entity::StereoAtom(StereoAtomId(0)),
    Entity::StereoAtom(StereoAtomId(2)),
    false
)]
#[case::stereo_bond_configuration(
    Entity::StereoBond(StereoBondId(0)),
    Entity::StereoBond(StereoBondId(1)),
    true
)]
#[case::entity_kind(Entity::Atom(AtomId(0)), Entity::Bond(BondId(0)), false)]
fn test_same_entity_color(
    initial_color_molecule: Molecule,
    #[case] entity: Entity,
    #[case] image: Entity,
    #[case] expected: bool,
) {
    let molecule = initial_color_molecule
        .normalize()
        .expect("fixed molecule normalizes");

    assert_eq!(
        same_entity_color(&molecule, entity, &molecule, image),
        expected
    );
}

#[rstest]
#[case::renumbered(
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [1 5 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[5 4 "1"] [4 3 "1"] [3 2 "1"] [2 1 "1"] [4 0 "1"]]}"#),
    true
)]
#[case::different_element(
    mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "C" "N"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
    false
)]
#[case::different_element_counts(
    mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "O" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
    false
)]
#[case::different_bond_order(
    mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "2"]]}"#),
    false
)]
#[case::electron_counts_moved(
    mol_dsl!(
        r#"{:atoms ["O" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [0 4 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["O" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [0 4 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,2,1,1,1]"}]}"#
    ),
    true
)]
#[case::different_electron_counts(
    mol_dsl!(
        r#"{:atoms ["O" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [0 4 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["O" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [0 4 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,1]"}]}"#
    ),
    false
)]
fn test_same_colors(#[case] source: Molecule, #[case] target: Molecule, #[case] expected: bool) {
    let source = source.normalize().unwrap();
    let target = target.normalize().unwrap();
    let source_incidence = source.incidence_graph(IncidenceLevel::Constitution);
    let target_incidence = target.incidence_graph(IncidenceLevel::Constitution);
    let source_colors = initial_colors(
        &source,
        &source_incidence,
        DescriptionLevel::Constitution,
        true,
    )
    .unwrap();
    let target_colors = initial_colors(
        &target,
        &target_incidence,
        DescriptionLevel::Constitution,
        true,
    )
    .unwrap();

    assert_eq!(
        same_colors(
            &source,
            &source_incidence,
            &source_colors,
            &target,
            &target_incidence,
            &target_colors,
        ),
        expected,
    );
}

#[rstest]
#[case::incomplete_order_one(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    }),
    Graph::new(2, &[[0, 1]]),
    vec![AutomorphismColor::Entity(0); 2],
    vec![
        SubdivisionNodeSource::Node(NodeId(0)),
        SubdivisionNodeSource::Node(NodeId(1)),
    ],
    vec![vec![NodeId(0), NodeId(1)]],
)]
#[case::modal_majority(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2).into_concrete()),
            (AtomId(1), AtomId(2), BondForm::from_order(2).into_concrete()),
            (AtomId(0), AtomId(2), BondForm::from_order(1).into_concrete()),
        ],
        ..Default::default()
    }),
    Graph::new(4, &[[0, 1], [1, 2], [0, 3], [3, 2]]),
    vec![
        AutomorphismColor::Entity(0),
        AutomorphismColor::Entity(0),
        AutomorphismColor::Entity(0),
        AutomorphismColor::Entity(1),
    ],
    vec![
        SubdivisionNodeSource::Node(NodeId(0)),
        SubdivisionNodeSource::Node(NodeId(1)),
        SubdivisionNodeSource::Node(NodeId(2)),
        SubdivisionNodeSource::Node(NodeId(3)),
    ],
    vec![vec![NodeId(0), NodeId(1), NodeId(2)], vec![NodeId(3)]],
)]
#[case::modal_tie(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2).into_concrete()),
            (AtomId(1), AtomId(2), BondForm::from_order(1).into_concrete()),
        ],
        ..Default::default()
    }),
    Graph::new(4, &[[0, 3], [3, 1], [1, 2]]),
    vec![
        AutomorphismColor::Entity(0),
        AutomorphismColor::Entity(0),
        AutomorphismColor::Entity(0),
        AutomorphismColor::Entity(2),
    ],
    vec![
        SubdivisionNodeSource::Node(NodeId(0)),
        SubdivisionNodeSource::Node(NodeId(1)),
        SubdivisionNodeSource::Node(NodeId(2)),
        SubdivisionNodeSource::Node(NodeId(3)),
    ],
    vec![vec![NodeId(0), NodeId(1), NodeId(2)], vec![NodeId(3)]],
)]
fn test_topology_search_input(
    #[case] molecule: Molecule,
    #[case] expected_graph: Graph,
    #[case] expected_colors: Vec<AutomorphismColor>,
    #[case] expected_node_sources: Vec<SubdivisionNodeSource>,
    #[case] expected_partition: Vec<Vec<NodeId>>,
) {
    let molecule = molecule.normalize().unwrap();
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Topology);
    let colors = initial_colors(
        &molecule,
        &incidence_graph,
        DescriptionLevel::Topology,
        true,
    )
    .unwrap();
    let input = topology_search_input(&molecule, &incidence_graph, &colors.entities);

    assert_eq!(input.adapter.graph, expected_graph);
    assert_eq!(input.adapter.colors, expected_colors);
    assert_eq!(input.adapter.node_sources, expected_node_sources);
    assert_eq!(
        input.adapter.source_node_count,
        incidence_graph.entity_count(EntityKind::Atom),
    );
    assert_eq!(input.initial_partition.cells, expected_partition);
}

#[rstest]
#[case::dative(Entity::DativeBond(DativeBondId(0)))]
#[case::aromatic(Entity::AromaticSystem(AromaticSystemId(0)))]
#[case::multicenter(Entity::MulticenterBond(MulticenterBondId(0)))]
#[case::noncovalent(Entity::NoncovalentBond(NoncovalentBondId(0)))]
#[case::stereo_atom(Entity::StereoAtom(StereoAtomId(0)))]
#[case::stereo_bond(Entity::StereoBond(StereoBondId(0)))]
fn test_remapping_from_order(initial_color_molecule: Molecule, #[case] excluded: Entity) {
    let incidence_graph = initial_color_molecule.incidence_graph(IncidenceLevel::Topology);
    let order = incidence_graph.graph().node_ids().collect::<Vec<_>>();
    let remapping = remapping_from_order(&initial_color_molecule, &incidence_graph, &order);

    assert_eq!(
        MoleculeCorrespondence::from(&remapping).right_of(excluded),
        Some(excluded)
    );
}

#[rstest]
#[case::identity(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
        ],
        ..Default::default()
    }),
    [0, 1, 2, 3, 4].map(NodeId).into(),
    molecule_remapping(&[
        vec![0, 1, 2],
        vec![0, 1],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    ])
)]
#[case::path_reversal(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
        ],
        ..Default::default()
    }),
    [2, 1, 0, 4, 3].map(NodeId).into(),
    molecule_remapping(&[
        vec![2, 1, 0],
        vec![1, 0],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    ])
)]
#[case::stereo_site_exchange(
    stereo_site_exchange_canonicalization_molecule(),
    [5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 14, 15, 16, 17, 10, 11, 12, 13, 19, 18].map(NodeId).into(),
    molecule_remapping(&[
        vec![5, 6, 7, 8, 9, 0, 1, 2, 3, 4],
        vec![4, 5, 6, 7, 0, 1, 2, 3],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        vec![1, 0],
        Vec::new(),
    ])
)]
fn test_remapping_from_generator(
    #[case] molecule: Molecule,
    #[case] generator: Vec<NodeId>,
    #[case] expected: MoleculeRemapping,
) {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);

    assert_eq!(
        remapping_from_generator(&molecule, &incidence_graph, &generator),
        expected
    );
}

#[rstest]
#[case::bond(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Topology,
    vec![
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::Bond(BondId(1)),
        Entity::Bond(BondId(0)),
    ],
)]
#[case::dative_bond(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 4],
        dative: vec![
            (vec![AtomId(0)], AtomId(2), DativeBondForm::from_order(1)),
            (vec![AtomId(1)], AtomId(3), DativeBondForm::from_order(1)),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Full,
    vec![
        Entity::Atom(AtomId(3)),
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::DativeBond(DativeBondId(1)),
        Entity::DativeBond(DativeBondId(0)),
    ],
)]
#[case::aromatic_system(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 4],
        aromatic: vec![
            (vec![AtomId(0), AtomId(1)], AromaticSystemForm::default()),
            (vec![AtomId(2), AtomId(3)], AromaticSystemForm::default()),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Full,
    vec![
        Entity::Atom(AtomId(3)),
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::AromaticSystem(AromaticSystemId(1)),
        Entity::AromaticSystem(AromaticSystemId(0)),
    ],
)]
#[case::multicenter_bond(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 4],
        multicenter: vec![
            (vec![AtomId(0), AtomId(1)], MulticenterBondForm::default()),
            (vec![AtomId(2), AtomId(3)], MulticenterBondForm::default()),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Full,
    vec![
        Entity::Atom(AtomId(3)),
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::MulticenterBond(MulticenterBondId(1)),
        Entity::MulticenterBond(MulticenterBondId(0)),
    ],
)]
#[case::noncovalent_bond(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 4],
        noncovalent: vec![
            ([AtomId(0), AtomId(1)], NoncovalentBondForm::default()),
            ([AtomId(2), AtomId(3)], NoncovalentBondForm::default()),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Full,
    vec![
        Entity::Atom(AtomId(3)),
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::NoncovalentBond(NoncovalentBondId(1)),
        Entity::NoncovalentBond(NoncovalentBondId(0)),
    ],
)]
#[case::stereo_atom(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        stereo_atoms: vec![
            (AtomId(0), Vec::new(), StereoAtomForm::default()),
            (AtomId(1), Vec::new(), StereoAtomForm::default()),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Full,
    vec![
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::StereoAtom(StereoAtomId(1)),
        Entity::StereoAtom(StereoAtomId(0)),
    ],
)]
#[case::stereo_bond(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 4],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2)),
            (AtomId(2), AtomId(3), BondForm::from_order(2)),
        ],
        stereo_bonds: vec![
            (
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                ],
                StereoBondForm::default(),
            ),
            (
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(3), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::LonePair),
                ],
                StereoBondForm::default(),
            ),
        ],
        ..Default::default()
    }),
    DescriptionLevel::Full,
    vec![
        Entity::Atom(AtomId(3)),
        Entity::Atom(AtomId(2)),
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::Bond(BondId(1)),
        Entity::Bond(BondId(0)),
        Entity::StereoBond(StereoBondId(1)),
        Entity::StereoBond(StereoBondId(0)),
    ],
)]
#[case::constitution_level(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 2],
        noncovalent: vec![([AtomId(0), AtomId(1)], NoncovalentBondForm::default())],
        stereo_atoms: vec![(AtomId(0), Vec::new(), StereoAtomForm::default())],
        ..Default::default()
    }),
    DescriptionLevel::Constitution,
    vec![
        Entity::Atom(AtomId(1)),
        Entity::Atom(AtomId(0)),
        Entity::NoncovalentBond(NoncovalentBondId(0)),
    ],
)]
fn test_entity_order(
    #[case] molecule: Molecule,
    #[case] level: DescriptionLevel,
    #[case] expected: Vec<Entity>,
) {
    let incidence_graph = molecule.incidence_graph(match level {
        DescriptionLevel::Topology => IncidenceLevel::Topology,
        DescriptionLevel::Constitution => IncidenceLevel::Constitution,
        DescriptionLevel::Structure | DescriptionLevel::Full => IncidenceLevel::Full,
    });
    let mut order = incidence_graph.graph().node_ids().collect::<Vec<_>>();
    order.reverse();

    assert_eq!(
        entity_order(&molecule, &incidence_graph, level, &order)
            .into_iter()
            .map(|node| incidence_graph.entity(node))
            .collect::<Vec<_>>(),
        expected,
    );
}

#[rstest]
#[case::localized_bond(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 2],
            bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
            ..Default::default()
        }),
        vec![Incidence::BondEndpoint, Incidence::BondEndpoint],
        0,
    )]
#[case::mixed_virtual_ligand_anchor(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![
                AtomForm::from_element(Element::C),
                AtomForm::from_element(Element::F),
                AtomForm::from_element(Element::Cl),
            ],
            bonds: vec![
                (AtomId(0), AtomId(1), BondForm::from_order(1)),
                (AtomId(0), AtomId(2), BondForm::from_order(1)),
            ],
            stereo_atoms: vec![(
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            )],
            ..Default::default()
        }),
        vec![
            Incidence::BondEndpoint,
            Incidence::BondEndpoint,
            Incidence::BondEndpoint,
            Incidence::BondEndpoint,
            Incidence::StereoSite,
            Incidence::StereoLigand(StereoLigandKind::Atom),
            Incidence::StereoLigand(StereoLigandKind::Atom),
            Incidence::StereoLigand(StereoLigandKind::ImplicitHydrogen),
            Incidence::StereoLigand(StereoLigandKind::LonePair),
        ],
        5,
    )]
fn test_automorphism_adapter_new(
    #[case] molecule: Molecule,
    #[case] expected_incidences: Vec<Incidence>,
    #[case] expected_occurrence_nodes: usize,
) {
    let molecule = molecule.normalize().unwrap();
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(
        &molecule,
        &incidence_graph,
        DescriptionLevel::Structure,
        true,
    )
    .unwrap();
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let source = incidence_graph.graph();

    assert_eq!(
        incidence_graph
            .incidences()
            .map(|(_, incidence)| incidence.clone())
            .collect::<Vec<_>>(),
        expected_incidences,
    );
    assert_eq!(
        adapter.graph().node_count(),
        source.node_count() + expected_occurrence_nodes,
    );
    assert_eq!(
        adapter.graph().edge_count(),
        source.edge_count() + expected_occurrence_nodes,
    );
    assert!(adapter.graph().is_simple());

    for node in source.node_ids() {
        assert_eq!(adapter.node_source(node), SubdivisionNodeSource::Node(node));
        assert_eq!(
            adapter.color(node),
            AutomorphismColor::Entity(colors.entities[node.index()]),
        );
    }
    for adapter_node in adapter.graph().node_ids().skip(source.node_count()) {
        let SubdivisionNodeSource::Edge(edge) = adapter.node_source(adapter_node) else {
            panic!("a node after the source nodes subdivides an incidence");
        };
        assert_eq!(
            adapter.color(adapter_node),
            AutomorphismColor::Incidence(colors.incidences[edge.index()]),
        );
        assert_eq!(adapter.graph().degree(adapter_node), 2);
    }
}

#[rstest]
fn test_automorphism_adapter_automorphisms() {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    });
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
    let colors = initial_colors(
        &molecule,
        &incidence_graph,
        DescriptionLevel::Constitution,
        true,
    )
    .unwrap();
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);

    assert_eq!(
        adapter.automorphisms(AutomorphismAlgorithm::Nauty),
        AutomorphismAdapterOutput {
            source_orbits: vec![NodeId(0), NodeId(0), NodeId(2)],
            source_canonical_labels: vec![NodeId(0), NodeId(1), NodeId(2)],
            source_generators: vec![vec![NodeId(1), NodeId(0), NodeId(2)]],
        },
    );
}

#[rstest]
#[case::whole_molecule(
    None,
    vec![
        [1, 2, 3, 4, 5, 0, 8, 6, 9, 10, 11, 7].map(NodeId).into(),
        [0, 5, 4, 3, 2, 1, 7, 6, 11, 10, 9, 8].map(NodeId).into(),
    ],
    vec![
        [1, 2, 3, 4, 5, 0, 8, 6, 9, 10, 11, 7].map(NodeId).into(),
        [0, 5, 4, 3, 2, 1, 7, 6, 11, 10, 9, 8].map(NodeId).into(),
    ]
)]
#[case::one_atom(
    Some(vec![AtomId(0)]),
    vec![
        [1, 2, 3, 4, 5, 0, 8, 6, 9, 10, 11, 7].map(NodeId).into(),
        [0, 5, 4, 3, 2, 1, 7, 6, 11, 10, 9, 8].map(NodeId).into(),
    ],
    vec![[0, 5, 4, 3, 2, 1, 7, 6, 11, 10, 9, 8].map(NodeId).into()]
)]
#[case::two_opposite_atoms(
    Some(vec![AtomId(0), AtomId(3)]),
    vec![
        [1, 2, 3, 4, 5, 0, 8, 6, 9, 10, 11, 7].map(NodeId).into(),
        [3, 4, 5, 0, 1, 2, 10, 9, 11, 7, 6, 8].map(NodeId).into(),
        [0, 5, 4, 3, 2, 1, 7, 6, 11, 10, 9, 8].map(NodeId).into(),
    ],
    vec![
        [3, 4, 5, 0, 1, 2, 10, 9, 11, 7, 6, 8].map(NodeId).into(),
        [0, 5, 4, 3, 2, 1, 7, 6, 11, 10, 9, 8].map(NodeId).into(),
    ]
)]
fn test_retain_preserving_generators(
    #[case] atoms: Option<Vec<AtomId>>,
    #[case] mut generators: Vec<Vec<NodeId>>,
    #[case] expected: Vec<Vec<NodeId>>,
) {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(5), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
            (AtomId(3), AtomId(4), BondForm::from_order(1)),
            (AtomId(4), AtomId(5), BondForm::from_order(1)),
        ],
        constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms,
            sum: NumForm::Lit(0),
        })
        .into(),
        ..Default::default()
    });
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);

    retain_preserving_generators(&molecule, &incidence_graph, &mut generators);

    assert_eq!(generators, expected);
}

#[rstest]
#[case::order_four(4)]
fn test_ordered_partition_refine_exhaustive_domain(
    canonicalize_context: CanonicalizeContext,
    #[case] atom_count: usize,
) {
    fn compare_refinement(
        partition: OrderedPartition,
        renumbered_partition: OrderedPartition,
        graph: &Graph,
        renumbered_graph: &Graph,
    ) {
        let reference = partition.clone().refine_by_all_cells(graph);
        let active = partition.clone().refine_by_active_cells(graph);
        let actual = partition.refine(graph);
        assert_eq!(active, reference);
        assert_eq!(actual, reference);
        let renumbered = renumbered_partition.refine(renumbered_graph);
        let node_count = graph.node_count();
        let transported = OrderedPartition {
            cells: renumbered
                .cells
                .iter()
                .map(|cell| {
                    let mut cell = cell
                        .iter()
                        .map(|node| NodeId((node_count - 1 - node.index()) as u32))
                        .collect::<Vec<_>>();
                    cell.sort_unstable();
                    cell
                })
                .collect(),
        };
        assert_eq!(actual, transported);

        let Some(cell_index) = actual.first_non_singleton_entity_cell(graph.node_count()) else {
            return;
        };
        for node in actual.cells[cell_index].clone() {
            let image = NodeId((node_count - 1 - node.index()) as u32);
            compare_refinement(
                actual.individualize(cell_index, node),
                renumbered.individualize(cell_index, image),
                graph,
                renumbered_graph,
            );
        }
    }

    let endpoint_pairs = (0..atom_count as u32)
        .flat_map(|first| ((first + 1)..atom_count as u32).map(move |second| [first, second]))
        .collect::<Vec<_>>();

    for complete_bond_form in [false, true] {
        for edge_mask in 0..(1_u64 << endpoint_pairs.len()) {
            let selected_pairs = endpoint_pairs
                .iter()
                .enumerate()
                .filter_map(|(position, &pair)| ((edge_mask >> position) & 1 == 1).then_some(pair))
                .collect::<Vec<_>>();
            let graph = Graph::new(atom_count, &selected_pairs);
            let renumbered_pairs = selected_pairs
                .iter()
                .map(|&[first, second]| {
                    [
                        atom_count as u32 - 1 - first,
                        atom_count as u32 - 1 - second,
                    ]
                })
                .collect::<Vec<_>>();
            let renumbered_graph = Graph::new(atom_count, &renumbered_pairs);
            compare_refinement(
                OrderedPartition::from_descriptors(&vec![0_u32; atom_count]),
                OrderedPartition::from_descriptors(&vec![0_u32; atom_count]),
                &graph,
                &renumbered_graph,
            );

            let molecule = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); atom_count],
                bonds: selected_pairs
                    .iter()
                    .map(|&[first, second]| {
                        let bond = BondForm::from_order(1);
                        (
                            AtomId(first),
                            AtomId(second),
                            if complete_bond_form {
                                bond.into_concrete()
                            } else {
                                bond
                            },
                        )
                    })
                    .collect(),
                ..Default::default()
            });
            let renumbered = molecule.clone().remap(&reverse_remapping(&molecule));
            let (expected, expected_remapping) =
                canonicalize_topology(molecule.clone(), &canonicalize_context).unwrap();
            let (actual, actual_remapping) =
                canonicalize_topology(renumbered.clone(), &canonicalize_context).unwrap();

            assert_eq!(actual, expected);
            assert!(molecule.framed_eq_under(&expected, &expected_remapping));
            assert!(renumbered.framed_eq_under(&actual, &actual_remapping));
        }
    }
}

#[rstest]
#[case::renumbered(
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [1 5 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[5 4 "1"] [4 3 "1"] [3 2 "1"] [2 1 "1"] [4 0 "1"]]}"#),
    true
)]
#[case::different_constitution(
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [1 5 "1"]]}"#),
    false
)]
#[case::different_size(
    mol_dsl!(r#"{:atoms ["C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"]]}"#),
    false
)]
fn test_same_cells(#[case] source: Molecule, #[case] target: Molecule, #[case] expected: bool) {
    let source_incidence = source.incidence_graph(IncidenceLevel::Constitution);
    let target_incidence = target.incidence_graph(IncidenceLevel::Constitution);
    let source_adapter = AutomorphismAdapter::new(
        &source_incidence,
        &initial_colors(
            &source,
            &source_incidence,
            DescriptionLevel::Constitution,
            true,
        )
        .unwrap(),
    );
    let target_adapter = AutomorphismAdapter::new(
        &target_incidence,
        &initial_colors(
            &target,
            &target_incidence,
            DescriptionLevel::Constitution,
            true,
        )
        .unwrap(),
    );
    let source_partition =
        OrderedPartition::from_descriptors(&source_adapter.colors).refine(source_adapter.graph());
    let target_partition =
        OrderedPartition::from_descriptors(&target_adapter.colors).refine(target_adapter.graph());

    assert_eq!(
        same_cells(
            &source_adapter,
            &source_partition,
            &target_adapter,
            &target_partition
        ),
        expected,
    );
}

#[rstest]
#[case::renumbered(
    mol_dsl!(r#"{:atoms ["N" "B"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
    mol_dsl!(r#"{:atoms ["B" "N"] :dative-bonds [{:donors [1] :acceptor 0 :attrs "1"}]}"#),
    vec![1, 0, 2],
    true
)]
#[case::donor_and_acceptor_exchanged(
    mol_dsl!(r#"{:atoms ["N" "B"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
    mol_dsl!(r#"{:atoms ["N" "B"] :dative-bonds [{:donors [1] :acceptor 0 :attrs "1"}]}"#),
    vec![0, 1, 2],
    false
)]
#[case::electron_counts(
    mol_dsl!(
        r#"{:atoms ["O" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [0 4 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["O" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [0 4 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#
    ),
    vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    true
)]
#[case::electron_counts_moved(
    mol_dsl!(
        r#"{:atoms ["O" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [0 4 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["O" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [0 4 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,2,1,1,1]"}]}"#
    ),
    vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
    false
)]
#[case::bond_moved(
    mol_dsl!(r#"{:atoms ["C" "O" "N"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "O" "N"] :bonds [[0 2 "1"] [2 1 "1"]]}"#),
    vec![0, 1, 2, 3, 4],
    false
)]
fn test_remapping_preserves_incidences(
    #[case] source: Molecule,
    #[case] target: Molecule,
    #[case] images: Vec<u32>,
    #[case] expected: bool,
) {
    let source_incidence = source.incidence_graph(IncidenceLevel::Constitution);
    let target_incidence = target.incidence_graph(IncidenceLevel::Constitution);
    let source_adapter = AutomorphismAdapter::new(
        &source_incidence,
        &initial_colors(
            &source,
            &source_incidence,
            DescriptionLevel::Constitution,
            true,
        )
        .unwrap(),
    );
    let target_adapter = AutomorphismAdapter::new(
        &target_incidence,
        &initial_colors(
            &target,
            &target_incidence,
            DescriptionLevel::Constitution,
            true,
        )
        .unwrap(),
    );
    let images = images.into_iter().map(NodeId).collect::<Vec<_>>();

    assert_eq!(
        remapping_preserves_incidences(&source_adapter, &target_adapter, &images),
        expected,
    );
    assert_eq!(
        source.framed_eq_under(
            &target,
            &remapping_from_generator(&source, &source_incidence, &images),
        ),
        expected,
    );
}

#[rstest]
fn test_canonical_search_color_labels() {
    let molecule = Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::O),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::N),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
        ],
        ..Default::default()
    });
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
    let colors = initial_colors(
        &molecule,
        &incidence_graph,
        DescriptionLevel::Constitution,
        true,
    )
    .unwrap();
    let adapter = AutomorphismAdapter::new(&incidence_graph, &colors);
    let mut relabeled = adapter.clone();
    relabeled.colors.iter_mut().for_each(|color| {
        *color = match *color {
            AutomorphismColor::Entity(value) => AutomorphismColor::Entity(u32::MAX - value),
            AutomorphismColor::Incidence(value) => AutomorphismColor::Incidence(u32::MAX - value),
        }
    });
    let descriptors = adapter.colors.clone();
    let leaf_candidate = |order: &[NodeId]| {
        let entity_order = entity_order(
            &molecule,
            &incidence_graph,
            DescriptionLevel::Constitution,
            order,
        );
        CanonicalCandidate {
            key: color_key(&molecule, &incidence_graph, &colors, &entity_order),
            entity_order,
        }
    };

    let expected = canonical_search(
        &adapter,
        &descriptors,
        AutomorphismAlgorithm::Nauty,
        &leaf_candidate,
    );
    let actual = canonical_search(
        &relabeled,
        &descriptors,
        AutomorphismAlgorithm::Nauty,
        &leaf_candidate,
    );

    assert_eq!(actual.candidate.key, expected.candidate.key);
}

#[rstest]
#[case::cells_differ(
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [1 5 "1"]]}"#),
    false
)]
#[case::single_node_cells(
    mol_dsl!(r#"{:atoms ["C" "O" "N"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["N" "O" "C"] :bonds [[2 1 "1"] [1 0 "1"]]}"#),
    true
)]
#[case::single_node_cells_check_fails(
    mol_dsl!(r#"{:atoms ["C" "O" "N"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "O" "N"] :bonds [[0 2 "1"] [2 1 "1"]]}"#),
    false
)]
#[case::first_path(
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 2 "1"] [2 4 "1"] [4 1 "1"] [1 3 "1"] [3 5 "1"] [5 0 "1"]]}"#),
    true
)]
#[case::search_finds_remapping(
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"]
                    [6 7 "1"] [7 0 "1"] [0 4 "1"] [1 3 "1"] [2 6 "1"] [5 7 "1"]]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[7 6 "1"] [6 5 "1"] [5 4 "1"] [4 3 "1"] [3 2 "1"] [2 1 "1"]
                    [1 0 "1"] [0 7 "1"] [7 3 "1"] [6 4 "1"] [5 1 "1"] [2 0 "1"]]}"#
    ),
    true
)]
#[case::search_finds_none(
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"] [4 5 "1"] [5 6 "1"]
                    [6 7 "1"] [7 4 "1"] [0 4 "1"] [1 5 "1"] [2 6 "1"] [3 7 "1"]]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"]
                    [6 7 "1"] [7 0 "1"] [0 4 "1"] [1 3 "1"] [2 6 "1"] [5 7 "1"]]}"#
    ),
    false
)]
fn test_remapping_exists(
    #[case] source: Molecule,
    #[case] target: Molecule,
    #[case] expected: bool,
) {
    let source_incidence = source.incidence_graph(IncidenceLevel::Topology);
    let target_incidence = target.incidence_graph(IncidenceLevel::Topology);
    let source_colors =
        initial_colors(&source, &source_incidence, DescriptionLevel::Topology, true).unwrap();
    let target_colors =
        initial_colors(&target, &target_incidence, DescriptionLevel::Topology, true).unwrap();
    let TopologySearchInput {
        adapter: source_adapter,
        initial_partition: source_partition,
    } = topology_search_input(&source, &source_incidence, &source_colors.entities);
    let TopologySearchInput {
        adapter: target_adapter,
        initial_partition: target_partition,
    } = topology_search_input(&target, &target_incidence, &target_colors.entities);

    assert_eq!(
        remapping_exists(
            &source_adapter,
            source_partition,
            &target_adapter,
            target_partition,
            AutomorphismAlgorithm::Nauty,
            &|images: &[NodeId]| {
                remapping_preserves_incidences(&source_adapter, &target_adapter, images)
            },
            false,
            &|_: &mut Vec<Vec<NodeId>>| {},
        ),
        expected,
    );
}

#[rstest]
#[case::mixed_bond_orders(Molecule::from_entries(MoleculeEntries {
    atoms: vec![AtomForm::from_element(Element::C); 4],
    bonds: vec![
        (AtomId(0), AtomId(1), BondForm::from_order(1).into_concrete()),
        (AtomId(1), AtomId(2), BondForm::from_order(2).into_concrete()),
        (AtomId(2), AtomId(3), BondForm::from_order(1).into_concrete()),
        (AtomId(3), AtomId(0), BondForm::from_order(3).into_concrete()),
    ],
    ..Default::default()
}))]
#[case::electronic_bond_forms(Molecule::from_entries(MoleculeEntries {
    atoms: vec![AtomForm::from_element(Element::C); 4],
    bonds: vec![
        (AtomId(0), AtomId(1), BondForm::from_order(1).into_concrete()),
        (
            AtomId(1),
            AtomId(2),
            BondForm::from_order(1).with_charge(1_i64).into_concrete(),
        ),
        (
            AtomId(2),
            AtomId(3),
            BondForm::from_order(1)
                .with_charge(0_i64)
                .with_unpaired_electrons((2_u8, 3_u8)),
        ),
    ],
    ..Default::default()
}))]
#[case::disconnected(Molecule::from_entries(MoleculeEntries {
    atoms: vec![
        AtomForm::from_element(Element::C),
        AtomForm::from_element(Element::N),
        AtomForm::from_element(Element::C),
        AtomForm::from_element(Element::N),
        AtomForm::from_element(Element::O),
    ],
    bonds: vec![
        (AtomId(0), AtomId(1), BondForm::from_order(1).into_concrete()),
        (AtomId(2), AtomId(3), BondForm::from_order(2).into_concrete()),
    ],
    ..Default::default()
}))]
#[case::topology_path(Molecule::from_entries(MoleculeEntries {
    atoms: vec![AtomForm::from_element(Element::C); 32],
    bonds: (0_u32..31)
        .map(|first| (
            AtomId(first),
            AtomId(first + 1),
            BondForm::from_order(1),
        ))
        .collect(),
    ..Default::default()
}))]
fn test_canonicalize_topology_relabeling(
    canonicalize_context: CanonicalizeContext,
    #[case] molecule: Molecule,
) {
    let expected = molecule
        .clone()
        .tracked_canonicalize_at(DescriptionLevel::Topology, &canonicalize_context)
        .map(|(canonical, _)| canonical)
        .unwrap();
    let renumbered = molecule.clone().remap(&reverse_remapping(&molecule));

    for source in [&molecule, &renumbered, &expected] {
        let (actual, remapping) =
            canonicalize_topology(source.clone(), &canonicalize_context).unwrap();

        assert_eq!(actual, expected);
        assert!(source.framed_eq_under(&actual, &remapping));
    }
}

#[rstest]
#[case::order_four(4)]
fn test_canonicalize_topology_exhaustive_domain(
    canonicalize_context: CanonicalizeContext,
    #[case] atom_count: usize,
) {
    let endpoint_pairs = (0..atom_count as u32)
        .flat_map(|first| ((first + 1)..atom_count as u32).map(move |second| [first, second]))
        .collect::<Vec<_>>();

    for complete_bond_form in [false, true] {
        for edge_mask in 0..(1_u64 << endpoint_pairs.len()) {
            let bonds = endpoint_pairs
                .iter()
                .enumerate()
                .filter_map(|(position, &[first, second])| {
                    if (edge_mask >> position) & 1 == 0 {
                        return None;
                    }
                    let bond = BondForm::from_order(1);
                    Some((
                        AtomId(first),
                        AtomId(second),
                        if complete_bond_form {
                            bond.into_concrete()
                        } else {
                            bond
                        },
                    ))
                })
                .collect::<Vec<_>>();
            let molecule = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm::from_element(Element::C); atom_count],
                bonds,
                ..Default::default()
            });

            let (canonical, remapping) =
                canonicalize_topology(molecule.clone(), &canonicalize_context).unwrap();

            assert!(
                molecule.framed_eq_under(&canonical, &remapping),
                "complete bond form {complete_bond_form}, edge mask {edge_mask:#08b}",
            );
            assert_eq!(canonical.check_integrity(), Ok(()));
            assert_eq!(
                canonical
                    .clone()
                    .tracked_canonicalize_at(DescriptionLevel::Topology, &canonicalize_context)
                    .map(|(canonical, _)| canonical),
                Ok(canonical.clone()),
                "complete bond form {complete_bond_form}, edge mask {edge_mask:#08b}",
            );

            for (index, atom_images) in permutations(atom_count).into_iter().enumerate() {
                let bond_count = molecule.bonds().count();
                let bond_images = if index % 2 == 0 {
                    (0..bond_count).collect::<Vec<_>>()
                } else {
                    (0..bond_count).rev().collect::<Vec<_>>()
                };
                let renumbering = molecule_remapping(&[
                    atom_images,
                    bond_images,
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                ]);
                let renumbered = molecule.clone().remap(&renumbering);

                assert_eq!(
                    renumbered.clone().tracked_canonicalize_at(DescriptionLevel::Topology, &canonicalize_context).map(|(canonical, _)| canonical),
                    Ok(canonical.clone()),
                    "complete bond form {complete_bond_form}, edge mask {edge_mask:#08b}, renumbering {index}",
                );
            }
        }
    }
}

#[rstest]
fn test_canonicalize_constitution_entity_kind(canonicalize_context: CanonicalizeContext) {
    let atoms = vec![AtomForm::from_element(Element::C); 4];
    let cases = [
        (
            "dative",
            Molecule::from_entries(MoleculeEntries {
                atoms: atoms.clone(),
                dative: vec![
                    (
                        vec![AtomId(0)],
                        AtomId(2),
                        DativeBondForm::new(NumForm::RangeFrom(1)),
                    ),
                    (
                        vec![AtomId(1)],
                        AtomId(3),
                        DativeBondForm::new(NumForm::RangeFrom(1)),
                    ),
                ],
                ..Default::default()
            }),
        ),
        (
            "aromatic",
            Molecule::from_entries(MoleculeEntries {
                atoms: atoms.clone(),
                aromatic: vec![
                    (
                        vec![AtomId(0), AtomId(1)],
                        AromaticSystemForm::from_electrons(vec![1, 2])
                            .with_charge(NumForm::var("q")),
                    ),
                    (
                        vec![AtomId(2), AtomId(3)],
                        AromaticSystemForm::from_electrons(vec![1, 2])
                            .with_charge(NumForm::var("q")),
                    ),
                ],
                ..Default::default()
            }),
        ),
        (
            "multicenter",
            Molecule::from_entries(MoleculeEntries {
                atoms: atoms.clone(),
                multicenter: vec![
                    (
                        vec![AtomId(0), AtomId(1)],
                        MulticenterBondForm::from_electrons(vec![2, 1]),
                    ),
                    (
                        vec![AtomId(2), AtomId(3)],
                        MulticenterBondForm::from_electrons(vec![2, 1]),
                    ),
                ],
                ..Default::default()
            }),
        ),
        (
            "noncovalent",
            Molecule::from_entries(MoleculeEntries {
                atoms,
                noncovalent: vec![
                    (
                        [AtomId(0), AtomId(1)],
                        NoncovalentBondForm::new(NoncovalentBondKindForm::Undetermined),
                    ),
                    (
                        [AtomId(2), AtomId(3)],
                        NoncovalentBondForm::new(NoncovalentBondKindForm::Undetermined),
                    ),
                ],
                ..Default::default()
            }),
        ),
    ];

    for (entity_kind, molecule) in cases {
        let (canonical, remapping) =
            canonicalize_constitution(molecule.clone(), &canonicalize_context).unwrap();

        assert!(
            molecule.framed_eq_under(&canonical, &remapping),
            "{entity_kind}"
        );
        assert_eq!(canonical.check_integrity(), Ok(()), "{entity_kind}");

        for (index, atom_images) in permutations(molecule.atoms().count())
            .into_iter()
            .enumerate()
        {
            let mut images = molecule
                .entity_counts()
                .map(|count| (0..count).collect::<Vec<_>>());
            images[0] = atom_images;
            if index % 2 == 1 {
                for entity_kind_images in &mut images[1..6] {
                    entity_kind_images.reverse();
                }
            }
            let renumbered = molecule.clone().remap(&molecule_remapping(&images));

            assert_eq!(
                renumbered
                    .clone()
                    .tracked_canonicalize_at(DescriptionLevel::Constitution, &canonicalize_context)
                    .map(|(canonical, _)| canonical),
                Ok(canonical.clone()),
                "{entity_kind}, renumbering {index}",
            );
        }
    }
}

#[rstest]
fn test_canonicalize_constitution_participant_order(canonicalize_context: CanonicalizeContext) {
    let left = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        aromatic: vec![
            (
                vec![AtomId(0), AtomId(1), AtomId(2)],
                AromaticSystemForm::from_electrons(vec![1, 2, 3]),
            ),
            (
                vec![AtomId(3), AtomId(4), AtomId(5)],
                AromaticSystemForm::from_electrons(vec![3, 2, 1]),
            ),
        ],
        multicenter: vec![
            (
                vec![AtomId(0), AtomId(3), AtomId(5)],
                MulticenterBondForm::from_electrons(vec![1, 2, 3]),
            ),
            (
                vec![AtomId(1), AtomId(2), AtomId(4)],
                MulticenterBondForm::from_electrons(vec![3, 2, 1]),
            ),
        ],
        ..Default::default()
    });
    let right = Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        aromatic: vec![
            (
                vec![AtomId(2), AtomId(0), AtomId(1)],
                AromaticSystemForm::from_electrons(vec![3, 1, 2]),
            ),
            (
                vec![AtomId(4), AtomId(5), AtomId(3)],
                AromaticSystemForm::from_electrons(vec![2, 1, 3]),
            ),
        ],
        multicenter: vec![
            (
                vec![AtomId(5), AtomId(0), AtomId(3)],
                MulticenterBondForm::from_electrons(vec![3, 1, 2]),
            ),
            (
                vec![AtomId(2), AtomId(4), AtomId(1)],
                MulticenterBondForm::from_electrons(vec![2, 1, 3]),
            ),
        ],
        ..Default::default()
    });

    assert_eq!(
        right
            .clone()
            .tracked_canonicalize_at(DescriptionLevel::Constitution, &canonicalize_context)
            .map(|(canonical, _)| canonical),
        left.clone()
            .tracked_canonicalize_at(DescriptionLevel::Constitution, &canonicalize_context)
            .map(|(canonical, _)| canonical),
    );
}

#[rstest]
#[case::atom_charge(
    r#"{:atoms ["C#c+" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C#c+" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::atom_charge_absent(
    r#"{:atoms ["C#c+" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::atom_charge_pair(
    r#"{:atoms ["C#c+" "C#c+" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C#c+" "C" "C" "C#c+" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::atom_constraint(
    r#"{:atoms ["C#v4" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C#v4" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::atom_constraint_absent(
    r#"{:atoms ["C#v4" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::atom_constraint_pair(
    r#"{:atoms ["C#v4" "C#v4" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C#v4" "C" "C" "C#v4" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_order(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "2"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_order_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_order_pair(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "2"] [1 2 "2"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "2"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::bond_constraint(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1#a"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1#a"] [4 5 "1"] [0 5 "1"]]}"#,
    true
)]
#[case::bond_constraint_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1#a"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
#[case::dative_bond(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [0] :acceptor 3 :attrs "1"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [1] :acceptor 4 :attrs "1"}]}"#,
    true
)]
#[case::dative_bond_distance(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [0] :acceptor 3 :attrs "1"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [0] :acceptor 2 :attrs "1"}]}"#,
    false
)]
#[case::dative_bond_role(
    r#"{:atoms ["N" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [0] :acceptor 3 :attrs "1"}]}"#,
    r#"{:atoms ["N" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :dative-bonds [{:donors [3] :acceptor 0 :attrs "1"}]}"#,
    false
)]
#[case::aromatic_electron_counts(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[2,1,1,1,1,1]"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,2,1,1]"}]}"#,
    true
)]
#[case::aromatic_electron_counts_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[2,1,1,1,1,1]"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]"}]}"#,
    false
)]
#[case::aromatic_electron_counts_pair(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[2,2,1,1,1,1]"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[2,1,1,2,1,1]"}]}"#,
    false
)]
#[case::aromatic_system_charge(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*#c-"} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:atoms [3 4 5] :attrs "*#c-"}]}"#,
    true
)]
#[case::aromatic_system_charge_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*#c-"} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::aromatic_system_constraint(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*#e2"} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:atoms [3 4 5] :attrs "*#e2"}]}"#,
    true
)]
#[case::aromatic_system_constraint_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*#e2"} {:atoms [3 4 5] :attrs "*"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [0 2 "1"] [3 4 "1"] [4 5 "1"] [3 5 "1"]]
          :aromatic-systems [{:atoms [0 1 2] :attrs "*"} {:atoms [3 4 5] :attrs "*"}]}"#,
    false
)]
#[case::multicenter_electron_counts(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:atoms [0 1 2] :attrs "[2,0,0]"}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:atoms [0 1 2] :attrs "[0,0,2]"}]}"#,
    true
)]
#[case::multicenter_electron_counts_distinct(
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:atoms [0 1 2] :attrs "[2,0,0]"}]}"#,
    r#"{:atoms ["B" "B" "B"]
          :multicenter-bonds [{:atoms [0 1 2] :attrs "[1,1,0]"}]}"#,
    false
)]
#[case::noncovalent_bond(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:atoms [0 3] :attrs "Hbd"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:atoms [1 4] :attrs "Hbd"}]}"#,
    true
)]
#[case::noncovalent_bond_distance(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:atoms [0 3] :attrs "Hbd"}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :noncovalent-bonds [{:atoms [0 2] :attrs "Hbd"}]}"#,
    false
)]
#[case::stereo_configuration(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    true
)]
#[case::stereo_configuration_distinct(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th1"}]}"#,
    false
)]
#[case::stereo_stereogenicity(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#g/"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0#g/"}]}"#,
    true
)]
#[case::stereo_stereogenicity_absent(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#g/"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::stereo_topicity(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,1)'"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0#o(0,1)'"}]}"#,
    true
)]
#[case::stereo_topicity_pair(
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,1)'"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    r#"{:atoms ["C" "F" "Cl" "Br" "C" "F" "Cl" "Br"]
          :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"] [4 5 "1"] [4 6 "1"] [4 7 "1"]]
          :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th0#o(0,2)'"}
                         {:site 4 :ligands [5 6 7 0] :attrs "Th0"}]}"#,
    false
)]
#[case::molecule_constraint(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0] :sum 0}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [3] :sum 0}}]}"#,
    true
)]
#[case::molecule_constraint_opposite(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0 3] :sum 0}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [1 4] :sum 0}}]}"#,
    true
)]
#[case::molecule_constraint_adjacent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0 1] :sum 0}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [3 2] :sum 0}}]}"#,
    true
)]
#[case::molecule_constraint_count(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0] :sum 0}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0 3] :sum 0}}]}"#,
    false
)]
#[case::molecule_constraint_pair(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0 1] :sum 0}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0 3] :sum 0}}]}"#,
    false
)]
#[case::molecule_constraint_absent(
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]
          :constraints [{:charge-sum {:atoms [0] :sum 0}}]}"#,
    r#"{:atoms ["C" "C" "C" "C" "C" "C"]
          :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [0 5 "1"]]}"#,
    false
)]
fn test_molecule_canonicalize_placement(
    canonicalize_context: CanonicalizeContext,
    #[case] lhs: &str,
    #[case] rhs: &str,
    #[case] expected: bool,
) {
    let [lhs, rhs] =
        [lhs, rhs].map(|source| Molecule::from_edn_str(source).expect("case molecule parses"));

    assert_eq!(
        lhs.clone().canonicalize(&canonicalize_context)
            == rhs.clone().canonicalize(&canonicalize_context),
        expected
    );
    assert_eq!(lhs.canonical_eq(&rhs, &canonicalize_context), expected);
    if expected {
        assert_eq!(
            lhs.canonical_hash(&canonicalize_context),
            rhs.canonical_hash(&canonicalize_context),
        );
    }
}

#[rstest]
fn test_molecule_canonicalize_stereo_frame(
    stereo_constraint_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let constraint = StereoAtomConstraintForm::Topicity(TopicityForm {
        pair: StereoLigandPair::new(0usize.into(), 2usize.into()),
        relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
    });
    let mut editor = stereo_constraint_canonicalization_molecule.edit();
    editor
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .constraints = constraint.clone().into();
    *editor.constraints_mut() =
        Constraint::StereoAtom(StereoAtomId(0), StereoKind::Tetrahedral, constraint).into();
    let source = editor.finish().unwrap();
    let reframed = source
        .clone()
        .reframe()
        .expect("integrity-valid molecule reframes");
    let canonical = source
        .clone()
        .canonicalize(&canonicalize_context)
        .expect("fixed molecule canonicalizes");

    assert_eq!(
        reframed.canonicalize(&canonicalize_context),
        Ok(canonical.clone()),
    );
    assert_eq!(
        canonical.clone().canonicalize(&canonicalize_context),
        Ok(canonical),
    );
}

#[rstest]
#[case::atom_charge(Molecule::from_entries(MoleculeEntries {
    atoms: vec![AtomForm::from_element(Element::C).with_charge(NumForm::lit_set(Vec::<i64>::new()))],
    ..Default::default()
}))]
#[case::atom_constraint(Molecule::from_entries(MoleculeEntries {
    atoms: vec![AtomForm::from_element(Element::C)
        .with_constraint(AtomConstraintForm::Valence(NumForm::lit_set(Vec::<i64>::new())))],
    ..Default::default()
}))]
#[case::dative_order(Molecule::from_entries(MoleculeEntries {
    atoms: vec![
        AtomForm::from_element(Element::N),
        AtomForm::from_element(Element::B),
    ],
    dative: vec![(
        vec![AtomId(0)],
        AtomId(1),
        DativeBondForm::new(NumForm::lit_set(Vec::<i64>::new())),
    )],
    ..Default::default()
}))]
#[case::stereo_coset(Molecule::from_entries(MoleculeEntries {
    atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
        .map(AtomForm::from_element)
        .into(),
    bonds: (1..=4)
        .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
        .collect(),
    stereo_atoms: vec![(
        AtomId(0),
        (1..=4)
            .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
            .collect(),
        StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::lit_set([])),
    )],
    ..Default::default()
}))]
#[case::electron_count_length(Molecule::from_entries(MoleculeEntries {
    atoms: vec![AtomForm::from_element(Element::C); 3],
    aromatic: vec![(
        vec![AtomId(0), AtomId(1), AtomId(2)],
        AromaticSystemForm::from_electrons(vec![2, 1]),
    )],
    ..Default::default()
}))]
fn test_molecule_canonicalize_contradiction(
    canonicalize_context: CanonicalizeContext,
    #[case] molecule: Molecule,
) {
    assert_eq!(
        molecule.canonicalize(&canonicalize_context),
        Err(MoleculeCanonicalizeError::Contradiction(Contradiction)),
    );
}

#[rstest]
fn test_molecule_canonical_eq(
    initial_color_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
) {
    let renumbered = initial_color_molecule
        .clone()
        .remap(&reverse_remapping(&initial_color_molecule));

    assert!(initial_color_molecule.canonical_eq(&renumbered, &canonicalize_context));
}

#[rstest]
#[case::ethanol_ethylamine(
    mol_dsl!(r#"{:atoms ["C" "C" "O"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "C" "N"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
    false
)]
#[case::hexane_2_methylpentane(
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [1 5 "1"]]}"#),
    false
)]
#[case::bromochlorofluoroiodomethane_renumbered(
    mol_dsl!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["I" "Br" "Cl" "F" "C"]
            :bonds [[4 0 "1"] [4 1 "1"] [4 2 "1"] [4 3 "1"]]
            :stereo-atoms [{:site 4 :ligands [3 2 1 0] :attrs "Th1"}]}"#
    ),
    true
)]
#[case::bromochlorofluoroiodomethane_mirror_image(
    mol_dsl!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["I" "Br" "Cl" "F" "C"]
            :bonds [[4 0 "1"] [4 1 "1"] [4 2 "1"] [4 3 "1"]]
            :stereo-atoms [{:site 4 :ligands [3 2 1 0] :attrs "Th0"}]}"#
    ),
    false
)]
#[case::decalin_renumbered(
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]
                    [5 6 "1"] [6 7 "1"] [7 8 "1"] [8 9 "1"] [9 4 "1"]]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[9 8 "1"] [8 7 "1"] [7 6 "1"] [6 5 "1"] [5 4 "1"] [4 9 "1"]
                    [4 3 "1"] [3 2 "1"] [2 1 "1"] [1 0 "1"] [0 5 "1"]]}"#
    ),
    true
)]
#[case::cyclohexane_two_cyclopropanes(
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "C" "C" "C" "C" "C"] :bonds [[0 1 "1"] [1 2 "1"] [2 0 "1"] [3 4 "1"] [4 5 "1"] [5 3 "1"]]}"#),
    false
)]
#[case::cuneane_renumbered(
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"]
                    [6 7 "1"] [7 0 "1"] [0 4 "1"] [1 3 "1"] [2 6 "1"] [5 7 "1"]]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[7 6 "1"] [6 5 "1"] [5 4 "1"] [4 3 "1"] [3 2 "1"] [2 1 "1"]
                    [1 0 "1"] [0 7 "1"] [7 3 "1"] [6 4 "1"] [5 1 "1"] [2 0 "1"]]}"#
    ),
    true
)]
#[case::cubane_cuneane(
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 0 "1"] [4 5 "1"] [5 6 "1"]
                    [6 7 "1"] [7 4 "1"] [0 4 "1"] [1 5 "1"] [2 6 "1"] [3 7 "1"]]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 6 "1"]
                    [6 7 "1"] [7 0 "1"] [0 4 "1"] [1 3 "1"] [2 6 "1"] [5 7 "1"]]}"#
    ),
    false
)]
#[case::meso_dichlorobutane_mirror_image(
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "Cl" "Cl"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [1 4 "1"] [2 5 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 4 2 [:h 1]] :attrs "Th1"}
                           {:site 2 :ligands [3 5 1 [:h 2]] :attrs "Th0"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "Cl" "Cl"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [1 4 "1"] [2 5 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 4 2 [:h 1]] :attrs "Th0"}
                           {:site 2 :ligands [3 5 1 [:h 2]] :attrs "Th1"}]}"#
    ),
    true
)]
#[case::meso_and_chiral_dichlorobutane(
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "Cl" "Cl"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [1 4 "1"] [2 5 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 4 2 [:h 1]] :attrs "Th1"}
                           {:site 2 :ligands [3 5 1 [:h 2]] :attrs "Th0"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "Cl" "Cl"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [1 4 "1"] [2 5 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 4 2 [:h 1]] :attrs "Th1"}
                           {:site 2 :ligands [3 5 1 [:h 2]] :attrs "Th1"}]}"#
    ),
    false
)]
#[case::dichlorobutane_enantiomers(
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "Cl" "Cl"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [1 4 "1"] [2 5 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 4 2 [:h 1]] :attrs "Th1"}
                           {:site 2 :ligands [3 5 1 [:h 2]] :attrs "Th1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C" "Cl" "Cl"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [1 4 "1"] [2 5 "1"]]
            :stereo-atoms [{:site 1 :ligands [0 4 2 [:h 1]] :attrs "Th0"}
                           {:site 2 :ligands [3 5 1 [:h 2]] :attrs "Th0"}]}"#
    ),
    false
)]
#[case::charge_sum_rotated(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
            (AtomId(3), AtomId(4), BondForm::from_order(1)),
            (AtomId(4), AtomId(5), BondForm::from_order(1)),
            (AtomId(5), AtomId(0), BondForm::from_order(1)),
        ],
        constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0), AtomId(1)]),
            sum: NumForm::Lit(0),
        })
        .into(),
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
            (AtomId(3), AtomId(4), BondForm::from_order(1)),
            (AtomId(4), AtomId(5), BondForm::from_order(1)),
            (AtomId(5), AtomId(0), BondForm::from_order(1)),
        ],
        constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(3), AtomId(4)]),
            sum: NumForm::Lit(0),
        })
        .into(),
        ..Default::default()
    }),
    true
)]
#[case::charge_sum_other_atoms(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
            (AtomId(3), AtomId(4), BondForm::from_order(1)),
            (AtomId(4), AtomId(5), BondForm::from_order(1)),
            (AtomId(5), AtomId(0), BondForm::from_order(1)),
        ],
        constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0), AtomId(1)]),
            sum: NumForm::Lit(0),
        })
        .into(),
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
            (AtomId(2), AtomId(3), BondForm::from_order(1)),
            (AtomId(3), AtomId(4), BondForm::from_order(1)),
            (AtomId(4), AtomId(5), BondForm::from_order(1)),
            (AtomId(5), AtomId(0), BondForm::from_order(1)),
        ],
        constraints: Constraint::Molecule(MoleculeConstraint::ChargeSum {
            atoms: Some(vec![AtomId(0), AtomId(2)]),
            sum: NumForm::Lit(0),
        })
        .into(),
        ..Default::default()
    }),
    false
)]
fn test_molecule_canonical_eq_compounds(
    canonicalize_context: CanonicalizeContext,
    #[case] left: Molecule,
    #[case] right: Molecule,
    #[case] expected: bool,
) {
    assert_eq!(left.canonical_eq(&right, &canonicalize_context), expected);
    assert_eq!(right.canonical_eq(&left, &canonicalize_context), expected);
}

#[rstest]
#[case::same_frame([4, 2, 1, 3], 0, true)]
#[case::odd_frame([2, 4, 1, 3], 1, true)]
#[case::even_frame([2, 1, 4, 3], 0, true)]
#[case::opposite_odd_frame([2, 4, 1, 3], 0, false)]
#[case::opposite_even_frame([2, 1, 4, 3], 1, false)]
fn test_molecule_canonical_eq_stereo_atom(
    stereo_atom_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
    #[case] ligand_ids: [u32; 4],
    #[case] coset: u32,
    #[case] expected: bool,
    #[values(false, true)] relabel: bool,
) {
    let mut elements = [Element::C, Element::F, Element::Cl, Element::Br, Element::I];
    if relabel {
        elements.swap(1, 2);
    }
    let other = Molecule::from_entries(MoleculeEntries {
        atoms: elements.into_iter().map(AtomForm::from_element).collect(),
        bonds: (1..=4)
            .map(|atom| (AtomId(0), AtomId(atom), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            ligand_ids
                .into_iter()
                .map(|atom| {
                    let atom = match (relabel, atom) {
                        (true, 1) => 2,
                        (true, 2) => 1,
                        _ => atom,
                    };
                    StereoLigand::new(AtomId(atom), StereoLigandKind::Atom)
                })
                .collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, coset),
        )],
        ..Default::default()
    });
    let source = stereo_atom_canonicalization_molecule;
    let canonical = source.clone().canonicalize(&canonicalize_context).unwrap();
    let other_canonical = other.clone().canonicalize(&canonicalize_context).unwrap();

    if expected {
        assert_eq!(canonical, other_canonical);
    } else {
        assert_ne!(canonical, other_canonical);
    }
    assert_eq!(source.canonical_eq(&other, &canonicalize_context), expected);
    assert_eq!(other.canonical_eq(&source, &canonicalize_context), expected);
}

#[rstest]
#[case::same_frame([3, 2, 5, 4], 1, true)]
#[case::one_end([2, 3, 5, 4], 0, true)]
#[case::both_ends([2, 3, 4, 5], 1, true)]
#[case::exchanged_ends([5, 4, 3, 2], 1, true)]
#[case::opposite([2, 3, 5, 4], 1, false)]
fn test_molecule_canonical_eq_stereo_bond(
    stereo_bond_canonicalization_molecule: Molecule,
    canonicalize_context: CanonicalizeContext,
    #[case] ligand_ids: [u32; 4],
    #[case] coset: u32,
    #[case] expected: bool,
    #[values(false, true)] relabel: bool,
) {
    let mut elements = [
        Element::C,
        Element::C,
        Element::F,
        Element::Cl,
        Element::Br,
        Element::I,
    ];
    if relabel {
        elements.swap(2, 3);
    }
    let other = Molecule::from_entries(MoleculeEntries {
        atoms: elements.into_iter().map(AtomForm::from_element).collect(),
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(1), AtomId(4), BondForm::from_order(1)),
            (AtomId(1), AtomId(5), BondForm::from_order(1)),
        ],
        stereo_bonds: vec![(
            BondId(0),
            ligand_ids
                .into_iter()
                .map(|atom| {
                    let atom = match (relabel, atom) {
                        (true, 2) => 3,
                        (true, 3) => 2,
                        _ => atom,
                    };
                    StereoLigand::new(AtomId(atom), StereoLigandKind::Atom)
                })
                .collect(),
            StereoBondForm::new(StereoKind::CisTrans, coset),
        )],
        ..Default::default()
    });
    let source = stereo_bond_canonicalization_molecule;
    let canonical = source.clone().canonicalize(&canonicalize_context).unwrap();
    let other_canonical = other.clone().canonicalize(&canonicalize_context).unwrap();

    if expected {
        assert_eq!(canonical, other_canonical);
    } else {
        assert_ne!(canonical, other_canonical);
    }
    assert_eq!(source.canonical_eq(&other, &canonicalize_context), expected);
    assert_eq!(other.canonical_eq(&source, &canonicalize_context), expected);
}

#[rstest]
#[case::renumbering(
    mol_dsl!(r#"{:atoms ["C" "N"] :bonds [[0 1 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["N" "C"] :bonds [[0 1 "1"]]}"#),
    true
)]
#[case::atom_form(
    mol_dsl!(r#"{:atoms ["C" "N"] :bonds [[0 1 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["F" "C"] :bonds [[0 1 "1"]]}"#),
    false
)]
#[case::ammonia_borane_renumbered(
    mol_dsl!(r#"{:atoms ["N" "B"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
    mol_dsl!(r#"{:atoms ["B" "N"] :dative-bonds [{:donors [1] :acceptor 0 :attrs "1"}]}"#),
    true
)]
#[case::methoxyamine_methylhydroxylamine(
    mol_dsl!(r#"{:atoms ["C" "O" "N"] :bonds [[0 1 "1"] [1 2 "1"]]}"#),
    mol_dsl!(r#"{:atoms ["C" "O" "N"] :bonds [[0 2 "1"] [2 1 "1"]]}"#),
    false
)]
#[case::bromochlorofluoroiodomethane_renumbered(
    mol_dsl!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["I" "Br" "Cl" "F" "C"]
            :bonds [[4 0 "1"] [4 1 "1"] [4 2 "1"] [4 3 "1"]]
            :stereo-atoms [{:site 4 :ligands [3 2 1 0] :attrs "Th1"}]}"#
    ),
    true
)]
#[case::bromochlorofluoroiodomethane_mirror_image(
    mol_dsl!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["I" "Br" "Cl" "F" "C"]
            :bonds [[4 0 "1"] [4 1 "1"] [4 2 "1"] [4 3 "1"]]
            :stereo-atoms [{:site 4 :ligands [3 2 1 0] :attrs "Th0"}]}"#
    ),
    false
)]
fn test_molecule_canonical_eq_remapping_witness(
    canonicalize_context: CanonicalizeContext,
    #[case] left: Molecule,
    #[case] right: Molecule,
    #[case] expected: bool,
) {
    assert!(left.clone().canonicalize(&canonicalize_context).is_ok());
    assert!(right.clone().canonicalize(&canonicalize_context).is_ok());
    assert_eq!(
        explicitly_dense_equivalent(left.entity_counts(), |remapping| {
            left.framed_eq_under(&right, remapping)
        }),
        expected
    );
    assert_eq!(left.canonical_eq(&right, &canonicalize_context), expected);
}

#[rstest]
#[case::hydrogen(StereoLigandKind::ImplicitHydrogen, false)]
#[case::lone_pair(StereoLigandKind::LonePair, false)]
#[case::virtual_ligand(StereoLigandKind::Virtual, true)]
fn test_molecule_canonical_eq_virtual_ligand(
    canonicalize_context: CanonicalizeContext,
    #[case] kind: StereoLigandKind,
    #[case] expected: bool,
) {
    let source = Molecule::from_entries(MoleculeEntries {
        atoms: [Element::C, Element::F, Element::Cl, Element::Br]
            .into_iter()
            .map(AtomForm::from_element)
            .collect(),
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::Virtual),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
        )],
        ..Default::default()
    });
    let other = Molecule::from_entries(MoleculeEntries {
        atoms: [Element::Br, Element::C, Element::F, Element::Cl]
            .into_iter()
            .map(AtomForm::from_element)
            .collect(),
        bonds: vec![
            (AtomId(1), AtomId(0), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(1)),
            (AtomId(1), AtomId(3), BondForm::from_order(1)),
        ],
        stereo_atoms: vec![(
            AtomId(1),
            vec![
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), kind),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
        )],
        ..Default::default()
    });

    assert_eq!(source.canonical_eq(&other, &canonicalize_context), expected);
    assert_eq!(other.canonical_eq(&source, &canonicalize_context), expected);
}

#[rstest]
#[case::one_fails(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([])),
        ],
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
        ],
        ..Default::default()
    }),
    false
)]
#[case::both_fail_identical(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([])),
        ],
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([])),
        ],
        ..Default::default()
    }),
    true
)]
#[case::both_fail_equal_normal_forms(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([1])),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![2, 1]),
        )],
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(1_i64),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![2, 1]),
        )],
        ..Default::default()
    }),
    true
)]
#[case::both_fail_different_frames(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(1_i64),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![2, 1]),
        )],
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(1_i64),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ],
        aromatic: vec![(
            vec![AtomId(1), AtomId(0), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![2, 1]),
        )],
        ..Default::default()
    }),
    false
)]
#[case::both_fail_different(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([])),
        ],
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::N).with_charge(NumForm::lit_set([])),
        ],
        ..Default::default()
    }),
    false
)]
#[case::both_fail_different_counts(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(NumForm::lit_set([])),
        ],
        ..Default::default()
    }),
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::N).with_charge(NumForm::lit_set([])),
            AtomForm::from_element(Element::O),
        ],
        ..Default::default()
    }),
    false
)]
fn test_molecule_canonical_eq_contradiction(
    canonicalize_context: CanonicalizeContext,
    #[case] left: Molecule,
    #[case] right: Molecule,
    #[case] expected: bool,
) {
    assert_eq!(left.canonical_eq(&right, &canonicalize_context), expected);
}
