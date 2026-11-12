//! Unit and bounded-exhaustive tests for molecule canonicalization.

use rstest::{fixture, rstest};
use umol_chem::element::Element;
use umol_perm::Permutation;

use super::super::tests::{
    canonicalize_context, initial_color_molecule, para_stereo_canonicalization_molecule,
    permutations, reverse_remapping, stereo_atom_canonicalization_molecule,
    stereo_bond_canonicalization_molecule, stereo_constraint_canonicalization_molecule,
};
use super::super::*;
use super::*;
use crate::ir::{
    AromaticSystemForm, AromaticSystemId, AtomConstraintForm, AtomForm, AtomId, BondForm, BondId,
    Constraint, DativeBondForm, DativeBondId, Entity, IncidenceLevel, MoleculeConstraint,
    MoleculeCorrespondence, MoleculeEntries, MulticenterBondForm, MulticenterBondId,
    NoncovalentBondForm, NoncovalentBondId, NoncovalentBondKind, NoncovalentBondKindForm, NumForm,
    StereoAtomConstraintForm, StereoAtomForm, StereoAtomId, StereoBondForm, StereoBondId,
    StereoConfigurationForm, StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
    StereoLigandPair, Topicity, TopicityForm, TopicityRelationForm,
};

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
    let colors = initial_colors(&molecule, &incidence_graph, DescriptionLevel::Structure)
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
    let colors = initial_colors(&molecule, &incidence_graph, DescriptionLevel::Structure)
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
fn test_initial_colors(
    initial_color_molecule: Molecule,
    #[case] lhs: Entity,
    #[case] rhs: Entity,
    #[case] expected: bool,
) {
    let molecule = initial_color_molecule
        .normalize()
        .expect("fixed molecule normalizes");
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(&molecule, &incidence_graph, DescriptionLevel::Structure)
        .expect("fixed molecule has colors");

    assert_eq!(
        colors.entities[incidence_graph.node_of(lhs).index()]
            == colors.entities[incidence_graph.node_of(rhs).index()],
        expected
    );
}

#[rstest]
fn test_initial_colors_incidence(initial_color_molecule: Molecule) {
    let molecule = initial_color_molecule
        .normalize()
        .expect("fixed molecule normalizes");
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let colors = initial_colors(&molecule, &incidence_graph, DescriptionLevel::Structure)
        .expect("fixed molecule has colors");
    let contributions = [1, 2, 2, 1, 1, 2, 1, 2, 2, 1];
    let incidences = incidence_graph
        .incidences()
        .map(|(edge, incidence)| {
            let contribution = edge
                .index()
                .checked_sub(14)
                .and_then(|position| contributions.get(position));
            ((incidence, contribution), colors.incidences[edge.index()])
        })
        .collect::<Vec<_>>();

    for (lhs, lhs_color) in &incidences {
        for (rhs, rhs_color) in &incidences {
            assert_eq!(lhs_color == rhs_color, lhs == rhs);
        }
    }
    for entity_color in &colors.entities {
        for (_, incidence_color) in &incidences {
            assert_ne!(entity_color, incidence_color);
        }
    }
}

#[rstest]
#[case::aromatic_electron_count_length(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 3],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            AromaticSystemForm::from_electrons(vec![2, 1]),
        )],
        ..Default::default()
    }),
    DescriptionLevel::Constitution
)]
#[case::multicenter_electron_count_length(
    Molecule::from_entries(MoleculeEntries {
        atoms: [Element::B, Element::H, Element::B].map(AtomForm::from_element).into(),
        multicenter: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2)],
            MulticenterBondForm::from_electrons(vec![1, 0, 1, 0]),
        )],
        ..Default::default()
    }),
    DescriptionLevel::Constitution
)]
#[case::stereo_atom_coset_range(
    Molecule::from_entries(MoleculeEntries {
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
            StereoAtomForm::new(StereoKind::Tetrahedral, 2u32),
        )],
        ..Default::default()
    }),
    DescriptionLevel::Structure
)]
fn test_initial_colors_error(#[case] molecule: Molecule, #[case] level: DescriptionLevel) {
    let incidence_graph = molecule.incidence_graph(match level {
        DescriptionLevel::Topology => IncidenceLevel::Topology,
        DescriptionLevel::Constitution => IncidenceLevel::Constitution,
        DescriptionLevel::Structure | DescriptionLevel::Full => IncidenceLevel::Full,
    });

    assert_eq!(
        initial_colors(&molecule, &incidence_graph, level),
        Err(Contradiction)
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
    let colors = initial_colors(&molecule, &incidence_graph, DescriptionLevel::Topology).unwrap();
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
    let colors = initial_colors(&molecule, &incidence_graph, DescriptionLevel::Structure).unwrap();
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
    let colors =
        initial_colors(&molecule, &incidence_graph, DescriptionLevel::Constitution).unwrap();
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
#[case::feature_free(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
        ],
        bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(1))],
        ..Default::default()
    }),
    vec![NodeId(1), NodeId(0), NodeId(2)],
    true,
)]
#[case::distinct_ligands(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
        ],
        bonds: (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    }),
    (0..10).map(NodeId::from).collect(),
    true,
)]
#[case::prochiral_explicit_hydrogens(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::H),
            AtomForm::from_element(Element::H),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
        ],
        bonds: (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    }),
    vec![
        NodeId(0), NodeId(2), NodeId(1), NodeId(3), NodeId(4),
        NodeId(6), NodeId(5), NodeId(7), NodeId(8), NodeId(9),
    ],
    false,
)]
#[case::stereo_site_exchange(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
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
    }),
    vec![
        NodeId(5), NodeId(6), NodeId(7), NodeId(8), NodeId(9),
        NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4),
        NodeId(14), NodeId(15), NodeId(16), NodeId(17),
        NodeId(10), NodeId(11), NodeId(12), NodeId(13),
        NodeId(19), NodeId(18),
    ],
    true,
)]
#[case::stereo_bond_endpoint_exchange(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
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
            [2, 3, 4, 5]
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .into(),
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    }),
    vec![
        NodeId(1), NodeId(0), NodeId(4), NodeId(5), NodeId(2), NodeId(3),
        NodeId(6), NodeId(9), NodeId(10), NodeId(7), NodeId(8), NodeId(11),
    ],
    true,
)]
#[case::stereo_site_exchange_one_topicity(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
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
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))
                    .with_constraint(StereoAtomConstraintForm::Topicity(TopicityForm {
                        pair: StereoLigandPair::new(0usize.into(), 1usize.into()),
                        relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
                    })),
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
    }),
    vec![
        NodeId(5), NodeId(6), NodeId(7), NodeId(8), NodeId(9),
        NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4),
        NodeId(14), NodeId(15), NodeId(16), NodeId(17),
        NodeId(10), NodeId(11), NodeId(12), NodeId(13),
        NodeId(19), NodeId(18),
    ],
    false,
)]
#[case::stereo_site_exchange_both_topicity(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C),
            AtomForm::from_element(Element::F),
            AtomForm::from_element(Element::Cl),
            AtomForm::from_element(Element::Br),
            AtomForm::from_element(Element::I),
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
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))
                    .with_constraint(StereoAtomConstraintForm::Topicity(TopicityForm {
                        pair: StereoLigandPair::new(0usize.into(), 1usize.into()),
                        relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
                    })),
            ),
            (
                AtomId(5),
                (6..=9)
                    .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                    .collect(),
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0))
                    .with_constraint(StereoAtomConstraintForm::Topicity(TopicityForm {
                        pair: StereoLigandPair::new(0usize.into(), 1usize.into()),
                        relation: TopicityRelationForm::Lit(Topicity::Enantiotopic),
                    })),
            ),
        ],
        ..Default::default()
    }),
    vec![
        NodeId(5), NodeId(6), NodeId(7), NodeId(8), NodeId(9),
        NodeId(0), NodeId(1), NodeId(2), NodeId(3), NodeId(4),
        NodeId(14), NodeId(15), NodeId(16), NodeId(17),
        NodeId(10), NodeId(11), NodeId(12), NodeId(13),
        NodeId(19), NodeId(18),
    ],
    true,
)]
fn test_retain_stereo_preserving_generators(
    #[case] molecule: Molecule,
    #[case] generator: Vec<NodeId>,
    #[case] expected_retained: bool,
) {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Full);
    let mut generators = vec![generator.clone()];
    let expected = if expected_retained {
        vec![generator]
    } else {
        Vec::new()
    };

    retain_stereo_preserving_generators(&molecule, &incidence_graph, &mut generators);

    assert_eq!(generators, expected);
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
    let colors =
        initial_colors(&molecule, &incidence_graph, DescriptionLevel::Constitution).unwrap();
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
