//! Unit tests for the symmetry, the colors and the remapping checks of a molecule.

use pretty_assertions::assert_eq;
use rstest::*;
use umol_chem::element::Element;

use super::*;
use crate::ir::atom::AtomForm;
use crate::ir::bond::BondForm;
use crate::ir::coloring::ConstitutionColoring;
use crate::ir::id::{AtomId, BondId, StereoAtomId, StereoBondId, StereoLigandPosition};
use crate::ir::stereo::{StereoAtomForm, StereoBondForm, StereoConfigurationForm, StereoKind};
use crate::ir::{
    ArithExpr, AromaticSystemForm, AromaticSystemId, AtomConstraintForm, DativeBondForm,
    DativeBondId, MulticenterBondForm, MulticenterBondId, NoncovalentBondForm, NoncovalentBondId,
    NoncovalentBondKind, NumForm, StereoAtomConstraintForm, StereoLigandPair, TopicityForm,
    TopicityRelationForm,
};
use crate::mol_dsl;

fn config() -> GraphSymmetryConfig<ConstitutionColoring> {
    GraphSymmetryConfig {
        coloring: ConstitutionColoring::full(),
        iterate_to_fixpoint: true,
        max_iterations: 16,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    }
}

#[fixture]
fn benzene_ring() -> Molecule {
    let atoms = vec![AtomForm::from_element(Element::C); 6];
    let bonds = (0..6)
        .map(|i| (AtomId(i), AtomId((i + 1) % 6), BondForm::from_order(1)))
        .collect();
    Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds,
        ..Default::default()
    })
}

// A tetrahedral center on atom 0 with the four given peripheral elements.
fn tetrahedral(peripherals: [Element; 4]) -> Molecule {
    let mut atoms = vec![AtomForm::from_element(Element::C)];
    atoms.extend(peripherals.into_iter().map(AtomForm::from_element));
    let bonds = (1..=4)
        .map(|i| (AtomId(0), AtomId(i), BondForm::from_order(1)))
        .collect();
    Molecule::from_entries(MoleculeEntries {
        atoms,
        bonds,
        stereo_atoms: vec![(
            AtomId(0),
            vec![
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ],
            StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    })
}

#[rstest]
fn test_molecule_graph_symmetry_constitutional(benzene_ring: Molecule) {
    let symmetry = benzene_ring.graph_symmetry(&config());
    assert_eq!(
        symmetry.automorphism_algorithm,
        AutomorphismAlgorithm::Nauty
    );
    assert!(symmetry.same_proper_orbit(AtomId(0), AtomId(3)));
    assert_eq!(
        symmetry.proper_orbit_of(AtomId(0)),
        vec![
            AtomId(0),
            AtomId(1),
            AtomId(2),
            AtomId(3),
            AtomId(4),
            AtomId(5),
        ],
    );
    assert!(!symmetry.is_chiral());
}

#[rstest]
fn test_molecule_graph_symmetry_chiral() {
    // Four distinct ligands ⇒ a genuine stereocenter, no compensating symmetry.
    let mol = tetrahedral([Element::F, Element::Cl, Element::Br, Element::I]);
    let symmetry = mol.graph_symmetry(&config());
    assert!(symmetry.is_chiral());
}

#[rstest]
fn test_molecule_graph_symmetry_enantiotopic() {
    // Two identical Cl ligands: enantiotopic — same star orbit, distinct proper.
    let mol = tetrahedral([Element::Cl, Element::Cl, Element::F, Element::Br]);
    let symmetry = mol.graph_symmetry(&config());
    assert!(!symmetry.is_chiral());
    assert!(!symmetry.same_proper_orbit(AtomId(1), AtomId(2)));
    assert!(symmetry.same_star_orbit(AtomId(1), AtomId(2)));
}

#[rstest]
#[case::atom(
    MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 5],
        bonds: (1..=4)
            .map(|id| (AtomId(0), AtomId(id), BondForm::from_order(1)))
            .collect(),
        stereo_atoms: vec![(
            AtomId(0),
            (1..=4)
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            StereoAtomForm::default(),
        )],
        ..Default::default()
    },
    AtomId(1),
    vec![AtomId(1), AtomId(2), AtomId(3), AtomId(4)],
)]
#[case::bond(
    MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(2)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(1), AtomId(4), BondForm::from_order(1)),
            (AtomId(1), AtomId(5), BondForm::from_order(1)),
        ],
        stereo_bonds: vec![(
            BondId(0),
            (2..=5)
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom))
                .collect(),
            StereoBondForm::default(),
        )],
        ..Default::default()
    },
    AtomId(2),
    vec![AtomId(2), AtomId(3), AtomId(4), AtomId(5)],
)]
fn test_molecule_graph_symmetry_kindless(
    #[case] entries: MoleculeEntries,
    #[case] atom: AtomId,
    #[case] expected_orbit: Vec<AtomId>,
) {
    let molecule = Molecule::try_from_entries(entries).expect("kindless stereo is admitted");
    for coloring in [
        ConstitutionColoring::entity_only(),
        ConstitutionColoring::full(),
    ] {
        let symmetry = molecule.graph_symmetry(&GraphSymmetryConfig {
            coloring,
            iterate_to_fixpoint: false,
            max_iterations: 1,
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        });

        assert_eq!(symmetry.proper_orbit_of(atom), expected_orbit);
        assert_eq!(symmetry.is_chiral(), false);
    }
}

#[rstest]
fn test_graph_symmetry_site_stabilizer(benzene_ring: Molecule) {
    let symmetry = benzene_ring.graph_symmetry(&config());
    // Every stabilizer generator fixes the distinguished site.
    let generators = symmetry.site_stabilizer(NodeId(0));
    assert!(generators.iter().all(|g| g[0] == NodeId(0)));
}

#[rstest]
fn test_graph_symmetry_is_chiral_no_stereo(benzene_ring: Molecule) {
    // No stereocenters ⇒ trivially achiral despite no improper generators.
    assert!(!benzene_ring.graph_symmetry(&config()).is_chiral());
}

#[rstest]
fn test_molecule_stereo_atom_symmetry_stereogenic() {
    let mol = tetrahedral([Element::F, Element::Cl, Element::Br, Element::I]);
    let gs = mol.graph_symmetry(&config());
    let stereo = mol.stereo_atom_symmetry(&gs, StereoAtomId(0));
    assert!(stereo.is_stereogenic());
    // Four distinct ligands ⇒ no local symmetry ⇒ pairwise diastereotopic.
    assert_eq!(
        stereo.topicity(StereoLigandPosition(0), StereoLigandPosition(1)),
        Topicity::Diastereotopic
    );
}

#[rstest]
fn test_molecule_graph_symmetry_out_of_range_coset() {
    // A malformed (out-of-range) coset must not panic graph_symmetry: the
    // element contributes no observable descriptor and the orbit partition
    // is still computed.
    let mut mol = tetrahedral([Element::F, Element::Cl, Element::Br, Element::I]);
    mol.stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 9);
    let gs = mol.graph_symmetry(&config());
    assert_eq!(gs.proper_orbit_of(AtomId(0)), vec![AtomId(0)]);
}

#[rstest]
fn test_molecule_observable_descriptor() {
    let mut mol = tetrahedral([Element::Cl, Element::Cl, Element::F, Element::Br]);
    mol.stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1);
    assert_eq!(
        mol.observable_descriptor(
            Entity::StereoAtom(StereoAtomId(0)),
            &[NodeId(0), NodeId(1), NodeId(1), NodeId(3), NodeId(4)],
        ),
        Some(0),
    );
}

#[rstest]
fn test_molecule_stereo_atom_symmetry_prochiral() {
    let mol = tetrahedral([Element::Cl, Element::Cl, Element::F, Element::Br]);
    let gs = mol.graph_symmetry(&config());
    let stereo = mol.stereo_atom_symmetry(&gs, StereoAtomId(0));
    // Two identical Cl ⇒ prochiral, not a genuine stereocenter.
    assert!(!stereo.is_stereogenic());
    assert_eq!(
        stereo.topicity(StereoLigandPosition(0), StereoLigandPosition(1)),
        Topicity::Enantiotopic
    ); // the two Cl
    assert_eq!(
        stereo.topicity(StereoLigandPosition(0), StereoLigandPosition(2)),
        Topicity::Diastereotopic
    ); // Cl vs F
}

#[rstest]
fn test_molecule_stereo_bond_symmetry() {
    // C0=C1 with four distinct substituents (F,Cl on C0; Br,I on C1): E/Z stereogenic.
    let mol = Molecule::from_entries(MoleculeEntries {
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
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
            ],
            StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
        )],
        ..Default::default()
    });
    let gs = mol.graph_symmetry(&config());
    let stereo = mol.stereo_bond_symmetry(&gs, StereoBondId(0));
    assert!(stereo.is_stereogenic());
    assert_eq!(
        stereo.topicity(StereoLigandPosition(0), StereoLigandPosition(1)),
        Topicity::Diastereotopic
    ); // F vs Cl on C0
}

#[rstest]
#[case::singleton(
    OrientedPermutationGroup::generate(4, &[]),
    StereoKind::Tetrahedral,
    StereoCoset::Lit(0),
    true,
)]
#[case::merged(
    OrientedPermutationGroup::generate(
        4,
        &[OrientedPermutation::proper(Permutation::from_image(&[1, 0, 2, 3]))],
    ),
    StereoKind::Tetrahedral,
    StereoCoset::Lit(0),
    false,
)]
#[case::invalid_generator(
    OrientedPermutationGroup::generate(3, &[]),
    StereoKind::Tetrahedral,
    StereoCoset::Lit(0),
    false,
)]
#[case::invalid_coset(
    OrientedPermutationGroup::generate(4, &[]),
    StereoKind::Tetrahedral,
    StereoCoset::Lit(2),
    false,
)]
#[case::undetermined(
    OrientedPermutationGroup::generate(4, &[]),
    StereoKind::Tetrahedral,
    StereoCoset::Undetermined,
    true,
)]
#[case::undetermined_merged(
    OrientedPermutationGroup::generate(4, &[
        OrientedPermutation::improper(Permutation::from_image(&[1, 0, 2, 3])),
    ]),
    StereoKind::Tetrahedral,
    StereoCoset::Undetermined,
    false,
)]
#[case::cis_trans_undetermined(
    OrientedPermutationGroup::generate(4, &[]),
    StereoKind::CisTrans,
    StereoCoset::Undetermined,
    true,
)]
#[case::cis_trans_merged(
    OrientedPermutationGroup::generate(4, &[
        OrientedPermutation::proper(Permutation::from_image(&[1, 0, 2, 3])),
    ]),
    StereoKind::CisTrans,
    StereoCoset::Undetermined,
    false,
)]
#[case::cis_trans_invalid_action(
    OrientedPermutationGroup::generate(4, &[
        OrientedPermutation::proper(Permutation::from_image(&[0, 2, 1, 3])),
    ]),
    StereoKind::CisTrans,
    StereoCoset::Undetermined,
    false,
)]
#[case::finite_domain(
    OrientedPermutationGroup::generate(4, &[]),
    StereoKind::Tetrahedral,
    StereoCoset::lit_set([0, 1]),
    true,
)]
#[case::finite_singleton(
    OrientedPermutationGroup::generate(4, &[]),
    StereoKind::CisTrans,
    StereoCoset::lit_set([1]),
    true,
)]
#[case::empty_domain(
    OrientedPermutationGroup::generate(4, &[]),
    StereoKind::Tetrahedral,
    StereoCoset::lit_set([]),
    false,
)]
#[case::invalid_domain(
    OrientedPermutationGroup::generate(4, &[]),
    StereoKind::Tetrahedral,
    StereoCoset::lit_set([0, 2]),
    false,
)]
#[case::nonbinary_literal(
    OrientedPermutationGroup::generate(4, &[]),
    StereoKind::SquarePlanar,
    StereoCoset::Lit(0),
    true,
)]
#[case::nonbinary_undetermined(
    OrientedPermutationGroup::generate(4, &[]),
    StereoKind::SquarePlanar,
    StereoCoset::Undetermined,
    false,
)]
fn test_stereo_symmetry_is_stereogenic(
    #[case] group: OrientedPermutationGroup,
    #[case] kind: StereoKind,
    #[case] coset: StereoCoset,
    #[case] expected: bool,
) {
    assert_eq!(
        StereoSymmetry { group, kind, coset }.is_stereogenic(),
        expected,
    );
}

#[rstest]
#[case::valid(
    &[NodeId(1), NodeId(0)],
    Some(Permutation::from_image(&[1, 0])),
)]
#[case::collision(&[NodeId(0), NodeId(0)], None)]
fn test_project_onto_ligands(#[case] generator: &[NodeId], #[case] expected: Option<Permutation>) {
    assert_eq!(
        project_onto_ligands(
            &[
                StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
            ],
            generator,
        ),
        expected,
    );
}

#[rstest]
#[case::hydrogen(StereoLigandKind::ImplicitHydrogen)]
#[case::lone_pair(StereoLigandKind::LonePair)]
#[case::virtual_ligand(StereoLigandKind::Virtual)]
fn test_project_onto_ligands_virtual(#[case] kind: StereoLigandKind) {
    let ligands = [
        StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
        StereoLigand::new(AtomId(0), kind),
    ];
    assert_eq!(
        project_onto_ligands(&ligands, &[NodeId(0), NodeId(2), NodeId(1), NodeId(3)]),
        Some(Permutation::from_image(&[1, 0, 2, 3])),
    );
}

#[rstest]
#[case::rotation(
    mol_dsl!(
        r#"{:atoms ["C" "F" "F" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "F" "F" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    vec![0, 2, 3, 1, 4, 6, 7, 5, 8, 9],
    true
)]
#[case::transposition(
    mol_dsl!(
        r#"{:atoms ["C" "F" "F" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "F" "F" "F" "Cl"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    vec![0, 2, 1, 3, 4, 6, 5, 7, 8, 9],
    false
)]
#[case::site_moved(
    mol_dsl!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#
    ),
    vec![1, 0, 2, 3, 4, 5, 6, 7, 8, 9],
    false
)]
#[case::renumbered(
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
    vec![4, 3, 2, 1, 0, 8, 7, 6, 5, 9],
    true
)]
#[case::renumbered_inverted(
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
    vec![4, 3, 2, 1, 0, 8, 7, 6, 5, 9],
    false
)]
#[case::stereo_bond(
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#
    ),
    vec![0, 1, 2, 3, 4, 5, 6, 7],
    true
)]
#[case::stereo_bond_inverted(
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]}"#
    ),
    mol_dsl!(
        r#"{:atoms ["C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "2"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct0"}]}"#
    ),
    vec![0, 1, 2, 3, 4, 5, 6, 7],
    false
)]
fn test_remapping_preserves_stereo(
    #[case] source: Molecule,
    #[case] target: Molecule,
    #[case] images: Vec<u32>,
    #[case] expected: bool,
) {
    let incidence_graph = source.incidence_graph(IncidenceLevel::Full);
    let images = images.into_iter().map(NodeId).collect::<Vec<_>>();

    assert_eq!(
        remapping_preserves_stereo(&source, &target, &incidence_graph, &images),
        expected,
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
#[case::bond_endpoint(
    mol_dsl!(r#"{:atoms ["C" "C"] :bonds [[0 1 "1"]]}"#),
    Entity::Bond(BondId(0)),
    AtomId(0),
    (Incidence::BondEndpoint, None)
)]
#[case::dative_donor(
    mol_dsl!(r#"{:atoms ["N" "B"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
    Entity::DativeBond(DativeBondId(0)),
    AtomId(0),
    (Incidence::DativeDonor, None)
)]
#[case::dative_acceptor(
    mol_dsl!(r#"{:atoms ["N" "B"] :dative-bonds [{:donors [0] :acceptor 1 :attrs "1"}]}"#),
    Entity::DativeBond(DativeBondId(0)),
    AtomId(1),
    (Incidence::DativeAcceptor, None)
)]
#[case::aromatic_atom_oxygen(
    mol_dsl!(
        r#"{:atoms ["O" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [0 4 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#
    ),
    Entity::AromaticSystem(AromaticSystemId(0)),
    AtomId(0),
    (Incidence::AromaticAtom, Some(2))
)]
#[case::aromatic_atom_carbon(
    mol_dsl!(
        r#"{:atoms ["O" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [0 4 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[2,1,1,1,1]"}]}"#
    ),
    Entity::AromaticSystem(AromaticSystemId(0)),
    AtomId(1),
    (Incidence::AromaticAtom, Some(1))
)]
#[case::aromatic_atom_without_counts(
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![AtomForm::from_element(Element::C); 6],
        aromatic: vec![(
            vec![AtomId(0), AtomId(1), AtomId(2), AtomId(3), AtomId(4), AtomId(5)],
            AromaticSystemForm::default(),
        )],
        ..Default::default()
    }),
    Entity::AromaticSystem(AromaticSystemId(0)),
    AtomId(0),
    (Incidence::AromaticAtom, None)
)]
#[case::multicenter_atom(
    mol_dsl!(
        r#"{:atoms ["B" "B" "H"] :multicenter-bonds [{:atoms [0 2 1] :attrs "[1,1,0]"}]}"#
    ),
    Entity::MulticenterBond(MulticenterBondId(0)),
    AtomId(1),
    (Incidence::MulticenterAtom, Some(0))
)]
fn test_incidence_color(
    #[case] molecule: Molecule,
    #[case] owner: Entity,
    #[case] atom: AtomId,
    #[case] expected: (Incidence, Option<i64>),
) {
    let incidence_graph = molecule.incidence_graph(IncidenceLevel::Constitution);
    let edge = incidence_graph
        .graph()
        .find_edge(
            incidence_graph.node_of(owner),
            incidence_graph.node_of(Entity::Atom(atom)),
        )
        .unwrap();
    let (incidence, electron_count) = incidence_color(&molecule, &incidence_graph, edge);

    assert_eq!((incidence.clone(), electron_count), expected);
}

#[fixture]
pub(crate) fn initial_color_molecule() -> Molecule {
    let normalized_three = NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
        ArithExpr::Lit(1),
        ArithExpr::Lit(2),
    ])));
    let normalized_one = NumForm::ArithExpr(Box::new(ArithExpr::Sum(vec![
        ArithExpr::Lit(0),
        ArithExpr::Lit(1),
    ])));
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            AtomForm::from_element(Element::C).with_charge(3_i64),
            AtomForm::from_element(Element::C)
                .with_charge(normalized_three)
                .with_constraint(AtomConstraintForm::Valence(NumForm::Lit(4))),
            AtomForm::from_element(Element::O).with_charge(3_i64),
            AtomForm::from_element(Element::C).with_charge(4_i64),
            AtomForm::from_element(Element::N),
            AtomForm::from_element(Element::H),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::new(normalized_one.clone())),
            (AtomId(2), AtomId(3), BondForm::from_order(2)),
            (AtomId(0), AtomId(5), BondForm::from_order(1)),
        ],
        dative: vec![
            (vec![AtomId(0)], AtomId(4), DativeBondForm::from_order(1)),
            (
                vec![AtomId(1)],
                AtomId(4),
                DativeBondForm::new(normalized_one),
            ),
            (vec![AtomId(2)], AtomId(4), DativeBondForm::from_order(2)),
        ],
        aromatic: vec![
            (
                vec![AtomId(0), AtomId(1)],
                AromaticSystemForm::from_electrons(vec![1, 2]),
            ),
            (
                vec![AtomId(2), AtomId(3)],
                AromaticSystemForm::from_electrons(vec![2, 1]),
            ),
            (
                vec![AtomId(4), AtomId(5)],
                AromaticSystemForm::from_electrons(vec![1, 2]).with_charge(1_i64),
            ),
        ],
        multicenter: vec![
            (
                vec![AtomId(0), AtomId(2)],
                MulticenterBondForm::from_electrons(vec![1, 2]),
            ),
            (
                vec![AtomId(1), AtomId(3)],
                MulticenterBondForm::from_electrons(vec![2, 1]),
            ),
        ],
        noncovalent: vec![
            (
                [AtomId(0), AtomId(5)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            ),
            (
                [AtomId(1), AtomId(5)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            ),
            (
                [AtomId(2), AtomId(5)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic),
            ),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(5), StereoLigandKind::Atom),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(0)),
            ),
            (
                AtomId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(1), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::Tetrahedral, StereoCoset::Lit(1)),
            ),
            (
                AtomId(2),
                vec![
                    StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                ],
                StereoAtomForm::new(StereoKind::SquarePlanar, StereoCoset::Lit(0)),
            ),
        ],
        stereo_bonds: vec![
            (
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(0)),
            ),
            (
                BondId(1),
                vec![
                    StereoLigand::new(AtomId(0), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(2), StereoLigandKind::LonePair),
                ],
                StereoBondForm::new(StereoKind::CisTrans, StereoCoset::Lit(1)),
            ),
        ],
        ..Default::default()
    })
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
