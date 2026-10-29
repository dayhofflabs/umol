//! Stored arrangement fixtures and an independent bounded assignment reference.

mod fixtures;
mod reference;

use fixtures::{
    cis_trans, evidenced_stereo, hydrogen_and_lone_pair, imine, implicit_hydrogen, lone_pair,
    partial_stereo, prochiral, remote_stereo, small_ring, square_planar, symmetric, tetrahedral,
};
use pretty_assertions::assert_eq;
use rstest::rstest;
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::StereoLigandKind::{Atom, ImplicitHydrogen, LonePair};
use umol_graph_ir::ir::{
    AtomId, BondId, CisTransStereoForm, ConstitutionColoring, ElementForm, Entity,
    GraphSymmetryConfig, Molecule, StereoAtomId, StereoBondId, StereoConfigurationForm,
    StereoCoset, StereoKind, StereoLigandKind, StereoTerm, Stereogenicity, TetrahedralStereoForm,
};

type StereoRecord = (
    Entity,
    Vec<(u32, StereoLigandKind)>,
    StereoConfigurationForm,
);

#[rstest]
#[case::tetrahedral(tetrahedral(), vec![
    (Entity::Atom(AtomId(0)), vec![(1, Atom), (2, Atom), (3, Atom), (4, Atom)],
     StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined)),
])]
#[case::prochiral(prochiral(), vec![
    (Entity::Atom(AtomId(0)), vec![(1, Atom), (2, Atom), (3, Atom), (4, Atom)],
     StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32)),
])]
#[case::symmetric(symmetric(), vec![
    (Entity::Atom(AtomId(0)), vec![(1, Atom), (2, Atom), (3, Atom), (4, Atom)],
     StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32)),
])]
#[case::implicit_hydrogen(implicit_hydrogen(), vec![
    (Entity::Atom(AtomId(0)), vec![(1, Atom), (2, Atom), (3, Atom), (0, ImplicitHydrogen)],
     StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined)),
])]
#[case::lone_pair(lone_pair(), vec![
    (Entity::Atom(AtomId(0)), vec![(1, Atom), (2, Atom), (3, Atom), (0, LonePair)],
     StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined)),
])]
#[case::hydrogen_and_lone_pair(hydrogen_and_lone_pair(), vec![
    (Entity::Atom(AtomId(0)), vec![(1, Atom), (2, Atom), (0, ImplicitHydrogen), (0, LonePair)],
     StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined)),
])]
#[case::cis_trans(cis_trans(), vec![
    (Entity::Bond(BondId(0)), vec![(2, Atom), (3, Atom), (4, Atom), (5, Atom)],
     StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Undetermined)),
])]
#[case::imine(imine(), vec![
    (Entity::Bond(BondId(0)), vec![(2, Atom), (0, ImplicitHydrogen), (3, Atom), (1, LonePair)],
     StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Undetermined)),
])]
#[case::remote_stereo(remote_stereo(), vec![
    (Entity::Atom(AtomId(0)), vec![(1, Atom), (2, Atom), (3, Atom), (4, Atom)],
     StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined)),
    (Entity::Atom(AtomId(1)), vec![(0, Atom), (5, Atom), (6, Atom), (7, Atom)],
     StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 0u32)),
    (Entity::Atom(AtomId(2)), vec![(0, Atom), (8, Atom), (9, Atom), (10, Atom)],
     StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 1u32)),
])]
#[case::small_ring(small_ring(), vec![
    (Entity::Bond(BondId(0)), vec![(5, Atom), (0, ImplicitHydrogen), (2, Atom), (1, ImplicitHydrogen)],
     StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Undetermined)),
])]
#[case::square_planar(square_planar(), vec![
    (Entity::Atom(AtomId(0)), vec![(1, Atom), (2, Atom), (3, Atom), (4, Atom)],
     StereoConfigurationForm::kinded(StereoKind::SquarePlanar, 0u32)),
])]
#[case::evidenced_stereo(evidenced_stereo(), vec![
    (Entity::Atom(AtomId(0)), vec![(1, Atom), (2, Atom), (3, Atom), (4, Atom)],
     StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined)),
    (Entity::Bond(BondId(4)), vec![(0, Atom), (1, ImplicitHydrogen), (6, Atom), (5, ImplicitHydrogen)],
     StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Undetermined)),
])]
#[case::partial_stereo(partial_stereo(), vec![
    (Entity::Atom(AtomId(0)), vec![(1, Atom), (2, Atom), (3, Atom), (4, Atom)],
     StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::lit_set([0, 1]))),
    (Entity::Bond(BondId(4)), vec![(0, Atom), (1, ImplicitHydrogen), (6, Atom), (5, ImplicitHydrogen)],
     StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::term(StereoTerm::var("configuration")))),
])]
fn test_molecule_edit_stereo_fixture(
    #[case] molecule: Molecule,
    #[case] expected: Vec<StereoRecord>,
) {
    let atoms = molecule.stereo_atoms().iter().map(|view| {
        (
            Entity::Atom(view.site_id()),
            view.ligand_ids()
                .iter()
                .map(|ligand| (ligand.atom_id.0, ligand.kind))
                .collect(),
            view.attributes().configuration.clone(),
        )
    });
    let bonds = molecule.stereo_bonds().iter().map(|view| {
        (
            Entity::Bond(view.site_id()),
            view.ligand_ids()
                .iter()
                .map(|ligand| (ligand.atom_id.0, ligand.kind))
                .collect(),
            view.attributes().configuration.clone(),
        )
    });
    assert_eq!(atoms.chain(bonds).collect::<Vec<_>>(), expected);
    assert_eq!(molecule.clone().edit().finish(), Ok(molecule));
}

#[rstest]
#[case::stereogenic(tetrahedral(), Stereogenicity::Stereogenic)]
#[case::prochiral(prochiral(), Stereogenicity::Prochiral)]
#[case::symmetric(symmetric(), Stereogenicity::Symmetric)]
fn test_stereo_symmetry_stereogenicity_literal(
    #[case] mut molecule: Molecule,
    #[case] expected: Stereogenicity,
    #[values(0u32, 1)] coset: u32,
) {
    molecule
        .stereo_atom_mut(StereoAtomId(0))
        .attributes_mut()
        .configuration = StereoConfigurationForm::kinded(StereoKind::Tetrahedral, coset);
    let graph = molecule.graph_symmetry(&GraphSymmetryConfig {
        coloring: ConstitutionColoring::full(),
        iterate_to_fixpoint: true,
        max_iterations: 16,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    });
    assert_eq!(
        molecule
            .stereo_atom_symmetry(&graph, StereoAtomId(0))
            .stereogenicity(),
        expected,
    );
}

#[rstest]
fn test_molecule_edit_partial_stereo(partial_stereo: Molecule) {
    assert_eq!(
        partial_stereo.atom(AtomId(7)).element(),
        &ElementForm::Undetermined
    );
    assert_eq!(partial_stereo.clone().edit().finish(), Ok(partial_stereo));
}

#[rstest]
fn test_absolute_assignments_stereo_evidence(evidenced_stereo: Molecule) {
    assert_eq!(
        evidenced_stereo
            .atom(AtomId(0))
            .constraints()
            .tetrahedral_stereo(),
        &TetrahedralStereoForm::Stereo(StereoCoset::Lit(1)),
    );
    assert_eq!(
        evidenced_stereo
            .bond(BondId(4))
            .constraints()
            .cis_trans_stereo(),
        &CisTransStereoForm::Stereo(StereoCoset::Lit(0)),
    );
    let source = [
        evidenced_stereo
            .stereo_atom(StereoAtomId(0))
            .attributes()
            .configuration
            .coset()
            .unwrap()
            .clone(),
        evidenced_stereo
            .stereo_bond(StereoBondId(0))
            .attributes()
            .configuration
            .coset()
            .unwrap()
            .clone(),
    ];
    assert_eq!(
        reference::absolute_assignments(&source, &[0, 1].into(), &[0, 1], |values| {
            values == [StereoCoset::Lit(1), StereoCoset::Lit(0)]
        }),
        Ok(vec![vec![StereoCoset::Lit(1), StereoCoset::Lit(0)]]),
    );
}
