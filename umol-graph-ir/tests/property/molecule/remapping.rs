//! Dense molecule-remapping properties.
//!
//! The generated domain contains every entity family, position-sensitive aromatic and
//! multicenter data, stereo frames, and constraints that reference every id family. Both generated
//! correspondences are complete non-identity cyclic permutations in every family. The success
//! properties use the asserted producer route; independently supplied coverage failures remain in
//! the exact unit suite for `try_remap`.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::Correspondence;

use crate::strategies::*;

fn identity<Id>(count: usize) -> Correspondence<Id>
where
    Id: Copy + Ord + From<usize>,
{
    Correspondence::from_images(&(0..count).map(Id::from).collect::<Vec<_>>(), count)
}

fn crossing<Id>(count: usize, seed: u64, family: u32) -> Correspondence<Id>
where
    Id: Copy + Ord + From<usize>,
{
    let shift = 1 + seed.rotate_right(family * 8) as usize % (count - 1);
    let images = (0..count)
        .map(|left| Id::from((left + shift) % count))
        .collect::<Vec<_>>();
    Correspondence::from_images(&images, count)
}

fn identity_correspondence(molecule: &Molecule) -> MoleculeCorrespondence {
    MoleculeCorrespondence::new(
        identity::<AtomId>(molecule.atoms().count()),
        identity::<BondId>(molecule.bonds().count()),
        identity::<DativeBondId>(molecule.dative_bonds().count()),
        identity::<AromaticSystemId>(molecule.aromatic_systems().count()),
        identity::<MulticenterBondId>(molecule.multicenter_bonds().count()),
        identity::<NoncovalentBondId>(molecule.noncovalent_bonds().count()),
        identity::<StereoAtomId>(molecule.stereo_atoms().count()),
        identity::<StereoBondId>(molecule.stereo_bonds().count()),
    )
}

fn crossing_correspondence(molecule: &Molecule, seed: u64) -> MoleculeCorrespondence {
    MoleculeCorrespondence::new(
        crossing::<AtomId>(molecule.atoms().count(), seed, 0),
        crossing::<BondId>(molecule.bonds().count(), seed, 1),
        crossing::<DativeBondId>(molecule.dative_bonds().count(), seed, 2),
        crossing::<AromaticSystemId>(molecule.aromatic_systems().count(), seed, 3),
        crossing::<MulticenterBondId>(molecule.multicenter_bonds().count(), seed, 4),
        crossing::<NoncovalentBondId>(molecule.noncovalent_bonds().count(), seed, 5),
        crossing::<StereoAtomId>(molecule.stereo_atoms().count(), seed, 6),
        crossing::<StereoBondId>(molecule.stereo_bonds().count(), seed, 7),
    )
}

fn atom(element: Element, charge: i64) -> AtomForm {
    let mut atom = AtomForm::from_element(element);
    atom.charge = NumForm::Lit(charge);
    atom
}

fn ligands(atoms: [u32; 4]) -> Vec<StereoLigand> {
    atoms
        .into_iter()
        .map(|atom| StereoLigand::new(AtomId(atom), StereoLigandKind::Atom))
        .collect()
}

fn remapping_molecule(atom_charge: i64, aromatic: i64, multicenter: i64) -> Molecule {
    Molecule::from_entries(MoleculeEntries {
        atoms: vec![
            atom(Element::C, atom_charge),
            atom(Element::N, atom_charge + 1),
            atom(Element::O, atom_charge + 2),
            atom(Element::F, atom_charge + 3),
            atom(Element::P, atom_charge + 4),
            atom(Element::S, atom_charge + 5),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(2)),
            (AtomId(2), AtomId(3), BondForm::from_order(3)),
            (AtomId(3), AtomId(4), BondForm::from_order(4)),
            (AtomId(4), AtomId(5), BondForm::from_order(5)),
            (AtomId(0), AtomId(5), BondForm::from_order(6)),
        ],
        dative: vec![
            (
                vec![AtomId(0), AtomId(2)],
                AtomId(1),
                DativeBondForm::from_order(1),
            ),
            (
                vec![AtomId(1), AtomId(3)],
                AtomId(4),
                DativeBondForm::from_order(2),
            ),
            (
                vec![AtomId(2), AtomId(4)],
                AtomId(5),
                DativeBondForm::from_order(3),
            ),
        ],
        aromatic: vec![
            (
                vec![AtomId(0), AtomId(2), AtomId(4)],
                AromaticSystemForm::from_electrons(vec![aromatic, aromatic + 1, aromatic + 2]),
            ),
            (
                vec![AtomId(1), AtomId(3), AtomId(5)],
                AromaticSystemForm::from_electrons(vec![aromatic + 3, aromatic + 4, aromatic + 5]),
            ),
            (
                vec![AtomId(0), AtomId(3), AtomId(5)],
                AromaticSystemForm::from_electrons(vec![aromatic + 6, aromatic + 7, aromatic + 8]),
            ),
        ],
        multicenter: vec![
            (
                vec![AtomId(0), AtomId(1), AtomId(4)],
                MulticenterBondForm::from_electrons(vec![
                    multicenter,
                    multicenter + 1,
                    multicenter + 2,
                ]),
            ),
            (
                vec![AtomId(1), AtomId(2), AtomId(5)],
                MulticenterBondForm::from_electrons(vec![
                    multicenter + 3,
                    multicenter + 4,
                    multicenter + 5,
                ]),
            ),
            (
                vec![AtomId(2), AtomId(3), AtomId(4)],
                MulticenterBondForm::from_electrons(vec![
                    multicenter + 6,
                    multicenter + 7,
                    multicenter + 8,
                ]),
            ),
        ],
        noncovalent: vec![
            (
                AtomId(0),
                AtomId(3),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            ),
            (
                AtomId(1),
                AtomId(4),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HalogenBond),
            ),
            (
                AtomId(2),
                AtomId(5),
                NoncovalentBondForm::from_kind(NoncovalentBondKind::Ionic),
            ),
        ],
        stereo_atoms: vec![
            (
                AtomId(0),
                ligands([1, 2, 3, 4]),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
            (
                AtomId(1),
                ligands([0, 2, 4, 5]),
                StereoAtomForm::new(StereoKind::Tetrahedral, 1u32),
            ),
            (
                AtomId(2),
                ligands([0, 1, 3, 5]),
                StereoAtomForm::new(StereoKind::Tetrahedral, 0u32),
            ),
        ],
        stereo_bonds: vec![
            (
                BondId(0),
                ligands([5, 2, 0, 1]),
                StereoBondForm::new(StereoKind::CisTrans, 0u32),
            ),
            (
                BondId(1),
                ligands([0, 3, 1, 2]),
                StereoBondForm::new(StereoKind::CisTrans, 1u32),
            ),
            (
                BondId(2),
                ligands([1, 4, 2, 3]),
                StereoBondForm::new(StereoKind::CisTrans, 0u32),
            ),
        ],
        constraints: Constraints::from_iter([
            Constraint::Relational(RelationalConstraint::DativeBondParallels {
                dative: DativeBondId(2),
                parallel: BondId(4),
            }),
            Constraint::Relational(RelationalConstraint::AromaticSystemContains {
                system: AromaticSystemId(2),
                atom: AtomId(5),
            }),
            Constraint::Relational(RelationalConstraint::MulticenterBondContains {
                bond: MulticenterBondId(2),
                atom: AtomId(4),
            }),
            Constraint::Relational(RelationalConstraint::NoncovalentBondContains {
                bond: NoncovalentBondId(2),
                atom: AtomId(5),
            }),
            Constraint::Relational(RelationalConstraint::StereoAtomSite {
                stereo_atom: StereoAtomId(2),
                atom: AtomId(2),
            }),
            Constraint::Relational(RelationalConstraint::StereoBondSite {
                stereo_bond: StereoBondId(2),
                bond: BondId(2),
            }),
        ]),
    })
}

fn remapping_scenario_strategy(
) -> impl Strategy<Value = (Molecule, MoleculeCorrespondence, MoleculeCorrespondence)> {
    (-2i64..=2, 0i64..=2, 0i64..=2, any::<u64>(), any::<u64>()).prop_map(
        |(atom_charge, aromatic, multicenter, first_seed, second_seed)| {
            let molecule = remapping_molecule(atom_charge, aromatic, multicenter);
            let first = crossing_correspondence(&molecule, first_seed);
            let second = crossing_correspondence(&molecule, second_seed);
            (molecule, first, second)
        },
    )
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_molecule_remap_equiv_under(
        (source, correspondence, _) in remapping_scenario_strategy(),
    ) {
        let remapped = source.remap(&correspondence);

        prop_assert!(source.equiv_under(&remapped, &correspondence));
    }

    #[test]
    fn test_molecule_remap_identity(
        (source, _, _) in remapping_scenario_strategy(),
    ) {
        let identity = identity_correspondence(&source);

        prop_assert_eq!(source.remap(&identity), source);
    }

    #[test]
    fn test_molecule_remap_inverse(
        (source, correspondence, _) in remapping_scenario_strategy(),
    ) {
        let remapped = source.remap(&correspondence);
        let restored = remapped.remap(&correspondence.reverse());

        prop_assert_eq!(restored, source);
    }

    #[test]
    fn test_molecule_remap_composition(
        (source, first, second) in remapping_scenario_strategy(),
    ) {
        let sequential = source.remap(&first).remap(&second);
        let direct = source.remap(&first.compose(&second));

        prop_assert_eq!(sequential, direct);
    }

    #[test]
    fn test_molecule_remap_integrity(
        (source, correspondence, _) in remapping_scenario_strategy(),
    ) {
        let remapped = source.remap(&correspondence);

        prop_assert_eq!(remapped.check_integrity(), Ok(()));
    }
}
