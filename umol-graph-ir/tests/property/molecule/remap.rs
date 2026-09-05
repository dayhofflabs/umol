//! Dense molecule-remapping properties.
//!
//! The generated domain contains every entity kind, position-sensitive aromatic and multicenter
//! data, stereo frames, and constraints that reference every id kind. Both generated
//! remappings are complete nonidentity cyclic permutations in every component. The success
//! properties use the asserted producer route; independently supplied coverage failures remain in
//! the exact unit suite for `try_remap`.

use std::fmt::Debug;

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::{EdgeId, GraphRemapping, NodeId, Remapping};
use umol_graph_ir::ir::MoleculeRemapping;

use crate::strategies::*;

fn crossing<Id>(count: usize, seed: u64, entity_kind: u32) -> Remapping<Id>
where
    Id: Copy + Debug + Into<usize> + From<usize>,
{
    let shift = 1 + seed.rotate_right(entity_kind * 8) as usize % (count - 1);
    let images = (0..count)
        .map(|left| Id::from((left + shift) % count))
        .collect::<Vec<_>>();
    Remapping::new(images).unwrap()
}

fn identity_remapping(molecule: &Molecule) -> MoleculeRemapping {
    MoleculeRemapping::new(
        GraphRemapping::identity(molecule.atoms().count(), molecule.bonds().count()),
        Remapping::identity(molecule.dative_bonds().count()),
        Remapping::identity(molecule.aromatic_systems().count()),
        Remapping::identity(molecule.multicenter_bonds().count()),
        Remapping::identity(molecule.noncovalent_bonds().count()),
        Remapping::identity(molecule.stereo_atoms().count()),
        Remapping::identity(molecule.stereo_bonds().count()),
    )
}

fn crossing_remapping(molecule: &Molecule, seed: u64) -> MoleculeRemapping {
    MoleculeRemapping::new(
        GraphRemapping::new(
            crossing::<NodeId>(molecule.atoms().count(), seed, 0),
            crossing::<EdgeId>(molecule.bonds().count(), seed, 1),
        ),
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
            atom(Element::Cl, atom_charge + 6),
            atom(Element::Br, atom_charge + 7),
            atom(Element::H, atom_charge + 8),
        ],
        bonds: vec![
            (AtomId(0), AtomId(1), BondForm::from_order(1)),
            (AtomId(1), AtomId(2), BondForm::from_order(2)),
            (AtomId(2), AtomId(3), BondForm::from_order(3)),
            (AtomId(3), AtomId(4), BondForm::from_order(4)),
            (AtomId(4), AtomId(5), BondForm::from_order(5)),
            (AtomId(0), AtomId(5), BondForm::from_order(6)),
            (AtomId(0), AtomId(2), BondForm::from_order(1)),
            (AtomId(0), AtomId(3), BondForm::from_order(1)),
            (AtomId(0), AtomId(4), BondForm::from_order(1)),
            (AtomId(1), AtomId(3), BondForm::from_order(1)),
            (AtomId(1), AtomId(4), BondForm::from_order(1)),
            (AtomId(1), AtomId(5), BondForm::from_order(1)),
            (AtomId(2), AtomId(5), BondForm::from_order(1)),
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
                vec![AtomId(6), AtomId(7), AtomId(8)],
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
                [AtomId(0), AtomId(3)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
            ),
            (
                [AtomId(1), AtomId(4)],
                NoncovalentBondForm::from_kind(NoncovalentBondKind::HalogenBond),
            ),
            (
                [AtomId(2), AtomId(5)],
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
                ligands([5, 2, 4, 3]),
                StereoBondForm::new(StereoKind::CisTrans, 0u32),
            ),
            (
                BondId(1),
                ligands([0, 4, 3, 5]),
                StereoBondForm::new(StereoKind::CisTrans, 1u32),
            ),
            (
                BondId(2),
                ligands([1, 5, 4, 0]),
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
) -> impl Strategy<Value = (Molecule, MoleculeRemapping, MoleculeRemapping)> {
    (-2i64..=2, 0i64..=2, 0i64..=2, any::<u64>(), any::<u64>()).prop_map(
        |(atom_charge, aromatic, multicenter, first_seed, second_seed)| {
            let molecule = remapping_molecule(atom_charge, aromatic, multicenter);
            let first = crossing_remapping(&molecule, first_seed);
            let second = crossing_remapping(&molecule, second_seed);
            (molecule, first, second)
        },
    )
}

fn inverse_images<Id: Copy + Into<usize> + From<usize>>(remapping: &Remapping<Id>) -> Vec<Id> {
    let mut pairs = (0..remapping.len())
        .map(Id::from)
        .map(|source| (Into::<usize>::into(remapping.map(source)), source))
        .collect::<Vec<_>>();
    pairs.sort_unstable_by_key(|&(target, _)| target);
    pairs.into_iter().map(|(_, source)| source).collect()
}

fn composed_images<Id: Copy + Into<usize> + From<usize>>(
    first: &Remapping<Id>,
    second: &Remapping<Id>,
) -> Vec<Id> {
    assert_eq!(first.len(), second.len());
    (0..first.len())
        .map(Id::from)
        .map(|source| second.map(first.map(source)))
        .collect()
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_molecule_remap_framed_eq_under(
        (source, remapping, _) in remapping_scenario_strategy(),
    ) {
        let remapped = source.remap(&remapping);

        prop_assert!(source.framed_eq_under(&remapped, &remapping));
    }

    #[test]
    fn test_molecule_remap_identity(
        (source, _, _) in remapping_scenario_strategy(),
    ) {
        let identity = identity_remapping(&source);

        prop_assert_eq!(source.remap(&identity), source);
    }

    #[test]
    fn test_molecule_remap_inverse(
        (source, remapping, _) in remapping_scenario_strategy(),
    ) {
        let remapped = source.remap(&remapping);
        let restored = remapped.remap(&MoleculeRemapping::new(
            GraphRemapping::new(
                Remapping::new(inverse_images(remapping.graph().nodes())).unwrap(),
                Remapping::new(inverse_images(remapping.graph().edges())).unwrap(),
            ),
            Remapping::new(inverse_images(remapping.dative_bonds())).unwrap(),
            Remapping::new(inverse_images(remapping.aromatic_systems())).unwrap(),
            Remapping::new(inverse_images(remapping.multicenter_bonds())).unwrap(),
            Remapping::new(inverse_images(remapping.noncovalent_bonds())).unwrap(),
            Remapping::new(inverse_images(remapping.stereo_atoms())).unwrap(),
            Remapping::new(inverse_images(remapping.stereo_bonds())).unwrap(),
        ));

        prop_assert_eq!(restored, source);
    }

    #[test]
    fn test_molecule_remap_composition(
        (source, first, second) in remapping_scenario_strategy(),
    ) {
        let sequential = source.remap(&first).remap(&second);
        let direct = source.remap(&MoleculeRemapping::new(
            GraphRemapping::new(
                Remapping::new(composed_images(
                    first.graph().nodes(),
                    second.graph().nodes(),
                ))
                .unwrap(),
                Remapping::new(composed_images(
                    first.graph().edges(),
                    second.graph().edges(),
                ))
                .unwrap(),
            ),
            Remapping::new(composed_images(first.dative_bonds(), second.dative_bonds())).unwrap(),
            Remapping::new(composed_images(
                first.aromatic_systems(),
                second.aromatic_systems(),
            ))
            .unwrap(),
            Remapping::new(composed_images(
                first.multicenter_bonds(),
                second.multicenter_bonds(),
            ))
            .unwrap(),
            Remapping::new(composed_images(
                first.noncovalent_bonds(),
                second.noncovalent_bonds(),
            ))
            .unwrap(),
            Remapping::new(composed_images(first.stereo_atoms(), second.stereo_atoms())).unwrap(),
            Remapping::new(composed_images(first.stereo_bonds(), second.stereo_bonds())).unwrap(),
        ));

        prop_assert_eq!(sequential, direct);
    }

}
