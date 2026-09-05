//! Molecule comparison properties.
//!
//! The identity-frame laws for `normalized_eq`, the correspondence laws for `framed_eq_under`,
//! and the complete comparison ladder deliberately use overlapping molecule domains. They
//! establish distinct relations: normalized equality in a shared frame, framed equality under an
//! explicit entity-id mapping, and canonical equality after entity-id renumbering. Successful and
//! intrinsically contradictory inputs have separate relation-law properties because only the
//! successful domain produces a canonical correspondence witness.

use std::fmt::Debug;

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::{AutomorphismAlgorithm, EdgeId, GraphRemapping, NodeId, Remapping};
use umol_graph_ir::ir::{
    AromaticSystemId, AtomId, Canonicalize, CanonicalizeContext, Contradiction, DativeBondId,
    Molecule, MoleculeCanonicalizeError, MoleculeEntries, MoleculeRemapping, MulticenterBondId,
    NoncovalentBondId, Normalize, NumForm, Reframe, StereoAtomId, StereoBondId,
};

use crate::strategies::{
    atom_form_strategy, intrinsic_contradiction_scenario_strategy, molecule_strategy,
    molecule_with_constraints_strategy, standardization_scenario_strategy,
    stereo_reframed_molecule_pair_strategy,
};

fn context() -> CanonicalizeContext {
    CanonicalizeContext {
        para_stereo: false,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    }
}

fn identity_remapping(molecule: &Molecule) -> MoleculeRemapping {
    fn identity<Id>(count: usize) -> Remapping<Id>
    where
        Id: Copy + Debug + Into<usize> + From<usize>,
    {
        let images: Vec<Id> = (0..count).map(Id::from).collect();
        Remapping::new(images).unwrap()
    }

    MoleculeRemapping::new(
        GraphRemapping::new(
            identity::<NodeId>(molecule.atoms().count()),
            identity::<EdgeId>(molecule.bonds().count()),
        ),
        identity::<DativeBondId>(molecule.dative_bonds().count()),
        identity::<AromaticSystemId>(molecule.aromatic_systems().count()),
        identity::<MulticenterBondId>(molecule.multicenter_bonds().count()),
        identity::<NoncovalentBondId>(molecule.noncovalent_bonds().count()),
        identity::<StereoAtomId>(molecule.stereo_atoms().count()),
        identity::<StereoBondId>(molecule.stereo_bonds().count()),
    )
}

fn atom_only_remapping(images: &[AtomId], _count: usize) -> MoleculeRemapping {
    fn empty<Id>() -> Remapping<Id>
    where
        Id: Copy + Debug + Into<usize> + From<usize>,
    {
        Remapping::default()
    }
    MoleculeRemapping::new(
        GraphRemapping::new(
            Remapping::new((images).iter().copied().map(NodeId::from).collect()).unwrap(),
            empty(),
        ),
        empty(),
        empty(),
        empty(),
        empty(),
        empty(),
        empty(),
    )
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
    fn test_molecule_normalized_eq_reflexive(molecule in molecule_with_constraints_strategy()) {
        prop_assert!(molecule.normalized_eq(&molecule));
    }

    #[test]
    fn test_molecule_normalized_eq_symmetric(
        left in molecule_with_constraints_strategy(),
        right in molecule_with_constraints_strategy(),
    ) {
        prop_assert_eq!(left.normalized_eq(&right), right.normalized_eq(&left));
    }

    #[test]
    fn test_molecule_framed_eq_under_composition(
        atoms in prop::collection::vec(atom_form_strategy(), 0..=5),
    ) {
        let count = atoms.len();
        let first_order = (0..count).collect::<Vec<_>>();
        let second_order = (0..count).rev().collect::<Vec<_>>();
        let mut third_order = first_order.clone();
        if count > 0 {
            third_order.rotate_left(1);
        }
        let molecule = |order: &[usize]| {
            Molecule::from_entries(MoleculeEntries {
                atoms: order.iter().map(|&index| atoms[index].clone()).collect(),
                ..Default::default()
            })
        };
        let correspondence = |left: &[usize], right: &[usize]| {
            let images = left
                .iter()
                .map(|original| {
                    AtomId::from(
                        right
                            .iter()
                            .position(|candidate| candidate == original)
                            .expect("orders contain the same indices"),
                    )
                })
                .collect::<Vec<_>>();
            atom_only_remapping(&images, count)
        };

        let first = molecule(&first_order);
        let second = molecule(&second_order);
        let third = molecule(&third_order);
        let first_second = correspondence(&first_order, &second_order);
        let second_third = correspondence(&second_order, &third_order);

        prop_assert!(first.framed_eq_under(&second, &first_second));
        prop_assert!(second.framed_eq_under(&third, &second_third));
        let composed = MoleculeRemapping::new(
            GraphRemapping::new(
                Remapping::new(composed_images(
                    first_second.graph().nodes(),
                    second_third.graph().nodes(),
                ))
                .unwrap(),
                Remapping::new(composed_images(
                    first_second.graph().edges(),
                    second_third.graph().edges(),
                ))
                .unwrap(),
            ),
            Remapping::new(composed_images(
                first_second.dative_bonds(),
                second_third.dative_bonds(),
            ))
            .unwrap(),
            Remapping::new(composed_images(
                first_second.aromatic_systems(),
                second_third.aromatic_systems(),
            ))
            .unwrap(),
            Remapping::new(composed_images(
                first_second.multicenter_bonds(),
                second_third.multicenter_bonds(),
            ))
            .unwrap(),
            Remapping::new(composed_images(
                first_second.noncovalent_bonds(),
                second_third.noncovalent_bonds(),
            ))
            .unwrap(),
            Remapping::new(composed_images(
                first_second.stereo_atoms(),
                second_third.stereo_atoms(),
            ))
            .unwrap(),
            Remapping::new(composed_images(
                first_second.stereo_bonds(),
                second_third.stereo_bonds(),
            ))
            .unwrap(),
        );
        prop_assert!(first.framed_eq_under(&third, &composed));
    }

    #[test]
    fn test_molecule_normalized_eq_agrees_with_equality_for_normalized_molecules(
        left in molecule_strategy(),
        right in molecule_strategy(),
    ) {
        prop_assert_eq!(left.normalized_eq(&right), left == right);
    }

    #[test]
    fn test_molecule_framed_eq_under_identity(
        molecule in molecule_with_constraints_strategy(),
    ) {
        let correspondence = identity_remapping(&molecule);
        let mut other = molecule.clone();
        if other.atoms().count() > 0 {
            other.atom_mut(AtomId(0)).attributes.charge = NumForm::Lit(99);
        }

        prop_assert_eq!(
            molecule.framed_eq_under(&other, &correspondence),
            molecule.framed_eq(&other),
        );
    }

    #[test]
    fn test_molecule_framed_eq_under_participant_frame(
        (left, right) in stereo_reframed_molecule_pair_strategy(),
    ) {
        let correspondence = identity_remapping(&left);

        prop_assert_eq!(
            left.framed_eq_under(&right, &correspondence),
            left.framed_eq(&right),
        );
    }

    #[test]
    fn test_molecule_framed_eq_under_inverse_correspondence(
        atoms in prop::collection::vec(atom_form_strategy(), 0..=5),
        change_mapped_atom in any::<bool>(),
    ) {
        let count = atoms.len();
        let left = Molecule::from_entries(MoleculeEntries {
            atoms: atoms.clone(),
            ..Default::default()
        });
        let mut right = Molecule::from_entries(MoleculeEntries {
            atoms: atoms.into_iter().rev().collect(),
            ..Default::default()
        });
        if change_mapped_atom && count > 0 {
            right.atom_mut(AtomId((count - 1) as u32)).attributes.charge = NumForm::Lit(99);
        }
        let images: Vec<AtomId> = (0..count).rev().map(AtomId::from).collect();
        let correspondence = MoleculeRemapping::new(
            GraphRemapping::new(
                Remapping::new(images.into_iter().map(NodeId::from).collect()).unwrap(),
                Remapping::default(),
            ),
            Remapping::default(),
            Remapping::default(),
            Remapping::default(),
            Remapping::default(),
            Remapping::default(),
            Remapping::default(),
        );

        let forward = left.framed_eq_under(&right, &correspondence);
        let reverse = right.framed_eq_under(&left, &correspondence);
        prop_assert_eq!(forward, reverse);
        prop_assert_eq!(forward, !change_mapped_atom || count == 0);
    }

    #[test]
    fn test_molecule_canonical_eq_standardization(
        scenario in standardization_scenario_strategy(),
    ) {
        let context = context();
        let source = scenario.molecule;
        let identical = source.clone();
        let normalized = source
            .clone()
            .normalize()
            .map_err(|_| TestCaseError::fail("generated molecule is intrinsically contradictory"))?;
        let normalized_again = normalized
            .clone()
            .normalize()
            .map_err(|_| TestCaseError::fail("normalized molecule became contradictory"))?;
        let reframed = source
            .clone()
            .reframe()
            .map_err(|_| TestCaseError::fail("generated molecule is intrinsically contradictory"))?;
        let canonical = source.clone().canonicalize(&context).map_err(|error| {
            TestCaseError::fail(format!("generated molecule did not canonicalize: {error}"))
        })?;
        let renumbered = source.remap(&scenario.remapping);

        prop_assert_eq!(&source, &identical);
        prop_assert!(source.normalized_eq(&identical));
        prop_assert!(source.framed_eq(&identical));
        prop_assert!(source.canonical_eq(&identical, &context));

        prop_assert!(source.normalized_eq(&normalized));
        prop_assert_eq!(
            source.normalized_eq(&normalized),
            normalized.normalized_eq(&source),
        );
        prop_assert!(normalized.normalized_eq(&normalized_again));
        prop_assert!(source.normalized_eq(&normalized_again));
        prop_assert!(source.framed_eq(&normalized));
        prop_assert!(normalized.framed_eq(&reframed));
        prop_assert!(source.framed_eq(&reframed));
        prop_assert_eq!(source.framed_eq(&reframed), reframed.framed_eq(&source),);

        prop_assert!(source.canonical_eq(&reframed, &context));
        prop_assert_eq!(
            source.canonical_eq(&reframed, &context),
            reframed.canonical_eq(&source, &context),
        );
        prop_assert!(reframed.canonical_eq(&renumbered, &context));
        prop_assert!(source.canonical_eq(&renumbered, &context));
        prop_assert!(renumbered.canonical_eq(&canonical, &context));
    }

    #[test]
    fn test_molecule_canonical_eq_contradiction(
        scenario in intrinsic_contradiction_scenario_strategy(),
    ) {
        let context = context();
        let [first, second, third] = scenario.molecules;

        prop_assert_ne!(&first, &second);
        prop_assert_ne!(&second, &third);
        prop_assert!(first.normalized_eq(&second));
        prop_assert!(second.normalized_eq(&third));
        prop_assert!(first.normalized_eq(&third));
        prop_assert!(first.framed_eq(&second));
        prop_assert!(second.framed_eq(&third));
        prop_assert!(first.framed_eq(&third));
        prop_assert!(first.canonical_eq(&second, &context));
        prop_assert!(second.canonical_eq(&third, &context));
        prop_assert!(first.canonical_eq(&third, &context));
        prop_assert_eq!(
            first.tracked_canonicalize(&context),
            Err(MoleculeCanonicalizeError::Contradiction(Contradiction)),
        );
    }
}
