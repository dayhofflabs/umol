//! Aggregate-canonicalization properties.
//!
//! The generated domain contains integrity-valid molecules with every entity kind and optional
//! constraints. Independently shuffled complete permutations in every entity namespace supply the
//! dense-remapping action. Exact fixtures and bounded exhaustive minima remain in the unit suite;
//! this module asserts the full normalization/reframe/canonicalize fixpoint and absorption matrix,
//! remapping invariance, equality, and canonical-hash laws without selecting a particular
//! symmetry-equivalent correspondence. Focused strategies separately cover feature-free,
//! constitution-bearing, and structure-bearing molecules, plus adjacent description-level
//! operands for binary equality.

use std::hash::{DefaultHasher, Hash, Hasher};

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_chem::element::Element;
use umol_graph_core::{AutomorphismAlgorithm, EdgeId, GraphRemapping, NodeId, Remapping};
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomForm, AtomId, BondForm, Canonicalize, CanonicalizeContext, Molecule,
    MoleculeEntries, MoleculeRemapping, NoncovalentBondForm, NoncovalentBondKind, Normalize,
    Reframe, StereoAtomForm, StereoKind, StereoLigand, StereoLigandKind,
};

use crate::strategies::{
    edge_set_strategy, element_strategy, molecule_dense_renumbering_strategy,
    standardization_scenario_strategy, stereo_reframed_molecule_pair_strategy,
};

fn context() -> CanonicalizeContext {
    CanonicalizeContext {
        para_stereo: false,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    }
}

fn structural_hash<T: Hash>(value: &T) -> u64 {
    let mut hasher = DefaultHasher::new();
    value.hash(&mut hasher);
    hasher.finish()
}

fn dense_renumbering_strategy(molecule: Molecule) -> BoxedStrategy<(Molecule, MoleculeRemapping)> {
    (
        Just(molecule.clone()),
        Just(molecule.atoms().ids().collect::<Vec<_>>()).prop_shuffle(),
        Just(molecule.bonds().ids().collect::<Vec<_>>()).prop_shuffle(),
        Just(molecule.dative_bonds().ids().collect::<Vec<_>>()).prop_shuffle(),
        Just(molecule.aromatic_systems().ids().collect::<Vec<_>>()).prop_shuffle(),
        Just(molecule.multicenter_bonds().ids().collect::<Vec<_>>()).prop_shuffle(),
        Just(molecule.noncovalent_bonds().ids().collect::<Vec<_>>()).prop_shuffle(),
        Just(molecule.stereo_atoms().ids().collect::<Vec<_>>()).prop_shuffle(),
        Just(molecule.stereo_bonds().ids().collect::<Vec<_>>()).prop_shuffle(),
    )
        .prop_map(
            |(
                molecule,
                atoms,
                bonds,
                dative,
                aromatic,
                multicenter,
                noncovalent,
                stereo_atoms,
                stereo_bonds,
            )| {
                (
                    molecule,
                    MoleculeRemapping::new(
                        GraphRemapping::new(
                            Remapping::new(atoms.iter().copied().map(NodeId::from).collect())
                                .unwrap(),
                            Remapping::new(bonds.iter().copied().map(EdgeId::from).collect())
                                .unwrap(),
                        ),
                        Remapping::new(dative).unwrap(),
                        Remapping::new(aromatic).unwrap(),
                        Remapping::new(multicenter).unwrap(),
                        Remapping::new(noncovalent).unwrap(),
                        Remapping::new(stereo_atoms).unwrap(),
                        Remapping::new(stereo_bonds).unwrap(),
                    ),
                )
            },
        )
        .boxed()
}

fn feature_free_dense_renumbering_strategy() -> BoxedStrategy<(Molecule, MoleculeRemapping)> {
    (2usize..=6)
        .prop_flat_map(|atom_count| {
            (
                prop::collection::vec(element_strategy(), atom_count),
                edge_set_strategy(atom_count),
            )
        })
        .prop_map(|(elements, edges)| {
            Molecule::from_entries(MoleculeEntries {
                atoms: elements.into_iter().map(AtomForm::from_element).collect(),
                bonds: edges
                    .into_iter()
                    .map(|[first, second]| (AtomId(first), AtomId(second), BondForm::from_order(1)))
                    .collect(),
                ..Default::default()
            })
        })
        .prop_flat_map(dense_renumbering_strategy)
        .boxed()
}

fn partially_featured_dense_renumbering_strategy() -> BoxedStrategy<(Molecule, MoleculeRemapping)> {
    prop_oneof![
        (prop::collection::vec(element_strategy(), 3..=6), 1u8..=3).prop_map(
            |(elements, bond_order)| {
                Molecule::from_entries(MoleculeEntries {
                    atoms: elements.into_iter().map(AtomForm::from_element).collect(),
                    bonds: vec![(AtomId(0), AtomId(1), BondForm::from_order(bond_order))],
                    noncovalent: vec![(
                        [AtomId(1), AtomId(2)],
                        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                    )],
                    ..Default::default()
                })
            }
        ),
        (0..StereoKind::Tetrahedral.count() as u32).prop_map(|coset| {
            Molecule::from_entries(MoleculeEntries {
                atoms: [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
                    .into_iter()
                    .map(AtomForm::from_element)
                    .collect(),
                bonds: (1..=4)
                    .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
                    .collect(),
                stereo_atoms: vec![(
                    AtomId(0),
                    (1..=4)
                        .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                        .collect(),
                    StereoAtomForm::new(StereoKind::Tetrahedral, coset),
                )],
                ..Default::default()
            })
        }),
    ]
    .prop_flat_map(dense_renumbering_strategy)
    .boxed()
}

fn canonical_eq_description_level_strategy(
) -> BoxedStrategy<(Molecule, Molecule, MoleculeRemapping)> {
    let atom_forms = || {
        [Element::C, Element::F, Element::Cl, Element::Br, Element::I]
            .into_iter()
            .map(AtomForm::from_element)
            .collect::<Vec<_>>()
    };
    let bonds = || {
        (1..=4)
            .map(|ligand| (AtomId(0), AtomId(ligand), BondForm::from_order(1)))
            .collect::<Vec<_>>()
    };
    let stereo_atoms = |coset| {
        vec![(
            AtomId(0),
            (1..=4)
                .map(|ligand| StereoLigand::new(AtomId(ligand), StereoLigandKind::Atom))
                .collect(),
            StereoAtomForm::new(StereoKind::Tetrahedral, coset),
        )]
    };

    prop_oneof![
        Just((
            Molecule::from_entries(MoleculeEntries {
                atoms: atom_forms(),
                bonds: bonds(),
                ..Default::default()
            }),
            Molecule::from_entries(MoleculeEntries {
                atoms: atom_forms(),
                bonds: bonds(),
                noncovalent: vec![(
                    [AtomId(1), AtomId(2)],
                    NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                )],
                ..Default::default()
            }),
        )),
        (0..StereoKind::Tetrahedral.count() as u32).prop_map(move |coset| {
            (
                Molecule::from_entries(MoleculeEntries {
                    atoms: atom_forms(),
                    bonds: bonds(),
                    noncovalent: vec![(
                        [AtomId(1), AtomId(2)],
                        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                    )],
                    ..Default::default()
                }),
                Molecule::from_entries(MoleculeEntries {
                    atoms: atom_forms(),
                    bonds: bonds(),
                    noncovalent: vec![(
                        [AtomId(1), AtomId(2)],
                        NoncovalentBondForm::from_kind(NoncovalentBondKind::HydrogenBond),
                    )],
                    stereo_atoms: stereo_atoms(coset),
                    ..Default::default()
                }),
            )
        }),
        (0..StereoKind::Tetrahedral.count() as u32).prop_map(move |coset| {
            let mut constrained_atoms = atom_forms();
            constrained_atoms[0].constraints = AtomConstraintForm::valence(4).into();
            (
                Molecule::from_entries(MoleculeEntries {
                    atoms: atom_forms(),
                    bonds: bonds(),
                    stereo_atoms: stereo_atoms(coset),
                    ..Default::default()
                }),
                Molecule::from_entries(MoleculeEntries {
                    atoms: constrained_atoms,
                    bonds: bonds(),
                    stereo_atoms: stereo_atoms(coset),
                    ..Default::default()
                }),
            )
        }),
    ]
    .prop_flat_map(|(lower, higher)| {
        dense_renumbering_strategy(higher)
            .prop_map(move |(higher, renumbering)| (lower.clone(), higher, renumbering))
    })
    .boxed()
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_molecule_canonicalize(
        (molecule, renumbering) in molecule_dense_renumbering_strategy(),
    ) {
        let context = context();
        let renumbered = molecule.remap(&renumbering);
        let canonical = molecule.clone().canonicalize(&context);
        let renumbered_canonical = renumbered.canonicalize(&context);

        prop_assert_eq!(&renumbered_canonical, &canonical);
        if let Ok(canonical) = canonical {
            let (with_correspondence, correspondence) = molecule
                .clone()
                .canonicalize_with_remapping(&context)
                .expect("successful canonicalization returns its correspondence");
            let reframed = molecule
                .remap(&correspondence)
                .reframe()
                .expect("a canonical correspondence preserves molecule integrity");

            prop_assert_eq!(&with_correspondence, &canonical);
            prop_assert_eq!(&reframed, &canonical);
            prop_assert!(molecule.framed_eq_under(&canonical, &correspondence));
        }
    }

    #[test]
    fn test_molecule_canonicalize_description_level(
        (molecule, renumbering) in prop_oneof![
            feature_free_dense_renumbering_strategy(),
            partially_featured_dense_renumbering_strategy(),
        ],
    ) {
        let context = context();
        let renumbered = molecule.remap(&renumbering);
        let (canonical, correspondence) = molecule
            .clone()
            .canonicalize_with_remapping(&context)
            .map_err(|error| {
                TestCaseError::fail(format!("generated molecule did not canonicalize: {error}"))
            })?;
        let (renumbered_canonical, renumbered_correspondence) = renumbered
            .clone()
            .canonicalize_with_remapping(&context)
            .map_err(|error| {
                TestCaseError::fail(format!(
                    "renumbered molecule did not canonicalize: {error}"
                ))
            })?;

        prop_assert_eq!(&renumbered_canonical, &canonical);
        prop_assert_eq!(molecule.clone().canonicalize(&context), Ok(canonical.clone()));
        prop_assert_eq!(renumbered.clone().canonicalize(&context), Ok(canonical.clone()));
        prop_assert_eq!(
            molecule.clone().canonical_hash(&context),
            renumbered.clone().canonical_hash(&context),
        );
        prop_assert_eq!(
            molecule.remap(&correspondence).reframe(),
            Ok(canonical.clone()),
        );
        prop_assert_eq!(
            renumbered.remap(&renumbered_correspondence).reframe(),
            Ok(canonical),
        );
    }

    #[test]
    fn test_molecule_canonicalize_standardization(
        scenario in standardization_scenario_strategy(),
    ) {
        let context = context();
        let source = scenario.molecule;
        let normalized = source.clone().normalize().map_err(|_| {
            TestCaseError::fail("generated molecule is intrinsically contradictory")
        })?;
        let reframed = source.clone().reframe().map_err(|_| {
            TestCaseError::fail("generated molecule is intrinsically contradictory")
        })?;
        let canonical = source.clone().canonicalize(&context).map_err(|error| {
            TestCaseError::fail(format!("generated molecule did not canonicalize: {error}"))
        })?;

        prop_assert_eq!(normalized.clone().normalize(), Ok(normalized.clone()));
        prop_assert_eq!(reframed.clone().reframe(), Ok(reframed.clone()));
        prop_assert_eq!(canonical.clone().canonicalize(&context), Ok(canonical.clone()));
        prop_assert_eq!(normalized.clone().reframe(), Ok(reframed.clone()));
        prop_assert_eq!(reframed.clone().normalize(), Ok(reframed.clone()));
        prop_assert_eq!(normalized.canonicalize(&context), Ok(canonical.clone()));
        prop_assert_eq!(reframed.canonicalize(&context), Ok(canonical.clone()));
        prop_assert_eq!(canonical.clone().normalize(), Ok(canonical.clone()));
        prop_assert_eq!(canonical.clone().reframe(), Ok(canonical));
    }

    #[test]
    fn test_molecule_canonical_hash(
        (molecule, renumbering) in molecule_dense_renumbering_strategy(),
    ) {
        let context = context();
        let renumbered = molecule.remap(&renumbering);

        prop_assert_eq!(
            molecule.clone().canonical_hash(&context),
            renumbered.clone().canonical_hash(&context),
        );
        if let Ok(canonical) = molecule.clone().canonicalize(&context) {
            prop_assert_eq!(
                molecule.canonical_hash(&context),
                Ok(structural_hash(&canonical)),
            );
        }
    }

    #[test]
    fn test_molecule_canonicalize_reframed(
        (left, right) in stereo_reframed_molecule_pair_strategy(),
    ) {
        let context = context();

        prop_assert_eq!(
            right.clone().canonicalize(&context),
            left.clone().canonicalize(&context),
        );
        prop_assert!(left.canonical_eq(&right, &context));
        prop_assert_eq!(right.canonical_hash(&context), left.canonical_hash(&context));
    }

    #[test]
    fn test_molecule_canonical_eq(
        (molecule, renumbering) in molecule_dense_renumbering_strategy(),
    ) {
        let context = context();
        let renumbered = molecule.remap(&renumbering);
        let canonical = molecule.clone().canonicalize(&context);

        prop_assert!(molecule.canonical_eq(&molecule, &context));
        prop_assert!(molecule.canonical_eq(&renumbered, &context));
        prop_assert_eq!(
            molecule.canonical_eq(&renumbered, &context),
            renumbered.canonical_eq(&molecule, &context),
        );
        if let Ok(canonical) = canonical {
            prop_assert!(renumbered.canonical_eq(&canonical, &context));
            prop_assert!(molecule.canonical_eq(&canonical, &context));
        }
    }

    #[test]
    fn test_molecule_canonical_eq_description_level(
        (lower, higher, renumbering) in canonical_eq_description_level_strategy(),
    ) {
        let context = context();
        let renumbered_higher = higher.remap(&renumbering);

        prop_assert!(!lower.canonical_eq(&renumbered_higher, &context));
        prop_assert!(!renumbered_higher.canonical_eq(&lower, &context));
    }

}
