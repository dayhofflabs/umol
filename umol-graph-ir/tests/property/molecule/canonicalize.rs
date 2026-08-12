//! Aggregate-canonicalization properties.
//!
//! The generated domain contains integrity-valid molecules with every entity family and optional
//! constraints. A complete reverse permutation in every entity namespace supplies the independent
//! dense-remapping action. Exact fixtures and bounded exhaustive minima remain in the unit suite;
//! this module asserts the public idempotence, remapping-invariance, and equality laws.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::{AutomorphismAlgorithm, Correspondence};
use umol_graph_ir::ir::{
    CanonicalizationContext, CanonicalizationLevel, Canonicalize, MoleculeCorrespondence,
};

use crate::strategies::*;

fn reverse<Id>(count: usize) -> Correspondence<Id>
where
    Id: Copy + Ord + From<usize>,
{
    Correspondence::from_images(&(0..count).rev().map(Id::from).collect::<Vec<_>>(), count)
}

fn reverse_correspondence(molecule: &Molecule) -> MoleculeCorrespondence {
    MoleculeCorrespondence::new(
        reverse::<AtomId>(molecule.atoms().count()),
        reverse::<BondId>(molecule.bonds().count()),
        reverse::<DativeBondId>(molecule.dative_bonds().count()),
        reverse::<AromaticSystemId>(molecule.aromatic_systems().count()),
        reverse::<MulticenterBondId>(molecule.multicenter_bonds().count()),
        reverse::<NoncovalentBondId>(molecule.noncovalent_bonds().count()),
        reverse::<StereoAtomId>(molecule.stereo_atoms().count()),
        reverse::<StereoBondId>(molecule.stereo_bonds().count()),
    )
}

fn context() -> CanonicalizationContext {
    CanonicalizationContext {
        para_stereo: false,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    }
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_molecule_canonicalize(molecule in molecule_with_constraints_strategy()) {
        let context = context();
        let renumbered = molecule.remap(&reverse_correspondence(&molecule));
        let canonical = molecule.clone().canonicalize(&context);
        let renumbered_canonical = renumbered.canonicalize(&context);

        prop_assert_eq!(&renumbered_canonical, &canonical);
        if let Ok(canonical) = canonical {
            prop_assert_eq!(canonical.clone().canonicalize(&context), Ok(canonical));
        }
    }

    #[test]
    fn test_molecule_canonical_eq(molecule in molecule_with_constraints_strategy()) {
        let context = context();
        let renumbered = molecule.remap(&reverse_correspondence(&molecule));
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
    fn test_molecule_canonical_eq_by(molecule in molecule_with_constraints_strategy()) {
        let context = context();
        let renumbered = molecule.remap(&reverse_correspondence(&molecule));

        for level in [
            CanonicalizationLevel::Topology,
            CanonicalizationLevel::Constitution,
            CanonicalizationLevel::Structure,
            CanonicalizationLevel::Full,
        ] {
            let canonical = molecule.clone().canonicalize_by(level, &context);

            prop_assert!(molecule.canonical_eq_by(&molecule, level, &context));
            prop_assert!(molecule.canonical_eq_by(&renumbered, level, &context));
            prop_assert_eq!(
                molecule.canonical_eq_by(&renumbered, level, &context),
                renumbered.canonical_eq_by(&molecule, level, &context),
            );
            if let Ok(canonical) = canonical {
                prop_assert!(renumbered.canonical_eq_by(&canonical, level, &context));
                prop_assert!(molecule.canonical_eq_by(&canonical, level, &context));
            }
        }
        prop_assert_eq!(
            molecule.canonical_eq_by(&renumbered, CanonicalizationLevel::Full, &context),
            molecule.canonical_eq(&renumbered, &context),
        );
    }
}
