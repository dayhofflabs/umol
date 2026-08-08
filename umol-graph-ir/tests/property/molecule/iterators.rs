//! Exact-size properties of molecule and entity-view iterators.

use proptest::prelude::*;
use proptest::test_runner::{FileFailurePersistence, TestCaseResult};
use umol_graph_ir::ir::{RingConfig, RingModel};

use super::REGRESSION_FILE;
use crate::strategies::*;

fn assert_exact_prefix<T>(
    mut iterator: impl ExactSizeIterator<Item = T>,
    expected_len: usize,
    prefix: usize,
) -> TestCaseResult {
    let prefix = prefix.min(expected_len);
    prop_assert_eq!(iterator.len(), expected_len);
    prop_assert_eq!(iterator.size_hint(), (expected_len, Some(expected_len)),);
    for consumed in 0..prefix {
        prop_assert!(iterator.next().is_some());
        let remaining = expected_len - consumed - 1;
        prop_assert_eq!(iterator.len(), remaining);
        prop_assert_eq!(iterator.size_hint(), (remaining, Some(remaining)));
    }
    for consumed in prefix..expected_len {
        prop_assert!(iterator.next().is_some());
        let remaining = expected_len - consumed - 1;
        prop_assert_eq!(iterator.len(), remaining);
        prop_assert_eq!(iterator.size_hint(), (remaining, Some(remaining)));
    }
    prop_assert!(iterator.next().is_none());
    prop_assert_eq!(iterator.len(), 0);
    Ok(())
}

proptest! {
    #![proptest_config(ProptestConfig {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(REGRESSION_FILE))),
        ..ProptestConfig::default()
    })]

    #[test]
    fn test_molecule_ast_view_iterators_exact_size(
        molecule in molecule_ast_strategy(),
        prefix in any::<usize>(),
    ) {
        let atom_count = molecule.atoms().count();
        assert_exact_prefix(molecule.atoms().ids(), atom_count, prefix)?;
        assert_exact_prefix(molecule.atoms().iter(), atom_count, prefix)?;

        let bond_count = molecule.bonds().count();
        assert_exact_prefix(molecule.bonds().ids(), bond_count, prefix)?;
        assert_exact_prefix(molecule.bonds().iter(), bond_count, prefix)?;

        let dative_count = molecule.dative_bonds().count();
        assert_exact_prefix(molecule.dative_bonds().ids(), dative_count, prefix)?;
        assert_exact_prefix(molecule.dative_bonds().iter(), dative_count, prefix)?;

        let aromatic_count = molecule.aromatic_systems().count();
        assert_exact_prefix(molecule.aromatic_systems().ids(), aromatic_count, prefix)?;
        assert_exact_prefix(molecule.aromatic_systems().iter(), aromatic_count, prefix)?;

        let multicenter_count = molecule.multicenter_bonds().count();
        assert_exact_prefix(
            molecule.multicenter_bonds().ids(),
            multicenter_count,
            prefix,
        )?;
        assert_exact_prefix(
            molecule.multicenter_bonds().iter(),
            multicenter_count,
            prefix,
        )?;

        let noncovalent_count = molecule.noncovalent_bonds().count();
        assert_exact_prefix(
            molecule.noncovalent_bonds().ids(),
            noncovalent_count,
            prefix,
        )?;
        assert_exact_prefix(
            molecule.noncovalent_bonds().iter(),
            noncovalent_count,
            prefix,
        )?;

        let stereo_atom_count = molecule.stereo_atoms().count();
        assert_exact_prefix(molecule.stereo_atoms().ids(), stereo_atom_count, prefix)?;
        assert_exact_prefix(molecule.stereo_atoms().iter(), stereo_atom_count, prefix)?;

        let stereo_bond_count = molecule.stereo_bonds().count();
        assert_exact_prefix(molecule.stereo_bonds().ids(), stereo_bond_count, prefix)?;
        assert_exact_prefix(molecule.stereo_bonds().iter(), stereo_bond_count, prefix)?;

        for atom in molecule.atoms().ids() {
            let neighbor_count = molecule.neighbors(atom).count();
            assert_exact_prefix(molecule.neighbors(atom), neighbor_count, prefix)?;

            let view = molecule.atom(atom);
            assert_exact_prefix(view.neighbors(), neighbor_count, prefix)?;
            assert_exact_prefix(view.bond_ids(), neighbor_count, prefix)?;

            let dative_incidence = molecule.dative_bonds().incident_ids(atom).count();
            assert_exact_prefix(
                molecule.dative_bonds().incident_ids(atom),
                dative_incidence,
                prefix,
            )?;
            assert_exact_prefix(
                molecule.dative_bonds().incident(atom),
                dative_incidence,
                prefix,
            )?;
            assert_exact_prefix(view.dative_bond_ids(), dative_incidence, prefix)?;
            assert_exact_prefix(view.dative_bonds(), dative_incidence, prefix)?;

            let aromatic_incidence = molecule.aromatic_systems().incident_ids(atom).count();
            assert_exact_prefix(
                molecule.aromatic_systems().incident_ids(atom),
                aromatic_incidence,
                prefix,
            )?;
            assert_exact_prefix(
                molecule.aromatic_systems().incident(atom),
                aromatic_incidence,
                prefix,
            )?;

            let multicenter_incidence = molecule.multicenter_bonds().incident_ids(atom).count();
            assert_exact_prefix(
                molecule.multicenter_bonds().incident_ids(atom),
                multicenter_incidence,
                prefix,
            )?;
            assert_exact_prefix(
                molecule.multicenter_bonds().incident(atom),
                multicenter_incidence,
                prefix,
            )?;
            assert_exact_prefix(
                view.multicenter_bond_ids(),
                multicenter_incidence,
                prefix,
            )?;
            assert_exact_prefix(view.multicenter_bonds(), multicenter_incidence, prefix)?;

            let noncovalent_incidence = molecule.noncovalent_bonds().incident_ids(atom).count();
            assert_exact_prefix(
                molecule.noncovalent_bonds().incident_ids(atom),
                noncovalent_incidence,
                prefix,
            )?;
            assert_exact_prefix(
                molecule.noncovalent_bonds().incident(atom),
                noncovalent_incidence,
                prefix,
            )?;
            assert_exact_prefix(
                view.noncovalent_bond_ids(),
                noncovalent_incidence,
                prefix,
            )?;
            assert_exact_prefix(view.noncovalent_bonds(), noncovalent_incidence, prefix)?;

            let stereo_atom_incidence = molecule.stereo_atoms().incident_ids(atom).count();
            assert_exact_prefix(
                molecule.stereo_atoms().incident_ids(atom),
                stereo_atom_incidence,
                prefix,
            )?;
            assert_exact_prefix(
                molecule.stereo_atoms().incident(atom),
                stereo_atom_incidence,
                prefix,
            )?;
        }

        for bond in molecule.bonds().iter() {
            assert_exact_prefix(bond.atoms(), 2, prefix)?;
        }

        for dative in molecule.dative_bonds().iter() {
            assert_exact_prefix(dative.donor_ids(), dative.donor_count(), prefix)?;
            assert_exact_prefix(dative.donors(), dative.donor_count(), prefix)?;
            assert_exact_prefix(dative.atom_ids(), dative.atom_count(), prefix)?;
            assert_exact_prefix(dative.atoms(), dative.atom_count(), prefix)?;
        }

        for aromatic in molecule.aromatic_systems().iter() {
            assert_exact_prefix(aromatic.atom_ids(), aromatic.atom_count(), prefix)?;
            assert_exact_prefix(aromatic.atoms(), aromatic.atom_count(), prefix)?;
        }

        for multicenter in molecule.multicenter_bonds().iter() {
            assert_exact_prefix(multicenter.atom_ids(), multicenter.atom_count(), prefix)?;
            assert_exact_prefix(multicenter.atoms(), multicenter.atom_count(), prefix)?;
        }

        for stereo_atom in molecule.stereo_atoms().iter() {
            assert_exact_prefix(
                stereo_atom.ligands(),
                stereo_atom.ligand_count(),
                prefix,
            )?;
        }

        for stereo_bond in molecule.stereo_bonds().iter() {
            assert_exact_prefix(
                stereo_bond.ligands(),
                stereo_bond.ligand_count(),
                prefix,
            )?;
        }

        let rings = molecule.rings(RingModel::default(), RingConfig::default());
        let ring_count = rings.count();
        assert_exact_prefix(rings.ids(), ring_count, prefix)?;
        assert_exact_prefix(rings.iter(), ring_count, prefix)?;
        let ring_set = rings.into_ring_set();
        assert_exact_prefix(ring_set.ids(), ring_count, prefix)?;
        assert_exact_prefix(ring_set.iter(), ring_count, prefix)?;
    }
}
