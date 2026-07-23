//! Molecule ring-enumeration properties.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_ast::ast::RingFamily;
use umol_graph_core::CycleEnumerationAlgorithm;

use crate::strategies::*;

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_molecule_ast_rings(ast in molecule_ast_strategy()) {
        let rings = ast
            .rings(CycleEnumerationAlgorithm::Vismara)
            .iter()
            .map(|ring| (ring.atoms().to_vec(), ring.bonds().to_vec()))
            .collect::<Vec<_>>();
        let explicit = ast
            .rings_with(
                RingFamily::Relevant,
                22,
                |_| true,
                CycleEnumerationAlgorithm::Vismara,
            )
            .iter()
            .map(|ring| (ring.atoms().to_vec(), ring.bonds().to_vec()))
            .collect::<Vec<_>>();

        prop_assert_eq!(rings, explicit);
    }

    #[test]
    fn test_molecule_ast_rings_with(
        ast in molecule_ast_strategy(),
        max_ring_size in 0usize..12,
        atom_cutoff in 0u32..16,
        relevant in any::<bool>(),
    ) {
        let family = if relevant {
            RingFamily::Relevant
        } else {
            RingFamily::Simple
        };
        let all = ast
            .rings_with(
                family,
                max_ring_size,
                |_| true,
                CycleEnumerationAlgorithm::Vismara,
            )
            .iter()
            .map(|ring| {
                let mut bonds = ring.bonds().to_vec();
                bonds.sort_unstable();
                bonds
            })
            .collect::<Vec<_>>();
        let filtered = ast.rings_with(
            family,
            max_ring_size,
            |atom| atom.0 < atom_cutoff,
            CycleEnumerationAlgorithm::Vismara,
        );

        for ring in filtered.iter() {
            let mut bonds = ring.bonds().to_vec();
            bonds.sort_unstable();

            prop_assert!(ring.len() <= max_ring_size);
            prop_assert!(ring.atoms().iter().all(|atom| atom.0 < atom_cutoff));
            prop_assert!(all.contains(&bonds));
        }
    }
}
