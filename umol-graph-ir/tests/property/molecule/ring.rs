//! Molecule ring-enumeration properties.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_ir::ir::{RingConfig, RingModel, RingSetKind, ValueAst};

use crate::strategies::*;

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_molecule_ast_rings(
        ast in molecule_ast_strategy(),
        max_ring_size in 0usize..12,
        relevant in any::<bool>(),
    ) {
        let kind = if relevant {
            RingSetKind::Relevant
        } else {
            RingSetKind::Simple
        };
        let rings = ast.rings(
            RingModel {
                kind,
                max_ring_size,
            },
            RingConfig::default(),
        );

        prop_assert_eq!(rings.kind(), kind);
        prop_assert_eq!(rings.max_ring_size(), max_ring_size);
        prop_assert_eq!(
            rings.ids().collect::<Vec<_>>(),
            rings.iter().map(|ring| ring.id).collect::<Vec<_>>()
        );
        for ring in rings.iter() {
            prop_assert!(ring.len() <= max_ring_size);
            prop_assert_eq!(rings.get(ring.id).map(|view| view.atoms()), Some(ring.atoms()));
            prop_assert_eq!(rings.get(ring.id).map(|view| view.bonds()), Some(ring.bonds()));
            for &atom in ring.atoms() {
                prop_assert!(rings.atom(atom).is_in_ring());
            }
            for &bond in ring.bonds() {
                prop_assert!(rings.bond(bond).is_in_ring());
            }
        }
        for atom in ast.atoms().iter() {
            let view = rings.atom(atom.id);
            let containing: Vec<_> = view.rings().collect();
            prop_assert_eq!(view.is_in_ring(), !containing.is_empty());
            prop_assert_eq!(view.ring_count(), ValueAst::Lit(containing.len() as i64));
            for size in 0..=max_ring_size.min(u8::MAX as usize) {
                let expected = containing
                    .iter()
                    .filter(|ring| ring.len() == size)
                    .count();
                prop_assert_eq!(
                    view.ring_size_count(size as u8),
                    ValueAst::Lit(expected as i64)
                );
            }
        }
        for bond in ast.bonds().iter() {
            let view = rings.bond(bond.id);
            let containing: Vec<_> = view.rings().collect();
            prop_assert_eq!(view.is_in_ring(), !containing.is_empty());
            prop_assert_eq!(view.ring_count(), ValueAst::Lit(containing.len() as i64));
            for size in 0..=max_ring_size.min(u8::MAX as usize) {
                let expected = containing
                    .iter()
                    .filter(|ring| ring.len() == size)
                    .count();
                prop_assert_eq!(
                    view.ring_size_count(size as u8),
                    ValueAst::Lit(expected as i64)
                );
            }
        }
    }

    #[test]
    fn test_molecule_ast_rings_reindexing(
        ast in molecule_ast_strategy(),
        max_ring_size in 3usize..12,
        relevant in any::<bool>(),
    ) {
        let kind = if relevant {
            RingSetKind::Relevant
        } else {
            RingSetKind::Simple
        };
        let model = RingModel {
            kind,
            max_ring_size,
        };
        let config = RingConfig::default();
        let mut expected: Vec<_> = ast
            .rings(model, config)
            .iter()
            .map(|ring| {
                let mut atoms = ring.atoms().to_vec();
                let mut bonds = ring.bonds().to_vec();
                atoms.sort_unstable();
                bonds.sort_unstable();
                (atoms, bonds)
            })
            .collect();
        expected.sort_unstable();

        let mut actual = Vec::new();
        for (component, correspondence) in ast.split() {
            for ring in component.rings(model, config).iter() {
                let mut atoms = ring
                    .atoms()
                    .iter()
                    .map(|&atom| {
                        correspondence
                            .atoms()
                            .right_of(atom)
                            .ok_or_else(|| {
                                TestCaseError::fail(format!(
                                    "component atom {atom:?} has no source atom"
                                ))
                            })
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let mut bonds = ring
                    .bonds()
                    .iter()
                    .map(|&bond| {
                        correspondence.bonds().right_of(bond).ok_or_else(|| {
                            TestCaseError::fail(format!(
                                "component bond {bond:?} has no source bond"
                            ))
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                atoms.sort_unstable();
                bonds.sort_unstable();
                actual.push((atoms, bonds));
            }
        }
        actual.sort_unstable();

        prop_assert_eq!(actual, expected);
    }
}
