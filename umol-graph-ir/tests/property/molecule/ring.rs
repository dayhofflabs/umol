//! Molecule ring-enumeration properties.
//! Reconstruction compares owning and borrowing routes on generated compatible
//! molecules; independent context cases check the documented panic freedom.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_ir::ir::{AtomId, BondId, NumForm, RingConfig, RingModel, RingSetKind, RingViews};

use crate::strategies::molecule_strategy;

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_molecule_rings(
        molecule in molecule_strategy(),
        max_ring_size in 0usize..12,
        relevant in any::<bool>(),
    ) {
        let kind = if relevant {
            RingSetKind::Relevant
        } else {
            RingSetKind::Simple
        };
        let rings = molecule.rings(
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
        for atom in molecule.atoms().iter() {
            let view = rings.atom(atom.id());
            let containing: Vec<_> = view.rings().collect();
            prop_assert_eq!(view.is_in_ring(), !containing.is_empty());
            prop_assert_eq!(view.ring_count(), NumForm::Lit(containing.len() as i64));
            for size in 0..=max_ring_size.min(u8::MAX as usize) {
                let expected = containing
                    .iter()
                    .filter(|ring| ring.len() == size)
                    .count();
                prop_assert_eq!(
                    view.ring_size_count(size as u8),
                    NumForm::Lit(expected as i64)
                );
            }
        }
        for bond in molecule.bonds().iter() {
            let view = rings.bond(bond.id());
            let containing: Vec<_> = view.rings().collect();
            prop_assert_eq!(view.is_in_ring(), !containing.is_empty());
            prop_assert_eq!(view.ring_count(), NumForm::Lit(containing.len() as i64));
            for size in 0..=max_ring_size.min(u8::MAX as usize) {
                let expected = containing
                    .iter()
                    .filter(|ring| ring.len() == size)
                    .count();
                prop_assert_eq!(
                    view.ring_size_count(size as u8),
                    NumForm::Lit(expected as i64)
                );
            }
        }
    }

    #[test]
    fn test_molecule_rings_reindexing(
        molecule in molecule_strategy(),
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
        let mut expected: Vec<_> = molecule
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
        for (component, correspondence) in molecule.tracked_split() {
            for ring in component.rings(model, config).iter() {
                let mut atoms = ring
                    .atoms()
                    .iter()
                    .map(|&atom| {
                        correspondence
                            .atoms()
                            .left_of(atom)
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
                        correspondence.bonds().left_of(bond).ok_or_else(|| {
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

    #[test]
    fn test_ring_views_from_ring_set(
        molecule in molecule_strategy(),
        max_ring_size in 0usize..12,
        relevant in any::<bool>(),
    ) {
        let model = RingModel {
            kind: if relevant { RingSetKind::Relevant } else { RingSetKind::Simple },
            max_ring_size,
        };
        let owned = molecule.rings(model, RingConfig::default());
        let ring_set = molecule.rings(model, RingConfig::default()).into_ring_set();
        let borrowed = RingViews::from_ring_set(&ring_set, &molecule);
        prop_assert_eq!(borrowed.kind(), owned.kind());
        prop_assert_eq!(borrowed.max_ring_size(), owned.max_ring_size());
        prop_assert_eq!(borrowed.count(), owned.count());
        prop_assert_eq!(borrowed.ids().collect::<Vec<_>>(), owned.ids().collect::<Vec<_>>());
        prop_assert_eq!(
            borrowed.iter().map(|ring| (ring.id, ring.atoms(), ring.bonds())).collect::<Vec<_>>(),
            owned.iter().map(|ring| (ring.id, ring.atoms(), ring.bonds())).collect::<Vec<_>>()
        );
        for id in molecule.atoms().ids() {
            let actual = borrowed.atom(id);
            let expected = owned.atom(id);
            prop_assert_eq!(actual.is_in_ring(), expected.is_in_ring());
            prop_assert_eq!(actual.rings().map(|ring| ring.id).collect::<Vec<_>>(), expected.rings().map(|ring| ring.id).collect::<Vec<_>>());
            prop_assert_eq!(actual.ring_count(), expected.ring_count());
            prop_assert_eq!(actual.smallest_ring_size(), expected.smallest_ring_size());
            prop_assert_eq!(actual.ring_degree(), expected.ring_degree());
            prop_assert_eq!(actual.ring_valence(), expected.ring_valence());
            for size in 0..=max_ring_size {
                prop_assert_eq!(actual.ring_size_count(size as u8), expected.ring_size_count(size as u8));
            }
        }
        for id in molecule.bonds().ids() {
            let actual = borrowed.bond(id);
            let expected = owned.bond(id);
            prop_assert_eq!(actual.is_in_ring(), expected.is_in_ring());
            prop_assert_eq!(actual.rings().map(|ring| ring.id).collect::<Vec<_>>(), expected.rings().map(|ring| ring.id).collect::<Vec<_>>());
            prop_assert_eq!(actual.ring_count(), expected.ring_count());
            prop_assert_eq!(actual.smallest_ring_size(), expected.smallest_ring_size());
            for size in 0..=max_ring_size {
                prop_assert_eq!(actual.ring_size_count(size as u8), expected.ring_size_count(size as u8));
            }
        }
        let extracted = borrowed.into_ring_set();
        prop_assert_eq!(extracted.kind(), ring_set.kind());
        prop_assert_eq!(extracted.max_ring_size(), ring_set.max_ring_size());
        prop_assert_eq!(
            extracted.iter().map(|ring| (ring.id, ring.atoms(), ring.bonds())).collect::<Vec<_>>(),
            ring_set.iter().map(|ring| (ring.id, ring.atoms(), ring.bonds())).collect::<Vec<_>>()
        );
        let extracted_graph = extracted.graph();
        let source_graph = ring_set.graph();
        prop_assert_eq!(extracted_graph.edges(), source_graph.edges());
    }

    #[test]
    fn test_ring_views_from_ring_set_context(
        source in molecule_strategy(),
        context in molecule_strategy(),
    ) {
        let ring_set = source.rings(RingModel::default(), RingConfig::default()).into_ring_set();
        let views = RingViews::from_ring_set(&ring_set, &context);
        for id in source.atoms().ids().chain([AtomId(context.atoms().count() as u32)]) {
            let atom = views.atom(id);
            let _ = (atom.ring_degree(), atom.ring_valence());
        }
        for id in source.bonds().ids().chain([BondId(context.bonds().count() as u32)]) {
            let bond = views.bond(id);
            let _ = (bond.is_in_ring(), bond.ring_count(), bond.ring_size_count(3), bond.smallest_ring_size());
            let _ = bond.rings().collect::<Vec<_>>();
        }
        prop_assert_eq!(
            views.iter().map(|ring| (ring.id, ring.atoms(), ring.bonds())).collect::<Vec<_>>(),
            ring_set.iter().map(|ring| (ring.id, ring.atoms(), ring.bonds())).collect::<Vec<_>>()
        );
    }
}
