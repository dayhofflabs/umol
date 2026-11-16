//! Reading normality, normalization invariance, and meet composition on raw forms.
//! Non-symmetry readings and stereo own assertions are compared across normalization.
//! Stereo derived readings use raw forms and compose with assertions through meet.
//! Numeric cases independently normalize sources before computing their meet.
//! Container/view equality checks own assertions except for kindless stereo entities,
//! whose assertions are independently required to be top despite their raw payload;
//! container unit cases independently cover raw contradictions and closure of absence.
//! Reframing invariance checks atom and bond stereo projections on molecules whose stored
//! participant orders change, independently of the fixed-frame normalization law.

use proptest::prelude::*;
use umol_chem::element::Element;
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomConstraintKey, AtomForm, AtomId, BondConstraintKey, BondForm, BondId,
    BooleanForm, GraphSymmetryConfig, Lattice, LigandPermutation, Molecule, MoleculeEntries,
    Normalize, NumForm, OrientedLigandPermutation, Reframe, RingConfig, RingModel,
    StereoAtomConstraintForm, StereoAtomConstraintKey, StereoAtomId, StereoBondConstraintForm,
    StereoBondConstraintKey, StereoBondId, StereoLigand, StereoLigandKind, StereoLigandPair,
    StereoLigandPosition, StereogenicityForm, TopicityRelationForm,
};
use umol_perm::{Orientation, Permutation};

use super::strategies::{
    molecule_with_constraints_strategy, raw_num_form_strategy, raw_stereo_atom_form_strategy,
    raw_stereo_bond_form_strategy,
};

macro_rules! check_readings {
    ($before:expr, $after:expr, $stored:expr, $key:expr) => {{
        let before = $before;
        let after = $after;
        let key = $key;
        let stored = $stored;
        let expected = (stored.asserted(key), stored.asserted_complete(key));
        prop_assert_eq!(before.asserted(key), expected.0);
        prop_assert_eq!(before.asserted_complete(key), expected.1);
        let original = [
            before.asserted(key),
            before.asserted_complete(key),
            before.derived(key),
            before.derived_complete(key),
            before.evidenced(key),
            before.evidenced_complete(key),
        ];
        let normalized = [
            after.asserted(key),
            after.asserted_complete(key),
            after.derived(key),
            after.derived_complete(key),
            after.evidenced(key),
            after.evidenced_complete(key),
        ];
        prop_assert_eq!(&original, &normalized);
        let closed = original[4].clone().and_then(|value| {
            if value.is_undetermined() {
                original[3].clone()
            } else {
                Some(value)
            }
        });
        prop_assert_eq!(&original[5], &closed);
        for value in original.into_iter().flatten() {
            prop_assert_eq!(value.clone().normalize(), Ok(value));
        }
    }};
}

proptest! {
    #[test]
    fn prop_constraint_readings_normalization(molecule in molecule_with_constraints_strategy()) {
        if let Ok(normal) = molecule.clone().normalize() {
            let rings = molecule.rings(RingModel::default(), RingConfig::default()).into_ring_set();
            for atom in molecule.atoms().iter() {
                let after = normal.atom(atom.id());
                for key in [AtomConstraintKey::Valence, AtomConstraintKey::DonatedPairs, AtomConstraintKey::AcceptedPairs, AtomConstraintKey::AromaticValence, AtomConstraintKey::MulticenterValence, AtomConstraintKey::TetrahedralStereo, AtomConstraintKey::Degree, AtomConstraintKey::TotalDegree, AtomConstraintKey::TotalHydrogens, AtomConstraintKey::TotalValence, AtomConstraintKey::RingValence] {
                    check_readings!(atom.constraints().with_rings(&rings), after.constraints().with_rings(&rings), &atom.attributes().constraints, key);
                }
                prop_assert_eq!(atom.valence(), after.valence());
                prop_assert_eq!(atom.covalence(), after.covalence());
                prop_assert_eq!(atom.heavy_atom_degree(), after.heavy_atom_degree());
                prop_assert_eq!(atom.heavy_atom_valence(), after.heavy_atom_valence());
            }
            for bond in molecule.bonds().iter() {
                for entry in bond.constraints().iter() {
                    check_readings!(bond.constraints().with_rings(&rings), normal.bond(bond.id()).constraints().with_rings(&rings), &bond.attributes().constraints, entry.key());
                }
            }
            for entity in molecule.dative_bonds().iter() {
                for entry in entity.constraints().iter() { check_readings!(entity.constraints(), normal.dative_bond(entity.id()).constraints(), &entity.attributes().constraints, entry.key()); }
            }
            for entity in molecule.aromatic_systems().iter() {
                for entry in entity.constraints().iter() { check_readings!(entity.constraints(), normal.aromatic_system(entity.id()).constraints(), &entity.attributes().constraints, entry.key()); }
            }
            for entity in molecule.multicenter_bonds().iter() {
                for entry in entity.constraints().iter() { check_readings!(entity.constraints(), normal.multicenter_bond(entity.id()).constraints(), &entity.attributes().constraints, entry.key()); }
            }
            for entity in molecule.noncovalent_bonds().iter() {
                for entry in entity.constraints().iter() { check_readings!(entity.constraints(), normal.noncovalent_bond(entity.id()).constraints(), &entity.attributes().constraints, entry.key()); }
            }
            for entity in molecule.stereo_atoms().iter() {
                let before = entity.constraints();
                let after = normal.stereo_atom(entity.id()).constraints();
                for entry in before.iter() {
                    let key = entry.key();
                    let stored = &entity.attributes().constraints;
                    let expected = if entity.attributes().configuration.kind().is_none() {
                        (Some(StereoAtomConstraintForm::undetermined(key)), Some(StereoAtomConstraintForm::undetermined(key)))
                    } else {
                        (stored.asserted(key), stored.asserted_complete(key))
                    };
                    prop_assert_eq!(before.asserted(key), expected.0.clone());
                    prop_assert_eq!(before.asserted_complete(key), expected.1.clone());
                    prop_assert_eq!(after.asserted(key), expected.0);
                    prop_assert_eq!(after.asserted_complete(key), expected.1);
                    for value in [before.asserted(key), before.asserted_complete(key)].into_iter().flatten() {
                        prop_assert_eq!(value.clone().normalize(), Ok(value));
                    }
                }
            }
            for entity in molecule.stereo_bonds().iter() {
                let before = entity.constraints();
                let after = normal.stereo_bond(entity.id()).constraints();
                for entry in before.iter() {
                    let key = entry.key();
                    let stored = &entity.attributes().constraints;
                    let expected = if entity.attributes().configuration.kind().is_none() {
                        (Some(StereoBondConstraintForm::undetermined(key)), Some(StereoBondConstraintForm::undetermined(key)))
                    } else {
                        (stored.asserted(key), stored.asserted_complete(key))
                    };
                    prop_assert_eq!(before.asserted(key), expected.0.clone());
                    prop_assert_eq!(before.asserted_complete(key), expected.1.clone());
                    prop_assert_eq!(after.asserted(key), expected.0);
                    prop_assert_eq!(after.asserted_complete(key), expected.1);
                    for value in [before.asserted(key), before.asserted_complete(key)].into_iter().flatten() {
                        prop_assert_eq!(value.clone().normalize(), Ok(value));
                    }
                }
            }
        }
    }

    #[test]
    fn prop_stereo_constraints_view_readings(
        atom in raw_stereo_atom_form_strategy(),
        bond in raw_stereo_bond_form_strategy(),
        include_constraints in any::<bool>(),
    ) {
        let degree = atom.configuration.kind().map_or(4, |kind| kind.degree());
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); (degree + 1).max(6)],
            bonds: (1..=degree).map(|other| (AtomId(0), AtomId(other as u32), BondForm::from_order(1)))
                .chain([(AtomId(1), AtomId(4), BondForm::from_order(1)),
                        (AtomId(1), AtomId(5), BondForm::from_order(1))]).collect(),
            stereo_atoms: vec![(AtomId(0), (1..=degree)
                .map(|id| StereoLigand::new(AtomId(id as u32), StereoLigandKind::Atom)).collect(), atom)],
            stereo_bonds: vec![(BondId(0), [2, 3, 4, 5].into_iter()
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).collect(), bond)],
            ..Default::default()
        });
        let original = molecule.clone();
        let symmetry = molecule.graph_symmetry(&GraphSymmetryConfig {
            include_constraints,
            ..Default::default()
        }).into_graph_symmetry();
        let stereo = molecule.stereo_atom(StereoAtomId(0));
        let view = stereo.constraints().with_symmetry(&symmetry);
        let keys = view.iter().map(|entry| entry.key()).chain([
            StereoAtomConstraintKey::Stereogenicity,
            StereoAtomConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))),
            StereoAtomConstraintKey::LigandSymmetry(OrientedLigandPermutation {
                permutation: LigandPermutation(Permutation::identity(stereo.ligand_count())),
                orientation: Orientation::Proper,
            }),
        ]);
        for key in keys {
            let derived = view.derived(key);
            if stereo.attributes().configuration.kind().is_none() {
                let top = Some(StereoAtomConstraintForm::undetermined(key));
                prop_assert_eq!(&derived, &top);
                prop_assert_eq!(view.asserted(key), top);
            } else {
                match (key, &derived) {
                    (StereoAtomConstraintKey::Stereogenicity, Some(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(_)))) => {},
                    (StereoAtomConstraintKey::Topicity(pair), Some(StereoAtomConstraintForm::Topicity(value))) => {
                        prop_assert_eq!(value.pair, pair);
                        prop_assert!(matches!(value.relation, TopicityRelationForm::Lit(_)));
                    },
                    (StereoAtomConstraintKey::LigandSymmetry(permutation), Some(StereoAtomConstraintForm::LigandSymmetry(value))) => {
                        prop_assert_eq!(value.permutation, permutation);
                        prop_assert!(matches!(value.invariant, BooleanForm::Lit(_)));
                    },
                    _ => prop_assert!(false, "derived reading must be literal at its key"),
                }
            }
            prop_assert_eq!(view.derived_complete(key), derived.clone());
            let expected = view.asserted(key).and_then(|asserted| derived.as_ref().and_then(|value| asserted.meet(value)));
            prop_assert_eq!(view.evidenced(key), expected.clone());
            prop_assert_eq!(view.evidenced_complete(key), expected);
            for value in [derived, view.evidenced(key)].into_iter().flatten() {
                prop_assert_eq!(value.clone().normalize(), Ok(value));
            }
        }
        let stereo = molecule.stereo_bond(StereoBondId(0));
        let view = stereo.constraints().with_symmetry(&symmetry);
        let keys = view.iter().map(|entry| entry.key()).chain([
            StereoBondConstraintKey::Stereogenicity,
            StereoBondConstraintKey::Topicity(StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1))),
            StereoBondConstraintKey::LigandSymmetry(OrientedLigandPermutation {
                permutation: LigandPermutation(Permutation::identity(stereo.ligand_count())),
                orientation: Orientation::Proper,
            }),
        ]);
        for key in keys {
            let derived = view.derived(key);
            if stereo.attributes().configuration.kind().is_none() {
                let top = Some(StereoBondConstraintForm::undetermined(key));
                prop_assert_eq!(&derived, &top);
                prop_assert_eq!(view.asserted(key), top);
            } else {
                match (key, &derived) {
                    (StereoBondConstraintKey::Stereogenicity, Some(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(_)))) => {},
                    (StereoBondConstraintKey::Topicity(pair), Some(StereoBondConstraintForm::Topicity(value))) => {
                        prop_assert_eq!(value.pair, pair);
                        prop_assert!(matches!(value.relation, TopicityRelationForm::Lit(_)));
                    },
                    (StereoBondConstraintKey::LigandSymmetry(permutation), Some(StereoBondConstraintForm::LigandSymmetry(value))) => {
                        prop_assert_eq!(value.permutation, permutation);
                        prop_assert!(matches!(value.invariant, BooleanForm::Lit(_)));
                    },
                    _ => prop_assert!(false, "derived reading must be literal at its key"),
                }
            }
            prop_assert_eq!(view.derived_complete(key), derived.clone());
            let expected = view.asserted(key).and_then(|asserted| derived.as_ref().and_then(|value| asserted.meet(value)));
            prop_assert_eq!(view.evidenced(key), expected.clone());
            prop_assert_eq!(view.evidenced_complete(key), expected);
            for value in [derived, view.evidenced(key)].into_iter().flatten() {
                prop_assert_eq!(value.clone().normalize(), Ok(value));
            }
        }
        prop_assert_eq!(molecule, original);
    }

    #[test]
    fn prop_constraint_readings_reframe(molecule in molecule_with_constraints_strategy()) {
        if let Ok(reframed) = molecule.clone().reframe() {
            for atom in molecule.atoms().iter() {
                let before = atom.constraints();
                let after = reframed.atom(atom.id()).constraints();
                let key = AtomConstraintKey::TetrahedralStereo;
                prop_assert_eq!(before.derived(key), after.derived(key));
                prop_assert_eq!(before.derived_complete(key), after.derived_complete(key));
                prop_assert_eq!(before.evidenced(key), after.evidenced(key));
                prop_assert_eq!(before.evidenced_complete(key), after.evidenced_complete(key));
            }
            for bond in molecule.bonds().iter() {
                let before = bond.constraints();
                let after = reframed.bond(bond.id()).constraints();
                let key = BondConstraintKey::CisTransStereo;
                prop_assert_eq!(before.derived(key), after.derived(key));
                prop_assert_eq!(before.derived_complete(key), after.derived_complete(key));
                prop_assert_eq!(before.evidenced(key), after.evidenced(key));
                prop_assert_eq!(before.evidenced_complete(key), after.evidenced_complete(key));
            }
        }
    }

    #[test]
    fn prop_atom_constraints_view_evidenced(
        assertion in raw_num_form_strategy(), order in raw_num_form_strategy(),
    ) {
        let mut atom = AtomForm::from_element(Element::C);
        atom.constraints.set(AtomConstraintForm::Valence(assertion.clone()));
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![atom, AtomForm::from_element(Element::C)],
            bonds: vec![(AtomId(0), AtomId(1), BondForm { order: order.clone(), ..Default::default() })],
            ..Default::default()
        });
        let normalized_order = order.normalize().unwrap();
        let derived = match normalized_order { NumForm::Lit(n) => NumForm::Lit(n), _ => NumForm::Undetermined };
        let expected = assertion.meet(&derived).map(AtomConstraintForm::Valence);
        let reading = molecule.atom(AtomId(0)).constraints();
        prop_assert_eq!(reading.evidenced(AtomConstraintKey::Valence), expected.clone());
        prop_assert_eq!(reading.evidenced_complete(AtomConstraintKey::Valence), expected);
    }
}
