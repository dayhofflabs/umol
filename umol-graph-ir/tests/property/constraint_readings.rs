//! Reading normality, normalization invariance, and meet composition on raw forms.
//! Molecule cases compare public readings across successful whole normalization;
//! numeric cases independently normalize sources before computing their meet.
//! Container/view equality checks own assertions except for kindless stereo entities,
//! whose assertions are independently required to be top despite their raw payload;
//! container unit cases independently cover raw contradictions and closure of absence.
//! Reframing invariance checks stereo readings on generated molecules whose stored
//! participant orders change, independently of the fixed-frame normalization law.

use proptest::prelude::*;
use umol_chem::element::Element;
use umol_graph_ir::ir::{
    AtomConstraintForm, AtomConstraintKey, AtomForm, AtomId, BondConstraintKey, BondForm, BondId,
    Lattice, Molecule, MoleculeEntries, Normalize, NumForm, Reframe, RingConfig, RingModel,
    StereoAtomConstraintForm, StereoAtomId, StereoBondConstraintForm, StereoBondId, StereoLigand,
    StereoLigandKind,
};

use super::strategies::{
    molecule_with_constraints_strategy, raw_num_form_strategy, raw_stereo_atom_form_strategy,
    raw_stereo_bond_form_strategy,
};

macro_rules! check_readings {
    ($before:expr, $after:expr, $stored:expr, $key:expr $(, $kindless:expr, $form:ty)?) => {{
        let before = $before;
        let after = $after;
        let key = $key;
        let stored = $stored;
        let expected = (stored.asserted(key), stored.asserted_complete(key));
        $(let expected = if $kindless {
            (Some(<$form>::undetermined(key)), Some(<$form>::undetermined(key)))
        } else { expected };)?
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
                for entry in entity.constraints().iter() { check_readings!(entity.constraints(), normal.stereo_atom(entity.id()).constraints(), &entity.attributes().constraints, entry.key(), entity.attributes().configuration.kind().is_none(), StereoAtomConstraintForm); }
            }
            for entity in molecule.stereo_bonds().iter() {
                for entry in entity.constraints().iter() { check_readings!(entity.constraints(), normal.stereo_bond(entity.id()).constraints(), &entity.attributes().constraints, entry.key(), entity.attributes().configuration.kind().is_none(), StereoBondConstraintForm); }
            }
        }
    }

    #[test]
    fn prop_stereo_constraints_view_normalization(
        atom in raw_stereo_atom_form_strategy(),
        bond in raw_stereo_bond_form_strategy(),
    ) {
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C); 6],
            bonds: [(0, 1), (0, 2), (0, 3), (0, 4), (1, 4), (1, 5)]
                .into_iter().map(|(a, b)| (AtomId(a), AtomId(b), BondForm::from_order(1))).collect(),
            stereo_atoms: vec![(AtomId(0), (1..=4)
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).collect(), atom)],
            stereo_bonds: vec![(BondId(0), [2, 3, 4, 5].into_iter()
                .map(|id| StereoLigand::new(AtomId(id), StereoLigandKind::Atom)).collect(), bond)],
            ..Default::default()
        });
        if let Ok(normal) = molecule.clone().normalize() {
            let atom = molecule.stereo_atom(StereoAtomId(0));
            let bond = molecule.stereo_bond(StereoBondId(0));
            for entry in atom.constraints().iter() {
                check_readings!(atom.constraints(), normal.stereo_atom(atom.id()).constraints(),
                    &atom.attributes().constraints, entry.key(),
                    atom.attributes().configuration.kind().is_none(), StereoAtomConstraintForm);
            }
            for entry in bond.constraints().iter() {
                check_readings!(bond.constraints(), normal.stereo_bond(bond.id()).constraints(),
                    &bond.attributes().constraints, entry.key(),
                    bond.attributes().configuration.kind().is_none(), StereoBondConstraintForm);
            }
        }
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
