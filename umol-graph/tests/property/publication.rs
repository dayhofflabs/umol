//! Representation-integrity preservation by chemistry-layer publishers.
//!
//! Resolved SMILES fixtures also check equality between consuming and borrowed
//! transformations; publication is checked through the editor's integrity gate.
//! Constructed single and coupled C/N/O/S rings with added valence assertions and
//! altered charges check Kekulizer's planned-atom conservation contract through
//! output validation and derived constraint readings. The rings have fixed donor
//! positions and admit matching; every failure must preserve the borrowed input.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_chem::element::Element;
use umol_graph::ingest::ingest_smiles;
use umol_graph::ops::model::{AromaticityModel, ChemistryModel, ValenceModel};
use umol_graph::ops::resolve::Resolver;
use umol_graph::ops::transform::{
    Aromatizer, DelocalizeCharge, KekulizeConfig, KekulizeError, Kekulizer,
    MaximumMatchingAlgorithm, Transformer,
};
use umol_graph::ops::validate::ValenceInvariantsValidator;
use umol_graph_ir::ir::{
    AromaticSystemForm, AtomConstraintForm, AtomConstraintKey, AtomForm, AtomId, BondForm, Lattice,
    Molecule, MoleculeEntries, NumForm,
};
use umol_utils::solution::Solution;

fn molecule_from_smiles(source: &str) -> Molecule {
    ingest_smiles(source).expect("publication-property SMILES fixture resolves")
}

fn kekule_molecule_strategy() -> impl Strategy<Value = Molecule> {
    prop::sample::select(vec!["C1=CC=CC=C1", "C1=CC=CN=C1"]).prop_map(molecule_from_smiles)
}

fn aromatic_molecule_strategy() -> impl Strategy<Value = Molecule> {
    prop::sample::select(vec!["c1ccccc1", "c1ccncc1"]).prop_map(molecule_from_smiles)
}

fn resolved_molecule_strategy() -> impl Strategy<Value = Molecule> {
    prop::sample::select(vec!["C", "CC", "C=C", "N", "O", "c1ccccc1"])
        .prop_map(molecule_from_smiles)
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_aromatizer_integrity_preservation(source in kekule_molecule_strategy()) {
        let transformer = Aromatizer::new(&AromaticityModel::daylight());
        let mut borrowed = source.clone();
        let published = transformer
            .transform(source)
            .map_err(|error| TestCaseError::fail(format!("aromatization failed: {error}")))?;
        transformer.transform_into(&mut borrowed)
            .map_err(|error| TestCaseError::fail(format!("aromatization failed: {error}")))?;

        prop_assert_eq!(&borrowed, &published);
        prop_assert_eq!(published.clone().edit().finish(), Ok(published));
    }

    #[test]
    fn test_delocalize_charge_integrity_preservation(source in aromatic_molecule_strategy()) {
        let mut borrowed = source.clone();
        let published = DelocalizeCharge
            .transform(source)
            .expect("delocalized-charge transformation is infallible");
        DelocalizeCharge.transform_into(&mut borrowed)
            .expect("delocalized-charge transformation is infallible");

        prop_assert_eq!(&borrowed, &published);
        prop_assert_eq!(published.clone().edit().finish(), Ok(published));
    }

    #[test]
    fn test_kekulizer_integrity_preservation(
        source in aromatic_molecule_strategy(),
        use_bipartite in any::<bool>(),
    ) {
        let algorithm = if use_bipartite {
            MaximumMatchingAlgorithm::HopcroftKarp
        } else {
            MaximumMatchingAlgorithm::Edmonds
        };
        let node_order = source.atoms().ids().collect::<Vec<AtomId>>();
        let transformer = Kekulizer::new(KekulizeConfig::new(algorithm), node_order);
        let mut borrowed = source.clone();
        let published = transformer
            .transform(source)
            .map_err(|error| TestCaseError::fail(format!("kekulization failed: {error}")))?;
        transformer.transform_into(&mut borrowed)
            .map_err(|error| TestCaseError::fail(format!("kekulization failed: {error}")))?;

        prop_assert_eq!(&borrowed, &published);
        prop_assert_eq!(published.clone().edit().finish(), Ok(published));
    }

    #[test]
    fn test_kekulizer_transform_valence(
        rings in prop::collection::vec(
            prop::sample::select(vec![Element::C, Element::N, Element::O, Element::S]), 1..=3),
        selected in any::<prop::sample::Index>(),
        valence in prop::option::of(0i64..=4),
        charge in -1i64..=1,
    ) {
        let mut entries = MoleculeEntries::default();
        let mut previous_bridge: Option<usize> = None;
        for element in rings {
            let start = entries.atoms.len();
            let size = if element == Element::C { 6 } else { 5 };
            for position in 0..size {
                let element = if position == 0 { element } else { Element::C };
                entries.atoms.push(AtomForm::from_element(element)
                    .with_implicit_hydrogens(i64::from(matches!(element, Element::C | Element::N)))
                    .with_lone_pairs(i64::from(matches!(element, Element::O | Element::S)))
                    .into_concrete());
                entries.bonds.push((AtomId((start + position) as u32),
                    AtomId((start + (position + 1) % size) as u32),
                    BondForm::from_order(1).into_concrete()));
            }
            let mut electrons = vec![1; size];
            if element != Element::C {
                electrons[0] = 2;
            }
            entries.aromatic.push(((start..start + size).map(|id| AtomId(id as u32)).collect(),
                AromaticSystemForm::from_electrons(electrons).into_concrete()));
            if let Some(previous) = previous_bridge {
                let current = start + 1;
                entries.bonds.push((AtomId(previous as u32), AtomId(current as u32),
                    BondForm::from_order(1).into_concrete()));
                entries.atoms[previous].implicit_hydrogens = NumForm::Lit(0);
                entries.atoms[current].implicit_hydrogens = NumForm::Lit(0);
            }
            previous_bridge = Some(start + size - 1);
        }
        let mut source = Molecule::from_entries(entries);
        prop_assert_eq!(ValenceInvariantsValidator.validate(&source), Ok(Solution::Determined(())));
        let members: Vec<_> = source.atoms().ids().collect();
        let id = members[selected.index(members.len())];
        {
            let mut atom = source.atom_mut(id);
            atom.attributes_mut().charge = NumForm::Lit(charge);
            if let Some(valence) = valence {
                atom.attributes_mut().constraints.set(AtomConstraintForm::valence(valence));
            }
        }
        let transformer = Kekulizer::new(KekulizeConfig::default(), source.atoms().ids().collect());
        let output = transformer.transform(source.clone());
        let mut borrowed = source.clone();
        match output {
            Ok(output) => {
                prop_assert_eq!(transformer.transform_into(&mut borrowed), Ok(()));
                prop_assert_eq!(&borrowed, &output);
                prop_assert_eq!(ValenceInvariantsValidator.validate(&output), Ok(Solution::Determined(())));
                for id in members {
                    let before = source.atom(id);
                    let after = output.atom(id);
                    prop_assert_eq!(after.charge(), before.charge());
                    let mut retained = before.attributes().constraints.clone();
                    retained.remove(AtomConstraintKey::AromaticValence);
                    prop_assert_eq!(&after.attributes().constraints, &retained);
                    let derived = after.constraints().derived_complete(AtomConstraintKey::Valence).unwrap();
                    let asserted = after.constraints().asserted(AtomConstraintKey::Valence).unwrap();
                    prop_assert!(derived.satisfies(&asserted));
                }
            }
            Err(error @ (KekulizeError::ContradictedAssertion { .. }
                | KekulizeError::ValenceInvariant { .. })) => {
                prop_assert_eq!(transformer.transform_into(&mut borrowed), Err(error));
                prop_assert_eq!(borrowed, source);
            }
            Err(error) => return Err(TestCaseError::fail(format!("unexpected planning error: {error}"))),
        }
    }

    #[test]
    fn test_resolver_integrity_preservation(mut source in resolved_molecule_strategy()) {
        let model = ChemistryModel {
            valence: ValenceModel::smiles(),
            ..ChemistryModel::default()
        };
        let resolver = Resolver::new(&model);
        let consumed = resolver.resolve(source.clone())
            .map_err(|error| TestCaseError::fail(format!("resolution failed: {error}")))?;
        let Solution::Determined(published) = consumed else {
            return Err(TestCaseError::fail("resolved input must remain determined"));
        };
        prop_assert_eq!(published.clone().edit().finish(), Ok(published));
        resolver
            .resolve_into(&mut source)
            .map_err(|error| TestCaseError::fail(format!("resolution failed: {error}")))?;

        prop_assert_eq!(source.clone().edit().finish(), Ok(source));
    }
}
