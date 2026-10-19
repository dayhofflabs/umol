//! Aromaticity transformer properties over generated ring representations.
//!
//! Raw C/N/O/S skeletons retain contradictions and undetermined fields. Stored
//! systems and saturated carbon skeletons provide independent successful inputs.
//! Publication checks compare consuming and in-place outcomes, representation
//! integrity, and unchanged inputs on rejection. Successful Kekule outputs also
//! exercise aromatization of localized double bonds.
//!
//! Kekulizer valence checks use matching-admitting rings with fixed donor positions.
//! Aromatizer additions must agree with open atom evidence and exclude non-aromatic
//! bonds. Independent carbon skeletons with disjoint localized double bonds and
//! hydrogens completing valence four exercise all three aromaticity rules, both
//! unasserted and carrying one varied atom and bond assertion.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_chem::element::Element;
use umol_graph::ops::model::{AromaticityModel, AromaticityRule, RingLimits};
use umol_graph::ops::transform::{
    Aromatizer, DelocalizeCharge, KekulizeConfig, KekulizeError, Kekulizer,
    MaximumMatchingAlgorithm, Transformer,
};
use umol_graph::ops::validate::ValenceInvariantsValidator;
use umol_graph_ir::ir::{
    AromaticSystemForm, AromaticValenceForm, AtomConstraintForm, AtomConstraintKey, AtomForm,
    AtomId, BondConstraintForm, BondForm, BondId, BooleanForm, ElectronCountsForm, Lattice,
    Molecule, MoleculeEntries, NumForm,
};
use umol_utils::solution::Solution;

use crate::strategies::{aromatic_skeleton, publication_molecule};

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            crate::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_aromatizer_transform_constraints(
        skeleton in aromatic_skeleton(),
        atom_index in any::<prop::sample::Index>(),
        bond_index in any::<prop::sample::Index>(),
        rule in 0usize..3,
    ) {
        let source = skeleton.molecule();
        let mut has_double = vec![false; source.atoms().count()];
        let bonds = source.bonds().iter().map(|bond| {
            let [a, b] = bond.atom_ids();
            let order = if !has_double[a.index()] && !has_double[b.index()] {
                has_double[a.index()] = true;
                has_double[b.index()] = true;
                2
            } else {
                1
            };
            (a, b, BondForm::from_order(order).into_concrete())
        }).collect();
        let unasserted = Molecule::from_entries(MoleculeEntries {
            atoms: source.atoms().iter().map(|atom| {
                AtomForm::from_element(Element::C).with_implicit_hydrogens(
                    4 - atom.neighbors().count() as i64 - i64::from(has_double[atom.id().index()])
                ).into_concrete()
            }).collect(),
            bonds,
            ..Default::default()
        });
        let mut asserted = unasserted.clone();
        let atom = AtomId(atom_index.index(source.atoms().count()) as u32);
        let bond = BondId(bond_index.index(source.bonds().count()) as u32);
        asserted.atom_mut(atom).attributes_mut().constraints.set(
            AtomConstraintForm::aromatic_valence(source.atom(atom).constraints().aromatic_valence().clone())
        );
        asserted.bond_mut(bond).attributes_mut().constraints.set(
            BondConstraintForm::Aromatic(*source.bond(bond).constraints().aromatic())
        );
        let transformer = Aromatizer::new(&AromaticityModel {
            rule: match rule {
                0 => AromaticityRule::Hueckel { ring_limits: RingLimits::default() },
                1 => AromaticityRule::Hmo { stabilization_threshold: 0.1 },
                _ => AromaticityRule::Clar,
            },
            ..AromaticityModel::daylight()
        });
        for input in [unasserted, asserted] {
            if let Ok(output) = transformer.transform(input.clone()) {
                for system in output.aromatic_systems().iter() {
                    let ElectronCountsForm::Lit(electrons) = &system.attributes().electrons else {
                        return Err(TestCaseError::fail("new system has nonliteral electron counts"));
                    };
                    let members: Vec<_> = system.atom_ids().collect();
                    prop_assert_eq!(electrons.len(), members.len());
                    for (&atom, &electrons) in members.iter().zip(electrons) {
                        let literal = AtomConstraintForm::aromatic_valence(
                            AromaticValenceForm::aromatic(electrons)
                        );
                        let evidence = input.atom(atom).constraints().evidenced(AtomConstraintKey::AromaticValence);
                        prop_assert_eq!(evidence.and_then(|value| value.meet(&literal)), Some(literal));
                    }
                    for bond in input.bonds().iter() {
                        if matches!(bond.constraints().aromatic(), BooleanForm::Lit(false)) {
                            prop_assert!(!bond.atom_ids().iter().all(|atom| members.contains(atom)));
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_transformer_transform_publication(
        source in publication_molecule(),
        use_bipartite in any::<bool>(),
    ) {
        let mut borrowed = source.clone();
        let published = DelocalizeCharge.transform(source.clone()).unwrap();
        prop_assert_eq!(DelocalizeCharge.transform_into(&mut borrowed), Ok(()));
        prop_assert_eq!(&borrowed, &published);
        prop_assert_eq!(published.clone().edit().finish(), Ok(published));

        let algorithm = if use_bipartite {
            MaximumMatchingAlgorithm::HopcroftKarp
        } else {
            MaximumMatchingAlgorithm::Edmonds
        };
        let transformer = Kekulizer::new(KekulizeConfig::new(algorithm), source.atoms().ids().collect());
        let mut borrowed = source.clone();
        let consumed = transformer.transform(source.clone());
        let in_place = transformer.transform_into(&mut borrowed);
        let mut kekule = vec![source.clone()];
        match consumed {
            Ok(published) => {
                prop_assert_eq!(in_place, Ok(()));
                prop_assert_eq!(&borrowed, &published);
                prop_assert_eq!(published.clone().edit().finish(), Ok(published.clone()));
                kekule.push(published);
            }
            Err(error) => {
                prop_assert_eq!(in_place, Err(error));
                prop_assert_eq!(borrowed, source);
            }
        }
        let transformer = Aromatizer::new(&AromaticityModel::daylight());
        for source in kekule {
            let mut borrowed = source.clone();
            let consumed = transformer.transform(source.clone());
            let in_place = transformer.transform_into(&mut borrowed);
            match consumed {
                Ok(published) => {
                    prop_assert_eq!(in_place, Ok(()));
                    prop_assert_eq!(&borrowed, &published);
                    prop_assert_eq!(published.clone().edit().finish(), Ok(published));
                }
                Err(error) => {
                    prop_assert_eq!(in_place, Err(error));
                    prop_assert_eq!(borrowed, source);
                }
            }
        }
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
}
