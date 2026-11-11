//! Participant-frame action properties.
//!
//! Local entity forms use independently generated compatible permutations. A raw aggregate
//! scenario with all six overlay kinds supplies two nonidentity, degree-compatible actions derived
//! from independent atom relabelings. The same actions exercise entity-kind span aggregates and the
//! three root aggregate shapes. Missing domains, wrong degrees, and covering supersets are checked
//! separately so successful action laws do not hide compatibility failures.
//!
//! The reframe-kind properties check `reframe_kind` against `reframe_by` on the local forms. A field
//! or a constraint whose kind is `Independent` is equal before and after a compatible action, and a
//! `Positional` field is reordered as the participants are. The stereo forms have no field of either
//! kind, so only their constraints are covered.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use strum::VariantArray;
use umol_graph_ir::ir::{
    AromaticSystemFieldKind, AromaticSystemForm, DativeBondFieldKind, DativeBondId,
    ElectronCountsForm, FrameTransport, MulticenterBondFieldKind, MulticenterBondForm,
    MulticenterBondId, NoncovalentBondFieldKind, Reaction, Reframe, ReframeKind, StereoKind,
};
use umol_perm::{DynPermutation, Permutation};

use crate::strategies::{
    aromatic_system_form_for, dative_bond_strategy, multicenter_bond_form_for,
    noncovalent_bond_form_strategy, standardization_scenario_strategy, stereo_atom_form_strategy,
    stereo_bond_form_strategy, stereo_frame_permutation_strategy,
};

fn dynamic_action_pair_strategy(
    degree: usize,
) -> impl Strategy<Value = (DynPermutation, DynPermutation)> {
    (
        Just((0..degree).collect::<Vec<_>>()).prop_shuffle(),
        Just((0..degree).collect::<Vec<_>>()).prop_shuffle(),
    )
        .prop_map(|(first, second)| {
            (
                DynPermutation::try_from(first).expect("a shuffled image is a permutation"),
                DynPermutation::try_from(second).expect("a shuffled image is a permutation"),
            )
        })
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            concat!(env!("CARGO_MANIFEST_DIR"), "/tests/property/frame.proptest-regressions"),
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_dative_bond_form_reframe_by(
        (form, first, second) in (1usize..=6).prop_flat_map(|degree| (
            dative_bond_strategy(),
            dynamic_action_pair_strategy(degree),
        )).prop_map(|(form, (first, second))| (form, first, second)),
    ) {
        let identity = DynPermutation::identity(first.degree());
        let inverse = first.inverse();
        let composite = first.compose(&second).expect("the actions have the same degree");

        prop_assert_eq!(form.clone().reframe_by(&identity), Some(form.clone()));
        prop_assert_eq!(
            form.clone().reframe_by(&first).and_then(|value| value.reframe_by(&inverse)),
            Some(form.clone()),
        );
        prop_assert_eq!(
            form.clone().reframe_by(&first).and_then(|value| value.reframe_by(&second)),
            form.reframe_by(&composite),
        );
    }

    #[test]
    fn test_aromatic_system_form_reframe_by(
        (form, first, second) in (3usize..=6).prop_flat_map(|degree| (
            aromatic_system_form_for(degree),
            dynamic_action_pair_strategy(degree),
        )).prop_map(|(form, (first, second))| (form, first, second)),
    ) {
        let identity = DynPermutation::identity(first.degree());
        let inverse = first.inverse();
        let composite = first.compose(&second).expect("the actions have the same degree");

        prop_assert_eq!(form.clone().reframe_by(&identity), Some(form.clone()));
        prop_assert_eq!(
            form.clone().reframe_by(&first).and_then(|value| value.reframe_by(&inverse)),
            Some(form.clone()),
        );
        prop_assert_eq!(
            form.clone().reframe_by(&first).and_then(|value| value.reframe_by(&second)),
            form.clone().reframe_by(&composite),
        );
        prop_assert_eq!(
            form.reframe_by(&DynPermutation::identity(first.degree() + 1)),
            None,
        );
    }

    #[test]
    fn test_multicenter_bond_form_reframe_by(
        (form, first, second) in (3usize..=6).prop_flat_map(|degree| (
            multicenter_bond_form_for(degree),
            dynamic_action_pair_strategy(degree),
        )).prop_map(|(form, (first, second))| (form, first, second)),
    ) {
        let identity = DynPermutation::identity(first.degree());
        let inverse = first.inverse();
        let composite = first.compose(&second).expect("the actions have the same degree");

        prop_assert_eq!(form.clone().reframe_by(&identity), Some(form.clone()));
        prop_assert_eq!(
            form.clone().reframe_by(&first).and_then(|value| value.reframe_by(&inverse)),
            Some(form.clone()),
        );
        prop_assert_eq!(
            form.clone().reframe_by(&first).and_then(|value| value.reframe_by(&second)),
            form.clone().reframe_by(&composite),
        );
        prop_assert_eq!(
            form.reframe_by(&DynPermutation::identity(first.degree() + 1)),
            None,
        );
    }

    #[test]
    fn test_noncovalent_bond_form_reframe_by(
        (form, first, second) in (
            noncovalent_bond_form_strategy(),
            dynamic_action_pair_strategy(2),
        ).prop_map(|(form, (first, second))| (form, first, second)),
    ) {
        let identity = DynPermutation::identity(first.degree());
        let inverse = first.inverse();
        let composite = first.compose(&second).expect("the actions have the same degree");

        prop_assert_eq!(form.clone().reframe_by(&identity), Some(form.clone()));
        prop_assert_eq!(
            form.clone().reframe_by(&first).and_then(|value| value.reframe_by(&inverse)),
            Some(form.clone()),
        );
        prop_assert_eq!(
            form.clone().reframe_by(&first).and_then(|value| value.reframe_by(&second)),
            form.clone().reframe_by(&composite),
        );
        prop_assert_eq!(form.reframe_by(&DynPermutation::identity(3)), None);
    }

    #[test]
    fn test_stereo_atom_form_reframe_by(
        (form, first, second) in stereo_atom_form_strategy().prop_flat_map(|form| {
            let kind = form.configuration.kind().expect("the strategy emits a kinded form");
            (
                Just(form),
                stereo_frame_permutation_strategy(kind),
                stereo_frame_permutation_strategy(kind),
            )
        }),
    ) {
        let identity = Permutation::identity(first.degree());
        let inverse = first.inverse();
        let composite = first.compose(second);

        prop_assert_eq!(form.clone().reframe_by(&identity), Some(form.clone()));
        prop_assert_eq!(
            form.clone().reframe_by(&first).and_then(|value| value.reframe_by(&inverse)),
            Some(form.clone()),
        );
        prop_assert_eq!(
            form.clone().reframe_by(&first).and_then(|value| value.reframe_by(&second)),
            form.clone().reframe_by(&composite),
        );
        prop_assert_eq!(
            form.reframe_by(&Permutation::identity(first.degree() - 1)),
            None,
        );
    }

    #[test]
    fn test_stereo_bond_form_reframe_by(
        (form, first, second) in (
            stereo_bond_form_strategy(),
            stereo_frame_permutation_strategy(StereoKind::CisTrans),
            stereo_frame_permutation_strategy(StereoKind::CisTrans),
        ),
    ) {
        let identity = Permutation::identity(first.degree());
        let inverse = first.inverse();
        let composite = first.compose(second);

        prop_assert_eq!(form.clone().reframe_by(&identity), Some(form.clone()));
        prop_assert_eq!(
            form.clone().reframe_by(&first).and_then(|value| value.reframe_by(&inverse)),
            Some(form.clone()),
        );
        prop_assert_eq!(
            form.clone().reframe_by(&first).and_then(|value| value.reframe_by(&second)),
            form.clone().reframe_by(&composite),
        );
        prop_assert_eq!(
            form.reframe_by(&Permutation::from_image(&[1, 2, 0, 3])),
            None,
        );
    }

    #[test]
    fn test_dative_bond_field_kind_reframe_kind(
        (form, action) in (1usize..=6).prop_flat_map(|degree| (
            dative_bond_strategy(),
            dynamic_action_pair_strategy(degree),
        )).prop_map(|(form, (action, _))| (form, action)),
    ) {
        let reframed = form.clone().reframe_by(&action).expect("a dative bond form takes any action");

        for &kind in DativeBondFieldKind::VARIANTS {
            let unchanged = match kind {
                DativeBondFieldKind::Order => reframed.order == form.order,
            };
            if kind.reframe_kind() == ReframeKind::Independent {
                prop_assert!(unchanged, "{:?}", kind);
            }
        }
    }

    #[test]
    fn test_dative_bond_constraint_key_reframe_kind(
        (form, action) in (1usize..=6).prop_flat_map(|degree| (
            dative_bond_strategy(),
            dynamic_action_pair_strategy(degree),
        )).prop_map(|(form, (action, _))| (form, action)),
    ) {
        let reframed = form.clone().reframe_by(&action).expect("a dative bond form takes any action");

        for constraint in form.constraints.iter() {
            let key = constraint.key();
            if key.reframe_kind() == ReframeKind::Independent {
                prop_assert_eq!(reframed.constraints.get(key), Some(constraint));
            }
        }
    }

    #[test]
    fn test_aromatic_system_field_kind_reframe_kind(
        (form, action) in (3usize..=6).prop_flat_map(|degree| (
            aromatic_system_form_for(degree),
            dynamic_action_pair_strategy(degree),
        )).prop_map(|(form, (action, _))| (form, action)),
    ) {
        let reframed = form.clone().reframe_by(&action).expect("the action has the degree of the form");
        let ElectronCountsForm::Lit(counts) = &form.electrons else {
            unreachable!("the strategy emits literal electron counts")
        };
        let permuted = AromaticSystemForm {
            electrons: ElectronCountsForm::Lit(
                action.act(counts).expect("the action has the degree of the form"),
            ),
            ..form.clone()
        };
        let same_field = |left: &AromaticSystemForm, right: &AromaticSystemForm, kind| match kind {
            AromaticSystemFieldKind::Electrons => left.electrons == right.electrons,
            AromaticSystemFieldKind::Charge => left.charge == right.charge,
            AromaticSystemFieldKind::UnpairedElectrons => {
                left.unpaired_electrons == right.unpaired_electrons
            }
        };

        for &kind in AromaticSystemFieldKind::VARIANTS {
            match kind.reframe_kind() {
                ReframeKind::Independent => {
                    prop_assert!(same_field(&reframed, &form, kind), "{:?}", kind);
                }
                ReframeKind::Positional => {
                    prop_assert!(same_field(&reframed, &permuted, kind), "{:?}", kind);
                }
                ReframeKind::Relative => {}
            }
        }
    }

    #[test]
    fn test_aromatic_system_constraint_key_reframe_kind(
        (form, action) in (3usize..=6).prop_flat_map(|degree| (
            aromatic_system_form_for(degree),
            dynamic_action_pair_strategy(degree),
        )).prop_map(|(form, (action, _))| (form, action)),
    ) {
        let reframed = form.clone().reframe_by(&action).expect("the action has the degree of the form");

        for constraint in form.constraints.iter() {
            let key = constraint.key();
            if key.reframe_kind() == ReframeKind::Independent {
                prop_assert_eq!(reframed.constraints.get(key), Some(constraint));
            }
        }
    }

    #[test]
    fn test_multicenter_bond_field_kind_reframe_kind(
        (form, action) in (3usize..=6).prop_flat_map(|degree| (
            multicenter_bond_form_for(degree),
            dynamic_action_pair_strategy(degree),
        )).prop_map(|(form, (action, _))| (form, action)),
    ) {
        let reframed = form.clone().reframe_by(&action).expect("the action has the degree of the form");
        let ElectronCountsForm::Lit(counts) = &form.electrons else {
            unreachable!("the strategy emits literal electron counts")
        };
        let permuted = MulticenterBondForm {
            electrons: ElectronCountsForm::Lit(
                action.act(counts).expect("the action has the degree of the form"),
            ),
            ..form.clone()
        };
        let same_field = |left: &MulticenterBondForm, right: &MulticenterBondForm, kind| match kind {
            MulticenterBondFieldKind::Electrons => left.electrons == right.electrons,
            MulticenterBondFieldKind::Charge => left.charge == right.charge,
            MulticenterBondFieldKind::UnpairedElectrons => {
                left.unpaired_electrons == right.unpaired_electrons
            }
        };

        for &kind in MulticenterBondFieldKind::VARIANTS {
            match kind.reframe_kind() {
                ReframeKind::Independent => {
                    prop_assert!(same_field(&reframed, &form, kind), "{:?}", kind);
                }
                ReframeKind::Positional => {
                    prop_assert!(same_field(&reframed, &permuted, kind), "{:?}", kind);
                }
                ReframeKind::Relative => {}
            }
        }
    }

    #[test]
    fn test_multicenter_bond_constraint_key_reframe_kind(
        (form, action) in (3usize..=6).prop_flat_map(|degree| (
            multicenter_bond_form_for(degree),
            dynamic_action_pair_strategy(degree),
        )).prop_map(|(form, (action, _))| (form, action)),
    ) {
        let reframed = form.clone().reframe_by(&action).expect("the action has the degree of the form");

        for constraint in form.constraints.iter() {
            let key = constraint.key();
            if key.reframe_kind() == ReframeKind::Independent {
                prop_assert_eq!(reframed.constraints.get(key), Some(constraint));
            }
        }
    }

    #[test]
    fn test_noncovalent_bond_field_kind_reframe_kind(
        (form, action) in (
            noncovalent_bond_form_strategy(),
            dynamic_action_pair_strategy(2),
        ).prop_map(|(form, (action, _))| (form, action)),
    ) {
        let reframed = form.clone().reframe_by(&action).expect("the action has degree two");

        for &kind in NoncovalentBondFieldKind::VARIANTS {
            let unchanged = match kind {
                NoncovalentBondFieldKind::Kind => reframed.kind == form.kind,
            };
            if kind.reframe_kind() == ReframeKind::Independent {
                prop_assert!(unchanged, "{:?}", kind);
            }
        }
    }

    #[test]
    fn test_noncovalent_bond_constraint_key_reframe_kind(
        (form, action) in (
            noncovalent_bond_form_strategy(),
            dynamic_action_pair_strategy(2),
        ).prop_map(|(form, (action, _))| (form, action)),
    ) {
        let reframed = form.clone().reframe_by(&action).expect("the action has degree two");

        for constraint in form.constraints.iter() {
            let key = constraint.key();
            if key.reframe_kind() == ReframeKind::Independent {
                prop_assert_eq!(reframed.constraints.get(key), Some(constraint));
            }
        }
    }

    #[test]
    fn test_stereo_atom_constraint_key_reframe_kind(
        (form, action) in stereo_atom_form_strategy().prop_flat_map(|form| {
            let kind = form.configuration.kind().expect("the strategy emits a kinded form");
            (Just(form), stereo_frame_permutation_strategy(kind))
        }),
    ) {
        let reframed = form.clone().reframe_by(&action).expect("the action belongs to the kind of the form");

        for constraint in form.constraints.iter() {
            let key = constraint.key();
            if key.reframe_kind() == ReframeKind::Independent {
                prop_assert_eq!(reframed.constraints.get(key), Some(constraint));
            }
        }
    }

    #[test]
    fn test_stereo_bond_constraint_key_reframe_kind(
        (form, action) in (
            stereo_bond_form_strategy(),
            stereo_frame_permutation_strategy(StereoKind::CisTrans),
        ),
    ) {
        let reframed = form.clone().reframe_by(&action).expect("the action belongs to the kind of the form");

        for constraint in form.constraints.iter() {
            let key = constraint.key();
            if key.reframe_kind() == ReframeKind::Independent {
                prop_assert_eq!(reframed.constraints.get(key), Some(constraint));
            }
        }
    }

    #[test]
    fn test_dative_bond_spans_reframe_by(scenario in standardization_scenario_strategy()) {
        let source = scenario.span.dative_bonds().clone();
        let first = scenario.first_action.dative_bonds();
        let second = scenario.second_action.dative_bonds();
        let inverse = first.inverse();
        let composite = first.compose(second).expect("the actions share one domain");

        prop_assert_eq!(source.clone().reframe_by(&first.identity()), Some(source.clone()));
        prop_assert_eq!(
            source.clone().reframe_by(first).and_then(|value| value.reframe_by(&inverse)),
            Some(source.clone()),
        );
        prop_assert_eq!(
            source.clone().reframe_by(first).and_then(|value| value.reframe_by(second)),
            source.reframe_by(&composite),
        );
    }

    #[test]
    fn test_aromatic_system_spans_reframe_by(scenario in standardization_scenario_strategy()) {
        let source = scenario.span.aromatic_systems().clone();
        let first = scenario.first_action.aromatic_systems();
        let second = scenario.second_action.aromatic_systems();
        let inverse = first.inverse();
        let composite = first.compose(second).expect("the actions share one domain");

        prop_assert_eq!(source.clone().reframe_by(&first.identity()), Some(source.clone()));
        prop_assert_eq!(
            source.clone().reframe_by(first).and_then(|value| value.reframe_by(&inverse)),
            Some(source.clone()),
        );
        prop_assert_eq!(
            source.clone().reframe_by(first).and_then(|value| value.reframe_by(second)),
            source.reframe_by(&composite),
        );
    }

    #[test]
    fn test_multicenter_bond_spans_reframe_by(scenario in standardization_scenario_strategy()) {
        let source = scenario.span.multicenter_bonds().clone();
        let first = scenario.first_action.multicenter_bonds();
        let second = scenario.second_action.multicenter_bonds();
        let inverse = first.inverse();
        let composite = first.compose(second).expect("the actions share one domain");

        prop_assert_eq!(source.clone().reframe_by(&first.identity()), Some(source.clone()));
        prop_assert_eq!(
            source.clone().reframe_by(first).and_then(|value| value.reframe_by(&inverse)),
            Some(source.clone()),
        );
        prop_assert_eq!(
            source.clone().reframe_by(first).and_then(|value| value.reframe_by(second)),
            source.reframe_by(&composite),
        );
    }

    #[test]
    fn test_noncovalent_bond_spans_reframe_by(scenario in standardization_scenario_strategy()) {
        let source = scenario.span.noncovalent_bonds().clone();
        let first = scenario.first_action.noncovalent_bonds();
        let second = scenario.second_action.noncovalent_bonds();
        let inverse = first.inverse();
        let composite = first.compose(second).expect("the actions share one domain");

        prop_assert_eq!(source.clone().reframe_by(&first.identity()), Some(source.clone()));
        prop_assert_eq!(
            source.clone().reframe_by(first).and_then(|value| value.reframe_by(&inverse)),
            Some(source.clone()),
        );
        prop_assert_eq!(
            source.clone().reframe_by(first).and_then(|value| value.reframe_by(second)),
            source.reframe_by(&composite),
        );
    }

    #[test]
    fn test_stereo_atom_spans_reframe_by(scenario in standardization_scenario_strategy()) {
        let source = scenario.span.stereo_atoms().clone();
        let first = scenario.first_action.stereo_atoms();
        let second = scenario.second_action.stereo_atoms();
        let inverse = first.inverse();
        let composite = first.compose(second).expect("the actions share one domain");

        prop_assert_eq!(source.clone().reframe_by(&first.identity()), Some(source.clone()));
        prop_assert_eq!(
            source.clone().reframe_by(first).and_then(|value| value.reframe_by(&inverse)),
            Some(source.clone()),
        );
        prop_assert_eq!(
            source.clone().reframe_by(first).and_then(|value| value.reframe_by(second)),
            source.reframe_by(&composite),
        );
    }

    #[test]
    fn test_stereo_bond_spans_reframe_by(scenario in standardization_scenario_strategy()) {
        let source = scenario.span.stereo_bonds().clone();
        let first = scenario.first_action.stereo_bonds();
        let second = scenario.second_action.stereo_bonds();
        let inverse = first.inverse();
        let composite = first.compose(second).expect("the actions share one domain");

        prop_assert_eq!(source.clone().reframe_by(&first.identity()), Some(source.clone()));
        prop_assert_eq!(
            source.clone().reframe_by(first).and_then(|value| value.reframe_by(&inverse)),
            Some(source.clone()),
        );
        prop_assert_eq!(
            source.clone().reframe_by(first).and_then(|value| value.reframe_by(second)),
            source.reframe_by(&composite),
        );
    }

    #[test]
    fn test_overlays_frame_action_compose(scenario in standardization_scenario_strategy()) {
        let first = &scenario.first_action;
        let second = &scenario.second_action;
        let inverse = first.inverse();
        let composite = first.compose(second).expect("the actions share one domain");
        let molecule = scenario.molecule;
        let reaction = Reaction::new(molecule.clone(), Default::default());
        let span = scenario.span;

        prop_assert_ne!(first, &first.identity());
        prop_assert_ne!(second, &second.identity());
        prop_assert_eq!(molecule.clone().reframe_by(&first.identity()), Some(molecule.clone()));
        prop_assert_eq!(
            molecule.clone().reframe_by(first).and_then(|value| value.reframe_by(&inverse)),
            Some(molecule.clone()),
        );
        prop_assert_eq!(
            molecule.clone().reframe_by(first).and_then(|value| value.reframe_by(second)),
            molecule.clone().reframe_by(&composite),
        );
        prop_assert_eq!(
            reaction.clone().reframe_by(first).and_then(|value| value.reframe_by(second)),
            reaction.clone().reframe_by(&composite),
        );
        prop_assert_eq!(
            span.clone().reframe_by(first).and_then(|value| value.reframe_by(second)),
            span.reframe_by(&composite),
        );
    }

    #[test]
    fn test_overlays_frame_action_domain(scenario in standardization_scenario_strategy()) {
        let missing = Reaction::default().representative_action();
        let reaction = Reaction::new(scenario.molecule.clone(), Default::default());

        prop_assert_eq!(
            scenario.molecule.clone().reframe_by(&scenario.covering_action),
            scenario.molecule.clone().reframe_by(&scenario.first_action),
        );
        prop_assert_eq!(scenario.molecule.clone().reframe_by(&missing), None);
        prop_assert_eq!(
            scenario.molecule.clone().reframe_by(&scenario.incompatible_action),
            None,
        );
        prop_assert_eq!(
            reaction.clone().reframe_by(&scenario.covering_action),
            reaction.clone().reframe_by(&scenario.first_action),
        );
        prop_assert_eq!(reaction.clone().reframe_by(&missing), None);
        prop_assert_eq!(reaction.reframe_by(&scenario.incompatible_action), None);
        prop_assert_eq!(
            scenario.span.clone().reframe_by(&scenario.covering_action),
            scenario.span.clone().reframe_by(&scenario.first_action),
        );
        prop_assert_eq!(scenario.span.clone().reframe_by(&missing), None);
        prop_assert_eq!(
            scenario
                .span
                .clone()
                .reframe_by(&scenario.incompatible_action),
            None,
        );
        prop_assert_eq!(scenario.first_action.compose(&missing), None);
        prop_assert_eq!(
            scenario.first_action.compose(&scenario.incompatible_action),
            None,
        );
        prop_assert!(scenario
            .covering_action
            .dative_bonds()
            .contains(DativeBondId(1)));
        prop_assert!(!scenario
            .span
            .dative_bonds()
            .contains(DativeBondId(1)));
        prop_assert!(scenario
            .reaction
            .representative_action()
            .multicenter_bonds()
            .contains(MulticenterBondId(7)));
    }
}
