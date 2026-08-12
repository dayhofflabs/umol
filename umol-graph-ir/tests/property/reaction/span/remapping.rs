//! Dense union-frame remapping properties for reaction spans.
//!
//! Generated spans include direct span entries and materializable reactions. Two independently
//! generated total correspondences exercise exact identity, inverse, composition, integrity, and
//! preservation of both compact side projections.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::Correspondence;
use umol_graph_ir::ir::{AtomId, EntitySpan, Molecule, MoleculeCorrespondence, ReactionSpan};

use super::reaction_span_scenario_strategy;

#[derive(Clone, Copy)]
enum Side {
    Left,
    Right,
}

fn is_present<T>(span: &EntitySpan<T>, side: Side) -> bool {
    match side {
        Side::Left => span.lhs().is_some(),
        Side::Right => span.rhs().is_some(),
    }
}

fn projected_atom_correspondence(
    span: &ReactionSpan,
    union: &Correspondence<AtomId>,
    side: Side,
) -> Correspondence<AtomId> {
    let source_present = span
        .atoms()
        .iter()
        .map(|attributes| is_present(attributes, side))
        .collect::<Vec<_>>();
    let mut target_present = vec![false; union.right_count()];
    for (source, present) in source_present.iter().copied().enumerate() {
        if present {
            let target = union
                .right_of(AtomId::from(source))
                .expect("union correspondence is total");
            target_present[target.0 as usize] = true;
        }
    }
    let mut target_projected = vec![None; union.right_count()];
    let mut target_count = 0;
    for (target, present) in target_present.into_iter().enumerate() {
        if present {
            target_projected[target] = Some(AtomId::from(target_count));
            target_count += 1;
        }
    }

    let mut pairs = Vec::new();
    let mut source_count = 0;
    for (source, present) in source_present.into_iter().enumerate() {
        if present {
            let target = union
                .right_of(AtomId::from(source))
                .expect("union correspondence is total");
            pairs.push((
                AtomId::from(source_count),
                target_projected[target.0 as usize]
                    .expect("present union atom has a projected target"),
            ));
            source_count += 1;
        }
    }
    Correspondence::new(pairs, source_count, target_count)
        .expect("projection of a union bijection is a bijection")
}

fn projection(span: &ReactionSpan, side: Side) -> Molecule {
    match side {
        Side::Left => span.lhs(),
        Side::Right => span.rhs(),
    }
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_reaction_span_remap_identity(
        scenario in reaction_span_scenario_strategy(),
    ) {
        let identity = scenario.first.compose(&scenario.first.reverse());

        prop_assert_eq!(scenario.span.remap(&identity), scenario.span);
    }

    #[test]
    fn test_reaction_span_remap_inverse(
        scenario in reaction_span_scenario_strategy(),
    ) {
        let remapped = scenario.span.remap(&scenario.first);
        let restored = remapped.remap(&scenario.first.reverse());

        prop_assert_eq!(restored, scenario.span);
    }

    #[test]
    fn test_reaction_span_remap_composition(
        scenario in reaction_span_scenario_strategy(),
    ) {
        let sequential = scenario.span.remap(&scenario.first).remap(&scenario.second);
        let direct = scenario.span.remap(&scenario.first.compose(&scenario.second));

        prop_assert_eq!(sequential, direct);
    }

    #[test]
    fn test_reaction_span_remap_integrity(
        scenario in reaction_span_scenario_strategy(),
    ) {
        let remapped = scenario.span.remap(&scenario.first);

        prop_assert_eq!(remapped.check_integrity(), Ok(()));
    }

    #[test]
    fn test_reaction_span_remap_projection(
        scenario in reaction_span_scenario_strategy(),
    ) {
        let remapped = scenario.span.remap(&scenario.first);
        for side in [Side::Left, Side::Right] {
            let source = projection(&scenario.span, side);
            let target = projection(&remapped, side);
            let atoms = projected_atom_correspondence(
                &scenario.span,
                scenario.first.atoms(),
                side,
            );
            let correspondence = MoleculeCorrespondence::induce(&source, &target, atoms)
                .expect("a remapped span induces its side correspondence");

            prop_assert!(correspondence.is_total());
            prop_assert_eq!(source.remap(&correspondence), target);
        }
    }
}
