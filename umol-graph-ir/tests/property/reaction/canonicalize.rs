//! Reaction canonicalization properties.
//!
//! Materializable reactions and independently renumbered reaction spans exercise exact canonical
//! forms, equality, normalization, reversal, and covariant application. Named defects separately
//! cover integrity failure, discontinuous deltas, and intrinsically contradictory forms. Nauty is
//! currently the only canonical-labeling selector; frozen canonical fixtures, rather than a
//! tautological second case, remain the compatibility target for future algorithms.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::{AutomorphismAlgorithm, Correspondence};
use umol_graph_ir::ir::canonicalize::canonicalize_reaction_with_correspondence;
use umol_graph_ir::ir::{
    AtomDelta, CanonicalizationContext, CanonicalizationLevel, Canonicalize, Contradiction, Delta,
    Deltas, Entity, EntitySpan, Molecule, MoleculeCorrespondence, MoleculeEntries, NumForm,
    Reaction, ReactionCanonicalizationError, ReactionDerivation, ReactionIntegrityError,
    ReactionSpan,
};

use super::span::reaction_span_scenario_strategy;
use crate::strategies::*;

#[derive(Clone, Copy)]
enum Side {
    Left,
    Right,
}

fn context() -> CanonicalizationContext {
    CanonicalizationContext {
        para_stereo: false,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    }
}

fn is_present<T>(value: &EntitySpan<T>, side: Side) -> bool {
    match side {
        Side::Left => value.lhs().is_some(),
        Side::Right => value.rhs().is_some(),
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
                .expect("canonical union correspondence is total");
            target_present[target.index()] = true;
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
                .expect("canonical union correspondence is total");
            pairs.push((
                AtomId::from(source_count),
                target_projected[target.index()]
                    .expect("present union atom has a projected target"),
            ));
            source_count += 1;
        }
    }
    Correspondence::new(pairs, source_count, target_count)
        .expect("projection of a union bijection is a bijection")
}

fn rhs_to_product_atoms(
    span: &ReactionSpan,
    pattern_to_host: &MoleculeCorrespondence,
    derivation: &ReactionDerivation,
) -> Correspondence<AtomId> {
    let side = span.correspondence();
    let mut pairs = Vec::new();
    for &(left, right) in side.atoms().matched_pairs() {
        let host = pattern_to_host
            .atoms()
            .right_of(left)
            .expect("application correspondence is total on the pattern");
        if let Some(product) = derivation.comap().atoms().right_of(host) {
            pairs.push((right, product));
        }
    }
    pairs.extend(
        side.atoms()
            .right_unmatched()
            .into_iter()
            .zip(derivation.comap().atoms().right_unmatched()),
    );
    Correspondence::new(
        pairs,
        side.atoms().right_count(),
        derivation.rhs().atoms().count(),
    )
    .expect("reaction rhs embeds injectively in its application product")
}

fn product_correspondence(
    source_span: &ReactionSpan,
    canonical_span: &ReactionSpan,
    union: &MoleculeCorrespondence,
    source_match: &MoleculeCorrespondence,
    canonical_match: &MoleculeCorrespondence,
    source: &ReactionDerivation,
    canonical: &ReactionDerivation,
) -> MoleculeCorrespondence {
    let rhs_action = projected_atom_correspondence(source_span, union.atoms(), Side::Right);
    let source_rhs_to_product = rhs_to_product_atoms(source_span, source_match, source);
    let canonical_rhs_to_product = rhs_to_product_atoms(canonical_span, canonical_match, canonical);
    let mut pairs = Vec::new();
    for host in (0..source.comap().atoms().left_count()).map(AtomId::from) {
        match (
            source.comap().atoms().right_of(host),
            canonical.comap().atoms().right_of(host),
        ) {
            (Some(left), Some(right)) => pairs.push((left, right)),
            (None, None) => {}
            _ => panic!("corresponding applications preserve the same host atoms"),
        }
    }
    for source_rhs in source_span.correspondence().atoms().right_unmatched() {
        let canonical_rhs = rhs_action
            .right_of(source_rhs)
            .expect("canonical action transports every rhs atom");
        pairs.push((
            source_rhs_to_product
                .right_of(source_rhs)
                .expect("source application contains every created atom"),
            canonical_rhs_to_product
                .right_of(canonical_rhs)
                .expect("canonical application contains every created atom"),
        ));
    }
    let atoms = Correspondence::new(
        pairs,
        source.rhs().atoms().count(),
        canonical.rhs().atoms().count(),
    )
    .expect("corresponding products have a bijective atom transport");
    MoleculeCorrespondence::induce(source.rhs(), canonical.rhs(), atoms)
        .expect("product atom transport induces all entity families")
}

fn canonicalization_error_strategy(
) -> impl Strategy<Value = (Reaction, ReactionCanonicalizationError)> {
    prop_oneof![
        (0u32..64).prop_map(|id| {
            (
                Reaction::new(
                    Molecule::default(),
                    [Delta::Atom(AtomDelta::Remove {
                        id: AtomId(id),
                        attributes: AtomForm::default(),
                    })]
                    .into_iter()
                    .collect(),
                ),
                ReactionCanonicalizationError::Integrity(
                    ReactionIntegrityError::InvalidReference {
                        entity: Entity::Atom(AtomId(id)),
                    },
                ),
            )
        }),
        discontinuous_atom_update_reaction_strategy().prop_map(|reaction| {
            (
                reaction,
                ReactionCanonicalizationError::Contradiction(Contradiction),
            )
        }),
        any::<i64>().prop_map(|charge| {
            let atom = AtomForm::default()
                .with_charge(NumForm::lit_set(Vec::<i64>::new()))
                .with_implicit_hydrogens(charge);
            (
                Reaction::new(
                    Molecule::from_entries(MoleculeEntries {
                        atoms: vec![atom],
                        ..Default::default()
                    }),
                    Deltas::new(),
                ),
                ReactionCanonicalizationError::Contradiction(Contradiction),
            )
        }),
    ]
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_reaction_canonicalize(reaction in materializable_reaction_strategy()) {
        let context = context();
        let canonical = reaction.clone().canonicalize(&context).map_err(|error| {
            TestCaseError::fail(format!("generated reaction did not canonicalize: {error}"))
        })?;

        prop_assert_eq!(canonical.check_integrity(), Ok(()));
        prop_assert_eq!(canonical.clone().canonicalize(&context), Ok(canonical));
    }

    #[test]
    fn test_reaction_canonicalize_by(scenario in reaction_span_scenario_strategy()) {
        let context = context();
        let source = scenario.span.to_reaction();
        let renumbered = scenario.span.remap(&scenario.first).to_reaction();

        for level in [
            CanonicalizationLevel::Topology,
            CanonicalizationLevel::Constitution,
            CanonicalizationLevel::Structure,
            CanonicalizationLevel::Full,
        ] {
            let left = source.clone().canonicalize_by(level, &context);
            let right = renumbered.clone().canonicalize_by(level, &context);
            match (left, right) {
                (Ok(left), Ok(right)) => {
                    prop_assert!(left.canonical_eq_by(&right, level, &context));
                    if level == CanonicalizationLevel::Full {
                        prop_assert_eq!(left, right);
                    }
                }
                (Err(left), Err(right)) => prop_assert_eq!(left, right),
                (left, right) => {
                    prop_assert!(false, "canonicalization mismatch: {left:?} != {right:?}")
                }
            }
        }
        prop_assert_eq!(
            source.clone().canonicalize_by(CanonicalizationLevel::Full, &context),
            source.canonicalize(&context),
        );
    }

    #[test]
    fn test_reaction_canonical_eq_by(scenario in reaction_span_scenario_strategy()) {
        let context = context();
        let source = scenario.span.to_reaction();
        let renumbered = scenario.span.remap(&scenario.first).to_reaction();

        for level in [
            CanonicalizationLevel::Topology,
            CanonicalizationLevel::Constitution,
            CanonicalizationLevel::Structure,
            CanonicalizationLevel::Full,
        ] {
            prop_assert!(source.canonical_eq_by(&source, level, &context));
            prop_assert!(source.canonical_eq_by(&renumbered, level, &context));
            prop_assert_eq!(
                source.canonical_eq_by(&renumbered, level, &context),
                renumbered.canonical_eq_by(&source, level, &context),
            );
        }
        prop_assert_eq!(
            source.canonical_eq_by(&renumbered, CanonicalizationLevel::Full, &context),
            source.canonical_eq(&renumbered, &context),
        );
    }

    #[test]
    fn test_reaction_canonicalize_roundtrip(reaction in materializable_reaction_strategy()) {
        let context = context();
        let normalized = reaction
            .to_reaction_span()
            .expect("generated reaction materializes")
            .to_reaction();

        prop_assert_eq!(
            normalized.canonicalize(&context),
            reaction.canonicalize(&context),
        );
    }

    #[test]
    fn test_reaction_canonicalize_reversal(reaction in materializable_reaction_strategy()) {
        let context = context();
        let canonical = reaction.clone().canonicalize(&context).map_err(|error| {
            TestCaseError::fail(format!("generated reaction did not canonicalize: {error}"))
        })?;
        let canonical_reversed = canonical
            .reverse()
            .map_err(|error| TestCaseError::fail(format!("canonical reversal failed: {error}")))?
            .canonicalize(&context);
        let reversed_canonical = reaction
            .reverse()
            .map_err(|error| TestCaseError::fail(format!("reaction reversal failed: {error}")))?
            .canonicalize(&context);

        prop_assert_eq!(canonical_reversed, reversed_canonical);
    }

    #[test]
    fn test_reaction_canonicalize_application(
        (reaction, host, correspondence) in reaction_application_strategy(),
    ) {
        let context = context();
        let source_span = reaction
            .to_reaction_span()
            .expect("generated reaction materializes");
        let (canonical, union) = canonicalize_reaction_with_correspondence(&reaction, &context)
            .map_err(|error| {
                TestCaseError::fail(format!("generated reaction did not canonicalize: {error}"))
            })?;
        let canonical_span = canonical
            .to_reaction_span()
            .expect("canonical reaction materializes");
        let lhs_atoms = projected_atom_correspondence(&source_span, union.atoms(), Side::Left);
        let lhs_action = MoleculeCorrespondence::induce(&reaction.lhs, &canonical.lhs, lhs_atoms)
            .expect("canonical union action induces its lhs action");
        let canonical_match = lhs_action.reverse().compose(&correspondence);
        let source_application = reaction.apply_at(&host, &correspondence);
        let canonical_application = canonical.apply_at(&host, &canonical_match);

        match (source_application, canonical_application) {
            // Application is a partial operation. Canonical relabeling preserves its domain, but
            // the first diagnostic need not be stable when several embeddings or stereo frames
            // fail and their entity order changes.
            (Err(_), Err(_)) => {}
            (Ok(left), Ok(right)) => {
                let products = product_correspondence(
                    &source_span,
                    &canonical_span,
                    &union,
                    &correspondence,
                    &canonical_match,
                    &left,
                    &right,
                );
                prop_assert!(products.is_total());
                prop_assert!(left.rhs().equiv_under(right.rhs(), &products));
            }
            (left, right) => prop_assert!(false, "application mismatch: {left:?} != {right:?}"),
        }
    }

    #[test]
    fn test_reaction_canonicalize_error(
        (reaction, expected) in canonicalization_error_strategy(),
    ) {
        prop_assert_eq!(reaction.canonicalize(&context()), Err(expected));
    }
}
