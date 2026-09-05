//! Reaction canonicalization properties.
//!
//! Materializable reactions and independently renumbered reaction spans exercise exact canonical
//! forms, the complete normalization/reframe/canonicalize law matrix, equality, canonical hashes,
//! reversal, and covariant application. Named defects separately cover discontinuous deltas and
//! intrinsically contradictory forms. Nauty is currently the only canonical-labeling selector;
//! frozen canonical fixtures, rather than a tautological second case, remain the compatibility
//! target for future algorithms.

use std::hash::{DefaultHasher, Hash, Hasher};

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::{
    AutomorphismAlgorithm, Correspondence, EdgeId, GraphRemapping, NodeId, Remapping,
};
use umol_graph_ir::ir::{
    Canonicalize, CanonicalizeContext, Contradiction, Deltas, EntitySpan, Molecule,
    MoleculeCorrespondence, MoleculeEntries, MoleculeRemapping, Normalize, NumForm, Reaction,
    ReactionCanonicalizeError, ReactionSpan, Reframe,
};

use super::span::reaction_span_scenario_strategy;
use crate::strategies::*;

#[derive(Clone, Copy)]
enum Side {
    Left,
    Right,
}

fn context() -> CanonicalizeContext {
    CanonicalizeContext {
        para_stereo: false,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    }
}

fn structural_hash<T: Hash>(value: &T) -> u64 {
    let mut hasher = DefaultHasher::new();
    value.hash(&mut hasher);
    hasher.finish()
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
    application: &(Molecule, MoleculeCorrespondence),
) -> Correspondence<AtomId> {
    let (product, correspondence) = application;
    let side = span.correspondence();
    let mut pairs = Vec::new();
    for &(left, right) in side.atoms().matched_pairs() {
        let host = pattern_to_host
            .atoms()
            .right_of(left)
            .expect("application correspondence is total on the pattern");
        if let Some(product) = correspondence.atoms().right_of(host) {
            pairs.push((right, product));
        }
    }
    pairs.extend(
        side.atoms()
            .right_unmatched()
            .into_iter()
            .zip(correspondence.atoms().right_unmatched()),
    );
    Correspondence::new(pairs, side.atoms().right_count(), product.atoms().count())
        .expect("reaction rhs embeds injectively in its application product")
}

fn product_correspondence(
    source_span: &ReactionSpan,
    canonical_span: &ReactionSpan,
    union: &MoleculeCorrespondence,
    source_match: &MoleculeCorrespondence,
    canonical_match: &MoleculeCorrespondence,
    source: &(Molecule, MoleculeCorrespondence),
    canonical: &(Molecule, MoleculeCorrespondence),
) -> MoleculeCorrespondence {
    let rhs_action = projected_atom_correspondence(source_span, union.atoms(), Side::Right);
    let source_rhs_to_product = rhs_to_product_atoms(source_span, source_match, source);
    let canonical_rhs_to_product = rhs_to_product_atoms(canonical_span, canonical_match, canonical);
    let mut pairs = Vec::new();
    for host in (0..source.1.atoms().left_count()).map(AtomId::from) {
        match (
            source.1.atoms().right_of(host),
            canonical.1.atoms().right_of(host),
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
    let atoms = Correspondence::new(pairs, source.0.atoms().count(), canonical.0.atoms().count())
        .expect("corresponding products have a bijective atom transport");
    MoleculeCorrespondence::induce(&source.0, &canonical.0, atoms)
        .expect("product atom transport induces all entity kinds")
}

fn canonicalization_error_strategy() -> impl Strategy<Value = (Reaction, ReactionCanonicalizeError)>
{
    prop_oneof![
        discontinuous_atom_update_reaction_strategy().prop_map(|reaction| {
            (
                reaction,
                ReactionCanonicalizeError::Contradiction(Contradiction),
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
                ReactionCanonicalizeError::Contradiction(Contradiction),
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
        let (with_correspondence, correspondence) = reaction
            .clone()
            .tracked_canonicalize(&context)
            .map_err(|error| {
                TestCaseError::fail(format!(
                    "generated reaction did not canonicalize with a correspondence: {error}"
                ))
            })?;
        let reconstructed = reaction
            .to_reaction_span()
            .map_err(|_| TestCaseError::fail("generated reaction did not materialize"))?
            .remap(&correspondence)
            .reframe()
            .map_err(|_| {
                TestCaseError::fail("canonical correspondence produced a contradiction")
            })?
            .to_reaction();

        prop_assert_eq!(
            Reaction::try_new(canonical.lhs().clone(), canonical.deltas().clone()),
            Ok(canonical.clone()),
        );
        prop_assert_eq!(&with_correspondence, &canonical);
        prop_assert_eq!(reconstructed, canonical);
    }

    #[test]
    fn test_reaction_canonicalize_standardization(
        scenario in standardization_scenario_strategy(),
    ) {
        let context = context();
        let source = scenario.reaction;
        let identical = source.clone();
        let normalized = source
            .clone()
            .normalize()
            .map_err(|_| TestCaseError::fail("generated reaction is intrinsically contradictory"))?;
        let normalized_again = normalized
            .clone()
            .normalize()
            .map_err(|_| TestCaseError::fail("normalized reaction became contradictory"))?;
        let reframed = source
            .clone()
            .reframe()
            .map_err(|_| TestCaseError::fail("generated reaction is intrinsically contradictory"))?;
        let canonical = source.clone().canonicalize(&context).map_err(|error| {
            TestCaseError::fail(format!("generated reaction did not canonicalize: {error}"))
        })?;
        let renumbered = scenario.span.remap(&scenario.remapping).to_reaction();

        prop_assert_eq!(normalized.clone().normalize(), Ok(normalized.clone()));
        prop_assert_eq!(reframed.clone().reframe(), Ok(reframed.clone()));
        prop_assert_eq!(
            canonical.clone().canonicalize(&context),
            Ok(canonical.clone())
        );
        prop_assert_eq!(normalized.clone().reframe(), Ok(reframed.clone()));
        prop_assert_eq!(reframed.clone().normalize(), Ok(reframed.clone()));
        prop_assert_eq!(
            normalized.clone().canonicalize(&context),
            Ok(canonical.clone())
        );
        prop_assert_eq!(
            reframed.clone().canonicalize(&context),
            Ok(canonical.clone())
        );
        prop_assert_eq!(canonical.clone().normalize(), Ok(canonical.clone()));
        prop_assert_eq!(canonical.clone().reframe(), Ok(canonical.clone()));

        prop_assert_eq!(&source, &identical);
        prop_assert!(source.normalized_eq(&identical));
        prop_assert!(source.framed_eq(&identical));
        prop_assert!(source.canonical_eq(&identical, &context));
        prop_assert!(source.normalized_eq(&normalized));
        prop_assert_eq!(
            source.normalized_eq(&normalized),
            normalized.normalized_eq(&source),
        );
        prop_assert!(normalized.normalized_eq(&normalized_again));
        prop_assert!(source.normalized_eq(&normalized_again));
        prop_assert!(source.framed_eq(&normalized));
        prop_assert!(normalized.framed_eq(&reframed));
        prop_assert!(source.framed_eq(&reframed));
        prop_assert_eq!(source.framed_eq(&reframed), reframed.framed_eq(&source));
        prop_assert!(source.canonical_eq(&reframed, &context));
        prop_assert_eq!(
            source.canonical_eq(&reframed, &context),
            reframed.canonical_eq(&source, &context),
        );
        prop_assert!(reframed.canonical_eq(&renumbered, &context));
        prop_assert!(source.canonical_eq(&renumbered, &context));
        prop_assert!(renumbered.canonical_eq(&canonical, &context));
    }

    #[test]
    fn test_reaction_canonical_hash(scenario in reaction_span_scenario_strategy()) {
        let context = context();
        let source = scenario.span.to_reaction();
        let renumbered = scenario.span.remap(&scenario.first).to_reaction();

        prop_assert_eq!(
            source.clone().canonical_hash(&context),
            renumbered.clone().canonical_hash(&context),
        );
        if let Ok(canonical) = source.clone().canonicalize(&context) {
            prop_assert_eq!(
                source.canonical_hash(&context),
                Ok(structural_hash(&canonical)),
            );
        }
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
        let (canonical, union) = reaction
            .clone()
            .tracked_canonicalize(&context)
            .map_err(|error| {
                TestCaseError::fail(format!("generated reaction did not canonicalize: {error}"))
            })?;
        let canonical_span = canonical
            .to_reaction_span()
            .expect("canonical reaction materializes");
        let union = MoleculeCorrespondence::from(&union);
        let lhs_atoms = projected_atom_correspondence(&source_span, union.atoms(), Side::Left);
        let lhs_action = MoleculeCorrespondence::induce(reaction.lhs(), canonical.lhs(), lhs_atoms)
            .expect("canonical union action induces its lhs action");
        let canonical_match = lhs_action.reverse().compose(&correspondence).unwrap();
        let source_application = reaction.tracked_apply_at(&host, &correspondence);
        let canonical_application = canonical.tracked_apply_at(&host, &canonical_match);

        match (source_application, canonical_application) {
            // Application is a partial operation. Canonical relabeling preserves its domain, but
            // the first diagnostic need not be stable when several embeddings or stereo frames
            // fail and their entity order changes.
            (Err(_), Err(_)) => {}
            (Ok(None), Ok(None)) => {}
            (Ok(Some(left)), Ok(Some(right))) => {
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
                let remapping = MoleculeRemapping::new(
                    GraphRemapping::new(
                        Remapping::new(
                            products
                                .atoms()
                                .matched_pairs()
                                .iter()
                                .map(|&(_, right)| NodeId::from(right))
                                .collect(),
                        )
                        .unwrap(),
                        Remapping::new(
                            products
                                .bonds()
                                .matched_pairs()
                                .iter()
                                .map(|&(_, right)| EdgeId::from(right))
                                .collect(),
                        )
                        .unwrap(),
                    ),
                    Remapping::new(
                        products
                            .dative_bonds()
                            .matched_pairs()
                            .iter()
                            .map(|&(_, right)| right)
                            .collect(),
                    )
                    .unwrap(),
                    Remapping::new(
                        products
                            .aromatic_systems()
                            .matched_pairs()
                            .iter()
                            .map(|&(_, right)| right)
                            .collect(),
                    )
                    .unwrap(),
                    Remapping::new(
                        products
                            .multicenter_bonds()
                            .matched_pairs()
                            .iter()
                            .map(|&(_, right)| right)
                            .collect(),
                    )
                    .unwrap(),
                    Remapping::new(
                        products
                            .noncovalent_bonds()
                            .matched_pairs()
                            .iter()
                            .map(|&(_, right)| right)
                            .collect(),
                    )
                    .unwrap(),
                    Remapping::new(
                        products
                            .stereo_atoms()
                            .matched_pairs()
                            .iter()
                            .map(|&(_, right)| right)
                            .collect(),
                    )
                    .unwrap(),
                    Remapping::new(
                        products
                            .stereo_bonds()
                            .matched_pairs()
                            .iter()
                            .map(|&(_, right)| right)
                            .collect(),
                    )
                    .unwrap(),
                );
                prop_assert!(left.0.framed_eq_under(&right.0, &remapping));
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

    #[test]
    fn test_reaction_canonical_eq_contradiction(
        scenario in intrinsic_contradiction_scenario_strategy(),
    ) {
        let context = context();
        let [first, second, third] = scenario.reactions;

        prop_assert_ne!(&first, &second);
        prop_assert_ne!(&second, &third);
        prop_assert!(first.normalized_eq(&second));
        prop_assert!(second.normalized_eq(&third));
        prop_assert!(first.normalized_eq(&third));
        prop_assert!(first.framed_eq(&second));
        prop_assert!(second.framed_eq(&third));
        prop_assert!(first.framed_eq(&third));
        prop_assert!(first.canonical_eq(&second, &context));
        prop_assert!(second.canonical_eq(&third, &context));
        prop_assert!(first.canonical_eq(&third, &context));
        prop_assert_eq!(
            first.tracked_canonicalize(&context),
            Err(ReactionCanonicalizeError::Contradiction(Contradiction)),
        );
    }
}
