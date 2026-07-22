//! Property tests for reaction composition and its DPO and canonicalization invariants.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_core::SubgraphIsomorphismAlgorithm;
use umol_utils::solution::Solution;

use crate::strategies::*;

const ALG: SubgraphIsomorphismAlgorithm = SubgraphIsomorphismAlgorithm::Vf2;

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    /// Every composite is itself a valid reaction: applying it at its own `lhs` reproduces its
    /// `right()`. Catches frame-algebra errors in the composite construction.
    #[test]
    fn test_reaction_ast_compose_well_formed(
        a in reaction_strategy(),
        b in reaction_strategy(),
    ) {
        for composite in a.compose(&b, CompositionScope::Full) {
            if let Ok(span) = composite.to_reaction_span() {
                let right = span.rhs();
                prop_assert!(composite
                    .apply(&composite.lhs, ALG)
                    .unwrap()
                    .any(|derivation| derivation.unwrap().rhs() == &right));
            }
        }
    }

    /// Soundness: every product of a composite applied to `A`'s reactant is also a product of
    /// applying B after A — `compose` invents no reactions.
    #[test]
    fn test_reaction_ast_compose_sound(
        a in reaction_strategy(),
        b in reaction_strategy(),
    ) {
        let host = a.lhs.clone();
        let composites = a.compose(&b, CompositionScope::Full);
        let composed: Vec<MoleculeAst> = composites
            .iter()
            .flat_map(|composite| {
                composite
                    .apply(&host, ALG)
                    .unwrap()
                    .map(Result::unwrap)
                    .collect::<Vec<_>>()
            })
            .map(|derivation| derivation.rhs().clone())
            .collect();

        let intermediates: Vec<MoleculeAst> = a
            .apply(&host, ALG)
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();
        let mut sequential: Vec<MoleculeAst> = Vec::new();
        for intermediate in &intermediates {
            sequential.extend(
                b.apply(intermediate, ALG)
                    .unwrap()
                    .map(Result::unwrap)
                    .map(|derivation| derivation.rhs().clone()),
            );
        }

        for product in &composed {
            prop_assert!(sequential.contains(product));
        }
    }

    /// Every composite of two overlay reactions is a valid reaction: applying it at its own `lhs`
    /// reproduces its `right()`. Catches overlay frame-algebra errors in the composite, and (the
    /// reason it once failed) `apply_at` removing multiple same-kind overlays: composites routinely
    /// remove ≥2 overlays of one kind, which the pre-batching single-id lowering mishandled.
    #[test]
    fn test_reaction_ast_compose_well_formed_overlay(
        a in overlay_reaction_strategy(),
        b in overlay_reaction_strategy(),
    ) {
        for composite in a.compose(&b, CompositionScope::Full) {
            if let Ok(span) = composite.to_reaction_span() {
                let right = span.rhs();
                prop_assert!(composite
                    .apply(&composite.lhs, ALG)
                    .unwrap()
                    .any(|derivation| derivation.unwrap().rhs() == &right));
            }
        }
    }

    /// Soundness with overlays: every product of a composite applied to A's reactant is also a
    /// product of applying B after A — compose invents no reactions, overlays included.
    #[test]
    fn test_reaction_ast_compose_sound_overlay(
        a in overlay_reaction_strategy(),
        b in overlay_reaction_strategy(),
    ) {
        let host = a.lhs.clone();
        let composed: Vec<MoleculeAst> = a
            .compose(&b, CompositionScope::Full)
            .iter()
            .flat_map(|composite| {
                composite
                    .apply(&host, ALG)
                    .unwrap()
                    .map(Result::unwrap)
                    .collect::<Vec<_>>()
            })
            .map(|derivation| derivation.rhs().clone())
            .collect();

        let intermediates: Vec<MoleculeAst> = a
            .apply(&host, ALG)
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();
        let mut sequential: Vec<MoleculeAst> = Vec::new();
        for intermediate in &intermediates {
            sequential.extend(
                b.apply(intermediate, ALG)
                    .unwrap()
                    .map(Result::unwrap)
                    .map(|derivation| derivation.rhs().clone()),
            );
        }

        for product in &composed {
            prop_assert!(sequential.contains(product));
        }
    }

    /// P1 completeness (`Full`): every sequential product (B applied after A) is also some
    /// composite's product. Together with `compose_sound_overlay` (composed ⊆ seq) this is set
    /// equality at `host = a.lhs`. Covers stereo: the reactants carry stereo overlays and the deltas
    /// stereo ops, glued and applied across ligand frames by `meet_pushout` / `apply_at`.
    #[test]
    fn test_reaction_ast_compose_complete_overlay(
        a in overlay_reaction_strategy(),
        b in overlay_reaction_strategy(),
    ) {
        let host = a.lhs.clone();
        let composed: Vec<MoleculeAst> = a
            .compose(&b, CompositionScope::Full)
            .iter()
            .flat_map(|composite| {
                composite
                    .apply(&host, ALG)
                    .unwrap()
                    .map(Result::unwrap)
                    .collect::<Vec<_>>()
            })
            .map(|derivation| derivation.rhs().clone())
            .collect();

        let intermediates: Vec<MoleculeAst> = a
            .apply(&host, ALG)
            .unwrap()
            .map(Result::unwrap)
            .map(|derivation| derivation.rhs().clone())
            .collect();
        let mut sequential: Vec<MoleculeAst> = Vec::new();
        for intermediate in &intermediates {
            sequential.extend(
                b.apply(intermediate, ALG)
                    .unwrap()
                    .map(Result::unwrap)
                    .map(|derivation| derivation.rhs().clone()),
            );
        }

        for product in &sequential {
            prop_assert!(
                composed.contains(product),
                "sequential product missing from composed set (P1 completeness)"
            );
        }
    }

    /// Every composite is DPO-valid — no deleted atom leaves a dangling bond or overlay. Confirms
    /// the compose during-check yields dangling-free composites (via the tier-2 `DpoValidator`).
    #[test]
    fn test_reaction_ast_compose_dangling_free(
        a in overlay_reaction_strategy(),
        b in overlay_reaction_strategy(),
    ) {
        for composite in a.compose(&b, CompositionScope::Full) {
            prop_assert_eq!(
                DpoValidator
                    .validate_reaction(&composite.lhs, &composite.deltas)
                    .unwrap(),
                Solution::Determined(())
            );
        }
    }

    /// `RcAnchored` is a sound filter: every RC-anchored composite is also a `Full` composite.
    #[test]
    fn test_reaction_ast_compose_rc_anchored_subset(
        a in overlay_reaction_strategy(),
        b in overlay_reaction_strategy(),
    ) {
        let full = a.compose(&b, CompositionScope::Full);
        for composite in a.compose(&b, CompositionScope::RcAnchored) {
            prop_assert!(full.contains(&composite));
        }
    }

    /// P4 — determinism: `compose` returns the identical `Vec` on repeated calls, and is invariant
    /// under pre-canonicalizing the inputs (compose canonicalizes the deltas itself).
    #[test]
    fn test_reaction_ast_compose_determinism(
        a in overlay_reaction_strategy(),
        b in overlay_reaction_strategy(),
    ) {
        prop_assert_eq!(
            a.compose(&b, CompositionScope::Full),
            a.compose(&b, CompositionScope::Full)
        );
        prop_assert_eq!(
            a.compose(&b, CompositionScope::RcAnchored),
            a.compose(&b, CompositionScope::RcAnchored)
        );
        if let (Ok(ac), Ok(bc)) = (a.clone().canonicalize(), b.clone().canonicalize()) {
            prop_assert_eq!(
                a.compose(&b, CompositionScope::Full),
                ac.compose(&bc, CompositionScope::Full)
            );
        }
    }

    /// P3 — every composite's deltas are in canonical normal form.
    #[test]
    fn test_reaction_ast_compose_canonical_deltas(
        a in overlay_reaction_strategy(),
        b in overlay_reaction_strategy(),
    ) {
        for c in a.compose(&b, CompositionScope::Full) {
            let canonical = c
                .deltas
                .clone()
                .canonicalize()
                .map_err(|e| TestCaseError::fail(format!("composite deltas not canonical: {e:?}")))?;
            prop_assert_eq!(canonical, c.deltas);
        }
    }

    /// P6 — no parallel overlays: within each kind a composite's overlays have distinct participant
    /// sets, so correspondence reuses an id and never duplicates (spec §4.1).
    #[test]
    fn test_reaction_ast_compose_distinct_overlays(
        a in overlay_reaction_strategy(),
        b in overlay_reaction_strategy(),
    ) {
        for c in a.compose(&b, CompositionScope::Full) {
            let m = &c.lhs;

            let mut dative: Vec<(Vec<AtomId>, AtomId)> = m
                .dative_bonds()
                .iter()
                .map(|x| {
                    let mut donors: Vec<AtomId> = x.donor_ids().collect();
                    donors.sort();
                    (donors, x.acceptor_id())
                })
                .collect();
            let dative_count = dative.len();
            dative.sort();
            dative.dedup();
            prop_assert_eq!(dative.len(), dative_count, "duplicate dative bonds");

            let mut aromatic: Vec<Vec<AtomId>> = m
                .aromatic_systems()
                .iter()
                .map(|x| {
                    let mut v: Vec<AtomId> = x.atom_ids().collect();
                    v.sort();
                    v
                })
                .collect();
            let aromatic_count = aromatic.len();
            aromatic.sort();
            aromatic.dedup();
            prop_assert_eq!(aromatic.len(), aromatic_count, "duplicate aromatic systems");

            let mut multicenter: Vec<Vec<AtomId>> = m
                .multicenter_bonds()
                .iter()
                .map(|x| {
                    let mut v: Vec<AtomId> = x.atom_ids().collect();
                    v.sort();
                    v
                })
                .collect();
            let multicenter_count = multicenter.len();
            multicenter.sort();
            multicenter.dedup();
            prop_assert_eq!(
                multicenter.len(),
                multicenter_count,
                "duplicate multicenter bonds"
            );

            let mut noncovalent: Vec<[AtomId; 2]> = m
                .noncovalent_bonds()
                .iter()
                .map(|x| {
                    let mut p = x.atom_ids();
                    p.sort();
                    p
                })
                .collect();
            let noncovalent_count = noncovalent.len();
            noncovalent.sort();
            noncovalent.dedup();
            prop_assert_eq!(
                noncovalent.len(),
                noncovalent_count,
                "duplicate noncovalent bonds"
            );
        }
    }
}
