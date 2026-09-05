//! Aggregate molecule frame-transport and reframe properties.
//!
//! The generated domain contains integrity-valid molecules with all overlay kinds, positional
//! payloads, and recursive constraints. Transport is checked through its action algebra, while
//! reframing is checked by agreement of its fused and witness-returning executions. These overlap
//! deliberately: the former isolates a supplied action and the latter validates representative
//! selection plus the normalization prefix. Pipeline fixpoint and absorption laws are checked in
//! the canonicalization module.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph_ir::ir::{Contradiction, FrameTransport, Normalize, Reframe};

use crate::strategies::{
    intrinsic_contradiction_scenario_strategy, molecule_with_constraints_strategy,
};

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_molecule_reframe_by(molecule in molecule_with_constraints_strategy()) {
        let action = molecule.representative_action();
        let identity = action.identity();
        let inverse = action.inverse();
        let composite = action
            .compose(&inverse)
            .expect("an action and its inverse have the same domain");

        prop_assert_eq!(
            molecule.clone().reframe_by(&identity),
            Some(molecule.clone()),
        );
        prop_assert_eq!(
            molecule
                .clone()
                .reframe_by(&action)
                .and_then(|transported| transported.reframe_by(&inverse)),
            molecule.clone().reframe_by(&composite),
        );
        prop_assert_eq!(molecule.clone().reframe_by(&composite), Some(molecule));
    }

    #[test]
    fn test_molecule_representative_action_contradiction(
        scenario in intrinsic_contradiction_scenario_strategy(),
    ) {
        for molecule in scenario.molecules {
            let action = molecule.representative_action();

            prop_assert_eq!(
                action.compose(&action.identity()),
                Some(action.clone()),
            );
            prop_assert_eq!(molecule.clone().normalize(), Err(Contradiction));
            prop_assert_eq!(molecule.reframe(), Err(Contradiction));
        }
    }

    #[test]
    fn test_molecule_tracked_reframe(molecule in molecule_with_constraints_strategy()) {
        let fused = molecule.clone().reframe().map_err(|_| {
            TestCaseError::fail("generated molecule is intrinsically contradictory")
        })?;
        let (witnessed, action) = molecule.clone().tracked_reframe().map_err(|_| {
            TestCaseError::fail("generated molecule is intrinsically contradictory")
        })?;
        let transported = molecule
            .normalize()
            .map_err(|_| {
                TestCaseError::fail("generated molecule is intrinsically contradictory")
            })?
            .reframe_by(&action)
            .ok_or_else(|| TestCaseError::fail("representative action did not cover its source"))?
            .normalize()
            .map_err(|_| {
                TestCaseError::fail("transported molecule is intrinsically contradictory")
            })?;
        let selected_action = witnessed.representative_action();

        prop_assert_eq!(fused, witnessed.clone());
        prop_assert_eq!(transported, witnessed);
        prop_assert_eq!(selected_action.clone(), selected_action.identity());
    }

    #[test]
    fn test_molecule_framed_eq(molecule in molecule_with_constraints_strategy()) {
        let normalized = molecule.clone().normalize().map_err(|_| {
            TestCaseError::fail("generated molecule is intrinsically contradictory")
        })?;
        let reframed = molecule.clone().reframe().map_err(|_| {
            TestCaseError::fail("generated molecule is intrinsically contradictory")
        })?;

        prop_assert!(molecule.normalized_eq(&normalized));
        prop_assert!(molecule.framed_eq(&normalized));
        prop_assert!(molecule.framed_eq(&reframed));
    }
}
