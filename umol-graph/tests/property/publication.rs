//! Representation-integrity preservation by chemistry-layer resolvers.
//!
//! Raw C/N/O/S ring skeletons retain contradictions and undetermined fields. Stored
//! first-ring skeletons supply literal contributions independently of perception.
//! Saturated carbon skeletons supply an independently concrete resolution domain.
//! Each operation checks consuming/in-place outcome agreement, publication through
//! the editor's integrity gate, and preservation of the input on rejection.

use proptest::prelude::*;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_graph::ops::model::{AromaticityModel, ChemistryModel, ValenceModel};
use umol_graph::ops::resolve::{
    AromaticityFailurePolicy, AromaticityResolveConfig, AromaticityResolver, IsotopePolicy,
    ResolveConfig, Resolver,
};
use umol_graph_ir::ir::{Normalize, RollbackMethod};
use umol_utils::solution::Solution;

use super::strategies::publication_molecule;

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::Direct(
            super::REGRESSION_FILE,
        ))),
        ..Config::default()
    })]

    #[test]
    fn test_aromaticity_resolver_resolve_publication(
        source in publication_molecule(),
        keep_valence in any::<bool>(),
        keep_system in any::<bool>(),
    ) {
        let resolver = AromaticityResolver::with_config(&AromaticityModel::daylight(),
            AromaticityResolveConfig {
                aromatic_valence_failure: if keep_valence { AromaticityFailurePolicy::Keep }
                    else { AromaticityFailurePolicy::Error },
                aromatic_system_failure: if keep_system { AromaticityFailurePolicy::Keep }
                    else { AromaticityFailurePolicy::Error },
                ..Default::default()
            });
        let mut borrowed = source.clone();
        let consumed = resolver.resolve(source.clone());
        let in_place = resolver.resolve_into(&mut borrowed);
        let mut projections = vec![source.clone()];
        match consumed {
            Ok(Solution::Determined(published)) => {
                prop_assert_eq!(in_place, Ok(Solution::Determined(())));
                prop_assert_eq!(&borrowed, &published);
                prop_assert_eq!(published.clone().edit().finish(), Ok(published.clone()));
                projections.push(published);
            }
            Ok(Solution::Underdetermined(())) => {
                prop_assert_eq!(in_place, Ok(Solution::Underdetermined(())));
                prop_assert_eq!(borrowed, source);
            }
            Ok(Solution::Contradictory(contradiction)) => {
                prop_assert_eq!(in_place, Ok(Solution::Contradictory(contradiction)));
                prop_assert_eq!(borrowed, source);
            }
            Err(error) => {
                prop_assert_eq!(in_place, Err(error));
                prop_assert_eq!(borrowed, source);
            }
        }
        for source in projections {
            let mut borrowed = source.clone();
            let consumed = resolver.project(source.clone());
            let in_place = resolver.project_into(&mut borrowed);
            match consumed {
                Ok(Solution::Determined(published)) => {
                    prop_assert_eq!(in_place, Ok(Solution::Determined(())));
                    prop_assert_eq!(&borrowed, &published);
                    prop_assert_eq!(published.clone().edit().finish(), Ok(published));
                }
                Ok(outcome) => return Err(TestCaseError::fail(format!(
                    "projection has no non-determined success: {outcome:?}"))),
                Err(error) => {
                    prop_assert_eq!(in_place, Err(error));
                    prop_assert!(borrowed.normalized_eq(&source));
                }
            }
        }
    }

    #[test]
    fn test_resolver_resolve_publication(
        source in publication_molecule(),
        rollback_method in prop::sample::select(vec![RollbackMethod::Clone, RollbackMethod::Journal]),
    ) {
        let model = ChemistryModel {
            valence: ValenceModel::smiles(),
            ..ChemistryModel::default()
        };
        let resolver = Resolver::with_config(&model, ResolveConfig {
            isotope: IsotopePolicy::Natural, ..Default::default()
        });
        let mut borrowed = source.clone();
        let consumed = resolver.resolve(source.clone());
        let in_place = resolver.resolve_into(&mut borrowed, rollback_method);
        match consumed {
            Ok(Solution::Determined(published)) => {
                prop_assert_eq!(in_place, Ok(Solution::Determined(())));
                prop_assert_eq!(&borrowed, &published);
                prop_assert_eq!(published.clone().edit().finish(), Ok(published));
            }
            Ok(Solution::Underdetermined(())) => {
                prop_assert_eq!(in_place, Ok(Solution::Underdetermined(())));
                prop_assert!(borrowed.normalized_eq(&source));
            }
            Ok(Solution::Contradictory(contradiction)) => {
                prop_assert_eq!(in_place, Ok(Solution::Contradictory(contradiction)));
                prop_assert!(borrowed.normalized_eq(&source));
            }
            Err(error) => {
                prop_assert_eq!(in_place, Err(error));
                prop_assert!(borrowed.normalized_eq(&source));
            }
        }
    }
}
