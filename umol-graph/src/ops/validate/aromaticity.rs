//! Aromaticity validator. Wraps [`AromaticityPerception`] and verifies that
//! aromatic constraints and systems agree with the selected model.

use thiserror::Error;
use umol_graph_ir::ir::MoleculeAst;
use umol_utils::solution::Solution;

use crate::ops::aromaticity::{
    AromaticityConfig, AromaticityContradiction, AromaticityError, AromaticityInconsistency,
    AromaticityPerception,
};
use crate::ops::model::AromaticityModel;

#[derive(Clone, Debug)]
pub struct AromaticityConformanceValidator {
    perception: AromaticityPerception,
    config: AromaticityConfig,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum AromaticityValidatorContradiction {
    #[error("perception rejected the input: {0}")]
    Perception(AromaticityContradiction),
    #[error(transparent)]
    Inconsistency(#[from] AromaticityInconsistency),
}

impl AromaticityConformanceValidator {
    pub fn new(model: &AromaticityModel) -> Self {
        Self::with_config(model, AromaticityConfig::default())
    }

    pub fn with_config(model: &AromaticityModel, config: AromaticityConfig) -> Self {
        Self {
            perception: AromaticityPerception::new(model),
            config,
        }
    }

    pub fn validate(
        &self,
        ast: &MoleculeAst,
    ) -> Result<Solution<(), AromaticityValidatorContradiction>, AromaticityError> {
        match self.perception.derive(ast, self.config)? {
            Solution::Determined(derivation) => {
                if let Some(&inconsistency) = derivation.inconsistencies.first() {
                    Ok(Solution::Contradictory(inconsistency.into()))
                } else {
                    Ok(Solution::Determined(()))
                }
            }
            Solution::Underdetermined(_) => Ok(Solution::Underdetermined(())),
            Solution::Contradictory(contradiction) => Ok(Solution::Contradictory(
                AromaticityValidatorContradiction::Perception(contradiction),
            )),
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_core::{
        ConnectedComponentsAlgorithm, MaximumIndependentSetAlgorithm,
        RelevantCycleEnumerationAlgorithm, SimpleCycleEnumerationAlgorithm,
    };
    use umol_graph_ir::ir::{AromaticSystemId, AtomId, BondId, RingConfig};
    use umol_graph_ir::mol_dsl;

    use super::*;
    use crate::ops::model::{ElementScope, RingLimits};

    #[rstest]
    fn test_aromaticity_conformance_validator_new() {
        let model = AromaticityModel::daylight();
        let molecule = mol_dsl!(
            r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"]
                    [4 5 "1#a"] [5 0 "1#a"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]"}]
        }"#
        );

        assert_eq!(
            AromaticityConformanceValidator::new(&model).validate(&molecule),
            AromaticityConformanceValidator::with_config(&model, AromaticityConfig::default())
                .validate(&molecule)
        );
    }

    #[rstest]
    fn test_aromaticity_conformance_validator_with_config() {
        let model = AromaticityModel::daylight();
        let molecule = mol_dsl!(
            r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"]
                    [4 5 "1#a"] [5 0 "1#a"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]"}]
        }"#
        );

        assert_eq!(
            AromaticityConformanceValidator::with_config(
                &model,
                AromaticityConfig {
                    ring_config: RingConfig {
                        simple_cycle_algorithm: SimpleCycleEnumerationAlgorithm::ReadTarjan,
                        relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                    },
                    connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
                    maximum_independent_set_algorithm:
                        MaximumIndependentSetAlgorithm::BranchAndBound,
                },
            )
            .validate(&molecule),
            Ok(Solution::Determined(()))
        );
    }

    #[rstest]
    #[case::assertion_without_system(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{:atoms ["C#a"]}"#),
        Solution::Contradictory(AromaticityValidatorContradiction::Inconsistency(
            AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) },
        )),
    )]
    #[case::participant_set_mismatch(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :type "[1,1,1,1,1]"}]
        }"#),
        Solution::Contradictory(AromaticityValidatorContradiction::Inconsistency(
            AromaticityInconsistency::AromaticSystemFailure {
                system: AromaticSystemId(0),
            },
        )),
    )]
    #[case::contribution_mismatch(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[2,1,1,1,1,1]"}]
        }"#),
        Solution::Contradictory(AromaticityValidatorContradiction::Inconsistency(
            AromaticityInconsistency::AromaticSystemFailure {
                system: AromaticSystemId(0),
            },
        )),
    )]
    #[case::bond_constraint_mismatch(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1#a!"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"]
                    [4 5 "1#a"] [5 0 "1#a"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]"}]
        }"#),
        Solution::Contradictory(AromaticityValidatorContradiction::Inconsistency(
            AromaticityInconsistency::AromaticBondConstraintMismatch {
                bond: BondId(0),
                system: AromaticSystemId(0),
            },
        )),
    )]
    #[case::aromatic_valence_mismatch(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"]
                    [4 5 "1#a"] [5 0 "1#a"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[2,0,1,1,1,1]"}]
        }"#),
        Solution::Contradictory(AromaticityValidatorContradiction::Inconsistency(
            AromaticityInconsistency::AromaticValenceMismatch {
                atom: AtomId(0),
                system: AromaticSystemId(0),
            },
        )),
    )]
    #[case::model_rejection(
        AromaticityModel::mdl(),
        mol_dsl!(r#"{
            :atoms ["O#n1" "C#h" "C#h" "C#h" "C#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4] :type "[2,1,1,1,1]"}]
        }"#),
        Solution::Contradictory(AromaticityValidatorContradiction::Inconsistency(
            AromaticityInconsistency::AromaticSystemFailure {
                system: AromaticSystemId(0),
            },
        )),
    )]
    #[case::perception_rejection(
        AromaticityModel::Clar {
            scope: ElementScope::Any,
            ring_limits: RingLimits::default(),
        },
        mol_dsl!(r#"{
            :atoms ["N#h#a2" "C#h#a" "C#h#a" "C#h#a" "C#h#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 0 "1"]]
        }"#),
        Solution::Contradictory(AromaticityValidatorContradiction::Perception(
            AromaticityContradiction::ClarNonBenzenoid(
                "Clar model requires benzenoid input but non-carbon aromatic atoms are present"
                    .to_string(),
            ),
        )),
    )]
    #[case::absent(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{
            :atoms ["C" "C" "C" "C" "C" "C"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]"}]
        }"#),
        Solution::Determined(()),
    )]
    #[case::vacuous(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{
            :atoms ["C#a*" "C#a*" "C#a*" "C#a*" "C#a*" "C#a*"]
            :bonds [[0 1 "1#a*"] [1 2 "1#a*"] [2 3 "1#a*"] [3 4 "1#a*"]
                    [4 5 "1#a*"] [5 0 "1#a*"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]"}]
        }"#),
        Solution::Determined(()),
    )]
    #[case::conformant(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{
            :atoms ["C#a" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1#a"] [1 2 "1#a"] [2 3 "1#a"] [3 4 "1#a"]
                    [4 5 "1#a"] [5 0 "1#a"]]
            :aromatic-systems [{:atoms [0 1 2 3 4 5] :type "[1,1,1,1,1,1]"}]
        }"#),
        Solution::Determined(()),
    )]
    #[case::non_ground(
        AromaticityModel::daylight(),
        mol_dsl!(r#"{
            :atoms ["C#a+" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        }"#),
        Solution::Underdetermined(()),
    )]
    fn test_aromaticity_conformance_validator_validate(
        #[case] model: AromaticityModel,
        #[case] molecule: MoleculeAst,
        #[case] expected: Solution<(), AromaticityValidatorContradiction>,
    ) {
        assert_eq!(
            AromaticityConformanceValidator::new(&model).validate(&molecule),
            Ok(expected)
        );
    }

    #[rstest]
    fn test_aromaticity_conformance_validator_validate_error() {
        let model = AromaticityModel::Hmo {
            scope: ElementScope::Any,
            stabilization_threshold: 0.0,
        };
        let molecule = mol_dsl!(
            r#"{
            :atoms ["C#a2" "C#a2" "C#a2" "C#a2" "C#a2" "C#a2"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        }"#
        );

        assert_eq!(
            AromaticityConformanceValidator::new(&model).validate(&molecule),
            Err(AromaticityError::HmoMissingParameters(
                "no Van-Catledge parameters for C with 2 pi-electrons".to_string(),
            ))
        );
    }
}
