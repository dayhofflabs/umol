//! Validators over a resolved `MoleculeAst`, grouped by tier:
//! - tier 1 (integrity): entity-structure shape, constraint cross-checks (umol-ast);
//! - tier 2 (invariants): electron-count and spin physics;
//! - tier 3 (conformance): valence table / atom-typing, aromaticity, stereo (model-carrying).
//!
//! Each validator returns `Result<Solution<(), C>, E>`: `Determined` and
//! `Underdetermined` are both successful outcomes, only `Contradictory(C)` is a
//! `Solution`-side failure; setup-level failures live in `Err(E)`. The composite
//! [`Validator`] carries the chemistry model, runs the tiers in order (stereo
//! last), and lifts each engine's contradiction/error into unions via `From`.
//! `validate_atom` runs only the per-atom invariants (no surrounding molecule).

pub mod aromaticity;
pub mod invariant;
pub mod spin;
pub mod stereo;
pub mod valence;

pub use aromaticity::{AromaticityConformanceValidator, AromaticityValidatorContradiction};
pub use invariant::{ValenceInvariantsError, ValenceInvariantsValidator};
pub use spin::{SpinInvariantsContradiction, SpinInvariantsError, SpinInvariantsValidator};
pub use stereo::{
    StereoConformanceValidator, StereoValidateConfig, StereoValidatorContradiction,
    StereoValidatorError,
};
use thiserror::Error;
use umol_ast::ast::{
    AtomAst, ConstraintContradiction, ConstraintError, ConstraintValidateConfig,
    ConstraintValidator, EntityStructureContradiction, EntityStructureError,
    EntityStructureValidator, MoleculeAst, SubstructureMatchAlgorithm,
};
use umol_graph_core::{
    ConnectedComponentsAlgorithm, RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
};
use umol_utils::solution::Solution;
pub use valence::{
    ValenceConformanceContradiction, ValenceConformanceError, ValenceConformanceValidator,
};

use crate::ops::aromaticity::{AromaticityConfig, AromaticityError};
use crate::ops::invariant::ValenceMismatch;
use crate::ops::model::ChemistryModel;

/// Operational configuration for composite molecule validation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ValidateConfig {
    /// Algorithms used by model-independent constraint validation.
    pub constraint: ConstraintValidateConfig,
    /// Algorithms used by aromaticity conformance validation.
    pub aromaticity: AromaticityConfig,
    /// Algorithms and iteration limits used by stereo conformance validation.
    pub stereo: StereoValidateConfig,
}

impl Default for ValidateConfig {
    fn default() -> Self {
        Self {
            constraint: ConstraintValidateConfig {
                relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
                substructure_match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
                subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2Rdkit,
            },
            aromaticity: AromaticityConfig::default(),
            stereo: StereoValidateConfig::default(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Validator<'a> {
    // Integrity validators: MoleculeAst internal data consistency.
    pub entity_structure: EntityStructureValidator,
    pub constraint: ConstraintValidator,
    // Invariants validators: physical, model-independent invariants.
    pub valence_invariants: ValenceInvariantsValidator,
    pub spin_invariants: SpinInvariantsValidator,
    // Conformance validators: chemistry model compliance.
    pub valence_conformance: ValenceConformanceValidator<'a>,
    pub aromaticity: AromaticityConformanceValidator,
    pub stereo: StereoConformanceValidator,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValidatorContradiction {
    #[error(transparent)]
    ValenceInvariant(#[from] ValenceMismatch),
    #[error(transparent)]
    SpinInvariants(#[from] SpinInvariantsContradiction),
    #[error(transparent)]
    Constraint(#[from] ConstraintContradiction),
    #[error(transparent)]
    EntityStructure(#[from] EntityStructureContradiction),
    #[error(transparent)]
    ValenceConformance(#[from] ValenceConformanceContradiction),
    #[error(transparent)]
    Aromaticity(#[from] AromaticityValidatorContradiction),
    #[error(transparent)]
    Stereo(#[from] StereoValidatorContradiction),
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValidatorError {
    #[error(transparent)]
    ValenceInvariant(#[from] ValenceInvariantsError),
    #[error(transparent)]
    SpinInvariants(#[from] SpinInvariantsError),
    #[error(transparent)]
    Constraint(#[from] ConstraintError),
    #[error(transparent)]
    EntityStructure(#[from] EntityStructureError),
    #[error(transparent)]
    ValenceConformance(#[from] ValenceConformanceError),
    #[error(transparent)]
    Aromaticity(#[from] AromaticityError),
    #[error(transparent)]
    Stereo(#[from] StereoValidatorError),
}

impl<'a> Validator<'a> {
    pub fn new(model: &'a ChemistryModel) -> Self {
        Self::with_config(model, ValidateConfig::default())
    }

    pub fn with_config(model: &'a ChemistryModel, config: ValidateConfig) -> Self {
        Self {
            entity_structure: EntityStructureValidator,
            constraint: ConstraintValidator::new(config.constraint),
            valence_invariants: ValenceInvariantsValidator,
            spin_invariants: SpinInvariantsValidator,
            valence_conformance: ValenceConformanceValidator::new(&model.valence),
            aromaticity: AromaticityConformanceValidator::with_config(
                &model.aromaticity,
                config.aromaticity,
            ),
            stereo: StereoConformanceValidator::with_config(&model.stereo, config.stereo),
        }
    }

    /// Integrity: entity-structure shape and constraint cross-checks.
    pub fn validate_integrity(
        &self,
        ast: &MoleculeAst,
    ) -> Result<Solution<(), ValidatorContradiction>, ValidatorError> {
        let mut any_undetermined = false;
        match self.entity_structure.validate(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        match self.constraint.validate(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        Ok(verdict(any_undetermined))
    }

    /// Invariants: electron count and spin coupling.
    pub fn validate_invariants(
        &self,
        ast: &MoleculeAst,
    ) -> Result<Solution<(), ValidatorContradiction>, ValidatorError> {
        let mut any_undetermined = false;
        match self.valence_invariants.validate(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        match self.spin_invariants.validate(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        Ok(verdict(any_undetermined))
    }

    /// Conformance: chemistry-model conformance — valence table / atom-typing,
    /// aromaticity, and stereo.
    pub fn validate_conformance(
        &self,
        ast: &MoleculeAst,
    ) -> Result<Solution<(), ValidatorContradiction>, ValidatorError> {
        let mut any_undetermined = false;
        match self.valence_conformance.validate(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        match self.aromaticity.validate(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        match self.stereo.validate(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        Ok(verdict(any_undetermined))
    }

    /// All validators in order: integrity → invariants → conformance.
    pub fn validate(
        &self,
        ast: &MoleculeAst,
    ) -> Result<Solution<(), ValidatorContradiction>, ValidatorError> {
        let mut any_undetermined = false;
        for outcome in [
            self.validate_integrity(ast)?,
            self.validate_invariants(ast)?,
            self.validate_conformance(ast)?,
        ] {
            match outcome {
                Solution::Determined(()) => {}
                Solution::Underdetermined(()) => any_undetermined = true,
                Solution::Contradictory(c) => return Ok(Solution::Contradictory(c)),
            }
        }
        Ok(verdict(any_undetermined))
    }

    pub fn validate_atom(
        &self,
        atom: &AtomAst,
    ) -> Result<Solution<(), ValidatorContradiction>, ValidatorError> {
        let mut any_undetermined = false;
        match self.valence_invariants.validate_atom(atom)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        match self.spin_invariants.validate_atom(atom)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        Ok(verdict(any_undetermined))
    }
}

fn verdict(any_undetermined: bool) -> Solution<(), ValidatorContradiction> {
    if any_undetermined {
        Solution::Underdetermined(())
    } else {
        Solution::Determined(())
    }
}

#[cfg(test)]
mod tests {
    use proptest::prelude::*;
    use rstest::rstest;
    use umol_ast::ast::{
        AtomAst, AtomConstraintAst, AtomId, Constraint, ElementAst, MoleculeAst,
        MoleculeConstraint, MoleculeParts, RingConfig, UnpairedElectronsAst, ValueAst,
    };
    use umol_ast::{mol_dsl, mol_dsl_ground};
    use umol_chem::element::Element;
    use umol_chem::error::SpinStateError;
    use umol_chem::spin::SpinMultiplicity;
    use umol_graph_core::{
        AutomorphismAlgorithm, ConnectedComponentsAlgorithm, MaximumIndependentSetAlgorithm,
    };

    use super::*;
    use crate::ops::model::ChemistryModel;

    #[rstest]
    fn test_validate_config_default() {
        assert_eq!(
            ValidateConfig::default(),
            ValidateConfig {
                constraint: ConstraintValidateConfig {
                    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                    connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
                    substructure_match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
                    subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2Rdkit,
                },
                aromaticity: AromaticityConfig {
                    ring_config: RingConfig::default(),
                    connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
                    maximum_independent_set_algorithm:
                        MaximumIndependentSetAlgorithm::BranchAndBound,
                },
                stereo: StereoValidateConfig {
                    automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                    max_iterations: 16,
                },
            }
        );
    }

    #[rstest]
    #[case::ground(mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#))]
    #[case::non_ground(mol_dsl!(r#"{:atoms ["C"] :bonds []}"#))]
    fn test_validator_new(#[case] molecule: MoleculeAst) {
        let model = ChemistryModel::default();
        assert_eq!(
            Validator::new(&model).validate(&molecule),
            Validator::with_config(&model, ValidateConfig::default()).validate(&molecule)
        );
    }

    #[rstest]
    fn test_validator_with_config() {
        let molecule = mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#);
        let model = ChemistryModel::default();
        let config = ValidateConfig {
            constraint: ConstraintValidateConfig {
                relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
                substructure_match_algorithm: SubstructureMatchAlgorithm::GraphAndOverlays,
                subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Vf2Rdkit,
            },
            aromaticity: AromaticityConfig {
                ring_config: RingConfig::default(),
                connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
                maximum_independent_set_algorithm: MaximumIndependentSetAlgorithm::BranchAndBound,
            },
            stereo: StereoValidateConfig {
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                max_iterations: 8,
            },
        };
        let validator = Validator::with_config(&model, config);

        assert_eq!(
            validator.aromaticity.validate(&molecule),
            AromaticityConformanceValidator::with_config(&model.aromaticity, config.aromaticity)
                .validate(&molecule)
        );
        assert_eq!(
            validator.stereo.validate(&molecule),
            StereoConformanceValidator::with_config(&model.stereo, config.stereo)
                .validate(&molecule)
        );
    }

    #[rstest]
    #[case::ground(mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#), Solution::Determined(()))]
    #[case::non_ground(mol_dsl!(r#"{:atoms ["C"] :bonds []}"#), Solution::Underdetermined(()))]
    #[case::invalid_spin(
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst {
                element: ElementAst::Lit(Element::C),
                charge: ValueAst::Lit(0),
                implicit_hydrogens: ValueAst::Lit(2),
                lone_pairs: ValueAst::Lit(0),
                unpaired_electrons: UnpairedElectronsAst::from((2_u8, 2_u8)),
                ..Default::default()
            }],
            ..Default::default()
        }),
        Solution::Contradictory(ValidatorContradiction::SpinInvariants(
            SpinInvariantsContradiction::MoleculeAtom {
                atom: AtomId(0),
                error: SpinStateError::Incompatible {
                    unpaired_electrons: 2,
                    multiplicity: SpinMultiplicity::DOUBLET,
                },
            },
        )),
    )]
    #[case::valid_coupling_not_yet_evaluated(
        MoleculeAst::from_parts(MoleculeParts {
            constraints: Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                atoms: None,
                unpaired_electrons: UnpairedElectronsAst::from((2_u8, 3_u8)),
            }).into(),
            ..Default::default()
        }),
        Solution::Underdetermined(()),
    )]
    #[case::invalid_coupling(
        MoleculeAst::from_parts(MoleculeParts {
            constraints: Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                atoms: None,
                unpaired_electrons: UnpairedElectronsAst::from((2_u8, 2_u8)),
            }).into(),
            ..Default::default()
        }),
        Solution::Contradictory(ValidatorContradiction::SpinInvariants(
            SpinInvariantsContradiction::UnpairedElectronCoupling {
                constraint_index: 0,
                error: SpinStateError::Incompatible {
                    unpaired_electrons: 2,
                    multiplicity: SpinMultiplicity::DOUBLET,
                },
            },
        )),
    )]
    fn test_validator_validate(
        #[case] molecule: MoleculeAst,
        #[case] expected: Solution<(), ValidatorContradiction>,
    ) {
        let model = ChemistryModel::default();
        assert_eq!(
            Validator::new(&model).validate(&molecule).unwrap(),
            expected
        );
    }

    #[rstest]
    fn test_validator_validate_integrity() {
        let molecule = mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#);
        let model = ChemistryModel::default();
        assert_eq!(
            Validator::new(&model)
                .validate_integrity(&molecule)
                .unwrap(),
            Solution::Determined(())
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ground(
        mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#),
        Solution::Determined(()),
    )]
    #[case::partial_spin(
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst {
                element: ElementAst::Lit(Element::C),
                charge: ValueAst::Lit(0),
                implicit_hydrogens: ValueAst::Lit(4),
                lone_pairs: ValueAst::Lit(0),
                unpaired_electrons: UnpairedElectronsAst {
                    count: ValueAst::Lit(0),
                    multiplicity: ValueAst::Undetermined,
                },
                ..Default::default()
            }],
            ..Default::default()
        }),
        Solution::Underdetermined(()),
    )]
    #[case::invalid_spin(
        MoleculeAst::from_parts(MoleculeParts {
            atoms: vec![AtomAst {
                element: ElementAst::Lit(Element::C),
                charge: ValueAst::Lit(0),
                implicit_hydrogens: ValueAst::Lit(2),
                lone_pairs: ValueAst::Lit(0),
                unpaired_electrons: UnpairedElectronsAst::from((2_u8, 2_u8)),
                ..Default::default()
            }],
            ..Default::default()
        }),
        Solution::Contradictory(ValidatorContradiction::SpinInvariants(
            SpinInvariantsContradiction::MoleculeAtom {
                atom: AtomId(0),
                error: SpinStateError::Incompatible {
                    unpaired_electrons: 2,
                    multiplicity: SpinMultiplicity::DOUBLET,
                },
            },
        )),
    )]
    #[case::valid_coupling_not_yet_evaluated(
        MoleculeAst::from_parts(MoleculeParts {
            constraints: Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                atoms: None,
                unpaired_electrons: UnpairedElectronsAst::from((2_u8, 3_u8)),
            }).into(),
            ..Default::default()
        }),
        Solution::Underdetermined(()),
    )]
    #[case::invalid_coupling(
        MoleculeAst::from_parts(MoleculeParts {
            constraints: Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                atoms: None,
                unpaired_electrons: UnpairedElectronsAst::from((2_u8, 2_u8)),
            }).into(),
            ..Default::default()
        }),
        Solution::Contradictory(ValidatorContradiction::SpinInvariants(
            SpinInvariantsContradiction::UnpairedElectronCoupling {
                constraint_index: 0,
                error: SpinStateError::Incompatible {
                    unpaired_electrons: 2,
                    multiplicity: SpinMultiplicity::DOUBLET,
                },
            },
        )),
    )]
    fn test_validator_validate_invariants(
        #[case] molecule: MoleculeAst,
        #[case] expected: Solution<(), ValidatorContradiction>,
    ) {
        let model = ChemistryModel::default();
        assert_eq!(
            Validator::new(&model)
                .validate_invariants(&molecule)
                .unwrap(),
            expected,
        );
    }

    #[rstest]
    fn test_validator_validate_conformance() {
        let molecule = mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#);
        let model = ChemistryModel::default();
        assert_eq!(
            Validator::new(&model)
                .validate_conformance(&molecule)
                .unwrap(),
            Solution::Determined(())
        );
    }

    #[rstest]
    #[case::methane(4, None, UnpairedElectronsAst::from((0_u8, 1_u8)), Solution::Determined(()))]
    #[case::with_valence_constraint(3, Some(1), UnpairedElectronsAst::from((0_u8, 1_u8)), Solution::Determined(()))]
    #[case::partial_spin(
        4,
        None,
        UnpairedElectronsAst { count: ValueAst::Lit(0), multiplicity: ValueAst::Undetermined },
        Solution::Underdetermined(()),
    )]
    #[case::invalid_spin(
        2,
        None,
        UnpairedElectronsAst::from((2_u8, 2_u8)),
        Solution::Contradictory(ValidatorContradiction::SpinInvariants(
            SpinInvariantsContradiction::Atom {
                error: SpinStateError::Incompatible {
                    unpaired_electrons: 2,
                    multiplicity: SpinMultiplicity::DOUBLET,
                },
            },
        )),
    )]
    fn test_validator_validate_atom(
        #[case] hydrogens: i64,
        #[case] valence: Option<i64>,
        #[case] unpaired_electrons: UnpairedElectronsAst,
        #[case] expected: Solution<(), ValidatorContradiction>,
    ) {
        let mut atom = AtomAst::from_element(Element::C);
        atom.charge = ValueAst::Lit(0);
        atom.lone_pairs = ValueAst::Lit(0);
        atom.implicit_hydrogens = ValueAst::Lit(hydrogens);
        atom.unpaired_electrons = unpaired_electrons;
        if let Some(v) = valence {
            atom.constraints
                .set(AtomConstraintAst::Valence(ValueAst::Lit(v)));
        }
        let model = ChemistryModel::default();
        assert_eq!(
            Validator::new(&model).validate_atom(&atom).unwrap(),
            expected,
        );
    }

    proptest! {
        #[test]
        fn test_validator_validate_atom_spin(
            count in 0_u8..5,
            multiplicity in 0_u8..8,
        ) {
            let atom = AtomAst {
                element: ElementAst::Lit(Element::C),
                charge: ValueAst::Lit(0),
                implicit_hydrogens: ValueAst::Lit(4 - i64::from(count)),
                lone_pairs: ValueAst::Lit(0),
                unpaired_electrons: UnpairedElectronsAst::from((count, multiplicity)),
                ..Default::default()
            };
            let expected = SpinInvariantsValidator
                .validate_atom(&atom)
                .unwrap()
                .map_contradiction(ValidatorContradiction::SpinInvariants);
            let model = ChemistryModel::default();

            prop_assert_eq!(Validator::new(&model).validate_atom(&atom).unwrap(), expected);
        }

        #[test]
        fn test_validator_validate_invariants_spin(
            count in 0_u8..5,
            multiplicity in 0_u8..8,
        ) {
            let molecule = MoleculeAst::from_parts(MoleculeParts {
                atoms: vec![AtomAst {
                    element: ElementAst::Lit(Element::C),
                    charge: ValueAst::Lit(0),
                    implicit_hydrogens: ValueAst::Lit(4 - i64::from(count)),
                    lone_pairs: ValueAst::Lit(0),
                    unpaired_electrons: UnpairedElectronsAst::from((count, multiplicity)),
                    ..Default::default()
                }],
                ..Default::default()
            });
            let expected = SpinInvariantsValidator
                .validate(&molecule)
                .unwrap()
                .map_contradiction(ValidatorContradiction::SpinInvariants);
            let model = ChemistryModel::default();

            prop_assert_eq!(
                Validator::new(&model).validate_invariants(&molecule).unwrap(),
                expected,
            );
        }
    }
}
