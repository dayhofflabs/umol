//! Validators over a resolved `Molecule`, grouped by tier:
//! - tier 2 (invariants): constraint satisfaction, electron count, and spin;
//! - tier 3 (conformance): valence table / atom-typing, aromaticity, stereo (model-carrying).
//!
//! Each validator returns `Result<Solution<(), C>, E>`: `Determined` and
//! `Underdetermined` are both successful outcomes, only `Contradictory(C)` is a
//! `Solution`-side failure; setup-level failures live in `Err(E)`. The composite
//! [`Validator`] carries the chemistry model, runs the tiers in order (stereo
//! last), and lifts each engine's contradiction/error into unions via `From`.

pub mod aromaticity;
pub mod connectivity;
pub mod constraint;
pub mod spin;
pub mod stereo;
pub mod valence;

pub use aromaticity::{AromaticityConformanceContradiction, AromaticityConformanceValidator};
pub use connectivity::{
    ConnectivityConformanceContradiction, ConnectivityConformanceError,
    ConnectivityConformanceValidator, ConnectivityModel,
};
pub use constraint::{
    ConstraintInvariantsContradiction, ConstraintInvariantsError, ConstraintInvariantsValidator,
    ConstraintValidateConfig, DerivedKind, IncidenceConstraintInvariantsContradiction,
    IncidenceConstraintInvariantsValidator, MoleculeConstraintInvariantsContradiction,
    MoleculeConstraintInvariantsValidator, RelationalConstraintInvariantsContradiction,
    RelationalConstraintInvariantsValidator, RingConstraintInvariantsContradiction,
    RingConstraintInvariantsValidator,
};
pub use spin::{SpinInvariantsContradiction, SpinInvariantsError, SpinInvariantsValidator};
pub use stereo::{
    StereoConformanceContradiction, StereoConformanceError, StereoConformanceValidator,
    StereoValidateConfig,
};
use thiserror::Error;
use umol_graph_core::{ConnectedComponentsAlgorithm, RelevantCycleEnumerationAlgorithm};
use umol_graph_ir::ir::Molecule;
use umol_utils::solution::Solution;
pub use valence::{
    ValenceConformanceContradiction, ValenceConformanceError, ValenceConformanceValidator,
};

use crate::ops::aromaticity::{AromaticityConfig, AromaticityError};
use crate::ops::invariant::ValenceMismatch;
pub use crate::ops::invariant::{ValenceInvariantsError, ValenceInvariantsValidator};
use crate::ops::model::ChemistryModel;

/// Operational configuration for composite molecule validation.
///
/// The model-independent defaults are relevant cycles via Vismara and connected components via
/// BFS.
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
                derived_kind: DerivedKind::DerivedComplete,
            },
            aromaticity: AromaticityConfig::default(),
            stereo: StereoValidateConfig::default(),
        }
    }
}

/// Composite tier-2 invariant and tier-3 chemistry-model validator.
#[derive(Clone, Debug)]
pub struct Validator<'a> {
    // Invariants validators: model-independent constraint and physical semantics.
    pub constraint: ConstraintInvariantsValidator,
    pub valence_invariants: ValenceInvariantsValidator,
    pub spin_invariants: SpinInvariantsValidator,
    // Conformance validators: chemistry model compliance.
    pub connectivity: ConnectivityConformanceValidator<'a>,
    pub valence_conformance: ValenceConformanceValidator<'a>,
    pub aromaticity: AromaticityConformanceValidator,
    pub stereo: StereoConformanceValidator,
}

/// Semantic contradiction returned by one component of [`Validator`].
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValidateContradiction {
    #[error(transparent)]
    ValenceInvariant(#[from] ValenceMismatch),
    #[error(transparent)]
    SpinInvariants(#[from] SpinInvariantsContradiction),
    #[error(transparent)]
    Constraint(#[from] ConstraintInvariantsContradiction),
    #[error(transparent)]
    Connectivity(#[from] ConnectivityConformanceContradiction),
    #[error(transparent)]
    ValenceConformance(#[from] ValenceConformanceContradiction),
    #[error(transparent)]
    Aromaticity(#[from] AromaticityConformanceContradiction),
    #[error(transparent)]
    Stereo(#[from] StereoConformanceContradiction),
}

/// Setup, reference, or unsupported-operation failure returned by one validator component.
#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValidateError {
    #[error(transparent)]
    ValenceInvariant(#[from] ValenceInvariantsError),
    #[error(transparent)]
    SpinInvariants(#[from] SpinInvariantsError),
    #[error(transparent)]
    Constraint(#[from] ConstraintInvariantsError),
    #[error(transparent)]
    Connectivity(#[from] ConnectivityConformanceError),
    #[error(transparent)]
    ValenceConformance(#[from] ValenceConformanceError),
    #[error(transparent)]
    Aromaticity(#[from] AromaticityError),
    #[error(transparent)]
    Stereo(#[from] StereoConformanceError),
}

impl<'a> Validator<'a> {
    /// Construct a validator with the operation defaults in [`ValidateConfig`].
    pub fn new(model: &'a ChemistryModel) -> Self {
        Self::with_config(model, ValidateConfig::default())
    }

    /// Construct a validator with explicit operational configuration.
    pub fn with_config(model: &'a ChemistryModel, config: ValidateConfig) -> Self {
        Self {
            constraint: ConstraintInvariantsValidator::new(config.constraint),
            valence_invariants: ValenceInvariantsValidator,
            spin_invariants: SpinInvariantsValidator,
            connectivity: ConnectivityConformanceValidator::new(&model.connectivity),
            valence_conformance: ValenceConformanceValidator::new(&model.valence),
            aromaticity: AromaticityConformanceValidator::with_config(
                &model.aromaticity,
                config.aromaticity,
            ),
            stereo: StereoConformanceValidator::with_config(&model.stereo, config.stereo),
        }
    }

    /// Model-independent constraint, electron-count, and spin invariants.
    pub fn validate_invariants(
        &self,
        molecule: &Molecule,
    ) -> Result<Solution<(), ValidateContradiction>, ValidateError> {
        let mut any_undetermined = false;
        match self.constraint.validate(molecule)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        match self.valence_invariants.validate(molecule)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        match self.spin_invariants.validate(molecule)? {
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
        molecule: &Molecule,
    ) -> Result<Solution<(), ValidateContradiction>, ValidateError> {
        let mut any_undetermined = false;
        match self.connectivity.validate(molecule)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        match self.valence_conformance.validate(molecule)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        match self.aromaticity.validate(molecule)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        match self.stereo.validate(molecule)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        Ok(verdict(any_undetermined))
    }

    /// All semantic validators in order: invariants, then model conformance.
    pub fn validate(
        &self,
        molecule: &Molecule,
    ) -> Result<Solution<(), ValidateContradiction>, ValidateError> {
        let mut any_undetermined = false;
        match self.validate_invariants(molecule)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c)),
        }
        match self.validate_conformance(molecule)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c)),
        }
        Ok(verdict(any_undetermined))
    }
}

fn verdict(any_undetermined: bool) -> Solution<(), ValidateContradiction> {
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
    use umol_chem::element::Element;
    use umol_chem::error::SpinStateError;
    use umol_chem::spin::SpinMultiplicity;
    use umol_graph_core::{
        AutomorphismAlgorithm, ConnectedComponentsAlgorithm, MaximumIndependentSetAlgorithm,
    };
    use umol_graph_ir::ir::{
        AtomConstraintForm, AtomForm, AtomId, Constraint, DativeBondId, ElementForm, Entity,
        Molecule, MoleculeConstraint, MoleculeEntries, NumForm, RelationalConstraint, RingConfig,
        RingScope, StereoAtomConstraintForm, StereoAtomId, StereoBondConstraintForm, StereoBondId,
        StereoKind, StereogenicityForm, UnpairedElectronsForm,
    };
    use umol_graph_ir::{mol_dsl, mol_dsl_ground};

    use super::*;
    use crate::ops::aromaticity::AromaticityInconsistency;
    use crate::ops::model::ChemistryModel;

    #[rstest]
    fn test_validate_config_default() {
        assert_eq!(
            ValidateConfig::default(),
            ValidateConfig {
                constraint: ConstraintValidateConfig {
                    relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm::Vismara,
                    connected_components_algorithm: ConnectedComponentsAlgorithm::Bfs,
                    derived_kind: DerivedKind::DerivedComplete,
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
    fn test_validator_new(#[case] molecule: Molecule) {
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
                derived_kind: DerivedKind::DerivedComplete,
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
            validator.constraint,
            ConstraintInvariantsValidator::new(config.constraint)
        );
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
    #[case::incidence(ConstraintInvariantsContradiction::Incidence(
        IncidenceConstraintInvariantsContradiction::Atom {
            atom: AtomId(0),
            constraint: AtomConstraintForm::valence(1),
        }
    ))]
    #[case::ring(ConstraintInvariantsContradiction::Ring(RingConstraintInvariantsContradiction::Atom {
        atom: AtomId(0),
        constraint: AtomConstraintForm::ring_membership(RingScope::All, 1),
    }))]
    #[case::relational(ConstraintInvariantsContradiction::Relational(
        RelationalConstraintInvariantsContradiction {
            constraint: RelationalConstraint::DativeBondDonor {
                bond: DativeBondId(0),
                atom: AtomId(0),
            },
        }
    ))]
    #[case::molecule(ConstraintInvariantsContradiction::Molecule(
        MoleculeConstraintInvariantsContradiction {
            constraint: MoleculeConstraint::Connected { atoms: None },
        }
    ))]
    #[case::logical(ConstraintInvariantsContradiction::Logical {
        constraint: Constraint::And(Vec::new()),
    })]
    #[case::stereo_atom(ConstraintInvariantsContradiction::StereoAtom {
        id: StereoAtomId(0),
        kind: StereoKind::Tetrahedral,
        constraint: StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
    })]
    #[case::stereo_bond(ConstraintInvariantsContradiction::StereoBond {
        id: StereoBondId(0),
        kind: StereoKind::CisTrans,
        constraint: StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Undetermined),
    })]
    fn test_validator_contradiction_from(#[case] input: ConstraintInvariantsContradiction) {
        assert_eq!(
            ValidateContradiction::from(input.clone()),
            ValidateContradiction::Constraint(input)
        );
    }

    #[rstest]
    #[case::invalid_reference(ConstraintInvariantsError::InvalidReference {
        entity: Entity::Atom(AtomId(0)),
    })]
    #[case::dative_ring(ConstraintInvariantsError::DativeBondRingMembershipUnsupported {
        bond: DativeBondId(0),
    })]
    fn test_validator_error_from(#[case] input: ConstraintInvariantsError) {
        assert_eq!(
            ValidateError::from(input.clone()),
            ValidateError::Constraint(input)
        );
    }

    #[rstest]
    #[case::ground(mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#), Solution::Determined(()))]
    #[case::non_ground(mol_dsl!(r#"{:atoms ["C"] :bonds []}"#), Solution::Underdetermined(()))]
    #[case::invalid_spin(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm {
                element: ElementForm::Lit(Element::C),
                charge: NumForm::Lit(0),
                implicit_hydrogens: NumForm::Lit(2),
                lone_pairs: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm::from((2_u8, 2_u8)),
                ..Default::default()
            }],
            ..Default::default()
        }),
        Solution::Contradictory(ValidateContradiction::SpinInvariants(
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
        Molecule::from_entries(MoleculeEntries {
            constraints: Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                atoms: None,
                unpaired_electrons: UnpairedElectronsForm::from((2_u8, 3_u8)),
            }).into(),
            ..Default::default()
        }),
        Solution::Underdetermined(()),
    )]
    #[case::invalid_coupling(
        Molecule::from_entries(MoleculeEntries {
            constraints: Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                atoms: None,
                unpaired_electrons: UnpairedElectronsForm::from((2_u8, 2_u8)),
            }).into(),
            ..Default::default()
        }),
        Solution::Contradictory(ValidateContradiction::SpinInvariants(
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
        #[case] molecule: Molecule,
        #[case] expected: Solution<(), ValidateContradiction>,
    ) {
        let model = ChemistryModel::default();
        assert_eq!(
            Validator::new(&model).validate(&molecule).unwrap(),
            expected
        );
    }

    #[rstest]
    fn test_validator_validate_invariants_error() {
        let mut molecule = Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm::from_element(Element::C)],
            ..Default::default()
        });
        molecule
            .constraints_mut()
            .push(Constraint::Atom(AtomId(1), AtomConstraintForm::valence(0)));
        let model = ChemistryModel::default();

        assert_eq!(
            Validator::new(&model).validate_invariants(&molecule),
            Err(ValidateError::Constraint(
                ConstraintInvariantsError::InvalidReference {
                    entity: Entity::Atom(AtomId(1)),
                }
            ))
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ground(
        mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#),
        Solution::Determined(()),
    )]
    #[case::partial_spin(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm {
                element: ElementForm::Lit(Element::C),
                charge: NumForm::Lit(0),
                implicit_hydrogens: NumForm::Lit(4),
                lone_pairs: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm {
                    count: NumForm::Lit(0),
                    multiplicity: NumForm::Undetermined,
                },
                ..Default::default()
            }],
            ..Default::default()
        }),
        Solution::Underdetermined(()),
    )]
    #[case::invalid_spin(
        Molecule::from_entries(MoleculeEntries {
            atoms: vec![AtomForm {
                element: ElementForm::Lit(Element::C),
                charge: NumForm::Lit(0),
                implicit_hydrogens: NumForm::Lit(2),
                lone_pairs: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm::from((2_u8, 2_u8)),
                ..Default::default()
            }],
            ..Default::default()
        }),
        Solution::Contradictory(ValidateContradiction::SpinInvariants(
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
        Molecule::from_entries(MoleculeEntries {
            constraints: Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                atoms: None,
                unpaired_electrons: UnpairedElectronsForm::from((2_u8, 3_u8)),
            }).into(),
            ..Default::default()
        }),
        Solution::Underdetermined(()),
    )]
    #[case::invalid_coupling(
        Molecule::from_entries(MoleculeEntries {
            constraints: Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                atoms: None,
                unpaired_electrons: UnpairedElectronsForm::from((2_u8, 2_u8)),
            }).into(),
            ..Default::default()
        }),
        Solution::Contradictory(ValidateContradiction::SpinInvariants(
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
        #[case] molecule: Molecule,
        #[case] expected: Solution<(), ValidateContradiction>,
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
    #[case::ground(
        mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#),
        Solution::Determined(()),
    )]
    #[case::aromaticity(
        mol_dsl!(r#"{:atoms ["C#a"]}"#),
        Solution::Contradictory(ValidateContradiction::Aromaticity(
            AromaticityConformanceContradiction::Inconsistency(
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) },
            ),
        )),
    )]
    #[case::stereo(
        mol_dsl_ground!(r#"{
            :atoms ["C#h0#T1" "C#h3" "C#h3" "C#h3" "C#h3"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        }"#),
        Solution::Determined(()),
    )]
    #[case::aromaticity_partial(
        mol_dsl!(r#"{
            :atoms ["C#a+" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        }"#),
        Solution::Underdetermined(()),
    )]
    fn test_validator_validate_conformance(
        #[case] molecule: Molecule,
        #[case] expected: Solution<(), ValidateContradiction>,
    ) {
        let original = molecule.clone();
        let model = ChemistryModel::default();

        assert_eq!(
            Validator::new(&model)
                .validate_conformance(&molecule)
                .unwrap(),
            expected
        );
        assert_eq!(molecule, original);
    }

    proptest! {
        #[test]
        fn test_validator_validate_invariants_spin(
            count in 0_u8..5,
            multiplicity in 0_u8..8,
        ) {
            let molecule = Molecule::from_entries(MoleculeEntries {
                atoms: vec![AtomForm {
                    element: ElementForm::Lit(Element::C),
                    charge: NumForm::Lit(0),
                    implicit_hydrogens: NumForm::Lit(4 - i64::from(count)),
                    lone_pairs: NumForm::Lit(0),
                    unpaired_electrons: UnpairedElectronsForm::from((count, multiplicity)),
                    ..Default::default()
                }],
                ..Default::default()
            });
            let expected = SpinInvariantsValidator
                .validate(&molecule)
                .unwrap()
                .map_contradiction(ValidateContradiction::SpinInvariants);
            let model = ChemistryModel::default();

            prop_assert_eq!(
                Validator::new(&model).validate_invariants(&molecule).unwrap(),
                expected,
            );
        }
    }
}
