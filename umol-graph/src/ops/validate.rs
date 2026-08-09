//! Validators over a resolved `MoleculeAst`, grouped by tier:
//! - tier 1 (integrity): entity-structure shape, constraint cross-checks (umol-graph-ir);
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
use umol_graph_core::{
    ConnectedComponentsAlgorithm, RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
};
use umol_graph_ir::ir::{
    AtomAst, ConstraintContradiction, ConstraintError, ConstraintValidateConfig,
    ConstraintValidator, EntityStructureContradiction, EntityStructureError,
    EntityStructureValidator, MoleculeAst, SubstructureMatchAlgorithm,
};
use umol_utils::solution::Solution;
pub use valence::{
    ValenceConformanceContradiction, ValenceConformanceError, ValenceConformanceValidator,
};

use crate::ops::aromaticity::{AromaticityConfig, AromaticityError};
use crate::ops::invariant::ValenceMismatch;
use crate::ops::model::ChemistryModel;

/// Operational configuration for composite molecule validation.
///
/// The model-independent defaults are relevant cycles via Vismara, connected components via BFS,
/// graph-and-overlays substructure matching, and VF2-RDKit subgraph isomorphism.
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
        match self.validate_integrity(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c)),
        }
        match self.validate_invariants(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c)),
        }
        match self.validate_conformance(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c)),
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
    use umol_chem::element::Element;
    use umol_chem::error::SpinStateError;
    use umol_chem::spin::SpinMultiplicity;
    use umol_graph_core::{
        AutomorphismAlgorithm, ConnectedComponentsAlgorithm, MaximumIndependentSetAlgorithm,
    };
    use umol_graph_ir::ir::{
        AtomAst, AtomConstraintAst, AtomId, Constraint, DativeBondId, ElementForm, Entity,
        IncidenceConstraintContradiction, MoleculeAst, MoleculeConstraint,
        MoleculeConstraintContradiction, MoleculeEntries, NumForm, RelationalConstraint,
        RelationalConstraintContradiction, RingConfig, RingConstraintContradiction, RingScope,
        StereoAtomConstraintAst, StereoAtomId, StereoBondConstraintAst, StereoBondId, StereoKind,
        StereogenicityAst, UnpairedElectronsForm,
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
                substructure_match_algorithm: SubstructureMatchAlgorithm::Incidence,
                subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm::Ullmann,
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
            ConstraintValidator::new(config.constraint)
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
    #[case::incidence(ConstraintContradiction::Incidence(
        IncidenceConstraintContradiction::Atom {
            atom: AtomId(0),
            constraint: AtomConstraintAst::valence(1),
        }
    ))]
    #[case::ring(ConstraintContradiction::Ring(RingConstraintContradiction::Atom {
        atom: AtomId(0),
        constraint: AtomConstraintAst::ring_membership(RingScope::All, 1),
    }))]
    #[case::relational(ConstraintContradiction::Relational(
        RelationalConstraintContradiction {
            constraint: RelationalConstraint::DativeBondDonor {
                bond: DativeBondId(0),
                atom: AtomId(0),
            },
        }
    ))]
    #[case::molecule(ConstraintContradiction::Molecule(
        MoleculeConstraintContradiction {
            constraint: MoleculeConstraint::Connected { atoms: None },
        }
    ))]
    #[case::logical(ConstraintContradiction::Logical {
        constraint: Constraint::And(Vec::new()),
    })]
    #[case::stereo_atom(ConstraintContradiction::StereoAtom {
        id: StereoAtomId(0),
        kind: StereoKind::Tetrahedral,
        constraint: StereoAtomConstraintAst::Stereogenicity(StereogenicityAst::Undetermined),
    })]
    #[case::stereo_bond(ConstraintContradiction::StereoBond {
        id: StereoBondId(0),
        kind: StereoKind::CisTrans,
        constraint: StereoBondConstraintAst::Stereogenicity(StereogenicityAst::Undetermined),
    })]
    fn test_validator_contradiction_from(#[case] input: ConstraintContradiction) {
        assert_eq!(
            ValidatorContradiction::from(input.clone()),
            ValidatorContradiction::Constraint(input)
        );
    }

    #[rstest]
    #[case::invalid_reference(ConstraintError::InvalidReference {
        entity: Entity::Atom(AtomId(0)),
    })]
    #[case::dative_ring(ConstraintError::DativeBondRingMembershipUnsupported {
        bond: DativeBondId(0),
    })]
    fn test_validator_error_from(#[case] input: ConstraintError) {
        assert_eq!(
            ValidatorError::from(input.clone()),
            ValidatorError::Constraint(input)
        );
    }

    #[rstest]
    #[case::ground(mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#), Solution::Determined(()))]
    #[case::non_ground(mol_dsl!(r#"{:atoms ["C"] :bonds []}"#), Solution::Underdetermined(()))]
    #[case::invalid_spin(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst {
                element: ElementForm::Lit(Element::C),
                charge: NumForm::Lit(0),
                implicit_hydrogens: NumForm::Lit(2),
                lone_pairs: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm::from((2_u8, 2_u8)),
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
        MoleculeAst::from_entries(MoleculeEntries {
            constraints: Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                atoms: None,
                unpaired_electrons: UnpairedElectronsForm::from((2_u8, 3_u8)),
            }).into(),
            ..Default::default()
        }),
        Solution::Underdetermined(()),
    )]
    #[case::invalid_coupling(
        MoleculeAst::from_entries(MoleculeEntries {
            constraints: Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                atoms: None,
                unpaired_electrons: UnpairedElectronsForm::from((2_u8, 2_u8)),
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
    #[case::entity_structure(
        mol_dsl!(r#"{:atoms ["C"] :bonds [[0 0 "1"]]}"#),
        Solution::Contradictory(ValidatorContradiction::EntityStructure(
            EntityStructureContradiction::BondSelfLoop { atom: AtomId(0) },
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
    #[case::empty(
        MoleculeAst::default(),
        Ok(Solution::Determined(())),
    )]
    #[case::without_constraints(
        mol_dsl!(r#"{:atoms ["C"] :bonds []}"#),
        Ok(Solution::Determined(())),
    )]
    #[case::contradiction(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst::from_element(Element::C)],
            constraints: Constraint::Atom(AtomId(0), AtomConstraintAst::valence(1)).into(),
            ..Default::default()
        }),
        Ok(Solution::Contradictory(ValidatorContradiction::Constraint(
            ConstraintContradiction::Incidence(IncidenceConstraintContradiction::Atom {
                atom: AtomId(0),
                constraint: AtomConstraintAst::valence(1),
            }),
        ))),
    )]
    #[case::error(
        {
            let mut molecule = MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![AtomAst::from_element(Element::C)],
                ..Default::default()
            });
            molecule
                .constraints_mut()
                .push(Constraint::Atom(AtomId(1), AtomConstraintAst::valence(0)));
            molecule
        },
        Err(ValidatorError::Constraint(ConstraintError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        })),
    )]
    fn test_validator_validate_integrity(
        #[case] molecule: MoleculeAst,
        #[case] expected: Result<Solution<(), ValidatorContradiction>, ValidatorError>,
    ) {
        let model = ChemistryModel::default();
        assert_eq!(
            Validator::new(&model).validate_integrity(&molecule),
            expected
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::ground(
        mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#),
        Solution::Determined(()),
    )]
    #[case::partial_spin(
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst {
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
        MoleculeAst::from_entries(MoleculeEntries {
            atoms: vec![AtomAst {
                element: ElementForm::Lit(Element::C),
                charge: NumForm::Lit(0),
                implicit_hydrogens: NumForm::Lit(2),
                lone_pairs: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm::from((2_u8, 2_u8)),
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
        MoleculeAst::from_entries(MoleculeEntries {
            constraints: Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                atoms: None,
                unpaired_electrons: UnpairedElectronsForm::from((2_u8, 3_u8)),
            }).into(),
            ..Default::default()
        }),
        Solution::Underdetermined(()),
    )]
    #[case::invalid_coupling(
        MoleculeAst::from_entries(MoleculeEntries {
            constraints: Constraint::Molecule(MoleculeConstraint::UnpairedElectronCoupling {
                atoms: None,
                unpaired_electrons: UnpairedElectronsForm::from((2_u8, 2_u8)),
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
    #[case::ground(
        mol_dsl_ground!(r#"{:atoms ["C #h4"] :bonds []}"#),
        Solution::Determined(()),
    )]
    #[case::aromaticity(
        mol_dsl!(r#"{:atoms ["C#a"]}"#),
        Solution::Contradictory(ValidatorContradiction::Aromaticity(
            AromaticityValidatorContradiction::Inconsistency(
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) },
            ),
        )),
    )]
    #[case::stereo(
        mol_dsl_ground!(r#"{
            :atoms ["C#h0#T1" "C#h3" "C#h3" "C#h3" "C#h3"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        }"#),
        Solution::Contradictory(ValidatorContradiction::Stereo(
            StereoValidatorContradiction::MissingStereoAtom { atom: AtomId(0) },
        )),
    )]
    #[case::aromaticity_partial(
        mol_dsl!(r#"{
            :atoms ["C#a+" "C#a" "C#a" "C#a" "C#a" "C#a"]
            :bonds [[0 1 "1"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        }"#),
        Solution::Underdetermined(()),
    )]
    fn test_validator_validate_conformance(
        #[case] molecule: MoleculeAst,
        #[case] expected: Solution<(), ValidatorContradiction>,
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

    #[rstest]
    #[case::methane(4, None, UnpairedElectronsForm::from((0_u8, 1_u8)), Solution::Determined(()))]
    #[case::with_valence_constraint(3, Some(1), UnpairedElectronsForm::from((0_u8, 1_u8)), Solution::Determined(()))]
    #[case::partial_spin(
        4,
        None,
        UnpairedElectronsForm { count: NumForm::Lit(0), multiplicity: NumForm::Undetermined },
        Solution::Underdetermined(()),
    )]
    #[case::invalid_spin(
        2,
        None,
        UnpairedElectronsForm::from((2_u8, 2_u8)),
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
        #[case] unpaired_electrons: UnpairedElectronsForm,
        #[case] expected: Solution<(), ValidatorContradiction>,
    ) {
        let mut atom = AtomAst::from_element(Element::C);
        atom.charge = NumForm::Lit(0);
        atom.lone_pairs = NumForm::Lit(0);
        atom.implicit_hydrogens = NumForm::Lit(hydrogens);
        atom.unpaired_electrons = unpaired_electrons;
        if let Some(v) = valence {
            atom.constraints
                .set(AtomConstraintAst::Valence(NumForm::Lit(v)));
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
                element: ElementForm::Lit(Element::C),
                charge: NumForm::Lit(0),
                implicit_hydrogens: NumForm::Lit(4 - i64::from(count)),
                lone_pairs: NumForm::Lit(0),
                unpaired_electrons: UnpairedElectronsForm::from((count, multiplicity)),
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
            let molecule = MoleculeAst::from_entries(MoleculeEntries {
                atoms: vec![AtomAst {
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
                .map_contradiction(ValidatorContradiction::SpinInvariants);
            let model = ChemistryModel::default();

            prop_assert_eq!(
                Validator::new(&model).validate_invariants(&molecule).unwrap(),
                expected,
            );
        }
    }
}
