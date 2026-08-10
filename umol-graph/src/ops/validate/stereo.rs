//! Stereo validator.

use thiserror::Error;
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    AsLit, AtomId, BondId, BooleanForm, ConstitutionColoring, GraphSymmetry, GraphSymmetryConfig,
    Lattice, LigandSymmetryForm, Molecule, StereoAtomId, StereoBondId, StereoKind,
    StereoLigandPair, StereoSymmetry, Stereogenicity, StereogenicityForm, Topicity, TopicityForm,
    TopicityRelationForm,
};
use umol_perm::OrientedPermutation;
use umol_utils::solution::Solution;

use crate::ops::model::StereoModel;
use crate::ops::stereo::{StereoInconsistency, StereoPerception};

/// Operational configuration for stereo-conformance validation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StereoValidateConfig {
    /// Backend used for molecule graph-automorphism calculations.
    pub automorphism_algorithm: AutomorphismAlgorithm,
    /// Maximum number of graph-symmetry refinement passes when the stereo
    /// model enables para-stereo perception.
    pub max_iterations: usize,
}

impl Default for StereoValidateConfig {
    fn default() -> Self {
        Self {
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            max_iterations: 16,
        }
    }
}

impl StereoValidateConfig {
    fn graph_symmetry_config(
        self,
        model: &StereoModel,
    ) -> GraphSymmetryConfig<ConstitutionColoring> {
        GraphSymmetryConfig {
            coloring: ConstitutionColoring::full(),
            iterate_to_fixpoint: model.para_stereo,
            max_iterations: self.max_iterations,
            automorphism_algorithm: self.automorphism_algorithm,
        }
    }
}

#[derive(Clone, Debug)]
pub struct StereoConformanceValidator {
    model: StereoModel,
    config: StereoValidateConfig,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum StereoValidatorContradiction {
    #[error(transparent)]
    Inconsistency(#[from] StereoInconsistency),
    #[error("tetrahedral stereo constraint at atom {atom:?} has no stereo atom")]
    MissingStereoAtom { atom: AtomId },
    #[error("cis-trans stereo constraint at bond {bond:?} has no stereo bond")]
    MissingStereoBond { bond: BondId },
    #[error("stereo {kind:?} coset {coset} out of range 0..{count}")]
    CosetOutOfRange {
        kind: StereoKind,
        coset: u32,
        count: usize,
    },
    #[error("stereo {kind:?} has {ligands} ligands, expected degree {degree}")]
    LigandArity {
        kind: StereoKind,
        ligands: usize,
        degree: usize,
    },
    #[error("improper (enantio-) relation asserted on achiral stereo {kind:?}")]
    ImproperOnAchiral { kind: StereoKind },
    #[error("stereogenicity {derived:?} contradicts asserted {asserted:?}")]
    StereogenicityMismatch {
        asserted: StereogenicityForm,
        derived: Stereogenicity,
    },
    #[error("ligand pair {pair:?} topicity {derived:?} contradicts asserted {asserted:?}")]
    TopicityMismatch {
        pair: StereoLigandPair,
        asserted: TopicityRelationForm,
        derived: Topicity,
    },
    #[error("asserted ligand symmetry {asserted:?} not satisfied by the derived group")]
    LigandSymmetryViolation { asserted: LigandSymmetryForm },
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum StereoValidatorError {}

impl StereoConformanceValidator {
    pub fn new(model: &StereoModel) -> Self {
        Self::with_config(model, StereoValidateConfig::default())
    }

    pub fn with_config(model: &StereoModel, config: StereoValidateConfig) -> Self {
        Self {
            model: model.clone(),
            config,
        }
    }

    /// Validate every stereo element against the molecule's graph symmetry.
    /// First `Contradictory` wins; a ground assertion the derived value leaves
    /// open contributes `Underdetermined`.
    pub fn validate(
        &self,
        molecule: &Molecule,
    ) -> Result<Solution<(), StereoValidatorContradiction>, StereoValidatorError> {
        let derivation = StereoPerception::new(&self.model).derive(molecule);
        for &inconsistency in &derivation.inconsistencies {
            let asserted = match inconsistency {
                StereoInconsistency::TetrahedralStereoFailure { .. }
                | StereoInconsistency::TetrahedralStereoMismatch { .. }
                | StereoInconsistency::CisTransStereoFailure { .. }
                | StereoInconsistency::CisTransStereoMismatch { .. } => true,
                StereoInconsistency::StereoAtomFailure { stereo_atom } => {
                    let atom = molecule.stereo_atom(stereo_atom).site_id();
                    molecule
                        .atom(atom)
                        .attributes
                        .constraints
                        .tetrahedral_stereo()
                        .is_some_and(|constraint| !constraint.is_undetermined())
                }
                StereoInconsistency::StereoBondFailure { stereo_bond } => {
                    let bond = molecule.stereo_bond(stereo_bond).site_id();
                    molecule
                        .bond(bond)
                        .attributes
                        .constraints
                        .cis_trans_stereo()
                        .is_some_and(|constraint| !constraint.is_undetermined())
                }
            };
            if asserted {
                return Ok(Solution::Contradictory(inconsistency.into()));
            }
        }

        for (atom, _, _) in &derivation.atoms {
            let asserted = molecule
                .atom(*atom)
                .attributes
                .constraints
                .tetrahedral_stereo()
                .is_some_and(|constraint| !constraint.is_undetermined());
            if asserted && !molecule.stereo_atoms().is_at(*atom) {
                return Ok(Solution::Contradictory(
                    StereoValidatorContradiction::MissingStereoAtom { atom: *atom },
                ));
            }
        }
        for (bond, _, _) in &derivation.bonds {
            let asserted = molecule
                .bond(*bond)
                .attributes
                .constraints
                .cis_trans_stereo()
                .is_some_and(|constraint| !constraint.is_undetermined());
            if asserted && !molecule.stereo_bonds().is_at(*bond) {
                return Ok(Solution::Contradictory(
                    StereoValidatorContradiction::MissingStereoBond { bond: *bond },
                ));
            }
        }

        let symmetry = molecule.graph_symmetry(&self.config.graph_symmetry_config(&self.model));
        let mut any_undetermined = false;

        for id in molecule.stereo_atoms().ids() {
            match self.validate_stereo_atom(molecule, id, &symmetry) {
                Solution::Determined(()) => {}
                Solution::Underdetermined(()) => any_undetermined = true,
                Solution::Contradictory(c) => return Ok(Solution::Contradictory(c)),
            }
        }
        for id in molecule.stereo_bonds().ids() {
            match self.validate_stereo_bond(molecule, id, &symmetry) {
                Solution::Determined(()) => {}
                Solution::Underdetermined(()) => any_undetermined = true,
                Solution::Contradictory(c) => return Ok(Solution::Contradictory(c)),
            }
        }

        Ok(if any_undetermined {
            Solution::Underdetermined(())
        } else {
            Solution::Determined(())
        })
    }

    fn validate_stereo_atom(
        &self,
        molecule: &Molecule,
        id: StereoAtomId,
        symmetry: &GraphSymmetry,
    ) -> Solution<(), StereoValidatorContradiction> {
        let view = molecule.stereo_atom(id);
        let constraints = view.constraints();
        let stereogenicity = constraints.stereogenicity();
        let topicities: Vec<TopicityForm> = constraints.topicities().cloned().collect();
        let ligand_symmetries: Vec<LigandSymmetryForm> =
            constraints.ligand_symmetries().cloned().collect();
        match self.validate_kind(
            view.kind(),
            view.coset().as_lit(),
            view.ligands().count(),
            &stereogenicity,
            &topicities,
        ) {
            Solution::Contradictory(c) => return Solution::Contradictory(c),
            Solution::Underdetermined(()) => return Solution::Underdetermined(()),
            Solution::Determined(()) => {}
        }
        let sym = molecule.stereo_atom_symmetry(symmetry, id);
        self.validate_symmetry(&sym, &stereogenicity, &topicities, &ligand_symmetries)
    }

    fn validate_stereo_bond(
        &self,
        molecule: &Molecule,
        id: StereoBondId,
        symmetry: &GraphSymmetry,
    ) -> Solution<(), StereoValidatorContradiction> {
        let view = molecule.stereo_bond(id);
        let constraints = view.constraints();
        let stereogenicity = constraints.stereogenicity();
        let topicities: Vec<TopicityForm> = constraints.topicities().cloned().collect();
        let ligand_symmetries: Vec<LigandSymmetryForm> =
            constraints.ligand_symmetries().cloned().collect();
        match self.validate_kind(
            view.kind(),
            view.coset().as_lit(),
            view.ligands().count(),
            &stereogenicity,
            &topicities,
        ) {
            Solution::Contradictory(c) => return Solution::Contradictory(c),
            Solution::Underdetermined(()) => return Solution::Underdetermined(()),
            Solution::Determined(()) => {}
        }
        let sym = molecule.stereo_bond_symmetry(symmetry, id);
        self.validate_symmetry(&sym, &stereogenicity, &topicities, &ligand_symmetries)
    }

    /// Ligand arity (= kind degree), the achiral-kind relation gate, and coset
    /// range; an absent coset leaves the element underdetermined.
    fn validate_kind(
        &self,
        kind: StereoKind,
        coset: Option<u32>,
        ligand_count: usize,
        stereogenicity: &StereogenicityForm,
        topicities: &[TopicityForm],
    ) -> Solution<(), StereoValidatorContradiction> {
        if ligand_count != kind.degree() {
            return Solution::Contradictory(StereoValidatorContradiction::LigandArity {
                kind,
                ligands: ligand_count,
                degree: kind.degree(),
            });
        }
        let has_improper = (!stereogenicity.is_undetermined()
            && stereogenicity.matches(&StereogenicityForm::Lit(Stereogenicity::Prochiral)))
            || topicities.iter().any(|t| {
                !t.relation.is_undetermined()
                    && t.relation
                        .matches(&TopicityRelationForm::Lit(Topicity::Enantiotopic))
            });
        if !kind.is_chiral_class() && has_improper {
            return Solution::Contradictory(StereoValidatorContradiction::ImproperOnAchiral {
                kind,
            });
        }
        match coset {
            Some(n) if (n as usize) >= kind.count() => {
                Solution::Contradictory(StereoValidatorContradiction::CosetOutOfRange {
                    kind,
                    coset: n,
                    count: kind.count(),
                })
            }
            Some(_) => Solution::Determined(()),
            None => Solution::Underdetermined(()),
        }
    }

    /// Compare derived symmetry against constraints.
    fn validate_symmetry(
        &self,
        sym: &StereoSymmetry,
        stereogenicity: &StereogenicityForm,
        topicities: &[TopicityForm],
        ligand_symmetries: &[LigandSymmetryForm],
    ) -> Solution<(), StereoValidatorContradiction> {
        let mut any_undetermined = false;

        if !stereogenicity.is_undetermined() {
            let derived = sym.stereogenicity();
            if !stereogenicity.matches(&StereogenicityForm::Lit(derived)) {
                return Solution::Contradictory(
                    StereoValidatorContradiction::StereogenicityMismatch {
                        asserted: stereogenicity.clone(),
                        derived,
                    },
                );
            }
            any_undetermined |= !stereogenicity.is_ground();
        }

        for t in topicities {
            if t.relation.is_undetermined() {
                continue;
            }
            let derived = sym.topicity(t.pair.first(), t.pair.second());
            if !t.relation.matches(&TopicityRelationForm::Lit(derived)) {
                return Solution::Contradictory(StereoValidatorContradiction::TopicityMismatch {
                    pair: t.pair,
                    asserted: t.relation.clone(),
                    derived,
                });
            }
            any_undetermined |= !t.relation.is_ground();
        }

        for ls in ligand_symmetries {
            let op =
                OrientedPermutation::new(ls.permutation.permutation.0, ls.permutation.orientation);
            let in_group = sym.group().contains(op);
            let holds = match ls.invariant {
                BooleanForm::Lit(true) => in_group,
                BooleanForm::Lit(false) => !in_group,
                BooleanForm::Undetermined => true,
            };
            if !holds {
                return Solution::Contradictory(
                    StereoValidatorContradiction::LigandSymmetryViolation { asserted: *ls },
                );
            }
        }

        if any_undetermined {
            Solution::Underdetermined(())
        } else {
            Solution::Determined(())
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph_ir::ir::{
        LigandPermutation, OrientedLigandPermutation, StereoAtomConstraintForm,
        StereoBondConstraintForm, StereoConfigurationForm, StereoKind, StereoLigandPosition,
        StereogenicityForm,
    };
    use umol_graph_ir::mol_dsl_ground;
    use umol_perm::{Orientation, Permutation};

    use super::*;

    // Tetrahedral stereocenter, four distinct halide ligands → trivial ligand
    // symmetry, genuinely stereogenic, every ligand pair diastereotopic.
    const CFCLBRI: &str = r#"{:atoms ["C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#;

    // Cis/trans double bond, each terminus bearing two distinct substituents.
    const BUTENE: &str = r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct1"}]}"#;

    #[rstest]
    fn test_stereo_validate_config_default() {
        assert_eq!(
            StereoValidateConfig::default(),
            StereoValidateConfig {
                automorphism_algorithm: AutomorphismAlgorithm::Nauty,
                max_iterations: 16,
            }
        );
    }

    #[rstest]
    #[case::no_para(false, false)]
    #[case::para(true, true)]
    fn test_stereo_validate_config_graph_symmetry_config(
        #[case] para_stereo: bool,
        #[case] expected_fixpoint: bool,
    ) {
        let model = StereoModel {
            para_stereo,
            ..StereoModel::default()
        };
        let config = StereoValidateConfig {
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            max_iterations: 8,
        }
        .graph_symmetry_config(&model);

        assert_eq!(config.automorphism_algorithm, AutomorphismAlgorithm::Nauty);
        assert_eq!(config.iterate_to_fixpoint, expected_fixpoint);
        assert_eq!(config.max_iterations, 8);
    }

    #[rstest]
    fn test_stereo_conformance_validator_new() {
        let molecule = mol_dsl_ground!(CFCLBRI);
        assert_eq!(
            StereoConformanceValidator::new(&StereoModel::default()).validate(&molecule),
            StereoConformanceValidator::with_config(
                &StereoModel::default(),
                StereoValidateConfig::default(),
            )
            .validate(&molecule)
        );
    }

    #[rstest]
    #[case::bare(CFCLBRI, (|_: &mut Molecule| {}) as fn(&mut Molecule))]
    #[case::matching_stereogenicity(
        CFCLBRI,
        (|molecule: &mut Molecule| {
            molecule.stereo_atom_mut(StereoAtomId(0)).attributes
                .constraints
                .set(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Stereogenic)));
        }) as fn(&mut Molecule)
    )]
    fn test_stereo_conformance_validator_validate(
        #[case] dsl: &str,
        #[case] mutate: fn(&mut Molecule),
    ) {
        let mut molecule = mol_dsl_ground!(dsl);
        mutate(&mut molecule);
        let solution = StereoConformanceValidator::new(&StereoModel::default())
            .validate(&molecule)
            .unwrap();
        assert_eq!(solution, Solution::Determined(()));
    }

    #[rstest]
    #[case::missing_stereo_atom(
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h#T1" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
        }"#),
        Solution::Contradictory(StereoValidatorContradiction::MissingStereoAtom {
            atom: AtomId(1),
        }),
    )]
    #[case::missing_stereo_bond(
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]
        }"#),
        Solution::Contradictory(StereoValidatorContradiction::MissingStereoBond {
            bond: BondId(1),
        }),
    )]
    #[case::unrealizable_tetrahedral_stereo(
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "S#h0#T1" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "1"]]
        }"#),
        Solution::Contradictory(StereoValidatorContradiction::Inconsistency(
            StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(1) },
        )),
    )]
    #[case::unrealizable_cis_trans_stereo(
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h2" "C#h"]
            :bonds [[0 1 "1"] [1 2 "2#C1"]]
        }"#),
        Solution::Contradictory(StereoValidatorContradiction::Inconsistency(
            StereoInconsistency::CisTransStereoFailure { bond: BondId(1) },
        )),
    )]
    #[case::unrealizable_stereo_atom(
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h#T1" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
            :stereo-atoms [{:site 1 :ligands [0] :attrs "Th1"}]
        }"#),
        Solution::Contradictory(StereoValidatorContradiction::Inconsistency(
            StereoInconsistency::StereoAtomFailure {
                stereo_atom: StereoAtomId(0),
            },
        )),
    )]
    #[case::unrealizable_stereo_bond(
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]
            :stereo-bonds [{:site 1 :ligands [0] :attrs "Ct1"}]
        }"#),
        Solution::Contradictory(StereoValidatorContradiction::Inconsistency(
            StereoInconsistency::StereoBondFailure {
                stereo_bond: StereoBondId(0),
            },
        )),
    )]
    #[case::tetrahedral_stereo_mismatch(
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h#T1" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
            :stereo-atoms [{
                :site 1
                :ligands [0 2 3 [:h 1]]
                :attrs "Th0"
            }]
        }"#),
        Solution::Contradictory(StereoValidatorContradiction::Inconsistency(
            StereoInconsistency::TetrahedralStereoMismatch {
                atom: AtomId(1),
                stereo_atom: StereoAtomId(0),
            },
        )),
    )]
    #[case::cis_trans_stereo_mismatch(
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]
            :stereo-bonds [{
                :site 1
                :ligands [0 [:h 1] 3 [:h 2]]
                :attrs "Ct0"
            }]
        }"#),
        Solution::Contradictory(StereoValidatorContradiction::Inconsistency(
            StereoInconsistency::CisTransStereoMismatch {
                bond: BondId(1),
                stereo_bond: StereoBondId(0),
            },
        )),
    )]
    #[case::conformant_stereo_atom(
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h#T1" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
            :stereo-atoms [{
                :site 1
                :ligands [0 2 3 [:h 1]]
                :attrs "Th1"
            }]
        }"#),
        Solution::Determined(()),
    )]
    #[case::conformant_stereo_bond(
        mol_dsl_ground!(r#"{
            :atoms ["C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]
            :stereo-bonds [{
                :site 1
                :ligands [0 [:h 1] 3 [:h 2]]
                :attrs "Ct1"
            }]
        }"#),
        Solution::Determined(()),
    )]
    #[case::undetermined_tetrahedral_stereo(
        mol_dsl_ground!(r#"{:atoms ["C#T*"]}"#),
        Solution::Determined(()),
    )]
    #[case::undetermined_cis_trans_stereo(
        mol_dsl_ground!(r#"{
            :atoms ["C#h2" "C#h2"]
            :bonds [[0 1 "2#C*"]]
        }"#),
        Solution::Determined(()),
    )]
    fn test_stereo_conformance_validator_validate_constraint(
        #[case] molecule: Molecule,
        #[case] expected: Solution<(), StereoValidatorContradiction>,
    ) {
        assert_eq!(
            StereoConformanceValidator::new(&StereoModel::default())
                .validate(&molecule)
                .unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::coset_out_of_range(
        CFCLBRI,
        (|molecule: &mut Molecule| {
            molecule.stereo_atom_mut(StereoAtomId(0)).attributes.configuration = StereoConfigurationForm::kinded(StereoKind::Tetrahedral, 9);
        }) as fn(&mut Molecule),
        StereoValidatorContradiction::CosetOutOfRange {
            kind: StereoKind::Tetrahedral,
            coset: 9,
            count: 2,
        }
    )]
    #[case::arity(
        CFCLBRI,
        (|molecule: &mut Molecule| {
            molecule.stereo_atom_mut(StereoAtomId(0)).attributes.configuration = StereoConfigurationForm::kinded(StereoKind::TrigonalBipyramidal, 0);
        }) as fn(&mut Molecule),
        StereoValidatorContradiction::LigandArity {
            kind: StereoKind::TrigonalBipyramidal,
            ligands: 4,
            degree: 5,
        }
    )]
    #[case::improper_on_achiral(
        BUTENE,
        (|molecule: &mut Molecule| {
            molecule.stereo_bond_mut(StereoBondId(0)).attributes
                .constraints
                .set(StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Prochiral)));
        }) as fn(&mut Molecule),
        StereoValidatorContradiction::ImproperOnAchiral {
            kind: StereoKind::CisTrans,
        }
    )]
    #[case::stereogenicity_mismatch(
        CFCLBRI,
        (|molecule: &mut Molecule| {
            molecule.stereo_atom_mut(StereoAtomId(0)).attributes
                .constraints
                .set(StereoAtomConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Symmetric)));
        }) as fn(&mut Molecule),
        StereoValidatorContradiction::StereogenicityMismatch {
            asserted: StereogenicityForm::Lit(Stereogenicity::Symmetric),
            derived: Stereogenicity::Stereogenic,
        }
    )]
    #[case::topicity_mismatch(
        CFCLBRI,
        (|molecule: &mut Molecule| {
            molecule.stereo_atom_mut(StereoAtomId(0)).attributes
                .constraints
                .set(StereoAtomConstraintForm::Topicity(TopicityForm {
                    pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
                    relation: TopicityRelationForm::Lit(Topicity::Homotopic),
                }));
        }) as fn(&mut Molecule),
        StereoValidatorContradiction::TopicityMismatch {
            pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
            asserted: TopicityRelationForm::Lit(Topicity::Homotopic),
            derived: Topicity::Diastereotopic,
        }
    )]
    #[case::ligand_symmetry_violation(
        CFCLBRI,
        (|molecule: &mut Molecule| {
            molecule.stereo_atom_mut(StereoAtomId(0)).attributes
                .constraints
                .set(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
                    permutation: OrientedLigandPermutation {
                        permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
                        orientation: Orientation::Proper,
                    },
                    invariant: BooleanForm::Lit(true),
                }));
        }) as fn(&mut Molecule),
        StereoValidatorContradiction::LigandSymmetryViolation {
            asserted: LigandSymmetryForm {
                permutation: OrientedLigandPermutation {
                    permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
                    orientation: Orientation::Proper,
                },
                invariant: BooleanForm::Lit(true),
            },
        }
    )]
    fn test_stereo_conformance_validator_validate_error(
        #[case] dsl: &str,
        #[case] mutate: fn(&mut Molecule),
        #[case] expected: StereoValidatorContradiction,
    ) {
        let mut molecule = mol_dsl_ground!(dsl);
        mutate(&mut molecule);
        let solution = StereoConformanceValidator::new(&StereoModel::default())
            .validate(&molecule)
            .unwrap();
        assert_eq!(solution, Solution::Contradictory(expected));
    }
}
