//! Stereo validator.

use thiserror::Error;
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{
    AsLit, BondId, BooleanForm, CisTransStereoForm, ConstitutionColoring, FrameTransport,
    GraphSymmetry, GraphSymmetryConfig, Lattice, LigandSymmetryForm, Molecule, RingConfig,
    RingModel, RingSetKind, StereoAtomId, StereoBondId, StereoKind, StereoLigandKind,
    StereoLigandPair, StereoSymmetry, Stereogenicity, StereogenicityForm, Topicity, TopicityForm,
    TopicityRelationForm,
};
use umol_perm::{OrientedPermutation, Permutation};
use umol_utils::solution::Solution;

use crate::ops::model::StereoModel;
use crate::ops::stereo::{
    StereoContradiction, StereoDerivation, StereoInconsistency, StereoPerceiver,
};

/// Operational configuration for stereo-conformance validation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StereoValidateConfig {
    /// Backend used for molecule graph-automorphism calculations.
    pub automorphism_algorithm: AutomorphismAlgorithm,
    pub para_stereo_max_iterations: u32,
    pub ring_config: RingConfig,
}

impl Default for StereoValidateConfig {
    fn default() -> Self {
        Self {
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
            para_stereo_max_iterations: 16,
            ring_config: RingConfig::default(),
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
            max_iterations: self.para_stereo_max_iterations,
            automorphism_algorithm: self.automorphism_algorithm,
        }
    }
}

/// Validates stored stereo entities and constraints against a stereo model.
#[derive(Clone, Debug)]
pub struct StereoConformanceValidator {
    model: StereoModel,
    config: StereoValidateConfig,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum StereoConformanceContradiction {
    #[error(transparent)]
    Inconsistency(#[from] StereoInconsistency),
    #[error(
        "bond {bond:?} has stereo in a ring of size {ring_size}, below minimum {minimum_ring_size}"
    )]
    RingTooSmall {
        bond: BondId,
        ring_size: u32,
        minimum_ring_size: u32,
    },
    #[error("bond {bond:?} has ring-relative E stereo in a ring of size {ring_size}, below minimum {minimum_ring_size}")]
    RingTooSmallForE {
        bond: BondId,
        ring_size: u32,
        minimum_ring_size: u32,
    },
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
    /// An asserted topicity names a position outside the ligand frame.
    #[error("ligand pair {pair:?} is out of range for degree {degree}")]
    TopicityPositionOutOfRange {
        pair: StereoLigandPair,
        degree: usize,
    },
    #[error("asserted ligand symmetry {asserted:?} not satisfied by the derived group")]
    LigandSymmetryViolation { asserted: LigandSymmetryForm },
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum StereoConformanceError {}

impl StereoConformanceValidator {
    /// Construct a validator with the default stereo operation configuration.
    pub fn new(model: &StereoModel) -> Self {
        Self::with_config(model, StereoValidateConfig::default())
    }

    /// Construct a validator with explicit stereo operation configuration.
    pub fn with_config(model: &StereoModel, config: StereoValidateConfig) -> Self {
        Self {
            model: model.clone(),
            config,
        }
    }

    /// Validate stereo descriptions against the model's ring-size limits and graph symmetry.
    /// First `Contradictory` wins; a ground assertion the derived value leaves
    /// open contributes `Underdetermined`.
    /// Unresolved ligand identification contributes `Underdetermined` after
    /// established contradictions are checked.
    /// Asserted topicity positions must lie within the ligand frame, and ligand-symmetry
    /// permutations must have its degree; disagreements produce `Contradictory`.
    pub fn validate(
        &self,
        molecule: &Molecule,
    ) -> Result<Solution<(), StereoConformanceContradiction>, StereoConformanceError> {
        let (derivation, mut any_undetermined) =
            match StereoPerceiver::new(&self.model).derive(molecule, self.config.ring_config) {
                Solution::Determined(derivation) => (derivation, false),
                Solution::Underdetermined(derivation) => (derivation, true),
                Solution::Contradictory(StereoContradiction::Inconsistency(inconsistency)) => {
                    return Ok(Solution::Contradictory(inconsistency.into()));
                }
            };
        for &inconsistency in &derivation.inconsistencies {
            let asserted = match inconsistency {
                StereoInconsistency::TetrahedralStereoFailure { .. }
                | StereoInconsistency::TetrahedralStereoMismatch { .. }
                | StereoInconsistency::CisTransStereoFailure { .. }
                | StereoInconsistency::CisTransStereoMismatch { .. } => true,
                StereoInconsistency::StereoAtomFailure { stereo_atom } => {
                    let atom = molecule.stereo_atom(stereo_atom).site_id();
                    !molecule
                        .atom(atom)
                        .attributes()
                        .constraints
                        .tetrahedral_stereo()
                        .is_undetermined()
                }
                StereoInconsistency::StereoBondFailure { stereo_bond } => {
                    let bond = molecule.stereo_bond(stereo_bond).site_id();
                    !molecule
                        .bond(bond)
                        .attributes()
                        .constraints
                        .cis_trans_stereo()
                        .is_undetermined()
                }
            };
            if asserted {
                return Ok(Solution::Contradictory(inconsistency.into()));
            }
        }

        match self.validate_rings(molecule, &derivation) {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c)),
        }

        let symmetry = molecule.graph_symmetry(&self.config.graph_symmetry_config(&self.model));
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

    fn validate_rings(
        &self,
        molecule: &Molecule,
        derivation: &StereoDerivation,
    ) -> Solution<(), StereoConformanceContradiction> {
        let minimum = self.model.stereo_bond_min_ring_size;
        let e_minimum = self.model.stereo_bond_e_min_ring_size;
        let maximum = minimum.max(e_minimum);
        if maximum <= 3 {
            return Solution::Determined(());
        }
        let mut rings = None;
        let mut any_undetermined = false;
        for bond in molecule.bonds().iter() {
            let relation = bond.stereo_bond();
            if bond.order().as_lit() != Some(2)
                || !(relation.is_some_and(|view| {
                    view.attributes().configuration.kind() == Some(StereoKind::CisTrans)
                }) || matches!(
                    bond.constraints().cis_trans_stereo(),
                    CisTransStereoForm::Stereo(_)
                ))
            {
                continue;
            }
            let rings = rings.get_or_insert_with(|| {
                molecule
                    .rings(
                        RingModel {
                            kind: RingSetKind::Simple,
                            max_ring_size: maximum as usize - 1,
                        },
                        self.config.ring_config,
                    )
                    .into_ring_set()
            });
            if let Some(size) = rings.bond_smallest_ring_size(bond.id()) {
                if size < minimum as usize {
                    return Solution::Contradictory(StereoConformanceContradiction::RingTooSmall {
                        bond: bond.id(),
                        ring_size: size as u32,
                        minimum_ring_size: minimum,
                    });
                }
            }
            let frame = relation
                .map(|view| (view.ligand_ids(), &view.attributes().configuration))
                .or_else(|| {
                    derivation
                        .stereo_bonds
                        .iter()
                        .find(|(site, _, _)| *site == bond.id())
                        .map(|(_, ligands, attributes)| {
                            (ligands.as_slice(), &attributes.configuration)
                        })
                });
            for ring in rings
                .iter()
                .filter(|ring| ring.bonds().contains(&bond.id()) && ring.len() < e_minimum as usize)
            {
                let Some((ligands, configuration)) = frame else {
                    any_undetermined = true;
                    continue;
                };
                let mut image = [0, 1, 2, 3];
                let mut complete = true;
                for side in image.chunks_exact_mut(2) {
                    let position = side.iter().position(|&index| {
                        let ligand = ligands[index];
                        ligand.kind == StereoLigandKind::Atom
                            && bond.atom_ids().iter().any(|&endpoint| {
                                molecule
                                    .bonds()
                                    .of_id(endpoint, ligand.atom_id)
                                    .is_some_and(|id| id != bond.id() && ring.bonds().contains(&id))
                            })
                    });
                    if let Some(position) = position {
                        side.swap(0, position);
                    } else {
                        complete = false;
                    }
                }
                if !complete {
                    any_undetermined = true;
                    continue;
                }
                let coset = configuration
                    .clone()
                    .reframe_by(&Permutation::from_image(&image))
                    .and_then(|configuration| configuration.coset().and_then(AsLit::as_lit));
                if coset == Some(1) {
                    return Solution::Contradictory(
                        StereoConformanceContradiction::RingTooSmallForE {
                            bond: bond.id(),
                            ring_size: ring.len() as u32,
                            minimum_ring_size: e_minimum,
                        },
                    );
                }
                any_undetermined |= coset.is_none();
            }
        }
        if any_undetermined {
            Solution::Underdetermined(())
        } else {
            Solution::Determined(())
        }
    }

    fn validate_stereo_atom(
        &self,
        molecule: &Molecule,
        id: StereoAtomId,
        symmetry: &GraphSymmetry,
    ) -> Solution<(), StereoConformanceContradiction> {
        let view = molecule.stereo_atom(id);
        let constraints = view.constraints();
        let stereogenicity = constraints.stereogenicity();
        if view
            .attributes()
            .configuration
            .coset()
            .and_then(AsLit::as_lit)
            .is_none()
        {
            return Solution::Underdetermined(());
        }
        let sym = molecule.stereo_atom_symmetry(symmetry, id);
        self.validate_symmetry(
            &sym,
            stereogenicity,
            constraints.topicities(),
            constraints.ligand_symmetries(),
        )
    }

    fn validate_stereo_bond(
        &self,
        molecule: &Molecule,
        id: StereoBondId,
        symmetry: &GraphSymmetry,
    ) -> Solution<(), StereoConformanceContradiction> {
        let view = molecule.stereo_bond(id);
        let constraints = view.constraints();
        let stereogenicity = constraints.stereogenicity();
        if view
            .attributes()
            .configuration
            .coset()
            .and_then(AsLit::as_lit)
            .is_none()
        {
            return Solution::Underdetermined(());
        }
        let sym = molecule.stereo_bond_symmetry(symmetry, id);
        self.validate_symmetry(
            &sym,
            stereogenicity,
            constraints.topicities(),
            constraints.ligand_symmetries(),
        )
    }

    /// Compare derived symmetry against constraints.
    fn validate_symmetry<'a>(
        &self,
        sym: &StereoSymmetry,
        stereogenicity: &StereogenicityForm,
        topicities: impl Iterator<Item = &'a TopicityForm>,
        ligand_symmetries: impl Iterator<Item = &'a LigandSymmetryForm>,
    ) -> Solution<(), StereoConformanceContradiction> {
        let mut any_undetermined = false;
        let degree = sym.group().degree();

        if !stereogenicity.is_undetermined() {
            let derived = sym.stereogenicity();
            if !stereogenicity.matches(&StereogenicityForm::Lit(derived)) {
                return Solution::Contradictory(
                    StereoConformanceContradiction::StereogenicityMismatch {
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
            if t.pair.first().index() >= degree || t.pair.second().index() >= degree {
                return Solution::Contradictory(
                    StereoConformanceContradiction::TopicityPositionOutOfRange {
                        pair: t.pair,
                        degree,
                    },
                );
            }
            let derived = sym.topicity(t.pair.first(), t.pair.second());
            if !t.relation.matches(&TopicityRelationForm::Lit(derived)) {
                return Solution::Contradictory(StereoConformanceContradiction::TopicityMismatch {
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
            if op.degree() != degree {
                return Solution::Contradictory(
                    StereoConformanceContradiction::LigandSymmetryViolation { asserted: *ls },
                );
            }
            let in_group = sym.group().contains(op);
            let holds = match ls.invariant {
                BooleanForm::Lit(true) => in_group,
                BooleanForm::Lit(false) => !in_group,
                BooleanForm::Undetermined => true,
            };
            if !holds {
                return Solution::Contradictory(
                    StereoConformanceContradiction::LigandSymmetryViolation { asserted: *ls },
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
    use std::collections::BTreeSet;

    use rstest::rstest;
    use umol_chem::element::Element;
    use umol_graph_core::{GraphRemapping, NodeId, Remapping};
    use umol_graph_ir::ir::{
        AtomForm, AtomId, BondConstraintForm, BondForm, BondId, CisTransStereoForm,
        LigandPermutation, MoleculeEntries, MoleculeRemapping, NumForm, OrientedLigandPermutation,
        StereoAtomConstraintForm, StereoBondConstraintForm, StereoBondForm,
        StereoConfigurationForm, StereoCoset, StereoKind, StereoLigand, StereoLigandKind,
        StereoLigandPosition, StereoTerm, StereogenicityForm,
    };
    use umol_graph_ir::mol_dsl_concrete;
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
                para_stereo_max_iterations: 16,
                ring_config: Default::default(),
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
            para_stereo_max_iterations: 8,
            ring_config: Default::default(),
        }
        .graph_symmetry_config(&model);

        assert_eq!(config.automorphism_algorithm, AutomorphismAlgorithm::Nauty);
        assert_eq!(config.iterate_to_fixpoint, expected_fixpoint);
        assert_eq!(config.max_iterations, 8);
    }

    #[rstest]
    fn test_stereo_conformance_validator_new() {
        let molecule = mol_dsl_concrete!(CFCLBRI);
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
    #[case::atom_hydrogen_and_lone_pair(
        r#"{:atoms ["P#h#n#T1" "F" "Cl"] :bonds [[0 1 "1"] [0 2 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 [:h 0] [:lp 0]] :attrs "Th1"}]}"#,
        (|_: &mut Molecule| {}) as fn(&mut Molecule),
    )]
    #[case::bond_hydrogen_and_lone_pair(
        r#"{:atoms ["N#h#n" "N#h#n"] :bonds [[0 1 "2#C1"]]
            :stereo-bonds [{:site 0 :ligands [[:h 0] [:lp 0] [:h 1] [:lp 1]] :attrs "Ct1"}]}"#,
        (|_: &mut Molecule| {}) as fn(&mut Molecule),
    )]
    #[case::matching_stereogenicity(
        CFCLBRI,
        (|molecule: &mut Molecule| {
            molecule
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .constraints
                .set(StereoAtomConstraintForm::Stereogenicity(
                    StereogenicityForm::Lit(Stereogenicity::Stereogenic),
                ));
        }) as fn(&mut Molecule)
    )]
    #[case::non_stereogenic_cis_trans(
        r#"{:atoms ["C" "C" "F" "F" "Br" "I"]
            :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct0"}]}"#,
        (|molecule: &mut Molecule| {
            molecule.stereo_bond_mut(StereoBondId(0)).attributes_mut().constraints.set(
                StereoBondConstraintForm::Stereogenicity(StereogenicityForm::Lit(Stereogenicity::Symmetric)),
            );
        }) as fn(&mut Molecule)
    )]
    fn test_stereo_conformance_validator_validate(
        #[case] dsl: &str,
        #[case] mutate: fn(&mut Molecule),
    ) {
        let mut molecule = mol_dsl_concrete!(dsl);
        mutate(&mut molecule);
        let solution = StereoConformanceValidator::new(&StereoModel::default())
            .validate(&molecule)
            .unwrap();
        assert_eq!(solution, Solution::Determined(()));
    }

    #[rstest]
    #[case::missing_stereo_atom(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h#T1" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
        }"#),
        Solution::Determined(()),
    )]
    #[case::missing_stereo_bond(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]
        }"#),
        Solution::Determined(()),
    )]
    #[case::unrealizable_tetrahedral_stereo(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "S#h0#T1" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "1"]]
        }"#),
        Solution::Contradictory(StereoConformanceContradiction::Inconsistency(
            StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(1) },
        )),
    )]
    #[case::unrealizable_cis_trans_stereo(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h2" "C#h"]
            :bonds [[0 1 "1"] [1 2 "2#C1"]]
        }"#),
        Solution::Contradictory(StereoConformanceContradiction::Inconsistency(
            StereoInconsistency::CisTransStereoFailure { bond: BondId(1) },
        )),
    )]
    #[case::tetrahedral_stereo_mismatch(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h#T1" "N#h2" "O#h"]
            :bonds [[0 1 "1"] [1 2 "1"] [1 3 "1"]]
            :stereo-atoms [{
                :site 1
                :ligands [0 2 3 [:h 1]]
                :attrs "Th0"
            }]
        }"#),
        Solution::Contradictory(StereoConformanceContradiction::Inconsistency(
            StereoInconsistency::TetrahedralStereoMismatch {
                atom: AtomId(1),
                stereo_atom: StereoAtomId(0),
            },
        )),
    )]
    #[case::cis_trans_stereo_mismatch(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h3" "C#h" "C#h" "C#h3"]
            :bonds [[0 1 "1"] [1 2 "2#C1"] [2 3 "1"]]
            :stereo-bonds [{
                :site 1
                :ligands [0 [:h 1] 3 [:h 2]]
                :attrs "Ct0"
            }]
        }"#),
        Solution::Contradictory(StereoConformanceContradiction::Inconsistency(
            StereoInconsistency::CisTransStereoMismatch {
                bond: BondId(1),
                stereo_bond: StereoBondId(0),
            },
        )),
    )]
    #[case::conformant_stereo_atom(
        mol_dsl_concrete!(r#"{
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
        mol_dsl_concrete!(r#"{
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
        mol_dsl_concrete!(r#"{:atoms ["C#T*"]}"#),
        Solution::Determined(()),
    )]
    #[case::undetermined_cis_trans_stereo(
        mol_dsl_concrete!(r#"{
            :atoms ["C#h2" "C#h2"]
            :bonds [[0 1 "2#C*"]]
        }"#),
        Solution::Determined(()),
    )]
    fn test_stereo_conformance_validator_validate_constraint(
        #[case] molecule: Molecule,
        #[case] expected: Solution<(), StereoConformanceContradiction>,
    ) {
        assert_eq!(
            StereoConformanceValidator::new(&StereoModel::default())
                .validate(&molecule)
                .unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::below_general(7, 8, 0, StereoCoset::Lit(0), Solution::Contradictory(
        StereoConformanceContradiction::RingTooSmall {
            bond: BondId(0), ring_size: 7, minimum_ring_size: 8,
        },
    ))]
    #[case::general_boundary(8, 8, 0, StereoCoset::Lit(1), Solution::Determined(()))]
    #[case::z_below_e(7, 0, 8, StereoCoset::Lit(0), Solution::Determined(()))]
    #[case::e_below_e(7, 0, 8, StereoCoset::Lit(1), Solution::Contradictory(
        StereoConformanceContradiction::RingTooSmallForE {
            bond: BondId(0), ring_size: 7, minimum_ring_size: 8,
        },
    ))]
    #[case::e_boundary(8, 0, 8, StereoCoset::Lit(1), Solution::Determined(()))]
    #[case::both_boundaries(8, 8, 8, StereoCoset::Lit(1), Solution::Determined(()))]
    #[case::unrestricted(4, 0, 0, StereoCoset::Lit(1), Solution::Determined(()))]
    #[case::small_thresholds(4, 2, 3, StereoCoset::Lit(1), Solution::Determined(()))]
    #[case::open_below_general(7, 8, 8, StereoCoset::Undetermined, Solution::Contradictory(
        StereoConformanceContradiction::RingTooSmall {
            bond: BondId(0), ring_size: 7, minimum_ring_size: 8,
        },
    ))]
    #[case::open_below_e(7, 0, 8, StereoCoset::Undetermined, Solution::Underdetermined(()))]
    #[case::variable_below_e(7, 0, 8, StereoCoset::term(StereoTerm::var("configuration")), Solution::Underdetermined(()))]
    #[case::large_below_general(24, 25, 0, StereoCoset::Lit(0), Solution::Contradictory(
        StereoConformanceContradiction::RingTooSmall {
            bond: BondId(0), ring_size: 24, minimum_ring_size: 25,
        },
    ))]
    #[case::large_below_e(24, 0, 25, StereoCoset::Lit(1), Solution::Contradictory(
        StereoConformanceContradiction::RingTooSmallForE {
            bond: BondId(0), ring_size: 24, minimum_ring_size: 25,
        },
    ))]
    #[case::large_boundary(24, 24, 24, StereoCoset::Lit(1), Solution::Determined(()))]
    fn test_stereo_conformance_validator_validate_ring(
        #[case] size: u32,
        #[case] stereo_bond_min_ring_size: u32,
        #[case] stereo_bond_e_min_ring_size: u32,
        #[case] coset: StereoCoset,
        #[case] expected: Solution<(), StereoConformanceContradiction>,
        #[values(false, true)] unresolved_atom: bool,
    ) {
        let mut molecule = Molecule::from_entries(MoleculeEntries {
            atoms: (0..size)
                .map(|id| {
                    AtomForm::from_element(Element::C)
                        .with_implicit_hydrogens(if id < 2 { 1 } else { 2 })
                        .into_concrete()
                })
                .collect(),
            bonds: (0..size)
                .map(|id| {
                    (
                        AtomId(id),
                        AtomId((id + 1) % size),
                        BondForm::from_order(if id == 0 { 2 } else { 1 }).into_concrete(),
                    )
                })
                .collect(),
            stereo_bonds: vec![(
                BondId(0),
                vec![
                    StereoLigand::new(AtomId(size - 1), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
                    StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                    StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
                ],
                StereoBondForm::new(StereoKind::CisTrans, coset),
            )],
            ..Default::default()
        });
        if unresolved_atom {
            let mut atom = mol_dsl_concrete!(CFCLBRI);
            atom.stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .configuration = StereoConfigurationForm::Undetermined;
            molecule.combine_from(&atom);
        }
        let expected = if unresolved_atom && expected == Solution::Determined(()) {
            Solution::Underdetermined(())
        } else {
            expected
        };
        let original = molecule.clone();
        assert_eq!(
            StereoConformanceValidator::new(&StereoModel {
                stereo_bond_min_ring_size,
                stereo_bond_e_min_ring_size,
                ..Default::default()
            })
            .validate(&molecule),
            Ok(expected),
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::unasserted(CisTransStereoForm::Undetermined, 8, 8, Solution::Determined(()))]
    #[case::negative(CisTransStereoForm::NotStereo, 8, 8, Solution::Determined(()))]
    #[case::positive_below_general(CisTransStereoForm::Stereo(StereoCoset::Undetermined), 8, 8,
        Solution::Contradictory(StereoConformanceContradiction::RingTooSmall {
            bond: BondId(0), ring_size: 6, minimum_ring_size: 8,
        }))]
    #[case::positive_below_e(CisTransStereoForm::Stereo(StereoCoset::Undetermined), 6, 8,
        Solution::Underdetermined(()))]
    #[case::positive_unrestricted(CisTransStereoForm::Stereo(StereoCoset::Undetermined), 0, 0,
        Solution::Determined(()))]
    #[case::z_below_e(CisTransStereoForm::Stereo(StereoCoset::Lit(0)), 6, 8,
        Solution::Determined(()))]
    #[case::e_below_e(CisTransStereoForm::Stereo(StereoCoset::Lit(1)), 6, 8,
        Solution::Contradictory(StereoConformanceContradiction::RingTooSmallForE {
            bond: BondId(0), ring_size: 6, minimum_ring_size: 8,
        }))]
    fn test_stereo_conformance_validator_validate_ring_constraint(
        #[case] assertion: CisTransStereoForm,
        #[case] stereo_bond_min_ring_size: u32,
        #[case] stereo_bond_e_min_ring_size: u32,
        #[case] expected: Solution<(), StereoConformanceContradiction>,
    ) {
        let mut molecule = mol_dsl_concrete!(
            r#"{
            :atoms ["C#h" "C#h" "C#h2" "C#h2" "C#h2" "C#h2"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
        }"#
        );
        molecule
            .bond_mut(BondId(0))
            .attributes_mut()
            .constraints
            .set(BondConstraintForm::CisTransStereo(assertion));
        let original = molecule.clone();
        assert_eq!(
            StereoConformanceValidator::new(&StereoModel {
                stereo_bond_min_ring_size,
                stereo_bond_e_min_ring_size,
                ..Default::default()
            })
            .validate(&molecule),
            Ok(expected),
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    fn test_stereo_conformance_validator_validate_ring_transport(
        #[values(0, 1)] coset: u32,
        #[values(false, true)] swap_first: bool,
        #[values(false, true)] swap_second: bool,
        #[values(false, true)] exchange_sides: bool,
        #[values(false, true)] reverse_atoms: bool,
    ) {
        let mut frame = [
            StereoLigand::new(AtomId(6), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
            StereoLigand::new(AtomId(1), StereoLigandKind::ImplicitHydrogen),
        ];
        if swap_first {
            frame.swap(0, 1);
        }
        if swap_second {
            frame.swap(2, 3);
        }
        if exchange_sides {
            frame.rotate_left(2);
        }
        let molecule = Molecule::from_entries(MoleculeEntries {
            atoms: (0..7)
                .map(|id| {
                    AtomForm::from_element(Element::C)
                        .with_implicit_hydrogens(if id < 2 { 1 } else { 2 })
                        .into_concrete()
                })
                .collect(),
            bonds: (0..7)
                .map(|id| {
                    (
                        AtomId(id),
                        AtomId((id + 1) % 7),
                        BondForm::from_order(if id == 0 { 2 } else { 1 }).into_concrete(),
                    )
                })
                .collect(),
            stereo_bonds: vec![(
                BondId(0),
                frame.to_vec(),
                StereoBondForm::new(
                    StereoKind::CisTrans,
                    StereoCoset::Lit(coset ^ u32::from(swap_first) ^ u32::from(swap_second)),
                ),
            )],
            ..Default::default()
        });
        let molecule = if reverse_atoms {
            molecule.remap(&MoleculeRemapping::new(
                GraphRemapping::new(
                    Remapping::new((0..7).rev().map(NodeId).collect()).unwrap(),
                    Remapping::identity(7),
                ),
                Remapping::identity(0),
                Remapping::identity(0),
                Remapping::identity(0),
                Remapping::identity(0),
                Remapping::identity(0),
                Remapping::identity(1),
            ))
        } else {
            molecule
        };
        let original = molecule.clone();
        assert_eq!(
            StereoConformanceValidator::new(&StereoModel {
                stereo_bond_min_ring_size: 0,
                stereo_bond_e_min_ring_size: 8,
                ..Default::default()
            })
            .validate(&molecule),
            Ok(if coset == 0 {
                Solution::Determined(())
            } else {
                Solution::Contradictory(StereoConformanceContradiction::RingTooSmallForE {
                    bond: BondId(0),
                    ring_size: 7,
                    minimum_ring_size: 8,
                })
            }),
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::smaller_z_larger_e(0, 6, Solution::Contradictory(
        StereoConformanceContradiction::RingTooSmallForE {
            bond: BondId(0), ring_size: 5, minimum_ring_size: 6,
        },
    ))]
    #[case::larger_boundary(0, 5, Solution::Determined(()))]
    #[case::smaller_e_larger_z(1, 5, Solution::Contradictory(
        StereoConformanceContradiction::RingTooSmallForE {
            bond: BondId(0), ring_size: 4, minimum_ring_size: 5,
        },
    ))]
    #[case::smaller_boundary(1, 4, Solution::Determined(()))]
    fn test_stereo_conformance_validator_validate_rings(
        #[case] coset: u32,
        #[case] stereo_bond_e_min_ring_size: u32,
        #[case] expected: Solution<(), StereoConformanceContradiction>,
    ) {
        let mut molecule = mol_dsl_concrete!(
            r#"{
            :atoms ["C" "C#h" "C#h2" "C#h" "C#h2" "C#h2"]
            :bonds [[0 1 "2"] [0 2 "1"] [2 3 "1"] [3 1 "1"]
                    [0 4 "1"] [4 5 "1"] [5 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 4 3 [:h 1]] :attrs "Ct0"}]
        }"#
        );
        molecule
            .stereo_bond_mut(StereoBondId(0))
            .attributes_mut()
            .configuration =
            StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(coset));
        let original = molecule.clone();
        assert_eq!(
            StereoConformanceValidator::new(&StereoModel {
                stereo_bond_min_ring_size: 0,
                stereo_bond_e_min_ring_size,
                ..Default::default()
            })
            .validate(&molecule),
            Ok(expected),
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    fn test_stereo_conformance_validator_validate_ring_partial(
        #[values(false, true)] known_first: bool,
    ) {
        let unknown = mol_dsl_concrete!(
            r#"{
            :atoms ["C#h" "C#h" "C#h2" "C#h2" "C#h2" "C#h2"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :stereo-bonds [{:site 0 :ligands [5 [:h 0] 2 [:h 1]] :attrs "Ct*"}]
        }"#
        );
        let mut known = unknown.clone();
        known
            .stereo_bond_mut(StereoBondId(0))
            .attributes_mut()
            .configuration =
            StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Lit(1));
        let molecule = if known_first {
            known.combine(&unknown)
        } else {
            unknown.combine(&known)
        };
        let original = molecule.clone();
        assert_eq!(
            StereoConformanceValidator::new(&StereoModel {
                stereo_bond_min_ring_size: 0,
                stereo_bond_e_min_ring_size: 8,
                ..Default::default()
            })
            .validate(&molecule),
            Ok(Solution::Contradictory(
                StereoConformanceContradiction::RingTooSmallForE {
                    bond: BondId(if known_first { 0 } else { 6 }),
                    ring_size: 6,
                    minimum_ring_size: 8,
                }
            )),
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::acyclic(BUTENE, false, Solution::Determined(()))]
    #[case::unresolved_virtual(
        r#"{:atoms ["C#h#n" "C#h#n" "C#h2" "C#h2" "C#h2" "C#h2"]
            :bonds [[0 1 "2"] [1 2 "1"] [2 3 "1"] [3 4 "1"] [4 5 "1"] [5 0 "1"]]
            :stereo-bonds [{:site 0 :ligands [5 [:virt 0] 2 [:virt 1]] :attrs "Ct1"}]}"#,
        true,
        Solution::Contradictory(StereoConformanceContradiction::RingTooSmallForE {
            bond: BondId(0), ring_size: 6, minimum_ring_size: 64,
        }),
    )]
    fn test_stereo_conformance_validator_validate_ring_frame(
        #[case] dsl: &str,
        #[case] unresolved: bool,
        #[case] expected: Solution<(), StereoConformanceContradiction>,
    ) {
        let mut molecule = mol_dsl_concrete!(dsl);
        if unresolved {
            for id in [AtomId(0), AtomId(1)] {
                let mut atom = molecule.atom_mut(id);
                atom.attributes_mut().implicit_hydrogens = NumForm::Undetermined;
                atom.attributes_mut().lone_pairs = NumForm::Undetermined;
            }
        }
        let original = molecule.clone();
        assert_eq!(
            StereoConformanceValidator::new(&StereoModel {
                stereo_bond_min_ring_size: 0,
                stereo_bond_e_min_ring_size: 64,
                ..Default::default()
            })
            .validate(&molecule),
            Ok(expected),
        );
        assert_eq!(molecule, original);
    }

    #[rstest]
    #[case::atom_configuration(
        CFCLBRI,
        (|molecule: &mut Molecule| {
            molecule.stereo_atom_mut(StereoAtomId(0)).attributes_mut().configuration =
                StereoConfigurationForm::Undetermined;
        }) as fn(&mut Molecule)
    )]
    #[case::bond_configuration(
        BUTENE,
        (|molecule: &mut Molecule| {
            molecule.stereo_bond_mut(StereoBondId(0)).attributes_mut().configuration =
                StereoConfigurationForm::Undetermined;
        }) as fn(&mut Molecule)
    )]
    #[case::typed_atom_frame(
        CFCLBRI,
        (|molecule: &mut Molecule| {
            molecule.atom_mut(AtomId(0)).attributes_mut().implicit_hydrogens =
                NumForm::var("hydrogens");
        }) as fn(&mut Molecule)
    )]
    #[case::typed_bond_frame(
        BUTENE,
        (|molecule: &mut Molecule| {
            molecule.atom_mut(AtomId(0)).attributes_mut().implicit_hydrogens =
                NumForm::var("hydrogens");
        }) as fn(&mut Molecule)
    )]
    #[case::virtual_atom_frame(
        r#"{:atoms ["C" "F" "Cl" "Br"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 [:virt 0]] :attrs "Th1"}]}"#,
        (|molecule: &mut Molecule| {
            let mut atom = molecule.atom_mut(AtomId(0));
            atom.attributes_mut().implicit_hydrogens = NumForm::Undetermined;
            atom.attributes_mut().lone_pairs = NumForm::Undetermined;
        }) as fn(&mut Molecule)
    )]
    #[case::virtual_bond_frame(
        r#"{:atoms ["C" "C" "F" "Cl"]
            :bonds [[0 1 "2"] [0 2 "1"] [1 3 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 [:virt 0] 3 [:virt 1]] :attrs "Ct1"}]}"#,
        (|molecule: &mut Molecule| {
            for id in [AtomId(0), AtomId(1)] {
                let mut atom = molecule.atom_mut(id);
                atom.attributes_mut().implicit_hydrogens = NumForm::Undetermined;
                atom.attributes_mut().lone_pairs = NumForm::Undetermined;
            }
        }) as fn(&mut Molecule)
    )]
    #[case::atom_undetermined_coset(
        CFCLBRI,
        (|molecule: &mut Molecule| {
            let mut atom = molecule.stereo_atom_mut(StereoAtomId(0));
            atom.attributes_mut().configuration =
                StereoConfigurationForm::kinded(StereoKind::Tetrahedral, StereoCoset::Undetermined);
            atom.attributes_mut().constraints.set(StereoAtomConstraintForm::Stereogenicity(
                StereogenicityForm::Lit(Stereogenicity::Stereogenic),
            ));
        }) as fn(&mut Molecule)
    )]
    #[case::bond_undetermined_coset(
        BUTENE,
        (|molecule: &mut Molecule| {
            let mut bond = molecule.stereo_bond_mut(StereoBondId(0));
            bond.attributes_mut().configuration =
                StereoConfigurationForm::kinded(StereoKind::CisTrans, StereoCoset::Undetermined);
            bond.attributes_mut().constraints.set(StereoBondConstraintForm::Stereogenicity(
                StereogenicityForm::Lit(Stereogenicity::Stereogenic),
            ));
        }) as fn(&mut Molecule)
    )]
    #[case::atom_topicity(
        CFCLBRI,
        (|molecule: &mut Molecule| {
            molecule.stereo_atom_mut(StereoAtomId(0)).attributes_mut().constraints.set(
                StereoAtomConstraintForm::Topicity(TopicityForm {
                    pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
                    relation: TopicityRelationForm::LitSet(BTreeSet::from([
                        Topicity::Homotopic, Topicity::Diastereotopic,
                    ])),
                }),
            );
        }) as fn(&mut Molecule)
    )]
    #[case::bond_topicity(
        BUTENE,
        (|molecule: &mut Molecule| {
            molecule.stereo_bond_mut(StereoBondId(0)).attributes_mut().constraints.set(
                StereoBondConstraintForm::Topicity(TopicityForm {
                    pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
                    relation: TopicityRelationForm::LitSet(BTreeSet::from([
                        Topicity::Homotopic, Topicity::Diastereotopic,
                    ])),
                }),
            );
        }) as fn(&mut Molecule)
    )]
    fn test_stereo_conformance_validator_validate_underdetermined(
        #[case] dsl: &str,
        #[case] mutate: fn(&mut Molecule),
    ) {
        let mut molecule = mol_dsl_concrete!(dsl);
        mutate(&mut molecule);
        assert_eq!(
            StereoConformanceValidator::new(&StereoModel::default()).validate(&molecule),
            Ok(Solution::Underdetermined(())),
        );
    }

    #[rstest]
    #[case::improper_on_achiral(
        BUTENE,
        (|molecule: &mut Molecule| {
            molecule
                .stereo_bond_mut(StereoBondId(0))
                .attributes_mut()
                .constraints
                .set(StereoBondConstraintForm::Stereogenicity(
                    StereogenicityForm::Lit(Stereogenicity::Prochiral),
                ));
        }) as fn(&mut Molecule),
        StereoConformanceContradiction::StereogenicityMismatch {
            asserted: StereogenicityForm::Lit(Stereogenicity::Prochiral),
            derived: Stereogenicity::Stereogenic,
        }
    )]
    #[case::atom_extra_hydrogen(
        r#"{:atoms ["C#h#T1" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#,
        (|_: &mut Molecule| {}) as fn(&mut Molecule),
        StereoConformanceContradiction::Inconsistency(
            StereoInconsistency::TetrahedralStereoFailure { atom: AtomId(0) },
        ),
    )]
    #[case::bond_extra_neighbor(
        r#"{:atoms ["C" "C" "F" "Cl" "Br" "I" "H"]
            :bonds [[0 1 "2#C1"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"] [1 6 "1"]]
            :stereo-bonds [{:site 0 :ligands [2 3 4 5] :attrs "Ct1"}]}"#,
        (|_: &mut Molecule| {}) as fn(&mut Molecule),
        StereoConformanceContradiction::Inconsistency(
            StereoInconsistency::CisTransStereoFailure { bond: BondId(0) },
        ),
    )]
    #[case::stereogenicity_mismatch(
        CFCLBRI,
        (|molecule: &mut Molecule| {
            molecule
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .constraints
                .set(StereoAtomConstraintForm::Stereogenicity(
                    StereogenicityForm::Lit(Stereogenicity::Symmetric),
                ));
        }) as fn(&mut Molecule),
        StereoConformanceContradiction::StereogenicityMismatch {
            asserted: StereogenicityForm::Lit(Stereogenicity::Symmetric),
            derived: Stereogenicity::Stereogenic,
        }
    )]
    #[case::topicity_mismatch(
        CFCLBRI,
        (|molecule: &mut Molecule| {
            molecule
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .constraints
                .set(StereoAtomConstraintForm::Topicity(TopicityForm {
                    pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
                    relation: TopicityRelationForm::Lit(Topicity::Homotopic),
                }));
        }) as fn(&mut Molecule),
        StereoConformanceContradiction::TopicityMismatch {
            pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
            asserted: TopicityRelationForm::Lit(Topicity::Homotopic),
            derived: Topicity::Diastereotopic,
        }
    )]
    #[case::topicity_complement_mismatch(
        CFCLBRI,
        (|molecule: &mut Molecule| {
            molecule.stereo_atom_mut(StereoAtomId(0)).attributes_mut().constraints.set(
                StereoAtomConstraintForm::Topicity(TopicityForm {
                    pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
                    relation: TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])),
                }),
            );
        }) as fn(&mut Molecule),
        StereoConformanceContradiction::TopicityMismatch {
            pair: StereoLigandPair::new(StereoLigandPosition(0), StereoLigandPosition(1)),
            asserted: TopicityRelationForm::NotSet(BTreeSet::from([Topicity::Diastereotopic])),
            derived: Topicity::Diastereotopic,
        }
    )]
    #[case::ligand_symmetry_violation(
        CFCLBRI,
        (|molecule: &mut Molecule| {
            molecule
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .constraints
                .set(StereoAtomConstraintForm::LigandSymmetry(LigandSymmetryForm {
                    permutation: OrientedLigandPermutation {
                        permutation: LigandPermutation(Permutation::from_image(&[1, 0, 2, 3])),
                        orientation: Orientation::Proper,
                    },
                    invariant: BooleanForm::Lit(true),
                }));
        }) as fn(&mut Molecule),
        StereoConformanceContradiction::LigandSymmetryViolation {
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
        #[case] expected: StereoConformanceContradiction,
    ) {
        let mut molecule = mol_dsl_concrete!(dsl);
        mutate(&mut molecule);
        let solution = StereoConformanceValidator::new(&StereoModel::default())
            .validate(&molecule)
            .unwrap();
        assert_eq!(solution, Solution::Contradictory(expected));
    }

    #[rstest]
    #[case::configuration_mismatch(
        r#"{:atoms ["C#T!" "F" "Cl" "Br" "I"]
            :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs "Th1"}]}"#,
        (|_: &mut Molecule| {}) as fn(&mut Molecule),
        StereoConformanceContradiction::Inconsistency(
            StereoInconsistency::TetrahedralStereoMismatch {
                atom: AtomId(5), stereo_atom: StereoAtomId(1),
            },
        ),
    )]
    #[case::stereogenicity_mismatch(
        CFCLBRI,
        (|molecule: &mut Molecule| {
            molecule.stereo_atom_mut(StereoAtomId(1)).attributes_mut().constraints.set(
                StereoAtomConstraintForm::Stereogenicity(
                    StereogenicityForm::Lit(Stereogenicity::Symmetric),
                ),
            );
        }) as fn(&mut Molecule),
        StereoConformanceContradiction::StereogenicityMismatch {
            asserted: StereogenicityForm::Lit(Stereogenicity::Symmetric),
            derived: Stereogenicity::Stereogenic,
        },
    )]
    fn test_stereo_conformance_validator_validate_partial(
        #[case] known: &str,
        #[case] mutate: fn(&mut Molecule),
        #[case] expected: StereoConformanceContradiction,
    ) {
        let mut molecule = mol_dsl_concrete!(CFCLBRI);
        molecule
            .atom_mut(AtomId(0))
            .attributes_mut()
            .implicit_hydrogens = NumForm::var("hydrogens");
        molecule.combine_from(&mol_dsl_concrete!(known));
        mutate(&mut molecule);
        assert_eq!(
            StereoConformanceValidator::new(&StereoModel::default()).validate(&molecule),
            Ok(Solution::Contradictory(expected))
        );
    }

    #[rstest]
    #[case::last_position(0, 3)]
    #[case::second_out_of_range(0, 4)]
    #[case::both_out_of_range(4, 5)]
    fn test_stereo_conformance_validator_validate_symmetry_undetermined(
        #[case] first: u32,
        #[case] second: u32,
    ) {
        let molecule = mol_dsl_concrete!(CFCLBRI);
        let graph_symmetry = molecule.graph_symmetry(&GraphSymmetryConfig {
            coloring: ConstitutionColoring::full(),
            iterate_to_fixpoint: false,
            max_iterations: 1,
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        });
        let symmetry = molecule.stereo_atom_symmetry(&graph_symmetry, StereoAtomId(0));
        let topicity = TopicityForm {
            pair: StereoLigandPair::new(StereoLigandPosition(first), StereoLigandPosition(second)),
            relation: TopicityRelationForm::Undetermined,
        };

        assert_eq!(
            StereoConformanceValidator::new(&StereoModel::default()).validate_symmetry(
                &symmetry,
                &StereogenicityForm::Undetermined,
                [topicity].iter(),
                [].iter(),
            ),
            Solution::Determined(()),
        );
    }

    #[rstest]
    #[case::second_out_of_range(0, 4)]
    #[case::both_out_of_range(4, 5)]
    fn test_stereo_conformance_validator_validate_symmetry_topicity_error(
        #[case] first: u32,
        #[case] second: u32,
    ) {
        let molecule = mol_dsl_concrete!(CFCLBRI);
        let graph_symmetry = molecule.graph_symmetry(&GraphSymmetryConfig {
            coloring: ConstitutionColoring::full(),
            iterate_to_fixpoint: false,
            max_iterations: 1,
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        });
        let symmetry = molecule.stereo_atom_symmetry(&graph_symmetry, StereoAtomId(0));
        let pair = StereoLigandPair::new(StereoLigandPosition(first), StereoLigandPosition(second));
        let topicity = TopicityForm {
            pair,
            relation: TopicityRelationForm::Lit(Topicity::Diastereotopic),
        };

        assert_eq!(
            StereoConformanceValidator::new(&StereoModel::default()).validate_symmetry(
                &symmetry,
                &StereogenicityForm::Undetermined,
                [topicity].iter(),
                [].iter(),
            ),
            Solution::Contradictory(StereoConformanceContradiction::TopicityPositionOutOfRange {
                pair,
                degree: 4,
            }),
        );
    }

    #[rstest]
    #[case::short_positive(3, BooleanForm::Lit(true))]
    #[case::short_negative(3, BooleanForm::Lit(false))]
    #[case::short_undetermined(3, BooleanForm::Undetermined)]
    #[case::long_positive(5, BooleanForm::Lit(true))]
    #[case::long_negative(5, BooleanForm::Lit(false))]
    #[case::long_undetermined(5, BooleanForm::Undetermined)]
    fn test_stereo_conformance_validator_validate_symmetry_ligand_symmetry_error(
        #[case] degree: usize,
        #[case] invariant: BooleanForm,
    ) {
        let molecule = mol_dsl_concrete!(CFCLBRI);
        let graph_symmetry = molecule.graph_symmetry(&GraphSymmetryConfig {
            coloring: ConstitutionColoring::full(),
            iterate_to_fixpoint: false,
            max_iterations: 1,
            automorphism_algorithm: AutomorphismAlgorithm::Nauty,
        });
        let symmetry = molecule.stereo_atom_symmetry(&graph_symmetry, StereoAtomId(0));
        let asserted = LigandSymmetryForm {
            permutation: OrientedLigandPermutation {
                permutation: LigandPermutation(Permutation::identity(degree)),
                orientation: Orientation::Proper,
            },
            invariant,
        };

        assert_eq!(
            StereoConformanceValidator::new(&StereoModel::default()).validate_symmetry(
                &symmetry,
                &StereogenicityForm::Undetermined,
                [].iter(),
                [asserted].iter(),
            ),
            Solution::Contradictory(StereoConformanceContradiction::LigandSymmetryViolation {
                asserted,
            }),
        );
    }
    #[rstest]
    #[case::second_position(0, 4)]
    #[case::both_positions(4, 5)]
    fn test_stereo_conformance_validator_validate_topicity_error(
        #[case] first: u32,
        #[case] second: u32,
        #[values(false, true)] bond: bool,
    ) {
        let mut molecule = mol_dsl_concrete!(if bond { BUTENE } else { CFCLBRI });
        let pair = StereoLigandPair::new(StereoLigandPosition(first), StereoLigandPosition(second));
        let assertion = TopicityForm {
            pair,
            relation: TopicityRelationForm::Lit(Topicity::Diastereotopic),
        };
        if bond {
            molecule
                .stereo_bond_mut(StereoBondId(0))
                .attributes_mut()
                .constraints
                .set(StereoBondConstraintForm::Topicity(assertion));
        } else {
            molecule
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .constraints
                .set(StereoAtomConstraintForm::Topicity(assertion));
        }
        assert_eq!(
            StereoConformanceValidator::new(&StereoModel::default()).validate(&molecule),
            Ok(Solution::Contradictory(
                StereoConformanceContradiction::TopicityPositionOutOfRange { pair, degree: 4 }
            ))
        );
    }

    #[rstest]
    fn test_stereo_conformance_validator_validate_ligand_symmetry_error(
        #[values(false, true)] invariant: bool,
        #[values(false, true)] bond: bool,
    ) {
        let mut molecule = mol_dsl_concrete!(if bond { BUTENE } else { CFCLBRI });
        let asserted = LigandSymmetryForm {
            permutation: OrientedLigandPermutation {
                permutation: LigandPermutation(Permutation::identity(3)),
                orientation: Orientation::Proper,
            },
            invariant: BooleanForm::Lit(invariant),
        };
        if bond {
            molecule
                .stereo_bond_mut(StereoBondId(0))
                .attributes_mut()
                .constraints
                .set(StereoBondConstraintForm::LigandSymmetry(asserted));
        } else {
            molecule
                .stereo_atom_mut(StereoAtomId(0))
                .attributes_mut()
                .constraints
                .set(StereoAtomConstraintForm::LigandSymmetry(asserted));
        }
        assert_eq!(
            StereoConformanceValidator::new(&StereoModel::default()).validate(&molecule),
            Ok(Solution::Contradictory(
                StereoConformanceContradiction::LigandSymmetryViolation { asserted }
            ))
        );
    }
}
