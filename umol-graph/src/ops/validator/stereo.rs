//! Stereo validator.

use thiserror::Error;
use umol_ast::ast::{
    AsLit, GraphSymmetry, Lattice, StereoLigandPair, LigandSymmetryAst, MemOp, MoleculeAst,
    StereoAtomId, StereoBondId, StereoKind, StereoSymmetry, Stereogenicity, StereogenicityAst,
    Topicity, TopicityAst, TopicityRelationAst,
};
use umol_perm::OrientedPermutation;

use crate::ops::model::StereoModel;
use crate::ops::solution::Solution;

#[derive(Clone, Debug)]
pub struct StereoValidator {
    model: StereoModel,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum StereoValidatorContradiction {
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
        asserted: StereogenicityAst,
        derived: Stereogenicity,
    },
    #[error("ligand pair {pair:?} topicity {derived:?} contradicts asserted {asserted:?}")]
    TopicityMismatch {
        pair: StereoLigandPair,
        asserted: TopicityRelationAst,
        derived: Topicity,
    },
    #[error("asserted ligand symmetry {asserted:?} not satisfied by the derived group")]
    LigandSymmetryViolation { asserted: LigandSymmetryAst },
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum StereoValidatorError {}

impl StereoValidator {
    pub fn new(model: &StereoModel) -> Self {
        Self {
            model: model.clone(),
        }
    }

    /// Validate every stereo element against the molecule's graph symmetry.
    /// First `Contradictory` wins; a ground assertion the derived value leaves
    /// open contributes `Underdetermined`.
    pub fn validate(
        &self,
        ast: &MoleculeAst,
    ) -> Result<Solution<(), StereoValidatorContradiction>, StereoValidatorError> {
        let symmetry = ast.graph_symmetry(&self.model.graph_symmetry_config());
        let mut any_undetermined = false;

        for id in ast.stereo_atoms().ids() {
            match self.validate_stereo_atom(ast, id, &symmetry) {
                Solution::Determined(()) => {}
                Solution::Underdetermined(()) => any_undetermined = true,
                Solution::Contradictory(c) => return Ok(Solution::Contradictory(c)),
            }
        }
        for id in ast.stereo_bonds().ids() {
            match self.validate_stereo_bond(ast, id, &symmetry) {
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
        ast: &MoleculeAst,
        id: StereoAtomId,
        symmetry: &GraphSymmetry,
    ) -> Solution<(), StereoValidatorContradiction> {
        let view = ast.stereo_atom(id);
        let constraints = view.constraints();
        let stereogenicity = constraints.stereogenicity();
        let topicities: Vec<TopicityAst> = constraints.topicities().cloned().collect();
        let ligand_symmetries: Vec<LigandSymmetryAst> =
            constraints.ligand_symmetry().copied().collect();
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
        let sym = ast.stereo_atom_symmetry(symmetry, id);
        self.validate_symmetry(&sym, &stereogenicity, &topicities, &ligand_symmetries)
    }

    fn validate_stereo_bond(
        &self,
        ast: &MoleculeAst,
        id: StereoBondId,
        symmetry: &GraphSymmetry,
    ) -> Solution<(), StereoValidatorContradiction> {
        let view = ast.stereo_bond(id);
        let constraints = view.constraints();
        let stereogenicity = constraints.stereogenicity();
        let topicities: Vec<TopicityAst> = constraints.topicities().cloned().collect();
        let ligand_symmetries: Vec<LigandSymmetryAst> =
            constraints.ligand_symmetry().copied().collect();
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
        let sym = ast.stereo_bond_symmetry(symmetry, id);
        self.validate_symmetry(&sym, &stereogenicity, &topicities, &ligand_symmetries)
    }

    /// Ligand arity (= kind degree), the achiral-kind relation gate, and coset
    /// range; an absent coset leaves the element underdetermined.
    fn validate_kind(
        &self,
        kind: StereoKind,
        coset: Option<u32>,
        ligand_count: usize,
        stereogenicity: &StereogenicityAst,
        topicities: &[TopicityAst],
    ) -> Solution<(), StereoValidatorContradiction> {
        if ligand_count != kind.degree() {
            return Solution::Contradictory(StereoValidatorContradiction::LigandArity {
                kind,
                ligands: ligand_count,
                degree: kind.degree(),
            });
        }
        let has_improper = (!stereogenicity.is_undetermined()
            && stereogenicity.matches(&StereogenicityAst::Lit(Stereogenicity::Prochiral)))
            || topicities.iter().any(|t| {
                !t.relation.is_undetermined()
                    && t.relation
                        .matches(&TopicityRelationAst::Lit(Topicity::Enantiotopic))
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
        stereogenicity: &StereogenicityAst,
        topicities: &[TopicityAst],
        ligand_symmetries: &[LigandSymmetryAst],
    ) -> Solution<(), StereoValidatorContradiction> {
        let mut any_undetermined = false;

        if !stereogenicity.is_undetermined() {
            let derived = sym.stereogenicity();
            if !stereogenicity.matches(&StereogenicityAst::Lit(derived)) {
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
            if !t.relation.matches(&TopicityRelationAst::Lit(derived)) {
                return Solution::Contradictory(StereoValidatorContradiction::TopicityMismatch {
                    pair: t.pair,
                    asserted: t.relation.clone(),
                    derived,
                });
            }
            any_undetermined |= !t.relation.is_ground();
        }

        for ls in ligand_symmetries {
            let op = OrientedPermutation::new(ls.permutation.permutation.0, ls.permutation.orientation);
            let in_group = sym.group().contains(op);
            let holds = match ls.member {
                MemOp::In => in_group,
                MemOp::NotIn => !in_group,
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
    use umol_ast::ast::{
        OrientedLigandPermutation, LigandPermutation, StereoAtomConstraint, StereoBondConstraint,
        StereoConfigurationAst, StereoKind, StereoLigandId, StereogenicityAst,
    };
    use umol_ast::mol_ground;
    use umol_perm::{Orientation, Permutation};

    use super::*;

    // Tetrahedral stereocenter, four distinct halide ligands → trivial ligand
    // symmetry, genuinely stereogenic, every ligand pair diastereotopic.
    const CFCLBRI: &str = r#"{:atoms ["C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]
        :stereo-atoms [{:site 0 :ligands [1 2 3 4] :type "Th1"}]}"#;

    // Cis/trans double bond, each terminus bearing two distinct substituents.
    const BUTENE: &str = r#"{:atoms ["C" "C" "F" "Cl" "Br" "I"]
        :bonds [[0 1 "2"] [0 2 "1"] [0 3 "1"] [1 4 "1"] [1 5 "1"]]
        :stereo-bonds [{:site 0 :ligands [2 3 4 5] :type "Ct1"}]}"#;

    #[rstest]
    #[case::bare(CFCLBRI, (|_: &mut MoleculeAst| {}) as fn(&mut MoleculeAst))]
    #[case::matching_stereogenicity(
        CFCLBRI,
        (|ast: &mut MoleculeAst| {
            ast.stereo_atom_mut(StereoAtomId(0))
                .constraints
                .add(StereoAtomConstraint::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Stereogenic)));
        }) as fn(&mut MoleculeAst)
    )]
    fn test_stereo_validator_validate(#[case] dsl: &str, #[case] mutate: fn(&mut MoleculeAst)) {
        let mut ast = mol_ground!(dsl);
        mutate(&mut ast);
        let solution = StereoValidator::new(&StereoModel::default())
            .validate(&ast)
            .unwrap();
        assert_eq!(solution, Solution::Determined(()));
    }

    #[rstest]
    #[case::coset_out_of_range(
        CFCLBRI,
        (|ast: &mut MoleculeAst| {
            ast.stereo_atom_mut(StereoAtomId(0)).configuration = StereoConfigurationAst::kinded(StereoKind::Tetrahedral, 9);
        }) as fn(&mut MoleculeAst),
        StereoValidatorContradiction::CosetOutOfRange {
            kind: StereoKind::Tetrahedral,
            coset: 9,
            count: 2,
        }
    )]
    #[case::arity(
        CFCLBRI,
        (|ast: &mut MoleculeAst| {
            ast.stereo_atom_mut(StereoAtomId(0)).configuration = StereoConfigurationAst::kinded(StereoKind::TrigonalBipyramidal, 0);
        }) as fn(&mut MoleculeAst),
        StereoValidatorContradiction::LigandArity {
            kind: StereoKind::TrigonalBipyramidal,
            ligands: 4,
            degree: 5,
        }
    )]
    #[case::improper_on_achiral(
        BUTENE,
        (|ast: &mut MoleculeAst| {
            ast.stereo_bond_mut(StereoBondId(0))
                .constraints
                .add(StereoBondConstraint::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Prochiral)));
        }) as fn(&mut MoleculeAst),
        StereoValidatorContradiction::ImproperOnAchiral {
            kind: StereoKind::CisTrans,
        }
    )]
    #[case::stereogenicity_mismatch(
        CFCLBRI,
        (|ast: &mut MoleculeAst| {
            ast.stereo_atom_mut(StereoAtomId(0))
                .constraints
                .add(StereoAtomConstraint::Stereogenicity(StereogenicityAst::Lit(Stereogenicity::Symmetric)));
        }) as fn(&mut MoleculeAst),
        StereoValidatorContradiction::StereogenicityMismatch {
            asserted: StereogenicityAst::Lit(Stereogenicity::Symmetric),
            derived: Stereogenicity::Stereogenic,
        }
    )]
    #[case::topicity_mismatch(
        CFCLBRI,
        (|ast: &mut MoleculeAst| {
            ast.stereo_atom_mut(StereoAtomId(0))
                .constraints
                .add(StereoAtomConstraint::Topicity(TopicityAst {
                    pair: StereoLigandPair::new(StereoLigandId(0), StereoLigandId(1)),
                    relation: TopicityRelationAst::Lit(Topicity::Homotopic),
                }));
        }) as fn(&mut MoleculeAst),
        StereoValidatorContradiction::TopicityMismatch {
            pair: StereoLigandPair::new(StereoLigandId(0), StereoLigandId(1)),
            asserted: TopicityRelationAst::Lit(Topicity::Homotopic),
            derived: Topicity::Diastereotopic,
        }
    )]
    #[case::ligand_symmetry_violation(
        CFCLBRI,
        (|ast: &mut MoleculeAst| {
            ast.stereo_atom_mut(StereoAtomId(0))
                .constraints
                .add(StereoAtomConstraint::LigandSymmetry(LigandSymmetryAst {
                    permutation: OrientedLigandPermutation {
                        permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])),
                        orientation: Orientation::Proper,
                    },
                    member: MemOp::In,
                }));
        }) as fn(&mut MoleculeAst),
        StereoValidatorContradiction::LigandSymmetryViolation {
            asserted: LigandSymmetryAst {
                permutation: OrientedLigandPermutation {
                    permutation: LigandPermutation(Permutation::from_image(4, &[1, 0, 2, 3])),
                    orientation: Orientation::Proper,
                },
                member: MemOp::In,
            },
        }
    )]
    fn test_stereo_validator_validate_error(
        #[case] dsl: &str,
        #[case] mutate: fn(&mut MoleculeAst),
        #[case] expected: StereoValidatorContradiction,
    ) {
        let mut ast = mol_ground!(dsl);
        mutate(&mut ast);
        let solution = StereoValidator::new(&StereoModel::default())
            .validate(&ast)
            .unwrap();
        assert_eq!(solution, Solution::Contradictory(expected));
    }
}
