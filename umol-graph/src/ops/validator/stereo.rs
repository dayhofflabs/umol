//! Stereo validator.

use thiserror::Error;
use umol_ast::ast::{
    AsLit, GraphSymmetry, Lattice, MoleculeAst, StereoAtomId, StereoBondId, StereoKind,
    Stereogenicity, StereogenicityRelationAst, Topicity, TopicityRelationAst,
};

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
        _symmetry: &GraphSymmetry,
    ) -> Solution<(), StereoValidatorContradiction> {
        let view = ast.stereo_atom(id);
        let constraints = view.constraints();
        let improper_topicity = constraints.topicities().any(|t| {
            !t.rel.is_undetermined()
                && t.rel
                    .matches(&TopicityRelationAst::Lit(Topicity::Enantiotopic))
        });
        let stereogenicity = constraints.stereogenicity();
        let improper_stereogenicity = !stereogenicity.is_undetermined()
            && stereogenicity.matches(&StereogenicityRelationAst::Lit(Stereogenicity::Prochiral));
        self.validate_kind(
            view.kind(),
            view.coset().as_lit(),
            view.ligands().count(),
            improper_topicity || improper_stereogenicity,
        )
    }

    fn validate_stereo_bond(
        &self,
        ast: &MoleculeAst,
        id: StereoBondId,
        _symmetry: &GraphSymmetry,
    ) -> Solution<(), StereoValidatorContradiction> {
        let view = ast.stereo_bond(id);
        let constraints = view.constraints();
        let improper_topicity = constraints.topicities().any(|t| {
            !t.rel.is_undetermined()
                && t.rel
                    .matches(&TopicityRelationAst::Lit(Topicity::Enantiotopic))
        });
        let stereogenicity = constraints.stereogenicity();
        let improper_stereogenicity = !stereogenicity.is_undetermined()
            && stereogenicity.matches(&StereogenicityRelationAst::Lit(Stereogenicity::Prochiral));
        self.validate_kind(
            view.kind(),
            view.coset().as_lit(),
            view.ligands().count(),
            improper_topicity || improper_stereogenicity,
        )
    }

    /// Validate the kind-dependent invariants of a stereo element:
    /// 1. Ligand arity must equal the kind's degree.
    /// 2. An achiral kind admits no improper (enantio-) relation.
    /// 3. The coset must be in range (`None` coset → underdetermined).
    fn validate_kind(
        &self,
        kind: StereoKind,
        coset: Option<u32>,
        ligand_count: usize,
        has_improper_relation: bool,
    ) -> Solution<(), StereoValidatorContradiction> {
        if ligand_count != kind.degree() {
            return Solution::Contradictory(StereoValidatorContradiction::LigandArity {
                kind,
                ligands: ligand_count,
                degree: kind.degree(),
            });
        }
        if !kind.is_chiral_class() && has_improper_relation {
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
}
