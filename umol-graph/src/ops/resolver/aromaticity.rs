//! Aromaticity resolver. Wraps [`AromaticityPerception`] with the
//! resolver-flow closure and Solution-typed return shape. Input contract:
//! AST atoms carry `AromaticValence::Aromatic(_)` hints (filled in by
//! atom-typing); the resolver runs perception against those, validates
//! them, and writes back aromatic system entries on success.

use umol_ast::ast::MoleculeAst;

use crate::ops::aromaticity::{
    electrons_from_aromatic_constraint, AromaticityContradiction, AromaticityError,
    AromaticityPerception,
};
use crate::ops::config::AromaticityModel;
use crate::ops::solution::Solution;

#[derive(Clone, Debug)]
pub struct AromaticityResolver {
    perception: AromaticityPerception,
}

impl AromaticityResolver {
    pub fn new(model: &AromaticityModel) -> Self {
        Self {
            perception: AromaticityPerception::new(model),
        }
    }

    pub fn resolve(
        &self,
        ast: &mut MoleculeAst,
    ) -> Result<Solution<(), AromaticityContradiction>, AromaticityError> {
        let outcome = self
            .perception
            .find_systems(ast, electrons_from_aromatic_constraint)?;
        match outcome {
            Solution::Determined(systems) => {
                self.perception.add_systems(ast, systems);
                Ok(Solution::Determined(()))
            }
            Solution::Underdetermined(_) => Ok(Solution::Underdetermined(())),
            Solution::Contradictory(c) => Ok(Solution::Contradictory(c)),
        }
    }
}
