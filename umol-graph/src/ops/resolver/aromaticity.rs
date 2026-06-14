//! Aromaticity resolver. Wraps [`AromaticityPerception`] with the
//! resolver-flow closure and Solution-typed return shape. Input contract:
//! AST atoms carry `AromaticValence::Aromatic(_)` hints (filled in by
//! atom-typing); the resolver runs perception against those, validates
//! them, and writes back aromatic system entries on success.

use umol_ast::ast::{AromaticValenceAst, MoleculeAst, ValueAst};

use crate::ops::aromaticity::{AromaticityContradiction, AromaticityError, AromaticityPerception};
use crate::ops::model::AromaticityModel;
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
        let outcome = self.perception.find_systems(ast, |v| {
            // Skip atoms already in an aromatic system so re-runs don't re-add it.
            if v.is_in_aromatic_system() {
                return None;
            }
            match v.ast.constraints.aromatic_valence() {
                AromaticValenceAst::Aromatic(ValueAst::Lit(n)) if n >= 0 => Some(n as u8),
                _ => None,
            }
        })?;
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
