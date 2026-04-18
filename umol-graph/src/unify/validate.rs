//! Validation engine: checks that a `MoleculeAst` satisfies the chemistry.

use crate::ast::molecule::MoleculeAst;
use crate::unify::chemistry::Chemistry;
use crate::unify::error::ValidationError;

pub struct Validator<'s> {
    chemistry: &'s Chemistry,
}

impl<'s> Validator<'s> {
    pub fn new(chemistry: &'s Chemistry) -> Self {
        Self { chemistry }
    }

    pub fn validate(&self, ast: &MoleculeAst) -> Result<(), ValidationError> {
        for i in 0..ast.atoms().count() {
            if !self.chemistry.valence.valence_balance(ast, i) {
                return Err(ValidationError::Contradictory);
            }
        }
        if ast.atoms().iter().all(|v| v.data.is_ground()) {
            Ok(())
        } else {
            Err(ValidationError::Underdetermined)
        }
    }
}
