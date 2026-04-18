//! Validation engine: checks that a `MoleculeAst` satisfies the chemistry.

use crate::ast::molecule::MoleculeAst;
use crate::unify::chemistry::Chemistry;
use crate::unify::error::ValidationError;
use crate::unify::propagate::ElectronInvariant;

pub struct Validator<'s> {
    chemistry: &'s Chemistry,
}

impl<'s> Validator<'s> {
    pub fn new(chemistry: &'s Chemistry) -> Self {
        Self { chemistry }
    }

    pub fn validate(&self, ast: &MoleculeAst) -> Result<(), ValidationError> {
        let electron_invariant = ElectronInvariant;
        for i in 0..ast.atoms().count() {
            if !electron_invariant.validate(ast, i) {
                return Err(ValidationError::Contradictory);
            }
            if !self.chemistry.valence.validate(ast, i) {
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
