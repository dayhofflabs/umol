//! Validation engine: checks that a `MoleculeAst` satisfies the chemistry and
//! any attached constraints.

use crate::api::molecule::Molecule;
use crate::ast::molecule::MoleculeAst;
use crate::ops::chemistry::Chemistry;
use crate::ops::error::ValidationError;
use crate::ops::propagate::ElectronInvariant;

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
        if !ast.atoms().iter().all(|v| v.data.is_ground()) {
            return Err(ValidationError::Underdetermined);
        }
        let mol =
            Molecule::new(ast.clone()).map_err(|_| ValidationError::Underdetermined)?;
        for c in mol.ast().constraints().iter() {
            if !c.evaluate(&mol) {
                return Err(ValidationError::Contradictory);
            }
        }
        Ok(())
    }
}
