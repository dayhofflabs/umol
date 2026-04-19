//! Validation engine: checks that a `MoleculeAst` satisfies the chemistry and
//! any attached constraints.

use crate::api::molecule::Molecule;
use crate::ast::molecule::MoleculeAst;
use crate::ops::chemistry::Chemistry;
use crate::ops::propagate::ElectronInvariant;
use crate::ops::solution::Solution;

pub struct Validator<'s> {
    chemistry: &'s Chemistry,
}

impl<'s> Validator<'s> {
    pub fn new(chemistry: &'s Chemistry) -> Self {
        Self { chemistry }
    }

    pub fn validate(&self, ast: &MoleculeAst) -> Solution<()> {
        let electron_invariant = ElectronInvariant;
        for i in 0..ast.atoms().count() {
            if !electron_invariant.validate(ast, i) {
                return Solution::Contradictory;
            }
            if !self.chemistry.valence.validate(ast, i) {
                return Solution::Contradictory;
            }
        }
        if !ast.atoms().iter().all(|v| v.data.is_ground()) {
            return Solution::Underdetermined(());
        }
        let Ok(mol) = Molecule::new(ast.clone()) else {
            return Solution::Underdetermined(());
        };
        for c in mol.ast().constraints().iter() {
            if !c.evaluate(&mol) {
                return Solution::Contradictory;
            }
        }
        Solution::Determined(())
    }
}
