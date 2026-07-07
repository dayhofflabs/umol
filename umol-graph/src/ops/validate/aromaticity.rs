//! Aromaticity validator. Wraps [`AromaticityPerception`] and verifies that
//! the aromatic systems already in the AST agree with what perception
//! independently finds. Input contract: AST atoms carry filled-in
//! `AromaticValence::Aromatic(Lit(n))` (atom-typing has run) and the AST
//! already carries one or more `AromaticSystemAst` entries.

use thiserror::Error;
use umol_ast::ast::{AromaticValenceAst, AtomId, MoleculeAst, ValueAst};
use umol_utils::solution::Solution;

use crate::ops::aromaticity::{AromaticityContradiction, AromaticityError, AromaticityPerception};
use crate::ops::model::AromaticityModel;

#[derive(Clone, Debug)]
pub struct AromaticityConformanceValidator {
    perception: AromaticityPerception,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum AromaticityValidatorContradiction {
    #[error("perception rejected the input: {0}")]
    Perception(AromaticityContradiction),
    #[error(
        "aromatic system count mismatch: AST has {ast_count}, perception found {perception_count}"
    )]
    SystemCountMismatch {
        ast_count: usize,
        perception_count: usize,
    },
    #[error("aromatic system atoms differ: AST {ast_atoms:?}, perception {perception_atoms:?}")]
    AtomsMismatch {
        ast_atoms: Vec<AtomId>,
        perception_atoms: Vec<AtomId>,
    },
}

impl AromaticityConformanceValidator {
    pub fn new(model: &AromaticityModel) -> Self {
        Self {
            perception: AromaticityPerception::new(model),
        }
    }

    pub fn validate(
        &self,
        ast: &MoleculeAst,
    ) -> Result<Solution<(), AromaticityValidatorContradiction>, AromaticityError> {
        let outcome = self.perception.find_systems(ast, |v| {
            match v
                .ast
                .constraints
                .aromatic_valence()
                .unwrap_or(&AromaticValenceAst::Undetermined)
            {
                AromaticValenceAst::Aromatic(ValueAst::Lit(n)) if *n >= 0 => Some(*n as u8),
                _ => None,
            }
        })?;
        let perception_systems = match outcome {
            Solution::Determined(systems) => systems,
            Solution::Underdetermined(_) => return Ok(Solution::Underdetermined(())),
            Solution::Contradictory(c) => {
                return Ok(Solution::Contradictory(
                    AromaticityValidatorContradiction::Perception(c),
                ));
            }
        };

        let ast_systems: Vec<Vec<AtomId>> = ast
            .aromatic_systems()
            .iter()
            .map(|view| {
                let mut atoms: Vec<AtomId> = view.atom_ids().collect();
                atoms.sort_unstable();
                atoms
            })
            .collect();
        let mut ast_systems_sorted = ast_systems;
        ast_systems_sorted.sort_by(|a, b| a.first().cmp(&b.first()));

        if ast_systems_sorted.len() != perception_systems.len() {
            return Ok(Solution::Contradictory(
                AromaticityValidatorContradiction::SystemCountMismatch {
                    ast_count: ast_systems_sorted.len(),
                    perception_count: perception_systems.len(),
                },
            ));
        }
        for (ast_atoms, (perception_atoms, _)) in
            ast_systems_sorted.iter().zip(perception_systems.iter())
        {
            if ast_atoms != perception_atoms {
                return Ok(Solution::Contradictory(
                    AromaticityValidatorContradiction::AtomsMismatch {
                        ast_atoms: ast_atoms.clone(),
                        perception_atoms: perception_atoms.clone(),
                    },
                ));
            }
        }
        Ok(Solution::Determined(()))
    }
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_ast::ast::{
        AromaticValenceAst, AtomAst, AtomConstraint, BondAst, MoleculeAst, SpinStateAst, ValueAst,
    };
    use umol_chem::element::Element;

    use super::*;
    use crate::ops::model::{ElementScope, RingLimits};
    use crate::ops::resolve::aromaticity::AromaticityResolver;

    #[fixture]
    fn carbon_only() -> AromaticityModel {
        AromaticityModel::HueckelRule {
            scope: ElementScope::AllowList(vec![Element::C]),
            ring_limits: RingLimits::default(),
        }
    }

    #[fixture]
    fn benzene() -> MoleculeAst {
        let atoms: Vec<AtomAst> = (0..6)
            .map(|_| {
                let mut atom = AtomAst::from_element(Element::C);
                atom.charge = ValueAst::Lit(0);
                atom.spin = SpinStateAst::closed_shell();
                atom.constraints.set(AtomConstraint::AromaticValence(
                    AromaticValenceAst::Aromatic(ValueAst::Lit(1)),
                ));
                atom
            })
            .collect();
        let bonds: Vec<_> = (0..6)
            .map(|i| (AtomId(i), AtomId((i + 1) % 6), BondAst::from_order(1)))
            .collect();
        MoleculeAst::from_atoms_and_bonds(atoms, bonds)
    }

    #[rstest]
    fn test_aromaticity_conformance_validator_validate(
        benzene: MoleculeAst,
        carbon_only: AromaticityModel,
    ) {
        let mut ast = benzene;
        AromaticityResolver::new(&carbon_only)
            .resolve(&mut ast)
            .unwrap();
        assert_eq!(
            AromaticityConformanceValidator::new(&carbon_only)
                .validate(&ast)
                .unwrap(),
            Solution::Determined(())
        );
    }

    #[rstest]
    fn test_aromaticity_conformance_validator_validate_error(
        benzene: MoleculeAst,
        carbon_only: AromaticityModel,
    ) {
        assert_eq!(
            AromaticityConformanceValidator::new(&carbon_only)
                .validate(&benzene)
                .unwrap(),
            Solution::Contradictory(AromaticityValidatorContradiction::SystemCountMismatch {
                ast_count: 0,
                perception_count: 1,
            })
        );
    }
}
