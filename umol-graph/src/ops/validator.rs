//! Validators: physics invariants (electron count, spin coupling),
//! constraint cross-checks, entity-structure shape checks, and aromaticity
//! verification.
//!
//! Each validator borrows a `MoleculeAst` (or `AtomAst`) and returns
//! `Result<Solution<(), C>, E>`. Determined and Underdetermined are both
//! successful outcomes; only `Contradictory(C)` is a failure on the `Solution`
//! side. Setup-level failures (parameter-table gaps, etc.) live in `Err(E)`;
//! tier-2 validators that have no setup use uninhabited error types.
//!
//! The composite [`Validator`] runs the four sub-validators in order and lifts
//! their per-engine `Contradiction` and `Error` types into unions via `From`
//! impls. `validate_atom` runs only those sub-validators that make sense
//! without a surrounding molecule (atom-typing registry use).
//! `AromaticityValidator` is configured separately because it carries a model.

pub mod aromaticity;
pub mod constraint;
pub mod entity;
pub mod invariant;
pub mod spin;

pub use aromaticity::{AromaticityValidator, AromaticityValidatorContradiction};
pub use constraint::{ConstraintContradiction, ConstraintError, ConstraintValidator};
pub use entity::{EntityStructureContradiction, EntityStructureError, EntityStructureValidator};
pub use invariant::{
    ElectronInvariantContradiction, ElectronInvariantError, ElectronInvariantValidator,
};
pub use spin::{SpinCouplingContradiction, SpinCouplingError, SpinCouplingValidator};
use thiserror::Error;
use umol_ast::ast::{AtomAst, MoleculeAst};

use crate::ops::solution::Solution;

#[derive(Clone, Copy, Debug, Default)]
pub struct Validator {
    pub electron_invariant: ElectronInvariantValidator,
    pub spin_coupling: SpinCouplingValidator,
    pub constraint: ConstraintValidator,
    pub entity_structure: EntityStructureValidator,
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValidatorContradiction {
    #[error(transparent)]
    ElectronInvariant(#[from] ElectronInvariantContradiction),
    #[error(transparent)]
    SpinCoupling(#[from] SpinCouplingContradiction),
    #[error(transparent)]
    Constraint(#[from] ConstraintContradiction),
    #[error(transparent)]
    EntityStructure(#[from] EntityStructureContradiction),
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ValidatorError {
    #[error(transparent)]
    ElectronInvariant(#[from] ElectronInvariantError),
    #[error(transparent)]
    SpinCoupling(#[from] SpinCouplingError),
    #[error(transparent)]
    Constraint(#[from] ConstraintError),
    #[error(transparent)]
    EntityStructure(#[from] EntityStructureError),
}

impl Validator {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn validate(
        &self,
        ast: impl AsRef<MoleculeAst>,
    ) -> Result<Solution<(), ValidatorContradiction>, ValidatorError> {
        let ast = ast.as_ref();
        let mut any_undetermined = false;

        // Run validators in order. First contradiction wins.
        match self.entity_structure.validate(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        match self.electron_invariant.validate(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        match self.spin_coupling.validate(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        match self.constraint.validate(ast)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }

        Ok(if any_undetermined {
            Solution::Underdetermined(())
        } else {
            Solution::Determined(())
        })
    }

    pub fn validate_atom(
        &self,
        atom: &AtomAst,
    ) -> Result<Solution<(), ValidatorContradiction>, ValidatorError> {
        let mut any_undetermined = false;

        match self.electron_invariant.validate_atom(atom)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }
        match self.spin_coupling.validate_atom(atom)? {
            Solution::Determined(()) => {}
            Solution::Underdetermined(()) => any_undetermined = true,
            Solution::Contradictory(c) => return Ok(Solution::Contradictory(c.into())),
        }

        Ok(if any_undetermined {
            Solution::Underdetermined(())
        } else {
            Solution::Determined(())
        })
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_ast::ast::{
        AtomAst, AtomConstraint, AtomId, BondAst, Constraints, ImplicitHydrogensAst, MoleculeAst,
        SpinStateAst, ValueAst,
    };
    use umol_shared::element::Element;

    use super::*;

    fn ground_methane_atom() -> AtomAst {
        let mut atom = AtomAst::from_element(Element::C);
        atom.charge = ValueAst::Lit(0);
        atom.lone_pairs = ValueAst::Lit(0);
        atom.implicit_hydrogens = ImplicitHydrogensAst::Lit(4);
        atom.spin = SpinStateAst::from((0_u8, 1_u8));
        atom
    }

    fn ethane() -> MoleculeAst {
        let mut ch3_a = AtomAst::from_element(Element::C);
        ch3_a.charge = ValueAst::Lit(0);
        ch3_a.lone_pairs = ValueAst::Lit(0);
        ch3_a.implicit_hydrogens = ImplicitHydrogensAst::Lit(3);
        ch3_a.spin = SpinStateAst::from((0_u8, 1_u8));
        let ch3_b = ch3_a.clone();
        MoleculeAst::from_parts(
            vec![ch3_a, ch3_b],
            vec![(AtomId(0), AtomId(1), BondAst::from_order(1))],
            vec![],
            vec![],
            vec![],
            vec![],
            Constraints::default(),
        )
    }

    #[rstest]
    fn test_validator_composite_validate_determined() {
        let v = Validator::new();
        let result = v.validate(ethane()).unwrap();
        assert!(matches!(result, Solution::Determined(())));
    }

    #[rstest]
    fn test_validator_composite_validate_atom_determined() {
        let v = Validator::new();
        let atom = ground_methane_atom();
        let result = v.validate_atom(&atom).unwrap();
        assert!(matches!(result, Solution::Determined(())));
    }

    #[rstest]
    fn test_validator_composite_validate_atom_with_constraint_only() {
        let v = Validator::new();
        let mut atom = AtomAst::from_element(Element::C);
        atom.charge = ValueAst::Lit(0);
        atom.lone_pairs = ValueAst::Lit(0);
        atom.implicit_hydrogens = ImplicitHydrogensAst::Lit(3);
        atom.spin = SpinStateAst::from((0_u8, 1_u8));
        atom.constraints
            .add(AtomConstraint::Valence(ValueAst::Lit(1)));
        let result = v.validate_atom(&atom).unwrap();
        assert!(matches!(result, Solution::Determined(())));
    }
}
