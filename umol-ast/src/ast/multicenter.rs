//! Multicenter bond AST.

use super::constraint::MulticenterBondConstraint;
use super::spin::SpinStateAst;
use super::value::ValueAst;

/// Multicenter bond: structural attributes of a bond spanning three or more
/// atoms. Electron count and spin live on the system itself; individual atom
/// participation is described by `MulticenterValence` atom constraints.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MulticenterBondAst {
    pub charge: ValueAst,
    pub spin: SpinStateAst,
    pub electrons: ValueAst,
    pub constraints: Vec<MulticenterBondConstraint>,
}

impl MulticenterBondAst {
    pub fn new(charge: ValueAst, spin: SpinStateAst, electrons: ValueAst) -> Self {
        Self {
            charge,
            spin,
            electrons,
            constraints: Vec::new(),
        }
    }

    pub fn is_ground(&self) -> bool {
        self.charge.is_ground() && self.spin.is_ground() && self.electrons.is_ground()
    }

    /// `self` (pattern) matches `target` iff every admissible assignment
    /// of `target` is also admissible by `self`, checked field-wise.
    pub fn matches(&self, target: &MulticenterBondAst) -> bool {
        self.charge.matches(&target.charge)
            && self.spin.matches(&target.spin)
            && self.electrons.matches(&target.electrons)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::default_all_undetermined(MulticenterBondAst::default(), false)]
    #[case::ground(MulticenterBondAst { charge: ValueAst::Lit(0), spin: SpinStateAst::new(0, 1),
        electrons: ValueAst::Lit(2), constraints: Vec::new() }, true)]
    fn test_multicenter_bond_ast_is_ground(
        #[case] ast: MulticenterBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_ground(), expected);
    }
}
