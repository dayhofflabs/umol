//! Multicenter bond AST.

use super::constraint::MulticenterBondConstraints;
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
    pub constraints: MulticenterBondConstraints,
}

impl MulticenterBondAst {
    pub fn new(charge: ValueAst, spin: SpinStateAst, electrons: ValueAst) -> Self {
        Self {
            charge,
            spin,
            electrons,
            constraints: MulticenterBondConstraints::new(),
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
    #[case::all_undetermined(MulticenterBondAst::default(), false)]
    #[case::charge_only(MulticenterBondAst::new(ValueAst::Lit(0), SpinStateAst::default(), ValueAst::Undetermined), false)]
    #[case::all_ground(MulticenterBondAst::new(ValueAst::Lit(0), SpinStateAst::new(0, 1), ValueAst::Lit(2)), true)]
    fn test_multicenter_bond_ast_is_ground(
        #[case] ast: MulticenterBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_ground(MulticenterBondAst::default(),
        MulticenterBondAst::new(ValueAst::Lit(0), SpinStateAst::new(0, 1), ValueAst::Lit(2)), true)]
    #[case::exact(MulticenterBondAst::new(ValueAst::Lit(0), SpinStateAst::new(0, 1), ValueAst::Lit(2)),
        MulticenterBondAst::new(ValueAst::Lit(0), SpinStateAst::new(0, 1), ValueAst::Lit(2)), true)]
    #[case::electrons_mismatch(MulticenterBondAst::new(ValueAst::Undetermined, SpinStateAst::default(), ValueAst::Lit(4)),
        MulticenterBondAst::new(ValueAst::Lit(0), SpinStateAst::new(0, 1), ValueAst::Lit(2)), false)]
    #[case::charge_mismatch(MulticenterBondAst::new(ValueAst::Lit(1), SpinStateAst::default(), ValueAst::Undetermined),
        MulticenterBondAst::new(ValueAst::Lit(0), SpinStateAst::default(), ValueAst::Undetermined), false)]
    fn test_multicenter_bond_ast_matches(
        #[case] pattern: MulticenterBondAst,
        #[case] target: MulticenterBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }
}
