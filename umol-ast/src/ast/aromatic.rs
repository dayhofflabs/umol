//! Aromatic system AST.

use super::constraint::AromaticSystemConstraints;
use super::spin::SpinStateAst;
use super::value::ValueAst;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AromaticSystemAst {
    pub charge: ValueAst,
    pub spin: SpinStateAst,
    pub electrons: ValueAst,
    pub constraints: AromaticSystemConstraints,
}

impl AromaticSystemAst {
    pub fn new(charge: ValueAst, spin: SpinStateAst, electrons: ValueAst) -> Self {
        Self {
            charge,
            spin,
            electrons,
            constraints: AromaticSystemConstraints::new(),
        }
    }

    pub fn is_ground(&self) -> bool {
        self.charge.is_ground() && self.spin.is_ground() && self.electrons.is_ground()
    }

    /// `self` (pattern) matches `target` iff every admissible assignment
    /// of `target` is also admissible by `self`, checked field-wise.
    pub fn matches(&self, target: &AromaticSystemAst) -> bool {
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
    #[case::all_undetermined(AromaticSystemAst::default(), false)]
    #[case::charge_only(AromaticSystemAst::new(ValueAst::Lit(0), SpinStateAst::default(), ValueAst::Undetermined), false)]
    #[case::all_ground(AromaticSystemAst::new(ValueAst::Lit(0), SpinStateAst::new(0, 1), ValueAst::Lit(6)), true)]
    fn test_aromatic_system_ast_is_ground(
        #[case] ast: AromaticSystemAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_ground(AromaticSystemAst::default(), AromaticSystemAst::new(ValueAst::Lit(0), SpinStateAst::new(0, 1), ValueAst::Lit(6)), true)]
    #[case::exact(AromaticSystemAst::new(ValueAst::Lit(0), SpinStateAst::new(0, 1), ValueAst::Lit(6)),
        AromaticSystemAst::new(ValueAst::Lit(0), SpinStateAst::new(0, 1), ValueAst::Lit(6)), true)]
    #[case::electrons_mismatch(AromaticSystemAst::new(ValueAst::Undetermined, SpinStateAst::default(), ValueAst::Lit(6)),
        AromaticSystemAst::new(ValueAst::Lit(0), SpinStateAst::new(0, 1), ValueAst::Lit(10)), false)]
    fn test_aromatic_system_ast_matches(
        #[case] pattern: AromaticSystemAst,
        #[case] target: AromaticSystemAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }
}
