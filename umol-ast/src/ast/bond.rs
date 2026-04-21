//! Bond-level AST fragments shared across crates.

use super::constraint::BondConstraint;
use super::spin::SpinStateAst;
use super::value::ValueAst;

/// Bond AST: structural representation of a bond plus bond-level constraints
/// (aromatic flag, ring membership).
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct BondAst {
    pub order: ValueAst,
    pub charge: ValueAst,
    pub spin: SpinStateAst,
    pub constraints: Vec<BondConstraint>,
}

impl BondAst {
    pub fn new(order: ValueAst) -> Self {
        Self {
            order,
            charge: ValueAst::default(),
            spin: SpinStateAst::default(),
            constraints: Vec::new(),
        }
    }

    pub fn from_order(order: u8) -> Self {
        Self::new(ValueAst::Lit(order as i64))
    }

    pub fn is_ground(&self) -> bool {
        self.order.is_ground() && self.charge.is_ground() && self.spin.is_ground()
    }

    /// `self` (pattern) matches `target` iff every admissible assignment
    /// of `target` is also admissible by `self`, checked field-wise.
    pub fn matches(&self, target: &BondAst) -> bool {
        self.order.matches(&target.order)
            && self.charge.matches(&target.charge)
            && self.spin.matches(&target.spin)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    fn ground() -> BondAst {
        BondAst {
            order: ValueAst::Lit(1),
            charge: ValueAst::Lit(0),
            spin: SpinStateAst::closed_shell(),
            constraints: Vec::new(),
        }
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::new(BondAst::new(ValueAst::Lit(2)),
        BondAst { order: ValueAst::Lit(2), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: Vec::new() })]
    #[case::from_order(BondAst::from_order(3),
        BondAst { order: ValueAst::Lit(3), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: Vec::new() })]
    fn test_bond_ast_new(#[case] actual: BondAst, #[case] expected: BondAst) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(BondAst::default(), false)]
    #[case::order_only(BondAst::from_order(1), false)]
    #[case::all_ground(ground(), true)]
    #[case::charge_undetermined(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::closed_shell(),
        constraints: Vec::new() }, false)]
    fn test_bond_ast_is_ground(#[case] ast: BondAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_ground(BondAst::default(), ground(), true)]
    #[case::same_order(BondAst::from_order(2), BondAst::from_order(2), true)]
    #[case::order_mismatch(BondAst::from_order(2), BondAst::from_order(1), false)]
    #[case::pattern_more_specific_than_target(BondAst::from_order(2), BondAst::default(), false)]
    #[case::charge_mismatch(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(0), spin: SpinStateAst::default(), constraints: Vec::new() },
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: Vec::new() }, false)]
    #[case::charge_wildcard_pattern(BondAst::from_order(1),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: Vec::new() }, true)]
    fn test_bond_ast_matches(
        #[case] pattern: BondAst,
        #[case] target: BondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }
}
