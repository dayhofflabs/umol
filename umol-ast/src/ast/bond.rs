//! Bond-level AST fragments shared across crates.

use std::mem;

use super::constraint::BondConstraints;
use super::spin::SpinStateAst;
use super::value::ValueAst;

/// Bond AST: structural representation of a bond plus bond-level constraints
/// (aromatic flag, ring membership).
#[derive(Clone, Debug, Default, PartialEq, Eq, Hash)]
pub struct BondAst {
    pub order: ValueAst,
    pub charge: ValueAst,
    pub spin: SpinStateAst,
    pub constraints: BondConstraints,
}

impl BondAst {
    pub fn new(order: ValueAst) -> Self {
        Self {
            order,
            charge: ValueAst::default(),
            spin: SpinStateAst::default(),
            constraints: BondConstraints::new(),
        }
    }

    pub fn from_order(order: u8) -> Self {
        Self::new(ValueAst::Lit(order as i64))
    }

    pub fn with_order(mut self, order: impl Into<ValueAst>) -> Self {
        self.order = order.into();
        self
    }

    pub fn with_charge(mut self, charge: impl Into<ValueAst>) -> Self {
        self.charge = charge.into();
        self
    }

    pub fn with_spin(mut self, spin: impl Into<SpinStateAst>) -> Self {
        self.spin = spin.into();
        self
    }

    pub fn with_constraints(mut self, constraints: impl Into<BondConstraints>) -> Self {
        self.constraints = constraints.into();
        self
    }

    /// Fill `Undetermined` value-bearing fields with zero defaults: charge
    /// to `Lit(0)`, spin to closed-shell singlet `(0, 1)`. Existing values
    /// and `constraints` are preserved. If `order` is `Undetermined`, the
    /// result is not ground.
    pub fn zeroed(mut self) -> Self {
        if self.charge.is_undetermined() {
            self.charge = ValueAst::Lit(0);
        }
        if self.spin.is_undetermined() {
            self.spin = SpinStateAst::from((0_u8, 1_u8));
        }
        self
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

    /// Simplify every value-bearing field in place.
    pub fn simplify_values(&mut self) {
        self.order = mem::take(&mut self.order).simplify();
        self.charge = mem::take(&mut self.charge).simplify();
        self.spin.simplify_values();
        self.constraints.simplify_each();
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
            constraints: BondConstraints::new(),
        }
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::new(BondAst::new(ValueAst::Lit(2)),
        BondAst { order: ValueAst::Lit(2), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::new() })]
    #[case::from_order(BondAst::from_order(3),
        BondAst { order: ValueAst::Lit(3), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::new() })]
    fn test_bond_ast_new(#[case] actual: BondAst, #[case] expected: BondAst) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(BondAst::default(), false)]
    #[case::order_only(BondAst::from_order(1), false)]
    #[case::all_ground(ground(), true)]
    #[case::charge_undetermined(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::closed_shell(),
        constraints: BondConstraints::new() }, false)]
    fn test_bond_ast_is_ground(#[case] ast: BondAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_ground(BondAst::default(), ground(), true)]
    #[case::same_order(BondAst::from_order(2), BondAst::from_order(2), true)]
    #[case::order_mismatch(BondAst::from_order(2), BondAst::from_order(1), false)]
    #[case::pattern_more_specific_than_target(BondAst::from_order(2), BondAst::default(), false)]
    #[case::charge_mismatch(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(0), spin: SpinStateAst::default(), constraints: BondConstraints::new() },
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: BondConstraints::new() }, false)]
    #[case::charge_wildcard_pattern(BondAst::from_order(1),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: BondConstraints::new() }, true)]
    fn test_bond_ast_matches(
        #[case] pattern: BondAst,
        #[case] target: BondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::from_ground_order(
        BondAst::from_order(1).zeroed(),
        BondAst {
            order: ValueAst::Lit(1),
            charge: ValueAst::Lit(0),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            constraints: BondConstraints::new(),
        },
    )]
    #[case::preserves_set_charge(
        BondAst::from_order(2).with_charge(1_i64).zeroed(),
        BondAst {
            order: ValueAst::Lit(2),
            charge: ValueAst::Lit(1),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            constraints: BondConstraints::new(),
        },
    )]
    fn test_bond_ast_zeroed(#[case] actual: BondAst, #[case] expected: BondAst) {
        assert_eq!(actual, expected);
    }
}
