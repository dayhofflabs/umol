//! Dative bond AST.

use std::mem;

use super::constraint::DativeBondConstraints;
use super::value::ValueAst;

/// Dative bond: variable-arity coordinative bond with a designated acceptor
/// atom. Donors and acceptor live together in the molecule's relation
/// participants array (sorted ascending by `NodeId`); `acceptor_slot` is the
/// index of the acceptor in that array. The two-atom case has a single donor;
/// larger participant sets describe a multi-donor / single-acceptor bond.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DativeBondAst {
    /// Index of the acceptor in the relation's sorted participants array.
    /// Always meaningful in the context of a `MoleculeAst`; meaningless on a
    /// standalone `DativeBondAst`.
    pub acceptor_slot: u8,
    /// Bond order — number of electron pairs donated.
    pub order: ValueAst,
    pub constraints: DativeBondConstraints,
}

impl DativeBondAst {
    pub fn new(order: ValueAst) -> Self {
        Self {
            acceptor_slot: 0,
            order,
            constraints: DativeBondConstraints::new(),
        }
    }

    pub fn from_order(order: u8) -> Self {
        Self::new(ValueAst::Lit(order as i64))
    }

    pub fn with_acceptor_slot(mut self, slot: u8) -> Self {
        self.acceptor_slot = slot;
        self
    }

    pub fn with_order(mut self, order: impl Into<ValueAst>) -> Self {
        self.order = order.into();
        self
    }

    pub fn with_constraints(mut self, constraints: impl Into<DativeBondConstraints>) -> Self {
        self.constraints = constraints.into();
        self
    }

    /// No-op: `DativeBondAst` has no value-bearing fields besides `order`,
    /// which is essential and never filled. Provided for API symmetry.
    pub fn into_ground(self) -> Self {
        self
    }

    /// Equivalent to `into_ground()`. `DativeBondAst` has no constraint
    /// defaults.
    pub fn into_zeroed(self) -> Self {
        self.into_ground()
    }

    pub fn is_ground(&self) -> bool {
        self.order.is_ground()
    }

    /// `self` (pattern) matches `target` iff every admissible assignment of
    /// `target` is also admissible by `self`, checked field-wise on `order`.
    /// `acceptor_slot` is structural — equality is enforced by the matching
    /// driver, not here.
    pub fn matches(&self, target: &DativeBondAst) -> bool {
        self.order.matches(&target.order)
    }

    /// Simplify every value-bearing field in place.
    pub fn simplify_values(&mut self) {
        self.order = mem::take(&mut self.order).simplify();
        self.constraints.simplify_each();
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::constraint::DativeBondConstraint;

    fn ground() -> DativeBondAst {
        DativeBondAst::from_order(1)
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(DativeBondAst::default(), false)]
    #[case::order_lit(ground(), true)]
    #[case::order_undetermined(DativeBondAst::new(ValueAst::Undetermined), false)]
    #[case::with_ground_constraint(DativeBondAst { acceptor_slot: 0, order: ValueAst::Lit(1),
        constraints: DativeBondConstraints::from_iter([DativeBondConstraint::RingSize(ValueAst::Lit(6))]) }, true)]
    fn test_dative_bond_ast_is_ground(#[case] ast: DativeBondAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rstest]
    #[case::same_order(DativeBondAst::from_order(1), DativeBondAst::from_order(1), true)]
    #[case::order_mismatch(DativeBondAst::from_order(1), DativeBondAst::from_order(2), false)]
    #[case::pattern_wildcard_matches(
        DativeBondAst::new(ValueAst::Undetermined),
        DativeBondAst::from_order(2),
        true,
    )]
    #[case::pattern_more_specific(
        DativeBondAst::from_order(2),
        DativeBondAst::new(ValueAst::Undetermined),
        false,
    )]
    fn test_dative_bond_ast_matches(
        #[case] pattern: DativeBondAst,
        #[case] target: DativeBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    fn test_dative_bond_ast_into_ground() {
        let bond = DativeBondAst::from_order(1);
        assert_eq!(bond.clone().into_ground(), bond);
    }

    #[rstest]
    fn test_dative_bond_ast_into_zeroed() {
        let bond = DativeBondAst::from_order(1);
        assert_eq!(bond.clone().into_zeroed(), bond);
    }
}
