//! Dative bond AST.

use std::mem;

use super::constraint::{DativeBondConstraint, DativeBondConstraints};
use super::traits::Lattice;
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

    /// Add a single constraint, replacing any existing entry of the same
    /// kind (last-wins per `DativeBondConstraints::add`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<DativeBondConstraint>) -> Self {
        self.constraints.add(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same kind (last-wins per `DativeBondConstraints::add`). Does
    /// not clear existing constraints; use `bond.constraints.clear()` or
    /// direct field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<DativeBondConstraint>,
    {
        for c in constraints {
            self.constraints.add(c.into());
        }
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

    /// Simplify every value-bearing field in place.
    pub fn simplify_values(&mut self) {
        self.order = mem::take(&mut self.order).simplify();
        self.constraints.simplify_each();
    }
}

impl Lattice for DativeBondAst {
    fn is_undetermined(&self) -> bool {
        self.order.is_undetermined() && self.constraints.is_undetermined()
    }

    fn is_ground(&self) -> bool {
        self.order.is_ground() && self.constraints.is_ground()
    }

    /// `acceptor_slot` is a structural anchor (u8 equality required).
    fn meet(&self, other: &Self) -> Option<Self> {
        if self.acceptor_slot != other.acceptor_slot {
            return None;
        }
        Some(Self {
            acceptor_slot: self.acceptor_slot,
            order: self.order.meet(&other.order)?,
            constraints: self.constraints.meet(&other.constraints)?,
        })
    }

    /// Mismatched `acceptor_slot` widens to `Self::default()` (entity-level
    /// undetermined).
    fn join(&self, other: &Self) -> Self {
        if self.acceptor_slot != other.acceptor_slot {
            return Self::default();
        }
        Self {
            acceptor_slot: self.acceptor_slot,
            order: self.order.join(&other.order),
            constraints: self.constraints.join(&other.constraints),
        }
    }

    /// `acceptor_slot` is a structural anchor: equality required.
    fn matches(&self, target: &Self) -> bool {
        self.acceptor_slot == target.acceptor_slot
            && self.order.matches(&target.order)
            && self.constraints.matches(&target.constraints)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::value::Expr;

    #[rustfmt::skip]
    #[rstest]
    #[case::new(DativeBondAst::new(ValueAst::Lit(2)),
        DativeBondAst { acceptor_slot: 0, order: ValueAst::Lit(2), constraints: DativeBondConstraints::new() })]
    #[case::from_order(DativeBondAst::from_order(3),
        DativeBondAst { acceptor_slot: 0, order: ValueAst::Lit(3), constraints: DativeBondConstraints::new() })]
    fn test_dative_bond_ast_new(#[case] actual: DativeBondAst, #[case] expected: DativeBondAst) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_acceptor_slot(DativeBondAst::from_order(1).with_acceptor_slot(2),
        DativeBondAst { acceptor_slot: 2, order: ValueAst::Lit(1), constraints: DativeBondConstraints::new() })]
    #[case::with_order(DativeBondAst::default().with_order(2),
        DativeBondAst { acceptor_slot: 0, order: ValueAst::Lit(2), constraints: DativeBondConstraints::new() })]
    #[case::with_constraint(DativeBondAst::from_order(1).with_constraint(DativeBondConstraint::Aromatic),
        DativeBondAst { acceptor_slot: 0, order: ValueAst::Lit(1),
            constraints: DativeBondConstraints::from(DativeBondConstraint::Aromatic) })]
    #[case::with_constraints_extends(
        DativeBondAst::from_order(1)
            .with_constraint(DativeBondConstraint::Aromatic)
            .with_constraints([DativeBondConstraint::ring_count(1), DativeBondConstraint::ring_size(6)]),
        DativeBondAst { acceptor_slot: 0, order: ValueAst::Lit(1),
            constraints: DativeBondConstraints::from_iter([
                DativeBondConstraint::Aromatic,
                DativeBondConstraint::ring_count(1),
                DativeBondConstraint::ring_size(6),
            ]) })]
    #[case::with_constraint_replaces_same_kind(
        DativeBondAst::from_order(1)
            .with_constraint(DativeBondConstraint::ring_size(5))
            .with_constraint(DativeBondConstraint::ring_size(6)),
        DativeBondAst { acceptor_slot: 0, order: ValueAst::Lit(1),
            constraints: DativeBondConstraints::from(DativeBondConstraint::ring_size(6)) })]
    fn test_dative_bond_ast_with_methods(#[case] actual: DativeBondAst, #[case] expected: DativeBondAst) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::from_order(DativeBondAst::from_order(1))]
    #[case::with_constraint(DativeBondAst::from_order(1).with_constraint(DativeBondConstraint::Aromatic))]
    fn test_dative_bond_ast_into_ground(#[case] bond: DativeBondAst) {
        assert_eq!(bond.clone().into_ground(), bond);
    }

    #[rstest]
    #[case::from_order(DativeBondAst::from_order(1))]
    #[case::with_constraint(DativeBondAst::from_order(1).with_constraint(DativeBondConstraint::Aromatic))]
    fn test_dative_bond_ast_into_zeroed(#[case] bond: DativeBondAst) {
        assert_eq!(bond.clone().into_zeroed(), bond.into_ground());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(DativeBondAst::default(), false)]
    #[case::order_lit(DativeBondAst::from_order(1), true)]
    #[case::order_undetermined(DativeBondAst::new(ValueAst::Undetermined), false)]
    #[case::ground_with_constraint(DativeBondAst { acceptor_slot: 0, order: ValueAst::Lit(1),
        constraints: DativeBondConstraints::from(DativeBondConstraint::ring_size(6)) }, true)]
    fn test_dative_bond_ast_is_ground(#[case] ast: DativeBondAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rstest]
    #[case::same_order(DativeBondAst::from_order(1), DativeBondAst::from_order(1), true)]
    #[case::order_mismatch(DativeBondAst::from_order(1), DativeBondAst::from_order(2), false)]
    #[case::pattern_wildcard_matches(
        DativeBondAst::new(ValueAst::Undetermined),
        DativeBondAst::from_order(2),
        true
    )]
    #[case::pattern_more_specific(
        DativeBondAst::from_order(2),
        DativeBondAst::new(ValueAst::Undetermined),
        false
    )]
    #[case::acceptor_slot_mismatch(
        DativeBondAst::from_order(1).with_acceptor_slot(0),
        DativeBondAst::from_order(1).with_acceptor_slot(1),
        false
    )]
    #[case::acceptor_slot_match_with_wildcard_order(
        DativeBondAst::new(ValueAst::Undetermined).with_acceptor_slot(2),
        DativeBondAst::from_order(1).with_acceptor_slot(2),
        true
    )]
    #[case::constraint_required_present(
        DativeBondAst::from_order(1).with_constraint(DativeBondConstraint::Aromatic),
        DativeBondAst::from_order(1).with_constraint(DativeBondConstraint::Aromatic),
        true
    )]
    #[case::constraint_required_absent(
        DativeBondAst::from_order(1).with_constraint(DativeBondConstraint::Aromatic),
        DativeBondAst::from_order(1),
        false
    )]
    fn test_dative_bond_ast_matches(
        #[case] pattern: DativeBondAst,
        #[case] target: DativeBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    fn test_dative_bond_ast_simplify_values() {
        let mut bond = DativeBondAst {
            acceptor_slot: 0,
            order: ValueAst::Expr(Box::new(Expr::Lit(2))),
            constraints: DativeBondConstraints::from(DativeBondConstraint::RingSize(
                ValueAst::Expr(Box::new(Expr::Lit(6))),
            )),
        };
        bond.simplify_values();
        assert_eq!(
            bond,
            DativeBondAst {
                acceptor_slot: 0,
                order: ValueAst::Lit(2),
                constraints: DativeBondConstraints::from(DativeBondConstraint::ring_size(6)),
            }
        );
    }

    #[rstest]
    #[case::both_default(
        DativeBondAst::default(),
        DativeBondAst::default(),
        Some(DativeBondAst::default())
    )]
    #[case::acceptor_slot_mismatch(
        DativeBondAst { acceptor_slot: 0, order: ValueAst::Lit(1), constraints: DativeBondConstraints::new() },
        DativeBondAst { acceptor_slot: 1, order: ValueAst::Lit(1), constraints: DativeBondConstraints::new() },
        None,
    )]
    #[case::narrows_order(
        DativeBondAst { acceptor_slot: 0, order: ValueAst::Undetermined, constraints: DativeBondConstraints::new() },
        DativeBondAst { acceptor_slot: 0, order: ValueAst::Lit(1), constraints: DativeBondConstraints::new() },
        Some(DativeBondAst { acceptor_slot: 0, order: ValueAst::Lit(1), constraints: DativeBondConstraints::new() }),
    )]
    fn test_dative_bond_ast_meet(
        #[case] a: DativeBondAst,
        #[case] b: DativeBondAst,
        #[case] expected: Option<DativeBondAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::acceptor_slot_mismatch_widens_to_default(
        DativeBondAst { acceptor_slot: 0, order: ValueAst::Lit(1), constraints: DativeBondConstraints::new() },
        DativeBondAst { acceptor_slot: 1, order: ValueAst::Lit(1), constraints: DativeBondConstraints::new() },
        DativeBondAst::default(),
    )]
    fn test_dative_bond_ast_join(
        #[case] a: DativeBondAst,
        #[case] b: DativeBondAst,
        #[case] expected: DativeBondAst,
    ) {
        assert_eq!(a.join(&b), expected);
    }
}
