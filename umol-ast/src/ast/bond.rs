//! Bond-level AST fragments shared across crates.

use umol_ast_macros::{Canonicalize, Lattice};

use super::constraint::{BondConstraint, BondConstraints};
use super::spin::SpinStateAst;
use super::traits::Lattice;
use super::value::ValueAst;

/// Bond AST: structural representation of a bond plus bond-level constraints
/// (aromatic flag, ring membership).
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Canonicalize, Lattice)]
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

    /// Add a single constraint, replacing any existing entry of the same
    /// kind (last-wins per `BondConstraints::add`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<BondConstraint>) -> Self {
        self.constraints.add(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same kind (last-wins per `BondConstraints::add`). Does not
    /// clear existing constraints; use `bond.constraints.clear()` or direct
    /// field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<BondConstraint>,
    {
        for c in constraints {
            self.constraints.add(c.into());
        }
        self
    }

    /// Overwrite with `other`, keeping existing values and constraints.
    pub fn update(&self, other: &BondAst) -> BondAst {
        let mut constraints = self.constraints.clone();
        for c in other.constraints.iter() {
            constraints.remove_by_key(c.key());
        }
        for c in other.constraints.iter() {
            if !c.is_undetermined() {
                constraints.add(c.clone());
            }
        }
        BondAst {
            order: if other.order.is_undetermined() {
                self.order.clone()
            } else {
                other.order.clone()
            },
            charge: if other.charge.is_undetermined() {
                self.charge.clone()
            } else {
                other.charge.clone()
            },
            spin: if other.spin.is_undetermined() {
                self.spin.clone()
            } else {
                other.spin.clone()
            },
            constraints,
        }
    }

    /// Fill `Undetermined` value-bearing struct fields with zero defaults:
    /// charge → `Lit(0)`, spin → closed-shell singlet `(0, 1)`. Existing
    /// values and `constraints` are preserved. The result is ground iff
    /// `order` is already ground.
    pub fn into_ground(mut self) -> Self {
        if self.charge.is_undetermined() {
            self.charge = ValueAst::Lit(0);
        }
        if self.spin.is_undetermined() {
            self.spin = SpinStateAst::from((0_u8, 1_u8));
        }
        self
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::constraint::RingScope;
    use crate::ast::error::Contradiction;
    use crate::ast::traits::{Canonicalize, Lattice};
    use crate::ast::BooleanAst;

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
    #[case::with_order(BondAst::default().with_order(2_i64),
        BondAst { order: ValueAst::Lit(2), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::new() })]
    #[case::with_charge(BondAst::from_order(1).with_charge(-1_i64),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(-1), spin: SpinStateAst::default(), constraints: BondConstraints::new() })]
    #[case::with_spin(BondAst::from_order(1).with_spin((0_u8, 1_u8)),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::closed_shell(), constraints: BondConstraints::new() })]
    #[case::with_constraint(BondAst::from_order(1).with_constraint(BondConstraint::Aromatic(BooleanAst::Lit(true))),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: BondConstraints::from(BondConstraint::Aromatic(BooleanAst::Lit(true))) })]
    #[case::with_constraints_extends(
        BondAst::from_order(1)
            .with_constraint(BondConstraint::Aromatic(BooleanAst::Lit(true)))
            .with_constraints([BondConstraint::ring_membership(RingScope::All, 1), BondConstraint::ring_membership(RingScope::Size(6), 1)]),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: BondConstraints::from_iter([
                BondConstraint::Aromatic(BooleanAst::Lit(true)),
                BondConstraint::ring_membership(RingScope::All, 1),
                BondConstraint::ring_membership(RingScope::Size(6), 1),
            ]) })]
    #[case::with_constraint_appends_same_scope(
        BondAst::from_order(1)
            .with_constraint(BondConstraint::ring_membership(RingScope::All, 1))
            .with_constraint(BondConstraint::ring_membership(RingScope::All, 2)),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: BondConstraints::from_iter([
                BondConstraint::ring_membership(RingScope::All, 1),
                BondConstraint::ring_membership(RingScope::All, 2),
            ]) })]
    #[case::with_constraint_appends_multi_valued_ring_size(
        BondAst::from_order(1)
            .with_constraint(BondConstraint::ring_membership(RingScope::Size(5), 1))
            .with_constraint(BondConstraint::ring_membership(RingScope::Size(6), 1)),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: BondConstraints::from_iter([
                BondConstraint::ring_membership(RingScope::Size(5), 1),
                BondConstraint::ring_membership(RingScope::Size(6), 1),
            ]) })]
    fn test_bond_ast_with_methods(#[case] actual: BondAst, #[case] expected: BondAst) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::overwrite(
        BondAst::from_order(1).with_charge(0_i64).with_constraint(BondConstraint::ring_membership(RingScope::Size(6), 1_i64)),
        BondAst::new(ValueAst::Undetermined).with_charge(-1_i64).with_constraint(BondConstraint::ring_membership(RingScope::Size(6), 2_i64)),
        BondAst::from_order(1).with_charge(-1_i64).with_constraint(BondConstraint::ring_membership(RingScope::Size(6), 2_i64)))]
    #[case::remove_constraint(
        BondAst::from_order(1).with_constraint(BondConstraint::ring_membership(RingScope::Size(6), 1_i64)),
        BondAst::new(ValueAst::Undetermined).with_constraint(BondConstraint::ring_membership(RingScope::Size(6), ValueAst::Undetermined)),
        BondAst::from_order(1))]
    fn test_bond_ast_update(#[case] lhs: BondAst, #[case] rhs: BondAst, #[case] expected: BondAst) {
        assert_eq!(lhs.update(&rhs), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::from_ground_order(
        BondAst::from_order(1).into_ground(),
        BondAst {
            order: ValueAst::Lit(1),
            charge: ValueAst::Lit(0),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            constraints: BondConstraints::new(),
        },
    )]
    #[case::preserves_set_charge(
        BondAst::from_order(2).with_charge(1_i64).into_ground(),
        BondAst {
            order: ValueAst::Lit(2),
            charge: ValueAst::Lit(1),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            constraints: BondConstraints::new(),
        },
    )]
    #[case::preserves_constraints(
        BondAst::from_order(1).with_constraint(BondConstraint::Aromatic(BooleanAst::Lit(true))).into_ground(),
        BondAst {
            order: ValueAst::Lit(1),
            charge: ValueAst::Lit(0),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            constraints: BondConstraints::from(BondConstraint::Aromatic(BooleanAst::Lit(true))),
        },
    )]
    fn test_bond_ast_into_ground(#[case] actual: BondAst, #[case] expected: BondAst) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(BondAst::default(), false)]
    #[case::order_only(BondAst::from_order(1), false)]
    #[case::all_ground(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(0), spin: SpinStateAst::closed_shell(),
        constraints: BondConstraints::new() }, true)]
    #[case::charge_undetermined(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::closed_shell(),
        constraints: BondConstraints::new() }, false)]
    #[case::ground_with_constraint(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(0), spin: SpinStateAst::closed_shell(),
        constraints: BondConstraints::from(BondConstraint::Aromatic(BooleanAst::Lit(true))) }, true)]
    fn test_bond_ast_is_ground(#[case] ast: BondAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_order(
        BondAst { order: ValueAst::lit_set([2]), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::new() },
        Ok(BondAst { order: ValueAst::Lit(2), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::new() }),
    )]
    #[case::order_empty_litset_contradiction(
        BondAst { order: ValueAst::lit_set(Vec::<i64>::new()), charge: ValueAst::Undetermined, spin: SpinStateAst::default(), constraints: BondConstraints::new() },
        Err(Contradiction),
    )]
    fn test_bond_ast_canonicalize(
        #[case] input: BondAst,
        #[case] expected: Result<BondAst, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_ground(BondAst::default(),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(0), spin: SpinStateAst::closed_shell(), constraints: BondConstraints::new() }, true)]
    #[case::same_order(BondAst::from_order(2), BondAst::from_order(2), true)]
    #[case::order_mismatch(BondAst::from_order(2), BondAst::from_order(1), false)]
    #[case::pattern_more_specific_than_target(BondAst::from_order(2), BondAst::default(), false)]
    #[case::charge_mismatch(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(0), spin: SpinStateAst::default(), constraints: BondConstraints::new() },
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: BondConstraints::new() }, false)]
    #[case::charge_wildcard_pattern(BondAst::from_order(1),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: BondConstraints::new() }, true)]
    #[case::spin_mismatch(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::closed_shell(), constraints: BondConstraints::new() },
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: (2_u8, 3_u8).into(), constraints: BondConstraints::new() }, false)]
    #[case::spin_wildcard_pattern(BondAst::from_order(1),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::closed_shell(), constraints: BondConstraints::new() }, true)]
    #[case::constraint_required_present(
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: BondConstraints::from(BondConstraint::Aromatic(BooleanAst::Lit(true))) },
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: BondConstraints::from(BondConstraint::Aromatic(BooleanAst::Lit(true))) }, true)]
    #[case::constraint_required_absent(
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: BondConstraints::from(BondConstraint::Aromatic(BooleanAst::Lit(true))) },
        BondAst::from_order(1), false)]
    fn test_bond_ast_matches(
        #[case] pattern: BondAst,
        #[case] target: BondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::both_default(BondAst::default(), BondAst::default(), Some(BondAst::default()))]
    #[case::narrows_field(
        BondAst::from_order(2),
        BondAst { order: ValueAst::Undetermined, charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: BondConstraints::new() },
        Some(BondAst { order: ValueAst::Lit(2), charge: ValueAst::Lit(1), spin: SpinStateAst::default(), constraints: BondConstraints::new() }),
    )]
    #[case::incompatible_order(BondAst::from_order(2), BondAst::from_order(3), None)]
    fn test_bond_ast_meet(
        #[case] a: BondAst,
        #[case] b: BondAst,
        #[case] expected: Option<BondAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::widens_to_set(BondAst::from_order(2), BondAst::from_order(3),
        BondAst { order: ValueAst::lit_set([2, 3]), charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
        constraints: BondConstraints::new() },
    )]
    fn test_bond_ast_join(#[case] a: BondAst, #[case] b: BondAst, #[case] expected: BondAst) {
        assert_eq!(a.join(&b), expected);
    }

    #[rstest]
    #[case::changed(
        BondAst::default(),
        BondAst::from_order(2),
        true,
        BondAst::from_order(2)
    )]
    #[case::no_change(
        BondAst::from_order(2),
        BondAst::from_order(2),
        false,
        BondAst::from_order(2)
    )]
    fn test_bond_ast_narrow_from(
        #[case] mut target: BondAst,
        #[case] source: BondAst,
        #[case] expected_changed: bool,
        #[case] expected_after: BondAst,
    ) {
        let changed = target.narrow_from(&source);
        assert_eq!(changed, expected_changed);
        assert_eq!(target, expected_after);
    }
}
