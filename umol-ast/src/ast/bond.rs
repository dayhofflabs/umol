//! Bond-level AST fragments shared across crates.

use umol_ast_macros::{Canonicalize, Lattice};

use super::constraint::{BondConstraintAst, BondConstraintsAst};
use super::spin::{UnpairedElectronsAst, UnpairedElectronsUpdate};
use super::traits::{Canonicalize, Lattice};
use super::value::ValueAst;

/// Bond AST: structural representation of a bond plus bond-level constraints
/// (aromatic flag, ring membership).
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Canonicalize, Lattice)]
pub struct BondAst {
    pub order: ValueAst,
    pub charge: ValueAst,
    pub spin: UnpairedElectronsAst,
    pub constraints: BondConstraintsAst,
}

/// Attribute update for a localized bond. Scalar fields are optional, spin is
/// updated independently by component, and undetermined constraints remove
/// their key.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BondUpdate {
    pub order: Option<ValueAst>,
    pub charge: Option<ValueAst>,
    pub spin: UnpairedElectronsUpdate,
    pub constraints: BondConstraintsAst,
}

impl BondAst {
    pub fn new(order: ValueAst) -> Self {
        Self {
            order,
            charge: ValueAst::default(),
            spin: UnpairedElectronsAst::default(),
            constraints: BondConstraintsAst::new(),
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

    pub fn with_spin(mut self, spin: impl Into<UnpairedElectronsAst>) -> Self {
        self.spin = spin.into();
        self
    }

    /// Add a single constraint, replacing any existing entry of the same
    /// key (last-wins per `BondConstraintsAst::set`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<BondConstraintAst>) -> Self {
        self.constraints.set(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same key (last-wins per `BondConstraintsAst::set`). Does not
    /// clear existing constraints; use `bond.constraints.clear()` or direct
    /// field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<BondConstraintAst>,
    {
        self.constraints
            .extend(constraints.into_iter().map(Into::into));
        self
    }

    /// Apply an attribute update, leaving omitted leaves and constraint keys unchanged.
    pub fn update(&self, update: &BondUpdate) -> BondAst {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        BondAst {
            order: update.order.clone().unwrap_or_else(|| self.order.clone()),
            charge: update.charge.clone().unwrap_or_else(|| self.charge.clone()),
            spin: self.spin.update(&update.spin),
            constraints,
        }
    }

    /// Derive the minimal canonical attribute update carrying `self` to `other`.
    pub fn difference_to(&self, other: &Self) -> BondUpdate {
        let mut constraints = BondConstraintsAst::new();
        for new in other.constraints.iter() {
            if self
                .constraints
                .get(new.key())
                .is_none_or(|old| !old.canonical_eq(new))
            {
                constraints.set(new.clone());
            }
        }
        for old in self.constraints.iter() {
            if other.constraints.get(old.key()).is_none() {
                constraints.set(old.as_undetermined());
            }
        }
        BondUpdate {
            order: (!self.order.canonical_eq(&other.order)).then(|| other.order.clone()),
            charge: (!self.charge.canonical_eq(&other.charge)).then(|| other.charge.clone()),
            spin: self.spin.difference_to(&other.spin),
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
            self.spin = UnpairedElectronsAst::from((0_u8, 1_u8));
        }
        self
    }
}

impl From<&str> for BondAst {
    fn from(s: &str) -> Self {
        s.parse().expect("invalid bond string")
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
    use crate::ast::{BooleanAst, CisTransStereoAst};

    #[rustfmt::skip]
    #[rstest]
    #[case::new(BondAst::new(ValueAst::Lit(2)),
        BondAst { order: ValueAst::Lit(2), charge: ValueAst::Undetermined, spin: UnpairedElectronsAst::default(), constraints: BondConstraintsAst::new() })]
    #[case::from_order(BondAst::from_order(3),
        BondAst { order: ValueAst::Lit(3), charge: ValueAst::Undetermined, spin: UnpairedElectronsAst::default(), constraints: BondConstraintsAst::new() })]
    fn test_bond_ast_new(#[case] actual: BondAst, #[case] expected: BondAst) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_order(BondAst::default().with_order(2_i64),
        BondAst { order: ValueAst::Lit(2), charge: ValueAst::Undetermined, spin: UnpairedElectronsAst::default(), constraints: BondConstraintsAst::new() })]
    #[case::with_charge(BondAst::from_order(1).with_charge(-1_i64),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(-1), spin: UnpairedElectronsAst::default(), constraints: BondConstraintsAst::new() })]
    #[case::with_spin(BondAst::from_order(1).with_spin((0_u8, 1_u8)),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: UnpairedElectronsAst::closed_shell(), constraints: BondConstraintsAst::new() })]
    #[case::with_constraint(BondAst::from_order(1).with_constraint(BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: UnpairedElectronsAst::default(),
            constraints: BondConstraintsAst::from(BondConstraintAst::Aromatic(BooleanAst::Lit(true))) })]
    #[case::with_constraints_extends(
        BondAst::from_order(1)
            .with_constraint(BondConstraintAst::Aromatic(BooleanAst::Lit(true)))
            .with_constraints([BondConstraintAst::ring_membership(RingScope::All, 1), BondConstraintAst::ring_membership(RingScope::Size(6), 1)]),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: UnpairedElectronsAst::default(),
            constraints: BondConstraintsAst::from_iter([
                BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
                BondConstraintAst::ring_membership(RingScope::All, 1),
                BondConstraintAst::ring_membership(RingScope::Size(6), 1),
            ]) })]
    #[case::with_constraint_appends_same_scope(
        BondAst::from_order(1)
            .with_constraint(BondConstraintAst::ring_membership(RingScope::All, 1))
            .with_constraint(BondConstraintAst::ring_membership(RingScope::All, 2)),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: UnpairedElectronsAst::default(),
            constraints: BondConstraintsAst::from_iter([
                BondConstraintAst::ring_membership(RingScope::All, 1),
                BondConstraintAst::ring_membership(RingScope::All, 2),
            ]) })]
    #[case::with_constraint_appends_multi_valued_ring_size(
        BondAst::from_order(1)
            .with_constraint(BondConstraintAst::ring_membership(RingScope::Size(5), 1))
            .with_constraint(BondConstraintAst::ring_membership(RingScope::Size(6), 1)),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: UnpairedElectronsAst::default(),
            constraints: BondConstraintsAst::from_iter([
                BondConstraintAst::ring_membership(RingScope::Size(5), 1),
                BondConstraintAst::ring_membership(RingScope::Size(6), 1),
            ]) })]
    fn test_bond_ast_with_methods(#[case] actual: BondAst, #[case] expected: BondAst) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::order(BondAst::from_order(1), BondUpdate { order: Some(ValueAst::Lit(2)), ..Default::default() }, BondAst::from_order(2))]
    #[case::order_undetermined(BondAst::from_order(1), BondUpdate { order: Some(ValueAst::Undetermined), ..Default::default() }, BondAst::default())]
    #[case::charge(BondAst::from_order(1).with_charge(0_i64), BondUpdate { charge: Some(ValueAst::Lit(-1)), ..Default::default() }, BondAst::from_order(1).with_charge(-1_i64))]
    #[case::charge_undetermined(BondAst::from_order(1).with_charge(-1_i64), BondUpdate { charge: Some(ValueAst::Undetermined), ..Default::default() }, BondAst::from_order(1))]
    #[case::spin_unpaired(BondAst::from_order(1).with_spin((2_u8, 3_u8)), BondUpdate { spin: UnpairedElectronsUpdate { count: Some(ValueAst::Lit(0)), multiplicity: None }, ..Default::default() }, BondAst::from_order(1).with_spin((0_u8, 3_u8)))]
    #[case::spin_multiplicity(BondAst::from_order(1).with_spin((2_u8, 3_u8)), BondUpdate { spin: UnpairedElectronsUpdate { count: None, multiplicity: Some(ValueAst::Lit(1)) }, ..Default::default() }, BondAst::from_order(1).with_spin((2_u8, 1_u8)))]
    #[case::spin_unpaired_undetermined(BondAst::from_order(1).with_spin((2_u8, 3_u8)), BondUpdate { spin: UnpairedElectronsUpdate { count: Some(ValueAst::Undetermined), multiplicity: None }, ..Default::default() }, BondAst::from_order(1).with_spin(UnpairedElectronsAst { count: ValueAst::Undetermined, multiplicity: ValueAst::Lit(3) }))]
    #[case::spin_multiplicity_undetermined(BondAst::from_order(1).with_spin((2_u8, 3_u8)), BondUpdate { spin: UnpairedElectronsUpdate { count: None, multiplicity: Some(ValueAst::Undetermined) }, ..Default::default() }, BondAst::from_order(1).with_spin(UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined }))]
    #[case::constraint_set(BondAst::from_order(1), BondUpdate { constraints: BondConstraintsAst::from(BondConstraintAst::Aromatic(BooleanAst::Lit(true))), ..Default::default() }, BondAst::from_order(1).with_constraint(BondConstraintAst::Aromatic(BooleanAst::Lit(true))))]
    #[case::constraint_replace(BondAst::from_order(1).with_constraint(BondConstraintAst::ring_membership(RingScope::Size(6), 1_i64)), BondUpdate { constraints: BondConstraintsAst::from(BondConstraintAst::ring_membership(RingScope::Size(6), 2_i64)), ..Default::default() }, BondAst::from_order(1).with_constraint(BondConstraintAst::ring_membership(RingScope::Size(6), 2_i64)))]
    #[case::constraint_remove(BondAst::from_order(1).with_constraint(BondConstraintAst::ring_membership(RingScope::Size(6), 1_i64)), BondUpdate { constraints: BondConstraintsAst::from(BondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined)), ..Default::default() }, BondAst::from_order(1))]
    fn test_bond_ast_update(#[case] bond: BondAst, #[case] update: BondUpdate, #[case] expected: BondAst) {
        assert_eq!(bond.update(&update), expected);
    }

    #[rstest]
    #[case::empty(BondAst::from_order(1).with_charge(-1_i64).with_spin((2_u8, 3_u8)).with_constraint(BondConstraintAst::Aromatic(BooleanAst::Lit(true))))]
    fn test_bond_ast_update_identity(#[case] bond: BondAst) {
        assert_eq!(bond.update(&BondUpdate::default()), bond);
    }

    #[rstest]
    fn test_bond_ast_difference_to() {
        let bond = BondAst::from_order(1)
            .with_charge(0_i64)
            .with_spin((2_u8, 3_u8))
            .with_constraints([
                BondConstraintAst::Aromatic(BooleanAst::Lit(true)),
                BondConstraintAst::ring_membership(RingScope::Size(6), 1_i64),
            ]);
        let other = BondAst::from_order(2)
            .with_spin((2_u8, 1_u8))
            .with_constraints([
                BondConstraintAst::Aromatic(BooleanAst::Lit(false)),
                BondConstraintAst::CisTransStereo(CisTransStereoAst::NotStereo),
            ]);
        assert_eq!(
            bond.difference_to(&other),
            BondUpdate {
                order: Some(ValueAst::Lit(2)),
                charge: Some(ValueAst::Undetermined),
                spin: UnpairedElectronsUpdate {
                    count: None,
                    multiplicity: Some(ValueAst::Lit(1)),
                },
                constraints: BondConstraintsAst::from_iter([
                    BondConstraintAst::Aromatic(BooleanAst::Lit(false)),
                    BondConstraintAst::CisTransStereo(CisTransStereoAst::NotStereo),
                    BondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined,),
                ]),
            }
        );
    }

    #[rstest]
    #[case::canonical(BondAst::from_order(1).with_charge(1_i64), BondAst::from_order(1).with_charge(ValueAst::lit_set([1])))]
    fn test_bond_ast_difference_to_identity(#[case] bond: BondAst, #[case] other: BondAst) {
        assert_eq!(bond.difference_to(&other), BondUpdate::default());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::from_ground_order(
        BondAst::from_order(1).into_ground(),
        BondAst {
            order: ValueAst::Lit(1),
            charge: ValueAst::Lit(0),
            spin: UnpairedElectronsAst::from((0_u8, 1_u8)),
            constraints: BondConstraintsAst::new(),
        },
    )]
    #[case::preserves_set_charge(
        BondAst::from_order(2).with_charge(1_i64).into_ground(),
        BondAst {
            order: ValueAst::Lit(2),
            charge: ValueAst::Lit(1),
            spin: UnpairedElectronsAst::from((0_u8, 1_u8)),
            constraints: BondConstraintsAst::new(),
        },
    )]
    #[case::preserves_constraints(
        BondAst::from_order(1).with_constraint(BondConstraintAst::Aromatic(BooleanAst::Lit(true))).into_ground(),
        BondAst {
            order: ValueAst::Lit(1),
            charge: ValueAst::Lit(0),
            spin: UnpairedElectronsAst::from((0_u8, 1_u8)),
            constraints: BondConstraintsAst::from(BondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        },
    )]
    fn test_bond_ast_into_ground(#[case] actual: BondAst, #[case] expected: BondAst) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(BondAst::default(), false)]
    #[case::order_only(BondAst::from_order(1), false)]
    #[case::all_ground(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(0), spin: UnpairedElectronsAst::closed_shell(),
        constraints: BondConstraintsAst::new() }, true)]
    #[case::charge_undetermined(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: UnpairedElectronsAst::closed_shell(),
        constraints: BondConstraintsAst::new() }, false)]
    #[case::ground_with_constraint(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(0), spin: UnpairedElectronsAst::closed_shell(),
        constraints: BondConstraintsAst::from(BondConstraintAst::Aromatic(BooleanAst::Lit(true))) }, true)]
    fn test_bond_ast_is_ground(#[case] ast: BondAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_order(
        BondAst { order: ValueAst::lit_set([2]), charge: ValueAst::Undetermined, spin: UnpairedElectronsAst::default(), constraints: BondConstraintsAst::new() },
        Ok(BondAst { order: ValueAst::Lit(2), charge: ValueAst::Undetermined, spin: UnpairedElectronsAst::default(), constraints: BondConstraintsAst::new() }),
    )]
    #[case::order_empty_litset_contradiction(
        BondAst { order: ValueAst::lit_set(Vec::<i64>::new()), charge: ValueAst::Undetermined, spin: UnpairedElectronsAst::default(), constraints: BondConstraintsAst::new() },
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
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(0), spin: UnpairedElectronsAst::closed_shell(), constraints: BondConstraintsAst::new() }, true)]
    #[case::same_order(BondAst::from_order(2), BondAst::from_order(2), true)]
    #[case::order_mismatch(BondAst::from_order(2), BondAst::from_order(1), false)]
    #[case::pattern_more_specific_than_target(BondAst::from_order(2), BondAst::default(), false)]
    #[case::charge_mismatch(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(0), spin: UnpairedElectronsAst::default(), constraints: BondConstraintsAst::new() },
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(1), spin: UnpairedElectronsAst::default(), constraints: BondConstraintsAst::new() }, false)]
    #[case::charge_wildcard_pattern(BondAst::from_order(1),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Lit(1), spin: UnpairedElectronsAst::default(), constraints: BondConstraintsAst::new() }, true)]
    #[case::spin_mismatch(BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: UnpairedElectronsAst::closed_shell(), constraints: BondConstraintsAst::new() },
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: (2_u8, 3_u8).into(), constraints: BondConstraintsAst::new() }, false)]
    #[case::spin_wildcard_pattern(BondAst::from_order(1),
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: UnpairedElectronsAst::closed_shell(), constraints: BondConstraintsAst::new() }, true)]
    #[case::constraint_required_present(
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: UnpairedElectronsAst::default(),
            constraints: BondConstraintsAst::from(BondConstraintAst::Aromatic(BooleanAst::Lit(true))) },
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: UnpairedElectronsAst::default(),
            constraints: BondConstraintsAst::from(BondConstraintAst::Aromatic(BooleanAst::Lit(true))) }, true)]
    #[case::constraint_required_absent(
        BondAst { order: ValueAst::Lit(1), charge: ValueAst::Undetermined, spin: UnpairedElectronsAst::default(),
            constraints: BondConstraintsAst::from(BondConstraintAst::Aromatic(BooleanAst::Lit(true))) },
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
        BondAst { order: ValueAst::Undetermined, charge: ValueAst::Lit(1), spin: UnpairedElectronsAst::default(), constraints: BondConstraintsAst::new() },
        Some(BondAst { order: ValueAst::Lit(2), charge: ValueAst::Lit(1), spin: UnpairedElectronsAst::default(), constraints: BondConstraintsAst::new() }),
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
        BondAst { order: ValueAst::lit_set([2, 3]), charge: ValueAst::Undetermined, spin: UnpairedElectronsAst::default(),
        constraints: BondConstraintsAst::new() },
    )]
    fn test_bond_ast_join(#[case] a: BondAst, #[case] b: BondAst, #[case] expected: BondAst) {
        assert_eq!(a.join(&b), Ok(expected));
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
