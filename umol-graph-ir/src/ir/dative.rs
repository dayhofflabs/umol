//! Dative bond AST.

use umol_graph_core::{BiRelationData, ParticipantPosition};
use umol_graph_ir_macros::{Canonicalize, Lattice};

use super::constraint::{DativeBondConstraintAst, DativeBondConstraintsAst};
use super::traits::Canonicalize;
use super::value::ValueAst;

/// Dative bond data: bond order (number of electron pairs donated) and
/// constraints. The acceptor and donor atoms are the participants of the
/// owning molecule's dative relation; this value holds no participant refs.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Canonicalize, Lattice)]
pub struct DativeBondAst {
    /// Bond order — number of electron pairs donated.
    pub order: ValueAst,
    pub constraints: DativeBondConstraintsAst,
}

/// Attribute update for a dative bond. The order is optional, and
/// undetermined constraints remove their key.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DativeBondUpdate {
    pub order: Option<ValueAst>,
    pub constraints: DativeBondConstraintsAst,
}

impl DativeBondAst {
    pub fn new(order: ValueAst) -> Self {
        Self {
            order,
            constraints: DativeBondConstraintsAst::new(),
        }
    }

    pub fn from_order(order: u8) -> Self {
        Self::new(ValueAst::Lit(order as i64))
    }

    pub fn with_order(mut self, order: impl Into<ValueAst>) -> Self {
        self.order = order.into();
        self
    }

    /// Add a single constraint, replacing any existing entry of the same
    /// key (last-wins per `DativeBondConstraintsAst::set`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<DativeBondConstraintAst>) -> Self {
        self.constraints.set(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same key (last-wins per `DativeBondConstraintsAst::set`). Does
    /// not clear existing constraints; use `bond.constraints.clear()` or
    /// direct field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<DativeBondConstraintAst>,
    {
        self.constraints
            .extend(constraints.into_iter().map(Into::into));
        self
    }

    /// No-op: `DativeBondAst` has no value-bearing fields besides `order`,
    /// which is essential and never filled. Provided for API symmetry.
    pub fn into_ground(self) -> Self {
        self
    }

    /// Apply an attribute update, leaving omitted leaves and constraint keys unchanged.
    pub fn update(&self, update: &DativeBondUpdate) -> DativeBondAst {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        DativeBondAst {
            order: update.order.clone().unwrap_or_else(|| self.order.clone()),
            constraints,
        }
    }

    /// Derive the minimal canonical attribute update carrying `self` to `other`.
    pub fn difference_to(&self, other: &Self) -> DativeBondUpdate {
        let mut constraints = DativeBondConstraintsAst::new();
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
        DativeBondUpdate {
            order: (!self.order.canonical_eq(&other.order)).then(|| other.order.clone()),
            constraints,
        }
    }
}

impl From<&str> for DativeBondAst {
    fn from(s: &str) -> Self {
        s.parse().expect("invalid dative bond string")
    }
}

impl BiRelationData for DativeBondAst {
    /// `order` is a scalar; neither the acceptor nor the donor factor is position-indexed.
    fn on_permutation(
        &mut self,
        _order_1: &[ParticipantPosition],
        _order_2: &[ParticipantPosition],
    ) {
    }

    fn is_permutation_invariant(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ir::boolean::BooleanAst;
    use crate::ir::constraint::RingScope;
    use crate::ir::error::Contradiction;
    use crate::ir::traits::{Canonicalize, Lattice};

    #[rustfmt::skip]
    #[rstest]
    #[case::new(DativeBondAst::new(ValueAst::Lit(2)),
        DativeBondAst { order: ValueAst::Lit(2), constraints: DativeBondConstraintsAst::new() })]
    #[case::from_order(DativeBondAst::from_order(3),
        DativeBondAst { order: ValueAst::Lit(3), constraints: DativeBondConstraintsAst::new() })]
    fn test_dative_bond_ast_new(#[case] actual: DativeBondAst, #[case] expected: DativeBondAst) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_order(DativeBondAst::default().with_order(2),
        DativeBondAst { order: ValueAst::Lit(2), constraints: DativeBondConstraintsAst::new() })]
    #[case::with_constraint(DativeBondAst::from_order(1).with_constraint(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        DativeBondAst { order: ValueAst::Lit(1),
            constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))) })]
    #[case::with_constraints_extends(
        DativeBondAst::from_order(1)
            .with_constraint(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)))
            .with_constraints([DativeBondConstraintAst::ring_membership(RingScope::All, 1), DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1)]),
        DativeBondAst { order: ValueAst::Lit(1),
            constraints: DativeBondConstraintsAst::from_iter([
                DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)),
                DativeBondConstraintAst::ring_membership(RingScope::All, 1),
                DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1),
            ]) })]
    #[case::with_constraint_appends_different_scopes(
        DativeBondAst::from_order(1)
            .with_constraint(DativeBondConstraintAst::ring_membership(RingScope::Size(5), 1))
            .with_constraint(DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1)),
        DativeBondAst { order: ValueAst::Lit(1),
            constraints: DativeBondConstraintsAst::from_iter([
                DativeBondConstraintAst::ring_membership(RingScope::Size(5), 1),
                DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1),
            ]) })]
    fn test_dative_bond_ast_with_methods(#[case] actual: DativeBondAst, #[case] expected: DativeBondAst) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::from_order(DativeBondAst::from_order(1))]
    #[case::with_constraint(DativeBondAst::from_order(1).with_constraint(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))))]
    fn test_dative_bond_ast_into_ground(#[case] bond: DativeBondAst) {
        assert_eq!(bond.clone().into_ground(), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::order(DativeBondAst::from_order(1), DativeBondUpdate { order: Some(ValueAst::Lit(2)), ..Default::default() }, DativeBondAst::from_order(2))]
    #[case::order_undetermined(DativeBondAst::from_order(1), DativeBondUpdate { order: Some(ValueAst::Undetermined), ..Default::default() }, DativeBondAst::default())]
    #[case::constraint_set(DativeBondAst::from_order(1), DativeBondUpdate { constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))), ..Default::default() }, DativeBondAst::from_order(1).with_constraint(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))))]
    #[case::constraint_replace(DativeBondAst::from_order(1).with_constraint(DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1_i64)), DativeBondUpdate { constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::Size(6), 2_i64)), ..Default::default() }, DativeBondAst::from_order(1).with_constraint(DativeBondConstraintAst::ring_membership(RingScope::Size(6), 2_i64)))]
    #[case::constraint_remove(DativeBondAst::from_order(1).with_constraint(DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1_i64)), DativeBondUpdate { constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined)), ..Default::default() }, DativeBondAst::from_order(1))]
    fn test_dative_bond_ast_update(
        #[case] bond: DativeBondAst,
        #[case] update: DativeBondUpdate,
        #[case] expected: DativeBondAst,
    ) {
        assert_eq!(bond.update(&update), expected);
    }

    #[rstest]
    #[case::empty(DativeBondAst::from_order(1).with_constraint(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))))]
    fn test_dative_bond_ast_update_identity(#[case] bond: DativeBondAst) {
        assert_eq!(bond.update(&DativeBondUpdate::default()), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fields_and_constraints(
        DativeBondAst::from_order(1).with_constraints([
            DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true)),
            DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1_i64),
        ]),
        DativeBondAst::from_order(2).with_constraints([
            DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false)),
            DativeBondConstraintAst::ring_membership(RingScope::All, 2_i64),
        ]),
        DativeBondUpdate {
            order: Some(ValueAst::Lit(2)),
            constraints: DativeBondConstraintsAst::from_iter([
                DativeBondConstraintAst::Aromatic(BooleanAst::Lit(false)),
                DativeBondConstraintAst::ring_membership(RingScope::All, 2_i64),
                DativeBondConstraintAst::ring_membership(RingScope::Size(6), ValueAst::Undetermined),
            ]),
        },
    )]
    fn test_dative_bond_ast_difference_to(
        #[case] bond: DativeBondAst,
        #[case] other: DativeBondAst,
        #[case] expected: DativeBondUpdate,
    ) {
        assert_eq!(bond.difference_to(&other), expected);
    }

    #[rstest]
    #[case::canonical(
        DativeBondAst::from_order(1),
        DativeBondAst::new(ValueAst::lit_set([1])),
    )]
    fn test_dative_bond_ast_difference_to_identity(
        #[case] bond: DativeBondAst,
        #[case] other: DativeBondAst,
    ) {
        assert_eq!(bond.difference_to(&other), DativeBondUpdate::default());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(DativeBondAst::default(), false)]
    #[case::order_lit(DativeBondAst::from_order(1), true)]
    #[case::order_undetermined(DativeBondAst::new(ValueAst::Undetermined), false)]
    #[case::ground_with_constraint(DativeBondAst { order: ValueAst::Lit(1),
        constraints: DativeBondConstraintsAst::from(DativeBondConstraintAst::ring_membership(RingScope::Size(6), 1)) }, true)]
    fn test_dative_bond_ast_is_ground(#[case] ast: DativeBondAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_order(
        DativeBondAst::new(ValueAst::lit_set([1])),
        Ok(DativeBondAst::from_order(1)),
    )]
    #[case::order_empty_litset_contradiction(
        DativeBondAst::new(ValueAst::lit_set(Vec::<i64>::new())),
        Err(Contradiction),
    )]
    fn test_dative_bond_ast_canonicalize(
        #[case] input: DativeBondAst,
        #[case] expected: Result<DativeBondAst, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
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
    #[case::constraint_required_present(
        DativeBondAst::from_order(1).with_constraint(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        DativeBondAst::from_order(1).with_constraint(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))),
        true
    )]
    #[case::constraint_required_absent(
        DativeBondAst::from_order(1).with_constraint(DativeBondConstraintAst::Aromatic(BooleanAst::Lit(true))),
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
    #[case::both_default(
        DativeBondAst::default(),
        DativeBondAst::default(),
        Some(DativeBondAst::default())
    )]
    #[case::order_mismatch(DativeBondAst::from_order(1), DativeBondAst::from_order(2), None)]
    #[case::narrows_order(
        DativeBondAst { order: ValueAst::Undetermined, constraints: DativeBondConstraintsAst::new() },
        DativeBondAst { order: ValueAst::Lit(1), constraints: DativeBondConstraintsAst::new() },
        Some(DativeBondAst { order: ValueAst::Lit(1), constraints: DativeBondConstraintsAst::new() }),
    )]
    fn test_dative_bond_ast_meet(
        #[case] a: DativeBondAst,
        #[case] b: DativeBondAst,
        #[case] expected: Option<DativeBondAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::order_widens_to_lit_set(
        DativeBondAst::from_order(1),
        DativeBondAst::from_order(2),
        DativeBondAst { order: ValueAst::lit_set([1, 2]), constraints: DativeBondConstraintsAst::new() },
    )]
    fn test_dative_bond_ast_join(
        #[case] a: DativeBondAst,
        #[case] b: DativeBondAst,
        #[case] expected: DativeBondAst,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }
}
