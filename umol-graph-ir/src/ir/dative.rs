//! Dative bond form.

use umol_graph_core::{BiRelationData, ParticipantPosition};
use umol_graph_ir_macros::{Lattice, Normalize};

use super::constraint::{DativeBondConstraintForm, DativeBondConstraintsForm};
use super::num::NumForm;
use super::traits::Equiv;

/// Dative bond data: bond order (number of electron pairs donated) and
/// constraints. The acceptor and donor atoms are the participants of the
/// owning molecule's dative relation; this value holds no participant refs.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Normalize, Lattice)]
pub struct DativeBondForm {
    /// Bond order — number of electron pairs donated.
    pub order: NumForm,
    pub constraints: DativeBondConstraintsForm,
}

/// Attribute update for a dative bond. The order is optional, and
/// undetermined constraints remove their key.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DativeBondUpdate {
    pub order: Option<NumForm>,
    pub constraints: DativeBondConstraintsForm,
}

impl DativeBondForm {
    pub fn new(order: NumForm) -> Self {
        Self {
            order,
            constraints: DativeBondConstraintsForm::new(),
        }
    }

    pub fn from_order(order: u8) -> Self {
        Self::new(NumForm::Lit(order as i64))
    }

    pub fn with_order(mut self, order: impl Into<NumForm>) -> Self {
        self.order = order.into();
        self
    }

    /// Add a single constraint, replacing any existing entry of the same
    /// key (last-wins per `DativeBondConstraintsForm::set`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<DativeBondConstraintForm>) -> Self {
        self.constraints.set(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same key (last-wins per `DativeBondConstraintsForm::set`). Does
    /// not clear existing constraints; use `bond.constraints.clear()` or
    /// direct field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<DativeBondConstraintForm>,
    {
        self.constraints
            .extend(constraints.into_iter().map(Into::into));
        self
    }

    /// No-op: `DativeBondForm` has no value-bearing fields besides `order`,
    /// which is essential and never filled. Provided for API symmetry.
    pub fn into_ground(self) -> Self {
        self
    }

    /// Apply an attribute update, leaving omitted leaves and constraint keys unchanged.
    pub fn update(&self, update: &DativeBondUpdate) -> DativeBondForm {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        DativeBondForm {
            order: update.order.clone().unwrap_or_else(|| self.order.clone()),
            constraints,
        }
    }

    /// Derive the minimal canonical attribute update carrying `self` to `other`.
    pub fn difference_to(&self, other: &Self) -> DativeBondUpdate {
        let mut constraints = DativeBondConstraintsForm::new();
        for new in other.constraints.iter() {
            if self
                .constraints
                .get(new.key())
                .is_none_or(|old| !old.equiv(new))
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
            order: (!self.order.equiv(&other.order)).then(|| other.order.clone()),
            constraints,
        }
    }
}

impl From<&str> for DativeBondForm {
    fn from(s: &str) -> Self {
        s.parse().expect("invalid dative bond string")
    }
}

impl BiRelationData for DativeBondForm {
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
    use crate::ir::boolean::BooleanForm;
    use crate::ir::constraint::RingScope;
    use crate::ir::error::Contradiction;
    use crate::ir::traits::{Lattice, Normalize};

    #[rustfmt::skip]
    #[rstest]
    #[case::new(DativeBondForm::new(NumForm::Lit(2)),
        DativeBondForm { order: NumForm::Lit(2), constraints: DativeBondConstraintsForm::new() })]
    #[case::from_order(DativeBondForm::from_order(3),
        DativeBondForm { order: NumForm::Lit(3), constraints: DativeBondConstraintsForm::new() })]
    fn test_dative_bond_form_new(#[case] actual: DativeBondForm, #[case] expected: DativeBondForm) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_order(DativeBondForm::default().with_order(2),
        DativeBondForm { order: NumForm::Lit(2), constraints: DativeBondConstraintsForm::new() })]
    #[case::with_constraint(DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        DativeBondForm { order: NumForm::Lit(1),
            constraints: DativeBondConstraintsForm::from(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))) })]
    #[case::with_constraints_extends(
        DativeBondForm::from_order(1)
            .with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)))
            .with_constraints([DativeBondConstraintForm::ring_membership(RingScope::All, 1), DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1)]),
        DativeBondForm { order: NumForm::Lit(1),
            constraints: DativeBondConstraintsForm::from_iter([
                DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
                DativeBondConstraintForm::ring_membership(RingScope::All, 1),
                DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
            ]) })]
    #[case::with_constraint_appends_different_scopes(
        DativeBondForm::from_order(1)
            .with_constraint(DativeBondConstraintForm::ring_membership(RingScope::Size(5), 1))
            .with_constraint(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1)),
        DativeBondForm { order: NumForm::Lit(1),
            constraints: DativeBondConstraintsForm::from_iter([
                DativeBondConstraintForm::ring_membership(RingScope::Size(5), 1),
                DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1),
            ]) })]
    fn test_dative_bond_form_with_methods(#[case] actual: DativeBondForm, #[case] expected: DativeBondForm) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::from_order(DativeBondForm::from_order(1))]
    #[case::with_constraint(DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    fn test_dative_bond_form_into_ground(#[case] bond: DativeBondForm) {
        assert_eq!(bond.clone().into_ground(), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::order(DativeBondForm::from_order(1), DativeBondUpdate { order: Some(NumForm::Lit(2)), ..Default::default() }, DativeBondForm::from_order(2))]
    #[case::order_undetermined(DativeBondForm::from_order(1), DativeBondUpdate { order: Some(NumForm::Undetermined), ..Default::default() }, DativeBondForm::default())]
    #[case::constraint_set(DativeBondForm::from_order(1), DativeBondUpdate { constraints: DativeBondConstraintsForm::from(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))), ..Default::default() }, DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    #[case::constraint_replace(DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1_i64)), DativeBondUpdate { constraints: DativeBondConstraintsForm::from(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 2_i64)), ..Default::default() }, DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 2_i64)))]
    #[case::constraint_remove(DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1_i64)), DativeBondUpdate { constraints: DativeBondConstraintsForm::from(DativeBondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined)), ..Default::default() }, DativeBondForm::from_order(1))]
    fn test_dative_bond_form_update(
        #[case] bond: DativeBondForm,
        #[case] update: DativeBondUpdate,
        #[case] expected: DativeBondForm,
    ) {
        assert_eq!(bond.update(&update), expected);
    }

    #[rstest]
    #[case::empty(DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    fn test_dative_bond_form_update_identity(#[case] bond: DativeBondForm) {
        assert_eq!(bond.update(&DativeBondUpdate::default()), bond);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fields_and_constraints(
        DativeBondForm::from_order(1).with_constraints([
            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true)),
            DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1_i64),
        ]),
        DativeBondForm::from_order(2).with_constraints([
            DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false)),
            DativeBondConstraintForm::ring_membership(RingScope::All, 2_i64),
        ]),
        DativeBondUpdate {
            order: Some(NumForm::Lit(2)),
            constraints: DativeBondConstraintsForm::from_iter([
                DativeBondConstraintForm::Aromatic(BooleanForm::Lit(false)),
                DativeBondConstraintForm::ring_membership(RingScope::All, 2_i64),
                DativeBondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined),
            ]),
        },
    )]
    fn test_dative_bond_form_difference_to(
        #[case] bond: DativeBondForm,
        #[case] other: DativeBondForm,
        #[case] expected: DativeBondUpdate,
    ) {
        assert_eq!(bond.difference_to(&other), expected);
    }

    #[rstest]
    #[case::canonical(
        DativeBondForm::from_order(1),
        DativeBondForm::new(NumForm::lit_set([1])),
    )]
    fn test_dative_bond_form_difference_to_identity(
        #[case] bond: DativeBondForm,
        #[case] other: DativeBondForm,
    ) {
        assert_eq!(bond.difference_to(&other), DativeBondUpdate::default());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(DativeBondForm::default(), false)]
    #[case::order_lit(DativeBondForm::from_order(1), true)]
    #[case::order_undetermined(DativeBondForm::new(NumForm::Undetermined), false)]
    #[case::ground_with_constraint(DativeBondForm { order: NumForm::Lit(1),
        constraints: DativeBondConstraintsForm::from(DativeBondConstraintForm::ring_membership(RingScope::Size(6), 1)) }, true)]
    fn test_dative_bond_form_is_ground(#[case] form: DativeBondForm, #[case] expected: bool) {
        assert_eq!(form.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_order(
        DativeBondForm::new(NumForm::lit_set([1])),
        Ok(DativeBondForm::from_order(1)),
    )]
    #[case::order_empty_litset_contradiction(
        DativeBondForm::new(NumForm::lit_set(Vec::<i64>::new())),
        Err(Contradiction),
    )]
    fn test_dative_bond_form_normalize(
        #[case] input: DativeBondForm,
        #[case] expected: Result<DativeBondForm, Contradiction>,
    ) {
        assert_eq!(input.normalize(), expected);
    }

    #[rstest]
    #[case::same_order(DativeBondForm::from_order(1), DativeBondForm::from_order(1), true)]
    #[case::order_mismatch(DativeBondForm::from_order(1), DativeBondForm::from_order(2), false)]
    #[case::pattern_wildcard_matches(
        DativeBondForm::new(NumForm::Undetermined),
        DativeBondForm::from_order(2),
        true
    )]
    #[case::pattern_more_specific(
        DativeBondForm::from_order(2),
        DativeBondForm::new(NumForm::Undetermined),
        false
    )]
    #[case::constraint_required_present(
        DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        true
    )]
    #[case::constraint_required_absent(
        DativeBondForm::from_order(1).with_constraint(DativeBondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        DativeBondForm::from_order(1),
        false
    )]
    fn test_dative_bond_form_matches(
        #[case] pattern: DativeBondForm,
        #[case] target: DativeBondForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::both_default(
        DativeBondForm::default(),
        DativeBondForm::default(),
        Some(DativeBondForm::default())
    )]
    #[case::order_mismatch(DativeBondForm::from_order(1), DativeBondForm::from_order(2), None)]
    #[case::narrows_order(
        DativeBondForm { order: NumForm::Undetermined, constraints: DativeBondConstraintsForm::new() },
        DativeBondForm { order: NumForm::Lit(1), constraints: DativeBondConstraintsForm::new() },
        Some(DativeBondForm { order: NumForm::Lit(1), constraints: DativeBondConstraintsForm::new() }),
    )]
    fn test_dative_bond_form_meet(
        #[case] a: DativeBondForm,
        #[case] b: DativeBondForm,
        #[case] expected: Option<DativeBondForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::order_widens_to_lit_set(
        DativeBondForm::from_order(1),
        DativeBondForm::from_order(2),
        DativeBondForm { order: NumForm::lit_set([1, 2]), constraints: DativeBondConstraintsForm::new() },
    )]
    fn test_dative_bond_form_join(
        #[case] a: DativeBondForm,
        #[case] b: DativeBondForm,
        #[case] expected: DativeBondForm,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }
}
