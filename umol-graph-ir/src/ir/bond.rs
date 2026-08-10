//! Bond-level AST fragments shared across crates.

use umol_graph_ir_macros::{Canonicalize, Lattice};

use super::constraint::{BondConstraintForm, BondConstraintsForm};
use super::num::NumForm;
use super::spin::{UnpairedElectronsForm, UnpairedElectronsUpdate};
use super::traits::{Canonicalize, Lattice};

/// Bond AST: structural representation of a bond plus bond-level constraints
/// (aromatic flag, ring membership).
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Canonicalize, Lattice)]
pub struct BondForm {
    pub order: NumForm,
    pub charge: NumForm,
    pub unpaired_electrons: UnpairedElectronsForm,
    pub constraints: BondConstraintsForm,
}

/// Attribute update for a localized bond. Scalar fields are optional, unpaired-electron components
/// are updated independently, and undetermined constraints remove their key.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BondUpdate {
    pub order: Option<NumForm>,
    pub charge: Option<NumForm>,
    pub unpaired_electrons: UnpairedElectronsUpdate,
    pub constraints: BondConstraintsForm,
}

impl BondForm {
    pub fn new(order: NumForm) -> Self {
        Self {
            order,
            charge: NumForm::default(),
            unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: BondConstraintsForm::new(),
        }
    }

    pub fn from_order(order: u8) -> Self {
        Self::new(NumForm::Lit(order as i64))
    }

    pub fn with_order(mut self, order: impl Into<NumForm>) -> Self {
        self.order = order.into();
        self
    }

    pub fn with_charge(mut self, charge: impl Into<NumForm>) -> Self {
        self.charge = charge.into();
        self
    }

    pub fn with_unpaired_electrons(
        mut self,
        unpaired_electrons: impl Into<UnpairedElectronsForm>,
    ) -> Self {
        self.unpaired_electrons = unpaired_electrons.into();
        self
    }

    /// Add a single constraint, replacing any existing entry of the same
    /// key (last-wins per `BondConstraintsForm::set`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<BondConstraintForm>) -> Self {
        self.constraints.set(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same key (last-wins per `BondConstraintsForm::set`). Does not
    /// clear existing constraints; use `bond.constraints.clear()` or direct
    /// field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<BondConstraintForm>,
    {
        self.constraints
            .extend(constraints.into_iter().map(Into::into));
        self
    }

    /// Apply an attribute update, leaving omitted leaves and constraint keys unchanged.
    pub fn update(&self, update: &BondUpdate) -> BondForm {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        BondForm {
            order: update.order.clone().unwrap_or_else(|| self.order.clone()),
            charge: update.charge.clone().unwrap_or_else(|| self.charge.clone()),
            unpaired_electrons: self.unpaired_electrons.update(&update.unpaired_electrons),
            constraints,
        }
    }

    /// Derive the minimal canonical attribute update carrying `self` to `other`.
    pub fn difference_to(&self, other: &Self) -> BondUpdate {
        let mut constraints = BondConstraintsForm::new();
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
            unpaired_electrons: self
                .unpaired_electrons
                .difference_to(&other.unpaired_electrons),
            constraints,
        }
    }

    /// Fill `Undetermined` value-bearing struct fields with zero defaults:
    /// charge → `Lit(0)`, unpaired electrons → closed-shell singlet `(0, 1)`. Existing
    /// values and `constraints` are preserved. The result is ground iff
    /// `order` is already ground.
    pub fn into_ground(mut self) -> Self {
        if self.charge.is_undetermined() {
            self.charge = NumForm::Lit(0);
        }
        if self.unpaired_electrons.is_undetermined() {
            self.unpaired_electrons = UnpairedElectronsForm::from((0_u8, 1_u8));
        }
        self
    }
}

impl From<&str> for BondForm {
    fn from(s: &str) -> Self {
        s.parse().expect("invalid bond string")
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ir::constraint::RingScope;
    use crate::ir::error::Contradiction;
    use crate::ir::traits::{Canonicalize, Lattice};
    use crate::ir::{BooleanForm, CisTransStereoForm};

    #[rustfmt::skip]
    #[rstest]
    #[case::new(BondForm::new(NumForm::Lit(2)),
        BondForm { order: NumForm::Lit(2), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() })]
    #[case::from_order(BondForm::from_order(3),
        BondForm { order: NumForm::Lit(3), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() })]
    fn test_bond_form_new(#[case] actual: BondForm, #[case] expected: BondForm) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_order(BondForm::default().with_order(2_i64),
        BondForm { order: NumForm::Lit(2), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() })]
    #[case::with_charge(BondForm::from_order(1).with_charge(-1_i64),
        BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(-1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() })]
    #[case::with_unpaired_electrons(BondForm::from_order(1).with_unpaired_electrons((0_u8, 1_u8)),
        BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::closed_shell(), constraints: BondConstraintsForm::new() })]
    #[case::with_constraint(BondForm::from_order(1).with_constraint(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: BondConstraintsForm::from(BondConstraintForm::Aromatic(BooleanForm::Lit(true))) })]
    #[case::with_constraints_extends(
        BondForm::from_order(1)
            .with_constraint(BondConstraintForm::Aromatic(BooleanForm::Lit(true)))
            .with_constraints([BondConstraintForm::ring_membership(RingScope::All, 1), BondConstraintForm::ring_membership(RingScope::Size(6), 1)]),
        BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: BondConstraintsForm::from_iter([
                BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
                BondConstraintForm::ring_membership(RingScope::All, 1),
                BondConstraintForm::ring_membership(RingScope::Size(6), 1),
            ]) })]
    #[case::with_constraint_appends_same_scope(
        BondForm::from_order(1)
            .with_constraint(BondConstraintForm::ring_membership(RingScope::All, 1))
            .with_constraint(BondConstraintForm::ring_membership(RingScope::All, 2)),
        BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: BondConstraintsForm::from_iter([
                BondConstraintForm::ring_membership(RingScope::All, 1),
                BondConstraintForm::ring_membership(RingScope::All, 2),
            ]) })]
    #[case::with_constraint_appends_multi_valued_ring_size(
        BondForm::from_order(1)
            .with_constraint(BondConstraintForm::ring_membership(RingScope::Size(5), 1))
            .with_constraint(BondConstraintForm::ring_membership(RingScope::Size(6), 1)),
        BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: BondConstraintsForm::from_iter([
                BondConstraintForm::ring_membership(RingScope::Size(5), 1),
                BondConstraintForm::ring_membership(RingScope::Size(6), 1),
            ]) })]
    fn test_bond_form_with_methods(#[case] actual: BondForm, #[case] expected: BondForm) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::order(BondForm::from_order(1), BondUpdate { order: Some(NumForm::Lit(2)), ..Default::default() }, BondForm::from_order(2))]
    #[case::order_undetermined(BondForm::from_order(1), BondUpdate { order: Some(NumForm::Undetermined), ..Default::default() }, BondForm::default())]
    #[case::charge(BondForm::from_order(1).with_charge(0_i64), BondUpdate { charge: Some(NumForm::Lit(-1)), ..Default::default() }, BondForm::from_order(1).with_charge(-1_i64))]
    #[case::charge_undetermined(BondForm::from_order(1).with_charge(-1_i64), BondUpdate { charge: Some(NumForm::Undetermined), ..Default::default() }, BondForm::from_order(1))]
    #[case::unpaired_electrons_count(BondForm::from_order(1).with_unpaired_electrons((2_u8, 3_u8)), BondUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Lit(0)), multiplicity: None }, ..Default::default() }, BondForm::from_order(1).with_unpaired_electrons((0_u8, 3_u8)))]
    #[case::unpaired_electrons_multiplicity(BondForm::from_order(1).with_unpaired_electrons((2_u8, 3_u8)), BondUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Lit(1)) }, ..Default::default() }, BondForm::from_order(1).with_unpaired_electrons((2_u8, 1_u8)))]
    #[case::unpaired_electrons_count_undetermined(BondForm::from_order(1).with_unpaired_electrons((2_u8, 3_u8)), BondUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: Some(NumForm::Undetermined), multiplicity: None }, ..Default::default() }, BondForm::from_order(1).with_unpaired_electrons(UnpairedElectronsForm { count: NumForm::Undetermined, multiplicity: NumForm::Lit(3) }))]
    #[case::unpaired_electrons_multiplicity_undetermined(BondForm::from_order(1).with_unpaired_electrons((2_u8, 3_u8)), BondUpdate { unpaired_electrons: UnpairedElectronsUpdate { count: None, multiplicity: Some(NumForm::Undetermined) }, ..Default::default() }, BondForm::from_order(1).with_unpaired_electrons(UnpairedElectronsForm { count: NumForm::Lit(2), multiplicity: NumForm::Undetermined }))]
    #[case::constraint_set(BondForm::from_order(1), BondUpdate { constraints: BondConstraintsForm::from(BondConstraintForm::Aromatic(BooleanForm::Lit(true))), ..Default::default() }, BondForm::from_order(1).with_constraint(BondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    #[case::constraint_replace(BondForm::from_order(1).with_constraint(BondConstraintForm::ring_membership(RingScope::Size(6), 1_i64)), BondUpdate { constraints: BondConstraintsForm::from(BondConstraintForm::ring_membership(RingScope::Size(6), 2_i64)), ..Default::default() }, BondForm::from_order(1).with_constraint(BondConstraintForm::ring_membership(RingScope::Size(6), 2_i64)))]
    #[case::constraint_remove(BondForm::from_order(1).with_constraint(BondConstraintForm::ring_membership(RingScope::Size(6), 1_i64)), BondUpdate { constraints: BondConstraintsForm::from(BondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined)), ..Default::default() }, BondForm::from_order(1))]
    fn test_bond_form_update(#[case] bond: BondForm, #[case] update: BondUpdate, #[case] expected: BondForm) {
        assert_eq!(bond.update(&update), expected);
    }

    #[rstest]
    #[case::empty(BondForm::from_order(1).with_charge(-1_i64).with_unpaired_electrons((2_u8, 3_u8)).with_constraint(BondConstraintForm::Aromatic(BooleanForm::Lit(true))))]
    fn test_bond_form_update_identity(#[case] bond: BondForm) {
        assert_eq!(bond.update(&BondUpdate::default()), bond);
    }

    #[rstest]
    fn test_bond_form_difference_to() {
        let bond = BondForm::from_order(1)
            .with_charge(0_i64)
            .with_unpaired_electrons((2_u8, 3_u8))
            .with_constraints([
                BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
                BondConstraintForm::ring_membership(RingScope::Size(6), 1_i64),
            ]);
        let other = BondForm::from_order(2)
            .with_unpaired_electrons((2_u8, 1_u8))
            .with_constraints([
                BondConstraintForm::Aromatic(BooleanForm::Lit(false)),
                BondConstraintForm::CisTransStereo(CisTransStereoForm::NotStereo),
            ]);
        assert_eq!(
            bond.difference_to(&other),
            BondUpdate {
                order: Some(NumForm::Lit(2)),
                charge: Some(NumForm::Undetermined),
                unpaired_electrons: UnpairedElectronsUpdate {
                    count: None,
                    multiplicity: Some(NumForm::Lit(1)),
                },
                constraints: BondConstraintsForm::from_iter([
                    BondConstraintForm::Aromatic(BooleanForm::Lit(false)),
                    BondConstraintForm::CisTransStereo(CisTransStereoForm::NotStereo),
                    BondConstraintForm::ring_membership(RingScope::Size(6), NumForm::Undetermined,),
                ]),
            }
        );
    }

    #[rstest]
    #[case::canonical(BondForm::from_order(1).with_charge(1_i64), BondForm::from_order(1).with_charge(NumForm::lit_set([1])))]
    fn test_bond_form_difference_to_identity(#[case] bond: BondForm, #[case] other: BondForm) {
        assert_eq!(bond.difference_to(&other), BondUpdate::default());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::from_ground_order(
        BondForm::from_order(1).into_ground(),
        BondForm {
            order: NumForm::Lit(1),
            charge: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            constraints: BondConstraintsForm::new(),
        },
    )]
    #[case::preserves_set_charge(
        BondForm::from_order(2).with_charge(1_i64).into_ground(),
        BondForm {
            order: NumForm::Lit(2),
            charge: NumForm::Lit(1),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            constraints: BondConstraintsForm::new(),
        },
    )]
    #[case::preserves_constraints(
        BondForm::from_order(1).with_constraint(BondConstraintForm::Aromatic(BooleanForm::Lit(true))).into_ground(),
        BondForm {
            order: NumForm::Lit(1),
            charge: NumForm::Lit(0),
            unpaired_electrons: UnpairedElectronsForm::from((0_u8, 1_u8)),
            constraints: BondConstraintsForm::from(BondConstraintForm::Aromatic(BooleanForm::Lit(true))),
        },
    )]
    fn test_bond_form_into_ground(#[case] actual: BondForm, #[case] expected: BondForm) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(BondForm::default(), false)]
    #[case::order_only(BondForm::from_order(1), false)]
    #[case::all_ground(BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::closed_shell(),
        constraints: BondConstraintsForm::new() }, true)]
    #[case::charge_undetermined(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::closed_shell(),
        constraints: BondConstraintsForm::new() }, false)]
    #[case::ground_with_constraint(BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::closed_shell(),
        constraints: BondConstraintsForm::from(BondConstraintForm::Aromatic(BooleanForm::Lit(true))) }, true)]
    fn test_bond_form_is_ground(#[case] form: BondForm, #[case] expected: bool) {
        assert_eq!(form.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_order(
        BondForm { order: NumForm::lit_set([2]), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() },
        Ok(BondForm { order: NumForm::Lit(2), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() }),
    )]
    #[case::order_empty_litset_contradiction(
        BondForm { order: NumForm::lit_set(Vec::<i64>::new()), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() },
        Err(Contradiction),
    )]
    fn test_bond_form_canonicalize(
        #[case] input: BondForm,
        #[case] expected: Result<BondForm, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_ground(BondForm::default(),
        BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::closed_shell(), constraints: BondConstraintsForm::new() }, true)]
    #[case::same_order(BondForm::from_order(2), BondForm::from_order(2), true)]
    #[case::order_mismatch(BondForm::from_order(2), BondForm::from_order(1), false)]
    #[case::pattern_more_specific_than_target(BondForm::from_order(2), BondForm::default(), false)]
    #[case::charge_mismatch(BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(0), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() },
        BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() }, false)]
    #[case::charge_wildcard_pattern(BondForm::from_order(1),
        BondForm { order: NumForm::Lit(1), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() }, true)]
    #[case::unpaired_electrons_mismatch(BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::closed_shell(), constraints: BondConstraintsForm::new() },
        BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: (2_u8, 3_u8).into(), constraints: BondConstraintsForm::new() }, false)]
    #[case::unpaired_electrons_wildcard_pattern(BondForm::from_order(1),
        BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::closed_shell(), constraints: BondConstraintsForm::new() }, true)]
    #[case::constraint_required_present(
        BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: BondConstraintsForm::from(BondConstraintForm::Aromatic(BooleanForm::Lit(true))) },
        BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: BondConstraintsForm::from(BondConstraintForm::Aromatic(BooleanForm::Lit(true))) }, true)]
    #[case::constraint_required_absent(
        BondForm { order: NumForm::Lit(1), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
            constraints: BondConstraintsForm::from(BondConstraintForm::Aromatic(BooleanForm::Lit(true))) },
        BondForm::from_order(1), false)]
    fn test_bond_form_matches(
        #[case] pattern: BondForm,
        #[case] target: BondForm,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::both_default(BondForm::default(), BondForm::default(), Some(BondForm::default()))]
    #[case::narrows_field(
        BondForm::from_order(2),
        BondForm { order: NumForm::Undetermined, charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() },
        Some(BondForm { order: NumForm::Lit(2), charge: NumForm::Lit(1), unpaired_electrons: UnpairedElectronsForm::default(), constraints: BondConstraintsForm::new() }),
    )]
    #[case::incompatible_order(BondForm::from_order(2), BondForm::from_order(3), None)]
    fn test_bond_form_meet(
        #[case] a: BondForm,
        #[case] b: BondForm,
        #[case] expected: Option<BondForm>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::widens_to_set(BondForm::from_order(2), BondForm::from_order(3),
        BondForm { order: NumForm::lit_set([2, 3]), charge: NumForm::Undetermined, unpaired_electrons: UnpairedElectronsForm::default(),
        constraints: BondConstraintsForm::new() },
    )]
    fn test_bond_form_join(#[case] a: BondForm, #[case] b: BondForm, #[case] expected: BondForm) {
        assert_eq!(a.join(&b), Ok(expected));
    }

    #[rstest]
    #[case::changed(
        BondForm::default(),
        BondForm::from_order(2),
        true,
        BondForm::from_order(2)
    )]
    #[case::no_change(
        BondForm::from_order(2),
        BondForm::from_order(2),
        false,
        BondForm::from_order(2)
    )]
    fn test_bond_form_narrow_from(
        #[case] mut target: BondForm,
        #[case] source: BondForm,
        #[case] expected_changed: bool,
        #[case] expected_after: BondForm,
    ) {
        let changed = target.narrow_from(&source);
        assert_eq!(changed, expected_changed);
        assert_eq!(target, expected_after);
    }
}
