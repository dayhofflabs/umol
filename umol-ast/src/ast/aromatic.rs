//! Aromatic system AST.

use umol_ast_macros::{Canonicalize, Lattice};
use umol_graph_core::{ParticipantPosition, RelationData};

use super::constraint::{AromaticSystemConstraintAst, AromaticSystemConstraintsAst};
use super::electrons::ElectronCountsAst;
use super::spin::{SpinStateAst, SpinStateUpdate};
use super::traits::{Canonicalize, Lattice};
use super::value::ValueAst;

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Canonicalize, Lattice)]
pub struct AromaticSystemAst {
    pub electrons: ElectronCountsAst,
    pub charge: ValueAst,
    pub spin: SpinStateAst,
    pub constraints: AromaticSystemConstraintsAst,
}

/// Attribute update for an aromatic system. Ordinary fields are optional,
/// spin is updated independently by component, and undetermined constraints
/// remove their key.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AromaticSystemUpdate {
    pub electrons: Option<ElectronCountsAst>,
    pub charge: Option<ValueAst>,
    pub spin: SpinStateUpdate,
    pub constraints: AromaticSystemConstraintsAst,
}

impl From<&str> for AromaticSystemAst {
    fn from(s: &str) -> Self {
        s.parse().expect("invalid aromatic system string")
    }
}

impl RelationData for AromaticSystemAst {
    /// The per-member electron counts are positional, so they follow a participant reorder.
    fn on_permutation(&mut self, order: &[ParticipantPosition]) {
        self.electrons.permute(order);
    }

    fn is_permutation_invariant(&self) -> bool {
        self.electrons.is_undetermined()
    }
}

impl AromaticSystemAst {
    pub fn new(electrons: ElectronCountsAst) -> Self {
        Self {
            electrons,
            ..Default::default()
        }
    }

    pub fn from_electrons(electrons: Vec<i64>) -> Self {
        Self::new(ElectronCountsAst::Lit(electrons))
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
    /// kind (last-wins per `AromaticSystemConstraintsAst::set`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<AromaticSystemConstraintAst>) -> Self {
        self.constraints.set(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same kind (last-wins per `AromaticSystemConstraintsAst::set`).
    /// Does not clear existing constraints; use `system.constraints.clear()`
    /// or direct field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<AromaticSystemConstraintAst>,
    {
        for c in constraints {
            self.constraints.set(c.into());
        }
        self
    }

    /// Fill `Undetermined` value-bearing struct fields with zero defaults:
    /// charge → `Lit(0)`, spin → closed-shell singlet `(0, 1)`. `electrons`
    /// and `constraints` are preserved. The result is ground iff `electrons`
    /// is already `Lit`.
    pub fn into_ground(mut self) -> Self {
        if self.charge.is_undetermined() {
            self.charge = ValueAst::Lit(0);
        }
        if self.spin.is_undetermined() {
            self.spin = SpinStateAst::from((0_u8, 1_u8));
        }
        self
    }

    /// Apply an attribute update, leaving omitted leaves and constraint keys unchanged.
    pub fn update(&self, update: &AromaticSystemUpdate) -> AromaticSystemAst {
        let mut constraints = self.constraints.clone();
        constraints.update(&update.constraints);
        AromaticSystemAst {
            electrons: update
                .electrons
                .clone()
                .unwrap_or_else(|| self.electrons.clone()),
            charge: update.charge.clone().unwrap_or_else(|| self.charge.clone()),
            spin: self.spin.update(&update.spin),
            constraints,
        }
    }

    /// Derive the minimal canonical attribute update carrying `self` to `other`.
    pub fn difference_to(&self, other: &Self) -> AromaticSystemUpdate {
        let mut constraints = AromaticSystemConstraintsAst::new();
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
        AromaticSystemUpdate {
            electrons: (!self.electrons.canonical_eq(&other.electrons))
                .then(|| other.electrons.clone()),
            charge: (!self.charge.canonical_eq(&other.charge)).then(|| other.charge.clone()),
            spin: self.spin.difference_to(&other.spin),
            constraints,
        }
    }

    /// Reorder the positional `electrons` by `order`, tracking a participant
    /// reordering; charge / spin / constraints are positionless and unchanged.
    pub fn permute(&mut self, order: &[ParticipantPosition]) {
        self.electrons.permute(order);
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::error::Contradiction;
    use crate::ast::traits::Canonicalize;

    #[rustfmt::skip]
    #[rstest]
    #[case::literal(AromaticSystemAst::new(ElectronCountsAst::Lit(vec![1; 3])),
        AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1; 3]),
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: AromaticSystemConstraintsAst::new() })]
    fn test_aromatic_system_ast_new(
        #[case] actual: AromaticSystemAst,
        #[case] expected: AromaticSystemAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::literal(AromaticSystemAst::from_electrons(vec![1, 1, 1]),
        AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1; 3]),
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: AromaticSystemConstraintsAst::new() })]
    fn test_aromatic_system_ast_from_electrons(
        #[case] actual: AromaticSystemAst,
        #[case] expected: AromaticSystemAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_charge(AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_charge(-1),
        AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]),
            charge: ValueAst::Lit(-1), spin: SpinStateAst::default(),
            constraints: AromaticSystemConstraintsAst::new() })]
    #[case::with_spin(AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_spin((0_u8, 1_u8)),
        AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]),
            charge: ValueAst::Undetermined, spin: SpinStateAst::closed_shell(),
            constraints: AromaticSystemConstraintsAst::new() })]
    #[case::with_constraint(
        AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_constraint(AromaticSystemConstraintAst::electron_count(6)),
        AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]),
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)) })]
    #[case::with_constraint_replaces_same_kind(
        AromaticSystemAst::default()
            .with_constraint(AromaticSystemConstraintAst::electron_count(2))
            .with_constraint(AromaticSystemConstraintAst::electron_count(6)),
        AromaticSystemAst { electrons: ElectronCountsAst::Undetermined,
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)) })]
    fn test_aromatic_system_ast_with_methods(
        #[case] actual: AromaticSystemAst,
        #[case] expected: AromaticSystemAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::from_ground_electrons(AromaticSystemAst::from_electrons(vec![1; 6]).into_ground(),
        AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1; 6]), charge: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AromaticSystemConstraintsAst::new() })]
    #[case::preserves_set_charge(AromaticSystemAst::from_electrons(vec![1; 6]).with_charge(1_i64).into_ground(),
        AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1; 6]), charge: ValueAst::Lit(1), spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AromaticSystemConstraintsAst::new() })]
    #[case::preserves_constraints(AromaticSystemAst::from_electrons(vec![1; 6]).with_constraint(AromaticSystemConstraintAst::electron_count(6)).into_ground(),
        AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1; 6]), charge: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6)) })]
    fn test_aromatic_system_ast_into_ground(
        #[case] actual: AromaticSystemAst,
        #[case] expected: AromaticSystemAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::electrons(AromaticSystemAst::from_electrons(vec![1, 1, 1]), AromaticSystemUpdate { electrons: Some(ElectronCountsAst::Lit(vec![2, 2, 2])), ..Default::default() }, AromaticSystemAst::from_electrons(vec![2, 2, 2]))]
    #[case::electrons_undetermined(AromaticSystemAst::from_electrons(vec![1, 1, 1]), AromaticSystemUpdate { electrons: Some(ElectronCountsAst::Undetermined), ..Default::default() }, AromaticSystemAst::default())]
    #[case::charge(AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_charge(0_i64), AromaticSystemUpdate { charge: Some(ValueAst::Lit(-1)), ..Default::default() }, AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_charge(-1_i64))]
    #[case::charge_undetermined(AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_charge(-1_i64), AromaticSystemUpdate { charge: Some(ValueAst::Undetermined), ..Default::default() }, AromaticSystemAst::from_electrons(vec![1, 1, 1]))]
    #[case::spin_unpaired(AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_spin((2_u8, 3_u8)), AromaticSystemUpdate { spin: SpinStateUpdate { unpaired: Some(ValueAst::Lit(0)), multiplicity: None }, ..Default::default() }, AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_spin((0_u8, 3_u8)))]
    #[case::spin_multiplicity(AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_spin((2_u8, 3_u8)), AromaticSystemUpdate { spin: SpinStateUpdate { unpaired: None, multiplicity: Some(ValueAst::Lit(1)) }, ..Default::default() }, AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_spin((2_u8, 1_u8)))]
    #[case::constraint_set(AromaticSystemAst::from_electrons(vec![1, 1, 1]), AromaticSystemUpdate { constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(6_i64)), ..Default::default() }, AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_constraint(AromaticSystemConstraintAst::electron_count(6_i64)))]
    #[case::constraint_replace(AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_constraint(AromaticSystemConstraintAst::electron_count(6_i64)), AromaticSystemUpdate { constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(4_i64)), ..Default::default() }, AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_constraint(AromaticSystemConstraintAst::electron_count(4_i64)))]
    #[case::constraint_remove(AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_constraint(AromaticSystemConstraintAst::electron_count(6_i64)), AromaticSystemUpdate { constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(ValueAst::Undetermined)), ..Default::default() }, AromaticSystemAst::from_electrons(vec![1, 1, 1]))]
    fn test_aromatic_system_ast_update(
        #[case] system: AromaticSystemAst,
        #[case] update: AromaticSystemUpdate,
        #[case] expected: AromaticSystemAst,
    ) {
        assert_eq!(system.update(&update), expected);
    }

    #[rstest]
    #[case::empty(AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_charge(-1_i64).with_spin((2_u8, 3_u8)).with_constraint(AromaticSystemConstraintAst::electron_count(6_i64)))]
    fn test_aromatic_system_ast_update_identity(#[case] system: AromaticSystemAst) {
        assert_eq!(system.update(&AromaticSystemUpdate::default()), system);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::fields_and_constraints(
        AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_charge(0_i64).with_spin((2_u8, 3_u8)).with_constraint(AromaticSystemConstraintAst::electron_count(6_i64)),
        AromaticSystemAst::from_electrons(vec![2, 2, 2]).with_spin((2_u8, 1_u8)),
        AromaticSystemUpdate {
            electrons: Some(ElectronCountsAst::Lit(vec![2, 2, 2])),
            charge: Some(ValueAst::Undetermined),
            spin: SpinStateUpdate { unpaired: None, multiplicity: Some(ValueAst::Lit(1)) },
            constraints: AromaticSystemConstraintsAst::from(AromaticSystemConstraintAst::electron_count(ValueAst::Undetermined)),
        },
    )]
    fn test_aromatic_system_ast_difference_to(
        #[case] system: AromaticSystemAst,
        #[case] other: AromaticSystemAst,
        #[case] expected: AromaticSystemUpdate,
    ) {
        assert_eq!(system.difference_to(&other), expected);
    }

    #[rstest]
    #[case::canonical(
        AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_charge(1_i64),
        AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_charge(ValueAst::lit_set([1])),
    )]
    fn test_aromatic_system_ast_difference_to_identity(
        #[case] system: AromaticSystemAst,
        #[case] other: AromaticSystemAst,
    ) {
        assert_eq!(
            system.difference_to(&other),
            AromaticSystemUpdate::default()
        );
    }

    #[rstest]
    #[case::three_members(
        AromaticSystemAst::from_electrons(vec![10, 20, 30]).with_charge(-1),
        vec![
            ParticipantPosition(2),
            ParticipantPosition(0),
            ParticipantPosition(1),
        ],
        AromaticSystemAst::from_electrons(vec![30, 10, 20]).with_charge(-1),
    )]
    fn test_aromatic_system_ast_permute(
        #[case] mut input: AromaticSystemAst,
        #[case] order: Vec<ParticipantPosition>,
        #[case] expected: AromaticSystemAst,
    ) {
        input.permute(&order);
        assert_eq!(input, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::all_undetermined(AromaticSystemAst::default(), false)]
    #[case::charge_only(AromaticSystemAst::new(ElectronCountsAst::Undetermined).with_charge(0), false)]
    #[case::ground_no_atoms(AromaticSystemAst::new(ElectronCountsAst::Lit(Vec::new())).with_charge(0).with_spin((0, 1)), true)]
    #[case::all_ground_six(AromaticSystemAst::new(ElectronCountsAst::Lit(vec![1; 6])).with_charge(0).with_spin((0, 1)), true)]
    #[case::ground_with_constraint(AromaticSystemAst::new(ElectronCountsAst::Lit(vec![1; 6])).with_charge(0).with_spin((0, 1)).with_constraint(AromaticSystemConstraintAst::electron_count(6)), true)]
    fn test_aromatic_system_ast_is_ground(
        #[case] ast: AromaticSystemAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_charge(
        AromaticSystemAst::default().with_charge(ValueAst::lit_set([0])),
        Ok(AromaticSystemAst::default().with_charge(0)),
    )]
    #[case::charge_empty_litset_contradiction(
        AromaticSystemAst::default().with_charge(ValueAst::lit_set(Vec::<i64>::new())),
        Err(Contradiction),
    )]
    fn test_aromatic_system_ast_canonicalize(
        #[case] input: AromaticSystemAst,
        #[case] expected: Result<AromaticSystemAst, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_default(AromaticSystemAst::default(), AromaticSystemAst::default(), true)]
    #[case::default_matches_ground(AromaticSystemAst::default(), AromaticSystemAst::new(ElectronCountsAst::Lit(Vec::new())).with_charge(0).with_spin((0, 1)), true)]
    #[case::exact(AromaticSystemAst::new(ElectronCountsAst::Lit(vec![1; 6])).with_charge(0).with_spin((0, 1)),
        AromaticSystemAst::new(ElectronCountsAst::Lit(vec![1; 6])).with_charge(0).with_spin((0, 1)), true)]
    #[case::electrons_length_mismatch(AromaticSystemAst::new(ElectronCountsAst::Lit(vec![1; 5])),
        AromaticSystemAst::new(ElectronCountsAst::Lit(vec![1; 6])).with_charge(0).with_spin((0, 1)), false)]
    #[case::electrons_value_mismatch(AromaticSystemAst::new(ElectronCountsAst::Lit(vec![2; 6])),
        AromaticSystemAst::new(ElectronCountsAst::Lit(vec![1; 6])).with_charge(0).with_spin((0, 1)), false)]
    #[case::pattern_undetermined_electron_matches_lit(AromaticSystemAst::new(ElectronCountsAst::Undetermined),
      AromaticSystemAst::new(ElectronCountsAst::Lit(vec![1; 6])).with_charge(0).with_spin((0, 1)), true)]
    #[case::charge_mismatch(AromaticSystemAst::new(ElectronCountsAst::Undetermined).with_charge(1),
        AromaticSystemAst::new(ElectronCountsAst::Undetermined).with_charge(0), false)]
    #[case::spin_mismatch(AromaticSystemAst::new(ElectronCountsAst::Undetermined).with_spin((2_u8, 3_u8)),
        AromaticSystemAst::new(ElectronCountsAst::Undetermined).with_spin((0_u8, 1_u8)), false)]
    #[case::constraint_required_present(
        AromaticSystemAst::new(ElectronCountsAst::Undetermined).with_constraint(AromaticSystemConstraintAst::electron_count(6)),
        AromaticSystemAst::new(ElectronCountsAst::Undetermined).with_constraint(AromaticSystemConstraintAst::electron_count(6)),
        true)]
    #[case::constraint_required_absent(
        AromaticSystemAst::new(ElectronCountsAst::Undetermined).with_constraint(AromaticSystemConstraintAst::electron_count(6)),
        AromaticSystemAst::new(ElectronCountsAst::Undetermined),
        false)]
    fn test_aromatic_system_ast_matches(
        #[case] pattern: AromaticSystemAst,
        #[case] target: AromaticSystemAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    #[case::both_default(
        AromaticSystemAst::default(),
        AromaticSystemAst::default(),
        Some(AromaticSystemAst::default())
    )]
    #[case::electrons_length_mismatch(
        AromaticSystemAst::from_electrons(vec![1; 6]),
        AromaticSystemAst::from_electrons(vec![1; 5]),
        None,
    )]
    #[case::narrows_electrons(
        AromaticSystemAst::new(ElectronCountsAst::Undetermined),
        AromaticSystemAst::from_electrons(vec![1; 3]),
        Some(AromaticSystemAst::from_electrons(vec![1; 3])),
    )]
    fn test_aromatic_system_ast_meet(
        #[case] a: AromaticSystemAst,
        #[case] b: AromaticSystemAst,
        #[case] expected: Option<AromaticSystemAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::electrons_length_mismatch_widens_to_default(
        AromaticSystemAst::from_electrons(vec![1; 6]),
        AromaticSystemAst::from_electrons(vec![1; 5]),
        AromaticSystemAst::default(),
    )]
    fn test_aromatic_system_ast_join(
        #[case] a: AromaticSystemAst,
        #[case] b: AromaticSystemAst,
        #[case] expected: AromaticSystemAst,
    ) {
        assert_eq!(a.join(&b), Ok(expected));
    }
}
