//! Aromatic system AST.

use umol_ast_macros::{Canonicalize, Lattice};
use umol_graph_core::ParticipantPosition;

use super::constraint::{AromaticSystemConstraint, AromaticSystemConstraints};
use super::electrons::ElectronCountsAst;
use super::spin::SpinStateAst;
use super::traits::Lattice;
use super::value::ValueAst;

#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Canonicalize, Lattice)]
pub struct AromaticSystemAst {
    pub electrons: ElectronCountsAst,
    pub charge: ValueAst,
    pub spin: SpinStateAst,
    pub constraints: AromaticSystemConstraints,
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
    /// kind (last-wins per `AromaticSystemConstraints::add`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<AromaticSystemConstraint>) -> Self {
        self.constraints.add(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same kind (last-wins per `AromaticSystemConstraints::add`).
    /// Does not clear existing constraints; use `system.constraints.clear()`
    /// or direct field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<AromaticSystemConstraint>,
    {
        for c in constraints {
            self.constraints.add(c.into());
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

    /// Overwrite with `other`, keeping existing values and constraints. Spin merges
    /// field-wise (unpaired / multiplicity independently).
    pub fn update(&self, other: &AromaticSystemAst) -> AromaticSystemAst {
        let mut constraints = self.constraints.clone();
        for c in other.constraints.iter() {
            constraints.remove_by_key(c.key());
        }
        for c in other.constraints.iter() {
            if !c.is_undetermined() {
                constraints.add(c.clone());
            }
        }
        AromaticSystemAst {
            electrons: if other.electrons.is_undetermined() {
                self.electrons.clone()
            } else {
                other.electrons.clone()
            },
            charge: if other.charge.is_undetermined() {
                self.charge.clone()
            } else {
                other.charge.clone()
            },
            spin: SpinStateAst {
                unpaired: if other.spin.unpaired.is_undetermined() {
                    self.spin.unpaired.clone()
                } else {
                    other.spin.unpaired.clone()
                },
                multiplicity: if other.spin.multiplicity.is_undetermined() {
                    self.spin.multiplicity.clone()
                } else {
                    other.spin.multiplicity.clone()
                },
            },
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
    #[case::new(AromaticSystemAst::new(ElectronCountsAst::Lit(vec![1; 3])),
        AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1; 3]),
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: AromaticSystemConstraints::new() })]
    #[case::from_electrons(AromaticSystemAst::from_electrons(vec![1, 1, 1]),
        AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1; 3]),
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: AromaticSystemConstraints::new() })]
    fn test_aromatic_system_ast_constructors(
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
            constraints: AromaticSystemConstraints::new() })]
    #[case::with_spin(AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_spin((0_u8, 1_u8)),
        AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]),
            charge: ValueAst::Undetermined, spin: SpinStateAst::closed_shell(),
            constraints: AromaticSystemConstraints::new() })]
    #[case::with_constraint(
        AromaticSystemAst::from_electrons(vec![1, 1, 1]).with_constraint(AromaticSystemConstraint::electron_count(6)),
        AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1, 1, 1]),
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6)) })]
    #[case::with_constraint_replaces_same_kind(
        AromaticSystemAst::default()
            .with_constraint(AromaticSystemConstraint::electron_count(2))
            .with_constraint(AromaticSystemConstraint::electron_count(6)),
        AromaticSystemAst { electrons: ElectronCountsAst::Undetermined,
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6)) })]
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
        constraints: AromaticSystemConstraints::new() })]
    #[case::preserves_set_charge(AromaticSystemAst::from_electrons(vec![1; 6]).with_charge(1_i64).into_ground(),
        AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1; 6]), charge: ValueAst::Lit(1), spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AromaticSystemConstraints::new() })]
    #[case::preserves_constraints(AromaticSystemAst::from_electrons(vec![1; 6]).with_constraint(AromaticSystemConstraint::electron_count(6)).into_ground(),
        AromaticSystemAst { electrons: ElectronCountsAst::Lit(vec![1; 6]), charge: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6)) })]
    fn test_aromatic_system_ast_into_ground(
        #[case] actual: AromaticSystemAst,
        #[case] expected: AromaticSystemAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::all_undetermined(AromaticSystemAst::default(), false)]
    #[case::charge_only(AromaticSystemAst::new(ElectronCountsAst::Undetermined).with_charge(0), false)]
    #[case::ground_no_atoms(AromaticSystemAst::new(ElectronCountsAst::Lit(Vec::new())).with_charge(0).with_spin((0, 1)), true)]
    #[case::all_ground_six(AromaticSystemAst::new(ElectronCountsAst::Lit(vec![1; 6])).with_charge(0).with_spin((0, 1)), true)]
    #[case::ground_with_constraint(AromaticSystemAst::new(ElectronCountsAst::Lit(vec![1; 6])).with_charge(0).with_spin((0, 1)).with_constraint(AromaticSystemConstraint::electron_count(6)), true)]
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
        AromaticSystemAst::new(ElectronCountsAst::Undetermined).with_constraint(AromaticSystemConstraint::electron_count(6)),
        AromaticSystemAst::new(ElectronCountsAst::Undetermined).with_constraint(AromaticSystemConstraint::electron_count(6)),
        true)]
    #[case::constraint_required_absent(
        AromaticSystemAst::new(ElectronCountsAst::Undetermined).with_constraint(AromaticSystemConstraint::electron_count(6)),
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
        assert_eq!(a.join(&b), expected);
    }

    #[rstest]
    fn test_aromatic_system_ast_permute() {
        let mut input = AromaticSystemAst::from_electrons(vec![10, 20, 30]).with_charge(-1);
        input.permute(&[
            ParticipantPosition(2),
            ParticipantPosition(0),
            ParticipantPosition(1),
        ]);
        assert_eq!(
            input,
            AromaticSystemAst::from_electrons(vec![30, 10, 20]).with_charge(-1)
        );
    }
}
