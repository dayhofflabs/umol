//! Multicenter bond AST.

use std::mem;

use super::constraint::{MulticenterBondConstraint, MulticenterBondConstraints};
use super::electrons::ElectronCountsAst;
use super::spin::SpinStateAst;
use super::traits::Lattice;
use super::value::ValueAst;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MulticenterBondAst {
    pub electrons: ElectronCountsAst,
    pub charge: ValueAst,
    pub spin: SpinStateAst,
    pub constraints: MulticenterBondConstraints,
}

impl MulticenterBondAst {
    pub fn new(electrons: ElectronCountsAst) -> Self {
        Self {
            electrons,
            ..Default::default()
        }
    }

    pub fn from_counts(electrons: Vec<i64>) -> Self {
        Self::new(ElectronCountsAst::Lit(electrons))
    }

    pub fn with_electrons(mut self, electrons: impl Into<ElectronCountsAst>) -> Self {
        self.electrons = electrons.into();
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
    /// kind (last-wins per `MulticenterBondConstraints::add`). Chainable.
    pub fn with_constraint(mut self, constraint: impl Into<MulticenterBondConstraint>) -> Self {
        self.constraints.add(constraint.into());
        self
    }

    /// Add each constraint from the iterator, replacing any existing entry
    /// of the same kind (last-wins per `MulticenterBondConstraints::add`).
    /// Does not clear existing constraints; use `bond.constraints.clear()`
    /// or direct field assignment for wipe-and-replace.
    pub fn with_constraints<I>(mut self, constraints: I) -> Self
    where
        I: IntoIterator,
        I::Item: Into<MulticenterBondConstraint>,
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

    /// Equivalent to `into_ground()`. `MulticenterBondAst` has no constraint
    /// defaults.
    pub fn into_zeroed(self) -> Self {
        self.into_ground()
    }

    pub fn simplify_values(&mut self) {
        self.charge = mem::take(&mut self.charge).simplify();
        self.spin.simplify_values();
        self.constraints.simplify_each();
    }
}

impl Lattice for MulticenterBondAst {
    fn is_undetermined(&self) -> bool {
        self.electrons.is_undetermined()
            && self.charge.is_undetermined()
            && self.spin.is_undetermined()
            && self.constraints.is_undetermined()
    }

    fn is_ground(&self) -> bool {
        self.electrons.is_ground()
            && self.charge.is_ground()
            && self.spin.is_ground()
            && self.constraints.is_ground()
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        Some(Self {
            electrons: self.electrons.meet(&other.electrons)?,
            charge: self.charge.meet(&other.charge)?,
            spin: self.spin.meet(&other.spin)?,
            constraints: self.constraints.meet(&other.constraints)?,
        })
    }

    fn join(&self, other: &Self) -> Self {
        Self {
            electrons: self.electrons.join(&other.electrons),
            charge: self.charge.join(&other.charge),
            spin: self.spin.join(&other.spin),
            constraints: self.constraints.join(&other.constraints),
        }
    }

    /// `self` (pattern) matches `target` iff `electrons` (whole-vector),
    /// `charge`, `spin`, and `constraints` all match field-wise.
    fn matches(&self, target: &Self) -> bool {
        self.electrons.matches(&target.electrons)
            && self.charge.matches(&target.charge)
            && self.spin.matches(&target.spin)
            && self.constraints.matches(&target.constraints)
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::ast::value::ValueTerm;

    #[rustfmt::skip]
    #[rstest]
    #[case::new(MulticenterBondAst::new(ElectronCountsAst::Lit(vec![1; 3])),
        MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1; 3]),
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: MulticenterBondConstraints::new() })]
    #[case::from_counts(MulticenterBondAst::from_counts(vec![1, 1, 1]),
        MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![1; 3]),
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: MulticenterBondConstraints::new() })]
    fn test_multicenter_bond_ast_new(
        #[case] actual: MulticenterBondAst,
        #[case] expected: MulticenterBondAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_electrons(MulticenterBondAst::default().with_electrons(ElectronCountsAst::Lit(vec![2; 3])),
        MulticenterBondAst { electrons: ElectronCountsAst::Lit(vec![2; 3]),
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: MulticenterBondConstraints::new() })]
    #[case::with_charge(MulticenterBondAst::default().with_charge(-1),
        MulticenterBondAst { electrons: ElectronCountsAst::Undetermined,
            charge: ValueAst::Lit(-1), spin: SpinStateAst::default(),
            constraints: MulticenterBondConstraints::new() })]
    #[case::with_spin(MulticenterBondAst::default().with_spin((0_u8, 1_u8)),
        MulticenterBondAst { electrons: ElectronCountsAst::Undetermined,
            charge: ValueAst::Undetermined, spin: SpinStateAst::closed_shell(),
            constraints: MulticenterBondConstraints::new() })]
    #[case::with_constraint(
        MulticenterBondAst::default().with_constraint(MulticenterBondConstraint::electron_count(2)),
        MulticenterBondAst { electrons: ElectronCountsAst::Undetermined,
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2)) })]
    #[case::with_constraints_extends(
        MulticenterBondAst::default()
            .with_constraints([MulticenterBondConstraint::electron_count(2)]),
        MulticenterBondAst { electrons: ElectronCountsAst::Undetermined,
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2)) })]
    #[case::with_constraint_replaces_same_kind(
        MulticenterBondAst::default()
            .with_constraint(MulticenterBondConstraint::electron_count(2))
            .with_constraint(MulticenterBondConstraint::electron_count(4)),
        MulticenterBondAst { electrons: ElectronCountsAst::Undetermined,
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(4)) })]
    fn test_multicenter_bond_ast_with_methods(
        #[case] actual: MulticenterBondAst,
        #[case] expected: MulticenterBondAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::from_ground_electrons(
        MulticenterBondAst::from_counts(vec![1; 3]).into_ground(),
        MulticenterBondAst {
            electrons: ElectronCountsAst::Lit(vec![1; 3]),
            charge: ValueAst::Lit(0),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            constraints: MulticenterBondConstraints::new(),
        },
    )]
    #[case::preserves_set_charge(
        MulticenterBondAst::from_counts(vec![1; 3]).with_charge(1_i64).into_ground(),
        MulticenterBondAst {
            electrons: ElectronCountsAst::Lit(vec![1; 3]),
            charge: ValueAst::Lit(1),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            constraints: MulticenterBondConstraints::new(),
        },
    )]
    #[case::preserves_constraints(
        MulticenterBondAst::from_counts(vec![1; 3])
            .with_constraint(MulticenterBondConstraint::electron_count(3))
            .into_ground(),
        MulticenterBondAst {
            electrons: ElectronCountsAst::Lit(vec![1; 3]),
            charge: ValueAst::Lit(0),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            constraints: MulticenterBondConstraints::from(
                MulticenterBondConstraint::electron_count(3),
            ),
        },
    )]
    fn test_multicenter_bond_ast_into_ground(
        #[case] actual: MulticenterBondAst,
        #[case] expected: MulticenterBondAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::from_counts(MulticenterBondAst::from_counts(vec![1; 3]))]
    #[case::with_constraint(
        MulticenterBondAst::from_counts(vec![1; 3])
            .with_constraint(MulticenterBondConstraint::electron_count(3))
    )]
    fn test_multicenter_bond_ast_into_zeroed(#[case] bond: MulticenterBondAst) {
        assert_eq!(bond.clone().into_zeroed(), bond.into_ground());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(MulticenterBondAst::default(), false)]
    #[case::charge_only(MulticenterBondAst::new(ElectronCountsAst::Undetermined).with_charge(0), false)]
    #[case::ground_no_atoms(MulticenterBondAst::new(ElectronCountsAst::Lit(Vec::new())).with_charge(0).with_spin((0, 1)), true)]
    #[case::all_ground_three(
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![1; 3])).with_charge(0).with_spin((0, 1)),
        true,
    )]
    #[case::ground_with_constraint(
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![1; 3]))
            .with_charge(0).with_spin((0, 1))
            .with_constraint(MulticenterBondConstraint::electron_count(3)),
        true,
    )]
    fn test_multicenter_bond_ast_is_ground(
        #[case] ast: MulticenterBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_default(MulticenterBondAst::default(), MulticenterBondAst::default(), true)]
    #[case::default_matches_ground(
        MulticenterBondAst::default(),
        MulticenterBondAst::new(ElectronCountsAst::Lit(Vec::new())).with_charge(0).with_spin((0, 1)),
        true,
    )]
    #[case::exact(
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![1; 3])).with_charge(0).with_spin((0, 1)),
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![1; 3])).with_charge(0).with_spin((0, 1)),
        true,
    )]
    #[case::electrons_length_mismatch(
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![1; 2])),
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![1; 3])).with_charge(0).with_spin((0, 1)),
        false,
    )]
    #[case::electrons_value_mismatch(
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![2; 3])),
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![1; 3])).with_charge(0).with_spin((0, 1)),
        false,
    )]
    #[case::charge_mismatch(
        MulticenterBondAst::new(ElectronCountsAst::Undetermined).with_charge(1),
        MulticenterBondAst::new(ElectronCountsAst::Undetermined).with_charge(0),
        false,
    )]
    #[case::spin_mismatch(
        MulticenterBondAst::new(ElectronCountsAst::Undetermined).with_spin((2_u8, 3_u8)),
        MulticenterBondAst::new(ElectronCountsAst::Undetermined).with_spin((0_u8, 1_u8)),
        false,
    )]
    #[case::constraint_required_present(
        MulticenterBondAst::new(ElectronCountsAst::Undetermined)
            .with_constraint(MulticenterBondConstraint::electron_count(3)),
        MulticenterBondAst::new(ElectronCountsAst::Undetermined)
            .with_constraint(MulticenterBondConstraint::electron_count(3)),
        true,
    )]
    #[case::constraint_required_absent(
        MulticenterBondAst::new(ElectronCountsAst::Undetermined)
            .with_constraint(MulticenterBondConstraint::electron_count(3)),
        MulticenterBondAst::new(ElectronCountsAst::Undetermined),
        false,
    )]
    fn test_multicenter_bond_ast_matches(
        #[case] pattern: MulticenterBondAst,
        #[case] target: MulticenterBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    fn test_multicenter_bond_ast_simplify_values() {
        let mut bond = MulticenterBondAst {
            electrons: ElectronCountsAst::Lit(vec![1]),
            charge: ValueAst::term(ValueTerm::Lit(0)),
            spin: SpinStateAst::default(),
            constraints: MulticenterBondConstraints::from(
                MulticenterBondConstraint::ElectronCount(ValueAst::term(ValueTerm::Lit(2))),
            ),
        };
        bond.simplify_values();
        assert_eq!(
            bond,
            MulticenterBondAst {
                electrons: ElectronCountsAst::Lit(vec![1]),
                charge: ValueAst::Lit(0),
                spin: SpinStateAst::default(),
                constraints: MulticenterBondConstraints::from(
                    MulticenterBondConstraint::electron_count(2),
                ),
            },
        );
    }

    #[rstest]
    #[case::both_default(
        MulticenterBondAst::default(),
        MulticenterBondAst::default(),
        Some(MulticenterBondAst::default())
    )]
    #[case::electrons_length_mismatch(
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![2; 3])),
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![2; 4])),
        None,
    )]
    #[case::narrows_electrons(
        MulticenterBondAst::new(ElectronCountsAst::Undetermined),
        MulticenterBondAst::from_counts(vec![1, 2]),
        Some(MulticenterBondAst::from_counts(vec![1, 2])),
    )]
    fn test_multicenter_bond_ast_meet(
        #[case] a: MulticenterBondAst,
        #[case] b: MulticenterBondAst,
        #[case] expected: Option<MulticenterBondAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rstest]
    #[case::electrons_length_mismatch_widens_to_default(
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![2; 3])),
        MulticenterBondAst::new(ElectronCountsAst::Lit(vec![2; 4])),
        MulticenterBondAst::default(),
    )]
    fn test_multicenter_bond_ast_join(
        #[case] a: MulticenterBondAst,
        #[case] b: MulticenterBondAst,
        #[case] expected: MulticenterBondAst,
    ) {
        assert_eq!(a.join(&b), expected);
    }
}
