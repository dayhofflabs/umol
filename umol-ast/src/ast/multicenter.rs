//! Multicenter bond AST.

use std::mem;

use super::constraint::{MulticenterBondConstraint, MulticenterBondConstraints};
use super::spin::SpinStateAst;
use super::traits::Lattice;
use super::value::ValueAst;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct MulticenterBondAst {
    pub electrons: Vec<ValueAst>,
    pub charge: ValueAst,
    pub spin: SpinStateAst,
    pub constraints: MulticenterBondConstraints,
}

impl MulticenterBondAst {
    pub fn new(electrons: Vec<ValueAst>) -> Self {
        Self {
            electrons,
            ..Default::default()
        }
    }

    pub fn from_electrons(electrons: Vec<u8>) -> Self {
        Self::new(
            electrons
                .into_iter()
                .map(|n| ValueAst::Lit(n as i64))
                .collect(),
        )
    }

    pub fn with_electrons(mut self, electrons: Vec<ValueAst>) -> Self {
        self.electrons = electrons;
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
    /// charge → `Lit(0)`, spin → closed-shell singlet `(0, 1)`. Per-atom
    /// `electrons` entries and `constraints` are preserved. The result is
    /// ground iff every `electrons` entry is already ground.
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

    pub fn is_ground(&self) -> bool {
        self.charge.is_ground()
            && self.spin.is_ground()
            && self.electrons.iter().all(|v| v.is_ground())
    }

    /// `self` (pattern) matches `target` iff per-atom electrons match
    /// position-wise (length-equality required) and `charge` / `spin` match
    /// field-wise.
    pub fn matches(&self, target: &MulticenterBondAst) -> bool {
        self.charge.matches(&target.charge)
            && self.spin.matches(&target.spin)
            && self.electrons.len() == target.electrons.len()
            && self
                .electrons
                .iter()
                .zip(&target.electrons)
                .all(|(p, t)| p.matches(t))
    }

    pub fn simplify_values(&mut self) {
        self.charge = mem::take(&mut self.charge).simplify();
        self.spin.simplify_values();
        for e in self.electrons.iter_mut() {
            *e = mem::take(e).simplify();
        }
        self.constraints.simplify_each();
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
    #[case::new(MulticenterBondAst::new(vec![ValueAst::Lit(1); 3]),
        MulticenterBondAst { electrons: vec![ValueAst::Lit(1); 3],
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: MulticenterBondConstraints::new() })]
    #[case::from_electrons(MulticenterBondAst::from_electrons(vec![1, 1, 1]),
        MulticenterBondAst { electrons: vec![ValueAst::Lit(1); 3],
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
    #[case::with_electrons(MulticenterBondAst::default().with_electrons(vec![ValueAst::Lit(2); 3]),
        MulticenterBondAst { electrons: vec![ValueAst::Lit(2); 3],
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: MulticenterBondConstraints::new() })]
    #[case::with_charge(MulticenterBondAst::default().with_charge(-1),
        MulticenterBondAst { electrons: Vec::new(),
            charge: ValueAst::Lit(-1), spin: SpinStateAst::default(),
            constraints: MulticenterBondConstraints::new() })]
    #[case::with_spin(MulticenterBondAst::default().with_spin((0_u8, 1_u8)),
        MulticenterBondAst { electrons: Vec::new(),
            charge: ValueAst::Undetermined, spin: SpinStateAst::closed_shell(),
            constraints: MulticenterBondConstraints::new() })]
    #[case::with_constraint(
        MulticenterBondAst::default().with_constraint(MulticenterBondConstraint::electron_count(2)),
        MulticenterBondAst { electrons: Vec::new(),
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2)) })]
    #[case::with_constraints_extends(
        MulticenterBondAst::default()
            .with_constraints([MulticenterBondConstraint::electron_count(2)]),
        MulticenterBondAst { electrons: Vec::new(),
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: MulticenterBondConstraints::from(MulticenterBondConstraint::electron_count(2)) })]
    #[case::with_constraint_replaces_same_kind(
        MulticenterBondAst::default()
            .with_constraint(MulticenterBondConstraint::electron_count(2))
            .with_constraint(MulticenterBondConstraint::electron_count(4)),
        MulticenterBondAst { electrons: Vec::new(),
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
        MulticenterBondAst::from_electrons(vec![1; 3]).into_ground(),
        MulticenterBondAst {
            electrons: vec![ValueAst::Lit(1); 3],
            charge: ValueAst::Lit(0),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            constraints: MulticenterBondConstraints::new(),
        },
    )]
    #[case::preserves_set_charge(
        MulticenterBondAst::from_electrons(vec![1; 3]).with_charge(1_i64).into_ground(),
        MulticenterBondAst {
            electrons: vec![ValueAst::Lit(1); 3],
            charge: ValueAst::Lit(1),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            constraints: MulticenterBondConstraints::new(),
        },
    )]
    #[case::preserves_constraints(
        MulticenterBondAst::from_electrons(vec![1; 3])
            .with_constraint(MulticenterBondConstraint::electron_count(3))
            .into_ground(),
        MulticenterBondAst {
            electrons: vec![ValueAst::Lit(1); 3],
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
    #[case::from_electrons(MulticenterBondAst::from_electrons(vec![1; 3]))]
    #[case::with_constraint(
        MulticenterBondAst::from_electrons(vec![1; 3])
            .with_constraint(MulticenterBondConstraint::electron_count(3))
    )]
    fn test_multicenter_bond_ast_into_zeroed(#[case] bond: MulticenterBondAst) {
        assert_eq!(bond.clone().into_zeroed(), bond.into_ground());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_(MulticenterBondAst::default(), false)]
    #[case::charge_only(MulticenterBondAst::new(Vec::new()).with_charge(0), false)]
    #[case::ground_no_atoms(MulticenterBondAst::new(Vec::new()).with_charge(0).with_spin((0, 1)), true)]
    #[case::all_ground_three(
        MulticenterBondAst::new(vec![ValueAst::Lit(1); 3]).with_charge(0).with_spin((0, 1)),
        true,
    )]
    #[case::one_undetermined_electron(
        MulticenterBondAst::new(vec![ValueAst::Lit(1), ValueAst::Undetermined, ValueAst::Lit(1)])
            .with_charge(0).with_spin((0, 1)),
        false,
    )]
    #[case::ground_with_constraint(
        MulticenterBondAst::new(vec![ValueAst::Lit(1); 3])
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
        MulticenterBondAst::new(Vec::new()).with_charge(0).with_spin((0, 1)),
        true,
    )]
    #[case::exact(
        MulticenterBondAst::new(vec![ValueAst::Lit(1); 3]).with_charge(0).with_spin((0, 1)),
        MulticenterBondAst::new(vec![ValueAst::Lit(1); 3]).with_charge(0).with_spin((0, 1)),
        true,
    )]
    #[case::electrons_length_mismatch(
        MulticenterBondAst::new(vec![ValueAst::Lit(1); 2]),
        MulticenterBondAst::new(vec![ValueAst::Lit(1); 3]).with_charge(0).with_spin((0, 1)),
        false,
    )]
    #[case::electrons_value_mismatch(
        MulticenterBondAst::new(vec![ValueAst::Lit(2); 3]),
        MulticenterBondAst::new(vec![ValueAst::Lit(1); 3]).with_charge(0).with_spin((0, 1)),
        false,
    )]
    #[case::charge_mismatch(
        MulticenterBondAst::new(vec![ValueAst::Undetermined; 3]).with_charge(1),
        MulticenterBondAst::new(vec![ValueAst::Undetermined; 3]).with_charge(0),
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
            electrons: vec![ValueAst::Expr(Box::new(Expr::Lit(1)))],
            charge: ValueAst::Expr(Box::new(Expr::Lit(0))),
            spin: SpinStateAst::default(),
            constraints: MulticenterBondConstraints::from(
                MulticenterBondConstraint::ElectronCount(ValueAst::Expr(Box::new(Expr::Lit(2)))),
            ),
        };
        bond.simplify_values();
        assert_eq!(
            bond,
            MulticenterBondAst {
                electrons: vec![ValueAst::Lit(1)],
                charge: ValueAst::Lit(0),
                spin: SpinStateAst::default(),
                constraints: MulticenterBondConstraints::from(
                    MulticenterBondConstraint::electron_count(2),
                ),
            },
        );
    }
}
