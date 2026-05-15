//! Aromatic system AST.

use std::mem;

use super::constraint::{AromaticSystemConstraint, AromaticSystemConstraints};
use super::spin::SpinStateAst;
use super::traits::Lattice;
use super::value::ValueAst;

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct AromaticSystemAst {
    pub electrons: Vec<ValueAst>,
    pub charge: ValueAst,
    pub spin: SpinStateAst,
    pub constraints: AromaticSystemConstraints,
}

impl AromaticSystemAst {
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

    /// Equivalent to `into_ground()`. `AromaticSystemAst` has no constraint
    /// defaults.
    pub fn into_zeroed(self) -> Self {
        self.into_ground()
    }

    /// `self` (pattern) matches `target` iff per-atom electrons match
    /// position-wise (length-equality required) and `charge` / `spin` match
    /// field-wise.
    pub fn matches(&self, target: &AromaticSystemAst) -> bool {
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

impl Lattice for AromaticSystemAst {
    fn is_undetermined(&self) -> bool {
        self.electrons.iter().all(|v| v.is_undetermined())
            && self.charge.is_undetermined()
            && self.spin.is_undetermined()
            && self.constraints.is_undetermined()
    }

    fn is_ground(&self) -> bool {
        self.electrons.iter().all(|v| v.is_ground())
            && self.charge.is_ground()
            && self.spin.is_ground()
            && self.constraints.is_ground()
    }

    fn meet(&self, other: &Self) -> Option<Self> {
        if self.electrons.len() != other.electrons.len() {
            return None;
        }
        let electrons: Option<Vec<ValueAst>> = self
            .electrons
            .iter()
            .zip(&other.electrons)
            .map(|(a, b)| a.meet(b))
            .collect();
        Some(Self {
            electrons: electrons?,
            charge: self.charge.meet(&other.charge)?,
            spin: self.spin.meet(&other.spin)?,
            constraints: self.constraints.meet(&other.constraints)?,
        })
    }

    fn join(&self, other: &Self) -> Self {
        if self.electrons.len() != other.electrons.len() {
            return Self::default();
        }
        Self {
            electrons: self
                .electrons
                .iter()
                .zip(&other.electrons)
                .map(|(a, b)| a.join(b))
                .collect(),
            charge: self.charge.join(&other.charge),
            spin: self.spin.join(&other.spin),
            constraints: self.constraints.join(&other.constraints),
        }
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
    #[case::new(AromaticSystemAst::new(vec![ValueAst::Lit(1); 3]),
        AromaticSystemAst { electrons: vec![ValueAst::Lit(1); 3],
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: AromaticSystemConstraints::new() })]
    #[case::from_electrons(AromaticSystemAst::from_electrons(vec![1, 1, 1]),
        AromaticSystemAst { electrons: vec![ValueAst::Lit(1); 3],
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: AromaticSystemConstraints::new() })]
    fn test_aromatic_system_ast_new(
        #[case] actual: AromaticSystemAst,
        #[case] expected: AromaticSystemAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::with_electrons(AromaticSystemAst::default().with_electrons(vec![ValueAst::Lit(2); 3]),
        AromaticSystemAst { electrons: vec![ValueAst::Lit(2); 3],
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: AromaticSystemConstraints::new() })]
    #[case::with_charge(AromaticSystemAst::default().with_charge(-1),
        AromaticSystemAst { electrons: Vec::new(),
            charge: ValueAst::Lit(-1), spin: SpinStateAst::default(),
            constraints: AromaticSystemConstraints::new() })]
    #[case::with_spin(AromaticSystemAst::default().with_spin((0_u8, 1_u8)),
        AromaticSystemAst { electrons: Vec::new(),
            charge: ValueAst::Undetermined, spin: SpinStateAst::closed_shell(),
            constraints: AromaticSystemConstraints::new() })]
    #[case::with_constraint(
        AromaticSystemAst::default().with_constraint(AromaticSystemConstraint::electron_count(6)),
        AromaticSystemAst { electrons: Vec::new(),
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6)) })]
    #[case::with_constraints_extends(
        AromaticSystemAst::default()
            .with_constraints([AromaticSystemConstraint::electron_count(6)]),
        AromaticSystemAst { electrons: Vec::new(),
            charge: ValueAst::Undetermined, spin: SpinStateAst::default(),
            constraints: AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6)) })]
    #[case::with_constraint_replaces_same_kind(
        AromaticSystemAst::default()
            .with_constraint(AromaticSystemConstraint::electron_count(2))
            .with_constraint(AromaticSystemConstraint::electron_count(6)),
        AromaticSystemAst { electrons: Vec::new(),
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
        AromaticSystemAst { electrons: vec![ValueAst::Lit(1); 6], charge: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AromaticSystemConstraints::new() })]
    #[case::preserves_set_charge(AromaticSystemAst::from_electrons(vec![1; 6]).with_charge(1_i64).into_ground(),
        AromaticSystemAst { electrons: vec![ValueAst::Lit(1); 6], charge: ValueAst::Lit(1), spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AromaticSystemConstraints::new() })]
    #[case::preserves_constraints(AromaticSystemAst::from_electrons(vec![1; 6]).with_constraint(AromaticSystemConstraint::electron_count(6)).into_ground(),
        AromaticSystemAst { electrons: vec![ValueAst::Lit(1); 6], charge: ValueAst::Lit(0), spin: SpinStateAst::from((0_u8, 1_u8)),
        constraints: AromaticSystemConstraints::from(AromaticSystemConstraint::electron_count(6)) })]
    fn test_aromatic_system_ast_into_ground(
        #[case] actual: AromaticSystemAst,
        #[case] expected: AromaticSystemAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::from_electrons(AromaticSystemAst::from_electrons(vec![1; 6]))]
    #[case::with_constraint(AromaticSystemAst::from_electrons(vec![1; 6]).with_constraint(AromaticSystemConstraint::electron_count(6)))]
    fn test_aromatic_system_ast_into_zeroed(#[case] system: AromaticSystemAst) {
        assert_eq!(system.clone().into_zeroed(), system.into_ground());
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::all_undetermined(AromaticSystemAst::default(), false)]
    #[case::charge_only(AromaticSystemAst::new(Vec::new()).with_charge(0), false)]
    #[case::ground_no_atoms(AromaticSystemAst::new(Vec::new()).with_charge(0).with_spin((0, 1)), true)]
    #[case::all_ground_six(AromaticSystemAst::new(vec![ValueAst::Lit(1); 6]).with_charge(0).with_spin((0, 1)), true)]
    #[case::one_undetermined_electron(AromaticSystemAst::new(vec![ValueAst::Lit(1), ValueAst::Undetermined, ValueAst::Lit(1)]).with_charge(0).with_spin((0, 1)), false)]
    #[case::ground_with_constraint(AromaticSystemAst::new(vec![ValueAst::Lit(1); 6]).with_charge(0).with_spin((0, 1)).with_constraint(AromaticSystemConstraint::electron_count(6)), true)]
    fn test_aromatic_system_ast_is_ground(
        #[case] ast: AromaticSystemAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_default(AromaticSystemAst::default(), AromaticSystemAst::default(), true)]
    #[case::default_matches_ground(AromaticSystemAst::default(), AromaticSystemAst::new(Vec::new()).with_charge(0).with_spin((0, 1)), true)]
    #[case::exact(AromaticSystemAst::new(vec![ValueAst::Lit(1); 6]).with_charge(0).with_spin((0, 1)),
        AromaticSystemAst::new(vec![ValueAst::Lit(1); 6]).with_charge(0).with_spin((0, 1)), true)]
    #[case::electrons_length_mismatch(AromaticSystemAst::new(vec![ValueAst::Lit(1); 5]),
        AromaticSystemAst::new(vec![ValueAst::Lit(1); 6]).with_charge(0).with_spin((0, 1)), false)]
    #[case::electrons_value_mismatch(AromaticSystemAst::new(vec![ValueAst::Lit(2); 6]),
        AromaticSystemAst::new(vec![ValueAst::Lit(1); 6]).with_charge(0).with_spin((0, 1)), false)]
    #[case::pattern_undetermined_electron_matches_lit(AromaticSystemAst::new(vec![ValueAst::Undetermined; 6]),
      AromaticSystemAst::new(vec![ValueAst::Lit(1); 6]).with_charge(0).with_spin((0, 1)), true)]
    fn test_aromatic_system_ast_matches(
        #[case] pattern: AromaticSystemAst,
        #[case] target: AromaticSystemAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rstest]
    fn test_aromatic_system_ast_simplify_values() {
        let mut system = AromaticSystemAst {
            electrons: vec![ValueAst::Expr(Box::new(Expr::Lit(1)))],
            charge: ValueAst::Expr(Box::new(Expr::Lit(0))),
            spin: SpinStateAst::default(),
            constraints: AromaticSystemConstraints::from(AromaticSystemConstraint::ElectronCount(
                ValueAst::Expr(Box::new(Expr::Lit(6))),
            )),
        };
        system.simplify_values();
        assert_eq!(
            system,
            AromaticSystemAst {
                electrons: vec![ValueAst::Lit(1)],
                charge: ValueAst::Lit(0),
                spin: SpinStateAst::default(),
                constraints: AromaticSystemConstraints::from(
                    AromaticSystemConstraint::electron_count(6),
                ),
            },
        );
    }

    #[rstest]
    #[case::both_default(
        AromaticSystemAst::default(),
        AromaticSystemAst::default(),
        Some(AromaticSystemAst::default()),
    )]
    #[case::electrons_length_mismatch(
        AromaticSystemAst::from_electrons(vec![1; 6]),
        AromaticSystemAst::from_electrons(vec![1; 5]),
        None,
    )]
    #[case::narrows_electrons(
        AromaticSystemAst::new(vec![ValueAst::Undetermined; 3]),
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
}
