//! Multicenter bond AST.

use std::mem;

use super::constraint::MulticenterBondConstraints;
use super::spin::SpinStateAst;
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
        Self::new(electrons.into_iter().map(|n| ValueAst::Lit(n as i64)).collect())
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

    pub fn with_constraints(mut self, constraints: impl Into<MulticenterBondConstraints>) -> Self {
        self.constraints = constraints.into();
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
    fn test_multicenter_bond_ast_into_ground(
        #[case] actual: MulticenterBondAst,
        #[case] expected: MulticenterBondAst,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    fn test_multicenter_bond_ast_into_zeroed() {
        let bond = MulticenterBondAst::from_electrons(vec![1; 3]);
        assert_eq!(bond.clone().into_zeroed(), bond.into_ground());
    }
}
