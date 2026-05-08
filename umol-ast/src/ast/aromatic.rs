//! Aromatic system AST.

use std::mem;

use super::constraint::AromaticSystemConstraints;
use super::spin::SpinStateAst;
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

    pub fn with_constraints(mut self, constraints: impl Into<AromaticSystemConstraints>) -> Self {
        self.constraints = constraints.into();
        self
    }

    /// Fill `Undetermined` value-bearing fields with zero defaults: charge
    /// to `Lit(0)`, spin to closed-shell singlet `(0, 1)`. Per-atom
    /// `electrons` entries and `constraints` are preserved. The result is
    /// ground iff every `electrons` entry is already ground.
    pub fn zeroed(mut self) -> Self {
        if self.charge.is_undetermined() {
            self.charge = ValueAst::Lit(0);
        }
        if self.spin.is_undetermined() {
            self.spin = SpinStateAst::from((0_u8, 1_u8));
        }
        self
    }

    pub fn is_ground(&self) -> bool {
        self.charge.is_ground()
            && self.spin.is_ground()
            && self.electrons.iter().all(|v| v.is_ground())
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

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::all_undetermined(AromaticSystemAst::default(), false)]
    #[case::charge_only(AromaticSystemAst::new(Vec::new()).with_charge(0), false)]
    #[case::ground_no_atoms(AromaticSystemAst::new(Vec::new()).with_charge(0).with_spin((0, 1)), true)]
    #[case::all_ground_six(
        AromaticSystemAst::new(vec![ValueAst::Lit(1); 6]).with_charge(0).with_spin((0, 1)),
        true,
    )]
    #[case::one_undetermined_electron(
        AromaticSystemAst::new(vec![ValueAst::Lit(1), ValueAst::Undetermined, ValueAst::Lit(1)])
            .with_charge(0).with_spin((0, 1)),
        false,
    )]
    fn test_aromatic_system_ast_is_ground(
        #[case] ast: AromaticSystemAst,
        #[case] expected: bool,
    ) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_matches_default(AromaticSystemAst::default(), AromaticSystemAst::default(), true)]
    #[case::default_matches_ground(
        AromaticSystemAst::default(),
        AromaticSystemAst::new(Vec::new()).with_charge(0).with_spin((0, 1)),
        true,
    )]
    #[case::exact(
        AromaticSystemAst::new(vec![ValueAst::Lit(1); 6]).with_charge(0).with_spin((0, 1)),
        AromaticSystemAst::new(vec![ValueAst::Lit(1); 6]).with_charge(0).with_spin((0, 1)),
        true,
    )]
    #[case::electrons_length_mismatch(
        AromaticSystemAst::new(vec![ValueAst::Lit(1); 5]),
        AromaticSystemAst::new(vec![ValueAst::Lit(1); 6]).with_charge(0).with_spin((0, 1)),
        false,
    )]
    #[case::electrons_value_mismatch(
        AromaticSystemAst::new(vec![ValueAst::Lit(2); 6]),
        AromaticSystemAst::new(vec![ValueAst::Lit(1); 6]).with_charge(0).with_spin((0, 1)),
        false,
    )]
    #[case::pattern_undetermined_electron_matches_lit(
        AromaticSystemAst::new(vec![ValueAst::Undetermined; 6]),
        AromaticSystemAst::new(vec![ValueAst::Lit(1); 6]).with_charge(0).with_spin((0, 1)),
        true,
    )]
    fn test_aromatic_system_ast_matches(
        #[case] pattern: AromaticSystemAst,
        #[case] target: AromaticSystemAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::from_ground_electrons(
        AromaticSystemAst::from_electrons(vec![1; 6]).zeroed(),
        AromaticSystemAst {
            electrons: vec![ValueAst::Lit(1); 6],
            charge: ValueAst::Lit(0),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            constraints: AromaticSystemConstraints::new(),
        },
    )]
    #[case::preserves_set_charge(
        AromaticSystemAst::from_electrons(vec![1; 6]).with_charge(1_i64).zeroed(),
        AromaticSystemAst {
            electrons: vec![ValueAst::Lit(1); 6],
            charge: ValueAst::Lit(1),
            spin: SpinStateAst::from((0_u8, 1_u8)),
            constraints: AromaticSystemConstraints::new(),
        },
    )]
    fn test_aromatic_system_ast_zeroed(
        #[case] actual: AromaticSystemAst,
        #[case] expected: AromaticSystemAst,
    ) {
        assert_eq!(actual, expected);
    }
}
