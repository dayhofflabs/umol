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
    pub fn new(electrons: Vec<ValueAst>, charge: ValueAst, spin: SpinStateAst) -> Self {
        Self {
            electrons,
            charge,
            spin,
            constraints: MulticenterBondConstraints::new(),
        }
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
    #[case::charge_only(MulticenterBondAst::new(Vec::new(), ValueAst::Lit(0), SpinStateAst::default()), false)]
    #[case::ground_no_atoms(MulticenterBondAst::new(Vec::new(), ValueAst::Lit(0), SpinStateAst::new(0, 1)), true)]
    #[case::all_ground_three(
        MulticenterBondAst::new(
            vec![ValueAst::Lit(1); 3],
            ValueAst::Lit(0),
            SpinStateAst::new(0, 1),
        ),
        true,
    )]
    #[case::one_undetermined_electron(
        MulticenterBondAst::new(
            vec![ValueAst::Lit(1), ValueAst::Undetermined, ValueAst::Lit(1)],
            ValueAst::Lit(0),
            SpinStateAst::new(0, 1),
        ),
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
        MulticenterBondAst::new(Vec::new(), ValueAst::Lit(0), SpinStateAst::new(0, 1)),
        true,
    )]
    #[case::exact(
        MulticenterBondAst::new(vec![ValueAst::Lit(1); 3], ValueAst::Lit(0), SpinStateAst::new(0, 1)),
        MulticenterBondAst::new(vec![ValueAst::Lit(1); 3], ValueAst::Lit(0), SpinStateAst::new(0, 1)),
        true,
    )]
    #[case::electrons_length_mismatch(
        MulticenterBondAst::new(vec![ValueAst::Lit(1); 2], ValueAst::Undetermined, SpinStateAst::default()),
        MulticenterBondAst::new(vec![ValueAst::Lit(1); 3], ValueAst::Lit(0), SpinStateAst::new(0, 1)),
        false,
    )]
    #[case::electrons_value_mismatch(
        MulticenterBondAst::new(vec![ValueAst::Lit(2); 3], ValueAst::Undetermined, SpinStateAst::default()),
        MulticenterBondAst::new(vec![ValueAst::Lit(1); 3], ValueAst::Lit(0), SpinStateAst::new(0, 1)),
        false,
    )]
    #[case::charge_mismatch(
        MulticenterBondAst::new(vec![ValueAst::Undetermined; 3], ValueAst::Lit(1), SpinStateAst::default()),
        MulticenterBondAst::new(vec![ValueAst::Undetermined; 3], ValueAst::Lit(0), SpinStateAst::default()),
        false,
    )]
    fn test_multicenter_bond_ast_matches(
        #[case] pattern: MulticenterBondAst,
        #[case] target: MulticenterBondAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }
}
