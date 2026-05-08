//! Spin-state AST.

use std::mem;

use umol_shared::spin::SpinState;

use super::value::ValueAst;

/// Spin state: unpaired-electron count and multiplicity as independent
/// `ValueAst` fields. Both may be `Undetermined`, a literal, or an
/// expression pattern.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct SpinStateAst {
    pub unpaired: ValueAst,
    pub multiplicity: ValueAst,
}

impl SpinStateAst {
    pub fn new(unpaired: u8, multiplicity: u8) -> Self {
        Self::from_values(
            ValueAst::Lit(unpaired as i64),
            ValueAst::Lit(multiplicity as i64),
        )
    }

    pub fn from_values(unpaired: ValueAst, multiplicity: ValueAst) -> Self {
        Self {
            unpaired,
            multiplicity,
        }
    }

    pub fn from_state(state: SpinState) -> Self {
        Self::new(state.unpaired(), state.multiplicity().multiplicity())
    }

    pub fn closed_shell() -> Self {
        Self::from_state(SpinState::closed_shell())
    }

    /// Both fields are literal. A ground spin state may still be physically
    /// inconsistent — parity of `(unpaired, multiplicity)` is a tier-2
    /// physics invariant enforced by a propagator in the solver, not here.
    pub fn is_ground(&self) -> bool {
        matches!(&self.unpaired, ValueAst::Lit(_)) && matches!(&self.multiplicity, ValueAst::Lit(_))
    }

    /// Both fields are `Undetermined` — the spin state asserts nothing.
    pub fn is_undetermined(&self) -> bool {
        self.unpaired.is_undetermined() && self.multiplicity.is_undetermined()
    }

    /// Pattern matches target iff `unpaired` and `multiplicity` each
    /// match field-wise under `ValueAst::matches`.
    pub fn matches(&self, target: &Self) -> bool {
        self.unpaired.matches(&target.unpaired) && self.multiplicity.matches(&target.multiplicity)
    }

    /// Simplify both `unpaired` and `multiplicity` in place via
    /// [`ValueAst::simplify`].
    pub fn simplify_values(&mut self) {
        self.unpaired = mem::take(&mut self.unpaired).simplify();
        self.multiplicity = mem::take(&mut self.multiplicity).simplify();
    }
}

impl Default for SpinStateAst {
    fn default() -> Self {
        Self {
            unpaired: ValueAst::Undetermined,
            multiplicity: ValueAst::Undetermined,
        }
    }
}

impl From<(u8, u8)> for SpinStateAst {
    fn from((unpaired, multiplicity): (u8, u8)) -> Self {
        Self::new(unpaired, multiplicity)
    }
}

impl From<SpinState> for SpinStateAst {
    fn from(state: SpinState) -> Self {
        Self::from_state(state)
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_shared::spin;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case(SpinStateAst::default(), false)]
    #[case(SpinStateAst::from_state(spin!("#u2")), true)]
    #[case(SpinStateAst::new(2, 3), true)]
    #[case(SpinStateAst::from_values(ValueAst::Lit(2), ValueAst::Undetermined), false)]
    fn test_spin_state_ast_is_ground(#[case] ast: SpinStateAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::both_undetermined(SpinStateAst::default(), SpinStateAst::new(2, 3), true)]
    #[case::pattern_specific_target_undetermined(SpinStateAst::new(2, 3), SpinStateAst::default(), false)]
    #[case::exact(SpinStateAst::new(2, 3), SpinStateAst::new(2, 3), true)]
    #[case::unpaired_mismatch(SpinStateAst::new(2, 3), SpinStateAst::new(0, 3), false)]
    #[case::partial_pattern(SpinStateAst::from_values(ValueAst::Undetermined, ValueAst::Lit(3)), SpinStateAst::new(2, 3), true)]
    fn test_spin_state_ast_matches(
        #[case] pattern: SpinStateAst,
        #[case] target: SpinStateAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[test]
    fn test_spin_state_ast_from_state() {
        let ast = SpinStateAst::from_state(spin!("#u2"));
        assert_eq!(ast.unpaired, ValueAst::Lit(2));
        assert_eq!(ast.multiplicity, ValueAst::Lit(3));
    }

    #[test]
    fn test_spin_state_ast_closed_shell() {
        let ast = SpinStateAst::closed_shell();
        assert_eq!(ast.unpaired, ValueAst::Lit(0));
        assert_eq!(ast.multiplicity, ValueAst::Lit(1));
    }

    #[test]
    fn test_spin_state_from_str_roundtrip() {
        let s = SpinState::from_str("#u0#s1").unwrap();
        assert_eq!(SpinStateAst::from_state(s), SpinStateAst::closed_shell());
    }
}
