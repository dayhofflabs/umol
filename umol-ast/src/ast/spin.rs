//! Spin-state AST.

use std::mem;

use umol_shared::spin::{SpinMultiplicity, SpinState};

use super::traits::{AsLit, Lattice};
use super::value::ValueAst;

/// Spin state: unpaired-electron count and multiplicity as independent
/// `ValueAst` fields. Both may be `Undetermined`, a literal, or an
/// expression pattern.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SpinStateAst {
    pub unpaired: ValueAst,
    pub multiplicity: ValueAst,
}

impl SpinStateAst {
    pub fn closed_shell() -> Self {
        SpinState::closed_shell().into()
    }

    /// Simplify both `unpaired` and `multiplicity` in place via
    /// [`ValueAst::simplify`].
    pub fn simplify_values(&mut self) {
        self.unpaired = mem::take(&mut self.unpaired).simplify();
        self.multiplicity = mem::take(&mut self.multiplicity).simplify();
    }
}

impl Lattice for SpinStateAst {
    /// Both fields are `Undetermined` — the spin state asserts nothing.
    fn is_undetermined(&self) -> bool {
        self.unpaired.is_undetermined() && self.multiplicity.is_undetermined()
    }

    /// Both fields are ground. A ground spin state may still be physically
    /// inconsistent — parity of `(unpaired, multiplicity)` is a tier-2
    /// physics invariant enforced by a propagator in the solver, not here.
    fn is_ground(&self) -> bool {
        self.unpaired.is_ground() && self.multiplicity.is_ground()
    }

    /// Field-wise meet — no physics validation. A meet that combines two
    /// individually-physics-consistent inputs may produce a parity-invalid
    /// pair; the solver's propagator catches this separately.
    fn meet(&self, other: &Self) -> Option<Self> {
        Some(Self {
            unpaired: self.unpaired.meet(&other.unpaired)?,
            multiplicity: self.multiplicity.meet(&other.multiplicity)?,
        })
    }

    fn join(&self, other: &Self) -> Self {
        Self {
            unpaired: self.unpaired.join(&other.unpaired),
            multiplicity: self.multiplicity.join(&other.multiplicity),
        }
    }

    /// Pattern matches target iff `unpaired` and `multiplicity` each
    /// match field-wise under `ValueAst::matches`.
    fn matches(&self, target: &Self) -> bool {
        self.unpaired.matches(&target.unpaired) && self.multiplicity.matches(&target.multiplicity)
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
        Self {
            unpaired: ValueAst::Lit(unpaired as i64),
            multiplicity: ValueAst::Lit(multiplicity as i64),
        }
    }
}

impl From<SpinState> for SpinStateAst {
    fn from(state: SpinState) -> Self {
        (state.unpaired(), u8::from(state.multiplicity())).into()
    }
}

impl AsLit for SpinStateAst {
    type Lit = SpinState;

    /// Concrete [`SpinState`] when both fields are `Lit` *and* the
    /// `(unpaired, multiplicity)` pair satisfies physics parity. Strictly
    /// narrower than [`is_ground`](SpinStateAst::is_ground), which checks
    /// only that both fields are literal.
    #[inline]
    fn as_lit(&self) -> Option<SpinState> {
        let (ValueAst::Lit(u), ValueAst::Lit(m)) = (&self.unpaired, &self.multiplicity) else {
            return None;
        };
        let mult = SpinMultiplicity::from_repr(*m as u8)?;
        SpinState::try_new(*u as u8, mult).ok()
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
    #[case(spin!("#u2").into(), true)]
    #[case((2, 3).into(), true)]
    #[case(SpinStateAst { unpaired: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined }, false)]
    fn test_spin_state_ast_is_ground(#[case] ast: SpinStateAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::both_undetermined(SpinStateAst::default(), (2_u8, 3_u8).into(), true)]
    #[case::pattern_specific_target_undetermined((2_u8, 3_u8).into(), SpinStateAst::default(), false)]
    #[case::exact((2_u8, 3_u8).into(), (2_u8, 3_u8).into(), true)]
    #[case::unpaired_mismatch((2_u8, 3_u8).into(), (0_u8, 3_u8).into(), false)]
    #[case::partial_pattern(
        SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(3) },
        (2_u8, 3_u8).into(),
        true,
    )]
    fn test_spin_state_ast_matches(
        #[case] pattern: SpinStateAst,
        #[case] target: SpinStateAst,
        #[case] expected: bool,
    ) {
        assert_eq!(pattern.matches(&target), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::both_undetermined(SpinStateAst::default(), None)]
    #[case::unpaired_lit_only(
        SpinStateAst { unpaired: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined },
        None,
    )]
    #[case::multiplicity_lit_only(
        SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(3) },
        None,
    )]
    #[case::valid_triplet((2_u8, 3_u8).into(), Some(spin!("#u2#s3")))]
    #[case::valid_closed_shell((0_u8, 1_u8).into(), Some(SpinState::closed_shell()))]
    #[case::parity_invalid((1_u8, 1_u8).into(), None)]
    fn test_spin_state_ast_as_lit(
        #[case] ast: SpinStateAst,
        #[case] expected: Option<SpinState>,
    ) {
        assert_eq!(ast.as_lit(), expected);
    }

    #[test]
    fn test_spin_state_ast_from_spin_state() {
        let ast: SpinStateAst = spin!("#u2").into();
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
        let ast: SpinStateAst = s.into();
        assert_eq!(ast, SpinStateAst::closed_shell());
    }
}
