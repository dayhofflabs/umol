//! Spin-state AST.

use umol_ast_macros::{Canonicalize, Lattice};
use umol_shared::spin::{SpinMultiplicity, SpinState};

use super::traits::AsLit;
use super::value::ValueAst;

/// Spin state: unpaired-electron count and multiplicity as independent `ValueAst` fields.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Lattice, Canonicalize)]
pub struct SpinStateAst {
    pub unpaired: ValueAst,
    pub multiplicity: ValueAst,
}

impl SpinStateAst {
    pub fn closed_shell() -> Self {
        SpinState::closed_shell().into()
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
    use crate::ast::error::Contradiction;
    use crate::ast::traits::{Canonicalize, Lattice};
    use crate::ast::value::ValueTerm;

    #[rustfmt::skip]
    #[rstest]
    #[case::folds_field(
        SpinStateAst { unpaired: ValueAst::term(ValueTerm::Sum(vec![ValueTerm::Lit(1), ValueTerm::Lit(1)])), multiplicity: ValueAst::Lit(3) },
        Ok(SpinStateAst { unpaired: ValueAst::Lit(2), multiplicity: ValueAst::Lit(3) }),
    )]
    #[case::parity_invalid_is_allowed(
        SpinStateAst { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Lit(1) },
        Ok(SpinStateAst { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Lit(1) }),
    )]
    fn test_spin_state_ast_canonicalize(
        #[case] input: SpinStateAst,
        #[case] expected: Result<SpinStateAst, Contradiction>,
    ) {
        assert_eq!(input.canonicalize(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::both_undetermined(SpinStateAst::default())]
    #[case::valid_triplet(SpinStateAst { unpaired: ValueAst::Lit(2), multiplicity: ValueAst::Lit(3) })]
    #[case::parity_invalid(SpinStateAst { unpaired: ValueAst::Lit(1), multiplicity: ValueAst::Lit(1) })]
    fn test_spin_state_ast_canonicalize_identity(#[case] input: SpinStateAst) {
        assert_eq!(input.clone().canonicalize(), Ok(input));
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::default_undetermined(SpinStateAst::default(), false)]
    #[case::from_dsl(spin!("#u2").into(), true)]
    #[case::from_pair((2, 3).into(), true)]
    #[case::partial(SpinStateAst { unpaired: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined }, false)]
    fn test_spin_state_ast_is_ground(#[case] ast: SpinStateAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::both_undetermined(SpinStateAst::default(), true)]
    #[case::ground((2_u8, 3_u8).into(), false)]
    #[case::partial(SpinStateAst { unpaired: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined }, false)]
    fn test_spin_state_ast_is_undetermined(#[case] ast: SpinStateAst, #[case] expected: bool) {
        assert_eq!(ast.is_undetermined(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und_ground(SpinStateAst::default(), (2_u8, 3_u8).into(), Some((2_u8, 3_u8).into()))]
    #[case::exact((2_u8, 3_u8).into(), (2_u8, 3_u8).into(), Some((2_u8, 3_u8).into()))]
    #[case::unpaired_conflict((2_u8, 3_u8).into(), (0_u8, 1_u8).into(), None)]
    #[case::field_wise(
        SpinStateAst { unpaired: ValueAst::Undetermined, multiplicity: ValueAst::Lit(3) },
        SpinStateAst { unpaired: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined },
        Some((2_u8, 3_u8).into()),
    )]
    fn test_spin_state_ast_meet(
        #[case] a: SpinStateAst,
        #[case] b: SpinStateAst,
        #[case] expected: Option<SpinStateAst>,
    ) {
        assert_eq!(a.meet(&b), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::und((2_u8, 3_u8).into(), SpinStateAst::default(), SpinStateAst::default())]
    #[case::exact((2_u8, 3_u8).into(), (2_u8, 3_u8).into(), (2_u8, 3_u8).into())]
    #[case::field_wise_widen(
        (2_u8, 3_u8).into(),
        SpinStateAst { unpaired: ValueAst::Lit(2), multiplicity: ValueAst::Lit(1) },
        SpinStateAst { unpaired: ValueAst::Lit(2), multiplicity: ValueAst::lit_set([1, 3]) },
    )]
    fn test_spin_state_ast_join(
        #[case] a: SpinStateAst,
        #[case] b: SpinStateAst,
        #[case] expected: SpinStateAst,
    ) {
        assert_eq!(a.join(&b), expected);
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

    #[rstest]
    fn test_spin_state_ast_from_spin_state() {
        let ast: SpinStateAst = spin!("#u2").into();
        assert_eq!(ast.unpaired, ValueAst::Lit(2));
        assert_eq!(ast.multiplicity, ValueAst::Lit(3));
    }

    #[rstest]
    fn test_spin_state_ast_closed_shell() {
        let ast = SpinStateAst::closed_shell();
        assert_eq!(ast.unpaired, ValueAst::Lit(0));
        assert_eq!(ast.multiplicity, ValueAst::Lit(1));
    }

    #[rstest]
    fn test_spin_state_ast_from_str_roundtrip() {
        let s = SpinState::from_str("#u0#s1").unwrap();
        let ast: SpinStateAst = s.into();
        assert_eq!(ast, SpinStateAst::closed_shell());
    }
}
