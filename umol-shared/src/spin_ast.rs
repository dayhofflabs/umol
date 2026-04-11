//! Spin-state AST: ground or pattern over (unpaired_electrons, multiplicity).

use crate::error::SpinStateError;
use crate::spin::{SpinMultiplicity, SpinState};
use crate::value_ast::ValueAst;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum SpinStateAst {
    Wildcard,
    Lit(SpinState),
    Pair {
        unpaired: Option<ValueAst>,
        multiplicity: Option<ValueAst>,
    },
}

impl SpinStateAst {
    pub fn is_ground(&self) -> bool {
        matches!(self, Self::Lit(_))
    }

    /// Build an `Option<SpinStateAst>` from an `(unpaired, multiplicity)` `ValueAst` pair.
    /// Returns `None` when both inputs are `None`.
    pub fn from_pair(unpaired: Option<ValueAst>, multiplicity: Option<ValueAst>) -> Option<Self> {
        match (unpaired, multiplicity) {
            (None, None) => None,
            (u, m) => Some(Self::Pair {
                unpaired: u,
                multiplicity: m,
            }),
        }
    }

    /// Decompose into `(unpaired, multiplicity)` `ValueAst` pair, the form
    /// expected by the existing spin lowering helpers.
    pub fn to_pair(&self) -> (Option<ValueAst>, Option<ValueAst>) {
        match self {
            Self::Wildcard => (Some(ValueAst::Wildcard), Some(ValueAst::Wildcard)),
            Self::Lit(s) => (
                Some(ValueAst::Lit(s.unpaired_electrons() as i64)),
                Some(ValueAst::Lit(s.multiplicity().multiplicity() as i64)),
            ),
            Self::Pair {
                unpaired,
                multiplicity,
            } => (unpaired.clone(), multiplicity.clone()),
        }
    }

    /// Test whether a concrete `SpinState` satisfies this AST.
    pub fn matches(&self, state: SpinState) -> bool {
        match self {
            Self::Wildcard => true,
            Self::Lit(s) => *s == state,
            Self::Pair {
                unpaired,
                multiplicity,
            } => {
                let u_ok = match unpaired {
                    None => true,
                    Some(v) => v.matches(state.unpaired_electrons() as i64),
                };
                let m_ok = match multiplicity {
                    None => true,
                    Some(v) => v.matches(state.multiplicity().multiplicity() as i64),
                };
                u_ok && m_ok
            }
        }
    }

    /// Collapse to a `SpinState` when both unpaired and multiplicity are concrete literals.
    ///
    /// Returns `Ok(None)` when the AST is not fully grounded; `Err` if the literal pair
    /// is out of range or fails [`SpinState::are_compatible`].
    pub fn try_into_ground(&self) -> Result<Option<SpinState>, SpinStateError> {
        match self {
            Self::Wildcard => Ok(None),
            Self::Lit(s) => Ok(Some(*s)),
            Self::Pair {
                unpaired,
                multiplicity,
            } => match (unpaired, multiplicity) {
                (Some(ValueAst::Lit(u)), Some(ValueAst::Lit(m))) => {
                    let u_u8 = u8::try_from(*u).map_err(|_| {
                        SpinStateError::UnpairedElectronsOutOfRange {
                            unpaired_electrons: u8::MAX,
                        }
                    })?;
                    let m_u8 = u8::try_from(*m).map_err(|_| {
                        SpinStateError::MultiplicityOutOfRange {
                            multiplicity: u8::MAX,
                        }
                    })?;
                    let mult = SpinMultiplicity::from_multiplicity(m_u8)
                        .ok_or(SpinStateError::MultiplicityOutOfRange { multiplicity: m_u8 })?;
                    SpinState::try_new(u_u8, mult).map(Some)
                }
                _ => Ok(None),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use pretty_assertions::assert_eq;
    use rstest::*;

    use super::*;
    use crate::spin;

    #[rstest]
    #[case(SpinStateAst::Wildcard, false)]
    #[case(SpinStateAst::Lit(spin!("#u2")), true)]
    #[case(SpinStateAst::Pair { unpaired: Some(ValueAst::Lit(2)), multiplicity: Some(ValueAst::Lit(3)) }, false)]
    fn test_spin_state_ast_is_ground(#[case] ast: SpinStateAst, #[case] expected: bool) {
        assert_eq!(ast.is_ground(), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard_any(SpinStateAst::Wildcard, spin!("#u2"), true)]
    #[case::lit_match(SpinStateAst::Lit(spin!("#u2")), spin!("#u2"), true)]
    #[case::lit_mismatch(SpinStateAst::Lit(spin!("#u2")), spin!("#u0"), false)]
    #[case::pair_both_match(SpinStateAst::Pair { unpaired: Some(ValueAst::Lit(2)), multiplicity: Some(ValueAst::Lit(3)) }, spin!("#u2"), true)]
    #[case::pair_unpaired_only(SpinStateAst::Pair { unpaired: Some(ValueAst::Lit(2)), multiplicity: None }, spin!("#u2#s1"), true)]
    #[case::pair_unpaired_mismatch(SpinStateAst::Pair { unpaired: Some(ValueAst::Lit(2)), multiplicity: None }, spin!("#u0"), false)]
    #[case::pair_multiplicity_only(SpinStateAst::Pair { unpaired: None, multiplicity: Some(ValueAst::Lit(3)) }, spin!("#u2"), true)]
    #[case::pair_both_wildcard(SpinStateAst::Pair { unpaired: Some(ValueAst::Wildcard), multiplicity: Some(ValueAst::Wildcard) }, spin!("#u2"), true)]
    fn test_spin_state_ast_matches(#[case] ast: SpinStateAst, #[case] state: SpinState, #[case] expected: bool) {
        assert_eq!(ast.matches(state), expected);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::wildcard(SpinStateAst::Wildcard, Ok(None))]
    #[case::lit(SpinStateAst::Lit(spin!("#u2")), Ok(Some(spin!("#u2"))))]
    #[case::pair_concrete(SpinStateAst::Pair { unpaired: Some(ValueAst::Lit(2)), multiplicity: Some(ValueAst::Lit(3)) }, Ok(Some(spin!("#u2#s3"))))]
    #[case::pair_partial(SpinStateAst::Pair { unpaired: Some(ValueAst::Lit(2)), multiplicity: None }, Ok(None))]
    #[case::pair_wildcard(SpinStateAst::Pair { unpaired: Some(ValueAst::Wildcard), multiplicity: Some(ValueAst::Wildcard) }, Ok(None))]
    #[case::pair_incompatible(SpinStateAst::Pair { unpaired: Some(ValueAst::Lit(0)), multiplicity: Some(ValueAst::Lit(3)) }, Err(SpinStateError::Incompatible { unpaired_electrons: 0, multiplicity: SpinMultiplicity::Triplet }))]
    fn test_spin_state_ast_try_into_ground(#[case] ast: SpinStateAst, #[case] expected: Result<Option<SpinState>, SpinStateError>) {
        assert_eq!(ast.try_into_ground(), expected);
    }
}
