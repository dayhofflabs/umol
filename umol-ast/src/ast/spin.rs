//! Spin-state AST.

use umol_shared::error::SpinStateError;
use umol_shared::spin::{SpinMultiplicity, SpinState, MAX_UNPAIRED_ELECTRONS};

use super::value::ValueAst;

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

    /// Validate that literal fields fit within physical ranges and, when both
    /// fields are literal, that the `(unpaired, multiplicity)` pair is
    /// chemically compatible.
    pub fn validate(&self) -> Result<(), SpinStateError> {
        let u_lit =
            match &self.unpaired {
                ValueAst::Lit(n) => Some(u8::try_from(*n).map_err(|_| {
                    SpinStateError::UnpairedElectronsOutOfRange { unpaired: u8::MAX }
                })?),
                _ => None,
            };
        if let Some(u) = u_lit {
            if u > MAX_UNPAIRED_ELECTRONS {
                return Err(SpinStateError::UnpairedElectronsOutOfRange { unpaired: u });
            }
        }
        let m_mult = match &self.multiplicity {
            ValueAst::Lit(n) => {
                let m = u8::try_from(*n).map_err(|_| SpinStateError::MultiplicityOutOfRange {
                    multiplicity: u8::MAX,
                })?;
                Some(
                    SpinMultiplicity::from_multiplicity(m)
                        .ok_or(SpinStateError::MultiplicityOutOfRange { multiplicity: m })?,
                )
            }
            _ => None,
        };
        if let (Some(u), Some(m)) = (u_lit, m_mult) {
            if !SpinState::are_compatible(u, m) {
                return Err(SpinStateError::Incompatible {
                    unpaired: u,
                    multiplicity: m,
                });
            }
        }
        Ok(())
    }

    pub fn is_ground(&self) -> bool {
        matches!(&self.unpaired, ValueAst::Lit(_))
            && matches!(&self.multiplicity, ValueAst::Lit(_))
            && self.validate().is_ok()
    }

    /// Pattern matches target iff `unpaired` and `multiplicity` each
    /// match field-wise under `ValueAst::matches`.
    pub fn matches(&self, target: &Self) -> bool {
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

    #[rustfmt::skip]
    #[rstest]
    #[case::ok(SpinStateAst::new(2, 3), Ok(()))]
    #[case::unpaired_ok_partial(SpinStateAst::from_values(ValueAst::Lit(9), ValueAst::Undetermined), Ok(()))]
    #[case::unpaired_out_of_range(SpinStateAst::from_values(ValueAst::Lit(10), ValueAst::Undetermined),
        Err(SpinStateError::UnpairedElectronsOutOfRange { unpaired: 10 }))]
    #[case::multiplicity_out_of_range(SpinStateAst::from_values(ValueAst::Undetermined, ValueAst::Lit(11)),
        Err(SpinStateError::MultiplicityOutOfRange { multiplicity: 11 }))]
    #[case::incompatible(SpinStateAst::new(0, 3), Err(SpinStateError::Incompatible { unpaired: 0, multiplicity: SpinMultiplicity::Triplet }))]
    #[case::negative_unpaired(SpinStateAst::from_values(ValueAst::Lit(-1), ValueAst::Undetermined),
        Err(SpinStateError::UnpairedElectronsOutOfRange { unpaired: u8::MAX }))]
    fn test_spin_state_ast_validate(
        #[case] ast: SpinStateAst,
        #[case] expected: Result<(), SpinStateError>,
    ) {
        assert_eq!(ast.validate(), expected);
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
