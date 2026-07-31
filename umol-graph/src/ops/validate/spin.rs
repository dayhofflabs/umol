//! Spin-state invariant validation. A complete literal unpaired-electron count
//! and multiplicity must form a physically valid [`SpinState`]; a pair with
//! either component still non-literal is underdetermined.

use thiserror::Error;
use umol_ast::ast::{AsLit, AtomAst, MoleculeAst, UnpairedElectronsAst};
use umol_chem::error::SpinStateError;
use umol_chem::spin::SpinState;
use umol_utils::solution::Solution;

#[derive(Clone, Copy, Debug, Default)]
pub struct SpinInvariantsValidator;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum SpinInvariantsContradiction {
    #[error("atom has invalid unpaired electrons: {error}")]
    Atom { error: SpinStateError },
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum SpinInvariantsError {}

impl SpinInvariantsValidator {
    pub fn validate(
        &self,
        _ast: impl AsRef<MoleculeAst>,
    ) -> Result<Solution<(), SpinInvariantsContradiction>, SpinInvariantsError> {
        Ok(Solution::Determined(()))
    }

    pub fn validate_atom(
        &self,
        atom: &AtomAst,
    ) -> Result<Solution<(), SpinInvariantsContradiction>, SpinInvariantsError> {
        Ok(validate_unpaired_electrons(&atom.unpaired_electrons)
            .map_contradiction(|error| SpinInvariantsContradiction::Atom { error }))
    }
}

fn validate_unpaired_electrons(
    unpaired_electrons: &UnpairedElectronsAst,
) -> Solution<(), SpinStateError> {
    let Some(unpaired_electrons) = unpaired_electrons.as_lit() else {
        return Solution::Underdetermined(());
    };
    match SpinState::try_from(unpaired_electrons) {
        Ok(_) => Solution::Determined(()),
        Err(error) => Solution::Contradictory(error),
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;
    use rstest::rstest;
    use umol_ast::ast::ValueAst;
    use umol_chem::spin::SpinMultiplicity;

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::closed_shell((0_u8, 1_u8).into(), Solution::Determined(()))]
    #[case::doublet((1_u8, 2_u8).into(), Solution::Determined(()))]
    #[case::open_shell_singlet((2_u8, 1_u8).into(), Solution::Determined(()))]
    #[case::triplet((2_u8, 3_u8).into(), Solution::Determined(()))]
    #[case::count_undetermined(
        UnpairedElectronsAst { count: ValueAst::Undetermined, multiplicity: ValueAst::Lit(1) },
        Solution::Underdetermined(()),
    )]
    #[case::multiplicity_undetermined(
        UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Undetermined },
        Solution::Underdetermined(()),
    )]
    #[case::count_negative(
        UnpairedElectronsAst { count: ValueAst::Lit(-1), multiplicity: ValueAst::Lit(1) },
        Solution::Contradictory(SpinInvariantsContradiction::Atom {
            error: SpinStateError::UnpairedElectronsOutOfRange { count: -1 },
        }),
    )]
    #[case::count_above_u8(
        UnpairedElectronsAst { count: ValueAst::Lit(256), multiplicity: ValueAst::Lit(1) },
        Solution::Contradictory(SpinInvariantsContradiction::Atom {
            error: SpinStateError::UnpairedElectronsOutOfRange { count: 256 },
        }),
    )]
    #[case::multiplicity_zero(
        UnpairedElectronsAst { count: ValueAst::Lit(0), multiplicity: ValueAst::Lit(0) },
        Solution::Contradictory(SpinInvariantsContradiction::Atom {
            error: SpinStateError::MultiplicityOutOfRange { multiplicity: 0 },
        }),
    )]
    #[case::multiplicity_above_u8(
        UnpairedElectronsAst { count: ValueAst::Lit(0), multiplicity: ValueAst::Lit(256) },
        Solution::Contradictory(SpinInvariantsContradiction::Atom {
            error: SpinStateError::MultiplicityOutOfRange { multiplicity: 256 },
        }),
    )]
    #[case::parity(
        UnpairedElectronsAst { count: ValueAst::Lit(2), multiplicity: ValueAst::Lit(2) },
        Solution::Contradictory(SpinInvariantsContradiction::Atom {
            error: SpinStateError::Incompatible {
                unpaired_electrons: 2,
                multiplicity: SpinMultiplicity::DOUBLET,
            },
        }),
    )]
    #[case::above_maximum(
        UnpairedElectronsAst { count: ValueAst::Lit(0), multiplicity: ValueAst::Lit(2) },
        Solution::Contradictory(SpinInvariantsContradiction::Atom {
            error: SpinStateError::Incompatible {
                unpaired_electrons: 0,
                multiplicity: SpinMultiplicity::DOUBLET,
            },
        }),
    )]
    fn test_spin_invariants_validator_validate_atom(
        #[case] unpaired_electrons: UnpairedElectronsAst,
        #[case] expected: Solution<(), SpinInvariantsContradiction>,
    ) {
        assert_eq!(
            SpinInvariantsValidator
                .validate_atom(&AtomAst {
                    unpaired_electrons,
                    ..Default::default()
                })
                .unwrap(),
            expected,
        );
    }
}
