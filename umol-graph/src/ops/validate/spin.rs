//! Per-entity spin-coupling parity check: a literal `(unpaired, multiplicity)`
//! pair must satisfy `multiplicity = unpaired - 2k + 1` for some `k ∈
//! 0..=unpaired/2`. Runs on any entity carrying a `SpinStateAst` (atom,
//! aromatic system, multicenter bond).
//!
//! Stub: always returns `Determined`. Implementation pending; complete literal
//! pairs are validated by conversion to `umol_chem::spin::SpinState`.

use thiserror::Error;
use umol_ast::ast::{AtomAst, MoleculeAst};
use umol_utils::solution::Solution;

#[derive(Clone, Copy, Debug, Default)]
pub struct SpinInvariantsValidator;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum SpinInvariantsContradiction {}

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
        _atom: &AtomAst,
    ) -> Result<Solution<(), SpinInvariantsContradiction>, SpinInvariantsError> {
        Ok(Solution::Determined(()))
    }
}
