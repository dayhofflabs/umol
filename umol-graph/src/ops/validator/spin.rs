//! Per-entity spin-coupling parity check: a literal `(unpaired, multiplicity)`
//! pair must satisfy `multiplicity = unpaired - 2k + 1` for some `k ∈
//! 0..=unpaired/2`. Runs on any entity carrying a `SpinStateAst` (atom,
//! aromatic system, multicenter bond).
//!
//! Stub: always returns `Determined`. Implementation pending; the parity
//! rule is in `umol_shared::spin::SpinState::are_compatible`.

use thiserror::Error;
use umol_ast::ast::{AtomAst, MoleculeAst};

use umol_shared::solution::Solution;

#[derive(Clone, Copy, Debug, Default)]
pub struct SpinCouplingValidator;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum SpinCouplingContradiction {}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum SpinCouplingError {}

impl SpinCouplingValidator {
    pub fn validate(
        &self,
        _ast: impl AsRef<MoleculeAst>,
    ) -> Result<Solution<(), SpinCouplingContradiction>, SpinCouplingError> {
        Ok(Solution::Determined(()))
    }

    pub fn validate_atom(
        &self,
        _atom: &AtomAst,
    ) -> Result<Solution<(), SpinCouplingContradiction>, SpinCouplingError> {
        Ok(Solution::Determined(()))
    }
}
