//! Cross-check between local atom constraints and topology-derived values
//! across all entity types, plus molecule-scope constraint evaluation
//! (`:connected`, `:total-charge`, etc.).
//!
//! Stub: always returns `Determined`. Filled in once the per-relation
//! constraint evaluators land.

use thiserror::Error;
use umol_ast::ast::MoleculeAst;

use umol_shared::solution::Solution;

#[derive(Clone, Copy, Debug, Default)]
pub struct ConstraintValidator;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ConstraintContradiction {}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ConstraintError {}

impl ConstraintValidator {
    pub fn validate(
        &self,
        _ast: impl AsRef<MoleculeAst>,
    ) -> Result<Solution<(), ConstraintContradiction>, ConstraintError> {
        // TODO: stub. Per-relation constraint evaluators not yet implemented.
        // Aromatic systems: `ElectronCount(#e) == sum(electrons) - system.charge`.
        // Multicenter bonds: analogous rule.
        // Rings: sum of ring size counts == total ring count.
        // Molecule-scope constraints: `:connected`, `:total-charge`, etc.
        Ok(Solution::Determined(()))
    }
}
