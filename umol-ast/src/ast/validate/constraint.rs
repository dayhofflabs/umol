//! Tier-1 constraint validator: cross-entity and molecule-scope constraint evaluation. Run at AST
//! construction/raise and available standalone; never consults a chemistry model.

use thiserror::Error;
use umol_utils::solution::Solution;

use super::super::molecule::MoleculeAst;

/// Cross-check between local atom constraints and topology-derived values across
/// all entity types, plus molecule-scope constraint evaluation (`:connected`,
/// `:total-charge`, etc.).
///
/// Stub: always returns `Determined`. Filled in once the per-relation constraint
/// evaluators land.
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
        // Entity constraint consistency: an entity's inline (entity-local) constraints
        // and the molecule-scope (`:constraints`) entries referencing it must be jointly
        // satisfiable — a same-kind conflict (e.g. inline `#v4` vs `{:atom [i {:valence 3}]}`)
        // is a contradiction.
        // Aromatic systems: `ElectronCount(#e) == sum(electrons) - system.charge`.
        // Multicenter bonds: analogous rule.
        // Rings: sum of ring size counts == total ring count.
        // Molecule-scope constraints: `:connected`, `:total-charge`, etc.
        Ok(Solution::Determined(()))
    }
}
