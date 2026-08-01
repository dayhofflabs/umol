//! Tier-1 constraint validator: cross-entity and molecule-scope constraint evaluation. Run at AST
//! construction/raise and available standalone; never consults a chemistry model.

use thiserror::Error;
use umol_graph_core::{
    ConnectedComponentsAlgorithm, RelevantCycleEnumerationAlgorithm, SubgraphIsomorphismAlgorithm,
};
use umol_utils::solution::Solution;

pub mod incidence;
pub mod relational;
pub mod ring;

pub use incidence::{IncidenceConstraintContradiction, IncidenceConstraintValidator};
pub use relational::{RelationalConstraintContradiction, RelationalConstraintValidator};
pub use ring::{RingConstraintContradiction, RingConstraintValidator};

use super::super::entity::Entity;
use super::super::id::DativeBondId;
use super::super::molecule::MoleculeAst;
use super::super::substructure::SubstructureMatchAlgorithm;

/// Algorithm selectors used by complete model-independent constraint validation.
///
/// Focused validators take only the selectors they require. This bundle has no
/// default at the AST layer so every algorithm choice remains explicit.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ConstraintValidateConfig {
    pub relevant_cycle_algorithm: RelevantCycleEnumerationAlgorithm,
    pub connected_components_algorithm: ConnectedComponentsAlgorithm,
    pub substructure_match_algorithm: SubstructureMatchAlgorithm,
    pub subgraph_isomorphism_algorithm: SubgraphIsomorphismAlgorithm,
}

/// Cross-check between local atom constraints and topology-derived values across
/// all entity types, plus molecule-scope constraint evaluation (`:connected`,
/// `:total-charge`, etc.).
///
/// Stub: always returns `Determined`. Filled in once the per-relation constraint
/// evaluators land.
#[derive(Clone, Copy, Debug, Default)]
pub struct ConstraintValidator;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ConstraintContradiction {
    #[error(transparent)]
    Incidence(#[from] IncidenceConstraintContradiction),
    #[error(transparent)]
    Ring(#[from] RingConstraintContradiction),
    #[error(transparent)]
    Relational(#[from] RelationalConstraintContradiction),
}

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum ConstraintError {
    #[error("constraint references unavailable {entity}")]
    InvalidReference { entity: Entity },
    /// Dative-bond ring topology is deferred to the coordination/haptic entity split in doc 117.
    #[error("ring membership is not defined for dative bond {bond:?}")]
    DativeBondRingMembershipUnsupported { bond: DativeBondId },
}

impl ConstraintValidator {
    pub fn validate(
        &self,
        _ast: &MoleculeAst,
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

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::ast::constraint::{
        AtomConstraintAst, BondConstraintAst, RelationalConstraint, RingScope,
    };
    use crate::ast::id::{AtomId, BondId, DativeBondId};

    #[rstest]
    #[case::incidence(
        IncidenceConstraintContradiction::Atom {
            atom: AtomId(2),
            constraint: AtomConstraintAst::valence(4),
        },
        ConstraintContradiction::Incidence(IncidenceConstraintContradiction::Atom {
            atom: AtomId(2),
            constraint: AtomConstraintAst::valence(4),
        })
    )]
    #[case::ring(
        RingConstraintContradiction::Bond {
            bond: BondId(3),
            constraint: BondConstraintAst::ring_membership(RingScope::Size(6), 1),
        },
        ConstraintContradiction::Ring(RingConstraintContradiction::Bond {
            bond: BondId(3),
            constraint: BondConstraintAst::ring_membership(RingScope::Size(6), 1),
        })
    )]
    #[case::relational(
        RelationalConstraintContradiction {
            constraint: RelationalConstraint::DativeBondDonor {
                bond: DativeBondId(1),
                atom: AtomId(2),
            },
        },
        ConstraintContradiction::Relational(RelationalConstraintContradiction {
            constraint: RelationalConstraint::DativeBondDonor {
                bond: DativeBondId(1),
                atom: AtomId(2),
            },
        })
    )]
    fn test_constraint_contradiction_from<T>(
        #[case] input: T,
        #[case] expected: ConstraintContradiction,
    ) where
        T: Into<ConstraintContradiction>,
    {
        assert_eq!(input.into(), expected);
    }

    #[rstest]
    #[case::invalid_reference(
        ConstraintError::InvalidReference {
            entity: Entity::Atom(AtomId(4)),
        },
        "constraint references unavailable atom 4"
    )]
    #[case::dative_ring(
        ConstraintError::DativeBondRingMembershipUnsupported {
            bond: DativeBondId(5),
        },
        "ring membership is not defined for dative bond DativeBondId(5)"
    )]
    fn test_constraint_error_display(#[case] input: ConstraintError, #[case] expected: &str) {
        assert_eq!(input.to_string(), expected);
    }
}
