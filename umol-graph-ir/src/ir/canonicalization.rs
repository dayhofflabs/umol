//! Aggregate canonicalization inputs and failures.

use thiserror::Error;
use umol_graph_core::AutomorphismAlgorithm;

use super::error::Contradiction;
use super::molecule::MoleculeIntegrityError;
use super::reaction_span::ReactionSpanIntegrityError;
use super::validate::ReactionIntegrityError;

/// Semantic and operational inputs to aggregate canonicalization.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CanonicalizationContext {
    /// Whether stereo-sensitive refinement is iterated to a para-stereo fixpoint.
    pub para_stereo: bool,
    /// Graph automorphism algorithm used during canonical-frame search.
    pub automorphism_algorithm: AutomorphismAlgorithm,
}

/// Structural level used to select or compare an aggregate's entity frame.
///
/// The variants form a nested hierarchy. [`Topology`](Self::Topology) contains atoms and localized
/// bonds. [`Constitution`](Self::Constitution) additionally contains dative bonds, aromatic
/// systems, multicenter bonds, and noncovalent bonds. [`Full`](Self::Full) additionally contains
/// stereo atoms and stereo bonds. Constraints do not contribute to any structural level.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CanonicalizationLevel {
    Topology,
    Constitution,
    Full,
}

/// Failure to construct a canonical [`Molecule`](super::Molecule).
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum MoleculeCanonicalizationError {
    /// The molecule does not satisfy its representation-integrity contract.
    #[error(transparent)]
    Integrity(#[from] MoleculeIntegrityError),
    /// Intrinsic normalization of a carried value reached a contradiction.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

/// Failure to construct a canonical [`ReactionSpan`](super::ReactionSpan).
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionSpanCanonicalizationError {
    /// The reaction span does not satisfy its representation-integrity contract.
    #[error(transparent)]
    Integrity(#[from] ReactionSpanIntegrityError),
    /// Intrinsic normalization of a carried value reached a contradiction.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

/// Failure to construct a canonical [`Reaction`](super::Reaction).
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionCanonicalizationError {
    /// The reaction does not satisfy its representation-integrity contract.
    #[error(transparent)]
    Integrity(#[from] ReactionIntegrityError),
    /// Intrinsic normalization or span materialization reached a contradiction.
    #[error(transparent)]
    Contradiction(#[from] Contradiction),
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;
    use crate::ir::{AtomId, Entity};

    #[rstest]
    #[case::integrity(
        MoleculeCanonicalizationError::from(MoleculeIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        }),
        MoleculeCanonicalizationError::Integrity(MoleculeIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        }),
    )]
    #[case::contradiction(
        MoleculeCanonicalizationError::from(Contradiction),
        MoleculeCanonicalizationError::Contradiction(Contradiction)
    )]
    fn test_molecule_canonicalization_error_from(
        #[case] actual: MoleculeCanonicalizationError,
        #[case] expected: MoleculeCanonicalizationError,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::integrity(
        ReactionSpanCanonicalizationError::from(ReactionSpanIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        }),
        ReactionSpanCanonicalizationError::Integrity(
            ReactionSpanIntegrityError::InvalidReference {
                entity: Entity::Atom(AtomId(1)),
            },
        ),
    )]
    #[case::contradiction(
        ReactionSpanCanonicalizationError::from(Contradiction),
        ReactionSpanCanonicalizationError::Contradiction(Contradiction)
    )]
    fn test_reaction_span_canonicalization_error_from(
        #[case] actual: ReactionSpanCanonicalizationError,
        #[case] expected: ReactionSpanCanonicalizationError,
    ) {
        assert_eq!(actual, expected);
    }

    #[rstest]
    #[case::integrity(
        ReactionCanonicalizationError::from(ReactionIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        }),
        ReactionCanonicalizationError::Integrity(ReactionIntegrityError::InvalidReference {
            entity: Entity::Atom(AtomId(1)),
        }),
    )]
    #[case::contradiction(
        ReactionCanonicalizationError::from(Contradiction),
        ReactionCanonicalizationError::Contradiction(Contradiction)
    )]
    fn test_reaction_canonicalization_error_from(
        #[case] actual: ReactionCanonicalizationError,
        #[case] expected: ReactionCanonicalizationError,
    ) {
        assert_eq!(actual, expected);
    }
}
