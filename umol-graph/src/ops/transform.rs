//! Transformations of determined molecular representations.

pub mod aromatizer;
pub mod delocalize_charge;
pub mod kekulizer;

pub use aromatizer::{AromatizeError, Aromatizer};
pub use delocalize_charge::DelocalizeCharge;
pub use kekulizer::{KekulizeConfig, KekulizeError, Kekulizer, MaximumMatchingAlgorithm};
use umol_graph_ir::ir::Molecule;

/// Rewrites a determined molecular representation.
///
/// # Semantic properties
///
/// - Successful `transform` and `transform_into` produce the same molecule.
/// - Failed `transform_into` restores the input under `Molecule::normalized_eq`,
///   preserving participant order.
/// - Every published result satisfies molecule representation integrity.
/// - `transform_iter` preserves its source and yields independently mutable molecules.
pub trait Transformer {
    /// Failure reported by this transformation.
    type Error;

    /// Consumes the input and returns its transformed value.
    ///
    /// # Errors
    ///
    /// Returns the transformation's error on rejection. The input is dropped.
    fn transform(&self, molecule: Molecule) -> Result<Molecule, Self::Error>;

    /// Transforms the borrowed molecule in place, restoring it on failure.
    ///
    /// # Errors
    ///
    /// Returns the same errors as `transform`; the input is restored under
    /// `Molecule::normalized_eq` before returning.
    fn transform_into(&self, molecule: &mut Molecule) -> Result<(), Self::Error>;

    /// Lazily yields the transformation's independent results.
    ///
    /// Candidate copying and transformation begin when iteration requests a result.
    /// Deterministic transformations yield one molecule on success and none on error;
    /// transformations with alternatives enumerate their results in implementation order.
    fn transform_iter<'a>(&'a self, molecule: &'a Molecule) -> impl Iterator<Item = Molecule> + 'a;
}
