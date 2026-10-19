//! Transformations of determined molecular representations.

pub mod aromaticity;
pub mod hydrogens;

pub use aromaticity::{
    AromatizeError, Aromatizer, DelocalizeCharge, KekulizeConfig, KekulizeError, Kekulizer,
    MaximumMatchingAlgorithm,
};
pub use hydrogens::{FoldHydrogensError, HydrogenFolder, HydrogenUnfolder, UnfoldHydrogensError};
use umol_graph_ir::ir::Molecule;

/// Rewrites a determined molecular representation.
///
/// Derives a complete plan before mutation, reporting every rejection during planning.
/// The plan applies to its unchanged input and preserves representation integrity;
/// a failed application or publication is an implementation defect.
///
/// # Semantic properties
///
/// - Successful `transform` and `transform_into` produce the same molecule.
/// - Rejected `transform_into` leaves the input unchanged.
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

    /// Plans and transforms the borrowed molecule in place without a transaction or rollback.
    ///
    /// # Errors
    ///
    /// Returns the same planning errors as `transform`; rejection leaves the input unchanged.
    fn transform_into(&self, molecule: &mut Molecule) -> Result<(), Self::Error>;

    /// Lazily yields the transformation's independent results.
    ///
    /// Candidate copying and transformation begin when iteration requests a result.
    /// Deterministic transformations yield one molecule on success and none on error;
    /// transformations with alternatives enumerate their results in implementation order.
    fn transform_iter<'a>(&'a self, molecule: &'a Molecule) -> impl Iterator<Item = Molecule> + 'a;
}
