//! Transformations of determined molecular representations.

pub mod aromaticity;
pub mod hydrogens;
pub mod stereo;

use std::convert::Infallible;

pub use aromaticity::{
    AromatizeError, Aromatizer, ChargeDelocalizer, KekulizeConfig, KekulizeError, Kekulizer,
    MaximumMatchingAlgorithm,
};
pub use hydrogens::{FoldHydrogensError, HydrogenFolder, HydrogenUnfolder, UnfoldHydrogensError};
pub use stereo::AddStereoConfig;
use umol_graph_ir::ir::{AtomId, Molecule};

use crate::ops::aromaticity::AromaticityConfig;
use crate::ops::model::AromaticityModel;

/// Rewrites a determined molecular representation.
///
/// Implemented by transformation engines. [`Transform`] provides the named methods on Molecule.
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

/// Named transformations on a molecular representation.
///
/// Consuming methods return the transformed value and drop the input on rejection.
/// Their `_into` counterparts mutate the receiver and leave it unchanged on rejection.
/// Models, configurations, and atom orders are supplied explicitly.
///
/// Each call constructs the corresponding [`Transformer`] and publishes independently.
/// In a sequence of in-place calls, a later rejection preserves earlier successful changes.
///
/// # Examples
///
/// ```
/// use umol_graph::ops::transform::Transform;
/// use umol_graph_ir::mol_dsl_concrete;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let explicit = mol_dsl_concrete!(r#"{:atoms ["C#h3" "H"] :bonds [[0 1 "1"]]}"#);
/// let mut molecule = explicit.fold_hydrogens()?;
/// assert_eq!(molecule, mol_dsl_concrete!(r#"{:atoms ["C#h4"]}"#));
///
/// molecule.unfold_hydrogens_into()?;
/// assert_eq!(molecule.atoms().count(), 5);
/// # Ok(())
/// # }
/// ```
///
/// # Semantic properties
///
/// - Consuming and in-place methods agree under structural equality for equal inputs and arguments.
/// - Each successful result satisfies molecule representation integrity.
/// - Each rejected in-place call preserves its entry state under structural equality.
pub trait Transform: Sized {
    /// Transfers localized π electrons into perceived aromatic systems using [`Aromatizer`].
    ///
    /// # Errors
    ///
    /// Reports missing model parameters, underdetermined input, or a contradiction in the
    /// proposed aromatic systems. Rejection drops the input.
    fn aromatize(
        self,
        model: &AromaticityModel,
        config: AromaticityConfig,
    ) -> Result<Self, AromatizeError>;

    /// Aromatizes the receiver in place using [`Aromatizer`].
    ///
    /// # Errors
    ///
    /// Reports the same errors as [`Self::aromatize`], leaving the receiver unchanged.
    fn aromatize_into(
        &mut self,
        model: &AromaticityModel,
        config: AromaticityConfig,
    ) -> Result<(), AromatizeError>;

    /// Localizes aromatic systems using [`Kekulizer`] and the supplied atom order.
    ///
    /// # Errors
    ///
    /// Reports incomplete or unsupported system data, an invalid atom order, unavailable
    /// matching, or a contradiction in the proposed localization. Rejection drops the input.
    fn kekulize(
        self,
        node_order: Vec<AtomId>,
        config: KekulizeConfig,
    ) -> Result<Self, KekulizeError>;

    /// Kekulizes the receiver in place using [`Kekulizer`] and the supplied atom order.
    ///
    /// # Errors
    ///
    /// Reports the same errors as [`Self::kekulize`], leaving the receiver unchanged.
    fn kekulize_into(
        &mut self,
        node_order: Vec<AtomId>,
        config: KekulizeConfig,
    ) -> Result<(), KekulizeError>;

    /// Delocalizes charge over homogeneous aromatic systems using [`ChargeDelocalizer`].
    ///
    /// This operation is infallible; systems without an applicable plan retain their state.
    fn delocalize_charge(self) -> Result<Self, Infallible>;

    /// Delocalizes charge in place using [`ChargeDelocalizer`].
    ///
    /// This operation is infallible; systems without an applicable plan retain their state.
    fn delocalize_charge_into(&mut self) -> Result<(), Infallible>;

    /// Folds eligible ordinary explicit hydrogens using [`HydrogenFolder`].
    ///
    /// # Errors
    ///
    /// Reports a negative affected implicit-hydrogen count or a stereo frame that cannot
    /// represent the folded hydrogens. Rejection drops the input.
    fn fold_hydrogens(self) -> Result<Self, FoldHydrogensError>;

    /// Folds eligible ordinary explicit hydrogens in place using [`HydrogenFolder`].
    ///
    /// # Errors
    ///
    /// Reports the same errors as [`Self::fold_hydrogens`], leaving the receiver unchanged.
    fn fold_hydrogens_into(&mut self) -> Result<(), FoldHydrogensError>;

    /// Expands literal implicit-hydrogen counts using [`HydrogenUnfolder`].
    ///
    /// # Errors
    ///
    /// Reports a negative count, a nonliteral count backing a virtual stereo hydrogen,
    /// or disagreement between a count and its stereo frame. Rejection drops the input.
    fn unfold_hydrogens(self) -> Result<Self, UnfoldHydrogensError>;

    /// Expands literal implicit-hydrogen counts in place using [`HydrogenUnfolder`].
    ///
    /// # Errors
    ///
    /// Reports the same errors as [`Self::unfold_hydrogens`], leaving the receiver unchanged.
    fn unfold_hydrogens_into(&mut self) -> Result<(), UnfoldHydrogensError>;
}

impl Transform for Molecule {
    fn aromatize(
        self,
        model: &AromaticityModel,
        config: AromaticityConfig,
    ) -> Result<Self, AromatizeError> {
        Aromatizer::with_config(model, config).transform(self)
    }

    fn aromatize_into(
        &mut self,
        model: &AromaticityModel,
        config: AromaticityConfig,
    ) -> Result<(), AromatizeError> {
        Aromatizer::with_config(model, config).transform_into(self)
    }

    fn kekulize(
        self,
        node_order: Vec<AtomId>,
        config: KekulizeConfig,
    ) -> Result<Self, KekulizeError> {
        Kekulizer::new(config, node_order).transform(self)
    }

    fn kekulize_into(
        &mut self,
        node_order: Vec<AtomId>,
        config: KekulizeConfig,
    ) -> Result<(), KekulizeError> {
        Kekulizer::new(config, node_order).transform_into(self)
    }

    fn delocalize_charge(self) -> Result<Self, Infallible> {
        ChargeDelocalizer.transform(self)
    }

    fn delocalize_charge_into(&mut self) -> Result<(), Infallible> {
        ChargeDelocalizer.transform_into(self)
    }

    fn fold_hydrogens(self) -> Result<Self, FoldHydrogensError> {
        HydrogenFolder.transform(self)
    }

    fn fold_hydrogens_into(&mut self) -> Result<(), FoldHydrogensError> {
        HydrogenFolder.transform_into(self)
    }

    fn unfold_hydrogens(self) -> Result<Self, UnfoldHydrogensError> {
        HydrogenUnfolder.transform(self)
    }

    fn unfold_hydrogens_into(&mut self) -> Result<(), UnfoldHydrogensError> {
        HydrogenUnfolder.transform_into(self)
    }
}

#[cfg(test)]
mod tests;
