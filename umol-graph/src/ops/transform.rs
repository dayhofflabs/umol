//! Pure transformations over a fully resolved `Molecule`. Distinct from
//! the resolvers in `ops/resolver`: a transformer rewrites a determined IR
//! into another determined IR without filling in undetermined values. Each
//! concrete transformer carries its own `Error` type via the trait's
//! associated type.

pub mod aromatizer;
pub mod delocalize_charge;
pub mod kekulizer;

pub use aromatizer::{Aromatizer, AromatizerError};
pub use delocalize_charge::DelocalizeCharge;
pub use kekulizer::{KekulizationConfig, Kekulizer, KekulizerError, MaximumMatchingAlgorithm};
use umol_graph_ir::ir::Molecule;

pub trait Transformer {
    type Error;

    fn transform_into(&self, molecule: &mut Molecule) -> Result<(), Self::Error>;

    fn transform(&self, molecule: &Molecule) -> Result<Molecule, Self::Error> {
        let mut out = molecule.clone();
        self.transform_into(&mut out)?;
        Ok(out)
    }

    /// Yields every result the transformer can produce. For deterministic
    /// transformers this is a single-element iterator; for non-deterministic
    /// ones this enumerates the alternatives. On error the iterator is
    /// empty.
    fn generate_all<'a>(
        &'a self,
        molecule: &'a Molecule,
    ) -> Box<dyn Iterator<Item = Molecule> + 'a>;
}
