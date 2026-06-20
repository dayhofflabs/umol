//! Pure transformations over a fully resolved `MoleculeAst`. Distinct from
//! the resolvers in `ops/resolver`: a transformer rewrites a determined AST
//! into another determined AST without filling in undetermined values. Each
//! concrete transformer carries its own `Error` type via the trait's
//! associated type.

pub mod aromatizer;
pub mod kekulizer;

pub use aromatizer::{Aromatizer, AromatizerError};
pub use kekulizer::{KekulizationModel, Kekulizer, KekulizerError};
use umol_ast::ast::MoleculeAst;

pub trait Transformer {
    type Error;

    fn transform_into(&self, ast: &mut MoleculeAst) -> Result<(), Self::Error>;

    fn transform(&self, ast: &MoleculeAst) -> Result<MoleculeAst, Self::Error> {
        let mut out = ast.clone();
        self.transform_into(&mut out)?;
        Ok(out)
    }

    /// Yields every result the transformer can produce. For deterministic
    /// transformers this is a single-element iterator; for non-deterministic
    /// ones this enumerates the alternatives. On error the iterator is
    /// empty.
    fn generate_all<'a>(
        &'a self,
        ast: &'a MoleculeAst,
    ) -> Box<dyn Iterator<Item = MoleculeAst> + 'a>;
}
