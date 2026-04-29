//! Pure transformations over a fully resolved `MoleculeAst`. Distinct from
//! the resolvers in `ops/`: a transformer rewrites a determined AST into
//! another determined AST without filling in undetermined values. Two
//! concrete implementations land in this module: [`Aromatize`] (Kekulé form
//! → aromatic-system form, via the same perception as resolution) and
//! [`Kekulize`] (aromatic-system form → Kekulé bond orders, via perfect
//! matching on the aromatic subgraph).

pub mod aromatize;
pub mod kekulize;

use thiserror::Error;
use umol_ast::ast::{AromaticSystemIdx, MoleculeAst};

use crate::ops::aromaticity::{AromaticityContradiction, AromaticityError};

pub use aromatize::Aromatize;
pub use kekulize::{Kekulize, KekulizationModel};

/// Flat error set covering every concrete `Transformer`. Keeps the trait
/// dyn-compatible (`Vec<Box<dyn Transformer>>` works) without forcing
/// callers to box the error.
#[derive(Debug, Error, Clone, PartialEq)]
pub enum TransformerError {
    #[error("aromatize setup error: {0}")]
    AromatizeSetup(#[from] AromaticityError),
    #[error("aromatize chemistry contradiction: {0}")]
    AromatizeContradiction(AromaticityContradiction),
    #[error("aromatize: input is underdetermined; resolve first")]
    AromatizeUnderdetermined,
    #[error("kekulize: no perfect matching exists for aromatic system {0:?}")]
    KekulizeNoMatching(AromaticSystemIdx),
}

pub trait Transformer {
    fn transform_into(&self, ast: &mut MoleculeAst) -> Result<(), TransformerError>;

    fn transform(&self, ast: &MoleculeAst) -> Result<MoleculeAst, TransformerError> {
        let mut out = ast.clone();
        self.transform_into(&mut out)?;
        Ok(out)
    }

    /// Yields every result the transformer can produce. For deterministic
    /// transformers (Aromatize) this is a single-element iterator; for
    /// non-deterministic ones (Kekulize when enumeration is wired up) this
    /// enumerates the alternatives. On error the iterator is empty.
    fn generate_all<'a>(
        &'a self,
        ast: &'a MoleculeAst,
    ) -> Box<dyn Iterator<Item = MoleculeAst> + 'a>;
}
