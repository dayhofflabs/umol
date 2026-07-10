//! umol DSL (surface) and AST (semantic) representations.

pub mod ast;
pub mod dsl;
pub mod macros;

/// The `mol!` visual-literal macro (desugars to an L2 `MoleculeSpec` build).
pub use umol_ast_macros::mol;
