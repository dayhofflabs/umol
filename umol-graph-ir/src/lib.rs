//! umol DSL (surface) and graph IR representations.

pub mod dsl;
pub mod ir;
pub mod macros;

/// The `mol!` visual-literal macro (desugars to `MoleculeSpec` build).
pub use umol_graph_ir_macros::{frag, mol};
