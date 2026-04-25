//! Graph-based molecular models.

// Disabled during the umol-ast wiring refactor:
// - `ast` and `dsl` have moved out into the `umol-ast` crate; the in-tree
//   copies are kept for reference until the rewrite lands.
// - `api` is on hold pending the data/engines/configs cleanup at the AST
//   layer.
// pub mod ast;
// pub mod dsl;
// pub mod api;

pub mod bond;
pub mod diagnostics;
pub mod io;
pub mod position;
pub mod ops;
pub mod span;
pub mod table_ir;
