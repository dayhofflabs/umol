//! Graph-based molecular models.

// Disabled during the umol-graph engines/config/data restructure (doc 92):
// - `ast` and `dsl` have moved into `umol-ast`; the in-tree copies are kept
//   for reference until the rewrite lands.
// - `api` is on hold pending the engine/config rebuild.
// - `io` parsers depend on `api::Molecule` and the old `ops::resolve::Resolver`
//   shape; re-enabled in phase 7 once the new resolver lands.
// pub mod ast;
// pub mod dsl;
// pub mod api;
// pub mod io;

pub mod diagnostics;
pub mod position;
pub mod ops;
pub mod span;
pub mod table_ir;
