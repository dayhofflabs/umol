//! Binding exceptions. Minimal for now: a single catchable `ParseError` for DSL
//! parse failures. The full three-tier hierarchy mirroring umol-ast's error tiers
//! (doc 137 "Error mapping") lands with the Rust-side error sweep.

use pyo3::create_exception;
use pyo3::exceptions::PyException;
use pyo3::PyErr;
use umol_ast::dsl::ParseError as AstParseError;

create_exception!(
    umol,
    ParseError,
    PyException,
    "Raised when a umol DSL string fails to parse."
);

/// Map an `umol_ast` parse error onto the catchable `umol.ParseError`.
pub(crate) fn parse_error(error: AstParseError) -> PyErr {
    ParseError::new_err(error.to_string())
}
