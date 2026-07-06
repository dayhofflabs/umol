//! Python bindings for umol: the `umol._native` extension module backing the
//! `umol` Python package. Types mirror the Rust API (see `umol-ast`).

use pyo3::prelude::*;

/// The native extension module. Wrapper types are registered here as the
/// binding grows; empty at scaffold stage.
#[pymodule]
fn _native(_module: &Bound<'_, PyModule>) -> PyResult<()> {
    Ok(())
}
