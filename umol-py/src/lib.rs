//! Python bindings for umol: the `umol._native` extension module backing the
//! `umol` Python package. Types mirror the Rust API (see `umol-ast`).

use pyo3::prelude::*;

#[cfg(feature = "graph")]
use crate::{
    element::Element,
    molecule::MoleculeAst,
    value::{MemOp, RelOp, ValueTerm},
};

#[cfg(feature = "graph")]
mod element;
#[cfg(feature = "graph")]
mod molecule;
#[cfg(feature = "graph")]
mod value;

/// The native extension module. Wrapper types are registered here as the
/// binding grows; the graph domain is gated behind the `graph` feature.
#[pymodule]
#[cfg_attr(not(feature = "graph"), allow(unused_variables))]
fn _native(module: &Bound<'_, PyModule>) -> PyResult<()> {
    #[cfg(feature = "graph")]
    {
        module.add_class::<Element>()?;
        module.add_class::<MoleculeAst>()?;
        module.add_class::<RelOp>()?;
        module.add_class::<MemOp>()?;
        module.add_class::<ValueTerm>()?;
    }
    Ok(())
}
