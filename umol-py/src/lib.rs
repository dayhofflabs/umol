//! Python bindings for umol: the `umol._native` extension module backing the
//! `umol` Python package. Types mirror the Rust API (see `umol-ast`).

use pyo3::prelude::*;

#[cfg(feature = "graph")]
use crate::{
    atom::{AtomAst, AtomId, AtomView, AtomViews, ElementAst, IsotopeMassAst, SpinStateAst},
    constraint::{
        AromaticValenceAst, AtomConstraintAst, AtomConstraintKey, AtomConstraintsAst,
        MulticenterValenceAst, RingMembershipAst, RingScope,
    },
    element::Element,
    molecule::MoleculeAst,
    stereo::{Permutation, StereoCosetAst, StereoTerm, TetrahedralStereoAst},
    value::{MemOp, RelOp, ValueAst, ValuePredicate, ValueTerm},
};

#[cfg(feature = "graph")]
mod atom;
#[cfg(feature = "graph")]
mod constraint;
#[cfg(feature = "graph")]
mod convert;
#[cfg(feature = "graph")]
mod element;
#[cfg(feature = "graph")]
mod molecule;
#[cfg(feature = "graph")]
mod stereo;
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
        module.add_class::<ValuePredicate>()?;
        module.add_class::<ValueAst>()?;
        module.add_class::<ElementAst>()?;
        module.add_class::<IsotopeMassAst>()?;
        module.add_class::<SpinStateAst>()?;
        module.add_class::<AtomAst>()?;
        module.add_class::<AtomId>()?;
        module.add_class::<AtomView>()?;
        module.add_class::<AtomViews>()?;
        module.add_class::<AromaticValenceAst>()?;
        module.add_class::<MulticenterValenceAst>()?;
        module.add_class::<RingScope>()?;
        module.add_class::<RingMembershipAst>()?;
        module.add_class::<AtomConstraintAst>()?;
        module.add_class::<AtomConstraintKey>()?;
        module.add_class::<AtomConstraintsAst>()?;
        module.add_class::<Permutation>()?;
        module.add_class::<StereoTerm>()?;
        module.add_class::<StereoCosetAst>()?;
        module.add_class::<TetrahedralStereoAst>()?;
    }
    Ok(())
}
