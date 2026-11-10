//! Stereo processing sites, selections, and inconsistency diagnostics.

use std::collections::BTreeSet;

use pyo3::basic::CompareOp;
use pyo3::prelude::*;
use umol_graph::ops::stereo::{
    StereoInconsistency as GraphStereoInconsistency, StereoSelection as GraphStereoSelection,
    StereoSite as GraphStereoSite,
};
use umol_graph_ir::ir::{
    AtomId as GraphIrAtomId, BondId as GraphIrBondId, StereoAtomId as GraphIrStereoAtomId,
    StereoBondId as GraphIrStereoBondId,
};

use super::{StereoKind, Stereogenicity};
use crate::convert::{hash_rust, FrozenSetItems};

/// An atom or bond topology ID, independent of a molecule and of stored stereo.
#[pyclass(frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoSite {
    Atom(u32),
    Bond(u32),
}

#[pymethods]
impl StereoSite {
    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    fn __richcmp__(&self, other: &Self, op: CompareOp) -> bool {
        op.matches(self.to_rust().cmp(&other.to_rust()))
    }

    fn __repr__(&self) -> String {
        match self {
            Self::Atom(id) => format!("StereoSite.Atom({id})"),
            Self::Bond(id) => format!("StereoSite.Bond({id})"),
        }
    }
}

impl StereoSite {
    pub(crate) fn from_rust(site: GraphStereoSite) -> Self {
        match site {
            GraphStereoSite::Atom(id) => Self::Atom(id.0),
            GraphStereoSite::Bond(id) => Self::Bond(id.0),
        }
    }

    pub(crate) fn to_rust(self) -> GraphStereoSite {
        match self {
            Self::Atom(id) => GraphStereoSite::Atom(GraphIrAtomId(id)),
            Self::Bond(id) => GraphStereoSite::Bond(GraphIrBondId(id)),
        }
    }
}

/// Intersecting site, kind, and stereogenicity filters.
///
/// None is unrestricted; an empty collection selects nothing. IDs are checked
/// against a molecule by the executing operation, not during construction.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StereoSelection(GraphStereoSelection);

#[pymethods]
impl StereoSelection {
    #[new]
    #[pyo3(signature = (*, sites=None, kinds=None, stereogenicities=None))]
    fn new(
        sites: Option<&Bound<'_, PyAny>>,
        kinds: Option<&Bound<'_, PyAny>>,
        stereogenicities: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        Ok(Self::from_rust(GraphStereoSelection {
            sites: sites
                .map(|values| {
                    values
                        .try_iter()?
                        .map(|value| Ok(value?.extract::<StereoSite>()?.to_rust()))
                        .collect::<PyResult<_>>()
                })
                .transpose()?,
            kinds: kinds
                .map(|values| {
                    values
                        .try_iter()?
                        .map(|value| Ok(value?.extract::<StereoKind>()?.to_rust()))
                        .collect::<PyResult<_>>()
                })
                .transpose()?,
            stereogenicities: stereogenicities
                .map(|values| {
                    values
                        .try_iter()?
                        .map(|value| Ok(value?.extract::<Stereogenicity>()?.to_rust()))
                        .collect::<PyResult<_>>()
                })
                .transpose()?,
        }))
    }

    #[getter]
    fn sites(&self) -> Option<FrozenSetItems<BTreeSet<StereoSite>>> {
        self.to_rust()
            .sites
            .as_ref()
            .map(|values| values.iter().copied().map(StereoSite::from_rust).collect())
    }

    #[getter]
    fn kinds(&self) -> Option<FrozenSetItems<BTreeSet<StereoKind>>> {
        self.to_rust()
            .kinds
            .as_ref()
            .map(|values| values.iter().copied().map(StereoKind::from_rust).collect())
    }

    #[getter]
    fn stereogenicities(&self) -> Option<FrozenSetItems<BTreeSet<Stereogenicity>>> {
        self.to_rust().stereogenicities.as_ref().map(|values| {
            values
                .iter()
                .copied()
                .map(Stereogenicity::from_rust)
                .collect()
        })
    }

    /// Test supplied values without looking up a molecule or calculating stereo.
    #[pyo3(signature = (site, kind=None, stereogenicity=None))]
    fn matches(
        &self,
        site: StereoSite,
        kind: Option<StereoKind>,
        stereogenicity: Option<Stereogenicity>,
    ) -> bool {
        self.to_rust().matches(
            site.to_rust(),
            kind.map(StereoKind::to_rust),
            stereogenicity.map(Stereogenicity::to_rust),
        )
    }

    pub(crate) fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "StereoSelection(sites={}, kinds={}, stereogenicities={})",
            self.sites()
                .into_pyobject(py)?
                .repr()?
                .extract::<String>()?,
            self.kinds()
                .into_pyobject(py)?
                .repr()?
                .extract::<String>()?,
            self.stereogenicities()
                .into_pyobject(py)?
                .repr()?
                .extract::<String>()?
        ))
    }
}

impl StereoSelection {
    pub(crate) fn from_rust(selection: GraphStereoSelection) -> Self {
        Self(selection)
    }
    pub(crate) fn to_rust(&self) -> &GraphStereoSelection {
        &self.0
    }
}

/// A stereo evidence or stored-entity inconsistency, with the Rust variant fields.
#[pyclass(frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum StereoInconsistency {
    TetrahedralStereoFailure { atom: u32 },
    StereoAtomFailure { stereo_atom: u32 },
    TetrahedralStereoMismatch { atom: u32, stereo_atom: u32 },
    CisTransStereoFailure { bond: u32 },
    StereoBondFailure { stereo_bond: u32 },
    CisTransStereoMismatch { bond: u32, stereo_bond: u32 },
}

#[pymethods]
impl StereoInconsistency {
    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    fn __richcmp__(&self, other: &Self, op: CompareOp) -> bool {
        op.matches(self.to_rust().cmp(&other.to_rust()))
    }

    fn __repr__(&self) -> String {
        match self {
            Self::TetrahedralStereoFailure { atom } => format!("StereoInconsistency.TetrahedralStereoFailure(atom={atom})"),
            Self::StereoAtomFailure { stereo_atom } => format!("StereoInconsistency.StereoAtomFailure(stereo_atom={stereo_atom})"),
            Self::TetrahedralStereoMismatch { atom, stereo_atom } => format!("StereoInconsistency.TetrahedralStereoMismatch(atom={atom}, stereo_atom={stereo_atom})"),
            Self::CisTransStereoFailure { bond } => format!("StereoInconsistency.CisTransStereoFailure(bond={bond})"),
            Self::StereoBondFailure { stereo_bond } => format!("StereoInconsistency.StereoBondFailure(stereo_bond={stereo_bond})"),
            Self::CisTransStereoMismatch { bond, stereo_bond } => format!("StereoInconsistency.CisTransStereoMismatch(bond={bond}, stereo_bond={stereo_bond})"),
        }
    }
    fn __str__(&self) -> String {
        self.to_rust().to_string()
    }
}

impl StereoInconsistency {
    pub(crate) fn from_rust(value: GraphStereoInconsistency) -> Self {
        match value {
            GraphStereoInconsistency::TetrahedralStereoFailure { atom } => {
                Self::TetrahedralStereoFailure { atom: atom.0 }
            }
            GraphStereoInconsistency::StereoAtomFailure { stereo_atom } => {
                Self::StereoAtomFailure {
                    stereo_atom: stereo_atom.0,
                }
            }
            GraphStereoInconsistency::TetrahedralStereoMismatch { atom, stereo_atom } => {
                Self::TetrahedralStereoMismatch {
                    atom: atom.0,
                    stereo_atom: stereo_atom.0,
                }
            }
            GraphStereoInconsistency::CisTransStereoFailure { bond } => {
                Self::CisTransStereoFailure { bond: bond.0 }
            }
            GraphStereoInconsistency::StereoBondFailure { stereo_bond } => {
                Self::StereoBondFailure {
                    stereo_bond: stereo_bond.0,
                }
            }
            GraphStereoInconsistency::CisTransStereoMismatch { bond, stereo_bond } => {
                Self::CisTransStereoMismatch {
                    bond: bond.0,
                    stereo_bond: stereo_bond.0,
                }
            }
        }
    }
    pub(crate) fn to_rust(self) -> GraphStereoInconsistency {
        match self {
            Self::TetrahedralStereoFailure { atom } => {
                GraphStereoInconsistency::TetrahedralStereoFailure {
                    atom: GraphIrAtomId(atom),
                }
            }
            Self::StereoAtomFailure { stereo_atom } => {
                GraphStereoInconsistency::StereoAtomFailure {
                    stereo_atom: GraphIrStereoAtomId(stereo_atom),
                }
            }
            Self::TetrahedralStereoMismatch { atom, stereo_atom } => {
                GraphStereoInconsistency::TetrahedralStereoMismatch {
                    atom: GraphIrAtomId(atom),
                    stereo_atom: GraphIrStereoAtomId(stereo_atom),
                }
            }
            Self::CisTransStereoFailure { bond } => {
                GraphStereoInconsistency::CisTransStereoFailure {
                    bond: GraphIrBondId(bond),
                }
            }
            Self::StereoBondFailure { stereo_bond } => {
                GraphStereoInconsistency::StereoBondFailure {
                    stereo_bond: GraphIrStereoBondId(stereo_bond),
                }
            }
            Self::CisTransStereoMismatch { bond, stereo_bond } => {
                GraphStereoInconsistency::CisTransStereoMismatch {
                    bond: GraphIrBondId(bond),
                    stereo_bond: GraphIrStereoBondId(stereo_bond),
                }
            }
        }
    }
}

#[cfg(test)]
mod tests;
