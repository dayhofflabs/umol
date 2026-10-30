//! Stereo processing sites, selections, and perception queries.

use std::collections::BTreeSet;

use pyo3::basic::CompareOp;
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use umol_graph::ops::stereo::{
    StereoDerivation as GraphStereoDerivation, StereoInconsistency as GraphStereoInconsistency,
    StereoPerception as GraphStereoPerception, StereoPerceptionError as GraphStereoPerceptionError,
    StereoSelection as GraphStereoSelection, StereoSite as GraphStereoSite,
};
use umol_graph_ir::ir::{
    AtomId as GraphIrAtomId, BondId as GraphIrBondId, StereoAtomId as GraphIrStereoAtomId,
    StereoBondId as GraphIrStereoBondId,
};

use super::{
    StereoAtomForm, StereoBondForm, StereoCoset, StereoKind, StereoLigand, Stereogenicity,
};
use crate::convert::{hash_rust, FrozenSetItems, TupleItems};
use crate::entity::EntityForm;
use crate::error::stereo_perception_error;
use crate::model::stereo::StereoModel;
use crate::molecule::Molecule;

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

/// Immutable stereo derivation with ordered rows and retained inconsistencies.
///
/// Properties return tuples and read-only copies of the entity forms. Construction
/// preserves the supplied rows without checking them against a molecule.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StereoDerivation(GraphStereoDerivation);

#[pymethods]
impl StereoDerivation {
    #[new]
    #[pyo3(signature = (*, stereo_atoms=Vec::new(), stereo_bonds=Vec::new(), skipped_stereo_bonds=Vec::new(), inconsistencies=Vec::new()))]
    fn new(
        py: Python<'_>,
        stereo_atoms: Vec<(u32, Vec<StereoLigand>, Py<StereoAtomForm>)>,
        stereo_bonds: Vec<(u32, Vec<StereoLigand>, Py<StereoBondForm>)>,
        skipped_stereo_bonds: Vec<u32>,
        inconsistencies: Vec<StereoInconsistency>,
    ) -> Self {
        Self::from_rust(GraphStereoDerivation {
            stereo_atoms: stereo_atoms
                .into_iter()
                .map(|(id, ligands, form)| {
                    (
                        GraphIrAtomId(id),
                        ligands.into_iter().map(StereoLigand::to_rust).collect(),
                        form.borrow(py).to_rust().clone(),
                    )
                })
                .collect(),
            stereo_bonds: stereo_bonds
                .into_iter()
                .map(|(id, ligands, form)| {
                    (
                        GraphIrBondId(id),
                        ligands.into_iter().map(StereoLigand::to_rust).collect(),
                        form.borrow(py).to_rust().clone(),
                    )
                })
                .collect(),
            skipped_stereo_bonds: skipped_stereo_bonds
                .into_iter()
                .map(GraphIrBondId)
                .collect(),
            inconsistencies: inconsistencies
                .into_iter()
                .map(StereoInconsistency::to_rust)
                .collect(),
        })
    }

    #[getter]
    fn stereo_atoms<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let values = self
            .to_rust()
            .stereo_atoms
            .iter()
            .map(|(id, ligands, form)| {
                Ok((
                    id.0,
                    TupleItems(
                        ligands
                            .iter()
                            .copied()
                            .map(StereoLigand::from_rust)
                            .collect::<Vec<_>>(),
                    ),
                    StereoAtomForm::new_readonly(py, form.clone())?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, values)
    }

    #[getter]
    fn stereo_bonds<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let values = self
            .to_rust()
            .stereo_bonds
            .iter()
            .map(|(id, ligands, form)| {
                Ok((
                    id.0,
                    TupleItems(
                        ligands
                            .iter()
                            .copied()
                            .map(StereoLigand::from_rust)
                            .collect::<Vec<_>>(),
                    ),
                    StereoBondForm::new_readonly(py, form.clone())?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        PyTuple::new(py, values)
    }

    #[getter]
    fn skipped_stereo_bonds(&self) -> TupleItems<Vec<u32>> {
        self.to_rust()
            .skipped_stereo_bonds
            .iter()
            .map(|id| id.0)
            .collect()
    }

    #[getter]
    fn inconsistencies(&self) -> TupleItems<Vec<StereoInconsistency>> {
        self.to_rust()
            .inconsistencies
            .iter()
            .copied()
            .map(StereoInconsistency::from_rust)
            .collect()
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!("StereoDerivation(stereo_atoms={}, stereo_bonds={}, skipped_stereo_bonds={}, inconsistencies={})",
            self.stereo_atoms(py)?.repr()?.extract::<String>()?,
            self.stereo_bonds(py)?.repr()?.extract::<String>()?,
            self.skipped_stereo_bonds().into_pyobject(py)?.repr()?.extract::<String>()?,
            self.inconsistencies().into_pyobject(py)?.repr()?.extract::<String>()?))
    }
}

impl StereoDerivation {
    pub(crate) fn from_rust(derivation: GraphStereoDerivation) -> Self {
        Self(derivation)
    }
    pub(crate) fn to_rust(&self) -> &GraphStereoDerivation {
        &self.0
    }
}

/// Stereo perception using a snapshot of the supplied stereo model.
#[pyclass(frozen, from_py_object)]
#[derive(Clone, Debug)]
pub struct StereoPerception(GraphStereoPerception);

#[pymethods]
impl StereoPerception {
    #[new]
    fn new(model: &StereoModel) -> Self {
        Self::from_rust(GraphStereoPerception::new(&model.to_rust()))
    }

    /// Derive source stereo evidence and check existing arrangements.
    fn derive(&self, molecule: &Molecule) -> StereoDerivation {
        StereoDerivation::from_rust(self.to_rust().derive(molecule.to_rust()))
    }

    /// Sites skipped by the model's minimum stereo-bond ring size.
    fn skipped_stereo_bonds(&self, molecule: &Molecule) -> FrozenSetItems<BTreeSet<u32>> {
        self.to_rust()
            .skipped_stereo_bonds(molecule.to_rust())
            .into_iter()
            .map(|id| id.0)
            .collect()
    }

    /// Derive a tetrahedral frame, or return None for an excluded or incompatible site.
    fn derive_stereo_atom(
        &self,
        py: Python<'_>,
        molecule: &Molecule,
        atom: u32,
        coset: &StereoCoset,
    ) -> PyResult<Option<(TupleItems<Vec<StereoLigand>>, StereoAtomForm)>> {
        if molecule
            .to_rust()
            .atoms()
            .get(GraphIrAtomId(atom))
            .is_none()
        {
            return Err(stereo_perception_error(
                py,
                GraphStereoPerceptionError::InvalidSite {
                    site: GraphStereoSite::Atom(GraphIrAtomId(atom)),
                },
            ));
        }
        Ok(self
            .to_rust()
            .derive_stereo_atom(molecule.to_rust(), GraphIrAtomId(atom), &coset.to_rust(py))
            .map(|(ligands, form)| {
                (
                    ligands.into_iter().map(StereoLigand::from_rust).collect(),
                    StereoAtomForm::from_rust(form),
                )
            }))
    }

    /// Derive a cis/trans frame, or return None for an excluded or incompatible site.
    fn derive_stereo_bond(
        &self,
        py: Python<'_>,
        molecule: &Molecule,
        bond: u32,
        coset: &StereoCoset,
    ) -> PyResult<Option<(TupleItems<Vec<StereoLigand>>, StereoBondForm)>> {
        if molecule
            .to_rust()
            .bonds()
            .get(GraphIrBondId(bond))
            .is_none()
        {
            return Err(stereo_perception_error(
                py,
                GraphStereoPerceptionError::InvalidSite {
                    site: GraphStereoSite::Bond(GraphIrBondId(bond)),
                },
            ));
        }
        Ok(self
            .to_rust()
            .derive_stereo_bond(molecule.to_rust(), GraphIrBondId(bond), &coset.to_rust(py))
            .map(|(ligands, form)| {
                (
                    ligands.into_iter().map(StereoLigand::from_rust).collect(),
                    StereoBondForm::from_rust(form),
                )
            }))
    }
}

impl StereoPerception {
    pub(crate) fn from_rust(perception: GraphStereoPerception) -> Self {
        Self(perception)
    }
    pub(crate) fn to_rust(&self) -> &GraphStereoPerception {
        &self.0
    }
}

#[cfg(test)]
mod tests;
