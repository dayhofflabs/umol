//! Exact and AST unpaired-electron values.
#![allow(clippy::absolute_paths)] // the `#[pyclass(hash)]` macro expands to absolute paths

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use umol_chem::spin::{SpinState as ChemSpinState, UnpairedElectrons as ChemUnpairedElectrons};
use umol_graph_ir::ir::{
    AsLit, UnpairedElectronsAst as GraphIrUnpairedElectronsAst,
    UnpairedElectronsUpdate as GraphIrUnpairedElectronsUpdate,
};

use crate::convert::{hash_rust, into_py_variant};
use crate::lattice::impl_py_lattice;
use crate::value::{ValueAst, ValueLike};

/// Exact unpaired-electron count and spin multiplicity without physical validation.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct UnpairedElectrons(ChemUnpairedElectrons);

#[pymethods]
impl UnpairedElectrons {
    #[new]
    fn new(count: i64, multiplicity: i64) -> Self {
        Self::from_rust(ChemUnpairedElectrons {
            count,
            multiplicity,
        })
    }

    #[getter]
    fn count(&self) -> i64 {
        self.0.count
    }

    #[getter]
    fn multiplicity(&self) -> i64 {
        self.0.multiplicity
    }

    fn __repr__(&self) -> String {
        let unpaired_electrons = self.to_rust();
        format!(
            "UnpairedElectrons(count={}, multiplicity={})",
            unpaired_electrons.count, unpaired_electrons.multiplicity,
        )
    }
}

impl UnpairedElectrons {
    pub(crate) fn from_rust(unpaired_electrons: ChemUnpairedElectrons) -> Self {
        Self(unpaired_electrons)
    }

    pub(crate) fn to_rust(self) -> ChemUnpairedElectrons {
        self.0
    }
}

/// A physically valid unpaired-electron count and spin multiplicity.
#[pyclass(eq, hash, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct SpinState(ChemSpinState);

#[pymethods]
impl SpinState {
    #[new]
    #[pyo3(signature = (*, unpaired_electrons, multiplicity))]
    fn new(unpaired_electrons: i64, multiplicity: i64) -> PyResult<Self> {
        ChemSpinState::try_from(ChemUnpairedElectrons {
            count: unpaired_electrons,
            multiplicity,
        })
        .map(Self::from_rust)
        .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    #[getter]
    fn unpaired_electrons(&self) -> u8 {
        self.0.unpaired_electrons()
    }

    #[getter]
    fn multiplicity(&self) -> u8 {
        self.0.multiplicity().into()
    }

    fn __repr__(&self) -> String {
        let spin_state = self.to_rust();
        format!(
            "SpinState(unpaired_electrons={}, multiplicity={})",
            spin_state.unpaired_electrons(),
            u8::from(spin_state.multiplicity()),
        )
    }
}

impl SpinState {
    pub(crate) fn from_rust(spin_state: ChemSpinState) -> Self {
        Self(spin_state)
    }

    pub(crate) fn to_rust(self) -> ChemSpinState {
        self.0
    }
}

/// Unpaired-electron count and multiplicity as independent value fields.
#[pyclass]
pub struct UnpairedElectronsAst {
    #[pyo3(get)]
    count: Py<ValueAst>,
    #[pyo3(get)]
    multiplicity: Py<ValueAst>,
}

#[pymethods]
impl UnpairedElectronsAst {
    #[new]
    fn new(py: Python<'_>, count: ValueLike, multiplicity: ValueLike) -> PyResult<Self> {
        Ok(UnpairedElectronsAst {
            count: count.to_py(py)?,
            multiplicity: multiplicity.to_py(py)?,
        })
    }

    /// Return the exact pair when both component expressions are literal.
    fn as_lit(&self, py: Python<'_>) -> Option<UnpairedElectrons> {
        Some(UnpairedElectrons::from_rust(self.to_rust(py).as_lit()?))
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "UnpairedElectronsAst({}, {})",
            self.count.bind(py).as_any().repr()?.extract::<String>()?,
            self.multiplicity
                .bind(py)
                .as_any()
                .repr()?
                .extract::<String>()?,
        ))
    }
}

impl UnpairedElectronsAst {
    pub(crate) fn from_rust(
        py: Python<'_>,
        ast: &GraphIrUnpairedElectronsAst,
    ) -> PyResult<UnpairedElectronsAst> {
        Ok(UnpairedElectronsAst {
            count: into_py_variant(py, ValueAst::from_rust(py, &ast.count)?)?,
            multiplicity: into_py_variant(py, ValueAst::from_rust(py, &ast.multiplicity)?)?,
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrUnpairedElectronsAst {
        GraphIrUnpairedElectronsAst {
            count: self.count.bind(py).borrow().to_rust(py),
            multiplicity: self.multiplicity.bind(py).borrow().to_rust(py),
        }
    }
}

impl_py_lattice!(
    UnpairedElectronsAst,
    GraphIrUnpairedElectronsAst,
    |value: &UnpairedElectronsAst, py: Python<'_>| -> PyResult<GraphIrUnpairedElectronsAst> {
        Ok(value.to_rust(py))
    },
    |py: Python<'_>, value: GraphIrUnpairedElectronsAst| -> PyResult<UnpairedElectronsAst> {
        UnpairedElectronsAst::from_rust(py, &value)
    }
);

/// Independent updates to the unpaired-electron count and multiplicity.
///
/// `None` leaves a component unchanged. A `ValueAst`, including
/// `ValueAst.Undetermined()`, sets that component exactly.
#[pyclass]
pub struct UnpairedElectronsUpdate {
    count: Option<Py<ValueAst>>,
    multiplicity: Option<Py<ValueAst>>,
}

#[pymethods]
impl UnpairedElectronsUpdate {
    #[new]
    #[pyo3(signature = (*, count=None, multiplicity=None))]
    fn new(
        py: Python<'_>,
        count: Option<ValueLike>,
        multiplicity: Option<ValueLike>,
    ) -> PyResult<Self> {
        Self::from_rust(
            py,
            &GraphIrUnpairedElectronsUpdate {
                count: count.map(|value| value.to_rust(py)),
                multiplicity: multiplicity.map(|value| value.to_rust(py)),
            },
        )
    }

    #[getter]
    fn count(&self, py: Python<'_>) -> Option<Py<ValueAst>> {
        self.count.as_ref().map(|value| value.clone_ref(py))
    }

    #[getter]
    fn multiplicity(&self, py: Python<'_>) -> Option<Py<ValueAst>> {
        self.multiplicity.as_ref().map(|value| value.clone_ref(py))
    }

    fn __eq__(&self, other: &Self, py: Python<'_>) -> bool {
        self.to_rust(py) == other.to_rust(py)
    }

    fn __hash__(&self, py: Python<'_>) -> u64 {
        hash_rust(&self.to_rust(py))
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let count = match &self.count {
            Some(value) => value.bind(py).as_any().repr()?.extract::<String>()?,
            None => "None".to_string(),
        };
        let multiplicity = match &self.multiplicity {
            Some(value) => value.bind(py).as_any().repr()?.extract::<String>()?,
            None => "None".to_string(),
        };
        Ok(format!(
            "UnpairedElectronsUpdate(count={count}, multiplicity={multiplicity})"
        ))
    }
}

impl UnpairedElectronsUpdate {
    pub(crate) fn from_rust(
        py: Python<'_>,
        update: &GraphIrUnpairedElectronsUpdate,
    ) -> PyResult<Self> {
        Ok(Self {
            count: update
                .count
                .as_ref()
                .map(|value| {
                    ValueAst::from_rust(py, value).and_then(|value| into_py_variant(py, value))
                })
                .transpose()?,
            multiplicity: update
                .multiplicity
                .as_ref()
                .map(|value| {
                    ValueAst::from_rust(py, value).and_then(|value| into_py_variant(py, value))
                })
                .transpose()?,
        })
    }

    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrUnpairedElectronsUpdate {
        GraphIrUnpairedElectronsUpdate {
            count: self
                .count
                .as_ref()
                .map(|value| value.bind(py).borrow().to_rust(py)),
            multiplicity: self
                .multiplicity
                .as_ref()
                .map(|value| value.bind(py).borrow().to_rust(py)),
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_chem::spin::{
        SpinMultiplicity, SpinState as ChemSpinState, UnpairedElectrons as ChemUnpairedElectrons,
    };
    use umol_graph_ir::ir::{
        UnpairedElectronsAst as GraphIrUnpairedElectronsAst,
        UnpairedElectronsUpdate as GraphIrUnpairedElectronsUpdate, ValueAst as GraphIrValueAst,
    };

    use super::*;

    #[rstest]
    #[case::closed_shell(ChemUnpairedElectrons { count: 0, multiplicity: 1 })]
    #[case::open_shell(ChemUnpairedElectrons { count: 2, multiplicity: 3 })]
    #[case::physics_invalid(ChemUnpairedElectrons { count: -1, multiplicity: 0 })]
    fn test_unpaired_electrons_roundtrip(#[case] unpaired_electrons: ChemUnpairedElectrons) {
        assert_eq!(
            UnpairedElectrons::from_rust(unpaired_electrons).to_rust(),
            unpaired_electrons,
        );
    }

    #[rstest]
    #[case::closed_shell(ChemSpinState::closed_shell())]
    #[case::doublet(ChemSpinState::new(1, SpinMultiplicity::DOUBLET).unwrap())]
    #[case::open_shell_singlet(ChemSpinState::new(2, SpinMultiplicity::SINGLET).unwrap())]
    #[case::triplet(ChemSpinState::new(2, SpinMultiplicity::TRIPLET).unwrap())]
    fn test_spin_state_roundtrip(#[case] spin_state: ChemSpinState) {
        assert_eq!(SpinState::from_rust(spin_state).to_rust(), spin_state);
    }

    #[rstest]
    #[case::complete(
        GraphIrUnpairedElectronsAst {
            count: GraphIrValueAst::Lit(2),
            multiplicity: GraphIrValueAst::Lit(3),
        },
        Some(ChemUnpairedElectrons { count: 2, multiplicity: 3 }),
    )]
    #[case::physics_invalid(
        GraphIrUnpairedElectronsAst {
            count: GraphIrValueAst::Lit(2),
            multiplicity: GraphIrValueAst::Lit(2),
        },
        Some(ChemUnpairedElectrons { count: 2, multiplicity: 2 }),
    )]
    #[case::count_partial(
        GraphIrUnpairedElectronsAst {
            count: GraphIrValueAst::Undetermined,
            multiplicity: GraphIrValueAst::Lit(3),
        },
        None,
    )]
    #[case::multiplicity_partial(
        GraphIrUnpairedElectronsAst {
            count: GraphIrValueAst::Lit(2),
            multiplicity: GraphIrValueAst::Undetermined,
        },
        None,
    )]
    fn test_unpaired_electrons_ast_as_lit(
        #[case] ast: GraphIrUnpairedElectronsAst,
        #[case] expected: Option<ChemUnpairedElectrons>,
    ) {
        Python::attach(|py| {
            let ast = UnpairedElectronsAst::from_rust(py, &ast).unwrap();
            assert_eq!(ast.as_lit(py).map(UnpairedElectrons::to_rust), expected);
        });
    }

    #[rstest]
    #[case::complete(GraphIrUnpairedElectronsAst {
        count: GraphIrValueAst::Lit(1),
        multiplicity: GraphIrValueAst::Lit(2),
    })]
    #[case::physics_invalid(GraphIrUnpairedElectronsAst {
        count: GraphIrValueAst::Lit(2),
        multiplicity: GraphIrValueAst::Lit(2),
    })]
    #[case::partial(GraphIrUnpairedElectronsAst {
        count: GraphIrValueAst::Undetermined,
        multiplicity: GraphIrValueAst::Undetermined,
    })]
    fn test_unpaired_electrons_ast_roundtrip(#[case] ast: GraphIrUnpairedElectronsAst) {
        Python::attach(|py| {
            assert_eq!(
                UnpairedElectronsAst::from_rust(py, &ast)
                    .unwrap()
                    .to_rust(py),
                ast
            );
        });
    }

    #[rstest]
    #[case::empty(GraphIrUnpairedElectronsUpdate::default())]
    #[case::single(GraphIrUnpairedElectronsUpdate {
        count: Some(GraphIrValueAst::Lit(2)),
        multiplicity: None,
    })]
    #[case::both(GraphIrUnpairedElectronsUpdate {
        count: Some(GraphIrValueAst::Lit(2)),
        multiplicity: Some(GraphIrValueAst::Lit(3)),
    })]
    #[case::explicit_undetermined(GraphIrUnpairedElectronsUpdate {
        count: Some(GraphIrValueAst::Undetermined),
        multiplicity: Some(GraphIrValueAst::Undetermined),
    })]
    fn test_unpaired_electrons_update_roundtrip(#[case] update: GraphIrUnpairedElectronsUpdate) {
        Python::attach(|py| {
            assert_eq!(
                UnpairedElectronsUpdate::from_rust(py, &update)
                    .unwrap()
                    .to_rust(py),
                update,
            );
        });
    }
}
