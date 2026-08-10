//! Electron-counts form: `Undetermined`, or a positional per-member-atom
//! count vector. A value leaf shared by the aromatic-system and multicenter-bond
//! bindings; the vector is positional (cell = member atom), aligned to the owning
//! entity's participant order.

use pyo3::prelude::*;
use umol_graph_ir::ir::{AsLit, ElectronCountsForm as GraphIrElectronCountsForm};

use crate::convert::{hash_rust, variant_repr};
use crate::lattice::impl_py_lattice;

/// A per-member-atom electron-count vector: undetermined, or a concrete list of
/// counts positionally aligned to the owning entity's atoms.
#[pyclass]
pub enum ElectronCountsForm {
    Undetermined(),
    Lit(Vec<i64>),
}

#[pymethods]
impl ElectronCountsForm {
    /// The concrete count vector, or `None` when undetermined.
    fn as_lit(&self) -> Option<Vec<i64>> {
        self.to_rust().as_lit()
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.to_rust() == other.to_rust()
    }

    fn __hash__(&self) -> u64 {
        hash_rust(&self.to_rust())
    }

    fn __repr__(slf: Py<Self>, py: Python<'_>) -> PyResult<String> {
        let (variant, arity) = match &*slf.bind(py).borrow() {
            ElectronCountsForm::Undetermined() => ("Undetermined", 0),
            ElectronCountsForm::Lit(_) => ("Lit", 1),
        };
        variant_repr(slf.bind(py).as_any(), "ElectronCountsForm", variant, arity)
    }
}

impl ElectronCountsForm {
    pub(crate) fn from_rust(ast: &GraphIrElectronCountsForm) -> Self {
        match ast {
            GraphIrElectronCountsForm::Undetermined => Self::Undetermined(),
            GraphIrElectronCountsForm::Lit(counts) => Self::Lit(counts.clone()),
        }
    }

    pub(crate) fn to_rust(&self) -> GraphIrElectronCountsForm {
        match self {
            Self::Undetermined() => GraphIrElectronCountsForm::Undetermined,
            Self::Lit(counts) => GraphIrElectronCountsForm::Lit(counts.clone()),
        }
    }
}

impl_py_lattice!(
    ElectronCountsForm,
    GraphIrElectronCountsForm,
    |value: &ElectronCountsForm, _py: Python<'_>| -> PyResult<GraphIrElectronCountsForm> {
        Ok(value.to_rust())
    },
    |_py: Python<'_>, value: GraphIrElectronCountsForm| -> PyResult<ElectronCountsForm> {
        Ok(ElectronCountsForm::from_rust(&value))
    }
);

/// Setter coercion for an electron-counts field: a Python `list[int]` → `Lit`, or an
/// `ElectronCountsForm` passthrough (matching `impl From<Vec<i64>>`).
#[derive(FromPyObject)]
pub(crate) enum ElectronCountsLike {
    Ast(Py<ElectronCountsForm>),
    Lit(Vec<i64>),
}

impl ElectronCountsLike {
    pub(crate) fn to_rust(&self, py: Python<'_>) -> GraphIrElectronCountsForm {
        match self {
            ElectronCountsLike::Lit(counts) => GraphIrElectronCountsForm::Lit(counts.clone()),
            ElectronCountsLike::Ast(a) => a.bind(py).borrow().to_rust(),
        }
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(GraphIrElectronCountsForm::Undetermined)]
    #[case(GraphIrElectronCountsForm::Lit(vec![1, 1, 1, 1, 1, 1]))]
    #[case(GraphIrElectronCountsForm::Lit(vec![]))]
    fn test_electron_counts_form_roundtrip(#[case] ast: GraphIrElectronCountsForm) {
        assert_eq!(ElectronCountsForm::from_rust(&ast).to_rust(), ast);
    }

    #[rstest]
    #[case(GraphIrElectronCountsForm::Undetermined, None)]
    #[case(GraphIrElectronCountsForm::Lit(vec![2, 0, 2]), Some(vec![2, 0, 2]))]
    fn test_electron_counts_form_as_lit(
        #[case] ast: GraphIrElectronCountsForm,
        #[case] expected: Option<Vec<i64>>,
    ) {
        assert_eq!(ElectronCountsForm::from_rust(&ast).as_lit(), expected);
    }

    #[rstest]
    fn test_electron_counts_like_to_rust() {
        Python::attach(|py| {
            // a bare list coerces to Lit
            assert_eq!(
                ElectronCountsLike::Lit(vec![1, 0, 1]).to_rust(py),
                GraphIrElectronCountsForm::Lit(vec![1, 0, 1])
            );
            // an ElectronCountsForm passes through
            let ast = Py::new(py, ElectronCountsForm::Lit(vec![2, 2])).unwrap();
            assert_eq!(
                ElectronCountsLike::Ast(ast).to_rust(py),
                GraphIrElectronCountsForm::Lit(vec![2, 2])
            );
        });
    }
}
