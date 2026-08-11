//! Python bindings for connectivity-model values.

use pyo3::prelude::*;
use umol_graph::ops::validate::ConnectivityModel as GraphConnectivityModel;

/// Connectivity requirements used by molecule conformance validation.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ConnectivityModel(GraphConnectivityModel);

#[pymethods]
impl ConnectivityModel {
    #[new]
    #[pyo3(signature = (*, allow_disconnected, allow_disconnected_dative, allow_disconnected_aromatic, allow_disconnected_multicenter, allow_disconnected_noncovalent, allow_disconnected_stereo_atom, allow_disconnected_stereo_bond, allow_disconnected_constraints))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        allow_disconnected: bool,
        allow_disconnected_dative: bool,
        allow_disconnected_aromatic: bool,
        allow_disconnected_multicenter: bool,
        allow_disconnected_noncovalent: bool,
        allow_disconnected_stereo_atom: bool,
        allow_disconnected_stereo_bond: bool,
        allow_disconnected_constraints: bool,
    ) -> Self {
        Self(GraphConnectivityModel {
            allow_disconnected,
            allow_disconnected_dative,
            allow_disconnected_aromatic,
            allow_disconnected_multicenter,
            allow_disconnected_noncovalent,
            allow_disconnected_stereo_atom,
            allow_disconnected_stereo_bond,
            allow_disconnected_constraints,
        })
    }

    #[staticmethod]
    fn default() -> Self {
        Self(GraphConnectivityModel::default())
    }

    #[getter]
    fn allow_disconnected(&self) -> bool {
        self.0.allow_disconnected
    }

    #[getter]
    fn allow_disconnected_dative(&self) -> bool {
        self.0.allow_disconnected_dative
    }

    #[getter]
    fn allow_disconnected_aromatic(&self) -> bool {
        self.0.allow_disconnected_aromatic
    }

    #[getter]
    fn allow_disconnected_multicenter(&self) -> bool {
        self.0.allow_disconnected_multicenter
    }

    #[getter]
    fn allow_disconnected_noncovalent(&self) -> bool {
        self.0.allow_disconnected_noncovalent
    }

    #[getter]
    fn allow_disconnected_stereo_atom(&self) -> bool {
        self.0.allow_disconnected_stereo_atom
    }

    #[getter]
    fn allow_disconnected_stereo_bond(&self) -> bool {
        self.0.allow_disconnected_stereo_bond
    }

    #[getter]
    fn allow_disconnected_constraints(&self) -> bool {
        self.0.allow_disconnected_constraints
    }

    pub(crate) fn __repr__(&self) -> String {
        if self.0 == GraphConnectivityModel::default() {
            return "ConnectivityModel.default()".to_owned();
        }
        format!(
            "ConnectivityModel(allow_disconnected={}, allow_disconnected_dative={}, allow_disconnected_aromatic={}, allow_disconnected_multicenter={}, allow_disconnected_noncovalent={}, allow_disconnected_stereo_atom={}, allow_disconnected_stereo_bond={}, allow_disconnected_constraints={})",
            python_bool(self.0.allow_disconnected),
            python_bool(self.0.allow_disconnected_dative),
            python_bool(self.0.allow_disconnected_aromatic),
            python_bool(self.0.allow_disconnected_multicenter),
            python_bool(self.0.allow_disconnected_noncovalent),
            python_bool(self.0.allow_disconnected_stereo_atom),
            python_bool(self.0.allow_disconnected_stereo_bond),
            python_bool(self.0.allow_disconnected_constraints),
        )
    }
}

impl ConnectivityModel {
    pub(crate) fn from_rust(model: &GraphConnectivityModel) -> Self {
        Self(*model)
    }

    pub(crate) fn to_rust(self) -> GraphConnectivityModel {
        self.0
    }
}

fn python_bool(value: bool) -> &'static str {
    if value {
        "True"
    } else {
        "False"
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::default(GraphConnectivityModel::default(), "ConnectivityModel.default()")]
    #[case::strict(
        GraphConnectivityModel {
            allow_disconnected: false,
            allow_disconnected_dative: false,
            allow_disconnected_aromatic: false,
            allow_disconnected_multicenter: false,
            allow_disconnected_noncovalent: false,
            allow_disconnected_stereo_atom: false,
            allow_disconnected_stereo_bond: false,
            allow_disconnected_constraints: false,
        },
        "ConnectivityModel(allow_disconnected=False, allow_disconnected_dative=False, allow_disconnected_aromatic=False, allow_disconnected_multicenter=False, allow_disconnected_noncovalent=False, allow_disconnected_stereo_atom=False, allow_disconnected_stereo_bond=False, allow_disconnected_constraints=False)",
    )]
    fn test_connectivity_model_from_rust(
        #[case] model: GraphConnectivityModel,
        #[case] expected_repr: &str,
    ) {
        let python = ConnectivityModel::from_rust(&model);

        assert_eq!(python.to_rust(), model);
        assert_eq!(python.__repr__(), expected_repr);
    }
}
