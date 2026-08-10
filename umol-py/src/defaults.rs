//! Python bindings for DSL-to-graph-IR construction defaults.

use pyo3::prelude::*;
use umol_graph_ir::dsl::{
    MoleculeDefaults as GraphIrMoleculeDefaults, ReactionDefaults as GraphIrReactionDefaults,
};

/// Defaults applied while constructing a molecule from its DSL.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoleculeDefaults(GraphIrMoleculeDefaults);

#[pymethods]
impl MoleculeDefaults {
    /// Preserve every IR value; omitted input values remain undetermined.
    #[new]
    pub(crate) fn new() -> Self {
        Self(GraphIrMoleculeDefaults::new())
    }

    /// Ground ordinary entity fields while leaving constraints required.
    #[staticmethod]
    pub(crate) fn ground() -> Self {
        Self(GraphIrMoleculeDefaults::ground())
    }

    fn __repr__(&self) -> &'static str {
        if self.0 == GraphIrMoleculeDefaults::ground() {
            "MoleculeDefaults.ground()"
        } else {
            "MoleculeDefaults()"
        }
    }
}

impl MoleculeDefaults {
    pub(crate) fn to_rust(&self) -> &GraphIrMoleculeDefaults {
        &self.0
    }
}

/// Defaults applied while constructing a reaction from its DSL.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReactionDefaults(GraphIrReactionDefaults);

#[pymethods]
impl ReactionDefaults {
    /// Preserve every IR value; omitted input values remain undetermined.
    #[new]
    pub(crate) fn new() -> Self {
        Self(GraphIrReactionDefaults::new())
    }

    /// Ground ordinary fields in the LHS and delta entity snapshots.
    #[staticmethod]
    pub(crate) fn ground() -> Self {
        Self(GraphIrReactionDefaults::ground())
    }

    fn __repr__(&self) -> &'static str {
        if self.0 == GraphIrReactionDefaults::ground() {
            "ReactionDefaults.ground()"
        } else {
            "ReactionDefaults()"
        }
    }
}

impl ReactionDefaults {
    pub(crate) fn to_rust(&self) -> &GraphIrReactionDefaults {
        &self.0
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    fn test_molecule_defaults_new() {
        let defaults = MoleculeDefaults::new();

        assert_eq!(defaults.to_rust(), &GraphIrMoleculeDefaults::new());
        assert_eq!(defaults.__repr__(), "MoleculeDefaults()");
    }

    #[rstest]
    fn test_molecule_defaults_ground() {
        let defaults = MoleculeDefaults::ground();

        assert_eq!(defaults.to_rust(), &GraphIrMoleculeDefaults::ground());
        assert_eq!(defaults.__repr__(), "MoleculeDefaults.ground()");
    }

    #[rstest]
    fn test_reaction_defaults_new() {
        let defaults = ReactionDefaults::new();

        assert_eq!(defaults.to_rust(), &GraphIrReactionDefaults::new());
        assert_eq!(defaults.__repr__(), "ReactionDefaults()");
    }

    #[rstest]
    fn test_reaction_defaults_ground() {
        let defaults = ReactionDefaults::ground();

        assert_eq!(defaults.to_rust(), &GraphIrReactionDefaults::ground());
        assert_eq!(defaults.__repr__(), "ReactionDefaults.ground()");
    }
}
