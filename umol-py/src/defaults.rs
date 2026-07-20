//! Python bindings for DSL-to-AST construction defaults.

use pyo3::prelude::*;
use umol_ast::dsl::{
    MoleculeDefaults as AstMoleculeDefaults, ReactionDefaults as AstReactionDefaults,
};

/// Defaults applied while constructing a molecule AST from its DSL.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MoleculeDefaults(AstMoleculeDefaults);

#[pymethods]
impl MoleculeDefaults {
    /// Preserve omitted values as undetermined.
    #[new]
    pub(crate) fn new() -> Self {
        Self(AstMoleculeDefaults::new())
    }

    /// Ground ordinary entity fields while leaving constraints required.
    #[staticmethod]
    pub(crate) fn ground() -> Self {
        Self(AstMoleculeDefaults::ground())
    }

    fn __repr__(&self) -> &'static str {
        if self.0 == AstMoleculeDefaults::ground() {
            "MoleculeDefaults.ground()"
        } else {
            "MoleculeDefaults()"
        }
    }
}

impl MoleculeDefaults {
    pub(crate) fn to_rust(&self) -> AstMoleculeDefaults {
        self.0.clone()
    }
}

/// Defaults applied while constructing a reaction AST from its DSL.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReactionDefaults(AstReactionDefaults);

#[pymethods]
impl ReactionDefaults {
    /// Preserve omitted values as undetermined.
    #[new]
    pub(crate) fn new() -> Self {
        Self(AstReactionDefaults::new())
    }

    /// Ground ordinary fields in the LHS and delta entity snapshots.
    #[staticmethod]
    pub(crate) fn ground() -> Self {
        Self(AstReactionDefaults::ground())
    }

    fn __repr__(&self) -> &'static str {
        if self.0 == AstReactionDefaults::ground() {
            "ReactionDefaults.ground()"
        } else {
            "ReactionDefaults()"
        }
    }
}

impl ReactionDefaults {
    pub(crate) fn to_rust(&self) -> AstReactionDefaults {
        self.0.clone()
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    fn test_molecule_defaults_new() {
        let defaults = MoleculeDefaults::new();

        assert_eq!(defaults.to_rust(), AstMoleculeDefaults::new());
        assert_eq!(defaults.__repr__(), "MoleculeDefaults()");
    }

    #[rstest]
    fn test_molecule_defaults_ground() {
        let defaults = MoleculeDefaults::ground();

        assert_eq!(defaults.to_rust(), AstMoleculeDefaults::ground());
        assert_eq!(defaults.__repr__(), "MoleculeDefaults.ground()");
    }

    #[rstest]
    fn test_reaction_defaults_new() {
        let defaults = ReactionDefaults::new();

        assert_eq!(defaults.to_rust(), AstReactionDefaults::new());
        assert_eq!(defaults.__repr__(), "ReactionDefaults()");
    }

    #[rstest]
    fn test_reaction_defaults_ground() {
        let defaults = ReactionDefaults::ground();

        assert_eq!(defaults.to_rust(), AstReactionDefaults::ground());
        assert_eq!(defaults.__repr__(), "ReactionDefaults.ground()");
    }
}
