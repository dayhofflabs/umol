//! Python binding for atom-type registries used by valence resolution.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use umol_ast::ast::{ElementAst as AstElementAst, ValueAst as AstValueAst};
use umol_chem::element::{Element as ChemElement, MAX_ATOMIC_NUMBER};
use umol_graph::ops::valence::AtomTypeRegistry as GraphAtomTypeRegistry;

use crate::atom::AtomAst;
use crate::element::Element;

/// An immutable collection of atom patterns used by atom-typing valence resolution.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AtomTypeRegistry(GraphAtomTypeRegistry);

#[pymethods]
impl AtomTypeRegistry {
    /// The built-in atom-type registry.
    #[staticmethod]
    fn default() -> Self {
        Self::from_rust(GraphAtomTypeRegistry::default_registry())
    }

    /// Construct a registry from atom patterns with literal elements and charges.
    #[staticmethod]
    fn from_atoms(py: Python<'_>, atoms: Vec<Py<AtomAst>>) -> PyResult<Self> {
        let mut rust_atoms = Vec::with_capacity(atoms.len());
        for (index, atom) in atoms.iter().enumerate() {
            let atom = atom.bind(py).borrow();
            let atom = atom.inner();
            if !matches!(&atom.element, AstElementAst::Lit(_)) {
                return Err(PyValueError::new_err(format!(
                    "atom type registry entry {index} must have a literal element"
                )));
            }
            let AstValueAst::Lit(charge) = &atom.charge else {
                return Err(PyValueError::new_err(format!(
                    "atom type registry entry {index} must have a literal charge"
                )));
            };
            i8::try_from(*charge).map_err(|_| {
                PyValueError::new_err(format!(
                    "atom type registry entry {index} charge {charge} is outside -128..=127"
                ))
            })?;
            rust_atoms.push(atom.clone());
        }
        Ok(Self(GraphAtomTypeRegistry::from_atoms(rust_atoms)))
    }

    /// Parse a registry from its TOML representation.
    #[staticmethod]
    fn from_toml(input: &str) -> PyResult<Self> {
        GraphAtomTypeRegistry::from_toml_str(input)
            .map(Self)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// The stable content hash of the registry's atom patterns.
    #[getter]
    fn content_hash(&self) -> u64 {
        self.0.content_hash()
    }

    /// The content hash as a fixed-width lowercase hexadecimal string.
    #[getter]
    fn content_hash_hex(&self) -> String {
        self.0.content_hash_hex()
    }

    /// All patterns registered for an element, detached from the registry.
    fn patterns_for_element(&self, element: Element) -> Vec<AtomAst> {
        self.0
            .patterns_for_element(ChemElement::from(&element))
            .iter()
            .cloned()
            .map(AtomAst::from_inner)
            .collect()
    }

    /// Patterns registered for an element and charge, detached from the registry.
    fn patterns_for_element_and_charge(&self, element: Element, charge: i8) -> Vec<AtomAst> {
        self.0
            .patterns_for_element_and_charge(ChemElement::from(&element), charge)
            .iter()
            .cloned()
            .map(AtomAst::from_inner)
            .collect()
    }

    fn __repr__(&self) -> String {
        if &self.0 == GraphAtomTypeRegistry::default_registry() {
            return "AtomTypeRegistry.default()".to_owned();
        }

        let atoms = (1..=MAX_ATOMIC_NUMBER)
            .filter_map(ChemElement::from_atomic_number)
            .flat_map(|element| self.0.patterns_for_element(element))
            .map(|atom| format!("AtomAst.parse({:?})", atom.to_string()))
            .collect::<Vec<_>>()
            .join(", ");
        format!("AtomTypeRegistry.from_atoms([{atoms}])")
    }
}

impl AtomTypeRegistry {
    pub(crate) fn from_rust(registry: &GraphAtomTypeRegistry) -> Self {
        Self(registry.clone())
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for valence resolution configuration"
    )]
    pub(crate) fn to_rust(&self) -> GraphAtomTypeRegistry {
        self.0.clone()
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use pyo3::exceptions::PyValueError;
    use rstest::rstest;
    use umol_ast::ast::AtomAst as AstAtomAst;
    use umol_graph::registry;

    use super::*;

    #[rstest]
    fn test_atom_type_registry_default() {
        assert_eq!(
            AtomTypeRegistry::default().0,
            *GraphAtomTypeRegistry::default_registry()
        );
    }

    #[rstest]
    fn test_atom_type_registry_from_atoms() {
        Python::attach(|py| {
            let atoms = vec![
                "C#c0#v4".parse::<AstAtomAst>().unwrap(),
                "O#c0#v2".parse::<AstAtomAst>().unwrap(),
            ];
            let python_atoms = atoms
                .iter()
                .cloned()
                .map(|atom| Py::new(py, AtomAst::from_inner(atom)).unwrap())
                .collect();

            assert_eq!(
                AtomTypeRegistry::from_atoms(py, python_atoms).unwrap().0,
                GraphAtomTypeRegistry::from_atoms(atoms)
            );
        });
    }

    #[rstest]
    #[case::element(
        AstAtomAst::default().with_charge(0),
        "atom type registry entry 0 must have a literal element"
    )]
    #[case::charge(
        AstAtomAst::from_element(ChemElement::C),
        "atom type registry entry 0 must have a literal charge"
    )]
    #[case::charge_range(
        AstAtomAst::from_element(ChemElement::C).with_charge(128),
        "atom type registry entry 0 charge 128 is outside -128..=127"
    )]
    fn test_atom_type_registry_from_atoms_error(#[case] atom: AstAtomAst, #[case] expected: &str) {
        Python::attach(|py| {
            let atom = Py::new(py, AtomAst::from_inner(atom)).unwrap();
            let error = AtomTypeRegistry::from_atoms(py, vec![atom]).unwrap_err();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected
            );
        });
    }

    #[rstest]
    fn test_atom_type_registry_from_toml() {
        let input = "[C]\n0 = [\"C#c0#v4\"]\n[O]\n0 = [\"O#c0#v2\"]";

        assert_eq!(
            AtomTypeRegistry::from_toml(input).unwrap().0,
            GraphAtomTypeRegistry::from_toml_str(input).unwrap()
        );
    }

    #[rstest]
    fn test_atom_type_registry_from_toml_error() {
        Python::attach(|py| {
            let error = AtomTypeRegistry::from_toml("[X]\n0 = [\"X#c0\"]").unwrap_err();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "invalid atom type registry: unknown element: X"
            );
        });
    }

    #[rstest]
    fn test_atom_type_registry_content_hash() {
        let registry = AtomTypeRegistry(registry!["C#c0#v4", "O#c0#v2"]);

        assert_eq!(registry.content_hash(), registry.0.content_hash());
        assert_eq!(registry.content_hash_hex(), registry.0.content_hash_hex());
        assert_eq!(
            registry.content_hash_hex(),
            format!("{:016x}", registry.content_hash())
        );
    }

    #[rstest]
    fn test_atom_type_registry_patterns_for_element() {
        let registry = AtomTypeRegistry(registry!["C#c0#v4", "C#c+#v3", "O#c0#v2"]);

        assert_eq!(
            registry
                .patterns_for_element(Element::from(ChemElement::C))
                .into_iter()
                .map(|atom| atom.inner().clone())
                .collect::<Vec<_>>(),
            registry.0.patterns_for_element(ChemElement::C)
        );
    }

    #[rstest]
    #[case::neutral(0, registry!["C#c0#v4"].patterns_for_element_and_charge(ChemElement::C, 0).to_vec())]
    #[case::positive(1, registry!["C#c+#v3"].patterns_for_element_and_charge(ChemElement::C, 1).to_vec())]
    #[case::missing(-1, Vec::new())]
    fn test_atom_type_registry_patterns_for_element_and_charge(
        #[case] charge: i8,
        #[case] expected: Vec<AstAtomAst>,
    ) {
        let registry = AtomTypeRegistry(registry!["C#c0#v4", "C#c+#v3", "O#c0#v2"]);

        assert_eq!(
            registry
                .patterns_for_element_and_charge(Element::from(ChemElement::C), charge)
                .into_iter()
                .map(|atom| atom.inner().clone())
                .collect::<Vec<_>>(),
            expected
        );
    }

    #[rstest]
    #[case::default(AtomTypeRegistry::default(), "AtomTypeRegistry.default()")]
    #[case::empty(
        AtomTypeRegistry(GraphAtomTypeRegistry::new()),
        "AtomTypeRegistry.from_atoms([])"
    )]
    #[case::custom(
        AtomTypeRegistry(registry!["C#c0#v4", "O#c0#v2"]),
        "AtomTypeRegistry.from_atoms([AtomAst.parse(\"C#i=#c0#h0#n0#u0#s#v4#d0#t0#a!#m!#T!\"), AtomAst.parse(\"O#i=#c0#h0#n0#u0#s#v2#d0#t0#a!#m!#T!\")])"
    )]
    fn test_atom_type_registry_repr(#[case] registry: AtomTypeRegistry, #[case] expected: &str) {
        assert_eq!(registry.__repr__(), expected);
    }

    #[rstest]
    #[case::borrowed(Cow::Borrowed(GraphAtomTypeRegistry::default_registry()))]
    #[case::owned(Cow::Owned(registry!["C#c0#v4", "O#c0#v2"]))]
    fn test_atom_type_registry_from_rust(#[case] registry: Cow<'static, GraphAtomTypeRegistry>) {
        assert_eq!(AtomTypeRegistry::from_rust(registry.as_ref()).0, *registry);
    }

    #[rstest]
    #[case::default(AtomTypeRegistry::default())]
    #[case::custom(AtomTypeRegistry(registry!["C#c0#v4", "O#c0#v2"]))]
    fn test_atom_type_registry_to_rust(#[case] registry: AtomTypeRegistry) {
        assert_eq!(registry.to_rust(), registry.0);
    }
}
