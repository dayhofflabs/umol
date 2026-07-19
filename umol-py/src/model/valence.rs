//! Python binding for atom-type registries used by valence resolution.

use std::borrow::Cow;
use std::collections::BTreeMap;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyList;
use umol_ast::ast::{ElementAst as AstElementAst, ValueAst as AstValueAst};
use umol_chem::element::{Element as ChemElement, MAX_ATOMIC_NUMBER};
use umol_graph::ops::model::ValenceModel as GraphValenceModel;
use umol_graph::ops::valence::{
    AtomTypeRegistry as GraphAtomTypeRegistry, ValenceEntry as GraphValenceEntry,
    ValenceTable as GraphValenceTable,
};

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
        reason = "Python-to-Rust conversion API for ChemistryModel configuration"
    )]
    pub(crate) fn to_rust(&self) -> GraphAtomTypeRegistry {
        self.0.clone()
    }
}

/// Valence states for one element in a counts-based valence table.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ValenceEntry(GraphValenceEntry);

#[pymethods]
impl ValenceEntry {
    #[new]
    #[pyo3(signature = (*, target_covalences=Vec::new(), aromatic_valences=Vec::new()))]
    fn new(mut target_covalences: Vec<u8>, mut aromatic_valences: Vec<u8>) -> Self {
        target_covalences.sort_unstable();
        aromatic_valences.sort_unstable();
        Self(GraphValenceEntry {
            target_covalences,
            aromatic_valences,
        })
    }

    /// Lewis/Langmuir saturation targets, in ascending order.
    #[getter]
    fn target_covalences<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        PyList::new(py, &self.0.target_covalences)
    }

    /// Admissible aromatic valences, in ascending order.
    #[getter]
    fn aromatic_valences<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {
        PyList::new(py, &self.0.aromatic_valences)
    }

    fn __repr__(&self) -> String {
        format!(
            "ValenceEntry(target_covalences={:?}, aromatic_valences={:?})",
            self.0.target_covalences, self.0.aromatic_valences
        )
    }
}

impl ValenceEntry {
    pub(crate) fn from_rust(entry: &GraphValenceEntry) -> Self {
        Self(entry.clone())
    }

    pub(crate) fn to_rust(&self) -> GraphValenceEntry {
        self.0.clone()
    }
}

/// An immutable element-to-valence-entry table used by counts-based resolution.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ValenceTable(GraphValenceTable);

#[pymethods]
impl ValenceTable {
    #[new]
    #[pyo3(signature = (*, entries))]
    fn new(entries: BTreeMap<Element, ValenceEntry>) -> Self {
        let mut table = GraphValenceTable::empty();
        for (element, entry) in entries {
            table.insert(ChemElement::from(&element), entry.to_rust());
        }
        Self(table)
    }

    /// The built-in valence table.
    #[staticmethod]
    fn default() -> Self {
        Self::from_rust(GraphValenceTable::default_table())
    }

    /// Parse a valence table from its TOML representation.
    #[staticmethod]
    fn from_toml(input: &str) -> PyResult<Self> {
        GraphValenceTable::from_toml_str(input)
            .map(Self)
            .map_err(|error| PyValueError::new_err(error.to_string()))
    }

    /// The stable content hash of the table entries.
    #[getter]
    fn content_hash(&self) -> u64 {
        self.0.content_hash()
    }

    /// The content hash as a fixed-width lowercase hexadecimal string.
    #[getter]
    fn content_hash_hex(&self) -> String {
        self.0.content_hash_hex()
    }

    /// The entry for an element, detached from the table.
    fn entry(&self, element: Element) -> Option<ValenceEntry> {
        self.0
            .entry(ChemElement::from(&element))
            .map(ValenceEntry::from_rust)
    }

    fn __repr__(&self) -> String {
        if &self.0 == GraphValenceTable::default_table() {
            return "ValenceTable.default()".to_owned();
        }

        let entries = (1..=MAX_ATOMIC_NUMBER)
            .filter_map(ChemElement::from_atomic_number)
            .filter_map(|element| {
                self.0.entry(element).map(|entry| {
                    format!(
                        "Element('{}'): ValenceEntry(target_covalences={:?}, aromatic_valences={:?})",
                        element.symbol(),
                        entry.target_covalences,
                        entry.aromatic_valences
                    )
                })
            })
            .collect::<Vec<_>>()
            .join(", ");
        format!("ValenceTable(entries={{{entries}}})")
    }
}

impl ValenceTable {
    pub(crate) fn from_rust(table: &GraphValenceTable) -> Self {
        Self(table.clone())
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for ChemistryModel configuration"
    )]
    pub(crate) fn to_rust(&self) -> GraphValenceTable {
        self.0.clone()
    }
}

/// Valence perception model and its owned model data.
#[pyclass(eq, frozen, from_py_object)]
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ValenceModel {
    /// Atom-typing valence model.
    #[pyo3(constructor = (*, registry))]
    AtomTyping { registry: AtomTypeRegistry },
    /// Counts-based valence model.
    #[pyo3(constructor = (*, table))]
    Counts { table: ValenceTable },
}

#[pymethods]
impl ValenceModel {
    pub(crate) fn __repr__(&self) -> String {
        match self {
            Self::AtomTyping { registry } => {
                format!("ValenceModel.AtomTyping(registry={})", registry.__repr__())
            }
            Self::Counts { table } => {
                format!("ValenceModel.Counts(table={})", table.__repr__())
            }
        }
    }
}

impl ValenceModel {
    pub(crate) fn from_rust(model: &GraphValenceModel) -> Self {
        match model {
            GraphValenceModel::AtomTyping { registry } => Self::AtomTyping {
                registry: AtomTypeRegistry::from_rust(registry.as_ref()),
            },
            GraphValenceModel::Counts { table } => Self::Counts {
                table: ValenceTable::from_rust(table.as_ref()),
            },
        }
    }

    #[allow(
        dead_code,
        reason = "Python-to-Rust conversion API for ChemistryModel configuration"
    )]
    pub(crate) fn to_rust(&self) -> GraphValenceModel {
        match self {
            Self::AtomTyping { registry } => GraphValenceModel::AtomTyping {
                registry: Cow::Owned(registry.to_rust()),
            },
            Self::Counts { table } => GraphValenceModel::Counts {
                table: Cow::Owned(table.to_rust()),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use pyo3::exceptions::PyValueError;
    use rstest::rstest;
    use umol_ast::ast::AtomAst as AstAtomAst;
    use umol_graph::{registry, valence_table};

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

    #[rstest]
    #[case::empty(Vec::new(), Vec::new(), GraphValenceEntry {
        target_covalences: Vec::new(),
        aromatic_valences: Vec::new(),
    })]
    #[case::multi_state(vec![6, 2, 4, 2], vec![4, 2], GraphValenceEntry {
        target_covalences: vec![2, 2, 4, 6],
        aromatic_valences: vec![2, 4],
    })]
    fn test_valence_entry_new(
        #[case] target_covalences: Vec<u8>,
        #[case] aromatic_valences: Vec<u8>,
        #[case] expected: GraphValenceEntry,
    ) {
        assert_eq!(
            ValenceEntry::new(target_covalences, aromatic_valences).0,
            expected
        );
    }

    #[rstest]
    #[case::empty(ValenceEntry::new(Vec::new(), Vec::new()), Vec::new())]
    #[case::multi_state(
        ValenceEntry::new(vec![6, 2, 4], vec![2]),
        vec![2, 4, 6]
    )]
    fn test_valence_entry_target_covalences(
        #[case] entry: ValenceEntry,
        #[case] expected: Vec<u8>,
    ) {
        Python::attach(|py| {
            assert_eq!(
                entry
                    .target_covalences(py)
                    .unwrap()
                    .extract::<Vec<u8>>()
                    .unwrap(),
                expected
            );
        });
    }

    #[rstest]
    #[case::empty(ValenceEntry::new(Vec::new(), Vec::new()), Vec::new())]
    #[case::multi_state(
        ValenceEntry::new(vec![4], vec![4, 2]),
        vec![2, 4]
    )]
    fn test_valence_entry_aromatic_valences(
        #[case] entry: ValenceEntry,
        #[case] expected: Vec<u8>,
    ) {
        Python::attach(|py| {
            assert_eq!(
                entry
                    .aromatic_valences(py)
                    .unwrap()
                    .extract::<Vec<u8>>()
                    .unwrap(),
                expected
            );
        });
    }

    #[rstest]
    #[case::empty(
        ValenceEntry::new(Vec::new(), Vec::new()),
        "ValenceEntry(target_covalences=[], aromatic_valences=[])"
    )]
    #[case::multi_state(
        ValenceEntry::new(vec![6, 2, 4], vec![4, 2]),
        "ValenceEntry(target_covalences=[2, 4, 6], aromatic_valences=[2, 4])"
    )]
    fn test_valence_entry_repr(#[case] entry: ValenceEntry, #[case] expected: &str) {
        assert_eq!(entry.__repr__(), expected);
    }

    #[rstest]
    #[case::empty(GraphValenceEntry {
        target_covalences: Vec::new(),
        aromatic_valences: Vec::new(),
    })]
    #[case::multi_state(GraphValenceEntry {
        target_covalences: vec![2, 4, 6],
        aromatic_valences: vec![2, 4],
    })]
    fn test_valence_entry_from_rust(#[case] entry: GraphValenceEntry) {
        assert_eq!(ValenceEntry::from_rust(&entry).0, entry);
    }

    #[rstest]
    #[case::empty(ValenceEntry::new(Vec::new(), Vec::new()))]
    #[case::multi_state(ValenceEntry::new(vec![6, 2, 4], vec![4, 2]))]
    fn test_valence_entry_to_rust(#[case] entry: ValenceEntry) {
        assert_eq!(entry.to_rust(), entry.0);
    }

    #[rstest]
    #[case::empty(BTreeMap::new(), GraphValenceTable::empty())]
    #[case::multi_state(
        BTreeMap::from([
            (
                Element::from(ChemElement::C),
                ValenceEntry::new(vec![6, 2, 4], vec![3, 2]),
            ),
            (
                Element::from(ChemElement::O),
                ValenceEntry::new(vec![2], Vec::new()),
            ),
        ]),
        GraphValenceTable::from_toml_str(
            "[C]\ntarget_covalences = [2, 4, 6]\naromatic_valences = [2, 3]\n[O]\ntarget_covalences = [2]"
        ).unwrap()
    )]
    fn test_valence_table_new(
        #[case] entries: BTreeMap<Element, ValenceEntry>,
        #[case] expected: GraphValenceTable,
    ) {
        assert_eq!(ValenceTable::new(entries).0, expected);
    }

    #[rstest]
    fn test_valence_table_default() {
        assert_eq!(
            ValenceTable::default().0,
            *GraphValenceTable::default_table()
        );
    }

    #[rstest]
    fn test_valence_table_from_toml() {
        let input = "[C]\ntarget_covalences = [6, 2, 4]\naromatic_valences = [3, 2]";

        assert_eq!(
            ValenceTable::from_toml(input).unwrap().0,
            GraphValenceTable::from_toml_str(input).unwrap()
        );
    }

    #[rstest]
    fn test_valence_table_from_toml_error() {
        Python::attach(|py| {
            let error = ValenceTable::from_toml("[X]\ntarget_covalences = [1]").unwrap_err();

            assert!(error.is_instance_of::<PyValueError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "invalid valence table: unknown element: X"
            );
        });
    }

    #[rstest]
    fn test_valence_table_content_hash() {
        let table = ValenceTable(valence_table![C => [4, 2], O => [2]]);

        assert_eq!(table.content_hash(), table.0.content_hash());
        assert_eq!(table.content_hash_hex(), table.0.content_hash_hex());
        assert_eq!(
            table.content_hash_hex(),
            format!("{:016x}", table.content_hash())
        );
    }

    #[rstest]
    #[case::present(
        ChemElement::C,
        Some(ValenceEntry::new(vec![2, 4], Vec::new()))
    )]
    #[case::missing(ChemElement::N, None)]
    fn test_valence_table_entry(
        #[case] element: ChemElement,
        #[case] expected: Option<ValenceEntry>,
    ) {
        let table = ValenceTable(valence_table![C => [4, 2], O => [2]]);

        assert_eq!(table.entry(Element::from(element)), expected);
    }

    #[rstest]
    #[case::default(ValenceTable::default(), "ValenceTable.default()")]
    #[case::empty(ValenceTable(GraphValenceTable::empty()), "ValenceTable(entries={})")]
    #[case::custom(
        ValenceTable(valence_table![C => [4, 2], O => [2]]),
        "ValenceTable(entries={Element('C'): ValenceEntry(target_covalences=[2, 4], aromatic_valences=[]), Element('O'): ValenceEntry(target_covalences=[2], aromatic_valences=[])})"
    )]
    fn test_valence_table_repr(#[case] table: ValenceTable, #[case] expected: &str) {
        assert_eq!(table.__repr__(), expected);
    }

    #[rstest]
    #[case::borrowed(Cow::Borrowed(GraphValenceTable::default_table()))]
    #[case::owned(Cow::Owned(valence_table![C => [4, 2], O => [2]]))]
    fn test_valence_table_from_rust(#[case] table: Cow<'static, GraphValenceTable>) {
        assert_eq!(ValenceTable::from_rust(table.as_ref()).0, *table);
    }

    #[rstest]
    #[case::default(ValenceTable::default())]
    #[case::custom(ValenceTable(valence_table![C => [4, 2], O => [2]]))]
    fn test_valence_table_to_rust(#[case] table: ValenceTable) {
        assert_eq!(table.to_rust(), table.0);
    }

    #[rstest]
    #[case::atom_typing(
        GraphValenceModel::AtomTyping {
            registry: Cow::Owned(registry!["C#c0#v4", "O#c0#v2"]),
        },
        ValenceModel::AtomTyping {
            registry: AtomTypeRegistry(registry!["C#c0#v4", "O#c0#v2"]),
        }
    )]
    #[case::counts(
        GraphValenceModel::Counts {
            table: Cow::Owned(valence_table![C => [4, 2], O => [2]]),
        },
        ValenceModel::Counts {
            table: ValenceTable(valence_table![C => [4, 2], O => [2]]),
        }
    )]
    fn test_valence_model_from_rust(
        #[case] model: GraphValenceModel,
        #[case] expected: ValenceModel,
    ) {
        assert_eq!(ValenceModel::from_rust(&model), expected);
    }

    #[rstest]
    #[case::atom_typing(
        ValenceModel::AtomTyping {
            registry: AtomTypeRegistry(registry!["C#c0#v4", "O#c0#v2"]),
        },
        GraphValenceModel::AtomTyping {
            registry: Cow::Owned(registry!["C#c0#v4", "O#c0#v2"]),
        }
    )]
    #[case::counts(
        ValenceModel::Counts {
            table: ValenceTable(valence_table![C => [4, 2], O => [2]]),
        },
        GraphValenceModel::Counts {
            table: Cow::Owned(valence_table![C => [4, 2], O => [2]]),
        }
    )]
    fn test_valence_model_to_rust(
        #[case] model: ValenceModel,
        #[case] expected: GraphValenceModel,
    ) {
        let rust = model.to_rust();
        assert_eq!(rust, expected);
        match rust {
            GraphValenceModel::AtomTyping {
                registry: Cow::Owned(_),
            }
            | GraphValenceModel::Counts {
                table: Cow::Owned(_),
            } => {}
            other => panic!("expected owned valence model data, got {other:?}"),
        }
    }

    #[rstest]
    #[case::atom_typing(
        ValenceModel::AtomTyping {
            registry: AtomTypeRegistry(registry!["C#c0#v4"]),
        },
        "ValenceModel.AtomTyping(registry=AtomTypeRegistry.from_atoms([AtomAst.parse(\"C#i=#c0#h0#n0#u0#s#v4#d0#t0#a!#m!#T!\")]))"
    )]
    #[case::counts(
        ValenceModel::Counts {
            table: ValenceTable(valence_table![C => [4, 2]]),
        },
        "ValenceModel.Counts(table=ValenceTable(entries={Element('C'): ValenceEntry(target_covalences=[2, 4], aromatic_valences=[])}))"
    )]
    fn test_valence_model_repr(#[case] model: ValenceModel, #[case] expected: &str) {
        assert_eq!(model.__repr__(), expected);
    }
}
