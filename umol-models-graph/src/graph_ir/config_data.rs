//! Configuration data for GraphIR.

use std::collections::HashMap;
use std::fs;
use std::path::Path;
use std::sync::LazyLock;

use serde::Deserialize;
use smallvec::SmallVec;
use umol_data::Element;

use super::atom_type::{AtomTypeQuery, AtomTypeSpec};
use super::error::ResolutionError;

/// Atom type registry for GraphIR.
///
/// Each spec is stored under both `(element, Some(charge))` and `(element, None)`,
/// enabling O(1) lookup for both charge-specific and element-only queries.
#[derive(Debug, Clone, Default)]
pub struct AtomTypeRegistry {
    atom_types: HashMap<(Element, Option<i8>), Vec<AtomTypeSpec>>,
}

/// Two-level TOML map: element symbol -> charge string -> specs list.
type AtomTypeRegistryToml = HashMap<String, HashMap<String, Vec<AtomTypeSpec>>>;

impl AtomTypeRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn default_registry() -> &'static Self {
        &DEFAULT_ATOM_TYPE_REGISTRY
    }

    pub fn from_specs(specs: impl IntoIterator<Item = AtomTypeSpec>) -> Self {
        let mut reg = Self::new();
        for spec in specs {
            reg.add(spec);
        }
        reg
    }

    pub fn from_toml_str(input: &str) -> Result<Self, ResolutionError> {
        let parsed: AtomTypeRegistryToml = toml::from_str(input)
            .map_err(|e| ResolutionError::InvalidAtomTypeRegistry(e.to_string()))?;
        let mut atom_types: HashMap<(Element, Option<i8>), Vec<AtomTypeSpec>> = HashMap::new();
        for (element_key, charges) in &parsed {
            let element: Element = element_key.parse().map_err(|_| {
                ResolutionError::InvalidAtomTypeRegistry(format!(
                    "unknown element: {}",
                    element_key
                ))
            })?;
            for (charge_key, specs) in charges {
                let charge: i8 = charge_key.parse().map_err(|_| {
                    ResolutionError::InvalidAtomTypeRegistry(format!(
                        "invalid charge '{}' for element {}",
                        charge_key, element_key
                    ))
                })?;
                atom_types
                    .entry((element, Some(charge)))
                    .or_default()
                    .extend(specs.iter().copied());
                atom_types
                    .entry((element, None))
                    .or_default()
                    .extend(specs.iter().copied());
            }
        }
        Ok(AtomTypeRegistry { atom_types })
    }

    pub fn from_toml_file(path: &Path) -> Result<Self, ResolutionError> {
        let input = fs::read_to_string(path)
            .map_err(|e| ResolutionError::InvalidAtomTypeRegistry(e.to_string()))?;
        Self::from_toml_str(&input)
    }

    pub fn add(&mut self, spec: AtomTypeSpec) {
        self.atom_types
            .entry((spec.element(), Some(spec.charge())))
            .or_default()
            .push(spec);
        self.atom_types
            .entry((spec.element(), None))
            .or_default()
            .push(spec);
    }

    pub fn specs_for_element(&self, element: Element) -> &[AtomTypeSpec] {
        self.atom_types
            .get(&(element, None))
            .map_or(&[], |v| v.as_slice())
    }

    pub fn specs_for_element_and_charge(&self, element: Element, charge: i8) -> &[AtomTypeSpec] {
        self.atom_types
            .get(&(element, Some(charge)))
            .map_or(&[], |v| v.as_slice())
    }

    pub fn candidates_for(&self, query: &AtomTypeQuery) -> SmallVec<[AtomTypeSpec; 4]> {
        self.atom_types
            .get(&(query.element, query.charge))
            .into_iter()
            .flatten()
            .filter(|spec| query.matches_spec(spec))
            .copied()
            .collect()
    }
}

/// Public shorthand for defining atom type registries from spec strings.
///
/// Takes a flat, comma-separated list of atom type spec literals.
/// Element and charge keys are derived from each spec automatically.
///
/// ```ignore
/// let reg = registry!["[H+0v1]", "[C+0v4]", "[C+1v3]"];
/// ```
#[macro_export]
macro_rules! registry {
    ($($spec:expr),* $(,)?) => {{
        let mut registry = $crate::graph_ir::config_data::AtomTypeRegistry::new();
        $(
            registry.add($crate::spec!($spec));
        )*
        registry
    }};
}

static DEFAULT_ATOM_TYPE_REGISTRY: LazyLock<AtomTypeRegistry> = LazyLock::new(|| {
    AtomTypeRegistry::from_toml_str(include_str!("../../config/default-registry.toml"))
        .expect("built-in default registry must be valid")
});

/// Per-element valence data for counts-based validation.
#[derive(Debug, Clone)]
pub struct ValenceEntry {
    pub outer_electrons: u8,
    pub allowed_valences: Vec<i8>,
}

/// Valence table for counts-based validation.
#[derive(Debug, Clone)]
pub struct ValenceTable {
    entries: HashMap<Element, ValenceEntry>,
}

#[derive(Deserialize)]
struct ValenceEntryToml {
    outer_electrons: u8,
    allowed_valences: Vec<i8>,
}

impl ValenceTable {
    pub fn empty() -> Self {
        ValenceTable {
            entries: HashMap::new(),
        }
    }

    pub fn insert(&mut self, element: Element, entry: ValenceEntry) {
        self.entries.insert(element, entry);
    }

    pub fn default_table() -> &'static Self {
        &DEFAULT_VALENCE_TABLE
    }

    pub fn from_toml_str(input: &str) -> Result<Self, ResolutionError> {
        let parsed: HashMap<String, ValenceEntryToml> = toml::from_str(input)
            .map_err(|e| ResolutionError::InvalidValenceTable(e.to_string()))?;
        let mut entries = HashMap::new();
        for (symbol, entry) in parsed {
            let element: Element = symbol.parse().map_err(|_| {
                ResolutionError::InvalidValenceTable(format!("unknown element: {}", symbol))
            })?;
            entries.insert(
                element,
                ValenceEntry {
                    outer_electrons: entry.outer_electrons,
                    allowed_valences: entry.allowed_valences,
                },
            );
        }
        Ok(ValenceTable { entries })
    }

    pub fn entry(&self, element: Element) -> Option<&ValenceEntry> {
        self.entries.get(&element)
    }

    /// Compute implicit hydrogen count using RDKit-style counts logic.
    ///
    /// Walks `allowed_valences` in order. `-1` means unconstrained: implicit H
    /// is `max(0, outer_electrons − charge − explicit_valence)`. Otherwise the
    /// first allowed value ≥ `explicit_valence` gives
    /// `implicit_h = allowed − explicit_valence`.
    ///
    /// Returns `None` when no valid valence state exists.
    pub fn compute_implicit_hydrogens(
        &self,
        element: Element,
        charge: i8,
        explicit_valence: u8,
    ) -> Option<u8> {
        let entry = self.entries.get(&element)?;
        let num_electrons = (entry.outer_electrons as i16) - (charge as i16);

        for &allowed in &entry.allowed_valences {
            if allowed == -1 {
                return Some((num_electrons - explicit_valence as i16).max(0) as u8);
            }
            if allowed as u8 >= explicit_valence {
                return Some(allowed as u8 - explicit_valence);
            }
        }
        None
    }
}

/// Public shorthand for defining valence tables.
///
/// ```ignore
/// let table = valence_table! {
///     C => 4, [4],
///     N => 5, [3],
///     S => 6, [2, 4, 6],
///     Fe => 8, [-1],
/// };
/// ```
#[macro_export]
macro_rules! valence_table {
    ($($el:ident => $outer:expr, [$($v:expr),* $(,)?]),* $(,)?) => {{
        let mut table = $crate::graph_ir::config_data::ValenceTable::empty();
        $(
            table.insert(
                <umol_data::Element as std::str::FromStr>::from_str(stringify!($el))
                    .expect("invalid element symbol in valence_table!"),
                $crate::graph_ir::config_data::ValenceEntry {
                    outer_electrons: $outer,
                    allowed_valences: vec![$($v),*],
                },
            );
        )*
        table
    }};
}

static DEFAULT_VALENCE_TABLE: LazyLock<ValenceTable> = LazyLock::new(|| {
    ValenceTable::from_toml_str(include_str!("../../config/default-valence-table.toml"))
        .expect("built-in default valence table must be valid")
});

#[cfg(test)]
mod tests {
    use umol_data::Element;

    use super::*;

    #[test]
    fn atom_type_registry_from_toml() {
        let input = r#"
[C]
0 = ["[C+0v4a0m0]"]

[O]
-1 = ["[O-1/3v1a0m0]"]
"#;
        let registry = AtomTypeRegistry::from_toml_str(input).unwrap();
        assert_eq!(registry.specs_for_element_and_charge(Element::C, 0).len(), 1);
        assert_eq!(registry.specs_for_element_and_charge(Element::O, -1).len(), 1);
    }

    #[test]
    fn default_registry_is_populated() {
        let reg = AtomTypeRegistry::default_registry();
        assert!(!reg.specs_for_element(Element::C).is_empty());
    }

    #[test]
    fn registry_macro_builds() {
        let reg = registry!["[C+0v4]", "[C+1^3v3]"];
        assert_eq!(reg.specs_for_element(Element::C).len(), 2);
    }

    #[test]
    fn default_table_loads() {
        let table = ValenceTable::default_table();
        assert!(table.entry(Element::C).is_some());
    }

    #[test]
    fn carbon_neutral_four_bonds() {
        let table = ValenceTable::default_table();
        assert_eq!(table.compute_implicit_hydrogens(Element::C, 0, 4), Some(0));
    }

    #[test]
    fn carbon_neutral_two_bonds() {
        let table = ValenceTable::default_table();
        assert_eq!(table.compute_implicit_hydrogens(Element::C, 0, 2), Some(2));
    }

    #[test]
    fn carbon_exceeds_all_allowed() {
        let table = ValenceTable::default_table();
        assert_eq!(table.compute_implicit_hydrogens(Element::C, 0, 5), None);
    }

    #[test]
    fn nitrogen_neutral_three_bonds() {
        let table = ValenceTable::default_table();
        assert_eq!(table.compute_implicit_hydrogens(Element::N, 0, 3), Some(0));
    }

    #[test]
    fn sulfur_allowed_list() {
        let table = ValenceTable::default_table();
        assert_eq!(table.compute_implicit_hydrogens(Element::S, 0, 2), Some(0));
        assert_eq!(table.compute_implicit_hydrogens(Element::S, 0, 4), Some(0));
        assert_eq!(table.compute_implicit_hydrogens(Element::S, 0, 3), Some(1));
    }

    #[test]
    fn iron_unconstrained() {
        let table = ValenceTable::default_table();
        // Fe has [-1], outer_electrons=8
        assert_eq!(table.compute_implicit_hydrogens(Element::Fe, 0, 3), Some(5));
        assert_eq!(table.compute_implicit_hydrogens(Element::Fe, 2, 3), Some(3));
    }

    #[test]
    fn missing_element_returns_none() {
        let table =
            ValenceTable::from_toml_str("[H]\nouter_electrons = 1\nallowed_valences = [1]\n")
                .unwrap();
        assert_eq!(table.compute_implicit_hydrogens(Element::C, 0, 0), None);
    }
}
