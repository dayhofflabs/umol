//! Configuration data for GraphIR.

use std::collections::BTreeMap;
use std::fmt::Write;
use std::fs;
use std::path::Path;
use std::sync::LazyLock;

use serde::Deserialize;
use smallvec::SmallVec;
use umol_data::Element;
use xxhash_rust::const_xxh3::xxh3_64;

use crate::graph_ir::atom_type::{AtomTypeQuery, AtomTypeSpec};
use crate::graph_ir::error::ResolutionError;

/// Atom type registry for GraphIR.
///
/// Each spec is stored under both `(element, Some(charge))` and `(element, None)`,
/// enabling O(1) lookup for both charge-specific and element-only queries.
#[derive(Debug, Clone)]
pub struct AtomTypeRegistry {
    atom_types: BTreeMap<(Element, Option<i8>), Vec<AtomTypeSpec>>,
    content_hash: u64,
}

/// Two-level TOML map: element symbol -> charge string -> specs list.
type AtomTypeRegistryToml = BTreeMap<String, BTreeMap<String, Vec<AtomTypeSpec>>>;

impl AtomTypeRegistry {
    pub fn new() -> Self {
        let mut reg = AtomTypeRegistry {
            atom_types: BTreeMap::new(),
            content_hash: 0,
        };
        reg.recompute_hash();
        reg
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

    pub fn content_hash(&self) -> u64 {
        self.content_hash
    }

    pub fn content_hash_hex(&self) -> String {
        format!("{:016x}", self.content_hash)
    }

    fn recompute_hash(&mut self) {
        let mut buf = String::new();
        for ((element, charge), specs) in &self.atom_types {
            let _ = write!(buf, "{},{:?}:", element, charge);
            let mut spec_strs: Vec<String> = specs.iter().map(|s| s.to_string()).collect();
            spec_strs.sort();
            for s in &spec_strs {
                let _ = write!(buf, "{},", s);
            }
            buf.push('\n');
        }
        self.content_hash = xxh3_64(buf.as_bytes());
    }

    pub fn from_toml_str(input: &str) -> Result<Self, ResolutionError> {
        let parsed: AtomTypeRegistryToml = toml::from_str(input)
            .map_err(|e| ResolutionError::InvalidAtomTypeRegistry(e.to_string()))?;
        let mut atom_types: BTreeMap<(Element, Option<i8>), Vec<AtomTypeSpec>> = BTreeMap::new();
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
        let mut reg = AtomTypeRegistry {
            atom_types,
            content_hash: 0,
        };
        reg.recompute_hash();
        Ok(reg)
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
        self.recompute_hash();
    }

    /// Returns true if all charge-specific registry specs satisfy atom invariants.
    /// This is an optional consistency check and is not enforced at load time.
    pub fn verify_specs(&self) -> bool {
        self.atom_types
            .iter()
            .filter(|((_, charge), _)| charge.is_some())
            .all(|(_, specs)| specs.iter().all(|spec| spec.check_invariants().is_ok()))
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
            .filter(|spec| query.matches(spec))
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
/// let reg = registry!["{H+0v1}", "{C+0v4}", "{C+1v3}"];
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
    pub allowed_valences: Vec<i8>,
    pub allowed_aromatic_valences: Vec<u8>,
}

/// Valence table for counts-based validation.
#[derive(Debug, Clone)]
pub struct ValenceTable {
    entries: BTreeMap<Element, ValenceEntry>,
    content_hash: u64,
}

#[derive(Deserialize)]
struct ValenceEntryToml {
    allowed_valences: Vec<i8>,
    #[serde(default)]
    allowed_aromatic_valences: Vec<u8>,
}

impl ValenceTable {
    pub fn empty() -> Self {
        let mut table = ValenceTable {
            entries: BTreeMap::new(),
            content_hash: 0,
        };
        table.recompute_hash();
        table
    }

    pub fn insert(&mut self, element: Element, entry: ValenceEntry) {
        self.entries.insert(element, entry);
        self.recompute_hash();
    }

    pub fn content_hash(&self) -> u64 {
        self.content_hash
    }

    pub fn content_hash_hex(&self) -> String {
        format!("{:016x}", self.content_hash)
    }

    fn recompute_hash(&mut self) {
        let mut buf = String::new();
        for (element, entry) in &self.entries {
            let _ = write!(
                buf,
                "{}:{:?}:{:?}\n",
                element,
                entry.allowed_valences,
                entry.allowed_aromatic_valences
            );
        }
        self.content_hash = xxh3_64(buf.as_bytes());
    }

    pub fn default_table() -> &'static Self {
        &DEFAULT_VALENCE_TABLE
    }

    pub fn from_toml_str(input: &str) -> Result<Self, ResolutionError> {
        let parsed: BTreeMap<String, ValenceEntryToml> = toml::from_str(input)
            .map_err(|e| ResolutionError::InvalidValenceTable(e.to_string()))?;
        let mut entries = BTreeMap::new();
        for (symbol, entry) in parsed {
            let element: Element = symbol.parse().map_err(|_| {
                ResolutionError::InvalidValenceTable(format!("unknown element: {}", symbol))
            })?;
            entries.insert(
                element,
                ValenceEntry {
                    allowed_valences: entry.allowed_valences,
                    allowed_aromatic_valences: entry.allowed_aromatic_valences,
                },
            );
        }
        let mut table = ValenceTable {
            entries,
            content_hash: 0,
        };
        table.recompute_hash();
        Ok(table)
    }

    pub fn entry(&self, element: Element) -> Option<&ValenceEntry> {
        self.entries.get(&element)
    }

    /// Compute implicit hydrogen count using RDKit-style counts logic.
    ///
    /// For charged atoms with specific allowed valences, looks up the
    /// isoelectronic element (atomic_number − charge) following RDKit's
    /// effective-atomic-number convention. Walks `allowed_valences` in order:
    /// `-1` means unconstrained (implicit H = max(0, valence_electrons − charge
    /// − explicit_valence)); otherwise the first allowed value ≥
    /// `explicit_valence` gives `implicit_h = allowed − explicit_valence`.
    ///
    /// Returns `None` when no valid valence state exists.
    pub fn compute_implicit_hydrogens(
        &self,
        element: Element,
        charge: i8,
        explicit_valence: u8,
    ) -> Option<u8> {
        if element == Element::H {
            return Some(0);
        }
        let entry = self.entries.get(&element)?;
        // Element metadata is the canonical source of valence electron counts.
        let num_electrons = (element.valence_electrons() as i16) - (charge as i16);

        let effective_valences = if charge != 0 && !entry.allowed_valences.contains(&-1) {
            let eff_entry = self.entries.get(&element.shift(-charge)?)?;
            &eff_entry.allowed_valences
        } else {
            &entry.allowed_valences
        };

        for &allowed in effective_valences {
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
///     C => [4],
///     N => [3],
///     S => [2, 4, 6],
///     Fe => [-1],
/// };
/// ```
#[macro_export]
macro_rules! valence_table {
    ($($el:ident => [$($v:expr),* $(,)?]),* $(,)?) => {{
        let mut table = $crate::graph_ir::config_data::ValenceTable::empty();
        $(
            table.insert(
                <umol_data::Element as std::str::FromStr>::from_str(stringify!($el))
                    .expect("invalid element symbol in valence_table!"),
                $crate::graph_ir::config_data::ValenceEntry {
                    allowed_valences: vec![$($v),*],
                    allowed_aromatic_valences: vec![],
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

/// Per-element default valence used for `ImplicitHydrogens::Normal`.
#[derive(Debug, Clone)]
pub struct NormalValenceTable {
    neutral: BTreeMap<Element, u8>,
}

impl NormalValenceTable {
    pub fn default_table() -> &'static Self {
        &DEFAULT_NORMAL_VALENCE_TABLE
    }

    pub fn from_toml_str(input: &str) -> Result<Self, ResolutionError> {
        let parsed: BTreeMap<String, u8> = toml::from_str(input)
            .map_err(|e| ResolutionError::InvalidValenceTable(e.to_string()))?;
        let mut neutral = BTreeMap::new();
        for (symbol, valence) in parsed {
            let element: Element = symbol.parse().map_err(|_| {
                ResolutionError::InvalidValenceTable(format!("unknown element: {}", symbol))
            })?;
            neutral.insert(element, valence);
        }
        Ok(Self { neutral })
    }

    pub fn normal_valence_for(&self, element: Element, charge: i8) -> Option<u8> {
        if charge == 0 {
            return self.neutral.get(&element).copied();
        }
        element
            .shift(-charge)
            .and_then(|isoelectronic| self.neutral.get(&isoelectronic).copied())
            .or_else(|| self.neutral.get(&element).copied())
    }
}

static DEFAULT_NORMAL_VALENCE_TABLE: LazyLock<NormalValenceTable> = LazyLock::new(|| {
    NormalValenceTable::from_toml_str(include_str!(
        "../../config/default-normal-valence-table.toml"
    ))
    .expect("built-in default normal valence table must be valid")
});

#[cfg(test)]
mod tests {
    use umol_data::Element;

    use super::*;

    #[test]
    fn test_atom_type_registry_from_toml() {
        let input = r#"
[C]
0 = ["{C+0v4a0m0}"]

[O]
-1 = ["{O-1/3v1a0m0}"]
"#;
        let registry = AtomTypeRegistry::from_toml_str(input).unwrap();
        assert_eq!(
            registry.specs_for_element_and_charge(Element::C, 0).len(),
            1
        );
        assert_eq!(
            registry.specs_for_element_and_charge(Element::O, -1).len(),
            1
        );
        assert_eq!(registry.content_hash_hex(), "cb498fccad4c230a");
    }

    #[test]
    fn test_default_registry() {
        let reg = AtomTypeRegistry::default_registry();
        assert!(!reg.specs_for_element(Element::C).is_empty());
    }

    #[test]
    fn test_registry_macro() {
        let reg = registry!["{C+0v4}", "{C+1^3v3}"];
        assert_eq!(reg.specs_for_element(Element::C).len(), 2);
    }

    #[test]
    fn test_valence_table_from_toml() {
        let table = ValenceTable::from_toml_str(
            r#"
        [H]
        allowed_valences = [1]
        "#,
        )
        .unwrap();
        assert!(table.entry(Element::H).is_some());
        assert_eq!(table.entry(Element::H).unwrap().allowed_valences, [1]);
        assert_eq!(table.content_hash_hex(), "5c81fe410a7e98a0");
    }

    #[test]
    fn test_default_valence_table() {
        let table = ValenceTable::default_table();
        assert!(table.entry(Element::C).is_some());
    }

    #[test]
    fn test_valence_table_allowed() {
        let table = ValenceTable::default_table();
        assert_eq!(table.compute_implicit_hydrogens(Element::S, 0, 2), Some(0));
        assert_eq!(table.compute_implicit_hydrogens(Element::S, 0, 4), Some(0));
        assert_eq!(table.compute_implicit_hydrogens(Element::S, 0, 3), Some(1));
    }

    #[test]
    fn test_valence_table_unconstrained() {
        let table = ValenceTable::default_table();
        // Fe has [-1], valence_electrons=8
        assert_eq!(table.compute_implicit_hydrogens(Element::Fe, 0, 3), Some(5));
        assert_eq!(table.compute_implicit_hydrogens(Element::Fe, 2, 3), Some(3));
    }

    #[test]
    fn test_valence_table_missing() {
        let table = ValenceTable::from_toml_str("[H]\nallowed_valences = [1]\n").unwrap();
        assert_eq!(table.compute_implicit_hydrogens(Element::C, 0, 0), None);
    }

    #[test]
    fn test_default_normal_valence_table() {
        let table = NormalValenceTable::default_table();
        assert_eq!(table.normal_valence_for(Element::C, 0), Some(4));
        assert_eq!(table.normal_valence_for(Element::N, 1), Some(4));
        assert_eq!(table.normal_valence_for(Element::O, -1), Some(1));
    }

    #[test]
    fn test_normal_valence_table_from_toml() {
        let table = NormalValenceTable::from_toml_str(
            r#"
            C = 4
            N = 3
            "#,
        )
        .unwrap();
        assert_eq!(table.normal_valence_for(Element::C, 0), Some(4));
        assert_eq!(table.normal_valence_for(Element::N, 0), Some(3));
    }
}
