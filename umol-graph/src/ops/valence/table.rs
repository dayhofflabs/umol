//! Per-element valence data for the Counts valence resolver.
//!
//! `ValenceTable` carries the allowed σ-bond-order valences (and aromatic
//! valences) per element; `NormalValenceTable` carries the per-element default
//! used to fill `ImplicitHydrogensAst::Normal` during raise.

use std::collections::BTreeMap;
use std::fmt::Write;
use std::sync::LazyLock;

use serde::Deserialize;
use umol_shared::element::Element;
use xxhash_rust::const_xxh3::xxh3_64;

use crate::ops::config::ConfigError;

#[derive(Debug, Clone)]
pub struct ValenceEntry {
    pub allowed_valences: Vec<i8>,
    pub allowed_aromatic_valences: Vec<u8>,
}

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
            let _ = writeln!(
                buf,
                "{}:{:?}:{:?}",
                element, entry.allowed_valences, entry.allowed_aromatic_valences
            );
        }
        self.content_hash = xxh3_64(buf.as_bytes());
    }

    pub fn default_table() -> &'static Self {
        &DEFAULT_VALENCE_TABLE
    }

    pub fn from_toml_str(input: &str) -> Result<Self, ConfigError> {
        let parsed: BTreeMap<String, ValenceEntryToml> = toml::from_str(input)
            .map_err(|e| ConfigError::InvalidValenceTable(e.to_string()))?;
        let mut entries = BTreeMap::new();
        for (symbol, entry) in parsed {
            let element: Element = symbol.parse().map_err(|_| {
                ConfigError::InvalidValenceTable(format!("unknown element: {}", symbol))
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

    /// Compute implicit-hydrogen count using RDKit-style counts logic.
    ///
    /// Walks `allowed_valences` for the element (or its isoelectronic
    /// counterpart if charged and the entry doesn't carry `-1`), returning
    /// `allowed - explicit_valence` for the first entry that fits. `-1` means
    /// unconstrained: `max(0, valence_electrons - charge - explicit)`.
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

/// Defines a `ValenceTable` from element-name keys.
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
        let mut table = $crate::ops::valence::ValenceTable::empty();
        $(
            table.insert(
                <::umol_shared::element::Element as ::std::str::FromStr>::from_str(stringify!($el))
                    .expect("invalid element symbol in valence_table!"),
                $crate::ops::valence::ValenceEntry {
                    allowed_valences: vec![$($v),*],
                    allowed_aromatic_valences: vec![],
                },
            );
        )*
        table
    }};
}

static DEFAULT_VALENCE_TABLE: LazyLock<ValenceTable> = LazyLock::new(|| {
    ValenceTable::from_toml_str(include_str!("../../../config/default-valence-table.toml"))
        .expect("built-in default valence table must be valid")
});

#[derive(Debug, Clone)]
pub struct NormalValenceTable {
    neutral: BTreeMap<Element, u8>,
}

impl NormalValenceTable {
    pub fn default_table() -> &'static Self {
        &DEFAULT_NORMAL_VALENCE_TABLE
    }

    pub fn from_toml_str(input: &str) -> Result<Self, ConfigError> {
        let parsed: BTreeMap<String, u8> = toml::from_str(input)
            .map_err(|e| ConfigError::InvalidValenceTable(e.to_string()))?;
        let mut neutral = BTreeMap::new();
        for (symbol, valence) in parsed {
            let element: Element = symbol.parse().map_err(|_| {
                ConfigError::InvalidValenceTable(format!("unknown element: {}", symbol))
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
        "../../../config/default-normal-valence-table.toml"
    ))
    .expect("built-in default normal valence table must be valid")
});

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_shared::element::Element;

    use super::*;

    #[test]
    fn test_valence_table_from_toml() {
        let table = ValenceTable::from_toml_str(
            r#"
        [H]
        allowed_valences = [1]
        "#,
        )
        .unwrap();
        assert_eq!(table.entry(Element::H).unwrap().allowed_valences, [1]);
        assert_eq!(table.content_hash_hex(), "5c81fe410a7e98a0");
    }

    #[test]
    fn test_valence_table_default_table() {
        let table = ValenceTable::default_table();
        assert!(table.entry(Element::C).is_some());
    }

    #[rstest]
    #[case::s_two(Element::S, 0, 2, Some(0))]
    #[case::s_three(Element::S, 0, 3, Some(1))]
    #[case::s_four(Element::S, 0, 4, Some(0))]
    fn test_valence_table_compute_implicit_hydrogens_allowed(
        #[case] element: Element,
        #[case] charge: i8,
        #[case] explicit: u8,
        #[case] expected: Option<u8>,
    ) {
        let table = ValenceTable::default_table();
        assert_eq!(
            table.compute_implicit_hydrogens(element, charge, explicit),
            expected
        );
    }

    #[rstest]
    #[case::fe_neutral(Element::Fe, 0, 3, Some(5))]
    #[case::fe_two_plus(Element::Fe, 2, 3, Some(3))]
    fn test_valence_table_compute_implicit_hydrogens_unconstrained(
        #[case] element: Element,
        #[case] charge: i8,
        #[case] explicit: u8,
        #[case] expected: Option<u8>,
    ) {
        let table = ValenceTable::default_table();
        assert_eq!(
            table.compute_implicit_hydrogens(element, charge, explicit),
            expected
        );
    }

    #[test]
    fn test_valence_table_compute_implicit_hydrogens_missing() {
        let table = ValenceTable::from_toml_str("[H]\nallowed_valences = [1]\n").unwrap();
        assert_eq!(table.compute_implicit_hydrogens(Element::C, 0, 0), None);
    }

    #[rstest]
    #[case::c_neutral(Element::C, 0, Some(4))]
    #[case::n_plus(Element::N, 1, Some(4))]
    #[case::o_minus(Element::O, -1, Some(1))]
    fn test_normal_valence_table_default_table(
        #[case] element: Element,
        #[case] charge: i8,
        #[case] expected: Option<u8>,
    ) {
        let table = NormalValenceTable::default_table();
        assert_eq!(table.normal_valence_for(element, charge), expected);
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
