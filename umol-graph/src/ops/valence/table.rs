//! Per-element valence data for the Counts valence resolver.
//!
//! `ValenceTable` carries target covalences and aromatic valences per element.

use std::collections::BTreeMap;
use std::fmt::Write;
use std::sync::LazyLock;

use serde::Deserialize;
use umol_chem::element::Element;
use xxhash_rust::const_xxh3::xxh3_64;

use crate::ops::model::ConfigError;

#[derive(Debug, Clone)]
pub struct ValenceEntry {
    /// Lewis/Langmuir saturation targets, sorted smallest to largest. Counts
    /// uses the first entry ≥ localized valence when `#h` is free. Literal `#h`
    /// overrides.
    pub target_covalences: Vec<u8>,
    /// Admissible aromatic valence counts when aromaticity is active. Empty
    /// means the element is not aromatic-capable.
    pub aromatic_valences: Vec<u8>,
}

#[derive(Debug, Clone)]
pub struct ValenceTable {
    entries: BTreeMap<Element, ValenceEntry>,
    content_hash: u64,
}

#[derive(Deserialize)]
struct ValenceEntryToml {
    #[serde(default)]
    target_covalences: Vec<u8>,
    #[serde(default)]
    aromatic_valences: Vec<u8>,
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
        let mut entry = entry;
        entry.target_covalences.sort_unstable();
        entry.aromatic_valences.sort_unstable();
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
                element, entry.target_covalences, entry.aromatic_valences
            );
        }
        self.content_hash = xxh3_64(buf.as_bytes());
    }

    pub fn default_table() -> &'static Self {
        &DEFAULT_VALENCE_TABLE
    }

    pub fn from_toml_str(input: &str) -> Result<Self, ConfigError> {
        let parsed: BTreeMap<String, ValenceEntryToml> =
            toml::from_str(input).map_err(|e| ConfigError::InvalidValenceTable(e.to_string()))?;
        let mut entries = BTreeMap::new();
        for (symbol, entry) in parsed {
            let element: Element = symbol.parse().map_err(|_| {
                ConfigError::InvalidValenceTable(format!("unknown element: {}", symbol))
            })?;
            let mut target_covalences = entry.target_covalences;
            target_covalences.sort_unstable();
            let mut aromatic_valences = entry.aromatic_valences;
            aromatic_valences.sort_unstable();
            entries.insert(
                element,
                ValenceEntry {
                    target_covalences,
                    aromatic_valences,
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
}

/// Defines a `ValenceTable` from element-name keys.
#[macro_export]
macro_rules! valence_table {
    ($($el:ident => [$($v:expr),* $(,)?]),* $(,)?) => {{
        let mut table = $crate::ops::valence::ValenceTable::empty();
        $(
            table.insert(
                <::umol_chem::element::Element as ::std::str::FromStr>::from_str(stringify!($el))
                    .expect("invalid element symbol in valence_table!"),
                $crate::ops::valence::ValenceEntry {
                    target_covalences: vec![$($v),*],
                    aromatic_valences: vec![],
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

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_chem::element::Element;

    use super::*;

    #[fixture]
    fn table() -> ValenceTable {
        ValenceTable::from_toml_str(
            r#"
        [H]
        target_covalences = [1]
        [S]
        target_covalences = [6, 2, 4]
        aromatic_valences = [2]
        "#,
        )
        .unwrap()
    }

    #[rstest]
    #[case::h(Element::H, vec![1], vec![])]
    #[case::s(Element::S, vec![2, 4, 6], vec![2])]
    fn test_valence_table_from_toml(
        #[case] element: Element,
        #[case] expected_targets: Vec<u8>,
        #[case] expected_aromatic: Vec<u8>,
    ) {
        let t = table();
        let entry = t.entry(element).unwrap();
        assert_eq!(entry.target_covalences, expected_targets);
        assert_eq!(entry.aromatic_valences, expected_aromatic);
    }
}
