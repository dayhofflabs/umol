//! Per-element valence data for the Counts valence resolver.
//!
//! `ValenceTable` carries the allowed covalences per element.

use std::collections::BTreeMap;
use std::fmt::Write;
use std::sync::LazyLock;

use serde::Deserialize;
use umol_shared::element::Element;
use xxhash_rust::const_xxh3::xxh3_64;

use crate::ops::model::ConfigError;

#[derive(Debug, Clone)]
pub struct ValenceEntry {
    /// Admissible covalences for neutral atoms.
    /// Covalence counts bonding electrons counting toward shell completions.
    /// (Langmuir covalency), `localized_valence + implicit_H + aromatic_increment`.
    /// Covalences for charged atoms are computed by shifting element by charge.
    pub covalence_set: Vec<u8>,
}

#[derive(Debug, Clone)]
pub struct ValenceTable {
    entries: BTreeMap<Element, ValenceEntry>,
    content_hash: u64,
}

#[derive(Deserialize)]
struct ValenceEntryToml {
    #[serde(default)]
    covalence_set: Vec<u8>,
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
            let _ = writeln!(buf, "{}:{:?}", element, entry.covalence_set);
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
            entries.insert(
                element,
                ValenceEntry {
                    covalence_set: entry.covalence_set,
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

/// Defines a `ValenceTable` from element-name keys. Empty list means
/// "no preferred valence" (transition metals, ionic cores).
///
/// ```ignore
/// let table = valence_table! {
///     C => [4],
///     N => [3],
///     S => [2, 4, 6],
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
                    covalence_set: vec![$($v),*],
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
    use umol_shared::element::Element;

    use super::*;

    #[test]
    fn test_valence_table_from_toml() {
        let table = ValenceTable::from_toml_str(
            r#"
        [H]
        covalence_set = [1]
        "#,
        )
        .unwrap();
        assert_eq!(table.entry(Element::H).unwrap().covalence_set, [1]);
    }

    #[test]
    fn test_valence_table_default_table() {
        let table = ValenceTable::default_table();
        assert!(table.entry(Element::C).is_some());
    }
}
