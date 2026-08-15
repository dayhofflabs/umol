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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValenceEntry {
    /// Lewis/Langmuir saturation targets, sorted smallest to largest. Counts
    /// uses the first entry ≥ localized valence when `#h` is free. Literal `#h`
    /// overrides.
    pub target_covalences: Vec<u8>,
    /// Admissible aromatic valence counts when aromaticity is active. Empty
    /// means the element is not aromatic-capable.
    pub aromatic_valences: Vec<u8>,
    /// Aromatic valence counts consulted only when `aromatic_valences` admits
    /// no candidate for the atom (carbon's zero-π exocyclic-carbonyl reading).
    pub fallback_aromatic_valences: Vec<u8>,
}

#[derive(Debug, Clone)]
pub struct ValenceTable {
    entries: BTreeMap<Element, ValenceEntry>,
    content_hash: u64,
}

impl PartialEq for ValenceTable {
    fn eq(&self, other: &Self) -> bool {
        self.entries == other.entries
    }
}

impl Eq for ValenceTable {}

#[derive(Deserialize)]
struct ValenceEntryToml {
    #[serde(default)]
    target_covalences: Vec<u8>,
    #[serde(default)]
    aromatic_valences: Vec<u8>,
    #[serde(default)]
    fallback_aromatic_valences: Vec<u8>,
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
        entry.fallback_aromatic_valences.sort_unstable();
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
                "{}:{:?}:{:?}:{:?}",
                element,
                entry.target_covalences,
                entry.aromatic_valences,
                entry.fallback_aromatic_valences
            );
        }
        self.content_hash = xxh3_64(buf.as_bytes());
    }

    pub fn default_table() -> &'static Self {
        &DEFAULT_VALENCE_TABLE
    }

    /// The frozen umol SMILES table behind `ValenceModel::smiles()`.
    pub fn smiles_table() -> &'static Self {
        &SMILES_VALENCE_TABLE
    }

    /// The frozen MDL/CTfile table behind `ValenceModel::mdl()`.
    pub fn mdl_table() -> &'static Self {
        &MDL_VALENCE_TABLE
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
            let mut fallback_aromatic_valences = entry.fallback_aromatic_valences;
            fallback_aromatic_valences.sort_unstable();
            entries.insert(
                element,
                ValenceEntry {
                    target_covalences,
                    aromatic_valences,
                    fallback_aromatic_valences,
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
                    fallback_aromatic_valences: vec![],
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

static SMILES_VALENCE_TABLE: LazyLock<ValenceTable> = LazyLock::new(|| {
    ValenceTable::from_toml_str(include_str!("../../../config/smiles-valence-table.toml"))
        .expect("built-in SMILES valence table must be valid")
});

static MDL_VALENCE_TABLE: LazyLock<ValenceTable> = LazyLock::new(|| {
    ValenceTable::from_toml_str(include_str!("../../../config/mdl-valence-table.toml"))
        .expect("built-in MDL valence table must be valid")
});

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_chem::element::Element;

    use super::*;

    #[rstest]
    fn test_valence_entry_eq() {
        assert_eq!(
            ValenceEntry {
                target_covalences: vec![2, 4, 6],
                aromatic_valences: vec![2],
                fallback_aromatic_valences: vec![],
            },
            ValenceEntry {
                target_covalences: vec![2, 4, 6],
                aromatic_valences: vec![2],
                fallback_aromatic_valences: vec![],
            },
        );
    }

    #[rstest]
    #[case::target_covalences(
        ValenceEntry {
            target_covalences: vec![2, 4],
            aromatic_valences: vec![2],
            fallback_aromatic_valences: vec![],
        },
        ValenceEntry {
            target_covalences: vec![2, 4, 6],
            aromatic_valences: vec![2],
            fallback_aromatic_valences: vec![],
        },
    )]
    #[case::aromatic_valences(
        ValenceEntry {
            target_covalences: vec![2, 4, 6],
            aromatic_valences: vec![2],
            fallback_aromatic_valences: vec![],
        },
        ValenceEntry {
            target_covalences: vec![2, 4, 6],
            aromatic_valences: vec![2, 4],
            fallback_aromatic_valences: vec![],
        },
    )]
    #[case::fallback_aromatic_valences(
        ValenceEntry {
            target_covalences: vec![4],
            aromatic_valences: vec![1],
            fallback_aromatic_valences: vec![],
        },
        ValenceEntry {
            target_covalences: vec![4],
            aromatic_valences: vec![1],
            fallback_aromatic_valences: vec![0],
        },
    )]
    fn test_valence_entry_eq_difference(#[case] left: ValenceEntry, #[case] right: ValenceEntry) {
        assert_ne!(left, right);
    }

    #[rstest]
    fn test_valence_table_eq() {
        let mut left = ValenceTable::empty();
        left.insert(
            Element::C,
            ValenceEntry {
                target_covalences: vec![4, 2],
                aromatic_valences: vec![3, 2],
                fallback_aromatic_valences: vec![],
            },
        );
        left.insert(
            Element::O,
            ValenceEntry {
                target_covalences: vec![2],
                aromatic_valences: vec![2],
                fallback_aromatic_valences: vec![],
            },
        );
        let mut right = ValenceTable::empty();
        right.insert(
            Element::O,
            ValenceEntry {
                target_covalences: vec![2],
                aromatic_valences: vec![2],
                fallback_aromatic_valences: vec![],
            },
        );
        right.insert(
            Element::C,
            ValenceEntry {
                target_covalences: vec![2, 4],
                aromatic_valences: vec![2, 3],
                fallback_aromatic_valences: vec![],
            },
        );

        assert_eq!(left, right);
    }

    #[rstest]
    #[case::missing_element(None)]
    #[case::target_covalences(Some(ValenceEntry {
        target_covalences: vec![2, 4],
        aromatic_valences: vec![2],
        fallback_aromatic_valences: vec![],
    }))]
    #[case::aromatic_valences(Some(ValenceEntry {
        target_covalences: vec![2],
        aromatic_valences: vec![2, 4],
        fallback_aromatic_valences: vec![],
    }))]
    fn test_valence_table_eq_difference(#[case] right_oxygen: Option<ValenceEntry>) {
        let mut left = ValenceTable::empty();
        left.insert(
            Element::C,
            ValenceEntry {
                target_covalences: vec![4],
                aromatic_valences: vec![3],
                fallback_aromatic_valences: vec![],
            },
        );
        left.insert(
            Element::O,
            ValenceEntry {
                target_covalences: vec![2],
                aromatic_valences: vec![2],
                fallback_aromatic_valences: vec![],
            },
        );
        let mut right = ValenceTable::empty();
        right.insert(
            Element::C,
            ValenceEntry {
                target_covalences: vec![4],
                aromatic_valences: vec![3],
                fallback_aromatic_valences: vec![],
            },
        );
        if let Some(entry) = right_oxygen {
            right.insert(Element::O, entry);
        }

        assert_ne!(left, right);
    }

    #[rstest]
    fn test_valence_table_eq_metadata() {
        let mut table = ValenceTable::empty();
        table.insert(
            Element::C,
            ValenceEntry {
                target_covalences: vec![4],
                aromatic_valences: vec![3],
                fallback_aromatic_valences: vec![],
            },
        );
        let mut different_hash = table.clone();
        different_hash.content_hash = table.content_hash().wrapping_add(1);

        assert_eq!(table, different_hash);
    }

    #[rstest]
    fn test_valence_table_from_toml() {
        let mut expected = ValenceTable::empty();
        expected.insert(
            Element::H,
            ValenceEntry {
                target_covalences: vec![1],
                aromatic_valences: vec![],
                fallback_aromatic_valences: vec![],
            },
        );
        expected.insert(
            Element::C,
            ValenceEntry {
                target_covalences: vec![4],
                aromatic_valences: vec![1],
                fallback_aromatic_valences: vec![0, 1],
            },
        );
        expected.insert(
            Element::S,
            ValenceEntry {
                target_covalences: vec![2, 4, 6],
                aromatic_valences: vec![2],
                fallback_aromatic_valences: vec![],
            },
        );

        assert_eq!(
            ValenceTable::from_toml_str(
                r#"
                [H]
                target_covalences = [1]
                [C]
                target_covalences = [4]
                aromatic_valences = [1]
                fallback_aromatic_valences = [1, 0]
                [S]
                target_covalences = [6, 2, 4]
                aromatic_valences = [2]
                "#,
            ),
            Ok(expected),
        );
    }
}
