//! Per-element conventional σ-valences (the `normal_valence` and
//! `aromatic_normal_valence` data) used to resolve
//! `ImplicitHydrogensAst::Normal`. Split out of `ValenceTable` so the
//! atom-typing resolver (which has no `allowed_valences` data) and the
//! counts resolver share the same H-inference table without forcing the
//! atom-typing path to carry an unused valence-table.

use std::collections::BTreeMap;
use std::fmt::Write;
use std::sync::LazyLock;

use serde::Deserialize;
use umol_shared::element::Element;
use xxhash_rust::const_xxh3::xxh3_64;

use crate::ops::config::ConfigError;

#[derive(Debug, Clone)]
pub struct NormalValenceEntry {
    /// Conventional σ-valence used when an atom asserts
    /// `ImplicitHydrogensAst::Normal`. `None` for elements without a
    /// well-defined normal valence (transition metals, ionic cores).
    pub normal_valence: Option<u8>,
    /// Conventional σ-valence used on an aromatic atom: implicit H =
    /// `aromatic_normal_valence - actual σ-valence`. `None` for elements
    /// without a defined aromatic normal valence.
    pub aromatic_normal_valence: Option<u8>,
}

#[derive(Debug, Clone)]
pub struct NormalValenceTable {
    entries: BTreeMap<Element, NormalValenceEntry>,
    content_hash: u64,
}

#[derive(Deserialize)]
struct NormalValenceEntryToml {
    #[serde(default)]
    normal_valence: Option<u8>,
    #[serde(default)]
    aromatic_normal_valence: Option<u8>,
}

impl NormalValenceTable {
    pub fn empty() -> Self {
        let mut table = Self {
            entries: BTreeMap::new(),
            content_hash: 0,
        };
        table.recompute_hash();
        table
    }

    pub fn insert(&mut self, element: Element, entry: NormalValenceEntry) {
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
                element, entry.normal_valence, entry.aromatic_normal_valence,
            );
        }
        self.content_hash = xxh3_64(buf.as_bytes());
    }

    pub fn default_table() -> &'static Self {
        &DEFAULT_NORMAL_VALENCE_TABLE
    }

    pub fn from_toml_str(input: &str) -> Result<Self, ConfigError> {
        let parsed: BTreeMap<String, NormalValenceEntryToml> =
            toml::from_str(input).map_err(|e| ConfigError::InvalidValenceTable(e.to_string()))?;
        let mut entries = BTreeMap::new();
        for (symbol, entry) in parsed {
            let element: Element = symbol.parse().map_err(|_| {
                ConfigError::InvalidValenceTable(format!("unknown element: {}", symbol))
            })?;
            entries.insert(
                element,
                NormalValenceEntry {
                    normal_valence: entry.normal_valence,
                    aromatic_normal_valence: entry.aromatic_normal_valence,
                },
            );
        }
        let mut table = Self {
            entries,
            content_hash: 0,
        };
        table.recompute_hash();
        Ok(table)
    }

    pub fn entry(&self, element: Element) -> Option<&NormalValenceEntry> {
        self.entries.get(&element)
    }

    /// Conventional σ-valence used when an atom asserts
    /// `ImplicitHydrogensAst::Normal`. For charged atoms, falls back to the
    /// isoelectronic neutral element; if neither has `normal_valence` set,
    /// returns `None`.
    pub fn normal_valence_for(&self, element: Element, charge: i8) -> Option<u8> {
        if charge == 0 {
            return self.entries.get(&element)?.normal_valence;
        }
        if let Some(iso) = element.shift(-charge) {
            if let Some(v) = self.entries.get(&iso).and_then(|e| e.normal_valence) {
                return Some(v);
            }
        }
        self.entries.get(&element)?.normal_valence
    }

    /// Conventional σ-valence used to fill `ImplicitHydrogensAst::Normal` on
    /// an aromatic atom: H = `aromatic_normal_valence - actual σ-valence`.
    /// Returns `None` for charged atoms (no charged-aromatic-normal data) or
    /// for elements without an aromatic normal valence set.
    pub fn aromatic_normal_valence_for(&self, element: Element, charge: i8) -> Option<u8> {
        if charge != 0 {
            return None;
        }
        self.entries.get(&element)?.aromatic_normal_valence
    }

    /// Implicit-hydrogen count for a `Normal`-asserting atom. Dispatches on
    /// aromaticity: aromatic atoms use `aromatic_normal_valence_for`, others
    /// use `normal_valence_for`. Returns `None` when the relevant entry is
    /// missing.
    pub fn implicit_hydrogens_for(
        &self,
        element: Element,
        charge: i8,
        explicit_valence: u8,
        is_aromatic: bool,
    ) -> Option<u8> {
        let normal = if is_aromatic {
            self.aromatic_normal_valence_for(element, charge)?
        } else {
            self.normal_valence_for(element, charge)?
        };
        Some(normal.saturating_sub(explicit_valence))
    }
}

static DEFAULT_NORMAL_VALENCE_TABLE: LazyLock<NormalValenceTable> = LazyLock::new(|| {
    NormalValenceTable::from_toml_str(include_str!(
        "../../../config/default-normal-valence-table.toml"
    ))
    .expect("built-in default normal-valence table must be valid")
});

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_shared::element::Element;

    use super::*;

    #[fixture]
    fn default_table() -> &'static NormalValenceTable {
        NormalValenceTable::default_table()
    }

    #[rstest]
    fn test_normal_valence_table_default_table(default_table: &NormalValenceTable) {
        assert_eq!(default_table.normal_valence_for(Element::C, 0), Some(4));
        assert_eq!(
            default_table.aromatic_normal_valence_for(Element::C, 0),
            Some(3)
        );
    }

    #[rstest]
    #[case::single_entry(
        r#"[C]
        normal_valence = 4
        aromatic_normal_valence = 3
        "#,
        Element::C,
        NormalValenceEntry { normal_valence: Some(4), aromatic_normal_valence: Some(3) },
    )]
    #[case::aromatic_only(
        r#"[B]
        aromatic_normal_valence = 2
        "#,
        Element::B,
        NormalValenceEntry { normal_valence: None, aromatic_normal_valence: Some(2) },
    )]
    fn test_normal_valence_table_from_toml_str(
        #[case] input: &str,
        #[case] element: Element,
        #[case] expected: NormalValenceEntry,
    ) {
        let table = NormalValenceTable::from_toml_str(input).unwrap();
        let entry = table.entry(element).unwrap();
        assert_eq!(entry.normal_valence, expected.normal_valence);
        assert_eq!(
            entry.aromatic_normal_valence,
            expected.aromatic_normal_valence
        );
    }

    #[rstest]
    #[case::c_neutral(Element::C, 0, Some(4))]
    #[case::n_plus(Element::N, 1, Some(4))]
    #[case::o_minus(Element::O, -1, Some(1))]
    #[case::unset_element(Element::Fe, 0, None)]
    fn test_normal_valence_table_normal_valence_for(
        default_table: &NormalValenceTable,
        #[case] element: Element,
        #[case] charge: i8,
        #[case] expected: Option<u8>,
    ) {
        assert_eq!(default_table.normal_valence_for(element, charge), expected);
    }

    #[rstest]
    #[case::c_neutral(Element::C, 0, Some(3))]
    #[case::n_neutral(Element::N, 0, Some(2))]
    #[case::no_aromatic_entry(Element::F, 0, None)]
    #[case::charged_returns_none(Element::C, 1, None)]
    #[case::unset_element(Element::Fe, 0, None)]
    fn test_normal_valence_table_aromatic_normal_valence_for(
        default_table: &NormalValenceTable,
        #[case] element: Element,
        #[case] charge: i8,
        #[case] expected: Option<u8>,
    ) {
        assert_eq!(
            default_table.aromatic_normal_valence_for(element, charge),
            expected
        );
    }

    #[rstest]
    #[case::c_neutral_aromatic(Element::C, 0, true, 1, Some(2))]
    #[case::c_neutral_aliphatic(Element::C, 0, false, 1, Some(3))]
    #[case::n_neutral_aromatic(Element::N, 0, true, 2, Some(0))]
    #[case::missing_returns_none(Element::Fe, 0, false, 0, None)]
    fn test_normal_valence_table_implicit_hydrogens_for(
        default_table: &NormalValenceTable,
        #[case] element: Element,
        #[case] charge: i8,
        #[case] is_aromatic: bool,
        #[case] explicit_valence: u8,
        #[case] expected: Option<u8>,
    ) {
        assert_eq!(
            default_table.implicit_hydrogens_for(element, charge, explicit_valence, is_aromatic),
            expected
        );
    }
}
