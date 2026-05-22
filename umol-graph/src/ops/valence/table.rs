//! Per-element valence data for the Counts valence resolver.
//!
//! `ValenceTable` carries the allowed σ-bond-order valences (and aromatic
//! valences) per element. The conventional `Normal`-resolving valences
//! (used to fill `ValueAst::Undetermined`) live in
//! [`crate::ops::valence::NormalImplicitHydrogensTable`].

use std::collections::BTreeMap;
use std::fmt::Write;
use std::sync::LazyLock;

use serde::Deserialize;
use umol_shared::element::Element;
use xxhash_rust::const_xxh3::xxh3_64;

use crate::ops::config::ConfigError;

#[derive(Debug, Clone)]
pub struct ValenceEntry {
    /// Admissible σ-valences for the neutral atom. Empty means "no preferred
    /// σ-valence" (transition metals, post-transition ions): the atom is
    /// accepted at any explicit valence and contributes no implicit
    /// hydrogens. For charged atoms, lookup falls back to the isoelectronic
    /// neutral entry.
    pub allowed_valences: Vec<u8>,
    pub allowed_aromatic_valences: Vec<u8>,
}

#[derive(Debug, Clone)]
pub struct ValenceTable {
    entries: BTreeMap<Element, ValenceEntry>,
    content_hash: u64,
}

#[derive(Deserialize)]
struct ValenceEntryToml {
    #[serde(default)]
    allowed_valences: Vec<u8>,
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
                element, entry.allowed_valences, entry.allowed_aromatic_valences,
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

    /// Compute implicit-hydrogen count from the table's `allowed_valences`.
    ///
    /// For a charged atom, the lookup transparently switches to the
    /// isoelectronic neutral entry (`element.shift(-charge)`), so e.g. Na+
    /// reads as Ne and resolves to 0 implicit H without the parent table
    /// having to special-case the cation.
    ///
    /// Then walks `allowed_valences` and returns `allowed -
    /// explicit_valence` for the first entry that satisfies `explicit ≤
    /// allowed ≤ num_electrons`. The upper bound rejects covalent valences
    /// the atom can't actually reach (Na+ valence 1: needs 1 electron,
    /// has 0).
    ///
    /// An empty `allowed_valences` means the table imposes no σ-valence
    /// preference (transition metals, ionic cores): returns `Some(0)`. A
    /// non-empty list with no entry that fits the constraints returns
    /// `None`.
    pub fn compute_implicit_hydrogens(
        &self,
        element: Element,
        charge: i8,
        explicit_valence: u8,
    ) -> Option<u8> {
        if element == Element::H {
            return Some(0);
        }
        let num_electrons = (element.valence_electrons() as i16) - (charge as i16);

        let effective_valences: &[u8] = if charge != 0 {
            let iso = element.shift(-charge)?;
            &self.entries.get(&iso)?.allowed_valences
        } else {
            &self.entries.get(&element)?.allowed_valences
        };

        if effective_valences.is_empty() {
            return Some(0);
        }
        for &allowed in effective_valences {
            if allowed >= explicit_valence && (allowed as i16) <= num_electrons {
                return Some(allowed - explicit_valence);
            }
        }
        None
    }
}

/// Defines a `ValenceTable` from element-name keys. Empty list means
/// "no preferred σ-valence" (transition metals, ionic cores).
///
/// ```ignore
/// let table = valence_table! {
///     C => [4],
///     N => [3],
///     S => [2, 4, 6],
///     Fe => [],
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
    #[case::fe_neutral(Element::Fe, 0, 3, Some(0))]
    #[case::fe_two_plus(Element::Fe, 2, 3, Some(0))]
    #[case::na_plus(Element::Na, 1, 0, Some(0))]
    #[case::mg_two_plus(Element::Mg, 2, 0, Some(0))]
    fn test_valence_table_compute_implicit_hydrogens_no_preference(
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
}
