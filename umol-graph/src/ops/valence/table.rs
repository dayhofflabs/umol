//! Per-element valence data for the Counts valence resolver.
//!
//! `ValenceTable` carries the allowed covalences (and aromatic valences) per
//! element.

use std::collections::BTreeMap;
use std::fmt::Write;
use std::sync::LazyLock;

use serde::Deserialize;
use umol_shared::element::Element;
use xxhash_rust::const_xxh3::xxh3_64;

use crate::ops::model::ConfigError;

#[derive(Debug, Clone)]
pub struct ValenceEntry {
    /// Admissible **covalences** for the neutral atom — the count of covalent
    /// bonds it forms (Langmuir covalency), `localized_valence + implicit_H +
    /// aromatic_increment`: the electrons gained by sharing, one per bond. This
    /// is neither the realized localized valence `v` (which excludes implicit
    /// H and π participation) nor `total_valence` (which also counts a *donated*
    /// aromatic lone pair). For a neutral atom it is the octet/duet-completion
    /// count (C 4, N 3, O 2, F 1); charged atoms read the isoelectronic neutral
    /// entry via `element.shift(-charge)`. Empty means "no preferred covalence"
    /// (transition metals, ionic cores): the atom is accepted at any explicit
    /// valence and contributes no implicit hydrogens.
    pub covalence_set: Vec<u8>,
    /// Admissible aromatic valences (`av`, the π contribution) when the atom is
    /// aromatic: 1 = standard ring atom, 2 = lone-pair donor, 0 = acceptor.
    /// Non-empty marks the element as aromaticity-capable.
    pub aromatic_valence_set: Vec<u8>,
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
    #[serde(default)]
    aromatic_valence_set: Vec<u8>,
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
                element, entry.covalence_set, entry.aromatic_valence_set,
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
                    covalence_set: entry.covalence_set,
                    aromatic_valence_set: entry.aromatic_valence_set,
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

    /// Compute implicit-hydrogen count from the table's `covalence_set`.
    ///
    /// For a charged atom, the lookup transparently switches to the
    /// isoelectronic neutral entry (`element.shift(-charge)`), so e.g. Na+
    /// reads as Ne and resolves to 0 implicit H without the parent table
    /// having to special-case the cation.
    ///
    /// Then walks `covalence_set` and returns `allowed -
    /// explicit_valence` for the first entry that satisfies `explicit ≤
    /// allowed ≤ num_electrons`. The upper bound rejects covalent valences
    /// the atom can't actually reach (Na+ valence 1: needs 1 electron,
    /// has 0).
    ///
    /// An empty `covalence_set` means the table imposes no valence
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
            &self.entries.get(&iso)?.covalence_set
        } else {
            &self.entries.get(&element)?.covalence_set
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
/// "no preferred valence" (transition metals, ionic cores).
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
                    covalence_set: vec![$($v),*],
                    aromatic_valence_set: vec![],
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
        let table = ValenceTable::from_toml_str("[H]\ncovalence_set = [1]\n").unwrap();
        assert_eq!(table.compute_implicit_hydrogens(Element::C, 0, 0), None);
    }
}
