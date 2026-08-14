//! Atom-type registry: a lookup of canonical atom patterns keyed by element
//! and (optionally) charge. Consumed by the AtomTyping valence resolver.
//!
//! TOML-loaded entries are parsed via `AtomDsl` and raised to ground `AtomForm`
//! values with `ground()` defaults plus the valence-relevant constraints zeroed
//! (valence, donated/accepted pairs, aromatic valence, multicenter valence); all
//! other constraints stay unconstrained. Stored under both `(element, Some(charge))`
//! and `(element, None)` for the two lookup modes.

use std::collections::BTreeMap;
use std::fmt::Write;
use std::fs;
use std::path::Path;
use std::sync::LazyLock;

use umol_chem::element::Element;
use umol_graph_ir::dsl::{
    AromaticValenceDefault, AtomDefaults, AtomDsl, MulticenterValenceDefault, NumDefault,
};
use umol_graph_ir::ir::{AtomForm, ElementForm, IntoIr, NumForm};
use xxhash_rust::const_xxh3::xxh3_64;

use crate::ops::model::ConfigError;

#[derive(Debug, Clone)]
pub struct AtomTypeRegistry {
    atom_types: BTreeMap<(Element, Option<i8>), Vec<AtomForm>>,
    content_hash: u64,
}

impl PartialEq for AtomTypeRegistry {
    fn eq(&self, other: &Self) -> bool {
        self.atom_types == other.atom_types
    }
}

impl Eq for AtomTypeRegistry {}

type AtomTypeRegistryToml = BTreeMap<String, BTreeMap<String, Vec<String>>>;

impl Default for AtomTypeRegistry {
    fn default() -> Self {
        Self::new()
    }
}

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

    pub fn from_atoms(atoms: impl IntoIterator<Item = AtomForm>) -> Self {
        let mut reg = Self::new();
        for atom in atoms {
            reg.add(atom);
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
        for ((element, charge), atoms) in &self.atom_types {
            let _ = write!(buf, "{},{:?}:", element, charge);
            let mut atom_strs: Vec<String> = atoms.iter().map(|a| format!("{:?}", a)).collect();
            atom_strs.sort();
            for s in &atom_strs {
                let _ = write!(buf, "{},", s);
            }
            buf.push('\n');
        }
        self.content_hash = xxh3_64(buf.as_bytes());
    }

    pub fn from_toml_str(input: &str) -> Result<Self, ConfigError> {
        let parsed: AtomTypeRegistryToml = toml::from_str(input)
            .map_err(|e| ConfigError::InvalidAtomTypeRegistry(e.to_string()))?;
        // Zero the valance constraints.
        let defaults = AtomDefaults {
            valence: NumDefault::Zero,
            donated_pairs: NumDefault::Zero,
            accepted_pairs: NumDefault::Zero,
            aromatic_valence: AromaticValenceDefault::NotAromatic,
            multicenter_valence: MulticenterValenceDefault::NotMulticenter,
            ..AtomDefaults::ground()
        };
        let mut atom_types: BTreeMap<(Element, Option<i8>), Vec<AtomForm>> = BTreeMap::new();
        for (element_key, charges) in &parsed {
            let element: Element = element_key.parse().map_err(|_| {
                ConfigError::InvalidAtomTypeRegistry(format!("unknown element: {}", element_key))
            })?;
            for (charge_key, sources) in charges {
                let charge: i8 = charge_key.parse().map_err(|_| {
                    ConfigError::InvalidAtomTypeRegistry(format!(
                        "invalid charge '{}' for element {}",
                        charge_key, element_key
                    ))
                })?;
                let atoms: Vec<AtomForm> = sources
                    .iter()
                    .map(|source| parse_entry(source, &defaults, element, charge))
                    .collect::<Result<_, _>>()?;
                atom_types
                    .entry((element, Some(charge)))
                    .or_default()
                    .extend(atoms.iter().cloned());
                atom_types.entry((element, None)).or_default().extend(atoms);
            }
        }
        let mut reg = AtomTypeRegistry {
            atom_types,
            content_hash: 0,
        };
        reg.recompute_hash();
        Ok(reg)
    }

    pub fn from_toml_file(path: &Path) -> Result<Self, ConfigError> {
        let input = fs::read_to_string(path)
            .map_err(|e| ConfigError::InvalidAtomTypeRegistry(e.to_string()))?;
        Self::from_toml_str(&input)
    }

    pub fn add(&mut self, atom: AtomForm) {
        let element = match &atom.element {
            ElementForm::Lit(e) => *e,
            _ => panic!("registry entries must have literal elements"),
        };
        let charge = match &atom.charge {
            NumForm::Lit(c) => *c as i8,
            _ => panic!("registry entries must have literal charges"),
        };
        self.atom_types
            .entry((element, Some(charge)))
            .or_default()
            .push(atom.clone());
        self.atom_types
            .entry((element, None))
            .or_default()
            .push(atom);
        self.recompute_hash();
    }

    pub fn patterns_for_element(&self, element: Element) -> &[AtomForm] {
        self.atom_types
            .get(&(element, None))
            .map_or(&[], |v| v.as_slice())
    }

    pub fn patterns_for_element_and_charge(&self, element: Element, charge: i8) -> &[AtomForm] {
        self.atom_types
            .get(&(element, Some(charge)))
            .map_or(&[], |v| v.as_slice())
    }

    pub fn lookup(&self, element: Element, charge: Option<i8>) -> &[AtomForm] {
        self.atom_types
            .get(&(element, charge))
            .map_or(&[], |v| v.as_slice())
    }
}

fn parse_entry(
    source: &str,
    defaults: &AtomDefaults,
    element: Element,
    charge: i8,
) -> Result<AtomForm, ConfigError> {
    let dsl: AtomDsl = source
        .parse()
        .map_err(|e| ConfigError::InvalidAtomTypeRegistry(format!("{}: {}", source, e)))?;
    let atom: AtomForm = dsl.into_ir(defaults);
    let &ElementForm::Lit(atom_element) = &atom.element else {
        return Err(ConfigError::InvalidAtomTypeRegistry(format!(
            "atom '{}' has non-literal element",
            source
        )));
    };
    if atom_element != element {
        return Err(ConfigError::InvalidAtomTypeRegistry(format!(
            "atom '{}' element {} does not match section element {}",
            source, atom_element, element
        )));
    }
    let &NumForm::Lit(atom_charge) = &atom.charge else {
        return Err(ConfigError::InvalidAtomTypeRegistry(format!(
            "atom '{}' has non-literal charge",
            source
        )));
    };
    let atom_charge = atom_charge as i8;
    if atom_charge != charge {
        return Err(ConfigError::InvalidAtomTypeRegistry(format!(
            "atom '{}' charge {} does not match section charge {}",
            source, atom_charge, charge
        )));
    }
    Ok(atom)
}

/// Defines an `AtomTypeRegistry` from a flat list of ground atom-DSL literals.
/// Element and charge are derived from each literal.
///
/// ```ignore
/// let reg = registry!["H#v", "C#v4", "C#c+#v3"];
/// ```
#[macro_export]
macro_rules! registry {
    ($($source:expr),* $(,)?) => {{
        let mut registry = $crate::ops::valence::AtomTypeRegistry::new();
        $(
            let dsl: ::umol_graph_ir::dsl::AtomDsl = $source
                .parse()
                .expect("invalid atom DSL");
            let atom: ::umol_graph_ir::ir::AtomForm = <_ as ::umol_graph_ir::ir::IntoIr<
                ::umol_graph_ir::ir::AtomForm,
            >>::into_ir(dsl, &::umol_graph_ir::dsl::AtomDefaults::zeroed());
            registry.add(atom);
        )*
        registry
    }};
}

static DEFAULT_ATOM_TYPE_REGISTRY: LazyLock<AtomTypeRegistry> = LazyLock::new(|| {
    AtomTypeRegistry::from_toml_str(include_str!("../../../config/default-registry.toml"))
        .expect("built-in default registry must be valid")
});

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    fn test_atom_type_registry_eq() {
        assert_eq!(
            registry!["C#c0#v4", "O#c0#v2"],
            registry!["C#c0#v4", "O#c0#v2"],
        );
    }

    #[rstest]
    #[case::pattern(registry!["C#c0#v4"], registry!["C#c0#v3"])]
    #[case::element_bucket(registry!["C#c0#v4"], registry!["N#c0#v4"])]
    #[case::charge_bucket(registry!["C#c0#v4"], registry!["C#c+#v4"])]
    fn test_atom_type_registry_eq_difference(
        #[case] left: AtomTypeRegistry,
        #[case] right: AtomTypeRegistry,
    ) {
        assert_ne!(left, right);
    }

    #[rstest]
    fn test_atom_type_registry_eq_metadata() {
        let registry = registry!["C#c0#v4"];
        let mut different_hash = registry.clone();
        different_hash.content_hash = registry.content_hash().wrapping_add(1);

        assert_eq!(registry, different_hash);
    }

    #[rstest]
    fn test_atom_type_registry_default_registry() {
        let expected =
            AtomTypeRegistry::from_toml_str(include_str!("../../../config/default-registry.toml"))
                .unwrap();

        assert_eq!(AtomTypeRegistry::default_registry(), &expected);
    }

    #[rstest]
    fn test_atom_type_registry_from_toml() {
        let input = r#"
[C]
0 = ["C#c0#v4#a0"]

[O]
-1 = ["O#c-#n3#v#a0"]
"#;
        let defaults = AtomDefaults {
            valence: NumDefault::Zero,
            donated_pairs: NumDefault::Zero,
            accepted_pairs: NumDefault::Zero,
            aromatic_valence: AromaticValenceDefault::NotAromatic,
            multicenter_valence: MulticenterValenceDefault::NotMulticenter,
            ..AtomDefaults::ground()
        };
        let expected = AtomTypeRegistry::from_atoms(
            ["C#c0#v4#a0", "O#c-#n3#v#a0"]
                .map(|source| source.parse::<AtomDsl>().unwrap().into_ir(&defaults)),
        );

        assert_eq!(AtomTypeRegistry::from_toml_str(input), Ok(expected));
    }

    #[rstest]
    #[case::wrong_element(
        "[C]\n0 = [\"O#c0\"]",
        ConfigError::InvalidAtomTypeRegistry(
            "atom 'O#c0' element O does not match section element C".to_owned(),
        ),
    )]
    #[case::wrong_charge(
        "[C]\n0 = [\"C#c+\"]",
        ConfigError::InvalidAtomTypeRegistry(
            "atom 'C#c+' charge 1 does not match section charge 0".to_owned(),
        ),
    )]
    fn test_atom_type_registry_from_toml_error(#[case] input: &str, #[case] expected: ConfigError) {
        assert_eq!(AtomTypeRegistry::from_toml_str(input), Err(expected));
    }

    #[rstest]
    fn test_registry_macro() {
        let defaults = AtomDefaults::zeroed();
        let expected = AtomTypeRegistry::from_atoms(
            ["C#c0#v4", "C#c+#h3"]
                .map(|source| source.parse::<AtomDsl>().unwrap().into_ir(&defaults)),
        );

        assert_eq!(registry!["C#c0#v4", "C#c+#h3"], expected);
    }
}
