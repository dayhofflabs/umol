//! Atom-type registry: a lookup of canonical atom patterns keyed by element
//! and (optionally) charge. Consumed by the AtomTyping valence resolver.
//!
//! TOML-loaded entries are parsed via `AtomDsl` and raised with the registry
//! raise defaults: concrete inherent fields except isotope, plus the valence-relevant constraints
//! at their zero values (valence, donated/accepted pairs, aromatic valence,
//! multicenter valence); all other constraints stay unconstrained. Stored under
//! both `(element, Some(charge))` and `(element, None)` for the two lookup modes.

use std::collections::BTreeMap;
use std::fmt::Write;
use std::fs;
use std::path::Path;
use std::sync::LazyLock;

use umol_chem::element::Element;
use umol_graph_ir::dsl::{
    AromaticValenceDefault, AtomDefaults, AtomDsl, IsotopeDefault, MulticenterValenceDefault,
    NumDefault,
};
use umol_graph_ir::ir::{AtomForm, ElementForm, IntoIr, IsotopeMassForm, NumForm};
use xxhash_rust::const_xxh3::xxh3_64;

use crate::ops::model::ConfigError;

/// Valence patterns indexed by literal element and charge.
///
/// Every entry has a literal element, a literal charge in the i8 range, and an
/// Undetermined isotope. Registry entries cannot assert isotope information.
/// Other fields and constraints are retained as supplied.
///
/// # Semantic properties
///
/// Construction and insertion preserve entry order and duplicates within each
/// lookup bucket. Checked and asserted producers establish the same invariant.
/// Failed checked insertion preserves both the entries and the content hash.
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

    /// Constructs a registry from patterns satisfying the entry invariant.
    ///
    /// # Panics
    ///
    /// Panics if an entry has a non-literal element or charge, a charge outside
    /// the i8 range, or an isotope other than Undetermined.
    pub fn from_atoms(atoms: impl IntoIterator<Item = AtomForm>) -> Self {
        Self::try_from_atoms(atoms).expect("invalid atom type registry entries")
    }

    /// Constructs a registry from independently supplied patterns.
    ///
    /// # Errors
    ///
    /// Returns ConfigError::InvalidAtomTypeRegistry if an entry has a non-literal
    /// element or charge, a charge outside the i8 range, or an isotope other than
    /// Undetermined.
    pub fn try_from_atoms(atoms: impl IntoIterator<Item = AtomForm>) -> Result<Self, ConfigError> {
        let mut reg = Self::new();
        for atom in atoms {
            reg.try_add(atom)?;
        }
        Ok(reg)
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

    /// The raise defaults for registry entries: concrete inherent fields except isotope, plus
    /// the valence-relevant constraints at their zero values; all other
    /// constraints stay unconstrained.
    pub fn raise_defaults() -> AtomDefaults {
        AtomDefaults {
            isotope: IsotopeDefault::Required,
            valence: NumDefault::Zero,
            donated_pairs: NumDefault::Zero,
            accepted_pairs: NumDefault::Zero,
            aromatic_valence: AromaticValenceDefault::NotAromatic,
            multicenter_valence: MulticenterValenceDefault::NotMulticenter,
            ..AtomDefaults::concrete()
        }
    }

    /// Loads atom patterns from element and charge sections in TOML.
    ///
    /// # Errors
    ///
    /// Returns ConfigError::InvalidAtomTypeRegistry for invalid TOML or atom DSL,
    /// an invalid registry entry, or an entry that disagrees with its section keys.
    pub fn from_toml_str(input: &str) -> Result<Self, ConfigError> {
        let parsed: AtomTypeRegistryToml = toml::from_str(input)
            .map_err(|e| ConfigError::InvalidAtomTypeRegistry(e.to_string()))?;
        let defaults = Self::raise_defaults();
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

    /// Loads a registry from a TOML file.
    ///
    /// # Errors
    ///
    /// Returns ConfigError::InvalidAtomTypeRegistry if reading the file or
    /// loading its contents fails.
    pub fn from_toml_file(path: &Path) -> Result<Self, ConfigError> {
        let input = fs::read_to_string(path)
            .map_err(|e| ConfigError::InvalidAtomTypeRegistry(e.to_string()))?;
        Self::from_toml_str(&input)
    }

    /// Inserts a pattern satisfying the entry invariant.
    ///
    /// # Panics
    ///
    /// Panics before mutation if the entry has a non-literal element or charge,
    /// a charge outside the i8 range, or an isotope other than Undetermined.
    pub fn add(&mut self, atom: AtomForm) {
        self.try_add(atom)
            .expect("invalid atom type registry entry");
    }

    /// Inserts an independently supplied pattern, preserving the registry on error.
    ///
    /// # Errors
    ///
    /// Returns ConfigError::InvalidAtomTypeRegistry if the entry has a non-literal
    /// element or charge, a charge outside the i8 range, or an isotope other than
    /// Undetermined.
    pub fn try_add(&mut self, atom: AtomForm) -> Result<(), ConfigError> {
        let (element, charge) = entry_key(&atom)?;
        self.atom_types
            .entry((element, Some(charge)))
            .or_default()
            .push(atom.clone());
        self.atom_types
            .entry((element, None))
            .or_default()
            .push(atom);
        self.recompute_hash();
        Ok(())
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

fn entry_key(atom: &AtomForm) -> Result<(Element, i8), ConfigError> {
    let ElementForm::Lit(element) = atom.element else {
        return Err(ConfigError::InvalidAtomTypeRegistry(
            "registry entries must have literal elements".to_owned(),
        ));
    };
    let NumForm::Lit(charge) = atom.charge else {
        return Err(ConfigError::InvalidAtomTypeRegistry(
            "registry entries must have literal charges".to_owned(),
        ));
    };
    let charge = i8::try_from(charge).map_err(|_| {
        ConfigError::InvalidAtomTypeRegistry(format!(
            "registry entry charge {charge} is outside -128..=127"
        ))
    })?;
    if !matches!(atom.isotope_mass, IsotopeMassForm::Undetermined) {
        return Err(ConfigError::InvalidAtomTypeRegistry(
            "registry entries must have undetermined isotopes".to_owned(),
        ));
    }
    Ok((element, charge))
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
    let (atom_element, atom_charge) = entry_key(&atom)?;
    if atom_element != element {
        return Err(ConfigError::InvalidAtomTypeRegistry(format!(
            "atom '{}' element {} does not match section element {}",
            source, atom_element, element
        )));
    }
    if atom_charge != charge {
        return Err(ConfigError::InvalidAtomTypeRegistry(format!(
            "atom '{}' charge {} does not match section charge {}",
            source, atom_charge, charge
        )));
    }
    Ok(atom)
}

/// Defines an `AtomTypeRegistry` from a flat list of atom-DSL patterns.
/// Element and charge are derived from each literal.
/// Isotope must be omitted or explicitly Undetermined.
///
/// # Panics
///
/// Panics if parsing fails or a pattern violates the registry entry invariant.
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
            >>::into_ir(dsl, &$crate::ops::valence::AtomTypeRegistry::raise_defaults());
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
    use std::env;
    use std::process;
    use std::thread;

    use rstest::rstest;
    use umol_graph_ir::atom_dsl;

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

        let registry = AtomTypeRegistry::default_registry();
        assert_eq!(registry, &expected);
        for rows in registry.atom_types.values() {
            for row in rows {
                assert_eq!(row.isotope_mass, IsotopeMassForm::Undetermined);
            }
        }
        let rebuilt = AtomTypeRegistry::try_from_atoms(
            registry
                .atom_types
                .iter()
                .filter(|((_, charge), _)| charge.is_none())
                .flat_map(|(_, rows)| rows.iter().cloned()),
        )
        .unwrap();
        assert_eq!(registry, &rebuilt);
        assert_eq!(registry.content_hash(), rebuilt.content_hash());
    }

    #[rstest]
    #[case::empty(vec![])]
    #[case::patterns(vec![atom_dsl!("C#c0#h4#D0"), atom_dsl!("C#c+#h3")])]
    #[case::duplicates(vec![atom_dsl!("C#c0#h4"), atom_dsl!("C#c0#h4")])]
    #[case::lowest_charge(vec![atom_dsl!("C#c-128")])]
    fn test_atom_type_registry_try_from_atoms(#[case] atoms: Vec<AtomForm>) {
        let registry = AtomTypeRegistry::try_from_atoms(atoms.clone()).unwrap();
        assert_eq!(registry.patterns_for_element(Element::C), atoms);
        assert_eq!(registry, AtomTypeRegistry::from_atoms(atoms.clone()));
        for charge in [-128, 0, 1] {
            let expected: Vec<_> = atoms
                .iter()
                .filter(|atom| atom.charge == NumForm::Lit(charge))
                .cloned()
                .collect();
            assert_eq!(
                registry.patterns_for_element_and_charge(Element::C, charge as i8),
                expected
            );
        }
    }

    #[rstest]
    #[case::element(atom_dsl!("*#c0"), "registry entries must have literal elements")]
    #[case::charge(atom_dsl!("C"), "registry entries must have literal charges")]
    #[case::charge_high(atom_dsl!("C#c128"), "registry entry charge 128 is outside -128..=127")]
    #[case::charge_low(atom_dsl!("C#c-129"), "registry entry charge -129 is outside -128..=127")]
    #[case::natural(atom_dsl!("C#i=#c0"), "registry entries must have undetermined isotopes")]
    #[case::mass(atom_dsl!("C#i13#c0"), "registry entries must have undetermined isotopes")]
    #[case::set(atom_dsl!("C#i{12,13}#c0"), "registry entries must have undetermined isotopes")]
    #[case::variable(atom_dsl!("C#i?mass#c0"), "registry entries must have undetermined isotopes")]
    #[case::restricted_variable(atom_dsl!("C#i?mass :: {12,13}#c0"), "registry entries must have undetermined isotopes")]
    #[case::empty_set(AtomForm { isotope_mass: IsotopeMassForm::lit_set([]), ..atom_dsl!("C#c0") }, "registry entries must have undetermined isotopes")]
    fn test_atom_type_registry_try_from_atoms_error(#[case] atom: AtomForm, #[case] message: &str) {
        assert_eq!(
            AtomTypeRegistry::try_from_atoms([atom_dsl!("O#c0"), atom]),
            Err(ConfigError::InvalidAtomTypeRegistry(message.to_owned())),
        );
    }

    #[rstest]
    #[case::natural(atom_dsl!("C#i=#c0"))]
    #[case::mass(atom_dsl!("C#i13#c0"))]
    #[case::set(atom_dsl!("C#i{12,13}#c0"))]
    #[case::variable(atom_dsl!("C#i?mass#c0"))]
    #[case::element(atom_dsl!("*#c0"))]
    #[case::charge(atom_dsl!("C"))]
    #[should_panic(expected = "invalid atom type registry entries")]
    fn test_atom_type_registry_from_atoms_error(#[case] atom: AtomForm) {
        AtomTypeRegistry::from_atoms([atom]);
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
            isotope: IsotopeDefault::Required,
            valence: NumDefault::Zero,
            donated_pairs: NumDefault::Zero,
            accepted_pairs: NumDefault::Zero,
            aromatic_valence: AromaticValenceDefault::NotAromatic,
            multicenter_valence: MulticenterValenceDefault::NotMulticenter,
            ..AtomDefaults::concrete()
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
    #[case::charge_range(
        "[C]\n0 = [\"C#c256\"]",
        ConfigError::InvalidAtomTypeRegistry("registry entry charge 256 is outside -128..=127".to_owned()),
    )]
    #[case::nonliteral_element(
        "[C]\n0 = [\"{C,N}#c0\"]",
        ConfigError::InvalidAtomTypeRegistry("registry entries must have literal elements".to_owned()),
    )]
    #[case::nonliteral_charge(
        "[C]\n0 = [\"C#c?charge\"]",
        ConfigError::InvalidAtomTypeRegistry("registry entries must have literal charges".to_owned()),
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
    #[case::omitted("[C]\n0 = [\"C#v4\"]", Ok(registry!["C#v4"]))]
    #[case::undetermined("[C]\n0 = [\"C#i*#v4\"]", Ok(registry!["C#v4"]))]
    #[case::natural("[C]\n0 = [\"C#i=\"]", Err(ConfigError::InvalidAtomTypeRegistry("registry entries must have undetermined isotopes".to_owned())))]
    #[case::mass("[C]\n0 = [\"C#i13\"]", Err(ConfigError::InvalidAtomTypeRegistry("registry entries must have undetermined isotopes".to_owned())))]
    #[case::set("[C]\n0 = [\"C#i{12,13}\"]", Err(ConfigError::InvalidAtomTypeRegistry("registry entries must have undetermined isotopes".to_owned())))]
    #[case::variable("[C]\n0 = [\"C#i?mass\"]", Err(ConfigError::InvalidAtomTypeRegistry("registry entries must have undetermined isotopes".to_owned())))]
    #[case::restricted_variable("[C]\n0 = [\"C#i?mass :: {12,13}\"]", Err(ConfigError::InvalidAtomTypeRegistry("registry entries must have undetermined isotopes".to_owned())))]
    fn test_atom_type_registry_from_toml_file(
        #[case] input: &str,
        #[case] expected: Result<AtomTypeRegistry, ConfigError>,
    ) {
        let path = env::temp_dir().join(format!(
            "umol-registry-{}-{:?}.toml",
            process::id(),
            thread::current().id()
        ));
        fs::write(&path, input).unwrap();
        let result = AtomTypeRegistry::from_toml_file(&path);
        fs::remove_file(path).unwrap();
        assert_eq!(result, expected);
        assert_eq!(AtomTypeRegistry::from_toml_str(input), expected);
    }

    #[rstest]
    #[case::new_element("O#h2#D0", "[C]\n0 = [\"C#h4\"]\n[O]\n0 = [\"O#h2#D0\"]")]
    #[case::new_charge("C#c+#h3", "[C]\n0 = [\"C#h4\"]\n1 = [\"C#c+#h3\"]")]
    #[case::duplicate("C#h4", "[C]\n0 = [\"C#h4\", \"C#h4\"]")]
    #[case::highest_charge("C#c127", "[C]\n0 = [\"C#h4\"]\n127 = [\"C#c127\"]")]
    fn test_atom_type_registry_try_add(#[case] source: &str, #[case] expected: &str) {
        let mut registry = AtomTypeRegistry::from_toml_str("[C]\n0 = [\"C#h4\"]").unwrap();
        let atom = source
            .parse::<AtomDsl>()
            .unwrap()
            .into_ir(&AtomTypeRegistry::raise_defaults());
        assert_eq!(registry.try_add(atom), Ok(()));
        let expected = AtomTypeRegistry::from_toml_str(expected).unwrap();
        assert_eq!(registry, expected);
        assert_eq!(registry.content_hash(), expected.content_hash());
    }

    #[rstest]
    #[case::element(atom_dsl!("*#c0"), "registry entries must have literal elements")]
    #[case::charge(atom_dsl!("C"), "registry entries must have literal charges")]
    #[case::charge_high(atom_dsl!("C#c128"), "registry entry charge 128 is outside -128..=127")]
    #[case::charge_low(atom_dsl!("C#c-129"), "registry entry charge -129 is outside -128..=127")]
    #[case::natural(atom_dsl!("C#i=#c0"), "registry entries must have undetermined isotopes")]
    #[case::mass(atom_dsl!("C#i13#c0"), "registry entries must have undetermined isotopes")]
    #[case::set(atom_dsl!("C#i{12,13}#c0"), "registry entries must have undetermined isotopes")]
    #[case::variable(atom_dsl!("C#i?mass#c0"), "registry entries must have undetermined isotopes")]
    #[case::restricted_variable(atom_dsl!("C#i?mass :: {12,13}#c0"), "registry entries must have undetermined isotopes")]
    #[case::empty_set(AtomForm { isotope_mass: IsotopeMassForm::lit_set([]), ..atom_dsl!("C#c0") }, "registry entries must have undetermined isotopes")]
    fn test_atom_type_registry_try_add_error(#[case] atom: AtomForm, #[case] message: &str) {
        let mut registry = registry!["C#v4", "O#v2"];
        let original = registry.clone();
        assert_eq!(
            registry.try_add(atom),
            Err(ConfigError::InvalidAtomTypeRegistry(message.to_owned()))
        );
        assert_eq!(registry, original);
        assert_eq!(registry.content_hash(), original.content_hash());
    }

    #[rstest]
    #[case::natural(atom_dsl!("C#i=#c0"))]
    #[case::mass(atom_dsl!("C#i13#c0"))]
    #[case::set(atom_dsl!("C#i{12,13}#c0"))]
    #[case::variable(atom_dsl!("C#i?mass#c0"))]
    #[case::element(atom_dsl!("*#c0"))]
    #[case::charge(atom_dsl!("C"))]
    #[should_panic(expected = "invalid atom type registry entry")]
    fn test_atom_type_registry_add_error(#[case] atom: AtomForm) {
        AtomTypeRegistry::new().add(atom);
    }

    #[rstest]
    fn test_registry_macro() {
        let expected = AtomTypeRegistry::from_atoms([
            atom_dsl!("C#c0#h0#n0#u0#s#v4#d0#t0#a!#m!"),
            atom_dsl!("C#c+#h3#n0#u0#s#v0#d0#t0#a!#m!"),
        ]);

        assert_eq!(registry!["C#c0#v4", "C#c+#h3"], expected);
    }

    #[rstest]
    #[case::natural("C#i=")]
    #[case::mass("C#i13")]
    #[case::set("C#i{12,13}")]
    #[case::variable("C#i?mass")]
    #[case::restricted_variable("C#i?mass :: {12,13}")]
    #[should_panic(expected = "registry entries must have undetermined isotopes")]
    fn test_registry_macro_error(#[case] source: &str) {
        registry![source];
    }
}
