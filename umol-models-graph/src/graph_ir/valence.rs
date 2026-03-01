//! Valence resolution strategies and atom typing support for GraphIR.

use std::collections::HashMap;
use std::fmt::{self, Display};
use std::fs;
use std::path::Path;
use std::str::FromStr;
use std::sync::LazyLock;

use serde::de::{self, Deserializer};
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;
use umol_data::{Element, SpinMultiplicity};

use super::error::ResolutionError;
use super::molecule::{AtomIndex, MoleculeBuilder};

/// Atom typing specification used by the atom-typing valence strategy.
///
/// String notation:
/// - `[El...]` where `El` is an element symbol.
/// - tokens are optional and can appear in any order:
///   - `+n` / `-n` charge (default 0, bare `+`/`-` means 1)
///   - `/n` lone pairs
///   - `^n` unpaired electrons
///   - `*n` multiplicity (default `unpaired + 1`)
///   - `Hn` hydrogens
///   - `vn` valence
///   - `>n` donated pairs
///   - `<n` accepted pairs
///   - `an` aromatic valence
///   - `mn` multicenter valence
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AtomTypeSpec {
    element: Element,
    charge: i8,
    hydrogens: u8,
    lone_pairs: u8,
    unpaired_electrons: u8,
    multiplicity: SpinMultiplicity,
    valence: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    aromatic_valence: u8,
    multicenter_valence: u8,
}

impl AtomTypeSpec {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        element: Element,
        charge: i8,
        hydrogens: u8,
        lone_pairs: u8,
        unpaired_electrons: u8,
        multiplicity: SpinMultiplicity,
        valence: u8,
        donated_pairs: u8,
        accepted_pairs: u8,
        aromatic_valence: u8,
        multicenter_valence: u8,
    ) -> Result<Self, ResolutionError> {
        let expected = unpaired_electrons + 1;
        if multiplicity.multiplicity() > expected {
            return Err(ResolutionError::InvalidAtomSpec(format!(
                "multiplicity {} exceeds unpaired_electrons+1 ({})",
                multiplicity.multiplicity(),
                expected
            )));
        }
        Ok(Self {
            element,
            charge,
            hydrogens,
            lone_pairs,
            unpaired_electrons,
            multiplicity,
            valence,
            donated_pairs,
            accepted_pairs,
            aromatic_valence,
            multicenter_valence,
        })
    }

    pub fn element(&self) -> Element {
        self.element
    }

    pub fn charge(&self) -> i8 {
        self.charge
    }

    pub fn hydrogens(&self) -> u8 {
        self.hydrogens
    }

    pub fn lone_pairs(&self) -> u8 {
        self.lone_pairs
    }

    pub fn unpaired_electrons(&self) -> u8 {
        self.unpaired_electrons
    }

    pub fn multiplicity(&self) -> SpinMultiplicity {
        self.multiplicity
    }

    pub fn valence(&self) -> u8 {
        self.valence
    }

    pub fn donated_pairs(&self) -> u8 {
        self.donated_pairs
    }

    pub fn accepted_pairs(&self) -> u8 {
        self.accepted_pairs
    }

    pub fn aromatic_valence(&self) -> u8 {
        self.aromatic_valence
    }

    pub fn multicenter_valence(&self) -> u8 {
        self.multicenter_valence
    }
}

impl Display for AtomTypeSpec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}", self.element)?;
        match self.charge {
            0 => {}
            1 => write!(f, "+")?,
            -1 => write!(f, "-")?,
            c if c < 0 => write!(f, "{}", c)?,
            c => write!(f, "+{}", c)?,
        }
        if self.lone_pairs > 0 {
            write!(f, "/{}", self.lone_pairs)?;
        }
        if self.unpaired_electrons > 0 {
            write!(f, "^{}", self.unpaired_electrons)?;
        }
        if self.multiplicity.multiplicity() != self.unpaired_electrons.saturating_add(1) {
            write!(f, "*{}", self.multiplicity.multiplicity())?;
        }
        if self.hydrogens > 0 {
            write!(f, "H{}", self.hydrogens)?;
        }
        if self.valence > 0 {
            write!(f, "v{}", self.valence)?;
        }
        if self.donated_pairs > 0 {
            write!(f, ">{}", self.donated_pairs)?;
        }
        if self.accepted_pairs > 0 {
            write!(f, "<{}", self.accepted_pairs)?;
        }
        if self.aromatic_valence > 0 {
            write!(f, "a{}", self.aromatic_valence)?;
        }
        if self.multicenter_valence > 0 {
            write!(f, "m{}", self.multicenter_valence)?;
        }
        write!(f, "]")
    }
}

impl FromStr for AtomTypeSpec {
    type Err = ResolutionError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if !s.starts_with('[') || !s.ends_with(']') {
            return Err(ResolutionError::InvalidAtomSpec(format!(
                "atom type spec must be bracketed: {}",
                s
            )));
        }
        let body = &s[1..s.len() - 1];
        let mut chars = body.chars().peekable();

        let first = chars
            .next()
            .ok_or_else(|| ResolutionError::InvalidAtomSpec("empty atom type spec".to_string()))?;
        if !first.is_ascii_uppercase() {
            return Err(ResolutionError::InvalidAtomSpec(format!(
                "invalid element in {}",
                s
            )));
        }
        let mut elem = String::new();
        elem.push(first);
        if chars.peek().is_some_and(|c| c.is_ascii_lowercase()) {
            elem.push(chars.next().unwrap());
        }
        let element: Element = elem
            .parse()
            .map_err(|_| ResolutionError::InvalidAtomSpec(format!("invalid element: {}", elem)))?;

        let mut charge = 0i8;
        let mut seen_charge = false;
        let mut lone_pairs = 0_u8;
        let mut multiplicity: Option<SpinMultiplicity> = None;
        let mut hydrogens = 0u8;
        let mut valence = 0u8;
        let mut donated_pairs = 0u8;
        let mut accepted_pairs = 0u8;
        let mut unpaired_electrons = 0u8;
        let mut aromatic_valence = 0u8;
        let mut multicenter_valence = 0u8;

        while let Some(token) = chars.next() {
            let mut number = String::new();
            while chars.peek().is_some_and(|c| c.is_ascii_digit()) {
                number.push(chars.next().unwrap());
            }
            let num_u8 = |default: u8| -> Result<u8, ResolutionError> {
                if number.is_empty() {
                    Ok(default)
                } else {
                    number.parse::<u8>().map_err(|_| {
                        ResolutionError::InvalidAtomSpec(format!(
                            "invalid numeric token '{}' in {}",
                            number, s
                        ))
                    })
                }
            };
            match token {
                '+' => {
                    if seen_charge {
                        return Err(ResolutionError::InvalidAtomSpec(
                            "duplicate charge token".to_string(),
                        ));
                    }
                    charge = num_u8(1)? as i8;
                    seen_charge = true;
                }
                '-' => {
                    if seen_charge {
                        return Err(ResolutionError::InvalidAtomSpec(
                            "duplicate charge token".to_string(),
                        ));
                    }
                    charge = -(num_u8(1)? as i8);
                    seen_charge = true;
                }
                '/' => lone_pairs = num_u8(1)?,
                '^' => unpaired_electrons = num_u8(1)?,
                '*' => {
                    let m = num_u8(1)?;
                    multiplicity =
                        Some(SpinMultiplicity::from_multiplicity(m).ok_or_else(|| {
                            ResolutionError::InvalidAtomSpec(format!(
                                "invalid multiplicity {} in {}",
                                m, s
                            ))
                        })?);
                }
                'H' => hydrogens = num_u8(1)?,
                'v' => valence = num_u8(1)?,
                '>' => donated_pairs = num_u8(1)?,
                '<' => accepted_pairs = num_u8(1)?,
                'a' => aromatic_valence = num_u8(1)?,
                'm' => multicenter_valence = num_u8(1)?,
                _ => {
                    return Err(ResolutionError::InvalidAtomSpec(format!(
                        "unknown token '{}' in {}",
                        token, s
                    )))
                }
            }
        }

        let multiplicity = multiplicity.unwrap_or_else(|| {
            SpinMultiplicity::from_multiplicity(unpaired_electrons.saturating_add(1))
                .unwrap_or(SpinMultiplicity::Singlet)
        });

        Self::new(
            element,
            charge,
            hydrogens,
            lone_pairs,
            unpaired_electrons,
            multiplicity,
            valence,
            donated_pairs,
            accepted_pairs,
            aromatic_valence,
            multicenter_valence,
        )
    }
}

impl Serialize for AtomTypeSpec {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for AtomTypeSpec {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let s = String::deserialize(deserializer)?;
        s.parse().map_err(de::Error::custom)
    }
}

/// Optional query constraints for matching atom type specs.
#[derive(Debug, Clone, Copy)]
pub struct AtomTypeQuery {
    pub element: Element,
    pub charge: Option<i8>,
    pub hydrogens: Option<u8>,
    pub lone_pairs: Option<u8>,
    pub unpaired_electrons: Option<u8>,
    pub multiplicity: Option<SpinMultiplicity>,
    pub valence: Option<u8>,
    pub donated_pairs: Option<u8>,
    pub accepted_pairs: Option<u8>,
    pub aromatic_valence: Option<u8>,
    pub multicenter_valence: Option<u8>,
}

impl AtomTypeQuery {
    pub fn unconstrained(element: Element) -> Self {
        Self {
            element,
            charge: None,
            hydrogens: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            valence: None,
            donated_pairs: None,
            accepted_pairs: None,
            aromatic_valence: None,
            multicenter_valence: None,
        }
    }

    pub fn from_builder_atom(builder: &MoleculeBuilder, atom_index: AtomIndex) -> Self {
        let atom = builder.atom(atom_index).expect("atom_index must be valid");
        let valence = builder.atom_bond_order_sum(atom_index);
        let (donated_pairs, accepted_pairs) = builder.atom_dative_bond_order_sums(atom_index);
        let aromatic_valence = match atom.aromatic_hint() {
            Some(false) => Some(0),
            _ => None,
        };
        let multicenter_valence = if builder.atom_has_multicenter_bonds(atom_index) {
            None
        } else {
            Some(0)
        };
        Self {
            element: atom.element(),
            charge: atom.charge(),
            hydrogens: atom.hydrogens(),
            lone_pairs: atom.lone_pairs(),
            unpaired_electrons: atom.unpaired_electrons(),
            multiplicity: atom.multiplicity(),
            valence: Some(valence),
            donated_pairs: Some(donated_pairs),
            accepted_pairs: Some(accepted_pairs),
            aromatic_valence,
            multicenter_valence,
        }
    }

    pub fn matches_spec(&self, spec: &AtomTypeSpec) -> bool {
        self.charge.is_none_or(|v| v == spec.charge())
            && self.hydrogens.is_none_or(|v| v == spec.hydrogens())
            && self.lone_pairs.is_none_or(|v| v == spec.lone_pairs())
            && self
                .unpaired_electrons
                .is_none_or(|v| v == spec.unpaired_electrons())
            && self.multiplicity.is_none_or(|v| v == spec.multiplicity())
            && self.valence.is_none_or(|v| v == spec.valence())
            && self.donated_pairs.is_none_or(|v| v == spec.donated_pairs())
            && self
                .accepted_pairs
                .is_none_or(|v| v == spec.accepted_pairs())
            && self
                .aromatic_valence
                .is_none_or(|v| v == spec.aromatic_valence())
            && self
                .multicenter_valence
                .is_none_or(|v| v == spec.multicenter_valence())
    }
}

/// Built-in atom type registry.
///
/// Each spec is stored under both `(element, Some(charge))` and `(element, None)`,
/// enabling O(1) lookup for both charge-specific and element-only queries.
#[derive(Debug, Clone, Default)]
pub struct AtomTypeRegistry {
    atom_types: HashMap<(Element, Option<i8>), Vec<AtomTypeSpec>>,
}

#[derive(Debug, Deserialize)]
struct AtomTypeRegistryToml {
    atom_types: Vec<AtomTypeSpec>,
}

impl AtomTypeRegistry {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn default_registry() -> &'static Self {
        &DEFAULT_ATOM_TYPE_REGISTRY
    }

    pub fn from_specs(specs: impl IntoIterator<Item = AtomTypeSpec>) -> Self {
        let mut reg = Self::new();
        for spec in specs {
            reg.push(spec);
        }
        reg
    }

    pub fn from_toml_str(input: &str) -> Result<Self, ResolutionError> {
        let parsed: AtomTypeRegistryToml = toml::from_str(input)
            .map_err(|e| ResolutionError::InvalidAtomTypeRegistry(e.to_string()))?;
        Ok(Self::from_specs(parsed.atom_types))
    }

    pub fn from_toml_file(path: &Path) -> Result<Self, ResolutionError> {
        let input = fs::read_to_string(path)
            .map_err(|e| ResolutionError::InvalidAtomTypeRegistry(e.to_string()))?;
        Self::from_toml_str(&input)
    }

    pub fn push(&mut self, spec: AtomTypeSpec) {
        self.atom_types
            .entry((spec.element(), Some(spec.charge())))
            .or_default()
            .push(spec);
        self.atom_types
            .entry((spec.element(), None))
            .or_default()
            .push(spec);
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
            .filter(|spec| query.matches_spec(spec))
            .copied()
            .collect()
    }
}

/// Public shorthand for parsing a single atom type specification.
#[macro_export]
macro_rules! spec {
    ($s:expr) => {{
        use std::str::FromStr;
        $crate::graph_ir::valence::AtomTypeSpec::from_str($s).unwrap()
    }};
}

/// Public shorthand for defining atom type registries.
#[macro_export]
macro_rules! registry {
    ($($spec:expr),* $(,)?) => {{
        let mut registry = $crate::graph_ir::valence::AtomTypeRegistry::new();
        $(
            registry.push($crate::spec!($spec));
        )*
        registry
    }};
}

static DEFAULT_ATOM_TYPE_REGISTRY: LazyLock<AtomTypeRegistry> = LazyLock::new(|| {
    AtomTypeRegistry::from_toml_str(include_str!("../../spec/default-registry.toml"))
        .expect("built-in default registry must be valid")
});

/// Atom-level valence validator.
#[derive(Debug, Clone)]
pub enum ValenceValidator {
    AtomTyping(AtomTypeRegistry),
    Counts,
}

impl ValenceValidator {
    pub fn candidates_for(
        &self,
        builder: &MoleculeBuilder,
        atom_index: AtomIndex,
    ) -> SmallVec<[AtomTypeSpec; 4]> {
        match self {
            ValenceValidator::AtomTyping(registry) => {
                registry.candidates_for(&AtomTypeQuery::from_builder_atom(builder, atom_index))
            }
            ValenceValidator::Counts => todo!("counts-based valence is not implemented yet"),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use umol_data::Element;

    use super::*;
    use crate::graph_ir::AtomBuilder;

    #[test]
    fn atom_type_spec_parse_extended_fields() {
        let spec = AtomTypeSpec::from_str("[N+1/1>1<0^0*1H0v3a1m0]").unwrap();
        assert_eq!(spec.element(), Element::N);
        assert_eq!(spec.charge(), 1);
        assert_eq!(spec.lone_pairs(), 1);
        assert_eq!(spec.valence(), 3);
        assert_eq!(spec.donated_pairs(), 1);
        assert_eq!(spec.accepted_pairs(), 0);
        assert_eq!(spec.aromatic_valence(), 1);
        assert_eq!(spec.multicenter_valence(), 0);
    }

    #[test]
    fn atom_type_spec_display_roundtrip() {
        let input = "[C-1/1^2*1H1v2a1m2]";
        let parsed = AtomTypeSpec::from_str(input).unwrap();
        let reparsed = AtomTypeSpec::from_str(&parsed.to_string()).unwrap();
        assert_eq!(parsed, reparsed);
    }

    #[test]
    fn atom_type_registry_from_toml() {
        let input = r#"
atom_types = ["[C+0v4a0m0]", "[O-1/3v1a0m0]"]
"#;
        let reg = AtomTypeRegistry::from_toml_str(input).unwrap();
        assert_eq!(reg.specs_for_element_and_charge(Element::C, 0).len(), 1);
        assert_eq!(reg.specs_for_element_and_charge(Element::O, -1).len(), 1);
    }

    #[test]
    fn default_registry_is_populated() {
        let reg = AtomTypeRegistry::default_registry();
        assert!(!reg.specs_for_element(Element::C).is_empty());
    }

    #[test]
    fn registry_macro_builds() {
        let reg = registry!["[C+0v4]", "[C+1^3v3]"];
        assert_eq!(reg.specs_for_element(Element::C).len(), 2);
    }

    #[test]
    fn dative_bond_order_sums_donor_and_acceptor() {
        use crate::graph_ir::DativeBond;

        let mut mb = MoleculeBuilder::new();
        let n = mb.add_atom(AtomBuilder::new(Element::N));
        let b = mb.add_atom(AtomBuilder::new(Element::B));
        mb.add_dative_bond(DativeBond::new(n, b, 2));

        assert_eq!(mb.atom_dative_bond_order_sums(n), (2, 0));
        assert_eq!(mb.atom_dative_bond_order_sums(b), (0, 2));
    }
}
