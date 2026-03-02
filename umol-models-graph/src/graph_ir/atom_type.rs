//! Atom typing specifications and queries for valence resolution.

use std::fmt::{self, Display};
use std::str::FromStr;

use serde::de::{self, Deserializer};
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};
use umol_data::{Element, SpinMultiplicity, SpinState, MAX_UNPAIRED_ELECTRONS};

use super::error::ResolutionError;
use super::molecule::{AtomIndex, MoleculeBuilder};

/// Atom typing specification for valence resolution.
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
    spin: SpinState,
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
        let spin = SpinState::try_new(unpaired_electrons, multiplicity).ok_or_else(|| {
            ResolutionError::InvalidAtomSpec(format!(
                "invalid spin state: {} unpaired electrons, {} multiplicity",
                unpaired_electrons,
                multiplicity.multiplicity()
            ))
        })?;
        Ok(Self {
            element,
            charge,
            hydrogens,
            lone_pairs,
            spin,
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

    pub fn spin(&self) -> SpinState {
        self.spin
    }

    pub fn unpaired_electrons(&self) -> u8 {
        self.spin.unpaired_electrons()
    }

    pub fn multiplicity(&self) -> SpinMultiplicity {
        self.spin.multiplicity()
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
        let n = self.spin.unpaired_electrons();
        let m = self.spin.multiplicity();
        if n > 0 {
            write!(f, "^{}", n)?;
        }
        if m.multiplicity() != n.saturating_add(1) {
            write!(f, "*{}", m.multiplicity())?;
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

        if unpaired_electrons > MAX_UNPAIRED_ELECTRONS {
            return Err(ResolutionError::InvalidAtomSpec(format!(
                "unpaired electrons {} exceeds max ({})",
                unpaired_electrons, MAX_UNPAIRED_ELECTRONS
            )));
        }

        let multiplicity = match multiplicity {
            Some(m) => m,
            None => SpinState::max_multiplicity(unpaired_electrons)
                .ok_or_else(|| {
                    ResolutionError::InvalidAtomSpec(format!(
                        "cannot derive multiplicity for {} unpaired electrons in {}",
                        unpaired_electrons, s
                    ))
                })?
                .multiplicity(),
        };

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

/// Public shorthand for parsing a single atom type specification.
#[macro_export]
macro_rules! spec {
    ($s:expr) => {{
        use std::str::FromStr;
        $crate::graph_ir::atom_type::AtomTypeSpec::from_str($s).unwrap()
    }};
}

#[cfg(test)]
mod tests {

    use std::str::FromStr;

    use umol_data::Element;

    use super::*;

    #[test]
    fn atom_type_spec_parse_extended_fields() {
        let spec = AtomTypeSpec::from_str("[N+1/1^0*1H0v3>1<0a1m0]").unwrap();
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
}
