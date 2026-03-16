//! Atom typing specifications and queries for valence resolution.

use std::fmt::{self, Display};
use std::str::FromStr;

use serde::de::{self, Deserializer};
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};
use umol_data::{Element, SpinMultiplicity, SpinState, MAX_UNPAIRED_ELECTRONS};

use crate::graph_ir::error::ResolutionError;
use crate::graph_ir::molecule::{AtomIndex, MoleculeBuilder};

/// Aromatic valence of an atom: either non-aromatic or contributing n >= 0
/// valence to a delocalized pi-system. Each atom can participate in at
/// most one aromatic system.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum AromaticValence {
    /// Non-aromatic atom.
    None,
    /// Aromatic atom contributing  valence `n` (n >= 0)
    Valence(u8),
}

impl AromaticValence {
    pub fn valence(&self) -> u8 {
        match self {
            AromaticValence::None => 0,
            AromaticValence::Valence(n) => *n,
        }
    }

    /// Atom is aromatic if it contributes valence (n >= 0) to an aromatic system
    pub fn is_aromatic(&self) -> bool {
        matches!(self, AromaticValence::Valence(_))
    }
}

impl Display for AromaticValence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AromaticValence::None => Ok(()),
            AromaticValence::Valence(n) => write!(f, "a{}", n),
        }
    }
}

impl FromStr for AromaticValence {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some(rest) = s.strip_prefix('a') {
            let n: u8 = rest
                .parse()
                .map_err(|_| format!("invalid aromatic valence: {}", s))?;
            Ok(AromaticValence::Valence(n))
        } else {
            Err(format!("expected 'a' prefix: {}", s))
        }
    }
}

/// Constraint for matching aromatic valence in atom type queries.
///
/// Used by `AtomTypeQuery` to filter candidates during valence resolution.
/// Unlike `AromaticValence` (which is a concrete value on a spec), this
/// expresses what range of aromatic valences is acceptable.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AromaticConstraint {
    /// Must be non-aromatic (`AromaticValence::None`).
    None,
    /// Must be aromatic with any pi-electron count (`AromaticValence::Valence(_)`).
    Any,
    /// Must be aromatic with exactly `n` pi-electrons (`AromaticValence::Valence(n)`).
    Valence(u8),
}

impl AromaticConstraint {
    pub fn matches(&self, av: AromaticValence) -> bool {
        match self {
            AromaticConstraint::None => av == AromaticValence::None,
            AromaticConstraint::Any => av.is_aromatic(),
            AromaticConstraint::Valence(n) => av == AromaticValence::Valence(*n),
        }
    }
}

/// Atom typing specification for valence resolution.
///
/// String notation:
/// - `{El...}` where `El` is an element symbol.
/// - tokens are optional and can appear in any order:
///   - `+n` / `-n` charge (default 0, bare `+`/`-` means 1)
///   - `Hn` hydrogens
///   - `/n` lone pairs
///   - `^n` unpaired electrons
///   - `xn` multiplicity (default `unpaired + 1`)
///   - `vn` valence
///   - `>n` donated pairs
///   - `<n` accepted pairs
///   - `an` aromatic valence (n >= 0) or none (non-aromatic)
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
    aromatic_valence: AromaticValence,
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
        aromatic_valence: AromaticValence,
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

    pub fn aromatic_valence(&self) -> AromaticValence {
        self.aromatic_valence
    }

    pub fn multicenter_valence(&self) -> u8 {
        self.multicenter_valence
    }

    pub fn is_aromatic(&self) -> bool {
        self.aromatic_valence.is_aromatic()
    }
}

impl Display for AtomTypeSpec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{{{}", self.element)?;
        match self.charge {
            0 => {}
            1 => write!(f, "+")?,
            -1 => write!(f, "-")?,
            c if c < 0 => write!(f, "{}", c)?,
            c => write!(f, "+{}", c)?,
        }
        if self.hydrogens > 0 {
            if self.hydrogens == 1 {
                write!(f, "H")?;
            } else {
                write!(f, "H{}", self.hydrogens)?;
            }
        }
        if self.lone_pairs > 0 {
            write!(f, "/{}", self.lone_pairs)?;
        }
        let n = self.spin.unpaired_electrons();
        let m = self.spin.multiplicity();
        if n > 0 {
            write!(f, "^{}", n)?;
        }
        if m.multiplicity() != n + 1 {
            write!(f, "x{}", m.multiplicity())?;
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
        write!(f, "{}", self.aromatic_valence)?;
        if self.multicenter_valence > 0 {
            write!(f, "m{}", self.multicenter_valence)?;
        }
        write!(f, "}}")
    }
}

impl FromStr for AtomTypeSpec {
    type Err = ResolutionError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if !s.starts_with('{') || !s.ends_with('}') {
            return Err(ResolutionError::InvalidAtomSpec(format!(
                "atom type spec must be braced: {}",
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
        if let Some(&c) = chars.peek() {
            if c.is_ascii_lowercase() {
                let mut two = String::new();
                two.push(first);
                two.push(c);
                if two.parse::<Element>().is_ok() {
                    elem.push(chars.next().unwrap());
                }
            }
        }
        let element: Element = elem
            .parse()
            .map_err(|_| ResolutionError::InvalidAtomSpec(format!("invalid element: {}", elem)))?;

        let mut charge = None;
        let mut hydrogens = 0u8;
        let mut lone_pairs = 0_u8;
        let mut multiplicity: Option<SpinMultiplicity> = None;
        let mut valence = 0u8;
        let mut donated_pairs = 0u8;
        let mut accepted_pairs = 0u8;
        let mut unpaired_electrons = 0u8;
        let mut aromatic_valence = AromaticValence::None;
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
                    if charge.is_some() {
                        return Err(ResolutionError::InvalidAtomSpec(
                            "duplicate charge token".to_string(),
                        ));
                    }
                    charge = Some(num_u8(1)? as i8);
                }
                '-' => {
                    if charge.is_some() {
                        return Err(ResolutionError::InvalidAtomSpec(
                            "duplicate charge token".to_string(),
                        ));
                    }
                    charge = Some(-(num_u8(1)? as i8));
                }
                'H' => hydrogens = num_u8(1)?,
                '/' => lone_pairs = num_u8(1)?,
                '^' => unpaired_electrons = num_u8(1)?,
                'x' => {
                    let m = num_u8(1)?;
                    multiplicity =
                        Some(SpinMultiplicity::from_multiplicity(m).ok_or_else(|| {
                            ResolutionError::InvalidAtomSpec(format!(
                                "invalid multiplicity {} in {}",
                                m, s
                            ))
                        })?);
                }
                'v' => valence = num_u8(1)?,
                '>' => donated_pairs = num_u8(1)?,
                '<' => accepted_pairs = num_u8(1)?,
                'a' => aromatic_valence = AromaticValence::Valence(num_u8(1)?),
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
            charge.unwrap_or(0),
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
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
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
    pub aromatic_valence: Option<AromaticConstraint>,
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
        let aromatic_valence = if builder.atom_aromatic_hint(atom_index) {
            Some(AromaticConstraint::Any)
        } else if atom.aromatic_hint() == Some(false) {
            Some(AromaticConstraint::None)
        } else {
            None
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

    pub fn matches(&self, spec: &AtomTypeSpec) -> bool {
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
                .is_none_or(|c| c.matches(spec.aromatic_valence()))
            && self
                .multicenter_valence
                .is_none_or(|v| v == spec.multicenter_valence())
    }

    pub fn is_aromatic(&self) -> bool {
        self.aromatic_valence
            .is_some_and(|c| matches!(c, AromaticConstraint::Any | AromaticConstraint::Valence(_)))
    }
}

impl Display for AtomTypeQuery {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "?{{{}", self.element)?;
        match self.charge {
            None => {}
            Some(0) => write!(f, "+0")?,
            Some(1) => write!(f, "+")?,
            Some(-1) => write!(f, "-")?,
            Some(c) if c < 0 => write!(f, "{}", c)?,
            Some(c) => write!(f, "+{}", c)?,
        }
        if let Some(h) = self.hydrogens {
            if h == 1 {
                write!(f, "H")?;
            } else {
                write!(f, "H{}", h)?;
            }
        }
        if let Some(lp) = self.lone_pairs {
            write!(f, "/{}", lp)?;
        }
        if let Some(n) = self.unpaired_electrons {
            write!(f, "^{}", n)?;
        }
        if let Some(m) = self.multiplicity {
            write!(f, "x{}", m.multiplicity())?;
        }
        if let Some(v) = self.valence {
            write!(f, "v{}", v)?;
        }
        if let Some(d) = self.donated_pairs {
            write!(f, ">{}", d)?;
        }
        if let Some(a) = self.accepted_pairs {
            write!(f, "<{}", a)?;
        }
        match self.aromatic_valence {
            Some(AromaticConstraint::None) => write!(f, "a!")?,
            Some(AromaticConstraint::Any) => write!(f, "a*")?,
            Some(AromaticConstraint::Valence(n)) => write!(f, "a{}", n)?,
            None => {}
        }
        if let Some(mv) = self.multicenter_valence {
            write!(f, "m{}", mv)?;
        }
        write!(f, "}}")
    }
}

impl FromStr for AtomTypeQuery {
    type Err = ResolutionError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if !s.starts_with("?{") || !s.ends_with('}') {
            return Err(ResolutionError::InvalidAtomSpec(format!(
                "atom type query must use ?{{...}} notation: {}",
                s
            )));
        }
        let body = &s[2..s.len() - 1];
        let mut chars = body.chars().peekable();

        let first = chars
            .next()
            .ok_or_else(|| ResolutionError::InvalidAtomSpec("empty atom type query".to_string()))?;
        if !first.is_ascii_uppercase() {
            return Err(ResolutionError::InvalidAtomSpec(format!(
                "invalid element in {}",
                s
            )));
        }
        let mut elem = String::new();
        elem.push(first);
        if let Some(&c) = chars.peek() {
            if c.is_ascii_lowercase() {
                let mut two = String::new();
                two.push(first);
                two.push(c);
                if two.parse::<Element>().is_ok() {
                    elem.push(chars.next().unwrap());
                }
            }
        }
        let element: Element = elem
            .parse()
            .map_err(|_| ResolutionError::InvalidAtomSpec(format!("invalid element: {}", elem)))?;

        let mut query = AtomTypeQuery::unconstrained(element);
        let mut seen_charge = false;

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
                    query.charge = Some(num_u8(1)? as i8);
                    seen_charge = true;
                }
                '-' => {
                    if seen_charge {
                        return Err(ResolutionError::InvalidAtomSpec(
                            "duplicate charge token".to_string(),
                        ));
                    }
                    query.charge = Some(-(num_u8(1)? as i8));
                    seen_charge = true;
                }
                'H' => query.hydrogens = Some(num_u8(1)?),
                '/' => query.lone_pairs = Some(num_u8(1)?),
                '^' => query.unpaired_electrons = Some(num_u8(1)?),
                'x' => {
                    let m = num_u8(1)?;
                    query.multiplicity =
                        Some(SpinMultiplicity::from_multiplicity(m).ok_or_else(|| {
                            ResolutionError::InvalidAtomSpec(format!(
                                "invalid multiplicity {} in {}",
                                m, s
                            ))
                        })?);
                }
                'v' => query.valence = Some(num_u8(1)?),
                '>' => query.donated_pairs = Some(num_u8(1)?),
                '<' => query.accepted_pairs = Some(num_u8(1)?),
                'a' => {
                    if chars.peek() == Some(&'*') {
                        chars.next();
                        query.aromatic_valence = Some(AromaticConstraint::Any);
                    } else if chars.peek() == Some(&'!') {
                        chars.next();
                        query.aromatic_valence = Some(AromaticConstraint::None);
                    } else {
                        query.aromatic_valence = Some(AromaticConstraint::Valence(num_u8(1)?));
                    }
                }
                'm' => query.multicenter_valence = Some(num_u8(1)?),
                _ => {
                    return Err(ResolutionError::InvalidAtomSpec(format!(
                        "unknown token '{}' in {}",
                        token, s
                    )))
                }
            }
        }

        Ok(query)
    }
}

impl Serialize for AtomTypeQuery {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for AtomTypeQuery {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let s = String::deserialize(deserializer)?;
        s.parse().map_err(de::Error::custom)
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

/// Public shorthand for parsing a single atom type query.
#[macro_export]
macro_rules! query {
    ($s:expr) => {{
        use std::str::FromStr;
        $crate::graph_ir::atom_type::AtomTypeQuery::from_str($s).unwrap()
    }};
}

#[cfg(test)]
mod tests {

    use std::str::FromStr;

    use pretty_assertions::assert_eq;
    use rstest::*;
    use umol_data::Element;

    use super::*;

    #[test]
    fn test_aromatic_valence_display() {
        assert_eq!(AromaticValence::None.to_string(), "");
        assert_eq!(AromaticValence::Valence(0).to_string(), "a0");
        assert_eq!(AromaticValence::Valence(1).to_string(), "a1");
    }

    #[test]
    fn test_aromatic_valence_from_str() {
        assert!(AromaticValence::from_str("").is_err());
        assert_eq!(
            AromaticValence::from_str("a0").unwrap(),
            AromaticValence::Valence(0)
        );
        assert_eq!(
            AromaticValence::from_str("a1").unwrap(),
            AromaticValence::Valence(1)
        );
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::atom("{N}", Element::N, 0, 0, 0, 0, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::None, 0)]
    #[case::charge_plus("{N+}", Element::N, 1, 0, 0, 0, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::None, 0)]
    #[case::charge_minus("{N-}", Element::N, -1, 0, 0, 0, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::None, 0)]
    #[case::charge_minus_1("{N-1}", Element::N, -1, 0, 0, 0, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::None, 0)]
    #[case::charge_plus_1("{N+1}", Element::N, 1, 0, 0, 0, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::None, 0)]
    #[case::hydrogen("{NH}", Element::N, 0, 1, 0, 0, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::None, 0)]
    #[case::hydrogen1("{NH1}", Element::N, 0, 1, 0, 0, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::None, 0)]
    #[case::lone_pairs("{N/1}", Element::N, 0, 0, 1, 0, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::None, 0)]
    #[case::unpaired_electrons("{N^1}", Element::N, 0, 0, 0, 1, SpinMultiplicity::Doublet, 0, 0, 0, AromaticValence::None, 0)]
    #[case::multiplicity("{Nx1}", Element::N, 0, 0, 0, 0, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::None, 0)]
    #[case::valence("{Nv1}", Element::N, 0, 0, 0, 0, SpinMultiplicity::Singlet, 1, 0, 0, AromaticValence::None, 0)]
    #[case::donated_pairs("{N>1}", Element::N, 0, 0, 0, 0, SpinMultiplicity::Singlet, 0, 1, 0, AromaticValence::None, 0)]
    #[case::accepted_pairs("{N<1}", Element::N, 0, 0, 0, 0, SpinMultiplicity::Singlet, 0, 0, 1, AromaticValence::None, 0)]
    #[case::aromatic_valence_0("{N+0a0}", Element::N, 0, 0, 0, 0, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::Valence(0), 0)]
    #[case::aromatic_valence_1("{N+0a1}", Element::N, 0, 0, 0, 0, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::Valence(1), 0)]
    #[case::multicenter_valence_0("{Nm0}", Element::N, 0, 0, 0, 0, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::None, 0)]
    #[case::multicenter_valence_1("{Nm1}", Element::N, 0, 0, 0, 0, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::None, 1)]
    #[case::complete("{N-H/1^2x1v2a1m2}", Element::N, -1, 1, 1, 2, SpinMultiplicity::Singlet, 2, 0, 0, AromaticValence::Valence(1), 2)]
    #[case::permuted("{N^2v2a1m2-H/1^2x1}", Element::N, -1, 1, 1, 2, SpinMultiplicity::Singlet, 2, 0, 0, AromaticValence::Valence(1), 2)]
    fn test_atom_type_spec_from_str(
        #[case] input: &str,
        #[case] element: Element,
        #[case] charge: i8,
        #[case] hydrogens: u8,
        #[case] lone_pairs: u8,
        #[case] unpaired_electrons: u8,
        #[case] multiplicity: SpinMultiplicity,
        #[case] valence: u8,
        #[case] donated_pairs: u8,
        #[case] accepted_pairs: u8,
        #[case] aromatic_valence: AromaticValence,
        #[case] multicenter_valence: u8,
    ) {
        let spec = AtomTypeSpec::from_str(input).unwrap();
        assert_eq!(spec.element(), element, "element mismatch for {}", input);
        assert_eq!(spec.charge(), charge, "charge mismatch for {}", input);
        assert_eq!(spec.hydrogens(), hydrogens, "hydrogens mismatch for {}", input);
        assert_eq!(spec.lone_pairs(), lone_pairs, "lone pairs mismatch for {}", input);
        assert_eq!(spec.unpaired_electrons(), unpaired_electrons, "unpaired electrons mismatch for {}", input);
        assert_eq!(spec.multiplicity(), multiplicity, "multiplicity mismatch for {}", input);
        assert_eq!(spec.valence(), valence, "valence mismatch for {}", input);
        assert_eq!(spec.donated_pairs(), donated_pairs, "donated pairs mismatch for {}", input);
        assert_eq!(spec.accepted_pairs(), accepted_pairs, "accepted pairs mismatch for {}", input);
        assert_eq!(spec.aromatic_valence(), aromatic_valence, "aromatic valence mismatch for {}", input);
        assert_eq!(spec.multicenter_valence(), multicenter_valence, "multicenter valence mismatch for {}", input);
    }

    #[rstest]
    #[case::aromatic_a2("{C-Hv2a2}")]
    #[case::aromatic_a0("{C+Hv2a0}")]
    #[case::non_aromatic("{CH3v1}")]
    #[case::multicenter_m2("{C-H/1^2x1v2m2}")]
    // TODO: Fix multicenter valence
    // #[case::multicenter_m0("{C-H/1^2x1v2m0}")]
    fn test_atom_type_spec_display_roundtrip(#[case] input: &str) {
        let parsed = AtomTypeSpec::from_str(input).unwrap();
        let formatted = parsed.to_string();
        assert_eq!(input, formatted);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::unconstrained("?{C}", Element::C, None, None, None, None, None, None, None, None)]
    #[case::constrained("?{C-H/1^2x1v2a1m2}", Element::C, Some(-1), Some(1), Some(1), Some(2), Some(SpinMultiplicity::Singlet), Some(2), Some(AromaticConstraint::Valence(1)), Some(2))]
    #[case::aromatic_any("?{Cv2a*}", Element::C, None, None, None, None, None, Some(2), Some(AromaticConstraint::Any), None)]
    #[case::aromatic_none("?{Cv2a!}", Element::C, None, None, None, None, None, Some(2), Some(AromaticConstraint::None), None)]
    fn test_atom_type_query_from_str(
        #[case] input: &str,
        #[case] element: Element,
        #[case] charge: Option<i8>,
        #[case] hydrogens: Option<u8>,
        #[case] lone_pairs: Option<u8>,
        #[case] unpaired_electrons: Option<u8>,
        #[case] multiplicity: Option<SpinMultiplicity>,
        #[case] valence: Option<u8>,
        #[case] aromatic_valence: Option<AromaticConstraint>,
        #[case] multicenter_valence: Option<u8>,
    ) {
        let query = AtomTypeQuery::from_str(input).unwrap();
        assert_eq!(query.element, element, "element mismatch for {}", input);
        assert_eq!(query.charge, charge, "charge mismatch for {}", input);
        assert_eq!(query.hydrogens, hydrogens, "hydrogens mismatch for {}", input);
        assert_eq!(query.lone_pairs, lone_pairs, "lone pairs mismatch for {}", input);
        assert_eq!(query.unpaired_electrons, unpaired_electrons, "unpaired electrons mismatch for {}", input);
        assert_eq!(query.multiplicity, multiplicity, "multiplicity mismatch for {}", input);
        assert_eq!(query.valence, valence, "valence mismatch for {}", input);
        assert_eq!(query.aromatic_valence, aromatic_valence, "aromatic valence mismatch for {}", input);
        assert_eq!(query.multicenter_valence, multicenter_valence, "multicenter valence mismatch for {}", input);
    }

    #[rstest]
    #[case::unconstrained("?{C}")]
    #[case::constrained("?{C-H/1^2x1v2a1m2}")]
    #[case::aromatic_any("?{Cv2a*}")]
    #[case::aromatic_none("?{Cv2a!}")]
    fn test_atom_type_query_display_roundtrip(#[case] input: &str) {
        let parsed = AtomTypeQuery::from_str(input).unwrap();
        let formatted = parsed.to_string();
        assert_eq!(input, formatted);
    }

    #[rstest]
    #[case::any_matches_a1(AromaticConstraint::Any, AromaticValence::Valence(1), true)]
    #[case::any_matches_a0(AromaticConstraint::Any, AromaticValence::Valence(0), true)]
    #[case::any_rejects_none(AromaticConstraint::Any, AromaticValence::None, false)]
    #[case::none_matches_none(AromaticConstraint::None, AromaticValence::None, true)]
    #[case::none_rejects_a1(AromaticConstraint::None, AromaticValence::Valence(1), false)]
    #[case::exact_matches(AromaticConstraint::Valence(2), AromaticValence::Valence(2), true)]
    #[case::exact_rejects_wrong(AromaticConstraint::Valence(2), AromaticValence::Valence(1), false)]
    #[case::exact_rejects_none(AromaticConstraint::Valence(1), AromaticValence::None, false)]
    fn test_aromatic_constraint_matches(
        #[case] constraint: AromaticConstraint,
        #[case] valence: AromaticValence,
        #[case] expected: bool,
    ) {
        assert_eq!(constraint.matches(valence), expected);
    }
}
