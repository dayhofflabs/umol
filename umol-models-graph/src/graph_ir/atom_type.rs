//! Atom typing specifications and queries for valence resolution.

use std::fmt::{self, Display};
use std::str::FromStr;

use serde::de::{self, Deserializer};
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};
use thiserror::Error;
use umol_data::{Element, SpinMultiplicity, SpinStateError};

use crate::atom::AromaticValence;
use crate::graph_ir::atom_pattern::HydrogenPattern;
use crate::graph_ir::error::ResolutionError;
use crate::graph_ir::molecule::{AtomIndex, MoleculeBuilder};
use crate::graph_ir::Atom;
use crate::table_ir::atom::ImplicitHydrogens;

// TODO: Incorporate relevant variants into ValidationError, remove this enum.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum AtomError {
    // Deprecated: legacy `?{...}` atom type query parser only.
    #[error("atom type query must use ?{{...}} notation")]
    InvalidQueryFormat,
    // Deprecated: legacy `?{...}` atom type query parser only.
    #[error("empty atom type query")]
    EmptyQuery,
    #[error("invalid element value: {0}")]
    InvalidElement(String),
    #[error("invalid atom tag: {0}")]
    InvalidTag(String),
    #[error("duplicate atom tag: {0}")]
    DuplicateTag(String),
    #[error("invalid charge value: {0}")]
    InvalidCharge(String),
    #[error("invalid implicit hydrogens value: {0}")]
    InvalidImplicitHydrogens(String),
    #[error("invalid lone-pairs value: {0}")]
    InvalidLonePairs(String),
    #[error("invalid unpaired-electrons value: {0}")]
    InvalidUnpairedElectrons(String),
    #[error("invalid multiplicity value: {0}")]
    InvalidMultiplicity(String),
    #[error("invalid valence value: {0}")]
    InvalidValence(String),
    #[error("invalid donated-pairs value: {0}")]
    InvalidDonatedPairs(String),
    #[error("invalid accepted-pairs value: {0}")]
    InvalidAcceptedPairs(String),
    #[error("invalid aromatic-valence value: {0}")]
    InvalidAromaticValence(String),
    #[error("invalid multicenter-valence value: {0}")]
    InvalidMulticenterValence(String),
    #[error(transparent)]
    SpinState(#[from] SpinStateError),
    #[error("unexpected atom tag: {0}")]
    UnexpectedTag(String),
    #[error("field '{field}' out of range: {value} not in [{min}, {max}]")]
    OutOfRange {
        field: &'static str,
        value: i64,
        min: i64,
        max: i64,
    },
    #[error("charge {charge} out of bounds for {element}: expected [{min_charge}, {max_charge}]")]
    ChargeOutOfBounds {
        element: Element,
        charge: i8,
        min_charge: i8,
        max_charge: i8,
    },
    #[error(
        "electron invariant mismatch for {element}: inv_o={orbital_invariant}, inv_e={electron_invariant}"
    )]
    ElectronInvariantMismatch {
        element: Element,
        orbital_invariant: i16,
        electron_invariant: i16,
    },
}

impl From<AtomError> for ResolutionError {
    fn from(value: AtomError) -> Self {
        ResolutionError::InvalidAtom(value.to_string())
    }
}

/// Constraint for matching implicit hydrogen information in atom type queries.
///
/// Query notation:
/// - `H` / `H1` / `Hn` => `Hydrogens(n)`
/// - `H*` => `Any`
/// - `H=` => `Normal`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum HydrogenConstraint {
    Hydrogens(u8),
    Normal,
    Any,
}

impl HydrogenConstraint {
    pub fn matches(&self, hydrogens: u8) -> bool {
        match self {
            HydrogenConstraint::Hydrogens(n) => *n == hydrogens,
            HydrogenConstraint::Normal => false,
            HydrogenConstraint::Any => true,
        }
    }
}

impl HydrogenConstraint {
    pub fn from_implicit_hydrogens(implicit_hydrogens: ImplicitHydrogens) -> Self {
        match implicit_hydrogens {
            ImplicitHydrogens::Hydrogens(h) => HydrogenConstraint::Hydrogens(h),
            ImplicitHydrogens::Normal => HydrogenConstraint::Normal,
        }
    }
}

/// Constraint for matching aromatic valence in atom type queries.
///
/// Variants:
/// - None: Non-aromatic
/// - Any: Aromatic (unknown valence)
/// - Valence(n): Aromatic, n >= 0 valence electrons
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AromaticConstraint {
    None,
    Any,
    Valence(u8),
}

impl AromaticConstraint {
    pub fn matches(&self, av: AromaticValence) -> bool {
        match self {
            AromaticConstraint::None => av == AromaticValence::NotAromatic,
            AromaticConstraint::Any => av.is_aromatic(),
            AromaticConstraint::Valence(n) => av == AromaticValence::Valence(*n),
        }
    }
}

/// Optional query constraints for matching atom type specs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct AtomTypeQuery {
    pub element: Element,
    pub charge: Option<i8>,
    pub implicit_hydrogens: Option<HydrogenConstraint>,
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
            implicit_hydrogens: None,
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
        let hydrogen_constraint = if builder.atom_has_normal_implicit_hydrogens(atom_index) {
            Some(HydrogenConstraint::Normal)
        } else {
            match &atom.implicit_hydrogens {
                HydrogenPattern::Is(h) => Some(HydrogenConstraint::Hydrogens(*h)),
                HydrogenPattern::Normal => Some(HydrogenConstraint::Normal),
                HydrogenPattern::Any => None,
            }
        };
        let aromatic_constraint = if builder.atom_aromatic_hint(atom_index) {
            Some(AromaticConstraint::Any)
        } else if builder.atom_explicit_aromatic_hint(atom_index) == Some(false) {
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
            charge: atom.charge.into_option(),
            implicit_hydrogens: hydrogen_constraint,
            lone_pairs: atom.lone_pairs.into_option(),
            unpaired_electrons: atom.unpaired_electrons.into_option(),
            multiplicity: atom.multiplicity.into_option(),
            valence: Some(valence),
            donated_pairs: Some(donated_pairs),
            accepted_pairs: Some(accepted_pairs),
            aromatic_valence: aromatic_constraint,
            multicenter_valence,
        }
    }

    pub fn matches_atom(&self, atom: &Atom) -> bool {
        self.charge.is_none_or(|v| v == atom.charge())
            && self
                .implicit_hydrogens
                .is_none_or(|v| v.matches(atom.implicit_hydrogens()))
            && self.lone_pairs.is_none_or(|v| v == atom.lone_pairs())
            && self
                .unpaired_electrons
                .is_none_or(|v| v == atom.unpaired_electrons())
            && self.multiplicity.is_none_or(|v| v == atom.multiplicity())
            && self.valence.is_none_or(|v| v == atom.valence())
            && self.donated_pairs.is_none_or(|v| v == atom.donated_pairs())
            && self
                .accepted_pairs
                .is_none_or(|v| v == atom.accepted_pairs())
            && self
                .aromatic_valence
                .is_none_or(|c| c.matches(atom.aromatic_valence()))
            && self
                .multicenter_valence
                .is_none_or(|v| v == atom.multicenter_valence())
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
        if let Some(h) = self.implicit_hydrogens {
            match h {
                HydrogenConstraint::Hydrogens(1) => write!(f, "H")?,
                HydrogenConstraint::Hydrogens(n) => write!(f, "H{}", n)?,
                HydrogenConstraint::Normal => write!(f, "H=")?,
                HydrogenConstraint::Any => write!(f, "H*")?,
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
    type Err = AtomError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let trimmed = s.trim();
        if !trimmed.starts_with("?{") || !trimmed.ends_with('}') {
            return Err(AtomError::InvalidQueryFormat);
        }
        let body = &trimmed[2..trimmed.len() - 1];
        let mut chars = body.chars().peekable();

        let first = chars.next().ok_or(AtomError::EmptyQuery)?;
        if !first.is_ascii_uppercase() {
            return Err(AtomError::InvalidElement(body.to_string()));
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
            .map_err(|_| AtomError::InvalidElement(elem.clone()))?;

        let mut query = AtomTypeQuery::unconstrained(element);
        let mut seen_charge = false;

        while let Some(token) = chars.next() {
            if token.is_ascii_whitespace() {
                continue;
            }
            let mut number = String::new();
            while chars.peek().is_some_and(|c| c.is_ascii_digit()) {
                number.push(chars.next().unwrap());
            }
            let num_u8 = |default: u8, tag: &str| -> Result<u8, AtomError> {
                if number.is_empty() {
                    Ok(default)
                } else {
                    number
                        .parse::<u8>()
                        .map_err(|_| AtomError::UnexpectedTag(tag.to_string()))
                }
            };
            match token {
                '+' => {
                    if seen_charge {
                        return Err(AtomError::DuplicateTag("+".to_string()));
                    }
                    query.charge = Some(num_u8(1, "+")? as i8);
                    seen_charge = true;
                }
                '-' => {
                    if seen_charge {
                        return Err(AtomError::DuplicateTag("-".to_string()));
                    }
                    query.charge = Some(-(num_u8(1, "-")? as i8));
                    seen_charge = true;
                }
                'H' => {
                    if chars.peek() == Some(&'*') {
                        chars.next();
                        query.implicit_hydrogens = Some(HydrogenConstraint::Any);
                    } else if chars.peek() == Some(&'=') {
                        chars.next();
                        query.implicit_hydrogens = Some(HydrogenConstraint::Normal);
                    } else {
                        query.implicit_hydrogens =
                            Some(HydrogenConstraint::Hydrogens(num_u8(1, "H")?));
                    }
                }
                '/' => query.lone_pairs = Some(num_u8(1, "/")?),
                '^' => query.unpaired_electrons = Some(num_u8(1, "^")?),
                'x' => {
                    let m = num_u8(1, "x")?;
                    query.multiplicity = Some(
                        SpinMultiplicity::from_multiplicity(m)
                            .ok_or_else(|| AtomError::InvalidMultiplicity(m.to_string()))?,
                    );
                }
                'v' => query.valence = Some(num_u8(1, "v")?),
                '>' => query.donated_pairs = Some(num_u8(1, ">")?),
                '<' => query.accepted_pairs = Some(num_u8(1, "<")?),
                'a' => {
                    if chars.peek() == Some(&'*') {
                        chars.next();
                        query.aromatic_valence = Some(AromaticConstraint::Any);
                    } else if chars.peek() == Some(&'!') {
                        chars.next();
                        query.aromatic_valence = Some(AromaticConstraint::None);
                    } else {
                        query.aromatic_valence = Some(AromaticConstraint::Valence(num_u8(1, "a")?));
                    }
                }
                'm' => query.multicenter_valence = Some(num_u8(1, "m")?),
                _ => {
                    return Err(AtomError::UnexpectedTag(token.to_string()));
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
    use umol_data::{Element, SpinMultiplicity};

    use super::*;

    #[rustfmt::skip]
    #[rstest]
    #[case::unconstrained("?{C}", Element::C, None, None, None, None, None, None, None, None)]
    #[case::outer_spaced("  ?{C}  ", Element::C, None, None, None, None, None, None, None, None)]
    #[case::hydrogen_any("?{CH*}", Element::C, None, Some(HydrogenConstraint::Any), None, None, None, None, None, None)]
    #[case::hydrogen_normal("?{CH=}", Element::C, None, Some(HydrogenConstraint::Normal), None, None, None, None, None, None)]
    #[case::spaced_aromatic_none("?{B a!}", Element::B, None, None, None, None, None, None, Some(AromaticConstraint::None), None)]
    #[case::constrained("?{C-H/1^2x1v2a1m2}", Element::C, Some(-1), Some(HydrogenConstraint::Hydrogens(1)), Some(1), Some(2), Some(SpinMultiplicity::Singlet), Some(2), Some(AromaticConstraint::Valence(1)), Some(2))]
    #[case::aromatic_any("?{Cv2a*}", Element::C, None, None, None, None, None, Some(2), Some(AromaticConstraint::Any), None)]
    #[case::aromatic_none("?{Cv2a!}", Element::C, None, None, None, None, None, Some(2), Some(AromaticConstraint::None), None)]
    fn test_atom_type_query_from_str(
        #[case] input: &str,
        #[case] element: Element,
        #[case] charge: Option<i8>,
        #[case] hydrogens: Option<HydrogenConstraint>,
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
        assert_eq!(query.implicit_hydrogens, hydrogens, "hydrogens mismatch for {}", input);
        assert_eq!(query.lone_pairs, lone_pairs, "lone pairs mismatch for {}", input);
        assert_eq!(query.unpaired_electrons, unpaired_electrons, "unpaired electrons mismatch for {}", input);
        assert_eq!(query.multiplicity, multiplicity, "multiplicity mismatch for {}", input);
        assert_eq!(query.valence, valence, "valence mismatch for {}", input);
        assert_eq!(query.aromatic_valence, aromatic_valence, "aromatic valence mismatch for {}", input);
        assert_eq!(query.multicenter_valence, multicenter_valence, "multicenter valence mismatch for {}", input);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::invalid_format("{C}", AtomError::InvalidQueryFormat)]
    #[case::empty_query("?{}", AtomError::EmptyQuery)]
    #[case::invalid_element("?{c}", AtomError::InvalidElement("c".to_string()))]
    #[case::duplicate_charge("?{C+-}", AtomError::DuplicateTag("-".to_string()))]
    #[case::invalid_multiplicity("?{Cx11}", AtomError::InvalidMultiplicity("11".to_string()))]
    #[case::unexpected_token("?{Cq1}", AtomError::UnexpectedTag("q".to_string()))]
    fn test_atom_type_query_from_str_error(#[case] input: &str, #[case] expected: AtomError) {
        assert_eq!(AtomTypeQuery::from_str(input).unwrap_err(), expected);
    }

    #[rstest]
    #[case::unconstrained("?{C}")]
    #[case::hydrogen_any("?{CH*}")]
    #[case::hydrogen_normal("?{CH=}")]
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
    #[case::any_rejects_none(AromaticConstraint::Any, AromaticValence::NotAromatic, false)]
    #[case::none_matches_none(AromaticConstraint::None, AromaticValence::NotAromatic, true)]
    #[case::none_rejects_a1(AromaticConstraint::None, AromaticValence::Valence(1), false)]
    #[case::exact_matches(AromaticConstraint::Valence(2), AromaticValence::Valence(2), true)]
    #[case::exact_rejects_wrong(AromaticConstraint::Valence(2), AromaticValence::Valence(1), false)]
    #[case::exact_rejects_none(AromaticConstraint::Valence(1), AromaticValence::NotAromatic, false)]
    fn test_aromatic_constraint_matches(
        #[case] constraint: AromaticConstraint,
        #[case] valence: AromaticValence,
        #[case] expected: bool,
    ) {
        assert_eq!(constraint.matches(valence), expected);
    }
}
