//! Atom typing specifications and queries for valence resolution.

use std::fmt::{self, Display};
use std::str::FromStr;

use serde::de::{self, Deserializer};
use serde::ser::Serializer;
use serde::{Deserialize, Serialize};
use thiserror::Error;
use umol_data::{Element, SpinMultiplicity, SpinState, SpinStateError, MAX_UNPAIRED_ELECTRONS};

use crate::atom::{AromaticValence, ImplicitHydrogens};
use crate::graph_ir::error::ResolutionError;
use crate::graph_ir::molecule::{AtomIndex, MoleculeBuilder};

#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum AtomError {
    #[error("atom type spec must use {{...}} notation")]
    InvalidSpecFormat,
    #[error("atom type query must use ?{{...}} notation")]
    InvalidQueryFormat,
    #[error("empty atom type spec")]
    EmptySpec,
    #[error("empty atom type query")]
    EmptyQuery,
    #[error("invalid element token")]
    InvalidElement,
    #[error("duplicate charge token")]
    DuplicateChargeToken,
    #[error("invalid charge token")]
    InvalidCharge,
    #[error("invalid implicit hydrogens token")]
    InvalidImplicitHydrogens,
    #[error("invalid lone-pairs token")]
    InvalidLonePairs,
    #[error("invalid unpaired-electrons token")]
    InvalidUnpairedElectrons,
    #[error("invalid multiplicity token")]
    InvalidMultiplicity,
    #[error("invalid valence token")]
    InvalidValence,
    #[error("invalid donated-pairs token")]
    InvalidDonatedPairs,
    #[error("invalid accepted-pairs token")]
    InvalidAcceptedPairs,
    #[error("invalid aromatic-valence token")]
    InvalidAromaticValence,
    #[error("invalid multicenter-valence token")]
    InvalidMulticenterValence,
    #[error("unpaired electrons exceed maximum: {unpaired_electrons} > {max_unpaired_electrons}")]
    UnpairedElectronsLiteralExceedMax {
        unpaired_electrons: u8,
        max_unpaired_electrons: u8,
    },
    #[error(transparent)]
    SpinState(#[from] SpinStateError),
    #[error("unexpected token '{token}'")]
    UnexpectedToken { token: char },
    #[error("charge {charge} out of bounds for {element}: expected [{min_charge}, {max_charge}]")]
    ChargeOutOfBounds {
        element: Element,
        charge: i8,
        min_charge: i8,
        max_charge: i8,
    },
    #[error("valence {valence} exceeds max {max_valence} for {element}")]
    ValenceExceedsMax {
        element: Element,
        valence: u8,
        max_valence: u8,
    },
    #[error(
        "unpaired electrons {unpaired_electrons} exceed max {max_unpaired_electrons} for {element}"
    )]
    UnpairedElectronsExceedMax {
        element: Element,
        unpaired_electrons: u8,
        max_unpaired_electrons: u8,
    },
    #[error(
        "implicit hydrogens {implicit_hydrogens} exceed max {max_implicit_hydrogens} for {element}"
    )]
    ImplicitHydrogensExceedMax {
        element: Element,
        implicit_hydrogens: u8,
        max_implicit_hydrogens: u8,
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
    implicit_hydrogens: u8,
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
        implicit_hydrogens: u8,
        lone_pairs: u8,
        unpaired_electrons: u8,
        multiplicity: SpinMultiplicity,
        valence: u8,
        donated_pairs: u8,
        accepted_pairs: u8,
        aromatic_valence: AromaticValence,
        multicenter_valence: u8,
    ) -> Result<Self, AtomError> {
        let spin = SpinState::try_new(unpaired_electrons, multiplicity)?;
        Ok(Self {
            element,
            charge,
            implicit_hydrogens,
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

    pub fn implicit_hydrogens(&self) -> u8 {
        self.implicit_hydrogens
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

    pub fn check_invariants(&self) -> Result<(), AtomError> {
        let (min_charge, max_charge) = self.element.charge_bounds();
        if self.charge < min_charge || self.charge > max_charge {
            return Err(AtomError::ChargeOutOfBounds {
                element: self.element,
                charge: self.charge,
                min_charge,
                max_charge,
            });
        }

        let max_valence = self.element.max_valence();
        if self.valence > max_valence {
            return Err(AtomError::ValenceExceedsMax {
                element: self.element,
                valence: self.valence,
                max_valence,
            });
        }

        let unpaired_electrons = self.spin.unpaired_electrons();
        let max_unpaired_electrons = self.element.max_unpaired_electrons();
        if unpaired_electrons > max_unpaired_electrons {
            return Err(AtomError::UnpairedElectronsExceedMax {
                element: self.element,
                unpaired_electrons,
                max_unpaired_electrons,
            });
        }

        let max_implicit_hydrogens = self.element.max_implicit_hydrogens();
        if self.implicit_hydrogens > max_implicit_hydrogens {
            return Err(AtomError::ImplicitHydrogensExceedMax {
                element: self.element,
                implicit_hydrogens: self.implicit_hydrogens,
                max_implicit_hydrogens,
            });
        }

        let aromatic_valence = self.aromatic_valence.valence() as i16;
        let aromatic_increment = aromatic_increment(self.aromatic_valence) as i16;
        let total_e_inv_o = unpaired_electrons as i16
            + (2 * self.lone_pairs as i16)
            + (2 * self.donated_pairs as i16)
            + (2 * self.accepted_pairs as i16)
            + (2 * self.implicit_hydrogens as i16)
            + (2 * self.valence as i16)
            + aromatic_valence
            + aromatic_increment
            + (self.multicenter_valence as i16);

        let total_e_inv_e = (self.element.valence_electrons() as i16) - (self.charge as i16)
            + (self.implicit_hydrogens as i16)
            + (self.valence as i16)
            + aromatic_increment
            + (self.multicenter_valence as i16)
            + (2 * self.accepted_pairs as i16);

        if total_e_inv_o != total_e_inv_e {
            return Err(AtomError::ElectronInvariantMismatch {
                element: self.element,
                orbital_invariant: total_e_inv_o,
                electron_invariant: total_e_inv_e,
            });
        }

        Ok(())
    }
}

fn aromatic_increment(aromatic_valence: AromaticValence) -> u8 {
    match aromatic_valence {
        AromaticValence::None => 0,
        AromaticValence::Valence(0) => 0,
        AromaticValence::Valence(1) => 1,
        AromaticValence::Valence(2) => 0,
        AromaticValence::Valence(_) => 0,
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
        if self.implicit_hydrogens > 0 {
            if self.implicit_hydrogens == 1 {
                write!(f, "H")?;
            } else {
                write!(f, "H{}", self.implicit_hydrogens)?;
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
    type Err = AtomError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let trimmed = s.trim();
        if !trimmed.starts_with('{') || !trimmed.ends_with('}') {
            return Err(AtomError::InvalidSpecFormat);
        }
        let body = &trimmed[1..trimmed.len() - 1];
        let mut chars = body.chars().peekable();

        let first = chars.next().ok_or(AtomError::EmptySpec)?;
        if !first.is_ascii_uppercase() {
            return Err(AtomError::InvalidElement);
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
        let element: Element = elem.parse().map_err(|_| AtomError::InvalidElement)?;

        let mut charge = None;
        let mut implicit_hydrogens = 0u8;
        let mut lone_pairs = 0_u8;
        let mut multiplicity: Option<SpinMultiplicity> = None;
        let mut valence = 0u8;
        let mut donated_pairs = 0u8;
        let mut accepted_pairs = 0u8;
        let mut unpaired_electrons = 0u8;
        let mut aromatic_valence = AromaticValence::None;
        let mut multicenter_valence = 0u8;

        while let Some(token) = chars.next() {
            if token.is_ascii_whitespace() {
                continue;
            }
            let mut number = String::new();
            while chars.peek().is_some_and(|c| c.is_ascii_digit()) {
                number.push(chars.next().unwrap());
            }
            let num_u8 = |default: u8, error: AtomError| -> Result<u8, AtomError> {
                if number.is_empty() {
                    Ok(default)
                } else {
                    number.parse::<u8>().map_err(|_| error)
                }
            };
            match token {
                '+' => {
                    if charge.is_some() {
                        return Err(AtomError::DuplicateChargeToken);
                    }
                    charge = Some(num_u8(1, AtomError::InvalidCharge)? as i8);
                }
                '-' => {
                    if charge.is_some() {
                        return Err(AtomError::DuplicateChargeToken);
                    }
                    charge = Some(-(num_u8(1, AtomError::InvalidCharge)? as i8));
                }
                'H' => implicit_hydrogens = num_u8(1, AtomError::InvalidImplicitHydrogens)?,
                '/' => lone_pairs = num_u8(1, AtomError::InvalidLonePairs)?,
                '^' => unpaired_electrons = num_u8(1, AtomError::InvalidUnpairedElectrons)?,
                'x' => {
                    let m = num_u8(1, AtomError::InvalidMultiplicity)?;
                    multiplicity = Some(
                        SpinMultiplicity::from_multiplicity(m)
                            .ok_or_else(|| AtomError::InvalidMultiplicity)?,
                    );
                }
                'v' => valence = num_u8(1, AtomError::InvalidValence)?,
                '>' => donated_pairs = num_u8(1, AtomError::InvalidDonatedPairs)?,
                '<' => accepted_pairs = num_u8(1, AtomError::InvalidAcceptedPairs)?,
                'a' => {
                    aromatic_valence =
                        AromaticValence::Valence(num_u8(1, AtomError::InvalidAromaticValence)?)
                }
                'm' => multicenter_valence = num_u8(1, AtomError::InvalidMulticenterValence)?,
                _ => {
                    return Err(AtomError::UnexpectedToken { token });
                }
            }
        }

        if unpaired_electrons > MAX_UNPAIRED_ELECTRONS {
            return Err(AtomError::UnpairedElectronsLiteralExceedMax {
                unpaired_electrons,
                max_unpaired_electrons: MAX_UNPAIRED_ELECTRONS,
            });
        }

        let multiplicity = match multiplicity {
            Some(m) => m,
            None => SpinState::max_multiplicity(unpaired_electrons)
                .ok_or(SpinStateError::Underdetermined)?
                .multiplicity(),
        };

        Self::new(
            element,
            charge.unwrap_or(0),
            implicit_hydrogens,
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
        let hydrogen_constraint = atom
            .implicit_hydrogens()
            .map(HydrogenConstraint::from_implicit_hydrogens);
        let aromatic_constraint = if builder.atom_aromatic_hint(atom_index) {
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
            implicit_hydrogens: hydrogen_constraint,
            lone_pairs: atom.lone_pairs(),
            unpaired_electrons: atom.unpaired_electrons(),
            multiplicity: atom.multiplicity(),
            valence: Some(valence),
            donated_pairs: Some(donated_pairs),
            accepted_pairs: Some(accepted_pairs),
            aromatic_valence: aromatic_constraint,
            multicenter_valence,
        }
    }

    pub fn matches(&self, spec: &AtomTypeSpec) -> bool {
        self.charge.is_none_or(|v| v == spec.charge())
            && self
                .implicit_hydrogens
                .is_none_or(|v| v.matches(spec.implicit_hydrogens()))
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
            return Err(AtomError::InvalidElement);
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
        let element: Element = elem.parse().map_err(|_| AtomError::InvalidElement)?;

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
            let num_u8 = |default: u8, error: AtomError| -> Result<u8, AtomError> {
                if number.is_empty() {
                    Ok(default)
                } else {
                    number.parse::<u8>().map_err(|_| error)
                }
            };
            match token {
                '+' => {
                    if seen_charge {
                        return Err(AtomError::DuplicateChargeToken);
                    }
                    query.charge = Some(num_u8(1, AtomError::InvalidCharge)? as i8);
                    seen_charge = true;
                }
                '-' => {
                    if seen_charge {
                        return Err(AtomError::DuplicateChargeToken);
                    }
                    query.charge = Some(-(num_u8(1, AtomError::InvalidCharge)? as i8));
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
                        query.implicit_hydrogens = Some(HydrogenConstraint::Hydrogens(num_u8(
                            1,
                            AtomError::InvalidImplicitHydrogens,
                        )?));
                    }
                }
                '/' => query.lone_pairs = Some(num_u8(1, AtomError::InvalidLonePairs)?),
                '^' => {
                    query.unpaired_electrons = Some(num_u8(1, AtomError::InvalidUnpairedElectrons)?)
                }
                'x' => {
                    let m = num_u8(1, AtomError::InvalidMultiplicity)?;
                    query.multiplicity = Some(
                        SpinMultiplicity::from_multiplicity(m)
                            .ok_or_else(|| AtomError::InvalidMultiplicity)?,
                    );
                }
                'v' => query.valence = Some(num_u8(1, AtomError::InvalidValence)?),
                '>' => query.donated_pairs = Some(num_u8(1, AtomError::InvalidDonatedPairs)?),
                '<' => query.accepted_pairs = Some(num_u8(1, AtomError::InvalidAcceptedPairs)?),
                'a' => {
                    if chars.peek() == Some(&'*') {
                        chars.next();
                        query.aromatic_valence = Some(AromaticConstraint::Any);
                    } else if chars.peek() == Some(&'!') {
                        chars.next();
                        query.aromatic_valence = Some(AromaticConstraint::None);
                    } else {
                        query.aromatic_valence = Some(AromaticConstraint::Valence(num_u8(
                            1,
                            AtomError::InvalidAromaticValence,
                        )?));
                    }
                }
                'm' => {
                    query.multicenter_valence =
                        Some(num_u8(1, AtomError::InvalidMulticenterValence)?)
                }
                _ => {
                    return Err(AtomError::UnexpectedToken { token });
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
    use umol_data::{Element, SpinStateError};

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
    #[case::outer_spaced("  {N}  ", Element::N, 0, 0, 0, 0, SpinMultiplicity::Singlet, 0, 0, 0, AromaticValence::None, 0)]
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
        assert_eq!(spec.implicit_hydrogens(), hydrogens, "hydrogens mismatch for {}", input);
        assert_eq!(spec.lone_pairs(), lone_pairs, "lone pairs mismatch for {}", input);
        assert_eq!(spec.unpaired_electrons(), unpaired_electrons, "unpaired electrons mismatch for {}", input);
        assert_eq!(spec.multiplicity(), multiplicity, "multiplicity mismatch for {}", input);
        assert_eq!(spec.valence(), valence, "valence mismatch for {}", input);
        assert_eq!(spec.donated_pairs(), donated_pairs, "donated pairs mismatch for {}", input);
        assert_eq!(spec.accepted_pairs(), accepted_pairs, "accepted pairs mismatch for {}", input);
        assert_eq!(spec.aromatic_valence(), aromatic_valence, "aromatic valence mismatch for {}", input);
        assert_eq!(spec.multicenter_valence(), multicenter_valence, "multicenter valence mismatch for {}", input);
    }

    #[rustfmt::skip]
    #[rstest]
    #[case::invalid_format("N", AtomError::InvalidSpecFormat)]
    #[case::empty_spec("{}", AtomError::EmptySpec)]
    #[case::invalid_element("{n}", AtomError::InvalidElement)]
    #[case::duplicate_charge("{N+-}", AtomError::DuplicateChargeToken)]
    #[case::invalid_multiplicity("{Nx11}", AtomError::InvalidMultiplicity)]
    #[case::invalid_spin("{N^2x2}", AtomError::SpinState(SpinStateError::Incompatible {
        unpaired_electrons: 2,
        multiplicity: SpinMultiplicity::Doublet,
    }))]
    #[case::unexpected_token("{Nq1}", AtomError::UnexpectedToken { token: 'q' })]
    fn test_atom_type_spec_from_str_error(#[case] input: &str, #[case] expected: AtomError) {
        assert_eq!(AtomTypeSpec::from_str(input).unwrap_err(), expected);
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
    #[case::valid("{C+0v4a0m0}", None)]
    #[case::charge_out_of_bounds("{C+5}", Some(AtomError::ChargeOutOfBounds { element: Element::C, charge: 5, min_charge: -4, max_charge: 1 }))]
    #[case::valence_exceeds_max("{Hv2}", Some(AtomError::ValenceExceedsMax { element: Element::H, valence: 2, max_valence: 1, }))]
    #[case::unpaired_exceeds_max("{O^3x2}", Some(AtomError::UnpairedElectronsExceedMax { element: Element::O, unpaired_electrons: 3, max_unpaired_electrons: 2 }))]
    #[case::implicit_hydrogens_exceed_max("{OH4}", Some(AtomError::ImplicitHydrogensExceedMax { element: Element::O, implicit_hydrogens: 4, max_implicit_hydrogens: 3 }))]
    #[case::electron_invariant_mismatch("{Cv1}", Some(AtomError::ElectronInvariantMismatch { element: Element::C, orbital_invariant: 2, electron_invariant: 5 }))]
    fn test_atom_type_spec_check_invariants(
        #[case] input: &str,
        #[case] expected: Option<AtomError>,
    ) {
        let spec = AtomTypeSpec::from_str(input).unwrap();
        assert_eq!(spec.check_invariants().err(), expected);
    }

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
    #[case::invalid_element("?{c}", AtomError::InvalidElement)]
    #[case::duplicate_charge("?{C+-}", AtomError::DuplicateChargeToken)]
    #[case::invalid_multiplicity("?{Cx11}", AtomError::InvalidMultiplicity)]
    #[case::unexpected_token("?{Cq1}", AtomError::UnexpectedToken { token: 'q' })]
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
