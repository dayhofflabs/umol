//! Valence states of atoms (for valence graph validation)
//!
//! Defines the `ValenceState` struct used for strict atom typing in molecular graphs.
//! This module also implicitly defines an internal string notation for these states,
//! primarily intended for easy definition within data files or code, not for general exchange.
//!
//! ## Internal Notation Format
//!
//! The notation resembles SMARTS atom primitives but with extensions for lone pairs (`/`)
//! and optional multiplicity (`*`), supporting arbitrary numbers of unpaired electrons (`^`).
//! The general format is: `[Element<Charge><LonePairs><Unpaired><Multiplicity><Valence>]`
//!
//! - **Element:** Standard atomic symbol (e.g., `C`, `Fe`). Required.
//! - **Properties:** Charge, Lone Pairs, Unpaired Electrons, Multiplicity, and Valence.
//!   These can appear in **any order** within the brackets after the element.
//!   If a property is omitted, a default value is assumed (see below).
//!
//! ### Property Details:
//! - **Charge:** `[+-]\d*` (e.g., `+1`, `-2`, `+`). If sign is present without a number (`+`, `-`), assumes `1`. Defaults to `0` if omitted.
//! - **Lone Pairs:** `/\d*` (e.g., `/1`, `/0`). If `/` is present without a number, assumes `1`. Defaults to `0` if omitted.
//! - **Unpaired Electrons:** `\^\d*` (e.g., `^2`, `^0`). If `^` is present without a number, assumes `1`. Defaults to `0` if omitted.
//! - **Multiplicity:** `\*\d+` (e.g., `*1`, `*3`). Must include a number if present. Defaults to (`unpaired + 1`, or `1` if unpaired is `0`) if omitted.
//! - **Valence:** `v\d+` (e.g., `v4`, `v6`). Must include a number if present. Defaults to `0` if omitted (use with caution, usually should be specified).
//!
//! ### Examples (Illustrative for Parsing):
//! - `[C]` -> C, charge=0, lp=0, unpaired=0, mult=1, valence=0
//! - `[C+v4]` -> C, charge=+1, lp=0, unpaired=0, mult=1, valence=4
//! - `[Fe+2/1^4*5v6]` -> Fe, charge=+2, lp=1, unpaired=4, mult=5, valence=6
//! - `[O-/1^1]` -> O, charge=-1, lp=1, unpaired=1, mult=2 (default), valence=0

use crate::Element;
use regex::Regex;
use std::fmt::{self, Display};
use std::str::FromStr;
use umol::error::DataError;
use umol::{Error, Result};

/// Valence state describes the fictitious state of an atom involved in bond formation
/// Defined by element, charge, number of lone pairs, unpaired electrons, and valence
/// Does not uniquely define an atomic term symbol
/// Valence is used to validate ValenceAtom objects according to:
/// `valence = unpaired_electrons + sum(bond_orders) + num_implicit_hydrogens`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ValenceState {
    element: Element,
    charge: i8,
    lone_pairs: u8,
    unpaired_electrons: u8,
    multiplicity: u8,
    valence: u8,
}

impl ValenceState {
    /// Creates a new ValenceState.
    pub fn new(
        element: Element,
        charge: i8,
        lone_pairs: u8,
        unpaired_electrons: u8,
        multiplicity: u8,
        valence: u8,
    ) -> Self {
        debug_assert!(multiplicity >= 1 && multiplicity <= unpaired_electrons + 1);
        debug_assert!((unpaired_electrons + 1 - multiplicity) % 2 == 0);
        debug_assert!(element.valence_electrons() as i8 - charge >= 0);
        Self {
            element,
            charge,
            lone_pairs,
            unpaired_electrons,
            multiplicity,
            valence,
        }
    }

    pub fn element(&self) -> Element {
        self.element
    }

    pub fn charge(&self) -> i8 {
        self.charge
    }

    pub fn lone_pairs(&self) -> u8 {
        self.lone_pairs
    }

    pub fn unpaired_electrons(&self) -> u8 {
        self.unpaired_electrons
    }

    pub fn electron_count(&self) -> u8 {
        2 * self.lone_pairs + self.unpaired_electrons
    }

    pub fn multiplicity(&self) -> u8 {
        self.multiplicity
    }

    /// Target valence number used in the check formula:
    /// `valence = unpaired_electrons + sum(bond_orders) + num_implicit_hydrogens`
    pub fn valence(&self) -> u8 {
        self.valence
    }

    /// Calculates the default multiplicity (unpaired + 1)
    pub fn default_multiplicity(&self) -> u8 {
        self.unpaired_electrons + 1
    }
}

impl Display for ValenceState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "[{}", self.element)?;
        match self.charge {
            0 => (),
            1 => write!(f, "+")?,
            -1 => write!(f, "-")?,
            c if c < 0 => write!(f, "-{}", -c)?,
            c if c > 0 => write!(f, "+{}", c)?,
            _ => unreachable!(),
        };
        if self.lone_pairs > 0 {
            write!(f, "/{}", self.lone_pairs)?;
        }
        if self.unpaired_electrons > 0 {
            write!(f, "^{}", self.unpaired_electrons)?;
        }
        // Only display multiplicity if it's non-default
        if self.multiplicity != self.default_multiplicity() {
            write!(f, "*{}", self.multiplicity)?;
        }
        if self.valence > 0 {
            write!(f, "v{}", self.valence)?;
        }
        write!(f, "]")
    }
}

impl FromStr for ValenceState {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        if s.is_empty() {
            return Err(DataError::InvalidValenceState("Empty string".to_string()).into());
        }

        let val_state_pattern = Regex::new(r"^\[([A-Z][a-z]?)((?:[-+/*v^]\d*)*)]$").unwrap();
        if !val_state_pattern.is_match(s) {
            return Err(
                DataError::InvalidValenceState(format!("Invalid valence state: {}", s)).into(),
            );
        }

        let caps = val_state_pattern.captures(s).unwrap();
        let element = caps[1].parse::<Element>()?;
        let properties = caps.get(2);
        if properties.is_none() {
            return Ok(ValenceState::new(element, 0, 0, 0, 1, 0));
        }
        let properties = properties.unwrap().as_str();

        let property_pattern = Regex::new(r"([-+/*v^])(\d*)").unwrap();
        let mut pos_charge = None;
        let mut neg_charge = None;
        let mut lone_pairs = None;
        let mut unpaired = None;
        let mut multiplicity = None;
        let mut valence = None;

        for cap in property_pattern.captures_iter(properties) {
            let property = &cap[1];
            let value = &cap[2];
            match property {
                "+" => {
                    if pos_charge.is_some() || neg_charge.is_some() {
                        return Err(
                            DataError::InvalidValenceState("Duplicate charge".to_string()).into(),
                        );
                    } else {
                        pos_charge = if value.is_empty() {
                            Some(1)
                        } else {
                            value.parse::<i8>().ok()
                        };
                    }
                }
                "-" => {
                    if neg_charge.is_some() || pos_charge.is_some() {
                        return Err(
                            DataError::InvalidValenceState("Duplicate charge".to_string()).into(),
                        );
                    } else {
                        neg_charge = if value.is_empty() {
                            Some(-1)
                        } else {
                            value.parse::<i8>().ok()
                        };
                    }
                }
                "/" => {
                    if lone_pairs.is_some() {
                        return Err(DataError::InvalidValenceState(
                            "Duplicate lone pair".to_string(),
                        )
                        .into());
                    } else {
                        lone_pairs = if value.is_empty() {
                            Some(1)
                        } else {
                            value.parse::<u8>().ok()
                        };
                    }
                }
                "^" => {
                    if unpaired.is_some() {
                        return Err(DataError::InvalidValenceState(
                            "Duplicate unpaired".to_string(),
                        )
                        .into());
                    } else {
                        unpaired = if value.is_empty() {
                            Some(1)
                        } else {
                            value.parse::<u8>().ok()
                        };
                    }
                }
                "*" => {
                    multiplicity = if multiplicity.is_some() {
                        return Err(DataError::InvalidValenceState(
                            "Duplicate multiplicity".to_string(),
                        )
                        .into());
                    } else {
                        if value.is_empty() {
                            Some(1)
                        } else {
                            value.parse::<u8>().ok()
                        }
                    }
                }
                "v" => {
                    valence = if valence.is_some() {
                        return Err(DataError::InvalidValenceState(
                            "Duplicate valence".to_string(),
                        )
                        .into());
                    } else {
                        if value.is_empty() {
                            Some(0)
                        } else {
                            value.parse::<u8>().ok()
                        }
                    }
                }
                _ => unreachable!(),
            }
        }

        let charge = pos_charge.unwrap_or(neg_charge.unwrap_or(0));
        let lone_pairs = lone_pairs.unwrap_or(0);
        let unpaired = unpaired.unwrap_or(0);
        // Multiplcity defaults to unpaired + 1
        let multiplicity = multiplicity.unwrap_or(unpaired + 1);
        let valence = valence.unwrap_or(0);

        Ok(ValenceState::new(
            element,
            charge,
            lone_pairs,
            unpaired,
            multiplicity,
            valence,
        ))
    }
}

/// Shorthand macro for valence state parsing
#[macro_export]
macro_rules! vs {
    ($s:expr) => {
        $s.parse::<ValenceState>().unwrap()
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::e;
    use rstest::*;

    #[rstest]
    #[case(e!(H), 0, 0, 1, 2, 1, "[H^1v1]")] // H radical
    #[case(e!(C), 0, 0, 0, 1, 4, "[Cv4]")] // Standard C
    #[case(e!(C), 0, 1, 0, 1, 2, "[C/1v2]")] // Singlet carbene C
    #[case(e!(C), 0, 1, 2, 3, 2, "[C/1^2v2]")] // Triplet carbene C
    #[case(e!(C), 1, 0, 0, 1, 3, "[C+v3]")] // Carbenium ion
    #[case(e!(C), -1, 1, 0, 1, 3, "[C-/1v3]")] // Carbanion
    #[case(e!(P), 0, 1, 0, 1, 3, "[P/1v3]")] // Trivalent
    #[case(e!(P), 0, 0, 0, 1, 5, "[Pv5]")] // Pentavalent P
    #[case(e!(Gd), 3, 0, 7, 8, 0, "[Gd+3^7]")] // Gd(+3) ion
    #[case(e!(N), 0, 1, 3, 2, 0, "[N/1^3*2]")] // N atom doublet state
    fn test_valence_state_display(
        #[case] element: Element,
        #[case] charge: i8,
        #[case] lone_pairs: u8,
        #[case] unpaired_electrons: u8,
        #[case] multiplicity: u8,
        #[case] valence: u8,
        #[case] expected: &str,
    ) {
        let vs = ValenceState::new(
            element,
            charge,
            lone_pairs,
            unpaired_electrons,
            multiplicity,
            valence,
        );
        assert_eq!(format!("{}", vs), expected);
    }

    #[rstest]
    #[case("[C]", ValenceState::new(e!(C), 0, 0, 0, 1, 0))] // No properties
    #[case("[Li+]", ValenceState::new(e!(Li), 1, 0, 0, 1, 0))] // Only positive charge
    #[case("[F-]", ValenceState::new(e!(F), -1, 0, 0, 1, 0))] // Only negative charge
    #[case("[C/1]", ValenceState::new(e!(C), 0, 1, 0, 1, 0))] // Only lone pair
    #[case("[Na^1]", ValenceState::new(e!(Na), 0, 0, 1, 2, 0))] // Unpaired, default multiplicity
    #[case("[N^3*4]", ValenceState::new(e!(N), 0, 0, 3, 4, 0))] // Unpaired, explicit multiplicity
    #[case("[Cv4]", ValenceState::new(e!(C), 0, 0, 0, 1, 4))] // Only valence
    #[case("[C+v3]", ValenceState::new(e!(C), 1, 0, 0, 1, 3))] // Charge +1 default
    #[case("[O-/2^1]", ValenceState::new(e!(O), -1, 2, 1, 2, 0))] // Charge -1 default, lp 2, unpaired 1, mult default 2
    #[case("[Fe+2/1^4*5v6]", ValenceState::new(e!(Fe), 2, 1, 4, 5, 6))] // Full spec
    #[case("[N*2/1^3]", ValenceState::new(e!(N), 0, 1, 3, 2, 0))] // Arbitrary order, explicit mult
    #[case("[P^v5]", ValenceState::new(e!(P), 0, 0, 1, 2, 5))] // Unpaired default 1, Valence 5
    #[case("[S/v3+]", ValenceState::new(e!(S), 1, 1, 0, 1, 3))] // Arbitrary order, lp default 1
    #[case("[Gd+3^7]", ValenceState::new(e!(Gd), 3, 0, 7, 8, 0))] // High unpaired, default mult
    #[case("[C^2v2]", ValenceState::new(e!(C), 0, 0, 2, 3, 2))] // Triplet carbene
    #[case("[C^2*1v2]", ValenceState::new(e!(C), 0, 0, 2, 1, 2))] // Singlet carbene
    #[case("[C+/0^0v3]", ValenceState::new(e!(C), 1, 0, 0, 1, 3))] // Explicit charge, lp, unpaired
    fn test_from_str_valid(#[case] s: &str, #[case] expected: ValenceState) {
        assert_eq!(ValenceState::from_str(s).unwrap(), expected);
    }

    #[rstest]
    // Invalid cases
    #[case("")] // Empty
    #[case("C")] // No brackets
    #[case("[C")] // Mismatched bracket
    #[case("[]")] // No element
    #[case("[Xx]")] // Invalid element
    #[case("[C+1+1]")] // Duplicate charge
    #[case("[C/1/1]")] // Duplicate lone pair
    #[case("[C^1^1]")] // Duplicate unpaired
    #[case("[C*1*1]")] // Duplicate multiplicity
    #[case("[Cv1v1]")] // Duplicate valence
    #[case("[C+?]")] // Invalid char in properties
    #[case("[C+1a]")] // Invalid char after property
    #[case("[C+1.0]")] // Non-integer
    fn test_from_str_invalid(#[case] s: &str) {
        assert!(ValenceState::from_str(s).is_err());
    }

    #[test]
    fn test_valence_state_macro() {
        assert_eq!(vs!("[C+0/1^2*3v4]"), ValenceState::new(e!(C), 0, 1, 2, 3, 4));
    }
}
