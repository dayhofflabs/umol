//! Atom spec
//!
//! Defines atom specs for strictly typed molecular valence graphs and an internal notation,
//! primarily intended for easy definition within data files or code, not for general exchange.
//!
//! ## Internal Notation Format for atom specs
//!
//! The notation resembles SMARTS atom primitives but with extensions for lone pairs (`/`),
//! donated pairs (`>`) and accepted pairs (`<`), optional multiplicity (`*`),
//! and supporting arbitrary numbers of unpaired electrons (`^`).
//! The format is: `[Element<Charge><LonePairs><Donated><Accepted><Unpaired><Multiplicity><ImplicitHs><Valence>]`
//!
//! - **Element:** Standard atomic symbol (e.g., `C`, `Fe`). Required.
//! - **Properties:** Charge, Lone Pairs, Donated Pairs, Accepted Pairs, Unpaired Electrons,
//!   Multiplicity, Implicit Hs, and Valence.
//!   These can appear in any order within the brackets after the element.
//!   If a property is omitted, a default value is assumed (see below).
//!   A property symbol without a number is assumed to be `1`.
//!
//! Note that all property definitions must be explicit, in contrast to the SMARTS notation
//! semantics. Specifically, `[CH4]` defines valence = 0 and is thus invalid (because the number
//! of implicit Hs > valence), unlike SMARTS. The correct notation is `[CH4v4]`.
//!
//! ### Property Details
//! - **Charge:** `[+-]\d*` (e.g., `+1`, `-2`, `+0`, `+`). Defaults to `0` if omitted.
//! - **Lone Pairs:** `/\d*` (e.g., `/1`, `/2`, `/0`, `/`). Defaults to `0` if omitted.
//! - **Donated Pairs:** `>\d*` (e.g., `>1`, `>0`, `>`). Defaults to `0` if omitted.
//! - **Accepted Pairs:** `<` (e.g., `<1`, `<0`, `<`). Defaults to `0` if omitted.
//! - **Unpaired Electrons:** `\^\d*` (e.g., `^2`, `^0`, `^`). Defaults to `0` if omitted.
//! - **Multiplicity:** `\*\d+` (e.g., `*1`, `*3`). Defaults to `unpaired + 1` if omitted.
//! - **Implicit Hs:** `H\d+` (e.g., `H4`, `H6`). Defaults to `0` if omitted.
//! - **Valence:** `v\d+` (e.g., `v4`, `v6`). Defaults to `0` if omitted.
//!
//! ### Examples
//! - `[C]` -> C, charge=0, lp=0, unpaired=0, mult=1, implicit_hydrogens=0, valence=0
//! - `[CH4v4]` -> C, charge=0, lp=0, unpaired=0, mult=1, implicit_hydrogens=4, valence=4
//! - `[C+v4]` -> C, charge=+1, lp=0, unpaired=0, mult=1, implicit_hydrogens=0, valence=4
//! - `[N+0/1>1v3]` -> N, charge=+0, lp=1, donated=1, unpaired=0, mult=2, implicit_hydrogens=0, valence=3
//! - `[Fe+2/1<4^4*5v6]` -> Fe, charge=+2, lp=1, donated=4, unpaired=4, mult=5, valence=6
//! - `[O-/1^1]` -> O, charge=-1, lp=1, unpaired=1, mult=2 (default), valence=0

use std::fmt::{self, Display};
use std::str::FromStr;

use regex::Regex;
use serde::{Deserialize, Serialize};
use umol::error::DataError;
use umol::{Error, Result};
use umol_data::Element;

/// AtomSpec is a generalization of the atomic valence state for modeling
/// covalent and dative bonding by molecular graphs.
/// Does not uniquely define an atomic term symbol for quantum chemical calculations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct AtomSpec {
    element: Element,
    charge: i8,
    lone_pairs: u8,
    donated_pairs: u8,
    accepted_pairs: u8,
    unpaired_electrons: u8,
    multiplicity: u8,
    implicit_hydrogens: u8,
    valence: u8,
}

impl AtomSpec {
    pub fn new(
        element: Element,
        charge: i8,
        lone_pairs: u8,
        donated_pairs: u8,
        accepted_pairs: u8,
        unpaired_electrons: u8,
        multiplicity: u8,
        implicit_hydrogens: u8,
        valence: u8,
    ) -> Self {
        debug_assert!(
            multiplicity >= 1 && multiplicity <= unpaired_electrons + 1,
            "Multiplicity must be between 1 and unpaired electrons + 1"
        );
        debug_assert!(
            (unpaired_electrons + 1 - multiplicity).is_multiple_of(2),
            "Multiplicity must be even"
        );
        debug_assert!(
            donated_pairs <= lone_pairs,
            "Donated pairs must be less than or equal to lone pairs"
        );
        debug_assert!(
            (element.valence_electrons() as i8 - charge) >= 0,
            "Charge must be less than or equal to valence electrons"
        );
        debug_assert!(
            implicit_hydrogens <= element.max_implicit_hydrogens(),
            "Implicit hydrogens must be less than or equal to max implicit hydrogens"
        );
        debug_assert!(
            implicit_hydrogens <= valence,
            "Implicit hydrogens must be less than or equal to valence"
        );
        debug_assert!(
            valence <= element.max_valence(),
            "Valence must be less than or equal to max valence"
        );
        Self {
            element,
            charge,
            lone_pairs,
            donated_pairs,
            accepted_pairs,
            unpaired_electrons,
            multiplicity,
            implicit_hydrogens,
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

    pub fn donated_pairs(&self) -> u8 {
        self.donated_pairs
    }

    pub fn accepted_pairs(&self) -> u8 {
        self.accepted_pairs
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

    pub fn implicit_hydrogens(&self) -> u8 {
        self.implicit_hydrogens
    }

    pub fn valence(&self) -> u8 {
        self.valence
    }
}

impl Display for AtomSpec {
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
            write!(f, "/")?;
            if self.lone_pairs > 1 {
                write!(f, "{}", self.lone_pairs)?;
            }
        }
        if self.donated_pairs > 0 {
            write!(f, ">")?;
            if self.donated_pairs > 1 {
                write!(f, "{}", self.donated_pairs)?;
            }
        }
        if self.accepted_pairs > 0 {
            write!(f, "<")?;
            if self.accepted_pairs > 1 {
                write!(f, "{}", self.accepted_pairs)?;
            }
        }
        if self.unpaired_electrons > 0 {
            write!(f, "^")?;
            if self.unpaired_electrons > 1 {
                write!(f, "{}", self.unpaired_electrons)?;
            }
        }
        // Only display multiplicity if it's non-default
        if self.multiplicity != self.unpaired_electrons + 1 {
            write!(f, "*{}", self.multiplicity)?;
        }

        if self.implicit_hydrogens > 0 {
            write!(f, "H{}", self.implicit_hydrogens)?;
        }

        if self.valence > 0 {
            write!(f, "v")?;
            if self.valence > 1 {
                write!(f, "{}", self.valence)?;
            }
        }
        write!(f, "]")
    }
}

impl TryFrom<&str> for AtomSpec {
    type Error = Error;

    fn try_from(s: &str) -> Result<Self> {
        if s.is_empty() {
            return Err(DataError::InvalidAtomSpec("Empty atom spec".to_string()).into());
        }

        let val_state_pattern = Regex::new(r"^\[([A-Z][a-z]?)((?:[-+/<>*Hv^]\d*)*)]$").unwrap();
        if !val_state_pattern.is_match(s) {
            return Err(DataError::InvalidAtomSpec(s.to_string()).into());
        }

        let caps = val_state_pattern.captures(s).unwrap();
        let element = caps[1].parse::<Element>()?;
        let properties = caps.get(2);
        if properties.is_none() {
            return Ok(AtomSpec::new(element, 0, 0, 0, 0, 0, 1, 0, 0));
        }
        let properties = properties.unwrap().as_str();

        let property_pattern = Regex::new(r"([-+/<>*Hv^])(\d*)").unwrap();
        let mut pos_charge = None;
        let mut neg_charge = None;
        let mut lone_pairs = None;
        let mut donated = None;
        let mut accepted = None;
        let mut unpaired = None;
        let mut multiplicity = None;
        let mut implicit_hydrogens = None;
        let mut valence = None;

        for cap in property_pattern.captures_iter(properties) {
            let property = &cap[1];
            let value = &cap[2];
            match property {
                "+" => {
                    if pos_charge.is_some() || neg_charge.is_some() {
                        return Err(DataError::InvalidAtomSpec(
                            "Duplicate charge definition".to_string(),
                        )
                        .into());
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
                        return Err(DataError::InvalidAtomSpec(
                            "Duplicate charge definition".to_string(),
                        )
                        .into());
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
                        return Err(DataError::InvalidAtomSpec(
                            "Duplicate lone pair definition".to_string(),
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
                ">" => {
                    if donated.is_some() {
                        return Err(DataError::InvalidAtomSpec(
                            "Duplicate donated pair definition".to_string(),
                        )
                        .into());
                    } else {
                        donated = if value.is_empty() {
                            Some(1)
                        } else {
                            value.parse::<u8>().ok()
                        };
                    }
                }
                "<" => {
                    if accepted.is_some() {
                        return Err(DataError::InvalidAtomSpec(
                            "Duplicate accepted pair definition".to_string(),
                        )
                        .into());
                    } else {
                        accepted = if value.is_empty() {
                            Some(1)
                        } else {
                            value.parse::<u8>().ok()
                        };
                    }
                }
                "^" => {
                    if unpaired.is_some() {
                        return Err(DataError::InvalidAtomSpec(
                            "Duplicate unpaired electron specification".to_string(),
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
                        return Err(DataError::InvalidAtomSpec(
                            "Duplicate multiplicity".to_string(),
                        )
                        .into());
                    } else if value.is_empty() {
                        Some(1)
                    } else {
                        value.parse::<u8>().ok()
                    }
                }
                "H" => {
                    if implicit_hydrogens.is_some() {
                        return Err(DataError::InvalidAtomSpec(
                            "Duplicate implicit hydrogen specification".to_string(),
                        )
                        .into());
                    } else {
                        implicit_hydrogens = if value.is_empty() {
                            Some(1)
                        } else {
                            value.parse::<u8>().ok()
                        };
                    }
                }
                "v" => {
                    valence = if valence.is_some() {
                        return Err(
                            DataError::InvalidAtomSpec("Duplicate valence".to_string()).into()
                        );
                    } else if value.is_empty() {
                        Some(0)
                    } else {
                        value.parse::<u8>().ok()
                    }
                }
                _ => unreachable!(),
            }
        }

        // Set default values
        let charge = pos_charge.unwrap_or(neg_charge.unwrap_or(0));
        let lone_pairs = lone_pairs.unwrap_or(0);
        let donated = donated.unwrap_or(0);
        let accepted = accepted.unwrap_or(0);
        let unpaired = unpaired.unwrap_or(0);
        // Multiplcity defaults to unpaired + 1
        let multiplicity = multiplicity.unwrap_or(unpaired + 1);
        let implicit_hydrogens = implicit_hydrogens.unwrap_or(0);
        let valence = valence.unwrap_or(0);

        // Check constraints
        if multiplicity == 0
            || multiplicity > unpaired + 1
            || (unpaired + 1 - multiplicity) % 2 != 0
        {
            return Err(DataError::InvalidAtomMultiplicity(format!("{}", multiplicity)).into());
        }
        if donated > lone_pairs {
            return Err(DataError::InvalidAtomDonatedPairs(format!("{}", donated)).into());
        }
        if (element.valence_electrons() as i8 - charge) < 0 {
            return Err(DataError::InvalidAtomCharge(format!("{}", charge)).into());
        }
        if implicit_hydrogens > element.max_implicit_hydrogens() {
            return Err(
                DataError::InvalidAtomImplicitHydrogens(format!("{}", implicit_hydrogens)).into(),
            );
        }
        if implicit_hydrogens > valence {
            return Err(DataError::InvalidAtomImplicitHydrogens(format!("{}", valence)).into());
        }
        if valence > element.max_valence() {
            return Err(DataError::InvalidAtomValence(format!("{}", valence)).into());
        }

        Ok(AtomSpec::new(
            element,
            charge,
            lone_pairs,
            donated,
            accepted,
            unpaired,
            multiplicity,
            implicit_hydrogens,
            valence,
        ))
    }
}

impl FromStr for AtomSpec {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        Self::try_from(s)
    }
}

/// Shorthand macro for atom spec parsing
#[macro_export]
macro_rules! a {
    ($s:expr) => {
        $s.parse::<AtomSpec>().unwrap()
    };
}

#[cfg(test)]
mod tests {
    use rstest::*;
    use umol_data::e;

    use super::*;

    #[rstest]
    #[case(e!(H), 0, 0, 0, 0, 1, 2, 0, 1, "[H^v]")] // H radical
    #[case(e!(C), 0, 0, 0, 0, 0, 1, 0, 4, "[Cv4]")] // Standard C
    #[case(e!(C), 0, 0, 0, 0, 0, 1, 4, 4, "[CH4v4]")] // Standard C with 4 implicit Hs
    #[case(e!(C), 0, 1, 0, 0, 0, 1, 0, 2, "[C/v2]")] // Singlet carbene C
    #[case(e!(C), 0, 1, 0, 0, 2, 3, 0, 2, "[C/^2v2]")] // Triplet carbene C
    #[case(e!(C), 1, 0, 0, 0, 0, 1, 0, 3, "[C+v3]")] // Carbenium ion
    #[case(e!(C), -1, 1, 0, 0, 0, 1, 0, 3, "[C-/v3]")] // Carbanion
    #[case(e!(P), 0, 1, 0, 0, 0, 1, 0, 3, "[P/v3]")] // Trivalent
    #[case(e!(P), 0, 0, 0, 0, 0, 1, 0, 5, "[Pv5]")] // Pentavalent P
    #[case(e!(Gd), 3, 0, 0, 0, 7, 8, 0, 0, "[Gd+3^7]")] // Gd(+3) ion
    #[case(e!(N), 0, 1, 0, 0, 3, 2, 0, 0, "[N/^3*2]")] // N atom doublet state
    fn test_atom_spec_display(
        #[case] element: Element,
        #[case] charge: i8,
        #[case] lone_pairs: u8,
        #[case] donated_pairs: u8,
        #[case] accepted_pairs: u8,
        #[case] unpaired_electrons: u8,
        #[case] multiplicity: u8,
        #[case] implicit_hydrogens: u8,
        #[case] valence: u8,
        #[case] expected: &str,
    ) {
        let vs = AtomSpec::new(
            element,
            charge,
            lone_pairs,
            donated_pairs,
            accepted_pairs,
            unpaired_electrons,
            multiplicity,
            implicit_hydrogens,
            valence,
        );
        assert_eq!(format!("{}", vs), expected);
    }

    #[rstest]
    #[case("[C]", AtomSpec::new(e!(C), 0, 0, 0, 0, 0, 1, 0, 0))] // No properties
    #[case("[Li+]", AtomSpec::new(e!(Li), 1, 0, 0, 0, 0, 1, 0, 0))] // Only positive charge
    #[case("[F-]", AtomSpec::new(e!(F), -1, 0, 0, 0, 0, 1, 0, 0))] // Only negative charge
    #[case("[C/1]", AtomSpec::new(e!(C), 0, 1, 0, 0, 0, 1, 0, 0))] // Only lone pair
    #[case("[Na^1]", AtomSpec::new(e!(Na), 0, 0, 0, 0, 1, 2, 0, 0))] // Unpaired, default multiplicity
    #[case("[N^3*4]", AtomSpec::new(e!(N), 0, 0, 0, 0, 3, 4, 0, 0))] // Unpaired, explicit multiplicity
    #[case("[Cv4]", AtomSpec::new(e!(C), 0, 0, 0, 0, 0, 1, 0, 4))] // Only valence
    #[case("[C+v3]", AtomSpec::new(e!(C), 1, 0, 0, 0, 0, 1, 0, 3))] // Charge +1 default
    #[case("[O-/2^1]", AtomSpec::new(e!(O), -1, 2, 0, 0, 1, 2, 0, 0))] // Charge -1 default, lp 2, unpaired 1, mult default 2
    #[case("[Fe+2/1^4*5v6]", AtomSpec::new(e!(Fe), 2, 1, 0, 0, 4, 5, 0, 6))] // Ion, no accepted pairs
    #[case("[Fe+2/1<4^4*5v6]", AtomSpec::new(e!(Fe), 2, 1, 0, 4, 4, 5, 0, 6))] // Ion, 4 accepted lone pairs
    #[case("[N*2/1^3]", AtomSpec::new(e!(N), 0, 1, 0, 0, 3, 2, 0, 0))] // Arbitrary order, explicit mult
    #[case("[N/1>1]", AtomSpec::new(e!(N), 0, 1, 1, 0, 0, 1, 0,0))] // Donated lone pair
    #[case("[P^v5]", AtomSpec::new(e!(P), 0, 0, 0, 0, 1, 2, 0, 5))] // Unpaired default 1, Valence 5
    #[case("[S/v3+]", AtomSpec::new(e!(S), 1, 1, 0, 0, 0, 1, 0, 3))] // Arbitrary order, lp default 1
    #[case("[Gd+3^7]", AtomSpec::new(e!(Gd), 3, 0, 0, 0, 7, 8, 0, 0))] // High unpaired, default mult
    #[case("[C^2v2]", AtomSpec::new(e!(C), 0, 0, 0, 0, 2, 3, 0, 2))] // Triplet carbene
    #[case("[C^2*1v2]", AtomSpec::new(e!(C), 0, 0, 0, 0, 2, 1, 0, 2))] // Singlet carbene
    #[case("[C+/0^0v3]", AtomSpec::new(e!(C), 1, 0, 0, 0, 0, 1, 0, 3))] // Explicit charge, lp, unpaired
    fn test_atom_spec_from_str(#[case] s: &str, #[case] expected: AtomSpec) {
        assert_eq!(AtomSpec::from_str(s).unwrap(), expected);
    }

    #[rstest]
    #[case("")] // Empty
    #[case("C")] // No brackets
    #[case("[C")] // Mismatched bracket
    #[case("[]")] // No element
    #[case("[Xx]")] // Invalid element
    #[case("[C++]")] // Duplicate charge definition
    #[case("[C//]")] // Duplicate lone pair definition
    #[case("[C>>]")] // Duplicate donated pair definition
    #[case("[C<1<1]")] // Duplicate accepted pair definition
    #[case("[C^^1]")] // Duplicate unpaired electron definition
    #[case("[C**]")] // Duplicate multiplicity definition
    #[case("[CHH]")] // Duplicate implicit hydrogen definition
    #[case("[Cv1v1]")] // Duplicate valence definition
    #[case("[C+?]")] // Invalid char in properties
    #[case("[C+1a]")] // Invalid char after property
    #[case("[C+1.0]")] // Non-integer
    #[case("[C^2*2]")] // Invalid multiplicity
    #[case("[CH4v3]")] // Invalid implicit hydrogen definition
    fn test_atom_spec_from_str_invalid(#[case] s: &str) {
        assert!(AtomSpec::from_str(s).is_err());
    }

    #[test]
    fn test_atom_spec_serialize() {
        let vs = AtomSpec::new(e!(C), 0, 1, 0, 0, 2, 3, 0, 4);
        let serialized = serde_json::to_string(&vs).unwrap();
        assert_eq!(
            serialized,
            r#"{"element":"C","charge":0,"lone_pairs":1,"donated_pairs":0,"accepted_pairs":0,"unpaired_electrons":2,"multiplicity":3,"implicit_hydrogens":0,"valence":4}"#
        );
    }

    #[test]
    fn test_atom_spec_deserialize() {
        let serialized = r#"{"element":"C","charge":0,"lone_pairs":1,"donated_pairs":0,"accepted_pairs":0,"unpaired_electrons":2,"multiplicity":3,"implicit_hydrogens":0,"valence":4}"#;
        let vs = serde_json::from_str::<AtomSpec>(serialized).unwrap();
        assert_eq!(vs, AtomSpec::new(e!(C), 0, 1, 0, 0, 2, 3, 0, 4));
    }

    #[test]
    fn test_atom_spec_macro() {
        assert_eq!(
            a!("[C+0/1^2*3H0v4]"),
            AtomSpec::new(e!(C), 0, 1, 0, 0, 2, 3, 0, 4)
        );
    }
}
