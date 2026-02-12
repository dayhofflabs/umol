//! Atom spec definitions for GraphIR.

use std::fmt;
use std::str::FromStr;

use regex::Regex;
use umol_data::Element;

use super::error::ResolutionError;

/// AtomSpec is a generalization of the atomic valence state for modeling
/// covalent and dative bonding by molecular graphs.
/// Does not uniquely define an atomic term symbol for quantum chemical calculations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AtomSpec {
    element: Element,
    charge: i32,
    lone_pairs: u32,
    donated_pairs: u32,
    accepted_pairs: u32,
    unpaired_electrons: u32,
    multiplicity: u32,
    implicit_hydrogens: u32,
    valence: u32,
}

impl AtomSpec {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        element: Element,
        charge: i32,
        lone_pairs: u32,
        donated_pairs: u32,
        accepted_pairs: u32,
        unpaired_electrons: u32,
        multiplicity: u32,
        implicit_hydrogens: u32,
        valence: u32,
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
            (element.valence_e() as i32 - charge) >= 0,
            "Charge must be less than or equal to valence electrons"
        );
        debug_assert!(
            implicit_hydrogens <= element.max_implicit_hydrogens() as u32,
            "Implicit hydrogens must be less than or equal to max implicit hydrogens"
        );
        debug_assert!(
            implicit_hydrogens <= valence,
            "Implicit hydrogens must be less than or equal to valence"
        );
        debug_assert!(
            valence <= element.max_valence() as u32,
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

    pub fn charge(&self) -> i32 {
        self.charge
    }

    pub fn lone_pairs(&self) -> u32 {
        self.lone_pairs
    }

    pub fn donated_pairs(&self) -> u32 {
        self.donated_pairs
    }

    pub fn accepted_pairs(&self) -> u32 {
        self.accepted_pairs
    }

    pub fn unpaired_electrons(&self) -> u32 {
        self.unpaired_electrons
    }

    pub fn electron_count(&self) -> u32 {
        2 * self.lone_pairs + self.unpaired_electrons
    }

    pub fn multiplicity(&self) -> u32 {
        self.multiplicity
    }

    pub fn implicit_hydrogens(&self) -> u32 {
        self.implicit_hydrogens
    }

    pub fn valence(&self) -> u32 {
        self.valence
    }
}

impl fmt::Display for AtomSpec {
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
    type Error = ResolutionError;

    fn try_from(s: &str) -> Result<Self, ResolutionError> {
        if s.is_empty() {
            return Err(ResolutionError::InvalidAtomSpec("Empty atom spec".into()));
        }

        let val_state_pattern = Regex::new(r"^\[([A-Z][a-z]?)((?:[-+/<>*Hv^]\d*)*)]$").unwrap();
        if !val_state_pattern.is_match(s) {
            return Err(ResolutionError::InvalidAtomSpec(s.to_string()));
        }

        let caps = val_state_pattern.captures(s).unwrap();
        let element = caps[1]
            .parse::<Element>()
            .map_err(|e| ResolutionError::InvalidAtomSpec(format!("Invalid element: {}", e)))?;
        let properties = caps.get(2);
        if properties.is_none() {
            return Ok(AtomSpec::new(element, 0, 0, 0, 0, 0, 1, 0, 0));
        }
        let properties = properties.unwrap().as_str();

        let property_pattern = Regex::new(r"([-+/<>*Hv^])(\d*)").unwrap();
        let mut pos_charge: Option<i32> = None;
        let mut neg_charge: Option<i32> = None;
        let mut lone_pairs: Option<u32> = None;
        let mut donated: Option<u32> = None;
        let mut accepted: Option<u32> = None;
        let mut unpaired: Option<u32> = None;
        let mut multiplicity: Option<u32> = None;
        let mut implicit_hydrogens: Option<u32> = None;
        let mut valence: Option<u32> = None;

        for cap in property_pattern.captures_iter(properties) {
            let property = &cap[1];
            let value = &cap[2];
            match property {
                "+" => {
                    if pos_charge.is_some() || neg_charge.is_some() {
                        return Err(ResolutionError::InvalidAtomSpec(
                            "Duplicate charge definition".to_string(),
                        )
                        .into());
                    } else {
                        pos_charge = if value.is_empty() {
                            Some(1)
                        } else {
                            value.parse::<i32>().ok()
                        };
                    }
                }
                "-" => {
                    if neg_charge.is_some() || pos_charge.is_some() {
                        return Err(ResolutionError::InvalidAtomSpec(
                            "Duplicate charge definition".to_string(),
                        )
                        .into());
                    } else {
                        neg_charge = if value.is_empty() {
                            Some(-1)
                        } else {
                            value.parse::<i32>().ok().map(|v| -v)
                        };
                    }
                }
                "/" => {
                    if lone_pairs.is_some() {
                        return Err(ResolutionError::InvalidAtomSpec(
                            "Duplicate lone pair definition".to_string(),
                        )
                        .into());
                    } else {
                        lone_pairs = if value.is_empty() {
                            Some(1)
                        } else {
                            value.parse::<u32>().ok()
                        };
                    }
                }
                ">" => {
                    if donated.is_some() {
                        return Err(ResolutionError::InvalidAtomSpec(
                            "Duplicate donated pair definition".to_string(),
                        )
                        .into());
                    } else {
                        donated = if value.is_empty() {
                            Some(1)
                        } else {
                            value.parse::<u32>().ok()
                        };
                    }
                }
                "<" => {
                    if accepted.is_some() {
                        return Err(ResolutionError::InvalidAtomSpec(
                            "Duplicate accepted pair definition".to_string(),
                        )
                        .into());
                    } else {
                        accepted = if value.is_empty() {
                            Some(1)
                        } else {
                            value.parse::<u32>().ok()
                        };
                    }
                }
                "^" => {
                    if unpaired.is_some() {
                        return Err(ResolutionError::InvalidAtomSpec(
                            "Duplicate unpaired electron specification".to_string(),
                        )
                        .into());
                    } else {
                        unpaired = if value.is_empty() {
                            Some(1)
                        } else {
                            value.parse::<u32>().ok()
                        };
                    }
                }
                "*" => {
                    multiplicity = if multiplicity.is_some() {
                        return Err(ResolutionError::InvalidAtomSpec(
                            "Duplicate multiplicity".to_string(),
                        )
                        .into());
                    } else if value.is_empty() {
                        Some(1)
                    } else {
                        value.parse::<u32>().ok()
                    }
                }
                "H" => {
                    if implicit_hydrogens.is_some() {
                        return Err(ResolutionError::InvalidAtomSpec(
                            "Duplicate implicit hydrogen specification".to_string(),
                        )
                        .into());
                    } else {
                        implicit_hydrogens = if value.is_empty() {
                            Some(1)
                        } else {
                            value.parse::<u32>().ok()
                        };
                    }
                }
                "v" => {
                    valence = if valence.is_some() {
                        return Err(ResolutionError::InvalidAtomSpec(
                            "Duplicate valence".to_string(),
                        )
                        .into());
                    } else if value.is_empty() {
                        return Err(ResolutionError::InvalidAtomSpec(
                            "Duplicate valence".to_string(),
                        )
                        .into());
                    } else {
                        value.parse::<u32>().ok()
                    }
                }
                _ => unreachable!(),
            }
        }

        let charge = pos_charge.unwrap_or(neg_charge.unwrap_or(0));
        let lone_pairs = lone_pairs.unwrap_or(0);
        let donated = donated.unwrap_or(0);
        let accepted = accepted.unwrap_or(0);
        let unpaired = unpaired.unwrap_or(0);
        let multiplicity = multiplicity.unwrap_or(unpaired + 1);
        let implicit_hydrogens = implicit_hydrogens.unwrap_or(0);
        let valence = valence.unwrap_or(0);

        if multiplicity == 0
            || multiplicity > unpaired + 1
            || (unpaired + 1 - multiplicity) % 2 != 0
        {
            return Err(ResolutionError::InvalidAtomSpec(format!(
                "Invalid multiplicity: {}",
                multiplicity
            )));
        }
        if donated > lone_pairs {
            return Err(ResolutionError::InvalidAtomSpec(format!(
                "Invalid donated pairs: {}",
                donated
            )));
        }
        if (element.valence_e() as i32 - charge) < 0 {
            return Err(ResolutionError::InvalidAtomSpec(format!(
                "Invalid charge: {}",
                charge
            )));
        }
        if implicit_hydrogens > element.max_implicit_hydrogens() as u32 {
            return Err(ResolutionError::InvalidAtomSpec(format!(
                "Invalid implicit hydrogens: {}",
                implicit_hydrogens
            )));
        }
        if implicit_hydrogens > valence {
            return Err(ResolutionError::InvalidAtomSpec(format!(
                "Invalid valence (H): {}",
                valence
            )));
        }
        if valence > element.max_valence() as u32 {
            return Err(ResolutionError::InvalidAtomSpec(format!(
                "Invalid valence: {}",
                valence
            )));
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

// impl FromStr for AtomSpec {
//     type Error = ResolutionError;

//     fn from_str(s: &str) -> Result<Self, Self::Error> {
//         let atom_spec = Self::try_from(s).map_err(|e| e.into())?;
//         Ok(atom_spec)
//     }
// }
