//! Electronic configurations of atoms and atomic ions
use std::cmp;
use std::collections::HashMap;
use std::fmt::{self, Display};
use std::ops::Deref;
use std::sync::LazyLock;

use map_macro::hash_map;

use crate::{e, occ, Element, Occupation, MAX_UNPAIRED_ELECTRONS};

/// Electronic configuration of atom or atomic ion
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Configuration {
    element: Element,
    core_element: Option<Element>,
    valence_occupation: Occupation,
}

impl Configuration {
    /// Create new configuration
    pub fn new(
        element: Element,
        core_element: Option<Element>,
        valence_occupation: Occupation,
    ) -> Self {
        Self {
            element,
            core_element,
            valence_occupation,
        }
    }

    /// Return element
    pub fn element(&self) -> Element {
        self.element
    }

    /// Return core element
    pub fn core_element(&self) -> Option<Element> {
        self.core_element
    }

    /// Return core occupation
    pub fn core_occupation(&self) -> Option<Occupation> {
        self.core_element.map(get_total_occupation)
    }

    /// Return computed valence occupation
    pub fn valence_occupation(&self) -> Occupation {
        self.valence_occupation
    }

    /// Return total occupation
    pub fn total_occupation(&self) -> Occupation {
        get_total_occupation(self.element)
    }

    /// Return core electron count
    pub fn core_electron_count(&self) -> u8 {
        self.core_element
            .map(|core| core.atomic_number())
            .unwrap_or(0)
    }

    /// Return valence electron count
    pub fn valence_electron_count(&self) -> u8 {
        let valence_occupation = self.valence_occupation();
        valence_occupation.s()
            + valence_occupation.p()
            + valence_occupation.d()
            + valence_occupation.f()
    }

    /// Return computed electron count
    pub fn electron_count(&self) -> u8 {
        self.core_electron_count() + self.valence_electron_count()
    }

    /// Return atomic number
    pub fn atomic_number(&self) -> u8 {
        self.element.atomic_number()
    }

    /// Return charge
    pub fn charge(&self) -> i8 {
        self.atomic_number() as i8 - self.electron_count() as i8
    }
}

impl Display for Configuration {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}{}",
            self.core_element()
                .map(|core| format!("[{}] ", core.symbol()))
                .unwrap_or("".to_string()),
            self.valence_occupation(),
        )
    }
}

// TODO: Implement TryFrom<&str>, FromStr for Configuration
// TODO: Implement Serialize, Deserialize for Configuration

/// Ground state configuration
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct GroundState(Configuration);

impl GroundState {
    pub fn new(element: Element) -> Self {
        if let Some(gs) = GROUND_STATE_EXCEPTIONS.get(&element) {
            *gs
        } else {
            GroundState(get_aufbau_configuration(element))
        }
    }
}

impl Display for GroundState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Deref for GroundState {
    type Target = Configuration;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl From<GroundState> for Configuration {
    fn from(gs: GroundState) -> Self {
        gs.0
    }
}

pub const MAX_N_QUANTUM_NUMBER: u8 = 7;
/// Maximum l quantum number
pub const MAX_L_QUANTUM_NUMBER: u8 = 3;

// Order of filling of atomic orbitals according to Madelung rule
// (n, l, orbital_count, closing subshell)
#[rustfmt::skip]
static MADELUNG_ORDER: [(u8, u8, u8, bool); 19] = [
    (1, 0, 2, true), (2, 0, 2, false), (2, 1, 6, true), (3, 0, 2, false), (3, 1, 6, true),
    (4, 0, 2, false), (3, 2, 10, false), (4, 1, 6, true), (5, 0, 2, false), (4, 2, 10, false),
    (5, 1, 6, true), (6, 0, 2, false), (4, 3, 14, false), (5, 2, 10, false), (6, 1, 6, true),
    (7, 0, 2, false), (5, 3, 14, false), (6, 2, 10, false), (7, 1, 6, true)];

/// Atomic ground state exceptions
static GROUND_STATE_EXCEPTIONS: LazyLock<HashMap<Element, GroundState>> = LazyLock::new(|| {
    hash_map! {
        Element::Cr => GroundState(Configuration::new(e!(Cr), Some(e!(Ar)), occ!(s1d5))),
        Element::Cu => GroundState(Configuration::new(e!(Cu), Some(e!(Ar)), occ!(s1d10))),
        Element::Nb => GroundState(Configuration::new(e!(Nb), Some(e!(Kr)), occ!(s1d4))),
        Element::Mo => GroundState(Configuration::new(e!(Mo), Some(e!(Kr)), occ!(s1d5))),
        Element::Ru => GroundState(Configuration::new(e!(Ru), Some(e!(Kr)), occ!(s1d7))),
        Element::Rh => GroundState(Configuration::new(e!(Rh), Some(e!(Kr)), occ!(s1d8))),
        Element::Pd => GroundState(Configuration::new(e!(Pd), Some(e!(Kr)), occ!(s1d10))),
        Element::Ag => GroundState(Configuration::new(e!(Ag), Some(e!(Kr)), occ!(s1d10))),
        Element::La => GroundState(Configuration::new(e!(La), Some(e!(Xe)), occ!(s2d1))),
        Element::Ce => GroundState(Configuration::new(e!(Ce), Some(e!(Xe)), occ!(s2d1f1))),
        Element::Gd => GroundState(Configuration::new(e!(Gd), Some(e!(Xe)), occ!(s2d1f7))),
        Element::Pt => GroundState(Configuration::new(e!(Pt), Some(e!(Xe)), occ!(s1d9f14))),
        Element::Au => GroundState(Configuration::new(e!(Au), Some(e!(Xe)), occ!(s1d10f14))),
        Element::Ac => GroundState(Configuration::new(e!(Ac), Some(e!(Rn)), occ!(s2d1))),
        Element::Th => GroundState(Configuration::new(e!(Th), Some(e!(Rn)), occ!(s2d2))),
        Element::Pa => GroundState(Configuration::new(e!(Pa), Some(e!(Rn)), occ!(s2d1f2))),
        Element::U => GroundState(Configuration::new(e!(U), Some(e!(Rn)), occ!(s2d1f3))),
        Element::Np => GroundState(Configuration::new(e!(Np), Some(e!(Rn)), occ!(s2d1f4))),
        Element::Cm => GroundState(Configuration::new(e!(Cm), Some(e!(Rn)), occ!(s2d1f7))),
    }
});

/// Get total occupation
fn get_total_occupation(element: Element) -> Occupation {
    let mut remaining = element.atomic_number();
    let mut occupations = [0; 4];
    for (_, l, capacity, _) in MADELUNG_ORDER {
        let occupation = cmp::min(remaining, capacity);
        occupations[l as usize] += occupation;
        remaining -= occupation;
    }
    Occupation::new(
        occupations[0],
        occupations[1],
        occupations[2],
        occupations[3],
    )
}
/// Compute configurations from Aufbau principle, Madelung rule, and Hund's rules
fn get_aufbau_configuration(element: Element) -> Configuration {
    if GROUND_STATE_EXCEPTIONS.contains_key(&element) {
        // TODO: Add warning when logging is implemented
        eprintln!("Element {} has exceptional configuration", element);
    }
    let mut remaining = element.atomic_number();
    let mut subshell_occupations = Vec::new();
    let mut closing_subshell = (0, 0, 0, false);

    for subshell @ (n, l, capacity, is_closing) in MADELUNG_ORDER {
        let occupation = cmp::min(remaining, capacity);
        subshell_occupations.push((n, l, capacity, occupation));
        remaining -= occupation;
        if is_closing && remaining != 0 {
            closing_subshell = subshell;
        }
        if remaining == 0 {
            break;
        }
    }

    // Add unpaired electrons from last partially filled subshell
    let mut unpaired = 0;
    if let Some(&(_, _, capacity, occupation)) = subshell_occupations.last() {
        if occupation < capacity {
            if occupation <= capacity / 2 {
                // Less than half filled shell
                unpaired = occupation;
            } else {
                // More than half filled shell
                unpaired = capacity - occupation;
            }
        }
    }
    debug_assert!(
        subshell_occupations
            .iter()
            .map(|counts| counts.3)
            .sum::<u8>()
            == element.atomic_number(),
        "Total occupation must be equal to atomic number"
    );
    debug_assert!(
        unpaired <= MAX_UNPAIRED_ELECTRONS,
        "Unpaired electrons must be less than or equal to {}",
        MAX_UNPAIRED_ELECTRONS
    );
    debug_assert!(
        closing_subshell.0 <= MAX_N_QUANTUM_NUMBER,
        "Closing subshell n must be less than or equal to {}",
        MAX_N_QUANTUM_NUMBER
    );
    debug_assert!(
        closing_subshell.1 <= MAX_L_QUANTUM_NUMBER,
        "Closing subshell l must be less than or equal to {}",
        MAX_L_QUANTUM_NUMBER
    );

    let valence_occupation = get_valence_occupation(subshell_occupations, closing_subshell);
    let core_element = get_core_element(element);
    Configuration::new(element, core_element, valence_occupation)
}

/// Compute core element
fn get_core_element(element: Element) -> Option<Element> {
    match element.atomic_number() {
        1..=2 => None,
        3..=10 => Some(Element::He),
        11..=18 => Some(Element::Ne),
        19..=36 => Some(Element::Ar),
        37..=54 => Some(Element::Kr),
        55..=86 => Some(Element::Xe),
        87..=118 => Some(Element::Rn),
        _ => unreachable!(),
    }
}

/// Compute valence occupation
fn get_valence_occupation(
    subshell_counts: Vec<(u8, u8, u8, u8)>,
    closing_subshell: (u8, u8, u8, bool),
) -> Occupation {
    let mut valence_counts: [u8; 4] = [0; 4];

    let mut in_valence = closing_subshell == (0, 0, 0, false);
    for (n, l, _, count) in subshell_counts {
        if count == 0 {
            break;
        }

        if in_valence {
            valence_counts[l as usize] += count;
        }

        if !in_valence && (n, l) == (closing_subshell.0, closing_subshell.1) {
            in_valence = true;
        }
    }

    Occupation::new(
        valence_counts[0],
        valence_counts[1],
        valence_counts[2],
        valence_counts[3],
    )
}

#[cfg(test)]
mod tests {
    use rstest::*;

    use super::*;
    use crate::{e, occ};

    #[rstest]
    #[case(e!(H), Configuration::new(e!(H), None, occ!(s1)))]
    #[case(e!(He), Configuration::new(e!(He), None, occ!(s2)))]
    #[case(e!(Li), Configuration::new(e!(Li), Some(e!(He)), occ!(s1)))]
    #[case(e!(Be), Configuration::new(e!(Be), Some(e!(He)), occ!(s2)))]
    #[case(e!(B), Configuration::new(e!(B), Some(e!(He)), occ!(s2p1)))]
    #[case(e!(C), Configuration::new(e!(C), Some(e!(He)), occ!(s2p2)))]
    #[case(e!(N), Configuration::new(e!(N), Some(e!(He)), occ!(s2p3)))]
    #[case(e!(O), Configuration::new(e!(O), Some(e!(He)), occ!(s2p4)))]
    #[case(e!(F), Configuration::new(e!(F), Some(e!(He)), occ!(s2p5)))]
    #[case(e!(Ne), Configuration::new(e!(Ne), Some(e!(He)), occ!(s2p6)))]
    #[case(e!(Na), Configuration::new(e!(Na), Some(e!(Ne)), occ!(s1)))]
    #[case(e!(Mg), Configuration::new(e!(Mg), Some(e!(Ne)), occ!(s2)))]
    #[case(e!(Al), Configuration::new(e!(Al), Some(e!(Ne)), occ!(s2p1)))]
    #[case(e!(Si), Configuration::new(e!(Si), Some(e!(Ne)), occ!(s2p2)))]
    #[case(e!(P), Configuration::new(e!(P), Some(e!(Ne)), occ!(s2p3)))]
    #[case(e!(S), Configuration::new(e!(S), Some(e!(Ne)), occ!(s2p4)))]
    #[case(e!(Cl), Configuration::new(e!(Cl), Some(e!(Ne)), occ!(s2p5)))]
    #[case(e!(Ar), Configuration::new(e!(Ar), Some(e!(Ne)), occ!(s2p6)))]
    #[case(e!(K), Configuration::new(e!(K), Some(e!(Ar)), occ!(s1)))]
    #[case(e!(Ca), Configuration::new(e!(Ca), Some(e!(Ar)), occ!(s2)))]
    #[case(e!(Sc), Configuration::new(e!(Sc), Some(e!(Ar)), occ!(s2d1)))]
    #[case(e!(Ti), Configuration::new(e!(Ti), Some(e!(Ar)), occ!(s2d2)))]
    #[case(e!(V), Configuration::new(e!(V), Some(e!(Ar)), occ!(s2d3)))]
    #[case(e!(Cr), Configuration::new(e!(Cr), Some(e!(Ar)), occ!(s2d4)))] // Exception
    #[case(e!(Mn), Configuration::new(e!(Mn), Some(e!(Ar)), occ!(s2d5)))]
    #[case(e!(Fe), Configuration::new(e!(Fe), Some(e!(Ar)), occ!(s2d6)))]
    #[case(e!(Co), Configuration::new(e!(Co), Some(e!(Ar)), occ!(s2d7)))]
    #[case(e!(Ni), Configuration::new(e!(Ni), Some(e!(Ar)), occ!(s2d8)))]
    #[case(e!(Cu), Configuration::new(e!(Cu), Some(e!(Ar)), occ!(s2d9)))] // Exception
    #[case(e!(Zn), Configuration::new(e!(Zn), Some(e!(Ar)), occ!(s2d10)))]
    #[case(e!(Ce), Configuration::new(e!(Ce), Some(e!(Xe)), occ!(s2f2)))] // Exception
    #[case(e!(Eu), Configuration::new(e!(Eu), Some(e!(Xe)), occ!(s2f7)))]
    #[case(e!(Gd), Configuration::new(e!(Gd), Some(e!(Xe)), occ!(s2f8)))] // Exception
    #[case(e!(Xe), Configuration::new(e!(Xe), Some(e!(Kr)), occ!(s2p6d10)))]
    #[case(e!(Pb), Configuration::new(e!(Pb), Some(e!(Xe)), occ!(s2p2d10f14)))]
    fn test_get_aufbau_configuration(#[case] element: Element, #[case] expected: Configuration) {
        assert_eq!(get_aufbau_configuration(element), expected);
    }

    #[rstest]
    #[case(Configuration::new(e!(H), None, occ!(s1)), None, occ!(s1))]
    #[case(Configuration::new(e!(Li), Some(e!(He)), occ!(s1)), Some(occ!(s2)), occ!(s1))]
    #[case(Configuration::new(e!(Be), Some(e!(He)), occ!(s2)), Some(occ!(s2)), occ!(s2))]
    #[case(Configuration::new(e!(C), Some(e!(He)), occ!(s2p2)), Some(occ!(s2)), occ!(s2p2))]
    #[case(Configuration::new(e!(Ne), Some(e!(He)), occ!(s2p6)), Some(occ!(s2)), occ!(s2p6))]
    #[case(Configuration::new(e!(Cr), Some(e!(Ar)), occ!(s1d5)), Some(occ!(s6p12)), occ!(s1d5))]
    #[case(Configuration::new(e!(Xe), Some(e!(Kr)), occ!(s2p6d10)), Some(occ!(s8p18d10)), occ!(s2p6d10))]
    #[case(Configuration::new(e!(Ce), Some(e!(Xe)), occ!(s2d1f1)), Some(occ!(s10p24d20)), occ!(s2d1f1))]
    fn test_configuration_properties(
        #[case] config: Configuration,
        #[case] expected_core_occupation: Option<Occupation>,
        #[case] expected_valence_occupation: Occupation,
    ) {
        assert_eq!(config.core_occupation(), expected_core_occupation);
        assert_eq!(config.valence_occupation(), expected_valence_occupation);
    }

    #[rstest]
    #[case(Configuration::new(e!(H), None, occ!(s1)), "s1")]
    #[case(Configuration::new(e!(He), None, occ!(s2)), "s2")]
    #[case(Configuration::new(e!(Li), Some(e!(He)), occ!(s1)), "[He] s1")]
    #[case(Configuration::new(e!(Be), Some(e!(He)), occ!(s2)), "[He] s2")]
    #[case(Configuration::new(e!(B), Some(e!(He)), occ!(s2p1)), "[He] s2p1")]
    fn test_configuration_display(#[case] config: Configuration, #[case] expected: &str) {
        assert_eq!(format!("{}", config), expected);
    }
    #[rstest]
    #[case(e!(H), GroundState(Configuration::new(e!(H), None, occ!(s1))))]
    #[case(e!(C), GroundState(Configuration::new(e!(C), Some(e!(He)), occ!(s2p2))))]
    #[case(e!(Cr), GroundState(Configuration::new(e!(Cr), Some(e!(Ar)), occ!(s1d5))))] // Exception
    #[case(e!(Mn), GroundState(Configuration::new(e!(Mn), Some(e!(Ar)), occ!(s2d5))))]
    #[case(e!(Cu), GroundState(Configuration::new(e!(Cu), Some(e!(Ar)), occ!(s1d10))))] // Exception
    #[case(e!(La), GroundState(Configuration::new(e!(La), Some(e!(Xe)), occ!(s2d1))))] // Exception
    #[case(e!(Eu), GroundState(Configuration::new(e!(Eu), Some(e!(Xe)), occ!(s2f7))))]
    #[case(e!(Gd), GroundState(Configuration::new(e!(Gd), Some(e!(Xe)), occ!(s2d1f7))))] // Exception
    #[case(e!(Xe), GroundState(Configuration::new(e!(Xe), Some(e!(Kr)), occ!(s2p6d10))))]
    #[case(e!(Pb), GroundState(Configuration::new(e!(Pb), Some(e!(Xe)), occ!(s2p2d10f14))))]
    fn test_ground_state_new(#[case] element: Element, #[case] expected: GroundState) {
        assert_eq!(GroundState::new(element), expected);
    }

    #[rstest]
    #[case(GroundState::new(Element::H), "s1")]
    #[case(GroundState::new(Element::He), "s2")]
    #[case(GroundState::new(Element::Li), "[He] s1")]
    #[case(GroundState::new(Element::Be), "[He] s2")]
    #[case(GroundState::new(Element::B), "[He] s2p1")]
    fn test_ground_state_display(#[case] gs: GroundState, #[case] expected: &str) {
        assert_eq!(format!("{}", gs), expected);
    }
}
