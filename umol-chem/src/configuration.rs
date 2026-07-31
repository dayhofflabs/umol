//! Electronic configurations of atoms and atomic ions
use std::cmp;
use std::collections::HashMap;
use std::fmt::{self, Display};
use std::ops::Deref;
use std::sync::LazyLock;

use map_macro::hash_map;

use crate::element::Element;
use crate::occupation::Occupation;
use crate::{e, occ};

/// Electronic configuration of atom or atomic ion
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Configuration {
    element: Element,
    core_element: Option<Element>,
    valence_occupation: Occupation,
}

// TODO: Restructure: Add Subshell struct instead of (n, l, capacity, is_closing) tuple
// Use subsell struct in configuration, make occupation derivable instead of stored
// Allow addition for occupations
// Add more structured way of determining the closing subshell for element (probably constructor on
// Subshell).
// Remove temporary allocations
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

    /// Compute configurations from Aufbau principle, Madelung rule, and Hund's rules
    pub fn aufbau(element: Element) -> Self {
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

        debug_assert!(
            subshell_occupations
                .iter()
                .map(|counts| counts.3)
                .sum::<u8>()
                == element.atomic_number(),
            "Total occupation must be equal to atomic number"
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

        Configuration::new(element, element.core(), valence_occupation)
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
            Self(Configuration::aufbau(element))
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
// TODO: Implement Subshell struct instead of tuple
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

/// Number of unpaired electrons when `electrons` occupy `orbitals` degenerate
/// spatial orbitals at maximum multiplicity (Hund's first rule): each orbital
/// is filled singly before any is doubly occupied. `electrons` must not exceed
/// `2 * orbitals`.
pub fn hund_rule_unpaired(electrons: u8, capacity: u8) -> u8 {
    debug_assert!(electrons <= capacity);
    if electrons <= capacity / 2 {
        electrons
    } else {
        capacity - electrons
    }
}

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
    #[case::h(e!(H), Configuration::new(e!(H), None, occ!(s1)))]
    #[case::he(e!(He), Configuration::new(e!(He), None, occ!(s2)))]
    #[case::li(e!(Li), Configuration::new(e!(Li), Some(e!(He)), occ!(s1)))]
    #[case::be(e!(Be), Configuration::new(e!(Be), Some(e!(He)), occ!(s2)))]
    #[case::b(e!(B), Configuration::new(e!(B), Some(e!(He)), occ!(s2p1)))]
    #[case::c(e!(C), Configuration::new(e!(C), Some(e!(He)), occ!(s2p2)))]
    #[case::n(e!(N), Configuration::new(e!(N), Some(e!(He)), occ!(s2p3)))]
    #[case::o(e!(O), Configuration::new(e!(O), Some(e!(He)), occ!(s2p4)))]
    #[case::f(e!(F), Configuration::new(e!(F), Some(e!(He)), occ!(s2p5)))]
    #[case::ne(e!(Ne), Configuration::new(e!(Ne), Some(e!(He)), occ!(s2p6)))]
    #[case::na(e!(Na), Configuration::new(e!(Na), Some(e!(Ne)), occ!(s1)))]
    #[case::mg(e!(Mg), Configuration::new(e!(Mg), Some(e!(Ne)), occ!(s2)))]
    #[case::al(e!(Al), Configuration::new(e!(Al), Some(e!(Ne)), occ!(s2p1)))]
    #[case::si(e!(Si), Configuration::new(e!(Si), Some(e!(Ne)), occ!(s2p2)))]
    #[case::p(e!(P), Configuration::new(e!(P), Some(e!(Ne)), occ!(s2p3)))]
    #[case::s(e!(S), Configuration::new(e!(S), Some(e!(Ne)), occ!(s2p4)))]
    #[case::cl(e!(Cl), Configuration::new(e!(Cl), Some(e!(Ne)), occ!(s2p5)))]
    #[case::ar(e!(Ar), Configuration::new(e!(Ar), Some(e!(Ne)), occ!(s2p6)))]
    #[case::k(e!(K), Configuration::new(e!(K), Some(e!(Ar)), occ!(s1)))]
    #[case::ca(e!(Ca), Configuration::new(e!(Ca), Some(e!(Ar)), occ!(s2)))]
    #[case::sc(e!(Sc), Configuration::new(e!(Sc), Some(e!(Ar)), occ!(s2d1)))]
    #[case::ti(e!(Ti), Configuration::new(e!(Ti), Some(e!(Ar)), occ!(s2d2)))]
    #[case::v(e!(V), Configuration::new(e!(V), Some(e!(Ar)), occ!(s2d3)))]
    #[case::cr(e!(Cr), Configuration::new(e!(Cr), Some(e!(Ar)), occ!(s2d4)))] // Exception
    #[case::mn(e!(Mn), Configuration::new(e!(Mn), Some(e!(Ar)), occ!(s2d5)))]
    #[case::fe(e!(Fe), Configuration::new(e!(Fe), Some(e!(Ar)), occ!(s2d6)))]
    #[case::co(e!(Co), Configuration::new(e!(Co), Some(e!(Ar)), occ!(s2d7)))]
    #[case::ni(e!(Ni), Configuration::new(e!(Ni), Some(e!(Ar)), occ!(s2d8)))]
    #[case::cu(e!(Cu), Configuration::new(e!(Cu), Some(e!(Ar)), occ!(s2d9)))] // Exception
    #[case::zn(e!(Zn), Configuration::new(e!(Zn), Some(e!(Ar)), occ!(s2d10)))]
    #[case::ce(e!(Ce), Configuration::new(e!(Ce), Some(e!(Xe)), occ!(s2f2)))] // Exception
    #[case::eu(e!(Eu), Configuration::new(e!(Eu), Some(e!(Xe)), occ!(s2f7)))]
    #[case::gd(e!(Gd), Configuration::new(e!(Gd), Some(e!(Xe)), occ!(s2f8)))] // Exception
    #[case::xe(e!(Xe), Configuration::new(e!(Xe), Some(e!(Kr)), occ!(s2p6d10)))]
    #[case::pb(e!(Pb), Configuration::new(e!(Pb), Some(e!(Xe)), occ!(s2p2d10f14)))]
    fn test_aufbau_configuration(#[case] element: Element, #[case] expected: Configuration) {
        assert_eq!(Configuration::aufbau(element), expected);
    }

    #[rstest]
    #[case::h(Configuration::new(e!(H), None, occ!(s1)), None, occ!(s1))]
    #[case::li(Configuration::new(e!(Li), Some(e!(He)), occ!(s1)), Some(occ!(s2)), occ!(s1))]
    #[case::be(Configuration::new(e!(Be), Some(e!(He)), occ!(s2)), Some(occ!(s2)), occ!(s2))]
    #[case::c(Configuration::new(e!(C), Some(e!(He)), occ!(s2p2)), Some(occ!(s2)), occ!(s2p2))]
    #[case::ne(Configuration::new(e!(Ne), Some(e!(He)), occ!(s2p6)), Some(occ!(s2)), occ!(s2p6))]
    #[case::cr(Configuration::new(e!(Cr), Some(e!(Ar)), occ!(s1d5)), Some(occ!(s6p12)), occ!(s1d5))]
    #[case::xe(Configuration::new(e!(Xe), Some(e!(Kr)), occ!(s2p6d10)), Some(occ!(s8p18d10)), occ!(s2p6d10))]
    #[case::ce(Configuration::new(e!(Ce), Some(e!(Xe)), occ!(s2d1f1)), Some(occ!(s10p24d20)), occ!(s2d1f1))]
    fn test_configuration_properties(
        #[case] config: Configuration,
        #[case] expected_core_occupation: Option<Occupation>,
        #[case] expected_valence_occupation: Occupation,
    ) {
        assert_eq!(config.core_occupation(), expected_core_occupation);
        assert_eq!(config.valence_occupation(), expected_valence_occupation);
    }

    #[rstest]
    #[case::h(Configuration::new(e!(H), None, occ!(s1)), "s1")]
    #[case::he(Configuration::new(e!(He), None, occ!(s2)), "s2")]
    #[case::li(Configuration::new(e!(Li), Some(e!(He)), occ!(s1)), "[He] s1")]
    #[case::be(Configuration::new(e!(Be), Some(e!(He)), occ!(s2)), "[He] s2")]
    #[case::b(Configuration::new(e!(B), Some(e!(He)), occ!(s2p1)), "[He] s2p1")]
    fn test_configuration_display(#[case] config: Configuration, #[case] expected: &str) {
        assert_eq!(format!("{}", config), expected);
    }
    #[rstest]
    #[case::h(e!(H), GroundState(Configuration::new(e!(H), None, occ!(s1))))]
    #[case::c(e!(C), GroundState(Configuration::new(e!(C), Some(e!(He)), occ!(s2p2))))]
    #[case::cr(e!(Cr), GroundState(Configuration::new(e!(Cr), Some(e!(Ar)), occ!(s1d5))))] // Exception
    #[case::mn(e!(Mn), GroundState(Configuration::new(e!(Mn), Some(e!(Ar)), occ!(s2d5))))]
    #[case::cu(e!(Cu), GroundState(Configuration::new(e!(Cu), Some(e!(Ar)), occ!(s1d10))))] // Exception
    #[case::la(e!(La), GroundState(Configuration::new(e!(La), Some(e!(Xe)), occ!(s2d1))))] // Exception
    #[case::eu(e!(Eu), GroundState(Configuration::new(e!(Eu), Some(e!(Xe)), occ!(s2f7))))]
    #[case::gd(e!(Gd), GroundState(Configuration::new(e!(Gd), Some(e!(Xe)), occ!(s2d1f7))))] // Exception
    #[case::xe(e!(Xe), GroundState(Configuration::new(e!(Xe), Some(e!(Kr)), occ!(s2p6d10))))]
    #[case::pb(e!(Pb), GroundState(Configuration::new(e!(Pb), Some(e!(Xe)), occ!(s2p2d10f14))))]
    fn test_ground_state_new(#[case] element: Element, #[case] expected: GroundState) {
        assert_eq!(GroundState::new(element), expected);
    }

    #[rstest]
    #[case::h(GroundState::new(Element::H), "s1")]
    #[case::he(GroundState::new(Element::He), "s2")]
    #[case::li(GroundState::new(Element::Li), "[He] s1")]
    #[case::be(GroundState::new(Element::Be), "[He] s2")]
    #[case::b(GroundState::new(Element::B), "[He] s2p1")]
    fn test_ground_state_display(#[case] gs: GroundState, #[case] expected: &str) {
        assert_eq!(format!("{}", gs), expected);
    }
}
