//! Electronic configuration
use crate::{Element, Occupation, SpinState, MAX_UNPAIRED_ELECTRONS};
use std::cmp;
use std::fmt::{self, Display};

/// Electronic configuration of atom or ion, including occupation and spin state
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Configuration {
    element: Element,
    core_element: Option<Element>,
    valence_occupation: Occupation,
    spin_state: SpinState,
}

impl Configuration {
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
        self.core_element
            .map(|core| compute_regular_configuration(core.atomic_number()).1)
    }

    /// Return computed valence occupation
    pub fn valence_occupation(&self) -> Occupation {
        self.valence_occupation
    }

    /// Return computed spin state
    pub fn spin_state(&self) -> SpinState {
        self.spin_state
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
            "{}: {}{} ({})",
            self.element.symbol(),
            self.core_element()
                .map(|core| format!("[{}] ", core.symbol()))
                .unwrap_or("".to_string()),
            self.valence_occupation(),
            self.spin_state().name()
        )
    }
}

/// Maximum n quantum number
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

/// Compute configurations from Aufbau principle, Madelung rule, and Hund's rules
fn compute_regular_configuration(electron_count: u8) -> (Option<Element>, Occupation, SpinState) {
    let mut remaining = electron_count;
    let mut subshell_occupations = Vec::new();
    let mut closing_subshell = (0, 0, 0, false);

    for subshell @ (n, l, capacity, is_closing) in MADELUNG_ORDER {
        println!("n: {}, l: {}, capacity: {}", n, l, capacity);

        let occupation = cmp::min(remaining, capacity);
        subshell_occupations.push((n, l, capacity, occupation));
        remaining -= occupation;
        if is_closing && remaining != 0 {
            closing_subshell = subshell;
        }
        println!("subshell_occupations: {:?}", subshell_occupations);
        println!("remaining: {}", remaining);

        if remaining == 0 {
            break;
        }
    }
    println!("closing_subshell: {:?}", closing_subshell);

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
    println!("unpaired: {}", unpaired);
    println!("subshell_occupations: {:?}", subshell_occupations);
    debug_assert!(
        subshell_occupations
            .iter()
            .map(|counts| counts.3)
            .sum::<u8>()
            == electron_count
    );
    debug_assert!(unpaired <= MAX_UNPAIRED_ELECTRONS);
    debug_assert!(closing_subshell.0 <= MAX_N_QUANTUM_NUMBER);
    debug_assert!(closing_subshell.1 <= MAX_L_QUANTUM_NUMBER);

    let valence_occupation = compute_valence_occupation(subshell_occupations, closing_subshell);
    let core_element = compute_core_element(closing_subshell);
    let spin_state = SpinState::from_unpaired_electrons(unpaired).unwrap();
    println!("core_element: {:?}", core_element);
    println!("valence_occupation: {:?}", valence_occupation);
    println!("spin_state: {:?}", spin_state);
    (core_element, valence_occupation, spin_state)
}

/// Compute core element
fn compute_core_element(closing_subshell: (u8, u8, u8, bool)) -> Option<Element> {
    match closing_subshell {
        (0, 0, _, false) => None,
        (1, 0, _, true) => Some(Element::He),
        (2, 1, _, true) => Some(Element::Ne),
        (3, 1, _, true) => Some(Element::Ar),
        (4, 1, _, true) => Some(Element::Kr),
        (5, 1, _, true) => Some(Element::Xe),
        (6, 1, _, true) => Some(Element::Rn),
        (7, 1, _, true) => Some(Element::Og),
        _ => unreachable!(),
    }
}

/// Compute valence occupation
fn compute_valence_occupation(
    subshell_counts: Vec<(u8, u8, u8, u8)>,
    closing_subshell: (u8, u8, u8, bool),
) -> Occupation {
    let mut valence_counts: [u8; 4] = [0; 4];

    let mut in_valence = closing_subshell == (0, 0, 0, false);
    println!(
        "before closing_subshell: {:?}, in_valence: {}",
        closing_subshell, in_valence
    );
    for (n, l, _, count) in subshell_counts {
        if count == 0 {
            break;
        }

        println!(
            "n: {}, l: {}, count: {}, in_valence: {}",
            n, l, count, in_valence
        );

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
    use super::*;
    use rstest::*;

    #[rstest]
    #[case(1, (None, Occupation::new(1, 0, 0, 0), SpinState::Doublet))]
    #[case(2, (None, Occupation::new(2, 0, 0, 0), SpinState::Singlet))]
    #[case(3, (Some(Element::He), Occupation::new(1, 0, 0, 0), SpinState::Doublet))]
    #[case(4, (Some(Element::He), Occupation::new(2, 0, 0, 0), SpinState::Singlet))]
    #[case(5, (Some(Element::He), Occupation::new(2, 1, 0, 0), SpinState::Doublet))]
    #[case(6, (Some(Element::He), Occupation::new(2, 2, 0, 0), SpinState::Triplet))]
    #[case(7, (Some(Element::He), Occupation::new(2, 3, 0, 0), SpinState::Quartet))]
    #[case(8, (Some(Element::He), Occupation::new(2, 4, 0, 0), SpinState::Triplet))]
    #[case(9, (Some(Element::He), Occupation::new(2, 5, 0, 0), SpinState::Doublet))]
    #[case(10, (Some(Element::He), Occupation::new(2, 6, 0, 0), SpinState::Singlet))]
    #[case(11, (Some(Element::Ne), Occupation::new(1, 0, 0, 0), SpinState::Doublet))]
    #[case(12, (Some(Element::Ne), Occupation::new(2, 0, 0, 0), SpinState::Singlet))]
    #[case(13, (Some(Element::Ne), Occupation::new(2, 1, 0, 0), SpinState::Doublet))]
    #[case(14, (Some(Element::Ne), Occupation::new(2, 2, 0, 0), SpinState::Triplet))]
    #[case(15, (Some(Element::Ne), Occupation::new(2, 3, 0, 0), SpinState::Quartet))]
    #[case(16, (Some(Element::Ne), Occupation::new(2, 4, 0, 0), SpinState::Triplet))]
    #[case(17, (Some(Element::Ne), Occupation::new(2, 5, 0, 0), SpinState::Doublet))]
    #[case(18, (Some(Element::Ne), Occupation::new(2, 6, 0, 0), SpinState::Singlet))]
    #[case(19, (Some(Element::Ar), Occupation::new(1, 0, 0, 0), SpinState::Doublet))]
    #[case(20, (Some(Element::Ar), Occupation::new(2, 0, 0, 0), SpinState::Singlet))]
    #[case(21, (Some(Element::Ar), Occupation::new(2, 0, 1, 0), SpinState::Doublet))]
    #[case(22, (Some(Element::Ar), Occupation::new(2, 0, 2, 0), SpinState::Triplet))]
    #[case(23, (Some(Element::Ar), Occupation::new(2, 0, 3, 0), SpinState::Quartet))]
    #[case(24, (Some(Element::Ar), Occupation::new(2, 0, 4, 0), SpinState::Quintet))] // Exception
    #[case(25, (Some(Element::Ar), Occupation::new(2, 0, 5, 0), SpinState::Sextet))]
    #[case(26, (Some(Element::Ar), Occupation::new(2, 0, 6, 0), SpinState::Quintet))]
    #[case(27, (Some(Element::Ar), Occupation::new(2, 0, 7, 0), SpinState::Quartet))]
    #[case(28, (Some(Element::Ar), Occupation::new(2, 0, 8, 0), SpinState::Triplet))]
    #[case(29, (Some(Element::Ar), Occupation::new(2, 0, 9, 0), SpinState::Doublet))] // Exception
    #[case(30, (Some(Element::Ar), Occupation::new(2, 0, 10, 0), SpinState::Singlet))]
    fn test_compute_regular_configuration(
        #[case] z: u8,
        #[case] expected: (Option<Element>, Occupation, SpinState),
    ) {
        assert_eq!(compute_regular_configuration(z), expected);
    }

    // #[rstest]
    // #[case(Configuration, "H s1 (doublet)")]
    // #[case(Configuration, "s2 (singlet)")]
    // #[case(Configuration, "[He] s1 (doublet)")]
    // #[case(Configuration, "[He] s2 (singlet)")]
    // #[case(Configuration, "[He] s2p1 (doublet)")]
    // fn test_display(#[case] config: Configuration, #[case] expected: &str) {
    //     assert_eq!(format!("{}", config), expected);
    // }
}
