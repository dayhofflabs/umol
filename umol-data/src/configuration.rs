//! Electronic configuration
use crate::{Element, Occupation, SpinState, MAX_ATOMIC_NUMBER, MAX_UNPAIRED_ELECTRONS};
use once_cell::sync::Lazy;
use std::cmp;
use std::collections::HashMap;
use std::fmt::{self, Display};
use std::sync::Mutex;

/// Electronic configuration of atom or ion, including occupation and spin state
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Configuration {
    /// Ground state configuration of neutral atom with nuclear charge Z
    GroundState(u8),
    /// Ground state of atomic ion with nuclear charge Z and charge q
    Ion(u8, i8),
}

impl Configuration {
    /// Return core element
    pub fn core_element(&self) -> Option<Element> {
        compute_property(*self).0
    }

    /// Return core occupation
    pub fn core_occupation(&self) -> Option<Occupation> {
        compute_property(*self)
            .0
            .map(|core| compute_property(Configuration::GroundState(core.atomic_number())).1)
    }

    /// Return computed valence occupation
    /// Uses memoization for performance
    pub fn valence_occupation(&self) -> Occupation {
        compute_property(*self).1
    }

    /// Return computed spin state
    /// Uses memoization for performance
    pub fn spin_state(&self) -> SpinState {
        compute_property(*self).2
    }

    /// Return computed electron count
    pub fn electron_count(&self) -> u8 {
        match self {
            Configuration::GroundState(z) => *z,
            Configuration::Ion(z, charge) => (*z as i8 - charge) as u8,
        }
    }

    /// Return atomic number
    pub fn atomic_number(&self) -> u8 {
        match self {
            Configuration::GroundState(z) => *z,
            Configuration::Ion(z, _) => *z,
        }
    }

    /// Return charge
    pub fn charge(&self) -> i8 {
        match self {
            Configuration::GroundState(_) => 0,
            Configuration::Ion(_, charge) => *charge,
        }
    }
}

impl Display for Configuration {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Configuration::GroundState(_) => write!(
                f,
                "{}{} ({})",
                self.core_element()
                    .map(|core| format!("[{}] ", core.symbol()))
                    .unwrap_or("".to_string()),
                self.valence_occupation(),
                self.spin_state().name()
            ),
            Configuration::Ion(z, charge) => write!(f, "{}-{}", z, charge),
        }
    }
}

/// Maximum n quantum number
pub const MAX_N_QUANTUM_NUMBER: u8 = 7;
/// Maximum l quantum number
pub const MAX_L_QUANTUM_NUMBER: u8 = 3;

/// Cached configuration properties
type CacheKey = (u8, i8);
type CacheValue = (Option<Element>, Occupation, SpinState);
static CACHED_PROPERTIES: Lazy<Mutex<HashMap<CacheKey, CacheValue>>> =
    Lazy::new(|| Mutex::new(HashMap::new()));

// Order of filling of atomic orbitals according to Madelung rule
// (n, l, capacity, closing subshell)
#[rustfmt::skip]
static MADELUNG_ORDER: [(u8, u8, u8, bool); 19] = [
    (1, 0, 2, true), (2, 0, 2, false), (2, 1, 6, true), (3, 0, 2, false), (3, 1, 6, true),
    (4, 0, 2, false), (3, 2, 10, false), (4, 1, 6, true), (5, 0, 2, false), (4, 2, 10, false),
    (5, 1, 6, true), (6, 0, 2, false), (4, 3, 14, false), (5, 2, 10, false), (6, 1, 6, true),
    (7, 0, 2, false), (5, 3, 14, false), (6, 2, 10, false), (7, 1, 6, true)];

/// Compute configuration properties
/// Uses memoization for performance
fn compute_property(config: Configuration) -> CacheValue {
    let key = (config.atomic_number(), config.charge());

    // Check cache
    {
        let cache = CACHED_PROPERTIES.lock().unwrap();
        if let Some(cached) = cache.get(&key) {
            return *cached;
        }
    }

    // Perform computation
    let value = compute_property_impl(key.0, key.1);

    // Insert into cache
    {
        let mut cache = CACHED_PROPERTIES.lock().unwrap();
        cache.insert(key, value);
    }

    value
}

/// Compute configuration properties using a set of rules + exceptions
fn compute_property_impl(z: u8, charge: i8) -> CacheValue {
    println!("compute_property_impl: z: {}, charge: {}", z, charge);
    debug_assert!(z <= MAX_ATOMIC_NUMBER && charge <= z as i8);

    let electron_count = (z as i8 - charge) as u8;

    // Process exceptions for neutral ground states
    if charge == 0 {
        match z {
            24 => {
                return (
                    Some(Element::Ar),
                    Occupation::new(1, 0, 5, 0),
                    SpinState::Septet,
                )
            }
            29 => {
                return (
                    Some(Element::Ar),
                    Occupation::new(1, 0, 10, 0),
                    SpinState::Doublet,
                )
            }
            _ => compute_property_from_rules(electron_count),
        }
    } else {
        compute_property_from_rules(electron_count)
    }
}

/// Compute configuration properties using Aufbau principle, Madelung rule, and
/// Hund's rules
fn compute_property_from_rules(electron_count: u8) -> CacheValue {
    let mut remaining_count = electron_count;
    let mut subshell_counts = HashMap::new();

    // Track highest energy subshell with electrons
    let mut last_filled_nl = (0, 0);
    // Track last filled p-subshell (for determining core element)
    let mut last_closing_nl = (0, 0);

    for (n, l, capacity, closing) in MADELUNG_ORDER {
        if remaining_count == 0 {
            break;
        }
        println!("n: {}, l: {}, capacity: {}", n, l, capacity);

        let fill = cmp::min(remaining_count, capacity);
        *subshell_counts.entry((n, l)).or_insert(0) += fill;
        remaining_count -= fill;
        if fill > 0 {
            last_filled_nl = (n, l);
            if closing && remaining_count != 0 {
                last_closing_nl = (n, l);
            }
        }
        println!("subshell_counts: {:?}", subshell_counts);
        println!("remaining_count: {}", remaining_count);
    }
    println!("last_filled_nl: {:?}", last_filled_nl);
    println!("last_closing_nl: {:?}", last_closing_nl);

    // Add unpaired electrons from last partially filled subshell
    let mut unpaired_electrons = 0;
    if let Some(&filled_count) = subshell_counts.get(&last_filled_nl) {
        let (_n, l) = last_filled_nl;
        let num_orbitals = 2 * l + 1;
        let capacity = 2 * num_orbitals;

        if filled_count < capacity {
            if filled_count <= num_orbitals {
                // Less than half filled shell
                unpaired_electrons = filled_count;
            } else {
                // More than half filled shell
                unpaired_electrons = capacity - filled_count;
            }
        }
    }
    println!("unpaired_electrons: {}", unpaired_electrons);
    println!("subshell_counts: {:?}", subshell_counts);
    debug_assert!(subshell_counts.values().sum::<u8>() == electron_count);
    debug_assert!(unpaired_electrons <= MAX_UNPAIRED_ELECTRONS);
    debug_assert!(last_closing_nl.0 <= MAX_N_QUANTUM_NUMBER);
    debug_assert!(last_closing_nl.1 <= MAX_L_QUANTUM_NUMBER);

    let (core_occupation, valence_occupation) =
        compute_core_valence_occupation(subshell_counts, last_closing_nl);
    let core_element = compute_core_element(core_occupation);
    let spin_state = SpinState::from_unpaired_electrons(unpaired_electrons).unwrap();
    println!("core_element: {:?}", core_element);
    println!("core_occupation: {:?}", core_occupation);
    println!("valence_occupation: {:?}", valence_occupation);
    println!("spin_state: {:?}", spin_state);
    (core_element, valence_occupation, spin_state)
}

/// Compute core element
fn compute_core_element(core_occupation: Occupation) -> Option<Element> {
    match (
        core_occupation.s(),
        core_occupation.p(),
        core_occupation.d(),
        core_occupation.f(),
    ) {
        (0, 0, 0, 0) => None,
        (2, 0, 0, 0) => Some(Element::He),
        (4, 6, 0, 0) => Some(Element::Ne),
        (6, 12, 0, 0) => Some(Element::Ar),
        (8, 18, 10, 0) => Some(Element::Kr),
        (10, 24, 20, 0) => Some(Element::Xe),
        (12, 30, 30, 14) => Some(Element::Rn),
        (14, 36, 40, 28) => Some(Element::Og),
        _ => unreachable!(),
    }
}

/// Compute core and valence occupation
fn compute_core_valence_occupation(
    subshell_counts: HashMap<(u8, u8), u8>,
    last_closing_nl: (u8, u8),
) -> (Occupation, Occupation) {
    let mut core_counts: [u8; 4] = [0; 4];
    let mut valence_counts: [u8; 4] = [0; 4];

    let mut in_valence = last_closing_nl == (0, 0);
    println!(
        "before last_closing_nl: {:?}, in_valence: {}",
        last_closing_nl, in_valence
    );
    for (n, l, _, _) in MADELUNG_ORDER {
        let count = subshell_counts.get(&(n, l)).unwrap_or(&0);
        println!(
            "n: {}, l: {}, count: {}, in_valence: {}",
            n, l, count, in_valence
        );
        if *count == 0 {
            break;
        }

        if in_valence {
            valence_counts[l as usize] += *count;
        } else {
            core_counts[l as usize] += *count;
        }

        if !in_valence && (n, l) == last_closing_nl {
            in_valence = true;
        }
    }
    (
        Occupation::new(
            core_counts[0],
            core_counts[1],
            core_counts[2],
            core_counts[3],
        ),
        Occupation::new(
            valence_counts[0],
            valence_counts[1],
            valence_counts[2],
            valence_counts[3],
        ),
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
    #[case(24, (Some(Element::Ar), Occupation::new(1, 0, 5, 0), SpinState::Septet))]
    #[case(25, (Some(Element::Ar), Occupation::new(2, 0, 5, 0), SpinState::Sextet))]
    #[case(26, (Some(Element::Ar), Occupation::new(2, 0, 6, 0), SpinState::Quintet))]
    #[case(27, (Some(Element::Ar), Occupation::new(2, 0, 7, 0), SpinState::Quartet))]
    #[case(28, (Some(Element::Ar), Occupation::new(2, 0, 8, 0), SpinState::Triplet))]
    #[case(29, (Some(Element::Ar), Occupation::new(1, 0, 10, 0), SpinState::Doublet))]
    #[case(30, (Some(Element::Ar), Occupation::new(2, 0, 10, 0), SpinState::Singlet))]
    fn test_compute_property_impl(#[case] z: u8, #[case] expected: CacheValue) {
        assert_eq!(compute_property_impl(z, 0), expected);
    }

    #[rstest]
    #[case(Configuration::GroundState(1), "s1 (doublet)")]
    #[case(Configuration::GroundState(2), "s2 (singlet)")]
    #[case(Configuration::GroundState(3), "[He] s1 (doublet)")]
    #[case(Configuration::GroundState(4), "[He] s2 (singlet)")]
    #[case(Configuration::GroundState(5), "[He] s2p1 (doublet)")]
    fn test_display(#[case] config: Configuration, #[case] expected: &str) {
        assert_eq!(format!("{}", config), expected);
    }
}
