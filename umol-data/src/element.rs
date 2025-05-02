//! Element data

use convert_case::{Case, Casing};
use map_macro::hash_map;
use once_cell::sync::Lazy;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt::{self, Display};
use std::str::FromStr;
use umol::error::DataError;
use umol::{Error, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, PartialOrd, Ord)]
#[rustfmt::skip]
pub enum Element {
    H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca,
    Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr,
    Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba,
    La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W,
    Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th, Pa, U,
    Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Ds,
    Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og,
}

// Element data:
// 0. atomic number
// 1. atomic mass
// 2. symbol
// 3. period
// 4. group
// 5. valence electrons (includes all electrons in valence subshells)
// 6. max valence
// 7. max, min charge
// 8. max unpaired electrons
// 9. max implicit hydrogens
static ELEMENT_DATA: Lazy<
    HashMap<Element, (u8, f64, &'static str, u8, u8, u8, u8, (i8, i8), u8, u8)>,
> = Lazy::new(|| {
    hash_map! {
        Element::H => (1, 1.008, "H", 1, 1, 1, 2, (-1, 1), 1, 0),
        Element::He => (2, 4.002602, "He", 1, 32, 0, 0, (0, 0), 0, 0), // No valence electrons
        Element::Li => (3, 6.94, "Li", 2, 1, 1, 2, (-1, 1), 1, 0),
        Element::Be => (4, 9.0121831, "Be", 2, 2, 2, 2, (-2, 2), 2, 0),
        Element::B => (5, 10.81, "B", 2, 3, 27, 8, (-3, 3), 3, 3),
        Element::C => (6, 12.011, "C", 2, 28, 4, 8, (-4, 4), 4, 4),
        Element::N => (7, 14.007, "N", 2, 29, 5, 8, (-3, 5), 3, 3),
        Element::O => (8, 15.999, "O", 2, 30, 6, 8, (-2, 6), 2, 2),
        Element::F => (9, 18.998403163, "F", 2, 31, 7, 8, (-1, 7), 1, 1),
        Element::Ne => (10, 20.1797, "Ne", 2, 32, 0, 0, (0, 0), 0, 0), // No valence electrons
        Element::Na => (11, 22.98976928, "Na", 3, 1, 1, 2, (-1, 1), 1, 0),
        Element::Mg => (12, 24.305, "Mg", 3, 2, 2, 2, (0, 2), 2, 0),
        Element::Al => (13, 26.9815385, "Al", 3, 27, 3, 8, (-1, 3), 3, 0),
        Element::Si => (14, 28.085, "Si", 3, 28, 4, 8, (-4, 4), 4, 4),
        Element::P => (15, 30.973761998, "P", 3, 29, 5, 8, (-3, 5), 3, 3),
        Element::S => (16, 32.06, "S", 3, 30, 6, 8, (-2, 6), 2, 2),
        Element::Cl => (17, 35.45, "Cl", 3, 31, 7, 8, (-1, 7), 1, 1),
        Element::Ar => (18, 39.948, "Ar", 3, 32, 0, 0, (0, 0), 0, 0), // No valence electrons
        Element::K => (19, 39.0983, "K", 4, 1, 1, 2, (-1, 1), 1, 0),
        Element::Ca => (20, 40.078, "Ca", 4, 2, 2, 2, (0, 2), 2, 0),
        Element::Sc => (21, 44.955908, "Sc", 4, 3, 3, 12, (0, 3), 3, 0),
        Element::Ti => (22, 47.867, "Ti", 4, 18, 4, 18, (0, 4), 4, 0),
        Element::V => (23, 50.9415, "V", 4, 19, 5, 18, (0, 5), 5, 0),
        Element::Cr => (24, 51.9961, "Cr", 4, 20, 6, 18, (0, 6), 6, 0),
        Element::Mn => (25, 54.938044, "Mn", 4, 21, 7, 18, (0, 7), 7, 0),
        Element::Fe => (26, 55.845, "Fe", 4, 22,  8, 18, (0, 6), 6, 0),
        Element::Co => (27, 58.933194, "Co", 4, 23, 9, 18, (0, 5), 5, 0),
        Element::Ni => (28, 58.6934, "Ni", 4, 24, 10, 18, (0, 4), 4, 0),
        Element::Cu => (29, 63.546, "Cu", 4, 25, 11, 18, (0, 3), 3, 0),
        Element::Zn => (30, 65.38, "Zn", 4, 26, 12, 18, (0, 2), 2, 0),
        Element::Ga => (31, 69.723, "Ga", 4, 27, 13, 18, (0, 3), 3, 0),
        Element::Ge => (32, 72.63, "Ge", 4, 28, 14, 18, (0, 4), 4, 0),
        Element::As => (33, 74.921595, "As", 4, 29, 15, 18, (0, 3), 3, 3),
        Element::Se => (34, 78.971, "Se", 4, 30,16, 18, (0, 2), 2, 2),
        Element::Br => (35, 79.904, "Br", 4, 31, 17, 18, (0, 1), 1, 1),
        Element::Kr => (36, 83.798, "Kr", 4, 32, 18, 18, (0, 0), 0, 0), // Allow Kr compounds
        Element::Rb => (37, 85.4678, "Rb", 5, 1, 1, 2, (-1, 1), 1, 0),
        Element::Sr => (38, 87.62, "Sr", 5, 2, 2, 2, (0, 2), 2, 0),
        Element::Y => (39, 88.90584, "Y", 5, 3, 3, 12, (0, 3), 3, 0),
        Element::Zr => (40, 91.224, "Zr", 5, 18, 4, 18, (0, 4), 4, 0),
        Element::Nb => (41, 92.90637, "Nb", 5, 19, 5, 18, (0, 5), 5, 0),
        Element::Mo => (42, 95.95, "Mo", 5, 20, 6, 18, (0, 6), 6, 0),
        Element::Tc => (43, 98.0, "Tc", 5, 21, 7, 18, (0, 7), 7, 0),
        Element::Ru => (44, 101.07, "Ru", 5, 22, 8, 18, (0, 8), 8, 0),
        Element::Rh => (45, 102.90550, "Rh", 5, 23, 9, 18, (0, 6), 6, 0),
        Element::Pd => (46, 106.42, "Pd", 5, 24, 10, 18, (0, 5), 5, 0),
        Element::Ag => (47, 107.8682, "Ag", 5, 25, 11, 18, (0, 3), 3, 0),
        Element::Cd => (48, 112.414, "Cd", 5, 26, 12, 18, (0, 2), 2, 0),
        Element::In => (49, 114.818, "In", 5, 27, 13, 18, (0, 3), 3, 0),
        Element::Sn => (50, 118.710, "Sn", 5, 28, 14, 18, (0, 4), 4, 4),
        Element::Sb => (51, 121.760, "Sb", 5, 29, 15, 18, (0, 3), 3, 3),
        Element::Te => (52, 127.60, "Te", 5, 30,  16, 18, (0, 2), 2, 2),
        Element::I => (53, 126.90447, "I", 5, 31, 17, 18, (0, 1), 1, 1),
        Element::Xe => (54, 131.293, "Xe", 5, 32, 18, 18, (0, 0), 0, 0), // Allow Xe compounds
        Element::Cs => (55, 132.90545196, "Cs", 6, 1, 1, 2, (-1, 1), 1, 0),
        Element::Ba => (56, 137.327, "Ba", 6, 2, 2, 2, (0, 2), 2, 0),
        Element::La => (57, 138.90547, "La", 6, 3, 3, 18, (0, 3), 3, 0),
        Element::Ce => (58, 140.116, "Ce", 6, 4, 4, 20, (0, 4), 4, 0),
        Element::Pr => (59, 140.90766, "Pr", 6, 5, 5, 20, (0, 4), 4, 0),
        Element::Nd => (60, 144.242, "Nd", 6, 6, 6, 22, (-1, 3), 3, 0),
        Element::Pm => (61, 145.0, "Pm", 6, 7, 7, 22, (-1, 3), 3, 0),
        Element::Sm => (62, 150.36, "Sm", 6, 8, 8, 24, (-1, 3), 3, 0),
        Element::Eu => (63, 151.964, "Eu", 6, 9, 9, 24, (-1, 3), 3, 0),
        Element::Gd => (64, 157.25, "Gd", 6, 10, 10, 26, (-1, 3), 3, 0),
        Element::Tb => (65, 158.92535, "Tb", 6, 11, 11, 26, (-1, 4), 4, 0),
        Element::Dy => (66, 162.500, "Dy", 6, 12, 12, 28, (-1, 3), 3, 0),
        Element::Ho => (67, 164.93033, "Ho", 6, 13, 13, 28, (-1, 3), 3, 0),
        Element::Er => (68, 167.259, "Er", 6, 14, 14, 30, (-1, 3), 3, 0),
        Element::Tm => (69, 168.93422, "Tm", 6, 15, 15, 30, (-1, 3), 3, 0),
        Element::Yb => (70, 173.045, "Yb", 6, 16, 16, 32, (-1, 3), 3, 0),
        Element::Lu => (71, 174.9668, "Lu", 6, 17, 17, 32, (-1, 3), 3, 0),
        Element::Hf => (72, 178.49, "Hf", 6, 18, 18, 32, (0, 4), 4, 0),
        Element::Ta => (73, 180.94788, "Ta", 6, 19, 19, 32, (0, 5), 5, 0),
        Element::W => (74, 183.84, "W", 6, 20, 20, 32, (0, 6), 6, 0),
        Element::Re => (75, 186.207, "Re", 6, 21, 21, 32, (0, 7), 7, 0),
        Element::Os => (76, 190.23, "Os", 6, 22, 22, 32, (0, 8), 8, 0),
        Element::Ir => (77, 192.217, "Ir", 6, 23, 23, 32, (0, 6), 6, 0),
        Element::Pt => (78, 195.084, "Pt", 6, 24, 24, 32, (0, 6), 4, 0),
        Element::Au => (79, 196.966569, "Au", 6, 25, 25, 32, (0, 5), 1, 0),
        Element::Hg => (80, 200.592, "Hg", 6, 26, 26, 32, (0, 2), 2, 0),
        Element::Tl => (81, 204.38, "Tl", 6, 27, 27, 32, (0, 3), 3, 0),
        Element::Pb => (82, 207.2, "Pb", 6, 28, 28, 32, (0, 4), 4, 0),
        Element::Bi => (83, 208.98040, "Bi", 6, 29, 29, 32, (0, 3), 3, 3),
        Element::Po => (84, 209.0, "Po", 6, 30, 30, 32, (0, 2), 2, 2),
        Element::At => (85, 210.0, "At", 6, 31, 31, 32, (0, 7), 1, 1),
        Element::Rn => (86, 222.0, "Rn", 6, 32, 32, 32, (0, 8), 2, 0),
        Element::Fr => (87, 223.0, "Fr", 7, 1, 1, 2, (0, 1), 1, 0),
        Element::Ra => (88, 226.0, "Ra", 7, 2, 2, 2, (0, 2), 2, 0),
        Element::Ac => (89, 227.0, "Ac", 7, 3, 3, 18, (0, 3), 3, 0),
        Element::Th => (90, 232.0377, "Th", 7, 4, 4, 20, (0, 4), 4, 0),
        Element::Pa => (91, 231.03588, "Pa", 7, 5, 5, 20, (0, 5), 3, 0),
        Element::U => (92, 238.02891, "U", 7, 6, 6, 22, (0, 6), 4, 0),
        Element::Np => (93, 237.0, "Np", 7, 7, 7, 22, (0, 7), 5, 0),
        Element::Pu => (94, 244.0, "Pu", 7, 8, 8, 24, (0, 8), 6, 0),
        Element::Am => (95, 243.0, "Am", 7, 9, 9, 24, (0, 7), 7, 0),
        Element::Cm => (96, 247.0, "Cm", 7, 10, 10, 26, (0, 6), 8, 0),
        Element::Bk => (97, 247.0, "Bk", 7, 11, 11, 26, (0, 5), 5, 0),
        Element::Cf => (98, 251.0, "Cf", 7, 12, 12, 28, (0, 5), 4, 0),
        Element::Es => (99, 252.0, "Es", 7, 13, 13, 28, (0, 4), 3, 0),
        Element::Fm => (100, 257.0, "Fm", 7, 14, 14, 30, (0, 3), 2, 0),
        Element::Md => (101, 258.0, "Md", 7, 15, 15, 30, (0, 3), 1, 0),
        Element::No => (102, 259.0, "No", 7, 16, 16, 32, (0, 3), 0, 0),
        Element::Lr => (103, 262.0, "Lr", 7, 17, 17, 32, (0, 3), 1, 0),
        Element::Rf => (104, 267.0, "Rf", 7, 18, 18, 32, (0, 4), 2, 0),
        Element::Db => (105, 270.0, "Db", 7, 19,19, 32, (0, 5), 3, 0),
        Element::Sg => (106, 271.0, "Sg", 7, 20, 20, 32, (0, 6), 4, 0),
        Element::Bh => (107, 270.0, "Bh", 7, 21, 21, 32, (0, 7), 5, 0),
        Element::Hs => (108, 277.0, "Hs", 7, 22, 22, 32, (0, 8), 6, 0),
        Element::Mt => (109, 276.0, "Mt", 7, 23, 23, 32, (0, 6), 5, 0),
        Element::Ds => (110, 281.0, "Ds", 7, 24, 24, 32, (0, 6), 4, 0),
        Element::Rg => (111, 280.0, "Rg", 7, 25, 25, 32, (0, 5), 3, 0),
        Element::Cn => (112, 285.0, "Cn", 7, 26, 26, 32, (0, 4), 2, 0),
        Element::Nh => (113, 284.0, "Nh", 7, 27, 27, 32, (0, 3), 1, 0),
        Element::Fl => (114, 289.0, "Fl", 7, 28, 28, 32, (0, 2), 0, 0),
        Element::Mc => (115, 288.0, "Mc", 7, 29, 29, 32, (0, 0), 0, 0),
        Element::Lv => (116, 293.0, "Lv", 7, 30, 30, 32, (0, 0), 0, 0),
        Element::Ts => (117, 294.0, "Ts", 7, 31, 31, 32, (0, 0), 0, 0),
        Element::Og => (118, 294.0, "Og", 7, 32, 32, 32, (0, 0), 0, 0),
    }
});

/// Maximum atomic number, group number, and period number  
pub const MAX_ATOMIC_NUMBER: u8 = 118;
pub const MAX_PERIOD_NUMBER: u8 = 7;
pub const MIN_PERIOD_NUMBER: [u8; 33] = [
    0, 1, 2, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2,
    1,
];
pub const MAX_GROUP_NUMBER: [u8; 8] = [0, 2, 8, 8, 18, 18, 32, 32];

/// Element iterator
#[rustfmt::skip]
pub const ELEMENTS: [Element; 118] = [
    Element::H, Element::He, Element::Li, Element::Be, Element::B, Element::C,
    Element::N, Element::O, Element::F, Element::Ne, Element::Na, Element::Mg,
    Element::Al, Element::Si, Element::P, Element::S, Element::Cl, Element::Ar,
    Element::K, Element::Ca, Element::Sc, Element::Ti, Element::V, Element::Cr,
    Element::Mn, Element::Fe, Element::Co, Element::Ni, Element::Cu, Element::Zn,
    Element::Ga, Element::Ge, Element::As, Element::Se, Element::Br, Element::Kr,
    Element::Rb, Element::Sr, Element::Y, Element::Zr, Element::Nb, Element::Mo,
    Element::Tc, Element::Ru, Element::Rh, Element::Pd, Element::Ag, Element::Cd,
    Element::In, Element::Sn, Element::Sb, Element::Te, Element::I, Element::Xe,
    Element::Cs, Element::Ba, Element::La, Element::Ce, Element::Pr, Element::Nd,
    Element::Pm, Element::Sm, Element::Eu, Element::Gd, Element::Tb, Element::Dy,
    Element::Ho, Element::Er, Element::Tm, Element::Yb, Element::Lu, Element::Hf,
    Element::Ta, Element::W, Element::Re, Element::Os, Element::Ir, Element::Pt,
    Element::Au, Element::Hg, Element::Tl, Element::Pb, Element::Bi, Element::Po,
    Element::At, Element::Rn, Element::Fr, Element::Ra, Element::Ac, Element::Th,
    Element::Pa, Element::U, Element::Np, Element::Pu, Element::Am, Element::Cm,
    Element::Bk, Element::Cf, Element::Es, Element::Fm, Element::Md, Element::No,
    Element::Lr, Element::Rf, Element::Db, Element::Sg, Element::Bh, Element::Hs,
    Element::Mt, Element::Ds, Element::Rg, Element::Cn, Element::Nh, Element::Fl,
    Element::Mc, Element::Lv, Element::Ts, Element::Og,
];

/// Last element
pub const LAST_ELEMENT: Element = Element::Og;

static SYMBOL_TO_ELEMENT: Lazy<HashMap<&'static str, Element>> = Lazy::new(|| {
    ELEMENT_DATA
        .iter()
        .map(|(element, data)| (data.2, *element))
        .collect()
});

static ATOMIC_NUMBER_TO_ELEMENT: Lazy<HashMap<u8, Element>> = Lazy::new(|| {
    ELEMENT_DATA
        .iter()
        .map(|(element, data)| (data.0, *element))
        .collect()
});

static PERIOD_GROUP_TO_ELEMENT: Lazy<HashMap<(u8, u8), Element>> = Lazy::new(|| {
    ELEMENT_DATA
        .iter()
        .map(|(element, data)| ((data.3, data.4), *element))
        .collect()
});

impl Element {
    // Get element from symbol
    pub fn from_symbol(symbol: &str) -> Option<Self> {
        if symbol.is_empty() {
            return None;
        }

        // Handle case-insensitive lookup by normalizing capitalization
        let normalized_symbol = symbol.to_owned().to_case(Case::Title);
        SYMBOL_TO_ELEMENT.get(normalized_symbol.as_str()).copied()
    }

    // Get element from atomic number
    pub fn from_atomic_number(number: u8) -> Option<Self> {
        ATOMIC_NUMBER_TO_ELEMENT.get(&number).copied()
    }

    // Get element from period and group.
    // NOTE: 32-group layout is used for unique group assignment.
    pub fn from_period_group(period: u8, group: u8) -> Option<Self> {
        PERIOD_GROUP_TO_ELEMENT.get(&(period, group)).copied()
    }

    // Get atomic number for element
    pub fn atomic_number(&self) -> u8 {
        ELEMENT_DATA.get(self).unwrap().0
    }

    // Get atomic mass for element
    pub fn atomic_mass(&self) -> f64 {
        ELEMENT_DATA.get(self).unwrap().1
    }

    // Get symbol for element
    pub fn symbol(&self) -> &'static str {
        ELEMENT_DATA.get(self).unwrap().2
    }

    // Get period for element
    pub fn period(&self) -> u8 {
        ELEMENT_DATA.get(self).unwrap().3
    }

    // Get group number for element.
    // NOTE: 32-group layout is used for unique group assignment.
    pub fn group(&self) -> u8 {
        ELEMENT_DATA.get(self).unwrap().4
    }

    // Get number of valence electrons for element
    pub fn valence_electrons(&self) -> u8 {
        ELEMENT_DATA.get(self).unwrap().5
    }

    // Get max valence for element
    pub fn max_valence(&self) -> u8 {
        ELEMENT_DATA.get(self).unwrap().6
    }

    // Get minimum and maximum charge for element
    pub fn charge_bounds(&self) -> (i8, i8) {
        ELEMENT_DATA.get(self).unwrap().7
    }

    // Get maximum number of unpaired electrons for element
    pub fn max_unpaired_electrons(&self) -> u8 {
        ELEMENT_DATA.get(self).unwrap().8
    }

    // Get maximum number of implicit hydrogens for element
    pub fn max_implicit_hydrogens(&self) -> u8 {
        ELEMENT_DATA.get(self).unwrap().9
    }

    // Get next element in the periodic table
    pub fn next(&self) -> Option<Self> {
        if self == &LAST_ELEMENT {
            None
        } else {
            Self::from_atomic_number(self.atomic_number() + 1)
        }
    }

    // Get previous element in the periodic table
    pub fn previous(&self) -> Option<Self> {
        if self == &Element::H {
            None
        } else {
            Self::from_atomic_number(self.atomic_number() - 1)
        }
    }

    // Get element in the next period
    pub fn next_period(&self) -> Option<Self> {
        if self.period() == MAX_PERIOD_NUMBER {
            None
        } else {
            Self::from_period_group(self.period() + 1, self.group())
        }
    }

    // Get element in the previous period
    pub fn previous_period(&self) -> Option<Self> {
        if self.period() == MIN_PERIOD_NUMBER[self.group() as usize] {
            None
        } else {
            Self::from_period_group(self.period() - 1, self.group())
        }
    }

    // Get element in the next group
    pub fn next_group(&self) -> Option<Self> {
        if self.group() == MAX_GROUP_NUMBER[self.period() as usize] {
            None
        } else {
            Self::from_period_group(self.period(), self.group() + 1)
        }
    }

    // Get element in the previous group
    pub fn previous_group(&self) -> Option<Self> {
        if self.group() == 1 {
            None
        } else {
            Self::from_period_group(self.period(), self.group() - 1)
        }
    }
}

impl TryFrom<&str> for Element {
    type Error = Error;

    fn try_from(s: &str) -> Result<Self> {
        Self::from_symbol(s).ok_or_else(|| DataError::InvalidElement(s.to_string()).into())
    }
}

impl TryFrom<u8> for Element {
    type Error = Error;

    fn try_from(number: u8) -> Result<Self> {
        Self::from_atomic_number(number)
            .ok_or_else(|| DataError::InvalidElement(number.to_string()).into())
    }
}

impl FromStr for Element {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self> {
        Self::try_from(s)
    }
}

impl Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.symbol())
    }
}

/// Shorthand macro for element access
#[macro_export]
macro_rules! e {
    ($elem:ident) => {
        Element::$elem
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::*;
    use serde_json;

    #[test]
    fn test_element_from_symbol() {
        assert_eq!(Element::from_symbol("H"), Some(Element::H));
        assert_eq!(Element::from_symbol("He"), Some(Element::He));
        assert_eq!(Element::from_symbol("C"), Some(Element::C));
        assert_eq!(Element::from_symbol("invalid"), None);
    }

    #[test]
    fn test_element_from_atomic_number() {
        assert_eq!(Element::from_atomic_number(1), Some(Element::H));
        assert_eq!(Element::from_atomic_number(2), Some(Element::He));
        assert_eq!(Element::from_atomic_number(6), Some(Element::C));
        assert_eq!(Element::from_atomic_number(119), None);
    }

    #[test]
    fn test_element_from_period_group() {
        assert_eq!(Element::from_period_group(1, 1), Some(Element::H));
        assert_eq!(Element::from_period_group(2, 28), Some(Element::C));
        assert_eq!(Element::from_period_group(1, 3), None);
    }

    #[rstest]
    #[case(Element::H, "H", 1, 1.008, 1, 1, 1, 2)]
    #[case(Element::He, "He", 2, 4.002602, 1, 32, 0, 0)]
    fn test_element_properties(
        #[case] element: Element,
        #[case] symbol: &str,
        #[case] atomic_number: u8,
        #[case] atomic_mass: f64,
        #[case] period: u8,
        #[case] group: u8,
        #[case] valence_electrons: u8,
        #[case] max_valence: u8,
    ) {
        assert_eq!(element.symbol(), symbol);
        assert_eq!(element.atomic_number(), atomic_number);
        assert!((element.atomic_mass() - atomic_mass).abs() < 1e-10);
        assert_eq!(element.period(), period);
        assert_eq!(element.group(), group);
        assert_eq!(element.valence_electrons(), valence_electrons);
        assert_eq!(element.max_valence(), max_valence);
    }

    #[rstest]
    #[case(Element::H, -1, 1, 1, 0)]
    #[case(Element::He, 0, 0, 0, 0)]
    #[case(Element::C, -4, 4, 4, 4)]
    #[case(Element::N, -3, 5, 3, 3)]
    fn test_element_bounds(
        #[case] element: Element,
        #[case] min_charge: i8,
        #[case] max_charge: i8,
        #[case] max_unpaired: u8,
        #[case] max_implicit_hydrogens: u8,
    ) {
        let (actual_min, actual_max) = element.charge_bounds();
        assert_eq!(actual_min, min_charge);
        assert_eq!(actual_max, max_charge);
        assert_eq!(element.max_unpaired_electrons(), max_unpaired);
        assert_eq!(element.max_implicit_hydrogens(), max_implicit_hydrogens);
    }

    #[test]
    fn test_element_display() {
        assert_eq!(Element::H.to_string(), "H");
        assert_eq!(Element::He.to_string(), "He");
        assert_eq!(Element::Li.to_string(), "Li");
        assert_eq!(Element::Be.to_string(), "Be");
    }

    #[test]
    fn test_all_elements_have_data() {
        for element in ELEMENT_DATA.keys() {
            let data = ELEMENT_DATA.get(element).unwrap();
            assert!(data.0 > 0); // atomic number
            assert!(data.1 > 0.0); // atomic mass
            assert!(!data.2.is_empty()); // symbol
            assert!(data.3 <= 7); // period
            assert!(data.4 <= 32); // group
            assert!(data.5 <= 32); // valence electrons
            assert!(data.6 <= 32); // max valence
            assert!(data.7 .0 <= data.7 .1); // charge bounds
            assert!(data.8 <= 10); // max unpaired electrons
            assert!(data.9 <= 4); // max implicit hydrogens
        }
    }

    #[test]
    fn test_symbol_to_element_to_symbol() {
        for (symbol, element) in SYMBOL_TO_ELEMENT.iter() {
            assert_eq!(element.symbol(), *symbol);
        }
    }

    #[test]
    fn test_atomic_number_to_element_to_atomic_number() {
        for (number, element) in ATOMIC_NUMBER_TO_ELEMENT.iter() {
            assert_eq!(element.atomic_number(), *number);
        }
    }

    #[test]
    fn test_element_ordering() {
        // Test basic ordering
        assert!(Element::H < Element::He);
        assert!(Element::He < Element::Li);
        assert!(Element::C < Element::N);
        assert!(Element::N < Element::O);

        // Test equality
        assert!(Element::H == Element::H);
        assert!(Element::C == Element::C);

        // Test partial ordering
        assert!(Element::H <= Element::H);
        assert!(Element::H <= Element::He);
        assert!(Element::He >= Element::H);
    }

    #[test]
    fn test_element_next_previous() {
        assert_eq!(Element::H.previous(), None);
        assert_eq!(Element::H.next(), Some(Element::He));
        assert_eq!(Element::He.previous(), Some(Element::H));
        assert_eq!(Element::Og.previous(), Some(Element::Ts));
        assert_eq!(Element::Og.next(), None);
    }

    #[test]
    fn test_element_next_previous_period() {
        assert_eq!(Element::H.previous_period(), None);
        assert_eq!(Element::H.next_period(), Some(Element::Li));
        assert_eq!(Element::C.previous_period(), None);
        assert_eq!(Element::C.next_period(), Some(Element::Si));
        // assert_eq!(Element::U.previous_period(), Some(Element::Nd));
        // assert_eq!(Element::U.next_period(), None);
    }

    #[test]
    fn test_element_next_previous_group() {
        assert_eq!(Element::Li.previous_group(), None);
        // assert_eq!(Element::Li.next_group(), Some(Element::Be));
        // assert_eq!(Element::Ne.previous_group(), Some(Element::F));
        // assert_eq!(Element::Ne.next_group(), None);
    }

    #[test]
    fn test_element_serialization() {
        // Test serialization of individual elements
        assert_eq!(serde_json::to_string(&Element::H).unwrap(), r#""H""#);
        assert_eq!(serde_json::to_string(&Element::He).unwrap(), r#""He""#);
        assert_eq!(serde_json::to_string(&Element::C).unwrap(), r#""C""#);

        // Test serialization of a vector of elements
        let elements = vec![Element::H, Element::He, Element::C];
        assert_eq!(
            serde_json::to_string(&elements).unwrap(),
            r#"["H","He","C"]"#
        );
    }

    #[test]
    fn test_element_deserialization() {
        // Test deserialization of individual elements
        assert_eq!(
            serde_json::from_str::<Element>(r#""H""#).unwrap(),
            Element::H
        );
        assert_eq!(
            serde_json::from_str::<Element>(r#""He""#).unwrap(),
            Element::He
        );
        assert_eq!(
            serde_json::from_str::<Element>(r#""C""#).unwrap(),
            Element::C
        );

        // Test deserialization of a vector of elements
        let elements: Vec<Element> = serde_json::from_str(r#"["H","He","C"]"#).unwrap();
        assert_eq!(elements, vec![Element::H, Element::He, Element::C]);

        // Test error handling for invalid element symbols
        assert!(serde_json::from_str::<Element>(r#""Invalid""#).is_err());
        assert!(serde_json::from_str::<Element>(r#""123""#).is_err());
    }

    #[test]
    fn test_element_roundtrip() {
        // Test roundtrip serialization/deserialization
        let elements = vec![Element::H, Element::He, Element::C, Element::N, Element::O];

        let serialized = serde_json::to_string(&elements).unwrap();
        let deserialized: Vec<Element> = serde_json::from_str(&serialized).unwrap();

        assert_eq!(elements, deserialized);
    }

    #[test]
    fn test_element_macro() {
        assert_eq!(e!(H), Element::H);
        assert_eq!(e!(Fe), Element::Fe);
    }
}
