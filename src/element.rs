// Element data and validation

use crate::core::{Error, Result};
use crate::error::DataError;
use map_macro::hash_map;
use once_cell::sync::Lazy;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt::{self, Display};

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

static ELEMENT_DATA: Lazy<HashMap<Element, (u8, f64, &'static str, (i8, i8), u8)>> =
    Lazy::new(|| {
        hash_map! {
            Element::H => (1, 1.008, "H", (-1, 1), 1),
            Element::He => (2, 4.002602, "He", (0, 0), 0),
            Element::Li => (3, 6.94, "Li", (-1, 1), 1),
            Element::Be => (4, 9.0121831, "Be", (-2, 2), 2),
            Element::B => (5, 10.81, "B", (-3, 3), 3),
            Element::C => (6, 12.011, "C", (-4, 4), 4),
            Element::N => (7, 14.007, "N", (-3, 5), 3),
            Element::O => (8, 15.999, "O", (-2, 2), 2),
            Element::F => (9, 18.998403163, "F", (-1, 1), 1),
            Element::Ne => (10, 20.1797, "Ne", (0, 0), 0),
            Element::Na => (11, 22.98976928, "Na", (-1, 1), 1),
            Element::Mg => (12, 24.305, "Mg", (0, 2), 2),
            Element::Al => (13, 26.9815385, "Al", (-1, 3), 3),
            Element::Si => (14, 28.085, "Si", (-4, 4), 4),
            Element::P => (15, 30.973761998, "P", (-3, 5), 3),
            Element::S => (16, 32.06, "S", (-2, 6), 2),
            Element::Cl => (17, 35.45, "Cl", (-1, 7), 1),
            Element::Ar => (18, 39.948, "Ar", (0, 0), 0),
            Element::K => (19, 39.0983, "K", (-1, 1), 1),
            Element::Ca => (20, 40.078, "Ca", (0, 2), 2),
            Element::Sc => (21, 44.955908, "Sc", (-1, 3), 3),
            Element::Ti => (22, 47.867, "Ti", (-1, 4), 4),
            Element::V => (23, 50.9415, "V", (-1, 5), 5),
            Element::Cr => (24, 51.9961, "Cr", (-1, 6), 6),
            Element::Mn => (25, 54.938044, "Mn", (0, 7), 7),
            Element::Fe => (26, 55.845, "Fe", (-1, 6), 6),
            Element::Co => (27, 58.933194, "Co", (-1, 5), 5),
            Element::Ni => (28, 58.6934, "Ni", (-1, 4), 4),
            Element::Cu => (29, 63.546, "Cu", (-1, 3), 3),
            Element::Zn => (30, 65.38, "Zn", (0, 2), 2),
            Element::Ga => (31, 69.723, "Ga", (-1, 3), 3),
            Element::Ge => (32, 72.63, "Ge", (-4, 4), 4),
            Element::As => (33, 74.921595, "As", (-3, 5), 3),
            Element::Se => (34, 78.971, "Se", (-2, 6), 2),
            Element::Br => (35, 79.904, "Br", (-1, 7), 1),
            Element::Kr => (36, 83.798, "Kr", (0, 8), 2),
            Element::Rb => (37, 85.4678, "Rb", (-1, 1), 1),
            Element::Sr => (38, 87.62, "Sr", (0, 2), 2),
            Element::Y => (39, 88.90584, "Y", (-1, 3), 3),
            Element::Zr => (40, 91.224, "Zr", (-1, 4), 4),
            Element::Nb => (41, 92.90637, "Nb", (-1, 5), 5),
            Element::Mo => (42, 95.95, "Mo", (-1, 6), 6),
            Element::Tc => (43, 98.0, "Tc", (-1, 7), 7),
            Element::Ru => (44, 101.07, "Ru", (-1, 8), 6),
            Element::Rh => (45, 102.90550, "Rh", (-1, 6), 3),
            Element::Pd => (46, 106.42, "Pd", (-1, 4), 4),
            Element::Ag => (47, 107.8682, "Ag", (-1, 3), 3),
            Element::Cd => (48, 112.414, "Cd", (0, 2), 2),
            Element::In => (49, 114.818, "In", (-1, 3), 3),
            Element::Sn => (50, 118.710, "Sn", (-2, 4), 4),
            Element::Sb => (51, 121.760, "Sb", (-3, 5), 3),
            Element::Te => (52, 127.60, "Te", (-2, 6), 4),
            Element::I => (53, 126.90447, "I", (-1, 7), 1),
            Element::Xe => (54, 131.293, "Xe", (0, 8), 2),
            Element::Cs => (55, 132.90545196, "Cs", (-1, 1), 1),
            Element::Ba => (56, 137.327, "Ba", (0, 2), 2),
            Element::La => (57, 138.90547, "La", (-1, 3), 3),
            Element::Ce => (58, 140.116, "Ce", (-1, 4), 4),
            Element::Pr => (59, 140.90766, "Pr", (-1, 4), 5),
            Element::Nd => (60, 144.242, "Nd", (-1, 3), 6),
            Element::Pm => (61, 145.0, "Pm", (-1, 3), 5),
            Element::Sm => (62, 150.36, "Sm", (-1, 3), 6),
            Element::Eu => (63, 151.964, "Eu", (-1, 3), 9),
            Element::Gd => (64, 157.25, "Gd", (-1, 3), 10),
            Element::Tb => (65, 158.92535, "Tb", (-1, 4), 9),
            Element::Dy => (66, 162.500, "Dy", (-1, 3), 6),
            Element::Ho => (67, 164.93033, "Ho", (-1, 3), 5),
            Element::Er => (68, 167.259, "Er", (-1, 3), 6),
            Element::Tm => (69, 168.93422, "Tm", (-1, 3), 1),
            Element::Yb => (70, 173.045, "Yb", (-1, 3), 2),
            Element::Lu => (71, 174.9668, "Lu", (-1, 3), 3),
            Element::Hf => (72, 178.49, "Hf", (0, 4), 4),
            Element::Ta => (73, 180.94788, "Ta", (-1, 5), 5),
            Element::W => (74, 183.84, "W", (-1, 6), 6),
            Element::Re => (75, 186.207, "Re", (-1, 7), 7),
            Element::Os => (76, 190.23, "Os", (-1, 8), 6),
            Element::Ir => (77, 192.217, "Ir", (-1, 6), 5),
            Element::Pt => (78, 195.084, "Pt", (-1, 6), 4),
            Element::Au => (79, 196.966569, "Au", (-1, 3), 4),
            Element::Hg => (80, 200.592, "Hg", (0, 2), 2),
            Element::Tl => (81, 204.38, "Tl", (-1, 3), 3),
            Element::Pb => (82, 207.2, "Pb", (-2, 4), 4),
            Element::Bi => (83, 208.98040, "Bi", (-3, 3), 3),
            Element::Po => (84, 209.0, "Po", (-2, 6), 2),
            Element::At => (85, 210.0, "At", (-1, 7), 1),
            Element::Rn => (86, 222.0, "Rn", (0, 8), 2),
            Element::Fr => (87, 223.0, "Fr", (-1, 1), 1),
            Element::Ra => (88, 226.0, "Ra", (0, 2), 2),
            Element::Ac => (89, 227.0, "Ac", (0, 3), 3),
            Element::Th => (90, 232.0377, "Th", (0, 4), 4),
            Element::Pa => (91, 231.03588, "Pa", (0, 5), 3),
            Element::U => (92, 238.02891, "U", (0, 6), 4),
            Element::Np => (93, 237.0, "Np", (0, 7), 5),
            Element::Pu => (94, 244.0, "Pu", (0, 8), 6),
            Element::Am => (95, 243.0, "Am", (0, 7), 7),
            Element::Cm => (96, 247.0, "Cm", (0, 6), 8),
            Element::Bk => (97, 247.0, "Bk", (0, 5), 5),
            Element::Cf => (98, 251.0, "Cf", (0, 5), 4),
            Element::Es => (99, 252.0, "Es", (0, 4), 3),
            Element::Fm => (100, 257.0, "Fm", (0, 3), 2),
            Element::Md => (101, 258.0, "Md", (0, 3), 1),
            Element::No => (102, 259.0, "No", (0, 3), 0),
            Element::Lr => (103, 262.0, "Lr", (0, 3), 1),
            Element::Rf => (104, 267.0, "Rf", (0, 4), 2),
            Element::Db => (105, 270.0, "Db", (0, 5), 3),
            Element::Sg => (106, 271.0, "Sg", (0, 6), 4),
            Element::Bh => (107, 270.0, "Bh", (0, 7), 5),
            Element::Hs => (108, 277.0, "Hs", (0, 8), 6),
            Element::Mt => (109, 276.0, "Mt", (0, 6), 5),
            Element::Ds => (110, 281.0, "Ds", (0, 6), 4),
            Element::Rg => (111, 280.0, "Rg", (0, 5), 3),
            Element::Cn => (112, 285.0, "Cn", (0, 4), 2),
            Element::Nh => (113, 284.0, "Nh", (0, 0), 0),
            Element::Fl => (114, 289.0, "Fl", (0, 0), 0),
            Element::Mc => (115, 288.0, "Mc", (0, 0), 0),
            Element::Lv => (116, 293.0, "Lv", (0, 0), 0),
            Element::Ts => (117, 294.0, "Ts", (0, 0), 0),
            Element::Og => (118, 294.0, "Og", (0, 0), 0),
        }
    });

static SYMBOL_TO_ELEMENT: Lazy<HashMap<&'static str, Element>> = Lazy::new(|| {
    ELEMENT_DATA
        .iter()
        .map(|(element, (_, _, symbol, _, _))| (*symbol, *element))
        .collect()
});

static ATOMIC_NUMBER_TO_ELEMENT: Lazy<HashMap<u8, Element>> = Lazy::new(|| {
    ELEMENT_DATA
        .iter()
        .map(|(element, (number, _, _, _, _))| (*number, *element))
        .collect()
});

impl Element {
    pub fn from_symbol(symbol: &str) -> Option<Self> {
        SYMBOL_TO_ELEMENT.get(symbol).copied()
    }

    pub fn from_atomic_number(number: u8) -> Option<Self> {
        ATOMIC_NUMBER_TO_ELEMENT.get(&number).copied()
    }

    pub fn symbol(&self) -> &'static str {
        ELEMENT_DATA.get(self).unwrap().2
    }

    pub fn atomic_number(&self) -> u8 {
        ELEMENT_DATA.get(self).unwrap().0
    }

    pub fn atomic_mass(&self) -> f64 {
        ELEMENT_DATA.get(self).unwrap().1
    }

    pub fn validate_charge(&self, charge: i8) -> Result<()> {
        let (min_charge, max_charge) = ELEMENT_DATA.get(self).unwrap().3;

        if charge < min_charge || charge > max_charge {
            return Err(Error::Data(DataError::InvalidCharge {
                element: *self,
                charge,
                min: min_charge,
                max: max_charge,
            }));
        }
        Ok(())
    }

    pub fn validate_unpaired_electrons(&self, unpaired: u8) -> Result<()> {
        let max_unpaired = ELEMENT_DATA.get(self).unwrap().4;

        if unpaired > max_unpaired {
            return Err(Error::Data(DataError::InvalidUnpairedElectrons {
                element: *self,
                unpaired,
                max: max_unpaired,
            }));
        }
        Ok(())
    }

    pub fn charge_bounds(&self) -> (i8, i8) {
        ELEMENT_DATA.get(self).unwrap().3
    }

    pub fn max_unpaired_electrons(&self) -> u8 {
        ELEMENT_DATA.get(self).unwrap().4
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
    fn test_element_properties() {
        let h = Element::H;
        assert_eq!(h.symbol(), "H");
        assert_eq!(h.atomic_number(), 1);
        assert!((h.atomic_mass() - 1.008).abs() < 1e-10);

        let he = Element::He;
        assert_eq!(he.symbol(), "He");
        assert_eq!(he.atomic_number(), 2);
        assert!((he.atomic_mass() - 4.002602).abs() < 1e-10);
    }

    #[rstest]
    #[case(Element::H, -1, 1, 1)]
    #[case(Element::He, 0, 0, 0)]
    #[case(Element::Li, -1, 1, 1)]
    #[case(Element::Be, -2, 2, 2)]
    fn test_element_bounds(
        #[case] element: Element,
        #[case] min_charge: i8,
        #[case] max_charge: i8,
        #[case] max_unpaired: u8,
    ) {
        let (actual_min, actual_max) = element.charge_bounds();
        assert_eq!(actual_min, min_charge);
        assert_eq!(actual_max, max_charge);
        assert_eq!(element.max_unpaired_electrons(), max_unpaired);
    }

    #[rstest]
    #[case(Element::H, 0, true)]
    #[case(Element::H, 1, true)]
    #[case(Element::H, -1, true)]
    #[case(Element::H, 2, false)]
    #[case(Element::H, -2, false)]
    fn test_charge_validation(
        #[case] element: Element,
        #[case] charge: i8,
        #[case] should_be_valid: bool,
    ) {
        let result = element.validate_charge(charge);
        assert_eq!(result.is_ok(), should_be_valid);
    }

    #[rstest]
    #[case(Element::H, 0, true)]
    #[case(Element::H, 1, true)]
    #[case(Element::H, 2, false)]
    #[case(Element::He, 0, true)]
    #[case(Element::He, 1, false)]
    fn test_unpaired_electrons_validation(
        #[case] element: Element,
        #[case] unpaired: u8,
        #[case] should_be_valid: bool,
    ) {
        let result = element.validate_unpaired_electrons(unpaired);
        assert_eq!(result.is_ok(), should_be_valid);
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
            assert!(data.3 .0 <= data.3 .1); // charge bounds
            assert!(data.4 <= 10); // max unpaired electrons
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
        assert_eq!(e!(C), Element::C);
        assert_eq!(e!(O), Element::O);
        assert_eq!(e!(Fe), Element::Fe);
    }
}
