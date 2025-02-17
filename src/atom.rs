// Core AtomSite trait and Element enum

use once_cell::sync::Lazy;
use std::collections::HashMap;

pub trait AtomSite: Sized {
    fn element(&self) -> Option<Element>;
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Element {
    H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca,
    Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr,
    Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba,
    La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W,
    Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th, Pa, U,
    Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Ds,
    Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og,
}

static ELEMENT_DATA: Lazy<HashMap<Element, (u8, f64, &'static str)>> = Lazy::new(|| {
    let mut m = HashMap::new();
    m.insert(Element::H, (1, 1.008, "H"));
    m.insert(Element::He, (2, 4.002602, "He"));
    m.insert(Element::Li, (3, 6.94, "Li"));
    m.insert(Element::Be, (4, 9.0121831, "Be"));
    m.insert(Element::B, (5, 10.81, "B"));
    m.insert(Element::C, (6, 12.011, "C"));
    m.insert(Element::N, (7, 14.007, "N"));
    m.insert(Element::O, (8, 15.999, "O"));
    m.insert(Element::F, (9, 18.998403163, "F"));
    m.insert(Element::Ne, (10, 20.1797, "Ne"));
    m.insert(Element::Na, (11, 22.98976928, "Na"));
    m.insert(Element::Mg, (12, 24.305, "Mg"));
    m.insert(Element::Al, (13, 26.9815385, "Al"));
    m.insert(Element::Si, (14, 28.085, "Si"));
    m.insert(Element::P, (15, 30.973761998, "P"));
    m.insert(Element::S, (16, 32.06, "S"));
    m.insert(Element::Cl, (17, 35.45, "Cl"));
    m.insert(Element::Ar, (18, 39.948, "Ar"));
    m.insert(Element::K, (19, 39.0983, "K"));
    m.insert(Element::Ca, (20, 40.078, "Ca"));
    m.insert(Element::Sc, (21, 44.955908, "Sc"));
    m.insert(Element::Ti, (22, 47.867, "Ti"));
    m.insert(Element::V, (23, 50.9415, "Va"));
    m.insert(Element::Cr, (24, 51.9961, "Cr"));
    m.insert(Element::Mn, (25, 54.938044, "Mn"));
    m.insert(Element::Fe, (26, 55.845, "Fe"));
    m.insert(Element::Co, (27, 58.933194, "Co"));
    m.insert(Element::Ni, (28, 58.6934, "Ni"));
    m.insert(Element::Cu, (29, 63.546, "Cu"));
    m.insert(Element::Zn, (30, 65.38, "Zn"));
    m.insert(Element::Ga, (31, 69.723, "Ga"));
    m.insert(Element::Ge, (32, 72.63, "Ge"));
    m.insert(Element::As, (33, 74.921595, "As"));
    m.insert(Element::Se, (34, 78.971, "Se"));
    m.insert(Element::Br, (35, 79.904, "Br"));
    m.insert(Element::Kr, (36, 83.798, "Kr"));
    m.insert(Element::Rb, (37, 85.4678, "Rb"));
    m.insert(Element::Sr, (38, 87.62, "Sr"));
    m.insert(Element::Y, (39, 88.90584, "Y"));
    m.insert(Element::Zr, (40, 91.224, "Zr"));
    m.insert(Element::Nb, (41, 92.90637, "Nb"));
    m.insert(Element::Mo, (42, 95.95, "Mo"));
    m.insert(Element::Tc, (43, 98.0, "Tc"));
    m.insert(Element::Ru, (44, 101.07, "Ru"));
    m.insert(Element::Rh, (45, 102.90550, "Rh"));
    m.insert(Element::Pd, (46, 106.42, "Pd"));
    m.insert(Element::Ag, (47, 107.8682, "Ag"));
    m.insert(Element::Cd, (48, 112.414, "Cd"));
    m.insert(Element::In, (49, 114.818, "In"));
    m.insert(Element::Sn, (50, 118.710, "Sn"));
    m.insert(Element::Sb, (51, 121.760, "Sb"));
    m.insert(Element::Te, (52, 127.60, "Te"));
    m.insert(Element::I, (53, 126.90447, "I"));
    m.insert(Element::Xe, (54, 131.293, "Xe"));
    m.insert(Element::Cs, (55, 132.90545196, "Cs"));
    m.insert(Element::Ba, (56, 137.327, "Ba"));
    m.insert(Element::La, (57, 138.90547, "La"));
    m.insert(Element::Ce, (58, 140.116, "Ce"));
    m.insert(Element::Pr, (59, 140.90766, "Pr"));
    m.insert(Element::Nd, (60, 144.242, "Nd"));
    m.insert(Element::Pm, (61, 145.0, "Pm"));
    m.insert(Element::Sm, (62, 150.36, "Sm"));
    m.insert(Element::Eu, (63, 151.964, "Eu"));
    m.insert(Element::Gd, (64, 157.25, "Gd"));
    m.insert(Element::Tb, (65, 158.92535, "Tb"));
    m.insert(Element::Dy, (66, 162.500, "Dy"));
    m.insert(Element::Ho, (67, 164.93033, "Ho"));
    m.insert(Element::Er, (68, 167.259, "Er"));
    m.insert(Element::Tm, (69, 168.93422, "Tm"));
    m.insert(Element::Yb, (70, 173.045, "Yb"));
    m.insert(Element::Lu, (71, 174.9668, "Lu"));
    m.insert(Element::Hf, (72, 178.49, "Hf"));
    m.insert(Element::Ta, (73, 180.94788, "Ta"));
    m.insert(Element::W, (74, 183.84, "W"));
    m.insert(Element::Re, (75, 186.207, "Re"));
    m.insert(Element::Os, (76, 190.23, "Os"));
    m.insert(Element::Ir, (77, 192.217, "Ir"));
    m.insert(Element::Pt, (78, 195.084, "Pt"));
    m.insert(Element::Au, (79, 196.966569, "Au"));
    m.insert(Element::Hg, (80, 200.592, "Hg"));
    m.insert(Element::Tl, (81, 204.38, "Tl"));
    m.insert(Element::Pb, (82, 207.2, "Pb"));
    m.insert(Element::Bi, (83, 208.98040, "Bi"));
    m.insert(Element::Po, (84, 209.0, "Po"));
    m.insert(Element::At, (85, 210.0, "At"));
    m.insert(Element::Rn, (86, 222.0, "Rn"));
    m.insert(Element::Fr, (87, 223.0, "Fr"));
    m.insert(Element::Ra, (88, 226.0, "Ra"));
    m.insert(Element::Ac, (89, 227.0, "Ac"));
    m.insert(Element::Th, (90, 232.0377, "Th"));
    m.insert(Element::Pa, (91, 231.03588, "Pa"));
    m.insert(Element::U, (92, 238.02891, "U"));
    m.insert(Element::Np, (93, 237.0, "Np"));
    m.insert(Element::Pu, (94, 244.0, "Pu"));
    m.insert(Element::Am, (95, 243.0, "Am"));
    m.insert(Element::Cm, (96, 247.0, "Cm"));
    m.insert(Element::Bk, (97, 247.0, "Bk"));
    m.insert(Element::Cf, (98, 251.0, "Cf"));
    m.insert(Element::Es, (99, 252.0, "Es"));
    m.insert(Element::Fm, (100, 257.0, "Fm"));
    m.insert(Element::Md, (101, 258.0, "Md"));
    m.insert(Element::No, (102, 259.0, "No"));
    m.insert(Element::Lr, (103, 262.0, "Lr"));
    m.insert(Element::Rf, (104, 267.0, "Rf"));
    m.insert(Element::Db, (105, 270.0, "Db"));
    m.insert(Element::Sg, (106, 271.0, "Sg"));
    m.insert(Element::Bh, (107, 270.0, "Bh"));
    m.insert(Element::Hs, (108, 277.0, "Hs"));
    m.insert(Element::Mt, (109, 276.0, "Mt"));
    m.insert(Element::Ds, (110, 281.0, "Ds"));
    m.insert(Element::Rg, (111, 280.0, "Rg"));
    m.insert(Element::Cn, (112, 285.0, "Cn"));
    m.insert(Element::Nh, (113, 284.0, "Nh"));
    m.insert(Element::Fl, (114, 289.0, "Fl"));
    m.insert(Element::Mc, (115, 288.0, "Mc"));
    m.insert(Element::Lv, (116, 293.0, "Lv"));
    m.insert(Element::Ts, (117, 294.0, "Ts"));
    m.insert(Element::Og, (118, 294.0, "Og"));
    m
});

static SYMBOL_TO_ELEMENT: Lazy<HashMap<&'static str, Element>> = Lazy::new(|| {
    ELEMENT_DATA.iter().map(|(element, (_, _, symbol))| (*symbol, *element)).collect()
});


static ATOMIC_NUMBER_TO_ELEMENT: Lazy<HashMap<u8, Element>> = Lazy::new(|| {
    ELEMENT_DATA
        .iter()
        .map(|(element, (number, _, _))| (*number, *element))
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
}