//! Element definitions and data

use serde::{Deserialize, Serialize};
use std::fmt::{self, Display};
use std::str::FromStr;
use umol::error::DataError;
use umol::{Error, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
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

// Element data array indexed by atomic number - 1
// Each tuple contains: (atomic_number, reference_isotope_mass, atomic_weight, symbol, period, group,
// valence_electrons, max_valence, (min_charge, max_charge), max_unpaired_electrons, max_implicit_hydrogens)
#[allow(clippy::type_complexity)]
static ELEMENT_DATA: [(u8, u32, f64, &str, u8, u8, u8, u8, (i8, i8), u8, u8); 118] = [
    (1, 1, 1.0080, "H", 1, 1, 1, 2, (-1, 1), 1, 0),    // H
    (2, 4, 4.0026, "He", 1, 32, 0, 0, (0, 0), 0, 0),   // He: no valence electrons
    (3, 7, 6.94, "Li", 2, 1, 1, 2, (-1, 1), 1, 0),     // Li
    (4, 9, 9.0122, "Be", 2, 2, 2, 2, (-2, 2), 2, 0),   // Be
    (5, 11, 10.81, "B", 2, 3, 27, 8, (-3, 3), 3, 3),   // B
    (6, 12, 12.011, "C", 2, 28, 4, 8, (-4, 4), 4, 4),  // C
    (7, 14, 14.007, "N", 2, 29, 5, 8, (-3, 5), 3, 3),  // N
    (8, 16, 15.999, "O", 2, 30, 6, 8, (-2, 6), 2, 2),  // O
    (9, 19, 18.998, "F", 2, 31, 7, 8, (-1, 7), 1, 1),  // F
    (10, 20, 20.180, "Ne", 2, 32, 0, 0, (0, 0), 0, 0), // Ne: no valence electrons
    (11, 23, 22.990, "Na", 3, 1, 1, 2, (-1, 1), 1, 0), // Na
    (12, 24, 24.305, "Mg", 3, 2, 2, 2, (0, 2), 2, 0),  // Mg
    (13, 27, 26.982, "Al", 3, 27, 3, 8, (-1, 3), 3, 0), // Al
    (14, 28, 28.085, "Si", 3, 28, 4, 8, (-4, 4), 4, 4), // Si
    (15, 31, 30.974, "P", 3, 29, 5, 8, (-3, 5), 3, 3), // P
    (16, 32, 32.06, "S", 3, 30, 6, 8, (-2, 6), 2, 2),  // S
    (17, 35, 35.45, "Cl", 3, 31, 7, 8, (-1, 7), 1, 1), // Cl
    (18, 40, 39.95, "Ar", 3, 32, 0, 0, (0, 0), 0, 0),  // Ar: no valence electrons
    (19, 39, 39.098, "K", 4, 1, 1, 2, (-1, 1), 1, 0),  // K
    (20, 40, 40.078, "Ca", 4, 2, 2, 2, (0, 2), 2, 0),  // Ca
    (21, 45, 44.956, "Sc", 4, 3, 3, 12, (0, 3), 3, 0), // Sc
    (22, 48, 47.867, "Ti", 4, 18, 4, 18, (0, 4), 4, 0), // Ti
    (23, 51, 50.942, "V", 4, 19, 5, 18, (0, 5), 5, 0), // V
    (24, 52, 51.996, "Cr", 4, 20, 6, 18, (0, 6), 6, 0), // Cr
    (25, 55, 54.938, "Mn", 4, 21, 7, 18, (0, 7), 7, 0), // Mn
    (26, 56, 55.845, "Fe", 4, 22, 8, 18, (0, 6), 6, 0), // Fe
    (27, 59, 58.933, "Co", 4, 23, 9, 18, (0, 5), 5, 0), // Co
    (28, 58, 58.693, "Ni", 4, 24, 10, 18, (0, 4), 4, 0), // Ni
    (29, 63, 63.546, "Cu", 4, 25, 11, 18, (0, 3), 3, 0), // Cu
    (30, 64, 65.38, "Zn", 4, 26, 12, 18, (0, 2), 2, 0), // Zn
    (31, 69, 69.723, "Ga", 4, 27, 13, 18, (0, 3), 3, 0), // Ga
    (32, 74, 72.630, "Ge", 4, 28, 14, 18, (0, 4), 4, 0), // Ge
    (33, 75, 74.922, "As", 4, 29, 15, 18, (0, 3), 3, 3), // As
    (34, 80, 78.971, "Se", 4, 30, 16, 18, (0, 2), 2, 2), // Se
    (35, 79, 79.904, "Br", 4, 31, 17, 18, (0, 1), 1, 1), // Br
    (36, 84, 83.798, "Kr", 4, 32, 18, 18, (0, 0), 0, 0), // Kr: allow Kr compounds
    (37, 85, 85.468, "Rb", 5, 1, 1, 2, (-1, 1), 1, 0), // Rb
    (38, 88, 87.62, "Sr", 5, 2, 2, 2, (0, 2), 2, 0),   // Sr
    (39, 89, 88.906, "Y", 5, 3, 3, 12, (0, 3), 3, 0),  // Y
    (40, 90, 91.224, "Zr", 5, 18, 4, 18, (0, 4), 4, 0), // Zr
    (41, 93, 92.906, "Nb", 5, 19, 5, 18, (0, 5), 5, 0), // Nb
    (42, 98, 95.95, "Mo", 5, 20, 6, 18, (0, 6), 6, 0), // Mo
    (43, 97, 98.0, "Tc", 5, 21, 7, 18, (0, 7), 7, 0), // Tc: using 97Tc as reference isotope based on PubChem data
    (44, 102, 101.07, "Ru", 5, 22, 8, 18, (0, 8), 8, 0), // Ru
    (45, 103, 102.91, "Rh", 5, 23, 9, 18, (0, 6), 6, 0), // Rh
    (46, 106, 106.42, "Pd", 5, 24, 10, 18, (0, 5), 5, 0), // Pd
    (47, 107, 107.87, "Ag", 5, 25, 11, 18, (0, 3), 3, 0), // Ag
    (48, 114, 112.41, "Cd", 5, 26, 12, 18, (0, 2), 2, 0), // Cd
    (49, 115, 114.82, "In", 5, 27, 13, 18, (0, 3), 3, 0), // In
    (50, 120, 118.71, "Sn", 5, 28, 14, 18, (0, 4), 4, 4), // Sn
    (51, 121, 121.76, "Sb", 5, 29, 15, 18, (0, 3), 3, 3), // Sb
    (52, 130, 127.60, "Te", 5, 30, 16, 18, (0, 2), 2, 2), // Te
    (53, 127, 126.90, "I", 5, 31, 17, 18, (0, 1), 1, 1), // I
    (54, 132, 131.29, "Xe", 5, 32, 18, 18, (0, 0), 0, 0), // Xe: allow Xe compounds
    (55, 133, 132.91, "Cs", 6, 1, 1, 2, (-1, 1), 1, 0), // Cs
    (56, 138, 137.33, "Ba", 6, 2, 2, 2, (0, 2), 2, 0), // Ba
    (57, 139, 138.91, "La", 6, 3, 3, 18, (0, 3), 3, 0), // La
    (58, 140, 140.12, "Ce", 6, 4, 4, 20, (0, 4), 4, 0), // Ce
    (59, 141, 140.91, "Pr", 6, 5, 5, 20, (0, 4), 4, 0), // Pr
    (60, 142, 144.24, "Nd", 6, 6, 6, 22, (-1, 3), 3, 0), // Nd
    (61, 145, 145.0, "Pm", 6, 7, 7, 22, (-1, 3), 3, 0), // Pm
    (62, 152, 150.36, "Sm", 6, 8, 8, 24, (-1, 3), 3, 0), // Sm
    (63, 153, 151.96, "Eu", 6, 9, 9, 24, (-1, 3), 3, 0), // Eu
    (64, 158, 157.25, "Gd", 6, 10, 10, 26, (-1, 3), 3, 0), // Gd
    (65, 159, 158.93, "Tb", 6, 11, 11, 26, (-1, 4), 4, 0), // Tb
    (66, 164, 162.50, "Dy", 6, 12, 12, 28, (-1, 3), 3, 0), // Dy
    (67, 165, 164.93, "Ho", 6, 13, 13, 28, (-1, 3), 3, 0), // Ho
    (68, 166, 167.26, "Er", 6, 14, 14, 30, (-1, 3), 3, 0), // Er
    (69, 169, 168.93, "Tm", 6, 15, 15, 30, (-1, 3), 3, 0), // Tm
    (70, 174, 173.05, "Yb", 6, 16, 16, 32, (-1, 3), 3, 0), // Yb
    (71, 175, 174.97, "Lu", 6, 17, 17, 32, (-1, 3), 3, 0), // Lu
    (72, 180, 178.49, "Hf", 6, 18, 18, 32, (0, 4), 4, 0), // Hf
    (73, 181, 180.95, "Ta", 6, 19, 19, 32, (0, 5), 5, 0), // Ta
    (74, 184, 183.84, "W", 6, 20, 20, 32, (0, 6), 6, 0), // W
    (75, 187, 186.21, "Re", 6, 21, 21, 32, (0, 7), 7, 0), // Re
    (76, 192, 190.23, "Os", 6, 22, 22, 32, (0, 8), 8, 0), // Os
    (77, 193, 192.22, "Ir", 6, 23, 23, 32, (0, 6), 6, 0), // Ir
    (78, 195, 195.08, "Pt", 6, 24, 24, 32, (0, 6), 4, 0), // Pt
    (79, 197, 196.97, "Au", 6, 25, 25, 32, (0, 5), 1, 0), // Au
    (80, 202, 200.59, "Hg", 6, 26, 26, 32, (0, 2), 2, 0), // Hg
    (81, 205, 204.38, "Tl", 6, 27, 27, 32, (0, 3), 3, 0), // Tl
    (82, 208, 207.2, "Pb", 6, 28, 28, 32, (0, 4), 4, 0), // Pb
    (83, 209, 208.98, "Bi", 6, 29, 29, 32, (0, 3), 3, 3), // Bi
    (84, 209, 209.0, "Po", 6, 30, 30, 32, (0, 2), 2, 2), // Po
    (85, 210, 210.0, "At", 6, 31, 31, 32, (0, 7), 1, 1), // At
    (86, 222, 222.0, "Rn", 6, 32, 32, 32, (0, 8), 2, 0), // Rn
    (87, 223, 223.0, "Fr", 7, 1, 1, 2, (0, 1), 1, 0), // Fr
    (88, 226, 226.0, "Ra", 7, 2, 2, 2, (0, 2), 2, 0), // Ra
    (89, 227, 227.0, "Ac", 7, 3, 3, 18, (0, 3), 3, 0), // Ac
    (90, 232, 232.04, "Th", 7, 4, 4, 20, (0, 4), 4, 0), // Th
    (91, 231, 231.04, "Pa", 7, 5, 5, 20, (0, 5), 3, 0), // Pa
    (92, 238, 238.03, "U", 7, 6, 6, 22, (0, 6), 4, 0), // U
    (93, 237, 237.0, "Np", 7, 7, 7, 22, (0, 7), 5, 0), // Np
    (94, 244, 244.0, "Pu", 7, 8, 8, 24, (0, 8), 6, 0), // Pu
    (95, 243, 243.0, "Am", 7, 9, 9, 24, (0, 7), 7, 0), // Am
    (96, 247, 247.0, "Cm", 7, 10, 10, 26, (0, 6), 8, 0), // Cm
    (97, 247, 247.0, "Bk", 7, 11, 11, 26, (0, 5), 5, 0), // Bk
    (98, 251, 251.0, "Cf", 7, 12, 12, 28, (0, 5), 4, 0), // Cf
    (99, 252, 252.0, "Es", 7, 13, 13, 28, (0, 4), 3, 0), // Es
    (100, 257, 257.0, "Fm", 7, 14, 14, 30, (0, 3), 2, 0), // Fm
    (101, 258, 258.0, "Md", 7, 15, 15, 30, (0, 3), 1, 0), // Md
    (102, 259, 259.0, "No", 7, 16, 16, 32, (0, 3), 0, 0), // No
    (103, 266, 226.0, "Lr", 7, 17, 17, 32, (0, 3), 1, 0), // Lr
    (104, 267, 267.0, "Rf", 7, 18, 18, 32, (0, 4), 2, 0), // Rf
    (105, 268, 268.0, "Db", 7, 19, 19, 32, (0, 5), 3, 0), // Db
    (106, 269, 269.0, "Sg", 7, 20, 20, 32, (0, 6), 4, 0), // Sg
    (107, 270, 270.0, "Bh", 7, 21, 21, 32, (0, 7), 5, 0), // Bh
    (108, 269, 269.0, "Hs", 7, 22, 22, 32, (0, 8), 6, 0), // Hs
    (109, 277, 277.0, "Mt", 7, 23, 23, 32, (0, 6), 5, 0), // Mt
    (110, 282, 282.0, "Ds", 7, 24, 24, 32, (0, 6), 4, 0), // Ds
    (111, 282, 282.0, "Rg", 7, 25, 25, 32, (0, 5), 3, 0), // Rg
    (112, 286, 286.0, "Cn", 7, 26, 26, 32, (0, 4), 2, 0), // Cn
    (113, 286, 286.0, "Nh", 7, 27, 27, 32, (0, 3), 1, 0), // Nh
    (114, 290, 290.0, "Fl", 7, 28, 28, 32, (0, 2), 0, 0), // Fl
    (115, 290, 290.0, "Mc", 7, 29, 29, 32, (0, 0), 0, 0), // Mc
    (116, 293, 293.0, "Lv", 7, 30, 30, 32, (0, 0), 0, 0), // Lv
    (117, 294, 294.0, "Ts", 7, 31, 31, 32, (0, 0), 0, 0), // Ts
    (118, 294, 294.0, "Og", 7, 32, 32, 32, (0, 0), 0, 0), // Og: PubChem uses 295 although only 294Og is known
];

/// Maximum atomic number, group number, and period number  
pub const MAX_ATOMIC_NUMBER: u8 = 118;
pub const MAX_PERIOD_NUMBER: u8 = 7;
pub const MIN_PERIOD_NUMBER: [u8; 33] = [
    0, 1, 2, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 2, 2, 2, 2,
    1,
];
pub const MAX_GROUP_NUMBER: [u8; 8] = [0, 2, 8, 8, 18, 18, 32, 32];

/// Element array indexed by atomic number - 1
#[rustfmt::skip]
#[allow(dead_code)]
static ELEMENTS: [Element; 118] = [
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

/// Lookup table for element symbols using direct indexing
/// Index calculation:
/// - 1-character elements: index = (first_char - 'A')
/// - 2-character elements: index = 26 + (first_char - 'A') * 26 + (second_char - 'a')
static ELEMENT_LOOKUP: [bool; 702] = {
    let mut table = [false; 702];

    // 1-character elements (indices 0-25)
    table[1] = true; // B
    table[2] = true; // C
    table[5] = true; // F
    table[7] = true; // H
    table[8] = true; // I
    table[10] = true; // K
    table[13] = true; // N
    table[14] = true; // O
    table[15] = true; // P
    table[18] = true; // S
    table[20] = true; // U
    table[21] = true; // V
    table[22] = true; // W
    table[24] = true; // Y

    // 2-character elements (indices 26-701)
    table[28] = true; // Ac
    table[32] = true; // Ag
    table[37] = true; // Al
    table[38] = true; // Am
    table[43] = true; // Ar
    table[44] = true; // As
    table[45] = true; // At
    table[46] = true; // Au
    table[52] = true; // Ba
    table[56] = true; // Be
    table[59] = true; // Bh
    table[60] = true; // Bi
    table[62] = true; // Bk
    table[69] = true; // Br
    table[78] = true; // Ca
    table[81] = true; // Cd
    table[82] = true; // Ce
    table[83] = true; // Cf
    table[89] = true; // Cl
    table[90] = true; // Cm
    table[91] = true; // Cn
    table[92] = true; // Co
    table[95] = true; // Cr
    table[96] = true; // Cs
    table[98] = true; // Cu
    table[105] = true; // Db
    table[124] = true; // Ds
    table[128] = true; // Dy
    table[147] = true; // Er
    table[148] = true; // Es
    table[150] = true; // Eu
    table[160] = true; // Fe
    table[167] = true; // Fl
    table[168] = true; // Fm
    table[173] = true; // Fr
    table[182] = true; // Ga
    table[185] = true; // Gd
    table[186] = true; // Ge
    table[212] = true; // He
    table[213] = true; // Hf
    table[214] = true; // Hg
    table[222] = true; // Ho
    table[226] = true; // Hs
    table[247] = true; // In
    table[251] = true; // Ir
    table[303] = true; // Kr
    table[312] = true; // La
    table[320] = true; // Li
    table[329] = true; // Lr
    table[332] = true; // Lu
    table[333] = true; // Lv
    table[340] = true; // Mc
    table[341] = true; // Md
    table[344] = true; // Mg
    table[351] = true; // Mn
    table[352] = true; // Mo
    table[357] = true; // Mt
    table[364] = true; // Na
    table[365] = true; // Nb
    table[367] = true; // Nd
    table[368] = true; // Ne
    table[371] = true; // Nh
    table[372] = true; // Ni
    table[378] = true; // No
    table[379] = true; // Np
    table[396] = true; // Og
    table[408] = true; // Os
    table[416] = true; // Pa
    table[417] = true; // Pb
    table[419] = true; // Pd
    table[428] = true; // Pm
    table[430] = true; // Po
    table[433] = true; // Pr
    table[435] = true; // Pt
    table[436] = true; // Pu
    table[468] = true; // Ra
    table[469] = true; // Rb
    table[472] = true; // Re
    table[473] = true; // Rf
    table[474] = true; // Rg
    table[475] = true; // Rh
    table[481] = true; // Rn
    table[488] = true; // Ru
    table[495] = true; // Sb
    table[496] = true; // Sc
    table[498] = true; // Se
    table[500] = true; // Sg
    table[502] = true; // Si
    table[506] = true; // Sm
    table[507] = true; // Sn
    table[511] = true; // Sr
    table[520] = true; // Ta
    table[521] = true; // Tb
    table[522] = true; // Tc
    table[524] = true; // Te
    table[527] = true; // Th
    table[528] = true; // Ti
    table[531] = true; // Tl
    table[532] = true; // Tm
    table[538] = true; // Ts
    table[628] = true; // Xe
    table[651] = true; // Yb
    table[663] = true; // Zn
    table[667] = true; // Zr

    table
};

/// Last element
pub const LAST_ELEMENT: Element = Element::Og;

/// Period/Group lookup table using direct 2D array indexing
/// [period-1][group-1] -> Option<Element>
static PERIOD_GROUP_TO_ELEMENT: [[Option<Element>; 33]; 8] = {
    let mut table = [[None; 33]; 8];

    // Period 1
    table[0][0] = Some(Element::H); // (1,1)
    table[0][31] = Some(Element::He); // (1,32)

    // Period 2
    table[1][0] = Some(Element::Li); // (2,1)
    table[1][1] = Some(Element::Be); // (2,2)
    table[1][2] = Some(Element::B); // (2,3)
    table[1][27] = Some(Element::C); // (2,28)
    table[1][28] = Some(Element::N); // (2,29)
    table[1][29] = Some(Element::O); // (2,30)
    table[1][30] = Some(Element::F); // (2,31)
    table[1][31] = Some(Element::Ne); // (2,32)

    // Period 3
    table[2][0] = Some(Element::Na); // (3,1)
    table[2][1] = Some(Element::Mg); // (3,2)
    table[2][26] = Some(Element::Al); // (3,27)
    table[2][27] = Some(Element::Si); // (3,28)
    table[2][28] = Some(Element::P); // (3,29)
    table[2][29] = Some(Element::S); // (3,30)
    table[2][30] = Some(Element::Cl); // (3,31)
    table[2][31] = Some(Element::Ar); // (3,32)

    // Period 4
    table[3][0] = Some(Element::K); // (4,1)
    table[3][1] = Some(Element::Ca); // (4,2)
    table[3][2] = Some(Element::Sc); // (4,3)
    table[3][17] = Some(Element::Ti); // (4,18)
    table[3][18] = Some(Element::V); // (4,19)
    table[3][19] = Some(Element::Cr); // (4,20)
    table[3][20] = Some(Element::Mn); // (4,21)
    table[3][21] = Some(Element::Fe); // (4,22)
    table[3][22] = Some(Element::Co); // (4,23)
    table[3][23] = Some(Element::Ni); // (4,24)
    table[3][24] = Some(Element::Cu); // (4,25)
    table[3][25] = Some(Element::Zn); // (4,26)
    table[3][26] = Some(Element::Ga); // (4,27)
    table[3][27] = Some(Element::Ge); // (4,28)
    table[3][28] = Some(Element::As); // (4,29)
    table[3][29] = Some(Element::Se); // (4,30)
    table[3][30] = Some(Element::Br); // (4,31)
    table[3][31] = Some(Element::Kr); // (4,32)

    // Period 5
    table[4][0] = Some(Element::Rb); // (5,1)
    table[4][1] = Some(Element::Sr); // (5,2)
    table[4][2] = Some(Element::Y); // (5,3)
    table[4][17] = Some(Element::Zr); // (5,18)
    table[4][18] = Some(Element::Nb); // (5,19)
    table[4][19] = Some(Element::Mo); // (5,20)
    table[4][20] = Some(Element::Tc); // (5,21)
    table[4][21] = Some(Element::Ru); // (5,22)
    table[4][22] = Some(Element::Rh); // (5,23)
    table[4][23] = Some(Element::Pd); // (5,24)
    table[4][24] = Some(Element::Ag); // (5,25)
    table[4][25] = Some(Element::Cd); // (5,26)
    table[4][26] = Some(Element::In); // (5,27)
    table[4][27] = Some(Element::Sn); // (5,28)
    table[4][28] = Some(Element::Sb); // (5,29)
    table[4][29] = Some(Element::Te); // (5,30)
    table[4][30] = Some(Element::I); // (5,31)
    table[4][31] = Some(Element::Xe); // (5,32)

    // Period 6
    table[5][0] = Some(Element::Cs); // (6,1)
    table[5][1] = Some(Element::Ba); // (6,2)
    table[5][2] = Some(Element::La); // (6,3)
    table[5][3] = Some(Element::Ce); // (6,4)
    table[5][4] = Some(Element::Pr); // (6,5)
    table[5][5] = Some(Element::Nd); // (6,6)
    table[5][6] = Some(Element::Pm); // (6,7)
    table[5][7] = Some(Element::Sm); // (6,8)
    table[5][8] = Some(Element::Eu); // (6,9)
    table[5][9] = Some(Element::Gd); // (6,10)
    table[5][10] = Some(Element::Tb); // (6,11)
    table[5][11] = Some(Element::Dy); // (6,12)
    table[5][12] = Some(Element::Ho); // (6,13)
    table[5][13] = Some(Element::Er); // (6,14)
    table[5][14] = Some(Element::Tm); // (6,15)
    table[5][15] = Some(Element::Yb); // (6,16)
    table[5][16] = Some(Element::Lu); // (6,17)
    table[5][17] = Some(Element::Hf); // (6,18)
    table[5][18] = Some(Element::Ta); // (6,19)
    table[5][19] = Some(Element::W); // (6,20)
    table[5][20] = Some(Element::Re); // (6,21)
    table[5][21] = Some(Element::Os); // (6,22)
    table[5][22] = Some(Element::Ir); // (6,23)
    table[5][23] = Some(Element::Pt); // (6,24)
    table[5][24] = Some(Element::Au); // (6,25)
    table[5][25] = Some(Element::Hg); // (6,26)
    table[5][26] = Some(Element::Tl); // (6,27)
    table[5][27] = Some(Element::Pb); // (6,28)
    table[5][28] = Some(Element::Bi); // (6,29)
    table[5][29] = Some(Element::Po); // (6,30)
    table[5][30] = Some(Element::At); // (6,31)
    table[5][31] = Some(Element::Rn); // (6,32)

    // Period 7
    table[6][0] = Some(Element::Fr); // (7,1)
    table[6][1] = Some(Element::Ra); // (7,2)
    table[6][2] = Some(Element::Ac); // (7,3)
    table[6][3] = Some(Element::Th); // (7,4)
    table[6][4] = Some(Element::Pa); // (7,5)
    table[6][5] = Some(Element::U); // (7,6)
    table[6][6] = Some(Element::Np); // (7,7)
    table[6][7] = Some(Element::Pu); // (7,8)
    table[6][8] = Some(Element::Am); // (7,9)
    table[6][9] = Some(Element::Cm); // (7,10)
    table[6][10] = Some(Element::Bk); // (7,11)
    table[6][11] = Some(Element::Cf); // (7,12)
    table[6][12] = Some(Element::Es); // (7,13)
    table[6][13] = Some(Element::Fm); // (7,14)
    table[6][14] = Some(Element::Md); // (7,15)
    table[6][15] = Some(Element::No); // (7,16)
    table[6][16] = Some(Element::Lr); // (7,17)
    table[6][17] = Some(Element::Rf); // (7,18)
    table[6][18] = Some(Element::Db); // (7,19)
    table[6][19] = Some(Element::Sg); // (7,20)
    table[6][20] = Some(Element::Bh); // (7,21)
    table[6][21] = Some(Element::Hs); // (7,22)
    table[6][22] = Some(Element::Mt); // (7,23)
    table[6][23] = Some(Element::Ds); // (7,24)
    table[6][24] = Some(Element::Rg); // (7,25)
    table[6][25] = Some(Element::Cn); // (7,26)
    table[6][26] = Some(Element::Nh); // (7,27)
    table[6][27] = Some(Element::Fl); // (7,28)
    table[6][28] = Some(Element::Mc); // (7,29)
    table[6][29] = Some(Element::Lv); // (7,30)
    table[6][30] = Some(Element::Ts); // (7,31)
    table[6][31] = Some(Element::Og); // (7,32)

    table
};

// Helper function to normalize element symbol bytes to title case on stack
// Returns the normalized buffer and the length of the symbol (1 or 2).
fn normalize_symbol_bytes(bytes: &[u8]) -> Option<([u8; 2], usize)> {
    match bytes.len() {
        1 => {
            if !bytes[0].is_ascii_alphabetic() {
                return None;
            }
            let upper_b = bytes[0].to_ascii_uppercase();
            Some(([upper_b, 0], 1)) // Second byte is padding, length is 1
        }
        2 => {
            if !bytes[0].is_ascii_alphabetic() || !bytes[1].is_ascii_alphabetic() {
                return None;
            }
            let b1_upper = bytes[0].to_ascii_uppercase();
            let b2_lower = bytes[1].to_ascii_lowercase();
            Some(([b1_upper, b2_lower], 2))
        }
        _ => None,
    }
}

impl Element {
    // Get element from symbol byte string (allocation-free)
    pub fn from_symbol_bytes(symbol: &[u8]) -> Option<Self> {
        if let Some((key_buf, len)) = normalize_symbol_bytes(symbol) {
            match &key_buf[..len] {
                b"H" => Some(Element::H),
                b"He" => Some(Element::He),
                b"Li" => Some(Element::Li),
                b"Be" => Some(Element::Be),
                b"B" => Some(Element::B),
                b"C" => Some(Element::C),
                b"N" => Some(Element::N),
                b"O" => Some(Element::O),
                b"F" => Some(Element::F),
                b"Ne" => Some(Element::Ne),
                b"Na" => Some(Element::Na),
                b"Mg" => Some(Element::Mg),
                b"Al" => Some(Element::Al),
                b"Si" => Some(Element::Si),
                b"P" => Some(Element::P),
                b"S" => Some(Element::S),
                b"Cl" => Some(Element::Cl),
                b"Ar" => Some(Element::Ar),
                b"K" => Some(Element::K),
                b"Ca" => Some(Element::Ca),
                b"Sc" => Some(Element::Sc),
                b"Ti" => Some(Element::Ti),
                b"V" => Some(Element::V),
                b"Cr" => Some(Element::Cr),
                b"Mn" => Some(Element::Mn),
                b"Fe" => Some(Element::Fe),
                b"Co" => Some(Element::Co),
                b"Ni" => Some(Element::Ni),
                b"Cu" => Some(Element::Cu),
                b"Zn" => Some(Element::Zn),
                b"Ga" => Some(Element::Ga),
                b"Ge" => Some(Element::Ge),
                b"As" => Some(Element::As),
                b"Se" => Some(Element::Se),
                b"Br" => Some(Element::Br),
                b"Kr" => Some(Element::Kr),
                b"Rb" => Some(Element::Rb),
                b"Sr" => Some(Element::Sr),
                b"Y" => Some(Element::Y),
                b"Zr" => Some(Element::Zr),
                b"Nb" => Some(Element::Nb),
                b"Mo" => Some(Element::Mo),
                b"Tc" => Some(Element::Tc),
                b"Ru" => Some(Element::Ru),
                b"Rh" => Some(Element::Rh),
                b"Pd" => Some(Element::Pd),
                b"Ag" => Some(Element::Ag),
                b"Cd" => Some(Element::Cd),
                b"In" => Some(Element::In),
                b"Sn" => Some(Element::Sn),
                b"Sb" => Some(Element::Sb),
                b"Te" => Some(Element::Te),
                b"I" => Some(Element::I),
                b"Xe" => Some(Element::Xe),
                b"Cs" => Some(Element::Cs),
                b"Ba" => Some(Element::Ba),
                b"La" => Some(Element::La),
                b"Ce" => Some(Element::Ce),
                b"Pr" => Some(Element::Pr),
                b"Nd" => Some(Element::Nd),
                b"Pm" => Some(Element::Pm),
                b"Sm" => Some(Element::Sm),
                b"Eu" => Some(Element::Eu),
                b"Gd" => Some(Element::Gd),
                b"Tb" => Some(Element::Tb),
                b"Dy" => Some(Element::Dy),
                b"Ho" => Some(Element::Ho),
                b"Er" => Some(Element::Er),
                b"Tm" => Some(Element::Tm),
                b"Yb" => Some(Element::Yb),
                b"Lu" => Some(Element::Lu),
                b"Hf" => Some(Element::Hf),
                b"Ta" => Some(Element::Ta),
                b"W" => Some(Element::W),
                b"Re" => Some(Element::Re),
                b"Os" => Some(Element::Os),
                b"Ir" => Some(Element::Ir),
                b"Pt" => Some(Element::Pt),
                b"Au" => Some(Element::Au),
                b"Hg" => Some(Element::Hg),
                b"Tl" => Some(Element::Tl),
                b"Pb" => Some(Element::Pb),
                b"Bi" => Some(Element::Bi),
                b"Po" => Some(Element::Po),
                b"At" => Some(Element::At),
                b"Rn" => Some(Element::Rn),
                b"Fr" => Some(Element::Fr),
                b"Ra" => Some(Element::Ra),
                b"Ac" => Some(Element::Ac),
                b"Th" => Some(Element::Th),
                b"Pa" => Some(Element::Pa),
                b"U" => Some(Element::U),
                b"Np" => Some(Element::Np),
                b"Pu" => Some(Element::Pu),
                b"Am" => Some(Element::Am),
                b"Cm" => Some(Element::Cm),
                b"Bk" => Some(Element::Bk),
                b"Cf" => Some(Element::Cf),
                b"Es" => Some(Element::Es),
                b"Fm" => Some(Element::Fm),
                b"Md" => Some(Element::Md),
                b"No" => Some(Element::No),
                b"Lr" => Some(Element::Lr),
                b"Rf" => Some(Element::Rf),
                b"Db" => Some(Element::Db),
                b"Sg" => Some(Element::Sg),
                b"Bh" => Some(Element::Bh),
                b"Hs" => Some(Element::Hs),
                b"Mt" => Some(Element::Mt),
                b"Ds" => Some(Element::Ds),
                b"Rg" => Some(Element::Rg),
                b"Cn" => Some(Element::Cn),
                b"Nh" => Some(Element::Nh),
                b"Fl" => Some(Element::Fl),
                b"Mc" => Some(Element::Mc),
                b"Lv" => Some(Element::Lv),
                b"Ts" => Some(Element::Ts),
                b"Og" => Some(Element::Og),
                _ => None,
            }
        } else {
            None
        }
    }

    // Get element from symbol string (allocation-free)
    pub fn from_symbol(symbol: &str) -> Option<Self> {
        Self::from_symbol_bytes(symbol.as_bytes())
    }

    // Get element from atomic number
    pub const fn from_atomic_number(number: u8) -> Option<Self> {
        if number >= 1 && number <= 118 {
            Some(ELEMENTS[(number - 1) as usize])
        } else {
            None
        }
    }

    // Get element from period and group coordinates
    // NOTE: 32-group layout is used for unique group assignment.
    pub const fn from_period_group(period: u8, group: u8) -> Option<Self> {
        if period == 0 || period > MAX_PERIOD_NUMBER || group == 0 || group > 32 {
            None
        } else {
            PERIOD_GROUP_TO_ELEMENT[period as usize - 1][group as usize - 1]
        }
    }

    // Get atomic number for element
    pub const fn atomic_number(&self) -> u8 {
        let index = *self as usize;
        ELEMENT_DATA[index].0
    }

    // Get mass number of the reference isotope for element
    pub const fn reference_mass_number(&self) -> u32 {
        let index = *self as usize;
        ELEMENT_DATA[index].1
    }

    // Get atomic mass (standard atomic weight) for element
    pub const fn mass(&self) -> f64 {
        let index = *self as usize;
        ELEMENT_DATA[index].2
    }

    // Get symbol for element
    pub const fn symbol(&self) -> &'static str {
        let index = *self as usize;
        ELEMENT_DATA[index].3
    }

    // Get period for element
    pub const fn period(&self) -> u8 {
        let index = *self as usize;
        ELEMENT_DATA[index].4
    }

    // Get group number for element.
    // NOTE: 32-group layout is used for unique group assignment.
    pub const fn group(&self) -> u8 {
        let index = *self as usize;
        ELEMENT_DATA[index].5
    }

    // Get number of valence electrons for element
    pub const fn valence_electrons(&self) -> u8 {
        let index = *self as usize;
        ELEMENT_DATA[index].6
    }

    // Get max valence for element
    pub const fn max_valence(&self) -> u8 {
        let index = *self as usize;
        ELEMENT_DATA[index].7
    }

    // Get minimum and maximum charge for element
    pub const fn charge_bounds(&self) -> (i8, i8) {
        let index = *self as usize;
        ELEMENT_DATA[index].8
    }

    // Get maximum number of unpaired electrons for element
    pub const fn max_unpaired_electrons(&self) -> u8 {
        let index = *self as usize;
        ELEMENT_DATA[index].9
    }

    // Get maximum number of implicit hydrogens for element
    pub const fn max_implicit_hydrogens(&self) -> u8 {
        let index = *self as usize;
        ELEMENT_DATA[index].10
    }

    // Get next element in the periodic table
    pub const fn next(&self) -> Option<Self> {
        let atomic_number = self.atomic_number();
        if atomic_number < MAX_ATOMIC_NUMBER {
            Self::from_atomic_number(atomic_number + 1)
        } else {
            None
        }
    }

    // Get previous element in the periodic table
    pub const fn previous(&self) -> Option<Self> {
        let atomic_number = self.atomic_number();
        if atomic_number > 1 {
            Self::from_atomic_number(atomic_number - 1)
        } else {
            None
        }
    }

    // Get element in the next period (same group)
    pub const fn next_period(&self) -> Option<Self> {
        if self.period() == MAX_PERIOD_NUMBER {
            None
        } else {
            Self::from_period_group(self.period() + 1, self.group())
        }
    }

    // Get element in the previous period (same group)
    pub const fn previous_period(&self) -> Option<Self> {
        if self.period() == MIN_PERIOD_NUMBER[self.group() as usize] {
            None
        } else {
            Self::from_period_group(self.period() - 1, self.group())
        }
    }

    // Get element in the next group (same period)
    pub const fn next_group(&self) -> Option<Self> {
        if self.group() == MAX_GROUP_NUMBER[self.period() as usize] {
            None
        } else {
            Self::from_period_group(self.period(), self.group() + 1)
        }
    }

    // Get element in the previous group (same period)
    pub const fn previous_group(&self) -> Option<Self> {
        if self.group() == 1 {
            None
        } else {
            Self::from_period_group(self.period(), self.group() - 1)
        }
    }

    /// Get all elements
    pub const fn all() -> &'static [Self] {
        &ELEMENTS
    }

    /// Check if byte string contains valid element symbol using direct indexing
    pub fn is_element_bytes(symbol: &[u8]) -> bool {
        if let Some((key_buf, len)) = normalize_symbol_bytes(symbol) {
            match len {
                1 => {
                    let first = key_buf[0];
                    #[allow(clippy::manual_range_contains)]
                    if first >= b'A' && first <= b'Z' {
                        let index = (first - b'A') as usize;
                        ELEMENT_LOOKUP[index]
                    } else {
                        false
                    }
                }
                2 => {
                    let first = key_buf[0];
                    let second = key_buf[1];
                    #[allow(clippy::manual_range_contains)]
                    if first >= b'A' && first <= b'Z' && second >= b'a' && second <= b'z' {
                        let index = 26 + (first - b'A') as usize * 26 + (second - b'a') as usize;
                        ELEMENT_LOOKUP[index]
                    } else {
                        false
                    }
                }
                _ => false,
            }
        } else {
            false
        }
    }

    // Check if string contains valid element symbol (case-insensitive)
    pub fn is_element(symbol: &str) -> bool {
        Self::is_element_bytes(symbol.as_bytes())
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
    use float_cmp::*;
    use rstest::*;
    use serde_json;

    #[test]
    fn test_element_from_symbol_bytes() {
        assert_eq!(Element::from_symbol_bytes(b"H"), Some(Element::H));
        assert_eq!(Element::from_symbol_bytes(b"h"), Some(Element::H));
        assert_eq!(Element::from_symbol_bytes(b"He"), Some(Element::He));
        assert_eq!(Element::from_symbol_bytes(b"he"), Some(Element::He));
        assert_eq!(Element::from_symbol_bytes(b"HE"), Some(Element::He));
        assert_eq!(Element::from_symbol_bytes(b"hE"), Some(Element::He));
        assert_eq!(Element::from_symbol_bytes(b"Li"), Some(Element::Li));
        assert_eq!(Element::from_symbol_bytes(b"Be"), Some(Element::Be));
        assert_eq!(Element::from_symbol_bytes(b"B"), Some(Element::B));
        assert_eq!(Element::from_symbol_bytes(b"C"), Some(Element::C));
        assert_eq!(Element::from_symbol_bytes(b"N"), Some(Element::N));
        assert_eq!(Element::from_symbol_bytes(b"O"), Some(Element::O));
        assert_eq!(Element::from_symbol_bytes(b"F"), Some(Element::F));
        assert_eq!(Element::from_symbol_bytes(b"Ne"), Some(Element::Ne));
        assert_eq!(Element::from_symbol_bytes(b"Na"), Some(Element::Na));
        assert_eq!(Element::from_symbol_bytes(b"Mg"), Some(Element::Mg));
        assert_eq!(Element::from_symbol_bytes(b"Al"), Some(Element::Al));
        assert_eq!(Element::from_symbol_bytes(b"Si"), Some(Element::Si));
        assert_eq!(Element::from_symbol_bytes(b"P"), Some(Element::P));
        assert_eq!(Element::from_symbol_bytes(b"S"), Some(Element::S));
        assert_eq!(Element::from_symbol_bytes(b"Cl"), Some(Element::Cl));
        assert_eq!(Element::from_symbol_bytes(b"Ar"), Some(Element::Ar));
        assert_eq!(Element::from_symbol_bytes(b"K"), Some(Element::K));
        assert_eq!(Element::from_symbol_bytes(b"Ca"), Some(Element::Ca));
        assert_eq!(Element::from_symbol_bytes(b"Sc"), Some(Element::Sc));
        assert_eq!(Element::from_symbol_bytes(b"Ti"), Some(Element::Ti));
        assert_eq!(Element::from_symbol_bytes(b"V"), Some(Element::V));
        assert_eq!(Element::from_symbol_bytes(b"Cr"), Some(Element::Cr));
        assert_eq!(Element::from_symbol_bytes(b"Mn"), Some(Element::Mn));
        assert_eq!(Element::from_symbol_bytes(b"Fe"), Some(Element::Fe));
        assert_eq!(Element::from_symbol_bytes(b"Co"), Some(Element::Co));
        assert_eq!(Element::from_symbol_bytes(b"Ni"), Some(Element::Ni));
        assert_eq!(Element::from_symbol_bytes(b"Cu"), Some(Element::Cu));
        assert_eq!(Element::from_symbol_bytes(b"Zn"), Some(Element::Zn));
        assert_eq!(Element::from_symbol_bytes(b"Ga"), Some(Element::Ga));
        assert_eq!(Element::from_symbol_bytes(b"Ge"), Some(Element::Ge));
        assert_eq!(Element::from_symbol_bytes(b"As"), Some(Element::As));
        assert_eq!(Element::from_symbol_bytes(b"Se"), Some(Element::Se));
        assert_eq!(Element::from_symbol_bytes(b"Br"), Some(Element::Br));
        assert_eq!(Element::from_symbol_bytes(b"Kr"), Some(Element::Kr));
        assert_eq!(Element::from_symbol_bytes(b"Rb"), Some(Element::Rb));
        assert_eq!(Element::from_symbol_bytes(b"Sr"), Some(Element::Sr));
        assert_eq!(Element::from_symbol_bytes(b"Y"), Some(Element::Y));
        assert_eq!(Element::from_symbol_bytes(b"Zr"), Some(Element::Zr));
        assert_eq!(Element::from_symbol_bytes(b"Nb"), Some(Element::Nb));
        assert_eq!(Element::from_symbol_bytes(b"Mo"), Some(Element::Mo));
        assert_eq!(Element::from_symbol_bytes(b"Tc"), Some(Element::Tc));
        assert_eq!(Element::from_symbol_bytes(b"Ru"), Some(Element::Ru));
        assert_eq!(Element::from_symbol_bytes(b"Rh"), Some(Element::Rh));
        assert_eq!(Element::from_symbol_bytes(b"Pd"), Some(Element::Pd));
        assert_eq!(Element::from_symbol_bytes(b"Ag"), Some(Element::Ag));
        assert_eq!(Element::from_symbol_bytes(b"Cd"), Some(Element::Cd));
        assert_eq!(Element::from_symbol_bytes(b"In"), Some(Element::In));
        assert_eq!(Element::from_symbol_bytes(b"Sn"), Some(Element::Sn));
        assert_eq!(Element::from_symbol_bytes(b"Sb"), Some(Element::Sb));
        assert_eq!(Element::from_symbol_bytes(b"Te"), Some(Element::Te));
        assert_eq!(Element::from_symbol_bytes(b"I"), Some(Element::I));
        assert_eq!(Element::from_symbol_bytes(b"Xe"), Some(Element::Xe));
        assert_eq!(Element::from_symbol_bytes(b"Cs"), Some(Element::Cs));
        assert_eq!(Element::from_symbol_bytes(b"Ba"), Some(Element::Ba));
        assert_eq!(Element::from_symbol_bytes(b"La"), Some(Element::La));
        assert_eq!(Element::from_symbol_bytes(b"Ce"), Some(Element::Ce));
        assert_eq!(Element::from_symbol_bytes(b"Pr"), Some(Element::Pr));
        assert_eq!(Element::from_symbol_bytes(b"Nd"), Some(Element::Nd));
        assert_eq!(Element::from_symbol_bytes(b"Pm"), Some(Element::Pm));
        assert_eq!(Element::from_symbol_bytes(b"Sm"), Some(Element::Sm));
        assert_eq!(Element::from_symbol_bytes(b"Eu"), Some(Element::Eu));
        assert_eq!(Element::from_symbol_bytes(b"Gd"), Some(Element::Gd));
        assert_eq!(Element::from_symbol_bytes(b"Tb"), Some(Element::Tb));
        assert_eq!(Element::from_symbol_bytes(b"Dy"), Some(Element::Dy));
        assert_eq!(Element::from_symbol_bytes(b"Ho"), Some(Element::Ho));
        assert_eq!(Element::from_symbol_bytes(b"Er"), Some(Element::Er));
        assert_eq!(Element::from_symbol_bytes(b"Tm"), Some(Element::Tm));
        assert_eq!(Element::from_symbol_bytes(b"Yb"), Some(Element::Yb));
        assert_eq!(Element::from_symbol_bytes(b"Lu"), Some(Element::Lu));
        assert_eq!(Element::from_symbol_bytes(b"Hf"), Some(Element::Hf));
        assert_eq!(Element::from_symbol_bytes(b"Ta"), Some(Element::Ta));
        assert_eq!(Element::from_symbol_bytes(b"W"), Some(Element::W));
        assert_eq!(Element::from_symbol_bytes(b"Re"), Some(Element::Re));
        assert_eq!(Element::from_symbol_bytes(b"Os"), Some(Element::Os));
        assert_eq!(Element::from_symbol_bytes(b"Ir"), Some(Element::Ir));
        assert_eq!(Element::from_symbol_bytes(b"Pt"), Some(Element::Pt));
        assert_eq!(Element::from_symbol_bytes(b"Au"), Some(Element::Au));
        assert_eq!(Element::from_symbol_bytes(b"Hg"), Some(Element::Hg));
        assert_eq!(Element::from_symbol_bytes(b"Tl"), Some(Element::Tl));
        assert_eq!(Element::from_symbol_bytes(b"Pb"), Some(Element::Pb));
        assert_eq!(Element::from_symbol_bytes(b"Bi"), Some(Element::Bi));
        assert_eq!(Element::from_symbol_bytes(b"Po"), Some(Element::Po));
        assert_eq!(Element::from_symbol_bytes(b"At"), Some(Element::At));
        assert_eq!(Element::from_symbol_bytes(b"Rn"), Some(Element::Rn));
        assert_eq!(Element::from_symbol_bytes(b"Fr"), Some(Element::Fr));
        assert_eq!(Element::from_symbol_bytes(b"Ra"), Some(Element::Ra));
        assert_eq!(Element::from_symbol_bytes(b"Ac"), Some(Element::Ac));
        assert_eq!(Element::from_symbol_bytes(b"Th"), Some(Element::Th));
        assert_eq!(Element::from_symbol_bytes(b"Pa"), Some(Element::Pa));
        assert_eq!(Element::from_symbol_bytes(b"U"), Some(Element::U));
        assert_eq!(Element::from_symbol_bytes(b"Np"), Some(Element::Np));
        assert_eq!(Element::from_symbol_bytes(b"Pu"), Some(Element::Pu));
        assert_eq!(Element::from_symbol_bytes(b"Am"), Some(Element::Am));
        assert_eq!(Element::from_symbol_bytes(b"Cm"), Some(Element::Cm));
        assert_eq!(Element::from_symbol_bytes(b"Bk"), Some(Element::Bk));
        assert_eq!(Element::from_symbol_bytes(b"Cf"), Some(Element::Cf));
        assert_eq!(Element::from_symbol_bytes(b"Es"), Some(Element::Es));
        assert_eq!(Element::from_symbol_bytes(b"Fm"), Some(Element::Fm));
        assert_eq!(Element::from_symbol_bytes(b"Md"), Some(Element::Md));
        assert_eq!(Element::from_symbol_bytes(b"No"), Some(Element::No));
        assert_eq!(Element::from_symbol_bytes(b"Lr"), Some(Element::Lr));
        assert_eq!(Element::from_symbol_bytes(b"Rf"), Some(Element::Rf));
        assert_eq!(Element::from_symbol_bytes(b"Db"), Some(Element::Db));
        assert_eq!(Element::from_symbol_bytes(b"Sg"), Some(Element::Sg));
        assert_eq!(Element::from_symbol_bytes(b"Bh"), Some(Element::Bh));
        assert_eq!(Element::from_symbol_bytes(b"Hs"), Some(Element::Hs));
        assert_eq!(Element::from_symbol_bytes(b"Mt"), Some(Element::Mt));
        assert_eq!(Element::from_symbol_bytes(b"Ds"), Some(Element::Ds));
        assert_eq!(Element::from_symbol_bytes(b"Rg"), Some(Element::Rg));
        assert_eq!(Element::from_symbol_bytes(b"Cn"), Some(Element::Cn));
        assert_eq!(Element::from_symbol_bytes(b"Nh"), Some(Element::Nh));
        assert_eq!(Element::from_symbol_bytes(b"Fl"), Some(Element::Fl));
        assert_eq!(Element::from_symbol_bytes(b"Mc"), Some(Element::Mc));
        assert_eq!(Element::from_symbol_bytes(b"Lv"), Some(Element::Lv));
        assert_eq!(Element::from_symbol_bytes(b"Ts"), Some(Element::Ts));
        assert_eq!(Element::from_symbol_bytes(b"Og"), Some(Element::Og));
        assert_eq!(Element::from_symbol_bytes(b"invalid"), None);
    }

    #[test]
    fn test_element_from_symbol() {
        assert_eq!(Element::from_symbol("H"), Some(Element::H));
        assert_eq!(Element::from_symbol("h"), Some(Element::H));
        assert_eq!(Element::from_symbol("He"), Some(Element::He));
        assert_eq!(Element::from_symbol("he"), Some(Element::He));
        assert_eq!(Element::from_symbol("HE"), Some(Element::He));
        assert_eq!(Element::from_symbol("hE"), Some(Element::He));
        assert_eq!(Element::from_symbol("Li"), Some(Element::Li));
        assert_eq!(Element::from_symbol("Be"), Some(Element::Be));
        assert_eq!(Element::from_symbol("B"), Some(Element::B));
        assert_eq!(Element::from_symbol("C"), Some(Element::C));
        assert_eq!(Element::from_symbol("N"), Some(Element::N));
        assert_eq!(Element::from_symbol("O"), Some(Element::O));
        assert_eq!(Element::from_symbol("F"), Some(Element::F));
        assert_eq!(Element::from_symbol("Ne"), Some(Element::Ne));
        assert_eq!(Element::from_symbol("Na"), Some(Element::Na));
        assert_eq!(Element::from_symbol("Mg"), Some(Element::Mg));
        assert_eq!(Element::from_symbol("Al"), Some(Element::Al));
        assert_eq!(Element::from_symbol("Si"), Some(Element::Si));
        assert_eq!(Element::from_symbol("P"), Some(Element::P));
        assert_eq!(Element::from_symbol("S"), Some(Element::S));
        assert_eq!(Element::from_symbol("Cl"), Some(Element::Cl));
        assert_eq!(Element::from_symbol("Ar"), Some(Element::Ar));
        assert_eq!(Element::from_symbol("K"), Some(Element::K));
        assert_eq!(Element::from_symbol("Ca"), Some(Element::Ca));
        assert_eq!(Element::from_symbol("Sc"), Some(Element::Sc));
        assert_eq!(Element::from_symbol("Ti"), Some(Element::Ti));
        assert_eq!(Element::from_symbol("V"), Some(Element::V));
        assert_eq!(Element::from_symbol("Cr"), Some(Element::Cr));
        assert_eq!(Element::from_symbol("Mn"), Some(Element::Mn));
        assert_eq!(Element::from_symbol("Fe"), Some(Element::Fe));
        assert_eq!(Element::from_symbol("Co"), Some(Element::Co));
        assert_eq!(Element::from_symbol("Ni"), Some(Element::Ni));
        assert_eq!(Element::from_symbol("Cu"), Some(Element::Cu));
        assert_eq!(Element::from_symbol("Zn"), Some(Element::Zn));
        assert_eq!(Element::from_symbol("Ga"), Some(Element::Ga));
        assert_eq!(Element::from_symbol("Ge"), Some(Element::Ge));
        assert_eq!(Element::from_symbol("As"), Some(Element::As));
        assert_eq!(Element::from_symbol("Se"), Some(Element::Se));
        assert_eq!(Element::from_symbol("Br"), Some(Element::Br));
        assert_eq!(Element::from_symbol("Kr"), Some(Element::Kr));
        assert_eq!(Element::from_symbol("Rb"), Some(Element::Rb));
        assert_eq!(Element::from_symbol("Sr"), Some(Element::Sr));
        assert_eq!(Element::from_symbol("Y"), Some(Element::Y));
        assert_eq!(Element::from_symbol("Zr"), Some(Element::Zr));
        assert_eq!(Element::from_symbol("Nb"), Some(Element::Nb));
        assert_eq!(Element::from_symbol("Mo"), Some(Element::Mo));
        assert_eq!(Element::from_symbol("Tc"), Some(Element::Tc));
        assert_eq!(Element::from_symbol("Ru"), Some(Element::Ru));
        assert_eq!(Element::from_symbol("Rh"), Some(Element::Rh));
        assert_eq!(Element::from_symbol("Pd"), Some(Element::Pd));
        assert_eq!(Element::from_symbol("Ag"), Some(Element::Ag));
        assert_eq!(Element::from_symbol("Cd"), Some(Element::Cd));
        assert_eq!(Element::from_symbol("In"), Some(Element::In));
        assert_eq!(Element::from_symbol("Sn"), Some(Element::Sn));
        assert_eq!(Element::from_symbol("Sb"), Some(Element::Sb));
        assert_eq!(Element::from_symbol("Te"), Some(Element::Te));
        assert_eq!(Element::from_symbol("I"), Some(Element::I));
        assert_eq!(Element::from_symbol("Xe"), Some(Element::Xe));
        assert_eq!(Element::from_symbol("Cs"), Some(Element::Cs));
        assert_eq!(Element::from_symbol("Ba"), Some(Element::Ba));
        assert_eq!(Element::from_symbol("La"), Some(Element::La));
        assert_eq!(Element::from_symbol("Ce"), Some(Element::Ce));
        assert_eq!(Element::from_symbol("Pr"), Some(Element::Pr));
        assert_eq!(Element::from_symbol("Nd"), Some(Element::Nd));
        assert_eq!(Element::from_symbol("Pm"), Some(Element::Pm));
        assert_eq!(Element::from_symbol("Sm"), Some(Element::Sm));
        assert_eq!(Element::from_symbol("Eu"), Some(Element::Eu));
        assert_eq!(Element::from_symbol("Gd"), Some(Element::Gd));
        assert_eq!(Element::from_symbol("Tb"), Some(Element::Tb));
        assert_eq!(Element::from_symbol("Dy"), Some(Element::Dy));
        assert_eq!(Element::from_symbol("Ho"), Some(Element::Ho));
        assert_eq!(Element::from_symbol("Er"), Some(Element::Er));
        assert_eq!(Element::from_symbol("Tm"), Some(Element::Tm));
        assert_eq!(Element::from_symbol("Yb"), Some(Element::Yb));
        assert_eq!(Element::from_symbol("Lu"), Some(Element::Lu));
        assert_eq!(Element::from_symbol("Hf"), Some(Element::Hf));
        assert_eq!(Element::from_symbol("Ta"), Some(Element::Ta));
        assert_eq!(Element::from_symbol("W"), Some(Element::W));
        assert_eq!(Element::from_symbol("Re"), Some(Element::Re));
        assert_eq!(Element::from_symbol("Os"), Some(Element::Os));
        assert_eq!(Element::from_symbol("Ir"), Some(Element::Ir));
        assert_eq!(Element::from_symbol("Pt"), Some(Element::Pt));
        assert_eq!(Element::from_symbol("Au"), Some(Element::Au));
        assert_eq!(Element::from_symbol("Hg"), Some(Element::Hg));
        assert_eq!(Element::from_symbol("Tl"), Some(Element::Tl));
        assert_eq!(Element::from_symbol("Pb"), Some(Element::Pb));
        assert_eq!(Element::from_symbol("Bi"), Some(Element::Bi));
        assert_eq!(Element::from_symbol("Po"), Some(Element::Po));
        assert_eq!(Element::from_symbol("At"), Some(Element::At));
        assert_eq!(Element::from_symbol("Rn"), Some(Element::Rn));
        assert_eq!(Element::from_symbol("Fr"), Some(Element::Fr));
        assert_eq!(Element::from_symbol("Ra"), Some(Element::Ra));
        assert_eq!(Element::from_symbol("Ac"), Some(Element::Ac));
        assert_eq!(Element::from_symbol("Th"), Some(Element::Th));
        assert_eq!(Element::from_symbol("Pa"), Some(Element::Pa));
        assert_eq!(Element::from_symbol("U"), Some(Element::U));
        assert_eq!(Element::from_symbol("Np"), Some(Element::Np));
        assert_eq!(Element::from_symbol("Pu"), Some(Element::Pu));
        assert_eq!(Element::from_symbol("Am"), Some(Element::Am));
        assert_eq!(Element::from_symbol("Cm"), Some(Element::Cm));
        assert_eq!(Element::from_symbol("Bk"), Some(Element::Bk));
        assert_eq!(Element::from_symbol("Cf"), Some(Element::Cf));
        assert_eq!(Element::from_symbol("Es"), Some(Element::Es));
        assert_eq!(Element::from_symbol("Fm"), Some(Element::Fm));
        assert_eq!(Element::from_symbol("Md"), Some(Element::Md));
        assert_eq!(Element::from_symbol("No"), Some(Element::No));
        assert_eq!(Element::from_symbol("Lr"), Some(Element::Lr));
        assert_eq!(Element::from_symbol("Rf"), Some(Element::Rf));
        assert_eq!(Element::from_symbol("Db"), Some(Element::Db));
        assert_eq!(Element::from_symbol("Sg"), Some(Element::Sg));
        assert_eq!(Element::from_symbol("Bh"), Some(Element::Bh));
        assert_eq!(Element::from_symbol("Hs"), Some(Element::Hs));
        assert_eq!(Element::from_symbol("Mt"), Some(Element::Mt));
        assert_eq!(Element::from_symbol("Ds"), Some(Element::Ds));
        assert_eq!(Element::from_symbol("Rg"), Some(Element::Rg));
        assert_eq!(Element::from_symbol("Cn"), Some(Element::Cn));
        assert_eq!(Element::from_symbol("Nh"), Some(Element::Nh));
        assert_eq!(Element::from_symbol("Fl"), Some(Element::Fl));
        assert_eq!(Element::from_symbol("Mc"), Some(Element::Mc));
        assert_eq!(Element::from_symbol("Lv"), Some(Element::Lv));
        assert_eq!(Element::from_symbol("Ts"), Some(Element::Ts));
        assert_eq!(Element::from_symbol("Og"), Some(Element::Og));
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
    #[case(Element::H, "H", 1, 1, 1.008, 1, 1, 1, 2)]
    #[case(Element::He, "He", 2, 4, 4.0026, 1, 32, 0, 0)]
    fn test_element_properties(
        #[case] element: Element,
        #[case] symbol: &str,
        #[case] atomic_number: u8,
        #[case] reference_atomic_mass: u32,
        #[case] atomic_mass: f64,
        #[case] period: u8,
        #[case] group: u8,
        #[case] valence_electrons: u8,
        #[case] max_valence: u8,
    ) {
        assert_eq!(element.symbol(), symbol);
        assert_eq!(element.atomic_number(), atomic_number);
        assert_eq!(element.reference_mass_number(), reference_atomic_mass);
        assert!(approx_eq!(f64, element.mass(), atomic_mass, ulps = 4));
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
        for data in &ELEMENT_DATA {
            assert!(data.0 > 0); // atomic number
            assert!(data.2 > 0.0); // atomic mass
            assert!(!data.3.is_empty()); // symbol
            assert!(data.4 <= 7); // period
            assert!(data.5 <= 32); // group
            assert!(data.6 <= 32); // valence electrons
            assert!(data.7 <= 32); // max valence
            assert!(data.8 .0 <= data.8 .1); // charge bounds
            assert!(data.9 <= 10); // max unpaired electrons
            assert!(data.10 <= 4); // max implicit hydrogens
        }
    }

    #[test]
    fn test_atomic_number_to_element_to_atomic_number() {
        for (i, element) in ELEMENTS.iter().enumerate() {
            let atomic_number = (i + 1) as u8;
            assert_eq!(element.atomic_number(), atomic_number);
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
    fn test_period_group_lookup() {
        assert_eq!(Element::from_period_group(1, 1), Some(Element::H));
        assert_eq!(Element::from_period_group(1, 32), Some(Element::He));
        assert_eq!(Element::from_period_group(4, 22), Some(Element::Fe));
        assert_eq!(Element::from_period_group(6, 24), Some(Element::Pt));
        assert_eq!(Element::from_period_group(7, 32), Some(Element::Og));

        assert_eq!(Element::from_period_group(0, 1), None);
        assert_eq!(Element::from_period_group(8, 1), None);
        assert_eq!(Element::from_period_group(1, 0), None);
        assert_eq!(Element::from_period_group(1, 33), None);

        assert_eq!(Element::from_period_group(2, 4), None);
        assert_eq!(Element::from_period_group(3, 4), None);
    }

    #[test]
    fn test_element_next_previous_period() {
        assert_eq!(Element::H.previous_period(), None);
        assert_eq!(Element::H.next_period(), Some(Element::Li));
        assert_eq!(Element::C.previous_period(), None);
        assert_eq!(Element::C.next_period(), Some(Element::Si));
        assert_eq!(Element::U.previous_period(), Some(Element::Nd));
        assert_eq!(Element::U.next_period(), None);
    }

    #[test]
    fn test_element_next_previous_group() {
        assert_eq!(Element::Li.previous_group(), None);
        assert_eq!(Element::Li.next_group(), Some(Element::Be));
        assert_eq!(Element::Ne.previous_group(), Some(Element::F));
        assert_eq!(Element::Ne.next_group(), None);
    }

    #[test]
    fn test_element_is_element_bytes() {
        assert!(Element::is_element_bytes(b"H"));
        assert!(Element::is_element_bytes(b"He"));
        assert!(Element::is_element_bytes(b"C"));
        assert!(Element::is_element_bytes(b"c"));
        assert!(Element::is_element_bytes(b"Cu"));
        assert!(Element::is_element_bytes(b"cu"));
        assert!(Element::is_element_bytes(b"CU"));
        assert!(Element::is_element_bytes(b"Ru"));
        assert!(Element::is_element_bytes(b"U"));
        assert!(Element::is_element_bytes(b"u"));
        assert!(!Element::is_element_bytes(b"R1"));
        assert!(!Element::is_element_bytes(b"X"));
    }

    #[test]
    fn test_element_is_element() {
        assert!(Element::is_element("H"));
        assert!(Element::is_element("He"));
        assert!(Element::is_element("C"));
        assert!(Element::is_element("c"));
        assert!(Element::is_element("Cu"));
        assert!(Element::is_element("cu"));
        assert!(Element::is_element("CU"));
        assert!(Element::is_element("Ru"));
        assert!(Element::is_element("U"));
        assert!(Element::is_element("u"));
        assert!(!Element::is_element("R1"));
        assert!(!Element::is_element("X"));
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
