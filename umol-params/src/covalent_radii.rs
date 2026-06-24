//! Covalent radii from Pyykko & Atsumi.
//!
//! Single-bond: Pyykko & Atsumi, Chem. Eur. J. 2009, 15, 186-197.
//! Double-bond: Pyykko & Atsumi, Chem. Eur. J. 2009, 15, 12770-12779.
//! Triple-bond: Pyykko & Atsumi, Chem. Eur. J. 2005, 11, 3511-3520.

use umol_chem::element::Element;
use umol_chem::units::length::Length;

/// Covalent radii for single, double, and triple bonds.
/// Values stored in atomic units (Bohr). `None` means no data available.
#[derive(Debug, Clone, Copy)]
pub struct CovalentRadii {
    pub single: Length,
    pub double: Option<Length>,
    pub triple: Option<Length>,
}

impl CovalentRadii {
    /// Returns the radius for a given bond order (1, 2, or 3).
    /// Falls back to the next available radius if the requested order is unavailable.
    pub fn for_order(&self, order: u8) -> Option<Length> {
        match order {
            1 => Some(self.single),
            2 => self.double.or(Some(self.single)),
            3 => self.triple.or(self.double).or(Some(self.single)),
            _ => None,
        }
    }
}

/// Look up covalent radii for an element.
pub fn covalent_radii(element: Element) -> CovalentRadii {
    let z = element.atomic_number() as usize;
    if z == 0 || z > RADII.len() {
        return CovalentRadii {
            single: Length::picometer(150.0),
            double: None,
            triple: None,
        };
    }
    RADII[z - 1]
}

const fn pm(v: u16) -> Length {
    Length::picometer(v as f64)
}

const fn cr(single: u16, double: i16, triple: i16) -> CovalentRadii {
    CovalentRadii {
        single: pm(single),
        double: if double >= 0 {
            Some(pm(double as u16))
        } else {
            None
        },
        triple: if triple >= 0 {
            Some(pm(triple as u16))
        } else {
            None
        },
    }
}

const NONE: i16 = -1;

/// Pyykko & Atsumi covalent radii, indexed by atomic number - 1.
/// Values in pm, converted to Bohr at compile time.
#[rustfmt::skip]
static RADII: [CovalentRadii; 118] = [
    //                  single  double  triple
    cr(  32, NONE, NONE),  //   1 H
    cr(  46, NONE, NONE),  //   2 He
    cr( 133,  124, NONE),  //   3 Li
    cr( 102,   90,   85),  //   4 Be
    cr(  85,   78,   73),  //   5 B
    cr(  75,   67,   60),  //   6 C
    cr(  71,   60,   54),  //   7 N
    cr(  63,   57,   53),  //   8 O
    cr(  64,   59,   53),  //   9 F
    cr(  67,   96, NONE),  //  10 Ne
    cr( 155,  160, NONE),  //  11 Na
    cr( 139,  132,  127),  //  12 Mg
    cr( 126,  113,  111),  //  13 Al
    cr( 116,  107,  102),  //  14 Si
    cr( 111,  102,   94),  //  15 P
    cr( 103,   94,   95),  //  16 S
    cr(  99,   95,   93),  //  17 Cl
    cr(  96,  107,   96),  //  18 Ar
    cr( 196,  193, NONE),  //  19 K
    cr( 171,  147,  133),  //  20 Ca
    cr( 148,  116,  114),  //  21 Sc
    cr( 136,  117,  108),  //  22 Ti
    cr( 134,  112,  106),  //  23 V
    cr( 122,  111,  103),  //  24 Cr
    cr( 119,  105,  103),  //  25 Mn
    cr( 116,  109,  102),  //  26 Fe
    cr( 111,  103,   96),  //  27 Co
    cr( 110,  101,  101),  //  28 Ni
    cr( 112,  115,  120),  //  29 Cu
    cr( 118,  120, NONE),  //  30 Zn
    cr( 124,  117,  121),  //  31 Ga
    cr( 121,  111,  114),  //  32 Ge
    cr( 121,  114,  106),  //  33 As
    cr( 116,  107,  107),  //  34 Se
    cr( 114,  109,  110),  //  35 Br
    cr( 117,  121,  108),  //  36 Kr
    cr( 210,  202, NONE),  //  37 Rb
    cr( 185,  157,  139),  //  38 Sr
    cr( 163,  130,  124),  //  39 Y
    cr( 154,  127,  121),  //  40 Zr
    cr( 147,  125,  116),  //  41 Nb
    cr( 138,  121,  113),  //  42 Mo
    cr( 128,  120,  110),  //  43 Tc
    cr( 125,  114,  103),  //  44 Ru
    cr( 125,  110,  106),  //  45 Rh
    cr( 120,  117,  112),  //  46 Pd
    cr( 128,  139,  137),  //  47 Ag
    cr( 136,  144, NONE),  //  48 Cd
    cr( 142,  136,  146),  //  49 In
    cr( 140,  130,  132),  //  50 Sn
    cr( 140,  133,  127),  //  51 Sb
    cr( 136,  128,  121),  //  52 Te
    cr( 133,  129,  125),  //  53 I
    cr( 131,  135,  122),  //  54 Xe
    cr( 232,  209, NONE),  //  55 Cs
    cr( 196,  161,  149),  //  56 Ba
    cr( 180,  139,  139),  //  57 La
    cr( 163,  137,  131),  //  58 Ce
    cr( 176,  138,  128),  //  59 Pr
    cr( 174,  137, NONE),  //  60 Nd
    cr( 173,  135, NONE),  //  61 Pm
    cr( 172,  134, NONE),  //  62 Sm
    cr( 168,  134, NONE),  //  63 Eu
    cr( 169,  135,  132),  //  64 Gd
    cr( 168,  135, NONE),  //  65 Tb
    cr( 167,  133, NONE),  //  66 Dy
    cr( 166,  133, NONE),  //  67 Ho
    cr( 165,  133, NONE),  //  68 Er
    cr( 164,  131, NONE),  //  69 Tm
    cr( 170,  129, NONE),  //  70 Yb
    cr( 162,  131,  131),  //  71 Lu
    cr( 152,  128,  122),  //  72 Hf
    cr( 146,  126,  119),  //  73 Ta
    cr( 137,  120,  115),  //  74 W
    cr( 131,  119,  110),  //  75 Re
    cr( 129,  116,  109),  //  76 Os
    cr( 122,  115,  107),  //  77 Ir
    cr( 123,  112,  110),  //  78 Pt
    cr( 124,  121,  123),  //  79 Au
    cr( 133,  142, NONE),  //  80 Hg
    cr( 144,  142,  150),  //  81 Tl
    cr( 144,  135,  137),  //  82 Pb
    cr( 151,  141,  135),  //  83 Bi
    cr( 145,  135,  129),  //  84 Po
    cr( 147,  138,  138),  //  85 At
    cr( 142,  145,  133),  //  86 Rn
    cr( 223,  218, NONE),  //  87 Fr
    cr( 201,  173,  159),  //  88 Ra
    cr( 186,  153,  140),  //  89 Ac
    cr( 175,  143,  136),  //  90 Th
    cr( 169,  138,  129),  //  91 Pa
    cr( 170,  134,  118),  //  92 U
    cr( 171,  136,  116),  //  93 Np
    cr( 172,  135, NONE),  //  94 Pu
    cr( 166,  135, NONE),  //  95 Am
    cr( 166,  136, NONE),  //  96 Cm
    cr( 168,  139, NONE),  //  97 Bk
    cr( 168,  140, NONE),  //  98 Cf
    cr( 165,  140, NONE),  //  99 Es
    cr( 167, NONE, NONE),  // 100 Fm
    cr( 173,  139, NONE),  // 101 Md
    cr( 176,  159, NONE),  // 102 No
    cr( 161,  141, NONE),  // 103 Lr
    cr( 157,  140,  131),  // 104 Rf
    cr( 149,  136,  126),  // 105 Db
    cr( 143,  128,  121),  // 106 Sg
    cr( 141,  128,  119),  // 107 Bh
    cr( 134,  125,  118),  // 108 Hs
    cr( 129,  125,  113),  // 109 Mt
    cr( 128,  116,  112),  // 110 Ds
    cr( 121,  116,  118),  // 111 Rg
    cr( 122,  137,  130),  // 112 Cn
    cr( 136, NONE, NONE),  // 113 Nh
    cr( 143, NONE, NONE),  // 114 Fl
    cr( 162, NONE, NONE),  // 115 Mc
    cr( 175, NONE, NONE),  // 116 Lv
    cr( 165, NONE, NONE),  // 117 Ts
    cr( 157, NONE, NONE),  // 118 Og
];

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case(Element::H, 32, None, None)]
    #[case(Element::C, 75, Some(67), Some(60))]
    #[case(Element::N, 71, Some(60), Some(54))]
    #[case(Element::O, 63, Some(57), Some(53))]
    #[case(Element::Fe, 116, Some(109), Some(102))]
    fn test_covalent_radii(
        #[case] element: Element,
        #[case] single_pm: u16,
        #[case] double_pm: Option<u16>,
        #[case] triple_pm: Option<u16>,
    ) {
        let r = covalent_radii(element);
        let expected = Length::picometer(single_pm as f64);
        assert!((r.single.as_bohr() - expected.as_bohr()).abs() < 1e-10);
        assert_eq!(
            r.double.map(|v| (v.as_angstrom() * 100.0).round() as u16),
            double_pm,
        );
        assert_eq!(
            r.triple.map(|v| (v.as_angstrom() * 100.0).round() as u16),
            triple_pm,
        );
    }

    #[rstest]
    #[case(Element::C, 1, 75)]
    #[case(Element::C, 2, 67)]
    #[case(Element::C, 3, 60)]
    #[case(Element::H, 2, 32)] // falls back to single
    #[case(Element::H, 3, 32)] // falls back to single
    fn test_covalent_radii_for_order(
        #[case] element: Element,
        #[case] order: u8,
        #[case] expected_pm: u16,
    ) {
        let r = covalent_radii(element);
        let v = r.for_order(order).unwrap();
        let expected = Length::picometer(expected_pm as f64);
        assert!((v.as_bohr() - expected.as_bohr()).abs() < 1e-10);
    }
}
