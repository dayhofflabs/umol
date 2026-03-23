//! Van-Catledge HMO heteroatom parameters derived from PPP calculations.
//!
//! F. A. Van-Catledge, J. Org. Chem. 1980, 45, 4801-4802.
//!
//! Hamiltonian definitions:
//!   alpha_X = alpha_0 + h_X * beta_0
//!   beta_XY = k_XY * beta_0

use umol_data::Element;

pub struct VanCatledgeParams;

impl VanCatledgeParams {
    const N: usize = 13;

    const ATOM_TYPES: [(Element, u8); Self::N] = [
        (Element::C, 1),
        (Element::B, 0),
        (Element::N, 1),
        (Element::N, 2),
        (Element::O, 1),
        (Element::O, 2),
        (Element::F, 2),
        (Element::Si, 1),
        (Element::P, 1),
        (Element::P, 2),
        (Element::S, 1),
        (Element::S, 2),
        (Element::Cl, 2),
    ];

    const H_X: [f64; Self::N] = [
        0.00,  // C
        -0.45, // B
        0.51,  // N1
        1.37,  // N2
        0.97,  // O1
        2.09,  // O2
        2.71,  // F
        0.00,  // Si
        0.19,  // P1
        0.75,  // P2
        0.46,  // S1
        1.11,  // S2
        1.48,  // Cl
    ];

    const FREE_VALENCE_REF: [f64; Self::N] = [
        1.732, // C
        1.705, // B
        1.393, // N1
        1.583, // N2
        0.909, // O1
        0.942, // O2
        0.179, // F
        1.732, // Si
        1.409, // P1
        1.666, // P2
        0.962, // S1
        1.229, // S2
        0.321, // Cl
    ];

    //        C     B     N1    N2    O1    O2    F     Si    P1    P2    S1    S2    Cl
    #[rustfmt::skip]
    const K_XY: [[f64; Self::N]; Self::N] = [
        [1.00, 0.73, 1.02, 0.89, 1.06, 0.66, 0.52, 0.75, 0.77, 0.76, 0.81, 0.69, 0.62], // C
        [0.73, 0.87, 0.66, 0.53, 0.60, 0.35, 0.26, 0.57, 0.53, 0.54, 0.51, 0.44, 0.41], // B
        [1.02, 0.66, 1.09, 0.99, 1.14, 0.80, 0.65, 0.72, 0.78, 0.81, 0.83, 0.78, 0.77], // N1
        [0.89, 0.53, 0.99, 0.98, 1.13, 0.89, 0.77, 0.43, 0.55, 0.64, 0.68, 0.73, 0.80], // N2
        [1.06, 0.60, 1.14, 1.13, 1.26, 1.02, 0.92, 0.65, 0.75, 0.82, 0.84, 0.85, 0.88], // O1
        [0.66, 0.35, 0.80, 0.89, 1.02, 0.95, 0.94, 0.24, 0.31, 0.39, 0.43, 0.54, 0.70], // O2
        [0.52, 0.26, 0.65, 0.77, 0.92, 0.94, 1.04, 0.17, 0.21, 0.22, 0.28, 0.32, 0.51], // F
        [0.75, 0.57, 0.72, 0.43, 0.65, 0.24, 0.17, 0.64, 0.62, 0.52, 0.61, 0.40, 0.34], // Si
        [0.77, 0.53, 0.78, 0.55, 0.75, 0.31, 0.21, 0.62, 0.63, 0.58, 0.65, 0.48, 0.35], // P1
        [0.76, 0.54, 0.81, 0.64, 0.82, 0.39, 0.22, 0.52, 0.58, 0.63, 0.65, 0.60, 0.55], // P2
        [0.81, 0.51, 0.83, 0.68, 0.84, 0.43, 0.28, 0.61, 0.65, 0.65, 0.68, 0.58, 0.52], // S1
        [0.69, 0.44, 0.78, 0.73, 0.85, 0.54, 0.32, 0.40, 0.48, 0.60, 0.58, 0.63, 0.59], // S2
        [0.62, 0.41, 0.77, 0.80, 0.88, 0.70, 0.51, 0.34, 0.35, 0.55, 0.52, 0.59, 0.68], // Cl
    ];

    fn atom_type_index(element: Element, pi_electrons: u8) -> Option<usize> {
        Self::ATOM_TYPES
            .iter()
            .position(|&(e, n)| e == element && n == pi_electrons)
    }

    /// Coulomb integral correction h_X for the given atom type.
    pub fn h_x(element: Element, pi_electrons: u8) -> Option<f64> {
        Self::atom_type_index(element, pi_electrons).map(|i| Self::H_X[i])
    }

    /// Resonance integral correction k_XY for the given pair of atom types.
    pub fn k_xy(x: (Element, u8), y: (Element, u8)) -> Option<f64> {
        let ix = Self::atom_type_index(x.0, x.1)?;
        let iy = Self::atom_type_index(y.0, y.1)?;
        Some(Self::K_XY[ix][iy])
    }

    /// Free valence reference F_X for the given atom type.
    pub fn free_valence_ref(element: Element, pi_electrons: u8) -> Option<f64> {
        Self::atom_type_index(element, pi_electrons).map(|i| Self::FREE_VALENCE_REF[i])
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::carbon(Element::C, 1, Some(0.0))]
    #[case::boron(Element::B, 0, Some(-0.45))]
    #[case::nitrogen_1(Element::N, 1, Some(0.51))]
    #[case::nitrogen_2(Element::N, 2, Some(1.37))]
    #[case::oxygen_1(Element::O, 1, Some(0.97))]
    #[case::sulfur_2(Element::S, 2, Some(1.11))]
    #[case::chlorine(Element::Cl, 2, Some(1.48))]
    #[case::unknown(Element::He, 0, None)]
    fn test_h_x(#[case] element: Element, #[case] pi_electrons: u8, #[case] expected: Option<f64>) {
        assert_eq!(VanCatledgeParams::h_x(element, pi_electrons), expected);
    }

    #[rstest]
    #[case::cc(Element::C, 1, Element::C, 1, Some(1.0))]
    #[case::cn1(Element::C, 1, Element::N, 1, Some(1.02))]
    #[case::co1(Element::C, 1, Element::O, 1, Some(1.06))]
    #[case::n1n2(Element::N, 1, Element::N, 2, Some(0.99))]
    #[case::cl_cl(Element::Cl, 2, Element::Cl, 2, Some(0.68))]
    #[case::unknown(Element::He, 0, Element::C, 1, None)]
    fn test_k_xy(
        #[case] ex: Element,
        #[case] nx: u8,
        #[case] ey: Element,
        #[case] ny: u8,
        #[case] expected: Option<f64>,
    ) {
        assert_eq!(VanCatledgeParams::k_xy((ex, nx), (ey, ny)), expected);
    }

    #[test]
    fn test_k_xy_symmetric() {
        for &(ex, nx) in &VanCatledgeParams::ATOM_TYPES {
            for &(ey, ny) in &VanCatledgeParams::ATOM_TYPES {
                let kxy = VanCatledgeParams::k_xy((ex, nx), (ey, ny)).unwrap();
                let kyx = VanCatledgeParams::k_xy((ey, ny), (ex, nx)).unwrap();
                assert!(
                    (kxy - kyx).abs() < 1e-10,
                    "k_XY not symmetric for ({:?},{}) vs ({:?},{}): {} != {}",
                    ex,
                    nx,
                    ey,
                    ny,
                    kxy,
                    kyx
                );
            }
        }
    }

    #[rstest]
    #[case::carbon(Element::C, 1, Some(1.732))]
    #[case::boron(Element::B, 0, Some(1.705))]
    #[case::fluorine(Element::F, 2, Some(0.179))]
    #[case::unknown(Element::He, 0, None)]
    fn test_free_valence_ref(
        #[case] element: Element,
        #[case] pi_electrons: u8,
        #[case] expected: Option<f64>,
    ) {
        assert_eq!(
            VanCatledgeParams::free_valence_ref(element, pi_electrons),
            expected
        );
    }
}
