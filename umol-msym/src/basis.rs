//! Symmetry-adapted basis definition and reduction.

use umol_msym_sys as ffi;

use crate::irrep::Irrep;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CartesianAxis {
    X,
    Y,
    Z,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BasisKind {
    /// s-type: atoms themselves (permutation representation). l=0.
    Atom,
    /// p-type: Cartesian displacement on an atom. l=1.
    Displacement(CartesianAxis),
    /// Real spherical harmonic Y_l^m centered on an atom.
    RealSphericalHarmonic,
    /// Cartesian harmonic (symmetric multinomials x^i y^j z^k) centered on an atom.
    CartesianHarmonic,
}

#[derive(Debug, Clone)]
pub struct BasisFunction {
    pub atom_index: usize,
    pub kind: BasisKind,
    /// Index among shells of the same l on the same atom (0 = first, 1 = second, ...).
    /// Translated to the principal quantum number n = l + 1 + shell_index at the FFI boundary.
    pub shell_index: u32,
    pub l: i32,
    pub m: i32,
}

impl BasisFunction {
    /// All 2l+1 real spherical harmonic components of a single shell.
    pub fn shell(atom_index: usize, shell_index: u32, l: i32) -> Vec<BasisFunction> {
        (-l..=l)
            .map(|m| BasisFunction {
                atom_index,
                kind: BasisKind::RealSphericalHarmonic,
                shell_index,
                l,
                m,
            })
            .collect()
    }

    pub(crate) fn from_ffi(
        bf: &ffi::msym_basis_function_t,
        elements_base: *const ffi::msym_element_t,
    ) -> Self {
        let atom_index = unsafe { (bf.element as *const ffi::msym_element_t).offset_from(elements_base) as usize };
        let rsh = unsafe { bf.f.rsh };
        let kind = if bf.type_ == ffi::MSYM_BASIS_TYPE_CARTESIAN {
            BasisKind::CartesianHarmonic
        } else if rsh.l == 1 {
            let axis = match rsh.m {
                1 => CartesianAxis::X,
                -1 => CartesianAxis::Y,
                _ => CartesianAxis::Z,
            };
            BasisKind::Displacement(axis)
        } else {
            BasisKind::RealSphericalHarmonic
        };
        let shell_index = (rsh.n - rsh.l - 1) as u32;
        BasisFunction {
            atom_index,
            kind,
            shell_index,
            l: rsh.l,
            m: rsh.m,
        }
    }

    /// Principal quantum number for libmsym FFI (n = l + 1 + shell_index).
    pub(crate) fn ffi_n(&self) -> i32 {
        self.l + 1 + self.shell_index as i32
    }

    /// Fixed-width orbital name in the format libmsym expects for
    /// `msym_basis_function_t.name` (e.g. "1s", "2px", "3d1+").
    pub(crate) fn libmsym_name_bytes(&self) -> [i8; 8] {
        let n = self.ffi_n();
        let s = match self.l {
            0 => format!("{n}s"),
            1 => {
                let axis = match self.m {
                    1 => "x",
                    -1 => "y",
                    0 => "z",
                    _ => "?",
                };
                format!("{n}p{axis}")
            }
            l => {
                let shell = (b'd' - 2
                    + l as u8
                    // skip 'e'
                    + if l >= 3 { 1 } else { 0 }
                    // skip 'j'
                    + if l >= 7 { 1 } else { 0 }
                    // skip 'o'
                    + if l >= 11 { 1 } else { 0 }
                    // skip 'q'
                    + if l >= 12 { 1 } else { 0 }) as char;
                let sign = if self.m > 0 {
                    "+"
                } else if self.m < 0 {
                    "-"
                } else {
                    ""
                };
                format!("{n}{shell}{}{sign}", self.m.unsigned_abs())
            }
        };
        let mut name = [0i8; 8];
        for (i, &b) in s.as_bytes().iter().take(7).enumerate() {
            name[i] = b as i8;
        }
        name
    }
}

/// One symmetry-adapted linear combination, sparse in the basis function indices.
#[derive(Debug, Clone)]
pub struct Salc {
    pub coefficients: Vec<(usize, f64)>,
}

/// SALCs belonging to a single irrep.
#[derive(Debug, Clone)]
pub struct IrrepBasis {
    pub irrep: Irrep,
    pub salcs: Vec<Salc>,
}

/// Full symmetry-adapted basis across all irreps.
#[derive(Debug, Clone)]
pub struct SalcBasis {
    pub basis_functions: Vec<BasisFunction>,
    pub irreps: Vec<IrrepBasis>,
}

impl SalcBasis {
    pub fn irrep_basis(&self, irrep: Irrep) -> Option<&IrrepBasis> {
        self.irreps.iter().find(|ib| ib.irrep == irrep)
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;

    use super::*;

    #[rstest]
    #[case::s(0, 0, "1s")]
    #[case::px(1, 1, "2px")]
    #[case::py(1, -1, "2py")]
    #[case::pz(1, 0, "2pz")]
    #[case::d_m0(2, 0, "3d0")]
    #[case::d_m2_plus(2, 2, "3d2+")]
    #[case::d_m2_minus(2, -2, "3d2-")]
    #[case::f_m0(3, 0, "4f0")]
    #[case::f_m3_plus(3, 3, "4f3+")]
    #[case::f_m3_minus(3, -3, "4f3-")]
    #[case::g(4, 0, "5g0")]
    #[case::h(5, 0, "6h0")]
    #[case::i(6, 0, "7i0")]
    #[case::k(7, 0, "8k0")]
    #[case::l(8, 0, "9l0")]
    #[case::m(9, 0, "10m0")]
    #[case::n(10, 0, "11n0")]
    #[case::p(11, 0, "12p0")]
    #[case::r(12, 0, "13r0")]
    fn test_basis_function_libmsym_name(#[case] l: i32, #[case] m: i32, #[case] expected: &str) {
        let bf = BasisFunction {
            atom_index: 0,
            kind: BasisKind::RealSphericalHarmonic,
            shell_index: 0,
            l,
            m,
        };
        let name = bf.libmsym_name_bytes();
        let s: String = name.iter().take_while(|&&b| b != 0).map(|&b| b as u8 as char).collect();
        assert_eq!(s, expected);
    }
}
