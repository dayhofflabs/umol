//! Symmetry-adapted basis definition and reduction.

use crate::point_group::Irrep;

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
                let shell = (b'f' - 3
                    + l as u8
                    + if l >= 5 { 1 } else { 0 }
                    + if l >= 10 { 1 } else { 0 }
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
