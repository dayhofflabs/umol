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
    /// Cartesian multinomial x^i y^j z^k centered on an atom.
    Cartesian,
}

#[derive(Debug, Clone)]
pub struct BasisFunction {
    pub atom_index: usize,
    pub kind: BasisKind,
    pub n: i32,
    pub l: i32,
    pub m: i32,
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
