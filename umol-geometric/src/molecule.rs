//! Born-Oppenheimer molecular model: N classical nuclei in 3D space.

use nalgebra::{DMatrix, Vector3};
use umol_data::element::Element;
use umol_data::spin::SpinMultiplicity;
use umol_data::units::{Angle, Length};

use crate::coordinates::Coordinates;
use crate::point_group::{PointGroup, C1};

/// 3D molecular geometry under the Born-Oppenheimer approximation.
///
/// Type parameter `G` represents the point group symmetry.
/// Defaults to `C1` (no symmetry) for asymmetric molecules.
/// Coordinates are stored internally in atomic units (Bohr).
pub struct Molecule<G: PointGroup = C1> {
    elements: Vec<Element>,
    coords: Coordinates<G>,
    charge: i32,
    multiplicity: SpinMultiplicity,
}

impl<G: PointGroup> Molecule<G> {
    /// Number of atoms.
    pub fn atom_count(&self) -> usize {
        self.elements.len()
    }

    /// Total number of electrons, derived from elements and charge.
    pub fn num_electrons(&self) -> u32 {
        let nuclear_charge: u32 = self.elements.iter().map(|e| e.atomic_number() as u32).sum();
        (nuclear_charge as i64 - self.charge as i64) as u32
    }

    /// Element of atom `i`.
    pub fn element(&self, i: usize) -> Element {
        self.elements[i]
    }

    /// Molecular charge.
    pub fn charge(&self) -> i32 {
        self.charge
    }

    /// Spin multiplicity.
    pub fn multiplicity(&self) -> SpinMultiplicity {
        self.multiplicity
    }

    /// Cartesian coordinates as a 3×N matrix (Bohr).
    pub fn cartesian_coords(&self) -> &DMatrix<f64> {
        match &self.coords {
            Coordinates::Cartesian(m) => m,
            Coordinates::Symmetric { full_coords, .. } => full_coords,
        }
    }

    /// Euclidean distance between atoms `i` and `j`.
    pub fn distance(&self, i: usize, j: usize) -> Length {
        Length::bohr(self.vec(i, j).norm())
    }

    /// Bond angle i-j-k (vertex at j).
    pub fn angle(&self, i: usize, j: usize, k: usize) -> Angle {
        let u = self.vec(j, i);
        let v = self.vec(j, k);
        Angle::radians(u.angle(&v))
    }

    /// Torsion (dihedral) angle i-j-k-l.
    ///
    /// Defined as the angle between the i-j-k and j-k-l planes,
    /// with sign by IUPAC convention: positive for clockwise rotation
    /// when viewed along j→k.
    pub fn torsion(&self, i: usize, j: usize, k: usize, l: usize) -> Angle {
        let b1 = self.vec(i, j);
        let b2 = self.vec(j, k);
        let b3 = self.vec(k, l);
        let n1 = b1.cross(&b2);
        let n2 = b2.cross(&b3);
        let m = n1.cross(&b2.normalize());
        Angle::radians(f64::atan2(m.dot(&n2), n1.dot(&n2)))
    }

    /// Wilson out-of-plane angle: angle of atom `i` out of the plane
    /// defined by atoms j-k-l (with `j` as the central vertex).
    pub fn out_of_plane(&self, i: usize, j: usize, k: usize, l: usize) -> Angle {
        let e_ji = self.vec(j, i).normalize();
        let e_jk = self.vec(j, k);
        let e_jl = self.vec(j, l);
        let n = e_jk.cross(&e_jl);
        let n_norm = n.norm();
        if n_norm < 1e-14 {
            return Angle::radians(0.0);
        }
        Angle::radians((e_ji.dot(&n) / n_norm).clamp(-1.0, 1.0).asin())
    }

    fn vec(&self, from: usize, to: usize) -> Vector3<f64> {
        let m = self.cartesian_coords();
        Vector3::new(
            m[(0, to)] - m[(0, from)],
            m[(1, to)] - m[(1, from)],
            m[(2, to)] - m[(2, from)],
        )
    }
}

impl Molecule<C1> {
    /// Create a molecule from elements and Cartesian coordinates in Angstroms.
    ///
    /// `coords` is a flat slice of length 3*N: [x0, y0, z0, x1, y1, z1, ...].
    /// Coordinates are converted to Bohr internally.
    pub fn from_cartesian_angstrom(
        elements: Vec<Element>,
        coords: &[f64],
        charge: i32,
        spin: SpinMultiplicity,
    ) -> Self {
        let n = elements.len();
        assert_eq!(coords.len(), 3 * n, "coords must have 3*N elements");
        let matrix = DMatrix::from_fn(3, n, |r, c| Length::angstrom(coords[3 * c + r]).as_bohr());
        Self {
            elements,
            coords: Coordinates::Cartesian(matrix),
            charge,
            multiplicity: spin,
        }
    }

    /// Create a molecule from elements and Cartesian coordinates in Bohr.
    ///
    /// `coords` is a flat slice of length 3*N: [x0, y0, z0, x1, y1, z1, ...].
    pub fn from_cartesian_bohr(
        elements: Vec<Element>,
        coords: &[f64],
        charge: i32,
        spin: SpinMultiplicity,
    ) -> Self {
        let n = elements.len();
        assert_eq!(coords.len(), 3 * n, "coords must have 3*N elements");
        let matrix = DMatrix::from_fn(3, n, |r, c| coords[3 * c + r]);
        Self {
            elements,
            coords: Coordinates::Cartesian(matrix),
            charge,
            multiplicity: spin,
        }
    }
}

#[cfg(test)]
mod tests {
    use float_cmp::approx_eq;
    use rstest::rstest;
    use umol_data::element::Element::*;

    use super::*;

    fn mol(elements: &[Element], coords: &[f64]) -> Molecule {
        Molecule::from_cartesian_angstrom(elements.to_vec(), coords, 0, SpinMultiplicity::Singlet)
    }

    // Water: O-H 0.96 Å, H-O-H 104.5°
    #[rustfmt::skip]
    fn water() -> Molecule {
        mol(&[O, H, H], &[
            0.000,  0.000,  0.000,
            0.960,  0.000,  0.000,
           -0.240,  0.930,  0.000,
        ])
    }

    // Ethane: C-C 1.54 Å, staggered (H-C-C-H torsion ≈ 60°)
    #[rustfmt::skip]
    fn ethane() -> Molecule {
        mol(&[C, C, H, H, H, H, H, H], &[
            0.000,  0.000,  0.000,
            1.540,  0.000,  0.000,
           -0.360, -0.510,  0.890,
           -0.360, -0.510, -0.890,
           -0.360,  1.020,  0.000,
            1.900,  0.510,  0.890,
            1.900,  0.510, -0.890,
            1.900, -1.020,  0.000,
        ])
    }

    // Ethylene: C=C 1.34 Å, planar (all atoms in z=0)
    #[rustfmt::skip]
    fn ethylene() -> Molecule {
        mol(&[C, C, H, H, H, H], &[
            0.000,  0.000,  0.000,
            1.340,  0.000,  0.000,
           -0.540,  0.930,  0.000,
           -0.540, -0.930,  0.000,
            1.880,  0.930,  0.000,
            1.880, -0.930,  0.000,
        ])
    }

    // Acetylene: C≡C 1.20 Å, linear
    #[rustfmt::skip]
    fn acetylene() -> Molecule {
        mol(&[C, C, H, H], &[
            0.000, 0.000, 0.000,
            1.200, 0.000, 0.000,
           -1.060, 0.000, 0.000,
            2.260, 0.000, 0.000,
        ])
    }

    #[rstest]
    #[case::oh_water(water(), 0, 1, 0.96)]
    #[case::cc_ethane(ethane(), 0, 1, 1.54)]
    #[case::cc_ethylene(ethylene(), 0, 1, 1.34)]
    #[case::cc_acetylene(acetylene(), 0, 1, 1.20)]
    fn test_molecule_distance(
        #[case] m: Molecule,
        #[case] i: usize,
        #[case] j: usize,
        #[case] expected_ang: f64,
    ) {
        approx_eq!(
            f64,
            m.distance(i, j).as_angstrom(),
            expected_ang,
            epsilon = 0.001
        );
    }

    #[rstest]
    #[case::hoh_water(water(), 1, 0, 2, 104.5)]
    #[case::hcc_ethane(ethane(), 2, 0, 1, 109.5)]
    #[case::hcc_ethylene(ethylene(), 2, 0, 1, 120.0)]
    #[case::hcc_acetylene(acetylene(), 2, 0, 1, 180.0)]
    fn test_molecule_angle(
        #[case] m: Molecule,
        #[case] i: usize,
        #[case] j: usize,
        #[case] k: usize,
        #[case] expected_deg: f64,
    ) {
        approx_eq!(
            f64,
            m.angle(i, j, k).as_degrees(),
            expected_deg,
            epsilon = 1.0
        );
    }

    #[rstest]
    #[case::staggered_ethane(ethane(), 2, 0, 1, 5, 60.0)]
    #[case::trans_ethylene(ethylene(), 2, 0, 1, 5, 180.0)]
    fn test_molecule_torsion(
        #[case] m: Molecule,
        #[case] i: usize,
        #[case] j: usize,
        #[case] k: usize,
        #[case] l: usize,
        #[case] expected_deg: f64,
    ) {
        approx_eq!(
            f64,
            m.torsion(i, j, k, l).as_degrees().abs(),
            expected_deg,
            epsilon = 2.0
        );
    }

    #[rstest]
    #[case::planar_ethylene(ethylene(), 2, 0, 1, 3, 0.0)]
    fn test_molecule_out_of_plane(
        #[case] m: Molecule,
        #[case] i: usize,
        #[case] j: usize,
        #[case] k: usize,
        #[case] l: usize,
        #[case] expected_deg: f64,
    ) {
        approx_eq!(
            f64,
            m.out_of_plane(i, j, k, l).as_degrees(),
            expected_deg,
            epsilon = 1.0
        );
    }
}
