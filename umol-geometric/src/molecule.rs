//! Born-Oppenheimer molecular model: N classical nuclei in 3D space.

use nalgebra::{DMatrix, Vector3};
use umol_data::element::Element;
use umol_data::spin::SpinMultiplicity;
use umol_data::units::{Angle, Length};
use umol_msym::{
    detect_symmetry, symmetrize as symmetrize_centers, EquivalenceSet, Error as SymmetryError,
    PointGroup, SymmetryCenter, SymmetryOp, Thresholds,
};

use crate::coordinates::Coordinates;

/// 3D molecular geometry under the Born-Oppenheimer approximation.
///
/// Every molecule carries point group symmetry data. Defaults to C1 (trivial).
/// Coordinates are stored internally in atomic units (Bohr).
pub struct Molecule {
    elements: Vec<Element>,
    coords: Coordinates,
    charge: i32,
    multiplicity: SpinMultiplicity,

    group: &'static PointGroup,
    /// Atom index orbits under symmetry operations. Each inner vec is one equivalence set.
    equivalence_sets: Vec<Vec<usize>>,
    /// One permutation per symmetry operation: atom_permutations[op][i] = j means
    /// operation `op` maps atom i to atom j.
    atom_permutations: Vec<Vec<usize>>,
}

impl Molecule {
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
        let Coordinates::Cartesian(ref m) = self.coords;
        m
    }

    /// Point group.
    pub fn point_group(&self) -> &'static PointGroup {
        self.group
    }

    /// Equivalence sets (atom index orbits under symmetry operations).
    pub fn equivalence_sets(&self) -> &[Vec<usize>] {
        &self.equivalence_sets
    }

    /// Atom permutations: one per symmetry operation.
    pub fn atom_permutations(&self) -> &[Vec<usize>] {
        &self.atom_permutations
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

    /// Create a molecule from elements and Cartesian coordinates in Angstroms.
    ///
    /// `coords` is a flat slice of length 3*N: [x0, y0, z0, x1, y1, z1, ...].
    /// Coordinates are converted to Bohr internally. Symmetry defaults to C1.
    pub fn from_cartesian_angstrom(
        elements: Vec<Element>,
        coords: &[f64],
        charge: i32,
        spin: SpinMultiplicity,
    ) -> Self {
        let n = elements.len();
        assert_eq!(coords.len(), 3 * n, "coords must have 3*N elements");
        let matrix = DMatrix::from_fn(3, n, |r, c| Length::angstrom(coords[3 * c + r]).as_bohr());
        Self::from_parts(elements, Coordinates::Cartesian(matrix), charge, spin)
    }

    /// Create a molecule from elements and Cartesian coordinates in Bohr.
    ///
    /// `coords` is a flat slice of length 3*N: [x0, y0, z0, x1, y1, z1, ...].
    /// Symmetry defaults to C1.
    pub fn from_cartesian_bohr(
        elements: Vec<Element>,
        coords: &[f64],
        charge: i32,
        spin: SpinMultiplicity,
    ) -> Self {
        let n = elements.len();
        assert_eq!(coords.len(), 3 * n, "coords must have 3*N elements");
        let matrix = DMatrix::from_fn(3, n, |r, c| coords[3 * c + r]);
        Self::from_parts(elements, Coordinates::Cartesian(matrix), charge, spin)
    }

    /// Detect point group symmetry. Returns a new molecule with the discovered
    /// group, equivalence sets, and atom permutations.
    pub fn perceive_symmetry(&self, thresholds: Thresholds) -> Result<Molecule, SymmetryError> {
        let result = detect_symmetry(&self.to_symmetry_centers(), thresholds)?;
        let coords = self.cartesian_coords();
        let eq_sets = equivalence_sets_as_indices(&result.equivalence_sets);
        let atom_permutations =
            compute_atom_permutations(coords, result.group.operations(), &self.elements);

        Ok(Molecule {
            elements: self.elements.clone(),
            coords: Coordinates::Cartesian(coords.clone()),
            charge: self.charge,
            multiplicity: self.multiplicity,
            group: result.group,
            equivalence_sets: eq_sets,
            atom_permutations,
        })
    }

    /// Symmetrize: perceive symmetry and snap coordinates to exact symmetry.
    /// Returns a new molecule with exact symmetry.
    pub fn symmetrize(&self, thresholds: Thresholds) -> Result<Molecule, SymmetryError> {
        let result = symmetrize_centers(&self.to_symmetry_centers(), thresholds)?;

        let n = result.centers.len();
        let matrix = DMatrix::from_fn(3, n, |r, c| {
            Length::angstrom(result.centers[c].position[r]).as_bohr()
        });

        let eq_sets = equivalence_sets_as_indices(&result.equivalence_sets);
        let atom_permutations =
            compute_atom_permutations(&matrix, result.group.operations(), &self.elements);

        Ok(Molecule {
            elements: self.elements.clone(),
            coords: Coordinates::Cartesian(matrix),
            charge: self.charge,
            multiplicity: self.multiplicity,
            group: result.group,
            equivalence_sets: eq_sets,
            atom_permutations,
        })
    }

    /// Convert molecule atoms to libmsym SymmetryCenter format (positions in Angstroms).
    fn to_symmetry_centers(&self) -> Vec<SymmetryCenter> {
        let m = self.cartesian_coords();
        self.elements
            .iter()
            .enumerate()
            .map(|(i, &elem)| SymmetryCenter {
                atomic_number: elem.atomic_number() as i32,
                mass: elem.mass(),
                position: [
                    Length::bohr(m[(0, i)]).as_angstrom(),
                    Length::bohr(m[(1, i)]).as_angstrom(),
                    Length::bohr(m[(2, i)]).as_angstrom(),
                ],
                name: String::new(),
            })
            .collect()
    }

    fn from_parts(
        elements: Vec<Element>,
        coords: Coordinates,
        charge: i32,
        multiplicity: SpinMultiplicity,
    ) -> Self {
        let n = elements.len();
        let equivalence_sets: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
        let atom_permutations = vec![(0..n).collect()];
        Self {
            elements,
            coords,
            charge,
            multiplicity,
            group: PointGroup::c1(),
            equivalence_sets,
            atom_permutations,
        }
    }
}

/// Convert EquivalenceSets (which contain full SymmetryCenter data) to atom
/// index orbits by matching positions back to the result element list.
fn equivalence_sets_as_indices(sets: &[EquivalenceSet]) -> Vec<Vec<usize>> {
    // EquivalenceSet.centers are ordered the same as the context center list,
    // so indices are sequential across sets.
    let mut idx = 0;
    sets.iter()
        .map(|es| {
            let indices: Vec<usize> = (idx..idx + es.centers.len()).collect();
            idx += es.centers.len();
            indices
        })
        .collect()
}

/// Compute atom permutations for each symmetry operation.
///
/// For operation R with 3×3 matrix M, atom i at position p is mapped to
/// the atom j whose position is closest to M·p.
fn compute_atom_permutations(
    coords: &DMatrix<f64>,
    operations: &[SymmetryOp],
    elements: &[Element],
) -> Vec<Vec<usize>> {
    let n = coords.ncols();
    operations
        .iter()
        .map(|op| {
            let mut perm = vec![0usize; n];
            for i in 0..n {
                let p = Vector3::new(coords[(0, i)], coords[(1, i)], coords[(2, i)]);
                let rp = op.matrix * p;

                // Find nearest atom of same element
                let mut best_j = 0;
                let mut best_d2 = f64::MAX;
                for j in 0..n {
                    if elements[j] != elements[i] {
                        continue;
                    }
                    let q = Vector3::new(coords[(0, j)], coords[(1, j)], coords[(2, j)]);
                    let d2 = (rp - q).norm_squared();
                    if d2 < best_d2 {
                        best_d2 = d2;
                        best_j = j;
                    }
                }
                perm[i] = best_j;
            }
            perm
        })
        .collect()
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

    // Water: symmetric geometry for perception tests (O at origin, H symmetric about xz)
    #[rustfmt::skip]
    fn symmetric_water() -> Molecule {
        mol(&[O, H, H], &[
            0.000,  0.000,  0.117,
            0.000,  0.757, -0.469,
            0.000, -0.757, -0.469,
        ])
    }

    // Methane: Td geometry for perception tests
    #[rustfmt::skip]
    fn symmetric_methane() -> Molecule {
        mol(&[C, H, H, H, H], &[
            0.000,  0.000,  0.000,
            0.629,  0.629,  0.629,
           -0.629, -0.629,  0.629,
           -0.629,  0.629, -0.629,
            0.629, -0.629, -0.629,
        ])
    }

    #[rstest]
    #[case(symmetric_water(), "C2v", 4, 2)]
    #[case(symmetric_methane(), "Td", 24, 2)]
    fn test_molecule_perceive_symmetry(
        #[case] m: Molecule,
        #[case] expected_group: &str,
        #[case] expected_ops: usize,
        #[case] expected_eq_sets: usize,
    ) {
        let thresholds = Thresholds::default();
        let sym = m.perceive_symmetry(thresholds).unwrap();
        assert_eq!(sym.point_group().to_string(), expected_group);
        assert_eq!(sym.point_group().order(), expected_ops);
        assert_eq!(sym.equivalence_sets().len(), expected_eq_sets);
        assert_eq!(sym.atom_permutations().len(), expected_ops);

        // Each permutation must be a valid bijection
        let n = sym.atom_count();
        for perm in sym.atom_permutations() {
            assert_eq!(perm.len(), n);
            let mut seen = vec![false; n];
            for &j in perm {
                assert!(j < n);
                seen[j] = true;
            }
            assert!(seen.iter().all(|&s| s), "permutation is not a bijection");
        }
    }

    #[rstest]
    fn test_molecule_symmetrize() {
        let thresholds = Thresholds::default();
        let m = symmetric_water();
        let sym = m.symmetrize(thresholds).unwrap();
        assert_eq!(sym.point_group().to_string(), "C2v");
        // Symmetrized coordinates should be more symmetric than input
        let coords = sym.cartesian_coords();
        // y-coordinates of the two H atoms should be exactly opposite
        let y_h1 = coords[(1, 1)];
        let y_h2 = coords[(1, 2)];
        assert!(
            (y_h1 + y_h2).abs() < 1e-10,
            "H atoms not symmetric: {y_h1} vs {y_h2}"
        );
    }
}
