//! Born-Oppenheimer molecular model: N classical nuclei in 3D space.

use std::collections::HashMap;

use nalgebra::{DMatrix, DVector, Vector3};
use umol_msym::{
    compute_salcs as compute_salcs_raw, detect_symmetry,
    generate_symmetry_images as generate_image_centers, group,
    lower_symmetry as lower_symmetry_centers, symmetrize as symmetrize_centers, BasisFunction,
    BasisKind, CartesianAxis, EquivalenceSet, Irrep, MatrixRep, MsymError, PointGroup, SalcBasis,
    SchoenfliesSymbol, SymmetryCenter, Thresholds,
};
use umol_shared::element::Element;
use umol_shared::spin::SpinMultiplicity;
use umol_shared::units::angle::Angle;
use umol_shared::units::length::Length;

use crate::coordinates::Coordinates;

/// Numerical zero for cross-product norms and degenerate geometry guards.
const NUMERICAL_ZERO: f64 = 1e-14;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoordinateKind {
    Translation,
    Rotation,
    Vibration,
}

pub struct SymmetryCoordinate {
    pub irrep: Irrep,
    pub kind: CoordinateKind,
    /// Per-atom displacement vectors: `atom_vectors[i]` = `[dx, dy, dz]` for atom `i`.
    pub atom_vectors: Vec<[f64; 3]>,
}

pub struct SymmetryCoordinates {
    pub gamma_total: Vec<(Irrep, u32)>,
    pub gamma_trans: Vec<(Irrep, u32)>,
    pub gamma_rot: Vec<(Irrep, u32)>,
    pub gamma_vib: Vec<(Irrep, u32)>,
    pub coordinates: Vec<SymmetryCoordinate>,
}

pub struct SymmetryDescentResult {
    pub molecule: Molecule,
    pub parent_group: &'static PointGroup,
    pub transform: nalgebra::Matrix3<f64>,
}

/// 3D molecular geometry under the Born-Oppenheimer approximation.
///
/// Every molecule carries point group symmetry data. Defaults to C1 (trivial).
/// Coordinates are stored internally in atomic units (Bohr).
pub struct Molecule {
    elements: Vec<Element>,
    coordinates: Coordinates,
    charge: i32,
    multiplicity: SpinMultiplicity,

    group: &'static PointGroup,
    /// Matrix realization of the point group in the molecule's coordinate frame.
    representation: MatrixRep,
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
    pub fn electron_count(&self) -> u32 {
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
    pub fn cartesian_coordinates(&self) -> &DMatrix<f64> {
        let Coordinates::Cartesian(ref m) = self.coordinates;
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
        if n_norm < NUMERICAL_ZERO {
            return Angle::radians(0.0);
        }
        Angle::radians((e_ji.dot(&n) / n_norm).clamp(-1.0, 1.0).asin())
    }

    fn vec(&self, from: usize, to: usize) -> Vector3<f64> {
        let m = self.cartesian_coordinates();
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
    pub fn perceive_symmetry(&self, thresholds: Thresholds) -> Result<Molecule, MsymError> {
        let result = detect_symmetry(&self.to_symmetry_centers(), thresholds)?;
        let coords = self.cartesian_coordinates();
        let eq_sets = equivalence_sets_as_indices(&result.equivalence_sets);
        let atom_permutations =
            compute_atom_permutations(coords, &result.representation, &self.elements);

        Ok(Molecule {
            elements: self.elements.clone(),
            coordinates: Coordinates::Cartesian(coords.clone()),
            charge: self.charge,
            multiplicity: self.multiplicity,
            group: result.group,
            representation: result.representation,
            equivalence_sets: eq_sets,
            atom_permutations,
        })
    }

    /// Symmetrize: perceive symmetry and snap coordinates to exact symmetry.
    /// Returns a new molecule with exact symmetry.
    pub fn symmetrize(&self, thresholds: Thresholds) -> Result<Molecule, MsymError> {
        let result = symmetrize_centers(&self.to_symmetry_centers(), thresholds)?;

        let n = result.centers.len();
        let matrix = DMatrix::from_fn(3, n, |r, c| {
            Length::angstrom(result.centers[c].position[r]).as_bohr()
        });

        let eq_sets = equivalence_sets_as_indices(&result.equivalence_sets);
        let atom_permutations =
            compute_atom_permutations(&matrix, &result.representation, &self.elements);

        Ok(Molecule {
            elements: self.elements.clone(),
            coordinates: Coordinates::Cartesian(matrix),
            charge: self.charge,
            multiplicity: self.multiplicity,
            group: result.group,
            representation: result.representation,
            equivalence_sets: eq_sets,
            atom_permutations,
        })
    }

    /// Lower the symmetry to a specified subgroup.
    /// Returns the molecule re-perceived under the child group, plus the
    /// parent-to-child orientation transform and correlation table.
    pub fn lower_symmetry(
        &self,
        target: SchoenfliesSymbol,
        thresholds: Thresholds,
    ) -> Result<SymmetryDescentResult, MsymError> {
        let result = lower_symmetry_centers(&self.to_symmetry_centers(), target, thresholds)?;
        let coords = self.cartesian_coordinates();
        let eq_sets = equivalence_sets_as_indices(&result.equivalence_sets);
        let atom_permutations =
            compute_atom_permutations(coords, &result.child_representation, &self.elements);

        Ok(SymmetryDescentResult {
            molecule: Molecule {
                elements: self.elements.clone(),
                coordinates: Coordinates::Cartesian(coords.clone()),
                charge: self.charge,
                multiplicity: self.multiplicity,
                group: result.child_group,
                representation: result.child_representation,
                equivalence_sets: eq_sets,
                atom_permutations,
            },
            parent_group: result.parent_group,
            transform: result.transform,
        })
    }

    /// Generate a full molecule from an asymmetric unit and a target point group.
    ///
    /// Positions are in Angstroms, relative to the molecular center (origin).
    pub fn generate_symmetry_images(
        label: SchoenfliesSymbol,
        elements: &[Element],
        positions_angstrom: &[[f64; 3]],
        thresholds: Thresholds,
    ) -> Result<Molecule, MsymError> {
        let centers: Vec<SymmetryCenter> = elements
            .iter()
            .zip(positions_angstrom)
            .map(|(&elem, pos)| SymmetryCenter {
                atomic_number: elem.atomic_number() as i32,
                mass: elem.mass(),
                position: Vector3::from(*pos),
                name: String::new(),
            })
            .collect();

        let result = generate_image_centers(label, &centers, thresholds)?;

        let gen_elements: Vec<Element> = result
            .centers
            .iter()
            .map(|c| Element::from_atomic_number(c.atomic_number as u8).unwrap())
            .collect();

        let n = result.centers.len();
        let matrix = DMatrix::from_fn(3, n, |r, c| {
            Length::angstrom(result.centers[c].position[r]).as_bohr()
        });

        let eq_sets = equivalence_sets_as_indices(&result.equivalence_sets);
        let atom_permutations =
            compute_atom_permutations(&matrix, &result.representation, &gen_elements);

        Ok(Molecule {
            elements: gen_elements,
            coordinates: Coordinates::Cartesian(matrix),
            charge: 0,
            multiplicity: SpinMultiplicity::Singlet,
            group: result.group,
            representation: result.representation,
            equivalence_sets: eq_sets,
            atom_permutations,
        })
    }

    /// Compute symmetry coordinates: decompose 3N Cartesian degrees of freedom
    /// into symmetry-classified translation, rotation, and vibration coordinates.
    ///
    /// Returns projected coordinate vectors grouped by irrep and classified by kind.
    /// Symmetry must have been perceived first.
    pub fn symmetry_coordinates(&self, thresholds: Thresholds) -> SymmetryCoordinates {
        let group = self.point_group();
        let n = self.atom_count();
        let dim3n = 3 * n;

        if group.is_linear() {
            return self.symmetry_coordinates_linear(thresholds);
        }

        let rep = &self.representation;
        let perms = self.atom_permutations();
        let h = group.order();
        assert_eq!(rep.order(), h);
        assert_eq!(perms.len(), h);

        // Step 1: Γ_3N characters (one per class)
        let n_classes = group.class_sizes().len();
        let mut gamma_3n = vec![0.0; n_classes];
        // Use first operation of each class
        let mut class_seen = vec![false; n_classes];
        for (k, op) in group.ops().into_iter().enumerate() {
            let c = op.class();
            if class_seen[c] {
                continue;
            }
            class_seen[c] = true;
            let n_fixed: usize = perms[k]
                .iter()
                .enumerate()
                .filter(|(i, &j)| *i == j)
                .count();
            gamma_3n[c] = n_fixed as f64 * rep.matrix(op).trace();
        }

        // Step 2: reduce
        let gamma_total = group
            .reduce(&gamma_3n)
            .expect("valid 3N representation characters");
        let gamma_trans = group.translation_irreps();
        let gamma_rot = group.rotation_irreps();
        let gamma_vib = subtract_irrep_reps(&gamma_total, &gamma_trans, &gamma_rot);

        // Step 3: build D_3N(R) matrices for each operation
        let d3n_mats: Vec<DMatrix<f64>> = group
            .ops()
            .into_iter()
            .zip(perms.iter())
            .map(|(op, perm)| {
                let mut d = DMatrix::zeros(dim3n, dim3n);
                let matrix = rep.matrix(op);
                for i in 0..n {
                    let j = perm[i];
                    // block (3*j, 3*i) = M_R
                    for r in 0..3 {
                        for c in 0..3 {
                            d[(3 * j + r, 3 * i + c)] = matrix[(r, c)];
                        }
                    }
                }
                d
            })
            .collect();

        // Step 4: project for each irrep
        let irreps = group.irreps();
        let mut coordinates = Vec::new();

        for irrep in &irreps {
            let dim = irrep.dimension() as f64;

            // P_μ = (l_μ / h) Σ_R χ_μ(R) · D_3N(R)
            let mut proj = DMatrix::zeros(dim3n, dim3n);
            for (k, op) in group.ops().into_iter().enumerate() {
                let chi = op.character(*irrep);
                proj += &d3n_mats[k] * (dim * chi / h as f64);
            }

            // SVD to extract nonzero columns
            let svd = proj.svd(true, false);
            let u = svd.u.unwrap();
            let vecs: Vec<DVector<f64>> = (0..dim3n)
                .filter(|&i| svd.singular_values[i] > thresholds.projection)
                .map(|i| u.column(i).into_owned())
                .collect();

            for v in vecs {
                coordinates.push(SymmetryCoordinate {
                    irrep: *irrep,
                    kind: CoordinateKind::Vibration, // classified below
                    atom_vectors: flat_to_atom_vectors(v.as_slice()),
                });
            }
        }

        // Step 5: classify as trans/rot/vib
        classify_coordinates(
            &mut coordinates,
            self.cartesian_coordinates(),
            n,
            &thresholds,
        );

        SymmetryCoordinates {
            gamma_total,
            gamma_trans,
            gamma_rot,
            gamma_vib,
            coordinates,
        }
    }

    fn symmetry_coordinates_linear(&self, thresholds: Thresholds) -> SymmetryCoordinates {
        let n = self.atom_count();
        let dim3n = 3 * n;
        let centers = self.to_symmetry_centers();

        let basis = self.displacement_basis();

        let salc_basis = compute_salcs_raw(&centers, &basis, thresholds)
            .expect("linear SALC computation failed");

        // Convert SALCs to dense coordinate vectors
        let mut coordinates = Vec::new();
        for ib in &salc_basis.irreps {
            for salc in &ib.salcs {
                let mut flat = vec![0.0; dim3n];
                for &(j, c) in &salc.coefficients {
                    flat[j] = c;
                }
                coordinates.push(SymmetryCoordinate {
                    irrep: ib.irrep,
                    kind: CoordinateKind::Vibration,
                    atom_vectors: flat_to_atom_vectors(&flat),
                });
            }
        }

        // Classify
        classify_coordinates(
            &mut coordinates,
            self.cartesian_coordinates(),
            n,
            &thresholds,
        );

        // Build Γ decompositions from the classified coordinates
        let gamma_total = count_irrep_reps(&coordinates, |_| true);
        let gamma_trans = count_irrep_reps(&coordinates, |c| {
            matches!(c.kind, CoordinateKind::Translation)
        });
        let gamma_rot =
            count_irrep_reps(&coordinates, |c| matches!(c.kind, CoordinateKind::Rotation));
        let gamma_vib = count_irrep_reps(&coordinates, |c| {
            matches!(c.kind, CoordinateKind::Vibration)
        });

        SymmetryCoordinates {
            gamma_total,
            gamma_trans,
            gamma_rot,
            gamma_vib,
            coordinates,
        }
    }

    /// Permutation basis (N functions, one per atom, l=0).
    ///
    /// Suitable for reducing the atomic permutation representation into irreps
    /// via `salc_basis`.
    pub fn permutation_basis(&self) -> Vec<BasisFunction> {
        (0..self.atom_count())
            .map(|i| BasisFunction {
                atom_index: i,
                kind: BasisKind::Atom,
                shell_index: 0,
                l: 0,
                m: 0,
            })
            .collect()
    }

    /// Cartesian displacement basis (3N functions) for this molecule.
    ///
    /// Produces x, y, z displacement functions on each atom. These are the l=1
    /// basis functions that decompose into translation, rotation, and vibration.
    pub fn displacement_basis(&self) -> Vec<BasisFunction> {
        let axes = [
            (CartesianAxis::X, 1),
            (CartesianAxis::Y, -1),
            (CartesianAxis::Z, 0),
        ];
        (0..self.atom_count())
            .flat_map(|i| {
                axes.iter().map(move |&(axis, m)| BasisFunction {
                    atom_index: i,
                    kind: BasisKind::Displacement(axis),
                    shell_index: 0,
                    l: 1,
                    m,
                })
            })
            .collect()
    }

    /// AO basis from per-element shell specifications.
    ///
    /// Each entry in `shell_spec` maps an element to its shells as `(l, count)` pairs.
    /// Example: `mol.ao_basis(&[(O, &[(0, 2), (1, 1)]), (H, &[(0, 1)])])`
    pub fn ao_basis(&self, shell_spec: &[(Element, &[(u32, u32)])]) -> Vec<BasisFunction> {
        self.elements
            .iter()
            .enumerate()
            .flat_map(|(atom_index, elem)| {
                let shells = shell_spec
                    .iter()
                    .find(|(e, _)| e == elem)
                    .unwrap_or_else(|| panic!("no shell spec for {elem}"))
                    .1;
                let mut shell_index_by_l = HashMap::<u32, u32>::new();
                shells.iter().flat_map(move |&(l, count)| {
                    let base = *shell_index_by_l.entry(l).or_insert(0);
                    shell_index_by_l.insert(l, base + count);
                    (0..count)
                        .flat_map(move |i| BasisFunction::shell(atom_index, base + i, l as i32))
                })
            })
            .collect()
    }

    /// Compute symmetry-adapted linear combinations (SALCs) for a set of basis functions.
    ///
    /// Each `BasisFunction::atom_index` must be a valid index into this molecule's atoms.
    /// Symmetry must have been perceived first (via `perceive_symmetry` or `symmetrize`).
    pub fn salc_basis(
        &self,
        basis: &[BasisFunction],
        thresholds: Thresholds,
    ) -> Result<SalcBasis, MsymError> {
        compute_salcs_raw(&self.to_symmetry_centers(), basis, thresholds)
    }

    /// Convert molecule atoms to libmsym SymmetryCenter format (positions in Angstroms).
    fn to_symmetry_centers(&self) -> Vec<SymmetryCenter> {
        let m = self.cartesian_coordinates();
        self.elements
            .iter()
            .enumerate()
            .map(|(i, &elem)| SymmetryCenter {
                atomic_number: elem.atomic_number() as i32,
                mass: elem.mass(),
                position: Vector3::new(
                    Length::bohr(m[(0, i)]).as_angstrom(),
                    Length::bohr(m[(1, i)]).as_angstrom(),
                    Length::bohr(m[(2, i)]).as_angstrom(),
                ),
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
        let group = group!(C1);
        let equivalence_sets: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();
        let atom_permutations = vec![(0..n).collect()];
        Self {
            elements,
            coordinates: coords,
            charge,
            multiplicity,
            group,
            representation: MatrixRep::identity_only(group),
            equivalence_sets,
            atom_permutations,
        }
    }
}

fn flat_to_atom_vectors(flat: &[f64]) -> Vec<[f64; 3]> {
    flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect()
}

fn atom_vectors_to_flat(vecs: &[[f64; 3]]) -> Vec<f64> {
    vecs.iter().flat_map(|v| v.iter().copied()).collect()
}

/// Subtract trans and rot irrep multiplicities from total to get vibrational.
fn subtract_irrep_reps(
    total: &[(Irrep, u32)],
    trans: &[(Irrep, u32)],
    rot: &[(Irrep, u32)],
) -> Vec<(Irrep, u32)> {
    let mut result = total.to_vec();
    for sub in [trans, rot] {
        for &(irrep, count) in sub {
            if let Some(entry) = result.iter_mut().find(|(ir, _)| *ir == irrep) {
                entry.1 = entry.1.saturating_sub(count);
            }
        }
    }
    result.retain(|&(_, n)| n > 0);
    result
}

/// Count irrep multiplicities from coordinates matching a predicate.
fn count_irrep_reps(
    coords: &[SymmetryCoordinate],
    pred: impl Fn(&SymmetryCoordinate) -> bool,
) -> Vec<(Irrep, u32)> {
    let mut result: Vec<(Irrep, u32)> = Vec::new();
    for c in coords {
        if !pred(c) {
            continue;
        }
        if let Some(entry) = result.iter_mut().find(|(ir, _)| *ir == c.irrep) {
            entry.1 += 1;
        } else {
            result.push((c.irrep, 1));
        }
    }
    result
}

/// Build translation trial vectors (normalized): uniform displacement along x, y, z.
fn translation_trial_vectors(n_atoms: usize) -> Vec<DVector<f64>> {
    let dim3n = 3 * n_atoms;
    let inv_sqrt_n = 1.0 / (n_atoms as f64).sqrt();
    (0..3)
        .map(|axis| {
            let mut v = DVector::zeros(dim3n);
            for i in 0..n_atoms {
                v[3 * i + axis] = inv_sqrt_n;
            }
            v
        })
        .collect()
}

/// Build rotation trial vectors (normalized): r_i × e_axis for each atom.
fn rotation_trial_vectors(
    coords_bohr: &DMatrix<f64>,
    n_atoms: usize,
    zero: f64,
) -> Vec<DVector<f64>> {
    let dim3n = 3 * n_atoms;
    (0..3)
        .filter_map(|axis| {
            let mut v = DVector::zeros(dim3n);
            for i in 0..n_atoms {
                let x = coords_bohr[(0, i)];
                let y = coords_bohr[(1, i)];
                let z = coords_bohr[(2, i)];
                let (dx, dy, dz) = match axis {
                    0 => (0.0, z, -y),
                    1 => (-z, 0.0, x),
                    _ => (y, -x, 0.0),
                };
                v[3 * i] = dx;
                v[3 * i + 1] = dy;
                v[3 * i + 2] = dz;
            }
            let norm = v.norm();
            if norm > zero {
                v /= norm;
                Some(v)
            } else {
                None
            }
        })
        .collect()
}

/// Classify symmetry coordinates as translation, rotation, or vibration.
///
/// For each irrep subspace, projects translation and rotation trial vectors
/// into the subspace, then rotates the basis to separate them from vibrations.
fn classify_coordinates(
    coordinates: &mut Vec<SymmetryCoordinate>,
    coords_bohr: &DMatrix<f64>,
    n_atoms: usize,
    thresholds: &Thresholds,
) {
    let trans_trials = translation_trial_vectors(n_atoms);
    let rot_trials = rotation_trial_vectors(coords_bohr, n_atoms, thresholds.zero);

    // Group coordinate indices by irrep
    let mut irrep_groups: Vec<(Irrep, Vec<usize>)> = Vec::new();
    for (idx, c) in coordinates.iter().enumerate() {
        if let Some(entry) = irrep_groups.iter_mut().find(|(ir, _)| *ir == c.irrep) {
            entry.1.push(idx);
        } else {
            irrep_groups.push((c.irrep, vec![idx]));
        }
    }

    let tol = thresholds.orthogonalization;

    // For each irrep subspace, extract trans/rot components by projection
    let mut new_coords: Vec<SymmetryCoordinate> = Vec::with_capacity(coordinates.len());

    for (irrep, indices) in &irrep_groups {
        let k = indices.len();
        if k == 0 {
            continue;
        }

        // Build subspace matrix U: columns are the coordinate vectors (flattened)
        let dim3n = 3 * n_atoms;
        let mut u = DMatrix::zeros(dim3n, k);
        for (col, &idx) in indices.iter().enumerate() {
            let flat = atom_vectors_to_flat(&coordinates[idx].atom_vectors);
            for row in 0..dim3n {
                u[(row, col)] = flat[row];
            }
        }

        // Collect orthonormal basis vectors classified as trans, rot, or vib
        let mut trans_basis: Vec<DVector<f64>> = Vec::new();
        let mut rot_basis: Vec<DVector<f64>> = Vec::new();

        // Project translation trial vectors into this subspace
        for t in &trans_trials {
            // p = U * U^T * t (projection into subspace)
            let coeffs = u.transpose() * t;
            let proj = &u * &coeffs;
            let norm = proj.norm();
            if norm > tol {
                let v = proj / norm;
                trans_basis.push(v);
            }
        }

        // Orthogonalize trans_basis
        gram_schmidt(&mut trans_basis, tol);

        // Project rotation trial vectors, then remove trans components
        for r in &rot_trials {
            let coeffs = u.transpose() * r;
            let mut proj = &u * &coeffs;
            // Remove translation components
            for tb in &trans_basis {
                let d = proj.dot(tb);
                proj -= d * tb;
            }
            let norm = proj.norm();
            if norm > tol {
                let v = proj / norm;
                rot_basis.push(v);
            }
        }

        // Orthogonalize rot_basis
        gram_schmidt(&mut rot_basis, tol);

        // The rest of the subspace is vibration: orthogonal complement
        let mut vib_basis: Vec<DVector<f64>> = Vec::new();
        for col in 0..k {
            let mut v = u.column(col).into_owned();
            for tb in &trans_basis {
                let d = v.dot(tb);
                v -= d * tb;
            }
            for rb in &rot_basis {
                let d = v.dot(rb);
                v -= d * rb;
            }
            let norm = v.norm();
            if norm > tol {
                v /= norm;
                vib_basis.push(v);
            }
        }
        gram_schmidt(&mut vib_basis, tol);

        for v in trans_basis {
            new_coords.push(SymmetryCoordinate {
                irrep: *irrep,
                kind: CoordinateKind::Translation,
                atom_vectors: flat_to_atom_vectors(v.as_slice()),
            });
        }
        for v in rot_basis {
            new_coords.push(SymmetryCoordinate {
                irrep: *irrep,
                kind: CoordinateKind::Rotation,
                atom_vectors: flat_to_atom_vectors(v.as_slice()),
            });
        }
        for v in vib_basis {
            new_coords.push(SymmetryCoordinate {
                irrep: *irrep,
                kind: CoordinateKind::Vibration,
                atom_vectors: flat_to_atom_vectors(v.as_slice()),
            });
        }
    }

    *coordinates = new_coords;
}

// TODO: Move into shared algorithm crate
fn gram_schmidt(vecs: &mut Vec<DVector<f64>>, tol: f64) {
    let mut i = 0;
    while i < vecs.len() {
        for j in 0..i {
            let d = vecs[i].dot(&vecs[j]);
            let vj = vecs[j].clone();
            vecs[i] -= d * &vj;
        }
        let norm = vecs[i].norm();
        if norm < tol {
            vecs.remove(i);
        } else {
            vecs[i] /= norm;
            i += 1;
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
    representation: &MatrixRep,
    elements: &[Element],
) -> Vec<Vec<usize>> {
    let n = coords.ncols();
    representation
        .matrices()
        .iter()
        .map(|matrix| {
            let mut perm = vec![0usize; n];
            for i in 0..n {
                let p = Vector3::new(coords[(0, i)], coords[(1, i)], coords[(2, i)]);
                let rp = matrix * p;

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
    use std::iter::repeat_n;

    use float_cmp::approx_eq;
    use rstest::rstest;
    use umol_shared::element::Element::*;

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
    #[case(
        symmetric_water(),
        3, &["A1", "A1", "B1"]
    )]
    #[case(
        symmetric_methane(),
        5, &["A1", "A1", "T2", "T2", "T2"]
    )]
    fn test_molecule_salc_basis(
        #[case] m: Molecule,
        #[case] expected_total: usize,
        #[case] expected_irreps: &[&str],
    ) {
        let sym = m.perceive_symmetry(Thresholds::default()).unwrap();
        let result = sym
            .salc_basis(&sym.permutation_basis(), Thresholds::default())
            .unwrap();

        let total: usize = result.irreps.iter().map(|ib| ib.salcs.len()).sum();
        assert_eq!(total, expected_total);

        let mut salc_symbols: Vec<String> = result
            .irreps
            .iter()
            .flat_map(|ib| repeat_n(ib.irrep.symbol().to_owned(), ib.salcs.len()))
            .collect();
        salc_symbols.sort();
        let mut expected: Vec<&str> = expected_irreps.to_vec();
        expected.sort();
        assert_eq!(salc_symbols, expected);
    }

    #[rstest]
    fn test_molecule_symmetrize() {
        let thresholds = Thresholds::default();
        let m = symmetric_water();
        let sym = m.symmetrize(thresholds).unwrap();
        assert_eq!(sym.point_group().to_string(), "C2v");
        // Symmetrized coordinates should be more symmetric than input
        let coords = sym.cartesian_coordinates();
        // y-coordinates of the two H atoms should be exactly opposite
        let y_h1 = coords[(1, 1)];
        let y_h2 = coords[(1, 2)];
        assert!(
            (y_h1 + y_h2).abs() < 1e-10,
            "H atoms not symmetric: {y_h1} vs {y_h2}"
        );
    }

    // HCl: C∞v, linear diatomic
    #[rustfmt::skip]
    fn hcl() -> Molecule {
        mol(&[Cl, H], &[
            0.000,  0.000,  0.000,
            0.000,  0.000,  1.275,
        ])
    }

    // CO₂: D∞h, linear triatomic
    #[rustfmt::skip]
    fn co2() -> Molecule {
        mol(&[O, C, O], &[
            0.000,  0.000, -1.160,
            0.000,  0.000,  0.000,
            0.000,  0.000,  1.160,
        ])
    }

    #[rstest]
    // Water (C2v): 3N=9, trans=3, rot=3, vib=3: 2A1 + B1
    #[case(
        symmetric_water(), "C2v", 9, 3, 3, 3,
        &[("A1", 2), ("B1", 1)]
    )]
    // Methane (Td): 3N=15, trans=3, rot=3, vib=9: A1+E+2T2
    #[case(
        symmetric_methane(), "Td", 15, 3, 3, 9,
        &[("A1", 1), ("E", 2), ("T2", 6)]
    )]
    // HCl (C∞v): 3N=6, trans=3, rot=2, vib=1: Σ+
    #[case(
        hcl(), "C∞v", 6, 3, 2, 1,
        &[("Σ+", 1)]
    )]
    // CO₂ (D∞h): 3N=9, trans=3, rot=2, vib=4: Σ+g + Σ+u + Πu(×2)
    #[case(
        co2(), "D∞h", 9, 3, 2, 4,
        &[("Πu", 2), ("Σ+g", 1), ("Σ+u", 1)]
    )]
    fn test_molecule_symmetry_coordinates(
        #[case] m: Molecule,
        #[case] expected_group: &str,
        #[case] expected_total: usize,
        #[case] expected_trans: usize,
        #[case] expected_rot: usize,
        #[case] expected_vib: usize,
        #[case] expected_vib_irreps: &[(&str, usize)],
    ) {
        let thresholds = Thresholds::default();
        let sym = m.perceive_symmetry(thresholds).unwrap();
        assert_eq!(sym.point_group().to_string(), expected_group);

        let sc = sym.symmetry_coordinates(thresholds);

        let total = sc.coordinates.len();
        assert_eq!(total, expected_total);

        let n_trans = sc
            .coordinates
            .iter()
            .filter(|c| c.kind == CoordinateKind::Translation)
            .count();
        let n_rot = sc
            .coordinates
            .iter()
            .filter(|c| c.kind == CoordinateKind::Rotation)
            .count();
        let n_vib = sc
            .coordinates
            .iter()
            .filter(|c| c.kind == CoordinateKind::Vibration)
            .count();
        assert_eq!(n_trans, expected_trans, "translation count");
        assert_eq!(n_rot, expected_rot, "rotation count");
        assert_eq!(n_vib, expected_vib, "vibration count");

        // Vibrational irrep decomposition
        let mut vib_irreps: Vec<(String, usize)> = Vec::new();
        for c in sc
            .coordinates
            .iter()
            .filter(|c| c.kind == CoordinateKind::Vibration)
        {
            if let Some(entry) = vib_irreps.iter_mut().find(|(s, _)| s == c.irrep.symbol()) {
                entry.1 += 1;
            } else {
                vib_irreps.push((c.irrep.symbol().to_owned(), 1));
            }
        }
        vib_irreps.sort();
        let expected_sorted: Vec<(String, usize)> = {
            let mut v: Vec<_> = expected_vib_irreps
                .iter()
                .map(|(s, n)| (s.to_string(), *n))
                .collect();
            v.sort();
            v
        };
        assert_eq!(vib_irreps, expected_sorted);

        // Orthonormality
        for (i, c1) in sc.coordinates.iter().enumerate() {
            let v1 = DVector::from_column_slice(&atom_vectors_to_flat(&c1.atom_vectors));
            for (j, c2) in sc.coordinates.iter().enumerate() {
                let v2 = DVector::from_column_slice(&atom_vectors_to_flat(&c2.atom_vectors));
                let dot = v1.dot(&v2);
                if i == j {
                    assert!(
                        (dot - 1.0).abs() < 1e-6,
                        "coordinate {i} not normalized: {dot}"
                    );
                } else {
                    assert!(
                        dot.abs() < 1e-6,
                        "coordinates {i},{j} not orthogonal: {dot}"
                    );
                }
            }
        }
    }

    #[allow(dead_code)]
    #[rustfmt::skip]
    fn methane() -> Molecule {
        mol(&[C, H, H, H, H], &[
            0.000,  0.000,  0.000,
            0.629,  0.629,  0.629,
           -0.629, -0.629,  0.629,
           -0.629,  0.629, -0.629,
            0.629, -0.629, -0.629,
        ])
    }

    #[rstest]
    #[case(SchoenfliesSymbol::Cs, "Cs")]
    #[case(SchoenfliesSymbol::Cn(2), "C2")]
    #[case(SchoenfliesSymbol::Cnv(2), "C2v")]
    #[case(SchoenfliesSymbol::Cn(1), "C1")]
    fn test_molecule_lower_symmetry_water(
        #[case] target: SchoenfliesSymbol,
        #[case] expected_name: &str,
    ) {
        let m = symmetric_water()
            .perceive_symmetry(Thresholds::default())
            .unwrap();
        assert_eq!(m.point_group().to_string(), "C2v");

        let result = m.lower_symmetry(target, Thresholds::default()).unwrap();
        assert_eq!(result.molecule.point_group().to_string(), expected_name);
        assert_eq!(result.molecule.atom_count(), 3);
        assert_eq!(result.parent_group.to_string(), "C2v");
    }

    #[rstest]
    #[case(SchoenfliesSymbol::Cnv(3), "C3v")]
    #[case(SchoenfliesSymbol::Cnv(2), "C2v")]
    #[case(SchoenfliesSymbol::Dnd(2), "D2d")]
    fn test_molecule_lower_symmetry_methane(
        #[case] target: SchoenfliesSymbol,
        #[case] expected_name: &str,
    ) {
        let m = symmetric_methane()
            .perceive_symmetry(Thresholds::default())
            .unwrap();
        assert_eq!(m.point_group().to_string(), "Td");

        let result = m.lower_symmetry(target, Thresholds::default()).unwrap();
        assert_eq!(result.molecule.point_group().to_string(), expected_name);
        assert_eq!(result.molecule.atom_count(), 5);
        assert_eq!(result.parent_group.to_string(), "Td");
    }
}
