//! Bond perception from 3D coordinates.
//!
//! Assigns bond orders to atom pairs based on interatomic distances
//! and valence constraints, using Lagrangian relaxation.

use umol_data::element::Element;
use umol_data::units::Length;
use umol_params::covalent_radii::covalent_radii;

use crate::algorithms::optimization::{
    lagrangian_relaxation, Constraint, LagrangianConfig, Variable,
};
use crate::molecule::Molecule;

/// Maximum bond order considered.
const MAX_ORDER: usize = 3;

/// Bond order likelihoods for a single atom pair.
/// `likelihoods[k]` = likelihood of bond order k, for k = 0..=MAX_ORDER.
pub type Likelihoods = [f64; MAX_ORDER + 1];

/// Distance → bond order likelihood model.
///
/// For bond orders k ≥ 1: Gaussian centered at μ_k = r_cov^(k)(a) + r_cov^(k)(b).
/// For k = 0 (no bond): sigmoid rising with distance beyond the single-bond length.
///
/// `sigma` and `sigmoid_delta` are lengths; `sigmoid_alpha` is in Bohr⁻¹.
pub struct BondDistanceModel {
    /// Gaussian width per bond order (σ_1, σ_2, σ_3).
    pub sigma: [Length; MAX_ORDER],
    /// Sigmoid steepness (Bohr⁻¹) for the no-bond likelihood.
    pub sigmoid_alpha: f64,
    /// Sigmoid midpoint shift from the single-bond distance.
    pub sigmoid_delta: Length,
}

impl Default for BondDistanceModel {
    fn default() -> Self {
        Self {
            sigma: [Length::angstrom(0.10), Length::angstrom(0.10), Length::angstrom(0.10)],
            // 15.0 Å⁻¹ converted to Bohr⁻¹
            sigmoid_alpha: 15.0 / Length::angstrom(1.0).as_bohr(),
            sigmoid_delta: Length::angstrom(0.30),
        }
    }
}

impl BondDistanceModel {
    /// Compute unnormalized likelihoods for each bond order.
    pub fn score(&self, d: Length, elem_a: Element, elem_b: Element) -> Likelihoods {
        let d = d.as_bohr();
        let ra = covalent_radii(elem_a);
        let rb = covalent_radii(elem_b);
        let mu_1 = (ra.single + rb.single).as_bohr();

        // p_0: sigmoid
        let delta = self.sigmoid_delta.as_bohr();
        let p0 = 1.0 / (1.0 + (-self.sigmoid_alpha * (d - mu_1 - delta)).exp());

        let mut out = [p0, 0.0, 0.0, 0.0];
        for (k, slot) in out.iter_mut().enumerate().skip(1) {
            if let (Some(r_a), Some(r_b)) = (ra.for_order(k as u8), rb.for_order(k as u8)) {
                let mu_k = (r_a + r_b).as_bohr();
                let sigma = self.sigma[k - 1].as_bohr();
                *slot = (-(d - mu_k).powi(2) / (2.0 * sigma * sigma)).exp();
            }
        }
        out
    }

    /// Generous distance cutoff beyond which p_0 ≈ 1.
    fn cutoff(&self, elem_a: Element, elem_b: Element) -> Length {
        let ra = covalent_radii(elem_a);
        let rb = covalent_radii(elem_b);
        (ra.single + rb.single) * 2.0
    }
}

/// Configuration for bond perception.
pub struct BondPerceptionConfig {
    /// Distance → bond order model.
    pub model: BondDistanceModel,
    /// Target valence per atom. `None` uses the octet-rule default.
    pub target_valences: Option<Vec<u8>>,
    /// Maximum Lagrangian iterations.
    pub max_iter: usize,
    /// Initial step size for subgradient updates.
    pub step_scale: f64,
}

impl Default for BondPerceptionConfig {
    fn default() -> Self {
        Self {
            model: BondDistanceModel::default(),
            target_valences: None,
            max_iter: 200,
            step_scale: 0.5,
        }
    }
}

/// Default valence for a neutral atom (octet rule).
fn default_valence(elem: Element) -> u8 {
    let ve = elem.valence_electrons();
    if ve <= 4 { ve } else { 8u8.saturating_sub(ve) }
}

/// Result of bond perception.
#[derive(Debug)]
pub struct BondPerceptionResult {
    /// Assigned bond orders: (atom_i, atom_j, order).
    pub bonds: Vec<(usize, usize, u8)>,
    /// Whether all valence constraints are satisfied.
    pub feasible: bool,
    /// Residual valence violation per atom (actual - target).
    pub valence_residuals: Vec<i32>,
}

/// Perceive bonds from a 3D molecular geometry using Lagrangian relaxation.
pub fn perceive_bonds(
    mol: &Molecule,
    config: &BondPerceptionConfig,
) -> BondPerceptionResult {
    let n = mol.atom_count();

    let target_valence: Vec<u8> = match &config.target_valences {
        Some(v) => v.clone(),
        None => (0..n).map(|i| default_valence(mol.element(i))).collect(),
    };

    // Build candidate bonds: one variable per atom pair within cutoff.
    // Each variable has domain {0, 1, ..., MAX_ORDER} and participates
    // in two constraints (one per endpoint atom).
    struct BondPair {
        i: usize,
        j: usize,
    }
    let mut pairs: Vec<BondPair> = Vec::new();
    let mut variables: Vec<Variable> = Vec::new();

    for i in 0..n {
        for j in (i + 1)..n {
            let d = mol.distance(i, j);
            if d > config.model.cutoff(mol.element(i), mol.element(j)) {
                continue;
            }
            let lik = config.model.score(d, mol.element(i), mol.element(j));
            let log_lik: Vec<f64> = lik
                .iter()
                .map(|&p| if p > 1e-30 { p.ln() } else { f64::NEG_INFINITY })
                .collect();

            variables.push(Variable {
                log_likelihoods: log_lik,
                constraints: vec![(i, 1.0), (j, 1.0)],
            });
            pairs.push(BondPair { i, j });
        }
    }

    // One constraint per atom: sum of incident bond orders = target valence
    let constraints: Vec<Constraint> = target_valence
        .iter()
        .map(|&v| Constraint { rhs: v as f64 })
        .collect();

    let result = lagrangian_relaxation(
        &variables,
        &constraints,
        &LagrangianConfig {
            max_iter: config.max_iter,
            step_scale: config.step_scale,
        },
    );

    // Collect results
    let mut bonds = Vec::new();
    for (idx, pair) in pairs.iter().enumerate() {
        let order = result.assignments[idx] as u8;
        if order > 0 {
            bonds.push((pair.i, pair.j, order));
        }
    }

    let valence_residuals: Vec<i32> = result.residuals.iter().map(|&r| r as i32).collect();
    let feasible = result.feasible;

    BondPerceptionResult {
        bonds,
        feasible,
        valence_residuals,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::{fixture, rstest};
    use umol_data::element::Element::*;
    use umol_data::spin::SpinMultiplicity;

    #[fixture]
    fn config() -> BondPerceptionConfig {
        BondPerceptionConfig::default()
    }

    fn mol(elements: &[Element], coords: &[f64]) -> Molecule {
        Molecule::from_cartesian_angstrom(elements.to_vec(), coords, 0, SpinMultiplicity::Singlet)
    }

    fn sorted_bonds(r: &BondPerceptionResult) -> Vec<(usize, usize, u8)> {
        let mut b = r.bonds.clone();
        b.sort();
        b
    }

    // Coordinates in Angstroms, from standard experimental geometries.
    //
    // Ethane:    C-C 1.54, C-H 1.09, staggered
    // Ethylene:  C=C 1.34, C-H 1.09, planar
    // Acetylene: C≡C 1.20, C-H 1.06, linear
    // Water:     O-H 0.96, H-O-H 104.5°
    // Benzene:   C-C 1.40, C-H 1.09, planar hexagon

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

    #[rustfmt::skip]
    fn acetylene() -> Molecule {
        mol(&[C, C, H, H], &[
            0.000, 0.000, 0.000,
            1.200, 0.000, 0.000,
           -1.060, 0.000, 0.000,
            2.260, 0.000, 0.000,
        ])
    }

    #[rustfmt::skip]
    fn water() -> Molecule {
        mol(&[O, H, H], &[
            0.000,  0.000,  0.000,
            0.960,  0.000,  0.000,
           -0.240,  0.930,  0.000,
        ])
    }

    #[rustfmt::skip]
    fn benzene() -> Molecule {
        mol(&[C, C, C, C, C, C, H, H, H, H, H, H], &[
             1.400,  0.000,  0.000,
             0.700,  1.212,  0.000,
            -0.700,  1.212,  0.000,
            -1.400,  0.000,  0.000,
            -0.700, -1.212,  0.000,
             0.700, -1.212,  0.000,
             2.490,  0.000,  0.000,
             1.245,  2.156,  0.000,
            -1.245,  2.156,  0.000,
            -2.490,  0.000,  0.000,
            -1.245, -2.156,  0.000,
             1.245, -2.156,  0.000,
        ])
    }

    #[rstest]
    #[case::ethane(ethane(),     7, 7,  0, 0)]  // 1 CC + 6 CH, all single
    #[case::ethylene(ethylene(), 5, 4,  1, 0)]  // 1 CC + 4 CH; 1 double, 0 triple
    #[case::acetylene(acetylene(), 3, 2, 0, 1)] // 1 CC + 2 CH; 0 double, 1 triple
    #[case::water(water(),       2, 2,  0, 0)]  // 2 OH, all single
    #[case::benzene(benzene(),  12, 9,  3, 0)]  // 6 CC + 6 CH; 3 double (Kekule), 0 triple
    fn test_perceive_bonds(
        config: BondPerceptionConfig,
        #[case] m: Molecule,
        #[case] expected_bonds: usize,
        #[case] expected_single: usize,
        #[case] expected_double: usize,
        #[case] expected_triple: usize,
    ) {
        let result = perceive_bonds(&m, &config);
        assert!(result.feasible, "valence residuals: {:?}", result.valence_residuals);

        let bonds = sorted_bonds(&result);
        assert_eq!(bonds.len(), expected_bonds, "bond count: {bonds:?}");

        let singles = bonds.iter().filter(|(_, _, o)| *o == 1).count();
        let doubles = bonds.iter().filter(|(_, _, o)| *o == 2).count();
        let triples = bonds.iter().filter(|(_, _, o)| *o == 3).count();
        assert_eq!(singles, expected_single, "single bonds: {bonds:?}");
        assert_eq!(doubles, expected_double, "double bonds: {bonds:?}");
        assert_eq!(triples, expected_triple, "triple bonds: {bonds:?}");
    }
}
