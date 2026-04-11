//! Hueckel Molecular Orbital (HMO) aromaticity model.
//!
//! Uses Van-Catledge parameters and aufbau filling to compute the delocalization energy and pi-bond orders.
//! Compares the delocalization energy per pi-electron to a threshold to determine if the system is aromatic.

use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};

use nalgebra::{DMatrix, SymmetricEigen};
use umol_shared::Element;
use umol_params::quantum::ppp::van_catledge::VanCatledgeParams;

use super::{AromaticContribution, AromaticSystem, AromaticityError};
use crate::graph_ir::config::ElementScope;
use crate::graph_ir::molecule::AtomIndex;
use crate::graph_ir::molecule_builder::MoleculeBuilder;
use crate::graph_ir::rings::{Ring, RingSet};

#[derive(Clone, Debug)]
pub struct HmoAromaticity {
    pub element_scope: ElementScope,
    pub stabilization_threshold: f64,
}

impl HmoAromaticity {
    pub fn new(element_scope: ElementScope, stabilization_threshold: f64) -> Self {
        Self {
            element_scope,
            stabilization_threshold,
        }
    }

    pub(crate) fn build_calculator(
        &self,
        builder: &MoleculeBuilder,
        pi_atoms: &[AtomIndex],
    ) -> Result<HmoCalculator, AromaticityError> {
        let atom_to_idx: HashMap<AtomIndex, usize> =
            pi_atoms.iter().enumerate().map(|(i, &a)| (a, i)).collect();

        let mut h_values = Vec::with_capacity(pi_atoms.len());
        let mut electron_count: u32 = 0;
        let mut atom_types: Vec<(Element, u8)> = Vec::with_capacity(pi_atoms.len());
        for &atom in pi_atoms {
            let atom_data = builder
                .atom(atom)
                .ok_or_else(|| AromaticityError::HmoMissingAtom("missing atom".to_string()))?;
            let element = atom_data.element();
            let valence = builder.atom_aromatic_valence(atom);
            let hx = VanCatledgeParams::h_x(element, valence).ok_or_else(|| {
                AromaticityError::HmoMissingParameters(format!(
                    "no Van-Catledge parameters for {:?} with {} pi-electrons",
                    element, valence
                ))
            })?;
            h_values.push(hx);
            atom_types.push((element, valence));
            electron_count += valence as u32;
        }

        let mut bonds = Vec::new();
        for (i, &atom_i) in pi_atoms.iter().enumerate() {
            for neighbor in builder.atom_neighbor_indices(atom_i) {
                if let Some(&j) = atom_to_idx.get(&neighbor) {
                    if j > i {
                        let k = VanCatledgeParams::k_xy(atom_types[i], atom_types[j]).ok_or_else(
                            || {
                                AromaticityError::HmoMissingParameters(format!(
                                    "no Van-Catledge k_XY for {:?}-{:?}",
                                    atom_types[i], atom_types[j]
                                ))
                            },
                        )?;
                        bonds.push((i, j, k));
                    }
                }
            }
        }

        HmoCalculator::new(pi_atoms.to_vec(), electron_count, h_values, bonds)
    }

    fn is_element_supported(&self, element: Element) -> bool {
        let has_params = (0..=2).any(|n| VanCatledgeParams::h_x(element, n).is_some());
        match &self.element_scope {
            ElementScope::Any => has_params,
            ElementScope::AllowList(list) => list.contains(&element) && has_params,
        }
    }

    pub fn find_from_rings(
        &self,
        builder: &MoleculeBuilder,
        rings: &RingSet,
    ) -> Result<Vec<AromaticSystem>, AromaticityError> {
        let pi_atoms: Vec<AtomIndex> = builder
            .atom_indices()
            .filter(|&atom| {
                let atom_data = match builder.atom(atom) {
                    Some(a) => a,
                    None => return false,
                };
                if !self.is_element_supported(atom_data.element()) {
                    return false;
                }
                builder.atom_has_aromatic_candidate(atom)
            })
            .collect();

        if pi_atoms.is_empty() {
            return Ok(Vec::new());
        }

        let pi_set: HashSet<AtomIndex> = pi_atoms.iter().copied().collect();

        let mut visited: HashSet<AtomIndex> = HashSet::new();
        let mut components: Vec<Vec<AtomIndex>> = Vec::new();
        for &atom in &pi_atoms {
            if visited.contains(&atom) {
                continue;
            }
            let mut component = Vec::new();
            let mut stack = vec![atom];
            visited.insert(atom);
            while let Some(current) = stack.pop() {
                component.push(current);
                for neighbor in builder.atom_neighbor_indices(current) {
                    if pi_set.contains(&neighbor) && visited.insert(neighbor) {
                        stack.push(neighbor);
                    }
                }
            }
            component.sort_unstable();
            components.push(component);
        }

        let mut candidates = Vec::new();
        for component in &components {
            let has_ring_atom = component.iter().any(|a| rings.is_ring_atom(*a));
            if !has_ring_atom || component.len() < 3 {
                continue;
            }

            let result = self.build_calculator(builder, component)?.solve();

            let de_per_electron = if result.electron_count > 0 {
                result.delocalization_energy / result.electron_count as f64
            } else {
                0.0
            };

            if de_per_electron >= self.stabilization_threshold {
                let contributions: Vec<AromaticContribution> = result
                    .atom_indices
                    .iter()
                    .map(|&atom| {
                        let valence = builder.atom_aromatic_valence(atom);
                        AromaticContribution::new(atom, valence)
                    })
                    .collect();

                let component_set: HashSet<AtomIndex> = component.iter().copied().collect();
                let rings: Vec<Ring> = rings
                    .ring_indices()
                    .filter_map(|i| rings.ring(i))
                    .filter(|cycle| cycle.atoms().iter().all(|a| component_set.contains(a)))
                    .cloned()
                    .collect();

                candidates.push(AromaticSystem::with_rings(contributions, rings));
            }
        }

        Ok(candidates)
    }
}

/// HMO calculator: parameter-agnostic Hueckel MO solver.
pub(crate) struct HmoCalculator {
    pi_atoms: Vec<AtomIndex>,
    electron_count: u32,
    h_values: Vec<f64>,
    bonds: Vec<(usize, usize, f64)>,
}

impl HmoCalculator {
    fn new(
        pi_atoms: Vec<AtomIndex>,
        electron_count: u32,
        h_values: Vec<f64>,
        bonds: Vec<(usize, usize, f64)>,
    ) -> Result<Self, AromaticityError> {
        if pi_atoms.is_empty() {
            return Err(AromaticityError::HmoInvalidInput(
                "empty pi-system for HMO".to_string(),
            ));
        }
        if electron_count == 0 {
            return Err(AromaticityError::HmoInvalidInput(
                "zero pi-electrons".to_string(),
            ));
        }
        if !electron_count.is_multiple_of(2) {
            return Err(AromaticityError::HmoInvalidInput(
                "open-shell pi-system (odd electron count) not supported by HMO".to_string(),
            ));
        }
        let orbital_count = (electron_count / 2) as usize;
        if orbital_count > pi_atoms.len() {
            return Err(AromaticityError::HmoInvalidInput(
                "more electron pairs than orbitals".to_string(),
            ));
        }
        Ok(Self {
            pi_atoms,
            electron_count,
            h_values,
            bonds,
        })
    }

    pub(crate) fn hamiltonian(&self) -> DMatrix<f64> {
        let n = self.pi_atoms.len();
        let mut h = DMatrix::zeros(n, n);
        for (i, &hx) in self.h_values.iter().enumerate() {
            h[(i, i)] = hx;
        }
        for &(i, j, k) in &self.bonds {
            h[(i, j)] = k;
            h[(j, i)] = k;
        }
        h
    }

    pub(crate) fn solve(&self) -> HmoOutput {
        let n = self.pi_atoms.len();
        let h = self.hamiltonian();
        let eigen = SymmetricEigen::new(h);
        let eigenvalues = eigen.eigenvalues;
        let eigenvectors = eigen.eigenvectors;

        // Sort by descending eigenvalue (most bonding first in HMO convention
        // where H is in units of beta with alpha=0).
        let mut sorted_indices: Vec<usize> = (0..n).collect();
        sorted_indices.sort_by(|&a, &b| {
            eigenvalues[b]
                .partial_cmp(&eigenvalues[a])
                .unwrap_or(Ordering::Equal)
        });

        let orbital_count = (self.electron_count / 2) as usize;

        let total_pi_energy: f64 = sorted_indices[..orbital_count]
            .iter()
            .map(|&i| eigenvalues[i])
            .sum::<f64>()
            * 2.0;
        let reference_energy = 2.0 * (self.electron_count / 2) as f64;
        let delocalization_energy = total_pi_energy - reference_energy;

        // Density matrix P_ij = 2 * sum_{occupied k} c_{ik} * c_{jk}.
        let mut density = DMatrix::zeros(n, n);
        for &k in &sorted_indices[..orbital_count] {
            let col = eigenvectors.column(k);
            for i in 0..n {
                for j in 0..n {
                    density[(i, j)] += 2.0 * col[i] * col[j];
                }
            }
        }

        let mut bond_orders: BTreeMap<(AtomIndex, AtomIndex), f64> = BTreeMap::new();
        for &(i, j, _) in &self.bonds {
            let (a, b) = if self.pi_atoms[i] < self.pi_atoms[j] {
                (self.pi_atoms[i], self.pi_atoms[j])
            } else {
                (self.pi_atoms[j], self.pi_atoms[i])
            };
            bond_orders.insert((a, b), density[(i, j)]);
        }

        HmoOutput {
            atom_indices: self.pi_atoms.clone(),
            delocalization_energy,
            electron_count: self.electron_count,
            bond_orders,
        }
    }
}

/// Output of an HMO calculation on a single pi-system.
#[derive(Debug)]
pub struct HmoOutput {
    pub atom_indices: Vec<AtomIndex>,
    pub delocalization_energy: f64,
    pub electron_count: u32,
    pub bond_orders: BTreeMap<(AtomIndex, AtomIndex), f64>,
}

#[cfg(test)]
mod tests {
    use float_cmp::*;
    use rstest::*;

    use super::*;
    use crate::graph_ir::bond_pattern::BondPattern;
    use crate::graph_ir::config::RingEnumerationStrategy;
    use crate::graph_ir::rings::{RingEnumerator, RingFamily};

    const C1: &str = "C#h#v2#a";

    fn make_ring(atom_specs: &[&str]) -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = atom_specs
            .iter()
            .map(|s| builder.add_resolved_atom(s.parse().unwrap()))
            .collect();
        let n = atoms.len();
        for i in 0..n {
            builder.add_bond_unchecked(atoms[i], atoms[(i + 1) % n], BondPattern::new(1));
        }
        builder
    }

    fn make_fused(atom_specs: &[&str], edges: &[(usize, usize)]) -> MoleculeBuilder {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = atom_specs
            .iter()
            .map(|s| builder.add_resolved_atom(s.parse().unwrap()))
            .collect();
        for &(a, b) in edges {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondPattern::new(1));
        }
        builder
    }

    #[fixture]
    fn hmo_model() -> HmoAromaticity {
        HmoAromaticity::new(
            ElementScope::AllowList(vec![Element::C, Element::N, Element::O, Element::S]),
            0.1,
        )
    }

    #[fixture]
    fn benzene() -> MoleculeBuilder {
        make_ring(&[C1; 6])
    }

    #[rustfmt::skip]
    #[fixture]
    fn naphthalene() -> MoleculeBuilder {
        make_fused(
            &[C1; 10],
            &[
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
                (3, 6), (6, 7), (7, 8), (8, 9), (9, 4),
            ],
        )
    }

    #[rustfmt::skip]
    #[fixture]
    fn azulene() -> MoleculeBuilder {
        make_fused(
            &[C1; 10],
            &[
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
                (0, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 4),
            ],
        )
    }

    #[fixture]
    fn pyridine() -> MoleculeBuilder {
        make_ring(&["N#n#v2#a", C1, C1, C1, C1, C1])
    }

    #[fixture]
    fn pyrrole() -> MoleculeBuilder {
        make_ring(&["N#h#v2#a2", C1, C1, C1, C1])
    }

    #[fixture]
    fn cyclobutadiene() -> MoleculeBuilder {
        make_ring(&[C1; 4])
    }

    fn solve(model: &HmoAromaticity, builder: &MoleculeBuilder) -> HmoOutput {
        let atoms: Vec<AtomIndex> = builder.atom_indices().collect();
        model.build_calculator(builder, &atoms).unwrap().solve()
    }

    #[rstest]
    #[case::benzene(benzene(), 2.0)]
    #[case::naphthalene(naphthalene(), 3.683)]
    #[case::azulene(azulene(), 3.364)]
    #[case::pyridine(pyridine(), 2.614)]
    #[case::pyrrole(pyrrole(), 2.200)]
    fn test_hmo_model_delocalization_energy(
        hmo_model: HmoAromaticity,
        #[case] builder: MoleculeBuilder,
        #[case] expected_de: f64,
    ) {
        let result = solve(&hmo_model, &builder);
        assert!(approx_eq!(
            f64,
            result.delocalization_energy,
            expected_de,
            epsilon = 0.005
        ));
    }

    #[rstest]
    #[case::benzene(benzene(), 6, 0.667, 0.667)]
    #[case::azulene(azulene(), 11, 0.401, 0.664)]
    fn test_hmo_model_bond_orders(
        hmo_model: HmoAromaticity,
        #[case] builder: MoleculeBuilder,
        #[case] expected_count: usize,
        #[case] expected_min: f64,
        #[case] expected_max: f64,
    ) {
        let result = solve(&hmo_model, &builder);
        assert_eq!(result.bond_orders.len(), expected_count);
        let min = result
            .bond_orders
            .values()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let max = result
            .bond_orders
            .values()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(approx_eq!(f64, min, expected_min, epsilon = 0.001));
        assert!(approx_eq!(f64, max, expected_max, epsilon = 0.001));
    }

    #[rstest]
    #[case::benzene(benzene(), 1, Some(6))]
    #[case::naphthalene(naphthalene(), 1, Some(10))]
    #[case::cyclobutadiene(cyclobutadiene(), 0, None)]
    fn test_hmo_model_find_from_rings(
        hmo_model: HmoAromaticity,
        #[case] builder: MoleculeBuilder,
        #[case] expected_systems: usize,
        #[case] expected_atoms: Option<usize>,
    ) {
        let ring_info =
            RingEnumerator::new(RingFamily::Simple, &RingEnumerationStrategy::default())
                .enumerate_builder(&builder);
        let systems = hmo_model.find_from_rings(&builder, &ring_info).unwrap();
        assert_eq!(systems.len(), expected_systems);
        assert_eq!(
            systems.first().map(|s| s.contributions().len()),
            expected_atoms
        );
    }

    #[rstest]
    fn test_hmo_output(hmo_model: HmoAromaticity, benzene: MoleculeBuilder) {
        let atoms: Vec<AtomIndex> = benzene.atom_indices().collect();
        let output = hmo_model
            .build_calculator(&benzene, &atoms)
            .unwrap()
            .solve();

        assert_eq!(output.bond_orders.len(), 6);
        let min_bo = output
            .bond_orders
            .values()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let max_bo = output
            .bond_orders
            .values()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(approx_eq!(f64, min_bo, 0.667, epsilon = 0.001));
        assert!(approx_eq!(f64, max_bo, 0.667, epsilon = 0.001));
        assert!(approx_eq!(
            f64,
            output.delocalization_energy,
            2.0,
            epsilon = 0.005
        ));
        assert_eq!(output.electron_count, 6);
    }

    #[rstest]
    fn test_hmo_model_hamiltonian(hmo_model: HmoAromaticity, pyridine: MoleculeBuilder) {
        let atoms: Vec<AtomIndex> = pyridine.atom_indices().collect();
        let calc = hmo_model.build_calculator(&pyridine, &atoms).unwrap();
        let h = calc.hamiltonian();
        assert_eq!(h.nrows(), 6);
        assert_eq!(h.ncols(), 6);
        assert!(approx_eq!(f64, h[(0, 0)], 0.51, epsilon = 0.01));
    }
}
