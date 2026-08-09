//! Hückel Molecular Orbital (HMO) aromaticity perception.
//!
//! Builds a π-only Hamiltonian from Van-Catledge parameters, fills it via aufbau,
//! and compares the delocalization energy per π-electron to a configured threshold.

use std::cmp::Ordering;
use std::collections::BTreeMap;

use nalgebra::{DMatrix, SymmetricEigen};
use thiserror::Error;
use umol_chem::element::Element;
use umol_graph_core::ConnectedComponentsAlgorithm;
use umol_graph_ir::ir::{
    AromaticSystemForm, AtomId, AtomView, ElementForm, Molecule, RingSet, UnpairedElectronsForm,
};
use umol_params::quantum::ppp::van_catledge::VanCatledgeParams;

use crate::ops::model::ElementScope;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum HmoError {
    /// Configuration gap: dispatcher surfaces as `Err(AromaticityError)`.
    #[error("hmo: missing parameters: {0}")]
    MissingParameters(String),
    /// Algorithm preconditions failed (empty pi-system, odd electron count,
    /// orbital count > atoms): dispatcher surfaces as `Solution::Contradictory`.
    #[error("hmo: invalid input: {0}")]
    InvalidInput(String),
    /// Atom data not ground enough to evaluate: dispatcher surfaces as
    /// `Solution::Underdetermined`.
    #[error("hmo: undetermined atom data: {0}")]
    UndeterminedAtom(String),
}

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

    fn is_element_supported(&self, element: Element) -> bool {
        let has_params = (0..=2).any(|n| VanCatledgeParams::h_x(element, n).is_some());
        match &self.element_scope {
            ElementScope::Any => has_params,
            ElementScope::AllowList(list) => list.contains(&element) && has_params,
        }
    }

    pub fn find_from_rings<F>(
        &self,
        ast: &Molecule,
        rings: &RingSet,
        connected_components_algorithm: ConnectedComponentsAlgorithm,
        electrons_at: &F,
    ) -> Result<Vec<(Vec<AtomId>, AromaticSystemForm)>, HmoError>
    where
        F: Fn(&AtomView<'_>) -> Option<u8>,
    {
        let pi_atoms: Vec<AtomId> = ast
            .atoms()
            .iter()
            .filter_map(|view| {
                let element = match view.ast.element {
                    ElementForm::Lit(e) => e,
                    _ => return None,
                };
                if !self.is_element_supported(element) {
                    return None;
                }
                electrons_at(&view).map(|_| view.id)
            })
            .collect();

        if pi_atoms.is_empty() {
            return Ok(Vec::new());
        }

        let subgraph = ast.induced_subgraph(&pi_atoms);
        let extracted = ast.extract(&subgraph);
        let mut sorted_host: Vec<AtomId> = subgraph
            .atoms()
            .matched_pairs()
            .iter()
            .map(|&(_, host)| host)
            .collect();
        sorted_host.sort_unstable();
        let components: Vec<Vec<AtomId>> = extracted
            .graph()
            .enumerate_connected_components(connected_components_algorithm)
            .into_iter()
            .map(|c| c.into_iter().map(|sub| sorted_host[sub.index()]).collect())
            .collect();

        let mut candidates = Vec::new();
        for component in &components {
            let has_ring_atom = component.iter().any(|a| rings.contains_atom(*a));
            if !has_ring_atom || component.len() < 3 {
                continue;
            }

            let result = self.build_calculator(ast, component, electrons_at)?.solve();

            let de_per_electron = if result.electron_count > 0 {
                result.delocalization_energy / result.electron_count as f64
            } else {
                0.0
            };

            if de_per_electron >= self.stabilization_threshold {
                let mut atoms = result.atom_indices.clone();
                atoms.sort_unstable();

                let electrons: Vec<i64> = atoms
                    .iter()
                    .map(|&atom| electrons_at(&ast.atom(atom)).unwrap_or(0) as i64)
                    .collect();

                candidates.push((
                    atoms,
                    AromaticSystemForm::from_electrons(electrons)
                        .with_charge(0)
                        .with_unpaired_electrons(UnpairedElectronsForm::closed_shell()),
                ));
            }
        }

        Ok(candidates)
    }

    pub(crate) fn build_calculator<F>(
        &self,
        ast: &Molecule,
        pi_atoms: &[AtomId],
        electrons_at: &F,
    ) -> Result<HmoCalculator, HmoError>
    where
        F: Fn(&AtomView<'_>) -> Option<u8>,
    {
        let subgraph = ast.induced_subgraph(pi_atoms);
        let n = pi_atoms.len();
        let mut hamiltonian = DMatrix::zeros(n, n);
        let mut electron_count: u32 = 0;
        let mut atom_types: Vec<(Element, u8)> = Vec::with_capacity(n);
        for (i, &atom) in pi_atoms.iter().enumerate() {
            let view = ast.atom(atom);
            let element = match view.ast.element {
                ElementForm::Lit(e) => e,
                _ => {
                    return Err(HmoError::UndeterminedAtom(
                        "undetermined element".to_string(),
                    ))
                }
            };
            let valence = electrons_at(&view).ok_or_else(|| {
                HmoError::UndeterminedAtom("undetermined aromatic valence".to_string())
            })?;
            let hx = VanCatledgeParams::h_x(element, valence).ok_or_else(|| {
                HmoError::MissingParameters(format!(
                    "no Van-Catledge parameters for {:?} with {} pi-electrons",
                    element, valence
                ))
            })?;
            hamiltonian[(i, i)] = hx;
            atom_types.push((element, valence));
            electron_count += valence as u32;
        }

        let mut bond_positions = Vec::with_capacity(subgraph.bonds().matched_pair_count());
        for &(_, bid) in subgraph.bonds().matched_pairs() {
            let [ha, hb] = ast.bond(bid).atom_ids();
            let i = subgraph.atoms().left_of(ha).unwrap().index();
            let j = subgraph.atoms().left_of(hb).unwrap().index();
            let k = VanCatledgeParams::k_xy(atom_types[i], atom_types[j]).ok_or_else(|| {
                HmoError::MissingParameters(format!(
                    "no Van-Catledge k_XY for {:?}-{:?}",
                    atom_types[i], atom_types[j]
                ))
            })?;
            hamiltonian[(i, j)] = k;
            hamiltonian[(j, i)] = k;
            bond_positions.push((i, j));
        }

        HmoCalculator::new(
            pi_atoms.to_vec(),
            electron_count,
            hamiltonian,
            bond_positions,
        )
    }
}

pub(crate) struct HmoCalculator {
    pi_atoms: Vec<AtomId>,
    electron_count: u32,
    hamiltonian: DMatrix<f64>,
    bond_positions: Vec<(usize, usize)>,
}

impl HmoCalculator {
    fn new(
        pi_atoms: Vec<AtomId>,
        electron_count: u32,
        hamiltonian: DMatrix<f64>,
        bond_positions: Vec<(usize, usize)>,
    ) -> Result<Self, HmoError> {
        if pi_atoms.is_empty() {
            return Err(HmoError::InvalidInput(
                "empty pi-system for HMO".to_string(),
            ));
        }
        if electron_count == 0 {
            return Err(HmoError::InvalidInput("zero pi-electrons".to_string()));
        }
        if !electron_count.is_multiple_of(2) {
            return Err(HmoError::InvalidInput(
                "open-shell pi-system (odd electron count) not supported by HMO".to_string(),
            ));
        }
        let orbital_count = (electron_count / 2) as usize;
        if orbital_count > pi_atoms.len() {
            return Err(HmoError::InvalidInput(
                "more electron pairs than orbitals".to_string(),
            ));
        }
        Ok(Self {
            pi_atoms,
            electron_count,
            hamiltonian,
            bond_positions,
        })
    }

    pub(crate) fn solve(&self) -> HmoOutput {
        let n = self.pi_atoms.len();
        let eigen = SymmetricEigen::new(self.hamiltonian.clone());
        let eigenvalues = eigen.eigenvalues;
        let eigenvectors = eigen.eigenvectors;

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

        let mut density = DMatrix::zeros(n, n);
        for &k in &sorted_indices[..orbital_count] {
            let col = eigenvectors.column(k);
            for i in 0..n {
                for j in 0..n {
                    density[(i, j)] += 2.0 * col[i] * col[j];
                }
            }
        }

        let mut bond_orders: BTreeMap<(AtomId, AtomId), f64> = BTreeMap::new();
        for &(i, j) in &self.bond_positions {
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

#[derive(Debug)]
pub struct HmoOutput {
    pub atom_indices: Vec<AtomId>,
    pub delocalization_energy: f64,
    pub electron_count: u32,
    pub bond_orders: BTreeMap<(AtomId, AtomId), f64>,
}

#[cfg(test)]
mod tests {
    use float_cmp::*;
    use rstest::*;
    use umol_chem::element::Element;
    use umol_graph_ir::ir::{
        AromaticValenceForm, AtomConstraintForm, AtomForm, AtomId, BondForm, Molecule,
        MoleculeEntries, NumForm, RingConfig, RingModel, RingSetKind,
    };

    use super::*;

    fn aromatic(element: Element, pi: i64) -> (AtomForm, Option<i64>) {
        (AtomForm::from_element(element), Some(pi))
    }

    fn apply_pi(specs: Vec<(AtomForm, Option<i64>)>) -> Vec<AtomForm> {
        specs
            .into_iter()
            .map(|(mut atom, pi)| {
                if let Some(n) = pi {
                    atom.constraints.set(AtomConstraintForm::AromaticValence(
                        AromaticValenceForm::Aromatic(NumForm::Lit(n)),
                    ));
                }
                atom
            })
            .collect()
    }

    fn make_ring(specs: Vec<(AtomForm, Option<i64>)>) -> Molecule {
        let n = specs.len();
        let atoms = apply_pi(specs);
        let bonds: Vec<_> = (0..n)
            .map(|i| {
                (
                    AtomId(i as u32),
                    AtomId(((i + 1) % n) as u32),
                    BondForm::from_order(1),
                )
            })
            .collect();
        Molecule::from_entries(MoleculeEntries {
            atoms,
            bonds,
            ..Default::default()
        })
    }

    fn make_fused(specs: Vec<(AtomForm, Option<i64>)>, edges: &[(usize, usize)]) -> Molecule {
        let atoms = apply_pi(specs);
        let bonds: Vec<_> = edges
            .iter()
            .map(|&(a, b)| (AtomId(a as u32), AtomId(b as u32), BondForm::from_order(1)))
            .collect();
        Molecule::from_entries(MoleculeEntries {
            atoms,
            bonds,
            ..Default::default()
        })
    }

    fn solve_hmo(model: &HmoAromaticity, ast: &Molecule) -> HmoOutput {
        let atoms: Vec<AtomId> = ast.atoms().ids().collect();
        model
            .build_calculator(ast, &atoms, &|v| match v
                .ast
                .constraints
                .aromatic_valence()
                .unwrap_or(&AromaticValenceForm::Undetermined)
            {
                AromaticValenceForm::Aromatic(NumForm::Lit(n)) if *n >= 0 => Some(*n as u8),
                _ => None,
            })
            .unwrap()
            .solve()
    }

    #[fixture]
    fn hmo_model() -> HmoAromaticity {
        HmoAromaticity::new(
            ElementScope::AllowList(vec![Element::C, Element::N, Element::O, Element::S]),
            0.1,
        )
    }

    #[fixture]
    fn benzene() -> Molecule {
        make_ring(vec![aromatic(Element::C, 1); 6])
    }

    #[rustfmt::skip]
    #[fixture]
    fn naphthalene() -> Molecule {
        make_fused(
            vec![aromatic(Element::C, 1); 10],
            &[
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
                (3, 6), (6, 7), (7, 8), (8, 9), (9, 4),
            ],
        )
    }

    #[rustfmt::skip]
    #[fixture]
    fn azulene() -> Molecule {
        make_fused(
            vec![aromatic(Element::C, 1); 10],
            &[
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
                (0, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 4),
            ],
        )
    }

    #[fixture]
    fn pyridine() -> Molecule {
        make_ring(vec![
            aromatic(Element::N, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
        ])
    }

    #[fixture]
    fn pyrrole() -> Molecule {
        make_ring(vec![
            aromatic(Element::N, 2),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
        ])
    }

    #[fixture]
    fn cyclobutadiene() -> Molecule {
        make_ring(vec![aromatic(Element::C, 1); 4])
    }

    #[rstest]
    #[case::benzene(benzene(), 2.0)]
    #[case::naphthalene(naphthalene(), 3.683)]
    #[case::azulene(azulene(), 3.364)]
    #[case::pyridine(pyridine(), 2.614)]
    #[case::pyrrole(pyrrole(), 2.200)]
    fn test_hmo_aromaticity_delocalization_energy(
        hmo_model: HmoAromaticity,
        #[case] ast: Molecule,
        #[case] expected_de: f64,
    ) {
        let result = solve_hmo(&hmo_model, &ast);
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
    fn test_hmo_aromaticity_bond_orders(
        hmo_model: HmoAromaticity,
        #[case] ast: Molecule,
        #[case] expected_count: usize,
        #[case] expected_min: f64,
        #[case] expected_max: f64,
    ) {
        let result = solve_hmo(&hmo_model, &ast);
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
    fn test_hmo_aromaticity_find_from_rings(
        hmo_model: HmoAromaticity,
        #[case] ast: Molecule,
        #[case] expected_systems: usize,
        #[case] expected_atoms: Option<usize>,
    ) {
        let ring_info = ast
            .rings(
                RingModel {
                    kind: RingSetKind::Relevant,
                    max_ring_size: 22,
                },
                RingConfig::default(),
            )
            .into_ring_set();
        let systems = hmo_model
            .find_from_rings(
                &ast,
                &ring_info,
                ConnectedComponentsAlgorithm::Bfs,
                &|v| match v
                    .ast
                    .constraints
                    .aromatic_valence()
                    .unwrap_or(&AromaticValenceForm::Undetermined)
                {
                    AromaticValenceForm::Aromatic(NumForm::Lit(n)) if *n >= 0 => Some(*n as u8),
                    _ => None,
                },
            )
            .unwrap();
        assert_eq!(systems.len(), expected_systems);
        assert_eq!(systems.first().map(|s| s.0.len()), expected_atoms);
    }

    #[rstest]
    fn test_hmo_aromaticity_output(hmo_model: HmoAromaticity, benzene: Molecule) {
        let output = solve_hmo(&hmo_model, &benzene);

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
    fn test_hmo_aromaticity_hamiltonian(hmo_model: HmoAromaticity, pyridine: Molecule) {
        let atoms: Vec<AtomId> = pyridine.atoms().ids().collect();
        let calc = hmo_model
            .build_calculator(&pyridine, &atoms, &|v| match v
                .ast
                .constraints
                .aromatic_valence()
                .unwrap_or(&AromaticValenceForm::Undetermined)
            {
                AromaticValenceForm::Aromatic(NumForm::Lit(n)) if *n >= 0 => Some(*n as u8),
                _ => None,
            })
            .unwrap();
        let h = &calc.hamiltonian;
        assert_eq!(h.nrows(), 6);
        assert_eq!(h.ncols(), 6);
        assert!(approx_eq!(f64, h[(0, 0)], 0.51, epsilon = 0.01));
    }
}
