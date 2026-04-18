//! Hueckel Molecular Orbital (HMO) aromaticity model.
//!
//! Uses Van-Catledge parameters and aufbau filling to compute the delocalization energy and pi-bond orders.
//! Compares the delocalization energy per pi-electron to a threshold to determine if the system is aromatic.

use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};

use nalgebra::{DMatrix, SymmetricEigen};
use umol_shared::element::Element;
use umol_shared::value_ast::ValueAst;
use umol_params::quantum::ppp::van_catledge::VanCatledgeParams;

use umol_shared::atom_ast::ElementAst;

use super::AromaticityError;
use super::ElementScope;
use crate::ast::AtomIdx;
use crate::ast::aromatic::AromaticSystem;
use crate::ast::constraint::{AromaticValenceConstraint, AtomConstraint, BondConstraint};
use crate::ast::aromatic::AromaticSystemAst;
use crate::ast::molecule::MoleculeAst;
use crate::ast::rings::RingSet;

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

    pub fn find_from_rings(
        &self,
        ast: &MoleculeAst,
        rings: &RingSet,
    ) -> Result<Vec<AromaticSystem>, AromaticityError> {
        let pi_atoms: Vec<AtomIdx> = ast
            .atoms()
            .iter()
            .filter_map(|view| {
                let element = match view.data.element {
                    ElementAst::Lit(e) => e,
                    _ => return None,
                };
                if !self.is_element_supported(element) {
                    return None;
                }
                ast.atom_aromatic_pi_electrons(view.idx).map(|_| view.idx)
            })
            .collect();

        if pi_atoms.is_empty() {
            return Ok(Vec::new());
        }

        let pi_set: HashSet<AtomIdx> = pi_atoms.iter().copied().collect();

        let mut visited: HashSet<AtomIdx> = HashSet::new();
        let mut components: Vec<Vec<AtomIdx>> = Vec::new();
        for &atom in &pi_atoms {
            if visited.contains(&atom) {
                continue;
            }
            let mut component = Vec::new();
            let mut stack = vec![atom];
            visited.insert(atom);
            while let Some(current) = stack.pop() {
                component.push(current);
                for neighbor in ast.neighbors(current) {
                    let n = neighbor.atom;
                    if pi_set.contains(&n) && visited.insert(n) {
                        stack.push(n);
                    }
                }
            }
            component.sort_unstable();
            components.push(component);
        }

        let mut candidates = Vec::new();
        for component in &components {
            let has_ring_atom = component.iter().any(|a| rings.contains_atom(*a));
            if !has_ring_atom || component.len() < 3 {
                continue;
            }

            let result = self.build_calculator(ast, component)?.solve();

            let de_per_electron = if result.electron_count > 0 {
                result.delocalization_energy / result.electron_count as f64
            } else {
                0.0
            };

            if de_per_electron >= self.stabilization_threshold {
                let mut atoms = result.atom_indices.clone();
                atoms.sort_unstable();

                let atom_constraints: Vec<AtomConstraint> = atoms
                    .iter()
                    .map(|&atom| {
                        let pi = ast.atom_aromatic_pi_electrons(atom).unwrap_or(0);
                        AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(
                            ValueAst::Lit(pi as i64),
                        ))
                    })
                    .collect();

                let bonds = ast.induced_bonds(&atoms);
                let bond_constraints = vec![BondConstraint::Aromatic; bonds.len()];

                candidates.push(AromaticSystem::new(
                    atoms,
                    bonds,
                    AromaticSystemAst::default(),
                    atom_constraints,
                    bond_constraints,
                ));
            }
        }

        Ok(candidates)
    }

    pub(crate) fn build_calculator(
        &self,
        ast: &MoleculeAst,
        pi_atoms: &[AtomIdx],
    ) -> Result<HmoCalculator, AromaticityError> {
        let atom_to_idx: HashMap<AtomIdx, usize> =
            pi_atoms.iter().enumerate().map(|(i, &a)| (a, i)).collect();

        let mut h_values = Vec::with_capacity(pi_atoms.len());
        let mut electron_count: u32 = 0;
        let mut atom_types: Vec<(Element, u8)> = Vec::with_capacity(pi_atoms.len());
        for &atom in pi_atoms {
            let atom_ast = ast.atom(atom);
            let element = match atom_ast.data.element {
                ElementAst::Lit(e) => e,
                _ => {
                    return Err(AromaticityError::HmoMissingAtom(
                        "undetermined element".to_string(),
                    ))
                }
            };
            let valence = ast.atom_aromatic_pi_electrons(atom).ok_or_else(|| {
                AromaticityError::HmoMissingAtom("undetermined aromatic valence".to_string())
            })?;
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
            for neighbor in ast.neighbors(atom_i) {
                let n = neighbor.atom;
                if let Some(&j) = atom_to_idx.get(&n) {
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
}

pub(crate) struct HmoCalculator {
    pi_atoms: Vec<AtomIdx>,
    electron_count: u32,
    h_values: Vec<f64>,
    bonds: Vec<(usize, usize, f64)>,
}

impl HmoCalculator {
    fn new(
        pi_atoms: Vec<AtomIdx>,
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

        let mut bond_orders: BTreeMap<(AtomIdx, AtomIdx), f64> = BTreeMap::new();
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

#[derive(Debug)]
pub struct HmoOutput {
    pub atom_indices: Vec<AtomIdx>,
    pub delocalization_energy: f64,
    pub electron_count: u32,
    pub bond_orders: BTreeMap<(AtomIdx, AtomIdx), f64>,
}

#[cfg(test)]
mod tests {
    use float_cmp::*;
    use rstest::*;
    use umol_shared::element::Element;
    use umol_shared::value_ast::ValueAst;

    use super::*;
    use crate::ast::AtomIdx;
    use crate::ast::atom::AtomAst;
    use crate::ast::bond::BondAst;
    use crate::ast::constraint::{AromaticValenceConstraint, AtomConstraint, MoleculeConstraint};
    use crate::ast::molecule::MoleculeAst;
    use crate::ast::rings::RingEnumerationStrategy;
    use crate::ast::rings::{RingEnumerator, RingFamily};

    fn aromatic(element: Element, pi: i64) -> (AtomAst, Option<i64>) {
        (AtomAst::from_element(element), Some(pi))
    }

    fn pi_constraints(specs: &[(AtomAst, Option<i64>)]) -> Vec<MoleculeConstraint> {
        specs
            .iter()
            .enumerate()
            .filter_map(|(i, (_, pi))| {
                pi.map(|n| {
                    MoleculeConstraint::AtomPred(
                        AtomIdx(i as u32),
                        AtomConstraint::AromaticValence(AromaticValenceConstraint::Value(
                            ValueAst::Lit(n),
                        )),
                    )
                })
            })
            .collect()
    }

    fn make_ring(specs: Vec<(AtomAst, Option<i64>)>) -> MoleculeAst {
        let n = specs.len();
        let constraints = pi_constraints(&specs);
        let atoms: Vec<AtomAst> = specs.into_iter().map(|(a, _)| a).collect();
        let bonds: Vec<_> = (0..n)
            .map(|i| {
                (
                    AtomIdx(i as u32),
                    AtomIdx(((i + 1) % n) as u32),
                    BondAst::from_order(1),
                )
            })
            .collect();
        MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], constraints)
    }

    fn make_fused(specs: Vec<(AtomAst, Option<i64>)>, edges: &[(usize, usize)]) -> MoleculeAst {
        let constraints = pi_constraints(&specs);
        let atoms: Vec<AtomAst> = specs.into_iter().map(|(a, _)| a).collect();
        let bonds: Vec<_> = edges
            .iter()
            .map(|&(a, b)| (AtomIdx(a as u32), AtomIdx(b as u32), BondAst::from_order(1)))
            .collect();
        MoleculeAst::new(atoms, bonds, vec![], vec![], vec![], vec![], constraints)
    }

    fn solve_hmo(model: &HmoAromaticity, ast: &MoleculeAst) -> HmoOutput {
        let atoms: Vec<AtomIdx> = (0..ast.atoms().count() as u32).map(AtomIdx).collect();
        model.build_calculator(ast, &atoms).unwrap().solve()
    }

    #[fixture]
    fn hmo_model() -> HmoAromaticity {
        HmoAromaticity::new(
            ElementScope::AllowList(vec![Element::C, Element::N, Element::O, Element::S]),
            0.1,
        )
    }

    #[fixture]
    fn benzene() -> MoleculeAst {
        make_ring(vec![aromatic(Element::C, 1); 6])
    }

    #[rustfmt::skip]
    #[fixture]
    fn naphthalene() -> MoleculeAst {
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
    fn azulene() -> MoleculeAst {
        make_fused(
            vec![aromatic(Element::C, 1); 10],
            &[
                (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
                (0, 5), (5, 6), (6, 7), (7, 8), (8, 9), (9, 4),
            ],
        )
    }

    #[fixture]
    fn pyridine() -> MoleculeAst {
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
    fn pyrrole() -> MoleculeAst {
        make_ring(vec![
            aromatic(Element::N, 2),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
            aromatic(Element::C, 1),
        ])
    }

    #[fixture]
    fn cyclobutadiene() -> MoleculeAst {
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
        #[case] ast: MoleculeAst,
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
        #[case] ast: MoleculeAst,
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
        #[case] ast: MoleculeAst,
        #[case] expected_systems: usize,
        #[case] expected_atoms: Option<usize>,
    ) {
        let ring_info =
            RingEnumerator::new(RingFamily::Simple, &RingEnumerationStrategy::default())
                .enumerate(&ast);
        let systems = hmo_model.find_from_rings(&ast, &ring_info).unwrap();
        assert_eq!(systems.len(), expected_systems);
        assert_eq!(
            systems.first().map(|s| s.atom_count()),
            expected_atoms
        );
    }

    #[rstest]
    fn test_hmo_output(hmo_model: HmoAromaticity, benzene: MoleculeAst) {
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
    fn test_hmo_hamiltonian(hmo_model: HmoAromaticity, pyridine: MoleculeAst) {
        let atoms: Vec<AtomIdx> = (0..pyridine.atoms().count() as u32).map(AtomIdx).collect();
        let calc = hmo_model.build_calculator(&pyridine, &atoms).unwrap();
        let h = calc.hamiltonian();
        assert_eq!(h.nrows(), 6);
        assert_eq!(h.ncols(), 6);
        assert!(approx_eq!(f64, h[(0, 0)], 0.51, epsilon = 0.01));
    }
}
