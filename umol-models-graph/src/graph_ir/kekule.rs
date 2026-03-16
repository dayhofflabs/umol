//! Kekulization: assignment of definite single/double bond orders to aromatic bonds.
//!
//! This module is independent of the resolution pipeline. It operates on
//! validated `AromaticSystem` objects and is invoked explicitly by the caller.
//! The algorithm uses backtracking DFS with an optional bound on backtrack steps.

use std::collections::{BTreeMap, HashMap, HashSet};

use crate::graph_ir::aromaticity::AromaticSystem;
use crate::graph_ir::molecule::{AtomIndex, BondIndex, MoleculeBuilder};

/// Result of Kekulization for one aromatic system.
#[derive(Clone, Debug)]
pub struct KekuleAssignment {
    /// Bond order for each bond in the aromatic system.
    /// Key: bond index. Value: 1 (single) or 2 (double).
    pub bond_orders: BTreeMap<BondIndex, u8>,
}

/// Configuration for the Kekulizer.
#[derive(Clone, Debug)]
pub struct KekuleConfig {
    /// Maximum number of backtrack steps before giving up.
    pub max_backtrack_steps: usize,
    /// Optional HMO bond order hints: bonds with higher pi-bond-order are
    /// tried first as double bonds.
    pub bond_order_hints: Option<BTreeMap<(AtomIndex, AtomIndex), f64>>,
}

impl Default for KekuleConfig {
    fn default() -> Self {
        Self {
            max_backtrack_steps: 100_000,
            bond_order_hints: None,
        }
    }
}

/// Attempt to find a valid Kekulé structure for an aromatic system.
///
/// A valid Kekulé structure assigns single/double bond orders such that
/// every atom in the aromatic system that contributes an odd number of
/// pi-electrons is incident to exactly one double bond within the system,
/// and atoms contributing an even number of pi-electrons are incident to
/// zero double bonds within the system (their electrons come from lone pairs).
pub fn kekulize(
    builder: &MoleculeBuilder,
    system: &AromaticSystem,
    config: &KekuleConfig,
) -> Option<KekuleAssignment> {
    let aromatic_atoms: HashSet<AtomIndex> = system.atoms().collect();
    if aromatic_atoms.is_empty() {
        return Some(KekuleAssignment {
            bond_orders: BTreeMap::new(),
        });
    }

    // Collect aromatic bonds (both endpoints in the aromatic system).
    let mut aromatic_bonds: Vec<(BondIndex, AtomIndex, AtomIndex)> = Vec::new();
    let mut atom_bonds: HashMap<AtomIndex, Vec<usize>> = HashMap::new();
    for &atom in &aromatic_atoms {
        atom_bonds.entry(atom).or_default();
    }

    for bond_idx in builder.bond_indices() {
        if let Some((a, b)) = builder.bond_atom_indices(bond_idx) {
            if aromatic_atoms.contains(&a) && aromatic_atoms.contains(&b) {
                let idx = aromatic_bonds.len();
                aromatic_bonds.push((bond_idx, a, b));
                atom_bonds.entry(a).or_default().push(idx);
                atom_bonds.entry(b).or_default().push(idx);
            }
        }
    }

    // Determine the "demand" for each atom: how many double bonds it needs.
    // Atoms contributing 1 pi-electron need exactly 1 double bond.
    // Atoms contributing 2 pi-electrons need 0 double bonds (lone pair donors).
    // Atoms contributing 0 pi-electrons need 0 double bonds (e.g., tropylium C+).
    let mut demand: HashMap<AtomIndex, u8> = HashMap::new();
    for contrib in system.contributions() {
        let d = if contrib.aromatic_valence() == 1 {
            1
        } else {
            0
        };
        demand.insert(contrib.atom(), d);
    }

    // Sort bonds by HMO hint (highest pi-bond-order first) for better search order.
    let mut bond_order: Vec<usize> = (0..aromatic_bonds.len()).collect();
    if let Some(hints) = &config.bond_order_hints {
        bond_order.sort_by(|&a, &b| {
            let (_, a1, a2) = aromatic_bonds[a];
            let (_, b1, b2) = aromatic_bonds[b];
            let key_a = if a1 < a2 { (a1, a2) } else { (a2, a1) };
            let key_b = if b1 < b2 { (b1, b2) } else { (b2, b1) };
            let ha = hints.get(&key_a).copied().unwrap_or(0.0);
            let hb = hints.get(&key_b).copied().unwrap_or(0.0);
            hb.partial_cmp(&ha).unwrap_or(std::cmp::Ordering::Equal)
        });
    }

    // Backtracking DFS.
    let mut assignment: Vec<u8> = vec![0; aromatic_bonds.len()]; // 0 = unassigned
    let mut current_demand: HashMap<AtomIndex, u8> = demand.clone();
    let mut backtracks: usize = 0;

    if dfs_kekulize(
        0,
        &bond_order,
        &aromatic_bonds,
        &atom_bonds,
        &mut assignment,
        &mut current_demand,
        &mut backtracks,
        config.max_backtrack_steps,
    ) {
        let mut result = BTreeMap::new();
        for (i, &order) in assignment.iter().enumerate() {
            let (bond_idx, _, _) = aromatic_bonds[i];
            result.insert(bond_idx, order);
        }
        Some(KekuleAssignment {
            bond_orders: result,
        })
    } else {
        None
    }
}

fn dfs_kekulize(
    depth: usize,
    bond_order: &[usize],
    aromatic_bonds: &[(BondIndex, AtomIndex, AtomIndex)],
    atom_bonds: &HashMap<AtomIndex, Vec<usize>>,
    assignment: &mut [u8],
    current_demand: &mut HashMap<AtomIndex, u8>,
    backtracks: &mut usize,
    max_backtracks: usize,
) -> bool {
    if *backtracks >= max_backtracks {
        return false;
    }

    if depth >= bond_order.len() {
        return current_demand.values().all(|&d| d == 0);
    }

    let bond_idx = bond_order[depth];
    let (_, a, b) = aromatic_bonds[bond_idx];
    let da = *current_demand.get(&a).unwrap_or(&0);
    let db = *current_demand.get(&b).unwrap_or(&0);

    // Try double bond (if both endpoints still need one).
    if da > 0 && db > 0 {
        assignment[bond_idx] = 2;
        *current_demand.get_mut(&a).unwrap() -= 1;
        *current_demand.get_mut(&b).unwrap() -= 1;
        if dfs_kekulize(
            depth + 1,
            bond_order,
            aromatic_bonds,
            atom_bonds,
            assignment,
            current_demand,
            backtracks,
            max_backtracks,
        ) {
            return true;
        }
        *current_demand.get_mut(&a).unwrap() += 1;
        *current_demand.get_mut(&b).unwrap() += 1;
        *backtracks += 1;
    }

    // Try single bond.
    assignment[bond_idx] = 1;
    if dfs_kekulize(
        depth + 1,
        bond_order,
        aromatic_bonds,
        atom_bonds,
        assignment,
        current_demand,
        backtracks,
        max_backtracks,
    ) {
        return true;
    }
    *backtracks += 1;

    assignment[bond_idx] = 0;
    false
}

#[cfg(test)]
mod tests {
    use smallvec::SmallVec;
    use umol_data::Element;

    use super::*;
    use crate::graph_ir::aromaticity::AromaticContribution;
    use crate::graph_ir::atom::AtomBuilder;
    use crate::graph_ir::bond::BondBuilder;

    fn carbon_aromatic_1() -> AtomBuilder {
        let spec = crate::spec!("{Cv2a1H}");
        let mut ab = AtomBuilder::new(Element::C);
        ab.set_candidates(SmallVec::from_elem(spec, 1));
        ab
    }

    fn make_benzene_system() -> (MoleculeBuilder, AromaticSystem) {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..6)
            .map(|_| builder.add_atom(carbon_aromatic_1()))
            .collect();
        for i in 0..6 {
            builder.add_bond_unchecked(atoms[i], atoms[(i + 1) % 6], BondBuilder::new(1, None));
        }
        let system = AromaticSystem::new(atoms.iter().map(|&a| AromaticContribution::new(a, 1)));
        (builder, system)
    }

    fn make_naphthalene_system() -> (MoleculeBuilder, AromaticSystem) {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<AtomIndex> = (0..10)
            .map(|_| builder.add_atom(carbon_aromatic_1()))
            .collect();
        let r1 = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)];
        for (a, b) in r1 {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondBuilder::new(1, None));
        }
        let r2 = [(3, 6), (6, 7), (7, 8), (8, 9), (9, 4)];
        for (a, b) in r2 {
            builder.add_bond_unchecked(atoms[a], atoms[b], BondBuilder::new(1, None));
        }
        let system = AromaticSystem::new(atoms.iter().map(|&a| AromaticContribution::new(a, 1)));
        (builder, system)
    }

    #[test]
    fn kekulize_benzene() {
        let (builder, system) = make_benzene_system();
        let config = KekuleConfig::default();
        let result = kekulize(&builder, &system, &config);
        assert!(result.is_some());
        let assignment = result.unwrap();
        assert_eq!(assignment.bond_orders.len(), 6);
        let doubles: usize = assignment.bond_orders.values().filter(|&&v| v == 2).count();
        let singles: usize = assignment.bond_orders.values().filter(|&&v| v == 1).count();
        assert_eq!(doubles, 3);
        assert_eq!(singles, 3);

        // Verify each atom is incident to exactly one double bond.
        let atoms: Vec<AtomIndex> = builder.atom_indices().collect();
        for &atom in &atoms {
            let mut double_count = 0;
            for bond_idx in builder.atom_bond_indices(atom) {
                if let Some(&order) = assignment.bond_orders.get(&bond_idx) {
                    if order == 2 {
                        double_count += 1;
                    }
                }
            }
            assert_eq!(
                double_count, 1,
                "atom {:?} should have exactly 1 double bond",
                atom
            );
        }
    }

    #[test]
    fn kekulize_naphthalene() {
        let (builder, system) = make_naphthalene_system();
        let config = KekuleConfig::default();
        let result = kekulize(&builder, &system, &config);
        assert!(result.is_some());
        let assignment = result.unwrap();
        let doubles: usize = assignment.bond_orders.values().filter(|&&v| v == 2).count();
        assert_eq!(doubles, 5); // Naphthalene has 5 double bonds in Kekule form.
    }
}
