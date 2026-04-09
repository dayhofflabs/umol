//! Kekulization: assignment of definite single/double bond orders to aromatic bonds.
//!
//! This module is independent of the resolution pipeline. It operates on
//! validated `AromaticSystem` objects and is invoked explicitly by the caller.
//! The algorithm uses backtracking DFS with an optional bound on backtrack steps.

use std::cmp::Ordering;
use std::collections::{BTreeMap, HashMap, HashSet};

use thiserror::Error;

use super::aromaticity::AromaticSystem;
use super::molecule::{AtomIndex, BondIndex};
use super::molecule_builder::MoleculeBuilder;

#[derive(Debug, Error, Clone, PartialEq, Eq)]
pub enum KekulizationError {
    #[error("no valid kekulization assignment found")]
    NoAssignment,
    #[error("kekulization backtrack limit exceeded ({limit})")]
    BacktrackLimitExceeded { limit: usize },
    #[error("unsupported kekulization algorithm: {0}")]
    UnsupportedAlgorithm(String),
    #[error("invalid aromatic system: {0}")]
    InvalidSystem(String),
}

/// Result of Kekulization for one aromatic system.
#[derive(Clone, Debug)]
pub struct KekuleAssignment {
    pub bond_orders: BTreeMap<BondIndex, u8>,
}

/// Configuration for the Kekulizer.
#[derive(Clone, Debug)]
pub struct KekuleConfig {
    pub max_backtrack_steps: usize,
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
) -> Result<KekuleAssignment, KekulizationError> {
    let aromatic_atoms: HashSet<AtomIndex> = system.atoms().collect();
    if aromatic_atoms.is_empty() {
        return Ok(KekuleAssignment {
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
            hb.partial_cmp(&ha).unwrap_or(Ordering::Equal)
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
        Ok(KekuleAssignment {
            bond_orders: result,
        })
    } else {
        if backtracks >= config.max_backtrack_steps {
            Err(KekulizationError::BacktrackLimitExceeded {
                limit: config.max_backtrack_steps,
            })
        } else {
            Err(KekulizationError::NoAssignment)
        }
    }
}

#[allow(clippy::too_many_arguments)]
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
