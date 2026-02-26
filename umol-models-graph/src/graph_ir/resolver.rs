//! Resolution of TableIR molecules to GraphIR molecules.

use std::collections::HashMap;

use super::config::{ResolveConfig, TopologyResolveFlags};
use super::error::ResolutionError;
use super::molecule::{AtomIndex, MoleculeBuilder};
use super::multicenter::{MulticenterBond, MulticenterContribution, MulticenterSet};
use super::{AtomBuilder, Bond, Molecule};
use crate::table_ir::Molecule as TableMolecule;

/// Resolve a TableIR molecule to a GraphIR molecule using default configuration.
pub fn resolve_molecule(molecule: &TableMolecule) -> Result<Molecule, ResolutionError> {
    resolve_molecule_with(molecule, &ResolveConfig::default())
}

/// Resolve a TableIR molecule to a GraphIR molecule with the given configuration.
///
/// Populates a `MoleculeBuilder` from the TableIR data, runs the resolution
/// phases in order, then builds the final `Molecule`.
///
/// 1. **Topology** — build graph from TableIR, validate indices, self-loops,
///    parallel edges, connectivity.
/// 2. **Valence** — match atoms against valid valence states; for aromatic atoms,
///    enumerate candidate states with σ-bond sum and π-contribution.
/// 3. **Aromaticity** — select from valence candidates by ring membership and
///    Kekulé/Hückel feasibility.
/// 4. **Stereochemistry** — validate chiral centers and bond stereo against the
///    resolved topology and bond orders.
pub fn resolve_molecule_with(
    molecule: &TableMolecule,
    config: &ResolveConfig,
) -> Result<Molecule, ResolutionError> {
    let mut builder = resolve_topology_with(molecule, config)?;
    resolve_valence_with(&mut builder, config)?;
    resolve_aromaticity_with(&mut builder, config)?;
    resolve_stereo_with(&mut builder, config)?;
    builder.build(config)
}

fn resolve_topology_with(
    molecule: &TableMolecule,
    config: &ResolveConfig,
) -> Result<MoleculeBuilder, ResolutionError> {
    if config.topology.enabled {
        if !config
            .topology
            .flags
            .contains(TopologyResolveFlags::DISCONNECTED_MOLECULES)
            && molecule.component_count() > 1
        {
            return Err(ResolutionError::TopologyDisconnected);
        }

        for (i, bond) in molecule.bonds.iter().enumerate() {
            if bond.atoms.first() == bond.atoms.second() {
                return Err(ResolutionError::TopologySelfLoop(i as u32));
            }
        }

        let mut atom_bonds: HashMap<(u32, u32), Vec<u32>> = HashMap::new();
        for (i, bond) in molecule.bonds.iter().enumerate() {
            let key = (bond.atoms.first(), bond.atoms.second());
            atom_bonds.entry(key).or_default().push(i as u32);
        }
        for (_pair, indices) in atom_bonds {
            if indices.len() >= 2 {
                return Err(ResolutionError::TopologyParallelEdges(
                    indices[0], indices[1],
                ));
            }
        }
    }

    let n = molecule.atom_count();
    let m = molecule.bond_count();
    let mut builder = MoleculeBuilder::with_capacity(n, m);
    let mut node_indices: Vec<AtomIndex> = Vec::with_capacity(n);
    for atom in &molecule.atoms {
        node_indices.push(builder.add_atom(AtomBuilder::from_table_atom(atom)));
    }
    for bond in &molecule.bonds {
        let a = bond.atoms.first();
        let b = bond.atoms.second();
        let graph_bond = Bond::from_table_bond(bond)?;
        builder
            .add_bond(
                node_indices[a as usize],
                node_indices[b as usize],
                graph_bond,
            )
            .ok_or_else(|| {
                ResolutionError::InvalidBondSpec("bond endpoint index out of range".to_string())
            })?;
    }

    for mc in &molecule.multicenter_bonds {
        let sets: Vec<MulticenterSet> = mc
            .contributions()
            .iter()
            .map(|contrib| {
                let contributions: Vec<MulticenterContribution> = contrib
                    .atoms()
                    .iter()
                    .map(|&idx| {
                        if (idx as usize) >= n {
                            Err(ResolutionError::InvalidAtomSpec(format!(
                                "multicenter bond references atom index {} out of range (0..{})",
                                idx, n
                            )))
                        } else {
                            Ok(MulticenterContribution::topology_only(
                                node_indices[idx as usize],
                            ))
                        }
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                Ok(MulticenterSet::topology_only(
                    contributions.iter().map(|c| c.atom()),
                ))
            })
            .collect::<Result<Vec<_>, _>>()?;
        builder.add_multicenter_bond(MulticenterBond::new(sets));
    }

    Ok(builder)
}

fn resolve_valence_with(
    _builder: &mut MoleculeBuilder,
    _config: &ResolveConfig,
) -> Result<(), ResolutionError> {
    todo!()
}

fn resolve_aromaticity_with(
    _builder: &mut MoleculeBuilder,
    _config: &ResolveConfig,
) -> Result<(), ResolutionError> {
    todo!()
}

fn resolve_stereo_with(
    _builder: &mut MoleculeBuilder,
    _config: &ResolveConfig,
) -> Result<(), ResolutionError> {
    todo!()
}
