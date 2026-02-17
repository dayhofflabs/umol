//! Resolution of TableIR molecules to GraphIR molecules.

use super::config::ResolveConfig;
use super::error::ResolutionError;
use super::Molecule;
use crate::table_ir::Molecule as TableMolecule;

/// Resolve a TableIR molecule to a GraphIR molecule using default configuration.
pub fn resolve_molecule(molecule: &TableMolecule) -> Result<Molecule, ResolutionError> {
    resolve_molecule_with(molecule, &ResolveConfig::default())
}

/// Resolve a TableIR molecule to a GraphIR molecule with the given configuration.
///
/// Resolution proceeds in four phases, each narrowing the set of valid
/// interpretations without backtracking:
///
/// 1. **Topology** — index validity, self-loops, parallel edges, connectivity.
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
    let resolved = resolve_topology_with(molecule, config)?;
    let resolved = resolve_valence_with(resolved, config)?;
    let resolved = resolve_aromaticity_with(resolved, config)?;
    let resolved = resolve_stereo_with(resolved, config)?;
    Ok(resolved)
}

fn resolve_topology_with(
    _molecule: &TableMolecule,
    _config: &ResolveConfig,
) -> Result<Molecule, ResolutionError> {
    todo!()
}

fn resolve_valence_with(
    _molecule: Molecule,
    _config: &ResolveConfig,
) -> Result<Molecule, ResolutionError> {
    todo!()
}

fn resolve_aromaticity_with(
    _molecule: Molecule,
    _config: &ResolveConfig,
) -> Result<Molecule, ResolutionError> {
    todo!()
}

fn resolve_stereo_with(
    _molecule: Molecule,
    _config: &ResolveConfig,
) -> Result<Molecule, ResolutionError> {
    todo!()
}
