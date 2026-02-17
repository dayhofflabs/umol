//! Resolution of TableIR molecules to GraphIR molecules.

use super::config::ResolveConfig;
use super::error::ResolutionError;
use super::molecule::MoleculeBuilder;
use super::Molecule;
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
    _molecule: &TableMolecule,
    _config: &ResolveConfig,
) -> Result<MoleculeBuilder, ResolutionError> {
    todo!()
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
