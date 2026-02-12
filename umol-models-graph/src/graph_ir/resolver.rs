//! Resolution of TableIR molecules to GraphIR molecules.

use super::error::ResolutionError;
use super::Molecule;
use crate::table_ir::Molecule as TableMolecule;

/// Resolve a TableIR molecule to a GraphIR molecule.
pub fn resolve_molecule(_table_molecule: &TableMolecule) -> Result<Molecule, ResolutionError> {
    todo!()
}
