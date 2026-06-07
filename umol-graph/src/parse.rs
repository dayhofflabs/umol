//! Format parsing into a resolved [`MoleculeAst`] — parse + raise + resolve.
//!
//! These wrappers combine the io parse layer with the `ops` resolver. Per doc 065 a
//! cross-module combiner returns `Box<dyn UmolError>`: a parse `ParseError`, a
//! `ResolverContradiction`, or `ResolveUnderdetermined`, each boxed at this boundary.

use umol_ast::ast::{MoleculeAst, TryIntoAst};
use umol_shared::error::UmolError;

use crate::io::ctfile::config::CtfileIoConfig;
use crate::io::ctfile::parser::parse_mol_bytes_to_table_ir_with;
use crate::io::smiles::config::SmilesIoConfig;
use crate::io::smiles::parser::parse_smiles_bytes_to_table_ir_with;
use crate::ops::model::ChemistryModel;
use crate::ops::resolver::{ResolveUnderdetermined, Resolver};
use crate::ops::solution::Solution;

/// Parse SMILES to a resolved [`MoleculeAst`] using default IO config and model.
pub fn parse_smiles(input: &str) -> Result<MoleculeAst, Box<dyn UmolError>> {
    parse_smiles_bytes(input.as_bytes())
}

/// Parse SMILES bytes to a resolved [`MoleculeAst`] using default IO config and model.
pub fn parse_smiles_bytes(input: &[u8]) -> Result<MoleculeAst, Box<dyn UmolError>> {
    parse_smiles_bytes_with(
        input,
        &SmilesIoConfig::basic_opensmiles(),
        &ChemistryModel::default(),
    )
}

/// Parse SMILES to a resolved [`MoleculeAst`] with explicit IO config and model.
pub fn parse_smiles_with(
    input: &str,
    io_config: &SmilesIoConfig,
    model: &ChemistryModel,
) -> Result<MoleculeAst, Box<dyn UmolError>> {
    parse_smiles_bytes_with(input.as_bytes(), io_config, model)
}

/// Parse SMILES bytes to a resolved [`MoleculeAst`] with explicit IO config and model.
pub fn parse_smiles_bytes_with(
    input: &[u8],
    io_config: &SmilesIoConfig,
    model: &ChemistryModel,
) -> Result<MoleculeAst, Box<dyn UmolError>> {
    let table_mol = parse_smiles_bytes_to_table_ir_with(input, io_config)?;
    let mut ast: MoleculeAst = (&table_mol)
        .try_into_ast(&())
        .expect("table_ir → MoleculeAst raise is currently infallible");
    match Resolver::new(model).resolve(&mut ast)? {
        Solution::Determined(()) => Ok(ast),
        Solution::Underdetermined(()) => Err(Box::new(ResolveUnderdetermined)),
        Solution::Contradictory(c) => Err(Box::new(c)),
    }
}

/// Parse MOL to a resolved [`MoleculeAst`] using default IO config and model.
pub fn parse_mol(input: &str) -> Result<MoleculeAst, Box<dyn UmolError>> {
    parse_mol_bytes(input.as_bytes())
}

/// Parse MOL bytes to a resolved [`MoleculeAst`] using default IO config and model.
pub fn parse_mol_bytes(input: &[u8]) -> Result<MoleculeAst, Box<dyn UmolError>> {
    parse_mol_bytes_with(input, &CtfileIoConfig::basic(), &ChemistryModel::default())
}

/// Parse MOL to a resolved [`MoleculeAst`] with explicit IO config and model.
pub fn parse_mol_with(
    input: &str,
    io_config: &CtfileIoConfig,
    model: &ChemistryModel,
) -> Result<MoleculeAst, Box<dyn UmolError>> {
    parse_mol_bytes_with(input.as_bytes(), io_config, model)
}

/// Parse MOL bytes to a resolved [`MoleculeAst`] with explicit IO config and model.
pub fn parse_mol_bytes_with(
    input: &[u8],
    io_config: &CtfileIoConfig,
    model: &ChemistryModel,
) -> Result<MoleculeAst, Box<dyn UmolError>> {
    let table_mol = parse_mol_bytes_to_table_ir_with(input, io_config)?;
    let mut ast: MoleculeAst = (&table_mol)
        .try_into_ast(&())
        .expect("table_ir → MoleculeAst raise is currently infallible");
    match Resolver::new(model).resolve(&mut ast)? {
        Solution::Determined(()) => Ok(ast),
        Solution::Underdetermined(()) => Err(Box::new(ResolveUnderdetermined)),
        Solution::Contradictory(c) => Err(Box::new(c)),
    }
}
