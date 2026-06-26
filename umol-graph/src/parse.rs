//! Format parsing into a resolved [`MoleculeAst`] — parse + raise + resolve.
//!
//! These wrappers combine the io parse layer with the `ops` resolver.

use umol_ast::ast::{MoleculeAst, TryIntoAst};
use umol_io::ctfile::config::CtfileIoConfig;
use umol_io::ctfile::parser::parse_mol_bytes_to_table_ir_with;
use umol_io::smiles::config::SmilesIoConfig;
use umol_io::smiles::parser::parse_smiles_bytes_to_table_ir_with;
use umol_utils::error::UmolError;
use umol_utils::solution::Solution;

use crate::ops::model::ChemistryModel;
use crate::ops::resolve::{ResolveUnderdetermined, Resolver};

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
    let mut ast: MoleculeAst = (&table_mol).try_into_ast(&())?;
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
    let mut ast: MoleculeAst = (&table_mol).try_into_ast(&())?;
    match Resolver::new(model).resolve(&mut ast)? {
        Solution::Determined(()) => Ok(ast),
        Solution::Underdetermined(()) => Err(Box::new(ResolveUnderdetermined)),
        Solution::Contradictory(c) => Err(Box::new(c)),
    }
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;

    use rstest::*;
    use umol_ast::ast::AtomId;
    use umol_chem::element::Element;
    use umol_io::ctfile::config::CtfileIoConfig;
    use umol_io::ctfile::parse_mol_to_ast;
    use umol_io::smiles::parse_smiles_to_ast;

    use super::parse_mol_bytes_with;
    use crate::ops::model::{
        AromaticityModel, ChemistryModel, CountsModel, ElementScope, RingLimits, StereoModel,
        ValenceModel,
    };
    use crate::ops::valence::{CountsValence, ValenceTable};

    const METHANE_MOL: &str = "Methane\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    1.2345    2.3456    3.4567 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";

    const BENZENE_AROMATIC_MOL: &str = "benzene\n\n\n  6  6  0  0  0  0  0  0  0  0999 V2000\n    0.0000    1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8660   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.8660    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  4  0  0  0  0\n  2  3  4  0  0  0  0\n  3  4  4  0  0  0  0\n  4  5  4  0  0  0  0\n  5  6  4  0  0  0  0\n  6  1  4  0  0  0  0\nM  END\n";

    #[fixture]
    fn counts_model() -> CountsModel {
        CountsModel {
            table: Cow::Borrowed(ValenceTable::default_table()),
        }
    }

    #[rstest]
    #[case::methane(METHANE_MOL, 1, "C#i=#c0#h4#n0#u0#s#v0#a!")]
    #[case::benzene(BENZENE_AROMATIC_MOL, 6, "C#i=#c0#h#n0#u0#s#v2#a")]
    fn test_parse_mol_to_ast_counts_resolve(
        counts_model: CountsModel,
        #[case] input: &str,
        #[case] atom_count: u32,
        #[case] expected_atom: &str,
    ) {
        let mut ast = parse_mol_to_ast(input).unwrap();
        CountsValence::new(&counts_model).resolve(&mut ast).unwrap();
        assert_eq!(ast.atoms().count(), atom_count as usize);
        for i in 0..atom_count {
            assert_eq!(ast.atom(AtomId(i)).ast.to_string(), expected_atom);
        }
    }

    #[rstest]
    fn test_parse_smiles_to_ast_methane_counts_resolve(counts_model: CountsModel) {
        let mut ast = parse_smiles_to_ast("C").unwrap();
        CountsValence::new(&counts_model).resolve(&mut ast).unwrap();
        assert_eq!(
            ast.atom(AtomId(0)).ast.to_string(),
            "C#i=#c0#h4#n0#u0#s#v0#a!"
        );
    }

    #[rstest]
    fn test_parse_mol_bytes_with_resolver_methane_determined(counts_model: CountsModel) {
        let model = ChemistryModel {
            valence: ValenceModel::Counts(counts_model),
            aromaticity: AromaticityModel::HueckelRule {
                scope: ElementScope::AllowList(vec![Element::C]),
                ring_limits: RingLimits::default(),
            },
            stereo: StereoModel::default(),
        };
        let ast =
            parse_mol_bytes_with(METHANE_MOL.as_bytes(), &CtfileIoConfig::basic(), &model).unwrap();
        assert_eq!(
            ast.atom(AtomId(0)).ast.to_string(),
            "C#i=#c0#h4#n0#u0#s#v0#a!"
        );
    }
}
