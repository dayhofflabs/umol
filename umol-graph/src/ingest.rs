//! Conversion of parsed external-format values into graph models.

use std::any::Any;

use thiserror::Error;
use umol_ast::ast::{MoleculeAst, TryIntoAst};
use umol_io::smiles::{ParseError as SmilesParseError, Smiles};
use umol_io::table_ir::raise::RaiseError;
use umol_utils::error::UmolError;
use umol_utils::solution::Solution;

use crate::ops::model::ChemistryModel;
use crate::ops::resolve::{ResolveUnderdetermined, Resolver, ResolverContradiction, ResolverError};

/// Convert a parsed external-format value into a graph model.
pub trait Ingest {
    type Output;
    type Error;

    /// Interpret this format value under `model`.
    fn ingest(&self, model: &ChemistryModel) -> Result<Self::Output, Self::Error>;
}

/// Failure while ingesting a parsed molecular representation.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum MoleculeIngestError {
    #[error("{0}")]
    ModelConversion(#[from] RaiseError),
    #[error("{0}")]
    Contradiction(#[from] ResolverContradiction),
    #[error("{0}")]
    Underdetermined(#[from] ResolveUnderdetermined),
    #[error("{0}")]
    Execution(#[from] ResolverError),
}

impl UmolError for MoleculeIngestError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

/// Failure while accepting SMILES text as a determined molecule.
#[derive(Clone, Debug, PartialEq, Error)]
pub enum SmilesInputError {
    #[error("{0}")]
    Syntax(#[from] SmilesParseError),
    #[error("{0}")]
    ModelConversion(#[from] RaiseError),
    #[error("{0}")]
    Contradiction(#[from] ResolverContradiction),
    #[error("{0}")]
    Underdetermined(#[from] ResolveUnderdetermined),
    #[error("{0}")]
    Execution(#[from] ResolverError),
}

impl From<MoleculeIngestError> for SmilesInputError {
    fn from(error: MoleculeIngestError) -> Self {
        match error {
            MoleculeIngestError::ModelConversion(error) => Self::ModelConversion(error),
            MoleculeIngestError::Contradiction(error) => Self::Contradiction(error),
            MoleculeIngestError::Underdetermined(error) => Self::Underdetermined(error),
            MoleculeIngestError::Execution(error) => Self::Execution(error),
        }
    }
}

impl UmolError for SmilesInputError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl Ingest for Smiles {
    type Output = MoleculeAst;
    type Error = MoleculeIngestError;

    fn ingest(&self, model: &ChemistryModel) -> Result<Self::Output, Self::Error> {
        let mut ast: MoleculeAst = self.as_table_ir().try_into_ast(&())?;
        match Resolver::new(model).resolve(&mut ast)? {
            Solution::Determined(()) => Ok(ast),
            Solution::Underdetermined(()) => Err(ResolveUnderdetermined.into()),
            Solution::Contradictory(error) => Err(error.into()),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::error::Error as _;

    use rstest::rstest;
    use umol_ast::mol_dsl;

    use super::*;
    use crate::ops::aromaticity::{AromaticityContradiction, AromaticityError};
    use crate::ops::model::AromaticityModel;

    #[rstest]
    #[case::model_conversion(
        MoleculeIngestError::ModelConversion(RaiseError::WedgeConflict { atom: 2 }),
        "inconsistent wedge bonds at atom 2"
    )]
    #[case::contradiction(
        MoleculeIngestError::Contradiction(ResolverContradiction::Aromaticity(
            AromaticityContradiction::HmoInvalidInput(String::from("invalid input")),
        )),
        "hmo: invalid input: invalid input"
    )]
    #[case::underdetermined(
        MoleculeIngestError::Underdetermined(ResolveUnderdetermined),
        "resolution underdetermined"
    )]
    #[case::execution(
        MoleculeIngestError::Execution(ResolverError::Aromaticity(
            AromaticityError::HmoMissingParameters(String::from("carbon")),
        )),
        "hmo: missing parameters: carbon"
    )]
    fn test_molecule_ingest_error(#[case] error: MoleculeIngestError, #[case] expected: &str) {
        assert_eq!(error.to_string(), expected);
        assert_eq!(
            error.source().map(ToString::to_string).as_deref(),
            Some(expected)
        );
    }

    #[rstest]
    #[case::syntax(
        SmilesInputError::Syntax(SmilesParseError::LeadingWhitespace),
        "Leading whitespace"
    )]
    #[case::model_conversion(
        SmilesInputError::ModelConversion(RaiseError::WedgeConflict { atom: 2 }),
        "inconsistent wedge bonds at atom 2"
    )]
    #[case::contradiction(
        SmilesInputError::Contradiction(ResolverContradiction::Aromaticity(
            AromaticityContradiction::HmoInvalidInput(String::from("invalid input")),
        )),
        "hmo: invalid input: invalid input"
    )]
    #[case::underdetermined(
        SmilesInputError::Underdetermined(ResolveUnderdetermined),
        "resolution underdetermined"
    )]
    #[case::execution(
        SmilesInputError::Execution(ResolverError::Aromaticity(
            AromaticityError::HmoMissingParameters(String::from("carbon")),
        )),
        "hmo: missing parameters: carbon"
    )]
    fn test_smiles_input_error(#[case] error: SmilesInputError, #[case] expected: &str) {
        assert_eq!(error.to_string(), expected);
        assert_eq!(
            error.source().map(ToString::to_string).as_deref(),
            Some(expected)
        );
    }

    #[rstest]
    #[case::model_conversion(
        MoleculeIngestError::ModelConversion(RaiseError::WedgeConflict { atom: 2 }),
        SmilesInputError::ModelConversion(RaiseError::WedgeConflict { atom: 2 })
    )]
    #[case::contradiction(
        MoleculeIngestError::Contradiction(ResolverContradiction::Aromaticity(
            AromaticityContradiction::HmoInvalidInput(String::from("invalid input")),
        )),
        SmilesInputError::Contradiction(ResolverContradiction::Aromaticity(
            AromaticityContradiction::HmoInvalidInput(String::from("invalid input")),
        ))
    )]
    #[case::underdetermined(
        MoleculeIngestError::Underdetermined(ResolveUnderdetermined),
        SmilesInputError::Underdetermined(ResolveUnderdetermined)
    )]
    #[case::execution(
        MoleculeIngestError::Execution(ResolverError::Aromaticity(
            AromaticityError::HmoMissingParameters(String::from("carbon")),
        )),
        SmilesInputError::Execution(ResolverError::Aromaticity(
            AromaticityError::HmoMissingParameters(String::from("carbon")),
        ))
    )]
    fn test_smiles_input_error_from(
        #[case] input: MoleculeIngestError,
        #[case] expected: SmilesInputError,
    ) {
        assert_eq!(SmilesInputError::from(input), expected);
    }

    #[rstest]
    fn test_smiles_ingest() {
        let smiles = Smiles::parse("C").unwrap();
        let expected = mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"]}"#);
        let default_model = ChemistryModel::default();
        let permissive_model = ChemistryModel {
            aromaticity: AromaticityModel::permissive(),
            ..ChemistryModel::default()
        };

        assert_eq!(smiles.ingest(&default_model), Ok(expected.clone()));
        assert_eq!(smiles.ingest(&permissive_model), Ok(expected));
    }
}
