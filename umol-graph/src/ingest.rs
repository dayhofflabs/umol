//! Conversion of parsed external-format values into graph models.

use std::any::Any;

use thiserror::Error;
use umol_ast::ast::{MoleculeAst, TryIntoAst};
use umol_io::smiles::{ParseError as SmilesParseError, Smiles, SmilesIoConfig};
use umol_io::table_ir::raise::RaiseError;
use umol_utils::error::UmolError;
use umol_utils::solution::Solution;

use crate::ops::model::ChemistryModel;
use crate::ops::resolve::{
    ResolveConfig, ResolveUnderdetermined, Resolver, ResolverContradiction, ResolverError,
};

/// Convert a parsed external-format value into a graph model.
pub trait Interpret {
    type Output;
    type Error;

    /// Interpret this format value under the semantic model and resolve policy.
    fn interpret(
        &self,
        model: &ChemistryModel,
        resolve_config: &ResolveConfig,
    ) -> Result<Self::Output, Self::Error>;
}

/// Failure while interpreting a parsed molecular representation.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum MoleculeInterpretationError {
    #[error("{0}")]
    ModelConversion(#[from] RaiseError),
    #[error("{0}")]
    Contradiction(#[from] ResolverContradiction),
    #[error("{0}")]
    Underdetermined(#[from] ResolveUnderdetermined),
    #[error("{0}")]
    Execution(#[from] ResolverError),
}

impl UmolError for MoleculeInterpretationError {
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

impl From<MoleculeInterpretationError> for SmilesInputError {
    fn from(error: MoleculeInterpretationError) -> Self {
        match error {
            MoleculeInterpretationError::ModelConversion(error) => Self::ModelConversion(error),
            MoleculeInterpretationError::Contradiction(error) => Self::Contradiction(error),
            MoleculeInterpretationError::Underdetermined(error) => Self::Underdetermined(error),
            MoleculeInterpretationError::Execution(error) => Self::Execution(error),
        }
    }
}

impl UmolError for SmilesInputError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl Interpret for Smiles {
    type Output = MoleculeAst;
    type Error = MoleculeInterpretationError;

    fn interpret(
        &self,
        model: &ChemistryModel,
        resolve_config: &ResolveConfig,
    ) -> Result<Self::Output, Self::Error> {
        let mut ast: MoleculeAst = self.as_table_ir().try_into_ast(&())?;
        match Resolver::with_config(model, *resolve_config).resolve(&mut ast)? {
            Solution::Determined(()) => Ok(ast),
            Solution::Underdetermined(()) => Err(ResolveUnderdetermined.into()),
            Solution::Contradictory(error) => Err(error.into()),
        }
    }
}

/// Ingest SMILES text with the OpenSMILES configuration and default model.
pub fn ingest_smiles(input: &str) -> Result<MoleculeAst, SmilesInputError> {
    ingest_smiles_bytes(input.as_bytes())
}

/// Ingest SMILES bytes with the OpenSMILES configuration and default model.
pub fn ingest_smiles_bytes(input: &[u8]) -> Result<MoleculeAst, SmilesInputError> {
    ingest_smiles_bytes_with(
        input,
        &SmilesIoConfig::opensmiles(),
        &ChemistryModel::default(),
        &ResolveConfig::default(),
    )
}

/// Ingest SMILES text with explicit IO, chemistry, and resolve configuration.
pub fn ingest_smiles_with(
    input: &str,
    io_config: &SmilesIoConfig,
    model: &ChemistryModel,
    resolve_config: &ResolveConfig,
) -> Result<MoleculeAst, SmilesInputError> {
    ingest_smiles_bytes_with(input.as_bytes(), io_config, model, resolve_config)
}

/// Ingest SMILES bytes with explicit IO, chemistry, and resolve configuration.
pub fn ingest_smiles_bytes_with(
    input: &[u8],
    io_config: &SmilesIoConfig,
    model: &ChemistryModel,
    resolve_config: &ResolveConfig,
) -> Result<MoleculeAst, SmilesInputError> {
    let smiles = Smiles::parse_bytes_with(input, io_config)?;
    smiles
        .interpret(model, resolve_config)
        .map_err(SmilesInputError::from)
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;
    use std::error::Error as _;

    use rstest::rstest;
    use umol_ast::ast::{
        AromaticSystemId, AromaticValenceAst, AtomId, StereoCosetAst, TetrahedralStereoAst,
        ValueAst,
    };
    use umol_ast::{atom_dsl, mol_dsl};

    use super::*;
    use crate::ops::aromaticity::{AromaticityContradiction, AromaticityError};
    use crate::ops::model::{AromaticityModel, ElementScope, RingLimits, ValenceModel};
    use crate::ops::resolve::{AromaticityResolveConfig, StereoResolveConfig};
    use crate::ops::valence::AtomTypeRegistry;

    #[rstest]
    #[case::model_conversion(
        MoleculeInterpretationError::ModelConversion(RaiseError::WedgeConflict { atom: 2 }),
        "inconsistent wedge bonds at atom 2"
    )]
    #[case::contradiction(
        MoleculeInterpretationError::Contradiction(ResolverContradiction::Aromaticity(
            AromaticityContradiction::HmoInvalidInput(String::from("invalid input")),
        )),
        "hmo: invalid input: invalid input"
    )]
    #[case::underdetermined(
        MoleculeInterpretationError::Underdetermined(ResolveUnderdetermined),
        "resolution underdetermined"
    )]
    #[case::execution(
        MoleculeInterpretationError::Execution(ResolverError::Aromaticity(
            AromaticityError::HmoMissingParameters(String::from("carbon")),
        )),
        "hmo: missing parameters: carbon"
    )]
    fn test_molecule_interpretation_error(
        #[case] error: MoleculeInterpretationError,
        #[case] expected: &str,
    ) {
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
        MoleculeInterpretationError::ModelConversion(RaiseError::WedgeConflict { atom: 2 }),
        SmilesInputError::ModelConversion(RaiseError::WedgeConflict { atom: 2 })
    )]
    #[case::contradiction(
        MoleculeInterpretationError::Contradiction(ResolverContradiction::Aromaticity(
            AromaticityContradiction::HmoInvalidInput(String::from("invalid input")),
        )),
        SmilesInputError::Contradiction(ResolverContradiction::Aromaticity(
            AromaticityContradiction::HmoInvalidInput(String::from("invalid input")),
        ))
    )]
    #[case::underdetermined(
        MoleculeInterpretationError::Underdetermined(ResolveUnderdetermined),
        SmilesInputError::Underdetermined(ResolveUnderdetermined)
    )]
    #[case::execution(
        MoleculeInterpretationError::Execution(ResolverError::Aromaticity(
            AromaticityError::HmoMissingParameters(String::from("carbon")),
        )),
        SmilesInputError::Execution(ResolverError::Aromaticity(
            AromaticityError::HmoMissingParameters(String::from("carbon")),
        ))
    )]
    fn test_smiles_input_error_from(
        #[case] input: MoleculeInterpretationError,
        #[case] expected: SmilesInputError,
    ) {
        assert_eq!(SmilesInputError::from(input), expected);
    }

    #[rstest]
    #[case::default(ChemistryModel::default(), ResolveConfig::default())]
    #[case::permissive(
        ChemistryModel {
            aromaticity: AromaticityModel::permissive(),
            ..ChemistryModel::default()
        },
        ResolveConfig::default()
    )]
    fn test_smiles_interpret(#[case] model: ChemistryModel, #[case] resolve_config: ResolveConfig) {
        let smiles = Smiles::parse("C").unwrap();
        let expected = mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"]}"#);

        assert_eq!(smiles.interpret(&model, &resolve_config), Ok(expected));
    }

    #[rstest]
    #[case::wildcard("*", ChemistryModel::default(), ResolveConfig::default())]
    fn test_smiles_interpret_error(
        #[case] input: &str,
        #[case] model: ChemistryModel,
        #[case] resolve_config: ResolveConfig,
    ) {
        let parsed = Smiles::parse(input).unwrap();

        assert_eq!(
            parsed.interpret(&model, &resolve_config),
            Err(MoleculeInterpretationError::Underdetermined(
                ResolveUnderdetermined
            ))
        );
    }

    #[rstest]
    #[case::methane("C")]
    #[case::benzene("c1ccccc1")]
    fn test_ingest_smiles(#[case] input: &str) {
        assert_eq!(
            ingest_smiles(input),
            ingest_smiles_with(
                input,
                &SmilesIoConfig::opensmiles(),
                &ChemistryModel::default(),
                &ResolveConfig::default(),
            )
        );
    }

    #[rstest]
    #[case::syntax(" C", SmilesInputError::Syntax(SmilesParseError::LeadingWhitespace))]
    #[case::model_conversion(
        "C[S@]C",
        SmilesInputError::ModelConversion(RaiseError::TetrahedralLigandCount {
            atom: 1,
            count: 2,
        })
    )]
    #[case::underdetermined("*", SmilesInputError::Underdetermined(ResolveUnderdetermined))]
    fn test_ingest_smiles_error(#[case] input: &str, #[case] expected: SmilesInputError) {
        assert_eq!(ingest_smiles(input), Err(expected));
    }

    #[rstest]
    #[case::methane(b"C")]
    #[case::benzene(b"c1ccccc1")]
    fn test_ingest_smiles_bytes(#[case] input: &[u8]) {
        assert_eq!(
            ingest_smiles_bytes(input),
            ingest_smiles_bytes_with(
                input,
                &SmilesIoConfig::opensmiles(),
                &ChemistryModel::default(),
                &ResolveConfig::default(),
            )
        );
    }

    #[rstest]
    #[case::delocalized(
        SmilesIoConfig::opensmiles(),
        ChemistryModel::default(),
        ResolveConfig::default(),
        vec![ValueAst::Lit(0); 3],
        ValueAst::Lit(1)
    )]
    #[case::localized(
        SmilesIoConfig::opensmiles(),
        ChemistryModel::default(),
        ResolveConfig {
            aromaticity: AromaticityResolveConfig {
                perception: Default::default(),
                delocalize_charge: false,
                reset_aromatic_valence: false,
            },
            stereo: StereoResolveConfig::default(),
        },
        vec![ValueAst::Lit(1), ValueAst::Lit(0), ValueAst::Lit(0)],
        ValueAst::Lit(0)
    )]
    fn test_ingest_smiles_with_charge(
        #[case] io_config: SmilesIoConfig,
        #[case] model: ChemistryModel,
        #[case] resolve_config: ResolveConfig,
        #[case] expected_atom_charges: Vec<ValueAst>,
        #[case] expected_system_charge: ValueAst,
    ) {
        let ast =
            ingest_smiles_with("[cH+]1[cH][cH]1", &io_config, &model, &resolve_config).unwrap();

        assert_eq!(
            ast.atoms()
                .iter()
                .map(|atom| atom.ast.charge.clone())
                .collect::<Vec<_>>(),
            expected_atom_charges
        );
        assert_eq!(
            ast.aromatic_system(AromaticSystemId(0)).ast.charge,
            expected_system_charge
        );
    }

    #[rstest]
    #[case::retained(
        SmilesIoConfig::opensmiles(),
        ChemistryModel::default(),
        ResolveConfig::default(),
        vec![Some(AromaticValenceAst::Aromatic(ValueAst::Lit(1))); 6]
    )]
    #[case::reset(
        SmilesIoConfig::opensmiles(),
        ChemistryModel::default(),
        ResolveConfig {
            aromaticity: AromaticityResolveConfig {
                perception: Default::default(),
                delocalize_charge: true,
                reset_aromatic_valence: true,
            },
            stereo: StereoResolveConfig::default(),
        },
        vec![None; 6]
    )]
    fn test_ingest_smiles_with_aromatic_valence(
        #[case] io_config: SmilesIoConfig,
        #[case] model: ChemistryModel,
        #[case] resolve_config: ResolveConfig,
        #[case] expected: Vec<Option<AromaticValenceAst>>,
    ) {
        let ast = ingest_smiles_with("c1ccccc1", &io_config, &model, &resolve_config).unwrap();

        assert_eq!(
            ast.atoms()
                .iter()
                .map(|atom| atom.ast.constraints.aromatic_valence().cloned())
                .collect::<Vec<_>>(),
            expected
        );
    }

    #[rstest]
    #[case::retained(
        SmilesIoConfig::opensmiles(),
        ChemistryModel::default(),
        ResolveConfig::default(),
        Some(TetrahedralStereoAst::Stereo(StereoCosetAst::Lit(0)))
    )]
    #[case::reset(
        SmilesIoConfig::opensmiles(),
        ChemistryModel::default(),
        ResolveConfig {
            aromaticity: AromaticityResolveConfig::default(),
            stereo: StereoResolveConfig {
                reset_stereo_constraints: true,
            },
        },
        None
    )]
    fn test_ingest_smiles_with_stereo(
        #[case] io_config: SmilesIoConfig,
        #[case] model: ChemistryModel,
        #[case] resolve_config: ResolveConfig,
        #[case] expected: Option<TetrahedralStereoAst>,
    ) {
        let ast = ingest_smiles_with("C[C@H](N)O", &io_config, &model, &resolve_config).unwrap();

        assert_eq!(
            ast.atom(AtomId(1))
                .ast
                .constraints
                .tetrahedral_stereo()
                .cloned(),
            expected
        );
        assert!(ast.stereo_atoms().is_at(AtomId(1)));
    }

    #[rstest]
    #[case::contradiction(
        "[nH]1cccc1",
        SmilesIoConfig::opensmiles(),
        ChemistryModel {
            aromaticity: AromaticityModel::Clar {
                scope: ElementScope::Any,
                ring_limits: RingLimits::default(),
            },
            ..ChemistryModel::default()
        },
        ResolveConfig::default(),
        SmilesInputError::Contradiction(ResolverContradiction::Aromaticity(
            AromaticityContradiction::ClarNonBenzenoid(String::from(
                "Clar model requires benzenoid input but non-carbon aromatic atoms are present",
            )),
        ))
    )]
    #[case::underdetermined(
        "C",
        SmilesIoConfig::opensmiles(),
        ChemistryModel {
            valence: ValenceModel::AtomTyping {
                registry: Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!("C#c0")])),
            },
            ..ChemistryModel::default()
        },
        ResolveConfig::default(),
        SmilesInputError::Underdetermined(ResolveUnderdetermined)
    )]
    fn test_ingest_smiles_with_error(
        #[case] input: &str,
        #[case] io_config: SmilesIoConfig,
        #[case] model: ChemistryModel,
        #[case] resolve_config: ResolveConfig,
        #[case] expected: SmilesInputError,
    ) {
        assert_eq!(
            ingest_smiles_with(input, &io_config, &model, &resolve_config,),
            Err(expected)
        );
    }
}
