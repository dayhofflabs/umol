//! Conversion of parsed external-format values into graph models.

use std::any::Any;

use thiserror::Error;
use umol_ast::ast::{MoleculeAst, ReactionAst, TryIntoAst};
use umol_graph_core::{Correspondence, NodeId};
use umol_io::smiles::{ParseError as SmilesParseError, ReactionSmiles, Smiles, SmilesIoConfig};
use umol_io::table_ir::raise::RaiseError;
use umol_io::table_ir::Molecule as TableMolecule;
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

/// Failure while interpreting a parsed reaction representation.
#[derive(Clone, Debug, PartialEq, Eq, Error)]
pub enum ReactionInterpretationError {
    #[error("reactants: {0}")]
    Reactants(#[source] MoleculeInterpretationError),
    #[error("products: {0}")]
    Products(#[source] MoleculeInterpretationError),
    #[error(
        "atom-map class {class} cannot be projected into one correspondence \
         (reactant atoms: {reactant_count}, product atoms: {product_count})"
    )]
    AmbiguousAtomMapClass {
        class: u32,
        reactant_count: usize,
        product_count: usize,
    },
    #[error("reaction agents cannot be represented in ReactionAst")]
    AgentsUnsupported,
}

impl UmolError for ReactionInterpretationError {
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

/// Failure while accepting reaction SMILES text as a determined reaction.
#[derive(Clone, Debug, PartialEq, Error)]
pub enum ReactionSmilesInputError {
    #[error("{0}")]
    Syntax(#[from] SmilesParseError),
    #[error("{0}")]
    Interpretation(#[from] ReactionInterpretationError),
}

impl UmolError for ReactionSmilesInputError {
    fn as_any(&self) -> &dyn Any {
        self
    }
}

fn interpret_molecule(
    molecule: &TableMolecule,
    model: &ChemistryModel,
    resolve_config: &ResolveConfig,
) -> Result<MoleculeAst, MoleculeInterpretationError> {
    let mut ast: MoleculeAst = molecule.try_into_ast(&())?;
    match Resolver::with_config(model, *resolve_config).resolve(&mut ast)? {
        Solution::Determined(()) => Ok(ast),
        Solution::Underdetermined(()) => Err(ResolveUnderdetermined.into()),
        Solution::Contradictory(error) => Err(error.into()),
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
        interpret_molecule(self.as_table_ir(), model, resolve_config)
    }
}

impl Interpret for ReactionSmiles {
    type Output = ReactionAst;
    type Error = ReactionInterpretationError;

    fn interpret(
        &self,
        model: &ChemistryModel,
        resolve_config: &ResolveConfig,
    ) -> Result<Self::Output, Self::Error> {
        let reaction = self.as_table_ir();
        if reaction.agents.atom_count() != 0 {
            return Err(ReactionInterpretationError::AgentsUnsupported);
        }

        let mut matched_pairs = Vec::new();
        for (&class, (reactants, products)) in &reaction.atom_mapping {
            if reactants.len() > 1 || products.len() > 1 {
                return Err(ReactionInterpretationError::AmbiguousAtomMapClass {
                    class,
                    reactant_count: reactants.len(),
                    product_count: products.len(),
                });
            }
            if let ([reactant], [product]) = (reactants.as_slice(), products.as_slice()) {
                matched_pairs.push((NodeId(*reactant), NodeId(*product)));
            }
        }

        let lhs = interpret_molecule(&reaction.reactants, model, resolve_config)
            .map_err(ReactionInterpretationError::Reactants)?;
        let rhs = interpret_molecule(&reaction.products, model, resolve_config)
            .map_err(ReactionInterpretationError::Products)?;
        let atom_correspondence =
            Correspondence::new(matched_pairs, lhs.atoms().count(), rhs.atoms().count());

        Ok(ReactionAst::from_sides(lhs, rhs, atom_correspondence))
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

/// Ingest reaction SMILES text with the OpenSMILES configuration and default model.
pub fn ingest_reaction_smiles(input: &str) -> Result<ReactionAst, ReactionSmilesInputError> {
    ingest_reaction_smiles_bytes(input.as_bytes())
}

/// Ingest reaction SMILES bytes with the OpenSMILES configuration and default model.
pub fn ingest_reaction_smiles_bytes(input: &[u8]) -> Result<ReactionAst, ReactionSmilesInputError> {
    ingest_reaction_smiles_bytes_with(
        input,
        &SmilesIoConfig::opensmiles(),
        &ChemistryModel::default(),
        &ResolveConfig::default(),
    )
}

/// Ingest reaction SMILES text with explicit IO, chemistry, and resolve configuration.
pub fn ingest_reaction_smiles_with(
    input: &str,
    io_config: &SmilesIoConfig,
    model: &ChemistryModel,
    resolve_config: &ResolveConfig,
) -> Result<ReactionAst, ReactionSmilesInputError> {
    ingest_reaction_smiles_bytes_with(input.as_bytes(), io_config, model, resolve_config)
}

/// Ingest reaction SMILES bytes with explicit IO, chemistry, and resolve configuration.
pub fn ingest_reaction_smiles_bytes_with(
    input: &[u8],
    io_config: &SmilesIoConfig,
    model: &ChemistryModel,
    resolve_config: &ResolveConfig,
) -> Result<ReactionAst, ReactionSmilesInputError> {
    let reaction_smiles = ReactionSmiles::parse_bytes_with(input, io_config)?;
    reaction_smiles
        .interpret(model, resolve_config)
        .map_err(ReactionSmilesInputError::from)
}

#[cfg(test)]
mod tests {
    use std::borrow::Cow;
    use std::error::Error as _;

    use rstest::rstest;
    use umol_ast::ast::{
        AromaticSystemId, AromaticValenceAst, AtomId, StereoCoset, TetrahedralStereoAst, ValueAst,
    };
    use umol_ast::{atom_dsl, mol_dsl};

    use super::*;
    use crate::ops::aromaticity::{AromaticityContradiction, AromaticityError};
    use crate::ops::model::{AromaticityModel, ElementScope, RingLimits, ValenceModel};
    use crate::ops::resolve::{AromaticityResolveConfig, InconsistencyPolicy, StereoResolveConfig};
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
    #[case::reactants(
        ReactionInterpretationError::Reactants(
            MoleculeInterpretationError::ModelConversion(
                RaiseError::WedgeConflict { atom: 2 },
            ),
        ),
        "reactants: inconsistent wedge bonds at atom 2",
        Some("inconsistent wedge bonds at atom 2"),
    )]
    #[case::products(
        ReactionInterpretationError::Products(MoleculeInterpretationError::Underdetermined(
            ResolveUnderdetermined
        ),),
        "products: resolution underdetermined",
        Some("resolution underdetermined")
    )]
    #[case::ambiguous_atom_map_class(
        ReactionInterpretationError::AmbiguousAtomMapClass {
            class: 7,
            reactant_count: 2,
            product_count: 1,
        },
        "atom-map class 7 cannot be projected into one correspondence \
         (reactant atoms: 2, product atoms: 1)",
        None,
    )]
    #[case::agents_unsupported(
        ReactionInterpretationError::AgentsUnsupported,
        "reaction agents cannot be represented in ReactionAst",
        None
    )]
    fn test_reaction_interpretation_error(
        #[case] error: ReactionInterpretationError,
        #[case] expected: &str,
        #[case] expected_source: Option<&str>,
    ) {
        assert_eq!(error.to_string(), expected);
        assert_eq!(
            error.source().map(ToString::to_string).as_deref(),
            expected_source
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
    #[case::syntax(
        ReactionSmilesInputError::Syntax(SmilesParseError::LeadingWhitespace),
        "Leading whitespace",
        vec!["Leading whitespace"],
    )]
    #[case::reactants(
        ReactionSmilesInputError::Interpretation(ReactionInterpretationError::Reactants(
            MoleculeInterpretationError::ModelConversion(
                RaiseError::WedgeConflict { atom: 2 },
            ),
        )),
        "reactants: inconsistent wedge bonds at atom 2",
        vec![
            "reactants: inconsistent wedge bonds at atom 2",
            "inconsistent wedge bonds at atom 2",
            "inconsistent wedge bonds at atom 2",
        ],
    )]
    #[case::products(
        ReactionSmilesInputError::Interpretation(ReactionInterpretationError::Products(
            MoleculeInterpretationError::Underdetermined(ResolveUnderdetermined),
        )),
        "products: resolution underdetermined",
        vec![
            "products: resolution underdetermined",
            "resolution underdetermined",
            "resolution underdetermined",
        ],
    )]
    #[case::ambiguous_atom_map_class(
        ReactionSmilesInputError::Interpretation(
            ReactionInterpretationError::AmbiguousAtomMapClass {
                class: 7,
                reactant_count: 2,
                product_count: 1,
            },
        ),
        "atom-map class 7 cannot be projected into one correspondence \
         (reactant atoms: 2, product atoms: 1)",
        vec![
            "atom-map class 7 cannot be projected into one correspondence \
             (reactant atoms: 2, product atoms: 1)",
        ],
    )]
    #[case::agents_unsupported(
        ReactionSmilesInputError::Interpretation(
            ReactionInterpretationError::AgentsUnsupported,
        ),
        "reaction agents cannot be represented in ReactionAst",
        vec!["reaction agents cannot be represented in ReactionAst"],
    )]
    fn test_reaction_smiles_input_error(
        #[case] error: ReactionSmilesInputError,
        #[case] expected: &str,
        #[case] expected_sources: Vec<&str>,
    ) {
        assert_eq!(error.to_string(), expected);

        let mut source = error.source();
        let mut actual_sources = Vec::new();
        while let Some(current) = source {
            actual_sources.push(current.to_string());
            source = current.source();
        }

        assert_eq!(
            actual_sources,
            expected_sources
                .into_iter()
                .map(String::from)
                .collect::<Vec<_>>()
        );
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
    #[case::mapped(
        "[CH4:1]>>[CH4:1]",
        r#"{:deltas [] :lhs {:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"] :bonds []}}"#.parse().unwrap(),
    )]
    #[case::atom_change(
        "[NH4+:1]>>[NH3:1]",
        r##"{:deltas [{:atom {:modify [0 "#c0#h3#n"]}}] :lhs {:atoms ["N#i=#c+#h4#n0#u0#s#v0#d0#t0#a!#m!"] :bonds []}}"##.parse().unwrap(),
    )]
    #[case::bond_change(
        "[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]",
        r##"{:deltas [{:atom {:modify [0 "#h3#v"]}} {:atom {:modify [1 "#h3#v"]}} {:bond {:modify [0 "1"]}}] :lhs {:atoms ["C#i=#c0#h2#n0#u0#s#v2#d0#t0#a!#m!" "C#i=#c0#h2#n0#u0#s#v2#d0#t0#a!#m!"] :bonds [[0 1 "2#c0#u0#s"]]}}"##.parse().unwrap(),
    )]
    #[case::reactant_only(
        "[CH4:1].[OH2:2]>>[CH4:1]",
        r#"{:deltas [{:atom {:remove 1}}] :lhs {:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!" "O#i=#c0#h2#n2#u0#s#v0#d0#t0#a!#m!"] :bonds []}}"#.parse().unwrap(),
    )]
    #[case::product_only(
        "[CH4:1]>>[CH4:1].[OH2:2]",
        r#"{:deltas [{:atom {:add "O#i=#c0#h2#n2#u0#s#v0#d0#t0#a!#m!"}}] :lhs {:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"] :bonds []}}"#.parse().unwrap(),
    )]
    #[case::unmapped(
        "C>>O",
        r#"{:deltas [{:atom {:remove 0}} {:atom {:add "O#i=#c0#h2#n2#u0#s#v0#d0#t0#a!#m!"}}] :lhs {:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!"] :bonds []}}"#.parse().unwrap(),
    )]
    #[case::reordered(
        "[CH4:1].[OH2:2]>>[OH2:2].[CH4:1]",
        r#"{:deltas [] :lhs {:atoms ["C#i=#c0#h4#n0#u0#s#v0#d0#t0#a!#m!" "O#i=#c0#h2#n2#u0#s#v0#d0#t0#a!#m!"] :bonds []}}"#.parse().unwrap(),
    )]
    fn test_reaction_smiles_interpret(#[case] input: &str, #[case] expected: ReactionAst) {
        let reaction = ReactionSmiles::parse(input).unwrap();

        assert_eq!(
            reaction.interpret(&ChemistryModel::default(), &ResolveConfig::default()),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::agents(
        "*>O>*",
        ChemistryModel::default(),
        ReactionInterpretationError::AgentsUnsupported
    )]
    #[case::ambiguous_reactants(
        "[*:1].[*:1]>>[*:1]",
        ChemistryModel::default(),
        ReactionInterpretationError::AmbiguousAtomMapClass {
            class: 1,
            reactant_count: 2,
            product_count: 1,
        },
    )]
    #[case::ambiguous_products(
        "[*:1]>>[*:1].[*:1]",
        ChemistryModel::default(),
        ReactionInterpretationError::AmbiguousAtomMapClass {
            class: 1,
            reactant_count: 1,
            product_count: 2,
        },
    )]
    #[case::ambiguous_both(
        "[*:1].[*:1]>>[*:1].[*:1]",
        ChemistryModel::default(),
        ReactionInterpretationError::AmbiguousAtomMapClass {
            class: 1,
            reactant_count: 2,
            product_count: 2,
        },
    )]
    #[case::reactants_model_conversion(
        "C[S@]C>>",
        ChemistryModel::default(),
        ReactionInterpretationError::Reactants(
            MoleculeInterpretationError::ModelConversion(
                RaiseError::TetrahedralLigandCount { atom: 1, count: 2 },
            ),
        ),
    )]
    #[case::products_model_conversion(
        ">>C[S@]C",
        ChemistryModel::default(),
        ReactionInterpretationError::Products(
            MoleculeInterpretationError::ModelConversion(
                RaiseError::TetrahedralLigandCount { atom: 1, count: 2 },
            ),
        ),
    )]
    #[case::reactants_underdetermined(
        "*>>",
        ChemistryModel::default(),
        ReactionInterpretationError::Reactants(MoleculeInterpretationError::Underdetermined(
            ResolveUnderdetermined
        ),)
    )]
    #[case::products_underdetermined(
        ">>*",
        ChemistryModel::default(),
        ReactionInterpretationError::Products(MoleculeInterpretationError::Underdetermined(
            ResolveUnderdetermined
        ),)
    )]
    #[case::reactants_contradiction(
        "[nH]1cccc1>>",
        ChemistryModel {
            aromaticity: AromaticityModel::Clar {
                scope: ElementScope::Any,
                ring_limits: RingLimits::default(),
            },
            ..ChemistryModel::default()
        },
        ReactionInterpretationError::Reactants(
            MoleculeInterpretationError::Contradiction(
                ResolverContradiction::Aromaticity(
                    AromaticityContradiction::ClarNonBenzenoid(String::from(
                        "Clar model requires benzenoid input but non-carbon aromatic atoms are present",
                    )),
                ),
            ),
        ),
    )]
    #[case::products_contradiction(
        ">>[nH]1cccc1",
        ChemistryModel {
            aromaticity: AromaticityModel::Clar {
                scope: ElementScope::Any,
                ring_limits: RingLimits::default(),
            },
            ..ChemistryModel::default()
        },
        ReactionInterpretationError::Products(
            MoleculeInterpretationError::Contradiction(
                ResolverContradiction::Aromaticity(
                    AromaticityContradiction::ClarNonBenzenoid(String::from(
                        "Clar model requires benzenoid input but non-carbon aromatic atoms are present",
                    )),
                ),
            ),
        ),
    )]
    #[case::reactants_execution(
        "c1ccccc1>>",
        ChemistryModel {
            valence: ValenceModel::AtomTyping {
                registry: Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!(
                    "C#i=#c0#h0#n0#u0#s#v2#a2"
                )])),
            },
            aromaticity: AromaticityModel::Hmo {
                scope: ElementScope::Any,
                stabilization_threshold: 0.5,
            },
            ..ChemistryModel::default()
        },
        ReactionInterpretationError::Reactants(
            MoleculeInterpretationError::Execution(
                ResolverError::Aromaticity(AromaticityError::HmoMissingParameters(
                    String::from("no Van-Catledge parameters for C with 2 pi-electrons"),
                )),
            ),
        ),
    )]
    #[case::products_execution(
        ">>c1ccccc1",
        ChemistryModel {
            valence: ValenceModel::AtomTyping {
                registry: Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!(
                    "C#i=#c0#h0#n0#u0#s#v2#a2"
                )])),
            },
            aromaticity: AromaticityModel::Hmo {
                scope: ElementScope::Any,
                stabilization_threshold: 0.5,
            },
            ..ChemistryModel::default()
        },
        ReactionInterpretationError::Products(
            MoleculeInterpretationError::Execution(
                ResolverError::Aromaticity(AromaticityError::HmoMissingParameters(
                    String::from("no Van-Catledge parameters for C with 2 pi-electrons"),
                )),
            ),
        ),
    )]
    fn test_reaction_smiles_interpret_error(
        #[case] input: &str,
        #[case] model: ChemistryModel,
        #[case] expected: ReactionInterpretationError,
    ) {
        let reaction = ReactionSmiles::parse(input).unwrap();

        assert_eq!(
            reaction.interpret(&model, &ResolveConfig::default()),
            Err(expected)
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
        Some(TetrahedralStereoAst::Stereo(StereoCoset::Lit(0)))
    )]
    #[case::reset(
        SmilesIoConfig::opensmiles(),
        ChemistryModel::default(),
        ResolveConfig {
            aromaticity: AromaticityResolveConfig::default(),
            stereo: StereoResolveConfig {
                reset_stereo_constraints: true,
                inconsistency: InconsistencyPolicy::Error,
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

    #[rstest]
    #[case::mapped("[CH4:1]>>[CH4:1]")]
    #[case::unmapped("C>>O")]
    fn test_ingest_reaction_smiles(#[case] input: &str) {
        assert_eq!(
            ingest_reaction_smiles(input),
            ingest_reaction_smiles_with(
                input,
                &SmilesIoConfig::opensmiles(),
                &ChemistryModel::default(),
                &ResolveConfig::default(),
            )
        );
    }

    #[rstest]
    #[case::syntax(
        " C>>C",
        ReactionSmilesInputError::Syntax(SmilesParseError::LeadingWhitespace)
    )]
    #[case::extended_bond(
        "C~C>>C.C",
        ReactionSmilesInputError::Syntax(SmilesParseError::InvalidToken { pos: 1 }),
    )]
    #[case::agents(
        "C>O>C",
        ReactionSmilesInputError::Interpretation(ReactionInterpretationError::AgentsUnsupported,)
    )]
    #[case::ambiguous_atom_map_class(
        "[C:1].[O:1]>>[C:1]",
        ReactionSmilesInputError::Interpretation(
            ReactionInterpretationError::AmbiguousAtomMapClass {
                class: 1,
                reactant_count: 2,
                product_count: 1,
            },
        ),
    )]
    #[case::underdetermined(
        "*>>C",
        ReactionSmilesInputError::Interpretation(ReactionInterpretationError::Reactants(
            MoleculeInterpretationError::Underdetermined(ResolveUnderdetermined),
        ),)
    )]
    fn test_ingest_reaction_smiles_error(
        #[case] input: &str,
        #[case] expected: ReactionSmilesInputError,
    ) {
        assert_eq!(ingest_reaction_smiles(input), Err(expected));
    }

    #[rstest]
    #[case::mapped(b"[CH4:1]>>[CH4:1]")]
    #[case::unmapped(b"C>>O")]
    fn test_ingest_reaction_smiles_bytes(#[case] input: &[u8]) {
        assert_eq!(
            ingest_reaction_smiles_bytes(input),
            ingest_reaction_smiles_bytes_with(
                input,
                &SmilesIoConfig::opensmiles(),
                &ChemistryModel::default(),
                &ResolveConfig::default(),
            )
        );
    }

    #[rstest]
    #[case::syntax(
        b" C>>C",
        ReactionSmilesInputError::Syntax(SmilesParseError::LeadingWhitespace)
    )]
    #[case::interpretation(
        b"C>O>C",
        ReactionSmilesInputError::Interpretation(ReactionInterpretationError::AgentsUnsupported,)
    )]
    fn test_ingest_reaction_smiles_bytes_error(
        #[case] input: &[u8],
        #[case] expected: ReactionSmilesInputError,
    ) {
        assert_eq!(ingest_reaction_smiles_bytes(input), Err(expected));
    }

    #[rstest]
    #[case::io(
        "C~C>>C.C",
        SmilesIoConfig::lenient(),
        ChemistryModel::default(),
        ResolveConfig::default(),
        Err(
            ReactionSmilesInputError::Interpretation(ReactionInterpretationError::Reactants(
                MoleculeInterpretationError::Underdetermined(ResolveUnderdetermined),
            ),)
        )
    )]
    #[case::chemistry(
        "[nH]1cccc1>>",
        SmilesIoConfig::opensmiles(),
        ChemistryModel {
            aromaticity: AromaticityModel::Clar {
                scope: ElementScope::Any,
                ring_limits: RingLimits::default(),
            },
            ..ChemistryModel::default()
        },
        ResolveConfig::default(),
        Err(ReactionSmilesInputError::Interpretation(
            ReactionInterpretationError::Reactants(
                MoleculeInterpretationError::Contradiction(
                    ResolverContradiction::Aromaticity(
                        AromaticityContradiction::ClarNonBenzenoid(String::from(
                            "Clar model requires benzenoid input but non-carbon aromatic atoms are present",
                        )),
                    ),
                ),
            ),
        )),
    )]
    #[case::resolve(
        "[cH+:1]1[cH:2][cH:3]1>>[cH+:1]1[cH:2][cH:3]1",
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
        Ok(r##"{:deltas [] :lhs {:aromatic-systems [{:atoms [0 1 2] :type "[0,1,1]#c0#u0#s"}] :atoms ["C#i=#c+#h#n0#u0#s#v2#d0#t0#a0#m!" "C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!" "C#i=#c0#h#n0#u0#s#v2#d0#t0#a#m!"] :bonds [[0 2 "1#c0#u0#s#a"] [0 1 "1#c0#u0#s#a"] [1 2 "1#c0#u0#s#a"]]}}"##.parse().unwrap()),
    )]
    fn test_ingest_reaction_smiles_with(
        #[case] input: &str,
        #[case] io_config: SmilesIoConfig,
        #[case] model: ChemistryModel,
        #[case] resolve_config: ResolveConfig,
        #[case] expected: Result<ReactionAst, ReactionSmilesInputError>,
    ) {
        assert_eq!(
            ingest_reaction_smiles_with(input, &io_config, &model, &resolve_config),
            expected
        );
    }

    #[rstest]
    #[case::lenient(
        "C~C>>C.C",
        SmilesIoConfig::lenient(),
        ChemistryModel::default(),
        ResolveConfig::default()
    )]
    fn test_ingest_reaction_smiles_bytes_with(
        #[case] input: &str,
        #[case] io_config: SmilesIoConfig,
        #[case] model: ChemistryModel,
        #[case] resolve_config: ResolveConfig,
    ) {
        assert_eq!(
            ingest_reaction_smiles_bytes_with(
                input.as_bytes(),
                &io_config,
                &model,
                &resolve_config,
            ),
            ingest_reaction_smiles_with(input, &io_config, &model, &resolve_config)
        );
    }
}
