//! Conversion of parsed external-format values into graph models.

use std::any::Any;

use thiserror::Error;
use umol_graph_core::Correspondence;
use umol_graph_ir::ir::{AtomId, Molecule, Reaction, TryIntoIr};
use umol_io::smiles::{ParseError as SmilesParseError, ReactionSmiles, Smiles, SmilesIoConfig};
use umol_io::table_ir::raise::RaiseError;
use umol_io::table_ir::Molecule as TableMolecule;
use umol_utils::error::UmolError;
use umol_utils::solution::Solution;

use crate::ops::model::{ChemistryModel, ValenceModel};
use crate::ops::resolve::{
    ResolveConfig, ResolveContradiction, ResolveError, ResolveUnderdetermined, Resolver,
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
    Contradiction(#[from] ResolveContradiction),
    #[error("{0}")]
    Underdetermined(#[from] ResolveUnderdetermined),
    #[error("{0}")]
    Execution(#[from] ResolveError),
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
    #[error("reaction agents cannot be represented in Reaction")]
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
    Contradiction(#[from] ResolveContradiction),
    #[error("{0}")]
    Underdetermined(#[from] ResolveUnderdetermined),
    #[error("{0}")]
    Execution(#[from] ResolveError),
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
) -> Result<Molecule, MoleculeInterpretationError> {
    let mut molecule: Molecule = molecule.try_into_ir(&())?;
    match Resolver::with_config(model, *resolve_config).resolve(&mut molecule)? {
        Solution::Determined(_) => Ok(molecule),
        Solution::Underdetermined(report) => Err(ResolveUnderdetermined { report }.into()),
        Solution::Contradictory(error) => Err(error.into()),
    }
}

impl Interpret for Smiles {
    type Output = Molecule;
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
    type Output = Reaction;
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
                matched_pairs.push((AtomId(*reactant), AtomId(*product)));
            }
        }

        let lhs = interpret_molecule(&reaction.reactants, model, resolve_config)
            .map_err(ReactionInterpretationError::Reactants)?;
        let rhs = interpret_molecule(&reaction.products, model, resolve_config)
            .map_err(ReactionInterpretationError::Products)?;
        let atom_correspondence =
            Correspondence::new(matched_pairs, lhs.atoms().count(), rhs.atoms().count())
                .expect("correspondence producer preserves partial-bijection invariants");

        Ok(Reaction::from_sides(lhs, rhs, atom_correspondence)
            .expect("interpreted reaction sides preserve unique entity incidence"))
    }
}

/// Ingest SMILES text with the OpenSMILES configuration and the SMILES
/// valence preset — the reader carries its format's convention.
pub fn ingest_smiles(input: &str) -> Result<Molecule, SmilesInputError> {
    ingest_smiles_bytes(input.as_bytes())
}

/// Ingest SMILES bytes with the OpenSMILES configuration and the SMILES
/// valence preset — the reader carries its format's convention.
pub fn ingest_smiles_bytes(input: &[u8]) -> Result<Molecule, SmilesInputError> {
    ingest_smiles_bytes_with(
        input,
        &SmilesIoConfig::opensmiles(),
        &ChemistryModel {
            valence: ValenceModel::smiles(),
            ..ChemistryModel::default()
        },
        &ResolveConfig::default(),
    )
}

/// Ingest SMILES text with explicit IO, chemistry, and resolve configuration.
pub fn ingest_smiles_with(
    input: &str,
    io_config: &SmilesIoConfig,
    model: &ChemistryModel,
    resolve_config: &ResolveConfig,
) -> Result<Molecule, SmilesInputError> {
    ingest_smiles_bytes_with(input.as_bytes(), io_config, model, resolve_config)
}

/// Ingest SMILES bytes with explicit IO, chemistry, and resolve configuration.
pub fn ingest_smiles_bytes_with(
    input: &[u8],
    io_config: &SmilesIoConfig,
    model: &ChemistryModel,
    resolve_config: &ResolveConfig,
) -> Result<Molecule, SmilesInputError> {
    let smiles = Smiles::parse_bytes_with(input, io_config)?;
    smiles
        .interpret(model, resolve_config)
        .map_err(SmilesInputError::from)
}

/// Ingest reaction SMILES text with the OpenSMILES configuration and the
/// SMILES valence preset — the reader carries its format's convention.
pub fn ingest_reaction_smiles(input: &str) -> Result<Reaction, ReactionSmilesInputError> {
    ingest_reaction_smiles_bytes(input.as_bytes())
}

/// Ingest reaction SMILES bytes with the OpenSMILES configuration and the
/// SMILES valence preset — the reader carries its format's convention.
pub fn ingest_reaction_smiles_bytes(input: &[u8]) -> Result<Reaction, ReactionSmilesInputError> {
    ingest_reaction_smiles_bytes_with(
        input,
        &SmilesIoConfig::opensmiles(),
        &ChemistryModel {
            valence: ValenceModel::smiles(),
            ..ChemistryModel::default()
        },
        &ResolveConfig::default(),
    )
}

/// Ingest reaction SMILES text with explicit IO, chemistry, and resolve configuration.
pub fn ingest_reaction_smiles_with(
    input: &str,
    io_config: &SmilesIoConfig,
    model: &ChemistryModel,
    resolve_config: &ResolveConfig,
) -> Result<Reaction, ReactionSmilesInputError> {
    ingest_reaction_smiles_bytes_with(input.as_bytes(), io_config, model, resolve_config)
}

/// Ingest reaction SMILES bytes with explicit IO, chemistry, and resolve configuration.
pub fn ingest_reaction_smiles_bytes_with(
    input: &[u8],
    io_config: &SmilesIoConfig,
    model: &ChemistryModel,
    resolve_config: &ResolveConfig,
) -> Result<Reaction, ReactionSmilesInputError> {
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
    use smallvec::smallvec;
    use umol_graph_ir::ir::{
        AromaticSystemId, AromaticValenceForm, AtomId, BooleanForm, Deltas, ElectronCountsForm,
        NumForm, TetrahedralStereoForm,
    };
    use umol_graph_ir::{atom_dsl, mol_dsl};

    use super::*;
    use crate::ops::aromaticity::{
        AromaticityContradiction, AromaticityError, AromaticityInconsistency,
    };
    use crate::ops::model::{
        AromaticityModel, AromaticityRule, AromaticityTieBreak, ElementScope, RingLimits,
        ValenceModel, ValenceTieBreak,
    };
    use crate::ops::resolve::{
        AromaticityFailurePolicy, AromaticityResolveConfig, StereoResolveConfig,
    };
    use crate::ops::valence::{AtomCompletions, AtomTypeRegistry, ResolveReport};

    #[rstest]
    #[case::model_conversion(
        MoleculeInterpretationError::ModelConversion(RaiseError::WedgeConflict { atom: 2 }),
        "inconsistent wedge bonds at atom 2"
    )]
    #[case::contradiction(
        MoleculeInterpretationError::Contradiction(ResolveContradiction::Aromaticity(
            AromaticityContradiction::HmoInvalidInput(String::from("invalid input")),
        )),
        "hmo: invalid input: invalid input"
    )]
    #[case::underdetermined(
        MoleculeInterpretationError::Underdetermined(ResolveUnderdetermined::default()),
        "resolution underdetermined"
    )]
    #[case::execution(
        MoleculeInterpretationError::Execution(ResolveError::Aromaticity(
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
            ResolveUnderdetermined::default()
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
        "reaction agents cannot be represented in Reaction",
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
        SmilesInputError::Contradiction(ResolveContradiction::Aromaticity(
            AromaticityContradiction::HmoInvalidInput(String::from("invalid input")),
        )),
        "hmo: invalid input: invalid input"
    )]
    #[case::underdetermined(
        SmilesInputError::Underdetermined(ResolveUnderdetermined::default()),
        "resolution underdetermined"
    )]
    #[case::execution(
        SmilesInputError::Execution(ResolveError::Aromaticity(
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
        MoleculeInterpretationError::Contradiction(ResolveContradiction::Aromaticity(
            AromaticityContradiction::HmoInvalidInput(String::from("invalid input")),
        )),
        SmilesInputError::Contradiction(ResolveContradiction::Aromaticity(
            AromaticityContradiction::HmoInvalidInput(String::from("invalid input")),
        ))
    )]
    #[case::underdetermined(
        MoleculeInterpretationError::Underdetermined(ResolveUnderdetermined::default()),
        SmilesInputError::Underdetermined(ResolveUnderdetermined::default())
    )]
    #[case::execution(
        MoleculeInterpretationError::Execution(ResolveError::Aromaticity(
            AromaticityError::HmoMissingParameters(String::from("carbon")),
        )),
        SmilesInputError::Execution(ResolveError::Aromaticity(
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
            MoleculeInterpretationError::Underdetermined(ResolveUnderdetermined::default()),
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
        "reaction agents cannot be represented in Reaction",
        vec!["reaction agents cannot be represented in Reaction"],
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
    #[case::default_aromaticity(
        ChemistryModel {
            valence: ValenceModel {
                tie_break: ValenceTieBreak::MostSaturated,
                ..ValenceModel::default()
            },
            ..ChemistryModel::default()
        },
        ResolveConfig::default()
    )]
    #[case::permissive(
        ChemistryModel {
            valence: ValenceModel {
                tie_break: ValenceTieBreak::MostSaturated,
                ..ValenceModel::default()
            },
            aromaticity: AromaticityModel::permissive(),
            ..ChemistryModel::default()
        },
        ResolveConfig::default()
    )]
    fn test_smiles_interpret(#[case] model: ChemistryModel, #[case] resolve_config: ResolveConfig) {
        let smiles = Smiles::parse("C").unwrap();
        let expected = mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s"]}"#);

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
                ResolveUnderdetermined::default()
            ))
        );
    }

    #[rstest]
    #[case::mapped(
        "[CH4:1]>>[CH4:1]",
        r#"{:deltas [] :lhs {:atoms ["C#i=#c0#h4#n0#u0#s"] :bonds []}}"#.parse().unwrap(),
    )]
    #[case::atom_change(
        "[NH4+:1]>>[NH3:1]",
        r##"{:deltas [{:atom {:modify [0 "#c0#h3#n"]}}] :lhs {:atoms ["N#i=#c+#h4#n0#u0#s"] :bonds []}}"##.parse().unwrap(),
    )]
    #[case::bond_change(
        "[CH2:1]=[CH2:2]>>[CH3:1][CH3:2]",
        r##"{:deltas [{:atom {:modify [0 "#h3"]}} {:atom {:modify [1 "#h3"]}} {:bond {:modify [0 "1"]}}] :lhs {:atoms ["C#i=#c0#h2#n0#u0#s" "C#i=#c0#h2#n0#u0#s"] :bonds [[0 1 "2#c0#u0#s"]]}}"##.parse().unwrap(),
    )]
    #[case::reactant_only(
        "[CH4:1].[OH2:2]>>[CH4:1]",
        r#"{:deltas [{:atom {:remove 1}}] :lhs {:atoms ["C#i=#c0#h4#n0#u0#s" "O#i=#c0#h2#n2#u0#s"] :bonds []}}"#.parse().unwrap(),
    )]
    #[case::product_only(
        "[CH4:1]>>[CH4:1].[OH2:2]",
        r#"{:deltas [{:atom {:add "O#i=#c0#h2#n2#u0#s"}}] :lhs {:atoms ["C#i=#c0#h4#n0#u0#s"] :bonds []}}"#.parse().unwrap(),
    )]
    #[case::unmapped(
        "C>>O",
        r#"{:deltas [{:atom {:remove 0}} {:atom {:add "O#i=#c0#h2#n2#u0#s"}}] :lhs {:atoms ["C#i=#c0#h4#n0#u0#s"] :bonds []}}"#.parse().unwrap(),
    )]
    #[case::reordered(
        "[CH4:1].[OH2:2]>>[OH2:2].[CH4:1]",
        r#"{:deltas [] :lhs {:atoms ["C#i=#c0#h4#n0#u0#s" "O#i=#c0#h2#n2#u0#s"] :bonds []}}"#.parse().unwrap(),
    )]
    fn test_reaction_smiles_interpret(#[case] input: &str, #[case] expected: Reaction) {
        let reaction = ReactionSmiles::parse(input).unwrap();

        let model = ChemistryModel {
            valence: ValenceModel {
                tie_break: ValenceTieBreak::MostSaturated,
                ..ValenceModel::default()
            },
            ..ChemistryModel::default()
        };

        assert_eq!(
            reaction.interpret(&model, &ResolveConfig::default()),
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
            ResolveUnderdetermined::default()
        ),)
    )]
    #[case::products_underdetermined(
        ">>*",
        ChemistryModel::default(),
        ReactionInterpretationError::Products(MoleculeInterpretationError::Underdetermined(
            ResolveUnderdetermined::default()
        ),)
    )]
    #[case::reactants_contradiction(
        "[nH]1cccc1>>",
        ChemistryModel {
            aromaticity: AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Clar { ring_limits: RingLimits::default() }, tie_break: AromaticityTieBreak::Strict },
            ..ChemistryModel::default()
        },
        ReactionInterpretationError::Reactants(
            MoleculeInterpretationError::Contradiction(
                ResolveContradiction::Aromaticity(
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
            aromaticity: AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Clar { ring_limits: RingLimits::default() }, tie_break: AromaticityTieBreak::Strict },
            ..ChemistryModel::default()
        },
        ReactionInterpretationError::Products(
            MoleculeInterpretationError::Contradiction(
                ResolveContradiction::Aromaticity(
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
            valence: ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!(
                    "C#i=#c0#h0#n0#u0#s#v2#a2"
                )]))),
            aromaticity: AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Hmo { stabilization_threshold: 0.5 }, tie_break: AromaticityTieBreak::Strict },
            ..ChemistryModel::default()
        },
        ReactionInterpretationError::Reactants(
            MoleculeInterpretationError::Execution(
                ResolveError::Aromaticity(AromaticityError::HmoMissingParameters(
                    String::from("no Van-Catledge parameters for C with 2 pi-electrons"),
                )),
            ),
        ),
    )]
    #[case::products_execution(
        ">>c1ccccc1",
        ChemistryModel {
            valence: ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!(
                    "C#i=#c0#h0#n0#u0#s#v2#a2"
                )]))),
            aromaticity: AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Hmo { stabilization_threshold: 0.5 }, tie_break: AromaticityTieBreak::Strict },
            ..ChemistryModel::default()
        },
        ReactionInterpretationError::Products(
            MoleculeInterpretationError::Execution(
                ResolveError::Aromaticity(AromaticityError::HmoMissingParameters(
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
                &ChemistryModel {
                    valence: ValenceModel::smiles(),
                    ..ChemistryModel::default()
                },
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
    #[case::underdetermined(
        "*",
        SmilesInputError::Underdetermined(ResolveUnderdetermined::default())
    )]
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
                &ChemistryModel {
                    valence: ValenceModel::smiles(),
                    ..ChemistryModel::default()
                },
                &ResolveConfig::default(),
            )
        );
    }

    #[rstest]
    #[case::localized(
        SmilesIoConfig::opensmiles(),
        ChemistryModel::default(),
        ResolveConfig::default(),
        vec![NumForm::Lit(1), NumForm::Lit(0), NumForm::Lit(0)],
        NumForm::Lit(0)
    )]
    fn test_ingest_smiles_with_charge(
        #[case] io_config: SmilesIoConfig,
        #[case] model: ChemistryModel,
        #[case] resolve_config: ResolveConfig,
        #[case] expected_atom_charges: Vec<NumForm>,
        #[case] expected_system_charge: NumForm,
    ) {
        let molecule =
            ingest_smiles_with("[cH+]1[cH][cH]1", &io_config, &model, &resolve_config).unwrap();

        assert_eq!(
            molecule
                .atoms()
                .iter()
                .map(|atom| atom.attributes.charge.clone())
                .collect::<Vec<_>>(),
            expected_atom_charges
        );
        assert_eq!(
            molecule
                .aromatic_system(AromaticSystemId(0))
                .attributes
                .charge,
            expected_system_charge
        );
    }

    #[rstest]
    #[case::retained(
        SmilesIoConfig::opensmiles(),
        ChemistryModel::default(),
        ResolveConfig::default(),
        vec![None; 6]
    )]
    #[case::reset(
        SmilesIoConfig::opensmiles(),
        ChemistryModel::default(),
        ResolveConfig {
            aromaticity: AromaticityResolveConfig {
                reset_aromatic_valence: true,
                ..AromaticityResolveConfig::default()
            },
            stereo: StereoResolveConfig::default(),
        },
        vec![None; 6]
    )]
    fn test_ingest_smiles_with_aromatic_valence(
        #[case] io_config: SmilesIoConfig,
        #[case] model: ChemistryModel,
        #[case] resolve_config: ResolveConfig,
        #[case] expected: Vec<Option<AromaticValenceForm>>,
    ) {
        let molecule = ingest_smiles_with("c1ccccc1", &io_config, &model, &resolve_config).unwrap();

        assert_eq!(
            molecule
                .atoms()
                .iter()
                .map(|atom| atom.attributes.constraints.aromatic_valence().cloned())
                .collect::<Vec<_>>(),
            expected
        );
    }

    #[rstest]
    #[case::mdl_benzene(
        "c1ccccc1",
        AromaticityModel::mdl(),
        vec![1, 1, 1, 1, 1, 1],
    )]
    #[case::mdl_pyridine(
        "n1ccccc1",
        AromaticityModel::mdl(),
        vec![1, 1, 1, 1, 1, 1],
    )]
    #[case::daylight_furan(
        "o1cccc1",
        AromaticityModel::daylight(),
        vec![2, 1, 1, 1, 1],
    )]
    #[case::daylight_thiophene(
        "s1cccc1",
        AromaticityModel::daylight(),
        vec![2, 1, 1, 1, 1],
    )]
    #[case::daylight_pyrrole(
        "[nH]1cccc1",
        AromaticityModel::daylight(),
        vec![2, 1, 1, 1, 1],
    )]
    fn test_ingest_smiles_with_aromaticity(
        #[case] input: &str,
        #[case] aromaticity: AromaticityModel,
        #[case] expected_electrons: Vec<i64>,
    ) {
        let model = ChemistryModel {
            aromaticity,
            ..ChemistryModel::default()
        };
        let molecule = ingest_smiles_with(
            input,
            &SmilesIoConfig::opensmiles(),
            &model,
            &ResolveConfig::default(),
        )
        .unwrap();
        let system = molecule.aromatic_system(AromaticSystemId(0));

        assert_eq!(
            system.atom_ids().collect::<Vec<_>>(),
            (0..expected_electrons.len())
                .map(|index| AtomId(u32::try_from(index).unwrap()))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            system.electrons(),
            &ElectronCountsForm::Lit(expected_electrons)
        );
        assert_eq!(
            molecule.aromatic_systems().ids().collect::<Vec<_>>(),
            vec![AromaticSystemId(0)]
        );
    }

    #[rstest]
    #[case::pyrrole(
        "c1cccn1",
        mol_dsl!(r##"{:aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,1,1,2]#c0#u0#s"}] :atoms ["C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "N#i=#c0#h#n0#u0#s"] :bonds [[0 4 "1#c0#u0#s"] [0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"] [3 4 "1#c0#u0#s"]]}"##)
    )]
    #[case::pyridine(
        "c1ccncc1",
        mol_dsl!(r##"{:aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]#c0#u0#s"}] :atoms ["C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "N#i=#c0#h0#n#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s"] :bonds [[0 5 "1#c0#u0#s"] [0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"] [3 4 "1#c0#u0#s"] [4 5 "1#c0#u0#s"]]}"##)
    )]
    #[case::imidazole(
        "c1cncn1",
        mol_dsl!(r##"{:aromatic-systems [{:atoms [0 1 2 3 4] :attrs "[1,1,2,1,1]#c0#u0#s"}] :atoms ["C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "N#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "N#i=#c0#h0#n#u0#s"] :bonds [[0 4 "1#c0#u0#s"] [0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"] [3 4 "1#c0#u0#s"]]}"##)
    )]
    #[case::quinoline(
        "c1ccc2ccccc2n1",
        mol_dsl!(r##"{:aromatic-systems [{:atoms [0 1 2 3 4 5 6 7 8 9] :attrs "[1,1,1,1,1,1,1,1,1,1]#c0#u0#s"}] :atoms ["C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h0#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h0#n0#u0#s" "N#i=#c0#h0#n#u0#s"] :bonds [[0 9 "1#c0#u0#s"] [0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"] [3 8 "1#c0#u0#s"] [3 4 "1#c0#u0#s"] [4 5 "1#c0#u0#s"] [5 6 "1#c0#u0#s"] [6 7 "1#c0#u0#s"] [7 8 "1#c0#u0#s"] [8 9 "1#c0#u0#s"]]}"##)
    )]
    #[case::isoquinoline(
        "c1ccc2cnccc2c1",
        mol_dsl!(r##"{:aromatic-systems [{:atoms [0 1 2 3 4 5 6 7 8 9] :attrs "[1,1,1,1,1,1,1,1,1,1]#c0#u0#s"}] :atoms ["C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h0#n0#u0#s" "C#i=#c0#h#n0#u0#s" "N#i=#c0#h0#n#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h0#n0#u0#s" "C#i=#c0#h#n0#u0#s"] :bonds [[0 9 "1#c0#u0#s"] [0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"] [3 8 "1#c0#u0#s"] [3 4 "1#c0#u0#s"] [4 5 "1#c0#u0#s"] [5 6 "1#c0#u0#s"] [6 7 "1#c0#u0#s"] [7 8 "1#c0#u0#s"] [8 9 "1#c0#u0#s"]]}"##)
    )]
    #[case::chloronium(
        "C1C[Cl+]1",
        mol_dsl!(r##"{:atoms ["C#i=#c0#h2#n0#u0#s" "C#i=#c0#h2#n0#u0#s" "Cl#i=#c+#h0#n2#u0#s"] :bonds [[0 2 "1#c0#u0#s"] [0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"]]}"##)
    )]
    #[case::chlorine_trifluoride(
        "FCl(F)F",
        mol_dsl!(r##"{:atoms ["F#i=#c0#h0#n3#u0#s" "Cl#i=#c0#h0#n2#u0#s" "F#i=#c0#h0#n3#u0#s" "F#i=#c0#h0#n3#u0#s"] :bonds [[0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [1 3 "1#c0#u0#s"]]}"##)
    )]
    fn test_ingest_smiles_resolution(#[case] input: &str, #[case] expected: Molecule) {
        assert_eq!(ingest_smiles(input).unwrap(), expected);
    }

    #[rstest]
    #[case::imidazole(
        "c1cncn1",
        SmilesInputError::Underdetermined(ResolveUnderdetermined {
            report: ResolveReport {
                unresolved: AtomCompletions::from_iter([2, 4].map(|atom| (
                    AtomId(atom),
                    smallvec![
                        atom_dsl!("N#i=#c0#h0#n#u0#s#v2#a"),
                        atom_dsl!("N#i=#c0#h#n0#u0#s#v2#a2"),
                    ],
                ))),
                tie_breaks: Vec::new(),
            },
        })
    )]
    fn test_ingest_smiles_with_tie_break(#[case] input: &str, #[case] expected: SmilesInputError) {
        // Both tautomeric assignments survive `Strict`: the report carries the
        // two nitrogen splits.
        let model = ChemistryModel {
            valence: ValenceModel {
                tie_break: ValenceTieBreak::Strict,
                ..ValenceModel::smiles()
            },
            ..ChemistryModel::default()
        };
        assert_eq!(
            ingest_smiles_with(
                input,
                &SmilesIoConfig::opensmiles(),
                &model,
                &ResolveConfig::default(),
            ),
            Err(expected)
        );
    }

    #[rstest]
    fn test_ingest_smiles_components() {
        // Five bridged triazine rings: fifteen flexible nitrogens whose
        // assignment product exceeds the per-component bound as a whole
        // molecule but not per candidate-ring component.
        let molecule =
            ingest_smiles("C(c1ncncn1)(c1ncncn1)(c1ncncn1)CC(c1ncncn1)c1ncncn1").unwrap();
        let systems: Vec<Vec<AtomId>> = molecule
            .aromatic_systems()
            .iter()
            .map(|system| system.atom_ids().collect())
            .collect();
        assert_eq!(systems.len(), 5);
        assert!(molecule.is_concrete());
    }

    #[rstest]
    #[case::mdl_furan("o1cccc1")]
    #[case::mdl_thiophene("s1cccc1")]
    #[case::mdl_pyrrole("[nH]1cccc1")]
    fn test_ingest_smiles_with_aromaticity_policy(#[case] input: &str) {
        let model = ChemistryModel {
            aromaticity: AromaticityModel::mdl(),
            ..ChemistryModel::default()
        };
        let molecule = ingest_smiles_with(
            input,
            &SmilesIoConfig::opensmiles(),
            &model,
            &ResolveConfig {
                aromaticity: AromaticityResolveConfig {
                    aromatic_valence_failure: AromaticityFailurePolicy::Keep,
                    ..AromaticityResolveConfig::default()
                },
                stereo: StereoResolveConfig::default(),
            },
        )
        .unwrap();

        assert_eq!(
            molecule
                .atoms()
                .iter()
                .map(|atom| atom.attributes.constraints.aromatic_valence().cloned())
                .collect::<Vec<_>>(),
            vec![Some(AromaticValenceForm::Aromatic(NumForm::Undetermined)); 5]
        );
        assert_eq!(
            molecule
                .bonds()
                .iter()
                .map(|bond| bond.attributes.constraints.aromatic())
                .collect::<Vec<_>>(),
            vec![BooleanForm::Lit(true); 5]
        );
        assert_eq!(
            molecule.aromatic_systems().ids().collect::<Vec<_>>(),
            Vec::<AromaticSystemId>::new()
        );
    }

    #[rstest]
    #[case::retained(
        SmilesIoConfig::opensmiles(),
        ChemistryModel {
            valence: ValenceModel {
                tie_break: ValenceTieBreak::MostSaturated,
                ..ValenceModel::default()
            },
            ..ChemistryModel::default()
        },
        ResolveConfig::default(),
        None
    )]
    #[case::reset(
        SmilesIoConfig::opensmiles(),
        ChemistryModel {
            valence: ValenceModel {
                tie_break: ValenceTieBreak::MostSaturated,
                ..ValenceModel::default()
            },
            ..ChemistryModel::default()
        },
        ResolveConfig {
            aromaticity: AromaticityResolveConfig::default(),
            stereo: StereoResolveConfig {
                reset_stereo_constraints: true,
                ..StereoResolveConfig::default()
            },
        },
        None
    )]
    fn test_ingest_smiles_with_stereo(
        #[case] io_config: SmilesIoConfig,
        #[case] model: ChemistryModel,
        #[case] resolve_config: ResolveConfig,
        #[case] expected: Option<TetrahedralStereoForm>,
    ) {
        let molecule =
            ingest_smiles_with("C[C@H](N)O", &io_config, &model, &resolve_config).unwrap();

        assert_eq!(
            molecule
                .atom(AtomId(1))
                .attributes
                .constraints
                .tetrahedral_stereo()
                .cloned(),
            expected
        );
        assert!(molecule.stereo_atoms().is_at(AtomId(1)));
    }

    #[rstest]
    #[case::contradiction(
        "[nH]1cccc1",
        SmilesIoConfig::opensmiles(),
        ChemistryModel {
            aromaticity: AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Clar { ring_limits: RingLimits::default() }, tie_break: AromaticityTieBreak::Strict },
            ..ChemistryModel::default()
        },
        ResolveConfig::default(),
        SmilesInputError::Contradiction(ResolveContradiction::Aromaticity(
            AromaticityContradiction::ClarNonBenzenoid(String::from(
                "Clar model requires benzenoid input but non-carbon aromatic atoms are present",
            )),
        ))
    )]
    #[case::underdetermined(
        "C",
        SmilesIoConfig::opensmiles(),
        ChemistryModel {
            valence: ValenceModel::atom_typing(Cow::Owned(AtomTypeRegistry::from_atoms([atom_dsl!("C#c0")]))),
            ..ChemistryModel::default()
        },
        ResolveConfig::default(),
        SmilesInputError::Underdetermined(ResolveUnderdetermined::default())
    )]
    #[case::mdl_furan(
        "o1cccc1",
        SmilesIoConfig::opensmiles(),
        ChemistryModel {
            aromaticity: AromaticityModel::mdl(),
            ..ChemistryModel::default()
        },
        ResolveConfig::default(),
        SmilesInputError::Contradiction(ResolveContradiction::Aromaticity(
            AromaticityContradiction::Inconsistency(
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) },
            ),
        )),
    )]
    #[case::mdl_thiophene(
        "s1cccc1",
        SmilesIoConfig::opensmiles(),
        ChemistryModel {
            aromaticity: AromaticityModel::mdl(),
            ..ChemistryModel::default()
        },
        ResolveConfig::default(),
        SmilesInputError::Contradiction(ResolveContradiction::Aromaticity(
            AromaticityContradiction::Inconsistency(
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) },
            ),
        )),
    )]
    #[case::mdl_pyrrole(
        "[nH]1cccc1",
        SmilesIoConfig::opensmiles(),
        ChemistryModel {
            aromaticity: AromaticityModel::mdl(),
            ..ChemistryModel::default()
        },
        ResolveConfig::default(),
        SmilesInputError::Contradiction(ResolveContradiction::Aromaticity(
            AromaticityContradiction::Inconsistency(
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) },
            ),
        )),
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
                &ChemistryModel {
                    valence: ValenceModel::smiles(),
                    ..ChemistryModel::default()
                },
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
            MoleculeInterpretationError::Underdetermined(ResolveUnderdetermined::default()),
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
                &ChemistryModel {
                    valence: ValenceModel::smiles(),
                    ..ChemistryModel::default()
                },
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
    #[case::io("C~C>>C.C")]
    fn test_ingest_reaction_smiles_with_underdetermined_report(#[case] input: &str) {
        let result = ingest_reaction_smiles_with(
            input,
            &SmilesIoConfig::lenient(),
            &ChemistryModel::default(),
            &ResolveConfig::default(),
        );
        let Err(ReactionSmilesInputError::Interpretation(ReactionInterpretationError::Reactants(
            MoleculeInterpretationError::Underdetermined(underdetermined),
        ))) = result
        else {
            panic!("expected an underdetermined reactants interpretation: {result:?}");
        };
        assert!(!underdetermined.report.unresolved.is_empty());
    }

    #[rstest]
    #[case::chemistry(
        "[nH]1cccc1>>",
        SmilesIoConfig::opensmiles(),
        ChemistryModel {
            aromaticity: AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Clar { ring_limits: RingLimits::default() }, tie_break: AromaticityTieBreak::Strict },
            ..ChemistryModel::default()
        },
        ResolveConfig::default(),
        Err(ReactionSmilesInputError::Interpretation(
            ReactionInterpretationError::Reactants(
                MoleculeInterpretationError::Contradiction(
                    ResolveContradiction::Aromaticity(
                        AromaticityContradiction::ClarNonBenzenoid(String::from(
                            "Clar model requires benzenoid input but non-carbon aromatic atoms are present",
                        )),
                    ),
                ),
            ),
        )),
    )]
    #[case::mdl_furan(
        "o1cccc1>>C",
        SmilesIoConfig::opensmiles(),
        ChemistryModel {
            aromaticity: AromaticityModel::mdl(),
            ..ChemistryModel::default()
        },
        ResolveConfig::default(),
        Err(ReactionSmilesInputError::Interpretation(
            ReactionInterpretationError::Reactants(
                MoleculeInterpretationError::Contradiction(
                    ResolveContradiction::Aromaticity(
                        AromaticityContradiction::Inconsistency(
                            AromaticityInconsistency::AromaticValenceFailure {
                                atom: AtomId(0),
                            },
                        ),
                    ),
                ),
            ),
        )),
    )]
    #[case::mdl_thiophene(
        "s1cccc1>>C",
        SmilesIoConfig::opensmiles(),
        ChemistryModel {
            aromaticity: AromaticityModel::mdl(),
            ..ChemistryModel::default()
        },
        ResolveConfig::default(),
        Err(ReactionSmilesInputError::Interpretation(
            ReactionInterpretationError::Reactants(
                MoleculeInterpretationError::Contradiction(
                    ResolveContradiction::Aromaticity(
                        AromaticityContradiction::Inconsistency(
                            AromaticityInconsistency::AromaticValenceFailure {
                                atom: AtomId(0),
                            },
                        ),
                    ),
                ),
            ),
        )),
    )]
    #[case::mdl_pyrrole(
        "[nH]1cccc1>>C",
        SmilesIoConfig::opensmiles(),
        ChemistryModel {
            aromaticity: AromaticityModel::mdl(),
            ..ChemistryModel::default()
        },
        ResolveConfig::default(),
        Err(ReactionSmilesInputError::Interpretation(
            ReactionInterpretationError::Reactants(
                MoleculeInterpretationError::Contradiction(
                    ResolveContradiction::Aromaticity(
                        AromaticityContradiction::Inconsistency(
                            AromaticityInconsistency::AromaticValenceFailure {
                                atom: AtomId(0),
                            },
                        ),
                    ),
                ),
            ),
        )),
    )]
    #[case::resolve(
        "[cH+:1]1[cH:2][cH:3]1>>[cH+:1]1[cH:2][cH:3]1",
        SmilesIoConfig::opensmiles(),
        ChemistryModel::default(),
        ResolveConfig::default(),
        Ok(r##"{:deltas [] :lhs {:aromatic-systems [{:atoms [0 1 2] :attrs "[0,1,1]#c0#u0#s"}] :atoms ["C#i=#c+#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s"] :bonds [[0 2 "1#c0#u0#s"] [0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"]]}}"##.parse().unwrap()),
    )]
    fn test_ingest_reaction_smiles_with(
        #[case] input: &str,
        #[case] io_config: SmilesIoConfig,
        #[case] model: ChemistryModel,
        #[case] resolve_config: ResolveConfig,
        #[case] expected: Result<Reaction, ReactionSmilesInputError>,
    ) {
        assert_eq!(
            ingest_reaction_smiles_with(input, &io_config, &model, &resolve_config),
            expected
        );
    }

    #[rstest]
    #[case::mdl_benzene(
        "[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1>>[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1",
        AromaticityModel::mdl(),
        vec![1, 1, 1, 1, 1, 1],
    )]
    #[case::mdl_pyridine(
        "[n:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1>>[n:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1",
        AromaticityModel::mdl(),
        vec![1, 1, 1, 1, 1, 1],
    )]
    #[case::daylight_furan(
        "[o:1]1[cH:2][cH:3][cH:4][cH:5]1>>[o:1]1[cH:2][cH:3][cH:4][cH:5]1",
        AromaticityModel::daylight(),
        vec![2, 1, 1, 1, 1],
    )]
    #[case::daylight_thiophene(
        "[s:1]1[cH:2][cH:3][cH:4][cH:5]1>>[s:1]1[cH:2][cH:3][cH:4][cH:5]1",
        AromaticityModel::daylight(),
        vec![2, 1, 1, 1, 1],
    )]
    #[case::daylight_pyrrole(
        "[nH:1]1[cH:2][cH:3][cH:4][cH:5]1>>[nH:1]1[cH:2][cH:3][cH:4][cH:5]1",
        AromaticityModel::daylight(),
        vec![2, 1, 1, 1, 1],
    )]
    fn test_ingest_reaction_smiles_with_aromaticity(
        #[case] input: &str,
        #[case] aromaticity: AromaticityModel,
        #[case] expected_electrons: Vec<i64>,
    ) {
        let model = ChemistryModel {
            aromaticity,
            ..ChemistryModel::default()
        };
        let reaction = ingest_reaction_smiles_with(
            input,
            &SmilesIoConfig::opensmiles(),
            &model,
            &ResolveConfig::default(),
        )
        .unwrap();
        let system = reaction.lhs.aromatic_system(AromaticSystemId(0));

        assert_eq!(
            system.atom_ids().collect::<Vec<_>>(),
            (0..expected_electrons.len())
                .map(|index| AtomId(u32::try_from(index).unwrap()))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            system.electrons(),
            &ElectronCountsForm::Lit(expected_electrons)
        );
        assert_eq!(
            reaction.lhs.aromatic_systems().ids().collect::<Vec<_>>(),
            vec![AromaticSystemId(0)]
        );
        assert_eq!(reaction.deltas, Deltas::new());
    }

    #[rstest]
    #[case::mdl_furan("[o:1]1[cH:2][cH:3][cH:4][cH:5]1>>[o:1]1[cH:2][cH:3][cH:4][cH:5]1")]
    #[case::mdl_thiophene("[s:1]1[cH:2][cH:3][cH:4][cH:5]1>>[s:1]1[cH:2][cH:3][cH:4][cH:5]1")]
    #[case::mdl_pyrrole("[nH:1]1[cH:2][cH:3][cH:4][cH:5]1>>[nH:1]1[cH:2][cH:3][cH:4][cH:5]1")]
    fn test_ingest_reaction_smiles_with_aromaticity_policy(#[case] input: &str) {
        let model = ChemistryModel {
            aromaticity: AromaticityModel::mdl(),
            ..ChemistryModel::default()
        };
        let reaction = ingest_reaction_smiles_with(
            input,
            &SmilesIoConfig::opensmiles(),
            &model,
            &ResolveConfig {
                aromaticity: AromaticityResolveConfig {
                    aromatic_valence_failure: AromaticityFailurePolicy::Keep,
                    ..AromaticityResolveConfig::default()
                },
                stereo: StereoResolveConfig::default(),
            },
        )
        .unwrap();

        assert_eq!(
            reaction
                .lhs
                .atoms()
                .iter()
                .map(|atom| atom.attributes.constraints.aromatic_valence().cloned())
                .collect::<Vec<_>>(),
            vec![Some(AromaticValenceForm::Aromatic(NumForm::Undetermined)); 5]
        );
        assert_eq!(
            reaction
                .lhs
                .bonds()
                .iter()
                .map(|bond| bond.attributes.constraints.aromatic())
                .collect::<Vec<_>>(),
            vec![BooleanForm::Lit(true); 5]
        );
        assert_eq!(
            reaction.lhs.aromatic_systems().ids().collect::<Vec<_>>(),
            Vec::<AromaticSystemId>::new()
        );
        assert_eq!(reaction.deltas, Deltas::new());
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
