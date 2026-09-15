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
    IsotopePolicy, ResolveConfig, ResolveContradiction, ResolveError, ResolveUnderdetermined,
    Resolver,
};

/// Convert a parsed external-format value into a graph model.
pub trait Interpret {
    type Output;
    type Error;

    /// Interpret this format value under the semantic model and resolve policy.
    ///
    /// Format-owned values survive resolution: SMILES with an omitted isotope
    /// raises to Natural composition, which Strict also preserves.
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
/// valence preset and Natural isotope policy.
pub fn ingest_smiles(input: &str) -> Result<Molecule, SmilesInputError> {
    ingest_smiles_bytes(input.as_bytes())
}

/// Ingest SMILES bytes with the OpenSMILES configuration and the SMILES
/// valence preset and Natural isotope policy.
pub fn ingest_smiles_bytes(input: &[u8]) -> Result<Molecule, SmilesInputError> {
    ingest_smiles_bytes_with(
        input,
        &SmilesIoConfig::opensmiles(),
        &ChemistryModel {
            valence: ValenceModel::smiles(),
            ..ChemistryModel::default()
        },
        &ResolveConfig {
            isotope: IsotopePolicy::Natural,
            ..Default::default()
        },
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
/// SMILES valence preset and Natural isotope policy.
pub fn ingest_reaction_smiles(input: &str) -> Result<Reaction, ReactionSmilesInputError> {
    ingest_reaction_smiles_bytes(input.as_bytes())
}

/// Ingest reaction SMILES bytes with the OpenSMILES configuration and the
/// SMILES valence preset and Natural isotope policy.
pub fn ingest_reaction_smiles_bytes(input: &[u8]) -> Result<Reaction, ReactionSmilesInputError> {
    ingest_reaction_smiles_bytes_with(
        input,
        &SmilesIoConfig::opensmiles(),
        &ChemistryModel {
            valence: ValenceModel::smiles(),
            ..ChemistryModel::default()
        },
        &ResolveConfig {
            isotope: IsotopePolicy::Natural,
            ..Default::default()
        },
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
    use umol_chem::element::Element;
    use umol_graph_core::{NodeId, Remapping};
    use umol_graph_ir::ir::{
        AromaticSystemId, AromaticValenceForm, AtomForm, AtomId, BondConstraintForm, BondId,
        BooleanForm, Canonicalize, CanonicalizeContext, Constraint, Deltas, ElectronCountsForm,
        Entity, MoleculeEntries, MoleculeIntegrityError, NumForm, StereoAtomId, StereoBondForm,
        StereoCoset, StereoKind, StereoLigand, StereoLigandKind, TetrahedralStereoForm,
    };
    use umol_graph_ir::{atom_dsl, mol_dsl, mol_dsl_concrete};
    use umol_io::table_ir::{AtomPair, BondConfiguration};

    use super::*;
    use crate::ops::aromaticity::{
        AromaticityContradiction, AromaticityError, AromaticityInconsistency,
    };
    use crate::ops::model::{
        AromaticityModel, AromaticityRule, AromaticityTieBreak, ElementScope, ValenceModel,
        ValenceTieBreak,
    };
    use crate::ops::resolve::{
        AromaticityFailurePolicy, AromaticityResolveConfig, DischargeContradiction,
        StereoContradiction, StereoResolveConfig, ValenceContradiction,
    };
    use crate::ops::stereo::StereoInconsistency;
    use crate::ops::valence::{AtomCompletions, AtomTypeRegistry, AtomTypingError, ResolveReport};

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
    #[case::omitted("[CH4]", mol_dsl!(r#"{:atoms ["C#i=#c0#h4#n0#u0#s"]}"#))]
    #[case::mass("[13CH4]", mol_dsl!(r#"{:atoms ["C#i13#c0#h4#n0#u0#s"]}"#))]
    fn test_smiles_interpret_isotope(
        #[values(ValenceModel::smiles(), ValenceModel::default())] valence: ValenceModel,
        #[values(IsotopePolicy::Strict, IsotopePolicy::Natural)] isotope: IsotopePolicy,
        #[case] source: &str,
        #[case] molecule: Molecule,
    ) {
        let model = ChemistryModel {
            valence,
            ..Default::default()
        };
        let config = ResolveConfig {
            isotope,
            ..Default::default()
        };
        let expected: Result<_, MoleculeInterpretationError> = Ok(molecule);
        assert_eq!(
            Smiles::parse(source).unwrap().interpret(&model, &config),
            expected
        );
        let expected = expected.map_err(SmilesInputError::from);
        assert_eq!(
            ingest_smiles_with(source, &SmilesIoConfig::opensmiles(), &model, &config),
            expected
        );
        assert_eq!(
            ingest_smiles_bytes_with(
                source.as_bytes(),
                &SmilesIoConfig::opensmiles(),
                &model,
                &config
            ),
            expected
        );
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
    #[case::omitted("[CH4:1]>>[CH4:1]", "C#i=#c0#h4#n0#u0#s")]
    #[case::mass("[13CH4:1]>>[13CH4:1]", "C#i13#c0#h4#n0#u0#s")]
    fn test_reaction_smiles_interpret_isotope(
        #[values(IsotopePolicy::Strict, IsotopePolicy::Natural)] isotope: IsotopePolicy,
        #[case] source: &str,
        #[case] atom: &str,
    ) {
        let model = ChemistryModel {
            valence: ValenceModel::smiles(),
            ..Default::default()
        };
        let config = ResolveConfig {
            isotope,
            ..Default::default()
        };
        let expected: Result<_, ReactionInterpretationError> = Ok(Reaction::new(
            Molecule::from_entries(MoleculeEntries {
                atoms: vec![atom.parse().unwrap()],
                ..Default::default()
            }),
            Deltas::default(),
        ));
        assert_eq!(
            ReactionSmiles::parse(source)
                .unwrap()
                .interpret(&model, &config),
            expected
        );
        let expected = expected.map_err(ReactionSmilesInputError::from);
        assert_eq!(
            ingest_reaction_smiles_with(source, &SmilesIoConfig::opensmiles(), &model, &config),
            expected
        );
        assert_eq!(
            ingest_reaction_smiles_bytes_with(
                source.as_bytes(),
                &SmilesIoConfig::opensmiles(),
                &model,
                &config
            ),
            expected
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
                RaiseError::MoleculeEntries(MoleculeIntegrityError::StereoLigandArity { entity: Entity::StereoAtom(StereoAtomId(0)), kind: StereoKind::Tetrahedral, expected: 4, actual: 2 }),
            ),
        ),
    )]
    #[case::products_model_conversion(
        ">>C[S@]C",
        ChemistryModel::default(),
        ReactionInterpretationError::Products(
            MoleculeInterpretationError::ModelConversion(
                RaiseError::MoleculeEntries(MoleculeIntegrityError::StereoLigandArity { entity: Entity::StereoAtom(StereoAtomId(0)), kind: StereoKind::Tetrahedral, expected: 4, actual: 2 }),
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
            aromaticity: AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Clar, tie_break: AromaticityTieBreak::Strict },
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
            aromaticity: AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Clar, tie_break: AromaticityTieBreak::Strict },
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
                    "C#c0#h0#n0#u0#s#v2#a2"
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
                    "C#c0#h0#n0#u0#s#v2#a2"
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
    #[case::opening("C[C@H]1CCCCO1")]
    #[case::closing("O1CCCC[C@@H]1C")]
    #[case::lone_pair("C[S@](=O)CC")]
    #[case::alkene("F/C=C/Cl")]
    #[case::ring_direction("C/C=C1CO\\1")]
    #[case::conjugated("C/C=C/C=C/C")]
    fn test_interpret_molecule_bond_storage(#[case] input: &str) {
        let table = Smiles::parse(input).unwrap().into_table_ir();
        let expected = ingest_smiles(input).unwrap();
        let model = ChemistryModel {
            valence: ValenceModel::smiles(),
            ..ChemistryModel::default()
        };
        let context = CanonicalizeContext {
            para_stereo: false,
            automorphism_algorithm: umol_graph_core::AutomorphismAlgorithm::Nauty,
        };
        for is_reverse in [true, false] {
            let mut reordered = table.clone();
            if is_reverse {
                reordered.bonds.reverse();
                for frame in &mut reordered.stereo_bonds {
                    frame.bond = reordered.bonds.len() as u32 - 1 - frame.bond;
                }
            } else {
                reordered.bonds.rotate_left(1);
                for frame in &mut reordered.stereo_bonds {
                    frame.bond = (frame.bond + reordered.bonds.len() as u32 - 1)
                        % reordered.bonds.len() as u32;
                }
            }
            let actual = interpret_molecule(&reordered, &model, &ResolveConfig::default()).unwrap();
            assert!(actual.canonical_eq(&expected, &context));
        }
    }

    #[rstest]
    #[case::one_direction("F/C=C/Cl", vec![1, 0, 2, 3])]
    #[case::branched("F/C(Cl)=C/Br", vec![4, 3, 2, 1, 0])]
    fn test_interpret_molecule_bond_endpoints(#[case] input: &str, #[case] atoms: Vec<usize>) {
        let mut table = Smiles::parse(input).unwrap().into_table_ir();
        let expected = ingest_smiles(input).unwrap();
        let atoms = Remapping::new(atoms.into_iter().map(NodeId::from).collect()).unwrap();
        for frame in &mut table.stereo_bonds {
            if let BondConfiguration::Framed { references, .. } = &mut frame.configuration {
                let pair = table.bonds[frame.bond as usize].atoms;
                *references = references.map(|atom| atoms.map(NodeId(atom)).0);
                if atoms.map(NodeId(pair.first())) > atoms.map(NodeId(pair.second())) {
                    references.swap(0, 1);
                }
            }
        }
        table.atoms = atoms.remap_vec(table.atoms);
        for bond in &mut table.bonds {
            let first = atoms.map(NodeId(bond.atoms.first())).0;
            let second = atoms.map(NodeId(bond.atoms.second())).0;
            bond.atoms = AtomPair::new(first, second);
        }
        let model = ChemistryModel {
            valence: ValenceModel::smiles(),
            ..ChemistryModel::default()
        };
        let actual = interpret_molecule(&table, &model, &ResolveConfig::default()).unwrap();
        let context = CanonicalizeContext {
            para_stereo: false,
            automorphism_algorithm: umol_graph_core::AutomorphismAlgorithm::Nauty,
        };
        assert!(actual.canonical_eq(&expected, &context));
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
        SmilesInputError::ModelConversion(RaiseError::MoleculeEntries(MoleculeIntegrityError::StereoLigandArity { entity: Entity::StereoAtom(StereoAtomId(0)), kind: StereoKind::Tetrahedral, expected: 4, actual: 2 }))
    )]
    #[case::underdetermined(
        "*",
        SmilesInputError::Underdetermined(ResolveUnderdetermined::default())
    )]
    #[case::trigonal_carbon(
        "[C@](F)(Cl)Br",
        SmilesInputError::Contradiction(ResolveContradiction::Stereo(
            StereoContradiction::Inconsistency(StereoInconsistency::StereoAtomFailure {
                stereo_atom: StereoAtomId(0),
            }),
        ))
    )]
    #[case::trigonal_carbocation(
        "[C@+](F)(Cl)Br",
        SmilesInputError::Contradiction(ResolveContradiction::Stereo(
            StereoContradiction::Inconsistency(StereoInconsistency::StereoAtomFailure {
                stereo_atom: StereoAtomId(0),
            }),
        ))
    )]
    #[case::two_carbon_substituents(
        "[C@](F)Cl",
        SmilesInputError::ModelConversion(RaiseError::MoleculeEntries(
            MoleculeIntegrityError::StereoLigandArity {
                entity: Entity::StereoAtom(StereoAtomId(0)),
                kind: StereoKind::Tetrahedral, expected: 4, actual: 2,
            },
        ))
    )]
    #[case::one_carbon_substituent(
        "[C@]F",
        SmilesInputError::ModelConversion(RaiseError::MoleculeEntries(
            MoleculeIntegrityError::StereoLigandArity {
                entity: Entity::StereoAtom(StereoAtomId(0)),
                kind: StereoKind::Tetrahedral, expected: 4, actual: 1,
            },
        ))
    )]
    #[case::five_ligands(
        "[C@H](F)(Cl)(Br)I",
        SmilesInputError::ModelConversion(RaiseError::MoleculeEntries(
            MoleculeIntegrityError::StereoLigandArity {
                entity: Entity::StereoAtom(StereoAtomId(0)),
                kind: StereoKind::Tetrahedral, expected: 4, actual: 5,
            },
        ))
    )]
    #[case::duplicate_bracket_h(
        "[C@H2](F)(Cl)Br",
        SmilesInputError::ModelConversion(RaiseError::MoleculeEntries(
            MoleculeIntegrityError::DuplicateStereoLigand {
                entity: Entity::StereoAtom(StereoAtomId(0)),
                ligand: StereoLigand::new(AtomId(0), StereoLigandKind::ImplicitHydrogen),
            },
        ))
    )]
    #[case::biphenyl_unwritten_single(
        // OpenSMILES: a single (nonaromatic) bond between two aromatic atoms
        // must be explicitly represented; the unwritten bond raises as an
        // aromatic assertion no system can claim.
        "c1ccc(cc1)c1ccccc1",
        SmilesInputError::Contradiction(ResolveContradiction::Discharge(
            DischargeContradiction::Assertion {
                constraint: Constraint::Bond(
                    BondId(6),
                    BondConstraintForm::Aromatic(BooleanForm::Lit(true)),
                ),
            },
        ))
    )]
    #[case::cyclobutadiene(
        "c1ccc1",
        SmilesInputError::Contradiction(ResolveContradiction::Aromaticity(
            AromaticityContradiction::Inconsistency(
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) },
            ),
        ))
    )]
    #[case::dangling_direction("F/C=C", SmilesInputError::Syntax(SmilesParseError::DanglingBondDirection { bond: 0 }))]
    #[case::conflicting_direction("F/C(\\Cl)=CF", SmilesInputError::Syntax(SmilesParseError::CisTransConflict { atom: 1 }))]
    #[case::ring_direction("C/1CC/1", SmilesInputError::Syntax(SmilesParseError::MismatchedRingBondDirections { pos: 6, open_pos: 2 }))]
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
    #[case::equivalent("C[C@H]1CCCCO1", "O1CCCC[C@@H]1C", true)]
    #[case::mirror("C[C@H]1CCCCO1", "O1CCCC[C@H]1C", false)]
    #[case::cis_dichlorocyclohexane("Cl[C@H]1CCCC[C@H]1Cl", "Cl[C@@H]1CCCC[C@@H]1Cl", true)]
    #[case::trans_dichlorocyclohexane("Cl[C@H]1CCCC[C@@H]1Cl", "Cl[C@@H]1CCCC[C@H]1Cl", false)]
    #[case::glucose_epimer(
        "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
        "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O",
        false
    )]
    #[case::branch_order("N[C@H](F)Cl", "N[C@@H](Cl)F", true)]
    #[case::root_hydrogen("[C@H](F)(Cl)Br", "F[C@@H](Cl)Br", true)]
    #[case::lone_pair("C[S@](=O)CC", "C[S@@](CC)=O", true)]
    #[case::ring_label("C[C@H]1CCCCO1", "C[C@H]%99CCCCO%99", true)]
    #[case::both_signs("F/C=C/Cl", "F\\C=C\\Cl", true)]
    #[case::opposite_alkene("F/C=C/Cl", "F/C=C\\Cl", false)]
    #[case::bond_traversal("F/C=C/Cl", "Cl/C=C/F", true)]
    #[case::marked_substituent("F/C(Cl)=C/Br", "FC(/Cl)=C/Br", true)]
    #[case::bond_branch_order("F/C(Cl)=C/Br", "Cl/C(F)=C\\Br", true)]
    #[case::ring_marker("C/C=C1CO\\1", "C/C=C/1CO1", true)]
    #[case::conjugated("C/C=C/C=C/C", "C\\C=C\\C=C\\C", true)]
    #[case::conjugated_opposite("C/C=C/C=C/C", "C/C=C/C=C\\C", false)]
    #[case::explicit_hydrogen("C[C@]1([H])CCCCO1", "C[C@@]1(CCCCO1)[H]", true)]
    #[case::hydrogen_representation("[C@H](F)(Cl)Br", "[C@]([H])(F)(Cl)Br", false)]
    #[case::later_root("C.[C@H](F)(Cl)Br", "C.F[C@@H](Cl)Br", true)]
    #[case::mixed_digits("O1CCC[C@]21CCNC2", "O1CCC[C@@]12CCNC2", true)]
    fn test_ingest_smiles_stereo(#[case] left: &str, #[case] right: &str, #[case] expected: bool) {
        let left = ingest_smiles(left).unwrap();
        let right = ingest_smiles(right).unwrap();
        let context = CanonicalizeContext {
            para_stereo: false,
            automorphism_algorithm: umol_graph_core::AutomorphismAlgorithm::Nauty,
        };
        assert_eq!(left.canonical_eq(&right, &context), expected);
    }

    #[rstest]
    #[case::e("F/C=C/Cl", vec![(1, StereoCoset::Lit(1))])]
    #[case::z("F/C=C\\Cl", vec![(1, StereoCoset::Lit(0))])]
    #[case::unmarked("FC=CCl", vec![])]
    #[case::one_sided("F/C=CCl", vec![])]
    #[case::conjugated("C/C=C/C=C/C", vec![(1, StereoCoset::Lit(1)), (3, StereoCoset::Lit(1))])]
    fn test_ingest_smiles_bond_stereo(
        #[case] input: &str,
        #[case] expected: Vec<(usize, StereoCoset)>,
    ) {
        let molecule = ingest_smiles(input).unwrap();
        assert_eq!(
            molecule
                .stereo_bonds()
                .iter()
                .map(|bond| (usize::from(bond.site_id()), bond.attributes.clone()))
                .collect::<Vec<_>>(),
            expected
                .into_iter()
                .map(|(site, coset)| (site, StereoBondForm::new(StereoKind::CisTrans, coset)))
                .collect::<Vec<_>>()
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
            isotope: IsotopePolicy::Strict,
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
    #[case::borepin(
        "b1cccccc1",
        mol_dsl!(r##"{:aromatic-systems [{:atoms [0 1 2 3 4 5 6] :attrs "[0,1,1,1,1,1,1]#c0#u0#s"}] :atoms ["B#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s"] :bonds [[0 6 "1#c0#u0#s"] [0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"] [3 4 "1#c0#u0#s"] [4 5 "1#c0#u0#s"] [5 6 "1#c0#u0#s"]]}"##)
    )]
    #[case::borazine(
        "b1nbnbn1",
        mol_dsl!(r##"{:aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[0,2,0,2,0,2]#c0#u0#s"}] :atoms ["B#i=#c0#h#n0#u0#s" "N#i=#c0#h#n0#u0#s" "B#i=#c0#h#n0#u0#s" "N#i=#c0#h#n0#u0#s" "B#i=#c0#h#n0#u0#s" "N#i=#c0#h#n0#u0#s"] :bonds [[0 5 "1#c0#u0#s"] [0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"] [3 4 "1#c0#u0#s"] [4 5 "1#c0#u0#s"]]}"##)
    )]
    #[case::azaborine(
        "b1ccccn1",
        mol_dsl!(r##"{:aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[0,1,1,1,1,2]#c0#u0#s"}] :atoms ["B#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "N#i=#c0#h#n0#u0#s"] :bonds [[0 5 "1#c0#u0#s"] [0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"] [3 4 "1#c0#u0#s"] [4 5 "1#c0#u0#s"]]}"##)
    )]
    fn test_ingest_smiles_with_element_scope(#[case] input: &str, #[case] expected: Molecule) {
        // Boron rings resolve once the scope admits the element; the refusal
        // twin below pins the same input dying on scope, not zero handling.
        let model = ChemistryModel {
            valence: ValenceModel::smiles(),
            aromaticity: AromaticityModel::permissive(),
            ..ChemistryModel::default()
        };
        assert_eq!(
            ingest_smiles_with(
                input,
                &SmilesIoConfig::opensmiles(),
                &model,
                &ResolveConfig::default(),
            )
            .unwrap(),
            expected
        );
    }

    #[rstest]
    #[case::borepin(
        "b1cccccc1",
        SmilesInputError::Contradiction(ResolveContradiction::Aromaticity(
            AromaticityContradiction::Inconsistency(
                AromaticityInconsistency::AromaticValenceFailure { atom: AtomId(0) },
            ),
        ))
    )]
    fn test_ingest_smiles_with_element_scope_error(
        #[case] input: &str,
        #[case] expected: SmilesInputError,
    ) {
        let model = ChemistryModel {
            valence: ValenceModel::smiles(),
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
    #[case::pteridine(
        "c1cnc2c(n1)cncn2",
        mol_dsl!(r##"{:aromatic-systems [{:atoms [0 1 2 3 4 5 6 7 8 9] :attrs "[1,1,1,1,1,1,1,1,1,1]#c0#u0#s"}] :atoms ["C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "N#i=#c0#h0#n#u0#s" "C#i=#c0#h0#n0#u0#s" "C#i=#c0#h0#n0#u0#s" "N#i=#c0#h0#n#u0#s" "C#i=#c0#h#n0#u0#s" "N#i=#c0#h0#n#u0#s" "C#i=#c0#h#n0#u0#s" "N#i=#c0#h0#n#u0#s"] :bonds [[0 5 "1#c0#u0#s"] [0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"] [3 9 "1#c0#u0#s"] [3 4 "1#c0#u0#s"] [4 5 "1#c0#u0#s"] [4 6 "1#c0#u0#s"] [6 7 "1#c0#u0#s"] [7 8 "1#c0#u0#s"] [8 9 "1#c0#u0#s"]]}"##)
    )]
    #[case::isoquinoline(
        "c1ccc2cnccc2c1",
        mol_dsl!(r##"{:aromatic-systems [{:atoms [0 1 2 3 4 5 6 7 8 9] :attrs "[1,1,1,1,1,1,1,1,1,1]#c0#u0#s"}] :atoms ["C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h0#n0#u0#s" "C#i=#c0#h#n0#u0#s" "N#i=#c0#h0#n#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h0#n0#u0#s" "C#i=#c0#h#n0#u0#s"] :bonds [[0 9 "1#c0#u0#s"] [0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"] [3 8 "1#c0#u0#s"] [3 4 "1#c0#u0#s"] [4 5 "1#c0#u0#s"] [5 6 "1#c0#u0#s"] [6 7 "1#c0#u0#s"] [7 8 "1#c0#u0#s"] [8 9 "1#c0#u0#s"]]}"##)
    )]
    #[case::biphenyl_explicit_single(
        "c1ccc(cc1)-c1ccccc1",
        mol_dsl!(r##"{:aromatic-systems [{:atoms [0 1 2 3 4 5] :attrs "[1,1,1,1,1,1]#c0#u0#s"} {:atoms [6 7 8 9 10 11] :attrs "[1,1,1,1,1,1]#c0#u0#s"}] :atoms ["C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h0#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h0#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s"] :bonds [[0 5 "1#c0#u0#s"] [0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"] [3 4 "1#c0#u0#s"] [4 5 "1#c0#u0#s"] [3 6 "1#c0#u0#s"] [6 11 "1#c0#u0#s"] [6 7 "1#c0#u0#s"] [7 8 "1#c0#u0#s"] [8 9 "1#c0#u0#s"] [9 10 "1#c0#u0#s"] [10 11 "1#c0#u0#s"]]}"##)
    )]
    #[case::pyridone_exocyclic(
        "O=c1cccc[nH]1",
        mol_dsl!(r##"{:aromatic-systems [{:atoms [1 2 3 4 5 6] :attrs "[0,1,1,1,1,2]#c0#u0#s"}] :atoms ["O#i=#c0#h0#n2#u0#s" "C#i=#c0#h0#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "N#i=#c0#h#n0#u0#s"] :bonds [[0 1 "2#c0#u0#s"] [1 6 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"] [3 4 "1#c0#u0#s"] [4 5 "1#c0#u0#s"] [5 6 "1#c0#u0#s"]]}"##)
    )]
    #[case::tropone_exocyclic(
        "O=c1cccccc1",
        mol_dsl!(r##"{:aromatic-systems [{:atoms [1 2 3 4 5 6 7] :attrs "[0,1,1,1,1,1,1]#c0#u0#s"}] :atoms ["O#i=#c0#h0#n2#u0#s" "C#i=#c0#h0#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s" "C#i=#c0#h#n0#u0#s"] :bonds [[0 1 "2#c0#u0#s"] [1 7 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"] [3 4 "1#c0#u0#s"] [4 5 "1#c0#u0#s"] [5 6 "1#c0#u0#s"] [6 7 "1#c0#u0#s"]]}"##)
    )]
    #[case::hydrogen_atom(
        "[H]",
        mol_dsl!(r##"{:atoms ["H#i=#c0#h0#n0#u#s2"]}"##)
    )]
    #[case::methyl_radical(
        "[CH3]",
        mol_dsl!(r##"{:atoms ["C#i=#c0#h3#n0#u#s2"]}"##)
    )]
    #[case::methane_bracket(
        "[CH4]",
        mol_dsl!(r##"{:atoms ["C#i=#c0#h4#n0#u0#s"]}"##)
    )]
    #[case::chloronium(
        "C1C[Cl+]1",
        mol_dsl!(r##"{:atoms ["C#i=#c0#h2#n0#u0#s" "C#i=#c0#h2#n0#u0#s" "Cl#i=#c+#h0#n2#u0#s"] :bonds [[0 2 "1#c0#u0#s"] [0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"]]}"##)
    )]
    #[case::chlorine_trifluoride(
        "FCl(F)F",
        mol_dsl!(r##"{:atoms ["F#i=#c0#h0#n3#u0#s" "Cl#i=#c0#h0#n2#u0#s" "F#i=#c0#h0#n3#u0#s" "F#i=#c0#h0#n3#u0#s"] :bonds [[0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [1 3 "1#c0#u0#s"]]}"##)
    )]
    #[case::ring_opening_frame(
        "C[C@H]1CCCCO1",
        mol_dsl!(r##"{:atoms ["C#i=#c0#h3#n0#u0#s" "C#i=#c0#h1#n0#u0#s" "C#i=#c0#h2#n0#u0#s" "C#i=#c0#h2#n0#u0#s" "C#i=#c0#h2#n0#u0#s" "C#i=#c0#h2#n0#u0#s" "O#i=#c0#h0#n2#u0#s"] :bonds [[0 1 "1#c0#u0#s"] [1 6 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"] [3 4 "1#c0#u0#s"] [4 5 "1#c0#u0#s"] [5 6 "1#c0#u0#s"]] :stereo-atoms [{:site 1 :ligands [0 [:h 1] 6 2] :attrs :ccw}]}"##)
    )]
    #[case::ring_closing_frame(
        "O1CCCC[C@@H]1C",
        mol_dsl!(r##"{:atoms ["O#i=#c0#h0#n2#u0#s" "C#i=#c0#h2#n0#u0#s" "C#i=#c0#h2#n0#u0#s" "C#i=#c0#h2#n0#u0#s" "C#i=#c0#h2#n0#u0#s" "C#i=#c0#h1#n0#u0#s" "C#i=#c0#h3#n0#u0#s"] :bonds [[0 5 "1#c0#u0#s"] [0 1 "1#c0#u0#s"] [1 2 "1#c0#u0#s"] [2 3 "1#c0#u0#s"] [3 4 "1#c0#u0#s"] [4 5 "1#c0#u0#s"] [5 6 "1#c0#u0#s"]] :stereo-atoms [{:site 5 :ligands [4 [:h 5] 0 6] :attrs :cw}]}"##)
    )]
    #[case::tetrahedral_four_atoms(
        "[C@](F)(Cl)(Br)I",
        mol_dsl!(r#"{
            :atoms ["C#i=#c0#h0#n0#u0#s" "F#i=#c0#h0#n3#u0#s"
                    "Cl#i=#c0#h0#n3#u0#s" "Br#i=#c0#h0#n3#u0#s" "I#i=#c0#h0#n3#u0#s"]
            :bonds [[0 1 "1#c0#u0#s"] [0 2 "1#c0#u0#s"] [0 3 "1#c0#u0#s"] [0 4 "1#c0#u0#s"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs :ccw}]
        }"#)
    )]
    #[case::tetrahedral_bracket_h(
        "[C@H](F)(Cl)Br",
        mol_dsl!(r#"{
            :atoms ["C#i=#c0#h1#n0#u0#s" "F#i=#c0#h0#n3#u0#s"
                    "Cl#i=#c0#h0#n3#u0#s" "Br#i=#c0#h0#n3#u0#s"]
            :bonds [[0 1 "1#c0#u0#s"] [0 2 "1#c0#u0#s"] [0 3 "1#c0#u0#s"]]
            :stereo-atoms [{:site 0 :ligands [[:h 0] 1 2 3] :attrs :ccw}]
        }"#)
    )]
    #[case::tetrahedral_explicit_h(
        "[C@]([H])(F)(Cl)Br",
        mol_dsl!(r#"{
            :atoms ["C#i=#c0#h0#n0#u0#s" "H#i=#c0#h0#n0#u0#s"
                    "F#i=#c0#h0#n3#u0#s" "Cl#i=#c0#h0#n3#u0#s" "Br#i=#c0#h0#n3#u0#s"]
            :bonds [[0 1 "1#c0#u0#s"] [0 2 "1#c0#u0#s"] [0 3 "1#c0#u0#s"] [0 4 "1#c0#u0#s"]]
            :stereo-atoms [{:site 0 :ligands [1 2 3 4] :attrs :ccw}]
        }"#)
    )]
    #[case::carbanion_lone_pair(
        "[C@-](F)(Cl)Br",
        mol_dsl!(r#"{
            :atoms ["C#i=#c-#h0#n1#u0#s" "F#i=#c0#h0#n3#u0#s"
                    "Cl#i=#c0#h0#n3#u0#s" "Br#i=#c0#h0#n3#u0#s"]
            :bonds [[0 1 "1#c0#u0#s"] [0 2 "1#c0#u0#s"] [0 3 "1#c0#u0#s"]]
            :stereo-atoms [{:site 0 :ligands [[:lp 0] 1 2 3] :attrs :ccw}]
        }"#)
    )]
    #[case::sulfoxide_lone_pair(
        "[S@](=O)(C)CC",
        mol_dsl!(r#"{
            :atoms ["S#i=#c0#h0#n1#u0#s" "O#i=#c0#h0#n2#u0#s"
                    "C#i=#c0#h3#n0#u0#s" "C#i=#c0#h2#n0#u0#s" "C#i=#c0#h3#n0#u0#s"]
            :bonds [[0 1 "2#c0#u0#s"] [0 2 "1#c0#u0#s"] [0 3 "1#c0#u0#s"] [3 4 "1#c0#u0#s"]]
            :stereo-atoms [{:site 0 :ligands [[:lp 0] 1 2 3] :attrs :ccw}]
        }"#)
    )]
    #[case::sulfoxide_zwitterion(
        "[S@+]([O-])(C)CC",
        mol_dsl!(r#"{
            :atoms ["S#i=#c+#h0#n1#u0#s" "O#i=#c-#h0#n3#u0#s"
                    "C#i=#c0#h3#n0#u0#s" "C#i=#c0#h2#n0#u0#s" "C#i=#c0#h3#n0#u0#s"]
            :bonds [[0 1 "1#c0#u0#s"] [0 2 "1#c0#u0#s"] [0 3 "1#c0#u0#s"] [3 4 "1#c0#u0#s"]]
            :stereo-atoms [{:site 0 :ligands [[:lp 0] 1 2 3] :attrs :ccw}]
        }"#)
    )]
    #[case::alkene_opposite(
        "F/C=C/Cl",
        mol_dsl!(r#"{
            :atoms ["F#i=#c0#h0#n3#u0#s" "C#i=#c0#h1#n0#u0#s"
                    "C#i=#c0#h1#n0#u0#s" "Cl#i=#c0#h0#n3#u0#s"]
            :bonds [[0 1 "1#c0#u0#s"] [1 2 "2#c0#u0#s"] [2 3 "1#c0#u0#s"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct1"}]
        }"#)
    )]
    #[case::alkene_same(
        "F/C=C\\Cl",
        mol_dsl!(r#"{
            :atoms ["F#i=#c0#h0#n3#u0#s" "C#i=#c0#h1#n0#u0#s"
                    "C#i=#c0#h1#n0#u0#s" "Cl#i=#c0#h0#n3#u0#s"]
            :bonds [[0 1 "1#c0#u0#s"] [1 2 "2#c0#u0#s"] [2 3 "1#c0#u0#s"]]
            :stereo-bonds [{:site 1 :ligands [0 [:h 1] 3 [:h 2]] :attrs "Ct0"}]
        }"#)
    )]
    fn test_ingest_smiles_resolution(#[case] input: &str, #[case] expected: Molecule) {
        assert_eq!(ingest_smiles(input).unwrap(), expected);
    }

    #[rstest]
    #[case::branched("[C]([CH3])([CH3])([OH])[F]", mol_dsl_concrete!(r#"{:atoms ["C" "C#h3" "C#h3" "O#h1#n2" "F#n3"] :bonds [[0 1 "1"] [0 2 "1"] [0 3 "1"] [0 4 "1"]]}"#))]
    #[case::ring("[CH2]1[CH2][O]1", mol_dsl_concrete!(r#"{:atoms ["C#h2" "C#h2" "O#n2"] :bonds [[0 2 "1"] [0 1 "1"] [1 2 "1"]]}"#))]
    #[case::nitrile("[CH3][C]#[N]", mol_dsl_concrete!(r#"{:atoms ["C#h3" "C" "N#n1"] :bonds [[0 1 "1"] [1 2 "3"]]}"#))]
    #[case::amide("[CH3][C](=[O])[NH2]", mol_dsl_concrete!(r#"{:atoms ["C#h3" "C" "O#n2" "N#h2#n1"] :bonds [[0 1 "1"] [1 2 "2"] [1 3 "1"]]}"#))]
    #[case::sulfoxide("[CH3][S](=[O])[CH3]", mol_dsl_concrete!(r#"{:atoms ["C#h3" "S#n1" "O#n2" "C#h3"] :bonds [[0 1 "1"] [1 2 "2"] [1 3 "1"]]}"#))]
    #[case::zwitterion("[NH3+][CH2][C](=[O])[O-]", mol_dsl_concrete!(r#"{:atoms ["N#c+#h3" "C#h2" "C" "O#n2" "O#c-#n3"] :bonds [[0 1 "1"] [1 2 "1"] [2 3 "2"] [2 4 "1"]]}"#))]
    #[case::components("[NH4+].[Cl-].[13CH3]", mol_dsl_concrete!(r#"{:atoms ["N#c+#h4" "Cl#c-#n4" "C#i13#h3#u1#s2"]}"#))]
    fn test_ingest_smiles_with_resolution(
        #[values(ValenceModel::smiles(), ValenceModel::default())] mut valence: ValenceModel,
        #[values(ValenceTieBreak::Strict, ValenceTieBreak::MostSaturated)]
        tie_break: ValenceTieBreak,
        #[case] input: &str,
        #[case] expected: Molecule,
    ) {
        valence.tie_break = tie_break;
        let model = ChemistryModel {
            valence,
            ..Default::default()
        };
        assert_eq!(
            ingest_smiles_with(
                input,
                &SmilesIoConfig::default(),
                &model,
                &ResolveConfig::default()
            ),
            Ok(expected)
        );
    }

    #[rstest]
    #[case::carbanion("[C@-](F)(Cl)Br", atom_dsl!("C#i=#c-#h0#n1#u0#s"))]
    #[case::nitrogen("[N@](F)(Cl)Br", atom_dsl!("N#i=#c0#h0#n1#u0#s"))]
    #[case::phosphorus("[P@](F)(Cl)Br", atom_dsl!("P#i=#c0#h0#n1#u0#s"))]
    #[case::sulfonium("[S@+](F)(Cl)Br", atom_dsl!("S#i=#c+#h0#n1#u0#s"))]
    fn test_ingest_smiles_with_lone_pair(
        #[case] input: &str,
        #[case] expected: AtomForm,
        #[values(ValenceModel::smiles(), ValenceModel::default())] valence: ValenceModel,
    ) {
        let model = ChemistryModel {
            valence,
            ..ChemistryModel::default()
        };
        let molecule = ingest_smiles_with(
            input,
            &SmilesIoConfig::opensmiles(),
            &model,
            &ResolveConfig::default(),
        )
        .unwrap();
        assert_eq!(molecule.atom(AtomId(0)).attributes, &expected);
        assert_eq!(
            molecule.stereo_atoms().ids().collect::<Vec<_>>(),
            vec![StereoAtomId(0)]
        );
        assert_eq!(
            molecule
                .stereo_atom(StereoAtomId(0))
                .ligands()
                .map(|ligand| StereoLigand::new(ligand.atom_id(), ligand.kind()))
                .collect::<Vec<_>>(),
            &[
                StereoLigand::new(AtomId(0), StereoLigandKind::LonePair),
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
            ]
        );
    }

    #[rstest]
    #[case::carbon("[C@](F)(Cl)(Br)I", atom_dsl!("C#i=#c0#h0#n0#u0#s"))]
    #[case::nitrogen_cation("[N@+](F)(Cl)(Br)I", atom_dsl!("N#i=#c+#h0#n0#u0#s"))]
    #[case::phosphorus_cation("[P@+](F)(Cl)(Br)I", atom_dsl!("P#i=#c+#h0#n0#u0#s"))]
    fn test_ingest_smiles_with_four_ligands(
        #[case] input: &str,
        #[case] expected: AtomForm,
        #[values(ValenceModel::smiles(), ValenceModel::default())] valence: ValenceModel,
    ) {
        let model = ChemistryModel {
            valence,
            ..ChemistryModel::default()
        };
        let molecule = ingest_smiles_with(
            input,
            &SmilesIoConfig::opensmiles(),
            &model,
            &ResolveConfig::default(),
        )
        .unwrap();
        assert_eq!(molecule.atom(AtomId(0)).attributes, &expected);
        assert_eq!(
            molecule.stereo_atoms().ids().collect::<Vec<_>>(),
            vec![StereoAtomId(0)]
        );
        assert_eq!(
            molecule
                .stereo_atom(StereoAtomId(0))
                .ligands()
                .map(|ligand| StereoLigand::new(ligand.atom_id(), ligand.kind()))
                .collect::<Vec<_>>(),
            &[
                StereoLigand::new(AtomId(1), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(2), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(3), StereoLigandKind::Atom),
                StereoLigand::new(AtomId(4), StereoLigandKind::Atom),
            ]
        );
    }

    #[rstest]
    #[case::carbon_anion("[C-](F)(Cl)(Br)I", Element::C, -1)]
    #[case::carbon_anion_stereo("[C@-](F)(Cl)(Br)I", Element::C, -1)]
    #[case::carbon_cation("[C+](F)(Cl)(Br)I", Element::C, 1)]
    #[case::carbon_cation_stereo("[C@+](F)(Cl)(Br)I", Element::C, 1)]
    #[case::sulfur_cation("[S@+](F)(Cl)(Br)I", Element::S, 1)]
    #[case::sulfur_anion("[S@-](F)(Cl)Br", Element::S, -1)]
    fn test_ingest_smiles_with_atom_typing_error(
        #[case] input: &str,
        #[case] element: Element,
        #[case] charge: i8,
    ) {
        assert_eq!(
            ingest_smiles_with(
                input,
                &SmilesIoConfig::opensmiles(),
                &ChemistryModel::default(),
                &ResolveConfig::default()
            ),
            Err(SmilesInputError::Contradiction(
                ResolveContradiction::Valence(ValenceContradiction::AtomTyping(
                    AtomTypingError::NoMatch {
                        atom_id: AtomId(0),
                        element,
                        charge: Some(charge),
                    }
                ),)
            )),
        );
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
                isotope: IsotopePolicy::Strict,
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
            isotope: IsotopePolicy::Strict,
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
            aromaticity: AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Clar, tie_break: AromaticityTieBreak::Strict },
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
    #[case::ring_equivalent("[CH3:1][C@H:2]1[CH2:3][CH2:4][CH2:5][CH2:6][O:7]1>>[O:7]1[CH2:6][CH2:5][CH2:4][CH2:3][C@@H:2]1[CH3:1]", true)]
    #[case::ring_mirror("[CH3:1][C@H:2]1[CH2:3][CH2:4][CH2:5][CH2:6][O:7]1>>[O:7]1[CH2:6][CH2:5][CH2:4][CH2:3][C@H:2]1[CH3:1]", false)]
    #[case::alkene_equivalent("[F:1]/[CH:2]=[CH:3]/[Cl:4]>>[Cl:4]/[CH:3]=[CH:2]/[F:1]", true)]
    #[case::alkene_opposite("[F:1]/[CH:2]=[CH:3]/[Cl:4]>>[F:1]/[CH:2]=[CH:3]\\[Cl:4]", false)]
    fn test_ingest_reaction_smiles_stereo(#[case] input: &str, #[case] expected: bool) {
        let reaction = ingest_reaction_smiles(input).unwrap();
        let span = reaction.to_reaction_span().unwrap();
        let context = CanonicalizeContext {
            para_stereo: false,
            automorphism_algorithm: umol_graph_core::AutomorphismAlgorithm::Nauty,
        };
        assert_eq!(span.lhs().canonical_eq(&span.rhs(), &context), expected);
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
    #[case::reactant_direction("F/C=C>>C", ReactionSmilesInputError::Syntax(SmilesParseError::DanglingBondDirection { bond: 0 }))]
    #[case::product_direction("C>>F/C=C", ReactionSmilesInputError::Syntax(SmilesParseError::DanglingBondDirection { bond: 0 }))]
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
            aromaticity: AromaticityModel { scope: ElementScope::Any, rule: AromaticityRule::Clar, tie_break: AromaticityTieBreak::Strict },
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
        let system = reaction.lhs().aromatic_system(AromaticSystemId(0));

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
            reaction.lhs().aromatic_systems().ids().collect::<Vec<_>>(),
            vec![AromaticSystemId(0)]
        );
        assert_eq!(reaction.deltas(), &Deltas::new());
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
                isotope: IsotopePolicy::Strict,
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
                .lhs()
                .atoms()
                .iter()
                .map(|atom| atom.attributes.constraints.aromatic_valence().cloned())
                .collect::<Vec<_>>(),
            vec![Some(AromaticValenceForm::Aromatic(NumForm::Undetermined)); 5]
        );
        assert_eq!(
            reaction
                .lhs()
                .bonds()
                .iter()
                .map(|bond| bond.attributes.constraints.aromatic())
                .collect::<Vec<_>>(),
            vec![BooleanForm::Lit(true); 5]
        );
        assert_eq!(
            reaction.lhs().aromatic_systems().ids().collect::<Vec<_>>(),
            Vec::<AromaticSystemId>::new()
        );
        assert_eq!(reaction.deltas(), &Deltas::new());
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
