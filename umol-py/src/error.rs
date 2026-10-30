//! Binding exceptions.

use pyo3::exceptions::{PyException, PyIndexError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::{create_exception, PyErr};
use umol_graph::export::{
    ConveyError as GraphConveyError, ReactionConveyError as GraphReactionConveyError,
    ReactionSmilesOutputError as GraphReactionSmilesOutputError,
    SmilesOutputError as GraphSmilesOutputError,
};
use umol_graph::fingerprint::FingerprintError as GraphFingerprintError;
use umol_graph::ingest::{
    MoleculeInterpretationError as GraphMoleculeInterpretationError,
    ReactionInterpretationError as GraphReactionInterpretationError,
    ReactionSmilesInputError as GraphReactionSmilesInputError,
    SmilesInputError as GraphSmilesInputError,
};
use umol_graph::ops::stereo::StereoPerceptionError as GraphStereoPerceptionError;
use umol_graph::ops::transform::{
    AddStereoError as GraphAddStereoError, AromatizeError as GraphAromatizeError,
    FoldHydrogensError as GraphFoldHydrogensError, KekulizeError as GraphKekulizeError,
    RemoveStereoError as GraphRemoveStereoError, UnfoldHydrogensError as GraphUnfoldHydrogensError,
};
use umol_graph_ir::dsl::{MetadataError as GraphIrMetadataError, ParseError as GraphIrParseError};
use umol_graph_ir::ir::{
    Contradiction as GraphIrContradiction, MoleculeApplyError as GraphIrMoleculeApplyError,
    MoleculeIntegrityError as GraphIrMoleculeIntegrityError,
    ReactionIntegrityError as GraphIrReactionIntegrityError,
    TransactionError as GraphIrTransactionError,
};

use crate::atom::AtomFieldKind;
use crate::stereo::{StereoInconsistency, StereoKind, StereoSite};

create_exception!(
    umol,
    ContradictionError,
    PyException,
    "Raised when a umol operation reaches a contradiction."
);

create_exception!(
    umol,
    NoJoinError,
    PyException,
    "Raised when the terms have no least upper bound."
);

create_exception!(
    umol,
    ParseError,
    PyException,
    "Raised when textual molecular input fails to parse."
);

create_exception!(
    umol,
    ModelConversionError,
    PyException,
    "Raised when a molecular representation cannot be converted to the requested model."
);

create_exception!(
    umol,
    InvalidStructureError,
    PyException,
    "Raised when a molecular value fails an operation's structural preconditions."
);

create_exception!(
    umol,
    MetadataError,
    PyException,
    "Raised when DSL metadata violates namespace or graph-IR coherence invariants."
);

create_exception!(
    umol,
    UnderdeterminedError,
    PyException,
    "Raised when an operation requires a determined molecular value."
);

create_exception!(
    umol,
    TransactionError,
    PyException,
    "Raised when an edit batch cannot be applied or a transaction has been aborted."
);

create_exception!(
    umol,
    ConsumedError,
    PyRuntimeError,
    "Raised when an operation accesses an object whose contents have been consumed."
);

create_exception!(
    umol,
    InvalidatedViewError,
    PyRuntimeError,
    "Raised when a view or iterator is invalidated by owner mutation or consumption."
);

/// Map an `umol_graph_ir` parse error onto the catchable `umol.ParseError`.
pub(crate) fn parse_error(error: GraphIrParseError) -> PyErr {
    ParseError::new_err(error.to_string())
}

/// Map an `umol_graph_ir` contradiction onto the catchable `umol.ContradictionError`.
pub(crate) fn contradiction_error(error: GraphIrContradiction) -> PyErr {
    ContradictionError::new_err(error.to_string())
}

/// Map an `umol_graph_ir` metadata error onto the catchable `umol.MetadataError`.
pub(crate) fn metadata_error(error: GraphIrMetadataError) -> PyErr {
    MetadataError::new_err(error.to_string())
}

/// Map an `umol_graph_ir` transaction failure onto the catchable `umol.TransactionError`.
pub(crate) fn transaction_error(error: GraphIrTransactionError) -> PyErr {
    TransactionError::new_err(error.to_string())
}

/// Map an `umol_graph_ir` molecule-integrity failure onto the public Python taxonomy.
pub(crate) fn molecule_integrity_error(error: GraphIrMoleculeIntegrityError) -> PyErr {
    InvalidStructureError::new_err(error.to_string())
}

/// Map an independently assembled reaction-integrity failure onto `ValueError`.
pub(crate) fn reaction_integrity_error(error: GraphIrReactionIntegrityError) -> PyErr {
    PyValueError::new_err(error.to_string())
}

/// Map molecule-application failures onto the existing Python exception taxonomy.
pub(crate) fn molecule_apply_error(error: GraphIrMoleculeApplyError) -> PyErr {
    match error {
        GraphIrMoleculeApplyError::Transaction(error) => transaction_error(error),
        GraphIrMoleculeApplyError::Integrity(error) => molecule_integrity_error(error),
    }
}

/// Map a SMILES ingestion error onto the public Python taxonomy.
pub(crate) fn smiles_input_error(error: GraphSmilesInputError) -> PyErr {
    match error {
        GraphSmilesInputError::Syntax(error) => ParseError::new_err(error.to_string()),
        GraphSmilesInputError::ModelConversion(error) => {
            ModelConversionError::new_err(error.to_string())
        }
        GraphSmilesInputError::Contradiction(error) => {
            ContradictionError::new_err(error.to_string())
        }
        GraphSmilesInputError::Underdetermined(error) => {
            UnderdeterminedError::new_err(error.to_string())
        }
        GraphSmilesInputError::Execution(error) => PyRuntimeError::new_err(error.to_string()),
    }
}

/// Map the resolved reaction-SMILES operation error onto the public Python taxonomy.
pub(crate) fn reaction_smiles_input_error(error: GraphReactionSmilesInputError) -> PyErr {
    match error {
        GraphReactionSmilesInputError::Syntax(error) => ParseError::new_err(error.to_string()),
        GraphReactionSmilesInputError::Interpretation(error) => {
            let message = error.to_string();
            match error {
                GraphReactionInterpretationError::Reactants(
                    GraphMoleculeInterpretationError::ModelConversion(_),
                )
                | GraphReactionInterpretationError::Products(
                    GraphMoleculeInterpretationError::ModelConversion(_),
                )
                | GraphReactionInterpretationError::AmbiguousAtomMapClass { .. }
                | GraphReactionInterpretationError::AgentsUnsupported => {
                    ModelConversionError::new_err(message)
                }
                GraphReactionInterpretationError::Reactants(
                    GraphMoleculeInterpretationError::Contradiction(_),
                )
                | GraphReactionInterpretationError::Products(
                    GraphMoleculeInterpretationError::Contradiction(_),
                ) => ContradictionError::new_err(message),
                GraphReactionInterpretationError::Reactants(
                    GraphMoleculeInterpretationError::Underdetermined(_),
                )
                | GraphReactionInterpretationError::Products(
                    GraphMoleculeInterpretationError::Underdetermined(_),
                ) => UnderdeterminedError::new_err(message),
                GraphReactionInterpretationError::Reactants(
                    GraphMoleculeInterpretationError::Execution(_),
                )
                | GraphReactionInterpretationError::Products(
                    GraphMoleculeInterpretationError::Execution(_),
                ) => PyRuntimeError::new_err(message),
            }
        }
    }
}

/// Select the Python exception for a molecule that cannot be conveyed to SMILES.
fn convey_error(error: &GraphConveyError, message: String) -> PyErr {
    match error {
        GraphConveyError::Contradiction(_) => ContradictionError::new_err(message),
        GraphConveyError::Underdetermined => UnderdeterminedError::new_err(message),
        GraphConveyError::Projection(_)
        | GraphConveyError::Value { .. }
        | GraphConveyError::Constraint(_)
        | GraphConveyError::Entity { .. }
        | GraphConveyError::StereoAtom { .. }
        | GraphConveyError::StereoBond { .. } => ModelConversionError::new_err(message),
    }
}

/// Map a SMILES export error onto the public Python taxonomy.
pub(crate) fn smiles_output_error(error: GraphSmilesOutputError) -> PyErr {
    let message = error.to_string();
    match &error {
        GraphSmilesOutputError::Convey(error) => convey_error(error, message),
        GraphSmilesOutputError::Render(_) => ModelConversionError::new_err(message),
    }
}

/// Map a reaction-SMILES export error onto the public Python taxonomy.
pub(crate) fn reaction_smiles_output_error(error: GraphReactionSmilesOutputError) -> PyErr {
    let message = error.to_string();
    match &error {
        GraphReactionSmilesOutputError::Convey(GraphReactionConveyError::Materialization(_)) => {
            ContradictionError::new_err(message)
        }
        GraphReactionSmilesOutputError::Convey(
            GraphReactionConveyError::Reactants(error) | GraphReactionConveyError::Products(error),
        ) => convey_error(error, message),
        GraphReactionSmilesOutputError::Render(_) => ModelConversionError::new_err(message),
    }
}

/// Map a fingerprint operation error onto the public Python taxonomy.
pub(crate) fn fingerprint_error(error: GraphFingerprintError) -> PyErr {
    match error {
        GraphFingerprintError::NotConcrete => {
            UnderdeterminedError::new_err("fingerprint requires a determined molecule")
        }
        GraphFingerprintError::Inconsistent => {
            ContradictionError::new_err("reaction fingerprint input is inconsistent")
        }
        GraphFingerprintError::ZeroWidth => PyValueError::new_err("width must be positive"),
        GraphFingerprintError::WidthMismatch { left, right } => {
            PyValueError::new_err(format!("fingerprint width mismatch: {left} != {right}"))
        }
    }
}

/// Map aromatization errors onto the public Python exception taxonomy.
pub(crate) fn aromatize_error(error: GraphAromatizeError) -> PyErr {
    let message = error.to_string();
    match error {
        GraphAromatizeError::Setup(_) => PyRuntimeError::new_err(message),
        GraphAromatizeError::Underdetermined => UnderdeterminedError::new_err(message),
        GraphAromatizeError::Contradiction(_)
        | GraphAromatizeError::ContradictedAssertion { .. }
        | GraphAromatizeError::ValenceInvariant { .. } => ContradictionError::new_err(message),
    }
}

/// Map kekulization errors onto the public Python exception taxonomy.
pub(crate) fn kekulize_error(error: GraphKekulizeError) -> PyErr {
    let message = error.to_string();
    match error {
        GraphKekulizeError::UndeterminedElectrons(_)
        | GraphKekulizeError::UndeterminedCharge(_)
        | GraphKekulizeError::UndeterminedUnpairedElectrons(_)
        | GraphKekulizeError::UndeterminedExposedAtomCharge { .. }
        | GraphKekulizeError::UndeterminedExposedAtomLonePairs { .. } => {
            UnderdeterminedError::new_err(message)
        }
        GraphKekulizeError::InvalidNodeOrder { .. } => PyValueError::new_err(message),
        GraphKekulizeError::ElectronCountMismatch { .. } => InvalidStructureError::new_err(message),
        GraphKekulizeError::ContradictedAssertion { .. }
        | GraphKekulizeError::ValenceInvariant { .. } => ContradictionError::new_err(message),
        GraphKekulizeError::UnsupportedElectronContribution { .. }
        | GraphKekulizeError::OpenShell(_)
        | GraphKekulizeError::UnsupportedCharge { .. }
        | GraphKekulizeError::MixedExposureDemand(_)
        | GraphKekulizeError::MatchingDeficiency { .. }
        | GraphKekulizeError::NoMatching(_)
        | GraphKekulizeError::NonBipartiteMatching(_) => ModelConversionError::new_err(message),
    }
}

/// Map hydrogen-folding errors onto the public Python exception taxonomy.
pub(crate) fn fold_hydrogens_error(error: GraphFoldHydrogensError) -> PyErr {
    let message = error.to_string();
    match error {
        GraphFoldHydrogensError::NegativeImplicitHydrogenCount { .. } => {
            InvalidStructureError::new_err(message)
        }
        GraphFoldHydrogensError::StereoHydrogenConflict { .. } => {
            ModelConversionError::new_err(message)
        }
    }
}

/// Map hydrogen-unfolding errors onto the public Python exception taxonomy.
pub(crate) fn unfold_hydrogens_error(error: GraphUnfoldHydrogensError) -> PyErr {
    let message = error.to_string();
    match error {
        GraphUnfoldHydrogensError::NonConcreteStereoHydrogenCount { .. } => {
            UnderdeterminedError::new_err(message)
        }
        GraphUnfoldHydrogensError::NegativeImplicitHydrogenCount { .. }
        | GraphUnfoldHydrogensError::StereoHydrogenCountMismatch { .. } => {
            InvalidStructureError::new_err(message)
        }
    }
}

/// Preserve stereo query diagnostics in the existing exception categories.
pub(crate) fn stereo_perception_error(py: Python<'_>, error: GraphStereoPerceptionError) -> PyErr {
    let message = error.to_string();
    let exception = match &error {
        GraphStereoPerceptionError::InvalidSite { .. } => PyIndexError::new_err(message),
        GraphStereoPerceptionError::MissingStereoEntity { .. } => {
            InvalidStructureError::new_err(message)
        }
        GraphStereoPerceptionError::UnsupportedKind { .. }
        | GraphStereoPerceptionError::UnsupportedConfiguration { .. }
        | GraphStereoPerceptionError::UnsupportedFluxionality { .. } => {
            ModelConversionError::new_err(message)
        }
        GraphStereoPerceptionError::UnderdeterminedAtomField { .. }
        | GraphStereoPerceptionError::UnderdeterminedBondOrder { .. } => {
            UnderdeterminedError::new_err(message)
        }
        GraphStereoPerceptionError::Inconsistency(_) => ContradictionError::new_err(message),
    };
    let attach = || -> PyResult<()> {
        let value = exception.value(py);
        match error {
            GraphStereoPerceptionError::InvalidSite { site }
            | GraphStereoPerceptionError::MissingStereoEntity { site }
            | GraphStereoPerceptionError::UnsupportedConfiguration { site } => {
                value.setattr("site", StereoSite::from_rust(site))?;
            }
            GraphStereoPerceptionError::UnsupportedKind { kind }
            | GraphStereoPerceptionError::UnsupportedFluxionality { kind } => {
                value.setattr("kind", StereoKind::from_rust(kind))?;
            }
            GraphStereoPerceptionError::UnderdeterminedAtomField { atom, field } => {
                value.setattr("atom", atom.0)?;
                value.setattr("field", AtomFieldKind::from_rust(field))?;
            }
            GraphStereoPerceptionError::UnderdeterminedBondOrder { bond } => {
                value.setattr("bond", bond.0)?;
            }
            GraphStereoPerceptionError::Inconsistency(inconsistency) => {
                value.setattr(
                    "inconsistency",
                    StereoInconsistency::from_rust(inconsistency),
                )?;
            }
        }
        Ok(())
    };
    match attach() {
        Ok(()) => exception,
        Err(error) => error,
    }
}

/// Preserve addition's shared perception or inconsistency diagnostic.
pub(crate) fn add_stereo_error(py: Python<'_>, error: GraphAddStereoError) -> PyErr {
    match error {
        GraphAddStereoError::Perception(error) => stereo_perception_error(py, error),
        GraphAddStereoError::Inconsistency(error) => {
            stereo_perception_error(py, GraphStereoPerceptionError::Inconsistency(error))
        }
    }
}

/// Preserve removal's shared perception diagnostic.
pub(crate) fn remove_stereo_error(py: Python<'_>, error: GraphRemoveStereoError) -> PyErr {
    match error {
        GraphRemoveStereoError::Perception(error) => stereo_perception_error(py, error),
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use umol_graph::ingest::ingest_smiles;
    use umol_graph::ops::aromaticity::{
        AromaticityContradiction as GraphAromaticityContradiction,
        AromaticityError as GraphAromaticityError,
    };
    use umol_graph::ops::invariant::ValenceMismatch;
    use umol_graph::ops::resolve::{
        ResolveContradiction as GraphResolveContradiction, ResolveError as GraphResolveError,
        ResolveUnderdetermined as GraphResolveUnderdetermined,
    };
    use umol_graph::ops::stereo::{
        StereoInconsistency as GraphStereoInconsistency, StereoSite as GraphStereoSite,
    };
    use umol_graph_ir::dsl::MetadataError as GraphIrMetadataError;
    use umol_graph_ir::ir::{
        AromaticSystemId, AtomConstraintForm, AtomFieldKind as GraphIrAtomFieldKind, AtomId,
        BondId, Entity, EntityKind, NumForm, StereoAtomId, StereoBondId,
        StereoKind as GraphIrStereoKind,
    };
    use umol_io::smiles::ParseError as SmilesParseError;
    use umol_io::table_ir::raise::RaiseError;

    use super::*;

    #[rstest]
    fn test_parse_error() {
        Python::attach(|py| {
            let error = parse_error(GraphIrParseError::ExpectedElement);
            assert!(error.is_instance_of::<ParseError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "expected atom element"
            );
        });
    }

    #[rstest]
    fn test_contradiction_error() {
        Python::attach(|py| {
            let error = contradiction_error(GraphIrContradiction);
            assert!(error.is_instance_of::<ContradictionError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "reached a contradiction"
            );
        });
    }

    #[rstest]
    #[case::duplicate_keyword(
        GraphIrMetadataError::DuplicateKeyword("site".to_string()),
        "duplicate keyword: site",
    )]
    #[case::duplicate_atom_alias(
        GraphIrMetadataError::DuplicateAtomAlias("carbon".to_string()),
        "atom DSL already has alias: carbon",
    )]
    #[case::entity_out_of_range(
        GraphIrMetadataError::EntityOutOfRange(Entity::Atom(AtomId(2))),
        "metadata entity is out of range: atom 2"
    )]
    #[case::entity_not_added(
        GraphIrMetadataError::EntityNotAdded(Entity::Bond(BondId(3))),
        "metadata entity is not introduced by an add delta: bond 3"
    )]
    fn test_metadata_error(#[case] input: GraphIrMetadataError, #[case] expected_message: &str) {
        Python::attach(|py| {
            let error = metadata_error(input);
            assert!(error.is_instance_of::<MetadataError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected_message
            );
        });
    }

    #[rstest]
    #[case::handle_out_of_range(
        GraphIrTransactionError::HandleOutOfRange {
            kind: EntityKind::Atom,
            index: 3,
            count: 2,
        },
        "atom handle 3 is out of range for 2 entries",
    )]
    #[case::handle_removed(
        GraphIrTransactionError::HandleRemoved {
            kind: EntityKind::Bond,
            index: 1,
        },
        "bond handle 1 refers to a removed entity",
    )]
    #[case::duplicate_removal(
        GraphIrTransactionError::DuplicateRemoval {
            kind: EntityKind::StereoAtom,
        },
        "duplicate stereo atom in removal batch",
    )]
    #[case::old_state_mismatch(
        GraphIrTransactionError::OldStateMismatch,
        "precondition failed: old state does not match current"
    )]
    #[case::missing_entry(
        GraphIrTransactionError::MissingEntry,
        "missing constraint entry on remove"
    )]
    #[case::malformed_edit(
        GraphIrTransactionError::MalformedEdit(
            "AddDativeBond requires at least one participant atom",
        ),
        "malformed edit: AddDativeBond requires at least one participant atom"
    )]
    #[case::aborted(GraphIrTransactionError::Aborted, "transaction has been aborted")]
    fn test_transaction_error(
        #[case] input: GraphIrTransactionError,
        #[case] expected_message: &str,
    ) {
        Python::attach(|py| {
            let error = transaction_error(input);
            assert!(error.is_instance_of::<TransactionError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected_message
            );
        });
    }

    #[rstest]
    #[case::syntax(
        ingest_smiles(" C").unwrap_err(),
        "ParseError",
        "Leading whitespace"
    )]
    #[case::model_conversion(
        GraphSmilesInputError::ModelConversion(
            RaiseError::TetrahedralLigandCount { atom: 1, count: 2 },
        ),
        "ModelConversionError",
        "tetrahedral stereo at atom 1 with 2 ligands, expected 3 or 4 ligands"
    )]
    #[case::incomplete_stereo(
        ingest_smiles("C[S@]C").unwrap_err(),
        "ContradictionError",
        "stereo inconsistency: stereo atom StereoAtomId(0) cannot be realized"
    )]
    #[case::contradiction(
        GraphSmilesInputError::Contradiction(GraphResolveContradiction::Aromaticity(
            GraphAromaticityContradiction::HmoInvalidInput(String::from("invalid input")),
        )),
        "ContradictionError",
        "hmo: invalid input: invalid input"
    )]
    #[case::underdetermined(
        ingest_smiles("*").unwrap_err(),
        "UnderdeterminedError",
        "resolution underdetermined"
    )]
    #[case::execution(
        GraphSmilesInputError::Execution(GraphResolveError::Aromaticity(
            GraphAromaticityError::HmoMissingParameters(String::from("carbon")),
        )),
        "RuntimeError",
        "hmo: missing parameters: carbon"
    )]
    fn test_smiles_input_error(
        #[case] input: GraphSmilesInputError,
        #[case] expected_type: &str,
        #[case] expected_message: &str,
    ) {
        Python::attach(|py| {
            let error = smiles_input_error(input);
            assert_eq!(error.get_type(py).name().unwrap(), expected_type);
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected_message
            );
        });
    }

    #[rstest]
    #[case::syntax(
        GraphReactionSmilesInputError::Syntax(SmilesParseError::LeadingWhitespace),
        "ParseError",
        "Leading whitespace"
    )]
    #[case::reactant_model_conversion(
        GraphReactionSmilesInputError::Interpretation(
            GraphReactionInterpretationError::Reactants(
                GraphMoleculeInterpretationError::ModelConversion(
                    RaiseError::TetrahedralLigandCount { atom: 1, count: 2 },
                ),
            ),
        ),
        "ModelConversionError",
        "reactants: tetrahedral stereo at atom 1 with 2 ligands, expected 3 or 4 ligands"
    )]
    #[case::product_model_conversion(
        GraphReactionSmilesInputError::Interpretation(
            GraphReactionInterpretationError::Products(
                GraphMoleculeInterpretationError::ModelConversion(
                    RaiseError::TetrahedralLigandCount { atom: 2, count: 3 },
                ),
            ),
        ),
        "ModelConversionError",
        "products: tetrahedral stereo at atom 2 with 3 ligands, expected 3 or 4 ligands"
    )]
    #[case::reactant_contradiction(
        GraphReactionSmilesInputError::Interpretation(
            GraphReactionInterpretationError::Reactants(
                GraphMoleculeInterpretationError::Contradiction(
                    GraphResolveContradiction::Aromaticity(
                        GraphAromaticityContradiction::HmoInvalidInput(String::from(
                            "invalid reactant",
                        )),
                    ),
                ),
            ),
        ),
        "ContradictionError",
        "reactants: hmo: invalid input: invalid reactant"
    )]
    #[case::product_contradiction(
        GraphReactionSmilesInputError::Interpretation(GraphReactionInterpretationError::Products(
            GraphMoleculeInterpretationError::Contradiction(
                GraphResolveContradiction::Aromaticity(
                    GraphAromaticityContradiction::HmoInvalidInput(String::from(
                        "invalid product",
                    )),
                ),
            ),
        ),),
        "ContradictionError",
        "products: hmo: invalid input: invalid product"
    )]
    #[case::reactant_underdetermined(
        GraphReactionSmilesInputError::Interpretation(
            GraphReactionInterpretationError::Reactants(
                GraphMoleculeInterpretationError::Underdetermined(GraphResolveUnderdetermined),
            ),
        ),
        "UnderdeterminedError",
        "reactants: resolution underdetermined"
    )]
    #[case::product_underdetermined(
        GraphReactionSmilesInputError::Interpretation(GraphReactionInterpretationError::Products(
            GraphMoleculeInterpretationError::Underdetermined(GraphResolveUnderdetermined),
        ),),
        "UnderdeterminedError",
        "products: resolution underdetermined"
    )]
    #[case::reactant_execution(
        GraphReactionSmilesInputError::Interpretation(
            GraphReactionInterpretationError::Reactants(
                GraphMoleculeInterpretationError::Execution(GraphResolveError::Aromaticity(
                    GraphAromaticityError::HmoMissingParameters(String::from("reactant atom",)),
                ),),
            ),
        ),
        "RuntimeError",
        "reactants: hmo: missing parameters: reactant atom"
    )]
    #[case::product_execution(
        GraphReactionSmilesInputError::Interpretation(GraphReactionInterpretationError::Products(
            GraphMoleculeInterpretationError::Execution(GraphResolveError::Aromaticity(
                GraphAromaticityError::HmoMissingParameters(String::from("product atom",)),
            ),),
        ),),
        "RuntimeError",
        "products: hmo: missing parameters: product atom"
    )]
    #[case::ambiguous_atom_map_class(
        GraphReactionSmilesInputError::Interpretation(
            GraphReactionInterpretationError::AmbiguousAtomMapClass {
                class: 7,
                reactant_count: 2,
                product_count: 1,
            },
        ),
        "ModelConversionError",
        "atom-map class 7 cannot be projected into one correspondence (reactant atoms: 2, product atoms: 1)"
    )]
    #[case::agents(
        GraphReactionSmilesInputError::Interpretation(
            GraphReactionInterpretationError::AgentsUnsupported,
        ),
        "ModelConversionError",
        "reaction agents cannot be represented in Reaction"
    )]
    fn test_reaction_smiles_input_error(
        #[case] input: GraphReactionSmilesInputError,
        #[case] expected_type: &str,
        #[case] expected_message: &str,
    ) {
        Python::attach(|py| {
            let error = reaction_smiles_input_error(input);
            assert_eq!(error.get_type(py).name().unwrap(), expected_type);
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected_message
            );
        });
    }

    #[rstest]
    #[case::not_ground(
        GraphFingerprintError::NotConcrete,
        "UnderdeterminedError",
        "fingerprint requires a determined molecule"
    )]
    #[case::inconsistent(
        GraphFingerprintError::Inconsistent,
        "ContradictionError",
        "reaction fingerprint input is inconsistent"
    )]
    #[case::zero_width(
        GraphFingerprintError::ZeroWidth,
        "ValueError",
        "width must be positive"
    )]
    #[case::width_mismatch(
        GraphFingerprintError::WidthMismatch { left: 64, right: 32 },
        "ValueError",
        "fingerprint width mismatch: 64 != 32"
    )]
    fn test_fingerprint_error(
        #[case] input: GraphFingerprintError,
        #[case] expected_type: &str,
        #[case] expected_message: &str,
    ) {
        Python::attach(|py| {
            let error = fingerprint_error(input);
            assert_eq!(error.get_type(py).name().unwrap(), expected_type);
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected_message
            );
        });
    }

    #[rstest]
    #[case::setup(
        GraphAromatizeError::Setup(GraphAromaticityError::HmoMissingParameters("carbon".into())),
        "RuntimeError",
    )]
    #[case::underdetermined(GraphAromatizeError::Underdetermined, "UnderdeterminedError")]
    #[case::contradiction(
        GraphAromatizeError::Contradiction(GraphAromaticityContradiction::ClarNonBenzenoid(
            "heterocycle".into(),
        )),
        "ContradictionError",
    )]
    #[case::assertion(
        GraphAromatizeError::ContradictedAssertion {
            aromatic_system: AromaticSystemId(0),
            atom: AtomId(2),
            constraint: AtomConstraintForm::valence(3),
        },
        "ContradictionError",
    )]
    #[case::valence(
        GraphAromatizeError::ValenceInvariant {
            aromatic_system: AromaticSystemId(0),
            mismatch: ValenceMismatch::OrbitalCount {
                atom_id: AtomId(2), orbital_count: 10, electron_count: 8,
            },
        },
        "ContradictionError",
    )]
    fn test_aromatize_error(#[case] input: GraphAromatizeError, #[case] expected_type: &str) {
        let expected_message = input.to_string();
        Python::attach(|py| {
            let error = aromatize_error(input);
            assert_eq!(error.get_type(py).name().unwrap(), expected_type);
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected_message
            );
        });
    }

    #[rstest]
    #[case::electrons(
        GraphKekulizeError::UndeterminedElectrons(AromaticSystemId(0)),
        "UnderdeterminedError"
    )]
    #[case::charge(
        GraphKekulizeError::UndeterminedCharge(AromaticSystemId(0)),
        "UnderdeterminedError"
    )]
    #[case::unpaired(
        GraphKekulizeError::UndeterminedUnpairedElectrons(AromaticSystemId(0)),
        "UnderdeterminedError"
    )]
    #[case::electron_count(
        GraphKekulizeError::ElectronCountMismatch {
            aromatic_system: AromaticSystemId(0), member_count: 6, electron_count: 5,
        },
        "InvalidStructureError",
    )]
    #[case::contribution(
        GraphKekulizeError::UnsupportedElectronContribution {
            aromatic_system: AromaticSystemId(0), atom: AtomId(1), contribution: 3,
        },
        "ModelConversionError",
    )]
    #[case::open_shell(
        GraphKekulizeError::OpenShell(AromaticSystemId(0)),
        "ModelConversionError"
    )]
    #[case::unsupported_charge(
        GraphKekulizeError::UnsupportedCharge { aromatic_system: AromaticSystemId(0), charge: 2 },
        "ModelConversionError",
    )]
    #[case::mixed_exposure(
        GraphKekulizeError::MixedExposureDemand(AromaticSystemId(0)),
        "ModelConversionError"
    )]
    #[case::order(
        GraphKekulizeError::InvalidNodeOrder {
            aromatic_system: AromaticSystemId(0), missing: vec![AtomId(5)], duplicates: vec![],
        },
        "ValueError",
    )]
    #[case::deficiency(
        GraphKekulizeError::MatchingDeficiency {
            aromatic_system: AromaticSystemId(0), expected: 1, actual: 3,
        },
        "ModelConversionError",
    )]
    #[case::exposed_charge(
        GraphKekulizeError::UndeterminedExposedAtomCharge {
            aromatic_system: AromaticSystemId(0), atom: AtomId(1),
        },
        "UnderdeterminedError",
    )]
    #[case::exposed_lone_pairs(
        GraphKekulizeError::UndeterminedExposedAtomLonePairs {
            aromatic_system: AromaticSystemId(0), atom: AtomId(1),
        },
        "UnderdeterminedError",
    )]
    #[case::assertion(
        GraphKekulizeError::ContradictedAssertion {
            aromatic_system: AromaticSystemId(0),
            atom: AtomId(2),
            constraint: AtomConstraintForm::valence(2),
        },
        "ContradictionError",
    )]
    #[case::valence(
        GraphKekulizeError::ValenceInvariant {
            aromatic_system: AromaticSystemId(0),
            mismatch: ValenceMismatch::OrbitalCount {
                atom_id: AtomId(2), orbital_count: 10, electron_count: 8,
            },
        },
        "ContradictionError",
    )]
    #[case::no_matching(
        GraphKekulizeError::NoMatching(AromaticSystemId(0)),
        "ModelConversionError"
    )]
    #[case::non_bipartite(
        GraphKekulizeError::NonBipartiteMatching(AromaticSystemId(0)),
        "ModelConversionError"
    )]
    fn test_kekulize_error(#[case] input: GraphKekulizeError, #[case] expected_type: &str) {
        let expected_message = input.to_string();
        Python::attach(|py| {
            let error = kekulize_error(input);
            assert_eq!(error.get_type(py).name().unwrap(), expected_type);
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected_message
            );
        });
    }

    #[rstest]
    #[case::negative(
        GraphFoldHydrogensError::NegativeImplicitHydrogenCount { atom: AtomId(0), count: -1 },
        "InvalidStructureError",
    )]
    #[case::stereo(
        GraphFoldHydrogensError::StereoHydrogenConflict {
            entity: Entity::StereoAtom(StereoAtomId(0)), hydrogen: AtomId(1),
        },
        "ModelConversionError",
    )]
    fn test_fold_hydrogens_error(
        #[case] input: GraphFoldHydrogensError,
        #[case] expected_type: &str,
    ) {
        let expected_message = input.to_string();
        Python::attach(|py| {
            let error = fold_hydrogens_error(input);
            assert_eq!(error.get_type(py).name().unwrap(), expected_type);
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected_message
            );
        });
    }

    #[rstest]
    #[case::negative(
        GraphUnfoldHydrogensError::NegativeImplicitHydrogenCount { atom: AtomId(0), count: -1 },
        "InvalidStructureError",
    )]
    #[case::underdetermined(
        GraphUnfoldHydrogensError::NonConcreteStereoHydrogenCount {
            entity: Entity::StereoBond(StereoBondId(0)),
            atom: AtomId(1),
            value: NumForm::Undetermined,
        },
        "UnderdeterminedError",
    )]
    #[case::stereo(
        GraphUnfoldHydrogensError::StereoHydrogenCountMismatch {
            entity: Entity::StereoAtom(StereoAtomId(0)),
            atom: AtomId(1),
            implicit_hydrogens: 0,
            virtual_hydrogens: 1,
        },
        "InvalidStructureError",
    )]
    fn test_unfold_hydrogens_error(
        #[case] input: GraphUnfoldHydrogensError,
        #[case] expected_type: &str,
    ) {
        let expected_message = input.to_string();
        Python::attach(|py| {
            let error = unfold_hydrogens_error(input);
            assert_eq!(error.get_type(py).name().unwrap(), expected_type);
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                expected_message
            );
        });
    }

    #[rstest]
    #[case::invalid_atom(GraphStereoPerceptionError::InvalidSite { site: GraphStereoSite::Atom(AtomId(2)) }, "IndexError")]
    #[case::invalid_bond(GraphStereoPerceptionError::InvalidSite { site: GraphStereoSite::Bond(BondId(3)) }, "IndexError")]
    #[case::missing(GraphStereoPerceptionError::MissingStereoEntity { site: GraphStereoSite::Atom(AtomId(1)) }, "InvalidStructureError")]
    #[case::kind(GraphStereoPerceptionError::UnsupportedKind { kind: GraphIrStereoKind::SquarePlanar }, "ModelConversionError")]
    #[case::configuration(GraphStereoPerceptionError::UnsupportedConfiguration { site: GraphStereoSite::Bond(BondId(1)) }, "ModelConversionError")]
    #[case::fluxionality(GraphStereoPerceptionError::UnsupportedFluxionality { kind: GraphIrStereoKind::Tetrahedral }, "ModelConversionError")]
    #[case::field(GraphStereoPerceptionError::UnderdeterminedAtomField { atom: AtomId(2), field: GraphIrAtomFieldKind::Charge }, "UnderdeterminedError")]
    #[case::order(GraphStereoPerceptionError::UnderdeterminedBondOrder { bond: BondId(2) }, "UnderdeterminedError")]
    #[case::inconsistency(GraphStereoPerceptionError::Inconsistency(GraphStereoInconsistency::TetrahedralStereoMismatch { atom: AtomId(1), stereo_atom: StereoAtomId(0) }), "ContradictionError")]
    fn test_stereo_perception_error(
        #[case] input: GraphStereoPerceptionError,
        #[case] exception: &str,
    ) {
        Python::attach(|py| {
            let error = stereo_perception_error(py, input.clone());
            assert_eq!(error.get_type(py).name().unwrap(), exception);
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                input.to_string()
            );
            let value = error.value(py);
            match input {
                GraphStereoPerceptionError::InvalidSite { site }
                | GraphStereoPerceptionError::MissingStereoEntity { site }
                | GraphStereoPerceptionError::UnsupportedConfiguration { site } => {
                    assert_eq!(
                        value
                            .getattr("site")
                            .unwrap()
                            .extract::<StereoSite>()
                            .unwrap()
                            .to_rust(),
                        site
                    );
                }
                GraphStereoPerceptionError::UnsupportedKind { kind }
                | GraphStereoPerceptionError::UnsupportedFluxionality { kind } => {
                    assert_eq!(
                        value
                            .getattr("kind")
                            .unwrap()
                            .extract::<StereoKind>()
                            .unwrap()
                            .to_rust(),
                        kind
                    );
                }
                GraphStereoPerceptionError::UnderdeterminedAtomField { atom, field } => {
                    assert_eq!(
                        value.getattr("atom").unwrap().extract::<u32>().unwrap(),
                        atom.0
                    );
                    assert_eq!(
                        value
                            .getattr("field")
                            .unwrap()
                            .extract::<AtomFieldKind>()
                            .unwrap()
                            .to_rust(),
                        field
                    );
                }
                GraphStereoPerceptionError::UnderdeterminedBondOrder { bond } => {
                    assert_eq!(
                        value.getattr("bond").unwrap().extract::<u32>().unwrap(),
                        bond.0
                    );
                }
                GraphStereoPerceptionError::Inconsistency(inconsistency) => {
                    assert_eq!(
                        value
                            .getattr("inconsistency")
                            .unwrap()
                            .extract::<StereoInconsistency>()
                            .unwrap()
                            .to_rust(),
                        inconsistency
                    );
                }
            }
        });
    }
}
