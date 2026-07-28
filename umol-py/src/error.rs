//! Binding exceptions.

use pyo3::exceptions::{PyException, PyRuntimeError, PyValueError};
use pyo3::{create_exception, PyErr};
use umol_ast::ast::Contradiction as AstContradiction;
use umol_ast::dsl::{MetadataError as AstMetadataError, ParseError as AstParseError};
use umol_graph::fingerprint::FingerprintError as GraphFingerprintError;
use umol_graph::ingest::SmilesInputError as GraphSmilesInputError;

create_exception!(
    umol,
    ContradictionError,
    PyException,
    "Raised when a umol operation reaches a contradiction."
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
    "Raised when DSL metadata violates namespace or AST-coherence invariants."
);

create_exception!(
    umol,
    UnderdeterminedError,
    PyException,
    "Raised when an operation requires a determined molecular value."
);

/// Map an `umol_ast` parse error onto the catchable `umol.ParseError`.
pub(crate) fn parse_error(error: AstParseError) -> PyErr {
    ParseError::new_err(error.to_string())
}

/// Map an `umol_ast` contradiction onto the catchable `umol.ContradictionError`.
pub(crate) fn contradiction_error(error: AstContradiction) -> PyErr {
    ContradictionError::new_err(error.to_string())
}

/// Map an `umol_ast` metadata error onto the catchable `umol.MetadataError`.
pub(crate) fn metadata_error(error: AstMetadataError) -> PyErr {
    MetadataError::new_err(error.to_string())
}

/// Map the resolved SMILES operation error onto the public Python taxonomy.
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

/// Map a fingerprint operation error onto the public Python taxonomy.
pub(crate) fn fingerprint_error(error: GraphFingerprintError) -> PyErr {
    match error {
        GraphFingerprintError::NotGround => {
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

#[cfg(test)]
mod tests {
    use pyo3::prelude::*;
    use rstest::rstest;
    use umol_ast::ast::{AtomId, BondId, Entity};
    use umol_ast::dsl::MetadataError as AstMetadataError;
    use umol_graph::ingest::ingest_smiles;
    use umol_graph::ops::aromaticity::{
        AromaticityContradiction as GraphAromaticityContradiction,
        AromaticityError as GraphAromaticityError,
    };
    use umol_graph::ops::resolve::{
        ResolverContradiction as GraphResolverContradiction, ResolverError as GraphResolverError,
    };

    use super::*;

    #[rstest]
    fn test_parse_error() {
        Python::attach(|py| {
            let error = parse_error(AstParseError::ExpectedElement);
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
            let error = contradiction_error(AstContradiction);
            assert!(error.is_instance_of::<ContradictionError>(py));
            assert_eq!(
                error.value(py).str().unwrap().extract::<String>().unwrap(),
                "reached a contradiction"
            );
        });
    }

    #[rstest]
    #[case::duplicate_keyword(
        AstMetadataError::DuplicateKeyword("site".to_string()),
        "duplicate keyword: site",
    )]
    #[case::duplicate_atom_alias(
        AstMetadataError::DuplicateAtomAlias("carbon".to_string()),
        "atom DSL already has alias: carbon",
    )]
    #[case::entity_out_of_range(
        AstMetadataError::EntityOutOfRange(Entity::Atom(AtomId(2))),
        "metadata entity is out of range: atom 2"
    )]
    #[case::entity_not_added(
        AstMetadataError::EntityNotAdded(Entity::Bond(BondId(3))),
        "metadata entity is not introduced by an add delta: bond 3"
    )]
    fn test_metadata_error(#[case] input: AstMetadataError, #[case] expected_message: &str) {
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
    #[case::syntax(
        ingest_smiles(" C").unwrap_err(),
        "ParseError",
        "Leading whitespace"
    )]
    #[case::model_conversion(
        ingest_smiles("C[S@]C").unwrap_err(),
        "ModelConversionError",
        "tetrahedral stereo at atom 1 with 2 ligands, expected 3 or 4 ligands"
    )]
    #[case::contradiction(
        GraphSmilesInputError::Contradiction(GraphResolverContradiction::Aromaticity(
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
        GraphSmilesInputError::Execution(GraphResolverError::Aromaticity(
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
    #[case::not_ground(
        GraphFingerprintError::NotGround,
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
}
