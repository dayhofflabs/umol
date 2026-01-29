//! Minimal property-based crash test coverage for m4 parser

use proptest::prelude::*;
use proptest::sample::select;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_models_graph::io::smiles::config::SmilesIoConfig;
use umol_models_graph::io::smiles::parser::parse_smiles_bytes_with;
use umol_models_graph::io::smiles::{parse_smiles_bytes, ParseError};

// Generate ASCII strings from a token-friendly alphabet to bias towards SMILES-like inputs.
// This is intentionally permissive; the property is "no panics".
fn smilesish() -> impl Strategy<Value = Vec<u8>> {
    // Common SMILES characters: letters, digits, bonds, ring, parens, brackets, slash/backslash, percent, dot
    const ALPHABET: &[u8] =
        b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-=#$:/\\().%[]";
    proptest::collection::vec(select(ALPHABET.to_vec()), 0..256)
}

// TODO: Add more fine-grained property tests
// - Chains of arbitrary length are allowed (single atom, multiple atoms, different bond orders, brackets)
// - Branches of arbitrary length are allowed
// - Groups of arbitrary length are allowed
// - Rings of arbitrary length are allowed
// - Arbitrary number of components is allowed
// - Arbitrary number of branches is allowed
// - Arbitrary number of groups is allowed
// - Arbitrary number of rings is allowed
// - Arbitrary nesting of branches is allowed
// - Arbitrary nesting of rings is allowed
// - Arbitrary nesting of branches and rings is allowed
// - Arbitrary number of chiral centers is allowed
// - Arbitrary number of stereogenic bonds is allowed

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::WithSource("proptest-regressions"))),
        cases: 10000,
        ..Config::default()}
    )]

    // Crash-only: parser should never panic on arbitrary ASCII up to length 256
    #[test]
    fn never_panics_on_ascii(input in smilesish()) {
        let config = SmilesIoConfig::basic_opensmiles();
        let _ = std::panic::catch_unwind(|| {
            let _ = parse_smiles_bytes_with(&input, &config);
        }).expect("parse_smiles panicked");
    }

    // Error spans must point within the input bounds
    #[test]
    fn error_positions_within_bounds(input in smilesish()) {
        let res = parse_smiles_bytes(&input);
        if let Err(err) = res {
            let len = input.len();
            let ok = match err {
                ParseError::LeadingWhitespace => true,
                | ParseError::InvalidToken { pos } => pos < len,
                | ParseError::InvalidElement { pos } => pos < len,
                | ParseError::UnbalancedOpenParen { pos } => pos < len,
                | ParseError::UnbalancedCloseParen { pos } => pos < len,
                | ParseError::EmptyBranch { pos } => pos < len,
                | ParseError::EmptyGroup { pos } => pos < len,
                | ParseError::NonfinalGroup { pos } => pos < len,
                | ParseError::TrailingBond { pos } => pos < len,
                | ParseError::ConsecutiveBonds { pos } => pos < len,
                | ParseError::LeadingBond { pos } => pos < len,
                | ParseError::InvalidRingIndex { pos } => pos < len,
                | ParseError::LeadingRing { pos } => pos < len,
                | ParseError::LeadingDot { pos } => pos < len,
                | ParseError::TrailingDot { pos } => pos < len,
                | ParseError::ConsecutiveDots { pos } => pos < len,
                | ParseError::UnbalancedOpenBracket { pos } => pos < len,
                | ParseError::UnbalancedCloseBracket { pos } => pos < len,
                | ParseError::InvalidBracket { pos } => pos < len,
                | ParseError::MissingClassIndex { pos } => pos < len,
                | ParseError::StrayBracketField { pos } => pos < len,
                | ParseError::DuplicateBracketField { pos } => pos < len,
                | ParseError::BracketHwithHcount { pos } => pos < len,
                | ParseError::ChiralityOutOfRange { pos } => pos < len,
                | ParseError::InvalidCxProperty { pos } => pos < len,
                | ParseError::UnbalancedRingIndex { open_pos } => open_pos < len,
                | ParseError::MismatchedRingBondOrders { pos, open_pos } => pos < len && open_pos < len,
                | ParseError::MismatchedRingBondDirs { pos, open_pos } => pos < len && open_pos < len,
                | ParseError::MismatchedRingBondDonations { pos, open_pos } => pos < len && open_pos < len,
                | ParseError::DotBeforeRing { pos } => pos < len,
                | ParseError::EmptyBracket { pos } => pos < len,
                | ParseError::MissingChiralityIndex { pos } => pos < len,
                | ParseError::IndexOutOfBounds(_) => true,
            };
            prop_assert!(ok, "error positions out of bounds: {:?}, len={}", err, len);
        }
    }

    // Bonds in successful parses must reference valid atom indices
    // Note: self-loop bonds (e.g., C11) are syntactically valid and checked during topology validation
    #[test]
    fn bonds_well_formed_on_success(input in smilesish()) {
        if let Ok(mol) = parse_smiles_bytes(&input) {
            let n = mol.atoms.len() as u32;
            for b in &mol.bonds {
                let sa = b.start_atom();
                let ea = b.end_atom();
                // Bonds use 0-based atom indices
                prop_assert!(sa < n && ea < n, "bond endpoints out of bounds: {}-{} / n={}", sa, ea, n);
            }
        }
    }
}
