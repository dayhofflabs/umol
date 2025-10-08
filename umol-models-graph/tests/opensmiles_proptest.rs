//! Minimal property-based crash test coverage for m4 parser

use proptest::prelude::*;
use proptest::sample::select;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_models_graph::io::smiles::config::SmilesParseFlags;
use umol_models_graph::io::smiles::parser::parse_smiles_with;
use umol_models_graph::io::smiles::{parse_smiles, ParseError};

// Generate ASCII strings from a token-friendly alphabet to bias towards SMILES-like inputs.
// This is intentionally permissive; the property is "no panics".
fn smilesish() -> impl Strategy<Value = Vec<u8>> {
    // Common SMILES characters: letters, digits, bonds, ring, parens, brackets, slash/backslash, percent, dot
    const ALPHABET: &[u8] =
        b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-=#$:/\\().%[]";
    proptest::collection::vec(select(ALPHABET.to_vec()), 0..256)
}

proptest! {
    #![proptest_config(Config {
        failure_persistence: Some(Box::new(FileFailurePersistence::WithSource("proptest-regressions"))),
        cases: 10000,
        ..Config::default()}
    )]

    // Crash-only: parser should never panic on arbitrary ASCII up to length 256
    #[test]
    fn never_panics_on_ascii_lenient(input in smilesish()) {
        let flags = SmilesParseFlags::EXTENDED_WS | SmilesParseFlags::ALLOWS_COMMENTS | SmilesParseFlags::EXPLICIT_EOI;
        let _ = std::panic::catch_unwind(|| {
            let _ = parse_smiles_with(&input, &flags);
        }).expect("parse_smiles(lenient) panicked");
    }

    // Error spans must point within the input bounds (M6 strict)
    #[test]
    fn error_positions_within_bounds(input in smilesish()) {
        let res = parse_smiles(&input);
        if let Err(err) = res {
            let len = input.len();
            let ok = match err {
                ParseError::InvalidWhitespace { pos }
                | ParseError::InvalidComment { pos }
                | ParseError::InvalidToken { pos }
                | ParseError::InvalidElement { pos }
                | ParseError::UnbalancedOpenParen { pos }
                | ParseError::UnbalancedCloseParen { pos }
                | ParseError::EmptyBranch { pos }
                | ParseError::EmptyGroup { pos }
                | ParseError::NonfinalGroup { pos }
                | ParseError::TrailingBond { pos }
                | ParseError::ConsecutiveBonds { pos }
                | ParseError::LeadingBond { pos }
                | ParseError::InvalidRingIndex { pos }
                | ParseError::LeadingRing { pos }
                | ParseError::LeadingDot { pos }
                | ParseError::TrailingDot { pos }
                | ParseError::ConsecutiveDots { pos }
                | ParseError::UnbalancedOpenBracket { pos }
                | ParseError::UnbalancedCloseBracket { pos }
                | ParseError::InvalidBracket { pos }
                | ParseError::MissingClassIndex { pos }
                | ParseError::StrayBracketField { pos }
                | ParseError::DuplicateBracketField { pos }
                | ParseError::BracketHwithHcount { pos }
                | ParseError::ChiralityOutOfRange { pos }
                | ParseError::UnterminatedBlockComment { pos }
                 => pos < len,
                | ParseError::UnbalancedRingIndex { open_pos } => open_pos < len,
                | ParseError::MismatchedRingBondDirs { pos, open_pos }
                | ParseError::MismatchedRingBondOrders { pos, open_pos } => pos < len && open_pos < len,
                | ParseError::DotBeforeRing { pos } => pos < len,
                | ParseError::EmptyBracket { pos } => pos < len,
                | ParseError::MissingChiralityIndex { pos } => pos < len,
            };
            prop_assert!(ok, "error positions out of bounds: {:?}, len={}", err, len);
        }
    }

    // Bonds in successful parses must reference valid, distinct atom indices
    #[test]
    fn bonds_well_formed_on_success(input in smilesish()) {
        if let Ok(mol) = parse_smiles(&input) {
            let n = mol.atoms.len() as u32;
            for b in &mol.bonds {
                let sa = b.start_atom;
                let ea = b.end_atom;
                prop_assert!(sa < n && ea < n, "bond endpoints out of bounds: {}-{} / n={}", sa, ea, n);
                prop_assert!(sa != ea, "self-loop bond unexpectedly present");
            }
        }
    }
}
