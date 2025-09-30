//! Minimal property-based crash test coverage for m4 parser

use proptest::prelude::*;
use proptest::sample::select;
use proptest::test_runner::{Config, FileFailurePersistence};
use umol_models_graph::io::config::SmilesParseFlags;
use umol_models_graph::io::smiles::parser::parse_smiles_inner;
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
    fn m6_never_panics_on_ascii_lenient(input in smilesish()) {
        let flags = SmilesParseFlags::INTERTOKEN_WS | SmilesParseFlags::COMMENTS | SmilesParseFlags::EXPLICIT_EOI;
        let _ = std::panic::catch_unwind(|| {
            let _ = parse_smiles_inner(&input, flags);
        }).expect("parse_smiles_m6(lenient) panicked");
    }

    // Error spans must point within the input bounds (M6 strict)
    #[test]
    fn m6_error_positions_within_bounds(input in smilesish()) {
        let res = parse_smiles(&input);
        if let Err(err) = res {
            let len = input.len();
            let ok = match err {
                ParseError::InvalidWhitespace { pos }
                | ParseError::InvalidComment { pos }
                | ParseError::UnsupportedToken { pos }
                | ParseError::UnbalancedBranchOpen { pos }
                | ParseError::UnbalancedBranchClose { pos }
                | ParseError::EmptyBranch { pos }
                | ParseError::EmptyGroup { pos }
                | ParseError::TopLevelGroupTrailing { pos }
                | ParseError::TrailingBond { pos }
                | ParseError::ConsecutiveBond { pos }
                | ParseError::LeadingBond { pos }
                | ParseError::RingIndexInvalid { pos }
                | ParseError::LeadingRing { pos }
                | ParseError::RingSelfLoop { pos }
                | ParseError::RingTwoMember { pos }
                | ParseError::RingMultipleRings { pos }
                | ParseError::LeadingDot { pos }
                | ParseError::TrailingDot { pos }
                | ParseError::ConsecutiveDot { pos }
                | ParseError::UnbalancedOpenBracket { pos }
                | ParseError::UnbalancedCloseBracket { pos }
                | ParseError::InvalidBracket { pos }
                | ParseError::BracketEmptyClass { pos }
                | ParseError::FieldOutsideBracket { pos }
                | ParseError::BracketDuplicateField { pos }
                | ParseError::BracketHOnH { pos }
                | ParseError::GroupLeadingConnector { pos }
                | ParseError::UnterminatedBlockComment { pos }
                 => pos < len,
                | ParseError::RingUnclosed { open_pos } => open_pos < len,
                ParseError::RingBondDirConflict { pos, open_pos }
                | ParseError::RingBondOrderConflict { pos, open_pos } => pos < len && open_pos < len,
            };
            prop_assert!(ok, "error positions out of bounds: {:?}, len={}", err, len);
        }
    }

    // Bonds in successful parses must reference valid, distinct atom indices
    #[test]
    fn m6_bonds_well_formed_on_success(input in smilesish()) {
        if let Ok(mol) = parse_smiles(&input) {
            let n = mol.atoms.len() as u32;
            for b in &mol.bonds {
                let sa = b.start_atom.expect("bond missing start");
                let ea = b.end_atom.expect("bond missing end");
                prop_assert!(sa < n && ea < n, "bond endpoints out of bounds: {}-{} / n={}", sa, ea, n);
                prop_assert!(sa != ea, "self-loop bond unexpectedly present");
            }
        }
    }
}
