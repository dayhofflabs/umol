//! Minimal property-based crash test coverage for m4 parser

use proptest::prelude::*;
use proptest::sample::select;
use umol_models_graph::io::smiles::{parse_smiles_m4, M4Error};
use umol_models_graph::io::smiles::{parse_smiles_m5, M5Error};

// Generate ASCII strings from a token-friendly alphabet to bias towards SMILES-like inputs.
// This is intentionally permissive; the property is "no panics".
fn smilesish() -> impl Strategy<Value = Vec<u8>> {
    // Common SMILES characters: letters, digits, bonds, ring, parens, brackets, slash/backslash, percent, dot
    const ALPHABET: &[u8] =
        b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-=#$:/\\().%[]";
    proptest::collection::vec(select(ALPHABET.to_vec()), 0..256)
}

proptest! {
    // Crash-only: parser should never panic on arbitrary ASCII up to length 256
    #[test]
    fn m4_never_panics_on_ascii(input in smilesish()) {
        let _ = std::panic::catch_unwind(|| {
            let _ = parse_smiles_m4(&input);
        }).expect("parse_smiles_m4 panicked");
    }
}

proptest! {
    // Crash-only for M5
    #[test]
    fn m5_never_panics_on_ascii(input in smilesish()) {
        let _ = std::panic::catch_unwind(|| {
            let _ = parse_smiles_m5(&input);
        }).expect("parse_smiles_m5 panicked");
    }
}

proptest! {
    // Error spans must point within the input bounds
    #[test]
    fn m4_error_positions_within_bounds(input in smilesish()) {
        let res = parse_smiles_m4(&input);
        if let Err(err) = res {
            let len = input.len();
            let ok = match err {
                M4Error::UnsupportedToken { pos }
                | M4Error::UnbalancedBranchOpen { pos }
                | M4Error::UnbalancedBranchClose { pos }
                | M4Error::EmptyBranch { pos }
                | M4Error::EmptyGroup { pos }
                | M4Error::TopLevelGroupTrailing { pos }
                | M4Error::TrailingBond { pos }
                | M4Error::ConsecutiveBond { pos }
                | M4Error::LeadingBond { pos }
                | M4Error::RingIndexInvalid { pos }
                | M4Error::LeadingRing { pos }
                | M4Error::RingSelfLoop { pos }
                | M4Error::RingTwoMember { pos } => pos < len,

                M4Error::RingBondDirConflict { pos, open_pos }
                | M4Error::RingBondOrderConflict { pos, open_pos } => pos < len && open_pos < len,

                M4Error::RingUnclosed { open_pos } => open_pos < len,

                M4Error::LeadingDot { pos }
                | M4Error::TrailingDot { pos }
                | M4Error::ConsecutiveDot { pos } => pos < len,
            };
            prop_assert!(ok, "error positions out of bounds: {:?}, len={}", err, len);
        }
    }
}

proptest! {
    // Error spans must point within the input bounds (M5)
    #[test]
    fn m5_error_positions_within_bounds(input in smilesish()) {
        let res = parse_smiles_m5(&input);
        if let Err(err) = res {
            let len = input.len();
            let ok = match err {
                M5Error::UnsupportedToken { pos }
                | M5Error::UnbalancedBranchOpen { pos }
                | M5Error::UnbalancedBranchClose { pos }
                | M5Error::EmptyBranch { pos }
                | M5Error::EmptyGroup { pos }
                | M5Error::TopLevelGroupTrailing { pos }
                | M5Error::TrailingBond { pos }
                | M5Error::ConsecutiveBond { pos }
                | M5Error::LeadingBond { pos }
                | M5Error::RingIndexInvalid { pos }
                | M5Error::LeadingRing { pos }
                | M5Error::RingSelfLoop { pos }
                | M5Error::RingTwoMember { pos }
                | M5Error::LeadingDot { pos }
                | M5Error::TrailingDot { pos }
                | M5Error::ConsecutiveDot { pos }
                | M5Error::UnclosedBracket { pos }
                | M5Error::InvalidBracket { pos }
                | M5Error::BracketHCountTwoDigits { pos }
                | M5Error::BracketEmptyClass { pos } => pos < len,

                M5Error::RingBondDirConflict { pos, open_pos }
                | M5Error::RingBondOrderConflict { pos, open_pos } => pos < len && open_pos < len,

                M5Error::RingUnclosed { open_pos } => open_pos < len,
            };
            prop_assert!(ok, "error positions out of bounds: {:?}, len={}", err, len);
        }
    }
}

proptest! {
    // Bonds in successful parses must reference valid, distinct atom indices
    #[test]
    fn m4_bonds_well_formed_on_success(input in smilesish()) {
        if let Ok(mol) = parse_smiles_m4(&input) {
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

proptest! {
    // Bonds well-formed for M5
    #[test]
    fn m5_bonds_well_formed_on_success(input in smilesish()) {
        if let Ok(mol) = parse_smiles_m5(&input) {
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

proptest! {
    // M5 is a superset of M4 on inputs without brackets: outcomes should match
    #[test]
    fn m5_matches_m4_without_brackets(input in smilesish()) {
        if input.iter().any(|&b| b == b'[' || b == b']') { return Ok(()); }
        let r4 = parse_smiles_m4(&input);
        let r5 = parse_smiles_m5(&input);
        match (r4, r5) {
            (Ok(m4), Ok(m5)) => prop_assert_eq!(m4, m5),
            (Err(e4), Err(e5)) => prop_assert_eq!(format!("{:?}", e4), format!("{:?}", e5)),
            (a, b) => prop_assert!(false, "mismatch: {:?} vs {:?}", a, b),
        }
    }
}
