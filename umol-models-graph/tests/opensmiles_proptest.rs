//! Minimal property-based crash test coverage for m4 parser

use proptest::prelude::*;
use proptest::sample::select;
use umol_models_graph::io::smiles::{parse_smiles_m4, M4Error};

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
