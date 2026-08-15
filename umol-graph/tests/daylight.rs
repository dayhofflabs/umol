//! Daylight-docs local instrument (discussion/194 S6b2): resolves every
//! candidate SMILES extracted from the Daylight theory manual and tutorial
//! examples under the smiles preset and writes the classification and
//! manifest for review. Every example in those pages is canonical by
//! definition, so any non-determined row is either a documented didactic
//! divergence or a dialect gap.
//!
//! The corpus is a local instrument, never committed (`materials/` is
//! git-ignored); fetch it with `scripts/fetch-daylight-smiles.sh`. The test
//! skips when the candidate list is absent. Run with
//! `cargo test -p umol-graph --test daylight -- --ignored --nocapture`.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

use umol_graph::ingest::{ingest_smiles, SmilesInputError};
use umol_io::smiles::Smiles;

#[test]
#[ignore = "local instrument; requires materials/formats/daylight/candidates.tsv"]
fn daylight_instrument() {
    let staging = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../materials/formats/daylight");
    let Ok(candidates) = fs::read_to_string(staging.join("candidates.tsv")) else {
        eprintln!("candidates.tsv absent; run scripts/fetch-daylight-smiles.sh");
        return;
    };

    let mut class_counts: BTreeMap<&str, usize> = BTreeMap::new();
    let mut manifest = String::new();
    let mut parsed = 0usize;
    let mut seen = BTreeSet::new();

    for line in candidates.lines() {
        let Some((page, candidate)) = line.split_once('\t') else {
            continue;
        };
        if !seen.insert(candidate.to_owned()) {
            continue;
        }
        // The parser self-curates the extraction: non-SMILES prose drops out.
        if Smiles::parse(candidate).is_err() {
            continue;
        }
        parsed += 1;
        let class = match ingest_smiles(candidate) {
            Ok(_) => "determined",
            Err(SmilesInputError::Underdetermined(_)) => "underdetermined",
            Err(SmilesInputError::Contradiction(_)) => "contradictory",
            Err(_) => "execution failure",
        };
        *class_counts.entry(class).or_default() += 1;
        if class != "determined" {
            let _ = writeln!(manifest, "{class} [{page}] {candidate}");
        }
    }

    let mut distribution = String::new();
    let _ = writeln!(distribution, "parsed candidates: {parsed}");
    for (class, count) in &class_counts {
        let _ = writeln!(distribution, "{class}: {count}");
    }
    fs::write(staging.join("distribution.txt"), &distribution)
        .expect("staging directory is writable");
    fs::write(staging.join("manifest.txt"), &manifest).expect("staging directory is writable");
    println!("{distribution}");
    println!("{manifest}");
}
