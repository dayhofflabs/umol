//! VEHICLe local instrument (discussion/196 A5b): classifies every corpus row
//! by its resolution outcome under both valence tie-breaks, asserts the
//! aromatic-system invariants on every resolved molecule, checks
//! tautomer-cluster class consistency, and writes the class distribution and
//! failure manifest into the staging directory for review.
//!
//! The corpus is a local development instrument, never committed
//! (`materials/` is git-ignored); fetch it with `scripts/fetch-vehicle.sh`.
//! The test skips when the CSV is absent. Run with
//! `cargo test -p umol-graph --release --test vehicle -- --ignored`.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

use umol_graph::ingest::{ingest_smiles_with, SmilesInputError};
use umol_graph::ops::aromaticity::{AromaticityConfig, AromaticityPerception};
use umol_graph::ops::model::{AromaticityModel, ChemistryModel, ValenceModel, ValenceTieBreak};
use umol_graph::ops::resolve::ResolveConfig;
use umol_graph_core::AutomorphismAlgorithm;
use umol_graph_ir::ir::{CanonicalizationContext, Canonicalize, Molecule};
use umol_io::smiles::SmilesIoConfig;
use umol_utils::solution::Solution;

/// Outcome-pair class per corpus row (discussion/196, "What it measures").
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum RowClass {
    /// Both tie-breaks determined, equal: one resolution, no discretion.
    Determined,
    /// Plural under `Strict`, the tie-break picks.
    TieBreakPicks,
    /// The tie survives the key.
    TieSurvives,
    /// Fails at every level.
    Contradictory,
    /// The row does not parse or raise.
    ParseFailure,
    /// An outcome combination the tie-break hierarchy forbids.
    HierarchyViolation,
}

impl RowClass {
    fn label(self) -> &'static str {
        match self {
            Self::Determined => "determined",
            Self::TieBreakPicks => "tie-break picks",
            Self::TieSurvives => "tie survives",
            Self::Contradictory => "contradictory",
            Self::ParseFailure => "parse failure",
            Self::HierarchyViolation => "hierarchy violation",
        }
    }
}

fn ingest(smiles: &str, tie_break: ValenceTieBreak) -> Result<Molecule, SmilesInputError> {
    ingest_smiles_with(
        smiles,
        &SmilesIoConfig::opensmiles(),
        &ChemistryModel {
            valence: ValenceModel {
                tie_break,
                ..ValenceModel::smiles()
            },
            ..ChemistryModel::default()
        },
        &ResolveConfig::default(),
    )
}

fn classify(
    strict: &Result<Molecule, SmilesInputError>,
    most_saturated: &Result<Molecule, SmilesInputError>,
) -> RowClass {
    match (strict, most_saturated) {
        (Ok(a), Ok(b)) if a == b => RowClass::Determined,
        (Err(SmilesInputError::Underdetermined(_)), Ok(_)) => RowClass::TieBreakPicks,
        (Err(SmilesInputError::Underdetermined(_)), Err(SmilesInputError::Underdetermined(_))) => {
            RowClass::TieSurvives
        }
        (Err(SmilesInputError::Contradiction(_)), Err(SmilesInputError::Contradiction(_))) => {
            RowClass::Contradictory
        }
        (
            Err(SmilesInputError::Syntax(_) | SmilesInputError::ModelConversion(_)),
            Err(SmilesInputError::Syntax(_) | SmilesInputError::ModelConversion(_)),
        ) => RowClass::ParseFailure,
        _ => RowClass::HierarchyViolation,
    }
}

/// The two per-row invariants: every emitted system re-validates under the
/// resolution model's aromaticity rule (perception `derive` on the resolved
/// molecule reassesses each stored system from its electron contributions),
/// and no two systems overlap on an atom.
fn system_violations(molecule: &Molecule, cell: &str, regid: &str) -> Vec<String> {
    let mut violations = Vec::new();
    let perception = AromaticityPerception::new(&AromaticityModel::daylight());
    match perception.derive(molecule, AromaticityConfig::default()) {
        Ok(Solution::Determined(derivation)) => {
            for inconsistency in derivation.inconsistencies {
                violations.push(format!("{regid} [{cell}] rule: {inconsistency:?}"));
            }
        }
        other => violations.push(format!("{regid} [{cell}] rule: {other:?}")),
    }
    let mut seen = BTreeMap::new();
    for system in molecule.aromatic_systems().iter() {
        for atom in system.atom_ids() {
            if let Some(previous) = seen.insert(atom, system.id) {
                violations.push(format!(
                    "{regid} [{cell}] overlap: atom {atom:?} in {previous:?} and {:?}",
                    system.id
                ));
            }
        }
    }
    violations
}

#[test]
#[ignore = "local instrument; requires materials/aromaticity/vehicle/VEHICLe.csv"]
fn vehicle_instrument() {
    let staging =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../materials/aromaticity/vehicle");
    let Ok(corpus) = fs::read_to_string(staging.join("VEHICLe.csv")) else {
        eprintln!("VEHICLe.csv absent; run scripts/fetch-vehicle.sh to stage the corpus");
        return;
    };

    let mut class_counts: BTreeMap<RowClass, usize> = BTreeMap::new();
    let mut manifest = String::new();
    let mut invariant_violations: Vec<String> = Vec::new();
    let mut hierarchy_rows: Vec<String> = Vec::new();
    // cluster id -> (regid, class, most-saturated result)
    let mut clusters: BTreeMap<String, Vec<(String, RowClass, Option<Molecule>)>> = BTreeMap::new();
    let mut row_count = 0usize;

    for line in corpus.lines().skip(1) {
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split(',').collect();
        let (regid, smiles) = (fields[0], fields[1]);
        let cluster = fields.get(5).copied().unwrap_or("").trim();
        row_count += 1;

        let strict = ingest(smiles, ValenceTieBreak::Strict);
        let most_saturated = ingest(smiles, ValenceTieBreak::MostSaturated);

        let class = classify(&strict, &most_saturated);
        *class_counts.entry(class).or_default() += 1;

        for (cell, outcome) in [("strict", &strict), ("most-saturated", &most_saturated)] {
            if let Ok(molecule) = outcome {
                invariant_violations.extend(system_violations(molecule, cell, regid));
            }
        }
        match class {
            RowClass::HierarchyViolation => {
                hierarchy_rows.push(format!(
                    "{regid} {smiles}: strict {strict:?}, most-saturated {most_saturated:?}"
                ));
            }
            RowClass::Contradictory | RowClass::ParseFailure => {
                let _ = writeln!(
                    manifest,
                    "{} {regid} {smiles}: {}",
                    class.label(),
                    strict
                        .as_ref()
                        .err()
                        .map_or_else(String::new, |e| e.to_string()),
                );
            }
            _ => {}
        }
        if !cluster.is_empty() {
            clusters.entry(cluster.to_owned()).or_default().push((
                regid.to_owned(),
                class,
                most_saturated.ok(),
            ));
        }
    }

    // Tautomer-cluster oracle: rows of one cluster share a class; among the
    // resolved most-saturated members, count distinct canonical forms.
    let context = CanonicalizationContext {
        para_stereo: false,
        automorphism_algorithm: AutomorphismAlgorithm::Nauty,
    };
    let mut cluster_class_splits = 0usize;
    let mut cluster_form_counts: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for (cluster, members) in &clusters {
        if members.len() < 2 {
            continue;
        }
        let classes: Vec<RowClass> = {
            let mut classes: Vec<RowClass> = members.iter().map(|(_, class, _)| *class).collect();
            classes.sort_unstable();
            classes.dedup();
            classes
        };
        if classes.len() > 1 {
            cluster_class_splits += 1;
            let _ = writeln!(
                manifest,
                "cluster {cluster} class split: {:?}",
                members
                    .iter()
                    .map(|(regid, class, _)| format!("{regid}={}", class.label()))
                    .collect::<Vec<_>>()
            );
        }
        let mut forms: Vec<Molecule> = Vec::new();
        let mut resolved = 0usize;
        for (_, _, molecule) in members {
            let Some(molecule) = molecule else { continue };
            resolved += 1;
            let canonical = molecule
                .clone()
                .canonicalize(&context)
                .expect("resolved molecule canonicalizes");
            if !forms.contains(&canonical) {
                forms.push(canonical);
            }
        }
        if resolved > 0 {
            *cluster_form_counts
                .entry((resolved, forms.len()))
                .or_default() += 1;
        }
    }

    let mut distribution = String::new();
    let _ = writeln!(distribution, "rows: {row_count}");
    for (class, count) in &class_counts {
        let _ = writeln!(
            distribution,
            "{}: {count} ({:.2}%)",
            class.label(),
            100.0 * *count as f64 / row_count as f64
        );
    }
    let _ = writeln!(distribution, "clusters (>=2 members): {}", clusters.len());
    let _ = writeln!(distribution, "cluster class splits: {cluster_class_splits}");
    let _ = writeln!(
        distribution,
        "cluster (resolved members, distinct canonical forms) counts: {cluster_form_counts:?}"
    );
    let _ = writeln!(
        distribution,
        "invariant violations: {}",
        invariant_violations.len()
    );
    let _ = writeln!(
        distribution,
        "hierarchy violations: {}",
        hierarchy_rows.len()
    );

    for violation in &invariant_violations {
        let _ = writeln!(manifest, "invariant {violation}");
    }
    for row in &hierarchy_rows {
        let _ = writeln!(manifest, "hierarchy {row}");
    }

    fs::write(staging.join("distribution.txt"), &distribution)
        .expect("staging directory is writable");
    fs::write(staging.join("manifest.txt"), &manifest).expect("staging directory is writable");
    println!("{distribution}");

    assert!(
        hierarchy_rows.is_empty(),
        "tie-break hierarchy violated; see manifest.txt"
    );
    assert!(
        invariant_violations.is_empty(),
        "aromatic-system invariants violated; see manifest.txt"
    );
}
