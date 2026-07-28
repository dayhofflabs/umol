//! Molecule properties grouped by operation.
//!
//! The parent owns one regression file for the whole subject. Every child uses
//! that file through `super::REGRESSION_FILE`, so reorganizing properties by
//! operation does not orphan minimized failures.

const REGRESSION_FILE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/property/molecule.proptest-regressions"
);

#[path = "molecule/compaction.rs"]
mod compaction;
#[path = "molecule/comparison.rs"]
mod comparison;
#[path = "molecule/correspondence.rs"]
mod correspondence;
#[path = "molecule/iterators.rs"]
mod iterators;
#[path = "molecule/meet_pushout.rs"]
mod meet_pushout;
#[path = "molecule/references.rs"]
mod references;
#[path = "molecule/ring.rs"]
mod ring;
#[path = "molecule/serialization.rs"]
mod serialization;
#[path = "molecule/structure.rs"]
mod structure;
