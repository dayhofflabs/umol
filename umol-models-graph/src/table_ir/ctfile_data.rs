//! CTFile-specific data for roundtripping.
//!
//! Contains data structures that are specific to CTFile formats (MOL, SDF).
//! This data is preserved for exact roundtripping, scheduled to be replaced
//! by semantically defined structures.

use std::collections::BTreeMap;

use super::rgroup::RGroup;
use super::sgroup::SGroup;

/// Legacy group abbreviation
#[derive(Debug, Clone, PartialEq)]
pub struct LegacyGroupAbbreviation {
    pub atom_index1: usize, // Atoms on this side are abbreviated
    pub atom_index2: usize, // Attachment point to main structure
    pub label: String,
}

/// CTFile-specific data container for ExtendedMolecule
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CtfileData {
    pub sgroups: BTreeMap<usize, SGroup>,
    pub rgroups: BTreeMap<usize, RGroup>,
    pub legacy_group_abbreviations: Vec<LegacyGroupAbbreviation>,
}
