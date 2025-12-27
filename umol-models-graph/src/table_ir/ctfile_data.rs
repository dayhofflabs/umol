//! CTFile-specific data for roundtripping.
//!
//! Contains data structures that are specific to CTFile formats (MOL, SDF).
//! This data is preserved for exact roundtripping, scheduled to be replaced
//! by semantically defined structures.

use std::collections::BTreeMap;

use super::rgroup::RGroup;
use super::sgroup::SGroup;

/// CTFile-specific data container for ExtendedMolecule
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CtfileData {
    pub sgroups: BTreeMap<usize, SGroup>,
    pub rgroups: BTreeMap<usize, RGroup>,
}
