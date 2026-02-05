//! CXSMILES annotation data for roundtripping.
//!
//! Contains format-specific data that doesn't have clean semantic representation
//! but is needed for faithful roundtripping of CXSMILES.

use std::collections::BTreeMap;

use super::sgroup::SGroup;

/// CXSMILES annotation data
#[derive(Clone, Debug, Default, PartialEq)]
pub struct CxAnnotationData {
    /// Enhanced stereo groups: index -> set of atoms with their mode
    pub stereo_groups: BTreeMap<u32, StereoSet>,

    /// Component groupings (atom indices per component)
    /// Used when explicit grouping differs from graph connectivity
    pub components: Option<Vec<Vec<u32>>>,

    /// S-groups from CXSMILES Sg/SgD/SgH tags
    pub sgroups: BTreeMap<u32, SGroup>,
}

/// A set of stereocenters with a common interpretation mode
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoSet {
    pub atoms: Vec<u32>,
    pub mode: StereoSetMode,
}

/// How to interpret a group of stereocenters
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StereoSetMode {
    /// Centers flip together (racemate-like)
    Correlated,
    /// Centers flip independently (mixture)
    Independent,
}
