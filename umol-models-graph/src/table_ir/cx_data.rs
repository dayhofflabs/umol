//! CXSMILES annotation data for roundtripping.
//!
//! Contains format-specific data that doesn't have clean semantic representation
//! but is needed for faithful roundtripping of CXSMILES.

/// CXSMILES annotation data
#[derive(Clone, Debug, Default, PartialEq)]
pub struct CxAnnotationData {
    /// Enhanced stereo: sets of atoms with their stereo mode
    pub stereo_sets: Vec<StereoSet>,

    /// Component groupings (atom indices per component)
    /// Used when explicit grouping differs from graph connectivity
    pub components: Option<Vec<Vec<u32>>>,
}

/// A set of stereocenters with a common interpretation mode
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StereoSet {
    pub atoms: Vec<u32>,
    pub mode: StereoMode,
}

/// How to interpret the stereochemistry of a set of centers
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum StereoMode {
    /// Exactly as drawn
    Absolute,
    /// Centers flip together (racemate-like); group number identifies correlated sets
    Correlated(u32),
    /// Centers flip independently (mixture); group number identifies the set
    Independent(u32),
}
