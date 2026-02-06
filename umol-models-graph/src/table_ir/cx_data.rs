//! CXSMILES annotation data for roundtripping.
//!
//! Contains format-specific data that doesn't have clean semantic representation
//! but is needed for faithful roundtripping of CXSMILES.

use std::collections::BTreeMap;

use super::atom::{BicycloStereo, Chirality};
use super::rgroup::RGroup;
use super::sgroup::SGroup;

/// Local parity entry (@: / @@:). Chiral center with ordered substituents.
#[derive(Clone, Debug, PartialEq)]
pub struct LocalParityEntry {
    pub center: u32,
    pub substituents: Vec<u32>,
    pub chirality: Chirality,
}

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

    /// R-groups from CXSMILES LOG: tag (label, occurrence, rgroup_or_h)
    pub rgroups: BTreeMap<u32, RGroup>,

    /// R-group member structures from RG: tag (label -> SMILES strings)
    pub rgroup_members: BTreeMap<u32, Vec<String>>,

    /// Local parity from @: / @@: (chiral center, ordered substituents, chirality)
    pub local_parity: Option<Vec<LocalParityEntry>>,

    /// Bicyclic stereo from THB: / TLB: / TEB:
    pub bicyclo_stereo: Option<Vec<BicycloStereo>>,
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
