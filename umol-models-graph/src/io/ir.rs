//! Intermediate representation for molecular structures

use crate::io::config::ParsingConfig;
use crate::io::ctab::bond::{BondStereo, BondTopology, BondType};
use crate::io::ctab::sgroup::SGroup;
use umol::Result;
use umol_data::Element;

/// Input molecular format
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SourceFormat {
    MOL,
}

/// Molecule IR
#[derive(Debug, Clone)]
pub struct RawMolecule {
    // TODO: Add non-MOL specific header
    // pub header: Header,
    pub atoms: Vec<RawAtom>,
    pub bonds: Vec<RawBond>,
    // TODO: Add non-MOL specific properties
    // pub properties: Vec<PropertyEntries>,
    pub sgroups: Vec<SGroup>,
}

/// Atom IR
#[derive(Debug, Clone)]
pub struct RawAtom {
    // Core properties (always present)
    pub element: Element,
    pub position: Option<(f64, f64, f64)>,
    pub formal_charge: i8,
    pub isotope: Option<u16>,

    // TODO: Add query system
    pub query_type: Option<String>,
    pub atom_list: Option<Vec<Element>>,
    pub attachment_point: Option<u8>,
    pub ring_bond_count: Option<u8>,
    pub substitution_count: Option<u8>,
    pub unsaturated: Option<bool>,

    pub source_format: SourceFormat,
    pub original_text: Option<String>,
}

/// Bond IR
#[derive(Debug, Clone)]
pub struct RawBond {
    pub atom_indices: (usize, usize),
    pub order: BondType,
    pub stereo: Option<BondStereo>,
    pub topology: Option<BondTopology>,
    pub reacting_center: Option<u8>,

    pub source_format: SourceFormat,
    pub original_text: Option<String>,
}

impl RawMolecule {
    /// Check if the molecule contains any extended features
    pub fn has_extended_features(&self) -> bool {
        self.atoms.iter().any(|a| a.has_extended_features())
            || self.bonds.iter().any(|b| b.has_extended_features())
            || !self.sgroups.is_empty()
    }

    /// Check if the molecule has 3D coordinates
    pub fn has_3d_coordinates(&self) -> bool {
        self.atoms.iter().any(|a| a.position.is_some())
    }
}

impl RawAtom {
    /// Check if this atom has any extended features
    pub fn has_extended_features(&self) -> bool {
        self.query_type.is_some()
            || self.atom_list.is_some()
            || self.attachment_point.is_some()
            || self.ring_bond_count.is_some()
            || self.substitution_count.is_some()
            || self.unsaturated.is_some()
    }
}

impl RawBond {
    /// Check if this bond has any extended features
    pub fn has_extended_features(&self) -> bool {
        self.order.is_bondlike()
        // TODO: Add topology and reacting center checks when those types support it
    }
}

/// Trait for types that can be constructed from parsed molecular data
pub trait ParseTarget: Sized {
    fn allows_query_features() -> bool;
    fn allows_sgroups() -> bool;
    fn allows_rgroups() -> bool;

    fn from_parsed_data(parsed: RawMolecule, config: &ParsingConfig) -> Result<Self>;
}

impl ParsingConfig {
    pub fn for_target<T: ParseTarget>() -> Self {
        Self::default()
    }
}
