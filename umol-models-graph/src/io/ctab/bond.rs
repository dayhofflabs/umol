//! Bond type for CTab format.

use bitflags::bitflags;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Bond order, mapping common MOL V2000 codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondType {
    Single,           // MOL code 1
    Double,           // MOL code 2
    Triple,           // MOL code 3
    Aromatic,         // MOL code 4
    SingleOrDouble,   // MOL code 5
    SingleOrAromatic, // MOL code 6
    DoubleOrAromatic, // MOL code 7
    Any,              // MOL code 8
    Zero,             // zero-order bond, via ZBO property
}

impl BondType {
    pub fn is_standard(&self) -> bool {
        matches!(self, BondType::Single | BondType::Double | BondType::Triple | BondType::Aromatic)
    }
}

/// Double bond stereochemistry specified in MOL V2000 files.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondStereo {
    Cis,    // MOL code 1
    Trans,  // MOL code 6
    Either, // MOL code 3
}

impl Default for BondStereo {
    fn default() -> Self {
        BondStereo::Either
    }
}

impl BondStereo {
    pub fn is_default(&self) -> bool {
        matches!(self, BondStereo::Either)
    }
}

/// Single bond wedging specified in MOL V2000 files.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondDir {
    Wedge,  // MOL code 1 (Up / Begin Wedge)
    Dash,   // MOL code 6 (Down / Begin Dash)
    Either, // MOL code 4 (Either)
}

impl Default for BondDir {
    fn default() -> Self {
        BondDir::Either
    }
}

impl BondDir {
    pub fn is_default(&self) -> bool {
        matches!(self, BondDir::Either)
    }
}

/// Bond topology (chain, ring, either), if specified.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondTopology {
    Chain,  // MOL code 2
    Ring,   // MOL code 1
    Either, // MOL code 0 (default/unspecified)
}

impl Default for BondTopology {
    fn default() -> Self {
        BondTopology::Either
    }
}

impl BondTopology {
    pub fn is_default(&self) -> bool {
        matches!(self, BondTopology::Either)
    }
}

bitflags! {
    /// Bond reacting center status
    ///
    /// - `0`: Unmarked (default/unspecified)
    /// - `1`: A reacting center
    /// - `-1`: Not a reacting center (exclusive, cannot be combined)
    /// - `2`: No change in the bond (exclusive, cannot be combined)
    /// - `4`: Bond is made or broken during the reaction
    /// - `8`: Bond order changes during the reaction
    ///
    /// Some combinations are possible.
    /// The `CENTER` flag can be combined with `MADE_BROKEN` and/or `ORDER_CHANGED`.
    /// `NOT_CENTER` and `NO_CHANGE` are exclusive.
    #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
    pub struct BondReactingCenter: i16 {
        const UNMARKED         = 0b00000000; // MOL code 0 (default/unspecified)
        const CENTER           = 0b00000001; // MOL code 1
        const NOT_CENTER       = 0b00000010; // MOL code -1 (exclusive, cannot be combined)
        const NO_CHANGE        = 0b00000100; // MOL code 2 (exclusive, cannot be combined)
        const MADE_BROKEN      = 0b00001000; // MOL code 4
        const ORDER_CHANGED    = 0b00010000; // MOL code 8

        const MADE_BROKEN_AND_ORDER_CHANGED = Self::MADE_BROKEN.bits() | Self::ORDER_CHANGED.bits();
        const CENTER_AND_MADE_BROKEN = Self::CENTER.bits() | Self::MADE_BROKEN.bits();
        const CENTER_AND_ORDER_CHANGED = Self::CENTER.bits() | Self::ORDER_CHANGED.bits();
        const CENTER_AND_MADE_BROKEN_AND_ORDER_CHANGED = Self::CENTER.bits() | Self::MADE_BROKEN.bits() | Self::ORDER_CHANGED.bits();
    }
}

impl Default for BondReactingCenter {
    fn default() -> Self {
        BondReactingCenter::UNMARKED
    }
}

impl BondReactingCenter {
    pub fn is_default(&self) -> bool {
        *self == BondReactingCenter::UNMARKED
    }
}

/// Bond
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BondStandard {
    pub bond_type: BondType,
    pub stereo: Option<BondStereo>,
    pub dir: Option<BondDir>,
    pub properties: HashMap<String, String>,
}

impl BondStandard {
    pub fn new(bond_type: BondType) -> Self {
        Self {
            bond_type,
            stereo: None,
            dir: None,
            properties: HashMap::new(),
        }
    }
}

/// Bond
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Bond {
    pub bond_type: BondType,
    pub stereo: Option<BondStereo>,
    pub dir: Option<BondDir>,
    pub topology: Option<BondTopology>,
    pub reacting_center: Option<BondReactingCenter>,
    pub properties: HashMap<String, String>,
}

impl Bond {
    pub fn new(bond_type: BondType) -> Self {
        Self {
            bond_type,
            stereo: None,
            dir: None,
            topology: None,
            reacting_center: None,
            properties: HashMap::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bond_standard_serialize() {
        let mut bond = BondStandard::new(BondType::Double);
        bond.properties
            .insert("test_key".to_string(), "test_value".to_string());

        // Test YAML serialization
        let yaml = serde_yaml::to_string(&bond).expect("Failed to serialize BondStandard to YAML");
        let deserialized: BondStandard =
            serde_yaml::from_str(&yaml).expect("Failed to deserialize BondStandard from YAML");
        assert_eq!(bond, deserialized);
    }

    #[test]
    fn test_bond_serialize() {
        let bond = Bond::new(BondType::Double);
        let yaml = serde_yaml::to_string(&bond).expect("Failed to serialize Bond to YAML");
        let deserialized: Bond =
            serde_yaml::from_str(&yaml).expect("Failed to deserialize Bond from YAML");
        assert_eq!(bond, deserialized);
    }
}
