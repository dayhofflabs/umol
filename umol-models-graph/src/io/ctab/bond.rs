//! Bond type for CTab format.

use std::collections::HashMap;

use bitflags::bitflags;
use serde::{Deserialize, Serialize};

/// Bond
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Bond {
    pub bond_type: BondType,
    pub stereo: Option<BondStereo>,
    pub dir: Option<BondDir>,
    pub properties: HashMap<String, String>,
}

impl Bond {
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
pub struct BondLike {
    pub bond_type: BondType,
    pub stereo: Option<BondStereo>,
    pub dir: Option<BondDir>,
    pub topology: Option<BondTopology>,
    pub reacting_center: Option<BondReactingCenter>,
    pub properties: HashMap<String, String>,
}

impl BondLike {
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

/// Bond order
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondType {
    Single,           // MOL code 1
    Double,           // MOL code 2
    Triple,           // MOL code 3
    Aromatic,         // MOL code 4
    SingleOrDouble,   // MOL code 5 (query)
    SingleOrAromatic, // MOL code 6 (query)
    DoubleOrAromatic, // MOL code 7 (query)
    Any,              // MOL code 8 (query)
    Zero,             // MOL code 0 (extended range)
    Quadruple,        // MOL code 9 (extended range)
    Quintuple,        // MOL code 10 (extended range)
    Sextuple,         // MOL code 11 (extended range)
}

impl BondType {
    pub fn is_bondlike(&self) -> bool {
        !matches!(
            self,
            BondType::Single
                | BondType::Double
                | BondType::Triple
                | BondType::Aromatic
                | BondType::Quadruple
                | BondType::Quintuple
                | BondType::Sextuple
                | BondType::Zero
        )
    }
    pub fn is_query(&self) -> bool {
        matches!(
            self,
            BondType::SingleOrDouble
                | BondType::SingleOrAromatic
                | BondType::DoubleOrAromatic
                | BondType::Any
        )
    }
    pub fn is_extended_range(&self) -> bool {
        matches!(
            self,
            BondType::Quadruple | BondType::Quintuple | BondType::Sextuple | BondType::Zero
        )
    }
}

/// Bond properties specified in the bond block

/// Double bond stereochemistry
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondStereo {
    Cis,   // MOL code 1
    Trans, // MOL code 6
    #[default]
    Either, // MOL code 3
}

impl BondStereo {
    pub fn is_default(&self) -> bool {
        matches!(self, s if *s == Default::default())
    }
}

/// Single bond wedging specified in MOL V2000 files.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondDir {
    Wedge,  // MOL code 1 (Up / Begin Wedge)
    Dash,   // MOL code 6 (Down / Begin Dash)
    #[default]
    Either, // MOL code 4 (Either)
}

impl BondDir {
    pub fn is_default(&self) -> bool {
        matches!(self, d if *d == Default::default())
    }
}

/// Bond topology (chain, ring, either), if specified.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondTopology {
    Chain,  // MOL code 2
    Ring,   // MOL code 1
    #[default]
    Either, // MOL code 0 (default/unspecified)
}

impl BondTopology {
    pub fn is_default(&self) -> bool {
        matches!(self, t if *t == Default::default())
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bond_serialize() {
        let mut bond = Bond::new(BondType::Double);
        bond.properties
            .insert("test_key".to_string(), "test_value".to_string());

        // Test YAML serialization
        let yaml = serde_yaml::to_string(&bond).expect("Failed to serialize Bond to YAML");
        let deserialized: Bond =
            serde_yaml::from_str(&yaml).expect("Failed to deserialize Bond from YAML");
        assert_eq!(bond, deserialized);
    }

    #[test]
    fn test_bondlike_serialize() {
        let bond = BondLike::new(BondType::Double);
        let yaml = serde_yaml::to_string(&bond).expect("Failed to serialize BondLike to YAML");
        let deserialized: BondLike =
            serde_yaml::from_str(&yaml).expect("Failed to deserialize BondLike from YAML");
        assert_eq!(bond, deserialized);
    }
}
