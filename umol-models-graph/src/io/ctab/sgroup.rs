//! SGroup types for CTab format.

use std::collections::HashMap;

use umol::error::DataError;
use umol::error::Result;

/// SGroup type
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SGroupType {
    Superatom,     // SUP
    MultipleGroup, // MUL
    RepeatingUnit, // SRU
    Monomer,       // MON
    Mer,           // MER
    Copolymer,     // COP
    Crosslink,     // CRO
    Modification,  // MOD
    Graft,         // GRA
    Component,     // COM
    Mixture,       // MIX
    Formulation,   // FOR
    Data,          // DAT
    AnyPolymer,    // ANY
    Generic,       // GEN
}

/// SGroup subtype for polymers
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SGroupSubtype {
    Alternating, // ALT
    Random,      // RAN
    Block,       // BLO
}

/// SGroup connectivity types
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SGroupConnectivity {
    HeadToHead,    // HH
    HeadToTail,    // HT
    EitherUnknown, // EU
}

/// SGroup multiplier
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SGroupMultiplier {
    Count(u32),
    N,
}

/// SGroup bracket coordinates
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SGroupBracketCoords {
    pub bracket1: (f64, f64),
    pub bracket2: (f64, f64),
}

/// SGroup connecting bond
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SGroupConnectingBond {
    pub bond_index: usize,
    pub bond_vector: (f64, f64),
}

/// SGroup bracket style
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SGroupBracketStyle {
    Default, // 0 = default brackets
    Curved,  // 1 = curved (parenthetic) brackets
}

/// SGroup data type
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum SGroupDataType {
    Formatted,
    Numeric,
    Text,
}
/// SGroup data
#[derive(Debug, Clone, PartialEq)]
pub struct SGroupData {
    pub field_type: SGroupDataType,
    pub field_units: Option<String>,
    pub query_identifier: Option<String>,
    pub data_query_operator: Option<String>,
    pub data_content: Option<Vec<String>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SGroupDataDisplayType {
    Attached,
    Detached,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SGroupDataDisplayPlacement {
    Absolute,
    Relative,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SGroupDataDisplayUnits {
    None,
    DisplayUnits,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SGroupDataDisplayChars {
    All,
    Number(u32),
}

/// SGroup (Superatom group)
#[derive(Debug, Clone)]
pub struct SGroup {
    pub label: Option<u32>,                       // Label for SUP, SRU, etc. (SLB)
    pub subscript: Option<String>,                // Subscript text (e.g., "n", "2") (SMT)
    pub group_type: SGroupType,                   // STY: SGroup type
    pub group_subtype: Option<SGroupSubtype>,     // SST: polymer subtype (ALT, RAN, BLO)
    pub connectivity: Option<SGroupConnectivity>, // SCN: connectivity (HH, HT, EU)
    pub expansion: bool,                          // SDS: expansion flag
    pub multiplier: Option<SGroupMultiplier>,     // SMT: multiplier for multiple groups
    pub atom_indices: Vec<usize>,                 // SAL: atoms in SGroup
    pub bond_indices: Vec<usize>,                 // SBL: bonds in SGroup
    pub parent_atom_indices: Option<Vec<usize>>,  // SPA: parent atoms for multiple groups
    pub correspondence: Option<Vec<usize>>,       // CRS: correspondence for crosslinks
    pub connecting_bond: Option<SGroupConnectingBond>, // CBV: connecting bond
    pub bracket_coords: Option<SGroupBracketCoords>, // SDI: bracket display info
    pub hierarchy_parent: Option<usize>,          // SPL: parent SGroup for hierarchies
    pub component_number: Option<u32>,            // SNC: component order number
    pub bracket_style: Option<SGroupBracketStyle>, // SBT: bracket display style
    pub data: HashMap<String, SGroupData>,        // SDT, SCD, SED: data for DAT SGroups
}

impl SGroup {
    /// Create new SGroup
    pub fn new(group_type: SGroupType) -> Self {
        Self {
            label: None,
            subscript: None,
            group_type,
            group_subtype: None,
            connectivity: None,
            expansion: false,
            multiplier: None,
            atom_indices: Vec::new(),
            bond_indices: Vec::new(),
            parent_atom_indices: None,
            correspondence: None,
            connecting_bond: None,
            bracket_coords: None,
            hierarchy_parent: None,
            component_number: None,
            bracket_style: None,
            data: HashMap::new(),
        }
    }

    /// Parse SGroup type string
    pub fn get_type(input: &str) -> Result<SGroupType> {
        match input {
            "SUP" => Ok(SGroupType::Superatom),
            "MUL" => Ok(SGroupType::MultipleGroup),
            "SRU" => Ok(SGroupType::RepeatingUnit),
            "MON" => Ok(SGroupType::Monomer),
            "MER" => Ok(SGroupType::Mer),
            "COP" => Ok(SGroupType::Copolymer),
            "CRO" => Ok(SGroupType::Crosslink),
            "MOD" => Ok(SGroupType::Modification),
            "GRA" => Ok(SGroupType::Graft),
            "COM" => Ok(SGroupType::Component),
            "MIX" => Ok(SGroupType::Mixture),
            "FOR" => Ok(SGroupType::Formulation),
            "DAT" => Ok(SGroupType::Data),
            "ANY" => Ok(SGroupType::AnyPolymer),
            "GEN" => Ok(SGroupType::Generic),
            _ => Err(DataError::InvalidFragment(format!("Unknown SGroup type: {}", input)).into()),
        }
    }

    /// Parse SGroup subtype string
    pub fn get_subtype(input: &str) -> Result<SGroupSubtype> {
        match input {
            "ALT" => Ok(SGroupSubtype::Alternating),
            "RAN" => Ok(SGroupSubtype::Random),
            "BLO" => Ok(SGroupSubtype::Block),
            _ => {
                Err(DataError::InvalidFragment(format!("Unknown SGroup subtype: {}", input)).into())
            }
        }
    }

    /// Parse SGroup connectivity string
    pub fn get_connectivity(input: &str) -> Result<SGroupConnectivity> {
        match input {
            "HH" => Ok(SGroupConnectivity::HeadToHead),
            "HT" => Ok(SGroupConnectivity::HeadToTail),
            "EU" => Ok(SGroupConnectivity::EitherUnknown),
            _ => Err(
                DataError::InvalidFragment(format!("Unknown SGroup connectivity: {}", input))
                    .into(),
            ),
        }
    }

    /// Parse SGroup multiplier string
    pub fn get_multiplier(input: &str) -> Result<SGroupMultiplier> {
        if input == "n" || input == "N" {
            Ok(SGroupMultiplier::N)
        } else {
            match input.parse::<u32>() {
                Ok(count) => Ok(SGroupMultiplier::Count(count)),
                Err(_) => Err(DataError::InvalidFragment(format!(
                    "Invalid SGroup multiplier: {}",
                    input
                ))
                .into()),
            }
        }
    }

    /// Parse SGroup data type
    pub fn get_data_type(input: &str) -> Result<SGroupDataType> {
        match input {
            "F" => Ok(SGroupDataType::Formatted),
            "N" => Ok(SGroupDataType::Numeric),
            "T" => Ok(SGroupDataType::Text),
            _ => Err(DataError::InvalidFragment(format!(
                "Unknown SGroup data type: {}",
                input
            ))
            .into()),
        }
    }

    /// Parse SGroup data display type
    pub fn get_data_display_type(input: &str) -> Result<SGroupDataDisplayType> {
        match input {
            "A" => Ok(SGroupDataDisplayType::Attached),
            "D" => Ok(SGroupDataDisplayType::Detached),
            _ => Err(DataError::InvalidFragment(format!(
                "Unknown SGroup data display type: {}",
                input
            ))
            .into()),
        }
    }

    /// Parse SGroup data display placement
    pub fn get_data_display_placement(input: &str) -> Result<SGroupDataDisplayPlacement> {
        match input {
            "A" => Ok(SGroupDataDisplayPlacement::Absolute),
            "R" => Ok(SGroupDataDisplayPlacement::Relative),
            _ => Err(DataError::InvalidFragment(format!(
                "Unknown SGroup data display placement: {}",
                input
            ))
            .into()),
        }
    }

    /// Parse SGroup data display units
    pub fn get_data_display_units(input: &str) -> Result<SGroupDataDisplayUnits> {
        match input {
            " " => Ok(SGroupDataDisplayUnits::None),
            "U" => Ok(SGroupDataDisplayUnits::DisplayUnits),
            _ => Err(DataError::InvalidFragment(format!(
                "Unknown SGroup data display units: {}",
                input
            ))
            .into()),
        }
    }

    /// Parse SGroup data display chars
    pub fn get_data_display_chars(input: &str) -> Result<SGroupDataDisplayChars> {
        match input {
            "ALL" => Ok(SGroupDataDisplayChars::All),
            _ => match input.parse::<u32>() {
                Ok(count) => Ok(SGroupDataDisplayChars::Number(count)),
                Err(_) => Err(DataError::InvalidFragment(format!(
                    "Invalid SGroup data display chars: {}",
                    input
                ))
                .into()),
            },
        }
    }
}
