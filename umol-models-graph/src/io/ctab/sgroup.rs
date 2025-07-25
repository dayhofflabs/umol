//! SGroup types for CTab format.

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

/// SGroup bracket style
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SGroupBracketStyle {
    Default, // 0 = default brackets
    Curved,  // 1 = curved (parenthetic) brackets
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
    pub multiplier: Option<String>,               // SMT: multiplier for multiple groups
    pub atom_indices: Vec<usize>,                 // SAL: atoms in SGroup
    pub bond_indices: Vec<usize>,                 // SBL: bonds in SGroup
    pub parent_atom_indices: Option<Vec<usize>>,  // SPA: parent atoms for multiple groups
    pub bracket_style: Option<SGroupBracketStyle>, // SBT: bracket display style
    pub data_field_name: Option<String>,          // SDT: data field name for DAT SGroups
    pub data_field_info: Option<String>,          // SDT: data field info for DAT SGroups
    pub data_content: Option<Vec<String>>,        // SCD/SED: actual data content
    pub hierarchy_parent: Option<usize>,          // SPL: parent SGroup for hierarchies
    pub component_number: Option<u32>,            // SNC: component order number
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
            bracket_style: None,
            data_field_name: None, 
            data_field_info: None,
            data_content: None,
            hierarchy_parent: None, 
            component_number: None,
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
}
