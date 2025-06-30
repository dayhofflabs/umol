//! SGroup (Superatom group)

use umol::error::Result;

/// SGroup type
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SGroupType {
    Generic,         // GEN
    MultipleGroup,   // MUL
    RepeatingUnit,   // SRU
    Superatom,       // SUP (Superatom/Abbreviation)
    Data,            // DAT
    Unknown(String), // Placeholder for unhandled types
}

/// SGroup connectivity types
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SGroupConnectivity {
    HeadToHead,      // HH
    HeadToTail,      // HT
    EitherUnknown,   // EU
}

/// SGroup subtype for polymers
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SGroupSubtype {
    Alternating,     // ALT
    Random,          // RAN
    Block,           // BLO
}

/// SGroup bracket style
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SGroupBracketStyle {
    Default,         // 0 = default brackets
    Curved,          // 1 = curved (parenthetic) brackets
}

/// SGroup (Superatom group)
#[derive(Debug, Clone)]
pub struct SGroup {
    pub id: usize,
    pub label: Option<String>,                           // Label for SUP, SRU, etc. (SMT)
    pub subscript: Option<String>,                       // Subscript text (e.g., "n", "2") (SMT)
    pub group_type: SGroupType,                          // STY: SGroup type
    pub subtype: Option<SGroupSubtype>,                  // SST: polymer subtype (ALT, RAN, BLO)
    pub connectivity: Option<SGroupConnectivity>,        // SCN: connectivity (HH, HT, EU)
    pub multiplier: Option<String>,                      // SMT: multiplier for multiple groups
    pub atom_indices: Vec<usize>,                        // SAL: atoms in SGroup
    pub bond_indices: Vec<usize>,                        // SBL: bonds in SGroup
    pub parent_atom_indices: Option<Vec<usize>>,         // SPA: parent atoms for multiple groups
    pub bracket_style: Option<SGroupBracketStyle>,       // SBT: bracket display style
    pub data_field_name: Option<String>,                 // SDT: data field name for DAT SGroups
    pub data_field_info: Option<String>,                 // SDT: data field info for DAT SGroups
    pub data_content: Option<Vec<String>>,               // SCD/SED: actual data content
    pub hierarchy_parent: Option<usize>,                 // SPL: parent SGroup for hierarchies
    pub component_number: Option<u32>,                   // SNC: component order number
}

impl SGroup {
    /// Create new SGroup
    pub fn new(id: usize, group_type: SGroupType) -> Self {
        Self {
            id,
            group_type,
            label: None,
            subscript: None,
            atom_indices: Vec::new(),
            bond_indices: Vec::new(),
            subtype: None,
            connectivity: None,
            multiplier: None,
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
            "DAT" => Ok(SGroupType::Data),
            "GEN" => Ok(SGroupType::Generic),
            _ => Ok(SGroupType::Unknown(input.to_string())),
        }
    }
}
