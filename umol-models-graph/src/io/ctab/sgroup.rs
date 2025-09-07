//! SGroup types for CTab format.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

/// SGroup type
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
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
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupSubtype {
    Alternating, // ALT
    Random,      // RAN
    Block,       // BLO
}

/// SGroup connectivity types
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupConnectivity {
    HeadToHead,    // HH
    HeadToTail,    // HT
    EitherUnknown, // EU
}

/// SGroup multiplier term (variable or integer)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupMultiplierTerm {
    Variable(char),
    Integer(u32),
}

/// SGroup multiplier arithmetic operator
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupMultiplierOp {
    Add,
    Sub,
    Mul,
    Div,
}

/// SGroup multiplier for repeating unit properties
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupMultiplier {
    Single(SGroupMultiplierTerm),
    Expression {
        left: SGroupMultiplierTerm,
        op: SGroupMultiplierOp,
        right: SGroupMultiplierTerm,
    },
}

/// SGroup bracket coordinates
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SGroupBracketCoords {
    pub bracket1: (f64, f64),
    pub bracket2: (f64, f64),
}

/// SGroup connecting bond
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SGroupConnectingBond {
    pub bond_index: usize,
    pub bond_vector: (f64, f64),
}

/// SGroup bracket style
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupBracketStyle {
    #[default]
    Default, // 0 = default brackets
    Curved,  // 1 = curved (parenthetic) brackets
}

/// SGroup data type
#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
pub enum SGroupDataType {
    Formatted,
    Numeric,
    #[default]
    Text,
}

/// SGroup data
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SGroupData {
    pub field_type: SGroupDataType,
    pub field_units: Option<String>,
    pub query_identifier: Option<String>,
    pub data_query_operator: Option<String>,
    pub data_content: Option<Vec<String>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupDataDisplayType {
    Attached,
    Detached,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupDataDisplayPlacement {
    Absolute,
    Relative,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupDataDisplayUnits {
    None,
    DisplayUnits,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SGroupDataDisplayChars {
    All,
    Number(u32),
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SGroupDataDisplay {
    pub coords: (f64, f64),
    pub display_type: SGroupDataDisplayType,
    pub display_placement: SGroupDataDisplayPlacement,
    pub display_units: SGroupDataDisplayUnits,
    pub display_chars: SGroupDataDisplayChars,
}

/// SGroup (Substance group)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
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
    pub data: BTreeMap<String, SGroupData>,       // SDT, SCD, SED: data for DAT SGroups
    pub display: Option<SGroupDataDisplay>,       // SDD: display info for DAT SGroups
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
            data: BTreeMap::new(),
            display: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use pretty_assertions::assert_eq;

    use super::*;

    #[test]
    fn test_sgroup_serialize() {
        let sgroup = SGroup::new(SGroupType::Superatom);
        let yaml = serde_yaml::to_string(&sgroup).expect("Failed to serialize SGroup to YAML");
        let deserialized: SGroup =
            serde_yaml::from_str(&yaml).expect("Failed to deserialize SGroup from YAML");
        assert_eq!(sgroup, deserialized);
    }
}
