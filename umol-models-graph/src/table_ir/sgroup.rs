//! SGroup for TableIR.
// These are temporary types inherited from CTFile
// TODO: Replace by semantically defined structures that SGroup type combines.

/// SGroup type
#[derive(Debug, Clone, Copy, PartialEq)]
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
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SGroupSubtype {
    Alternating, // ALT
    Random,      // RAN
    Block,       // BLO
}

/// SGroup connectivity types
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SGroupConnectivity {
    HeadToHead,    // HH
    HeadToTail,    // HT
    EitherUnknown, // EU
}

/// SGroup multiplier term (variable or integer)
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SGroupMultiplierTerm {
    Variable(char),
    Integer(u32),
}

/// SGroup multiplier arithmetic operator
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SGroupMultiplierOp {
    Add,
    Sub,
    Mul,
    Div,
}

/// SGroup multiplier for repeating unit properties
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SGroupMultiplier {
    Single(SGroupMultiplierTerm),
    Expression {
        left: SGroupMultiplierTerm,
        op: SGroupMultiplierOp,
        right: SGroupMultiplierTerm,
    },
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
    pub bond_index: u32,
    pub bond_vector: (f64, f64),
}

/// SGroup bracket style
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SGroupBracketStyle {
    Default, // 0 = default brackets
    Curved,  // 1 = curved (parenthetic) brackets
}

/// SGroup data type
#[derive(Debug, Clone, PartialEq)]
pub enum SGroupDataType {
    Formatted,
    Numeric,
    Text,
}

/// SGroup data
#[derive(Debug, Clone, PartialEq)]
pub struct SGroupData {
    pub field_type: SGroupDataType,
    pub field_name: String,
    pub field_units: Option<String>,
    pub query_identifier: Option<String>,
    pub data_query_operator: Option<String>,
    pub data_content: Option<Vec<String>>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SGroupDataDisplayType {
    Attached,
    Detached,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SGroupDataDisplayPlacement {
    Absolute,
    Relative,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SGroupDataDisplayUnits {
    None,
    DisplayUnits,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SGroupDataDisplayChars {
    All,
    Number(u32),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SGroupDataDisplay {
    pub coords: (f64, f64),
    pub display_type: SGroupDataDisplayType,
    pub display_placement: SGroupDataDisplayPlacement,
    pub display_units: SGroupDataDisplayUnits,
    pub display_chars: SGroupDataDisplayChars,
}

/// SGroup (Substance group)
#[derive(Debug, Clone, PartialEq)]
pub struct SGroup {
    pub label: Option<u32>,
    pub subscript: Option<String>,
    pub group_type: SGroupType,
    pub group_subtype: Option<SGroupSubtype>,
    pub connectivity: Option<SGroupConnectivity>,
    pub expansion: bool,
    pub atom_indices: Vec<u32>,
    pub bond_indices: Vec<u32>,
    pub parent_atom_indices: Option<Vec<u32>>,
    pub correspondence: Option<Vec<u32>>,
    pub connecting_bond: Option<SGroupConnectingBond>,
    pub hierarchy_parent: Option<u32>,
    pub component_number: Option<u32>,
    pub bracket_style: Option<SGroupBracketStyle>,
    pub data: Option<SGroupData>,
    pub multiplier: Option<SGroupMultiplier>,
    pub bracket_coords: Option<SGroupBracketCoords>,
    pub display: Option<SGroupDataDisplay>,
}

impl SGroup {
    pub fn new(group_type: SGroupType) -> Self {
        Self {
            label: None,
            subscript: None,
            group_type,
            group_subtype: None,
            connectivity: None,
            expansion: false,
            atom_indices: Vec::new(),
            bond_indices: Vec::new(),
            parent_atom_indices: None,
            correspondence: None,
            connecting_bond: None,
            hierarchy_parent: None,
            component_number: None,
            bracket_style: None,
            data: None,
            multiplier: None,
            bracket_coords: None,
            display: None,
        }
    }
}
