//! SGroup (Superatom group)

/// SGroup type
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SGroupType {
    Generic,         // GEN
    MultipleGroup,   // MUL
    SRU,             // SRU (Structural Repeating Unit)
    Superatom,       // SUP (Superatom/Abbreviation)
    Data,            // DAT
    Unknown(String), // Placeholder for unhandled types
}

/// SGroup (Superatom group)
#[derive(Debug, Clone)]
pub struct SGroup {
    pub id: usize, // The original index from the file (sss)
    pub label: Option<String>, // Label for SUP, SRU, etc.
    pub subscript: Option<String>, // Subscript text (e.g., "n", "2")
    pub group_type: SGroupType,
    pub atom_indices: Vec<usize>, // 1-based MOL indices of atoms in the group
    // Storing bond endpoints (atom pairs) for now
    pub bond_endpoint_pairs: Vec<(usize, usize)>,
    // TODO: Add fields for other SGroup properties as needed (subtype, label, connectivity, etc.)
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
            bond_endpoint_pairs: Vec::new(),
        }
    }
}
