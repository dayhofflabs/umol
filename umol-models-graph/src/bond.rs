//! Bond type for the molecular graph model.

use std::collections::HashMap;

/// Represents bond order, mapping common MOL V2000 codes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondType {
    Single,   // MOL code 1
    Double,   // MOL code 2
    Triple,   // MOL code 3
    Aromatic, // MOL code 4
    Other,    // Placeholder for less common types
    SingleOrDouble, // MOL code 5
    SingleOrAromatic, // MOL code 6
    DoubleOrAromatic, // MOL code 7
    Any, // MOL code 8
}

/// Represents double bond stereochemistry specified in MOL V2000 files.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondStereo {
    Cis,   // MOL code 1
    Trans, // MOL code 6
    Either, // MOL code 3
}

/// Represents single bond wedging specified in MOL V2000 files.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondDir {
    Wedge, // MOL code 1 (Up / Begin Wedge)
    Dash,  // MOL code 6 (Down / Begin Dash)
}

/// Represents a bond in a graph-based molecular model, mirroring key RDKit properties.
#[derive(Debug, Clone)]
pub struct Bond {
    /// The type/order of the bond.
    pub bond_type: BondType,
    /// Double bond stereochemistry (`Cis`/`Trans`), if specified.
    pub stereo: Option<BondStereo>,
    /// Single bond directionality (wedge/dash for depiction), if specified.
    pub dir: Option<BondDir>,
    /// Generic string-based properties.
    pub properties: HashMap<String, String>,
}

impl Bond {
    /// Creates a new Bond with default properties for the given BondType.
    pub fn new(bond_type: BondType) -> Self {
        Self {
            bond_type,
            stereo: None,
            dir: None,
            properties: HashMap::new(),
        }
    }
} 