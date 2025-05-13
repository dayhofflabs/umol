//! Atom type for the molecular graph model.

use std::collections::HashMap;
use umol_data::Element;

/// Tetrahedral chirality specified in MOL files.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomStereoParity {
    /// Corresponds to MOL code 1, RDKit `CHI_TETRAHEDRAL_CW` (Clockwise / R).
    Odd,
    /// Corresponds to MOL code 2, RDKit `CHI_TETRAHEDRAL_CCW` (Counter-Clockwise / S).
    Even,
    /// Corresponds to MOL code 3, RDKit `CHI_UNSPECIFIED`.
    Either,
}

/// Atom
#[derive(Debug, Clone)]
pub struct Atom {
    /// Element
    pub element: Element,
    /// Charge
    pub charge: i8,
    /// Isotope mass number
    pub isotope_mass: Option<u32>,
    /// Tetrahedral chirality
    pub stereo_parity: Option<AtomStereoParity>,
    /// Hydrogen count
    pub hydrogen_count: Option<u8>,
    /// Valence
    pub valence: Option<u8>,
    /// Atom mapping number
    pub atom_map_num: Option<u32>,
    /// Radical flag
    pub radical: Option<u8>,
    /// Generic string-based properties.
    pub properties: HashMap<String, String>,
}

impl Atom {
    /// Create new Atom with default properties for given element
    pub fn new(element: Element) -> Self {
        Self {
            element,
            charge: 0,
            isotope_mass: None,
            stereo_parity: None,
            hydrogen_count: None,
            valence: None,
            atom_map_num: None,
            radical: None,
            properties: HashMap::new(),
        }
    }
}
