//! Atom type for the molecular graph model.

use std::collections::HashMap;
use umol_data::Element;

/// Represents tetrahedral chirality specified in MOL files.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomStereoParity {
    /// Corresponds to MOL code 1, RDKit `CHI_TETRAHEDRAL_CW` (Clockwise / R).
    Odd,
    /// Corresponds to MOL code 2, RDKit `CHI_TETRAHEDRAL_CCW` (Counter-Clockwise / S).
    Even,
    /// Corresponds to MOL code 3, RDKit `CHI_UNSPECIFIED`.
    Either,
}

/// Represents an atom in a graph-based molecular model, mirroring key RDKit properties.
#[derive(Debug, Clone)]
pub struct Atom {
    /// The chemical element.
    pub element: Element,
    /// Formal charge on the atom.
    pub formal_charge: i8,
    /// Isotope mass difference relative to the standard atomic weight for the element.
    /// `None` or `Some(0)` represents the default isotope.
    pub mass_difference: Option<i8>,
    /// Tetrahedral chirality, if specified.
    pub stereo_parity: Option<AtomStereoParity>,
    /// Number of explicit hydrogens attached, if specified (e.g., from HCOUNT).
    pub explicit_hydrogens: Option<u8>,
    /// Valence specified directly in the input (e.g., MOL VAL field).
    pub valence: Option<u8>,
    /// Atom mapping number, often used in reactions.
    pub atom_map_num: Option<u32>,
    /// Radical status: 1=singlet, 2=doublet, 3=triplet. `None` or `Some(0)` means non-radical.
    pub radical: Option<u8>,
    /// Generic string-based properties.
    pub properties: HashMap<String, String>,
}

impl Atom {
    /// Creates a new Atom with default properties for the given element.
    pub fn new(element: Element) -> Self {
        Self {
            element,
            formal_charge: 0,
            mass_difference: None,
            stereo_parity: None,
            explicit_hydrogens: None,
            valence: None,
            atom_map_num: None,
            radical: None,
            properties: HashMap::new(),
        }
    }
}
