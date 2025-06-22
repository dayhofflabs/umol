//! Atom type for the molecular graph model.

use std::collections::HashMap;
use umol_data::{Element, NamedIsotope};

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

/// Atom list (for query molecules in MOL files)
#[derive(Debug, Clone)]
pub(crate) struct AtomList {
    pub(crate) elements: Vec<Element>,
}

/// Generalized atom kind (for atom-like objects in MOL files)
#[derive(Debug, Clone)]
pub(crate) enum AtomSymbol {
    Element(Element),
    NamedIsotope(NamedIsotope),
    AtomList(AtomList),
    Unspecified(char),
    LonePair,
    RGroup(usize),
}

/// Atom
#[derive(Debug, Clone)]
pub struct Atom {
    pub element: Element,
    pub charge: i8,
    pub isotope_mass: Option<u32>,
    pub stereo_parity: Option<AtomStereoParity>,
    pub hydrogen_count: Option<u8>,
    pub valence: Option<u8>,
    pub atom_map_num: Option<u32>,
    pub radical: Option<u8>,
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

/// Generalized atom symbol (for atom-like objects in MOL files)
#[derive(Debug, Clone)]
pub struct AtomLike {
    pub symbol: AtomSymbol,
    pub charge: i8,
    pub isotope_mass: Option<u32>,
    pub stereo_parity: Option<AtomStereoParity>,
    pub hydrogen_count: Option<u8>,
    pub valence: Option<u8>,
    pub atom_map_num: Option<u32>,
}

impl AtomLike {
    pub fn new(symbol: AtomSymbol) -> Self {
        Self {
            symbol,
            charge: 0,
            isotope_mass: None,
            stereo_parity: None,
            hydrogen_count: None,
            valence: None,
            atom_map_num: None,
        }
    }
}

impl From<Atom> for AtomLike {
    fn from(atom: Atom) -> Self {
        Self {
            symbol: AtomSymbol::Element(atom.element),
            charge: atom.charge,
            isotope_mass: atom.isotope_mass,
            stereo_parity: atom.stereo_parity,
            hydrogen_count: atom.hydrogen_count,
            valence: atom.valence,
            atom_map_num: atom.atom_map_num,
        }
    }
}
