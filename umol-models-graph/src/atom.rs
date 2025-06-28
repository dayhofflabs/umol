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

/// Include stereochemistry in query atoms
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomStereoCare {
    /// None corresponds to MOL code 0, RDKit `molStereoCare = 0`
    /// Corresponds to MOL code 1, RDKit `molStereoCare = 1`
    Care,
}

/// Inversion/retention flag
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomInversionRetention {
    /// Corresponds to MOL code 1, RDKit `molInversionFlag = 4`
    Inverted,
    /// Corresponds to MOL code 2, RDKit `molInversionFlag = 8`
    Retained,
}

/// Atom exact change flag
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomExactChange {
    /// Corresponds to MOL code 1, RDKit `molExactChangeFlag = 4`
    Match,
}

/// Atom list (for query molecules in MOL files)
#[derive(Debug, Clone, PartialEq)]
pub struct AtomList {
    pub elements: Vec<Element>,
}

/// Generalized atom kind (for atom-like objects in MOL files)
#[derive(Debug, Clone, PartialEq)]
pub enum AtomSymbol {
    Element(Element),
    NamedIsotope(NamedIsotope),
    AtomList(AtomList),
    Unspecified(char),
    LonePair,
    RGroup(usize),
}

/// Atom
#[derive(Debug, Clone)]
pub struct AtomStandard {
    pub element: Element,
    pub charge: i8,
    pub radical: Option<u8>,
    pub isotope_mass: Option<u32>,
    pub stereo_parity: Option<AtomStereoParity>,
    pub hydrogen_count: Option<u8>,
    pub valence: Option<u8>,
    pub atom_map_num: Option<u32>,
    pub properties: HashMap<String, String>,
}

impl AtomStandard {
    /// Create new Atom with default properties for given element
    pub fn new(element: Element) -> Self {
        Self {
            element,
            charge: 0,
            radical: None,
            isotope_mass: None,
            stereo_parity: None,
            hydrogen_count: None,
            valence: None,
            atom_map_num: None,
            properties: HashMap::new(),
        }
    }
}

/// Link atom specification for LIN property
#[derive(Debug, Clone, PartialEq)]
pub struct LinkAtomSpec {
    pub repeat_count: u8,  // vvv >= 2
    pub bond1: usize,      // bbb (can be 0)
    pub bond2: usize,      // ccc (can be 0)
}

/// Generalized atom symbol (for atom-like objects in MOL files)
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub symbol: AtomSymbol,
    pub charge: i8,
    pub radical: Option<u8>,
    pub isotope_mass: Option<u32>,
    pub stereo_parity: Option<AtomStereoParity>,
    pub stereo_care: Option<AtomStereoCare>,
    pub hydrogen_count: Option<u8>,
    pub valence: Option<u8>,
    pub atom_map_num: Option<u32>,
    pub inversion_retention: Option<AtomInversionRetention>,
    pub exact_change: Option<AtomExactChange>,
    pub properties: HashMap<String, String>,
    
    // Query-specific properties
    pub attachment_point: Option<u8>,                    // APO: 0=none, 1=first, 2=second, 3=both
    pub attachment_order: Option<Vec<(usize, u8)>>,      // AAL: neighbor-order pairs
    pub ring_bond_count: Option<i8>,                     // RBC: -2=r*, -1=r0, 0=off, 2+=count
    pub substitution_count: Option<i8>,                  // SUB: -2=s*, -1=s0, 0=off, 1+=count
    pub unsaturated: Option<bool>,                       // UNS: Some(true)=on, None=off
    pub link_atom: Option<LinkAtomSpec>,                 // LIN: complex spec
}

impl Atom {
    pub fn new(symbol: AtomSymbol) -> Self {
        Self {
            symbol,
            charge: 0,
            radical: None,
            isotope_mass: None,
            stereo_parity: None,
            hydrogen_count: None,
            stereo_care: None,
            valence: None,
            atom_map_num: None,
            inversion_retention: None,
            exact_change: None,
            properties: HashMap::new(),
            
            // Query-specific properties - all default to None
            attachment_point: None,
            attachment_order: None,
            ring_bond_count: None,
            substitution_count: None,
            unsaturated: None,
            link_atom: None,
        }
    }
}



impl From<AtomStandard> for Atom {
    fn from(atom: AtomStandard) -> Self {
        Self {
            symbol: AtomSymbol::Element(atom.element),
            charge: atom.charge,
            radical: atom.radical,
            isotope_mass: atom.isotope_mass,
            stereo_parity: atom.stereo_parity,
            stereo_care: None,
            hydrogen_count: atom.hydrogen_count,
            valence: atom.valence,
            atom_map_num: atom.atom_map_num,
            inversion_retention: None,
            exact_change: None,
            properties: atom.properties,
            
            // Query-specific properties - not applicable to AtomStandard
            attachment_point: None,
            attachment_order: None,
            ring_bond_count: None,
            substitution_count: None,
            unsaturated: None,
            link_atom: None,
        }
    }
}
