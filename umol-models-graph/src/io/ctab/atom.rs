//! Atom type for CTab format.

use crate::io::ctab::query::QueryAtom;
use crate::io::ctab::rgroup::RGroup;
use std::collections::HashMap;
use umol_data::{Element, NamedIsotope};

/// Tetrahedral chirality specified in MOL files.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomStereoParity {
    Odd,    // sss = 1, RDKit `CHI_TETRAHEDRAL_CW` (Clockwise / R)
    Even,   // sss = 2, RDKit `CHI_TETRAHEDRAL_CCW` (Counter-Clockwise / S)
    Either, // sss = 3, RDKit `CHI_UNSPECIFIED`
}

/// Include stereochemistry in query atoms
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomStereoCare {
    Care, // hhh = 1, RDKit `molStereoCare = 1`
}

/// Inversion/retention flag
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomInversionRetention {
    Inverted, // hhh = 1, RDKit `molInversionFlag = 4`
    Retained, // hhh = 2, RDKit `molInversionFlag = 8`
}

/// Atom exact change flag
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomExactChange {
    Match, // hhh = 1, RDKit `molExactChangeFlag = 4`
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
    Query(QueryAtom),
    LonePair,
    RGroup(RGroup),
}

/// Atom
#[derive(Debug, Clone)]
pub struct AtomStandard {
    pub element: Element,
    pub charge: i8,
    pub radical: Option<AtomRadical>,
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

/// Radical type for RAD property
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AtomRadical {
    Singlet, // RAD vvv = 1
    Doublet, // RAD vvv = 2
    Triplet, // RAD vvv = 3
}

/// Attachment point type for APO property
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AttachmentPointType {
    First,  // APO vvv = 1
    Second, // APO vvv = 2
    Both,   // APO vvv = 3
}

/// Link atom for LIN property
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct LinkAtom {
    pub repeat_count: u8, // LIN vvv >= 2
    pub bond1: usize,     // LIN bbb (can be 0)
    pub bond2: usize,     // LIN ccc (can be 0)
}

/// Generalized atom symbol (for atom-like objects in MOL files)
#[derive(Debug, Clone, PartialEq)]
pub struct Atom {
    pub symbol: AtomSymbol,
    pub charge: i8,
    pub radical: Option<AtomRadical>,
    pub isotope_mass: Option<u32>,
    pub stereo_parity: Option<AtomStereoParity>,
    pub stereo_care: Option<AtomStereoCare>,
    pub hydrogen_count: Option<u8>,
    pub valence: Option<u8>,
    pub atom_map_num: Option<u32>,
    pub inversion_retention: Option<AtomInversionRetention>,
    pub exact_change: Option<AtomExactChange>,
    pub attachment_point: Option<AttachmentPointType>,
    pub attachment_order: Option<Vec<(usize, u8)>>,
    pub ring_bond_count: Option<i8>,
    pub substitution_count: Option<i8>,
    pub unsaturated: Option<bool>,
    pub link_atom: Option<LinkAtom>,
    pub properties: HashMap<String, String>,
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
            attachment_point: None,
            attachment_order: None,
            ring_bond_count: None,
            substitution_count: None,
            unsaturated: None,
            link_atom: None,
            properties: HashMap::new(),
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
            attachment_point: None,
            attachment_order: None,
            ring_bond_count: None,
            substitution_count: None,
            unsaturated: None,
            link_atom: None,
            properties: atom.properties,
        }
    }
}
