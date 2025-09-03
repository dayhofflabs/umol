//! Atom type for CTab format.

use crate::io::ctab::query::QueryAtom;
use crate::io::ctab::rgroup::RGroup;
use nalgebra;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use umol_data::{Element, NamedIsotope};

/// Atom
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Atom {
    pub element: Element,
    pub charge: i8,
    pub radical: Option<AtomRadical>,
    pub isotope_mass: Option<u32>,
    pub stereo_parity: Option<AtomStereoParity>,
    pub hydrogen_count: Option<u8>,
    pub valence: Option<u8>,
    pub atom_map_num: Option<u32>,
    pub position: Option<Point3D>,
    pub properties: HashMap<String, String>,
}

impl Atom {
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
            position: None,
            properties: HashMap::new(),
        }
    }
}

/// Atom-like structure
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AtomLike {
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
    pub ring_bond_count: Option<RingBondCount>,
    pub substitution_count: Option<SubstitutionCount>,
    pub unsaturated: Option<UnsaturatedAtom>,
    pub link_atom: Option<LinkAtom>,
    pub position: Option<Point3D>,
    pub properties: HashMap<String, String>,
}

impl AtomLike {
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
            position: None,
            properties: HashMap::new(),
        }
    }
}

impl From<Atom> for AtomLike {
    fn from(atom: Atom) -> Self {
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
            position: atom.position,
            properties: atom.properties,
        }
    }
}

/// Atom symbol includes both elements and atom-like structures
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum AtomSymbol {
    Element(Element),
    NamedIsotope(NamedIsotope),
    AtomList(AtomList),
    Query(QueryAtom),
    LonePair,
    RGroup(RGroup),
    Pseudoatom(String)
}

impl AtomSymbol {
    pub fn is_atomlike(&self) -> bool {
        !matches!(self, AtomSymbol::Element(_) | AtomSymbol::NamedIsotope(_))
    }
}

/// 3D coordinate type
pub type Point3D = nalgebra::Point3<f64>;

/// Atom properties specified in the atom block

/// Tetrahedral chirality specified in MOL files.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AtomStereoParity {
    Odd,    // sss = 1, RDKit `CHI_TETRAHEDRAL_CW` (Clockwise / R)
    Even,   // sss = 2, RDKit `CHI_TETRAHEDRAL_CCW` (Counter-Clockwise / S)
    Either, // sss = 3, RDKit `CHI_UNSPECIFIED`
}

impl Default for AtomStereoParity {
    fn default() -> Self {
        AtomStereoParity::Either
    }
}

impl AtomStereoParity {
    pub fn is_default(&self) -> bool {
        matches!(self, AtomStereoParity::Either)
    }
}

/// Include stereochemistry in query atoms
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AtomStereoCare {
    Care, // hhh = 1, RDKit `molStereoCare = 1`
}

impl AtomStereoCare {
    pub fn is_default(&self) -> bool {
        false
    }
}

/// Inversion/retention flag
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AtomInversionRetention {
    Inverted, // hhh = 1, RDKit `molInversionFlag = 4`
    Retained, // hhh = 2, RDKit `molInversionFlag = 8`
}

impl AtomInversionRetention {
    pub fn is_default(&self) -> bool {
        false
    }
}

/// Atom exact change flag
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AtomExactChange {
    Match, // hhh = 1, RDKit `molExactChangeFlag = 4`
}

impl AtomExactChange {
    pub fn is_default(&self) -> bool {
        false
    }
}

/// Atom properties specified in the properties block

/// Radical type for RAD property
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AtomRadical {
    Singlet, // RAD vvv = 1
    Doublet, // RAD vvv = 2
    Triplet, // RAD vvv = 3
}

impl AtomRadical {
    pub fn is_default(&self) -> bool {
        false
    }
}

/// Attachment point type for APO property
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AttachmentPointType {
    First,  // APO vvv = 1
    Second, // APO vvv = 2
    Both,   // APO vvv = 3
}

/// Ring bond count for RBC property
/// -2 = as drawn (r*), -1 = no ring bonds (r0), 0 = off, 2 = r2, 3 = r3, 4 = r4+
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum RingBondCount {
    AsDrawn,     // r*
    NoRingBonds, // r0
    R2,          // 2
    R3,          // 3
    R4Plus,      // 4+
}

/// Substitution count for SUB property
/// -2 = as drawn (s*), -1 = no substitution (s0), 0 = off, 1-5 = s1-s5, 6 = s6+
/// Extended range: 6-10 = s6-s10
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SubstitutionCount {
    AsDrawn,        // s*
    NoSubstitution, // s0
    S1,             // 1
    S2,             // 2
    S3,             // 3
    S4,             // 4
    S5,             // 5
    S6Plus,         // 6+
    S6,             // 6
    S7,             // 7
    S8,             // 8
    S9,             // 9
    S10,            // 10
}

/// Unsaturated atom for UNS property
/// 0 = off, 1 = on
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct UnsaturatedAtom;

/// Link atom for LIN property
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct LinkAtom {
    pub repeat_count: u8,           // LIN vvv >= 2
    pub subs_index1: usize,         // LIN bbb
    pub subs_index2: Option<usize>, // LIN ccc (optional)
}

/// Atom list
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AtomList {
    pub elements: Vec<Element>,
    pub exclusion: bool,
}

impl Default for AtomList {
    fn default() -> Self {
        Self {
            elements: vec![],
            exclusion: false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use umol_data::Element;

    #[test]
    fn test_atom_serialize() {
        let mut atom = Atom::new(Element::C);
        atom.charge = -1;
        atom.radical = Some(AtomRadical::Doublet);
        atom.hydrogen_count = Some(3);
        atom.position = Some(nalgebra::Point3::new(1.5, -2.3, 0.8));

        let yaml = serde_yaml::to_string(&atom).expect("Failed to serialize Atom to YAML");
        let deserialized: Atom =
            serde_yaml::from_str(&yaml).expect("Failed to deserialize Atom from YAML");
        assert_eq!(atom, deserialized);
    }

    #[test]
    fn test_atomlike_serialize() {
        let atom = AtomLike::new(AtomSymbol::Element(Element::C));
        let yaml = serde_yaml::to_string(&atom).expect("Failed to serialize AtomLike to YAML");
        let deserialized: AtomLike =
            serde_yaml::from_str(&yaml).expect("Failed to deserialize AtomLike from YAML");
        assert_eq!(atom, deserialized);
    }
}
