//! Intermediate representation for molecular structures

use std::fmt;

use serde::{Deserialize, Serialize};
use umol_data::{Element, NamedIsotope};

/// Input molecular format
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum SourceFormat {
    MOL,
    SMILES,
    SMARTS,
    #[default]
    UNKNOWN,
}

pub mod builder;

/// Molecule IR
///
/// Preserves all information from parsing without chemical interpretation.
/// Can be converted to validated molecular types through ParseTarget trait.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Molecule {
    // Core structure
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub ring_events: Vec<Ring>,

    // Fragment/Link architecture for structural organization
    pub fragments: Vec<Fragment>,
    pub links: Vec<Link>,

    // Properties
    pub properties: Vec<Property>,

    // Metadata
    pub comments: Vec<String>,
    pub source_format: SourceFormat,
}

impl fmt::Display for Molecule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Molecule (atoms: {}, bonds: {})",
            self.atoms.len(),
            self.bonds.len()
        )
    }
}

/// Unified atom representation
#[derive(Debug, Default, Clone, PartialEq)]
pub enum AtomSymbol {
    Element(Element),
    NamedIsotope(NamedIsotope),
    Query(QueryAtom),
    Variable(Variable),
    Pseudoatom(String),
    #[default]
    Unknown,
}

/// Variable atom
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Variable {}

/// Atom IR
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Atom {
    // Core atomic properties (common to all formats)
    pub symbol: AtomSymbol,
    pub position: Option<Point3D>,
    pub charge: Option<i32>,
    pub isotope: Option<u32>,
    pub radical: Option<AtomRadical>,
    pub hydrogen_count: Option<u32>,
    pub implicit_h: bool,

    pub aromatic: Option<bool>,
    pub chirality: Option<Chirality>,
    pub class: Option<u32>,

    // Metadata
    pub span_start: Option<u32>,
    pub span_end: Option<u32>,
    pub source_format: SourceFormat,
}

impl Atom {
    /// Create a new atom from an element
    pub fn from_element(element: Element) -> Self {
        Self {
            symbol: AtomSymbol::Element(element),
            ..Default::default()
        }
    }

    /// Create a new atom with aromatic flag set to false
    pub fn from_aliphatic_atom(element: Element) -> Self {
        Self {
            symbol: AtomSymbol::Element(element),
            aromatic: Some(false),
            ..Default::default()
        }
    }

    /// Create a new atom with aromatic flag set to true
    pub fn from_aromatic_atom(element: Element) -> Self {
        Self {
            symbol: AtomSymbol::Element(element),
            aromatic: Some(true),
            ..Default::default()
        }
    }
}

/// Extended query atom types (superset of MOL + SMARTS)
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub enum QueryAtom {
    Any,           // * = any atom
    Heavy,         // A = all except H
    Heteroatom,    // Q = any heteroatom (all except H, C)
    Halogen,       // X = F, Cl, Br, I
    Metal,         // M = any metal
    HeavyOrH,      // AH = any atom (CXSMILES extension)
    HeteroatomOrH, // QH = Q or H (CXSMILES extension)
    HalogenOrH,    // XH = X or H (CXSMILES extension)
    MetalOrH,      // MH = M or H (CXSMILES extension)
    #[default]
    Unknown,
}

/// Radical type
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AtomRadical {
    Singlet,
    Doublet,
    Triplet,
}

/// Bond IR
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Bond {
    // Core properties
    pub start_atom: u32,
    pub end_atom: u32,
    pub ring: Option<u32>,
    pub symbol: BondSymbol,
    pub stereo: Option<BondStereo>,
    pub direction: Option<BondDir>,

    // Metadata
    pub span_start: Option<u32>,
    pub span_end: Option<u32>,
    pub source_format: SourceFormat,
}

impl Bond {
    pub fn from_order(order: BondOrder) -> Self {
        Self {
            symbol: BondSymbol::Bond(order),
            ..Default::default()
        }
    }

    pub fn up() -> Self {
        Self {
            symbol: BondSymbol::Bond(BondOrder::Single),
            direction: Some(BondDir::Up),
            ..Default::default()
        }
    }

    pub fn down() -> Self {
        Self {
            symbol: BondSymbol::Bond(BondOrder::Single),
            direction: Some(BondDir::Down),
            ..Default::default()
        }
    }
}

/// Unified bond representation (concrete + queries)
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub enum BondSymbol {
    Bond(BondOrder),
    Query(QueryBond),
    #[default]
    Unknown,
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Quadruple,
    Aromatic,
    #[default]
    Unknown,
}

impl BondOrder {
    pub fn symbol(&self) -> &str {
        match self {
            BondOrder::Single => "-",
            BondOrder::Double => "=",
            BondOrder::Triple => "#",
            BondOrder::Quadruple => "$",
            BondOrder::Aromatic => ":",
            BondOrder::Unknown => "?",
        }
    }
}

#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub enum QueryBond {
    SingleOrDouble,
    SingleOrAromatic,
    DoubleOrAromatic,
    Any,
    #[default]
    Unknown,
}

/// Bond direction/wedging information
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub enum BondDir {
    Up,
    Down,
    Either,
    #[default]
    Unknown,
}

/// Double-bond stereochemistry (E/Z) annotation in IR
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondStereo {
    Cis,
    Trans,
    #[default]
    Either,
}

/// Chirality
#[derive(Debug, Default, Clone, Copy, PartialEq)]
pub enum Chirality {
    Clockwise,
    CounterClockwise,
    Tetrahedral {
        arr: u32,
    },
    Allenal {
        arr: u32,
    },
    SquarePlanar {
        arr: u32,
    },
    TrigonalBipyramidal {
        arr: u32,
    },
    Octahedral {
        arr: u32,
    },
    #[default]
    Unknown,
}

/// Ring closure
#[derive(Debug, Default, Clone, PartialEq)]
pub struct RingBond {
    pub index: Option<u32>,
    pub bond: Option<Bond>,
    pub ring: Option<u32>,
}

/// Ring
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ring {
    pub ring_idx: u32,
    pub atom_a: Option<u32>,
    pub atom_b: Option<u32>,
    pub open_start: Option<u32>,
    pub open_end: Option<u32>,
    pub close_start: Option<u32>,
    pub close_end: Option<u32>,
}

/// Fragment
#[derive(Debug, Clone, PartialEq)]
pub struct Fragment {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
}

/// Link
#[derive(Debug, Clone, PartialEq)]
pub struct Link {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
}

/// Property
#[derive(Debug, Clone, PartialEq)]
pub struct Property {
    pub name: String,
    pub value: String,
}

/// 3D coordinate type
#[derive(Debug, Clone, PartialEq)]
pub struct Point3D {
    x: f64,
    y: f64,
    z: f64,
}

impl Default for Point3D {
    fn default() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }
}
