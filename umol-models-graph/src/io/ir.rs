//! Intermediate representation for molecular structures

use umol_data::{Element, NamedIsotope};

use crate::io::ctab::atom::AtomRadical;
use crate::io::ctab::bond::{BondReactingCenter, BondStereo, BondTopology};

/// Input molecular format
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum SourceFormat {
    MOL,
    SMILES,
    SMARTS,
    #[default]
    UNKNOWN,
}

/// Molecule IR
///
/// Preserves all information from parsing without chemical interpretation.
/// Can be converted to validated molecular types through ParseTarget trait.
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Molecule {
    // Core structure
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,

    // Fragment/Link architecture for structural organization
    pub fragments: Vec<Fragment>,
    pub links: Vec<Link>,

    // Properties
    pub properties: Vec<Property>,

    // Metadata
    pub comments: Vec<String>,
    pub source_format: SourceFormat,
}

/// Unified atom representation
#[derive(Debug, Default, Clone, PartialEq)]
pub enum AtomSymbol {
    Element(Element),
    NamedIsotope(NamedIsotope),
    Query(QueryAtom),
    Variable(Variable),
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

    pub aromatic: Option<bool>,
    pub chirality: Option<Chirality>,
    pub class: Option<u32>,

    // Metadata
    pub source_format: SourceFormat,
    pub repr: Option<String>,
}

impl Atom {
    /// Create a new atom from an element
    pub fn from_element(element: Element) -> Self {
        Self {
            symbol: AtomSymbol::Element(element),
            repr: Some(element.to_string()),
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

/// Bond IR
#[derive(Debug, Default, Clone, PartialEq)]
pub struct Bond {
    pub atom_indices: (usize, usize),

    // Core properties
    pub symbol: BondSymbol,
    pub stereo: Option<BondStereo>,
    pub direction: Option<BondDirection>,

    // Query properties
    pub topology: Option<BondTopology>,
    pub reacting_center: Option<BondReactingCenter>,

    // Metadata
    pub source_format: SourceFormat,
    pub repr: Option<String>,
}

/// Extended query atom types (superset of MOL + SMARTS)
#[derive(Debug, Clone, PartialEq)]
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
}

/// Unified bond representation (concrete + queries)
#[derive(Debug, Default, Clone, PartialEq)]
pub enum BondSymbol {
    Bond(BondOrder),
    Query(QueryBond),
    #[default]
    Unknown,
}

#[derive(Debug, Clone, PartialEq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
}

#[derive(Debug, Clone, PartialEq)]
pub enum QueryBond {
    SingleOrDouble,
    SingleOrAromatic,
    DoubleOrAromatic,
    Any,
}

/// Bond direction/wedging information
#[derive(Debug, Clone, PartialEq)]
pub enum BondDirection {
    Up,
    Down,
    Either,
}

/// Chirality
#[derive(Debug, Default, Clone, PartialEq)]
pub enum Chirality {
    Clockwise,
    CounterClockwise,
    Tetrahedral {
        index: u32,
    },
    Allenal {
        index: u32,
    },
    SquarePlanar {
        index: u32,
    },
    TrigonalBipyramidal {
        index: u32,
    },
    Octahedral {
        index: u32,
    },
    #[default]
    Unknown,
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
