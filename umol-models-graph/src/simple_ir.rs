//! Simple IR for atom/bond-based molecular models.

use serde::{Deserialize, Serialize};
use umol_data::{Element, NamedIsotope};

use crate::span::Span;

/// Input molecular format
#[derive(Debug, Default, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum SourceFormat {
    MOL,
    SMILES,
    SMARTS,
    #[default]
    UNKNOWN,
}

pub mod builder;

/// Simple molecule IR
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Molecule {
    // Core structure
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub rings: Vec<Ring>,
    pub electrons: Option<u32>,

    // Fragments/links for substructures
    pub fragments: Vec<Fragment>,
    pub links: Vec<Link>,

    // Properties
    pub properties: Vec<Property>,

    // Metadata
    pub comments: Vec<String>,
    pub source_format: SourceFormat,
}

/// Atom IR
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Atom {
    pub symbol: AtomSymbol,
    pub position: Option<Point3D>,
    pub charge: Option<i32>,
    pub isotope: Option<u32>,
    pub unpaired_e: Option<u32>,
    pub hydrogens: Option<u32>,
    pub implicit_h: bool,
    pub aromatic: Option<bool>,
    pub chirality: Option<Chirality>,
    pub class: Option<u32>,

    // Metadata
    pub span: Option<Span>,
}

impl Atom {
    /// Create new aliphatic atom (aromatic flag false)
    pub fn from_aliphatic_atom(element: Element) -> Self {
        Self {
            symbol: AtomSymbol::Element(element),
            aromatic: Some(false),
            implicit_h: true,
            ..Default::default()
        }
    }

    /// Create new aliphatic atom including span
    pub fn from_aliphatic_atom_with_span(
        element: Element,
        span_start: Option<u32>,
        span_end: Option<u32>,
    ) -> Self {
        Self {
            symbol: AtomSymbol::Element(element),
            aromatic: Some(false),
            span: span_start.zip(span_end).map(|(s, e)| Span::bytes(s, e)),
            implicit_h: true,
            ..Default::default()
        }
    }

    /// Create new aromatic atom (aromatic flag true)
    pub fn from_aromatic_atom(element: Element) -> Self {
        Self {
            symbol: AtomSymbol::Element(element),
            aromatic: Some(true),
            implicit_h: true,
            ..Default::default()
        }
    }

    /// Create new aromatic atom including span
    pub fn from_aromatic_atom_with_span(
        element: Element,
        span_start: Option<u32>,
        span_end: Option<u32>,
    ) -> Self {
        Self {
            symbol: AtomSymbol::Element(element),
            aromatic: Some(true),
            span: span_start.zip(span_end).map(|(s, e)| Span::bytes(s, e)),
            implicit_h: true,
            ..Default::default()
        }
    }
}

/// Atom symbol
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum AtomSymbol {
    Element(Element),
    NamedIsotope(NamedIsotope),
    Query(QueryAtom),
    Variable(Variable),
    // TODO: Add internal structure
    Pseudoatom(String),
    #[default]
    Unknown,
}

/// Variable atom
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Variable {}

/// Extended query atom types (superset of MOL + SMARTS)
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
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

/// Bond IR
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Bond {
    pub start_atom: u32,
    pub end_atom: u32,
    pub symbol: BondSymbol,
    pub ring: Option<u32>,
    pub stereo: Option<BondStereo>,
    pub direction: Option<BondDir>,

    // Metadata
    pub span: Option<Span>,
}

impl Bond {
    pub fn from_order(start_atom: u32, end_atom: u32, order: BondOrder) -> Self {
        Self {
            start_atom,
            end_atom,
            symbol: BondSymbol::Bond(order),
            ..Default::default()
        }
    }

    pub fn from_order_with_span(
        start_atom: u32,
        end_atom: u32,
        order: BondOrder,
        span_start: Option<u32>,
        span_end: Option<u32>,
    ) -> Self {
        Self {
            start_atom,
            end_atom,
            symbol: BondSymbol::Bond(order),
            span: span_start.zip(span_end).map(|(s, e)| Span::bytes(s, e)),
            ..Default::default()
        }
    }
}

/// Bond symbol
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum BondSymbol {
    Bond(BondOrder),
    Query(QueryBond),
    #[default]
    Unknown,
}

/// Discrete bond order
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum BondOrder {
    Zero,
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
            BondOrder::Zero => ".",
            BondOrder::Single => "-",
            BondOrder::Double => "=",
            BondOrder::Triple => "#",
            BondOrder::Quadruple => "$",
            BondOrder::Aromatic => ":",
            BondOrder::Unknown => "?",
        }
    }
}

/// Query bond
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum QueryBond {
    SingleOrDouble,
    SingleOrAromatic,
    DoubleOrAromatic,
    Any,
    #[default]
    Unknown,
}

/// Bond direction/wedging information
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum BondDir {
    Up,
    Down,
    Either,
    #[default]
    Unknown,
}

/// Double-bond stereochemistry (E/Z) annotation in IR
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum BondStereo {
    Cis,
    Trans,
    #[default]
    Either,
}

/// Chirality
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
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

/// Ring
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Ring {
    pub ring_idx: u32,
    pub start_atom: Option<u32>,
    pub end_atom: Option<u32>,
    pub open_start: Option<u32>,
    pub open_end: Option<u32>,
    pub close_start: Option<u32>,
    pub close_end: Option<u32>,
}

/// Fragment
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Fragment {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
}

/// Link
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Link {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
}

/// Property
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Property {
    pub name: String,
    pub value: String,
}

/// 3D coordinates
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
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
