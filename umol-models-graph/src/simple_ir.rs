//! Simple IR for atom/bond-based molecular models.

use std::mem;

use serde::{Deserialize, Serialize};
use umol_data::{Element, NamedIsotope};

use crate::position::Point3D;
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
            span: Span::from_bytes_opt(span_start, span_end),
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
            span: Span::from_bytes_opt(span_start, span_end),
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
            span: Span::from_bytes_opt(span_start, span_end),
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
    pub open_span: Option<Span>,
    pub close_span: Option<Span>,
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

/// Atom data supplied by parsers before conversion into IR atoms.
pub struct AtomData {
    pub element: Element,
    pub isotope: Option<u32>,
    pub charge: Option<i32>,
    pub hydrogen_count: Option<u32>,
    pub class: Option<u32>,
    pub aromatic: bool,
    pub implicit_h: bool,
    pub chirality: Option<Chirality>,
    pub unknown_symbol: bool,
    pub span: Option<Span>,
}

/// Bond data supplied by parsers before conversion into IR bonds.
pub struct BondData {
    pub order: BondOrder,
    pub dir: Option<BondDir>,
    pub span: Option<Span>,
}

/// Builder used by tokenizers to incrementally assemble SIR molecules.
pub struct MoleculeBuilder {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    rings: Vec<Ring>,
    molecules: Vec<Molecule>,
}

impl MoleculeBuilder {
    pub fn with_capacity(approx_atoms: usize, approx_bonds: usize) -> Self {
        Self {
            atoms: Vec::with_capacity(approx_atoms),
            bonds: Vec::with_capacity(approx_bonds),
            rings: Vec::new(),
            molecules: Vec::new(),
        }
    }

    pub fn clear_reuse(&mut self) {
        self.atoms.clear();
        self.bonds.clear();
        self.rings.clear();
        self.molecules.clear();
    }

    #[inline]
    pub fn on_atom(&mut self, a: AtomData) -> u32 {
        let span = a.span;
        let atom = if a.unknown_symbol {
            Atom {
                symbol: AtomSymbol::Unknown,
                position: None,
                charge: a.charge,
                isotope: a.isotope,
                unpaired_e: None,
                hydrogens: a.hydrogen_count,
                implicit_h: a.implicit_h,
                aromatic: Some(a.aromatic),
                chirality: a.chirality,
                class: a.class,
                span,
            }
        } else {
            Atom {
                symbol: AtomSymbol::Element(a.element),
                position: None,
                isotope: a.isotope,
                unpaired_e: None,
                charge: a.charge,
                hydrogens: a.hydrogen_count,
                implicit_h: a.implicit_h,
                aromatic: Some(a.aromatic),
                chirality: a.chirality,
                class: a.class,
                span,
            }
        };

        self.atoms.push(atom);
        self.atoms.len() as u32
    }

    #[inline]
    pub fn on_bond(&mut self, start: u32, end: u32, b: BondData) {
        let span = b.span;
        let bond = Bond {
            start_atom: start,
            end_atom: end,
            symbol: BondSymbol::Bond(b.order),
            direction: b.dir,
            ring: None,
            stereo: None,
            span,
        };
        self.bonds.push(bond);
    }

    #[inline]
    pub fn on_atom_fast(
        &mut self,
        element: Element,
        aromatic: bool,
        span_start: Option<u32>,
        span_end: Option<u32>,
    ) -> u32 {
        let span = Span::from_bytes_opt(span_start, span_end);
        let atom = if aromatic {
            let mut a = Atom::from_aromatic_atom(element);
            a.span = span;
            a
        } else {
            let mut a = Atom::from_aliphatic_atom(element);
            a.span = span;
            a
        };
        self.atoms.push(atom);
        self.atoms.len() as u32
    }

    #[inline]
    pub fn on_bond_single_fast(
        &mut self,
        start_atom: u32,
        end_atom: u32,
        span_start: Option<u32>,
        span_end: Option<u32>,
    ) {
        let span = Span::from_bytes_opt(span_start, span_end);
        let mut bond = Bond::from_order(start_atom, end_atom, BondOrder::Single);
        bond.span = span;
        self.bonds.push(bond);
    }

    pub fn on_component_end(&mut self) {
        if self.atoms.is_empty() {
            return;
        }
        let mut mol = Molecule::default();
        mol.atoms = mem::take(&mut self.atoms);
        mol.bonds = mem::take(&mut self.bonds);
        mol.rings = mem::take(&mut self.rings);
        self.molecules.push(mol);
    }

    pub fn finish(&mut self) -> Vec<Molecule> {
        if !self.atoms.is_empty() || !self.bonds.is_empty() {
            self.on_component_end();
        }
        mem::take(&mut self.molecules)
    }

    #[inline]
    pub fn annotate_last_atom_span(&mut self, start: u32) {
        if let Some(a) = self.atoms.last_mut() {
            a.span = Some(match a.span {
                Some(span) => span.with_start(start),
                None => Span::bytes(start, start),
            });
        }
    }

    #[inline]
    pub fn annotate_last_bond_span(&mut self, start: u32) {
        if let Some(b) = self.bonds.last_mut() {
            b.span = Some(match b.span {
                Some(span) => span.with_start(start),
                None => Span::bytes(start, start),
            });
        }
    }

    #[inline]
    pub fn on_ring_open(
        &mut self,
        ring_idx: u32,
        start: Option<u32>,
        end: Option<u32>,
        atom_idx: Option<u32>,
    ) {
        self.rings.push(Ring {
            ring_idx,
            open_span: Span::from_bytes_opt(start, end),
            close_span: None,
            start_atom: atom_idx,
            end_atom: None,
        });
    }

    #[inline]
    pub fn on_ring_close(
        &mut self,
        ring_idx: u32,
        start: Option<u32>,
        end: Option<u32>,
        atom_idx: Option<u32>,
    ) {
        for ev in self.rings.iter_mut().rev() {
            if ev.ring_idx == ring_idx && ev.close_span.is_none() {
                ev.close_span = Span::from_bytes_opt(start, end);
                ev.end_atom = atom_idx;
                return;
            }
        }

        self.rings.push(Ring {
            ring_idx,
            open_span: None,
            close_span: Span::from_bytes_opt(start, end),
            start_atom: None,
            end_atom: atom_idx,
        });
    }
}
