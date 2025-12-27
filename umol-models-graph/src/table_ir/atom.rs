//! Atom IR for TableIR.

use std::collections::HashMap;

use umol_data::{Element, NamedIsotope};

use super::error::ConversionError;
use super::rgroup::RGroup;
use crate::span::Span;

/// Basic Atom IR
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Atom {
    pub element: Element,
    pub charge: Option<i8>,
    pub isotope_mass: Option<u32>,
    pub hydrogens: Option<u8>,
    pub implicit_h: bool,
    pub unpaired_e: Option<u8>,
    pub aromatic: Option<bool>,
    pub chirality: Option<Chirality>,
    pub class: Option<u32>,
    pub span: Option<Span>,
}

impl Atom {
    /// Create new aliphatic atom (aromatic flag false)
    pub fn aliphatic_atom(element: Element) -> Self {
        Self {
            element,
            charge: None,
            isotope_mass: None,
            hydrogens: None,
            implicit_h: true,
            unpaired_e: None,
            aromatic: Some(false),
            chirality: None,
            class: None,
            span: None,
        }
    }

    /// Create new aliphatic atom including span
    pub fn aliphatic_atom_with_span(element: Element, span: Span) -> Self {
        Self {
            element,
            charge: None,
            isotope_mass: None,
            hydrogens: None,
            implicit_h: true,
            unpaired_e: None,
            aromatic: Some(false),
            chirality: None,
            class: None,
            span: Some(span),
        }
    }

    /// Create new aromatic atom (aromatic flag true)
    pub fn aromatic_atom(element: Element) -> Self {
        Self {
            element,
            charge: None,
            isotope_mass: None,
            hydrogens: None,
            implicit_h: true,
            unpaired_e: None,
            aromatic: Some(true),
            chirality: None,
            class: None,
            span: None,
        }
    }

    /// Create new aromatic atom including span
    pub fn aromatic_atom_with_span(element: Element, span: Span) -> Self {
        Self {
            element,
            charge: None,
            isotope_mass: None,
            hydrogens: None,
            implicit_h: true,
            unpaired_e: None,
            aromatic: Some(true),
            chirality: None,
            class: None,
            span: Some(span),
        }
    }
}

/// Atom symbol (superset of CTFile and SMILES)
#[derive(Clone, Debug, PartialEq)]
pub enum AtomSymbol {
    Element(Element),
    NamedIsotope(NamedIsotope),
    GenericAtom(GenericAtom),
    AtomList(AtomList),
    RGroup(RGroup),
    Pseudoatom(String),
    LonePair,
}

impl AtomSymbol {
    /// Returns true if this is an extended atom structure (not a simple element or isotope)
    pub fn is_extended(&self) -> bool {
        !matches!(self, AtomSymbol::Element(_) | AtomSymbol::NamedIsotope(_))
    }
}

/// Generic atom types (SMILES and CXSMILES)
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GenericAtom {
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

impl GenericAtom {
    pub fn symbol(&self) -> &str {
        match self {
            GenericAtom::Any => "*",
            GenericAtom::Heavy => "A",
            GenericAtom::Heteroatom => "Q",
            GenericAtom::Halogen => "X",
            GenericAtom::Metal => "M",
            GenericAtom::HeavyOrH => "AH",
            GenericAtom::HeteroatomOrH => "QH",
            GenericAtom::HalogenOrH => "XH",
            GenericAtom::MetalOrH => "MH",
        }
    }

    pub fn from_symbol_bytes(s: &[u8]) -> Option<GenericAtom> {
        match s {
            b"*" => Some(GenericAtom::Any),
            b"A" => Some(GenericAtom::Heavy),
            b"Q" => Some(GenericAtom::Heteroatom),
            b"X" => Some(GenericAtom::Halogen),
            b"M" => Some(GenericAtom::Metal),
            b"AH" => Some(GenericAtom::HeavyOrH),
            b"QH" => Some(GenericAtom::HeteroatomOrH),
            b"XH" => Some(GenericAtom::HalogenOrH),
            b"MH" => Some(GenericAtom::MetalOrH),
            _ => None,
        }
    }

    pub fn from_symbol(s: &str) -> Option<GenericAtom> {
        Self::from_symbol_bytes(s.as_bytes())
    }
}

/// Atom list (inclusion or exclusion list of elements)
#[derive(Debug, Clone, PartialEq)]
pub struct AtomList {
    pub elements: Vec<Element>,
    pub exclusion: bool,
}

impl AtomList {
    pub fn empty() -> Self {
        Self {
            elements: Vec::new(),
            exclusion: false,
        }
    }
}

/// Chirality
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Chirality {
    Clockwise,
    CounterClockwise,
    Tetrahedral { arr: u32 },
    Allenal { arr: u32 },
    SquarePlanar { arr: u32 },
    TrigonalBipyramidal { arr: u32 },
    Octahedral { arr: u32 },
}

/// Extended atom IR
/// Temporary container for atomic features of generalized molecules.
/// Includes extended atom fields from CTFile and SMILES parsers extensions.
// TODO: Split into multiple semantically defined structures.
#[derive(Clone, Debug, PartialEq)]
pub struct ExtendedAtom {
    pub symbol: AtomSymbol,
    pub charge: Option<i8>,
    pub isotope_mass: Option<u32>,
    pub hydrogens: Option<u8>,
    pub implicit_h: bool,
    pub unpaired_e: Option<u8>,
    pub aromatic: Option<bool>,
    pub chirality: Option<Chirality>,
    pub class: Option<u32>,
    pub span: Option<Span>,
    pub stereo_parity: Option<AtomStereoParity>,
    pub stereo_care: Option<AtomStereoCare>,
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
    pub properties: HashMap<String, String>,
}
// End of TODO

impl ExtendedAtom {
    pub fn from_element(element: Element) -> Self {
        Self {
            symbol: AtomSymbol::Element(element),
            charge: None,
            isotope_mass: None,
            hydrogens: None,
            implicit_h: true,
            unpaired_e: None,
            aromatic: None,
            chirality: None,
            class: None,
            span: None,
            stereo_parity: None,
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

    /// Check if this atom has extended features that would be lost in conversion to basic Atom.
    pub fn has_extended_features(&self) -> bool {
        self.stereo_parity.is_some()
            || self.stereo_care.is_some()
            || self.valence.is_some()
            || self.atom_map_num.is_some()
            || self.inversion_retention.is_some()
            || self.exact_change.is_some()
            || self.attachment_point.is_some()
            || self.attachment_order.is_some()
            || self.ring_bond_count.is_some()
            || self.substitution_count.is_some()
            || self.unsaturated.is_some()
            || self.link_atom.is_some()
            || !self.properties.is_empty()
            || self.symbol.is_extended()
    }
}

impl From<Atom> for ExtendedAtom {
    fn from(atom: Atom) -> Self {
        Self {
            symbol: AtomSymbol::Element(atom.element),
            charge: atom.charge,
            isotope_mass: atom.isotope_mass,
            hydrogens: atom.hydrogens,
            implicit_h: atom.implicit_h,
            unpaired_e: atom.unpaired_e,
            aromatic: atom.aromatic,
            chirality: atom.chirality,
            class: atom.class,
            span: atom.span,
            stereo_parity: None,
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

impl TryFrom<ExtendedAtom> for Atom {
    type Error = ConversionError;

    fn try_from(extended: ExtendedAtom) -> Result<Self, Self::Error> {
        // Extract element from symbol (must be Element or NamedIsotope)
        let element = match extended.symbol {
            AtomSymbol::Element(e) => e,
            AtomSymbol::NamedIsotope(ni) => ni.element(),
            _ => {
                return Err(ConversionError::HasExtendedFeatures);
            }
        };

        // Check for extended features
        if extended.has_extended_features() {
            return Err(ConversionError::HasExtendedFeatures);
        }

        Ok(Self {
            element,
            charge: extended.charge,
            isotope_mass: extended.isotope_mass,
            hydrogens: extended.hydrogens,
            implicit_h: extended.implicit_h,
            unpaired_e: extended.unpaired_e,
            aromatic: extended.aromatic,
            chirality: extended.chirality,
            class: extended.class,
            span: extended.span,
        })
    }
}

/// Atom stereo parity
/// CTFile extension
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AtomStereoParity {
    Odd,  // Clockwise / R
    Even, // Counter-Clockwise / S
    Either,
}

/// Atom stereo care
/// CTFile extension
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AtomStereoCare {
    Care, // Stereo should be considered
}

/// Atom inversion/retention
/// CTFile extension
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AtomInversionRetention {
    Inverted,
    Retained,
}

/// Atom exact change flag
/// CTFile extension
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AtomExactChange {
    Match,
}

/// Attachment point type for R-groups
/// CTFile extension
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AttachmentPointType {
    First,
    Second,
    Both,
}

/// Ring bond count query
/// CTFile extension
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RingBondCount {
    AsDrawn,
    NoRingBonds,
    R2,
    R3,
    R4Plus,
}

/// Substitution count query
/// CTFile extension
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SubstitutionCount {
    AsDrawn,
    NoSubstitution,
    S1,
    S2,
    S3,
    S4,
    S5,
    S6,
    S6Plus,
    S7,
    S8,
    S9,
    S10,
}

/// Unsaturated atom flag (query feature)
/// CTFile extension
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UnsaturatedAtom;

/// Link atom (for polymers)
/// CTFile extension
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinkAtom {
    pub repeat_count: u8,
    pub subs_index1: usize,
    pub subs_index2: Option<usize>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atom_aliphatic() {
        let atom = Atom::aliphatic_atom(Element::C);
        assert_eq!(atom.element, Element::C);
        assert_eq!(atom.aromatic, Some(false));
        assert!(atom.implicit_h);
    }

    #[test]
    fn test_atom_aromatic() {
        let atom = Atom::aromatic_atom(Element::C);
        assert_eq!(atom.element, Element::C);
        assert_eq!(atom.aromatic, Some(true));
        assert!(atom.implicit_h);
    }

    #[test]
    fn test_extended_atom_from_element() {
        let extended = ExtendedAtom::from_element(Element::N);
        assert_eq!(extended.symbol, AtomSymbol::Element(Element::N));
        assert!(extended.implicit_h);
        assert!(extended.aromatic.is_none());
        assert!(extended.properties.is_empty());
        assert!(!extended.has_extended_features());
    }

    #[test]
    fn test_from_atom_to_extended_atom() {
        let atom = Atom {
            element: Element::C,
            charge: Some(1),
            isotope_mass: Some(13),
            hydrogens: Some(3),
            implicit_h: false,
            unpaired_e: Some(1),
            aromatic: Some(true),
            chirality: Some(Chirality::Clockwise),
            class: Some(5),
            span: None,
        };

        let extended: ExtendedAtom = atom.into();
        assert_eq!(extended.symbol, AtomSymbol::Element(Element::C));
        assert_eq!(extended.charge, Some(1));
        assert_eq!(extended.isotope_mass, Some(13));
        assert_eq!(extended.hydrogens, Some(3));
        assert!(!extended.implicit_h);
        assert_eq!(extended.unpaired_e, Some(1));
        assert_eq!(extended.aromatic, Some(true));
        assert_eq!(extended.chirality, Some(Chirality::Clockwise));
        assert_eq!(extended.class, Some(5));
        assert!(!extended.has_extended_features());
    }

    #[test]
    fn test_try_from_extended_atom_to_atom_element() {
        let extended = ExtendedAtom {
            symbol: AtomSymbol::Element(Element::O),
            charge: Some(-1),
            isotope_mass: None,
            hydrogens: Some(1),
            implicit_h: true,
            unpaired_e: None,
            aromatic: Some(false),
            chirality: None,
            class: None,
            span: None,
            stereo_parity: None,
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
        };

        let atom: Atom = extended.try_into().unwrap();
        assert_eq!(atom.element, Element::O);
        assert_eq!(atom.charge, Some(-1));
        assert_eq!(atom.hydrogens, Some(1));
        assert_eq!(atom.aromatic, Some(false));
    }

    #[test]
    fn test_try_from_extended_atom_to_atom_named_isotope() {
        let extended = ExtendedAtom {
            symbol: AtomSymbol::NamedIsotope(NamedIsotope::D),
            charge: None,
            isotope_mass: None,
            hydrogens: None,
            implicit_h: true,
            unpaired_e: None,
            aromatic: None,
            chirality: None,
            class: None,
            span: None,
            stereo_parity: None,
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
        };

        let atom: Atom = extended.try_into().unwrap();
        assert_eq!(atom.element, Element::H);
    }

    #[test]
    fn test_try_from_extended_atom_to_atom_invalid() {
        let extended = ExtendedAtom {
            symbol: AtomSymbol::GenericAtom(GenericAtom::Any),
            charge: None,
            isotope_mass: None,
            hydrogens: None,
            implicit_h: true,
            unpaired_e: None,
            aromatic: None,
            chirality: None,
            class: None,
            span: None,
            stereo_parity: None,
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
        };

        let result: Result<Atom, _> = extended.try_into();
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            ConversionError::HasExtendedFeatures
        ));
    }

    #[test]
    fn test_has_extended_features_basic() {
        let extended = ExtendedAtom {
            symbol: AtomSymbol::Element(Element::C),
            charge: Some(1),
            isotope_mass: Some(13),
            hydrogens: Some(2),
            implicit_h: false,
            unpaired_e: Some(1),
            aromatic: Some(true),
            chirality: Some(Chirality::Clockwise),
            class: Some(5),
            span: None,
            stereo_parity: None,
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
        };

        assert!(!extended.has_extended_features());
    }

    #[test]
    fn test_has_extended_features_extended() {
        let extended = ExtendedAtom {
            symbol: AtomSymbol::Element(Element::C),
            charge: None,
            isotope_mass: None,
            hydrogens: None,
            implicit_h: true,
            unpaired_e: None,
            aromatic: None,
            chirality: None,
            class: None,
            span: None,
            stereo_parity: Some(AtomStereoParity::Even),
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
        };

        assert!(extended.has_extended_features());
    }

    #[test]
    fn test_roundtrip_atom_to_extended_atom_to_atom() {
        let atom = Atom {
            element: Element::N,
            charge: Some(-1),
            isotope_mass: Some(15),
            hydrogens: Some(2),
            implicit_h: false,
            unpaired_e: None,
            aromatic: Some(false),
            chirality: Some(Chirality::CounterClockwise),
            class: Some(10),
            span: None,
        };

        let extended: ExtendedAtom = atom.into();
        let atom2: Atom = extended.try_into().unwrap();

        assert_eq!(atom.element, atom2.element);
        assert_eq!(atom.charge, atom2.charge);
        assert_eq!(atom.isotope_mass, atom2.isotope_mass);
        assert_eq!(atom.hydrogens, atom2.hydrogens);
        assert_eq!(atom.implicit_h, atom2.implicit_h);
        assert_eq!(atom.unpaired_e, atom2.unpaired_e);
        assert_eq!(atom.aromatic, atom2.aromatic);
        assert_eq!(atom.chirality, atom2.chirality);
        assert_eq!(atom.class, atom2.class);
        assert_eq!(atom.span, atom2.span);
    }
}
