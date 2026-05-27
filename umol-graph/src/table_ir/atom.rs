//! Atom types for TableIR.

use std::collections::HashMap;

use umol_shared::element::Element;
use umol_shared::isotope::NamedIsotope;
use umol_shared::spin::SpinMultiplicity;

use super::error::ConversionError;
use super::rgroup::RGroup;
use crate::span::Span;

/// Basic Atom IR
#[derive(Clone, Debug, PartialEq)]
pub struct Atom {
    pub element: Element,
    pub isotope_mass: Option<u32>,
    pub charge: Option<i8>,
    pub implicit_hydrogens: Option<u8>,
    pub valence: Option<u8>,
    pub lone_pairs: Option<u8>,
    pub unpaired_electrons: Option<u8>,
    pub multiplicity: Option<SpinMultiplicity>,
    pub aromatic: Option<bool>,
    pub chirality: Option<Chirality>,
    pub class: Option<u32>,
    pub label: Option<String>,
    pub value: Option<String>,
    pub span: Option<Span>,
}

impl Atom {
    /// Create new atom from element (default for MOL/CTFile)
    pub fn from_element(element: Element) -> Self {
        Self {
            element,
            isotope_mass: None,
            charge: None,
            implicit_hydrogens: None,
            valence: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic: None,
            chirality: None,
            class: None,
            label: None,
            value: None,
            span: None,
        }
    }
    /// Create new aliphatic atom (aromatic flag false, infer implicit hydrogens)
    pub fn aliphatic_atom(element: Element) -> Self {
        Self {
            element,
            isotope_mass: None,
            charge: None,
            implicit_hydrogens: None,
            valence: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic: Some(false),
            chirality: None,
            class: None,
            label: None,
            value: None,
            span: None,
        }
    }

    /// Create new aliphatic atom including span
    pub fn aliphatic_atom_with_span(element: Element, span: Span) -> Self {
        Self {
            element,
            isotope_mass: None,
            charge: None,
            implicit_hydrogens: None,
            valence: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic: Some(false),
            chirality: None,
            class: None,
            label: None,
            value: None,
            span: Some(span),
        }
    }

    /// Create new aromatic atom (aromatic flag true, infer implicit hydrogens)
    pub fn aromatic_atom(element: Element) -> Self {
        Self {
            element,
            isotope_mass: None,
            charge: None,
            implicit_hydrogens: None,
            valence: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic: Some(true),
            chirality: None,
            class: None,
            label: None,
            value: None,
            span: None,
        }
    }

    /// Create new aromatic atom including span
    pub fn aromatic_atom_with_span(element: Element, span: Span) -> Self {
        Self {
            element,
            isotope_mass: None,
            charge: None,
            implicit_hydrogens: None,
            valence: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic: Some(true),
            chirality: None,
            class: None,
            label: None,
            value: None,
            span: Some(span),
        }
    }
}

/// Atom symbol (superset of CTFile and SMILES)
#[derive(Clone, Debug, PartialEq)]
pub enum AtomSymbol {
    Element(Element),
    NamedIsotope(NamedIsotope),
    WildcardAtom(WildcardAtom),
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

/// Wildcard atom types (CTFile, SMILES, CXSMILES)
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum WildcardAtom {
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

impl WildcardAtom {
    pub fn symbol(&self) -> &str {
        match self {
            WildcardAtom::Any => "*",
            WildcardAtom::Heavy => "A",
            WildcardAtom::Heteroatom => "Q",
            WildcardAtom::Halogen => "X",
            WildcardAtom::Metal => "M",
            WildcardAtom::HeavyOrH => "AH",
            WildcardAtom::HeteroatomOrH => "QH",
            WildcardAtom::HalogenOrH => "XH",
            WildcardAtom::MetalOrH => "MH",
        }
    }

    pub fn from_symbol_bytes(s: &[u8]) -> Option<WildcardAtom> {
        match s {
            b"*" => Some(WildcardAtom::Any),
            b"A" => Some(WildcardAtom::Heavy),
            b"Q" => Some(WildcardAtom::Heteroatom),
            b"X" => Some(WildcardAtom::Halogen),
            b"M" => Some(WildcardAtom::Metal),
            b"AH" => Some(WildcardAtom::HeavyOrH),
            b"QH" => Some(WildcardAtom::HeteroatomOrH),
            b"XH" => Some(WildcardAtom::HalogenOrH),
            b"MH" => Some(WildcardAtom::MetalOrH),
            _ => None,
        }
    }

    pub fn from_symbol(s: &str) -> Option<WildcardAtom> {
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
    Unspecified,
    Tetrahedral { arr: u32 },
    Allenal { arr: u32 },
    SquarePlanar { arr: u32 },
    TrigonalBipyramidal { arr: u32 },
    Octahedral { arr: u32 },
}

/// Bicyclic bridge stereo (CXSMILES THB/TLB/TEB).
/// Describes the configuration at a bridgehead in a bicyclic system.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BicycloStereoData {
    pub ligand_atom: u32,
    pub connection_atom: u32,
    pub lower_bridge_atoms: Vec<u32>,
    pub higher_bridge_atoms: Vec<u32>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BicycloStereo {
    TowardsHigherBridge(BicycloStereoData),
    TowardsLowerBridge(BicycloStereoData),
    TowardsEitherBridge(BicycloStereoData),
}

/// Extended atom IR
/// Temporary container for atomic features of generalized molecules.
/// Includes extended atom fields from CTFile and SMILES parsers extensions.
// TODO: Split into multiple semantically defined structures.
#[derive(Clone, Debug, PartialEq)]
pub struct ExtendedAtom {
    pub symbol: AtomSymbol,
    pub isotope_mass: Option<u32>,
    pub charge: Option<i8>,
    pub implicit_hydrogens: Option<u8>,
    pub valence: Option<u8>,
    pub lone_pairs: Option<u8>,
    pub unpaired_electrons: Option<u8>,
    pub multiplicity: Option<SpinMultiplicity>,
    pub aromatic: Option<bool>,
    pub chirality: Option<Chirality>,
    pub class: Option<u32>,
    pub label: Option<String>,
    pub value: Option<String>,
    pub pattern: Option<String>,
    pub stereo_care: Option<AtomStereoCare>,
    pub inversion_retention: Option<AtomInversionRetention>,
    pub exact_change: Option<AtomExactChange>,
    pub attachment_point: Option<AttachmentPointType>,
    pub attachment_order: Option<Vec<(u32, u8)>>,
    pub ligand_order: Option<Vec<(u32, u8)>>,
    pub ring_bond_count: Option<RingBondCount>,
    pub substitution_count: Option<SubstitutionCount>,
    pub unsaturated: Option<UnsaturatedAtom>,
    pub link_atom: Option<LinkAtom>,
    pub properties: HashMap<String, String>,
    pub span: Option<Span>,
}
// End of TODO

impl ExtendedAtom {
    /// Create new extended atom from symbol (default for MOL/CTFile, no implicit hydrogens)
    pub fn from_atom_symbol(symbol: AtomSymbol) -> Self {
        Self {
            symbol,
            isotope_mass: None,
            charge: None,
            implicit_hydrogens: None,
            valence: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic: None,
            chirality: None,
            class: None,
            label: None,
            value: None,
            pattern: None,
            stereo_care: None,
            inversion_retention: None,
            exact_change: None,
            attachment_point: None,
            attachment_order: None,
            ligand_order: None,
            ring_bond_count: None,
            substitution_count: None,
            unsaturated: None,
            link_atom: None,
            properties: HashMap::new(),
            span: None,
        }
    }

    pub fn from_element(element: Element) -> Self {
        Self::from_atom_symbol(AtomSymbol::Element(element))
    }

    pub fn from_named_isotope(isotope: NamedIsotope) -> Self {
        let element = isotope.element();
        let mass_number = isotope.mass_number();
        let mut atom = Self::from_element(element);
        atom.isotope_mass = Some(mass_number);
        atom
    }

    /// Check if this atom has extended features that would be lost in conversion to basic Atom.
    pub fn has_extended_features(&self) -> bool {
        self.symbol.is_extended()
            || self.pattern.is_some()
            || self.stereo_care.is_some()
            || self.inversion_retention.is_some()
            || self.exact_change.is_some()
            || self.attachment_point.is_some()
            || self.attachment_order.is_some()
            || self.ligand_order.is_some()
            || self.ring_bond_count.is_some()
            || self.substitution_count.is_some()
            || self.unsaturated.is_some()
            || self.link_atom.is_some()
            || has_extended_properties(&self.properties)
    }
}

/// Check if properties HashMap contains extended properties (excluding basic alias/value).
fn has_extended_properties(properties: &HashMap<String, String>) -> bool {
    properties
        .keys()
        .any(|k| k != "molFileAlias" && k != "molFileValue")
}

impl From<Atom> for ExtendedAtom {
    fn from(atom: Atom) -> Self {
        Self {
            symbol: AtomSymbol::Element(atom.element),
            isotope_mass: atom.isotope_mass,
            charge: atom.charge,
            implicit_hydrogens: atom.implicit_hydrogens,
            valence: atom.valence,
            lone_pairs: atom.lone_pairs,
            unpaired_electrons: atom.unpaired_electrons,
            multiplicity: atom.multiplicity,
            aromatic: atom.aromatic,
            chirality: atom.chirality,
            class: atom.class,
            label: atom.label,
            value: atom.value,
            pattern: None,
            stereo_care: None,
            inversion_retention: None,
            exact_change: None,
            attachment_point: None,
            attachment_order: None,
            ligand_order: None,
            ring_bond_count: None,
            substitution_count: None,
            unsaturated: None,
            link_atom: None,
            properties: HashMap::new(),
            span: atom.span,
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
            isotope_mass: extended.isotope_mass,
            charge: extended.charge,
            implicit_hydrogens: extended.implicit_hydrogens,
            valence: extended.valence,
            lone_pairs: extended.lone_pairs,
            unpaired_electrons: extended.unpaired_electrons,
            multiplicity: extended.multiplicity,
            aromatic: extended.aromatic,
            chirality: extended.chirality,
            class: extended.class,
            label: extended.label,
            value: extended.value,
            span: extended.span,
        })
    }
}

/// Atom stereo parity
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
    pub min_repeat: u8,
    pub repeat_count: u8,
    pub subs_index1: u32,
    pub subs_index2: Option<u32>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atom_aliphatic() {
        let atom = Atom::aliphatic_atom(Element::C);
        assert_eq!(atom.element, Element::C);
        assert_eq!(atom.aromatic, Some(false));
        assert_eq!(atom.implicit_hydrogens, None);
    }

    #[test]
    fn test_atom_aromatic() {
        let atom = Atom::aromatic_atom(Element::C);
        assert_eq!(atom.element, Element::C);
        assert_eq!(atom.aromatic, Some(true));
        assert_eq!(atom.implicit_hydrogens, None);
    }

    #[test]
    fn test_extended_atom_from_element() {
        let extended = ExtendedAtom::from_element(Element::N);
        assert_eq!(extended.symbol, AtomSymbol::Element(Element::N));
        assert_eq!(extended.implicit_hydrogens, None);
        assert!(extended.aromatic.is_none());
        assert!(extended.properties.is_empty());
        assert!(!extended.has_extended_features());
    }

    #[test]
    fn test_from_atom_to_extended_atom() {
        let atom = Atom {
            element: Element::C,
            isotope_mass: Some(13),
            charge: Some(1),
            implicit_hydrogens: Some(3),
            valence: Some(3),
            lone_pairs: None,
            unpaired_electrons: Some(1),
            multiplicity: Some(SpinMultiplicity::Doublet),
            aromatic: Some(true),
            chirality: Some(Chirality::Clockwise),
            class: Some(5),
            label: None,
            value: None,
            span: None,
        };

        let extended: ExtendedAtom = atom.into();
        assert_eq!(extended.symbol, AtomSymbol::Element(Element::C));
        assert_eq!(extended.isotope_mass, Some(13));
        assert_eq!(extended.charge, Some(1));
        assert_eq!(
            extended.implicit_hydrogens,
            Some(3)
        );
        assert_eq!(extended.valence, Some(3));
        assert_eq!(extended.lone_pairs, None);
        assert_eq!(extended.unpaired_electrons, Some(1));
        assert_eq!(extended.multiplicity, Some(SpinMultiplicity::Doublet));
        assert_eq!(extended.aromatic, Some(true));
        assert_eq!(extended.chirality, Some(Chirality::Clockwise));
        assert_eq!(extended.class, Some(5));
        assert!(!extended.has_extended_features());
    }

    #[test]
    fn test_try_from_extended_atom_to_atom_element() {
        let extended = ExtendedAtom {
            symbol: AtomSymbol::Element(Element::O),
            isotope_mass: None,
            charge: Some(-1),
            implicit_hydrogens: None,
            valence: Some(1),
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic: Some(false),
            chirality: None,
            class: None,
            label: None,
            value: None,
            pattern: None,
            stereo_care: None,
            inversion_retention: None,
            exact_change: None,
            attachment_point: None,
            attachment_order: None,
            ligand_order: None,
            ring_bond_count: None,
            substitution_count: None,
            unsaturated: None,
            link_atom: None,
            properties: HashMap::new(),
            span: None,
        };

        let atom: Atom = extended.try_into().unwrap();
        assert_eq!(atom.element, Element::O);
        assert_eq!(atom.isotope_mass, None);
        assert_eq!(atom.charge, Some(-1));
        assert_eq!(atom.implicit_hydrogens, None);
        assert_eq!(atom.valence, Some(1));
        assert_eq!(atom.lone_pairs, None);
        assert!(atom.unpaired_electrons.is_none());
        assert!(atom.multiplicity.is_none());
        assert_eq!(atom.aromatic, Some(false));
    }

    #[test]
    fn test_try_from_extended_atom_to_atom_named_isotope() {
        let named_isotope = NamedIsotope::D;
        let extended = ExtendedAtom::from_named_isotope(named_isotope);
        let atom: Atom = extended.try_into().unwrap();
        assert_eq!(atom.element, Element::H);
        assert_eq!(atom.isotope_mass, Some(2));
        assert_eq!(atom.charge, None);
        assert_eq!(atom.implicit_hydrogens, None);
        assert_eq!(atom.valence, None);
        assert_eq!(atom.lone_pairs, None);
        assert!(atom.unpaired_electrons.is_none());
        assert_eq!(atom.aromatic, None);
    }

    #[test]
    fn test_try_from_extended_atom_to_atom_invalid() {
        let extended = ExtendedAtom {
            symbol: AtomSymbol::WildcardAtom(WildcardAtom::Any),
            isotope_mass: None,
            charge: None,
            implicit_hydrogens: None,
            valence: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic: None,
            chirality: None,
            class: None,
            label: None,
            value: None,
            pattern: None,
            stereo_care: None,
            inversion_retention: None,
            exact_change: None,
            attachment_point: None,
            attachment_order: None,
            ligand_order: None,
            ring_bond_count: None,
            substitution_count: None,
            unsaturated: None,
            link_atom: None,
            properties: HashMap::new(),
            span: None,
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
            isotope_mass: Some(13),
            charge: Some(1),
            implicit_hydrogens: Some(2),
            valence: None,
            lone_pairs: None,
            unpaired_electrons: Some(1),
            multiplicity: Some(SpinMultiplicity::Doublet),
            aromatic: Some(true),
            chirality: Some(Chirality::Clockwise),
            class: Some(5),
            label: None,
            value: None,
            pattern: None,
            stereo_care: None,
            inversion_retention: None,
            exact_change: None,
            attachment_point: None,
            attachment_order: None,
            ligand_order: None,
            ring_bond_count: None,
            substitution_count: None,
            unsaturated: None,
            link_atom: None,
            properties: HashMap::new(),
            span: None,
        };

        assert!(!extended.has_extended_features());
    }

    #[test]
    fn test_has_extended_features_extended() {
        let extended = ExtendedAtom {
            symbol: AtomSymbol::Element(Element::C),
            isotope_mass: None,
            charge: None,
            implicit_hydrogens: None,
            valence: None,
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic: None,
            chirality: None,
            class: None,
            label: None,
            value: None,
            pattern: None,
            stereo_care: Some(AtomStereoCare::Care),
            inversion_retention: None,
            exact_change: None,
            attachment_point: None,
            attachment_order: None,
            ligand_order: None,
            ring_bond_count: None,
            substitution_count: None,
            unsaturated: None,
            link_atom: None,
            properties: HashMap::new(),
            span: None,
        };

        assert!(extended.has_extended_features());
    }

    #[test]
    fn test_roundtrip_atom_to_extended_atom_to_atom() {
        let atom = Atom {
            element: Element::N,
            isotope_mass: Some(15),
            charge: Some(-1),
            implicit_hydrogens: Some(2),
            valence: Some(3),
            lone_pairs: None,
            unpaired_electrons: None,
            multiplicity: None,
            aromatic: Some(false),
            chirality: Some(Chirality::CounterClockwise),
            span: None,
            label: None,
            value: None,
            class: Some(10),
        };

        let extended: ExtendedAtom = atom.clone().into();
        let atom2: Atom = extended.try_into().unwrap();

        assert_eq!(atom.element, atom2.element);
        assert_eq!(atom.isotope_mass, atom2.isotope_mass);
        assert_eq!(atom.charge, atom2.charge);
        assert_eq!(atom.implicit_hydrogens, atom2.implicit_hydrogens);
        assert_eq!(atom.valence, atom2.valence);
        assert_eq!(atom.lone_pairs, atom2.lone_pairs);
        assert_eq!(atom.unpaired_electrons, atom2.unpaired_electrons);
        assert_eq!(atom.multiplicity, atom2.multiplicity);
        assert_eq!(atom.aromatic, atom2.aromatic);
        assert_eq!(atom.chirality, atom2.chirality);
        assert_eq!(atom.class, atom2.class);
        assert_eq!(atom.label, atom2.label);
        assert_eq!(atom.value, atom2.value);
        assert_eq!(atom.span, atom2.span);
    }
}
