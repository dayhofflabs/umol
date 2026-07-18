//! Atom types for TableIR.

use std::collections::HashMap;

use umol_chem::element::Element;
use umol_chem::isotope::NamedIsotope;
use umol_chem::spin::SpinMultiplicity;

use super::error::ConversionError;
use super::rgroup::RGroup;
use super::span::Span;

/// Basic Atom IR
#[derive(Clone, Debug, PartialEq)]
pub struct Atom {
    /// Concrete element, or `None` for the OpenSMILES `*` wildcard.
    pub element: Option<Element>,
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
            element: Some(element),
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

    /// Create an OpenSMILES `*` wildcard atom.
    pub fn wildcard() -> Self {
        Self {
            element: None,
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

    /// Create an OpenSMILES `*` wildcard atom including its source span.
    pub fn wildcard_with_span(span: Span) -> Self {
        Self {
            span: Some(span),
            ..Self::wildcard()
        }
    }

    /// Create new aliphatic atom (aromatic flag false, infer implicit hydrogens)
    pub fn aliphatic_atom(element: Element) -> Self {
        Self {
            element: Some(element),
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
            element: Some(element),
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
            element: Some(element),
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
            element: Some(element),
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
        !matches!(
            self.symbol,
            AtomSymbol::Element(_)
                | AtomSymbol::NamedIsotope(_)
                | AtomSymbol::WildcardAtom(WildcardAtom::Any)
        ) || self.pattern.is_some()
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
            symbol: match atom.element {
                Some(element) => AtomSymbol::Element(element),
                None => AtomSymbol::WildcardAtom(WildcardAtom::Any),
            },
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
        if extended.has_extended_features() {
            return Err(ConversionError::HasExtendedFeatures);
        }

        let element = match extended.symbol {
            AtomSymbol::Element(element) => Some(element),
            AtomSymbol::NamedIsotope(isotope) => Some(isotope.element()),
            AtomSymbol::WildcardAtom(WildcardAtom::Any) => None,
            _ => unreachable!("extended symbols were rejected above"),
        };

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
    use rstest::rstest;

    use super::*;

    #[rstest]
    fn test_atom_aliphatic_atom() {
        let atom = Atom::aliphatic_atom(Element::C);
        assert_eq!(atom.element, Some(Element::C));
        assert_eq!(atom.aromatic, Some(false));
        assert_eq!(atom.implicit_hydrogens, None);
    }

    #[rstest]
    fn test_atom_aromatic_atom() {
        let atom = Atom::aromatic_atom(Element::C);
        assert_eq!(atom.element, Some(Element::C));
        assert_eq!(atom.aromatic, Some(true));
        assert_eq!(atom.implicit_hydrogens, None);
    }

    #[rstest]
    fn test_atom_wildcard() {
        assert_eq!(
            Atom::wildcard(),
            Atom {
                element: None,
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
        );
    }

    #[rstest]
    fn test_atom_wildcard_with_span() {
        assert_eq!(
            Atom::wildcard_with_span(Span::bytes(2, 3)),
            Atom {
                element: None,
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
                span: Some(Span::bytes(2, 3)),
            }
        );
    }

    #[rstest]
    fn test_extended_atom_from_element() {
        let extended = ExtendedAtom::from_element(Element::N);
        assert_eq!(extended.symbol, AtomSymbol::Element(Element::N));
        assert_eq!(extended.implicit_hydrogens, None);
        assert!(extended.aromatic.is_none());
        assert!(extended.properties.is_empty());
        assert!(!extended.has_extended_features());
    }

    #[rstest]
    fn test_extended_atom_from_atom() {
        let atom = Atom {
            element: Some(Element::C),
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
        assert_eq!(extended.implicit_hydrogens, Some(3));
        assert_eq!(extended.valence, Some(3));
        assert_eq!(extended.lone_pairs, None);
        assert_eq!(extended.unpaired_electrons, Some(1));
        assert_eq!(extended.multiplicity, Some(SpinMultiplicity::Doublet));
        assert_eq!(extended.aromatic, Some(true));
        assert_eq!(extended.chirality, Some(Chirality::Clockwise));
        assert_eq!(extended.class, Some(5));
        assert!(!extended.has_extended_features());
    }

    #[rstest]
    fn test_atom_try_from_element() {
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
        assert_eq!(atom.element, Some(Element::O));
        assert_eq!(atom.isotope_mass, None);
        assert_eq!(atom.charge, Some(-1));
        assert_eq!(atom.implicit_hydrogens, None);
        assert_eq!(atom.valence, Some(1));
        assert_eq!(atom.lone_pairs, None);
        assert!(atom.unpaired_electrons.is_none());
        assert!(atom.multiplicity.is_none());
        assert_eq!(atom.aromatic, Some(false));
    }

    #[rstest]
    fn test_atom_try_from_named_isotope() {
        let named_isotope = NamedIsotope::D;
        let extended = ExtendedAtom::from_named_isotope(named_isotope);
        let atom: Atom = extended.try_into().unwrap();
        assert_eq!(atom.element, Some(Element::H));
        assert_eq!(atom.isotope_mass, Some(2));
        assert_eq!(atom.charge, None);
        assert_eq!(atom.implicit_hydrogens, None);
        assert_eq!(atom.valence, None);
        assert_eq!(atom.lone_pairs, None);
        assert!(atom.unpaired_electrons.is_none());
        assert_eq!(atom.aromatic, None);
    }

    #[rstest]
    fn test_atom_try_from_wildcard() {
        let atom = Atom::try_from(ExtendedAtom::from_atom_symbol(AtomSymbol::WildcardAtom(
            WildcardAtom::Any,
        )))
        .unwrap();
        assert_eq!(atom, Atom::wildcard());
    }

    #[rstest]
    #[case::heavy(WildcardAtom::Heavy)]
    #[case::heteroatom(WildcardAtom::Heteroatom)]
    #[case::halogen(WildcardAtom::Halogen)]
    #[case::metal(WildcardAtom::Metal)]
    #[case::heavy_or_h(WildcardAtom::HeavyOrH)]
    #[case::heteroatom_or_h(WildcardAtom::HeteroatomOrH)]
    #[case::halogen_or_h(WildcardAtom::HalogenOrH)]
    #[case::metal_or_h(WildcardAtom::MetalOrH)]
    fn test_atom_try_from_error(#[case] wildcard: WildcardAtom) {
        let extended = ExtendedAtom::from_atom_symbol(AtomSymbol::WildcardAtom(wildcard));
        assert_eq!(
            Atom::try_from(extended),
            Err(ConversionError::HasExtendedFeatures)
        );
    }

    #[rstest]
    fn test_extended_atom_has_extended_features_basic_fields() {
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

    #[rstest]
    #[case::element(AtomSymbol::Element(Element::C), false)]
    #[case::named_isotope(AtomSymbol::NamedIsotope(NamedIsotope::D), false)]
    #[case::wildcard(AtomSymbol::WildcardAtom(WildcardAtom::Any), false)]
    #[case::heavy(AtomSymbol::WildcardAtom(WildcardAtom::Heavy), true)]
    #[case::atom_list(AtomSymbol::AtomList(AtomList::empty()), true)]
    fn test_extended_atom_has_extended_features_symbol(
        #[case] symbol: AtomSymbol,
        #[case] expected: bool,
    ) {
        assert_eq!(
            ExtendedAtom::from_atom_symbol(symbol).has_extended_features(),
            expected
        );
    }

    #[rstest]
    fn test_extended_atom_has_extended_features() {
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

    #[rstest]
    fn test_atom_try_from_roundtrip() {
        let atom = Atom {
            element: Some(Element::N),
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

        assert_eq!(atom2, atom);
    }
}
